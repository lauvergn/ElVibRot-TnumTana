!===========================================================================
!===========================================================================
!This file is part of ElVibRot.
!
!    ElVibRot is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ElVibRot is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ElVibRot.  If not, see <http://www.gnu.org/licenses/>.
!
!    Copyright 2015 David Lauvergnat [1]
!      with contributions of
!        Josep Maria Luis (optimization) [2]
!        Ahai Chen (MPI) [1,4]
!        Lucien Dupuy (CRP) [5]
!
![1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France
![2]: Institut de Química Computacional and Departament de Química,
!        Universitat de Girona, Catalonia, Spain
![3]: Department of Chemistry, Aarhus University, DK-8000 Aarhus C, Denmark
![4]: Maison de la Simulation USR 3441, CEA Saclay, France
![5]: Laboratoire Univers et Particule de Montpellier, UMR 5299,
!         Université de Montpellier, France
!
!    ElVibRot includes:
!        - Tnum-Tana under the GNU LGPL3 license
!        - Somme subroutines of John Burkardt under GNU LGPL license
!             http://people.sc.fsu.edu/~jburkardt/
!        - Somme subroutines of SHTOOLS written by Mark A. Wieczorek under BSD license
!             http://shtools.ipgp.fr
!        - Some subroutine of QMRPack (see cpyrit.doc) Roland W. Freund and Noel M. Nachtigal:
!             https://www.netlib.org/linalg/qmr/
!
!===========================================================================
!===========================================================================
MODULE mod_basis_BtoG_GtoB_SGType4_MPI
USE mod_system
USE mod_nDindex
USE mod_basis_set_alloc
USE mod_param_SGType2
USE mod_basis_RCVec_SGType4,ONLY:typervec,typecvec,alloc_typervec,alloc_typecvec,      &
                                 dealloc_typervec,dealloc_typecvec
IMPLICIT NONE

PRIVATE
PUBLIC  PackedBasis_TO_tabR_index_MPI,tabPackedBasis_TO_tabRpacked_MPI
PUBLIC  tabPackedBasis_TO_tabR_MPI,tabR_TO_tabPackedBasis_MPI
PUBLIC  Mapping_table_allocate_MPI,Mapping_table_MPI
PUBLIC  ini_iGs_MPI,set_iGs_MPI_mc,Set_scheme_MPI


CONTAINS
!=======================================================================================
!> @breif allocate mapping table for different MPI scheme
!---------------------------------------------------------------------------------------
SUBROUTINE Mapping_table_allocate_MPI(basis_SG,Max_Srep)
  USE mod_system
  USE mod_basis_set_alloc,ONLY:basis
  USE mod_MPI_aux
  IMPLICIT NONE

  TYPE(basis),                       intent(inout) :: basis_SG
  Integer,                           intent(in)    :: Max_Srep

#if(run_MPI)

  IF(MPI_id==0) THEN
    IF(MPI_scheme/=1) THEN
      CALL allocate_array(basis_SG%para_SGType2%tab_iB_OF_SRep_TO_iB,1,Max_Srep)
    ELSE
      CALL allocate_array(basis_SG%para_SGType2%tab_iB_OF_SRep_TO_iB,                  &
                          bounds_MPI(1,MPI_id),bounds_MPI(2,MPI_id))
    ENDIF
  ELSE
    CALL allocate_array(basis_SG%para_SGType2%tab_iB_OF_SRep_TO_iB,                    &
                        bounds_MPI(1,MPI_id),bounds_MPI(2,MPI_id))
  ENDIF

#endif
ENDSUBROUTINE Mapping_table_allocate_MPI
!=======================================================================================


!=======================================================================================
!> mapping table for MPI
!> 
!> @breif setup for different MPI scheme
!> scheme 1,4: keep on each threads the required mapping table 
!> scheme 2: keep mapping table only on root threads
!> scheme 3: keep required mapping table on levels2 threads
!---------------------------------------------------------------------------------------
SUBROUTINE Mapping_table_MPI(basis_SG,Max_Srep)
  USE mod_system
  USE mod_basis_set_alloc,ONLY:basis
  USE mod_MPI_aux
  IMPLICIT NONE

  TYPE(basis),                       intent(inout) :: basis_SG
  Integer,                           intent(in)    :: Max_Srep

#if(run_MPI)

  ! note, no boartcast here
  IF(MPI_scheme/=1) CALL MPI_combine_array(basis_SG%para_SGType2%tab_iB_OF_SRep_TO_iB)

  ! for scheme=2
  IF(MPI_id/=0 .AND. MPI_scheme/=1 .AND. MPI_scheme/=4) THEN
    IF(allocated(basis_SG%para_SGType2%tab_iB_OF_SRep_TO_iB))                          &
      deallocate(basis_SG%para_SGType2%tab_iB_OF_SRep_TO_iB)
  ENDIF

  ! for scheme=3
  IF(MPI_scheme==3) THEN
    IF(MPI_id==0) THEN
      DO i_MPI=1,MPI_np-1
        IF(mod(i_MPI,n_level2)==0) THEN
          CALL MPI_Send(basis_SG%para_SGType2%tab_iB_OF_SRep_TO_iB,Max_Srep,           &
                        MPI_Integer4,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
        ENDIF
      ENDDO
    ELSEIF(mod(MPI_id,n_level2)==0) THEN
      CALL allocate_array(basis_SG%para_SGType2%tab_iB_OF_SRep_TO_iB,1,Max_Srep)
      CALL MPI_Recv(basis_SG%para_SGType2%tab_iB_OF_SRep_TO_iB,Max_Srep,MPI_Integer4,  &
                    root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
    ENDIF ! MPI_id==0
  ENDIF ! MPI_scheme==3

#endif
END SUBROUTINE Mapping_table_MPI
!=======================================================================================


!=======================================================================================
!> transfer from compact rep to SRep
!---------------------------------------------------------------------------------------
SUBROUTINE tabPackedBasis_TO_tabR_MPI(PsiR,all_RvecB_temp,iG,SGType2,nDI_index,        &
                                 reduce_Vlength,size_psi,Max_nDI_ib0,nDI_index_list)
  USE mod_system
  USE mod_basis_set_alloc
  USE mod_param_SGType2
  USE mod_nDindex
  USE mod_MPI_aux
  IMPLICIT NONE

  TYPE (TypeRVec),allocatable,       intent(inout) :: PsiR(:)
  TYPE (param_SGType2),              intent(inout) :: SGType2
  Real(kind=Rkind),                     intent(in) :: all_RvecB_temp(:)
  Integer(kind=Ikind),allocatable,      intent(in) :: nDI_index(:)
  Integer(kind=Ikind),allocatable,   intent(inout) :: nDI_index_list(:)
  integer(kind=MPI_INTEGER_KIND),       intent(in) :: size_psi
  integer(kind=MPI_INTEGER_KIND),       intent(in) :: reduce_Vlength  
  integer,                              intent(in) :: iG
  integer,                              intent(in) :: Max_nDI_ib0

  integer,allocatable                              :: temp_list(:)  
  integer                                          :: temp_length 
  
  integer :: ib0,nb_AT_iG,iB_ib0,nDI_ib0,iBSRep,iB,nDI
  integer :: index_find,ii,itab
  Logical :: once1

#if(run_MPI)

  temp_int=SGType2%tab_nb_OF_SRep(iG)*SGType2%nb0
  DO itab=1,size_psi
    CALL allocate_array(PsiR(itab)%V,1,temp_int)
    !IF(allocated(PsiR(itab)%V)) deallocate(PsiR(itab)%V)
    !allocate(PsiR(itab)%V(temp_int))
  ENDDO
  CALL allocate_array(temp_list,1,temp_int)
  once1=.TRUE.
  
  DO itab=1,size_psi
    temp_length=0
    DO ib0=1,SGType2%nb0
      iBSRep=SGType2%tab_Sum_nb_OF_SRep(iG)-SGType2%tab_nb_OF_SRep(iG)
      DO iB_ib0=1,SGType2%tab_nb_OF_SRep(iG)
        iBSRep=iBSRep+1
        nDI_ib0=SGType2%tab_iB_OF_SRep_TO_iB(iBSRep)
        IF (nDI_ib0>0 .AND. nDI_ib0<=Max_nDI_ib0) THEN
          nDI=(ib0-1)*Max_nDI_ib0               +nDI_ib0
          iB =(ib0-1)*SGType2%tab_nb_OF_SRep(iG)+iB_ib0

          temp_length=temp_length+1
          IF(once1) THEN
            SGType2%V_allcount=SGType2%V_allcount+1
            temp_list(temp_length)=nDI_index_list(SGType2%V_allcount)
          ENDIF
          temp_int=(itab-1)*reduce_Vlength+temp_list(temp_length)
          PsiR(itab)%V(iB)=all_RvecB_temp(temp_int)
          !tabR_iG(iB) = tabR(nDI)
        END IF
      END DO
    END DO  
    once1=.FALSE.
  ENDDO
  deallocate(temp_list)

#endif
END SUBROUTINE tabPackedBasis_TO_tabR_MPI

!=======================================================================================

!=======================================================================================
!> subroutine for ccounting the length of compact basis to be send to each threads
!
!> reduce_index_mpi: count the overall length for each threads
!> nDI_index: temp index for all possible nDI for each threads
!> nDI_index_list: record the relevant position of each elements in nDI_index
!
!> Warning: time consuming and memory consuming, need improvements
!---------------------------------------------------------------------------------------
SUBROUTINE PackedBasis_TO_tabR_index_MPI(iG,SGType2,reduce_index_mpi,nDI_index,        &
                                         Max_nDI_ib0,nDI_index_list)
  USE mod_system
  USE mod_basis_set_alloc
  USE mod_param_SGType2
  USE mod_nDindex
  USE mod_MPI_aux
  IMPLICIT NONE

  TYPE(param_SGType2),            intent(inout) :: SGType2
  Integer,                           intent(in) :: iG
  Integer,                           intent(in) :: Max_nDI_ib0
  Integer(kind=MPI_INTEGER_KIND), intent(inout) :: reduce_index_mpi
  Integer(kind=Ikind),allocatable,intent(inout) :: nDI_index(:)
  Integer(kind=Ikind),allocatable,intent(inout) :: nDI_index_list(:)
  
  Integer(kind=Ikind),allocatable               :: nDI_index_temp(:)
  Integer                                       :: iBSRep,iB,nDI
  Integer                                       :: ib0,nb_AT_iG,iB_ib0,nDI_ib0  
  Integer                                       :: ii

#if(run_MPI)

  !Max_nDI_ib0 = size(tabR)/SGType2%nb0
  !size_SRep=SGType2%tab_nb_OF_SRep(iG)*SGType2%nb0
  
  DO ib0=1,SGType2%nb0
    iBSRep = SGType2%tab_Sum_nb_OF_SRep(iG)-SGType2%tab_nb_OF_SRep(iG)
    DO iB_ib0=1,SGType2%tab_nb_OF_SRep(iG)
      iBSRep  = iBSRep + 1
      nDI_ib0 = SGType2%tab_iB_OF_SRep_TO_iB(iBSRep)
      IF (nDI_ib0 > 0 .AND. nDI_ib0 <= Max_nDI_ib0 ) THEN
        nDI=(ib0-1)*Max_nDI_ib0               +nDI_ib0
        iB =(ib0-1)*SGType2%tab_nb_OF_SRep(iG)+iB_ib0  
        ! counts actual number of V
        SGType2%V_allcount=SGType2%V_allcount+1

        !> record the mapping list on non-root threads according to nDI-----------------
        IF(MPI_id/=0) THEN
          DO ii=1,MAX(reduce_index_mpi,1)
            IF(nDI_index(ii)==nDI) THEN
              nDI_index_list(SGType2%V_allcount)=ii
              EXIT
            ENDIF
          ENDDO
        ENDIF ! for MPI_id/=0-----------------------------------------------------------
        
        ! if nDI is not on the current list, include the new nDI in the "nDI_index"
        IF(ii>MAX(reduce_index_mpi,1)) THEN
          reduce_index_mpi=reduce_index_mpi+1
          nDI_index(reduce_index_mpi)=nDI
          IF(MPI_id/=0) nDI_index_list(SGType2%V_allcount)=reduce_index_mpi
          
          !> if exceed the size of current array----------------------------------------
          IF(reduce_index_mpi==SGType2%num_nDI_index) THEN
            !> for a banlance of new allocation and memory waste
            SGType2%num_nDI_index=SGType2%num_nDI_index+Max(SGType2%num_nDI_index/5,500) 
            !> nDI_index_temp -> nDI_index
            ! nDI_index_temp is deallocated automatically
            CALL allocate_array(nDI_index_temp,1,SGType2%num_nDI_index)
            nDI_index_temp(1:reduce_index_mpi)=nDI_index
            CALL move_alloc(nDI_index_temp,nDI_index) 
          ENDIF ! for reduce_index_mpi==SGType2%num_nDI_index---------------------------
          
        ENDIF ! for ii>MAX(reduce_index_mpi,1)
      ENDIF ! for nDI_ib0 > 0 .AND. nDI_ib0 <= Max_nDI_ib0
    ENDDO
  ENDDO

#endif
END SUBROUTINE PackedBasis_TO_tabR_index_MPI 
!=======================================================================================

!=======================================================================================
! pack the vectors to send to each threads (not currently used)
!---------------------------------------------------------------------------------------

SUBROUTINE tabPackedBasis_TO_tabRpacked_MPI(all_RvecB_temp,Psi_RvecB,nDI_index,length)
  IMPLICIT NONE
  
  Real(kind=Rkind), allocatable,intent(inout) :: all_RvecB_temp(:)
  Real(kind=Rkind),             intent(in)    :: Psi_RvecB(:)
  Integer,allocatable,          intent(in)    :: nDI_index(:)
  Integer,                      intent(in)    :: length
  Integer ii  

#if(run_MPI)

  Do ii=1,length
    all_RvecB_temp(ii)=Psi_RvecB(nDI_index(ii))
  ENDDO

#endif 
END SUBROUTINE tabPackedBasis_TO_tabRpacked_MPI

!=======================================================================================

!=======================================================================================
! transfer from SRep to compact basis
!---------------------------------------------------------------------------------------
SUBROUTINE tabR_TO_tabPackedBasis_MPI(all_RvecB_temp2,PsiR,iG,SGType2,WeightiG,        &
                           nDI_index,reduce_Vlength,size_psi,Max_nDI_ib0,nDI_index_list)
  USE mod_system
  USE mod_basis_set_alloc
  USE mod_param_SGType2
  USE mod_nDindex
  USE mod_MPI_aux
  IMPLICIT NONE

  TYPE(TypeRVec),allocatable,         intent(inout) :: PsiR(:)
  TYPE(param_SGType2),                intent(inout) :: SGType2
  Real(kind=Rkind),                   intent(inout) :: all_RvecB_temp2(:)
  Integer(kind=Ikind),allocatable,       intent(in) :: nDI_index(:)
  Integer(kind=Ikind),allocatable,    intent(inout) :: nDI_index_list(:)
  integer,                               intent(in) :: iG
  integer(kind=MPI_INTEGER_KIND),        intent(in) :: size_psi
  integer,                               intent(in) :: Max_nDI_ib0
  integer(kind=MPI_INTEGER_KIND)        ,intent(in) :: reduce_Vlength  
  real(kind=Rkind),                      intent(in) :: WeightiG

  integer,allocatable                               :: temp_list(:)  
  integer                                           :: temp_length 
   
  integer :: ii,itab
  integer :: iBSRep,iB,nDI
  integer :: ib0,nb_AT_iG,iB_ib0,nDI_ib0
  Logical :: once1

#if(run_MPI)

  IF(size(PsiR) == 0) STOP 'ERROR in tabR_TO_tabPackedBasis_MPI'
  
  temp_int=SGType2%tab_nb_OF_SRep(iG)*SGType2%nb0
  CALL allocate_array(temp_list,1,temp_int)
  once1=.TRUE.
  
  DO itab=1,size_psi
    temp_length=0
    DO ib0=1,SGType2%nb0
      iBSRep=SGType2%tab_Sum_nb_OF_SRep(iG)-SGType2%tab_nb_OF_SRep(iG)
      DO iB_ib0=1,SGType2%tab_nb_OF_SRep(iG)
        iBSRep =iBSRep + 1
        nDI_ib0=SGType2%tab_iB_OF_SRep_TO_iB(iBSRep)
        IF(nDI_ib0>0 .AND. nDI_ib0<=Max_nDI_ib0) THEN
          nDI=(ib0-1)*Max_nDI_ib0               +nDI_ib0
          iB =(ib0-1)*SGType2%tab_nb_OF_SRep(iG)+iB_ib0
          temp_length=temp_length+1
          IF(once1) THEN
            SGType2%V_allcount2=SGType2%V_allcount2+1
            temp_list(temp_length)=nDI_index_list(SGType2%V_allcount2)
          ENDIF
          temp_int=(itab-1)*reduce_Vlength+temp_list(temp_length)
          all_RvecB_temp2(temp_int)=all_RvecB_temp2(temp_int)+WeightiG*PsiR(itab)%V(iB)
          !tabR(nDI)=tabR(nDI)+WeightiG*tabR_iG(iB)
        END IF
      END DO
    END DO
    once1=.FALSE.
  ENDDO
  deallocate(temp_list)

#endif
END SUBROUTINE tabR_TO_tabPackedBasis_MPI
!=======================================================================================

!=======================================================================================
SUBROUTINE Set_scheme_MPI_old(basis_SG)
  USE mod_basis_set_alloc,ONLY:basis
  IMPLICIT NONE

  TYPE(basis),                   intent(inout) :: basis_SG
  Integer                                      :: ave_iGs

#if(run_MPI)

  IF(MPI_scheme==0) THEN
    ave_iGs=( sum(basis_SG%para_SGType2%tab_nb_OF_SRep(:))                         &
                 *basis_SG%para_SGType2%nb0 )/MPI_np
    IF(ave_iGs<1000000) THEN
      MPI_scheme=4
    ELSE
      IF(MPI_np<Num_L2) THEN
        MPI_scheme=2
      ELSE
        MPI_scheme=3
      ENDIF
    ENDIF
  ENDIF

  ! initialize for MPI_scheme 3
  IF(MPI_scheme==3) THEN
    n_level2=MPI_np/(MPI_np/Num_L2+1)
    IF(mod(MPI_np,MPI_np/Num_L2+1)/=0) n_level2=n_level2+1
    ! write(out_unitp,*) 'MPI_scheme3 check: Num_L2=',Num_L2,'n_level2=',n_level2
  ENDIF

  ! initialize iGs distribution
  CALL ini_iGs_MPI(basis_SG,.FALSE.)

  ! if keeper compact vec on current thread
  IF(MPI_id==0 .OR. MPI_scheme==1) keep_MPI=.TRUE.

#endif
END SUBROUTINE Set_scheme_MPI_old
!=======================================================================================

!=======================================================================================
SUBROUTINE Set_scheme_MPI(basis_SG,lMax_Srep)
  USE mod_basis_set_alloc,ONLY:basis
  IMPLICIT NONE

  TYPE(basis),                   intent(inout) :: basis_SG
  Integer(kind=ILkind),          intent(in)    :: lMax_Srep
  Integer                                      :: ave_iGs

#if(run_MPI)

  IF(MPI_scheme==0) THEN
    ! do selection depend on the memory
    ! if do have big enough memory, set MPI_scheme=1
    IF(basis_SG%nDindB%Max_nDI*MPI_nb_WP<lMax_Srep/2) THEN
      MPI_scheme=1
    ELSE
      IF(MPI_np<Num_L2) THEN
        MPI_scheme=2
      ELSE
        MPI_scheme=3
      ENDIF
    ENDIF
  ENDIF

  ! initialize for MPI_scheme 3
  IF(MPI_scheme==3) THEN
    n_level2=MPI_np/(MPI_np/Num_L2+1)
    IF(mod(MPI_np,MPI_np/Num_L2+1)/=0) n_level2=n_level2+1
    ! write(out_unitp,*) 'MPI_scheme3 check: Num_L2=',Num_L2,'n_level2=',n_level2
  ENDIF

  ! initialize iGs distribution
  CALL ini_iGs_MPI(basis_SG,.FALSE.)

  ! if keeper compact vec on current thread
  IF(MPI_id==0 .OR. MPI_scheme==1) keep_MPI=.TRUE.

#endif
END SUBROUTINE Set_scheme_MPI
!=======================================================================================

!=======================================================================================
!> initialize iGs_MPI for each threads
!> note to .false. BasisnD%para_SGType2%once_action to perform once
!=======================================================================================
  SUBROUTINE ini_iGs_MPI(BasisnD,once)
    USE mod_basis_set_alloc,ONLY:basis
    USE mod_MPI_aux
    IMPLICIT NONE

    TYPE(basis),                    intent(inout) :: BasisnD
    Logical,                           intent(in) :: once

#if(run_MPI)

    nb_per_MPI=BasisnD%para_SGType2%nb_SG/MPI_np
    nb_rem_MPI=mod(BasisnD%para_SGType2%nb_SG,MPI_np) 

    CALL allocate_array(iGs_MPI,1,2,0,MPI_np-1)
    CALL allocate_array(iGs_MPI0,1,2,0,MPI_np-1)

    DO i_MPI=0,MPI_np-1
      bound1_MPI=i_MPI*nb_per_MPI+1+MIN(i_MPI,nb_rem_MPI)
      bound2_MPI=(i_MPI+1)*nb_per_MPI+MIN(i_MPI,nb_rem_MPI)+merge(1,0,nb_rem_MPI>i_MPI)
      iGs_MPI(1,i_MPI)=bound1_MPI
      iGs_MPI(2,i_MPI)=bound2_MPI
    ENDDO

    iGs_MPI0=iGs_MPI

    CALL set_iGs_MPI_mc(BasisnD,initialize=.TRUE.)

    IF(once) BasisnD%para_SGType2%once_action=.FALSE.

#endif
  ENDSUBROUTINE ini_iGs_MPI  
!=======================================================================================


!=======================================================================================
!> @brief 
!> set iGs_MPI_mc for reducing memory usage in MPI action
!> reduce level depends on MPI_mc
!> 
!> @param BasisnD [in]
!> @param initialize [in], optional, for initializing iGs_MPI_mc
!> @param iGs_MPI [in] 
!> @param iGs_MPI_mc [inout]
!=======================================================================================
SUBROUTINE set_iGs_MPI_mc(BasisnD,initialize)
  USE mod_basis_set_alloc,ONLY:basis
  USE mod_MPI
  IMPLICIT NONE

  TYPE(basis),                       intent(in) :: BasisnD
  Logical,optional,                  intent(in) :: initialize
!  Integer,                        intent(inout) :: iGs_MPI_mc(2,MPI_mc,0:MPI_np-1)

  Integer                                       :: ii
  Integer                                       :: d1
  Integer                                       :: d2
  Integer                                       :: temp_int

#if(run_MPI)

  IF(present(initialize)) THEN
    IF(initialize) THEN
      allocate(iGs_MPI_mc(1:2,1:MPI_mc,0:MPI_np-1))
      allocate(iGs_MPI_mc0(1:2,1:MPI_mc,0:MPI_np-1))
    ENDIF
  ENDIF

  iGs_MPI_mc(:,:,:)=0

  DO i_MPI=0,MPI_np-1
    d1=iGs_MPI(1,i_MPI)
    d2=iGs_MPI(2,i_MPI)
    nb_per_MPI=(d2-d1+1)/MPI_mc
    nb_rem_MPI=mod(d2-d1+1,MPI_mc) 
    DO ii=1,MPI_mc
      bound1_MPI=(ii-1)*nb_per_MPI+d1+MIN(ii-1,nb_rem_MPI)
      bound2_MPI=ii*nb_per_MPI+(d1-1)+MIN(ii-1,nb_rem_MPI)+merge(1,0,nb_rem_MPI>(ii-1))
      iGs_MPI_mc(1,ii,i_MPI)=bound1_MPI
      iGs_MPI_mc(2,ii,i_MPI)=bound2_MPI
    ENDDO
  ENDDO

  IF(present(initialize)) THEN
    IF(initialize) iGs_MPI_mc0=iGs_MPI_mc
  ENDIF

#endif
ENDSUBROUTINE set_iGs_MPI_mc 
!=======================================================================================

END MODULE mod_basis_BtoG_GtoB_SGType4_MPI
