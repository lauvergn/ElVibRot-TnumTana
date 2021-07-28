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
MODULE mod_OpPsi_SG4_MPI
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: sub_TabOpPsi_FOR_SGtype4_MPI
  PUBLIC :: sub_TabOpPsi_FOR_SGtype4_SRB_MPI,sub_TabOpPsi_FOR_SGtype4_SRG_MPI

  CONTAINS

!=======================================================================================
!> OpPsi is assigned in this routine with OpPsi%Vec=0 initially
!> after this subroutine, OpPsi will be ready
!======================================================================================= 
SUBROUTINE sub_TabOpPsi_FOR_SGtype4_MPI(Psi,OpPsi,para_Op)
  USE mod_system
  USE mod_nDindex
  USE mod_Coord_KEO,              ONLY:CoordType
  USE mod_SymAbelian,             ONLY:Calc_symab1_EOR_symab2
  USE mod_basis_set_alloc,        ONLY:basis
  USE mod_basis_BtoG_GtoB_SGType4,ONLY:tabPackedBasis_TO_tabR_AT_iG,                   &
                                       tabR_AT_iG_TO_tabPackedBasis,                   &
                                       TypeRVec,dealloc_TypeRVec
  USE mod_psi,                    ONLY:param_psi,ecri_psi,                             &
                                       Set_symab_OF_psiBasisRep
  USE mod_SetOp,                  ONLY:param_Op,write_param_Op
  IMPLICIT NONE

  TYPE(param_psi),                       intent(in)    :: Psi(:)
  TYPE(param_psi),                       intent(inout) :: OpPsi(:)
  TYPE(param_Op),                        intent(inout) :: para_Op

!  TYPE(CoordType),pointer                              :: mole
  TYPE(basis),pointer                                   :: BasisnD  
  Integer                                               :: ii
  Integer                                               :: itab
  Integer                                               :: iterm00
  Integer                                               :: OpPsi_symab
  Integer,allocatable                                   :: tab_l(:)
  Character(len=*),parameter                  :: name_sub='sub_TabOpPsi_FOR_SGtype4_MPI'
  
#if(run_MPI)

!  mole    => para_Op%mole
  BasisnD => para_Op%BasisnD

  IF(MPI_id==0) THEN
    IF(size(Psi)==0) THEN
      write(out_unitp,*) 'ERROR in ',name_sub
      write(out_unitp,*) '  The size of Psi(:) is zero!!'
      write(out_unitp,*) '  => Check the fortran.'
      STOP ' ERROR in sub_TabOpPsi_FOR_SGtype4_MPI: size(Psi) = 0'
    END IF
    IF (Psi(1)%cplx) THEN
      write(out_unitp,*) 'ERROR in ',name_sub
      write(out_unitp,*) '  Psi(1) is complex !!'
      write(out_unitp,*) '  => Check the fortran.'
      STOP ' ERROR in sub_TabOpPsi_FOR_SGtype4_MPI: Psi(1) is complex'
    END IF
  ENDIF
  !IF (para_Op%nb_bie /= 1) STOP 'nb_bie /= 1 in sub_TabOpPsi_FOR_SGtype4'

  IF(.NOT. allocated(iGs_MPI)) STOP 'error in sub_TabOpPsi_FOR_SGtype4_MPI,            &
                                     iGs_MPI not initialized'

  IF(keep_MPI) THEN
    DO itab=1,size(Psi)
      OpPsi(itab)         =Psi(itab) ! for the allocation. It has to be changed!
      OpPsi(itab)%RvecB(:)=ZERO
    END DO
  ENDIF

  CALL Action_MPI(Psi,OpPsi,BasisnD,para_Op)

  IF (print_level > 0 .AND. BasisnD%para_SGType2%nb_SG > 10**4) THEN
    write(out_unitp,'(a)',ADVANCE='yes') '-]'
  END IF
  CALL flush_perso(out_unitp)

  iterm00 = para_Op%derive_term_TO_iterm(0,0)
  IF (associated(para_Op%OpGrid)) THEN
    para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_done =                  &
                      para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid
    para_Op%OpGrid(iterm00)%para_FileGrid%Save_FileGrid_done =                 &
                      para_Op%OpGrid(iterm00)%para_FileGrid%Save_FileGrid
    ! grid calculated and saved in first action
    para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_done=.True. 
  END IF

  IF(keep_MPI) THEN
    DO ii=1,size(OpPsi)
      OpPsi_symab = Calc_symab1_EOR_symab2(para_Op%symab,Psi(ii)%symab)
      CALL Set_symab_OF_psiBasisRep(OpPsi(ii),OpPsi_symab)
      !write(out_unitp,*) 'para_Op,psi symab ',i,para_Op%symab,Psi(i)%symab
      !write(out_unitp,*) 'OpPsi_symab',i,OpPsi(i)%symab
    END DO
  ENDIF

#endif
  END SUBROUTINE sub_TabOpPsi_FOR_SGtype4_MPI
!=======================================================================================


!=======================================================================================
!> @brief operator action using MPI
!=======================================================================================
  SUBROUTINE Action_MPI(Psi,OpPsi,BasisnD,para_Op)
    USE mod_system
    USE mod_nDindex
    USE mod_psi,            ONLY:param_psi
    USE mod_SetOp,          ONLY:param_Op
    USE mod_basis_set_alloc,ONLY:basis
    USE mod_MPI_aux
    IMPLICIT NONE

    TYPE(param_psi),                       intent(in)    :: Psi(:)
    TYPE(param_psi),                       intent(inout) :: OpPsi(:)
    TYPE(basis),pointer,                   intent(inout) :: BasisnD
    TYPE(param_Op),                        intent(inout) :: para_Op

    Integer(kind=MPI_INTEGER_KIND),pointer               :: size_ST(:)
    Integer(kind=MPI_INTEGER_KIND),pointer               :: size_ST_mc(:,:)

#if(run_MPI)

    CALL system_clock(time_point1,time_rate,time_max)
    IF(BasisnD%para_SGType2%once_action) THEN 
      CALL time_perso('MPI loop in action begin')

      !CALL ini_iGs_MPI(BasisnD,.FALSE.)

      CALL allocate_array(time_MPI_local,0,MPI_np-1)
      CALL allocate_array(time_MPI_calcu,0,MPI_np-1)
      CALL allocate_array(BasisnD%para_SGType2%size_ST,0,MPI_np-1)
      CALL allocate_array(BasisnD%para_SGType2%size_ST_mc,1,MPI_mc,0,MPI_np-1)

      write(out_unitp,*) 'nb_Term,nb_SG',para_Op%nb_Term,BasisnD%para_SGType2%nb_SG
      IF(MPI_scheme/=1) THEN
        CALL allocate_array(para_Op%OpGrid(1)%para_FileGrid%Save_MemGrid_iG,           &
                            1,BasisnD%para_SGType2%nb_SG)
      ENDIF
    ENDIF

    ! size of Smolyak terms for each thread
    size_ST   =>BasisnD%para_SGType2%size_ST
    size_ST_mc=>BasisnD%para_SGType2%size_ST_mc
    IF(BasisnD%para_SGType2%once_action) THEN
      CALL get_size_ST(BasisnD,size_ST,size_ST_mc)
      write(out_unitp,*) 'size_ST:',size_ST(MPI_id),'from',MPI_id
      !write(out_unitp,*) 'size_ST_mc:',size_ST_mc(:,MPI_id),'from',MPI_id
    ENDIF
    
    !-----------------------------------------------------------------------------------
    !> scheme 1: keep compact basis on all threads, best for small compact basis
    !> scheme 2: most general MPI, keep compact basis on master only, 
    !!           provide control of Smolyak terms' memory, 
    !!           auto balance Smolyak terms on different threads
    !> scheme 3: 2-level Smolyak terms distribution based on scheme 2
    !> scheme 4: for more cores and more memory, need improvement
    !-----------------------------------------------------------------------------------
    SELECT CASE(MPI_scheme)
      CASE (1)
        CALL Action_MPI_S1(Psi,OpPsi,BasisnD,para_Op)
      CASE (2)
        CALL Action_MPI_S2(Psi,OpPsi,BasisnD,para_Op)
      CASE (3)
        CALL Action_MPI_S3(Psi,OpPsi,BasisnD,para_Op)
      CASE (4)
        CALL Action_MPI_S4(Psi,OpPsi,BasisnD,para_Op)
      CASE Default
        STOP 'error in MPI scheme for operator action'
    END SELECT
    !-----------------------------------------------------------------------------------

    CALL system_clock(time_point2,time_rate,time_max)
    IF(BasisnD%para_SGType2%once_action) CALL allocate_array(time_MPI_act_all,0,MPI_np-1)
    time_MPI_act_all(MPI_id)=merge(time_point2-time_point1,                            &
                              time_point2-time_point1+time_max,time_point2>=time_point1)
    time_MPI_action=time_MPI_action+time_MPI_act_all(MPI_id)
    
    IF(BasisnD%para_SGType2%once_action .AND. MPI_id==0)                               &
                   write(out_unitp,*) 'time MPI comm check: ',time_comm,' from ', MPI_id
    IF(BasisnD%para_SGType2%once_action) CALL time_perso('MPI loop in action end') 
    BasisnD%para_SGType2%once_action=.FALSE.

    write(out_unitp,*) 'time used in action on processor :', MPI_id,                   &    
                    real(time_MPI_act_all(MPI_id),kind=Rkind)/real(time_rate,kind=Rkind)
#endif
  END SUBROUTINE Action_MPI
!=======================================================================================


!#if(run_MPI)
!!=======================================================================================  
!!> action with MPI: scheme 1
!!======================================================================================= 
!  SUBROUTINE Action_MPI_S1(Psi,OpPsi,BasisnD,para_Op)
!    USE mod_system
!    USE mod_nDindex
!    USE mod_Coord_KEO,                  ONLY:CoordType
!    USE mod_basis_set_alloc,            ONLY:basis
!    USE mod_OpPsi_SG4,                  ONLY:sub_TabOpPsi_OF_ONEDP_FOR_SGtype4
!    USE mod_basis_BtoG_GtoB_SGType4,    ONLY:tabPackedBasis_TO_tabR_AT_iG,            &
!                                             tabR_AT_iG_TO_tabPackedBasis,            &
!                                             TypeRVec,dealloc_TypeRVec
!    USE mod_basis_BtoG_GtoB_SGType4_MPI,ONLY:set_iGs_MPI_mc
!    USE mod_psi,                        ONLY:param_psi
!    USE mod_SetOp,                      ONLY:param_Op
!    USE mod_MPI_aux
!    IMPLICIT NONE
!
!    TYPE(param_psi),                       intent(in)    :: Psi(:)
!    TYPE(param_psi),                       intent(inout) :: OpPsi(:)
!    TYPE(basis),pointer,                   intent(inout) :: BasisnD
!    TYPE(param_Op),                        intent(inout) :: para_Op
!
!    TYPE(TypeRVec),allocatable                           :: PsiR(:)
!    Real(kind=Rkind),allocatable                         :: PsiR_temp(:) 
!    Real(kind=Rkind),allocatable                         :: Psi_ST(:) 
!    Integer(kind=MPI_INTEGER_KIND)                       :: PsiR_temp_length(0:MPI_np-1)
!    Integer(kind=MPI_INTEGER_KIND)                       :: Psi_ST_length(0:MPI_np-1)
!    Integer(kind=MPI_INTEGER_KIND),pointer               :: size_psi
!    Integer(kind=MPI_INTEGER_KIND),pointer               :: size_ST(:) 
!    Integer                                              :: iG
!    Integer                                              :: ST_index1
!    Integer                                              :: ST_index2
!    Integer                                              :: ST_index_iG
!    Integer                                              :: size_ST_iG
!    Integer                                              :: itab
!    Integer                                              :: time_auto1
!    Integer                                              :: time_auto2
!    Integer                                              :: d1
!    Integer                                              :: d2
!    Logical,pointer                                      :: once_action
!
!    once_action => BasisnD%para_SGType2%once_action
!    size_ST     => BasisnD%para_SGType2%size_ST
!    
!    IF(once_action) THEN 
!      write(out_unitp,*) 'action with MPI: Scheme 1'
!    ELSE
!      ! auto distribute Smolyak terms to different threads 
!      write(out_unitp,*) '------------------------------------------------'
!      write(out_unitp,*) 'iGs_MPI at current step: ',iGs_MPI
!
!      IF(MPI_iGs_auto .AND. MPI_np>1) CALL auto_iGs_MPI(para_Op)
!      write(out_unitp,*) 'iGs_MPI at new step    : ',iGs_MPI
!      write(out_unitp,*) '------------------------------------------------'
!      
!      ! get new size_ST according to new iGs_MPI
!      IF(MPI_iGs_auto) THEN
!        Do i_MPI=0,MPI_np-1
!          size_ST(i_MPI)=0;
!          DO iG=iGs_MPI(1,i_MPI),iGs_MPI(2,i_MPI)
!            temp_int=BasisnD%para_SGType2%tab_nb_OF_SRep(iG)*BasisnD%para_SGType2%nb0
!            size_ST(i_MPI)=size_ST(i_MPI)+temp_int
!          ENDDO
!        ENDDO
!        
!        !CALL set_iGs_MPI_mc(BasisnD)
!      ENDIF
!    ENDIF
!
!    size_psi => BasisnD%para_SGType2%size_psi  
!    If(MPI_id==0) size_psi=INT(size(Psi),MPI_INTEGER_KIND) ! for itab
!    CALL MPI_BCAST(size_psi,size1_MPI,Int_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
!
!    IF(allocated(PsiR)) deallocate(PsiR)
!    allocate(PsiR(size_psi))
!    
!    time_MPI_local(MPI_id)=0
!
!    !-----------------------------------------------------------------------------------
!    IF(MPI_id==0) THEN
!      !-prepare PsiR(itab)%V to be send to other threads--------------------------------
!      CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,1)
!      DO i_MPI=1,MPI_np-1
!        CALL allocate_array(Psi_ST,1,size_psi*size_ST(i_MPI))
!        ST_index1=0
!        DO iG=iGs_MPI(1,i_MPI),iGs_MPI(2,i_MPI)
!          DO itab=1,size_psi
!            CALL tabPackedBasis_TO_tabR_AT_iG(PsiR(itab)%V,psi(itab)%RvecB,iG,         &
!                                              BasisnD%para_SGType2)
!            ST_index2=ST_index1+size(PsiR(itab)%V)
!            Psi_ST(ST_index1+1:ST_index2)=PsiR(itab)%V
!            ST_index1=ST_index2
!          ENDDO
!        ENDDO ! for iG
!        Psi_ST_length(i_MPI)=ST_index1
!        IF(once_action) write(out_unitp,*) 'length check:',i_MPI,size_ST(i_MPI),       &
!                 size_psi,size(psi(1)%RvecB),sum(BasisnD%para_SGType2%tab_nb_OF_SRep(:))
!
!        ! double check
!        IF(Abs(size_psi*size_ST(i_MPI)-Psi_ST_length(i_MPI))>0) THEN
!          write(out_unitp,*) 'error in MPI action part,check length'
!          STOP
!        ENDIF
!
!        CALL time_record(time_comm,time_temp1,time_temp2,1)
!        CALL MPI_Send(Psi_ST,size_psi*size_ST(i_MPI),Real_MPI,i_MPI,i_MPI,             &
!                      MPI_COMM_WORLD,MPI_err)
!        CALL time_record(time_comm,time_temp1,time_temp2,2)
!      ENDDO
!      CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,2)
!
!      !---------------------------------------------------------------------------------
!
!      !-calculations on MPI_id=0--------------------------------------------------------
!      DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)
!        !-transfore to SRep-------------------------------------------------------------
!        DO itab=1,size_psi
!          CALL tabPackedBasis_TO_tabR_AT_iG(PsiR(itab)%V,psi(itab)%RvecB,iG,           &
!                                            BasisnD%para_SGType2)
!        ENDDO
!
!        !-main calculation--------------------------------------------------------------
!        CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4(PsiR,iG,                                &
!                          BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval(:,iG),para_Op)  
!
!        !-back to compact form----------------------------------------------------------
!        DO itab=1,size_psi
!          CALL tabR_AT_iG_TO_tabPackedBasis(OpPsi(itab)%RvecB,PsiR(itab)%V,iG,         &
!                                            BasisnD%para_SGType2,BasisnD%WeightSG(iG))
!        ENDDO
!      ENDDO ! main loop of iG for calcuation on master
!
!      !-receive results from other threads----------------------------------------------
!      CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,1)
!      DO i_MPI=1,MPI_np-1
!        CALL time_record(time_comm,time_temp1,time_temp2,1)
!        CALL allocate_array(Psi_ST,1,size_psi*size_ST(i_MPI))
!        CALL MPI_Recv(Psi_ST,size_psi*size_ST(i_MPI),MPI_REAL8,i_MPI,                  &
!                      i_MPI,MPI_COMM_WORLD,MPI_stat,MPI_err)
!        CALL time_record(time_comm,time_temp1,time_temp2,2)
!        ST_index1=0
!        DO iG=iGs_MPI(1,i_MPI),iGs_MPI(2,i_MPI)
!          size_ST_iG=BasisnD%para_SGType2%tab_nb_OF_SRep(iG)*BasisnD%para_SGType2%nb0
!          DO itab=1,size_psi
!            ST_index2=ST_index1+size_ST_iG
!            d1=ST_index1+1
!            d2=ST_index2
!            CALL tabR_AT_iG_TO_tabPackedBasis(OpPsi(itab)%RvecB,Psi_ST(d1:d2),iG,      &
!                                              BasisnD%para_SGType2,BasisnD%WeightSG(iG))
!            ST_index1=ST_index2
!          ENDDO
!        ENDDO
!      ENDDO  ! for i_MPI=1,MPI_np-1  
!      CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,2)
!    ENDIF ! for MPI_id==0
!    
!    !-calculation on other threads------------------------------------------------------
!    IF(MPI_id/=0) THEN
!      ! get PsiR(itab)%V form root thread
!      CALL allocate_array(Psi_ST,1,size_psi*size_ST(MPI_id))
!      CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,1)
!      CALL MPI_Recv(Psi_ST,size_psi*size_ST(MPI_id),MPI_REAL8,root_MPI,MPI_id,         &
!                    MPI_COMM_WORLD,MPI_stat,MPI_err)
!      CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,2)
!
!      !-loop for main calculation-------------------------------------------------------
!      ST_index1=0
!      ST_index_iG=0
!      DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)
!        size_ST_iG=BasisnD%para_SGType2%tab_nb_OF_SRep(iG)*BasisnD%para_SGType2%nb0
!
!        !-extract SRep from PsiR_temp---------------------------------------------------
!        DO itab=1,size_psi
!          CALL allocate_array(PsiR(itab)%V,1,size_ST_iG)
!          ST_index2=ST_index1+size_ST_iG
!          PsiR(itab)%V=Psi_ST(ST_index1+1:ST_index2)
!          ST_index1=ST_index2
!        ENDDO
!
!        !-main calculation--------------------------------------------------------------
!        CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4(PsiR,iG,                                &
!                        BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval(:,iG),para_Op)  
!
!        !-pack PsiR(itab)%V-------------------------------------------------------------
!        Do itab=1,size_psi
!          d1=ST_index_iG+1
!          d2=ST_index_iG+size(PsiR(itab)%V)
!          Psi_ST(d1:d2)=PsiR(itab)%V
!          ST_index_iG=d2
!        ENDDO
!      ENDDO ! for iG
!      
!      CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,1)
!      !-send back PsiR(itab)%V----------------------------------------------------------
!      CALL MPI_Send(Psi_ST,size_psi*size_ST(MPI_id),MPI_REAL8,root_MPI,MPI_id,         &
!                    MPI_COMM_WORLD,MPI_err)
!      CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,2)
!    ENDIF ! for MPI_id/=0
!    !-----------------------------------------------------------------------------------
!
!    If(allocated(Psi_ST)) deallocate(Psi_ST)
!    DO itab=1,size_psi
!      CALL dealloc_TypeRVec(PsiR(itab))
!    END DO
!
!  END SUBROUTINE Action_MPI_S1
!!======================================================================================= 
!#endif

!=======================================================================================  
!> @brief action with MPI: scheme 1
!> @note keep compact RvecB on all threads
!> 
!> @param Psi [in] |\psi>
!> @param OpPsi [inout] O|\psi>
!> @param BasisnD [in] base infortmation
!> @param para_Op [in] Operator
!======================================================================================= 
  SUBROUTINE Action_MPI_S1(Psi,OpPsi,BasisnD,para_Op)
    USE mod_system
    USE mod_nDindex
    USE mod_psi,                        ONLY:param_psi
    USE mod_SetOp,                      ONLY:param_Op
    USE mod_Coord_KEO,                  ONLY:CoordType
    USE mod_basis_set_alloc,            ONLY:basis
    USE mod_OpPsi_SG4,                  ONLY:sub_TabOpPsi_OF_ONEDP_FOR_SGtype4
    USE mod_basis_BtoG_GtoB_SGType4,    ONLY:tabPackedBasis_TO_tabR_AT_iG,             &
                                             tabR_AT_iG_TO_tabPackedBasis,             &
                                             TypeRVec,dealloc_TypeRVec
    USE mod_basis_BtoG_GtoB_SGType4_MPI,ONLY:set_iGs_MPI_mc
    USE mod_MPI_aux
    IMPLICIT NONE

    TYPE(param_psi),                       intent(in)    :: Psi(:)
    TYPE(param_psi),                       intent(inout) :: OpPsi(:)
    TYPE(basis),pointer,                   intent(inout) :: BasisnD
    TYPE(param_Op),                        intent(inout) :: para_Op

    TYPE(TypeRVec),allocatable                           :: PsiR(:)
    Real(kind=Rkind),allocatable                         :: Vec(:)
    Integer(kind=MPI_INTEGER_KIND),pointer               :: size_psi
    Integer(kind=MPI_INTEGER_KIND),pointer               :: size_ST(:)
    Integer                                              :: size_RvecB
    Integer                                              :: iG
    Integer                                              :: itab
    Integer                                              :: d1
    Integer                                              :: d2
    Logical,pointer                                      :: once_action

#if(run_MPI)

    size_psi    => BasisnD%para_SGType2%size_psi
    size_ST     => BasisnD%para_SGType2%size_ST
    once_action => BasisnD%para_SGType2%once_action

    size_psi=INT(size(Psi),MPI_INTEGER_KIND) 
    size_RvecB=size(psi(1)%RvecB)

    IF(once_action) THEN 
      write(out_unitp,*) 'action with MPI: Scheme 1'
      write(out_unitp,*) '------------------------------------------------'
      !write(out_unitp,*) 'iGs_MPI: ',iGs_MPI
      !write(out_unitp,*) '------------------------------------------------'
      write(out_unitp,*) 'length check:',MPI_id,size_psi,size_RvecB,size_ST(MPI_id),   &
                                         sum(BasisnD%para_SGType2%tab_nb_OF_SRep(:))
    ENDIF ! once_action

    IF(allocated(PsiR)) deallocate(PsiR)
    allocate(PsiR(size_psi))

    CALL allocate_array(vec,1,size_RvecB*size_psi)

    !-calculation on each thread--------------------------------------------------------
    DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)

      !-transfore to SRep---------------------------------------------------------------
      DO itab=1,size_psi
        CALL tabPackedBasis_TO_tabR_AT_iG(PsiR(itab)%V,psi(itab)%RvecB,iG,             &
                                          BasisnD%para_SGType2)
      ENDDO

      !-main calculation----------------------------------------------------------------
      CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4(PsiR,iG,                                  &
                      BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval(:,iG),para_Op)  

      !-back to compact form------------------------------------------------------------
      DO itab=1,size_psi
        d1=size_RvecB*(itab-1)+1
        d2=size_RvecB* itab
        CALL tabR_AT_iG_TO_tabPackedBasis(vec(d1:d2),PsiR(itab)%V,iG,                  &
                                          BasisnD%para_SGType2,BasisnD%WeightSG(iG))
!        CALL tabR_AT_iG_TO_tabPackedBasis(OpPsi(itab)%RvecB,PsiR(itab)%V,iG,           &
!                                          BasisnD%para_SGType2,BasisnD%WeightSG(iG))
      ENDDO
    ENDDO

    !-collect result to root and boardcast----------------------------------------------
!    IF(MPI_np>1) THEN
!      CALL allocate_array(vec,1,size_RvecB*size_psi)
!
!      DO itab=1,size_psi
!        d1=size_RvecB*(itab-1)+1
!        d2=size_RvecB* itab
!        vec(d1:d2)=OpPsi(itab)%RvecB
!      ENDDO
!
!      CALL time_record(time_comm,time_temp1,time_temp2,1)
!      CALL MPI_Reduce_sum_Bcast(vec,size_RvecB*size_psi)
!      CALL time_record(time_comm,time_temp1,time_temp2,2)
!
!      DO itab=1,size_psi
!        d1=size_RvecB*(itab-1)+1
!        d2=size_RvecB* itab
!        OpPsi(itab)%RvecB=Vec(d1:d2)
!      ENDDO
!    ENDIF

    IF(MPI_np>1) THEN
      CALL time_record(time_comm,time_temp1,time_temp2,1)
      CALL MPI_Reduce_sum_Bcast(vec,size_RvecB*size_psi)
      CALL time_record(time_comm,time_temp1,time_temp2,2)
    ENDIF

    DO itab=1,size_psi
      d1=size_RvecB*(itab-1)+1
      d2=size_RvecB* itab
      OpPsi(itab)%RvecB=Vec(d1:d2)
    ENDDO

    IF(allocated(Vec)) deallocate(Vec)
    DO itab=1,size_psi
      CALL dealloc_TypeRVec(PsiR(itab))
    ENDDO

#endif
  END SUBROUTINE Action_MPI_S1
!======================================================================================= 

!=======================================================================================  
!> @brief action with MPI: scheme 2
!> @note imcrease MPI_mc to reduce memory usage  
!> @note set MPI_S2_L2=t to use more than one thread to distribute smolyak terms. 
!> require the preloading of mapping table for certain threads

!> @param Psi [in] |\psi>
!> @param OpPsi [inout] O|\psi>
!> @param BasisnD [in] base infortmation
!> @param para_Op [in] Operator
!======================================================================================= 
  SUBROUTINE Action_MPI_S2(Psi,OpPsi,BasisnD,para_Op)
    USE mod_system
    USE mod_nDindex
    USE mod_psi,                             ONLY:param_psi
    USE mod_SetOp,                           ONLY:param_Op
    USE mod_Coord_KEO,                       ONLY:CoordType
    USE mod_basis_set_alloc,                 ONLY:basis
    USE mod_OpPsi_SG4,                       ONLY:sub_TabOpPsi_OF_ONEDP_FOR_SGtype4
    USE mod_basis_BtoG_GtoB_SGType4,         ONLY:tabPackedBasis_TO_tabR_AT_iG,        &
                                                  tabR_AT_iG_TO_tabPackedBasis,        &
                                                  TypeRVec,dealloc_TypeRVec
    USE mod_basis_BtoG_GtoB_SGType4_MPI,     ONLY:set_iGs_MPI_mc
    USE mod_MPI_aux
    IMPLICIT NONE
    
    TYPE(param_psi),               intent(in)    :: Psi(:)
    TYPE(param_psi),               intent(inout) :: OpPsi(:)
    TYPE(basis),pointer,           intent(inout) :: BasisnD
    TYPE(param_Op),                intent(inout) :: para_Op
    
    TYPE(TypeRVec),allocatable                   :: PsiR(:)
    Real(kind=Rkind),allocatable                 :: Psi_ST(:)
    Real(kind=Rkind),allocatable                 :: RvecB_all(:)
    Real(kind=Rkind),allocatable                 :: Vec(:)
    Integer(kind=MPI_INTEGER_KIND)               :: Psi_ST_length(0:MPI_np-1)
    Integer(kind=MPI_INTEGER_KIND),pointer       :: size_ST(:) 
    Integer(kind=MPI_INTEGER_KIND),pointer       :: size_ST_mc(:,:) 
    Integer(kind=MPI_INTEGER_KIND),pointer       :: size_psi
    Integer                                      :: size_RvecB
    Integer                                      :: length_ST_mc
    Integer                                      :: ST_index1
    Integer                                      :: ST_index2
    Integer                                      :: ST_index_iG
    Integer                                      :: size_ST_iG
    Integer                                      :: iG
    Integer                                      :: itab
    Integer                                      :: time_auto1
    Integer                                      :: time_auto2
    Integer                                      :: d1
    Integer                                      :: d2
    Integer                                      :: ii
    Integer                                      :: i_mc
    Logical,pointer                              :: once_action
    Logical                                      :: iGs_change
 
#if(run_MPI)

    size_ST     => BasisnD%para_SGType2%size_ST
    size_ST_mc  => BasisnD%para_SGType2%size_ST_mc
    size_psi    => BasisnD%para_SGType2%size_psi
    once_action => BasisnD%para_SGType2%once_action

    IF(once_action) THEN
      write(out_unitp,*) 'action with MPI: Scheme 2'
      IF(MPI_S2_L2) THEN 
        write(out_unitp,*) 'enable 2-level distribution of Smolyak terms'
        IF(MPI_np<Num_L2) write(out_unitp,*) 'Warning: 2-level distribution requires more cores'
      ENDIF
      write(out_unitp,*) '------------------------------------------------'
      write(out_unitp,*) 'iGs_MPI initial: ',iGs_MPI
      write(out_unitp,*) '------------------------------------------------'
    ELSE
      IF(MPI_iGs_auto) THEN
        iGs_change=.FALSE.
        ! auto balance Smolyak terms assigned to different threads
        IF(MPI_np>1) CALL auto_iGs_MPI(para_Op,iGs_change)

        IF(iGs_change) THEN
          write(out_unitp,*) '------------------------------------------------'
          write(out_unitp,*) 'iGs_MPI changed, current: ',iGs_MPI
          write(out_unitp,*) '------------------------------------------------'
      
          ! get new size_ST according to new iGs_MPI
          CALL set_iGs_MPI_mc(BasisnD)
          CALL get_size_ST(BasisnD,size_ST,size_ST_mc)
        ENDIF
      ENDIF ! MPI_iGs_auto
    ENDIF

    ! initialize time count
    time_MPI_local=0
    time_MPI_calcu=0

    If(MPI_id==0) THEN
      size_psi=INT(size(Psi),MPI_INTEGER_KIND) ! for itab
      size_RvecB=size(psi(1)%RvecB)
    ENDIF
    CALL MPI_BCAST(size_psi,  size1_MPI,Int_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
    CALL MPI_BCAST(size_RvecB,size1_MPI,Int_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
    
    IF(allocated(PsiR)) deallocate(PsiR)
    allocate(PsiR(size_psi))

    !-----------------------------------------------------------------------------------
    IF(MPI_id==0) THEN
      CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,1)
      ! IF enable 2-level distributor
      !-1/Num_L2 cores will be used as second level distributor-------------------------
      ! n_level2 will >> MPI_np when MPI_S2_L2=f
      DO i_MPI=0,MPI_np-1,n_level2
        IF(i_MPI>0) THEN
          CALL allocate_array(RvecB_all,1,size_RvecB*size_psi)
          DO itab=1,size_psi
            d1=size_RvecB*(itab-1)+1
            d2=size_RvecB* itab
            RvecB_all(d1:d2)=psi(itab)%RvecB
          ENDDO
        
          CALL time_record(time_comm,time_temp1,time_temp2,1)
          CALL MPI_Send(RvecB_all,size_RvecB*size_psi,Real_MPI,i_MPI,                  &
                        i_MPI,MPI_COMM_WORLD,MPI_err)
          CALL time_record(time_comm,time_temp1,time_temp2,2)
        ENDIF
      ENDDO ! i_MPI=1,n_level2
      IF(allocated(RvecB_all)) deallocate(RvecB_all)
      CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,2)

      !-distributing smolyak terms from root threads------------------------------------
      DO i_mc=1,MPI_mc
        CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,1)
        DO i_MPI=1,MIN(MPI_np-1,n_level2-1)
          length_ST_mc=size_psi*size_ST_mc(i_mc,i_MPI)
          CALL allocate_array(Psi_ST,1,length_ST_mc)
          ST_index1=0
          DO iG=iGs_MPI_mc(1,i_mc,i_MPI),iGs_MPI_mc(2,i_mc,i_MPI)
            DO itab=1,size_psi
              CALL tabPackedBasis_TO_tabR_AT_iG(PsiR(itab)%V,psi(itab)%RvecB,iG,       &
                                                BasisnD%para_SGType2)
              ST_index2=ST_index1+size(PsiR(itab)%V)
              Psi_ST(ST_index1+1:ST_index2)=PsiR(itab)%V
              ST_index1=ST_index2
            ENDDO
          ENDDO ! for iG
          Psi_ST_length(i_MPI)=ST_index1

          IF(once_action .AND. i_mc==1) THEN 
            write(out_unitp,*) 'length check:',i_MPI,size_psi,size_RvecB,              &
                                size_ST(i_MPI),size_ST_mc(i_mc,i_MPI),sum(size_ST),    &
                                sum(BasisnD%para_SGType2%tab_nb_OF_SRep(:))
          ENDIF

          ! double check
          IF(Abs(length_ST_mc-Psi_ST_length(i_MPI))>0) THEN
            write(out_unitp,*) 'ST length check:',length_ST_mc,Psi_ST_length(i_MPI)
            STOP 'error in MPI action part,check length'
          ENDIF

          CALL time_record(time_comm,time_temp1,time_temp2,1)
          CALL MPI_Send(Psi_ST,length_ST_mc,Real_MPI,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
          CALL time_record(time_comm,time_temp1,time_temp2,2)
        ENDDO ! i_MPI=1,MIN(MPI_np-1,n_level2-1)
        CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,2)

        IF(allocated(Psi_ST)) deallocate(Psi_ST)

        !-calculations on MPI_id=0------------------------------------------------------
        CALL time_record(time_MPI_calcu(MPI_id),time_auto1,time_auto2,1)
        DO iG=iGs_MPI_mc(1,i_mc,MPI_id),iGs_MPI_mc(2,i_mc,MPI_id)
          !-transfore to SRep-----------------------------------------------------------
          DO itab=1,size_psi
            CALL tabPackedBasis_TO_tabR_AT_iG(PsiR(itab)%V,psi(itab)%RvecB,iG,         &
                                              BasisnD%para_SGType2)
          ENDDO

          !-main calculation------------------------------------------------------------
          CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4(PsiR,iG,                              &
                          BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval(:,iG),para_Op)

          !-back to compact form--------------------------------------------------------
          DO itab=1,size_psi
            CALL tabR_AT_iG_TO_tabPackedBasis(OpPsi(itab)%RvecB,PsiR(itab)%V,iG,       &
                                              BasisnD%para_SGType2,BasisnD%WeightSG(iG))
          ENDDO
        ENDDO ! main loop of iG for calcuation on root
        CALL time_record(time_MPI_calcu(MPI_id),time_auto1,time_auto2,2)

        !-receive results from other threads----------------------------------------------
        CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,1)
        DO i_MPI=1,MPI_np-1
          IF(mod(i_MPI,n_level2)==0) THEN !-level2 threads--------------------------------
            CALL allocate_array(Vec,1,size_RvecB*size_psi)
            CALL time_record(time_comm,time_temp1,time_temp2,1)
            CALL MPI_Recv(Vec,size_RvecB*size_psi,Real_MPI,i_MPI,                        &
                          i_MPI,MPI_COMM_WORLD,MPI_stat,MPI_err)
            CALL time_record(time_comm,time_temp1,time_temp2,2)

            DO itab=1,size_psi
              d1=size_RvecB*(itab-1)+1
              d2=size_RvecB* itab
              OpPsi(itab)%RvecB=OpPsi(itab)%RvecB+Vec(d1:d2)
            ENDDO
            IF(allocated(Vec)) deallocate(Vec)
          ELSE !-normal threads---------------------------------------------------------
            length_ST_mc=size_psi*size_ST_mc(i_mc,i_MPI)
            CALL allocate_array(Psi_ST,1,length_ST_mc)

            CALL time_record(time_comm,time_temp1,time_temp2,1)
            CALL MPI_Recv(Psi_ST,length_ST_mc,Real_MPI,i_MPI,i_MPI,                    &
                          MPI_COMM_WORLD,MPI_stat,MPI_err)
            CALL time_record(time_comm,time_temp1,time_temp2,2)

            ST_index1=0
            DO iG=iGs_MPI_mc(1,i_mc,i_MPI),iGs_MPI_mc(2,i_mc,i_MPI)
              size_ST_iG=BasisnD%para_SGType2%tab_nb_OF_SRep(iG)                       &
                            *BasisnD%para_SGType2%nb0
              DO itab=1,size_psi
                d1=ST_index1+1
                d2=ST_index1+size_ST_iG
                CALL tabR_AT_iG_TO_tabPackedBasis(OpPsi(itab)%RvecB,Psi_ST(d1:d2),iG,  &
                                              BasisnD%para_SGType2,BasisnD%WeightSG(iG))
                ST_index1=d2
              ENDDO
            ENDDO ! iG
            IF(allocated(Psi_ST)) deallocate(Psi_ST)
          ENDIF ! mod(i_MPI,n_level2)==0
        ENDDO ! for i_MPI=1,MPI_np-1  
        CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,2)
      ENDDO ! i_mc
    ENDIF ! for MPI_id==0

    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    IF(MPI_id/=0) THEN
      !-level 2 threads-----------------------------------------------------------------
      IF(mod(MPI_id,n_level2)==0) THEN
        CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,1)
        ! receive psi(itab)%RvecB from root threads
        CALL allocate_array(RvecB_all,1,size_RvecB*size_psi)
        CALL MPI_Recv(RvecB_all,size_RvecB*size_psi,Real_MPI,root_MPI,                 &
                      MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)

        !-send Smolyak terms to the other threads---------------------------------------
        DO i_mc=1,MPI_mc
          DO i_MPI=MPI_id+1,MIN(MPI_np-1,MPI_id+n_level2-1)
            length_ST_mc=size_psi*size_ST_mc(i_mc,i_MPI)
            CALL allocate_array(Psi_ST,1,length_ST_mc)
            ST_index1=0
            DO iG=iGs_MPI_mc(1,i_mc,i_MPI),iGs_MPI_mc(2,i_mc,i_MPI)
              DO itab=1,size_psi
                d1=size_RvecB*(itab-1)+1
                d2=size_RvecB* itab
                CALL tabPackedBasis_TO_tabR_AT_iG(PsiR(itab)%V,RvecB_all(d1:d2),iG,      &
                                                  BasisnD%para_SGType2)
                ST_index2=ST_index1+size(PsiR(itab)%V)
                Psi_ST(ST_index1+1:ST_index2)=PsiR(itab)%V
                ST_index1=ST_index2
              ENDDO
            ENDDO ! iG
            Psi_ST_length(i_MPI)=ST_index1

            ! double check
            IF(Abs(length_ST_mc-Psi_ST_length(i_MPI))>0) THEN
              write(out_unitp,*) 'error in MPI action part,check Psi_ST_length'
              STOP
            ENDIF

            CALL time_record(time_comm,time_temp1,time_temp2,1)
            CALL MPI_Send(Psi_ST,length_ST_mc,Real_MPI,i_MPI,i_MPI,                    &
                          MPI_COMM_WORLD,MPI_err)
            CALL time_record(time_comm,time_temp1,time_temp2,2)
          ENDDO ! i_MPI=MPI_id+1,MIN(MPI_np-1,MPI_id+n_level2-1)
          CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,2)

          IF(allocated(Psi_ST)) deallocate(Psi_ST)

          !-calculation-------------------------------------------------------------------
          CALL time_record(time_MPI_calcu(MPI_id),time_auto1,time_auto2,1)
          CALL allocate_array(Vec,1,size_RvecB*size_psi)
          DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)
            DO itab=1,size_psi
              d1=size_RvecB*(itab-1)+1
              d2=size_RvecB* itab
              CALL tabPackedBasis_TO_tabR_AT_iG(PsiR(itab)%V,RvecB_all(d1:d2),iG,      &
                                                BasisnD%para_SGType2)
            ENDDO

            CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4(PsiR,iG,                            &
                          BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval(:,iG),para_Op)

            DO itab=1,size_psi
              d1=size_RvecB*(itab-1)+1
              d2=size_RvecB* itab
              CALL tabR_AT_iG_TO_tabPackedBasis(Vec(d1:d2),PsiR(itab)%V,iG,            &
                                              BasisnD%para_SGType2,BasisnD%WeightSG(iG))
            ENDDO
          ENDDO
          CALL time_record(time_MPI_calcu(MPI_id),time_auto1,time_auto2,2)

          CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,1)
          CALL MPI_Send(Vec,size_RvecB*size_psi,Real_MPI,root_MPI,                     &
                        MPI_id,MPI_COMM_WORLD,MPI_err)
          CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,2)
          IF(allocated(RvecB_all)) deallocate(RvecB_all)
        ENDDO ! i_mc=1,MPI_mc
        
      ELSE !-other threads--------------------------------------------------------------

        DO i_mc=1,MPI_mc
          !-get PsiR(itab)%V from root or second level threads--------------------------
          length_ST_mc=size_psi*size_ST_mc(i_mc,MPI_id)
          CALL allocate_array(Psi_ST,1,length_ST_mc)

          CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,1)
          CALL MPI_Recv(Psi_ST,length_ST_mc,Real_MPI,MPI_id/n_level2*n_level2,MPI_id,  &
                        MPI_COMM_WORLD,MPI_stat,MPI_err)
          CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,2)
          
          !-loop for main calculation---------------------------------------------------
          CALL time_record(time_MPI_calcu(MPI_id),time_auto1,time_auto2,1)
          ST_index1=0
          ST_index_iG=0
          DO iG=iGs_MPI_mc(1,i_mc,MPI_id),iGs_MPI_mc(2,i_mc,MPI_id)
            size_ST_iG=BasisnD%para_SGType2%tab_nb_OF_SRep(iG)*BasisnD%para_SGType2%nb0

            !-extract SRep from PsiR_temp-------------------------------------------------
            DO itab=1,size_psi
              CALL allocate_array(PsiR(itab)%V,1,size_ST_iG)
              ST_index2=ST_index1+size_ST_iG
              PsiR(itab)%V=Psi_ST(ST_index1+1:ST_index2)
              ST_index1=ST_index2
            ENDDO

            !-main calculation----------------------------------------------------------
            CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4(PsiR,iG,                            &
                          BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval(:,iG),para_Op)

            !-pack and PsiR(itab)%V-----------------------------------------------------
            Do itab=1,size_psi
              d1=ST_index_iG+1
              d2=ST_index_iG+size(PsiR(itab)%V)
              Psi_ST(d1:d2)=PsiR(itab)%V
              ST_index_iG=d2
            ENDDO
          ENDDO ! for iG
          CALL time_record(time_MPI_calcu(MPI_id),time_auto1,time_auto2,2)

          !-send back PsiR(itab)%V------------------------------------------------------
          CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,1)
          CALL MPI_Send(Psi_ST,length_ST_mc,Real_MPI,root_MPI,MPI_id,                  &
                        MPI_COMM_WORLD,MPI_err)
          CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,2)
        ENDDO ! i_mc=1,MPI_mc
      ENDIF
    ENDIF ! for MPI_id/=0

    IF(allocated(Vec)) deallocate(Vec)
    IF(allocated(Psi_ST)) deallocate(Psi_ST)
    IF(allocated(RvecB_all)) deallocate(RvecB_all)
    DO itab=1,size_psi
      CALL dealloc_TypeRVec(PsiR(itab))
    END DO

#endif
  END SUBROUTINE Action_MPI_S2
!======================================================================================= 

!=======================================================================================  
!> @brief action with MPI: scheme 3
!> @note works on more than one nodes.
!> @note imcrease MPI_mc to reduce memory usage  
!> require the preloading of mapping table for certain threads
!
!> @param Psi [in] |\psi>
!> @param OpPsi [inout] O|\psi>
!> @param BasisnD [in] base infortmation
!> @param para_Op [in] Operator
! consider to merge with scheme 2
!======================================================================================= 
  SUBROUTINE Action_MPI_S3(Psi,OpPsi,BasisnD,para_Op)
    USE mod_system
    USE mod_nDindex
    USE mod_psi,                             ONLY:param_psi
    USE mod_SetOp,                           ONLY:param_Op
    USE mod_Coord_KEO,                       ONLY:CoordType
    USE mod_basis_set_alloc,                 ONLY:basis
    USE mod_OpPsi_SG4,                       ONLY:sub_TabOpPsi_OF_ONEDP_FOR_SGtype4
    USE mod_basis_BtoG_GtoB_SGType4,         ONLY:tabPackedBasis_TO_tabR_AT_iG,        &
                                                  tabR_AT_iG_TO_tabPackedBasis,        &
                                                  TypeRVec,dealloc_TypeRVec
    USE mod_basis_BtoG_GtoB_SGType4_MPI,     ONLY:set_iGs_MPI_mc
    USE mod_MPI_aux
    IMPLICIT NONE
    
    TYPE(param_psi),               intent(in)    :: Psi(:)
    TYPE(param_psi),               intent(inout) :: OpPsi(:)
    TYPE(basis),pointer,           intent(inout) :: BasisnD
    TYPE(param_Op),                intent(inout) :: para_Op
    
    TYPE(TypeRVec),allocatable                   :: PsiR(:)
    Real(kind=Rkind),allocatable                 :: Psi_ST(:)
    Real(kind=Rkind),allocatable                 :: RvecB_all(:)
    Real(kind=Rkind),allocatable                 :: Vec(:)
    Integer(kind=MPI_INTEGER_KIND)               :: Psi_ST_length(0:MPI_np-1)
    Integer(kind=MPI_INTEGER_KIND),pointer       :: size_ST(:) 
    Integer(kind=MPI_INTEGER_KIND),pointer       :: size_ST_mc(:,:) 
    Integer(kind=MPI_INTEGER_KIND),pointer       :: size_psi
    Integer                                      :: size_RvecB
    Integer                                      :: length_ST_mc
    Integer                                      :: ST_index1
    Integer                                      :: ST_index2
    Integer                                      :: ST_index_iG
    Integer                                      :: size_ST_iG
    Integer                                      :: iG
    Integer                                      :: itab
    Integer                                      :: time_auto1
    Integer                                      :: time_auto2
    Integer                                      :: d1
    Integer                                      :: d2
    Integer                                      :: ii
    Integer                                      :: i_mc
    Logical,pointer                              :: once_action
    Logical                                      :: iGs_change

#if(run_MPI)

    size_ST     => BasisnD%para_SGType2%size_ST
    size_ST_mc  => BasisnD%para_SGType2%size_ST_mc
    size_psi    => BasisnD%para_SGType2%size_psi
    once_action => BasisnD%para_SGType2%once_action

    IF(once_action) THEN
      write(out_unitp,*) 'action with MPI: Scheme 3'
      write(out_unitp,*) '------------------------------------------------'
      write(out_unitp,*) 'iGs_MPI initial: ',iGs_MPI
      write(out_unitp,*) '------------------------------------------------'
    ELSE
      IF(MPI_iGs_auto) THEN
        iGs_change=.FALSE.
        ! auto balance Smolyak terms assigned to different threads
        IF(MPI_np>1) CALL auto_iGs_MPI(para_Op,iGs_change)

        IF(iGs_change) THEN
          write(out_unitp,*) '------------------------------------------------'
          write(out_unitp,*) 'iGs_MPI changed, current: ',iGs_MPI
          write(out_unitp,*) '------------------------------------------------'
      
          ! get new size_ST according to new iGs_MPI
          CALL set_iGs_MPI_mc(BasisnD)
          CALL get_size_ST(BasisnD,size_ST,size_ST_mc)
        ENDIF
      ENDIF ! MPI_iGs_auto
    ENDIF

    ! initialize time count
    time_MPI_local=0
    time_MPI_calcu=0

    If(MPI_id==0) THEN
      size_psi=INT(size(Psi),MPI_INTEGER_KIND) ! for itab
      size_RvecB=size(psi(1)%RvecB)
    ENDIF
    CALL MPI_BCAST(size_psi,  size1_MPI,Int_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
    CALL MPI_BCAST(size_RvecB,size1_MPI,Int_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)

    IF(allocated(PsiR)) deallocate(PsiR)
    allocate(PsiR(size_psi))

    IF(MPI_id==0 .OR. MPI_nodes_p0) THEN

      DO i_mc=1,MPI_mc
        CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,1)
        DO i_MPI=MPI_sub_id(1),MPI_sub_id(2)
          length_ST_mc=size_psi*size_ST_mc(i_mc,i_MPI)
          CALL allocate_array(Psi_ST,1,length_ST_mc)
          ST_index1=0
          DO iG=iGs_MPI_mc(1,i_mc,i_MPI),iGs_MPI_mc(2,i_mc,i_MPI)
            DO itab=1,size_psi
              CALL tabPackedBasis_TO_tabR_AT_iG(PsiR(itab)%V,psi(itab)%RvecB,iG,       &
                                                BasisnD%para_SGType2)
              ST_index2=ST_index1+size(PsiR(itab)%V)
              Psi_ST(ST_index1+1:ST_index2)=PsiR(itab)%V
              ST_index1=ST_index2
            ENDDO
          ENDDO ! for iG
          Psi_ST_length(i_MPI)=ST_index1

          IF(once_action .AND. i_mc==1) THEN 
            write(out_unitp,*) 'length check:',i_MPI,size_psi,size_RvecB,              &
                                size_ST(i_MPI),size_ST_mc(i_mc,i_MPI),sum(size_ST),    &
                                sum(BasisnD%para_SGType2%tab_nb_OF_SRep(:))
          ENDIF

          ! double check
          IF(Abs(length_ST_mc-Psi_ST_length(i_MPI))>0) THEN
            write(out_unitp,*) 'ST length check:',length_ST_mc,Psi_ST_length(i_MPI)
            STOP 'error in MPI action part,check length'
          ENDIF

          CALL time_record(time_comm,time_temp1,time_temp2,1)
          CALL MPI_Send(Psi_ST,length_ST_mc,Real_MPI,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
          CALL time_record(time_comm,time_temp1,time_temp2,2)
        ENDDO ! i_MPI
        CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,2)
        IF(allocated(Psi_ST)) deallocate(Psi_ST)

        !-calculations on MPI_node_p0---------------------------------------------------
        CALL time_record(time_MPI_calcu(MPI_id),time_auto1,time_auto2,1)
        DO iG=iGs_MPI_mc(1,i_mc,MPI_id),iGs_MPI_mc(2,i_mc,MPI_id)
          !-transfore to SRep-----------------------------------------------------------
          DO itab=1,size_psi
            CALL tabPackedBasis_TO_tabR_AT_iG(PsiR(itab)%V,psi(itab)%RvecB,iG,         &
                                              BasisnD%para_SGType2)
          ENDDO

          !-main calculation------------------------------------------------------------
          CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4(PsiR,iG,                              &
                          BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval(:,iG),para_Op)

          !-back to compact form--------------------------------------------------------
          DO itab=1,size_psi
            CALL tabR_AT_iG_TO_tabPackedBasis(OpPsi(itab)%RvecB,PsiR(itab)%V,iG,       &
                                              BasisnD%para_SGType2,BasisnD%WeightSG(iG))
          ENDDO
        ENDDO ! main loop of iG for calcuation on root
        CALL time_record(time_MPI_calcu(MPI_id),time_auto1,time_auto2,2)

        !-receive results from other threads--------------------------------------------
        CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,1)
        DO i_MPI=MPI_sub_id(1),MPI_sub_id(2)
          length_ST_mc=size_psi*size_ST_mc(i_mc,i_MPI)
          CALL allocate_array(Psi_ST,1,length_ST_mc)

          CALL time_record(time_comm,time_temp1,time_temp2,1)
          CALL MPI_Recv(Psi_ST,length_ST_mc,Real_MPI,i_MPI,i_MPI,                      &
                        MPI_COMM_WORLD,MPI_stat,MPI_err)
          CALL time_record(time_comm,time_temp1,time_temp2,2)

          ST_index1=0
          DO iG=iGs_MPI_mc(1,i_mc,i_MPI),iGs_MPI_mc(2,i_mc,i_MPI)
            size_ST_iG=BasisnD%para_SGType2%tab_nb_OF_SRep(iG)                         &
                      *BasisnD%para_SGType2%nb0
            DO itab=1,size_psi
              d1=ST_index1+1
              d2=ST_index1+size_ST_iG
              CALL tabR_AT_iG_TO_tabPackedBasis(OpPsi(itab)%RvecB,Psi_ST(d1:d2),iG,    &
                                              BasisnD%para_SGType2,BasisnD%WeightSG(iG))
              ST_index1=d2
            ENDDO
          ENDDO ! iG
          IF(allocated(Psi_ST)) deallocate(Psi_ST)
        ENDDO ! i_MPI
      ENDDO ! i_mc

      ! collect result from each node-------------------------------------------------
      IF(MPI_id==0) THEN
        IF(MPI_nodes_num>1) THEN
          CALL allocate_array(RvecB_all,1,size_RvecB*size_psi)
          DO ii=1,MPI_nodes_num-1
            CALL MPI_Recv(RvecB_all,size_RvecB*size_psi,Real_MPI,MPI_nodes_p00(ii),    &
                          ii,MPI_COMM_WORLD,MPI_stat,MPI_err)

            DO itab=1,size_psi
              d1=size_RvecB*(itab-1)+1
              d2=size_RvecB* itab
              OpPsi(itab)%RvecB=OpPsi(itab)%RvecB+RvecB_all(d1:d2)
            ENDDO
          ENDDO

          ! share result with node_p0
          DO itab=1,size_psi
            d1=size_RvecB*(itab-1)+1
            d2=size_RvecB* itab
            RvecB_all(d1:d2)=OpPsi(itab)%RvecB
          ENDDO

          DO ii=1,MPI_nodes_num-1
            CALL MPI_Send(RvecB_all,size_RvecB*size_psi,Real_MPI,MPI_nodes_p00(ii),    &
                          MPI_nodes_p00(ii),MPI_COMM_WORLD,MPI_err)
          ENDDO
        ENDIF

      ELSEIF(MPI_nodes_p0) THEN

        CALL allocate_array(RvecB_all,1,size_RvecB*size_psi)
        DO itab=1,size_psi
          d1=size_RvecB*(itab-1)+1
          d2=size_RvecB* itab
          RvecB_all(d1:d2)=OpPsi(itab)%RvecB
        ENDDO

        CALL MPI_Send(RvecB_all,size_RvecB*size_psi,Real_MPI,root_MPI,MPI_node_id,   &
                      MPI_COMM_WORLD,MPI_err)
        CALL MPI_Recv(RvecB_all,size_RvecB*size_psi,Real_MPI,root_MPI,MPI_id,        &
                      MPI_COMM_WORLD,MPI_stat,MPI_err)
        DO itab=1,size_psi
          d1=size_RvecB*(itab-1)+1
          d2=size_RvecB* itab
          OpPsi(itab)%RvecB=RvecB_all(d1:d2)
        ENDDO

        IF(allocated(RvecB_all)) deallocate(RvecB_all)
      ENDIF
      CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,2)
    ENDIF ! MPI_id==0 .OR. MPI_nodes_p0

    !-----------------------------------------------------------------------------------
    !-----------------------------------------------------------------------------------
    IF(.NOT. (MPI_id==0 .OR. MPI_nodes_p0)) THEN
      DO i_mc=1,MPI_mc
        length_ST_mc=size_psi*size_ST_mc(i_mc,MPI_id)
        CALL allocate_array(Psi_ST,1,length_ST_mc)

        CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,1)
        CALL MPI_Recv(Psi_ST,length_ST_mc,Real_MPI,MPI_node_p0_id,MPI_id,              &
                      MPI_COMM_WORLD,MPI_stat,MPI_err)
        CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,2)

        !-loop for main calculation-----------------------------------------------------
        CALL time_record(time_MPI_calcu(MPI_id),time_auto1,time_auto2,1)
        ST_index1=0
        ST_index_iG=0
        DO iG=iGs_MPI_mc(1,i_mc,MPI_id),iGs_MPI_mc(2,i_mc,MPI_id)
          size_ST_iG=BasisnD%para_SGType2%tab_nb_OF_SRep(iG)*BasisnD%para_SGType2%nb0

          !-extract SRep from PsiR_temp-------------------------------------------------
          DO itab=1,size_psi
            CALL allocate_array(PsiR(itab)%V,1,size_ST_iG)
            ST_index2=ST_index1+size_ST_iG
            PsiR(itab)%V=Psi_ST(ST_index1+1:ST_index2)
            ST_index1=ST_index2
          ENDDO

          !-main calculation----------------------------------------------------------
          CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4(PsiR,iG,                            &
                        BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval(:,iG),para_Op)

          !-pack and PsiR(itab)%V-------------------------------------------------------
          Do itab=1,size_psi
            d1=ST_index_iG+1
            d2=ST_index_iG+size(PsiR(itab)%V)
            Psi_ST(d1:d2)=PsiR(itab)%V
            ST_index_iG=d2
          ENDDO
        ENDDO ! for iG
        CALL time_record(time_MPI_calcu(MPI_id),time_auto1,time_auto2,2)

        !-send back PsiR(itab)%V--------------------------------------------------------
        CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,1)
        CALL MPI_Send(Psi_ST,length_ST_mc,Real_MPI,MPI_node_p0_id,MPI_id,              &
                      MPI_COMM_WORLD,MPI_err)
        CALL time_record(time_MPI_local(MPI_id),time_auto1,time_auto2,2)
      ENDDO ! i_mc
    ENDIF ! .NOT. (MPI_id==0 .OR. MPI_nodes_p0)

#endif
  END SUBROUTINE Action_MPI_S3

!=======================================================================================  
!> @brief action with MPI: scheme 4
!> compress the terms to be send. 
!> not very effeicient, consider to remove it later
!=======================================================================================  
  SUBROUTINE Action_MPI_S4(Psi,OpPsi,BasisnD,para_Op)
    USE mod_system
    USE mod_nDindex
    USE mod_Coord_KEO,                  ONLY:CoordType
    USE mod_basis_set_alloc,            ONLY:basis
    USE mod_OpPsi_SG4,                  ONLY:sub_TabOpPsi_OF_ONEDP_FOR_SGtype4
    USE mod_basis_BtoG_GtoB_SGType4,    ONLY:tabPackedBasis_TO_tabR_AT_iG,             &
                                             tabR_AT_iG_TO_tabPackedBasis,             &
                                             TypeRVec,dealloc_TypeRVec
    USE mod_basis_BtoG_GtoB_SGType4_MPI,ONLY:PackedBasis_TO_tabR_index_MPI,            &
                                             tabR_TO_tabPackedBasis_MPI,               &
                                             tabPackedBasis_TO_tabR_MPI
    USE mod_psi,                        ONLY:param_psi
    USE mod_SetOp,                      ONLY:param_Op
    USE mod_MPI_aux
    IMPLICIT NONE
    
    TYPE(param_psi),                       intent(in)    :: Psi(:)
    TYPE(param_psi),                       intent(inout) :: OpPsi(:)
    TYPE(basis),pointer,                   intent(inout) :: BasisnD
    TYPE(param_Op),                        intent(inout) :: para_Op

    TYPE(TypeRVec),allocatable                           :: PsiR(:)
    Real(kind=Rkind),allocatable                         :: all_RvecB_temp(:)
    Real(kind=Rkind),allocatable                         :: all_RvecB_temp2(:)
    Integer(kind=MPI_INTEGER_KIND),pointer               :: size_psi
    Integer(kind=MPI_INTEGER_KIND),pointer               :: size_ST(:) 
    Integer(kind=MPI_INTEGER_KIND),pointer               :: reduce_Vlength_master(:)
    Integer(kind=MPI_INTEGER_KIND),pointer               :: reduce_Vlength
    Integer,pointer                                      :: Max_nDI_ib0
    Integer,pointer                                      :: V_allcount
    Integer,pointer                                      :: V_allcount2
    Integer                                              :: iG
    Integer                                              :: ii
    Integer                                              :: itab
    Logical,pointer                                      :: once_action

#if(run_MPI)
    
    size_psi => BasisnD%para_SGType2%size_psi
    size_ST  => BasisnD%para_SGType2%size_ST
    
    If(MPI_id==0) size_psi=INT(size(Psi),MPI_INTEGER_KIND) ! for itab
    CALL MPI_BCAST(size_psi,size1_MPI,Int_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)

    Max_nDI_ib0    => BasisnD%para_SGType2%Max_nDI_ib0
    reduce_Vlength => BasisnD%para_SGType2%reduce_Vlength
    V_allcount     => BasisnD%para_SGType2%V_allcount
    V_allcount2    => BasisnD%para_SGType2%V_allcount2
    once_action    => BasisnD%para_SGType2%once_action
    
    ! calculate total length of vectors for each threads--------------------------------
    IF(once_action .AND. MPI_id==0) THEN
      write(out_unitp,*) 'action with MPI: Scheme 4'
      allocate(BasisnD%para_SGType2%nDI_index_master(0:MPI_np-1))
      allocate(BasisnD%para_SGType2%reduce_Vlength_master(0:MPI_np-1))
    ENDIF
    reduce_Vlength_master=>BasisnD%para_SGType2%reduce_Vlength_master

    If(MPI_id==0) Max_nDI_ib0=size(psi(1)%RvecB)/BasisnD%para_SGType2%nb0
    CALL MPI_BCAST(Max_nDI_ib0,size1_MPI,Int_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
         
    ! clean PsiR     
    IF(allocated(PsiR)) deallocate(PsiR)
    allocate(PsiR(size_psi))
    
    ! works on the other threads
    !-----------------------------------------------------------------------------------
    ! all information for reduce_Vlength_MPI is keeped on master
    ! each thread keep its own    
    IF(MPI_id/=0) THEN
      !-generate index,  once only------------------------------------------------------
      IF(once_action) THEN
        ! initialize the size for index for pack psi on each threads
        BasisnD%para_SGType2%num_nDI_index=size_ST(MPI_id)/20
        reduce_Vlength=0
        V_allcount=0
        CALL allocate_array(BasisnD%para_SGType2%nDI_index,                            &
                            1,BasisnD%para_SGType2%num_nDI_index)
        CALL allocate_array(BasisnD%para_SGType2%nDI_index_list,1,size_ST(MPI_id))

        DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)
          ! note size(psi(i)%RvecB) same for all i 
          CALL PackedBasis_TO_tabR_index_MPI(iG,BasisnD%para_SGType2,reduce_Vlength,   &
                                         BasisnD%para_SGType2%nDI_index,Max_nDI_ib0,   &
                                         BasisnD%para_SGType2%nDI_index_list)
        ENDDO
        write(out_unitp,*) 'V_allcount check',V_allcount,size_ST(MPI_id),          &
                                              reduce_Vlength,' from ',MPI_id
        CALL MPI_Send(reduce_Vlength,size1_MPI,Int_MPI,root_MPI,MPI_id,                &
                      MPI_COMM_WORLD,MPI_err)
        CALL MPI_Send(BasisnD%para_SGType2%nDI_index(1:reduce_Vlength),                &
                      reduce_Vlength,Int_MPI,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_err)
      ENDIF
      
      !> check length of integer
      IF(Int(reduce_Vlength,8)*Int(size_psi,8)>huge(0_4)                               &
         .AND. MPI_INTEGER_KIND==4) THEN
        STOP 'integer exceed 32-bit MPI, use 64-bit MPI instead'
      ENDIF
      
      !-wait for master-----------------------------------------------------------------
      CALL allocate_array(all_RvecB_temp,1,reduce_Vlength*size_psi)
      CALL allocate_array(all_RvecB_temp2,1,reduce_Vlength*size_psi)
      Call MPI_Recv(all_RvecB_temp,reduce_Vlength*size_psi,MPI_REAL8,root_MPI,         &
                    MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)

      !-calculation on threads----------------------------------------------------------
      V_allcount=0
      V_allcount2=0   

      DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)
        ! all_RvecB_temp --> PsiR
        CALL tabPackedBasis_TO_tabR_MPI(PsiR,all_RvecB_temp,iG,BasisnD%para_SGType2,   &
                        BasisnD%para_SGType2%nDI_index,reduce_Vlength,size_psi,        &
                        Max_nDI_ib0,BasisnD%para_SGType2%nDI_index_list)

        !-main calculation--------------------------------------------------------------
        CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4(PsiR,iG,                                &
                        BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval(:,iG),para_Op)

        !-pack PsiR in the slave threads------------------------------------------------
        CALL tabR_TO_tabPackedBasis_MPI(all_RvecB_temp2,PsiR,iG,                       &
                         BasisnD%para_SGType2,BasisnD%WeightSG(iG),                    &
                         BasisnD%para_SGType2%nDI_index,reduce_Vlength,size_psi,       &
                         Max_nDI_ib0,BasisnD%para_SGType2%nDI_index_list)
      ENDDO 
      
      !-send packed PsiR to master------------------------------------------------------
      CALL MPI_Send(all_RvecB_temp2,reduce_Vlength*size_psi,MPI_REAL8,root_MPI,        &
                    MPI_id,MPI_COMM_WORLD,MPI_err)
                   
    ENDIF ! for MPI_id/=0  

    ! works on the master threads
    !-----------------------------------------------------------------------------------
    ! save index information for all threads on master & send vectors to other threads
    IF(MPI_id==0) THEN
      ! works for the other threads-----------------------------------------------------
      Do i_MPI=1,MPI_np-1
        !-wait for other threads--------------------------------------------------------
        CALL time_record(time_comm,time_temp1,time_temp2,1)
        IF(once_action) THEN
          CALL MPI_Recv(reduce_Vlength,size1_MPI,Int_MPI,i_MPI,i_MPI,                  &
                        MPI_COMM_WORLD,MPI_stat,MPI_err)
          reduce_Vlength_master(i_MPI)=reduce_Vlength
          allocate(BasisnD%para_SGType2%nDI_index_master(i_MPI)%array(reduce_Vlength))
          CALL MPI_Recv(BasisnD%para_SGType2%nDI_index_master(i_MPI)%array,            &
                        reduce_Vlength,Int_MPI,i_MPI,i_MPI,MPI_COMM_WORLD,             &
                        MPI_stat,MPI_err)
          write(out_unitp,*) 'length of comm list:',reduce_Vlength_master(i_MPI),      &
                             'from',i_MPI
        ENDIF ! for once_action
        CALL time_record(time_comm,time_temp1,time_temp2,2)

        ! pack vectores to send
        CALL allocate_array(all_RvecB_temp,1,reduce_Vlength_master(i_MPI)*size_psi)
        DO itab=1,size_psi
          DO ii=1,reduce_Vlength_master(i_MPI)
            temp_int=(itab-1)*reduce_Vlength_master(i_MPI)+ii
            all_RvecB_temp(temp_int)=psi(itab)%RvecB(                                  &
                                 BasisnD%para_SGType2%nDI_index_master(i_MPI)%array(ii))
          ENDDO
        ENDDO

        CALL time_record(time_comm,time_temp1,time_temp2,1)
        CALL MPI_Send(all_RvecB_temp,size_psi*reduce_Vlength_master(i_MPI),            &
                      MPI_REAL8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
        CALL time_record(time_comm,time_temp1,time_temp2,2)

      ENDDO ! for i_MPI=1,MPI_np-1

      !-calculation on master-----------------------------------------------------------
      DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)
        DO itab=1,size_psi
          CALL tabPackedBasis_TO_tabR_AT_iG(PsiR(itab)%V,psi(itab)%RvecB,              &
                                            iG,BasisnD%para_SGType2)
        ENDDO 

        !-main calculation--------------------------------------------------------------
        CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4(PsiR,iG,                                &
                        BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval(:,iG),para_Op) 
              
        DO itab=1,size_psi
          CALL tabR_AT_iG_TO_tabPackedBasis(OpPsi(itab)%RvecB,PsiR(itab)%V,iG,         &
                                            BasisnD%para_SGType2,BasisnD%WeightSG(iG))
        ENDDO 

      ENDDO ! for iG

      !-receive results from slave threads----------------------------------------------
      Do i_MPI=1,MPI_np-1
        CALL allocate_array(all_RvecB_temp,1,size_psi*reduce_Vlength_master(i_MPI))
        CALL time_record(time_comm,time_temp1,time_temp2,1)
        CALL MPI_Recv(all_RvecB_temp,size_psi*reduce_Vlength_master(i_MPI),Real_MPI,   &
                      i_MPI,i_MPI,MPI_COMM_WORLD,MPI_stat,MPI_err)
        CALL time_record(time_comm,time_temp1,time_temp2,2)

        !-extract results from other threads--------------------------------------------
        DO itab=1,size_psi
          Do ii=1,reduce_Vlength_master(i_MPI)
            temp_int=BasisnD%para_SGType2%nDI_index_master(i_MPI)%array(ii)
            OpPsi(itab)%RvecB(temp_int)=OpPsi(itab)%RvecB(temp_int)                    &
                               +all_RvecB_temp((itab-1)*reduce_Vlength_master(i_MPI)+ii)
          ENDDO 
        ENDDO
      ENDDO ! for i_MPI=1,MPI_np-1
    ENDIF ! for MPI_id==0
                  
    IF(allocated(all_RvecB_temp))  deallocate(all_RvecB_temp)
    IF(allocated(all_RvecB_temp2)) deallocate(all_RvecB_temp2)
    DO itab=1,size_psi
      CALL dealloc_TypeRVec(PsiR(itab))
    END DO

#endif
  END SUBROUTINE Action_MPI_S4
!=======================================================================================  


!=======================================================================================
!> sub_TabOpPsi_FOR_SGtype4 for working on full Smolyak rep.
!======================================================================================= 
  SUBROUTINE sub_TabOpPsi_FOR_SGtype4_SRB_MPI(Psi,OpPsi,para_Op)
    USE mod_system
    USE mod_nDindex
    USE mod_psi_set_alloc,  ONLY:param_psi
    USE mod_basis_set_alloc,ONLY:basis
    USE mod_SetOp,          ONLY:param_Op
    USE mod_MPI_aux
    IMPLICIT NONE

    TYPE(param_psi),                intent(in)    :: psi(:)
    TYPE(param_psi),                intent(inout) :: OpPsi(:)
    TYPE(param_Op),                 intent(inout) :: para_Op
    
    TYPE(basis),pointer                           :: BasisnD
    Integer                                       :: iG
    Integer                                       :: itab
    Integer                                       :: iterm00
    Integer                                       :: size_psi
    Integer                                       :: ii
    
#if(run_MPI)

    BasisnD => para_Op%BasisnD

    nb_mult_OpPsi=0
    size_psi=size(psi)

    ! initialize OpPsi
    DO ii=1,size_psi
      OpPsi(ii)=Psi(ii)
      OpPsi(ii)=ZERO
    ENDDO

    !DO iG=1,BasisnD%para_SGType2%nb_SG
    DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)
      CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4_SRB_MPI(iG,psi,OpPsi,                     &
                          BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval(:,iG),para_Op)
    ENDDO

    iterm00 = para_Op%derive_term_TO_iterm(0,0)
    IF (associated(para_Op%OpGrid)) THEN
      para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_done =                        &
                        para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid
      para_Op%OpGrid(iterm00)%para_FileGrid%Save_FileGrid_done =                       &
                        para_Op%OpGrid(iterm00)%para_FileGrid%Save_FileGrid
      para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_done=.True. 
    END IF

    ! set symab
  !  IF(MPI_id==0) THEN
  !    DO i=1,size(OpPsi)
  !      OpPsi_symab=Calc_symab1_EOR_symab2(para_Op%symab,Psi(i)%symab)
  !      CALL Set_symab_OF_psiBasisRep(OpPsi(i),OpPsi_symab)
  !    END DO
  !  ENDIF

#endif
  ENDSUBROUTINE sub_TabOpPsi_FOR_SGtype4_SRB_MPI
!=======================================================================================

!=======================================================================================
!> sub_TabOpPsi_OF_ONEDP_FOR_SGtype4 for working on fulll Smolyak rep. on basis
!=======================================================================================
  SUBROUTINE sub_TabOpPsi_OF_ONEDP_FOR_SGtype4_SRB_MPI(iG,psi,OpPsi,tab_l,para_Op)
  USE mod_system
  USE mod_nDindex
  USE mod_Coord_KEO,              ONLY:CoordType
  USE mod_basis_set_alloc,        ONLY:basis
  USE mod_OpPsi_SG4,              ONLY:get_OpGrid_type0_OF_ONEDP_FOR_SG4,              &
                                       get_OpGrid_type10_OF_ONEDP_FOR_SG4,             &
                                       get_OpGrid_type1_OF_ONEDP_FOR_SG4_new2
  USE mod_basis_BtoG_GtoB_SGType4,ONLY:TypeRVec,DerivOp_TO_RDP_OF_SmolaykRep,          &
                                       getbis_tab_nq,getbis_tab_nb,                    &
                                       BDP_TO_GDP_OF_SmolyakRep,                       &
                                       GDP_TO_BDP_OF_SmolyakRep
  USE mod_SetOp,                  ONLY:param_Op
  USE mod_psi_set_alloc,          ONLY:param_psi
  USE mod_MPI_aux
  IMPLICIT NONE


  TYPE(param_Op),                 intent(inout) :: para_Op
  TYPE(param_psi),                intent(in)    :: psi(:)
  TYPE(param_psi),                intent(inout) :: OpPsi(:)
  Integer,                        intent(in)    :: iG
  Integer,                        intent(in)    :: tab_l(:)
  
  TYPE(CoordType),pointer                       :: mole
  TYPE(basis),pointer                           :: BasisnD
  Real(kind=Rkind),allocatable                  :: Op_Psi(:,:) ! size(nq,nb0)
  Real(kind=Rkind),allocatable                  :: Psi_ch(:,:) ! size(nq,nb0)
  Real(kind=Rkind),allocatable                  :: PsiR(:)
  Real(kind=Rkind),allocatable                  :: V(:,:,:)
  Real(kind=Rkind),allocatable                  :: PsiRj(:,:)
  Real(kind=Rkind),allocatable                  :: PsiRi(:)
  Real(kind=Rkind),allocatable                  :: OpPsiR(:)
  Real(kind=Rkind),allocatable                  :: VPsi(:,:)
  Real(kind=Rkind),allocatable                  :: GGiq(:,:,:)
  Real(kind=Rkind),allocatable                  :: GridOp(:,:,:,:)
  Real(kind=Rkind),allocatable                  :: sqRhoOVERJac(:)
  Real(kind=Rkind),allocatable                  :: Jac(:)
  Integer,allocatable                           :: tab_nq(:)
  Integer,allocatable                           :: tab_nb(:)
  Integer                                       :: derive_termQdyn(2)
  Integer                                       :: Dim
  Integer                                       :: nb
  Integer                                       :: nq
  Integer                                       :: nb0
  Integer                                       :: ib0
  Integer                                       :: jb0
  
  Integer                                       :: iqi
  Integer                                       :: iqf
  Integer                                       :: jqi
  Integer                                       :: jqf
  
  Integer                                       :: d1
  Integer                                       :: d2
  Integer                                       :: num
  Integer                                       :: ii
  Integer                                       :: i
  Integer                                       :: j
  Integer                                       :: itab
  Integer                                       :: iterm
  
  Character(len=*),parameter             :: name_sub='sub_TabOpPsi_OF_ONEDP_FOR_SGtype4'

#if(run_MPI)

  num=1
  IF(psi(1)%cplx) num=2
  
  mole    => para_Op%mole
  BasisnD => para_Op%BasisnD

  d1=psi(1)%SR_B_index(iG)
  d2=psi(1)%SR_B_index(iG+1)-1
  CALL alloc_NParray(PsiR,(/d2-d1+1/),'PsiR',name_sub)

  Dim=size(tab_l)
  CALL alloc_NParray(tab_nb,(/Dim/),'tab_nb',name_sub)
  CALL alloc_NParray(tab_nq,(/Dim/),'tab_nq',name_sub)
  
  tab_nq(:)=getbis_tab_nq(tab_l,BasisnD%tab_basisPrimSG)
  tab_nb(:)=getbis_tab_nb(tab_l,BasisnD%tab_basisPrimSG)
  
  nb =BasisnD%para_SGType2%tab_nb_OF_SRep(iG)
  nq =BasisnD%para_SGType2%tab_nq_OF_SRep(iG)
  nb0=BasisnD%para_SGType2%nb0

  ! NOTE, need further optimize here for matrix manipulate
  SELECT CASE (para_Op%type_Op)
  CASE (0) !-0:Scalar-------------------------------------------------------------------

    CALL get_OpGrid_type0_OF_ONEDP_FOR_SG4(iG,tab_l,para_Op,V)
    CALL alloc_NParray(Op_Psi,(/nq,nb0/),'Op_Psi',name_sub)
    
    DO itab=1,size(psi)
      DO ii=1,num
        PsiR=psi(itab)%SR_B(d1:d2,ii)
        ! partial B to G
        CALL BDP_TO_GDP_OF_SmolyakRep(PsiR,BasisnD%tab_basisPrimSG,                    &
                                      tab_l,tab_nq,tab_nb,nb0)

        ! VPsi calculation, it has to be done before because V is not diagonal
        Op_Psi(:,:)=ZERO
        DO ib0=1,nb0
          jqi=1
          jqf=nq
          DO jb0=1,nb0
            ! was PsiR(itab)%V(jqi:jqf)
            Op_Psi(:,ib0)=Op_Psi(:,ib0)+V(:,ib0,jb0)*PsiR(jqi:jqf)
            jqi=jqf+1
            jqf=jqf+nq
          END DO
        END DO
        
        PsiR=reshape(Op_Psi,shape=(/nq*nb0/))
        CALL GDP_TO_BDP_OF_SmolyakRep(PsiR,BasisnD%tab_basisPrimSG,&
                                      tab_l,tab_nq,tab_nb,nb0)
        OpPsi(itab)%SR_B(d1:d2,ii)=PsiR
      ENDDO ! for ii
    ENDDO ! for itab

    IF(allocated(V)) CALL dealloc_NParray(V,'V',name_sub)
    IF(allocated(Op_Psi)) CALL dealloc_NParray(Op_Psi,'Op_Psi',name_sub)

  CASE (1) !-1: H: F2.d^2 + F1.d^1 + V--------------------------------------------------

    CALL get_OpGrid_type1_OF_ONEDP_FOR_SG4_new2(iG,tab_l,para_Op,GridOp)
    CALL alloc_NParray(Op_Psi,(/nq,nb0/),'Op_Psi',name_sub)
    CALL alloc_NParray(Psi_ch,(/nq,nb0/),'Psi_ch',name_sub)

    DO itab=1,size(psi)
      DO ii=1,num
        PsiR=psi(itab)%SR_B(d1:d2,ii)
        CALL BDP_TO_GDP_OF_SmolyakRep(PsiR,BasisnD%tab_basisPrimSG,                    &
                                      tab_l,tab_nq,tab_nb,nb0)

        Op_Psi(:,:)=ZERO
        DO iterm=1,para_Op%nb_Term
          IF(para_Op%OpGrid(iterm)%grid_zero) CYCLE
          Psi_ch(:,:)=reshape(PsiR,shape=(/nq,nb0/))
          DO ib0=1,nb0
            CALL DerivOp_TO_RDP_OF_SmolaykRep(Psi_ch(:,ib0),BasisnD%tab_basisPrimSG,   &
                                          tab_l,tab_nq,para_Op%derive_termQdyn(:,iterm))
          ENDDO
          ! GridOp(:,1,1,iterm)Psi_ch calculation
          DO ib0=1,nb0
            DO jb0=1,nb0
              Op_Psi(:,ib0)=Op_Psi(:,ib0)+GridOp(:,ib0,jb0,iterm)*Psi_ch(:,jb0)
            END DO
          END DO
        ENDDO ! for iterm=1,para_Op%nb_Term
        PsiR=reshape(Op_Psi,shape=(/nq*nb0/))
        CALL GDP_TO_BDP_OF_SmolyakRep(PsiR,BasisnD%tab_basisPrimSG,                    &
                                      tab_l,tab_nq,tab_nb,nb0)
        OpPsi(itab)%SR_B(d1:d2,ii)=PsiR
      ENDDO ! for ii
    ENDDO ! for itab

    IF(allocated(GridOp)) CALL dealloc_NParray(GridOp,'GridOp',name_sub)
    IF(allocated(Op_Psi)) CALL dealloc_NParray(Op_Psi,'Op_Psi',name_sub)
    IF(allocated(Psi_ch)) CALL dealloc_NParray(Psi_ch,'Psi_ch',name_sub)

  CASE (10) !-10: H: d^1 G d^1 +V-------------------------------------------------------

    CALL alloc_NParray(PsiRj,(/nq,mole%nb_act1 /),'PsiRj', name_sub)
    CALL alloc_NParray(PsiRi,(/nq /),             'PsiRi', name_sub)
    CALL alloc_NParray(VPsi, (/nq,nb0/),          'VPsi',  name_sub)
    CALL alloc_NParray(OpPsiR,(/nq/),             'OpPsiR',name_sub)

    CALL get_OpGrid_type10_OF_ONEDP_FOR_SG4(iG,tab_l,para_Op,V,GGiq,sqRhoOVERJac,Jac)

    DO itab=1,size(psi)
      DO ii=1,num
        PsiR=psi(itab)%SR_B(d1:d2,ii)
        CALL BDP_TO_GDP_OF_SmolyakRep(PsiR,BasisnD%tab_basisPrimSG,&
                                      tab_l,tab_nq,tab_nb,nb0)

        ! VPsi calculation, it has to be done before because V is not diagonal
        VPsi(:,:)=ZERO
        DO ib0=1,nb0
          jqi=1
          jqf=nq
          DO jb0=1,nb0
            VPsi(:,ib0)=VPsi(:,ib0)+V(:,ib0,jb0)*PsiR(jqi:jqf)
            jqi=jqf+1
            jqf=jqf+nq
          ENDDO
        ENDDO
        
        iqi=1
        iqf=nq
        DO ib0=1,nb0
          ! multiplication by sqRhoOVERJac
          PsiR(jqi:jqf)=PsiR(jqi:jqf)*sqRhoOVERJac(:)
          
          ! derivative with respect to Qj
          DO j=1,mole%nb_act1
            derive_termQdyn(:)=(/mole%liste_QactTOQdyn(j),0/)

            PsiRj(:,j)=PsiR(jqi:jqf)
            CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiRj(:,j),BasisnD%tab_basisPrimSG,      &
                                              tab_l,tab_nq,derive_termQdyn)
          ENDDO ! for j=1,mole%nb_act1
          
          OpPsiR(:)=ZERO
          DO i=1,mole%nb_act1
            PsiRi(:)=ZERO
            DO j=1,mole%nb_act1
              PsiRi(:)=PsiRi(:)+GGiq(:,j,i)*PsiRj(:,j)
            ENDDO
            PsiRi(:)=PsiRi(:)*Jac(:)

            derive_termQdyn(:)=(/mole%liste_QactTOQdyn(i),0/)

            CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiRi(:),BasisnD%tab_basisPrimSG,        &
                                              tab_l,tab_nq,derive_termQdyn)
            OpPsiR(:)=OpPsiR(:)+PsiRi(:)
          ENDDO ! for i=1,mole%nb_act1
          
          OpPsiR(:)=-HALF*OpPsiR(:)/(Jac(:)*sqRhoOVERJac(:))+VPsi(:,ib0)
          PsiR(jqi:jqf)=OpPsiR(:)
          
          iqi=iqf+1
          iqf=iqf+nq
        ENDDO ! for ib0=1,nb0
        
        CALL GDP_TO_BDP_OF_SmolyakRep(PsiR,BasisnD%tab_basisPrimSG,&
                                    tab_l,tab_nq,tab_nb,nb0)
        OpPsi(itab)%SR_B(d1:d2,ii)=PsiR
      ENDDO ! for ii
    ENDDO ! for itab

    IF(allocated(PsiRi))        CALL dealloc_NParray(PsiRi,'PsiRi',name_sub)
    IF(allocated(PsiRj))        CALL dealloc_NParray(PsiRj,'PsiRj',name_sub)
    IF(allocated(GGiq))         CALL dealloc_NParray(GGiq, 'GGiq', name_sub)
    IF(allocated(sqRhoOVERJac)) CALL dealloc_NParray(sqRhoOVERJac,'sqRhoOVERJac',name_sub)
    IF(allocated(Jac))          CALL dealloc_NParray(Jac,  'Jac',  name_sub)
    IF(allocated(V))            CALL dealloc_NParray(V,    'V',    name_sub)
    IF(allocated(VPsi))         CALL dealloc_NParray(VPsi, 'VPsi', name_sub)

  CASE Default !------------------------------------------------------------------------
    STOP 'error: case in sub_TabOpPsi_OF_ONEDP_FOR_SGtype4_SRB_MPI'
  END SELECT

  IF(allocated(tab_nb)) CALL dealloc_NParray(tab_nb,'tab_nb',name_sub)
  IF(allocated(tab_nq)) CALL dealloc_NParray(tab_nq,'tab_nq',name_sub)
  IF(allocated(PsiR))   CALL dealloc_NParray(PsiR,'tab_nq',name_sub)

#endif

  END SUBROUTINE sub_TabOpPsi_OF_ONEDP_FOR_SGtype4_SRB_MPI
!=======================================================================================

!=======================================================================================
!> sub_TabOpPsi_OF_ONEDP_FOR_SGtype4 for working on fulll Smolyak rep.
!=======================================================================================
  SUBROUTINE sub_TabOpPsi_OF_ONEDP_FOR_SGtype4_SRG_MPI(iG,psi,OpPsi,tab_l,para_Op)
  USE mod_system
  USE mod_nDindex
  USE mod_Coord_KEO,              ONLY:CoordType
  USE mod_basis_set_alloc,        ONLY:basis
  USE mod_OpPsi_SG4,              ONLY:get_OpGrid_type0_OF_ONEDP_FOR_SG4,              &
                                       get_OpGrid_type10_OF_ONEDP_FOR_SG4,             &
                                       get_OpGrid_type1_OF_ONEDP_FOR_SG4_new2
  USE mod_basis_BtoG_GtoB_SGType4,ONLY:TypeRVec,DerivOp_TO_RDP_OF_SmolaykRep,          &
                                       getbis_tab_nq,getbis_tab_nb
  USE mod_SetOp,                  ONLY:param_Op
  USE mod_psi_set_alloc,          ONLY:param_psi
  USE mod_MPI_aux
  IMPLICIT NONE


  TYPE(param_Op),                 intent(inout) :: para_Op
  TYPE(param_psi),                intent(in)    :: psi(:)
  TYPE(param_psi),                intent(inout) :: OpPsi(:)
  Integer,                        intent(in)    :: iG
  Integer,                        intent(in)    :: tab_l(:)
  
  TYPE(CoordType),pointer                       :: mole
  TYPE(basis),pointer                           :: BasisnD
  Real(kind=Rkind),allocatable                  :: Op_Psi(:,:) ! size(nq,nb0)
  Real(kind=Rkind),allocatable                  :: Psi_ch(:,:) ! size(nq,nb0)
  Real(kind=Rkind),allocatable                  :: PsiR(:)
  Real(kind=Rkind),allocatable                  :: V(:,:,:)
  Real(kind=Rkind),allocatable                  :: PsiRj(:,:)
  Real(kind=Rkind),allocatable                  :: PsiRi(:)
  Real(kind=Rkind),allocatable                  :: OpPsiR(:)
  Real(kind=Rkind),allocatable                  :: VPsi(:,:)
  Real(kind=Rkind),allocatable                  :: GGiq(:,:,:)
  Real(kind=Rkind),allocatable                  :: GridOp(:,:,:,:)
  Real(kind=Rkind),allocatable                  :: sqRhoOVERJac(:)
  Real(kind=Rkind),allocatable                  :: Jac(:)
  Integer,allocatable                           :: tab_nq(:)
  Integer,allocatable                           :: tab_nb(:)
  Integer                                       :: derive_termQdyn(2)
  Integer                                       :: Dim
  Integer                                       :: nb
  Integer                                       :: nq
  Integer                                       :: nb0
  Integer                                       :: ib0
  Integer                                       :: jb0
  
  Integer                                       :: iqi
  Integer                                       :: iqf
  Integer                                       :: jqi
  Integer                                       :: jqf
  
  Integer                                       :: d1
  Integer                                       :: d2
  Integer                                       :: num
  Integer                                       :: ii
  Integer                                       :: i
  Integer                                       :: j
  Integer                                       :: itab
  Integer                                       :: iterm
  
  Character(len=*),parameter             :: name_sub='sub_TabOpPsi_OF_ONEDP_FOR_SGtype4'

#if(run_MPI)

  num=1
  IF(psi(1)%cplx) num=2
  
  mole    => para_Op%mole
  BasisnD => para_Op%BasisnD

  d1=psi(1)%SR_G_index(iG)
  d2=psi(1)%SR_G_index(iG+1)-1
  CALL alloc_NParray(PsiR,(/d2-d1+1/),'PsiR',name_sub)

  Dim=size(tab_l)
  CALL alloc_NParray(tab_nb,(/Dim/),'tab_nb',name_sub)
  CALL alloc_NParray(tab_nq,(/Dim/),'tab_nq',name_sub)
  
  tab_nq(:)=getbis_tab_nq(tab_l,BasisnD%tab_basisPrimSG)
  tab_nb(:)=getbis_tab_nb(tab_l,BasisnD%tab_basisPrimSG)
  
  nb =BasisnD%para_SGType2%tab_nb_OF_SRep(iG)
  nq =BasisnD%para_SGType2%tab_nq_OF_SRep(iG)
  nb0=BasisnD%para_SGType2%nb0

  ! NOTE, need further optimize here for matrix manipulate
  SELECT CASE (para_Op%type_Op)
  CASE (0) !-0:Scalar-------------------------------------------------------------------
    CALL get_OpGrid_type0_OF_ONEDP_FOR_SG4(iG,tab_l,para_Op,V)

    !CALL allocate_array(Op_Psi,1,nq,1,nb0)
    CALL alloc_NParray(Op_Psi,(/nq,nb0/),'Op_Psi',name_sub)
    
    DO itab=1,size(psi)
      DO ii=1,num
        Op_Psi(:,:)=ZERO
      
        ! partial B to G
        !CALL BDP_TO_GDP_OF_SmolyakRep(PsiR(itab)%V,BasisnD%tab_basisPrimSG,&
        !                              tab_l,tab_nq,tab_nb,nb0)
        ! PsiR(itab)%V --> psi%SR_G(psi%SR_G_index(iG):psi%SR_G_index(iG+1)-1)

        ! VPsi calculation, it has to be done before because V is not diagonal
        DO ib0=1,nb0
          jqi=1
          jqf=nq
          DO jb0=1,nb0
            ! was PsiR(itab)%V(jqi:jqf)
            Op_Psi(:,ib0)=Op_Psi(:,ib0)+V(:,ib0,jb0)*psi(itab)%SR_G(d1-1+jqi:d1-1+jqf,ii)
            jqi=jqf+1
            jqf=jqf+nq
          END DO
        END DO
        OpPsi(itab)%SR_G(d1:d2,ii)=reshape(Op_Psi,shape=(/nq*nb0/))
      ENDDO ! for ii
    ENDDO ! for itab

    IF(allocated(V)) CALL dealloc_NParray(V,'V',name_sub)
    IF(allocated(Op_Psi)) CALL dealloc_NParray(Op_Psi,'Op_Psi',name_sub)

  CASE (1) !-1: H: F2.d^2 + F1.d^1 + V--------------------------------------------------
    CALL get_OpGrid_type1_OF_ONEDP_FOR_SG4_new2(iG,tab_l,para_Op,GridOp)

    CALL alloc_NParray(Op_Psi,(/nq,nb0/),'Op_Psi',name_sub)
    CALL alloc_NParray(Psi_ch,(/nq,nb0/),'Psi_ch',name_sub)

    DO itab=1,size(psi)
      DO ii=1,num
        Op_Psi(:,:)=ZERO
        DO iterm=1,para_Op%nb_Term
          IF(para_Op%OpGrid(iterm)%grid_zero) CYCLE
          Psi_ch(:,:)=reshape(psi(itab)%SR_G(d1:d2,ii),shape=(/nq,nb0/))

          DO ib0=1,nb0
            CALL DerivOp_TO_RDP_OF_SmolaykRep(Psi_ch(:,ib0),BasisnD%tab_basisPrimSG,   &
                                          tab_l,tab_nq,para_Op%derive_termQdyn(:,iterm))
          ENDDO

          ! GridOp(:,1,1,iterm)Psi_ch calculation
          DO ib0=1,nb0
            DO jb0=1,nb0
              Op_Psi(:,ib0)=Op_Psi(:,ib0)+GridOp(:,ib0,jb0,iterm)*Psi_ch(:,jb0)
            END DO
          END DO
        ENDDO ! for iterm=1,para_Op%nb_Term

        OpPsi(itab)%SR_G(d1:d2,ii)= reshape(Op_Psi,shape=(/nq*nb0/))
      ENDDO ! for ii
    ENDDO ! for itab

    IF(allocated(GridOp)) CALL dealloc_NParray(GridOp,'GridOp',name_sub)
    IF(allocated(Op_Psi)) CALL dealloc_NParray(Op_Psi,'Op_Psi',name_sub)
    IF(allocated(Psi_ch)) CALL dealloc_NParray(Psi_ch,'Psi_ch',name_sub)

  CASE (10) !-10: H: d^1 G d^1 +V-------------------------------------------------------

    CALL alloc_NParray(PsiRj,(/nq,mole%nb_act1 /),'PsiRj', name_sub)
    CALL alloc_NParray(PsiRi,(/nq /),             'PsiRi', name_sub)
    CALL alloc_NParray(VPsi, (/nq,nb0/),          'VPsi',  name_sub)
    CALL alloc_NParray(OpPsiR,(/nq/),             'OpPsiR',name_sub)

    CALL get_OpGrid_type10_OF_ONEDP_FOR_SG4(iG,tab_l,para_Op,V,GGiq,sqRhoOVERJac,Jac)

    DO itab=1,size(psi)
      DO ii=1,num
        PsiR=psi(itab)%SR_G(d1:d2,ii)
        ! VPsi calculation, it has to be done before because V is not diagonal
        VPsi(:,:)=ZERO
        DO ib0=1,nb0
          jqi=1
          jqf=nq
          DO jb0=1,nb0
            VPsi(:,ib0)=VPsi(:,ib0)+V(:,ib0,jb0)*PsiR(jqi:jqf)
            jqi=jqf+1
            jqf=jqf+nq
          ENDDO
        ENDDO
        
        iqi=1
        iqf=nq
        DO ib0=1,nb0
          ! multiplication by sqRhoOVERJac
          PsiR(jqi:jqf)=PsiR(jqi:jqf)*sqRhoOVERJac(:)
          
          ! derivative with respect to Qj
          DO j=1,mole%nb_act1
            derive_termQdyn(:)=(/mole%liste_QactTOQdyn(j),0/)

            PsiRj(:,j)=PsiR(jqi:jqf)
            CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiRj(:,j),BasisnD%tab_basisPrimSG,      &
                                              tab_l,tab_nq,derive_termQdyn)
          ENDDO ! for j=1,mole%nb_act1
          
          OpPsiR(:)=ZERO
          DO i=1,mole%nb_act1
            PsiRi(:)=ZERO
            DO j=1,mole%nb_act1
              PsiRi(:)=PsiRi(:)+GGiq(:,j,i)*PsiRj(:,j)
            ENDDO
            PsiRi(:)=PsiRi(:)*Jac(:)

            derive_termQdyn(:)=(/mole%liste_QactTOQdyn(i),0/)

            CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiRi(:),BasisnD%tab_basisPrimSG,        &
                                              tab_l,tab_nq,derive_termQdyn)
            OpPsiR(:)=OpPsiR(:)+PsiRi(:)
          ENDDO ! for i=1,mole%nb_act1
          
          OpPsiR(:)=-HALF*OpPsiR(:)/(Jac(:)*sqRhoOVERJac(:))+VPsi(:,ib0)
          PsiR(jqi:jqf)=OpPsiR(:)
          
          iqi=iqf+1
          iqf=iqf+nq
        ENDDO ! for ib0=1,nb0
        OpPsi(itab)%SR_G(d1:d2,ii)=PsiR
      ENDDO ! for ii
    ENDDO ! for itab

    IF(allocated(PsiRi))        CALL dealloc_NParray(PsiRi,'PsiRi',name_sub)
    IF(allocated(PsiRj))        CALL dealloc_NParray(PsiRj,'PsiRj',name_sub)
    IF(allocated(GGiq))         CALL dealloc_NParray(GGiq, 'GGiq', name_sub)
    IF(allocated(sqRhoOVERJac)) CALL dealloc_NParray(sqRhoOVERJac,'sqRhoOVERJac',name_sub)
    IF(allocated(Jac))          CALL dealloc_NParray(Jac,  'Jac',  name_sub)
    IF(allocated(V))            CALL dealloc_NParray(V,    'V',    name_sub)
    IF(allocated(VPsi))         CALL dealloc_NParray(VPsi, 'VPsi', name_sub)

  CASE Default !------------------------------------------------------------------------
    STOP 'error: case in sub_TabOpPsi_OF_ONEDP_FOR_SGtype4_SRG_MPI'
  END SELECT

  IF(allocated(tab_nb)) CALL dealloc_NParray(tab_nb,'tab_nb',name_sub)
  IF(allocated(tab_nq)) CALL dealloc_NParray(tab_nq,'tab_nq',name_sub)
  IF(allocated(PsiR))   CALL dealloc_NParray(PsiR,'tab_nq',name_sub)

#endif
  END SUBROUTINE sub_TabOpPsi_OF_ONEDP_FOR_SGtype4_SRG_MPI
!=======================================================================================


!=======================================================================================
!> sub_TabOpPsi_FOR_SGtype4 for working on full Smolyak rep.
!======================================================================================= 
  SUBROUTINE sub_TabOpPsi_FOR_SGtype4_SRG_MPI(Psi,OpPsi,para_Op)
    USE mod_system
    USE mod_nDindex
    USE mod_psi_set_alloc,  ONLY:param_psi
    USE mod_basis_set_alloc,ONLY:basis
    USE mod_SetOp,          ONLY:param_Op
    USE mod_MPI_aux
    IMPLICIT NONE


    TYPE(param_psi),                intent(in)    :: psi(:)
    TYPE(param_psi),                intent(inout) :: OpPsi(:)
    TYPE(param_Op),                 intent(inout) :: para_Op
    
    TYPE(basis),pointer                           :: BasisnD
    Integer                                       :: iG
    Integer                                       :: itab
    Integer                                       :: iterm00
    Integer                                       :: size_psi
    Integer                                       :: ii

#if(run_MPI)

    BasisnD => para_Op%BasisnD

    nb_mult_OpPsi=0
    size_psi=size(psi)

    ! initialize OpPsi
    DO ii=1,size_psi
      OpPsi(ii)=Psi(ii)
      OpPsi(ii)=ZERO
    ENDDO

    !DO iG=1,BasisnD%para_SGType2%nb_SG
    DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)
      CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4_SRG_MPI(iG,psi,OpPsi,                     &
                          BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval(:,iG),para_Op)
    ENDDO

    iterm00 = para_Op%derive_term_TO_iterm(0,0)
    IF (associated(para_Op%OpGrid)) THEN
      para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_done =                        &
                        para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid
      para_Op%OpGrid(iterm00)%para_FileGrid%Save_FileGrid_done =                       &
                        para_Op%OpGrid(iterm00)%para_FileGrid%Save_FileGrid
      para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_done=.True. 
    END IF

  ! Warning, set symab when converge back to packed basis
!  IF(MPI_id==0) THEN
!    DO i=1,size(OpPsi)
!      OpPsi_symab=Calc_symab1_EOR_symab2(para_Op%symab,Psi(i)%symab)
!      CALL Set_symab_OF_psiBasisRep(OpPsi(i),OpPsi_symab)
!    END DO
!  ENDIF

#endif
  ENDSUBROUTINE sub_TabOpPsi_FOR_SGtype4_SRG_MPI
!=======================================================================================


!=======================================================================================
!> auto adjust the distribution of iG on each threads according to 
!> time used in the previous action
!=======================================================================================
  SUBROUTINE auto_iGs_MPI_old(para_Op)
    USE mod_system
    USE mod_SetOp,ONLY:param_Op
    USE mod_MPI_aux
    IMPLICIT NONE
    
    TYPE(param_Op),                 intent(in)    :: para_Op
    
    Real                                          :: ave_time
    Real                                          :: rest_time
    Real                                          :: temp_time
    Real                                          :: ave_iGtime_threads(0:MPI_np-1)
    Real                                          :: master_commu_time
    Integer                                       :: iG_threads(0:MPI_np-1)
    Integer                                       :: iG_change
    Integer                                       :: iG_total
    Integer                                       :: ii
    Integer                                       :: i_MPI2
    Integer                                       :: current_i_MPI
    Integer                                       :: num_iG

#if(run_MPI)

    CALL MPI_collect_info(time_MPI_act_all) ! time_MPI_act_all infor on master

    IF(MPI_id==0) THEN 
      write(out_unitp,*) 'time used in action for each threads:',time_MPI_act_all
    ENDIF
    
    IF(MPI_id==0) THEN
      iG_total=para_Op%BasisnD%para_SGType2%nb_SG

      ! get the average time spend on each threads, 
      ! which is the ideal time for each threads
      ave_time=Real(SUM(time_MPI_act_all),kind=Rkind)/Real(MPI_np,kind=Rkind)
      
      ! get the average time for each iG on different threads
      DO i_MPI=1,MPI_np-1
        ave_iGtime_threads(i_MPI)=Real(time_MPI_act_all(i_MPI),kind=Rkind)             &
                                   /Real(iGs_MPI(2,i_MPI)-iGs_MPI(1,i_MPI)+1,kind=Rkind)
      ENDDO
      ! master takes more time on communcation. 
      ! assume the ave time for iGs is same as threads=1
      ave_iGtime_threads(0)=ave_iGtime_threads(1)
      master_commu_time=MAX(0.,Real(time_MPI_act_all(0)-time_MPI_act_all(1),kind=Rkind))

      ! may consider to use inter Trapezoid
      current_i_MPI=0
      DO i_MPI=0,MPI_np-2
        rest_time=0
        temp_time=0
        IF(i_MPI==0) THEN
          rest_time=Real(time_MPI_act_all(i_MPI),kind=Rkind)-master_commu_time
          num_iG=iGs_MPI(2,i_MPI)-iGs_MPI(1,i_MPI)+1
          temp_time=master_commu_time
        ENDIF
        
        iG_threads(i_MPI)=0
        DO i_MPI2=current_i_MPI,MPI_np-1
          IF(temp_time+rest_time<=ave_time) THEN
            temp_time=temp_time+rest_time
            iG_threads(i_MPI)=iG_threads(i_MPI)+num_iG
            
            ! go to next block
            IF(i_MPI2<=MPI_np-2) THEN
              rest_time=Real(time_MPI_act_all(i_MPI2+1),kind=Rkind)
              num_iG=iGs_MPI(2,i_MPI2+1)-iGs_MPI(1,i_MPI2+1)+1
            ENDIF
          ELSE
            ! make sure there is at least 1 Smolyak term for each threads
            temp_int=MAX(1,FLOOR((ave_time-temp_time)/ave_iGtime_threads(i_MPI2)))
            iG_threads(i_MPI)=iG_threads(i_MPI)+temp_int

            ! record the current block
            current_i_MPI=i_MPI2
            ! record the rest time and iGs in current block
            rest_time=Real(time_MPI_act_all(i_MPI2),kind=Rkind)-(ave_time-temp_time)
            num_iG=iGs_MPI(2,i_MPI2)-iGs_MPI(1,i_MPI2)+1-temp_int
            EXIT
          ENDIF
        ENDDO ! for i_MPI2=current_i_MPI,MPI_np-1
      ENDDO ! for i_MPI=0,MPI_np-2
      ! this make sure all iGs are accounted for 
      iG_threads(MPI_np-1)=iG_total-SUM(iG_threads(0:MPI_np-2))

      !> set up new iGs arrangement
      iGs_MPI(1,0)=1
      iGs_MPI(2,0)=iGs_MPI(1,0)+iG_threads(0)-1
      DO i_MPI=1,MPI_np-1
        iGs_MPI(1,i_MPI)=iGs_MPI(2,i_MPI-1)+1
        iGs_MPI(2,i_MPI)=iGs_MPI(1,i_MPI)+iG_threads(i_MPI)-1
      ENDDO
    ENDIF ! for MPI_id==0

    ! share to all threads
    CALL MPI_Bcast_matrix(iGs_MPI(1:2,0:MPI_np-1),1,2,0,MPI_np-1,root_MPI,shift2=1)

#endif
  ENDSUBROUTINE auto_iGs_MPI_old
!=======================================================================================


!=======================================================================================
!> auto adjust the distribution of iG on each threads according to 
!> time used in the previous action
!=======================================================================================
  SUBROUTINE auto_iGs_MPI_old2(para_Op)
    USE mod_system
    USE mod_SetOp,ONLY:param_Op
    USE mod_MPI_aux
    IMPLICIT NONE
    
    TYPE(param_Op),                 intent(in)    :: para_Op
    
    Real                                          :: ave_time
    Real                                          :: rest_time
    Real                                          :: temp_time
    Real                                          :: ave_iGtime_threads(0:MPI_np-1)
    Real                                          :: master_commu_time
    Integer                                       :: iG_threads(0:MPI_np-1)
    Integer                                       :: iG_change
    Integer                                       :: iG_total
    Integer                                       :: ii
    Integer                                       :: i_MPI2
    Integer                                       :: current_i_MPI
    Integer                                       :: num_iG
    Logical                                       :: returnall

#if(run_MPI)

    CALL MPI_collect_info(time_MPI_act_all) ! time_MPI_act_all infor on master

    returnall=.FALSE.
    IF(MPI_id==0) THEN 
      write(out_unitp,*) 'time used in action for each threads:',time_MPI_act_all
      ave_time=Real(SUM(time_MPI_act_all),kind=Rkind)/Real(MPI_np,kind=Rkind)
      IF(ALL(ABS(ave_time-time_MPI_act_all)/ave_time<0.12)) returnall=.TRUE.
    ENDIF

    CALL MPI_BCAST(returnall,size1_MPI,MPI_Logical,root_MPI,MPI_COMM_WORLD,MPI_err)
    IF(returnall) RETURN
    
    IF(MPI_id==0) THEN
      iG_total=para_Op%BasisnD%para_SGType2%nb_SG
      
      ! get the average time for each iG on different threads
      DO i_MPI=0,MPI_np-1
        ave_iGtime_threads(i_MPI)=Real(time_MPI_act_all(i_MPI),kind=Rkind)             &
                                   /Real(iGs_MPI(2,i_MPI)-iGs_MPI(1,i_MPI)+1,kind=Rkind)
      ENDDO
      master_commu_time=Real(time_MPI_act_all(0),kind=Rkind)                           &
                        -MAX(0.,Real(iGs_MPI(2,0)-iGs_MPI(1,0)+1,kind=Rkind)           &
                                *ave_iGtime_threads(MPI_np-1))

      ! get the average time spend on each threads, 
      ! which is the ideal time for each threads
      ave_time=(Real(SUM(time_MPI_act_all),kind=Rkind)-master_commu_time)              &
               /Real(MPI_np,kind=Rkind)

      DO ii=1,MPI_np
        ! may consider to use inter Trapezoid
        current_i_MPI=0
        rest_time=0
        DO i_MPI=0,MPI_np-2
          temp_time=0
          IF(i_MPI==0) THEN
            rest_time=Real(time_MPI_act_all(i_MPI),kind=Rkind)
            num_iG=MAX(1,iGs_MPI(2,i_MPI)-iGs_MPI(1,i_MPI)+1)
          ENDIF
          
          iG_threads(i_MPI)=0
          DO i_MPI2=current_i_MPI,MPI_np-1
            IF(temp_time+rest_time<=ave_time) THEN
              temp_time=temp_time+rest_time
              iG_threads(i_MPI)=iG_threads(i_MPI)+num_iG
              
              ! go to next block
              IF(i_MPI2<=MPI_np-2) THEN
                rest_time=Real(time_MPI_act_all(i_MPI2+1),kind=Rkind)
                num_iG=MAX(1,iGs_MPI(2,i_MPI2+1)-iGs_MPI(1,i_MPI2+1)+1)
              ENDIF
            ELSE
              ! make sure there is at least 1 Smolyak term for each threads
              temp_int=MAX(1,FLOOR((ave_time-temp_time)/ave_iGtime_threads(i_MPI2)))
              iG_threads(i_MPI)=iG_threads(i_MPI)+temp_int
              
              ! record the current block
              current_i_MPI=i_MPI2
              ! record the rest time and iGs in current block
              rest_time=Real(time_MPI_act_all(i_MPI2),kind=Rkind)-(ave_time-temp_time)
              num_iG=MAX(1,iGs_MPI(2,i_MPI2)-iGs_MPI(1,i_MPI2)+1-temp_int)
              EXIT
            ENDIF
          ENDDO ! for i_MPI2=current_i_MPI,MPI_np-1
        ENDDO ! for i_MPI=0,MPI_np-2
        
        IF(SUM(iG_threads(0:MPI_np-2))<iG_total) THEN
          EXIT
        ELSE
          ave_time=ave_time*(1.0-1.0/MPI_np)
        ENDIF
        
        IF(ii==MPI_np) STOP 'error in auto_iGs_MPI'
      ENDDO

      ! this make sure all iGs are accounted for
      iG_threads(MPI_np-1)=iG_total-SUM(iG_threads(0:MPI_np-2))

      !> set up new iGs arrangement
      iGs_MPI(1,0)=1
      iGs_MPI(2,0)=iGs_MPI(1,0)+iG_threads(0)-1
      DO i_MPI=1,MPI_np-1
        iGs_MPI(1,i_MPI)=iGs_MPI(2,i_MPI-1)+1
        iGs_MPI(2,i_MPI)=iGs_MPI(1,i_MPI)+iG_threads(i_MPI)-1
      ENDDO
    ENDIF ! for MPI_id==0

    ! share to all threads
    CALL MPI_Bcast_matrix(iGs_MPI(1:2,0:MPI_np-1),1,2,0,MPI_np-1,root_MPI,shift2=1)

#endif
  ENDSUBROUTINE auto_iGs_MPI_old2
!=======================================================================================


!=======================================================================================
!> auto adjust the distribution of iG on each threads according to 
!> time used in the previous action
!=======================================================================================
  SUBROUTINE auto_iGs_MPI_old3(para_Op)
    USE mod_system
    USE mod_SetOp,ONLY:param_Op
    USE mod_MPI_aux
    IMPLICIT NONE
    
    TYPE(param_Op),                 intent(in)    :: para_Op
    
    Real                                          :: ave_iGtime_threads(0:MPI_np-1)
    Real                                          :: ave_time
    Real                                          :: rest_time
    Real                                          :: count_time
    Integer                                       :: iG_threads(0:MPI_np-1)
    Integer                                       :: iG_total
    Integer                                       :: i_MPI2
    Integer                                       :: current_i_MPI
    Integer                                       :: num_iG
    Integer                                       :: ii
    Logical                                       :: returnall

#if(run_MPI)

    CALL MPI_collect_info(time_MPI_act_all) ! time_MPI_act_all infor on master
    CALL MPI_collect_info(time_MPI_local)

    returnall=.FALSE.
    IF(MPI_id==0) THEN 
      write(out_unitp,*) 'time used in action for each threads:',time_MPI_act_all
      !write(out_unitp,*) 'commu time used in action for each threads:',time_MPI_local
      ave_time=Real(SUM(time_MPI_act_all),kind=Rkind)/Real(MPI_np,kind=Rkind)
      IF(ALL(ABS(ave_time-time_MPI_act_all)/ave_time<0.1)) returnall=.TRUE.
      write(out_unitp,*) 'time balance:',ABS(ave_time-time_MPI_act_all)/ave_time<0.1
    ENDIF

    CALL MPI_BCAST(returnall,size1_MPI,MPI_Logical,root_MPI,MPI_COMM_WORLD,MPI_err)
    IF(returnall) RETURN
    
    IF(MPI_id==0) THEN
      time_MPI_calcu=time_MPI_act_all-time_MPI_local
      iG_total=para_Op%BasisnD%para_SGType2%nb_SG
      
      ! get the average time for each iG on different threads
      DO i_MPI=0,MPI_np-1
        ave_iGtime_threads(i_MPI)=Real(time_MPI_calcu(i_MPI),kind=Rkind)               &
                                   /Real(iGs_MPI(2,i_MPI)-iGs_MPI(1,i_MPI)+1,kind=Rkind)
      ENDDO

      ! get the average time spend on each threads, 
      ! which is the ideal time for each threads
      ave_time=Real(SUM(time_MPI_act_all),kind=Rkind)/Real(MPI_np,kind=Rkind)

      DO ii=1,MPI_np
        ! may consider to use inter Trapezoid
        current_i_MPI=0
        rest_time=0.
        DO i_MPI=0,MPI_np-2
          count_time=time_MPI_local(i_MPI)
          IF(i_MPI==0) THEN
            rest_time=Real(time_MPI_calcu(i_MPI),kind=Rkind)
            num_iG=MAX(1,iGs_MPI(2,i_MPI)-iGs_MPI(1,i_MPI)+1)
          ENDIF
          
          iG_threads(i_MPI)=0
          DO i_MPI2=current_i_MPI,MPI_np-1
            IF(rest_time+count_time<=ave_time) THEN
              count_time=count_time+rest_time
              iG_threads(i_MPI)=iG_threads(i_MPI)+num_iG

              ! go to next block
              IF(i_MPI2<=MPI_np-2) THEN
                rest_time=Real(time_MPI_calcu(i_MPI2+1),kind=Rkind)
                num_iG=MAX(1,iGs_MPI(2,i_MPI2+1)-iGs_MPI(1,i_MPI2+1)+1)
              ENDIF
            ELSE
              ! make sure there is at least 1 Smolyak term for each threads
              temp_int=MAX(1,FLOOR((ave_time-count_time)/ave_iGtime_threads(i_MPI2)))
              iG_threads(i_MPI)=iG_threads(i_MPI)+temp_int

              ! record the current block
              current_i_MPI=i_MPI2
              ! record the rest time and iGs in current block
              rest_time=Real(time_MPI_calcu(i_MPI2),kind=Rkind)-(ave_time-count_time)
              num_iG=MAX(1,iGs_MPI(2,i_MPI2)-iGs_MPI(1,i_MPI2)+1-temp_int)
              EXIT
            ENDIF
          ENDDO ! for i_MPI2=current_i_MPI,MPI_np-1
        ENDDO ! for i_MPI=0,MPI_np-2

        IF(SUM(iG_threads(0:MPI_np-2))<iG_total) THEN
          EXIT
        ELSE
          ave_time=ave_time*(1.0-1.0/MPI_np/10)
        ENDIF
        IF(ii==MPI_np) STOP 'error in auto_iGs_MPI'
      ENDDO

      ! this make sure all iGs are accounted for
      iG_threads(MPI_np-1)=iG_total-SUM(iG_threads(0:MPI_np-2))

      !> set up new iGs arrangement
      iGs_MPI(1,0)=1
      iGs_MPI(2,0)=iGs_MPI(1,0)+iG_threads(0)-1
      DO i_MPI=1,MPI_np-1
        iGs_MPI(1,i_MPI)=iGs_MPI(2,i_MPI-1)+1
        iGs_MPI(2,i_MPI)=iGs_MPI(1,i_MPI)+iG_threads(i_MPI)-1
      ENDDO
    ENDIF ! for MPI_id==0

    ! share to all threads
    CALL MPI_Bcast_matrix(iGs_MPI(1:2,0:MPI_np-1),1,2,0,MPI_np-1,root_MPI,shift2=1)

#endif
  ENDSUBROUTINE auto_iGs_MPI_old3
!=======================================================================================


!=======================================================================================
!> auto adjust the distribution of iG on each threads according to 
!> time used in the previous action
!=======================================================================================
  SUBROUTINE auto_iGs_MPI_old4(para_Op)
    USE mod_system
    USE mod_SetOp,ONLY:param_Op
    USE mod_MPI_aux
    IMPLICIT NONE
    
    TYPE(param_Op),                 intent(in)    :: para_Op
    
    Real                                          :: ave_iGtime_threads(0:MPI_np-1)
    Real                                          :: ave_time
    Real                                          :: rest_time
    Real                                          :: count_time
    Integer                                       :: iG_threads(0:MPI_np-1)
    Integer                                       :: iG_total
    Integer                                       :: i_MPI2
    Integer                                       :: current_i_MPI
    Integer                                       :: num_iG
    Integer                                       :: ii
    Logical                                       :: returnall

#if(run_MPI)

    CALL MPI_collect_info(time_MPI_act_all) ! time_MPI_act_all infor on master
    CALL MPI_collect_info(time_MPI_local)

    returnall=.FALSE.
    IF(MPI_id==0) THEN 
      write(out_unitp,*) 'time used in action for each threads:',time_MPI_act_all
      !write(out_unitp,*) 'commu time used in action for each threads:',time_MPI_local
      ave_time=Real(SUM(time_MPI_act_all),kind=Rkind)/Real(MPI_np,kind=Rkind)
      IF(ALL(ABS(ave_time-time_MPI_act_all)/ave_time<0.1)) returnall=.TRUE.
      write(out_unitp,*) 'time balance:',ABS(ave_time-time_MPI_act_all)/ave_time<0.1
    ENDIF

    CALL MPI_BCAST(returnall,size1_MPI,MPI_Logical,root_MPI,MPI_COMM_WORLD,MPI_err)
    IF(returnall) RETURN
    
    IF(MPI_id==0) THEN
      time_MPI_calcu=time_MPI_act_all-time_MPI_local
      iG_total=para_Op%BasisnD%para_SGType2%nb_SG
      
      ! get the average time for each iG on different threads
      DO i_MPI=0,MPI_np-1
        ave_iGtime_threads(i_MPI)=Real(time_MPI_calcu(i_MPI),kind=Rkind)               &
                                 /Real(iGs_MPI(2,i_MPI)-iGs_MPI(1,i_MPI)+1,kind=Rkind)
      ENDDO

      ! get the average time spend on each threads, 
      ! which is the ideal time for each threads
      ave_time=Real(SUM(time_MPI_act_all),kind=Rkind)/Real(MPI_np,kind=Rkind)

      DO ii=1,MPI_np
        ! may consider to use inter Trapezoid
        current_i_MPI=0
        rest_time=0.
        DO i_MPI=0,MPI_np-2
          count_time=time_MPI_local(i_MPI)
          IF(i_MPI==0) THEN
            rest_time=Real(time_MPI_calcu(i_MPI),kind=Rkind)
            num_iG=MAX(1,iGs_MPI(2,i_MPI)-iGs_MPI(1,i_MPI)+1)
          ENDIF
          
          iG_threads(i_MPI)=0
          DO i_MPI2=current_i_MPI,MPI_np-1
            IF(rest_time+count_time<=ave_time) THEN
              count_time=count_time+rest_time
              iG_threads(i_MPI)=iG_threads(i_MPI)+num_iG

              ! go to next block
              IF(i_MPI2<=MPI_np-2) THEN
                rest_time=Real(time_MPI_calcu(i_MPI2+1),kind=Rkind)
                num_iG=MAX(1,iGs_MPI(2,i_MPI2+1)-iGs_MPI(1,i_MPI2+1)+1)
              ENDIF
            ELSE
              ! make sure there is at least 1 Smolyak term for each threads
              temp_int=MAX(1,FLOOR((ave_time-count_time)/ave_iGtime_threads(i_MPI2)))
              iG_threads(i_MPI)=iG_threads(i_MPI)+temp_int

              ! record the current block
              current_i_MPI=i_MPI2
              ! record the rest time and iGs in current block
              rest_time=Real(time_MPI_calcu(i_MPI2),kind=Rkind)-(ave_time-count_time)
              num_iG=MAX(1,iGs_MPI(2,i_MPI2)-iGs_MPI(1,i_MPI2)+1-temp_int)
              EXIT
            ENDIF
          ENDDO ! for i_MPI2=current_i_MPI,MPI_np-1
        ENDDO ! for i_MPI=0,MPI_np-2

        IF(SUM(iG_threads(0:MPI_np-2))<iG_total) THEN
          EXIT
        ELSE
          ave_time=ave_time*(1.0-1.0/MPI_np/10)
        ENDIF
        IF(ii==MPI_np) STOP 'error in auto_iGs_MPI'
      ENDDO

      ! this make sure all iGs are accounted for
      iG_threads(MPI_np-1)=iG_total-SUM(iG_threads(0:MPI_np-2))

      !> set up new iGs arrangement
      iGs_MPI(1,0)=1
      iGs_MPI(2,0)=iGs_MPI(1,0)+iG_threads(0)-1
      DO i_MPI=1,MPI_np-1
        iGs_MPI(1,i_MPI)=iGs_MPI(2,i_MPI-1)+1
        iGs_MPI(2,i_MPI)=iGs_MPI(1,i_MPI)+iG_threads(i_MPI)-1
      ENDDO
    ENDIF ! for MPI_id==0

    ! share to all threads
    CALL MPI_Bcast_matrix(iGs_MPI(1:2,0:MPI_np-1),1,2,0,MPI_np-1,root_MPI,shift2=1)

#endif
  ENDSUBROUTINE auto_iGs_MPI_old4
!=======================================================================================

!=======================================================================================
!> auto adjust the distribution of iG on each threads according to 
!> time used in the previous action
!=======================================================================================
  SUBROUTINE auto_iGs_MPI_old5(para_Op)
    USE mod_system
    USE mod_SetOp,ONLY:param_Op
    USE mod_MPI_aux
    IMPLICIT NONE
    
    TYPE(param_Op),                 intent(in)    :: para_Op
    
    Real                                          :: ave_iGtime_threads(0:MPI_np-1)
    Real                                          :: ave_time
    Real                                          :: rest_time
    Real                                          :: count_time
    Real                                          :: master_local_time
    Integer                                       :: iG_threads(0:MPI_np-1)
    Integer                                       :: iG_total
    Integer                                       :: i_MPI2
    Integer                                       :: current_i_MPI
    Integer                                       :: num_iG
    Integer                                       :: ii
    Logical                                       :: returnall

#if(run_MPI)

    CALL MPI_collect_info(time_MPI_act_all) ! time_MPI_act_all infor on master
    CALL MPI_collect_info(time_MPI_calcu)
    !CALL MPI_collect_info(time_MPI_local)
    CALL MPI_Reduce_sum_Bcast(time_MPI_local(0:MPI_np-1),MPI_np)
    master_local_time=sum(time_MPI_local(1:MPI_np-1))
    !IF(MPI_id==0) time_MPI_local(0)=sum(time_MPI_local(1:MPI_np-1))
    !CALL MPI_Bcast_(time_MPI_local(0),size1_MPI,root_MPI)

    returnall=.FALSE.
    IF(MPI_id==0) THEN 
      write(out_unitp,*) 'time used in action for each processer:',time_MPI_act_all
      write(out_unitp,*) 'local comm time used for each processer:',time_MPI_local
      ave_time=Real(SUM(time_MPI_act_all),kind=Rkind)/Real(MPI_np,kind=Rkind)
      IF(ALL(ABS(ave_time-time_MPI_act_all)/ave_time<0.1)) returnall=.TRUE.
      write(out_unitp,*) 'time balance:',ABS(ave_time-time_MPI_act_all)/ave_time<0.1
    ENDIF

    CALL MPI_BCAST(returnall,size1_MPI,MPI_Logical,root_MPI,MPI_COMM_WORLD,MPI_err)
    IF(returnall) RETURN
    
    IF(MPI_id==0) THEN
      !time_MPI_calcu=time_MPI_act_all-time_MPI_local
      time_MPI_calcu=time_MPI_calcu+time_MPI_local
      time_MPI_calcu(0)=time_MPI_act_all(0)-time_MPI_local(0)
      iG_total=para_Op%BasisnD%para_SGType2%nb_SG
      
      ! get the average time for each iG on different threads
      DO i_MPI=0,MPI_np-1
        ave_iGtime_threads(i_MPI)=Real(time_MPI_calcu(i_MPI),kind=Rkind)               &
                                   /Real(iGs_MPI(2,i_MPI)-iGs_MPI(1,i_MPI)+1,kind=Rkind)
      ENDDO

      ! get the average time spend on each threads, 
      ! which is the ideal time for each threads
      ! ave_time=Real(SUM(time_MPI_act_all),kind=Rkind)/Real(MPI_np,kind=Rkind)
      ave_time=Real(SUM(time_MPI_calcu)+master_local_time,kind=Rkind)    &
              /Real(MPI_np,kind=Rkind)

      DO ii=1,MPI_np
        ! may consider to use inter Trapezoid
        current_i_MPI=0
        rest_time=0.
        DO i_MPI=0,MPI_np-2
          count_time=0
          IF(i_MPI==0) THEN
            rest_time=Real(time_MPI_calcu(i_MPI)+master_local_time,kind=Rkind)
            num_iG=MAX(1,iGs_MPI(2,i_MPI)-iGs_MPI(1,i_MPI)+1)
          ENDIF
          
          iG_threads(i_MPI)=0
          DO i_MPI2=current_i_MPI,MPI_np-1
            IF(rest_time+count_time<=ave_time) THEN
              count_time=count_time+rest_time
              iG_threads(i_MPI)=iG_threads(i_MPI)+num_iG

              ! go to next block
              IF(i_MPI2<=MPI_np-2) THEN
                rest_time=Real(time_MPI_calcu(i_MPI2+1),kind=Rkind)
                num_iG=MAX(1,iGs_MPI(2,i_MPI2+1)-iGs_MPI(1,i_MPI2+1)+1)
              ENDIF
            ELSE
              ! make sure there is at least 1 Smolyak term for each threads
              temp_int=MAX(1,FLOOR((ave_time-count_time)/ave_iGtime_threads(i_MPI2)))
              IF(i_MPI==0) temp_int=MAX(1,FLOOR((ave_time-count_time-master_local_time)&
                                   /ave_iGtime_threads(i_MPI2)))
              iG_threads(i_MPI)=iG_threads(i_MPI)+temp_int

              ! record the current block
              current_i_MPI=i_MPI2
              ! record the rest time and iGs in current block
              rest_time=Real(time_MPI_calcu(i_MPI2),kind=Rkind)-(ave_time-count_time)
              num_iG=MAX(1,iGs_MPI(2,i_MPI2)-iGs_MPI(1,i_MPI2)+1-temp_int)
              EXIT
            ENDIF
          ENDDO ! for i_MPI2=current_i_MPI,MPI_np-1
        ENDDO ! for i_MPI=0,MPI_np-2

        IF(SUM(iG_threads(0:MPI_np-2))>iG_total) THEN
          ave_time=ave_time*(1.0-1.0/MPI_np/10)
        ELSEIF(Real(SUM(iG_threads(0:MPI_np-2)))/Real(iG_total)<0.7) THEN
          ave_time=ave_time*(1.0+1.0/MPI_np/10)
        ELSE
          EXIT ! normal
        ENDIF

        IF(ii==MPI_np) STOP 'error in auto_iGs_MPI'
      ENDDO

      ! this make sure all iGs are accounted for
      iG_threads(MPI_np-1)=iG_total-SUM(iG_threads(0:MPI_np-2))

      !> set up new iGs arrangement
      iGs_MPI(1,0)=1
      iGs_MPI(2,0)=iGs_MPI(1,0)+iG_threads(0)-1
      DO i_MPI=1,MPI_np-1
        iGs_MPI(1,i_MPI)=iGs_MPI(2,i_MPI-1)+1
        iGs_MPI(2,i_MPI)=iGs_MPI(1,i_MPI)+iG_threads(i_MPI)-1
      ENDDO
    ENDIF ! for MPI_id==0

    ! share to all threads
    CALL MPI_Bcast_matrix(iGs_MPI(1:2,0:MPI_np-1),1,2,0,MPI_np-1,root_MPI,shift2=1)

#endif
  ENDSUBROUTINE auto_iGs_MPI_old5
!=======================================================================================

!=======================================================================================
!> auto adjust the distribution of iG on each threads according to 
!> time used in the previous action
!---------------------------------------------------------------------------------------
  SUBROUTINE auto_iGs_MPI(para_Op,iGs_change)
    USE mod_system
    USE mod_SetOp,ONLY:param_Op
    USE mod_MPI_aux
    IMPLICIT NONE
    
    TYPE(param_Op),                intent(in)    :: para_Op
    Logical,                       intent(inout) :: iGs_change
    
    Real                                         :: ave_iGtime_threads(0:MPI_np-1)
    Real                                         :: ave_time
    Real                                         :: rest_time
    Real                                         :: count_time
    Real                                         :: master_local_time
    Integer                                      :: time_MPI_real_comm(0:MPI_np-1)
    Integer                                      :: iG_threads(0:MPI_np-1)
    Integer                                      :: iG_total
    Integer                                      :: i_MPI2
    Integer                                      :: current_i_MPI
    Integer                                      :: num_iG
    Integer                                      :: ii
    Integer                                      :: jj
    Logical                                      :: returnall

#if(run_MPI)

    CALL MPI_collect_info(time_MPI_act_all) ! time_MPI_act_all infor on master
    CALL MPI_collect_info(time_MPI_calcu)
    CALL MPI_collect_info(time_MPI_local)

    iGs_change=.TRUE.
    IF(MPI_id==0) THEN 
      write(out_unitp,*) 'time used in action for each processer:',time_MPI_act_all
      write(out_unitp,*) 'local comm time used for each processer:',time_MPI_local
      write(out_unitp,*) 'calculation time used for each processer:',time_MPI_calcu
      ave_time=Real(SUM(time_MPI_act_all),kind=Rkind)/Real(MPI_np,kind=Rkind)
      IF(ALL(ABS(ave_time-time_MPI_act_all)/ave_time<0.1)) iGs_change=.FALSE.
      write(out_unitp,*) 'time balance:',ABS(ave_time-time_MPI_act_all)/ave_time<0.1
    ENDIF

    CALL MPI_BCAST(iGs_change,size1_MPI,MPI_Logical,root_MPI,MPI_COMM_WORLD,MPI_err)
    IF(.NOT. iGs_change) RETURN

    ! get the average time spend on each threads, 
    ! which is the ideal time for each threads
    ave_time=Real(SUM(time_MPI_local(0:MPI_np-1))+SUM(time_MPI_calcu(0:MPI_np-1)))     &
            /Real(MPI_np)

    IF(MPI_id==0) THEN
      iG_total=para_Op%BasisnD%para_SGType2%nb_SG
      
      ! get the average time for each iG on different threads
      DO i_MPI=0,MPI_np-1
        ave_iGtime_threads(i_MPI)=Real(time_MPI_calcu(i_MPI),kind=Rkind)               &
                                 /Real(iGs_MPI(2,i_MPI)-iGs_MPI(1,i_MPI)+1,kind=Rkind)
      ENDDO

      DO ii=1,MPI_np
        ! may consider to use inter Trapezoid
        current_i_MPI=0
        rest_time=0.
        DO i_MPI=0,MPI_np-2
          count_time=time_MPI_local(i_MPI)
          IF(i_MPI==0) THEN
            rest_time=Real(time_MPI_calcu(i_MPI),kind=Rkind)
            num_iG=MAX(1,iGs_MPI(2,i_MPI)-iGs_MPI(1,i_MPI)+1)
          ENDIF

          iG_threads(i_MPI)=0
          DO i_MPI2=current_i_MPI,MPI_np-1
            IF(rest_time+count_time<=ave_time) THEN
              count_time=count_time+rest_time
              iG_threads(i_MPI)=iG_threads(i_MPI)+num_iG

              ! go to next block
              IF(i_MPI2<=MPI_np-2) THEN
                rest_time=Real(time_MPI_calcu(i_MPI2+1),kind=Rkind)
                num_iG=MAX(1,iGs_MPI(2,i_MPI2+1)-iGs_MPI(1,i_MPI2+1)+1)
              ENDIF
            ELSE
              ! make sure there is at least 1 Smolyak term for each threads
              temp_int=MAX(1,FLOOR((ave_time-count_time)/ave_iGtime_threads(i_MPI2)))
              iG_threads(i_MPI)=iG_threads(i_MPI)+temp_int

              ! record the current block
              current_i_MPI=i_MPI2
              ! record the rest time and iGs in current block
              rest_time=Real(time_MPI_calcu(i_MPI2),kind=Rkind)-(ave_time-count_time)
              num_iG=MAX(1,iGs_MPI(2,i_MPI2)-iGs_MPI(1,i_MPI2)+1-temp_int)
              EXIT
            ENDIF
          ENDDO ! for i_MPI2=current_i_MPI,MPI_np-1
        ENDDO ! for i_MPI=0,MPI_np-2

        IF(SUM(iG_threads(0:MPI_np-2))>iG_total) THEN
          ave_time=ave_time*(1.0-1.0/MPI_np/10)
        ELSEIF(MPI_np>5 .AND. Real(SUM(iG_threads(0:MPI_np-2)))/Real(iG_total)<0.7) THEN
          ave_time=ave_time*(1.0+1.0/MPI_np/10)
        ELSE
          EXIT ! normal
        ENDIF

        IF(ii==MPI_np) THEN
          write(out_unitp,*) 'Warning: error in auto_iGs_MPI, back to initial'
          iGs_MPI=iGs_MPI0
          iGs_MPI_mc=iGs_MPI_mc0
          RETURN
        ENDIF
      ENDDO

      ! this make sure all iGs are accounted for
      iG_threads(MPI_np-1)=iG_total-SUM(iG_threads(0:MPI_np-2))

      !> set up new iGs arrangement
      iGs_MPI(1,0)=1
      iGs_MPI(2,0)=iGs_MPI(1,0)+iG_threads(0)-1
      DO i_MPI=1,MPI_np-1
        iGs_MPI(1,i_MPI)=iGs_MPI(2,i_MPI-1)+1
        iGs_MPI(2,i_MPI)=iGs_MPI(1,i_MPI)+iG_threads(i_MPI)-1
      ENDDO
    ENDIF ! for MPI_id==0

    ! share to all threads
    CALL MPI_Bcast_matrix(iGs_MPI(1:2,0:MPI_np-1),1,2,0,MPI_np-1,root_MPI,shift2=1)

#endif
  ENDSUBROUTINE auto_iGs_MPI
!=======================================================================================

!=======================================================================================
  SUBROUTINE get_size_ST(BasisnD,size_ST,size_ST_mc)
    USE mod_basis_set_alloc,ONLY:basis
    USE mod_MPI
    IMPLICIT NONE

    TYPE(basis),                              intent(in) :: BasisnD
    Integer(kind=MPI_INTEGER_KIND),pointer,intent(inout) :: size_ST(:)
    Integer(kind=MPI_INTEGER_KIND),pointer,intent(inout) :: size_ST_mc(:,:)

    Integer                                              :: iG
    Integer                                              :: ii
    Integer                                              :: temp_int

#if(run_MPI)

    Do i_MPI=0,MPI_np-1
      size_ST(i_MPI)=0
      DO iG=iGs_MPI(1,i_MPI),iGs_MPI(2,i_MPI)
        temp_int=BasisnD%para_SGType2%tab_nb_OF_SRep(iG)*BasisnD%para_SGType2%nb0
        size_ST(i_MPI)=size_ST(i_MPI)+temp_int
      ENDDO

      IF(MPI_mc>1) THEN
        DO ii=1,MPI_mc
          size_ST_mc(ii,i_MPI)=0
          DO iG=iGs_MPI_mc(1,ii,i_MPI),iGs_MPI_mc(2,ii,i_MPI)
            temp_int=BasisnD%para_SGType2%tab_nb_OF_SRep(iG)*BasisnD%para_SGType2%nb0
            size_ST_mc(ii,i_MPI)=size_ST_mc(ii,i_MPI)+temp_int
          ENDDO
        ENDDO
      ELSEIF(MPI_mc==1) THEN
        size_ST_mc(1,0:MPI_np-1)=size_ST(0:MPI_np-1)
      ELSE
        STOP 'MPI_mc error in get_size_ST' 
      ENDIF
    ENDDO ! for i_MPI=0,MPI_np-1

#endif
  ENDSUBROUTINE get_size_ST
!=======================================================================================

END MODULE mod_OpPsi_SG4_MPI
