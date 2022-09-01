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


!=======================================================================================
! subroutine not used anymore:
! Overlap_psi_Hpsi_matrix_MPI
! Overlap_psi_Hpsi_matrix_MPI2
!=======================================================================================

MODULE mod_psi_Op_MPI
  USE mod_basis
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: calculate_overlap_MPI,calculate_overlap1D_MPI
  PUBLIC :: distribute_psi_MPI
  PUBLIC :: Overlap_HS_matrix_MPI3,Overlap_H_matrix_MPI4
  PUBLIC :: Overlap_psi1_psi2_SRB_MPI,Overlap_psi1_psi2_SRG_MPI
  PUBLIC :: Set_symab_OF_psiBasisRep_MPI
  PUBLIC :: norm2_psi_SR_MPI,sub_LCpsi_TO_psi_MPI

  CONTAINS

!=======================================================================================
!> @brief MPI of Set_symab_OF_psiBasisRep when each threads contain part of vec
!> for symmetrization (with abelian group) of psi in BasisRep
!> for MPI_scheme 1 or the calculation of residual g
!
!  be careful with bounds_MPI. It should be consistant with the *vecB in psi
!---------------------------------------------------------------------------------------
  SUBROUTINE Set_symab_OF_psiBasisRep_MPI(psi,symab,changes)
    USE mod_system
    USE mod_psi_set_alloc
    USE mod_MPI_aux
    IMPLICIT NONE

    TYPE(param_psi),      intent(inout) :: psi
    Integer,optional,        intent(in) :: symab
    Logical,optional,     intent(inout) :: changes

    Integer                             :: loc_symab
    Integer                             :: ib
    Integer                             :: ib_temp
    Integer                             :: d1
    Integer                             :: d2
    Integer                             :: Get_symabOFSymAbelianOFBasis_AT_ib ! function

#if(run_MPI)

    IF(present(changes)) changes=.FALSE.

    d1=bounds_MPI(1,MPI_id)
    d2=bounds_MPI(2,MPI_id)

    IF(psi%BasisRep)THEN
      IF(psi%nb_bi==1 .AND. psi%nb_be==1) THEN
        IF(present(symab)) THEN
          loc_symab=symab
        ELSE
          ! find the symmtry (symab of the largest coef)
          IF (psi%cplx) THEN
            ib=maxloc(abs(psi%CvecB(d1:d2)),dim=1)
          ELSE
            ib=maxloc(abs(psi%RvecB(d1:d2)),dim=1)
          END IF

          CALL MPI_Reduce_max_Bcast(ib)

          loc_symab=Get_symabOFSymAbelianOFBasis_AT_ib(psi%BasisnD,ib)
          !write(out_unitp,*) 'maxloc,loc_symab check',ib,loc_symab,MPI_id
        ENDIF
      ELSE
        loc_symab=-1
      END IF
      psi%symab=loc_symab
    ELSE
      psi%symab=-1
    END IF

    IF(MPI_scheme==1 .OR. MPI_id==0) THEN
      d1=1
      d2=psi%nb_tot
    ENDIF

    IF(psi%symab>=0 .AND. psi%symab<=7) THEN
      IF(psi%cplx .AND. allocated(psi%CvecB)) THEN
        DO ib=d1,d2
          IF(psi%symab/=Get_symabOFSymAbelianOFBasis_AT_ib(psi%BasisnD,ib)) THEN
            IF(ABS(psi%CvecB(ib))>1.0E-20 .AND. present(changes)) changes=.TRUE.
            psi%CvecB(ib)=CZERO
          ENDIF
        ENDDO
      ELSEIF (.NOT. psi%cplx .AND. allocated(psi%RvecB)) THEN
        DO ib=d1,d2
          IF(psi%symab/=Get_symabOFSymAbelianOFBasis_AT_ib(psi%BasisnD,ib)) THEN
            IF(ABS(psi%RvecB(ib))>1.0E-20 .AND. present(changes)) changes=.TRUE.
            psi%RvecB(ib)=ZERO
          ENDIF
        ENDDO
      ENDIF
    ENDIF

    ! share information of changes
    IF(present(changes)) THEN
      CALL MPI_Reduce(changes,temp_logi,size1_MPI,MPI_Logical,MPI_LOR,                 &
                    root_MPI,MPI_COMM_WORLD,MPI_err)
      IF(MPI_id==0) changes=temp_logi
      CALL MPI_BCAST(changes,size1_MPI,MPI_Logical,root_MPI,MPI_COMM_WORLD,MPI_err)
    ENDIF

#endif
  END SUBROUTINE Set_symab_OF_psiBasisRep_MPI

!=======================================================================================
!> calculate overlap: <psi1|psi2> on Smolyak rep on Basis.
!---------------------------------------------------------------------------------------
  SUBROUTINE Overlap_psi1_psi2_SRB_MPI(Overlap,psi1,psi2)
    USE mod_system
    USE mod_psi_set_alloc
    USE mod_MPI_aux
    IMPLICIT NONE

    TYPE(param_psi),                intent(in)    :: psi1
    TYPE(param_psi),                intent(in)    :: psi2
    Complex(kind=Rkind),            intent(inout) :: Overlap

    Integer                                       :: iG
    Integer                                       :: d1
    Integer                                       :: d2

#if(run_MPI)

    Overlap=CMPLX(ZERO,ZERO,kind=Rkind)
    DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)
      d1=psi1%SR_B_index(iG)
      d2=psi1%SR_B_index(iG+1)-1
      IF(psi1%cplx) THEN
        Overlap=Overlap+psi1%BasisnD%WeightSG(iG)                                  &
              *dot_product(CMPLX(psi1%SR_B(d1:d2,1),psi1%SR_B(d1:d2,2),kind=Rkind),&
                           CMPLX(psi2%SR_B(d1:d2,1),psi2%SR_B(d1:d2,2),kind=Rkind))
      ELSE
        Overlap=Overlap+psi1%BasisnD%WeightSG(iG)                                  &
                *dot_product(psi1%SR_B(d1:d2,1),psi2%SR_B(d1:d2,1))
      ENDIF
    ENDDO

    !> reduce sum and boardcast
    CALL MPI_Reduce_sum_Bcast(Overlap)

#endif
  ENDSUBROUTINE Overlap_psi1_psi2_SRB_MPI
!=======================================================================================


!=======================================================================================
!> calculate overlap: <psi1|psi2> on Smolyak rep.
!---------------------------------------------------------------------------------------
  SUBROUTINE Overlap_psi1_psi2_SRG_MPI(Overlap,psi1,psi2)
    USE mod_system
    USE mod_psi_set_alloc
    USE mod_param_SGType2
    USE mod_MPI_aux
    USE mod_basis_BtoG_GtoB_SGType4,ONLY:getbis_tab_nq,getbis_tab_nb,              &
                                         GDP_TO_BDP_OF_SmolyakRep
    IMPLICIT NONE


    TYPE(param_psi),                intent(in)    :: psi1
    TYPE(param_psi),                intent(in)    :: psi2
    Complex(kind=Rkind),            intent(inout) :: Overlap
    
    TYPE(param_SGType2),pointer                   :: SGType2
    Real(kind=Rkind),allocatable                  :: SR_B1_R(:)
    Real(kind=Rkind),allocatable                  :: SR_B1_C(:)
    Real(kind=Rkind),allocatable                  :: SR_B2_R(:)
    Real(kind=Rkind),allocatable                  :: SR_B2_C(:)
    Integer,allocatable                           :: tab_nq(:)
    Integer,allocatable                           :: tab_nb(:)
    Integer                                       :: d1
    Integer                                       :: d2
    Integer                                       :: dim
    Integer                                       :: iG

#if(run_MPI)

    SGType2 => psi1%BasisnD%para_SGType2
    Overlap=CMPLX(ZERO,ZERO,kind=Rkind)

    DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)
      d1=psi1%SR_G_index(iG)
      d2=psi1%SR_G_index(iG+1)-1
      
      dim=size(SGType2%nDind_SmolyakRep%Tab_nDval(:,iG))
      IF(allocated(tab_nb)) deallocate(tab_nb)
      IF(allocated(tab_nq)) deallocate(tab_nq)
      allocate(tab_nb(dim))
      allocate(tab_nq(dim))
      tab_nq(:)=getbis_tab_nq(SGType2%nDind_SmolyakRep%Tab_nDval(:,iG),                &
                              psi1%BasisnD%tab_basisPrimSG)
      tab_nb(:)=getbis_tab_nb(SGType2%nDind_SmolyakRep%Tab_nDval(:,iG),                &
                              psi1%BasisnD%tab_basisPrimSG)

      IF(psi1%cplx) THEN
!            Overlap=Overlap+psi1%BasisnD%WeightSG(iG)                                  &
!                  *dot_product(CMPLX(psi1%SR_G(d1:d2,1),psi1%SR_G(d1:d2,2)),           &
!                               CMPLX(psi2%SR_G(d1:d2,1),psi2%SR_G(d1:d2,2)))
        CALL allocate_array(SR_B1_R,1,d2-d1+1)
        CALL allocate_array(SR_B1_C,1,d2-d1+1)
        CALL allocate_array(SR_B2_R,1,d2-d1+1)
        CALL allocate_array(SR_B2_C,1,d2-d1+1)
        SR_B1_R=psi1%SR_G(d1:d2,1)
        SR_B1_C=psi1%SR_G(d1:d2,2)
        SR_B2_R=psi2%SR_G(d1:d2,1)
        SR_B2_C=psi2%SR_G(d1:d2,2)

        CALL GDP_TO_BDP_OF_SmolyakRep(SR_B1_R,psi1%BasisnD%tab_basisPrimSG,            &
                 SGType2%nDind_SmolyakRep%Tab_nDval(:,iG),tab_nq,tab_nb,SGType2%nb0)
        CALL GDP_TO_BDP_OF_SmolyakRep(SR_B1_C,psi1%BasisnD%tab_basisPrimSG,            &
                 SGType2%nDind_SmolyakRep%Tab_nDval(:,iG),tab_nq,tab_nb,SGType2%nb0)
        CALL GDP_TO_BDP_OF_SmolyakRep(SR_B2_R,psi1%BasisnD%tab_basisPrimSG,            &
                 SGType2%nDind_SmolyakRep%Tab_nDval(:,iG),tab_nq,tab_nb,SGType2%nb0)
        CALL GDP_TO_BDP_OF_SmolyakRep(SR_B2_C,psi1%BasisnD%tab_basisPrimSG,            &
                 SGType2%nDind_SmolyakRep%Tab_nDval(:,iG),tab_nq,tab_nb,SGType2%nb0)

        Overlap=Overlap+psi1%BasisnD%WeightSG(iG)                                      &
                 *dot_product(CMPLX(SR_B1_R,SR_B1_C),CMPLX(SR_B2_R,SR_B2_C))
      ELSE
!            Overlap=Overlap+psi1%BasisnD%WeightSG(iG)                                  &
!                           *dot_product(psi1%SR_G(d1:d2,1),psi2%SR_G(d1:d2,1))
        CALL allocate_array(SR_B1_R,1,d2-d1+1)
        CALL allocate_array(SR_B2_R,1,d2-d1+1)
        SR_B1_R=psi1%SR_G(d1:d2,1)
        SR_B2_R=psi2%SR_G(d1:d2,1)
        
        CALL GDP_TO_BDP_OF_SmolyakRep(SR_B1_R,psi1%BasisnD%tab_basisPrimSG,            &
                 SGType2%nDind_SmolyakRep%Tab_nDval(:,iG),tab_nq,tab_nb,SGType2%nb0)
        CALL GDP_TO_BDP_OF_SmolyakRep(SR_B2_R,psi1%BasisnD%tab_basisPrimSG,            &
                 SGType2%nDind_SmolyakRep%Tab_nDval(:,iG),tab_nq,tab_nb,SGType2%nb0)

        Overlap=Overlap+psi1%BasisnD%WeightSG(iG)*dot_product(SR_B1_R,SR_B2_R)
      ENDIF
    ENDDO

    !> reduce sum and boardcast
    CALL MPI_Reduce_sum_Bcast(Overlap)

#endif
  ENDSUBROUTINE Overlap_psi1_psi2_SRG_MPI
!=======================================================================================


!=======================================================================================
!> normalize psi in Smolyak rep. 
!> scheme=1: calculate the normalization constant and normalize 
!> scheme=2: just calculate the normalization constant 
!> scheme=3: normalize with existing normalization constant
!
! improve later
!---------------------------------------------------------------------------------------
  SUBROUTINE norm2_psi_SR_MPI(psi,scheme)
    USE mod_system
    USE mod_psi_set_alloc
    IMPLICIT NONE
    
    TYPE(param_psi),                intent(inout) :: psi
    Integer,                        intent(in)    :: scheme

    Complex(kind=Rkind)                           :: Overlap

#if(run_MPI)

    !> calcualte the normalization constant 
    IF(scheme==1 .OR. scheme==2) THEN
      IF(psi%SRG_MPI) THEN
        CALL Overlap_psi1_psi2_SRG_MPI(Overlap,psi,psi)
      ELSE IF(psi%SRB_MPI) THEN
        CALL Overlap_psi1_psi2_SRB_MPI(Overlap,psi,psi)
      ENDIF
      psi%norm2=Real(Overlap,kind=Rkind)
    ENDIF

    IF(scheme==1 .OR. scheme==3) THEN
      IF(psi%SRG_MPI) THEN
        psi%SR_G=psi%SR_G/sqrt(psi%norm2)
      ELSE IF(psi%SRB_MPI) THEN
        psi%SR_B=psi%SR_B/sqrt(psi%norm2)
      ENDIF
      psi%norm2=1.0
    ENDIF

#endif
  ENDSUBROUTINE norm2_psi_SR_MPI
!=======================================================================================


!=======================================================================================
! subroutine for calculation of matrix H_overlap(i,j) for <psi(i)|Hpsi(j)> 
! MPI version, takes too much memory
!---------------------------------------------------------------------------------------
  SUBROUTINE Overlap_psi_Hpsi_matrix_MPI(H_overlap,S_overlap,psi,Hpsi,ndim,With_Grid)
    USE mod_system
    USE mod_psi_set_alloc
    USE mod_psi_Op,ONLY:Overlap_psi1_psi2
    USE mod_MPI_aux
    IMPLICIT NONE

    TYPE(param_psi), intent(inout)              :: psi(:)  !< inout only non-root threats
    TYPE(param_psi), intent(inout)              :: Hpsi(:) !< inout only non-root threats
    Logical,optional,intent(in)                 :: With_Grid
    Integer,         intent(in)                 :: ndim
    
    Real(kind=Rkind),allocatable,intent(inout)  :: H_overlap(:,:)
    Real(kind=Rkind),allocatable,intent(inout)  :: S_overlap(:,:)

    Character(len=*),parameter                  :: name_sub='Overlap_psi_Hpsi_matrix_MPI'
    !Real(kind=Rkind),allocatable                :: H_flat(:)
    !Real(kind=Rkind),allocatable                :: S_flat(:)
    Complex(kind=Rkind)                         :: Overlap
    Integer                                     :: num
    Integer                                     :: ii
    Integer                                     :: jj

#if(run_MPI)

    IF(allocated(H_overlap)) deallocate(H_overlap)
    IF(allocated(S_overlap)) deallocate(S_overlap)
    CALL alloc_NParray(H_overlap,[ndim,ndim],"H",name_sub)
    CALL alloc_NParray(S_overlap,[ndim,ndim],"S",name_sub)

    !allocate(H_flat(ndim*ndim))
    !allocate(S_flat(ndim*ndim))
    
    ! calculate on master without MPI
    IF(MPI_id==0) THEN
      num=0
      DO ii=1,ndim
        DO jj=1,ndim
          num=num+1
          CALL Overlap_psi1_psi2(Overlap,psi(jj),Hpsi(ii),With_Grid=With_Grid)
          H_overlap(jj,ii)=real(Overlap,kind=Rkind)
          
          CALL Overlap_psi1_psi2(Overlap,psi(jj), psi(ii),With_Grid=With_Grid)
          S_overlap(jj,ii)=real(Overlap,kind=Rkind)
          
          !H_flat(num)=H_overlap(jj,ii)
          !S_flat(num)=S_overlap(jj,ii)
        ENDDO
      ENDDO
    ENDIF ! for MPI_id==0

    CALL MPI_Bcast_matrix(H_overlap,1,ndim,1,ndim,root_MPI)
    CALL MPI_Bcast_matrix(S_overlap,1,ndim,1,ndim,root_MPI)

  !  CALL MPI_Bcast(H_flat,ndim*ndim,MPI_Real8,root_MPI,MPI_COMM_WORLD,MPI_err)
  !  CALL MPI_Bcast(S_flat,ndim*ndim,MPI_Real8,root_MPI,MPI_COMM_WORLD,MPI_err)
  !
  !  IF(MPI_id/=0) THEN
  !    num=0
  !    DO ii=1,ndim
  !      DO jj=1,ndim
  !        num=num+1
  !        H_overlap(jj,ii)=H_flat(num)
  !        S_overlap(jj,ii)=S_flat(num)
  !      ENDDO
  !    ENDDO
  !  ENDIF ! for MPI_id/=0
    
    ! Using MPI, takes too much memory'
    ! psi(1:ndim) are share on all threads, Hpsi(:) are distributed to each threads
    !CALL Overlap_psi1_psi2_MPI(H_overlap,S_overlap,psi,Hpsi,ndim,With_Grid)

#endif
  END SUBROUTINE Overlap_psi_Hpsi_matrix_MPI
!=======================================================================================


!=======================================================================================
! subroutine for calculation of matrix H_overlap(i,j) for <psi(i)|Hpsi(j)> 
! MPI version, less memory
!---------------------------------------------------------------------------------------
  SUBROUTINE Overlap_psi_Hpsi_matrix_MPI2(H_overlap,S_overlap,psi,Hpsi,ndim,With_Grid)
    USE mod_system
    USE mod_psi_set_alloc
    USE mod_MPI_aux
    IMPLICIT NONE

    TYPE(param_psi), intent(inout)              :: psi(:)  !< inout only non-root threats
    TYPE(param_psi), intent(inout)              :: Hpsi(:) !< inout only non-root threats
    Logical,optional,intent(in)                 :: With_Grid
    Integer,         intent(in)                 :: ndim

    Real(kind=Rkind),allocatable,intent(inout)  :: H_overlap(:,:)
    Real(kind=Rkind),allocatable,intent(inout)  :: S_overlap(:,:)

    Character(len=*),parameter                  :: name_sub='Overlap_psi_Hpsi_matrix_MPI2'
    Complex(kind=Rkind)                         :: Overlap
    Integer                                     :: ii
    Integer                                     :: jj

#if(run_MPI)

    IF(allocated(H_overlap)) deallocate(H_overlap)
    IF(allocated(S_overlap)) deallocate(S_overlap)
    CALL alloc_NParray(H_overlap,[ndim,ndim],"H",name_sub)
    CALL alloc_NParray(S_overlap,[ndim,ndim],"S",name_sub)

    !-------------------------------------------------------------------------------------
    ! calculate on master without MPI
  !  IF(MPI_id==0) THEN
  !    DO ii=1,ndim
  !      DO jj=1,ndim
  !        CALL Overlap_psi1_psi2(Overlap,psi(jj),Hpsi(ii),With_Grid=With_Grid)
  !        H_overlap(jj,ii)=real(Overlap,kind=Rkind)
  !
  !        CALL Overlap_psi1_psi2(Overlap,psi(jj), psi(ii),With_Grid=With_Grid)
  !        S_overlap(jj,ii)=real(Overlap,kind=Rkind)
  !
  !      ENDDO
  !    ENDDO
  !  ENDIF ! for MPI_id==0
  !  
  !  CALL MPI_Bcast_matrix(H_overlap,1,ndim,1,ndim,root_MPI)
  !  CALL MPI_Bcast_matrix(S_overlap,1,ndim,1,ndim,root_MPI)

    !-------------------------------------------------------------------------------------
    ! calculate with MPI
    CALL Overlap_psi1_psi2_MPI2(H_overlap,S_overlap,psi,Hpsi,ndim,With_Grid)

#endif
  END SUBROUTINE Overlap_psi_Hpsi_matrix_MPI2
!=======================================================================================

!=======================================================================================
!> Subroutine for the calculation of matrix H_overlap(i,j) for <psi(i)|Hpsi(j)> 
!>                                      and S_overlap(i,j) for <psi(i)| psi(j)> 
!> MPI V3

!> overlap is cut into submatrix for MPI parallel   
!> note, this requre new  Overlap_psipsi_MPI3 routine
!
!> the submatrix of psi and Hpsi are keeped also for the calculation 
!> of Residual g in MakeResidual_Davidson_MPI2
!---------------------------------------------------------------------------------------
  SUBROUTINE Overlap_HS_matrix_MPI3(H_overlap,S_overlap,psi,Hpsi,ndim0,ndim,With_Grid)
    USE mod_system
    USE mod_psi_set_alloc
    USE mod_MPI_aux
    IMPLICIT NONE
    
    TYPE(param_psi), intent(inout)              :: psi(:)  !< inout only non-root threats
    TYPE(param_psi), intent(inout)              :: Hpsi(:) !< inout only non-root threats
    Logical,optional,intent(in)                 :: With_Grid
    Integer,         intent(in)                 :: ndim0
    Integer,         intent(in)                 :: ndim
    
    Real(kind=Rkind),allocatable,intent(inout)  :: H_overlap(:,:)
    Real(kind=Rkind),allocatable,intent(inout)  :: S_overlap(:,:)
    
    Real(kind=Rkind),allocatable                :: H0_overlap(:,:)
    Real(kind=Rkind),allocatable                :: S0_overlap(:,:)
    
    Real(kind=Rkind),allocatable                :: H_overlapp(:,:)
    Real(kind=Rkind),allocatable                :: S_overlapp(:,:)

    Character(len=*),parameter                  :: name_sub='Overlap_HS_matrix_MPI3'
    Complex(kind=Rkind)                         :: Overlap
    Integer                                     :: ii
    Integer                                     :: jj

#if(run_MPI)

    CALL increase_martix(H_overlap,name_sub,ndim0,ndim)
    CALL increase_martix(S_overlap,name_sub,ndim0,ndim)

!  IF(allocated(H_overlap)) THEN
!!    CALL alloc_NParray(H0_overlap,[ndim,ndim],"H0",name_sub)
!!    H0_overlap(1:ndim0,1:ndim0)=H_overlap(1:ndim0,1:ndim0)
!!    deallocate(H_overlap)
!!    CALL alloc_NParray(H_overlap,[ndim,ndim],"H0",name_sub)
!!    H_overlap(1:ndim0,1:ndim0)=H0_overlap(1:ndim0,1:ndim0)
!    CALL alloc_NParray(H0_overlap,[ndim,ndim],"H0",name_sub)
!    H0_overlap(1:ndim0,1:ndim0)=H_overlap(1:ndim0,1:ndim0)
!    CALL move_alloc(H0_overlap,H_overlap) ! moves the allocation from H0_overlap to ...
!  ELSE
!    CALL alloc_NParray(H_overlap,[ndim,ndim],"H",name_sub)
!  ENDIF
!  
!  IF(allocated(S_overlap)) THEN
!!    CALL alloc_NParray(S0_overlap,[ndim,ndim],"S0",name_sub)
!!    S0_overlap(1:ndim0,1:ndim0)=S_overlap(1:ndim0,1:ndim0)
!!    deallocate(S_overlap)
!!    CALL alloc_NParray(S_overlap,[ndim,ndim],"S",name_sub)
!!    S_overlap(1:ndim0,1:ndim0)=S0_overlap(1:ndim0,1:ndim0)
!    CALL alloc_NParray(S0_overlap,[ndim,ndim],"S0",name_sub)
!    S0_overlap(1:ndim0,1:ndim0)=S_overlap(1:ndim0,1:ndim0)
!    CALL move_alloc(S0_overlap,S_overlap) ! moves the allocation from S0_overlap to ...
!  ELSE
!    CALL alloc_NParray(S_overlap,[ndim,ndim],"H",name_sub)
!  ENDIF

  !-------------------------------------------------------------------------------------
  ! calculate on master without MPI
!  CALL alloc_NParray(H_overlapp,[ndim,ndim],"H",name_sub)
!  CALL alloc_NParray(S_overlapp,[ndim,ndim],"S",name_sub)
!  IF(MPI_id==0) THEN
!    DO ii=1,ndim
!      DO jj=1,ndim
!        CALL Overlap_psi1_psi2(Overlap,psi(jj),Hpsi(ii),With_Grid=With_Grid)
!        H_overlapp(jj,ii)=real(Overlap,kind=Rkind)
!
!        CALL Overlap_psi1_psi2(Overlap,psi(jj), psi(ii),With_Grid=With_Grid)
!        S_overlapp(jj,ii)=real(Overlap,kind=Rkind)
!
!      ENDDO
!    ENDDO
!  ENDIF ! for MPI_id==0
!  
!  CALL MPI_Bcast_matrix(H_overlapp,1,ndim,1,ndim,root_MPI)
!  CALL MPI_Bcast_matrix(S_overlapp,1,ndim,1,ndim,root_MPI)

  !-------------------------------------------------------------------------------------
    ! calculate with MPI
    IF(MPI_scheme==1) THEN
      CALL calculate_overlap_S1_MPI(psi,1,ndim,With_Grid,Hpsi=Hpsi,                    &
                                    S_overlap=S_overlap,H_overlap=H_overlap)
    ELSE
      CALL Overlap_psi1_psi2_MPI5(H_overlap,S_overlap,psi,Hpsi,0,ndim,With_Grid)
    ENDIF
!  write(*,*) 'overlapp check',MAXVAL(ABS(H_overlapp-H_overlap)),                       &
!                              MAXVAL(ABS(S_overlapp-S_overlap)),                       &
!                              MAXVAL(ABS(H_overlap)),MAXVAL(ABS(S_overlap)),           &
!                              MAXVAL(ABS(H_overlapp)),MAXVAL(ABS(S_overlapp))

#endif
  END SUBROUTINE Overlap_HS_matrix_MPI3
!=======================================================================================

!=======================================================================================
!> MPI Subroutine for the calculation of matrix H_overlap(i,j) for <psi(i)|Hpsi(j)> 

!> overlap is cut into submatrix for MPI parallel   
!> note, this requre new Overlap_psipsi_MPI3 routine
!
!> the submatrix of psi have been distribed in sub_NewVec_Davidson
!> Hpsi are distribed and keeped also for the calculation
!> of Residual g in MakeResidual_Davidson_MPI3
!---------------------------------------------------------------------------------------
  SUBROUTINE Overlap_H_matrix_MPI4(H_overlap,psi,Hpsi,ndim0,ndim,With_Grid)
    USE mod_system
    USE mod_psi_set_alloc
    USE mod_MPI_aux
    IMPLICIT NONE
    
    TYPE(param_psi), intent(inout)              :: psi(:)  !< inout only non-root threats
    TYPE(param_psi), intent(inout)              :: Hpsi(:) !< inout only non-root threats
    Logical,optional,intent(in)                 :: With_Grid
    Integer,         intent(in)                 :: ndim0
    Integer,         intent(in)                 :: ndim
    
    Real(kind=Rkind),allocatable,intent(inout)  :: H_overlap(:,:)
    Real(kind=Rkind),allocatable                :: H0_overlap(:,:)
  !  Real(kind=Rkind),allocatable                :: H_overlapp(:,:)

    Character(len=*),parameter                  :: name_sub='Overlap_H_matrix_MPI4'
    Complex(kind=Rkind)                         :: Overlap
    Integer                                     :: ii
    Integer                                     :: jj

#if(run_MPI)

    CALL increase_martix(H_overlap,name_sub,ndim0,ndim)

!  IF(allocated(H_overlap)) THEN
!!    CALL alloc_NParray(H0_overlap,[ndim,ndim],"H0",name_sub)
!!    H0_overlap(1:ndim0,1:ndim0)=H_overlap(1:ndim0,1:ndim0)
!!    deallocate(H_overlap)
!!    CALL alloc_NParray(H_overlap,[ndim,ndim],"H0",name_sub)
!!    H_overlap(1:ndim0,1:ndim0)=H0_overlap(1:ndim0,1:ndim0)
!    CALL alloc_NParray(H0_overlap,[ndim,ndim],"H0",name_sub)
!    H0_overlap(1:ndim0,1:ndim0)=H_overlap(1:ndim0,1:ndim0)
!    CALL move_alloc(H0_overlap,H_overlap) ! moves the allocation from H0_overlap to ...
!  ELSE
!    CALL alloc_NParray(H_overlap,[ndim,ndim],"H",name_sub)
!  ENDIF

  !-------------------------------------------------------------------------------------
  ! calculate on master without MPI
!  CALL alloc_NParray(H_overlapp,[ndim,ndim],"H",name_sub)
!  IF(MPI_id==0) THEN
!    DO ii=1,ndim
!      DO jj=1,ndim
!        CALL Overlap_psi1_psi2(Overlap,psi(jj),Hpsi(ii),With_Grid=With_Grid)
!        H_overlapp(jj,ii)=real(Overlap,kind=Rkind)
!      ENDDO
!    ENDDO
!  ENDIF ! for MPI_id==0
!  
!  CALL MPI_Bcast_matrix(H_overlapp,1,ndim,1,ndim,root_MPI)

    !-----------------------------------------------------------------------------------
    ! calculate with MPI
    IF(MPI_scheme==1) THEN
      CALL calculate_overlap_S1_MPI(psi,ndim0+1,ndim,With_Grid=With_Grid,Hpsi=Hpsi,    &
                                    H_overlap=H_overlap)
    ELSE
      CALL Overlap_psi_Hpsi_MPI(H_overlap,psi,Hpsi,ndim0,ndim,With_Grid)
    ENDIF
!  write(*,*) 'H_overlapp check',MAXVAL(ABS(H_overlapp-H_overlap)),    &
!                                MAXVAL(ABS(H_overlap)),MAXVAL(ABS(H_overlapp))

#endif
  END SUBROUTINE Overlap_H_matrix_MPI4

!=======================================================================================

!=======================================================================================
!> MPI Subroutine for the calculation of matrix S_overlap(i,j) for <psi(i)|psi(j)> 

!> overlap is cut into submatrix for MPI parallel   
!> note, this requre new Overlap_psipsi_MPI3 routine
!
!> psi are distribed and keeped also for the calculation
!> of H_overlap in sub_MakeHPsi_Davidson 
!> and Residual g in MakeResidual_Davidson_MPI3
!---------------------------------------------------------------------------------------
  SUBROUTINE Overlap_S_matrix_MPI4(S_overlap,psi,ndim0,ndim,With_Grid)
    USE mod_system
    USE mod_psi_set_alloc
    USE mod_MPI_aux
    IMPLICIT NONE
    
    TYPE(param_psi), intent(inout)              :: psi(:)  !< inout only non-root threats
    Logical,optional,intent(in)                 :: With_Grid
    Integer,         intent(in)                 :: ndim0
    Integer,         intent(in)                 :: ndim
    
    Real(kind=Rkind),allocatable,intent(inout)  :: S_overlap(:,:)
    Real(kind=Rkind),allocatable                :: S0_overlap(:,:)
    !Real(kind=Rkind),allocatable                :: S_overlapp(:,:)

    Character(len=*),parameter                  :: name_sub='Overlap_S_matrix_MPI4'
    Complex(kind=Rkind)                         :: Overlap
    Integer                                     :: ii
    Integer                                     :: jj

#if(run_MPI)

    CALL increase_martix(S_overlap,name_sub,ndim0,ndim)
    
!  IF(allocated(S_overlap)) THEN
!!    CALL alloc_NParray(S0_overlap,[ndim,ndim],"S0",name_sub)
!!    S0_overlap(1:ndim0,1:ndim0)=S_overlap(1:ndim0,1:ndim0)
!!    deallocate(S_overlap)
!!    CALL alloc_NParray(S_overlap,[ndim,ndim],"S",name_sub)
!!    S_overlap(1:ndim0,1:ndim0)=S0_overlap(1:ndim0,1:ndim0)
!    CALL alloc_NParray(S0_overlap,[ndim,ndim],"S0",name_sub)
!    S0_overlap(1:ndim0,1:ndim0)=S_overlap(1:ndim0,1:ndim0)
!    CALL move_alloc(S0_overlap,S_overlap) ! moves the allocation from S0_overlap to ...
!  ELSE
!    CALL alloc_NParray(S_overlap,[ndim,ndim],"H",name_sub)
!  ENDIF
  
  !-------------------------------------------------------------------------------------
  ! calculate on master without MPI
!  CALL alloc_NParray(S_overlapp,[ndim,ndim],"S",name_sub)
!  IF(MPI_id==0) THEN
!    DO ii=1,ndim
!      DO jj=1,ndim
!        CALL Overlap_psi1_psi2(Overlap,psi(jj), psi(ii),With_Grid=With_Grid)
!        S_overlapp(jj,ii)=real(Overlap,kind=Rkind)
!      ENDDO
!    ENDDO
!  ENDIF ! for MPI_id==0
!  
!  CALL MPI_Bcast_matrix(S_overlapp,1,ndim,1,ndim,root_MPI)

    !-----------------------------------------------------------------------------------
    ! calculate with MPI
    !CALL Overlap_psi_psi_MPI(S_overlap,psi,ndim0,ndim,With_Grid)
    ! no distribution  of psi
    CALL calculate_overlap_MPI(psi,ndim0+1,ndim,With_Grid=With_Grid,S_overlap=S_overlap)

!  write(*,*) 'S_overlapp check',MAXVAL(ABS(S_overlapp-S_overlap)),    &
!                                MAXVAL(ABS(S_overlap)),MAXVAL(ABS(S_overlapp))

#endif
  END SUBROUTINE Overlap_S_matrix_MPI4
!=======================================================================================

!=======================================================================================
!> Overlap_S_matrix_MPI without the distribution of psi
!---------------------------------------------------------------------------------------
  SUBROUTINE Overlap_S0_matrix_MPI4(S_overlap,psi,ndim0,ndim,With_Grid)
    USE mod_system
    USE mod_psi_set_alloc
    USE mod_MPI_aux
    IMPLICIT NONE

    TYPE(param_psi), intent(inout)              :: psi(:)  !< inout only non-root threats
    Real(kind=Rkind),allocatable,intent(inout)  :: S_overlap(:,:)
    Logical,optional,intent(in)                 :: With_Grid
    Integer,         intent(in)                 :: ndim0
    Integer,         intent(in)                 :: ndim
    
    Character(len=*),parameter                  :: name_sub='Overlap_S_matrix_MPI4'
    Complex(kind=Rkind)                         :: Overlap
    Integer                                     :: ii
    Integer                                     :: jj

#if(run_MPI)

    CALL increase_martix(S_overlap,name_sub,ndim0,ndim)
    
    ! calculate without the distribution of psi
    CALL calculate_overlap_MPI(psi,ndim0+1,ndim,With_Grid=With_Grid,S_overlap=S_overlap)

#endif
  END SUBROUTINE Overlap_S0_matrix_MPI4
!=======================================================================================

!=======================================================================================
!> subroutine for the calculation of Overlap_psi_Hpsi with MPI 
!> overlap is cut into submatrix for MPI parallel   
!> note, this requre new Overlap_psi1_psi2 routine
!
!> the submatrix of Hpsi are keeped also for the calculation 
!> of Residual g in MakeResidual_Davidson_MPI2
!---------------------------------------------------------------------------------------
  SUBROUTINE Overlap_psi_Hpsi_MPI(H_overlap,psi,Hpsi,ndim0,ndim,With_Grid)
    USE mod_system
    USE mod_psi_set_alloc
    USE mod_MPI_aux
    IMPLICIT NONE

    TYPE(param_psi), intent(inout)          :: psi(:)  !< on non-root threads allocated
    TYPE(param_psi), intent(inout)          :: Hpsi(:) !< on non-root threads allocated
    Real(kind=Rkind),intent(inout)          :: H_overlap(:,:)
    Logical,optional,intent(in)             :: With_Grid
    Integer,         intent(in)             :: ndim0
    Integer,         intent(in)             :: ndim

    complex(kind=Rkind)                     :: Overlap
    logical                                 :: With_Grid_loc
    Integer                                 :: i
    Integer                                 :: j
    Logical                                 :: root_jobs

#if(run_MPI)

    With_Grid_loc=.FALSE.
    IF (present(With_Grid)) With_Grid_loc = With_Grid

    ! only on master: psi%CvecG,psi%RvecG,psi%CvecB,psi%RvecB
    CALL distribute_psi_pack_MPI(Hpsi,ndim0+1,ndim,With_Grid=With_Grid_loc)
    CALL calculate_overlap_MPI(psi,ndim0+1,ndim,With_Grid=With_Grid_loc,Hpsi=Hpsi,       &
                               H_overlap=H_overlap)

#endif
  END SUBROUTINE Overlap_psi_Hpsi_MPI
!=======================================================================================

!=======================================================================================
!> subroutine for the calculation of S_overlap(Overlap_psi_psi) with MPI 
!> overlap is cut into submatrix for MPI parallel   
!> note, this requre new Overlap_psi1_psi2 routine
!
!> the submatrix of psi are keeped also for the calculation 
!> of H_overlap and Residual g in MakeResidual_Davidson_MPI3
!---------------------------------------------------------------------------------------
  SUBROUTINE Overlap_psi_psi_MPI(S_overlap,psi,ndim0,ndim,With_Grid)
    USE mod_system
    USE mod_psi_set_alloc
    USE mod_MPI_aux
    IMPLICIT NONE

    TYPE(param_psi), intent(inout)          :: psi(:)  !< on non-root threads allocated
    Real(kind=Rkind),intent(inout)          :: S_overlap(:,:)
    Logical,optional,intent(in)             :: With_Grid
    Integer,         intent(in)             :: ndim0
    Integer,         intent(in)             :: ndim

    logical                                 :: With_Grid_loc
    Integer                                 :: i,j
    Logical                                 :: root_jobs

#if(run_MPI)

    With_Grid_loc=.FALSE.
    IF(present(With_Grid)) With_Grid_loc=With_Grid

    ! only on master: psi%CvecG,psi%RvecG,psi%CvecB,psi%RvecB
    CALL distribute_psi_pack_MPI(psi,ndim0+1,ndim,With_Grid=With_Grid_loc)
    CALL calculate_overlap_MPI(psi,ndim0+1,ndim,With_Grid=With_Grid_loc,                 &
                               S_overlap=S_overlap)

#endif
  END SUBROUTINE Overlap_psi_psi_MPI
!=======================================================================================

!=======================================================================================
!> subroutine for distribute psi(ndim1:ndim2) to different threads 
!> the distribution is depends on the length of vec (RvecB, CvecB, RvecG, or CvecG)
!---------------------------------------------------------------------------------------
  SUBROUTINE distribute_psi_MPI(psi,ndim1,ndim2,With_Grid)
    USE mod_system
    USE mod_psi_set_alloc
    USE mod_MPI_aux
    IMPLICIT NONE
    
    TYPE(param_psi), intent(inout)          :: psi(:) 
    Integer,         intent(in)             :: ndim1
    Integer,         intent(in)             :: ndim2
    Logical,optional,intent(in)             :: With_Grid
   
    Integer                                 :: i,j
    Logical                                 :: With_Grid_loc

#if(run_MPI)

    With_Grid_loc=.FALSE.
    IF (present(With_Grid)) With_Grid_loc=With_Grid

    IF(MPI_id==0) THEN
      DO i=ndim1,ndim2
        IF(.NOT. With_Grid_loc) THEN
          nb_per_MPI=psi(i)%nb_tot/MPI_np
          nb_rem_MPI=mod(psi(i)%nb_tot,MPI_np) !remainder jobs
        ELSE
          nb_per_MPI=psi(i)%nb_qaie/MPI_np
          nb_rem_MPI=mod(psi(i)%nb_qaie,MPI_np) !remainder jobs
        ENDIF
        
        ! Send array
        DO i_MPI=1,MPI_np-1
          bound1_MPI=i_MPI*nb_per_MPI+1+MIN(i_MPI,nb_rem_MPI)
          bound2_MPI=(i_MPI+1)*nb_per_MPI+MIN(i_MPI,nb_rem_MPI)                          &
                                         +merge(1,0,nb_rem_MPI>i_MPI)
          IF (.NOT. With_Grid_loc) THEN
            IF (psi(i)%cplx) THEN
              CALL MPI_Send(psi(i)%CvecB(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1, &
                            MPI_Complex8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
            ELSE
              CALL MPI_Send(psi(i)%RvecB(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1, &
                            MPI_Real8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
            ENDIF
          ELSE
            IF (psi(i)%cplx) THEN
              CALL MPI_Send(psi(i)%CvecG(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1, &
                            MPI_Complex8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)        
            ELSE
              CALL MPI_Send(psi(i)%RvecG(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1, &
                            MPI_Real8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
            ENDIF
          ENDIF ! for .NOT. With_Grid_loc 
        ENDDO ! for i_MPI=1,MPI_np-1
      ENDDO ! for i=ndim1,ndim2
    ENDIF ! for MPI_id==0
    
    ! MPI/=0------------------------------------------------------------------------------
    IF(MPI_id/=0) THEN
      DO i=ndim1,ndim2
        IF(.NOT. With_Grid_loc) THEN
          nb_per_MPI=psi(i)%nb_tot/MPI_np
          nb_rem_MPI=mod(psi(i)%nb_tot,MPI_np) !remainder jobs
        ELSE
          nb_per_MPI=psi(i)%nb_qaie/MPI_np
          nb_rem_MPI=mod(psi(i)%nb_qaie,MPI_np) !remainder jobs
        ENDIF
        bound1_MPI=MPI_id*nb_per_MPI+1+MIN(MPI_id,nb_rem_MPI)
        bound2_MPI=(MPI_id+1)*nb_per_MPI+MIN(MPI_id,nb_rem_MPI)                          &
                                        +merge(1,0,nb_rem_MPI>MPI_id)
        
        ! receive array
        IF (.NOT. With_Grid_loc) THEN
          IF (psi(i)%cplx) THEN
            IF(.NOT. allocated(psi(i)%CvecB))                                            &
                 allocate(psi(i)%CvecB(bound1_MPI:bound2_MPI))
            CALL MPI_Recv(psi(i)%CvecB(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,   &
                          MPI_Complex8,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
          ELSE
            IF(.NOT. allocated(psi(i)%RvecB))                                            &
                 allocate(psi(i)%RvecB(bound1_MPI:bound2_MPI))
            CALL MPI_Recv( psi(i)%RvecB(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,  &
                          MPI_Real8,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
          ENDIF
        ELSE
          IF (psi(i)%cplx) THEN
            IF(.NOT. allocated(psi(i)%CvecG))                                            &
                 allocate(psi(i)%CvecG(bound1_MPI:bound2_MPI))
            CALL MPI_Recv(psi(i)%CvecG(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,   &
                          MPI_Complex8,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
          ELSE
            IF(.NOT. allocated(psi(i)%RvecG))                                            &
                 allocate(psi(i)%RvecG(bound1_MPI:bound2_MPI))
            CALL MPI_Recv(psi(i)%RvecG(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,   &
                          MPI_Real8,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
          ENDIF
        ENDIF
      ENDDO ! for i=ndim1,ndim2
    ENDIF ! for MPI_id/=0  

#endif
  END SUBROUTINE distribute_psi_MPI
!=======================================================================================

!=======================================================================================
!> subroutine for distribute psi(ndim1:ndim2) to different threads 
!> the distribution is depends on the length of vec (RvecB, CvecB, RvecG, or CvecG)
!> vector is packed on master and unpacked on theards to reduce comm. time.
!---------------------------------------------------------------------------------------
  SUBROUTINE distribute_psi_pack_MPI(psi,ndim1,ndim2,With_Grid)
    USE mod_system
    USE mod_psi_set_alloc
    USE mod_MPI_aux
    IMPLICIT NONE
    
    TYPE(param_psi), intent(inout)          :: psi(:) 
    Integer,         intent(in)             :: ndim1
    Integer,         intent(in)             :: ndim2
    Logical,optional,intent(in)             :: With_Grid
   
    Complex(kind=Rkind),allocatable         :: Cvec(:)
    Real(kind=Rkind),allocatable            :: Rvec(:)
    Integer                                 :: length(MPI_np)
    Integer                                 :: count
    Integer                                 :: i,j
    Logical                                 :: With_Grid_loc

#if(run_MPI)

    With_Grid_loc=.FALSE.
    IF (present(With_Grid)) With_Grid_loc=With_Grid

    DO i_MPI=1,MPI_np-1
      length(i_MPI)=0
      DO i=ndim1,ndim2
        IF(.NOT. With_Grid_loc) THEN
          nb_per_MPI=psi(i)%nb_tot/MPI_np
          nb_rem_MPI=mod(psi(i)%nb_tot,MPI_np) 
        ELSE
          nb_per_MPI=psi(i)%nb_qaie/MPI_np
          nb_rem_MPI=mod(psi(i)%nb_qaie,MPI_np) 
        ENDIF
        bound1_MPI=i_MPI*nb_per_MPI+1+MIN(i_MPI,nb_rem_MPI)
        bound2_MPI=(i_MPI+1)*nb_per_MPI+MIN(i_MPI,nb_rem_MPI)                          &
                   +merge(1,0,nb_rem_MPI>i_MPI)

        length(i_MPI)=length(i_MPI)+(bound2_MPI-bound1_MPI+1)
      ENDDO
    ENDDO

    ! IF(MPI_id==0) THEN
    !   DO i_MPI=1,MPI_np-1
    IF(MPI_id==0 .OR. MPI_nodes_p0) THEN ! for S2 & S3
      DO i_MPI=MPI_sub_id(1),MPI_sub_id(2)

        IF(psi(1)%cplx) THEN
          allocate(Cvec(length(i_MPI)))
        ELSE
          allocate(Rvec(length(i_MPI)))
        ENDIF
        count=1

        DO i=ndim1,ndim2
          IF(.NOT. With_Grid_loc) THEN
            nb_per_MPI=psi(i)%nb_tot/MPI_np
            nb_rem_MPI=mod(psi(i)%nb_tot,MPI_np) !remainder jobs
          ELSE
            nb_per_MPI=psi(i)%nb_qaie/MPI_np
            nb_rem_MPI=mod(psi(i)%nb_qaie,MPI_np) !remainder jobs
          ENDIF
          bound1_MPI=i_MPI*nb_per_MPI+1+MIN(i_MPI,nb_rem_MPI)
          bound2_MPI=(i_MPI+1)*nb_per_MPI+MIN(i_MPI,nb_rem_MPI)                        &
                                         +merge(1,0,nb_rem_MPI>i_MPI)

          ! pack array
          IF (.NOT. With_Grid_loc) THEN
            IF (psi(i)%cplx) THEN
              Cvec(count:bound2_MPI-bound1_MPI+count)=psi(i)%CvecB(bound1_MPI:bound2_MPI)
            ELSE
              Rvec(count:bound2_MPI-bound1_MPI+count)=psi(i)%RvecB(bound1_MPI:bound2_MPI)
            ENDIF
          ELSE
            IF (psi(i)%cplx) THEN
              Cvec(count:bound2_MPI-bound1_MPI+count)=psi(i)%CvecG(bound1_MPI:bound2_MPI)
            ELSE
              Rvec(count:bound2_MPI-bound1_MPI+count)=psi(i)%RvecG(bound1_MPI:bound2_MPI)
            ENDIF
          ENDIF ! for .NOT. With_Grid_loc 
          count=bound2_MPI-bound1_MPI+count+1
          
        ENDDO ! for i=ndim1,ndim2 

        ! send array
        IF(psi(1)%cplx) THEN
          CALL MPI_Send(Cvec,count-1,MPI_Complex8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
          deallocate(Cvec)
        ELSE
          CALL MPI_Send(Rvec,count-1,MPI_Real8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
          deallocate(Rvec)
        ENDIF
      ENDDO ! for i_MPI
    ENDIF ! for MPI_id==0 .OR. MPI_nodes_p0
    
    !-----------------------------------------------------------------------------------
    ! IF(MPI_id/=0) THEN
    IF(.NOT.(MPI_id==0 .OR. MPI_nodes_p0)) THEN ! for S2 & S3
      ! receive array
      IF(psi(1)%cplx) THEN
        allocate(Cvec(length(MPI_id)))
        ! CALL MPI_Recv(Cvec,length(MPI_id),MPI_Complex8,root_MPI,MPI_id,                &
        !               MPI_COMM_WORLD,MPI_stat,MPI_err)
        CALL MPI_Recv(Cvec,length(MPI_id),MPI_Complex8,MPI_node_p0_id,MPI_id,          &
                      MPI_COMM_WORLD,MPI_stat,MPI_err)
      ELSE
        allocate(Rvec(length(MPI_id)))
        ! CALL MPI_Recv(Rvec,length(MPI_id),MPI_Real8,root_MPI,MPI_id,                   &
        !               MPI_COMM_WORLD,MPI_stat,MPI_err)
        CALL MPI_Recv(Rvec,length(MPI_id),MPI_Real8,MPI_node_p0_id,MPI_id,             &
                      MPI_COMM_WORLD,MPI_stat,MPI_err)
      ENDIF

      count=1
      DO i=ndim1,ndim2
        IF(.NOT. With_Grid_loc) THEN
          nb_per_MPI=psi(i)%nb_tot/MPI_np
          nb_rem_MPI=mod(psi(i)%nb_tot,MPI_np) !remainder jobs
        ELSE
          nb_per_MPI=psi(i)%nb_qaie/MPI_np
          nb_rem_MPI=mod(psi(i)%nb_qaie,MPI_np) !remainder jobs
        ENDIF
        bound1_MPI=MPI_id*nb_per_MPI+1+MIN(MPI_id,nb_rem_MPI)
        bound2_MPI=(MPI_id+1)*nb_per_MPI+MIN(MPI_id,nb_rem_MPI)                        &
                                        +merge(1,0,nb_rem_MPI>MPI_id)

        ! unpack
        IF (.NOT. With_Grid_loc) THEN
          IF (psi(i)%cplx) THEN
            IF(.NOT. allocated(psi(i)%CvecB))                                          &
                      allocate(psi(i)%CvecB(bound1_MPI:bound2_MPI))
            psi(i)%CvecB(bound1_MPI:bound2_MPI)=Cvec(count:count+bound2_MPI-bound1_MPI)
          ELSE
            IF(.NOT. allocated(psi(i)%RvecB))                                          &
                      allocate(psi(i)%RvecB(bound1_MPI:bound2_MPI))
            psi(i)%RvecB(bound1_MPI:bound2_MPI)=Rvec(count:count+bound2_MPI-bound1_MPI)
          ENDIF
        ELSE
          IF (psi(i)%cplx) THEN
            IF(.NOT. allocated(psi(i)%CvecG))                                          &
                      allocate(psi(i)%CvecG(bound1_MPI:bound2_MPI))
            psi(i)%CvecG(bound1_MPI:bound2_MPI)=Cvec(count:count+bound2_MPI-bound1_MPI)
          ELSE
            IF(.NOT. allocated(psi(i)%RvecG))                                          &
                      allocate(psi(i)%RvecG(bound1_MPI:bound2_MPI))
            psi(i)%RvecG(bound1_MPI:bound2_MPI)=Rvec(count:count+bound2_MPI-bound1_MPI)
          ENDIF
        ENDIF
        count=bound2_MPI-bound1_MPI+count+1
      ENDDO ! for i=ndim1,ndim2
      
      IF(psi(1)%cplx) THEN
        deallocate(Cvec)
      ELSE
        deallocate(Rvec)
      ENDIF
    ENDIF ! for .NOT.(MPI_id==0 .OR. MPI_nodes_p0)

#endif
  END SUBROUTINE distribute_psi_pack_MPI
!=======================================================================================

!=======================================================================================
!> @brief calculate the overlap of <psi|psi> or <psi|Hpsi> with MPI
!> working for ndim1=<i<=ndim2 
!
!>   NOTE: the pre-distribution of vectors (RvecB, CvecB, RvecG, or CvecG) are required
!>         i.e. subroutine distribute_psi_MPI
!>   NOTE: works according to the vectors stored on different threads 
!>         limited by bound1_MPI and bound2_MPI, see "Overlap_psipsi_MPI3"
!---------------------------------------------------------------------------------------
  SUBROUTINE calculate_overlap_MPI(psi,ndim1,ndim2,With_Grid,Hpsi,S_overlap,H_overlap)
    USE mod_system
    USE mod_psi_set_alloc
    USE mod_MPI_aux
    IMPLICIT NONE
    
    TYPE(param_psi),          intent(in)    :: psi(:) 
    TYPE(param_psi),optional, intent(in)    :: Hpsi(:) 
    Real(kind=Rkind),optional,intent(inout) :: H_overlap(:,:)
    Real(kind=Rkind),optional,intent(inout) :: S_overlap(:,:)
    Logical,optional,         intent(in)    :: With_Grid
    Integer,                  intent(in)    :: ndim1
    Integer,                  intent(in)    :: ndim2
    
    Complex(kind=Rkind)                     :: Overlap
    Integer                                 :: i,j
    Logical                                 :: With_Grid_loc

#if(run_MPI)

    With_Grid_loc=.FALSE.
    IF(present(With_Grid)) With_Grid_loc=With_Grid

    IF((present(Hpsi) .AND. (.NOT. present(H_overlap))) .OR.                           &
       ((.NOT. present(H_overlap)) .AND. (.NOT. present(S_overlap)))) THEN
      STOP 'variable presented error in calculate_overlap_MPI' 
    ENDIF

    !-overlap <psi|Hpsi>----------------------------------------------------------------
    IF(present(Hpsi) .AND. present(H_overlap)) THEN 
      DO i=1,ndim1-1
        DO j=ndim1,ndim2
          CALL Overlap_psipsi_MPI3(Overlap,psi(j),Hpsi(i),With_Grid=With_Grid)
          H_overlap(j,i)=real(Overlap,kind=Rkind) 
        ENDDO 
      END DO
      
      DO i=ndim1,ndim2
        DO j=1,ndim2
          CALL Overlap_psipsi_MPI3(Overlap,psi(j),Hpsi(i),With_Grid=With_Grid)
          H_overlap(j,i)=real(Overlap,kind=Rkind) 
        ENDDO 
      ENDDO

      ! collect and broadcast result
      CALL MPI_Reduce_sum_matrix(H_overlap,ndim1,ndim2,1,ndim1-1,root_MPI)
      CALL MPI_Bcast_matrix     (H_overlap,ndim1,ndim2,1,ndim1-1,root_MPI)
      CALL MPI_Reduce_sum_matrix(H_overlap,1,ndim2,ndim1,ndim2,root_MPI)
      CALL MPI_Bcast_matrix     (H_overlap,1,ndim2,ndim1,ndim2,root_MPI)
    ENDIF 

    !-overlap <psi|psi>-----------------------------------------------------------------
    ! calculate half of the matrix. Note the conjucate for complex case 
    IF(present(S_overlap)) THEN 
      DO i=1,ndim1-1
        DO j=ndim1,ndim2
          IF(j<=i) THEN
            CALL Overlap_psipsi_MPI3(Overlap,psi(j),psi(i),With_Grid=With_Grid)
            S_overlap(j,i)=real(Overlap,kind=Rkind)
            S_overlap(i,j)=S_overlap(j,i)
          ENDIF
        ENDDO 
      ENDDO 
      
      DO i=ndim1,ndim2
        DO j=1,ndim2
          IF(j<=i) THEN
            CALL Overlap_psipsi_MPI3(Overlap,psi(j),psi(i),With_Grid=With_Grid)
            S_overlap(j,i)=real(Overlap,kind=Rkind)
            S_overlap(i,j)=S_overlap(j,i)
          ENDIF
        ENDDO 
      ENDDO

      ! collect and broadcast result 
      CALL MPI_Reduce_sum_matrix(S_overlap,ndim1,ndim2,1,ndim1-1,root_MPI)
      CALL MPI_Bcast_matrix     (S_overlap,ndim1,ndim2,1,ndim1-1,root_MPI)
      CALL MPI_Reduce_sum_matrix(S_overlap,1,ndim2,ndim1,ndim2,root_MPI)
      CALL MPI_Bcast_matrix     (S_overlap,1,ndim2,ndim1,ndim2,root_MPI)

    ENDIF

#endif
  END SUBROUTINE calculate_overlap_MPI
!=======================================================================================


!=======================================================================================
!> @brief calculate <psi|psi> and/or <psi|H|psi> in MPI scheme 1 
!---------------------------------------------------------------------------------------
  SUBROUTINE calculate_overlap_S1_MPI(psi,ndim1,ndim2,With_Grid,Hpsi,S_overlap,H_overlap)
    USE mod_system
    USE mod_psi_set_alloc
    USE mod_psi_Op,ONLY:Overlap_psi1_psi2
    USE mod_MPI_aux
    IMPLICIT NONE

    TYPE(param_psi),          intent(inout) :: psi(:)
    TYPE(param_psi),optional, intent(inout) :: Hpsi(:)
    Real(kind=Rkind),optional,intent(inout) :: H_overlap(:,:)
    Real(kind=Rkind),optional,intent(inout) :: S_overlap(:,:)
    Logical,optional,            intent(in) :: With_Grid
    Integer,                     intent(in) :: ndim1
    Integer,                     intent(in) :: ndim2

    Complex(kind=Rkind)                     :: Overlap
    Real(kind=Rkind),allocatable            :: vec(:)
    Integer                                 :: ii
    Integer                                 :: jj
    Integer                                 :: nn
    Integer                                 :: num_all
    Integer,allocatable                     :: kk(:,:)
    Integer                                 :: b_MPI(2,0:MPI_np-1)

#if(run_MPI)

    IF((present(Hpsi) .AND. (.NOT. present(H_overlap))) .OR.                           &
       ((.NOT. present(H_overlap)) .AND. (.NOT. present(S_overlap)))) THEN
      STOP 'variable presented error in calculate_overlap_S1_MPI'  
    ENDIF

    num_all=(ndim2-ndim1+1)*(ndim1-1)+(ndim2-ndim1+1)*ndim2
    CALL allocate_array(kk,1,2,1,num_all)

    nn=0
    DO ii=1,ndim1-1
      DO jj=ndim1,ndim2
        nn=nn+1
        kk(1,nn)=ii
        kk(2,nn)=jj
      ENDDO
    ENDDO

    DO ii=ndim1,ndim2
      DO jj=1,ndim2
        nn=nn+1
        kk(1,nn)=ii
        kk(2,nn)=jj
      ENDDO
    ENDDO

    nb_per_MPI=num_all/MPI_np
    nb_rem_MPI=mod(num_all,MPI_np) 
    DO i_MPI=0,MPI_np-1
      b_MPI(1,i_MPI)=i_MPI*nb_per_MPI+1+MIN(i_MPI,nb_rem_MPI)
      b_MPI(2,i_MPI)=(i_MPI+1)*nb_per_MPI+MIN(i_MPI,nb_rem_MPI)                        &
                                         +merge(1,0,nb_rem_MPI>i_MPI)
    ENDDO

    !-calculate H_overlap---------------------------------------------------------------
    IF(present(Hpsi) .AND. present(H_overlap)) THEN

      IF(MPI_id/=0) CALL allocate_array(vec,b_MPI(1,MPI_id),b_MPI(2,MPI_id))
      DO ii=b_MPI(1,MPI_id),b_MPI(2,MPI_id)
        CALL Overlap_psi1_psi2(Overlap,psi(kk(2,ii)),Hpsi(kk(1,ii)),With_Grid=With_Grid)
        IF(MPI_id==0) THEN
          H_overlap(kk(2,ii),kk(1,ii))=real(Overlap,kind=Rkind)
        ELSE
          vec(ii)=real(Overlap,kind=Rkind)
        ENDIF
      ENDDO

      ! collect result
      IF(MPI_id/=0) THEN
        CALL MPI_Send(vec(b_MPI(1,MPI_id):b_MPI(2,MPI_id)),                            &
                          b_MPI(2,MPI_id)-b_MPI(1,MPI_id)+1,Real_MPI,root_MPI,MPI_id,  &
                          MPI_COMM_WORLD,MPI_err)
      ELSE
        DO i_MPI=1,MPI_np-1
          CALL allocate_array(vec,b_MPI(1,i_MPI),b_MPI(2,i_MPI))
          CALL MPI_Recv(vec(b_MPI(1,i_MPI):b_MPI(2,i_MPI)),                            &
                            b_MPI(2,i_MPI)-b_MPI(1,i_MPI)+1,Real_MPI,i_MPI,i_MPI,      &
                            MPI_COMM_WORLD,MPI_stat,MPI_err)
          DO ii=b_MPI(1,i_MPI),b_MPI(2,i_MPI)
            H_overlap(kk(2,ii),kk(1,ii))=vec(ii)
          ENDDO
        ENDDO ! for i_MPI=1,MPI_np-1
      ENDIF ! for MPI_id

      IF(ndim1>1) CALL MPI_Bcast_matrix(H_overlap,ndim1,ndim2,1,ndim1-1,root_MPI)
      CALL MPI_Bcast_matrix(H_overlap,1,ndim2,ndim1,ndim2,root_MPI)

    ENDIF ! for present(Hpsi) .AND. present(H_overlap)


    ! improve for S sym later
    !-overlap <psi|psi>-----------------------------------------------------------------
    IF(present(S_overlap)) THEN
      IF(MPI_id/=0) CALL allocate_array(vec,b_MPI(1,MPI_id),b_MPI(2,MPI_id))
      DO ii=b_MPI(1,MPI_id),b_MPI(2,MPI_id)
        CALL Overlap_psi1_psi2(Overlap,psi(kk(2,ii)),psi(kk(1,ii)),With_Grid=With_Grid)
        IF(MPI_id==0) THEN
          S_overlap(kk(2,ii),kk(1,ii))=real(Overlap,kind=Rkind)
        ELSE
          vec(ii)=real(Overlap,kind=Rkind)
        ENDIF
      ENDDO

      ! collect result
      IF(MPI_id/=0) THEN
        CALL MPI_Send(vec(b_MPI(1,MPI_id):b_MPI(2,MPI_id)),                            &
                          b_MPI(2,MPI_id)-b_MPI(1,MPI_id)+1,Real_MPI,root_MPI,MPI_id,  &
                          MPI_COMM_WORLD,MPI_err)
      ELSE
        DO i_MPI=1,MPI_np-1
          CALL allocate_array(vec,b_MPI(1,i_MPI),b_MPI(2,i_MPI))
          CALL MPI_Recv(vec(b_MPI(1,i_MPI):b_MPI(2,i_MPI)),                            &
                            b_MPI(2,i_MPI)-b_MPI(1,i_MPI)+1,Real_MPI,i_MPI,i_MPI,      &
                            MPI_COMM_WORLD,MPI_stat,MPI_err)
          DO ii=b_MPI(1,i_MPI),b_MPI(2,i_MPI)
            S_overlap(kk(2,ii),kk(1,ii))=vec(ii)
          ENDDO
        ENDDO ! for i_MPI=1,MPI_np-1
      ENDIF ! for MPI_id

      IF(ndim1>1) CALL MPI_Bcast_matrix(S_overlap,ndim1,ndim2,1,ndim1-1,root_MPI)
      CALL MPI_Bcast_matrix(S_overlap,1,ndim2,ndim1,ndim2,root_MPI)

    ENDIF ! for present(S_overlap)

    IF(allocated(vec)) deallocate(vec)

#endif
  ENDSUBROUTINE calculate_overlap_S1_MPI
!=======================================================================================


!=======================================================================================
! calculate overlap for <psi(n+1)|psi(i)> or/and <Hpsi(n+1)|psi(i)> i=1...n
!---------------------------------------------------------------------------------------
  SUBROUTINE calculate_overlap1D_MPI(psi,ndim,With_Grid,Hpsi,S_overlap1D,H_overlap1D)
    USE mod_system
    USE mod_psi_set_alloc
    USE mod_MPI_aux
    IMPLICIT NONE
    
    TYPE(param_psi),          intent(in)    :: psi(:) 
    TYPE(param_psi),optional, intent(in)    :: Hpsi(:) 
    Integer,                  intent(in)    :: ndim
    Logical,optional,         intent(in)    :: With_Grid
    Real(kind=Rkind),optional,intent(inout) :: H_overlap1D(:)
    Real(kind=Rkind),optional,intent(inout) :: S_overlap1D(:)
    
    Complex(kind=Rkind)                     :: Overlap
    Real(kind=Rkind)                        :: Overlap1D_temp(ndim)
    Integer                                 :: i,j
    Logical                                 :: With_Grid_loc

#if(run_MPI)

    With_Grid_loc=.FALSE.
    IF(present(With_Grid)) With_Grid_loc=With_Grid

    IF((present(Hpsi) .AND. (.NOT. present(H_overlap1D))) .OR.                         &
       ((.NOT. present(S_overlap1D)) .AND. (.NOT. present(H_overlap1D)))) THEN
      write(out_unitp,*) 'variable presented error in calculate_overlap_MPI'  
      STOP
    ENDIF

    !-overlap <psi|Hpsi>----------------------------------------------------------------
    IF(present(Hpsi) .AND. present(H_overlap1D)) THEN 
      DO i=1,ndim
        CALL Overlap_psipsi_MPI3(Overlap,psi(i),Hpsi(ndim),With_Grid=With_Grid)
        H_overlap1D(i)=real(Overlap,kind=Rkind) 
      END DO

      ! collect and broadcast result
!      CALL MPI_Reduce(H_overlap1D,Overlap1D_temp,ndim,MPI_Real8,MPI_SUM,root_MPI,      &
!                      MPI_COMM_WORLD,MPI_err)
!      IF(MPI_id==0) H_overlap1D=Overlap1D_temp
!      CALL MPI_BCAST(H_overlap1D,ndim,MPI_Real8,root_MPI,MPI_COMM_WORLD,MPI_err)
      CALL MPI_Reduce_sum_Bcast(H_overlap1D,ndim)
    ENDIF 

    !-overlap <psi|psi>-----------------------------------------------------------------  
    IF(present(S_overlap1D)) THEN 
      DO i=1,ndim
        CALL Overlap_psipsi_MPI3(Overlap,psi(i),psi(ndim),With_Grid=With_Grid)
        S_overlap1D(i)=real(Overlap,kind=Rkind)
      ENDDO 

      ! collect and broadcast result
!      CALL MPI_Reduce(S_overlap1D(1:ndim),Overlap1D_temp,ndim,Real_MPI,MPI_SUM,        &
!                      root_MPI,MPI_COMM_WORLD,MPI_err)
!      IF(MPI_id==0) S_overlap1D(1:ndim)=Overlap1D_temp(1:ndim)
!      CALL MPI_BCAST(S_overlap1D,ndim,MPI_Real8,root_MPI,MPI_COMM_WORLD,MPI_err)
      CALL MPI_Reduce_sum_Bcast(S_overlap1D,ndim)
    ENDIF

#endif
  END SUBROUTINE calculate_overlap1D_MPI
!=======================================================================================


!=======================================================================================
!> subroutine for the calculation of Overlap_psi1_psi2 with MPI 
!> overlap is cut into submatrix for MPI parallel   
!> note, this requre new  Overlap_psi1_psi2 routine
!
!> the submatrix of psi and Hpsi are keeped also for the calculation 
!> of Residual g in MakeResidual_Davidson_MPI2
!---------------------------------------------------------------------------------------
  SUBROUTINE Overlap_psi1_psi2_MPI5(H_overlap,S_overlap,psi,Hpsi,ndim0,ndim,With_Grid)
    USE mod_system
    USE mod_psi_set_alloc
    USE mod_MPI_aux
    IMPLICIT NONE

    TYPE(param_psi), intent(inout)          :: psi(:)  !< on non-root threads allocated
    TYPE(param_psi), intent(inout)          :: Hpsi(:) !< on non-root threads allocated
    Real(kind=Rkind),intent(inout)          :: H_overlap(:,:)
    Real(kind=Rkind),intent(inout)          :: S_overlap(:,:)
    Logical,optional,intent(in)             :: With_Grid
    Integer,         intent(in)             :: ndim0
    Integer,         intent(in)             :: ndim

    complex(kind=Rkind)                     :: Overlap
    logical                                 :: With_Grid_loc
    Integer                                 :: i
    Integer                                 :: j
    Logical                                 :: root_jobs

#if(run_MPI)

    With_Grid_loc=.FALSE.
    IF (present(With_Grid)) With_Grid_loc = With_Grid

    ! only on master: psi%CvecG,psi%RvecG,psi%CvecB,psi%RvecB
    CALL distribute_psi_pack_MPI( psi,ndim0+1,ndim,With_Grid)
    CALL distribute_psi_pack_MPI(Hpsi,ndim0+1,ndim,With_Grid)

    CALL calculate_overlap_MPI(psi,ndim0+1,ndim,With_Grid,Hpsi=Hpsi,                   &
                               S_overlap=S_overlap,H_overlap=H_overlap)

    ! note <psi|Hpsi> should be completely renew 
  !  CALL calculate_overlap_MPI(psi,ndim0+1,ndim,With_Grid,S_overlap=S_overlap)
  !  CALL calculate_overlap_MPI(psi,ndim0+1,ndim,With_Grid,Hpsi=Hpsi,H_overlap=H_overlap)

#endif
  END SUBROUTINE Overlap_psi1_psi2_MPI5
!=======================================================================================


!=======================================================================================
!> subroutine for the calculation of Overlap_psi1_psi2 with MPI 
!> overlap is cut into submatrix for MPI parallel   
!> note, this requre new  Overlap_psi1_psi2 routine
!
!> the submatrix of psi and Hpsi are keeped also for the calculation 
!> of Residual g in MakeResidual_Davidson_MPI2
!---------------------------------------------------------------------------------------
  SUBROUTINE Overlap_psi1_psi2_MPI3(H_overlap,S_overlap,psi,Hpsi,ndim0,ndim,With_Grid)
    USE mod_system
    USE mod_psi_set_alloc
    USE mod_MPI_aux
    IMPLICIT NONE
    
    TYPE(param_psi), intent(inout)          :: psi(:)  !< on non-root threads allocated
    TYPE(param_psi), intent(inout)          :: Hpsi(:) !< on non-root threads allocated
    Real(kind=Rkind),intent(inout)          :: H_overlap(:,:)
    Real(kind=Rkind),intent(inout)          :: S_overlap(:,:)
    Logical,optional,intent(in)             :: With_Grid
    Integer,         intent(in)             :: ndim0
    Integer,         intent(in)             :: ndim

    complex(kind=Rkind)                     :: Overlap
    logical                                 :: With_Grid_loc
    Integer                                 :: i
    Integer                                 :: j
    Logical                                 :: root_jobs

#if(run_MPI)

    With_Grid_loc=.FALSE.
    IF (present(With_Grid)) With_Grid_loc = With_Grid

    ! only on master: psi%CvecG,psi%RvecG,psi%CvecB,psi%RvecB
    IF(MPI_id==0) THEN
      DO i=ndim0+1,ndim
        
        IF(.NOT. With_Grid_loc) THEN
          nb_per_MPI=Hpsi(i)%nb_tot/MPI_np
          nb_rem_MPI=mod(Hpsi(i)%nb_tot,MPI_np) !remainder jobs
        ELSE
          nb_per_MPI=Hpsi(i)%nb_qaie/MPI_np
          nb_rem_MPI=mod(Hpsi(i)%nb_qaie,MPI_np) !remainder jobs
        ENDIF
        
        ! Send array
        DO i_MPI=1,MPI_np-1
          bound1_MPI=i_MPI*nb_per_MPI+1+MIN(i_MPI,nb_rem_MPI)
          bound2_MPI=(i_MPI+1)*nb_per_MPI+MIN(i_MPI,nb_rem_MPI)+merge(1,0,nb_rem_MPI>i_MPI)
          IF (.NOT. With_Grid_loc) THEN
            IF (psi(i)%cplx) THEN
              CALL MPI_Send(Hpsi(i)%CvecB(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,&
                            MPI_complex8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
              CALL MPI_Send(psi(i)%CvecB(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1, &
                            MPI_Complex8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
            ELSE
              CALL MPI_Send(Hpsi(i)%RvecB(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,&
                            MPI_Real8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
              CALL MPI_Send(psi(i)%RvecB(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1, &
                            MPI_Real8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
            
            ENDIF
          ELSE
            IF (psi(i)%cplx) THEN
              CALL MPI_Send(Hpsi(i)%CvecG(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,&
                            MPI_complex8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
              CALL MPI_Send(psi(i)%CvecG(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1, &
                            MPI_Complex8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)        
            ELSE
              CALL MPI_Send(Hpsi(i)%RvecG(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,&
                            MPI_Real8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
              CALL MPI_Send(psi(i)%RvecG(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1, &
                            MPI_Real8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
            
            ENDIF
          ENDIF ! for .NOT. With_Grid_loc 
        ENDDO ! for i_MPI=1,MPI_np-1
      ENDDO ! for i=ndim0+1,ndim

      !> calcuation on master
      DO i=1,ndim0
        DO j=ndim0+1,ndim
          CALL Overlap_psipsi_MPI3(Overlap,psi(j),Hpsi(i),With_Grid=With_Grid)
          H_overlap(j,i)=real(Overlap,kind=Rkind) 
          CALL Overlap_psipsi_MPI3(Overlap,psi(j), psi(i),With_Grid=With_Grid)
          S_overlap(j,i)=real(Overlap,kind=Rkind)
        ENDDO ! for j=ndim0+1,ndim
      END DO ! for i=ndim0+1,ndim
      
      DO i=ndim0+1,ndim
        DO j=1,ndim
          CALL Overlap_psipsi_MPI3(Overlap,psi(j),Hpsi(i),With_Grid=With_Grid)
          H_overlap(j,i)=real(Overlap,kind=Rkind) 
          CALL Overlap_psipsi_MPI3(Overlap,psi(j), psi(i),With_Grid=With_Grid)
          S_overlap(j,i)=real(Overlap,kind=Rkind)
        ENDDO ! for j=ndim0+1,ndim
      END DO ! for i=ndim0+1,ndim

    ENDIF ! for MPI_id==0

    ! MPI/=0------------------------------------------------------------------------------
    IF(MPI_id/=0) THEN
      DO i=ndim0+1,ndim
        IF(.NOT. With_Grid_loc) THEN
          nb_per_MPI=Hpsi(i)%nb_tot/MPI_np
          nb_rem_MPI=mod(Hpsi(i)%nb_tot,MPI_np) !remainder jobs
        ELSE
          nb_per_MPI=Hpsi(i)%nb_qaie/MPI_np
          nb_rem_MPI=mod(Hpsi(i)%nb_qaie,MPI_np) !remainder jobs
        ENDIF
        bound1_MPI=MPI_id*nb_per_MPI+1+MIN(MPI_id,nb_rem_MPI)
        bound2_MPI=(MPI_id+1)*nb_per_MPI+MIN(MPI_id,nb_rem_MPI)+merge(1,0,nb_rem_MPI>MPI_id)
        
        ! receive array
        IF (.NOT. With_Grid_loc) THEN
          IF (psi(i)%cplx) THEN
            IF(.NOT. allocated(Hpsi(i)%CvecB))                                           &
                allocate(Hpsi(i)%CvecB(bound1_MPI:bound2_MPI))
            CALL MPI_Recv(Hpsi(i)%CvecB(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,  &
                          MPI_complex8,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
            IF(.NOT. allocated(psi(i)%CvecB))                                            &
                 allocate(psi(i)%CvecB(bound1_MPI:bound2_MPI))
            CALL MPI_Recv( psi(i)%CvecB(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,  &
                          MPI_Complex8,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
          ELSE
            IF(.NOT. allocated(Hpsi(i)%RvecB))                                           &
                 allocate(Hpsi(i)%RvecB(bound1_MPI:bound2_MPI))
            CALL MPI_Recv(Hpsi(i)%RvecB(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,  &
                          MPI_Real8,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
          
            IF(.NOT. allocated(psi(i)%RvecB))                                            &
                 allocate(psi(i)%RvecB(bound1_MPI:bound2_MPI))
            CALL MPI_Recv( psi(i)%RvecB(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,  &
                          MPI_Real8,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
            
          ENDIF
        ELSE
          IF (psi(i)%cplx) THEN
            IF(.NOT. allocated(Hpsi(i)%CvecG))                                           &
                 allocate(Hpsi(i)%CvecG(bound1_MPI:bound2_MPI))
            CALL MPI_Recv(Hpsi(i)%CvecG(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,  &
                          MPI_complex8,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
            IF(.NOT. allocated(psi(i)%CvecG))                                            &
                 allocate(psi(i)%CvecG(bound1_MPI:bound2_MPI))
            CALL MPI_Recv( psi(i)%CvecG(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,  &
                          MPI_Complex8,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
          ELSE
            IF(.NOT. allocated(Hpsi(i)%RvecB))                                           &
                 allocate(Hpsi(i)%RvecG(bound1_MPI:bound2_MPI))
            CALL MPI_Recv(Hpsi(i)%RvecB(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,  &
                          MPI_Real8,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
            IF(.NOT. allocated(psi(i)%RvecB))                                            &
                 allocate(psi(i)%RvecG(bound1_MPI:bound2_MPI))
            CALL MPI_Recv( psi(i)%RvecB(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,  &
                          MPI_Real8,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
          ENDIF
        ENDIF
      ENDDO ! for i=ndim0+1,ndim

      !> calcuation on each thread
      DO i=1,ndim0
        DO j=ndim0+1,ndim
          CALL Overlap_psipsi_MPI3(Overlap,psi(j),Hpsi(i),With_Grid=With_Grid)
          H_overlap(j,i)=real(Overlap,kind=Rkind) 
          CALL Overlap_psipsi_MPI3(Overlap,psi(j), psi(i),With_Grid=With_Grid)
          S_overlap(j,i)=real(Overlap,kind=Rkind)
        ENDDO ! for j=ndim0+1,ndim
      END DO ! for i=ndim0+1,ndim
      
      DO i=ndim0+1,ndim
        DO j=1,ndim
          CALL Overlap_psipsi_MPI3(Overlap,psi(j),Hpsi(i),With_Grid=With_Grid)
          H_overlap(j,i)=real(Overlap,kind=Rkind) 
          CALL Overlap_psipsi_MPI3(Overlap,psi(j), psi(i),With_Grid=With_Grid)
          S_overlap(j,i)=real(Overlap,kind=Rkind)
        ENDDO ! for j=ndim0+1,ndim
      END DO ! for i=ndim0+1,ndim

    ENDIF ! for MPI_id/=0  

    CALL MPI_Reduce_sum_matrix(H_overlap,ndim0+1,ndim,1,ndim0,root_MPI)
    CALL MPI_Bcast_matrix     (H_overlap,ndim0+1,ndim,1,ndim0,root_MPI)
    CALL MPI_Reduce_sum_matrix(H_overlap,1,ndim,ndim0+1,ndim,root_MPI)
    CALL MPI_Bcast_matrix     (H_overlap,1,ndim,ndim0+1,ndim,root_MPI)

    CALL MPI_Reduce_sum_matrix(S_overlap,ndim0+1,ndim,1,ndim0,root_MPI)
    CALL MPI_Bcast_matrix     (S_overlap,ndim0+1,ndim,1,ndim0,root_MPI)
    CALL MPI_Reduce_sum_matrix(S_overlap,1,ndim,ndim0+1,ndim,root_MPI)
    CALL MPI_Bcast_matrix     (S_overlap,1,ndim,ndim0+1,ndim,root_MPI)

#endif
  END SUBROUTINE Overlap_psi1_psi2_MPI3
!=======================================================================================

!=======================================================================================
! subroutine for the calculation of Overlap_psi1_psi2 with MPI 
! ONLY the overlap of i=j are calculated on the other threads to reduce memory 
! psi and Hpsi are keeped also for the calculation 
! of Residual g in MakeResidual_Davidson_MPI2
!---------------------------------------------------------------------------------------
  SUBROUTINE Overlap_psi1_psi2_MPI2(H_overlap,S_overlap,psi,Hpsi,ndim,With_Grid)
    USE mod_system
    USE mod_psi_set_alloc
    USE mod_MPI_aux
    IMPLICIT NONE
    
    TYPE(param_psi), intent(inout)          :: psi(:)  !< on non-root threads allocated
    TYPE(param_psi), intent(inout)          :: Hpsi(:) !< on non-root threads allocated
    Real(kind=Rkind),intent(inout)          :: H_overlap(:,:)
    Real(kind=Rkind),intent(inout)          :: S_overlap(:,:)
    Logical,optional,intent(in)             :: With_Grid
    Integer,         intent(in)             :: ndim

    complex(kind=Rkind)                     :: Overlap
    logical                                 :: With_Grid_loc
    Integer                                 :: i
    Integer                                 :: j
    Logical                                 :: send_once(MPI_np-1)
    Logical                                 :: root_jobs

#if(run_MPI)

    With_Grid_loc=.FALSE.
    IF (present(With_Grid)) With_Grid_loc = With_Grid

    ! only on master: psi%CvecG,psi%RvecG,psi%CvecB,psi%RvecB
    nb_per_MPI=(ndim)/MPI_np
    nb_rem_MPI=mod(ndim,MPI_np) !remainder jobs
    IF(MPI_id==0) THEN
      !> send Hpsi and psi
      DO i_MPI=1,MPI_np-1
        bound1_MPI=i_MPI*nb_per_MPI+1+MIN(i_MPI,nb_rem_MPI)
        bound2_MPI=(i_MPI+1)*nb_per_MPI+MIN(i_MPI,nb_rem_MPI)+merge(1,0,nb_rem_MPI>i_MPI)
        DO i=1,ndim
          IF(i>=bound1_MPI .AND. i<=bound2_MPI) THEN
            IF (.NOT. With_Grid_loc) THEN
              IF (psi(i)%cplx) THEN
                CALL MPI_Send(Hpsi(i)%CvecB,Hpsi(i)%nb_tot,MPI_complex8,i_MPI,         &
                              i_MPI,MPI_COMM_WORLD,MPI_err)
                CALL MPI_Send(psi(i)%CvecB,psi(i)%nb_tot,MPI_Complex8,i_MPI,           &
                              i_MPI,MPI_COMM_WORLD,MPI_err)

              ELSE
                CALL MPI_Send(Hpsi(i)%RvecB,Hpsi(i)%nb_tot,MPI_Real8,i_MPI,            &
                              i_MPI,MPI_COMM_WORLD,MPI_err)
                CALL MPI_Send(psi(i)%RvecB,psi(i)%nb_tot,MPI_Real8,i_MPI,              &
                              i_MPI,MPI_COMM_WORLD,MPI_err)
              ENDIF
            ELSE
              IF (psi(i)%cplx) THEN
                CALL MPI_Send(Hpsi(i)%CvecG,Hpsi(i)%nb_qaie,MPI_complex8,i_MPI,        &
                              i_MPI,MPI_COMM_WORLD,MPI_err)
                CALL MPI_Send(psi(i)%CvecG,psi(i)%nb_qaie,MPI_Complex8,i_MPI,          &
                              i_MPI,MPI_COMM_WORLD,MPI_err)
              ELSE
                CALL MPI_Send(Hpsi(i)%RvecG,Hpsi(i)%nb_qaie,MPI_Real8,i_MPI,           &
                              i_MPI,MPI_COMM_WORLD,MPI_err)
                CALL MPI_Send(psi(i)%RvecG,psi(i)%nb_qaie,MPI_Real8,i_MPI,             &
                              i_MPI,MPI_COMM_WORLD,MPI_err)
              ENDIF
            ENDIF ! for .NOT. With_Grid_loc 
          ENDIF ! for i>=bound1_MPI .AND. i<=bound2_MPI
        ENDDO ! for i=1,ndim
      ENDDO ! for i_MPI=1,MPI_np-1

      !> calcuation on master
      DO i=1,ndim
        DO j=1,ndim
          root_jobs=.TRUE.
          DO i_MPI=1,MPI_np-1
            bound1_MPI=i_MPI*nb_per_MPI+1+MIN(i_MPI,nb_rem_MPI)
            bound2_MPI=(i_MPI+1)*nb_per_MPI+MIN(i_MPI,nb_rem_MPI)                      &
                                           +merge(1,0,nb_rem_MPI>i_MPI)
            IF((i>=bound1_MPI .AND. i<=bound2_MPI) .AND.                               &
               (j>=bound1_MPI .AND. j<=bound2_MPI)) THEN
               root_jobs=.FALSE.
            ENDIF
          ENDDO ! i_MPI=1,MPI_np-1
          
          IF(root_jobs) THEN
            CALL Overlap_psipsi_MPI(Overlap,psi(j),Hpsi(i),With_Grid=With_Grid)
            H_overlap(j,i)=real(Overlap,kind=Rkind) 
            CALL Overlap_psipsi_MPI(Overlap,psi(j), psi(i),With_Grid=With_Grid)
            S_overlap(j,i)=real(Overlap,kind=Rkind)
          ENDIF
        ENDDO ! for j=1,ndim
      END DO ! for i=1,ndim

      DO i_MPI=1,MPI_np-1
        bound1_MPI=i_MPI*nb_per_MPI+1+MIN(i_MPI,nb_rem_MPI)
        bound2_MPI=(i_MPI+1)*nb_per_MPI+MIN(i_MPI,nb_rem_MPI)+merge(1,0,nb_rem_MPI>i_MPI)
        CALL MPI_Recv_matrix(H_overlap,bound1_MPI,bound2_MPI,bound1_MPI,bound2_MPI,    &
                             i_MPI,i_MPI)
        CALL MPI_Recv_matrix(S_overlap,bound1_MPI,bound2_MPI,bound1_MPI,bound2_MPI,    &
                             i_MPI,i_MPI)
      ENDDO 

    ENDIF ! for MPI_id==0

    !-------------------------------------------------------------------------------------
    bound1_MPI=MPI_id*nb_per_MPI+1+MIN(MPI_id,nb_rem_MPI)
    bound2_MPI=(MPI_id+1)*nb_per_MPI+MIN(MPI_id,nb_rem_MPI)+merge(1,0,nb_rem_MPI>MPI_id)
    IF(MPI_id/=0) THEN
      DO i=1,ndim
        IF(i>=bound1_MPI .AND. i<=bound2_MPI) THEN
          IF (.NOT. With_Grid_loc) THEN
            IF (psi(i)%cplx) THEN
              IF(.NOT. allocated(Hpsi(i)%CvecB)) CALL alloc_NParray(Hpsi(i)%CvecB,     &
                                            [Hpsi(i)%nb_tot],'Hpsi%CvecB','alloc_psi')
              CALL MPI_Recv(Hpsi(i)%CvecB,Hpsi(i)%nb_tot,MPI_complex8,root_MPI,MPI_id, &
                            MPI_COMM_WORLD,MPI_stat,MPI_err)
              IF(.NOT. allocated(psi(i)%CvecB))  CALL alloc_NParray(psi(i)%CvecB,      &
                                              [psi(i)%nb_tot],'psi%CvecB','alloc_psi')
              CALL MPI_Recv( psi(i)%CvecB, psi(i)%nb_tot,MPI_Complex8,root_MPI,MPI_id, &
                            MPI_COMM_WORLD,MPI_stat,MPI_err)
            ELSE
              IF(.NOT. allocated(Hpsi(i)%RvecB)) CALL alloc_NParray(Hpsi(i)%RvecB,     &
                                            [Hpsi(i)%nb_tot],'Hpsi%RvecB','alloc_psi')
              CALL MPI_Recv(Hpsi(i)%RvecB,Hpsi(i)%nb_tot,MPI_Real8,root_MPI,MPI_id,    &
                            MPI_COMM_WORLD,MPI_stat,MPI_err)
            
              IF(.NOT. allocated(psi(i)%RvecB))  CALL alloc_NParray(psi(i)%RvecB,      &
                                              [psi(i)%nb_tot],'psi%RvecB','alloc_psi')
              CALL MPI_Recv( psi(i)%RvecB, psi(i)%nb_tot,MPI_Real8,root_MPI,MPI_id,    &
                            MPI_COMM_WORLD,MPI_stat,MPI_err)
            ENDIF
          ELSE
            IF (psi(i)%cplx) THEN
              IF(.NOT. allocated(Hpsi(i)%CvecG)) CALL alloc_NParray(Hpsi(i)%CvecG,     &
                                           [Hpsi(i)%nb_qaie],'Hpsi%CvecG','alloc_psi')
              CALL MPI_Recv(Hpsi(i)%CvecG,Hpsi(i)%nb_qaie,MPI_complex8,root_MPI,MPI_id,&
                            MPI_COMM_WORLD,MPI_stat,MPI_err)
              IF(.NOT. allocated(psi(i)%CvecG))  CALL alloc_NParray(psi(i)%CvecG,      &
                                             [psi(i)%nb_qaie],'psi%CvecG','alloc_psi')
              CALL MPI_Recv( psi(i)%CvecG, psi(i)%nb_qaie,MPI_Complex8,root_MPI,MPI_id,&
                            MPI_COMM_WORLD,MPI_stat,MPI_err)
            ELSE
              IF(.NOT. allocated(Hpsi(i)%RvecB)) CALL alloc_NParray(Hpsi(i)%RvecG,     &
                                            [Hpsi(i)%nb_tot],'Hpsi%RvecG','alloc_psi')
              CALL MPI_Recv(Hpsi(i)%RvecB,Hpsi(i)%nb_tot,MPI_Real8,root_MPI,MPI_id,    &
                            MPI_COMM_WORLD,MPI_stat,MPI_err)
              IF(.NOT. allocated(psi(i)%RvecB))  CALL alloc_NParray(psi(i)%RvecG,      &
                                              [psi(i)%nb_tot],'psi%RvecG','alloc_psi')
              CALL MPI_Recv( psi(i)%RvecB, psi(i)%nb_tot,MPI_Real8,root_MPI,MPI_id,    &
                            MPI_COMM_WORLD,MPI_stat,MPI_err)
            ENDIF
          ENDIF ! for .NOT. With_Grid_loc
        ENDIF ! for i>=bound1_MPI .AND. i<=bound2_MPI
      ENDDO ! for i=1,ndim

      ! main calculation
      DO i=1,ndim
        DO j=1,ndim
          IF((i>=bound1_MPI .AND. i<=bound2_MPI) .AND.                                 &
             (j>=bound1_MPI .AND. j<=bound2_MPI)) THEN
            IF (.NOT. With_Grid_loc) THEN
              IF (psi(j)%symab > -1 .AND. Hpsi(i)%symab > -1                           &
                                    .AND. Hpsi(i)%symab /= Hpsi(i)%symab) THEN
              ELSE
                CALL Overlap_psipsi_MPI(Overlap,psi(j),Hpsi(i),With_Grid=With_Grid)
                H_overlap(j,i)=real(Overlap,kind=Rkind)
                CALL Overlap_psipsi_MPI(Overlap,psi(j), psi(i),With_Grid=With_Grid)
                S_overlap(j,i)=real(Overlap,kind=Rkind)
              ENDIF
            ENDIF
          ENDIF    
        ENDDO
      ENDDO

      CALL MPI_Send_matrix(H_overlap,bound1_MPI,bound2_MPI,bound1_MPI,bound2_MPI,      &
                           root_MPI,MPI_id)
      CALL MPI_Send_matrix(S_overlap,bound1_MPI,bound2_MPI,bound1_MPI,bound2_MPI,      &
                           root_MPI,MPI_id)
    ENDIF ! for MPI_id/=0  

  !  DO i_MPI=1,MPI_np-1
  !    bound1_MPI=i_MPI*nb_per_MPI+1+MIN(i_MPI,nb_rem_MPI)
  !    bound2_MPI=(i_MPI+1)*nb_per_MPI+MIN(i_MPI,nb_rem_MPI)+merge(1,0,nb_rem_MPI>i_MPI)
  !    CALL MPI_Bcast_matrix(H_overlap,bound1_MPI,bound2_MPI,bound1_MPI,bound2_MPI,i_MPI)
  !    CALL MPI_Bcast_matrix(S_overlap,bound1_MPI,bound2_MPI,bound1_MPI,bound2_MPI,i_MPI)
  !  ENDDO

    ! a bit waste of comm. time
    CALL MPI_Bcast_matrix(H_overlap,1,ndim,1,ndim,root_MPI)
    CALL MPI_Bcast_matrix(S_overlap,1,ndim,1,ndim,root_MPI)

#endif
  END SUBROUTINE Overlap_psi1_psi2_MPI2
!=======================================================================================

!=======================================================================================
! subroutine for the calculation of Overlap_psi1_psi2 with MPI 
! in loop: 
! DO i=i_l,i_u
! DO j=j_l,j_u
!---------------------------------------------------------------------------------------
  SUBROUTINE Overlap_psi1_psi2_MPI(H_overlap,S_overlap,psi,Hpsi,ndim,With_Grid)
    USE mod_system
    USE mod_psi_set_alloc
    USE mod_MPI_aux
    IMPLICIT NONE

    TYPE(param_psi), intent(inout)          :: psi(:)  !< on non-root threads allocated
    TYPE(param_psi), intent(inout)          :: Hpsi(:) !< on non-root threads allocated
    Real(kind=Rkind),intent(inout)          :: H_overlap(:,:)
    Real(kind=Rkind),intent(inout)          :: S_overlap(:,:)
    Logical,optional,intent(in)             :: With_Grid
    Integer,intent(in)                      :: ndim

    complex(kind=Rkind)                     :: Overlap
    logical                                 :: With_Grid_loc
    Integer                                 :: i
    Integer                                 :: j
    Logical                                 :: send_once(MPI_np-1)

#if(run_MPI)

    With_Grid_loc=.FALSE.
    IF (present(With_Grid)) With_Grid_loc = With_Grid

    ! only on master: psi%CvecG,psi%RvecG,psi%CvecB,psi%RvecB
    nb_per_MPI=(ndim)/MPI_np
    nb_rem_MPI=mod(ndim,MPI_np) !remainder jobs
    
    IF(MPI_id==0) THEN
      !> send Hpsi and psi
      DO i_MPI=1,MPI_np-1
        bound1_MPI=i_MPI*nb_per_MPI+1+MIN(i_MPI,nb_rem_MPI)
        bound2_MPI=(i_MPI+1)*nb_per_MPI+MIN(i_MPI,nb_rem_MPI)                          &
                                       +merge(1,0,nb_rem_MPI>i_MPI)
        DO i=1,ndim
          IF(.NOT. With_Grid_loc) THEN
            IF(psi(i)%cplx) THEN
              IF(i>=bound1_MPI .AND. i<=bound2_MPI) CALL MPI_Send(Hpsi(i)%CvecB,       &
                           Hpsi(i)%nb_tot,Cplx_MPI,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
              CALL MPI_Send(psi(i)%CvecB,psi(i)%nb_tot,MPI_Complex8,i_MPI,             &
                            i_MPI,MPI_COMM_WORLD,MPI_err)
            ELSE
              IF(i>=bound1_MPI .AND. i<=bound2_MPI) CALL MPI_Send(Hpsi(i)%RvecB,       &
                             Hpsi(i)%nb_tot,Real_MPI,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
              CALL MPI_Send(psi(i)%RvecB,psi(i)%nb_tot,MPI_Real8,i_MPI,                &
                            i_MPI,MPI_COMM_WORLD,MPI_err)
            ENDIF
          ELSE
            IF(psi(i)%cplx) THEN
              IF(i>=bound1_MPI .AND. i<=bound2_MPI) CALL MPI_Send(Hpsi(i)%CvecG,       &
                          Hpsi(i)%nb_qaie,Cplx_MPI,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
              CALL MPI_Send(psi(i)%CvecG,psi(i)%nb_qaie,MPI_Complex8,i_MPI,            &
                            i_MPI,MPI_COMM_WORLD,MPI_err)
            ELSE
              IF(i>=bound1_MPI .AND. i<=bound2_MPI) CALL MPI_Send(Hpsi(i)%RvecG,       &
                            Hpsi(i)%nb_qaie,Real_MPI,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
              CALL MPI_Send(psi(i)%RvecG,psi(i)%nb_qaie,MPI_Real8,i_MPI,               &
                            i_MPI,MPI_COMM_WORLD,MPI_err)
            ENDIF
          ENDIF ! for .NOT. With_Grid_loc  
        ENDDO ! for i=1,ndim
      ENDDO ! for i_MPI=1,MPI_np-1
      
      !> calcuation on master
      i_MPI=0
      DO i=1,ndim
        bound1_MPI=i_MPI*nb_per_MPI+1+MIN(i_MPI,nb_rem_MPI)
        bound2_MPI=(i_MPI+1)*nb_per_MPI+MIN(i_MPI,nb_rem_MPI)+merge(1,0,nb_rem_MPI>i_MPI)
        IF(i>=bound1_MPI .AND. i<=bound2_MPI) THEN
          DO j=1,ndim
            ! main claculation on master
            CALL Overlap_psipsi_MPI(Overlap,psi(j),Hpsi(i),With_Grid=With_Grid)
            H_overlap(j,i)=real(Overlap,kind=Rkind) 
            CALL Overlap_psipsi_MPI(Overlap,psi(j), psi(i),With_Grid=With_Grid)
            S_overlap(j,i)=real(Overlap,kind=Rkind)
          ENDDO ! for j=1,ndim
        ENDIF ! for i>=bound1_MPI .AND. i<=bound2_MPI
      END DO ! for i=1,ndim
      
    ENDIF ! for MPI_id==0
      
    !-------------------------------------------------------------------------------------
    bound1_MPI=MPI_id*nb_per_MPI+1+MIN(MPI_id,nb_rem_MPI)
    bound2_MPI=(MPI_id+1)*nb_per_MPI+MIN(MPI_id,nb_rem_MPI)+merge(1,0,nb_rem_MPI>MPI_id)
    IF(MPI_id/=0) THEN
      DO i=1,ndim
        IF (.NOT. With_Grid_loc) THEN
          IF (psi(i)%cplx) THEN
            IF(i>=bound1_MPI .AND. i<=bound2_MPI) THEN
              IF(.NOT. allocated(Hpsi(i)%CvecB)) CALL alloc_NParray(Hpsi(i)%CvecB,     &
                                            [Hpsi(i)%nb_tot],'Hpsi%CvecB','alloc_psi')
              CALL MPI_Recv(Hpsi(i)%CvecB,Hpsi(i)%nb_tot,MPI_complex8,root_MPI,MPI_id, &
                            MPI_COMM_WORLD,MPI_stat,MPI_err)
            ENDIF          
            IF(.NOT. allocated(psi(i)%CvecB)) CALL alloc_NParray(psi(i)%CvecB,         &
                                              [psi(i)%nb_tot],'psi%CvecB','alloc_psi')
            CALL MPI_Recv( psi(i)%CvecB, psi(i)%nb_tot,MPI_Complex8,root_MPI,MPI_id,   &
                          MPI_COMM_WORLD,MPI_stat,MPI_err)
          ELSE
            IF(i>=bound1_MPI .AND. i<=bound2_MPI) THEN
              IF(.NOT. allocated(Hpsi(i)%RvecB)) CALL alloc_NParray(Hpsi(i)%RvecB,     &
                                            [Hpsi(i)%nb_tot],'Hpsi%RvecB','alloc_psi')
              CALL MPI_Recv(Hpsi(i)%RvecB,Hpsi(i)%nb_tot,MPI_Real8,root_MPI,MPI_id,    &
                            MPI_COMM_WORLD,MPI_stat,MPI_err)
            ENDIF
            IF(.NOT. allocated(psi(i)%RvecB))  CALL alloc_NParray(psi(i)%RvecB,        &
                                              [psi(i)%nb_tot],'psi%RvecB','alloc_psi')
            CALL MPI_Recv( psi(i)%RvecB, psi(i)%nb_tot,MPI_Real8,root_MPI,MPI_id,      &
                          MPI_COMM_WORLD,MPI_stat,MPI_err)
          ENDIF
        ELSE
          IF (psi(i)%cplx) THEN
            IF(i>=bound1_MPI .AND. i<=bound2_MPI) THEN
              IF(.NOT. allocated(Hpsi(i)%CvecG)) CALL alloc_NParray(Hpsi(i)%CvecG,     &
                                           [Hpsi(i)%nb_qaie],'Hpsi%CvecG','alloc_psi')
              CALL MPI_Recv(Hpsi(i)%CvecG,Hpsi(i)%nb_qaie,MPI_complex8,root_MPI,MPI_id,&
                            MPI_COMM_WORLD,MPI_stat,MPI_err)
            ENDIF
            IF(.NOT. allocated(psi(i)%CvecG))  CALL alloc_NParray(psi(i)%CvecG,        &
                                             [psi(i)%nb_qaie],'psi%CvecG','alloc_psi')
            CALL MPI_Recv( psi(i)%CvecG, psi(i)%nb_qaie,MPI_Complex8,root_MPI,MPI_id,  &
                          MPI_COMM_WORLD,MPI_stat,MPI_err)
          ELSE
            IF(i>=bound1_MPI .AND. i<=bound2_MPI) THEN
              IF(.NOT. allocated(Hpsi(i)%RvecB)) CALL alloc_NParray(Hpsi(i)%RvecB,     &
                                            [Hpsi(i)%nb_tot],'Hpsi%RvecB','alloc_psi')
              CALL MPI_Recv(Hpsi(i)%RvecB,Hpsi(i)%nb_tot,MPI_Real8,root_MPI,MPI_id,    &
                            MPI_COMM_WORLD,MPI_stat,MPI_err)
            ENDIF
            IF(.NOT. allocated(psi(i)%RvecB))  CALL alloc_NParray(psi(i)%RvecB,        &
                                              [psi(i)%nb_tot],'psi%RvecB','alloc_psi')
            CALL MPI_Recv( psi(i)%RvecB, psi(i)%nb_tot,MPI_Real8,root_MPI,MPI_id,      &
                          MPI_COMM_WORLD,MPI_stat,MPI_err)
          ENDIF
        ENDIF ! for .NOT. With_Grid_loc
      ENDDO ! for i=1,ndim
      
      ! main calculation
      DO i=1,ndim
        IF(i>=bound1_MPI .AND. i<=bound2_MPI) THEN
          DO j=1,ndim
            IF (.NOT. With_Grid_loc) THEN
              IF (psi(j)%symab > -1 .AND. Hpsi(i)%symab > -1                           &
                                    .AND. Hpsi(i)%symab /= Hpsi(i)%symab) THEN
              ELSE
                CALL Overlap_psipsi_MPI(Overlap,psi(j),Hpsi(i),With_Grid=With_Grid)
                H_overlap(j,i)=real(Overlap,kind=Rkind)
                CALL Overlap_psipsi_MPI(Overlap,psi(j), psi(i),With_Grid=With_Grid)
                S_overlap(j,i)=real(Overlap,kind=Rkind)
              ENDIF
            ENDIF
          ENDDO
        ENDIF
      ENDDO
    ENDIF ! for MPI_id/=0  

    DO i_MPI=0,MPI_np-1
      bound1_MPI=i_MPI*nb_per_MPI+1+MIN(i_MPI,nb_rem_MPI)
      bound2_MPI=(i_MPI+1)*nb_per_MPI+MIN(i_MPI,nb_rem_MPI)+merge(1,0,nb_rem_MPI>i_MPI)
      CALL MPI_Bcast_matrix(H_overlap,1,ndim,bound1_MPI,bound2_MPI,i_MPI)
      CALL MPI_Bcast_matrix(S_overlap,1,ndim,bound1_MPI,bound2_MPI,i_MPI)
    ENDDO

#endif
  END SUBROUTINE Overlap_psi1_psi2_MPI
!=======================================================================================

!=======================================================================================
!> subroutine calculating overlap of psi1 and psi2 with MPI
!>  be careful with the way distributing array, 
!>  which is ready in Overlap_psi1_psi2_MPI3
!---------------------------------------------------------------------------------------
  SUBROUTINE Overlap_psipsi_MPI3(Overlap,psi1,psi2,With_Grid,Channel_ie)
    USE mod_system
    USE mod_psi_set_alloc
    USE mod_MPI_aux
    IMPLICIT NONE

    !-variables for the WP--------------------------------------------------------------
    TYPE(param_psi),intent(in)                :: psi1
    TYPE(param_psi),intent(in)                :: psi2
    Complex(kind=Rkind)                       :: Overlap
    Logical,optional,intent(in)               :: With_Grid
    Integer,optional,intent(in)               :: Channel_ie

    !-working variables-----------------------------------------------------------------
    Logical                                   :: With_Grid_loc
    Complex(kind=Rkind)                       :: temp
    Real(kind=Rkind)                          :: WrhonD
    Real(kind=Rkind)                          :: Roverlap
    Real(kind=Rkind)                          :: Rtemp
    Real(kind=Rkind),allocatable              :: wrho(:)
    Integer                                   :: locChannel_ie
    Integer                                   :: i_qa
    Integer                                   :: i_qaie
    Integer                                   :: i_be
    Integer                                   :: i_bi
    Integer                                   :: i_ba
    Integer                                   :: i_baie
    Integer                                   :: f_baie
    Integer                                   :: iie
    Integer                                   :: fie
    Integer                                   :: iie_MPI
    Integer                                   :: fie_MPI

#if(run_MPI)

    !-for debuging----------------------------------------------------------------------
    character (len=*), parameter :: name_sub='Overlap_psipsi_MPI3'
    logical,parameter :: debug = .FALSE.
    ! logical,parameter :: debug = .TRUE.

    !-----------------------------------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) 'psi1'
      CALL ecri_psi(psi=psi1)

      write(out_unitp,*) 'psi2'
      CALL ecri_psi(psi=psi2)
      write(out_unitp,*) 'GridRep,BasisRep ?'
      IF (present(With_Grid)) write(out_unitp,*) 'With_Grid',With_Grid
      IF (present(Channel_ie)) write(out_unitp,*) 'Channel_ie',Channel_ie
    ENDIF
    !-----------------------------------------------------------------------------------

    With_Grid_loc = .FALSE.

    IF(present(With_Grid)) With_Grid_loc=With_Grid

    locChannel_ie = 0
    IF(present(Channel_ie)) locChannel_ie=Channel_ie

    IF (psi1%nb_baie>psi1%nb_tot) THEN
      With_Grid_loc = .FALSE.
    ENDIF

    ! get bound1_MPI and bound2_MPI
    IF(.NOT. With_Grid_loc) THEN
      nb_per_MPI=psi1%nb_tot/MPI_np
      nb_rem_MPI=mod(psi1%nb_tot,MPI_np)
    ELSE
      nb_per_MPI=psi1%nb_qaie/MPI_np
      nb_rem_MPI=mod(psi1%nb_qaie,MPI_np) 
    ENDIF
    bound1_MPI=MPI_id*nb_per_MPI+1+MIN(MPI_id,nb_rem_MPI)
    bound2_MPI=(MPI_id+1)*nb_per_MPI+MIN(MPI_id,nb_rem_MPI)+merge(1,0,nb_rem_MPI>MPI_id)

    ! With_Grid_loc: F
    IF(With_Grid_loc) THEN
      IF(psi1%cplx .AND. allocated(psi1%CvecG) .AND. allocated(psi2%CvecG)) THEN
      ELSE IF(.NOT. psi1%cplx .AND. allocated(psi1%RvecG) .AND.                        &
               allocated(psi2%RvecG)) THEN
      ELSE
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' impossible to calculate the GridRep overlap'
        write(out_unitp,*) ' With_Grid_loc=t but problem with the allocation GridRep'
        write(out_unitp,*) 'allocated(psi1%CvecG)',allocated(psi1%CvecG)
        write(out_unitp,*) 'allocated(psi2%CvecG)',allocated(psi2%CvecG)
        write(out_unitp,*) 'allocated(psi1%RvecG)',allocated(psi1%RvecG)
        write(out_unitp,*) 'allocated(psi2%RvecG)',allocated(psi2%RvecG)
        write(out_unitp,*) ' psi1'
        CALL ecri_psi(psi=psi1,ecri_GridRep=.TRUE.)
        write(out_unitp,*) ' psi2'
        CALL ecri_psi(psi=psi2,ecri_GridRep=.TRUE.)
        STOP
      ENDIF
    ELSE
      IF(psi1%cplx .AND.allocated(psi1%CvecB) .AND. allocated(psi2%CvecB)) THEN
      ELSE IF(.NOT. psi1%cplx .AND. allocated(psi1%RvecB) .AND.                        &
              allocated(psi2%RvecB)) THEN
      ELSE
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' impossible to calculate the BasisRep overlap'
        write(out_unitp,*) ' With_Grid_loc=f (on basis) but problem with the allocation of BasisRep'
        write(out_unitp,*) 'allocated(psi1%CvecB)',allocated(psi1%CvecB)
        write(out_unitp,*) 'allocated(psi2%CvecB)',allocated(psi2%CvecB)
        write(out_unitp,*) 'allocated(psi1%RvecB)',allocated(psi1%RvecB)
        write(out_unitp,*) 'allocated(psi2%RvecB)',allocated(psi2%RvecB)
        write(out_unitp,*) ' psi1'
        CALL ecri_psi(psi=psi1,ecri_BasisRep=.TRUE.)
        write(out_unitp,*) ' psi2'
        CALL ecri_psi(psi=psi2,ecri_BasisRep=.TRUE.)
        STOP
      ENDIF
    ENDIF

    Overlap = cmplx(ZERO,ZERO,kind=Rkind)
    IF(.NOT. With_Grid_loc) THEN
      i_baie=1
      f_baie=psi1%nb_tot
      IF(psi1%nb_tot==psi1%nb_baie .AND.  locChannel_ie>0 .AND.           &
         locChannel_ie <= psi1%nb_bi*psi1%nb_be) THEN
        i_baie = 1 + (locChannel_ie-1)*psi1%nb_ba
        f_baie = i_baie-1 + psi1%nb_ba
      END IF

      IF(bound2_MPI>i_baie .AND. bound1_MPI<f_baie) THEN
        f_baie=MIN(f_baie,bound2_MPI)
        i_baie=MAX(i_baie,bound1_MPI)
        
        IF(psi1%symab>-1 .AND. psi2%symab>-1 .AND. psi1%symab/=psi2%symab) THEN
          !Overlap = cmplx(ZERO,ZERO,kind=Rkind)
        ELSE
          IF(psi1%cplx) THEN
            Overlap=dot_product(psi1%CvecB(i_baie:f_baie),psi2%CvecB(i_baie:f_baie))
          ELSE
            ROverlap=dot_product(psi1%RvecB(i_baie:f_baie),psi2%RvecB(i_baie:f_baie))
            Overlap=cmplx(ROverlap,ZERO,kind=Rkind)
          ENDIF
        ENDIF
      ENDIF ! bound2_MPI>i_baie .AND. bound1_MPI<f_baie

    ELSE ! With_Grid_loc
    
      CALL alloc_NParray(wrho,[psi1%nb_qa],"wrho",name_sub)
      DO i_qa=1,psi1%nb_qa
        wrho(i_qa) = Rec_WrhonD(psi1%BasisnD,i_qa)
      ENDDO

      IF(psi1%cplx) THEN
        iie=1
        fie=psi1%nb_qa
        DO i_be=1,psi1%nb_be
        DO i_bi=1,psi1%nb_bi
          IF(bound2_MPI>iie .AND. bound1_MPI<fie) THEN
            fie_MPI=MIN(fie,bound2_MPI)
            iie_MPI=MAX(iie,bound1_MPI)
            Overlap=Overlap+dot_product(psi1%CvecG(iie_MPI:fie_MPI),                   &
                            wrho(Mod(iie_MPI,psi1%nb_qa):Mod(fie_MPI,psi1%nb_qa))      &
                            *psi2%CvecG(iie_MPI:fie_MPI))
          ENDIF
          iie=iie+psi1%nb_qa
          fie=fie+psi1%nb_qa
        ENDDO
        ENDDO
      ELSE
        iie=1
        fie=psi1%nb_qa
        DO i_be=1,psi1%nb_be
        DO i_bi=1,psi1%nb_bi
          IF(bound2_MPI>iie .AND. bound1_MPI<fie) THEN
            fie_MPI=MIN(fie,bound2_MPI)
            iie_MPI=MAX(iie,bound1_MPI)
            ROverlap=ROverlap+dot_product(psi1%RvecG(iie_MPI:fie_MPI),                 &
                              wrho(Mod(iie_MPI,psi1%nb_qa):Mod(fie_MPI,psi1%nb_qa))    &
                              *psi2%RvecG(iie_MPI:fie_MPI))
          ENDIF
          iie=iie + psi1%nb_qa
          fie=fie + psi1%nb_qa
        ENDDO
        ENDDO
        Overlap=cmplx(ROverlap,ZERO,kind=Rkind)
      ENDIF

      CALL dealloc_NParray(wrho,"wrho",name_sub)

    ENDIF

    !-----------------------------------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'Overlap : ',Overlap
      write(out_unitp,*) 'END ',name_sub
    ENDIF
    !-----------------------------------------------------------------------------------

#endif
  END SUBROUTINE Overlap_psipsi_MPI3
!=======================================================================================

!=======================================================================================
! this is a temp subroutine for running Overlap_psi1_psi2 on all threads with MPI
! it will be replaced by the original 'Overlap_psi1_psi2' 
! and thus completely removed later
!---------------------------------------------------------------------------------------
  SUBROUTINE Overlap_psipsi_MPI(Overlap,psi1,psi2,With_Grid,Channel_ie)
    USE mod_system
    USE mod_psi_set_alloc
    IMPLICIT NONE

!----- variables for the WP ----------------------------------------
    TYPE (param_psi), intent(in)    :: psi1,psi2
    complex (kind=Rkind)            :: Overlap
    logical, optional, intent(in)   :: With_Grid
    integer, optional, intent(in)   :: Channel_ie

!------ working variables ---------------------------------
    logical              :: With_Grid_loc
    integer              :: locChannel_ie
    integer              :: i_qa,i_qaie
    integer              :: i_be,i_bi,i_ba
    integer              :: i_baie,f_baie
    integer              :: i_modif_q
    real (kind=Rkind)    :: WrhonD
    complex (kind=Rkind) :: temp
    real (kind=Rkind)    :: Roverlap,Rtemp
    integer              :: iie,fie
    real (kind=Rkind), allocatable :: wrho(:)

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Overlap_psipsi_MPI'
    logical,parameter :: debug = .FALSE.
!     logical,parameter :: debug = .TRUE.

#if(run_MPI)

!-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) 'psi1'
      CALL ecri_psi(psi=psi1)

      write(out_unitp,*) 'psi2'
      CALL ecri_psi(psi=psi2)
      write(out_unitp,*) 'GridRep,BasisRep ?'
      IF (present(With_Grid)) write(out_unitp,*) 'With_Grid',With_Grid
      IF (present(Channel_ie)) write(out_unitp,*) 'Channel_ie',Channel_ie
    END IF
!-----------------------------------------------------------

    With_Grid_loc = .FALSE.

    IF (present(With_Grid)) With_Grid_loc = With_Grid

    locChannel_ie = 0
    IF (present(Channel_ie)) locChannel_ie = Channel_ie

    IF (psi1%nb_baie > psi1%nb_tot) THEN
      With_Grid_loc = .FALSE.
    END IF

    ! With_Grid_loc: F
    IF (With_Grid_loc) THEN
      IF (psi1%cplx .AND.                                             &
       allocated(psi1%CvecG) .AND. allocated(psi2%CvecG)) THEN
      ELSE IF (.NOT. psi1%cplx .AND.                                  &
       allocated(psi1%RvecG) .AND. allocated(psi2%RvecG)) THEN
      ELSE
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' impossible to calculate the GridRep overlap'
        write(out_unitp,*) ' With_Grid_loc=t but problem with the allocation GridRep'
        write(out_unitp,*) 'allocated(psi1%CvecG)',allocated(psi1%CvecG)
        write(out_unitp,*) 'allocated(psi2%CvecG)',allocated(psi2%CvecG)
        write(out_unitp,*) 'allocated(psi1%RvecG)',allocated(psi1%RvecG)
        write(out_unitp,*) 'allocated(psi2%RvecG)',allocated(psi2%RvecG)
        write(out_unitp,*) ' psi1'
        CALL ecri_psi(psi=psi1,ecri_GridRep=.TRUE.)
        write(out_unitp,*) ' psi2'
        CALL ecri_psi(psi=psi2,ecri_GridRep=.TRUE.)
        STOP
      END IF
    ELSE
      IF (psi1%cplx .AND.                                             &
       allocated(psi1%CvecB) .AND. allocated(psi2%CvecB)) THEN
      ELSE IF (.NOT. psi1%cplx .AND.                                  &
       allocated(psi1%RvecB) .AND. allocated(psi2%RvecB)) THEN
      ELSE
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' impossible to calculate the BasisRep overlap'
        write(out_unitp,*) ' With_Grid_loc=f (on basis) but problem with the allocation of BasisRep'
        write(out_unitp,*) 'allocated(psi1%CvecB)',allocated(psi1%CvecB)
        write(out_unitp,*) 'allocated(psi2%CvecB)',allocated(psi2%CvecB)
        write(out_unitp,*) 'allocated(psi1%RvecB)',allocated(psi1%RvecB)
        write(out_unitp,*) 'allocated(psi2%RvecB)',allocated(psi2%RvecB)
        write(out_unitp,*) ' psi1'
        CALL ecri_psi(psi=psi1,ecri_BasisRep=.TRUE.)
        write(out_unitp,*) ' psi2'
        CALL ecri_psi(psi=psi2,ecri_BasisRep=.TRUE.)
        STOP
      END IF
    END IF

    IF (.NOT. With_Grid_loc) THEN
      i_baie=1
      f_baie=psi1%nb_tot
      IF (psi1%nb_tot == psi1%nb_baie .AND.  locChannel_ie > 0 .AND.  &
                              locChannel_ie <= psi1%nb_bi*psi1%nb_be) THEN
        i_baie = 1 + (locChannel_ie-1)*psi1%nb_ba
        f_baie = i_baie-1 + psi1%nb_ba
      END IF
      IF (psi1%symab > -1 .AND. psi2%symab > -1 .AND. psi1%symab /= psi2%symab) THEN
        Overlap = cmplx(ZERO,ZERO,kind=Rkind)
      ELSE
        IF (psi1%cplx) THEN
          Overlap = dot_product( psi1%CvecB(i_baie:f_baie) ,          &
                                 psi2%CvecB(i_baie:f_baie) )
        ELSE
          ROverlap = dot_product( psi1%RvecB(i_baie:f_baie) ,         &
                                  psi2%RvecB(i_baie:f_baie) )
          Overlap = cmplx(ROverlap,ZERO,kind=Rkind)
        END IF
      END IF

    ELSE

!       - initialization ----------------------------------
      Overlap = cmplx(ZERO,ZERO,kind=Rkind)

      CALL alloc_NParray(wrho,[psi1%nb_qa],"wrho",name_sub)
      DO i_qa=1,psi1%nb_qa
        wrho(i_qa) = Rec_WrhonD(psi1%BasisnD,i_qa)
      END DO

      IF (psi1%cplx) THEN
        iie = 1
        fie = psi1%nb_qa
        DO i_be=1,psi1%nb_be
        DO i_bi=1,psi1%nb_bi
          Overlap = Overlap + dot_product(                            &
            psi1%CvecG(iie:fie),wrho*psi2%CvecG(iie:fie))
          iie = iie + psi1%nb_qa
          fie = fie + psi1%nb_qa
        END DO
        END DO
      ELSE
        iie = 1
        fie = psi1%nb_qa
        DO i_be=1,psi1%nb_be
        DO i_bi=1,psi1%nb_bi
          Overlap = Overlap + cmplx(dot_product(                      &
            psi1%RvecG(iie:fie),wrho*psi2%RvecG(iie:fie)) ,kind=Rkind)
          iie = iie + psi1%nb_qa
          fie = fie + psi1%nb_qa
        END DO
        END DO
      END IF

      CALL dealloc_NParray(wrho,"wrho",name_sub)

    END IF

!----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'Overlap : ',Overlap
      write(out_unitp,*) 'END ',name_sub
    END IF
!----------------------------------------------------------

#endif
  END SUBROUTINE Overlap_psipsi_MPI
!=======================================================================================


!=======================================================================================
! Save vectors
!---------------------------------------------------------------------------------------
  SUBROUTINE sub_LCpsi_TO_psi_MPI(psi,Vec,ndim,nb_save)
    USE mod_system
    USE mod_psi_set_alloc
    USE mod_MPI_aux
    IMPLICIT NONE

    TYPE(param_psi),            intent(inout) :: psi(ndim)
    Real(kind=Rkind),           intent(in)    :: Vec(ndim,ndim)
    Integer,                    intent(in)    :: ndim
    Integer,                    intent(in)    :: nb_save

    Real(kind=Rkind),allocatable              :: PsiRk(:)
    Integer                                   :: symab_psi_old(ndim)
    Integer                                   :: isym
    Integer                                   :: ii
    Integer                                   :: kk

#if(run_MPI)

    IF(sum(abs(Vec))>ONETENTH**10) THEN
      IF(psi(1)%BasisRep) THEN
        DO ii=1,ndim
          symab_psi_old(ii)=psi(ii)%symab
        END DO

        CALL allocate_array(PsiRk,1,ndim)

        DO kk=bounds_MPI(1,MPI_id),bounds_MPI(2,MPI_id)
          DO ii=1,ndim
            PsiRk(ii)=psi(ii)%RvecB(kk)
          END DO

          PsiRk(:)=matmul(PsiRk,Vec)

          DO ii=1,ndim
            psi(ii)%RvecB(kk)=PsiRk(ii)
          END DO
        END DO

        IF(allocated(PsiRk)) deallocate(PsiRk)

        DO ii=1,ndim
          CALL MPI_combine_array(psi(ii)%RvecB,MS=MPI_scheme)
        ENDDO

        DO ii=1,ndim
          isym=maxloc(abs(Vec(:,ii)),dim=1)
          CALL Set_symab_OF_psiBasisRep_MPI(psi(ii),symab_psi_old(isym))
        END DO

      ELSE !----------------------------------------------------------------------------

        CALL allocate_array(PsiRk,1,ndim)

        DO kk=bounds_MPI(1,MPI_id),bounds_MPI(2,MPI_id)
          DO ii=1,ndim
            PsiRk(ii)=psi(ii)%RvecG(kk)
          END DO

          PsiRk(:)=matmul(PsiRk,Vec)

          DO ii=1,ndim
            psi(ii)%RvecG(kk)=PsiRk(ii)
          END DO
        END DO

        IF(allocated(PsiRk)) deallocate(PsiRk)

        DO ii=1,ndim
          CALL MPI_combine_array(psi(ii)%RvecG,MS=MPI_scheme)
        ENDDO

      ENDIF ! for psi(1)%BasisRep

    ELSE !------------------------------------------------------------------------------
      CONTINUE ! nothing!!!
    END IF

#endif
  END SUBROUTINE sub_LCpsi_TO_psi_MPI
!=======================================================================================

ENDMODULE mod_psi_Op_MPI

