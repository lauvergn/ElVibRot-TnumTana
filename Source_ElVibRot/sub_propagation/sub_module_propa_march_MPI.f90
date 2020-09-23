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

MODULE mod_march_MPI
  USE mod_system
  USE mod_psi,  ONLY:param_psi,alloc_NParray,dealloc_NParray,dealloc_psi
  USE mod_propa,ONLY:param_propa,Calc_AutoCorr,Write_AutoCorr
  USE mod_march_SG4
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: march_SIL_MPI

  CONTAINS
!=======================================================================================
!> @brief MPI version of march_SIL
!> according to Smolyak rep. 
!> Note, this works for massive cluster with big memmory
!
! require a new action subroutine, SR vector ready in V
! rewrite the atcion function to make it easier for callimng from different subroutines
!=======================================================================================  
  SUBROUTINE march_SIL_MPI(TT,no,psi,psi0,para_H,para_propa)
    USE mod_Op,                         ONLY:param_Op,sub_OpPsi
    USE mod_psi,                        ONLY:norm2_psi
    USE mod_psi_Op_MPI,                 ONLY:norm2_psi_SR_MPI,Overlap_psi1_psi2_SRB_MPI
    USE mod_psi_set_alloc,              ONLY:param_psi,psi_times_SR_MPI
    USE mod_OpPsi_MPI,                  ONLY:sub_scaledOpPsi_SR_MPI
    USE mod_propa_MPI,                  ONLY:Calc_AutoCorr_SR_MPI 
    USE mod_basis_BtoG_GtoB_SGType4,    ONLY:TypeRVec
    USE mod_basis_BtoG_GtoB_SGType4_MPI,ONLY:ini_iGs_MPI
    IMPLICIT NONE

    !-variables for namelist minimum----------------------------------------------------
    TYPE(param_Op),                 intent(in)    :: para_H

    !-variables for the WP propagation--------------------------------------------------
    TYPE(param_propa),              intent(inout) :: para_propa
    TYPE(param_psi),                intent(inout) :: psi0
    TYPE(param_psi),                intent(inout) :: psi
    Real(kind=Rkind),               intent(in)    :: TT
    Integer,                        intent(in)    :: no
    
    TYPE(param_psi),allocatable                   :: tab_KrylovSpace(:)
    TYPE(param_psi)                               :: w1
    Complex(kind=Rkind)      :: H(para_propa%para_poly%npoly,para_propa%para_poly%npoly)
    Complex(kind=Rkind)               :: psi0_psiKrylovSpace(para_propa%para_poly%npoly)
    Complex(kind=Rkind),allocatable               :: UPsiOnKrylov(:)
    Complex(kind=Rkind),allocatable               :: Vec(:,:)
    Complex(kind=Rkind)                           :: Overlap
    Complex(kind=Rkind)                           :: cdot
    Real(kind=Rkind),allocatable                  :: Eig(:)
    Real(kind=Rkind)                              :: E0
    Real(kind=Rkind)                              :: phase
    Real(kind=Rkind)                              :: micro_deltaT
    Real(kind=Rkind)                              :: micro_T
    Real(kind=Rkind)                              :: micro_phase
    Integer                                       :: n 
    Integer                                       :: k

    Character(len=*),parameter                    :: name_sub='march_SIL_MPI'

#if(run_MPI)

    !> for Krylov Space
    allocate(tab_KrylovSpace(para_propa%para_poly%npoly+1)) ! param_psi 
    
    !> initialize iGs_MPI
    IF(para_H%BasisnD%para_SGType2%once_action) CALL ini_iGs_MPI(para_H%BasisnD,.TRUE.)
    
    !> Extract compact psi to each threads and transfer to gird rep.
    CALL SmolyakR_distribute_SRB_MPI(psi ,para_Op=para_H)
    IF(para_propa%nb_micro>1) CALL SmolyakR_distribute_SRB_MPI(psi0,para_Op=para_H)
    
    !> set first Krylov Space
    tab_KrylovSpace(1)=psi 
    IF(para_propa%nb_micro>1) psi0_psiKrylovSpace(:)=CZERO

    E0=para_H%E0
    H(:,:)=CZERO ! initialize H
    DO k=1,para_propa%para_poly%npoly
      IF(para_propa%nb_micro>1) THEN
        ! <psi(t)|v_k'> ~ <v_0|v_k'>
        psi0_psiKrylovSpace(k)=Calc_AutoCorr_SR_MPI(psi0,tab_KrylovSpace(k),           &
                                                    para_propa,TT,Write_AC=.FALSE.)
      END IF

      ! |w1>=H|v_k>
      CALL sub_OpPsi(Psi=tab_KrylovSpace(k),OpPsi=w1,para_Op=para_H)
      
      ! Warning:require this two line to get same result as normal case. 
      ! waiting for new idea
      !CALL SmolyakR_to_packedB_SRB_MPI(w1,para_Op=para_H)
      !CALL SmolyakR_distribute_SRB_MPI(w1 ,para_Op=para_H)

      ! case k=1, alpha_1=<v_1|H|v_1>
      IF(k==1) THEN 
        !> Energy shift, E0, calculation for the first iteration. E0=<psi |H| psi> 
        !> This shift is important to improve the stapility
        !> Note phase need to be taking into account at the end of the iterations
        CALL Overlap_psi1_psi2_SRB_MPI(Overlap,tab_KrylovSpace(1),w1)
        !CALL MPI_Reduce_sum_Bcast(Overlap) ! reduce sum and boardcast
        E0=real(Overlap,kind=Rkind)
      ENDIF
      ! scaling with E0 
      CALL sub_scaledOpPsi_SR_MPI(Psi=tab_KrylovSpace(k),OpPsi=w1,E0=E0,Esc=ONE)

      ! v_{k}=w_{k-1}/||w_{k-1}||=w_{k-1}/beta_{k-1}
      ! k=1: w_{k}=H|v_{k}>-alpha_k|v_{k}>
      ! k>1: w_{k}=H|v_{k}>-alpha_k|v_{k}>-beta_{k-1}|v_{k-1}>
      ! beta_{k}=<v_k|H|v_{k+1}>=<v_{k+1}|H|v_{k}>=H(k,k+1)=H(k+1,k)
      IF(k>1) THEN
        !w1=w1-H(k,k-1)*tab_KrylovSpace(k-1) !> H|v_{k}>-beta_{k-1}|v_{k-1}>
        CALL w1_minus_Cv_SR_MPI(w1,H(k,k-1),tab_KrylovSpace(k-1))
      ENDIF
      ! alpha_k=<v_k|H|v_k>
      CALL Overlap_psi1_psi2_SRB_MPI(Overlap,w1,tab_KrylovSpace(k))
      !CALL MPI_Reduce_sum_Bcast(Overlap)
      H(k,k)=Overlap ! alpha_k

      !> diagonalization then exit when:
      !> 1. k reaches npoly 2. the scheme converges
      ! psi(t+dt)=U|psi(t)>=sum_{l=1}^n <vec_l|psi(t)> e^{-i*Ei*dt}|vec_l>
      !                    =sum_{k=1}^n a_k|v_k>
      !                 a_k=sum_{l=1}^n<vec_l|v0>e^{-i*Ei*dt}<v_s|vec_l>
      ! works on all the threads
      CALL UPsi_spec_MPI(UPsiOnKrylov,H(1:k,1:k),Vec,Eig,                              &
                                                para_propa%WPdeltaT,k,With_diago=.TRUE.)

      IF(abs(UPsiOnKrylov(k))<para_propa%para_poly%poly_tol .OR.                       &
         k==para_propa%para_poly%npoly) THEN
        n=k
        write(out_unitp,*) n,'abs(UPsiOnKrylov(n)',abs(UPsiOnKrylov(n))
        EXIT
      END IF

      !> w_k=(H|v_{k}>-beta_{k-1}|v_{k-1}>)-alpha_k|v_{k}>="w1"-alpha_k|v_k>
      ! w1=w1-Overlap*tab_KrylovSpace(k) !< now "w1" is w_k
      CALL w1_minus_Cv_SR_MPI(w1,Overlap,tab_KrylovSpace(k))

      !> beta_{k}=||w_k||=sqrt(<w_k|w_k>)
      CALL norm2_psi_SR_MPI(w1,2) !< 2: just calculate the normalization constant 
      !CALL MPI_Reduce_sum_Bcast(w1%norm2)
      H(k+1,k)=sqrt(w1%norm2)
      H(k,k+1)=conjg(H(k+1,k))

      !> assign v_{k+1}=w_k/beta_k
      CALL norm2_psi_SR_MPI(w1,3) !< 3: normalize SR_B with existing norm. constant
      tab_KrylovSpace(k+1)=w1    
    ENDDO ! for k=1,para_propa%para_poly%npoly
    !CALL SRB_to_packB_write(w1,call_from='endofSILloop')

    psi=ZERO
    DO k=1,n
      ! psi(t+dt)=sum_{k=1}^n a_k|v_k>
      ! psi=psi+UPsiOnKrylov(k)*tab_KrylovSpace(k)
      CALL w1_minus_Cv_SR_MPI(psi,-UPsiOnKrylov(k),tab_KrylovSpace(k))
    END DO

    !> check normalization
    CALL norm2_psi_SR_MPI(psi,2) !< 2: just calculate the normalization constant 
    !CALL norm2_psi(psi)

    IF(psi%norm2 > para_propa%max_norm2) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' STOP propagation: norm > max_norm',psi%norm2
      para_propa%march_error  =.TRUE.
      para_propa%test_max_norm=.TRUE.
      STOP
    ENDIF

    ! Phase Shift
    phase=E0*para_propa%WPdeltaT  
    !psi=psi*exp(-EYE*phase)
    CALL psi_times_SR_MPI(psi,exp(-EYE*phase))
    
    ! transfer back to packed basis, valid on master only
    CALL SmolyakR_to_packedB_SRB_MPI(psi,para_Op=para_H)

    ! autocorelation
    IF(para_propa%nb_micro>1) THEN
      micro_deltaT=para_propa%WPdeltaT/real(para_propa%nb_micro,kind=Rkind)
      micro_phase =phase/real(para_propa%nb_micro,kind=Rkind)
    
      phase =ZERO
      micro_T=ZERO
      
      DO k=1,para_propa%nb_micro
        micro_T=micro_T+micro_deltaT
        phase=phase+micro_phase
        
        CALL UPsi_spec_MPI(UPsiOnKrylov,H,Vec,Eig,micro_T,n,With_diago=.FALSE.)
        ! sum(a_k <v0|v_k'>)
        cdot=sum(UPsiOnKrylov(1:n)*psi0_psiKrylovSpace(1:n)) ! we cannot use dot_product
        cdot=cdot*exp(-EYE*phase)
        IF(MPI_id==0) CALL Write_AutoCorr(no,TT+micro_T,cdot)
      ENDDO 
      CALL flush_perso(no)
    ELSE
      !cdot=Calc_AutoCorr_SR_MPI(psi0,psi,para_propa,TT,Write_AC=.FALSE.)
      IF(MPI_id==0) THEN
        cdot=Calc_AutoCorr(psi0,psi,para_propa,TT,Write_AC=.FALSE.)
        CALL Write_AutoCorr(no,TT+para_propa%WPdeltaT,cdot)
        CALL flush_perso(no)
      ENDIF
    ENDIF

    ! transfer back to packed basis, valid on master only
    !CALL SmolyakR_to_packedB_SRB_MPI(psi,para_Op=para_H)
    
    ! deallocation  
    DO k=1,size(tab_KrylovSpace)
       CALL dealloc_psi(tab_KrylovSpace(k)) ! deallocation of SR_B included
    ENDDO
    deallocate(tab_KrylovSpace)

    IF(allocated(Vec))          CALL dealloc_NParray(Vec,'Vec',name_sub)
    IF(allocated(Eig))          CALL dealloc_NParray(Eig,'Eig',name_sub)
    IF(allocated(UPsiOnKrylov)) CALL dealloc_NParray(UPsiOnKrylov,'UPsiOnKrylov',name_sub)
    CALL dealloc_psi(w1)

#endif
  END SUBROUTINE march_SIL_MPI
!=======================================================================================


!=======================================================================================
!> transfer from Smolyak basis to packed basis
!=======================================================================================
  SUBROUTINE SRB_to_packB_write(psi,call_from)
    USE mod_Op,                     ONLY:param_Op
    USE mod_psi_set_alloc,          ONLY:param_psi
    USE mod_basis_BtoG_GtoB_SGType4,ONLY:tabR_AT_iG_TO_tabPackedBasis
    USE mod_psi_Op_MPI,             ONLY:Overlap_psi1_psi2_SRB_MPI
    USE mod_MPI_aux
    IMPLICIT NONE

    TYPE(param_psi),                intent(inout) :: psi
    Character(len=*),               intent(in)    :: call_from
    
    Complex(kind=Rkind)                           :: Overlap
    Real(kind=Rkind),allocatable                  :: SRB(:)
    Real(kind=Rkind),allocatable                  :: CvecR(:)
    Real(kind=Rkind),allocatable                  :: CvecC(:)
    Integer                                       :: iG
    Integer                                       :: d1
    Integer                                       :: d2

#if(run_MPI)

    CALL allocate_array(CvecR,1,size(psi%CvecB))
    CALL allocate_array(CvecC,1,size(psi%CvecB))
    DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)
      d1=psi%SR_B_index(iG)
      d2=psi%SR_B_index(iG+1)-1
      CALL allocate_array(SRB,1,d2-d1+1)
      SRB=psi%SR_B(d1:d2,1) 
      CALL tabR_AT_iG_TO_tabPackedBasis(CvecR,SRB,iG,psi%BasisnD%para_SGType2,         &
                                        psi%BasisnD%WeightSG(iG)) 
      SRB=psi%SR_B(d1:d2,2)
      CALL tabR_AT_iG_TO_tabPackedBasis(CvecC,SRB,iG,psi%BasisnD%para_SGType2,         &
                                        psi%BasisnD%WeightSG(iG)) 
    ENDDO
    
    CALL Overlap_psi1_psi2_SRB_MPI(Overlap,psi,psi)
    
    write(out_unitp,*) call_from, ' Real part: ',CvecR
    write(out_unitp,*) call_from, ' imag part: ',CvecC
    write(out_unitp,*) call_from, ' overlap  : ',Overlap

#endif
  END SUBROUTINE SRB_to_packB_write
!=======================================================================================


!=======================================================================================
!> @brief MPI version of march_SIL
!> according to Smolyak rep. 
!> Note, this works for massive cluster with big memmory
!
! require a new action subroutine, SR vector ready in V
! rewrite the atcion function to make it easier for callimng from different subroutines
!=======================================================================================  
  SUBROUTINE march_SIL_MPI_old(TT,no,psi,psi0,para_H,para_propa)
    USE mod_Op,                         ONLY:param_Op,sub_OpPsi
    USE mod_OpPsi_MPI,                  ONLY:sub_scaledOpPsi_SR_MPI
    USE mod_psi_set_alloc,              ONLY:param_psi,psi_times_SR_MPI
    USE mod_psi_Op_MPI,                 ONLY:norm2_psi_SR_MPI,Overlap_psi1_psi2_SRG_MPI
    USE mod_propa_MPI,                  ONLY:Calc_AutoCorr_SR_MPI 
    USE mod_basis_BtoG_GtoB_SGType4,    ONLY:TypeRVec
    USE mod_basis_BtoG_GtoB_SGType4_MPI,ONLY:ini_iGs_MPI
    IMPLICIT NONE

    !-variables for namelist minimum----------------------------------------------------
    TYPE(param_Op),                 intent(in)    :: para_H

    !-variables for the WP propagation--------------------------------------------------
    TYPE(param_propa),              intent(inout) :: para_propa
    TYPE(param_psi),                intent(inout) :: psi0
    TYPE(param_psi),                intent(inout) :: psi
    Real(kind=Rkind),               intent(in)    :: TT
    Integer,                        intent(in)    :: no
    
    TYPE(param_psi),allocatable                   :: tab_KrylovSpace(:)
    TYPE(param_psi)                               :: w1
    Complex(kind=Rkind)      :: H(para_propa%para_poly%npoly,para_propa%para_poly%npoly)
    Complex(kind=Rkind)               :: psi0_psiKrylovSpace(para_propa%para_poly%npoly)
    Complex(kind=Rkind),allocatable               :: UPsiOnKrylov(:)
    Complex(kind=Rkind),allocatable               :: Vec(:,:)
    Complex(kind=Rkind)                           :: Overlap
    Complex(kind=Rkind)                           :: cdot
    Real(kind=Rkind),allocatable                  :: Eig(:)
    Real(kind=Rkind)                              :: E0
    Real(kind=Rkind)                              :: phase
    Real(kind=Rkind)                              :: micro_deltaT
    Real(kind=Rkind)                              :: micro_T
    Real(kind=Rkind)                              :: micro_phase
    Integer                                       :: n 
    Integer                                       :: k

    Character(len=*),parameter                    :: name_sub='march_SIL_MPI'

#if(run_MPI)

    !> for Krylov Space
    allocate(tab_KrylovSpace(para_propa%para_poly%npoly+1)) ! param_psi 
    
    !> initialize iGs_MPI
    IF(para_H%BasisnD%para_SGType2%once_action) CALL ini_iGs_MPI(para_H%BasisnD,.TRUE.)
    
    !> Extract compact psi to each threads and transfer to gird rep.
    CALL SmolyakR_distribute_SRG_MPI(psi ,para_Op=para_H)
    CALL SmolyakR_distribute_SRG_MPI(psi0,para_Op=para_H)
    
    !> set first Krylov Space
    tab_KrylovSpace(1)=psi 
    IF(para_propa%nb_micro>1) psi0_psiKrylovSpace(:)=CZERO
    
    E0=para_H%E0
    H(:,:)=CZERO ! initialize H
    DO k=1,para_propa%para_poly%npoly
    
      IF(para_propa%nb_micro>1) THEN
        ! <psi(t)|v_k'> ~ <v_0|v_k'>
        psi0_psiKrylovSpace(k)=Calc_AutoCorr_SR_MPI(psi0,tab_KrylovSpace(k),           &
                                                    para_propa,TT,Write_AC=.FALSE.)
      END IF

      ! |w1>=H|v_k>
      CALL sub_OpPsi(Psi=tab_KrylovSpace(k),OpPsi=w1,para_Op=para_H)

      ! case k=1, alpha_1=<v_1|H|v_1>
      IF(k==1) THEN 
        !> Energy shift, E0, calculation for the first iteration. E0=<psi |H| psi> 
        !> This shift is important to improve the stapility
        !> Note phase need to be taking into account at the end of the iterations
        CALL Overlap_psi1_psi2_SRG_MPI(Overlap,tab_KrylovSpace(1),w1)
        !CALL MPI_Reduce_sum_Bcast(Overlap) ! reduce sum and boardcast
        E0=real(Overlap,kind=Rkind)
      ENDIF

      ! scaling with E0 
      CALL sub_scaledOpPsi_SR_MPI(Psi=tab_KrylovSpace(k),OpPsi=w1,E0=E0,Esc=ONE)
      
      ! v_{k}=w_{k-1}/||w_{k-1}||=w_{k-1}/beta_{k-1}
      ! k=1: w_{k}=H|v_{k}>-alpha_k|v_{k}>
      ! k>1: w_{k}=H|v_{k}>-alpha_k|v_{k}>-beta_{k-1}|v_{k-1}>
      ! beta_{k}=<v_k|H|v_{k+1}>=<v_{k+1}|H|v_{k}>=H(k,k+1)=H(k+1,k)
      IF(k>1) THEN
        !w1=w1-H(k,k-1)*tab_KrylovSpace(k-1) !> H|v_{k}>-beta_{k-1}|v_{k-1}>
        CALL w1_minus_Cv_SR_MPI(w1,H(k,k-1),tab_KrylovSpace(k-1))
      ENDIF

      ! alpha_k=<v_k|H|v_k>
      CALL Overlap_psi1_psi2_SRG_MPI(Overlap,w1,tab_KrylovSpace(k))
      !CALL MPI_Reduce_sum_Bcast(Overlap)
      H(k,k)=Overlap ! alpha_k

      !> diagonalization then exit when:
      !> 1. k reaches npoly 2. the scheme converges
      ! psi(t+dt)=U|psi(t)>=sum_{l=1}^n <vec_l|psi(t)> e^{-i*Ei*dt}|vec_l>
      !                    =sum_{k=1}^n a_k|v_k>
      !                 a_k=sum_{l=1}^n<vec_l|v0>e^{-i*Ei*dt}<v_s|vec_l>
      ! works on all the threads
      CALL UPsi_spec_MPI(UPsiOnKrylov,H(1:k,1:k),Vec,Eig,                              &
                                                para_propa%WPdeltaT,k,With_diago=.TRUE.)
      IF(abs(UPsiOnKrylov(k))<para_propa%para_poly%poly_tol .OR.                       &
         k==para_propa%para_poly%npoly) THEN
        n=k
        write(out_unitp,*) n,'abs(UPsiOnKrylov(n)',abs(UPsiOnKrylov(n))
        EXIT
      END IF

      !> w_k=(H|v_{k}>-beta_{k-1}|v_{k-1}>)-alpha_k|v_{k}>="w1"-alpha_k|v_k>
      ! w1=w1-Overlap*tab_KrylovSpace(k) !< now "w1" is w_k
      CALL w1_minus_Cv_SR_MPI(w1,Overlap,tab_KrylovSpace(k))

      !> beta_{k}=||w_k||=sqrt(<w_k|w_k>)
      CALL norm2_psi_SR_MPI(w1,2) !< 2: just calculate the normalization constant 
      !CALL MPI_Reduce_sum_Bcast(w1%norm2)
      H(k+1,k)=sqrt(w1%norm2)
      H(k,k+1)=H(k+1,k)

      !> assign v_{k+1}=w_k/beta_k
      CALL norm2_psi_SR_MPI(w1,3) !< 3: normalize SR_G with existing norm. constant
      tab_KrylovSpace(k+1)=w1    
    ENDDO ! for k=1,para_propa%para_poly%npoly
    write(out_unitp,*) 'abs(UPsiOnKrylov)',abs(UPsiOnKrylov(1:n)), 'from ',MPI_id

    psi=ZERO
    DO k=1,n
      ! psi(t+dt)=sum_{k=1}^n a_k|v_k>
      ! psi=psi+UPsiOnKrylov(k)*tab_KrylovSpace(k)
      CALL psi_times_SR_MPI(tab_KrylovSpace(k),UPsiOnKrylov(k),psi)
    END DO

    !> check normalization
    CALL norm2_psi_SR_MPI(psi,2) !< 2: just calculate the normalization constant 
    IF(psi%norm2 > para_propa%max_norm2) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' STOP propagation: norm > max_norm',psi%norm2
      para_propa%march_error  =.TRUE.
      para_propa%test_max_norm=.TRUE.
      STOP
    ENDIF

    ! Phase Shift
    phase=E0*para_propa%WPdeltaT  
    ! psi=psi*exp(-EYE*phase)
    CALL psi_times_SR_MPI(psi,exp(-EYE*phase))

    ! autocorelation
    IF(para_propa%nb_micro>1) THEN
      micro_deltaT=para_propa%WPdeltaT/real(para_propa%nb_micro,kind=Rkind)
      micro_phase =phase/real(para_propa%nb_micro,kind=Rkind)
    
      phase =ZERO
      micro_T=ZERO
      
      DO k=1,para_propa%nb_micro
        micro_T=micro_T+micro_deltaT
        phase=phase+micro_phase
        
        CALL UPsi_spec_MPI(UPsiOnKrylov,H,Vec,Eig,micro_T,n,With_diago=.FALSE.)
        ! sum(a_k <v0|v_k'>)
        cdot=sum(UPsiOnKrylov(1:n)*psi0_psiKrylovSpace(1:n)) ! we cannot use dot_product
        cdot=cdot*exp(-EYE*phase)
        CALL Write_AutoCorr(no,TT+micro_T,cdot)
      ENDDO 
      CALL flush_perso(no)
    ELSE
      cdot=Calc_AutoCorr_SR_MPI(psi0,psi,para_propa,TT,Write_AC=.FALSE.)
      CALL Write_AutoCorr(no,TT+para_propa%WPdeltaT,cdot)
      CALL flush_perso(no)
    ENDIF

    ! deallocation  
    DO k=1,size(tab_KrylovSpace)
       CALL dealloc_psi(tab_KrylovSpace(k)) ! deallocation of SR_G included
    ENDDO
    deallocate(tab_KrylovSpace)

    IF(allocated(Vec))          CALL dealloc_NParray(Vec,'Vec',name_sub)
    IF(allocated(Eig))          CALL dealloc_NParray(Eig,'Eig',name_sub)
    IF(allocated(UPsiOnKrylov)) CALL dealloc_NParray(UPsiOnKrylov,'UPsiOnKrylov',name_sub)
    CALL dealloc_psi(w1)

    ! update 
    CALL SmolyakR_to_packedB_SRG_MPI(psi,para_Op=para_H)

#endif
  END SUBROUTINE march_SIL_MPI_old
!=======================================================================================


!=======================================================================================
!> @brief distribute Smolyak terms to threads
!> Only the real part of the packed basis are converged to Smolyak rep. 
!> 1. transfer from packed basis to basis in Smolyak rep.
!> 2. transfer from basis Re. to grid rep. 
!> 3. grid Smolyak rep. are reserved for the current time step
!
!> if psi%clpx, the final SR_B will be (:,2) dim matrix, with real and imag. part resp.
!> otherwise, SR_B will be (:,1)
!=======================================================================================
  SUBROUTINE SmolyakR_distribute_SRB_MPI(psi,para_Op)
    USE mod_basis_set_alloc
    USE mod_param_SGType2
    USE mod_nDindex
    USE mod_psi_set_alloc,          ONLY:param_psi
    USE mod_Op,                     ONLY:param_Op
    USE mod_basis_BtoG_GtoB_SGType4,ONLY:tabPackedBasis_TO_tabR_AT_iG,                 &
                                         BDP_TO_GDP_OF_SmolyakRep,                     &
                                         getbis_tab_nq,getbis_tab_nb
    USE mod_MPI_aux
    IMPLICIT NONE
    
    TYPE(param_psi),                intent(inout) :: psi
    TYPE(param_Op),                 intent(in)    :: para_Op
    
    TYPE(param_SGType2),pointer                   :: SGType2
    Real(kind=Rkind),allocatable                  :: SR_B0(:,:)
    Integer,pointer                               :: tab_l(:,:)
    Integer,allocatable                           :: tab_nq(:)
    Integer,allocatable                           :: tab_nb(:)
    Integer                                       :: iG
    Integer                                       :: nb0
    Integer                                       :: dim
    
    Character(len=*),parameter             :: name_sub='SmolyakR_distribute_SR_MPI'

#if(run_MPI)

    psi%SRB_MPI=.TRUE.

    SGType2 => para_Op%BasisnD%para_SGType2
    tab_l   => SGType2%nDind_SmolyakRep%Tab_nDval(:,:)
    nb0=SGType2%nb0

    !> allocate grid SR, the relevant index record the vector position for each iG
    !> psi%SR_B(psi%SR_B_index(iG):psi%SR_B_index(iG+1)-1) => SRB(iG)
    CALL allocate_array(psi%SR_B_length,0,MPI_np-1)
    CALL allocate_array(psi%SR_B_index,iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)+1)
    psi%SR_B_index(iGs_MPI(1,MPI_id))=1

    DO i_MPI=0,MPI_np-1
      psi%SR_B_length(i_MPI)=0
      DO iG=iGs_MPI(1,i_MPI),iGs_MPI(2,i_MPI)
        psi%SR_B_length(i_MPI)=psi%SR_B_length(i_MPI)+SGType2%tab_nb_OF_SRep(iG)*nb0
        IF(i_MPI==MPI_id) psi%SR_B_index(iG+1)=psi%SR_B_length(MPI_id)+1 ! for each threads
      ENDDO
    ENDDO

!    DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)+1
!      write(*,*) 'iGs_MPI check: iG=',iG,psi%SR_B_index(iG),iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)+1
!    ENDDO

    !> allocate SR_B for each threads
    IF(psi%cplx) THEN
      CALL allocate_array(psi%SR_B,1,psi%SR_B_length(MPI_id),1,2)
    ELSE 
      CALL allocate_array(psi%SR_B,1,psi%SR_B_length(MPI_id),1,1)
    ENDIF

    !> distribute SR on Basis to each threads
    IF(MPI_id==0) THEN
      DO i_MPI=1,MPI_np-1
      
        !> allocate temprary array for sending
        IF(psi%cplx) THEN
          CALL allocate_array(SR_B0,1,psi%SR_B_length(i_MPI),1,2)
        ELSE
          CALL allocate_array(SR_B0,1,psi%SR_B_length(i_MPI),1,1)
        ENDIF
        
        !> loop to send SR_B
        DO iG=iGs_MPI(1,i_MPI),iGs_MPI(2,i_MPI)
          IF(psi%cplx) THEN
            ! Real part
            CALL pack_Basis_to_SR_Basis_MPI(iG,Real(psi%CvecB,kind=Rkind),SR_B0(:,1),  &
                                            psi%SR_B_index,para_Op)
            ! Imaginary part
            CALL pack_Basis_to_SR_Basis_MPI(iG,aimag(psi%CvecB),SR_B0(:,2),            &
                                            psi%SR_B_index,para_Op)
          ELSE
            CALL pack_Basis_to_SR_Basis_MPI(iG,psi%RvecB,SR_B0(:,1),                   &
                                            psi%SR_B_index,para_Op)
          ENDIF
        ENDDO
      
        ! send SR_B
        IF(psi%cplx) THEN
          CALL MPI_Send(SR_B0(1:psi%SR_B_length(i_MPI),1),psi%SR_B_length(i_MPI),      &
                              Real_MPI,root_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
          CALL MPI_Send(SR_B0(1:psi%SR_B_length(i_MPI),2),psi%SR_B_length(i_MPI),      &
                              Real_MPI,root_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
        ELSE
          CALL MPI_Send(SR_B0(1:psi%SR_B_length(i_MPI),1),psi%SR_B_length(i_MPI),      &
                              Real_MPI,root_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
        ENDIF
      ENDDO
      
      ! Smolyak rep. on MPI_id=0
      DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)
        IF(psi%cplx) THEN
          ! Real part
          CALL pack_Basis_to_SR_Basis_MPI(iG,Real(psi%CvecB,kind=Rkind),psi%SR_B(:,1), &
                                          psi%SR_B_index,para_Op)
          ! Imaginary part
          CALL pack_Basis_to_SR_Basis_MPI(iG,aimag(psi%CvecB),psi%SR_B(:,2),           &
                                          psi%SR_B_index,para_Op)
        ELSE
          CALL pack_Basis_to_SR_Basis_MPI(iG,psi%RvecB,psi%SR_B(:,1),                  &
                                          psi%SR_B_index,para_Op)
        ENDIF
      ENDDO
    ENDIF ! for MPI_id==0
    
    !> receive SR_B at each threads
    IF(MPI_id/=0) THEN
      IF(psi%cplx) THEN
        CALL MPI_Recv(psi%SR_B(1:psi%SR_B_length(MPI_id),1),psi%SR_B_length(MPI_id),   &
                       Real_MPI,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
        CALL MPI_Recv(psi%SR_B(1:psi%SR_B_length(MPI_id),2),psi%SR_B_length(MPI_id),   &
                       Real_MPI,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
      ELSE
        CALL MPI_Recv(psi%SR_B(1:psi%SR_B_length(MPI_id),1),psi%SR_B_length(MPI_id),   &
                       Real_MPI,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
      ENDIF
    ENDIF

    ! free space
    IF(allocated(SR_B0)) deallocate(SR_B0)

#endif
  END SUBROUTINE SmolyakR_distribute_SRB_MPI
!=======================================================================================


!=======================================================================================
!> @brief distribute Smolyak terms to threads
!> Only the real part of the packed basis are converged to Smolyak rep. 
!> 1. transfer from packed basis to basis in Smolyak rep.
!> 2. transfer from basis Re. to grid rep. 
!> 3. grid Smolyak rep. are reserved for the current time step
!
!> if psi%clpx, the final SR_G will be (:,2) dim matrix, with real and imag. part resp.
!> otherwise, SR_G will be (:,1)
!=======================================================================================
  SUBROUTINE SmolyakR_distribute_SRG_MPI(psi,para_Op)
    USE mod_basis_set_alloc
    USE mod_param_SGType2
    USE mod_nDindex
    USE mod_psi_set_alloc,          ONLY:param_psi
    USE mod_Op,                     ONLY:param_Op
    USE mod_basis_BtoG_GtoB_SGType4,ONLY:tabPackedBasis_TO_tabR_AT_iG,                 &
                                         BDP_TO_GDP_OF_SmolyakRep,                     &
                                         getbis_tab_nq,getbis_tab_nb
    USE mod_MPI_aux
    IMPLICIT NONE
    
    TYPE(param_psi),                intent(inout) :: psi
    TYPE(param_Op),                 intent(in)    :: para_Op
    
    TYPE(param_SGType2),pointer                   :: SGType2
    Real(kind=Rkind),allocatable                  :: SR_G0(:,:)
    Integer,pointer                               :: tab_l(:,:)
    Integer,allocatable                           :: tab_nq(:)
    Integer,allocatable                           :: tab_nb(:)
    Integer                                       :: iG
    Integer                                       :: nb0
    Integer                                       :: dim

    Character(len=*),parameter             :: name_sub='SmolyakR_distribute_SR_MPI'

#if(run_MPI)

    psi%SRG_MPI=.TRUE.

    SGType2 => para_Op%BasisnD%para_SGType2
    tab_l   => SGType2%nDind_SmolyakRep%Tab_nDval(:,:)
    
    dim=size(tab_l,1)
    CALL alloc_NParray(tab_nb,(/dim/),'tab_nb',name_sub)
    CALL alloc_NParray(tab_nq,(/dim/),'tab_nq',name_sub)

    !> allocate grid SR, the relevant index record the vector position for each iG
    !> psi%SR_G(psi%SR_G_index(iG):psi%SR_G_index(iG+1)-1) => SRG(iG)
    nb0=SGType2%nb0
    CALL allocate_array(psi%SR_G_length,0,MPI_np-1)
    CALL allocate_array(psi%SR_G_index,iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)+1)
    psi%SR_G_index(iGs_MPI(1,MPI_id))=1
    
    !IF(MPI_id==0) CALL allocate_array(psi%SR_G_index0,1,SGType2%nb_SG)
    
    DO i_MPI=0,MPI_np-1
      psi%SR_G_length(i_MPI)=0
      !IF(MPI_id==0) psi%SR_G_index0(iGs_MPI(1,i_MPI))=1
      DO iG=iGs_MPI(1,i_MPI),iGs_MPI(2,i_MPI)
        tab_nq(:)=getbis_tab_nq(tab_l(:,iG),para_Op%BasisnD%tab_basisPrimSG)
        psi%SR_G_length(i_MPI)=psi%SR_G_length(i_MPI)+product(tab_nq)*nb0
        !IF(MPI_id==0) psi%SR_G_index0(iG+1)=psi%SR_G_length(i_MPI)+1
        IF(i_MPI==MPI_id) psi%SR_G_index(iG+1)=psi%SR_G_length(MPI_id)+1 ! for each threads
      ENDDO
    ENDDO
    
!    DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)+1
!      write(*,*) 'iGs_MPI check: iG=',iG,psi%SR_G_index(iG),iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)+1
!    ENDDO
    
    !> allocate SR_G for each threads
    IF(psi%cplx) THEN
      CALL allocate_array(psi%SR_G,1,psi%SR_G_length(MPI_id),1,2)
    ELSE 
      CALL allocate_array(psi%SR_G,1,psi%SR_G_length(MPI_id),1,1)
    ENDIF

    !> distribute SR on grid to each threads
    IF(MPI_id==0) THEN
      DO i_MPI=1,MPI_np-1
      
        !> allocate temprary array for sending
        IF(psi%cplx) THEN
          CALL allocate_array(SR_G0,1,psi%SR_G_length(i_MPI),1,2)
        ELSE
          CALL allocate_array(SR_G0,1,psi%SR_G_length(i_MPI),1,1)
        ENDIF
        
        !> loop to send SR_G
        DO iG=iGs_MPI(1,i_MPI),iGs_MPI(2,i_MPI)
          tab_nq(:)=getbis_tab_nq(tab_l(:,iG),para_Op%BasisnD%tab_basisPrimSG)
          tab_nb(:)=getbis_tab_nb(tab_l(:,iG),para_Op%BasisnD%tab_basisPrimSG)

!          ! get SR on Basis
!          IF(psi%cplx) THEN
!            !-Real part-----------------------------------------------------------------
!            CALL tabPackedBasis_TO_tabR_AT_iG(psi%SR_B,Real(psi%CvecB,kind=Rkind),iG,  &
!                                              SGType2%tab_nb_OF_SRep(iG)*nb0)
!            ! transfer to SR on grid
!            CALL BDP_TO_GDP_OF_SmolyakRep(psi%SR_B,para_Op%BasisnD%tab_basisPrimSG,    &
!                                          tab_l(:,iG),tab_nq,tab_nb,nb0)
!            !< psi%SR_B is now actually psi%SR_G at iG
!            IF(size(psi%SR_B)-(psi%SR_G_index(iG+1)-psi%SR_G_index(iG))/=0)            &
!                              STOP 'error in SmolyakR_distribute_MPI, check psi%SR_B size'
!            SR_G0(psi%SR_G_index(iG):psi%SR_G_index(iG+1)-1,1)=psi%SR_B
!            
!            !-Imaginary part------------------------------------------------------------
!            CALL tabPackedBasis_TO_tabR_AT_iG(psi%SR_B,aimag(psi%CvecB),iG,            &
!                                              SGType2%tab_nb_OF_SRep(iG)*nb0)
!            ! transfer to SR on grid
!            CALL BDP_TO_GDP_OF_SmolyakRep(psi%SR_B,para_Op%BasisnD%tab_basisPrimSG,    &
!                                          tab_l(:,iG),tab_nq,tab_nb,nb0)
!            SR_G0(psi%SR_G_index(iG):psi%SR_G_index(iG+1)-1,2)=psi%SR_B
!          ELSE
!            CALL tabPackedBasis_TO_tabR_AT_iG(psi%SR_B,psi%RvecB,iG,                   &
!                                              SGType2%tab_nb_OF_SRep(iG)*nb0)
!            ! transfer to SR on grid
!            CALL BDP_TO_GDP_OF_SmolyakRep(psi%SR_B,para_Op%BasisnD%tab_basisPrimSG,    &
!                                        tab_l(:,iG),tab_nq,tab_nb,nb0)
!            SR_G0(psi%SR_G_index(iG):psi%SR_G_index(iG+1)-1,1)=psi%SR_B
!          ENDIF
          IF(psi%cplx) THEN
            ! Real part
            CALL pack_Basis_to_SR_grid_MPI(iG,Real(psi%CvecB,kind=Rkind),SR_G0(:,1),   &
                                       psi%SR_G_index,tab_nq,tab_nb,tab_l(:,iG),para_Op)
            ! Imaginary part
            CALL pack_Basis_to_SR_grid_MPI(iG,aimag(psi%CvecB),SR_G0(:,2),             &
                                       psi%SR_G_index,tab_nq,tab_nb,tab_l(:,iG),para_Op)
          ELSE
            CALL pack_Basis_to_SR_grid_MPI(iG,psi%RvecB,SR_G0(:,1),                    &
                                       psi%SR_G_index,tab_nq,tab_nb,tab_l(:,iG),para_Op)
          ENDIF
        ENDDO
      
        ! send SR_G
        IF(psi%cplx) THEN
          CALL MPI_Send(SR_G0(1:psi%SR_G_length(i_MPI),1),psi%SR_G_length(i_MPI),      &
                              Real_MPI,root_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
          CALL MPI_Send(SR_G0(1:psi%SR_G_length(i_MPI),2),psi%SR_G_length(i_MPI),      &
                              Real_MPI,root_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
        ELSE
          CALL MPI_Send(SR_G0(1:psi%SR_G_length(i_MPI),1),psi%SR_G_length(i_MPI),      &
                              Real_MPI,root_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
        ENDIF
      ENDDO
      
      ! Smolyak rep. on MPI_id=0
      DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)
        tab_nq(:)=getbis_tab_nq(tab_l(:,iG),para_Op%BasisnD%tab_basisPrimSG)
        tab_nb(:)=getbis_tab_nb(tab_l(:,iG),para_Op%BasisnD%tab_basisPrimSG)
        
        IF(psi%cplx) THEN
          ! Real part
          CALL pack_Basis_to_SR_grid_MPI(iG,Real(psi%CvecB,kind=Rkind),psi%SR_G(:,1),  &
                                       psi%SR_G_index,tab_nq,tab_nb,tab_l(:,iG),para_Op)
          ! Imaginary part
          CALL pack_Basis_to_SR_grid_MPI(iG,aimag(psi%CvecB),psi%SR_G(:,2),            &
                                       psi%SR_G_index,tab_nq,tab_nb,tab_l(:,iG),para_Op)
        ELSE
          CALL pack_Basis_to_SR_grid_MPI(iG,psi%RvecB,psi%SR_G(:,1),                   &
                                       psi%SR_G_index,tab_nq,tab_nb,tab_l(:,iG),para_Op)
        ENDIF
      ENDDO
    ENDIF ! for MPI_id==0
    
    !> receive SR_G at each threads
    IF(MPI_id/=0) THEN
      IF(psi%cplx) THEN
        CALL MPI_Recv(psi%SR_G(1:psi%SR_G_length(MPI_id),1),psi%SR_G_length(MPI_id),   &
                       Real_MPI,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
        CALL MPI_Recv(psi%SR_G(1:psi%SR_G_length(MPI_id),2),psi%SR_G_length(MPI_id),   &
                       Real_MPI,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
      ELSE
        CALL MPI_Recv(psi%SR_G(1:psi%SR_G_length(MPI_id),1),psi%SR_G_length(MPI_id),   &
                       Real_MPI,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
      ENDIF
    ENDIF

    ! free space
    IF(allocated(psi%SR_B)) deallocate(psi%SR_B)
    IF(allocated(SR_G0)) deallocate(SR_G0)

#endif
  END SUBROUTINE SmolyakR_distribute_SRG_MPI
!=======================================================================================


!=======================================================================================
!> transfer Smolyak rep. back to packed Basisi rep.
!=======================================================================================
  SUBROUTINE SmolyakR_to_packedB_SRB_MPI(psi,para_Op)
    USE mod_basis_set_alloc
    USE mod_param_SGType2
    USE mod_nDindex
    USE mod_SymAbelian,             ONLY:Calc_symab1_EOR_symab2
    USE mod_psi_set_alloc,          ONLY:param_psi
    USE mod_psi_Op,                 ONLY:Set_symab_OF_psiBasisRep
    USE mod_Op,                     ONLY:param_Op
    USE mod_basis_BtoG_GtoB_SGType4,ONLY:tabR_AT_iG_TO_tabPackedBasis,                 &
                                         GDP_TO_BDP_OF_SmolyakRep
    USE mod_MPI_aux
    IMPLICIT NONE
    
    TYPE(param_psi),                intent(inout) :: psi
    TYPE(param_Op),                 intent(in)    :: para_Op
    
    TYPE(param_SGType2),pointer                   :: SGType2
    Real(kind=Rkind),allocatable                  :: SR_B0(:)
    Real(kind=Rkind),allocatable                  :: Rvec(:,:)
    Integer                                       :: nb0
    Integer                                       :: psi_symab
    Integer                                       :: iG

    Character(len=*),parameter             :: name_sub='SmolyakR_to_packedB_SRB_MPI'

#if(run_MPI)

    SGType2 => para_Op%BasisnD%para_SGType2
    nb0=SGType2%nb0

    ! send 
    IF(MPI_id/=0) THEN
      IF(psi%cplx) THEN
        CALL MPI_Send(psi%SR_B(1:psi%SR_B_length(MPI_id),1),psi%SR_B_length(MPI_id),   &
                      Real_MPI,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_err)
        CALL MPI_Send(psi%SR_B(1:psi%SR_B_length(MPI_id),2),psi%SR_B_length(MPI_id),   &
                      Real_MPI,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_err)
      ELSE
        CALL MPI_Send(psi%SR_B(1:psi%SR_B_length(MPI_id),1),psi%SR_B_length(MPI_id),   &
                      Real_MPI,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_err)
      ENDIF
    ENDIF

    IF(MPI_id==0) THEN
      ! MPI_id=0 prat
      IF(psi%cplx) THEN
        CALL allocate_array(Rvec,1,size(psi%CvecB),1,2)
      ENDIF

      DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)
        IF(psi%cplx) THEN
          CALL SR_Basis_to_pack_Basis_MPI(iG,Rvec(:,1),psi%SR_B(:,1),psi%SR_B_index,   &
                                          para_Op)
          CALL SR_Basis_to_pack_Basis_MPI(iG,Rvec(:,2),psi%SR_B(:,2),psi%SR_B_index,   &
                                          para_Op)
        ELSE
          CALL SR_Basis_to_pack_Basis_MPI(iG,psi%RvecB,psi%SR_B(:,1),psi%SR_B_index,   &
                                          para_Op)
        ENDIF
      ENDDO
      IF(psi%cplx) psi%CvecB=cmplx(Rvec(:,1),Rvec(:,2),kind=Rkind)

      ! receive SRG from the other threads
      DO i_MPI=1,MPI_np-1
        CALL allocate_array(SR_B0,1,psi%SR_B_length(i_MPI))
        IF(psi%cplx) THEN
          CALL MPI_Recv(SR_B0(1:psi%SR_B_length(i_MPI)),psi%SR_B_length(i_MPI),      &
                        Real_MPI,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_stat,MPI_err)
          DO iG=iGs_MPI(1,i_MPI),iGs_MPI(2,i_MPI)
            CALL SR_Basis_to_pack_Basis_MPI(iG,Rvec(:,1),SR_B0,psi%SR_B_index,          &
                                            para_Op)
          ENDDO
          
          CALL MPI_Recv(SR_B0(1:psi%SR_B_length(i_MPI)),psi%SR_B_length(i_MPI),      &
                        Real_MPI,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_stat,MPI_err)
          DO iG=iGs_MPI(1,i_MPI),iGs_MPI(2,i_MPI)
            CALL SR_Basis_to_pack_Basis_MPI(iG,Rvec(:,2),SR_B0,psi%SR_B_index,          &
                                            para_Op)
          ENDDO
          psi%CvecB=cmplx(Rvec(:,1),Rvec(:,2),kind=Rkind)
        ELSE
          CALL MPI_Recv(SR_B0(1:psi%SR_B_length(i_MPI)),psi%SR_B_length(i_MPI),      &
                        Real_MPI,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_stat,MPI_err)
          DO iG=iGs_MPI(1,i_MPI),iGs_MPI(2,i_MPI)
            CALL SR_Basis_to_pack_Basis_MPI(iG,psi%RvecB,SR_B0,psi%SR_B_index,para_Op)
          ENDDO
        ENDIF ! psi%cplx
      ENDDO ! for i_MPI=1,MPI_np-1

      psi_symab=Calc_symab1_EOR_symab2(para_Op%symab,psi%symab)
      CALL Set_symab_OF_psiBasisRep(psi,psi_symab)
    ENDIF ! for MPI_id==0

    IF(allocated(Rvec))      deallocate(Rvec)
    IF(allocated(SR_B0))     deallocate(SR_B0)
    IF(allocated(psi%SR_B))  deallocate(psi%SR_B)
    
    psi%SRB_MPI=.FALSE.

#endif
  ENDSUBROUTINE SmolyakR_to_packedB_SRB_MPI
!=======================================================================================


!=======================================================================================
!> transfer Smolyak rep. back to packed Basisi rep.
!=======================================================================================
  SUBROUTINE SmolyakR_to_packedB_SRG_MPI(psi,para_Op)
    USE mod_basis_set_alloc
    USE mod_param_SGType2
    USE mod_nDindex
    USE mod_SymAbelian,             ONLY:Calc_symab1_EOR_symab2
    USE mod_psi_set_alloc,          ONLY:param_psi
    USE mod_psi_Op,                 ONLY:Set_symab_OF_psiBasisRep
    USE mod_Op,                     ONLY:param_Op
    USE mod_basis_BtoG_GtoB_SGType4,ONLY:tabR_AT_iG_TO_tabPackedBasis,                 &
                                         GDP_TO_BDP_OF_SmolyakRep,                     &
                                         getbis_tab_nq,getbis_tab_nb
    USE mod_MPI_aux
    IMPLICIT NONE
    
    TYPE(param_psi),                intent(inout) :: psi
    TYPE(param_Op),                 intent(in)    :: para_Op
    
    TYPE(param_SGType2),pointer                   :: SGType2
    Real(kind=Rkind),allocatable                  :: SR_G0(:)
    Real(kind=Rkind),allocatable                  :: Rvec(:,:)
    Integer,pointer                               :: tab_l(:,:)
    Integer,allocatable                           :: tab_nq(:)
    Integer,allocatable                           :: tab_nb(:)
    Integer                                       :: nb0
    Integer                                       :: psi_symab
    Integer                                       :: iG
    Integer                                       :: dim

    Character(len=*),parameter             :: name_sub='SmolyakR_to_packedB_SR_MPI'

#if(run_MPI)

    SGType2 => para_Op%BasisnD%para_SGType2
    tab_l   => SGType2%nDind_SmolyakRep%Tab_nDval(:,:)
    nb0=SGType2%nb0
    
    dim=size(tab_l,1)
    CALL alloc_NParray(tab_nb,(/dim/),'tab_nb',name_sub)
    CALL alloc_NParray(tab_nq,(/dim/),'tab_nq',name_sub)

    ! send 
    IF(MPI_id/=0) THEN
      IF(psi%cplx) THEN
        CALL MPI_Send(psi%SR_G(1:psi%SR_G_length(MPI_id),1),psi%SR_G_length(MPI_id),   &
                      Real_MPI,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_err)
        CALL MPI_Send(psi%SR_G(1:psi%SR_G_length(MPI_id),2),psi%SR_G_length(MPI_id),   &
                      Real_MPI,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_err)
      ELSE
        CALL MPI_Send(psi%SR_G(1:psi%SR_G_length(MPI_id),1),psi%SR_G_length(MPI_id),   &
                      Real_MPI,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_err)
      ENDIF
    ENDIF

    IF(MPI_id==0) THEN
      ! MPI_id=0 prat
      IF(psi%cplx) THEN
        CALL allocate_array(Rvec,1,size(psi%CvecB),1,2)
      ENDIF

      DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)
        tab_nq(:)=getbis_tab_nq(tab_l(:,iG),para_Op%BasisnD%tab_basisPrimSG)
        tab_nb(:)=getbis_tab_nb(tab_l(:,iG),para_Op%BasisnD%tab_basisPrimSG)

        IF(psi%cplx) THEN
          CALL SR_grid_to_pack_Basis_MPI(iG,Rvec(:,1),psi%SR_G(:,1),psi%SR_G_index,    &
                                         tab_nq,tab_nb,tab_l(:,iG),para_Op)
          CALL SR_grid_to_pack_Basis_MPI(iG,Rvec(:,2),psi%SR_G(:,2),psi%SR_G_index,    &
                                         tab_nq,tab_nb,tab_l(:,iG),para_Op)
        ELSE
          CALL SR_grid_to_pack_Basis_MPI(iG,psi%RvecB,psi%SR_G(:,1),psi%SR_G_index,    &
                                         tab_nq,tab_nb,tab_l(:,iG),para_Op)
        ENDIF

!        IF(psi%cplx) THEN
!          CALL GDP_TO_BDP_OF_SmolyakRep(psi%SR_G(:,1),para_Op%BasisnD%tab_basisPrimSG,&
!                                        tab_l,tab_nq,tab_nb,nb0)
!          CALL GDP_TO_BDP_OF_SmolyakRep(psi%SR_G(:,2),para_Op%BasisnD%tab_basisPrimSG,&
!                                        tab_l,tab_nq,tab_nb,nb0)
!          CALL tabR_AT_iG_TO_tabPackedBasis(Rvec(:,1),psi%SR_G(:,1),iG,SGType2,       &
!                                            para_Op%BasisnD%WeightSG(iG))
!          CALL tabR_AT_iG_TO_tabPackedBasis(Rvec(:,2),psi%SR_G(:,2),iG,SGType2,       &
!                                            para_Op%BasisnD%WeightSG(iG))
!          psi%CvecB=cmplx(Rvec(:,1),Rvec(:,2))
!        ELSE
!          CALL GDP_TO_BDP_OF_SmolyakRep(psi%SR_G(:,1),para_Op%BasisnD%tab_basisPrimSG,&
!                                        tab_l,tab_nq,tab_nb,nb0)
!          CALL tabR_AT_iG_TO_tabPackedBasis(psi%RvecB,psi%SR_G(:,1),iG,SGType2,        &
!                                            para_Op%BasisnD%WeightSG(iG))
!        ENDIF
      ENDDO
      IF(psi%cplx) psi%CvecB=cmplx(Rvec(:,1),Rvec(:,2),kind=Rkind)

      ! receive SRG from the other threads
      DO i_MPI=1,MPI_np-1
        CALL allocate_array(SR_G0,1,psi%SR_G_length(i_MPI))
        IF(psi%cplx) THEN
          CALL MPI_Recv(SR_G0(1:psi%SR_G_length(i_MPI)),psi%SR_G_length(i_MPI),      &
                        Real_MPI,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_stat,MPI_err)
          DO iG=iGs_MPI(1,i_MPI),iGs_MPI(2,i_MPI)
            CALL SR_grid_to_pack_Basis_MPI(iG,Rvec(:,1),SR_G0,psi%SR_G_index,          &
                                           tab_nq,tab_nb,tab_l(:,iG),para_Op)
          ENDDO
          
          CALL MPI_Recv(SR_G0(1:psi%SR_G_length(i_MPI)),psi%SR_G_length(i_MPI),      &
                        Real_MPI,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_stat,MPI_err)
          DO iG=iGs_MPI(1,i_MPI),iGs_MPI(2,i_MPI)
            CALL SR_grid_to_pack_Basis_MPI(iG,Rvec(:,2),SR_G0,psi%SR_G_index,          &
                                           tab_nq,tab_nb,tab_l(:,iG),para_Op)
          ENDDO
          psi%CvecB=cmplx(Rvec(:,1),Rvec(:,2),kind=Rkind)
        ELSE
          CALL MPI_Recv(SR_G0(1:psi%SR_G_length(i_MPI)),psi%SR_G_length(i_MPI),      &
                        Real_MPI,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_stat,MPI_err)
          DO iG=iGs_MPI(1,i_MPI),iGs_MPI(2,i_MPI)
            CALL SR_grid_to_pack_Basis_MPI(iG,psi%RvecB,SR_G0,psi%SR_G_index,          &
                                           tab_nq,tab_nb,tab_l(:,iG),para_Op)
          ENDDO
        ENDIF ! psi%cplx
      ENDDO ! for i_MPI=1,MPI_np-1

      psi_symab=Calc_symab1_EOR_symab2(para_Op%symab,psi%symab)
      CALL Set_symab_OF_psiBasisRep(psi,psi_symab)
    ENDIF ! for MPI_id==0

    IF(allocated(Rvec)) deallocate(Rvec)
    IF(allocated(SR_G0)) deallocate(SR_G0)
    IF(allocated(tab_nq)) deallocate(tab_nq)
    IF(allocated(tab_nb)) deallocate(tab_nb)

#endif
  ENDSUBROUTINE SmolyakR_to_packedB_SRG_MPI
!=======================================================================================


!=======================================================================================
!> @brief get Smolyak rep on Basis from packed basis
!=======================================================================================
  SUBROUTINE pack_Basis_to_SR_Basis_MPI(iG,packB,SR_B,SR_B_index,para_Op)
    USE mod_basis_set_alloc
    USE mod_basis_BtoG_GtoB_SGType4,ONLY:tabPackedBasis_TO_tabR_AT_iG
    USE mod_Op,                     ONLY:param_Op
    USE mod_MPI_aux
    IMPLICIT NONE

    TYPE(param_Op),                 intent(in)    :: para_Op
    Real(kind=Rkind),               intent(in)    :: packB(:)
    Real(kind=Rkind),               intent(inout) :: SR_B(:)
    Integer,                        intent(in)    :: SR_B_index(:)
    Integer,                        intent(in)    :: iG

    Real(kind=Rkind),allocatable                  :: SRB(:)

#if(run_MPI)

    ! transfer to SR on Basis
    CALL tabPackedBasis_TO_tabR_AT_iG(SRB,packB,iG,para_Op%BasisnD%para_SGType2)

    !< store in psi%SR_B
    IF(size(SRB)-(SR_B_index(iG+1)-SR_B_index(iG))/=0)                                &
                         STOP 'error in SmolyakR_distribute_SR_MPI, check psi%SR_B size'
    SR_B(SR_B_index(iG):SR_B_index(iG+1)-1)=SRB

#endif
  ENDSUBROUTINE pack_Basis_to_SR_Basis_MPI
!=======================================================================================


!=======================================================================================
!> @brief get Smolyak rep on grid from packed basis
!=======================================================================================
  SUBROUTINE pack_Basis_to_SR_grid_MPI(iG,packB,SR_G,SR_G_index,                       &
                                       tab_nq,tab_nb,tab_l,para_Op)
    USE mod_basis_set_alloc
    USE mod_param_SGType2
    USE mod_basis_BtoG_GtoB_SGType4,ONLY:tabPackedBasis_TO_tabR_AT_iG,                 &
                                         BDP_TO_GDP_OF_SmolyakRep
    USE mod_Op,                     ONLY:param_Op
    USE mod_MPI_aux
    IMPLICIT NONE

    TYPE(param_Op),                 intent(in)    :: para_Op
    Real(kind=Rkind),               intent(inout) :: SR_G(:)
    Real(kind=Rkind),               intent(in)    :: packB(:)
    Integer,                        intent(in)    :: SR_G_index(:)
    Integer,allocatable,            intent(in)    :: tab_nq(:)
    Integer,allocatable,            intent(in)    :: tab_nb(:)
    Integer,                        intent(in)    :: tab_l(:)
    Integer,                        intent(in)    :: iG

    Real(kind=Rkind),allocatable                  :: SR_B(:)
    TYPE(param_SGType2),pointer                   :: SGType2

#if(run_MPI)

    SGType2 => para_Op%BasisnD%para_SGType2
  
    CALL tabPackedBasis_TO_tabR_AT_iG(SR_B,packB,iG,SGType2)
    ! transfer to SR on grid
    CALL BDP_TO_GDP_OF_SmolyakRep(SR_B,para_Op%BasisnD%tab_basisPrimSG,                &
                                  tab_l,tab_nq,tab_nb,SGType2%nb0)
    !< psi%SR_B is now actually psi%SR_G at iG
    IF(size(SR_B)-(SR_G_index(iG+1)-SR_G_index(iG))/=0)                                &
                         STOP 'error in SmolyakR_distribute_SR_MPI, check psi%SR_B size'
    SR_G(SR_G_index(iG):SR_G_index(iG+1)-1)=SR_B

#endif
  ENDSUBROUTINE pack_Basis_to_SR_grid_MPI
!=======================================================================================


!=======================================================================================
!> @brief transfer Smolyak rep on grid to packed basis
!=======================================================================================
  SUBROUTINE SR_Basis_to_pack_Basis_MPI(iG,packB,SR_B,SR_B_index,para_Op)
    USE mod_basis_set_alloc
    USE mod_param_SGType2
    USE mod_basis_BtoG_GtoB_SGType4,ONLY:tabR_AT_iG_TO_tabPackedBasis
    USE mod_Op,                     ONLY:param_Op
    USE mod_MPI_aux
    IMPLICIT NONE

    TYPE(param_Op),                 intent(in)    :: para_Op
    Real(kind=Rkind),               intent(in)    :: SR_B(:)
    Real(kind=Rkind),               intent(inout) :: packB(:)
    Integer,                        intent(in)    :: SR_B_index(:)
    Integer,                        intent(in)    :: iG
    
    Real(kind=Rkind),allocatable                  :: SRB(:)
    Integer                                       :: nb0

#if(run_MPI)

    nb0=para_Op%BasisnD%para_SGType2%nb0

    ! the SR_B at iG is stored in SRB
    CALL allocate_array(SRB,1,SR_B_index(iG+1)-SR_B_index(iG))
    SRB=SR_B(SR_B_index(iG):SR_B_index(iG+1)-1)

    ! now SR_B stores the basis rep. 
    CALL tabR_AT_iG_TO_tabPackedBasis(packB,SRB,iG,para_Op%BasisnD%para_SGType2,       &
                                      para_Op%BasisnD%WeightSG(iG))
    IF(allocated(SRB)) deallocate(SRB)

#endif
  ENDSUBROUTINE SR_Basis_to_pack_Basis_MPI
!=======================================================================================


!=======================================================================================
!> @brief transfer Smolyak rep on grid to packed basis
!=======================================================================================
  SUBROUTINE SR_grid_to_pack_Basis_MPI(iG,packB,SR_G,SR_G_index,                       &
                                       tab_nq,tab_nb,tab_l,para_Op)
    USE mod_basis_set_alloc
    USE mod_param_SGType2
    USE mod_basis_BtoG_GtoB_SGType4,ONLY:tabR_AT_iG_TO_tabPackedBasis,                 &
                                         GDP_TO_BDP_OF_SmolyakRep
    USE mod_Op,                     ONLY:param_Op
    USE mod_MPI_aux
    IMPLICIT NONE

    TYPE(param_Op),                 intent(in)    :: para_Op
    Real(kind=Rkind),               intent(in)    :: SR_G(:)
    Real(kind=Rkind),               intent(inout) :: packB(:)
    Integer,                        intent(in)    :: SR_G_index(:)
    Integer,allocatable,            intent(in)    :: tab_nq(:)
    Integer,allocatable,            intent(in)    :: tab_nb(:)
    Integer,                        intent(in)    :: tab_l(:)
    Integer,                        intent(in)    :: iG
    
    Real(kind=Rkind),allocatable                  :: SR_temp(:)
    Integer                                       :: nb0

#if(run_MPI)

    nb0=para_Op%BasisnD%para_SGType2%nb0
    
    ! the SR_G at iG is stored in SR_B
    CALL allocate_array(SR_temp,1,SR_G_index(iG+1)-SR_G_index(iG))
    SR_temp=SR_G(SR_G_index(iG):SR_G_index(iG+1)-1)
    
    CALL GDP_TO_BDP_OF_SmolyakRep(SR_temp,para_Op%BasisnD%tab_basisPrimSG,             &
                                  tab_l,tab_nq,tab_nb,nb0)
    ! now SR_B stores the basis rep. 
    CALL tabR_AT_iG_TO_tabPackedBasis(packB,SR_temp,iG,para_Op%BasisnD%para_SGType2,   &
                                      para_Op%BasisnD%WeightSG(iG))
    IF(allocated(SR_temp)) deallocate(SR_temp)

#endif
  ENDSUBROUTINE
!=======================================================================================


!=======================================================================================
!> calculate psi1=psi1-C*psi2
!> e.g. w1=w1-H(k,k-1)*tab_KrylovSpace(k-1) !> H|v_{k}>-beta_{k-1}|v_{k-1}>
!>      w1=w1-Overlap*tab_KrylovSpace(k)
!=======================================================================================
   SUBROUTINE w1_minus_Cv_SR_MPI(psi1,Cplx,psi2)
     USE mod_psi_set_alloc,          ONLY:param_psi
     IMPLICIT NONE
     
     TYPE(param_psi),                intent(inout) :: psi1
     TYPE(param_psi),                intent(in)    :: psi2
     Complex(kind=Rkind),            intent(in)    :: Cplx
     
     Real(kind=Rkind)                              :: R
     Real(kind=Rkind)                              :: C

#if(run_MPI)

     R=Real(Cplx,kind=Rkind)
     C=aimag(Cplx)

     IF(psi1%SRG_MPI) THEN
       psi1%SR_G(:,1)=psi1%SR_G(:,1)-(psi2%SR_G(:,1)*R-psi2%SR_G(:,2)*C)
       psi1%SR_G(:,2)=psi1%SR_G(:,2)-(psi2%SR_G(:,1)*C+psi2%SR_G(:,2)*R)
     ELSE IF(psi1%SRB_MPI) THEN
       psi1%SR_B(:,1)=psi1%SR_B(:,1)-(psi2%SR_B(:,1)*R-psi2%SR_B(:,2)*C)
       psi1%SR_B(:,2)=psi1%SR_B(:,2)-(psi2%SR_B(:,1)*C+psi2%SR_B(:,2)*R)
     ENDIF

     IF(psi1%symab/=psi2%symab) THEN
       IF(psi1%symab==-2) THEN
         psi1%symab=psi2%symab
       ELSE
         psi1%symab=-1
       ENDIF
     ENDIF

#endif
   ENDSUBROUTINE w1_minus_Cv_SR_MPI
!=======================================================================================

!=======================================================================================
  SUBROUTINE UPsi_spec_MPI(UPsiOnKrylov,H,Vec,Eig,deltaT,n,With_diago)
    IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
    integer,                           intent(in)     :: n
    complex (kind=Rkind),              intent(in)     :: H(n,n)
    complex (kind=Rkind), allocatable, intent(inout)  :: Vec(:,:)
    real (kind=Rkind),    allocatable, intent(inout)  :: Eig(:)
    complex (kind=Rkind), allocatable, intent(inout)  :: UPsiOnKrylov(:)
    real (kind=Rkind),                 intent(in)     :: deltaT
    logical,                           intent(in)     :: With_diago

!------ working variables ---------------------------------
    complex (kind=Rkind) :: coef_i
    integer              :: i

!----- for debuging --------------------------------------------------
    integer, parameter :: nmax = 12
    logical, parameter :: debug=.FALSE.
    !logical, parameter :: debug=.TRUE.
    character (len=*), parameter :: name_sub='UPsi_spec'
!-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) 'deltaT',deltaT
      write(out_unitp,*) 'n',n
      write(out_unitp,*) 'With_diago',With_diago
      IF (With_diago .AND. n <= nmax) THEN
        write(out_unitp,*) 'H'
        CALL Write_Mat(H,out_unitp,6)
      END IF
    END IF
!-----------------------------------------------------------

    IF (With_diago) THEN
      IF (allocated(Vec)) CALL dealloc_NParray(Vec,'Vec',name_sub)
      CALL alloc_NParray(Vec,[n,n],'Vec',name_sub)

      IF (allocated(Eig)) CALL dealloc_NParray(Eig,'Eig',name_sub)
      CALL alloc_NParray(Eig,[n],  'Eig',name_sub)

      IF (allocated(UPsiOnKrylov))                                    &
            CALL dealloc_NParray(UPsiOnKrylov,'UPsiOnKrylov',name_sub)
      CALL alloc_NParray(UPsiOnKrylov,[n],'UPsiOnKrylov',name_sub)

      CALL diagonalization_HerCplx(H,Eig,Vec,n,3,1,.TRUE.)
    END IF

    ! loop on the eigenvectors
    UPsiOnKrylov = CZERO
    DO i=1,n
      coef_i = Vec(1,i) ! just (1,i) because psi on the Krylov space is [1,0,0,...0]
      coef_i = coef_i * exp(-EYE*Eig(i)*deltaT) ! spectral propa

      UPsiOnKrylov(:) = UPsiOnKrylov(:) + conjg(Vec(:,i))*coef_i ! update U.psi on the Krylov space

    END DO

  END SUBROUTINE UPsi_spec_MPI
!=======================================================================================

END MODULE mod_march_MPI
