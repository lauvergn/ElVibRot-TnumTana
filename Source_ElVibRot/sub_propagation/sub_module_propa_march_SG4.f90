!===========================================================================
!===========================================================================
!This file is part of ElVibRot.
!
!    ElVibRot is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ElVibRot is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ElVibRot.  If not, see <http://www.gnu.org/licenses/>.
!
!    Copyright 2015  David Lauvergnat
!      with contributions of Mamadou Ndong, Josep Maria Luis
!
!    ElVibRot includes:
!        - Tnum-Tana under the GNU LGPL3 license
!        - Somme subroutines of John Burkardt under GNU LGPL license
!             http://people.sc.fsu.edu/~jburkardt/
!        - Somme subroutines of SHTOOLS written by Mark A. Wieczorek under BSD license
!             http://shtools.ipgp.fr
!===========================================================================
!===========================================================================

 MODULE mod_march_SG4
 USE mod_system
 USE mod_file
 USE mod_psi_set_alloc, ONLY : param_psi
 USE mod_propa,         ONLY : param_propa,cof,Calc_AutoCorr,Write_AutoCorr
 IMPLICIT NONE

 PRIVATE
 PUBLIC :: march_noD_SG4_BasisRep,march_noD_SG4_GridRep

 CONTAINS

 SUBROUTINE march_noD_SG4_BasisRep(T,no,psi,psi0,para_H,para_propa)
 USE mod_system
!$ USE omp_lib, only : OMP_GET_THREAD_NUM
 USE mod_nDindex
 USE mod_Coord_KEO,                ONLY : zmatrix
 USE mod_basis_set_alloc,          ONLY : basis
 USE mod_basis_RCVec_SGType4
 USE mod_basis_BtoG_GtoB_SGType4,  ONLY : tabPackedBasis_TO_tabR_AT_iG, &
                  tabR_AT_iG_TO_tabPackedBasis,TypeRVec,alloc_TypeRVec,dealloc_TypeRVec, &
                  BDP_TO_GDP_OF_SmolyakRep,GDP_TO_BDP_OF_SmolyakRep, &
                  getbis_tab_nq,getbis_tab_nb

 USE mod_Op,              ONLY : param_Op,write_param_Op
 USE mod_OpPsi,           ONLY : sub_TabOpPsi_OF_ONEDP_FOR_SGtype4
 USE mod_psi_set_alloc,   ONLY : param_psi,copy_psi2TOpsi1,alloc_psi,dealloc_psi,ecri_psi
 USE mod_psi_Op,          ONLY : norme_psi,Overlap_psi1_psi2
 USE mod_psi_SimpleOp,    ONLY : operator (*),operator (+),operator (-),assignment (=)
 IMPLICIT NONE


!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_psi)   :: psi0,psi

!------ working variables ---------------------------------
      complex (kind=Rkind) :: E,cdot
      complex (kind=Rkind) :: psi0Hkpsi0(0:para_propa%para_poly%npoly)
      real (kind=Rkind)    :: microT,T,microdeltaT,phase,microphase
      integer              :: no

      complex (kind=Rkind) :: rt,rtj,rt_tmp
      integer              :: it,j,i_qaie_corr,j_exit
      real (kind=Rkind)    :: E0,rg,norm2_w2


 ! local variables
 TYPE (zmatrix), pointer :: mole
 TYPE (basis),   pointer :: BasisnD


 TYPE (TypeCVec)    :: PsiCVec,Psi0Cvec,Cw1,Cw2(1)
 TYPE (param_psi)   :: MarchPsi

 integer :: ib,i,iG,nb_thread,itab,nq,nb,D,ith,err_sub

  integer,            allocatable       :: tab_nq(:)
  integer,            allocatable       :: tab_nb(:)
  integer,            allocatable       :: tab_l(:)

!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='march_noD_SG4_BasisRep'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'deltaT',para_propa%WPdeltaT
        write(out_unitp,*) 'E0,Esc',para_H%E0,para_H%Esc
        write(out_unitp,*) 'nOD',para_propa%para_poly%npoly
      END IF
!-----------------------------------------------------------
STOP 'march_noD_SG4_BasisRep'
!
!  mole    => para_H%mole
!  BasisnD => para_H%BasisnD
!
!  D = BasisnD%nb_basis
!
!
!  IF (para_H%nb_bie /= 1) STOP 'nb_bie /= 1'
!
!  CALL alloc_NParray(tab_l ,(/ D /),'tab_l', name_sub)
!  CALL alloc_NParray(tab_nb,(/ D /),'tab_nb',name_sub)
!  CALL alloc_NParray(tab_nq,(/ D /),'tab_nq',name_sub)
!
!
!  MarchPsi = RPsi
!  MarchPsi = ZERO
!
! IF (print_level > 0 .AND. BasisnD%para_SGType2%nb_SG > 10**5 ) THEN
!   write(out_unitp,'(a)')              'OpPsi SG4 (%): [--0-10-20-30-40-50-60-70-80-90-100]'
!   write(out_unitp,'(a)',ADVANCE='no') 'OpPsi SG4 (%): ['
!   CALL flush_perso(out_unitp)
! END IF
! psi0Hkpsi0(:) = CZERO
! E0            = para_H%E0
!
! 21 format(a,100(x,e12.5))
!
!!to be sure to have the correct number of threads
!nb_thread = BasisnD%para_SGType2%nb_threads
!nb_thread = 1
!
!
!!--------------------------------------------------------------
!!-- For the initialization of tab_l(:) and the use of ADD_ONE_TO_nDindex in the parallel loop
!ith = 0
!!$ ith = OMP_GET_THREAD_NUM()
!tab_l(:) = BasisnD%para_SGType2%nDval_init(:,ith+1)
!!--------------------------------------------------------------

!! we are not using the parallel do, to be able to use the correct initialized tab_l with nDval_init
!DO iG=BasisnD%para_SGType2%iG_th(ith+1),BasisnD%para_SGType2%fG_th(ith+1)
!
!   CALL ADD_ONE_TO_nDindex(BasisnD%para_SGType2%nDind_SmolyakRep,tab_l,iG=iG,err_sub=err_sub)
!
!   tab_nq(:) = getbis_tab_nq(tab_l,BasisnD%tab_basisPrimSG)
!   tab_nb(:) = getbis_tab_nb(tab_l,BasisnD%tab_basisPrimSG)
!
!   !write(6,*) 'iG',iG
!   !transfert part of the psi%RvecB(:) to PsiCvec%V and psi0%CvecB(:) to Psi0Cvec%V
!   CALL tabPackedBasis_TO_tabC_AT_iG(PsiCVec%V, Psi%CvecB, iG,BasisnD%para_SGType2)
!   CALL tabPackedBasis_TO_tabC_AT_iG(Psi0CVec%V,Psi0%CvecB,iG,BasisnD%para_SGType2)
!   nb = size(PsiRvec%V)
!
!   !auto-c is done on the basis
!   psi0Hkpsi0(0) = psi0Hkpsi0(0) + BasisnD%WeightSG(iG) * dot_product(Psi0Cvec%V,PsiCvec%V)
!
!   CALL alloc_TypeCVec(w1,nb)
!   CALL alloc_TypeCVec(w2(1),nb)
!
!
!   ! March_SG4
!   w1%V(:)     = Psivec%V
!   rtj         = CONE
!
!   DO j=1,para_propa%para_poly%npoly
!
!
!     w2(1)%V(:)  = w1%V
!
!     CALL sub_TabOpPsiC_OF_ONEDP_FOR_SGtype4(w2,iG,tab_l,para_H) ! in w2, we have H.Rw1
!
!     !E0 = dot_product(w1%V,w2(1)%V)
!
!        !write(6,21) 'w1',w1%V
!        !write(6,21) 'w2',w2(1)%V
!
!     w2(1)%V(:) = w2(1)%V - E0 * w1%V ! equivalent sub_scaledOpPsi
!
!        !write(6,21) 'w2',w2(1)%V
!
!     w1%V(:)    = w2(1)%V
!
!     !CALL Overlap_psi1_psi2(psi0Hkpsi0(j),psi0,w1)
!     psi0Hkpsi0(j) = psi0Hkpsi0(j) + BasisnD%WeightSG(iG) * dot_product(Psi0Cvec%V,w1%V)
!
!
!     rtj = rtj * cmplx(ZERO,-para_propa%WPdeltaT/real(j,kind=Rkind),kind=Rkind)
!
!     Rw2(1)%V(:) = Rw1%V * rtj
!
!        !write(6,21) 'w2*rtj',w2(1)%V
!
!     PsiCvec%V(:) = PsiCvec%V + w2(1)%V
!
!        !write(6,21) 'Rpsi',PsiRvec%V
!        !write(6,21) 'Ipsi',PsiIvec%V
!
!     norm2_w2 = dot_product(w2(1)%V,w2(1)%V)
!
!
!     IF (debug) write(out_unitp,*) 'iG,j,norm2_w2',iG,j,norm2_w2
!
!     IF (norm2_w2 > TEN**15) THEN
!       write(out_unitp,*) ' ERROR in ',name_sub
!       write(out_unitp,*) ' Norm of the vector is TOO large (> 10^15)',j,    &
!                                             norm2_w2
!       write(out_unitp,*) ' => Reduce the time step !!'
!       STOP
!     END IF
!     j_exit = j
!     IF (norm2_w2 < para_propa%para_poly%poly_tol) EXIT
!
!
!   END DO
!
!   CALL tabC_AT_iG_TO_tabPackedBasis(MarchPsi%CvecB,PsiCvec%V,         &
!                   iG,BasisnD%para_SGType2,BasisnD%WeightSG(iG))
!
!
!   write(out_unitp,*) 'iG,j,norm2_w2',iG,j,norm2_w2
!   !deallocate PsiRvec, Psi0Rvec
!   CALL dealloc_TypeCVec(w1)
!   CALL dealloc_TypeCVec(Cw2(1))
!
!   IF (print_level > 0  .AND. BasisnD%para_SGType2%nb_SG > 10**5 .AND. &
!       mod(iG,max(1,BasisnD%para_SGType2%nb_SG/10)) == 0) THEN
!     write(out_unitp,'(a)',ADVANCE='no') '---'
!     CALL flush_perso(out_unitp)
!   END IF
!
!   !write(6,*) 'iG done:',iG ; flush(6)
! END DO
!
! Psi%CvecB(:) = MarchPsi%CvecB
!
! CALL dealloc_psi(MarchRPsi)
! CALL dealloc_psi(MarchIPsi)
!
! CALL dealloc_NParray(tab_l ,'tab_l', name_sub)
! CALL dealloc_NParray(tab_nb,'tab_nb',name_sub)
! CALL dealloc_NParray(tab_nq,'tab_nq',name_sub)
!
! IF (print_level > 0 .AND. BasisnD%para_SGType2%nb_SG > 10**5) THEN
!   write(out_unitp,'(a)',ADVANCE='yes') '----]'
! END IF
! CALL flush_perso(out_unitp)
!
! IF (abs(norm2_w2) > para_propa%para_poly%poly_tol) THEN
!   write(out_unitp,*) ' ERROR in ',name_sub
!   write(out_unitp,*) ' Norm of the last vector is TOO large',norm2_w2
!   write(out_unitp,*) ' poly_tol: ',para_propa%para_poly%poly_tol
!   write(out_unitp,*) ' => npoly or max_poly are TOO small',               &
!                 para_propa%para_poly%npoly
!   STOP
! END IF
!
! !write(out_unitp,*) ' Psi before phase shift '
! !CALL ecri_psi(psi=psi)
! !- Phase Shift -----------------------------------
! phase = para_H%E0*para_propa%WPdeltaT
! psi   = psi * exp(-cmplx(ZERO,phase,kind=Rkind))
!
! !write(out_unitp,*) ' Psi after phase shift '
! !CALL ecri_psi(psi=psi)
!
!
! !- check norm ------------------
! CALL norme_psi(psi,GridRep=.FALSE.,BasisRep=.TRUE.,Renorm=.FALSE.)
! IF ( psi%norme > psi%max_norme) THEN
!   T  = T + para_propa%WPdeltaT
!   write(out_unitp,*) ' ERROR in ',name_sub
!   write(out_unitp,*) ' STOP propagation: norm > max_norm',psi%norme
!   para_propa%march_error   = .TRUE.
!   para_propa%test_max_norm = .TRUE.
!   STOP
! END IF
!
! CALL Overlap_psi1_psi2(cdot,psi0,psi)
! CALL Write_AutoCorr(no,T+para_propa%WPdeltaT,cdot)

!-----------------------------------------------------------
      IF (debug) THEN
        CALL norme_psi(psi)
        write(out_unitp,*) 'norm psi',psi%norme
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

 END SUBROUTINE march_noD_SG4_BasisRep
 SUBROUTINE march_noD_SG4_BasisRep_v0(T,no,psi,psi0,para_H,para_propa)
 USE mod_system
!$ USE omp_lib, only : OMP_GET_THREAD_NUM
 USE mod_nDindex
 USE mod_Coord_KEO,                ONLY : zmatrix
 USE mod_basis_set_alloc,          ONLY : basis
 USE mod_basis_BtoG_GtoB_SGType4,  ONLY : tabPackedBasis_TO_tabR_AT_iG, &
                  tabR_AT_iG_TO_tabPackedBasis,TypeRVec,alloc_TypeRVec,dealloc_TypeRVec, &
                  BDP_TO_GDP_OF_SmolyakRep,GDP_TO_BDP_OF_SmolyakRep, &
                  getbis_tab_nq,getbis_tab_nb

 USE mod_Op,              ONLY : param_Op,write_param_Op
 USE mod_OpPsi,           ONLY : sub_TabOpPsi_OF_ONEDP_FOR_SGtype4
 USE mod_psi_set_alloc,   ONLY : param_psi,copy_psi2TOpsi1,alloc_psi,dealloc_psi,ecri_psi
 USE mod_psi_Op,          ONLY : norme_psi,Overlap_psi1_psi2
 USE mod_psi_SimpleOp,    ONLY : operator (*),operator (+),operator (-),assignment (=)
 IMPLICIT NONE


!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_psi)   :: psi0,psi

!------ working variables ---------------------------------
      complex (kind=Rkind) :: E,cdot
      complex (kind=Rkind) :: psi0Hkpsi0(0:para_propa%para_poly%npoly)
      real (kind=Rkind)    :: microT,T,microdeltaT,phase,microphase
      integer              :: no

      complex (kind=Rkind) :: rt,rtj,rt_tmp
      integer              :: it,j,i_qaie_corr,j_exit
      real (kind=Rkind)    :: E0,rg,norm2_w2


 ! local variables
 TYPE (zmatrix), pointer :: mole
 TYPE (basis),   pointer :: BasisnD

 TYPE (TypeRVec)    :: PsiRVec,PsiIVec,Psi0Rvec,Psi0Ivec,Rw1,Rw2(1),Iw1,Iw2(1)
 TYPE (param_psi)   :: Rpsi,IPsi,Rpsi0,IPsi0,MarchRpsi,MarchIpsi


 !TYPE (TypeCVec)    :: PsiCVec,Psi0Cvec,Cw1,Cw2(1)

 integer :: ib,i,iG,nb_thread,itab,nq,nb,D,ith,err_sub

  integer,            allocatable       :: tab_nq(:)
  integer,            allocatable       :: tab_nb(:)
  integer,            allocatable       :: tab_l(:)

!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='march_noD_SG4_BasisRep_v0'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'deltaT',para_propa%WPdeltaT
        write(out_unitp,*) 'E0,Esc',para_H%E0,para_H%Esc
        write(out_unitp,*) 'nOD',para_propa%para_poly%npoly
      END IF
!-----------------------------------------------------------


  mole    => para_H%mole
  BasisnD => para_H%BasisnD

  D = BasisnD%nb_basis


  IF (para_H%nb_bie /= 1) STOP 'nb_bie /= 1'

  IF (OpPsi_omp == 0) THEN
    nb_thread = 1
  ELSE
    nb_thread = OpPsi_maxth
  END IF

  CALL alloc_NParray(tab_l ,(/ D /),'tab_l', name_sub)
  CALL alloc_NParray(tab_nb,(/ D /),'tab_nb',name_sub)
  CALL alloc_NParray(tab_nq,(/ D /),'tab_nq',name_sub)

  CALL copy_psi2TOpsi1(RPsi,Psi,alloc=.FALSE.)
  RPsi%cplx = .FALSE.
  CALL alloc_psi(RPsi,BasisRep=.TRUE.)
  IPsi  = RPsi
  RPsi0 = RPsi
  IPsi0 = RPsi
  MarchRpsi = RPsi ; MarchRpsi = ZERO
  MarchIpsi = IPsi ; MarchIpsi = ZERO

  RPsi%RvecB(:) = Real(Psi%CvecB(:),kind=Rkind)
  IPsi%RvecB(:) = Aimag(Psi%CvecB(:))

  RPsi0%RvecB(:) = Real(Psi0%CvecB(:),kind=Rkind)
  IPsi0%RvecB(:) = Aimag(Psi0%CvecB(:))

 IF (print_level > 0 .AND. BasisnD%para_SGType2%nb_SG > 10**5 ) THEN
   write(out_unitp,'(a)')              'OpPsi SG4 (%): [--0-10-20-30-40-50-60-70-80-90-100]'
   write(out_unitp,'(a)',ADVANCE='no') 'OpPsi SG4 (%): ['
   CALL flush_perso(out_unitp)
 END IF
 psi0Hkpsi0(:) = CZERO
 E0            = para_H%E0

 21 format(a,100(x,e12.5))

!to be sure to have the correct number of threads
nb_thread = BasisnD%para_SGType2%nb_threads
nb_thread = 1


 !--------------------------------------------------------------
 !-- For the initialization of tab_l(:) and the use of ADD_ONE_TO_nDindex in the parallel loop
 ith = 0
 !$ ith = OMP_GET_THREAD_NUM()
 tab_l(:) = BasisnD%para_SGType2%nDval_init(:,ith+1)
 !--------------------------------------------------------------

 ! we are not using the parallel do, to be able to use the correct initialized tab_l with nDval_init
 DO iG=BasisnD%para_SGType2%iG_th(ith+1),BasisnD%para_SGType2%fG_th(ith+1)

   CALL ADD_ONE_TO_nDindex(BasisnD%para_SGType2%nDind_SmolyakRep,tab_l,iG=iG,err_sub=err_sub)

   tab_nq(:) = getbis_tab_nq(tab_l,BasisnD%tab_basisPrimSG)
   tab_nb(:) = getbis_tab_nb(tab_l,BasisnD%tab_basisPrimSG)

   !write(6,*) 'iG',iG
   !transfert part of the psi%RvecB(:) to PsiRvec%R and psi0%RvecB(:) to Psi0Rvec%R
   ! the real and imaginary part are splited
   CALL tabPackedBasis_TO_tabR_AT_iG(PsiRvec%V, RPsi%RvecB, iG,BasisnD%para_SGType2)
   CALL tabPackedBasis_TO_tabR_AT_iG(PsiIvec%V, IPsi%RvecB, iG,BasisnD%para_SGType2)
   CALL tabPackedBasis_TO_tabR_AT_iG(Psi0Rvec%V,RPsi0%RvecB,iG,BasisnD%para_SGType2)
   CALL tabPackedBasis_TO_tabR_AT_iG(Psi0Ivec%V,IPsi0%RvecB,iG,BasisnD%para_SGType2)
   nb = size(PsiRvec%V)

   !auto-c is done on the basis
   psi0Hkpsi0(0) = psi0Hkpsi0(0) + BasisnD%WeightSG(iG) *     cmplx(        &
     dot_product(Psi0Rvec%V,PsiRvec%V)+dot_product(Psi0Ivec%V,PsiIvec%V) ,  &
     dot_product(Psi0Rvec%V,PsiIvec%V)-dot_product(Psi0Ivec%V,PsiRvec%V), kind=Rkind)

   CALL alloc_TypeRVec(Rw1,nb)
   CALL alloc_TypeRVec(Iw1,nb)
   CALL alloc_TypeRVec(Rw2(1),nb)
   CALL alloc_TypeRVec(Iw2(1),nb)


   ! March_SG4


   Rw1%V(:)     = PsiRvec%V
   Iw1%V(:)     = PsiIvec%V

       !write(6,21) 'Rw1',Rw1%V
       !write(6,21) 'Iw1',Iw1%V

   rtj          = CONE

   DO j=1,para_propa%para_poly%npoly


     Rw2(1)%V(:)  = Rw1%V
     Iw2(1)%V(:)  = Iw1%V

     CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4(Rw2,iG,tab_l,para_H) ! in w2, we have H.Rw1
     CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4(Iw2,iG,tab_l,para_H) ! in w2, we have H.Iw1

     !E0 = dot_product(Rw1%V,Rw2(1)%V)+dot_product(Iw1%V,Iw2(1)%V)

        !write(6,21) 'Rw1',Rw1%V
        !write(6,21) 'Iw1',Iw1%V
        !write(6,21) 'Rw2',Rw2(1)%V
        !write(6,21) 'Iw2',Iw2(1)%V

     Rw2(1)%V(:) = Rw2(1)%V - E0 * Rw1%V ! equivalent sub_scaledOpPsi
     Iw2(1)%V(:) = Iw2(1)%V - E0 * Iw1%V ! equivalent sub_scaledOpPsi

        !write(6,21) 'Rw2',Rw2(1)%V
        !write(6,21) 'Iw2',Iw2(1)%V

     Rw1%V(:)    = Rw2(1)%V
     Iw1%V(:)    = Iw2(1)%V

     !CALL Overlap_psi1_psi2(psi0Hkpsi0(j),psi0,w1)
     psi0Hkpsi0(j) = psi0Hkpsi0(j) + BasisnD%WeightSG(iG) * cmplx(      &
       dot_product(Psi0Rvec%V,Rw1%V)+dot_product(Psi0Ivec%V,Iw1%V) ,    &
       dot_product(Psi0Rvec%V,Iw1%V)-dot_product(Psi0Ivec%V,Rw1%V), kind=Rkind)


     rtj = rtj * cmplx(ZERO,-para_propa%WPdeltaT/real(j,kind=Rkind),kind=Rkind)

     Rw2(1)%V(:) = Rw1%V * real(rtj,kind=Rkind) - Iw1%V * Aimag(rtj)
     Iw2(1)%V(:) = Rw1%V * Aimag(rtj)           + Iw1%V * real(rtj,kind=Rkind)

        !write(6,21) 'Rw2*rtj',Rw2(1)%V
        !write(6,21) 'Iw2*rtj',Iw2(1)%V

     PsiRvec%V(:) = PsiRvec%V + Rw2(1)%V
     PsiIvec%V(:) = PsiIvec%V + Iw2(1)%V

        !write(6,21) 'Rpsi',PsiRvec%V
        !write(6,21) 'Ipsi',PsiIvec%V

     norm2_w2 = dot_product(Rw2(1)%V,Rw2(1)%V) + dot_product(Iw2(1)%V,Iw2(1)%V)


     IF (debug) write(out_unitp,*) 'iG,j,norm2_w2',iG,j,norm2_w2

     IF (norm2_w2 > TEN**15) THEN
       write(out_unitp,*) ' ERROR in ',name_sub
       write(out_unitp,*) ' Norm of the vector is TOO large (> 10^15)',j,    &
                                             norm2_w2
       write(out_unitp,*) ' => Reduce the time step !!'
       STOP
     END IF
     j_exit = j
     IF (norm2_w2 < para_propa%para_poly%poly_tol) EXIT


   END DO

   CALL tabR_AT_iG_TO_tabPackedBasis(MarchRPsi%RvecB,PsiRvec%V,         &
                   iG,BasisnD%para_SGType2,BasisnD%WeightSG(iG))
   CALL tabR_AT_iG_TO_tabPackedBasis(MarchIPsi%RvecB,PsiIvec%V,         &
                   iG,BasisnD%para_SGType2,BasisnD%WeightSG(iG))


   write(out_unitp,*) 'iG,j,norm2_w2',iG,j,norm2_w2
   !deallocate PsiRvec, Psi0Rvec
   CALL dealloc_TypeRVec(PsiRvec)
   CALL dealloc_TypeRVec(PsiIvec)
   CALL dealloc_TypeRVec(Psi0Rvec)
   CALL dealloc_TypeRVec(Psi0Ivec)
   CALL dealloc_TypeRVec(Rw1)
   CALL dealloc_TypeRVec(Iw1)
   CALL dealloc_TypeRVec(Rw2(1))
   CALL dealloc_TypeRVec(Iw2(1))

   IF (print_level > 0  .AND. BasisnD%para_SGType2%nb_SG > 10**5 .AND. &
       mod(iG,max(1,BasisnD%para_SGType2%nb_SG/10)) == 0) THEN
     write(out_unitp,'(a)',ADVANCE='no') '---'
     CALL flush_perso(out_unitp)
   END IF

   !write(6,*) 'iG done:',iG ; flush(6)
 END DO

 Psi%CvecB(:) = cmplx(MarchRPsi%RvecB,MarchIPsi%RvecB,kind=Rkind)

 CALL dealloc_psi(RPsi)
 CALL dealloc_psi(IPsi)
 CALL dealloc_psi(RPsi0)
 CALL dealloc_psi(IPsi0)
 CALL dealloc_psi(MarchRPsi)
 CALL dealloc_psi(MarchIPsi)

 CALL dealloc_NParray(tab_l ,'tab_l', name_sub)
 CALL dealloc_NParray(tab_nb,'tab_nb',name_sub)
 CALL dealloc_NParray(tab_nq,'tab_nq',name_sub)

 IF (print_level > 0 .AND. BasisnD%para_SGType2%nb_SG > 10**5) THEN
   write(out_unitp,'(a)',ADVANCE='yes') '----]'
 END IF
 CALL flush_perso(out_unitp)

 IF (abs(norm2_w2) > para_propa%para_poly%poly_tol) THEN
   write(out_unitp,*) ' ERROR in ',name_sub
   write(out_unitp,*) ' Norm of the last vector is TOO large',norm2_w2
   write(out_unitp,*) ' poly_tol: ',para_propa%para_poly%poly_tol
   write(out_unitp,*) ' => npoly or max_poly are TOO small',               &
                 para_propa%para_poly%npoly
   STOP
 END IF

 !write(out_unitp,*) ' Psi before phase shift '
 !CALL ecri_psi(psi=psi)
 !- Phase Shift -----------------------------------
 phase = para_H%E0*para_propa%WPdeltaT
 psi = psi * exp(-cmplx(ZERO,phase,kind=Rkind))

 !write(out_unitp,*) ' Psi after phase shift '
 !CALL ecri_psi(psi=psi)


 !- check norm ------------------
 CALL norme_psi(psi,GridRep=.FALSE.,BasisRep=.TRUE.,Renorm=.FALSE.)
 IF ( psi%norme > psi%max_norme) THEN
   T  = T + para_propa%WPdeltaT
   write(out_unitp,*) ' ERROR in ',name_sub
   write(out_unitp,*) ' STOP propagation: norm > max_norm',psi%norme
   para_propa%march_error   = .TRUE.
   para_propa%test_max_norm = .TRUE.
   STOP
 END IF

 CALL Overlap_psi1_psi2(cdot,psi0,psi)
 CALL Write_AutoCorr(no,T+para_propa%WPdeltaT,cdot)



!-----------------------------------------------------------
      IF (debug) THEN
        CALL norme_psi(psi)
        write(out_unitp,*) 'norm psi',psi%norme
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE march_noD_SG4_BasisRep_v0
 SUBROUTINE march_noD_SG4_GridRep(T,no,psi,psi0,para_H,para_propa)
 USE mod_system
 USE mod_nDindex

 USE mod_Coord_KEO,                ONLY : zmatrix, get_Qact, get_d0GG

 USE mod_basis_set_alloc,          ONLY : basis
 USE mod_basis,                    ONLY : Rec_WrhonD,Rec_Qact_SG4
 USE mod_basis_BtoG_GtoB_SGType4,  ONLY : tabPackedBasis_TO_tabR_AT_iG, &
                                          tabR_AT_iG_TO_tabPackedBasis, &
                                          tabR2grid_TO_tabR1_AT_iG,     &
                              TypeRVec,alloc_TypeRVec,dealloc_TypeRVec, &
                              TypeCVec,alloc_TypeCVec,dealloc_TypeCVec, &
                     BDP_TO_GDP_OF_SmolyakRep,GDP_TO_BDP_OF_SmolyakRep, &
                              getbis_tab_nq,getbis_tab_nb

 USE mod_Op,              ONLY : param_Op,write_param_Op
 USE mod_OpPsi,           ONLY : sub_TabOpPsi_OF_ONEGDP_WithOp_FOR_SGtype4
 USE mod_psi_set_alloc,   ONLY : param_psi,copy_psi2TOpsi1,alloc_psi,dealloc_psi,ecri_psi
 USE mod_psi_Op,          ONLY : norme_psi
 USE mod_psi_SimpleOp,    ONLY : operator (*),operator (+),operator (-),assignment (=)
 IMPLICIT NONE


!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_psi)   :: psi0,psi

!------ working variables ---------------------------------
      complex (kind=Rkind) :: E,cdot
      complex (kind=Rkind) :: psi0Hkpsi0(0:para_propa%para_poly%npoly)
      real (kind=Rkind)    :: microT,T,microdeltaT,phase,microphase
      integer              :: no

      complex (kind=Rkind) :: rt,rtj,rt_tmp
      integer              :: it,j,i_qaie_corr,j_exit
      real (kind=Rkind)    :: rg,norm2_w2



 ! local variables
 TYPE (zmatrix), pointer :: mole
 TYPE (basis),   pointer :: BasisnD

 TYPE (TypeCVec)    :: PsiRVec,PsiIVec,Psi0Rvec,Psi0Ivec,Rw1,Rw2(1),Iw1,Iw2(1)
 TYPE (TypeRVec)    :: Wrho

 TYPE (param_psi)   :: Marchpsi

 TYPE (param_psi)   :: Rpsi,IPsi,Rpsi0,IPsi0,MarchRpsi,MarchIpsi

 integer :: ib,i,iG,nb_thread,itab,nq,nb,D,iq,iq_ib,L,iterm00,err_sub

  integer,            allocatable       :: tab_nq(:)
  integer,            allocatable       :: tab_iq(:)
  integer,            allocatable       :: tab_nb(:)
  integer,            allocatable       :: tab_l(:)

  real (kind=Rkind),  allocatable       :: V(:)
  real (kind=Rkind),  allocatable       :: Qact(:)
  TYPE (Type_nDindex)                   :: nDind_DPG    ! multidimensional DP index
  real (kind=Rkind),  allocatable       :: GG(:,:,:)
  real (kind=Rkind),  allocatable       :: sqRhoOVERJac(:),Jac(:)
  real (kind=Rkind)                     :: Rho

!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='march_noD_SG4_GridRep'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'deltaT',para_propa%WPdeltaT
        write(out_unitp,*) 'E0,Esc',para_H%E0,para_H%Esc
        write(out_unitp,*) 'nOD',para_propa%para_poly%npoly
      END IF
!-----------------------------------------------------------


  mole    => para_H%mole
  BasisnD => para_H%BasisnD

  D = size(BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval,dim=1)


  IF (para_H%nb_bie /= 1) STOP 'nb_bie /= 1'

  IF (OpPsi_omp == 0) THEN
    nb_thread = 1
  ELSE
    nb_thread = OpPsi_maxth
  END IF

  CALL alloc_NParray(tab_l ,(/ D /),'tab_l', name_sub)
  CALL alloc_NParray(tab_nb,(/ D /),'tab_nb',name_sub)
  CALL alloc_NParray(tab_nq,(/ D /),'tab_nq',name_sub)
  CALL alloc_NParray(tab_iq,(/ D /),'tab_iq',name_sub)

 Marchpsi = Psi

  CALL copy_psi2TOpsi1(RPsi,Psi,alloc=.FALSE.)
  RPsi%cplx = .FALSE.
  CALL alloc_psi(RPsi,BasisRep=.TRUE.)
  IPsi  = RPsi
  RPsi0 = RPsi
  IPsi0 = RPsi
  MarchRpsi = RPsi ; MarchRpsi = ZERO
  MarchIpsi = IPsi ; MarchIpsi = ZERO

  RPsi%RvecB(:) = Real(Psi%CvecB(:),kind=Rkind)
  IPsi%RvecB(:) = Aimag(Psi%CvecB(:))

  RPsi0%RvecB(:) = Real(Psi0%CvecB(:),kind=Rkind)
  IPsi0%RvecB(:) = Aimag(Psi0%CvecB(:))

 IF (print_level > 0 .AND. BasisnD%para_SGType2%nb_SG > 10**5 ) THEN
   write(out_unitp,'(a)')              'OpPsi SG4 (%): [--0-10-20-30-40-50-60-70-80-90-100]'
   write(out_unitp,'(a)',ADVANCE='no') 'OpPsi SG4 (%): ['
   CALL flush_perso(out_unitp)
 END IF
 psi0Hkpsi0(:) = CZERO

 21 format(a,100(x,f8.6))

 DO iG=1,BasisnD%para_SGType2%nb_SG

   tab_l(:)  = BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval(:,iG)
   tab_nq(:) = getbis_tab_nq(tab_l,BasisnD%tab_basisPrimSG)
   tab_nb(:) = getbis_tab_nb(tab_l,BasisnD%tab_basisPrimSG)

   !write(6,*) 'iG',iG
   !transfert part of the psi%RvecB(:) to PsiRvec%R and psi0%RvecB(:) to Psi0Rvec%V
   ! the real and imaginary part are splited
STOP 'not yet !!!!'
!   CALL tabPackedBasis_TO_tabR_AT_iG(PsiRvec%R, RPsi%RvecB, iG,BasisnD%para_SGType2)
!   CALL tabPackedBasis_TO_tabR_AT_iG(PsiIvec%R, IPsi%RvecB, iG,BasisnD%para_SGType2)
!   CALL tabPackedBasis_TO_tabR_AT_iG(Psi0Rvec%R,RPsi0%RvecB,iG,BasisnD%para_SGType2)
!   CALL tabPackedBasis_TO_tabR_AT_iG(Psi0Ivec%R,IPsi0%RvecB,iG,BasisnD%para_SGType2)
!
!   ! here, auto-c is done on the basis
!   psi0Hkpsi0(0) = psi0Hkpsi0(0) + BasisnD%WeightSG(iG) *     cmplx(        &
!     dot_product(Psi0Rvec%R,PsiRvec%R)+dot_product(Psi0Ivec%R,PsiIvec%R) ,  &
!     dot_product(Psi0Rvec%R,PsiIvec%R)-dot_product(Psi0Ivec%R,PsiRvec%R), kind=Rkind)
!
!   ! partial B to G
!   CALL BDP_TO_GDP_OF_SmolyakRep(PsiRvec%R,BasisnD%tab_basisPrimSG,tab_l,tab_nq,tab_nb)
!   CALL BDP_TO_GDP_OF_SmolyakRep(PsiIvec%R,BasisnD%tab_basisPrimSG,tab_l,tab_nq,tab_nb)
!
!   CALL BDP_TO_GDP_OF_SmolyakRep(Psi0Rvec%R,BasisnD%tab_basisPrimSG,tab_l,tab_nq,tab_nb)
!   CALL BDP_TO_GDP_OF_SmolyakRep(Psi0Ivec%R,BasisnD%tab_basisPrimSG,tab_l,tab_nq,tab_nb)
!
!   nq = size(PsiRvec%R)
!   CALL alloc_TypeRVec(Rw1,nq)
!   CALL alloc_TypeRVec(Iw1,nq)
!   CALL alloc_TypeRVec(Rw2(1),nq)
!   CALL alloc_TypeRVec(Iw2(1),nq)
!   CALL alloc_TypeRVec(Wrho,nq)
!
!
!   !transfert part of the potential
!   iterm00 = para_H%derive_term_TO_iterm(0,0)
!
!   IF (para_H%OpGrid(iterm00)%grid_zero) THEN
!     CALL alloc_NParray(V,(/ nq /),'V',name_sub)
!     V(:) = ZERO
!   ELSE
!     IF (para_H%OpGrid(iterm00)%grid_cte) THEN
!       CALL alloc_NParray(V,(/ nq /),'V',name_sub)
!       V(:) = para_H%OpGrid(iterm00)%Mat_cte(1,1)
!     ELSE
!       CALL tabR2grid_TO_tabR1_AT_iG(V,para_H%OpGrid(iterm00)%Grid(:,1,1),&
!                                     iG,BasisnD%para_SGType2)
!     END IF
!   END IF
!
!   ! G calculation
!   CALL alloc_NParray(sqRhoOVERJac,(/ nq /),      'sqRhoOVERJac',name_sub)
!   CALL alloc_NParray(Jac,         (/ nq /),      'Jac',         name_sub)
!   CALL alloc_NParray(GG,(/nq,mole%nb_act1,mole%nb_act1/),'GGiq',name_sub)
!   CALL alloc_NParray(Qact,       (/mole%nb_var/),'Qact',        name_sub)
!   CALL init_nDindexPrim(nDind_DPG,BasisnD%nb_basis,tab_nq,type_OF_nDindex=-1)
!   DO iq=1,nq
!     CALL get_Qact(Qact,mole%ActiveTransfo) ! rigid, flexible coordinates
!     CALL Rec_Qact_SG4(Qact,BasisnD%tab_basisPrimSG,tab_l,nDind_DPG,iq,mole,err_sub)
!
!     CALL get_d0GG(Qact,para_H%para_Tnum,mole,d0GG=GG(iq,:,:),          &
!                                         def=.TRUE.,Jac=Jac(iq),Rho=Rho)
!     sqRhoOVERJac(iq) = sqrt(Rho/Jac(iq))
!
!   END DO
!
!   DO iq=1,nq
!     !Wrho%R(iq) = BasisnD%WeightSG(iG)
!     Wrho%R(iq) = ONE
!     CALL calc_nDindex(nDind_DPG,iq,tab_iq)
!     DO ib=1,BasisnD%nb_basis
!       L     = tab_l(ib)
!       iq_ib = tab_iq(ib)
!       Wrho%R(iq) = Wrho%R(iq) * Rec_WrhonD(BasisnD%tab_basisPrimSG(L,ib),iq_ib)
!     END DO
!   END DO
!   CALL dealloc_nDindex(nDind_DPG)
!   ! March_SG4
!   Rw1%R(:)     = PsiRvec%R
!   Iw1%R(:)     = PsiIvec%R
!
!   rtj          = CONE
!
!   DO j=1,para_propa%para_poly%npoly
!
!     Rw2(1)%R(:)  = Rw1%R
!     Iw2(1)%R(:)  = Iw1%R
!     !CALL sub_TabOpPsi_OF_ONEGDP_FOR_SGtype4(Rw2,iG,para_H) ! in Rw2, we have H.Rw1
!     !CALL sub_TabOpPsi_OF_ONEGDP_FOR_SGtype4(Iw2,iG,para_H) ! in Iw2, we have H.Iw1
!
!     CALL sub_TabOpPsi_OF_ONEGDP_WithOp_FOR_SGtype4(Rw2,iG,para_H,      &
!                                                    V,GG,sqRhoOVERJac,Jac)
!     CALL sub_TabOpPsi_OF_ONEGDP_WithOp_FOR_SGtype4(Iw2,iG,para_H,      &
!                                                    V,GG,sqRhoOVERJac,Jac)
!
!        !write(6,21) 'Rw1',Rw1%R
!        !write(6,21) 'Iw1',Iw1%R
!        !write(6,21) 'Rw2',Rw2(1)%R
!        !write(6,21) 'Iw2',Iw2(1)%R
!
!     Rw2(1)%R(:) = Rw2(1)%R - para_H%E0 * Rw1%R ! equivalent sub_scaledOpPsi
!     Iw2(1)%R(:) = Iw2(1)%R - para_H%E0 * Iw1%R ! equivalent sub_scaledOpPsi
!
!        !write(6,21) 'Rw2',Rw2(1)%R
!        !write(6,21) 'Iw2',Iw2(1)%R
!
!     Rw1%R(:)    = Rw2(1)%R
!     Iw1%R(:)    = Iw2(1)%R
!
!     !CALL Overlap_psi1_psi2(psi0Hkpsi0(j),psi0,w1)
!     !wrong to be adapted for the grid
!     psi0Hkpsi0(j) = psi0Hkpsi0(j) + BasisnD%WeightSG(iG) * cmplx(      &
!       dot_product(Psi0Rvec%R,Wrho%R*Rw1%R)+dot_product(Psi0Ivec%R,Wrho%R*Iw1%R) ,    &
!       dot_product(Psi0Rvec%R,Wrho%R*Iw1%R)-dot_product(Psi0Ivec%R,Wrho%R*Rw1%R), kind=Rkind)
!
!     rtj = rtj * cmplx(ZERO,-para_propa%WPdeltaT/real(j,kind=Rkind),kind=Rkind)
!
!     Rw2(1)%R(:) = Rw1%R * real(rtj,kind=Rkind) - Iw1%R * Aimag(rtj)
!     Iw2(1)%R(:) = Rw1%R * Aimag(rtj)           + Iw1%R * real(rtj,kind=Rkind)
!
!     PsiRvec%R(:) = PsiRvec%R + Rw2(1)%R
!     PsiIvec%R(:) = PsiIvec%R + Iw2(1)%R
!
!     norm2_w2 = dot_product(Rw2(1)%R,Wrho%R*Rw2(1)%R) + dot_product(Iw2(1)%R,Wrho%R*Iw2(1)%R)
!
!
!     IF (debug) write(out_unitp,*) 'iG,j,norm2_w2',iG,j,norm2_w2
!
!     IF (norm2_w2 > TEN**15) THEN
!       write(out_unitp,*) ' ERROR in ',name_sub
!       write(out_unitp,*) ' Norm of the vector is TOO large (> 10^15)',j,    &
!                                             norm2_w2
!       write(out_unitp,*) ' => Reduce the time step !!'
!       STOP
!     END IF
!     j_exit = j
!     IF (norm2_w2 < para_propa%para_poly%poly_tol) EXIT
!
!
!   END DO
!   write(out_unitp,*) 'iG,j,norm2_w2',iG,j,norm2_w2
!
!   CALL GDP_TO_BDP_OF_SmolyakRep(PsiRvec%R,BasisnD%tab_basisPrimSG,tab_l,tab_nq,tab_nb)
!   CALL GDP_TO_BDP_OF_SmolyakRep(PsiIvec%R,BasisnD%tab_basisPrimSG,tab_l,tab_nq,tab_nb)
!
!
!   !transfert back, PsiRvec%R to the psi%RvecB(:)
!   CALL tabR_AT_iG_TO_tabPackedBasis(MarchRPsi%RvecB,PsiRvec%R,         &
!                   iG,BasisnD%para_SGType2,BasisnD%WeightSG(iG))
!   CALL tabR_AT_iG_TO_tabPackedBasis(MarchIPsi%RvecB,PsiIvec%R,         &
!                   iG,BasisnD%para_SGType2,BasisnD%WeightSG(iG))
!
!   !deallocate PsiRvec, Psi0Rvec
!   CALL dealloc_TypeRVec(PsiRvec)
!   CALL dealloc_TypeRVec(PsiIvec)
!   CALL dealloc_TypeRVec(Psi0Rvec)
!   CALL dealloc_TypeRVec(Psi0Ivec)
!   CALL dealloc_TypeRVec(Rw1)
!   CALL dealloc_TypeRVec(Iw1)
!   CALL dealloc_TypeRVec(Rw2(1))
!   CALL dealloc_TypeRVec(Iw2(1))
!   CALL dealloc_TypeRVec(Wrho)
!   CALL dealloc_NParray(V,           'V',           name_sub)
!   CALL dealloc_NParray(sqRhoOVERJac,'sqRhoOVERJac',name_sub)
!   CALL dealloc_NParray(Jac,         'Jac',         name_sub)
!   CALL dealloc_NParray(GG,          'GG',          name_sub)
!   CALL dealloc_NParray(Qact,        'Qact',        name_sub)
!
!
!   IF (print_level > 0  .AND. BasisnD%para_SGType2%nb_SG > 10**5 .AND. &
!       mod(iG,max(1,BasisnD%para_SGType2%nb_SG/10)) == 0) THEN
!     write(out_unitp,'(a)',ADVANCE='no') '---'
!     CALL flush_perso(out_unitp)
!   END IF
!
!   !write(6,*) 'iG done:',iG ; flush(6)
 END DO
!
! Psi%CvecB(:) = cmplx(MarchRPsi%RvecB,MarchIPsi%RvecB,kind=Rkind)
!
!
! CALL dealloc_psi(RPsi)
! CALL dealloc_psi(IPsi)
! CALL dealloc_psi(RPsi0)
! CALL dealloc_psi(IPsi0)
! CALL dealloc_psi(MarchRPsi)
! CALL dealloc_psi(MarchIPsi)
!
!
! CALL dealloc_NParray(tab_l ,'tab_l', name_sub)
! CALL dealloc_NParray(tab_nb,'tab_nb',name_sub)
! CALL dealloc_NParray(tab_nq,'tab_nq',name_sub)
! CALL dealloc_NParray(tab_iq,'tab_nq',name_sub)
!
! IF (print_level > 0 .AND. BasisnD%para_SGType2%nb_SG > 10**5) THEN
!   write(out_unitp,'(a)',ADVANCE='yes') '----]'
! END IF
! CALL flush_perso(out_unitp)
!
! IF (abs(norm2_w2) > para_propa%para_poly%poly_tol) THEN
!   write(out_unitp,*) ' ERROR in ',name_sub
!   write(out_unitp,*) ' Norm of the last vector is TOO large',norm2_w2
!   write(out_unitp,*) ' poly_tol: ',para_propa%para_poly%poly_tol
!   write(out_unitp,*) ' => npoly or max_poly are TOO small',               &
!                 para_propa%para_poly%npoly
!   STOP
! END IF
!
! !- Phase Shift -----------------------------------
! phase = para_H%E0*para_propa%WPdeltaT
! psi = psi * exp(-cmplx(ZERO,phase,kind=Rkind))
!
! !- check norm ------------------
! CALL norme_psi(psi,GridRep=.FALSE.,BasisRep=.TRUE.,Renorm=.FALSE.)
! IF ( psi%norme > psi%max_norme) THEN
!   T  = T + para_propa%WPdeltaT
!   write(out_unitp,*) ' ERROR in ',name_sub
!   write(out_unitp,*) ' STOP propagation: norm > max_norm',psi%norme
!   para_propa%march_error   = .TRUE.
!   para_propa%test_max_norm = .TRUE.
!   STOP
! END IF
!
!
!
! microdeltaT = para_propa%WPdeltaT/                                &
!               real(para_propa%nb_micro,kind=Rkind)
! microphase = phase/                                               &
!               real(para_propa%nb_micro,kind=Rkind)
!
! phase = ZERO
! microT = ZERO
!
! DO it=1,para_propa%nb_micro
!
!   microT = microT + microdeltaT
!   phase = phase + microphase
!
!   rtj = cmplx(ONE,ZERO,kind=Rkind)
!   cdot = psi0Hkpsi0(0)
!   DO j=1,j_exit
!
!     rt_tmp =  cmplx(ZERO,-microT/real(j,kind=Rkind),kind=Rkind)
!     rtj = rt_tmp * rtj
!     cdot = cdot + psi0Hkpsi0(j) * rtj
!   END DO
!   cdot = cdot * cmplx( cos(phase),-sin(phase),kind=Rkind)
!   CALL Write_AutoCorr(no,T+microT,cdot)
! END DO
!
!

!-----------------------------------------------------------
      IF (debug) THEN
        CALL norme_psi(psi)
        write(out_unitp,*) 'norm psi',psi%norme
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------


 END SUBROUTINE march_noD_SG4_GridRep
 SUBROUTINE march_noD_SG4_GridRep_v0(T,no,psi,psi0,para_H,para_propa)
 USE mod_system
 USE mod_nDindex

 USE mod_Coord_KEO,                     ONLY : zmatrix, get_Qact, get_d0GG

 USE mod_basis_set_alloc,          ONLY : basis
 USE mod_basis,                    ONLY : Rec_WrhonD,Rec_Qact_SG4
 USE mod_basis_BtoG_GtoB_SGType4,  ONLY : tabPackedBasis_TO_tabR_AT_iG, &
                                          tabR_AT_iG_TO_tabPackedBasis, &
                                          tabR2grid_TO_tabR1_AT_iG,     &
                              TypeRVec,alloc_TypeRVec,dealloc_TypeRVec, &
                     BDP_TO_GDP_OF_SmolyakRep,GDP_TO_BDP_OF_SmolyakRep, &
                              getbis_tab_nq,getbis_tab_nb

 USE mod_Op,              ONLY : param_Op,write_param_Op
 USE mod_OpPsi,           ONLY : sub_TabOpPsi_OF_ONEGDP_WithOp_FOR_SGtype4
 USE mod_psi_set_alloc,   ONLY : param_psi,copy_psi2TOpsi1,alloc_psi,dealloc_psi,ecri_psi
 USE mod_psi_Op,          ONLY : norme_psi
 USE mod_psi_SimpleOp,    ONLY : operator (*),operator (+),operator (-),assignment (=)
 IMPLICIT NONE


!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_psi)   :: psi0,psi

!------ working variables ---------------------------------
      complex (kind=Rkind) :: E,cdot
      complex (kind=Rkind) :: psi0Hkpsi0(0:para_propa%para_poly%npoly)
      real (kind=Rkind)    :: microT,T,microdeltaT,phase,microphase
      integer              :: no

      complex (kind=Rkind) :: rt,rtj,rt_tmp
      integer              :: it,j,i_qaie_corr,j_exit
      real (kind=Rkind)    :: rg,norm2_w2



 ! local variables
 TYPE (zmatrix), pointer :: mole
 TYPE (basis),   pointer :: BasisnD

 TYPE (TypeRVec)    :: PsiRVec,PsiIVec,Psi0Rvec,Psi0Ivec,Rw1,Rw2(1),Iw1,Iw2(1),Wrho
 TYPE (param_psi)   :: Rpsi,IPsi,Rpsi0,IPsi0,MarchRpsi,MarchIpsi

 integer :: ib,i,iG,nb_thread,itab,nq,nb,D,iq,iq_ib,L,iterm00,err_sub

  integer,            allocatable       :: tab_nq(:)
  integer,            allocatable       :: tab_iq(:)
  integer,            allocatable       :: tab_nb(:)
  integer,            allocatable       :: tab_l(:)

  real (kind=Rkind),  allocatable       :: V(:)
  real (kind=Rkind),  allocatable       :: Qact(:)
  TYPE (Type_nDindex)                   :: nDind_DPG    ! multidimensional DP index
  real (kind=Rkind),  allocatable       :: GG(:,:,:)
  real (kind=Rkind),  allocatable       :: sqRhoOVERJac(:),Jac(:)
  real (kind=Rkind)                     :: Rho

!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='march_noD_SG4_GridRep_v0'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'deltaT',para_propa%WPdeltaT
        write(out_unitp,*) 'E0,Esc',para_H%E0,para_H%Esc
        write(out_unitp,*) 'nOD',para_propa%para_poly%npoly
      END IF
!-----------------------------------------------------------


  mole    => para_H%mole
  BasisnD => para_H%BasisnD

  D = size(BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval,dim=1)


  IF (para_H%nb_bie /= 1) STOP 'nb_bie /= 1'

  IF (OpPsi_omp == 0) THEN
    nb_thread = 1
  ELSE
    nb_thread = OpPsi_maxth
  END IF

  CALL alloc_NParray(tab_l ,(/ D /),'tab_l', name_sub)
  CALL alloc_NParray(tab_nb,(/ D /),'tab_nb',name_sub)
  CALL alloc_NParray(tab_nq,(/ D /),'tab_nq',name_sub)
  CALL alloc_NParray(tab_iq,(/ D /),'tab_iq',name_sub)

  CALL copy_psi2TOpsi1(RPsi,Psi,alloc=.FALSE.)
  RPsi%cplx = .FALSE.
  CALL alloc_psi(RPsi,BasisRep=.TRUE.)
  IPsi  = RPsi
  RPsi0 = RPsi
  IPsi0 = RPsi
  MarchRpsi = RPsi ; MarchRpsi = ZERO
  MarchIpsi = IPsi ; MarchIpsi = ZERO

  RPsi%RvecB(:) = Real(Psi%CvecB(:),kind=Rkind)
  IPsi%RvecB(:) = Aimag(Psi%CvecB(:))

  RPsi0%RvecB(:) = Real(Psi0%CvecB(:),kind=Rkind)
  IPsi0%RvecB(:) = Aimag(Psi0%CvecB(:))

 IF (print_level > 0 .AND. BasisnD%para_SGType2%nb_SG > 10**5 ) THEN
   write(out_unitp,'(a)')              'OpPsi SG4 (%): [--0-10-20-30-40-50-60-70-80-90-100]'
   write(out_unitp,'(a)',ADVANCE='no') 'OpPsi SG4 (%): ['
   CALL flush_perso(out_unitp)
 END IF
 psi0Hkpsi0(:) = CZERO

 21 format(a,100(x,f8.6))

 DO iG=1,BasisnD%para_SGType2%nb_SG

   tab_l(:)  = BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval(:,iG)
   tab_nq(:) = getbis_tab_nq(tab_l,BasisnD%tab_basisPrimSG)
   tab_nb(:) = getbis_tab_nb(tab_l,BasisnD%tab_basisPrimSG)

   !write(6,*) 'iG',iG
   !transfert part of the psi%RvecB(:) to PsiRvec%V and psi0%RvecB(:) to Psi0Rvec%V
   ! the real and imaginary part are splited
   CALL tabPackedBasis_TO_tabR_AT_iG(PsiRvec%V, RPsi%RvecB, iG,BasisnD%para_SGType2)
   CALL tabPackedBasis_TO_tabR_AT_iG(PsiIvec%V, IPsi%RvecB, iG,BasisnD%para_SGType2)
   CALL tabPackedBasis_TO_tabR_AT_iG(Psi0Rvec%V,RPsi0%RvecB,iG,BasisnD%para_SGType2)
   CALL tabPackedBasis_TO_tabR_AT_iG(Psi0Ivec%V,IPsi0%RvecB,iG,BasisnD%para_SGType2)

   ! here, auto-c is done on the basis
   psi0Hkpsi0(0) = psi0Hkpsi0(0) + BasisnD%WeightSG(iG) *     cmplx(        &
     dot_product(Psi0Rvec%V,PsiRvec%V)+dot_product(Psi0Ivec%V,PsiIvec%V) ,  &
     dot_product(Psi0Rvec%V,PsiIvec%V)-dot_product(Psi0Ivec%V,PsiRvec%V), kind=Rkind)

   ! partial B to G
   CALL BDP_TO_GDP_OF_SmolyakRep(PsiRvec%V,BasisnD%tab_basisPrimSG,tab_l,tab_nq,tab_nb)
   CALL BDP_TO_GDP_OF_SmolyakRep(PsiIvec%V,BasisnD%tab_basisPrimSG,tab_l,tab_nq,tab_nb)

   CALL BDP_TO_GDP_OF_SmolyakRep(Psi0Rvec%V,BasisnD%tab_basisPrimSG,tab_l,tab_nq,tab_nb)
   CALL BDP_TO_GDP_OF_SmolyakRep(Psi0Ivec%V,BasisnD%tab_basisPrimSG,tab_l,tab_nq,tab_nb)

   nq = size(PsiRvec%V)
   CALL alloc_TypeRVec(Rw1,nq)
   CALL alloc_TypeRVec(Iw1,nq)
   CALL alloc_TypeRVec(Rw2(1),nq)
   CALL alloc_TypeRVec(Iw2(1),nq)
   CALL alloc_TypeRVec(Wrho,nq)


   !transfert part of the potential
   iterm00 = para_H%derive_term_TO_iterm(0,0)

   IF (para_H%OpGrid(iterm00)%grid_zero) THEN
     CALL alloc_NParray(V,(/ nq /),'V',name_sub)
     V(:) = ZERO
   ELSE
     IF (para_H%OpGrid(iterm00)%grid_cte) THEN
       CALL alloc_NParray(V,(/ nq /),'V',name_sub)
       V(:) = para_H%OpGrid(iterm00)%Mat_cte(1,1)
     ELSE
       CALL tabR2grid_TO_tabR1_AT_iG(V,para_H%OpGrid(iterm00)%Grid(:,1,1),&
                                     iG,BasisnD%para_SGType2)
     END IF
   END IF


   ! G calculation
   CALL alloc_NParray(sqRhoOVERJac,(/ nq /),      'sqRhoOVERJac',name_sub)
   CALL alloc_NParray(Jac,         (/ nq /),      'Jac',         name_sub)
   CALL alloc_NParray(GG,(/nq,mole%nb_act1,mole%nb_act1/),'GGiq',name_sub)
   CALL alloc_NParray(Qact,       (/mole%nb_var/),'Qact',        name_sub)
   CALL init_nDindexPrim(nDind_DPG,BasisnD%nb_basis,tab_nq,type_OF_nDindex=-1)

   DO iq=1,nq
     CALL get_Qact(Qact,mole%ActiveTransfo) ! rigid, flexible coordinates
     CALL Rec_Qact_SG4(Qact,BasisnD%tab_basisPrimSG,tab_l,nDind_DPG,iq,mole,err_sub)

     CALL get_d0GG(Qact,para_H%para_Tnum,mole,d0GG=GG(iq,:,:),          &
                                         def=.TRUE.,Jac=Jac(iq),Rho=Rho)
     sqRhoOVERJac(iq) = sqrt(Rho/Jac(iq))

   END DO


   DO iq=1,nq
     !Wrho%V(iq) = BasisnD%WeightSG(iG)
     Wrho%V(iq) = ONE
     CALL calc_nDindex(nDind_DPG,iq,tab_iq)
     DO ib=1,BasisnD%nb_basis
       L     = tab_l(ib)
       iq_ib = tab_iq(ib)
       Wrho%V(iq) = Wrho%V(iq) * Rec_WrhonD(BasisnD%tab_basisPrimSG(L,ib),iq_ib)
     END DO
   END DO
   CALL dealloc_nDindex(nDind_DPG)

   ! March_SG4
   Rw1%V(:)     = PsiRvec%V
   Iw1%V(:)     = PsiIvec%V

   rtj          = CONE

   DO j=1,para_propa%para_poly%npoly

     Rw2(1)%V(:)  = Rw1%V
     Iw2(1)%V(:)  = Iw1%V
     !CALL sub_TabOpPsi_OF_ONEGDP_FOR_SGtype4(Rw2,iG,para_H) ! in Rw2, we have H.Rw1
     !CALL sub_TabOpPsi_OF_ONEGDP_FOR_SGtype4(Iw2,iG,para_H) ! in Iw2, we have H.Iw1

     CALL sub_TabOpPsi_OF_ONEGDP_WithOp_FOR_SGtype4(Rw2,iG,para_H,      &
                                                    V,GG,sqRhoOVERJac,Jac)
     CALL sub_TabOpPsi_OF_ONEGDP_WithOp_FOR_SGtype4(Iw2,iG,para_H,      &
                                                    V,GG,sqRhoOVERJac,Jac)

        !write(6,21) 'Rw1',Rw1%V
        !write(6,21) 'Iw1',Iw1%V
        !write(6,21) 'Rw2',Rw2(1)%V
        !write(6,21) 'Iw2',Iw2(1)%V

     Rw2(1)%V(:) = Rw2(1)%V - para_H%E0 * Rw1%V ! equivalent sub_scaledOpPsi
     Iw2(1)%V(:) = Iw2(1)%V - para_H%E0 * Iw1%V ! equivalent sub_scaledOpPsi

        !write(6,21) 'Rw2',Rw2(1)%V
        !write(6,21) 'Iw2',Iw2(1)%V

     Rw1%V(:)    = Rw2(1)%V
     Iw1%V(:)    = Iw2(1)%V

     !CALL Overlap_psi1_psi2(psi0Hkpsi0(j),psi0,w1)
     !wrong to be adapted for the grid
     psi0Hkpsi0(j) = psi0Hkpsi0(j) + BasisnD%WeightSG(iG) * cmplx(      &
       dot_product(Psi0Rvec%V,Wrho%V*Rw1%V)+dot_product(Psi0Ivec%V,Wrho%V*Iw1%V) ,    &
       dot_product(Psi0Rvec%V,Wrho%V*Iw1%V)-dot_product(Psi0Ivec%V,Wrho%V*Rw1%V), kind=Rkind)

     rtj = rtj * cmplx(ZERO,-para_propa%WPdeltaT/real(j,kind=Rkind),kind=Rkind)

     Rw2(1)%V(:) = Rw1%V * real(rtj,kind=Rkind) - Iw1%V * Aimag(rtj)
     Iw2(1)%V(:) = Rw1%V * Aimag(rtj)           + Iw1%V * real(rtj,kind=Rkind)

     PsiRvec%V(:) = PsiRvec%V + Rw2(1)%V
     PsiIvec%V(:) = PsiIvec%V + Iw2(1)%V

     norm2_w2 = dot_product(Rw2(1)%V,Wrho%V*Rw2(1)%V) + dot_product(Iw2(1)%V,Wrho%V*Iw2(1)%V)


     IF (debug) write(out_unitp,*) 'iG,j,norm2_w2',iG,j,norm2_w2

     IF (norm2_w2 > TEN**15) THEN
       write(out_unitp,*) ' ERROR in ',name_sub
       write(out_unitp,*) ' Norm of the vector is TOO large (> 10^15)',j,    &
                                             norm2_w2
       write(out_unitp,*) ' => Reduce the time step !!'
       STOP
     END IF
     j_exit = j
     IF (norm2_w2 < para_propa%para_poly%poly_tol) EXIT


   END DO
   write(out_unitp,*) 'iG,j,norm2_w2',iG,j,norm2_w2

   CALL GDP_TO_BDP_OF_SmolyakRep(PsiRvec%V,BasisnD%tab_basisPrimSG,tab_l,tab_nq,tab_nb)
   CALL GDP_TO_BDP_OF_SmolyakRep(PsiIvec%V,BasisnD%tab_basisPrimSG,tab_l,tab_nq,tab_nb)


   !transfert back, PsiRvec%V to the psi%RvecB(:)
   CALL tabR_AT_iG_TO_tabPackedBasis(MarchRPsi%RvecB,PsiRvec%V,         &
                   iG,BasisnD%para_SGType2,BasisnD%WeightSG(iG))
   CALL tabR_AT_iG_TO_tabPackedBasis(MarchIPsi%RvecB,PsiIvec%V,         &
                   iG,BasisnD%para_SGType2,BasisnD%WeightSG(iG))

   !deallocate PsiRvec, Psi0Rvec
   CALL dealloc_TypeRVec(PsiRvec)
   CALL dealloc_TypeRVec(PsiIvec)
   CALL dealloc_TypeRVec(Psi0Rvec)
   CALL dealloc_TypeRVec(Psi0Ivec)
   CALL dealloc_TypeRVec(Rw1)
   CALL dealloc_TypeRVec(Iw1)
   CALL dealloc_TypeRVec(Rw2(1))
   CALL dealloc_TypeRVec(Iw2(1))
   CALL dealloc_TypeRVec(Wrho)
   CALL dealloc_NParray(V,           'V',           name_sub)
   CALL dealloc_NParray(sqRhoOVERJac,'sqRhoOVERJac',name_sub)
   CALL dealloc_NParray(Jac,         'Jac',         name_sub)
   CALL dealloc_NParray(GG,          'GG',          name_sub)
   CALL dealloc_NParray(Qact,        'Qact',        name_sub)


   IF (print_level > 0  .AND. BasisnD%para_SGType2%nb_SG > 10**5 .AND. &
       mod(iG,max(1,BasisnD%para_SGType2%nb_SG/10)) == 0) THEN
     write(out_unitp,'(a)',ADVANCE='no') '---'
     CALL flush_perso(out_unitp)
   END IF

   !write(6,*) 'iG done:',iG ; flush(6)
 END DO

 Psi%CvecB(:) = cmplx(MarchRPsi%RvecB,MarchIPsi%RvecB,kind=Rkind)


 CALL dealloc_psi(RPsi)
 CALL dealloc_psi(IPsi)
 CALL dealloc_psi(RPsi0)
 CALL dealloc_psi(IPsi0)
 CALL dealloc_psi(MarchRPsi)
 CALL dealloc_psi(MarchIPsi)


 CALL dealloc_NParray(tab_l ,'tab_l', name_sub)
 CALL dealloc_NParray(tab_nb,'tab_nb',name_sub)
 CALL dealloc_NParray(tab_nq,'tab_nq',name_sub)
 CALL dealloc_NParray(tab_iq,'tab_nq',name_sub)

 IF (print_level > 0 .AND. BasisnD%para_SGType2%nb_SG > 10**5) THEN
   write(out_unitp,'(a)',ADVANCE='yes') '----]'
 END IF
 CALL flush_perso(out_unitp)

 IF (abs(norm2_w2) > para_propa%para_poly%poly_tol) THEN
   write(out_unitp,*) ' ERROR in ',name_sub
   write(out_unitp,*) ' Norm of the last vector is TOO large',norm2_w2
   write(out_unitp,*) ' poly_tol: ',para_propa%para_poly%poly_tol
   write(out_unitp,*) ' => npoly or max_poly are TOO small',               &
                 para_propa%para_poly%npoly
   STOP
 END IF

 !- Phase Shift -----------------------------------
 phase = para_H%E0*para_propa%WPdeltaT
 psi = psi * exp(-cmplx(ZERO,phase,kind=Rkind))

 !- check norm ------------------
 CALL norme_psi(psi,GridRep=.FALSE.,BasisRep=.TRUE.,Renorm=.FALSE.)
 IF ( psi%norme > psi%max_norme) THEN
   T  = T + para_propa%WPdeltaT
   write(out_unitp,*) ' ERROR in ',name_sub
   write(out_unitp,*) ' STOP propagation: norm > max_norm',psi%norme
   para_propa%march_error   = .TRUE.
   para_propa%test_max_norm = .TRUE.
   STOP
 END IF



 microdeltaT = para_propa%WPdeltaT/                                &
               real(para_propa%nb_micro,kind=Rkind)
 microphase = phase/                                               &
               real(para_propa%nb_micro,kind=Rkind)

 phase = ZERO
 microT = ZERO

 DO it=1,para_propa%nb_micro

   microT = microT + microdeltaT
   phase = phase + microphase

   rtj = cmplx(ONE,ZERO,kind=Rkind)
   cdot = psi0Hkpsi0(0)
   DO j=1,j_exit

     rt_tmp =  cmplx(ZERO,-microT/real(j,kind=Rkind),kind=Rkind)
     rtj = rt_tmp * rtj
     cdot = cdot + psi0Hkpsi0(j) * rtj
   END DO
   cdot = cdot * cmplx( cos(phase),-sin(phase),kind=Rkind)
   CALL Write_AutoCorr(no,T+microT,cdot)
 END DO



!-----------------------------------------------------------
      IF (debug) THEN
        CALL norme_psi(psi)
        write(out_unitp,*) 'norm psi',psi%norme
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------


 END SUBROUTINE march_noD_SG4_GridRep_v0

 END MODULE mod_march_SG4

