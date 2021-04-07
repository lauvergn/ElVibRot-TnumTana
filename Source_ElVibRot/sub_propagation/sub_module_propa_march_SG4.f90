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
 MODULE mod_march_SG4
 USE mod_system
 USE mod_psi,    ONLY : param_psi
 USE mod_propa,  ONLY : param_propa,Calc_AutoCorr,Write_AutoCorr
 IMPLICIT NONE

 INTERFACE march_noD_SG4
   MODULE PROCEDURE march_noD_SG4_BasisRep_v2
 END INTERFACE

 PRIVATE
 PUBLIC :: march_noD_SG4



 CONTAINS


 SUBROUTINE march_noD_SG4_BasisRep_v2(T,no,psi,psi0,para_H,para_propa)
 USE mod_system
!$ USE omp_lib, only : OMP_GET_THREAD_NUM
 USE mod_nDindex
 USE mod_Coord_KEO,                ONLY : CoordType
 USE mod_basis_set_alloc,          ONLY : basis
 USE mod_basis_BtoG_GtoB_SGType4,  ONLY : tabPackedBasis_TO_tabR_AT_iG, &
                  tabR_AT_iG_TO_tabPackedBasis,TypeRVec,alloc_TypeRVec,dealloc_TypeRVec, &
                  BDP_TO_GDP_OF_SmolyakRep,GDP_TO_BDP_OF_SmolyakRep, &
                  getbis_tab_nq,getbis_tab_nb

 USE mod_psi,      ONLY : param_psi,copy_psi2TOpsi1,alloc_psi,   &
                          dealloc_psi,ecri_psi,Overlap_psi1_psi2,&
                          norm2_psi,ReNorm_psi

 USE mod_Op,              ONLY : param_Op,write_param_Op,sub_OpPsi
 USE mod_OpPsi_SG4,       ONLY : sub_TabOpPsi_OF_ONEDP_FOR_SGtype4
 USE mod_MPI_aux
 IMPLICIT NONE


!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_psi)   :: psi0,psi

!------ working variables ---------------------------------
      complex (kind=Rkind) :: E,cdot
      real (kind=Rkind)    :: T
      integer              :: no



 ! local variables
 TYPE (CoordType), pointer :: mole
 TYPE (basis),   pointer :: BasisnD

 TYPE (TypeRVec)    :: PsiRVec,PsiIVec,Psi0Rvec,Psi0Ivec
 TYPE (param_psi)   :: Rpsi,IPsi,Rpsi0,IPsi0,MarchRpsi,MarchIpsi,Oppsi

 integer :: ib,i,iG,nb_thread,itab,nq,nb,D,ith,err_sub

  integer,            allocatable       :: tab_l(:)

!----- for debuging --------------------------------------------------
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
  character (len=*), parameter :: name_sub='march_noD_SG4_BasisRep_v2'
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

  IF(openmpi) THEN
    CALL sub_OpPsi(psi,Oppsi,para_H)
    psi=Oppsi
  ELSE

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

  !RPsi0%RvecB(:) = Real(Psi0%CvecB(:),kind=Rkind)
  !IPsi0%RvecB(:) = Aimag(Psi0%CvecB(:))

 IF (print_level > 0 .AND. BasisnD%para_SGType2%nb_SG > 10**5 ) THEN
   write(out_unitp,'(a)')              'OpPsi SG4 (%): [--0-10-20-30-40-50-60-70-80-90-100]'
   write(out_unitp,'(a)',ADVANCE='no') 'OpPsi SG4 (%): ['
   CALL flush_perso(out_unitp)
 END IF

 !to be sure to have the correct number of threads, we use
 !   BasisnD%para_SGType2%nb_threads
 !$OMP parallel                                                  &
 !$OMP default(none)                                             &
 !$OMP shared(RPsi,IPsi,MarchRpsi,MarchIpsi)                     &
 !$OMP shared(T,para_H,BasisnD,para_propa,print_level,out_unitp) &
 !$OMP shared(MPI_id)                                            &
 !$OMP private(itab,iG,tab_l,ith,err_sub,PsiRvec,PsiIvec)        &
 !$OMP num_threads(BasisnD%para_SGType2%nb_threads)

 CALL alloc_NParray(tab_l ,(/ BasisnD%para_SGType2%nDind_SmolyakRep%ndim /),'tab_l', name_sub)

 !--------------------------------------------------------------
 !-- For the initialization of tab_l(:) and the use of ADD_ONE_TO_nDindex in the parallel loop
 ith = 0
 !$ ith = OMP_GET_THREAD_NUM()
 tab_l(:) = BasisnD%para_SGType2%nDval_init(:,ith+1)
 !--------------------------------------------------------------

 ! we are not using the parallel do, to be able to use the correct initialized tab_l with nDval_init
 DO iG=BasisnD%para_SGType2%iG_th(ith+1),BasisnD%para_SGType2%fG_th(ith+1)

   CALL ADD_ONE_TO_nDindex(BasisnD%para_SGType2%nDind_SmolyakRep,tab_l,iG=iG,err_sub=err_sub)

   !write(out_unitp,*) 'iG',iG
   !transfert part of the psi%RvecB(:) to PsiRvec%R and psi0%RvecB(:) to Psi0Rvec%R
   ! the real and imaginary part are splited

   CALL tabPackedBasis_TO_tabR_AT_iG(PsiRvec%V, RPsi%RvecB, iG,BasisnD%para_SGType2)
   CALL tabPackedBasis_TO_tabR_AT_iG(PsiIvec%V, IPsi%RvecB, iG,BasisnD%para_SGType2)
   !CALL tabPackedBasis_TO_tabR_AT_iG(Psi0Rvec%V,RPsi0%RvecB,iG,BasisnD%para_SGType2)
   !CALL tabPackedBasis_TO_tabR_AT_iG(Psi0Ivec%V,IPsi0%RvecB,iG,BasisnD%para_SGType2)

   CALL march_noD_ONE_DP_SG4(T,PsiRVec,PsiIVec,iG,tab_l,para_H,para_propa)

   CALL tabR_AT_iG_TO_tabPackedBasis(MarchRPsi%RvecB,PsiRvec%V,         &
                           iG,BasisnD%para_SGType2,BasisnD%WeightSG(iG))
   CALL tabR_AT_iG_TO_tabPackedBasis(MarchIPsi%RvecB,PsiIvec%V,         &
                           iG,BasisnD%para_SGType2,BasisnD%WeightSG(iG))

   !deallocate PsiRvec, Psi0Rvec
   CALL dealloc_TypeRVec(PsiRvec)
   CALL dealloc_TypeRVec(PsiIvec)
   !CALL dealloc_TypeRVec(Psi0Rvec)
   !CALL dealloc_TypeRVec(Psi0Ivec)

   IF (print_level > 0  .AND. BasisnD%para_SGType2%nb_SG > 10**5 .AND. &
       mod(iG,max(1,BasisnD%para_SGType2%nb_SG/10)) == 0 .AND. MPI_id==0) THEN
     write(out_unitp,'(a)',ADVANCE='no') '---'
     CALL flush_perso(out_unitp)
   END IF

   !write(out_unitp,*) 'iG done:',iG ; flush(out_unitp)
 END DO
 CALL dealloc_NParray(tab_l,'tabl_l',name_sub)
 !$OMP   END PARALLEL

 Psi%CvecB(:) = cmplx(MarchRPsi%RvecB,MarchIPsi%RvecB,kind=Rkind)

 CALL dealloc_psi(RPsi)
 CALL dealloc_psi(IPsi)
 !CALL dealloc_psi(RPsi0)
 !CALL dealloc_psi(IPsi0)
 CALL dealloc_psi(MarchRPsi)
 CALL dealloc_psi(MarchIPsi)

  ENDIF ! for if(openmpi)

 IF (print_level > 0 .AND. BasisnD%para_SGType2%nb_SG > 10**5) THEN
   IF(MPI_id==0) write(out_unitp,'(a)',ADVANCE='yes') '----]'
 END IF
 CALL flush_perso(out_unitp)

 !- check norm ------------------
 IF(keep_MPI) CALL norm2_psi(psi,GridRep=.FALSE.,BasisRep=.TRUE.)
 IF(openmpi .AND. MPI_scheme/=1)  CALL MPI_Bcast_(psi%norm2,size1_MPI,root_MPI)

 IF (debug) write(out_unitp,*) 'norm^2',psi%norm2

 IF ( psi%norm2 > para_propa%max_norm2) THEN
   T  = T + para_propa%WPdeltaT
   write(out_unitp,*) ' ERROR in ',name_sub
   write(out_unitp,*) ' STOP propagation: norm > max_norm',psi%norm2
   para_propa%march_error   = .TRUE.
   para_propa%test_max_norm = .TRUE.
   STOP
 END IF
 IF(keep_MPI) CALL Renorm_psi(psi)
 IF(openmpi .AND. MPI_scheme/=1) CALL MPI_Bcast_(psi%norm2,size1_MPI,root_MPI)

 IF(MPI_id==0) THEN
   CALL Overlap_psi1_psi2(cdot,psi0,psi)
   CALL Write_AutoCorr(no,T+para_propa%WPdeltaT,cdot)
 ENDIF

!-----------------------------------------------------------
 IF (debug) THEN
   CALL norm2_psi(psi)
   write(out_unitp,*) 'norm psi',psi%norm2
   write(out_unitp,*) 'END ',name_sub
 END IF
!-----------------------------------------------------------

 END SUBROUTINE march_noD_SG4_BasisRep_v2

 SUBROUTINE march_noD_ONE_DP_SG4(T,PsiRVec,PsiIVec,iG,tab_l,para_H,para_propa)
 USE mod_system
 USE mod_basis_set_alloc,          ONLY : basis
 USE mod_basis_BtoG_GtoB_SGType4,  ONLY : TypeRVec,alloc_TypeRVec,dealloc_TypeRVec
 USE mod_OpPsi_SG4,                ONLY : sub_TabOpPsi_OF_ONEDP_FOR_SGtype4
 USE mod_Op,                       ONLY : param_Op,write_param_Op
 IMPLICIT NONE


 TYPE (TypeRVec),      intent(inout) :: PsiRVec,PsiIVec
 real (kind=Rkind),    intent(in)    :: T
 integer,              intent(in)    :: iG
 integer,              intent(in)    :: tab_l(:)

!----- variables pour la namelist minimum ----------------------------
 TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ----------------------------
 TYPE (param_propa) :: para_propa


!------ working variables ---------------------------------
 real (kind=Rkind)    :: phase
 complex (kind=Rkind) :: rtj
 integer              :: nb,j,j_exit
 real (kind=Rkind)    :: E0,norm2_w2,norm2_AT_iG
 TYPE (TypeRVec)      :: RIw2(2),RIw1(2)


!----- for debuging --------------------------------------------------
 logical, parameter :: debug=.FALSE.
 !logical, parameter :: debug=.TRUE.
 character (len=*), parameter :: name_sub='march_noD_ONE_DP_SG4'
!-----------------------------------------------------------
 IF (debug) THEN
   write(out_unitp,*) 'BEGINNING ',name_sub
   write(out_unitp,*) 'T',T
   write(out_unitp,*) 'deltaT',para_propa%WPdeltaT
   write(out_unitp,*) 'nOD',para_propa%para_poly%npoly
 END IF
!-----------------------------------------------------------

   nb = size(PsiRvec%V)

   CALL alloc_TypeRVec(RIw1(1),nb)
   CALL alloc_TypeRVec(RIw1(2),nb)
   CALL alloc_TypeRVec(RIw2(1),nb)
   CALL alloc_TypeRVec(RIw2(2),nb)


   ! March_SG4

   RIw1(1)%V(:)     = PsiRvec%V
   RIw1(2)%V(:)     = PsiIvec%V
   norm2_AT_iG = dot_product(PsiRvec%V,PsiRvec%V) + dot_product(PsiIvec%V,PsiIvec%V)
   IF (debug) write(out_unitp,*) 'norm2_AT_iG',norm2_AT_iG


   ! initial energy of the given DP
   RIw2(1)%V(:)  = RIw1(1)%V
   RIw2(2)%V(:)  = RIw1(2)%V

   CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4(RIw2,iG,tab_l,para_H) ! in RIw2, we have H.RIw1(1) and H.RIw1(2)

   E0 = (dot_product(RIw1(1)%V,RIw2(1)%V)+dot_product(RIw1(2)%V,RIw2(2)%V))/norm2_AT_iG

   rtj          = CONE

   DO j=1,para_propa%para_poly%npoly


     RIw2(1)%V(:)  = RIw1(1)%V
     RIw2(2)%V(:)  = RIw1(2)%V

     CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4(RIw2,iG,tab_l,para_H) ! in RIw2, we have H.RIw1(1) and H.RIw1(2)

     RIw2(1)%V(:) = RIw2(1)%V - E0 * RIw1(1)%V ! equivalent sub_scaledOpPsi
     RIw2(2)%V(:) = RIw2(2)%V - E0 * RIw1(2)%V ! equivalent sub_scaledOpPsi

     RIw1(1)%V(:)    = RIw2(1)%V
     RIw1(2)%V(:)    = RIw2(2)%V


     rtj = rtj * cmplx(ZERO,-para_propa%WPdeltaT/real(j,kind=Rkind),kind=Rkind)

     RIw2(1)%V(:) = RIw1(1)%V * real(rtj,kind=Rkind) - RIw1(2)%V * Aimag(rtj)
     RIw2(2)%V(:) = RIw1(1)%V * Aimag(rtj)           + RIw1(2)%V * real(rtj,kind=Rkind)


     PsiRvec%V(:) = PsiRvec%V + RIw2(1)%V
     PsiIvec%V(:) = PsiIvec%V + RIw2(2)%V

     norm2_w2 = (dot_product(RIw2(1)%V,RIw2(1)%V) + dot_product(RIw2(2)%V,RIw2(2)%V)) / norm2_AT_iG


     IF (debug) write(out_unitp,*) 'iG,j,norm2_w2/norm2_AT_iG',iG,j,norm2_w2

     IF (norm2_w2 > TEN**15) THEN
       write(out_unitp,*) ' ERROR in ',name_sub
       write(out_unitp,*) ' Norm^2 of the vector is TOO large (> 10^15)',j,norm2_w2
       write(out_unitp,*) ' => Reduce the time step !!'
       STOP
     END IF
     j_exit = j
     IF (norm2_w2 < para_propa%para_poly%poly_tol) EXIT


   END DO

   !- Phase Shift -----------------------------------
   phase = E0*para_propa%WPdeltaT
   RIw1(1)%V(:) =  PsiRvec%V*cos(phase) + PsiIvec%V * sin(phase)
   RIw1(2)%V(:) = -PsiRvec%V*sin(phase) + PsiIvec%V * cos(phase)

   PsiRvec%V = RIw1(1)%V
   PsiIvec%V = RIw1(2)%V

   CALL dealloc_TypeRVec(RIw1(1))
   CALL dealloc_TypeRVec(RIw1(2))
   CALL dealloc_TypeRVec(RIw2(1))
   CALL dealloc_TypeRVec(RIw2(2))

 !-----------------------------------------------------------
 IF (debug) THEN
   write(out_unitp,*) 'iG,j_exit,norm2_w2',iG,j_exit,norm2_w2
   write(out_unitp,*) 'END ',name_sub
 END IF
 !-----------------------------------------------------------

END SUBROUTINE march_noD_ONE_DP_SG4
 SUBROUTINE march_noD_SG4_BasisRep_v1(T,no,psi,psi0,para_H,para_propa)
 USE mod_system
!$ USE omp_lib, only : OMP_GET_THREAD_NUM
 USE mod_nDindex
 USE mod_Coord_KEO,                ONLY : CoordType
 USE mod_basis_set_alloc,          ONLY : basis
 USE mod_basis_BtoG_GtoB_SGType4,  ONLY : tabPackedBasis_TO_tabR_AT_iG, &
                  tabR_AT_iG_TO_tabPackedBasis,TypeRVec,alloc_TypeRVec,dealloc_TypeRVec, &
                  BDP_TO_GDP_OF_SmolyakRep,GDP_TO_BDP_OF_SmolyakRep, &
                  getbis_tab_nq,getbis_tab_nb

 USE mod_psi,             ONLY : param_psi,copy_psi2TOpsi1,alloc_psi,   &
                                 dealloc_psi,ecri_psi,Overlap_psi1_psi2,&
                                 norm2_psi,ReNorm_psi

 USE mod_Op,              ONLY : param_Op,write_param_Op
 USE mod_OpPsi_SG4,       ONLY : sub_TabOpPsi_OF_ONEDP_FOR_SGtype4

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
      real (kind=Rkind)    :: E0,rg,norm2_w2,norm2_AT_iG


 ! local variables
 TYPE (CoordType), pointer :: mole
 TYPE (basis),   pointer :: BasisnD

 TYPE (TypeRVec)    :: PsiRVec,PsiIVec,Psi0Rvec,Psi0Ivec,Rw1,Iw1,RIw2(2)
 TYPE (param_psi)   :: Rpsi,IPsi,Rpsi0,IPsi0,MarchRpsi,MarchIpsi

 integer :: ib,i,iG,nb_thread,itab,nq,nb,D,ith,err_sub

  integer,            allocatable       :: tab_nq(:)
  integer,            allocatable       :: tab_nb(:)
  integer,            allocatable       :: tab_l(:)

!----- for debuging --------------------------------------------------
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
  character (len=*), parameter :: name_sub='march_noD_SG4_BasisRep_v1'
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

   !write(out_unitp,*) 'iG',iG
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
   CALL alloc_TypeRVec(RIw2(1),nb)
   CALL alloc_TypeRVec(RIw2(2),nb)


   ! March_SG4

   Rw1%V(:)     = PsiRvec%V
   Iw1%V(:)     = PsiIvec%V
   norm2_AT_iG = dot_product(PsiRvec%V,PsiRvec%V) + dot_product(PsiIvec%V,PsiIvec%V)
   IF (debug) write(out_unitp,*) 'norm2_AT_iG',norm2_AT_iG


   ! initial energy of the given DP
   RIw2(1)%V(:)  = Rw1%V
   RIw2(2)%V(:)  = Iw1%V

   CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4(RIw2,iG,tab_l,para_H) ! in RIw2, we have H.Rw1 and H.Iw1

   E0 = (dot_product(Rw1%V,RIw2(1)%V)+dot_product(Iw1%V,RIw2(2)%V))/norm2_AT_iG

   rtj          = CONE

   DO j=1,para_propa%para_poly%npoly


     RIw2(1)%V(:)  = Rw1%V
     RIw2(2)%V(:)  = Iw1%V

     CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4(RIw2,iG,tab_l,para_H) ! in RIw2, we have H.Rw1 and H.Iw1

     RIw2(1)%V(:) = RIw2(1)%V - E0 * Rw1%V ! equivalent sub_scaledOpPsi
     RIw2(2)%V(:) = RIw2(2)%V - E0 * Iw1%V ! equivalent sub_scaledOpPsi

     Rw1%V(:)    = RIw2(1)%V
     Iw1%V(:)    = RIw2(2)%V

     !CALL Overlap_psi1_psi2(psi0Hkpsi0(j),psi0,w1)
     psi0Hkpsi0(j) = psi0Hkpsi0(j) + BasisnD%WeightSG(iG) * cmplx(      &
       dot_product(Psi0Rvec%V,Rw1%V)+dot_product(Psi0Ivec%V,Iw1%V) ,    &
       dot_product(Psi0Rvec%V,Iw1%V)-dot_product(Psi0Ivec%V,Rw1%V), kind=Rkind)


     rtj = rtj * cmplx(ZERO,-para_propa%WPdeltaT/real(j,kind=Rkind),kind=Rkind)

     RIw2(1)%V(:) = Rw1%V * real(rtj,kind=Rkind) - Iw1%V * Aimag(rtj)
     RIw2(2)%V(:) = Rw1%V * Aimag(rtj)           + Iw1%V * real(rtj,kind=Rkind)


     PsiRvec%V(:) = PsiRvec%V + RIw2(1)%V
     PsiIvec%V(:) = PsiIvec%V + RIw2(2)%V

     norm2_w2 = (dot_product(RIw2(1)%V,RIw2(1)%V) + dot_product(RIw2(2)%V,RIw2(2)%V)) / norm2_AT_iG


     IF (debug) write(out_unitp,*) 'iG,j,norm2_w2/norm2_AT_iG',iG,j,norm2_w2

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

   !- Phase Shift -----------------------------------
   phase = E0*para_propa%WPdeltaT
   RIw2(1)%V(:) =  PsiRvec%V*cos(phase) + PsiIvec%V * sin(phase)
   RIw2(2)%V(:) = -PsiRvec%V*sin(phase) + PsiIvec%V * cos(phase)

   CALL tabR_AT_iG_TO_tabPackedBasis(MarchRPsi%RvecB,RIw2(1)%V,         &
                   iG,BasisnD%para_SGType2,BasisnD%WeightSG(iG))
   CALL tabR_AT_iG_TO_tabPackedBasis(MarchIPsi%RvecB,RIw2(2)%V,         &
                   iG,BasisnD%para_SGType2,BasisnD%WeightSG(iG))


   IF (debug)  write(out_unitp,*) 'iG,j,norm2_w2',iG,j,norm2_w2
   !deallocate PsiRvec, Psi0Rvec
   CALL dealloc_TypeRVec(PsiRvec)
   CALL dealloc_TypeRVec(PsiIvec)
   CALL dealloc_TypeRVec(Psi0Rvec)
   CALL dealloc_TypeRVec(Psi0Ivec)
   CALL dealloc_TypeRVec(Rw1)
   CALL dealloc_TypeRVec(Iw1)
   CALL dealloc_TypeRVec(RIw2(1))
   CALL dealloc_TypeRVec(RIw2(2))

   IF (print_level > 0  .AND. BasisnD%para_SGType2%nb_SG > 10**5 .AND. &
       mod(iG,max(1,BasisnD%para_SGType2%nb_SG/10)) == 0 .AND. MPI_id==0) THEN
     write(out_unitp,'(a)',ADVANCE='no') '---'
     CALL flush_perso(out_unitp)
   END IF

   !write(out_unitp,*) 'iG done:',iG ; flush(out_unitp)
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
   IF(MPI_id==0) write(out_unitp,'(a)',ADVANCE='yes') '----]'
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

 !- check norm ------------------
 CALL norm2_psi(psi,GridRep=.FALSE.,BasisRep=.TRUE.)
 IF (debug) write(out_unitp,*) 'norm^2',psi%norm2

 IF ( psi%norm2 > para_propa%max_norm2) THEN
   T  = T + para_propa%WPdeltaT
   write(out_unitp,*) ' ERROR in ',name_sub
   write(out_unitp,*) ' STOP propagation: norm > max_norm',psi%norm2
   para_propa%march_error   = .TRUE.
   para_propa%test_max_norm = .TRUE.
   STOP
 END IF
 !CALL Renorm_psi(psi)



 CALL Overlap_psi1_psi2(cdot,psi0,psi)
 CALL Write_AutoCorr(no,T+para_propa%WPdeltaT,cdot)



!-----------------------------------------------------------
      IF (debug) THEN
        CALL norm2_psi(psi)
        write(out_unitp,*) 'norm psi',psi%norm2
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

 END SUBROUTINE march_noD_SG4_BasisRep_v1
 SUBROUTINE march_noD_SG4_BasisRep_v0(T,no,psi,psi0,para_H,para_propa)
 USE mod_system
!$ USE omp_lib, only : OMP_GET_THREAD_NUM
 USE mod_nDindex
 USE mod_Coord_KEO,                ONLY : CoordType
 USE mod_basis_set_alloc,          ONLY : basis
 USE mod_basis_BtoG_GtoB_SGType4,  ONLY : tabPackedBasis_TO_tabR_AT_iG, &
                  tabR_AT_iG_TO_tabPackedBasis,TypeRVec,alloc_TypeRVec,dealloc_TypeRVec, &
                  BDP_TO_GDP_OF_SmolyakRep,GDP_TO_BDP_OF_SmolyakRep, &
                  getbis_tab_nq,getbis_tab_nb

 USE mod_psi,             ONLY : param_psi,copy_psi2TOpsi1,alloc_psi,   &
                                 dealloc_psi,ecri_psi,Overlap_psi1_psi2,&
                                 norm2_psi,ReNorm_psi

 USE mod_Op,              ONLY : param_Op,write_param_Op
 USE mod_OpPsi_SG4,       ONLY : sub_TabOpPsi_OF_ONEDP_FOR_SGtype4
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
      real (kind=Rkind)    :: E0,rg,norm2_w2,norm2_AT_iG


 ! local variables
 TYPE (CoordType), pointer :: mole
 TYPE (basis),   pointer :: BasisnD

 TYPE (TypeRVec)    :: PsiRVec,PsiIVec,Psi0Rvec,Psi0Ivec,Rw1,Rw2(1),Iw1,Iw2(1)
 TYPE (param_psi)   :: Rpsi,IPsi,Rpsi0,IPsi0,MarchRpsi,MarchIpsi

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

   !write(out_unitp,*) 'iG',iG
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
   norm2_AT_iG = dot_product(PsiRvec%V,PsiRvec%V) + dot_product(PsiIvec%V,PsiIvec%V)
   IF (debug) write(out_unitp,*) 'norm2_AT_iG',norm2_AT_iG


   ! initial energy of the given DP
   Rw2(1)%V(:)  = Rw1%V
   Iw2(1)%V(:)  = Iw1%V

   CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4(Rw2,iG,tab_l,para_H) ! in w2, we have H.Rw1
   CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4(Iw2,iG,tab_l,para_H) ! in w2, we have H.Iw1

   E0 = (dot_product(Rw1%V,Rw2(1)%V)+dot_product(Iw1%V,Iw2(1)%V))/norm2_AT_iG

       !write(out_unitp,21) 'Rw1',Rw1%V
       !write(out_unitp,21) 'Iw1',Iw1%V

   rtj          = CONE

   DO j=1,para_propa%para_poly%npoly


     Rw2(1)%V(:)  = Rw1%V
     Iw2(1)%V(:)  = Iw1%V

     CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4(Rw2,iG,tab_l,para_H) ! in w2, we have H.Rw1
     CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4(Iw2,iG,tab_l,para_H) ! in w2, we have H.Iw1

        !write(out_unitp,21) 'Rw1',Rw1%V
        !write(out_unitp,21) 'Iw1',Iw1%V
        !write(out_unitp,21) 'Rw2',Rw2(1)%V
        !write(out_unitp,21) 'Iw2',Iw2(1)%V

     Rw2(1)%V(:) = Rw2(1)%V - E0 * Rw1%V ! equivalent sub_scaledOpPsi
     Iw2(1)%V(:) = Iw2(1)%V - E0 * Iw1%V ! equivalent sub_scaledOpPsi

        !write(out_unitp,21) 'Rw2',Rw2(1)%V
        !write(out_unitp,21) 'Iw2',Iw2(1)%V

     Rw1%V(:)    = Rw2(1)%V
     Iw1%V(:)    = Iw2(1)%V

     !CALL Overlap_psi1_psi2(psi0Hkpsi0(j),psi0,w1)
     psi0Hkpsi0(j) = psi0Hkpsi0(j) + BasisnD%WeightSG(iG) * cmplx(      &
       dot_product(Psi0Rvec%V,Rw1%V)+dot_product(Psi0Ivec%V,Iw1%V) ,    &
       dot_product(Psi0Rvec%V,Iw1%V)-dot_product(Psi0Ivec%V,Rw1%V), kind=Rkind)


     rtj = rtj * cmplx(ZERO,-para_propa%WPdeltaT/real(j,kind=Rkind),kind=Rkind)

     Rw2(1)%V(:) = Rw1%V * real(rtj,kind=Rkind) - Iw1%V * Aimag(rtj)
     Iw2(1)%V(:) = Rw1%V * Aimag(rtj)           + Iw1%V * real(rtj,kind=Rkind)

        !write(out_unitp,21) 'Rw2*rtj',Rw2(1)%V
        !write(out_unitp,21) 'Iw2*rtj',Iw2(1)%V

     PsiRvec%V(:) = PsiRvec%V + Rw2(1)%V
     PsiIvec%V(:) = PsiIvec%V + Iw2(1)%V

        !write(out_unitp,21) 'Rpsi',PsiRvec%V
        !write(out_unitp,21) 'Ipsi',PsiIvec%V

     norm2_w2 = (dot_product(Rw2(1)%V,Rw2(1)%V) + dot_product(Iw2(1)%V,Iw2(1)%V)) / norm2_AT_iG


     IF (debug) write(out_unitp,*) 'iG,j,norm2_w2/norm2_AT_iG',iG,j,norm2_w2

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

   !- Phase Shift -----------------------------------
   phase = E0*para_propa%WPdeltaT
   Rw2(1)%V(:) =  PsiRvec%V*cos(phase) + PsiIvec%V * sin(phase)
   Iw2(1)%V(:) = -PsiRvec%V*sin(phase) + PsiIvec%V * cos(phase)

   !Rw2(1)%V(:) =  PsiRvec%V
   !Iw2(1)%V(:) =  PsiIvec%V

   CALL tabR_AT_iG_TO_tabPackedBasis(MarchRPsi%RvecB,Rw2(1)%V,         &
                   iG,BasisnD%para_SGType2,BasisnD%WeightSG(iG))
   CALL tabR_AT_iG_TO_tabPackedBasis(MarchIPsi%RvecB,Iw2(1)%V,         &
                   iG,BasisnD%para_SGType2,BasisnD%WeightSG(iG))


   IF (debug)  write(out_unitp,*) 'iG,j,norm2_w2',iG,j,norm2_w2
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
       mod(iG,max(1,BasisnD%para_SGType2%nb_SG/10)) == 0 .AND. MPI_id==0) THEN
     write(out_unitp,'(a)',ADVANCE='no') '---'
     CALL flush_perso(out_unitp)
   END IF

   !write(out_unitp,*) 'iG done:',iG ; flush(out_unitp)
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
   IF(MPI_id==0) write(out_unitp,'(a)',ADVANCE='yes') '----]'
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
 !phase = para_H%E0*para_propa%WPdeltaT
 !psi = psi * exp(-cmplx(ZERO,phase,kind=Rkind))

 !write(out_unitp,*) ' Psi after phase shift '
 !CALL ecri_psi(psi=psi)


 !- check norm ------------------
 CALL norm2_psi(psi,GridRep=.FALSE.,BasisRep=.TRUE.)
 IF (debug) write(out_unitp,*) 'norm^2',psi%norm2

 IF ( psi%norm2 > para_propa%max_norm2) THEN
   T  = T + para_propa%WPdeltaT
   write(out_unitp,*) ' ERROR in ',name_sub
   write(out_unitp,*) ' STOP propagation: norm > max_norm',psi%norm2
   para_propa%march_error   = .TRUE.
   para_propa%test_max_norm = .TRUE.
   STOP
 END IF
 CALL Renorm_psi(psi)



 CALL Overlap_psi1_psi2(cdot,psi0,psi)
 CALL Write_AutoCorr(no,T+para_propa%WPdeltaT,cdot)



!-----------------------------------------------------------
      IF (debug) THEN
        CALL norm2_psi(psi)
        write(out_unitp,*) 'norm psi',psi%norm2
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

 END SUBROUTINE march_noD_SG4_BasisRep_v0

 END MODULE mod_march_SG4

