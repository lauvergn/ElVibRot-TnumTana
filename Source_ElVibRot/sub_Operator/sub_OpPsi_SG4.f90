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


MODULE mod_OpPsi_SG4

!PRIVATE
!PUBLIC :: sub_TabOpPsi_OF_ONEDP_FOR_SGtype4,sub_TabOpPsi_OF_ONEGDP_FOR_SGtype4
!PUBLIC :: sub_TabOpPsi_OF_ONEGDP_WithOp_FOR_SGtype4

CONTAINS

 SUBROUTINE sub_OpPsi_FOR_SGtype4(Psi,OpPsi,para_Op)

 USE mod_system
!$ USE omp_lib, only : OMP_GET_THREAD_NUM
 USE mod_nDindex
 USE mod_Coord_KEO,                ONLY : zmatrix

 USE mod_basis_set_alloc,          ONLY : basis
 USE mod_basis_RCVec_SGType4,      ONLY : TypeRVec,dealloc_TypeRVec
 USE mod_basis_BtoG_GtoB_SGType4,  ONLY : tabPackedBasis_TO_tabR_AT_iG, &
                                          tabR_AT_iG_TO_tabPackedBasis

 USE mod_psi_set_alloc,            ONLY : param_psi,ecri_psi,assignment (=)

 USE mod_Op,                       ONLY : param_Op,write_param_Op
 IMPLICIT NONE

 TYPE (param_psi), intent(in)      :: Psi
 TYPE (param_psi), intent(inout)   :: OpPsi

 TYPE (param_Op),  intent(inout)   :: para_Op

 ! local variables
 TYPE (zmatrix), pointer :: mole
 TYPE (basis),   pointer :: BasisnD

 TYPE (TypeRVec)         :: PsiR

 integer       :: ib,i,iG,iiG,nb_thread,ith,iterm00
 integer       :: tab_l(para_Op%BasisnD%nb_basis)
 logical       :: not_init

 !----- for debuging ----------------------------------------------
 character (len=*), parameter :: name_sub='sub_OpPsi_FOR_SGtype4'
 logical, parameter :: debug = .FALSE.
 !logical, parameter :: debug = .TRUE.
 !-----------------------------------------------------------------
 IF (debug) THEN
   write(out_unitp,*) 'BEGINNING ',name_sub
   write(out_unitp,*) 'nb_bie,nb_baie',para_Op%nb_bie,para_Op%nb_baie
   write(out_unitp,*) 'nb_act1',para_Op%mole%nb_act1
   write(out_unitp,*) 'nb_var',para_Op%mole%nb_var
   !CALL flush_perso(out_unitp)
   !CALL write_param_Op(para_Op)
   write(out_unitp,*)
   write(out_unitp,*) 'PsiBasisRep'
   write(out_unitp,*) 'Psi%RvecB',Psi%RvecB
   CALL flush_perso(out_unitp)
 END IF
 !-----------------------------------------------------------------
  !CALL Check_mem()

  mole    => para_Op%mole
  BasisnD => para_Op%BasisnD

  IF (Psi%cplx) STOP 'cplx'
  IF (para_Op%nb_bie /= 1) STOP 'nb_bie /= 1'

  IF (OpPsi_omp == 0) THEN
    nb_thread = 1
  ELSE
    nb_thread = OpPsi_maxth
  END IF


 nb_mult_OpPsi   = 0
 OpPsi           = Psi ! for the allocation. It has to be changed!
 OpPsi%RvecB(:)  = ZERO

 IF (print_level > 0 .AND. BasisnD%para_SGType2%nb_SG > 10**5 ) THEN
   write(out_unitp,'(a)')              'OpPsi SG4 (%): [--0-10-20-30-40-50-60-70-80-90-100]'
   write(out_unitp,'(a)',ADVANCE='no') 'OpPsi SG4 (%): ['
   CALL flush_perso(out_unitp)
 END IF

 IF (OpPsi_omp == 0) THEN
   DO iG=1,BasisnD%para_SGType2%nb_SG

     !write(6,*) 'iG',iG
     !transfert part of the psi%RvecB(:) to PsiR%V
     CALL tabPackedBasis_TO_tabR_AT_iG(PsiR%V,psi%RvecB,    &
                                         iG,BasisnD%para_SGType2)
     !write(6,*) 'iG,PsiR',iG,PsiR%V

     CALL sub_OpPsi_OF_ONEDP_FOR_SGtype4_ompGrid(PsiR,iG,para_Op)

     !transfert back, PsiR%V (HPsi) to the psi%RvecB(:)
     CALL tabR_AT_iG_TO_tabPackedBasis(OpPsi%RvecB,PsiR%V,  &
                     iG,BasisnD%para_SGType2,BasisnD%WeightSG(iG))

     !deallocate PsiR
     CALL dealloc_TypeRVec(PsiR)


     IF (print_level > 0  .AND. BasisnD%para_SGType2%nb_SG > 10**5 .AND. &
         mod(iG,max(1,BasisnD%para_SGType2%nb_SG/10)) == 0) THEN
       write(out_unitp,'(a)',ADVANCE='no') '---'
       CALL flush_perso(out_unitp)
     END IF

     !write(6,*) 'iG done:',iG ; flush(6)
   END DO
 ELSE

   !to be sure to have the correct number of threads
   nb_thread = BasisnD%para_SGType2%nb_threads

   !$OMP parallel                                                &
   !$OMP default(none)                                           &
   !$OMP shared(Psi,OpPsi)                                       &
   !$OMP shared(para_Op,BasisnD,print_level,out_unitp)           &
   !$OMP private(iG,iiG,tab_l,PsiR,ith)                          &
   !$OMP num_threads(BasisnD%para_SGType2%nb_threads)

   !--------------------------------------------------------------
   !-- For the initialization of tab_l(:) and the use of ADD_ONE_TO_nDindex in the parallel loop
   ith = 0
   !$ ith = OMP_GET_THREAD_NUM()
   tab_l(:) = BasisnD%para_SGType2%nDval_init(:,ith+1)
   !--------------------------------------------------------------

   ! we are not using the parallel do, to be able to use the correct initialized tab_l with nDval_init
   DO iG=BasisnD%para_SGType2%iG_th(ith+1),BasisnD%para_SGType2%fG_th(ith+1)

     CALL ADD_ONE_TO_nDindex(BasisnD%para_SGType2%nDind_SmolyakRep,tab_l,iG=iG)

     !write(6,*) 'iG',iG ; flush(6)
     !transfert part of the psi%RvecB(:) to PsiR%V
     CALL tabPackedBasis_TO_tabR_AT_iG(PsiR%V,psi%RvecB,iG,BasisnD%para_SGType2)
     !write(6,*) 'iG,PsiR',iG,PsiR%V ; flush(6)

     CALL sub_OpPsi_OF_ONEDP_FOR_SGtype4(PsiR,iG,tab_l,para_Op)

     !transfert back, PsiR%V (HPsi) to the psi%RvecB(:)
     CALL tabR_AT_iG_TO_tabPackedBasis(OpPsi%RvecB,PsiR%V,              &
                           iG,BasisnD%para_SGType2,BasisnD%WeightSG(iG))

     !deallocate PsiR(:)
     CALL dealloc_TypeRVec(PsiR)


     IF (print_level > 0  .AND. BasisnD%para_SGType2%nb_SG > 10**5 .AND. &
         mod(iG,max(1,BasisnD%para_SGType2%nb_SG/10)) == 0) THEN
       write(out_unitp,'(a)',ADVANCE='no') '---'
       CALL flush_perso(out_unitp)
     END IF

     !write(6,*) 'iG done:',iG ; flush(6)
   END DO
  !$OMP   END PARALLEL
END IF

 IF (print_level > 0 .AND. BasisnD%para_SGType2%nb_SG > 10**5) THEN
   write(out_unitp,'(a)',ADVANCE='yes') '----]'
 END IF
 CALL flush_perso(out_unitp)

 iterm00 = para_Op%derive_term_TO_iterm(0,0)
 IF (associated(para_Op%OpGrid)) THEN
   para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_done =            &
                     para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid
   para_Op%OpGrid(iterm00)%para_FileGrid%Save_FileGrid_done =           &
                     para_Op%OpGrid(iterm00)%para_FileGrid%Save_FileGrid
 END IF

!-----------------------------------------------------------
 IF (debug) THEN
   write(out_unitp,*) 'OpPsiBasisRep'
   write(out_unitp,*) 'OpPsi%RvecB',OpPsi%RvecB
   write(out_unitp,*)
   write(out_unitp,*) 'END ',name_sub
 END IF
!-----------------------------------------------------------

  END SUBROUTINE sub_OpPsi_FOR_SGtype4

  SUBROUTINE sub_OpPsi_OF_ONEDP_FOR_SGtype4(PsiR,iG,tab_l,para_Op)
   USE mod_system
!$ USE omp_lib, only : OMP_GET_THREAD_NUM
   USE mod_nDindex

   USE mod_Coord_KEO,               ONLY : zmatrix, get_Qact, get_d0GG

   USE mod_basis_set_alloc,         ONLY : basis
   USE mod_basis,                   ONLY : Rec_Qact_SG4_with_Tab_iq
   USE mod_basis_RCVec_SGType4,      ONLY : TypeRVec
   USE mod_basis_BtoG_GtoB_SGType4, ONLY : BDP_TO_GDP_OF_SmolyakRep,    &
                                              GDP_TO_BDP_OF_SmolyakRep, &
                                          DerivOp_TO_RDP_OF_SmolaykRep, &
                                          tabR2grid_TO_tabR1_AT_iG,     &
                                          getbis_tab_nq,getbis_tab_nb
  USE mod_Op,                      ONLY : param_Op,write_param_Op
  IMPLICIT NONE

  TYPE (TypeRVec),                    intent(inout)    :: PsiR
  integer,                            intent(in)       :: iG,tab_l(:)

  TYPE (param_Op),                    intent(inout)       :: para_Op

  !local variables
  TYPE (zmatrix), pointer :: mole
  TYPE(basis),    pointer :: BasisnD

  integer :: iq,nb,nq,i,j,D

  integer                               :: derive_termQdyn(2)

  real (kind=Rkind),  allocatable       :: V(:)
  real (kind=Rkind),  allocatable       :: Qact(:)
  real (kind=Rkind),  allocatable       :: GGiq(:,:,:)
  real (kind=Rkind),  allocatable       :: sqRhoOVERJac(:),Jac(:)
  real (kind=Rkind)                     :: Rho

  real (kind=Rkind),  allocatable       :: PsiRj(:,:),PsiRi(:),OpPsiR(:)

  integer,            allocatable       :: tab_nq(:)
  integer,            allocatable       :: tab_nb(:)
  integer,            allocatable       :: tab_iq(:)
  integer :: iterm00
  integer :: err_sub


  !----- for debuging ----------------------------------------------
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub='sub_OpPsi_OF_ONEDP_FOR_SGtype4'
  logical, parameter :: debug = .FALSE.
  !logical, parameter :: debug = .TRUE.
  !-----------------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*) 'nb_bie,nb_baie',para_Op%nb_bie,para_Op%nb_baie
    write(out_unitp,*) 'nb_act1',para_Op%mole%nb_act1
    write(out_unitp,*) 'nb_var',para_Op%mole%nb_var
    CALL flush_perso(out_unitp)
  END IF
  !-----------------------------------------------------------------

   !write(6,*) '================================' ; flush(6)
   !write(6,*) '============ START =============' ; flush(6)

   mole    => para_Op%mole
   BasisnD => para_Op%BasisnD

   D = size(tab_l)

   CALL alloc_NParray(Qact,(/mole%nb_var/),'Qact',name_sub)
   CALL alloc_NParray(tab_nb,(/ D /),'tab_nb',name_sub)
   CALL alloc_NParray(tab_nq,(/ D /),'tab_nq',name_sub)
   CALL alloc_NParray(tab_iq,(/ D /),'tab_iq',name_sub)

   tab_nq(:) = getbis_tab_nq(tab_l,BasisnD%tab_basisPrimSG)
   tab_nb(:) = getbis_tab_nb(tab_l,BasisnD%tab_basisPrimSG)


   nb = size(PsiR%V) ! because PsiR is in basis Rep
   nq = BasisnD%para_SGType2%tab_nq_OF_SRep(iG)

   CALL alloc_NParray(PsiRj,(/ nq,mole%nb_act1 /),'PsiRj',       name_sub)
   CALL alloc_NParray(PsiRi,       (/ nq /),      'PsiRi',       name_sub)


   CALL get_OpGrid_OF_ONEDP_FOR_SG4(iG,tab_l,para_Op,V,GGiq,sqRhoOVERJac,Jac)

   !write(6,*) ' GGiq:',GGiq
   !write(6,*) ' Jac :',Jac
   !write(6,*) ' sqRhoOVERJac :',sqRhoOVERJac
   !write(6,*) 'd0GG ... done' ; flush(6)
   CALL alloc_NParray(OpPsiR,(/nq/),'OpPsiR',name_sub)

   !write(6,*) ' V Basis of PsiR :',PsiR%V

   ! partial B to G
   CALL BDP_TO_GDP_OF_SmolyakRep(PsiR%V,BasisnD%tab_basisPrimSG,        &
                                   tab_l,tab_nq,tab_nb)
   ! now PsiR is on the grid
   !write(6,*) ' R Grid of PsiR :',PsiR%V
   !write(6,*) ' R Grid of PsiR : done',OMP_GET_THREAD_NUM() ; flush(6)

   ! multiplication by sqRhoOVERJac
   PsiR%V(:) = PsiR%V(:) * sqRhoOVERJac(:)

   !write(6,*) ' R*sq Grid of PsiR :',PsiR%V
   !write(6,*) ' R*sq ... Grid of PsiR : done',OMP_GET_THREAD_NUM() ; flush(6)


   ! derivative with respect to Qj
   DO j=1,mole%nb_act1
     derive_termQdyn(:) = (/ mole%liste_QactTOQsym(j),0 /)

     PsiRj(:,j) = PsiR%V(:)
     CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiRj(:,j),BasisnD%tab_basisPrimSG, &
                                         tab_l,tab_nq,derive_termQdyn)

     !write(6,*) 'dQj Grid of PsiR',j,PsiRj(:,j) ; flush(6)
   END DO
   !write(6,*) ' dQj R*sq ... Grid of PsiR : done',OMP_GET_THREAD_NUM() ; flush(6)

   OpPsiR(:) = ZERO
   DO i=1,mole%nb_act1

     PsiRi(:) = ZERO
     DO j=1,mole%nb_act1
       PsiRi(:) = PsiRi(:) + GGiq(:,j,i) * PsiRj(:,j)
     END DO
     PsiRi(:) = PsiRi(:) * Jac(:)
     !write(6,*) ' Jac Gij* ... Grid of SRep : done',iG,i,OMP_GET_THREAD_NUM() ; flush(6)

     derive_termQdyn(:) = (/ mole%liste_QactTOQsym(i),0 /)

     CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiRi(:),BasisnD%tab_basisPrimSG,&
                                         tab_l,tab_nq,derive_termQdyn)

     !write(6,*) 'shape OpPsiR,PsiRi ',shape(OpPsiR),shape(PsiRi)
     !write(6,*) 'DerivOp_TO_RDP_OF_SmolaykRep : done',OMP_GET_THREAD_NUM() ; flush(6)

     OpPsiR(:) = OpPsiR(:) + PsiRi(:)
   END DO

   !write(6,*) ' dQi Jac Gij* ... Grid of PsiR : done',OMP_GET_THREAD_NUM() ; flush(6)

   OpPsiR(:) = (  -HALF*OpPsiR(:)/Jac(:) + PsiR%V(:)*V(:) ) / sqRhoOVERJac(:)
   !write(6,*) ' OpPsiR Grid of PsiR',OpPsiR
   !write(6,*) ' OpPsiR Grid of PsiR : done',OMP_GET_THREAD_NUM() ; flush(6)


   !tranfert OpPsiR (on the grid) to PsiR
   PsiR%V(:) = OpPsiR(:)

   CALL dealloc_NParray(OpPsiR,'OpPsiR', name_sub)


   CALL GDP_TO_BDP_OF_SmolyakRep(PsiR%V,BasisnD%tab_basisPrimSG,&
                                 tab_l,tab_nq,tab_nb)
   ! now PsiR%V is on the Basis
   !write(6,*) ' V Basis of OpPsiR :',PsiR%V

   IF (allocated(PsiRi))        CALL dealloc_NParray(PsiRi,       'PsiRi',       name_sub)
   IF (allocated(PsiRj))        CALL dealloc_NParray(PsiRj,       'PsiRj',       name_sub)
   IF (allocated(GGiq))         CALL dealloc_NParray(GGiq,        'GGiq',        name_sub)
   IF (allocated(tab_nb))       CALL dealloc_NParray(tab_nb,      'tab_nb',      name_sub)
   IF (allocated(tab_nq))       CALL dealloc_NParray(tab_nq,      'tab_nq',      name_sub)
   IF (allocated(tab_iq))       CALL dealloc_NParray(tab_iq,      'tab_iq',      name_sub)
   IF (allocated(Qact))         CALL dealloc_NParray(Qact,        'Qact',        name_sub)
   IF (allocated(sqRhoOVERJac)) CALL dealloc_NParray(sqRhoOVERJac,'sqRhoOVERJac',name_sub)
   IF (allocated(Jac))          CALL dealloc_NParray(Jac,         'Jac',         name_sub)
   IF (allocated(V))            CALL dealloc_NParray(V,           'V',           name_sub)

   !write(6,*) '============ END ===============' ; flush(6)
   !write(6,*) '================================' ; flush(6)

  !-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'END ',name_sub
    CALL flush_perso(out_unitp)
  END IF
  !-----------------------------------------------------------
  !CALL UnCheck_mem() ; stop
  END SUBROUTINE sub_OpPsi_OF_ONEDP_FOR_SGtype4

 SUBROUTINE sub_TabOpPsi_FOR_SGtype4(Psi,OpPsi,para_Op)
 USE mod_system
!$ USE omp_lib, only : OMP_GET_THREAD_NUM
 USE mod_nDindex

 USE mod_Coord_KEO,                     ONLY : zmatrix

 USE mod_basis_set_alloc,          ONLY : basis
 USE mod_basis_BtoG_GtoB_SGType4,  ONLY : tabPackedBasis_TO_tabR_AT_iG, &
                  tabR_AT_iG_TO_tabPackedBasis,TypeRVec,dealloc_TypeRVec

 USE mod_psi_set_alloc,            ONLY : param_psi,ecri_psi

 USE mod_Op,                       ONLY : param_Op,write_param_Op
 IMPLICIT NONE

 TYPE (param_psi), intent(in)      :: Psi(:)
 TYPE (param_psi), intent(inout)   :: OpPsi(:)

 TYPE (param_Op),  intent(inout)   :: para_Op

 ! local variables
 TYPE (zmatrix), pointer :: mole
 TYPE (basis),   pointer :: BasisnD

 TYPE (TypeRVec),   allocatable    :: PsiR(:)

 integer       :: ib,i,iG,iiG,nb_thread,itab,ith,iterm00
 integer       :: tab_l(para_Op%BasisnD%nb_basis)
 logical       :: not_init

 !----- for debuging ----------------------------------------------
 character (len=*), parameter :: name_sub='sub_TabOpPsi_FOR_SGtype4'
 logical, parameter :: debug = .FALSE.
 !logical, parameter :: debug = .TRUE.
 !-----------------------------------------------------------------
 IF (debug) THEN
   write(out_unitp,*) 'BEGINNING ',name_sub
   write(out_unitp,*) 'nb_bie,nb_baie',para_Op%nb_bie,para_Op%nb_baie
   write(out_unitp,*) 'nb_act1',para_Op%mole%nb_act1
   write(out_unitp,*) 'nb_var',para_Op%mole%nb_var
   !CALL flush_perso(out_unitp)
   !CALL write_param_Op(para_Op)
   write(out_unitp,*)
   write(out_unitp,*) 'PsiBasisRep'
   !DO itab=1,size(Psi)
     !CALL ecri_psi(Psi=Psi(itab))
   !END DO
   CALL flush_perso(out_unitp)
 END IF
 !-----------------------------------------------------------------
  mole    => para_Op%mole
  BasisnD => para_Op%BasisnD

  IF (Psi(1)%cplx) STOP 'cplx'
  IF (para_Op%nb_bie /= 1) STOP 'nb_bie /= 1'

  IF (OpPsi_omp == 0) THEN
    nb_thread = 1
  ELSE
    nb_thread = OpPsi_maxth
  END IF

 nb_mult_OpPsi = 0
 DO itab=1,size(Psi)
   OpPsi(itab)           = Psi(itab) ! for the allocation. It has to be changed!
   OpPsi(itab)%RvecB(:)  = ZERO
 END DO

 IF (print_level > 0 .AND. BasisnD%para_SGType2%nb_SG > 10**5 ) THEN
   write(out_unitp,'(a)')              'OpPsi SG4 (%): [--0-10-20-30-40-50-60-70-80-90-100]'
   write(out_unitp,'(a)',ADVANCE='no') 'OpPsi SG4 (%): ['
   CALL flush_perso(out_unitp)
 END IF

 IF (OpPsi_omp == 0) THEN
   DO iG=1,BasisnD%para_SGType2%nb_SG

     !write(6,*) 'iG',iG ;    CALL flush_perso(out_unitp)
     !transfert part of the psi(:)%RvecB(:) to PsiR(itab)%V
     allocate(PsiR(size(Psi)))
     DO itab=1,size(Psi)
       CALL tabPackedBasis_TO_tabR_AT_iG(PsiR(itab)%V,psi(itab)%RvecB,    &
                                         iG,BasisnD%para_SGType2)
       !write(6,*) 'iG,itab,PsiR',iG,itab,PsiR(itab)%V  ;    CALL flush_perso(out_unitp)
     END DO
     !write(6,*) 'iG',iG, ' unpack done' ;    CALL flush_perso(out_unitp)

     CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4_ompGrid(PsiR,iG,para_Op)
     !write(6,*) 'iG',iG, ' TabOpPsi done' ;    CALL flush_perso(out_unitp)

     !transfert back, PsiR(itab)%V (HPsi) to the psi(:)%RvecB(:)
     DO itab=1,size(Psi)
       CALL tabR_AT_iG_TO_tabPackedBasis(OpPsi(itab)%RvecB,PsiR(itab)%V,  &
                     iG,BasisnD%para_SGType2,BasisnD%WeightSG(iG))
     END DO
     !write(6,*) 'iG',iG, ' pack done' ;    CALL flush_perso(out_unitp)

     !deallocate PsiR(:)
     DO itab=1,size(Psi)
       CALL dealloc_TypeRVec(PsiR(itab))
     END DO
     deallocate(PsiR)

     IF (print_level > 0  .AND. BasisnD%para_SGType2%nb_SG > 10**5 .AND. &
         mod(iG,max(1,BasisnD%para_SGType2%nb_SG/10)) == 0) THEN
       write(out_unitp,'(a)',ADVANCE='no') '---'
       CALL flush_perso(out_unitp)
     END IF

     !write(6,*) 'iG done:',iG ; flush(6)
   END DO
 ELSE

   !to be sure to have the correct number of threads
   nb_thread = BasisnD%para_SGType2%nb_threads

   !$OMP parallel                                                &
   !$OMP default(none)                                           &
   !$OMP shared(Psi,OpPsi)                                       &
   !$OMP shared(para_Op,BasisnD,print_level,out_unitp)           &
   !$OMP private(itab,iG,iiG,tab_l,PsiR,ith)                     &
   !$OMP num_threads(BasisnD%para_SGType2%nb_threads)

   !--------------------------------------------------------------
   !-- For the initialization of tab_l(:) and the use of ADD_ONE_TO_nDindex in the parallel loop
   ith = 0
   !$ ith = OMP_GET_THREAD_NUM()
   tab_l(:) = BasisnD%para_SGType2%nDval_init(:,ith+1)
   !--------------------------------------------------------------

   ! we are not using the parallel do, to be able to use the correct initialized tab_l with nDval_init
   DO iG=BasisnD%para_SGType2%iG_th(ith+1),BasisnD%para_SGType2%fG_th(ith+1)

     CALL ADD_ONE_TO_nDindex(BasisnD%para_SGType2%nDind_SmolyakRep,tab_l,iG=iG)

     !write(6,*) 'iG',iG
     !transfert part of the psi(:)%RvecB(:) to PsiR(itab)%V
     allocate(PsiR(size(Psi)))
     DO itab=1,size(Psi)
       CALL tabPackedBasis_TO_tabR_AT_iG(PsiR(itab)%V,psi(itab)%RvecB,    &
                                         iG,BasisnD%para_SGType2)
       !write(6,*) 'iG,itab,PsiR',iG,itab,PsiR(itab)%V
     END DO

     CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4(PsiR,iG,tab_l,para_Op)

     !transfert back, PsiR(itab)%V (HPsi) to the psi(:)%RvecB(:)
     DO itab=1,size(Psi)
       CALL tabR_AT_iG_TO_tabPackedBasis(OpPsi(itab)%RvecB,PsiR(itab)%V,  &
                     iG,BasisnD%para_SGType2,BasisnD%WeightSG(iG))
     END DO

     !deallocate PsiR(:)
     DO itab=1,size(Psi)
       CALL dealloc_TypeRVec(PsiR(itab))
     END DO
     deallocate(PsiR)

     IF (print_level > 0  .AND. BasisnD%para_SGType2%nb_SG > 10**5 .AND. &
         mod(iG,max(1,BasisnD%para_SGType2%nb_SG/10)) == 0) THEN
       write(out_unitp,'(a)',ADVANCE='no') '---'
       CALL flush_perso(out_unitp)
     END IF

     !write(6,*) 'iG done:',iG ; flush(6)
   END DO
  !$OMP   END PARALLEL
END IF

 IF (print_level > 0 .AND. BasisnD%para_SGType2%nb_SG > 10**5) THEN
   write(out_unitp,'(a)',ADVANCE='yes') '----]'
 END IF
 CALL flush_perso(out_unitp)

 iterm00 = para_Op%derive_term_TO_iterm(0,0)
 IF (associated(para_Op%OpGrid)) THEN
   para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_done =            &
                     para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid
   para_Op%OpGrid(iterm00)%para_FileGrid%Save_FileGrid_done =           &
                     para_Op%OpGrid(iterm00)%para_FileGrid%Save_FileGrid
 END IF

!-----------------------------------------------------------
 IF (debug) THEN
   write(out_unitp,*) 'OpPsiGridRep'
   !DO itab=1,size(Psi)
     !CALL ecri_psi(Psi=OpPsi(itab))
   !END DO
   write(out_unitp,*)
   write(out_unitp,*) 'END ',name_sub
 END IF
!-----------------------------------------------------------

  END SUBROUTINE sub_TabOpPsi_FOR_SGtype4

  SUBROUTINE sub_TabOpPsi_OF_ONEDP_FOR_SGtype4(PsiR,iG,tab_l,para_Op)
  USE mod_system
  USE mod_nDindex

  USE mod_Coord_KEO,               ONLY : zmatrix, get_Qact, get_d0GG

  USE mod_basis_set_alloc,         ONLY : basis
  USE mod_basis,                   ONLY : Rec_Qact_SG4_with_Tab_iq
  USE mod_basis_BtoG_GtoB_SGType4, ONLY : TypeRVec,                     &
                     BDP_TO_GDP_OF_SmolyakRep,GDP_TO_BDP_OF_SmolyakRep, &
                  DerivOp_TO_RDP_OF_SmolaykRep,tabR2grid_TO_tabR1_AT_iG,&
                     getbis_tab_nq,getbis_tab_nb
  USE mod_Op,                      ONLY : param_Op,write_param_Op
  IMPLICIT NONE

  TYPE (TypeRVec),                    intent(inout)    :: PsiR(:)
  integer,                            intent(in)       :: iG,tab_l(:)

  TYPE (param_Op),                    intent(inout)    :: para_Op

  !local variables
  TYPE (zmatrix), pointer :: mole
  TYPE(basis),    pointer :: BasisnD

  integer :: iq,nb,nq,i,j,itab,D

  integer                               :: derive_termQdyn(2)

  real (kind=Rkind),  allocatable       :: V(:)
  real (kind=Rkind),  allocatable       :: Qact(:)
  real (kind=Rkind),  allocatable       :: GGiq(:,:,:)
  real (kind=Rkind),  allocatable       :: sqRhoOVERJac(:),Jac(:)
  real (kind=Rkind)                     :: Rho

  real (kind=Rkind),  allocatable       :: PsiRj(:,:),PsiRi(:),OpPsiR(:)

  integer,            allocatable       :: tab_nq(:)
  integer,            allocatable       :: tab_nb(:)
  integer,            allocatable       :: tab_iq(:)
  integer :: iterm00
  integer :: err_sub


  !----- for debuging ----------------------------------------------
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub='sub_TabOpPsi_OF_ONEDP_FOR_SGtype4'
  logical, parameter :: debug = .FALSE.
  !logical, parameter :: debug = .TRUE.
  !-----------------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*) 'nb_bie,nb_baie',para_Op%nb_bie,para_Op%nb_baie
    write(out_unitp,*) 'nb_act1',para_Op%mole%nb_act1
    write(out_unitp,*) 'nb_var',para_Op%mole%nb_var
    CALL flush_perso(out_unitp)
  END IF
  !-----------------------------------------------------------------
  !CALL Check_mem()

   !write(6,*) '================================' ; flush(6)
   !write(6,*) '============ START =============' ; flush(6)

   mole    => para_Op%mole
   BasisnD => para_Op%BasisnD

   D = size(tab_l)

   CALL alloc_NParray(Qact,(/mole%nb_var/),'Qact',name_sub)
   CALL alloc_NParray(tab_nb,(/ D /),'tab_nb',name_sub)
   CALL alloc_NParray(tab_nq,(/ D /),'tab_nq',name_sub)
   CALL alloc_NParray(tab_iq,(/ D /),'tab_iq',name_sub)

   tab_nq(:) = getbis_tab_nq(tab_l,BasisnD%tab_basisPrimSG)
   tab_nb(:) = getbis_tab_nb(tab_l,BasisnD%tab_basisPrimSG)


   nb = size(PsiR(1)%V) ! because PsiR is in basis Rep
   nq = BasisnD%para_SGType2%tab_nq_OF_SRep(iG)

   CALL alloc_NParray(PsiRj,(/ nq,mole%nb_act1 /),'PsiRj',       name_sub)
   CALL alloc_NParray(PsiRi,       (/ nq /),      'PsiRi',       name_sub)


   CALL get_OpGrid_OF_ONEDP_FOR_SG4(iG,tab_l,para_Op,V,GGiq,sqRhoOVERJac,Jac)


   DO itab=1,size(PsiR)

     CALL alloc_NParray(OpPsiR,(/nq/),'OpPsiR',name_sub)

     ! partial B to G
     CALL BDP_TO_GDP_OF_SmolyakRep(PsiR(itab)%V,BasisnD%tab_basisPrimSG,&
                                   tab_l,tab_nq,tab_nb)
     ! now PsiR is on the grid
     !write(6,*) ' R Grid of PsiR(itab) :',PsiR(itab)%V
     !write(6,*) ' R Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(6)

     ! multiplication by sqRhoOVERJac
     PsiR(itab)%V(:) = PsiR(itab)%V(:) * sqRhoOVERJac(:)

     !write(6,*) ' R*sq Grid of PsiR(itab) :',PsiR
     !write(6,*) ' R*sq ... Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(6)


     ! derivative with respect to Qj
     DO j=1,mole%nb_act1
       derive_termQdyn(:) = (/ mole%liste_QactTOQsym(j),0 /)

       PsiRj(:,j) = PsiR(itab)%V(:)
       CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiRj(:,j),BasisnD%tab_basisPrimSG, &
                                         tab_l,tab_nq,derive_termQdyn)

       !write(6,*) 'dQj Grid of PsiR(itab)',j,PsiRj(:,j) ; flush(6)
     END DO
     !write(6,*) ' dQj R*sq ... Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(6)

     OpPsiR(:) = ZERO
     DO i=1,mole%nb_act1

       PsiRi(:) = ZERO
       DO j=1,mole%nb_act1
         PsiRi(:) = PsiRi(:) + GGiq(:,j,i) * PsiRj(:,j)
       END DO
       PsiRi(:) = PsiRi(:) * Jac(:)
       !write(6,*) ' Jac Gij* ... Grid of SRep : done',iG,itab,i,OMP_GET_THREAD_NUM() ; flush(6)

       derive_termQdyn(:) = (/ mole%liste_QactTOQsym(i),0 /)

       CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiRi(:),BasisnD%tab_basisPrimSG,&
                                         tab_l,tab_nq,derive_termQdyn)

       OpPsiR(:) = OpPsiR(:) + PsiRi(:)

     END DO
     !write(6,*) ' dQi Jac Gij* ... Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(6)

     OpPsiR(:) = (  -HALF*OpPsiR(:)/Jac(:) + PsiR(itab)%V(:)*V(:) ) / sqRhoOVERJac(:)
     !write(6,*) ' OpPsiR Grid of PsiR(itab)',OpPsiR
     !write(6,*) ' OpPsiR Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(6)

     !tranfert OpPsiR (on the grid) to PsiR(itab)
     PsiR(itab)%V(:) = OpPsiR(:)



     CALL dealloc_NParray(OpPsiR,'OpPsiR', name_sub)


     CALL GDP_TO_BDP_OF_SmolyakRep(PsiR(itab)%V,BasisnD%tab_basisPrimSG,&
                                   tab_l,tab_nq,tab_nb)
     ! now PsiR(itab)%V is on the Basis


   END DO

   IF (allocated(PsiRi))        CALL dealloc_NParray(PsiRi,       'PsiRi',       name_sub)
   IF (allocated(PsiRj))        CALL dealloc_NParray(PsiRj,       'PsiRj',       name_sub)
   IF (allocated(GGiq))         CALL dealloc_NParray(GGiq,        'GGiq',        name_sub)
   IF (allocated(tab_nb))       CALL dealloc_NParray(tab_nb,      'tab_nb',      name_sub)
   IF (allocated(tab_nq))       CALL dealloc_NParray(tab_nq,      'tab_nq',      name_sub)
   IF (allocated(tab_iq))       CALL dealloc_NParray(tab_iq,      'tab_iq',      name_sub)
   IF (allocated(Qact))         CALL dealloc_NParray(Qact,        'Qact',        name_sub)
   IF (allocated(sqRhoOVERJac)) CALL dealloc_NParray(sqRhoOVERJac,'sqRhoOVERJac',name_sub)
   IF (allocated(Jac))          CALL dealloc_NParray(Jac,         'Jac',         name_sub)
   IF (allocated(V))            CALL dealloc_NParray(V,           'V',           name_sub)

   !write(6,*) '============ END ===============' ; flush(6)
   !write(6,*) '================================' ; flush(6)

  !-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'END ',name_sub
    CALL flush_perso(out_unitp)
  END IF
  !-----------------------------------------------------------
  !CALL UnCheck_mem() ; stop
  END SUBROUTINE sub_TabOpPsi_OF_ONEDP_FOR_SGtype4

  ! with this subroutine, PsiR(:)%V(:) are on the grid
  SUBROUTINE sub_TabOpPsi_OF_ONEGDP_FOR_SGtype4(PsiR,iG,para_Op)
  USE mod_system
  USE mod_nDindex

  USE mod_Coord_KEO,               ONLY : zmatrix, get_Qact, get_d0GG

  USE mod_basis_set_alloc,         ONLY : basis
  USE mod_basis,                   ONLY : Rec_Qact_SG4_with_Tab_iq
  USE mod_basis_BtoG_GtoB_SGType4, ONLY : TypeRVec,                     &
                     BDP_TO_GDP_OF_SmolyakRep,GDP_TO_BDP_OF_SmolyakRep, &
                   DerivOp_TO_RDP_OF_SmolaykRep,tabR2grid_TO_tabR1_AT_iG,&
                     getbis_tab_nq,getbis_tab_nb
  USE mod_Op,                      ONLY : param_Op,write_param_Op
  IMPLICIT NONE

  TYPE (TypeRVec),                    intent(inout)    :: PsiR(:)
  integer,                            intent(in)       :: iG

  TYPE (param_Op),                    intent(inout)    :: para_Op


  !local variables
  TYPE (zmatrix), pointer :: mole
  TYPE(basis),    pointer :: BasisnD

  integer :: iq,nb,nq,i,j,itab,D

  integer                               :: derive_termQdyn(2)

  real (kind=Rkind),  allocatable       :: V(:)
  real (kind=Rkind),  allocatable       :: Qact(:)
  real (kind=Rkind),  allocatable       :: GGiq(:,:,:)
  real (kind=Rkind),  allocatable       :: sqRhoOVERJac(:),Jac(:)
  real (kind=Rkind)                     :: Rho

  real (kind=Rkind),  allocatable       :: PsiRj(:,:),PsiRi(:),OpPsiR(:)

  integer,            allocatable       :: tab_nq(:)
  integer,            allocatable       :: tab_nb(:)
  integer,            allocatable       :: tab_l(:)
  integer,            allocatable       :: tab_iq(:)

  integer :: iterm00
  integer :: err_sub


  !----- for debuging ----------------------------------------------
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub='sub_TabOpPsi_OF_ONEGDP_FOR_SGtype4'
  !logical, parameter :: debug = .FALSE.
  logical, parameter :: debug = .TRUE.
  !-----------------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*) 'nb_bie,nb_baie',para_Op%nb_bie,para_Op%nb_baie
    write(out_unitp,*) 'nb_act1',para_Op%mole%nb_act1
    write(out_unitp,*) 'nb_var',para_Op%mole%nb_var
    CALL flush_perso(out_unitp)
  END IF
  !-----------------------------------------------------------------
  !CALL Check_mem()

   !write(6,*) '================================' ; flush(6)
   !write(6,*) '============ START =============' ; flush(6)

   mole    => para_Op%mole
   BasisnD => para_Op%BasisnD

   D = size(BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval,dim=1)

   CALL alloc_NParray(Qact,(/mole%nb_var/),'Qact',name_sub)
   CALL alloc_NParray(tab_l ,(/ D /),'tab_l', name_sub)
   CALL alloc_NParray(tab_nb,(/ D /),'tab_nb',name_sub)
   CALL alloc_NParray(tab_nq,(/ D /),'tab_nq',name_sub)
   CALL alloc_NParray(tab_iq,(/ D /),'tab_iq',name_sub)

   tab_l(:)  = BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval(:,iG)
   tab_nq(:) = getbis_tab_nq(tab_l,BasisnD%tab_basisPrimSG)
   tab_nb(:) = getbis_tab_nb(tab_l,BasisnD%tab_basisPrimSG)


   nb = size(PsiR(1)%V) ! because PsiR is in basis Rep
   nq = BasisnD%para_SGType2%tab_nq_OF_SRep(iG)

   CALL alloc_NParray(PsiRj,(/ nq,mole%nb_act1 /),'PsiRj',       name_sub)
   CALL alloc_NParray(PsiRi,       (/ nq /),      'PsiRi',       name_sub)

   CALL get_OpGrid_OF_ONEDP_FOR_SG4(iG,tab_l,para_Op,V,GGiq,sqRhoOVERJac,Jac)

   DO itab=1,size(PsiR)

     CALL alloc_NParray(OpPsiR,(/nq/),'OpPsiR',name_sub)

     ! partial B to G
     !!!!! no because PsiR(itab)%V(:) is already on the grid


     ! multiplication by sqRhoOVERJac
     PsiR(itab)%V(:) = PsiR(itab)%V(:) * sqRhoOVERJac(:)

     !write(6,*) ' R*sq Grid of PsiR(itab) :',PsiR
     !write(6,*) ' R*sq ... Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(6)


     ! derivative with respect to Qj
     DO j=1,mole%nb_act1
       derive_termQdyn(:) = (/ mole%liste_QactTOQsym(j),0 /)

       PsiRj(:,j) = PsiR(itab)%V(:)
       CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiRj(:,j),BasisnD%tab_basisPrimSG, &
                                         tab_l,tab_nq,derive_termQdyn)

       !write(6,*) 'dQj Grid of PsiR(itab)',j,PsiRj(:,j) ; flush(6)
     END DO
     !write(6,*) ' dQj R*sq ... Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(6)

     OpPsiR(:) = ZERO
     DO i=1,mole%nb_act1

       PsiRi(:) = ZERO
       DO j=1,mole%nb_act1
         PsiRi(:) = PsiRi(:) + GGiq(:,j,i) * PsiRj(:,j)
       END DO
       PsiRi(:) = PsiRi(:) * Jac(:)
       !write(6,*) ' Jac Gij* ... Grid of SRep : done',iG,itab,i,OMP_GET_THREAD_NUM() ; flush(6)

       derive_termQdyn(:) = (/ mole%liste_QactTOQsym(i),0 /)

       CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiRi(:),BasisnD%tab_basisPrimSG,&
                                         tab_l,tab_nq,derive_termQdyn)

       OpPsiR(:) = OpPsiR(:) + PsiRi(:)

     END DO
     !write(6,*) ' dQi Jac Gij* ... Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(6)

     OpPsiR(:) = (  -HALF*OpPsiR(:)/Jac(:) + PsiR(itab)%V(:)*V(:) ) / sqRhoOVERJac(:)
     !write(6,*) ' OpPsiR Grid of PsiR(itab)',OpPsiR
     !write(6,*) ' OpPsiR Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(6)

     !tranfert OpPsiR (on the grid) to PsiR(itab)
     PsiR(itab)%V(:) = OpPsiR(:)

     CALL dealloc_NParray(OpPsiR,'OpPsiR', name_sub)

     ! we kept PsiR(itab)%V on the grid


   END DO

   IF (allocated(PsiRi))        CALL dealloc_NParray(PsiRi,       'PsiRi',       name_sub)
   IF (allocated(PsiRj))        CALL dealloc_NParray(PsiRj,       'PsiRj',       name_sub)
   IF (allocated(GGiq))         CALL dealloc_NParray(GGiq,        'GGiq',        name_sub)
   IF (allocated(tab_nb))       CALL dealloc_NParray(tab_nb,      'tab_nb',      name_sub)
   IF (allocated(tab_nq))       CALL dealloc_NParray(tab_nq,      'tab_nq',      name_sub)
   IF (allocated(tab_iq))       CALL dealloc_NParray(tab_iq,      'tab_iq',      name_sub)
   IF (allocated(tab_l))        CALL dealloc_NParray(tab_l,       'tab_l',       name_sub)
   IF (allocated(Qact))         CALL dealloc_NParray(Qact,        'Qact',        name_sub)
   IF (allocated(sqRhoOVERJac)) CALL dealloc_NParray(sqRhoOVERJac,'sqRhoOVERJac',name_sub)
   IF (allocated(Jac))          CALL dealloc_NParray(Jac,         'Jac',         name_sub)
   IF (allocated(V))            CALL dealloc_NParray(V,           'V',           name_sub)

   !write(6,*) '============ END ===============' ; flush(6)
   !write(6,*) '================================' ; flush(6)

  !-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'END ',name_sub
    CALL flush_perso(out_unitp)
  END IF
  !-----------------------------------------------------------
  !CALL UnCheck_mem() ; stop
  END SUBROUTINE sub_TabOpPsi_OF_ONEGDP_FOR_SGtype4

  SUBROUTINE sub_TabOpPsi_OF_ONEGDP_WithOp_FOR_SGtype4(PsiR,iG,para_Op, &
                                                   V,GG,sqRhoOVERJac,Jac)
  USE mod_system

  USE mod_Coord_KEO,               ONLY : zmatrix, get_Qact

  USE mod_basis_set_alloc,         ONLY : basis
  USE mod_basis_BtoG_GtoB_SGType4, ONLY : TypeRVec,                     &
                     BDP_TO_GDP_OF_SmolyakRep,GDP_TO_BDP_OF_SmolyakRep, &
                   DerivOp_TO_RDP_OF_SmolaykRep,tabR2grid_TO_tabR1_AT_iG,&
                     getbis_tab_nq,getbis_tab_nb
  USE mod_Op,                      ONLY : param_Op,write_param_Op
  IMPLICIT NONE

  TYPE (TypeRVec),                    intent(inout)    :: PsiR(:)
  integer,                            intent(in)       :: iG

  TYPE (param_Op),                    intent(inout)    :: para_Op
  real (kind=Rkind),                  intent(in)       :: V(:)
  real (kind=Rkind),                  intent(in)       :: GG(:,:,:)
  real (kind=Rkind),                  intent(in)       :: sqRhoOVERJac(:),Jac(:)

  !local variables
  TYPE (zmatrix), pointer :: mole
  TYPE(basis),    pointer :: BasisnD

  integer :: iq,nb,nq,i,j,itab,D

  integer                               :: derive_termQdyn(2)


  real (kind=Rkind)                     :: Rho

  real (kind=Rkind),  allocatable       :: PsiRj(:,:),PsiRi(:),OpPsiR(:)

  integer,            allocatable       :: tab_nq(:)
  integer,            allocatable       :: tab_nb(:)
  integer,            allocatable       :: tab_l(:)
  integer :: err_sub


  !----- for debuging ----------------------------------------------
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub='sub_TabOpPsi_OF_ONEGDP_WithOp_FOR_SGtype4'
  logical, parameter :: debug = .FALSE.
  !logical, parameter :: debug = .TRUE.
  !-----------------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*) 'nb_bie,nb_baie',para_Op%nb_bie,para_Op%nb_baie
    write(out_unitp,*) 'nb_act1',para_Op%mole%nb_act1
    write(out_unitp,*) 'nb_var',para_Op%mole%nb_var
    CALL flush_perso(out_unitp)
  END IF
  !-----------------------------------------------------------------
  !CALL Check_mem()

   !write(6,*) '================================' ; flush(6)
   !write(6,*) '============ START =============' ; flush(6)

   mole    => para_Op%mole
   BasisnD => para_Op%BasisnD

   D = size(BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval,dim=1)

   CALL alloc_NParray(tab_l ,(/ D /),'tab_l', name_sub)
   CALL alloc_NParray(tab_nb,(/ D /),'tab_nb',name_sub)
   CALL alloc_NParray(tab_nq,(/ D /),'tab_nq',name_sub)

   tab_l(:)  = BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval(:,iG)
   tab_nq(:) = getbis_tab_nq(tab_l,BasisnD%tab_basisPrimSG)
   tab_nb(:) = getbis_tab_nb(tab_l,BasisnD%tab_basisPrimSG)


   nb = size(PsiR(1)%V) ! because PsiR is in basis Rep
   nq = BasisnD%para_SGType2%tab_nq_OF_SRep(iG)

   CALL alloc_NParray(PsiRj,(/ nq,mole%nb_act1 /),'PsiRj',       name_sub)
   CALL alloc_NParray(PsiRi,       (/ nq /),      'PsiRi',       name_sub)

   DO itab=1,size(PsiR)

     CALL alloc_NParray(OpPsiR,(/nq/),'OpPsiR',name_sub)

     ! partial B to G
     !!!!! no because PsiR(itab)%V(:) is already on the grid


     ! multiplication by sqRhoOVERJac
     PsiR(itab)%V(:) = PsiR(itab)%V(:) * sqRhoOVERJac(:)

     !write(6,*) ' R*sq Grid of PsiR(itab) :',PsiR
     !write(6,*) ' R*sq ... Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(6)


     ! derivative with respect to Qj
     DO j=1,mole%nb_act1
       derive_termQdyn(:) = (/ mole%liste_QactTOQsym(j),0 /)

       PsiRj(:,j) = PsiR(itab)%V(:)
       CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiRj(:,j),BasisnD%tab_basisPrimSG, &
                                         tab_l,tab_nq,derive_termQdyn)

       !write(6,*) 'dQj Grid of PsiR(itab)',j,PsiRj(:,j) ; flush(6)
     END DO
     !write(6,*) ' dQj R*sq ... Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(6)

     OpPsiR(:) = ZERO
     DO i=1,mole%nb_act1

       PsiRi(:) = ZERO
       DO j=1,mole%nb_act1
         PsiRi(:) = PsiRi(:) + GG(:,j,i) * PsiRj(:,j)
       END DO
       PsiRi(:) = PsiRi(:) * Jac(:)
       !write(6,*) ' Jac Gij* ... Grid of SRep : done',iG,itab,i,OMP_GET_THREAD_NUM() ; flush(6)

       derive_termQdyn(:) = (/ mole%liste_QactTOQsym(i),0 /)

       CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiRi(:),BasisnD%tab_basisPrimSG,&
                                         tab_l,tab_nq,derive_termQdyn)

       OpPsiR(:) = OpPsiR(:) + PsiRi(:)

     END DO
     !write(6,*) ' dQi Jac Gij* ... Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(6)

     OpPsiR(:) = (  -HALF*OpPsiR(:)/Jac(:) + PsiR(itab)%V(:)*V(:) ) / sqRhoOVERJac(:)
     !write(6,*) ' OpPsiR Grid of PsiR(itab)',OpPsiR
     !write(6,*) ' OpPsiR Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(6)

     !tranfert OpPsiR (on the grid) to PsiR(itab)
     PsiR(itab)%V(:) = OpPsiR(:)



     CALL dealloc_NParray(OpPsiR,'OpPsiR', name_sub)

     ! we kept PsiR(itab)%V on the grid


   END DO

   IF (allocated(PsiRi))        CALL dealloc_NParray(PsiRi,       'PsiRi',       name_sub)
   IF (allocated(PsiRj))        CALL dealloc_NParray(PsiRj,       'PsiRj',       name_sub)
   IF (allocated(tab_nb))       CALL dealloc_NParray(tab_nb,      'tab_nb',      name_sub)
   IF (allocated(tab_nq))       CALL dealloc_NParray(tab_nq,      'tab_nq',      name_sub)
   IF (allocated(tab_l))        CALL dealloc_NParray(tab_l,       'tab_l',       name_sub)

   !write(6,*) '============ END ===============' ; flush(6)
   !write(6,*) '================================' ; flush(6)

  !-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'END ',name_sub
    CALL flush_perso(out_unitp)
  END IF
  !-----------------------------------------------------------
  !CALL UnCheck_mem() ; stop
  END SUBROUTINE sub_TabOpPsi_OF_ONEGDP_WithOp_FOR_SGtype4

  SUBROUTINE sub_TabOpPsi_OF_ONEDP_FOR_SGtype4_ompGrid(PsiR,iG,para_Op)
  USE mod_system
  USE mod_nDindex

  USE mod_Coord_KEO,               ONLY : zmatrix, get_Qact, get_d0GG

  USE mod_basis_set_alloc,         ONLY : basis
  USE mod_basis,                   ONLY : Rec_Qact_SG4
  USE mod_basis_BtoG_GtoB_SGType4, ONLY : TypeRVec,                     &
                     BDP_TO_GDP_OF_SmolyakRep,GDP_TO_BDP_OF_SmolyakRep, &
                   DerivOp_TO_RDP_OF_SmolaykRep,tabR2grid_TO_tabR1_AT_iG,&
                     getbis_tab_nq,getbis_tab_nb
  USE mod_Op,                      ONLY : param_Op,write_param_Op
  IMPLICIT NONE

  TYPE (TypeRVec),                    intent(inout)    :: PsiR(:)
  integer,                            intent(in)       :: iG

  TYPE (param_Op),                    intent(inout)    :: para_Op


  !local variables
  TYPE (zmatrix), pointer :: mole
  TYPE(basis),    pointer :: BasisnD

  integer :: iq,nb,nq,i,j,itab,D

  integer                               :: derive_termQdyn(2)

  real (kind=Rkind),  allocatable       :: V(:)
  real (kind=Rkind),  allocatable       :: Qact(:)
  real (kind=Rkind),  allocatable       :: GGiq(:,:,:)
  TYPE (Type_nDindex)                   :: nDind_DPG    ! multidimensional DP index
  real (kind=Rkind),  allocatable       :: sqRhoOVERJac(:),Jac(:)
  real (kind=Rkind)                     :: Rho

  real (kind=Rkind),  allocatable       :: PsiRj(:,:),PsiRi(:),OpPsiR(:)

  integer,            allocatable       :: tab_nq(:)
  integer,            allocatable       :: tab_nb(:)
  integer,            allocatable       :: tab_l(:)
  integer :: iterm00
  integer :: err_sub,nb_thread


  !----- for debuging ----------------------------------------------
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub='sub_TabOpPsi_OF_ONEDP_FOR_SGtype4_ompGrid'
  logical, parameter :: debug = .FALSE.
  !logical, parameter :: debug = .TRUE.
  !-----------------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*) 'nb_bie,nb_baie',para_Op%nb_bie,para_Op%nb_baie
    write(out_unitp,*) 'nb_act1',para_Op%mole%nb_act1
    write(out_unitp,*) 'nb_var',para_Op%mole%nb_var
    CALL flush_perso(out_unitp)
  END IF
  !-----------------------------------------------------------------
  !CALL Check_mem()

   !write(6,*) '================================' ; flush(6)
   !write(6,*) '============ START =============' ; flush(6)


  IF (Grid_omp == 0 .AND. OpPsi_omp > 0) THEN
    nb_thread = 1
  ELSE
    nb_thread = Grid_maxth
  END IF


   mole    => para_Op%mole
   BasisnD => para_Op%BasisnD

   D = size(BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval,dim=1)

   CALL alloc_NParray(tab_l ,(/ D /),'tab_l', name_sub)
   CALL alloc_NParray(tab_nb,(/ D /),'tab_nb',name_sub)
   CALL alloc_NParray(tab_nq,(/ D /),'tab_nq',name_sub)

   tab_l(:)  = BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval(:,iG)
   tab_nq(:) = getbis_tab_nq(tab_l,BasisnD%tab_basisPrimSG)
   tab_nb(:) = getbis_tab_nb(tab_l,BasisnD%tab_basisPrimSG)


   nb = size(PsiR(1)%V) ! because PsiR is in basis Rep
   nq = BasisnD%para_SGType2%tab_nq_OF_SRep(iG)

   CALL alloc_NParray(PsiRj,(/ nq,mole%nb_act1 /),'PsiRj',       name_sub)
   CALL alloc_NParray(PsiRi,       (/ nq /),      'PsiRi',       name_sub)
   CALL alloc_NParray(sqRhoOVERJac,(/ nq /),      'sqRhoOVERJac',name_sub)
   CALL alloc_NParray(Jac,         (/ nq /),      'Jac',         name_sub)


   CALL get_OpGrid_OF_ONEDP_FOR_SG4(iG,tab_l,para_Op,V)

   ! G calculation
   CALL alloc_NParray(GGiq,(/nq,mole%nb_act1,mole%nb_act1/),'GGiq',name_sub)
   CALL init_nDindexPrim(nDind_DPG,BasisnD%nb_basis,tab_nq,type_OF_nDindex=-1)

 !$OMP parallel                                               &
 !$OMP default(none)                                          &
 !$OMP shared(mole,BasisnD,para_Op,nq,tab_l,nDind_DPG,iG)     &
 !$OMP shared(sqRhoOVERJac,Jac,GGiq)                          &
 !$OMP private(iq,Qact,err_sub,Rho)                           &
 !$OMP num_threads(nb_thread)
   CALL alloc_NParray(Qact,(/mole%nb_var/),'Qact',name_sub)
 !$OMP   DO SCHEDULE(STATIC)
   DO iq=1,nq
     CALL get_Qact(Qact,mole%ActiveTransfo) ! rigid, flexible coordinates
     CALL Rec_Qact_SG4(Qact,BasisnD%tab_basisPrimSG,tab_l,nDind_DPG,iq,mole,err_sub)

     CALL get_d0GG(Qact,para_Op%para_Tnum,mole,d0GG=GGiq(iq,:,:),        &
                                         def=.TRUE.,Jac=Jac(iq),Rho=Rho)
     sqRhoOVERJac(iq) = sqrt(Rho/Jac(iq))

   END DO
!$OMP   END DO
   CALL dealloc_NParray(Qact,'Qact',name_sub)
!$OMP   END PARALLEL

   CALL dealloc_nDindex(nDind_DPG)


   DO itab=1,size(PsiR)

     CALL alloc_NParray(OpPsiR,(/nq/),'OpPsiR',name_sub)

     ! partial B to G
     CALL BDP_TO_GDP_OF_SmolyakRep(PsiR(itab)%V,BasisnD%tab_basisPrimSG,&
                                   tab_l,tab_nq,tab_nb)
     ! now PsiR is on the grid
     !write(6,*) ' R Grid of PsiR(itab) :',PsiR(itab)%V
     !write(6,*) ' R Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(6)

     ! multiplication by sqRhoOVERJac
     PsiR(itab)%V(:) = PsiR(itab)%V(:) * sqRhoOVERJac(:)

     !write(6,*) ' R*sq Grid of PsiR(itab) :',PsiR
     !write(6,*) ' R*sq ... Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(6)


     ! derivative with respect to Qj
     DO j=1,mole%nb_act1
       derive_termQdyn(:) = (/ mole%liste_QactTOQsym(j),0 /)

       PsiRj(:,j) = PsiR(itab)%V(:)
       CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiRj(:,j),BasisnD%tab_basisPrimSG, &
                                         tab_l,tab_nq,derive_termQdyn)

       !write(6,*) 'dQj Grid of PsiR(itab)',j,PsiRj(:,j) ; flush(6)
     END DO
     !write(6,*) ' dQj R*sq ... Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(6)

     OpPsiR(:) = ZERO
     DO i=1,mole%nb_act1

       PsiRi(:) = ZERO
       DO j=1,mole%nb_act1
         PsiRi(:) = PsiRi(:) + GGiq(:,j,i) * PsiRj(:,j)
       END DO
       PsiRi(:) = PsiRi(:) * Jac(:)
       !write(6,*) ' Jac Gij* ... Grid of SRep : done',iG,itab,i,OMP_GET_THREAD_NUM() ; flush(6)

       derive_termQdyn(:) = (/ mole%liste_QactTOQsym(i),0 /)

       CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiRi(:),BasisnD%tab_basisPrimSG,&
                                         tab_l,tab_nq,derive_termQdyn)

       OpPsiR(:) = OpPsiR(:) + PsiRi(:)

     END DO
     !write(6,*) ' dQi Jac Gij* ... Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(6)

     OpPsiR(:) = (  -HALF*OpPsiR(:)/Jac(:) + PsiR(itab)%V(:)*V(:) ) / sqRhoOVERJac(:)
     !write(6,*) ' OpPsiR Grid of PsiR(itab)',OpPsiR
     !write(6,*) ' OpPsiR Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(6)

     !tranfert OpPsiR (on the grid) to PsiR(itab)
     PsiR(itab)%V(:) = OpPsiR(:)



     CALL dealloc_NParray(OpPsiR,'OpPsiR', name_sub)


     CALL GDP_TO_BDP_OF_SmolyakRep(PsiR(itab)%V,BasisnD%tab_basisPrimSG,&
                                   tab_l,tab_nq,tab_nb)
     ! now PsiR(itab)%V is on the Basis


   END DO

   IF (allocated(PsiRi))        CALL dealloc_NParray(PsiRi,       'PsiRi',       name_sub)
   IF (allocated(PsiRj))        CALL dealloc_NParray(PsiRj,       'PsiRj',       name_sub)
   IF (allocated(GGiq))         CALL dealloc_NParray(GGiq,        'GGiq',        name_sub)
   IF (allocated(tab_nb))       CALL dealloc_NParray(tab_nb,      'tab_nb',      name_sub)
   IF (allocated(tab_nq))       CALL dealloc_NParray(tab_nq,      'tab_nq',      name_sub)
   IF (allocated(tab_l))        CALL dealloc_NParray(tab_l,       'tab_l',       name_sub)
   IF (allocated(sqRhoOVERJac)) CALL dealloc_NParray(sqRhoOVERJac,'sqRhoOVERJac',name_sub)
   IF (allocated(Jac))          CALL dealloc_NParray(Jac,         'Jac',         name_sub)
   IF (allocated(V))            CALL dealloc_NParray(V,           'V',           name_sub)

   !write(6,*) '============ END ===============' ; flush(6)
   !write(6,*) '================================' ; flush(6)

  !-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'END ',name_sub
    CALL flush_perso(out_unitp)
  END IF
  !-----------------------------------------------------------
  !CALL UnCheck_mem() ; stop
  END SUBROUTINE sub_TabOpPsi_OF_ONEDP_FOR_SGtype4_ompGrid

  SUBROUTINE sub_OpPsi_OF_ONEDP_FOR_SGtype4_ompGrid(PsiR,iG,para_Op)
  USE mod_system
  USE mod_nDindex

  USE mod_Coord_KEO,               ONLY : zmatrix, get_Qact, get_d0GG

  USE mod_basis_set_alloc,         ONLY : basis
  USE mod_basis,                   ONLY : Rec_Qact_SG4
  USE mod_basis_BtoG_GtoB_SGType4, ONLY : TypeRVec,                     &
                     BDP_TO_GDP_OF_SmolyakRep,GDP_TO_BDP_OF_SmolyakRep, &
                   DerivOp_TO_RDP_OF_SmolaykRep,tabR2grid_TO_tabR1_AT_iG,&
                     getbis_tab_nq,getbis_tab_nb
  USE mod_Op,                      ONLY : param_Op,write_param_Op

  IMPLICIT NONE

  TYPE (TypeRVec),                    intent(inout)    :: PsiR
  integer,                            intent(in)       :: iG

  TYPE (param_Op),                    intent(inout)    :: para_Op


  !local variables
  TYPE (zmatrix), pointer :: mole
  TYPE(basis),    pointer :: BasisnD

  integer :: iq,nb,nq,i,j,D

  integer                               :: derive_termQdyn(2)

  real (kind=Rkind),  allocatable       :: V(:)
  real (kind=Rkind),  allocatable       :: Qact(:)
  real (kind=Rkind),  allocatable       :: GGiq(:,:,:)
  TYPE (Type_nDindex)                   :: nDind_DPG    ! multidimensional DP index
  real (kind=Rkind),  allocatable       :: sqRhoOVERJac(:),Jac(:)
  real (kind=Rkind)                     :: Rho

  real (kind=Rkind),  allocatable       :: PsiRj(:,:),PsiRi(:),OpPsiR(:)

  integer,            allocatable       :: tab_nq(:)
  integer,            allocatable       :: tab_nb(:)
  integer,            allocatable       :: tab_l(:)
  integer :: iterm00
  integer :: err_sub,nb_thread


  !----- for debuging ----------------------------------------------
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub='sub_OpPsi_OF_ONEDP_FOR_SGtype4_ompGrid'
  logical, parameter :: debug = .FALSE.
  !logical, parameter :: debug = .TRUE.
  !-----------------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*) 'nb_bie,nb_baie',para_Op%nb_bie,para_Op%nb_baie
    write(out_unitp,*) 'nb_act1',para_Op%mole%nb_act1
    write(out_unitp,*) 'nb_var',para_Op%mole%nb_var
    CALL flush_perso(out_unitp)
  END IF
  !-----------------------------------------------------------------
  !CALL Check_mem()

   !write(6,*) '================================' ; flush(6)
   !write(6,*) '============ START =============' ; flush(6)


  IF (Grid_omp == 0 .AND. OpPsi_omp > 0) THEN
    nb_thread = 1
  ELSE
    nb_thread = Grid_maxth
  END IF


   mole    => para_Op%mole
   BasisnD => para_Op%BasisnD

   D = size(BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval,dim=1)

   CALL alloc_NParray(tab_l ,(/ D /),'tab_l', name_sub)
   CALL alloc_NParray(tab_nb,(/ D /),'tab_nb',name_sub)
   CALL alloc_NParray(tab_nq,(/ D /),'tab_nq',name_sub)

   tab_l(:)  = BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval(:,iG)
   tab_nq(:) = getbis_tab_nq(tab_l,BasisnD%tab_basisPrimSG)
   tab_nb(:) = getbis_tab_nb(tab_l,BasisnD%tab_basisPrimSG)


   nb = size(PsiR%V) ! because PsiR is in basis Rep
   nq = BasisnD%para_SGType2%tab_nq_OF_SRep(iG)

   CALL alloc_NParray(PsiRj,(/ nq,mole%nb_act1 /),'PsiRj',       name_sub)
   CALL alloc_NParray(PsiRi,       (/ nq /),      'PsiRi',       name_sub)
   CALL alloc_NParray(sqRhoOVERJac,(/ nq /),      'sqRhoOVERJac',name_sub)
   CALL alloc_NParray(Jac,         (/ nq /),      'Jac',         name_sub)

   !transfert part of the potential
   CALL get_OpGrid_OF_ONEDP_FOR_SG4(iG,tab_l,para_Op,V)


   ! G calculation
   CALL alloc_NParray(GGiq,(/nq,mole%nb_act1,mole%nb_act1/),'GGiq',name_sub)
   CALL init_nDindexPrim(nDind_DPG,BasisnD%nb_basis,tab_nq,type_OF_nDindex=-1)

 !$OMP parallel                                               &
 !$OMP default(none)                                          &
 !$OMP shared(mole,BasisnD,para_Op,nq,tab_l,nDind_DPG,iG)     &
 !$OMP shared(sqRhoOVERJac,Jac,GGiq)                          &
 !$OMP private(iq,Qact,err_sub,Rho)                           &
 !$OMP num_threads(nb_thread)
   CALL alloc_NParray(Qact,(/mole%nb_var/),'Qact',name_sub)
 !$OMP   DO SCHEDULE(STATIC)
   DO iq=1,nq
     CALL get_Qact(Qact,mole%ActiveTransfo) ! rigid, flexible coordinates
     CALL Rec_Qact_SG4(Qact,BasisnD%tab_basisPrimSG,tab_l,nDind_DPG,iq,mole,err_sub)

     CALL get_d0GG(Qact,para_Op%para_Tnum,mole,d0GG=GGiq(iq,:,:),        &
                                         def=.TRUE.,Jac=Jac(iq),Rho=Rho)
     sqRhoOVERJac(iq) = sqrt(Rho/Jac(iq))

   END DO
!$OMP   END DO
   CALL dealloc_NParray(Qact,'Qact',name_sub)
!$OMP   END PARALLEL

   CALL dealloc_nDindex(nDind_DPG)

   CALL alloc_NParray(OpPsiR,(/nq/),'OpPsiR',name_sub)

   ! partial B to G
   CALL BDP_TO_GDP_OF_SmolyakRep(PsiR%V,BasisnD%tab_basisPrimSG,tab_l,tab_nq,tab_nb)
   ! now PsiR is on the grid
   !write(6,*) ' R Grid of PsiR :',PsiR%V
   !write(6,*) ' R Grid of PsiR : done',OMP_GET_THREAD_NUM() ; flush(6)

   ! multiplication by sqRhoOVERJac
   PsiR%V(:) = PsiR%V(:) * sqRhoOVERJac(:)

   !write(6,*) ' R*sq Grid of PsiR :',PsiR
   !write(6,*) ' R*sq ... Grid of PsiR : done',OMP_GET_THREAD_NUM() ; flush(6)


   ! derivative with respect to Qj
   DO j=1,mole%nb_act1
     derive_termQdyn(:) = (/ mole%liste_QactTOQsym(j),0 /)

     PsiRj(:,j) = PsiR%V(:)
     CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiRj(:,j),BasisnD%tab_basisPrimSG, &
                                         tab_l,tab_nq,derive_termQdyn)

     !write(6,*) 'dQj Grid of PsiR',j,PsiRj(:,j) ; flush(6)
   END DO
   !write(6,*) ' dQj R*sq ... Grid of PsiR : done',OMP_GET_THREAD_NUM() ; flush(6)

   OpPsiR(:) = ZERO
   DO i=1,mole%nb_act1

     PsiRi(:) = ZERO
     DO j=1,mole%nb_act1
       PsiRi(:) = PsiRi(:) + GGiq(:,j,i) * PsiRj(:,j)
     END DO
     PsiRi(:) = PsiRi(:) * Jac(:)
     !write(6,*) ' Jac Gij* ... Grid of SRep : done',iG,i,OMP_GET_THREAD_NUM() ; flush(6)

     derive_termQdyn(:) = (/ mole%liste_QactTOQsym(i),0 /)

     CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiRi(:),BasisnD%tab_basisPrimSG,&
                                       tab_l,tab_nq,derive_termQdyn)

     OpPsiR(:) = OpPsiR(:) + PsiRi(:)

   END DO
   !write(6,*) ' dQi Jac Gij* ... Grid of PsiR : done',OMP_GET_THREAD_NUM() ; flush(6)

   OpPsiR(:) = (  -HALF*OpPsiR(:)/Jac(:) + PsiR%V(:)*V(:) ) / sqRhoOVERJac(:)
   !write(6,*) ' OpPsiR Grid of PsiR',OpPsiR
   !write(6,*) ' OpPsiR Grid of PsiR : done',OMP_GET_THREAD_NUM() ; flush(6)

   !tranfert OpPsiR (on the grid) to PsiR
   PsiR%V(:) = OpPsiR(:)



   CALL dealloc_NParray(OpPsiR,'OpPsiR', name_sub)


   CALL GDP_TO_BDP_OF_SmolyakRep(PsiR%V,BasisnD%tab_basisPrimSG,&
                                   tab_l,tab_nq,tab_nb)
   ! now PsiR%V is on the Basis

   IF (allocated(PsiRi))        CALL dealloc_NParray(PsiRi,       'PsiRi',       name_sub)
   IF (allocated(PsiRj))        CALL dealloc_NParray(PsiRj,       'PsiRj',       name_sub)
   IF (allocated(GGiq))         CALL dealloc_NParray(GGiq,        'GGiq',        name_sub)
   IF (allocated(tab_nb))       CALL dealloc_NParray(tab_nb,      'tab_nb',      name_sub)
   IF (allocated(tab_nq))       CALL dealloc_NParray(tab_nq,      'tab_nq',      name_sub)
   IF (allocated(tab_l))        CALL dealloc_NParray(tab_l,       'tab_l',       name_sub)
   IF (allocated(sqRhoOVERJac)) CALL dealloc_NParray(sqRhoOVERJac,'sqRhoOVERJac',name_sub)
   IF (allocated(Jac))          CALL dealloc_NParray(Jac,         'Jac',         name_sub)
   IF (allocated(V))            CALL dealloc_NParray(V,           'V',           name_sub)

   !write(6,*) '============ END ===============' ; flush(6)
   !write(6,*) '================================' ; flush(6)

  !-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'END ',name_sub
    CALL flush_perso(out_unitp)
  END IF
  !-----------------------------------------------------------
  !CALL UnCheck_mem() ; stop
  END SUBROUTINE sub_OpPsi_OF_ONEDP_FOR_SGtype4_ompGrid


  SUBROUTINE get_OpGrid_OF_ONEDP_FOR_SG4(iG,tab_l,para_Op,V,GGiq,sqRhoOVERJac,Jac)
  USE mod_system
  USE mod_nDindex

  USE mod_Coord_KEO,               ONLY : zmatrix, get_Qact, get_d0GG
  use mod_PrimOp,                  only: param_d0matop, init_d0matop,   &
                                         param_typeop, set_d0matop_at_qact, &
                                         dealloc_tab_of_d0matop

  USE mod_basis_set_alloc,         ONLY : basis
  USE mod_basis,                   ONLY : Rec_Qact_SG4_with_Tab_iq
  USE mod_basis_BtoG_GtoB_SGType4, ONLY : TypeRVec,                     &
                     BDP_TO_GDP_OF_SmolyakRep,GDP_TO_BDP_OF_SmolyakRep, &
                  DerivOp_TO_RDP_OF_SmolaykRep,tabR2grid_TO_tabR1_AT_iG,&
                     getbis_tab_nq,getbis_tab_nb
  USE mod_Op,                      ONLY : param_Op,write_param_Op


  IMPLICIT NONE

  real (kind=Rkind),  allocatable,    intent(inout)    :: V(:)
  real (kind=Rkind),  allocatable,    intent(inout), optional    :: GGiq(:,:,:)
  real (kind=Rkind),  allocatable,    intent(inout), optional    :: sqRhoOVERJac(:),Jac(:)
  integer,                            intent(in)       :: iG,tab_l(:)

  TYPE (param_Op),                    intent(inout)    :: para_Op

  !local variables
  TYPE (zmatrix), pointer :: mole
  TYPE(basis),    pointer :: BasisnD

  integer :: iq,nq,D,iOp

  integer                               :: derive_termQdyn(2)
  real (kind=Rkind),  allocatable       :: Qact(:)
  real (kind=Rkind)                     :: Rho

  integer,            allocatable       :: tab_nq(:)
  integer,            allocatable       :: tab_iq(:)
  integer :: iterm00,itabR,nR
  integer :: err_sub
  logical :: KEO
  logical :: lformatted=.TRUE.

  TYPE (param_d0MatOp), allocatable :: d0MatOp(:)
  character (len=Line_len) :: FileName_RV


  !----- for debuging ----------------------------------------------
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub='get_OpGrid_OF_ONEDP_FOR_SG4'
  logical, parameter :: debug = .FALSE.
  !logical, parameter :: debug = .TRUE.
  !-----------------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*) 'nb_bie,nb_baie',para_Op%nb_bie,para_Op%nb_baie
    write(out_unitp,*) 'nb_act1',para_Op%mole%nb_act1
    write(out_unitp,*) 'nb_var',para_Op%mole%nb_var
    CALL flush_perso(out_unitp)
  END IF
  !-----------------------------------------------------------------

   !write(6,*) '================================' ; flush(6)
   !write(6,*) '============ START =============' ; flush(6)

   mole    => para_Op%mole
   BasisnD => para_Op%BasisnD

   KEO = present(GGiq) .AND.present(sqRhoOVERJac) .AND. present(Jac) .AND. para_Op%n_op == 0
   D = size(tab_l)

   CALL alloc_NParray(Qact,(/mole%nb_var/),'Qact',name_sub)
   CALL alloc_NParray(tab_nq,(/ D /),'tab_nq',name_sub)
   CALL alloc_NParray(tab_iq,(/ D /),'tab_iq',name_sub)

   tab_nq(:) = getbis_tab_nq(tab_l,BasisnD%tab_basisPrimSG)

   nq = BasisnD%para_SGType2%tab_nq_OF_SRep(iG)

   IF (KEO) THEN
     CALL alloc_NParray(sqRhoOVERJac,(/ nq /),      'sqRhoOVERJac',name_sub)
     CALL alloc_NParray(Jac,         (/ nq /),      'Jac',         name_sub)
   END IF

   !transfert part of the scalar part of the potential
   iterm00 = para_Op%derive_term_TO_iterm(0,0)

   lformatted = para_Op%OpGrid(iterm00)%para_FileGrid%Formatted_FileGrid

   IF (debug) THEN
     write(out_unitp,*) iG,'Save_MemGrid,Save_MemGrid_done',            &
                 para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid,    &
                 para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_done
     write(out_unitp,*) iG,'Save_FileGrid,Save_FileGrid_done',          &
                para_Op%OpGrid(iterm00)%para_FileGrid%Save_FileGrid,    &
                para_Op%OpGrid(iterm00)%para_FileGrid%Save_FileGrid_done
   END IF

   IF (para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_done) THEN
     IF (para_Op%OpGrid(iterm00)%grid_zero) THEN
       CALL alloc_NParray(V,(/ nq /),'V',name_sub)
       V(:) = ZERO
     ELSE
       IF (para_Op%OpGrid(iterm00)%grid_cte) THEN
         CALL alloc_NParray(V,(/ nq /),'V',name_sub)
         V(:) = para_Op%OpGrid(iterm00)%Mat_cte(1,1)
       ELSE
         IF (associated(para_Op%OpGrid)) THEN
         IF (associated(para_Op%OpGrid(iterm00)%Grid)) THEN

           CALL tabR2grid_TO_tabR1_AT_iG(V,para_Op%OpGrid(iterm00)%Grid(:,1,1),&
                                         iG,BasisnD%para_SGType2)
         END IF
         END IF
       END IF
     END IF
   ELSE IF (para_Op%OpGrid(iterm00)%para_FileGrid%Save_FileGrid_done .OR. &
            para_Op%OpGrid(iterm00)%para_FileGrid%Read_FileGrid) THEN
    !$OMP  CRITICAL (get_OpGrid_OF_ONEDP_FOR_SG4_CRIT2)
     FileName_RV = trim(para_Op%OpGrid(iterm00)%file_Grid%name) // '_SGterm' // int_TO_char(iG)

     IF (iG == 1 .AND. debug) THEN
       write(out_unitp,*) 'iG,nq, FileName_RV',iG,nq,FileName_RV
       CALL flush_perso(out_unitp)
     END IF

     CALL sub_ReadRV(V,FileName_RV,lformatted,err_sub)
     IF (err_sub /= 0) STOP 'error while reading the grid'
     IF (nq /= size(V)) THEN
       write(out_unitp,*) ' ERROR in ',name_sub
       write(out_unitp,*) ' the size of the subroutine and file grids are different'
       write(out_unitp,*) ' nq (subroutine) and nq (file):',nq,size(V)
       write(out_unitp,*) ' file name: ',FileName_RV
       STOP
     END IF
     IF (iG == 1 .AND. debug) THEN
       write(out_unitp,*) 'iG,nq, read V',iG,nq
       CALL flush_perso(out_unitp)
     END IF

    !$OMP  END CRITICAL (get_OpGrid_OF_ONEDP_FOR_SG4_CRIT2)

   ELSE IF (.NOT. para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_done) THEN
     CALL alloc_NParray(V,(/ nq /),'V',name_sub)
     allocate(d0MatOp(para_Op%para_PES%nb_scalar_Op+2))
     DO iOp=1,size(d0MatOp)
       CALL Init_d0MatOp(d0MatOp(iOp),para_Op%param_TypeOp,para_Op%para_PES%nb_elec)
     END DO
   END IF

   ! G calculation
   CALL alloc_NParray(GGiq,(/nq,mole%nb_act1,mole%nb_act1/),'GGiq',name_sub)

   tab_iq(:) = 1 ; tab_iq(1) = 0

   DO iq=1,nq

     CALL ADD_ONE_TO_nDval_m1(tab_iq,tab_nq)

     CALL get_Qact(Qact,mole%ActiveTransfo) ! rigid, flexible coordinates
     CALL Rec_Qact_SG4_with_Tab_iq(Qact,BasisnD%tab_basisPrimSG,tab_l,tab_iq,mole,err_sub)

     IF (KEO) THEN
       CALL get_d0GG(Qact,para_Op%para_Tnum,mole,d0GG=GGiq(iq,:,:),     &
                                         def=.TRUE.,Jac=Jac(iq),Rho=Rho)

       sqRhoOVERJac(iq) = sqrt(Rho/Jac(iq))
     END IF

     IF (.NOT. para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_done  .AND. &
         .NOT. para_Op%OpGrid(iterm00)%para_FileGrid%Save_FileGrid_done .AND. &
         .NOT. para_Op%OpGrid(iterm00)%para_FileGrid%Read_FileGrid) THEN

        CALL Set_d0MatOp_AT_Qact(Qact,d0MatOp,mole,                     &
                                 para_Op%para_Tnum,para_Op%para_PES)
        V(iq) = d0MatOp(iterm00)%ReVal(1,1,1)

        IF (iG == 1 .AND. debug) THEN
          write(out_unitp,*) 'iG,nq,V',iG,nq,V
          CALL flush_perso(out_unitp)
        END IF

     END IF

   END DO


   IF (para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid .AND.       &
     .NOT. para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_done) THEN
   IF (associated(para_Op%OpGrid(iterm00)%Grid)) THEN
     itabR = BasisnD%para_SGType2%tab_Sum_nq_OF_SRep(iG)
     nR    = BasisnD%para_SGType2%tab_nq_OF_SRep(iG)
     para_Op%OpGrid(iterm00)%Grid(itabR-nR+1:itabR,1,1) = V

     IF (iG == 1 .AND. debug) THEN
       write(out_unitp,*) 'iG,nq, save mem V',iG,nq
       CALL flush_perso(out_unitp)
     END IF

   END IF
   END IF

   IF (para_Op%OpGrid(iterm00)%para_FileGrid%Save_FileGrid .AND.           &
      .NOT. para_Op%OpGrid(iterm00)%para_FileGrid%Save_FileGrid_done .AND. &
      .NOT. para_Op%OpGrid(iterm00)%para_FileGrid%Read_FileGrid) THEN

!    !$OMP  CRITICAL (get_OpGrid_OF_ONEDP_FOR_SG4_CRIT1)
     FileName_RV = trim(para_Op%OpGrid(iterm00)%file_Grid%name) // '_SGterm' // int_TO_char(iG)
     CALL sub_WriteRV(V,FileName_RV,lformatted,err_sub)

     IF (iG == 1 .AND. debug) THEN
       write(out_unitp,*) 'iG,nq, save file V',iG,nq
       CALL flush_perso(out_unitp)
     END IF

!    !$OMP  END CRITICAL (get_OpGrid_OF_ONEDP_FOR_SG4_CRIT1)
   END IF

  IF (allocated(tab_nq))       CALL dealloc_NParray(tab_nq,      'tab_nq',      name_sub)
  IF (allocated(tab_iq))       CALL dealloc_NParray(tab_iq,      'tab_iq',      name_sub)
  IF (allocated(Qact))         CALL dealloc_NParray(Qact,        'Qact',        name_sub)

  IF (allocated(d0MatOp)) THEN
    CALL dealloc_Tab_OF_d0MatOp(d0MatOp)
    deallocate(d0MatOp)
  END IF

   !write(6,*) '============ END ===============' ; flush(6)
   !write(6,*) '================================' ; flush(6)

  !-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'END ',name_sub
    CALL flush_perso(out_unitp)
  END IF
  !-----------------------------------------------------------
  END SUBROUTINE get_OpGrid_OF_ONEDP_FOR_SG4

END MODULE mod_OpPsi_SG4
