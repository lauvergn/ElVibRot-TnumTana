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
MODULE mod_OpPsi_SG4

PRIVATE
PUBLIC :: sub_OpPsi_FOR_SGtype4,sub_TabOpPsi_FOR_SGtype4
PUBLIC :: sub_TabOpPsi_OF_ONEDP_FOR_SGtype4
PUBLIC :: sub_TabOpPsi_OF_ONEGDP_WithOp_FOR_SGtype4
PUBLIC :: sub_OpPsi_OF_ONEDP_FOR_SGtype4
PUBLIC :: get_OpGrid_type0_OF_ONEDP_FOR_SG4
PUBLIC :: get_OpGrid_type10_OF_ONEDP_FOR_SG4

PUBLIC :: get_OpGrid_type1_OF_ONEDP_FOR_SG4_new2

CONTAINS

 SUBROUTINE sub_OpPsi_FOR_SGtype4(Psi,OpPsi,para_Op)

 USE mod_system
!$ USE omp_lib, only : OMP_GET_THREAD_NUM
 USE mod_nDindex
 USE mod_Coord_KEO,                ONLY : CoordType

 USE mod_basis_set_alloc,          ONLY : basis
 USE mod_basis_RCVec_SGType4,      ONLY : TypeRVec,dealloc_TypeRVec
 USE mod_basis_BtoG_GtoB_SGType4,  ONLY : tabPackedBasis_TO_tabR_AT_iG, &
                                          tabR_AT_iG_TO_tabPackedBasis

 USE mod_psi,                      ONLY : param_psi,ecri_psi

 USE mod_SetOp,                    ONLY : param_Op,write_param_Op
 IMPLICIT NONE

 TYPE (param_psi), intent(in)      :: Psi
 TYPE (param_psi), intent(inout)   :: OpPsi

 TYPE (param_Op),  intent(inout)   :: para_Op

 ! local variables
 TYPE (CoordType), pointer :: mole
 TYPE (basis),   pointer :: BasisnD

 TYPE (TypeRVec)         :: PsiR

 integer                :: ib,i,iG,iiG,nb_thread,ith,iterm00
 integer, allocatable   :: tab_l(:)
 logical                :: not_init

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

  IF (Psi%cplx) STOP 'cplx in sub_OpPsi_FOR_SGtype4'
  IF (para_Op%nb_bie /= 1) STOP 'nb_bie /= 1 in sub_OpPsi_FOR_SGtype4'

  IF (SG4_omp == 0) THEN
    nb_thread = 1
  ELSE
    nb_thread = SG4_maxth
  END IF


 nb_mult_OpPsi   = 0
 OpPsi           = Psi ! for the allocation. It has to be changed!
 OpPsi%RvecB(:)  = ZERO

 IF (print_level > 0 .AND. BasisnD%para_SGType2%nb_SG > 10**4 ) THEN
   write(out_unitp,'(a)')              'OpPsi SG4 (%): [-10-20-30-40-50-60-70-80-90-100]'
   write(out_unitp,'(a)',ADVANCE='no') 'OpPsi SG4 (%): ['
   CALL flush_perso(out_unitp)
 END IF

 IF (OpPsi_omp == 0) THEN
   DO iG=1,BasisnD%para_SGType2%nb_SG

     !write(out_unitp,*) 'iG',iG
     !transfert part of the psi%RvecB(:) to PsiR%V
     CALL tabPackedBasis_TO_tabR_AT_iG(PsiR%V,psi%RvecB,    &
                                         iG,BasisnD%para_SGType2)
     !write(out_unitp,*) 'iG,PsiR',iG,PsiR%V

     CALL sub_OpPsi_OF_ONEDP_FOR_SGtype4_ompGrid(PsiR,iG,para_Op)

     !transfert back, PsiR%V (HPsi) to the psi%RvecB(:)
     CALL tabR_AT_iG_TO_tabPackedBasis(OpPsi%RvecB,PsiR%V,  &
                     iG,BasisnD%para_SGType2,BasisnD%WeightSG(iG))

     !deallocate PsiR
     CALL dealloc_TypeRVec(PsiR)


     IF (print_level > 0  .AND. BasisnD%para_SGType2%nb_SG > 10**4 .AND. &
         mod(iG,max(1,BasisnD%para_SGType2%nb_SG/10)) == 0) THEN
       write(out_unitp,'(a)',ADVANCE='no') '---'
       CALL flush_perso(out_unitp)
     END IF

     !write(out_unitp,*) 'iG done:',iG ; flush(out_unitp)
   END DO
 ELSE

   !to be sure to have the correct number of threads
   nb_thread = BasisnD%para_SGType2%nb_threads

   !$OMP parallel                                                &
   !$OMP default(none)                                           &
   !$OMP shared(Psi,OpPsi)                                       &
   !$OMP shared(para_Op,BasisnD,print_level,out_unitp,MPI_id)    &
   !$OMP private(iG,iiG,tab_l,PsiR,ith)                          &
   !$OMP num_threads(BasisnD%para_SGType2%nb_threads)

   !--------------------------------------------------------------
   !-- For the initialization of tab_l(:) and the use of ADD_ONE_TO_nDindex in the parallel loop
   CALL alloc_NParray(tab_l,(/ BasisnD%para_SGType2%nDind_SmolyakRep%ndim /),'tabl_l',name_sub)

   ith = 0
   !$ ith = OMP_GET_THREAD_NUM()
   tab_l(:) = BasisnD%para_SGType2%nDval_init(:,ith+1)
   !--------------------------------------------------------------

   ! we are not using the parallel do, to be able to use the correct initialized tab_l with nDval_init
   DO iG=BasisnD%para_SGType2%iG_th(ith+1),BasisnD%para_SGType2%fG_th(ith+1)

     CALL ADD_ONE_TO_nDindex(BasisnD%para_SGType2%nDind_SmolyakRep,tab_l,iG=iG)

     !write(out_unitp,*) 'iG',iG ; flush(out_unitp)
     !transfert part of the psi%RvecB(:) to PsiR%V
     CALL tabPackedBasis_TO_tabR_AT_iG(PsiR%V,psi%RvecB,iG,BasisnD%para_SGType2)
     !write(out_unitp,*) 'iG,PsiR',iG,PsiR%V ; flush(out_unitp)

     CALL sub_OpPsi_OF_ONEDP_FOR_SGtype4(PsiR,iG,tab_l,para_Op)

     !transfert back, PsiR%V (HPsi) to the psi%RvecB(:)
     CALL tabR_AT_iG_TO_tabPackedBasis(OpPsi%RvecB,PsiR%V,              &
                           iG,BasisnD%para_SGType2,BasisnD%WeightSG(iG))

     !deallocate PsiR(:)
     CALL dealloc_TypeRVec(PsiR)


     IF (print_level > 0  .AND. BasisnD%para_SGType2%nb_SG > 10**4 .AND. &
         mod(iG,max(1,BasisnD%para_SGType2%nb_SG/10)) == 0) THEN
       write(out_unitp,'(a)',ADVANCE='no') '---'
       CALL flush_perso(out_unitp)
     END IF

     !write(out_unitp,*) 'iG done:',iG ; flush(out_unitp)
   END DO
   CALL dealloc_NParray(tab_l,'tabl_l',name_sub)
  !$OMP   END PARALLEL
END IF

 IF (print_level > 0 .AND. BasisnD%para_SGType2%nb_SG > 10**4) THEN
   write(out_unitp,'(a)',ADVANCE='yes') '-]'
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

 USE mod_Coord_KEO,               ONLY : CoordType

 USE mod_basis_set_alloc,         ONLY : basis
 !USE mod_basis,                   ONLY : Rec_Qact_SG4_with_Tab_iq
 USE mod_basis_RCVec_SGType4,     ONLY : TypeRVec
 USE mod_basis_BtoG_GtoB_SGType4, ONLY : BDP_TO_GDP_OF_SmolyakRep,    &
                                            GDP_TO_BDP_OF_SmolyakRep, &
                                        DerivOp_TO_RDP_OF_SmolaykRep, &
                                        tabR2grid_TO_tabR1_AT_iG,     &
                                        getbis_tab_nq,getbis_tab_nb
 USE mod_SetOp,                    ONLY : param_Op,write_param_Op
 IMPLICIT NONE

 TYPE (TypeRVec),                    intent(inout)    :: PsiR
 integer,                            intent(in)       :: iG,tab_l(:)

 TYPE (param_Op),                    intent(inout)       :: para_Op

 !local variables
 TYPE (CoordType), pointer :: mole
 TYPE(basis),    pointer :: BasisnD

 integer :: iq,nb,nq,i,j,D
 integer :: ib0,jb0,nb0,iqi,iqf,jqi,jqf

 integer                               :: derive_termQdyn(2)

 real (kind=Rkind),  allocatable       :: V(:,:,:)
 real (kind=Rkind),  allocatable       :: Qact(:)
 real (kind=Rkind),  allocatable       :: GGiq(:,:,:)
 real (kind=Rkind),  allocatable       :: sqRhoOVERJac(:),Jac(:)
 real (kind=Rkind)                     :: Rho

 real (kind=Rkind),  allocatable       :: PsiRj(:,:),PsiRi(:),OpPsiR(:)
 real (kind=Rkind),  allocatable       :: VPsi(:,:) ! size(nq,nb0)

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

  !write(out_unitp,*) '================================' ; flush(out_unitp)
  !write(out_unitp,*) '============ START =============' ; flush(out_unitp)

  mole    => para_Op%mole
  BasisnD => para_Op%BasisnD

 D = size(tab_l)

 CALL alloc_NParray(Qact,(/mole%nb_var/),'Qact',name_sub)
 CALL alloc_NParray(tab_nb,(/ D /),'tab_nb',name_sub)
 CALL alloc_NParray(tab_nq,(/ D /),'tab_nq',name_sub)
 CALL alloc_NParray(tab_iq,(/ D /),'tab_iq',name_sub)

 tab_nq(:) = getbis_tab_nq(tab_l,BasisnD%tab_basisPrimSG)
 tab_nb(:) = getbis_tab_nb(tab_l,BasisnD%tab_basisPrimSG)


 nb  = BasisnD%para_SGType2%tab_nb_OF_SRep(iG)
 nq  = BasisnD%para_SGType2%tab_nq_OF_SRep(iG)
 nb0 = BasisnD%para_SGType2%nb0

 CALL alloc_NParray(PsiRj,(/ nq,mole%nb_act1 /),'PsiRj',       name_sub)
 CALL alloc_NParray(PsiRi,       (/ nq /),      'PsiRi',       name_sub)


 CALL get_OpGrid_type10_OF_ONEDP_FOR_SG4(iG,tab_l,para_Op,V,GGiq,sqRhoOVERJac,Jac)

 !write(out_unitp,*) ' GGiq:',GGiq
 !write(out_unitp,*) ' Jac :',Jac
 !write(out_unitp,*) ' sqRhoOVERJac :',sqRhoOVERJac
 !write(out_unitp,*) 'd0GG ... done' ; flush(out_unitp)
 CALL alloc_NParray(OpPsiR,(/nq/),'OpPsiR',name_sub)

 !write(out_unitp,*) ' V Basis of PsiR :',PsiR%V

 ! partial B to G
 CALL BDP_TO_GDP_OF_SmolyakRep(PsiR%V,BasisnD%tab_basisPrimSG,        &
                                 tab_l,tab_nq,tab_nb,nb0)
 ! now PsiR is on the grid
 !write(out_unitp,*) ' R Grid of PsiR :',PsiR%V
 !write(out_unitp,*) ' R Grid of PsiR : done',OMP_GET_THREAD_NUM() ; flush(out_unitp)

 ! VPsi calculation, it has to be done before because V is not diagonal
 CALL alloc_NParray(VPsi,(/nq,nb0/),'VPsi',name_sub)
 VPsi(:,:) = ZERO
 DO ib0=1,nb0
   jqi = 1
   jqf = nq
   DO jb0=1,nb0
     VPsi(:,ib0) = VPsi(:,ib0) + V(:,ib0,jb0) * PsiR%V(jqi:jqf)
     jqi = jqf+1
     jqf = jqf+nq
   END DO
 END DO

 iqi = 1
 iqf = nq
 DO ib0=1,nb0

   ! multiplication by sqRhoOVERJac
   PsiR%V(iqi:iqf) = PsiR%V(iqi:iqf) * sqRhoOVERJac(:)

   !write(out_unitp,*) ' R*sq Grid of PsiR :',PsiR%V
   !write(out_unitp,*) ' R*sq ... Grid of PsiR : done',OMP_GET_THREAD_NUM() ; flush(out_unitp)


   ! derivative with respect to Qj
   DO j=1,mole%nb_act1
     derive_termQdyn(:) = (/ mole%liste_QactTOQdyn(j),0 /)

     PsiRj(:,j) = PsiR%V(iqi:iqf)
     CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiRj(:,j),BasisnD%tab_basisPrimSG, &
                                         tab_l,tab_nq,derive_termQdyn)

     !write(out_unitp,*) 'dQj Grid of PsiR',j,PsiRj(:,j) ; flush(out_unitp)
   END DO
   !write(out_unitp,*) ' dQj R*sq ... Grid of PsiR : done',OMP_GET_THREAD_NUM() ; flush(out_unitp)

   OpPsiR(:) = ZERO
   DO i=1,mole%nb_act1

     PsiRi(:) = ZERO
     DO j=1,mole%nb_act1
       PsiRi(:) = PsiRi(:) + GGiq(:,j,i) * PsiRj(:,j)
     END DO
     PsiRi(:) = PsiRi(:) * Jac(:)
     !write(out_unitp,*) ' Jac Gij* ... Grid of SRep : done',iG,i,OMP_GET_THREAD_NUM() ; flush(out_unitp)

     derive_termQdyn(:) = (/ mole%liste_QactTOQdyn(i),0 /)

     CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiRi(:),BasisnD%tab_basisPrimSG,&
                                         tab_l,tab_nq,derive_termQdyn)

     !write(out_unitp,*) 'shape OpPsiR,PsiRi ',shape(OpPsiR),shape(PsiRi)
     !write(out_unitp,*) 'DerivOp_TO_RDP_OF_SmolaykRep : done',OMP_GET_THREAD_NUM() ; flush(out_unitp)

     OpPsiR(:) = OpPsiR(:) + PsiRi(:)
   END DO

   !write(out_unitp,*) ' dQi Jac Gij* ... Grid of PsiR : done',OMP_GET_THREAD_NUM() ; flush(out_unitp)

   OpPsiR(:) =   -HALF*OpPsiR(:) / ( Jac(:)*sqRhoOVERJac(:) ) + VPsi(:,ib0)
   !write(out_unitp,*) ' OpPsiR Grid of PsiR',OpPsiR
   !write(out_unitp,*) ' OpPsiR Grid of PsiR : done',OMP_GET_THREAD_NUM() ; flush(out_unitp)

   !tranfert OpPsiR (on the grid) to PsiR
   PsiR%V(iqi:iqf) = OpPsiR(:)

   iqi = iqf+1
   iqf = iqf+nq
 END DO

 IF (allocated(OpPsiR)) CALL dealloc_NParray(OpPsiR,'OpPsiR', name_sub)
 IF (allocated(VPsi))   CALL dealloc_NParray(VPsi,  'VPsi',   name_sub)

 CALL GDP_TO_BDP_OF_SmolyakRep(PsiR%V,BasisnD%tab_basisPrimSG,&
                               tab_l,tab_nq,tab_nb,nb0)
 ! now PsiR%V is on the Basis
 !write(out_unitp,*) ' V Basis of OpPsiR :',PsiR%V

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

 !write(out_unitp,*) '============ END ===============' ; flush(out_unitp)
 !write(out_unitp,*) '================================' ; flush(out_unitp)

 !-----------------------------------------------------------
 IF (debug) THEN
   write(out_unitp,*) 'END ',name_sub
   CALL flush_perso(out_unitp)
 END IF
 !-----------------------------------------------------------
 !CALL UnCheck_mem() ; stop
 END SUBROUTINE sub_OpPsi_OF_ONEDP_FOR_SGtype4
  SUBROUTINE sub_OpPsi_OF_ONEDP_FOR_SGtype4_ompGrid(PsiR,iG,para_Op)
  USE mod_system
  USE mod_nDindex

  USE mod_Coord_KEO,               ONLY : CoordType, get_d0GG

  USE mod_basis_set_alloc,         ONLY : basis
  USE mod_basis,                   ONLY : Rec_Qact_SG4
  USE mod_basis_BtoG_GtoB_SGType4, ONLY : TypeRVec,                     &
                     BDP_TO_GDP_OF_SmolyakRep,GDP_TO_BDP_OF_SmolyakRep, &
                   DerivOp_TO_RDP_OF_SmolaykRep,tabR2grid_TO_tabR1_AT_iG,&
                     getbis_tab_nq,getbis_tab_nb
  USE mod_SetOp,                      ONLY : param_Op,write_param_Op

  IMPLICIT NONE

  TYPE (TypeRVec),                    intent(inout)    :: PsiR
  integer,                            intent(in)       :: iG
  TYPE (param_Op),                    intent(inout)    :: para_Op


  !local variables
  TYPE (CoordType), pointer :: mole
  TYPE(basis),    pointer :: BasisnD

  integer :: iq,nb,nq,i,j,D

  integer                               :: derive_termQdyn(2)

  real (kind=Rkind),  allocatable       :: V(:,:,:)
  real (kind=Rkind),  allocatable       :: Qact(:)
  real (kind=Rkind),  allocatable       :: GGiq(:,:,:)
  TYPE (Type_nDindex)                   :: nDind_DPG    ! multidimensional DP index
  real (kind=Rkind),  allocatable       :: sqRhoOVERJac(:),Jac(:)
  real (kind=Rkind)                     :: Rho

  real (kind=Rkind),  allocatable       :: PsiRj(:,:),PsiRi(:),OpPsiR(:)
  real (kind=Rkind),  allocatable       :: VPsi(:,:) ! size(nq,nb0)

  integer,            allocatable       :: tab_nq(:)
  integer,            allocatable       :: tab_nb(:)
  integer,            allocatable       :: tab_l(:)
  integer :: iterm00
  integer :: ib0,jb0,nb0,iqi,iqf,jqi,jqf
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

   !write(out_unitp,*) '================================' ; flush(out_unitp)
   !write(out_unitp,*) '============ START =============' ; flush(out_unitp)


  IF (Grid_omp == 0 .AND. OpPsi_omp > 0) THEN
    nb_thread = 1
  ELSE
    nb_thread = Grid_maxth
  END IF


   mole    => para_Op%mole
   BasisnD => para_Op%BasisnD

   D   = size(BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval,dim=1)

   CALL alloc_NParray(tab_l ,(/ D /),'tab_l', name_sub)
   CALL alloc_NParray(tab_nb,(/ D /),'tab_nb',name_sub)
   CALL alloc_NParray(tab_nq,(/ D /),'tab_nq',name_sub)

   tab_l(:)  = BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval(:,iG)
   tab_nq(:) = getbis_tab_nq(tab_l,BasisnD%tab_basisPrimSG)
   tab_nb(:) = getbis_tab_nb(tab_l,BasisnD%tab_basisPrimSG)


   nb  = BasisnD%para_SGType2%tab_nb_OF_SRep(iG) ! because PsiR is in basis Rep
   nq  = BasisnD%para_SGType2%tab_nq_OF_SRep(iG)
   nb0 = BasisnD%para_SGType2%nb0

   CALL alloc_NParray(PsiRj,(/ nq,mole%nb_act1 /),'PsiRj',       name_sub)
   CALL alloc_NParray(PsiRi,       (/ nq /),      'PsiRi',       name_sub)
   CALL alloc_NParray(sqRhoOVERJac,(/ nq /),      'sqRhoOVERJac',name_sub)
   CALL alloc_NParray(Jac,         (/ nq /),      'Jac',         name_sub)

   !transfert part of the potential
   CALL get_OpGrid_type0_OF_ONEDP_FOR_SG4(iG,tab_l,para_Op,V)


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
   CALL BDP_TO_GDP_OF_SmolyakRep(PsiR%V,BasisnD%tab_basisPrimSG,  &
                                 tab_l,tab_nq,tab_nb,nb0)
   ! now PsiR is on the grid
   !write(out_unitp,*) ' R Grid of PsiR :',PsiR%V
   !write(out_unitp,*) ' R Grid of PsiR : done',OMP_GET_THREAD_NUM() ; flush(out_unitp)

 ! VPsi calculation, it has to be done before because V is not diagonal
 CALL alloc_NParray(VPsi,(/nq,nb0/),'VPsi',name_sub)
 VPsi(:,:) = ZERO
 DO ib0=1,nb0
   jqi = 1
   jqf = nq
   DO jb0=1,nb0
     VPsi(:,ib0) = VPsi(:,ib0) + V(:,ib0,jb0) * PsiR%V(jqi:jqf)
     jqi = jqf+1
     jqf = jqf+nq
   END DO
 END DO

 iqi = 1
 iqf = nq
 DO ib0=1,nb0

   ! multiplication by sqRhoOVERJac
   PsiR%V(iqi:iqf) = PsiR%V(iqi:iqf) * sqRhoOVERJac(:)

   !write(out_unitp,*) ' R*sq Grid of PsiR :',PsiR
   !write(out_unitp,*) ' R*sq ... Grid of PsiR : done',OMP_GET_THREAD_NUM() ; flush(out_unitp)


   ! derivative with respect to Qj
   DO j=1,mole%nb_act1
     derive_termQdyn(:) = (/ mole%liste_QactTOQdyn(j),0 /)

     PsiRj(:,j) = PsiR%V(iqi:iqf)
     CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiRj(:,j),BasisnD%tab_basisPrimSG, &
                                         tab_l,tab_nq,derive_termQdyn)

     !write(out_unitp,*) 'dQj Grid of PsiR',j,PsiRj(:,j) ; flush(out_unitp)
   END DO
   !write(out_unitp,*) ' dQj R*sq ... Grid of PsiR : done',OMP_GET_THREAD_NUM() ; flush(out_unitp)

   OpPsiR(:) = ZERO
   DO i=1,mole%nb_act1

     PsiRi(:) = ZERO
     DO j=1,mole%nb_act1
       PsiRi(:) = PsiRi(:) + GGiq(:,j,i) * PsiRj(:,j)
     END DO
     PsiRi(:) = PsiRi(:) * Jac(:)
     !write(out_unitp,*) ' Jac Gij* ... Grid of SRep : done',iG,i,OMP_GET_THREAD_NUM() ; flush(out_unitp)

     derive_termQdyn(:) = (/ mole%liste_QactTOQdyn(i),0 /)

     CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiRi(:),BasisnD%tab_basisPrimSG,&
                                       tab_l,tab_nq,derive_termQdyn)

     OpPsiR(:) = OpPsiR(:) + PsiRi(:)

   END DO
   !write(out_unitp,*) ' dQi Jac Gij* ... Grid of PsiR : done',OMP_GET_THREAD_NUM() ; flush(out_unitp)

   OpPsiR(:) =   -HALF*OpPsiR(:) / ( Jac(:)*sqRhoOVERJac(:) ) + VPsi(:,ib0)
   !write(out_unitp,*) ' OpPsiR Grid of PsiR',OpPsiR
   !write(out_unitp,*) ' OpPsiR Grid of PsiR : done',OMP_GET_THREAD_NUM() ; flush(out_unitp)

   !tranfert OpPsiR (on the grid) to PsiR
   PsiR%V(iqi:iqf) = OpPsiR(:)

   iqi = iqf+1
   iqf = iqf+nq
 END DO


 IF (allocated(OpPsiR)) CALL dealloc_NParray(OpPsiR,'OpPsiR', name_sub)
 IF (allocated(VPsi))   CALL dealloc_NParray(VPsi,  'VPsi',   name_sub)


 CALL GDP_TO_BDP_OF_SmolyakRep(PsiR%V,BasisnD%tab_basisPrimSG,          &
                                   tab_l,tab_nq,tab_nb,nb0)
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

 !write(out_unitp,*) '============ END ===============' ; flush(out_unitp)
 !write(out_unitp,*) '================================' ; flush(out_unitp)

 !-----------------------------------------------------------
 IF (debug) THEN
   write(out_unitp,*) 'END ',name_sub
   CALL flush_perso(out_unitp)
 END IF
 !-----------------------------------------------------------
 !CALL UnCheck_mem() ; stop
 END SUBROUTINE sub_OpPsi_OF_ONEDP_FOR_SGtype4_ompGrid

!=======================================================================================
!> OpPsi is assigned in this routine with OpPsi%Vec=0 initially
!> after this subroutine, OpPsi will be ready
!=======================================================================================
SUBROUTINE sub_TabOpPsi_FOR_SGtype4(Psi,OpPsi,para_Op)
  USE mod_system
  !$ USE omp_lib, only : OMP_GET_THREAD_NUM
  USE mod_nDindex

  USE mod_Coord_KEO,                ONLY : CoordType

  USE mod_SymAbelian,               ONLY : Calc_symab1_EOR_symab2

  USE mod_basis_set_alloc,          ONLY : basis
  USE mod_basis_BtoG_GtoB_SGType4,  ONLY : tabPackedBasis_TO_tabR_AT_iG, &
                                           tabR_AT_iG_TO_tabPackedBasis, &
                                           TypeRVec,dealloc_TypeRVec
  USE mod_psi,                      ONLY : param_psi,ecri_psi,          &
                                           Set_symab_OF_psiBasisRep

  USE mod_SetOp,                    ONLY : param_Op,write_param_Op
  IMPLICIT NONE

  TYPE (param_psi), intent(in)      :: Psi(:)
  TYPE (param_psi), intent(inout)   :: OpPsi(:)
  TYPE (param_Op),  intent(inout)   :: para_Op

  ! local variables
  TYPE (CoordType), pointer :: mole
  TYPE (basis),   pointer :: BasisnD

  TYPE (TypeRVec),   allocatable       :: PsiR(:)

  Real (kind=Rkind), allocatable       :: PsiR_temp(:)
  integer                              :: ib,i,iG,iiG,nb_thread
  integer                              :: itab,ith,iterm00,packet_size,OpPsi_symab
  integer, allocatable                 :: tab_l(:)
  logical                              :: not_init


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
  !-----------------------------------------------------------------------------
  mole    => para_Op%mole
  BasisnD => para_Op%BasisnD

  IF(MPI_id==0) THEN
    IF (size(Psi) == 0) THEN
      write(out_unitp,*) 'ERROR in ',name_sub
      write(out_unitp,*) '  The size of Psi(:) is zero!!'
      write(out_unitp,*) '  => Check the fortran.'
      STOP ' ERROR in sub_TabOpPsi_FOR_SGtype4: size(Psi) = 0'
    END IF
    IF (Psi(1)%cplx) THEN
      write(out_unitp,*) 'ERROR in ',name_sub
      write(out_unitp,*) '  Psi(1) is complex !!'
      write(out_unitp,*) '  => Check the fortran.'
      STOP ' ERROR in sub_TabOpPsi_FOR_SGtype4: Psi(1) is complex'
    END IF
  ENDIF
  !IF (para_Op%nb_bie /= 1) STOP 'nb_bie /= 1 in sub_TabOpPsi_FOR_SGtype4'

  !-----------------------------------------------------------------------------
  IF (SG4_omp == 0) THEN
    nb_thread = 1
  ELSE
    nb_thread = SG4_maxth
  END IF
  !-----------------------------------------------------------------------------

  nb_mult_OpPsi = 0
  IF(MPI_id==0) THEN
    DO itab=1,size(Psi)
      OpPsi(itab)           = Psi(itab) ! for the allocation. It has to be changed!
      OpPsi(itab)%RvecB(:)  = ZERO
    END DO
  ENDIF

  IF (print_level > 0 .AND. BasisnD%para_SGType2%nb_SG > 10**4 ) THEN
    write(out_unitp,'(a)')              'OpPsi SG4 (%): [-10-20-30-40-50-60-70-80-90-100]'
    write(out_unitp,'(a)',ADVANCE='no') 'OpPsi SG4 (%): ['
    CALL flush_perso(out_unitp)
  END IF

!  IF(openmpi) THEN
!    CALL Action_MPI(Psi,OpPsi,BasisnD,para_Op)
!  ELSE
!    !---------------------------------------------------------------------------
    ! for direct_KEO case
    IF (SG4_omp == 0 .AND. para_Op%direct_KEO) THEN

      DO iG=1,BasisnD%para_SGType2%nb_SG

        !transfert part of the psi(:)%RvecB(:) to PsiR(itab)%V
        allocate(PsiR(size(Psi)))
        ! transfer form basis to grid-------------------------------------------
        DO itab=1,size(Psi)
          CALL tabPackedBasis_TO_tabR_AT_iG(PsiR(itab)%V,psi(itab)%RvecB,    &
                                            iG,BasisnD%para_SGType2)
          !write(out_unitp,*) 'iG,itab,PsiR',iG,itab,PsiR(itab)%V  ;    CALL flush_perso(out_unitp)
        END DO
        !write(out_unitp,*) 'iG',iG, ' unpack done' ;    CALL flush_perso(out_unitp)

        CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4_ompGrid(PsiR,iG,para_Op)
        !write(out_unitp,*) 'iG',iG, ' TabOpPsi done' ;    CALL flush_perso(out_unitp)

        !transfert back, PsiR(itab)%V (HPsi) to the psi(:)%RvecB(:)
        DO itab=1,size(Psi)
          CALL tabR_AT_iG_TO_tabPackedBasis(OpPsi(itab)%RvecB,PsiR(itab)%V,  &
                         iG,BasisnD%para_SGType2,BasisnD%WeightSG(iG))
        END DO
        !write(out_unitp,*) 'iG',iG, ' pack done' ;    CALL flush_perso(out_unitp)

        !deallocate PsiR(:)
        DO itab=1,size(Psi)
          CALL dealloc_TypeRVec(PsiR(itab))
        END DO
        deallocate(PsiR)

        IF (print_level > 0  .AND. BasisnD%para_SGType2%nb_SG > 10**4 .AND.    &
            mod(iG,max(1,BasisnD%para_SGType2%nb_SG/10)) == 0 .AND. MPI_id==0) THEN
          write(out_unitp,'(a)',ADVANCE='no') '---'
          CALL flush_perso(out_unitp)
        END IF

      END DO
    ELSE IF (allocated(BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval)) THEN

      packet_size=max(1,BasisnD%para_SGType2%nb_SG/SG4_maxth/10)
      !$OMP   PARALLEL DEFAULT(NONE)                              &
      !$OMP   SHARED(Psi,OpPsi)                                   &
      !$OMP   SHARED(para_Op,BasisnD)                             &
      !$OMP   SHARED(print_level,out_unitp,SG4_maxth,packet_size) &
      !$OMP   SHARED(MPI_id)                                      &
      !$OMP   PRIVATE(iG,itab,PsiR)                               &
      !$OMP   NUM_THREADS(SG4_maxth)
      allocate(PsiR(size(Psi)))
      !!$OMP   DO SCHEDULE(DYNAMIC,packet_size)
      !!$OMP   DO SCHEDULE(GUIDED)
      !$OMP   DO SCHEDULE(STATIC)
      DO iG=1,BasisnD%para_SGType2%nb_SG
        !transfert part of the psi(:)%RvecB(:) to PsiR(itab)%V
        DO itab=1,size(Psi)
          CALL tabPackedBasis_TO_tabR_AT_iG(PsiR(itab)%V,psi(itab)%RvecB,  &
                                            iG,BasisnD%para_SGType2)
        END DO

        CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4(PsiR,iG,                    &
             BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval(:,iG),para_Op)

        !transfert back, PsiR(itab)%V (HPsi) to the psi(:)%RvecB(:)
        DO itab=1,size(Psi)
          CALL tabR_AT_iG_TO_tabPackedBasis(OpPsi(itab)%RvecB,PsiR(itab)%V,&
                             iG,BasisnD%para_SGType2,BasisnD%WeightSG(iG))
        END DO

        !deallocate PsiR(:)
        DO itab=1,size(Psi)
          CALL dealloc_TypeRVec(PsiR(itab))
        END DO

        IF (print_level > 0  .AND. BasisnD%para_SGType2%nb_SG > 10**4 .AND. &
            mod(iG,max(1,BasisnD%para_SGType2%nb_SG/10)) == 0 .AND. MPI_id==0) THEN
          write(out_unitp,'(a)',ADVANCE='no') '---'
          CALL flush_perso(out_unitp)
        END IF

        !write(out_unitp,*) 'iG done:',iG ; flush(out_unitp)
      END DO
      !$OMP   END DO
      deallocate(PsiR) ! with ifx it has to be removed !!
      !$OMP   END PARALLEL

    ELSE IF (BasisnD%para_SGType2%nb_tasks /= BasisnD%para_SGType2%nb_threads) THEN ! version 2

      !$OMP parallel do                                             &
      !$OMP default(none)                                           &
      !$OMP shared(Psi,OpPsi,para_Op,BasisnD)                       &
      !$OMP private(i)                                              &
      !$OMP num_threads(SG4_maxth)
      DO i=1,BasisnD%para_SGType2%nb_tasks
        CALL sub_TabOpPsi_OF_SeveralDP_FOR_SGtype4(Psi,OpPsi,para_Op,          &
                      BasisnD%para_SGType2%nDval_init(:,i),                    &
               BasisnD%para_SGType2%iG_th(i),BasisnD%para_SGType2%fG_th(i))
      END DO
      !$OMP   END PARALLEL DO
    ELSE
      !to be sure to have the correct number of threads, we use
      !   BasisnD%para_SGType2%nb_threads

      !$OMP parallel                                                &
      !$OMP default(none)                                           &
      !$OMP shared(Psi,OpPsi)                                       &
      !$OMP shared(para_Op,BasisnD,print_level,out_unitp,MPI_id)    &
      !$OMP private(itab,iG,iiG,tab_l,PsiR,ith)                     &
      !$OMP num_threads(BasisnD%para_SGType2%nb_threads)

      !--------------------------------------------------------------
      !-- For the initialization of tab_l(:) and the use of ADD_ONE_TO_nDindex in the parallel loop
      CALL alloc_NParray(tab_l,(/ BasisnD%para_SGType2%nDind_SmolyakRep%ndim /),'tabl_l',name_sub)
      ith = 0
      !$ ith = OMP_GET_THREAD_NUM()
      tab_l(:) = BasisnD%para_SGType2%nDval_init(:,ith+1)
      !--------------------------------------------------------------

      !write(out_unitp,*) 'ith,tab_l(:)',ith,':',tab_l

      ! we are not using the parallel do, to be able to use the correct initialized tab_l with nDval_init
      DO iG=BasisnD%para_SGType2%iG_th(ith+1),BasisnD%para_SGType2%fG_th(ith+1)
        !write(out_unitp,*) 'iG',iG ; flush(out_unitp)

        CALL ADD_ONE_TO_nDindex(BasisnD%para_SGType2%nDind_SmolyakRep,tab_l,iG=iG)

        !transfert part of the psi(:)%RvecB(:) to PsiR(itab)%V
        allocate(PsiR(size(Psi)))
        DO itab=1,size(Psi)
          CALL tabPackedBasis_TO_tabR_AT_iG(PsiR(itab)%V,psi(itab)%RvecB,    &
                                            iG,BasisnD%para_SGType2)
          !write(out_unitp,*) 'iG,itab,PsiR',iG,itab,PsiR(itab)%V
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

        IF (print_level > 0  .AND. BasisnD%para_SGType2%nb_SG > 10**4 .AND. &
            mod(iG,max(1,BasisnD%para_SGType2%nb_SG/10)) == 0 .AND. MPI_id==0) THEN
          write(out_unitp,'(a)',ADVANCE='no') '---'
          CALL flush_perso(out_unitp)
        END IF

        !write(out_unitp,*) 'iG done:',iG ; flush(out_unitp)
      END DO
      CALL dealloc_NParray(tab_l,'tabl_l',name_sub)
      !$OMP   END PARALLEL
    END IF
!  ENDIF ! for openmpi

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
    para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_done=.True.
  END IF

  IF(MPI_id==0) THEN
    DO i=1,size(OpPsi)
      OpPsi_symab = Calc_symab1_EOR_symab2(para_Op%symab,Psi(i)%symab)
      CALL Set_symab_OF_psiBasisRep(OpPsi(i),OpPsi_symab)

      !write(out_unitp,*) 'para_Op,psi symab ',i,para_Op%symab,Psi(i)%symab
      !write(out_unitp,*) 'OpPsi_symab',i,OpPsi(i)%symab
    END DO
  ENDIF

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
!=======================================================================================

 SUBROUTINE sub_TabOpPsi_OF_SeveralDP_FOR_SGtype4(Psi,OpPsi,para_Op,    &
                                                  tab_l_init,initG,endG)
 USE mod_system
 USE mod_nDindex
 USE mod_Coord_KEO,                ONLY : CoordType
 USE mod_basis_set_alloc,          ONLY : basis
 USE mod_basis_BtoG_GtoB_SGType4,  ONLY : tabPackedBasis_TO_tabR_AT_iG, &
                                          tabR_AT_iG_TO_tabPackedBasis, &
                                          TypeRVec,dealloc_TypeRVec

 USE mod_psi,                      ONLY : param_psi,ecri_psi
 USE mod_SetOp,                    ONLY : param_Op,write_param_Op
 IMPLICIT NONE

 TYPE (param_psi), intent(in)      :: Psi(:)
 TYPE (param_psi), intent(inout)   :: OpPsi(:)

 TYPE (param_Op),  intent(inout)   :: para_Op
 integer,          intent(in)      :: initG,endG,tab_l_init(:)

 ! local variables
 TYPE (basis),      pointer        :: BasisnD
 TYPE (TypeRVec),   allocatable    :: PsiR(:)
 integer                           :: iG,itab
 integer,           allocatable    :: tab_l(:)

 !----- for debuging ----------------------------------------------
 character (len=*), parameter :: name_sub='sub_TabOpPsi_OF_SeveralDP_FOR_SGtype4'
 logical, parameter :: debug = .FALSE.
 !logical, parameter :: debug = .TRUE.
 !-----------------------------------------------------------------
 IF (debug) THEN
   write(out_unitp,*) 'BEGINNING ',name_sub
   CALL flush_perso(out_unitp)
 END IF
 !-----------------------------------------------------------------
  BasisnD => para_Op%BasisnD



   !--------------------------------------------------------------
   !-- For the initialization of tab_l(:) and the use of ADD_ONE_TO_nDindex in the parallel loop
   CALL alloc_NParray(tab_l,(/ BasisnD%para_SGType2%nDind_SmolyakRep%ndim /),&
                     'tab_l',name_sub)
   tab_l(:) = tab_l_init(:)
   !--------------------------------------------------------------

   DO iG=initG,endG

     CALL ADD_ONE_TO_nDindex(BasisnD%para_SGType2%nDind_SmolyakRep,tab_l,iG=iG)

     !write(out_unitp,*) 'iG',iG
     !transfert part of the psi(:)%RvecB(:) to PsiR(itab)%V
     allocate(PsiR(size(Psi)))
     DO itab=1,size(Psi)
       CALL tabPackedBasis_TO_tabR_AT_iG(PsiR(itab)%V,psi(itab)%RvecB,    &
                                         iG,BasisnD%para_SGType2)
       !write(out_unitp,*) 'iG,itab,PsiR',iG,itab,PsiR(itab)%V
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

     IF (print_level > 0  .AND. BasisnD%para_SGType2%nb_SG > 10**4 .AND. &
         mod(iG,max(1,BasisnD%para_SGType2%nb_SG/10)) == 0 .AND. MPI_id==0) THEN
       write(out_unitp,'(a)',ADVANCE='no') '---'
       CALL flush_perso(out_unitp)
     END IF

   END DO

!-----------------------------------------------------------
 IF (debug) THEN
   write(out_unitp,*)
   write(out_unitp,*) 'END ',name_sub
 END IF
!-----------------------------------------------------------

  END SUBROUTINE sub_TabOpPsi_OF_SeveralDP_FOR_SGtype4


  SUBROUTINE sub_TabOpPsi_OF_ONEDP_FOR_SGtype4_old(PsiR,iG,tab_l,para_Op,PsiROnGrid)
  USE mod_system
  USE mod_nDindex

  USE mod_Coord_KEO,               ONLY : CoordType

  USE mod_basis_set_alloc,         ONLY : basis
!  USE mod_basis,                   ONLY : Rec_Qact_SG4_with_Tab_iq
  USE mod_basis_BtoG_GtoB_SGType4, ONLY : TypeRVec,                     &
                     BDP_TO_GDP_OF_SmolyakRep,GDP_TO_BDP_OF_SmolyakRep, &
                  DerivOp_TO_RDP_OF_SmolaykRep,tabR2grid_TO_tabR1_AT_iG,&
                     getbis_tab_nq,getbis_tab_nb
  USE mod_SetOp,                      ONLY : param_Op,write_param_Op
  IMPLICIT NONE

  TYPE (TypeRVec),                    intent(inout)    :: PsiR(:)
  logical,          optional,         intent(in)       :: PsiROnGrid
  integer,                            intent(in)       :: iG,tab_l(:)

  TYPE (param_Op),                    intent(inout)    :: para_Op

  !local variables
  TYPE (CoordType), pointer :: mole
  TYPE(basis),    pointer :: BasisnD

  integer :: iq,nb,nq,i,j,itab,D
  logical :: PsiROnBasis_loc

  integer                               :: derive_termQdyn(2)

  real (kind=Rkind),  allocatable       :: V(:,:,:)
  real(kind=Rkind),   allocatable       :: GridOp(:,:,:,:)
  real (kind=Rkind),  allocatable       :: GGiq(:,:,:)
  real (kind=Rkind),  allocatable       :: sqRhoOVERJac(:),Jac(:)
  real (kind=Rkind)                     :: Rho

  real (kind=Rkind),  allocatable       :: PsiRj(:,:),PsiRi(:),OpPsiR(:)
  real (kind=Rkind),  allocatable       :: PsiTemp1(:),PsiTemp2(:)

  integer,            allocatable       :: tab_nq(:)
  integer,            allocatable       :: tab_nb(:)
  integer :: iterm,iterm00
  integer :: err_sub


  !----- for debuging ----------------------------------------------
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub='sub_TabOpPsi_OF_ONEDP_FOR_SGtype4_old'
  logical, parameter :: debug = .FALSE.
  !logical, parameter :: debug = .TRUE.
  !-----------------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*) 'nb_bie,nb_baie',para_Op%nb_bie,para_Op%nb_baie
    write(out_unitp,*) 'nb_act1',para_Op%mole%nb_act1
    write(out_unitp,*) 'nb_var',para_Op%mole%nb_var
    write(out_unitp,*) 'nb psi',size(PsiR)
    write(out_unitp,*) 'iG, tab_l(:)',iG,':',tab_l(:)
    write(out_unitp,*) 'nq',para_Op%BasisnD%para_SGType2%tab_nq_OF_SRep(iG)

    CALL flush_perso(out_unitp)
  END IF
  !-----------------------------------------------------------------
  !CALL Check_mem()

   !write(out_unitp,*) '================================' ; flush(out_unitp)
   !write(out_unitp,*) '============ START =============' ; flush(out_unitp)

   mole    => para_Op%mole
   BasisnD => para_Op%BasisnD

   IF (present(PsiROnGrid)) THEN
     PsiROnBasis_loc = .NOT. PsiROnGrid
   ELSE
     PsiROnBasis_loc = .TRUE.
   END IF
   D = size(tab_l)

   CALL alloc_NParray(tab_nb,(/ D /),'tab_nb',name_sub)
   CALL alloc_NParray(tab_nq,(/ D /),'tab_nq',name_sub)

   tab_nq(:) = getbis_tab_nq(tab_l,BasisnD%tab_basisPrimSG)
   tab_nb(:) = getbis_tab_nb(tab_l,BasisnD%tab_basisPrimSG)

   nb = size(PsiR(1)%V) ! because PsiR is in basis Rep
   nq = BasisnD%para_SGType2%tab_nq_OF_SRep(iG)

   SELECT CASE (para_Op%type_Op)
   CASE (0) ! 0 : Scalar

     CALL get_OpGrid_type0_OF_ONEDP_FOR_SG4(iG,tab_l,para_Op,V)

     DO itab=1,size(PsiR)

       IF (PsiROnBasis_loc) THEN
         ! partial B to G
         CALL BDP_TO_GDP_OF_SmolyakRep(PsiR(itab)%V,BasisnD%tab_basisPrimSG,&
                                       tab_l,tab_nq,tab_nb)
         PsiR(itab)%V(:) = PsiR(itab)%V(:) * V(:,1,1)

         ! partial B to G
         CALL GDP_TO_BDP_OF_SmolyakRep(PsiR(itab)%V,BasisnD%tab_basisPrimSG,&
                                       tab_l,tab_nq,tab_nb)
         ! now PsiR(itab)%V is on the Basis
       ELSE
         PsiR(itab)%V(:) = PsiR(itab)%V(:) * V(:,1,1)
       END IF

     END DO

     IF (allocated(V)) CALL dealloc_NParray(V,'V',name_sub)

   CASE (1) ! 1 : H: F2.d^2 + F1.d^1 + V
     IF (debug) write(out_unitp,*) 'nb_Term',para_Op%nb_Term

     CALL get_OpGrid_type1_OF_ONEDP_FOR_SG4_new2(iG,tab_l,para_Op,GridOp)

     CALL alloc_NParray(PsiTemp1,       (/ nq /),      'PsiTemp1',       name_sub)
     CALL alloc_NParray(PsiTemp2,       (/ nq /),      'PsiTemp2',       name_sub)

     DO itab=1,size(PsiR)

       IF (PsiROnBasis_loc) THEN
         ! partial B to G
         CALL BDP_TO_GDP_OF_SmolyakRep(PsiR(itab)%V,BasisnD%tab_basisPrimSG,&
                                       tab_l,tab_nq,tab_nb)

         PsiTemp1(:) = ZERO
         DO iterm=1,para_Op%nb_Term

           IF (para_Op%OpGrid(iterm)%grid_zero) CYCLE

           PsiTemp2(:) = PsiR(itab)%V(:)

           CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiTemp2,BasisnD%tab_basisPrimSG, &
                               tab_l,tab_nq,para_Op%derive_termQdyn(:,iterm))

           PsiTemp1(:) = PsiTemp1(:) + PsiTemp2 * GridOp(:,1,1,iterm)

         END DO

         PsiR(itab)%V(:) = PsiTemp1(:)

         ! partial B to G
         CALL GDP_TO_BDP_OF_SmolyakRep(PsiR(itab)%V,BasisnD%tab_basisPrimSG,&
                                       tab_l,tab_nq,tab_nb)
         ! now PsiR(itab)%V is on the Basis
       ELSE

         PsiTemp1(:) = ZERO
         DO iterm=1,para_Op%nb_Term

           IF (para_Op%OpGrid(iterm)%grid_zero) CYCLE

           PsiTemp2(:) = PsiR(itab)%V(:)

           CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiTemp2,BasisnD%tab_basisPrimSG, &
                               tab_l,tab_nq,para_Op%derive_termQdyn(:,iterm))

           PsiTemp1(:) = PsiTemp1(:) + PsiTemp2 * GridOp(:,1,1,iterm)

         END DO

         PsiR(itab)%V(:) = PsiTemp1(:)

       END IF

     END DO

     IF (allocated(GridOp))   CALL dealloc_NParray(GridOp,  'GridOp',  name_sub)
     IF (allocated(PsiTemp1)) CALL dealloc_NParray(PsiTemp1,'PsiTemp1',name_sub)
     IF (allocated(PsiTemp2)) CALL dealloc_NParray(PsiTemp2,'PsiTemp2',name_sub)

   CASE (10) !10: H: d^1 G d^1 +V

     CALL alloc_NParray(PsiRj,(/ nq,mole%nb_act1 /),'PsiRj',       name_sub)
     CALL alloc_NParray(PsiRi,       (/ nq /),      'PsiRi',       name_sub)


     CALL get_OpGrid_type10_OF_ONEDP_FOR_SG4(iG,tab_l,para_Op,V,GGiq,sqRhoOVERJac,Jac)


     DO itab=1,size(PsiR)

       CALL alloc_NParray(OpPsiR,(/nq/),'OpPsiR',name_sub)

       IF (PsiROnBasis_loc) THEN
         ! partial B to G
         CALL BDP_TO_GDP_OF_SmolyakRep(PsiR(itab)%V,BasisnD%tab_basisPrimSG,&
                                       tab_l,tab_nq,tab_nb)
       END IF

       ! now PsiR is on the grid
       !write(out_unitp,*) ' R Grid of PsiR(itab) :',PsiR(itab)%V
       !write(out_unitp,*) ' R Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(out_unitp)

       ! multiplication by sqRhoOVERJac
       PsiR(itab)%V(:) = PsiR(itab)%V(:) * sqRhoOVERJac(:)

       !write(out_unitp,*) ' R*sq Grid of PsiR(itab) :',PsiR
       !write(out_unitp,*) ' R*sq ... Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(out_unitp)


       ! derivative with respect to Qj
       DO j=1,mole%nb_act1
         derive_termQdyn(:) = (/ mole%liste_QactTOQdyn(j),0 /)

         PsiRj(:,j) = PsiR(itab)%V(:)
         CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiRj(:,j),BasisnD%tab_basisPrimSG, &
                                           tab_l,tab_nq,derive_termQdyn)

         !write(out_unitp,*) 'dQj Grid of PsiR(itab)',j,PsiRj(:,j) ; flush(out_unitp)
       END DO
       !write(out_unitp,*) ' dQj R*sq ... Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(out_unitp)

       OpPsiR(:) = ZERO
       DO i=1,mole%nb_act1

         PsiRi(:) = ZERO
         DO j=1,mole%nb_act1
           PsiRi(:) = PsiRi(:) + GGiq(:,j,i) * PsiRj(:,j)
         END DO
         PsiRi(:) = PsiRi(:) * Jac(:)
         !write(out_unitp,*) ' Jac Gij* ... Grid of SRep : done',iG,itab,i,OMP_GET_THREAD_NUM() ; flush(out_unitp)

         derive_termQdyn(:) = (/ mole%liste_QactTOQdyn(i),0 /)

         CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiRi(:),BasisnD%tab_basisPrimSG,&
                                           tab_l,tab_nq,derive_termQdyn)

         OpPsiR(:) = OpPsiR(:) + PsiRi(:)

       END DO
       !write(out_unitp,*) ' dQi Jac Gij* ... Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(out_unitp)

       OpPsiR(:) = (  -HALF*OpPsiR(:)/Jac(:) + PsiR(itab)%V(:)*V(:,1,1) ) / sqRhoOVERJac(:)
       !write(out_unitp,*) ' OpPsiR Grid of PsiR(itab)',OpPsiR
       !write(out_unitp,*) ' OpPsiR Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(out_unitp)

       !tranfert OpPsiR (on the grid) to PsiR(itab)
       PsiR(itab)%V(:) = OpPsiR(:)

       CALL dealloc_NParray(OpPsiR,'OpPsiR', name_sub)

       IF (PsiROnBasis_loc) THEN
         ! partial B to G
         CALL GDP_TO_BDP_OF_SmolyakRep(PsiR(itab)%V,BasisnD%tab_basisPrimSG,&
                                     tab_l,tab_nq,tab_nb)
         ! now PsiR(itab)%V is on the Basis
       END IF


     END DO

     IF (allocated(PsiRi))        CALL dealloc_NParray(PsiRi,       'PsiRi',       name_sub)
     IF (allocated(PsiRj))        CALL dealloc_NParray(PsiRj,       'PsiRj',       name_sub)
     IF (allocated(GGiq))         CALL dealloc_NParray(GGiq,        'GGiq',        name_sub)
     IF (allocated(sqRhoOVERJac)) CALL dealloc_NParray(sqRhoOVERJac,'sqRhoOVERJac',name_sub)
     IF (allocated(Jac))          CALL dealloc_NParray(Jac,         'Jac',         name_sub)
     IF (allocated(V))            CALL dealloc_NParray(V,           'V',           name_sub)


   CASE Default
     STOP 'no default'
   END SELECT

   IF (allocated(tab_nb)) CALL dealloc_NParray(tab_nb,'tab_nb',name_sub)
   IF (allocated(tab_nq)) CALL dealloc_NParray(tab_nq,'tab_nq',name_sub)
   !write(out_unitp,*) '============ END ===============' ; flush(out_unitp)
   !write(out_unitp,*) '================================' ; flush(out_unitp)
  !-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'END ',name_sub
    CALL flush_perso(out_unitp)
  END IF
  !-----------------------------------------------------------
  !CALL UnCheck_mem() ; stop
  END SUBROUTINE sub_TabOpPsi_OF_ONEDP_FOR_SGtype4_old

!=======================================================================================
  SUBROUTINE sub_TabOpPsi_OF_ONEDP_FOR_SGtype4(PsiR,iG,tab_l,para_Op,PsiROnGrid)
  USE mod_system
  USE mod_nDindex

  USE mod_Coord_KEO,               ONLY : CoordType

  USE mod_basis_set_alloc,         ONLY : basis
  !USE mod_basis,                   ONLY : Rec_Qact_SG4_with_Tab_iq
  USE mod_basis_BtoG_GtoB_SGType4, ONLY : TypeRVec,                     &
                     BDP_TO_GDP_OF_SmolyakRep,GDP_TO_BDP_OF_SmolyakRep, &
                  DerivOp_TO_RDP_OF_SmolaykRep,tabR2grid_TO_tabR1_AT_iG,&
                     getbis_tab_nq,getbis_tab_nb
  USE mod_SetOp,                      ONLY : param_Op,write_param_Op
  IMPLICIT NONE

  TYPE (TypeRVec),                    intent(inout)    :: PsiR(:)
  logical,          optional,         intent(in)       :: PsiROnGrid
  integer,                            intent(in)       :: iG,tab_l(:)

  TYPE (param_Op),                    intent(inout)    :: para_Op

  !local variables
  TYPE (CoordType), pointer :: mole
  TYPE(basis),    pointer :: BasisnD

  integer :: iq,nb,nq,i,j,itab,D
  integer :: ib0,jb0,nb0,iqi,iqf,jqi,jqf,ibi,ibf

  logical :: PsiROnBasis_loc

  integer                               :: derive_termQdyn(2)

  real (kind=Rkind),  allocatable       :: V(:,:,:)
  real (kind=Rkind),  allocatable       :: GridOp(:,:,:,:)
  real (kind=Rkind),  allocatable       :: GGiq(:,:,:)
  real (kind=Rkind),  allocatable       :: sqRhoOVERJac(:),Jac(:)
  real (kind=Rkind)                     :: Rho

  real (kind=Rkind),  allocatable       :: PsiRj(:,:),PsiRi(:),OpPsiR(:),VPsi(:,:)
  real (kind=Rkind),  allocatable       :: OpPsi(:,:) ! size(nq,nb0)
  real (kind=Rkind),  allocatable       :: Psi_ch(:,:) ! size(nq,nb0)

  integer,            allocatable       :: tab_nq(:)
  integer,            allocatable       :: tab_nb(:)
  integer :: iterm,iterm00
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
    write(out_unitp,*) 'nb psi',size(PsiR)
    write(out_unitp,*) 'iG, tab_l(:)',iG,':',tab_l(:)
    write(out_unitp,*) 'nb0',para_Op%BasisnD%para_SGType2%nb0
    write(out_unitp,*) 'nq',para_Op%BasisnD%para_SGType2%tab_nq_OF_SRep(iG)
    write(out_unitp,*) 'nb',para_Op%BasisnD%para_SGType2%tab_nb_OF_SRep(iG)
    IF (present(PsiROnGrid)) write(out_unitp,*) 'PsiROnGrid',PsiROnGrid
    CALL flush_perso(out_unitp)
  END IF
  !-----------------------------------------------------------------
  !CALL Check_mem()

   !write(out_unitp,*) '================================' ; flush(out_unitp)
   !write(out_unitp,*) '============ START =============' ; flush(out_unitp)

   mole    => para_Op%mole
   BasisnD => para_Op%BasisnD

   IF (present(PsiROnGrid)) THEN
     PsiROnBasis_loc = .NOT. PsiROnGrid
   ELSE
     PsiROnBasis_loc = .TRUE.
   END IF
   D = size(tab_l)

   CALL alloc_NParray(tab_nb,(/ D /),'tab_nb',name_sub)
   CALL alloc_NParray(tab_nq,(/ D /),'tab_nq',name_sub)

   tab_nq(:) = getbis_tab_nq(tab_l,BasisnD%tab_basisPrimSG)
   tab_nb(:) = getbis_tab_nb(tab_l,BasisnD%tab_basisPrimSG)

   nb  = BasisnD%para_SGType2%tab_nb_OF_SRep(iG)
   nq  = BasisnD%para_SGType2%tab_nq_OF_SRep(iG)
   nb0 = BasisnD%para_SGType2%nb0

   SELECT CASE (para_Op%type_Op)
   CASE (0) ! 0 : Scalar

     CALL get_OpGrid_type0_OF_ONEDP_FOR_SG4(iG,tab_l,para_Op,V)

     CALL alloc_NParray(OpPsi,(/nq,nb0/),'OpPsi',name_sub)

     DO itab=1,size(PsiR)

       OpPsi(:,:) = ZERO

       IF (PsiROnBasis_loc) THEN
         ! partial B to G
         CALL BDP_TO_GDP_OF_SmolyakRep(PsiR(itab)%V,BasisnD%tab_basisPrimSG,&
                                       tab_l,tab_nq,tab_nb,nb0)
       END IF

       ! VPsi calculation, it has to be done before because V is not diagonal
       DO ib0=1,nb0
         jqi = 1
         jqf = nq
         DO jb0=1,nb0
           OpPsi(:,ib0) = OpPsi(:,ib0) + V(:,ib0,jb0) * PsiR(itab)%V(jqi:jqf)
           jqi = jqf+1
           jqf = jqf+nq
         END DO
       END DO


       PsiR(itab)%V(:) = reshape(OpPsi,shape=(/nq*nb0/))

       IF (PsiROnBasis_loc) THEN
         ! partial B to G
         CALL GDP_TO_BDP_OF_SmolyakRep(PsiR(itab)%V,BasisnD%tab_basisPrimSG,&
                                       tab_l,tab_nq,tab_nb,nb0)
         ! now PsiR(itab)%V is on the Basis
       END IF

     END DO

     IF (allocated(V))     CALL dealloc_NParray(V,'V',name_sub)
     IF (allocated(OpPsi)) CALL dealloc_NParray(OpPsi,'OpPsi',name_sub)

   CASE (1) ! 1 : H: F2.d^2 + F1.d^1 + V
     IF (debug) write(out_unitp,*) 'nb_Term',para_Op%nb_Term

     CALL get_OpGrid_type1_OF_ONEDP_FOR_SG4_new2(iG,tab_l,para_Op,GridOp)

     CALL alloc_NParray(OpPsi,(/nq,nb0/),'OpPsi',name_sub)
     CALL alloc_NParray(Psi_ch,(/nq,nb0/),'Psi_ch',name_sub)


     DO itab=1,size(PsiR)

       IF (PsiROnBasis_loc) THEN
         ! partial B to G
           CALL BDP_TO_GDP_OF_SmolyakRep(PsiR(itab)%V,                  &
                                         BasisnD%tab_basisPrimSG,       &
                                         tab_l,tab_nq,tab_nb,nb0)
       END IF

       OpPsi(:,:) = ZERO

       DO iterm=1,para_Op%nb_Term
         IF (para_Op%OpGrid(iterm)%grid_zero) CYCLE

         Psi_ch(:,:) = reshape(PsiR(itab)%V,shape=(/nq,nb0/))
         DO ib0=1,nb0
           CALL DerivOp_TO_RDP_OF_SmolaykRep(Psi_ch(:,ib0),BasisnD%tab_basisPrimSG, &
                               tab_l,tab_nq,para_Op%derive_termQdyn(:,iterm))
         END DO


         ! GridOp(:,1,1,iterm)Psi_ch calculation
         DO ib0=1,nb0
         DO jb0=1,nb0
           OpPsi(:,ib0) = OpPsi(:,ib0) + GridOp(:,ib0,jb0,iterm) * Psi_ch(:,jb0)
         END DO
         END DO

       END DO


       PsiR(itab)%V(:) = reshape(OpPsi,shape=(/nq*nb0/))



       IF (PsiROnBasis_loc) THEN
         ! partial B to G
           CALL GDP_TO_BDP_OF_SmolyakRep(PsiR(itab)%V,                  &
                                         BasisnD%tab_basisPrimSG,       &
                                         tab_l,tab_nq,tab_nb,nb0)
         ! now PsiR(itab)%V is on the basis
       END IF

     END DO

     IF (allocated(GridOp))   CALL dealloc_NParray(GridOp,  'GridOp',  name_sub)
     IF (allocated(OpPsi))    CALL dealloc_NParray(OpPsi,   'OpPsi',name_sub)
     IF (allocated(Psi_ch))   CALL dealloc_NParray(Psi_ch,  'Psi_ch',name_sub)

   CASE (10) !10: H: d^1 G d^1 +V
     CALL alloc_NParray(PsiRj, (/ nq,mole%nb_act1 /),'PsiRj',       name_sub)
     CALL alloc_NParray(PsiRi,       (/ nq /),       'PsiRi',       name_sub)
     CALL alloc_NParray(VPsi,  (/nq,nb0/),           'VPsi',        name_sub)
     CALL alloc_NParray(OpPsiR,(/nq/),               'OpPsiR',      name_sub)


     CALL get_OpGrid_type10_OF_ONEDP_FOR_SG4(iG,tab_l,para_Op,V,GGiq,sqRhoOVERJac,Jac)


     DO itab=1,size(PsiR)

       IF (PsiROnBasis_loc) THEN
         ! partial B to G
         CALL BDP_TO_GDP_OF_SmolyakRep(PsiR(itab)%V,BasisnD%tab_basisPrimSG,&
                                       tab_l,tab_nq,tab_nb,nb0)
       END IF

       ! VPsi calculation, it has to be done before because V is not diagonal
       VPsi(:,:) = ZERO
       DO ib0=1,nb0
         jqi = 1
         jqf = nq
         DO jb0=1,nb0
           VPsi(:,ib0) = VPsi(:,ib0) + V(:,ib0,jb0) * PsiR(itab)%V(jqi:jqf)
           jqi = jqf+1
           jqf = jqf+nq
         END DO
       END DO

       iqi = 1
       iqf = nq
       DO ib0=1,nb0

         ! multiplication by sqRhoOVERJac
         PsiR(itab)%V(iqi:iqf) = PsiR(itab)%V(iqi:iqf) * sqRhoOVERJac(:)

         !write(out_unitp,*) ' R*sq Grid of PsiR :',PsiR%V
         !write(out_unitp,*) ' R*sq ... Grid of PsiR : done',OMP_GET_THREAD_NUM() ; flush(out_unitp)


         ! derivative with respect to Qj
         DO j=1,mole%nb_act1
           derive_termQdyn(:) = (/ mole%liste_QactTOQdyn(j),0 /)

           PsiRj(:,j) = PsiR(itab)%V(iqi:iqf)
           CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiRj(:,j),BasisnD%tab_basisPrimSG, &
                                             tab_l,tab_nq,derive_termQdyn)
           !write(out_unitp,*) 'dQj Grid of PsiR',j,PsiRj(:,j) ; flush(out_unitp)
         END DO
         !write(out_unitp,*) ' dQj R*sq ... Grid of PsiR : done',OMP_GET_THREAD_NUM() ; flush(out_unitp)

         OpPsiR(:) = ZERO
         DO i=1,mole%nb_act1

           PsiRi(:) = ZERO
           DO j=1,mole%nb_act1
             PsiRi(:) = PsiRi(:) + GGiq(:,j,i) * PsiRj(:,j)
           END DO
           PsiRi(:) = PsiRi(:) * Jac(:)
           !write(out_unitp,*) ' Jac Gij* ... Grid of SRep : done',iG,i,OMP_GET_THREAD_NUM() ; flush(out_unitp)

           derive_termQdyn(:) = (/ mole%liste_QactTOQdyn(i),0 /)

           CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiRi(:),BasisnD%tab_basisPrimSG,&
                                               tab_l,tab_nq,derive_termQdyn)

           !write(out_unitp,*) 'shape OpPsiR,PsiRi ',shape(OpPsiR),shape(PsiRi)
           !write(out_unitp,*) 'DerivOp_TO_RDP_OF_SmolaykRep : done',OMP_GET_THREAD_NUM() ; flush(out_unitp)

           OpPsiR(:) = OpPsiR(:) + PsiRi(:)
         END DO

         !write(out_unitp,*) ' dQi Jac Gij* ... Grid of PsiR : done',OMP_GET_THREAD_NUM() ; flush(out_unitp)

         OpPsiR(:) =   -HALF*OpPsiR(:) / ( Jac(:)*sqRhoOVERJac(:) ) + VPsi(:,ib0)
         !write(out_unitp,*) ' OpPsiR Grid of PsiR',OpPsiR
         !write(out_unitp,*) ' OpPsiR Grid of PsiR : done',OMP_GET_THREAD_NUM() ; flush(out_unitp)

         !tranfert OpPsiR (on the grid) to PsiR
         PsiR(itab)%V(iqi:iqf) = OpPsiR(:)

         iqi = iqf+1
         iqf = iqf+nq
       END DO

       IF (PsiROnBasis_loc) THEN
         ! partial B to G
         CALL GDP_TO_BDP_OF_SmolyakRep(PsiR(itab)%V,BasisnD%tab_basisPrimSG,&
                                     tab_l,tab_nq,tab_nb,nb0)
         ! now PsiR(itab)%V is on the Basis
       END IF


     END DO

     IF (allocated(PsiRi))        CALL dealloc_NParray(PsiRi,       'PsiRi',       name_sub)
     IF (allocated(PsiRj))        CALL dealloc_NParray(PsiRj,       'PsiRj',       name_sub)
     IF (allocated(GGiq))         CALL dealloc_NParray(GGiq,        'GGiq',        name_sub)
     IF (allocated(sqRhoOVERJac)) CALL dealloc_NParray(sqRhoOVERJac,'sqRhoOVERJac',name_sub)
     IF (allocated(Jac))          CALL dealloc_NParray(Jac,         'Jac',         name_sub)
     IF (allocated(V))            CALL dealloc_NParray(V,           'V',           name_sub)
     IF (allocated(VPsi))         CALL dealloc_NParray(VPsi        ,'VPsi',        name_sub)


   CASE Default
     STOP 'no default'
   END SELECT

   IF (allocated(tab_nb)) CALL dealloc_NParray(tab_nb,'tab_nb',name_sub)
   IF (allocated(tab_nq)) CALL dealloc_NParray(tab_nq,'tab_nq',name_sub)
   !write(out_unitp,*) '============ END ===============' ; flush(out_unitp)
   !write(out_unitp,*) '================================' ; flush(out_unitp)
  !-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'END ',name_sub
    CALL flush_perso(out_unitp)
  END IF
  !-----------------------------------------------------------
  !CALL UnCheck_mem() ; stop
  END SUBROUTINE sub_TabOpPsi_OF_ONEDP_FOR_SGtype4
!=======================================================================================

!=======================================================================================
!> sub_TabOpPsi_OF_ONEDP_FOR_SGtype4 for working on fulll Smolyak rep. on basis
!=======================================================================================
  SUBROUTINE sub_TabOpPsi_OF_ONEDP_FOR_SGtype4_SRB_MPI(iG,psi,OpPsi,tab_l,para_Op)
  USE mod_system
  USE mod_nDindex
  USE mod_Coord_KEO,              ONLY:CoordType
  USE mod_basis_set_alloc,        ONLY:basis
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

END SUBROUTINE sub_TabOpPsi_OF_ONEDP_FOR_SGtype4_SRG_MPI
!=======================================================================================

  SUBROUTINE sub_TabOpPsi_OF_ONEGDP_WithOp_FOR_SGtype4(PsiR,iG,para_Op, &
                                                   V,GG,sqRhoOVERJac,Jac)
  USE mod_system

  USE mod_Coord_KEO,               ONLY : CoordType

  USE mod_basis_set_alloc,         ONLY : basis
  USE mod_basis_BtoG_GtoB_SGType4, ONLY : TypeRVec,                     &
                     BDP_TO_GDP_OF_SmolyakRep,GDP_TO_BDP_OF_SmolyakRep, &
                   DerivOp_TO_RDP_OF_SmolaykRep,tabR2grid_TO_tabR1_AT_iG,&
                     getbis_tab_nq,getbis_tab_nb
  USE mod_SetOp,                      ONLY : param_Op,write_param_Op
  IMPLICIT NONE

  TYPE (TypeRVec),                    intent(inout)    :: PsiR(:)
  integer,                            intent(in)       :: iG

  TYPE (param_Op),                    intent(inout)    :: para_Op
  real (kind=Rkind),                  intent(in)       :: V(:)
  real (kind=Rkind),                  intent(in)       :: GG(:,:,:)
  real (kind=Rkind),                  intent(in)       :: sqRhoOVERJac(:),Jac(:)

  !local variables
  TYPE (CoordType), pointer :: mole
  TYPE(basis),    pointer :: BasisnD

  integer :: iq,nb,nq,i,j,itab,D
  integer :: ib0,jb0,nb0,iqi,iqf,jqi,jqf

  integer                               :: derive_termQdyn(2)


  real (kind=Rkind)                     :: Rho

  real (kind=Rkind),  allocatable       :: PsiRj(:,:),PsiRi(:),OpPsiR(:)
  real (kind=Rkind),  allocatable       :: VPsi(:,:) ! size(nq,nb0)

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

   !write(out_unitp,*) '================================' ; flush(out_unitp)
   !write(out_unitp,*) '============ START =============' ; flush(out_unitp)

   mole    => para_Op%mole
   BasisnD => para_Op%BasisnD

   D = size(BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval,dim=1)

   CALL alloc_NParray(tab_l ,(/ D /),'tab_l', name_sub)
   CALL alloc_NParray(tab_nb,(/ D /),'tab_nb',name_sub)
   CALL alloc_NParray(tab_nq,(/ D /),'tab_nq',name_sub)

   tab_l(:)  = BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval(:,iG)
   tab_nq(:) = getbis_tab_nq(tab_l,BasisnD%tab_basisPrimSG)
   tab_nb(:) = getbis_tab_nb(tab_l,BasisnD%tab_basisPrimSG)


   nb = BasisnD%para_SGType2%tab_nb_OF_SRep(iG)
   nq = BasisnD%para_SGType2%tab_nq_OF_SRep(iG)
   nb0 = BasisnD%para_SGType2%nb0

   CALL alloc_NParray(PsiRj,(/ nq,mole%nb_act1 /),'PsiRj',       name_sub)
   CALL alloc_NParray(PsiRi,       (/ nq /),      'PsiRi',       name_sub)

   DO itab=1,size(PsiR)

     CALL alloc_NParray(OpPsiR,(/nq/),'OpPsiR',name_sub)

     ! partial B to G
     !!!!! no because PsiR(itab)%V(:) is already on the grid


     ! multiplication by sqRhoOVERJac
     PsiR(itab)%V(:) = PsiR(itab)%V(:) * sqRhoOVERJac(:)

     !write(out_unitp,*) ' R*sq Grid of PsiR(itab) :',PsiR
     !write(out_unitp,*) ' R*sq ... Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(out_unitp)


     ! derivative with respect to Qj
     DO j=1,mole%nb_act1
       derive_termQdyn(:) = (/ mole%liste_QactTOQdyn(j),0 /)

       PsiRj(:,j) = PsiR(itab)%V(:)
       CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiRj(:,j),BasisnD%tab_basisPrimSG, &
                                         tab_l,tab_nq,derive_termQdyn)

       !write(out_unitp,*) 'dQj Grid of PsiR(itab)',j,PsiRj(:,j) ; flush(out_unitp)
     END DO
     !write(out_unitp,*) ' dQj R*sq ... Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(out_unitp)

     OpPsiR(:) = ZERO
     DO i=1,mole%nb_act1

       PsiRi(:) = ZERO
       DO j=1,mole%nb_act1
         PsiRi(:) = PsiRi(:) + GG(:,j,i) * PsiRj(:,j)
       END DO
       PsiRi(:) = PsiRi(:) * Jac(:)
       !write(out_unitp,*) ' Jac Gij* ... Grid of SRep : done',iG,itab,i,OMP_GET_THREAD_NUM() ; flush(out_unitp)

       derive_termQdyn(:) = (/ mole%liste_QactTOQdyn(i),0 /)

       CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiRi(:),BasisnD%tab_basisPrimSG,&
                                         tab_l,tab_nq,derive_termQdyn)

       OpPsiR(:) = OpPsiR(:) + PsiRi(:)

     END DO
     !write(out_unitp,*) ' dQi Jac Gij* ... Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(out_unitp)

     OpPsiR(:) = (  -HALF*OpPsiR(:)/Jac(:) + PsiR(itab)%V(:)*V(:) ) / sqRhoOVERJac(:)
     !write(out_unitp,*) ' OpPsiR Grid of PsiR(itab)',OpPsiR
     !write(out_unitp,*) ' OpPsiR Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(out_unitp)

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

   !write(out_unitp,*) '============ END ===============' ; flush(out_unitp)
   !write(out_unitp,*) '================================' ; flush(out_unitp)

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

  USE mod_Coord_KEO,               ONLY : CoordType, get_d0GG

  USE mod_basis_set_alloc,         ONLY : basis
  USE mod_basis,                   ONLY : Rec_Qact_SG4
  USE mod_basis_BtoG_GtoB_SGType4, ONLY : TypeRVec,                     &
                     BDP_TO_GDP_OF_SmolyakRep,GDP_TO_BDP_OF_SmolyakRep, &
                   DerivOp_TO_RDP_OF_SmolaykRep,tabR2grid_TO_tabR1_AT_iG,&
                     getbis_tab_nq,getbis_tab_nb
  USE mod_SetOp,                      ONLY : param_Op,write_param_Op
  IMPLICIT NONE

  TYPE (TypeRVec),                    intent(inout)    :: PsiR(:)
  integer,                            intent(in)       :: iG

  TYPE (param_Op),                    intent(inout)    :: para_Op


  !local variables
  TYPE (CoordType), pointer :: mole
  TYPE(basis),    pointer :: BasisnD

  integer :: iq,nb,nq,i,j,itab,D
  integer :: ib0,jb0,nb0,iqi,iqf,jqi,jqf

  integer                               :: derive_termQdyn(2)

  real (kind=Rkind),  allocatable       :: V(:,:,:)
  real (kind=Rkind),  allocatable       :: Qact(:)
  real (kind=Rkind),  allocatable       :: GGiq(:,:,:)
  TYPE (Type_nDindex)                   :: nDind_DPG    ! multidimensional DP index
  real (kind=Rkind),  allocatable       :: sqRhoOVERJac(:),Jac(:)
  real (kind=Rkind)                     :: Rho

  real (kind=Rkind),  allocatable       :: PsiRj(:,:),PsiRi(:),OpPsiR(:)
  real (kind=Rkind),  allocatable       :: VPsi(:,:) ! size(nq,nb0)

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

   !write(out_unitp,*) '================================' ; flush(out_unitp)
   !write(out_unitp,*) '============ START =============' ; flush(out_unitp)


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


   nb  = BasisnD%para_SGType2%tab_nb_OF_SRep(iG)
   nq  = BasisnD%para_SGType2%tab_nq_OF_SRep(iG)
   nb0 = BasisnD%para_SGType2%nb0

   CALL alloc_NParray(PsiRj,(/ nq,mole%nb_act1 /),'PsiRj',       name_sub)
   CALL alloc_NParray(PsiRi,       (/ nq /),      'PsiRi',       name_sub)
   CALL alloc_NParray(sqRhoOVERJac,(/ nq /),      'sqRhoOVERJac',name_sub)
   CALL alloc_NParray(Jac,         (/ nq /),      'Jac',         name_sub)


   CALL get_OpGrid_type10_OF_ONEDP_FOR_SG4(iG,tab_l,para_Op,V)

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
     !write(out_unitp,*) ' R Grid of PsiR(itab) :',PsiR(itab)%V
     !write(out_unitp,*) ' R Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(out_unitp)

     ! multiplication by sqRhoOVERJac
     PsiR(itab)%V(:) = PsiR(itab)%V(:) * sqRhoOVERJac(:)

     !write(out_unitp,*) ' R*sq Grid of PsiR(itab) :',PsiR
     !write(out_unitp,*) ' R*sq ... Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(out_unitp)


     ! derivative with respect to Qj
     DO j=1,mole%nb_act1
       derive_termQdyn(:) = (/ mole%liste_QactTOQdyn(j),0 /)

       PsiRj(:,j) = PsiR(itab)%V(:)
       CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiRj(:,j),BasisnD%tab_basisPrimSG, &
                                         tab_l,tab_nq,derive_termQdyn)

       !write(out_unitp,*) 'dQj Grid of PsiR(itab)',j,PsiRj(:,j) ; flush(out_unitp)
     END DO
     !write(out_unitp,*) ' dQj R*sq ... Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(out_unitp)

     OpPsiR(:) = ZERO
     DO i=1,mole%nb_act1

       PsiRi(:) = ZERO
       DO j=1,mole%nb_act1
         PsiRi(:) = PsiRi(:) + GGiq(:,j,i) * PsiRj(:,j)
       END DO
       PsiRi(:) = PsiRi(:) * Jac(:)
       !write(out_unitp,*) ' Jac Gij* ... Grid of SRep : done',iG,itab,i,OMP_GET_THREAD_NUM() ; flush(out_unitp)

       derive_termQdyn(:) = (/ mole%liste_QactTOQdyn(i),0 /)

       CALL DerivOp_TO_RDP_OF_SmolaykRep(PsiRi(:),BasisnD%tab_basisPrimSG,&
                                         tab_l,tab_nq,derive_termQdyn)

       OpPsiR(:) = OpPsiR(:) + PsiRi(:)

     END DO
     !write(out_unitp,*) ' dQi Jac Gij* ... Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(out_unitp)

     OpPsiR(:) = (  -HALF*OpPsiR(:)/Jac(:) + PsiR(itab)%V(:)*V(:,1,1) ) / sqRhoOVERJac(:)
     !write(out_unitp,*) ' OpPsiR Grid of PsiR(itab)',OpPsiR
     !write(out_unitp,*) ' OpPsiR Grid of PsiR(itab) : done',itab,OMP_GET_THREAD_NUM() ; flush(out_unitp)

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

   !write(out_unitp,*) '============ END ===============' ; flush(out_unitp)
   !write(out_unitp,*) '================================' ; flush(out_unitp)

  !-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'END ',name_sub
    CALL flush_perso(out_unitp)
  END IF
  !-----------------------------------------------------------
  !CALL UnCheck_mem() ; stop
  END SUBROUTINE sub_TabOpPsi_OF_ONEDP_FOR_SGtype4_ompGrid

  SUBROUTINE get_OpGrid_type10_OF_ONEDP_FOR_SG4(iG,tab_l,para_Op,V,GGiq,sqRhoOVERJac,Jac)
  USE mod_system
  USE mod_nDindex

  USE mod_Coord_KEO,               ONLY: CoordType, get_d0GG
  use mod_PrimOp,                  only: param_d0matop, init_d0matop, Get_iOp_FROM_n_Op,   &
                                         param_typeop, get_d0MatOp_AT_Qact, &
                                         dealloc_tab_of_d0matop

  USE mod_basis_set_alloc,         ONLY : basis
  USE mod_basis,                   ONLY : Rec_Qact_SG4_with_Tab_iq
  USE mod_basis_BtoG_GtoB_SGType4, ONLY : TypeRVec,                     &
                     BDP_TO_GDP_OF_SmolyakRep,GDP_TO_BDP_OF_SmolyakRep, &
                  DerivOp_TO_RDP_OF_SmolaykRep,tabR2grid_TO_tabR1_AT_iG,&
                     getbis_tab_nq,getbis_tab_nb,tabR2gridbis_TO_tabR1_AT_iG
  USE mod_SetOp,                      ONLY : param_Op,write_param_Op


  IMPLICIT NONE

  real (kind=Rkind),  allocatable,    intent(inout)              :: V(:,:,:)
  real (kind=Rkind),  allocatable,    intent(inout), optional    :: GGiq(:,:,:)
  real (kind=Rkind),  allocatable,    intent(inout), optional    :: sqRhoOVERJac(:),Jac(:)
  integer,                            intent(in)                 :: iG,tab_l(:)

  TYPE (param_Op),                    intent(inout)    :: para_Op

  !local variables
  TYPE (CoordType), pointer :: mole
  TYPE(basis),    pointer :: BasisnD

  integer :: iq,nq,D

  integer                               :: derive_termQdyn(2)
  real (kind=Rkind),  allocatable       :: Qact(:)
  real (kind=Rkind)                     :: Rho

  integer,            allocatable       :: tab_nq(:)
  integer,            allocatable       :: tab_iq(:)
  integer :: iterm00,iOp,itabR,nR,nb0,i,j
  integer :: err_sub
  logical :: KEO
  logical :: lformatted=.TRUE.

  TYPE (param_d0MatOp), allocatable :: d0MatOp(:)
  character (len=Line_len)          :: FileName_RV


  !----- for debuging ----------------------------------------------
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub='get_OpGrid_type10_OF_ONEDP_FOR_SG4'
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

   !write(out_unitp,*) '================================' ; flush(out_unitp)
   !write(out_unitp,*) '============ START =============' ; flush(out_unitp)

   mole    => para_Op%mole
   BasisnD => para_Op%BasisnD

   KEO = present(GGiq) .AND. present(sqRhoOVERJac) .AND.                &
         present(Jac)  .AND. para_Op%n_op == 0
   D = size(tab_l)

   CALL alloc_NParray(Qact,(/mole%nb_var/),'Qact',name_sub)
   CALL alloc_NParray(tab_nq,(/ D /),'tab_nq',name_sub)
   CALL alloc_NParray(tab_iq,(/ D /),'tab_iq',name_sub)

   tab_nq(:) = getbis_tab_nq(tab_l,BasisnD%tab_basisPrimSG)

   nq  = BasisnD%para_SGType2%tab_nq_OF_SRep(iG)
   nb0 = BasisnD%para_SGType2%nb0

   IF (KEO) THEN
     CALL alloc_NParray(sqRhoOVERJac,(/ nq /),      'sqRhoOVERJac',name_sub)
     CALL alloc_NParray(Jac,         (/ nq /),      'Jac',         name_sub)
   END IF

   !transfert part of the scalar part of the potential
   iterm00 = para_Op%derive_term_TO_iterm(0,0)
   iOp     = Get_iOp_FROM_n_Op(para_Op%n_Op) ! 2,1,2 + n_Op

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
       CALL alloc_NParray(V,(/ nq,nb0,nb0 /),'V',name_sub)
       V(:,:,:) = ZERO
     ELSE
       IF (para_Op%OpGrid(iterm00)%grid_cte) THEN
         CALL alloc_NParray(V,(/ nq,nb0,nb0 /),'V',name_sub)
         DO i=1,nb0
         DO j=1,nb0
           V(:,j,i) = para_Op%OpGrid(iterm00)%Mat_cte(j,i)
         END DO
         END DO
       ELSE
         CALL alloc_NParray(V,(/ nq,nb0,nb0 /),'V',name_sub)
         IF (associated(para_Op%OpGrid)) THEN
         IF (associated(para_Op%OpGrid(iterm00)%Grid)) THEN
         DO i=1,nb0
         DO j=1,nb0
           CALL tabR2gridbis_TO_tabR1_AT_iG(V(:,j,i),                   &
                                   para_Op%OpGrid(iterm00)%Grid(:,j,i), &
                                            iG,BasisnD%para_SGType2)
         END DO
         END DO
         END IF
         END IF
       END IF
     END IF

   ELSE IF (para_Op%OpGrid(iterm00)%para_FileGrid%Save_FileGrid_done .OR. &
            para_Op%OpGrid(iterm00)%para_FileGrid%Read_FileGrid) THEN
    STOP 'nb0>1 not yet'
    IF (nb0 /= 1) STOP 'nb0>1 not yet'

    !$OMP  CRITICAL (get_OpGrid_type10_OF_ONEDP_FOR_SG4_CRIT2)
     FileName_RV = trim(para_Op%OpGrid(iterm00)%file_Grid%name) // '_SGterm' // int_TO_char(iG)

     IF (iG == 1 .AND. debug) THEN
       write(out_unitp,*) 'iG,nq, FileName_RV',iG,nq,FileName_RV
       CALL flush_perso(out_unitp)
     END IF

     ! this subroutine does not work because V as 3 dim (before 1)
     !CALL sub_ReadRV(V,FileName_RV,lformatted,err_sub)
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

    !$OMP  END CRITICAL (get_OpGrid_type10_OF_ONEDP_FOR_SG4_CRIT2)

   ELSE IF (.NOT. para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_done) THEN
     CALL alloc_NParray(V,(/ nq,nb0,nb0 /),'V',name_sub)
     allocate(d0MatOp(para_Op%para_ReadOp%nb_scalar_Op+2))
     DO i=1,size(d0MatOp)
       CALL Init_d0MatOp(d0MatOp(i),para_Op%param_TypeOp,para_Op%para_ReadOp%nb_elec)
     END DO
     !was
!     DO iOp=1,size(d0MatOp) ! iOp value is change here
!       CALL Init_d0MatOp(d0MatOp(iOp),para_Op%param_TypeOp,para_Op%para_ReadOp%nb_elec)
!     END DO
   END IF

   ! G calculation
   IF (KEO) CALL alloc_NParray(GGiq,(/nq,mole%nb_act1,mole%nb_act1/),'GGiq',name_sub)

   tab_iq(:) = 1 ; tab_iq(1) = 0

   DO iq=1,nq

     CALL ADD_ONE_TO_nDval_m1(tab_iq,tab_nq)

     CALL Rec_Qact_SG4_with_Tab_iq(Qact,BasisnD%tab_basisPrimSG,tab_l,tab_iq,mole,err_sub)

     IF (KEO) THEN
       CALL get_d0GG(Qact,para_Op%para_Tnum,mole,d0GG=GGiq(iq,:,:),     &
                                         def=.TRUE.,Jac=Jac(iq),Rho=Rho)

       sqRhoOVERJac(iq) = sqrt(Rho/Jac(iq))
     END IF

     IF (.NOT. para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_done  .AND. &
         .NOT. para_Op%OpGrid(iterm00)%para_FileGrid%Save_FileGrid_done .AND. &
         .NOT. para_Op%OpGrid(iterm00)%para_FileGrid%Read_FileGrid) THEN

        CALL get_d0MatOp_AT_Qact(Qact,d0MatOp,mole,                     &
                                 para_Op%para_Tnum,para_Op%para_ReadOp%PrimOp_t)

        DO i=1,nb0
        DO j=1,nb0
          V(iq,j,i) = d0MatOp(iOp)%ReVal(j,i,iterm00)
          !was V(iq,j,i) = d0MatOp(iterm00)%ReVal(j,i,1)
        END DO
        END DO

        IF (iG == 1 .AND. debug) THEN
          write(out_unitp,*) 'iG,nq,V',iG,nq,V
          CALL flush_perso(out_unitp)
        END IF

     END IF

   END DO
   !IF (debug) write(out_unitp,*) 'V(:,:,:)',V(:,:,:)


   IF (para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid .AND.       &
     .NOT. para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_done) THEN
   IF (associated(para_Op%OpGrid(iterm00)%Grid)) THEN
     itabR = BasisnD%para_SGType2%tab_Sum_nq_OF_SRep(iG)
     nR    = BasisnD%para_SGType2%tab_nq_OF_SRep(iG)
     para_Op%OpGrid(iterm00)%Grid(itabR-nR+1:itabR,:,:) = V

     IF (iG == 1 .AND. debug) THEN
       write(out_unitp,*) 'iG,nq, save mem V',iG,nq
       CALL flush_perso(out_unitp)
     END IF

   END IF
   END IF

   IF (para_Op%OpGrid(iterm00)%para_FileGrid%Save_FileGrid .AND.           &
      .NOT. para_Op%OpGrid(iterm00)%para_FileGrid%Save_FileGrid_done .AND. &
      .NOT. para_Op%OpGrid(iterm00)%para_FileGrid%Read_FileGrid) THEN
    STOP 'nb0>1 not yet'
    IF (nb0 /= 1) STOP 'nb0>1 not yet'

!    !$OMP  CRITICAL (get_OpGrid_type10_OF_ONEDP_FOR_SG4_CRIT1)
     FileName_RV = trim(para_Op%OpGrid(iterm00)%file_Grid%name) // '_SGterm' // int_TO_char(iG)

     ! this subroutine does not work because V as 3 dim (before 1)
     !CALL sub_WriteRV(V,FileName_RV,lformatted,err_sub)

     IF (iG == 1 .AND. debug) THEN
       write(out_unitp,*) 'iG,nq, save file V',iG,nq
       CALL flush_perso(out_unitp)
     END IF

!    !$OMP  END CRITICAL (get_OpGrid_type10_OF_ONEDP_FOR_SG4_CRIT1)
   END IF

  IF (allocated(tab_nq))       CALL dealloc_NParray(tab_nq,      'tab_nq',      name_sub)
  IF (allocated(tab_iq))       CALL dealloc_NParray(tab_iq,      'tab_iq',      name_sub)
  IF (allocated(Qact))         CALL dealloc_NParray(Qact,        'Qact',        name_sub)

  IF (allocated(d0MatOp)) THEN
    CALL dealloc_Tab_OF_d0MatOp(d0MatOp)
    deallocate(d0MatOp)
  END IF

   !write(out_unitp,*) '============ END ===============' ; flush(out_unitp)
   !write(out_unitp,*) '================================' ; flush(out_unitp)

  !-----------------------------------------------------------
  IF (debug) THEN
    DO i=1,nb0
    DO j=1,nb0
      write(out_unitp,*) 'V(:,j,i) ',j,i,V(:,j,i)
    END DO
    END DO
    write(out_unitp,*) 'GGiq(iq,:,:)'
    DO iq=1,nq
      CALL Write_Mat(GGiq(iq,:,:),out_unitp,6)
    END DO
    write(out_unitp,*) 'END ',name_sub
    CALL flush_perso(out_unitp)
  END IF
  !-----------------------------------------------------------
  END SUBROUTINE get_OpGrid_type10_OF_ONEDP_FOR_SG4

  SUBROUTINE get_OpGrid_type1_OF_ONEDP_FOR_SG4_new(iG,tab_l,para_Op,GridOp)
  USE mod_system
  USE mod_nDindex

  USE mod_Coord_KEO,               ONLY : CoordType
  use mod_PrimOp,                  only: param_d0matop, init_d0matop, Get_iOp_FROM_n_Op,  &
                                         param_typeop, TnumKEO_TO_tab_d0H, get_d0MatOp_AT_Qact, &
                                         dealloc_tab_of_d0matop

  USE mod_basis_set_alloc,         ONLY : basis
  USE mod_basis,                   ONLY : Rec_Qact_SG4_with_Tab_iq
  USE mod_basis_BtoG_GtoB_SGType4, ONLY : TypeRVec,                     &
                     BDP_TO_GDP_OF_SmolyakRep,GDP_TO_BDP_OF_SmolyakRep, &
                  DerivOp_TO_RDP_OF_SmolaykRep,tabR2grid_TO_tabR1_AT_iG,&
                     getbis_tab_nq,getbis_tab_nb,tabR2gridbis_TO_tabR1_AT_iG
  USE mod_SetOp,                      ONLY : param_Op,write_param_Op


  IMPLICIT NONE

  real (kind=Rkind),  allocatable,    intent(inout)              :: GridOp(:,:,:,:)
  integer,                            intent(in)                 :: iG,tab_l(:)

  TYPE (param_Op),                    intent(inout)    :: para_Op

  !local variables
  TYPE (CoordType), pointer :: mole
  TYPE(basis),    pointer :: BasisnD


  integer                               :: derive_termQdyn(2)
  real (kind=Rkind),  allocatable       :: Qact(:)

  integer,            allocatable       :: tab_nq(:)
  integer,            allocatable       :: tab_iq(:)

  integer :: iq,nq,D,iOp,iterm,iterm00,itabR,nR,nb0,i,j
  integer :: err_sub

  logical :: Op_term_done(para_Op%nb_Term)
  TYPE (param_d0MatOp), allocatable :: d0MatOp(:)
  character (len=Line_len) :: FileName_RV


  !----- for debuging ----------------------------------------------
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub='get_OpGrid_type1_OF_ONEDP_FOR_SG4_new'
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

 IF (iG == 1 .AND. debug) write(out_unitp,*) '================================' ; flush(out_unitp)
 IF (iG == 1 .AND. debug) write(out_unitp,*) '============ START =============' ; flush(out_unitp)

   mole    => para_Op%mole
   BasisnD => para_Op%BasisnD


   D = size(tab_l)

   CALL alloc_NParray(Qact,(/mole%nb_var/),'Qact',name_sub)
   CALL alloc_NParray(tab_nq,(/ D /),'tab_nq',name_sub)
   CALL alloc_NParray(tab_iq,(/ D /),'tab_iq',name_sub)

   tab_nq(:) = getbis_tab_nq(tab_l,BasisnD%tab_basisPrimSG)

   nq  = BasisnD%para_SGType2%tab_nq_OF_SRep(iG)
   nb0 = BasisnD%para_SGType2%nb0

   CALL alloc_NParray(GridOp,(/ nq,nb0,nb0,para_Op%nb_Term /),'GridOp',name_sub)
   GridOp(:,:,:,:) = ZERO
   Op_term_done(:) = .FALSE.


   !transfert part of the scalar part of the potential
   iterm00 = para_Op%derive_term_TO_iterm(0,0)
   iOp     = Get_iOp_FROM_n_Op(para_Op%n_Op)
   !write(out_unitp,*) name_sub,' iOp',iOp

   IF (debug) THEN
     write(out_unitp,*) iG,'Save_MemGrid,Save_MemGrid_done',            &
                 para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid,    &
                 para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_done
   END IF

   IF (para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_done) THEN
      IF (associated(para_Op%OpGrid)) THEN
        DO iterm=1,para_Op%nb_Term
          IF (para_Op%OpGrid(iterm)%grid_zero) THEN
            IF (debug) write(out_unitp,*) ' Op term: zero',iG,nq

            GridOp(:,:,:,iterm) = ZERO
            Op_term_done(iterm) = .TRUE.
          ELSE
            IF (para_Op%OpGrid(iterm)%grid_cte) THEN
              IF (debug) write(out_unitp,*) ' Op term: cte',iG,iterm,nq
              DO i=1,nb0
              DO j=1,nb0
                GridOp(:,j,i,iterm) = para_Op%OpGrid(iterm)%Mat_cte(j,i)
              END DO
              END DO
              Op_term_done(iterm) = .TRUE.
            ELSE
              IF (associated(para_Op%OpGrid(iterm)%Grid)) THEN
                IF (debug) write(out_unitp,*) ' Op term: from memory',iG,iterm,nq
                Op_term_done(iterm) = .TRUE.
                DO i=1,nb0
                DO j=1,nb0
                  CALL tabR2gridbis_TO_tabR1_AT_iG(GridOp(:,j,i,iterm),    &
                                        para_Op%OpGrid(iterm)%Grid(:,j,i), &
                                                  iG,BasisnD%para_SGType2)
                END DO
                END DO
              ELSE
                IF (debug) write(out_unitp,*) ' Op term: from memory SRep',iG,iterm,nq
                Op_term_done(iterm) = .TRUE.
                GridOp(:,:,:,iterm) = reshape(para_Op%OpGrid(iterm)%SRep%SmolyakRep(iG)%V,shape=[nq,nb0,nb0])
              END IF
            END IF
          END IF
        END DO
      END IF

   ELSE IF (.NOT. para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_done) THEN
     allocate(d0MatOp(para_Op%para_ReadOp%nb_scalar_Op+2))
     DO i=1,size(d0MatOp)
       CALL Init_d0MatOp(d0MatOp(i),para_Op%param_TypeOp,para_Op%para_ReadOp%nb_elec)
     END DO

     tab_iq(:) = 1 ; tab_iq(1) = 0

     DO iq=1,nq

       CALL ADD_ONE_TO_nDval_m1(tab_iq,tab_nq)

       CALL Rec_Qact_SG4_with_Tab_iq(Qact,BasisnD%tab_basisPrimSG,tab_l,tab_iq,mole,err_sub)

       CALL get_d0MatOp_AT_Qact(Qact,d0MatOp,mole,para_Op%para_Tnum,    &
                                para_Op%para_ReadOp%PrimOp_t)

       IF (para_Op%n_Op == 0) THEN ! H
         CALL TnumKEO_TO_tab_d0H(Qact,d0MatOp(iOp),mole,para_Op%para_Tnum) ! here the vep is added to the potential

         DO iterm=1,para_Op%nb_Term ! just the diagonal elements
           IF (iterm == iterm00) CYCLE ! without the vep+potential
           DO i=1,nb0
             GridOp(iq,i,i,iterm) = d0MatOp(iOp)%ReVal(i,i,iterm)
           END DO
         END DO
       END IF

       ! now the saclar part (diagonal+off diagonal elements)
       ! when Op=H, potential+vep
       GridOp(iq,:,:,iterm00) = d0MatOp(iOp)%ReVal(:,:,iterm00)

     END DO

     IF (iG == 1 .AND. debug) THEN
       write(out_unitp,*) 'iG,nq,GridOp: calc',iG,nq,GridOp(:,:,:,iterm00)
       CALL flush_perso(out_unitp)
     END IF
   END IF

   DO iterm=1,para_Op%nb_Term
     IF (para_Op%OpGrid(iterm)%para_FileGrid%Save_MemGrid .AND.           &
        .NOT. para_Op%OpGrid(iterm)%para_FileGrid%Save_MemGrid_done .AND. &
        Op_term_done(iterm) ) THEN

       IF (associated(para_Op%OpGrid(iterm)%Grid)) THEN
         itabR = BasisnD%para_SGType2%tab_Sum_nq_OF_SRep(iG)
         nR    = BasisnD%para_SGType2%tab_nq_OF_SRep(iG)

         IF (iG == 1 .AND. debug) write(out_unitp,*) 'iG,nq,GridOp: save mem Grid',iterm


         para_Op%OpGrid(iterm00)%Grid(itabR-nR+1:itabR,:,:) = GridOp(:,:,:,iterm)
         para_Op%OpGrid(iterm)%grid_zero = .FALSE.
         para_Op%OpGrid(iterm)%grid_cte  = .FALSE.

       ELSE IF (allocated(para_Op%OpGrid(iterm)%SRep%SmolyakRep)) THEN
         IF (iG == 1 .AND. debug) write(out_unitp,*) 'iG,nq,GridOp SRep: save mem Grid',iterm

         para_Op%OpGrid(iterm)%SRep%SmolyakRep(iG)%V = reshape(GridOp(:,:,:,iterm),shape=[nq*nb0**2])
         para_Op%OpGrid(iterm)%grid_zero = .FALSE.
         para_Op%OpGrid(iterm)%grid_cte  = .FALSE.
       END IF

       IF (iG == 1 .AND. debug) THEN
         write(out_unitp,*) 'iG,nq, save mem GridOp',iG,nq
         CALL flush_perso(out_unitp)
       END IF

     END IF

   END DO


  IF (allocated(tab_nq))       CALL dealloc_NParray(tab_nq,      'tab_nq',      name_sub)
  IF (allocated(tab_iq))       CALL dealloc_NParray(tab_iq,      'tab_iq',      name_sub)
  IF (allocated(Qact))         CALL dealloc_NParray(Qact,        'Qact',        name_sub)

  IF (allocated(d0MatOp)) THEN
    CALL dealloc_Tab_OF_d0MatOp(d0MatOp)
    deallocate(d0MatOp)
  END IF

 IF (iG == 1 .AND. debug)  write(out_unitp,*) '============ END ===============' ; flush(out_unitp)
 IF (iG == 1 .AND. debug)  write(out_unitp,*) '================================' ; flush(out_unitp)

  !-----------------------------------------------------------
  IF (debug) THEN
    DO iterm=1,para_Op%nb_Term
      write(out_unitp,*) 'GridOp(:,:,:,iterm)',iterm,GridOp(:,:,:,iterm)
    END DO
    write(out_unitp,*) 'END ',name_sub
    CALL flush_perso(out_unitp)
  END IF
  !-----------------------------------------------------------
  END SUBROUTINE get_OpGrid_type1_OF_ONEDP_FOR_SG4_new

!=======================================================================================
! new2 get_OpGrid_type1_OF_ONEDP_FOR_SG4 added SmolyakRep part
  SUBROUTINE get_OpGrid_type1_OF_ONEDP_FOR_SG4_new2(iG,tab_l,para_Op,GridOp)
  USE mod_system
  USE mod_nDindex

  USE mod_Coord_KEO,               ONLY : CoordType
  use mod_PrimOp,                  ONLY : param_d0matop, init_d0matop, Get_iOp_FROM_n_Op,  &
                                          param_typeop, TnumKEO_TO_tab_d0H, get_d0MatOp_AT_Qact, &
                                          dealloc_tab_of_d0matop

  USE mod_basis_set_alloc,         ONLY : basis
  USE mod_basis,                   ONLY : Rec_Qact_SG4_with_Tab_iq
  USE mod_basis_BtoG_GtoB_SGType4, ONLY : TypeRVec,                     &
                     BDP_TO_GDP_OF_SmolyakRep,GDP_TO_BDP_OF_SmolyakRep, &
                  DerivOp_TO_RDP_OF_SmolaykRep,tabR2grid_TO_tabR1_AT_iG,&
                     getbis_tab_nq,getbis_tab_nb,tabR2gridbis_TO_tabR1_AT_iG
  USE mod_SetOp,                      ONLY : param_Op,write_param_Op
  USE mod_MPI_aux
  IMPLICIT NONE

  real (kind=Rkind),  allocatable,    intent(inout)              :: GridOp(:,:,:,:)
  integer,                            intent(in)                 :: iG,tab_l(:)

  TYPE (param_Op),                    intent(inout)    :: para_Op

  !local variables
  TYPE (CoordType), pointer :: mole
  TYPE(basis),    pointer :: BasisnD


  integer                               :: derive_termQdyn(2)
  real (kind=Rkind),  allocatable       :: Qact(:)

  integer,            allocatable       :: tab_nq(:)
  integer,            allocatable       :: tab_iq(:)

  integer :: iq,nq,D,iOp,iterm,iterm00,itabR,nR,nb0,i,j
  integer :: err_sub

  logical                           :: KEO_done,Op_term_done(para_Op%nb_Term)
  Logical                           :: Save_MemGrid_done
  TYPE (param_d0MatOp), allocatable :: d0MatOp(:)
  character (len=Line_len) :: FileName_RV


  !----- for debuging ----------------------------------------------
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub='get_OpGrid_type1_OF_ONEDP_FOR_SG4_new2'
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

 IF (iG == 1 .AND. debug) write(out_unitp,*) '================================' ; flush(out_unitp)
 IF (iG == 1 .AND. debug) write(out_unitp,*) '============ START =============' ; flush(out_unitp)

   mole    => para_Op%mole
   BasisnD => para_Op%BasisnD


   D = size(tab_l)

   CALL alloc_NParray(Qact,(/mole%nb_var/),'Qact',name_sub)
   CALL alloc_NParray(tab_nq,(/ D /),'tab_nq',name_sub)
   CALL alloc_NParray(tab_iq,(/ D /),'tab_iq',name_sub)

   tab_nq(:) = getbis_tab_nq(tab_l,BasisnD%tab_basisPrimSG)

   nq  = BasisnD%para_SGType2%tab_nq_OF_SRep(iG)
   nb0 = BasisnD%para_SGType2%nb0

   CALL alloc_NParray(GridOp,(/ nq,nb0,nb0,para_Op%nb_Term /),'GridOp',name_sub)
   GridOp(:,:,:,:) = ZERO
   Op_term_done(:) = .FALSE.


   !transfert part of the scalar part of the potential
   iterm00 = para_Op%derive_term_TO_iterm(0,0)
   iOp     = Get_iOp_FROM_n_Op(para_Op%n_Op)
   !write(out_unitp,*) name_sub,' iOp',iOp

   IF (iG == 1 .AND. debug) THEN
     write(out_unitp,*) iG,'Save_MemGrid,Save_MemGrid_done',            &
                 para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid,    &
                 para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_done
   END IF

   !first deal with constant or zero term.
   KEO_done = .TRUE.
   DO iterm=1,para_Op%nb_Term

     IF (para_Op%OpGrid(iterm)%grid_zero) THEN
       IF (iG == 1 .AND. debug) write(out_unitp,*) ' Op term: zero',iG,nq

       GridOp(:,:,:,iterm) = ZERO
       Op_term_done(iterm) = .TRUE.
     ELSE IF (para_Op%OpGrid(iterm)%grid_cte) THEN
       IF (iG == 1 .AND. debug) write(out_unitp,*) ' Op term: cte',iG,iterm,nq
       DO i=1,nb0
       DO j=1,nb0
           GridOp(:,j,i,iterm) = para_Op%OpGrid(iterm)%Mat_cte(j,i)
       END DO
         END DO
         Op_term_done(iterm) = .TRUE.
     END IF
     IF (iterm00 /= iterm) KEO_done = KEO_done .AND. Op_term_done(iterm)
   END DO
   ! Now KEO_done=.TRUE., when all KEO terms have been treated

   Save_MemGrid_done=para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_done
   IF(openmpi .AND. MPI_scheme/=1)                                                     &
   Save_MemGrid_done=para_Op%OpGrid(1)%para_FileGrid%Save_MemGrid_iG(iG)

   IF (Save_MemGrid_done) THEN

      IF (associated(para_Op%OpGrid)) THEN
        DO iterm=1,para_Op%nb_Term
          IF (Op_term_done(iterm)) CYCLE

          IF (associated(para_Op%OpGrid(iterm)%Grid)) THEN
            IF (iG == 1 .AND. debug) write(out_unitp,*) ' Op term: from memory',iG,iterm,nq
            Op_term_done(iterm) = .TRUE.
            DO i=1,nb0
            DO j=1,nb0
              CALL tabR2gridbis_TO_tabR1_AT_iG(GridOp(:,j,i,iterm),     &
                                     para_Op%OpGrid(iterm)%Grid(:,j,i), &
                                                iG,BasisnD%para_SGType2)
            END DO
            END DO
          ELSE
            IF (iG == 1 .AND. debug) write(out_unitp,*) ' Op term: from memory SRep',iG,iterm,nq
            Op_term_done(iterm) = .TRUE.
            GridOp(:,:,:,iterm) = reshape(para_Op%OpGrid(iterm)%SRep%SmolyakRep(iG)%V,shape=[nq,nb0,nb0])
          END IF

          IF (iG == 1 .AND. debug) THEN
            write(out_unitp,*) 'iG,nq,GridOp: from mem',iG,nq,GridOp(:,:,:,iterm)
            CALL flush_perso(out_unitp)
          END IF

        END DO
      END IF

   ELSE IF (.NOT. Save_MemGrid_done) THEN

     allocate(d0MatOp(para_Op%para_ReadOp%nb_scalar_Op+2))
     DO i=1,size(d0MatOp)
       CALL Init_d0MatOp(d0MatOp(i),para_Op%param_TypeOp,para_Op%para_ReadOp%nb_elec)
     END DO

     tab_iq(:) = 1 ; tab_iq(1) = 0

     DO iq=1,nq

       CALL ADD_ONE_TO_nDval_m1(tab_iq,tab_nq)

       CALL Rec_Qact_SG4_with_Tab_iq(Qact,BasisnD%tab_basisPrimSG,tab_l,tab_iq,mole,err_sub)

       CALL get_d0MatOp_AT_Qact(Qact,d0MatOp,mole,para_Op%para_Tnum,    &
                                para_Op%para_ReadOp%PrimOp_t)

       IF (para_Op%n_Op == 0 .AND. .NOT. KEO_done) THEN ! H
         CALL TnumKEO_TO_tab_d0H(Qact,d0MatOp(iOp),mole,para_Op%para_Tnum) ! here the vep is added to the potential

         DO iterm=1,para_Op%nb_Term ! just the diagonal elements
           IF (iterm == iterm00) CYCLE ! without the vep+potential
           DO i=1,nb0
             GridOp(iq,i,i,iterm) = d0MatOp(iOp)%ReVal(i,i,iterm)
           END DO
         END DO
       END IF

       ! now the saclar part (diagonal+off diagonal elements)
       ! when Op=H, potential+vep
       GridOp(iq,:,:,iterm00) = d0MatOp(iOp)%ReVal(:,:,iterm00)

     END DO

     Op_term_done(:) = .TRUE.


     IF (iG == 1 .AND. debug) THEN
       write(out_unitp,*) 'iG,nq,GridOp: calc',iG,nq,GridOp(:,:,:,iterm00)
       CALL flush_perso(out_unitp)
     END IF
   END IF

   DO iterm=1,para_Op%nb_Term
     Save_MemGrid_done=para_Op%OpGrid(iterm)%para_FileGrid%Save_MemGrid_done
     IF(openmpi .AND. MPI_scheme/=1)                                                   &
     Save_MemGrid_done=para_Op%OpGrid(1)%para_FileGrid%Save_MemGrid_iG(iG)

     IF (para_Op%OpGrid(iterm)%para_FileGrid%Save_MemGrid .AND.                        &
        .NOT. Save_MemGrid_done .AND. Op_term_done(iterm) ) THEN

       ! tag: para_Op%BasisnD%SparseGrid_type==4
       IF (allocated(para_Op%OpGrid(iterm)%SRep%SmolyakRep)) THEN
         IF (iG == 1 .AND. debug) write(out_unitp,*) 'iG,nq,GridOp SRep: save mem Grid',iterm

         para_Op%OpGrid(iterm)%SRep%SmolyakRep(iG)%V = reshape(GridOp(:,:,:,iterm),shape=[nq*nb0**2])
         para_Op%OpGrid(iterm)%grid_zero = .FALSE.
         para_Op%OpGrid(iterm)%grid_cte  = .FALSE.

         IF(openmpi .AND. MPI_scheme/=1 .AND. iterm==para_Op%nb_Term)                  &
                               para_Op%OpGrid(1)%para_FileGrid%Save_MemGrid_iG(iG)=.TRUE.

         IF (iG == 1 .AND. debug) THEN
           write(out_unitp,*) 'iG,nq,GridOp: save mem Grid',iG,nq,para_Op%OpGrid(iterm)%SRep%SmolyakRep(iG)%V
           CALL flush_perso(out_unitp)
         END IF

       ELSE IF (associated(para_Op%OpGrid(iterm)%Grid)) THEN
         itabR = BasisnD%para_SGType2%tab_Sum_nq_OF_SRep(iG)
         nR    = BasisnD%para_SGType2%tab_nq_OF_SRep(iG)

         IF (iG == 1 .AND. debug) write(out_unitp,*) 'iG,nq,GridOp: save mem Grid',iterm

         para_Op%OpGrid(iterm)%Grid(itabR-nR+1:itabR,:,:) = GridOp(:,:,:,iterm)
         para_Op%OpGrid(iterm)%grid_zero = .FALSE.
         para_Op%OpGrid(iterm)%grid_cte  = .FALSE.

         IF (iG == 1 .AND. debug) THEN
           write(out_unitp,*) 'iG,nq,GridOp: save mem Grid',iG,nq,para_Op%OpGrid(iterm)%Grid(itabR-nR+1:itabR,:,:)
           CALL flush_perso(out_unitp)
         END IF

       END IF

     END IF

   END DO


  IF (allocated(tab_nq))       CALL dealloc_NParray(tab_nq,      'tab_nq',      name_sub)
  IF (allocated(tab_iq))       CALL dealloc_NParray(tab_iq,      'tab_iq',      name_sub)
  IF (allocated(Qact))         CALL dealloc_NParray(Qact,        'Qact',        name_sub)

  IF (allocated(d0MatOp)) THEN
    CALL dealloc_Tab_OF_d0MatOp(d0MatOp)
    deallocate(d0MatOp)
  END IF

 IF (iG == 1 .AND. debug)  write(out_unitp,*) '============ END ===============' ; flush(out_unitp)
 IF (iG == 1 .AND. debug)  write(out_unitp,*) '================================' ; flush(out_unitp)

  !-----------------------------------------------------------
  IF (debug) THEN
    !DO iterm=1,para_Op%nb_Term
    !  write(out_unitp,*) 'GridOp(:,:,:,iterm)',iterm,GridOp(:,:,:,iterm)
    !END DO
    write(out_unitp,*) 'END ',name_sub
    CALL flush_perso(out_unitp)
  END IF
  !-----------------------------------------------------------
  END SUBROUTINE get_OpGrid_type1_OF_ONEDP_FOR_SG4_new2
!=======================================================================================

  SUBROUTINE get_OpGrid_type1_OF_ONEDP_FOR_SG4(iG,tab_l,para_Op,GridOp)
  USE mod_system
  USE mod_nDindex

  USE mod_Coord_KEO
  use mod_PrimOp

  USE mod_basis_set_alloc
  USE mod_basis
  USE mod_basis_BtoG_GtoB_SGType4
  USE mod_SetOp


  IMPLICIT NONE

  real(kind=Rkind), allocatable,      intent(inout)    :: GridOp(:,:,:,:)
  integer,                            intent(in)       :: iG,tab_l(:)
  TYPE (param_Op),                    intent(inout)    :: para_Op

  !local variables
  TYPE (CoordType), pointer :: mole
  TYPE(basis),    pointer :: BasisnD

  integer :: iq,nq,D,iOp

  integer                               :: derive_termQdyn(2)
  real (kind=Rkind),  allocatable       :: Qact(:)

  integer,            allocatable       :: tab_nq(:)
  integer,            allocatable       :: tab_iq(:)
  integer :: iterm,iterm00,itabR,nR,nb0,i,j
  integer :: err_sub
  logical :: KEO
  logical :: Op_term_done(para_Op%nb_Term)

  TYPE (param_d0MatOp), allocatable :: d0MatOp(:)


  !----- for debuging ----------------------------------------------
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub='get_OpGrid_type1_OF_ONEDP_FOR_SG4'
  logical, parameter :: debug = .FALSE.
  !logical, parameter :: debug = .TRUE.
  !-----------------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*) 'nb_bie,nb_baie',para_Op%nb_bie,para_Op%nb_baie
    write(out_unitp,*) 'nb_act1',para_Op%mole%nb_act1
    write(out_unitp,*) 'nb_Term',para_Op%nb_Term
    write(out_unitp,*) 'asso OpGrid',associated(para_Op%OpGrid)

    CALL flush_perso(out_unitp)
  END IF
  !-----------------------------------------------------------------

   !write(out_unitp,*) '================================' ; flush(out_unitp)
   !write(out_unitp,*) '============ START =============' ; flush(out_unitp)

   mole    => para_Op%mole
   BasisnD => para_Op%BasisnD

   D = size(tab_l)

   nq  = BasisnD%para_SGType2%tab_nq_OF_SRep(iG)
   nb0 = BasisnD%para_SGType2%nb0
   iterm00 = para_Op%derive_term_TO_iterm(0,0)

   IF (debug) THEN
     write(out_unitp,*) iG,'Save_MemGrid,Save_MemGrid_done',            &
                       para_Op%OpGrid(1)%para_FileGrid%Save_MemGrid,    &
                       para_Op%OpGrid(1)%para_FileGrid%Save_MemGrid_done
   END IF

   CALL alloc_NParray(GridOp,(/ nq,nb0,nb0,para_Op%nb_Term /),'GridOp',name_sub)
   Op_term_done(:) = .FALSE.

   IF (associated(para_Op%OpGrid)) THEN
     DO iterm=1,para_Op%nb_Term
       IF (para_Op%OpGrid(iterm)%grid_zero) THEN
         IF (debug) write(out_unitp,*) ' Op term: zero',iG,nq

         GridOp(:,:,:,iterm) = ZERO
         Op_term_done(iterm) = .TRUE.
       ELSE
         IF (para_Op%OpGrid(iterm)%grid_cte) THEN
           IF (debug) write(out_unitp,*) ' Op term: cte',iG,nq
           DO i=1,nb0
           DO j=1,nb0
             GridOp(:,j,i,iterm) = para_Op%OpGrid(iterm)%Mat_cte(j,i)
           END DO
           END DO
           Op_term_done(iterm) = .TRUE.
         ELSE
           IF (associated(para_Op%OpGrid(iterm)%Grid)) THEN
           IF (debug) write(out_unitp,*) ' Op term: from memory',iG,nq

             Op_term_done(iterm) = .TRUE.
             DO i=1,nb0
             DO j=1,nb0
               CALL tabR2gridbis_TO_tabR1_AT_iG(GridOp(:,j,i,iterm),    &
                                     para_Op%OpGrid(iterm)%Grid(:,j,i), &
                                               iG,BasisnD%para_SGType2)
             END DO
             END DO
           END IF
         END IF
       END IF
     END DO
   END IF

   !Op_term_done(iterm00) = .TRUE.
   IF (.NOT. Op_term_done(iterm00)) THEN
     IF (debug) write(out_unitp,*) ' Op term (potential): calculated',iG,nq
     write(out_unitp,*) ' Op term (potential): calculated',iG,nq

     CALL alloc_NParray(Qact,(/mole%nb_var/),'Qact',name_sub)
     CALL alloc_NParray(tab_nq,(/ D /),'tab_nq',name_sub)
     CALL alloc_NParray(tab_iq,(/ D /),'tab_iq',name_sub)

     tab_nq(:) = getbis_tab_nq(tab_l,BasisnD%tab_basisPrimSG)

     allocate(d0MatOp(para_Op%para_ReadOp%nb_scalar_Op+2))
     DO iOp=1,size(d0MatOp)
       CALL Init_d0MatOp(d0MatOp(iOp),para_Op%param_TypeOp,para_Op%para_ReadOp%nb_elec)
     END DO

     tab_iq(:) = 1 ; tab_iq(1) = 0

     DO iq=1,nq

       CALL ADD_ONE_TO_nDval_m1(tab_iq,tab_nq)

       CALL Rec_Qact_SG4_with_Tab_iq(Qact,BasisnD%tab_basisPrimSG,tab_l,tab_iq,mole,err_sub)

!       IF (KEO) THEN
!         CALL get_d0GG(Qact,para_Op%para_Tnum,mole,d0GG=GGiq(iq,:,:),  &
!                                        def=.TRUE.,Jac=Jac(iq),Rho=Rho)
!       END IF

          CALL get_d0MatOp_AT_Qact(Qact,d0MatOp,mole,para_Op%para_Tnum, &
                                   para_Op%para_ReadOp%PrimOp_t)
          GridOp(iq,:,:,iterm00) = d0MatOp(iterm00)%ReVal(:,:,1)

     END DO

     IF (allocated(tab_nq))  CALL dealloc_NParray(tab_nq,'tab_nq',name_sub)
     IF (allocated(tab_iq))  CALL dealloc_NParray(tab_iq,'tab_iq',name_sub)
     IF (allocated(Qact))    CALL dealloc_NParray(Qact,  'Qact',  name_sub)

     IF (allocated(d0MatOp)) THEN
       CALL dealloc_Tab_OF_d0MatOp(d0MatOp)
       deallocate(d0MatOp)
     END IF

   END IF


  !write(out_unitp,*) '============ END ===============' ; flush(out_unitp)
  !write(out_unitp,*) '================================' ; flush(out_unitp)

  !-----------------------------------------------------------
  IF (debug) THEN
    DO iterm=1,para_Op%nb_Term
      write(out_unitp,*) 'GridOp(:,:,:,iterm)',iterm,GridOp(:,:,:,iterm)
    END DO
    write(out_unitp,*) 'END ',name_sub
    CALL flush_perso(out_unitp)
  END IF
  !-----------------------------------------------------------
  END SUBROUTINE get_OpGrid_type1_OF_ONEDP_FOR_SG4
!=======================================================================================

  SUBROUTINE get_OpGrid_type0_OF_ONEDP_FOR_SG4(iG,tab_l,para_Op,V)
  USE mod_system
  USE mod_nDindex

  USE mod_Coord_KEO,               ONLY: CoordType
  use mod_PrimOp,                  only: param_d0matop, init_d0matop,   &
                                         param_typeop, get_d0MatOp_AT_Qact, &
                                         dealloc_tab_of_d0matop

  USE mod_basis_set_alloc,         ONLY : basis
  USE mod_basis,                   ONLY : Rec_Qact_SG4_with_Tab_iq
  USE mod_basis_BtoG_GtoB_SGType4, ONLY : TypeRVec,                     &
                     BDP_TO_GDP_OF_SmolyakRep,GDP_TO_BDP_OF_SmolyakRep, &
                  DerivOp_TO_RDP_OF_SmolaykRep,tabR2grid_TO_tabR1_AT_iG,&
                     getbis_tab_nq,getbis_tab_nb,tabR2gridbis_TO_tabR1_AT_iG
  USE mod_SetOp,                   ONLY : param_Op,write_param_Op


  IMPLICIT NONE

  real (kind=Rkind),  allocatable,    intent(inout)    :: V(:,:,:)
  integer,                            intent(in)       :: iG,tab_l(:)

  TYPE (param_Op),                    intent(inout)    :: para_Op

  !local variables
  TYPE (CoordType), pointer :: mole
  TYPE(basis),    pointer :: BasisnD

  integer :: iq,nq,D,iOp

  integer                               :: derive_termQdyn(2)
  real (kind=Rkind),  allocatable       :: Qact(:)
  real (kind=Rkind)                     :: Rho

  integer,            allocatable       :: tab_nq(:)
  integer,            allocatable       :: tab_iq(:)
  integer :: iterm00,itabR,nR,nb0,i,j
  integer :: err_sub
  logical :: lformatted=.TRUE.

  TYPE (param_d0MatOp), allocatable :: d0MatOp(:)
  character (len=Line_len) :: FileName_RV


  !----- for debuging ----------------------------------------------
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub='get_OpGrid_type0_OF_ONEDP_FOR_SG4'
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

   !write(out_unitp,*) '================================' ; flush(out_unitp)
   !write(out_unitp,*) '============ START =============' ; flush(out_unitp)

   mole    => para_Op%mole
   BasisnD => para_Op%BasisnD

   D = size(tab_l)

   CALL alloc_NParray(Qact,(/mole%nb_var/),'Qact',name_sub)
   CALL alloc_NParray(tab_nq,(/ D /),'tab_nq',name_sub)
   CALL alloc_NParray(tab_iq,(/ D /),'tab_iq',name_sub)

   tab_nq(:) = getbis_tab_nq(tab_l,BasisnD%tab_basisPrimSG)

   nq  = BasisnD%para_SGType2%tab_nq_OF_SRep(iG)
   nb0 = BasisnD%para_SGType2%nb0


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
       CALL alloc_NParray(V,(/ nq,nb0,nb0 /),'V',name_sub)
       V(:,:,:) = ZERO
     ELSE
       IF (para_Op%OpGrid(iterm00)%grid_cte) THEN
         CALL alloc_NParray(V,(/ nq,nb0,nb0 /),'V',name_sub)
         DO i=1,nb0
         DO j=1,nb0
            V(:,j,i) = para_Op%OpGrid(iterm00)%Mat_cte(j,i)
         END DO
         END DO
       ELSE
         IF (associated(para_Op%OpGrid)) THEN
         IF (associated(para_Op%OpGrid(iterm00)%Grid)) THEN
         CALL alloc_NParray(V,(/ nq,nb0,nb0 /),'V',name_sub)
           DO i=1,nb0
           DO j=1,nb0
             CALL tabR2gridbis_TO_tabR1_AT_iG(V(:,j,i),                 &
                                   para_Op%OpGrid(iterm00)%Grid(:,j,i), &
                                               iG,BasisnD%para_SGType2)
           END DO
           END DO
         END IF
         END IF
       END IF
     END IF
   ELSE IF (para_Op%OpGrid(iterm00)%para_FileGrid%Save_FileGrid_done .OR. &
            para_Op%OpGrid(iterm00)%para_FileGrid%Read_FileGrid) THEN
    STOP 'nb0>1 not yet'

     !$OMP  CRITICAL (get_OpGrid_type0_OF_ONEDP_FOR_SG4_CRIT2)
     FileName_RV = trim(para_Op%OpGrid(iterm00)%file_Grid%name) // '_SGterm' // int_TO_char(iG)

     IF (iG == 1 .AND. debug) THEN
       write(out_unitp,*) 'iG,nq, FileName_RV',iG,nq,FileName_RV
       CALL flush_perso(out_unitp)
     END IF

     ! the rank of V is 3 now, this subroutine work with rank 1
     !CALL sub_ReadRV(V,FileName_RV,lformatted,err_sub)
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
     !$OMP  END CRITICAL (get_OpGrid_type0_OF_ONEDP_FOR_SG4_CRIT2)

   ELSE IF (.NOT. para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_done) THEN
     CALL alloc_NParray(V,(/ nq,nb0,nb0 /),'V',name_sub)
     allocate(d0MatOp(para_Op%para_ReadOp%nb_scalar_Op+2))
     DO iOp=1,size(d0MatOp)
       CALL Init_d0MatOp(d0MatOp(iOp),para_Op%param_TypeOp,para_Op%para_ReadOp%nb_elec)
     END DO
   END IF

   tab_iq(:) = 1 ; tab_iq(1) = 0

   DO iq=1,nq

     CALL ADD_ONE_TO_nDval_m1(tab_iq,tab_nq)

     CALL Rec_Qact_SG4_with_Tab_iq(Qact,BasisnD%tab_basisPrimSG,tab_l,tab_iq,mole,err_sub)

     IF (.NOT. para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_done  .AND. &
         .NOT. para_Op%OpGrid(iterm00)%para_FileGrid%Save_FileGrid_done .AND. &
         .NOT. para_Op%OpGrid(iterm00)%para_FileGrid%Read_FileGrid) THEN

        CALL get_d0MatOp_AT_Qact(Qact,d0MatOp,mole,para_Op%para_Tnum,   &
                                 para_Op%para_ReadOp%PrimOp_t)
        V(iq,:,:) = d0MatOp(iterm00)%ReVal(:,:,1)

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
     para_Op%OpGrid(iterm00)%Grid(itabR-nR+1:itabR,:,:) = V(:,:,:)

     IF (iG == 1 .AND. debug) THEN
       write(out_unitp,*) 'iG,nq, save mem V',iG,nq
       CALL flush_perso(out_unitp)
     END IF

   END IF
   END IF

   IF (para_Op%OpGrid(iterm00)%para_FileGrid%Save_FileGrid .AND.           &
      .NOT. para_Op%OpGrid(iterm00)%para_FileGrid%Save_FileGrid_done .AND. &
      .NOT. para_Op%OpGrid(iterm00)%para_FileGrid%Read_FileGrid) THEN
    STOP 'nb0>1 not yet'

!    !$OMP  CRITICAL (get_OpGrid_type0_OF_ONEDP_FOR_SG4_CRIT1)
     FileName_RV = trim(para_Op%OpGrid(iterm00)%file_Grid%name) // '_SGterm' // int_TO_char(iG)

     ! the rank of V is 3 now, this subroutine work with rank 1
     !CALL sub_WriteRV(V,FileName_RV,lformatted,err_sub)

     IF (iG == 1 .AND. debug) THEN
       write(out_unitp,*) 'iG,nq, save file V',iG,nq
       CALL flush_perso(out_unitp)
     END IF

!    !$OMP  END CRITICAL (get_OpGrid_type0_OF_ONEDP_FOR_SG4_CRIT1)
   END IF

  IF (allocated(tab_nq))       CALL dealloc_NParray(tab_nq,      'tab_nq',      name_sub)
  IF (allocated(tab_iq))       CALL dealloc_NParray(tab_iq,      'tab_iq',      name_sub)
  IF (allocated(Qact))         CALL dealloc_NParray(Qact,        'Qact',        name_sub)

  IF (allocated(d0MatOp)) THEN
    CALL dealloc_Tab_OF_d0MatOp(d0MatOp)
    deallocate(d0MatOp)
  END IF

   !write(out_unitp,*) '============ END ===============' ; flush(out_unitp)
   !write(out_unitp,*) '================================' ; flush(out_unitp)

  !-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'END ',name_sub
    CALL flush_perso(out_unitp)
  END IF
  !-----------------------------------------------------------
  END SUBROUTINE get_OpGrid_type0_OF_ONEDP_FOR_SG4

END MODULE mod_OpPsi_SG4
