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

#if(run_MPI)
PUBLIC :: get_OpGrid_type1_OF_ONEDP_FOR_SG4_MPI
PUBLIC :: sub_TabOpPsi_FOR_SGtype4_SRB_MPI,sub_TabOpPsi_FOR_SGtype4_SRG_MPI
PUBLIC :: ini_iGs_MPI,auto_iGs_MPI
#endif

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
 USE mod_MPI
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

  IF (Psi%cplx) STOP 'cplx'
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
         mod(iG,max(1,BasisnD%para_SGType2%nb_SG/10)) == 0 .AND. MPI_id==0) THEN
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
         mod(iG,max(1,BasisnD%para_SGType2%nb_SG/10)) == 0 .AND. MPI_id==0) THEN
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

 USE mod_Coord_KEO,               ONLY : CoordType, get_Qact, get_d0GG

 USE mod_basis_set_alloc,         ONLY : basis
 USE mod_basis,                   ONLY : Rec_Qact_SG4_with_Tab_iq
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

  USE mod_Coord_KEO,               ONLY : CoordType, get_Qact, get_d0GG

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
  USE mod_MPI_Aid
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

#if(run_MPI)
  integer(kind=MPI_INTEGER_KIND),pointer :: size_PsiR_V(:)
#endif

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
    IF (Psi(1)%cplx) STOP 'cplx'
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

#if(run_MPI)  
  IF(openmpi) THEN
    CALL system_clock(time_point1,time_rate,time_max)
    IF(BasisnD%para_SGType2%once_action) THEN 
      CALL time_perso('MPI loop in action begin')
      CALL ini_iGs_MPI(para_Op,.FALSE.)
      allocate(BasisnD%para_SGType2%size_PsiR_V(0:MPI_np-1))
      DO iterm00=1,para_Op%nb_Term
        CALL allocate_array(para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_iG,     &
                            1,BasisnD%para_SGType2%nb_SG) ! FALSE initially
      ENDDO
      
      allocate(time_MPI_commu(0:MPI_np-1))
      allocate(time_MPI_calcu(0:MPI_np-1))
    ENDIF

!    ! jobs equally assigned to different threads
!    ! the remainder jobs are assigned to each thread from 0
!    nb_per_MPI=BasisnD%para_SGType2%nb_SG/MPI_np
!    !If(mod(BasisnD%para_SGType2%nb_SG,MPI_np)/=0) nb_per_MPI=nb_per_MPI+1
!    nb_rem_MPI=mod(BasisnD%para_SGType2%nb_SG,MPI_np) !remainder jobs 

!    IF(BasisnD%para_SGType2%once_action) THEN
!      allocate(BasisnD%para_SGType2%size_PsiR_V(0:MPI_np-1))
!      !size of %V for each thread
!      Do i_mpi=0,MPI_np-1
!        BasisnD%para_SGType2%size_PsiR_V(i_mpi)=0;
!        !DO iG=i_mpi*nb_per_MPI+1,MIN((i_mpi+1)*nb_per_MPI,BasisnD%para_SGType2%nb_SG)
!        bound1_MPI=i_mpi*nb_per_MPI+1+MIN(i_mpi,nb_rem_MPI)
!        bound2_MPI=(i_mpi+1)*nb_per_MPI+MIN(i_mpi,nb_rem_MPI)+merge(1,0,nb_rem_MPI>i_mpi)
!        DO iG=bound1_MPI,bound2_MPI
!          temp_int=BasisnD%para_SGType2%tab_nb_OF_SRep(iG)*BasisnD%para_SGType2%nb0
!          BasisnD%para_SGType2%size_PsiR_V(i_mpi)=BasisnD%para_SGType2                 &
!                                                         %size_PsiR_V(i_mpi)+temp_int
!        ENDDO
!      ENDDO
!      write(out_unitp,*) 'size_PsiR_V:',BasisnD%para_SGType2%size_PsiR_V(MPI_id),      &
!                         'from',MPI_id
!    ENDIF

    !size of %V for each thread
    size_PsiR_V=>BasisnD%para_SGType2%size_PsiR_V
    IF(BasisnD%para_SGType2%once_action) THEN
      Do i_MPI=0,MPI_np-1
        size_PsiR_V(i_mpi)=0;
        DO iG=iGs_MPI(1,i_MPI),iGs_MPI(2,i_MPI)
          temp_int=BasisnD%para_SGType2%tab_nb_OF_SRep(iG)*BasisnD%para_SGType2%nb0
          size_PsiR_V(i_mpi)=size_PsiR_V(i_mpi)+temp_int
        ENDDO
      ENDDO
      write(out_unitp,*) 'size_PsiR_V:',BasisnD%para_SGType2%size_PsiR_V(MPI_id),      &
                         'from',MPI_id
    ENDIF
    !-----------------------------------------------------------------------------------
    !> action with MPI scheme 1: for more cores, need improvement
    !> action with MPI scheme 2: for less cores
    !> action with MPI scheme 3: simple MPI, available with auto iG distribution
    SELECT CASE(MPI_scheme)
      CASE (1)
        CALL Action_MPI_S1(Psi,OpPsi,BasisnD,para_Op,size_PsiR_V)
      CASE (2)
        CALL Action_MPI_S2(Psi,OpPsi,BasisnD,para_Op,size_PsiR_V)
      CASE (3)
        CALL Action_MPI_S3(Psi,OpPsi,BasisnD,para_Op,size_PsiR_V)
      CASE Default
        IF(size_PsiR_V(0)<1000000) THEN  !< 1000000 according to a few effeiciency test
          CALL Action_MPI_S1(Psi,OpPsi,BasisnD,para_Op,size_PsiR_V)
        ELSE
          CALL Action_MPI_S3(Psi,OpPsi,BasisnD,para_Op,size_PsiR_V)
        ENDIF
    END SELECT
    !-----------------------------------------------------------------------------------

    CALL system_clock(time_point2,time_rate,time_max)
    IF(BasisnD%para_SGType2%once_action) allocate(time_MPI_act_all(0:MPI_np-1))
    time_MPI_act_all(MPI_id)=merge(time_point2-time_point1,                            &
                              time_point2-time_point1+time_max,time_point2>=time_point1)
    time_MPI_action=time_MPI_action+time_MPI_act_all(MPI_id)
    
    IF(BasisnD%para_SGType2%once_action .AND. MPI_id==0)                               &
                   write(out_unitp,*) 'time MPI comm check: ',time_comm,' from ', MPI_id
    IF(BasisnD%para_SGType2%once_action) CALL time_perso('MPI loop in action end') 
    BasisnD%para_SGType2%once_action=.FALSE.

  ELSE
#endif
    ! non openmpi cases
    !---------------------------------------------------------------------------
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
      deallocate(PsiR)
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
#if(run_MPI)
  ENDIF ! for openmpi
#endif

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

#if(run_MPI)
!=======================================================================================
!> sub_TabOpPsi_FOR_SGtype4 for working on full Smolyak rep.
!======================================================================================= 
SUBROUTINE sub_TabOpPsi_FOR_SGtype4_SRB_MPI(Psi,OpPsi,para_Op)
  USE mod_system
  USE mod_nDindex
  USE mod_psi_set_alloc,  ONLY:param_psi
  USE mod_basis_set_alloc,ONLY:basis
  USE mod_SetOp,          ONLY:param_Op
  USE mod_MPI_Aid
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
    CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4_SRB_MPI(iG,psi,OpPsi,                       &
                          BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval(:,iG),para_Op)
  ENDDO

  iterm00 = para_Op%derive_term_TO_iterm(0,0)
  IF (associated(para_Op%OpGrid)) THEN
    para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_done =                          &
                      para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid
    para_Op%OpGrid(iterm00)%para_FileGrid%Save_FileGrid_done =                         &
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

ENDSUBROUTINE sub_TabOpPsi_FOR_SGtype4_SRB_MPI
!=======================================================================================
#endif

#if(run_MPI)
!=======================================================================================
!> sub_TabOpPsi_FOR_SGtype4 for working on full Smolyak rep.
!======================================================================================= 
SUBROUTINE sub_TabOpPsi_FOR_SGtype4_SRG_MPI(Psi,OpPsi,para_Op)
  USE mod_system
  USE mod_nDindex
  USE mod_psi_set_alloc,  ONLY:param_psi
  USE mod_basis_set_alloc,ONLY:basis
  USE mod_SetOp,          ONLY:param_Op
  USE mod_MPI_Aid
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
    CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4_SRG_MPI(iG,psi,OpPsi,                       &
                          BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval(:,iG),para_Op)
  ENDDO

  iterm00 = para_Op%derive_term_TO_iterm(0,0)
  IF (associated(para_Op%OpGrid)) THEN
    para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_done =                          &
                      para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid
    para_Op%OpGrid(iterm00)%para_FileGrid%Save_FileGrid_done =                         &
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

ENDSUBROUTINE sub_TabOpPsi_FOR_SGtype4_SRG_MPI
!=======================================================================================
#endif

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

  USE mod_Coord_KEO,               ONLY : CoordType, get_Qact, get_d0GG

  USE mod_basis_set_alloc,         ONLY : basis
  USE mod_basis,                   ONLY : Rec_Qact_SG4_with_Tab_iq
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

  USE mod_Coord_KEO,               ONLY : CoordType, get_Qact, get_d0GG

  USE mod_basis_set_alloc,         ONLY : basis
  USE mod_basis,                   ONLY : Rec_Qact_SG4_with_Tab_iq
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
  USE mod_MPI_Aid
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
  USE mod_MPI_Aid
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

  USE mod_Coord_KEO,               ONLY : CoordType, get_Qact

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

  USE mod_Coord_KEO,               ONLY : CoordType, get_Qact, get_d0GG

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

  USE mod_Coord_KEO,               ONLY: CoordType, get_Qact, get_d0GG
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
  USE mod_MPI


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
  character (len=Line_len) :: FileName_RV


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
     allocate(d0MatOp(para_Op%para_PES%nb_scalar_Op+2))
     DO i=1,size(d0MatOp) 
       CALL Init_d0MatOp(d0MatOp(i),para_Op%param_TypeOp,para_Op%para_PES%nb_elec)
     END DO
     !was
!     DO iOp=1,size(d0MatOp) ! iOp value is change here
!       CALL Init_d0MatOp(d0MatOp(iOp),para_Op%param_TypeOp,para_Op%para_PES%nb_elec)
!     END DO
   END IF

   ! G calculation
   IF (KEO) CALL alloc_NParray(GGiq,(/nq,mole%nb_act1,mole%nb_act1/),'GGiq',name_sub)

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

        CALL get_d0MatOp_AT_Qact(Qact,d0MatOp,mole,                     &
                                 para_Op%para_Tnum,para_Op%para_PES)

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

  USE mod_Coord_KEO,               ONLY : CoordType, get_Qact
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
     allocate(d0MatOp(para_Op%para_PES%nb_scalar_Op+2))
     DO i=1,size(d0MatOp)
       CALL Init_d0MatOp(d0MatOp(i),para_Op%param_TypeOp,para_Op%para_PES%nb_elec)
     END DO

     tab_iq(:) = 1 ; tab_iq(1) = 0

     DO iq=1,nq

       CALL ADD_ONE_TO_nDval_m1(tab_iq,tab_nq)

       CALL get_Qact(Qact,mole%ActiveTransfo) ! rigid, flexible coordinates
       CALL Rec_Qact_SG4_with_Tab_iq(Qact,BasisnD%tab_basisPrimSG,tab_l,tab_iq,mole,err_sub)

       CALL get_d0MatOp_AT_Qact(Qact,d0MatOp,mole,                     &
                                   para_Op%para_Tnum,para_Op%para_PES)

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

  USE mod_Coord_KEO,               ONLY : CoordType, get_Qact
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
  USE mod_MPI_Aid
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

  logical :: KEO_done,Op_term_done(para_Op%nb_Term)
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

#if(run_MPI)
   IF(para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_iG(iG)) THEN
#else
   IF (para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_done) THEN
#endif
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

#if(run_MPI)
   ELSE IF(.NOT. para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_iG(iG)) THEN
#else
   ELSE IF (.NOT. para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_done) THEN
#endif
     allocate(d0MatOp(para_Op%para_PES%nb_scalar_Op+2))
     DO i=1,size(d0MatOp)
       CALL Init_d0MatOp(d0MatOp(i),para_Op%param_TypeOp,para_Op%para_PES%nb_elec)
     END DO

     tab_iq(:) = 1 ; tab_iq(1) = 0

     DO iq=1,nq

       CALL ADD_ONE_TO_nDval_m1(tab_iq,tab_nq)

       CALL get_Qact(Qact,mole%ActiveTransfo) ! rigid, flexible coordinates
       CALL Rec_Qact_SG4_with_Tab_iq(Qact,BasisnD%tab_basisPrimSG,tab_l,tab_iq,mole,err_sub)

       CALL get_d0MatOp_AT_Qact(Qact,d0MatOp,mole,                     &
                                   para_Op%para_Tnum,para_Op%para_PES)

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
#if(run_MPI)
     IF(para_Op%OpGrid(iterm)%para_FileGrid%Save_MemGrid .AND.                         &
        .NOT. para_Op%OpGrid(iterm)%para_FileGrid%Save_MemGrid_iG(iG) .AND.            &
        Op_term_done(iterm)) THEN
#else
     IF (para_Op%OpGrid(iterm)%para_FileGrid%Save_MemGrid .AND.           &
        .NOT. para_Op%OpGrid(iterm)%para_FileGrid%Save_MemGrid_done .AND. &
        Op_term_done(iterm) ) THEN
#endif
       IF (associated(para_Op%OpGrid(iterm)%Grid)) THEN
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

       ELSE IF (allocated(para_Op%OpGrid(iterm)%SRep%SmolyakRep)) THEN
         IF (iG == 1 .AND. debug) write(out_unitp,*) 'iG,nq,GridOp SRep: save mem Grid',iterm

         !> @warning array allocation upon assignment, only valid from Fortran 2003  
         para_Op%OpGrid(iterm)%SRep%SmolyakRep(iG)%V = reshape(GridOp(:,:,:,iterm),shape=[nq*nb0**2])
         para_Op%OpGrid(iterm)%grid_zero = .FALSE.
         para_Op%OpGrid(iterm)%grid_cte  = .FALSE.
#if(run_MPI)
         para_Op%OpGrid(iterm)%para_FileGrid%Save_MemGrid_iG(iG)=.TRUE.
#endif
         IF (iG == 1 .AND. debug) THEN
           write(out_unitp,*) 'iG,nq,GridOp: save mem Grid',iG,nq,para_Op%OpGrid(iterm)%SRep%SmolyakRep(iG)%V
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

     allocate(d0MatOp(para_Op%para_PES%nb_scalar_Op+2))
     DO iOp=1,size(d0MatOp)
       CALL Init_d0MatOp(d0MatOp(iOp),para_Op%param_TypeOp,para_Op%para_PES%nb_elec)
     END DO

     tab_iq(:) = 1 ; tab_iq(1) = 0

     DO iq=1,nq

       CALL ADD_ONE_TO_nDval_m1(tab_iq,tab_nq)

       CALL get_Qact(Qact,mole%ActiveTransfo) ! rigid, flexible coordinates
       CALL Rec_Qact_SG4_with_Tab_iq(Qact,BasisnD%tab_basisPrimSG,tab_l,tab_iq,mole,err_sub)

!       IF (KEO) THEN
!         CALL get_d0GG(Qact,para_Op%para_Tnum,mole,d0GG=GGiq(iq,:,:),  &
!                                        def=.TRUE.,Jac=Jac(iq),Rho=Rho)
!       END IF

          CALL get_d0MatOp_AT_Qact(Qact,d0MatOp,mole,                   &
                                   para_Op%para_Tnum,para_Op%para_PES)
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

  USE mod_Coord_KEO,               ONLY: CoordType, get_Qact
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
     allocate(d0MatOp(para_Op%para_PES%nb_scalar_Op+2))
     DO iOp=1,size(d0MatOp)
       CALL Init_d0MatOp(d0MatOp(iOp),para_Op%param_TypeOp,para_Op%para_PES%nb_elec)
     END DO
   END IF

   tab_iq(:) = 1 ; tab_iq(1) = 0

   DO iq=1,nq

     CALL ADD_ONE_TO_nDval_m1(tab_iq,tab_nq)

     CALL get_Qact(Qact,mole%ActiveTransfo) ! rigid, flexible coordinates
     CALL Rec_Qact_SG4_with_Tab_iq(Qact,BasisnD%tab_basisPrimSG,tab_l,tab_iq,mole,err_sub)

     IF (.NOT. para_Op%OpGrid(iterm00)%para_FileGrid%Save_MemGrid_done  .AND. &
         .NOT. para_Op%OpGrid(iterm00)%para_FileGrid%Save_FileGrid_done .AND. &
         .NOT. para_Op%OpGrid(iterm00)%para_FileGrid%Read_FileGrid) THEN

        CALL get_d0MatOp_AT_Qact(Qact,d0MatOp,mole,                     &
                                 para_Op%para_Tnum,para_Op%para_PES)
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
  

#if(run_MPI)
!=======================================================================================  
!> action with MPI: scheme 1
!=======================================================================================  
  SUBROUTINE Action_MPI_S1(Psi,OpPsi,BasisnD,para_Op,size_PsiR_V)
    USE mod_system
    USE mod_nDindex
    USE mod_Coord_KEO,                ONLY : CoordType
    USE mod_basis_set_alloc,          ONLY : basis
    USE mod_basis_BtoG_GtoB_SGType4,  ONLY : tabPackedBasis_TO_tabR_AT_iG, &
                                             tabR_AT_iG_TO_tabPackedBasis, &
                                             TypeRVec,dealloc_TypeRVec,    &
                                             PackedBasis_TO_tabR_index_MPI,&
                                             tabR_TO_tabPackedBasis_MPI,   &
                                             tabPackedBasis_TO_tabR_MPI
    USE mod_psi,                      ONLY : param_psi
    USE mod_SetOp,                    ONLY : param_Op
    USE mod_MPI
    USE mod_MPI_Aid
    IMPLICIT NONE
    
    TYPE(param_psi),                       intent(in)    :: Psi(:)
    TYPE(param_psi),                       intent(inout) :: OpPsi(:)
    TYPE(basis),pointer,                   intent(inout) :: BasisnD
    Integer(kind=MPI_INTEGER_KIND),pointer,intent(inout) :: size_PsiR_V(:) 
    TYPE(param_Op),                        intent(inout) :: para_Op

    TYPE(TypeRVec),allocatable                           :: PsiR(:)
    Real(kind=Rkind),allocatable                         :: all_RvecB_temp(:)
    Real(kind=Rkind),allocatable                         :: all_RvecB_temp2(:)
    Integer(kind=MPI_INTEGER_KIND),pointer               :: Psi_size_MPI0
    Integer(kind=MPI_INTEGER_KIND),pointer               :: reduce_Vlength_master(:)
    Integer(kind=MPI_INTEGER_KIND),pointer               :: reduce_Vlength
    Integer,pointer                                      :: Max_nDI_ib0
    Integer,pointer                                      :: V_allcount
    Integer,pointer                                      :: V_allcount2
    Integer                                              :: iG
    integer                                              :: iG_MPI
    Integer                                              :: ii
    Integer                                              :: itab
    Logical,pointer                                      :: once_action


    Psi_size_MPI0  => BasisnD%para_SGType2%Psi_size_MPI0  
    If(MPI_id==0) Psi_size_MPI0=INT(size(Psi),MPI_INTEGER_KIND) ! for itab
    CALL MPI_BCAST(Psi_size_MPI0,size1_MPI,MPI_int_def,root_MPI,MPI_COMM_WORLD,MPI_err)

    Max_nDI_ib0    => BasisnD%para_SGType2%Max_nDI_ib0
    reduce_Vlength => BasisnD%para_SGType2%reduce_Vlength
    V_allcount     => BasisnD%para_SGType2%V_allcount
    V_allcount2    => BasisnD%para_SGType2%V_allcount2
    once_action    => BasisnD%para_SGType2%once_action
    
    ! calculate total length of vectors for each threads--------------------------------
    IF(once_action .AND. MPI_id==0) THEN
      write(out_unitp,*) 'action with MPI: Scheme 1'
      allocate(BasisnD%para_SGType2%nDI_index_master(0:MPI_np-1))
      allocate(BasisnD%para_SGType2%reduce_Vlength_master(0:MPI_np-1))
    ENDIF
    reduce_Vlength_master=>BasisnD%para_SGType2%reduce_Vlength_master

    If(MPI_id==0) Max_nDI_ib0=size(psi(1)%RvecB)/BasisnD%para_SGType2%nb0
    CALL MPI_BCAST(Max_nDI_ib0,size1_MPI,MPI_int_fortran,root_MPI,                     &
                   MPI_COMM_WORLD,MPI_err)
         
    ! clean PsiR     
    IF(allocated(PsiR)) deallocate(PsiR)
    allocate(PsiR(Psi_size_MPI0))
    
    ! works on the other threads
    !-----------------------------------------------------------------------------------
    ! all information for reduce_Vlength_MPI is keeped on master
    ! each thread keep its own    
    IF(MPI_id/=0) THEN
      !-generate index,  once only------------------------------------------------------
      IF(once_action) THEN
        ! initialize the size for index for pack psi on each threads
        BasisnD%para_SGType2%num_nDI_index=size_PsiR_V(MPI_id)/20
        reduce_Vlength=0
        V_allcount=0
        CALL allocate_array(BasisnD%para_SGType2%nDI_index,                            &
                            1,BasisnD%para_SGType2%num_nDI_index)
        CALL allocate_array(BasisnD%para_SGType2%nDI_index_list,1,size_PsiR_V(MPI_id))

        DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)
          ! note size(psi(i)%RvecB) same for all i 
          CALL PackedBasis_TO_tabR_index_MPI(iG,BasisnD%para_SGType2,reduce_Vlength,   &
                                         BasisnD%para_SGType2%nDI_index,Max_nDI_ib0,   &
                                         BasisnD%para_SGType2%nDI_index_list)
        ENDDO
        write(out_unitp,*) 'V_allcount check',V_allcount,size_PsiR_V(MPI_id),          &
                                              reduce_Vlength,' from ',MPI_id
        CALL MPI_Send(reduce_Vlength,size1_MPI,MPI_int_fortran,root_MPI,MPI_id,        &
                      MPI_COMM_WORLD,MPI_err)
        CALL MPI_Send(BasisnD%para_SGType2%nDI_index(1:reduce_Vlength),                &
                      reduce_Vlength,MPI_int_fortran,root_MPI,MPI_id,                  &
                      MPI_COMM_WORLD,MPI_err)
      ENDIF
      
      !> check length of integer
      IF(Int(reduce_Vlength,8)*Int(Psi_size_MPI0,8)>huge(0_4)                          &
         .AND. MPI_INTEGER_KIND==4) THEN
        STOP 'integer exceed 32-bit MPI, use 64-bit MPI instead'
      ENDIF
      
      !-wait for master-----------------------------------------------------------------
      CALL allocate_array(all_RvecB_temp,1,reduce_Vlength*Psi_size_MPI0)
      CALL allocate_array(all_RvecB_temp2,1,reduce_Vlength*Psi_size_MPI0)
      Call MPI_Recv(all_RvecB_temp,reduce_Vlength*Psi_size_MPI0,MPI_REAL8,root_MPI,    &
                    MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)

      !-calculation on threads----------------------------------------------------------
      V_allcount=0
      V_allcount2=0   

      DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)
        ! all_RvecB_temp --> PsiR
        CALL tabPackedBasis_TO_tabR_MPI(PsiR,all_RvecB_temp,iG,BasisnD%para_SGType2,   &
                        BasisnD%para_SGType2%nDI_index,reduce_Vlength,Psi_size_MPI0,   &
                        Max_nDI_ib0,BasisnD%para_SGType2%nDI_index_list)

        !-main calculation--------------------------------------------------------------
        CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4(PsiR,iG,                                &
                        BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval(:,iG),para_Op)

        !-pack PsiR in the slave threads------------------------------------------------
        CALL tabR_TO_tabPackedBasis_MPI(all_RvecB_temp2,PsiR,iG,                       &
                         BasisnD%para_SGType2,BasisnD%WeightSG(iG),                    &
                         BasisnD%para_SGType2%nDI_index,reduce_Vlength,Psi_size_MPI0,  &
                         Max_nDI_ib0,BasisnD%para_SGType2%nDI_index_list)
      ENDDO 
      
      !-send packed PsiR to master------------------------------------------------------
      CALL MPI_Send(all_RvecB_temp2,reduce_Vlength*Psi_size_MPI0,MPI_REAL8,root_MPI,   &
                    MPI_id,MPI_COMM_WORLD,MPI_err)
                   
    ENDIF ! for MPI_id/=0  

    ! works on the master threads
    !-----------------------------------------------------------------------------------
    ! save index information for all threads on master & send vectors to other threads
    IF(MPI_id==0) THEN
      ! works for the other threads-----------------------------------------------------
      Do i_mpi=1,MPI_np-1
        !-wait for other threads--------------------------------------------------------
        CALL time_record(time_comm,time_temp1,time_temp2,1)
        IF(once_action) THEN
          CALL MPI_Recv(reduce_Vlength,size1_MPI,MPI_int_fortran,i_mpi,i_mpi,          &
                        MPI_COMM_WORLD,MPI_stat,MPI_err)
          reduce_Vlength_master(i_mpi)=reduce_Vlength
          allocate(BasisnD%para_SGType2%nDI_index_master(i_mpi)%array(reduce_Vlength))
          CALL MPI_Recv(BasisnD%para_SGType2%nDI_index_master(i_mpi)%array,            &
                        reduce_Vlength,MPI_int_fortran,i_mpi,i_mpi,                    &
                        MPI_COMM_WORLD,MPI_stat,MPI_err)
          write(out_unitp,*) 'length of comm list:',reduce_Vlength_master(i_mpi),      &
                             'from',i_mpi
        ENDIF ! for once_action
        CALL time_record(time_comm,time_temp1,time_temp2,2)

        ! pack vectores to send
        CALL allocate_array(all_RvecB_temp,1,reduce_Vlength_master(i_mpi)*Psi_size_MPI0)
        DO itab=1,Psi_size_MPI0
          DO ii=1,reduce_Vlength_master(i_mpi)
            temp_int=(itab-1)*reduce_Vlength_master(i_mpi)+ii
            all_RvecB_temp(temp_int)=psi(itab)%RvecB(                                  &
                                 BasisnD%para_SGType2%nDI_index_master(i_mpi)%array(ii))
          ENDDO
        ENDDO

        CALL time_record(time_comm,time_temp1,time_temp2,1)
        CALL MPI_Send(all_RvecB_temp,Psi_size_MPI0*reduce_Vlength_master(i_mpi),       &
                      MPI_REAL8,i_mpi,i_mpi,MPI_COMM_WORLD,MPI_err)
        CALL time_record(time_comm,time_temp1,time_temp2,2)

      ENDDO ! for i_mpi=1,MPI_np-1

      !-calculation on master-----------------------------------------------------------
      DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)
        DO itab=1,Psi_size_MPI0
          CALL tabPackedBasis_TO_tabR_AT_iG(PsiR(itab)%V,psi(itab)%RvecB,              &
                                            iG,BasisnD%para_SGType2)
        ENDDO 

        !-main calculation--------------------------------------------------------------
        CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4(PsiR,iG,                                &
                        BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval(:,iG),para_Op) 
              
        DO itab=1,Psi_size_MPI0
          CALL tabR_AT_iG_TO_tabPackedBasis(OpPsi(itab)%RvecB,PsiR(itab)%V,iG,         &
                                            BasisnD%para_SGType2,BasisnD%WeightSG(iG))
        ENDDO 

      ENDDO ! for iG

      !-receive results from slave threads----------------------------------------------
      Do i_mpi=1,MPI_np-1
        CALL allocate_array(all_RvecB_temp,1,Psi_size_MPI0*reduce_Vlength_master(i_mpi))
        CALL time_record(time_comm,time_temp1,time_temp2,1)
        CALL MPI_Recv(all_RvecB_temp,Psi_size_MPI0*reduce_Vlength_master(i_mpi),       &
                      MPI_REAL8,i_mpi,i_mpi,MPI_COMM_WORLD,MPI_stat,MPI_err)
        CALL time_record(time_comm,time_temp1,time_temp2,2)

        !-extract results from other threads--------------------------------------------
        DO itab=1,Psi_size_MPI0
          Do ii=1,reduce_Vlength_master(i_mpi)
            temp_int=BasisnD%para_SGType2%nDI_index_master(i_mpi)%array(ii)
            OpPsi(itab)%RvecB(temp_int)=OpPsi(itab)%RvecB(temp_int)                    &
                               +all_RvecB_temp((itab-1)*reduce_Vlength_master(i_mpi)+ii)
          ENDDO 
        ENDDO
      ENDDO ! for i_mpi=1,MPI_np-1
    ENDIF ! for MPI_id==0
                  
    IF(allocated(all_RvecB_temp))  deallocate(all_RvecB_temp)
    IF(allocated(all_RvecB_temp2)) deallocate(all_RvecB_temp2)
    DO itab=1,Psi_size_MPI0
      CALL dealloc_TypeRVec(PsiR(itab))
    END DO
    
  END SUBROUTINE Action_MPI_S1
!=======================================================================================  
#endif

#if(run_MPI)
!=======================================================================================  
!> action with MPI: scheme 3
!======================================================================================= 
  SUBROUTINE Action_MPI_S3(Psi,OpPsi,BasisnD,para_Op,size_PsiR_V)
    USE mod_system
    USE mod_nDindex
    USE mod_Coord_KEO,                ONLY : CoordType
    USE mod_basis_set_alloc,          ONLY : basis
    USE mod_basis_BtoG_GtoB_SGType4,  ONLY : tabPackedBasis_TO_tabR_AT_iG, &
                                             tabR_AT_iG_TO_tabPackedBasis, &
                                             TypeRVec,dealloc_TypeRVec,    &
                                             tabR_TO_tabPackedBasis_MPI,   &
                                             tabPackedBasis_TO_tabR_MPI

    USE mod_psi,                      ONLY : param_psi
    USE mod_SetOp,                    ONLY : param_Op
    USE mod_MPI_Aid
    IMPLICIT NONE

    TYPE(param_psi),                       intent(in)    :: Psi(:)
    TYPE(param_psi),                       intent(inout) :: OpPsi(:)
    TYPE(basis),pointer,                   intent(inout) :: BasisnD
    Integer(kind=MPI_INTEGER_KIND),pointer,intent(inout) :: size_PsiR_V(:) 
    TYPE(param_Op),                        intent(inout) :: para_Op

    TYPE(TypeRVec),allocatable                           :: PsiR(:)
    Real(kind=Rkind),allocatable                         :: PsiR_temp(:) 
    Integer(kind=MPI_INTEGER_KIND)                       :: PsiR_temp_length(0:MPI_np-1)
    Integer(kind=MPI_INTEGER_KIND),pointer               :: Psi_size_MPI0
    Integer                                              :: iG
    integer                                              :: iG_MPI
    Integer                                              :: PsiR_count1
    Integer                                              :: PsiR_count2
    Integer                                              :: PsiR_count_iG
    Integer                                              :: PsiR_V_iG_size
    Integer                                              :: itab
    Integer                                              :: time_auto1
    Integer                                              :: time_auto2
    Logical,pointer                                      :: once_action


    once_action => BasisnD%para_SGType2%once_action

    ! auto distribute Smolyak terms to different threads 
    IF(.NOT. once_action .AND. allocated(iGs_MPI)) THEN
      write(out_unitp,*) '------------------------------------------------'
      write(out_unitp,*) 'iGs_MPI at current step: ',iGs_MPI
      IF(MPI_np>1) CALL auto_iGs_MPI(para_Op)
      write(out_unitp,*) 'iGs_MPI at new step    : ',iGs_MPI
      write(out_unitp,*) '------------------------------------------------'
      
      ! get new size_PsiR_V according to new iGs_MPI
      Do i_MPI=0,MPI_np-1
        size_PsiR_V(i_mpi)=0;
        DO iG=iGs_MPI(1,i_MPI),iGs_MPI(2,i_MPI)
          temp_int=BasisnD%para_SGType2%tab_nb_OF_SRep(iG)*BasisnD%para_SGType2%nb0
          size_PsiR_V(i_mpi)=size_PsiR_V(i_mpi)+temp_int
        ENDDO
      ENDDO
    ENDIF

    IF(once_action) THEN
      IF(MPI_id==0) write(out_unitp,*) 'action with MPI: Scheme 3'
    ENDIF

    ! set iGs_MPI
!    IF(once_action) THEN
!      DO i_mpi=0,MPI_np-1
!        bound1_MPI=i_mpi*nb_per_MPI+1+MIN(i_mpi,nb_rem_MPI)
!        bound2_MPI=(i_mpi+1)*nb_per_MPI+MIN(i_mpi,nb_rem_MPI)                          &
!                                       +merge(1,0,nb_rem_MPI>i_mpi)
!        iGs_MPI(1,i_mpi)=bound1_MPI
!        iGs_MPI(2,i_mpi)=bound2_MPI
!      ENDDO
!    ENDIF

    Psi_size_MPI0 => BasisnD%para_SGType2%Psi_size_MPI0  
    If(MPI_id==0) Psi_size_MPI0=INT(size(Psi),MPI_INTEGER_KIND) ! for itab
    CALL MPI_BCAST(Psi_size_MPI0,size1_MPI,MPI_int_def,root_MPI,MPI_COMM_WORLD,MPI_err)

    IF(allocated(PsiR)) deallocate(PsiR)
    allocate(PsiR(Psi_size_MPI0))
    !-----------------------------------------------------------------------------------
    IF(MPI_id==0) THEN
      !-prepare PsiR(itab)%V to be send to other threads--------------------------------
      time_MPI_commu(MPI_id)=0
      CALL time_record(time_MPI_commu(MPI_id),time_auto1,time_auto2,1)
      DO i_mpi=1,MPI_np-1
        CALL allocate_array(PsiR_temp,1,Psi_size_MPI0*size_PsiR_V(i_mpi))
        PsiR_count1=0
        DO iG_MPI=iGs_MPI(1,i_mpi),iGs_MPI(2,i_mpi)
          DO itab=1,Psi_size_MPI0
            CALL tabPackedBasis_TO_tabR_AT_iG(PsiR(itab)%V,psi(itab)%RvecB,iG_MPI,     &
                                              BasisnD%para_SGType2)
            PsiR_count2=PsiR_count1+size(PsiR(itab)%V)
            PsiR_temp(PsiR_count1+1:PsiR_count2)=PsiR(itab)%V
            PsiR_count1=PsiR_count2
          ENDDO
        ENDDO ! for iG_MPI
        PsiR_temp_length(i_mpi)=PsiR_count1
        write(out_unitp,*) 'length check:',size_PsiR_V(i_mpi),size(psi(1)%RvecB)
        
        ! double check
        IF(Abs(Psi_size_MPI0*size_PsiR_V(i_mpi)-PsiR_temp_length(i_mpi))>0) THEN
          write(out_unitp,*) 'error in MPI action part,check length'
          STOP
        ENDIF
        
        CALL time_record(time_comm,time_temp1,time_temp2,1)
        CALL MPI_Send(PsiR_temp,Psi_size_MPI0*size_PsiR_V(i_mpi),MPI_real_fortran,     &
                      i_mpi,i_mpi,MPI_COMM_WORLD,MPI_err)
        CALL time_record(time_comm,time_temp1,time_temp2,2)
      ENDDO
      CALL time_record(time_MPI_commu(MPI_id),time_auto1,time_auto2,2)

      !---------------------------------------------------------------------------------

      !-calculations on MPI_id=0--------------------------------------------------------
      DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)
        !-transfore to SRep-------------------------------------------------------------
        DO itab=1,Psi_size_MPI0
          CALL tabPackedBasis_TO_tabR_AT_iG(PsiR(itab)%V,psi(itab)%RvecB,iG,           &
                                            BasisnD%para_SGType2)
        ENDDO

        !-main calculation--------------------------------------------------------------
        CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4(PsiR,iG,                                &
                          BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval(:,iG),para_Op)  

        !-back to compact form----------------------------------------------------------
        DO itab=1,Psi_size_MPI0
          CALL tabR_AT_iG_TO_tabPackedBasis(OpPsi(itab)%RvecB,PsiR(itab)%V,iG,         &
                                            BasisnD%para_SGType2,BasisnD%WeightSG(iG))
        ENDDO
      ENDDO ! main loop of iG for calcuation on master

      !-receive results from other threads----------------------------------------------
      CALL time_record(time_MPI_commu(MPI_id),time_auto1,time_auto2,1)
      DO i_mpi=1,MPI_np-1
        CALL time_record(time_comm,time_temp1,time_temp2,1)
        If(allocated(PsiR_temp)) deallocate(PsiR_temp)
        allocate(PsiR_temp(Psi_size_MPI0*size_PsiR_V(i_mpi)))
        !CALL allocate_array(PsiR_temp,1,Psi_size_MPI0*size_PsiR_V(i_mpi))
        CALL MPI_Recv(PsiR_temp,Psi_size_MPI0*size_PsiR_V(i_mpi),MPI_REAL8,i_mpi,      &
                      i_mpi,MPI_COMM_WORLD,MPI_stat,MPI_err)
        CALL time_record(time_comm,time_temp1,time_temp2,2)

        PsiR_count1=0
        DO iG_MPI=iGs_MPI(1,i_mpi),iGs_MPI(2,i_mpi)
          PsiR_V_iG_size=BasisnD%para_SGType2%tab_nb_OF_SRep(iG_MPI)                   &
                        *BasisnD%para_SGType2%nb0
          DO itab=1,Psi_size_MPI0
            PsiR_count2=PsiR_count1+PsiR_V_iG_size
            IF(allocated(PsiR(itab)%V)) deallocate(PsiR(itab)%V)
            allocate(PsiR(itab)%V(PsiR_V_iG_size))
            PsiR(itab)%V=PsiR_temp(PsiR_count1+1:PsiR_count2)
            PsiR_count1=PsiR_count2
            CALL tabR_AT_iG_TO_tabPackedBasis(OpPsi(itab)%RvecB,PsiR(itab)%V,          &
                                   iG_MPI,BasisnD%para_SGType2,BasisnD%WeightSG(iG_MPI))
          ENDDO
        ENDDO
      ENDDO  ! for i_mpi=1,MPI_np-1  
      CALL time_record(time_MPI_commu(MPI_id),time_auto1,time_auto2,2)
    ENDIF ! for MPI_id==0
    
    !-calculation on other threads------------------------------------------------------
    IF(MPI_id/=0) THEN
      time_MPI_commu(MPI_id)=0
      allocate(PsiR_temp(Psi_size_MPI0*size_PsiR_V(MPI_id)))
      ! waiting for master, get PsiR(itab)%V
      CALL time_record(time_MPI_commu(MPI_id),time_auto1,time_auto2,1)
      CALL MPI_Recv(PsiR_temp,Psi_size_MPI0*size_PsiR_V(MPI_id),MPI_REAL8,root_MPI,    &
                    MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
      CALL time_record(time_MPI_commu(MPI_id),time_auto1,time_auto2,2)
                              
      !-loop for main calculation-------------------------------------------------------
      PsiR_count1=0
      DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)
        PsiR_count_iG=PsiR_count1
        PsiR_V_iG_size=BasisnD%para_SGType2%tab_nb_OF_SRep(iG)*BasisnD%para_SGType2%nb0

        !-extract SRep from PsiR_temp---------------------------------------------------
        DO itab=1,Psi_size_MPI0
          PsiR_count2=PsiR_count1+PsiR_V_iG_size
          PsiR(itab)%V=PsiR_temp(PsiR_count1+1:PsiR_count2)
          PsiR_count1=PsiR_count2
        ENDDO

        !-main calculation--------------------------------------------------------------
        CALL sub_TabOpPsi_OF_ONEDP_FOR_SGtype4(PsiR,iG,                                &
                        BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval(:,iG),para_Op)  

        !-pack and PsiR(itab)%V---------------------------------------------------------
        Do itab=1,Psi_size_MPI0
          PsiR_temp(PsiR_count_iG+1:PsiR_count_iG+size(PsiR(itab)%V))=PsiR(itab)%V
          PsiR_count_iG=PsiR_count_iG+size(PsiR(itab)%V)
        ENDDO
      ENDDO ! for iG
      
      CALL time_record(time_MPI_commu(MPI_id),time_auto1,time_auto2,1)
      !-send back PsiR(itab)%V----------------------------------------------------------
      CALL MPI_Send(PsiR_temp,Psi_size_MPI0*size_PsiR_V(MPI_id),MPI_REAL8,             &
                    root_MPI,MPI_id,MPI_COMM_WORLD,MPI_err)
      CALL time_record(time_MPI_commu(MPI_id),time_auto1,time_auto2,2)
    ENDIF ! for MPI_id/=0
    !-----------------------------------------------------------------------------------

    If(allocated(PsiR_temp)) deallocate(PsiR_temp)
    DO itab=1,Psi_size_MPI0
      CALL dealloc_TypeRVec(PsiR(itab))
    END DO
       
  END SUBROUTINE Action_MPI_S3
!======================================================================================= 
#endif

#if(run_MPI)
!=======================================================================================  
!> action with MPI: scheme 2
!======================================================================================= 
  SUBROUTINE Action_MPI_S2(Psi,OpPsi,BasisnD,para_Op,size_PsiR_V)
    USE mod_system
    USE mod_nDindex
    USE mod_Coord_KEO,                ONLY : CoordType
    USE mod_basis_set_alloc,          ONLY : basis
    USE mod_basis_BtoG_GtoB_SGType4,  ONLY : tabPackedBasis_TO_tabR_AT_iG, &
                                             tabR_AT_iG_TO_tabPackedBasis, &
                                             TypeRVec,dealloc_TypeRVec,    &
                                             tabR_TO_tabPackedBasis_MPI,   &
                                             tabPackedBasis_TO_tabR_MPI

    USE mod_psi,                      ONLY : param_psi
    USE mod_SetOp,                    ONLY : param_Op
    USE mod_MPI_Aid
    IMPLICIT NONE
    
    TYPE(param_psi),                       intent(in)    :: Psi(:)
    TYPE(param_psi),                       intent(inout) :: OpPsi(:)
    TYPE(basis),pointer,                   intent(inout) :: BasisnD
    Integer(kind=MPI_INTEGER_KIND),pointer,intent(inout) :: size_PsiR_V(:) 
    TYPE(param_Op),                        intent(inout) :: para_Op
    
    TYPE(TypeRVec),allocatable                           :: PsiR(:)
    Real(kind=Rkind),allocatable                         :: PsiR_temp(:) 
    Integer(kind=MPI_INTEGER_KIND)                       :: PsiR_temp_length(0:MPI_np-1)
    Integer(kind=MPI_INTEGER_KIND),pointer               :: Psi_size_MPI0
    Integer                                              :: iG
    integer                                              :: iG_MPI
    Integer                                              :: itab
    Integer                                              :: size_RvecB
    Logical,pointer                                      :: once_action
!    
!    
!    once_action   => BasisnD%para_SGType2%once_action
!    Psi_size_MPI0 => BasisnD%para_SGType2%Psi_size_MPI0  
!
!    IF(once_action) THEN
!      IF(MPI_id==0) write(out_unitp,*) 'action with MPI: Scheme 2'
!    ENDIF
!    
!    If(MPI_id==0) THEN
!      Psi_size_MPI0=INT(size(Psi),MPI_INTEGER_KIND) ! for itab
!      size_RvecB=size(psi(1)%RvecB)
!    ENDIF
!    CALL MPI_BCAST(Psi_size_MPI0,size1_MPI,MPI_int_def,root_MPI,MPI_COMM_WORLD,MPI_err)
!    CALL MPI_BCAST(size_RvecB,   size1_MPI,MPI_int_def,root_MPI,MPI_COMM_WORLD,MPI_err)
!
!    IF(MPI_id==0) THEN
!      DO i_MPI=1,MPI_np-1
!        DO itab=1,Psi_size_MPI0
!          CALL MPI_Send(psi(itab)%RvecB,size_RvecB,)
!          CALL MPI_Send(reduce_Vlength,size1_MPI,MPI_int_fortran,root_MPI,MPI_id,        &
!                        MPI_COMM_WORLD,MPI_err)
!        ENDDO
!      ENDDO
!    ENDIF
!    
!    IF(MPI_id/=0) THEN
!    ENDIF
!    
  END SUBROUTINE Action_MPI_S2
!======================================================================================= 
#endif

#if(run_MPI)
!=======================================================================================
!> initialize iGs_MPI for each threads
!> note to .false. para_Op%BasisnD%para_SGType2%once_action to perform once
!=======================================================================================
  SUBROUTINE ini_iGs_MPI(para_Op,once)
    USE mod_SetOp,ONLY:param_Op
    USE mod_MPI_Aid
    IMPLICIT NONE
    
    TYPE(param_Op),                 intent(in)    :: para_Op
    Logical,                        intent(in)    :: once
    
    nb_per_MPI=para_Op%BasisnD%para_SGType2%nb_SG/MPI_np
    nb_rem_MPI=mod(para_Op%BasisnD%para_SGType2%nb_SG,MPI_np) 
    
    CALL allocate_array(iGs_MPI,1,2,0,MPI_np-1)
    
    DO i_MPI=0,MPI_np-1
      bound1_MPI=i_MPI*nb_per_MPI+1+MIN(i_MPI,nb_rem_MPI)
      bound2_MPI=(i_MPI+1)*nb_per_MPI+MIN(i_MPI,nb_rem_MPI)+merge(1,0,nb_rem_MPI>i_MPI)
      iGs_MPI(1,i_MPI)=bound1_MPI
      iGs_MPI(2,i_MPI)=bound2_MPI
    ENDDO
    IF(once) para_Op%BasisnD%para_SGType2%once_action=.FALSE.
    
  ENDSUBROUTINE ini_iGs_MPI  
!=======================================================================================
#endif

#if(run_MPI)
!=======================================================================================
!> auto adjust the distribution of iG on each threads according to 
!> time used in the previous action
!=======================================================================================
  SUBROUTINE auto_iGs_MPI_old(para_Op)
    USE mod_system
    USE mod_SetOp,ONLY:param_Op
    USE mod_MPI_Aid
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


    CALL MPI_collect_info(time_MPI_act_all) ! time_MPI_act_all infor on master

    IF(MPI_id==0) THEN 
      write(*,*) 'time used in action for each threads:',time_MPI_act_all
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
      master_commu_time=MAX(0.,Real(time_MPI_act_all(0)-time_MPI_act_all(1),        &
                                                                            kind=Rkind))
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
            write(*,*) 'iG_threads checkcheck',temp_int,iG_threads(i_MPI)
            
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

  ENDSUBROUTINE auto_iGs_MPI_old
!=======================================================================================
#endif

#if(run_MPI)
!=======================================================================================
!> auto adjust the distribution of iG on each threads according to 
!> time used in the previous action
!=======================================================================================
  SUBROUTINE auto_iGs_MPI_old2(para_Op)
    USE mod_system
    USE mod_SetOp,ONLY:param_Op
    USE mod_MPI_Aid
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


    CALL MPI_collect_info(time_MPI_act_all) ! time_MPI_act_all infor on master

    returnall=.FALSE.
    IF(MPI_id==0) THEN 
      write(*,*) 'time used in action for each threads:',time_MPI_act_all
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

  ENDSUBROUTINE auto_iGs_MPI_old2
!=======================================================================================
#endif

#if(run_MPI)
!=======================================================================================
!> auto adjust the distribution of iG on each threads according to 
!> time used in the previous action
!=======================================================================================
  SUBROUTINE auto_iGs_MPI(para_Op)
    USE mod_system
    USE mod_SetOp,ONLY:param_Op
    USE mod_MPI_Aid
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


    CALL MPI_collect_info(time_MPI_act_all) ! time_MPI_act_all infor on master
    CALL MPI_collect_info(time_MPI_commu)
     
    returnall=.FALSE.
    IF(MPI_id==0) THEN 
      write(out_unitp,*) 'time used in action for each threads:',time_MPI_act_all
      write(out_unitp,*) 'commu time used in action for each threads:',time_MPI_commu
      ave_time=Real(SUM(time_MPI_act_all),kind=Rkind)/Real(MPI_np,kind=Rkind)
      IF(ALL(ABS(ave_time-time_MPI_act_all)/ave_time<0.12)) returnall=.TRUE.
      write(out_unitp,*) 'time balance:',ABS(ave_time-time_MPI_act_all)/ave_time<0.12
    ENDIF

    CALL MPI_BCAST(returnall,size1_MPI,MPI_Logical,root_MPI,MPI_COMM_WORLD,MPI_err)
    IF(returnall) RETURN
    
    IF(MPI_id==0) THEN
      time_MPI_calcu=time_MPI_act_all-time_MPI_commu
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
          count_time=time_MPI_commu(i_MPI)
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

  ENDSUBROUTINE auto_iGs_MPI
!=======================================================================================
#endif

#if(run_MPI)
!=======================================================================================
!> @brief subroutine for recording time 
!> @param time_sum should be initialized before calling this function
!=======================================================================================
  SUBROUTINE time_record(time_sum,time1,time2,point)
    USE mod_MPI
    IMPLICIT NONE

    Integer,                        intent(inout) :: time_sum
    Integer,                        intent(inout) :: time1
    Integer,                        intent(inout) :: time2
    Integer,                           intent(in) :: point

    IF(point==1) THEN
      CALL system_clock(time1,time_rate,time_max)
    ELSEIF(point==2) THEN
      CALL system_clock(time2,time_rate,time_max)
      time_sum=time_sum+merge(time2-time1,time2-time1+time_max,time2>=time1)
    ELSE
      STOP 'error when calling time_record'
    ENDIF
  ENDSUBROUTINE
#endif


END MODULE mod_OpPsi_SG4
