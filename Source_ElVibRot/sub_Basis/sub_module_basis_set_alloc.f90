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
 MODULE mod_basis_set_alloc
  use mod_system

  use mod_dnSVM,   only: type_dnmat,type_dncplxmat,type_intvec,          &
                         alloc_array, alloc_dncplxmat, alloc_dnmat,      &
                         dealloc_dnmat, dealloc_dncplxmat, dealloc_array,&
                         dealloc_intvec, sub_intvec1_to_intvec2,         &
                         write_dnsvm, write_dnmat, write_dncplxmat
  use mod_nDindex, only: type_ndindex, write_ndindex, dealloc_ndindex,  &
                         alloc_array, alloc_nparray,                    &
                         dealloc_nparray, dealloc_array

      use mod_RotBasis_Param ! all
      use mod_Basis_Grid_Param
      USE mod_SymAbelian
      USE mod_param_SGType2
      USE mod_Basis_L_TO_n
      IMPLICIT NONE

      PRIVATE

        TYPE basis
          logical           :: active                = .FALSE.    ! (F)
          logical           :: print_info_OF_basisDP = .TRUE.     ! (T)

          integer           :: ndim        = 0                    !  dimension of the basis
                                                                  !  i.e. nb of variables
          integer, allocatable :: iQdyn(:)                        !  dim : ndim
          integer, allocatable :: Tabder_Qdyn_TO_Qbasis(:)        ! Tabder_Qdyn_TO_Qbasis(0:nb_var)

          TYPE (Type_SymAbelian), pointer :: P_SymAbelian => null() ! for the abelian symmetry

          integer, allocatable :: tab_ndim_index(:,:)   ! tab_ndim_index(ndim,nb). This table is defined only
                                                        ! when the basis is packed

          integer                        :: nb = 0                            !  nb of basis functions
          integer                        :: nb_init = 0                       !  nb of basis functions (before contraction)
          TYPE (Basis_Grid_Param)        :: Basis_Grid_Para

          real (kind=Rkind), allocatable :: EneH0(:)     ! EeneH0(nb) : EeneH0(ib)=<d0b(:,ib) I H0 I d0b(:,ib)>

          TYPE (Type_dnMat)              :: dnRGB     ! basis functions d0b(nq,nb) ....
          TYPE (Type_dnCplxMat)          :: dnCGB     ! basis functions d0cb(nq,nb)
          TYPE (Type_dnMat)              :: dnRBG     ! basis functions td0b(nb,nq)  (transpose)
          TYPE (Type_dnCplxMat)          :: dnCBG     ! basis functions td0cb(nb,nq) (transpose)
          TYPE (Type_dnMat)              :: dnRBGwrho ! basis functions td0b(nq,nb)  (transpose+wrho)
          TYPE (Type_dnCplxMat)          :: dnCBGwrho ! basis functions td0cb(nq,nb) (transpose+wrho)


          ! the projection of d1b and d2b => d0b
          logical                        :: dnBBRep      = .FALSE. ! (F) If we calculate the BasisRep representation
          logical                        :: dnBBRep_done = .FALSE. ! (F) If we calculate the BasisRep representation

          TYPE (Type_dnMat)              :: dnRBB ! matrices which enables to tranform d./dQ d2./dQ2 on the basis
          TYPE (Type_dnCplxMat)          :: dnCBB ! matrices which enables to tranform d./dQ d2./dQ2 on the basis

          logical                        :: dnGGRep      = .FALSE. ! (T) If we calculate the dnRGG
          logical                        :: dnGGRep_done = .FALSE. ! (F) If we calculate the dnRGG
          TYPE (Type_dnMat)              :: dnRGG ! matrices which enables to tranform d./dQ d2./dQ2 on the grid

          logical                        :: cplx = .FALSE. !  .T. if the basis set is complex (def .F.)

          integer                        :: nq_max_Nested = -1        ! Value to calculate Nested grid point (with Nested=1)
                                                                     ! With the value -1, the value is automatically defined
          integer                        :: Nested        =  0        ! When value > 0, we use Nested Grid
          real (kind=Rkind), allocatable :: x(:,:)  !  grid points x(ndim,nq)
          real (kind=Rkind), allocatable :: w(:)    !  weight w(nq)
          real (kind=Rkind), allocatable :: rho(:)  !  rho : wrho(nq)
          real (kind=Rkind), allocatable :: wrho(:) ! weight * rho : wrho(nq)
          integer,           allocatable :: nrho(:) ! nrho(dim), to define the volume element for Tnum

          logical                        :: check_basis            = .TRUE.   ! if T, the basis set is checked (ortho ...)
          logical                        :: check_nq_OF_basis      = .TRUE.   ! if T, the nq is adapted to nb (nq>= nb)
          logical                        :: packed                 = .FALSE.  ! packed=.T. if the basis set is packed (true nD basis)
          logical                        :: packed_done            = .FALSE.  ! packed_done=.T., if the basis has been packed
                                                                          ! Remark: A nD-contracted basis is allways packed

          logical                        :: primitive              = .FALSE.  ! IF True, the basis set is a primitive basis
          logical                        :: primitive_done         = .FALSE.  ! This parameter is used to check if primitive basis-sets or set-up
                                                                              ! for direct-product basis-set

          logical                        :: auto_basis             = .FALSE.  ! it is done automatically

          logical                        :: contrac                = .FALSE.  !  .T. if the basis set is contracted
          logical                        :: auto_contrac           = .FALSE.  ! it is done automatically
          logical                        :: contrac_analysis       = .FALSE.  ! perform contraction, but only for the analysis (RD)
          real (kind=Rkind)              :: max_ene_contrac        = ONETENTH ! maximal energy to select nbc (in ua)
          integer                        :: max_nbc                = 0
          integer                        :: min_nbc                = 0
          integer                        :: auto_contrac_type1_TO  = 100      ! type 1 coordinates of other basis are changed
          integer                        :: auto_contrac_type21_TO = 100      ! type 1 coordinates of other basis are changed
          logical                        :: POGridRep              = .FALSE.  ! used PO-GridRep (auto_contrac=t)
          logical                        :: POGridRep_polyortho    = .FALSE.  ! used PO-GridRep with othonormal poly (auto_contrac=t)
          integer                        :: nqPLUSnbc_TO_nqc       = 0        ! grid points add to nbc to get nqc (POGridRep)
          logical                        :: xPOGridRep_done        = .FALSE.  ! (F). If T construct a basis with a x already present
          integer                        :: nbc                    = 0        ! nb of basis functions after contraction
          integer                        :: nqc                    = 0        !  nb of grid points after contraction
          logical                        :: make_cubature          = .FALSE.
          logical                        :: Restart_make_cubature  = .FALSE.
          logical                        :: read_contrac_file      = .FALSE.  ! .T. if the basis set is contracted
          TYPE(param_file)               :: file_contrac                      ! file for read contraction coef
          real (kind=Rkind), allocatable :: Rvec(:,:)                     ! real eigenvectors for the contraction

          integer                        :: type      = 0     ! basis type
          character (len=Name_len)       :: name      = "0"   ! name of the basis set
          real (kind=Rkind), allocatable :: A(:)              ! Range [A:B]
          real (kind=Rkind), allocatable :: B(:)              ! Range [A:B]
          real (kind=Rkind), allocatable :: Q0(:)             ! scaling parameters
          real (kind=Rkind), allocatable :: scaleQ(:)         ! scaling parameters
          integer, allocatable           :: opt_A(:)          ! integer for the optimisation
          integer, allocatable           :: opt_B(:)          ! integer for the optimisation
          integer, allocatable           :: opt_Q0(:)         ! integer for the optimisation
          integer, allocatable           :: opt_scaleQ(:)     ! integer for the optimisation
          integer                        :: nb_transfo = 0    ! basis type
          real (kind=Rkind), allocatable :: cte_transfo(:,:)  ! for nb_transfo
          integer                        :: opt_param = 0     ! number of parameters to be optimized



          integer                        :: nb_basis                 = 0
          logical                        :: tab_basis_done           = .FALSE. ! When it TRUE the tab_Pbasis(:) are set up
          logical                        :: tab_basis_linked         = .FALSE. ! When it TRUE the tab_Pbasis(:)%Pbasis points to other basis

          TYPE (P_basis), pointer        :: tab_Pbasis(:)            => null() !  tab_Pbasis(nb_basis)
          TYPE (Type_IntVec), pointer    :: Tab_OF_Tabnb2(:)         => null() ! Tab_OF_Tabnb2(nb_basis), for SparseBasis or Pruned basis
          TYPE (Type_nDindex)            :: nDindG                             ! enable to use multidimensional index for the grid
          TYPE (Type_nDindex), pointer   :: nDindB                   => null() ! enable to use multidimensional index for the basis functions
          TYPE (Type_nDindex), pointer   :: nDindB_uncontracted      => null() ! enable to use multidimensional index for the uncontracted basis functions
          integer                        :: Type_OF_nDindB           = 1       ! enable to chose the initialization of nDindB
          integer                        :: nDinit_OF_nDindB         = 1       ! enable to chose nDinit of nDindB
          real (kind=Rkind)              :: Norm_OF_nDindB           = huge(1) ! Norm for the initialization of nDindB
          real (kind=Rkind)              :: weight_OF_nDindB         = ONE     ! weight for the initialization of nDindB
          integer                        :: MaxCoupling_OF_nDindB    = -1      ! number of coupling modes (default all)
          integer                        :: nb_OF_MinNorm_OF_nDindB  = 1
          integer                        :: Div_nb_TO_Norm_OF_nDindB = 1
          logical                        :: contrac_WITH_nDindB      = .FALSE.

          integer                        :: SparseGrid_type          = 0 ! 0 no sparse grid
                                                                         ! 1 old several nD-grids
                                                                         ! 2 new sparse grid with ONE nD-grids

          logical                        :: SparseGrid_With_Cuba    = .TRUE. ! When 2 or more are true, the program choses the optimal one
          logical                        :: SparseGrid_With_Smolyak = .TRUE. ! When only one is true, the program tries to use only one
          logical                        :: SparseGrid_With_DP      = .TRUE. ! Remark: when only SparseGrid_With_Cuba=T, and the grid does not exit the program stops

          TYPE (Basis_L_TO_n)         :: L_TO_nq
          TYPE (Basis_L_TO_n)         :: L_TO_nb

          integer                     :: L_SparseGrid             = -1      ! parameter for the number of points of SparseGrid

          logical                     :: nqSG_SMALLER_nqDP        = .TRUE.

          integer                     :: L_SparseBasis            = -1      ! parameter for the number of points of SparseGrid

          logical                     :: With_L                   = .FALSE.

          integer                     :: nb_SG                    = 0       ! number of direct-product grids (SparseGrid)
          real (kind=Rkind), allocatable :: WeightSG(:)                     ! WeightSG(nb_SG)
          TYPE (P_basis), pointer     :: tab_PbasisSG(:)          => null() ! tab_PbasisSG(nb_SG) fort SGtype=1
          TYPE (basis), pointer       :: tab_basisPrimSG(:,:)     => null() ! tab_basis(nb_basis,0:Lmax)

          TYPE (param_SGType2)        :: para_SGType2

          TYPE (RotBasis_Param)       :: RotBasis


        END TYPE basis

        TYPE P_basis
          TYPE (basis), pointer :: Pbasis => null()
        END TYPE P_basis

        TYPE param_AllBasis
          logical                     :: alloc      = .FALSE.

          integer                     :: nb_be      = 1
          integer                     :: nb_bi      = 1

          TYPE (basis), pointer       :: BasisnD    => null()
          TYPE (basis), pointer       :: Basis2n    => null()
          TYPE (basis), pointer       :: BasisElec  => null()
          TYPE (basis), pointer       :: BasisRot   => null()
        END TYPE param_AllBasis

      INTERFACE alloc_array
        MODULE PROCEDURE alloc_array_OF_Basisdim0,alloc_array_OF_Basisdim1
        MODULE PROCEDURE alloc_array_OF_Basisdim2

        MODULE PROCEDURE alloc_array_OF_P_Basisdim1,alloc_array_OF_P_Basisdim2

      END INTERFACE

      INTERFACE dealloc_array
        MODULE PROCEDURE dealloc_array_OF_Basisdim0,dealloc_array_OF_Basisdim1
        MODULE PROCEDURE dealloc_array_OF_Basisdim2

        MODULE PROCEDURE dealloc_array_OF_P_Basisdim1,dealloc_array_OF_P_Basisdim2
      END INTERFACE

      PUBLIC  basis, alloc_init_basis, dealloc_basis, clean_basis, &
              basis2TObasis1, RecWrite_basis, RecWriteMini_basis,  &
              alloc_dnb_OF_basis, dealloc_dnb_OF_basis,            &
              alloc_xw_OF_basis, dealloc_xw_OF_basis
      PUBLIC  get_x_OF_basis, get_x_AT_iq_OF_basis
      PUBLIC  get_rho_OF_basis, get_rho_AT_iq_OF_basis
      PUBLIC  get_w_OF_basis, get_w_AT_iq_OF_basis
      PUBLIC  get_wrho_OF_basis, get_wrho_AT_iq_OF_basis
      PUBLIC  Get_MatdnRGG, Get_MatdnRGB, Get_MatdnRBB
      PUBLIC  Get_MatdnCGB, Get_MatdnCBB
      PUBLIC  Get2_MatdnRGB, Get2_MatdnRBG
      PUBLIC  Set_nq_OF_basis, get_nq_FROM_basis, get_nqa_FROM_basis, get_nb_FROM_basis
      PUBLIC  get_tab_nq_OF_Qact, get_nb_bi_FROM_AllBasis
      PUBLIC  get_nb_be_FROM_basis

      PUBLIC  P_basis, alloc_tab_Pbasis_OF_basis
      PUBLIC  param_AllBasis, alloc_AllBasis, dealloc_AllBasis

      PUBLIC alloc_array, dealloc_array

      CONTAINS

!      ==========================================================
!       alloc memory of the TYPE basis
!      ==========================================================
      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
       SUBROUTINE alloc_init_basis(basis_set)
       TYPE (basis),         intent(inout) :: basis_set

       integer      :: i
       integer      :: err_mem,memory
       character (len=*), parameter :: name_sub='alloc_init_basis'


!       IF (basis_set%ndim < 1 .AND. .NOT. NewBasisEl) THEN
!         write(out_unitp,*) ' ERROR in ',name_sub
!         write(out_unitp,*) '  WRONG parameter values: ndim',basis_set%ndim
!         write(out_unitp,*) '  CHECK the fortran !!'
!         STOP
!       END IF

       IF (.NOT. associated(basis_set%nDindB)) THEN
         CALL alloc_array(basis_set%nDindB,"basis_set%nDindB",name_sub)
       END IF

       IF (basis_set%ndim > 0) THEN
         IF (basis_set%ndim > 0 .AND. .NOT. allocated(basis_set%iQdyn) ) THEN
           CALL alloc_NParray(basis_set%iQdyn,(/ basis_set%ndim /),       &
                           "basis_set%iQdyn",name_sub)
           basis_set%iQdyn(:) = 0
         END IF
         IF (basis_set%nb_Transfo > 0 .AND. .NOT. allocated(basis_set%cte_Transfo) ) THEN
           CALL alloc_NParray(basis_set%cte_Transfo,(/ 20,basis_set%nb_Transfo /),&
                             "basis_set%cte_Transfo",name_sub)
           basis_set%cte_Transfo(:,:) = ZERO
         END IF

         IF (basis_set%ndim > 0 .AND. .NOT. allocated(basis_set%A) ) THEN
           CALL alloc_NParray(basis_set%A,(/ basis_set%ndim /),           &
                           "basis_set%A",name_sub)
           basis_set%A(:) = ZERO
         END IF
         IF (basis_set%ndim > 0 .AND. .NOT. allocated(basis_set%B) ) THEN
           CALL alloc_NParray(basis_set%B,(/ basis_set%ndim /),           &
                           "basis_set%B",name_sub)
           basis_set%B(:) = ZERO
         END IF
         IF (basis_set%ndim > 0 .AND. .NOT. allocated(basis_set%Q0) ) THEN
           CALL alloc_NParray(basis_set%Q0,(/ basis_set%ndim /),          &
                           "basis_set%Q0",name_sub)
           basis_set%Q0(:) = ZERO
         END IF
         IF (basis_set%ndim > 0 .AND. .NOT. allocated(basis_set%scaleQ) ) THEN
           CALL alloc_NParray(basis_set%scaleQ,(/ basis_set%ndim /),      &
                           "basis_set%scaleQ",name_sub)
           basis_set%scaleQ(:) = ONE
         END IF

         IF (basis_set%ndim > 0 .AND. .NOT. allocated(basis_set%opt_A) ) THEN
           CALL alloc_NParray(basis_set%opt_A,(/ basis_set%ndim /),       &
                           "basis_set%opt_A",name_sub)
           basis_set%opt_A(:) = 0
         END IF
         IF (basis_set%ndim > 0 .AND. .NOT. allocated(basis_set%opt_B) ) THEN
           CALL alloc_NParray(basis_set%opt_B,(/ basis_set%ndim /),       &
                           "basis_set%opt_B",name_sub)
           basis_set%opt_B(:) = 0
         END IF
         IF (basis_set%ndim > 0 .AND. .NOT. allocated(basis_set%opt_Q0) ) THEN
           CALL alloc_NParray(basis_set%opt_Q0,(/ basis_set%ndim /),      &
                           "basis_set%opt_Q0",name_sub)
           basis_set%opt_Q0(:) = 0
         END IF
         IF (basis_set%ndim > 0 .AND. .NOT. allocated(basis_set%opt_scaleQ) ) THEN
           CALL alloc_NParray(basis_set%opt_scaleQ,(/ basis_set%ndim /),  &
                           "basis_set%opt_scaleQ",name_sub)
           basis_set%opt_scaleQ(:) = 0
         END IF

         IF (basis_set%ndim > 0 .AND. .NOT. allocated(basis_set%nrho) ) THEN
           CALL alloc_NParray(basis_set%nrho,(/ basis_set%ndim /),        &
                           "basis_set%nrho",name_sub)
           basis_set%nrho(:) = 1 ! with nrho(i)=1, the voume element is dT = dQi
         END IF
       END IF


       END SUBROUTINE alloc_init_basis

       SUBROUTINE alloc_dnb_OF_basis(basis_set)
         TYPE(basis), intent(inout) :: basis_set
         integer           :: i,nq,nderiv_loc


         !----- for debuging --------------------------------------------------
         integer :: err_mem,memory
         character (len=*), parameter :: name_sub='alloc_dnb_OF_basis'
         logical,parameter :: debug=.FALSE.
         !logical,parameter :: debug=.TRUE.
         !----- for debuging --------------------------------------------------

         IF (debug) write(out_unitp,*) 'BEGINNING ',name_sub
         CALL flush_perso(out_unitp)

         nq = get_nq_FROM_basis(basis_set)

         nderiv_loc = 2
         IF (basis_set%ndim == 0) THEN
           nderiv_loc = 0
           write(out_unitp,*) ' WARNING in ',name_sub
           write(out_unitp,*) '  ndim = 0'
         END IF

         IF (basis_set%ndim < 0 .OR. basis_set%nb < 1 .OR. nq < 1) THEN
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) '  WRONG paramter values: ndim or nb or nq',      &
                       basis_set%ndim,basis_set%nb,nq
           write(out_unitp,*) '  CHECK the fortran !!'
           STOP
         END IF
!         IF (basis_set%ndim < 1 .OR. basis_set%nb < 1 .OR. nq < 1) THEN
!           write(out_unitp,*) ' ERROR in ',name_sub
!           write(out_unitp,*) '  WRONG paramter values: ndim or nb or nq',      &
!                       basis_set%ndim,basis_set%nb,nq
!           write(out_unitp,*) '  CHECK the fortran !!'
!           STOP
!         END IF
         ! first deallocation
         CALL dealloc_dnb_OF_basis(basis_set)

         ! then allocation
         IF (basis_set%ndim > 0) THEN
           CALL alloc_NParray(basis_set%tab_ndim_index,(/basis_set%ndim,basis_set%nb /), &
                         "basis_set%tab_ndim_index",name_sub)
           basis_set%tab_ndim_index(:,:) = 0
         END IF

         !write(out_unitp,*) 'cplx,nb,nq,ndim',basis_set%cplx,basis_set%nb,nq,basis_set%ndim
         CALL flush_perso(6)

         IF (basis_set%cplx) THEN

           CALL alloc_dnCplxMat(basis_set%dnCGB,                        &
                                nq,basis_set%nb,basis_set%ndim,nderiv=nderiv_loc)

           CALL alloc_dnCplxMat(basis_set%dnCBG,                        &
                                basis_set%nb,nq,basis_set%ndim,nderiv=0)

         ELSE

           CALL alloc_dnMat(basis_set%dnRGB,                            &
                                nq,basis_set%nb,basis_set%ndim,nderiv=nderiv_loc)

           CALL alloc_dnMat(basis_set%dnRBG,                            &
                                basis_set%nb,nq,basis_set%ndim,nderiv=0)


         END IF
         IF (debug) write(out_unitp,*) 'END ',name_sub
         CALL flush_perso(out_unitp)

       END SUBROUTINE alloc_dnb_OF_basis
       SUBROUTINE dealloc_dnb_OF_basis(basis_set)
         TYPE(basis), intent(inout) :: basis_set

         !----- for debuging --------------------------------------------------
         integer :: err_mem,memory
         character (len=*), parameter :: name_sub='dealloc_dnb_OF_basis'
         logical,parameter :: debug=.FALSE.
         !logical,parameter :: debug=.TRUE.
         !----- for debuging --------------------------------------------------

         IF (debug) write(out_unitp,*) 'BEGINNING ',name_sub
         CALL flush_perso(out_unitp)

         IF (allocated(basis_set%tab_ndim_index)) THEN
           CALL dealloc_NParray(basis_set%tab_ndim_index,                 &
                             'basis_set%tab_ndim_index',name_sub)
         END IF

         CALL dealloc_dnMat(basis_set%dnRGB)
         CALL dealloc_dnMat(basis_set%dnRBG)
         CALL dealloc_dnMat(basis_set%dnRBGwrho)
         CALL dealloc_dnMat(basis_set%dnRBB)
         CALL dealloc_dnMat(basis_set%dnRGG)

         CALL dealloc_dnCplxMat(basis_set%dnCGB)
         CALL dealloc_dnCplxMat(basis_set%dnCBG)
         CALL dealloc_dnCplxMat(basis_set%dnCBGwrho)
         CALL dealloc_dnCplxMat(basis_set%dnCBB)
         !CALL dealloc_dnCplxMat(basis_set%dnCGG) !not yet in the type

         IF (debug) write(out_unitp,*) 'END ',name_sub
         CALL flush_perso(out_unitp)

       END SUBROUTINE dealloc_dnb_OF_basis
       SUBROUTINE alloc_xw_OF_basis(basis_set)
         TYPE(basis), intent(inout) :: basis_set
         integer :: nq,Lmax
         integer :: err_mem,memory
         character (len=*), parameter :: name_sub='alloc_xw_OF_basis'

         nq = get_nq_FROM_basis(basis_set)

         Lmax = max(0,basis_set%L_SparseBasis)

         IF (basis_set%ndim < 1 .OR. nq < 1) THEN
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) '  WRONG parameter values: ndim or nq',basis_set%ndim,nq
           write(out_unitp,*) '  CHECK the fortran !!'
           STOP
         END IF

         IF (basis_set%xPOGridRep_done) THEN
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) '  xPOGridRep_done = .TRUE. is impossible'
           write(out_unitp,*) '  CHECK the fortran !!'
           STOP
         END IF

         IF (allocated(basis_set%x))  THEN
           CALL dealloc_NParray(basis_set%x,"basis_set%x",name_sub)
         END IF
         CALL alloc_NParray(basis_set%x,(/ basis_set%ndim,nq /),          &
                         "basis_set%x",name_sub)
         basis_set%x     = ZERO

         IF (allocated(basis_set%w))  THEN
           CALL dealloc_NParray(basis_set%w,"basis_set%w",name_sub)
         END IF
         CALL alloc_NParray(basis_set%w,(/ nq /),"basis_set%w",name_sub)
         basis_set%w     = ZERO

         IF (allocated(basis_set%wrho))  THEN
           CALL dealloc_NParray(basis_set%wrho,"basis_set%wrho",name_sub)
         END IF
         CALL alloc_NParray(basis_set%wrho,(/ nq /),"basis_set%wrho",name_sub)
         basis_set%wrho     = ZERO

         IF (allocated(basis_set%rho))  THEN
           CALL dealloc_NParray(basis_set%rho,"basis_set%rho",name_sub)
         END IF
         CALL alloc_NParray(basis_set%rho,(/ nq /),"basis_set%rho",name_sub)
         basis_set%rho     = ZERO

       END SUBROUTINE alloc_xw_OF_basis
       SUBROUTINE dealloc_xw_OF_basis(basis_set)
         TYPE(basis), intent(inout) :: basis_set

         integer :: err_mem,memory
         character (len=*), parameter :: name_sub='dealloc_xw_OF_basis'

         IF (allocated(basis_set%x))  THEN
           CALL dealloc_NParray(basis_set%x,"basis_set%x",name_sub)
         END IF

         IF (allocated(basis_set%w))  THEN
           CALL dealloc_NParray(basis_set%w,"basis_set%w",name_sub)
         END IF

         IF (allocated(basis_set%wrho))  THEN
           CALL dealloc_NParray(basis_set%wrho,"basis_set%wrho",name_sub)
         END IF

         IF (allocated(basis_set%rho))  THEN
           CALL dealloc_NParray(basis_set%rho,"basis_set%rho",name_sub)
         END IF

       END SUBROUTINE dealloc_xw_OF_basis
       FUNCTION get_x_OF_basis(basis_set)
         TYPE(basis), intent(in) :: basis_set
         real (kind=Rkind), allocatable :: get_x_OF_basis(:,:)

         integer :: err_mem,memory
         character (len=*), parameter :: name_sub='get_x_OF_basis'

         IF (allocated(basis_set%x)) get_x_OF_basis = basis_set%x

       END FUNCTION get_x_OF_basis
       FUNCTION get_x_AT_iq_OF_basis(basis_set,iq)
         TYPE(basis), intent(in) :: basis_set
         integer,     intent(in) :: iq
         real (kind=Rkind), allocatable :: get_x_AT_iq_OF_basis(:)

         integer :: err_mem,memory
         character (len=*), parameter :: name_sub='get_x_AT_iq_OF_basis'

         IF (allocated(basis_set%x)) get_x_AT_iq_OF_basis = basis_set%x(:,iq)

       END FUNCTION get_x_AT_iq_OF_basis
       FUNCTION get_rho_OF_basis(basis_set)
         TYPE(basis), intent(in) :: basis_set
         real (kind=Rkind), allocatable :: get_rho_OF_basis(:)

         integer :: err_mem,memory
         character (len=*), parameter :: name_sub='get_rho_OF_basis'

         IF (allocated(basis_set%rho)) get_rho_OF_basis = basis_set%rho

       END FUNCTION get_rho_OF_basis
       FUNCTION get_rho_AT_iq_OF_basis(basis_set,iq)
         TYPE(basis), intent(in) :: basis_set
         integer,     intent(in) :: iq
         real (kind=Rkind) :: get_rho_AT_iq_OF_basis

         integer :: err_mem,memory
         character (len=*), parameter :: name_sub='get_rho_AT_iq_OF_basis'

         IF (allocated(basis_set%rho)) THEN
            get_rho_AT_iq_OF_basis = basis_set%rho(iq)
         ELSE
           STOP 'basis_set%rho(:) is not allocated. STOP in get_rho_AT_iq_OF_basis'
         END IF

       END FUNCTION get_rho_AT_iq_OF_basis
       FUNCTION get_w_OF_basis(basis_set,L)
         TYPE(basis),       intent(in) :: basis_set
         integer, optional, intent(in) :: L
         real (kind=Rkind), allocatable :: get_w_OF_basis(:)

         integer :: err_mem,memory
         character (len=*), parameter :: name_sub='get_w_OF_basis'

         IF (present(L)) THEN
           STOP 'not yet in get_w_OF_basis'
         ELSE
           IF (allocated(basis_set%w)) get_w_OF_basis = basis_set%w
         END IF

       END FUNCTION get_w_OF_basis
       FUNCTION get_w_AT_iq_OF_basis(basis_set,iq,L)
         TYPE(basis),       intent(in) :: basis_set
         integer,           intent(in) :: iq
         integer, optional, intent(in) :: L
         real (kind=Rkind) :: get_w_AT_iq_OF_basis

         integer :: err_mem,memory
         character (len=*), parameter :: name_sub='get_w_AT_iq_OF_basis'

         IF (present(L)) THEN
           STOP 'not yet in get_w_AT_iq_OF_basis'
         ELSE
           IF (allocated(basis_set%w)) THEN
              get_w_AT_iq_OF_basis = basis_set%w(iq)
           ELSE
             STOP 'basis_set%w(:) is not allocated. STOP in get_w_AT_iq_OF_basis'
           END IF
         END IF

       END FUNCTION get_w_AT_iq_OF_basis
       FUNCTION get_wrho_OF_basis(basis_set,L)
         TYPE(basis),       intent(in) :: basis_set
         integer, optional, intent(in) :: L
         real (kind=Rkind), allocatable :: get_wrho_OF_basis(:)

         integer :: err_mem,memory
         character (len=*), parameter :: name_sub='get_wrho_OF_basis'

         IF (present(L)) THEN
           STOP 'not yet in get_w_OF_basis'
         ELSE
           IF (allocated(basis_set%wrho)) THEN
             CALL alloc_NParray(get_wrho_OF_basis,shape(basis_set%wrho),&
                               'get_wrho_OF_basis',name_sub)
             get_wrho_OF_basis = basis_set%wrho
           ELSE
             CALL alloc_NParray(get_wrho_OF_basis,(/ basis_set%nb /),&
                               'get_wrho_OF_basis',name_sub)
             get_wrho_OF_basis = ONE
           END IF
         END IF

       END FUNCTION get_wrho_OF_basis
       FUNCTION get_wrho_AT_iq_OF_basis(basis_set,iq,L)
         TYPE(basis),       intent(in) :: basis_set
         integer,           intent(in) :: iq
         integer, optional, intent(in) :: L
         real (kind=Rkind) :: get_wrho_AT_iq_OF_basis

         integer :: err_mem,memory
         character (len=*), parameter :: name_sub='get_wrho_AT_iq_OF_basis'

         IF (present(L)) THEN
           STOP 'L present: not yet. STOP in get_wrho_AT_iq_OF_basis'
         ELSE
           IF (allocated(basis_set%wrho)) THEN
              get_wrho_AT_iq_OF_basis = basis_set%wrho(iq)
           ELSE
             STOP 'basis_set%wrho(:) is not allocated. STOP in get_wrho_AT_iq_OF_basis'
           END IF
         END IF

       END FUNCTION get_wrho_AT_iq_OF_basis
       SUBROUTINE alloc_tab_Pbasis_OF_basis(basis_set)
         TYPE(basis), intent(inout) :: basis_set

         integer :: i
         integer :: err_mem,memory
         character (len=*), parameter :: name_sub='alloc_tab_Pbasis_OF_basis'

         IF (basis_set%nb_basis < 0) THEN
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) '  WRONG paramter values: nb_basis',basis_set%nb_basis
           write(out_unitp,*) '  CHECK the fortran !!'
           STOP
         END IF

         IF (associated(basis_set%tab_Pbasis))  THEN
           DO i=1,size(basis_set%tab_Pbasis)
             CALL dealloc_basis(basis_set%tab_Pbasis(i)%Pbasis)

             CALL dealloc_array(basis_set%tab_Pbasis(i)%Pbasis,         &
                               'basis_set%tab_Pbasis(i)%Pbasis',name_sub)
           END DO
           CALL dealloc_array(basis_set%tab_Pbasis,                     &
                             'basis_set%tab_Pbasis',name_sub)
         END IF

         CALL alloc_array(basis_set%tab_Pbasis,(/ basis_set%nb_basis /),&
                         'basis_set%tab_Pbasis',name_sub)

         DO i=1,basis_set%nb_basis
           CALL alloc_array(basis_set%tab_Pbasis(i)%Pbasis,             &
                           'basis_set%tab_Pbasis(i)%Pbasis',name_sub)
         END DO

         IF (associated(basis_set%Tab_OF_Tabnb2)) THEN
           DO i=1,size(basis_set%Tab_OF_Tabnb2)
             CALL dealloc_IntVec(basis_set%Tab_OF_Tabnb2(i))
           END DO
           CALL dealloc_array(basis_set%Tab_OF_Tabnb2,                  &
                             'basis_set%Tab_OF_Tabnb2',name_sub)
         END IF
         CALL alloc_array(basis_set%Tab_OF_Tabnb2,(/ basis_set%nb_basis /),&
                         'basis_set%Tab_OF_Tabnb2',name_sub)


       END SUBROUTINE alloc_tab_Pbasis_OF_basis

!================================================================
!       dealloc memory of the TYPE basis
!================================================================
      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
       RECURSIVE SUBROUTINE dealloc_basis(basis_set,Basis_FOR_SG,keep_Rvec)
         USE mod_param_RD
         USE mod_MPI
         IMPLICIT NONE

         TYPE (basis)      :: basis_set
         logical, optional :: Basis_FOR_SG,keep_Rvec

         logical           :: Basis_FOR_SG_loc,keep_Rvec_loc
         integer           :: i,j

         integer :: err_mem,memory
         character (len=*), parameter :: name_sub='dealloc_basis'
         
!         character(14) :: name_subp='dealloc_basis'
!         character(2)  :: name_int
!         character(16) :: name_all
!         i=MPI_id
!         write(name_int, '(I2)') i
!         name_all=name_subp//name_int
         
         Basis_FOR_SG_loc = .FALSE.
         IF (present(Basis_FOR_SG)) Basis_FOR_SG_loc = Basis_FOR_SG

         keep_Rvec_loc = .FALSE.
         IF (present(keep_Rvec)) keep_Rvec_loc = keep_Rvec

         basis_set%active                = .FALSE.
         basis_set%print_info_OF_basisDP = .TRUE.

         IF (.NOT. keep_Rvec_loc) THEN
           basis_set%ndim         = 0
           IF (allocated(basis_set%iQdyn))  THEN
             CALL dealloc_NParray(basis_set%iQdyn,"basis_set%iQdyn",name_sub)
           END IF
           IF (allocated(basis_set%Tabder_Qdyn_TO_Qbasis) )  THEN
             CALL dealloc_NParray(basis_set%Tabder_Qdyn_TO_Qbasis,       &
                               "basis_set%Tabder_Qdyn_TO_Qbasis",name_sub)
           END IF
         END IF

         IF (allocated(basis_set%nrho) )  THEN
           CALL dealloc_NParray(basis_set%nrho,"basis_set%nrho",name_sub)
         END IF


         IF (.NOT. Basis_FOR_SG_loc) THEN
           CALL dealloc_SymAbelian(basis_set%P_SymAbelian)
         ELSE
           nullify(basis_set%P_SymAbelian)
         END IF

         basis_set%check_basis            = .TRUE.
         basis_set%check_nq_OF_basis      = .TRUE.
         basis_set%packed                 = .FALSE.
         basis_set%packed_done            = .FALSE.
         basis_set%primitive              = .FALSE.
         basis_set%primitive_done         = .FALSE.

         IF (.NOT. keep_Rvec_loc) THEN
           basis_set%auto_basis             = .FALSE.
           basis_set%contrac                = .FALSE.
           basis_set%auto_contrac           = .FALSE.
           basis_set%contrac_analysis       = .FALSE.
           basis_set%max_ene_contrac        = ONETENTH
           basis_set%POGridRep              = .FALSE.
           basis_set%POGridRep_polyortho    = .FALSE.
           basis_set%xPOGridRep_done        = .FALSE.
           basis_set%nbc                    = 0
           basis_set%nqc                    = 0
           basis_set%max_nbc                = 0
           basis_set%min_nbc                = 0
           basis_set%auto_contrac_type1_TO  = 100
           basis_set%auto_contrac_type21_TO = 100
           basis_set%nqPLUSnbc_TO_nqc       = 0
           basis_set%read_contrac_file      = .FALSE.
           IF (allocated(basis_set%Rvec))  THEN
             CALL dealloc_NParray(basis_set%Rvec,"basis_set%Rvec",name_sub)
           END IF
           CALL file_dealloc(basis_set%file_contrac)
         END IF

         basis_set%type  = 0
         basis_set%name  = "0"

         CALL dealloc_dnb_OF_basis(basis_set)

         CALL dealloc_dnMat(basis_set%dnRGG)
         basis_set%dnGGRep      = .FALSE.
         basis_set%dnGGRep_done = .FALSE.
         IF (.NOT. keep_Rvec_loc) THEN
           basis_set%nb      = 0
           basis_set%nb_init = 0
         END IF

         CALL Set_nq_OF_basis(basis_set,nq=0)
         CALL Basis_Grid_ParamTOBasis_Grid_Param_init(basis_set%Basis_Grid_Para)

         IF (allocated(basis_set%EneH0))    THEN
           CALL dealloc_NParray(basis_set%EneH0,                        &
                               "basis_set%EneH0",name_sub)
         END IF

         basis_set%dnBBRep      = .FALSE.
         basis_set%dnBBRep_done = .FALSE.
         CALL dealloc_dnMat(basis_set%dnRBB)
         CALL dealloc_dnCplxMat(basis_set%dnCBB)

         basis_set%cplx  = .FALSE.

         basis_set%Nested        = 0
         basis_set%nq_max_Nested = -1

         CALL dealloc_xw_OF_basis(basis_set)

         IF (allocated(basis_set%cte_Transfo))         THEN
           CALL dealloc_NParray(basis_set%cte_Transfo,"basis_set%cte_Transfo",name_sub)
         END IF

         basis_set%nb_Transfo = 0
         IF (allocated(basis_set%A))         THEN
           CALL dealloc_NParray(basis_set%A,"basis_set%A",name_sub)
         END IF
         IF (allocated(basis_set%B))         THEN
           CALL dealloc_NParray(basis_set%B,"basis_set%B",name_sub)
         END IF
         IF (allocated(basis_set%Q0))        THEN
           CALL dealloc_NParray(basis_set%Q0,"basis_set%Q0",name_sub)
         END IF
         IF (allocated(basis_set%scaleQ))    THEN
           CALL dealloc_NParray(basis_set%scaleQ,"basis_set%scaleQ",name_sub)
         END IF

         IF (allocated(basis_set%opt_A) ) THEN
           CALL dealloc_NParray(basis_set%opt_A,"basis_set%opt_A",name_sub)
         END IF
         IF (allocated(basis_set%opt_B) ) THEN
           CALL dealloc_NParray(basis_set%opt_B,"basis_set%opt_B",name_sub)
         END IF
         IF (allocated(basis_set%opt_Q0) ) THEN
           CALL dealloc_NParray(basis_set%opt_Q0,"basis_set%opt_Q0",name_sub)
         END IF
         IF (allocated(basis_set%opt_scaleQ) ) THEN
           CALL dealloc_NParray(basis_set%opt_scaleQ,"basis_set%opt_scaleQ",name_sub)
         END IF


         CALL dealloc_Basis_L_TO_n(basis_set%L_TO_nb)
         CALL dealloc_Basis_L_TO_n(basis_set%L_TO_nq)

         basis_set%L_SparseGrid             = -1
         basis_set%L_SparseBasis            = -1
         basis_set%With_L                   = .FALSE.

         basis_set%SparseGrid_type          = 0

         basis_set%SparseGrid_With_Cuba     = .TRUE.
         basis_set%SparseGrid_With_Smolyak  = .TRUE.
         basis_set%SparseGrid_With_DP       = .TRUE.
         basis_set%nqSG_SMALLER_nqDP        = .TRUE.
         basis_set%nb_basis                 = 0
         basis_set%nb_SG                    = 0

         IF (associated(basis_set%Tab_OF_Tabnb2)) THEN
           DO i=1,size(basis_set%Tab_OF_Tabnb2)
             CALL dealloc_IntVec(basis_set%Tab_OF_Tabnb2(i))
           END DO
           CALL dealloc_array(basis_set%Tab_OF_Tabnb2,                  &
                             'basis_set%Tab_OF_Tabnb2',name_sub)
         END IF

         IF (allocated(basis_set%WeightSG))    THEN
           CALL dealloc_NParray(basis_set%WeightSG,                       &
                             "basis_set%WeightSG",name_sub)
         END IF

         IF (associated(basis_set%tab_PbasisSG))  THEN
           DO i=1,size(basis_set%tab_PbasisSG)
             CALL dealloc_basis(basis_set%tab_PbasisSG(i)%Pbasis,Basis_FOR_SG=.TRUE.)

             CALL dealloc_array(basis_set%tab_PbasisSG(i)%Pbasis,       &
                               'basis_set%tab_PbasisSG(i)%Pbasis',name_sub)
           END DO
           CALL dealloc_array(basis_set%tab_PbasisSG,                   &
                             'basis_set%tab_PbasisSG',name_sub)
         END IF
         IF (associated(basis_set%tab_basisPrimSG))  THEN
           DO i=lbound(basis_set%tab_basisPrimSG,dim=1),ubound(basis_set%tab_basisPrimSG,dim=1)
           DO j=lbound(basis_set%tab_basisPrimSG,dim=2),ubound(basis_set%tab_basisPrimSG,dim=2)
             CALL dealloc_basis(basis_set%tab_basisPrimSG(i,j))
           END DO
           END DO
           CALL dealloc_array(basis_set%tab_basisPrimSG,                &
                             'basis_set%tab_basisPrimSG',name_sub)
         END IF


         IF (keep_Rvec_loc) THEN
           IF (associated(basis_set%tab_Pbasis))  THEN
             DO i=1,size(basis_set%tab_Pbasis)
               IF (.NOT. basis_set%tab_basis_linked) THEN
                 CALL dealloc_basis(basis_set%tab_Pbasis(i)%Pbasis,keep_Rvec=.TRUE.)
               END IF
             END DO
           END IF
         ELSE
           IF (associated(basis_set%tab_Pbasis))  THEN
             DO i=1,size(basis_set%tab_Pbasis)
               IF (.NOT. basis_set%tab_basis_linked .AND.               &
                       associated(basis_set%tab_Pbasis(i)%Pbasis) ) THEN
                 CALL dealloc_basis(basis_set%tab_Pbasis(i)%Pbasis)
                 CALL dealloc_array(basis_set%tab_Pbasis(i)%Pbasis,     &
                                   'basis_set%tab_Pbasis(i)%Pbasis',name_sub)
               END IF
             END DO
             CALL dealloc_array(basis_set%tab_Pbasis,                   &
                               'basis_set%tab_Pbasis',name_sub)
             basis_set%tab_basis_linked = .FALSE.
           END IF
           basis_set%tab_basis_done = .FALSE.

         END IF
         CALL dealloc_nDindex(basis_set%nDindG)

         IF (.NOT. keep_Rvec_loc) THEN
           IF (associated(basis_set%nDindB) .AND. .NOT. Basis_FOR_SG_loc) THEN
             CALL dealloc_nDindex(basis_set%nDindB)
             CALL dealloc_array(basis_set%nDindB,                       &
                               'basis_set%nDindB',name_sub)
           END IF
           basis_set%Type_OF_nDindB        = 1
           basis_set%Norm_OF_nDindB        = huge(1)
           basis_set%weight_OF_nDindB      = ONE
           basis_set%MaxCoupling_OF_nDindB = -1

           IF (associated(basis_set%nDindB_uncontracted) .AND. .NOT. Basis_FOR_SG_loc) THEN
             CALL dealloc_nDindex(basis_set%nDindB_uncontracted)
             CALL dealloc_array(basis_set%nDindB_uncontracted,          &
                               'basis_set%nDindB_uncontracted',name_sub)
           END IF
         END IF

         CALL dealloc_SGType2(basis_set%para_SGType2)
         CALL dealloc_RotBasis_Param(basis_set%RotBasis)

       END SUBROUTINE dealloc_basis

      SUBROUTINE alloc_array_OF_Basisdim0(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (basis), pointer, intent(inout) :: tab

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=0
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_Basisdim0'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (associated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       memory = 1 ! true pointer
       allocate(tab,stat=err_mem) ! in alloc_array
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'basis')

      END SUBROUTINE alloc_array_OF_Basisdim0
      SUBROUTINE dealloc_array_OF_Basisdim0(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (basis), pointer, intent(inout) :: tab
      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_Basisdim0'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN
       IF (.NOT. associated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = 1
       deallocate(tab,stat=err_mem) ! in alloc_array
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'basis')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_Basisdim0

      SUBROUTINE alloc_array_OF_Basisdim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      TYPE (basis), pointer, intent(inout) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=1
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_Basisdim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (associated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)

         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1)),stat=err_mem) ! in alloc_array
       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1)),stat=err_mem) ! in alloc_array
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'basis')

      END SUBROUTINE alloc_array_OF_Basisdim1
      SUBROUTINE dealloc_array_OF_Basisdim1(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (basis), pointer, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_Basisdim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN
       IF (.NOT. associated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab)
       deallocate(tab,stat=err_mem) ! in alloc_array
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'basis')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_Basisdim1
      SUBROUTINE alloc_array_OF_Basisdim2(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      TYPE (basis), pointer, intent(inout) :: tab(:,:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub


      integer, parameter :: ndim=2
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_Basisdim2'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       IF (associated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)

         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1),                              & ! in alloc_array
                      tab_lb(2):tab_ub(2)),stat=err_mem)
       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1),tab_ub(2)),stat=err_mem)  ! in alloc_array
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'basis')

      END SUBROUTINE alloc_array_OF_Basisdim2
      SUBROUTINE dealloc_array_OF_Basisdim2(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (basis), pointer, intent(inout) :: tab(:,:)
      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_Basisdim2'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN

       IF (.NOT. associated(tab))                                       &
                 CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab)
       deallocate(tab,stat=err_mem)  ! in alloc_array
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'basis')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_Basisdim2


      SUBROUTINE alloc_array_OF_P_Basisdim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      TYPE (P_basis), pointer, intent(inout) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=1
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_P_Basisdim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (associated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)

         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1)),stat=err_mem) ! in alloc_array
       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1)),stat=err_mem) ! in alloc_array
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'P_basis')

      END SUBROUTINE alloc_array_OF_P_Basisdim1
      SUBROUTINE dealloc_array_OF_P_Basisdim1(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (P_basis), pointer, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_P_Basisdim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN
       IF (.NOT. associated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab)
       deallocate(tab,stat=err_mem) ! in alloc_array
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'P_basis')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_P_Basisdim1

      SUBROUTINE alloc_array_OF_P_Basisdim2(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      TYPE (P_basis), pointer, intent(inout) :: tab(:,:)
      integer,                 intent(in)    :: tab_ub(:)
      integer, optional,       intent(in)    :: tab_lb(:)
      character (len=*),       intent(in)    :: name_var,name_sub


      integer, parameter :: ndim=2
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_P_Basisdim2'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       IF (associated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)

         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1),                              & ! in alloc_array
                      tab_lb(2):tab_ub(2)),stat=err_mem)
       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1),tab_ub(2)),stat=err_mem) ! in alloc_array
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'P_basis')

      END SUBROUTINE alloc_array_OF_P_Basisdim2
      SUBROUTINE dealloc_array_OF_P_Basisdim2(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (P_basis), pointer, intent(inout) :: tab(:,:)
      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_P_Basisdim2'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN

       IF (.NOT. associated(tab))                                       &
                 CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab)
       deallocate(tab,stat=err_mem) ! in alloc_array
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'P_basis')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_P_Basisdim2

!================================================================
!     clean basis when the basis is pack and for dnBBRep
!================================================================
      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
       SUBROUTINE clean_basis(basis_set)
         IMPLICIT NONE

         TYPE (basis) :: basis_set
         integer      :: i,j

         integer :: err_mem,memory
         character (len=*), parameter :: name_sub='clean_basis'

         IF (.NOT. basis_set%packed_done) RETURN

         IF (basis_set%dnBBRep_done) THEN

           write(out_unitp,*) 'clean d1b, d2b, d1cb, d2cb'
           IF (associated(basis_set%dnRGB%d1)) THEN
             CALL dealloc_array(basis_set%dnRGB%d1,'basis_set%dnRGB%d1',name_sub)
             basis_set%dnRGB%nderiv = 0
           END IF

           IF (associated(basis_set%dnRGB%d2)) THEN
             CALL dealloc_array(basis_set%dnRGB%d2,'basis_set%dnRGB%d2',name_sub)
             basis_set%dnRGB%nderiv = 0
           END IF

           IF (associated(basis_set%dnCGB%d1)) THEN
             CALL dealloc_array(basis_set%dnCGB%d1,'basis_set%dnCGB%d1',name_sub)
             basis_set%dnCGB%nderiv = 0
           END IF

           IF (associated(basis_set%dnCGB%d2)) THEN
             CALL dealloc_array(basis_set%dnCGB%d2,'basis_set%dnCGB%d2',name_sub)
             basis_set%dnCGB%nderiv = 0
           END IF

         END IF
         IF (print_level > -1) write(out_unitp,*) 'clean all the tab_basis'

         IF (associated(basis_set%Tab_OF_Tabnb2)) THEN
           DO i=1,size(basis_set%Tab_OF_Tabnb2)
             CALL dealloc_IntVec(basis_set%Tab_OF_Tabnb2(i))
           END DO
           CALL dealloc_array(basis_set%Tab_OF_Tabnb2,                  &
                             'basis_set%Tab_OF_Tabnb2',name_sub)
         END IF

         IF (allocated(basis_set%WeightSG))  THEN
           CALL dealloc_NParray(basis_set%WeightSG,"basis_set%WeightSG",name_sub)
         END IF

         IF (associated(basis_set%tab_PbasisSG))  THEN
           DO i=1,size(basis_set%tab_PbasisSG)
             CALL dealloc_basis(basis_set%tab_PbasisSG(i)%Pbasis)
             CALL dealloc_array(basis_set%tab_PbasisSG(i)%Pbasis,       &
                               'basis_set%tab_PbasisSG(i)%Pbasis',name_sub)
           END DO
           CALL dealloc_array(basis_set%tab_PbasisSG,                   &
                             'basis_set%tab_PbasisSG',name_sub)
         END IF
         IF (associated(basis_set%tab_basisPrimSG))  THEN
           DO i=lbound(basis_set%tab_basisPrimSG,dim=1),ubound(basis_set%tab_basisPrimSG,dim=1)
           DO j=lbound(basis_set%tab_basisPrimSG,dim=2),ubound(basis_set%tab_basisPrimSG,dim=2)
             CALL dealloc_basis(basis_set%tab_basisPrimSG(i,j))
           END DO
           END DO
           CALL dealloc_array(basis_set%tab_basisPrimSG,                &
                             'basis_set%tab_basisPrimSG',name_sub)
         END IF


       END SUBROUTINE clean_basis

!     ===========================================================
!      basis2=basis1
!     ===========================================================
      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      RECURSIVE SUBROUTINE basis2TObasis1(basis_set1,basis_set2,init_only,with_SG,Basis_FOR_SG)
        USE mod_param_RD
        IMPLICIT NONE
        TYPE (basis), intent(in)    :: basis_set2
        TYPE (basis), intent(inout) :: basis_set1
        logical, optional           :: init_only,with_SG,Basis_FOR_SG


        logical                     :: init_only_loc,with_SG_loc,Basis_FOR_SG_loc
        integer                     :: i,j,li,ui,lj,uj,iopt
        integer                     :: n1,n2,n3,n4

        integer :: err_mem,memory
         character (len=*), parameter :: name_sub='basis2TObasis1'

        !CALL RecWrite_basis(basis_set2,write_all=.TRUE.)

        CALL dealloc_basis(basis_set1)

        init_only_loc = .FALSE.
        IF (present(init_only)) init_only_loc = init_only

        with_SG_loc = .FALSE.
        IF (present(with_SG)) with_SG_loc = with_SG

        Basis_FOR_SG_loc = .FALSE.
        IF (present(Basis_FOR_SG)) Basis_FOR_SG_loc = Basis_FOR_SG

        basis_set1%active                = basis_set2%active
        basis_set1%print_info_OF_basisDP = basis_set2%print_info_OF_basisDP

        basis_set1%ndim              = basis_set2%ndim

        basis_set1%nb_Transfo        = basis_set2%nb_Transfo

        CALL alloc_init_basis(basis_set1)

        IF (init_only_loc) THEN
          basis_set1%nb              = basis_set2%nb_init
          basis_set1%Basis_Grid_Para = basis_set2%Basis_Grid_Para
          CALL Basis_Grid_Param_initTOBasis_Grid_Param(basis_set1%Basis_Grid_Para)
        ELSE
          basis_set1%nb              = basis_set2%nb
          basis_set1%Basis_Grid_Para = basis_set2%Basis_Grid_Para

        END IF
        basis_set1%nb_init           = basis_set2%nb_init

        basis_set1%nq_max_Nested     = basis_set2%nq_max_Nested
        basis_set1%Nested            = basis_set2%Nested
        basis_set1%cplx              = basis_set2%cplx

        basis_set1%check_basis       = basis_set2%check_basis
        basis_set1%check_nq_OF_basis = basis_set2%check_nq_OF_basis
        basis_set1%packed            = basis_set2%packed
        basis_set1%primitive         = basis_set2%primitive
        IF (init_only_loc) THEN
          basis_set1%primitive_done       = .FALSE.
          basis_set1%packed_done          = .FALSE.
        ELSE
          basis_set1%primitive_done       = basis_set2%primitive_done
          basis_set1%packed_done          = basis_set2%packed_done
        END IF

        basis_set1%auto_basis             = basis_set2%auto_basis

        basis_set1%contrac                = basis_set2%contrac
        basis_set1%auto_contrac           = basis_set2%auto_contrac
        basis_set1%contrac_analysis       = basis_set2%contrac_analysis
        basis_set1%max_ene_contrac        = basis_set2%max_ene_contrac
        basis_set1%make_cubature          = basis_set2%make_cubature
        basis_set1%Restart_make_cubature  = basis_set2%Restart_make_cubature

        basis_set1%max_nbc                = basis_set2%max_nbc
        basis_set1%min_nbc                = basis_set2%min_nbc
        basis_set1%auto_contrac_type1_TO  = basis_set2%auto_contrac_type1_TO
        basis_set1%auto_contrac_type21_TO = basis_set2%auto_contrac_type21_TO

        basis_set1%POGridRep             = basis_set2%POGridRep
        basis_set1%POGridRep_polyortho   = basis_set2%POGridRep_polyortho
        basis_set1%xPOGridRep_done       = basis_set2%xPOGridRep_done
        basis_set1%nbc               = basis_set2%nbc
        basis_set1%nqc               = basis_set2%nqc
        basis_set1%read_contrac_file = basis_set2%read_contrac_file
        basis_set1%file_contrac      = basis_set2%file_contrac

        basis_set1%type              = basis_set2%type
        basis_set1%name              = basis_set2%name

        basis_set1%nb_basis          = basis_set2%nb_basis

        IF (allocated(basis_set2%iQdyn)) THEN
           basis_set1%iQdyn = basis_set2%iQdyn
        END IF

        IF (allocated(basis_set2%nrho) ) THEN
           basis_set1%nrho = basis_set2%nrho
        END IF

        IF (.NOT. init_only_loc .AND. .NOT. Basis_FOR_SG_loc) THEN
          CALL SymAbelian1_TO_SymAbelian2(basis_set2%P_SymAbelian,      &
                                                basis_set1%P_SymAbelian)
        END IF

        IF (allocated(basis_set2%tab_ndim_index) .AND. .NOT. init_only_loc) THEN
          CALL alloc_NParray(basis_set1%tab_ndim_index,                   &
                                  (/ basis_set2%ndim,basis_set2%nb /),  &
                          "basis_set1%tab_ndim_index",name_sub)
          basis_set1%tab_ndim_index = basis_set2%tab_ndim_index
        END IF

        basis_set1%opt_param = basis_set2%opt_param
        IF (allocated(basis_set2%cte_Transfo)) THEN
          IF (allocated(basis_set1%cte_Transfo)) THEN
            CALL dealloc_NParray(basis_set1%cte_Transfo,                  &
                                "basis_set1%cte_Transfo",name_sub)
          END IF
          CALL alloc_NParray(basis_set1%cte_Transfo,                      &
                                        shape(basis_set2%cte_Transfo),  &
                            "basis_set1%cte_Transfo",name_sub)

          basis_set1%cte_Transfo = basis_set2%cte_Transfo
        END IF
        IF (allocated(basis_set2%A)) THEN
          basis_set1%A      = basis_set2%A
        END IF
        IF (allocated(basis_set2%B)) THEN
          basis_set1%B      = basis_set2%B
        END IF
        IF (allocated(basis_set2%Q0)) THEN
          basis_set1%Q0      = basis_set2%Q0
        END IF
        IF (allocated(basis_set2%scaleQ)) THEN
          basis_set1%scaleQ  = basis_set2%scaleQ
        END IF
        IF (allocated(basis_set2%opt_A)) THEN
          basis_set1%opt_A      = basis_set2%opt_A
        END IF
        IF (allocated(basis_set2%opt_B)) THEN
          basis_set1%opt_B      = basis_set2%opt_B
        END IF
        IF (allocated(basis_set2%opt_Q0)) THEN
          basis_set1%opt_Q0      = basis_set2%opt_Q0
        END IF
        IF (allocated(basis_set2%opt_scaleQ)) THEN
          basis_set1%opt_scaleQ  = basis_set2%opt_scaleQ
        END IF


        IF (allocated(basis_set2%Tabder_Qdyn_TO_Qbasis)) THEN
          n1 = ubound(basis_set2%Tabder_Qdyn_TO_Qbasis,dim=1)
          CALL alloc_NParray(basis_set1%Tabder_Qdyn_TO_Qbasis,(/ n1 /),   &
                          "basis_set1%Tabder_Qdyn_TO_Qbasis",name_sub, (/ 0 /))
          basis_set1%Tabder_Qdyn_TO_Qbasis = basis_set2%Tabder_Qdyn_TO_Qbasis
        END IF

        IF (basis_set2%nb_basis > 0) THEN
          IF (init_only_loc) THEN
            basis_set1%tab_basis_done   = .FALSE.
          ELSE
            basis_set1%tab_basis_done   = basis_set2%tab_basis_done
          END IF
          basis_set1%tab_basis_linked = basis_set2%tab_basis_linked
          IF (basis_set2%tab_basis_linked) THEN
            IF (associated(basis_set2%tab_Pbasis) ) THEN
              CALL alloc_array(basis_set1%tab_Pbasis,(/ basis_set2%nb_basis /),&
                              'basis_set1%tab_Pbasis',name_sub)
              DO i=1,basis_set2%nb_basis
                basis_set1%tab_Pbasis(i)%Pbasis => basis_set2%tab_Pbasis(i)%Pbasis
              END DO
            END IF
            !STOP 'pb with tab_basis_linked'
          ELSE
            IF (associated(basis_set2%tab_Pbasis) ) THEN
              CALL alloc_array(basis_set1%tab_Pbasis,(/ basis_set2%nb_basis /),&
                              'basis_set1%tab_Pbasis',name_sub)

              DO i=1,basis_set2%nb_basis
                CALL alloc_array(basis_set1%tab_Pbasis(i)%Pbasis,       &
                                'basis_set1%tab_Pbasis(i)%Pbasis',name_sub)
                CALL basis2TObasis1(basis_set1%tab_Pbasis(i)%Pbasis,    &
                                    basis_set2%tab_Pbasis(i)%Pbasis,    &
                                    init_only_loc)
              END DO
            END IF
          END IF
          IF (associated(basis_set2%Tab_OF_Tabnb2)) THEN
            CALL alloc_array(basis_set1%Tab_OF_Tabnb2,                  &
                                              (/ basis_set2%nb_basis /),&
                            'basis_set1%Tab_OF_Tabnb2',name_sub)
            DO i=1,basis_set2%nb_basis
              CALL sub_IntVec1_TO_IntVec2(basis_set2%Tab_OF_Tabnb2(i),  &
                                          basis_set1%Tab_OF_Tabnb2(i))
            END DO
          END IF
        END IF
        basis_set1%nDindG = basis_set2%nDindG

        IF (.NOT. Basis_FOR_SG_loc) THEN
          IF (.NOT. associated(basis_set2%nDindB)) THEN
            CALL RecWrite_basis(basis_set2)
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' nDindB of basis_set2 is not associated'
            STOP
          END IF
          basis_set1%nDindB = basis_set2%nDindB
          IF (associated(basis_set2%nDindB_uncontracted)) THEN
            CALL alloc_array(basis_set1%nDindB_uncontracted,            &
                            'basis_set1%nDindB_uncontracted',name_sub)
            basis_set1%nDindB_uncontracted = basis_set2%nDindB_uncontracted
          END IF
        END IF

        basis_set1%Type_OF_nDindB          = basis_set2%Type_OF_nDindB
        basis_set1%Norm_OF_nDindB          = basis_set2%Norm_OF_nDindB
        basis_set1%weight_OF_nDindB        = basis_set2%weight_OF_nDindB
        basis_set1%nDinit_OF_nDindB        = basis_set2%nDinit_OF_nDindB
        basis_set1%MaxCoupling_OF_nDindB   = basis_set2%MaxCoupling_OF_nDindB
        basis_set1%nb_OF_MinNorm_OF_nDindB = basis_set2%nb_OF_MinNorm_OF_nDindB
        basis_set1%Div_nb_TO_Norm_OF_nDindB= basis_set2%Div_nb_TO_Norm_OF_nDindB
        basis_set1%contrac_WITH_nDindB     = basis_set2%contrac_WITH_nDindB

        basis_set1%dnBBRep                 = basis_set2%dnBBRep
        basis_set1%dnGGRep                 = basis_set2%dnGGRep

        IF (.NOT. init_only_loc) THEN

          IF (allocated(basis_set2%Rvec)) THEN
            CALL alloc_NParray(basis_set1%Rvec,shape(basis_set2%Rvec),    &
                            "basis_set1%Rvec",name_sub)
            basis_set1%Rvec  = basis_set2%Rvec
          END IF

          IF (basis_set2%dnRGB%alloc) basis_set1%dnRGB = basis_set2%dnRGB
          IF (basis_set2%dnCGB%alloc) basis_set1%dnCGB = basis_set2%dnCGB

          IF (basis_set2%dnRBG%alloc) basis_set1%dnRBG = basis_set2%dnRBG
          IF (basis_set2%dnCBG%alloc) basis_set1%dnCBG = basis_set2%dnCBG


          basis_set1%dnBBRep_done = basis_set2%dnBBRep_done
          IF (basis_set2%dnCBB%alloc) basis_set1%dnCBB = basis_set2%dnCBB
          IF (basis_set2%dnRBB%alloc) basis_set1%dnRBB = basis_set1%dnRBB

          basis_set1%dnGGRep_done = basis_set2%dnGGRep_done
          IF (basis_set2%dnRGG%alloc) basis_set1%dnRGG = basis_set2%dnRGG

          IF (basis_set2%dnRBGwrho%alloc) basis_set1%dnRBGwrho = basis_set2%dnRBGwrho
          IF (basis_set2%dnCBGwrho%alloc) basis_set1%dnCBGwrho = basis_set2%dnCBGwrho

         IF (allocated(basis_set1%EneH0))    THEN
           CALL dealloc_NParray(basis_set1%EneH0,                       &
                               "basis_set1%EneH0",name_sub)
         END IF
         IF (allocated(basis_set2%EneH0))    THEN
           CALL alloc_NParray(basis_set1%EneH0,shape(basis_set2%EneH0), &
                             "basis_se1t%EneH0",name_sub)
           basis_set1%EneH0(:) = basis_set2%EneH0(:)
         END IF



          IF (allocated(basis_set2%x)) THEN
            CALL alloc_NParray(basis_set1%x,shape(basis_set2%x),          &
                            "basis_set1%x",name_sub)
            basis_set1%x      = basis_set2%x
          END IF
          IF (allocated(basis_set2%w)) THEN
            CALL alloc_NParray(basis_set1%w,shape(basis_set2%w),          &
                            "basis_set1%w",name_sub)
            basis_set1%w      = basis_set2%w
          END IF
          IF (allocated(basis_set2%wrho)) THEN
            CALL alloc_NParray(basis_set1%wrho,shape(basis_set2%wrho),    &
                            "basis_set1%wrho",name_sub)
            basis_set1%wrho      = basis_set2%wrho
          END IF
          IF (allocated(basis_set2%rho)) THEN
            CALL alloc_NParray(basis_set1%rho,shape(basis_set2%rho),      &
                            "basis_set1%rho",name_sub)
            basis_set1%rho      = basis_set2%rho
          END IF

        END IF

        ! for sparse basis and grid -----------------------
        basis_set1%L_TO_nb     = basis_set2%L_TO_nb
        basis_set1%L_TO_nq     = basis_set2%L_TO_nq

        basis_set1%With_L                   = basis_set2%With_L
        basis_set1%L_SparseBasis            = basis_set2%L_SparseBasis
        basis_set1%L_SparseGrid             = basis_set2%L_SparseGrid
        basis_set1%SparseGrid_type          = basis_set2%SparseGrid_type
        basis_set1%SparseGrid_With_Cuba     = basis_set2%SparseGrid_With_Cuba
        basis_set1%SparseGrid_With_Smolyak  = basis_set2%SparseGrid_With_Smolyak
        basis_set1%SparseGrid_With_DP       = basis_set2%SparseGrid_With_DP
        basis_set1%nqSG_SMALLER_nqDP        = basis_set2%nqSG_SMALLER_nqDP
        basis_set1%nb_SG                    = basis_set2%nb_SG


        IF (init_only_loc) basis_set1%nb_SG = 0
        IF ( basis_set2%nb_SG > 0 .AND. .NOT. init_only_loc) THEN
          IF (allocated(basis_set2%WeightSG)) THEN
            CALL alloc_NParray(basis_set1%WeightSG,shape(basis_set2%WeightSG), &
                            "basis_set1%WeightSG",name_sub)
            basis_set1%WeightSG = basis_set2%WeightSG
          END IF

          IF (with_SG_loc .AND. associated(basis_set2%tab_PbasisSG)) THEN
            CALL alloc_array(basis_set1%tab_PbasisSG,(/ basis_set2%nb_SG /), &
                            'basis_set1%tab_PbasisSG',name_sub)
            DO i=1,basis_set2%nb_SG
              CALL alloc_array(basis_set1%tab_PbasisSG(i)%Pbasis,       &
                              'basis_set1%tab_PbasisSG(i)%Pbasis',name_sub)

              CALL basis2TObasis1(basis_set1%tab_PbasisSG(i)%Pbasis,         &
                                  basis_set2%tab_PbasisSG(i)%Pbasis,Basis_FOR_SG=.TRUE.)

              basis_set1%tab_PbasisSG(i)%Pbasis%nDindB     => basis_set1%nDindB
              IF (associated(basis_set1%nDindB_uncontracted)) THEN
                basis_set1%tab_PbasisSG(i)%Pbasis%nDindB_uncontracted  => basis_set1%nDindB_uncontracted
              END IF
              basis_set1%tab_PbasisSG(i)%Pbasis%P_SymAbelian => basis_set1%P_SymAbelian
            END DO
          END IF

          IF (associated(basis_set2%tab_basisPrimSG)) THEN
            li=lbound(basis_set2%tab_basisPrimSG,dim=1)
            ui=ubound(basis_set2%tab_basisPrimSG,dim=1)
            lj=lbound(basis_set2%tab_basisPrimSG,dim=2)
            uj=ubound(basis_set2%tab_basisPrimSG,dim=2)
            CALL alloc_array(basis_set1%tab_basisPrimSG,(/ ui,uj /),    &
                            'basis_set1%tab_basisPrimSG',name_sub, (/ li,lj /))

            DO i=li,uj
            DO j=lj,uj
              CALL basis2TObasis1(basis_set1%tab_basisPrimSG(i,j),           &
                                  basis_set2%tab_basisPrimSG(i,j))
            END DO
            END DO
         END IF
        END IF

        basis_set1%para_SGType2 = basis_set2%para_SGType2

        basis_set1%RotBasis     = basis_set2%RotBasis

      END SUBROUTINE basis2TObasis1

      SUBROUTINE Get_MatdnRGG(basis_set,MatRGG,dnba_ind)
      USE mod_system
      IMPLICIT NONE

      TYPE (basis),intent(in)        :: basis_set
      real(kind=Rkind),intent(inout) :: MatRGG(:,:)
      integer, intent(in)            :: dnba_ind(2)

      integer :: nq
!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='Get_MatdnRGG'
!-----------------------------------------------------------
      IF (.NOT. basis_set%packed) RETURN
      nq = get_nq_FROM_basis(basis_set)

      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'ndim',basis_set%ndim
        write(out_unitp,*) 'dnba_ind',dnba_ind
        write(out_unitp,*) 'nq',nq
        write(out_unitp,*) 'shape(MatRGG)',shape(MatRGG)
        write(out_unitp,*) 'alloc dnRGG',basis_set%dnRGG%alloc
        CALL write_dnSVM(basis_set%dnRGG)
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

      IF (.NOT. basis_set%dnRGG%alloc) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'basis_set%dnRGG is not allocated!!'
        write(out_unitp,*) 'CHECK the fortran source'
        STOP
      END IF

      IF (dnba_ind(1) == 0 .AND. dnba_ind(2) == 0) THEN
        CALL mat_id(MatRGG(:,:),nq,nq)
      ELSE IF (dnba_ind(1) == 0) THEN ! first derivative
        MatRGG(:,:) = basis_set%dnRGG%d1(:,:,dnba_ind(2))
      ELSE IF (dnba_ind(2) == 0) THEN ! first derivative
        MatRGG(:,:) = basis_set%dnRGG%d1(:,:,dnba_ind(1))
      ELSE ! 2d derivative
        MatRGG(:,:) = basis_set%dnRGG%d2(:,:,dnba_ind(1),dnba_ind(2))
      END IF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'MatRGG',dnba_ind
        CALL write_Mat(MatRGG,out_unitp,5)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------
      END SUBROUTINE Get_MatdnRGG

      SUBROUTINE Get_MatdnRGB(basis_set,RMatdnb,dnba_ind)
      USE mod_system
      IMPLICIT NONE

      TYPE (basis),intent(in)        :: basis_set
      real(kind=Rkind),intent(inout) :: RMatdnb(:,:)
      integer, intent(in)            :: dnba_ind(2)


      integer :: nq
!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='Get_MatdnRGB'
!-----------------------------------------------------------
      IF (.NOT. basis_set%packed) RETURN
      nq = get_nq_FROM_basis(basis_set)
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nq',nq
      END IF
!-----------------------------------------------------------

    IF (NewBasisEl .AND. basis_set%ndim == 0) THEN
      CALL mat_id(RMatdnb,basis_set%nb,basis_set%nb)
    ELSE
      IF (dnba_ind(1) == 0 .AND. dnba_ind(2) == 0) THEN
        RMatdnb(:,:) = basis_set%dnRGB%d0(:,:)
      ELSE IF (dnba_ind(1) == 0) THEN ! first derivative
        RMatdnb(:,:) = basis_set%dnRGB%d1(:,:,dnba_ind(2))
      ELSE IF (dnba_ind(2) == 0) THEN ! first derivative
        RMatdnb(:,:) = basis_set%dnRGB%d1(:,:,dnba_ind(1))
      ELSE ! 2d derivative
        RMatdnb(:,:) = basis_set%dnRGB%d2(:,:,dnba_ind(1),dnba_ind(2))
      END IF
    END IF


!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------
      END SUBROUTINE Get_MatdnRGB
      FUNCTION Get2_MatdnRGB(basis_set,dnba_ind) RESULT (RMatdnb)
      USE mod_system
      IMPLICIT NONE

      TYPE (basis),intent(in)        :: basis_set
      real(kind=Rkind), allocatable  :: RMatdnb(:,:)
      integer, intent(in)            :: dnba_ind(2)


      integer :: nq
!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='Get2_MatdnRGB'
!-----------------------------------------------------------
      IF (.NOT. basis_set%packed) RETURN
      nq = get_nq_FROM_basis(basis_set)
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nq',nq
      END IF
!-----------------------------------------------------------

    IF (NewBasisEl .AND. basis_set%ndim == 0) THEN
      allocate(RMatdnb(basis_set%nb,basis_set%nb))
      CALL mat_id(RMatdnb,basis_set%nb,basis_set%nb)
    ELSE
      IF (dnba_ind(1) == 0 .AND. dnba_ind(2) == 0) THEN
        RMatdnb = basis_set%dnRGB%d0(:,:)
      ELSE IF (dnba_ind(1) == 0) THEN ! first derivative
        RMatdnb = basis_set%dnRGB%d1(:,:,dnba_ind(2))
      ELSE IF (dnba_ind(2) == 0) THEN ! first derivative
        RMatdnb = basis_set%dnRGB%d1(:,:,dnba_ind(1))
      ELSE ! 2d derivative
        RMatdnb = basis_set%dnRGB%d2(:,:,dnba_ind(1),dnba_ind(2))
      END IF
    END IF


!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------
      END FUNCTION Get2_MatdnRGB

      FUNCTION Get2_MatdnRBG(basis_set) RESULT (RMatdnb)
      USE mod_system
      IMPLICIT NONE

      TYPE (basis),intent(in)        :: basis_set
      real(kind=Rkind), allocatable  :: RMatdnb(:,:)


      integer :: nq
!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='Get2_MatdnRBG'
!-----------------------------------------------------------
      IF (.NOT. basis_set%packed) RETURN
      nq = get_nq_FROM_basis(basis_set)
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nq',nq
      END IF
!-----------------------------------------------------------
    IF (NewBasisEl .AND. basis_set%ndim == 0) THEN
      allocate(RMatdnb(basis_set%nb,basis_set%nb))
      CALL mat_id(RMatdnb,basis_set%nb,basis_set%nb)
    ELSE
      RMatdnb = basis_set%dnRBG%d0(:,:)
    END IF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------
      END FUNCTION Get2_MatdnRBG

      SUBROUTINE Get_MatdnCGB(basis_set,CMatdnb,dnba_ind)
      USE mod_system
      IMPLICIT NONE

      TYPE (basis),intent(in)           :: basis_set
      complex(kind=Rkind),intent(inout) :: CMatdnb(:,:)
      integer, intent(in)               :: dnba_ind(2)


      integer :: nq
!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='Get_MatdnCGB'
!-----------------------------------------------------------
      IF (.NOT. basis_set%packed) RETURN
      nq = get_nq_FROM_basis(basis_set)
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nq',nq
      END IF
!-----------------------------------------------------------
    IF (NewBasisEl .AND. basis_set%ndim == 0) THEN
      CALL Cplx_mat_id(CMatdnb,basis_set%nb,basis_set%nb)
    ELSE
      IF (dnba_ind(1) == 0 .AND. dnba_ind(2) == 0) THEN
        CMatdnb(:,:) = basis_set%dnCGB%d0(:,:)
      ELSE IF (dnba_ind(1) == 0) THEN ! first derivative
        CMatdnb(:,:) = basis_set%dnCGB%d1(:,:,dnba_ind(2))
      ELSE IF (dnba_ind(2) == 0) THEN ! first derivative
        CMatdnb(:,:) = basis_set%dnCGB%d1(:,:,dnba_ind(1))
      ELSE ! 2d derivative
        CMatdnb(:,:) = basis_set%dnCGB%d2(:,:,dnba_ind(1),dnba_ind(2))
      END IF
    END IF


!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------
      END SUBROUTINE Get_MatdnCGB

      SUBROUTINE Get_MatdnRBB(basis_set,MatRBB,dnba_ind)
      USE mod_system
      IMPLICIT NONE

      TYPE (basis),intent(in)        :: basis_set
      real(kind=Rkind),intent(inout) :: MatRBB(:,:)
      integer, intent(in)            :: dnba_ind(2)

      integer :: nb

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='Get_MatdnRBB'
!-----------------------------------------------------------
      IF (.NOT. basis_set%packed) RETURN
      nb = basis_set%nb
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb',nb
        write(out_unitp,*) 'shape(MatRBB)',shape(MatRBB)
        write(out_unitp,*) 'alloc dnRBB',basis_set%dnRBB%alloc
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

      IF (.NOT. basis_set%dnRBB%alloc) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'basis_set%dnRBB is not allocated!!'
        write(out_unitp,*) 'CHECK the fortran source'
        STOP
      END IF

      IF (dnba_ind(1) == 0 .AND. dnba_ind(2) == 0) THEN
        CALL mat_id(MatRBB(:,:),nb,nb)
      ELSE IF (dnba_ind(1) == 0) THEN ! first derivative
        MatRBB(:,:) = basis_set%dnRBB%d1(:,:,dnba_ind(2))
      ELSE IF (dnba_ind(2) == 0) THEN ! first derivative
        MatRBB(:,:) = basis_set%dnRBB%d1(:,:,dnba_ind(1))
      ELSE ! 2d derivative
        MatRBB(:,:) = basis_set%dnRBB%d2(:,:,dnba_ind(1),dnba_ind(2))
      END IF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'MatRBB',dnba_ind
        CALL write_Mat(MatRBB,out_unitp,5)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------
      END SUBROUTINE Get_MatdnRBB
      SUBROUTINE Get_MatdnCBB(basis_set,MatCBB,dnba_ind)
      USE mod_system
      IMPLICIT NONE

      TYPE (basis),intent(in)           :: basis_set
      complex(kind=Rkind),intent(inout) :: MatCBB(:,:)
      integer, intent(in)               :: dnba_ind(2)

      integer :: nb

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='Get_MatdnCBB'
!-----------------------------------------------------------
      IF (.NOT. basis_set%packed) RETURN
      nb = basis_set%nb
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb',nb
        write(out_unitp,*) 'shape(MatCBB)',shape(MatCBB)
        write(out_unitp,*) 'alloc dnCBB',basis_set%dnCBB%alloc
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

      IF (.NOT. basis_set%dnCBB%alloc) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'basis_set%dnCBB is not allocated!!'
        write(out_unitp,*) 'CHECK the fortran source'
        STOP
      END IF

      IF (dnba_ind(1) == 0 .AND. dnba_ind(2) == 0) THEN
        CALL Cplx_mat_id(MatCBB(:,:),nb,nb)
      ELSE IF (dnba_ind(1) == 0) THEN ! first derivative
        MatCBB(:,:) = basis_set%dnCBB%d1(:,:,dnba_ind(2))
      ELSE IF (dnba_ind(2) == 0) THEN ! first derivative
        MatCBB(:,:) = basis_set%dnCBB%d1(:,:,dnba_ind(1))
      ELSE ! 2d derivative
        MatCBB(:,:) = basis_set%dnCBB%d2(:,:,dnba_ind(1),dnba_ind(2))
      END IF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'MatCBB',dnba_ind
        CALL write_Mat(MatCBB,out_unitp,5)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------
      END SUBROUTINE Get_MatdnCBB
!================================================================
!       write the type basis
!================================================================
      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      RECURSIVE SUBROUTINE RecWrite_basis(basis_set,write_all)
      USE mod_system
      USE mod_MPI
      IMPLICIT NONE

       TYPE (basis) :: basis_set
       logical, intent(in), optional :: write_all

       logical                      :: write_all_loc
       integer                      :: i,j,L,ib,val,li,ui,lj,uj,nb_basis,n1,nq
       integer,allocatable              :: tab_nq(:)
       character(len=:), allocatable     :: Rec_line

       integer :: err_mem,memory
       character (len=*), parameter :: Rec_tab = '===='
       integer, save                :: iRec = -1


       iRec = iRec + 1

       Rec_line = String_TO_String(Rec_tab)
       DO i=2,iRec
         Rec_line = String_TO_String(Rec_line // Rec_tab)
       END DO
       Rec_line = String_TO_String(Rec_line // '=Rec' // int_TO_char(iRec) // '=')

       write_all_loc = .FALSE.
       IF (present(write_all)) write_all_loc = write_all

       nb_basis = basis_set%nb_basis
       write(out_unitp,*) Rec_line,'BEGINNING RecWrite_basis'
       write(out_unitp,*) Rec_line,'write_all: ',write_all_loc
       write(out_unitp,*)

       write(out_unitp,*) Rec_line,'ndim',basis_set%ndim
       write(out_unitp,*) Rec_line,'nb,nq',basis_set%nb,          &
                                            get_nq_FROM_basis(basis_set)
       write(out_unitp,*) Rec_line,'nb_init,nq_init (before contraction)',&
                     basis_set%nb_init,get_nq_FROM_basis(basis_set,init=.TRUE.)

       write(out_unitp,*) Rec_line,'Nested,nq_max_Nested',basis_set%Nested,basis_set%nq_max_Nested

       CALL Write_Basis_Grid_Param(basis_set%Basis_Grid_Para,Rec_line)

       write(out_unitp,*) Rec_line,'print_info_OF_basisDP',basis_set%print_info_OF_basisDP
       IF (allocated(basis_set%nrho)) write(out_unitp,*) Rec_line,'nrho_OF_Qbasis',basis_set%nrho(:)

       write(out_unitp,*)
       CALL Write_SymAbelian(basis_set%P_SymAbelian)

       write(out_unitp,*) Rec_line,'check_basis',basis_set%check_basis
       write(out_unitp,*) Rec_line,'check_nq_OF_basis',basis_set%check_nq_OF_basis
       write(out_unitp,*) Rec_line,'packed,packed_done',basis_set%packed,basis_set%packed_done
       write(out_unitp,*) Rec_line,'primitive',basis_set%primitive
       write(out_unitp,*) Rec_line,'primitive_done',basis_set%primitive_done

       write(out_unitp,*) Rec_line,'contrac',basis_set%contrac
       write(out_unitp,*) Rec_line,'auto_contrac',basis_set%auto_contrac
       write(out_unitp,*) Rec_line,'contrac_analysis',basis_set%contrac_analysis
       write(out_unitp,*) Rec_line,'POGridRep,POGridRep_polyortho',basis_set%POGridRep,basis_set%POGridRep_polyortho
       write(out_unitp,*) Rec_line,'nqPLUSnbc_TO_nqc',basis_set%nqPLUSnbc_TO_nqc
       write(out_unitp,*) Rec_line,'max_ene_contrac',basis_set%max_ene_contrac
       write(out_unitp,*) Rec_line,'min_nbc,max_nbc',basis_set%min_nbc,basis_set%max_nbc
       write(out_unitp,*) Rec_line,'auto_contrac_type1_TO',basis_set%auto_contrac_type1_TO
       write(out_unitp,*) Rec_line,'auto_contrac_type21_TO',basis_set%auto_contrac_type21_TO
       write(out_unitp,*) Rec_line,'nbc,nqc',basis_set%nbc,basis_set%nqc
       write(out_unitp,*) Rec_line,'read_contrac_file',basis_set%read_contrac_file
       IF (basis_set%read_contrac_file)                                      &
          write(out_unitp,*) Rec_line,'name_contract_file',basis_set%file_contrac%name
       write(out_unitp,*)

       write(out_unitp,*) Rec_line,'type,name',basis_set%type,basis_set%name
       write(out_unitp,*)
       write(out_unitp,*) Rec_line,'cplx ',basis_set%cplx
       write(out_unitp,*) Rec_line,'dnBBRep,dnBBRep_done',basis_set%dnBBRep,basis_set%dnBBRep_done

       write(out_unitp,*)
       IF ( allocated(basis_set%iQdyn) ) THEN
         write(out_unitp,*) Rec_line,'iQdyn',basis_set%iQdyn(1:basis_set%ndim)
         write(out_unitp,*)
       END IF
       IF (allocated(basis_set%Tabder_Qdyn_TO_Qbasis) ) THEN
         write(out_unitp,*) Rec_line,'Tabder_Qdyn_TO_Qbasis',basis_set%Tabder_Qdyn_TO_Qbasis(:)
       END IF

       write(out_unitp,*) Rec_line,'alloc tab_ndim_index',allocated(basis_set%tab_ndim_index)
       IF (allocated(basis_set%tab_ndim_index) ) THEN
         DO ib=1,basis_set%nb
           write(out_unitp,*) 'ib,index',ib,basis_set%tab_ndim_index(:,ib)
         END DO
       END IF

       write(out_unitp,*) Rec_line,'------------------------------------------------'
       write(out_unitp,*) Rec_line,'------------------------------------------------'
       write(out_unitp,*) Rec_line,'------------------------------------------------'
       write(out_unitp,*) Rec_line,'SparseGrid_type        :',basis_set%SparseGrid_type
       write(out_unitp,*) Rec_line,'SparseGrid_With_Cuba   :',basis_set%SparseGrid_With_Cuba
       write(out_unitp,*) Rec_line,'SparseGrid_With_Smolyak:',basis_set%SparseGrid_With_Smolyak
       write(out_unitp,*) Rec_line,'SparseGrid_With_DP     :',basis_set%SparseGrid_With_DP
       write(out_unitp,*) Rec_line,'nb_basis',basis_set%nb_basis

       write(out_unitp,*) Rec_line,'------------------------------------------------'
       write(out_unitp,*) Rec_line,'- nDindG ---------------------------------------'
       CALL Write_nDindex(basis_set%nDindG,name_info=Rec_line)

       write(out_unitp,*) Rec_line,'- nDindB ---------------------------------------'
       write(out_unitp,*) Rec_line,'Type_OF_nDindB',basis_set%Type_OF_nDindB
       write(out_unitp,*) Rec_line,'Norm_OF_nDindB',basis_set%Norm_OF_nDindB
       write(out_unitp,*) Rec_line,'weight_OF_nDindB',basis_set%weight_OF_nDindB
       write(out_unitp,*) Rec_line,'nDinit_OF_nDindB',basis_set%nDinit_OF_nDindB
       write(out_unitp,*) Rec_line,'MaxCoupling_OF_nDindB',basis_set%MaxCoupling_OF_nDindB
       write(out_unitp,*) Rec_line,'nb_OF_MinNorm_OF_nDindB',basis_set%nb_OF_MinNorm_OF_nDindB
       write(out_unitp,*) Rec_line,'Div_nb_TO_Norm_OF_nDindB',basis_set%Div_nb_TO_Norm_OF_nDindB
       write(out_unitp,*) Rec_line,'contrac_WITH_nDindB',basis_set%contrac_WITH_nDindB
       write(out_unitp,*) Rec_line,'asso nDindB',associated(basis_set%nDindB)

       IF (associated(basis_set%nDindB)) CALL Write_nDindex(basis_set%nDindB,name_info=Rec_line)


       write(out_unitp,*) Rec_line,'tab_basis_done',basis_set%tab_basis_done
       write(out_unitp,*) Rec_line,' tab_basis_linked',basis_set%tab_basis_linked
       IF (basis_set%nb_basis > 0 .AND. .NOT. basis_set%packed) THEN
         write(out_unitp,*) Rec_line,'------------------------------------------------'
         write(out_unitp,*) Rec_line,'- tab_Pbasis for direct_product or SparseBasis -'
         write(out_unitp,*) Rec_line,'------------------------------------------------'
         DO i=1,basis_set%nb_basis
           IF (write_all_loc) THEN
             write(out_unitp,*) Rec_line,'  tab_Pbasis:',i
             CALL RecWrite_basis(basis_set%tab_Pbasis(i)%Pbasis,write_all_loc)
           ELSE
             write(out_unitp,*) Rec_line,'  tab_Pbasis:',i,' is not written'
           END IF
         END DO
         write(out_unitp,*) Rec_line,'------------------------------------------------'
         write(out_unitp,*) Rec_line,'- END tab_Pbasis --------------------------------'
         write(out_unitp,*) Rec_line,'------------------------------------------------'

         write(out_unitp,*) Rec_line,'------------------------------------------------'
         write(out_unitp,*) Rec_line,'- Tab_OF_Tabnb2 for SparseBasis ----------------'
         write(out_unitp,*) Rec_line,'------------------------------------------------'
         IF (associated(basis_set%Tab_OF_Tabnb2)) THEN
           DO i=1,basis_set%nb_basis
             IF (associated(basis_set%Tab_OF_Tabnb2(i)%vec) ) THEN
               write(out_unitp,*) Rec_line,'ibasis,nbb',i,basis_set%Tab_OF_Tabnb2(i)%nb_var_vec
               write(out_unitp,*) Rec_line,'Tab_OF_Tabnb2(ibasis)',basis_set%Tab_OF_Tabnb2(i)%vec
             END IF
           END DO
         END IF
         write(out_unitp,*) Rec_line,'------------------------------------------------'
         write(out_unitp,*) Rec_line,'- END Tab_OF_Tabnb2-----------------------------'
         write(out_unitp,*) Rec_line,'------------------------------------------------'
       END IF

       write(out_unitp,*) Rec_line,'------------------------------------------------'
       write(out_unitp,*) Rec_line,'- for Sparse Grid (Basis) ----------------------'
       write(out_unitp,*) Rec_line,'------------------------------------------------'
       write(out_unitp,*) Rec_line,' SparseGrid_type',basis_set%SparseGrid_type
       write(out_unitp,*) Rec_line,' With_L',basis_set%With_L
       write(out_unitp,*) Rec_line,' L_SparseGrid ',basis_set%L_SparseGrid
       write(out_unitp,*) Rec_line,' L_SparseBasis',basis_set%L_SparseBasis
       write(out_unitp,*) Rec_line,' Parameters: nq(L) = A + B * L**expo'
       CALL Write_Basis_L_TO_n(basis_set%L_TO_nq,Rec_line)
       write(out_unitp,*) Rec_line,' Parameters: nb(L) = A + B * L**expo'
       CALL Write_Basis_L_TO_n(basis_set%L_TO_nb,Rec_line)

       IF (basis_set%SparseGrid_type == 1) THEN

         CALL alloc_NParray(tab_nq,(/ nb_basis /),"tab_nq","RecWrite_basis")
         DO i=1,basis_set%nb_SG
           DO j=1,nb_basis
             tab_nq(j) =                                                &
              get_nq_FROM_basis(basis_set%tab_PbasisSG(i)%Pbasis%tab_Pbasis(j)%Pbasis)
           END DO
           write(out_unitp,*) Rec_line,'  tab_PbasisSG + WeightSG:',&
                                                  i,basis_set%WeightSG(i)
           write(out_unitp,*) Rec_line,'  nq_tot,nq(:):',         &
                   get_nq_FROM_basis(basis_set%tab_PbasisSG(i)%Pbasis), &
                                                         ' : ',tab_nq(:)
           write(out_unitp,*) Rec_line,'  tab_PbasisSG:',i,' not written'
           !CALL RecWrite_basis(basis_set%tab_PbasisSG(i)%Pbasis,write_all_loc)
         END DO
         CALL dealloc_NParray(tab_nq,"tab_nq","RecWrite_basis")

         IF (associated(basis_set%tab_basisPrimSG)) THEN
         write(out_unitp,*) Rec_line,'------------------------------------------------'
         write(out_unitp,*) Rec_line,'- tab_basisPrimSG ------------------------------'
         write(out_unitp,*) Rec_line,'------------------------------------------------'
         DO L=0,ubound(basis_set%tab_basisPrimSG,dim=2)
         DO i=1,ubound(basis_set%tab_basisPrimSG,dim=1)
           IF (write_all_loc) THEN
             write(out_unitp,*) Rec_line,'  tab_basisPrimSG(i,L): ',i,L
             CALL RecWrite_basis(basis_set%tab_basisPrimSG(i,L),write_all_loc)
           ELSE
             write(out_unitp,*) Rec_line,'  tab_basisPrimSG(i,L):',i,L,' is not written'
           END IF
         END DO
         END DO
         write(out_unitp,*) Rec_line,'------------------------------------------------'
         write(out_unitp,*) Rec_line,'- END tab_basisPrimSG --------------------------'
         write(out_unitp,*) Rec_line,'------------------------------------------------'
         END IF

         write(out_unitp,*) Rec_line,'------------------------------------------------'
         write(out_unitp,*) Rec_line,'- END Sparse Grid (Basis) ----------------------'
         write(out_unitp,*) Rec_line,'------------------------------------------------'
       END IF

       IF (basis_set%nb_Transfo > 0 .AND. allocated(basis_set%cte_Transfo)) THEN
         DO i=1,basis_set%nb_Transfo
           write(out_unitp,*) Rec_line,'cte_Transfo',i,basis_set%cte_Transfo(:,i)
         END DO
       END IF

       IF (basis_set%ndim > 0 ) THEN
         IF ( allocated(basis_set%A) .AND. allocated(basis_set%opt_A)) THEN
           write(out_unitp,*) Rec_line,'A    ',basis_set%A
           write(out_unitp,*) Rec_line,'opt_A',basis_set%opt_A
         END IF
         IF ( allocated(basis_set%B) .AND. allocated(basis_set%opt_B)) THEN
           write(out_unitp,*) Rec_line,'B    ',basis_set%B
           write(out_unitp,*) Rec_line,'opt_B',basis_set%opt_B
         END IF
         IF ( allocated(basis_set%Q0) .AND. allocated(basis_set%opt_Q0)) THEN
           write(out_unitp,*) Rec_line,'Q0    ',basis_set%Q0
           write(out_unitp,*) Rec_line,'opt_Q0',basis_set%opt_Q0
         END IF
         IF ( allocated(basis_set%scaleQ) .AND. allocated(basis_set%opt_scaleQ)) THEN
           write(out_unitp,*) Rec_line,'scaleQ    ',basis_set%scaleQ
           write(out_unitp,*) Rec_line,'opt_scaleQ',basis_set%opt_scaleQ
         END IF
       END IF
       IF(MPI_id==0) write(out_unitp,*) 'Parameter(s) to be optimized?: ',basis_set%opt_param


       nq = get_nq_FROM_basis(basis_set)
       IF (nq > 0 .AND.  basis_set%ndim > 0  .AND. write_all_loc) THEN
         IF ( allocated(basis_set%x) ) THEN
           write(out_unitp,*) Rec_line,'x'
           CALL Write_Mat(basis_set%x,out_unitp,8,name_info=Rec_line)
           write(out_unitp,*)
         END IF
         IF ( allocated(basis_set%w) ) THEN
           write(out_unitp,*) Rec_line,'w'
           CALL Write_Vec(basis_set%w,out_unitp,8,name_info=Rec_line)
           write(out_unitp,*)
         END IF
         IF ( allocated(basis_set%wrho) ) THEN
           write(out_unitp,*) Rec_line,'w*rho'
           CALL Write_Vec(basis_set%wrho,out_unitp,8,name_info=Rec_line)
           write(out_unitp,*)
         END IF
         IF ( allocated(basis_set%rho) ) THEN
           write(out_unitp,*) Rec_line,'rho'
           CALL Write_Vec(basis_set%rho,out_unitp,8,name_info=Rec_line)
           write(out_unitp,*)
         END IF

       END IF

       IF (basis_set%nb > 0 .AND. nq > 0 .AND. basis_set%ndim > 0       &
                                         .AND. write_all_loc) THEN

         write(out_unitp,*) Rec_line,'dnRGB'
         CALL Write_dnSVM(basis_set%dnRGB)
         write(out_unitp,*) Rec_line,'dnRBG'
         CALL Write_dnSVM(basis_set%dnRBG)

         write(out_unitp,*) Rec_line,'dnCGB'
         CALL Write_dnSVM(basis_set%dnCGB)
         write(out_unitp,*) Rec_line,'dnCBG'
         CALL Write_dnSVM(basis_set%dnCBG)

         write(out_unitp,*) Rec_line,'dnGGRep:      ',basis_set%dnGGRep
         write(out_unitp,*) Rec_line,'dnGGRep_done: ',basis_set%dnGGRep_done
         write(out_unitp,*) Rec_line,'alloc dnRGG:  ',basis_set%dnRGG%alloc
         CALL Write_dnMat(basis_set%dnRGG)

         write(out_unitp,*) Rec_line,'dnRBB'
         CALL Write_dnMat(basis_set%dnRBB)
         write(out_unitp,*) Rec_line,'dnCBB'
         CALL Write_dnCplxMat(basis_set%dnCBB)

       END IF



         IF (allocated(basis_set%EneH0)) THEN
           IF(MPI_id==0) write(out_unitp,*) Rec_line,'EneH0 = <d0b(:,ib) I H0 I d0b(:,ib)>'
           CALL Write_Vec(basis_set%EneH0,out_unitp,8,name_info=Rec_line)
         END IF

       write(out_unitp,*) Rec_line,'END RecWrite_basis'
       iRec = iRec - 1
       deallocate(Rec_line)

       END SUBROUTINE RecWrite_basis

      RECURSIVE SUBROUTINE RecWriteMini_basis(basis_set)
      USE mod_system
      IMPLICIT NONE

       TYPE (basis) :: basis_set

       integer             :: i,j,ib,val,li,ui,lj,uj,nb_basis,n1,L,nq,nq_init

       character(len=:), allocatable     :: Rec_line

       character (len=*), parameter :: Rec_tab = '===='
       integer, save                :: iRec = -1

       iRec = iRec + 1

       Rec_line = String_TO_String(Rec_tab)
       DO i=2,iRec
         Rec_line = String_TO_String(Rec_line // Rec_tab)
       END DO
       Rec_line = String_TO_String(Rec_line // '=RecMini' // int_TO_char(iRec) // '=')


       nq      = get_nq_FROM_basis(basis_set)
       nq_init = get_nq_FROM_basis(basis_set,init=.TRUE.)

       nb_basis = basis_set%nb_basis
       write(out_unitp,*) Rec_line,'BEGINNING RecWriteMini_basis'
       write(out_unitp,*)
       write(out_unitp,*) Rec_line,'ndim',basis_set%ndim
       write(out_unitp,*) Rec_line,'nb,nq',basis_set%nb,nq
       write(out_unitp,*) Rec_line,'nb,nq (before contraction)',basis_set%nb_init,nq_init

       write(out_unitp,*) Rec_line,'Nested,nq_max_Nested',basis_set%Nested,basis_set%nq_max_Nested
       CALL Write_Basis_Grid_Param(basis_set%Basis_Grid_Para,Rec_line)



       write(out_unitp,*) Rec_line,'print_info_OF_basisDP',basis_set%print_info_OF_basisDP

       write(out_unitp,*) Rec_line,'alloc nrho_OF_Qbasis',allocated(basis_set%nrho)
       IF (allocated(basis_set%nrho)) write(out_unitp,*) Rec_line,'nrho_OF_Qbasis',basis_set%nrho(:)
       write(out_unitp,*)
       write(out_unitp,*) Rec_line,'asso SymAbelian',associated(basis_set%P_SymAbelian)


       write(out_unitp,*) Rec_line,'check_basis',basis_set%check_basis
       write(out_unitp,*) Rec_line,'check_nq_OF_basis',basis_set%check_nq_OF_basis
       write(out_unitp,*) Rec_line,'packed,packed_done',basis_set%packed,basis_set%packed_done
       write(out_unitp,*) Rec_line,'primitive',basis_set%primitive
       write(out_unitp,*) Rec_line,'primitive_done',basis_set%primitive_done

       write(out_unitp,*) Rec_line,'contrac',basis_set%contrac
       write(out_unitp,*) Rec_line,'auto_contrac',basis_set%auto_contrac
       write(out_unitp,*) Rec_line,'contrac_analysis',basis_set%contrac_analysis
       write(out_unitp,*) Rec_line,'POGridRep,POGridRep_polyortho',basis_set%POGridRep,basis_set%POGridRep_polyortho
       write(out_unitp,*) Rec_line,'nqPLUSnbc_TO_nqc',basis_set%nqPLUSnbc_TO_nqc
       write(out_unitp,*) Rec_line,'max_ene_contrac',basis_set%max_ene_contrac
       write(out_unitp,*) Rec_line,'min_nbc,max_nbc',basis_set%min_nbc,basis_set%max_nbc
       write(out_unitp,*) Rec_line,'auto_contrac_type1_TO',basis_set%auto_contrac_type1_TO
       write(out_unitp,*) Rec_line,'auto_contrac_type21_TO',basis_set%auto_contrac_type21_TO
       write(out_unitp,*) Rec_line,'nbc,nqc',basis_set%nbc,basis_set%nqc
       write(out_unitp,*) Rec_line,'read_contrac_file',basis_set%read_contrac_file
       IF (basis_set%read_contrac_file)                                      &
          write(out_unitp,*) Rec_line,'name_contract_file',basis_set%file_contrac%name
       write(out_unitp,*)

       write(out_unitp,*) Rec_line,'type,name',basis_set%type,basis_set%name
       write(out_unitp,*)
       write(out_unitp,*) Rec_line,'cplx ',basis_set%cplx
       write(out_unitp,*) Rec_line,'dnBBRep,dnBBRep_done',basis_set%dnBBRep,basis_set%dnBBRep_done

       write(out_unitp,*)
       write(out_unitp,*) Rec_line,'alloc iQdyn',allocated(basis_set%iQdyn)
       IF ( allocated(basis_set%iQdyn) ) THEN
         write(out_unitp,*) Rec_line,'iQdyn',basis_set%iQdyn(1:basis_set%ndim)
         write(out_unitp,*)
       END IF
       write(out_unitp,*) Rec_line,'alloc Tabder_Qdyn_TO_Qbasis',allocated(basis_set%Tabder_Qdyn_TO_Qbasis)
       IF (allocated(basis_set%Tabder_Qdyn_TO_Qbasis) ) THEN
         write(out_unitp,*) Rec_line,'Tabder_Qdyn_TO_Qbasis',basis_set%Tabder_Qdyn_TO_Qbasis(:)
       END IF
       write(out_unitp,*) Rec_line,'alloc tab_ndim_index',allocated(basis_set%tab_ndim_index)

       IF (basis_set%nb_Transfo > 0 .AND. allocated(basis_set%cte_Transfo)) THEN
         DO i=1,basis_set%nb_Transfo
           write(out_unitp,*) Rec_line,'cte_Transfo',i,basis_set%cte_Transfo(:,i)
         END DO
       END IF

       write(out_unitp,*) Rec_line,'alloc A,B,Q0,scaleQ',         &
                       allocated(basis_set%A),allocated(basis_set%B),   &
                   allocated(basis_set%Q0),allocated(basis_set%scaleQ)
       write(out_unitp,*) Rec_line,'alloc opt A,B,Q0,scaleQ',     &
               allocated(basis_set%opt_A),allocated(basis_set%opt_B),   &
            allocated(basis_set%opt_Q0),allocated(basis_set%opt_scaleQ)

       write(out_unitp,*) Rec_line,'alloc x,w,wrho',allocated(basis_set%x), &
                       allocated(basis_set%w),allocated(basis_set%wrho)

       write(out_unitp,*) Rec_line,'dnGGRep:      ',basis_set%dnGGRep
       write(out_unitp,*) Rec_line,'dnGGRep_done: ',basis_set%dnGGRep_done
       write(out_unitp,*) Rec_line,'alloc dnRGG:  ',basis_set%dnRGG%alloc
       write(out_unitp,*) Rec_line,'alloc dnRGB:  ',basis_set%dnRGB%alloc
       write(out_unitp,*) Rec_line,'alloc dnCGB:  ',basis_set%dnCGB%alloc
       write(out_unitp,*) Rec_line,'alloc dnRBB:  ',basis_set%dnRBB%alloc
       write(out_unitp,*) Rec_line,'alloc dnCBB:  ',basis_set%dnCBB%alloc



       write(out_unitp,*) Rec_line,'------------------------------------------------'
       write(out_unitp,*) Rec_line,'------------------------------------------------'
       write(out_unitp,*) Rec_line,'------------------------------------------------'
       write(out_unitp,*) Rec_line,'SparseGrid_type        :',basis_set%SparseGrid_type
       write(out_unitp,*) Rec_line,'SparseGrid_With_Cuba   :',basis_set%SparseGrid_With_Cuba
       write(out_unitp,*) Rec_line,'SparseGrid_With_Smolyak:',basis_set%SparseGrid_With_Smolyak
       write(out_unitp,*) Rec_line,'SparseGrid_With_DP     :',basis_set%SparseGrid_With_DP
       write(out_unitp,*) Rec_line,'nb_basis',basis_set%nb_basis

       write(out_unitp,*) Rec_line,'------------------------------------------------'
       write(out_unitp,*) Rec_line,'- nDindG ---------------------------------------'
       CALL Write_nDindex(basis_set%nDindG,name_info=Rec_line)

       write(out_unitp,*) Rec_line,'- nDindB ---------------------------------------'
       write(out_unitp,*) Rec_line,'Type_OF_nDindB',basis_set%Type_OF_nDindB
       write(out_unitp,*) Rec_line,'Norm_OF_nDindB',basis_set%Norm_OF_nDindB
       write(out_unitp,*) Rec_line,'weight_OF_nDindB',basis_set%weight_OF_nDindB
       write(out_unitp,*) Rec_line,'nDinit_OF_nDindB',basis_set%nDinit_OF_nDindB
       write(out_unitp,*) Rec_line,'MaxCoupling_OF_nDindB',basis_set%MaxCoupling_OF_nDindB
       write(out_unitp,*) Rec_line,'nb_OF_MinNorm_OF_nDindB',basis_set%nb_OF_MinNorm_OF_nDindB
       write(out_unitp,*) Rec_line,'Div_nb_TO_Norm_OF_nDindB',basis_set%Div_nb_TO_Norm_OF_nDindB
       write(out_unitp,*) Rec_line,'contrac_WITH_nDindB',basis_set%contrac_WITH_nDindB
       CALL Write_nDindex(basis_set%nDindB,name_info=Rec_line)


       write(out_unitp,*)
       write(out_unitp,*) Rec_line,'------------------------------------------------'
       write(out_unitp,*) Rec_line,'- tab_Pbasis for direct_product or SparseBasis -'
       write(out_unitp,*) Rec_line,'------------------------------------------------'
       write(out_unitp,*) Rec_line,' tab_basis_done',basis_set%tab_basis_done
       write(out_unitp,*) Rec_line,' tab_basis_linked',basis_set%tab_basis_linked

       write(out_unitp,*) Rec_line,' asso tab_Pbasis:',associated(basis_set%tab_Pbasis)
       IF ( associated(basis_set%tab_Pbasis) ) THEN
         DO i=1,basis_set%nb_basis
           write(out_unitp,*) Rec_line,'  tab_Pbasis:',i
           CALL RecWriteMini_basis(basis_set%tab_Pbasis(i)%Pbasis)
         END DO
       END IF
       write(out_unitp,*) Rec_line,'------------------------------------------------'
       write(out_unitp,*) Rec_line,'- END tab_Pbasis --------------------------------'
       write(out_unitp,*) Rec_line,'------------------------------------------------'

       write(out_unitp,*) Rec_line,'------------------------------------------------'
       write(out_unitp,*) Rec_line,'- Tab_OF_Tabnb2 for SparseBasis ----------------'
       write(out_unitp,*) Rec_line,'------------------------------------------------'
       write(out_unitp,*) Rec_line,'asso Tab_OF_Tabnb2',associated(basis_set%Tab_OF_Tabnb2)
       IF (associated(basis_set%Tab_OF_Tabnb2)) THEN
         DO i=1,basis_set%nb_basis
           write(out_unitp,*) Rec_line,'asso Tab_OF_Tabnb2(i)%vec',i,associated(basis_set%Tab_OF_Tabnb2(i)%vec)
           write(out_unitp,*) Rec_line,'ibasis,nbb',i,basis_set%Tab_OF_Tabnb2(i)%nb_var_vec
         END DO
       END IF
       write(out_unitp,*) Rec_line,'------------------------------------------------'
       write(out_unitp,*) Rec_line,'- END Tab_OF_Tabnb2-----------------------------'
       write(out_unitp,*) Rec_line,'------------------------------------------------'

       write(out_unitp,*) Rec_line,'------------------------------------------------'
       write(out_unitp,*) Rec_line,'- for Sparse Grid (Basis) ----------------------'
       write(out_unitp,*) Rec_line,'------------------------------------------------'
       write(out_unitp,*) Rec_line,' With_L',basis_set%With_L
       write(out_unitp,*) Rec_line,' Parameters: nq(L) = A + B * L**expo'
       CALL Write_Basis_L_TO_n(basis_set%L_TO_nq,Rec_line)
       write(out_unitp,*) Rec_line,' Parameters: nb(L) = A + B * L**expo'
       CALL Write_Basis_L_TO_n(basis_set%L_TO_nb,Rec_line)

       write(out_unitp,*) Rec_line,' asso basis_set%tab_PbasisSG',associated(basis_set%tab_PbasisSG)
       IF (associated(basis_set%tab_PbasisSG)) THEN
         DO i=1,basis_set%nb_SG
           write(out_unitp,*) Rec_line,'  tab_PbasisSG + WeightSG:',i,basis_set%WeightSG(i)
           CALL RecWriteMini_basis(basis_set%tab_PbasisSG(i)%Pbasis)
         END DO

       END IF
       write(out_unitp,*) Rec_line,'------------------------------------------------'
       write(out_unitp,*) Rec_line,'- END Sparse Grid (Basis) ----------------------'
       write(out_unitp,*) Rec_line,'------------------------------------------------'

       write(out_unitp,*) Rec_line,'END RecWriteMini_basis'
       iRec = iRec - 1

       deallocate(Rec_line)
       END SUBROUTINE RecWriteMini_basis

!=======================================================================
!
!      para_AllBasis subroutines
!
!=======================================================================
      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      SUBROUTINE alloc_AllBasis(para_AllBasis)

         TYPE (param_AllBasis), intent(inout) :: para_AllBasis

         integer      :: i,err_mem,memory

         IF (para_AllBasis%alloc) THEN
           write(out_unitp,*) ' ERROR in alloc_AllBasis'
           write(out_unitp,*) 'para_AllBasis is already allocated!'
           write(out_unitp,*) ' CHECK the source!!!!!'
           STOP
         END IF

         para_AllBasis%alloc = .TRUE.

         para_AllBasis%nb_be = 1
         para_AllBasis%nb_bi = 1

         CALL alloc_array(para_AllBasis%BasisnD,'para_AllBasis%BasisnD','alloc_AllBasis')
         CALL alloc_array(para_AllBasis%Basis2n,'para_AllBasis%Basis2n','alloc_AllBasis')


      END SUBROUTINE alloc_AllBasis

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      SUBROUTINE dealloc_AllBasis(para_AllBasis)

         TYPE(param_AllBasis), intent(inout) :: para_AllBasis
         integer      :: i,err_mem,memory

         para_AllBasis%alloc = .FALSE.

         para_AllBasis%nb_be = 1
         para_AllBasis%nb_bi = 1

         CALL dealloc_basis(para_AllBasis%BasisnD)
         CALL dealloc_array(para_AllBasis%BasisnD,'para_AllBasis%BasisnD','alloc_AllBasis')


         CALL dealloc_basis(para_AllBasis%Basis2n)
         CALL dealloc_array(para_AllBasis%Basis2n,'para_AllBasis%Basis2n','alloc_AllBasis')


         nullify(para_AllBasis%BasisElec)
         nullify(para_AllBasis%BasisRot)

      END SUBROUTINE dealloc_AllBasis

      SUBROUTINE Set_nq_OF_basis(basis_set,nq)

      TYPE (basis), intent(inout) :: basis_set
      integer, intent (in)        :: nq

!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Set_nq_OF_basis'
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------

      basis_set%Basis_Grid_Para%nq = nq

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE Set_nq_OF_basis

      FUNCTION get_nq_FROM_basis(basis_set,LGrid,init) RESULT(nq)

      TYPE (basis), intent(in) :: basis_set
      integer, optional :: LGrid
      logical, optional :: init
      integer           :: nq

!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='get_nq_FROM_basis'
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL RecWriteMini_basis(basis_set)
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------

      IF (basis_set%ndim == 0) THEN ! for the electronic basis
        nq = get_nb_FROM_basis(basis_set)
      ELSE
        IF (present(init)) THEN
          nq  = basis_set%Basis_Grid_Para%nq_init
        ELSE
          nq  = basis_set%Basis_Grid_Para%nq
        END IF
      END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'nq',nq
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------

      END FUNCTION get_nq_FROM_basis

      RECURSIVE FUNCTION get_nqa_FROM_basis(basis_set) RESULT(nqa)

      TYPE (basis), intent(in) :: basis_set
      integer           :: nqa


      integer           :: ib,L,nqa_SG,i_SG

!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='get_nqa_FROM_basis'
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL RecWriteMini_basis(basis_set)
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------
    IF (basis_set%ndim == 0 .AND. NewBasisEl) THEN
      nqa = 1 ! Correct value ???
    ELSE
      IF (basis_set%packed_done) THEN
        nqa = get_nq_FROM_basis(basis_set)
      ELSE
        IF (basis_set%nb_basis == 0 ) STOP ' ERROR with packed!!!'

        SELECT CASE (basis_set%SparseGrid_type)
        CASE (0) ! Direct product
          nqa = 1
          DO ib=1,basis_set%nb_basis
            nqa = nqa * get_nqa_FROM_basis(basis_set%tab_Pbasis(ib)%Pbasis)
          END DO

        CASE (1) ! Sparse basis (Smolyak 1st implementation)
          ! Just the grid point number for the SG:
          nqa = 0
          DO i_SG=1,basis_set%nb_SG
            nqa = nqa + get_nqa_FROM_basis(basis_set%tab_PbasisSG(i_SG)%Pbasis)
          END DO

        CASE (2,4) ! Sparse basis (Smolyak 2d  implementation)
          nqa = 0
          DO i_SG=1,basis_set%nb_SG

            nqa_SG = 1
            DO ib=1,basis_set%nb_basis
              L = basis_set%para_SGType2%nDind_SmolyakRep%Tab_nDval(ib,i_SG)
              nqa_SG = nqa_SG * get_nqa_FROM_basis(basis_set%tab_basisPrimSG(L,ib))
            END DO

            nqa = nqa + nqa_SG
          END DO

        CASE DEFAULT
          write(out_unitp,*) ' ERROR in',name_sub
          write(out_unitp,*) ' WRONG SparseGrid_type',basis_set%SparseGrid_type
          write(out_unitp,*) ' The possibilities are: 0, 1, 2, 4'
          STOP
        END SELECT
      END IF
    END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'nqa',nqa
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------

      END FUNCTION get_nqa_FROM_basis

      RECURSIVE SUBROUTINE get_tab_nq_OF_Qact(tab_nq,basis_set)

      TYPE (basis), intent(in) :: basis_set
      integer, intent (inout)  :: tab_nq(basis_set%ndim)

      integer                  :: tab1_nq(basis_set%ndim)
      integer                  :: idim,ndim,ib,L,i_SG
!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='get_tab_nq_OF_Qact'
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL RecWriteMini_basis(basis_set)
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------

      tab_nq(:) = 0
      IF (basis_set%primitive) THEN
        DO idim=1,basis_set%ndim
          tab_nq(idim) = get_nq_FROM_basis(basis_set)
        END DO

      ELSE
        IF (basis_set%nb_basis == 0 ) STOP ' ERROR with packed!!!'

        SELECT CASE (basis_set%SparseGrid_type)
        CASE (0) ! Direct product
          idim = 0
          DO ib=1,basis_set%nb_basis
            ndim = basis_set%tab_Pbasis(ib)%Pbasis%ndim
            CALL get_tab_nq_OF_Qact(tab_nq(idim+1:idim+ndim),           &
                                        basis_set%tab_Pbasis(ib)%Pbasis)
            idim = idim+ndim
          END DO
        CASE (1) ! Sparse basis (Smolyak 1st implementation)
          idim = 0
          DO ib=1,basis_set%nb_basis
            ndim = basis_set%tab_basisPrimSG(0,ib)%ndim
            DO L=0,basis_set%L_SparseGrid
               CALL get_tab_nq_OF_Qact(tab1_nq(idim+1:idim+ndim),       &
                                          basis_set%tab_basisPrimSG(L,ib))
               tab_nq(idim+1:idim+ndim) = tab_nq(idim+1:idim+ndim) +    &
                                               tab1_nq(idim+1:idim+ndim)
            END DO
            idim = idim+ndim
          END DO
        CASE (2,4) ! Sparse basis (Smolyak 2d  implementation)
          idim = 0
          DO ib=1,basis_set%nb_basis
            ndim = basis_set%tab_basisPrimSG(0,ib)%ndim
            DO L=0,basis_set%L_SparseGrid
               CALL get_tab_nq_OF_Qact(tab1_nq(idim+1:idim+ndim),       &
                                          basis_set%tab_basisPrimSG(L,ib))
               tab_nq(idim+1:idim+ndim) = tab_nq(idim+1:idim+ndim) +    &
                                               tab1_nq(idim+1:idim+ndim)
            END DO
            idim = idim+ndim
          END DO

        CASE DEFAULT
          write(out_unitp,*) ' ERROR in',name_sub
          write(out_unitp,*) ' WRONG SparseGrid_type',basis_set%SparseGrid_type
          write(out_unitp,*) ' The possibilities are: 0, 1, 2, 4'
          STOP
        END SELECT
      END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------

      END SUBROUTINE get_tab_nq_OF_Qact

      integer FUNCTION get_nb_FROM_basis(basis_set)

      TYPE (basis), intent(in) :: basis_set

!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='get_nb_FROM_basis'
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL RecWriteMini_basis(basis_set)
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------

      IF (associated(basis_set%nDindB)) THEN
        get_nb_FROM_basis = basis_set%nDindB%max_nDI
      ELSE
        get_nb_FROM_basis = basis_set%nb
      END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------

      END FUNCTION get_nb_FROM_basis

      integer FUNCTION get_nb_bi_FROM_AllBasis(AllBasis)

      TYPE (param_AllBasis), intent(in) :: AllBasis

!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='get_nb_bi_FROM_AllBasis'
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------

      IF (associated(AllBasis%Basis2n)) THEN
        get_nb_bi_FROM_AllBasis = get_nb_FROM_basis(AllBasis%Basis2n)
      ELSE
        get_nb_bi_FROM_AllBasis = 1
      END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------

      END FUNCTION get_nb_bi_FROM_AllBasis

  FUNCTION get_nb_be_FROM_basis(basis_set)
  USE mod_system
  IMPLICIT NONE

  integer :: get_nb_be_FROM_basis
  !----- variables for the construction of H ---------------------------
  TYPE (basis), intent(in) :: basis_set

  integer :: ib
  !----- for debuging --------------------------------------------------
      !logical, parameter :: debug = .TRUE.
       logical, parameter :: debug = .FALSE.
  !-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING get_nb_be_FROM_basis'
  END IF

  get_nb_be_FROM_basis = -1
  IF (NewBasisEl) THEN
    IF (basis_set%nb_basis > 0) THEN
      DO ib=1,size(basis_set%tab_Pbasis)
        IF (basis_set%tab_Pbasis(ib)%Pbasis%type == 2) THEN
          get_nb_be_FROM_basis = basis_set%tab_Pbasis(ib)%Pbasis%nb
          EXIT
        END IF
      END DO
    ELSE
      IF (basis_set%type == 2) THEN
        get_nb_be_FROM_basis = basis_set%nb
      END IF
    END IF
  END IF

  IF (debug) THEN
    write(out_unitp,*) 'END get_nb_be_FROM_basis'
  END IF
  END FUNCTION get_nb_be_FROM_basis

END MODULE mod_basis_set_alloc
