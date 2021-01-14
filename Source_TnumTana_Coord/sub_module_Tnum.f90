!===========================================================================
!===========================================================================
!This file is part of Tnum-Tana.
!
!    Tnum-Tana is a free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    Tnum-Tana is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ElVibRot.  If not, see <http://www.gnu.org/licenses/>.
!
!    Copyright 2015  David Lauvergnat
!      Tnum is written David Lauvergnat [1]
!      Tana is written by Mamadou Ndong [1] and David Lauvergnat [1]
!         with contributions
!          Emil Lund klinting (coupling with MidasCpp) [3]'
!
![1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France
![3]: Department of Chemistry, Aarhus University, DK-8000 Aarhus C, Denmark
!
!===========================================================================
!===========================================================================
MODULE mod_Tnum
      USE mod_system
      USE mod_nDFit,            only: param_nDFit
      USE mod_QTransfo,         only: type_qtransfo, write_qtransfo,    &
                                      dealloc_qtransfo, dealloc_array,  &
                                      alloc_array, read_qtransfo,       &
                                      sub_type_name_of_qin,             &
                                      sub_check_lineartransfo,          &
                                      qtransfo1toqtransfo2
      USE mod_LinearNMTransfo,  only: type_nmtransfo, dealloc_array,    &
                                      alloc_array, read_lineartransfo,  &
                                      read_nmtransfo, alloc_lineartransfo
      USE mod_RPHTransfo,       only: type_rphtransfo, write_rphtransfo,&
                                      dealloc_array, alloc_array,       &
                                      rphtransfo1torphtransfo2,         &
                                      set_rphtransfo, dealloc_rphtransfo
      USE CurviRPH_mod,         only: curvirph_type, dealloc_curvirph,  &
                                      curvirph1_to_curvirph2, get_CurviRPH
      USE mod_ActiveTransfo
      USE mod_CartesianTransfo
      USE mod_Tana_Sum_OpnD

      IMPLICIT NONE

        PRIVATE

!-------------------------------------------------------------------------------
        TYPE param_PES_FromTnum

          integer           :: nb_elec              = 1       ! nb of electronic PES

          logical           :: deriv_WITH_FiniteDiff= .FALSE. ! IF true, force to use finite difference scheme

          real (kind=Rkind) :: pot0                 = ZERO    ! minimum of the PES
          real (kind=Rkind) :: min_pot              = HUGE(ONE) ! minimum of the PES
          real (kind=Rkind) :: max_pot              = -HUGE(ONE) ! minimum of the PES

          real (kind=Rkind) :: stepOp               = ONETENTH**2

          logical           :: calc_scalar_Op       = .FALSE.

          logical           :: opt                  = .FALSE. ! If we minimize the PES (not used yet)
          logical           :: pot_cplx             = .FALSE. ! complex PES
          logical           :: HarD                 = .TRUE.  ! Harmonic Domain for the PES
          integer           :: pot_itQtransfo       = -1      ! for new Qtransfo (default nb_QTransfo, dyn. coordinates)
          integer           :: nb_scalar_Op         = 0       ! nb of Operator
          integer           :: nb_CAP               = 0       ! nb of Operator

          ! parameters for on-the-fly calculations
          logical           :: OnTheFly             = .FALSE. ! On-the-fly calculation f PES and dipole
          logical           :: Read_OnTheFly_only   = .FALSE. ! Read On-the-fly calculation f PES and dipole only
                                                              ! without ab initio calculation
           integer :: charge       = 0
           integer :: multiplicity = -1

          character (len=Name_len)      :: ab_initio_prog     = ' '  ! gaussian, gamess

          character (len=Name_longlen)  :: ab_initio_methEne     = ' '
          character (len=Name_longlen)  :: ab_initio_methDip     = ' '
          character (len=Name_longlen)  :: ab_initio_basisEne    = ' '
          character (len=Name_longlen)  :: ab_initio_basisDip    = ' '


          character (len=Line_len)      :: commande_unix      = ' '
          logical                       :: header             = .FALSE.
          logical                       :: footer             = .FALSE.
          character (len=Name_len)      :: file_name_OTF      = ' '
          character (len=Name_len)      :: file_name_fchk     = ' '

          ! parameters for the Quantum Model Lib (ECAM), KEO+PES
          logical                       :: QMLib = .FALSE.

          ! parameters for the PES obtained from internal fit program
          logical                               :: nDfit_Op            = .FALSE.
          character (len=Line_len)              :: BaseName_nDfit_file = ' '
          character (len=Line_len)              :: nDFit_V_name_Fit    = ' '
          character (len=Line_len) ,allocatable :: nDFit_Scalar_Op_name_Fit(:)
          !TYPE (param_nDFit)              :: para_nDFit_V
          !TYPE (param_nDFit), allocatable :: para_nDFit_Scalar_Op(:)

        END TYPE param_PES_FromTnum
!-------------------------------------------------------------------------------
        TYPE CoordType
          logical            :: WriteCC    = .FALSE.

          real (kind=Rkind)  :: stepQ = ONETENTH**4
          logical            :: num_x = .FALSE.


          integer :: nb_act           = 0
          integer :: nb_var           = 0         ! nb_var= 3*nat-6 + nb_extra_Coord
          integer :: ndimG            = 0
          integer :: nb_extra_Coord   = 0
          integer :: ncart            = 0
          integer :: ncart_act        = 0
          integer :: nat0             = 0
          integer :: nat              = 0
          integer :: nat_act          = 0

          ! for the ab initio calculation
          integer                          :: charge       = 0
          integer                          :: multiplicity = -1
          integer                          :: nb_elec      = -1      ! here it is the number of electrons

          integer,                 pointer :: Z(:)         => null() ! true pointer
          character (len=Name_len),pointer :: name_Qdyn(:) => null() ! true pointer
          character (len=Name_len),pointer :: name_Qact(:) => null() ! true pointer
          character (len=Name_len),pointer :: symbole(:)   => null() ! true pointer

          logical                          :: cos_th       = .FALSE.  ! T => coordinate (valence angle) => cos(th)
                                                                      ! F => coordinate (valence angle) => th

          logical                          :: Without_Rot         = .FALSE.
          logical                          :: Centered_ON_CoM     = .TRUE.
          logical                          :: With_VecCOM         = .FALSE.

          logical                          :: Old_Qtransfo        = .FALSE.
          logical                          :: Cart_transfo        = .FALSE.
          logical                          :: Rot_Dip_with_EC     = .FALSE.



          integer                          :: nb_Qtransfo         = -1
          integer                          :: opt_param           =  0
          integer, allocatable             :: opt_Qdyn(:)

          integer                            :: itPrim               = -1
          integer                            :: itNM                 = -1
          integer                            :: itRPH                = -1
          TYPE (Type_ActiveTransfo), pointer, public :: ActiveTransfo=> null() ! it'll point on tab_Qtransfo
          TYPE (Type_NMTransfo),     pointer :: NMTransfo            => null() ! it'll point on tab_Qtransfo
          TYPE (Type_RPHTransfo),    pointer :: RPHTransfo           => null() ! it'll point on tab_Qtransfo
          TYPE (Type_Qtransfo),      pointer :: tab_Qtransfo(:)      => null()
          TYPE (Type_Qtransfo),      pointer :: tab_Cart_transfo(:)  => null()
          TYPE (Type_RPHTransfo),    pointer :: RPHTransfo_inact2n   => null() ! For the inactive coordinates (type 21)

          TYPE (CurviRPH_type)               :: CurviRPH

          integer, pointer :: liste_QactTOQsym(:) => null()   ! true pointer
          integer, pointer :: liste_QsymTOQact(:) => null()   ! true pointer
          integer, pointer :: liste_QactTOQdyn(:) => null()   ! true pointer
          integer, pointer :: liste_QdynTOQact(:) => null()   ! true pointer
          integer, pointer :: nrho_OF_Qact(:)     => null()   ! enables to define the volume element
          integer, pointer :: nrho_OF_Qdyn(:)     => null()   ! enables to define the volume element


          integer :: nb_act1   =0
          integer :: nb_inact2n=0,nb_inact21=0,nb_inact22=0
          integer :: nb_inact20=0,nb_inact=0
          integer :: nb_inact31=0
          integer :: nb_rigid0 =0,nb_rigid100=0,nb_rigid=0


          real (kind=Rkind), pointer :: masses(:) => null()        ! true pointer
          integer,           pointer :: active_masses(:) => null() ! for partial hessian (PVSCF)

          real (kind=Rkind), pointer :: d0sm(:) => null()
          real (kind=Rkind)          :: Mtot = ZERO
          real (kind=Rkind)          :: Mtot_inv = ZERO
        CONTAINS
          PROCEDURE, PRIVATE, PASS(mole1) :: CoordType2_TO_CoordType1
          GENERIC,   PUBLIC  :: assignment(=) => CoordType2_TO_CoordType1
        END TYPE CoordType
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!       The extended type is defined to be able to use the TYPE zmatrix,
!       which is not a well-chosen name (it was for the first Tnum versions).
!-------------------------------------------------------------------------------
        TYPE, EXTENDS (CoordType) :: zmatrix
        END TYPE zmatrix
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
        !!@description: TODO
        !!@param: TODO
!-------------------------------------------------------------------------------
        TYPE Tnum
          real (kind=Rkind)          :: stepT   = ONETENTH**4
          logical                    :: num_GG  = .FALSE.
          logical                    :: num_g   = .FALSE.
          logical                    :: num_x   = .FALSE.

          logical                    :: Gdiago  = .FALSE.
          logical                    :: Gcte    = .FALSE.
          integer                    :: NonGcteRange(2) = 0
          logical                    :: Inertia = .TRUE.
          real (kind=Rkind), pointer :: Gref(:,:) => null() ! reference value if Gcte=.true.

          integer                    :: nrho              = 2
          integer                    :: vep_type          = -1 ! The default depends on the GG metric
                                                               ! 0   : without vep
                                                               ! 1   : normal vep
                                                               ! 100 : normal vep for the full metric tensor
                                                               !         (not the reduced one), for coord type=100
                                                               ! -100: vep type100 with a bug on dng (do not use)
          integer                    :: JJ                = 0
          logical                    :: WriteT            = .FALSE.
          logical                    :: With_Cart_Transfo = .TRUE.

          logical                    :: Tana              = .FALSE.
          logical                    :: MidasCppForm      = .FALSE.
          logical                    :: MCTDHForm         = .FALSE.
          logical                    :: LaTeXForm         = .FALSE.
          logical                    :: VSCFForm          = .FALSE.

          logical                    :: f2f1_ana          = .FALSE.

          integer                    :: KEO_TalyorOFQinact2n = -1

          TYPE (param_PES_FromTnum)  :: para_PES_FromTnum

          TYPE (sum_opnd)            :: TWOxKEO,ExpandTWOxKEO
        CONTAINS
          PROCEDURE, PRIVATE, PASS(Tnum1) :: Tnum2_TO_Tnum1
          GENERIC,   PUBLIC  :: assignment(=) => Tnum2_TO_Tnum1
        END TYPE Tnum
!-------------------------------------------------------------------------------

      PUBLIC :: param_PES_FromTnum, Tnum, dealloc_Tnum

      PUBLIC :: Write_f2f1vep, Write_TcorTrot

      PUBLIC :: zmatrix
      PUBLIC :: Write_mole, Read_mole, dealloc_zmat, mole1TOmole2

      PUBLIC :: CoordType
      PUBLIC :: Write_CoordType, Read_CoordType,                        &
                               dealloc_CoordType, CoordType1TOCoordType2
      PUBLIC :: check_charge, type_var_analysis_OF_CoordType
      PUBLIC :: Sub_paraRPH_TO_CoordType, Sub_CoordType_TO_paraRPH,     &
                                            Sub_CoordType_TO_paraRPH_new
      PUBLIC :: CoordTypeRPH_TO_CoordTypeFlex
      PUBLIC :: Set_OptimizationPara_FROM_CoordType

      PUBLIC :: get_CurviRPH



      CONTAINS
!================================================================
! Transfert data from the 1st Qtransfo to the CoordType
!================================================================
  SUBROUTINE Set_masses_Z_TO_CoordType(mole,Qtransfo)
        TYPE (CoordType),     intent(inout) :: mole
        TYPE (Type_Qtransfo), intent(in)    :: Qtransfo

        character (len=Name_len) :: name_transfo

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory,err_io
      character (len=*), parameter :: name_sub = "Set_masses_Z_TO_CoordType"
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      !-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
      END IF
      !-----------------------------------------------------------


        name_transfo = Qtransfo%name_transfo
        CALL string_uppercase_TO_lowercase(name_transfo)

        SELECT CASE (name_transfo)

        CASE ('zmat') ! It should be one of the first transfo read

          mole%nat0      = Qtransfo%ZmatTransfo%nat0
          mole%nat       = Qtransfo%ZmatTransfo%nat

          mole%nb_var    = Qtransfo%ZmatTransfo%nb_var

          mole%ncart     = Qtransfo%ZmatTransfo%ncart
          mole%nat_act   = Qtransfo%ZmatTransfo%nat_act
          mole%ncart_act = Qtransfo%ZmatTransfo%ncart_act


          mole%Z       => Qtransfo%ZmatTransfo%Z
          mole%symbole => Qtransfo%ZmatTransfo%symbole
          mole%masses  => Qtransfo%ZmatTransfo%masses


        CASE ('bunch','bunch_poly') ! It should one of the first transfo

          mole%nat0      = Qtransfo%BunchTransfo%nat0
          mole%nat       = Qtransfo%BunchTransfo%nat
          mole%nb_var    = Qtransfo%BunchTransfo%nb_var
          mole%ncart     = Qtransfo%BunchTransfo%ncart
          mole%nat_act   = Qtransfo%BunchTransfo%nat_act
          mole%ncart_act = Qtransfo%BunchTransfo%ncart_act

          mole%Z       => Qtransfo%BunchTransfo%Z
          mole%symbole => Qtransfo%BunchTransfo%symbole
          mole%masses  => Qtransfo%BunchTransfo%masses

        CASE ('qtox_ana')

          mole%nat0      = Qtransfo%QTOXanaTransfo%nat0
          mole%nat       = Qtransfo%QTOXanaTransfo%nat
          mole%nb_var    = Qtransfo%QTOXanaTransfo%nb_var
          mole%ncart     = Qtransfo%QTOXanaTransfo%ncart
          mole%nat_act   = Qtransfo%QTOXanaTransfo%nat_act
          mole%ncart_act = Qtransfo%QTOXanaTransfo%ncart_act

          mole%Z       => Qtransfo%QTOXanaTransfo%Z
          mole%symbole => Qtransfo%QTOXanaTransfo%symbole
          mole%masses  => Qtransfo%QTOXanaTransfo%masses

        END SELECT

      !-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
      !-----------------------------------------------------------

  END SUBROUTINE Set_masses_Z_TO_CoordType

!================================================================
!       Write the "mole" (old with zmatrix type)
!       Write the CoordType
!================================================================
  SUBROUTINE Write_mole(mole,print_mole)
  TYPE (zmatrix), intent(in)            :: mole
  logical,        intent(in), optional  :: print_mole

  IF (present(print_mole)) THEN
    CALL  Write_CoordType(mole%CoordType,print_mole)
  ELSE
    CALL  Write_CoordType(mole%CoordType)
  END IF

  END SUBROUTINE Write_mole

  SUBROUTINE Write_CoordType(mole,print_mole)
    TYPE (CoordType), intent(in)             :: mole
    logical,          intent(in), optional   :: print_mole

      integer i,j,ni,it
      logical :: print_mole_loc


      IF (.NOT. present(print_mole)) THEN
        IF (print_CoordType_done .AND. print_level <= 1) RETURN
      END IF
      print_CoordType_done = .TRUE. ! from sub_module_system

        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING Write_CoordType'
        write(out_unitp,*)
        write(out_unitp,*)
        write(out_unitp,*) 'Without_Rot:      ',mole%Without_Rot
        write(out_unitp,*) 'Centered_ON_CoM   ',mole%Centered_ON_CoM
        write(out_unitp,*) 'Cart_transfo:     ',mole%Cart_transfo
        write(out_unitp,*) 'stepQ,num_x:      ',mole%stepQ,mole%num_x
        write(out_unitp,*) 'Parameter(s) to be optimized?:     ',mole%opt_param
        write(out_unitp,*) 'Coordinates Qdyn to be optimized?: ',mole%opt_Qdyn(:)

!       -------------------------------------------------------
!       -------------------------------------------------------
        write(out_unitp,*)
        write(out_unitp,*) '------------------------------------------'
        write(out_unitp,*) '--- Number of variables ------------------'
        write(out_unitp,*)
        write(out_unitp,*) 'nb_act,nb_var',mole%nb_act,mole%nb_var
        write(out_unitp,*) 'nb_extra_Coord',mole%nb_extra_Coord
        write(out_unitp,*) 'Number of variables of type...'
        write(out_unitp,*) 'type  1: active                       (nb_act1):',  &
                    mole%nb_act1
        write(out_unitp,*) 'type 31: for gaussian WP           (nb_inact31):',  &
                    mole%nb_inact31
        write(out_unitp,*) 'type 21: inactive harmonic         (nb_inact21):',  &
                    mole%nb_inact21
        write(out_unitp,*) 'type 22: inactive non harmonic     (nb_inact22):',  &
                    mole%nb_inact22
        write(out_unitp,*) 'type 20: adiabatically constrained (nb_inact20):',  &
                    mole%nb_inact20
        write(out_unitp,*) 'type  0: rigid                      (nb_rigid0):',  &
                    mole%nb_rigid0
        write(out_unitp,*) 'type100: rigid                    (nb_rigid100):',  &
                    mole%nb_rigid100
        write(out_unitp,*)
        write(out_unitp,*) 'nb_inact2n = nb_inact21 + nb_inact22  :',   &
                    mole%nb_inact2n
        write(out_unitp,*) 'nb_act= nb_act1+nb_inact31+nb_inact2n :',   &
                    mole%nb_act
        write(out_unitp,*) 'nb_inact   = nb_inact20 + nb_inact2n  :',   &
                    mole%nb_inact
        write(out_unitp,*) 'nb_rigid   = nb_rigid0 + nb_rigid100  :',   &
                    mole%nb_rigid
        write(out_unitp,*) '------------------------------------------'
        write(out_unitp,*) 'ndimG = nb_act+6 + nb_act+3 or nb_act :',   &
                    mole%ndimG
        write(out_unitp,*) '------------------------------------------'
        write(out_unitp,*)
        write(out_unitp,*)
        write(out_unitp,*) '------------------------------------------'
        write(out_unitp,*) '------------------------------------------'
        write(out_unitp,*) 'tab_Qtransfo'
        write(out_unitp,*) '  itPrim: ',mole%itPrim
        write(out_unitp,*) '  itNM:   ',mole%itNM
        write(out_unitp,*) '  itRPH:  ',mole%itRPH

        write(out_unitp,*) '------------------------------------------'
        DO it=1,mole%nb_Qtransfo
          CALL Write_Qtransfo(mole%tab_Qtransfo(it))
        END DO
        write(out_unitp,*) '------------------------------------------'

        IF (associated(mole%tab_Cart_transfo)) THEN
          write(out_unitp,*) '------------------------------------------'
          write(out_unitp,*) '------------------------------------------'
          write(out_unitp,*) 'tab_Cart_transfo'
          write(out_unitp,*) '------------------------------------------'
          DO it=1,size(mole%tab_Cart_transfo)
            CALL Write_Qtransfo(mole%tab_Cart_transfo(it))
          END DO
          write(out_unitp,*) '------------------------------------------'
        END IF

!       -------------------------------------------------------
        write(out_unitp,*) 'ActiveTransfo from mole'
        CALL Write_ActiveTransfo(mole%ActiveTransfo)

        write(out_unitp,*) 'nrho_OF_Qdyn',mole%nrho_OF_Qdyn
        write(out_unitp,*) 'nrho_OF_Qact',mole%nrho_OF_Qact

        write(out_unitp,*) '------------------------------------------'
        write(out_unitp,*) '------------------------------------------'
        write(out_unitp,*)
!       -------------------------------------------------------

        IF (associated(mole%RPHTransfo_inact2n)) THEN
          write(out_unitp,*) '------------------------------------------'
          write(out_unitp,*) '------------------------------------------'
          write(out_unitp,*) '----------- RPHTransfo_inact2n -----------'

          CALL Write_RPHTransfo(mole%RPHTransfo_inact2n)

          write(out_unitp,*) '------------------------------------------'
          write(out_unitp,*) '------------------------------------------'
        END IF


!       -------------------------------------------------------
!       Mtot_inv and sqrt(masses(i)) masses(i)
!       -------------------------------------------------------
        write(out_unitp,*)
        write(out_unitp,*) '------------------------------------------'
        write(out_unitp,*) '--- Masses ... ---------------------------'
        write(out_unitp,*)
        write(out_unitp,*) 'ncart,ncart_act',mole%ncart,mole%ncart_act
        write(out_unitp,*) 'ATOM  mass    sqrt(mass) ...'
        DO i=1,mole%ncart_act,3

          write(out_unitp,91) i,mole%masses(i),mole%d0sm(i)
 91       format('at: ',i3,f18.6,f15.6)

        END DO
        IF (mole%Mtot_inv /= ZERO)                                      &
             write(out_unitp,*) 'Mtot_inv,Mtot',mole%Mtot_inv,ONE/mole%Mtot_inv
        write(out_unitp,*) '------------------------------------------'
        write(out_unitp,*)
!       -------------------------------------------------------
!       -------------------------------------------------------

        write(out_unitp,*) 'END Write_CoordType'


  END SUBROUTINE Write_CoordType


!================================================================
!
!     dealloc zmatrix type or CoordType
!
!================================================================
  SUBROUTINE dealloc_zmat(mole)
  TYPE (zmatrix), intent(inout) :: mole

    CALL dealloc_CoordType(mole%CoordType)

  END SUBROUTINE dealloc_zmat

  SUBROUTINE dealloc_CoordType(mole)
  TYPE (CoordType), intent(inout) :: mole

       integer        :: it

      character (len=*), parameter :: name_sub='dealloc_CoordType'


      !write(out_unitp,*) 'BEGINNING ',name_sub
      !CALL flush_perso(out_unitp)

        mole%WriteCC      = .FALSE.

        mole%stepQ        = ONETENTH**4
        mole%num_x        = .FALSE.

        mole%nb_act       = 0
        mole%nb_var       = 0
        mole%nb_extra_Coord = 0

        mole%ndimG        = 0
        mole%ncart        = 0
        mole%ncart_act    = 0
        mole%nat          = 0
        mole%nat0         = 0
        mole%nat_act      = 0


        mole%charge       = 0
        mole%multiplicity = -1
        mole%nb_elec      = -1

        nullify(mole%name_Qdyn) ! true pointer
        nullify(mole%name_Qact) ! true pointer
        nullify(mole%Z)         ! true pointer
        nullify(mole%symbole)   ! true pointer
        nullify(mole%masses)    ! true pointer

        mole%Without_Rot     = .FALSE.
        mole%Centered_ON_CoM = .TRUE.
        mole%cos_th          = .FALSE.

        mole%Old_Qtransfo    = .FALSE.
        mole%Cart_transfo    = .FALSE.

        mole%nb_Qtransfo     = -1
        mole%itNM            = -1
        mole%itRPH           = -1
        mole%itPrim          = -1
        mole%opt_param       = 0

        IF (allocated(mole%opt_Qdyn)) THEN
          CALL dealloc_NParray(mole%opt_Qdyn,"mole%opt_Qdyn",name_sub)
        END IF

        nullify(mole%ActiveTransfo)  ! true pointer
        nullify(mole%NMTransfo)      ! true pointer
        nullify(mole%RPHTransfo)     ! true pointer

        CALL dealloc_CurviRPH(mole%CurviRPH)

        IF (associated(mole%tab_Qtransfo)) THEN
          DO it=lbound(mole%tab_Qtransfo,dim=1),ubound(mole%tab_Qtransfo,dim=1)
            CALL dealloc_Qtransfo(mole%tab_Qtransfo(it))
          END DO
          CALL dealloc_array(mole%tab_Qtransfo,                         &
                            "mole%tab_Qtransfo",name_sub)
        END IF

        IF (associated(mole%tab_Cart_transfo)) THEN
          DO it=lbound(mole%tab_Cart_transfo,dim=1),ubound(mole%tab_Cart_transfo,dim=1)
            CALL dealloc_Qtransfo(mole%tab_Cart_transfo(it))
          END DO
          CALL dealloc_array(mole%tab_Cart_transfo,                     &
                            "mole%tab_Cart_transfo",name_sub)
        END IF

      IF (associated(mole%RPHTransfo_inact2n)) THEN
        CALL dealloc_array(mole%RPHTransfo_inact2n,                     &
                          'mole%RPHTransfo_inact2n',name_sub)
      END IF

        nullify(mole%liste_QactTOQsym)   ! true pointer
        nullify(mole%liste_QsymTOQact)   ! true pointer
        nullify(mole%liste_QactTOQdyn)   ! true pointer
        nullify(mole%liste_QdynTOQact)   ! true pointer

        IF (associated(mole%nrho_OF_Qact))  THEN
          CALL dealloc_array(mole%nrho_OF_Qact,                         &
                            "mole%nrho_OF_Qact",name_sub)
        END IF
        IF (associated(mole%nrho_OF_Qdyn))  THEN
          CALL dealloc_array(mole%nrho_OF_Qdyn,                         &
                            "mole%nrho_OF_Qdyn",name_sub)
        END IF

        mole%nb_act1     = 0
        mole%nb_inact2n  = 0
        mole%nb_inact21  = 0
        mole%nb_inact22  = 0
        mole%nb_inact20  = 0
        mole%nb_inact    = 0
        mole%nb_inact31  = 0
        mole%nb_rigid0   = 0
        mole%nb_rigid100 = 0
        mole%nb_rigid    = 0

        IF (associated(mole%active_masses))  THEN
          CALL dealloc_array(mole%active_masses,                        &
                            "mole%active_masses",name_sub)
        END IF

        IF (associated(mole%d0sm))  THEN
          CALL dealloc_array(mole%d0sm,"mole%d0sm",name_sub)
        END IF
        mole%Mtot         = ZERO
        mole%Mtot_inv     = ZERO

      !write(out_unitp,*) 'END ',name_sub

  END SUBROUTINE dealloc_CoordType
!===============================================================================
!     read zmatrix type or CoordType parameters for Tnum
!===============================================================================
  !!@description: read zmatrix type parameters for Tnum
  SUBROUTINE Read_mole(mole,para_Tnum,const_phys)
  USE mod_constant
  IMPLICIT NONE

  TYPE (zmatrix),   intent(inout)   :: mole
  TYPE (Tnum),      intent(inout)   :: para_Tnum
  TYPE (constant),  intent(inout)   :: const_phys

    CALL Read_CoordType(mole%CoordType,para_Tnum,const_phys)

  END SUBROUTINE Read_mole

  !!@description: read CoordType parameters for Tnum
  !!@param: TODO
  SUBROUTINE Read_CoordType(mole,para_Tnum,const_phys)
      USE mod_ActiveTransfo,    only : Read_ActiveTransfo
      USE mod_ZmatTransfo,      only : Read_ZmatTransfo
      USE mod_CartesianTransfo, only : Write_CartesianTransfo
      USE mod_constant
      USE mod_MPI
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType), intent(inout)   :: mole
      TYPE (Tnum),      intent(inout)   :: para_Tnum
      TYPE (constant),  intent(inout)   :: const_phys


!----- physical and mathematical constants ----------------------------
!      - for the definition of d0d1d2d3_Qeq -------------
      logical :: Qeq_sym,Q0_sym

!-------------------------------------------------------
      logical :: WriteCC ! write cartesian coordinates
      logical :: WriteT  ! write T, GG, g

!-------------------------------------------------------
!     - for the CoordType --------------------------------
      logical :: zmat,bunch,cart,cos_th

      integer :: nb_act,nb_extra_Coord
      integer :: nat,nb_var

!     - for Q_transfo ----------------------------------
      logical :: Old_Qtransfo,Cart_transfo,Rot_Dip_with_EC
      integer :: nb_Qtransfo,nb_Qin
      character (len=Name_len) :: name_transfo

!     - for NM_TO_sym ----------------------------------
      logical :: NM,NM_TO_sym,hessian_old,hessian_cart,hessian_onthefly
      logical :: purify_hess,k_Half
      character (len=Name_len)      :: file_hessian

!     - for rotational constraints ---------------------
!     Without_Rot = .TRUE. => calculation without rotation (constraints on the rotation)
      logical :: Without_Rot
!     Centered_ON_CoM = .TRUE. => the molecule is recentered on the center of mass
      logical :: Centered_ON_CoM
      integer :: JJ  ! JJ=0, the rotation is not printed. JJ>0 the rotation is printed

!     - for new orientation (for zmat only) ---------------------
      logical :: New_Orient
      real (kind=Rkind) :: vAt1(3),vAt2(3),vAt3(3)


!     - for the symmetry -------------------------------
!     sym : .TRUE. if we use the symmetry
      logical :: sym,check_sym

!     - for the molecule -------------------------------
      integer :: charge,multiplicity
      logical :: header,footer
      character (len=Name_len)      :: file_name_OTF,file_name_fchk
      character (len=Line_len)      :: commande_unix
      character (len=Name_longlen)  :: ab_initio_meth,ab_initio_basis
      character (len=Name_longlen)  :: ab_initio_methEne,ab_initio_basisEne
      character (len=Name_longlen)  :: ab_initio_methDip,ab_initio_basisdip

      character (len=Name_len)      :: ab_initio_prog

!     - for the files -----------------------------------------------
      integer :: nio


!     - for Tnum or Tana ----------------------------------------------
      integer           :: nrho,vep_type,NonGcteRange(2)
      logical           :: num_GG,num_g,num_x,Gdiago,Gcte,With_VecCOM
      logical           :: With_Tab_dnQflex
      logical           :: Tana,MidasCppForm,MCTDHForm,LaTeXForm,VSCFForm,f2f1_ana
      real (kind=Rkind) :: stepT,stepOp
      integer           :: KEO_TalyorOFQinact2n ! taylor epxansion along coordinate 2n (21) types
!     - end for the CoordType ----------------------------


!     - divers ------------------------------------------
      integer :: i,n,it,iat,i_Q,iQout,iQin

      NAMELIST /variables/nat,zmat,bunch,cos_th,                        &
                     Without_Rot,Centered_ON_CoM,JJ,                    &
                     New_Orient,vAt1,vAt2,vAt3,With_VecCOM,             &
                     nb_var,nb_act,With_Tab_dnQflex,nb_extra_Coord,     &
                     Old_Qtransfo,nb_Qtransfo,Cart_transfo,             &
                     Rot_Dip_with_EC,sym,check_sym,                     &
                     NM,NM_TO_sym,hessian_old,purify_hess,k_Half,       &
                     hessian_cart,hessian_onthefly,file_hessian,stepOp, &
                     stepT,num_GG,num_g,num_x,nrho,Tana,vep_type,       &
                     Gdiago,Gcte,NonGcteRange,                          &
                     MidasCppForm,MCTDHForm,LaTeXForm,VSCFForm,         &
                     KEO_TalyorOFQinact2n,f2f1_ana,                     &
                     charge,multiplicity,                               &
                     header,footer,file_name_OTF,file_name_fchk,        &
                     commande_unix,ab_initio_prog,                      &
                     ab_initio_meth,ab_initio_basis,                    &
                     ab_initio_methEne,ab_initio_basisEne,              &
                     ab_initio_methDip,ab_initio_basisDip,              &
                     WriteCC,WriteT

!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Read_CoordType'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      !IF (debug) THEN
        write(out_unitp,*) '================================================='
        write(out_unitp,*) '================================================='
        write(out_unitp,*) 'BEGINNING ',name_sub
      !END IF

!-------------------------------------------------------
!     initializations
!-------------------------------------------------------
      IF (.NOT. const_phys%constant_done) THEN
        CALL sub_constantes(const_phys,Read_Namelist=.FALSE.)
      END IF


      stepOp               = ZERO
      stepT                = ONETENTH**4
      num_GG                = .FALSE.
      num_g                = .FALSE.
      num_x                = .FALSE.
      Gdiago               = .FALSE.
      Gcte                 = .FALSE.
      NonGcteRange(:)      = 0
      Tana                 = .FALSE.
      MidasCppForm         = .FALSE.
      MCTDHForm            = .FALSE.
      LaTeXForm            = .FALSE.
      VSCFForm             = .FALSE.

      nrho                 = 1
      vep_type             = -1
      KEO_TalyorOFQinact2n = -1
      f2f1_ana             = .FALSE.

      With_VecCOM          = .FALSE.
      Without_Rot          = .FALSE.
      Centered_ON_CoM      = .TRUE.
      JJ                   = 0
      New_Orient           = .FALSE.
      vAt1(:)              = ZERO
      vAt2(:)              = ZERO
      vAt3(:)              = ZERO

      zmat                 = .TRUE.
      bunch                = .FALSE.
      cos_th               = .TRUE.
      cart                 = .FALSE.

      Cart_transfo         = .FALSE.
      Rot_Dip_with_EC      = .FALSE.
      Old_Qtransfo         = .TRUE.
      nb_Qtransfo          = -1

      With_Tab_dnQflex     = .FALSE.

      NM                   = .FALSE.
      NM_TO_sym            = .FALSE.
      purify_hess          = .FALSE.
      k_Half               = .FALSE.
      file_hessian         = 'xx_freq.fchk'
      hessian_old          = .TRUE.
      hessian_cart         = .TRUE.
      hessian_onthefly     = .FALSE.

      WriteCC              = .FALSE.
      WriteT               = .FALSE.

      nat                  = 0
      nb_var               = 0
      nb_act               = 0
      nb_extra_Coord       = 0
      sym                  = .FALSE.
      check_sym            = .TRUE.
      charge               = 0
      multiplicity         = -1
      header               = .FALSE.
      footer               = .FALSE.
      file_name_fchk       = ''
      file_name_OTF        = 'xx'
      ab_initio_meth       = ' hf '
      ab_initio_basis      = ' sto-3g '
      ab_initio_methEne    = ''
      ab_initio_basisEne   = ''
      ab_initio_methDip    = ''
      ab_initio_basisDip   = ''
      ab_initio_prog       = 'g03'
      commande_unix        = 'g03.run xx >err'

      read(in_unitp,variables,IOSTAT=err_read)
      IF (err_read < 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' End-of-file or End-of-record'
        write(out_unitp,*) ' The namelist "variables" is probably absent'
        write(out_unitp,*) ' check your data!'
        write(out_unitp,*) ' ERROR in ',name_sub
        STOP
      ELSE IF (err_read > 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' Some parameter name of the namelist "variables" are probaly wrong'
        write(out_unitp,*) ' check your data!'
        write(out_unitp,variables)
        write(out_unitp,*) ' ERROR in ',name_sub
        STOP
      END IF
      IF (print_level > 1) write(out_unitp,variables)

      IF (stepT  == ZERO)  stepT  = ONETENTH**4
      IF (stepOp == ZERO)  stepOp = stepT
      IF (Rot_Dip_with_EC) Cart_transfo = .TRUE.

      para_Tnum%JJ                   = JJ
      para_Tnum%With_Cart_Transfo    = (JJ>0) .AND. Cart_transfo

      para_Tnum%stepT                = stepT
      para_Tnum%num_x                = num_x
      para_Tnum%num_GG               = num_GG
      para_Tnum%num_g                = num_g

      para_Tnum%WriteT               = WriteT
      para_Tnum%Gdiago               = Gdiago
      para_Tnum%Gcte                 = Gcte

      IF ((nrho == 20 .OR. nrho == 10) .AND. vep_type /= -1 .AND. vep_type /= 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' nrho    ',nrho
        write(out_unitp,*) ' vep_type',vep_type
        write(out_unitp,*) ' nrho = 20 or 10                       => old way to force the vep to zero.'
        write(out_unitp,*) ' vep_type is defined (/= -1) and /= 0  => vep is not forced to zero'
        write(out_unitp,*) ' You have two possibilities:'
        write(out_unitp,*) '   keep the nrho values (20 or 10)  and change vep_type to 0 or -1'
        write(out_unitp,*) '   keep the vep_type values         and change nrho to 1 or 2'
        write(out_unitp,*) ' check your data!'
        write(out_unitp,variables)
        write(out_unitp,*) ' ERROR in ',name_sub
        STOP ' ERROR in Read_CoordType: inconsistent nrho and vep_types values.'
      END IF
      IF (nrho == 0)                  vep_type = 0 ! no vep for euclidean normalization
      IF (nrho == 10 .OR. nrho == 20) vep_type = 0 ! old way to force the vep to zero
      IF (Gcte .AND. vep_type == -1)  vep_type = 0 ! it is not always the case [dTau=sin(th)dth ]
      IF (vep_type == -1)             vep_type = 1
      para_Tnum%vep_type            = vep_type
      IF (nrho == 10)                 nrho = 1
      IF (nrho == 20)                 nrho = 2
      para_Tnum%nrho                = nrho

      IF (para_Tnum%vep_type == -100) THEN
        write(out_unitp,*) '======================================================'
        write(out_unitp,*) '======================================================'
        write(out_unitp,*) '======================================================'
        write(out_unitp,*) '======================================================'
        write(out_unitp,*) '======================================================'
        write(out_unitp,*) '======================================================'
        write(out_unitp,*) ' WARNING in ',name_sub
        write(out_unitp,*) ' vep_type',para_Tnum%vep_type
        write(out_unitp,*)
        write(out_unitp,*) '  You are using vep_type to recover the vep values with a bug.'
        write(out_unitp,*) '  Relevent, only when the coordinate types are 100.'
        write(out_unitp,*)
        write(out_unitp,*) '  DO NOT USE IT!!'
        write(out_unitp,*) ' WARNING in ',name_sub
        write(out_unitp,*) '======================================================'
        write(out_unitp,*) '======================================================'
        write(out_unitp,*) '======================================================'
        write(out_unitp,*) '======================================================'
        write(out_unitp,*) '======================================================'
        write(out_unitp,*) '======================================================'
      END IF

      para_Tnum%Tana                 = Tana

      para_Tnum%LaTeXForm            = (Tana .AND. LaTeXForm)
      para_Tnum%MCTDHForm            = (Tana .AND. MCTDHForm)
      para_Tnum%VSCFForm             = (Tana .AND. VSCFForm)
      para_Tnum%MidasCppForm         = (Tana .AND. MidasCppForm)

      para_Tnum%KEO_TalyorOFQinact2n = KEO_TalyorOFQinact2n
      para_Tnum%f2f1_ana             = f2f1_ana

      mole%stepQ        = stepT
      mole%num_x        = num_x


      mole%Old_Qtransfo    = Old_Qtransfo
      mole%nb_Qtransfo     = nb_Qtransfo
      mole%Cart_transfo    = Cart_transfo
      mole%Rot_Dip_with_EC = Rot_Dip_with_EC
      mole%Without_Rot     = Without_Rot
      mole%Centered_ON_CoM = Centered_ON_CoM
      mole%With_VecCOM     = With_VecCOM

      mole%nb_extra_Coord  = nb_extra_Coord

!=======================================================================
!===== New Coordinate transformations ==================================
!=======================================================================
      IF (mole%nb_Qtransfo > 1) THEN
        IF(MPI_id==0) THEN
          write(out_unitp,*) '================================================='
          write(out_unitp,*) '================================================='
          write(out_unitp,*) 'New Coordinate transformations',mole%nb_Qtransfo
          write(out_unitp,*) '================================================='
        ENDIF
        CALL alloc_array(mole%tab_Qtransfo,(/mole%nb_Qtransfo/),        &
                        "mole%tab_Qtransfo",name_sub)
        nb_Qin          = 0
        mole%opt_param  = 0
        DO it=1,nb_Qtransfo
          IF (debug) write(out_unitp,*) 'New Qtransfo',it

          mole%tab_Qtransfo(it)%num_transfo = it

          IF (it > 0) nb_Qin = mole%nb_var

          IF (it > 1) THEN ! not cartessian coordinates (Qout)
            mole%tab_Qtransfo(it)%type_Qout => mole%tab_Qtransfo(it-1)%type_Qin
            mole%tab_Qtransfo(it)%name_Qout => mole%tab_Qtransfo(it-1)%name_Qin
            mole%tab_Qtransfo(it)%nb_Qout    = mole%tab_Qtransfo(it-1)%nb_Qin
            IF (mole%tab_Qtransfo(it)%nb_Qout < 1) THEN
              write(out_unitp,*) ' ERROR in ',name_sub
              write(out_unitp,*) '  it:',it
              write(out_unitp,*) '  nb_Qout < 1 !!',mole%tab_Qtransfo(it)%nb_Qout
              write(out_unitp,*) ' Check the fortran !!'
              STOP
            END IF
          END IF

          CALL flush_perso(out_unitp)
          IF (debug) write(out_unitp,*) 'read Qtransfo',it

          CALL read_Qtransfo(mole%tab_Qtransfo(it),nb_Qin,mole%nb_extra_Coord,  &
                             const_phys%mendeleev)
          mole%tab_Qtransfo(it)%BeforeActive = (it == nb_Qtransfo-1)

          mole%opt_param = mole%opt_param + mole%tab_Qtransfo(it)%opt_param

          IF (mole%tab_Qtransfo(it)%Primitive_Coord) mole%itPrim = it

          name_transfo = mole%tab_Qtransfo(it)%name_transfo
          CALL string_uppercase_TO_lowercase(name_transfo)

          CALL Set_masses_Z_TO_CoordType(mole,mole%tab_Qtransfo(it))

          SELECT CASE (name_transfo)
          CASE ('bunch','bunch_poly')
            ! because we need BunchTransfo for Poly transfo
            mole%tab_Qtransfo(it+1)%BunchTransfo => mole%tab_Qtransfo(it)%BunchTransfo

          CASE ("nm")
            mole%itNM  = it
            IF (associated(mole%NMTransfo)) THEN
              write(out_unitp,*) ' ERROR in ',name_sub
              write(out_unitp,*) '  TWO NM transformations have been read'
              write(out_unitp,*) '  Only one is possible'
              write(out_unitp,*) ' Check your data !!'
              STOP
            ELSE
              mole%NMTransfo => mole%tab_Qtransfo(it)%NMTransfo
            END IF

          CASE ("rph")
            IF (it /= nb_Qtransfo-1) THEN
               write(out_unitp,*) ' WARNNING in ',name_sub
               write(out_unitp,*) ' The RPH (Reaction Path Hamiltonian) transfortmation MUST be just before "active".'
               write(out_unitp,*) 'it,name_transfo: ',it,mole%tab_Qtransfo(it)%name_transfo
               !STOP
            END IF
            IF (mole%itRPH /= -1) THEN
              write(out_unitp,*) ' ERROR in ',name_sub
              write(out_unitp,*) '  TWO RPH transformations have been read'
              write(out_unitp,*) '  Only one is possible'
              write(out_unitp,*) ' Check your data !!'
              STOP
            ELSE
              mole%itRPH      = it
              mole%RPHTransfo => mole%tab_Qtransfo(it)%RPHTransfo
            END IF

          CASE ('active')
            IF (it /= nb_Qtransfo) THEN
               write(out_unitp,*) ' ERROR in ',name_sub
               write(out_unitp,*) ' The active transfortmation MUST be the last one'
               write(out_unitp,*) 'it,name_transfo: ',it,mole%tab_Qtransfo(it)%name_transfo
               STOP
            END IF
            IF (associated(mole%ActiveTransfo)) THEN
              write(out_unitp,*) ' ERROR in ',name_sub
              write(out_unitp,*) '  TWO active transformations have been read'
              write(out_unitp,*) '  Only one is possible'
              write(out_unitp,*) ' Check your data !!'
              STOP
            ELSE
              mole%ActiveTransfo => mole%tab_Qtransfo(it)%ActiveTransfo
            END IF

            mole%ActiveTransfo%With_Tab_dnQflex = With_Tab_dnQflex


          CASE default
            CONTINUE
          END SELECT

          IF (debug) write(out_unitp,*) 'END: New Qtransfo',it
          CALL flush_perso(out_unitp)
        END DO
        write(out_unitp,*) 'Parameter(s) to be optimized?: ',mole%opt_param

        !=======================================================================
        ! analyzis of the transformations:
        name_transfo = mole%tab_Qtransfo(1)%name_transfo
        CALL string_uppercase_TO_lowercase(name_transfo)
        IF (name_transfo /= 'zmat'  .AND. name_transfo /= 'bunch' .AND. &
            name_transfo /= 'bunch_poly' .AND.                          &
            name_transfo /= 'qtox_ana') THEN
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) ' The first transfortmation MUST be ',    &
                      '"zmat" or "bunch" or "bunch_poly" or "QTOX_ana".'
           write(out_unitp,*) 'name_transfo: ',mole%tab_Qtransfo(1)%name_transfo
           STOP
        END IF
        !=======================================================================

        !=======================================================================
        name_transfo = mole%tab_Qtransfo(nb_Qtransfo)%name_transfo
        CALL string_uppercase_TO_lowercase(name_transfo)
        IF (name_transfo /= 'active') THEN
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) ' The last transfortmation MUST be "active".'
           write(out_unitp,*) 'name_transfo: ',mole%tab_Qtransfo(nb_Qtransfo)%name_transfo
           STOP
        END IF
        !=======================================================================


        mole%liste_QactTOQsym => mole%ActiveTransfo%list_QactTOQdyn
        mole%liste_QactTOQdyn => mole%ActiveTransfo%list_QactTOQdyn
        mole%liste_QsymTOQact => mole%ActiveTransfo%list_QdynTOQact
        mole%liste_QdynTOQact => mole%ActiveTransfo%list_QdynTOQact

        mole%name_Qdyn        => mole%tab_Qtransfo(nb_Qtransfo)%name_Qout

        CALL alloc_array(mole%nrho_OF_Qact,(/mole%nb_var/),             &
                        "mole%nrho_OF_Qact",name_sub)
        mole%nrho_OF_Qact(:) = 0
        CALL alloc_array(mole%nrho_OF_Qdyn,(/mole%nb_var/),             &
                        "mole%nrho_OF_Qdyn",name_sub)
        mole%nrho_OF_Qdyn(:) = 0

        IF (mole%tab_Qtransfo(nb_Qtransfo)%nb_Qin /= mole%nb_var) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nb_var /= nb_Qin',mole%nb_var,mole%tab_Qtransfo(nb_Qtransfo)%nb_Qin
          STOP
        END IF

        CALL type_var_analysis_OF_CoordType(mole)

        !=======================================================================
        !RPH transformation and type 21 or 22 coordinates are not compatible anymore.
        IF (associated(mole%RPHTransfo) .AND. mole%nb_inact2n > 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' asso mole%RPHTransfo ',associated(mole%RPHTransfo)
          write(out_unitp,*) ' mole%nb_inact2n      ',mole%nb_inact2n
          write(out_unitp,*) ' RPHTransfo and type 21 or 22 coordinates are not compatible anymore'
          STOP
        END IF

        mole%name_Qact => mole%tab_Qtransfo(nb_Qtransfo)%name_Qin

        IF (mole%nb_act < 1) THEN
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) ' nb_act < 1 !!',nb_act
           STOP
        END IF
        mole%tab_Qtransfo(:)%nb_act = mole%nb_act

        IF (debug) THEN
          DO it=1,nb_Qtransfo
            CALL Write_Qtransfo(mole%tab_Qtransfo(it),force_print=.TRUE.)
          END DO
        END IF

        IF(MPI_id==0) THEN
          write(out_unitp,*) '================================================='
          write(out_unitp,*) '================================================='
          write(out_unitp,*) 'END New Coordinates transformation'
          write(out_unitp,*) '================================================='
          write(out_unitp,*) '================================================='
          CALL flush_perso(out_unitp)
        ENDIF
!=======================================================================
!===== Old Coordinates transformation ==================================
!=======================================================================
      ELSE IF (zmat) THEN
        mole%nb_Qtransfo = 2  ! zmat + active
        IF (sym) mole%nb_Qtransfo = mole%nb_Qtransfo + 1
        IF (NM .OR. NM_TO_sym) mole%nb_Qtransfo = mole%nb_Qtransfo + 1

        IF(MPI_id==0) THEN
          write(out_unitp,*) '================================================='
          write(out_unitp,*) '================================================='
          write(out_unitp,*) 'Old Coordinates transformation',            &
                               mole%nb_Qtransfo
          write(out_unitp,*) '================================================='
        ENDIF

        CALL alloc_array(mole%tab_Qtransfo,(/mole%nb_Qtransfo/),        &
                        "mole%tab_Qtransfo",name_sub)

        !===============================================================
        !===== FIRST TRANSFO : zmat ====================================
        !===============================================================
        it = 1
        mole%itPrim = it

        mole%tab_Qtransfo(it)%name_transfo     = 'zmat'
        mole%tab_Qtransfo(it)%inTOout          = .TRUE.
        mole%tab_Qtransfo(it)%num_transfo      = it
        mole%tab_Qtransfo(it)%Primitive_Coord  = .TRUE.

        write(out_unitp,*) ' transfo: ',mole%tab_Qtransfo(it)%name_transfo
        write(out_unitp,*) ' num_transfo',mole%tab_Qtransfo(it)%num_transfo

        mole%tab_Qtransfo(it)%ZmatTransfo%cos_th = cos_th
        mole%tab_Qtransfo(it)%ZmatTransfo%nat0   = nat
        mole%tab_Qtransfo(it)%ZmatTransfo%nat    = nat + 1
        mole%tab_Qtransfo(it)%ZmatTransfo%nb_var = max(1,3*nat-6)
        mole%tab_Qtransfo(it)%ZmatTransfo%ncart  = 3*(nat+1)
        mole%tab_Qtransfo(it)%nb_Qin             = max(1,3*nat-6)
        mole%tab_Qtransfo(it)%nb_Qout            = 3*(nat+1)

        IF (debug) write(out_unitp,*) 'nat0,nat,nb_var,ncart',          &
                     mole%tab_Qtransfo(it)%ZmatTransfo%nat0,            &
                     mole%tab_Qtransfo(it)%ZmatTransfo%nat,             &
                     mole%tab_Qtransfo(it)%ZmatTransfo%nb_var,          &
                     mole%tab_Qtransfo(it)%ZmatTransfo%ncart

        mole%nat0            = mole%tab_Qtransfo(it)%ZmatTransfo%nat0
        mole%nat             = mole%tab_Qtransfo(it)%ZmatTransfo%nat
        mole%nb_var          = mole%tab_Qtransfo(it)%ZmatTransfo%nb_var
        mole%ncart           = mole%tab_Qtransfo(it)%ZmatTransfo%ncart
        mole%tab_Qtransfo(it)%Primitive_Coord = .TRUE.

        mole%tab_Qtransfo(it)%ZmatTransfo%New_Orient = New_Orient
        mole%tab_Qtransfo(it)%ZmatTransfo%vAt1       = vAt1
        mole%tab_Qtransfo(it)%ZmatTransfo%vAt2       = vAt2
        mole%tab_Qtransfo(it)%ZmatTransfo%vAt3       = vAt3

        CALL sub_Type_Name_OF_Qin(mole%tab_Qtransfo(it),"Qzmat")
        mole%tab_Qtransfo(it)%ZmatTransfo%type_Qin => mole%tab_Qtransfo(it)%type_Qin
        mole%tab_Qtransfo(it)%ZmatTransfo%name_Qin => mole%tab_Qtransfo(it)%name_Qin

        CALL Read_ZmatTransfo(mole%tab_Qtransfo(it)%ZmatTransfo,const_phys%mendeleev)

        CALL Set_masses_Z_TO_CoordType(mole,mole%tab_Qtransfo(it))


        ! for Qout type, name ....
        CALL alloc_array(mole%tab_Qtransfo(it)%type_Qout,(/3*(nat+1)/), &
                        "mole%tab_Qtransfo(it)%type_Qout",name_sub)

        CALL alloc_array(mole%tab_Qtransfo(it)%name_Qout,(/3*(nat+1)/), &
                        "mole%tab_Qtransfo(it)%name_Qout",name_sub)
        mole%tab_Qtransfo(it)%type_Qout(:) = 1 ! cartesian type

        DO i=1,mole%tab_Qtransfo(it)%nb_Qout
          iat = (i-1)/3 +1
          IF (mod(i,3) == 1) CALL make_nameQ(mole%tab_Qtransfo(it)%name_Qout(i),"X",iat,0)
          IF (mod(i,3) == 2) CALL make_nameQ(mole%tab_Qtransfo(it)%name_Qout(i),"Y",iat,0)
          IF (mod(i,3) == 0) CALL make_nameQ(mole%tab_Qtransfo(it)%name_Qout(i),"Z",iat,0)
        END DO

        IF (print_level > 1) CALL Write_QTransfo(mole%tab_Qtransfo(it))
        write(out_unitp,*) '------------------------------------------'
        CALL flush_perso(out_unitp)
        !===============================================================
        !===== SYM TRANSFO : linear ====================================
        !===============================================================
        IF (sym) THEN
          it = it + 1

          ! for Qout type, name ....
          mole%tab_Qtransfo(it)%type_Qout => mole%tab_Qtransfo(it-1)%type_Qin
          mole%tab_Qtransfo(it)%name_Qout => mole%tab_Qtransfo(it-1)%name_Qin

          mole%tab_Qtransfo(it)%name_transfo = 'linear'
          mole%tab_Qtransfo(it)%inTOout      = .TRUE.
          mole%tab_Qtransfo(it)%num_transfo  = it
          write(out_unitp,*) ' transfo: ',mole%tab_Qtransfo(it)%name_transfo
          write(out_unitp,*) ' num_transfo',mole%tab_Qtransfo(it)%num_transfo

          mole%tab_Qtransfo(it)%nb_Qin  = mole%nb_var
          mole%tab_Qtransfo(it)%nb_Qout = mole%nb_var
          mole%tab_Qtransfo(it)%LinearTransfo%inv = .FALSE.
          CALL flush_perso(out_unitp)

          CALL Read_LinearTransfo(mole%tab_Qtransfo(it)%LinearTransfo,  &
                                  mole%nb_var)
          CALL sub_Type_Name_OF_Qin(mole%tab_Qtransfo(it),"Qlinear")

          CALL Sub_Check_LinearTransfo(mole%tab_Qtransfo(it))

          IF (print_level > 1) CALL Write_QTransfo(mole%tab_Qtransfo(it))
          write(out_unitp,*) '------------------------------------------'
          CALL flush_perso(out_unitp)
        END IF

        !===============================================================
        !===== SYM TRANSFO : NM ========================================
        !===============================================================
        IF (NM .OR. NM_TO_sym) THEN
          it = it + 1
          mole%itNM = it

          ! for Qout type, name ....
          mole%tab_Qtransfo(it)%type_Qout => mole%tab_Qtransfo(it-1)%type_Qin
          mole%tab_Qtransfo(it)%name_Qout => mole%tab_Qtransfo(it-1)%name_Qin

          mole%tab_Qtransfo(it)%name_transfo = 'NM'
          mole%tab_Qtransfo(it)%inTOout      = .TRUE.
          mole%tab_Qtransfo(it)%num_transfo  = it
          write(out_unitp,*) ' transfo: ',mole%tab_Qtransfo(it)%name_transfo
          write(out_unitp,*) ' num_transfo',mole%tab_Qtransfo(it)%num_transfo

          mole%tab_Qtransfo(it)%nb_Qin  = mole%nb_var
          mole%tab_Qtransfo(it)%nb_Qout = mole%nb_var
          mole%tab_Qtransfo(it)%LinearTransfo%inv = .FALSE.
          CALL flush_perso(out_unitp)

          CALL alloc_array(mole%tab_Qtransfo(it)%NMTransfo,             &
                          "mole%tab_Qtransfo(it)%NMTransfo",name_sub)
          mole%NMTransfo => mole%tab_Qtransfo(it)%NMTransfo

          mole%tab_Qtransfo(it)%NMTransfo%purify_hess      = purify_hess
          mole%tab_Qtransfo(it)%NMTransfo%k_Half           = k_Half
          mole%tab_Qtransfo(it)%NMTransfo%hessian_old      = hessian_old
          mole%tab_Qtransfo(it)%NMTransfo%hessian_onthefly = hessian_onthefly
          mole%tab_Qtransfo(it)%NMTransfo%hessian_cart     = hessian_cart

          mole%tab_Qtransfo(it)%NMTransfo%file_hessian%name      = trim(file_hessian)
          mole%tab_Qtransfo(it)%NMTransfo%file_hessian%unit      = 0
          mole%tab_Qtransfo(it)%NMTransfo%file_hessian%formatted = .TRUE.
          mole%tab_Qtransfo(it)%NMTransfo%file_hessian%append    = .FALSE.
          mole%tab_Qtransfo(it)%NMTransfo%file_hessian%old       = hessian_old

          CALL Read_NMTransfo(mole%tab_Qtransfo(it)%NMTransfo,mole%nb_var)
          CALL sub_Type_Name_OF_Qin(mole%tab_Qtransfo(it),"QNM")

          CALL alloc_LinearTransfo(mole%tab_Qtransfo(it)%LinearTransfo,mole%nb_var)
          CALL mat_id(mole%tab_Qtransfo(it)%LinearTransfo%mat,mole%nb_var,mole%nb_var)
          CALL mat_id(mole%tab_Qtransfo(it)%LinearTransfo%mat_inv,mole%nb_var,mole%nb_var)
          CALL Sub_Check_LinearTransfo(mole%tab_Qtransfo(it))

          IF (print_level > 1) CALL Write_QTransfo(mole%tab_Qtransfo(it))
          write(out_unitp,*) '------------------------------------------'
          CALL flush_perso(out_unitp)
        END IF

        !===============================================================
        !===== LAST TRANSFO : active ===================================
        !===============================================================
        it = it + 1

        mole%tab_Qtransfo(it)%name_transfo = 'active'
        mole%tab_Qtransfo(it)%inTOout      = .TRUE.
        mole%tab_Qtransfo(it)%num_transfo  = it
        write(out_unitp,*) ' transfo: ',mole%tab_Qtransfo(it)%name_transfo
        write(out_unitp,*) ' num_transfo',mole%tab_Qtransfo(it)%num_transfo

        IF (it /= mole%nb_Qtransfo) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' The last read transformation is not the active one!'
          write(out_unitp,*) ' it,mole%nb_Qtransfo',it,mole%nb_Qtransfo
          STOP
        END IF
        mole%tab_Qtransfo(it)%nb_Qin  = mole%nb_var
        mole%tab_Qtransfo(it)%nb_Qout = mole%nb_var

        ! for Qout type, name ....
        mole%tab_Qtransfo(it)%type_Qout => mole%tab_Qtransfo(it-1)%type_Qin
        mole%tab_Qtransfo(it)%name_Qout => mole%tab_Qtransfo(it-1)%name_Qin


        CALL alloc_array(mole%tab_Qtransfo(it)%ActiveTransfo,           &
                        "mole%tab_Qtransfo(it)%ActiveTransfo",name_sub)
        CALL Read_ActiveTransfo(mole%tab_Qtransfo(it)%ActiveTransfo,    &
                                mole%nb_var)

        mole%ActiveTransfo => mole%tab_Qtransfo(it)%ActiveTransfo

        mole%ActiveTransfo%With_Tab_dnQflex = With_Tab_dnQflex

        mole%liste_QactTOQsym => mole%ActiveTransfo%list_QactTOQdyn
        mole%liste_QactTOQdyn => mole%ActiveTransfo%list_QactTOQdyn
        mole%liste_QsymTOQact => mole%ActiveTransfo%list_QdynTOQact
        mole%liste_QdynTOQact => mole%ActiveTransfo%list_QdynTOQact

        mole%name_Qdyn        => mole%tab_Qtransfo(it)%name_Qout

        CALL alloc_array(mole%nrho_OF_Qact,(/mole%nb_var/),             &
                        "mole%nrho_OF_Qact",name_sub)
        mole%nrho_OF_Qact(:) = 0

        CALL alloc_array(mole%nrho_OF_Qdyn,(/mole%nb_var/),             &
                        "mole%nrho_OF_Qdyn",name_sub)
        mole%nrho_OF_Qdyn(:) = 0

        CALL type_var_analysis_OF_CoordType(mole)
        mole%name_Qact        => mole%tab_Qtransfo(it)%name_Qin

        CALL flush_perso(out_unitp)


        IF (print_level > 1) CALL Write_QTransfo(mole%tab_Qtransfo(it))
        write(out_unitp,*) '------------------------------------------'

        write(out_unitp,*) '================================================='
        write(out_unitp,*) 'END Old Coordinates transformation'
        write(out_unitp,*) '================================================='
        write(out_unitp,*) '================================================='
      ELSE
        write(out_unitp,*) '================================================='
        write(out_unitp,*) '================================================='
        write(out_unitp,*) ' Old Coordinates transformation'
        write(out_unitp,*) '   with calc_f2_f1Q_ana (zmat=f)'
        write(out_unitp,*) '================================================='
        write(out_unitp,*) '================================================='

        mole%nb_var  = nb_var
        mole%nb_act  = nb_var
        mole%nb_act1 = nb_var

        allocate(mole%ActiveTransfo)
        CALL alloc_ActiveTransfo(mole%ActiveTransfo,nb_var)

        mole%liste_QactTOQsym => mole%ActiveTransfo%list_QactTOQdyn
        mole%liste_QactTOQdyn => mole%ActiveTransfo%list_QactTOQdyn
        mole%liste_QsymTOQact => mole%ActiveTransfo%list_QdynTOQact
        mole%liste_QdynTOQact => mole%ActiveTransfo%list_QdynTOQact

        mole%ActiveTransfo%list_QactTOQdyn(:) = 1

        IF (.NOT. associated(mole%name_Qdyn)) THEN
          CALL alloc_array(mole%name_Qdyn,(/nb_var/),"mole%name_Qdyn",name_sub)
        END IF

        DO i=1,nb_var
          mole%ActiveTransfo%list_QactTOQdyn(i) = i
          mole%ActiveTransfo%list_QdynTOQact(i) = i
          CALL make_nameQ(mole%name_Qdyn(i),'Qf2ana',i,0)
        END DO

        write(out_unitp,*) '================================================='
        write(out_unitp,*) 'END Old Coordinates transformation'
        write(out_unitp,*) '   with calc_f2_f1Q_ana (zmat=f)'
        write(out_unitp,*) '================================================='
        write(out_unitp,*) '================================================='
      END IF
      CALL flush_perso(out_unitp)

!=======================================================================
!===== set up: Mtot, Mtot_inv, d0sm =================================
!=======================================================================
      IF (mole%nb_Qtransfo /= -1) THEN
        CALL alloc_array(mole%d0sm,(/mole%ncart/),"mole%d0sm",name_sub)

        CALL alloc_array(mole%active_masses,(/mole%ncart/),             &
                        "mole%active_masses",name_sub)
        mole%active_masses(:) = 1

        mole%d0sm(:)   = sqrt(mole%masses(:))

        mole%Mtot      = sum(mole%masses(1:mole%ncart_act:3))
        mole%Mtot_inv  = ONE/mole%Mtot
      END IF

!=======================================================================
!===== Cart Coordinates transformation =================================
!=======================================================================
      IF (mole%Cart_transfo) THEN
        write(out_unitp,*) '================================================='
        write(out_unitp,*) '================================================='
        write(out_unitp,*) ' Read Cartesian transformation(s)'
        write(out_unitp,*) '================================================='

        CALL alloc_array(mole%tab_Cart_transfo,(/1/),                   &
                        "mole%tab_Cart_transfo",name_sub)

        CALL alloc_array(mole%tab_Cart_transfo(1)%CartesianTransfo%d0sm,&
                                                    (/mole%ncart_act/), &
                        "mole%tab_Cart_transfo(1)%CartesianTransfo%d0sm",name_sub)
        mole%tab_Cart_transfo(1)%CartesianTransfo%d0sm = mole%d0sm(1:mole%ncart_act)

        CALL alloc_array(mole%tab_Cart_transfo(1)%CartesianTransfo%masses_at,   &
                                                    (/ mole%nat_act /), &
                        "mole%tab_Cart_transfo(1)%CartesianTransfo%masses_at",name_sub)
        mole%tab_Cart_transfo(1)%CartesianTransfo%masses_at(:) = mole%masses(1:mole%ncart_act:3)
        mole%tab_Cart_transfo(1)%CartesianTransfo%nat_act      = mole%nat_act


        mole%tab_Cart_transfo(1)%num_transfo = 1
        CALL read_Qtransfo(mole%tab_Cart_transfo(1),mole%ncart_act,             &
                           mole%nb_extra_Coord,const_phys%mendeleev)

        IF (print_level > 1) CALL Write_CartesianTransfo(mole%tab_Cart_transfo(1)%CartesianTransfo)

        write(out_unitp,*) '================================================='
        write(out_unitp,*) '================================================='
      END IF

      mole%WriteCC         = WriteCC

      IF (NonGcteRange(1) > 0 .AND. NonGcteRange(1) > mole%nb_act1) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' Wrong value for NonGcteRange(:)',NonGcteRange(:)
        STOP
      END IF
      IF (NonGcteRange(1) > 0 .AND. NonGcteRange(2) > 0 .AND. NonGcteRange(2) > mole%nb_act1) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' Wrong value for NonGcteRange(:)',NonGcteRange(:)
        STOP
      END IF
      IF (NonGcteRange(1) > 0 .AND. NonGcteRange(2) < NonGcteRange(1)) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' Wrong value for NonGcteRange(:)',NonGcteRange(:)
        STOP
      END IF
      para_Tnum%NonGcteRange(:)   = NonGcteRange(:)

!=======================================================================
!=======================================================================
!=======================================================================
     !check is Tana is possible : nb_Qtransfo = 3
     para_Tnum%Tana = para_Tnum%Tana .AND. mole%nb_Qtransfo == 3
     !check is Tana is possible : 2st transfo poly
     name_transfo = mole%tab_Qtransfo(2)%name_transfo
     CALL string_uppercase_TO_lowercase(name_transfo)
     para_Tnum%Tana = para_Tnum%Tana .AND. name_transfo == "poly"
     ! we don't need to check the 1st and the last Qtransfo

!=======================================================================
!=======================================================================

!=======================================================================
!=======================================================================
!    special case if mole%ActiveTransfo%list_act_OF_Qdyn(i) = 200
     IF (count(mole%ActiveTransfo%list_act_OF_Qdyn(:) == 200) > 0) THEN
       para_Tnum%num_GG = .TRUE.
       para_Tnum%num_g  = .FALSE.
       para_Tnum%num_x  = .FALSE.
       write(out_unitp,*) ' WARNING: ActiveTransfo%list_act_OF_Qdyn(i) = 200'
       write(out_unitp,*) ' We are using numerical derivatives !!'
       write(out_unitp,*) ' num_GG,num_g,num_x',para_Tnum%num_GG,  &
                                      para_Tnum%num_g,para_Tnum%num_x
       write(out_unitp,*) ' stepT',para_Tnum%stepT
       write(out_unitp,*)
     END IF
!=======================================================================
!=======================================================================

!=======================================================================
!===== Ab initio parameters ============================================
!=======================================================================
      mole%charge          = charge
      mole%multiplicity    = multiplicity


     para_Tnum%para_PES_FromTnum%stepOp             = stepOp
     para_Tnum%para_PES_FromTnum%charge             = charge
     para_Tnum%para_PES_FromTnum%multiplicity       = multiplicity

     para_Tnum%para_PES_FromTnum%ab_initio_prog     = ab_initio_prog
     para_Tnum%para_PES_FromTnum%commande_unix      = commande_unix
     para_Tnum%para_PES_FromTnum%header             = header
     para_Tnum%para_PES_FromTnum%footer             = footer
     para_Tnum%para_PES_FromTnum%file_name_OTF      = file_name_OTF
     para_Tnum%para_PES_FromTnum%file_name_fchk     = file_name_fchk

     IF (len_trim(ab_initio_methEne) == 0) THEN
       ab_initio_methEne = ab_initio_meth
     END IF
     IF (len_trim(ab_initio_basisEne) == 0) THEN
       ab_initio_basisEne = ab_initio_basis
     END IF

     IF (len_trim(ab_initio_methDip) == 0) THEN
       ab_initio_methDip = ab_initio_meth
     END IF
     IF (len_trim(ab_initio_basisDip) == 0) THEN
       ab_initio_basisDip = ab_initio_basis
     END IF

     para_Tnum%para_PES_FromTnum%ab_initio_methEne  = ab_initio_methEne
     para_Tnum%para_PES_FromTnum%ab_initio_methDip  = ab_initio_methDip
     para_Tnum%para_PES_FromTnum%ab_initio_basisEne = ab_initio_basisEne
     para_Tnum%para_PES_FromTnum%ab_initio_basisDip = ab_initio_basisDip

!=======================================================================
!=======================================================================
      IF (debug) THEN
         write(out_unitp,*)
         DO it=1,mole%nb_Qtransfo
           write(out_unitp,*) '========================='
           write(out_unitp,*) ' Qtransfo,name_transfo',it," : ",        &
                                trim(mole%tab_Qtransfo(it)%name_transfo)

           write(out_unitp,*) 'asso name_Qout and type_Qout',           &
               associated(mole%tab_Qtransfo(it)%name_Qout),             &
                            associated(mole%tab_Qtransfo(it)%type_Qout)
           IF (associated(mole%tab_Qtransfo(it)%name_Qout) .AND.        &
                       associated(mole%tab_Qtransfo(it)%type_Qout)) THEN
             DO i_Q=1,mole%tab_Qtransfo(it)%nb_Qout
               write(out_unitp,*) 'i_Q,name_Qout,type_Qout',i_Q," ",    &
                           trim(mole%tab_Qtransfo(it)%name_Qout(i_Q)),  &
                                    mole%tab_Qtransfo(it)%type_Qout(i_Q)
               CALL flush_perso(out_unitp)
             END DO
           END IF

           write(out_unitp,*) 'asso name_Qin and type_Qin',             &
               associated(mole%tab_Qtransfo(it)%name_Qin),              &
                            associated(mole%tab_Qtransfo(it)%type_Qin)
           IF (associated(mole%tab_Qtransfo(it)%name_Qin) .AND.         &
                       associated(mole%tab_Qtransfo(it)%type_Qin)) THEN
             DO i_Q=1,mole%tab_Qtransfo(it)%nb_Qin
               write(out_unitp,*) 'i_Q,name_Qin,type_Qin',i_Q," ",      &
                           trim(mole%tab_Qtransfo(it)%name_Qin(i_Q)),   &
                                    mole%tab_Qtransfo(it)%type_Qin(i_Q)
               CALL flush_perso(out_unitp)
             END DO
           END IF

           write(out_unitp,*) '========================='
         END DO
      END IF

      write(out_unitp,*) 'END ',name_sub
      write(out_unitp,*) '================================================='
      write(out_unitp,*) '================================================='

      CALL flush_perso(out_unitp)
      END SUBROUTINE Read_CoordType
!===============================================================================

!================================================================
!       Copy two zmatrix type variables or two CoordType variables
!================================================================
  SUBROUTINE CoordType1TOCoordType2(mole1,mole2)

  ! for the CoordType and Tnum --------------------------------------
  TYPE (CoordType), intent(in)    :: mole1
  TYPE (CoordType), intent(inout) :: mole2

      integer :: it
      character (len=Name_len) :: name_transfo

      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      character (len=*), parameter :: name_sub='CoordType1TOCoordType2'

      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL flush_perso(out_unitp)
      END IF


      CALL dealloc_CoordType(mole2)

      IF (mole1%nat < 3 .OR. mole1%nb_var < 1 .OR.                      &
          mole1%nb_act < 1 .OR. mole1%ncart < 9) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' mole1 is probably not allocated !!'
        write(out_unitp,*) 'nat',mole1%nat
        write(out_unitp,*) 'nb_var',mole1%nb_var
        write(out_unitp,*) 'nb_act',mole1%nb_act
        write(out_unitp,*) 'ncart',mole1%ncart
        write(out_unitp,*) ' Check the Fortran source !!!'
        STOP
      END IF

      mole2%stepQ        = mole1%stepQ
      mole2%num_x        = mole1%num_x


      mole2%nb_act          = mole1%nb_act
      mole2%nb_var          = mole1%nb_var
      mole2%nb_extra_Coord  = mole1%nb_extra_Coord

      mole2%ndimG        = mole1%ndimG
      mole2%ncart        = mole1%ncart
      mole2%ncart_act    = mole1%ncart_act
      mole2%nat          = mole1%nat
      mole2%nat0         = mole1%nat0
      mole2%nat_act      = mole1%nat_act

      ! for the ab initio calculation
      mole2%charge       = mole1%charge
      mole2%multiplicity = mole1%multiplicity
      mole2%nb_elec      = mole1%nb_elec

      mole2%WriteCC      = mole1%WriteCC

      mole2%With_VecCOM     = mole1%With_VecCOM
      mole2%Without_Rot     = mole1%Without_Rot
      mole2%Centered_ON_CoM = mole1%Centered_ON_CoM
      mole2%cos_th          = mole1%cos_th


      mole2%Old_Qtransfo     = mole1%Old_Qtransfo
      mole2%Cart_transfo     = mole1%Cart_transfo
      mole2%nb_Qtransfo      = mole1%nb_Qtransfo
      mole2%itNM             = mole1%itNM
      mole2%itRPH            = mole1%itRPH
      mole2%itPrim           = mole1%itPrim
      mole2%opt_param        = mole1%opt_param

      IF (allocated(mole1%opt_Qdyn)) THEN
        CALL alloc_NParray(mole2%opt_Qdyn,shape(mole1%opt_Qdyn),        &
                          "mole2%opt_Qdyn",name_sub)
        mole2%opt_Qdyn(:) = mole1%opt_Qdyn(:)
      END IF


      CALL alloc_array(mole2%tab_Qtransfo,(/mole2%nb_Qtransfo/),        &
                      "mole2%tab_Qtransfo",name_sub)
      DO it=1,mole2%nb_Qtransfo
        !write(out_unitp,*) 'it',it ; flush(out_unitp)
        CALL Qtransfo1TOQtransfo2(mole1%tab_Qtransfo(it),               &
                                  mole2%tab_Qtransfo(it))

        name_transfo = mole2%tab_Qtransfo(it)%name_transfo
        CALL string_uppercase_TO_lowercase(name_transfo)
        SELECT CASE (name_transfo)
        CASE ("nm")
          mole2%NMTransfo => mole2%tab_Qtransfo(it)%NMTransfo
        CASE ("rph")
          mole2%RPHTransfo => mole2%tab_Qtransfo(it)%RPHTransfo
        CASE ('active')
          mole2%ActiveTransfo => mole2%tab_Qtransfo(it)%ActiveTransfo
        CASE ('zmat') ! it can be one of the last one
          mole2%Z       => mole2%tab_Qtransfo(it)%ZmatTransfo%Z
          mole2%symbole => mole2%tab_Qtransfo(it)%ZmatTransfo%symbole
          mole2%masses  => mole2%tab_Qtransfo(it)%ZmatTransfo%masses
        CASE ('bunch','bunch_poly') ! it has to be one of the last one
          mole2%Z       => mole2%tab_Qtransfo(it)%BunchTransfo%Z
          mole2%symbole => mole2%tab_Qtransfo(it)%BunchTransfo%symbole
          mole2%masses  => mole2%tab_Qtransfo(it)%BunchTransfo%masses
        CASE ('qtox_ana') ! it has to be one of the last one
          mole2%Z       => mole2%tab_Qtransfo(it)%QTOXanaTransfo%Z
          mole2%symbole => mole2%tab_Qtransfo(it)%QTOXanaTransfo%symbole
          mole2%masses  => mole2%tab_Qtransfo(it)%QTOXanaTransfo%masses
        CASE default
          CONTINUE
        END SELECT

        IF (it > 1) THEN
          mole2%tab_Qtransfo(it)%type_Qout => mole2%tab_Qtransfo(it-1)%type_Qin
          mole2%tab_Qtransfo(it)%name_Qout => mole2%tab_Qtransfo(it-1)%name_Qin
        END IF

      END DO



      IF (associated(mole1%RPHTransfo_inact2n)) THEN
        CALL alloc_array(mole2%RPHTransfo_inact2n,                      &
                        'mole2%RPHTransfo_inact2n',name_sub)
        CALL RPHTransfo1TORPHTransfo2(mole1%RPHTransfo_inact2n,         &
                                      mole2%RPHTransfo_inact2n)
      END IF

      CALL CurviRPH1_TO_CurviRPH2(mole1%CurviRPH,mole2%CurviRPH)

      !IF (mole1%Cart_transfo .AND. associated(mole1%tab_Cart_transfo)) THEN
      IF (associated(mole1%tab_Cart_transfo)) THEN

        CALL alloc_array(mole2%tab_Cart_transfo,                        &
                                         shape(mole1%tab_Cart_transfo), &
                        "mole2%tab_Cart_transfo",name_sub)

        DO it=1,size(mole1%tab_Cart_transfo)
          !write(out_unitp,*) 'it',it ; flush(out_unitp)
          CALL Qtransfo1TOQtransfo2(mole1%tab_Cart_transfo(it),         &
                                    mole2%tab_Cart_transfo(it))
        END DO
      END IF


      ! the 6 tables are true pointers
      mole2%liste_QactTOQsym => mole2%ActiveTransfo%list_QactTOQdyn
      mole2%liste_QactTOQdyn => mole2%ActiveTransfo%list_QactTOQdyn
      mole2%liste_QsymTOQact => mole2%ActiveTransfo%list_QdynTOQact
      mole2%liste_QdynTOQact => mole2%ActiveTransfo%list_QdynTOQact

      mole2%name_Qact        => mole2%tab_Qtransfo(mole2%nb_Qtransfo)%name_Qin
      mole2%name_Qdyn        => mole2%tab_Qtransfo(mole2%nb_Qtransfo)%name_Qout

      CALL alloc_array(mole2%nrho_OF_Qact,(/mole2%nb_var/),             &
                      "mole2%nrho_OF_Qact",name_sub)
      mole2%nrho_OF_Qact(:)  = mole1%nrho_OF_Qact(:)

      CALL alloc_array(mole2%nrho_OF_Qdyn,(/mole2%nb_var/),             &
                      "mole2%nrho_OF_Qdyn",name_sub)
      mole2%nrho_OF_Qdyn(:)  = mole1%nrho_OF_Qdyn(:)

      mole2%nb_act1          = mole1%nb_act1
      mole2%nb_inact2n       = mole1%nb_inact2n
      mole2%nb_inact21       = mole1%nb_inact21
      mole2%nb_inact22       = mole1%nb_inact22
      mole2%nb_inact20       = mole1%nb_inact20
      mole2%nb_inact         = mole1%nb_inact
      mole2%nb_inact31       = mole1%nb_inact31
      mole2%nb_rigid0        = mole1%nb_rigid0
      mole2%nb_rigid100      = mole1%nb_rigid100
      mole2%nb_rigid         = mole1%nb_rigid

      CALL alloc_array(mole2%active_masses,(/mole2%ncart/),             &
                      "mole2%active_masses",name_sub)
      mole2%active_masses    = mole1%active_masses

      CALL alloc_array(mole2%d0sm,(/mole2%ncart/),                      &
                      "mole2%d0sm",name_sub)
      mole2%d0sm             = mole1%d0sm
      mole2%Mtot             = mole1%Mtot
      mole2%Mtot_inv         = mole1%Mtot_inv


      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF

  END SUBROUTINE CoordType1TOCoordType2
  SUBROUTINE mole1TOmole2(mole1,mole2)
  TYPE (zmatrix), intent(in)    :: mole1
  TYPE (zmatrix), intent(inout) :: mole2

    CALL CoordType1TOCoordType2(mole1%CoordType,mole2%CoordType)

  END SUBROUTINE mole1TOmole2
  !!@description: Copy two CoordType variables
  !!@param: TODO
  SUBROUTINE CoordType2_TO_CoordType1(mole1,mole2)

  ! for the CoordType and Tnum --------------------------------------
  CLASS (CoordType), intent(inout) :: mole1
  TYPE (CoordType),  intent(in)    :: mole2

      integer :: it
      character (len=Name_len) :: name_transfo

      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      character (len=*), parameter :: name_sub='CoordType2_TO_CoordType1'

      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL flush_perso(out_unitp)
      END IF


      CALL dealloc_CoordType(mole1)

      IF (mole2%nat < 3 .OR. mole2%nb_var < 1 .OR.                      &
          mole2%nb_act < 1 .OR. mole2%ncart < 9) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' mole2 is probably not allocated !!'
        write(out_unitp,*) 'nat',mole2%nat
        write(out_unitp,*) 'nb_var',mole2%nb_var
        write(out_unitp,*) 'nb_act',mole2%nb_act
        write(out_unitp,*) 'ncart',mole2%ncart
        write(out_unitp,*) ' Check the Fortran source !!!'
        STOP
      END IF

      mole1%stepQ        = mole2%stepQ
      mole1%num_x        = mole2%num_x


      mole1%nb_act          = mole2%nb_act
      mole1%nb_var          = mole2%nb_var
      mole1%nb_extra_Coord  = mole2%nb_extra_Coord

      mole1%ndimG        = mole2%ndimG
      mole1%ncart        = mole2%ncart
      mole1%ncart_act    = mole2%ncart_act
      mole1%nat          = mole2%nat
      mole1%nat0         = mole2%nat0
      mole1%nat_act      = mole2%nat_act

      ! for the ab initio calculation
      mole1%charge       = mole2%charge
      mole1%multiplicity = mole2%multiplicity
      mole1%nb_elec      = mole2%nb_elec

      mole1%WriteCC      = mole2%WriteCC

      mole1%With_VecCOM     = mole2%With_VecCOM
      mole1%Without_Rot     = mole2%Without_Rot
      mole1%Centered_ON_CoM = mole2%Centered_ON_CoM
      mole1%cos_th          = mole2%cos_th


      mole1%Old_Qtransfo     = mole2%Old_Qtransfo
      mole1%Cart_transfo     = mole2%Cart_transfo
      mole1%nb_Qtransfo      = mole2%nb_Qtransfo
      mole1%itNM             = mole2%itNM
      mole1%itRPH            = mole2%itRPH
      mole1%itPrim           = mole2%itPrim
      mole1%opt_param        = mole2%opt_param

      IF (allocated(mole2%opt_Qdyn)) THEN
        CALL alloc_NParray(mole1%opt_Qdyn,shape(mole2%opt_Qdyn),        &
                          "mole1%opt_Qdyn",name_sub)
        mole1%opt_Qdyn(:) = mole2%opt_Qdyn(:)
      END IF


      CALL alloc_array(mole1%tab_Qtransfo,(/mole1%nb_Qtransfo/),        &
                      "mole1%tab_Qtransfo",name_sub)
      DO it=1,mole1%nb_Qtransfo
        !write(out_unitp,*) 'it',it ; flush(out_unitp)
        CALL Qtransfo1TOQtransfo2(mole2%tab_Qtransfo(it),               &
                                  mole1%tab_Qtransfo(it))

        name_transfo = mole1%tab_Qtransfo(it)%name_transfo
        CALL string_uppercase_TO_lowercase(name_transfo)
        SELECT CASE (name_transfo)
        CASE ("nm")
          mole1%NMTransfo => mole1%tab_Qtransfo(it)%NMTransfo
        CASE ("rph")
          mole1%RPHTransfo => mole1%tab_Qtransfo(it)%RPHTransfo
        CASE ('active')
          mole1%ActiveTransfo => mole1%tab_Qtransfo(it)%ActiveTransfo
        CASE ('zmat') ! it can be one of the last one
          mole1%Z       => mole1%tab_Qtransfo(it)%ZmatTransfo%Z
          mole1%symbole => mole1%tab_Qtransfo(it)%ZmatTransfo%symbole
          mole1%masses  => mole1%tab_Qtransfo(it)%ZmatTransfo%masses
        CASE ('bunch','bunch_poly') ! it has to be one of the last one
          mole1%Z       => mole1%tab_Qtransfo(it)%BunchTransfo%Z
          mole1%symbole => mole1%tab_Qtransfo(it)%BunchTransfo%symbole
          mole1%masses  => mole1%tab_Qtransfo(it)%BunchTransfo%masses
        CASE ('qtox_ana') ! it has to be one of the last one
          mole1%Z       => mole1%tab_Qtransfo(it)%QTOXanaTransfo%Z
          mole1%symbole => mole1%tab_Qtransfo(it)%QTOXanaTransfo%symbole
          mole1%masses  => mole1%tab_Qtransfo(it)%QTOXanaTransfo%masses
        CASE default
          CONTINUE
        END SELECT

        IF (it > 1) THEN
          mole1%tab_Qtransfo(it)%type_Qout => mole1%tab_Qtransfo(it-1)%type_Qin
          mole1%tab_Qtransfo(it)%name_Qout => mole1%tab_Qtransfo(it-1)%name_Qin
        END IF

      END DO



      IF (associated(mole2%RPHTransfo_inact2n)) THEN
        CALL alloc_array(mole1%RPHTransfo_inact2n,                      &
                        'mole1%RPHTransfo_inact2n',name_sub)
        CALL RPHTransfo1TORPHTransfo2(mole2%RPHTransfo_inact2n,         &
                                      mole1%RPHTransfo_inact2n)
      END IF

      CALL CurviRPH1_TO_CurviRPH2(mole2%CurviRPH,mole1%CurviRPH)

      !IF (mole2%Cart_transfo .AND. associated(mole2%tab_Cart_transfo)) THEN
      IF (associated(mole2%tab_Cart_transfo)) THEN

        CALL alloc_array(mole1%tab_Cart_transfo,                        &
                                         shape(mole2%tab_Cart_transfo), &
                        "mole1%tab_Cart_transfo",name_sub)

        DO it=1,size(mole2%tab_Cart_transfo)
          !write(out_unitp,*) 'it',it ; flush(out_unitp)
          CALL Qtransfo1TOQtransfo2(mole2%tab_Cart_transfo(it),         &
                                    mole1%tab_Cart_transfo(it))
        END DO
      END IF


      ! the 6 tables are true pointers
      mole1%liste_QactTOQsym => mole1%ActiveTransfo%list_QactTOQdyn
      mole1%liste_QactTOQdyn => mole1%ActiveTransfo%list_QactTOQdyn
      mole1%liste_QsymTOQact => mole1%ActiveTransfo%list_QdynTOQact
      mole1%liste_QdynTOQact => mole1%ActiveTransfo%list_QdynTOQact

      mole1%name_Qact        => mole1%tab_Qtransfo(mole1%nb_Qtransfo)%name_Qin
      mole1%name_Qdyn        => mole1%tab_Qtransfo(mole1%nb_Qtransfo)%name_Qout

      CALL alloc_array(mole1%nrho_OF_Qact,(/mole1%nb_var/),             &
                      "mole1%nrho_OF_Qact",name_sub)
      mole1%nrho_OF_Qact(:)  = mole2%nrho_OF_Qact(:)

      CALL alloc_array(mole1%nrho_OF_Qdyn,(/mole1%nb_var/),             &
                      "mole1%nrho_OF_Qdyn",name_sub)
      mole1%nrho_OF_Qdyn(:)  = mole2%nrho_OF_Qdyn(:)

      mole1%nb_act1          = mole2%nb_act1
      mole1%nb_inact2n       = mole2%nb_inact2n
      mole1%nb_inact21       = mole2%nb_inact21
      mole1%nb_inact22       = mole2%nb_inact22
      mole1%nb_inact20       = mole2%nb_inact20
      mole1%nb_inact         = mole2%nb_inact
      mole1%nb_inact31       = mole2%nb_inact31
      mole1%nb_rigid0        = mole2%nb_rigid0
      mole1%nb_rigid100      = mole2%nb_rigid100
      mole1%nb_rigid         = mole2%nb_rigid

      CALL alloc_array(mole1%active_masses,(/mole1%ncart/),             &
                      "mole1%active_masses",name_sub)
      mole1%active_masses    = mole2%active_masses

      CALL alloc_array(mole1%d0sm,(/mole1%ncart/),                      &
                      "mole1%d0sm",name_sub)
      mole1%d0sm             = mole2%d0sm
      mole1%Mtot             = mole2%Mtot
      mole1%Mtot_inv         = mole2%Mtot_inv


      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF

     END SUBROUTINE CoordType2_TO_CoordType1

  SUBROUTINE Tnum2_TO_Tnum1(Tnum1,Tnum2)

  CLASS (Tnum), intent(inout) :: Tnum1
  TYPE (Tnum),  intent(in)    :: Tnum2

  Tnum1%stepT                = Tnum2%stepT
  Tnum1%num_GG               = Tnum2%num_GG
  Tnum1%num_g                = Tnum2%num_g
  Tnum1%num_x                = Tnum2%num_x
  Tnum1%Gdiago               = Tnum2%Gdiago
  Tnum1%Gcte                 = Tnum2%Gcte
  Tnum1%NonGcteRange         = Tnum2%NonGcteRange
  Tnum1%Inertia              = Tnum2%Inertia
  Tnum1%nrho                 = Tnum2%nrho
  Tnum1%vep_type             = Tnum2%vep_type
  Tnum1%JJ                   = Tnum2%JJ
  Tnum1%WriteT               = Tnum2%WriteT
  Tnum1%With_Cart_Transfo    = Tnum2%With_Cart_Transfo

  Tnum1%Tana                 = Tnum2%Tana
  Tnum1%MidasCppForm         = Tnum2%MidasCppForm
  Tnum1%MCTDHForm            = Tnum2%MCTDHForm
  Tnum1%LaTeXForm            = Tnum2%LaTeXForm
  Tnum1%VSCFForm             = Tnum2%VSCFForm
  Tnum1%f2f1_ana             = Tnum2%f2f1_ana
  Tnum1%KEO_TalyorOFQinact2n = Tnum2%KEO_TalyorOFQinact2n

  Tnum1%para_PES_FromTnum    = Tnum2%para_PES_FromTnum ! OK no pointers

  Tnum1%TWOxKEO              = Tnum2%TWOxKEO
  Tnum1%ExpandTWOxKEO        = Tnum2%ExpandTWOxKEO

  IF (associated(Tnum2%Gref)) THEN
    CALL alloc_array(Tnum1%Gref,shape(Tnum2%Gref),'Tnum1%Gref','Tnum2_TO_Tnum1')
    Tnum1%Gref(:,:) = Tnum2%Gref
  END IF

  END SUBROUTINE Tnum2_TO_Tnum1
  SUBROUTINE dealloc_Tnum(para_Tnum)

  CLASS (Tnum), intent(inout) :: para_Tnum

  para_Tnum%stepT   = ONETENTH**4
  para_Tnum%num_GG  = .FALSE.
  para_Tnum%num_g   = .FALSE.
  para_Tnum%num_x   = .FALSE.

  para_Tnum%Gdiago  = .FALSE.
  para_Tnum%Gcte    = .FALSE.
  para_Tnum%NonGcteRange(2) = 0
  para_Tnum%Inertia = .TRUE.

  para_Tnum%nrho              = 2
  para_Tnum%vep_type          = -1
  para_Tnum%JJ                = 0
  para_Tnum%WriteT            = .FALSE.
  para_Tnum%With_Cart_Transfo = .TRUE.

  para_Tnum%Tana              = .FALSE.
  para_Tnum%MidasCppForm      = .FALSE.
  para_Tnum%MCTDHForm         = .FALSE.
  para_Tnum%LaTeXForm         = .FALSE.
  para_Tnum%VSCFForm          = .FALSE.

  para_Tnum%f2f1_ana          = .FALSE.

  para_Tnum%KEO_TalyorOFQinact2n = -1

  !CALL dealloc_PES_Frompara_Tnum(para_Tnum%para_PES_FromTnum)
  IF (allocated(para_Tnum%para_PES_FromTnum%nDFit_Scalar_Op_name_Fit)) THEN
    deallocate(para_Tnum%para_PES_FromTnum%nDFit_Scalar_Op_name_Fit)
  END IF

  CALL delete_op(para_Tnum%TWOxKEO)
  CALL delete_op(para_Tnum%ExpandTWOxKEO)

  IF (associated(para_Tnum%Gref)) THEN
    CALL dealloc_array(para_Tnum%Gref,'para_Tnum%Gref','dealloc_Tnum')
  END IF

  END SUBROUTINE dealloc_Tnum
!================================================================
!      check charge multiplicity
!================================================================
      SUBROUTINE check_charge(mole)
       TYPE (CoordType), intent(inout) :: mole

       IF (print_level > 1) THEN
         write(out_unitp,*)
         write(out_unitp,*) ' - check charge multiplicity -------------'
       END IF

       mole%nb_elec = sum(mole%Z(1:mole%nat0)) + mole%charge
       IF (mole%multiplicity <0) THEN
         mole%multiplicity = 1 + mod(mole%nb_elec,2)
       END IF

       IF (mod(mole%nb_elec+1,2) /= mod(mole%multiplicity,2)) THEN
         write(out_unitp,*) 'ERROR in check_charge'
         write(out_unitp,*) 'nb_elec,charge,multiplicity INCOMPATIBLE'
         write(out_unitp,*) mole%nb_elec,mole%charge,mole%multiplicity
         STOP
       END IF

       IF (print_level > 1) THEN
         write(out_unitp,*) 'nb_elec,charge,multiplicity',                &
                   mole%nb_elec,mole%charge,mole%multiplicity
         write(out_unitp,*) 'Z:',mole%Z(1:mole%nat0)
         write(out_unitp,*) ' - End check charge multiplicity ----------'
       END IF

      END SUBROUTINE check_charge

!================================================================
!      analysis of the variable type
!
!       analysis of list_act_OF_Qdyn(.) to define ActiveTransfo%list_QactTOQdyn(.) and ActiveTransfo%list_QdynTOQact(.)
!================================================================
  SUBROUTINE type_var_analysis_OF_CoordType(mole)
  TYPE (CoordType), intent(inout) :: mole


      integer :: iv_inact20,iv_rigid0,iv_inact21,iv_act1,iv_inact22
      integer :: iv_inact31,iv_rigid100
      integer :: n_test
      integer :: i,ii,i_Q,iQin,iQout
      integer :: nb_Qtransfo,it

      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'type_var_analysis_OF_CoordType'
      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.

      IF (debug) write(out_unitp,*) 'BEGINNING ',name_sub
      IF (print_level > 0 .OR. debug) THEN
         write(out_unitp,*) '-analysis of the variable type ---'
         write(out_unitp,*) '  mole%list_act_OF_Qdyn',                  &
                                     mole%ActiveTransfo%list_act_OF_Qdyn
      END IF

!      ! first transfert list_act_OF_Qdyn to RPHtransfo
!      IF (associated(mole%RPHTransfo)) THEN
!      IF (mole%RPHTransfo%option == 0) THEN
!
!        ! This transformation is done only if option==0.
!        !    For option=1, everything is already done
!
!        IF (.NOT. (mole%RPHTransfo%init_Qref .OR. mole%RPHTransfo%init)) THEN
!          CALL Set_RPHTransfo(mole%RPHTransfo,mole%ActiveTransfo%list_act_OF_Qdyn)
!        END IF
!
!        DO i=1,mole%nb_var
!          IF (mole%ActiveTransfo%list_act_OF_Qdyn(i) == 21)             &
!                              mole%ActiveTransfo%list_act_OF_Qdyn(i) = 1
!        END DO
!      END IF
!      END IF

      mole%nb_act1    =                                                 &
       count(abs(mole%ActiveTransfo%list_act_OF_Qdyn(1:mole%nb_var))==1)
      mole%nb_inact20 = count(mole%ActiveTransfo%list_act_OF_Qdyn==20) +&
                        count(mole%ActiveTransfo%list_act_OF_Qdyn==200)
      mole%nb_inact21 = count(mole%ActiveTransfo%list_act_OF_Qdyn==21) +&
                        count(mole%ActiveTransfo%list_act_OF_Qdyn==210)
      mole%nb_inact31 = count(mole%ActiveTransfo%list_act_OF_Qdyn==31)
      mole%nb_inact22 = count(mole%ActiveTransfo%list_act_OF_Qdyn==22)
      mole%nb_rigid0  =                                                 &
            count(mole%ActiveTransfo%list_act_OF_Qdyn(1:mole%nb_var)==0)
      mole%nb_rigid100= count(mole%ActiveTransfo%list_act_OF_Qdyn==100)



      n_test = mole%nb_act1    +                                        &
               mole%nb_inact31 +                                        &
               mole%nb_inact20 + mole%nb_inact21 + mole%nb_inact22 +    &
               mole%nb_rigid0  + mole%nb_rigid100

      IF (n_test /= mole%nb_var) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' type of variable impossible',              &
                                    mole%ActiveTransfo%list_act_OF_Qdyn
        write(out_unitp,*) 'nb_act1+nb_inact+nb_rigid... and nb_var',   &
                              n_test,mole%nb_var
        STOP
      END IF
!---------------------------------------------------------------



!--------------------------------------------------------------
!----- variable list in order: active, inactive harmo ...
!      =>  ActiveTransfo%list_QactTOQdyn

       iv_act1    = 0
       iv_inact31 = iv_act1    + mole%nb_act1
       iv_inact22 = iv_inact31 + mole%nb_inact31
       iv_inact21 = iv_inact22 + mole%nb_inact22
       iv_rigid100= iv_inact21 + mole%nb_inact21
       iv_inact20 = iv_rigid100+ mole%nb_rigid100
       iv_rigid0  = iv_inact20 + mole%nb_inact20


       mole%ActiveTransfo%list_QactTOQdyn(:) = 0

       DO i=1,mole%nb_var

        SELECT CASE (mole%ActiveTransfo%list_act_OF_Qdyn(i))
        CASE (1,-1)
           iv_act1 = iv_act1 + 1
           mole%ActiveTransfo%list_QactTOQdyn(iv_act1) = i
        CASE (31)
           iv_inact31 = iv_inact31 + 1
           mole%ActiveTransfo%list_QactTOQdyn(iv_inact31) = i
        CASE (21,210)
           iv_inact21 = iv_inact21 + 1
           mole%ActiveTransfo%list_QactTOQdyn(iv_inact21) = i
        CASE (22)
           iv_inact22 = iv_inact22 + 1
           mole%ActiveTransfo%list_QactTOQdyn(iv_inact22) = i
        CASE (20,200)
           iv_inact20 = iv_inact20 + 1
           mole%ActiveTransfo%list_QactTOQdyn(iv_inact20) = i
        CASE (0)
           iv_rigid0 = iv_rigid0 + 1
           mole%ActiveTransfo%list_QactTOQdyn(iv_rigid0) = i
        CASE (100)
           iv_rigid100 = iv_rigid100 + 1
           mole%ActiveTransfo%list_QactTOQdyn(iv_rigid100) = i
         CASE default
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) ' act_OF_Qdyn',mole%ActiveTransfo%list_act_OF_Qdyn(i),' impossible'
           STOP
         END SELECT

       END DO

!     variable list in CoordType order ------------------------
!     =>  ActiveTransfo%list_QdynTOQact
!
      mole%ActiveTransfo%list_QdynTOQact(:) = 0
      DO i=1,size(mole%ActiveTransfo%list_QdynTOQact)
        mole%ActiveTransfo%list_QdynTOQact(mole%ActiveTransfo%list_QactTOQdyn(i)) = i
      END DO

      mole%nb_inact2n = mole%nb_inact21 + mole%nb_inact22
      mole%nb_act     = mole%nb_act1    + mole%nb_inact31 + mole%nb_inact2n
      mole%nb_inact   = mole%nb_inact20 + mole%nb_inact2n
      mole%nb_rigid   = mole%nb_rigid0  + mole%nb_rigid100

      mole%ndimG      = mole%nb_act     + 3
      IF (mole%With_VecCOM) mole%ndimG = mole%ndimG + 3 ! with the COM
      IF (mole%Without_Rot) mole%ndimG  = mole%nb_act


      nb_Qtransfo = mole%nb_Qtransfo
      mole%ActiveTransfo%nb_var      = mole%nb_var

      mole%ActiveTransfo%nb_act      = mole%nb_act
      mole%ActiveTransfo%nb_act1     = mole%nb_act1

      mole%ActiveTransfo%nb_inact2n  = mole%nb_inact2n
      mole%ActiveTransfo%nb_inact21  = mole%nb_inact21
      mole%ActiveTransfo%nb_inact22  = mole%nb_inact22
      mole%ActiveTransfo%nb_inact20  = mole%nb_inact20
      mole%ActiveTransfo%nb_inact    = mole%nb_inact

      mole%ActiveTransfo%nb_inact31  = mole%nb_inact31

      mole%ActiveTransfo%nb_rigid0   = mole%nb_rigid0
      mole%ActiveTransfo%nb_rigid100 = mole%nb_rigid100
      mole%ActiveTransfo%nb_rigid    = mole%nb_rigid

      ! list_QactTOQdyn has been changed (not the same Qact order),
      !   then ActiveTransfo%Qact0 has to be changed as well.
      DO iQin=1,mole%nb_var
         iQout = mole%ActiveTransfo%list_QactTOQdyn(iQin)
         mole%ActiveTransfo%Qact0(iQin) = mole%ActiveTransfo%Qdyn0(iQout)
      END DO


      ! set nrho_OF_Qact in function of nrho_OF_Qdyn and list_QactTOQdyn
      CALL sub_Type_Name_OF_Qin(mole%tab_Qtransfo(nb_Qtransfo),"Qact")
      DO iQin=1,mole%nb_var

         iQout = mole%ActiveTransfo%list_QactTOQdyn(iQin)
         !write(out_unitp,*) 'iQin,iQout,type_Qout',iQin,iQout,mole%tab_Qtransfo(nb_Qtransfo)%type_Qout(iQout)
         !flush(out_unitp)
         mole%tab_Qtransfo(nb_Qtransfo)%type_Qin(iQin) =                &
                        mole%tab_Qtransfo(nb_Qtransfo)%type_Qout(iQout)
         mole%tab_Qtransfo(nb_Qtransfo)%name_Qin(iQin) =                &
                         mole%tab_Qtransfo(nb_Qtransfo)%name_Qout(iQout)
         mole%nrho_OF_Qact(iQin) = mole%nrho_OF_Qdyn(iQout)

      END DO

      mole%tab_Qtransfo(:)%nb_act = mole%nb_act

      IF (print_level > 0 .OR. debug) THEN
        write(out_unitp,*) 'mole%nrho_OF_Qact(:)    ',mole%nrho_OF_Qact(:)
        write(out_unitp,*) 'mole%...%list_QactTOQdyn',mole%ActiveTransfo%list_QactTOQdyn
        write(out_unitp,*) 'mole%...%list_QdynTOQact',mole%ActiveTransfo%list_QdynTOQact
        write(out_unitp,*)
        write(out_unitp,*) 'mole%...%type_Qout(:)',mole%tab_Qtransfo(nb_Qtransfo)%type_Qout(:)
        write(out_unitp,*) 'mole%...%type_Qin(:)',mole%tab_Qtransfo(nb_Qtransfo)%type_Qin(:)
        write(out_unitp,*)
        write(out_unitp,*) '- End analysis of the variable type ---'
      END IF
      IF (debug) write(out_unitp,*) 'END ',name_sub

  END SUBROUTINE type_var_analysis_OF_CoordType

!================================================================
! ++    write the hamiltonian f2 f1
!================================================================
      SUBROUTINE Write_f2f1vep(f2ij,f1i,vep,rho,nb_act)
      USE mod_system
      IMPLICIT NONE

       integer nb_act
       real (kind=Rkind) ::  f2ij(nb_act,nb_act),f1i(nb_act),vep,rho

       integer i

       write(out_unitp,*) 'nb_act',nb_act
       write(out_unitp,11) vep,rho
 11    format(' vep rho = ',2f30.15)
       write(out_unitp,*)
       CALL Write_Vec(f1i,out_unitp,10,name_info=' f1i = ')
       write(out_unitp,*)
       write(out_unitp,*) ' f2ij'
       CALL Write_Mat(f2ij,out_unitp,4)
       write(out_unitp,*)

       end subroutine Write_f2f1vep
!================================================================
! ++    write the hamiltonian Tcor Trot
!================================================================
      SUBROUTINE Write_TcorTrot(Tcor2,Tcor1,Trot,nb_act)
      USE mod_system
      IMPLICIT NONE

       integer :: nb_act
       real (kind=Rkind) ::  Tcor2(nb_act,3),Tcor1(3),Trot(3,3)

       integer :: i

       write(out_unitp,*) 'nb_act',nb_act
       write(out_unitp,*)
       write(out_unitp,*) 'Tcor2 * i* Jx'
       write(out_unitp,21) (Tcor2(i,1),i=1,nb_act)
       write(out_unitp,*) 'Tcor2 * i* Jy'
       write(out_unitp,21) (Tcor2(i,2),i=1,nb_act)
       write(out_unitp,*) 'Tcor2 * i* Jz'
       write(out_unitp,21) (Tcor2(i,3),i=1,nb_act)
 21    format(' Tcor2 = ',2x,10(f18.10,3x))

       write(out_unitp,*)
       write(out_unitp,*) 'Tcor1*i*Jx, Tcor1*i*Jy, Tcor1*i*Jz'
       write(out_unitp,31) Tcor1(1),Tcor1(2),Tcor1(3)
 31    format(' Tcor1 = ',2x,3(f18.10,3x))

       write(out_unitp,*)
       write(out_unitp,*) 'Trot'
       CALL Write_Mat(Trot,out_unitp,3)
       write(out_unitp,*)

       end subroutine Write_TcorTrot

  SUBROUTINE Sub_paraRPH_TO_CoordType(mole)
  TYPE (CoordType), intent(inout) :: mole

      integer :: it
      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.

      IF (.NOT. associated(mole%RPHTransfo)) RETURN
      mole%RPHTransfo%option = 0 ! If option/=0, we set up to option 0 to be able to perform
                           ! the transformations Sub_paraRPH_TO_mole and Sub_mole_TO_paraRPH

      IF (debug) write(out_unitp,*) 'BEGINNING Sub_paraRPH_TO_CoordType'

      mole%ActiveTransfo%list_act_OF_Qdyn(:) = mole%RPHTransfo%list_act_OF_Qdyn(:)
      mole%nb_act1                           = mole%RPHTransfo%nb_act1
      mole%nb_inact21                        = mole%RPHTransfo%nb_inact21


      mole%nb_inact2n = mole%nb_inact21 + mole%nb_inact22
      mole%nb_act     = mole%nb_act1 + mole%nb_inact31 + mole%nb_inact2n
      mole%nb_inact   = mole%nb_inact20 + mole%nb_inact2n
      mole%nb_rigid   = mole%nb_rigid0 + mole%nb_rigid100
      mole%ndimG      = mole%nb_act+3
      IF (mole%With_VecCOM) mole%ndimG = mole%ndimG + 3 ! with the COM
      IF (mole%Without_Rot) mole%ndimG  = mole%nb_act


      mole%ActiveTransfo%nb_var      = mole%nb_var

      mole%ActiveTransfo%nb_act      = mole%nb_act
      mole%ActiveTransfo%nb_act1     = mole%nb_act1

      mole%ActiveTransfo%nb_inact2n  = mole%nb_inact2n
      mole%ActiveTransfo%nb_inact21  = mole%nb_inact21
      mole%ActiveTransfo%nb_inact22  = mole%nb_inact22
      mole%ActiveTransfo%nb_inact20  = mole%nb_inact20
      mole%ActiveTransfo%nb_inact    = mole%nb_inact

      mole%ActiveTransfo%nb_inact31  = mole%nb_inact31

      mole%ActiveTransfo%nb_rigid0   = mole%nb_rigid0
      mole%ActiveTransfo%nb_rigid100 = mole%nb_rigid100
      mole%ActiveTransfo%nb_rigid    = mole%nb_rigid

      IF (debug) write(out_unitp,*) 'END Sub_paraRPH_TO_CoordType'


  END SUBROUTINE Sub_paraRPH_TO_CoordType
  SUBROUTINE Sub_CoordType_TO_paraRPH(mole)
  TYPE (CoordType), intent(inout) :: mole

      integer :: i
      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.

      IF (.NOT. associated(mole%RPHTransfo)) RETURN
      IF (mole%RPHTransfo%option /= 0) RETURN
      ! This transformation is done only if option==0.
      !    For option=1, everything is already done

      IF (debug) write(out_unitp,*) 'BEGINNING Sub_CoordType_TO_paraRPH'


      DO i=1,mole%nb_var
        IF (mole%ActiveTransfo%list_act_OF_Qdyn(i) == 21)               &
        mole%ActiveTransfo%list_act_OF_Qdyn(i) = 1
      END DO
      mole%nb_act1 = mole%nb_act1 + mole%nb_inact21
      mole%nb_inact21 = 0



      mole%nb_inact2n = mole%nb_inact21 + mole%nb_inact22
      mole%nb_act     = mole%nb_act1 + mole%nb_inact31 + mole%nb_inact2n
      mole%nb_inact   = mole%nb_inact20 + mole%nb_inact2n
      mole%nb_rigid   = mole%nb_rigid0 + mole%nb_rigid100
      mole%ndimG      = mole%nb_act+3
      IF (mole%With_VecCOM) mole%ndimG = mole%ndimG + 3 ! with the COM
      IF (mole%Without_Rot) mole%ndimG      = mole%nb_act


      mole%ActiveTransfo%nb_var      = mole%nb_var

      mole%ActiveTransfo%nb_act      = mole%nb_act
      mole%ActiveTransfo%nb_act1     = mole%nb_act1

      mole%ActiveTransfo%nb_inact2n  = mole%nb_inact2n
      mole%ActiveTransfo%nb_inact21  = mole%nb_inact21
      mole%ActiveTransfo%nb_inact22  = mole%nb_inact22
      mole%ActiveTransfo%nb_inact20  = mole%nb_inact20
      mole%ActiveTransfo%nb_inact    = mole%nb_inact

      mole%ActiveTransfo%nb_inact31  = mole%nb_inact31

      mole%ActiveTransfo%nb_rigid0   = mole%nb_rigid0
      mole%ActiveTransfo%nb_rigid100 = mole%nb_rigid100
      mole%ActiveTransfo%nb_rigid    = mole%nb_rigid

      IF (debug) write(out_unitp,*) 'END Sub_CoordType_TO_paraRPH'


      END SUBROUTINE Sub_CoordType_TO_paraRPH

      SUBROUTINE Sub_CoordType_TO_paraRPH_new(mole)

      TYPE (CoordType) :: mole

      integer :: i
      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.

      IF (.NOT. associated(mole%RPHTransfo)) RETURN
      IF (mole%RPHTransfo%option /= 0) RETURN
      ! This transformation is done only if option==0.
      !    For option=1, everything is already done

      IF (debug) write(out_unitp,*) 'BEGINNING Sub_CoordType_TO_paraRPH_new'


      mole%RPHTransfo%list_act_OF_Qdyn = mole%ActiveTransfo%list_act_OF_Qdyn

      DO i=1,mole%nb_var
        IF (mole%ActiveTransfo%list_act_OF_Qdyn(i) == 21)               &
                              mole%ActiveTransfo%list_act_OF_Qdyn(i) = 1
      END DO

      IF (debug) write(out_unitp,*) 'END Sub_CoordType_TO_paraRPH_new'


      END SUBROUTINE Sub_CoordType_TO_paraRPH_new

      ! this subroutine eanbles to perform automatic contraction
  SUBROUTINE CoordTypeRPH_TO_CoordTypeFlex(mole)
  USE mod_FlexibleTransfo, only : Read_FlexibleTransfo
  IMPLICIT NONE

  TYPE (CoordType), intent(inout) :: mole

      integer :: nb_Qtransfo,nb_Qin
      integer, allocatable :: list_flex(:)
      character (len=Name_len) :: name_transfo


!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = "CoordTypeRPH_TO_CoordTypeFlex"
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.



      IF (.NOT. associated(mole%RPHTransfo)) RETURN
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
       END IF
!-----------------------------------------------------------

      IF (mole%itRPH == -1) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' mole%itRPH = -1 is impossible for the RPH transformation'
        write(out_unitp,*) '  Check the fortran source!!'
        STOP
      END IF
      IF (debug) write(out_unitp,*) 'The RPH transformation is at:',mole%itRPH

      nb_Qin = mole%tab_Qtransfo(mole%itRPH)%nb_Qin


      ! modification of the transformation
      IF (debug) write(out_unitp,*) 'mole%RPHTransfo%list_act_OF_Qdyn', &
                                    mole%RPHTransfo%list_act_OF_Qdyn(:)
      CALL alloc_NParray(list_flex,(/nb_Qin/),"list_flex",name_sub)
      list_flex = mole%RPHTransfo%list_act_OF_Qdyn(:)
      WHERE (list_flex == 21)
        list_flex = 20
      END WHERE
      IF (debug) write(out_unitp,*) 'list_flex',list_flex(:)
      mole%tab_Qtransfo(mole%itRPH)%name_transfo = 'flexible'

      CALL Read_FlexibleTransfo(mole%tab_Qtransfo(mole%itRPH)%FlexibleTransfo,  &
                                nb_Qin,list_flex)

      CALL dealloc_NParray(list_flex,"list_flex",name_sub)

      CALL sub_Type_Name_OF_Qin(mole%tab_Qtransfo(mole%itRPH),"Qflex")

      ! remove the RPH transformation and the pointer
      CALL dealloc_RPHTransfo(mole%tab_Qtransfo(mole%itRPH)%RPHTransfo)
      nullify(mole%RPHTransfo)
      mole%itRPH = -1

       IF (debug) THEN
         write(out_unitp,*) 'END ',name_sub
       END IF

  END SUBROUTINE CoordTypeRPH_TO_CoordTypeFlex


  SUBROUTINE Set_OptimizationPara_FROM_CoordType(mole,Set_Val)
  USE mod_system
  IMPLICIT NONE

  !----- for the CoordType and Tnum --------------------------------------
  TYPE (CoordType), intent(in) :: mole
  integer,          intent(in) :: Set_Val

      integer :: i,i1,i2,it,itone,nopt,iQdyn
!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug = .FALSE.
!      logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'Set_OptimizationPara_FROM_CoordType'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
      END IF
!---------------------------------------------------------------------

      IF (mole%opt_param > 0) THEN

        DO it=1,mole%nb_Qtransfo
          IF (mole%tab_Qtransfo(it)%opt_param > 0) THEN
            IF (associated(mole%tab_Qtransfo(it)%oneDTransfo)) THEN

              DO itone=1,mole%tab_Qtransfo(it)%nb_transfo
                nopt = count(mole%tab_Qtransfo(it)%oneDTransfo(itone)%opt_cte /= 0)
                i1 = para_FOR_optimization%i_OptParam+1
                i2 = para_FOR_optimization%i_OptParam+nopt
                IF (debug) write(out_unitp,*) 'it,itone,nopt',it,itone,nopt
                IF (nopt > 0) THEN
                  para_FOR_optimization%nb_OptParam =                   &
                                para_FOR_optimization%nb_OptParam + nopt

                  IF (Set_Val == -1) THEN
                    DO i=1,nopt
                      IF (mole%tab_Qtransfo(it)%oneDTransfo(itone)%opt_cte(i) /= 0) &
                        para_FOR_optimization%Val_RVec(i1+i-1) =        &
                            mole%tab_Qtransfo(it)%oneDTransfo(itone)%cte(i)
                    END DO
                  ELSE IF (Set_Val == 1) THEN
                    DO i=1,nopt
                      IF (mole%tab_Qtransfo(it)%oneDTransfo(itone)%opt_cte(i) /= 0) &
                        mole%tab_Qtransfo(it)%oneDTransfo(itone)%cte(i) = &
                                    para_FOR_optimization%Val_RVec(i1+i-1)
                    END DO
                  END IF
                  para_FOR_optimization%i_OptParam = i2
                END IF
              END DO
            ELSE
              STOP 'problem with opt_para !!'
            END IF
          END IF
        END DO

      END IF

      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'nb_OptParam ',para_FOR_optimization%nb_OptParam
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF

  END SUBROUTINE Set_OptimizationPara_FROM_CoordType
END MODULE mod_Tnum
