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
   MODULE mod_PrimOp_def
   USE mod_system
   USE mod_OTF_def, only: param_OTF, dealloc_OTF
   use mod_nDFit,   only: param_nDFit, dealloc_nDFit
   use mod_CAP
   use mod_HStep
   IMPLICIT NONE
   PRIVATE

     TYPE PrimOp_t
        real (kind=Rkind) :: pot0                 = ZERO       ! reference energy of the PES (given)
        real (kind=Rkind) :: pot_Qref             = ZERO       ! Energy at the reference geometry (calculated)
        real (kind=Rkind) :: min_pot              =  HUGE(ONE) ! minimum of the PES
        real (kind=Rkind) :: max_pot              = -HUGE(ONE) ! minimum of the PES

        integer           :: nb_elec              = 1       ! nb of electronic PES

        logical           :: calc_scalar_Op       = .FALSE.

        logical           :: opt                  = .FALSE. ! If we minimize the PES (not used yet)
        logical           :: pot_cplx             = .FALSE. ! complex PES

        logical           :: Read_OnTheFly_only   = .FALSE. ! Read On-the-fly calculation of PES and dipole only
                                                            ! without ab initio calculation
        logical           :: HarD                 = .TRUE.  ! Harmonic Domain for the PES
        logical           :: deriv_WITH_FiniteDiff= .FALSE. ! IF true, force to use finite difference scheme
        real (kind=Rkind) :: stepOp               = ONETENTH**2
        integer           :: pot_itQtransfo       = -1      ! for new Qtransfo (default nb_QTransfo, dyn. coordinates)
        integer           :: nb_scalar_Op         = 0       ! nb of Operators
        integer           :: nb_CAP               = 0       ! nb of CAP Operators
        integer           :: nb_FluxOp            = 0       ! nb of flux Operators

        integer           :: Type_HamilOp         = 1       ! 1:  F2.d^2 + F1.d^1 + V
                                                            ! 10: d^1 G d^1 + V
        logical           :: direct_KEO           = .FALSE. ! IF true, the metric tensor is recomputed
        logical           :: direct_ScalOp        = .FALSE. ! IF true, the scalar op is recomputed

        ! For CAP: complex absorbing potential (when nb_CAP > 0)
        TYPE (CAP_t), allocatable     :: tab_CAP(:)
        ! For HStep: Heaviside step function (when nb_FluxOp > 0)
        TYPE (HStep_t), allocatable     :: tab_HStep(:)

        ! For on-the-fly calculations
        logical                       :: OnTheFly             = .FALSE. ! On-the-fly calculation of PES and dipole
        TYPE (param_OTF)              :: para_OTF
        TYPE (param_OTF)              :: para_OTF_Dip
        logical                       :: levelEne_EQ_levelDip = .TRUE.

        ! parameters for the Quantum Model Lib (ECAM), KEO+PES
        logical                       :: QMLib   = .FALSE.
        logical                       :: QMLib_G = .FALSE.
        integer, allocatable          :: Qit_TO_QQMLib(:)


        ! For the nDfit
        logical                         :: nDfit_Op = .FALSE.
        TYPE (param_nDFit)              :: para_nDFit_V
        TYPE (param_nDFit), allocatable :: para_nDFit_Scalar_Op(:)
      CONTAINS
        PROCEDURE, PRIVATE, PASS(PrimOp1) :: PrimOp2_TO_PrimOp1
        GENERIC,   PUBLIC  :: assignment(=) => PrimOp2_TO_PrimOp1
      END TYPE PrimOp_t

      TYPE, EXTENDS(PrimOp_t) :: param_PES
       ! nothing, just to keep the name
      END TYPE param_PES

    PUBLIC :: param_PES, write_param_PES ! for PVSCF
    PUBLIC :: PrimOp_t, write_PrimOp, dealloc_PrimOp, Sub_PES_FromTnum_TO_PrimOp

  CONTAINS

  SUBROUTINE write_param_PES(para_PES)
  IMPLICIT NONE
  TYPE (param_PES) :: para_PES
    CALL write_PrimOp(para_PES%PrimOp_t)
  END SUBROUTINE write_param_PES

  SUBROUTINE write_PrimOp(PrimOp)
  USE mod_OTF_def, only: write_OTF
  USE mod_nDFit,   only: write_ndfit
  IMPLICIT NONE

      TYPE (PrimOp_t) :: PrimOp

      integer :: i

      write(out_unitp,*) ' BEGINNING write_PrimOp (write_param_PES)'


      write(out_unitp,*) 'PrimOp%stepOp',PrimOp%stepOp

      write(out_unitp,*) 'PrimOp%opt',PrimOp%opt

      write(out_unitp,*) 'PrimOp%pot0',PrimOp%pot0
      write(out_unitp,*) 'PrimOp%pot_Qref',PrimOp%pot_Qref
      write(out_unitp,*) 'PrimOp%min_pot',PrimOp%min_pot
      write(out_unitp,*) 'PrimOp%max_pot',PrimOp%max_pot
      write(out_unitp,*) 'PrimOp%HarD',PrimOp%HarD
      write(out_unitp,*) 'PrimOp%nb_elec',PrimOp%nb_elec
      write(out_unitp,*) 'PrimOp%pot_cplx',PrimOp%pot_cplx
      write(out_unitp,*) 'PrimOp%OnTheFly',PrimOp%OnTheFly
      write(out_unitp,*) 'PrimOp%pot_itQtransfo',PrimOp%pot_itQtransfo

      write(out_unitp,*) 'PrimOp%calc_scalar_Op',PrimOp%calc_scalar_Op
      write(out_unitp,*) 'PrimOp%nb_scalar_Op',PrimOp%nb_scalar_Op

      write(out_unitp,*) 'PrimOp%nb_CAP',PrimOp%nb_CAP
      IF (allocated(PrimOp%tab_CAP)) then
        DO i=1,size(PrimOp%tab_CAP)
          CALL Write_CAP(PrimOp%tab_CAP(i))
        END DO
      END IF

      write(out_unitp,*) 'PrimOp%nb_FluxOp',PrimOp%nb_FluxOp
      IF (allocated(PrimOp%tab_HStep)) then
        DO i=1,size(PrimOp%tab_HStep)
          CALL Write_HStep(PrimOp%tab_HStep(i))
        END DO
      END IF

      write(out_unitp,*) 'PrimOp%deriv_WITH_FiniteDiff',PrimOp%deriv_WITH_FiniteDiff
      write(out_unitp,*) 'PrimOp%nDfit_Op',PrimOp%nDfit_Op

      ! parameters for the Quantum Model Lib (ECAM), KEO+PES
      write(out_unitp,*) 'PrimOp%QMLib',PrimOp%QMLib
      IF (allocated(PrimOp%Qit_TO_QQMLib)) THEN
        write(out_unitp,*) 'PrimOp%Qit_TO_QQMLib',PrimOp%Qit_TO_QQMLib
      END IF

      write(out_unitp,*) 'PrimOp%Type_HamilOp',PrimOp%Type_HamilOp
      write(out_unitp,*) 'PrimOp%direct_KEO',PrimOp%direct_KEO
      write(out_unitp,*) 'PrimOp%direct_ScalOp',PrimOp%direct_ScalOp

      write(out_unitp,*)
      IF (PrimOp%OnTheFly) THEN
        write(out_unitp,*) 'PrimOp%levelEne_EQ_levelDip',PrimOp%levelEne_EQ_levelDip

        write(out_unitp,*) 'for PES'
        CALL write_OTF(PrimOp%para_OTF)

        write(out_unitp,*) 'for saclar operators (dipole ...)'
        CALL write_OTF(PrimOp%para_OTF_Dip)

      END IF

      IF (PrimOp%nDfit_Op) THEN
        write(out_unitp,*) 'Write para_nDFit_V'
        CALL Write_nDFit(PrimOp%para_nDFit_V)
        IF (allocated(PrimOp%para_nDFit_Scalar_Op)) THEN
          DO i=1,size(PrimOp%para_nDFit_Scalar_Op)
            write(out_unitp,*) 'Write para_nDFit_Scalar_Op',i
            CALL Write_nDFit(PrimOp%para_nDFit_Scalar_Op(i))
          END DO
        END IF
      END IF


    write(out_unitp,*) ' END write_PrimOp'
    CALL flush_perso(out_unitp)
  END SUBROUTINE write_PrimOp
  SUBROUTINE PrimOp2_TO_PrimOp1(PrimOp1,PrimOp2)
  IMPLICIT NONE
      CLASS (PrimOp_t), intent(inout) :: PrimOp1
      TYPE (PrimOp_t),  intent(in)    :: PrimOp2

      integer :: i

      !write(out_unitp,*) ' BEGINNING PrimOp2_TO_PrimOp1'

      PrimOp1%pot0                  = PrimOp2%pot0
      PrimOp1%pot_Qref              = PrimOp2%pot_Qref
      PrimOp1%min_pot               = PrimOp2%min_pot
      PrimOp1%max_pot               = PrimOp2%max_pot
      PrimOp1%nb_elec               = PrimOp2%nb_elec
      PrimOp1%calc_scalar_Op        = PrimOp2%calc_scalar_Op
      PrimOp1%opt                   = PrimOp2%opt
      PrimOp1%pot_cplx              = PrimOp2%pot_cplx
      PrimOp1%OnTheFly              = PrimOp2%OnTheFly
      PrimOp1%Read_OnTheFly_only    = PrimOp2%Read_OnTheFly_only
      PrimOp1%HarD                  = PrimOp2%HarD
      PrimOp1%deriv_WITH_FiniteDiff = PrimOp2%deriv_WITH_FiniteDiff
      PrimOp1%stepOp                = PrimOp2%stepOp
      PrimOp1%pot_itQtransfo        = PrimOp2%pot_itQtransfo
      PrimOp1%nb_scalar_Op          = PrimOp2%nb_scalar_Op
      PrimOp1%Type_HamilOp          = PrimOp2%Type_HamilOp
      PrimOp1%direct_KEO            = PrimOp2%direct_KEO
      PrimOp1%direct_ScalOp         = PrimOp2%direct_ScalOp

      PrimOp1%para_OTF              = PrimOp2%para_OTF
      PrimOp1%para_OTF_Dip          = PrimOp2%para_OTF_Dip
      PrimOp1%levelEne_EQ_levelDip  = PrimOp2%levelEne_EQ_levelDip

      PrimOp1%QMLib                 = PrimOp2%QMLib
      IF (allocated(PrimOp2%Qit_TO_QQMLib)) THEN
        PrimOp1%Qit_TO_QQMLib         = PrimOp2%Qit_TO_QQMLib
      END IF

      PrimOp1%nDfit_Op              = PrimOp2%nDfit_Op
      PrimOp1%para_nDFit_V          = PrimOp2%para_nDFit_V
      IF (allocated(PrimOp2%para_nDFit_Scalar_Op)) THEN
        allocate(PrimOp1%para_nDFit_Scalar_Op(size(PrimOp2%para_nDFit_Scalar_Op)))
        DO i=1,size(PrimOp2%para_nDFit_Scalar_Op)
          PrimOp1%para_nDFit_Scalar_Op(i) = PrimOp2%para_nDFit_Scalar_Op(i)
        END DO
      END IF

      PrimOp1%nb_CAP              = PrimOp2%nb_CAP
      IF (allocated(PrimOp1%tab_CAP)) STOP 'PrimOp2_TO_PrimOp1: tab_CAP is allocated'
      IF (allocated(PrimOp2%tab_CAP)) THEN
        allocate(PrimOp1%tab_CAP(size(PrimOp2%tab_CAP)))
        DO i=1,size(PrimOp2%tab_CAP)
          PrimOp1%tab_CAP(i) = PrimOp2%tab_CAP(i)
        END DO
      END IF

      PrimOp1%nb_FluxOp           = PrimOp2%nb_FluxOp
      IF (allocated(PrimOp1%tab_HStep)) STOP 'PrimOp2_TO_PrimOp1: tab_HStep is allocated'
      IF (allocated(PrimOp2%tab_HStep)) THEN
        allocate(PrimOp1%tab_HStep(size(PrimOp2%tab_HStep)))
        DO i=1,size(PrimOp2%tab_HStep)
          PrimOp1%tab_HStep(i) = PrimOp2%tab_HStep(i)
        END DO
      END IF

     !write(out_unitp,*) ' END PrimOp2_TO_PrimOp1'
     !CALL flush_perso(out_unitp)
  END SUBROUTINE PrimOp2_TO_PrimOp1
  SUBROUTINE dealloc_PrimOp(PrimOp)
  IMPLICIT NONE
      CLASS (PrimOp_t), intent(inout) :: PrimOp

      integer :: i

      !write(out_unitp,*) ' BEGINNING dealloc_PrimOp'

      PrimOp%pot0                  = ZERO
      PrimOp%pot_Qref              = ZERO
      PrimOp%min_pot               = HUGE(ONE)
      PrimOp%max_pot               = -HUGE(ONE)
      PrimOp%nb_elec               = 1
      PrimOp%calc_scalar_Op        = .FALSE.
      PrimOp%opt                   = .FALSE.
      PrimOp%pot_cplx              = .FALSE.
      PrimOp%OnTheFly              = .FALSE.
      PrimOp%Read_OnTheFly_only    = .FALSE.
      PrimOp%HarD                  = .TRUE.
      PrimOp%deriv_WITH_FiniteDiff = .FALSE.
      PrimOp%stepOp                = ONETENTH**2
      PrimOp%pot_itQtransfo        = -1
      PrimOp%nb_scalar_Op          = 0
      PrimOp%nb_CAP                = 0
      PrimOp%nb_FluxOp             = 0
      PrimOp%Type_HamilOp          = 1
      PrimOp%direct_KEO            = .FALSE.
      PrimOp%direct_ScalOp         = .FALSE.

      IF (allocated(PrimOp%tab_CAP)) then
        DO i=1,size(PrimOp%tab_CAP)
          CALL dealloc_CAP(PrimOp%tab_CAP(i))
        END DO
        deallocate(PrimOp%tab_CAP)
      END IF

      IF (allocated(PrimOp%tab_HStep)) then
        DO i=1,size(PrimOp%tab_HStep)
          CALL dealloc_HStep(PrimOp%tab_HStep(i))
        END DO
        deallocate(PrimOp%tab_HStep)
      END IF

      CALL dealloc_OTF(PrimOp%para_OTF)
      CALL dealloc_OTF(PrimOp%para_OTF_Dip)

      PrimOp%levelEne_EQ_levelDip  =  .TRUE.

      PrimOp%QMLib                 = .FALSE.
      IF (allocated(PrimOp%Qit_TO_QQMLib)) deallocate(PrimOp%Qit_TO_QQMLib)

      PrimOp%nDfit_Op              = .FALSE.
      CALL dealloc_nDFit(PrimOp%para_nDFit_V)
      IF (allocated(PrimOp%para_nDFit_Scalar_Op)) THEN
        DO i=1,size(PrimOp%para_nDFit_Scalar_Op)
          CALL dealloc_nDFit(PrimOp%para_nDFit_Scalar_Op(i))
        END DO
        deallocate(PrimOp%para_nDFit_Scalar_Op)
      END IF

     !write(out_unitp,*) ' END dealloc_PrimOp'
     !CALL flush_perso(out_unitp)
  END SUBROUTINE dealloc_PrimOp
  SUBROUTINE Sub_PES_FromTnum_TO_PrimOp(PrimOp,para_PES_FromTnum)
      USE mod_OTF_def,   only: init_OTF
      USE mod_Coord_KEO, only: param_PES_FromTnum
      IMPLICIT NONE

      TYPE (PrimOp_t)           :: PrimOp
      TYPE (param_PES_FromTnum) :: para_PES_FromTnum

      character (len=Name_len)      :: ab_initio_prog
      integer :: i
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Sub_PES_FromTnum_TO_PrimOp'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         CALL flush_perso(out_unitp)
       END IF
!-----------------------------------------------------------

        PrimOp%stepOp                = para_PES_FromTnum%stepOp
        PrimOp%opt                   = para_PES_FromTnum%opt
        PrimOp%pot0                  = para_PES_FromTnum%pot0
        PrimOp%HarD                  = para_PES_FromTnum%HarD
        PrimOp%nb_elec               = para_PES_FromTnum%nb_elec
        PrimOp%pot_cplx              = para_PES_FromTnum%pot_cplx
        PrimOp%OnTheFly              = para_PES_FromTnum%OnTheFly
        PrimOp%pot_itQtransfo        = para_PES_FromTnum%pot_itQtransfo
        PrimOp%nb_scalar_Op          = para_PES_FromTnum%nb_scalar_Op
        PrimOp%nb_CAP                = para_PES_FromTnum%nb_CAP

        PrimOp%deriv_WITH_FiniteDiff = para_PES_FromTnum%deriv_WITH_FiniteDiff
        PrimOp%nDfit_Op              = para_PES_FromTnum%nDfit_Op

        PrimOp%QMLib                  = para_PES_FromTnum%QMLib
        PrimOp%QMLib_G                = para_PES_FromTnum%QMLib_G

        PrimOp%para_OTF%charge       = para_PES_FromTnum%charge
        PrimOp%para_OTF%multiplicity = para_PES_FromTnum%multiplicity
        IF (PrimOp%OnTheFly) THEN
          PrimOp%nb_scalar_Op          = 3
        END IF

        CALL init_OTF(PrimOp%para_OTF)
        CALL init_OTF(PrimOp%para_OTF_DIP)

        ab_initio_prog = para_PES_FromTnum%ab_initio_prog
        CALL string_uppercase_TO_lowercase(ab_initio_prog)

        SELECT CASE (ab_initio_prog)
        CASE ('g03','g98','gaussian')
          PrimOp%para_OTF%file_data%name=trim(para_PES_FromTnum%file_name_OTF)   // '.com'
          PrimOp%para_OTF%file_log%name=trim(para_PES_FromTnum%file_name_OTF)    // '.log'
          PrimOp%para_OTF%file_header%name=trim(para_PES_FromTnum%file_name_OTF) // '.header'
          PrimOp%para_OTF%file_footer%name=trim(para_PES_FromTnum%file_name_OTF) // '.footer'
          IF (len_trim(para_PES_FromTnum%file_name_fchk) == 0) THEN
            PrimOp%para_OTF%file_FChk%name = 'Test.FChk'
          ELSE
            PrimOp%para_OTF%file_FChk%name = para_PES_FromTnum%file_name_fchk
          END IF
        CASE ('gamess')
          PrimOp%para_OTF%file_data%name=trim(para_PES_FromTnum%file_name_OTF)   // '.inpg'
          PrimOp%para_OTF%file_log%name=trim(para_PES_FromTnum%file_name_OTF)    // '.outg'
          PrimOp%para_OTF%file_header%name=trim(para_PES_FromTnum%file_name_OTF) // '.header'
          PrimOp%para_OTF%file_footer%name=trim(para_PES_FromTnum%file_name_OTF) // '.footer'
          PrimOp%para_OTF%file_pun%name=trim(para_PES_FromTnum%file_name_OTF)    // '.dat'
        CASE ('gamess2014')
          PrimOp%para_OTF%file_data%name=trim(para_PES_FromTnum%file_name_OTF)   // '.inp'
          PrimOp%para_OTF%file_log%name=trim(para_PES_FromTnum%file_name_OTF)    // '.log'
          PrimOp%para_OTF%file_header%name=trim(para_PES_FromTnum%file_name_OTF) // '.header'
          PrimOp%para_OTF%file_footer%name=trim(para_PES_FromTnum%file_name_OTF) // '.footer'
          PrimOp%para_OTF%file_pun%name=trim(para_PES_FromTnum%file_name_OTF)    // '.pun'
        CASE ('molpro')
          PrimOp%para_OTF%file_data%name=trim(para_PES_FromTnum%file_name_OTF)   // '.inp'
          PrimOp%para_OTF%file_log%name=trim(para_PES_FromTnum%file_name_OTF)    // '.outm'
          PrimOp%para_OTF%file_header%name=trim(para_PES_FromTnum%file_name_OTF) // '.header'
          PrimOp%para_OTF%file_footer%name=trim(para_PES_FromTnum%file_name_OTF) // '.footer'
          PrimOp%para_OTF%file_pun%name=trim(para_PES_FromTnum%file_name_OTF)    // '.log'
        CASE ('generic')
          PrimOp%para_OTF%file_data%name=trim(para_PES_FromTnum%file_name_OTF)   // '.evrti'
          PrimOp%para_OTF%file_log%name=trim(para_PES_FromTnum%file_name_OTF)    // '.evrto'
          PrimOp%para_OTF%file_header%name=trim(para_PES_FromTnum%file_name_OTF) // '.header'
          PrimOp%para_OTF%file_footer%name=trim(para_PES_FromTnum%file_name_OTF) // '.footer'
          PrimOp%para_OTF%file_pun%name=trim(para_PES_FromTnum%file_name_OTF)    // '.dat'
        CASE ('g09') ! particular case because the keyword formchk is obsolet
          PrimOp%para_OTF%file_data%name=trim(para_PES_FromTnum%file_name_OTF)   // '.com'
          PrimOp%para_OTF%file_log%name=trim(para_PES_FromTnum%file_name_OTF)    // '.log'
          PrimOp%para_OTF%file_header%name=trim(para_PES_FromTnum%file_name_OTF) // '.header'
          PrimOp%para_OTF%file_footer%name=trim(para_PES_FromTnum%file_name_OTF) // '.footer'
          IF (len_trim(para_PES_FromTnum%file_name_fchk) == 0) THEN
           PrimOp%para_OTF%file_FChk%name = trim(para_PES_FromTnum%file_name_OTF)// '.fchk'
          ELSE
            PrimOp%para_OTF%file_FChk%name = para_PES_FromTnum%file_name_fchk
          END IF
        CASE default ! ERROR: wrong program !
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' The ab initio program is UNKNOWN ',para_PES_FromTnum%ab_initio_prog
          STOP
        END SELECT
        IF (print_level > 1) THEN
          write(out_unitp,*) ' Files for the OTF'
          write(out_unitp,*) trim(PrimOp%para_OTF%file_data%name)
          write(out_unitp,*) trim(PrimOp%para_OTF%file_log%name)
          write(out_unitp,*) trim(PrimOp%para_OTF%file_Fchk%name)
          write(out_unitp,*) trim(PrimOp%para_OTF%file_pun%name)
          write(out_unitp,*) ' Program for the OTF ',para_PES_FromTnum%ab_initio_prog
          write(out_unitp,*) ' Unix script for the OTF ',para_PES_FromTnum%commande_unix
        END IF

        PrimOp%para_OTF%header          = para_PES_FromTnum%header
        PrimOp%para_OTF%footer          = para_PES_FromTnum%footer
        PrimOp%para_OTF%file_name       = para_PES_FromTnum%file_name_OTF
        PrimOp%para_OTF%ab_initio_prog  = para_PES_FromTnum%ab_initio_prog
        PrimOp%para_OTF%commande_unix   = para_PES_FromTnum%commande_unix

        PrimOp%para_OTF_DIP                 = PrimOp%para_OTF

        PrimOp%para_OTF%ab_initio_meth      = para_PES_FromTnum%ab_initio_methEne
        PrimOp%para_OTF%ab_initio_basis     = para_PES_FromTnum%ab_initio_basisEne
        PrimOp%para_OTF_Dip%ab_initio_meth  = para_PES_FromTnum%ab_initio_methDip
        PrimOp%para_OTF_Dip%ab_initio_basis = para_PES_FromTnum%ab_initio_basisDip

        PrimOp%levelEne_EQ_levelDip         =                         &
          (para_PES_FromTnum%ab_initio_methEne  .EQ. para_PES_FromTnum%ab_initio_methDip) .AND.&
          (para_PES_FromTnum%ab_initio_basisEne .EQ. para_PES_FromTnum%ab_initio_basisDip)


      PrimOp%nDfit_Op = para_PES_FromTnum%nDfit_Op
      IF (PrimOp%nDfit_Op) THEN
        PrimOp%para_nDFit_V%Param_Fit_file%name =                       &
                                      para_PES_FromTnum%nDFit_V_name_Fit

        PrimOp%para_nDFit_V%name_Fit = para_PES_FromTnum%nDFit_V_name_Fit

        IF (allocated(para_PES_FromTnum%nDFit_Scalar_Op_name_Fit)) THEN

          allocate(PrimOp%para_nDFit_Scalar_Op(PrimOp%nb_scalar_Op))
          DO i=1,PrimOp%nb_scalar_Op
            PrimOp%para_nDFit_Scalar_Op(i)%Param_Fit_file%name =        &
                           para_PES_FromTnum%nDFit_Scalar_Op_name_Fit(i)
          END DO

          deallocate(para_PES_FromTnum%nDFit_Scalar_Op_name_Fit) ! we don't need anymore.

        END IF
      END IF


      IF (debug) THEN
        CALL write_PrimOp(PrimOp)
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF

  END SUBROUTINE Sub_PES_FromTnum_TO_PrimOp

  END MODULE mod_PrimOp_def
