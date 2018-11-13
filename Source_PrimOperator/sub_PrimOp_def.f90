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

   MODULE mod_PrimOp_def
   USE mod_system
   USE mod_OTF_def, only: param_OTF
   use mod_nDFit,   only: param_ndfit
   IMPLICIT NONE
     PRIVATE
       TYPE param_PES
          logical           :: init_done            = .FALSE. ! to check if the initalization is done
          real (kind=Rkind) :: pot0                 = ZERO    ! minimum of the PES
          real (kind=Rkind) :: min_pot              = HUGE(ONE) ! minimum of the PES
          real (kind=Rkind) :: max_pot              = -HUGE(ONE) ! minimum of the PES

          integer           :: nb_elec              = 1       ! nb of electronic PES

          logical           :: calc_scalar_Op       = .FALSE.

          logical           :: opt                  = .FALSE. ! If we minimize the PES (not used yet)
          logical           :: pot_cplx             = .FALSE. ! complex PES
          logical           :: OnTheFly             = .FALSE. ! On-the-fly calculation f PES and dipole
          logical           :: Read_OnTheFly_only   = .FALSE. ! Read On-the-fly calculation f PES and dipole only
                                                              ! without ab initio calculation
          logical           :: pot_act              = .TRUE.  ! active variables for the PES
          logical           :: pot_cart             = .FALSE. ! the potential is in cartesian coordinates
          logical           :: HarD                 = .TRUE.  ! Harmonic Domain for the PES
          logical           :: deriv_WITH_FiniteDiff= .FALSE. ! IF true, force to use finite difference scheme
          real (kind=Rkind) :: stepOp               = ONETENTH**2
          integer           :: pot_itQtransfo       = -1      ! for new Qtransfo (default nb_QTransfo, dyn. coordinates)
          integer           :: nb_scalar_Op         = 0       ! nb of Operator
          integer           :: Type_HamilOp         = 1       ! 1:  F2.d^2 + F1.d^1 + V
                                                              ! 10: d^1 G d^1 + V
          logical           :: direct_KEO           = .FALSE. ! IF true the metric tensor is recomputed

          ! For on-the-fly calculations
          TYPE (param_OTF)              :: para_OTF
          TYPE (param_OTF)              :: para_OTF_Dip
          logical                       :: levelEne_EQ_levelDip = .TRUE.

          ! For the nDfit
          logical                       :: nDfit_Op = .FALSE.
          TYPE (param_nDFit)            :: para_nDFit_V
          TYPE (param_nDFit), pointer   :: para_nDFit_Scalar_Op(:) => null()

        END TYPE param_PES

     PUBLIC param_PES, write_param_PES, Sub_PES_FromTnum_TO_PES

   CONTAINS

      SUBROUTINE write_param_PES(para_PES)
      USE mod_OTF_def, only: write_OTF
      USE mod_nDFit,   only: write_ndfit

      TYPE (param_PES) :: para_PES

      integer :: i



      write(out_unitp,*) ' BEGINNING write_param_PES'


      write(out_unitp,*) 'para_PES%init_done',para_PES%init_done
      write(out_unitp,*) 'para_PES%stepOp',para_PES%stepOp

      write(out_unitp,*) 'para_PES%opt',para_PES%opt

      write(out_unitp,*) 'para_PES%pot0',para_PES%pot0
      write(out_unitp,*) 'para_PES%min_pot',para_PES%min_pot
      write(out_unitp,*) 'para_PES%max_pot',para_PES%max_pot
      write(out_unitp,*) 'para_PES%HarD',para_PES%HarD
      write(out_unitp,*) 'para_PES%pot_cart',para_PES%pot_cart
      write(out_unitp,*) 'para_PES%nb_elec',para_PES%nb_elec
      write(out_unitp,*) 'para_PES%pot_cplx',para_PES%pot_cplx
      write(out_unitp,*) 'para_PES%OnTheFly',para_PES%OnTheFly
      write(out_unitp,*) 'para_PES%pot_itQtransfo',para_PES%pot_itQtransfo

      write(out_unitp,*) 'para_PES%calc_scalar_Op',para_PES%calc_scalar_Op
      write(out_unitp,*) 'para_PES%nb_scalar_Op',para_PES%nb_scalar_Op

      write(out_unitp,*) 'para_PES%deriv_WITH_FiniteDiff',para_PES%deriv_WITH_FiniteDiff
      write(out_unitp,*) 'para_PES%nDfit_Op',para_PES%nDfit_Op


      write(out_unitp,*) 'para_PES%Type_HamilOp',para_PES%Type_HamilOp

      write(out_unitp,*)
      IF (para_PES%OnTheFly) THEN
        write(out_unitp,*) 'para_PES%levelEne_EQ_levelDip',para_PES%levelEne_EQ_levelDip

        write(out_unitp,*) 'for PES'
        CALL write_OTF(para_PES%para_OTF)

        write(out_unitp,*) 'for saclar operators (dipole ...)'
        CALL write_OTF(para_PES%para_OTF_Dip)

      END IF

      IF (para_PES%nDfit_Op) THEN
        write(out_unitp,*) 'Write para_nDFit_V'
        CALL Write_nDFit(para_PES%para_nDFit_V)
        IF (associated(para_PES%para_nDFit_Scalar_Op)) THEN
          DO i=1,size(para_PES%para_nDFit_Scalar_Op)
            write(out_unitp,*) 'Write para_nDFit_Scalar_Op',i
            CALL Write_nDFit(para_PES%para_nDFit_Scalar_Op(i))
          END DO
        END IF
      END IF


      write(out_unitp,*) ' END write_param_PES'
      CALL flush_perso(out_unitp)
      END SUBROUTINE write_param_PES

      SUBROUTINE Sub_PES_FromTnum_TO_PES(para_PES,para_PES_FromTnum)
      USE mod_OTF_def,   only: init_G03_OTF
      USE mod_Coord_KEO, only: param_PES_FromTnum
      IMPLICIT NONE

      TYPE (param_PES)          :: para_PES
      TYPE (param_PES_FromTnum) :: para_PES_FromTnum

      integer :: i
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Sub_PES_FromTnum_TO_PES'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         CALL flush_perso(out_unitp)
       END IF
!-----------------------------------------------------------

        para_PES%stepOp                = para_PES_FromTnum%stepOp
        para_PES%opt                   = para_PES_FromTnum%opt
        para_PES%pot0                  = para_PES_FromTnum%pot0
        para_PES%pot_act               = para_PES_FromTnum%pot_act
        para_PES%pot_cart              = para_PES_FromTnum%pot_cart
        para_PES%HarD                  = para_PES_FromTnum%HarD
        para_PES%nb_elec               = para_PES_FromTnum%nb_elec
        para_PES%pot_cplx              = para_PES_FromTnum%pot_cplx
        para_PES%OnTheFly              = para_PES_FromTnum%OnTheFly
        para_PES%pot_itQtransfo        = para_PES_FromTnum%pot_itQtransfo
        para_PES%nb_scalar_Op          = para_PES_FromTnum%nb_scalar_Op
        para_PES%deriv_WITH_FiniteDiff = para_PES_FromTnum%deriv_WITH_FiniteDiff
        para_PES%nDfit_Op              = para_PES_FromTnum%nDfit_Op

      !IF (para_PES%OnTheFly) THEN
        para_PES%para_OTF%charge       = para_PES_FromTnum%charge
        para_PES%para_OTF%multiplicity = para_PES_FromTnum%multiplicity
        IF (para_PES%OnTheFly) THEN
          para_PES%nb_scalar_Op          = 3
        END IF


        CALL init_G03_OTF(para_PES%para_OTF)
        CALL init_G03_OTF(para_PES%para_OTF_DIP)

        SELECT CASE (para_PES_FromTnum%ab_initio_prog)
        CASE ('g03','g98','G98','G03','gaussian')
          para_PES%para_OTF%file_data%name=trim(para_PES_FromTnum%file_name_OTF)   // '.com'
          para_PES%para_OTF%file_log%name=trim(para_PES_FromTnum%file_name_OTF)    // '.log'
          para_PES%para_OTF%file_header%name=trim(para_PES_FromTnum%file_name_OTF) // '.header'
          para_PES%para_OTF%file_footer%name=trim(para_PES_FromTnum%file_name_OTF) // '.footer'
          IF (len_trim(para_PES_FromTnum%file_name_fchk) == 0) THEN
            para_PES%para_OTF%file_FChk%name = 'Test.FChk'
          ELSE
            para_PES%para_OTF%file_FChk%name = para_PES_FromTnum%file_name_fchk
          END IF
        CASE ('gamess','GAMESS')
          para_PES%para_OTF%file_data%name=trim(para_PES_FromTnum%file_name_OTF)   // '.inpg'
          para_PES%para_OTF%file_log%name=trim(para_PES_FromTnum%file_name_OTF)    // '.outg'
          para_PES%para_OTF%file_header%name=trim(para_PES_FromTnum%file_name_OTF) // '.header'
          para_PES%para_OTF%file_footer%name=trim(para_PES_FromTnum%file_name_OTF) // '.footer'
          para_PES%para_OTF%file_pun%name=trim(para_PES_FromTnum%file_name_OTF)    // '.dat'
        CASE ('gamess2014')
          para_PES%para_OTF%file_data%name=trim(para_PES_FromTnum%file_name_OTF)   // '.inp'
          para_PES%para_OTF%file_log%name=trim(para_PES_FromTnum%file_name_OTF)    // '.log'
          para_PES%para_OTF%file_header%name=trim(para_PES_FromTnum%file_name_OTF) // '.header'
          para_PES%para_OTF%file_footer%name=trim(para_PES_FromTnum%file_name_OTF) // '.footer'
          para_PES%para_OTF%file_pun%name=trim(para_PES_FromTnum%file_name_OTF)    // '.pun'
        CASE ('generic')
          para_PES%para_OTF%file_data%name=trim(para_PES_FromTnum%file_name_OTF)   // '.evrti'
          para_PES%para_OTF%file_log%name=trim(para_PES_FromTnum%file_name_OTF)    // '.evrto'
          para_PES%para_OTF%file_header%name=trim(para_PES_FromTnum%file_name_OTF) // '.header'
          para_PES%para_OTF%file_footer%name=trim(para_PES_FromTnum%file_name_OTF) // '.footer'
          para_PES%para_OTF%file_pun%name=trim(para_PES_FromTnum%file_name_OTF)    // '.dat'
        CASE ('g09','G09') ! particular case because the keyword formchk is obsolet
          para_PES%para_OTF%file_data%name=trim(para_PES_FromTnum%file_name_OTF)   // '.com'
          para_PES%para_OTF%file_log%name=trim(para_PES_FromTnum%file_name_OTF)    // '.log'
          para_PES%para_OTF%file_header%name=trim(para_PES_FromTnum%file_name_OTF) // '.header'
          para_PES%para_OTF%file_footer%name=trim(para_PES_FromTnum%file_name_OTF) // '.footer'
          IF (len_trim(para_PES_FromTnum%file_name_fchk) == 0) THEN
           para_PES%para_OTF%file_FChk%name = trim(para_PES_FromTnum%file_name_OTF)// '.fchk'
          ELSE
            para_PES%para_OTF%file_FChk%name = para_PES_FromTnum%file_name_fchk
          END IF
        CASE default ! ERROR: wrong program !
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' The ab initio program is UNKNOWN ',para_PES_FromTnum%ab_initio_prog
          STOP
        END SELECT
        IF (print_level > 1) THEN
          write(out_unitp,*) ' Files for the OTF'
          write(out_unitp,*) trim(para_PES%para_OTF%file_data%name)
          write(out_unitp,*) trim(para_PES%para_OTF%file_log%name)
          write(out_unitp,*) trim(para_PES%para_OTF%file_Fchk%name)
          write(out_unitp,*) trim(para_PES%para_OTF%file_pun%name)
          write(out_unitp,*) ' Program for the OTF ',para_PES_FromTnum%ab_initio_prog
          write(out_unitp,*) ' Unix script for the OTF ',para_PES_FromTnum%commande_unix
        END IF
        para_PES%para_OTF%header          = para_PES_FromTnum%header
        para_PES%para_OTF%footer          = para_PES_FromTnum%footer
        para_PES%para_OTF%file_name       = para_PES_FromTnum%file_name_OTF
        para_PES%para_OTF%ab_initio_prog  = para_PES_FromTnum%ab_initio_prog
        para_PES%para_OTF%commande_unix   = para_PES_FromTnum%commande_unix

        para_PES%para_OTF_DIP = para_PES%para_OTF

        para_PES%para_OTF%ab_initio_meth      = para_PES_FromTnum%ab_initio_methEne
        para_PES%para_OTF%ab_initio_basis     = para_PES_FromTnum%ab_initio_basisEne
        para_PES%para_OTF_Dip%ab_initio_meth  = para_PES_FromTnum%ab_initio_methDip
        para_PES%para_OTF_Dip%ab_initio_basis = para_PES_FromTnum%ab_initio_basisDip

        para_PES%levelEne_EQ_levelDip         =                         &
          (para_PES_FromTnum%ab_initio_methEne  .EQ. para_PES_FromTnum%ab_initio_methDip) .AND.&
          (para_PES_FromTnum%ab_initio_basisEne .EQ. para_PES_FromTnum%ab_initio_basisDip)

      !END IF

      para_PES%nDfit_Op = para_PES_FromTnum%nDfit_Op
      IF (para_PES%nDfit_Op) THEN
        para_PES%para_nDFit_V%Param_Fit_file%name =                     &
                      para_PES_FromTnum%para_nDFit_V%Param_Fit_file%name
        IF (associated(para_PES_FromTnum%para_nDFit_Scalar_Op)) THEN

          allocate(para_PES%para_nDFit_Scalar_Op(para_PES%nb_scalar_Op))
          DO i=1,para_PES%nb_scalar_Op
            para_PES%para_nDFit_Scalar_Op(i)%Param_Fit_file%name =      &
                para_PES_FromTnum%para_nDFit_Scalar_Op(i)%Param_Fit_file%name
          END DO

        END IF
      END IF

      IF (debug) THEN
        CALL write_param_PES(para_PES)
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF

      END SUBROUTINE Sub_PES_FromTnum_TO_PES

   END MODULE mod_PrimOp_def

