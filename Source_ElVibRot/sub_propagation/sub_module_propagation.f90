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

      MODULE mod_propa
      USE mod_system
      USE mod_Constant,      only: get_conv_au_to_unit, real_wu, convrwu_to_r
      USE mod_field,         ONLY : param_field
      USE mod_psi_set_alloc, ONLY : param_psi
      USE mod_param_WP0,     ONLY : param_WP0
      USE mod_type_ana_psi,  ONLY : param_ana_psi
      IMPLICIT NONE

PRIVATE
PUBLIC :: param_poly,param_control,param_Davidson,param_propa
PUBLIC :: read_propagation,read_davidson
PUBLIC :: dealloc_param_propa,sub_analyze_WP_OpWP,sub_analyze_mini_WP_OpWP
PUBLIC :: Read_AutoCorr,Write_AutoCorr,Calc_AutoCorr
PUBLIC :: SaveWP_restart,ReadWP_restart
PUBLIC :: initialisation1_poly,cof

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
        TYPE param_poly

        logical           :: init_done = .FALSE.

        real (kind=Rkind) :: DHmax,Hmax,Hmin        ! parameters
        real (kind=Rkind) :: deltaE,Esc,E0          ! parameters
        real (kind=Rkind) :: alpha                  ! parameters
        real (kind=Rkind) :: poly_tol
        integer           :: max_poly               ! max number of polynomials
        real (kind=Rkind), pointer :: coef_poly(:) => null()  !
        integer           :: npoly                  ! number polynomials
        integer           :: npoly_Opt              ! optimal number polynomials

        END TYPE param_poly
        TYPE param_control

        logical          :: restart            ! restart

        integer          :: nb_WP              ! number of initial WP and WP target's
        integer          :: nb_WPba            ! number of basis function of each WP
        integer          :: max_iter           ! number of iterations for the control
        real (kind=Rkind)    :: conv               ! convergence for the objective
        real (kind=Rkind)    :: alpha,gamma        ! parameters for the convergence
        real (kind=Rkind)    :: Max_alpha          ! parameters for the convergence
        logical          :: Krotov             ! Krotov Algorithm
        logical          :: Turinici           ! Turinici Algorithm
        logical          :: envelopp           ! use an envelopp (s(t)) for the field
        real (kind=Rkind)    :: Tenvelopp          ! parameter of the envelopp sin(t/T*pi)
        logical          :: Obj_TO_alpha       ! the objectifs modify alpha

        integer, pointer :: tab_WP0(:) => null()         ! numero of the intial WP (nb_WP)
        integer, pointer :: tab_WPt(:) => null()        ! numero of the target WP (nb_WP)

        logical          :: lo_WP_save_T
        integer          :: nb_DT
        integer, pointer :: tab_WP_save_T(:,:) => null() ! (nb_WP,0:nb_DT) IF enough memory, save the time-dependent WP


        logical          :: post_control       ! T => one propagation for
                                               ! the WP analysis
        logical          :: gate               ! T => read a matrix for the
                                               ! definition of the gate
        logical          :: cplx_gate          ! T => if Mgate0 and Mgatet are complex
        complex (kind=Rkind), pointer :: Mgate0(:,:) => null()  ! gate Matrix (for WP0)
        complex (kind=Rkind), pointer :: Mgatet(:,:) => null()  ! gate Matrix (fort WPt)

        END TYPE param_control
        TYPE param_Davidson
         integer :: num_resetH,num_checkS
         logical :: one_by_one
         integer :: residual_max_nb
         integer :: max_it,nb_WP,max_WP,num_LowestWP,nb_WP0

         logical :: read_WP,read_listWP,precond,formatted_file_readWP
         integer :: nb_readWP,nb_readWP_OF_List
         character (len=Line_len) :: name_file_readWP
         real (kind=Rkind) :: precond_tol

         logical :: save_all
         real (kind=Rkind) :: save_max_ene,scaled_max_ene
         integer :: save_max_nb
         character (len=Line_len) :: name_file_saveWP
         logical :: formatted_file_WP

         logical :: all_lower_states,lower_states,project_WP0,Op_Transfo
         logical :: Hmin_propa,Hmax_propa
         real (kind=Rkind) :: thresh_project = 0.8_Rkind

         real (kind=Rkind) :: Max_ene
         integer           :: symab


         real (kind=Rkind) :: RMS_ene,conv_ene
         real (kind=Rkind) :: RMS_resi,conv_resi
         integer :: conv_hermitian
         integer :: NewVec_type = 2
         logical :: With_Grid   = .FALSE.

         ! parameters specific to the filter diagonalization
         real (kind=Rkind) :: E0_filter     = ZERO   ! the center of the windowin ua (but read in cm-1)
         real (kind=Rkind) :: W_filter      = ZERO   ! the center of the windowin ua (but read in cm-1)
         real (kind=Rkind) :: LambdaMin     = ZERO   ! The lower energy of the window
         real (kind=Rkind) :: LambdaMax     = ZERO   ! The higher energy of the window
         integer           :: L_filter      = 50     ! the size of filter basis
         integer           :: Lmax_filter   = 1000   ! the maximal size of filter basis
         integer           :: M_filter      = 100    ! the number of terms in the Chebyshev expansion
         integer           :: DeltaM_filter = 0      ! the number of terms in the Chebyshev expansion
         integer           :: Mmax_filter   = 10000  ! the number of terms in the Chebyshev expansion

        END TYPE param_Davidson

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
        TYPE param_propa


        real (kind=Rkind)   ::  Hmax,Hmin     !
        logical             ::  once_control_Hmin=.TRUE.!< control the calculation of 
                                                        !! Hmin once at the first action
        logical             ::  auto_Hmax     !  .TRUE. => Hmax is obtained with a propagation
                                              !            with imaginary time (type_WPpropa=-3)
                                              !            (default .FALSE.)
        real (kind=Rkind)   ::  WPTmax        ! propagation time in au
        real (kind=Rkind)   ::  WPT0          ! initial time in au
        logical             ::  restart       ! restart the propagation

        real (kind=Rkind)   ::  WPdeltaT      ! time step in au
        integer             ::  nb_micro      ! nb of microiterations

        integer             ::  max_ana       ! nb of wavefunctions calculated with imaginary propagation

        integer             ::  type_WPpropa  ! propagation type :
                                              !  0 from the spectrum (not working)
                                              !  1 Chebychef (default)
                                              !  2 nOD
                                              !  3 imaginary propagation
                                              ! 22 Time dependant
                                              ! 221,222,223 Time dependant (with Mux or Muy or Muz)
                                              ! 24 Time dependant (for the control)
                                              ! 241,242,243 Time dependant (with Mux or Muy or Muz)

        logical                  :: With_field = .FALSE.   ! propagation with a field
        character (len=Name_len) :: name_WPpropa  !  spectral (not working)
                                              !  cheby (default)
                                              !  nOD
                                              !  Hmin,Hmax (imaginary propagation: nOD)
                                              !  TDnOD

        logical         ::       spectral     ! Def (f). If T => spectral propagatio


        integer         ::       type_corr    !  kind of correlation function
                                              !  default 0 => <psi0(T)|psi(T)>
                                              !  1 =>  psi(T,i_qa_corr,...) (for a grid point)
        integer         ::       i_qa_corr    !  grid point for type_corr = 1
        integer         ::       channel_ie_corr !   Channel for type_corr = 1
        integer         ::       Op_corr      !  defined the Operator when type_corr=2



        integer         ::       n_WPecri        !  write WP every n_WPecri time steps
        logical         ::       WPpsi           !  write WP
        logical         ::       WPpsi2          !  write WP2
        logical         ::       write_iter      !  write iteration
        logical         ::       write_GridRep   !  write WP of the grid
        logical         ::       write_BasisRep  !  write WP of the basis
        logical         ::       write_WPAdia    !  if T, write WP on the adiabatic PES (otherwise on the diabatic ones)

        TYPE (param_ana_psi)  :: ana_psi         ! to control the WP analysis


        logical         ::       control         !  use control (need type_WPpropa 24)
        logical         ::       test_max_norm   ! IF .TRUE., Error in the propagation (norm too large)
        logical         ::       march_error     ! IF .TRUE., Error in the propagation
        integer         ::       num_Op = 0      !  Operator used for the propagation (H)

!----- variables for the WP propagation ----------------------------
        TYPE (param_psi), pointer :: work_WP(:) => null()

        TYPE (param_file)    :: file_WP,file_autocorr,file_spectrum
        TYPE (param_file)    :: file_WP_restart ! files to be able to write the propagation

        TYPE(param_WP0)      :: para_WP0
        TYPE(param_control)  :: para_control
        TYPE(param_Davidson) :: para_Davidson
        TYPE(param_poly)     :: para_poly
        TYPE(param_field)    :: para_field

        integer           :: TFnexp2
        real (kind=Rkind) :: TFmaxE,TFminE

        END TYPE param_propa

      CONTAINS
!================================================================
!
!    alloc / dealloc param_poly
!
!================================================================
      !!@description: alloc param_poly
      !!@param: para_poly
      SUBROUTINE alloc_param_poly(para_poly)
      IMPLICIT NONE
      TYPE (param_poly), intent(inout) :: para_poly

        IF (para_poly%max_poly >= 0 .AND.                               &
            .NOT. associated(para_poly%coef_poly)) THEN

          CALL alloc_array(para_poly%coef_poly,(/para_poly%max_poly/),  &
                          "para_poly%coef_poly","alloc_param_poly")
        END IF

      END SUBROUTINE alloc_param_poly

      !!@description: dealloc param_poly
      !!@param: para_poly
      SUBROUTINE dealloc_param_poly(para_poly)
      IMPLICIT NONE
      TYPE (param_poly), intent(inout) :: para_poly

        IF (associated(para_poly%coef_poly)) THEN
          CALL dealloc_array(para_poly%coef_poly,                       &
                          "para_poly%coef_poly","alloc_param_poly")
        END IF

      END SUBROUTINE dealloc_param_poly

!================================================================
!
!    alloc / dealloc param_control
!
!================================================================
      !!@description: dealloc param_control
      !!@param: para_control
      SUBROUTINE dealloc_param_control(para_control)
      IMPLICIT NONE
      TYPE (param_control), intent(inout) :: para_control

        IF (associated(para_control%tab_WP0)) THEN
          CALL dealloc_array(para_control%tab_WP0,                      &
                            "para_control%tab_WP0","dealloc_param_control")
        END IF

        IF (associated(para_control%tab_WPt)) THEN
          CALL dealloc_array(para_control%tab_WPt,                      &
                            "para_control%tab_WPt","dealloc_param_control")
        END IF

        IF (associated(para_control%tab_WP_save_T)) THEN
          CALL dealloc_array(para_control%tab_WP_save_T,                &
                            "para_control%tab_WP_save_T","dealloc_param_control")
        END IF

        IF (associated(para_control%Mgate0)) THEN
          CALL dealloc_array(para_control%Mgate0,                       &
                            "para_control%Mgate0","dealloc_param_control")
        END IF

        IF (associated(para_control%Mgatet)) THEN
          CALL dealloc_array(para_control%Mgatet,                       &
                            "para_control%Mgatet","dealloc_param_control")
        END IF

      END SUBROUTINE dealloc_param_control


!================================================================
!
!    alloc / dealloc param_propa
!
!================================================================
      SUBROUTINE dealloc_param_propa(para_propa)
      USE mod_system
      USE mod_field,         ONLY : dealloc_param_field
      USE mod_param_WP0,     ONLY : dealloc_param_WP0
      USE mod_psi_set_alloc, ONLY : dealloc_psi,dealloc_array
      USE mod_type_ana_psi,  ONLY : dealloc_ana_psi
      IMPLICIT NONE
      TYPE (param_propa), intent(inout) :: para_propa

      integer :: i

      CALL dealloc_ana_psi(para_propa%ana_psi)

      CALL dealloc_param_WP0(para_propa%para_WP0)
      CALL dealloc_param_control(para_propa%para_control)
      CALL dealloc_param_poly(para_propa%para_poly)
      CALL dealloc_param_field(para_propa%para_field)

      IF (associated(para_propa%work_WP)) THEN
        DO i=0,size(para_propa%work_WP)-1
          CALL dealloc_psi(para_propa%work_WP(i))
        END DO
        CALL dealloc_array(para_propa%work_WP,                          &
                          "para_propa%work_WP","dealloc_param_propa")
      END IF

      END SUBROUTINE dealloc_param_propa

      SUBROUTINE SaveWP_restart(T,WP,file_restart)
      USE mod_system
      USE mod_psi_set_alloc
      USE mod_file
      IMPLICIT NONE


!----- variables for the WP propagation ----------------------------
      TYPE (param_file), intent(inout) :: file_restart
      TYPE (param_psi),  intent(in)    :: WP(:)
      real (kind=Rkind), intent(in)    :: T


!------ working parameters --------------------------------
      integer       :: i,no_restart


!----- for debuging --------------------------------------------------
      character (len=*), parameter ::name_sub='SaveWP_restart'
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) ' nb_psi',size(WP)
        write(out_unitp,*)
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------


      CALL file_open(file_restart,no_restart)
      write(no_restart,*) T
      IF(MPI_id==0) THEN
        DO i=1,size(WP)
          write(no_restart,*) WP(i)%CvecB
        END DO
      ENDIF
      close(no_restart)

!----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
       END IF
!----------------------------------------------------------

      END SUBROUTINE SaveWP_restart
      SUBROUTINE ReadWP_restart(T,WP,file_restart)
      USE mod_system
      USE mod_psi_set_alloc
      USE mod_file
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_file), intent(inout) :: file_restart
      TYPE (param_psi),  intent(inout) :: WP(:)
      real (kind=Rkind), intent(inout) :: T

!------ working parameters --------------------------------
      integer       :: i,no_restart,err_read

!----- for debuging --------------------------------------------------
      character (len=*), parameter ::name_sub='ReadWP_restart'
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) ' nb_psi',size(WP)
        write(out_unitp,*)
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

      CALL file_open(file_restart,no_restart,old=.TRUE.)

      T = ZERO
      !err_read = 0
      read(no_restart,*,IOSTAT=err_read) T
      IF (err_read /= 0) THEN
        write(out_unitp,*) ' WARNING in ',name_sub
        write(out_unitp,*) ' T (time) is not present in the restart file'
        write(out_unitp,*) ' file name:',trim(file_restart%name)
        write(out_unitp,*) ' => No restart, T0=0'
        T = ZERO
      ELSE
        write(out_unitp,*) 'T0 for the restart:',T
        DO i=1,size(WP)
          read(no_restart,*) WP(i)%CvecB
        END DO
      END IF

      close(no_restart)

!----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
       END IF
!----------------------------------------------------------

      END SUBROUTINE ReadWP_restart
!==============================================================
!
!     Calculation of the autocorrelation function
!
!==============================================================

      FUNCTION Calc_AutoCorr(psi0,psi,para_propa,T,Write_AC)
      USE mod_system
      USE mod_psi_set_alloc,   ONLY : param_psi,ecri_psi
      USE mod_psi_B_TO_G,      ONLY : sub_PsiBasisRep_TO_GridRep
      USE mod_psi_Op,          ONLY : Overlap_psi1_psi2
      USE mod_MPI
      IMPLICIT NONE


!----- variables for the WP propagation ----------------------------
      TYPE (param_propa)   :: para_propa
      TYPE (param_psi)     :: psi0,psi

      real (kind=Rkind)    :: T
      logical, optional    :: Write_AC
      complex (kind=Rkind) :: Calc_AutoCorr

      complex (kind=Rkind) :: cdot
      logical              :: Write_AC_loc
      integer              :: i_qaie_corr


      IF (present(Write_AC)) THEN
        Write_AC_loc = Write_AC
      ELSE
        Write_AC_loc = .FALSE.
      END IF


      IF (para_propa%type_corr == 0) THEN
        CALL Overlap_psi1_psi2(cdot,psi0,psi)
      ELSE IF (para_propa%type_corr == 1) THEN
        CALL sub_PsiBasisRep_TO_GridRep(psi)
        i_qaie_corr = para_propa%i_qa_corr +                            &
                    (para_propa%channel_ie_corr-1)*psi%nb_qa
        cdot = psi%CvecG(i_qaie_corr)
      ELSE IF (para_propa%type_corr == 2) THEN
        STOP 'test new auto'
      ELSE
        STOP 'This autocorrelation procedure is not defined'
      END IF

      IF (Write_AC_loc) THEN
        CALL Write_AutoCorr(para_propa%file_autocorr%unit,T,cdot)
      END IF

      Calc_AutoCorr = cdot

      END FUNCTION Calc_AutoCorr

!==============================================================
!
!     write/read the autocorrelation function
!
!==============================================================
      SUBROUTINE Write_AutoCorr(no,T,cdot)
      USE mod_system
      IMPLICIT NONE

      integer              :: no
      real (kind=Rkind)    :: T
      complex (kind=Rkind) :: cdot

      write(no,*) 'AutoCor ',T,real(cdot,kind=Rkind),aimag(cdot),abs(cdot)

      END SUBROUTINE Write_AutoCorr
      SUBROUTINE Read_AutoCorr(no,T,cdot)
      USE mod_system
      IMPLICIT NONE

      integer                :: no,err_io
      real (kind=Rkind)      :: T
      complex (kind=Rkind)   :: cdot

      character (len=Name_len) :: name
      real (kind=Rkind)      :: a,b,c

      read(no,*,IOSTAT=err_io) name,T,a,b,c
      IF (err_io == 0) THEN
       cdot = cmplx(a,b,kind=Rkind)
      ELSE
       cdot = czero
      END IF
!     write(out_unitp,*) name,T,cdot

      END SUBROUTINE Read_AutoCorr


!==============================================================
!
!     initialisation for the Chebychev and nOD propagation
!
!==============================================================
      SUBROUTINE initialisation1_poly(para_poly,deltaT,type_propa)
      USE mod_system
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_poly) :: para_poly
      real (kind=Rkind) :: deltaT
      integer           :: type_propa



      integer       :: npoly_cheby

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='initialisation1_poly'
      logical, parameter :: debug =.FALSE.
      !logical, parameter :: debug =.TRUE.
!-----------------------------------------------------------
       IF (para_poly%init_done) RETURN
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'deltaT',deltaT
         write(out_unitp,*) 'Hmax',para_poly%Hmax
         write(out_unitp,*) 'DHmax',para_poly%DHmax
         write(out_unitp,*) 'Hmin',para_poly%Hmin
         CALL flush_perso(out_unitp)
       END IF
!-----------------------------------------------------------

      IF (para_poly%Hmin > para_poly%Hmax) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) 'Hmax',para_poly%Hmax
         write(out_unitp,*) 'DHmax',para_poly%DHmax
         write(out_unitp,*) 'Hmin',para_poly%Hmin
         write(out_unitp,*) ' Hmin > Hmax '
         STOP " ERROR in " // name_sub // " : Hmin > Hmax"
      END IF

!-----------------------------------------------------------
!     E0 is the zero of energy for energy scaling (i.e. the centre of the range).
!     ESC is the energy scaling parameter ( so that scaled energy lies between -1 and 1)

      para_poly%deltaE = para_poly%Hmax - para_poly%Hmin
      para_poly%E0     = para_poly%Hmin + HALF * para_poly%deltaE
      para_poly%alpha  = HALF * para_poly%deltaE * deltaT
      para_poly%Esc    = ONE

      IF (type_propa == 1) para_poly%Esc = HALF * para_poly%deltaE
      IF (type_propa == 33) para_poly%E0     = para_poly%Hmin
      IF (type_propa == 33) para_poly%alpha  = para_poly%deltaE *       &
                                                  deltaT


      IF (debug) THEN
        write(out_unitp,*) 'Hmin,Hmax',para_poly%Hmin,para_poly%Hmax
        write(out_unitp,*) 'deltaE',para_poly%deltaE
        write(out_unitp,*) 'E0,Esc',para_poly%E0,para_poly%Esc
        write(out_unitp,*) 'alpha (r)',para_poly%alpha
        CALL flush_perso(out_unitp)
      END IF
!     ------------------------------------------------------

!     ------------------------------------------------------
      IF (type_propa == 1) THEN
!       - Chebychev coeficients calculation ------------------
        CALL cofcheb(para_poly%alpha,npoly_cheby,                       &
                     para_poly%coef_poly,                               &
                     para_poly%max_poly,para_poly%poly_tol)
        IF (npoly_cheby > para_poly%npoly) para_poly%npoly = npoly_cheby

      ELSE
        CALL cof_nOD(para_poly%alpha,DeltaT,para_poly%npoly,            &
                     para_poly%coef_poly,                               &
                     para_poly%max_poly,para_poly%poly_tol)
      END IF
!     ------------------------------------------------------

      para_poly%init_done = .TRUE.

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'npoly',para_poly%npoly
         write(out_unitp,*) 'deltaT',deltaT
         write(out_unitp,*) 'END ',name_sub
       END IF
!-----------------------------------------------------------
       CALL flush_perso(out_unitp)

       END SUBROUTINE initialisation1_poly


!==============================================================
!
!     Chebychev coeficients calculation
!     the number of polynomia ncheb is optimized
!
!==============================================================
      !!@description: Calculation of the Chebychev coefficients
      !!              the number of polynomia ncheb is optimized
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      SUBROUTINE cofcheb(r,ncheb,cf,nchmx,tol)
      USE mod_system
      IMPLICIT NONE

      integer       :: ncheb,nchmx
      real (kind=Rkind) :: cf(:)
      real (kind=Rkind) :: r,tol,ctst,chk
      integer       :: i,nset,k6,ier


!----- for debuging --------------------------------------------------
      logical, parameter :: debug =.FALSE.
      !logical, parameter :: debug =.TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING cofcheb'
         write(out_unitp,*) 'r,tol',r,tol
         write(out_unitp,*) 'nchmx',nchmx
         CALL flush_perso(out_unitp)
       END IF
!-----------------------------------------------------------

!----------------------------------------------
      ncheb = max(ONE,r)
      IF (debug) write(out_unitp,*) 'r, ncheb : ',r,ncheb
      IF(MPI_id==0) CALL flush_perso(out_unitp)

      IF (ncheb > nchmx) THEN
         write(out_unitp,*) ' ERROR in cofcheb'
         write(out_unitp,*) '*** nchmx too small '
         write(out_unitp,*) 'nchmx < ncheb',nchmx,ncheb
         STOP
      END IF

      CALL cof(r,ncheb,cf)

!     -----------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Chebychev coeficients',ncheb
        write(out_unitp,2009) (i,cf(i),i=1,ncheb)
        CALL flush_perso(out_unitp)
      END IF
!     -----------------------------------------
!----------------------------------------------

!----------------------------------------------
! Set nch by using tolerance parameter, tol

  30  ctst = abs(cf(ncheb)/cf(1))
      if(ctst > tol)then
        ncheb=ncheb+10
        IF (ncheb > nchmx) THEN
          write(out_unitp,*) ' ERROR in cofcheb'
          write(out_unitp,*) '*** nchmx too small '
          write(out_unitp,*) 'nchmx < ncheb',nchmx,ncheb
          STOP
        END IF
        call cof(r,ncheb,cf)
      go to 30
      else
        nset=ncheb
        DO  k6=1,nset
            chk=abs(cf(k6)/cf(1))
            if(chk <= tol)then
               ncheb=k6
               go to 45
            endif
        END DO
      endif
   45 continue
!----------------------------------------------

      write(out_unitp,2007) ncheb
 2007 format(1x,'Ncheb after optimization = ',i6)


!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Chebychev coeficients',ncheb
        write(out_unitp,2009) (i,cf(i),i=1,ncheb)
 2009   format(3(i6,e16.5))
        write(out_unitp,*) 'END cofcheb'
      END IF
!-----------------------------------------------------------
      IF(MPI_id==0) CALL flush_perso(out_unitp)

      END SUBROUTINE cofcheb
!==============================================================
!
!     Chebychev coeficients calculation
!
!==============================================================
      !!@description: Chebychev coefficients calculation
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      SUBROUTINE cof(r,ncheb,cf)
      USE mod_system
      IMPLICIT NONE

      integer ncheb
      real (kind=Rkind) :: cf(:)
      real (kind=Rkind) :: r,r1

      integer       :: i,ier

      r1 = r
      CALL mmbsjn(r1,ncheb,cf,ier)

      !write(6,*) 'Chebychev coefficients'
      !write(6,*) '  with r=',r1

      DO i=2,ncheb
        cf(i) = cf(i) + cf(i)
        !write(6,*) 'i,cf',i,cf(i)
      END DO
!stop
      END SUBROUTINE cof

!==============================================================
!
!     Taylor coeficients calculation (nOD)
!
!==============================================================
      !!@description: Taylor coeficients calculation (nOD)
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      SUBROUTINE cof_nOD(alpha,DeltaT,nOD,coef,Max_nOD,epsi)
      USE mod_system
      IMPLICIT NONE

      integer       :: nOD,Max_nOD
      real (kind=Rkind) :: alpha,DeltaT,epsi
      real (kind=Rkind) :: coef(:)



      real (kind=Rkind) :: reste,xi
      integer       :: i

!----- for debuging --------------------------------------------------
      logical, parameter :: debug =.FALSE.
!     logical, parameter :: debug =.TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING cof_nOD'
         write(out_unitp,*) 'alpha,DeltaT,epsi',alpha,DeltaT,epsi
         write(out_unitp,*) 'nOD,Max_nOD',nOD,Max_nOD
       END IF
!-----------------------------------------------------------

!----------------------------------------------

      !- determination of the order: nOD
      IF (nOD < 1) THEN
        i = 1
        xi = alpha
        reste = ONE +alpha
        DO WHILE (abs(log(reste)-alpha) > epsi .AND. i <= Max_nOD)
          i = i+1
          xi = xi * alpha / real(i,kind=Rkind)
          reste = reste +xi
          !write(out_unitp,*) 'i,xi,reste',i,xi,log(reste)-alpha
        END DO
        nOD  = i
        write(out_unitp,*) 'Optimal nOD1',nOD,log(reste)-alpha
        DO WHILE (xi > epsi .AND. i <= Max_nOD)
          i = i+1
          xi = xi * alpha / real(i,kind=Rkind)
          !write(out_unitp,*) 'i,xi',i,xi
        END DO
        nOD  = i
        IF (debug) write(out_unitp,*) 'Optimal nOD2',nOD,xi
      ELSE
        IF (debug) write(out_unitp,*) 'Fixed   nOD ',nOD
      END IF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END cof_nOD'
      END IF
!-----------------------------------------------------------

      END SUBROUTINE cof_nOD

!=======================================================================================
SUBROUTINE sub_analyze_WP_OpWP(T,WP,nb_WP,para_H,para_propa,adia,para_field)
  USE mod_system
  USE mod_Op,              ONLY : param_Op,sub_PsiOpPsi,sub_PsiDia_TO_PsiAdia_WITH_MemGrid
  USE mod_field,           ONLY : param_field,sub_dnE
  USE mod_ExactFact

  USE mod_psi_set_alloc,   ONLY : param_psi,ecri_psi,alloc_psi,dealloc_psi
  USE mod_ana_psi,         ONLY : sub_analyze_psi,norm2_psi
  USE mod_psi_B_TO_G,      ONLY : sub_PsiBasisRep_TO_GridRep
  USE mod_psi_SimpleOp,    ONLY : operator (*),operator (+),operator (-),assignment (=)
  USE mod_MPI
  IMPLICIT NONE

  real (kind=Rkind) :: T      ! time

!----- variables pour la namelist minimum ----------------------------
  TYPE (param_Op)   :: para_H

!- variables for the WP propagation ----------------------------
  TYPE (param_propa) :: para_propa
  TYPE (param_field), optional :: para_field
  logical,            optional :: adia

  integer            :: nb_WP
  TYPE (param_psi)   :: WP(:)

  !logical            :: ana_mini = .TRUE.
  logical            :: ana_mini = .FALSE.

!-- working parameters --------------------------------
  TYPE (param_psi)   :: w1,w2

  integer       :: j,i,i_bi,i_be,i_bie
  complex (kind=Rkind) :: ET  ! energy
  character (len=30)   :: info
  logical :: BasisRep,GridRep,adia_loc,Write_psi2_Grid,Write_psi_Grid

!- for the field --------------------------------------------------
  real (kind=Rkind)    :: dnE(3)
!- for the field --------------------------------------------------

!- for debuging --------------------------------------------------
  character (len=*), parameter :: name_sub='sub_analyze_WP_OpWP'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
!-------------------------------------------------------
  IF (debug) THEN
   write(out_unitp,*) 'BEGINNING ',name_sub
   write(out_unitp,*) 'Tmax,deltaT',para_propa%WPTmax,para_propa%WPdeltaT
   write(out_unitp,*)
   write(out_unitp,*) 'nb_ba,nb_qa',WP(1)%nb_ba,WP(1)%nb_qa
   write(out_unitp,*) 'nb_bi',WP(1)%nb_bi
   write(out_unitp,*)

   DO i=1,nb_WP
     CALL norm2_psi(WP(i),.FALSE.,.TRUE.,.FALSE.)
     write(out_unitp,*) 'normepsi0 BasisRep',i,WP(1)%norme

     write(out_unitp,*) 'WP(i)%BasisRep',i
     CALL ecri_psi(T=ZERO,psi=WP(1),                               &
                   ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)
     write(out_unitp,*) 'WP(i)%GridRep',i
     CALL ecri_psi(T=ZERO,psi=WP(1),                               &
                   ecri_GridRep=.TRUE.,ecri_BasisRep=.FALSE.)
   END DO
  END IF
!-----------------------------------------------------------

  IF (ana_mini) THEN
    IF (present(adia)) THEN
      IF (present(para_field)) THEN
        CALL sub_analyze_mini_WP_OpWP(T,WP,nb_WP,para_H,para_propa,adia,para_field)
      ELSE
        CALL sub_analyze_mini_WP_OpWP(T,WP,nb_WP,para_H,para_propa,adia)
      END IF
    ELSE
      IF (present(para_field)) THEN
        CALL sub_analyze_mini_WP_OpWP(T,WP,nb_WP,para_H,para_propa,para_field=para_field)
      ELSE
        CALL sub_analyze_mini_WP_OpWP(T,WP,nb_WP,para_H,para_propa)
      END IF
    END IF
    RETURN
  END IF

  Write_psi2_Grid = para_propa%ana_psi%Write_psi2_Grid
  Write_psi_Grid  = para_propa%ana_psi%Write_psi_Grid

  IF(MPI_id==0) THEN
    BasisRep = WP(1)%BasisRep
    GridRep  = WP(1)%GridRep
  ENDIF

  para_propa%ana_psi%T = T

  IF (present(adia)) THEN
    adia_loc = adia
  ELSE
    adia_loc = para_propa%Write_WPAdia
  END IF

  IF (.NOT. para_H%para_ReadOp%para_FileGrid%Save_MemGrid_done) adia_loc = .FALSE.

  !-----------------------------------------------------------
  ! => the WPs on the Grid
  IF (.NOT. para_propa%ana_psi%GridDone) THEN
    CALL time_perso('sub_PsiBasisRep_TO_GridRep ini')
    DO i=1,nb_WP
      IF(MPI_id==0) CALL sub_PsiBasisRep_TO_GridRep(WP(i))
    END DO
    CALL time_perso('sub_PsiBasisRep_TO_GridRep end')
  END IF
  para_propa%ana_psi%GridDone = .TRUE.
  !-----------------------------------------------------------

  !-----------------------------------------------------------
   IF (present(para_field)) THEN
     CALL sub_dnE(dnE,0,T,para_field)
     para_propa%ana_psi%With_field = .TRUE.
     para_propa%ana_psi%field      = dnE
   ELSE
     para_propa%ana_psi%With_field = .FALSE.
     para_propa%ana_psi%field      = (/ZERO,ZERO,ZERO/)
   END IF
  !-----------------------------------------------------------

  !-----------------------------------------------------------
  w2 = WP(1) ! just for the initialization
  IF(MPI_id==0) CALL alloc_psi(w2,GridRep=.TRUE.)

  DO i=1,nb_WP
    !-----------------------------------------------------------
    ! => Analysis for diabatic potential (always done)

    ! =>first the energy
    IF(MPI_id==0) THEN
      w1 = WP(i)
      CALL norm2_psi(w1,.FALSE.,.TRUE.,.FALSE.)
    ENDIF
    CALL sub_PsiOpPsi(ET,w1,w2,para_H)
    
    IF(MPI_id==0) THEN
      WP(i)%CAvOp = ET/w1%norme

      para_propa%ana_psi%num_psi = i
      para_propa%ana_psi%Ene     = real(WP(i)%CAvOp,kind=Rkind)

      ! => The analysis (diabatic)
      para_propa%ana_psi%adia = .FALSE.
      para_propa%ana_psi%Write_psi2_Grid = Write_psi2_Grid
      para_propa%ana_psi%Write_psi_Grid  = Write_psi_Grid
      CALL sub_analyze_psi(WP(i),para_propa%ana_psi)

      IF (para_propa%ana_psi%ExactFact > 0) THEN
        write(out_unitp,*) i,'Exact Factorization analysis at ',T, ' ua'
        IF (present(para_field)) THEN
          CALL sub_ExactFact_analysis(T,WP(i),para_propa%ana_psi,para_H,  &
                         para_propa%WPTmax,para_propa%WPdeltaT,para_field)
        ELSE
          CALL sub_ExactFact_analysis(T,WP(i),para_propa%ana_psi,para_H,  &
                                    para_propa%WPTmax,para_propa%WPdeltaT)
        END IF
      END IF


!    ! => The analysis (adiabatic)
!    IF (adia_loc) THEN
!      w1 = WP(i)
!      para_propa%ana_psi%adia = .TRUE.
!      para_propa%ana_psi%Write_psi2_Grid = Write_psi2_Grid
!      para_propa%ana_psi%Write_psi_Grid  = .FALSE.
!      CALL sub_PsiDia_TO_PsiAdia_WITH_MemGrid(w1,para_H)
!      CALL sub_analyze_psi(w1,para_propa%ana_psi)
!    END IF
      CALL flush_perso(out_unitp)

      para_propa%ana_psi%adia = .FALSE.
      CALL alloc_psi(WP(i),BasisRep=BasisRep,GridRep=GridRep)
    ENDIF ! for MPI_id==0
  END DO

  CALL dealloc_psi(w1,delete_all=.TRUE.)
  CALL dealloc_psi(w2,delete_all=.TRUE.)


  para_propa%ana_psi%GridDone = .FALSE.

!----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'END ',name_sub
  END IF
!-------------------------------------------------------
END SUBROUTINE sub_analyze_WP_OpWP
!=======================================================================================

SUBROUTINE sub_analyze_mini_WP_OpWP(T,WP,nb_WP,para_H,para_propa,adia,para_field)
  USE mod_system
  USE mod_Op,              ONLY : param_Op,sub_PsiOpPsi,sub_PsiDia_TO_PsiAdia_WITH_MemGrid
  USE mod_field,           ONLY : param_field,sub_dnE

  USE mod_psi_set_alloc,   ONLY : param_psi,ecri_psi,alloc_psi,dealloc_psi
  USE mod_ana_psi,         ONLY : sub_analyze_psi,norm2_psi,Channel_weight
  USE mod_psi_B_TO_G,      ONLY : sub_PsiBasisRep_TO_GridRep
  USE mod_psi_SimpleOp,    ONLY : operator (*),operator (+),operator (-),assignment (=)
  USE mod_MPI
  IMPLICIT NONE

  real (kind=Rkind) :: T      ! time

!----- variables pour la namelist minimum ----------------------------
  TYPE (param_Op)   :: para_H

!- variables for the WP propagation ----------------------------
  TYPE (param_propa) :: para_propa
  TYPE (param_field), optional :: para_field
  logical,            optional :: adia

  integer            :: nb_WP
  TYPE (param_psi)   :: WP(:)

!-- working parameters --------------------------------
  TYPE (param_psi)   :: w2

  !integer       :: j,i,i_bi,i_be,i_bie

  integer       :: i,i_bi,i_be

  complex (kind=Rkind)              :: ET  ! energy
  real (kind=Rkind)                 :: E
  TYPE(REAL_WU)                     :: RWU_E
  character(len=:), allocatable     :: psi_line


  real (kind=Rkind), allocatable    :: tab_WeightChannels(:,:)
  real (kind=Rkind)                 :: Psi_norm2
!- for the field --------------------------------------------------
  real (kind=Rkind)    :: dnE(3)
!- for the field --------------------------------------------------

!- for debuging --------------------------------------------------
  character (len=*), parameter :: name_sub='sub_analyze_mini_WP_OpWP'
  logical, parameter :: debug=.FALSE.
! logical, parameter :: debug=.TRUE.
!-------------------------------------------------------
  IF (debug) THEN
   write(out_unitp,*) 'BEGINNING ',name_sub
   write(out_unitp,*) 'Tmax,deltaT',para_propa%WPTmax,para_propa%WPdeltaT
   write(out_unitp,*)
   write(out_unitp,*) 'nb_ba,nb_qa',WP(1)%nb_ba,WP(1)%nb_qa
   write(out_unitp,*) 'nb_bi',WP(1)%nb_bi
   write(out_unitp,*)

   DO i=1,nb_WP
     CALL norm2_psi(WP(i),.FALSE.,.TRUE.,.FALSE.)
     write(out_unitp,*) 'normepsi0 BasisRep',i,WP(1)%norme

     write(out_unitp,*) 'WP(i)%BasisRep',i
     CALL ecri_psi(T=ZERO,psi=WP(1),                               &
                   ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)
     write(out_unitp,*) 'WP(i)%GridRep',i
     CALL ecri_psi(T=ZERO,psi=WP(1),                               &
                   ecri_GridRep=.TRUE.,ecri_BasisRep=.FALSE.)
   END DO
  END IF
!-----------------------------------------------------------


  !-----------------------------------------------------------
  IF(MPI_id==0) w2 = WP(1) ! just for the initialization
  IF(MPI_id==0) CALL alloc_psi(w2)

  DO i=1,nb_WP
    !-----------------------------------------------------------
    ! => Analysis for diabatic potential (always done)

    ! =>first the energy
    IF(MPI_id==0) CALL norm2_psi(WP(i))
    CALL sub_PsiOpPsi(ET,WP(i),w2,para_H)
    
    IF(MPI_id==0) THEN
      WP(i)%CAvOp = ET/WP(i)%norme

      RWU_E  = REAL_WU(real(WP(i)%CAvOp,kind=Rkind),'au','E')
      E      = convRWU_TO_R(RWU_E ,WorkingUnit=.FALSE.)

      CALL Channel_weight(tab_WeightChannels,WP(i),GridRep=.FALSE.,BasisRep=.TRUE.)
      Psi_norm2 = sum(tab_WeightChannels)

      ! add the psi number + the time
      psi_line = 'norm^2-WP #WP ' // int_TO_char(i) // ' ' // real_TO_char(T,Rformat='f12.2')

      ! add the energy
      psi_line = psi_line // ' ' // real_TO_char(E,Rformat='f8.5')

      ! add the field (if necessary)
      IF (present(para_field)) THEN
        CALL sub_dnE(dnE,0,T,para_field)
        psi_line = psi_line // ' ' // real_TO_char(dnE(1),Rformat='f8.5')
        psi_line = psi_line // ' ' // real_TO_char(dnE(2),Rformat='f8.5')
        psi_line = psi_line // ' ' // real_TO_char(dnE(3),Rformat='f8.5')
      END IF

      psi_line = psi_line // ' ' // real_TO_char(Psi_norm2,Rformat='f10.7')
      DO i_be=1,WP(i)%nb_be
      DO i_bi=1,WP(i)%nb_bi
        psi_line = psi_line // ' ' // real_TO_char(tab_WeightChannels(i_bi,i_be),Rformat='f10.7')
      END DO
      END DO

      write(out_unitp,*) psi_line
    ENDIF ! for MPI_id==0

  END DO

  CALL dealloc_psi(w2,delete_all=.TRUE.)

  IF (allocated(tab_WeightChannels)) deallocate(tab_WeightChannels)
  IF (allocated(psi_line))           deallocate(psi_line)

!----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'END ',name_sub
  END IF
!-------------------------------------------------------
END SUBROUTINE sub_analyze_mini_WP_OpWP
!================================================================
! ++    read the namelist for the wavepacket propagation
!
! Tmax  : time of the propagation
! deltaT : step of time
!
! n_WPecri : write the WP every n_ecri*deltaT
!       WPpsi2 = .TRUE. => writing of psi*psi
!       WPpsi  = .TRUE. => writing of Im(psi) and Real(psi)
!
!       type_propa : type of propagation
!                    0 from the spectrum
!                    1 Chebychef (default)
!
!       type_Hpsi  : type of Hpsi
!                    0 built H in BasisRep (default)
!
!       definition of WP0 (not normalized):
!       for each variable Qi : exp[-((Q-Qeq)/sigma)2+i*imp_k*(Q-Qeq)]
!
!================================================================
      SUBROUTINE read_propagation(para_propa,nb_act1,nb_vp_spec_out)
      USE mod_system
      USE mod_psi_set_alloc, ONLY : alloc_array
      USE mod_param_WP0,     ONLY : alloc_param_WP0
      USE mod_Constant
      USE mod_MPI
      IMPLICIT NONE

!------ parameter for the propagation ---------------------
      TYPE (param_propa) :: para_propa
      integer            :: nb_act1,nb_vp_spec_out

!----- physical and mathematical constants ---------------------------
      real (kind=Rkind) :: auTOcm_inv


!----- time parameters -----------------------------------------------
      !real (kind=Rkind) :: WPTmax ,WPdeltaT
      TYPE (REAL_WU)    :: WPTmax ,WPdeltaT

      integer           :: nb_micro


!     - for propagation with polynomials -------------------
      logical           :: auto_Hmax
      real (kind=Rkind) :: DHmax
      integer           :: max_poly,npoly
      real (kind=Rkind) :: poly_tol


      logical                  :: spectral
      integer                  :: nb_vp_spec
      integer                  :: type_WPpropa
      character (len=Name_len) :: name_WPpropa

      integer            :: n_WPecri
      logical            :: WPpsi2 ,WPpsi,write_DVR,write_FBR,write_iter,write_WPAdia
      character (len=Line_len) :: file_autocorr,file_WP,file_spectrum,file_restart


      integer         ::       type_corr       !  kind of correlation function
                                               !  default 0 => <psi0(T)|psi(T)>
                                               !  1 =>  psi(T,i_qa_corr,...) (for a grid point)
      integer         ::       i_qa_corr       !  grid point for type_corr = 1
      integer         ::       channel_ie_corr !  Channel for type_corr = 1
      integer         ::       Op_corr         !  Op used with type_corr=2
      logical         ::       restart         !  restart the propagation


      logical      :: New_Read_WP0,read_listWP0,read_file
      logical      :: WP0restart,WP0cplx
      logical      :: lect_WP0DVR ,lect_WP0FBR,lect_WP0FBRall,WP0FBR
      integer      :: WP0n_h,WP0nb_elec,WP0_DIP,WP0nrho,nb_WP0
      integer      :: WP0_nb_CleanChannel
      character (len=Line_len) :: file_WP0
      real (kind=Rkind) :: th_WP0
      integer           :: TFnexp2
      TYPE (REAL_WU)    :: TFmaxE,TFminE

!------ initial WP definition -----------------------------
!     for each variable Qi : exp[-((Q-Qeq)/sigma)2+i*imp_k*(Q-Qeq)+i*phase]
      real (kind=Rkind) :: sigma,imp_k,Qeq,phase

!------ initial WP for the control -----------------------
      integer       :: nb_WP,nb_WPba,max_iter
      real (kind=Rkind) :: conv,alpha,Max_alpha,gamma,Tenvelopp
      logical       :: post_control,gate,cplx_gate
      logical       :: Krotov,Turinici,envelopp,Obj_TO_alpha
      real (kind=Rkind),pointer    :: MRgate(:,:)
      complex (kind=Rkind),pointer :: MCgate(:,:)
      integer       :: err
!------ working variables ---------------------------------
      integer   :: i

      NAMELIST /propa/WPTmax,WPdeltaT,nb_micro,restart,                 &
                      type_WPpropa,spectral,nb_vp_spec,                 &
                      write_iter,n_WPecri,                              &
                  WPpsi2,WPpsi,file_WP,write_DVR,write_FBR,write_WPAdia,&
                      lect_WP0DVR,lect_WP0FBR,lect_WP0FBRall,           &
                      WP0FBR,                                           &
                      WP0n_h,WP0nb_elec,WP0_DIP,th_WP0,WP0nrho,         &
                      WP0_nb_CleanChannel,                              &
                      WP0restart,file_WP0,WP0cplx,                      &
                      nb_WP0,New_Read_WP0,read_listWP0,read_file,       &
                      DHmax,auto_Hmax,                                  &
                      max_poly,poly_tol,npoly,                          &
                      type_corr,Op_corr,i_qa_corr,channel_ie_corr,      &
                      file_autocorr,file_spectrum,file_restart,         &
                      TFnexp2,TFmaxE,TFminE

      NAMELIST /defWP0/sigma,Qeq,imp_k,phase
      NAMELIST /control/nb_WP,nb_WPba,max_iter,conv,alpha,Max_alpha,    &
                        gamma,                                          &
                        Krotov,Turinici,envelopp,Tenvelopp,Obj_TO_alpha,&
                        post_control,gate,cplx_gate

!----- for debuging ----------------------------------------
      character (len=*), parameter :: name_sub='read_propagation'
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

      IF(MPI_id==0) write(out_unitp,*) ' PROPAGATION PARAMETERS: propa, defWP0, control'
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')

!--------- memory allocation -----------------------------
        IF (nb_act1 < 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nb_act1 is < 1 !!',nb_act1
          STOP
        END IF

!------- read the namelist -----------------------------

        WPdeltaT            = REAL_WU(TEN,'au','t')  ! 10 ua (time)
        WPTmax              = REAL_WU(HUNDRED,'au','t')  ! 100 ua (time)


        nb_micro            = 1
        restart             = .FALSE.

        DHmax               = -TEN
        auto_Hmax           = .FALSE.
        max_poly            = 5000
        npoly               = 0
        poly_tol            = ZERO

        type_WPpropa        = 1
        name_WPpropa        = 'cheby'
        spectral            = .FALSE.
        nb_vp_spec          = 0

        write_iter          = .TRUE.
        n_WPecri            = 1
        WPpsi2              = .FALSE.
        WPpsi               = .FALSE.
        write_DVR           = .TRUE.
        write_FBR           = .FALSE.
        write_WPAdia        = .FALSE.
        file_autocorr       = 'file_auto'
        file_spectrum       = 'file_spectrum'
        file_WP             = 'file_WP'
        file_restart        = 'file_WP_restart'

        nb_WP0              = 1
        read_file           = .FALSE.
        New_Read_WP0        = .FALSE.
        read_listWP0        = .FALSE.
        WP0restart          = .FALSE.
        WP0cplx             = .TRUE.
        file_WP0            = 'file_WP0'
        lect_WP0DVR         = .FALSE.
        lect_WP0FBR         = .FALSE.
        lect_WP0FBRall      = .TRUE.
        WP0FBR              = .TRUE.
        WP0n_h              = 1
        WP0nb_elec          = 1
        WP0_DIP             = 0
        th_WP0              = ZERO
        WP0_nb_CleanChannel = 0
        WP0nrho             = -1

        TFnexp2             = 15
        TFminE              = REAL_WU(0,'cm-1','E')  ! 0 cm-1
        TFmaxE              = REAL_WU(TEN**4,'cm-1','E')  ! 10000 cm-1
        type_corr           = 0
        i_qa_corr           = 0
        Op_corr             = 0
        channel_ie_corr     = 0

        read(in_unitp,propa)

        IF (n_WPecri < 0) THEN
          write(out_unitp,*) ' WARNING: n_WPecri should be > 0'
          write(out_unitp,*) '   => n_WPecri=1'
          n_WPecri      = 1
        END IF

        IF (print_level > 0) write(out_unitp,propa)


        para_propa%WPTmax       = convRWU_WorkingUnit_TO_R(WPTmax)
        para_propa%WPdeltaT     = convRWU_WorkingUnit_TO_R(WPdeltaT)

        para_propa%nb_micro     = nb_micro

        para_propa%para_poly%max_poly     = max_poly
        para_propa%para_poly%npoly        = npoly

        CALL alloc_param_poly(para_propa%para_poly)

        CALL alloc_array(para_propa%work_WP,(/npoly+6/),                &
                        "para_propa%work_WP",name_sub,(/0/))

        para_propa%spectral      = spectral
        nb_vp_spec_out           = nb_vp_spec
        para_propa%with_field    = .FALSE.
        IF (type_WPpropa .EQ.1) THEN
!         Chebychev
          IF (poly_tol .EQ. ZERO) poly_tol = ONETENTH**8
          IF (DHmax .EQ. -TEN) DHmax = HALF
          name_WPpropa  = 'cheby'
        ELSE IF (type_WPpropa .EQ.2) THEN
!         Taylor
          IF (poly_tol .EQ. ZERO) poly_tol = ONETENTH**20
          IF (DHmax .EQ. -TEN) DHmax = ZERO
          name_WPpropa  = 'nOD'
        ELSE IF (type_WPpropa .EQ. 3 .OR. type_WPpropa .EQ. -3) THEN
!         im SOD
          IF (poly_tol .EQ. ZERO) poly_tol = ONETENTH**8
          IF (DHmax .EQ. -TEN) DHmax = ZERO
          IF (type_WPpropa .EQ.  3) name_WPpropa  = 'Emin'
          IF (type_WPpropa .EQ. -3) name_WPpropa  = 'Emax'
        ELSE IF (type_WPpropa .EQ. 33 .OR. type_WPpropa .EQ. -33) THEN
!         Davidson
          IF (poly_tol .EQ. ZERO) poly_tol = ONETENTH**5
          IF (DHmax .EQ. -TEN) DHmax = ZERO
          name_WPpropa  = 'Davidson'
        ELSE IF (type_WPpropa .EQ. 34) THEN
!         im SOD
          IF (poly_tol .EQ. ZERO) poly_tol = ONETENTH**5
          IF (DHmax .EQ. -TEN) DHmax = ZERO
          name_WPpropa  = 'Conjugated Gradient'
        ELSE IF (type_WPpropa==22 .OR. type_WPpropa==221 .OR.           &
                 type_WPpropa==222 .OR. type_WPpropa==223) THEN
!         nOD with a time dependant pulse in Hamiltonian (W(t))
          IF (poly_tol .EQ. ZERO) poly_tol = ONETENTH**20
          IF (DHmax .EQ. -TEN) DHmax = ZERO
          name_WPpropa  = 'TDH nOD'
          para_propa%with_field    = .TRUE.
        ELSE IF (type_WPpropa==50 .OR. type_WPpropa==54 .OR.            &
                 type_WPpropa==52 ) THEN
!         RK4 with a time dependant pulse in Hamiltonian (W(t))
          IF (poly_tol .EQ. ZERO) poly_tol = ONETENTH**20
          IF (DHmax .EQ. -TEN) DHmax = ZERO
          name_WPpropa  = 'TDH_RK4'
          IF (type_WPpropa==52) name_WPpropa  = 'TDH_RK2'
          para_propa%with_field    = .TRUE.
        ELSE IF (type_WPpropa==5) THEN
!         RK4 without a time dependant pulse in Hamiltonian (W(t))
          IF (poly_tol .EQ. ZERO) poly_tol = ONETENTH**20
          IF (DHmax .EQ. -TEN) DHmax = ZERO
          name_WPpropa  = 'RK4'
          para_propa%with_field    = .FALSE.

        ELSE IF (type_WPpropa==6) THEN
!         ModMidPoint without a time dependant pulse in Hamiltonian (W(t))
          IF (poly_tol .EQ. ZERO) poly_tol = ONETENTH**20
          IF (DHmax .EQ. -TEN) DHmax = ZERO
          name_WPpropa  = 'ModMidPoint'
          para_propa%with_field    = .FALSE.

        ELSE IF (type_WPpropa==7) THEN
!         Bulirsch-Stoer without a time dependant pulse in Hamiltonian (W(t))
          IF (poly_tol .EQ. ZERO) poly_tol = ONETENTH**20
          IF (DHmax .EQ. -TEN) DHmax = ZERO
          name_WPpropa  = 'Bulirsch-Stoer'
          para_propa%with_field    = .FALSE.

        ELSE IF (type_WPpropa==24 .OR. type_WPpropa==241 .OR.           &
                 type_WPpropa==242 .OR. type_WPpropa==243) THEN
!         nOD with a time dependant pulse in Hamiltonian (W(t))
!         for the control only
          IF (poly_tol .EQ. ZERO) poly_tol = ONETENTH**20
          IF (DHmax .EQ. -TEN) DHmax = ZERO
          name_WPpropa  = 'TDnOD'
          para_propa%with_field    = .TRUE.
        ELSE IF (type_WPpropa == 100) THEN
!         test
          IF (poly_tol .EQ. ZERO) poly_tol = ONETENTH**20
          IF (DHmax .EQ. -TEN) DHmax = ZERO
          name_WPpropa  = 'TDH test'
          para_propa%with_field    = .TRUE.
        ELSE
           write(out_unitp,*) 'ERROR in ',name_sub
           write(out_unitp,*) ' type of propagation(',type_WPpropa,')'
           write(out_unitp,*) ' is not possible'
           STOP
        ENDIF


        para_propa%para_poly%poly_tol     = poly_tol
        para_propa%para_poly%DHmax        = DHmax
        para_propa%auto_Hmax              = auto_Hmax


        para_propa%type_WPpropa = type_WPpropa
        para_propa%name_WPpropa = name_WPpropa

        para_propa%write_iter         = write_iter
        para_propa%n_WPecri           = n_WPecri
        para_propa%WPpsi2             = WPpsi2
        para_propa%WPpsi              = WPpsi
        para_propa%write_GridRep      = write_DVR
        para_propa%write_BasisRep     = write_FBR
        para_propa%write_WPAdia       = write_WPAdia

        para_propa%file_autocorr%name    = make_FileName(file_autocorr)
        IF (err_file_name(para_propa%file_autocorr%name,name_sub='read_propagation') /= 0) THEN
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) '  the file_autocorr file name is empty'
          write(out_unitp,*) '  => check your data!!'
          STOP
        END IF

        para_propa%file_spectrum%name    = make_FileName(file_spectrum)
        IF (err_file_name(para_propa%file_spectrum%name,name_sub='read_propagation') /= 0) THEN
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) '  the file_spectrum file name is empty'
          write(out_unitp,*) '  => check your data!!'
          STOP
        END IF

        para_propa%file_WP%name          = make_FileName(file_WP)
        IF (err_file_name(para_propa%file_WP%name,name_sub='read_propagation') /= 0) THEN
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) '  the file_WP file name is empty'
          write(out_unitp,*) '  => check your data!!'
          STOP
        END IF

        para_propa%file_WP_restart%name  = make_FileName(file_restart)
        IF (restart .AND. .NOT. check_file_exist_WITH_file_name(para_propa%file_WP_restart%name)) THEN
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) '  the restart file does not exist or its file name is empty'
          write(out_unitp,*) '  file_restart: ',file_restart
          write(out_unitp,*) '  => check your data!!'
          STOP
        END IF

        para_propa%para_WP0%nb_WP0        = nb_WP0
        para_propa%para_WP0%New_Read_WP0  = New_Read_WP0
        para_propa%para_WP0%read_file     = read_file
        para_propa%para_WP0%read_listWP0  = read_listWP0
        para_propa%para_WP0%lect_WP0GridRep   = lect_WP0DVR
        para_propa%para_WP0%lect_WP0BasisRep   = lect_WP0FBR
        para_propa%para_WP0%lect_WP0BasisRepall= lect_WP0FBRall
        para_propa%para_WP0%WP0BasisRep        = WP0FBR
        para_propa%para_WP0%WP0n_h        = WP0n_h
        para_propa%para_WP0%WP0nb_elec    = WP0nb_elec
        para_propa%para_WP0%WP0_DIP       = WP0_DIP
        para_propa%para_WP0%th_WP0        = th_WP0
        para_propa%para_WP0%WP0nrho       = WP0nrho
        para_propa%para_WP0%WP0restart    = WP0restart
        para_propa%para_WP0%WP0cplx       = WP0cplx
        para_propa%para_WP0%file_WP0%name = make_FileName(file_WP0)
        IF (read_file .AND. .NOT.                                       &
            check_file_exist_WITH_file_name(para_propa%para_WP0%file_WP0%name)) THEN
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) '  the file_WP0 file does not exist or its file name is empty'
          write(out_unitp,*) '  file_WP0: ',file_WP0
          write(out_unitp,*) '  => check your data!!'
          STOP
        END IF
        para_propa%para_WP0%WP0_nb_CleanChannel = WP0_nb_CleanChannel
        IF (WP0_nb_CleanChannel > 0) THEN
          write(out_unitp,*)
          write(out_unitp,*) "===================================="
          write(out_unitp,*) "== WP0_CleanChannel ================"
          CALL alloc_param_WP0(para_propa%para_WP0,                     &
                        WP0Grid_Gaussian=.FALSE.,WP0_CleanChannel=.TRUE.)
          read(in_unitp,*) para_propa%para_WP0%WP0_CleanChannellist
          write(out_unitp,*) " list: ",para_propa%para_WP0%WP0_CleanChannellist
          write(out_unitp,*)
          write(out_unitp,*) "===================================="
        END IF

        IF (type_corr < 0 .AND. type_corr > 2) THEN
          write(out_unitp,*) ' ERROR in read_propa'
          write(out_unitp,*) ' type_corr < 0 and type_corr >2'
          STOP
        END IF
        IF (type_corr == 1 .AND. (i_qa_corr <1 .OR. channel_ie_corr<1)) THEN
          write(out_unitp,*) ' ERROR in read_propa'
          write(out_unitp,*) ' type_corr == 1 and ',                             &
                      'wrong i_qa_corr or channel_ie_corr'
          write(out_unitp,*) 'type_corr,i_qa_corr,channel_ie_corr',             &
                      type_corr,i_qa_corr,channel_ie_corr
          STOP
        END IF
        para_propa%type_corr       = type_corr
        para_propa%Op_corr         = Op_corr
        para_propa%i_qa_corr       = i_qa_corr
        para_propa%channel_ie_corr = channel_ie_corr
        para_propa%restart         = restart


        para_propa%TFnexp2        = TFnexp2


        para_propa%TFmaxE       = convRWU_WorkingUnit_TO_R(TFmaxE)
        para_propa%TFminE       = convRWU_WorkingUnit_TO_R(TFminE)

        IF (para_propa%control) THEN
          para_propa%with_field    = .TRUE.
          nb_WP        = 1
          nb_WPba      = 0
          max_iter     = 0
          conv         = ONE - ONETENTH**3
          alpha        = ONE
          Max_alpha    = TEN**3
          gamma        = ONETENTH
          Krotov       = .TRUE.
          Turinici     = .FALSE.
          envelopp     = .TRUE.
          Tenvelopp    = ZERO
          Obj_TO_alpha = .FALSE.
          post_control = .TRUE.
          gate         = .FALSE.
          cplx_gate    = .FALSE.
          read(in_unitp,control)

          IF (Tenvelopp == ZERO) Tenvelopp = convRWU_TO_R(WPTmax)
          IF (print_level > 0) write(out_unitp,control)

          IF (nb_WP < 1) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' nb_WP < 1 !!!',nb_WP
            STOP
          END IF
          IF (nb_WPba < 1) nb_WPba = nb_WP

          para_propa%para_control%nb_WP           = nb_WP
          para_propa%para_control%nb_WPba         = nb_WPba
          para_propa%para_control%max_iter        = max_iter
          para_propa%para_control%conv            = conv
          para_propa%para_control%alpha           = alpha
          para_propa%para_control%Max_alpha       = Max_alpha
          para_propa%para_control%gamma           = gamma

          para_propa%para_control%Krotov          = Krotov
          para_propa%para_control%Turinici        = Turinici
          para_propa%para_control%envelopp        = envelopp
          para_propa%para_control%Tenvelopp       = Tenvelopp
          para_propa%para_control%Obj_TO_alpha    = Obj_TO_alpha

          para_propa%para_control%post_control    = post_control
          para_propa%para_control%gate            = gate
          para_propa%para_control%cplx_gate       = cplx_gate
          IF (gate) THEN
            CALL alloc_array(para_propa%para_control%tab_WP0,           &
                                                         (/nb_WPba/),   &
                            "para_propa%para_control%tab_WP0",name_sub)
            CALL alloc_array(para_propa%para_control%Mgate0,            &
                                                   (/nb_WP,nb_WPba/),   &
                            "para_propa%para_control%Mgate0",name_sub)
            CALL alloc_array(para_propa%para_control%Mgatet,            &
                                                   (/nb_WP,nb_WPba/),   &
                            "para_propa%para_control%Mgatet",name_sub)
            read(in_unitp,*) para_propa%para_control%tab_WP0
            write(out_unitp,*) 'tab_WP0',para_propa%para_control%tab_WP0

            IF (cplx_gate) THEN
              nullify(MCgate)
              CALL alloc_array(MCgate,(/nb_WP,nb_WPba/),"MCgate",name_sub)

              CALL Read_Mat(MCgate,in_unitp,nb_WPba,err)
              IF (err /= 0) THEN
                write(out_unitp,*) 'ERROR in ',name_sub
                write(out_unitp,*) ' reading the matrix "Mgate0"'
                STOP
              END IF
              para_propa%para_control%Mgate0(:,:) = MCgate(:,:)

              CALL Read_Mat(MCgate,in_unitp,nb_WPba,err)
              IF (err /= 0) THEN
                write(out_unitp,*) 'ERROR in ',name_sub
                write(out_unitp,*) ' reading the matrix "Mgatet"'
                STOP
              END IF
              para_propa%para_control%Mgatet(:,:) = MCgate(:,:)

              CALL dealloc_array(MCgate,"MCgate",name_sub)
            ELSE
              nullify(MRgate)
              CALL alloc_array(MRgate,(/nb_WP,nb_WPba/),"MRgate",name_sub)

              CALL Read_Mat(MRgate,in_unitp,nb_WPba,err)
              IF (err /= 0) THEN
                write(out_unitp,*) 'ERROR in ',name_sub
                write(out_unitp,*) ' reading the matrix "Mgate0"'
                STOP
              END IF
              para_propa%para_control%Mgate0(:,:) = MRgate(:,:)

              CALL Read_Mat(MRgate,in_unitp,nb_WPba,err)
              IF (err /= 0) THEN
                write(out_unitp,*) 'ERROR in ',name_sub
                write(out_unitp,*) ' reading the matrix "Mgatet"'
                STOP
              END IF
              para_propa%para_control%Mgatet(:,:) = MRgate(:,:)

              CALL dealloc_array(MRgate,"MRgate",name_sub)
            END IF

            IF (debug) THEN
              write(out_unitp,*) 'Gate Matrix (WP0)',nb_WP
              CALL Write_Mat(para_propa%para_control%Mgate0,out_unitp,4)
              write(out_unitp,*) 'Gate Matrix (WPt)',nb_WP
              CALL Write_Mat(para_propa%para_control%Mgatet,out_unitp,4)
            END IF

            DO i=1,nb_WP
              CALL Write_Vec(para_propa%para_control%Mgate0(i,:),       &
                             out_unitp,4,name_info="#WP0 " // int_TO_char(i) )
              CALL Write_Vec(para_propa%para_control%Mgatet(i,:),       &
                             out_unitp,4,name_info="#WPt " // int_TO_char(i))
            END DO
          ELSE
            CALL alloc_array(para_propa%para_control%tab_WP0,(/nb_WP/), &
                            "para_propa%para_control%tab_WP0",name_sub)
            CALL alloc_array(para_propa%para_control%tab_WPt,(/nb_WP/), &
                            "para_propa%para_control%tab_WPt",name_sub)

            read(in_unitp,*) para_propa%para_control%tab_WP0
            read(in_unitp,*) para_propa%para_control%tab_WPt
            write(out_unitp,*) 'tab_WP0',para_propa%para_control%tab_WP0
            write(out_unitp,*) 'tab_WPt',para_propa%para_control%tab_WPt
          END IF

        ELSE
          IF (lect_WP0DVR) THEN
            write(out_unitp,*) ' read the WP0 in GridRep'
          ELSE IF (lect_WP0FBR) THEN
            write(out_unitp,*) ' read the WP0 in BasisRep (old way)'
          ELSE IF (New_Read_WP0) THEN
            write(out_unitp,*) ' read the WP0 in BasisRep (new way)'
          ELSE

            para_propa%para_WP0%nb_act1 = nb_act1

            CALL alloc_param_WP0(para_propa%para_WP0,                   &
                       WP0Grid_Gaussian=.TRUE.,WP0_CleanChannel=.FALSE.)

            write(out_unitp,*)
            IF (print_level > 0)                                        &
             write(out_unitp,*) 'WP0(Q)=exp[-((Q-Qeq)/sigma)^2+i*imp_k*(Q-Qeq)+i*phase]'
            IF (print_level > 0) write(out_unitp,*) 'WP0sigma WP0Qeq WP0imp_k WP0phase'
            DO i=1,nb_act1
              sigma    = ONETENTH
              Qeq      = ZERO
              imp_k    = ZERO
              phase    = ZERO
              read(in_unitp,defWP0)
              para_propa%para_WP0%WP0sigma(i) = sigma
              para_propa%para_WP0%WP0Qeq(i)   = Qeq
              para_propa%para_WP0%WP0imp_k(i) = imp_k
              para_propa%para_WP0%WP0phase(i) = phase

              IF (print_level > 0) write(out_unitp,*) i,sigma,Qeq,imp_k,phase
            END DO

          END IF
        END IF

      END SUBROUTINE read_propagation
!=======================================================
!
!     read parameters for the Davidson
!
!=======================================================
      SUBROUTINE read_davidson(para_Davidson,para_propa)
      USE mod_system
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_Davidson) :: para_Davidson
      TYPE (param_propa)    :: para_propa

!----- for the namelist -------------------------------------
      integer :: num_resetH,num_checkS
      integer :: max_it,nb_WP,max_WP,num_LowestWP

      logical :: read_WP,read_listWP,precond,formatted_file_readWP
      integer :: nb_readWP,nb_readWP_OF_List
      real (kind=Rkind) :: precond_tol
      character (len=Line_len) :: name_file_readWP

      logical :: save_all
       real (kind=Rkind) :: scaled_max_ene
      integer :: save_max_nb
      character (len=Line_len) :: name_file_saveWP

      logical :: lower_states,all_lower_states,project_wp0
      logical :: one_residue,one_by_one,With_Grid
      TYPE (REAL_WU)    :: Max_ene

      integer           :: conv_hermitian
      TYPE (REAL_WU)    :: conv_ene
      TYPE (REAL_WU)    :: conv_resi

      integer           :: symab
      integer           :: NewVec_type,residual_max_nb

      TYPE(REAL_WU)     :: E0_filter ! the center of the window
      TYPE(REAL_WU)     :: W_filter  ! the width of the window

      integer           :: L_filter      = 50   ! the size of filter basis
      integer           :: M_filter      = 100  ! the number of terms in the Chebyshev expansion
      integer           :: DeltaM_filter = 0    ! the number of terms in the Chebyshev expansion
      integer           :: Mmax_filter   = 0
      logical           :: auto_Hmax,Hmin_propa,Hmax_propa
      real (kind=Rkind) :: DHmax,Hmin,Hmax
      integer           :: max_poly
      real (kind=Rkind) :: poly_tol
      logical           :: Op_Transfo
      integer           :: err_io

      namelist /davidson/num_resetH,num_checkS,                         &
                        read_WP,nb_readWP,read_listWP,nb_readWP_OF_List,&
                         precond,precond_tol,NewVec_type,               &
                         name_file_readWP,formatted_file_readWP,        &
                         lower_states,one_residue,project_wp0,          &
                         residual_max_nb,one_by_one,                    &
                         all_lower_states,Max_ene,num_LowestWP,         &
                         conv_ene,conv_resi,symab,                      &
                         max_it,nb_WP,max_WP,conv_hermitian,With_Grid,  &
                   scaled_max_ene,save_all,save_max_nb,name_file_saveWP,&
                         Op_Transfo,E0_filter,W_filter,                 &
                         L_filter,M_filter,DeltaM_filter,Mmax_filter,   &
                       DHmax,auto_Hmax,Hmin,Hmax,Hmin_propa,Hmax_propa, &
                         max_poly,poly_tol

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='read_Davidson'
      logical, parameter :: debug =.FALSE.
!      logical, parameter :: debug =.TRUE.
!-----------------------------------------------------------
      write(out_unitp,*) ' DAVIDSON PARAMETERS'
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
      END IF
!-----------------------------------------------------------
      num_resetH            = -1 ! no reset for Davidson ortherwise restart davidson with nb_diago WP
      num_checkS            = 10 ! check overlap matrix

      nb_WP                 = 0
      max_WP                = 0
      read_WP               = .FALSE. ! T => read the the WP in file_WP
      read_listWP           = .FALSE. ! T => read a list of WP (nb_readWP)
      nb_readWP_OF_List     = 0    ! 0 => all
      nb_readWP             = -1      ! if -1 => nb_diago, ortherwise more WP
      num_LowestWP          = -1        ! With this number, the WP with lower energy are excluded
      name_file_readWP      = 'file_WP'
      formatted_file_readWP = para_Davidson%formatted_file_readWP ! to restor the default defined in ini_data

      precond           = .FALSE. ! precondition
      precond_tol       = ONETENTH

      lower_states      = .FALSE. ! T => converge on the lower states
      all_lower_states  = .FALSE. ! T => converge all the lower states (E<Max_ene)
      project_WP0       = .FALSE. ! T => project on WP0 (revelant if lower_states=f)
      one_residue       = .FALSE. ! T => increase H with ONE residue (largest)
      one_by_one        = .FALSE. ! T => converge vector one by one
      residual_max_nb   = huge(1) ! the number maximal of residual

      Max_ene           = REAL_WU(FIVE*TEN**3,'cm-1','E')  ! 5000 cm-1

      symab             = -1
      conv_ene          = REAL_WU(ONETENTH**7,'au','E')
      conv_resi         = REAL_WU(FIVE*ONETENTH**5,'au','E')
      max_it            = 100
      conv_hermitian    = 0 ! convergence if max_ene < 10**conv_hermitian * "non_hermitic"
      NewVec_type       = 2
      With_Grid         = .FALSE.

      scaled_max_ene    = 1.1_Rkind ! save_max_ene = scaled_max_ene * max_ene
      save_all          = .FALSE.
      save_max_nb       = -1 ! If -1, we did not use this number (default)
      name_file_saveWP  = 'file_WP'

      ! for the filter diagonalization
      E0_filter         = REAL_WU(ZERO,      'cm-1','E')  ! the center of the window in ua (but read in cm-1)
      W_filter          = REAL_WU(200._Rkind,'cm-1','E')  ! the window in ua (but read in cm-1)

      L_filter          = 100   ! the size of filter basis set
      M_filter          = 300  ! it is equivalent to npoly in the propgation subroutines
      DeltaM_filter     = 0
      Mmax_filter       = 30000

      DHmax             = -TEN
      Hmin              =  huge(ONE)
      Hmax              = -huge(ONE)
      auto_Hmax         = .FALSE. ! to calculate automatically Hmin end Hmax
      Hmin_propa        = .FALSE. ! with Davidson
      Hmax_propa        = .FALSE. ! with Davidson

      max_poly          = 500
      poly_tol          = ZERO
      Op_Transfo        = .FALSE.

      read(in_unitp,davidson,IOSTAT=err_io)
      IF (print_level > 0) write(out_unitp,davidson)

      IF (.NOT. lower_states .AND. .NOT. all_lower_states .AND. &
          .NOT. project_WP0) lower_states = .TRUE.


      IF (residual_max_nb == huge(1) .AND. one_residue) residual_max_nb = 1

      ! for the filter diagonalization
      IF (DeltaM_filter == 0) DeltaM_filter = L_filter
      para_Davidson%E0_filter     = convRWU_TO_R(E0_filter)
      para_Davidson%W_filter      = convRWU_TO_R(W_filter)

      para_Davidson%L_filter      = L_filter
      para_Davidson%M_filter      = M_filter
      para_Davidson%DeltaM_filter = DeltaM_filter
      para_Davidson%Mmax_filter   = Mmax_filter

      IF (poly_tol == ZERO) poly_tol = ONETENTH**8
      IF (DHmax == -TEN)    DHmax = HALF
      para_propa%para_poly%poly_tol     = poly_tol
      para_propa%para_poly%max_poly     = max_poly
      para_propa%para_poly%npoly        = M_filter
      para_propa%para_poly%DHmax        = DHmax
      para_propa%auto_Hmax              = auto_Hmax
      para_propa%Hmax                   = Hmax
      para_propa%Hmin                   = Hmin

      para_Davidson%Hmin_propa          = Hmin_propa
      para_Davidson%Hmax_propa          = Hmax_propa
      IF (Hmin_propa .OR. Hmax_propa) THEN
        nb_WP            = 0
        lower_states     = .TRUE.
        project_WP0      = .FALSE.
        all_lower_states = .FALSE.
        max_it           = 50
        num_resetH       = max_it
        num_checkS       = max_it
        NewVec_type      = 2
      END IF

      para_Davidson%Op_Transfo          = Op_Transfo


      IF (project_WP0 .EQV. (lower_states .OR. all_lower_states)) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' You have to chose between ',                       &
             '"lower_states" or "all_lower_states" or "project_WP0"',   &
                            lower_states,all_lower_states,project_WP0
        write(out_unitp,*) ' CHECK your data'
        STOP
      END IF

      IF (all_lower_states) THEN
        lower_states = .TRUE.
      END IF

      para_Davidson%num_resetH            = num_resetH
      para_Davidson%num_checkS            = num_checkS

      para_Davidson%num_LowestWP          = num_LowestWP

      para_Davidson%read_WP               = read_WP
      para_Davidson%nb_readWP             = nb_readWP
      para_Davidson%nb_readWP_OF_List     = nb_readWP_OF_List
      para_Davidson%read_listWP           = read_listWP
      para_Davidson%precond               = precond
      para_Davidson%precond_tol           = precond_tol
      para_Davidson%name_file_readWP      = make_FileName(name_file_readWP)
      para_Davidson%formatted_file_readWP = formatted_file_readWP
      IF (read_WP .AND. .NOT. check_file_exist_WITH_file_name(para_Davidson%name_file_readWP)) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) '  the name_file_readWP does not exist or its file name is empty'
        write(out_unitp,*) '  name_file_readWP: ',name_file_readWP
        write(out_unitp,*) '  => check your data!!'
        STOP
      END IF

      para_Davidson%lower_states          = lower_states
      para_Davidson%all_lower_states      = all_lower_states
      para_Davidson%Max_ene               = convRWU_TO_R(Max_ene)

      para_Davidson%project_WP0           = project_WP0
      para_Davidson%residual_max_nb       = residual_max_nb
      para_Davidson%one_by_one            = one_by_one
      para_Davidson%NewVec_type           = NewVec_type
      para_Davidson%With_Grid             = With_Grid


      para_Davidson%symab             = symab
      para_Davidson%conv_ene          = convRWU_TO_R(conv_ene)
      para_Davidson%conv_resi         = convRWU_TO_R(conv_resi)

      para_Davidson%max_it            = max_it

      para_Davidson%nb_WP             = nb_WP
      para_Davidson%max_WP            = max_WP
      para_Davidson%conv_hermitian    = conv_hermitian

      para_Davidson%save_all          = save_all
      para_Davidson%scaled_max_ene    = scaled_max_ene
      para_Davidson%save_max_nb       = save_max_nb
      para_Davidson%name_file_saveWP  = make_FileName(name_file_saveWP)
      IF (err_file_name(para_Davidson%name_file_saveWP,name_sub='read_davidson') /= 0) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) '  the name_file_saveWP file name is empty'
        write(out_unitp,*) '  => check your data!!'
        STOP
      END IF

      IF(MPI_id==0) THEN
        write(out_unitp,*) 'para_Davidson%Max_ene       ',para_Davidson%Max_ene
        write(out_unitp,*) 'para_Davidson%scaled_max_ene',para_Davidson%scaled_max_ene
        write(out_unitp,*) 'para_Davidson%conv_ene      ',para_Davidson%conv_ene
        write(out_unitp,*) 'para_Davidson%conv_resi     ',para_Davidson%conv_resi
      ENDIF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE read_davidson

      END MODULE mod_propa

