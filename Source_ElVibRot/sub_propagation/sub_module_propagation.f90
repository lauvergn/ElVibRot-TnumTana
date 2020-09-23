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
      MODULE mod_propa
      USE mod_system
      USE mod_Constant,  ONLY : get_conv_au_to_unit,real_wu,            &
                                convRWU_TO_R_WITH_WorkingUnit,          &
                                convRWU_TO_R_WITH_WritingUnit,RWU_WriteUnit
      USE mod_psi,       ONLY : param_WP0,param_ana_psi
      USE mod_field,     ONLY : param_field
      IMPLICIT NONE

PRIVATE
PUBLIC :: param_poly,param_control,param_Davidson,param_propa
PUBLIC :: read_propagation,read_davidson
PUBLIC :: dealloc_param_propa,sub_analyze_WP_OpWP,sub_analyze_mini_WP_OpWP
PUBLIC :: Read_AutoCorr,Write_AutoCorr,Calc_AutoCorr
PUBLIC :: SaveWP_restart,ReadWP_restart

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
        TYPE param_poly

        logical                    :: init_done       = .FALSE.
        real (kind=Rkind)          :: DHmax           = -TEN
        real (kind=Rkind)          :: Hmax            = ZERO
        real (kind=Rkind)          :: Hmin            = ZERO
        real (kind=Rkind)          :: deltaE          = ZERO
        real (kind=Rkind)          :: Esc             = ONE
        real (kind=Rkind)          :: E0              = ZERO
        real (kind=Rkind)          :: alpha           = ZERO
        real (kind=Rkind)          :: poly_tol        = ZERO
        integer                    :: max_poly        = 500 ! max number of polynomials
        real (kind=Rkind), pointer :: coef_poly(:)    => null()  !
        integer                    :: npoly           = 0 ! number polynomials
        integer                    :: npoly_Opt       = 0 ! optimal number polynomials

        END TYPE param_poly
        TYPE param_control

        logical          :: restart            ! restart

        integer          :: nb_WP              ! number of initial WP and WP target's
        integer          :: nb_WPba            ! number of basis function of each WP
        integer          :: max_iter           ! number of iterations for the control
        real (kind=Rkind):: conv               ! convergence for the objective
        real (kind=Rkind):: alpha,gamma        ! parameters for the convergence
        real (kind=Rkind):: Max_alpha          ! parameters for the convergence
        logical          :: Krotov             ! Krotov Algorithm
        logical          :: Turinici           ! Turinici Algorithm
        logical          :: envelopp           ! use an envelopp (s(t)) for the field
        real (kind=Rkind):: Tenvelopp          ! parameter of the envelopp sin(t/T*pi)
        logical          :: Obj_TO_alpha       ! the objectifs modify alpha

        integer, pointer :: tab_WP0(:) => null()        ! the intial WPs (nb_WP)
        integer, pointer :: tab_WPt(:) => null()        ! the target WPs (nb_WP)

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

         integer :: num_resetH =-1 ! -1 no reset
         integer :: num_checkS = 10
         integer :: max_it     = 100

         integer :: nb_WP        = 0
         integer :: max_WP       = 0
         integer :: num_LowestWP = -1        ! With this number, the WP with lower energy are excluded
         integer :: nb_WP0       = 0


         logical :: read_WP                            = .FALSE. ! T => read the the WP in file_WP
         logical :: read_listWP                        = .FALSE. ! T => read a list of WP (nb_readWP)
         logical :: formatted_file_readWP              = .TRUE.
         integer :: nb_readWP                          = -1      ! if -1 => nb_diago, ortherwise more WP
         integer :: nb_readWP_OF_List                  = 0     ! 0 => all
         character (len=Line_len) :: name_file_readWP  = 'file_WP'

         logical :: precond                            = .FALSE. ! precondition (not used)
         real (kind=Rkind) :: precond_tol              = ONETENTH

         logical :: save_all                           = .FALSE.
         integer :: save_interal                       = 1 ! for MPI ???
         real (kind=Rkind) :: save_max_ene             = -huge(ONE)
         real (kind=Rkind) :: scaled_max_ene           = 1.1_Rkind ! save_max_ene = scaled_max_ene * max_ene
         integer :: save_max_nb                        = -1 ! If -1, we did not use this number (default)
         character (len=Line_len) :: name_file_saveWP  = 'file_WP'
         logical :: formatted_file_WP                  = .TRUE.

         logical :: all_lower_states  = .FALSE. ! T => converge all the lower states (E<Max_ene)
         logical :: lower_states      = .FALSE. ! T => converge on the lower states
         logical :: project_WP0       = .FALSE. ! T => project on WP0 (revelant if lower_states=f)
         logical :: one_by_one        = .FALSE. ! T => converge vector one by one
         integer :: residual_max_nb   = huge(1) ! the number maximal of residual

         logical :: Op_Transfo        = .FALSE.


         logical :: Hmin_propa = .FALSE.  ! with Davidson
         logical :: Hmax_propa = .FALSE.  ! with Davidson
         real (kind=Rkind) :: thresh_project = 0.8_Rkind

         real (kind=Rkind) :: Max_ene = -HUGE(ONE)
         integer           :: symab = -1 ! -1 nosym


         real (kind=Rkind) :: RMS_ene   = ZERO ! not used
         real (kind=Rkind) :: conv_ene  = ONETENTH**7
         real (kind=Rkind) :: RMS_resi  = ZERO ! not used
         real (kind=Rkind) :: conv_resi = FIVE*ONETENTH**5

         integer :: conv_hermitian    = 0 ! convergence if max_ene < 10**conv_hermitian * "non_hermitic"
         integer :: NewVec_type       = 2
         logical :: With_Grid         = .FALSE.

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

        real (kind=Rkind)   ::  Hmax              = ZERO
        real (kind=Rkind)   ::  Hmin              = ZERO
        logical             ::  once_Hmin         =.TRUE.!< control the calculation of 
                                                         !< Hmin once at the first action
        logical             ::  auto_Hmax         = .FALSE.    !  .TRUE. => Hmax is obtained with a propagation
                                              !            with imaginary time (type_WPpropa=-3)
                                              !            (default .FALSE.)
        real (kind=Rkind)   ::  WPTmax       = HUNDRED ! propagation time in au
        real (kind=Rkind)   ::  WPT0          ! initial time in au
        logical             ::  restart      = .FALSE. ! restart the propagation

        real (kind=Rkind)   ::  WPdeltaT      = TEN ! time step in au
        integer             ::  nb_micro      = 1 ! nb of microiterations
        logical             ::  One_Iteration = .FALSE. ! Propagation stops after one time step.



        integer             ::  max_ana      = 0 ! nb of wavefunctions calculated with imaginary propagation

        integer             ::  type_WPpropa = 0 ! propagation type :
                                              !  1 Chebychev (default)
                                              !  2 nOD (Taylor)
                                              !  5 RK4
                                              !  6 ModmidPoint
                                              !  7 Bulirsch-Stoer
                                              !  8 SIL
                                              !  9 SIP
                                              ! 10 Spectral

                                              !  3 (Emin) or  -3 (Emax) imaginary propagation
                                              ! 33 (Emin) or -33 (Emax) Davidson
                                              ! 34 (Emin) Conjugated Gradient

                                              ! 22 or 24 Time dependant Hamiltonian nOD
                                              ! 52 Time dependant Hamiltonian RK2
                                              ! 54 Time dependant Hamiltonian RK4


        logical                  :: With_field = .FALSE.   ! propagation with a field
        character (len=Name_len) :: name_WPpropa = ' ' !  spectral (not working)
                                              !  cheby (default)
                                              !  nOD
                                              !  Hmin,Hmax (imaginary propagation: nOD)
                                              !  TDnOD

        logical         ::       spectral    = .FALSE. ! Def (f). If T => spectral propagatio


        integer         ::       type_corr       = 0 !  kind of correlation function
                                                     !  default 0 => <psi0(T)|psi(T)>
                                                     !  1 =>  psi(T,i_qa_corr,...) (for a grid point)
        integer         ::       i_qa_corr       = 0 !  grid point for type_corr = 1
        integer         ::       channel_ie_corr = 0 !   Channel for type_corr = 1
        integer         ::       Op_corr         = 0 !  defined the Operator when type_corr=2

        integer         ::       n_WPecri        = 1 !  write WP every n_WPecri time steps
        logical         ::       WPpsi           = .FALSE. !  write WP
        logical         ::       WPpsi2          = .FALSE. !  write WP2
        logical         ::       write_iter      = .TRUE. !  write iteration
        logical         ::       write_GridRep   = .TRUE. !  write WP of the grid
        logical         ::       write_BasisRep  = .FALSE.!  write WP of the basis
        logical         ::       write_WPAdia    = .FALSE. !  if T, write WP on the adiabatic PES (otherwise on the diabatic ones)

        TYPE (param_ana_psi)  :: ana_psi         ! to control the WP analysis


        logical           ::       control       = .FALSE.  !  use control (need type_WPpropa 24)
        logical           ::       test_max_norm = .FALSE.  ! IF .TRUE., Error in the propagation (norm too large)
        real (kind=Rkind) ::       max_norm2 = 1.3_Rkind ! ! if wp%norm^2 > max_norm^2 => stop
        logical           ::       march_error = .FALSE.    ! IF .TRUE., Error in the propagation
        integer           ::       num_Op = 0      !  Operator used for the propagation (H)

!----- variables for the WP propagation ----------------------------

        TYPE (param_file)    :: file_WP,file_autocorr,file_spectrum
        TYPE (param_file)    :: file_WP_restart ! files to be able to write the propagation

        TYPE(param_WP0)      :: para_WP0
        TYPE(param_control)  :: para_control
        TYPE(param_Davidson) :: para_Davidson
        TYPE(param_poly)     :: para_poly
        TYPE(param_field)    :: para_field

        integer           :: TFnexp2 = 15
        real (kind=Rkind) :: TFmaxE  = HUGE(ONE)
        real (kind=Rkind) :: TFminE  = ZERO

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
      USE mod_field, ONLY : dealloc_param_field
      USE mod_psi,   ONLY : param_WP0,dealloc_param_WP0,dealloc_ana_psi,&
                            dealloc_psi,dealloc_array
      IMPLICIT NONE

      TYPE (param_propa), intent(inout) :: para_propa

      integer :: i

      CALL dealloc_ana_psi(para_propa%ana_psi)

      CALL dealloc_param_WP0(para_propa%para_WP0)
      CALL dealloc_param_control(para_propa%para_control)
      CALL dealloc_param_poly(para_propa%para_poly)
      CALL dealloc_param_field(para_propa%para_field)

      END SUBROUTINE dealloc_param_propa

      SUBROUTINE SaveWP_restart(T,WP,file_restart)
      USE mod_system
      USE mod_psi,    ONLY : param_psi
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
        write(out_unitp,*) ' T=',T
        write(out_unitp,*) ' nb_WP,size WP',size(WP),size(WP(1)%CvecB)
        write(out_unitp,*)
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------


      CALL file_open(file_restart,no_restart)

      write(no_restart,*) T,size(WP),size(WP(1)%CvecB)
      IF(keep_MPI) THEN
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
      USE mod_psi,    ONLY : param_psi
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_file), intent(inout) :: file_restart
      TYPE (param_psi),  intent(inout) :: WP(:)
      real (kind=Rkind), intent(inout) :: T

!------ working parameters --------------------------------
      integer       :: i,no_restart,err_read
      integer       :: nb_WP_file,size_WP_file

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
      read(no_restart,*,IOSTAT=err_read) T,nb_WP_file,size_WP_file
      IF (err_read /= 0) THEN
        write(out_unitp,*) ' WARNING in ',name_sub
        write(out_unitp,*) ' T (time) is not present in the restart file'
        write(out_unitp,*) ' file name:',trim(file_restart%name)
        write(out_unitp,*) ' => No restart, T0=0'
        T = ZERO
      ELSE
        write(out_unitp,*) 'T0 for the restart:',T
        write(out_unitp,*) 'ReadWP_restart: nb_WP,size WP            ',size(WP),size(WP(1)%CvecB)
        write(out_unitp,*) 'ReadWP_restart: nb_WP,size WP (from file)',nb_WP_file,size_WP_file

        IF (nb_WP_file /= size(WP) .OR. size_WP_file /= size(WP(1)%CvecB)) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' Inconsistent WP size values!'
          STOP ' ERROR in ReadWP_restart:  Inconsistent WP size values.'
        END IF

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
      USE mod_psi,    ONLY : param_psi,ecri_psi,Overlap_psi1_psi2,      &
                             sub_PsiBasisRep_TO_GridRep
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
        IF(MPI_id==0) CALL Write_AutoCorr(para_propa%file_autocorr%unit,T,cdot)
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


!=======================================================================================
SUBROUTINE sub_analyze_WP_OpWP(T,WP,nb_WP,para_H,para_propa,adia,para_field)
  USE mod_system
  USE mod_Op,              ONLY : param_Op,sub_PsiOpPsi,sub_psiHitermPsi, &
                                  sub_PsiDia_TO_PsiAdia_WITH_MemGrid
  USE mod_field,           ONLY : param_field,sub_dnE
  USE mod_ExactFact

  USE mod_psi, ONLY : param_psi,ecri_psi,alloc_psi,dealloc_psi,         &
                      sub_analyze_psi,norm2_psi,alloc_psi,modif_ana_psi,&
                      sub_PsiBasisRep_TO_GridRep,Write_ana_psi
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

  logical            :: ana_mini = .TRUE.  ! turn off further analysis
!  logical            :: ana_mini = .FALSE.
  logical            :: G,G2,B,B2,With_field

!-- working parameters --------------------------------
  TYPE (param_psi)   :: w1,w2

  integer                            :: j,i,i_bi,i_be,i_bie
  complex (kind=Rkind)               :: ET  ! energy
  character (len=:),    allocatable  :: info

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
   write(out_unitp,*)
   write(out_unitp,*) 'nb_ba,nb_qa',WP(1)%nb_ba,WP(1)%nb_qa
   write(out_unitp,*) 'nb_bi',WP(1)%nb_bi
   write(out_unitp,*)
   CALL Write_ana_psi(para_propa%ana_psi)

   DO i=1,nb_WP
     CALL norm2_psi(WP(i),.FALSE.,.TRUE.,.FALSE.)
     write(out_unitp,*) 'norm2psi0 BasisRep',i,WP(1)%norm2

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
        CALL sub_analyze_mini_WP_OpWP(T,WP,nb_WP,para_H,adia=adia,para_field=para_field)
      ELSE
        CALL sub_analyze_mini_WP_OpWP(T,WP,nb_WP,para_H,adia=adia)
      END IF
    ELSE
      IF (present(para_field)) THEN
        CALL sub_analyze_mini_WP_OpWP(T,WP,nb_WP,para_H,para_field=para_field)
      ELSE
        CALL sub_analyze_mini_WP_OpWP(T,WP,nb_WP,para_H)
      END IF
    END IF
    RETURN
  END IF

  Write_psi2_Grid = para_propa%ana_psi%Write_psi2_Grid
  Write_psi_Grid  = para_propa%ana_psi%Write_psi_Grid

  IF(keep_MPI) THEN
    BasisRep = WP(1)%BasisRep
    GridRep  = WP(1)%GridRep
  ENDIF


  IF (present(adia)) THEN
    adia_loc = adia
  ELSE
    adia_loc = para_propa%Write_WPAdia .OR. para_propa%ana_psi%adia
  END IF
  IF (.NOT. para_H%para_ReadOp%para_FileGrid%Save_MemGrid_done) adia_loc = .FALSE.

  !-----------------------------------------------------------
  ! => the WPs on the Grid
!  IF (.NOT. para_propa%ana_psi%GridDone) THEN
!    IF(openmpi) CALL time_perso('sub_PsiBasisRep_TO_GridRep ini')
!    DO i=1,nb_WP
!      IF(keep_MPI) CALL sub_PsiBasisRep_TO_GridRep(WP(i))
!    END DO
!    IF(openmpi) CALL time_perso('sub_PsiBasisRep_TO_GridRep end')
!  END IF
!  para_propa%ana_psi%GridDone = .TRUE.
  !-----------------------------------------------------------

  !-----------------------------------------------------------
   With_field = present(para_field)
   IF (present(para_field)) THEN
     CALL sub_dnE(dnE,0,T,para_field)
   ELSE
     dnE = (/ZERO,ZERO,ZERO/)
   END IF
  !-----------------------------------------------------------

  !-----------------------------------------------------------
  w2 = WP(1) ! just for the initialization
  IF(keep_MPI) CALL alloc_psi(w2,GridRep=.TRUE.)

  DO i=1,nb_WP
    !-----------------------------------------------------------
    ! => Analysis for diabatic potential (always done)

    ! =>first the energy
    IF(keep_MPI) THEN
      w1 = WP(i)
      CALL norm2_psi(w1,GridRep=.FALSE.,BasisRep=.TRUE.)
    ENDIF
    CALL sub_PsiOpPsi(ET,w1,w2,para_H)

    IF (para_H%spectral_done) THEN
      w2 = WP(i) ! save Spectral rep
      ! WP is back to the initial basis to be able to make the analysis
      IF (para_H%cplx) THEN
        WP(i)%CvecB(:) = matmul(conjg(para_H%Cvp),WP(i)%CvecB)
      ELSE
        WP(i)%CvecB(:) = matmul(para_H%Rvp,WP(i)%CvecB)
      END IF
    END IF

    IF(keep_MPI) CALL sub_PsiBasisRep_TO_GridRep(WP(i))


    IF(keep_MPI) THEN
      WP(i)%CAvOp = ET/w1%norm2

      G  = (para_propa%WPpsi  .AND.para_propa%write_GridRep)
      G2 = (para_propa%WPpsi2 .AND.para_propa%write_GridRep)
      B  = (para_propa%WPpsi  .AND.para_propa%write_BasisRep)
      B2 = (para_propa%WPpsi2 .AND.para_propa%write_BasisRep)

      CALL modif_ana_psi(para_propa%ana_psi,                            &
                         T=T,num_psi=i,Ene=real(WP(i)%CAvOp,kind=Rkind),&
                         With_field=With_field,field=dnE,               &
                         Write_psi_Grid=G,  Write_psi2_Grid=G2,         &
                         Write_psi_Basis=B, Write_psi2_Basis=B2)

      ! => The analysis (diabatic)
      CALL sub_analyze_psi(WP(i),para_propa%ana_psi,adia=.FALSE.)

      IF (para_propa%ana_psi%ExactFact > 0) THEN
        w1 = WP(i)
        write(out_unitp,*) i,'Exact Factorization analysis at ',T, ' ua'
        IF (present(para_field)) THEN
          CALL sub_ExactFact_analysis(T,w1,para_propa%ana_psi,para_H,   &
                       para_propa%WPTmax,para_propa%WPdeltaT,para_field)
        ELSE
          CALL sub_ExactFact_analysis(T,w1,para_propa%ana_psi,para_H,  &
                                 para_propa%WPTmax,para_propa%WPdeltaT)
        END IF
      END IF

      IF (para_propa%ana_psi%AvHiterm) THEN
        w1   = WP(i)
        info = real_TO_char(T,Rformat='f12.2')
        CALL sub_psiHitermPsi(w1,i,info,para_H)
      END IF

      ! => The analysis (adiabatic)
      IF (adia_loc) THEN
        w1 = WP(i)
        CALL sub_PsiDia_TO_PsiAdia_WITH_MemGrid(w1,para_H)
        para_propa%ana_psi%GridDone = .TRUE.
        CALL sub_analyze_psi(w1,para_propa%ana_psi,adia=.TRUE.)
      END IF
      CALL flush_perso(out_unitp)

      CALL alloc_psi(WP(i),BasisRep=BasisRep,GridRep=GridRep)
    ENDIF ! for keep_MPI

    IF (para_H%spectral_done) THEN
      WP(i) = w2 ! restore Spectral rep
    END IF

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

SUBROUTINE sub_analyze_mini_WP_OpWP(T,WP,nb_WP,para_H,adia,para_field)
  USE mod_system
  USE mod_Op,    ONLY : param_Op,sub_PsiOpPsi,sub_PsiDia_TO_PsiAdia_WITH_MemGrid
  USE mod_field, ONLY : param_field,sub_dnE

  USE mod_psi,   ONLY : param_psi,ecri_psi,alloc_psi,dealloc_psi,       &
                      sub_analyze_psi,norm2_psi,alloc_psi,modif_ana_psi,&
                      Channel_weight,sub_PsiBasisRep_TO_GridRep
  USE mod_MPI_aux
  IMPLICIT NONE

  real (kind=Rkind) :: T      ! time

!----- variables pour la namelist minimum ----------------------------
  TYPE (param_Op)   :: para_H

!- variables for the WP propagation ----------------------------
  TYPE (param_field), optional :: para_field
  logical,            optional :: adia

  integer            :: nb_WP
  TYPE (param_psi)   :: WP(:)

!-- working parameters --------------------------------
  TYPE (param_psi)   :: w1,w2

  integer       :: i,i_bi,i_be

  complex (kind=Rkind)              :: ET  ! energy
  real (kind=Rkind)                 :: E
  TYPE(REAL_WU)                     :: RWU_E
  character(len=:), allocatable     :: psi_line


  real (kind=Rkind), allocatable    :: tab_WeightChannels(:,:)
  real (kind=Rkind)                 :: Psi_norm2
  logical                           :: With_ENE = .TRUE.
  !logical                           :: With_ENE = .FALSE.

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
   write(out_unitp,*)
   write(out_unitp,*) 'nb_ba,nb_qa',WP(1)%nb_ba,WP(1)%nb_qa
   write(out_unitp,*) 'nb_bi',WP(1)%nb_bi
   write(out_unitp,*)

   DO i=1,nb_WP
     CALL norm2_psi(WP(i),.FALSE.,.TRUE.,.FALSE.)
     write(out_unitp,*) 'norm2psi0 BasisRep',i,WP(1)%norm2

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
  w2 = WP(1) ! just for the initialization
  CALL alloc_psi(w2)

  DO i=1,nb_WP

    w1 = WP(i)
    !-----------------------------------------------------------
    ! => Analysis for diabatic potential (always done)

    IF(keep_MPI) CALL Channel_weight(tab_WeightChannels,w1,GridRep=.FALSE.,BasisRep=.TRUE.)
    Psi_norm2 = sum(tab_WeightChannels)
    IF(openmpi .AND. MPI_scheme/=1) CALL MPI_Bcast_(Psi_norm2,size1_MPI,root_MPI)

    ! add the psi number + the time
    psi_line = 'norm^2-WP #WP ' // int_TO_char(i) // ' ' // real_TO_char(T,Rformat='f12.2')

    IF (With_ENE) THEN
      ! =>first the energy
      CALL sub_PsiOpPsi(ET,w1,w2,para_H)
      IF(openmpi .AND. MPI_scheme/=1) CALL MPI_Bcast_(ET,size1_MPI,root_MPI)
      w1%CAvOp = ET/Psi_norm2

      RWU_E  = REAL_WU(real(w1%CAvOp,kind=Rkind),'au','E')
      E      = convRWU_TO_R_WITH_WritingUnit(RWU_E)

      ! add the energy
      psi_line = psi_line // ' ' // real_TO_char(E,Rformat='f8.5')

    ELSE
      ! add the energy
      psi_line = psi_line // ' ' // 'xxxxxxxx'
    END IF

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


  END DO

  CALL dealloc_psi(w2,delete_all=.TRUE.)
  CALL dealloc_psi(w1,delete_all=.TRUE.)

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
      USE mod_psi,     ONLY : alloc_param_WP0
      IMPLICIT NONE

!------ parameter for the propagation ---------------------
      TYPE (param_propa) :: para_propa
      integer            :: nb_act1,nb_vp_spec_out

!----- physical and mathematical constants ---------------------------
      real (kind=Rkind) :: auTOcm_inv


!----- time parameters -----------------------------------------------
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
      logical       :: post_control,gate,cplx_gate,One_Iteration
      logical       :: Krotov,Turinici,envelopp,Obj_TO_alpha
      real (kind=Rkind),pointer    :: MRgate(:,:)
      complex (kind=Rkind),pointer :: MCgate(:,:)
      integer       :: err
!------ working variables ---------------------------------
      integer   :: i

      NAMELIST /propa/WPTmax,WPdeltaT,nb_micro,restart,                 &
                      name_WPpropa,type_WPpropa,spectral,nb_vp_spec,    &
                      One_Iteration,                                    &
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

      write(out_unitp,*) ' PROPAGATION PARAMETERS: propa, defWP0, control'
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
        One_Iteration       = .FALSE.

        DHmax               = -TEN
        auto_Hmax           = .FALSE.
        max_poly            = 5000
        npoly               = 0
        poly_tol            = ZERO

        type_WPpropa        = 0
        name_WPpropa        = ''
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


        para_propa%WPTmax                 = convRWU_TO_R_WITH_WorkingUnit(WPTmax)
        para_propa%WPdeltaT               = convRWU_TO_R_WITH_WorkingUnit(WPdeltaT)
        para_propa%nb_micro               = nb_micro
        para_propa%One_Iteration          = One_Iteration
        para_propa%para_poly%max_poly     = max_poly
        para_propa%para_poly%npoly        = npoly

        CALL alloc_param_poly(para_propa%para_poly)

        para_propa%spectral      = spectral
        nb_vp_spec_out           = nb_vp_spec
        para_propa%with_field    = .FALSE.

  IF (type_WPpropa > 0 .AND. name_WPpropa /= '') THEN
    write(out_unitp,*) 'ERROR in ',name_sub
    write(out_unitp,*) '  type_WPpropa and name_WPpropa are defined.'
    write(out_unitp,*) '  type_WPpropa: ',type_WPpropa
    write(out_unitp,*) '  name_WPpropa: ',trim(name_WPpropa)
    write(out_unitp,*) '  You have to chose only one.'
    write(out_unitp,*) '  => check your data!!'
    STOP
  END IF
  IF (type_WPpropa == 0) THEN
     CALL string_uppercase_TO_lowercase(name_WPpropa)
     SELECT CASE (name_WPpropa)
       CASE ('cheby','chebychev')
         type_WPpropa = 1
       CASE ('nod','taylor')
         type_WPpropa = 2
       CASE ('emin','emin-relax','relax')
         type_WPpropa = 3
       CASE ('emax')
         type_WPpropa = -3
       CASE ('davidson','emin-davidson')
         type_WPpropa = 33
       CASE ('emax-davidson')
         type_WPpropa = -33
       CASE ('cg','conjugated gradient')
         type_WPpropa = 34
       CASE ('tdh-nod','tdh-taylor')
         type_WPpropa = 22
       CASE ('tdh-rk2')
         type_WPpropa = 52
       CASE ('tdh-rk4')
         type_WPpropa = 54
       CASE ('rk4')
         type_WPpropa = 5

       CASE ('modmidpoint')
         type_WPpropa = 6
       CASE ('bulirsch-stoer')
         type_WPpropa = 7
       CASE ('sil')
         type_WPpropa = 8
       CASE ('sip')
         type_WPpropa = 9
       CASE ('spectral')
         type_WPpropa = 10

       CASE ('control','opt-control')
         type_WPpropa = 24
       CASE ('test')
         type_WPpropa = 100
       END SELECT
  END IF
  SELECT CASE (type_WPpropa)
  CASE (1) !         Chebychev
          IF (poly_tol .EQ. ZERO) poly_tol = ONETENTH**12
          IF (DHmax .EQ. -TEN) DHmax = HALF
          name_WPpropa  = 'cheby'
  CASE (2)!         Taylor
          IF (poly_tol .EQ. ZERO) poly_tol = ONETENTH**20
          IF (DHmax .EQ. -TEN) DHmax = ZERO
          name_WPpropa  = 'nOD'

  CASE (5) !         RK4 without a time dependant pulse in Hamiltonian (W(t))
          IF (poly_tol .EQ. ZERO) poly_tol = ONETENTH**20
          IF (DHmax .EQ. -TEN) DHmax = ZERO
          name_WPpropa  = 'RK4'
          para_propa%with_field    = .FALSE.

  CASE (6) !         ModMidPoint without a time dependant pulse in Hamiltonian (W(t))
          IF (poly_tol .EQ. ZERO) poly_tol = ONETENTH**20
          IF (DHmax .EQ. -TEN) DHmax = ZERO
          name_WPpropa  = 'ModMidPoint'
          para_propa%with_field    = .FALSE.

  CASE (7) !         Bulirsch-Stoer without a time dependant pulse in Hamiltonian (W(t))
          IF (poly_tol .EQ. ZERO) poly_tol = ONETENTH**20
          IF (DHmax .EQ. -TEN) DHmax = ZERO
          name_WPpropa  = 'Bulirsch-Stoer'
          para_propa%with_field    = .FALSE.

  CASE (8)!         Short Interative Propagation (Lanczos/Davidson)
          IF (poly_tol .EQ. ZERO) poly_tol = ONETENTH**20
          IF (DHmax .EQ. -TEN) DHmax = ZERO
          name_WPpropa  = 'SIL'
          para_propa%with_field    = .FALSE.

  CASE (9)!         Short Interative Propagation (Lanczos/Davidson)
          IF (poly_tol .EQ. ZERO) poly_tol = ONETENTH**20
          IF (DHmax .EQ. -TEN) DHmax = ZERO
          name_WPpropa  = 'SIP'
          para_propa%with_field    = .FALSE.

  CASE (10)!         Spectral representation Propagation
          poly_tol                 = ZERO
          DHmax                    = ZERO
          name_WPpropa             = 'Spectral'
          para_propa%with_field    = .FALSE.


  CASE (3,-3) !         im Relax
          IF (poly_tol .EQ. ZERO) poly_tol = ONETENTH**5
          IF (DHmax .EQ. -TEN) DHmax = ZERO
          IF (type_WPpropa ==  3) name_WPpropa  = 'Emin'
          IF (type_WPpropa == -3) name_WPpropa  = 'Emax'
  CASE (33,-33) !         Davidson
          IF (poly_tol .EQ. ZERO) poly_tol = ONETENTH**5
          IF (DHmax .EQ. -TEN) DHmax = ZERO
          name_WPpropa  = 'Davidson'
  CASE (34) ! Conjugated Gradient
          IF (poly_tol .EQ. ZERO) poly_tol = ONETENTH**5
          IF (DHmax .EQ. -TEN) DHmax = ZERO
          name_WPpropa  = 'Conjugated Gradient'


  CASE (22,221,222,223,24,241,242,243) !         nOD with a time dependant pulse in Hamiltonian (W(t))
          IF (poly_tol .EQ. ZERO) poly_tol = ONETENTH**20
          IF (DHmax .EQ. -TEN) DHmax = ZERO
          name_WPpropa  = 'TDH-nOD'
          para_propa%with_field    = .TRUE.
  CASE (50,54,52) !         RK4 with a time dependant pulse in Hamiltonian (W(t))
          poly_tol                 = ZERO
          IF (DHmax .EQ. -TEN) DHmax = ZERO
          name_WPpropa  = 'TDH_RK4'
          IF (type_WPpropa==52) name_WPpropa  = 'TDH_RK2'
          para_propa%with_field    = .TRUE.


  CASE (100) !         test
          IF (poly_tol .EQ. ZERO) poly_tol = ONETENTH**20
          IF (DHmax .EQ. -TEN) DHmax = ZERO
          name_WPpropa  = 'test'
          para_propa%with_field    = .TRUE.
  CASE Default
           write(out_unitp,*) 'ERROR in ',name_sub
           write(out_unitp,*) ' type of propagation(',type_WPpropa,')'
           write(out_unitp,*) ' is not possible'
           STOP
  END SELECT


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


        para_propa%TFmaxE       = convRWU_TO_R_WITH_WorkingUnit(TFmaxE)
        para_propa%TFminE       = convRWU_TO_R_WITH_WorkingUnit(TFminE)

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

          IF (Tenvelopp == ZERO) Tenvelopp = convRWU_TO_R_WITH_WorkingUnit(WPTmax)
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
      integer :: save_interal !< save WP every 'save_interal' step
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
                         max_poly,poly_tol,save_interal

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
      save_interal      = 1

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

      read(in_unitp,davidson)
      IF (print_level > 0) write(out_unitp,davidson)

      IF (.NOT. lower_states .AND. .NOT. all_lower_states .AND. &
          .NOT. project_WP0) lower_states = .TRUE.


      IF (residual_max_nb == huge(1) .AND. one_residue) residual_max_nb = 1

      ! for the filter diagonalization
      IF (DeltaM_filter == 0) DeltaM_filter = L_filter
      para_Davidson%E0_filter     = convRWU_TO_R_WITH_WorkingUnit(E0_filter)
      para_Davidson%W_filter      = convRWU_TO_R_WITH_WorkingUnit(W_filter)

      para_Davidson%L_filter      = L_filter
      para_Davidson%M_filter      = M_filter
      para_Davidson%DeltaM_filter = DeltaM_filter
      para_Davidson%Mmax_filter   = Mmax_filter

      IF (poly_tol == ZERO)    poly_tol = ONETENTH**8
      IF (DHmax    == -TEN)       DHmax = HALF
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
      para_Davidson%Max_ene               = convRWU_TO_R_WITH_WorkingUnit(Max_ene)

      para_Davidson%project_WP0           = project_WP0
      para_Davidson%residual_max_nb       = residual_max_nb
      para_Davidson%one_by_one            = one_by_one
      para_Davidson%NewVec_type           = NewVec_type
      para_Davidson%With_Grid             = With_Grid


      para_Davidson%symab             = symab
      para_Davidson%conv_ene          = convRWU_TO_R_WITH_WorkingUnit(conv_ene)
      para_Davidson%conv_resi         = convRWU_TO_R_WITH_WorkingUnit(conv_resi)

      para_Davidson%max_it            = max_it

      para_Davidson%nb_WP             = nb_WP
      para_Davidson%max_WP            = max_WP
      para_Davidson%conv_hermitian    = conv_hermitian

      para_Davidson%save_all          = save_all
      para_Davidson%save_interal      = save_interal
      para_Davidson%scaled_max_ene    = scaled_max_ene
      para_Davidson%save_max_nb       = save_max_nb
      para_Davidson%name_file_saveWP  = make_FileName(name_file_saveWP)

      IF(nb_WP/=0) THEN
        MPI_nb_WP=nb_WP
      ELSE 
        MPI_nb_WP=100
      ENDIF

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
