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

      MODULE mod_analysis
      USE mod_system
      USE mod_file
      IMPLICIT NONE
        TYPE param_ana
          integer :: max_ana           = -1        ! nb of level to be analyzed. If max_ana=-1, all level are analyzed.
          real (kind=Rkind) :: max_ene = TEN**4    ! (cm-1)
          integer, allocatable :: Qtransfo_type(:) ! type of the transformation
          logical :: ana = .TRUE.                  ! analysis will be done

          integer :: print_psi                     ! nb of level write on the grid
          logical :: psi2                          ! print psi^2 instead of psi
          logical :: Rho1D,Rho2D                   ! reduced densities (1D or 2D) along coordinates
          integer, allocatable :: Weight_Rho(:)       ! enable to use a weight (0=>constant=1, +/-1=>step ...)
          real (kind=Rkind), allocatable :: Qana_Weight(:) ! geometry (Qact order) for the analysis (use with Weight_Rho)

          logical :: psi1D_Q0,psi2D_Q0             ! reduced densities (1D or 2D) along coordinates
          real (kind=Rkind), allocatable :: Qana(:)  ! geometry (Qact order) for the analysis

          logical :: Read_zpe                      ! (default F), T is Ezpe is read
          real (kind=Rkind) :: Ezpe                ! ground state reference energy

          logical :: davidson                      ! flag for Davidson diagonalization
          logical :: arpack                        ! flag for arpack diagonalization
          logical :: filter                        ! flag for filter diagonalization

          logical :: VibRot                        ! flag for RoVibrational levels (from the J=0 levels)
          integer :: JJmax = 10                    ! J values when VibRot=t

          logical :: intensity                     ! flag for intensity (life time) calculations
          logical :: NLO                           ! flag for NLO calculations
          logical :: Psi_ScalOp                    ! flag for NLO calculations
          integer :: CRP  = 0                      ! CRP=1, to CRP calculation
          real (kind=Rkind) :: CRP_Ene = ZERO      ! Total energy for CRP
          real (kind=Rkind) :: CRP_DEne = ZERO     ! Energy increment for the CRP
          integer :: nb_CRP_Ene  = 1               ! Number of CRP calculation

          real (kind=Rkind) :: Temp                ! temperature (K) (for intensity)

          logical :: print                         ! flag for priting H, Vec matrices
          logical :: Spectral_ScalOp               ! IF T, calculate the spectral representation of scalar Op
          integer :: MaxWP_TO_Write_MatOp = 100    ! when the WP is larger than this parameter,
                                                   ! the spectral representation of operator is not printed
          logical :: propa                         ! flag for WP propagation
          logical :: control                       ! flag for the optimal control

          character (len=Line_len) :: name_file_spectralWP = 'file_WPspectral'
          logical                  :: formatted_file_WP    = .TRUE.

        END TYPE param_ana

        TYPE param_intensity
          logical :: l_Aif        = .FALSE. ! Einstein coef (s-1)
          logical :: l_Int        = .TRUE.  ! intensity
          logical :: l_CrossSec   = .FALSE. ! Cross Section (m^2)
          logical :: l_IntVR      = .FALSE. ! intensity vib+rot
          logical :: l_Tau        = .FALSE. ! life time (s)

          logical :: pola_xyz(3)  = (/.TRUE.,.TRUE.,.TRUE./) ! use some of component of the dipole moment

          logical :: l_md         = .FALSE.
          logical :: l_md2        = .FALSE.

          logical :: dip_Qact     = .FALSE.

          logical                    :: Inertia   = .FALSE.
          real (kind=Rkind)          :: Temp      = 298_Rkind           ! temperature (K) (for intensity)
          real (kind=Rkind)          :: ABCref(3) = (/ZERO,ZERO,ZERO/)  ! rotational constant (au)
          real (kind=Rkind), allocatable :: ABC(:,:)

          integer           :: JJmax = 10
          real (kind=Rkind) :: Emin,Emax      ! in cm-1
          real (kind=Rkind) :: Ewidth         ! 1. in cm-1
          real (kind=Rkind) :: nEstep         ! 10
          logical :: l_lorentz = .TRUE.      ! (T) if F => gaussian

          TYPE (param_file) :: file_spectrum,file_intensity
          TYPE (param_file) :: file_resart_int

        END TYPE param_intensity
      CONTAINS
      SUBROUTINE read_analyse(para_ana,Qana)

      USE mod_system
      IMPLICIT NONE

!----- variables pour la namelist analyse ----------------------------
      TYPE (param_ana)  :: para_ana
      real (kind=Rkind) :: Qana(:)


      integer       :: nb_harm_ana,max_ana,print_psi,MaxWP_TO_Write_MatOp,JJmax
      logical       :: ana,print,propa,intensity,Psi_ScalOp
      logical       :: control,davidson,arpack,filter,NLO,VibRot
      integer       :: CRP,nb_CRP_Ene
      logical       :: Rho1D,Rho2D,psi2,psi1D_Q0,psi2D_Q0,Wheight_rho
      logical       :: QTransfo,formatted_file_WP
      logical        :: Spectral_ScalOp
      TYPE (REAL_WU) :: CRP_Ene,CRP_DEne,ene0,Ezpe,max_ene
      real (kind=Rkind) :: Temp

      character (len=Line_len) :: name_file_spectralWP
      character (len=name_len) :: name_dum

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = "read_analyse"
      logical, parameter :: debug=.FALSE.
!      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------

      NAMELIST /analyse/ana,print,nb_harm_ana,max_ana,max_ene,          &
                        propa,                                          &
                        print_psi,psi2,psi1D_Q0,psi2D_Q0,QTransfo,      &
                        Rho1D,Rho2D,Wheight_rho,                        &
                        intensity,NLO,CRP,CRP_Ene,CRP_DEne,nb_CRP_Ene,  &
                        Psi_ScalOp,VibRot,JJmax,                        &
                        ene0,Ezpe,Temp,                                 &
                name_file_spectralWP,formatted_file_WP,FilePsiVersion,  &
                        control,davidson,arpack,filter,Spectral_ScalOp, &
                        MaxWP_TO_Write_MatOp


      write(out_unitp,*) ' ANALYSIS PARAMETERS'

      print                = .FALSE.
      psi2                 = .FALSE.
      QTransfo             = .FALSE.
      Rho1D                = .FALSE.
      Rho2D                = .FALSE.
      Wheight_rho          = .FALSE.
      psi1D_Q0             = .FALSE.
      psi2D_Q0             = .FALSE.
      propa                = .FALSE.
      intensity            = .FALSE.
      Psi_ScalOp           = .FALSE.
      NLO                  = .FALSE.
      CRP                  = 0
      CRP_Ene              = REAL_WU(ZERO,'cm-1','E')
      CRP_DEne             = REAL_WU(ZERO,'cm-1','E')
      nb_CRP_Ene           = 1

      VibRot               = .FALSE.
      JJmax                = -1
      control              = .FALSE.
      davidson             = .FALSE.
      filter               = .FALSE.
      arpack               = .FALSE.
      Spectral_ScalOp      = .FALSE.
      MaxWP_TO_Write_MatOp = 100
      name_file_spectralWP = 'file_WPspectral'
      formatted_file_WP    = .TRUE.
      FilePsiVersion       = 0

      ana             = .TRUE.
      print_psi       = 0
      nb_harm_ana     = 1
      max_ana         = -1
      Temp            = -ONE
      ene0            = REAL_WU(huge(ONE),'cm-1','E')
      Ezpe            = REAL_WU(huge(ONE),'cm-1','E')
      max_ene         = REAL_WU(TEN**4,   'cm-1','E') ! 10 000 cm-1
      read(in_unitp,analyse)
      IF (print_level > 0) write(out_unitp,analyse)
      write(out_unitp,*)
      IF (davidson .AND. arpack) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' davidson=t and arpack=t'
        write(out_unitp,*) ' You HAVE to select only one!'
        STOP
      END IF
      IF (.NOT. formatted_file_WP) FilePsiVersion = 1
      write(out_unitp,*) 'name_file_spectralWP,formatted_file_WP: ', &
                          name_file_spectralWP,formatted_file_WP

      IF (Qtransfo) THEN
        CALL alloc_NParray(para_ana%Qtransfo_type,shape(Qana),"para_ana%Qtransfo_type",name_sub)
        read(in_unitp,*) name_dum,para_ana%Qtransfo_type(:)
      END IF

      IF (Wheight_rho) THEN
        CALL alloc_NParray(para_ana%Qana_Weight,shape(Qana),"para_ana%Qana_Weight",name_sub)
        read(in_unitp,*) name_dum,para_ana%Qana_Weight(:)

        CALL alloc_NParray(para_ana%Weight_Rho,shape(Qana),"para_ana%Weight_Rho",name_sub)
        read(in_unitp,*) name_dum,para_ana%Weight_Rho(:)
      END IF

      Spectral_ScalOp = Spectral_ScalOp .OR. intensity .OR. NLO .OR. Psi_ScalOp



      IF (ene0%val /= huge(ONE) .AND. Ezpe%val == huge(ONE) ) Ezpe = ene0
      para_ana%Read_zpe = (Ezpe%val /= huge(ONE))
      IF (max_ana < 0) max_ana = huge(1)

      para_ana%print                = print
      para_ana%propa                = propa
      para_ana%control              = control
      para_ana%davidson             = davidson
      para_ana%arpack               = arpack
      para_ana%filter               = filter
      para_ana%Spectral_ScalOp      = Spectral_ScalOp
      para_ana%MaxWP_TO_Write_MatOp = MaxWP_TO_Write_MatOp
      para_ana%name_file_spectralWP = name_file_spectralWP
      para_ana%formatted_file_WP    = formatted_file_WP

      para_ana%max_ene         = convRWU_TO_R(max_ene)
      para_ana%max_ana         = max_ana

      para_ana%ana             = ana
      para_ana%psi2            = psi2
      para_ana%Rho1D           = Rho1D
      para_ana%Rho2D           = Rho2D
      para_ana%psi1D_Q0        = psi1D_Q0
      para_ana%psi2D_Q0        = psi2D_Q0

      CALL alloc_NParray(para_ana%Qana,shape(Qana),"para_ana%Qana",name_sub)
      para_ana%Qana(:) = Qana(:)

      para_ana%print_psi       = print_psi
      para_ana%intensity       = intensity
      para_ana%Psi_ScalOp      = Psi_ScalOp
      para_ana%NLO             = NLO

      para_ana%CRP             = CRP
      para_ana%CRP_Ene         = convRWU_TO_R(CRP_Ene)
      para_ana%CRP_DEne        = convRWU_TO_R(CRP_DEne)
      para_ana%nb_CRP_Ene      = nb_CRP_Ene


      write(out_unitp,*) 'CRP,E,DE,nb_E   : ',para_ana%CRP,             &
                  para_ana%CRP_Ene,para_ana%CRP_DEne,para_ana%nb_CRP_Ene


      para_ana%Ezpe            = convRWU_TO_R(Ezpe)
      para_ana%Temp            = Temp

      para_ana%VibRot          = VibRot
      para_ana%JJmax           = JJmax
      IF (.NOT. VibRot) para_ana%JJmax = -1


      write(out_unitp,*) 'Ezpe   : ',RWU_Write(Ezpe,WithUnit=.TRUE.,WorkingUnit=.FALSE.)
      write(out_unitp,*) 'max_ene: ',RWU_Write(max_ene,WithUnit=.TRUE.,WorkingUnit=.FALSE.)


      END SUBROUTINE read_analyse

      SUBROUTINE dealloc_para_ana(para_ana)

      USE mod_system
      IMPLICIT NONE

!----- variables pour la namelist analyse ----------------------------
      TYPE (param_ana), intent(inout) :: para_ana

      IF (allocated(para_ana%Qtransfo_type)) THEN
        CALL dealloc_NParray(para_ana%Qtransfo_type,"para_ana%Qtransfo_type","dealloc_para_ana")
      END IF

      IF (allocated(para_ana%Qana_Weight)) THEN
        CALL dealloc_NParray(para_ana%Qana_Weight,"para_ana%Qana_Weight","dealloc_para_ana")
      END IF
      IF (allocated(para_ana%Weight_Rho)) THEN
        CALL dealloc_NParray(para_ana%Weight_Rho,"para_ana%Weight_Rho","dealloc_para_ana")
      END IF
      IF (allocated(para_ana%Qana)) THEN
        CALL dealloc_NParray(para_ana%Qana,"para_ana%Qana","dealloc_para_ana")
      END IF


      END SUBROUTINE dealloc_para_ana


      SUBROUTINE read_intensity(para_intensity)
      USE mod_system
      IMPLICIT NONE

      TYPE (param_intensity) :: para_intensity
      real (kind=Rkind) :: auTOcm_inv

      logical       :: l_Aif,l_Int,l_Tau,l_md,l_md2,l_IntVR,l_CrossSec
      logical       :: dip_Qact
      logical       :: pola_xyz(3)
      real (kind=Rkind) :: Temp,ABC(3)
      integer       :: JJmax

      logical       :: l_lorentz,Inertia
      TYPE(REAL_WU) :: Emin,Emax,Ewidth

      integer       :: nEstep

      character (len=Line_len) :: file_spectrum,file_intensity
      character (len=Line_len) :: file_resart_int

      integer       :: nb_t

      NAMELIST /intensity/Temp,l_Aif,l_Int,l_IntVR,l_Tau,l_md,l_md2,    &
                          l_CrossSec,                                   &
                          pola_xyz,                                     &
                          ABC,JJmax,Inertia,                            &
                          dip_Qact,                                     &
                          Emin,Emax,Ewidth,nEstep,l_lorentz,            &
                          file_spectrum,file_intensity,file_resart_int


!--------------------------------------------------------------------
      write(out_unitp,*) ' INTENSITY PARAMETERS'
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')

        l_Int       = .TRUE.
        l_CrossSec  = .FALSE.
        l_IntVR     = .FALSE.
        l_Aif       = .FALSE.
        l_Tau       = .FALSE.
        l_md        = .FALSE.
        l_md2       = .FALSE.
        dip_Qact    = .FALSE.
        pola_xyz(:) = .TRUE.
        Temp        = 298.15_Rkind
        ABC(:)      = ZERO
        JJmax       = 100
        Inertia     = .FALSE.

        l_lorentz   = .TRUE.
        Emin        = REAL_WU(ZERO,'cm-1','E')
        Emax        = REAL_WU(ZERO,'cm-1','E')
        Ewidth      = REAL_WU(ONE ,'cm-1','E')
        nEstep      = 10

        file_spectrum   = 'spectrum'
        file_intensity  = 'intensity'
        file_resart_int = 'restart.int'

        read(in_unitp,intensity)
        nb_t = 0
        IF (l_Int) nb_t = nb_t + 1
        IF (l_CrossSec) nb_t = nb_t + 1
        IF (l_IntVR) nb_t = nb_t + 1
        IF (l_Aif) nb_t = nb_t + 1
        IF (l_Tau) nb_t = nb_t + 1
        IF (nb_t /= 1) THEN
          write(out_unitp,*) ' ERROR in read_intensity'
          write(out_unitp,*) ' You have to chose ONE option:'
          write(out_unitp,*) '     l_Int=t or'
          write(out_unitp,*) '     l_CrossSec=t or'
          write(out_unitp,*) '     l_IntVR=t or'
          write(out_unitp,*) '     l_Aif=t or'
          write(out_unitp,*) '     l_Tau=t'
          write(out_unitp,*) ' but l_Int,l_IntVR,l_Aif,l_Tau:',         &
                                   l_Int,l_CrossSec,l_IntVR,l_Aif,l_Tau
          STOP
        END IF

        IF (print_level > 0) write(out_unitp,intensity)
        write(out_unitp,*)

        para_intensity%l_Int       = l_Int
        para_intensity%l_CrossSec  = l_CrossSec
        para_intensity%l_IntVR     = l_IntVR
        para_intensity%l_Aif       = l_Aif
        para_intensity%l_Tau       = l_Tau
        para_intensity%l_md        = l_md
        para_intensity%l_md2       = l_md2
        para_intensity%pola_xyz(:) = pola_xyz(:)

        para_intensity%Temp      = Temp
        para_intensity%ABCref(:) = ABC(:) / auTOcm_inv
        para_intensity%JJmax     = JJmax
        para_intensity%Inertia   = Inertia

        para_intensity%dip_Qact  = dip_Qact

        para_intensity%l_lorentz   = l_lorentz
        para_intensity%Emin        = convRWU_TO_R(Emin)
        para_intensity%Emax        = convRWU_TO_R(Emax)
        para_intensity%Ewidth      = convRWU_TO_R(Ewidth)
        para_intensity%nEstep      = nEstep


        para_intensity%file_spectrum%name   = file_spectrum
        para_intensity%file_intensity%name  = file_intensity
        para_intensity%file_resart_int%name = file_resart_int

      END SUBROUTINE read_intensity

      END MODULE mod_analysis

