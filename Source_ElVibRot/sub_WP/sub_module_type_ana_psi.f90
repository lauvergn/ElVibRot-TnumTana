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

      MODULE mod_type_ana_psi
      USE mod_system
      IMPLICIT NONE

        TYPE param_ana_psi

          logical                        :: ana           = .TRUE.  ! analyze Psi
          integer                        :: num_psi       = 0       ! The numbering of psi (or wp)
          logical                        :: GridDone      = .FALSE. ! flag for the WP on the grid

          ! Boltzman population analysis
          logical                        :: Boltzmann_pop  = .TRUE.       ! population analysis
          real (kind=Rkind)              :: Temp           = -ONE         ! temperature
          real (kind=Rkind)              :: Part_func      = -ONE         ! partition function
          real (kind=Rkind)              :: ZPE            = -Huge(ONE)   ! ZPE
          real (kind=Rkind)              :: Ene            = -Huge(ONE)   ! Ene

          ! population analysis + reduced density
          TYPE (param_file)              :: file_PsiRho
          real (kind=Rkind), allocatable :: max_RedDensity(:)         ! the gaussian weighted (almost the last basis function) reduced density
          real (kind=Rkind), allocatable :: tab_WeightChannels(:,:)   ! tab_WeightChannels(nb_bi,nb_be)
          real (kind=Rkind)              :: Psi_norm2     = -ONE      ! norm^2 of psi
          logical                        :: adia          = .FALSE.   ! To perform special analysis with adiabatic states
          logical                        :: Rho1D         = .FALSE.   ! reduced densities (1D) along coordinate
          logical                        :: Rho2D         = .FALSE.   ! reduced densities (2D) along coordinate
          integer,           allocatable :: Weight_Rho(:)             ! enable to use a weight (0=>constant=1, +/-1=>step ...)
          real (kind=Rkind), allocatable :: Qana_weight(:)            ! geometry (Qact order???) for the analysis (use with Weight_Rho)
          integer                        :: Rho_type      = 2         ! Change the weight of rho along the coordiantes:
                                                                      ! Rho_type=0, without rho(Q).w(Q) => dT=dQ
                                                                      ! Rho_type=1, with    rho(Q).w(Q) => dT=rho(Q).w(Q)dQ (=>sum Rho(Q)=1)
                                                                      ! Rho_type=2, without w(Q);       => dT=w(Q)dQ


          logical                        :: AvQ           = .FALSE.   ! Average Values (Qact ...)
          integer,           allocatable :: Qtransfo_type(:)          ! ???type of the transformation


          logical                        :: AvScalOp      = .FALSE.   ! Average value of the scalar operators


          ! For quantum coherence Mij = Int [rho_i(Q)*rho_j(Q)/rho(Q) dQ]
          integer                        :: coherence = 0       ! default, 0: no coherence calculation
                                                                !          1: coherence calculation
          real (kind=Rkind)              :: coherence_epsi = ONETENTH**6  ! To avoid numerical trouble when rho(Q) is almost zero

          ! For the exact factorization
          integer                        :: ExactFact = 0       ! default, 0: no exact factorization analysis
          TYPE (param_file)              :: file_ExactFactPotCut


          ! 1D and 2D cut of psi at Qana
          TYPE (param_file)              :: file_PsiCut
          logical                        :: psi1D_Q0 = .FALSE.  ! reduced densities (1D or 2D) along coordinates
          logical                        :: psi2D_Q0 = .FALSE.  ! reduced densities (1D or 2D) along coordinates
          real (kind=Rkind), allocatable :: Qana(:)             ! geometry (Qact order) for the analysis

          ! write the psi
          TYPE (param_file)              :: file_Psi
          logical                        :: Write_psi        = .FALSE. ! write psi (or psi^2 ...)
          logical                        :: Write_psi2_Grid  = .FALSE. ! write the density on the grid
          logical                        :: Write_psi2_Basis = .FALSE. ! write the density on the basis functions
          logical                        :: Write_psi_Grid  = .FALSE.  ! write psi on the grid
          logical                        :: Write_psi_Basis = .FALSE.  ! write psi on the basis

          ! Propagation parameters
          logical                        :: propa         = .FALSE.      ! To perform special analysis of WP
          logical                        :: With_field    = .FALSE.      ! when the field is present
          real (kind=Rkind)              :: T             = ZERO         ! time in au
          real (kind=Rkind)              :: max_Psi_norm2 =  1.3_Rkind   ! if norm^2 >max_norm^2 => stop the WP propagation
          real (kind=Rkind)              :: field(3)= (/ZERO,ZERO,ZERO/) ! Electric field

        END TYPE param_ana_psi

        INTERFACE assignment (=)
          MODULE PROCEDURE ana_psi2_TO_ana_psi1
        END INTERFACE

        CONTAINS

      SUBROUTINE init_ana_psi(ana_psi,ana,num_psi,adia,           &
                              Write_psi2_Grid,Write_psi2_Basis,   &
                              Write_psi_Grid,Write_psi_Basis,     &
                              AvQ,Qtransfo_type,                  &
                              AvScalOp,                           &
                              coherence,coherence_epsi,           &
                              ExactFact,                          &
                              Rho1D,Rho2D,Weight_Rho,Qana_Weight, &
                              Rho_type,                           &
                              psi1D_Q0,psi2D_Q0,Qana,             &
                              propa,T,                            &
                              Ene,Boltzmann_pop,Temp,Part_func,ZPE)
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_ana_psi), intent(inout)      :: ana_psi
      logical,                        optional :: ana
      integer,                        optional :: num_psi
      logical,                        optional :: adia

      logical,                        optional :: Write_psi2_Grid,Write_psi2_Basis
      logical,                        optional :: Write_psi_Grid,Write_psi_Basis

      logical,                        optional :: AvQ,AvScalOp
      integer,           allocatable, optional :: Qtransfo_type(:)     ! type of the transformation

      integer,                        optional :: coherence         ! coherence_tyep (0 non calculation)
      real (kind=Rkind),              optional :: coherence_epsi    ! To avoid numerical trouble when rho(Q) is almost zero

      integer,                        optional ::ExactFact          ! for the exact factorization analysis

      logical,                        optional :: Rho1D,Rho2D
      integer,           allocatable, optional :: Weight_Rho(:)        ! enable to use a weight (0=>constant=1, +/-1=>step ...)
      real (kind=Rkind), allocatable, optional :: Qana_Weight(:)       ! geometry (Qact order) for the analysis (use with Weight_Rho)
      integer,                        optional :: Rho_type

      logical,                        optional :: psi1D_Q0,psi2D_Q0
      real (kind=Rkind), allocatable, optional :: Qana(:)              ! geometry (Qact order) for the analysis

      real (kind=Rkind),              optional :: Ene

      logical,                        optional :: propa
      real (kind=Rkind),              optional :: T ! time


      logical,                        optional :: Boltzmann_pop
      real (kind=Rkind),              optional :: Temp,Part_func,ZPE

      logical :: ana_weight

      character (len=*), parameter :: name_sub='init_ana_psi'



      !------------------------------------------------------------
      ana_psi%ana           = .TRUE. ! analyze Psi
      IF (present(ana))         ana_psi%ana         = ana
      ana_psi%num_psi       = 0      ! The numbering of psi (or wp)
      IF (present(num_psi))     ana_psi%num_psi     = num_psi
      !------------------------------------------------------------

      !------------------------------------------------------------
      ! write the psi
      ana_psi%Write_psi2_Grid  = .FALSE. ! write the density on the grid
      ana_psi%Write_psi2_Basis = .FALSE. ! write the density on the basis functions
      ana_psi%Write_psi_Grid   = .FALSE.  ! write psi on the grid
      ana_psi%Write_psi_Basis  = .FALSE.  ! write psi on the basis
      IF (present(Write_psi2_Grid))       ana_psi%Write_psi2_Grid       = Write_psi2_Grid
      IF (present(Write_psi2_Basis))      ana_psi%Write_psi2_Basis      = Write_psi2_Basis
      IF (present(Write_psi_Grid))        ana_psi%Write_psi_Grid        = Write_psi_Grid
      IF (present(Write_psi_Basis))       ana_psi%Write_psi_Basis       = Write_psi_Basis
      ana_psi%Write_psi = ana_psi%Write_psi2_Grid .OR. ana_psi%Write_psi2_Basis .OR.    &
                          ana_psi%Write_psi_Grid  .OR. ana_psi%Write_psi_Basis

      !------------------------------------------------------------
      ! for the reduced denstity analysis
      ana_psi%adia     = .FALSE.     ! To perform special analysis with adiabatic states
      IF (present(adia))        ana_psi%adia        = adia
      ana_psi%Rho1D    = .FALSE.
      IF (present(Rho1D))       ana_psi%Rho1D       = Rho1D
      ana_psi%Rho2D    = .FALSE.
      IF (present(Rho2D))       ana_psi%Rho2D       = Rho2D
      ana_psi%Rho_type = 2
      IF (present(Rho_type))       ana_psi%Rho_type       = Rho_type
      IF (ana_psi%Rho_type < 0 .OR. ana_psi%Rho_type > 2) ana_psi%Rho_type = 2

      IF (present(Qana_Weight) .AND. present(Weight_Rho)) THEN
        ana_weight = allocated(Qana_Weight) .AND. allocated(Weight_Rho)
      ELSE
        ana_weight = .FALSE.
      END IF

      IF (ana_weight) THEN

        IF (size(Qana_Weight) /= size(Weight_Rho) ) THEN
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) ' The size of Qana_Weight is different from the size of Weight_Rho'
          write(out_unitp,*) ' It is not possible!'
          write(out_unitp,*) ' Check the fortran!!'
          STOP
        END IF

        CALL alloc_NParray(ana_psi%Qana_Weight,shape(Qana_Weight),        &
                          "ana_psi%Qana_Weight",name_sub)
        ana_psi%Qana_Weight(:) = Qana_Weight(:)

        CALL alloc_NParray(ana_psi%Weight_Rho,shape(Weight_Rho),          &
                          "ana_psi%Weight_Rho",name_sub)
        ana_psi%Weight_Rho(:) = Weight_Rho(:)

      END IF
      !------------------------------------------------------------

      !------------------------------------------------------------
      ! analysis with 1D or 2D cuts
      ana_psi%psi1D_Q0 = .FALSE.
      ana_psi%psi2D_Q0 = .FALSE.
      IF (present(psi1D_Q0))    ana_psi%psi1D_Q0    = psi1D_Q0
      IF (present(psi2D_Q0))    ana_psi%psi2D_Q0    = psi2D_Q0
      IF (present(Qana)) THEN
        IF (allocated(Qana)) THEN
          CALL alloc_NParray(ana_psi%Qana,shape(Qana),"ana_psi%Qana",name_sub)
          ana_psi%Qana(:) = Qana(:)
        END IF
      END IF
      !------------------------------------------------------------

      !------------------------------------------------------------
      ! Average of Qi ...
      ana_psi%AvQ = .FALSE.
      IF (present(AvQ))         ana_psi%AvQ         = AvQ
      IF (present(Qtransfo_type)) THEN
      IF (allocated(Qtransfo_type)) THEN
        CALL alloc_NParray(ana_psi%Qtransfo_type,shape(Qtransfo_type),  &
                            "ana_psi%Qtransfo_type",name_sub)
        ana_psi%Qtransfo_type(:) = Qtransfo_type(:)
      END IF
      END IF
      !------------------------------------------------------------


      !------------------------------------------------------------
      ! Average of the Scalar Operators
      ana_psi%AvScalOp = .FALSE.
      IF (present(AvScalOp))         ana_psi%AvScalOp         = AvScalOp
      !------------------------------------------------------------


      !------------------------------------------------------------
      ! coherence calculation
      ana_psi%coherence = 0       ! default, 0: no coherence calculation
      IF (present(coherence)) ana_psi%coherence = coherence

      ana_psi%coherence_epsi = ONETENTH**6
      IF (present(coherence_epsi)) ana_psi%coherence_epsi = coherence_epsi
      !------------------------------------------------------------

      !------------------------------------------------------------
      ! exact factorization analysis
      ana_psi%ExactFact = 0       ! default, 0: no exact factorization analysis
      IF (present(ExactFact)) ana_psi%ExactFact = ExactFact
      !------------------------------------------------------------

      !------------------------------------------------------------
      ! Boltzman population analysis
      ana_psi%Ene       = -Huge(ONE) ! Ene
      IF (present(Ene))         ana_psi%Ene         = Ene

      ana_psi%Boltzmann_pop       = .TRUE. ! Boltzmann_pop
      IF (present(Boltzmann_pop))  ana_psi%Boltzmann_pop  = Boltzmann_pop
      ana_psi%Temp      = -ONE
      IF (present(Temp))        ana_psi%Temp        = Temp
      ana_psi%Part_func = -ONE
      IF (present(Part_func))   ana_psi%Part_func   = Part_func
      ana_psi%ZPE       = -Huge(ONE) ! ZPE
      IF (present(ZPE))         ana_psi%ZPE         = ZPE
      !------------------------------------------------------------

      !------------------------------------------------------------
      ana_psi%propa               = .FALSE.             ! To perform special analysis of WP
      ana_psi%T                   = ZERO                ! time in au
      ana_psi%Psi_norm2           = -ONE                ! norm^2 of psi
      ana_psi%max_Psi_norm2       =  1.3_Rkind          ! if norm^2 >max_norm^2 => stop the W
      ana_psi%field(:)            = (/ZERO,ZERO,ZERO/)  ! Electric field
      ana_psi%With_field          = .FALSE.             ! Propagation with a field
      !------------------------------------------------------------


      ana_psi%GridDone = .FALSE.

      END SUBROUTINE init_ana_psi
      SUBROUTINE dealloc_ana_psi(ana_psi)
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_ana_psi)        :: ana_psi


      IF (allocated(ana_psi%Qana_Weight)) THEN
        CALL dealloc_NParray(ana_psi%Qana_Weight,"ana_psi%Qana_Weight","dealloc_ana_psi")
      END IF
      IF (allocated(ana_psi%Weight_Rho)) THEN
        CALL dealloc_NParray(ana_psi%Weight_Rho,"ana_psi%Weight_Rho","dealloc_ana_psi")
      END IF
      IF (allocated(ana_psi%Qana)) THEN
        CALL dealloc_NParray(ana_psi%Qana,"ana_psi%Qana","dealloc_ana_psi")
      END IF

      IF (allocated(ana_psi%Qtransfo_type)) THEN
        CALL dealloc_NParray(ana_psi%Qtransfo_type,"ana_psi%Qtransfo_type","dealloc_ana_psi")
      END IF

      IF (allocated(ana_psi%max_RedDensity)) THEN
        CALL dealloc_NParray(ana_psi%max_RedDensity,"ana_psi%max_RedDensity","dealloc_ana_psi")
      END IF

     IF (allocated(ana_psi%tab_WeightChannels)) THEN
        CALL dealloc_NParray(ana_psi%tab_WeightChannels,                &
                            "ana_psi%tab_WeightChannels","dealloc_ana_psi")
      END IF

      CALL init_ana_psi(ana_psi)
      CALL file_dealloc(ana_psi%file_PsiRho)
      CALL file_dealloc(ana_psi%file_PsiCut)
      CALL file_dealloc(ana_psi%file_Psi)
      CALL file_dealloc(ana_psi%file_ExactFactPotCut)


      END SUBROUTINE dealloc_ana_psi

      SUBROUTINE ana_psi2_TO_ana_psi1(ana_psi1,ana_psi2)
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_ana_psi), intent(inout) :: ana_psi1
      TYPE (param_ana_psi), intent(in)    :: ana_psi2

     ana_psi1%ana           = ana_psi2%ana
     ana_psi1%num_psi       = ana_psi2%num_psi
     ana_psi1%GridDone      = ana_psi2%GridDone

     ana_psi1%Write_psi           = ana_psi2%Write_psi
     ana_psi1%Write_psi2_Grid     = ana_psi2%Write_psi2_Grid
     ana_psi1%Write_psi2_Basis    = ana_psi2%Write_psi2_Basis
     ana_psi1%Write_psi_Grid      = ana_psi2%Write_psi_Grid
     ana_psi1%Write_psi_Basis     = ana_psi2%Write_psi_Basis

     ana_psi1%Coherence_epsi      = ana_psi2%Coherence_epsi
     ana_psi1%Coherence           = ana_psi2%Coherence

     ana_psi1%ExactFact            = ana_psi2%ExactFact
     ana_psi1%file_ExactFactPotCut = ana_psi2%file_ExactFactPotCut

     ana_psi1%AvQ = ana_psi2%AvQ
     IF (allocated(ana_psi2%Qtransfo_type)) THEN
       CALL alloc_NParray(ana_psi1%Qtransfo_type,                       &
                                         shape(ana_psi2%Qtransfo_type), &
                         "ana_psi1%Qtransfo_type","ana_psi2_TO_ana_psi1")
       ana_psi1%Qtransfo_type(:) = ana_psi2%Qtransfo_type(:)
     END IF

     ana_psi1%AvScalOp      = ana_psi2%AvScalOp

     ana_psi1%adia          = ana_psi2%adia
     ana_psi1%Rho1D         = ana_psi2%Rho1D
     ana_psi1%Rho2D         = ana_psi2%Rho2D
     ana_psi1%Psi_norm2     = ana_psi2%Psi_norm2
     ana_psi1%file_PsiRho   = ana_psi2%file_PsiRho
     ana_psi1%Rho_type      = ana_psi2%Rho_type

     IF (allocated(ana_psi2%tab_WeightChannels)) THEN
        CALL alloc_NParray(ana_psi1%tab_WeightChannels,                      &
                                         shape(ana_psi2%tab_WeightChannels), &
                          "ana_psi1%tab_WeightChannels","ana_psi2_TO_ana_psi1")
        ana_psi1%tab_WeightChannels(:,:) = ana_psi2%tab_WeightChannels(:,:)
      END IF

      IF (allocated(ana_psi2%Qana_Weight)) THEN
        CALL alloc_NParray(ana_psi1%Qana_Weight,shape(ana_psi2%Qana_Weight),&
                          "ana_psi1%Qana_Weight","ana_psi2_TO_ana_psi1")
        ana_psi1%Qana_Weight(:) = ana_psi2%Qana_Weight(:)
      END IF
      IF (allocated(ana_psi2%Weight_Rho)) THEN
        CALL alloc_NParray(ana_psi1%Weight_Rho,shape(ana_psi2%Weight_Rho),&
                          "ana_psi1%Weight_Rho","ana_psi2_TO_ana_psi1")
        ana_psi1%Weight_Rho(:) = ana_psi2%Weight_Rho(:)
      END IF
     IF (allocated(ana_psi2%max_RedDensity)) THEN
        CALL alloc_NParray(ana_psi1%max_RedDensity,                     &
                                         shape(ana_psi2%max_RedDensity),&
                          "ana_psi1%max_RedDensity","ana_psi2_TO_ana_psi1")
        ana_psi1%max_RedDensity(:) = ana_psi2%max_RedDensity(:)
     END IF


      ana_psi1%psi1D_Q0     = ana_psi2%psi1D_Q0
      ana_psi1%psi2D_Q0     = ana_psi2%psi2D_Q0
      ana_psi1%file_PsiCut  = ana_psi2%file_PsiCut
      IF (allocated(ana_psi2%Qana)) THEN
        CALL alloc_NParray(ana_psi1%Qana,shape(ana_psi2%Qana),            &
                          "ana_psi1%Qana","ana_psi2_TO_ana_psi1")
        ana_psi1%Qana(:) = ana_psi2%Qana(:)
      END IF


      ana_psi1%Boltzmann_pop   = ana_psi2%Boltzmann_pop
      ana_psi1%Ene             = ana_psi2%Ene
      ana_psi1%Temp            = ana_psi2%Temp
      ana_psi1%Part_func       = ana_psi2%Part_func
      ana_psi1%ZPE             = ana_psi2%ZPE


     ana_psi1%propa         = ana_psi2%propa
     ana_psi1%T             = ana_psi2%T
     ana_psi1%max_Psi_norm2 = ana_psi2%max_Psi_norm2
     ana_psi1%field(:)      = ana_psi2%field(:)
     ana_psi1%With_field    = ana_psi2%With_field
     ana_psi1%file_Psi      = ana_psi2%file_Psi

      END SUBROUTINE ana_psi2_TO_ana_psi1

      SUBROUTINE Write_ana_psi(ana_psi)
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_ana_psi), intent(in)      :: ana_psi

      character (len=*), parameter :: name_sub='Write_ana_psi'

      write(out_unitp,*) 'BEGINNING ',name_sub

      write(out_unitp,*) 'ana',ana_psi%ana
      write(out_unitp,*) 'num_psi',ana_psi%num_psi
      write(out_unitp,*) 'GridDone',ana_psi%GridDone
       write(out_unitp,*)
      write(out_unitp,*) 'Boltzmann population:'
      write(out_unitp,*) 'Boltzmann_pop',ana_psi%Boltzmann_pop
      write(out_unitp,*) 'Temp',ana_psi%Temp
      write(out_unitp,*) 'Part_func',ana_psi%Part_func
      write(out_unitp,*) 'ZPE',ana_psi%ZPE
      write(out_unitp,*) 'Ene',ana_psi%Ene

      write(out_unitp,*)
      write(out_unitp,*) 'population analysis + reduced density:'
      write(out_unitp,*) 'file_PsiRho:'
      CALL file_Write(ana_psi%file_PsiRho)
      IF (allocated(ana_psi%max_RedDensity))                            &
              write(out_unitp,*) 'max_RedDensity',ana_psi%max_RedDensity
      IF (allocated(ana_psi%tab_WeightChannels))                        &
       write(out_unitp,*) 'tab_WeightChannels',ana_psi%tab_WeightChannels

      write(out_unitp,*) 'Psi_norm2',ana_psi%Psi_norm2
      write(out_unitp,*) 'adia',ana_psi%adia
      write(out_unitp,*) 'Rho1D',ana_psi%Rho1D
      write(out_unitp,*) 'Rho2D',ana_psi%Rho2D
      IF (allocated(ana_psi%Weight_Rho))                                &
                      write(out_unitp,*) 'Weight_Rho',ana_psi%Weight_Rho
      IF (allocated(ana_psi%Qana_weight))                               &
                   write(out_unitp,*) 'Qana_weight',ana_psi%Qana_weight
      write(out_unitp,*) 'Rho_type',ana_psi%Rho_type

      write(out_unitp,*)
      write(out_unitp,*) 'Average over coordinates:'
      write(out_unitp,*) 'AvQ',ana_psi%AvQ
      IF (allocated(ana_psi%Qtransfo_type))                             &
                write(out_unitp,*) 'Qtransfo_type',ana_psi%Qtransfo_type

      write(out_unitp,*)
      write(out_unitp,*) 'Average over scalar operators:'
      write(out_unitp,*) 'AvScalOp',ana_psi%AvScalOp

      write(out_unitp,*)
      write(out_unitp,*) 'Coherence:?'
      write(out_unitp,*) 'Coherence (Coherence_type)',ana_psi%Coherence
      write(out_unitp,*) 'Coherence_epsi',ana_psi%Coherence_epsi

      write(out_unitp,*)
      write(out_unitp,*) 'Exact Factorisation analysis:?'
      write(out_unitp,*) 'ExactFact',ana_psi%ExactFact


      write(out_unitp,*)
      write(out_unitp,*) '1D and 2D cut of psi at Qana:'
      write(out_unitp,*) 'file_PsiCut:'
      CALL file_Write(ana_psi%file_PsiCut)
      write(out_unitp,*) 'psi1D_Q0',ana_psi%psi1D_Q0
      write(out_unitp,*) 'psi2D_Q0',ana_psi%psi2D_Q0
      IF (allocated(ana_psi%Qana)) write(out_unitp,*) 'Qana',ana_psi%Qana


      write(out_unitp,*)
      write(out_unitp,*) 'Write_psi:',ana_psi%Write_psi
      write(out_unitp,*) 'file_Psi:'
      CALL file_Write(ana_psi%file_Psi)
      write(out_unitp,*) 'Write_psi2_Grid',ana_psi%Write_psi2_Grid
      write(out_unitp,*) 'Write_psi2_Basis',ana_psi%Write_psi2_Basis
      write(out_unitp,*) 'Write_psi_Grid',ana_psi%Write_psi_Grid
      write(out_unitp,*) 'Write_psi_Basis',ana_psi%Write_psi_Basis

      write(out_unitp,*)
      write(out_unitp,*) 'Propagation parameters:'
      write(out_unitp,*) 'propa',ana_psi%propa
      write(out_unitp,*) 'With_field',ana_psi%With_field
      write(out_unitp,*) 'T (time)',ana_psi%T
      write(out_unitp,*) 'max_Psi_norm2',ana_psi%max_Psi_norm2
      write(out_unitp,*) 'field(:)',ana_psi%field(:)

      END SUBROUTINE Write_ana_psi

      END MODULE mod_type_ana_psi

