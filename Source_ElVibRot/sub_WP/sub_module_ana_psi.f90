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

      MODULE mod_ana_psi
      USE mod_system
      IMPLICIT NONE

        TYPE param_ana_psi

          logical                        :: ana    = .TRUE.     ! analyze Psi
          logical                        :: AvQ    = .FALSE.    ! Average Values (Qact ...)
          logical                        :: Rho1D  = .FALSE.    ! reduced densities (1D or 2D) along coordinates
          logical                        :: Rho2D  = .FALSE.    ! reduced densities (1D or 2D) along coordinates
          integer,           allocatable :: Weight_Rho(:)       ! enable to use a weight (0=>constant=1, +/-1=>step ...)
          real (kind=Rkind), allocatable :: Qana_weight(:)      ! geometry (Qact order) for the analysis (use with Weight_Rho)

          integer,           allocatable :: Qtransfo_type(:)    ! type of the transformation

          logical                        :: psi1D_Q0 = .FALSE.  ! reduced densities (1D or 2D) along coordinates
          logical                        :: psi2D_Q0 = .FALSE.  ! reduced densities (1D or 2D) along coordinates
          real (kind=Rkind), allocatable :: Qana(:)             ! geometry (Qact order) for the analysis
          integer                        :: iwp    = 0

          real (kind=Rkind)              :: Temp = -ONE         ! temperature
          real (kind=Rkind)              :: Part_func = -ONE    ! partition function
          real (kind=Rkind)              :: ZPE = -Huge(ONE)    ! ZPE
          real (kind=Rkind)              :: Ene = -Huge(ONE)    ! Ene

          real (kind=Rkind), allocatable :: max_RedDensity(:)          ! the gaussian weighted (almost the last basis function) reduced density

          real (kind=Rkind)              :: Psi_norm2     = -ONE       ! norm^2 of psi
          real (kind=Rkind)              :: max_Psi_norm2 =  1.3_Rkind ! if norm^2 >max_norm^2 => stop the WP propagation
          real (kind=Rkind), allocatable :: tab_WeightChannels(:,:)    ! tab_WeightChannels(nb_bi,nb_be)

        END TYPE param_ana_psi

        INTERFACE assignment (=)
          MODULE PROCEDURE ana_psi2_TO_ana_psi1
        END INTERFACE

        CONTAINS

      SUBROUTINE init_ana_psi(ana_psi,ana,iwp,AvQ,                      &
                              Rho1D,Rho2D,Weight_Rho,Qana_Weight,       &
                              Qtransfo_type,                            &
                              psi1D_Q0,psi2D_Q0,Qana,                   &
                              Ene,Temp,Part_func,ZPE)
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_ana_psi)                     :: ana_psi
      logical,                        optional :: AvQ,ana,Rho1D,Rho2D
      integer,           allocatable, optional :: Weight_Rho(:)        ! enable to use a weight (0=>constant=1, +/-1=>step ...)
      real (kind=Rkind), allocatable, optional :: Qana_Weight(:)       ! geometry (Qact order) for the analysis (use with Weight_Rho)
      integer,           allocatable, optional :: Qtransfo_type(:)     ! type of the transformation
      logical,                        optional :: psi1D_Q0,psi2D_Q0
      real (kind=Rkind), allocatable, optional :: Qana(:)              ! geometry (Qact order) for the analysis
      integer,                        optional :: iwp
      real (kind=Rkind),              optional :: Ene,Temp,Part_func,ZPE

      logical :: ana_weight

      character (len=*), parameter :: name_sub='init_ana_psi'

      ana_psi%ana = .TRUE.
      ana_psi%AvQ = .FALSE.
      ana_psi%iwp = 0

      ana_psi%Rho1D    = .FALSE.
      ana_psi%Rho2D    = .FALSE.
      ana_psi%psi1D_Q0 = .FALSE.
      ana_psi%psi2D_Q0 = .FALSE.

      ana_psi%Ene       = -Huge(ONE) ! Ene
      ana_psi%Temp      = -ONE
      ana_psi%Part_func = -ONE
      ana_psi%ZPE       = -Huge(ONE) ! Ene


      IF (present(ana))         ana_psi%ana         = ana
      IF (present(iwp))         ana_psi%iwp         = iwp
      IF (present(Rho1D))       ana_psi%Rho1D       = Rho1D
      IF (present(Rho2D))       ana_psi%Rho2D       = Rho2D

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

      IF (present(AvQ))         ana_psi%AvQ         = AvQ
      IF (present(psi1D_Q0))    ana_psi%psi1D_Q0    = psi1D_Q0
      IF (present(psi2D_Q0))    ana_psi%psi2D_Q0    = psi2D_Q0
      IF (present(Qana)) THEN
        IF (allocated(Qana)) THEN
          CALL alloc_NParray(ana_psi%Qana,shape(Qana),"ana_psi%Qana",name_sub)
          ana_psi%Qana(:) = Qana(:)
        END IF
      END IF


      IF (present(Qtransfo_type)) THEN
        IF (allocated(Qtransfo_type)) THEN
          CALL alloc_NParray(ana_psi%Qtransfo_type,shape(Qtransfo_type),  &
                            "ana_psi%Qtransfo_type",name_sub)
          ana_psi%Qtransfo_type(:) = Qtransfo_type(:)
        END IF
      END IF

      IF (present(Ene))         ana_psi%Ene         = Ene

      IF (present(Temp))        ana_psi%Temp        = Temp
      IF (present(Part_func))   ana_psi%Part_func   = Part_func
      IF (present(ZPE))         ana_psi%ZPE         = ZPE

      ana_psi%Psi_norm2           = -ONE       ! norm^2 of psi
      ana_psi%max_Psi_norm2       =  1.3_Rkind ! if norm^2 >max_norm^2 => stop the W

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

      END SUBROUTINE dealloc_ana_psi

      SUBROUTINE ana_psi2_TO_ana_psi1(ana_psi1,ana_psi2)
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_ana_psi), intent(inout) :: ana_psi1
      TYPE (param_ana_psi), intent(in)    :: ana_psi2


      ana_psi1%ana = ana_psi2%ana
      ana_psi1%iwp = ana_psi2%iwp
      ana_psi1%AvQ = ana_psi2%AvQ

      ana_psi1%Rho1D    = ana_psi2%Rho1D
      ana_psi1%Rho2D    = ana_psi2%Rho2D
      ana_psi1%psi1D_Q0 = ana_psi2%psi1D_Q0
      ana_psi1%psi2D_Q0 = ana_psi2%psi2D_Q0

      ana_psi1%Ene       = ana_psi2%Ene
      ana_psi1%Temp      = ana_psi2%Temp
      ana_psi1%Part_func = ana_psi2%Part_func
      ana_psi1%ZPE       = ana_psi2%ZPE

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

      IF (allocated(ana_psi2%Qana)) THEN
        CALL alloc_NParray(ana_psi1%Qana,shape(ana_psi2%Qana),            &
                          "ana_psi1%Qana","ana_psi2_TO_ana_psi1")
        ana_psi1%Qana(:) = ana_psi2%Qana(:)
      END IF

      IF (allocated(ana_psi2%Qtransfo_type)) THEN
        CALL alloc_NParray(ana_psi1%Qtransfo_type,                        &
                                         shape(ana_psi2%Qtransfo_type), &
                          "ana_psi1%Qtransfo_type","ana_psi2_TO_ana_psi1")
        ana_psi1%Qtransfo_type(:) = ana_psi2%Qtransfo_type(:)
      END IF

     IF (allocated(ana_psi2%max_RedDensity)) THEN
        CALL alloc_NParray(ana_psi1%max_RedDensity,                     &
                                         shape(ana_psi2%max_RedDensity),&
                          "ana_psi1%max_RedDensity","ana_psi2_TO_ana_psi1")
        ana_psi1%max_RedDensity(:) = ana_psi2%max_RedDensity(:)
     END IF


     ana_psi1%Psi_norm2     = ana_psi2%Psi_norm2
     ana_psi1%max_Psi_norm2 = ana_psi2%max_Psi_norm2
     IF (allocated(ana_psi2%tab_WeightChannels)) THEN
        CALL alloc_NParray(ana_psi1%tab_WeightChannels,                      &
                                         shape(ana_psi2%tab_WeightChannels), &
                          "ana_psi1%tab_WeightChannels","ana_psi2_TO_ana_psi1")
        ana_psi1%tab_WeightChannels(:,:) = ana_psi2%tab_WeightChannels(:,:)
      END IF

      END SUBROUTINE ana_psi2_TO_ana_psi1

      END MODULE mod_ana_psi

