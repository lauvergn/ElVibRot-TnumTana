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

 MODULE mod_march
 USE mod_system
 USE mod_file
 USE mod_field,         ONLY : param_field
 USE mod_psi_set_alloc, ONLY : param_psi
 USE mod_propa,         ONLY : param_propa,cof,Calc_AutoCorr,Write_AutoCorr,SaveWP_restart
 USE mod_march_SG4
 use mod_Constant,      only: get_conv_au_to_unit
 IMPLICIT NONE

 PRIVATE
 PUBLIC :: march_gene,march_nOD_im,march_FOD_Opti_im

      CONTAINS
!================================================================
!
!    march_gene  with or without field
!
!================================================================
      SUBROUTINE march_gene(T,WP,WP0,nb_WP,print_Op,                    &
                            para_H,para_propa,tab_Op,para_field)
      USE mod_system
      USE mod_Op,              ONLY : param_Op
      USE mod_psi_set_alloc,   ONLY : param_psi,ecri_psi
      USE mod_ana_psi,         ONLY : norm2_psi
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)           :: para_H
      TYPE (param_Op), optional :: tab_Op(:)

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa)           :: para_propa
      TYPE (param_field), optional :: para_field

      integer          :: nb_WP
      TYPE (param_psi) :: WP(:),WP0(:)

!----- for printing --------------------------------------------------
      logical           :: print_Op
      real (kind=Rkind) :: T      ! time
      integer           :: i,no

!------ working parameters --------------------------------
      logical :: SGtype4 = .FALSE.
      !logical :: SGtype4 = .TRUE.

!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='march_gene'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Tmax,deltaT',para_propa%WPTmax,para_propa%WPdeltaT
        write(out_unitp,*) 'Hmin,Hmax',para_propa%para_poly%Hmin,       &
                                       para_propa%para_poly%Hmax
        write(out_unitp,*)
        write(out_unitp,*) 'nb_ba,nb_qa',WP%nb_ba,WP%nb_qa
        write(out_unitp,*) 'nb_bi',WP%nb_bi
        write(out_unitp,*)
        write(out_unitp,*) 'type_WPpropa',para_propa%type_WPpropa
        write(out_unitp,*)

        DO i=1,nb_WP
          CALL norm2_psi(WP(i),.FALSE.,.TRUE.,.FALSE.)
          write(out_unitp,*) 'norm WP',i,WP(i)%norme

          write(out_unitp,*) 'WP BasisRep',i
          CALL ecri_psi(T=ZERO,psi=WP(i),                               &
                        ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)
          write(out_unitp,*) 'WP GridRep',i
          CALL ecri_psi(T=ZERO,psi=WP(i),                               &
                        ecri_GridRep=.TRUE.,ecri_BasisRep=.FALSE.)
        END DO

      END IF
!-----------------------------------------------------------

      IF (para_propa%With_field .AND. .NOT. present(para_field)         &
          .AND. .NOT. present(tab_Op)) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'With_field=T and para_field and tab_Op are not present'
        write(out_unitp,*) 'check the source!!'
        STOP
      END IF

!-----------------------------------------------------------
      para_propa%march_error = .FALSE.
      no = para_propa%file_autocorr%unit

      SGtype4    = SGtype4 .AND. (para_H%BasisnD%SparseGrid_type == 4)

      SELECT CASE (para_propa%type_WPpropa)

      CASE (1)
        !CALL march_cheby_old(T,no,WP(1),WP0(1),para_H,para_propa)
        CALL march_cheby(T,no,WP(1),WP0(1),para_H,para_propa)

      CASE (2)
        !IF (SGtype4 .AND. direct_KEO) THEN
        IF (SGtype4) THEN
          CALL  march_noD_SG4(T,no,WP(1),WP0(1),para_H,para_propa)
        ELSE
          CALL  march_noD(T,no,WP(1),WP0(1),para_H,para_propa)
        END IF

      CASE (5)
        CALL march_RK4(T,no,WP(1),WP0(1),para_H,para_propa)

      CASE (6)
        CALL march_ModMidPoint(T,no,WP(1),WP0(1),para_H,para_propa)

      CASE (7)
        CALL march_BS(T,no,WP(1),WP0(1),para_H,para_propa)

      CASE (22,24) ! nOD
!       CALL march_split_field(T,WP(:),nb_WP,print_Op,
!    *                       para_H,tab_Op(3:5),
!    *                       para_field,para_propa)
!       CALL march_new_noD_field(T,WP(:),nb_WP,print_Op,
!    *                       para_H,tab_Op(3:5),
!    *                       para_field,para_propa)
        CALL march_noD_field(T,WP(:),nb_WP,print_Op,                    &
                             para_H,tab_Op(1:3),                        &
                             para_field,para_propa)

      CASE (50,54) ! RK4
        CALL march_RK4_field(T,WP(:),nb_WP,print_Op,                    &
                             para_H,tab_Op(1:3),                        &
                             para_field,para_propa)

      CASE (52) ! RK2
        CALL march_RK2_field(T,WP(:),nb_WP,print_Op,                    &
                             para_H,tab_Op(1:3),                        &
                             para_field,para_propa)

      CASE DEFAULT
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'The type_WPpropa is unknown',para_propa%type_WPpropa
        write(out_unitp,*) 'check your data !!'
        STOP

      END SELECT

      CALL SaveWP_restart(T+para_propa%WPdeltaT,WP,para_propa%file_WP_restart)

!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Psi at T+DT'
        DO i=1,nb_WP
          CALL ecri_psi(psi=WP(i),                               &
                        ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)
        END DO
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------
      END SUBROUTINE march_gene
!================================================================
!
!    march RK2
!         a time dependant pulse in Hamiltonian (W(t))
!         H is the square matrix (dimension n)
!
!================================================================
      !!@description: march RK2  a time dependant pulse
      !!             in Hamiltonian (W(t)) H is the square matrix (dimension n)
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      SUBROUTINE march_RK2_field(T,WP,nb_WP,print_Op,                   &
                                 para_H,para_Dip,para_field,para_propa)
      USE mod_system
      USE mod_Op,              ONLY : param_Op
      USE mod_psi_set_alloc,   ONLY : param_psi,ecri_psi
      USE mod_ana_psi,         ONLY : norm2_psi
      USE mod_psi_SimpleOp
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H
      TYPE (param_Op)   :: para_Dip(3)

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_field) :: para_field

      integer          :: nb_WP
      TYPE (param_psi) :: WP(:)
      TYPE (param_psi), pointer :: w1,w2,w3,w4,w5,w6

!----- for printing --------------------------------------------------
      logical :: print_Op




!------ working parameters --------------------------------


      integer :: i
      integer           :: ip
      real (kind=Rkind) :: T      ! time
      real (kind=Rkind) :: T_DT
      real (kind=Rkind) :: phase

!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='march_RK2_field'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Tmax,deltaT',para_propa%WPTmax,para_propa%WPdeltaT
        write(out_unitp,*) 'Hmin,Hmax',para_propa%para_poly%Hmin,               &
                                para_propa%para_poly%Hmax
        write(out_unitp,*)
        write(out_unitp,*) 'nb_ba,nb_qa',WP%nb_ba,WP%nb_qa
        write(out_unitp,*) 'nb_bi',WP%nb_bi
        write(out_unitp,*)

        DO i=1,nb_WP
          CALL norm2_psi(WP(i),.FALSE.,.TRUE.,.FALSE.)
          write(out_unitp,*) 'norm WP',i,WP(i)%norme

          write(out_unitp,*) 'WP BasisRep',i
          CALL ecri_psi(T=ZERO,psi=WP(i),                              &
                        ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)
          write(out_unitp,*) 'WP GridRep',i
          CALL ecri_psi(T=ZERO,psi=WP(i),                              &
                        ecri_GridRep=.TRUE.,ecri_BasisRep=.FALSE.)
        END DO
       END IF
!-----------------------------------------------------------
      w1 => para_propa%work_WP(1)
      w2 => para_propa%work_WP(2)
      w3 => para_propa%work_WP(3)
      w4 => para_propa%work_WP(4)
      w5 => para_propa%work_WP(5)
      w6 => para_propa%work_WP(6)

      DO i=1,nb_WP
!-----------------------------------------------------------
!      w1 = WPdeltaT * (-i)(H - dnE(T).Dip).WP
!      w1 = WPdeltaT * fcn(WP,T)

       CALL fcn_field(T,WP(i),w1,para_H,para_Dip,para_field,w5)
       w1 = w1 * para_propa%WPdeltaT

!-----------------------------------------------------------
!      w6 = WP + w1
!      w2 = WPdeltaT * (-i)(H - dnE(T+DT).Dip).w6
!      w2 = WPdeltaT * fcn(WP+w1,T+DT)
       T_DT = T + para_propa%WPdeltaT
       w6 = WP(i) + w1

       CALL fcn_field(T_DT,w6,w2,para_H,para_Dip,para_field,w5)
       w2 = w2 * para_propa%WPdeltaT

       w5 = w1+w2
       WP(i) = WP(i) + w5*HALF

!      - Phase Shift -----------------
       phase = para_H%E0*para_propa%WPdeltaT
       WP(i) = WP(i) * exp(-cmplx(ZERO,phase,kind=Rkind))

!     - check norm ------------------
      CALL norm2_psi(WP(i),GridRep=.FALSE.,BasisRep=.TRUE.)
      IF ( WP(i)%norme > WP(i)%max_norme) THEN
        T  = T + para_propa%WPdeltaT
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' STOP propagation: norme > max_norme',              &
                         WP%norme
        para_propa%march_error   = .TRUE.
        para_propa%test_max_norm = .TRUE.
        STOP
      END IF

      END DO

!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------


      END SUBROUTINE march_RK2_field
!================================================================
!
!    march RK4
!         a time dependant pulse in Hamiltonian (W(t))
!         H is the square matrix (dimension n)
!
!================================================================
      !!@description: ch RK4
      !!         a time dependant pulse in Hamiltonian (W(t))
      !!         H is the square matrix (dimension n)
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      SUBROUTINE march_RK4_field(T,WP,nb_WP,print_Op,                   &
                                 para_H,para_Dip,para_field,para_propa)
      USE mod_system
      USE mod_Op,              ONLY : param_Op
      USE mod_psi_set_alloc,   ONLY : param_psi,ecri_psi
      USE mod_ana_psi,          ONLY : norm2_psi
      USE mod_psi_SimpleOp
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H
      TYPE (param_Op)   :: para_Dip(3)

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_field) :: para_field

      integer          :: nb_WP
      TYPE (param_psi) :: WP(:)
      TYPE (param_psi), pointer :: w1,w2,w3,w4,w5,w6

!----- for printing --------------------------------------------------
      logical ::print_Op



!------ working parameters --------------------------------


      integer           :: i,ip
      real (kind=Rkind) :: T      ! time
      real (kind=Rkind) :: T_DT,T_DT_HALF
      real (kind=Rkind) :: DT,DTo2,DTo6
      real (kind=Rkind) :: phase

      integer  ::   nioWP

!----- for the field --------------------------------------------------
      real (kind=Rkind)    :: dnE(3)
      real (kind=Rkind)    :: ww(3)
      logical :: make_field
!----- for the field --------------------------------------------------


!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='march_RK4_field'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Tmax,deltaT',para_propa%WPTmax,para_propa%WPdeltaT
        write(out_unitp,*) 'Hmin,Hmax',para_propa%para_poly%Hmin,               &
                                para_propa%para_poly%Hmax
        write(out_unitp,*)
        write(out_unitp,*) 'nb_ba,nb_qa',WP%nb_ba,WP%nb_qa
        write(out_unitp,*) 'nb_bi',WP%nb_bi
        write(out_unitp,*)

        DO i=1,nb_WP
          CALL norm2_psi(WP(i),.FALSE.,.TRUE.,.FALSE.)
          write(out_unitp,*) 'norm WP',i,WP(i)%norme

          write(out_unitp,*) 'WP BasisRep',i
          CALL ecri_psi(T=ZERO,psi=WP(i),                               &
                        ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)
          write(out_unitp,*) 'WP GridRep',i
          CALL ecri_psi(T=ZERO,psi=WP(i),                               &
                        ecri_GridRep=.TRUE.,ecri_BasisRep=.FALSE.)
        END DO

      END IF
!-----------------------------------------------------------
      T_DT      = T + para_propa%WPdeltaT
      T_DT_HALF = T + para_propa%WPdeltaT*HALF
      DT        = para_propa%WPdeltaT
      DTo2      = para_propa%WPdeltaT/TWO
      DTo6      = para_propa%WPdeltaT/SIX

      DO i=1,6
        para_propa%work_WP(i) = WP(1)
      END DO
      w1 => para_propa%work_WP(1)
      w2 => para_propa%work_WP(2)
      w3 => para_propa%work_WP(3)
      w4 => para_propa%work_WP(4)
      w5 => para_propa%work_WP(5)
      w6 => para_propa%work_WP(6)
!-----------------------------------------------------------

      DO i=1,nb_WP

!     w1 = -i*H(T)*WP
      CALL fcn_field(T,WP(i),w1,para_H,para_Dip,para_field,w5)

!     w2 = WP + w1 * (DT/2)
      w2 = WP(i) + w1 * DTo2

!     w3 = -iH(T+DT/2)*w2
      IF (para_field%type == 'grid') THEN
        CALL fcn_field(T,w2,w3,para_H,para_Dip,para_field,w5)
      ELSE
        CALL fcn_field(T_DT_HALF,w2,w3,para_H,para_Dip,para_field,w5)
      END IF

!     w2 = WP + w3 * (DT/2)
      w2 = WP(i) + w3 * DTo2

!     w4 = -iH(T+DT/2)*w2
      IF (para_field%type == 'grid') THEN
        CALL fcn_field(T,w2,w4,para_H,para_Dip,para_field,w5)
      ELSE
        CALL fcn_field(T_DT_HALF,w2,w4,para_H,para_Dip,para_field,w5)
      END IF

!     w2 = WP + w4 * DT
      w2 = WP(i) + w4 * DT

!     w4 = w3 + w4
      w4 = w3 + w4

!     w3 = -iH(T)*w2
      IF (para_field%type == 'grid') THEN
        CALL fcn_field(T,w2,w3,para_H,para_Dip,para_field,w5)
      ELSE
        CALL fcn_field(T_DT,w2,w3,para_H,para_Dip,para_field,w5)
      END IF

!     WP = WP + (DT/6)*(w1+w3+2.D0*w4)
      w5 = w1 + w3
      w5 = w5 + w4*TWO
      WP(i) = WP(i) + w5*DTo6


!     - Phase Shift -----------------
      phase = para_H%E0*para_propa%WPdeltaT
      WP(i) = WP(i) * exp(-cmplx(ZERO,phase,kind=Rkind))

!     - check norm ------------------
      CALL norm2_psi(WP(i),GridRep=.FALSE.,BasisRep=.TRUE.)
      IF ( WP(i)%norme > WP(i)%max_norme) THEN
        T  = T + para_propa%WPdeltaT
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' STOP propagation: norme > max_norme',              &
                         WP%norme
        para_propa%march_error   = .TRUE.
        para_propa%test_max_norm = .TRUE.
        STOP
      END IF

      END DO

!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------

      END SUBROUTINE march_RK4_field
      SUBROUTINE march_RK4(T,no,WP,WP0,para_H,para_propa)
      USE mod_system
      USE mod_Op,              ONLY : param_Op
      USE mod_psi_set_alloc,   ONLY : param_psi,ecri_psi
      USE mod_ana_psi,          ONLY : norm2_psi
      USE mod_psi_SimpleOp
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa

      TYPE (param_psi) :: WP,WP0
      TYPE (param_psi), pointer :: w1,w2,w3,w4,w5,w6

      integer :: no
!----- for printing --------------------------------------------------
      logical ::print_Op



!------ working parameters --------------------------------

      complex (kind=Rkind) :: cdot

      integer           :: i,ip
      real (kind=Rkind) :: T      ! time
      real (kind=Rkind) :: T_DT,T_DT_HALF
      real (kind=Rkind) :: DT,DTo2,DTo6
      real (kind=Rkind) :: phase

      integer  ::   nioWP

!----- for the field --------------------------------------------------
      real (kind=Rkind)    :: dnE(3)
      real (kind=Rkind)    :: ww(3)
      logical :: make_field
!----- for the field --------------------------------------------------


!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='march_RK4'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Tmax,deltaT',para_propa%WPTmax,para_propa%WPdeltaT
        write(out_unitp,*) 'Hmin,Hmax',para_propa%para_poly%Hmin,               &
                                para_propa%para_poly%Hmax
        write(out_unitp,*)
        write(out_unitp,*) 'nb_ba,nb_qa',WP%nb_ba,WP%nb_qa
        write(out_unitp,*) 'nb_bi',WP%nb_bi
        write(out_unitp,*)

          CALL norm2_psi(WP,.FALSE.,.TRUE.,.FALSE.)
          write(out_unitp,*) 'norm WP',i,WP%norme

          write(out_unitp,*) 'WP BasisRep'
          CALL ecri_psi(T=ZERO,psi=WP,                               &
                        ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)
          write(out_unitp,*) 'WP GridRep'
          CALL ecri_psi(T=ZERO,psi=WP,                               &
                        ecri_GridRep=.TRUE.,ecri_BasisRep=.FALSE.)

      END IF
!-----------------------------------------------------------
      T_DT      = T + para_propa%WPdeltaT
      T_DT_HALF = T + para_propa%WPdeltaT*HALF
      DT        = para_propa%WPdeltaT
      DTo2      = para_propa%WPdeltaT/TWO
      DTo6      = para_propa%WPdeltaT/SIX

      DO i=1,6
        para_propa%work_WP = WP
      END DO
      w1 => para_propa%work_WP(1)
      w2 => para_propa%work_WP(2)
      w3 => para_propa%work_WP(3)
      w4 => para_propa%work_WP(4)
      w5 => para_propa%work_WP(5)
      w6 => para_propa%work_WP(6)
!-----------------------------------------------------------

      !w1 = -i.H.WP
      CALL fcn(WP,w1,para_H)

      !w2 = WP + w1 * (DT/2)
      w2 = WP + w1 * DTo2

      !w3 = -iH*w2
      CALL fcn(w2,w3,para_H)


      !w2 = WP + w3 * (DT/2)
      w2 = WP + w3 * DTo2

      !w4 = -iH*w2
      CALL fcn(w2,w4,para_H)

      !w2 = WP + w4 * DT
      w2 = WP + w4 * DT

      !w4 = w3 + w4
      w4 = w3 + w4

      !w3 = -iH(T)*w2
      CALL fcn(w2,w3,para_H)


      !WP = WP + (DT/6)*(w1+w3+2.D0*w4)
      w5 = w1 + w3
      w5 = w5 + w4*TWO
      WP = WP + w5*DTo6


      !- Phase Shift -----------------
      phase = para_H%E0*para_propa%WPdeltaT
      WP = WP * exp(-cmplx(ZERO,phase,kind=Rkind))

      !- check norm ------------------
      CALL norm2_psi(WP,GridRep=.FALSE.,BasisRep=.TRUE.)
      IF ( WP%norme > WP%max_norme) THEN
        T  = T + para_propa%WPdeltaT
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' STOP propagation: norme > max_norme',              &
                         WP%norme
        para_propa%march_error   = .TRUE.
        para_propa%test_max_norm = .TRUE.
        STOP
      END IF

      cdot = Calc_AutoCorr(WP0,WP,para_propa,T,Write_AC=.FALSE.)
      CALL Write_AutoCorr(no,T + para_propa%WPdeltaT,cdot)
      CALL flush_perso(no)


!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------

      END SUBROUTINE march_RK4
      SUBROUTINE march_BS(T,no,WP,WP0,para_H,para_propa)
      USE mod_system
      USE mod_Op,              ONLY : param_Op
      USE mod_psi_set_alloc,   ONLY : param_psi,ecri_psi,dealloc_psi
      USE mod_ana_psi,          ONLY : norm2_psi
      USE mod_psi_SimpleOp
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa

      TYPE (param_psi) :: WP,WP0



      TYPE (param_psi), allocatable :: yt0(:),yt1(:)
      TYPE (param_psi) :: yerr,ycorr



      integer :: no
!----- for printing --------------------------------------------------
      logical ::print_Op



!------ working parameters --------------------------------

      complex (kind=Rkind) :: cdot
      real(kind=Rkind)     :: x,err0,err1
      integer              :: i,j,n,m
      integer              :: ip, order
      real (kind=Rkind)    :: T      ! time
      real (kind=Rkind)    :: DT
      real (kind=Rkind)    :: phase

      integer  ::   nioWP




!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='march_BS'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Tmax,deltaT',para_propa%WPTmax,para_propa%WPdeltaT
        write(out_unitp,*) 'Hmin,Hmax',para_propa%para_poly%Hmin,               &
                                para_propa%para_poly%Hmax
        write(out_unitp,*)
        write(out_unitp,*) 'nb_ba,nb_qa',WP%nb_ba,WP%nb_qa
        write(out_unitp,*) 'nb_bi',WP%nb_bi
        write(out_unitp,*)

          CALL norm2_psi(WP,.FALSE.,.TRUE.,.FALSE.)
          write(out_unitp,*) 'norm WP',i,WP%norme

          write(out_unitp,*) 'WP BasisRep'
          CALL ecri_psi(T=ZERO,psi=WP,                               &
                        ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)
          write(out_unitp,*) 'WP GridRep'
          CALL ecri_psi(T=ZERO,psi=WP,                               &
                        ecri_GridRep=.TRUE.,ecri_BasisRep=.FALSE.)

      END IF
!-----------------------------------------------------------

  order = para_propa%para_poly%npoly
  IF (order < 1) THEN
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) ' order < 1  ...'
    write(out_unitp,*) ' npoly has a wrong value',para_propa%para_poly%npoly
  END IF

  DT        = para_propa%WPdeltaT

  yerr = WP
  yerr = ZERO

  allocate(yt0(0:0))
  m = 2
  yt0(0) = WP

  CALL march_ModMidPoint(T,no,yt0(0),WP0,para_H,para_propa,m)
  err0 = huge(ONE)

  DO j=1,order
    allocate(yt1(0:j))
    m = 2*j+2
    yt1(0) = WP
    CALL march_ModMidPoint(T,no,yt1(0),WP0,para_H,para_propa,m)

    !extrapolation
    DO i=1,j
      x = real(2*j+2,kind=8)/real(2*(j-i)+2,kind=8)
      x = ONE/(x**2-ONE)
      yerr = yt1(i-1)-yt0(i-1)
      yt1(i) = yt1(i-1) + yerr*x
    END DO

    yerr = yt1(j) - yt0(j-1)
    CALL norm2_psi(yerr)
    err1 = sqrt(yerr%norme)

    DO i=lbound(yt0,dim=1),ubound(yt0,dim=1)
      CALL dealloc_psi(yt0(i),delete_all=.TRUE.)
    END DO
    deallocate(yt0)
    allocate(yt0(0:j))
    DO i=lbound(yt1,dim=1),ubound(yt1,dim=1)
      yt0(i) = yt1(i)
      CALL dealloc_psi(yt1(i),delete_all=.TRUE.)
    END DO
    deallocate(yt1)
    !write(6,*) 'march_bs',j,err1

    !IF (err1 < 1.d-6 .OR. err1 > err0) EXIT
    IF (err1 < 1.d-10) EXIT

    err0 = err1
  END DO
  write(6,*) 'end march_bs',min(j,order),err1

  WP = yt0(min(j,order))
  DO i=lbound(yt0,dim=1),ubound(yt0,dim=1)
    CALL dealloc_psi(yt0(i),delete_all=.TRUE.)
  END DO
  deallocate(yt0)
  CALL dealloc_psi(yerr,delete_all=.TRUE.)


      !- Phase Shift -----------------
      phase = para_H%E0*para_propa%WPdeltaT
      WP = WP * exp(-cmplx(ZERO,phase,kind=Rkind))

      !- check norm ------------------
      CALL norm2_psi(WP,GridRep=.FALSE.,BasisRep=.TRUE.)
      IF ( WP%norme > WP%max_norme) THEN
        T  = T + para_propa%WPdeltaT
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' STOP propagation: norme > max_norme',              &
                         WP%norme
        para_propa%march_error   = .TRUE.
        para_propa%test_max_norm = .TRUE.
        STOP
      END IF

      cdot = Calc_AutoCorr(WP0,WP,para_propa,T,Write_AC=.FALSE.)
      CALL Write_AutoCorr(no,T + para_propa%WPdeltaT,cdot)
      CALL flush_perso(no)


!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------

      END SUBROUTINE march_BS
      SUBROUTINE march_ModMidPoint(T,no,WP,WP0,para_H,para_propa,order)
      USE mod_system
      USE mod_Op,              ONLY : param_Op
      USE mod_psi_set_alloc,   ONLY : param_psi,ecri_psi,dealloc_psi
      USE mod_ana_psi,          ONLY : norm2_psi
      USE mod_psi_SimpleOp
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      integer, optional  :: order
      TYPE (param_psi) :: WP,WP0


      TYPE (param_psi) :: dWP,zkm,zk,zkp

      integer :: no
!----- for printing --------------------------------------------------
      logical ::print_Op



!------ working parameters --------------------------------

      complex (kind=Rkind) :: cdot

      integer           :: i,k,ip,order_loc
      real (kind=Rkind) :: T      ! time
      real (kind=Rkind) :: DTT
      real (kind=Rkind) :: phase

      integer  ::   nioWP

!----- for the field --------------------------------------------------
      real (kind=Rkind)    :: dnE(3)
      real (kind=Rkind)    :: ww(3)
      logical :: make_field
!----- for the field --------------------------------------------------


!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='march_ModMidPoint'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Tmax,deltaT',para_propa%WPTmax,para_propa%WPdeltaT
        write(out_unitp,*) 'Hmin,Hmax',para_propa%para_poly%Hmin,               &
                                para_propa%para_poly%Hmax
        write(out_unitp,*)
        write(out_unitp,*) 'nb_ba,nb_qa',WP%nb_ba,WP%nb_qa
        write(out_unitp,*) 'nb_bi',WP%nb_bi
        write(out_unitp,*)

          CALL norm2_psi(WP,.FALSE.,.TRUE.,.FALSE.)
          write(out_unitp,*) 'norm WP',i,WP%norme

          write(out_unitp,*) 'WP BasisRep'
          CALL ecri_psi(T=ZERO,psi=WP,                               &
                        ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)
          write(out_unitp,*) 'WP GridRep'
          CALL ecri_psi(T=ZERO,psi=WP,                               &
                        ecri_GridRep=.TRUE.,ecri_BasisRep=.FALSE.)

      END IF
!-----------------------------------------------------------
      IF (present(order)) THEN
        order_loc = order
      ELSE
        order_loc = para_propa%para_poly%npoly
      END IF
      IF (order_loc < 1) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' order_loc < 1  ...'
        write(out_unitp,*) ' npoly or order have wrong values'
        write(out_unitp,*) ' npoly: ',para_propa%para_poly%npoly
        IF (present(order))  &
          write(out_unitp,*) ' order: ',order
        STOP
      END IF
!write(6,*) 'order_loc',order_loc
      DTT        = para_propa%WPdeltaT / real(order_loc,kind=Rkind)

      zkm = WP
      CALL fcn(WP,dWP,para_H)
      zk  = WP + dWP * DTT


      DO k=1,order_loc-1

        CALL fcn(zk,dWP,para_H)
        zkp = zkm + TWO*DTT * dWP
        zkm = zk
        zk  = zkp

      END DO

      CALL fcn(zk,zkp,para_H)
      dWP = DTT * zkp

      zkp = zk + zkm

      zk =  zkp + dWP
      WP = HALF * zk


      CALL dealloc_psi(dWP,delete_all=.TRUE.)
      CALL dealloc_psi(zkm,delete_all=.TRUE.)
      CALL dealloc_psi(zk ,delete_all=.TRUE.)
      CALL dealloc_psi(zkp,delete_all=.TRUE.)

      IF (present(order)) RETURN

      !- Phase Shift -----------------
      phase = para_H%E0*para_propa%WPdeltaT
      WP = WP * exp(-cmplx(ZERO,phase,kind=Rkind))

      !- check norm ------------------
      CALL norm2_psi(WP,GridRep=.FALSE.,BasisRep=.TRUE.)
      IF ( WP%norme > WP%max_norme) THEN
        T  = T + para_propa%WPdeltaT
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' STOP propagation: norme > max_norme',              &
                         WP%norme
        para_propa%march_error   = .TRUE.
        para_propa%test_max_norm = .TRUE.
        STOP
      END IF

      cdot = Calc_AutoCorr(WP0,WP,para_propa,T,Write_AC=.FALSE.)
      CALL Write_AutoCorr(no,T + para_propa%WPdeltaT,cdot)
      CALL flush_perso(no)


!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------

      END SUBROUTINE march_ModMidPoint

!================================================================
!
!   fcn =>  d psi /dt = fcn(psi,t)
!      fcn : -i (H -mu.e(t))|psi>
!
!================================================================
      !!@description: n =>  d psi /dt = fcn(psi,t)
      !!      fcn : -i (H -mu.e(t))|psi>
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      SUBROUTINE fcn_field(T,WP,dWP,para_H,para_Dip,                    &
                           para_field,w1)
      USE mod_system

      USE mod_psi_set_alloc,   ONLY : param_psi,ecri_psi
      USE mod_ana_psi,          ONLY : norm2_psi
      USE mod_psi_SimpleOp

      USE mod_field,           ONLY : param_field,sub_dnE
      USE mod_Op,              ONLY : param_Op, sub_OpPsi,sub_scaledOpPsi
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H
      TYPE (param_Op)   :: para_Dip(3)

!----- variables for the WP propagation ----------------------------
      TYPE (param_field) :: para_field

      TYPE (param_psi) :: WP,dWP,w1

!----- for printing --------------------------------------------------
      logical ::print_Op



!------ working parameters --------------------------------


      integer           :: ip
      real (kind=Rkind) :: T      ! time
      real (kind=Rkind) :: T_DT

      integer  ::   nioWP

!----- for the field --------------------------------------------------
      real (kind=Rkind)    :: dnE(3)
      real (kind=Rkind)    :: ww(3)
      logical :: make_field
!----- for the field --------------------------------------------------


!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='fcn_field'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_ba,nb_qa',WP%nb_ba,WP%nb_qa
        write(out_unitp,*) 'nb_bi',WP%nb_bi
        write(out_unitp,*)

        CALL norm2_psi(WP,.FALSE.,.TRUE.,.FALSE.)
        write(out_unitp,*) 'norm WP BasisRep',WP%norme

        write(out_unitp,*) 'WP'
        CALL ecri_psi(T=T,psi=WP)
       END IF
!-----------------------------------------------------------

!-----------------------------------------------------------
!      dWP = (-i)(H - dnE(T).Dip).WP

       CALL sub_OpPsi(WP,dWP,para_H)
       CALL sub_scaledOpPsi(WP,dWP,para_H%E0,ONE)
       CALL sub_dnE(dnE,0,T,para_field)
       DO ip=1,3
         IF (.NOT. para_field%pola_xyz(ip)) CYCLE

         w1 = ZERO
         CALL sub_OpPsi(WP,w1,para_Dip(ip))
         dWP = dWP - w1 * dnE(ip)

       END DO
       dWP = dWP * (-EYE)

      IF (debug) THEN
        write(out_unitp,*) 'dWP'
        CALL ecri_psi(T=T,psi=dWP)
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------


      END SUBROUTINE fcn_field
      SUBROUTINE fcn(WP,dWP,para_H)
      USE mod_system

      USE mod_psi_set_alloc,   ONLY : param_psi,ecri_psi
      USE mod_ana_psi,          ONLY : norm2_psi
      USE mod_psi_SimpleOp
      USE mod_Op,              ONLY : param_Op, sub_OpPsi,sub_scaledOpPsi
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ------------------------------
      TYPE (param_psi) :: WP,dWP

!----- for printing --------------------------------------------------
      logical ::print_Op


      integer  ::   nioWP


!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='fcn'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_ba,nb_qa',WP%nb_ba,WP%nb_qa
        write(out_unitp,*) 'nb_bi',WP%nb_bi
        write(out_unitp,*)

        CALL norm2_psi(WP,.FALSE.,.TRUE.,.FALSE.)
        write(out_unitp,*) 'norm WP BasisRep',WP%norme

        write(out_unitp,*) 'WP'
        CALL ecri_psi(T=ZERO,psi=WP)
       END IF
!-----------------------------------------------------------

!-----------------------------------------------------------
!      dWP = -i.H.WP

       CALL sub_OpPsi(WP,dWP,para_H)
       CALL sub_scaledOpPsi(WP,dWP,para_H%E0,ONE)
       dWP = dWP * (-EYE)

      IF (debug) THEN
        write(out_unitp,*) 'dWP'
        CALL ecri_psi(T=ZERO,psi=dWP)
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------


      END SUBROUTINE fcn
!================================================================
!
!    march nOD propagation with
!         a time dependant pulse in Hamiltonian (W(t))
!
!================================================================
      !!@description: arch nOD propagation with
      !!         a time dependant pulse in Hamiltonian (W(t))
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      SUBROUTINE march_noD_field(T,WP,nb_WP,print_Op,                   &
                                 para_H,para_Dip,                       &
                                 para_field,para_propa)
      USE mod_system

      USE mod_field,           ONLY : param_field,sub_dnE
      USE mod_Op,              ONLY : param_Op, sub_OpPsi,sub_scaledOpPsi

      USE mod_psi_set_alloc,   ONLY : param_psi,ecri_psi
      USE mod_ana_psi,          ONLY : norm2_psi
      USE mod_psi_SimpleOp
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H
      TYPE (param_Op)   :: para_Dip(3)

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_field) :: para_field

      integer          :: nb_WP
      TYPE (param_psi) :: WP(:)
      TYPE (param_psi), pointer :: w1,w2

!----- for printing --------------------------------------------------
      logical ::print_Op



!------ working parameters --------------------------------

      complex (kind=Rkind)   :: Overlap,cdot,rt,rti,fac_k,norm
      real (kind=Rkind) :: phase,limit

      integer       :: it,it_max,no,max_der,nb_der
      integer       :: i,j,k,jt,ip,iq
      real (kind=Rkind) :: T      ! time
      real (kind=Rkind) :: T_Delta! time+deltaT
      integer       :: max_ecri

      integer  ::   nioWP

!----- for the field --------------------------------------------------
      real (kind=Rkind), pointer :: tab_dnE(:,:)
      real (kind=Rkind)    :: dnE(3)
      real (kind=Rkind)    :: ww(3)
      logical :: make_field
!----- for the field --------------------------------------------------


!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug=.FALSE.
!      logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='march_nOd_field'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Tmax,deltaT',para_propa%WPTmax,para_propa%WPdeltaT
        write(out_unitp,*) 'Hmin,Hmax',para_propa%para_poly%Hmin,               &
                                para_propa%para_poly%Hmax
        write(out_unitp,*)
        write(out_unitp,*) 'nb_ba,nb_qa',WP(1)%nb_ba,WP(1)%nb_qa
        write(out_unitp,*) 'nb_bi',WP(1)%nb_bi
        write(out_unitp,*)

        DO i=1,nb_WP
          CALL norm2_psi(WP(i),.FALSE.,.TRUE.,.FALSE.)
          write(out_unitp,*) 'norm WP',i,WP(i)%norme

          write(out_unitp,*) 'WP BasisRep',i
          CALL ecri_psi(T=ZERO,psi=WP(i),                               &
                        ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)
          write(out_unitp,*) 'WP GridRep',i
          CALL ecri_psi(T=ZERO,psi=WP(i),                               &
                        ecri_GridRep=.TRUE.,ecri_BasisRep=.FALSE.)
        END DO
      END IF
!-----------------------------------------------------------


!-----------------------------------------------------------
      w1 => para_propa%work_WP(para_propa%para_poly%npoly+1)
      w2 => para_propa%work_WP(para_propa%para_poly%npoly+2)

      max_der = 0
      nb_der = para_propa%para_poly%npoly-1
      IF (para_field%max_der >= 0)                                      &
           nb_der = min(para_field%max_der,para_propa%para_poly%npoly-1)

      nullify(tab_dnE)
      CALL alloc_array(tab_dnE,(/nb_der,3/),"tab_dnE","march_noD_field",(/0,1/))
      tab_dnE(:,:) = ZERO
!-----------------------------------------------------------

!-----------------------------------------------------------

!     - propagation ----------------
      fac_k = cmplx(para_propa%WPdeltaT,kind=Rkind)
      CALL sub_dnE(dnE,0,T,para_field)
      tab_dnE(0,:) = -dnE(:)
      DO k=1,nb_der
        CALL sub_dnE(dnE,k,T,para_field)
        tab_dnE(k,:) = -dnE(:)
        fac_k = fac_k * cmplx(para_propa%WPdeltaT/real(k,kind=Rkind),kind=Rkind)
        limit = maxval(abs(dnE(:)))*fac_k
        IF (limit > para_propa%para_poly%poly_tol) max_der=k
!       write(out_unitp,*) 'k,max_der,limit',k,max_der,limit
      END DO
!     max_der=para_propa%para_poly%npoly-1
!     write(out_unitp,*) 'max_der',max_der

      DO j=1,nb_WP
        para_propa%work_WP(0) = WP(j)
        DO i=0,para_propa%para_poly%npoly-1

          CALL sub_OpPsi(para_propa%work_WP(i),w2,para_H)
          CALL sub_scaledOpPsi(para_propa%work_WP(i),w2,para_H%E0,ONE)

          para_propa%work_WP(i+1) = w2  ! H.psi(i)


          fac_k = cmplx(ZERO,-para_propa%WPdeltaT /                     &
                           real(i+1,kind=Rkind),kind=Rkind)

          para_propa%work_WP(i+1) = para_propa%work_WP(i+1) * fac_k

          DO ip=1,3
            IF (.NOT. para_field%pola_xyz(ip)) CYCLE
            IF (maxval(abs(tab_dnE(:,ip))) <                            &
                para_propa%para_poly%poly_tol) CYCLE
            fac_k = cmplx(ZERO,-para_propa%WPdeltaT /                   &
                                 real(i+1,kind=Rkind),kind=Rkind)

            w1 = para_propa%work_WP(i) * (cmplx(tab_dnE(0,ip),kind=Rkind) * fac_k)

            DO k=1,min(i,max_der)
               fac_k = fac_k * cmplx(para_propa%WPdeltaT/real(k,kind=Rkind),kind=Rkind)
               rt = cmplx(tab_dnE(k,ip),kind=Rkind)*fac_k
               w1 = w1 + para_propa%work_WP(i-k) * rt
            END DO

            CALL sub_OpPsi(w1,w2,para_Dip(ip))

            para_propa%work_WP(i+1) = para_propa%work_WP(i+1) + w2
          END DO
          WP(j) = WP(j) + para_propa%work_WP(i+1)

          CALL norm2_psi(para_propa%work_WP(i+1),.FALSE.,.TRUE.,.FALSE.)
!         write(out_unitp,*) 'norme psi i+1',i+1,para_propa%work_WP(i+1)%norme
          IF (para_propa%work_WP(i+1)%norme <                           &
                        para_propa%para_poly%poly_tol) EXIT
          IF (para_propa%work_WP(i+1)%norme > TEN**15) THEN
             write(out_unitp,*) ' ERROR in ',name_sub
             write(out_unitp,*) ' The norme of dnpsi(i+1) is huge !',           &
                                          para_propa%work_WP(i+1)%norme
             write(out_unitp,*) ' iteration i:',i
             write(out_unitp,*) ' Reduce the time step, WPDeltaT:',             &
                             para_propa%WPdeltaT
             para_propa%march_error = .TRUE.
             STOP
          END IF
        END DO



        IF (print_Op) write(out_unitp,*) 'T,max_der,ipoly,norme',               &
                                  T,max_der,i,                          &
        para_propa%work_WP(min(i+1,para_propa%para_poly%npoly))%norme

        para_propa%work_WP(0) =                                         &
        para_propa%work_WP(min(i+1,para_propa%para_poly%npoly))

!       - Phase Shift -----------------
        phase = para_H%E0*para_propa%WPdeltaT
        WP(j) = WP(j) * exp(-cmplx(ZERO,phase,kind=Rkind))

!       - check norm ------------------
        CALL norm2_psi(WP(j),GridRep=.FALSE.,BasisRep=.TRUE.)
        IF ( WP(j)%norme > WP(j)%max_norme) THEN
          T  = T + para_propa%WPdeltaT
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' STOP propagation: norm > max_norm',              &
                         WP(j)%norme
          para_propa%march_error   = .TRUE.
          para_propa%test_max_norm = .TRUE.
          STOP
        END IF

      END DO

      CALL dealloc_array(tab_dnE,"tab_dnE","march_noD_field")


!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------


      END SUBROUTINE march_noD_field
!================================================================
!
!     march Taylor
!
!================================================================
      !!@description: march Taylor
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      SUBROUTINE march_noD(T,no,psi,psi0,para_H,para_propa)
      USE mod_system
      USE mod_Op,              ONLY : param_Op, sub_PsiOpPsi, sub_OpPsi,sub_scaledOpPsi

      USE mod_psi_set_alloc,   ONLY : param_psi,ecri_psi
      USE mod_ana_psi,          ONLY : norm2_psi
      USE mod_psi_SimpleOp
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_psi)   :: psi0,psi

!------ working variables ---------------------------------
      TYPE (param_psi), pointer :: w1,w2
      complex (kind=Rkind) :: E,cdot
      complex (kind=Rkind) :: psi0Hkpsi0(0:para_propa%para_poly%npoly)
      real (kind=Rkind)    :: microT,T,microdeltaT,phase,microphase
      integer              :: no

      complex (kind=Rkind) :: rt,rtj,rt_tmp
      integer              :: it,i,j,i_qaie_corr,j_exit
      real (kind=Rkind)    :: rg

!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='march_nOd'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'deltaT',para_propa%WPdeltaT
        write(out_unitp,*) 'E0,Esc',para_H%E0,para_H%Esc
        write(out_unitp,*) 'nOD',para_propa%para_poly%npoly
      END IF
!-----------------------------------------------------------

      w1 => para_propa%work_WP(1)
      w2 => para_propa%work_WP(2)

      w2 = psi
      CALL sub_PsiOpPsi(E,psi,w2,para_H)
!     para_H%E0 = real(E,kind=Rkind)
!     write(out_unitp,*) 'E0,Esc',para_H%E0,para_H%Esc

      psi0Hkpsi0(:) = cmplx(ZERO,ZERO,kind=Rkind)

      psi0Hkpsi0(0) = Calc_AutoCorr(psi0,psi,para_propa,T,Write_AC=.FALSE.)

      w1 = psi
      rtj = cmplx(ONE,ZERO,kind=Rkind)

 21 format(a,100(x,e12.5))

      !write(6,21) 'Rw1',Real(w1%CvecB,kind=Rkind)
      !write(6,21) 'Iw1',AImag(w1%CvecB)

      DO j=1,para_propa%para_poly%npoly

        !write(6,21) 'Rw1',Real(w1%CvecB,kind=Rkind)
        !write(6,21) 'Iw1',AImag(w1%CvecB)

        IF (j > 1) CALL sub_OpPsi(w1,w2,para_H) ! already done in PsiHPsi
        !write(6,21) 'Rw2',Real(w2%CvecB,kind=Rkind)
        !write(6,21) 'Iw2',AImag(w2%CvecB)

        CALL sub_scaledOpPsi(w1,w2,para_H%E0,ONE)

        !write(6,21) 'Rw2',Real(w2%CvecB,kind=Rkind)
        !write(6,21) 'Iw2',AImag(w2%CvecB)

        w1 = w2

        psi0Hkpsi0(j) = Calc_AutoCorr(psi0,w1,para_propa,T,Write_AC=.FALSE.)

        !CALL Overlap_psi1_psi2(psi0Hkpsi0(j),psi0,w1)

        rtj = rtj *                                                     &
          cmplx(ZERO,-para_propa%WPdeltaT/real(j,kind=Rkind),kind=Rkind)
        w2 = w1 * rtj

        !write(6,21) 'Rw2*rtj',Real(w2%CvecB,kind=Rkind)
        !write(6,21) 'Iw2*rtj',AImag(w2%CvecB)

        psi = psi + w2

        !write(6,21) 'Rpsi',Real(psi%CvecB,kind=Rkind)
        !write(6,21) 'Ipsi',AImag(psi%CvecB)

        CALL norm2_psi(w2)

        IF (debug) write(out_unitp,*) 'j,norme w2',j,w2%norme

        IF (w2%norme > TEN**15) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' Norm of the vector is TOO large (> 10^15)',j,    &
                                                w2%norme
          write(out_unitp,*) ' => Reduce the time step !!'
          STOP
        END IF
        j_exit = j
        IF (w2%norme < para_propa%para_poly%poly_tol) EXIT

      END DO
      write(out_unitp,*) 'j_exit,norms',j_exit,abs(w2%norme)

!     write(out_unitp,*) 'j,norme',j,abs(w1%norme*rtj)
      IF (abs(w2%norme) > para_propa%para_poly%poly_tol) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' Norm of the last vector is TOO large',w2%norme
        write(out_unitp,*) ' poly_tol: ',para_propa%para_poly%poly_tol
        write(out_unitp,*) ' => npoly or max_poly are TOO small',               &
                      para_propa%para_poly%npoly
        STOP
      END IF

 !write(out_unitp,*) ' Psi before phase shift '
 !CALL ecri_psi(psi=psi)
!     - Phase Shift -----------------------------------
      phase = para_H%E0*para_propa%WPdeltaT
      psi = psi * exp(-cmplx(ZERO,phase,kind=Rkind))

 !write(out_unitp,*) ' Psi after phase shift '
 !CALL ecri_psi(psi=psi)

!    - check norm ------------------
      CALL norm2_psi(psi,GridRep=.FALSE.,BasisRep=.TRUE.)
      IF ( psi%norme > psi%max_norme) THEN
        T  = T + para_propa%WPdeltaT
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' STOP propagation: norm > max_norm',                &
                         psi%norme
        para_propa%march_error   = .TRUE.
        para_propa%test_max_norm = .TRUE.
        STOP
      END IF



      microdeltaT = para_propa%WPdeltaT/                                &
                    real(para_propa%nb_micro,kind=Rkind)
      microphase = phase/                                               &
                    real(para_propa%nb_micro,kind=Rkind)

      phase = ZERO
      microT = ZERO

      DO it=1,para_propa%nb_micro

        microT = microT + microdeltaT
        phase = phase + microphase

        rtj = cmplx(ONE,ZERO,kind=Rkind)
        cdot = psi0Hkpsi0(0)
        DO j=1,j_exit

          rt_tmp =  cmplx(ZERO,-microT/real(j,kind=Rkind),kind=Rkind)
          rtj = rt_tmp * rtj
          cdot = cdot + psi0Hkpsi0(j) * rtj
        END DO
        cdot = cdot * cmplx( cos(phase),-sin(phase),kind=Rkind)
        CALL Write_AutoCorr(no,T+microT,cdot)
      END DO



!-----------------------------------------------------------
      IF (debug) THEN
        CALL norm2_psi(psi)
        write(out_unitp,*) 'norm psi',psi%norme
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------
!STOP

 END SUBROUTINE march_noD
!================================================================
!
!     march cheby
!
!================================================================
      !!@description: march cheby
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      SUBROUTINE march_cheby(T,no,psi,psi0,para_H,para_propa)
      USE mod_system
      USE mod_Op,              ONLY : param_Op, sub_OpPsi,sub_scaledOpPsi

      USE mod_psi_set_alloc,   ONLY : param_psi,ecri_psi
      USE mod_ana_psi,         ONLY : norm2_psi
      USE mod_psi_SimpleOp
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_psi)   :: psi0,psi

!------ working variables ---------------------------------
      TYPE (param_psi), pointer :: w1,w2,w3
      TYPE (param_psi)          :: psi_save

      complex (kind=Rkind) :: cdot
      complex (kind=Rkind) :: psi0Hkpsi0(0:para_propa%para_poly%npoly)
      real (kind=Rkind)    :: microT,T,microdeltaT,phase,microphase
      integer              :: no

      complex (kind=Rkind) :: rt,rt2,rtj,rt_tmp
      integer              :: it,j,jt,jt_exit
      real (kind=Rkind)    :: r,rg
      integer              :: i,i_qaie_corr,icheb
      real (kind=Rkind)    :: norm_exit
      integer              :: max_cheby=10


!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING march_cheby'
        write(out_unitp,*) 'deltaT',para_propa%WPdeltaT
        write(out_unitp,*) 'deltaE',para_propa%para_poly%deltaE
        write(out_unitp,*) 'E0',para_H%E0
        write(out_unitp,*) 'Esc',para_H%Esc
        write(out_unitp,*) 'ncheby',para_propa%para_poly%npoly

        write(out_unitp,*) 'psi%BasisRep psi%GridRep',psi%BasisRep,psi%GridRep

        !write(out_unitp,*) 'psi'
        !CALL ecri_psi(T=T,psi=psi,ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)

      END IF
!-----------------------------------------------------------
    w1 => para_propa%work_WP(1)
    w2 => para_propa%work_WP(2)
    w3 => para_propa%work_WP(3)

    psi_save = psi

    DO icheb=1,max_cheby
      para_propa%march_error = .FALSE.

      para_propa%para_poly%deltaE = para_propa%para_poly%Hmax -          &
                                                para_propa%para_poly%Hmin
      para_propa%para_poly%E0     = para_propa%para_poly%Hmin +          &
                                       HALF * para_propa%para_poly%deltaE
      para_propa%para_poly%alpha  = HALF * para_propa%para_poly%deltaE * &
                                                     para_propa%WPdeltaT
      para_propa%para_poly%Esc    = HALF * para_propa%para_poly%deltaE
      para_H%E0                   = para_propa%para_poly%E0
      para_H%Esc                  = para_propa%para_poly%Esc

      IF (debug) THEN
        write(out_unitp,*) 'deltaT',para_propa%WPdeltaT
        write(out_unitp,*) 'deltaE',para_propa%para_poly%deltaE
        write(out_unitp,*) 'alpha ',para_propa%para_poly%alpha
        write(out_unitp,*) 'E0    ',para_H%E0
        write(out_unitp,*) 'Esc   ',para_H%Esc
        write(out_unitp,*) 'ncheby',para_propa%para_poly%npoly
      END IF

      psi0Hkpsi0(:) = CZERO

      r = max(ONE,HALF * para_propa%para_poly%deltaE * para_propa%WPdeltaT)
      CALL cof(r,para_propa%para_poly%npoly,para_propa%para_poly%coef_poly)

      IF (debug) THEN
        write(out_unitp,*) 'r,deltaE,WPdeltaT',r,para_propa%para_poly%deltaE,para_propa%WPdeltaT
        write(out_unitp,*) 'npoly',para_propa%para_poly%npoly
        !write(out_unitp,*) 'coef_poly',para_propa%para_poly%coef_poly(1:para_propa%para_poly%npoly)
      END IF

      psi0Hkpsi0(0) = Calc_AutoCorr(psi0,psi,para_propa,T,Write_AC=.FALSE.)

      rt  = cmplx(ZERO,-ONE,kind=Rkind)
      rt2 = cmplx(ZERO,-TWO,kind=Rkind)

!     - The first term of the expansion ------------------
      w1  = psi
      psi = psi * para_propa%para_poly%coef_poly(1)


      psi0Hkpsi0(1) =  psi0Hkpsi0(0)

!     - The second term of the expansion -----------------
      write(out_unitp,'(a)',advance='no') 'cheby rec:'
      CALL sub_OpPsi(w1,w2,para_H)
      CALL sub_scaledOpPsi(w1,w2,para_H%E0,para_H%Esc)

      w2  = w2 * rt
      psi = psi + w2 * para_propa%para_poly%coef_poly(2)

      psi0Hkpsi0(2) = Calc_AutoCorr(psi0,w2,para_propa,T,Write_AC=.FALSE.)


!     - The higher terms of the expansion ----------------

      DO jt=3,para_propa%para_poly%npoly

        CALL sub_OpPsi(w2,w3,para_H)
        CALL sub_scaledOpPsi(w2,w3,para_H%E0,para_H%Esc)
        IF (mod(jt,100) == 0) write(out_unitp,'(a)',advance='no') '.'

!        Recurrence relations of the Chebychev expansion:
         w3 = w1 + w3 * rt2


         w1 = w2
         w2 = w3

         psi = psi + w2 * para_propa%para_poly%coef_poly(jt)

         CALL norm2_psi(w2)

         norm_exit = abs(w2%norme*para_propa%para_poly%coef_poly(jt))
         jt_exit = jt
         IF (debug) write(out_unitp,*) 'jt,norms',jt,norm_exit

         IF (norm_exit > TEN**15) THEN
           write(out_unitp,*) ' WARNING: Norm^2 of the vector is TOO large (> 10^15)',jt,norm_exit
           para_propa%march_error = .TRUE.
           EXIT
         END IF

         psi0Hkpsi0(jt) = Calc_AutoCorr(psi0,w2,para_propa,T,Write_AC=.FALSE.)

         IF (norm_exit < para_propa%para_poly%poly_tol) EXIT

      END DO
      write(out_unitp,*) 'jt_exit,norms',jt_exit,abs(w2%norme),norm_exit
      para_propa%para_poly%npoly_Opt = jt_exit

      IF (.NOT. para_propa%march_error .AND. norm_exit > para_propa%para_poly%poly_tol) THEN
        write(out_unitp,*) ' ERROR in march_cheby'
        write(out_unitp,*) ' WARNING: Norm^2 of the last vector TOO large',jt_exit,norm_exit
        para_propa%march_error = .TRUE.
      END IF

      IF (para_propa%march_error) THEN
        write(out_unitp,*) ' Old Hmin,Hmax',para_propa%para_poly%Hmin,para_propa%para_poly%Hmax

        para_propa%para_poly%Hmax = para_propa%para_poly%Hmax +         &
                                  para_propa%para_poly%deltaE * ONETENTH
        !para_propa%para_poly%Hmin = para_propa%para_poly%Hmin -         &
        !                          para_propa%para_poly%deltaE * ONETENTH

        write(out_unitp,*) ' New Hmin,Hmax',para_propa%para_poly%Hmin,para_propa%para_poly%Hmax

        psi = psi_save
      ELSE
        EXIT
      END IF

    END DO

    IF (para_propa%march_error) THEN
      write(out_unitp,*) ' ERROR in march_cheby'
      write(out_unitp,*) ' It cannot converge ...'
      write(out_unitp,*) ' => Reduce the time step !!'
      STOP
    END IF

    CALL dealloc_psi(psi_save,delete_all=.TRUE.)


!     - Phase Shift -----------------------------------
      phase = para_H%E0*para_propa%WPdeltaT
      psi = psi * exp(-cmplx(ZERO,phase,kind=Rkind))



      microdeltaT = para_propa%WPdeltaT/                                &
                    real(para_propa%nb_micro,kind=Rkind)
      microphase = phase/real(para_propa%nb_micro,kind=Rkind)

      phase = ZERO
      microT = ZERO

      !write(out_unitp,*) 'para_propa%nb_micro',para_propa%nb_micro
      !write(out_unitp,*) 'microdeltaT',microdeltaT
      !write(out_unitp,*) 'microphase',microphase


      DO it=1,para_propa%nb_micro

        microT = microT + microdeltaT
        phase = phase + microphase

        r = HALF * para_propa%para_poly%deltaE * microT
        CALL cof(r,jt_exit,para_propa%para_poly%coef_poly)

        cdot = cmplx(ZERO,ZERO,kind=Rkind)
        !write(out_unitp,*) 'jt,cdot',0,cdot
        DO jt=1,jt_exit

          cdot = cdot + psi0Hkpsi0(jt) *                                &
            cmplx(para_propa%para_poly%coef_poly(jt),ZERO,kind=Rkind)
          !write(out_unitp,*) 'jt,psi0Hkpsi0(jt),coef_poly(jt),cdot',jt,psi0Hkpsi0(jt),para_propa%para_poly%coef_poly(jt),cdot

        END DO
        cdot = cdot * cmplx( cos(phase),-sin(phase),kind=Rkind)
        !write(out_unitp,*) 'T,cdot,no,type_corr',T+microT,cdot,no,para_propa%type_corr
        CALL Write_AutoCorr(no,T+microT,cdot)
      END DO
      CALL flush_perso(no)




!-----------------------------------------------------------
      IF (debug) THEN
        !write(out_unitp,*) 'psi'
        !CALL ecri_psi(T=T+para_propa%WPdeltaT,psi=psi,ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)
        write(out_unitp,*) 'END march_cheby'
      END IF
!-----------------------------------------------------------


      END SUBROUTINE march_cheby
      SUBROUTINE march_cheby_old(T,no,psi,psi0,para_H,para_propa)
      USE mod_system
      USE mod_Op,              ONLY : param_Op, sub_OpPsi,sub_scaledOpPsi

      USE mod_psi_set_alloc,   ONLY : param_psi,ecri_psi
      USE mod_ana_psi,         ONLY : norm2_psi
      USE mod_psi_SimpleOp
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_psi)   :: psi0,psi

!------ working variables ---------------------------------
      TYPE (param_psi), pointer :: w1,w2,w3

      complex (kind=Rkind) :: cdot
      complex (kind=Rkind) :: psi0Hkpsi0(0:para_propa%para_poly%npoly)
      real (kind=Rkind)    :: microT,T,microdeltaT,phase,microphase
      integer              :: no

      complex (kind=Rkind) :: rt,rt2,rtj,rt_tmp
      integer              :: it,j,jt,jt_exit
      real (kind=Rkind)    :: r,rg
      integer              :: i,i_qaie_corr
      real (kind=Rkind)    :: norm_exit


!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING march_cheby_old'
        write(out_unitp,*) 'deltaT',para_propa%WPdeltaT
        write(out_unitp,*) 'deltaE',para_propa%para_poly%deltaE
        write(out_unitp,*) 'E0',para_H%E0
        write(out_unitp,*) 'Esc',para_H%Esc
        write(out_unitp,*) 'ncheby',para_propa%para_poly%npoly

        write(out_unitp,*) 'psi%BasisRep psi%GridRep',psi%BasisRep,psi%GridRep

        !write(out_unitp,*) 'psi'
        !CALL ecri_psi(T=T,psi=psi,ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)

      END IF
!-----------------------------------------------------------

      w1 => para_propa%work_WP(1)
      w2 => para_propa%work_WP(2)
      w3 => para_propa%work_WP(3)

      psi0Hkpsi0(:) = cmplx(ZERO,ZERO,kind=Rkind)

      r = HALF * para_propa%para_poly%deltaE * para_propa%WPdeltaT
      CALL cof(r,para_propa%para_poly%npoly,                            &
                 para_propa%para_poly%coef_poly)
      !write(out_unitp,*) 'r,deltaE,WPdeltaT',r,para_propa%para_poly%deltaE,para_propa%WPdeltaT
      !write(out_unitp,*) 'npoly,coef_poly',para_propa%para_poly%npoly,  &
      !  para_propa%para_poly%coef_poly(1:para_propa%para_poly%npoly)

      psi0Hkpsi0(0) = Calc_AutoCorr(psi0,psi,para_propa,T,Write_AC=.FALSE.)

      rt  = cmplx(ZERO,-ONE,kind=Rkind)
      rt2 = cmplx(ZERO,-TWO,kind=Rkind)

!     - The first term of the expansion ------------------
      w1  = psi
      psi = psi * para_propa%para_poly%coef_poly(1)


      psi0Hkpsi0(1) =  psi0Hkpsi0(0)

!     - The second term of the expansion -----------------
      write(out_unitp,'(a)',advance='no') 'cheby rec:'
      CALL sub_OpPsi(w1,w2,para_H)
      CALL sub_scaledOpPsi(w1,w2,para_H%E0,para_H%Esc)

      w2  = w2 * rt
      psi = psi + w2 * para_propa%para_poly%coef_poly(2)

      psi0Hkpsi0(2) = Calc_AutoCorr(psi0,w2,para_propa,T,Write_AC=.FALSE.)


!     - The higher terms of the expansion ----------------

      DO jt=3,para_propa%para_poly%npoly

        CALL sub_OpPsi(w2,w3,para_H)
        CALL sub_scaledOpPsi(w2,w3,para_H%E0,para_H%Esc)
        IF (mod(jt,100) == 0) write(out_unitp,'(a)',advance='no') '.'

!        Recurrence relations of the Chebychev expansion:
         w3 = w1 + w3 * rt2


         w1 = w2
         w2 = w3

         psi = psi + w2 * para_propa%para_poly%coef_poly(jt)

         CALL norm2_psi(w2)

         norm_exit = abs(w2%norme*para_propa%para_poly%coef_poly(jt))
         jt_exit = jt
         IF (debug) write(out_unitp,*) 'jt,norms',jt,norm_exit

          IF (norm_exit > TEN**15) THEN
            write(out_unitp,*) ' ERROR in march_cheby'
            write(out_unitp,*) ' Norm of the vector is TOO large (> 10^15)',    &
                                                jt,norm_exit
           write(out_unitp,*) ' => Reduce the time step !!'
           STOP
         END IF
         IF (norm_exit < para_propa%para_poly%poly_tol) EXIT

         psi0Hkpsi0(jt) = Calc_AutoCorr(psi0,w2,para_propa,T,Write_AC=.FALSE.)

      END DO
      write(out_unitp,*) 'jt_exit,norms',jt_exit,abs(w2%norme),norm_exit

      IF (norm_exit > para_propa%para_poly%poly_tol) THEN
        write(out_unitp,*) ' ERROR in march_cheby'
        write(out_unitp,*) ' Norm of the last vector TOO large',norm_exit
        write(out_unitp,*) ' poly_tol: ',para_propa%para_poly%poly_tol
        write(out_unitp,*) ' => npoly TOO small',para_propa%para_poly%npoly
        STOP
      END IF

!     - Phase Shift -----------------------------------
      phase = para_H%E0*para_propa%WPdeltaT
      psi = psi * exp(-cmplx(ZERO,phase,kind=Rkind))



      microdeltaT = para_propa%WPdeltaT/                                &
                    real(para_propa%nb_micro,kind=Rkind)
      microphase = phase/real(para_propa%nb_micro,kind=Rkind)

      phase = ZERO
      microT = ZERO

      !write(out_unitp,*) 'para_propa%nb_micro',para_propa%nb_micro
      !write(out_unitp,*) 'microdeltaT',microdeltaT
      !write(out_unitp,*) 'microphase',microphase


      DO it=1,para_propa%nb_micro

        microT = microT + microdeltaT
        phase = phase + microphase

        r = HALF * para_propa%para_poly%deltaE * microT
        CALL cof(r,jt_exit,para_propa%para_poly%coef_poly)

        cdot = cmplx(ZERO,ZERO,kind=Rkind)
        !write(out_unitp,*) 'jt,cdot',0,cdot
        DO jt=1,jt_exit

          cdot = cdot + psi0Hkpsi0(jt) *                                &
            cmplx(para_propa%para_poly%coef_poly(jt),ZERO,kind=Rkind)
          !write(out_unitp,*) 'jt,psi0Hkpsi0(jt),coef_poly(jt),cdot',jt,psi0Hkpsi0(jt),para_propa%para_poly%coef_poly(jt),cdot

        END DO
        cdot = cdot * cmplx( cos(phase),-sin(phase),kind=Rkind)
        !write(out_unitp,*) 'T,cdot,no,type_corr',T+microT,cdot,no,para_propa%type_corr
        CALL Write_AutoCorr(no,T+microT,cdot)
      END DO
      CALL flush_perso(no)




!-----------------------------------------------------------
      IF (debug) THEN
        !write(out_unitp,*) 'psi'
        !CALL ecri_psi(T=T+para_propa%WPdeltaT,psi=psi,ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)
        write(out_unitp,*) 'END march_cheby_old'
      END IF
!-----------------------------------------------------------


      END SUBROUTINE march_cheby_old
!==============================================================
!
! nOD propagator in imagynary time
!
!==============================================================
      !!@description: nOD propagator in imagynary time
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      SUBROUTINE march_nOD_im(T,no,psi,psi0,w1,w2,para_H,para_propa)
      USE mod_system
      USE mod_Op,              ONLY : param_Op, sub_OpPsi,sub_scaledOpPsi

      USE mod_psi_set_alloc,   ONLY : param_psi,ecri_psi
      USE mod_ana_psi,         ONLY : norm2_psi
      USE mod_psi_SimpleOp
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_psi)   :: psi0,psi,w1,w2


!------ working variables ---------------------------------
      integer              :: no
      real (kind=Rkind)    :: rtj
      integer              :: j
      real  (kind=Rkind)   :: T,signDT
      complex (kind=Rkind) :: E0


!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING march_nOD_im'
        write(out_unitp,*) 'deltaT',para_propa%WPdeltaT
        write(out_unitp,*) 'nOD',para_propa%para_poly%npoly
        write(out_unitp,*) 'type_WPpropa',para_propa%type_WPpropa
        write(out_unitp,*) 'para_H%E0',para_H%E0
      END IF
!-----------------------------------------------------------



      signDT = -ONE
      IF (para_propa%type_WPpropa < 0) signDT = ONE


      w1 = psi
      DO j=1,para_propa%para_poly%max_poly

        rtj =  signDT*para_propa%WPdeltaT/real(j,kind=Rkind)

        CALL sub_OpPsi(w1,w2,para_H)
        CALL sub_scaledOpPsi(w1,w2,para_H%E0,ONE)

        w1  = w2 * rtj
        psi = psi + w1

        CALL norm2_psi(w1)
        CALL norm2_psi(psi)
        !write(out_unitp,*) 'j,wi%n/psi%n',j,w1%norme/psi%norme
        IF (w1%norme/psi%norme < para_propa%para_poly%poly_tol) EXIT

      END DO
      IF (para_propa%write_iter .OR. debug)                             &
                  write(out_unitp,*) 'j,wi%n/psi%n',j,w1%norme/psi%norme

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END march_nOD_im'
      END IF
!-----------------------------------------------------------

      END SUBROUTINE march_nOD_im
!================================================================
!
!     march for relaxation with optimal DeltaT (first Order)
!
!================================================================
      !!@description: march for relaxation with optimal DeltaT (first Order)
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      SUBROUTINE march_FOD_Opti_im(psi,Ene0,T,it,para_H,para_propa)
      USE mod_system
      USE mod_Op,              ONLY : param_Op, sub_OpPsi,sub_scaledOpPsi
      USE mod_psi_set_alloc,   ONLY : param_psi,ecri_psi,dealloc_psi
      USE mod_psi_Op,          ONLY : Overlap_psi1_psi2
      USE mod_ana_psi,         ONLY : norm2_psi

      USE mod_psi_SimpleOp
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)        :: para_H

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa)     :: para_propa
      TYPE (param_psi)       :: psi

      TYPE (param_psi)       :: Hpsi,H2psi,H3psi,g

      real (kind=Rkind)      :: Ene0
      logical                :: cplx
      complex (kind=Rkind)   :: CS
      real (kind=Rkind)      :: RS

!------ working parameters --------------------------------
      complex (kind=Rkind)   :: cdot
      integer                :: it,it_all,no
      integer                :: i,ii,j,iqa
      integer                :: max_ecri
      real (kind=Rkind)      :: DeltaT,T      ! time
      real (kind=Rkind)      :: DeltaE,Deltapsi,epsi,sign,normeg
      real (kind=Rkind)      :: avH1,avH2,avH3,A,B,C,D,DT1,DT2,S,E1,E2
      complex (kind=Rkind)   :: CavH1,CavH2,CavH3

      integer                ::   nioWP


!----- for debuging --------------------------------------------------
       logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING march_FOD_Opti_im'
        write(out_unitp,*) 'Tmax,deltaT',para_propa%WPTmax,para_propa%WPdeltaT
        write(out_unitp,*) 'Hmin,Hmax',para_propa%para_poly%Hmin,               &
                                para_propa%para_poly%Hmax
        write(out_unitp,*)
       END IF
!-----------------------------------------------------------

      IF (para_H%cplx) THEN
        write(out_unitp,*) ' ERROR in march_FOD_Opti_im'
        write(out_unitp,*) ' This subroutine cannot be use with complex Hmiltonian'
        write(out_unitp,*) ' You can use "march_nOD_im" instead !'
        STOP
      END IF

!------ initialization -------------------------------------
      epsi = para_propa%para_poly%poly_tol


      Hpsi  = psi
      H2psi = psi
      g     = psi
      IF (psi%cplx) H3psi = psi
      !- Optimal DeltaT and energy ----------------------------

      ! Hpsi = H.psi
      CALL sub_OpPsi(psi,Hpsi,para_H)
      CALL sub_scaledOpPsi(psi,Hpsi,para_H%E0,ONE)

      ! H2psi = H.H.psi = H.Hpsi
      CALL sub_OpPsi(Hpsi,H2psi,para_H)
      CALL sub_scaledOpPsi(Hpsi,H2psi,para_H%E0,ONE)

      ! H3psi = H.H.H.psi = H.H2psi
      ! <psi|H3|psi> = <Hspi|H2psi> when psi is real
      CALL Overlap_psi1_psi2(CavH1,psi,Hpsi)
      CALL Overlap_psi1_psi2(CavH2,psi,H2psi)
      IF (psi%cplx) THEN
        CALL sub_OpPsi(H2psi,H3psi,para_H)
        CALL sub_scaledOpPsi(H2psi,H3psi,para_H%E0,ONE)
        CALL Overlap_psi1_psi2(CavH3,psi,H3psi)
      ELSE
        CALL Overlap_psi1_psi2(CavH3,Hpsi,H2psi)
      END IF
      avH1 = CavH1
      avH2 = CavH2
      avH3 = CavH3

      IF (debug) write(out_unitp,*) 'avH..',avH1,avH2,avH3


      A =  avH1*avH1 - avH2
      B = -avH1*avH2 + avH3
      C =  avH2*avH2 - avH3*avH1
      D = B*B - FOUR*A*C
      DT1 = (-B+sqrt(D))/(C+C)
      DT2 = (-B-sqrt(D))/(C+C)


      S = ONE-TWO*DT1*avH1+DT1*DT1*avH2
      E1= avH1-TWO*DT1*avH2+DT1*DT1*avH3
      E1= E1/S
      !write(out_unitp,*) 'DT1',DT1,S,E1*get_Conv_au_TO_unit('E','cm-1')
      S = ONE-TWO*DT2*avH1+DT2*DT2*avH2
      E2= avH1-TWO*DT2*avH2+DT2*DT2*avH3
      E2= E2/S
      !write(out_unitp,*) 'DT2',DT2,S,E2*get_Conv_au_TO_unit('E','cm-1')

      IF (E1<E2) THEN
        DeltaT = DT1
      ELSE
        DeltaT = DT2
      END IF

      S = ONE - TWO*DeltaT*avH1 + DeltaT**2*avH2
      DO WHILE (abs(S) < ONETENTH**6)
        write(out_unitp,*) 'WARNING: S<',ONETENTH**6
        write(out_unitp,*) 'WP it, DeltaT,S',it,DeltaT,S
        DeltaT = DeltaT * HALF
        S = ONE - TWO*DeltaT*avH1 + DeltaT**2*avH2
      END DO

      Ene0 = avH1 - TWO*DeltaT*avH2 + DeltaT**2*avH3
      Ene0 = Ene0/S
      DeltaE = Ene0 - avH1

      !- Optimal DeltaT and energy ----------------------------

      !- propagation ------------------------------------------
      g = psi         *(-Ene0)
      g = g +    Hpsi *(ONE + Ene0*DeltaT)
      g = g +   H2psi *(          -DeltaT)
      CALL norm2_psi(g)
      normeg = sqrt(g%norme)

      H2psi = psi + Hpsi * (-DeltaT)
      psi = H2psi * (ONE/sqrt(S))

      Deltapsi = -A/S * DeltaT**2


      !- propagation ------------------------------------------
      IF (para_propa%write_iter .OR. debug) THEN
        write(out_unitp,*) 'WP it sqrt(norme g)',it,                &
            sqrt(g%norme)*get_Conv_au_TO_unit('E','cm-1')
      END IF

      IF (debug) THEN
        write(out_unitp,*) 'WP it, E',it,                               &
                                  Ene0*get_Conv_au_TO_unit('E','cm-1'), &
                                DeltaE*get_Conv_au_TO_unit('E','cm-1')
        write(out_unitp,*) 'WP it sqrt(norme g)',it,                &
            sqrt(g%norme)*get_Conv_au_TO_unit('E','cm-1')
        write(out_unitp,*) 'WP it, DeltaT,S',it,DeltaT,S
      END IF

      T = T + DeltaT

      CALL dealloc_psi(Hpsi)
      CALL dealloc_psi(H2psi)
      CALL dealloc_psi(H3psi)
      CALL dealloc_psi(g)

!----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'END march_FOD_Opti_im'
       END IF
!----------------------------------------------------------


      END SUBROUTINE march_FOD_Opti_im

!================================================================
!
!    march nOD propagation with
!         a time dependant pulse in Hamiltonian (W(t))
!
!================================================================
      SUBROUTINE march_new0_noD_field(T,WP,nb_WP,print_Op,              &
                                      para_H,para_Dip,                  &
                                      para_field,para_propa)
      USE mod_system

      USE mod_field,           ONLY : param_field,sub_dnE
      USE mod_Op,              ONLY : param_Op, sub_OpPsi,sub_scaledOpPsi

      USE mod_psi_set_alloc,   ONLY : param_psi,ecri_psi
      USE mod_ana_psi,          ONLY : norm2_psi
      USE mod_psi_SimpleOp
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H
      TYPE (param_Op)   :: para_Dip(3)

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_field) :: para_field

      integer          :: nb_WP
      TYPE (param_psi) :: WP(nb_WP)
      TYPE (param_psi), pointer :: w1,w2

!----- for printing --------------------------------------------------
      logical ::print_Op



!------ working parameters --------------------------------
      real (kind=Rkind) :: fac_k,phase

      complex (kind=Rkind)   :: Overlap,cdot,rt,rti,norm,coef,coef_ip(3)
      real (kind=Rkind) :: limit

      integer       :: it,it_max,no,max_der,nb_der
      integer       :: i,j,k,m,jt,ip,iq
      real (kind=Rkind) :: T      ! time
      real (kind=Rkind) :: T_Delta! time+deltaT
      integer       :: max_ecri

      integer  ::   nioWP

!     logical :: taylor=.TRUE.
      logical :: taylor=.FALSE.
!----- for the field --------------------------------------------------
      real (kind=Rkind), allocatable :: tab_dnE(:,:)
      real (kind=Rkind)    :: dnE(3)
      real (kind=Rkind)    :: ww(3)
      logical :: make_field
!----- for the field --------------------------------------------------

!----- function -------------------------------------------------------
      real (kind=Rkind)    :: binomial,factor


!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='march_nOd_field'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Tmax,deltaT',para_propa%WPTmax,para_propa%WPdeltaT
        write(out_unitp,*) 'Hmin,Hmax',para_propa%para_poly%Hmin,               &
                                para_propa%para_poly%Hmax
        write(out_unitp,*)
        write(out_unitp,*) 'nb_ba,nb_qa',WP(1)%nb_ba,WP(1)%nb_qa
        write(out_unitp,*) 'nb_bi',WP(1)%nb_bi
        write(out_unitp,*)

        DO i=1,nb_WP
          CALL norm2_psi(WP(i),.FALSE.,.TRUE.,.FALSE.)
          write(out_unitp,*) 'norm WP',i,WP(i)%norme

          write(out_unitp,*) 'WP BasisRep',i
          CALL ecri_psi(T=ZERO,psi=WP(i),                               &
                        ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)
          write(out_unitp,*) 'WP GridRep',i
          CALL ecri_psi(T=ZERO,psi=WP(i),                               &
                        ecri_GridRep=.TRUE.,ecri_BasisRep=.FALSE.)
        END DO
      END IF
!-----------------------------------------------------------


!-----------------------------------------------------------
      w1 => para_propa%work_WP(para_propa%para_poly%npoly+1)
      w2 => para_propa%work_WP(para_propa%para_poly%npoly+2)
!-----------------------------------------------------------

!-----------------------------------------------------------
      max_der = 0
      nb_der = 3*para_propa%para_poly%npoly
      IF (para_field%max_der >= 0)                                      &
           nb_der = min(para_field%max_der,para_propa%para_poly%npoly-1)

      CALL alloc_NParray(tab_dnE,(/nb_der,3/),                            &
                      "tab_dnE","march_new0_noD_field",(/0,1/))
      tab_dnE(:,:) = ZERO

      DO k=0,nb_der
        CALL sub_dnE(dnE,k,T,para_field)
        tab_dnE(k,:) = -dnE(:)
      END DO
      max_der = nb_der
      write(out_unitp,*) 'max_der,npoly',max_der,para_propa%para_poly%npoly

!-----------------------------------------------------------
      DO j=1,nb_WP
!       Derivative calculation
        fac_k = ONE
        DO k=0,para_propa%para_poly%npoly-1

          IF (k == 0) THEN
            para_propa%work_WP(k) = WP(j)
          ELSE
!           initialization with the H contribution
            CALL sub_OpPsi(para_propa%work_WP(k-1),w1,para_H)
            CALL sub_scaledOpPsi(para_propa%work_WP(k-1),w1,            &
                                 para_H%E0,ONE)
!           ici pas de binome car m=0
            coef = -EYE
            para_propa%work_WP(k) = w1 * coef
            DO m=0,k-1
!             Dip contributions
              DO ip=1,3
                IF (.NOT. para_field%pola_xyz(ip)) CYCLE
                coef = -EYE*binomial(k-1,m) * tab_dnE(m,ip)
                CALL sub_OpPsi(para_propa%work_WP(k-1-m),w1,            &
                              para_Dip(ip))
                para_propa%work_WP(k) = para_propa%work_WP(k) +         &
                                         coef * w1
              END DO
            END DO
          END IF


          IF (taylor) THEN
!           Just the Taylor serie
            IF (k == 0) CYCLE ! because WP(j) is already the zero-order
            fac_k = fac_k * para_propa%WPdeltaT / real(k,kind=Rkind)
            w1 = para_propa%work_WP(k) * fac_k
          ELSE
!           Infinite sum of the derivative of the field
            CALL sub_OpPsi(para_propa%work_WP(k),w2,para_H)
            CALL sub_scaledOpPsi(para_propa%work_WP(k),w2,para_H%E0,ONE)
            coef = para_propa%WPdeltaT**(k+1) / factor(k+1)
            w1 = w2 * coef
            coef_ip(:) = ZERO
!           DO m=k+1,max_der
            DO m=k+1,max_der
!             coef = binomial(m-1,k)*para_propa%WPdeltaT**m/factor(m)
              coef = para_propa%WPdeltaT**m /                           &
                  (real(m,kind=Rkind) * factor(k) * factor(m-1-k) )
              coef_ip(:) = coef_ip(:) + coef*tab_dnE(m-1-k,:)
!             write(out_unitp,*) 'k,m',k,m,abs(coef)
              IF (abs(coef) <                                           &
                      para_propa%para_poly%poly_tol*ONETENTH**0) EXIT
            END DO
            DO ip=1,3
              IF (.NOT. para_field%pola_xyz(ip)) CYCLE
              CALL sub_OpPsi(para_propa%work_WP(k),w2,para_Dip(ip))
              w1 = w1 + w2 * coef_ip(ip)
            END DO
          END IF


          WP(j) = WP(j) + w1 * (-EYE)
          CALL norm2_psi(w1)
          write(out_unitp,*) 'norme dkpsi k',k,w1%norme
          IF (w1%norme < para_propa%para_poly%poly_tol) EXIT

          IF (w1%norme > TEN**15) THEN
             write(out_unitp,*) ' ERROR in ',name_sub
             write(out_unitp,*) ' The norme of dnpsi(k) is huge !',w1%norme
             write(out_unitp,*) ' iteration k:',k
             write(out_unitp,*) ' Reduce the time step, WPDeltaT:',             &
                             para_propa%WPdeltaT
             para_propa%march_error = .TRUE.
             STOP
           END IF

        END DO

!       - Phase Shift -----------------
        phase = para_H%E0*para_propa%WPdeltaT
        WP(j) = WP(j) * exp(-EYE*phase)

!       - check norm ------------------
        CALL norm2_psi(WP(j),GridRep=.FALSE.,BasisRep=.TRUE.)
        IF ( WP(j)%norme > WP(j)%max_norme) THEN
           T  = T + para_propa%WPdeltaT
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) ' STOP propagation: norm > max_norm',             &
                         WP(j)%norme
           para_propa%march_error   = .TRUE.
           para_propa%test_max_norm = .TRUE.
           STOP
        END IF

      END DO

      CALL dealloc_NParray(tab_dnE,"tab_dnE","march_new0_noD_field")

!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------


      END SUBROUTINE march_new0_noD_field
!================================================================
!
!    march nOD propagation with
!         a time dependant pulse in Hamiltonian (W(t))
!
!================================================================
      SUBROUTINE march_new_noD_field(T,WP,nb_WP,print_Op,               &
                                     para_H,para_Dip,para_field,para_propa)
      USE mod_system

      USE mod_field,           ONLY : param_field,sub_dnE
      USE mod_Op,              ONLY : param_Op, sub_OpPsi,sub_scaledOpPsi

      USE mod_psi_set_alloc,   ONLY : param_psi,ecri_psi
      USE mod_ana_psi,         ONLY : norm2_psi
      USE mod_psi_SimpleOp
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H
      TYPE (param_Op)   :: para_Dip(3)

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_field) :: para_field

      integer          :: nb_WP
      TYPE (param_psi) :: WP(nb_WP)
      TYPE (param_psi), pointer :: w1,w2

!----- for printing --------------------------------------------------
      logical ::print_Op



!------ working parameters --------------------------------

      complex (kind=Rkind)   :: Overlap,cdot,rt,rti,fac_k,norm
      real (kind=Rkind) :: limit,phase

      integer       :: it,it_max,no,max_der,nb_der
      integer       :: i,j,k,jt,ip,iq
      real (kind=Rkind) :: T      ! time
      real (kind=Rkind) :: T_Delta! time+deltaT
      integer       :: max_ecri
      real (kind=Rkind) :: IntE(3)
      real (kind=Rkind) :: max_norm

      integer  ::   nioWP
      real (kind=Rkind)    :: fac
      real (kind=Rkind)    :: binomial,factor ! functions

      logical, parameter :: taylor = .TRUE.
!     logical, parameter :: taylor = .FALSE.

!----- for the field --------------------------------------------------
      real (kind=Rkind), allocatable :: tab_dnE(:,:)
      real (kind=Rkind)    :: dnE(3)
      real (kind=Rkind)    :: ww(3)
      logical :: make_field
!----- for the field --------------------------------------------------


!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='march_new_nOd_field'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Tmax,deltaT',para_propa%WPTmax,para_propa%WPdeltaT
        write(out_unitp,*) 'Hmin,Hmax',para_propa%para_poly%Hmin,               &
                                para_propa%para_poly%Hmax
        write(out_unitp,*)
        write(out_unitp,*) 'nb_ba,nb_qa',WP(1)%nb_ba,WP(1)%nb_qa
        write(out_unitp,*) 'nb_bi',WP(1)%nb_bi
        write(out_unitp,*)

        DO i=1,nb_WP
          CALL norm2_psi(WP(i),.FALSE.,.TRUE.,.FALSE.)
          write(out_unitp,*) 'norm WP',i,WP(i)%norme

          write(out_unitp,*) 'WP BasisRep',i
          CALL ecri_psi(T=ZERO,psi=WP(i),                               &
                        ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)
          write(out_unitp,*) 'WP GridRep',i
          CALL ecri_psi(T=ZERO,psi=WP(i),                               &
                        ecri_GridRep=.TRUE.,ecri_BasisRep=.FALSE.)
        END DO
      END IF
!-----------------------------------------------------------


!-----------------------------------------------------------
      w1 => para_propa%work_WP(para_propa%para_poly%npoly+1)
      w2 => para_propa%work_WP(para_propa%para_poly%npoly+2)

      max_der = 0
      nb_der = para_propa%para_poly%npoly-1
      IF (para_field%max_der >= 0)                                      &
           nb_der = min(para_field%max_der,para_propa%para_poly%npoly-1)

      CALL alloc_NParray(tab_dnE,(/nb_der,3/),                            &
                      "tab_dnE",name_sub,(/0,1/))
      tab_dnE(:,:) = ZERO

!     - propagation ----------------
      fac_k = para_propa%WPdeltaT
      CALL sub_dnE(dnE,0,T,para_field)
      tab_dnE(0,:) = -dnE(:)
      DO k=1,nb_der
        CALL sub_dnE(dnE,k,T,para_field)
        tab_dnE(k,:) = -dnE(:)
        fac_k = fac_k * para_propa%WPdeltaT/real(k,kind=Rkind)
        limit = maxval(abs(dnE(:)))*fac_k
        IF (limit > para_propa%para_poly%poly_tol) max_der=k
!       write(out_unitp,*) 'k,max_der,limit',k,max_der,limit
      END DO
      max_der=para_propa%para_poly%npoly-1
      write(out_unitp,*) 'max_der',max_der

      DO j=1,nb_WP
        para_propa%work_WP(0) = WP(j)
        max_norm = ZERO
        DO i=0,para_propa%para_poly%npoly-1

          CALL sub_OpPsi(para_propa%work_WP(i),w2,para_H)
          CALL sub_scaledOpPsi(para_propa%work_WP(i),w2,para_H%E0,ONE)

          para_propa%work_WP(i+1) = w2  ! H.psi(i)


          fac_k = cmplx(ZERO,-para_propa%WPdeltaT /                     &
                            real(i+1,kind=Rkind),kind=Rkind)


          para_propa%work_WP(i+1) = para_propa%work_WP(i+1) * fac_k

          DO ip=1,3
            IF (.NOT. para_field%pola_xyz(ip)) CYCLE
            fac_k = cmplx(ZERO,-para_propa%WPdeltaT /                   &
                                 real(i+1,kind=Rkind),kind=Rkind)

            w1 = para_propa%work_WP(i) * (tab_dnE(0,ip)* fac_k)

            DO k=1,min(i,max_der)
              fac_k = fac_k * para_propa%WPdeltaT/real(k,kind=Rkind)
              rt = tab_dnE(k,ip)*fac_k
              w1 = w1 + para_propa%work_WP(i-k) * rt
            END DO

            CALL sub_OpPsi(w1,w2,para_Dip(ip))

            para_propa%work_WP(i+1) = para_propa%work_WP(i+1) + w2
          END DO

          IF (taylor) THEN
            w1 = para_propa%work_WP(i+1)
          ELSE
!           H contribution
            CALL sub_OpPsi(para_propa%work_WP(i),w2,para_H)
            CALL sub_scaledOpPsi(para_propa%work_WP(i),w2,para_H%E0,ONE)
            w1 = w2 * (para_propa%WPdeltaT/real(i+1,kind=Rkind))

!           Dip.field contributions
            IntE(:) = ZERO
            DO k=i+1,max_der
              fac = para_propa%WPdeltaT**(k-i) /                        &
                   (real(k,kind=Rkind)*factor(k-i-1))
              IntE(:) = IntE(:) + fac * tab_dnE(k-1-i,:)
              IF (abs(fac) <                                            &
                  para_propa%para_poly%poly_tol*ONETENTH**2) EXIT
            END DO
            DO ip=1,3
              IF (.NOT. para_field%pola_xyz(ip)) CYCLE
              CALL sub_OpPsi(para_propa%work_WP(i),w2,para_Dip(ip))
              w1 = w1 + w2 * IntE(ip)
            END DO

            w1 = w1 * (-EYE)

          END IF

          WP(j) = WP(j) + w1

          CALL norm2_psi(w1)
          IF (w1%norme > max_norm) max_norm = w1%norme
          write(out_unitp,*) 'norms of w1(i),max_norm',i,w1%norme,max_norm

          IF (w1%norme < para_propa%para_poly%poly_tol) EXIT
          IF (w1%norme > TEN**15) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' The norme of dnpsi(i+1) is huge !',w1%norme
            write(out_unitp,*) ' iteration i:',i
            write(out_unitp,*) ' Reduce the time step, WPDeltaT:',              &
                             para_propa%WPdeltaT
            para_propa%march_error = .TRUE.
            STOP
          END IF
        END DO



        IF (print_Op) write(out_unitp,*) 'T,max_der,ipoly,norme',               &
                                  T,max_der,i,w1%norme

        para_propa%work_WP(0) =                                         &
           para_propa%work_WP(min(i+1,para_propa%para_poly%npoly))

!       - Phase Shift -----------------
        phase = para_H%E0*para_propa%WPdeltaT
        WP(j) = WP(j) * exp(-EYE*phase)

!       - check norm ------------------
        CALL norm2_psi(WP(j),GridRep=.FALSE.,BasisRep=.TRUE.)
        IF ( WP(j)%norme > WP(j)%max_norme) THEN
          T  = T + para_propa%WPdeltaT
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' STOP propagation: norm > max_norm',              &
                         WP(j)%norme
          para_propa%march_error   = .TRUE.
          para_propa%test_max_norm = .TRUE.
          STOP
        END IF

      END DO

      CALL dealloc_NParray(tab_dnE,"tab_dnE",name_sub)

!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------


      END SUBROUTINE march_new_noD_field
!================================================================
!
!    march nOD propagation with
!         a time dependant pulse in Hamiltonian (W(t))
!
!================================================================
      SUBROUTINE march_split_field(T,WP,nb_WP,print_Op,                 &
                                 para_H,para_Dip,                       &
                                 para_field,para_propa)
      USE mod_system

      USE mod_field,           ONLY : param_field,sub_dnE
      USE mod_Op,              ONLY : param_Op, sub_OpPsi,sub_OpiPsi,sub_scaledOpPsi

      USE mod_psi_B_TO_G,      ONLY : sub_PsiBasisRep_TO_GridRep,sub_PsiGridRep_TO_BasisRep
      USE mod_psi_set_alloc,   ONLY : param_psi,ecri_psi
      USE mod_ana_psi,         ONLY : norm2_psi
      USE mod_psi_SimpleOp
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H
      TYPE (param_Op)   :: para_Dip(3)

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_field) :: para_field

      integer          :: nb_WP
      TYPE (param_psi) :: WP(nb_WP)
      TYPE (param_psi), pointer :: w1,w2,w3
      TYPE (param_psi), pointer :: expT,expV

!----- for printing --------------------------------------------------
      logical ::print_Op



!------ working parameters --------------------------------

      complex (kind=Rkind)   :: Overlap,cdot,rt,rti,fac_k,norm
      real (kind=Rkind) :: limit

      integer       :: it,it_max,no,max_der,nb_der
      integer       :: i,j,k,jt,ip,iq
      real (kind=Rkind) :: T      ! time
      real (kind=Rkind) :: T_Delta! time+deltaT
      integer       :: max_ecri
      real (kind=Rkind) :: max_norm,phase

!----- for the field --------------------------------------------------
      real (kind=Rkind), allocatable :: tab_dnE(:,:)
      real (kind=Rkind)    :: dnE(3)
      real (kind=Rkind)    :: ww(3)
      logical :: make_field
!----- for the field --------------------------------------------------


!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='march_split_field'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Tmax,deltaT',para_propa%WPTmax,para_propa%WPdeltaT
        write(out_unitp,*) 'Hmin,Hmax',para_propa%para_poly%Hmin,               &
                                para_propa%para_poly%Hmax
        write(out_unitp,*)
        write(out_unitp,*) 'nb_ba,nb_qa',WP(1)%nb_ba,WP(1)%nb_qa
        write(out_unitp,*) 'nb_bi',WP(1)%nb_bi
        write(out_unitp,*)

        DO i=1,nb_WP
          CALL norm2_psi(WP(i),.FALSE.,.TRUE.,.FALSE.)
          write(out_unitp,*) 'norm WP',i,WP(i)%norme

          write(out_unitp,*) 'WP BasisRep',i
          CALL ecri_psi(T=ZERO,psi=WP(i),                               &
                        ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)
          write(out_unitp,*) 'WP GridRep',i
          CALL ecri_psi(T=ZERO,psi=WP(i),                               &
                        ecri_GridRep=.TRUE.,ecri_BasisRep=.FALSE.)
        END DO
      END IF
!-----------------------------------------------------------

!     psi(t+dt) = exp(-i.V/2.Dt).exp(-iT.Dt).exp(-i.V/2.Dt).psi(t)


!-----------------------------------------------------------
      w1 => para_propa%work_WP(para_propa%para_poly%npoly+1)
      w2 => para_propa%work_WP(para_propa%para_poly%npoly+2)
      w3 => para_propa%work_WP(para_propa%para_poly%npoly+3)
      expT => para_propa%work_WP(para_propa%para_poly%npoly+4)
      expV => para_propa%work_WP(para_propa%para_poly%npoly+5)

      w1 = WP(1)
      w1 = ONE
!     expT = exp(-idt/2 T)
      CALL sub_OpiPsi(w1,expT,para_H,3)
      expT = (-EYE*para_propa%WPdeltaT*HALF)*expT
      expT%CvecB(:)=exp(expT%CvecB(:))

      ip=3
      CALL sub_PsiBasisRep_TO_GridRep(w1)
      w1%CvecG(:)=ONE
      CALL sub_PsiGridRep_TO_BasisRep(w1)
      CALL sub_OpiPsi(w1,w2,para_H,1)
      CALL sub_OpPsi(w1,w3,para_Dip(ip))
      w1 = w2-dnE(ip)*w3
      expV = (-EYE*para_propa%WPdeltaT)*w1
      CALL sub_PsiBasisRep_TO_GridRep(expV)
      expV%CvecG(:)=exp(expV%CvecG(:))


      DO j=1,nb_WP

!       w3 = exp(-idt/2 T).WP
        w3%CvecB(:) = expT%CvecB(:) * WP(j)%CvecB(:)

!       w3 = exp(-idt (V-E.mu)).w3
        CALL sub_PsiBasisRep_TO_GridRep(w3)
        w3%CvecG(:) = expV%CvecG(:) * w3%CvecG(:)
        CALL sub_PsiGridRep_TO_BasisRep(w3)

!       WP(j) = exp(-idt/2 T).w3
        WP(j)%CvecB(:) = expT%CvecB(:) * w3%CvecB(:)

!       - Phase Shift -----------------
        phase = para_H%E0*para_propa%WPdeltaT
        WP(j) = WP(j) * exp(-EYE*phase)

!       - check norm ------------------
        CALL norm2_psi(WP(j),GridRep=.FALSE.,BasisRep=.TRUE.)
        IF ( WP(j)%norme > WP(j)%max_norme) THEN
          T  = T + para_propa%WPdeltaT
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' STOP propagation: norm > max_norm',              &
                         WP(j)%norme
          para_propa%march_error   = .TRUE.
          para_propa%test_max_norm = .TRUE.
          STOP
        END IF

      END DO

!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------


      END SUBROUTINE march_split_field

      END MODULE mod_march

