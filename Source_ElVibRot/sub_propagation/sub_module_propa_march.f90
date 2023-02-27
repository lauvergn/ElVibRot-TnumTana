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
 MODULE mod_march
 USE mod_system
 use mod_Constant,      ONLY : get_conv_au_to_unit
 USE mod_psi,           ONLY : param_psi,alloc_NParray,dealloc_NParray,dealloc_psi
 USE mod_propa,         ONLY : param_propa,param_poly,Calc_AutoCorr,    &
                               Write_AutoCorr,SaveWP_restart,           &
                               sub_analyze_mini_WP_OpWP
 USE mod_field,         ONLY : param_field
 USE mod_march_MPI,     ONLY : march_SIL_MPI
 USE mod_march_SG4
 IMPLICIT NONE

 PRIVATE
 PUBLIC :: march_gene,march_nOD_im,march_FOD_Opti_im,initialisation1_poly

      CONTAINS
!=======================================================================================
!    march_gene  with or without field
!=======================================================================================
      SUBROUTINE march_gene(T,WP,WP0,nb_WP,print_Op,                    &
                            para_H,para_propa,tab_Op,para_field)
      USE mod_system
      USE mod_psi,   ONLY : param_psi,ecri_psi,norm2_psi
      USE mod_Op,    ONLY : param_Op
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
        IF(.NOT. openmpi) write(out_unitp,*) 'Hmin,Hmax',para_propa%para_poly%Hmin,    &
                                                         para_propa%para_poly%Hmax
        write(out_unitp,*)
        write(out_unitp,*) 'nb_ba,nb_qa',WP%nb_ba,WP%nb_qa
        write(out_unitp,*) 'nb_bi',WP%nb_bi
        write(out_unitp,*)
        write(out_unitp,*) 'type_WPpropa',para_propa%type_WPpropa
        write(out_unitp,*)

        DO i=1,nb_WP
          CALL norm2_psi(WP(i),.FALSE.,.TRUE.,.FALSE.)
          write(out_unitp,*) 'norm WP',i,WP(i)%norm2

          write(out_unitp,*) 'WP BasisRep',i
          CALL ecri_psi(T=ZERO,psi=WP(i),                               &
                        ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)
          write(out_unitp,*) 'WP GridRep',i
          CALL ecri_psi(T=ZERO,psi=WP(i),                               &
                        ecri_GridRep=.TRUE.,ecri_BasisRep=.FALSE.)
        END DO
      END IF
!-----------------------------------------------------------

      IF (para_propa%One_Iteration) THEN
        IF(MPI_id==0) THEN
          write(out_unitp,*) '================================================='
          write(out_unitp,*) ' VIB: BEGINNING ',name_sub
          CALL time_perso('march_gene')
          write(out_unitp,*)
        ENDIF
      END IF

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

      ! para_propa%type_WPpropa 1
      ! SGtype4 F
      SELECT CASE (para_propa%type_WPpropa)

      CASE (1) ! Chebyshev
        !CALL march_cheby_old(T,no,WP(1),WP0(1),para_H,para_propa)
        CALL march_cheby(T,no,WP(1),WP0(1),para_H,para_propa)

      CASE (2) ! n-Order Taylor expansion
        !IF (SGtype4 .AND. direct_KEO) THEN
        IF (SGtype4) THEN
          CALL  march_noD_SG4(T,no,WP(1),WP0(1),para_H,para_propa)
        ELSE
          CALL  march_noD(T,no,WP(1),WP0(1),para_H,para_propa)
        END IF

      CASE (5) ! Runge-Kunta
        IF (para_propa%para_poly%npoly == 2) THEN
          CALL march_RK2(T,no,WP(1),WP0(1),para_H,para_propa)
        else
          CALL march_RK4(T,no,WP(1),WP0(1),para_H,para_propa)
        END IF

      CASE (6) !
        CALL march_ModMidPoint(T,no,WP(1),WP0(1),para_H,para_propa)

      CASE (7)
        CALL march_BS(T,no,WP(1),WP0(1),para_H,para_propa)

      CASE (8) ! Short Iterative Lanczos
        IF(openmpi .AND. SRep_MPI) THEN
          CALL march_SIL_MPI(T,no,WP(1),WP0(1),para_H,para_propa)
        ELSE
          CALL march_SIL(T,no,WP(1),WP0(1),para_H,para_propa)
        ENDIF

      CASE (9) ! Short Iterative Propagation (Lanczos with complete orthogonalisation)
        !CALL march_SIP(T,no,WP(1),WP0(1),para_H,para_propa)
        !CALL march_SIP2(T,no,WP(1),WP0(1),para_H,para_propa)
        IF (SGtype4) THEN
          !CALL march_SIP2_Srep(T,no,WP(1),WP0(1),para_H,para_propa)
          CALL march_SIP2_Srep2(T,no,WP(1),WP0(1),para_H,para_propa)
        ELSE
          CALL march_SIP(T,no,WP(1),WP0(1),para_H,para_propa)
          !CALL march_SIP2(T,no,WP(1),WP0(1),para_H,para_propa)
        END IF

      CASE (10) ! Spectral (bug???)
        CALL march_Spectral(T,no,WP(1),WP0(1),para_H,para_propa)

      CASE (22,24) ! n-Order Taylor expansion (with Field)
!       CALL march_split_field(T,WP(:),nb_WP,print_Op,
!    *                       para_H,tab_Op(3:5),
!    *                       para_field,para_propa)
!       CALL march_new_noD_field(T,WP(:),nb_WP,print_Op,
!    *                       para_H,tab_Op(3:5),
!    *                       para_field,para_propa)
        CALL march_noD_field(T,WP(:),nb_WP,print_Op,                    &
                             para_H,tab_Op(1:3),                        &
                             para_field,para_propa)

      CASE (50,54) ! 4th-Order Runge-Kunta (with Field)
        CALL march_RK4_field(T,WP(:),nb_WP,print_Op,                    &
                             para_H,tab_Op(1:3),                        &
                             para_field,para_propa)

      CASE (52)  ! 2d-Order Runge-Kunta (with Field)
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

      IF (para_propa%One_Iteration) THEN
        CALL sub_analyze_mini_WP_OpWP(T+para_propa%WPdeltaT,WP,1,para_H,        &
                                                         para_propa%ana_psi)

        IF(MPI_id==0) THEN
          write(out_unitp,*)
          CALL time_perso('march_gene')
          write(out_unitp,*)
          write(out_unitp,*) ' VIB: END march_gene'
          write(out_unitp,*) '================================================='
          write(out_unitp,*) 'Propagation stops after one time step.'
          write(out_unitp,*) ' ElVibRot-Tnum AU REVOIR!!!'
          write(out_unitp,*) '================================================'
        ENDIF
        STOP 'Propagation stops after one time step.'
      END IF

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
!=======================================================================================

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
      USE mod_psi,   ONLY : param_psi,ecri_psi,norm2_psi
      USE mod_Op,    ONLY : param_Op
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H
      TYPE (param_Op)   :: para_Dip(3)

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_field) :: para_field

      integer          :: nb_WP
      TYPE (param_psi) :: WP(:)

!----- for printing --------------------------------------------------
      logical :: print_Op




!------ working parameters --------------------------------
      TYPE (param_psi)  :: w1,w2,w3,w4,w5,w6
      integer           :: i
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
          write(out_unitp,*) 'norm WP',i,WP(i)%norm2

          write(out_unitp,*) 'WP BasisRep',i
          CALL ecri_psi(T=ZERO,psi=WP(i),                              &
                        ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)
          write(out_unitp,*) 'WP GridRep',i
          CALL ecri_psi(T=ZERO,psi=WP(i),                              &
                        ecri_GridRep=.TRUE.,ecri_BasisRep=.FALSE.)
        END DO
       END IF
!-----------------------------------------------------------

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
      IF ( WP(i)%norm2 > para_propa%max_norm2) THEN
        T  = T + para_propa%WPdeltaT
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' STOP propagation: norm2 > max_norm2',WP%norm2
        para_propa%march_error   = .TRUE.
        para_propa%test_max_norm = .TRUE.
        STOP
      END IF

      END DO


       CALL dealloc_psi(w1)
       CALL dealloc_psi(w2)
       CALL dealloc_psi(w3)
       CALL dealloc_psi(w4)
       CALL dealloc_psi(w5)
       CALL dealloc_psi(w6)

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
      USE mod_psi,   ONLY : param_psi,ecri_psi,norm2_psi
      USE mod_Op,    ONLY : param_Op
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H
      TYPE (param_Op)   :: para_Dip(3)

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_field) :: para_field

      integer          :: nb_WP
      TYPE (param_psi) :: WP(:)

!----- for printing --------------------------------------------------
      logical ::print_Op



!------ working parameters --------------------------------
      TYPE (param_psi)  :: w1,w2,w3,w4,w5,w6
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
          write(out_unitp,*) 'norm WP',i,WP(i)%norm2

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
      IF ( WP(i)%norm2 > para_propa%max_norm2) THEN
        T  = T + para_propa%WPdeltaT
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' STOP propagation: norm2 > max_norm2',WP%norm2
        para_propa%march_error   = .TRUE.
        para_propa%test_max_norm = .TRUE.
        STOP
      END IF

      END DO

       CALL dealloc_psi(w1)
       CALL dealloc_psi(w2)
       CALL dealloc_psi(w3)
       CALL dealloc_psi(w4)
       CALL dealloc_psi(w5)
       CALL dealloc_psi(w6)

!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------

      END SUBROUTINE march_RK4_field
      SUBROUTINE march_Euler(T,no,WP,WP0,para_H,para_propa)
      USE mod_system
      USE mod_psi,   ONLY : param_psi,ecri_psi,norm2_psi
      USE mod_Op,    ONLY : param_Op
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa

      TYPE (param_psi) :: WP,WP0

      integer :: no
!----- for printing --------------------------------------------------
      logical ::print_Op



!------ working parameters --------------------------------
      TYPE (param_psi)     :: w1,w2,w6
      complex (kind=Rkind) :: cdot
      integer              :: i,ip
      real (kind=Rkind)    :: T      ! time
      real (kind=Rkind)    :: T_DT
      real (kind=Rkind)    :: phase

      integer  ::   nioWP

!----- for the field --------------------------------------------------
      real (kind=Rkind)    :: dnE(3)
      real (kind=Rkind)    :: ww(3)
      logical :: make_field
!----- for the field --------------------------------------------------


!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='march_Euler'
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
          write(out_unitp,*) 'norm WP',i,WP%norm2

          write(out_unitp,*) 'WP BasisRep'
          CALL ecri_psi(T=ZERO,psi=WP,                               &
                        ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)
          write(out_unitp,*) 'WP GridRep'
          CALL ecri_psi(T=ZERO,psi=WP,                               &
                        ecri_GridRep=.TRUE.,ecri_BasisRep=.FALSE.)

      END IF
!-----------------------------------------------------------

       !-----------------------------------------------------------
       !w1 = WPdeltaT * fcn(WP,T)

       CALL fcn(WP,w1,para_H)
       w1 = w1 * para_propa%WPdeltaT
       WP = WP + w1

       !- Phase Shift -----------------
       phase = para_H%E0*para_propa%WPdeltaT
       WP = WP * exp(-cmplx(ZERO,phase,kind=Rkind))

      !- check norm ------------------
      CALL norm2_psi(WP,GridRep=.FALSE.,BasisRep=.TRUE.)
      IF ( WP%norm2 > para_propa%max_norm2) THEN
        T  = T + para_propa%WPdeltaT
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' STOP propagation: norm2 > max_norm2',WP%norm2
        para_propa%march_error   = .TRUE.
        para_propa%test_max_norm = .TRUE.
        STOP
      END IF

      cdot = Calc_AutoCorr(WP0,WP,para_propa,T,Write_AC=.FALSE.)
      CALL Write_AutoCorr(no,T + para_propa%WPdeltaT,cdot)
      flush(no)

      CALL dealloc_psi(w1)
      CALL dealloc_psi(w2)
      CALL dealloc_psi(w6)
!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------

END SUBROUTINE march_Euler
      SUBROUTINE march_RK2(T,no,WP,WP0,para_H,para_propa)
      USE mod_system
      USE mod_psi,   ONLY : param_psi,ecri_psi,norm2_psi
      USE mod_Op,    ONLY : param_Op
      USE mod_MPI_aux
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa

      TYPE (param_psi) :: WP,WP0

      integer :: no
!----- for printing --------------------------------------------------
      logical ::print_Op



!------ working parameters --------------------------------
      TYPE (param_psi)     :: w1,w2,w6
      complex (kind=Rkind) :: cdot
      integer              :: i,ip
      real (kind=Rkind)    :: T      ! time
      real (kind=Rkind)    :: T_DT
      real (kind=Rkind)    :: phase

      integer  ::   nioWP

!----- for the field --------------------------------------------------
      real (kind=Rkind)    :: dnE(3)
      real (kind=Rkind)    :: ww(3)
      logical :: make_field
!----- for the field --------------------------------------------------


!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='march_RK2'
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
          write(out_unitp,*) 'norm WP',i,WP%norm2

          write(out_unitp,*) 'WP BasisRep'
          CALL ecri_psi(T=ZERO,psi=WP,                               &
                        ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)
          write(out_unitp,*) 'WP GridRep'
          CALL ecri_psi(T=ZERO,psi=WP,                               &
                        ecri_GridRep=.TRUE.,ecri_BasisRep=.FALSE.)

      END IF
!-----------------------------------------------------------

      !-----------------------------------------------------------
       !w1 = WPdeltaT * fcn(WP,T)

      CALL fcn(WP,w1,para_H)
      IF(keep_MPI) w1 = w1 * para_propa%WPdeltaT

      !-----------------------------------------------------------
      !w6 = WP + w1
      !w2 = WPdeltaT * fcn(WP+w1,T+DT)
      T_DT = T + para_propa%WPdeltaT
      IF(keep_MPI) w6 = WP + w1
      CALL fcn(w6,w2,para_H)

      IF(keep_MPI) w2 = w2 * para_propa%WPdeltaT

      IF(keep_MPI) THEN
        w6     = w1+w2
        WP     = WP + w6*HALF

        !- Phase Shift -----------------
        phase = para_H%E0*para_propa%WPdeltaT
        WP = WP * exp(-cmplx(ZERO,phase,kind=Rkind))
      ENDIF

      !- check norm ------------------
      IF(keep_MPI) CALL norm2_psi(WP,GridRep=.FALSE.,BasisRep=.TRUE.)
      IF(openmpi)  CALL MPI_Bcast_(WP%norm2,size1_MPI,root_MPI)

      IF ( WP%norm2 > para_propa%max_norm2) THEN
        T  = T + para_propa%WPdeltaT
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' STOP propagation: norm2 > max_norm2',WP%norm2
        para_propa%march_error   = .TRUE.
        para_propa%test_max_norm = .TRUE.
        STOP
      END IF

      IF(MPI_id==0) THEN
        cdot = Calc_AutoCorr(WP0,WP,para_propa,T,Write_AC=.FALSE.)
        CALL Write_AutoCorr(no,T + para_propa%WPdeltaT,cdot)
        flush(no)
      ENDIF

      IF(keep_MPI) THEN
        CALL dealloc_psi(w1)
        CALL dealloc_psi(w2)
        CALL dealloc_psi(w6)
      ENDIF
!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------

END SUBROUTINE march_RK2
SUBROUTINE march_RK4(T,no,WP,WP0,para_H,para_propa)
      USE mod_system
      USE mod_psi,   ONLY : param_psi,ecri_psi,norm2_psi
      USE mod_Op,    ONLY : param_Op
      USE mod_MPI_aux
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa

      TYPE (param_psi) :: WP,WP0

      integer :: no
!----- for printing --------------------------------------------------
      logical ::print_Op



!------ working parameters --------------------------------
      TYPE (param_psi)     :: w1,w2,w3,w4,w
      complex (kind=Rkind) :: cdot
      integer              :: i,ip
      real (kind=Rkind)    :: T      ! time
      real (kind=Rkind)    :: T_DT,T_DT_HALF
      real (kind=Rkind)    :: DT,DTo2,DTo6
      real (kind=Rkind)    :: phase

      integer  ::   nioWP

!----- for the field --------------------------------------------------
      real (kind=Rkind)    :: dnE(3)
      real (kind=Rkind)    :: ww(3)
      logical :: make_field
!----- for the field --------------------------------------------------


!----- for debuging --------------------------------------------------
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
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
          write(out_unitp,*) 'norm WP',i,WP%norm2

          write(out_unitp,*) 'WP BasisRep',T
          CALL ecri_psi(T=T,psi=WP,                               &
                        ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)
          write(out_unitp,*) 'WP GridRep',T
          CALL ecri_psi(T=T,psi=WP,                               &
                        ecri_GridRep=.TRUE.,ecri_BasisRep=.FALSE.)

      END IF
!-----------------------------------------------------------
      T_DT      = T + para_propa%WPdeltaT
      T_DT_HALF = T + para_propa%WPdeltaT*HALF
      DT        = para_propa%WPdeltaT
      DTo2      = para_propa%WPdeltaT/TWO
      DTo6      = para_propa%WPdeltaT/SIX


!-----------------------------------------------------------

      !w1 = -i.H.WP
      CALL fcn(WP,w1,para_H)

      IF (debug) THEN
        write(out_unitp,*) 'k1 BasisRep'
        CALL ecri_psi(T=T,psi=w1)
      END IF

      !w = WP + w1 * (DT/2)
      IF(keep_MPI) w = WP + w1 * DTo2

      !w2 = -iH*w
      CALL fcn(w,w2,para_H)
      IF (debug) THEN
        write(out_unitp,*) 'k2 BasisRep'
        CALL ecri_psi(T=T,psi=w2)
      END IF


      !w = WP + w2 * (DT/2)
      IF(keep_MPI) w = WP + w2 * DTo2

      !w3 = -iH*w2
      CALL fcn(w,w3,para_H)
      IF (debug) THEN
        write(out_unitp,*) 'k3 BasisRep'
        CALL ecri_psi(T=T,psi=w3)
      END IF



      !w = WP + w3 * DT
      IF(keep_MPI) w = WP + w3 * DT


      !w4 = -iH(T)*w
      CALL fcn(w,w4,para_H)
      IF (debug) THEN
        write(out_unitp,*) 'k4 BasisRep'
        CALL ecri_psi(T=T,psi=w4)
      END IF


      !WP = WP + (DT/6)*(w1 + 2*w2 + 2*w3 + w4)
      IF(keep_MPI) THEN
        WP = WP + (w1+w4)*DTo6+(w2+w3)*(TWO*DTo6)
      ENDIF

      !- Phase Shift -----------------
      phase = para_H%E0*para_propa%WPdeltaT
      IF(keep_MPI) WP = WP * exp(-cmplx(ZERO,phase,kind=Rkind))

      !- check norm ------------------
      IF(keep_MPI) CALL norm2_psi(WP,GridRep=.FALSE.,BasisRep=.TRUE.)
      IF(openmpi)  CALL MPI_Bcast_(WP%norm2,size1_MPI,root_MPI)

      IF ( WP%norm2 > para_propa%max_norm2) THEN
        T  = T + para_propa%WPdeltaT
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' STOP propagation: norm2 > max_norm2',WP%norm2
        para_propa%march_error   = .TRUE.
        para_propa%test_max_norm = .TRUE.
        STOP
      END IF

      IF(MPI_id==0) THEN
        cdot = Calc_AutoCorr(WP0,WP,para_propa,T,Write_AC=.FALSE.)
        CALL Write_AutoCorr(no,T + para_propa%WPdeltaT,cdot)
        flush(no)
      ENDIF

      IF(keep_MPI) THEN
        CALL dealloc_psi(w1)
        CALL dealloc_psi(w2)
        CALL dealloc_psi(w3)
        CALL dealloc_psi(w4)
        CALL dealloc_psi(w)
      ENDIF
!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'WP BasisRep',T_DT
        CALL ecri_psi(T=T_DT,psi=WP,                               &
                      ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)
        write(out_unitp,*) 'WP GridRep',T_DT
        CALL ecri_psi(T=T_DT,psi=WP,                               &
                      ecri_GridRep=.TRUE.,ecri_BasisRep=.FALSE.)

        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------

END SUBROUTINE march_RK4
SUBROUTINE march_RK4_old(T,no,WP,WP0,para_H,para_propa)
      USE mod_system
      USE mod_psi,   ONLY : param_psi,ecri_psi,norm2_psi
      USE mod_Op,    ONLY : param_Op
      USE mod_MPI_aux
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa

      TYPE (param_psi) :: WP,WP0

      integer :: no
!----- for printing --------------------------------------------------
      logical ::print_Op



!------ working parameters --------------------------------
      TYPE (param_psi)     :: w1,w2,w3,w4,w5,w6
      complex (kind=Rkind) :: cdot
      integer              :: i,ip
      real (kind=Rkind)    :: T      ! time
      real (kind=Rkind)    :: T_DT,T_DT_HALF
      real (kind=Rkind)    :: DT,DTo2,DTo6
      real (kind=Rkind)    :: phase

      integer  ::   nioWP

!----- for the field --------------------------------------------------
      real (kind=Rkind)    :: dnE(3)
      real (kind=Rkind)    :: ww(3)
      logical :: make_field
!----- for the field --------------------------------------------------


!----- for debuging --------------------------------------------------
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='march_RK4_old'
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
          write(out_unitp,*) 'norm WP',i,WP%norm2

          write(out_unitp,*) 'WP BasisRep',T
          CALL ecri_psi(T=T,psi=WP,                               &
                        ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)
          write(out_unitp,*) 'WP GridRep',T
          CALL ecri_psi(T=T,psi=WP,                               &
                        ecri_GridRep=.TRUE.,ecri_BasisRep=.FALSE.)

      END IF
!-----------------------------------------------------------
      T_DT      = T + para_propa%WPdeltaT
      T_DT_HALF = T + para_propa%WPdeltaT*HALF
      DT        = para_propa%WPdeltaT
      DTo2      = para_propa%WPdeltaT/TWO
      DTo6      = para_propa%WPdeltaT/SIX


!-----------------------------------------------------------

      !w1 = -i.H.WP
      CALL fcn(WP,w1,para_H)

      IF (debug) THEN
        write(out_unitp,*) 'k1 BasisRep'
        CALL ecri_psi(T=T,psi=WP)
      END IF

      !w2 = WP + w1 * (DT/2)
      IF(keep_MPI) w2 = WP + w1 * DTo2

      !w3 = -iH*w2
      CALL fcn(w2,w3,para_H)
      IF (debug) THEN
        write(out_unitp,*) 'k2 BasisRep'
        CALL ecri_psi(T=T,psi=WP)
      END IF


      !w2 = WP + w3 * (DT/2)
      IF(keep_MPI) w2 = WP + w3 * DTo2

      !w4 = -iH*w2
      CALL fcn(w2,w4,para_H)
      IF (debug) THEN
        write(out_unitp,*) 'k3 BasisRep'
        CALL ecri_psi(T=T,psi=WP)
      END IF



      !w2 = WP + w4 * DT
      IF(keep_MPI) w2 = WP + w4 * DT

      !w4 = w3 + w4
      IF(keep_MPI) w4 = w3 + w4

      !w3 = -iH(T)*w2
      CALL fcn(w2,w3,para_H)
      IF (debug) THEN
        write(out_unitp,*) 'k4 BasisRep'
        CALL ecri_psi(T=T,psi=WP)
      END IF


      !WP = WP + (DT/6)*(w1+w3+2.D0*w4)
      IF(keep_MPI) THEN
        w5 = w1 + w3
        w5 = w5 + w4*TWO
        WP = WP + w5*DTo6
      ENDIF

      !- Phase Shift -----------------
      phase = para_H%E0*para_propa%WPdeltaT
      IF(keep_MPI) WP = WP * exp(-cmplx(ZERO,phase,kind=Rkind))

      !- check norm ------------------
      IF(keep_MPI) CALL norm2_psi(WP,GridRep=.FALSE.,BasisRep=.TRUE.)
      IF(openmpi)  CALL MPI_Bcast_(WP%norm2,size1_MPI,root_MPI)

      IF ( WP%norm2 > para_propa%max_norm2) THEN
        T  = T + para_propa%WPdeltaT
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' STOP propagation: norm2 > max_norm2',WP%norm2
        para_propa%march_error   = .TRUE.
        para_propa%test_max_norm = .TRUE.
        STOP
      END IF

      IF(MPI_id==0) THEN
        cdot = Calc_AutoCorr(WP0,WP,para_propa,T,Write_AC=.FALSE.)
        CALL Write_AutoCorr(no,T + para_propa%WPdeltaT,cdot)
        flush(no)
      ENDIF

      IF(keep_MPI) THEN
        CALL dealloc_psi(w1)
        CALL dealloc_psi(w2)
        CALL dealloc_psi(w3)
        CALL dealloc_psi(w4)
        CALL dealloc_psi(w5)
        CALL dealloc_psi(w6)
      ENDIF
!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'WP BasisRep',T_DT
        CALL ecri_psi(T=T_DT,psi=WP,                               &
                      ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)
        write(out_unitp,*) 'WP GridRep',T_DT
        CALL ecri_psi(T=T_DT,psi=WP,                               &
                      ecri_GridRep=.TRUE.,ecri_BasisRep=.FALSE.)

        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------

END SUBROUTINE march_RK4_old
      SUBROUTINE march_BS(T,no,WP,WP0,para_H,para_propa)
      USE mod_system
      USE mod_psi,   ONLY : param_psi,ecri_psi,norm2_psi,dealloc_psi
      USE mod_Op,    ONLY : param_Op
      USE mod_MPI_aux
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
          write(out_unitp,*) 'norm WP',i,WP%norm2

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

  IF(keep_MPI) THEN
    yerr = WP
    yerr = ZERO
  ENDIF

  allocate(yt0(0:0))
  m = 2
  IF(keep_MPI) yt0(0) = WP

  CALL march_ModMidPoint(T,no,yt0(0),WP0,para_H,para_propa,m)
  err0 = huge(ONE)

  DO j=1,order
    allocate(yt1(0:j))
    m = 2*j+2
    IF(keep_MPI) yt1(0) = WP
    CALL march_ModMidPoint(T,no,yt1(0),WP0,para_H,para_propa,m)

    !extrapolation
    DO i=1,j
      x = real(2*j+2,kind=Rkind)/real(2*(j-i)+2,kind=Rkind)
      x = ONE/(x**2-ONE)
      IF(keep_MPI) yerr = yt1(i-1)-yt0(i-1)
      IF(keep_MPI) yt1(i) = yt1(i-1) + yerr*x
    END DO

    IF(keep_MPI) yerr = yt1(j) - yt0(j-1)
    IF(keep_MPI) CALL norm2_psi(yerr)
    IF(openmpi) CALL MPI_Bcast_(yerr%norm2,size1_MPI,root_MPI)
    err1 = sqrt(yerr%norm2)

    DO i=lbound(yt0,dim=1),ubound(yt0,dim=1)
      IF(keep_MPI) CALL dealloc_psi(yt0(i),delete_all=.TRUE.)
    END DO

    deallocate(yt0)
    allocate(yt0(0:j))
    DO i=lbound(yt1,dim=1),ubound(yt1,dim=1)
      IF(keep_MPI) yt0(i) = yt1(i)
      IF(keep_MPI) CALL dealloc_psi(yt1(i),delete_all=.TRUE.)
    END DO
    deallocate(yt1)
    !write(out_unitp,*) 'march_bs',j,err1

    !IF (err1 < 1.d-6 .OR. err1 > err0) EXIT
    IF (err1 < 1.d-10) EXIT

    err0 = err1
  END DO
  write(out_unitp,*) 'end march_bs',min(j,order),err1

  IF(keep_MPI) WP = yt0(min(j,order))
  DO i=lbound(yt0,dim=1),ubound(yt0,dim=1)
    IF(keep_MPI) CALL dealloc_psi(yt0(i),delete_all=.TRUE.)
  END DO
  deallocate(yt0)
  IF(keep_MPI) CALL dealloc_psi(yerr,delete_all=.TRUE.)

      !- Phase Shift -----------------
      phase = para_H%E0*para_propa%WPdeltaT
      IF(keep_MPI) WP = WP * exp(-cmplx(ZERO,phase,kind=Rkind))

      !- check norm ------------------
      IF(keep_MPI) CALL norm2_psi(WP,GridRep=.FALSE.,BasisRep=.TRUE.)
      IF(openmpi) CALL MPI_Bcast_(WP%norm2,size1_MPI,root_MPI)
      IF ( WP%norm2 > para_propa%max_norm2) THEN
        T  = T + para_propa%WPdeltaT
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' STOP propagation: norm2 > max_norm2',WP%norm2
        para_propa%march_error   = .TRUE.
        para_propa%test_max_norm = .TRUE.
        STOP
      END IF

      IF(MPI_id==0) THEN
        cdot = Calc_AutoCorr(WP0,WP,para_propa,T,Write_AC=.FALSE.)
        CALL Write_AutoCorr(no,T + para_propa%WPdeltaT,cdot)
        flush(no)
      ENDIF

!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------

      END SUBROUTINE march_BS
      SUBROUTINE march_ModMidPoint(T,no,WP,WP0,para_H,para_propa,order)
      USE mod_system
      USE mod_psi,   ONLY : param_psi,ecri_psi,norm2_psi,dealloc_psi
      USE mod_Op,    ONLY : param_Op
      USE mod_MPI_aux
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
          write(out_unitp,*) 'norm WP',i,WP%norm2

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
!write(out_unitp,*) 'order_loc',order_loc
      DTT        = para_propa%WPdeltaT / real(order_loc,kind=Rkind)

      IF(keep_MPI) zkm = WP
      CALL fcn(WP,dWP,para_H)
      IF(keep_MPI) zk  = WP + dWP * DTT


      DO k=1,order_loc-1

        CALL fcn(zk,dWP,para_H)
        IF(keep_MPI) THEN
          zkp = zkm + TWO*DTT * dWP
          zkm = zk
          zk  = zkp
        ENDIF
      END DO

      CALL fcn(zk,zkp,para_H)
      IF(keep_MPI) THEN
        dWP = DTT * zkp

        zkp = zk + zkm

        zk =  zkp + dWP
        WP = HALF * zk
      ENDIF

      IF(keep_MPI) THEN
        CALL dealloc_psi(dWP,delete_all=.TRUE.)
        CALL dealloc_psi(zkm,delete_all=.TRUE.)
        CALL dealloc_psi(zk ,delete_all=.TRUE.)
        CALL dealloc_psi(zkp,delete_all=.TRUE.)
      ENDIF

      IF (present(order)) RETURN

      !- Phase Shift -----------------
      phase = para_H%E0*para_propa%WPdeltaT
      IF(keep_MPI) WP = WP * exp(-cmplx(ZERO,phase,kind=Rkind))

      !- check norm ------------------
      IF(keep_MPI) CALL norm2_psi(WP,GridRep=.FALSE.,BasisRep=.TRUE.)
      IF(openmpi)  CALL MPI_Bcast_(WP%norm2,size1_MPI,root_MPI)

      IF ( WP%norm2 > para_propa%max_norm2) THEN
        T  = T + para_propa%WPdeltaT
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' STOP propagation: norm2 > max_norm2',WP%norm2
        para_propa%march_error   = .TRUE.
        para_propa%test_max_norm = .TRUE.
        STOP
      END IF

      IF(MPI_id==0) THEN
        cdot = Calc_AutoCorr(WP0,WP,para_propa,T,Write_AC=.FALSE.)
        CALL Write_AutoCorr(no,T + para_propa%WPdeltaT,cdot)
        flush(no)
      ENDIF

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
      USE mod_psi,   ONLY : param_psi,ecri_psi,norm2_psi
      USE mod_Op,    ONLY : param_Op,sub_OpPsi,sub_scaledOpPsi
      USE mod_field, ONLY : param_field,sub_dnE
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
        write(out_unitp,*) 'norm WP BasisRep',WP%norm2

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

         !w1 = ZERO
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
      USE mod_basis, ONLY : get_nb_TDParam_FROM_basis

      USE mod_psi,   ONLY : param_psi,ecri_psi,norm2_psi
      USE mod_Op,    ONLY : param_Op, sub_OpPsi,sub_scaledOpPsi
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ------------------------------
      TYPE (param_psi) :: WP,dWP

!----- for printing --------------------------------------------------
      logical ::print_Op


      integer  ::   nioWP
      integer  ::   nb_TDParam


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
        write(out_unitp,*) 'norm WP BasisRep',WP%norm2

        write(out_unitp,*) 'WP'
        CALL ecri_psi(T=ZERO,psi=WP)
       END IF
!-----------------------------------------------------------

       !nb_TDParam = get_nb_TDParam_FROM_basis(WP%BasisnD)
       IF(keep_MPI) write(out_unitp,*) 'In fcn: nb_TDParam',WP%nb_TDParam

!-----------------------------------------------------------
!      dWP = -i.H.WP

       CALL sub_OpPsi(WP,dWP,para_H)
       CALL sub_scaledOpPsi(WP,dWP,para_H%E0,ONE)
       IF(keep_MPI) dWP = dWP * (-EYE)

       IF (debug) THEN
        write(out_unitp,*) 'dWP'
        CALL ecri_psi(T=ZERO,psi=dWP)
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------


      END SUBROUTINE fcn
      SUBROUTINE Make_SMatrix_WITH_TDParam(S,WP,para_H)
      USE mod_system
      USE mod_basis, ONLY : get_nb_TDParam_FROM_basis

      USE mod_psi,   ONLY : param_psi,ecri_psi,norm2_psi
      USE mod_Op,    ONLY : param_Op, sub_OpPsi,sub_scaledOpPsi
      IMPLICIT NONE

      complex(kind=Rkind),  INTENT(INOUT) :: S(:,:)
      TYPE (param_Op),      intent(in)    :: para_H
      TYPE (param_psi),     intent(in)    :: WP


      integer  ::   nb_TDParam,dim_S


!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='Make_SMatrix_WITH_TDParam'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_tot',WP%nb_tot
        write(out_unitp,*) 'nb_TDParam',WP%nb_TDParam
        write(out_unitp,*)
        !write(out_unitp,*) 'WP'
        !CALL ecri_psi(T=ZERO,psi=WP)
       END IF
!-----------------------------------------------------------

       dim_S = WP%nb_tot+WP%nb_TDParam
       IF (dim_S /= size(S(:,1))) STOP 'ERROR in Make_SMatrix_WITH_TDParam: wrong size.'

!-----------------------------------------------------------
       !CALL Cplx_mat_id(S, WP%nb_tot,dim_S)
       S(1:WP%nb_tot,1:WP%nb_tot) = Identity_Mat(WP%nb_tot)
       STOP 'IN Make_SMatrix_WITH_TDParam : a vérifer !!'

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'S with TDParam'
        CALL Write_Mat(S,out_unitp,6)
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------

END SUBROUTINE Make_SMatrix_WITH_TDParam

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
      USE mod_psi,   ONLY : param_psi,ecri_psi,norm2_psi
      USE mod_Op,    ONLY : param_Op,sub_OpPsi,sub_scaledOpPsi
      USE mod_field, ONLY : param_field,sub_dnE
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H
      TYPE (param_Op)   :: para_Dip(3)

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_field) :: para_field

      integer          :: nb_WP
      TYPE (param_psi) :: WP(:)

!----- for printing --------------------------------------------------
      logical ::print_Op



!------ working parameters --------------------------------
      TYPE (param_psi)              :: w1,w2
      TYPE (param_psi) ,allocatable :: work_WP(:)

      complex (kind=Rkind)   :: Overlap,cdot,rt,rti,fac_k,norm
      real (kind=Rkind)      :: phase,limit
      integer                :: it,it_max,no,max_der,nb_der
      integer                :: i,j,k,jt,ip,iq
      real (kind=Rkind)      :: T      ! time
      real (kind=Rkind)      :: T_Delta! time+deltaT
      integer                :: max_ecri

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
          write(out_unitp,*) 'norm WP',i,WP(i)%norm2

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
      max_der = 0
      nb_der = para_propa%para_poly%npoly-1
      IF (para_field%max_der >= 0)                                      &
           nb_der = min(para_field%max_der,para_propa%para_poly%npoly-1)

      nullify(tab_dnE)
      CALL alloc_array(tab_dnE,[nb_der,3],"tab_dnE",name_sub,[0,1])
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

      CALL alloc_NParray(work_WP,[para_propa%para_poly%npoly-1],        &
                        'work_WP',name_sub,[0])
      DO j=1,nb_WP
        work_WP(0) = WP(j)

        DO i=0,para_propa%para_poly%npoly-1

          CALL sub_OpPsi(work_WP(i),w2,para_H)
          CALL sub_scaledOpPsi(work_WP(i),w2,para_H%E0,ONE)

          !work_WP(i+1) = w2  ! H.psi(i)


          fac_k = cmplx(ZERO,-para_propa%WPdeltaT /                     &
                           real(i+1,kind=Rkind),kind=Rkind)


          work_WP(i+1) = w2 * fac_k

          DO ip=1,3
            IF (.NOT. para_field%pola_xyz(ip)) CYCLE
            IF (maxval(abs(tab_dnE(:,ip))) <                            &
                para_propa%para_poly%poly_tol) CYCLE
            fac_k = cmplx(ZERO,-para_propa%WPdeltaT /                   &
                                 real(i+1,kind=Rkind),kind=Rkind)

            w1 = work_WP(i) * (cmplx(tab_dnE(0,ip),kind=Rkind) * fac_k)

            DO k=1,min(i,max_der)
               fac_k = fac_k * cmplx(para_propa%WPdeltaT/real(k,kind=Rkind),kind=Rkind)
               rt = cmplx(tab_dnE(k,ip),kind=Rkind)*fac_k
               w1 = w1 + work_WP(i-k) * rt
            END DO

            CALL sub_OpPsi(w1,w2,para_Dip(ip))

            work_WP(i+1) = work_WP(i+1) + w2

          END DO
          WP(j) = WP(j) + work_WP(i+1)

          CALL norm2_psi(work_WP(i+1))
          !write(out_unitp,*) 'norm2 psi i+1',i+1,work_WP(i+1)%norm2
          IF (work_WP(i+1)%norm2 < para_propa%para_poly%poly_tol) EXIT
          IF (work_WP(i+1)%norm2 > TEN**15) THEN
             write(out_unitp,*) ' ERROR in ',name_sub
             write(out_unitp,*) ' The norm2 of dnpsi(i+1) is huge !',           &
                                          work_WP(i+1)%norm2
             write(out_unitp,*) ' iteration i:',i
             write(out_unitp,*) ' Reduce the time step, WPDeltaT:',             &
                             para_propa%WPdeltaT
             para_propa%march_error = .TRUE.
             STOP
          END IF
        END DO



        IF (print_Op) write(out_unitp,*) 'T,max_der,ipoly,norm2',T,max_der,i,&
                          work_WP(min(i+1,para_propa%para_poly%npoly))%norm2

        work_WP(0) = work_WP(min(i+1,para_propa%para_poly%npoly))

!       - Phase Shift -----------------
        phase = para_H%E0*para_propa%WPdeltaT
        WP(j) = WP(j) * exp(-cmplx(ZERO,phase,kind=Rkind))

!       - check norm ------------------
        CALL norm2_psi(WP(j),GridRep=.FALSE.,BasisRep=.TRUE.)
        IF (WP(j)%norm2 > para_propa%max_norm2) THEN
          T  = T + para_propa%WPdeltaT
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' STOP propagation: norm^2 > max_norm^2',WP(j)%norm2
          write(out_unitp,*) ' norm^2:     ',WP(j)%norm2
          write(out_unitp,*) ' max_norm^2: ',para_propa%max_norm2
          para_propa%march_error   = .TRUE.
          para_propa%test_max_norm = .TRUE.
          STOP
        END IF

      END DO

      CALL dealloc_array(tab_dnE,"tab_dnE","march_noD_field")
      CALL dealloc_psi(w1)
      CALL dealloc_psi(w2)
      CALL dealloc_NParray(work_WP,'work_WP',name_sub)

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
      USE mod_psi,   ONLY : param_psi,ecri_psi,norm2_psi
      USE mod_Op,    ONLY : param_Op,sub_PsiOpPsi,sub_OpPsi,sub_scaledOpPsi
      USE mod_MPI_aux
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_psi)   :: psi0,psi

!------ working variables ---------------------------------
      TYPE (param_psi)     :: w1,w2
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

      IF(keep_MPI) w2 = psi
      CALL sub_PsiOpPsi(E,psi,w2,para_H)
      IF(openmpi .AND. MPI_scheme/=1) CALL MPI_Bcast_(E,size1_MPI,root_MPI)
!     para_H%E0 = real(E,kind=Rkind)
!     write(out_unitp,*) 'E0,Esc',para_H%E0,para_H%Esc

      psi0Hkpsi0(:) = cmplx(ZERO,ZERO,kind=Rkind)

      IF(keep_MPI) psi0Hkpsi0(0) = Calc_AutoCorr(psi0,psi,para_propa,T,Write_AC=.FALSE.)

      IF(keep_MPI) w1 = psi
      rtj = cmplx(ONE,ZERO,kind=Rkind)

 21 format(a,100(x,e12.5))

      !write(out_unitp,21) 'Rw1',Real(w1%CvecB,kind=Rkind)
      !write(out_unitp,21) 'Iw1',AImag(w1%CvecB)

      DO j=1,para_propa%para_poly%npoly

        !write(out_unitp,21) 'Rw1',Real(w1%CvecB,kind=Rkind)
        !write(out_unitp,21) 'Iw1',AImag(w1%CvecB)

        IF (j > 1) CALL sub_OpPsi(w1,w2,para_H) ! already done in PsiHPsi
        !write(out_unitp,21) 'Rw2',Real(w2%CvecB,kind=Rkind)
        !write(out_unitp,21) 'Iw2',AImag(w2%CvecB)

        CALL sub_scaledOpPsi(w1,w2,para_H%E0,ONE)

        !write(out_unitp,21) 'Rw2',Real(w2%CvecB,kind=Rkind)
        !write(out_unitp,21) 'Iw2',AImag(w2%CvecB)

        IF(keep_MPI) w1 = w2

        IF(keep_MPI) psi0Hkpsi0(j) = Calc_AutoCorr(psi0,w1,para_propa,T,Write_AC=.FALSE.)

        !CALL Overlap_psi1_psi2(psi0Hkpsi0(j),psi0,w1)

        rtj = rtj *                                                     &
          cmplx(ZERO,-para_propa%WPdeltaT/real(j,kind=Rkind),kind=Rkind)
        IF(keep_MPI) w2 = w1 * rtj

        !write(out_unitp,21) 'Rw2*rtj',Real(w2%CvecB,kind=Rkind)
        !write(out_unitp,21) 'Iw2*rtj',AImag(w2%CvecB)

        IF(keep_MPI) psi = psi + w2

        !write(out_unitp,21) 'Rpsi',Real(psi%CvecB,kind=Rkind)
        !write(out_unitp,21) 'Ipsi',AImag(psi%CvecB)

        IF(keep_MPI) CALL norm2_psi(w2)
        IF(openmpi .AND. MPI_scheme/=1)  CALL MPI_Bcast_(w2%norm2,size1_MPI,root_MPI)

        IF (debug) write(out_unitp,*) 'j,norm2 w2',j,w2%norm2

        IF (w2%norm2 > TEN**15) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' Norm of the vector is TOO large (> 10^15)',j,    &
                                                w2%norm2
          write(out_unitp,*) ' => Reduce the time step !!'
          STOP
        END IF
        j_exit = j
        IF (w2%norm2 < para_propa%para_poly%poly_tol) EXIT

      END DO
      write(out_unitp,*) 'j_exit,norms',j_exit,abs(w2%norm2)

!     write(out_unitp,*) 'j,norm2',j,abs(w1%norm2*rtj)
      IF (abs(w2%norm2) > para_propa%para_poly%poly_tol) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' Norm of the last vector is TOO large',w2%norm2
        write(out_unitp,*) ' poly_tol: ',para_propa%para_poly%poly_tol
        write(out_unitp,*) ' => npoly or max_poly are TOO small',               &
                      para_propa%para_poly%npoly
        STOP
      END IF

 !write(out_unitp,*) ' Psi before phase shift '
 !CALL ecri_psi(psi=psi)
!     - Phase Shift -----------------------------------
      phase = para_H%E0*para_propa%WPdeltaT
      IF(keep_MPI) psi   = psi * exp(-cmplx(ZERO,phase,kind=Rkind))

 !write(out_unitp,*) ' Psi after phase shift '
 !CALL ecri_psi(psi=psi)

!    - check norm ------------------
      IF(keep_MPI) CALL norm2_psi(psi,GridRep=.FALSE.,BasisRep=.TRUE.)
      IF(openmpi .AND. MPI_scheme/=1)  CALL MPI_Bcast_(psi%norm2,size1_MPI,root_MPI)

      IF ( psi%norm2 > para_propa%max_norm2) THEN
        T  = T + para_propa%WPdeltaT
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' STOP propagation: norm > max_norm',psi%norm2
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

      IF(MPI_id==0) THEN
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
      ENDIF

      IF(keep_MPI) THEN
        CALL dealloc_psi(w1)
        CALL dealloc_psi(w2)
      ENDIF
!-----------------------------------------------------------
      IF (debug) THEN
        CALL norm2_psi(psi)
        write(out_unitp,*) 'norm psi',psi%norm2
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

 END SUBROUTINE march_noD
!================================================================
!
!     march Short Iteration Progration (Lanczos ...)
!
!================================================================
      !!@description: march SIP
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
  SUBROUTINE march_SIP(T,no,psi,psi0,para_H,para_propa)
      USE mod_system
      USE mod_psi,   ONLY : param_psi,ecri_psi,norm2_psi,renorm_psi,Overlap_psi1_psi2
      USE mod_Op,    ONLY : param_Op, sub_PsiOpPsi, sub_OpPsi,sub_scaledOpPsi
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)                      :: para_H

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa),    intent(inout)  :: para_propa
      TYPE (param_psi),      intent(in)     :: psi0
      TYPE (param_psi),      intent(inout)  :: psi
      integer,               intent(in)     :: no
      real (kind=Rkind),     intent(in)     :: T

!------ working variables ---------------------------------
      complex (kind=Rkind)      :: cdot
      real (kind=Rkind)         :: E0,microT,microdeltaT,phase,microphase
      complex (kind=Rkind)      :: psi0_psiKrylovSpace(para_propa%para_poly%npoly)

      integer                            :: n
      complex (kind=Rkind)               :: H(para_propa%para_poly%npoly,para_propa%para_poly%npoly)
      complex (kind=Rkind), allocatable  :: UPsiOnKrylov(:)
      complex (kind=Rkind), allocatable  :: Vec(:,:)
      real (kind=Rkind),    allocatable  :: Eig(:)

      TYPE (param_psi),     allocatable  :: tab_KrylovSpace(:)
      complex (kind=Rkind)               :: Overlap


      integer              :: it,i,j

!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='march_SIP'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'deltaT',para_propa%WPdeltaT
        write(out_unitp,*) 'E0,Esc',para_H%E0,para_H%Esc
        write(out_unitp,*) 'order',para_propa%para_poly%npoly
      END IF
!-----------------------------------------------------------

      allocate(tab_KrylovSpace(para_propa%para_poly%npoly+1))

      tab_KrylovSpace(1) = psi


      IF (para_propa%nb_micro > 1) THEN
        psi0_psiKrylovSpace(:) = CZERO
        psi0_psiKrylovSpace(1) = Calc_AutoCorr(psi0,tab_KrylovSpace(1),         &
                                         para_propa,T,Write_AC=.FALSE.)
      END IF

      E0 = para_H%E0
      H(:,:) = CZERO
      ! loop for H|psi>, H^2|psi>, H^3|psi>...
      DO j=2,para_propa%para_poly%npoly+1
        IF (debug) write(out_unitp,*) 'in ',name_sub,' it:',j-1
        CALL sub_OpPsi(Psi  =tab_KrylovSpace(j-1),                              &
                       OpPsi=tab_KrylovSpace(j),para_Op=para_H)

        IF (j == 2) THEN
          ! Energy shift, E0, calculation for the first iteration.
          ! since E0=<psi |H| psi> = <tab_KrylovSpace(1) |H|tab_KrylovSpace(1)> =
          ! ..  <tab_KrylovSpace(1) | tab_KrylovSpace(2)>
          ! This shift is important to improve the stability.
          ! => But the phase need to be taking into account at the end of the iterations.
          CALL Overlap_psi1_psi2(Overlap,tab_KrylovSpace(1),tab_KrylovSpace(2))
          E0 = real(Overlap,kind=Rkind)

        END IF
        CALL sub_scaledOpPsi(Psi  =tab_KrylovSpace(j-1),                &
                             OpPsi=tab_KrylovSpace(j),E0=E0,Esc=ONE)

        !make part of the H matrix (related to the vector j-1)
        ! OpPsi of the vector j-1 is in tab_KrylovSpace(j).
        ! That why we need a loop util npoly+1.
        DO i=1,j-1
          CALL Overlap_psi1_psi2(Overlap,tab_KrylovSpace(i),tab_KrylovSpace(j))
          H(i,j-1) = Overlap
          H(j-1,i) = conjg(Overlap)
        END DO

        ! psi(t+dt)=sum_{i}^{n} <Vec|psi_0>exp(-i*Ei*dt) |Vec>
        ! n=j-1
        CALL UPsi_spec(UPsiOnKrylov,H(1:j-1,1:j-1),Vec,Eig,             &
                              para_propa%WPdeltaT,j-1,With_diago=.TRUE.)
        IF (debug) write(out_unitp,*) j-1,'abs(UPsiOnKrylov(j-1)',abs(UPsiOnKrylov(j-1))
        IF (abs(UPsiOnKrylov(j-1)) < para_propa%para_poly%poly_tol .OR. &
            j == para_propa%para_poly%npoly+1) THEN
          n = j-1
          write(out_unitp,*) n,'abs(UPsiOnKrylov(n)',abs(UPsiOnKrylov(j-1))
          EXIT
        END IF

        ! ortho, the Krylov space: new vector
        !2d version (more stable than the 1st and the 3d ones)
        CALL renorm_psi(tab_KrylovSpace(j))
        DO i=1,j-1
          CALL Overlap_psi1_psi2(Overlap,tab_KrylovSpace(i),tab_KrylovSpace(j))
          IF (abs(Overlap) == ZERO) CYCLE
          tab_KrylovSpace(j) = tab_KrylovSpace(j) - tab_KrylovSpace(i) * Overlap
          CALL renorm_psi(tab_KrylovSpace(j))
        END DO
        DO i=1,j-1
          CALL Overlap_psi1_psi2(Overlap,tab_KrylovSpace(i),tab_KrylovSpace(j))
          IF (abs(Overlap) == ZERO) CYCLE
          tab_KrylovSpace(j) = tab_KrylovSpace(j) - tab_KrylovSpace(i) * Overlap
          CALL renorm_psi(tab_KrylovSpace(j))
        END DO


        IF (para_propa%nb_micro > 1) THEN
          psi0_psiKrylovSpace(j) = Calc_AutoCorr(psi0,tab_KrylovSpace(j),&
                                          para_propa,T,Write_AC=.FALSE.)
        END IF

      END DO


      IF (abs(UPsiOnKrylov(n)) > para_propa%para_poly%poly_tol) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' The last vector, UPsiOnKrylov(n), coefficient is TOO large'
        write(out_unitp,*) '    abs(UPsiOnKrylov(n)',abs(UPsiOnKrylov(n))
        write(out_unitp,*) '    poly_tol: ',para_propa%para_poly%poly_tol
        write(out_unitp,*) ' => npoly is TOO small',para_propa%para_poly%npoly
        write(out_unitp,*) ' or'
        write(out_unitp,*) ' => Reduce the time step !!'
        STOP
      END IF

      IF (debug) THEN
        write(out_unitp,*) 'Eig',Eig(1:n)
        write(out_unitp,*) 'abs(UPsiOnKrylov)',abs(UPsiOnKrylov(1:n))
      END IF

      Psi = ZERO
      DO i=1,n
        Psi = Psi + UPsiOnKrylov(i)*tab_KrylovSpace(i)
      END DO

      !- check norm ------------------
      CALL norm2_psi(psi)
      IF ( psi%norm2 > para_propa%max_norm2) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' STOP propagation: norm > max_norm',psi%norm2
        para_propa%march_error   = .TRUE.
        para_propa%test_max_norm = .TRUE.
        STOP
      END IF

      !- Phase Shift -----------------
      phase = E0*para_propa%WPdeltaT
      psi   = psi * exp(-EYE*phase)



      !- autocorelation -----------------
      IF (para_propa%nb_micro > 1) THEN
        microdeltaT = para_propa%WPdeltaT/real(para_propa%nb_micro,kind=Rkind)
        microphase  = phase/real(para_propa%nb_micro,kind=Rkind)

        phase  = ZERO
        microT = ZERO
        DO it=1,para_propa%nb_micro

          microT = microT + microdeltaT
          phase  = phase   + microphase

          CALL UPsi_spec(UPsiOnKrylov,H,Vec,Eig,microT,n,With_diago=.FALSE.)

          cdot = sum(UPsiOnKrylov(1:n)*psi0_psiKrylovSpace(1:n)) ! we cannot use dot_product
          cdot = cdot * exp(-EYE*phase)
          CALL Write_AutoCorr(no,T+microT,cdot)
        END DO
        flush(no)
      ELSE
        cdot = Calc_AutoCorr(psi0,psi,para_propa,T,Write_AC=.FALSE.)
        CALL Write_AutoCorr(no,T + para_propa%WPdeltaT,cdot)
        flush(no)
      END IF


      ! deallocation
      DO i=1,size(tab_KrylovSpace)
         CALL dealloc_psi(tab_KrylovSpace(i))
      END DO
      deallocate(tab_KrylovSpace)
      IF (allocated(Vec))          CALL dealloc_NParray(Vec,'Vec',name_sub)
      IF (allocated(Eig))          CALL dealloc_NParray(Eig,'Eig',name_sub)
      IF (allocated(UPsiOnKrylov)) CALL dealloc_NParray(UPsiOnKrylov,'UPsiOnKrylov',name_sub)


!-----------------------------------------------------------
      IF (debug) THEN
        CALL norm2_psi(psi)
        write(out_unitp,*) 'norm psi',psi%norm2
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

  END SUBROUTINE march_SIP
  SUBROUTINE march_SIP2(T,no,psi,psi0,para_H,para_propa)
      USE mod_system
      USE mod_psi,   ONLY : param_psi,ecri_psi,norm2_psi,renorm_psi,Overlap_psi1_psi2
      USE mod_Op,    ONLY : param_Op, sub_PsiOpPsi, sub_OpPsi,sub_scaledOpPsi
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)                      :: para_H

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa),    intent(inout)  :: para_propa
      TYPE (param_psi),      intent(in)     :: psi0
      TYPE (param_psi),      intent(inout)  :: psi
      integer,               intent(in)     :: no
      real (kind=Rkind),     intent(in)     :: T

!------ working variables ---------------------------------
      complex (kind=Rkind)      :: cdot
      real (kind=Rkind)         :: E0,microT,microdeltaT,phase,microphase
      complex (kind=Rkind)      :: psi0_psiKrylovSpace(para_propa%para_poly%npoly)

      integer                            :: n
      complex (kind=Rkind)               :: H(para_propa%para_poly%npoly,para_propa%para_poly%npoly)
      complex (kind=Rkind)               :: S(para_propa%para_poly%npoly+1,para_propa%para_poly%npoly+1)
      complex (kind=Rkind), allocatable  :: UPsiOnKrylov(:)
      complex (kind=Rkind), allocatable  :: Vec(:,:)
      real (kind=Rkind),    allocatable  :: Eig(:)

      TYPE (param_psi),     allocatable  :: tab_KrylovSpace(:)
      complex (kind=Rkind)               :: Overlap,UPspectral_n


      integer              :: it,i,j

!----- for debuging --------------------------------------------------
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='march_SIP2'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'deltaT',para_propa%WPdeltaT
        write(out_unitp,*) 'E0,Esc',para_H%E0,para_H%Esc
        write(out_unitp,*) 'order',para_propa%para_poly%npoly
        write(out_unitp,*) 'tol',para_propa%para_poly%poly_tol
      END IF
!-----------------------------------------------------------

      allocate(tab_KrylovSpace(para_propa%para_poly%npoly+1))

      tab_KrylovSpace(1) = psi

      E0 = para_H%E0
      H(:,:) = CZERO
      S(:,:) = CZERO
      S(1,1) = CONE
      ! loop for H|psi>, H^2|psi>, H^3|psi>...
      DO j=2,para_propa%para_poly%npoly+1
        IF (debug) write(out_unitp,*) 'in ',name_sub,' it:',j-1
        CALL sub_OpPsi(Psi  =tab_KrylovSpace(j-1),                      &
                       OpPsi=tab_KrylovSpace(j),para_Op=para_H)
!        IF (j == 2) THEN
!          ! Energy shift, E0, calculation for the first iteration.
!          ! since E0=<psi |H| psi> = <tab_KrylovSpace(1) |H|tab_KrylovSpace(1)> =
!          ! ..  <tab_KrylovSpace(1) | tab_KrylovSpace(2)>
!          ! This shift is important to improve the stability.
!          ! => But the phase need to be taking into account at the end of the iterations.
!          CALL Overlap_psi1_psi2(Overlap,tab_KrylovSpace(1),tab_KrylovSpace(2))
!          E0 = real(Overlap,kind=Rkind)
!        END IF
        CALL sub_scaledOpPsi(Psi  =tab_KrylovSpace(j-1),                &
                             OpPsi=tab_KrylovSpace(j),E0=E0,Esc=ONE)

        !make part of the S matrix (related to the vector j)
        ! and then H (part of S)
        DO i=1,j
          CALL Overlap_psi1_psi2(Overlap,tab_KrylovSpace(i),tab_KrylovSpace(j))
          S(i,j) = Overlap
          S(j,i) = conjg(Overlap)
        END DO
        H(1:j-1,1:j-1) = S(1:j-1,2:j)

        write(out_unitp,*)
        write(out_unitp,*) 'Real(S)'
        CALL Write_Mat(real(S(1:j,1:j),kind=Rkind),out_unitp,6,Rformat='(f18.14)')
        write(out_unitp,*) 'im(S)'
        CALL Write_Mat(aimag(S(1:j,1:j)),out_unitp,6,Rformat='(f18.14)')


        ! psi(t+dt)=sum_{i}^{n} <Vec|psi_0>exp(-i*Ei*dt) |Vec>
        ! n=j-1
        CALL UGPsi_spec(UPsiOnKrylov,UPspectral_n,H(1:j-1,1:j-1),       &
                        S(1:j-1,1:j-1),Vec,Eig,para_propa%WPdeltaT,     &
                        j-1,With_diago=.TRUE.)
        !IF (debug) write(out_unitp,*) j-1,'abs(UPspectral_n)',abs(UPspectral_n)
        IF (abs(UPspectral_n) < para_propa%para_poly%poly_tol .OR. &
            j == para_propa%para_poly%npoly+1) THEN
          n = j-1
          write(out_unitp,*) n,'abs(UPspectral_n)',abs(UPspectral_n)
          EXIT
        END IF

      END DO

      IF (abs(UPspectral_n) > para_propa%para_poly%poly_tol) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' The last vector, UPsiOnKrylov(n), coefficient is TOO large'
        write(out_unitp,*) '    abs(UPspectral_n)',abs(UPspectral_n)
        write(out_unitp,*) '    poly_tol: ',para_propa%para_poly%poly_tol
        write(out_unitp,*) ' => npoly is TOO small',para_propa%para_poly%npoly
        write(out_unitp,*) ' or'
        write(out_unitp,*) ' => Reduce the time step !!'
        STOP
      END IF

      IF (debug) THEN
        write(out_unitp,*) 'Eig',Eig(1:n)
        write(out_unitp,*) 'abs(UPsiOnKrylov)',abs(UPsiOnKrylov(1:n))
      END IF

      Psi = ZERO
      DO i=1,n
        Psi = Psi + UPsiOnKrylov(i)*tab_KrylovSpace(i)
      END DO

      !- check norm ------------------
      CALL norm2_psi(psi)
      IF ( psi%norm2 > para_propa%max_norm2) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' STOP propagation: norm > max_norm',psi%norm2
        para_propa%march_error   = .TRUE.
        para_propa%test_max_norm = .TRUE.
        STOP
      END IF

      !- Phase Shift -----------------
      phase = E0*para_propa%WPdeltaT
      psi   = psi * exp(-EYE*phase)
      cdot = Calc_AutoCorr(psi0,psi,para_propa,T,Write_AC=.FALSE.)
      CALL Write_AutoCorr(no,T + para_propa%WPdeltaT,cdot)
      flush(no)


      ! deallocation
      DO i=1,size(tab_KrylovSpace)
         CALL dealloc_psi(tab_KrylovSpace(i))
      END DO
      deallocate(tab_KrylovSpace)
      IF (allocated(Vec))          CALL dealloc_NParray(Vec,'Vec',name_sub)
      IF (allocated(Eig))          CALL dealloc_NParray(Eig,'Eig',name_sub)
      IF (allocated(UPsiOnKrylov)) CALL dealloc_NParray(UPsiOnKrylov,'UPsiOnKrylov',name_sub)

!-----------------------------------------------------------
      IF (debug) THEN
        CALL norm2_psi(psi)
        write(out_unitp,*) 'norm psi',psi%norm2
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

  END SUBROUTINE march_SIP2
  SUBROUTINE march_SIP2_SRep(T,no,psi,psi0,para_H,para_propa)
      USE mod_system
      USE mod_psi
      USE mod_Op
      USE mod_OpPsi_SG4
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)                      :: para_H

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa),    intent(inout)  :: para_propa
      TYPE (param_psi),      intent(in)     :: psi0
      TYPE (param_psi),      intent(inout)  :: psi
      integer,               intent(in)     :: no
      real (kind=Rkind),     intent(in)     :: T

!------ working variables ---------------------------------
      complex (kind=Rkind)      :: cdot
      real (kind=Rkind)         :: E0,phase

      integer                             :: n
      complex (kind=Rkind)                :: H(para_propa%para_poly%npoly,para_propa%para_poly%npoly)
      complex (kind=Rkind)                :: S(para_propa%para_poly%npoly+1,para_propa%para_poly%npoly+1)
      complex (kind=Rkind), allocatable   :: UPsiOnKrylov(:)
      complex (kind=Rkind)                :: UPspectral_n

      complex (kind=Rkind), allocatable   :: Vec(:,:)
      real (kind=Rkind),    allocatable   :: Eig(:)

      real (kind=Rkind),    allocatable   :: S11,S12,S21,S22
      integer                             :: iG,tab_l(psi%BasisnD%nb_basis)

      TYPE (param_psi)                    :: RCPsi(2),RCpsi0(2)
      TYPE (Type_SmolyakRep)              :: PsiSRep(2),Psi0SRep(2)
      TYPE (Type_SmolyakRep), allocatable :: tab_KrylovSpace(:,:)

      TYPE (Basis),  pointer    :: BasisnD       => null()   ! .TRUE. pointer

      integer              :: it,i,j

!----- for debuging --------------------------------------------------
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='march_SIP2_SRep'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'deltaT',para_propa%WPdeltaT
        write(out_unitp,*) 'E0,Esc',para_H%E0,para_H%Esc
        write(out_unitp,*) 'order',para_propa%para_poly%npoly
      END IF
!-----------------------------------------------------------
     BasisnD => psi%BasisnD


      RCPsi  = psi
      !RCPsi0 = psi0

      CALL tabPackedBasis_TO_SmolyakRepBasis(PsiSRep(1),RCPsi(1)%RVecB, &
                                BasisnD%tab_basisPrimSG,BasisnD%nDindB, &
                                BasisnD%para_SGType2)
      CALL tabPackedBasis_TO_SmolyakRepBasis(PsiSRep(2),RCPsi(2)%RVecB, &
                                BasisnD%tab_basisPrimSG,BasisnD%nDindB, &
                                BasisnD%para_SGType2)

      allocate(tab_KrylovSpace(2,para_propa%para_poly%npoly+1)) ! the 2 for the Real and imaginary part.

      tab_KrylovSpace(1,1) = PsiSRep(1)
      tab_KrylovSpace(2,1) = PsiSRep(2)


      E0 = para_H%E0
      H(:,:) = CZERO
      S(:,:) = CZERO
      S(1,1) = CONE

      ! loop for H|psi>, H^2|psi>, H^3|psi>...
      DO j=2,para_propa%para_poly%npoly+1
        IF (debug) write(out_unitp,*) 'in ',name_sub,' it:',j-1

        ! loop on the Smolyak terms
        tab_l(:) = BasisnD%para_SGType2%nDval_init(:,1)
        DO iG=1,BasisnD%para_SGType2%nb_SG
          CALL ADD_ONE_TO_nDindex(BasisnD%para_SGType2%nDind_SmolyakRep,tab_l,iG=iG)

          ! H psi (we must copy because sub_OpPsi_OF_ONEDP_FOR_SGtype4 works on the one vector: psi=Hpsi)
          tab_KrylovSpace(1,j) = tab_KrylovSpace(1,j-1)
          tab_KrylovSpace(2,j) = tab_KrylovSpace(2,j-1)

          CALL sub_OpPsi_OF_ONEDP_FOR_SGtype4(tab_KrylovSpace(1,j)%SmolyakRep(iG),  &
                                              iG,tab_l,para_H)
          CALL sub_OpPsi_OF_ONEDP_FOR_SGtype4(tab_KrylovSpace(2,j)%SmolyakRep(iG),  &
                                              iG,tab_l,para_H)

          !shifting HPsi (E0 is real)
          tab_KrylovSpace(1,j)%SmolyakRep(iG)%V = tab_KrylovSpace(1,j  )%SmolyakRep(iG)%V - &
                                               E0*tab_KrylovSpace(1,j-1)%SmolyakRep(iG)%V
          tab_KrylovSpace(2,j)%SmolyakRep(iG)%V = tab_KrylovSpace(2,j  )%SmolyakRep(iG)%V - &
                                               E0*tab_KrylovSpace(2,j-1)%SmolyakRep(iG)%V

          ! update the S matrix (related to the vector j)
          DO i=1,j
            !CALL Overlap_psi1_psi2(Overlap,tab_KrylovSpace(i),tab_KrylovSpace(j))
            !Overlap=< V(1,i)-i*V(2,i) | V(1,j)+i*V(2,j)> = (S11+S22)+i*(S12-S21)
            S11 = dot_product(tab_KrylovSpace(1,i)%SmolyakRep(iG)%V,tab_KrylovSpace(1,j)%SmolyakRep(iG)%V)
            S22 = dot_product(tab_KrylovSpace(2,i)%SmolyakRep(iG)%V,tab_KrylovSpace(2,j)%SmolyakRep(iG)%V)
            S12 = dot_product(tab_KrylovSpace(1,i)%SmolyakRep(iG)%V,tab_KrylovSpace(2,j)%SmolyakRep(iG)%V)
            S21 = dot_product(tab_KrylovSpace(2,i)%SmolyakRep(iG)%V,tab_KrylovSpace(1,j)%SmolyakRep(iG)%V)

            S(i,j) = S(i,j) + BasisnD%WeightSG(iG)*cmplx(S11+S22, S12-S21,kind=Rkind)

            S(j,i) = conjg(S(i,j))
            !S11 = dot_product(tab_KrylovSpace(1,j)%SmolyakRep(iG)%V,tab_KrylovSpace(1,i)%SmolyakRep(iG)%V)
            !S22 = dot_product(tab_KrylovSpace(2,j)%SmolyakRep(iG)%V,tab_KrylovSpace(2,i)%SmolyakRep(iG)%V)
            !S12 = dot_product(tab_KrylovSpace(1,j)%SmolyakRep(iG)%V,tab_KrylovSpace(2,i)%SmolyakRep(iG)%V)
            !S21 = dot_product(tab_KrylovSpace(2,j)%SmolyakRep(iG)%V,tab_KrylovSpace(1,i)%SmolyakRep(iG)%V)
            !S(j,i) = S(j,i) + BasisnD%WeightSG(iG)*cmplx(S11+S22, S12-S21,kind=Rkind)

          END DO

        END DO
        ! End of the Smolyak term loop

        write(out_unitp,*)
        write(out_unitp,*) 'Real(S)'
        CALL Write_Mat(real(S(1:j,1:j),kind=Rkind),out_unitp,6,Rformat='(f18.14)')
        write(out_unitp,*) 'im(S)'
        CALL Write_Mat(aimag(S(1:j,1:j)),out_unitp,6,Rformat='(f18.14)')

        ! The matrix H is a part of S
        H(1:j-1,1:j-1) = S(1:j-1,2:j)
        !H(1:j-1,1:j-1) = HALF*(H(1:j-1,1:j-1) + conjg(transpose(H(1:j-1,1:j-1))))

        ! psi(t+dt)=sum_{i}^{n} <Vec|psi_0>exp(-i*Ei*dt) |Vec>
        ! n=j-1
        CALL UGPsi_spec(UPsiOnKrylov,UPspectral_n,H(1:j-1,1:j-1),       &
                        S(1:j-1,1:j-1),Vec,Eig,para_propa%WPdeltaT,     &
                        j-1,With_diago=.TRUE.)
        !IF (debug) write(out_unitp,*) j-1,'abs(UPspectral_n)',abs(UPspectral_n)
        IF (abs(UPspectral_n) < para_propa%para_poly%poly_tol .OR. &
            j == para_propa%para_poly%npoly+1) THEN
          n = j-1
          write(out_unitp,*) n,'abs(UPspectral_n)',abs(UPspectral_n)
          EXIT
        END IF

      END DO

      IF (abs(UPspectral_n) > para_propa%para_poly%poly_tol) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' The last vector, UPsiOnKrylov(n), coefficient is TOO large'
        write(out_unitp,*) '    abs(UPspectral_n)',abs(UPspectral_n)
        write(out_unitp,*) '    poly_tol: ',para_propa%para_poly%poly_tol
        write(out_unitp,*) ' => npoly is TOO small',para_propa%para_poly%npoly
        write(out_unitp,*) ' or'
        write(out_unitp,*) ' => Reduce the time step !!'
        STOP
      END IF

      IF (debug) THEN
        write(out_unitp,*) 'Eig',Eig(1:n)
        write(out_unitp,*) 'abs(UPsiOnKrylov)',abs(UPsiOnKrylov(1:n))
      END IF

      !Psi = ZERO
      PsiSRep(1) = ZERO
      PsiSRep(2) = ZERO
      DO i=1,n
        !Psi = Psi + UPsiOnKrylov(i)*tab_KrylovSpace(i)
        PsiSRep(1) = PsiSRep(1) +                                       &
                real(UPsiOnKrylov(i),kind=Rkind)*tab_KrylovSpace(1,i) - &
                aimag(UPsiOnKrylov(i))*tab_KrylovSpace(2,i)

        PsiSRep(2) = PsiSRep(2) +                                       &
                real(UPsiOnKrylov(i),kind=Rkind)*tab_KrylovSpace(2,i) + &
                aimag(UPsiOnKrylov(i))*tab_KrylovSpace(1,i)

      END DO

      CALL SmolyakRepBasis_TO_tabPackedBasis(PsiSRep(1),RCPsi(1)%RVecB, &
                                BasisnD%tab_basisPrimSG,BasisnD%nDindB, &
                                BasisnD%para_SGType2,BasisnD%WeightSG)
      CALL SmolyakRepBasis_TO_tabPackedBasis(PsiSRep(2),RCPsi(2)%RVecB, &
                                BasisnD%tab_basisPrimSG,BasisnD%nDindB, &
                                BasisnD%para_SGType2,BasisnD%WeightSG)

      psi = RCPsi


      !- check norm ------------------
      CALL norm2_psi(psi)
      IF ( psi%norm2 > para_propa%max_norm2) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' STOP propagation: norm > max_norm',psi%norm2
        para_propa%march_error   = .TRUE.
        para_propa%test_max_norm = .TRUE.
        STOP
      END IF

      !- Phase Shift -----------------
      phase = E0*para_propa%WPdeltaT
      psi   = psi * exp(-EYE*phase)
      cdot = Calc_AutoCorr(psi0,psi,para_propa,T,Write_AC=.FALSE.)
      CALL Write_AutoCorr(no,T + para_propa%WPdeltaT,cdot)
      flush(no)


      ! deallocation
      DO i=1,size(tab_KrylovSpace(1,:))
         CALL dealloc_SmolyakRep(tab_KrylovSpace(1,i))
         CALL dealloc_SmolyakRep(tab_KrylovSpace(2,i))
      END DO
      deallocate(tab_KrylovSpace)
       CALL dealloc_psi(RCPsi(1))
       CALL dealloc_psi(RCPsi(2))

      IF (allocated(Vec))          CALL dealloc_NParray(Vec,'Vec',name_sub)
      IF (allocated(Eig))          CALL dealloc_NParray(Eig,'Eig',name_sub)
      IF (allocated(UPsiOnKrylov)) CALL dealloc_NParray(UPsiOnKrylov,'UPsiOnKrylov',name_sub)

!-----------------------------------------------------------
      IF (debug) THEN
        CALL norm2_psi(psi)
        write(out_unitp,*) 'norm psi',psi%norm2
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

  END SUBROUTINE march_SIP2_SRep
  SUBROUTINE march_SIP2_SRep2(T,no,psi,psi0,para_H,para_propa)
      USE mod_system
      USE mod_psi
      USE mod_Op
      USE mod_OpPsi_SG4
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)                      :: para_H

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa),    intent(inout)  :: para_propa
      TYPE (param_psi),      intent(in)     :: psi0
      TYPE (param_psi),      intent(inout)  :: psi
      integer,               intent(in)     :: no
      real (kind=Rkind),     intent(in)     :: T

!------ working variables ---------------------------------
      complex (kind=Rkind)      :: cdot
      real (kind=Rkind)         :: E0,phase

      integer                             :: n
      complex (kind=Rkind)                :: H(para_propa%para_poly%npoly,para_propa%para_poly%npoly)
      complex (kind=Rkind)                :: S(para_propa%para_poly%npoly+1,para_propa%para_poly%npoly+1)
      complex (kind=Rkind), allocatable   :: UPsiOnKrylov(:)
      complex (kind=Rkind)                :: UPspectral_n

      complex (kind=Rkind), allocatable   :: Vec(:,:)
      real (kind=Rkind),    allocatable   :: Eig(:)

      real (kind=Rkind),    allocatable   :: S11,S12,S21,S22
      integer                             :: iG,tab_l(psi%BasisnD%nb_basis)

      TYPE (param_psi)                    :: RCPsi(2),RCpsi0(2)
      TYPE (Type_SmolyakRep)              :: PsiSRep(2),Psi0SRep(2)
      TYPE (Type_SmolyakRep), allocatable :: tab_KrylovSpace(:,:)
      TYPE (Type_SmolyakRep)              :: HPsiSRep(2)

      TYPE (Basis),  pointer    :: BasisnD       => null()   ! .TRUE. pointer

      integer              :: it,i,j

!----- for debuging --------------------------------------------------
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='march_SIP2_SRep2'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'deltaT',para_propa%WPdeltaT
        write(out_unitp,*) 'E0,Esc',para_H%E0,para_H%Esc
        write(out_unitp,*) 'order',para_propa%para_poly%npoly
      END IF
!-----------------------------------------------------------
     BasisnD => psi%BasisnD

      RCPsi  = psi
      !RCPsi0 = psi0

      CALL tabPackedBasis_TO_SmolyakRepBasis(PsiSRep(1),RCPsi(1)%RVecB, &
                                BasisnD%tab_basisPrimSG,BasisnD%nDindB, &
                                BasisnD%para_SGType2)
      CALL tabPackedBasis_TO_SmolyakRepBasis(PsiSRep(2),RCPsi(2)%RVecB, &
                                BasisnD%tab_basisPrimSG,BasisnD%nDindB, &
                                BasisnD%para_SGType2)

      allocate(tab_KrylovSpace(2,para_propa%para_poly%npoly+1)) ! the 2 for the Real and imaginary part.

      tab_KrylovSpace(1,1) = PsiSRep(1)
      tab_KrylovSpace(2,1) = PsiSRep(2)


      E0 = para_H%E0
      H(:,:) = CZERO
      S(:,:) = CZERO
      S(1,1) = CONE

      ! loop for H|psi>, H^2|psi>, H^3|psi>...
      DO j=2,para_propa%para_poly%npoly+1
        IF (debug) write(out_unitp,*) 'in ',name_sub,' it:',j-1

        ! loop on the Smolyak terms
        tab_l(:) = BasisnD%para_SGType2%nDval_init(:,1)
        DO iG=1,BasisnD%para_SGType2%nb_SG
          CALL ADD_ONE_TO_nDindex(BasisnD%para_SGType2%nDind_SmolyakRep,tab_l,iG=iG)

          ! H psi (we must copy because sub_OpPsi_OF_ONEDP_FOR_SGtype4 works on the one vector: psi=Hpsi)
          tab_KrylovSpace(1,j) = tab_KrylovSpace(1,j-1)
          tab_KrylovSpace(2,j) = tab_KrylovSpace(2,j-1)

          CALL sub_OpPsi_OF_ONEDP_FOR_SGtype4(tab_KrylovSpace(1,j)%SmolyakRep(iG),  &
                                              iG,tab_l,para_H)
          CALL sub_OpPsi_OF_ONEDP_FOR_SGtype4(tab_KrylovSpace(2,j)%SmolyakRep(iG),  &
                                              iG,tab_l,para_H)

          !shifting HPsi (E0 is real)
          tab_KrylovSpace(1,j)%SmolyakRep(iG)%V = tab_KrylovSpace(1,j  )%SmolyakRep(iG)%V - &
                                               E0*tab_KrylovSpace(1,j-1)%SmolyakRep(iG)%V
          tab_KrylovSpace(2,j)%SmolyakRep(iG)%V = tab_KrylovSpace(2,j  )%SmolyakRep(iG)%V - &
                                               E0*tab_KrylovSpace(2,j-1)%SmolyakRep(iG)%V

          ! update the S matrix (related to the vector j)
          DO i=1,j
            !CALL Overlap_psi1_psi2(Overlap,tab_KrylovSpace(i),tab_KrylovSpace(j))
            !Overlap=< V(1,i)-i*V(2,i) | V(1,j)+i*V(2,j)> = (S11+S22)+i*(S12-S21)
            S11 = dot_product(tab_KrylovSpace(1,i)%SmolyakRep(iG)%V,tab_KrylovSpace(1,j)%SmolyakRep(iG)%V)
            S22 = dot_product(tab_KrylovSpace(2,i)%SmolyakRep(iG)%V,tab_KrylovSpace(2,j)%SmolyakRep(iG)%V)
            S12 = dot_product(tab_KrylovSpace(1,i)%SmolyakRep(iG)%V,tab_KrylovSpace(2,j)%SmolyakRep(iG)%V)
            S21 = dot_product(tab_KrylovSpace(2,i)%SmolyakRep(iG)%V,tab_KrylovSpace(1,j)%SmolyakRep(iG)%V)

            S(i,j) = S(i,j) + BasisnD%WeightSG(iG)*cmplx(S11+S22, S12-S21,kind=Rkind)

            IF (j /= i) S(j,i) = conjg(S(i,j))
            !S11 = dot_product(tab_KrylovSpace(1,j)%SmolyakRep(iG)%V,tab_KrylovSpace(1,i)%SmolyakRep(iG)%V)
            !S22 = dot_product(tab_KrylovSpace(2,j)%SmolyakRep(iG)%V,tab_KrylovSpace(2,i)%SmolyakRep(iG)%V)
            !S12 = dot_product(tab_KrylovSpace(1,j)%SmolyakRep(iG)%V,tab_KrylovSpace(2,i)%SmolyakRep(iG)%V)
            !S21 = dot_product(tab_KrylovSpace(2,j)%SmolyakRep(iG)%V,tab_KrylovSpace(1,i)%SmolyakRep(iG)%V)
            !S(j,i) = S(j,i) + BasisnD%WeightSG(iG)*cmplx(S11+S22, S12-S21,kind=Rkind)

          END DO


          ! H psi (we must copy because sub_OpPsi_OF_ONEDP_FOR_SGtype4 works on the one vector: psi=Hpsi)
          ! => for j-1
          HPsiSRep(1) = tab_KrylovSpace(1,j-1)
          HPsiSRep(2) = tab_KrylovSpace(2,j-1)

          CALL sub_OpPsi_OF_ONEDP_FOR_SGtype4(HPsiSRep(1)%SmolyakRep(iG),iG,tab_l,para_H)
          CALL sub_OpPsi_OF_ONEDP_FOR_SGtype4( HPsiSRep(2)%SmolyakRep(iG),iG,tab_l,para_H)

          !shifting HPsi (E0 is real)
          HPsiSRep(1)%SmolyakRep(iG)%V = HPsiSRep(1)%SmolyakRep(iG)%V - &
                           E0*tab_KrylovSpace(1,j-1)%SmolyakRep(iG)%V
          HPsiSRep(2)%SmolyakRep(iG)%V = HPsiSRep(2)%SmolyakRep(iG)%V - &
                           E0*tab_KrylovSpace(2,j-1)%SmolyakRep(iG)%V

          DO i=1,j-1
            !CALL Overlap_psi1_psi2(Overlap,tab_KrylovSpace(i),tab_KrylovSpace(j))
            !Overlap=< V(1,i)-i*V(2,i) | V(1,j)+i*V(2,j)> = (S11+S22)+i*(S12-S21)
            S11 = dot_product(tab_KrylovSpace(1,i)%SmolyakRep(iG)%V,HPsiSRep(1)%SmolyakRep(iG)%V)
            S22 = dot_product(tab_KrylovSpace(2,i)%SmolyakRep(iG)%V,HPsiSRep(2)%SmolyakRep(iG)%V)
            S12 = dot_product(tab_KrylovSpace(1,i)%SmolyakRep(iG)%V,HPsiSRep(2)%SmolyakRep(iG)%V)
            S21 = dot_product(tab_KrylovSpace(2,i)%SmolyakRep(iG)%V,HPsiSRep(1)%SmolyakRep(iG)%V)

            H(i,j-1) = H(i,j-1) + BasisnD%WeightSG(iG)*cmplx(S11+S22, S12-S21,kind=Rkind)

            IF (i /= (j-1) ) H(j-1,i) = conjg(H(i,j-1)) !!????

          END DO

        END DO
        ! End of the Smolyak term loop

        write(out_unitp,*)
        write(out_unitp,*) 'Real(S)'
        CALL Write_Mat(real(S(1:j,1:j),kind=Rkind),out_unitp,6,Rformat='(f18.14)')
        write(out_unitp,*) 'im(S)'
        CALL Write_Mat(aimag(S(1:j,1:j)),out_unitp,6,Rformat='(f18.14)')

        ! The matrix H is a part of S
        !H(1:j-1,1:j-1) = S(1:j-1,2:j)
        !H(1:j-1,1:j-1) = HALF*(H(1:j-1,1:j-1) + conjg(transpose(H(1:j-1,1:j-1))))

         !DO i=1,j-1
         !  H(i,i) = real(H(i,i),kind=Rkind)
         !END DO

        ! psi(t+dt)=sum_{i}^{n} <Vec|psi_0>exp(-i*Ei*dt) |Vec>
        ! n=j-1
        CALL UGPsi_spec(UPsiOnKrylov,UPspectral_n,H(1:j-1,1:j-1),     &
                        S(1:j-1,1:j-1),Vec,Eig,para_propa%WPdeltaT,     &
                        j-1,With_diago=.TRUE.)
        !IF (debug) write(out_unitp,*) j-1,'abs(UPspectral_n)',abs(UPspectral_n)
        IF (abs(UPspectral_n) < para_propa%para_poly%poly_tol .OR. &
            j == para_propa%para_poly%npoly+1) THEN
          n = j-1
          write(out_unitp,*) n,'abs(UPspectral_n)',abs(UPspectral_n)
          EXIT
        END IF

      END DO

      IF (abs(UPspectral_n) > para_propa%para_poly%poly_tol) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' The last vector, UPsiOnKrylov(n), coefficient is TOO large'
        write(out_unitp,*) '    abs(UPspectral_n)',abs(UPspectral_n)
        write(out_unitp,*) '    poly_tol: ',para_propa%para_poly%poly_tol
        write(out_unitp,*) ' => npoly is TOO small',para_propa%para_poly%npoly
        write(out_unitp,*) ' or'
        write(out_unitp,*) ' => Reduce the time step !!'
        !STOP
      END IF

      IF (debug) THEN
        write(out_unitp,*) 'Eig',Eig(1:n)
        write(out_unitp,*) 'abs(UPsiOnKrylov)',abs(UPsiOnKrylov(1:n))
      END IF

      !Psi = ZERO
      PsiSRep(1) = ZERO
      PsiSRep(2) = ZERO
      DO i=1,n
        !Psi = Psi + UPsiOnKrylov(i)*tab_KrylovSpace(i)
        PsiSRep(1) = PsiSRep(1) +                                       &
                real(UPsiOnKrylov(i),kind=Rkind)*tab_KrylovSpace(1,i) - &
                aimag(UPsiOnKrylov(i))*tab_KrylovSpace(2,i)

        PsiSRep(2) = PsiSRep(2) +                                       &
                real(UPsiOnKrylov(i),kind=Rkind)*tab_KrylovSpace(2,i) + &
                aimag(UPsiOnKrylov(i))*tab_KrylovSpace(1,i)

      END DO

      CALL SmolyakRepBasis_TO_tabPackedBasis(PsiSRep(1),RCPsi(1)%RVecB, &
                                BasisnD%tab_basisPrimSG,BasisnD%nDindB, &
                                BasisnD%para_SGType2,BasisnD%WeightSG)
      CALL SmolyakRepBasis_TO_tabPackedBasis(PsiSRep(2),RCPsi(2)%RVecB, &
                                BasisnD%tab_basisPrimSG,BasisnD%nDindB, &
                                BasisnD%para_SGType2,BasisnD%WeightSG)

      psi = RCPsi


      !- check norm ------------------
      CALL norm2_psi(psi)
      IF ( psi%norm2 > para_propa%max_norm2) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' STOP propagation: norm > max_norm',psi%norm2
        para_propa%march_error   = .TRUE.
        para_propa%test_max_norm = .TRUE.
        STOP
      END IF

      !- Phase Shift -----------------
      phase = E0*para_propa%WPdeltaT
      psi   = psi * exp(-EYE*phase)
      cdot = Calc_AutoCorr(psi0,psi,para_propa,T,Write_AC=.FALSE.)
      CALL Write_AutoCorr(no,T + para_propa%WPdeltaT,cdot)
      flush(no)


      ! deallocation
      DO i=1,size(tab_KrylovSpace(1,:))
         CALL dealloc_SmolyakRep(tab_KrylovSpace(1,i))
         CALL dealloc_SmolyakRep(tab_KrylovSpace(2,i))
      END DO
      deallocate(tab_KrylovSpace)

      CALL dealloc_SmolyakRep(HPsiSRep(1))
      CALL dealloc_SmolyakRep(HPsiSRep(2))

       CALL dealloc_psi(RCPsi(1))
       CALL dealloc_psi(RCPsi(2))

      IF (allocated(Vec))          CALL dealloc_NParray(Vec,'Vec',name_sub)
      IF (allocated(Eig))          CALL dealloc_NParray(Eig,'Eig',name_sub)
      IF (allocated(UPsiOnKrylov)) CALL dealloc_NParray(UPsiOnKrylov,'UPsiOnKrylov',name_sub)

!-----------------------------------------------------------
      IF (debug) THEN
        CALL norm2_psi(psi)
        write(out_unitp,*) 'norm psi',psi%norm2
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

  END SUBROUTINE march_SIP2_SRep2
  SUBROUTINE march_SIL(T,no,psi,psi0,para_H,para_propa)
      USE mod_system
      USE mod_Op,    ONLY : param_Op, sub_PsiOpPsi, sub_OpPsi,sub_scaledOpPsi
      USE mod_psi,   ONLY : param_psi,ecri_psi,norm2_psi,renorm_psi,    &
                            renorm_psi_With_norm2,Overlap_psi1_psi2
      USE mod_MPI_aux
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)     :: para_H

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa),    intent(inout)  :: para_propa
      TYPE (param_psi),      intent(in)     :: psi0
      TYPE (param_psi),      intent(inout)  :: psi
      integer,               intent(in)     :: no
      real (kind=Rkind),     intent(in)     :: T

!------ working variables ---------------------------------
      complex (kind=Rkind)      :: cdot
      real (kind=Rkind)         :: E0,microT,microdeltaT,phase,microphase
      complex (kind=Rkind)      :: psi0_psiKrylovSpace(para_propa%para_poly%npoly)

      integer                            :: n
      complex (kind=Rkind)               :: H(para_propa%para_poly%npoly,para_propa%para_poly%npoly)
      complex (kind=Rkind), allocatable  :: UPsiOnKrylov(:)
      complex (kind=Rkind), allocatable  :: Vec(:,:)
      real (kind=Rkind),    allocatable  :: Eig(:)

      TYPE (param_psi),     allocatable  :: tab_KrylovSpace(:)
      complex (kind=Rkind)               :: Overlap
      TYPE (param_psi)                   :: w1


      integer              :: k,it
      Logical              :: MPI_exit=.FALSE.

!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='march_SIL'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'deltaT',para_propa%WPdeltaT
        write(out_unitp,*) 'E0,Esc',para_H%E0,para_H%Esc
        write(out_unitp,*) 'order',para_propa%para_poly%npoly
      END IF
!-----------------------------------------------------------

      MPI_exit=.FALSE.

      allocate(tab_KrylovSpace(para_propa%para_poly%npoly+1))

      IF(keep_MPI) tab_KrylovSpace(1) = psi

      IF (para_propa%nb_micro > 1) THEN
        IF(keep_MPI) psi0_psiKrylovSpace(:) = CZERO
        IF(keep_MPI) psi0_psiKrylovSpace(1) = Calc_AutoCorr(psi0,tab_KrylovSpace(1), &
                                         para_propa,T,Write_AC=.FALSE.)
      END IF

      E0 = para_H%E0
      IF(keep_MPI) H(:,:) = CZERO
      DO k=1,para_propa%para_poly%npoly
        IF (debug) write(out_unitp,*) 'in ',name_sub,' it:',k

        IF (para_propa%nb_micro > 1) THEN
          IF(keep_MPI) psi0_psiKrylovSpace(k) = Calc_AutoCorr(psi0,tab_KrylovSpace(k),&
                                          para_propa,T,Write_AC=.FALSE.)
        END IF

        !alpha_k calculation. Rq: |w1> = H| tab_KrylovSpace(k)>
        CALL sub_OpPsi(Psi=tab_KrylovSpace(k),OpPsi=w1,para_Op=para_H)
        IF (k == 1) THEN
          ! Energy shift, E0, calculation for the first iteration.
          ! since E0=<psi |H| psi> = <tab_KrylovSpace(1) |H|tab_KrylovSpace(1)> =
          ! ..  <tab_KrylovSpace(1) | w1>
          ! This shift is important to improve the stapility.
          ! => But the phase need to be taking into account at the end of the iterations.
          IF(keep_MPI) CALL Overlap_psi1_psi2(Overlap,tab_KrylovSpace(1),w1)
          CALL MPI_Bcast_(Overlap,size1_MPI,root_MPI)
          E0 = real(Overlap,kind=Rkind)
        END IF
        CALL sub_scaledOpPsi(Psi=tab_KrylovSpace(k),OpPsi=w1,E0=E0,Esc=ONE)

        ! |w_k> = H| tab_KrylovSpace(k)> - beta_k-1 * |tab_KrylovSpace(k-1)>
        ! |w_k> = |w1>                   -MatH(k,k-1)*|tab_KrylovSpace(k-1)>
        IF (k == 1) THEN
          ! w1 = w1
        ELSE
          IF(keep_MPI) w1 = w1 - H(k,k-1) * tab_KrylovSpace(k-1)
        END IF

        ! alpha_k = <w_k | tab_KrylovSpace(k)> = <w1 | tab_KrylovSpace(k)>
        IF(keep_MPI) THEN
          CALL Overlap_psi1_psi2(Overlap,w1,tab_KrylovSpace(k))
          H(k,k) = Overlap

          !diagonalization then exit when:
          !  (i) k reaches npoly (ii) the scheme converges
          CALL UPsi_spec(UPsiOnKrylov,H(1:k,1:k),Vec,Eig,                 &
                                  para_propa%WPdeltaT,k,With_diago=.TRUE.)
          !write(out_unitp,*) k,'abs(UPsiOnKrylov(k)',abs(UPsiOnKrylov(k))
          IF (abs(UPsiOnKrylov(k)) < para_propa%para_poly%poly_tol .OR. &
            k == para_propa%para_poly%npoly) THEN
            n = k
            write(out_unitp,*) n,'abs(UPsiOnKrylov(n)',abs(UPsiOnKrylov(n))
            !EXIT
            MPI_exit=.TRUE.
          END IF
        ENDIF

        IF(openmpi) CALL MPI_Bcast_(MPI_exit,size1_MPI,root_MPI)
        IF(MPI_exit) EXIT

        ! orthogonalisation: |v_k> = |w_k> - alpha_k * |tab_KrylovSpace(k)>
        !                     |w1> = |w1>  - Overlap * |tab_KrylovSpace(k)>
        IF(keep_MPI) w1 = w1 - Overlap*tab_KrylovSpace(k)

        !beta_k = <v_k | v_k> = <w1 | w1>
        !       => beta_k = sqrt(norm2)
        IF(keep_MPI) THEN
          CALL norm2_psi(w1)
          H(k+1,k) = sqrt(w1%norm2)
          H(k,k+1) = H(k+1,k)
          CALL renorm_psi_With_norm2(w1)
          tab_KrylovSpace(k+1) = w1
        ENDIF

      END DO

      IF (debug) write(out_unitp,*) 'abs(UPsiOnKrylov)',abs(UPsiOnKrylov(1:n))

      IF(keep_MPI) THEN
        Psi = ZERO
        DO k=1,n
          Psi = Psi + UPsiOnKrylov(k)*tab_KrylovSpace(k)
        END DO
      ENDIF

      !- check norm ------------------
      IF(keep_MPI) CALL norm2_psi(psi)
      IF(openmpi)  CALL MPI_Bcast_(psi%norm2,size1_MPI,root_MPI)
      IF ( psi%norm2 > para_propa%max_norm2) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' STOP propagation: norm > max_norm',                &
                         psi%norm2
        para_propa%march_error   = .TRUE.
        para_propa%test_max_norm = .TRUE.
        STOP
      END IF

      !- Phase Shift -----------------
      phase = E0*para_propa%WPdeltaT
      IF(keep_MPI) psi   = psi * exp(-EYE*phase)

      !- autocorelation -----------------
      IF(MPI_id==0) THEN
      IF (para_propa%nb_micro > 1) THEN
        microdeltaT = para_propa%WPdeltaT/real(para_propa%nb_micro,kind=Rkind)
        microphase  = phase/real(para_propa%nb_micro,kind=Rkind)

        phase  = ZERO
        microT = ZERO
        DO it=1,para_propa%nb_micro

          microT = microT + microdeltaT
          phase  = phase   + microphase

          CALL UPsi_spec(UPsiOnKrylov,H,Vec,Eig,microT,n,With_diago=.FALSE.)

          cdot = sum(UPsiOnKrylov(1:n)*psi0_psiKrylovSpace(1:n)) ! we cannot use dot_product
          cdot = cdot * exp(-EYE*phase)
          CALL Write_AutoCorr(no,T+microT,cdot)
        END DO
        flush(no)
      ELSE
        cdot = Calc_AutoCorr(psi0,psi,para_propa,T,Write_AC=.FALSE.)
        CALL Write_AutoCorr(no,T + para_propa%WPdeltaT,cdot)
        flush(no)
      END IF
      ENDIF

      ! deallocation
      DO k=1,size(tab_KrylovSpace)
        IF(keep_MPI) CALL dealloc_psi(tab_KrylovSpace(k))
      END DO
      deallocate(tab_KrylovSpace)
      IF (allocated(Vec))          CALL dealloc_NParray(Vec,'Vec',name_sub)
      IF (allocated(Eig))          CALL dealloc_NParray(Eig,'Eig',name_sub)
      IF (allocated(UPsiOnKrylov)) CALL dealloc_NParray(UPsiOnKrylov,'UPsiOnKrylov',name_sub)
      IF(keep_MPI) CALL dealloc_psi(w1)
!-----------------------------------------------------------
      IF (debug) THEN
        CALL norm2_psi(psi)
        write(out_unitp,*) 'norm psi',psi%norm2
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

  END SUBROUTINE march_SIL
  SUBROUTINE UPsi_spec(UPsiOnKrylov,H,Vec,Eig,deltaT,n,With_diago)
      USE mod_system
      IMPLICIT NONE



!----- variables for the WP propagation ----------------------------
      integer,                           intent(in)     :: n
      complex (kind=Rkind),              intent(in)     :: H(n,n)
      complex (kind=Rkind), allocatable, intent(inout)  :: Vec(:,:)
      real (kind=Rkind),    allocatable, intent(inout)  :: Eig(:)
      complex (kind=Rkind), allocatable, intent(inout)  :: UPsiOnKrylov(:)
      real (kind=Rkind),                 intent(in)     :: deltaT
      logical,                           intent(in)     :: With_diago

!------ working variables ---------------------------------
      complex (kind=Rkind) :: coef_i
      integer              :: i

!----- for debuging --------------------------------------------------
      integer, parameter :: nmax = 12
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='UPsi_spec'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'deltaT',deltaT
        write(out_unitp,*) 'n',n
        write(out_unitp,*) 'With_diago',With_diago
        IF (With_diago .AND. n <= nmax) THEN
          write(out_unitp,*) 'H'
          CALL Write_Mat(H,out_unitp,6)
        END IF
      END IF
!-----------------------------------------------------------

      IF (With_diago) THEN
        IF (allocated(Vec)) CALL dealloc_NParray(Vec,'Vec',name_sub)
        CALL alloc_NParray(Vec,[n,n],'Vec',name_sub)

        IF (allocated(Eig)) CALL dealloc_NParray(Eig,'Eig',name_sub)
        CALL alloc_NParray(Eig,[n],  'Eig',name_sub)

        IF (allocated(UPsiOnKrylov))                                    &
              CALL dealloc_NParray(UPsiOnKrylov,'UPsiOnKrylov',name_sub)
        CALL alloc_NParray(UPsiOnKrylov,[n],'UPsiOnKrylov',name_sub)

        CALL diagonalization(H,Eig,Vec,3,1,.TRUE.)
      END IF

      ! loop on the eigenvectors
      UPsiOnKrylov = CZERO
      DO i=1,n
        coef_i = Vec(1,i) ! just (1,i) because psi on the Krylov space is [1,0,0,...0]
        coef_i = coef_i * exp(-EYE*Eig(i)*deltaT) ! spectral propa

        UPsiOnKrylov(:) = UPsiOnKrylov(:) + conjg(Vec(:,i))*coef_i ! update U.psi on the Krylov space

      END DO

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'abs(UPsiOnKrylov(:))',abs(UPsiOnKrylov)
        write(out_unitp,*) 'minval(abs(UPsiOnKrylov(:)))',minval(abs(UPsiOnKrylov))
        write(out_unitp,*) 'norm2 UPsiOnKrylov',dot_product(UPsiOnKrylov,UPsiOnKrylov)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------


  END SUBROUTINE UPsi_spec
  SUBROUTINE UGPsi_spec(UPsiOnKrylov,UPspectral_n,H,S,Vec,Eig,deltaT,n,With_diago)
      USE mod_system
      IMPLICIT NONE



!----- variables for the WP propagation ----------------------------
      integer,                           intent(in)     :: n
      complex (kind=Rkind),              intent(in)     :: H(:,:)
      complex (kind=Rkind),              intent(in)     :: S(:,:)
      complex (kind=Rkind), allocatable, intent(inout)  :: Vec(:,:)
      real (kind=Rkind),    allocatable, intent(inout)  :: Eig(:)
      complex (kind=Rkind), allocatable, intent(inout)  :: UPsiOnKrylov(:)
      real (kind=Rkind),                 intent(in)     :: deltaT
      logical,                           intent(in)     :: With_diago
      complex (kind=Rkind),              intent(inout)  :: UPspectral_n
!------ working variables ---------------------------------
      complex (kind=Rkind) :: coef_i
      integer              :: i,j

      complex (kind=Rkind)     :: Ho(n,n)
      complex (kind=Rkind)     :: So(n,n),SoDiag(n)
      complex (kind=Rkind)     :: C1(n,n)
      complex (kind=Rkind)     :: Sii,Sij
      real (kind=Rkind)        :: max_diff

      real (kind=Rkind)        :: tol_ortho = ONETENTH**10

!----- for debuging --------------------------------------------------
      integer, parameter :: nmax = 12
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='UGPsi_spec'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'deltaT',deltaT
        write(out_unitp,*) 'n',n
        write(out_unitp,*) 'With_diago',With_diago
        IF (With_diago .AND. n <= nmax) THEN
          write(out_unitp,*)
          write(out_unitp,*) 'H'
          write(out_unitp,*) 'Real(H)'
          CALL Write_Mat(real(H,kind=Rkind),out_unitp,6)
          write(out_unitp,*) 'im(H)'
          CALL Write_Mat(aimag(H),out_unitp,6)
          write(out_unitp,*)
          write(out_unitp,*) 'Real(S)'
          CALL Write_Mat(real(S,kind=Rkind),out_unitp,6)
          write(out_unitp,*) 'im(S)'
          CALL Write_Mat(aimag(S),out_unitp,6)
        END IF
      END IF
!-----------------------------------------------------------

      IF (With_diago) THEN
        IF (allocated(Vec)) CALL dealloc_NParray(Vec,'Vec',name_sub)
        CALL alloc_NParray(Vec,[n,n],'Vec',name_sub)

        IF (allocated(Eig)) CALL dealloc_NParray(Eig,'Eig',name_sub)
        CALL alloc_NParray(Eig,[n],  'Eig',name_sub)

        CALL diagonalization(S,Eig,Vec,3,1,.TRUE.)
        IF (debug) write(out_unitp,*) 'Eig of S',Eig
        IF (debug) write(out_unitp,*) 'min(abs(Eig)) of S',minval(abs(Eig))

        !first transfo S -> I (ortho Smidt)
        So = S
        C1 = CZERO
        DO i=1,n
          Sii = S(i,i)
          IF (abs(Sii) <tol_ortho) THEN
            write(out_unitp,*) 'ERROR in ',name_sub
            write(out_unitp,*) 'S(i,i) is too small',abs(Sii)
            STOP 'ERROR in UGPsi_spec: S(i,i) is too small'
          END IF
          C1(i,i) = ONE / sqrt(abs(Sii))
        END DO

        DO i=1,n
          Sii = dot_product(C1(:,i),matmul(S,C1(:,i)))
          IF (abs(Sii) < tol_ortho) THEN
            write(out_unitp,*) 'ERROR in ',name_sub
            write(out_unitp,*) ' S(i,i) is too small',abs(Sii)
            STOP 'ERROR in UGPsi_spec: S(i,i) is too small'
          END IF
          C1(:,i) = C1(:,i) / sqrt(abs(Sii))
          Sii = dot_product(C1(:,i),matmul(S,C1(:,i)))
          DO j=i+1,n
            Sij = dot_product(C1(:,i),matmul(S,C1(:,j)))
            C1(:,j) = C1(:,j) - (Sij/Sii)*C1(:,i)
          END DO
        END DO
        Ho = matmul(transpose(conjg(C1)),matmul(H,C1))
        So = matmul(transpose(conjg(C1)),matmul(S,C1))

        IF (debug)  write(out_unitp,*) ' Max non-Hermitian Ho :',maxval(abs(Ho-conjg(transpose(Ho))))

        IF (debug .AND. n <= nmax) THEN

          write(out_unitp,*)
          write(out_unitp,*) 'Ho'
          write(out_unitp,*) 'Real(Ho)'
          CALL Write_Mat(real(Ho,kind=Rkind),out_unitp,6)
          write(out_unitp,*) 'im(Ho)'
          CALL Write_Mat(aimag(Ho),out_unitp,6)

          write(out_unitp,*)
          write(out_unitp,*) 'So: Identity matrix'
          write(out_unitp,*) 'Real(So)'
          CALL Write_Mat(real(So,kind=Rkind),out_unitp,6)
          write(out_unitp,*) 'im(So)'
          CALL Write_Mat(aimag(So),out_unitp,6)
          write(out_unitp,*)
        END IF


        DO i=1,n
          SoDiag(i) = So(i,i)
          So(i,i)   = abs(So(i,i)) - CONE
        END DO
        max_diff = maxval(abs(So))
        IF (debug)  write(out_unitp,*) ' Max_diff on (So-Idmatrix) :',max_diff

        IF (max_diff > ONETENTH**6) THEN
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) ' So is not the identity matrix. max_diff:',max_diff
          write(out_unitp,*)
          write(out_unitp,*) 'Real(So)'
          CALL Write_Mat(real(So,kind=Rkind),out_unitp,6)
          write(out_unitp,*) 'im(So)'
          CALL Write_Mat(aimag(So),out_unitp,6)
          write(out_unitp,*)
          STOP 'ERROR in UGPsi_spec: So is not the identity matrix.'
        END IF

        IF (allocated(Vec)) CALL dealloc_NParray(Vec,'Vec',name_sub)
        CALL alloc_NParray(Vec,[n,n],'Vec',name_sub)

        IF (allocated(Eig)) CALL dealloc_NParray(Eig,'Eig',name_sub)
        CALL alloc_NParray(Eig,[n],  'Eig',name_sub)

        IF (allocated(UPsiOnKrylov))                                    &
              CALL dealloc_NParray(UPsiOnKrylov,'UPsiOnKrylov',name_sub)
        CALL alloc_NParray(UPsiOnKrylov,[n],'UPsiOnKrylov',name_sub)

        CALL diagonalization(Ho,Eig,Vec,3,1,.TRUE.)
        IF (debug) THEN
          write(out_unitp,*) 'Eig',Eig

          write(out_unitp,*)
          write(out_unitp,*) 'Real(Vec)'
          CALL Write_Mat(real(Vec,kind=Rkind),out_unitp,6)
          write(out_unitp,*) 'im(Vec)'
          CALL Write_Mat(aimag(Vec),out_unitp,6)
          write(out_unitp,*)
        END IF

      END IF

      ! loop on the eigenvectors
      UPsiOnKrylov = CZERO

      DO i=1,n
        IF (debug) write(out_unitp,*) 'norm2 Vec(:,i)',i,dot_product(Vec(:,i),Vec(:,i))

        coef_i = Vec(1,i) ! just (1,i) because psi on the Krylov subspace is [1,0,0,...0]
                          ! that why we use Schmidt orthogonalisation
                          ! We have to multiply by So(i,i)=SoDiag(i), the values can be -1.
        coef_i = coef_i * exp(-EYE*Eig(i)*SoDiag(i)*deltaT) ! spectral propa

        UPsiOnKrylov(:) = UPsiOnKrylov(:) + conjg(Vec(:,i))*coef_i ! update U.psi on the Krylov space
      END DO
      IF (debug) write(out_unitp,*) 'norm2 UPsiOnKrylov',dot_product(UPsiOnKrylov,UPsiOnKrylov)

      UPspectral_n = UPsiOnKrylov(n)

      UPsiOnKrylov(:) = matmul(C1,UPsiOnKrylov)

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'abs(UPspectral_n)',abs(UPspectral_n)
        write(out_unitp,*) 'abs(UPsiOnKrylov)',abs(UPsiOnKrylov)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------


  END SUBROUTINE UGPsi_spec
  SUBROUTINE UGPsi_taylor(UPsiOnKrylov,UPspectral_n,H,S,Vec,Eig,deltaT,n,With_diago)
      USE mod_system
      IMPLICIT NONE



!----- variables for the WP propagation ----------------------------
      integer,                           intent(in)     :: n
      complex (kind=Rkind),              intent(in)     :: H(:,:)
      complex (kind=Rkind),              intent(in)     :: S(:,:)
      complex (kind=Rkind), allocatable, intent(inout)  :: Vec(:,:)
      real (kind=Rkind),    allocatable, intent(inout)  :: Eig(:)
      complex (kind=Rkind), allocatable, intent(inout)  :: UPsiOnKrylov(:)
      real (kind=Rkind),                 intent(in)     :: deltaT
      logical,                           intent(in)     :: With_diago
      complex (kind=Rkind),              intent(inout)  :: UPspectral_n
!------ working variables ---------------------------------
      complex (kind=Rkind) :: coef_i
      integer              :: i,j

      complex (kind=Rkind)     :: Ho(n,n)
      real (kind=Rkind)        :: E0
      complex (kind=Rkind)     :: So(n,n),SoDiag(n)
      complex (kind=Rkind)     :: C1(n,n)
      complex (kind=Rkind)     :: Sii,Sij
      real (kind=Rkind)        :: max_diff
      complex (kind=Rkind)     :: Vtemp(n)


      real (kind=Rkind)        :: tol_ortho = ONETENTH**10
      integer                  :: max_it_taylor = 1000
      real (kind=Rkind)        :: tol_taylor = ONETENTH**20
      real (kind=Rkind)        :: max_norm_UPsiOnKrylov
!----- for debuging --------------------------------------------------
      integer, parameter :: nmax = 12
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='UGPsi_taylor'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'deltaT',deltaT
        write(out_unitp,*) 'n',n
        write(out_unitp,*) 'With_diago',With_diago
        IF (With_diago .AND. n <= nmax) THEN
          write(out_unitp,*)
          write(out_unitp,*) 'H'
          write(out_unitp,*) 'Real(H)'
          CALL Write_Mat(real(H,kind=Rkind),out_unitp,6)
          write(out_unitp,*) 'im(H)'
          CALL Write_Mat(aimag(H),out_unitp,6)
          write(out_unitp,*)
          write(out_unitp,*) 'Real(S)'
          CALL Write_Mat(real(S,kind=Rkind),out_unitp,6)
          write(out_unitp,*) 'im(S)'
          CALL Write_Mat(aimag(S),out_unitp,6)
        END IF
      END IF
!-----------------------------------------------------------

      IF (With_diago) THEN
        IF (allocated(Vec)) CALL dealloc_NParray(Vec,'Vec',name_sub)
        CALL alloc_NParray(Vec,[n,n],'Vec',name_sub)

        IF (allocated(Eig)) CALL dealloc_NParray(Eig,'Eig',name_sub)
        CALL alloc_NParray(Eig,[n],  'Eig',name_sub)

        CALL diagonalization(S,Eig,Vec,3,1,.TRUE.)
        IF (debug) write(out_unitp,*) 'Eig of S',Eig
        IF (debug) write(out_unitp,*) 'min(abs(Eig)) of S',minval(abs(Eig))

        !first transfo S -> I (ortho Smidt)
        So = S
        C1 = CZERO
        DO i=1,n
          Sii = S(i,i)
          IF (abs(Sii) <tol_ortho) THEN
            write(out_unitp,*) 'ERROR in ',name_sub
            write(out_unitp,*) 'S(i,i) is too small',abs(Sii)
            STOP 'ERROR in UGPsi_spec: S(i,i) is too small'
          END IF
          C1(i,i) = ONE / sqrt(abs(Sii))
        END DO

        DO i=1,n
          Sii = dot_product(C1(:,i),matmul(S,C1(:,i)))
          IF (abs(Sii) < tol_ortho) THEN
            write(out_unitp,*) 'ERROR in ',name_sub
            write(out_unitp,*) ' S(i,i) is too small',abs(Sii)
            STOP 'ERROR in UGPsi_spec: S(i,i) is too small'
          END IF
          C1(:,i) = C1(:,i) / sqrt(abs(Sii))
          Sii = dot_product(C1(:,i),matmul(S,C1(:,i)))
          DO j=i+1,n
            Sij = dot_product(C1(:,i),matmul(S,C1(:,j)))
            C1(:,j) = C1(:,j) - (Sij/Sii)*C1(:,i)
          END DO
        END DO
        Ho = matmul(transpose(conjg(C1)),matmul(H,C1))
        So = matmul(transpose(conjg(C1)),matmul(S,C1))

        IF (debug .AND. n <= nmax) THEN

          write(out_unitp,*)
          write(out_unitp,*) 'Ho'
          write(out_unitp,*) 'Real(Ho)'
          CALL Write_Mat(real(Ho,kind=Rkind),out_unitp,6)
          write(out_unitp,*) 'im(Ho)'
          CALL Write_Mat(aimag(Ho),out_unitp,6)

          write(out_unitp,*)
          write(out_unitp,*) 'So: Identity matrix'
          write(out_unitp,*) 'Real(So)'
          CALL Write_Mat(real(So,kind=Rkind),out_unitp,6)
          write(out_unitp,*) 'im(So)'
          CALL Write_Mat(aimag(So),out_unitp,6)
          write(out_unitp,*)
        END IF


        DO i=1,n
          SoDiag(i) = So(i,i)
          So(i,i)   = abs(So(i,i)) - CONE
        END DO
        max_diff = maxval(abs(So))
        IF (debug)  write(out_unitp,*) ' Max_diff on (So-Idmatrix) :',max_diff

        IF (max_diff > ONETENTH**6) THEN
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) ' So is not the identity matrix. max_diff:',max_diff
          write(out_unitp,*)
          write(out_unitp,*) 'Real(So)'
          CALL Write_Mat(real(So,kind=Rkind),out_unitp,6)
          write(out_unitp,*) 'im(So)'
          CALL Write_Mat(aimag(So),out_unitp,6)
          write(out_unitp,*)
          STOP 'ERROR in UGPsi_spec: So is not the identity matrix.'
        END IF

        !change Ho to take into account the sign of the diagonal of So (we assume, So is diagonal matrix)
        !  -> Ho = So^-1 Ho
        DO i=1,n
          Ho(i,:) = Ho(i,:)/SoDiag(i)
        END DO
        E0 = real(H(1,1),kind=Rkind)
        Ho = -EYE*deltaT*(Ho-E0)

        IF (allocated(UPsiOnKrylov))                                    &
              CALL dealloc_NParray(UPsiOnKrylov,'UPsiOnKrylov',name_sub)
        CALL alloc_NParray(UPsiOnKrylov,[n],'UPsiOnKrylov',name_sub)

      END IF

      ! loop on the eigenvectors
      UPsiOnKrylov    = CZERO
      UPsiOnKrylov(1) = CONE
      Vtemp           = UPsiOnKrylov ! here we have the 0 order term.
      !write(out_unitp,*) '0 abs(UPsiOnKrylov)',abs(UPsiOnKrylov)
      !write(out_unitp,*) '0 Vtemp',abs(Vtemp)
      max_norm_UPsiOnKrylov = ONE

      DO i=1,max_it_taylor
        Vtemp = matmul(Ho,Vtemp) / cmplx(i,kind=Rkind)
        !write(out_unitp,*) i,' Vtemp',abs(Vtemp)
        UPsiOnKrylov = UPsiOnKrylov + Vtemp
        !write(out_unitp,*) i,'abs(UPsiOnKrylov)',abs(UPsiOnKrylov)

        IF (sqrt(abs(dot_product(UPsiOnKrylov,UPsiOnKrylov))) > max_norm_UPsiOnKrylov) &
           max_norm_UPsiOnKrylov = sqrt(abs(dot_product(UPsiOnKrylov,UPsiOnKrylov)))
        IF (max_norm_UPsiOnKrylov > TEN**10) THEN
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) ' UPsiOnKrylov is too large at the Taylor iteration:',i
          write(out_unitp,*) ' max_norm_UPsiOnKrylov',max_norm_UPsiOnKrylov
          write(out_unitp,*) 'abs(UPsiOnKrylov)',abs(UPsiOnKrylov)
          write(out_unitp,*)
          STOP ' ERROR in Taylor: UPsiOnKrylov too large'
        END IF

        IF (sqrt(abs(dot_product(Vtemp,Vtemp))) < tol_taylor) EXIT

      END DO
      IF (debug) write(out_unitp,*) 'Taylor order',i,sqrt(abs(dot_product(Vtemp,Vtemp)))
      IF (debug) write(out_unitp,*) 'max_norm_UPsiOnKrylov',max_norm_UPsiOnKrylov

      IF (debug) write(out_unitp,*) 'norm2 UPsiOnKrylov',dot_product(UPsiOnKrylov,UPsiOnKrylov)

      UPspectral_n = UPsiOnKrylov(n)*exp(EYE*deltaT*E0)

      UPsiOnKrylov(:) = matmul(C1,conjg(UPsiOnKrylov)) ! need the conjg ???

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'abs(UPspectral_n)',abs(UPspectral_n)
        write(out_unitp,*) 'abs(UPsiOnKrylov)',abs(UPsiOnKrylov)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------


  END SUBROUTINE UGPsi_taylor
  SUBROUTINE UPsi_spec_v1(UPsiOnKrylov,H,deltaT,n)
      USE mod_system
      IMPLICIT NONE



!----- variables for the WP propagation ----------------------------
      integer,               intent(in)     :: n
      complex (kind=Rkind),  intent(in)     :: H(n,n)
      complex (kind=Rkind),  intent(inout)  :: UPsiOnKrylov(n)
      real (kind=Rkind),     intent(in)     :: deltaT

!------ working variables ---------------------------------
      complex (kind=Rkind) :: Vec(n,n)
      real (kind=Rkind)    :: Eig(n)
      complex (kind=Rkind) :: coef_i
      integer              :: i

!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='UPsi_spec_v1'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'deltaT',deltaT
      END IF
!-----------------------------------------------------------

      CALL diagonalization(H,Eig,Vec,3,1,.TRUE.)

      ! loop on the eigenvectors
      UPsiOnKrylov = CZERO
      DO i=1,n
        coef_i = Vec(1,i) ! just (1,i) because psi on the Krylov space is [1,0,0,...0]
        coef_i = coef_i * exp(-EYE*Eig(i)*deltaT) ! spectral propa

        UPsiOnKrylov(:) = UPsiOnKrylov(:) + conjg(Vec(:,i))*coef_i ! update U.psi on the Krylov space

      END DO

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'abs(UPsiOnKrylov)',abs(UPsiOnKrylov)
        write(out_unitp,*) 'norm2 UPsiOnKrylov',dot_product(UPsiOnKrylov,UPsiOnKrylov)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------


  END SUBROUTINE UPsi_spec_v1
  SUBROUTINE march_SIP_v1(T,no,psi,psi0,para_H,para_propa)
      USE mod_system
      USE mod_psi,   ONLY : param_psi,ecri_psi,norm2_psi,renorm_psi,Overlap_psi1_psi2
      USE mod_Op,    ONLY : param_Op, sub_PsiOpPsi, sub_OpPsi,sub_scaledOpPsi
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)    :: para_H

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa),    intent(inout)  :: para_propa
      TYPE (param_psi),      intent(in)     :: psi0
      TYPE (param_psi),      intent(inout)  :: psi
      integer,               intent(in)     :: no
      real (kind=Rkind),     intent(in)     :: T

!------ working variables ---------------------------------
      TYPE (param_psi)          :: w1,w2
      complex (kind=Rkind)      :: cdot
      real (kind=Rkind)         :: E0,microT,microdeltaT,phase,microphase
      integer                   :: n


      complex (kind=Rkind) :: H(para_propa%para_poly%npoly,para_propa%para_poly%npoly)
      complex (kind=Rkind) :: UPsiOnKrylov(para_propa%para_poly%npoly)

      TYPE (param_psi), allocatable :: tab_KrylovSpace(:)
      complex (kind=Rkind) :: Overlap,ReNor,coef_i

      logical, parameter :: test = .FALSE.
      !logical, parameter :: test = .TRUE.

      integer              :: it,i,j

!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='march_SIP_v1'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'deltaT',para_propa%WPdeltaT
        write(out_unitp,*) 'E0,Esc',para_H%E0,para_H%Esc
        write(out_unitp,*) 'order',para_propa%para_poly%npoly
      END IF
!-----------------------------------------------------------

      allocate(tab_KrylovSpace(para_propa%para_poly%npoly+1))

      tab_KrylovSpace(1) = psi

      E0 = para_H%E0
      H(:,:) = CZERO
      DO j=2,para_propa%para_poly%npoly+1
        CALL sub_OpPsi(Psi  =tab_KrylovSpace(j-1),                      &
                       OpPsi=tab_KrylovSpace(j),para_Op=para_H)
        IF (j == 2) THEN
          CALL Overlap_psi1_psi2(Overlap,tab_KrylovSpace(1),tab_KrylovSpace(2))
          E0 = real(Overlap,kind=Rkind)
        END IF
        CALL sub_scaledOpPsi(Psi  =tab_KrylovSpace(j-1),                &
                             OpPsi=tab_KrylovSpace(j),E0=E0,Esc=ONE)

        !make part of the H matrix (related to the vector j-1)
        ! OpPsi of the vector j-1 is in tab_KrylovSpace(j).
        ! That why we need a loop util npoly+1.
        DO i=1,j-1
          CALL Overlap_psi1_psi2(Overlap,tab_KrylovSpace(i),tab_KrylovSpace(j))
          H(i,j-1) = Overlap
          H(j-1,i) = conjg(Overlap)
        END DO

        CALL UPsi_spec_v1(UPsiOnKrylov(1:j-1),H(1:j-1,1:j-1),para_propa%WPdeltaT,j-1)
        !write(out_unitp,*) j-1,'abs(UPsiOnKrylov(j-1)',abs(UPsiOnKrylov(j-1))
        IF (abs(UPsiOnKrylov(j-1)) < para_propa%para_poly%poly_tol .OR. &
            j == para_propa%para_poly%npoly+1) THEN
          n = j-1
          write(out_unitp,*) n,'abs(UPsiOnKrylov(n)',abs(UPsiOnKrylov(n))
          EXIT
        END IF

        ! ortho, the Krylov space: new vector
!        !1st version (not very stable)
!        CALL renorm_psi(tab_KrylovSpace(j))
!        DO i=1,j-1
!          CALL Overlap_psi1_psi2(Overlap,tab_KrylovSpace(j),tab_KrylovSpace(i))
!          IF (abs(Overlap) == ZERO) CYCLE
!          ReNor = ONE/sqrt(ONE-abs(Overlap)**2)
!          w1                 = tab_KrylovSpace(j)*ReNor
!          w2                 = tab_KrylovSpace(i) * (Overlap*ReNor)
!          tab_KrylovSpace(j) = w1-w2
!        END DO

        !2e version (more stable than the 1st)
        CALL renorm_psi(tab_KrylovSpace(j))
        DO i=1,j-1
          CALL Overlap_psi1_psi2(Overlap,tab_KrylovSpace(i),tab_KrylovSpace(j))
          IF (abs(Overlap) == ZERO) CYCLE
          tab_KrylovSpace(j) = tab_KrylovSpace(j) - tab_KrylovSpace(i) * Overlap
          CALL renorm_psi(tab_KrylovSpace(j))
        END DO

!        !3e version (less stable than the 2e)
!        DO i=1,j-1
!          CALL Overlap_psi1_psi2(Overlap,tab_KrylovSpace(j),tab_KrylovSpace(i))
!          IF (abs(Overlap) == ZERO) CYCLE
!          CALL norm2_psi(tab_KrylovSpace(j))
!          w1                 = tab_KrylovSpace(j)
!          w2                 = tab_KrylovSpace(i) * (conjg(Overlap)/tab_KrylovSpace(j)%norm2)
!          tab_KrylovSpace(j) = w1-w2
!        END DO
!        CALL renorm_psi(tab_KrylovSpace(j))

      END DO

      IF (debug) write(out_unitp,*) 'abs(UPsiOnKrylov)',abs(UPsiOnKrylov(1:n))

      Psi = ZERO
      DO i=1,n
        Psi = Psi + UPsiOnKrylov(i)*tab_KrylovSpace(i)
      END DO

      !- check norm ------------------
      CALL norm2_psi(psi)
      IF ( psi%norm2 > para_propa%max_norm2) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' STOP propagation: norm > max_norm',                &
                         psi%norm2
        para_propa%march_error   = .TRUE.
        para_propa%test_max_norm = .TRUE.
        STOP
      END IF

      !- Phase Shift -----------------
      phase = E0*para_propa%WPdeltaT
      psi = psi * exp(-EYE*phase)

!      microdeltaT = para_propa%WPdeltaT/real(para_propa%nb_micro,kind=Rkind)
!      microphase  = phase/real(para_propa%nb_micro,kind=Rkind)
!
!
!      !write(out_unitp,*) 'para_propa%nb_micro',para_propa%nb_micro
!      !write(out_unitp,*) 'microdeltaT',microdeltaT
!      !write(out_unitp,*) 'microphase',microphase
!
!      phase  = ZERO
!      microT = ZERO
!      DO it=1,para_propa%nb_micro
!
!        microT = microT + microdeltaT
!        phase = phase + microphase
!
!        rtj = cmplx(ONE,ZERO,kind=Rkind)
!        cdot = psi0Hkpsi0(0)
!        DO j=1,j_exit
!
!          rt_tmp =  cmplx(ZERO,-microT/real(j,kind=Rkind),kind=Rkind)
!          rtj = rt_tmp * rtj
!          cdot = cdot + psi0Hkpsi0(j) * rtj
!        END DO
!        cdot = cdot * cmplx( cos(phase),-sin(phase),kind=Rkind)
!        CALL Write_AutoCorr(no,T+microT,cdot)
!      END DO
!      flush(no)

      cdot = Calc_AutoCorr(psi0,psi,para_propa,T,Write_AC=.FALSE.)
      CALL Write_AutoCorr(no,T + para_propa%WPdeltaT,cdot)
      flush(no)


      ! deallocation
      DO i=1,size(tab_KrylovSpace)
         CALL dealloc_psi(tab_KrylovSpace(i))
      END DO
      deallocate(tab_KrylovSpace)
      CALL dealloc_psi(w1)
      CALL dealloc_psi(w2)



!-----------------------------------------------------------
      IF (debug) THEN
        CALL norm2_psi(psi)
        write(out_unitp,*) 'norm psi',psi%norm2
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      IF (debug) STOP


  END SUBROUTINE march_SIP_v1
  SUBROUTINE march_SIP_v0(T,no,psi,psi0,para_H,para_propa)
      USE mod_system
      USE mod_psi,   ONLY : param_psi,ecri_psi,norm2_psi,renorm_psi,Overlap_psi1_psi2
      USE mod_Op,    ONLY : param_Op, sub_PsiOpPsi, sub_OpPsi,sub_scaledOpPsi
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)     :: para_H

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa),    intent(in)     :: para_propa
      TYPE (param_psi),      intent(in)     :: psi0
      TYPE (param_psi),      intent(inout)  :: psi
      integer,               intent(in)     :: no
      real (kind=Rkind),     intent(in)     :: T

!------ working variables ---------------------------------
      TYPE (param_psi)          :: w1,w2
      complex (kind=Rkind)      :: cdot
      real (kind=Rkind)         :: E0,microT,microdeltaT,phase,microphase


      complex (kind=Rkind) :: H(para_propa%para_poly%npoly,para_propa%para_poly%npoly)
      complex (kind=Rkind) :: S(para_propa%para_poly%npoly,para_propa%para_poly%npoly)
      complex (kind=Rkind) :: Vec(para_propa%para_poly%npoly,para_propa%para_poly%npoly)
      real (kind=Rkind)    :: Eig(para_propa%para_poly%npoly)
      complex (kind=Rkind) :: PsiOnVec(para_propa%para_poly%npoly)
      complex (kind=Rkind) :: HPsiOnVec(para_propa%para_poly%npoly)
      complex (kind=Rkind) :: HPsiOnKrylov(para_propa%para_poly%npoly)

      TYPE (param_psi), allocatable :: tab_KrylovSpace(:)
      complex (kind=Rkind) :: Overlap,ReNor

      logical, parameter :: test = .FALSE.
      !logical, parameter :: test = .TRUE.

      integer              :: it,i,j

!----- for debuging --------------------------------------------------
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='march_SIP_v0'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'deltaT',para_propa%WPdeltaT
        write(out_unitp,*) 'E0,Esc',para_H%E0,para_H%Esc
        write(out_unitp,*) 'order',para_propa%para_poly%npoly
      END IF
!-----------------------------------------------------------

      allocate(tab_KrylovSpace(para_propa%para_poly%npoly+1))

      tab_KrylovSpace(1) = psi

      E0 = para_H%E0
      H(:,:) = CZERO
      DO j=2,para_propa%para_poly%npoly+1
        CALL sub_OpPsi(Psi  =tab_KrylovSpace(j-1),                      &
                       OpPsi=tab_KrylovSpace(j),para_Op=para_H)
        IF (j == 2) THEN
          CALL Overlap_psi1_psi2(Overlap,tab_KrylovSpace(1),tab_KrylovSpace(2))
          E0 = real(Overlap,kind=Rkind)
        END IF
        CALL sub_scaledOpPsi(Psi  =tab_KrylovSpace(j-1),                &
                             OpPsi=tab_KrylovSpace(j),E0=E0,Esc=ONE)

        !make part of the H matrix (related to the vector j-1)
        ! OpPsi of the vector j-1 is in tab_KrylovSpace(j).
        ! That why we need a loop util npoly+1.
        DO i=1,j-1
          CALL Overlap_psi1_psi2(Overlap,tab_KrylovSpace(i),tab_KrylovSpace(j))
          H(i,j-1) = Overlap
          H(j-1,i) = conjg(Overlap)
        END DO

        IF (j == para_propa%para_poly%npoly+1) EXIT
        ! ortho, the Krylov space: new vector
        CALL renorm_psi(tab_KrylovSpace(j))
!        DO i=1,j-1
!          CALL Overlap_psi1_psi2(Overlap,tab_KrylovSpace(j),tab_KrylovSpace(i))
!          IF (abs(Overlap) == ZERO) CYCLE
!          ReNor = ONE/sqrt(ONE-abs(Overlap)**2)
!          w1                 = tab_KrylovSpace(j)*ReNor
!          w2                 = tab_KrylovSpace(i) * (Overlap*ReNor)
!          tab_KrylovSpace(j) = w1-w2
!        END DO
        DO i=1,j-1
          CALL Overlap_psi1_psi2(Overlap,tab_KrylovSpace(j),tab_KrylovSpace(i))
          IF (abs(Overlap) == ZERO) CYCLE
          !ReNor = ONE/sqrt(ONE-abs(Overlap)**2)
          w1                 = tab_KrylovSpace(j)
          w2                 = tab_KrylovSpace(i) * conjg(Overlap)
          tab_KrylovSpace(j) = w1-w2
          CALL renorm_psi(tab_KrylovSpace(j))

        END DO

      END DO

      IF (debug) THEN
        write(out_unitp,*) 'H'
        CALL Write_Mat(H,out_unitp,6)
      END IF

      CALL diagonalization(H,Eig,Vec,3,1,.TRUE.)

      IF (debug) THEN
        write(out_unitp,*) 'Vec'
        CALL Write_Mat(Vec,out_unitp,6)
        write(out_unitp,*) 'Eig',Eig
      END IF

      IF (test) THEN
        !overlapp matrix S(i,j)
        DO j=1,para_propa%para_poly%npoly
        DO i=1,para_propa%para_poly%npoly
          CALL Overlap_psi1_psi2(Overlap,tab_KrylovSpace(i),tab_KrylovSpace(j))
          S(i,j) = Overlap
        END DO
        END DO
        write(out_unitp,*) 'S'
        CALL Write_Mat(S,out_unitp,6)

        write(out_unitp,*) 'H'
        CALL Write_Mat(H,out_unitp,6)

        write(out_unitp,*) 'Eig1',dot_product(Vec(:,1),matmul(H,Vec(:,1)))

        H = CZERO
        DO i=1,para_propa%para_poly%npoly
          H(i,i) = Eig(i)
        END DO
        H = matmul(conjg(Vec),matmul(H,transpose(Vec)))
        !H = matmul(transpose(conjg(Vec)),matmul(H,Vec)) ! wrong expression

        write(out_unitp,*) 'H'
        CALL Write_Mat(H,out_unitp,6)
        STOP
      END IF

      ! since the first vector of the orthonormalized Krylov space is psi
      HPsiOnKrylov(:) = CZERO  ; HPsiOnKrylov(1) = CONE
      PsiOnVec    = matmul(transpose(Vec),HPsiOnKrylov)
      !DO i=1,para_propa%para_poly%npoly
      !  PsiOnVec(i) = Vec(1,i)
      !END DO
      IF (debug) write(out_unitp,*) 'PsiOnVec',PsiOnVec

      HPsiOnVec   = exp(-EYE*Eig*para_propa%WPdeltaT)*PsiOnVec
      IF (debug) write(out_unitp,*) 'HPsiOnVec',HPsiOnVec


      HPsiOnKrylov  = matmul(conjg(Vec),HPsiOnVec)
      !HPsiOnKrylov = CZERO
      !DO i=1,para_propa%para_poly%npoly
      !  HPsiOnKrylov(:) = HPsiOnKrylov(:) + conjg(Vec(:,i))*HPsiOnVec(i)
      !END DO
      IF (debug) write(out_unitp,*) 'HPsiOnKrylov',HPsiOnKrylov
      IF (debug) write(out_unitp,*) 'abs(HPsiOnKrylov)',abs(HPsiOnKrylov)

      Psi = ZERO
      DO i=1,para_propa%para_poly%npoly
        Psi = Psi + HPsiOnKrylov(i)*tab_KrylovSpace(i)
      END DO

      ! deallocation
      DO i=1,size(tab_KrylovSpace)
         CALL dealloc_psi(tab_KrylovSpace(i))
      END DO
      deallocate(tab_KrylovSpace)
      CALL dealloc_psi(w1)
      CALL dealloc_psi(w2)

      !- Phase Shift -----------------
      phase = E0*para_propa%WPdeltaT
      psi = psi * exp(-EYE*phase)

      cdot = Calc_AutoCorr(psi0,psi,para_propa,T,Write_AC=.FALSE.)
      CALL Write_AutoCorr(no,T + para_propa%WPdeltaT,cdot)
      flush(no)


!-----------------------------------------------------------
      IF (debug) THEN
        CALL norm2_psi(psi)
        write(out_unitp,*) 'norm psi',psi%norm2
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      IF (debug) STOP


  END SUBROUTINE march_SIP_v0
!================================================================
!
!     march with spectral representation
!
!================================================================
      !!@description: march Spectral
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
  SUBROUTINE march_Spectral(T,no,psi,psi0,para_H,para_propa)
      USE mod_system
      USE mod_psi,   ONLY : param_psi,ecri_psi,norm2_psi,renorm_psi,Overlap_psi1_psi2
      USE mod_Op,    ONLY : param_Op, sub_PsiOpPsi, sub_OpPsi,sub_scaledOpPsi
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)     :: para_H

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa),    intent(in)     :: para_propa
      TYPE (param_psi),      intent(in)     :: psi0
      TYPE (param_psi),      intent(inout)  :: psi
      integer,               intent(in)     :: no
      real (kind=Rkind),     intent(in)     :: T

!------ working variables ---------------------------------
      !logical :: WithSpectralWP = .FALSE.
      logical :: WithSpectralWP = .TRUE.

      !TYPE (param_psi), pointer :: w1,w2
      complex (kind=Rkind) :: cdot
      real (kind=Rkind)    :: microT,microdeltaT,phase,microphase



      complex (kind=Rkind), allocatable :: PsiOnSpectral(:)
      complex (kind=Rkind), allocatable :: HPsiOnSpectral(:)

      integer              :: it,i,j

!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='march_Spectral'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'deltaT',para_propa%WPdeltaT
        write(out_unitp,*) 'E0,Esc',para_H%E0,para_H%Esc
        write(out_unitp,*) 'order',para_propa%para_poly%npoly
      END IF
!-----------------------------------------------------------

      IF (.NOT. para_H%spectral_done .AND. para_H%spectral_Op /= para_H%n_Op) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' H is not in a spectral representation or'
        write(out_unitp,*) ' The spectral representation not done on H.'
        write(out_unitp,*) ' spectral_done',para_H%spectral_done
        write(out_unitp,*) ' spectral_Op',para_H%spectral_Op
        STOP
      END IF

      microdeltaT = para_propa%WPdeltaT/                                &
                    real(para_propa%nb_micro,kind=Rkind)

      CALL alloc_NParray(PsiOnSpectral, [para_H%nb_tot],'PsiOnSpectral', name_sub)
      CALL alloc_NParray(HPsiOnSpectral,[para_H%nb_tot],'HPsiOnSpectral',name_sub)
      IF (WithSpectralWP) THEN

        microT = ZERO
        DO it=1,para_propa%nb_micro
          microT = microT + microdeltaT

          IF (para_H%cplx) THEN
            ! 2d: propagation with spectral representation
            Psi%CvecB   = exp(-EYE*para_H%Cdiag*microdeltaT)*Psi%CvecB
          ELSE
            ! 2d: propagation with spectral representation
            ! warning: para_H%Rdiag does not exist when "direct=4"
            Psi%CvecB   = exp(-EYE*para_H%Rdiag*microdeltaT)*Psi%CvecB
          END IF

          cdot = Calc_AutoCorr(psi0,psi,para_propa,T+microT,Write_AC=.FALSE.)
          CALL Write_AutoCorr(no,T+microT,cdot)
          flush(no)

        END DO
      ELSE

        IF (para_H%cplx) THEN
          ! 1st: project psi on the spectral basis
          PsiOnSpectral(:) = matmul(transpose(para_H%Cvp),Psi%CvecB)

          ! 2d: propagation with spectral representation
          HPsiOnSpectral   = exp(-EYE*para_H%Cdiag*para_propa%WPdeltaT)*PsiOnSpectral

          ! 3d: back to the initial basis
          Psi%CvecB = matmul(conjg(para_H%Cvp),HPsiOnSpectral)

        ELSE
          ! 1st: project psi on the spectral basis
          PsiOnSpectral(:) = matmul(transpose(para_H%Rvp),Psi%CvecB)

          ! 2d: propagation with spectral representation
          HPsiOnSpectral   = exp(-EYE*para_H%Rdiag*para_propa%WPdeltaT)*PsiOnSpectral

          ! 3d: back to the initial basis
          Psi%CvecB = matmul(para_H%Rvp,HPsiOnSpectral)

          !write(6,*) 'E',sum(conjg(HPsiOnSpectral)*HPsiOnSpectral*para_H%Rdiag)

        END IF

        CALL dealloc_NParray(PsiOnSpectral, 'PsiOnSpectral', name_sub)
        CALL dealloc_NParray(HPsiOnSpectral,'HPsiOnSpectral',name_sub)
        cdot = Calc_AutoCorr(psi0,psi,para_propa,T,Write_AC=.FALSE.)
        CALL Write_AutoCorr(no,T + para_propa%WPdeltaT,cdot)
        flush(no)
      END IF




!-----------------------------------------------------------
      IF (debug) THEN
        CALL norm2_psi(psi)
        write(out_unitp,*) 'norm psi',psi%norm2
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------
!STOP

 END SUBROUTINE march_Spectral

!=======================================================================================
!     march cheby
!=======================================================================================
      !!@description: march cheby
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      SUBROUTINE march_cheby(T,no,psi,psi0,para_H,para_propa)
      USE mod_system
      USE mod_psi,   ONLY : param_psi,ecri_psi,norm2_psi
      USE mod_Op,    ONLY : param_Op, sub_OpPsi,sub_scaledOpPsi
      USE mod_propa
      USE mod_MPI_aux
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_psi)   :: psi0,psi

!------ working variables ---------------------------------
      TYPE (param_psi)          :: w1,w2,w3
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
      Logical              :: exitall ! for MPI
      Logical              :: Type_FileGrid4

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
      exitall=.FALSE. !< MPI: for controling the exit of all threads
      Type_FileGrid4=(para_H%para_ReadOp%para_FileGrid%Type_FileGrid==4)

      IF(keep_MPI) psi_save = psi

      DO icheb=1,max_cheby

        !IF(openmpi) THEN
        IF(Type_FileGrid4) THEN
          IF(keep_MPI) w1  = psi
          IF(para_propa%once_Hmin) THEN
            CALL sub_OpPsi(w1,w2,para_H) ! calculate once for Hmax
            CALL sub_Hmax(para_propa,para_H)

            para_propa%once_Hmin=.FALSE.
            para_propa%Hmax = para_propa%Hmax + para_propa%para_poly%DHmax
            para_propa%para_poly%Hmin = para_propa%Hmin
            para_propa%para_poly%Hmax = para_propa%Hmax

            CALL initialisation1_poly(para_propa%para_poly,                            &
                                      para_propa%WPdeltaT,                             &
                                      para_propa%type_WPpropa)

            para_H%scaled = .TRUE.
            para_H%E0     = para_propa%para_poly%E0
            para_H%Esc    = para_propa%para_poly%Esc

            !write(out_unitp,*) 'Hmin,Hmax check:', para_propa%Hmin,para_propa%Hmax
          ENDIF
        ENDIF

        para_propa%march_error      = .FALSE.
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

        ! calculate cofficients for Chebychev in para_propa%para_poly%coef_poly
        !r = max(ONE,HALF * para_propa%para_poly%deltaE * para_propa%WPdeltaT)
        ! When R=ONE the calculation is wrong, because deltaE should changed
        r = HALF * para_propa%para_poly%deltaE * para_propa%WPdeltaT
        CALL cof(r,para_propa%para_poly%npoly,para_propa%para_poly%coef_poly)

        IF (debug) THEN
          write(out_unitp,*) 'r,deltaE,WPdeltaT',r,para_propa%para_poly%deltaE,para_propa%WPdeltaT
          write(out_unitp,*) 'npoly',para_propa%para_poly%npoly
          !write(out_unitp,*) 'coef_poly',para_propa%para_poly%coef_poly(1:para_propa%para_poly%npoly)
        END IF

        psi0Hkpsi0(0) = Calc_AutoCorr(psi0,psi,para_propa,T,Write_AC=.FALSE.)
        IF(openmpi .AND. MPI_scheme/=1) CALL MPI_Bcast_(psi0Hkpsi0(0),size1_MPI,root_MPI)

        rt  = cmplx(ZERO,-ONE,kind=Rkind)
        rt2 = cmplx(ZERO,-TWO,kind=Rkind)

!     - The first term of the expansion ------------------
        ! IF(.NOT. openmpi) w1  = psi
        IF(.NOT. Type_FileGrid4) w1  = psi
        IF(keep_MPI) psi = psi * para_propa%para_poly%coef_poly(1)

        psi0Hkpsi0(1) =  psi0Hkpsi0(0)

!     - The second term of the expansion -----------------
        write(out_unitp,'(a)',advance='no') 'cheby rec:'

        !IF(.NOT. para_propa%once_Hmin .OR. .NOT. openmpi) CALL sub_OpPsi(w1,w2,para_H)
!        IF((Type_FileGrid4 .AND. .NOT. para_propa%once_Hmin) .OR. (.NOT. Type_FileGrid4)) THEN
!          CALL sub_OpPsi(w1,w2,para_H)
!        ENDIF
        IF(Type_FileGrid4) THEN
          IF(.NOT. para_propa%once_Hmin) CALL sub_OpPsi(w1,w2,para_H)
        ELSE
          CALL sub_OpPsi(w1,w2,para_H)
        ENDIF
        CALL sub_scaledOpPsi(w1,w2,para_H%E0,para_H%Esc) ! limited to MPI_id==0 in subroutine

        IF(keep_MPI) THEN
          w2  = w2 * rt
          psi = psi + w2 * para_propa%para_poly%coef_poly(2)

          psi0Hkpsi0(2) = Calc_AutoCorr(psi0,w2,para_propa,T,Write_AC=.FALSE.)
        ENDIF

!     - The higher terms of the expansion ----------------

        DO jt=3,para_propa%para_poly%npoly

          CALL sub_OpPsi(w2,w3,para_H)
          CALL sub_scaledOpPsi(w2,w3,para_H%E0,para_H%Esc)
          IF (mod(jt,100) == 0) write(out_unitp,'(a)',advance='no') '.'

          IF(keep_MPI) THEN
!           Recurrence relations of the Chebychev expansion:
            w3 = w1 + w3 * rt2
            w1 = w2
            w2 = w3
            psi = psi + w2 * para_propa%para_poly%coef_poly(jt)
          ENDIF

          IF(keep_MPI) THEN
            CALL norm2_psi(w2)
            norm_exit = abs(w2%norm2*para_propa%para_poly%coef_poly(jt))
          ENDIF
          IF(openmpi) CALL MPI_Bcast_(norm_exit,size1_MPI,root_MPI)

          jt_exit = jt
          IF (debug) write(out_unitp,*) 'jt,norms',jt,norm_exit

          IF (norm_exit > TEN**15) THEN
            write(out_unitp,*) ' WARNING: Norm^2 of the vector is TOO large (> 10^15)',jt,norm_exit
            IF(keep_MPI) para_propa%march_error = .TRUE.
            EXIT
          END IF

          IF(keep_MPI) THEN
            psi0Hkpsi0(jt) = Calc_AutoCorr(psi0,w2,para_propa,T,Write_AC=.FALSE.)
!#if(run_MPI)
!            IF (norm_exit < para_propa%para_poly%poly_tol) exitall=.TRUE.
!#else
!            IF (norm_exit < para_propa%para_poly%poly_tol) EXIT
!#endif
            IF(openmpi) THEN
              IF (norm_exit < para_propa%para_poly%poly_tol) exitall=.TRUE.
            ELSE
              IF (norm_exit < para_propa%para_poly%poly_tol) EXIT
            ENDIF
          ENDIF

          !> MPI: the other threads are waiting for master here
          ! CALL MPI_Bcast(exitall,size1_MPI,Real_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
          IF(openmpi .AND. MPI_scheme/=1) THEN
            CALL MPI_Bcast_(exitall,size1_MPI,root_MPI)
            CALL MPI_Bcast_(w2%norm2,size1_MPI,root_MPI)
          ENDIF
          IF(exitall) EXIT

        END DO ! for jt=3,para_propa%para_poly%npoly

        write(out_unitp,*) 'jt_exit,norms',jt_exit,abs(w2%norm2),norm_exit
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

          IF(keep_MPI) psi = psi_save
        ELSE
          EXIT
        END IF

      END DO ! for max_cheby

      IF (para_propa%march_error) THEN
        write(out_unitp,*) ' ERROR in march_cheby'
        write(out_unitp,*) ' It cannot converge ...'
        write(out_unitp,*) ' => Reduce the time step !!'
        STOP
      END IF

      CALL dealloc_psi(psi_save,delete_all=.TRUE.)
      CALL dealloc_psi(w1)
      CALL dealloc_psi(w2)
      CALL dealloc_psi(w3)

!     - Phase Shift -----------------------------------
      IF(keep_MPI) THEN
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
        IF(MPI_id==0) CALL Write_AutoCorr(no,T+microT,cdot)
      END DO
      IF(MPI_id==0) flush(no)
      ENDIF ! keep_MPI

!-----------------------------------------------------------
      IF (debug) THEN
        !write(out_unitp,*) 'psi'
        !CALL ecri_psi(T=T+para_propa%WPdeltaT,psi=psi,ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)
        write(out_unitp,*) 'END march_cheby'
      END IF
!-----------------------------------------------------------

      END SUBROUTINE march_cheby
!=======================================================================================

      SUBROUTINE march_cheby_old(T,no,psi,psi0,para_H,para_propa)
      USE mod_system
      USE mod_psi,   ONLY : param_psi,ecri_psi,norm2_psi
      USE mod_Op,    ONLY : param_Op, sub_OpPsi,sub_scaledOpPsi
      USE mod_MPI_aux
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_psi)   :: psi0,psi

!------ working variables ---------------------------------
      TYPE (param_psi)     :: w1,w2,w3

      complex (kind=Rkind) :: cdot
      complex (kind=Rkind) :: psi0Hkpsi0(0:para_propa%para_poly%npoly)
      real (kind=Rkind)    :: microT,T,microdeltaT,phase,microphase
      integer              :: no

      complex (kind=Rkind) :: rt,rt2,rtj,rt_tmp
      integer              :: it,j,jt,jt_exit
      real (kind=Rkind)    :: r,rg
      integer              :: i,i_qaie_corr
      real (kind=Rkind)    :: norm_exit
      Logical              :: Type_FileGrid4


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
      Type_FileGrid4=(para_H%para_ReadOp%para_FileGrid%Type_FileGrid==4)

      !IF(openmpi) THEN
      IF(Type_FileGrid4) THEN
        w1=psi
        IF(para_propa%once_Hmin) THEN
          CALL sub_OpPsi(w1,w2,para_H)
          CALL sub_Hmax(para_propa,para_H)
          para_propa%once_Hmin=.FALSE.

          para_propa%Hmax = para_propa%Hmax + para_propa%para_poly%DHmax
          para_propa%para_poly%Hmin = para_propa%Hmin
          para_propa%para_poly%Hmax = para_propa%Hmax

          CALL initialisation1_poly(para_propa%para_poly,                              &
                                    para_propa%WPdeltaT,                               &
                                    para_propa%type_WPpropa)

          para_H%scaled = .TRUE.
          para_H%E0     = para_propa%para_poly%E0
          para_H%Esc    = para_propa%para_poly%Esc

          !write(out_unitp,*) 'Hmin,Hmax check:', para_propa%Hmin,para_propa%Hmax
        ENDIF
      ENDIF ! direct4

      IF(keep_MPI) psi0Hkpsi0(:) = cmplx(ZERO,ZERO,kind=Rkind)

      r = HALF * para_propa%para_poly%deltaE * para_propa%WPdeltaT
      CALL cof(r,para_propa%para_poly%npoly,                            &
                 para_propa%para_poly%coef_poly)
      !write(out_unitp,*) 'r,deltaE,WPdeltaT',r,para_propa%para_poly%deltaE,para_propa%WPdeltaT
      !write(out_unitp,*) 'npoly,coef_poly',para_propa%para_poly%npoly,  &
      !  para_propa%para_poly%coef_poly(1:para_propa%para_poly%npoly)

      IF(keep_MPI) psi0Hkpsi0(0) = Calc_AutoCorr(psi0,psi,para_propa,T,Write_AC=.FALSE.)

!#if(run_MPI)
!      CALL MPI_Bcast(psi0Hkpsi0(0),size1_MPI,Cplx_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
!#endif
      IF(openmpi .AND. MPI_scheme/=1) CALL MPI_Bcast_(psi0Hkpsi0(0),size1_MPI,root_MPI)

      rt  = cmplx(ZERO,-ONE,kind=Rkind)
      rt2 = cmplx(ZERO,-TWO,kind=Rkind)

!     - The first term of the expansion ------------------
      IF(keep_MPI) THEN
        !IF(.NOT. openmpi) w1  = psi
        IF(.NOT. Type_FileGrid4) w1  = psi
        psi = psi * para_propa%para_poly%coef_poly(1)
        psi0Hkpsi0(1) =  psi0Hkpsi0(0)
      ENDIF

!     - The second term of the expansion -----------------
      write(out_unitp,'(a)',advance='no') 'cheby rec:'

      !IF(.NOT. para_propa%once_Hmin .OR. .NOT. openmpi) CALL sub_OpPsi(w1,w2,para_H)
!      IF((Type_FileGrid4 .AND. .NOT. para_propa%once_Hmin) .OR. (.NOT. Type_FileGrid4)) THEN
!        CALL sub_OpPsi(w1,w2,para_H)
!      ENDIF
      IF(Type_FileGrid4) THEN
        IF(.NOT. para_propa%once_Hmin) CALL sub_OpPsi(w1,w2,para_H)
      ELSE
        CALL sub_OpPsi(w1,w2,para_H)
      ENDIF
      CALL sub_scaledOpPsi(w1,w2,para_H%E0,para_H%Esc)

      IF(keep_MPI) THEN
        w2  = w2 * rt
        psi = psi + w2 * para_propa%para_poly%coef_poly(2)

        psi0Hkpsi0(2) = Calc_AutoCorr(psi0,w2,para_propa,T,Write_AC=.FALSE.)
      ENDIF

!     - The higher terms of the expansion ----------------

      DO jt=3,para_propa%para_poly%npoly

        CALL sub_OpPsi(w2,w3,para_H)
        CALL sub_scaledOpPsi(w2,w3,para_H%E0,para_H%Esc)
        IF (mod(jt,100) == 0) write(out_unitp,'(a)',advance='no') '.'

        IF(keep_MPI) THEN
!        Recurrence relations of the Chebychev expansion:
         w3 = w1 + w3 * rt2

         w1 = w2
         w2 = w3

         psi = psi + w2 * para_propa%para_poly%coef_poly(jt)

         CALL norm2_psi(w2)

         norm_exit = abs(w2%norm2*para_propa%para_poly%coef_poly(jt))
         jt_exit = jt
         IF (debug) write(out_unitp,*) 'jt,norms',jt,norm_exit

          IF (norm_exit > TEN**15) THEN
            write(out_unitp,*) ' ERROR in march_cheby'
            write(out_unitp,*) ' Norm of the vector is TOO large (> 10^15)',    &
                                                jt,norm_exit
           write(out_unitp,*) ' => Reduce the time step !!'
           STOP
         END IF
         psi0Hkpsi0(jt) = Calc_AutoCorr(psi0,w2,para_propa,T,Write_AC=.FALSE.)
        ENDIF ! for keep_MPI

        IF(openmpi .AND. MPI_scheme/=1) CALL MPI_Bcast_(norm_exit,size1_MPI,root_MPI)

        IF (norm_exit < para_propa%para_poly%poly_tol) EXIT

      END DO
      write(out_unitp,*) 'jt_exit,norms',jt_exit,abs(w2%norm2),norm_exit
      IF (debug) write(out_unitp,*) 'psi0Hkpsi0',psi0Hkpsi0(0:jt_exit)

      IF (norm_exit > para_propa%para_poly%poly_tol) THEN
        write(out_unitp,*) ' ERROR in march_cheby'
        write(out_unitp,*) ' Norm of the last vector TOO large',norm_exit
        write(out_unitp,*) ' poly_tol: ',para_propa%para_poly%poly_tol
        write(out_unitp,*) ' => npoly TOO small',para_propa%para_poly%npoly
        STOP
      END IF

!     - Phase Shift -----------------------------------
      IF(keep_MPI) THEN
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
        IF(MPI_id==0) CALL Write_AutoCorr(no,T+microT,cdot)
      END DO
      IF(MPI_id==0) flush(no)
      ENDIF ! keep_MPI

      CALL dealloc_psi(w1)
      CALL dealloc_psi(w2)
      CALL dealloc_psi(w3)


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
      USE mod_psi,   ONLY : param_psi,ecri_psi,norm2_psi
      USE mod_Op,    ONLY : param_Op, sub_OpPsi,sub_scaledOpPsi
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

        write(out_unitp,'(a)',advance='no') 'nOD_im rec:'

      w1 = psi
      DO j=1,para_propa%para_poly%max_poly
        IF (mod(j,10) == 0) write(out_unitp,'(a)',advance='no') '.'

        rtj =  signDT*para_propa%WPdeltaT/real(j,kind=Rkind)

        CALL sub_OpPsi(w1,w2,para_H)
        CALL sub_scaledOpPsi(w1,w2,para_H%E0,ONE)

        w1  = w2 * rtj
        psi = psi + w1

        CALL norm2_psi(w1)
        CALL norm2_psi(psi)
        !write(out_unitp,*) 'j,wi%n/psi%n',j,w1%norm2/psi%norm2
        IF (w1%norm2/psi%norm2 < para_propa%para_poly%poly_tol) EXIT

      END DO
      write(out_unitp,*) 'jt_exit,rap norm',j,w1%norm2/psi%norm2

      IF (para_propa%write_iter .OR. debug)                             &
                  write(out_unitp,*) 'j,wi%n/psi%n',j,w1%norm2/psi%norm2

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
      USE mod_psi,   ONLY : param_psi,ecri_psi,dealloc_psi,Overlap_psi1_psi2,norm2_psi
      USE mod_Op,    ONLY : param_Op, sub_OpPsi,sub_scaledOpPsi
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
      real (kind=Rkind)      :: DeltaE,Deltapsi,epsi,sign,norm2g
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
      norm2g = sqrt(g%norm2)

      H2psi = psi + Hpsi * (-DeltaT)
      psi = H2psi * (ONE/sqrt(S))

      Deltapsi = -A/S * DeltaT**2


      !- propagation ------------------------------------------
      IF (para_propa%write_iter .OR. debug) THEN
        write(out_unitp,*) 'WP it sqrt(norm2 g)',it,                &
            sqrt(g%norm2)*get_Conv_au_TO_unit('E','cm-1')
      END IF

      IF (debug) THEN
        write(out_unitp,*) 'WP it, E',it,                               &
                                  Ene0*get_Conv_au_TO_unit('E','cm-1'), &
                                DeltaE*get_Conv_au_TO_unit('E','cm-1')
        write(out_unitp,*) 'WP it sqrt(norm2 g)',it,                &
            sqrt(g%norm2)*get_Conv_au_TO_unit('E','cm-1')
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
      USE mod_psi,   ONLY : param_psi,ecri_psi,norm2_psi
      USE mod_field, ONLY : param_field,sub_dnE
      USE mod_Op,    ONLY : param_Op, sub_OpPsi,sub_scaledOpPsi
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H
      TYPE (param_Op)   :: para_Dip(3)

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_field) :: para_field

      integer          :: nb_WP
      TYPE (param_psi) :: WP(nb_WP)

!----- for printing --------------------------------------------------
      logical ::print_Op



!------ working parameters --------------------------------
      TYPE (param_psi)              :: w1,w2
      TYPE (param_psi) ,allocatable :: work_WP(:)

      real (kind=Rkind)      :: fac_k,phase
      complex (kind=Rkind)   :: Overlap,cdot,rt,rti,norm,coef,coef_ip(3)
      real (kind=Rkind)      :: limit

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
          write(out_unitp,*) 'norm WP',i,WP(i)%norm2

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
      max_der = 0
      nb_der = 3*para_propa%para_poly%npoly
      IF (para_field%max_der >= 0)                                      &
           nb_der = min(para_field%max_der,para_propa%para_poly%npoly-1)

      CALL alloc_NParray(tab_dnE,[nb_der,3],                            &
                      "tab_dnE","march_new0_noD_field",[0,1])
      tab_dnE(:,:) = ZERO

      DO k=0,nb_der
        CALL sub_dnE(dnE,k,T,para_field)
        tab_dnE(k,:) = -dnE(:)
      END DO
      max_der = nb_der
      write(out_unitp,*) 'max_der,npoly',max_der,para_propa%para_poly%npoly

      CALL alloc_NParray(work_WP,[para_propa%para_poly%npoly-1],        &
                        'work_WP',name_sub,[0])
!-----------------------------------------------------------
      DO j=1,nb_WP
!       Derivative calculation
        fac_k = ONE
        DO k=0,para_propa%para_poly%npoly-1

          IF (k == 0) THEN
            work_WP(k) = WP(j)
          ELSE
!           initialization with the H contribution
            CALL sub_OpPsi(work_WP(k-1),w1,para_H)
            CALL sub_scaledOpPsi(work_WP(k-1),w1,para_H%E0,ONE)
!           ici pas de binome car m=0
            coef = -EYE
            work_WP(k) = w1 * coef
            DO m=0,k-1
!             Dip contributions
              DO ip=1,3
                IF (.NOT. para_field%pola_xyz(ip)) CYCLE
                coef = -EYE*binomial(k-1,m) * tab_dnE(m,ip)
                CALL sub_OpPsi(work_WP(k-1-m),w1,para_Dip(ip))
                work_WP(k) = work_WP(k) + coef * w1
              END DO
            END DO
          END IF


          IF (taylor) THEN
!           Just the Taylor serie
            IF (k == 0) CYCLE ! because WP(j) is already the zero-order
            fac_k = fac_k * para_propa%WPdeltaT / real(k,kind=Rkind)
            w1 = work_WP(k) * fac_k
          ELSE
!           Infinite sum of the derivative of the field
            CALL sub_OpPsi(work_WP(k),w2,para_H)
            CALL sub_scaledOpPsi(work_WP(k),w2,para_H%E0,ONE)
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
              CALL sub_OpPsi(work_WP(k),w2,para_Dip(ip))
              w1 = w1 + w2 * coef_ip(ip)
            END DO
          END IF


          WP(j) = WP(j) + w1 * (-EYE)
          CALL norm2_psi(w1)
          write(out_unitp,*) 'norm2 dkpsi k',k,w1%norm2
          IF (w1%norm2 < para_propa%para_poly%poly_tol) EXIT

          IF (w1%norm2 > TEN**15) THEN
             write(out_unitp,*) ' ERROR in ',name_sub
             write(out_unitp,*) ' The norm2 of dnpsi(k) is huge !',w1%norm2
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
        IF ( WP(j)%norm2 > para_propa%max_norm2) THEN
           T  = T + para_propa%WPdeltaT
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) ' STOP propagation: norm > max_norm',             &
                         WP(j)%norm2
           para_propa%march_error   = .TRUE.
           para_propa%test_max_norm = .TRUE.
           STOP
        END IF

      END DO

      CALL dealloc_NParray(tab_dnE,"tab_dnE",name_sub)
      CALL dealloc_psi(w1)
      CALL dealloc_psi(w2)
      CALL dealloc_NParray(work_WP,'work_WP',name_sub)

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
      USE mod_psi,   ONLY : param_psi,ecri_psi,norm2_psi
      USE mod_field, ONLY : param_field,sub_dnE
      USE mod_Op,    ONLY : param_Op, sub_OpPsi,sub_scaledOpPsi
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H
      TYPE (param_Op)   :: para_Dip(3)

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_field) :: para_field

      integer          :: nb_WP
      TYPE (param_psi) :: WP(nb_WP)


!----- for printing --------------------------------------------------
      logical ::print_Op



!------ working parameters --------------------------------
      TYPE (param_psi)              :: w1,w2
      TYPE (param_psi), allocatable :: work_WP(:)
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
          write(out_unitp,*) 'norm WP',i,WP(i)%norm2

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

      max_der = 0
      nb_der = para_propa%para_poly%npoly-1
      IF (para_field%max_der >= 0)                                      &
           nb_der = min(para_field%max_der,para_propa%para_poly%npoly-1)

      CALL alloc_NParray(tab_dnE,[nb_der,3],                          &
                      "tab_dnE",name_sub,[0,1])
      tab_dnE(:,:) = ZERO

      CALL alloc_NParray(work_WP,[para_propa%para_poly%npoly-1],        &
                        'work_WP',name_sub,[0])

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
        work_WP(0) = WP(j)
        max_norm = ZERO
        DO i=0,para_propa%para_poly%npoly-1

          CALL sub_OpPsi(work_WP(i),w2,para_H)
          CALL sub_scaledOpPsi(work_WP(i),w2,para_H%E0,ONE)

          work_WP(i+1) = w2  ! H.psi(i)


          fac_k = cmplx(ZERO,-para_propa%WPdeltaT /                     &
                            real(i+1,kind=Rkind),kind=Rkind)


          work_WP(i+1) = work_WP(i+1) * fac_k

          DO ip=1,3
            IF (.NOT. para_field%pola_xyz(ip)) CYCLE
            fac_k = cmplx(ZERO,-para_propa%WPdeltaT /                   &
                                 real(i+1,kind=Rkind),kind=Rkind)

            w1 = work_WP(i) * (tab_dnE(0,ip)* fac_k)

            DO k=1,min(i,max_der)
              fac_k = fac_k * para_propa%WPdeltaT/real(k,kind=Rkind)
              rt = tab_dnE(k,ip)*fac_k
              w1 = w1 + work_WP(i-k) * rt
            END DO

            CALL sub_OpPsi(w1,w2,para_Dip(ip))

            work_WP(i+1) = work_WP(i+1) + w2
          END DO

          IF (taylor) THEN
            w1 = work_WP(i+1)
          ELSE
!           H contribution
            CALL sub_OpPsi(work_WP(i),w2,para_H)
            CALL sub_scaledOpPsi(work_WP(i),w2,para_H%E0,ONE)
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
              CALL sub_OpPsi(work_WP(i),w2,para_Dip(ip))
              w1 = w1 + w2 * IntE(ip)
            END DO

            w1 = w1 * (-EYE)

          END IF

          WP(j) = WP(j) + w1

          CALL norm2_psi(w1)
          IF (w1%norm2 > max_norm) max_norm = w1%norm2
          write(out_unitp,*) 'norms of w1(i),max_norm',i,w1%norm2,max_norm

          IF (w1%norm2 < para_propa%para_poly%poly_tol) EXIT
          IF (w1%norm2 > TEN**15) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' The norm2 of dnpsi(i+1) is huge !',w1%norm2
            write(out_unitp,*) ' iteration i:',i
            write(out_unitp,*) ' Reduce the time step, WPDeltaT:',              &
                             para_propa%WPdeltaT
            para_propa%march_error = .TRUE.
            STOP
          END IF
        END DO



        IF (print_Op) write(out_unitp,*) 'T,max_der,ipoly,norm2',               &
                                  T,max_der,i,w1%norm2

        work_WP(0) = work_WP(min(i+1,para_propa%para_poly%npoly))

!       - Phase Shift -----------------
        phase = para_H%E0*para_propa%WPdeltaT
        WP(j) = WP(j) * exp(-EYE*phase)

!       - check norm ------------------
        CALL norm2_psi(WP(j),GridRep=.FALSE.,BasisRep=.TRUE.)
        IF ( WP(j)%norm2 > para_propa%max_norm2) THEN
          T  = T + para_propa%WPdeltaT
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' STOP propagation: norm > max_norm',              &
                         WP(j)%norm2
          para_propa%march_error   = .TRUE.
          para_propa%test_max_norm = .TRUE.
          STOP
        END IF

      END DO

      CALL dealloc_NParray(tab_dnE,"tab_dnE",name_sub)
      CALL dealloc_psi(w1)
      CALL dealloc_psi(w2)
      CALL dealloc_NParray(work_WP,'work_WP',name_sub)

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
      USE mod_psi,   ONLY : param_psi,ecri_psi,norm2_psi,               &
                            sub_PsiBasisRep_TO_GridRep,                 &
                            sub_PsiGridRep_TO_BasisRep
      USE mod_field, ONLY : param_field,sub_dnE
      USE mod_Op,    ONLY : param_Op, sub_OpPsi,sub_OpiPsi,sub_scaledOpPsi

      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H
      TYPE (param_Op)   :: para_Dip(3)

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_field) :: para_field

      integer          :: nb_WP
      TYPE (param_psi) :: WP(nb_WP)


!----- for printing --------------------------------------------------
      logical ::print_Op



!------ working parameters --------------------------------
      TYPE (param_psi)       :: w1,w2,w3
      TYPE (param_psi)       :: expT,expV
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
          write(out_unitp,*) 'norm WP',i,WP(i)%norm2

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
        IF ( WP(j)%norm2 > para_propa%max_norm2) THEN
          T  = T + para_propa%WPdeltaT
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' STOP propagation: norm > max_norm',              &
                         WP(j)%norm2
          para_propa%march_error   = .TRUE.
          para_propa%test_max_norm = .TRUE.
          STOP
        END IF

      END DO

      CALL dealloc_psi(w1)
      CALL dealloc_psi(w2)
      CALL dealloc_psi(w3)
      CALL dealloc_psi(expT)
      CALL dealloc_psi(expV)

!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------


      END SUBROUTINE march_split_field

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
         flush(out_unitp)
       END IF
!-----------------------------------------------------------

      IF (para_poly%Hmin > para_poly%Hmax) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) 'Hmax',para_poly%Hmax
         write(out_unitp,*) 'DHmax',para_poly%DHmax
         write(out_unitp,*) 'Hmin',para_poly%Hmin
         write(out_unitp,*) ' Hmin > Hmax '
         !STOP " ERROR in " // name_sub // " : Hmin > Hmax"
         STOP
      END IF

!-----------------------------------------------------------
!     E0 is the zero of energy for energy scaling (i.e. the centre of the range).
!     ESC is the energy scaling parameter ( so that scaled energy lies between -1 and 1)

      para_poly%deltaE = para_poly%Hmax - para_poly%Hmin
      para_poly%E0     = para_poly%Hmin + HALF * para_poly%deltaE
      para_poly%alpha  = HALF * para_poly%deltaE * deltaT
      para_poly%Esc    = ONE

      IF (type_propa == 1)  para_poly%Esc    = HALF * para_poly%deltaE
      IF (type_propa == 33) para_poly%E0     = para_poly%Hmin
      IF (type_propa == 33) para_poly%alpha  = para_poly%deltaE * deltaT


      IF (debug) THEN
        write(out_unitp,*) 'Hmin,Hmax',para_poly%Hmin,para_poly%Hmax
        write(out_unitp,*) 'deltaE',para_poly%deltaE
        write(out_unitp,*) 'E0,Esc',para_poly%E0,para_poly%Esc
        write(out_unitp,*) 'alpha (r)',para_poly%alpha
        flush(out_unitp)
      END IF
!     ------------------------------------------------------

!     ------------------------------------------------------
      IF (type_propa == 1) THEN
!       - Chebychev coefficients calculation ------------------
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
       flush(out_unitp)

       END SUBROUTINE initialisation1_poly


!==============================================================
!
!     Chebychev coefficients calculation
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
         flush(out_unitp)
       END IF
!-----------------------------------------------------------

!----------------------------------------------
      ncheb = max(ONE,r)
      IF (debug) write(out_unitp,*) 'r, ncheb : ',r,ncheb
      flush(out_unitp)

      IF (ncheb > nchmx) THEN
         write(out_unitp,*) ' ERROR in cofcheb'
         write(out_unitp,*) '*** nchmx too small '
         write(out_unitp,*) 'nchmx < ncheb',nchmx,ncheb
         STOP
      END IF

      CALL cof(r,ncheb,cf)

!     -----------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Chebychev coefficients',ncheb
        write(out_unitp,2009) (i,cf(i),i=1,ncheb)
        flush(out_unitp)
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
        write(out_unitp,*) 'Chebychev coefficients',ncheb
        write(out_unitp,2009) (i,cf(i),i=1,ncheb)
 2009   format(3(i6,e16.5))
        write(out_unitp,*) 'END cofcheb'
      END IF
!-----------------------------------------------------------
      flush(out_unitp)

      END SUBROUTINE cofcheb
!==============================================================
!
!     Chebychev coefficients calculation
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

      !write(out_unitp,*) 'Chebychev coefficients'
      !write(out_unitp,*) '  with r=',r1

      DO i=2,ncheb
        cf(i) = cf(i) + cf(i)
        !write(out_unitp,*) 'i,cf',i,cf(i)
      END DO
!stop
      END SUBROUTINE cof

!==============================================================
!
!     Taylor coefficients calculation (nOD)
!
!==============================================================
      !!@description: Taylor coefficients calculation (nOD)
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      SUBROUTINE cof_nOD(alpha,DeltaT,nOD,coef,Max_nOD,epsi)
      USE mod_system
      IMPLICIT NONE

      integer           :: nOD,Max_nOD
      real (kind=Rkind) :: alpha,DeltaT,epsi
      real (kind=Rkind) :: coef(:)



      real (kind=Rkind) :: reste,xi
      integer           :: i

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

      END MODULE mod_march
