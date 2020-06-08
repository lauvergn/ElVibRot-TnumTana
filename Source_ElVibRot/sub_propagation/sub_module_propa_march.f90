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
#if(run_MPI)
#else
        write(out_unitp,*) 'Hmin,Hmax',para_propa%para_poly%Hmin,       &
                                       para_propa%para_poly%Hmax
#endif
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

      CASE (5) ! 4th-Order Runge-Kunta
        CALL march_RK4(T,no,WP(1),WP0(1),para_H,para_propa)

      CASE (6) !
        CALL march_ModMidPoint(T,no,WP(1),WP0(1),para_H,para_propa)

      CASE (7)
        CALL march_BS(T,no,WP(1),WP0(1),para_H,para_propa)

      CASE (8) ! Short Iterative Lanczos
#if(run_MPI)
        IF(SRep_MPI) THEN
          CALL march_SIL_MPI(T,no,WP(1),WP0(1),para_H,para_propa)
        ELSE
#endif
          CALL march_SIL(T,no,WP(1),WP0(1),para_H,para_propa)
#if(run_MPI)
        ENDIF
#endif
      CASE (9) ! Short Iterative Propagation (Lanczos with complete orthogonalisation)
        CALL march_SIP(T,no,WP(1),WP0(1),para_H,para_propa)

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
        CALL sub_analyze_mini_WP_OpWP(T+para_propa%WPdeltaT,WP,1,para_H)

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
      SUBROUTINE march_RK4(T,no,WP,WP0,para_H,para_propa)
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
          write(out_unitp,*) 'norm WP',i,WP%norm2

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
      CALL flush_perso(no)

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

      END SUBROUTINE march_RK4
      SUBROUTINE march_BS(T,no,WP,WP0,para_H,para_propa)
      USE mod_system
      USE mod_psi,   ONLY : param_psi,ecri_psi,norm2_psi,dealloc_psi
      USE mod_Op,    ONLY : param_Op
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
      x = real(2*j+2,kind=Rkind)/real(2*(j-i)+2,kind=Rkind)
      x = ONE/(x**2-ONE)
      yerr = yt1(i-1)-yt0(i-1)
      yt1(i) = yt1(i-1) + yerr*x
    END DO

    yerr = yt1(j) - yt0(j-1)
    CALL norm2_psi(yerr)
    err1 = sqrt(yerr%norm2)

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
    !write(out_unitp,*) 'march_bs',j,err1

    !IF (err1 < 1.d-6 .OR. err1 > err0) EXIT
    IF (err1 < 1.d-10) EXIT

    err0 = err1
  END DO
  write(out_unitp,*) 'end march_bs',min(j,order),err1

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
      CALL flush_perso(no)


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
      CALL alloc_array(tab_dnE,(/nb_der,3/),"tab_dnE",name_sub,(/0,1/))
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

      w2 = psi
      CALL sub_PsiOpPsi(E,psi,w2,para_H)
!     para_H%E0 = real(E,kind=Rkind)
!     write(out_unitp,*) 'E0,Esc',para_H%E0,para_H%Esc

      psi0Hkpsi0(:) = cmplx(ZERO,ZERO,kind=Rkind)

      psi0Hkpsi0(0) = Calc_AutoCorr(psi0,psi,para_propa,T,Write_AC=.FALSE.)

      w1 = psi
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

        w1 = w2

        psi0Hkpsi0(j) = Calc_AutoCorr(psi0,w1,para_propa,T,Write_AC=.FALSE.)

        !CALL Overlap_psi1_psi2(psi0Hkpsi0(j),psi0,w1)

        rtj = rtj *                                                     &
          cmplx(ZERO,-para_propa%WPdeltaT/real(j,kind=Rkind),kind=Rkind)
        w2 = w1 * rtj

        !write(out_unitp,21) 'Rw2*rtj',Real(w2%CvecB,kind=Rkind)
        !write(out_unitp,21) 'Iw2*rtj',AImag(w2%CvecB)

        psi = psi + w2

        !write(out_unitp,21) 'Rpsi',Real(psi%CvecB,kind=Rkind)
        !write(out_unitp,21) 'Ipsi',AImag(psi%CvecB)

        CALL norm2_psi(w2)

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
      psi = psi * exp(-cmplx(ZERO,phase,kind=Rkind))

 !write(out_unitp,*) ' Psi after phase shift '
 !CALL ecri_psi(psi=psi)

!    - check norm ------------------
      CALL norm2_psi(psi,GridRep=.FALSE.,BasisRep=.TRUE.)
      IF ( psi%norm2 > para_propa%max_norm2) THEN
        T  = T + para_propa%WPdeltaT
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' STOP propagation: norm > max_norm',                &
                         psi%norm2
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

      CALL dealloc_psi(w1)
      CALL dealloc_psi(w2)

!-----------------------------------------------------------
      IF (debug) THEN
        CALL norm2_psi(psi)
        write(out_unitp,*) 'norm psi',psi%norm2
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------
!STOP

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
        psi0_psiKrylovSpace(1) = Calc_AutoCorr(psi0,tab_KrylovSpace(1), &
                                         para_propa,T,Write_AC=.FALSE.)
      END IF

      E0 = para_H%E0
      H(:,:) = CZERO
      ! loop for H|psi>, H^2|psi>, H^3|psi>...
      DO j=2,para_propa%para_poly%npoly+1
        IF (debug) write(out_unitp,*) 'in ',name_sub,' it:',j-1
        CALL sub_OpPsi(Psi  =tab_KrylovSpace(j-1),                      &
                       OpPsi=tab_KrylovSpace(j),para_Op=para_H)
        IF (j == 2) THEN
          ! Energy shift, E0, calculation for the first iteration.
          ! since E0=<psi |H| psi> = <tab_KrylovSpace(1) |H|tab_KrylovSpace(1)> =
          ! ..  <tab_KrylovSpace(1) | tab_KrylovSpace(2)>
          ! This shift is important to improve the stapility.
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
        !write(out_unitp,*) j-1,'abs(UPsiOnKrylov(j-1)',abs(UPsiOnKrylov(j-1))
        IF (abs(UPsiOnKrylov(j-1)) < para_propa%para_poly%poly_tol .OR. &
            j == para_propa%para_poly%npoly+1) THEN
          n = j-1
          write(out_unitp,*) n,'abs(UPsiOnKrylov(n)',abs(UPsiOnKrylov(n))
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

        IF (para_propa%nb_micro > 1) THEN
          psi0_psiKrylovSpace(j) = Calc_AutoCorr(psi0,tab_KrylovSpace(j),&
                                          para_propa,T,Write_AC=.FALSE.)
        END IF

      END DO

      IF (abs(UPsiOnKrylov(n)) > para_propa%para_poly%poly_tol) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' The last vector, UPsiOnKrylov(n), coeficient is TOO large'
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
        CALL flush_perso(no)
      ELSE
        cdot = Calc_AutoCorr(psi0,psi,para_propa,T,Write_AC=.FALSE.)
        CALL Write_AutoCorr(no,T + para_propa%WPdeltaT,cdot)
        CALL flush_perso(no)
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
  SUBROUTINE march_SIL(T,no,psi,psi0,para_H,para_propa)
      USE mod_system
      USE mod_Op,    ONLY : param_Op, sub_PsiOpPsi, sub_OpPsi,sub_scaledOpPsi
      USE mod_psi,   ONLY : param_psi,ecri_psi,norm2_psi,renorm_psi,    &
                            renorm_psi_With_norm2,Overlap_psi1_psi2

USE mod_basis_BtoG_GtoB_SGType4,ONLY:tabPackedBasis_TO_tabR_AT_iG

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


      allocate(tab_KrylovSpace(para_propa%para_poly%npoly+1))

      tab_KrylovSpace(1) = psi

      IF (para_propa%nb_micro > 1) THEN
        psi0_psiKrylovSpace(:) = CZERO
!        psi0_psiKrylovSpace(1) = Calc_AutoCorr(psi0,tab_KrylovSpace(1), &
!                                         para_propa,T,Write_AC=.FALSE.)
      END IF

      E0 = para_H%E0
      H(:,:) = CZERO
      DO k=1,para_propa%para_poly%npoly
        IF (debug) write(out_unitp,*) 'in ',name_sub,' it:',k

        IF (para_propa%nb_micro > 1) THEN
          psi0_psiKrylovSpace(k) = Calc_AutoCorr(psi0,tab_KrylovSpace(k),&
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
          CALL Overlap_psi1_psi2(Overlap,tab_KrylovSpace(1),w1)
          E0 = real(Overlap,kind=Rkind)
        END IF
        CALL sub_scaledOpPsi(Psi=tab_KrylovSpace(k),OpPsi=w1,E0=E0,Esc=ONE)

        ! |w_k> = H| tab_KrylovSpace(k)> - beta_k-1 * |tab_KrylovSpace(k-1)>
        ! |w_k> = |w1>                   -MatH(k,k-1)*|tab_KrylovSpace(k-1)>
        IF (k == 1) THEN
          ! w1 = w1
        ELSE
          w1 = w1 - H(k,k-1) * tab_KrylovSpace(k-1)
        END IF

        ! alpha_k = <w_k | tab_KrylovSpace(k)> = <w1 | tab_KrylovSpace(k)>
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
          EXIT
        END IF

        ! orthogonalisation: |v_k> = |w_k> - alpha_k * |tab_KrylovSpace(k)>
        !                     |w1> = |w1>  - Overlap * |tab_KrylovSpace(k)>
        w1 = w1 - Overlap*tab_KrylovSpace(k)

        !beta_k = <v_k | v_k> = <w1 | w1>
        !       => beta_k = sqrt(norm2)
        CALL norm2_psi(w1)
        H(k+1,k) = sqrt(w1%norm2)
        H(k,k+1) = H(k+1,k)
        CALL renorm_psi_With_norm2(w1)
        tab_KrylovSpace(k+1) = w1

      END DO

      IF (debug) write(out_unitp,*) 'abs(UPsiOnKrylov)',abs(UPsiOnKrylov(1:n))

      Psi = ZERO
      DO k=1,n
        Psi = Psi + UPsiOnKrylov(k)*tab_KrylovSpace(k)
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
        CALL flush_perso(no)
      ELSE
        cdot = Calc_AutoCorr(psi0,psi,para_propa,T,Write_AC=.FALSE.)
        CALL Write_AutoCorr(no,T + para_propa%WPdeltaT,cdot)
        CALL flush_perso(no)
      END IF


      ! deallocation
      DO k=1,size(tab_KrylovSpace)
         CALL dealloc_psi(tab_KrylovSpace(k))
      END DO
      deallocate(tab_KrylovSpace)
      IF (allocated(Vec))          CALL dealloc_NParray(Vec,'Vec',name_sub)
      IF (allocated(Eig))          CALL dealloc_NParray(Eig,'Eig',name_sub)
      IF (allocated(UPsiOnKrylov)) CALL dealloc_NParray(UPsiOnKrylov,'UPsiOnKrylov',name_sub)
      CALL dealloc_psi(w1)
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

        CALL diagonalization_HerCplx(H,Eig,Vec,n,3,1,.TRUE.)
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
        write(out_unitp,*) 'abs(UPsiOnKrylov)',abs(UPsiOnKrylov)
        write(out_unitp,*) 'norm2 UPsiOnKrylov',dot_product(UPsiOnKrylov,UPsiOnKrylov)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------


  END SUBROUTINE UPsi_spec
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

      CALL diagonalization_HerCplx(H,Eig,Vec,n,3,1,.TRUE.)

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
!      CALL flush_perso(no)

      cdot = Calc_AutoCorr(psi0,psi,para_propa,T,Write_AC=.FALSE.)
      CALL Write_AutoCorr(no,T + para_propa%WPdeltaT,cdot)
      CALL flush_perso(no)


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
          E0 = real(Overlap)
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

      CALL diagonalization_HerCplx(H,Eig,Vec,para_propa%para_poly%npoly,&
                                   3,1,.TRUE.)

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
      CALL flush_perso(no)


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
      TYPE (param_psi), pointer :: w1,w2
      complex (kind=Rkind) :: cdot
      !real (kind=Rkind)    :: microT,microdeltaT,phase,microphase



      complex (kind=Rkind), allocatable :: PsiOnSpectral(:)
      complex (kind=Rkind), allocatable :: HPsiOnSpectral(:)

      TYPE (param_psi), allocatable :: tab_KrylovSpace(:)
      complex (kind=Rkind) :: Overlap,ReNor

      logical, parameter :: test = .FALSE.
      !logical, parameter :: test = .TRUE.

      integer              :: it,i,j

!----- for debuging --------------------------------------------------
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
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

      CALL alloc_NParray(PsiOnSpectral, [para_H%nb_tot],'PsiOnSpectral', name_sub)
      CALL alloc_NParray(HPsiOnSpectral,[para_H%nb_tot],'HPsiOnSpectral',name_sub)

!      ! 1st: project psi on the spectral basis
!      IF (para_H%cplx) THEN
!        PsiOnSpectral(:) = matmul(transpose(para_H%Cvp),Psi%CvecB)
!      ELSE
!        PsiOnSpectral(:) = matmul(transpose(para_H%Rvp),Psi%CvecB)
!      END IF
!
!      ! 2d: propagation with spectral representation
!      IF (para_H%cplx) THEN
!        HPsiOnSpectral   = exp(-EYE*para_H%Cdiag*para_propa%WPdeltaT)*PsiOnSpectral
!      ELSE
!        HPsiOnSpectral   = exp(-EYE*para_H%Rdiag*para_propa%WPdeltaT)*PsiOnSpectral
!      END IF
!
!      ! 3d: back to the initial basis
!      IF (para_H%cplx) THEN
!        Psi%CvecB = matmul(conjg(para_H%Cvp),HPsiOnSpectral)
!      ELSE
!        Psi%CvecB = matmul(para_H%Rvp,HPsiOnSpectral)
!      END IF


      ! 1st: project psi on the spectral basis
      IF (para_H%cplx) THEN
        PsiOnSpectral(:) = matmul(transpose(para_H%Cvp),Psi%CvecB)
      ELSE
        DO i=1,para_H%nb_tot
          PsiOnSpectral(i) = dot_product(para_H%Rvp(:,i),Psi%CvecB)
        END DO
      END IF

      ! 2d: propagation with spectral representation
      IF (para_H%cplx) THEN
        HPsiOnSpectral   = exp(-EYE*para_H%Cdiag*para_propa%WPdeltaT)*PsiOnSpectral
      ELSE
        HPsiOnSpectral   = exp(-EYE*para_H%Rdiag*para_propa%WPdeltaT)*PsiOnSpectral
      END IF

      ! 3d: back to the initial basis
      IF (para_H%cplx) THEN
        Psi%CvecB = matmul(conjg(para_H%Cvp),HPsiOnSpectral)
      ELSE
        Psi%CvecB = CZERO
        DO i=1,para_H%nb_tot
          Psi%CvecB(:) = Psi%CvecB(:) + HPsiOnSpectral(i) * para_H%Rvp(:,i)
        END DO
      END IF


      CALL dealloc_NParray(PsiOnSpectral, 'PsiOnSpectral', name_sub)
      CALL dealloc_NParray(HPsiOnSpectral,'HPsiOnSpectral',name_sub)


      cdot = Calc_AutoCorr(psi0,psi,para_propa,T,Write_AC=.FALSE.)
      CALL Write_AutoCorr(no,T + para_propa%WPdeltaT,cdot)
      CALL flush_perso(no)

!-----------------------------------------------------------
      IF (debug) THEN
        CALL norm2_psi(psi)
        write(out_unitp,*) 'norm psi',psi%norm2
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------
!STOP

 END SUBROUTINE march_Spectral

#if(run_MPI)
!=======================================================================================
!> @brief MPI version of march_SIL
!> according to Smolyak rep. 
!> Note, this works for massive cluster with big memmory
!
! require a new action subroutine, SR vector ready in V
! rewrite the atcion function to make it easier for callimng from different subroutines
!=======================================================================================  
  SUBROUTINE march_SIL_MPI(TT,no,psi,psi0,para_H,para_propa)
    USE mod_system
    USE mod_Op,                     ONLY:param_Op,sub_OpPsi
    USE mod_OpPsi,                  ONLY:sub_scaledOpPsi_SR_MPI
    USE mod_psi_set_alloc,          ONLY:param_psi,psi_times_SR_MPI
    USE mod_psi_Op,                 ONLY:norm2_psi_SR_MPI,Overlap_psi1_psi2_SRB_MPI
    USE mod_psi,                    ONLY:norm2_psi
    USE mod_propa,                  ONLY:Calc_AutoCorr_SR_MPI 
    USE mod_OpPsi_SG4,              ONLY:ini_iGs_MPI
    USE mod_basis_BtoG_GtoB_SGType4,ONLY:TypeRVec
    IMPLICIT NONE

    !-variables for namelist minimum----------------------------------------------------
    TYPE(param_Op),                 intent(in)    :: para_H

    !-variables for the WP propagation--------------------------------------------------
    TYPE(param_propa),              intent(inout) :: para_propa
    TYPE(param_psi),                intent(inout) :: psi0
    TYPE(param_psi),                intent(inout) :: psi
    Real(kind=Rkind),               intent(in)    :: TT
    Integer,                        intent(in)    :: no
    
    TYPE(param_psi),allocatable                   :: tab_KrylovSpace(:)
    TYPE(param_psi)                               :: w1
    Complex(kind=Rkind)      :: H(para_propa%para_poly%npoly,para_propa%para_poly%npoly)
    Complex(kind=Rkind)               :: psi0_psiKrylovSpace(para_propa%para_poly%npoly)
    Complex(kind=Rkind),allocatable               :: UPsiOnKrylov(:)
    Complex(kind=Rkind),allocatable               :: Vec(:,:)
    Complex(kind=Rkind)                           :: Overlap
    Complex(kind=Rkind)                           :: cdot
    Real(kind=Rkind),allocatable                  :: Eig(:)
    Real(kind=Rkind)                              :: E0
    Real(kind=Rkind)                              :: phase
    Real(kind=Rkind)                              :: micro_deltaT
    Real(kind=Rkind)                              :: micro_T
    Real(kind=Rkind)                              :: micro_phase
    Integer                                       :: n 
    Integer                                       :: k

    Character(len=*),parameter                    :: name_sub='march_SIL_MPI'


    !> for Krylov Space
    allocate(tab_KrylovSpace(para_propa%para_poly%npoly+1)) ! param_psi 
    
    !> initialize iGs_MPI
    IF(para_H%BasisnD%para_SGType2%once_action) CALL ini_iGs_MPI(para_H,.TRUE.)
    
    !> Extract compact psi to each threads and transfer to gird rep.
    CALL SmolyakR_distribute_SRB_MPI(psi ,para_Op=para_H)
    IF(para_propa%nb_micro>1) CALL SmolyakR_distribute_SRB_MPI(psi0,para_Op=para_H)
    
    !> set first Krylov Space
    tab_KrylovSpace(1)=psi 
    IF(para_propa%nb_micro>1) psi0_psiKrylovSpace(:)=CZERO

    E0=para_H%E0
    H(:,:)=CZERO ! initialize H
    DO k=1,para_propa%para_poly%npoly
      IF(para_propa%nb_micro>1) THEN
        ! <psi(t)|v_k'> ~ <v_0|v_k'>
        psi0_psiKrylovSpace(k)=Calc_AutoCorr_SR_MPI(psi0,tab_KrylovSpace(k),           &
                                                    para_propa,TT,Write_AC=.FALSE.)
      END IF

      ! |w1>=H|v_k>
      CALL sub_OpPsi(Psi=tab_KrylovSpace(k),OpPsi=w1,para_Op=para_H)
      
      ! Warning:require this two line to get same result as normal case. 
      ! waiting for new idea
      !CALL SmolyakR_to_packedB_SRB_MPI(w1,para_Op=para_H)
      !CALL SmolyakR_distribute_SRB_MPI(w1 ,para_Op=para_H)

      ! case k=1, alpha_1=<v_1|H|v_1>
      IF(k==1) THEN 
        !> Energy shift, E0, calculation for the first iteration. E0=<psi |H| psi> 
        !> This shift is important to improve the stapility
        !> Note phase need to be taking into account at the end of the iterations
        CALL Overlap_psi1_psi2_SRB_MPI(Overlap,tab_KrylovSpace(1),w1)
        !CALL MPI_Reduce_sum_Bcast(Overlap) ! reduce sum and boardcast
        E0=real(Overlap,kind=Rkind)
      ENDIF
      ! scaling with E0 
      CALL sub_scaledOpPsi_SR_MPI(Psi=tab_KrylovSpace(k),OpPsi=w1,E0=E0,Esc=ONE)

      ! v_{k}=w_{k-1}/||w_{k-1}||=w_{k-1}/beta_{k-1}
      ! k=1: w_{k}=H|v_{k}>-alpha_k|v_{k}>
      ! k>1: w_{k}=H|v_{k}>-alpha_k|v_{k}>-beta_{k-1}|v_{k-1}>
      ! beta_{k}=<v_k|H|v_{k+1}>=<v_{k+1}|H|v_{k}>=H(k,k+1)=H(k+1,k)
      IF(k>1) THEN
        !w1=w1-H(k,k-1)*tab_KrylovSpace(k-1) !> H|v_{k}>-beta_{k-1}|v_{k-1}>
        CALL w1_minus_Cv_SR_MPI(w1,H(k,k-1),tab_KrylovSpace(k-1))
      ENDIF
      ! alpha_k=<v_k|H|v_k>
      CALL Overlap_psi1_psi2_SRB_MPI(Overlap,w1,tab_KrylovSpace(k))
      !CALL MPI_Reduce_sum_Bcast(Overlap)
      H(k,k)=Overlap ! alpha_k

      !> diagonalization then exit when:
      !> 1. k reaches npoly 2. the scheme converges
      ! psi(t+dt)=U|psi(t)>=sum_{l=1}^n <vec_l|psi(t)> e^{-i*Ei*dt}|vec_l>
      !                    =sum_{k=1}^n a_k|v_k>
      !                 a_k=sum_{l=1}^n<vec_l|v0>e^{-i*Ei*dt}<v_s|vec_l>
      ! works on all the threads
      CALL UPsi_spec(UPsiOnKrylov,H(1:k,1:k),Vec,Eig,                                  &
                                                para_propa%WPdeltaT,k,With_diago=.TRUE.)

      IF(abs(UPsiOnKrylov(k))<para_propa%para_poly%poly_tol .OR.                       &
         k==para_propa%para_poly%npoly) THEN
        n=k
        write(out_unitp,*) n,'abs(UPsiOnKrylov(n)',abs(UPsiOnKrylov(n))
        EXIT
      END IF

      !> w_k=(H|v_{k}>-beta_{k-1}|v_{k-1}>)-alpha_k|v_{k}>="w1"-alpha_k|v_k>
      ! w1=w1-Overlap*tab_KrylovSpace(k) !< now "w1" is w_k
      CALL w1_minus_Cv_SR_MPI(w1,Overlap,tab_KrylovSpace(k))

      !> beta_{k}=||w_k||=sqrt(<w_k|w_k>)
      CALL norm2_psi_SR_MPI(w1,2) !< 2: just calculate the normalization constant 
      !CALL MPI_Reduce_sum_Bcast(w1%norm2)
      H(k+1,k)=sqrt(w1%norm2)
      H(k,k+1)=conjg(H(k+1,k))

      !> assign v_{k+1}=w_k/beta_k
      CALL norm2_psi_SR_MPI(w1,3) !< 3: normalize SR_B with existing norm. constant
      tab_KrylovSpace(k+1)=w1    
    ENDDO ! for k=1,para_propa%para_poly%npoly
    !CALL SRB_to_packB_write(w1,call_from='endofSILloop')

    psi=ZERO
    DO k=1,n
      ! psi(t+dt)=sum_{k=1}^n a_k|v_k>
      ! psi=psi+UPsiOnKrylov(k)*tab_KrylovSpace(k)
      CALL w1_minus_Cv_SR_MPI(psi,-UPsiOnKrylov(k),tab_KrylovSpace(k))
    END DO

    !> check normalization
    CALL norm2_psi_SR_MPI(psi,2) !< 2: just calculate the normalization constant 
    !CALL norm2_psi(psi)

    IF(psi%norm2 > para_propa%max_norm2) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' STOP propagation: norm > max_norm',psi%norm2
      para_propa%march_error  =.TRUE.
      para_propa%test_max_norm=.TRUE.
      STOP
    ENDIF

    ! Phase Shift
    phase=E0*para_propa%WPdeltaT  
    !psi=psi*exp(-EYE*phase)
    CALL psi_times_SR_MPI(psi,exp(-EYE*phase))
    
    ! transfer back to packed basis, valid on master only
    CALL SmolyakR_to_packedB_SRB_MPI(psi,para_Op=para_H)

    ! autocorelation
    IF(para_propa%nb_micro>1) THEN
      micro_deltaT=para_propa%WPdeltaT/real(para_propa%nb_micro,kind=Rkind)
      micro_phase =phase/real(para_propa%nb_micro,kind=Rkind)
    
      phase =ZERO
      micro_T=ZERO
      
      DO k=1,para_propa%nb_micro
        micro_T=micro_T+micro_deltaT
        phase=phase+micro_phase
        
        CALL UPsi_spec(UPsiOnKrylov,H,Vec,Eig,micro_T,n,With_diago=.FALSE.)
        ! sum(a_k <v0|v_k'>)
        cdot=sum(UPsiOnKrylov(1:n)*psi0_psiKrylovSpace(1:n)) ! we cannot use dot_product
        cdot=cdot*exp(-EYE*phase)
        IF(MPI_id==0) CALL Write_AutoCorr(no,TT+micro_T,cdot)
      ENDDO 
      CALL flush_perso(no)
    ELSE
      !cdot=Calc_AutoCorr_SR_MPI(psi0,psi,para_propa,TT,Write_AC=.FALSE.)
      IF(MPI_id==0) THEN
        cdot=Calc_AutoCorr(psi0,psi,para_propa,TT,Write_AC=.FALSE.)
        CALL Write_AutoCorr(no,TT+para_propa%WPdeltaT,cdot)
        CALL flush_perso(no)
      ENDIF
    ENDIF

    ! transfer back to packed basis, valid on master only
    !CALL SmolyakR_to_packedB_SRB_MPI(psi,para_Op=para_H)
    
    ! deallocation  
    DO k=1,size(tab_KrylovSpace)
       CALL dealloc_psi(tab_KrylovSpace(k)) ! deallocation of SR_B included
    ENDDO
    deallocate(tab_KrylovSpace)

    IF(allocated(Vec))          CALL dealloc_NParray(Vec,'Vec',name_sub)
    IF(allocated(Eig))          CALL dealloc_NParray(Eig,'Eig',name_sub)
    IF(allocated(UPsiOnKrylov)) CALL dealloc_NParray(UPsiOnKrylov,'UPsiOnKrylov',name_sub)
    CALL dealloc_psi(w1)

  END SUBROUTINE march_SIL_MPI
!=======================================================================================
#endif

#if(run_MPI)
!=======================================================================================
!> transfer from Smolyak basis to packed basis
!=======================================================================================
  SUBROUTINE SRB_to_packB_write(psi,call_from)
    USE mod_system
    USE mod_Op,                     ONLY:param_Op
    USE mod_psi_set_alloc,          ONLY:param_psi
    USE mod_basis_BtoG_GtoB_SGType4,ONLY:tabR_AT_iG_TO_tabPackedBasis
    USE mod_psi_Op,                 ONLY:Overlap_psi1_psi2_SRB_MPI
    USE mod_MPI_Aid
    IMPLICIT NONE

    TYPE(param_psi),                intent(inout) :: psi
    Character(len=*),               intent(in)    :: call_from
    
    Complex(kind=Rkind)                           :: Overlap
    Real(kind=Rkind),allocatable                  :: SRB(:)
    Real(kind=Rkind),allocatable                  :: CvecR(:)
    Real(kind=Rkind),allocatable                  :: CvecC(:)
    Integer                                       :: iG
    Integer                                       :: d1
    Integer                                       :: d2

    CALL allocate_array(CvecR,1,size(psi%CvecB))
    CALL allocate_array(CvecC,1,size(psi%CvecB))
    DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)
      d1=psi%SR_B_index(iG)
      d2=psi%SR_B_index(iG+1)-1
      CALL allocate_array(SRB,1,d2-d1+1)
      SRB=psi%SR_B(d1:d2,1) 
      CALL tabR_AT_iG_TO_tabPackedBasis(CvecR,SRB,iG,psi%BasisnD%para_SGType2,         &
                                        psi%BasisnD%WeightSG(iG)) 
      SRB=psi%SR_B(d1:d2,2)
      CALL tabR_AT_iG_TO_tabPackedBasis(CvecC,SRB,iG,psi%BasisnD%para_SGType2,         &
                                        psi%BasisnD%WeightSG(iG)) 
    ENDDO
    
    CALL Overlap_psi1_psi2_SRB_MPI(Overlap,psi,psi)
    
    write(*,*) call_from, ' Real part: ',CvecR
    write(*,*) call_from, ' imag part: ',CvecC
    write(*,*) call_from, ' overlap  : ',Overlap
  END SUBROUTINE SRB_to_packB_write
!=======================================================================================
#endif

#if(run_MPI)
!=======================================================================================
!> @brief MPI version of march_SIL
!> according to Smolyak rep. 
!> Note, this works for massive cluster with big memmory
!
! require a new action subroutine, SR vector ready in V
! rewrite the atcion function to make it easier for callimng from different subroutines
!=======================================================================================  
  SUBROUTINE march_SIL_MPI_old(TT,no,psi,psi0,para_H,para_propa)
    USE mod_system
    USE mod_Op,                     ONLY:param_Op,sub_OpPsi
    USE mod_OpPsi,                  ONLY:sub_scaledOpPsi_SR_MPI
    USE mod_psi_set_alloc,          ONLY:param_psi,psi_times_SR_MPI
    USE mod_psi_Op,                 ONLY:norm2_psi_SR_MPI,Overlap_psi1_psi2_SRG_MPI
    USE mod_propa,                  ONLY:Calc_AutoCorr_SR_MPI 
    USE mod_OpPsi_SG4,              ONLY:ini_iGs_MPI
    USE mod_basis_BtoG_GtoB_SGType4,ONLY:TypeRVec
    IMPLICIT NONE

    !-variables for namelist minimum----------------------------------------------------
    TYPE(param_Op),                 intent(in)    :: para_H

    !-variables for the WP propagation--------------------------------------------------
    TYPE(param_propa),              intent(inout) :: para_propa
    TYPE(param_psi),                intent(inout) :: psi0
    TYPE(param_psi),                intent(inout) :: psi
    Real(kind=Rkind),               intent(in)    :: TT
    Integer,                        intent(in)    :: no
    
    TYPE(param_psi),allocatable                   :: tab_KrylovSpace(:)
    TYPE(param_psi)                               :: w1
    Complex(kind=Rkind)      :: H(para_propa%para_poly%npoly,para_propa%para_poly%npoly)
    Complex(kind=Rkind)               :: psi0_psiKrylovSpace(para_propa%para_poly%npoly)
    Complex(kind=Rkind),allocatable               :: UPsiOnKrylov(:)
    Complex(kind=Rkind),allocatable               :: Vec(:,:)
    Complex(kind=Rkind)                           :: Overlap
    Complex(kind=Rkind)                           :: cdot
    Real(kind=Rkind),allocatable                  :: Eig(:)
    Real(kind=Rkind)                              :: E0
    Real(kind=Rkind)                              :: phase
    Real(kind=Rkind)                              :: micro_deltaT
    Real(kind=Rkind)                              :: micro_T
    Real(kind=Rkind)                              :: micro_phase
    Integer                                       :: n 
    Integer                                       :: k

    Character(len=*),parameter                    :: name_sub='march_SIL_MPI'

    !> for Krylov Space
    allocate(tab_KrylovSpace(para_propa%para_poly%npoly+1)) ! param_psi 
    
    !> initialize iGs_MPI
    IF(para_H%BasisnD%para_SGType2%once_action) CALL ini_iGs_MPI(para_H,.TRUE.)
    
    !> Extract compact psi to each threads and transfer to gird rep.
    CALL SmolyakR_distribute_SRG_MPI(psi ,para_Op=para_H)
    CALL SmolyakR_distribute_SRG_MPI(psi0,para_Op=para_H)
    
    !> set first Krylov Space
    tab_KrylovSpace(1)=psi 
    IF(para_propa%nb_micro>1) psi0_psiKrylovSpace(:)=CZERO
    
    E0=para_H%E0
    H(:,:)=CZERO ! initialize H
    DO k=1,para_propa%para_poly%npoly
    
      IF(para_propa%nb_micro>1) THEN
        ! <psi(t)|v_k'> ~ <v_0|v_k'>
        psi0_psiKrylovSpace(k)=Calc_AutoCorr_SR_MPI(psi0,tab_KrylovSpace(k),           &
                                                    para_propa,TT,Write_AC=.FALSE.)
      END IF

      ! |w1>=H|v_k>
      CALL sub_OpPsi(Psi=tab_KrylovSpace(k),OpPsi=w1,para_Op=para_H)

      ! case k=1, alpha_1=<v_1|H|v_1>
      IF(k==1) THEN 
        !> Energy shift, E0, calculation for the first iteration. E0=<psi |H| psi> 
        !> This shift is important to improve the stapility
        !> Note phase need to be taking into account at the end of the iterations
        CALL Overlap_psi1_psi2_SRG_MPI(Overlap,tab_KrylovSpace(1),w1)
        !CALL MPI_Reduce_sum_Bcast(Overlap) ! reduce sum and boardcast
        E0=real(Overlap,kind=Rkind)
      ENDIF

      ! scaling with E0 
      CALL sub_scaledOpPsi_SR_MPI(Psi=tab_KrylovSpace(k),OpPsi=w1,E0=E0,Esc=ONE)
      
      ! v_{k}=w_{k-1}/||w_{k-1}||=w_{k-1}/beta_{k-1}
      ! k=1: w_{k}=H|v_{k}>-alpha_k|v_{k}>
      ! k>1: w_{k}=H|v_{k}>-alpha_k|v_{k}>-beta_{k-1}|v_{k-1}>
      ! beta_{k}=<v_k|H|v_{k+1}>=<v_{k+1}|H|v_{k}>=H(k,k+1)=H(k+1,k)
      IF(k>1) THEN
        !w1=w1-H(k,k-1)*tab_KrylovSpace(k-1) !> H|v_{k}>-beta_{k-1}|v_{k-1}>
        CALL w1_minus_Cv_SR_MPI(w1,H(k,k-1),tab_KrylovSpace(k-1))
      ENDIF

      ! alpha_k=<v_k|H|v_k>
      CALL Overlap_psi1_psi2_SRG_MPI(Overlap,w1,tab_KrylovSpace(k))
      !CALL MPI_Reduce_sum_Bcast(Overlap)
      H(k,k)=Overlap ! alpha_k

      !> diagonalization then exit when:
      !> 1. k reaches npoly 2. the scheme converges
      ! psi(t+dt)=U|psi(t)>=sum_{l=1}^n <vec_l|psi(t)> e^{-i*Ei*dt}|vec_l>
      !                    =sum_{k=1}^n a_k|v_k>
      !                 a_k=sum_{l=1}^n<vec_l|v0>e^{-i*Ei*dt}<v_s|vec_l>
      ! works on all the threads
      CALL UPsi_spec(UPsiOnKrylov,H(1:k,1:k),Vec,Eig,                                  &
                                                para_propa%WPdeltaT,k,With_diago=.TRUE.)
      IF(abs(UPsiOnKrylov(k))<para_propa%para_poly%poly_tol .OR.                       &
         k==para_propa%para_poly%npoly) THEN
        n=k
        write(out_unitp,*) n,'abs(UPsiOnKrylov(n)',abs(UPsiOnKrylov(n))
        EXIT
      END IF

      !> w_k=(H|v_{k}>-beta_{k-1}|v_{k-1}>)-alpha_k|v_{k}>="w1"-alpha_k|v_k>
      ! w1=w1-Overlap*tab_KrylovSpace(k) !< now "w1" is w_k
      CALL w1_minus_Cv_SR_MPI(w1,Overlap,tab_KrylovSpace(k))

      !> beta_{k}=||w_k||=sqrt(<w_k|w_k>)
      CALL norm2_psi_SR_MPI(w1,2) !< 2: just calculate the normalization constant 
      !CALL MPI_Reduce_sum_Bcast(w1%norm2)
      H(k+1,k)=sqrt(w1%norm2)
      H(k,k+1)=H(k+1,k)

      !> assign v_{k+1}=w_k/beta_k
      CALL norm2_psi_SR_MPI(w1,3) !< 3: normalize SR_G with existing norm. constant
      tab_KrylovSpace(k+1)=w1    
    ENDDO ! for k=1,para_propa%para_poly%npoly
    write(out_unitp,*) 'abs(UPsiOnKrylov)',abs(UPsiOnKrylov(1:n)), 'from ',MPI_id

    psi=ZERO
    DO k=1,n
      ! psi(t+dt)=sum_{k=1}^n a_k|v_k>
      ! psi=psi+UPsiOnKrylov(k)*tab_KrylovSpace(k)
      CALL psi_times_SR_MPI(tab_KrylovSpace(k),UPsiOnKrylov(k),psi)
    END DO

    !> check normalization
    CALL norm2_psi_SR_MPI(psi,2) !< 2: just calculate the normalization constant 
    IF(psi%norm2 > para_propa%max_norm2) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' STOP propagation: norm > max_norm',psi%norm2
      para_propa%march_error  =.TRUE.
      para_propa%test_max_norm=.TRUE.
      STOP
    ENDIF

    ! Phase Shift
    phase=E0*para_propa%WPdeltaT  
    ! psi=psi*exp(-EYE*phase)
    CALL psi_times_SR_MPI(psi,exp(-EYE*phase))

    ! autocorelation
    IF(para_propa%nb_micro>1) THEN
      micro_deltaT=para_propa%WPdeltaT/real(para_propa%nb_micro,kind=Rkind)
      micro_phase =phase/real(para_propa%nb_micro,kind=Rkind)
    
      phase =ZERO
      micro_T=ZERO
      
      DO k=1,para_propa%nb_micro
        micro_T=micro_T+micro_deltaT
        phase=phase+micro_phase
        
        CALL UPsi_spec(UPsiOnKrylov,H,Vec,Eig,micro_T,n,With_diago=.FALSE.)
        ! sum(a_k <v0|v_k'>)
        cdot=sum(UPsiOnKrylov(1:n)*psi0_psiKrylovSpace(1:n)) ! we cannot use dot_product
        cdot=cdot*exp(-EYE*phase)
        CALL Write_AutoCorr(no,TT+micro_T,cdot)
      ENDDO 
      CALL flush_perso(no)
    ELSE
      cdot=Calc_AutoCorr_SR_MPI(psi0,psi,para_propa,TT,Write_AC=.FALSE.)
      CALL Write_AutoCorr(no,TT+para_propa%WPdeltaT,cdot)
      CALL flush_perso(no)
    ENDIF

    ! deallocation  
    DO k=1,size(tab_KrylovSpace)
       CALL dealloc_psi(tab_KrylovSpace(k)) ! deallocation of SR_G included
    ENDDO
    deallocate(tab_KrylovSpace)

    IF(allocated(Vec))          CALL dealloc_NParray(Vec,'Vec',name_sub)
    IF(allocated(Eig))          CALL dealloc_NParray(Eig,'Eig',name_sub)
    IF(allocated(UPsiOnKrylov)) CALL dealloc_NParray(UPsiOnKrylov,'UPsiOnKrylov',name_sub)
    CALL dealloc_psi(w1)

    ! update 
    CALL SmolyakR_to_packedB_SRG_MPI(psi,para_Op=para_H)

  END SUBROUTINE march_SIL_MPI_old
!=======================================================================================
#endif

#if(run_MPI)
!=======================================================================================
!> @brief distribute Smolyak terms to threads
!> Only the real part of the packed basis are converged to Smolyak rep. 
!> 1. transfer from packed basis to basis in Smolyak rep.
!> 2. transfer from basis Re. to grid rep. 
!> 3. grid Smolyak rep. are reserved for the current time step
!
!> if psi%clpx, the final SR_B will be (:,2) dim matrix, with real and imag. part resp.
!> otherwise, SR_B will be (:,1)
!=======================================================================================
  SUBROUTINE SmolyakR_distribute_SRB_MPI(psi,para_Op)
    USE mod_system
    USE mod_basis_set_alloc
    USE mod_param_SGType2
    USE mod_nDindex
    USE mod_psi_set_alloc,          ONLY:param_psi
    USE mod_Op,                     ONLY:param_Op
    USE mod_basis_BtoG_GtoB_SGType4,ONLY:tabPackedBasis_TO_tabR_AT_iG,                 &
                                         BDP_TO_GDP_OF_SmolyakRep,                     &
                                         getbis_tab_nq,getbis_tab_nb
    USE mod_MPI_Aid
    IMPLICIT NONE
    
    TYPE(param_psi),                intent(inout) :: psi
    TYPE(param_Op),                 intent(in)    :: para_Op
    
    TYPE(param_SGType2),pointer                   :: SGType2
    Real(kind=Rkind),allocatable                  :: SR_B0(:,:)
    Integer,pointer                               :: tab_l(:,:)
    Integer,allocatable                           :: tab_nq(:)
    Integer,allocatable                           :: tab_nb(:)
    Integer                                       :: iG
    Integer                                       :: nb0
    Integer                                       :: dim
    
    Character(len=*),parameter             :: name_sub='SmolyakR_distribute_SR_MPI'
    
    psi%SRB_MPI=.TRUE.

    SGType2 => para_Op%BasisnD%para_SGType2
    tab_l   => SGType2%nDind_SmolyakRep%Tab_nDval(:,:)
    nb0=SGType2%nb0

    !> allocate grid SR, the relevant index record the vector position for each iG
    !> psi%SR_B(psi%SR_B_index(iG):psi%SR_B_index(iG+1)-1) => SRB(iG)
    CALL allocate_array(psi%SR_B_length,0,MPI_np-1)
    CALL allocate_array(psi%SR_B_index,iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)+1)
    psi%SR_B_index(iGs_MPI(1,MPI_id))=1

    DO i_MPI=0,MPI_np-1
      psi%SR_B_length(i_MPI)=0
      DO iG=iGs_MPI(1,i_MPI),iGs_MPI(2,i_MPI)
        psi%SR_B_length(i_MPI)=psi%SR_B_length(i_MPI)+SGType2%tab_nb_OF_SRep(iG)*nb0
        IF(i_MPI==MPI_id) psi%SR_B_index(iG+1)=psi%SR_B_length(MPI_id)+1 ! for each threads
      ENDDO
    ENDDO

!    DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)+1
!      write(*,*) 'iGs_MPI check: iG=',iG,psi%SR_B_index(iG),iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)+1
!    ENDDO

    !> allocate SR_B for each threads
    IF(psi%cplx) THEN
      CALL allocate_array(psi%SR_B,1,psi%SR_B_length(MPI_id),1,2)
    ELSE 
      CALL allocate_array(psi%SR_B,1,psi%SR_B_length(MPI_id),1,1)
    ENDIF

    !> distribute SR on Basis to each threads
    IF(MPI_id==0) THEN
      DO i_MPI=1,MPI_np-1
      
        !> allocate temprary array for sending
        IF(psi%cplx) THEN
          CALL allocate_array(SR_B0,1,psi%SR_B_length(i_MPI),1,2)
        ELSE
          CALL allocate_array(SR_B0,1,psi%SR_B_length(i_MPI),1,1)
        ENDIF
        
        !> loop to send SR_B
        DO iG=iGs_MPI(1,i_MPI),iGs_MPI(2,i_MPI)
          IF(psi%cplx) THEN
            ! Real part
            CALL pack_Basis_to_SR_Basis_MPI(iG,Real(psi%CvecB,kind=Rkind),SR_B0(:,1),  &
                                            psi%SR_B_index,para_Op)
            ! Imaginary part
            CALL pack_Basis_to_SR_Basis_MPI(iG,aimag(psi%CvecB),SR_B0(:,2),            &
                                            psi%SR_B_index,para_Op)
          ELSE
            CALL pack_Basis_to_SR_Basis_MPI(iG,psi%RvecB,SR_B0(:,1),                   &
                                            psi%SR_B_index,para_Op)
          ENDIF
        ENDDO
      
        ! send SR_B
        IF(psi%cplx) THEN
          CALL MPI_Send(SR_B0(1:psi%SR_B_length(i_MPI),1),psi%SR_B_length(i_MPI),      &
                              MPI_real_fortran,root_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
          CALL MPI_Send(SR_B0(1:psi%SR_B_length(i_MPI),2),psi%SR_B_length(i_MPI),      &
                              MPI_real_fortran,root_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
        ELSE
          CALL MPI_Send(SR_B0(1:psi%SR_B_length(i_MPI),1),psi%SR_B_length(i_MPI),      &
                              MPI_real_fortran,root_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
        ENDIF
      ENDDO
      
      ! Smolyak rep. on MPI_id=0
      DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)
        IF(psi%cplx) THEN
          ! Real part
          CALL pack_Basis_to_SR_Basis_MPI(iG,Real(psi%CvecB,kind=Rkind),psi%SR_B(:,1), &
                                          psi%SR_B_index,para_Op)
          ! Imaginary part
          CALL pack_Basis_to_SR_Basis_MPI(iG,aimag(psi%CvecB),psi%SR_B(:,2),           &
                                          psi%SR_B_index,para_Op)
        ELSE
          CALL pack_Basis_to_SR_Basis_MPI(iG,psi%RvecB,psi%SR_B(:,1),                  &
                                          psi%SR_B_index,para_Op)
        ENDIF
      ENDDO
    ENDIF ! for MPI_id==0
    
    !> receive SR_B at each threads
    IF(MPI_id/=0) THEN
      IF(psi%cplx) THEN
        CALL MPI_Recv(psi%SR_B(1:psi%SR_B_length(MPI_id),1),psi%SR_B_length(MPI_id),   &
                       MPI_real_fortran,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
        CALL MPI_Recv(psi%SR_B(1:psi%SR_B_length(MPI_id),2),psi%SR_B_length(MPI_id),   &
                       MPI_real_fortran,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
      ELSE
        CALL MPI_Recv(psi%SR_B(1:psi%SR_B_length(MPI_id),1),psi%SR_B_length(MPI_id),   &
                       MPI_real_fortran,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
      ENDIF
    ENDIF

    ! free space
    IF(allocated(SR_B0)) deallocate(SR_B0)

  END SUBROUTINE SmolyakR_distribute_SRB_MPI
!=======================================================================================
#endif

#if(run_MPI)
!=======================================================================================
!> @brief distribute Smolyak terms to threads
!> Only the real part of the packed basis are converged to Smolyak rep. 
!> 1. transfer from packed basis to basis in Smolyak rep.
!> 2. transfer from basis Re. to grid rep. 
!> 3. grid Smolyak rep. are reserved for the current time step
!
!> if psi%clpx, the final SR_G will be (:,2) dim matrix, with real and imag. part resp.
!> otherwise, SR_G will be (:,1)
!=======================================================================================
  SUBROUTINE SmolyakR_distribute_SRG_MPI(psi,para_Op)
    USE mod_system
    USE mod_basis_set_alloc
    USE mod_param_SGType2
    USE mod_nDindex
    USE mod_psi_set_alloc,          ONLY:param_psi
    USE mod_Op,                     ONLY:param_Op
    USE mod_basis_BtoG_GtoB_SGType4,ONLY:tabPackedBasis_TO_tabR_AT_iG,                 &
                                         BDP_TO_GDP_OF_SmolyakRep,                     &
                                         getbis_tab_nq,getbis_tab_nb
    USE mod_MPI_Aid
    IMPLICIT NONE
    
    TYPE(param_psi),                intent(inout) :: psi
    TYPE(param_Op),                 intent(in)    :: para_Op
    
    TYPE(param_SGType2),pointer                   :: SGType2
    Real(kind=Rkind),allocatable                  :: SR_G0(:,:)
    Integer,pointer                               :: tab_l(:,:)
    Integer,allocatable                           :: tab_nq(:)
    Integer,allocatable                           :: tab_nb(:)
    Integer                                       :: iG
    Integer                                       :: nb0
    Integer                                       :: dim
    
    Character(len=*),parameter             :: name_sub='SmolyakR_distribute_SR_MPI'
    
    psi%SRG_MPI=.TRUE.

    SGType2 => para_Op%BasisnD%para_SGType2
    tab_l   => SGType2%nDind_SmolyakRep%Tab_nDval(:,:)
    
    dim=size(tab_l,1)
    CALL alloc_NParray(tab_nb,(/dim/),'tab_nb',name_sub)
    CALL alloc_NParray(tab_nq,(/dim/),'tab_nq',name_sub)

    !> allocate grid SR, the relevant index record the vector position for each iG
    !> psi%SR_G(psi%SR_G_index(iG):psi%SR_G_index(iG+1)-1) => SRG(iG)
    nb0=SGType2%nb0
    CALL allocate_array(psi%SR_G_length,0,MPI_np-1)
    CALL allocate_array(psi%SR_G_index,iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)+1)
    psi%SR_G_index(iGs_MPI(1,MPI_id))=1
    
    !IF(MPI_id==0) CALL allocate_array(psi%SR_G_index0,1,SGType2%nb_SG)
    
    DO i_MPI=0,MPI_np-1
      psi%SR_G_length(i_MPI)=0
      !IF(MPI_id==0) psi%SR_G_index0(iGs_MPI(1,i_MPI))=1
      DO iG=iGs_MPI(1,i_MPI),iGs_MPI(2,i_MPI)
        tab_nq(:)=getbis_tab_nq(tab_l(:,iG),para_Op%BasisnD%tab_basisPrimSG)
        psi%SR_G_length(i_MPI)=psi%SR_G_length(i_MPI)+product(tab_nq)*nb0
        !IF(MPI_id==0) psi%SR_G_index0(iG+1)=psi%SR_G_length(i_MPI)+1
        IF(i_MPI==MPI_id) psi%SR_G_index(iG+1)=psi%SR_G_length(MPI_id)+1 ! for each threads
      ENDDO
    ENDDO
    
!    DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)+1
!      write(*,*) 'iGs_MPI check: iG=',iG,psi%SR_G_index(iG),iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)+1
!    ENDDO
    
    !> allocate SR_G for each threads
    IF(psi%cplx) THEN
      CALL allocate_array(psi%SR_G,1,psi%SR_G_length(MPI_id),1,2)
    ELSE 
      CALL allocate_array(psi%SR_G,1,psi%SR_G_length(MPI_id),1,1)
    ENDIF

    !> distribute SR on grid to each threads
    IF(MPI_id==0) THEN
      DO i_MPI=1,MPI_np-1
      
        !> allocate temprary array for sending
        IF(psi%cplx) THEN
          CALL allocate_array(SR_G0,1,psi%SR_G_length(i_MPI),1,2)
        ELSE
          CALL allocate_array(SR_G0,1,psi%SR_G_length(i_MPI),1,1)
        ENDIF
        
        !> loop to send SR_G
        DO iG=iGs_MPI(1,i_MPI),iGs_MPI(2,i_MPI)
          tab_nq(:)=getbis_tab_nq(tab_l(:,iG),para_Op%BasisnD%tab_basisPrimSG)
          tab_nb(:)=getbis_tab_nb(tab_l(:,iG),para_Op%BasisnD%tab_basisPrimSG)

!          ! get SR on Basis
!          IF(psi%cplx) THEN
!            !-Real part-----------------------------------------------------------------
!            CALL tabPackedBasis_TO_tabR_AT_iG(psi%SR_B,Real(psi%CvecB,kind=Rkind),iG,  &
!                                              SGType2%tab_nb_OF_SRep(iG)*nb0)
!            ! transfer to SR on grid
!            CALL BDP_TO_GDP_OF_SmolyakRep(psi%SR_B,para_Op%BasisnD%tab_basisPrimSG,    &
!                                          tab_l(:,iG),tab_nq,tab_nb,nb0)
!            !< psi%SR_B is now actually psi%SR_G at iG
!            IF(size(psi%SR_B)-(psi%SR_G_index(iG+1)-psi%SR_G_index(iG))/=0)            &
!                              STOP 'error in SmolyakR_distribute_MPI, check psi%SR_B size'
!            SR_G0(psi%SR_G_index(iG):psi%SR_G_index(iG+1)-1,1)=psi%SR_B
!            
!            !-Imaginary part------------------------------------------------------------
!            CALL tabPackedBasis_TO_tabR_AT_iG(psi%SR_B,aimag(psi%CvecB),iG,            &
!                                              SGType2%tab_nb_OF_SRep(iG)*nb0)
!            ! transfer to SR on grid
!            CALL BDP_TO_GDP_OF_SmolyakRep(psi%SR_B,para_Op%BasisnD%tab_basisPrimSG,    &
!                                          tab_l(:,iG),tab_nq,tab_nb,nb0)
!            SR_G0(psi%SR_G_index(iG):psi%SR_G_index(iG+1)-1,2)=psi%SR_B
!          ELSE
!            CALL tabPackedBasis_TO_tabR_AT_iG(psi%SR_B,psi%RvecB,iG,                   &
!                                              SGType2%tab_nb_OF_SRep(iG)*nb0)
!            ! transfer to SR on grid
!            CALL BDP_TO_GDP_OF_SmolyakRep(psi%SR_B,para_Op%BasisnD%tab_basisPrimSG,    &
!                                        tab_l(:,iG),tab_nq,tab_nb,nb0)
!            SR_G0(psi%SR_G_index(iG):psi%SR_G_index(iG+1)-1,1)=psi%SR_B
!          ENDIF
          IF(psi%cplx) THEN
            ! Real part
            CALL pack_Basis_to_SR_grid_MPI(iG,Real(psi%CvecB,kind=Rkind),SR_G0(:,1),   &
                                       psi%SR_G_index,tab_nq,tab_nb,tab_l(:,iG),para_Op)
            ! Imaginary part
            CALL pack_Basis_to_SR_grid_MPI(iG,aimag(psi%CvecB),SR_G0(:,2),             &
                                       psi%SR_G_index,tab_nq,tab_nb,tab_l(:,iG),para_Op)
          ELSE
            CALL pack_Basis_to_SR_grid_MPI(iG,psi%RvecB,SR_G0(:,1),                    &
                                       psi%SR_G_index,tab_nq,tab_nb,tab_l(:,iG),para_Op)
          ENDIF
        ENDDO
      
        ! send SR_G
        IF(psi%cplx) THEN
          CALL MPI_Send(SR_G0(1:psi%SR_G_length(i_MPI),1),psi%SR_G_length(i_MPI),      &
                              MPI_real_fortran,root_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
          CALL MPI_Send(SR_G0(1:psi%SR_G_length(i_MPI),2),psi%SR_G_length(i_MPI),      &
                              MPI_real_fortran,root_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
        ELSE
          CALL MPI_Send(SR_G0(1:psi%SR_G_length(i_MPI),1),psi%SR_G_length(i_MPI),      &
                              MPI_real_fortran,root_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
        ENDIF
      ENDDO
      
      ! Smolyak rep. on MPI_id=0
      DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)
        tab_nq(:)=getbis_tab_nq(tab_l(:,iG),para_Op%BasisnD%tab_basisPrimSG)
        tab_nb(:)=getbis_tab_nb(tab_l(:,iG),para_Op%BasisnD%tab_basisPrimSG)
        
        IF(psi%cplx) THEN
          ! Real part
          CALL pack_Basis_to_SR_grid_MPI(iG,Real(psi%CvecB,kind=Rkind),psi%SR_G(:,1),  &
                                       psi%SR_G_index,tab_nq,tab_nb,tab_l(:,iG),para_Op)
          ! Imaginary part
          CALL pack_Basis_to_SR_grid_MPI(iG,aimag(psi%CvecB),psi%SR_G(:,2),            &
                                       psi%SR_G_index,tab_nq,tab_nb,tab_l(:,iG),para_Op)
        ELSE
          CALL pack_Basis_to_SR_grid_MPI(iG,psi%RvecB,psi%SR_G(:,1),                   &
                                       psi%SR_G_index,tab_nq,tab_nb,tab_l(:,iG),para_Op)
        ENDIF
      ENDDO
    ENDIF ! for MPI_id==0
    
    !> receive SR_G at each threads
    IF(MPI_id/=0) THEN
      IF(psi%cplx) THEN
        CALL MPI_Recv(psi%SR_G(1:psi%SR_G_length(MPI_id),1),psi%SR_G_length(MPI_id),   &
                       MPI_real_fortran,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
        CALL MPI_Recv(psi%SR_G(1:psi%SR_G_length(MPI_id),2),psi%SR_G_length(MPI_id),   &
                       MPI_real_fortran,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
      ELSE
        CALL MPI_Recv(psi%SR_G(1:psi%SR_G_length(MPI_id),1),psi%SR_G_length(MPI_id),   &
                       MPI_real_fortran,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
      ENDIF
    ENDIF

    ! free space
    IF(allocated(psi%SR_B)) deallocate(psi%SR_B)
    IF(allocated(SR_G0)) deallocate(SR_G0)

  END SUBROUTINE SmolyakR_distribute_SRG_MPI
!=======================================================================================
#endif

#if(run_MPI)
!=======================================================================================
!> transfer Smolyak rep. back to packed Basisi rep.
!=======================================================================================
  SUBROUTINE SmolyakR_to_packedB_SRB_MPI(psi,para_Op)
    USE mod_system
    USE mod_basis_set_alloc
    USE mod_param_SGType2
    USE mod_nDindex
    USE mod_SymAbelian,             ONLY:Calc_symab1_EOR_symab2
    USE mod_psi_set_alloc,          ONLY:param_psi
    USE mod_psi_Op,                 ONLY:Set_symab_OF_psiBasisRep
    USE mod_Op,                     ONLY:param_Op
    USE mod_basis_BtoG_GtoB_SGType4,ONLY:tabR_AT_iG_TO_tabPackedBasis,                 &
                                         GDP_TO_BDP_OF_SmolyakRep
    USE mod_MPI_Aid
    IMPLICIT NONE
    
    TYPE(param_psi),                intent(inout) :: psi
    TYPE(param_Op),                 intent(in)    :: para_Op
    
    TYPE(param_SGType2),pointer                   :: SGType2
    Real(kind=Rkind),allocatable                  :: SR_B0(:)
    Real(kind=Rkind),allocatable                  :: Rvec(:,:)
    Integer                                       :: nb0
    Integer                                       :: psi_symab
    Integer                                       :: iG

    Character(len=*),parameter             :: name_sub='SmolyakR_to_packedB_SRB_MPI'


    SGType2 => para_Op%BasisnD%para_SGType2
    nb0=SGType2%nb0

    ! send 
    IF(MPI_id/=0) THEN
      IF(psi%cplx) THEN
        CALL MPI_Send(psi%SR_B(1:psi%SR_B_length(MPI_id),1),psi%SR_B_length(MPI_id),   &
                      MPI_real_fortran,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_err)
        CALL MPI_Send(psi%SR_B(1:psi%SR_B_length(MPI_id),2),psi%SR_B_length(MPI_id),   &
                      MPI_real_fortran,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_err)
      ELSE
        CALL MPI_Send(psi%SR_B(1:psi%SR_B_length(MPI_id),1),psi%SR_B_length(MPI_id),   &
                      MPI_real_fortran,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_err)
      ENDIF
    ENDIF

    IF(MPI_id==0) THEN
      ! MPI_id=0 prat
      IF(psi%cplx) THEN
        CALL allocate_array(Rvec,1,size(psi%CvecB),1,2)
      ENDIF

      DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)
        IF(psi%cplx) THEN
          CALL SR_Basis_to_pack_Basis_MPI(iG,Rvec(:,1),psi%SR_B(:,1),psi%SR_B_index,   &
                                          para_Op)
          CALL SR_Basis_to_pack_Basis_MPI(iG,Rvec(:,2),psi%SR_B(:,2),psi%SR_B_index,   &
                                          para_Op)
        ELSE
          CALL SR_Basis_to_pack_Basis_MPI(iG,psi%RvecB,psi%SR_B(:,1),psi%SR_B_index,   &
                                          para_Op)
        ENDIF
      ENDDO
      IF(psi%cplx) psi%CvecB=cmplx(Rvec(:,1),Rvec(:,2),kind=Rkind)

      ! receive SRG from the other threads
      DO i_MPI=1,MPI_np-1
        CALL allocate_array(SR_B0,1,psi%SR_B_length(i_MPI))
        IF(psi%cplx) THEN
          CALL MPI_Recv(SR_B0(1:psi%SR_B_length(i_MPI)),psi%SR_B_length(i_MPI),      &
                        MPI_real_fortran,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_stat,MPI_err)
          DO iG=iGs_MPI(1,i_MPI),iGs_MPI(2,i_MPI)
            CALL SR_Basis_to_pack_Basis_MPI(iG,Rvec(:,1),SR_B0,psi%SR_B_index,          &
                                            para_Op)
          ENDDO
          
          CALL MPI_Recv(SR_B0(1:psi%SR_B_length(i_MPI)),psi%SR_B_length(i_MPI),      &
                        MPI_real_fortran,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_stat,MPI_err)
          DO iG=iGs_MPI(1,i_MPI),iGs_MPI(2,i_MPI)
            CALL SR_Basis_to_pack_Basis_MPI(iG,Rvec(:,2),SR_B0,psi%SR_B_index,          &
                                            para_Op)
          ENDDO
          psi%CvecB=cmplx(Rvec(:,1),Rvec(:,2),kind=Rkind)
        ELSE
          CALL MPI_Recv(SR_B0(1:psi%SR_B_length(i_MPI)),psi%SR_B_length(i_MPI),      &
                        MPI_real_fortran,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_stat,MPI_err)
          DO iG=iGs_MPI(1,i_MPI),iGs_MPI(2,i_MPI)
            CALL SR_Basis_to_pack_Basis_MPI(iG,psi%RvecB,SR_B0,psi%SR_B_index,          &
                                            para_Op)
          ENDDO
        ENDIF ! psi%cplx
      ENDDO ! for i_MPI=1,MPI_np-1

      psi_symab=Calc_symab1_EOR_symab2(para_Op%symab,psi%symab)
      CALL Set_symab_OF_psiBasisRep(psi,psi_symab)
    ENDIF ! for MPI_id==0

    IF(allocated(Rvec))      deallocate(Rvec)
    IF(allocated(SR_B0))     deallocate(SR_B0)
    IF(allocated(psi%SR_B))  deallocate(psi%SR_B)
    
    psi%SRB_MPI=.FALSE.

  ENDSUBROUTINE SmolyakR_to_packedB_SRB_MPI
!=======================================================================================
#endif

#if(run_MPI)
!=======================================================================================
!> transfer Smolyak rep. back to packed Basisi rep.
!=======================================================================================
  SUBROUTINE SmolyakR_to_packedB_SRG_MPI(psi,para_Op)
    USE mod_system
    USE mod_basis_set_alloc
    USE mod_param_SGType2
    USE mod_nDindex
    USE mod_SymAbelian,             ONLY:Calc_symab1_EOR_symab2
    USE mod_psi_set_alloc,          ONLY:param_psi
    USE mod_psi_Op,                 ONLY:Set_symab_OF_psiBasisRep
    USE mod_Op,                     ONLY:param_Op
    USE mod_basis_BtoG_GtoB_SGType4,ONLY:tabR_AT_iG_TO_tabPackedBasis,                 &
                                         GDP_TO_BDP_OF_SmolyakRep,                     &
                                         getbis_tab_nq,getbis_tab_nb
    USE mod_MPI_Aid
    IMPLICIT NONE
    
    TYPE(param_psi),                intent(inout) :: psi
    TYPE(param_Op),                 intent(in)    :: para_Op
    
    TYPE(param_SGType2),pointer                   :: SGType2
    Real(kind=Rkind),allocatable                  :: SR_G0(:)
    Real(kind=Rkind),allocatable                  :: Rvec(:,:)
    Integer,pointer                               :: tab_l(:,:)
    Integer,allocatable                           :: tab_nq(:)
    Integer,allocatable                           :: tab_nb(:)
    Integer                                       :: nb0
    Integer                                       :: psi_symab
    Integer                                       :: iG
    Integer                                       :: dim

    Character(len=*),parameter             :: name_sub='SmolyakR_to_packedB_SR_MPI'


    SGType2 => para_Op%BasisnD%para_SGType2
    tab_l   => SGType2%nDind_SmolyakRep%Tab_nDval(:,:)
    nb0=SGType2%nb0
    
    dim=size(tab_l,1)
    CALL alloc_NParray(tab_nb,(/dim/),'tab_nb',name_sub)
    CALL alloc_NParray(tab_nq,(/dim/),'tab_nq',name_sub)

    ! send 
    IF(MPI_id/=0) THEN
      IF(psi%cplx) THEN
        CALL MPI_Send(psi%SR_G(1:psi%SR_G_length(MPI_id),1),psi%SR_G_length(MPI_id),   &
                      MPI_real_fortran,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_err)
        CALL MPI_Send(psi%SR_G(1:psi%SR_G_length(MPI_id),2),psi%SR_G_length(MPI_id),   &
                      MPI_real_fortran,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_err)
      ELSE
        CALL MPI_Send(psi%SR_G(1:psi%SR_G_length(MPI_id),1),psi%SR_G_length(MPI_id),   &
                      MPI_real_fortran,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_err)
      ENDIF
    ENDIF

    IF(MPI_id==0) THEN
      ! MPI_id=0 prat
      IF(psi%cplx) THEN
        CALL allocate_array(Rvec,1,size(psi%CvecB),1,2)
      ENDIF

      DO iG=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)
        tab_nq(:)=getbis_tab_nq(tab_l(:,iG),para_Op%BasisnD%tab_basisPrimSG)
        tab_nb(:)=getbis_tab_nb(tab_l(:,iG),para_Op%BasisnD%tab_basisPrimSG)

        IF(psi%cplx) THEN
          CALL SR_grid_to_pack_Basis_MPI(iG,Rvec(:,1),psi%SR_G(:,1),psi%SR_G_index,    &
                                         tab_nq,tab_nb,tab_l(:,iG),para_Op)
          CALL SR_grid_to_pack_Basis_MPI(iG,Rvec(:,2),psi%SR_G(:,2),psi%SR_G_index,    &
                                         tab_nq,tab_nb,tab_l(:,iG),para_Op)
        ELSE
          CALL SR_grid_to_pack_Basis_MPI(iG,psi%RvecB,psi%SR_G(:,1),psi%SR_G_index,    &
                                         tab_nq,tab_nb,tab_l(:,iG),para_Op)
        ENDIF

!        IF(psi%cplx) THEN
!          CALL GDP_TO_BDP_OF_SmolyakRep(psi%SR_G(:,1),para_Op%BasisnD%tab_basisPrimSG,&
!                                        tab_l,tab_nq,tab_nb,nb0)
!          CALL GDP_TO_BDP_OF_SmolyakRep(psi%SR_G(:,2),para_Op%BasisnD%tab_basisPrimSG,&
!                                        tab_l,tab_nq,tab_nb,nb0)
!          CALL tabR_AT_iG_TO_tabPackedBasis(Rvec(:,1),psi%SR_G(:,1),iG,SGType2,       &
!                                            para_Op%BasisnD%WeightSG(iG))
!          CALL tabR_AT_iG_TO_tabPackedBasis(Rvec(:,2),psi%SR_G(:,2),iG,SGType2,       &
!                                            para_Op%BasisnD%WeightSG(iG))
!          psi%CvecB=cmplx(Rvec(:,1),Rvec(:,2))
!        ELSE
!          CALL GDP_TO_BDP_OF_SmolyakRep(psi%SR_G(:,1),para_Op%BasisnD%tab_basisPrimSG,&
!                                        tab_l,tab_nq,tab_nb,nb0)
!          CALL tabR_AT_iG_TO_tabPackedBasis(psi%RvecB,psi%SR_G(:,1),iG,SGType2,        &
!                                            para_Op%BasisnD%WeightSG(iG))
!        ENDIF
      ENDDO
      IF(psi%cplx) psi%CvecB=cmplx(Rvec(:,1),Rvec(:,2),kind=Rkind)

      ! receive SRG from the other threads
      DO i_MPI=1,MPI_np-1
        CALL allocate_array(SR_G0,1,psi%SR_G_length(i_MPI))
        IF(psi%cplx) THEN
          CALL MPI_Recv(SR_G0(1:psi%SR_G_length(i_MPI)),psi%SR_G_length(i_MPI),      &
                        MPI_real_fortran,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_stat,MPI_err)
          DO iG=iGs_MPI(1,i_MPI),iGs_MPI(2,i_MPI)
            CALL SR_grid_to_pack_Basis_MPI(iG,Rvec(:,1),SR_G0,psi%SR_G_index,          &
                                           tab_nq,tab_nb,tab_l(:,iG),para_Op)
          ENDDO
          
          CALL MPI_Recv(SR_G0(1:psi%SR_G_length(i_MPI)),psi%SR_G_length(i_MPI),      &
                        MPI_real_fortran,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_stat,MPI_err)
          DO iG=iGs_MPI(1,i_MPI),iGs_MPI(2,i_MPI)
            CALL SR_grid_to_pack_Basis_MPI(iG,Rvec(:,2),SR_G0,psi%SR_G_index,          &
                                           tab_nq,tab_nb,tab_l(:,iG),para_Op)
          ENDDO
          psi%CvecB=cmplx(Rvec(:,1),Rvec(:,2),kind=Rkind)
        ELSE
          CALL MPI_Recv(SR_G0(1:psi%SR_G_length(i_MPI)),psi%SR_G_length(i_MPI),      &
                        MPI_real_fortran,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_stat,MPI_err)
          DO iG=iGs_MPI(1,i_MPI),iGs_MPI(2,i_MPI)
            CALL SR_grid_to_pack_Basis_MPI(iG,psi%RvecB,SR_G0,psi%SR_G_index,          &
                                           tab_nq,tab_nb,tab_l(:,iG),para_Op)
          ENDDO
        ENDIF ! psi%cplx
      ENDDO ! for i_MPI=1,MPI_np-1

      psi_symab=Calc_symab1_EOR_symab2(para_Op%symab,psi%symab)
      CALL Set_symab_OF_psiBasisRep(psi,psi_symab)
    ENDIF ! for MPI_id==0

    IF(allocated(Rvec)) deallocate(Rvec)
    IF(allocated(SR_G0)) deallocate(SR_G0)
    IF(allocated(tab_nq)) deallocate(tab_nq)
    IF(allocated(tab_nb)) deallocate(tab_nb)

  ENDSUBROUTINE SmolyakR_to_packedB_SRG_MPI
!=======================================================================================
#endif

#if(run_MPI)
!=======================================================================================
!> @brief get Smolyak rep on Basis from packed basis
!=======================================================================================
  SUBROUTINE pack_Basis_to_SR_Basis_MPI(iG,packB,SR_B,SR_B_index,para_Op)
    USE mod_system
    USE mod_basis_set_alloc
    USE mod_basis_BtoG_GtoB_SGType4,ONLY:tabPackedBasis_TO_tabR_AT_iG
    USE mod_Op,                     ONLY:param_Op
    USE mod_MPI_Aid
    IMPLICIT NONE

    TYPE(param_Op),                 intent(in)    :: para_Op
    Real(kind=Rkind),               intent(in)    :: packB(:)
    Real(kind=Rkind),               intent(inout) :: SR_B(:)
    Integer,                        intent(in)    :: SR_B_index(:)
    Integer,                        intent(in)    :: iG

    Real(kind=Rkind),allocatable                  :: SRB(:)

    ! transfer to SR on Basis
    CALL tabPackedBasis_TO_tabR_AT_iG(SRB,packB,iG,para_Op%BasisnD%para_SGType2)

    !< store in psi%SR_B
    IF(size(SRB)-(SR_B_index(iG+1)-SR_B_index(iG))/=0)                                &
                         STOP 'error in SmolyakR_distribute_SR_MPI, check psi%SR_B size'
    SR_B(SR_B_index(iG):SR_B_index(iG+1)-1)=SRB

  ENDSUBROUTINE pack_Basis_to_SR_Basis_MPI
!=======================================================================================
#endif

#if(run_MPI)
!=======================================================================================
!> @brief get Smolyak rep on grid from packed basis
!=======================================================================================
  SUBROUTINE pack_Basis_to_SR_grid_MPI(iG,packB,SR_G,SR_G_index,                  &
                                       tab_nq,tab_nb,tab_l,para_Op)
    USE mod_system
    USE mod_basis_set_alloc
    USE mod_param_SGType2
    USE mod_basis_BtoG_GtoB_SGType4,ONLY:tabPackedBasis_TO_tabR_AT_iG,                 &
                                         BDP_TO_GDP_OF_SmolyakRep
    USE mod_Op,                     ONLY:param_Op
    USE mod_MPI_Aid
    IMPLICIT NONE

    TYPE(param_Op),                 intent(in)    :: para_Op
    Real(kind=Rkind),               intent(inout) :: SR_G(:)
    Real(kind=Rkind),               intent(in)    :: packB(:)
    Integer,                        intent(in)    :: SR_G_index(:)
    Integer,allocatable,            intent(in)    :: tab_nq(:)
    Integer,allocatable,            intent(in)    :: tab_nb(:)
    Integer,                        intent(in)    :: tab_l(:)
    Integer,                        intent(in)    :: iG

    Real(kind=Rkind),allocatable                  :: SR_B(:)
    TYPE(param_SGType2),pointer                   :: SGType2

    SGType2 => para_Op%BasisnD%para_SGType2
  
    CALL tabPackedBasis_TO_tabR_AT_iG(SR_B,packB,iG,SGType2)
    ! transfer to SR on grid
    CALL BDP_TO_GDP_OF_SmolyakRep(SR_B,para_Op%BasisnD%tab_basisPrimSG,                &
                                  tab_l,tab_nq,tab_nb,SGType2%nb0)
    !< psi%SR_B is now actually psi%SR_G at iG
    IF(size(SR_B)-(SR_G_index(iG+1)-SR_G_index(iG))/=0)                                &
                         STOP 'error in SmolyakR_distribute_SR_MPI, check psi%SR_B size'
    SR_G(SR_G_index(iG):SR_G_index(iG+1)-1)=SR_B

  ENDSUBROUTINE pack_Basis_to_SR_grid_MPI
!=======================================================================================
#endif

#if(run_MPI)
!=======================================================================================
!> @brief transfer Smolyak rep on grid to packed basis
!=======================================================================================
  SUBROUTINE SR_Basis_to_pack_Basis_MPI(iG,packB,SR_B,SR_B_index,para_Op)
    USE mod_system
    USE mod_basis_set_alloc
    USE mod_param_SGType2
    USE mod_basis_BtoG_GtoB_SGType4,ONLY:tabR_AT_iG_TO_tabPackedBasis
    USE mod_Op,                     ONLY:param_Op
    USE mod_MPI_Aid
    IMPLICIT NONE
    
    
    TYPE(param_Op),                 intent(in)    :: para_Op
    Real(kind=Rkind),               intent(in)    :: SR_B(:)
    Real(kind=Rkind),               intent(inout) :: packB(:)
    Integer,                        intent(in)    :: SR_B_index(:)
    Integer,                        intent(in)    :: iG
    
    Real(kind=Rkind),allocatable                  :: SRB(:)
    Integer                                       :: nb0


    nb0=para_Op%BasisnD%para_SGType2%nb0

    ! the SR_B at iG is stored in SRB
    CALL allocate_array(SRB,1,SR_B_index(iG+1)-SR_B_index(iG))
    SRB=SR_B(SR_B_index(iG):SR_B_index(iG+1)-1)

    ! now SR_B stores the basis rep. 
    CALL tabR_AT_iG_TO_tabPackedBasis(packB,SRB,iG,para_Op%BasisnD%para_SGType2,       &
                                      para_Op%BasisnD%WeightSG(iG))
    IF(allocated(SRB)) deallocate(SRB)

  ENDSUBROUTINE SR_Basis_to_pack_Basis_MPI
!=======================================================================================
#endif

#if(run_MPI)
!=======================================================================================
!> @brief transfer Smolyak rep on grid to packed basis
!=======================================================================================
  SUBROUTINE SR_grid_to_pack_Basis_MPI(iG,packB,SR_G,SR_G_index,                       &
                                       tab_nq,tab_nb,tab_l,para_Op)
    USE mod_system
    USE mod_basis_set_alloc
    USE mod_param_SGType2
    USE mod_basis_BtoG_GtoB_SGType4,ONLY:tabR_AT_iG_TO_tabPackedBasis,                 &
                                         GDP_TO_BDP_OF_SmolyakRep
    USE mod_Op,                     ONLY:param_Op
    USE mod_MPI_Aid
    IMPLICIT NONE
    
    
    TYPE(param_Op),                 intent(in)    :: para_Op
    Real(kind=Rkind),               intent(in)    :: SR_G(:)
    Real(kind=Rkind),               intent(inout) :: packB(:)
    Integer,                        intent(in)    :: SR_G_index(:)
    Integer,allocatable,            intent(in)    :: tab_nq(:)
    Integer,allocatable,            intent(in)    :: tab_nb(:)
    Integer,                        intent(in)    :: tab_l(:)
    Integer,                        intent(in)    :: iG
    
    Real(kind=Rkind),allocatable                  :: SR_temp(:)
    Integer                                       :: nb0

    
    nb0=para_Op%BasisnD%para_SGType2%nb0
    
    ! the SR_G at iG is stored in SR_B
    CALL allocate_array(SR_temp,1,SR_G_index(iG+1)-SR_G_index(iG))
    SR_temp=SR_G(SR_G_index(iG):SR_G_index(iG+1)-1)
    
    CALL GDP_TO_BDP_OF_SmolyakRep(SR_temp,para_Op%BasisnD%tab_basisPrimSG,             &
                                  tab_l,tab_nq,tab_nb,nb0)
    ! now SR_B stores the basis rep. 
    CALL tabR_AT_iG_TO_tabPackedBasis(packB,SR_temp,iG,para_Op%BasisnD%para_SGType2,   &
                                      para_Op%BasisnD%WeightSG(iG))
    IF(allocated(SR_temp)) deallocate(SR_temp)
  ENDSUBROUTINE
!=======================================================================================
#endif

#if(run_MPI)
!=======================================================================================
!> calculate psi1=psi1-C*psi2
!> e.g. w1=w1-H(k,k-1)*tab_KrylovSpace(k-1) !> H|v_{k}>-beta_{k-1}|v_{k-1}>
!>      w1=w1-Overlap*tab_KrylovSpace(k)
!=======================================================================================
   SUBROUTINE w1_minus_Cv_SR_MPI(psi1,Cplx,psi2)
     USE mod_psi_set_alloc,          ONLY:param_psi
     IMPLICIT NONE
     
     TYPE(param_psi),                intent(inout) :: psi1
     TYPE(param_psi),                intent(in)    :: psi2
     Complex(kind=Rkind),            intent(in)    :: Cplx
     
     Real(kind=Rkind)                              :: R
     Real(kind=Rkind)                              :: C


     R=Real(Cplx,kind=Rkind)
     C=aimag(Cplx)

     IF(psi1%SRG_MPI) THEN
       psi1%SR_G(:,1)=psi1%SR_G(:,1)-(psi2%SR_G(:,1)*R-psi2%SR_G(:,2)*C)
       psi1%SR_G(:,2)=psi1%SR_G(:,2)-(psi2%SR_G(:,1)*C+psi2%SR_G(:,2)*R)
     ELSE IF(psi1%SRB_MPI) THEN
       psi1%SR_B(:,1)=psi1%SR_B(:,1)-(psi2%SR_B(:,1)*R-psi2%SR_B(:,2)*C)
       psi1%SR_B(:,2)=psi1%SR_B(:,2)-(psi2%SR_B(:,1)*C+psi2%SR_B(:,2)*R)
     ENDIF

     IF(psi1%symab/=psi2%symab) THEN
       IF(psi1%symab==-2) THEN
         psi1%symab=psi2%symab
       ELSE
         psi1%symab=-1
       ENDIF
     ENDIF
     
   ENDSUBROUTINE w1_minus_Cv_SR_MPI
!=======================================================================================
#endif

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

      psi_save = psi

      DO icheb=1,max_cheby

#if(run_MPI)
        w1  = psi
        IF(para_propa%once_Hmin) THEN
          CALL sub_OpPsi(w1,w2,para_H) ! calculate once for Hmax
          CALL sub_Hmax(para_propa,para_H)
          para_propa%once_Hmin=.FALSE.
                
          para_propa%Hmax = para_propa%Hmax + para_propa%para_poly%DHmax
          para_propa%para_poly%Hmin = para_propa%Hmin
          para_propa%para_poly%Hmax = para_propa%Hmax

          CALL initialisation1_poly(para_propa%para_poly,                                &
                                    para_propa%WPdeltaT,                                 &
                                    para_propa%type_WPpropa)
          
          para_H%scaled = .TRUE.
          para_H%E0     = para_propa%para_poly%E0
          para_H%Esc    = para_propa%para_poly%Esc
          write(out_unitp,*) 'Hmin,Hmax check:', para_propa%Hmin,para_propa%Hmax
        ENDIF
#endif

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

        ! calculate cofficients for Chebychev in para_propa%para_poly%coef_poly
        !r = max(ONE,HALF * para_propa%para_poly%deltaE * para_propa%WPdeltaT)
        ! possible a problem here for very small system
        r = HALF * para_propa%para_poly%deltaE * para_propa%WPdeltaT
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
#if(run_MPI)
        !w1  = psi
#else
        w1  = psi
#endif
        psi = psi * para_propa%para_poly%coef_poly(1)

        psi0Hkpsi0(1) =  psi0Hkpsi0(0)

!     - The second term of the expansion -----------------
        write(out_unitp,'(a)',advance='no') 'cheby rec:'
        CALL sub_OpPsi(w1,w2,para_H)
        CALL sub_scaledOpPsi(w1,w2,para_H%E0,para_H%Esc) ! limited to MPI_id==0 in subroutine

        IF(MPI_id==0) THEN
          w2  = w2 * rt
          psi = psi + w2 * para_propa%para_poly%coef_poly(2)
        ENDIF
        
        psi0Hkpsi0(2) = Calc_AutoCorr(psi0,w2,para_propa,T,Write_AC=.FALSE.)

!     - The higher terms of the expansion ----------------

        DO jt=3,para_propa%para_poly%npoly

          CALL sub_OpPsi(w2,w3,para_H)
          CALL sub_scaledOpPsi(w2,w3,para_H%E0,para_H%Esc)
          IF (mod(jt,100) == 0) write(out_unitp,'(a)',advance='no') '.'
          
          IF(MPI_id==0) THEN
!           Recurrence relations of the Chebychev expansion:
            w3 = w1 + w3 * rt2
            w1 = w2
            w2 = w3
            psi = psi + w2 * para_propa%para_poly%coef_poly(jt)

            CALL norm2_psi(w2)

            norm_exit = abs(w2%norm2*para_propa%para_poly%coef_poly(jt))
          ENDIF
#if(run_MPI)
          CALL MPI_Bcast(norm_exit,size1_MPI,MPI_Real8,root_MPI,MPI_COMM_WORLD,MPI_err)
#endif
          jt_exit = jt
          IF (debug) write(out_unitp,*) 'jt,norms',jt,norm_exit

          IF (norm_exit > TEN**15) THEN
            write(out_unitp,*) ' WARNING: Norm^2 of the vector is TOO large (> 10^15)',jt,norm_exit
            IF(MPI_id==0) para_propa%march_error = .TRUE.
            EXIT
          END IF

          IF(MPI_id==0) THEN
            psi0Hkpsi0(jt) = Calc_AutoCorr(psi0,w2,para_propa,T,Write_AC=.FALSE.)
#if(run_MPI)
            IF (norm_exit < para_propa%para_poly%poly_tol) exitall=.TRUE.
#else
            IF (norm_exit < para_propa%para_poly%poly_tol) EXIT
#endif
          ENDIF
#if(run_MPI)
          !> MPI the other threads are waiting for master here
          CALL MPI_Bcast(exitall,size1_MPI,MPI_Real8,root_MPI,MPI_COMM_WORLD,MPI_err)
          IF(exitall) EXIT
#endif
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

          psi = psi_save
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
!=======================================================================================      

      SUBROUTINE march_cheby_old(T,no,psi,psi0,para_H,para_propa)
      USE mod_system
      USE mod_psi,   ONLY : param_psi,ecri_psi,norm2_psi
      USE mod_Op,    ONLY : param_Op, sub_OpPsi,sub_scaledOpPsi
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
         IF (norm_exit < para_propa%para_poly%poly_tol) EXIT

         psi0Hkpsi0(jt) = Calc_AutoCorr(psi0,w2,para_propa,T,Write_AC=.FALSE.)

      END DO
      write(out_unitp,*) 'jt_exit,norms',jt_exit,abs(w2%norm2),norm_exit

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


      w1 = psi
      DO j=1,para_propa%para_poly%max_poly

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

      CALL alloc_NParray(tab_dnE,(/nb_der,3/),                            &
                      "tab_dnE","march_new0_noD_field",(/0,1/))
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

      CALL alloc_NParray(tab_dnE,(/nb_der,3/),                          &
                      "tab_dnE",name_sub,(/0,1/))
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
         CALL flush_perso(out_unitp)
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

