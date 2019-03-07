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

MODULE mod_FullPropa
USE mod_Constant
IMPLICIT NONE

PRIVATE
PUBLIC :: sub_propagation,sub_propagation3

CONTAINS

!================================================================
!     WP propagation
!
! para_AllOp : table of operator
!
!     para_H    = para_AllOp%tab_Op(1)
!     para_Dipx = para_AllOp%tab_Op(2)
!     para_Dipy = para_AllOp%tab_Op(3)
!     para_Dipz = para_AllOp%tab_Op(4)
!
!     0 : spectral propagation (WP is projected on eigenvectors)
!     1 : Cheby    propagation
!     2 : nOD      propagation
!     22: nOD      propagation + field
!     3 : relaxation
!     33: relaxation (n states)
!
!
!================================================================
      SUBROUTINE sub_propagation(WP0,para_AllOp,para_propa)
      USE mod_system
      USE mod_Op
      USE mod_propa
      USE mod_psi_set_alloc
      USE mod_field
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa)           :: para_propa
      TYPE (param_psi), intent(in) :: WP0(1)

!----- variables for the operators ----------------------------
      TYPE (param_AllOp), target :: para_AllOp

!----- variables for H ---------------------------------------------
      TYPE (param_Op), pointer   :: para_H
      complex (kind=Rkind)       :: Et
      integer                    :: nb_diago,max_diago
      TYPE (param_psi),pointer   :: psi(:)
      real (kind=Rkind),pointer  :: Ene0(:)

!----- for printing --------------------------------------------------
      logical :: print_Op


!------ working variables ---------------------------------
      integer             :: i,j,jt
      real (kind=Rkind)   :: dnE,T
      real (kind=Rkind)   :: w,wmin,wmax,stepw
      integer             :: iOp,n_Op,iDip
      TYPE (param_psi)    :: WP(1)

      !logical :: SGtype4 = .FALSE.
      logical :: SGtype4 = .TRUE.
      logical :: direct_KEO

!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
      para_H => para_AllOp%tab_Op(1)
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING propagation'
        write(out_unitp,*) 'n',WP0(1)%nb_tot
        write(out_unitp,*)
        write(out_unitp,*) 'WP0 BasisRep'
        CALL ecri_psi(psi=WP0(1))
        write(out_unitp,*)
        CALL write_param_Op(para_H)
      END IF
!-----------------------------------------------------------


      para_propa%Hmax = para_propa%Hmax + para_propa%para_poly%DHmax

      para_propa%para_poly%Hmin = para_propa%Hmin
      para_propa%para_poly%Hmax = para_propa%Hmax



      write(out_unitp,*) 'Tmax,DeltaT (ua)=> ',                         &
             para_propa%WPTmax,para_propa%WPdeltaT
      write(out_unitp,*) 'Tmax,DeltaT (fs)=> ',                         &
               para_propa%WPTmax*get_Conv_au_TO_unit('t','fs'),         &
             para_propa%WPdeltaT*get_Conv_au_TO_unit('t','fs')
      write(out_unitp,*) 'Tmax,DeltaT (ps)=> ',                         &
               para_propa%WPTmax*get_Conv_au_TO_unit('t','ps'),         &
             para_propa%WPdeltaT*get_Conv_au_TO_unit('t','ps')
      write(out_unitp,*) '... DeltaE,Emax (cm-1)',                      &
          TWO*pi/para_propa%WPTmax   * get_Conv_au_TO_unit('E','cm-1'), &
          TWO*pi/para_propa%WPdeltaT * get_Conv_au_TO_unit('E','cm-1')

      direct_KEO = para_H%para_ReadOp%para_FileGrid%Save_MemGrid
      direct_KEO = direct_KEO .AND. para_H%BasisnD%dnGGRep
      direct_KEO = direct_KEO .AND. (para_H%type_Op == 10)
      direct_KEO = direct_KEO .AND. para_H%direct_KEO

      SGtype4    = SGtype4 .AND. (para_H%BasisnD%SparseGrid_type == 4)

      IF (abs(para_propa%type_WPpropa) /= 33 .AND. abs(para_propa%type_WPpropa) /= 34) THEN
        WP(1) = WP0(1)
      END IF

      IF (para_propa%n_WPecri < 1) para_propa%n_WPecri = 2 * int(para_propa%WPTmax/para_propa%WPdeltaT) !nothing written

      SELECT CASE (para_propa%type_WPpropa)

      CASE (1,2,5,6,7)

        IF (SGtype4 .AND. direct_KEO) THEN
          !CALL sub_propagation11_SG4(WP0,WP,1,para_H,para_propa)
          CALL sub_propagation11(WP0,WP,1,para_H,para_propa)
        ELSE
          CALL sub_propagation11(WP0,WP,1,para_H,para_propa)
        END IF

        CALL TF_autocorr(para_propa)

      CASE (-3,3)

        CALL sub_propagation3(Et,WP0(1),WP(1),para_H,para_propa)

      CASE (34)

        nb_diago = min(para_propa%max_ana,para_H%nb_tot)
        nullify(psi)
        CALL alloc_array(psi,(/nb_diago/),"psi","sub_propagation")
        nullify(Ene0)
        CALL alloc_array(Ene0,(/nb_diago/),"Ene0","sub_propagation")

        CALL sub_propagation34(psi,Ene0,nb_diago,                       &
                               para_H,para_propa)

        CALL dealloc_array(psi,"psi","sub_propagation")
        CALL dealloc_array(Ene0,"Ene0","sub_propagation")

      CASE (22,24,50,52,54,221,222,223)

!       - for initialization of field variables -----------
!       - enable to do loop of the field ------------------
        print_Op=.TRUE.
        CALL init0_field(para_propa%para_field,para_propa%WPTmax)
        CALL read_field(para_propa%para_field)

!       - dipole moment on the BasisRep basis ---------------------
        iOp = 3
        DO i=iOp,iOp+2
          iDip = para_AllOp%tab_Op(i)%n_Op
          IF (para_propa%para_field%pola_xyz(iDip)) THEN
            write(out_unitp,*) 'Propagation with ',                     &
                   trim(para_AllOp%tab_Op(i)%name_Op),                  &
                   para_AllOp%tab_Op(i)%n_Op
          END IF
        END DO

        IF (para_propa%para_field%stepw == 0) THEN
          write(out_unitp,*) 'propagation without scan in w'
          WP(1) = WP0(1)
          CALL sub_propagation24(WP,1,print_Op,                         &
                                 para_propa%para_field,.FALSE.,         &
                                 para_H,para_AllOp%tab_Op(iOp:iOp+2),   &
                                 para_propa)
        ELSE

          write(out_unitp,*) 'propagation with scan in w'
          wmin  = para_propa%para_field%wmin
          wmax  = para_propa%para_field%wmax
          stepw = para_propa%para_field%stepw
          w     = wmin - stepw
          print_Op = .FALSE.
          !print_Op = .TRUE.

          DO WHILE (w < wmax)
            w = w + stepw

            para_propa%para_field%w(:,:) = w

            write(out_unitp,*) 'propagation with w =',w
            WP = WP0
            CALL sub_propagation24(WP,1,print_Op,                       &
                                   para_propa%para_field,.FALSE.,       &
                                   para_H,para_AllOp%tab_Op(iOp:iOp+2), &
                                   para_propa)

          END DO
        END IF

      CASE (100)

        CALL sub_propagation100(WP0,1,print_Op,                         &
                                para_propa%para_field,.FALSE.,          &
                                para_AllOp,para_propa)

      CASE DEFAULT

        write(out_unitp,*) ' sub_propagation : NO spectral propagation'
        STOP

      END SELECT
      write(out_unitp,*) 'Number of Hamiltonian operations (H I psi >)',para_H%nb_OpPsi



      CALL dealloc_psi(WP(1),delete_all=.TRUE.)
      nullify(para_H)

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END propagation'
      END IF
!-----------------------------------------------------------
      end subroutine sub_propagation
!================================================================
!
!     3 : propagation in imaginary time => ground state
!
!================================================================
      SUBROUTINE sub_propagation3(E0,psi0,psi,para_H,para_propa)
      USE mod_system
      USE mod_psi_B_TO_G
      USE mod_Op
      USE mod_propa
      USE mod_march
      USE mod_psi_set_alloc
      USE mod_ana_psi
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H
      TYPE (param_psi), intent(in)   :: psi0

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa)   :: para_propa
      TYPE (param_psi)     :: psi,w1,w2
      TYPE (param_ana_psi) :: ana_psi

!------ working parameters --------------------------------
      complex (kind=Rkind) :: cdot

      integer       :: it,no
      integer       :: i,j
      integer       :: max_ecri
      real (kind=Rkind) :: T      ! time
      real (kind=Rkind) :: T_Delta! time+deltaT
      complex (kind=Rkind) :: E0,E1
      real (kind=Rkind) :: DeltaE,epsi,RE0

      logical       :: FOD
      integer  ::   nioWP
      logical :: BasisRep,GridRep



!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING sub_propagation3'
        write(out_unitp,*) 'Tmax,deltaT',para_propa%WPTmax,para_propa%WPdeltaT
        write(out_unitp,*) 'Hmin,Hmax',para_propa%para_poly%Hmin,               &
                                para_propa%para_poly%Hmax
        write(out_unitp,*)
        write(out_unitp,*) 'nb_ba,nb_qa',psi%nb_ba,psi%nb_qa
        write(out_unitp,*) 'nb_bi',psi%nb_bi
        write(out_unitp,*)

        write(out_unitp,*) 'psiBasisRep'
        CALL ecri_psi(psi=psi,ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)
        write(out_unitp,*) 'psiGridRep'
        CALL ecri_psi(psi=psi,ecri_GridRep=.TRUE.,ecri_BasisRep=.FALSE.)

       END IF
!-----------------------------------------------------------

      BasisRep = psi0%BasisRep
      GridRep  = psi0%GridRep

!-----------------------------------------------------------
!     copy psi0 in psi
      w1  = psi0
      w2  = psi0
      psi = psi0
!-----------------------------------------------------------

!-----------------------------------------------------------
      write(out_unitp,*) ' vib : propagation: ',para_propa%name_WPpropa

!     - parameters for poly (cheby and nOD) ... ------------
      CALL initialisation1_poly(para_propa%para_poly,                   &
                                para_propa%WPdeltaT,                    &
                                para_propa%type_WPpropa)

!     - scaling of H ---------------------------------------
      para_H%scaled = .TRUE.
      para_H%E0     = para_propa%para_poly%E0
      para_H%Esc    = para_propa%para_poly%Esc
!-----------------------------------------------------------


!-----------------------------------------------------------
      epsi = para_propa%para_poly%poly_tol
      T    = ZERO
      it   = 0

      E1 = ZERO
      CALL sub_PsiOpPsi(E0,psi,w1,para_H)
      E1 = E0
      DeltaE = abs(para_H%Hmax - para_H%Hmin)
      FOD = .NOT. (para_propa%type_WPpropa == -3) .AND. .NOT. para_H%cplx
      IF (para_propa%write_iter .OR. debug) THEN
        write(out_unitp,21) T,E0*get_Conv_au_TO_unit('E','cm-1')
      END IF
 21   format('ImTimeProp ',f12.2,' (',2(1x,f18.4),')')

!-----------------------------------------------------------

!------- propagation loop ---------------------------------
      DO WHILE (T <= para_propa%WPTmax .AND. DeltaE > epsi)

!       --------------------------------------------------------
        IF (FOD) THEN
          IF (para_propa%write_iter .OR. debug) write(out_unitp,*) 'march FOD'
          para_H%E0     = ZERO
          E1 = E0
          CALL march_FOD_Opti_im(psi,RE0,T,it,para_H,para_propa)
!         E0 = RE0
!  normalement, RE0 doit etre l'energie....

          CALL renorm_psi(psi)
          E1 = E0
          CALL sub_PsiOpPsi(E0,psi,w1,para_H)
          DeltaE = abs(E1-E0)
          FOD = (DeltaE > ONETENTH**6)  ! about 3 cm-1
        ELSE
          IF (para_propa%write_iter .OR. debug) write(out_unitp,*) 'march nOD'
          para_H%E0     = para_propa%para_poly%E0
          IF (para_propa%type_WPpropa == 3) para_H%E0     = ZERO
          CALL march_nOD_im(T,no,psi,psi0,w1,w2,para_H,para_propa)

          CALL renorm_psi(psi)
          E1 = E0
          CALL sub_PsiOpPsi(E0,psi,w1,para_H)
          DeltaE = E0-E1
          IF (DeltaE*real(para_propa%type_WPpropa,kind=Rkind) > ZERO) THEN
            para_propa%WPdeltaT = para_propa%WPdeltaT/TWO
          END IF

          DeltaE = abs(DeltaE)
          T  = T + para_propa%WPdeltaT
        END IF

        IF (para_propa%write_iter .OR. debug) THEN
           write(out_unitp,21) T,E0*get_Conv_au_TO_unit('E','cm-1')
           !write(out_unitp,21) T,E0
        END IF
        it = it + 1

        CALL flush_perso(out_unitp)

      END DO
!----------------------------------------------------------


!----------------------------------------------------------
!     - write the final WP---------------------------------
      IF (T > para_propa%WPTmax)                                        &
          write(out_unitp,*) ' WARNING : the WP is not fully relaxed'

      para_propa%ana_psi%Write_psi2_Grid = .FALSE.
      para_propa%ana_psi%Write_psi_Grid  = .TRUE.
      para_propa%ana_psi%Write_psi_Basis = .TRUE.

      !CALL Write_ana_psi(para_propa%ana_psi)
      CALL sub_analyze_WP_OpWP(T,(/ psi /),1,para_H,para_propa,adia=.FALSE.)

      IF (debug .OR. psi%nb_tot < 1000) THEN
        write(out_unitp,*) 'WP (BasisRep) at T=',T
        CALL ecri_psi(T=T,psi=psi,ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)
      END IF

!----------------------------------------------------------
      write(out_unitp,*)
      write(out_unitp,*) 'Number of Hamiltonian operation',para_H%nb_OpPsi
      write(out_unitp,*)
      write(out_unitp,*) 'relaxed E/pot0 (ua)',E0
      write(out_unitp,*) 'relaxed E (ua)',E0+para_H%pot0
      write(out_unitp,*) 'relaxed En/pot0 (cm-1)',                              &
                                              E0*get_Conv_au_TO_unit('E','cm-1')
      write(out_unitp,*)
      write(out_unitp,*) 'DHmax (ua)',E0 -                                      &
               para_propa%para_poly%Hmax + para_propa%para_poly%DHmax
!----------------------------------------------------------


      CALL dealloc_psi(w1,delete_all=.TRUE.)
      CALL dealloc_psi(w2,delete_all=.TRUE.)

!----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'END sub_propagation3'
       END IF
!----------------------------------------------------------


      END SUBROUTINE sub_propagation3

      SUBROUTINE sub_propagation34(psi,Ene0,nb_diago,                   &
                                   para_H,para_propa)
      USE mod_system
      USE mod_Op
      !USE mod_psi
      USE mod_psi_set_alloc
      USE mod_psi_B_TO_G
      USE mod_psi_SimpleOp
      USE mod_ana_psi
      USE mod_psi_Op
      USE mod_propa
      USE mod_march
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa

      integer            :: nb_diago,nb_diagoR
      TYPE (param_psi)   :: psi(nb_diago)
      real (kind=Rkind)      :: Ene0(nb_diago)
      TYPE (param_psi)   :: Hpsi,H2psi,g
      logical            :: cplx
      complex (kind=Rkind)   :: CS
      real (kind=Rkind)      :: RS


!------ working parameters --------------------------------
      complex (kind=Rkind) :: CEne0,cdot

      integer       :: it,it_all,no
      integer       :: i,ii,j,iqa
      integer       :: max_ecri
      real (kind=Rkind) :: DeltaT,T      ! time
      real (kind=Rkind) :: DeltaE,Deltapsi,epsi,normeg,th
      real (kind=Rkind) :: avH1,avH2,avH3,A,B,C,D,DT1,DT2,S,E1,E2
      real (kind=Rkind) :: Qact(para_H%mole%nb_act1)
      real (kind=Rkind) :: psi_q(nb_diago)
      complex (kind=Rkind)   :: CavH1,CavH2,CavH3

      logical       :: lect = .FALSE.

      integer  ::   nioWP


!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING sub_propagation34'
        write(out_unitp,*) 'Tmax,deltaT',para_propa%WPTmax,para_propa%WPdeltaT
        write(out_unitp,*) 'Hmin,Hmax',para_propa%para_poly%Hmin,               &
                                para_propa%para_poly%Hmax
        write(out_unitp,*)
       END IF
!-----------------------------------------------------------

!-----------------------------------------------------------
      write(out_unitp,*) ' propagation34: ',para_propa%name_WPpropa
      write(out_unitp,*) ' nb_diago',nb_diago
      write(out_unitp,*)

!     - scaling of H ---------------------------------------
      para_H%E0     = ZERO
      para_H%Esc    = ONE
      para_H%scaled = .TRUE.
!-----------------------------------------------------------


!------ initialization -------------------------------------
      epsi = para_propa%para_poly%poly_tol

      cplx = para_H%cplx
      CALL init_psi(Hpsi,para_H,cplx)
      CALL alloc_psi(Hpsi)
      CALL init_psi(H2psi,para_H,cplx)
      CALL alloc_psi(H2psi)
      CALL init_psi(g,para_H,cplx)
      CALL alloc_psi(g)


      IF (lect) THEN
        read(99,*) nb_diagoR
        write(out_unitp,*) ' nb_diagoR',nb_diagoR
        nb_diago = min(nb_diagoR,nb_diago)
        write(out_unitp,*) ' new nb_diago',nb_diago
        DO i=1,nb_diago
          CALL init_psi(psi(i),para_H,cplx)
          psi(i)%GridRep=.TRUE.
          CALL alloc_psi(psi(i))
        END DO

        IF (psi(1)%cplx) THEN
          DO iqa=1,psi(1)%nb_qa
            read(99,*) Qact,psi_q(1:nb_diago)
            DO i=1,nb_diago
              psi(i)%CvecG(iqa) = psi_q(i)
            END DO
          END DO
        ELSE
          DO iqa=1,psi(1)%nb_qa
            read(99,*) Qact,psi_q(1:nb_diago)
            DO i=1,nb_diago
              psi(i)%RvecG(iqa) = psi_q(i)
           END DO
          END DO
        END IF


        DO i=1,nb_diago
          CALL sub_PsiGridRep_TO_BasisRep(psi(i))
        END DO


      ELSE
        DO i=1,nb_diago
          CALL init_psi(psi(i),para_H,cplx)
          CALL alloc_psi(psi(i))
          IF (psi(i)%cplx) THEN
            psi(i)%CvecB(:) = ZERO
            psi(i)%CvecB(i) = ONE
!           DO j=1,psi(i)%nb_baie
!             CALL random_number(a)
!             psi(i)%CvecB(j) = cmplx(a-HALF,ZERO,kind=Rkind)
!           END DO
          ELSE
            psi(i)%RvecB(:) = ZERO
            psi(i)%RvecB(i) = ONE
!           DO j=1,psi(i)%nb_baie
!             CALL random_number(a)
!             psi(i)%RvecB(j) = a-HALF
!           END DO
          END IF

        END DO

      END IF

      DO i=1,nb_diago
        CALL renorm_psi(psi(i))
        CALL sub_PsiOpPsi(CEne0,psi(i),Hpsi,para_H)
        Ene0(i) = CEne0
      END DO

      CALL trie_psi(psi,Ene0,nb_diago)
!------ initialization -------------------------------------


!------- propagation loop ---------------------------------
      CALL file_open(para_propa%file_WP,nioWP)
      it_all = 0
      DO i=1,nb_diago

        T    = ZERO
        it   = 0

!       - Schmidt ortho ------------------------------------
        DO j=1,i-1
          CALL Overlap_psi1_psi2(CS,psi(i),psi(j))
          RS = - real(CS,kind=Rkind)
          psi(i) = psi(i) + psi(j) * RS
        END DO
        CALL renorm_psi(psi(i))
!         - Schmidt ortho ------------------------------------

!       Hpsi = H.psi(i)
        CALL sub_PsiOpPsi(CEne0,psi(i),Hpsi,para_H)
        Ene0(i) = CEne0
        DeltaE = para_H%Hmax-para_H%Hmin
        normeg = para_H%Hmax-para_H%Hmin

        write(out_unitp,*) '-------------------------------------------'
        write(out_unitp,*) 'WP i, E',i,                                 &
                        Ene0(i)*get_Conv_au_TO_unit('E','cm-1'),DeltaE

!       DO WHILE (T .LE. para_propa%WPTmax .AND. abs(DeltaE) .GT. epsi)
        DO WHILE (T <= para_propa%WPTmax .AND. normeg > ONETENTH**5)

!         - Schmidt ortho ------------------------------------
          DO j=1,i-1
            CALL Overlap_psi1_psi2(CS,psi(i),psi(j))
            RS = - real(CS,kind=Rkind)
            psi(i) = psi(i) + psi(j) * RS
          END DO
          CALL renorm_psi(psi(i))
!         - Schmidt ortho ------------------------------------

!         - CG minimization ------------------------------------
!         Hpsi = H.psi
          CALL sub_OpPsi(psi(i),Hpsi,para_H)
          CALL sub_scaledOpPsi(psi(i),Hpsi,para_H%E0,ONE)

!         H2psi = H.H.psi = H.Hpsi
          CALL sub_OpPsi(Hpsi,H2psi,para_H)
          CALL sub_scaledOpPsi(Hpsi,H2psi,para_H%E0,ONE)

!         H3psi = H.H.H.psi = H.H2psi
!         <psi|H3|psi> = <Hspi|H2psi> when psi is real
          CALL Overlap_psi1_psi2(CavH1,psi(i),Hpsi)
          CALL Overlap_psi1_psi2(CavH2,psi(i),H2psi)
          CALL Overlap_psi1_psi2(CavH3,Hpsi,H2psi)
          avH1 = CavH1
          avH2 = CavH2
          avH3 = CavH3


          IF (debug) write(out_unitp,*) 'avH..',avH1,avH2,avH3

!         |g> = H |psi> - E|psi> with E = <psi|H|psi> (E=avH1)
!         Rq : <g|psi> = 0
!         GC :  we optimize |psi+> = cos(th)*|psi> + sin(th)*|g>
!               such the energy, E(th), is minimal
!               E(th) = cos2(th)*A + sin2(th)*B + 2cos(th)sin(th)*C
!               => tan(2th) = -2C/(B-A)
          A = avH1
          B = avH3-TWO*avH2*avH1+avH1**3
          C = avH2 - avH1**2
          B = B/C
          C = sqrt(C)

          write(out_unitp,*) 'A,B,C',A,B,C

          th = HALF*atan( -TWO*C/(B-A) )
          E1 = A*cos(th)**2 + B*sin(th)**2 + TWO*C*cos(th)*sin(th)


          write(out_unitp,*) 'it,th E',it,th,E1*get_Conv_au_TO_unit('E','cm-1')


          Ene0(i) = E1
          DeltaE = Ene0(i) - avH1

!         - Optimal DeltaT and energy ----------------------------

          IF (i == psi(i)%nb_baie) DeltaE = ZERO

!         - propagation ------------------------------------------
          g = Hpsi + psi(i) * (-avH1)
          CALL norm2_psi(g)
          normeg = sqrt(g%norme)
          CALL renorm_psi_WITH_norm2(g)
          write(out_unitp,*) 'it normeg C',it,normeg,C


          H2psi  = psi(i) * cos(th)
          psi(i) = H2psi + g * sin(th)


          CALL norm2_psi(psi(i))
          write(out_unitp,*) 'it norme npsi',it,psi(i)%norme
          CALL renorm_psi_WITH_norm2(psi(i))
          CALL sub_PsiOpPsi(CEne0,psi(i),Hpsi,para_H)
          Ene0(i) = CEne0
          write(out_unitp,*) 'it E with npsi',it,                       &
                                Ene0(i)*get_Conv_au_TO_unit('E','cm-1')


!         - propagation ------------------------------------------

          IF (debug) THEN
            write(out_unitp,*) 'WP it, E',it,                           &
                       Ene0(i)*get_Conv_au_TO_unit('E','cm-1'),         &
                        DeltaE*get_Conv_au_TO_unit('E','cm-1')
            write(out_unitp,*) 'WP it sqrt(norme g)',it,                &
                         sqrt(g%norme)*get_Conv_au_TO_unit('E','cm-1')
            write(out_unitp,*) 'WP it, DeltaT,S',it,DeltaT,S
          END IF

          T = T + DeltaT
          it = it + 1

        END DO
        IF (T> para_propa%WPTmax) write(out_unitp,*) ' WARNING: not converged!'
        write(out_unitp,*) 'WP i, E',i,                                 &
                         Ene0(i)*get_Conv_au_TO_unit('E','cm-1'),DeltaE
        write(out_unitp,*) 'it,T,Deltapsi',it,T,Deltapsi
        it_all = it_all + it
      END DO
      write(out_unitp,*) '-------------------------------------------'
      write(out_unitp,*) 'iterations for eigenvalues: ',it_all
      write(out_unitp,*) '-------------------------------------------'
!----------------------------------------------------------



!----------------------------------------------------------
!     - write the final WP---------------------------------
!     DO i=1,nb_diago
!       write(out_unitp,*) 'WP (i) (BasisRep) ',i
!       CALL ecri_psi(T,psi(i))
!     END DO
!----------------------------------------------------------

      close(nioWP)

      CALL dealloc_psi(Hpsi)
      CALL dealloc_psi(H2psi)
      CALL dealloc_psi(g)

!----------------------------------------------------------
      write(out_unitp,21) T,DeltaE*get_Conv_au_TO_unit('E','cm-1'),     &
                     real(Ene0(:))*get_Conv_au_TO_unit('E','cm-1')
 21   format('ImTimeProp ',f12.2,f18.6,1x,20(1x,f18.4))
      write(out_unitp,*)
!----------------------------------------------------------


!----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'END sub_propagation34'
       END IF
!----------------------------------------------------------


      end subroutine sub_propagation34
!================================================================
!
!     11 : Cheby or nOD propagation
!          without H  (we use directly OpPsi)
!          Hmin and Hmax are the parameter to scale H
!          psi is the intial WP
!          w1,w2,w3 are working WP
!
!================================================================
      SUBROUTINE sub_propagation11(psi0,psi,nb_WP,para_H,para_propa)
      USE mod_system
      USE mod_Op
      USE mod_psi_set_alloc
      USE mod_psi_B_TO_G
      USE mod_ana_psi
      USE mod_propa
      USE mod_march
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ----------------------------
      integer                      :: nb_WP
      TYPE (param_propa)           :: para_propa
      TYPE (param_psi), intent(in) :: psi0(nb_WP)
      TYPE (param_psi)             :: psi(nb_WP)

!------ working parameters --------------------------------
      complex (kind=Rkind) :: cdot

      integer       :: it,itmax,no,no_restart
      integer       :: i,j,i_bie,i_qaie_corr
      integer       :: max_ecri,Write_restart
      character (len=Name_len) :: name_dum
      real (kind=Rkind) :: T      ! time

      integer  ::   nioWP
      TYPE (param_psi)             :: wp_Adia
      logical :: BasisRep,GridRep,test
      character (len=len_trim(para_propa%file_WP_restart%name)) :: Restart_file_name

!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (nb_WP > 1) STOP
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING sub_propagation11'
        write(out_unitp,*) 'Tmax,DeltaT',para_propa%WPTmax,para_propa%WPdeltaT
        write(out_unitp,*) 'Hmin,Hmax',para_propa%para_poly%Hmin,               &
                                para_propa%para_poly%Hmax
        write(out_unitp,*)
        write(out_unitp,*) 'nb_ba,nb_qa',psi(1)%nb_ba,psi(1)%nb_qa
        write(out_unitp,*) 'nb_bi',psi(1)%nb_bi
        write(out_unitp,*)

        write(out_unitp,*) 'psiBasisRep'
        CALL ecri_psi(T=ZERO,psi=psi(1),ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)
        write(out_unitp,*) 'psiGridRep'
        CALL ecri_psi(T=ZERO,psi=psi(1),ecri_GridRep=.TRUE.,ecri_BasisRep=.FALSE.)
      END IF
!-----------------------------------------------------------


!-----------------------------------------------------------
      BasisRep = psi0(1)%BasisRep
      GridRep  = psi0(1)%GridRep

      write(out_unitp,*) ' vib : propagation',para_propa%name_WPpropa

!     - parameters for poly (cheby and nOD) ... ------------
      CALL initialisation1_poly(para_propa%para_poly,                   &
                                para_propa%WPdeltaT,                    &
                                para_propa%type_WPpropa)

!     - scaling of H ---------------------------------------
      para_H%scaled = .TRUE.
      para_H%E0     = para_propa%para_poly%E0
      para_H%Esc    = para_propa%para_poly%Esc
!-----------------------------------------------------------


!-----------------------------------------------------------

!------- propagation loop ---------------------------------

      T = ZERO
      IF (para_propa%restart) THEN
        CALL file_open(para_propa%file_autocorr,no,append=.TRUE.)
        backspace(no)
        read(no,*) name_dum,T
        close(no)
        write(out_unitp,*) 'T0 for the restart:',T

        CALL file_open(para_propa%file_WP_restart,no_restart)
        read(no_restart,*) psi(1)%CvecB
        close(no_restart)

        CALL file_open(para_propa%file_autocorr,no,append=.TRUE.)

      ELSE
        CALL file_open(para_propa%file_autocorr,no)
        cdot = Calc_AutoCorr(psi0(1),psi(1),para_propa,T,Write_AC=.TRUE.)
      END IF

      Restart_file_name = para_propa%file_WP_restart%name
      it            = 0
      itmax         = (para_propa%WPTmax-T)/para_propa%WPdeltaT

      DO WHILE ( (T - (para_propa%WPTmax-para_propa%WPdeltaT) <         &
                 para_propa%WPdeltaT/TEN**5) .AND. psi(1)%norme < psi(1)%max_norme)

           para_propa%ana_psi%Write_psi2_Grid = (mod(it,para_propa%n_WPecri) == 0) .AND. para_propa%WPpsi2
           para_propa%ana_psi%Write_psi_Grid  = (mod(it,para_propa%n_WPecri) == 0) .AND. para_propa%WPpsi

           CALL sub_analyze_WP_OpWP(T,psi,1,para_H,para_propa)

           !para_propa%file_WP_restart%name = Restart_file_name // '_it' // int_TO_char(mod(it,max(1,itmax/100)))
           !CALL file_open(para_propa%file_WP_restart,no_restart)
           !write(no_restart,*) psi(1)%CvecB
           !close(no_restart)

         CALL march_gene(T,psi(1:1),psi0(1:1),1,.FALSE.,para_H,para_propa)

         it = it + 1
         T  = T + para_propa%WPdeltaT

      END DO

!----------------------------------------------------------


!----------------------------------------------------------
!     - write the final WP---------------------------------
      write(out_unitp,*) 'WP (BasisRep) at T=',T
      IF (debug) CALL ecri_psi(T=T,psi=psi(1),ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)

      para_propa%ana_psi%Write_psi2_Grid = para_propa%WPpsi2
      para_propa%ana_psi%Write_psi_Grid  = para_propa%WPpsi
      CALL sub_analyze_WP_OpWP(T,psi,1,para_H,para_propa)

      para_propa%file_WP_restart%name = Restart_file_name
      write(out_unitp,*) ' Last restart file: ',trim(para_propa%file_WP_restart%name)
      CALL file_open(para_propa%file_WP_restart,no_restart)
      write(no_restart,*) psi(1)%CvecB
      close(no_restart)
!----------------------------------------------------------

      CALL file_close(para_propa%file_autocorr)
      CALL dealloc_psi(WP_adia,delete_all=.TRUE.)
      IF (psi(1)%norme >= psi(1)%max_norme) STOP


!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END sub_propagation11'
      END IF
!----------------------------------------------------------

      end subroutine sub_propagation11
!================================================================
!
!    24 : nOD propagation with
!         a time dependant pulse in Hamiltonian (W(t))
!         H is the square matrix (dimension n)
!         Hmin and Hmax are the parameter to scale H
!
!================================================================
      SUBROUTINE sub_propagation24(WP,nb_WP,print_Op,                   &
                                   para_field_new,make_field,           &
                                   para_H,para_Dip,para_propa)
      USE mod_system
      USE mod_psi_B_TO_G
      USE mod_Op
      USE mod_propa
      USE mod_march
      !USE mod_psi
      USE mod_psi_set_alloc
      USE mod_ana_psi
      USE mod_field
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H
      TYPE (param_Op)   :: para_Dip(3)

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_field) :: para_field_new
      real (kind=Rkind)      :: E_new

      integer            :: nb_WP
      TYPE (param_psi)   :: WP(nb_WP)

!----- for printing --------------------------------------------------
      logical ::print_Op
      logical ::print_Op_loc


!------ working parameters --------------------------------
      complex (kind=Rkind) :: E0(nb_WP),avE(nb_WP)  ! energy

      integer       :: it,it_max,i
      real (kind=Rkind) :: T      ! time

!----- for the field --------------------------------------------------
      real (kind=Rkind)    :: ww(3)
      real (kind=Rkind)    :: dnE(3)
      logical :: make_field
!----- for the field --------------------------------------------------


!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING sub_propagation24'
        write(out_unitp,*) 'Tmax,deltaT',para_propa%WPTmax,para_propa%WPdeltaT
        write(out_unitp,*) 'Hmin,Hmax',para_propa%para_poly%Hmin,               &
                                para_propa%para_poly%Hmax
        write(out_unitp,*)
        write(out_unitp,*) 'nb_ba,nb_qa',WP(1)%nb_ba,WP(1)%nb_qa
        write(out_unitp,*) 'nb_bi',WP(1)%nb_bi
        write(out_unitp,*)

        write(out_unitp,*) 'WP(1)%BasisRep'
        CALL ecri_psi(T=ZERO,psi=WP(1),ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)
        write(out_unitp,*) 'WP(1)%GridRep'
        CALL ecri_psi(T=ZERO,psi=WP(1),ecri_GridRep=.TRUE.,ecri_BasisRep=.FALSE.)
       END IF
!-----------------------------------------------------------

!-----------------------------------------------------------
      IF (print_Op) THEN
        write(out_unitp,*) ' Propagation ',para_propa%name_WPpropa
        IF (para_propa%para_field%pola_xyz(1)) write(out_unitp,*) 'with Dipx'
        IF (para_propa%para_field%pola_xyz(2)) write(out_unitp,*) 'with Dipy'
        IF (para_propa%para_field%pola_xyz(3)) write(out_unitp,*) 'with Dipz'
      END IF

!     - parameters for poly (cheby and nOD) ... ------------
      CALL initialisation1_poly(para_propa%para_poly,                   &
                                para_propa%WPdeltaT,                    &
                                para_propa%type_WPpropa)

!-----------------------------------------------------------

!     - scaling of H ---------------------------------------
!     para_H%scaled = .FALSE.
!     para_H%E0     = ZERO
!     para_H%Esc    = ONE

      para_H%scaled = .TRUE.
      para_H%E0     = para_propa%para_poly%E0
      para_H%Esc    = para_propa%para_poly%Esc

      para_propa%para_poly%phase = para_H%E0*para_propa%WPdeltaT
      write(out_unitp,*) 'para_H%E0,para_H%Esc',para_H%E0,para_H%Esc
      write(out_unitp,*) 'phase',para_propa%para_poly%phase

!-----------------------------------------------------------
      T = ZERO
      IF (para_propa%WPdeltaT < 0) T = para_propa%WPTmax
      it     = 0
      it_max = para_propa%WPTmax/abs(para_propa%WPdeltaT)-1
      avE(:) = ZERO

      ww = para_propa%para_field%w(:,1)
      IF (print_Op) THEN
        write(out_unitp,*) 'For the first pulse'
        write(out_unitp,*) '  nb_pulse',para_propa%para_field%nb_pulse
        write(out_unitp,*) 'ww (ua),FieldE0 (ua)',ww,                   &
                                           para_propa%para_field%E0(:,1)
        write(out_unitp,*) 'ww (cm-1),FieldE0 (V cm-1)',                &
                                  ww*get_Conv_au_TO_unit('E','cm-1'),   &
            para_propa%para_field%E0(:,1)*get_Conv_au_TO_unit('Electric Field','V.cm-1')
        write(out_unitp,11) para_propa%para_field%E0(:,1)**2 *          &
                           get_Conv_au_TO_unit('EF intensity','W.cm-2')
 11     format('G (W cm-2):  ',d16.8)
        write(out_unitp,*)

      END IF
!-----------------------------------------------------------


!------- propagation loop ---------------------------------
      DO WHILE (it <= it_max)

        para_propa%ana_psi%Write_psi2_Grid = (mod(it,para_propa%n_WPecri) == 0) .AND. para_propa%WPpsi2
        para_propa%ana_psi%Write_psi_Grid  = (mod(it,para_propa%n_WPecri) == 0) .AND. para_propa%WPpsi

        CALL sub_analyze_WP_OpWP(T,WP,nb_WP,para_H,para_propa,          &
                                       para_field=para_propa%para_field)

        IF (it ==0) E0(:)  = WP(:)%CAvOp
        avE(:) = avE(:) + WP(:)%CAvOp-E0(:)


         print_Op_loc = print_Op .AND. mod(it,para_propa%n_WPecri) == 0
         CALL march_gene(T,WP(:),WP(:),nb_WP,print_Op_loc,              &
                         para_H,para_propa,                             &
                         para_Dip,para_propa%para_field)


        T = T + para_propa%WPdeltaT
        it = it + 1

       END DO ! loop on the Time iteration
       CALL flush_perso(out_unitp)

!----------------------------------------------------------
!     - write the final WP---------------------------------
      T  = para_propa%WPTmax
      write(out_unitp,*) '=================================='
      write(out_unitp,*) 'Final normalized WP'
      DO i=1,nb_WP
        CALL renorm_psi(WP(i),GridRep=.FALSE.,BasisRep=.TRUE.)
      END DO

      para_propa%ana_psi%Write_psi2_Grid = para_propa%WPpsi2
      para_propa%ana_psi%Write_psi_Grid  = para_propa%WPpsi
      CALL sub_analyze_WP_OpWP(T,WP,nb_WP,para_H,para_propa,            &
                                 para_field=para_propa%para_field)
      avE(:) = avE(:) + WP(:)%CAvOp-E0(:)

      DO i=1,nb_WP
        avE(i) = avE(i)/real(it_max,kind=Rkind)
        write(out_unitp,31) i,ww,real(avE(i),kind=Rkind)
 31     format('average absorbed energy at w ',i3,4(1x,f20.10))
      END DO ! j loop (nb_WP)

!----------------------------------------------------------


!----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'END sub_propagation24'
       END IF
!----------------------------------------------------------


      end subroutine sub_propagation24
!================================================================
!
!    1OO : subroutine to test things ...
!
!================================================================
      SUBROUTINE sub_propagation100(WP,nb_WP,print_Op,                  &
                                   para_field_new,make_field,           &
                                   para_AllOp,para_propa)
      USE mod_system
      USE mod_Op
      USE mod_propa
      USE mod_march
      !USE mod_psi
      USE mod_psi_set_alloc
      USE mod_psi_io
      USE mod_ana_psi
      USE mod_field
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_AllOp)  :: para_AllOp
      integer             :: iOp
      logical             :: print_Op

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_field) :: para_field_new
      logical :: make_field

      integer            :: nb_WP
      TYPE (param_psi)   :: WP(nb_WP)
      TYPE (param_psi)   :: OpWP(nb_WP)


!----- for printing --------------------------------------------------


!------ working parameters --------------------------------
      TYPE (param_psi)   :: w1,w2,w3

      integer              :: i
      logical              :: cplx
      complex (kind=Rkind) :: ET  ! energy

!----- for the field --------------------------------------------------


!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING sub_propagation100'
        write(out_unitp,*) 'Tmax,deltaT',para_propa%WPTmax,para_propa%WPdeltaT
        write(out_unitp,*) 'Hmin,Hmax',para_propa%para_poly%Hmin,               &
                                para_propa%para_poly%Hmax
        write(out_unitp,*)
        write(out_unitp,*) 'nb_ba,nb_qa',WP(1)%nb_ba,WP(1)%nb_qa
        write(out_unitp,*) 'nb_bi',WP(1)%nb_bi
        write(out_unitp,*)

        write(out_unitp,*) 'WP(1)'
        CALL ecri_psi(T=ZERO,psi=WP(1))
       END IF
!-----------------------------------------------------------
      cplx = .FALSE.
      CALL init_psi(w1,para_AllOp%tab_Op(1),cplx)
      CALL alloc_psi(w1)
      w1%RvecB(:) = (/ (i,i=1,w1%nb_tot) /)
      !CALL sub_PsiBasisRep_TO_GridRep(w1)
      w2 = w1

      write(out_unitp,*) 'debut'
      DO i=1,1
         !write(out_unitp,*) 'i',i ; CALL flush_perso(out_unitp)
         CALL sub_OpPsi(w1,w2,para_AllOp%tab_Op(1))
         !CALL sub_d0d1d2PsiBasisRep_TO_GridRep(w1,(/ 1,2 /))
         !CALL sub_PsiBasisRep_TO_GridRep(w1)
         !CALL sub_PsiGridRep_TO_BasisRep(w1)
      END DO
      write(out_unitp,*) 'end'



      RETURN





      !=================================================
      DO iOp=1,size(para_AllOp%tab_Op)
        IF (para_propa%num_Op == para_AllOp%tab_Op(iOp)%n_Op) EXIT
      END DO
      IF (iOp > size(para_AllOp%tab_Op)) THEN
        write(out_unitp,*) ' ERROR in sub_propagation100'
        write(out_unitp,*) ' The Operator "',para_propa%num_Op,'" is not in the list!'
        write(out_unitp,*) ' Change the "num_Op" value in "propagation namelist"'
        write(out_unitp,*) ' Possible values are:'
        DO iOp=1,size(para_AllOp%tab_Op)
          write(out_unitp,*) 'iOp,num_Op,Name_Op',iOp,para_AllOp%tab_Op(iOp)%n_Op, &
                      para_AllOp%tab_Op(iOp)%name_Op
        END DO
        STOP
      END IF
      write(out_unitp,*) 'num_Op',iOp,para_propa%num_Op

      DO i=1,nb_WP
         write(out_unitp,*) 'WP,i',i
         !CALL ecri_psiBasisRepnotall_nD(WP(i),out_unitp,ONETENTH**10,.TRUE.,i)
         !CALL ecri_psi(Psi=WP(i))
         CALL sub_OpPsi(WP(i),OpWP(i),para_AllOp%tab_Op(iOp))
         CALL Overlap_psi1_psi2(ET,WP(i),OpWP(i))

         write(out_unitp,*) 'OpWP,i',i,ET*get_Conv_au_TO_unit('E','cm-1')

         !CALL ecri_psi(Psi=OpWP(i))
         WP(i) = OpWP(i)
         CALL dealloc_psi(OpWP(i))
      END DO

      write(out_unitp,*) 'file_WP0: ',para_propa%para_WP0%file_WP0%name

      CALL sub_save_psi(WP,nb_WP,para_propa%para_WP0%file_WP0)

      RETURN
      !=================================================


      CALL dealloc_psi(w1)
      CALL dealloc_psi(w2)

!----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'END sub_propagation100'
       END IF
!----------------------------------------------------------

      end subroutine sub_propagation100
END MODULE mod_FullPropa
