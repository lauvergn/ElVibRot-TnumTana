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

!==============================================================
!     Wave packet initialisation
!
!     reading WP0:
!       lect_WP0GridRep = .TRUE.
!   or  lect_WP0BasisRep = .TRUE.
!
!     calculating WP0 on the grid (GridRep)
!          WP0(qi) =  exp[-((Q-Qeq)/sigma)2]*exp[i*imp_k*(Q-Qeq)]
!
!
!     If WP0BasisRep=.TRUE. => WP0 is BasisRep and nWP0=dim(WP0BasisRep)
!               .FALSE.=> WP0 is GridRep and nWP0=dim(WP0GridRep)
!
!     If WP0_DIP=1,2,3     => WP0 = WP0*dip(i) [i=1,2,3 => x,y,z]
!             =0 (default) => WP0 = WP0
!==============================================================
      SUBROUTINE psi0(WP0,para_WP0,mole)
      USE mod_system
      USE mod_Coord_KEO
      USE mod_basis
      USE mod_psi_set_alloc
      USE mod_psi_Op
      USE mod_psi_B_TO_G
      USE mod_param_WP0
      USE mod_ana_psi
      USE mod_psi_io
#IF(run_MPI)      
      USE mod_MPI
      USE mod_MPI_Aid
#ENDIF      
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_WP0) :: para_WP0
      TYPE (param_psi), intent(inout) :: WP0(1)

!----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

!------ working variables ---------------------------------
      integer  :: nio,nb_WPdum
      integer  :: ecri_numi,ecri_nume

!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING psi0'
        write(out_unitp,*) 'WP0n_h,WP0nb_elec',para_WP0%WP0n_h,para_WP0%WP0nb_elec
        write(out_unitp,*)
        CALL ecri_init_psi(WP0(1))

        write(out_unitp,*)
        write(out_unitp,*) 'WP0BasisRep',para_WP0%WP0BasisRep
        write(out_unitp,*) 'lect_WP0GridRep,lect_WP0BasisRep',         &
                  para_WP0%lect_WP0GridRep,para_WP0%lect_WP0BasisRep
        write(out_unitp,*)
        write(out_unitp,*) 'nb_basis_act1',WP0(1)%BasisnD%nb_basis
        CALL RecWrite_basis(WP0(1)%BasisnD)
        write(out_unitp,*)
        CALL flush_perso(out_unitp)
      END IF

!-----------------------------------------------------------
      IF (para_WP0%New_Read_WP0) THEN
        IF(MPI_id==0) CALL sub_read_psi0(WP0,para_WP0,max_WP=1)
      ELSE
        ecri_numi = para_WP0%WP0n_h
        ecri_nume = para_WP0%WP0nb_elec


        ! alloc the WP0 ----------------------------------------------------------------
        ! only allocated for MPI_id=0 inside the subroutine
        ! for psi%CvecB,psi%RvecB,psi%CvecG,psi%RvecG      
        CALL alloc_psi(WP0(1))

        IF (.NOT. para_WP0%lect_WP0BasisRep .AND.                         &
            .NOT. para_WP0%WP0restart) THEN

          CALL alloc_psi(WP0(1),GridRep=.TRUE.)

          IF (para_WP0%lect_WP0GridRep) THEN
            !- read WP0 on the grid ----------------------------
            CALL lect_psiBasisRep(WP0(1)%CvecG,para_WP0%WP0cplx,          &
                             WP0(1)%nb_qa,WP0(1)%nb_bi,WP0(1)%nb_be,      &
                             para_WP0%WP0n_h,para_WP0%WP0nb_elec)
          ELSE
            !- initialisation with a gaussian WP ---------------
            CALL psi0_gaussGridRep(WP0(1),para_WP0,mole)
          END IF

          IF (debug) THEN
            write(out_unitp,*) 'psiGridRep ini'
            CALL ecri_psi(ZERO,WP0(1),out_unitp,.TRUE.,.FALSE.)
          END IF

          CALL norm2_psi(WP0(1),GridRep=.TRUE.)
          IF(MPI_id==0) write(out_unitp,*) 'normeWP GridRep',WP0(1)%norme
          CALL flush_perso(out_unitp)
          CALL renorm_psi_WITH_norm2(WP0(1),GridRep=.TRUE.)
          IF(MPI_id==0) write(out_unitp,*) 'normeWP GridRep',WP0(1)%norme
          CALL flush_perso(out_unitp)

          IF (debug) THEN
            write(out_unitp,*) 'psiGridRep normalized'
            CALL ecri_psi(ZERO,WP0(1),out_unitp,.TRUE.,.FALSE.)
          END IF

          !- GridRep=>BasisRep -------------------------------------------------
          IF (para_WP0%WP0BasisRep) THEN

            CALL sub_PsiGridRep_TO_BasisRep(WP0(1))

            CALL norm2_psi(WP0(1),BasisRep=.TRUE.)

            write(out_unitp,*) 'normeWP BasisRep',WP0(1)%norme

            IF (abs(ONE-WP0(1)%norme) >= ONETENTH**5) THEN
              write(out_unitp,*) ' WARNNIG in psi0'
              write(out_unitp,*) ' the transformation GridRep to BasisRep is NOT exact'
              write(out_unitp,*) ' => used more basis functions'
              !CALL ecri_psi(psi=WP0(1),                                      &
              !              ecri_BasisRep=.TRUE.,ecri_GridRep=.TRUE.)
            END IF
            CALL renorm_psi(WP0(1),BasisRep=.TRUE.)

          END IF

        ELSE IF (para_WP0%WP0restart ) THEN
          !write(out_unitp,*) 'WP0%cplx',WP0%cplx
          CALL file_open(para_WP0%file_WP0,nio)
          read(nio,*) nb_WPdum
          CALL lect_psiBasisRepnotall_nD(WP0(1),nio,WP0(1)%cplx,para_WP0%file_WP0%formatted)
          close(nio)

        ELSE IF (para_WP0%lect_WP0BasisRep .AND. para_WP0%lect_WP0BasisRepall) THEN

          IF (WP0(1)%nb_baie > WP0(1)%nb_tot) THEN
          CALL lect_psiBasisRep(WP0(1)%CvecB,para_WP0%WP0cplx,                  &
                           WP0(1)%nb_tot,1,1,                                   &
                           para_WP0%WP0n_h,para_WP0%WP0nb_elec)
          ELSE
          CALL lect_psiBasisRep(WP0(1)%CvecB,para_WP0%WP0cplx,                  &
                           WP0(1)%nb_ba,WP0(1)%nb_bi,WP0(1)%nb_be,                 &
                           para_WP0%WP0n_h,para_WP0%WP0nb_elec)
          END IF

        ELSE IF (para_WP0%lect_WP0BasisRep                                     &
                                .AND. .NOT. para_WP0%lect_WP0BasisRepall) THEN

          CALL lect_psiBasisRepnotall(WP0(1),para_WP0%WP0cplx)
        ELSE
          write(out_unitp,*) ' ERROR in psi0'
          write(out_unitp,*) ' I do not what to do!!!'
          STOP
        END IF ! for .NOT. para_WP0%lect_WP0BasisRep

!     IF (para_WP0%WP0_DIP .GT. 0) THEN
!       IF (para_WP0%WP0BasisRep) THEN
!       ELSE
!       ENDIF
!     END IF

      END IF ! for para_WP0%New_Read_WP0

      IF(MPI_id==0) CALL renorm_psi(WP0(1),BasisRep=.TRUE.)
      IF(MPI_id==0) write(out_unitp,*) 'normeWP BasisRep',WP0(1)%norme

      !- clear WP0%...GridRep, if not need ------------------
      IF (para_WP0%WP0BasisRep) THEN
        WP0(1)%GridRep = .FALSE.
        CALL alloc_psi(WP0(1)) ! deallocate the grid representation
      END IF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'WP0BasisRep',WP0(1)%norme
        IF (para_WP0%WP0BasisRep) THEN
          CALL ecri_psi(ZERO,WP0(1),out_unitp,.FALSE.,.TRUE.)
        ELSE
          CALL ecri_psi(ZERO,WP0(1),out_unitp,.TRUE.,.FALSE.)
        END IF
        write(out_unitp,*) 'END psi0'
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

      END SUBROUTINE psi0
!==============================================================
!     Wave packet initialisation
!
!     for each active dimension :
!          psi0(qi) =  exp[-((Q-Qeq)/sigma)^2]*exp[i*imp_k*(Q-Qeq)+i*phase]
!
!     If WP0_DIP=1,2,3     => WP0 = WP0*dip(i) [i=1,2,3 => x,y,z]
!             =0 (default) => WP0 = WP0
!
!==============================================================
      SUBROUTINE psi0_gaussGridRep(WP0,para_WP0,mole)
      USE mod_system
      USE mod_Coord_KEO
      USE mod_psi_set_alloc
      USE mod_param_WP0
      IMPLICIT NONE


!----- variables for the WP propagation ----------------------------
      TYPE (param_WP0) :: para_WP0
      TYPE (param_psi)   :: WP0
      integer            :: WP0nb_elec,WP0n_h

!----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole


      real (kind=Rkind)      :: Qact(WP0%nb_act1)



      real (kind=Rkind)      :: z,ze,zk
      real (kind=Rkind)      :: czk,szk
      integer                :: i_act
      integer                :: i_qa,i_qaie

!------ dipole moment ----------------------------------------------------
      real (kind=Rkind)      :: rhonD

!----- for debuging --------------------------------------------------
      logical,parameter :: debug = .FALSE.
      !logical,parameter :: debug = .TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING psi0_gaussGridRep'
        write(out_unitp,*) 'nb_act1',WP0%nb_act1
        write(out_unitp,*) 'nq_a,n_h,nb_elec',WP0%nb_qa,WP0%nb_bi,WP0%nb_be
        write(out_unitp,*) 'WP0n_h,WP0nb_elec',para_WP0%WP0n_h,para_WP0%WP0nb_elec
        write(out_unitp,*) 'sigma,Qeq,imp_k',                                   &
           para_WP0%WP0sigma,para_WP0%WP0Qeq,para_WP0%WP0imp_k
        write(out_unitp,*) 'WP0_DIP',para_WP0%WP0_DIP
        write(out_unitp,*)
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------


      WP0n_h = para_WP0%WP0n_h
      WP0nb_elec = para_WP0%WP0nb_elec
      IF (WP0n_h <1) WP0n_h = 1
      IF (WP0nb_elec <1) WP0nb_elec = 1



!     - initisalisation ----------------------------------
      IF (WP0%cplx) THEN
        WP0%CvecG(:) = cmplx(ZERO,ZERO,kind=Rkind)
      ELSE
        WP0%RvecG(:) = ZERO
      END IF

      DO i_qa=1,WP0%nb_qa

!       - calculation of Qact -------------------------------
        CALL Rec_Qact(Qact,WP0%BasisnD,i_qa,mole)
!       - calculation of WrhonD ------------------------------
        rhonD = Rec_rhonD(WP0%BasisnD,i_qa)

!       ------------------------------------------------
        ze = ZERO
        zk = ZERO
        DO i_act=1,WP0%nb_act1
           z = (Qact(i_act)-para_WP0%WP0Qeq(i_act))/                    &
                                     para_WP0%WP0sigma(i_act)
           ze = ze + z*z
           zk = zk + para_WP0%WP0phase(i_act) +                         &
             (Qact(i_act)-para_WP0%WP0Qeq(i_act))*para_WP0%WP0imp_k(i_act)
        END DO
        zk = mod(zk,TWO*pi)
        ze = exp(-ze)/sqrt(sqrt(pi/TWO))**WP0%nb_act1
        ze = ze /sqrt(product(para_WP0%WP0sigma))
        IF (para_WP0%WP0nrho == 1) ze = ze /sqrt(rhonD)
        czk = cos(zk)
        szk = sin(zk)

        i_qaie = i_qa + ( (WP0n_h-1) +                                  &
               (WP0nb_elec-1)*WP0%nb_bi ) * WP0%nb_qa

        IF (WP0%cplx) THEN
          WP0%CvecG(i_qaie) = cmplx(czk*ze,szk*ze,kind=Rkind)
          IF (debug) write(out_unitp,11) Qact,WP0%CvecG(i_qaie)
 11       format(10f15.5)
        ELSE
          WP0%RvecG(i_qaie) = ze
          IF (debug) write(out_unitp,11) Qact,WP0%RvecG(i_qaie)
        END IF

      END DO

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'END psi0_gaussGridRep'
         CALL flush_perso(out_unitp)
       END IF
!-----------------------------------------------------------


      END SUBROUTINE psi0_gaussGridRep
