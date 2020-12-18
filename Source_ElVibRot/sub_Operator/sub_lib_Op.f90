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

!=====================================================================
!
!  Op(Q).d0b(Q,i).W(Q) and d0b(Q,i) calculations for the nD quadrature point k
!
!=====================================================================
      SUBROUTINE calc_td0b_OpRVd0bW(iq,k,td0b,d0MatOpd0bWrho,WnD,kmem,  &
                                    d0MatOp,para_Op,BasisnD)
      USE mod_system
      USE mod_PrimOp, only: Write_d0MatOp
      USE mod_basis
      USE mod_SetOp
      IMPLICIT NONE

!----- for the basis set ----------------------------------------------
      TYPE (Basis) :: BasisnD

!----- Operator variables --------------------------------------------
      TYPE (param_Op)      :: para_Op
      integer              :: k,kmem,iq
      real (kind=Rkind)    :: td0b(para_Op%nb_ba,kmem)
      TYPE (param_d0MatOp) :: d0MatOpd0bWrho(kmem,para_Op%nb_ba)


      TYPE (param_d0MatOp) :: d0MatOp

!------ quadrature weights ----------------------------------------
      real (kind=Rkind) :: WnD


!------ working variables ---------------------------------
      integer  :: ib
      integer  :: j_act,i_act,i,j

      real (kind=Rkind) :: d0bnD
      real (kind=Rkind) :: d1bnD(para_Op%nb_act1)
      real (kind=Rkind) :: d2bnD(para_Op%nb_act1,para_Op%nb_act1)

      real (kind=Rkind) :: dnbnD

      integer ::  i_term,i_Op

!----- for debuging --------------------------------------------------
       character (len=*), parameter :: name_sub = 'calc_td0b_OpRVd0bW'
       logical, parameter :: debug = .FALSE.
       !logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'iq,k,kmem',iq,k,kmem
         write(out_unitp,*) 'nb_ba',para_Op%nb_ba
         write(out_unitp,*) 'nb_act1',para_Op%nb_act1
         write(out_unitp,*)
         write(out_unitp,*) 'WnD',WnD
         ! CALL write_param_Op(para_Op)
          write(out_unitp,*) 'd0MatOp:'
          CALL Write_d0MatOp(d0MatOp)
         write(out_unitp,*)
       END IF
!-----------------------------------------------------------


!-----------------------------------------------------------
!     test if the basis is real
!     ------------------------------------------------------
      IF (BasisnD%cplx) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' You are using the REAL subroutine '
         write(out_unitp,*) ' for the d0b calculation, but the basis is COMPLEX'
         STOP
      END IF
!-----------------------------------------------------------

!-----------------------------------------------------------
      IF (d0MatOp%cplx) THEN
        d0MatOp%ImVal(:,:)    = d0MatOp%ImVal(:,:) * WnD
      END IF
      d0MatOp%ReVal(:,:,:)    = d0MatOp%ReVal(:,:,:) * WnD

      DO ib=1,para_Op%nb_ba

        CALL d0d1d2bnDQact(d0bnD,d1bnD,d2bnD,BasisnD,iq,ib,para_Op%mole)
        td0b(ib,k)            = d0bnD

        !- initialisation ----
        d0MatOpd0bWrho(k,ib)%ReVal(:,:,:) = ZERO

        DO i_term=1,d0MatOp%nb_term
          i_act = d0MatOp%derive_termQact(1,i_term)
          j_act = d0MatOp%derive_termQact(2,i_term)

          !- 2d order derivatives ------------------
          IF ( j_act > 0 .AND. i_act > 0 ) THEN
            dnbnD = d2bnD(j_act,i_act)

          !- first order derivatives ---------------
          ELSE IF ( j_act > 0 .AND. i_act <= 0 ) THEN
            dnbnD = d1bnD(j_act)

          ELSE IF ( j_act <= 0 .AND. i_act > 0 ) THEN
            dnbnD = d1bnD(i_act)

          !- no derivative  of deformation part ----------------------
          ELSE
            dnbnD = d0bnD

          END IF

          i    = min(0,i_act)
          j    = min(0,j_act)
          i_Op = d0MatOpd0bWrho(k,ib)%derive_term_TO_iterm(j,i)

          d0MatOpd0bWrho(k,ib)%ReVal(:,:,i_Op) =                        &
                                d0MatOpd0bWrho(k,ib)%ReVal(:,:,i_Op) +  &
                                       dnbnD * d0MatOp%ReVal(:,:,i_term)
        END DO

        IF (d0MatOp%cplx) THEN
          d0MatOpd0bWrho(k,ib)%ImVal(:,:) = d0bnD * d0MatOp%ImVal(:,:)
        END IF

      END DO
      !-----------------------------------------------------------------


!-----------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) ' td0b(:,k)',k,para_Op%nb_ba
        CALL Write_Vec(td0b(:,k),out_unitp,8)
        write(out_unitp,*)
        write(out_unitp,*) ' d0MatOpd0bWrho(:,:)'
        DO i=1,ubound(d0MatOpd0bWrho,dim=2)
          write(out_unitp,*) 'k,i',k,i
          CALL Write_d0MatOp(d0MatOpd0bWrho(k,i))
        END DO
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF

      END SUBROUTINE calc_td0b_OpRVd0bW
