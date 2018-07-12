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
!=====================================================================
!=====================================================================
!
!      hermiticity analysis
!      and symetrisation
!
!=====================================================================
      SUBROUTINE sub_hermitic_H(H,nb_bases,non_hermitic,sym)
      USE mod_system
      IMPLICIT NONE

!------ active Matrix H ------------------------------------------

      integer   nb_bases
      real (kind=Rkind) ::    H(nb_bases,nb_bases)
      real (kind=Rkind) ::    non_hermitic
      logical           ::    sym


!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING sub_hermitic_H'
      END IF
!-----------------------------------------------------------

      non_hermitic = maxval(abs(H-transpose(H)))*HALF
      IF (sym) H = (H+transpose(H))*HALF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'non_hermitique H',non_hermitic
        write(out_unitp,*) 'END sub_hermitic_H'
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_hermitic_H
!=====================================================================
!
!      hermiticity analysis
!      and symetrisation
!
!=====================================================================
      SUBROUTINE sub_hermitic_cplxH(H,nb_bases,non_hermitic,sym)
      USE mod_system
      IMPLICIT NONE

!------ active Matrix H ------------------------------------------

      integer       nb_bases
      complex (kind=Rkind) ::    H(nb_bases,nb_bases)
      real (kind=Rkind)    ::    val,non_hermitic
      logical              ::    sym


      integer   i,j

!----- for debuging --------------------------------------------------
      logical debug
      parameter (debug=.FALSE.)
!     parameter (debug=.TRUE.)
!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING sub_hermitic_cplxH'
      END IF
!-----------------------------------------------------------

      non_hermitic = ZERO
      DO i=1,nb_bases
        DO j=i+1,nb_bases

          val = abs( H(i,j) - H(j,i) )
          IF (sym) THEN
            H(i,j) = (H(i,j) + H(j,i) )*cmplx(HALF,ZERO,kind=Rkind)
            H(j,i) = H(i,j)
          END IF

          IF ( val > non_hermitic) non_hermitic = val
!         write(out_unitp,*) 'sub_hermitic_cplxH',val,non_hermitic
        END DO
      END DO

      non_hermitic = non_hermitic * HALF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'non_hermitique cplxH',non_hermitic
        write(out_unitp,*) 'END sub_hermitic_cplxH'
      END IF
!-----------------------------------------------------------

      RETURN
      end subroutine sub_hermitic_cplxH
!=====================================================================
!
!       analysis of the overlap matrix (upper part)
!
!=====================================================================
      SUBROUTINE sub_ana_S(S,nb_bases,max_Sii,max_Sij,write_maxS)
      USE mod_system
      IMPLICIT NONE

!------ active Matrix H ------------------------------------------

      integer   nb_bases
      real (kind=Rkind) :: S(nb_bases,nb_bases)
      real (kind=Rkind) :: max_Sii,max_Sij
      logical           :: write_maxS


      integer   i,j

!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING sub_ana_S'
        CALL Write_Mat(S,out_unitp,5)
      END IF
!-----------------------------------------------------------

!     - analysis of the overlap matrix
      max_Sii = ZERO
      max_Sij = ZERO
      DO i=1,nb_bases
        max_Sii=max(max_Sii,abs(S(i,i)-ONE))
        DO j=i+1,nb_bases
          max_Sij = max(max_Sij,abs(S(i,j)))
        END DO
      END DO

      IF (write_maxS) THEN
         write(out_unitp,"(' Max Overlap:',2e11.3)") max_Sii,max_Sij
         CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) ' Max Overlap:',max_Sii,max_Sij
        write(out_unitp,*) 'END sub_ana_S'
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_ana_S
!=====================================================================
!
!       analysis of the overlap matrix (upper part)
!
!=====================================================================
      SUBROUTINE sub_ana_cplxS(S,nb_bases,max_Sii,max_Sij,write_maxS)
      USE mod_system
      IMPLICIT NONE

!------ active Matrix H ------------------------------------------

      integer              :: nb_bases
      complex (kind=Rkind) :: S(nb_bases,nb_bases)
      real (kind=Rkind)    :: max_Sii,max_Sij
      logical              :: write_maxS


      integer   :: i,j

!----- for debuging --------------------------------------------------
      logical debug
      parameter (debug=.FALSE.)
!     parameter (debug=.TRUE.)
!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING sub_ana_cplxS'
      END IF
!-----------------------------------------------------------

!     - analysis of the overlap matrix
      max_Sii = ZERO
      max_Sij = ZERO
      DO i=1,nb_bases
        max_Sii=max(max_Sii,abs(S(i,i)-CONE))
        DO j=i+1,nb_bases
          max_Sij = max(max_Sij,abs(S(i,j)))
        END DO
      END DO

      IF (write_maxS) THEN
         write(out_unitp,"(' Max Overlap:',2e11.3)") max_Sii,max_Sij
      END IF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) max_Sii,max_Sij
        write(out_unitp,*) 'sub_ana_cplxS'
      END IF
!-----------------------------------------------------------

      RETURN
      end subroutine sub_ana_cplxS

