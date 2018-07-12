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

!=============================================================
!
!      determination des tous les Ln(xi)=serie_fourier(n,i)
!      + les derivees 1er et 2d.
!
!=============================================================
      SUBROUTINE sub_quadra_fourier(base,nosym,nstep,nb_shift,tab_shift)

      USE mod_system
      USE mod_basis
      IMPLICIT NONE

!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis) :: base
      logical       :: num,nosym
      real (kind=Rkind) :: step

      integer  :: nstep,nb_shift
      integer  :: tab_shift(nb_shift)
      integer  :: ib,ishift,nb_i
!---------------------------------------------------------------------
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      logical :: deriv
      real (kind=Rkind) :: d0,d1,d2,d3,dx
      integer       :: i,ii,k,l,nq
!---------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       nq = get_nq_FROM_basis(base)
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING sub_quadra_fourier'
         write(out_unitp,*) 'nb,nq',base%nb,nq
       END IF
!-----------------------------------------------------------

       IF (base%nb <= 0) STOP 'ERROR nb<=0'
      base%packed            = .TRUE.
      base%packed_done       = .TRUE.

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!      test sur nb_fourier et nb_quadra
!      "cos 1 0" or "cos": all terms
!      "cos 2 0" : only sine
!      "cos 2 -1": only cosine
!----------------------------------------------------------------------------
      IF (base%check_nq_OF_basis) THEN
        IF (print_level > -1) write(out_unitp,*) '    Basis: Fourier series'
        IF (print_level > -1) write(out_unitp,*) '      old nb_fourier',base%nb
        IF (mod(base%nb,nb_shift) == 0) THEN
          nb_i = base%nb/nb_shift
        ELSE
          nb_i = base%nb/nb_shift + 1
        END IF
        base%nb = nb_i * nb_shift
        IF (print_level > -1) write(out_unitp,*) '      nb_fourier',base%nb
        IF (print_level > -1) write(out_unitp,*) '      nb_shift:',nb_shift,'tab',tab_shift
      ELSE
        nstep    = 1
        nb_i     = base%nb
        nb_shift = 1
      END IF

      IF (.NOT. base%xPOGridRep_done) THEN
        IF (base%check_nq_OF_basis) THEN
          IF (print_level > -1) write(out_unitp,*) '      old nb_quadra',nq
          IF ( nq < nb_i*nstep ) nq = nb_i*nstep
          IF (nosym .AND. nq == nb_i*nstep .AND. mod(nq,2) == 0) nq = nq + 1
          IF (print_level > -1) write(out_unitp,*) '      new nb_quadra',nq
        END IF

        IF (print_level > -1) write(out_unitp,*) 'fourier: nb,nq',base%nb,nq
        CALL Set_nq_OF_basis(base,nq)
        CALL alloc_xw_OF_basis(base)
        CALL gauss_fourier(base%x(1,:),base%w,nq)
        IF (nosym) THEN
          dx = base%x(1,2)-base%x(1,1)
          base%x = base%x + dx*HALF
        END IF
        base%wrho(:) = base%w(:)
      END IF
      base%rho(:)  = ONE

!     calcul des valeurs de la serie de fourier et des derivees en chaque
!     point de la quadrature.
      deriv = .TRUE.


      CALL alloc_dnb_OF_basis(base)
      ib = 0
      DO ishift=1,nb_shift

        DO i=1,nb_i
          ib = ib + 1
          ii = nstep*i+tab_shift(ishift)
          IF (print_level > -1 .AND. base%print_info_OF_basisDP)  THEN
            IF (mod(ii,2) == 0 .AND. ii < 11) THEN
              write(out_unitp,*) 'basis function: fourrier',ii,'or sin',int(ii/2)
            ELSE IF (mod(ii,2) == 1 .AND. ii < 11) THEN
              write(out_unitp,*) 'basis function: fourrier',ii,'or cos',int(ii/2)
            END IF
          END IF
          base%tab_ndim_index(1,ib) = ii
          DO k=1,nq
            CALL d0d1d2d3fourier(base%x(1,k),d0,d1,d2,d3,ii)
            base%dnRGB%d0(k,ib)     = d0
            base%dnRGB%d1(k,ib,1)   = d1
            base%dnRGB%d2(k,ib,1,1) = d2
!           write(out_unitp,*) ib,k,ii,d0,d1,d2
          END DO
        END DO

      END DO

      IF (base%nb == nq .AND. mod(base%nb,2) == 0) THEN
        base%dnRGB%d0(:,nq) = base%dnRGB%d0(:,nq) / sqrt(TWO)
        base%dnRGB%d1(:,nq,:) = base%dnRGB%d1(:,nq,:) / sqrt(TWO)
        base%dnRGB%d2(:,nq,:,:) = base%dnRGB%d2(:,nq,:,:) / sqrt(TWO)
      END IF

!-----------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(base)
        write(out_unitp,*) 'END sub_quadra_fourier'
      END IF
!-----------------------------------------------------------

      RETURN
      end subroutine sub_quadra_fourier

