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
!      Particule-in-a-box basis set
!
!=============================================================
      SUBROUTINE sub_quadra_box(base,nosym)
      USE mod_system
      USE mod_basis
      IMPLICIT NONE

!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis)      :: base
      logical           :: nosym


      integer           :: ib
!---------------------------------------------------------------------
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      logical            :: deriv
      real (kind=Rkind)  :: d0,d1,d2,d3
      integer            :: i,ii,k,nq
!---------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_quadra_box'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
       nq = get_nq_FROM_basis(base)
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'nb,nq',base%nb,nq
       END IF
!-----------------------------------------------------------

       IF (base%nb <= 0) STOP 'ERROR nb<=0'

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!      test sur sub_quadra_box et nb_quadra
!----------------------------------------------------------------------------
      base%packed            = .TRUE.
      base%packed_done       = .TRUE.

      IF (base%check_nq_OF_basis) THEN
        write(out_unitp,*) '    Basis: ',name_sub
        write(out_unitp,*) '      nb_box',base%nb
      END IF

      IF (.NOT. base%xPOGridRep_done) THEN
        IF (base%check_nq_OF_basis) THEN
          write(out_unitp,*) '      old nb_quadra',nq
          IF ( nq < base%nb ) nq = base%nb
          write(out_unitp,*) '      new nb_quadra',nq
        END IF
        CALL Set_nq_OF_basis(base,nq)
        CALL alloc_xw_OF_basis(base)

        IF (nosym) THEN
          CALL gauss_box_nosym(base%x(1,:),base%w,nq)
        ELSE
          CALL gauss_box(base%x(1,:),base%w,nq)
        END IF

        base%wrho(:) = base%w(:)
        base%rho(:)  = ONE
      END IF
      base%rho(:)  = ONE


      CALL alloc_dnb_OF_basis(base)


      DO ib=1,base%nb
        base%tab_ndim_index(1,ib) = ib
        IF (debug) write(out_unitp,*) 'basis, particle in a box[0,Pi]:',ib
        DO k=1,nq
          CALL d0d1d2d3box(base%x(1,k),d0,d1,d2,d3,ib)
          base%dnRGB%d0(k,ib)     = d0
          base%dnRGB%d1(k,ib,1)   = d1
          base%dnRGB%d2(k,ib,1,1) = d2
         !write(out_unitp,*) ib,k,d0,d1,d2
        END DO
    END DO

      IF (base%nb == nq) THEN
        base%dnRGB%d0(:,nq) = base%dnRGB%d0(:,nq) / sqrt(TWO)
        base%dnRGB%d1(:,nq,:) = base%dnRGB%d1(:,nq,:) / sqrt(TWO)
        base%dnRGB%d2(:,nq,:,:) = base%dnRGB%d2(:,nq,:,:) / sqrt(TWO)
      END IF

!-----------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(base,write_all=.TRUE.)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_quadra_box

