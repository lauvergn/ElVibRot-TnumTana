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

!=============================================================
!
!      determination des tous les Ln(xi)=poly_legendre(n,i)
!      + les derivees 1er et 2d.
!
!      paire = 0  => paire
!      paire = 1  => impaire
!      paire =-1  => tout
!
!=============================================================
      SUBROUTINE sub_quadra_legendre(base,paire)
      USE mod_system
      USE mod_basis
      IMPLICIT NONE

!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis)  :: base
      integer       :: paire
!---------------------------------------------------------------------
!---------------------------------------------------------------------

      logical           :: num
      real (kind=Rkind) :: step
      logical           :: deriv
      integer           :: i,j,k,nq
!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       nq = get_nq_FROM_basis(base)
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING sub_quadra_legendre'
         write(out_unitp,*) 'nb,nq',base%nb,nq
         write(out_unitp,*) 'ndim',base%ndim
         write(out_unitp,*) 'paire',paire
       END IF
!-----------------------------------------------------------
       deriv = .TRUE.
       num   = .FALSE.

       IF (base%nb <= 0) STOP 'ERROR nb<=0'

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!      test sur nb_legendre et nb_quadra
!      nb_quadra > nb_legendre
!----------------------------------------------------------------------------
      IF (base%check_nq_OF_basis .AND. print_level > -1) THEN
        write(out_unitp,*) '    Basis: Legendre polynomia'
        write(out_unitp,*) '      nb_legendre',base%nb
      END IF

      base%packed            = .TRUE.
      base%packed_done       = .TRUE.


      IF (.NOT. base%xPOGridRep_done) THEN
        IF (base%check_nq_OF_basis) THEN
          IF (print_level > -1) write(out_unitp,*) '      old nb_quadra',nq
          IF (paire == -1) THEN
            IF ( nq < base%nb ) nq = base%nb
          ELSE
            IF ( nq < 2*base%nb ) nq = 2*base%nb+1
          END IF
          IF (print_level > -1) write(out_unitp,*) '      new nb_quadra',nq
        END IF
        CALL Set_nq_OF_basis(base,nq)

        CALL alloc_xw_OF_basis(base)
        !CALL gauleg(-ONE,ONE,base%x(1,:),base%w,nq)
        CALL gauleg128(-ONE,ONE,base%x(1,:),base%w,nq)
        base%wrho(:) = base%w(:)
      END IF
      base%rho(:)  = ONE

!     calcul des valeurs des polynomes de legendre et des derivees en chaque
!     point de la quadrature.
      deriv = .TRUE.
      CALL alloc_dnb_OF_basis(base)

      IF (paire == 0) THEN
        IF (base%print_info_OF_basisDP .AND. print_level > -1)          &
                       write(out_unitp,*) '     even Legendre polynomia'
        base%tab_ndim_index(1,:) = (/ (2*i-1,i=1,base%nb) /)
        CALL d0d1d2Plm_0_grid(base%x(1,:),base%dnRGB%d0,base%dnRGB%d1,base%dnRGB%d2,   &
                              base%nb,nq,deriv,num,step)

      ELSE IF (paire == 1) THEN
        IF (base%print_info_OF_basisDP .AND. print_level > -1)          &
                        write(out_unitp,*) '     odd Legendre polynomia'
        base%tab_ndim_index(1,:) = (/ (2*i,i=1,base%nb) /)
        CALL d0d1d2Plm_1_grid(base%x(1,:),base%dnRGB%d0,base%dnRGB%d1,base%dnRGB%d2,   &
                              base%nb,nq,deriv,num,step)
      ELSE
        IF (base%print_info_OF_basisDP .AND. print_level > -1)          &
                        write(out_unitp,*) '     All Legendre polynomia'
        base%tab_ndim_index(1,:) = (/ (i,i=1,base%nb) /)
        CALL d0d1d2Plm_grid(base%x(1,:),base%dnRGB%d0,base%dnRGB%d1,base%dnRGB%d2,     &
                            base%nb,nq,deriv,num,step)
      END IF

!-----------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(base)
        write(out_unitp,*) 'END sub_quadra_legendre'
      END IF
!-----------------------------------------------------------
      END SUBROUTINE sub_quadra_legendre
!=============================================================
!
!      determination des tous les Ln(xi)=poly_legendre(n,i)
!      + les derivees 1er et 2d.
!
!=============================================================
      SUBROUTINE transfo_cosTOangle(base)
      USE mod_system
      USE mod_basis
      IMPLICIT NONE

!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis) :: base


      integer i,j
      real(kind=Rkind)  c,s,ss


      DO i=1,get_nq_FROM_basis(base)
        c  = base%x(1,i)
        ss = ONE - c*c
        s  = sqrt(ss)
        base%x(1,i) = acos(c)
        base%w(i)   = base%w(i)/s
        base%rho(i) = s
!       the wrho(:) is unchanged  !!
!       but nrho(1) has to be changed
        base%nrho(1) = 2  ! it means, the variable, x, is substituted by cos(theta)

        DO j=1,base%nb
          base%dnRGB%d2(i,j,1,1) = -c * base%dnRGB%d1(i,j,1) + ss * base%dnRGB%d2(i,j,1,1)
          base%dnRGB%d1(i,j,1)   = -s * base%dnRGB%d1(i,j,1)
        END DO
      END DO

      end subroutine transfo_cosTOangle

