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
!      determination des tous les Hn(xi)=poly_h(n,i)
!      + les derivees 1er et 2d.
!
!=============================================================
      SUBROUTINE sub_quadra_laguerre(base)

      USE mod_system
      USE mod_basis
      IMPLICIT NONE

!---------- variables passees en argument ----------------------------
      TYPE (basis)     :: base

!---------------------------------------------------------------------
!---------------------------------------------------------------------
      integer  :: i,nq

      integer  :: nb_Laguerre

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_quadra_laguerre'
!      logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
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
!      test sur nb_herm et nb_quadra
!      nb_quadra > nb_herm
!----------------------------------------------------------------------------
      IF (base%check_nq_OF_basis .AND. print_level > -1) THEN
        write(out_unitp,*) '    Basis: Laguerre polynomia'
        write(out_unitp,*) '      nb_Laguerre',base%nb
      END IF
      base%packed            = .TRUE.
      base%packed_done       = .TRUE.


      IF (.NOT. base%xPOGridRep_done) THEN
        IF (base%check_nq_OF_basis) THEN
          IF (print_level > -1) write(out_unitp,*) '      old nb_quadra',nq
          IF ( nq < base%nb ) nq = base%nb + 1
          IF (print_level > -1) write(out_unitp,*) '      new nb_quadra',nq
        END IF
        CALL Set_nq_OF_basis(base,nq)

        CALL alloc_xw_OF_basis(base)

        !CALL hercom(nq,base%x(1,:),base%w)
        base%wrho(:) = base%w(:)
        !base%xPOGridRep_done = .TRUE.
      END IF
      base%rho(:)  = ONE
STOP 'laguerre'

!     calcul des valeurs des polynomes de nb_Laguerre et des derivees en chaque
!     point de la quadrature.
      CALL alloc_dnb_OF_basis(base)

      IF (base%print_info_OF_basisDP .AND. print_level > -1)            &
                        write(out_unitp,*) '      All Laguerre polynomia'
        base%tab_ndim_index(1,:) = (/ (i,i=1,base%nb) /)
        !CALL d0d1d2poly_Laguerre_exp_grille(                            &
        !                     base%x(1,:),                               &
        !                     base%dnRGB%d0(:,:),                             &
        !                     base%dnRGB%d1(:,:,1),                           &
        !                     base%dnRGB%d2(:,:,1,1),                         &
        !                     base%nb,nq)

!-----------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(base)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_quadra_laguerre
