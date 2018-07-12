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
!      determination des tous les Hn(xi)=poly_h(n,i)
!      + les derivees 1er et 2d.
!
!=============================================================
      SUBROUTINE sub_quadra_hermite(base,paire)

      USE mod_system
      USE mod_basis
      IMPLICIT NONE

!---------- variables passees en argument ----------------------------
      TYPE (basis)     :: base
      integer          :: paire
      logical          :: num
      real(kind=Rkind) :: step
!---------------------------------------------------------------------
!---------------------------------------------------------------------
      logical  :: deriv
      integer  :: i,nq,nb_nosym

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_quadra_hermite'
      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
!-----------------------------------------------------------
       nq = get_nq_FROM_basis(base)

       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING sub_quadra_hermite'
         write(out_unitp,*) 'L_SparseBasis',base%L_SparseBasis
         write(out_unitp,*) 'nb',base%nb
         write(out_unitp,*) 'nq',nq
       END IF
!-----------------------------------------------------------
       deriv = .TRUE.
       num   = .FALSE.

       IF (base%nb <= 0) THEN
         base%nb = Get_nb_FROM_l_OF_PrimBasis(base%L_SparseBasis,base)
       END IF

       IF (base%nb <= 0) THEN
         write(out_unitp,*) 'ERROR in sub_quadra_hermite'
         write(out_unitp,*) 'nb<=0',base%nb
         STOP 'ERROR nb<=0'
       END IF

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!      test sur nb_herm et nb_quadra
!      nb_quadra > nb_herm
!----------------------------------------------------------------------------
      IF (base%check_nq_OF_basis .AND. print_level > -1) THEN
        write(out_unitp,*) '    Basis: Hermite polynomia'
        write(out_unitp,*) '      nb_hermite',base%nb
      END IF
      IF (print_level > -1) write(out_unitp,*) '      nb_hermite',base%nb

      base%packed            = .TRUE.
      base%packed_done       = .TRUE.


      IF (.NOT. base%xPOGridRep_done) THEN
        IF (base%check_nq_OF_basis) THEN
          IF (print_level > -1) write(out_unitp,*) '      nb_hermite',base%nb
          IF (print_level > -1) write(out_unitp,*) '      old nb_quadra',nq
          IF (paire == -1) THEN
            nb_nosym = base%nb
          ELSE IF (paire == 1) THEN ! odd
            nb_nosym = 2*base%nb - 1
          ELSE ! even
            nb_nosym = 2*base%nb
          END IF
        END IF

        SELECT CASE (base%Nested)
        CASE(1)

          IF (mod(base%nq_max_Nested,2) == 0) base%nq_max_Nested = base%nq_max_Nested + 1

          IF (base%check_nq_OF_basis .AND. mod(nq,2) == 0)    nq = nq + 1 ! here nq must be odd : 2n-1
          IF (base%check_nq_OF_basis .AND. nq < 2*nb_nosym-1) nq = 2*nb_nosym-1
          CALL Set_nq_OF_basis(base,nq)
          IF (print_level > -1) write(out_unitp,*) '      new nb_quadra',nq

          CALL alloc_xw_OF_basis(base)
          CALL grid_HermiteNested1(base%x,base%w,nq,base%nq_max_Nested)

        CASE(2) ! not yet
          IF (mod(base%nq_max_Nested,2) == 0) base%nq_max_Nested = base%nq_max_Nested + 1

          IF (base%check_nq_OF_basis .AND. mod(nq,2) == 0)    nq = nq + 1 ! here nq must be odd : 2n-1
          IF (base%check_nq_OF_basis .AND. nq < 2*nb_nosym-1) nq = 2*nb_nosym-1
          CALL Set_nq_OF_basis(base,nq)
          IF (print_level > -1) write(out_unitp,*) '      new nb_quadra',nq

          CALL alloc_xw_OF_basis(base)

          STOP 'not yet nested2'
        CASE Default

          IF (base%check_nq_OF_basis .AND. nq < nb_nosym) nq = nb_nosym + 1
          CALL Set_nq_OF_basis(base,nq)
          IF (print_level > -1) write(out_unitp,*) '      new nb_quadra',nq

          CALL alloc_xw_OF_basis(base)
          CALL hercom(nq,base%x(1,:),base%w)

        END SELECT
        !base%xPOGridRep_done = .TRUE.
        base%wrho(:) = base%w(:)
        base%rho(:)  = ONE
      END IF


!     calcul des valeurs des polynomes de hermites et des derivees en chaque
!     point de la quadrature.
      deriv = .TRUE.
      CALL alloc_dnb_OF_basis(base)

      IF (paire == 0) THEN
        IF (base%print_info_OF_basisDP .AND. print_level > -1)          &
                       write(out_unitp,*) '      even Hermite polynomia'
        base%tab_ndim_index(1,:) = (/ (2*i-1,i=1,base%nb) /)
        CALL d0d1d2poly_Hermite_0_exp_grille(                           &
                             base%x(1,:),                               &
                             base%dnRGB%d0(:,:),                             &
                             base%dnRGB%d1(:,:,1),                           &
                             base%dnRGB%d2(:,:,1,1),                         &
                             base%nb,nq,deriv,num,step)
      ELSE IF (paire == 1) THEN
        IF (base%print_info_OF_basisDP .AND. print_level > -1)          &
                         write(out_unitp,*) '      odd Hermite polynomia'
        base%tab_ndim_index(1,:) = (/ (2*i,i=1,base%nb) /)
        CALL d0d1d2poly_Hermite_1_exp_grille(                           &
                             base%x(1,:),                               &
                             base%dnRGB%d0(:,:),                             &
                             base%dnRGB%d1(:,:,1),                           &
                             base%dnRGB%d2(:,:,1,1),                         &
                             base%nb,nq,deriv,num,step)
      ELSE
        IF (base%print_info_OF_basisDP .AND. print_level > -1)          &
                        write(out_unitp,*) '      All Hermite polynomia'
        base%tab_ndim_index(1,:) = (/ (i,i=1,base%nb) /)
        CALL d0d1d2poly_Hermite_exp_grille(                             &
                             base%x(1,:),                               &
                             base%dnRGB%d0(:,:),                             &
                             base%dnRGB%d1(:,:,1),                           &
                             base%dnRGB%d2(:,:,1,1),                         &
                             base%nb,nq,deriv,num,step)
      END IF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'base%dnRGB'
        CALL write_dnSVM(base%dnRGB)
        !CALL RecWrite_basis(base,write_all=.FALSE.)
        write(out_unitp,*) 'END sub_quadra_hermite'
      END IF
!-----------------------------------------------------------

      end subroutine sub_quadra_hermite
!=============================================================
!
!      determination des tous les Hn(xi)=poly_h(n,i)
!      + les derivees 1er et 2d.
!
!=============================================================
      SUBROUTINE sub_quadra_hermitebox(base,paire)

      USE mod_system
      USE mod_basis
      IMPLICIT NONE

!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis) :: base
      integer paire
      logical num
      real(kind=Rkind)  step
!---------------------------------------------------------------------
!---------------------------------------------------------------------
      logical  :: deriv
      real(kind=Rkind) :: A,B,xeq,scale,dx
      integer  :: nq,i,j,nqo

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_quadra_hermitebox'
!      logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
!-----------------------------------------------------------
       nqo = get_nq_FROM_basis(base)
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING sub_quadra_hermitebox'
         write(out_unitp,*) 'nb,nq',base%nb,nqo
         write(out_unitp,*) 'num,step',num,step
       END IF
!-----------------------------------------------------------

       deriv = .TRUE.
       num   = .FALSE.

       IF (base%nb <= 0) STOP 'ERROR nb<=0'

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!      test sur nb_herm et nb_quadra
!      nb_quadra > nb_herm
!----------------------------------------------------------------------------
      IF (base%print_info_OF_basisDP .AND. print_level > -1)            &
                 write(out_unitp,*) '    Basis: Hermite polynomia+BoxAB'
      base%packed            = .TRUE.
      base%packed_done       = .TRUE.


      IF (base%xPOGridRep_done) THEN
        write(out_unitp,*) 'ERROR in sub_quadra_hermitebox'
        write(out_unitp,*) 'xPOGridRep_done=t is not possible for this basis!'
        write(out_unitp,*) 'CHECK your data'
        STOP
      END IF
      IF (base%check_nq_OF_basis) THEN
        IF (print_level > -1) write(out_unitp,*) '      nb_hermite',base%nb
        IF (print_level > -1) write(out_unitp,*) '      old nb_quadra',nqo
        IF (paire .EQ. -1) THEN
          IF ( nqo < base%nb ) nqo = base%nb + 1
        ELSE
          IF ( nqo < 2*base%nb ) nqo = 2*base%nb+1
        END IF
        IF (print_level > -1) write(out_unitp,*) '      new nb_quadra',nqo
      END IF
      CALL Set_nq_OF_basis(base,nqo)

      CALL alloc_xw_OF_basis(base)
!     - The true Gauss-Hermite grid -----------------------------------------
!       with a small number of grid points (nb)
      nq = base%nb
      CALL hercom(nq,base%x(1,:),base%w)
!     - The true Gauss-Hermite grid -----------------------------------------


!     - We add grid points between A and B (regularly spaced) ---------------
!       without weight (0.)
!     ------ for the domain [A,B] -------------------------------------------
      A = base%x(1,1)
      B = base%x(1,nq)
      IF (base%print_info_OF_basisDP .AND. print_level > -1)            &
                                   write(out_unitp,*) 'HermBox: A,B',A,B
      dx = (B-A)/real(nqo-nq,kind=Rkind)
      base%x(1,nq+1:nqo) =                                              &
         (/ (A+dx*(real(i,kind=Rkind)-HALF),i=1,nqo-nq) /)
      base%w(nq+1:nqo) = ZERO


!     - sort the grid points --------------------------------------------------
      DO i=1,nqo
      DO j=i+1,nqo
        IF (base%x(1,i) > base%x(1,j) ) THEN
          A           = base%x(1,j)
          base%x(1,j) = base%x(1,i)
          base%x(1,i) = A
          A         = base%w(j)
          base%w(j) = base%w(i)
          base%w(i) = A
        END IF
      END DO
      END DO
!     - We add grid points between A and B (regularly spaced) ---------------


      base%wrho(:) = base%w(:)
      base%rho(:)  = ONE

!     calcul des valeurs des polynomes de hermites et des derivees en chaque
!     point de la quadrature.
      deriv = .TRUE.
      CALL alloc_dnb_OF_basis(base)

      IF (paire == 0) THEN ! even
        IF (base%print_info_OF_basisDP .AND. print_level > -1)          &
                       write(out_unitp,*) '      even Hermite polynomia'
        base%tab_ndim_index(1,:) = (/ (2*i-1,i=1,base%nb) /)
        CALL d0d1d2poly_Hermite_0_exp_grille(                           &
                             base%x(1,:),                               &
                             base%dnRGB%d0(:,:),                             &
                             base%dnRGB%d1(:,:,1),                           &
                             base%dnRGB%d2(:,:,1,1),                         &
                             base%nb,nqo,deriv,num,step)
      ELSE IF (paire == 1) THEN ! odd
        IF (base%print_info_OF_basisDP .AND. print_level > -1)          &
                        write(out_unitp,*) '      odd Hermite polynomia'
        base%tab_ndim_index(1,:) = (/ (2*i,i=1,base%nb) /)
        CALL d0d1d2poly_Hermite_1_exp_grille(                           &
                             base%x(1,:),                               &
                             base%dnRGB%d0(:,:),                             &
                             base%dnRGB%d1(:,:,1),                           &
                             base%dnRGB%d2(:,:,1,1),                         &
                             base%nb,nqo,deriv,num,step)
      ELSE
        IF (base%print_info_OF_basisDP .AND. print_level > -1)          &
                        write(out_unitp,*) '      All Hermite polynomia'
        base%tab_ndim_index(1,:) = (/ (i,i=1,base%nb) /)
        CALL d0d1d2poly_Hermite_exp_grille(                             &
                             base%x(1,:),                               &
                             base%dnRGB%d0(:,:),                             &
                             base%dnRGB%d1(:,:,1),                           &
                             base%dnRGB%d2(:,:,1,1),                         &
                             base%nb,nqo,deriv,num,step)
      END IF


!-----------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(base)
        write(out_unitp,*) 'END sub_quadra_hermite'
      END IF
!-----------------------------------------------------------

      end subroutine sub_quadra_hermitebox


      SUBROUTINE grid_HermiteNested1(x,w,nq,nqmax)

      USE mod_system
      USE mod_basis
      IMPLICIT NONE

      integer, intent(in) :: nq,nqmax
      real(kind=Rkind), intent(inout) :: x(nq)
      real(kind=Rkind), intent(inout) :: w(nq)

!---------------------------------------------------------------------
!---------- local variables ----------------------------
      TYPE (basis)     :: base
      logical          :: num,deriv
      real(kind=Rkind) :: step
!---------------------------------------------------------------------
!---------------------------------------------------------------------
      !real(kind=Rkind) :: A,B,xeq,scale,dx
      !integer  :: nq,iq,i,j,nb_nosym,nq0,dnq,nqo
      real(kind=Rkind), allocatable :: x_loc(:)
      real(kind=Rkind), allocatable :: w_loc(:)
      integer  :: dnb


!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='grid_HermiteNested1'
!      logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nq,nqmax',nq,nqmax
      END IF
!-----------------------------------------------------------

      deriv = .TRUE.
      num   = .FALSE.

      IF (nq > nqmax .OR. mod(nqmax-nq,2) /= 0) THEN
        write(out_unitp,*) 'nq,nqmax',nq,nqmax
        STOP 'ERROR in grid_HermiteNested1'
      END IF


      CALL alloc_NParray(x_loc,(/ nqmax /),"x_loc",name_sub)
      CALL alloc_NParray(w_loc,(/ nqmax /),"w_loc",name_sub)
      CALL hercom(nqmax,x_loc(:),w_loc(:))
      !write(out_unitp,*) 'old w(:)',w_loc(:)
      !write(out_unitp,*) 'old x(:)',x_loc(:)

      ! new grid
      dnb = (nqmax-nq)/2
      x(:) = x_loc(1+dnb:nqmax-dnb)


      !weight
      CALL Set_nq_OF_basis(base,nq)
      base%nb   = (nq+1)/2
      base%ndim = 1
      CALL alloc_dnb_OF_basis(base)
      CALL d0d1d2poly_Hermite_exp_grille(                               &
                             x(:),                               &
                             base%dnRGB%d0(:,:),                             &
                             base%dnRGB%d1(:,:,1),                           &
                             base%dnRGB%d2(:,:,1,1),                         &
                             base%nb,nq,deriv,num,step)

      CALL Weight_OF_grid(w,base%dnRGB%d0,base%nb,nq)

      CALL dealloc_dnb_OF_basis(base)

      CALL dealloc_NParray(x_loc,"x_loc",name_sub)
      CALL dealloc_NParray(w_loc,"w_loc",name_sub)

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!      test sur nb_herm et nb_quadra
!      nb_quadra > nb_herm
!----------------------------------------------------------------------------

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE grid_HermiteNested1

      SUBROUTINE sub_quadra_HermiteNested2(base,paire)

      USE mod_system
      USE mod_basis
      IMPLICIT NONE

!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis) :: base
      integer paire
      logical num
      real(kind=Rkind)  step
!---------------------------------------------------------------------
!---------------------------------------------------------------------
      logical  :: deriv
      real(kind=Rkind) :: A,B,xeq,scale,dx
      integer  :: nb,nq,iq,i,j,nb_nosym,nq0,dnq,nq1,nq2,nqo
      real(kind=Rkind), allocatable :: x_loc(:)
      real(kind=Rkind), allocatable :: w_loc(:)


!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_quadra_HermiteNested2'
!      logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
!-----------------------------------------------------------
       nqo = get_nq_FROM_basis(base)
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb,nq',base%nb,nqo
        write(out_unitp,*) 'num,step',num,step
      END IF
!-----------------------------------------------------------


      deriv = .TRUE.
      num   = .FALSE.

      IF (base%nb <= 0) STOP 'ERROR nb<=0'

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!      test sur nb_herm et nb_quadra
!      nb_quadra > nb_herm
!----------------------------------------------------------------------------
      IF (base%print_info_OF_basisDP) write(out_unitp,*) '    Basis: Hermite polynomia+Nested'
      base%packed            = .TRUE.
      base%packed_done       = .TRUE.


      IF (base%xPOGridRep_done) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'xPOGridRep_done=t is not possible for this basis!'
        write(out_unitp,*) 'CHECK your data'
        STOP
      END IF
      IF (base%check_nq_OF_basis) THEN
        write(out_unitp,*) '      nb_hermite',base%nb
        write(out_unitp,*) '      old nb_quadra',nqo
        IF (paire == -1) THEN
          nb_nosym = base%nb
        ELSE IF (paire == 1) THEN ! even
          nb_nosym = 2*base%nb - 1
        ELSE ! no sym
          nb_nosym = 2*base%nb
        END IF
        IF (mod(nqo,2) == 0) nqo = nqo + 1
        IF (nqo < 2*nb_nosym-1) THEN
          nqo = 2*nb_nosym-1
        END IF
        write(out_unitp,*) '      new nb_quadra',nqo
      END IF

      IF (base%L_SparseGrid < 0) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'base%L_SparseGrid MUST larger than -1',base%L_SparseGrid
        write(out_unitp,*) 'CHECK the fortran !!'
        STOP
      END IF
      IF (base%nq_max_Nested < 1) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'nq_max_Nested < 1',base%nq_max_Nested
        write(out_unitp,*) 'CHECK your data or the fortran!!'
        STOP
      END IF
      IF (mod(base%nq_max_Nested,2) == 0) base%nq_max_Nested = base%nq_max_Nested + 1

      ! nested Hermite Grid (here without weight)
      IF (base%L_SparseGrid <= base%nq_max_Nested/2) THEN
        nqo      = base%L_SparseGrid * 2 + 1
        nb_nosym = base%L_SparseGrid + 1
        nq0      = base%nq_max_Nested/2 + 1
        dnq      = nqo/2
        nq1      = nq0-dnq
        nq2      = nq0+dnq
      ELSE
        nqo      = base%nq_max_Nested
        nb_nosym = base%nq_max_Nested
        nq1      = 1
        nq2      = base%nq_max_Nested
      END IF
      write(out_unitp,*) 'base%nq_max_Nested',base%nq_max_Nested
      write(out_unitp,*) 'base%L_SparseGrid,nq,nb_nosym,nq1,nq2',base%L_SparseGrid,nqo,nb_nosym,nq1,nq2
        CALL Set_nq_OF_basis(base,nqo)

      CALL alloc_xw_OF_basis(base)

      CALL alloc_NParray(x_loc,(/ base%nq_max_Nested /),                  &
                        "x_loc","sub_quadra_HermiteNested2")
      CALL alloc_NParray(w_loc,(/ base%nq_max_Nested /),                  &
                        "w_loc","sub_quadra_HermiteNested2")

      CALL hercom(base%nq_max_Nested,x_loc(:),w_loc(:))
      !write(out_unitp,*) 'old w(:)',w_loc(:)
      !write(out_unitp,*) 'old x(:)',x_loc(:)
      base%x(1,:) = x_loc(nq1:nq2)
      write(out_unitp,*) 'new x(:)',base%x(1,:)



      base%rho(:)  = ONE

      ! weights calculation
      nb = base%nb ! save nb
      base%nb = nb_nosym
      deriv = .FALSE.

      CALL alloc_dnb_OF_basis(base)
      CALL d0d1d2poly_Hermite_exp_grille(                               &
                             base%x(1,:),                               &
                             base%dnRGB%d0(:,:),                             &
                             base%dnRGB%d1(:,:,1),                           &
                             base%dnRGB%d2(:,:,1,1),                         &
                             base%nb,nqo,deriv,num,step)

      CALL Weight_OF_grid_basis(base)
      CALL dealloc_dnb_OF_basis(base)
      base%nb = nb
      write(out_unitp,*) 'new w',base%w

!     calcul des valeurs des polynomes de hermites et des derivees en chaque
!     point de la quadrature.
      deriv = .TRUE.
      CALL alloc_dnb_OF_basis(base)

      IF (paire == 0) THEN
        IF (base%print_info_OF_basisDP) write(out_unitp,*) '      even Hermite polynomia'
        base%tab_ndim_index(1,:) = (/ (2*i-1,i=1,base%nb) /)
        CALL d0d1d2poly_Hermite_0_exp_grille(                           &
                             base%x(1,:),                               &
                             base%dnRGB%d0(:,:),                             &
                             base%dnRGB%d1(:,:,1),                           &
                             base%dnRGB%d2(:,:,1,1),                         &
                             base%nb,nqo,deriv,num,step)
      ELSE IF (paire == 1) THEN
        IF (base%print_info_OF_basisDP) write(out_unitp,*) '      odd Hermite polynomia'
        base%tab_ndim_index(1,:) = (/ (2*i,i=1,base%nb) /)
        CALL d0d1d2poly_Hermite_1_exp_grille(                           &
                             base%x(1,:),                               &
                             base%dnRGB%d0(:,:),                             &
                             base%dnRGB%d1(:,:,1),                           &
                             base%dnRGB%d2(:,:,1,1),                         &
                             base%nb,nqo,deriv,num,step)
      ELSE
        IF (base%print_info_OF_basisDP) write(out_unitp,*) '      All Hermite polynomia'
        base%tab_ndim_index(1,:) = (/ (i,i=1,base%nb) /)
        CALL d0d1d2poly_Hermite_exp_grille(                             &
                             base%x(1,:),                               &
                             base%dnRGB%d0(:,:),                             &
                             base%dnRGB%d1(:,:,1),                           &
                             base%dnRGB%d2(:,:,1,1),                         &
                             base%nb,nqo,deriv,num,step)
      END IF

      CALL dealloc_NParray(x_loc,"x_loc","sub_quadra_HermiteNested2")
      CALL dealloc_NParray(w_loc,"w_loc","sub_quadra_HermiteNested2")

!-----------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(base)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_quadra_HermiteNested2

     SUBROUTINE sub_quadra_hermite_cuba(base)

      USE mod_system
      USE mod_basis
      IMPLICIT NONE

!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis)     :: base
      logical          :: num = .FALSE.
      real(kind=Rkind) :: step = ZERO
!---------------------------------------------------------------------
!---------------------------------------------------------------------
      logical          :: deriv = .TRUE.
      integer          :: p,d,option
      integer          :: tab_nb(base%ndim)
      real(kind=Rkind) :: x(base%ndim)
      real(kind=Rkind) :: d0(base%ndim)
      real(kind=Rkind) :: d1(base%ndim)
      real(kind=Rkind) :: d2(base%ndim)
      integer          :: i,iB,iQ,id,jd,nq,LB,LG


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_quadra_hermite_cuba'
      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
!-----------------------------------------------------------
       nq = get_nq_FROM_basis(base)
       IF (base%L_SparseBasis < 0) THEN
         LB = base%nb-1
       ELSE
         LB = base%L_SparseBasis
       END IF

       IF (base%L_SparseGrid < 0) THEN
         LG = nq-1
       ELSE
         LG = base%L_SparseGrid
       END IF

       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING sub_quadra_hermite_cuba'
         write(out_unitp,*) 'LB,LG',LB,LG
         write(out_unitp,*) 'ndim',base%ndim
       END IF
!-----------------------------------------------------------
       IF (base%ndim < 1) STOP ' Problem with ndim'

       IF (base%ndim == 1) STOP ' You MUST use "HO" basis set'

       ! ndim is the "n" parameter of Burkardt subroutine
       ! d is the largest degree of the polynomials: d=2*(nq-1) = 2*LG
       IF (LB < 0) STOP 'ERROR LB<=0'

       ! test if nq >= nb
      IF (base%check_nq_OF_basis) THEN
        write(out_unitp,*) '    Basis: Hermite polynomia, cubature'
        write(out_unitp,*) '      LB,LG',LB,LG
      END IF
      base%packed            = .TRUE.
      base%packed_done       = .TRUE.


      IF (.NOT. base%xPOGridRep_done) THEN
        IF (base%check_nq_OF_basis) THEN
          IF (print_level > -1) write(out_unitp,*) '      old LG',LG
          IF ( LG < LB ) LG = LB + 1
          IF (print_level > -1) write(out_unitp,*) '      new LG',LG
          IF (print_level > -1) write(out_unitp,*) '           d',2*LG
        END IF

        d = 2*LG ! degree
        ! calculation of the grid point number (with the EN_02... size subroutine)
        SELECT CASE (d)
        CASE(0)
          IF (print_level > -1) write(out_unitp,*) ' en_r2_01_1_size'
          CALL en_r2_01_1_size(base%ndim,nq)
        CASE(2)
          IF (print_level > -1) write(out_unitp,*) ' en_r2_02_xiu_size'
          CALL en_r2_02_xiu_size(base%ndim,nq)
        CASE(4)
          IF (print_level > -1) write(out_unitp,*) ' en_r2_05_1_size'
          option = 1
          CALL en_r2_05_1_size(base%ndim,option,nq)
        CASE(6)
          IF (print_level > -1) write(out_unitp,*) ' en_r2_07_3_size'
          option = 1
          CALL en_r2_07_3_size(base%ndim,option,nq)
        CASE(8)
          IF (print_level > -1) write(out_unitp,*) ' en_r2_09_1_size'
          option = 1
          CALL en_r2_09_1_size(base%ndim,option,nq)
        CASE(10)
         IF (print_level > -1)  write(out_unitp,*) ' en_r2_11_1_size'
          option = 1
          CALL en_r2_11_1_size(base%ndim,option,nq)
        CASE DEFAULT
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' we cannot use cubature with d > 10'
          write(out_unitp,*) 'd',d
          write(out_unitp,*) 'nb,nq',base%nb,nq
          write(out_unitp,*) 'ndim',base%ndim
          nq = -1
          !RETURN
        END SELECT
        IF (base%print_info_OF_basisDP .AND. print_level > -1) THEN
          write(out_unitp,*) 'maximum degree:',d
          write(out_unitp,*) 'Cubature point number:',nq
        END IF
        CALL Set_nq_OF_basis(base,nq)


        IF (nq < 1) RETURN
        CALL alloc_xw_OF_basis(base)


        ! the grid points
        SELECT CASE (d)
        CASE(0)
          CALL en_r2_01_1(base%ndim,nq,base%x,base%w)
        CASE(2)
          CALL en_r2_02_xiu(base%ndim,nq,base%x,base%w)
        CASE(4)
          option = 1
          CALL en_r2_05_1(base%ndim,option,nq,base%x,base%w)
        CASE(6)
          option = 1
          CALL en_r2_07_3(base%ndim,option,nq,base%x,base%w)
        CASE(8)
          option = 1
          CALL en_r2_09_1(base%ndim,option,nq,base%x,base%w)
        CASE(10)
          option = 1
          CALL en_r2_11_1(base%ndim,option,nq,base%x,base%w)
        CASE DEFAULT
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' we cannot use cubature with d > 10'
          write(out_unitp,*) 'nb,nq',base%nb,nq
          write(out_unitp,*) 'ndim',base%ndim
          nq = -1
          STOP
        END SELECT

        DO iQ=1,nq
          x(:) = base%x(:,iQ)
          base%w(iQ) = base%w(iQ) / exp(-dot_product(x,x))
        END DO

        IF (base%print_info_OF_basisDP .AND. print_level > -1) THEN
          write(out_unitp,*) 'cubature',nq
          DO iQ=1,nq
            write(6,*) iQ,base%x(:,iQ),base%w(iQ)
          END DO
        END IF

        base%wrho(:) = base%w(:)
        !base%xPOGridRep_done = .TRUE.
      END IF
      base%rho(:)  = ONE



      ! number of basis functions:

      CALL dealloc_nDindex(base%nDindB)
      CALL alloc_nDindex(base%nDindB,ndim=base%ndim)

      base%nDindB%packed = .TRUE.
      tab_nb(:) = LB+1
      CALL init_nDindexPrim(base%nDindB,ndim=base%ndim,                 &
                            Type_OF_nDindex=0,                          &
                            MaxNorm=real(LB,kind=Rkind),                &
                            nDsize=tab_nb)
      IF (debug) CALL Write_nDindex(base%nDindB,name_sub)
      base%nb = base%nDindB%max_nDI
      IF (base%print_info_OF_basisDP .AND. print_level > -1) THEN
        write(out_unitp,*) 'Basis function number:',base%nb
      END IF


!     calcul des valeurs des polynomes de hermites et des derivees en chaque
!     point de la quadrature.
      deriv = .TRUE.
      CALL alloc_dnb_OF_basis(base)


      DO iQ=1,nq
        x(:) = base%x(:,iQ)
        DO iB=1,base%nb
          CALL calc_nDindex(base%nDindB,iB,tab_nb)

          base%tab_ndim_index(:,iB) = tab_nb(:)


          DO i=1,base%ndim
            CALL d0d1d2poly_Hermite_exp(x(i),tab_nb(i)-1,               &
                                       d0(i),d1(i),d2(i),deriv,num,step)
          END DO
          base%dnRGB%d0(iQ,iB) = product(d0)

          DO id=1,base%ndim
            base%dnRGB%d1(iQ,iB,id)    = d1(id)
            base%dnRGB%d2(iQ,iB,id,id) = d2(id)
            DO i=1,base%ndim
              IF (i == id) CYCLE
              base%dnRGB%d1(iQ,iB,id)    = base%dnRGB%d1(iQ,iB,id)    * d0(i)
              base%dnRGB%d2(iQ,iB,id,id) = base%dnRGB%d2(iQ,iB,id,id) * d0(i)
            END DO
          END DO

          DO id=1,base%ndim
          DO jd=id+1,base%ndim

            base%dnRGB%d2(iQ,iB,id,jd) = d1(id) * d1(jd)
            DO i=1,base%ndim
              IF (i == id .OR. i ==jd) CYCLE
              base%dnRGB%d2(iQ,iB,id,jd) = base%dnRGB%d2(iQ,iB,id,jd) * d0(i)
            END DO
            base%dnRGB%d2(iQ,iB,jd,id) = base%dnRGB%d2(iQ,iB,id,jd)
          END DO
          END DO

        END DO
      END DO

!      CALL sort_basis(base)
!      base%check_basis = .TRUE.
!      write(out_unitp,*) 'coucou ',name_sub
!      CALL check_ortho_basis(base,.FALSE.)
!      base%check_basis = .FALSE.
!      write(out_unitp,*) 'old w',base%w
!      CALL Weight_OF_grid_basis(base)
!      write(out_unitp,*) 'new w',base%w

!-----------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(base)
        write(out_unitp,*) 'END sub_quadra_hermite_cuba'
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_quadra_hermite_cuba
      SUBROUTINE sub_quadra_hermite_cuba_DML(base,err_grid)

      USE mod_system
      USE mod_basis
      IMPLICIT NONE

!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis)     :: base
      logical          :: num = .FALSE.
      real(kind=Rkind) :: step = ZERO
!---------------------------------------------------------------------
!---------------------------------------------------------------------
      logical          :: deriv = .TRUE.
      integer          :: p,d,option
      integer          :: tab_nb(base%ndim)
      real(kind=Rkind) :: x(base%ndim)
      real(kind=Rkind) :: d0(base%ndim)
      real(kind=Rkind) :: d1(base%ndim)
      real(kind=Rkind) :: d2(base%ndim)
      integer          :: i,iB,iQ,id,jd,idum,nq


      integer                  :: LB,LG
      character (len=Name_len) :: name_i,name_j
      integer                  :: nio,err_io
      TYPE (param_file)        :: cubature_file
      logical                  :: err_grid

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_quadra_hermite_cuba_DML'
      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
!-----------------------------------------------------------
       nq = get_nq_FROM_basis(base)
       IF (base%L_SparseBasis < 0) THEN
         LB = base%nb-1
       ELSE
         LB = base%L_SparseBasis
       END IF

       IF (base%L_SparseGrid < 0) THEN
         LG = nq-1
       ELSE
         LG = base%L_SparseGrid
       END IF

       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) '    Basis: Hermite polynomia, cubture'
         write(out_unitp,*) '      LB,LG',LB,LG
       END IF
!-----------------------------------------------------------
       err_Grid = .FALSE.

       IF (base%ndim < 1) STOP ' Problem with ndim'

       IF (base%ndim == 1) STOP ' You MUST use "HO" basis set'

       ! ndim is the "n" parameter of Burkardt subroutine
       ! d is the largest degree of the polynomials: d=2*(nq-1)
       IF (LB < 0) STOP 'ERROR LB<0'

       ! test if nq >= nb
      IF (base%check_nq_OF_basis) THEN
        write(out_unitp,*) '    Basis: Hermite polynomia, cubature'
        write(out_unitp,*) '      LB,LG',LB,LG
      END IF
      base%packed            = .TRUE.
      base%packed_done       = .TRUE.


      IF (.NOT. base%xPOGridRep_done) THEN
        IF (base%check_nq_OF_basis) THEN
          IF (print_level > -1) write(out_unitp,*) '      old LG',LG
          IF (LG < LB) LG = LB + 1
          IF (print_level > -1) write(out_unitp,*) '      new LG',LG
        END IF

        IF (LG /= 0) THEN
          CALL Write_int_IN_char(base%ndim,name_i)
          CALL Write_int_IN_char(LG,       name_j)
          cubature_file%name = trim(EVRT_path) //                       &
                                          '/Internal_data/HermCuba/' // &
                                 trim(name_i) // 'D_deg' // trim(name_j)


          write(out_unitp,*) 'cubature_file%name: ',cubature_file%name
          CALL file_open(cubature_file,nio,old=.TRUE.,err_file=err_io)

          IF (err_io == 0) THEN
            read(nio,*,iostat=err_io) nq
          ELSE
            write(out_unitp,*) 'WARNNING in ',name_sub
            write(out_unitp,*) 'cannot find cubature file: ',cubature_file%name
            err_Grid = .TRUE.
          END IF
          write(out_unitp,*) '      cubature nq',nq
          close(nio)
          CALL Set_nq_OF_basis(base,nq)


          err_Grid = err_Grid .OR. (nq < 1)
          IF (err_Grid) RETURN



          CALL alloc_xw_OF_basis(base)

          CALL file_open(cubature_file,nio,old=.TRUE.,err_file=err_io)
          read(nio,*,iostat=err_io)

          DO iq=1,nq
            read(nio,*,iostat=err_io) idum,base%x(:,iq),base%w(iq)
          END DO
          close(nio)
        ELSE ! LG = 0 (degree zero => only one point (0., 0., ...0.)
          nq = 1
          CALL Set_nq_OF_basis(base,nq)

          CALL alloc_xw_OF_basis(base)
          base%x(:,1) = ZERO
          base%w(1)   = sqrt(pi)**base%ndim
        END IF
        base%wrho(:) = base%w(:)
        !base%xPOGridRep_done = .TRUE.
      END IF
      base%rho(:)  = ONE


      ! number of basis functions:

      CALL dealloc_nDindex(base%nDindB)
      CALL alloc_nDindex(base%nDindB,ndim=base%ndim)

      base%nDindB%packed = .TRUE.
      tab_nb(:) = LB+1
      CALL init_nDindexPrim(base%nDindB,ndim=base%ndim,                 &
                            Type_OF_nDindex=0,                          &
                            MaxNorm=real(LB,kind=Rkind),                &
                            nDsize=tab_nb)
      IF (debug) CALL Write_nDindex(base%nDindB,name_sub)
      base%nb = base%nDindB%max_nDI
      IF (base%print_info_OF_basisDP) write(out_unitp,*) 'Basis function number:',base%nb


!     calcul des valeurs des polynomes de hermites et des derivees en chaque
!     point de la quadrature.
      deriv = .TRUE.
      CALL alloc_dnb_OF_basis(base)


      DO iQ=1,nq
        x(:) = base%x(:,iQ)
        DO iB=1,base%nb
          CALL calc_nDindex(base%nDindB,iB,tab_nb)

          base%tab_ndim_index(:,iB) = tab_nb(:)


          DO i=1,base%ndim
            CALL d0d1d2poly_Hermite_exp(x(i),tab_nb(i)-1,               &
                                       d0(i),d1(i),d2(i),deriv,num,step)
          END DO
          base%dnRGB%d0(iQ,iB) = product(d0)

          DO id=1,base%ndim
            base%dnRGB%d1(iQ,iB,id)    = d1(id)
            base%dnRGB%d2(iQ,iB,id,id) = d2(id)
            DO i=1,base%ndim
              IF (i == id) CYCLE
              base%dnRGB%d1(iQ,iB,id)    = base%dnRGB%d1(iQ,iB,id)    * d0(i)
              base%dnRGB%d2(iQ,iB,id,id) = base%dnRGB%d2(iQ,iB,id,id) * d0(i)
            END DO
          END DO

          DO id=1,base%ndim
          DO jd=id+1,base%ndim

            base%dnRGB%d2(iQ,iB,id,jd) = d1(id) * d1(jd)
            DO i=1,base%ndim
              IF (i == id .OR. i ==jd) CYCLE
              base%dnRGB%d2(iQ,iB,id,jd) = base%dnRGB%d2(iQ,iB,id,jd) * d0(i)
            END DO
            base%dnRGB%d2(iQ,iB,jd,id) = base%dnRGB%d2(iQ,iB,id,jd)
          END DO
          END DO

        END DO
      END DO

!      CALL sort_basis(base)
!      base%check_basis = .TRUE.
!      write(out_unitp,*) 'coucou ',name_sub
!      CALL check_ortho_basis(base,.FALSE.)
!      base%check_basis = .FALSE.
!      write(out_unitp,*) 'old w',base%w
!      CALL Weight_OF_grid_basis(base)
!      write(out_unitp,*) 'new w',base%w

!-----------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(base)
        write(out_unitp,*) 'END sub_quadra_hermite_cuba_DML'
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_quadra_hermite_cuba_DML

      SUBROUTINE transfo_Q_TO_tQ(base)
      USE mod_system
      USE mod_basis
      USE mod_dnSVM
      IMPLICIT NONE

!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis) :: base


      TYPE (Type_dnS)    :: dntQ,dntQ_inv
      integer            :: i,dnErr
      real(kind=Rkind)   :: q,d1,d2
      integer, parameter :: itype = 2

      CALL alloc_dnS(dntQ,nb_var_deriv=1,nderiv=3)
      CALL alloc_dnS(dntQ_inv,nb_var_deriv=1,nderiv=3)

      IF (.NOT. allocated(base%cte_Transfo) ) THEN
        CALL alloc_NParray(base%cte_Transfo, (/ 20,1 /),'base%cte_Transfo','transfo_Q_TO_tQ')
      END IF

      DO i=1,get_nq_FROM_basis(base)
        q  = base%x(1,i)
        CALL sub_dntf(-itype,dntQ_inv,q,base%cte_Transfo(:,1),dnErr)
        IF (dnErr /= 0) THEN
          write(out_unitp,*) ' ERROR in transfo_Q_TO_tQ'
          write(out_unitp,*) '   ERROR in the sub_dntf call for the coordinates, iQdyn:',base%iQdyn(1)
          STOP 'ERROR in sub_dntf called from transfo_Q_TO_tQ'
        END IF

        base%x(1,i) = dntQ_inv%d0

        CALL sub_dntf(itype,dntQ,dntQ_inv%d0,base%cte_Transfo(:,1),dnErr)
        IF (dnErr /= 0) THEN
          write(out_unitp,*) ' ERROR in transfo_Q_TO_tQ'
          write(out_unitp,*) '   ERROR in the sub_dntf call for the coordinates, iQdyn:',base%iQdyn(1)
          STOP 'ERROR in sub_dntf called from transfo_Q_TO_tQ'
        END IF

        d1 = dntQ%d1(1)
        d2 = dntQ%d2(1,1)


        base%w(i)   = base%w(i)*base%rho(i)/d1
        base%rho(i) = d1
!       the wrho(:) is unchanged  !!
!       but nrho(1) has to be changed
        base%nrho(1) = itype  ! it means, the variable, x, is substituted by cos(theta)

        base%dnRGB%d2(i,:,1,1) = d2 * base%dnRGB%d1(i,:,1) + d1*d1 * base%dnRGB%d2(i,:,1,1)
        base%dnRGB%d1(i,:,1)   = d1 * base%dnRGB%d1(i,:,1)
      END DO

      CALL dealloc_dnS(dntQ)
      CALL dealloc_dnS(dntQ_inv)


      END SUBROUTINE transfo_Q_TO_tQ
