!===========================================================================
!===========================================================================
!This file is part of Tnum-Tana.
!
!    Tnum-Tana is a free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    Tnum-Tana is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ElVibRot.  If not, see <http://www.gnu.org/licenses/>.
!
!    Copyright 2015  David Lauvergnat
!      with contributions of Mamadou Ndong
!
!===========================================================================
!===========================================================================
      MODULE mod_dnM
      use mod_system, only: rkind, alloc_array, zero, out_unitp, &
                            dealloc_array, write_error_not_null, &
                            sub_test_tab_ub, sub_test_tab_lb,    &
                            error_memo_allo, write_error_null,   &
                            czero, write_vecmat, one, two
      use mod_dnS, only: alloc_array, dealloc_array, type_dns,   &
                         check_alloc_dns, alloc_dns, write_dns
      use mod_dnV, only: alloc_array, dealloc_array, type_dnvec, &
                         check_alloc_dnvec, alloc_dnvec
      IMPLICIT NONE

      PRIVATE

      TYPE Type_dnMat
          logical                     :: alloc=.FALSE.

          integer                     :: nderiv       = 0
          integer                     :: nb_var_deriv = 0
          integer                     :: nb_var_Matl  = 0
          integer                     :: nb_var_Matc  = 0

          real (kind=Rkind), pointer  :: d0(:,:)      => null()
          real (kind=Rkind), pointer  :: d1(:,:,:)    => null()
          real (kind=Rkind), pointer  :: d2(:,:,:,:)  => null()
          real (kind=Rkind), pointer  :: d3(:,:,:,:,:)=> null()

      END TYPE Type_dnMat

      TYPE Type_dnCplxMat
          logical                     :: alloc=.FALSE.

          integer                     :: nderiv       = 0
          integer                     :: nb_var_deriv = 0
          integer                     :: nb_var_Matl  = 0
          integer                     :: nb_var_Matc  = 0

          complex (kind=Rkind), pointer  :: d0(:,:)      => null()
          complex (kind=Rkind), pointer  :: d1(:,:,:)    => null()
          complex (kind=Rkind), pointer  :: d2(:,:,:,:)  => null()
          complex (kind=Rkind), pointer  :: d3(:,:,:,:,:)=> null()

      END TYPE Type_dnCplxMat

        INTERFACE assignment (=)
          MODULE PROCEDURE sub_dnMat2_TO_dnMat1,sub_dnCplxMat2_TO_dnCplxMat1
        END INTERFACE

      INTERFACE alloc_array
        MODULE PROCEDURE alloc_array_OF_dnMatdim1
        MODULE PROCEDURE alloc_array_OF_dnCplxMatdim1
      END INTERFACE
      INTERFACE dealloc_array
        MODULE PROCEDURE dealloc_array_OF_dnMatdim1
        MODULE PROCEDURE dealloc_array_OF_dnCplxMatdim1
      END INTERFACE

      PUBLIC :: Type_dnMat, alloc_dnMat, dealloc_dnMat, check_alloc_dnMat, Write_dnMat
      PUBLIC :: Type_dnCplxMat, alloc_dnCplxMat, dealloc_dnCplxMat, check_alloc_dnCplxMat, Write_dnCplxMat

      PUBLIC :: assignment (=), alloc_array, dealloc_array
      PUBLIC :: sub_dnMat1_TO_dnMat2, sub_dnMat1_TO_LargerdnMat2, sub_dnMat1_TO_dnMat2_partial
      PUBLIC :: dnVec_TO_dnMat, sub_dnMat_TO_dnS, sub_dnS_TO_dnMat
      PUBLIC :: sub_ZERO_TO_dnMat, sub_ZERO_TO_dnCplxMat
      PUBLIC :: dnVec1_wPLUS_dnMat2_TO_dnMat3,dnMat1_PLUS_dnMat2_TO_dnMat3
      PUBLIC :: dnMat1_MUL_dnMat2_TO_dnMat3, dnVec1_MUL_dnMat2_TO_dnVec3,dnMat1_MUL_dnVec2_TO_dnVec3
      PUBLIC :: TRANS_dnMat1_TO_dnMat2,INV_dnMat1_TO_dnMat2, Det_OF_dnMat_TO_dnS

      CONTAINS
!
!================================================================
!
!     allocation
!
!================================================================


      SUBROUTINE alloc_dnMat(dnMat,                                     &
                           nb_var_Matl,nb_var_Matc,nb_var_deriv,nderiv)
        TYPE (Type_dnMat) :: dnMat
        integer, optional :: nb_var_Matl,nb_var_Matc,nb_var_deriv,nderiv
        integer :: nd,nml,nmc
        integer :: err_mem

        IF (present(nderiv)) dnMat%nderiv = nderiv
        IF (present(nb_var_deriv)) dnMat%nb_var_deriv = nb_var_deriv
        IF (present(nb_var_Matl)) dnMat%nb_var_Matl = nb_var_Matl
        IF (present(nb_var_Matc)) dnMat%nb_var_Matc = nb_var_Matc

        IF (dnMat%nb_var_deriv == 0) dnMat%nderiv = 0

        nd = dnMat%nb_var_deriv
        nml = dnMat%nb_var_Matl
        nmc = dnMat%nb_var_Matc

        IF (dnMat%alloc) RETURN
        dnMat%alloc = .TRUE.


        IF (nml > 0 .AND. nmc > 0) THEN
          CALL alloc_array(dnMat%d0,(/ nml,nmc /),'dnMat%d0','alloc_dnMat')
          dnMat%d0(:,:) = ZERO

          IF (dnMat%nderiv >= 1) THEN
            CALL alloc_array(dnMat%d1,(/ nml,nmc,nd /),'dnMat%d1','alloc_dnMat')
            dnMat%d1(:,:,:) = ZERO
          END IF

          IF (dnMat%nderiv >= 2) THEN
            CALL alloc_array(dnMat%d2,(/ nml,nmc,nd,nd /),'dnMat%d2','alloc_dnMat')
            dnMat%d2(:,:,:,:) = ZERO
          END IF

          IF (dnMat%nderiv >= 3) THEN
            CALL alloc_array(dnMat%d3,(/ nml,nmc,nd,nd,nd /),'dnMat%d3','alloc_dnMat')
            dnMat%d3(:,:,:,:,:) = ZERO
          END IF

          IF (dnMat%nderiv > 3) THEN
            write(out_unitp,*) ' ERROR in alloc_dnMat'
            write(out_unitp,*) ' nderiv MUST be < 4',dnMat%nderiv
            STOP
          END IF
        ELSE
          write(out_unitp,*) ' ERROR in alloc_dnMat'
          !write(out_unitp,*) ' nb_var_deriv MUST be > 0',nd
          write(out_unitp,*) ' AND nb_var_matl MUST be > 0',nml
          write(out_unitp,*) ' AND nb_var_matc MUST be > 0',nmc
          STOP
        END IF

      END SUBROUTINE alloc_dnMat

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE dealloc_dnMat(dnMat)
        TYPE (Type_dnMat) :: dnMat
        integer :: nd,nml,nmc
        integer :: err_mem

        nd = dnMat%nb_var_deriv
        nml = dnMat%nb_var_Matl
        nmc = dnMat%nb_var_Matc

        IF (associated(dnMat%d0)) THEN
          CALL dealloc_array(dnMat%d0,'dnMat%d0','dealloc_dnMat')
        END IF

        IF (associated(dnMat%d1)) THEN
          CALL dealloc_array(dnMat%d1,'dnMat%d1','dealloc_dnMat')
        END IF

        IF (associated(dnMat%d2)) THEN
          CALL dealloc_array(dnMat%d2,'dnMat%d2','dealloc_dnMat')
        END IF

        IF (associated(dnMat%d3)) THEN
          CALL dealloc_array(dnMat%d3,'dnMat%d3','dealloc_dnMat')
        END IF

        dnMat%alloc    = .FALSE.

        dnMat%nderiv       = 0
        dnMat%nb_var_deriv = 0
        dnMat%nb_var_Matl  = 0
        dnMat%nb_var_Matc  = 0

      END SUBROUTINE dealloc_dnMat

      SUBROUTINE alloc_array_OF_dnMatdim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      TYPE (Type_dnMat), pointer, intent(out) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=1
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_dnMatdim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (associated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)

         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1)),stat=err_mem)
       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Type_dnMat')

      END SUBROUTINE alloc_array_OF_dnMatdim1
      SUBROUTINE dealloc_array_OF_dnMatdim1(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (Type_dnMat), pointer, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_dnMatdim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN
       IF (.NOT. associated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Type_dnMat')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_dnMatdim1

      SUBROUTINE alloc_dnCplxMat(dnMat,                                 &
                           nb_var_Matl,nb_var_Matc,nb_var_deriv,nderiv)
        TYPE (Type_dnCplxMat) :: dnMat
        integer, optional :: nb_var_Matl,nb_var_Matc,nb_var_deriv,nderiv
        integer :: nd,nml,nmc
        integer :: err_mem

        IF (present(nderiv)) dnMat%nderiv = nderiv
        IF (present(nb_var_deriv)) dnMat%nb_var_deriv = nb_var_deriv
        IF (present(nb_var_Matl)) dnMat%nb_var_Matl = nb_var_Matl
        IF (present(nb_var_Matc)) dnMat%nb_var_Matc = nb_var_Matc

        IF (dnMat%nb_var_deriv == 0) dnMat%nderiv = 0

        nd = dnMat%nb_var_deriv
        nml = dnMat%nb_var_Matl
        nmc = dnMat%nb_var_Matc

        IF (dnMat%alloc) RETURN
        dnMat%alloc = .TRUE.


        IF (nml > 0 .AND. nmc > 0) THEN
          CALL alloc_array(dnMat%d0,(/ nml,nmc /),'dnMat%d0','alloc_dnCplxMat')
          dnMat%d0(:,:) = CZERO

          IF (dnMat%nderiv >= 1) THEN
            CALL alloc_array(dnMat%d1,(/ nml,nmc,nd /),'dnMat%d1','alloc_dnCplxMat')
            dnMat%d1(:,:,:) = CZERO
          END IF

          IF (dnMat%nderiv >= 2) THEN
            CALL alloc_array(dnMat%d2,(/ nml,nmc,nd,nd /),'dnMat%d2','alloc_dnCplxMat')
            dnMat%d2(:,:,:,:) = CZERO
          END IF

          IF (dnMat%nderiv >= 3) THEN
            CALL alloc_array(dnMat%d3,(/ nml,nmc,nd,nd,nd /),'dnMat%d3','alloc_dnCplxMat')
            dnMat%d3(:,:,:,:,:) = CZERO
          END IF

          IF (dnMat%nderiv > 3) THEN
            write(out_unitp,*) ' ERROR in alloc_dnCplxMat'
            write(out_unitp,*) ' nderiv MUST be < 4',dnMat%nderiv
            STOP
          END IF
        ELSE
          write(out_unitp,*) ' ERROR in alloc_dnCplxMat'
          write(out_unitp,*) ' AND nb_var_matl MUST be > 0',nml
          write(out_unitp,*) ' AND nb_var_matc MUST be > 0',nmc
          STOP
        END IF

      END SUBROUTINE alloc_dnCplxMat

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE dealloc_dnCplxMat(dnMat)
        TYPE (Type_dnCplxMat) :: dnMat
        integer :: nd,nml,nmc
        integer :: err_mem

        nd = dnMat%nb_var_deriv
        nml = dnMat%nb_var_Matl
        nmc = dnMat%nb_var_Matc

        IF (associated(dnMat%d0)) THEN
          CALL dealloc_array(dnMat%d0,'dnMat%d0','dealloc_dnCplxMat')
        END IF

        IF (associated(dnMat%d1)) THEN
          CALL dealloc_array(dnMat%d1,'dnMat%d1','dealloc_dnCplxMat')
        END IF

        IF (associated(dnMat%d2)) THEN
          CALL dealloc_array(dnMat%d2,'dnMat%d2','dealloc_dnCplxMat')
        END IF

        IF (associated(dnMat%d3)) THEN
          CALL dealloc_array(dnMat%d3,'dnMat%d3','dealloc_dnCplxMat')
        END IF

        dnMat%alloc    = .FALSE.

        dnMat%nderiv       = 0
        dnMat%nb_var_deriv = 0
        dnMat%nb_var_Matl  = 0
        dnMat%nb_var_Matc  = 0

      END SUBROUTINE dealloc_dnCplxMat

      SUBROUTINE alloc_array_OF_dnCplxMatdim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      TYPE (Type_dnCplxMat), pointer, intent(out) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=1
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_dnCplxMatdim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (associated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)

         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1)),stat=err_mem)
       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Type_dnCplxMat')

      END SUBROUTINE alloc_array_OF_dnCplxMatdim1
      SUBROUTINE dealloc_array_OF_dnCplxMatdim1(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (Type_dnCplxMat), pointer, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_dnCplxMatdim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN
       IF (.NOT. associated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Type_dnCplxMat')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_dnCplxMatdim1


!================================================================
!
!     check if alloc has been done
!
!================================================================

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE check_alloc_dnMat(A,name_A,name_sub)
        TYPE (Type_dnMat), intent(in) :: A
        character (len=*), intent(in) :: name_A
        character (len=*), intent(in) :: name_sub

        IF ( .NOT. A%alloc) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) name_A,' has NOT been allocated with "alloc_dnMat"'
          write(out_unitp,*) ' CHECK the source!!!!!'
          STOP
        END IF
      END SUBROUTINE check_alloc_dnMat

      SUBROUTINE check_alloc_dnCplxMat(A,name_A,name_sub)
        TYPE (Type_dnCplxMat), intent(in) :: A
        character (len=*), intent(in) :: name_A
        character (len=*), intent(in) :: name_sub

        IF ( .NOT. A%alloc) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) name_A,' has NOT been allocated with "alloc_dnCplxMat"'
          write(out_unitp,*) ' CHECK the source!!!!!'
          STOP
        END IF
      END SUBROUTINE check_alloc_dnCplxMat

!================================================================
!        write the derived type
!================================================================

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE Write_dnMat(dnMat,nderiv)
        TYPE (Type_dnMat) :: dnMat
        integer, optional :: nderiv
        integer :: nderiv_loc
        integer :: i,j,k
        integer :: nl,nc

        IF (.NOT. dnMat%alloc) THEN
          write(out_unitp,*) 'BEGINNING Write_dnMat'
          write(out_unitp,*) 'dnMat is not allocated',dnMat%alloc
          write(out_unitp,*) 'END Write_dnMat'
          RETURN
        END IF
        !CALL check_alloc_dnMat(dnMat,'dnMat','Write_dnMat')

        nderiv_loc = dnMat%nderiv
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        write(out_unitp,*) 'BEGINNING Write_dnMat'
        write(out_unitp,*) 'nderiv,nb_var_deriv',dnMat%nderiv,dnMat%nb_var_deriv

        nl = dnMat%nb_var_Matl
        nc = dnMat%nb_var_Matc

        IF (nderiv_loc >= 0 .AND. associated(dnMat%d0)) THEN
          write(out_unitp,*) 'd0'
          CALL Write_VecMat(dnMat%d0,out_unitp,5)
        END IF
        IF (nderiv_loc > 0 .AND. associated(dnMat%d1)) THEN
          DO i=1,dnMat%nb_var_deriv
            write(out_unitp,*) 'd1',i
            CALL Write_VecMat(dnMat%d1(:,:,i),out_unitp,5)
          END DO
        END IF
        IF (nderiv_loc > 1 .AND. associated(dnMat%d2)) THEN
          DO i=1,dnMat%nb_var_deriv
          DO j=i,dnMat%nb_var_deriv
            write(out_unitp,*) 'd2',i,j
            CALL Write_VecMat(dnMat%d2(:,:,i,j),out_unitp,5)
          END DO
          END DO
        END IF
        IF (nderiv_loc > 2 .AND. associated(dnMat%d3)) THEN
          DO i=1,dnMat%nb_var_deriv
          DO j=i,dnMat%nb_var_deriv
          DO k=j,dnMat%nb_var_deriv
            write(out_unitp,*) 'd3',i,j,k
            CALL Write_VecMat(dnMat%d3(:,:,i,j,k),out_unitp,5)
          END DO
          END DO
          END DO
        END IF

        write(out_unitp,*) 'END Write_dnMat'


      END SUBROUTINE Write_dnMat

      SUBROUTINE Write_dnCplxMat(dnMat,nderiv)
        TYPE (Type_dnCplxMat) :: dnMat
        integer, optional :: nderiv
        integer :: nderiv_loc
        integer :: i,j,k
        integer :: nl,nc

        IF (.NOT. dnMat%alloc) THEN
          write(out_unitp,*) 'BEGINNING Write_dnCplxMat'
          write(out_unitp,*) 'dnMat is not allocated',dnMat%alloc
          write(out_unitp,*) 'END Write_dnCplxMat'
          RETURN
        END IF
        !CALL check_alloc_dnCplxMat(dnMat,'dnMat','Write_dnCplxMat')

        nderiv_loc = dnMat%nderiv
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        write(out_unitp,*) 'BEGINNING Write_dnCplxMat'
        write(out_unitp,*) 'nderiv,nb_var_deriv',dnMat%nderiv,dnMat%nb_var_deriv

        nl = dnMat%nb_var_Matl
        nc = dnMat%nb_var_Matc

        IF (nderiv_loc >= 0 .AND. associated(dnMat%d0)) THEN
          write(out_unitp,*) 'd0'
          CALL Write_VecMat(dnMat%d0,out_unitp,5)
        END IF
        IF (nderiv_loc > 0 .AND. associated(dnMat%d1)) THEN
          DO i=1,dnMat%nb_var_deriv
            write(out_unitp,*) 'd1',i
            CALL Write_VecMat(dnMat%d1(:,:,i),out_unitp,5)
          END DO
        END IF
        IF (nderiv_loc > 1 .AND. associated(dnMat%d2)) THEN
          DO i=1,dnMat%nb_var_deriv
          DO j=i,dnMat%nb_var_deriv
            write(out_unitp,*) 'd2',i,j
            CALL Write_VecMat(dnMat%d2(:,:,i,j),out_unitp,5)
          END DO
          END DO
        END IF
        IF (nderiv_loc > 2 .AND. associated(dnMat%d3)) THEN
          DO i=1,dnMat%nb_var_deriv
          DO j=i,dnMat%nb_var_deriv
          DO k=j,dnMat%nb_var_deriv
            write(out_unitp,*) 'd3',i,j,k
            CALL Write_VecMat(dnMat%d3(:,:,i,j,k),out_unitp,5)
          END DO
          END DO
          END DO
        END IF

        write(out_unitp,*) 'END Write_dnCplxMat'


      END SUBROUTINE Write_dnCplxMat

!================================================================
!        dnS2 = dnS1 , dnVec2 = dnVec1 ...
!        transfer Vec(iVec) => R or R => Vec(iVec)
!================================================================
      SUBROUTINE sub_dnCplxMat2_TO_dnCplxMat1(dnMat1,dnMat2)
        TYPE (Type_dnCplxMat), intent(inout) :: dnMat1
        TYPE (Type_dnCplxMat), intent(in)    :: dnMat2

        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='sub_dnCplxMat2_TO_dnCplxMat1'

        CALL check_alloc_dnCplxMat(dnMat2,'dnMat2',name_sub)
        IF (dnMat1%alloc) THEN
          CALL dealloc_dnCplxMat(dnMat1)
        END IF
        CALL alloc_dnCplxMat(dnMat1,dnMat2%nb_var_Matl,dnMat2%nb_var_Matc,&
                                      dnMat2%nb_var_deriv,dnMat2%nderiv)

        nderiv_loc = min(dnMat1%nderiv,dnMat2%nderiv)

        IF (dnMat2%nb_var_deriv /= dnMat1%nb_var_deriv) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nb_var_deriv in dnMat1 and dnMat2 are different!',&
                    dnMat1%nb_var_deriv,dnMat2%nb_var_deriv
          STOP
        END IF
        IF (dnMat1%nb_var_Matl /= dnMat2%nb_var_Matl) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nb_var_Matl in dnMat1 and dnMat2 are different!',&
                    dnMat1%nb_var_Matl,dnMat2%nb_var_Matl
          STOP
        END IF
        IF (dnMat1%nb_var_Matc /= dnMat2%nb_var_Matc) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nb_var_Matc in dnMat1 and dnMat2 are different!',&
                    dnMat1%nb_var_Matc,dnMat2%nb_var_Matc
          STOP
        END IF

        IF (nderiv_loc == 0) THEN
           dnMat1%d0 = dnMat2%d0
        ELSE IF (nderiv_loc == 1) THEN
           dnMat1%d0 = dnMat2%d0
           dnMat1%d1 = dnMat2%d1
        ELSE IF (nderiv_loc == 2) THEN
           dnMat1%d0 = dnMat2%d0
           dnMat1%d1 = dnMat2%d1
           dnMat1%d2 = dnMat2%d2
        ELSE IF (nderiv_loc == 3) THEN
           dnMat1%d0 = dnMat2%d0
           dnMat1%d1 = dnMat2%d1
           dnMat1%d2 = dnMat2%d2
           dnMat1%d3 = dnMat2%d3
        ELSE
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nderiv > 4 is NOT possible',nderiv_loc
          write(out_unitp,*) 'It souhld never append! Check the source'
          STOP
        END IF
      END SUBROUTINE sub_dnCplxMat2_TO_dnCplxMat1
      SUBROUTINE sub_dnMat2_TO_dnMat1(dnMat1,dnMat2)
        TYPE (Type_dnMat), intent(inout) :: dnMat1
        TYPE (Type_dnMat), intent(in)    :: dnMat2

        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='sub_dnMat2_TO_dnMat1'

        CALL check_alloc_dnMat(dnMat2,'dnMat2',name_sub)
        IF (dnMat1%alloc) THEN
          CALL dealloc_dnMat(dnMat1)
        END IF
        CALL alloc_dnMat(dnMat1,dnMat2%nb_var_Matl,dnMat2%nb_var_Matc,  &
                                      dnMat2%nb_var_deriv,dnMat2%nderiv)

        nderiv_loc = min(dnMat1%nderiv,dnMat2%nderiv)

        IF (dnMat2%nb_var_deriv /= dnMat1%nb_var_deriv) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nb_var_deriv in dnMat1 and dnMat2 are different!',&
                    dnMat1%nb_var_deriv,dnMat2%nb_var_deriv
          STOP
        END IF
        IF (dnMat1%nb_var_Matl /= dnMat2%nb_var_Matl) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nb_var_Matl in dnMat1 and dnMat2 are different!',&
                    dnMat1%nb_var_Matl,dnMat2%nb_var_Matl
          STOP
        END IF
        IF (dnMat1%nb_var_Matc /= dnMat2%nb_var_Matc) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nb_var_Matc in dnMat1 and dnMat2 are different!',&
                    dnMat1%nb_var_Matc,dnMat2%nb_var_Matc
          STOP
        END IF

        IF (nderiv_loc == 0) THEN
           dnMat1%d0 = dnMat2%d0
        ELSE IF (nderiv_loc == 1) THEN
           dnMat1%d0 = dnMat2%d0
           dnMat1%d1 = dnMat2%d1
        ELSE IF (nderiv_loc == 2) THEN
           dnMat1%d0 = dnMat2%d0
           dnMat1%d1 = dnMat2%d1
           dnMat1%d2 = dnMat2%d2
        ELSE IF (nderiv_loc == 3) THEN
           dnMat1%d0 = dnMat2%d0
           dnMat1%d1 = dnMat2%d1
           dnMat1%d2 = dnMat2%d2
           dnMat1%d3 = dnMat2%d3
        ELSE
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nderiv > 4 is NOT possible',nderiv_loc
          write(out_unitp,*) 'It souhld never append! Check the source'
          STOP
        END IF
      END SUBROUTINE sub_dnMat2_TO_dnMat1

      SUBROUTINE sub_dnMat1_TO_dnMat2(dnMat1,dnMat2,nderiv)
        TYPE (Type_dnMat) :: dnMat1,dnMat2
        integer, optional :: nderiv

        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='sub_dnMat1_TO_dnMat2'

        CALL check_alloc_dnMat(dnMat1,'dnMat1',name_sub)
        IF (.NOT. dnMat2%alloc) THEN
          CALL alloc_dnMat(dnMat2,dnMat1%nb_var_Matl,dnMat1%nb_var_Matc,&
                                      dnMat1%nb_var_deriv,dnMat1%nderiv)
        END IF

        nderiv_loc = min(dnMat1%nderiv,dnMat2%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (dnMat1%nb_var_deriv /= dnMat2%nb_var_deriv) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nb_var_deriv in dnMat1 and dnMat2 are different!',&
                    dnMat1%nb_var_deriv,dnMat2%nb_var_deriv
          STOP
        END IF
        IF (dnMat1%nb_var_Matl /= dnMat2%nb_var_Matl) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nb_var_Matl in dnMat1 and dnMat2 are different!',&
                    dnMat1%nb_var_Matl,dnMat2%nb_var_Matl
          STOP
        END IF
        IF (dnMat1%nb_var_Matc /= dnMat2%nb_var_Matc) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nb_var_Matc in dnMat1 and dnMat2 are different!',&
                    dnMat1%nb_var_Matc,dnMat2%nb_var_Matc
          STOP
        END IF

        IF (nderiv_loc == 0) THEN
           dnMat2%d0 = dnMat1%d0
        ELSE IF (nderiv_loc == 1) THEN
           dnMat2%d0 = dnMat1%d0
           dnMat2%d1 = dnMat1%d1
        ELSE IF (nderiv_loc == 2) THEN
           dnMat2%d0 = dnMat1%d0
           dnMat2%d1 = dnMat1%d1
           dnMat2%d2 = dnMat1%d2
        ELSE IF (nderiv_loc == 3) THEN
           dnMat2%d0 = dnMat1%d0
           dnMat2%d1 = dnMat1%d1
           dnMat2%d2 = dnMat1%d2
           dnMat2%d3 = dnMat1%d3
        ELSE
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nderiv > 4 is NOT possible',nderiv_loc
          write(out_unitp,*) 'It souhld never append! Check the source'
          STOP
        END IF
      END SUBROUTINE sub_dnMat1_TO_dnMat2

      SUBROUTINE sub_dnMat1_TO_LargerdnMat2(dnMat1,dnMat2,add_l,add_c,nderiv)
        TYPE (Type_dnMat) :: dnMat1,dnMat2
        integer :: add_l,add_c
        integer, optional :: nderiv

        integer :: nderiv_loc,l1,c1
        character (len=*), parameter :: name_sub='sub_dnMat1_TO_LargerdnMat2'

        CALL check_alloc_dnMat(dnMat1,'dnMat1',name_sub)
        IF (.NOT. dnMat2%alloc) THEN
          CALL alloc_dnMat(dnMat2,                                      &
                      dnMat1%nb_var_Matl+add_l,dnMat1%nb_var_Matc+add_c,&
                                      dnMat1%nb_var_deriv,dnMat1%nderiv)
        ELSE
          IF (dnMat1%nb_var_Matl+add_l /= dnMat2%nb_var_Matl .OR.       &
              dnMat1%nb_var_Matc+add_c /= dnMat2%nb_var_Matc ) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) '   Incompatible nb_var_Matl and add_l'
            write(out_unitp,*) '  dnMat1%nb_var_Matl,add_l',dnMat1%nb_var_Matl,add_l
            write(out_unitp,*) '  dnMat2%nb_var_Matl',dnMat2%nb_var_Matl
            write(out_unitp,*) ' OR Incompatible nb_var_Matc and add_c'
            write(out_unitp,*) '  dnMat1%nb_var_Matc,add_c',dnMat1%nb_var_Matc,add_c
            write(out_unitp,*) '  dnMat2%nb_var_Matc',dnMat2%nb_var_Matc
            STOP
          END IF
          IF (dnMat1%nb_var_deriv /= dnMat2%nb_var_deriv) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' nb_var_deriv in dnMat1 and dnMat2 are different!',&
                    dnMat1%nb_var_deriv,dnMat2%nb_var_deriv
            STOP
          END IF
        END IF

        nderiv_loc = min(dnMat1%nderiv,dnMat2%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        c1 = dnMat1%nb_var_Matc
        l1 = dnMat1%nb_var_Matl

        CALL sub_ZERO_TO_dnMat(dnMat2,nderiv_loc)

        IF (nderiv_loc == 0) THEN
           dnMat2%d0(1:l1,1:c1) = dnMat1%d0
        ELSE IF (nderiv_loc == 1) THEN
           dnMat2%d0(1:l1,1:c1)   = dnMat1%d0
           dnMat2%d1(1:l1,1:c1,:) = dnMat1%d1
        ELSE IF (nderiv_loc == 2) THEN
           dnMat2%d0(1:l1,1:c1)     = dnMat1%d0
           dnMat2%d1(1:l1,1:c1,:)   = dnMat1%d1
           dnMat2%d2(1:l1,1:c1,:,:) = dnMat1%d2
        ELSE IF (nderiv_loc == 3) THEN
           dnMat2%d0(1:l1,1:c1)       = dnMat1%d0
           dnMat2%d1(1:l1,1:c1,:)     = dnMat1%d1
           dnMat2%d2(1:l1,1:c1,:,:)   = dnMat1%d2
           dnMat2%d3(1:l1,1:c1,:,:,:) = dnMat1%d3
        ELSE
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nderiv > 4 is NOT possible',nderiv_loc
          write(out_unitp,*) 'It souhld never append! Check the source'
          STOP
        END IF
      END SUBROUTINE sub_dnMat1_TO_LargerdnMat2

      SUBROUTINE sub_dnMat1_TO_dnMat2_partial(dnMat1,dnMat2,nderiv)
        TYPE (Type_dnMat) :: dnMat1,dnMat2
        integer, optional :: nderiv

        integer :: nderiv_loc,nd
        character (len=*), parameter :: name_sub='sub_dnMat1_TO_dnMat2_partial'

        CALL check_alloc_dnMat(dnMat1,'dnMat1',name_sub)
        CALL check_alloc_dnMat(dnMat2,'dnMat2',name_sub)

        nderiv_loc = min(dnMat1%nderiv,dnMat2%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        nd = min(dnMat1%nb_var_deriv,dnMat2%nb_var_deriv)

        IF (dnMat1%nb_var_Matl /= dnMat2%nb_var_Matl) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nb_var_Matl in dnMat1 and dnMat2 are different!',&
                    dnMat1%nb_var_Matl,dnMat2%nb_var_Matl
          STOP
        END IF
        IF (dnMat1%nb_var_Matc /= dnMat2%nb_var_Matc) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nb_var_Matc in dnMat1 and dnMat2 are different!',&
                    dnMat1%nb_var_Matc,dnMat2%nb_var_Matc
          STOP
        END IF

        IF (nderiv_loc == 0) THEN
           dnMat2%d0                     = dnMat1%d0
        ELSE IF (nderiv_loc == 1) THEN
           dnMat2%d0                     = dnMat1%d0
           dnMat2%d1(:,:,1:nd)           = dnMat1%d1(:,:,1:nd)
        ELSE IF (nderiv_loc == 2) THEN
           dnMat2%d0                     = dnMat1%d0
           dnMat2%d1(:,:,1:nd)           = dnMat1%d1(:,:,1:nd)
           dnMat2%d2(:,:,1:nd,1:nd)      = dnMat1%d2(:,:,1:nd,1:nd)
        ELSE IF (nderiv_loc == 3) THEN
           dnMat2%d0                     = dnMat1%d0
           dnMat2%d1(:,:,1:nd)           = dnMat1%d1(:,:,1:nd)
           dnMat2%d2(:,:,1:nd,1:nd)      = dnMat1%d2(:,:,1:nd,1:nd)
           dnMat2%d3(:,:,1:nd,1:nd,1:nd) = dnMat1%d3(:,:,1:nd,1:nd,1:nd)
        ELSE
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nderiv > 4 is NOT possible',nderiv_loc
          write(out_unitp,*) 'It souhld never append! Check the source'
          STOP
        END IF
      END SUBROUTINE sub_dnMat1_TO_dnMat2_partial
      SUBROUTINE dnVec_TO_dnMat(dnVec,dnMat)
        TYPE (Type_dnMat), intent (inout) :: dnMat
        TYPE (Type_dnVec), intent (in)    :: dnVec


        character (len=*), parameter :: name_sub='dnVec_TO_dnMat'

        !write(out_unitp,*) 'In ',name_sub,': dnVec'
        !CALL Write_dnVec(dnVec)


        CALL check_alloc_dnVec(dnVec,'dnVec',name_sub)

        IF (.NOT. dnMat%alloc) THEN
          CALL alloc_dnMat(dnMat,                                       &
                           nb_var_Matl=dnVec%nb_var_vec,                &
                           nb_var_Matc=dnVec%nb_var_deriv,              &
                           nb_var_deriv=dnVec%nb_var_deriv,             &
                           nderiv=dnVec%nderiv-1)
        END IF

        IF (dnVec%nderiv == 0) RETURN

        IF (dnVec%nderiv > 0) THEN
          dnMat%d0 = dnVec%d1
        END IF

        IF (dnVec%nderiv > 1) THEN
          dnMat%d1 = dnVec%d2
        END IF

        IF (dnVec%nderiv > 2) THEN
          dnMat%d2 = dnVec%d3
        END IF

        !write(out_unitp,*) 'In ',name_sub,': dnMat'
        !CALL Write_dnMat(dnMat)

      END SUBROUTINE dnVec_TO_dnMat
      SUBROUTINE sub_dnMat_TO_dnS(dnMat,ilin,icol,dnS,nderiv)
        TYPE (Type_dnMat) :: dnMat
        integer           :: ilin,icol
        TYPE (Type_dnS)   :: dnS
        integer, optional :: nderiv

        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='sub_dnMat_TO_dnS'

        CALL check_alloc_dnMat(dnMat,'dnMat',name_sub)
        CALL check_alloc_dnS(dnS,'dnS',name_sub)

        nderiv_loc = min(dnMat%nderiv,dnS%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)


        IF (nderiv_loc == 0) THEN
           dnS%d0 = dnMat%d0(ilin,icol)
        ELSE IF (nderiv_loc == 1) THEN
           dnS%d0 = dnMat%d0(ilin,icol)
           dnS%d1 = dnMat%d1(ilin,icol,:)
        ELSE IF (nderiv_loc == 2) THEN
           dnS%d0 = dnMat%d0(ilin,icol)
           dnS%d1 = dnMat%d1(ilin,icol,:)
           dnS%d2 = dnMat%d2(ilin,icol,:,:)
        ELSE IF (nderiv_loc == 3) THEN
           dnS%d0 = dnMat%d0(ilin,icol)
           dnS%d1 = dnMat%d1(ilin,icol,:)
           dnS%d2 = dnMat%d2(ilin,icol,:,:)
           dnS%d3 = dnMat%d3(ilin,icol,:,:,:)
        ELSE
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nderiv > 4 is NOT possible',nderiv
          write(out_unitp,*) 'It souhld never append! Check the source'
          STOP
        END IF
      END SUBROUTINE sub_dnMat_TO_dnS
      SUBROUTINE sub_dnS_TO_dnMat(dnS,dnMat,ilin,icol,nderiv)
        TYPE (Type_dnMat) :: dnMat
        integer           :: ilin,icol
        TYPE (Type_dnS)   :: dnS
        integer, optional :: nderiv

        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='sub_dnS_TO_dnMat'

        CALL check_alloc_dnMat(dnMat,'dnMat',name_sub)
        CALL check_alloc_dnS(dnS,'dnS',name_sub)

        nderiv_loc = min(dnMat%nderiv,dnS%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)


        IF (nderiv_loc == 0) THEN
           dnMat%d0(ilin,icol) = dnS%d0
        ELSE IF (nderiv_loc == 1) THEN
           dnMat%d0(ilin,icol)   = dnS%d0
           dnMat%d1(ilin,icol,:) = dnS%d1
        ELSE IF (nderiv_loc == 2) THEN
           dnMat%d0(ilin,icol)     = dnS%d0
           dnMat%d1(ilin,icol,:)   = dnS%d1
           dnMat%d2(ilin,icol,:,:) = dnS%d2
        ELSE IF (nderiv_loc == 3) THEN
           dnMat%d0(ilin,icol)       = dnS%d0
           dnMat%d1(ilin,icol,:)     = dnS%d1
           dnMat%d2(ilin,icol,:,:)   = dnS%d2
           dnMat%d3(ilin,icol,:,:,:) = dnS%d3
        ELSE
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nderiv > 4 is NOT possible',nderiv
          write(out_unitp,*) 'It souhld never append! Check the source'
          STOP
        END IF
      END SUBROUTINE sub_dnS_TO_dnMat
      SUBROUTINE dnVec1_wPLUS_dnMat2_TO_dnMat3(dnVec1,iVec1,w1,         &
                                               dnMat2,ilin2,icol2,w2,   &
                                               dnMat3,ilin3,icol3,nderiv)
        TYPE (Type_dnMat)   :: dnMat2,dnMat3
        TYPE (Type_dnVec)   :: dnVec1
        integer, intent(in) :: iVec1,ilin2,icol2,ilin3,icol3
        integer, optional   :: nderiv
        real (kind=Rkind)   :: w1,w2


        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='dnVec1_wPLUS_dnMat2_TO_dnMat3'

        CALL check_alloc_dnMat(dnMat2,'dnMat2',name_sub)
        CALL check_alloc_dnMat(dnMat3,'dnMat3',name_sub)
        CALL check_alloc_dnVec(dnVec1,'dnVec1',name_sub)

        nderiv_loc = min(dnMat2%nderiv,dnMat3%nderiv,dnVec1%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (dnMat2%nb_var_deriv /= dnMat3%nb_var_deriv) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nb_var_deriv in dnMat2 and dnMat3 are different!',&
                    dnMat2%nb_var_deriv,dnMat3%nb_var_deriv
          STOP
        END IF
        IF (dnMat2%nb_var_deriv /= dnVec1%nb_var_deriv) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nb_var_deriv in dnMat2 and dnVec1 are different!',&
                    dnMat2%nb_var_deriv,dnVec1%nb_var_deriv
          STOP
        END IF

        IF (iVec1 <0 .OR. iVec1 > dnVec1%nb_var_vec) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' iVec1 is not compatible dnVec1',          &
                                                iVec1,dnVec1%nb_var_vec
          STOP
        END IF

        IF (ilin2 <0 .OR. ilin2 > dnMat2%nb_var_Matl .OR.               &
            icol2 <0 .OR. icol2 > dnMat2%nb_var_Matc) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' ilin2 is not compatible dnMat2',          &
                                                ilin2,dnMat2%nb_var_Matl
         write(out_unitp,*) ' OR icol2 is not compatible dnMat2',       &
                                                icol2,dnMat2%nb_var_Matc
          STOP
        END IF

        IF (ilin3 <0 .OR. ilin3 > dnMat3%nb_var_Matl .OR.               &
            icol3 <0 .OR. icol3 > dnMat3%nb_var_Matc) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' ilin3 is not compatible dnMat3',          &
                                                ilin3,dnMat3%nb_var_Matl
         write(out_unitp,*) ' OR icol3 is not compatible dnMat3',       &
                                                icol3,dnMat3%nb_var_Matc
          STOP
        END IF

        IF (nderiv_loc >= 0) THEN
           dnMat3%d0(ilin3,icol3) = w2*dnMat2%d0(ilin2,icol2)+w1*dnVec1%d0(iVec1)
        END IF
        IF (nderiv_loc >= 1) THEN
           dnMat3%d1(ilin3,icol3,:) = w2*dnMat2%d1(ilin2,icol2,:)+w1*dnVec1%d1(iVec1,:)
        END IF
        IF (nderiv_loc >= 2) THEN
           dnMat3%d2(ilin3,icol3,:,:) = w2*dnMat2%d2(ilin2,icol2,:,:)+w1*dnVec1%d2(iVec1,:,:)
        END IF
        IF (nderiv_loc >= 3) THEN
           dnMat3%d3(ilin3,icol3,:,:,:) = w2*dnMat2%d3(ilin2,icol2,:,:,:)+w1*dnVec1%d3(iVec1,:,:,:)
        END IF

        IF (nderiv_loc >= 4) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nderiv > 4 is NOT possible',nderiv
          write(out_unitp,*) 'It souhld never append! Check the source'
          STOP
        END IF
      END SUBROUTINE dnVec1_wPLUS_dnMat2_TO_dnMat3

      SUBROUTINE dnMat1_PLUS_dnMat2_TO_dnMat3(dnMat1,dnMat2,dnMat3,w1,w2,nderiv)
        TYPE (Type_dnMat)            :: dnMat1,dnMat2,dnMat3
        real(kind=Rkind), optional   :: w1,w2
        integer,          optional   :: nderiv


        real(kind=Rkind)   :: w1_loc,w2_loc

        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='dnMat1_PLUS_dnMat2_TO_dnMat3'

        CALL check_alloc_dnMat(dnMat1,'dnMat1',name_sub)
        CALL check_alloc_dnMat(dnMat2,'dnMat2',name_sub)
        IF (.NOT. dnMat3%alloc) THEN
          CALL alloc_dnMat(dnMat3,dnMat1%nb_var_Matl,dnMat2%nb_var_Matc,&
                                      dnMat1%nb_var_deriv,dnMat1%nderiv)
        END IF

        nderiv_loc = min(dnMat2%nderiv,dnMat3%nderiv,dnMat1%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (dnMat2%nb_var_deriv /= dnMat3%nb_var_deriv) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nb_var_deriv in dnMat2 and dnMat3 are different!',&
                    dnMat2%nb_var_deriv,dnMat3%nb_var_deriv
          STOP
        END IF
        IF (dnMat2%nb_var_deriv /= dnMat1%nb_var_deriv) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nb_var_deriv in dnMat2 and dnMat1 are different!',&
                    dnMat2%nb_var_deriv,dnMat1%nb_var_deriv
          STOP
        END IF

        IF (present(w1)) THEN
          w1_loc = w1
        ELSE
          w1_loc = ONE
        END IF
        IF (present(w2)) THEN
          w2_loc = w2
        ELSE
          w2_loc = ONE
        END IF


          IF (nderiv_loc >= 0) THEN
            dnMat3%d0 = w1_loc*dnMat1%d0 + w2_loc*dnMat2%d0
          END IF
          IF (nderiv_loc >= 1) THEN
            dnMat3%d1 = w1_loc*dnMat1%d1 + w2_loc*dnMat2%d1
          END IF
          IF (nderiv_loc >= 2) THEN
            dnMat3%d2 = w1_loc*dnMat1%d2 + w2_loc*dnMat2%d2

          END IF
          IF (nderiv_loc >= 3) THEN
            dnMat3%d3 = w1_loc*dnMat1%d3 + w2_loc*dnMat2%d3
          END IF

        IF (nderiv_loc >= 4) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nderiv > 4 is NOT possible',nderiv
          write(out_unitp,*) 'It souhld never append! Check the source'
          STOP
        END IF
      END SUBROUTINE dnMat1_PLUS_dnMat2_TO_dnMat3

      SUBROUTINE dnMat1_MUL_dnMat2_TO_dnMat3(dnMat1,dnMat2,dnMat3,nderiv)
        TYPE (Type_dnMat)   :: dnMat1,dnMat2,dnMat3
        integer, optional   :: nderiv

        integer :: id,jd,kd


        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='dnMat1_MUL_dnMat2_TO_dnMat3'

        CALL check_alloc_dnMat(dnMat1,'dnMat1',name_sub)
        CALL check_alloc_dnMat(dnMat2,'dnMat2',name_sub)
        IF (.NOT. dnMat3%alloc) THEN
          CALL alloc_dnMat(dnMat3,dnMat1%nb_var_Matl,dnMat2%nb_var_Matc,&
                                      dnMat1%nb_var_deriv,dnMat1%nderiv)
        END IF

        nderiv_loc = min(dnMat2%nderiv,dnMat3%nderiv,dnMat1%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (dnMat2%nb_var_deriv /= dnMat3%nb_var_deriv) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nb_var_deriv in dnMat2 and dnMat3 are different!',&
                    dnMat2%nb_var_deriv,dnMat3%nb_var_deriv
          STOP
        END IF
        IF (dnMat2%nb_var_deriv /= dnMat1%nb_var_deriv) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nb_var_deriv in dnMat2 and dnMat1 are different!',&
                    dnMat2%nb_var_deriv,dnMat1%nb_var_deriv
          STOP
        END IF



          IF (nderiv_loc >= 0) THEN
            dnMat3%d0 = matmul(dnMat1%d0,dnMat2%d0)
          END IF
          IF (nderiv_loc >= 1) THEN
            DO id=1,dnMat1%nb_var_deriv
              dnMat3%d1(:,:,id) = matmul(dnMat1%d1(:,:,id),dnMat2%d0) + &
                                  matmul(dnMat1%d0,dnMat2%d1(:,:,id))
            END DO
          END IF
          IF (nderiv_loc >= 2) THEN
            DO id=1,dnMat1%nb_var_deriv
            DO jd=1,dnMat1%nb_var_deriv

              dnMat3%d2(:,:,id,jd) = matmul(dnMat1%d2(:,:,id,jd),dnMat2%d0) + &
                                matmul(dnMat1%d1(:,:,id),dnMat2%d1(:,:,jd)) + &
                                matmul(dnMat1%d1(:,:,jd),dnMat2%d1(:,:,id)) + &
                                     matmul(dnMat1%d0,dnMat2%d2(:,:,id,jd))

            END DO
            END DO

          END IF
          IF (nderiv_loc >= 3) THEN
            DO id=1,dnMat1%nb_var_deriv
            DO jd=1,dnMat1%nb_var_deriv
            DO kd=1,dnMat1%nb_var_deriv

              dnMat3%d3(:,:,id,jd,kd) =                                 &
                       matmul(dnMat1%d3(:,:,id,jd,kd),dnMat2%d0(:,:)) + &

                       matmul(dnMat1%d2(:,:,id,jd),dnMat2%d1(:,:,kd)) + &
                       matmul(dnMat1%d2(:,:,id,kd),dnMat2%d1(:,:,jd)) + &
                       matmul(dnMat1%d2(:,:,jd,kd),dnMat2%d1(:,:,id)) + &

                       matmul(dnMat1%d1(:,:,id),dnMat2%d2(:,:,jd,kd)) + &
                       matmul(dnMat1%d1(:,:,jd),dnMat2%d2(:,:,id,kd)) + &
                       matmul(dnMat1%d1(:,:,kd),dnMat2%d2(:,:,id,jd)) + &

                       matmul(dnMat1%d0(:,:),dnMat2%d3(:,:,id,jd,kd))
            END DO
            END DO
            END DO
          END IF

        IF (nderiv_loc >= 4) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nderiv > 4 is NOT possible',nderiv
          write(out_unitp,*) 'It souhld never append! Check the source'
          STOP
        END IF
      END SUBROUTINE dnMat1_MUL_dnMat2_TO_dnMat3

     SUBROUTINE dnVec1_MUL_dnMat2_TO_dnVec3(dnVec1,dnMat2,dnVec3,nderiv)
        TYPE (Type_dnMat)   :: dnMat2
        TYPE (Type_dnVec)   :: dnVec1,dnVec3

        integer, optional   :: nderiv

        integer :: id,jd,kd


        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='dnVec1_MUL_dnMat2_TO_dnVec3'

        CALL check_alloc_dnVec(dnVec1,'dnVec1',name_sub)
        CALL check_alloc_dnMat(dnMat2,'dnMat2',name_sub)
        nderiv_loc = min(dnMat2%nderiv,dnVec1%nderiv)
        IF (.NOT. dnVec3%alloc) THEN
          CALL alloc_dnVec(dnVec3,dnMat2%nb_var_Matc,dnMat2%nb_var_deriv,nderiv_loc)
        END IF

        nderiv_loc = min(dnMat2%nderiv,dnVec3%nderiv,dnVec1%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (dnMat2%nb_var_deriv /= dnVec3%nb_var_deriv) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nb_var_deriv in dnMat2 and dnVec3 are different!',&
                    dnMat2%nb_var_deriv,dnVec3%nb_var_deriv
          STOP
        END IF
        IF (dnMat2%nb_var_deriv /= dnVec1%nb_var_deriv) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nb_var_deriv in dnMat2 and dnVec1 are different!',&
                    dnMat2%nb_var_deriv,dnVec1%nb_var_deriv
          STOP
        END IF



          IF (nderiv_loc >= 0) THEN
            dnVec3%d0 = matmul(dnVec1%d0,dnMat2%d0)
          END IF
          IF (nderiv_loc >= 1) THEN
            DO id=1,dnVec1%nb_var_deriv
              dnVec3%d1(:,id) = matmul(dnVec1%d1(:,id),dnMat2%d0) + &
                                  matmul(dnVec1%d0,dnMat2%d1(:,:,id))
            END DO
          END IF
          IF (nderiv_loc >= 2) THEN
            DO id=1,dnVec1%nb_var_deriv
            DO jd=1,dnVec1%nb_var_deriv

              dnVec3%d2(:,id,jd) = matmul(dnVec1%d2(:,id,jd),dnMat2%d0) + &
                                matmul(dnVec1%d1(:,id),dnMat2%d1(:,:,jd)) + &
                                matmul(dnVec1%d1(:,jd),dnMat2%d1(:,:,id)) + &
                                     matmul(dnVec1%d0,dnMat2%d2(:,:,id,jd))

            END DO
            END DO

          END IF
          IF (nderiv_loc >= 3) THEN
            DO id=1,dnVec1%nb_var_deriv
            DO jd=1,dnVec1%nb_var_deriv
            DO kd=1,dnVec1%nb_var_deriv

              dnVec3%d3(:,id,jd,kd) =                                 &
                       matmul(dnVec1%d3(:,id,jd,kd),dnMat2%d0(:,:)) + &

                       matmul(dnVec1%d2(:,id,jd),dnMat2%d1(:,:,kd)) + &
                       matmul(dnVec1%d2(:,id,kd),dnMat2%d1(:,:,jd)) + &
                       matmul(dnVec1%d2(:,jd,kd),dnMat2%d1(:,:,id)) + &

                       matmul(dnVec1%d1(:,id),dnMat2%d2(:,:,jd,kd)) + &
                       matmul(dnVec1%d1(:,jd),dnMat2%d2(:,:,id,kd)) + &
                       matmul(dnVec1%d1(:,kd),dnMat2%d2(:,:,id,jd)) + &

                       matmul(dnVec1%d0(:),dnMat2%d3(:,:,id,jd,kd))
            END DO
            END DO
            END DO
          END IF

        IF (nderiv_loc >= 4) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nderiv > 4 is NOT possible',nderiv
          write(out_unitp,*) 'It souhld never append! Check the source'
          STOP
        END IF
      END SUBROUTINE dnVec1_MUL_dnMat2_TO_dnVec3

      SUBROUTINE dnMat1_MUL_dnVec2_TO_dnVec3(dnMat1,dnVec2,dnVec3,nderiv)
        TYPE (Type_dnMat)   :: dnMat1
        TYPE (Type_dnVec)   :: dnVec2,dnVec3

        integer, optional   :: nderiv

        integer :: id,jd,kd


        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='dnMat1_MUL_dnVec2_TO_dnVec3'

        CALL check_alloc_dnMat(dnMat1,'dnMat1',name_sub)
        CALL check_alloc_dnvec(dnvec2,'dnVec2',name_sub)
        nderiv_loc = min(dnVec2%nderiv,dnMat1%nderiv)
        IF (.NOT. dnVec3%alloc) THEN
          CALL alloc_dnVec(dnVec3,dnMat1%nb_var_Matl,dnMat1%nb_var_deriv,nderiv_loc)
        END IF

        nderiv_loc = min(dnVec2%nderiv,dnVec3%nderiv,dnMat1%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (dnVec2%nb_var_deriv /= dnVec3%nb_var_deriv) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nb_var_deriv in dnVec2 and dnVec3 are different!',&
                    dnVec2%nb_var_deriv,dnVec3%nb_var_deriv
          STOP
        END IF
        IF (dnVec2%nb_var_deriv /= dnMat1%nb_var_deriv) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nb_var_deriv in dnVec2 and dnMat1 are different!',&
                    dnVec2%nb_var_deriv,dnMat1%nb_var_deriv
          STOP
        END IF



          IF (nderiv_loc >= 0) THEN
            dnVec3%d0 = matmul(dnMat1%d0,dnVec2%d0)
          END IF
          IF (nderiv_loc >= 1) THEN
            DO id=1,dnMat1%nb_var_deriv
              dnVec3%d1(:,id) = matmul(dnMat1%d1(:,:,id),dnVec2%d0) +   &
                                matmul(dnMat1%d0,dnVec2%d1(:,id))
            END DO
          END IF
          IF (nderiv_loc >= 2) THEN
            DO id=1,dnMat1%nb_var_deriv
            DO jd=1,dnMat1%nb_var_deriv

              dnVec3%d2(:,id,jd) = matmul(dnMat1%d2(:,:,id,jd),dnVec2%d0) + &
                                matmul(dnMat1%d1(:,:,id),dnVec2%d1(:,jd)) + &
                                matmul(dnMat1%d1(:,:,jd),dnVec2%d1(:,id)) + &
                                     matmul(dnMat1%d0,dnVec2%d2(:,id,jd))

            END DO
            END DO

          END IF
          IF (nderiv_loc >= 3) THEN
            DO id=1,dnMat1%nb_var_deriv
            DO jd=1,dnMat1%nb_var_deriv
            DO kd=1,dnMat1%nb_var_deriv

              dnVec3%d3(:,id,jd,kd) =                                 &
                       matmul(dnMat1%d3(:,:,id,jd,kd),dnVec2%d0(:)) + &

                       matmul(dnMat1%d2(:,:,id,jd),dnVec2%d1(:,kd)) + &
                       matmul(dnMat1%d2(:,:,id,kd),dnVec2%d1(:,jd)) + &
                       matmul(dnMat1%d2(:,:,jd,kd),dnVec2%d1(:,id)) + &

                       matmul(dnMat1%d1(:,:,id),dnVec2%d2(:,jd,kd)) + &
                       matmul(dnMat1%d1(:,:,jd),dnVec2%d2(:,id,kd)) + &
                       matmul(dnMat1%d1(:,:,kd),dnVec2%d2(:,id,jd)) + &

                       matmul(dnMat1%d0(:,:),dnVec2%d3(:,id,jd,kd))
            END DO
            END DO
            END DO
          END IF

        IF (nderiv_loc >= 4) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nderiv > 4 is NOT possible',nderiv
          write(out_unitp,*) 'It souhld never append! Check the source'
          STOP
        END IF
      END SUBROUTINE dnMat1_MUL_dnVec2_TO_dnVec3



      SUBROUTINE TRANS_dnMat1_TO_dnMat2(dnMat1,dnMat2,nderiv)
        TYPE (Type_dnMat)   :: dnMat1,dnMat2
        integer, optional   :: nderiv

        integer :: id,jd,kd


        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='TRANS_dnMat1_TO_dnMat2'

        CALL check_alloc_dnMat(dnMat1,'dnMat1',name_sub)
        IF (.NOT. dnMat2%alloc) THEN
          CALL alloc_dnMat(dnMat2,dnMat1%nb_var_Matc,dnMat1%nb_var_Matl,&
                                      dnMat1%nb_var_deriv,dnMat1%nderiv)
        END IF

        nderiv_loc = min(dnMat2%nderiv,dnMat1%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (dnMat2%nb_var_deriv /= dnMat1%nb_var_deriv) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nb_var_deriv in dnMat2 and dnMat1 are different!',&
                    dnMat2%nb_var_deriv,dnMat1%nb_var_deriv
          STOP
        END IF



          IF (nderiv_loc >= 0) THEN
            dnMat2%d0 = transpose(dnMat1%d0)
          END IF

          IF (nderiv_loc >= 1) THEN
            DO id=1,dnMat1%nb_var_deriv
              dnMat2%d1(:,:,id) = transpose(dnMat1%d1(:,:,id))
            END DO
          END IF

          IF (nderiv_loc >= 2) THEN
            DO id=1,dnMat1%nb_var_deriv
            DO jd=1,dnMat1%nb_var_deriv
              dnMat2%d2(:,:,id,jd) = transpose(dnMat1%d2(:,:,id,jd))
            END DO
            END DO

          END IF
          IF (nderiv_loc >= 3) THEN
            DO id=1,dnMat1%nb_var_deriv
            DO jd=1,dnMat1%nb_var_deriv
            DO kd=1,dnMat1%nb_var_deriv
              dnMat2%d3(:,:,id,jd,kd) = transpose(dnMat1%d3(:,:,id,jd,kd))
            END DO
            END DO
            END DO
          END IF

        IF (nderiv_loc >= 4) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nderiv > 4 is NOT possible',nderiv
          write(out_unitp,*) 'It souhld never append! Check the source'
          STOP
        END IF
      END SUBROUTINE TRANS_dnMat1_TO_dnMat2

      SUBROUTINE INV_dnMat1_TO_dnMat2(dnMat1,dnMat2,nderiv)
        TYPE (Type_dnMat)   :: dnMat1,dnMat2
        integer, optional   :: nderiv

        integer :: i,j,k


        integer :: nderiv_loc
        integer :: err_mem,memory
!----- for debuging --------------------------------------------------
       character (len=*), parameter :: name_sub='INV_dnMat1_TO_dnMat2'
       logical, parameter :: debug = .FALSE.
       !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) ' BEGINNING ',name_sub
         write(out_unitp,*)
         write(out_unitp,*) 'dnMat1'
         CALL Write_dnMat(dnMat1)
         flush(out_unitp)
       END IF
!-----------------------------------------------------------

        CALL check_alloc_dnMat(dnMat1,'dnMat1',name_sub)
        IF (.NOT. dnMat2%alloc) THEN
          CALL alloc_dnMat(dnMat2,dnMat1%nb_var_Matl,dnMat1%nb_var_Matc,&
                                      dnMat1%nb_var_deriv,dnMat1%nderiv)
        END IF

        nderiv_loc = min(dnMat2%nderiv,dnMat1%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (dnMat2%nb_var_deriv /= dnMat1%nb_var_deriv) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nb_var_deriv in dnMat2 and dnMat1 are different!',&
                    dnMat2%nb_var_deriv,dnMat1%nb_var_deriv
          STOP
        END IF


        CALL inv_m1_TO_m2(dnMat1%d0,dnMat2%d0,dnMat1%nb_var_Matl,0,ZERO) ! not SVD

        IF (nderiv_loc >= 1) THEN
          DO i=1,dnMat1%nb_var_deriv
            dnMat2%d1(:,:,i) =  -matmul(                                &
                                  matmul(dnMat2%d0,dnMat1%d1(:,:,i))   ,&
                                        dnMat2%d0)
          END DO
        END IF

        IF (nderiv_loc >= 2) THEN
          DO i=1,dnMat1%nb_var_deriv

            dnMat2%d2(:,:,i,i) = -matmul(                               &
                      TWO*matmul(dnMat2%d1(:,:,i),dnMat1%d1(:,:,i)) +   &
                          matmul(dnMat2%d0,dnMat1%d2(:,:,i,i))         ,&
                                          dnMat2%d0)

          END DO

          DO i=1,dnMat1%nb_var_deriv
          DO j=i+1,dnMat1%nb_var_deriv

            dnMat2%d2(:,:,i,j) = -matmul(                               &
                           matmul(dnMat2%d1(:,:,i),dnMat1%d1(:,:,j)) +  &
                           matmul(dnMat2%d1(:,:,j),dnMat1%d1(:,:,i)) +  &
                           matmul(dnMat2%d0(:,:),dnMat1%d2(:,:,i,j))   ,&
                                          dnMat2%d0)

            dnMat2%d2(:,:,j,i) = dnMat2%d2(:,:,i,j)

          END DO
          END DO
        END IF

        IF (nderiv_loc >= 3) THEN
          write(out_unitp,*) ' not yet nderiv=3',name_sub
          STOP
        END IF

        IF (nderiv_loc >= 4) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nderiv > 4 is NOT possible',nderiv
          write(out_unitp,*) 'It should never append! Check the source'
          STOP
        END IF

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'dnMat2 = dnMat1^-1'
         CALL Write_dnMat(dnMat2)
         write(out_unitp,*) ' END ',name_sub
         flush(out_unitp)
       END IF
!-----------------------------------------------------------

      END SUBROUTINE INV_dnMat1_TO_dnMat2

      SUBROUTINE Det_OF_dnMat_TO_dnS(dnMat,dnS,nderiv,dnMat_inv)
        TYPE (Type_dnMat), intent(in)            :: dnMat
        TYPE (Type_dnS), intent(inout)           :: dnS
        TYPE (Type_dnMat), intent(in), optional  :: dnMat_inv

        integer, optional   :: nderiv

        integer :: i,j,k,l

        TYPE (Type_dnMat)               :: dnMat_inv_loc
        real (kind=Rkind), allocatable  :: mat_temp(:,:)
        real (kind=Rkind)  :: trace

        integer :: nderiv_loc
        integer :: err_mem,memory
!----- for debuging --------------------------------------------------
       character (len=*), parameter :: name_sub='Det_OF_dnMat_TO_dnS'
       logical, parameter :: debug = .FALSE.
       !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) ' BEGINNING ',name_sub
         write(out_unitp,*)
         write(out_unitp,*) 'dnMat'
         CALL Write_dnMat(dnMat)
         IF (present(dnMat_inv)) THEN
           write(out_unitp,*) 'dnMat_inv'
           CALL Write_dnMat(dnMat_inv)
         END IF
         flush(out_unitp)
       END IF
!-----------------------------------------------------------

        CALL check_alloc_dnMat(dnMat,'dnMat',name_sub)
        IF (present(dnMat_inv)) CALL check_alloc_dnMat(dnMat_inv,'dnMat_inv',name_sub)

        IF (.NOT. dnS%alloc) THEN
          CALL alloc_dnS(dnS,dnMat%nb_var_deriv,dnMat%nderiv)
        END IF

        nderiv_loc = dnMat%nderiv
        IF (present(dnMat_inv)) nderiv_loc = min(nderiv_loc,dnMat_inv%nderiv)
        IF (present(nderiv))    nderiv_loc = min(nderiv_loc,nderiv)

        IF (present(dnMat_inv)) THEN
        IF (dnMat%nb_var_deriv /= dnMat_inv%nb_var_deriv) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nb_var_deriv in dnMat and dnMat_inv are different!',&
                    dnMat%nb_var_deriv,dnMat_inv%nb_var_deriv
          STOP
        END IF
        END IF

        !CALL alloc_NParray(mat_temp,shape(dnMat%d0),'mat_temp',name_sub)

        IF (present(dnMat_inv)) THEN
          CALL Det_OF_m1(dnMat%d0,dnS%d0,dnMat%nb_var_Matl)

          IF (nderiv_loc >= 1) THEN
            DO i=1,dnMat%nb_var_deriv
              dnS%d1(i) = ZERO
              DO k=1,dnMat%nb_var_Matl
              DO l=1,dnMat%nb_var_Matl
                dnS%d1(i) = dnS%d1(i) + dnMat_inv%d0(k,l)*dnMat%d1(l,k,i)
              END DO
              END DO

              !mat_temp = matmul(dnMat_inv%d0,dnMat%d1(:,:,i))
              !DO k=1,dnMat%nb_var_Matl
              !  dnS%d1(i) = dnS%d1(i) + mat_temp(k,k)
              !END DO
            END DO
          END IF

          IF (nderiv_loc >= 2) THEN
            DO i=1,dnMat%nb_var_deriv
            DO j=i,dnMat%nb_var_deriv
              dnS%d2(i,j) = dnS%d1(i) * dnS%d1(j)

              DO k=1,dnMat%nb_var_Matl
              DO l=1,dnMat%nb_var_Matl
                dnS%d2(i,j) = dnS%d2(i,j) + dnMat_inv%d0(k,l)*dnMat%d2(l,k,i,j) + &
                                            dnMat_inv%d1(k,l,j)*dnMat%d1(l,k,i)
              END DO
              END DO

              !mat_temp = matmul(dnMat_inv%d0,dnMat%d2(:,:,i,j)) +       &
              !           matmul(dnMat_inv%d1(:,:,j),dnMat%d1(:,:,i))
              !DO k=1,dnMat%nb_var_Matl
              !  dnS%d2(i,j) = dnS%d2(i,j) + mat_temp(k,k)
              !END DO

              dnS%d2(j,i) = dnS%d2(i,j)
            END DO
            END DO
            dnS%d2(:,:) = dnS%d2(:,:) * dnS%d0

          END IF

          IF (nderiv_loc >= 1) THEN
            dnS%d1(:) = dnS%d1(:) * dnS%d0
          END IF

        ELSE
          CALL INV_dnMat1_TO_dnMat2(dnMat,dnMat_inv_loc,nderiv_loc)
          CALL Det_OF_m1(dnMat%d0,dnS%d0,dnMat%nb_var_Matl)


          IF (nderiv_loc >= 1) THEN
            DO i=1,dnMat%nb_var_deriv
              dnS%d1(i) = ZERO
              DO k=1,dnMat%nb_var_Matl
              DO l=1,dnMat%nb_var_Matl
                dnS%d1(i) = dnS%d1(i) + dnMat_inv_loc%d0(k,l)*dnMat%d1(l,k,i)
              END DO
              END DO

              !mat_temp = matmul(dnMat_inv_loc%d0,dnMat%d1(:,:,i))
              !DO k=1,dnMat%nb_var_Matl
              !  dnS%d1(i) = dnS%d1(i) + mat_temp(k,k)
              !END DO
            END DO
          END IF


          IF (nderiv_loc >= 2) THEN
            DO i=1,dnMat%nb_var_deriv
            DO j=i,dnMat%nb_var_deriv
              dnS%d2(i,j) = dnS%d1(i) * dnS%d1(j)

              DO k=1,dnMat%nb_var_Matl
              DO l=1,dnMat%nb_var_Matl
                dnS%d2(i,j) = dnS%d2(i,j) + dnMat_inv_loc%d0(k,l)*dnMat%d2(l,k,i,j) + &
                                            dnMat_inv_loc%d1(k,l,j)*dnMat%d1(l,k,i)
              END DO
              END DO

              !mat_temp = matmul(dnMat_inv_loc%d0,dnMat%d2(:,:,i,j)) +   &
              !           matmul(dnMat_inv_loc%d1(:,:,j),dnMat%d1(:,:,i))
              !DO k=1,dnMat%nb_var_Matl
              !  dnS%d2(i,j) = dnS%d2(i,j) + mat_temp(k,k)
              !END DO

              dnS%d2(j,i) = dnS%d2(i,j)
            END DO
            END DO
            dnS%d2(:,:) = dnS%d2(:,:) * dnS%d0


          END IF

          IF (nderiv_loc >= 1) THEN
            dnS%d1(:) = dnS%d1(:) * dnS%d0
          END IF
        END IF
        !CALL dealloc_NParray(mat_temp,'mat_temp',name_sub)



!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'dnS = dndet(Mat)'
         CALL Write_dnS(dnS)
         write(out_unitp,*) ' END ',name_sub
         flush(out_unitp)
       END IF
!-----------------------------------------------------------

      END SUBROUTINE Det_OF_dnMat_TO_dnS

!
!================================================================
!
!     dnS = 0
!     dnVec = 0
!     dnMat = 0
!
!================================================================
!
      !!@description: TODO
      !!@param: TODO
      SUBROUTINE sub_ZERO_TO_dnMat(dnMat,nderiv)
        TYPE (Type_dnMat) :: dnMat
        integer, optional :: nderiv
        integer :: nderiv_loc

        CALL check_alloc_dnMat(dnMat,'dnMat','sub_ZERO_TO_dnMat')

        nderiv_loc = dnMat%nderiv
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

          dnMat%d0(:,:) = ZERO
          IF (nderiv_loc == 1) THEN
            dnMat%d1(:,:,:) = ZERO
          ELSE IF (nderiv_loc == 2) THEN
            dnMat%d1(:,:,:) = ZERO
            dnMat%d2(:,:,:,:) = ZERO
          ELSE IF (nderiv_loc == 3) THEN
            dnMat%d1(:,:,:) = ZERO
            dnMat%d2(:,:,:,:) = ZERO
            dnMat%d3(:,:,:,:,:) = ZERO
          ELSE IF (nderiv_loc > 3) THEN
            write(out_unitp,*) ' ERROR in sub_ZERO_TO_dnMat'
            write(out_unitp,*) ' nderiv_loc MUST be < 4',nderiv_loc
            STOP
          END IF

      END SUBROUTINE sub_ZERO_TO_dnMat

      SUBROUTINE sub_ZERO_TO_dnCplxMat(dnMat,nderiv)
        TYPE (Type_dnCplxMat) :: dnMat
        integer, optional :: nderiv
        integer :: nderiv_loc

        CALL check_alloc_dnCplxMat(dnMat,'dnMat','sub_ZERO_TO_dnCplxMat')

        nderiv_loc = dnMat%nderiv
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

          dnMat%d0(:,:) = CZERO
          IF (nderiv_loc == 1) THEN
            dnMat%d1(:,:,:) = CZERO
          ELSE IF (nderiv_loc == 2) THEN
            dnMat%d1(:,:,:) = CZERO
            dnMat%d2(:,:,:,:) = CZERO
          ELSE IF (nderiv_loc == 3) THEN
            dnMat%d1(:,:,:) = CZERO
            dnMat%d2(:,:,:,:) = CZERO
            dnMat%d3(:,:,:,:,:) = CZERO
          ELSE IF (nderiv_loc > 3) THEN
            write(out_unitp,*) ' ERROR in sub_ZERO_TO_dnCplxMat'
            write(out_unitp,*) ' nderiv_loc MUST be < 4',nderiv_loc
            STOP
          END IF

      END SUBROUTINE sub_ZERO_TO_dnCplxMat
      END MODULE mod_dnM

