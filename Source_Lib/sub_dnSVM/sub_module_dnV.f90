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

MODULE mod_dnV
      USE mod_system
      IMPLICIT NONE

      PRIVATE

      !==============================================
      !!@description: TODO
      !!@param: TODO
      TYPE Type_dnVec
          logical                     :: alloc=.FALSE.

          integer                     :: nderiv       = 0
          integer                     :: nb_var_deriv = 0
          integer                     :: nb_var_vec   = 0

          real (kind=Rkind), pointer  :: d0(:)        => null()
          real (kind=Rkind), pointer  :: d1(:,:)      => null()
          real (kind=Rkind), pointer  :: d2(:,:,:)    => null()
          real (kind=Rkind), pointer  :: d3(:,:,:,:)  => null()
      CONTAINS
        PROCEDURE, PRIVATE, PASS(dnVec1) :: sub_dnVec2_TO_dnVec1
        GENERIC,   PUBLIC  :: assignment(=) => sub_dnVec2_TO_dnVec1
      END TYPE Type_dnVec

      INTERFACE alloc_array
        MODULE PROCEDURE alloc_array_OF_dnVecdim1,alloc_array_OF_dnVecdim2
      END INTERFACE
      INTERFACE dealloc_array
        MODULE PROCEDURE dealloc_array_OF_dnVecdim1,dealloc_array_OF_dnVecdim2
      END INTERFACE

      PUBLIC :: Type_dnVec, alloc_array, dealloc_array
      PUBLIC :: alloc_dnVec, dealloc_dnVec, check_alloc_dnVec, Write_dnVec, sub_Normalize_dnVec

      PUBLIC :: get_nderiv_FROM_dnVec,get_nb_var_deriv_FROM_dnVec


      PUBLIC :: sub_dnVec_TO_dnS, sub_dnS_TO_dnVec, sub_dnVec1_TO_dnVec2
      PUBLIC :: sub_PartdnVec1_TO_PartdnVec2, sub_dnVec1_TO_dnVec2_WithIvec, sub_dnVec1_TO_dnVec2_partial
      PUBLIC :: sub_dnVec1_wTO_dnVec2

      PUBLIC :: sub_dnVec1_wADDTO_dnVec2, sub_dnVec1_PLUS_dnVec2_TO_dnVec3, dnVec2_wPLUS_dnVec3_TO_dnVec1
      PUBLIC :: sub_ZERO_TO_dnVec, test_ZERO_OF_dnVec, Vec_wADDTO_dnVec2_ider

      PUBLIC :: sub_dot_product_dnVec1_dnVec2_TO_dnS, Sub_crossproduct_dnVec1_dnVec2_TO_dnVec3
      PUBLIC :: dnVec2_O_dnVec1_TO_dnVec3
      PUBLIC :: sub_dnVec1_PROD_dnS2_TO_dnVec3

      ! with new type: dnS_t (QML)
      !PUBLIC :: sub_dnVec_TO_dnSt,sub_dnSt_TO_dnVec

      CONTAINS
!
!================================================================
!
!     allocation
!
!================================================================

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE alloc_dnVec(dnVec,nb_var_vec,nb_var_deriv,nderiv)
        TYPE (Type_dnVec) :: dnVec
        integer, optional :: nb_var_vec,nb_var_deriv,nderiv
        integer :: nd,nv,err_mem

        IF (present(nderiv)) dnVec%nderiv = nderiv
        IF (present(nb_var_deriv)) dnVec%nb_var_deriv = nb_var_deriv
        IF (present(nb_var_vec)) dnVec%nb_var_vec = nb_var_vec

        IF (dnVec%nb_var_deriv == 0) dnVec%nderiv = 0

        nd = dnVec%nb_var_deriv
        nv = dnVec%nb_var_vec

        IF (dnVec%alloc) RETURN
        dnVec%alloc = .TRUE.

!        write(out_unitp,*) 'BEGINNING alloc_dnVec'
!        write(out_unitp,*) 'nderiv',dnVec%nderiv
!        write(out_unitp,*) 'nb_var_deriv',dnVec%nb_var_deriv
!        write(out_unitp,*) 'nb_var_vec',dnVec%nb_var_vec


        IF (nv > 0) THEN
          CALL alloc_array(dnVec%d0,(/ nv /),'dnVec%d0','alloc_dnVec')
          dnVec%d0(:) = ZERO

          IF (dnVec%nderiv >= 1) THEN
            CALL alloc_array(dnVec%d1,(/ nv,nd /),'dnVec%d1','alloc_dnVec')
            dnVec%d1(:,:) = ZERO
          END IF

          IF (dnVec%nderiv >= 2) THEN
            CALL alloc_array(dnVec%d2,(/ nv,nd,nd /),'dnVec%d2','alloc_dnVec')
            dnVec%d2(:,:,:) = ZERO
          END IF

          IF (dnVec%nderiv >= 3) THEN
            CALL alloc_array(dnVec%d3,(/ nv,nd,nd,nd /),'dnVec%d3','alloc_dnVec')
            dnVec%d3(:,:,:,:) = ZERO
          END IF

          IF (dnVec%nderiv > 3) THEN
            write(out_unitp,*) ' ERROR in alloc_dnVec'
            write(out_unitp,*) ' nderiv MUST be < 4',dnVec%nderiv
            STOP
          END IF
        ELSE
          write(out_unitp,*) ' ERROR in alloc_dnVec'
          write(out_unitp,*) ' nb_var_vec MUST be > 0',nv
          STOP
        END IF

!       write(out_unitp,*) 'associated(d0)',associated(dnVec%d0)
!       write(out_unitp,*) 'associated(d1)',associated(dnVec%d1)
!       write(out_unitp,*) 'associated(d2)',associated(dnVec%d2)
!       write(out_unitp,*) 'associated(d3)',associated(dnVec%d3)
!        write(out_unitp,*) 'END alloc_dnVec'

      END SUBROUTINE alloc_dnVec

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE dealloc_dnVec(dnVec)
        TYPE (Type_dnVec) :: dnVec

        integer :: nd,nv,err_mem

        nd = dnVec%nb_var_deriv
        nv = dnVec%nb_var_vec

!        write(out_unitp,*) 'BEGINNING alloc_dnVec'
!        write(out_unitp,*) 'nderiv',dnVec%nderiv
!        write(out_unitp,*) 'nb_var_deriv',dnVec%nb_var_deriv
!        write(out_unitp,*) 'nb_var_vec',dnVec%nb_var_vec

        IF (associated(dnVec%d0)) THEN
          CALL dealloc_array(dnVec%d0,'dnVec%d0','dealloc_dnVec')
        END IF

        IF (associated(dnVec%d1)) THEN
          CALL dealloc_array(dnVec%d1,'dnVec%d1','dealloc_dnVec')
        END IF

        IF (associated(dnVec%d2)) THEN
          CALL dealloc_array(dnVec%d2,'dnVec%d2','dealloc_dnVec')
        END IF

        IF (associated(dnVec%d3)) THEN
          CALL dealloc_array(dnVec%d3,'dnVec%d3','dealloc_dnVec')
        END IF


        dnVec%alloc    = .FALSE.


        dnVec%nderiv       = 0
        dnVec%nb_var_deriv = 0
        dnVec%nb_var_vec   = 0
        !write(out_unitp,*) 'END dealloc_dnVec'

      END SUBROUTINE dealloc_dnVec

      SUBROUTINE alloc_array_OF_dnVecdim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      TYPE (Type_dnVec), pointer, intent(inout) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=1
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_dnVecdim1'
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
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Type_dnVec')

      END SUBROUTINE alloc_array_OF_dnVecdim1
      SUBROUTINE dealloc_array_OF_dnVecdim1(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (Type_dnVec), pointer, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_dnVecdim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN
       IF (.NOT. associated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Type_dnVec')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_dnVecdim1

      SUBROUTINE alloc_array_OF_dnVecdim2(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      TYPE (Type_dnVec), pointer, intent(inout) :: tab(:,:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub


      integer, parameter :: ndim=2
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_dnVecdim2'
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
         allocate(tab(tab_lb(1):tab_ub(1),                              &
                      tab_lb(2):tab_ub(2)),stat=err_mem)
       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1),tab_ub(2)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Type_dnVec')

      END SUBROUTINE alloc_array_OF_dnVecdim2
      SUBROUTINE dealloc_array_OF_dnVecdim2(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (Type_dnVec), pointer, intent(inout) :: tab(:,:)
      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_dnVecdim2'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN

       IF (.NOT. associated(tab))                                       &
                 CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Type_dnVec')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_dnVecdim2

!================================================================
!
!     check if alloc has been done
!
!================================================================
      !!@description: TODO
      !!@param: TODO
      SUBROUTINE check_alloc_dnVec(A,name_A,name_sub)
        TYPE (Type_dnVec), intent(in) :: A
        character (len=*), intent(in) :: name_A
        character (len=*), intent(in) :: name_sub

        IF ( .NOT. A%alloc) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) name_A,' has NOT been allocated with "alloc_dnVec"'
          write(out_unitp,*) ' CHECK the source!!!!!'
          STOP
        END IF
      END SUBROUTINE check_alloc_dnVec

  FUNCTION get_nderiv_FROM_dnVec(Vec) RESULT(nderiv)

    integer                          :: nderiv
    TYPE (Type_dnVec), intent(in)    :: Vec

    nderiv = Vec%nderiv

    IF (.NOT. associated(Vec%d0)) THEN
      nderiv = -1
    ELSE IF (.NOT. associated(Vec%d1)) THEN
      nderiv = 0
    ELSE IF (.NOT. associated(Vec%d2)) THEN
      nderiv = 1
    ELSE IF (.NOT. associated(Vec%d3)) THEN
      nderiv = 2
    ELSE
      nderiv = 3
    END IF

    IF (Vec%nderiv /= nderiv) THEN
      write(out_unitp,*) ' ERROR in get_nderiv_FROM_dnVec'
      write(out_unitp,*) '  Problem with nderiv in Vec'
      CALL Write_dnVec(Vec)
      STOP 'ERROR in get_nderiv_FROM_dnVec'
    END IF

    END FUNCTION get_nderiv_FROM_dnVec
  FUNCTION get_nb_var_deriv_FROM_dnVec(Vec) RESULT(nb_var_deriv)

    integer                          :: nb_var_deriv
    TYPE (Type_dnVec), intent(in)    :: Vec

    nb_var_deriv = Vec%nb_var_deriv

    IF (.NOT. associated(Vec%d1)) THEN
      nb_var_deriv = 0
    ELSE
      nb_var_deriv = size(Vec%d1,dim=2)
    END IF

    IF (Vec%nb_var_deriv /= nb_var_deriv .AND. Vec%nderiv > 0) THEN
      write(out_unitp,*) ' ERROR in get_nb_var_deriv_FROM_dnVec'
      write(out_unitp,*) '  Problem with nb_var_deriv in Vec'
      CALL Write_dnVec(Vec)
      STOP 'ERROR in get_nb_var_deriv_FROM_dnVec'
    END IF

    END FUNCTION get_nb_var_deriv_FROM_dnVec


!================================================================
!        write the derived type
!================================================================
      !!@description: write the derived type
      !!@param: TODO
      SUBROUTINE Write_dnVec(dnVec,nderiv)
        TYPE (Type_dnVec) :: dnVec
        integer, optional :: nderiv
        integer :: i,j,k
        integer :: nderiv_loc

        CALL check_alloc_dnVec(dnVec,'dnVec','Write_dnVec')

        nderiv_loc = dnVec%nderiv
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        write(out_unitp,*) 'BEGINNING Write dnVec'
        write(out_unitp,*) 'nderiv',dnVec%nderiv
        write(out_unitp,*) 'nb_var_vec,nb_var_deriv',                           &
                          dnVec%nb_var_vec,dnVec%nb_var_deriv
        IF (nderiv_loc >= 0 .AND. associated(dnVec%d0)) THEN
          write(out_unitp,*) 'd0'
          !write(out_unitp,*) dnVec%d0
          CALL Write_VecMat(dnVec%d0,out_unitp,10)
        END IF
        IF (nderiv_loc > 0 .AND. associated(dnVec%d1)) THEN
          DO i=1,dnVec%nb_var_deriv
            write(out_unitp,*) 'd1',i
            !write(out_unitp,*) dnVec%d1(:,i)
            CALL Write_VecMat(dnVec%d1(:,i),out_unitp,10)
          END DO
        END IF
        IF (nderiv_loc > 1 .AND. associated(dnVec%d2)) THEN
          DO i=1,dnVec%nb_var_deriv
          DO j=i,dnVec%nb_var_deriv
            write(out_unitp,*) 'd2',i,j
            !write(out_unitp,*) dnVec%d2(:,i,j)
            CALL Write_VecMat(dnVec%d2(:,i,j),out_unitp,10)
          END DO
          END DO
        END IF
        IF (nderiv_loc > 2 .AND. associated(dnVec%d3)) THEN
          DO i=1,dnVec%nb_var_deriv
          DO j=i,dnVec%nb_var_deriv
          DO k=j,dnVec%nb_var_deriv
            write(out_unitp,*) 'd3',i,j,k
            !write(out_unitp,*) dnVec%d3(:,i,j,k)
            CALL Write_VecMat(dnVec%d3(:,i,j,k),out_unitp,10)
          END DO
          END DO
          END DO
        END IF

        write(out_unitp,*) 'END Write dnVec'


      END SUBROUTINE Write_dnVec
      SUBROUTINE Vec_wADDTO_dnVec2_ider(Vec,w,dnVec2,ider,nderiv)
        real (kind=Rkind),  intent(in)            :: Vec(:)
        TYPE (Type_dnVec),  intent(inout)         :: dnVec2
        integer,            intent(in),  optional :: ider(:)
        integer,            intent(in),  optional :: nderiv
        real (kind=Rkind),  intent(in)            :: w

        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='Vec_wADDTO_dnVec2_ider'

        CALL check_alloc_dnVec(dnVec2,'dnVec2',name_sub)

        nderiv_loc = dnVec2%nderiv
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (size(Vec) /= dnVec2%nb_var_Vec) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' size(Vec) and nb_var_Vec must be equal',&
                                           size(Vec),dnVec2%nb_var_Vec
          write(out_unitp,*) ' CHECK the fortran source!!'
          STOP
        END IF
        IF (present(ider)) THEN
          IF (size(ider) > nderiv_loc) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' size(ider) cannot be > and nderiv_loc.'
            write(out_unitp,*) ' size(ider)',size(ider)
            write(out_unitp,*) ' dnVec2%nderiv',dnVec2%nderiv
            IF (present(nderiv)) write(out_unitp,*) ' nderiv',nderiv
            write(out_unitp,*) ' CHECK the fortran source!!'
            STOP
          END IF
          IF (any(ider < 1) .OR. any(ider > dnVec2%nb_var_deriv)) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' Some ider(:) values are out-of-range.'
            write(out_unitp,*) ' ider(:)',ider
            write(out_unitp,*) ' derivative range [1:',dnVec2%nderiv,']'
            write(out_unitp,*) ' CHECK the fortran source!!'
            STOP
          END IF
        END IF


        IF (present(ider)) THEN
          SELECT CASE (size(ider))
          CASE (3)
            dnVec2%d3(:,ider(1),ider(2),ider(3)) = w*Vec +              &
                                    dnVec2%d3(:,ider(1),ider(2),ider(3))
          CASE (2)
            dnVec2%d2(:,ider(1),ider(2)) = w*Vec + dnVec2%d2(:,ider(1),ider(2))
          CASE (1)
            dnVec2%d1(:,ider(1)) = w*Vec + dnVec2%d1(:,ider(1))
          CASE (0)
            dnVec2%d0(:) = w*Vec + dnVec2%d0
          CASE Default
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' nderiv_loc > 3 is NOT possible',nderiv_loc
            write(out_unitp,*) 'It should never append! Check the source'
            STOP
          END SELECT
        ELSE
            dnVec2%d0(:) = w*Vec + dnVec2%d0
        END IF

      END SUBROUTINE Vec_wADDTO_dnVec2_ider
!================================================================
!        dnS2 = dnS1 , dnVec2 = dnVec1 ...
!        transfer Vec(iVec) => R or R => Vec(iVec)
!================================================================
      !!@description: TODO
      !!@param: TODO
!      SUBROUTINE sub_dnVec_TO_dnSt(dnVec,dnS,iVec)
!      USE mod_QML_dnS
!
!        TYPE (Type_dnVec) :: dnVec
!        TYPE (dnS_t)      :: dnS
!        integer           :: iVec
!
!        character (len=*), parameter :: name_sub='sub_dnVec_TO_dnSt'
!
!        CALL check_alloc_dnVec(dnVec,'dnVec',name_sub)
!
!        IF (iVec < 1 .OR. iVec > dnVec%nb_var_Vec) THEN
!          write(out_unitp,*) ' ERROR in ',name_sub
!          write(out_unitp,*) ' iVec < 1 or iVec > dnVec%nb_var_Vec',            &
!                    dnVec%nb_var_vec,iVec
!          STOP
!        END IF
!
!        SELECT CASE (dnVec%nderiv)
!        CASE (0)
!          CALL QML_set_dnS(dnS,dnVec%d0(iVec))
!        CASE (1)
!          CALL QML_set_dnS(dnS,dnVec%d0(iVec),     &
!                               dnVec%d1(iVec,:))
!        CASE (2)
!          CALL QML_set_dnS(dnS,dnVec%d0(iVec),     &
!                               dnVec%d1(iVec,:),   &
!                               dnVec%d2(iVec,:,:))
!        CASE (3)
!          CALL QML_set_dnS(dnS,dnVec%d0(iVec),     &
!                               dnVec%d1(iVec,:),   &
!                               dnVec%d2(iVec,:,:), &
!                               dnVec%d3(iVec,:,:,:))
!        END SELECT
!
!      END SUBROUTINE sub_dnVec_TO_dnSt
!      SUBROUTINE sub_dnSt_TO_dnVec(dnS,dnVec,iVec)
!      USE mod_QML_dnS
!
!        TYPE (Type_dnVec) :: dnVec
!        TYPE (dnS_t)      :: dnS
!        integer           :: iVec
!
!        integer           :: nderiv,nb_var_deriv
!        character (len=*), parameter :: name_sub='sub_dnSt_TO_dnVec'
!
!        nb_var_deriv = QML_get_ndim_FROM_dnS(dnS)
!        nderiv       = QML_get_nderiv_FROM_dnS(dnS)
!
!        CALL check_alloc_dnVec(dnVec,'dnVec',name_sub)
!
!        nderiv = min(dnVec%nderiv,nderiv)
!
!        IF (nderiv > 0) THEN ! it has to be test because ndim of dnS is zero when nderiv=0
!        IF (dnVec%nb_var_deriv /= nb_var_deriv) THEN
!          write(out_unitp,*) ' ERROR in ',name_sub
!          write(out_unitp,*) ' nb_var_deriv in dnVec and dnS are different!',   &
!                    dnVec%nb_var_deriv,nb_var_deriv
!          STOP
!        END IF
!        END IF
!
!        IF (iVec < 1 .OR. iVec > dnVec%nb_var_Vec) THEN
!          write(out_unitp,*) ' ERROR in ',name_sub
!          write(out_unitp,*) ' iVec < 1 or iVec > dnVec%nb_var_Vec',            &
!                    dnVec%nb_var_vec,iVec
!          STOP
!        END IF
!
!        SELECT CASE (nderiv)
!        CASE (0)
!          CALL QML_sub_get_dn_FROM_dnS(dnS,dnVec%d0(iVec))
!        CASE (1)
!          CALL QML_sub_get_dn_FROM_dnS(dnS,dnVec%d0(iVec),     &
!                                           dnVec%d1(iVec,:))
!        CASE (2)
!          CALL QML_sub_get_dn_FROM_dnS(dnS,dnVec%d0(iVec),     &
!                                           dnVec%d1(iVec,:),   &
!                                           dnVec%d2(iVec,:,:))
!        CASE (3)
!          CALL QML_sub_get_dn_FROM_dnS(dnS,dnVec%d0(iVec),     &
!                                           dnVec%d1(iVec,:),   &
!                                           dnVec%d2(iVec,:,:), &
!                                           dnVec%d3(iVec,:,:,:))
!        END SELECT
!
!      END SUBROUTINE sub_dnSt_TO_dnVec
      SUBROUTINE sub_dnVec_TO_dnS(dnVec,dnS,iVec,nderiv)
      use mod_dnS, only: type_dns, alloc_dns, check_alloc_dns, write_dns

        TYPE (Type_dnVec) :: dnVec
        TYPE (Type_dnS)   :: dnS
        integer :: iVec
        integer, optional :: nderiv

        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='sub_dnVec_TO_dnS'

        CALL check_alloc_dnVec(dnVec,'dnVec',name_sub)
        !CALL check_alloc_dnS(dnS,'dnS',name_sub)
        IF (.NOT. dnS%alloc) THEN
          CALL alloc_dnS(dnS,dnVec%nb_var_deriv,dnVec%nderiv)
        END IF

        nderiv_loc = min(dnVec%nderiv,dnS%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (dnVec%nb_var_deriv /= dnS%nb_var_deriv) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nb_var_deriv in dnVec and dnS are different!',   &
                    dnVec%nb_var_deriv,dnS%nb_var_deriv
          STOP
        END IF
        IF (iVec < 1 .OR. iVec > dnVec%nb_var_Vec) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' iVec < 1 or iVec > dnVec%nb_var_Vec',            &
                    dnVec%nb_var_vec,iVec
          STOP
        END IF
        IF (nderiv_loc == 0) THEN
           dnS%d0 = dnVec%d0(iVec)
        ELSE IF (nderiv_loc == 1) THEN
           dnS%d0 = dnVec%d0(iVec)
           dnS%d1(:) = dnVec%d1(iVec,:)
        ELSE IF (nderiv_loc == 2) THEN
           dnS%d0 = dnVec%d0(iVec)
           dnS%d1(:) = dnVec%d1(iVec,:)
           dnS%d2(:,:) = dnVec%d2(iVec,:,:)
        ELSE IF (nderiv_loc == 3) THEN
           dnS%d0 = dnVec%d0(iVec)
           dnS%d1(:) = dnVec%d1(iVec,:)
           dnS%d2(:,:) = dnVec%d2(iVec,:,:)
           dnS%d3(:,:,:) = dnVec%d3(iVec,:,:,:)
        ELSE
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nderiv_loc > 4 is NOT possible',nderiv_loc
          write(out_unitp,*) 'It souhld never append! Check the source'
          STOP
        END IF
      END SUBROUTINE sub_dnVec_TO_dnS

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE sub_dnS_TO_dnVec(dnS,dnVec,iVec,nderiv)
      use mod_dnS, only: type_dns, check_alloc_dns, write_dns

        TYPE (Type_dnVec) :: dnVec
        TYPE (Type_dnS)   :: dnS
        integer :: iVec
        integer, optional :: nderiv

        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='sub_dnS_TO_dnVec'

        CALL check_alloc_dnVec(dnVec,'dnVec',name_sub)
        CALL check_alloc_dnS(dnS,'dnS',name_sub)

        nderiv_loc = min(dnVec%nderiv,dnS%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (dnVec%nb_var_deriv /= dnS%nb_var_deriv) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nb_var_deriv in dnVec and dnS are different!',   &
                    dnVec%nb_var_deriv,dnS%nb_var_deriv
          STOP
        END IF
        IF (iVec < 1 .OR. iVec > dnVec%nb_var_Vec) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' iVec < 1 or iVec > dnVec%nb_var_Vec',            &
                    dnVec%nb_var_vec,iVec
          STOP
        END IF
        IF (nderiv_loc == 0) THEN
           dnVec%d0(iVec) = dnS%d0
        ELSE IF (nderiv_loc == 1) THEN
           dnVec%d0(iVec) = dnS%d0
           dnVec%d1(iVec,:) = dnS%d1(:)
        ELSE IF (nderiv_loc == 2) THEN
           dnVec%d0(iVec) = dnS%d0
           dnVec%d1(iVec,:) = dnS%d1(:)
           dnVec%d2(iVec,:,:) = dnS%d2(:,:)
        ELSE IF (nderiv_loc == 3) THEN
           dnVec%d0(iVec) = dnS%d0
           dnVec%d1(iVec,:) = dnS%d1(:)
           dnVec%d2(iVec,:,:) = dnS%d2(:,:)
           dnVec%d3(iVec,:,:,:) = dnS%d3(:,:,:)
        ELSE
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nderiv_loc > 4 is NOT possible',nderiv_loc
          write(out_unitp,*) 'It souhld never append! Check the source'
          STOP
        END IF
      END SUBROUTINE sub_dnS_TO_dnVec

      SUBROUTINE sub_dot_product_dnVec1_dnVec2_TO_dnS(dnVec1,dnVec2,dnS,nderiv)
      use mod_dnS, only: type_dns, alloc_dns, check_alloc_dns, write_dns

        TYPE (Type_dnVec) :: dnVec1,dnVec2
        TYPE (Type_dnS)   :: dnS
        integer, optional :: nderiv

        integer :: id,jd,kd,nderiv_loc
        character (len=*), parameter :: name_sub='sub_dot_product_dnVec1_dnVec2_TO_dnS'

        CALL check_alloc_dnVec(dnVec1,'dnVec',name_sub)
        CALL check_alloc_dnVec(dnVec2,'dnVec',name_sub)
        IF (.NOT. dnS%alloc) THEN
          CALL alloc_dnS(dnS,dnVec1%nb_var_deriv,dnVec1%nderiv)
        END IF

        nderiv_loc = min(dnVec1%nderiv,dnVec2%nderiv,dnS%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (dnVec1%nb_var_deriv /= dnS%nb_var_deriv .OR.                &
            dnVec2%nb_var_deriv /= dnS%nb_var_deriv .OR.                &
            dnVec1%nb_var_deriv /= dnVec2%nb_var_deriv) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nb_var_deriv in dnVec1 and dnVec2 and dnS are different!',   &
                    dnVec1%nb_var_deriv,dnVec2%nb_var_deriv,dnS%nb_var_deriv
          STOP
        END IF

        IF (nderiv_loc == 0) THEN
           dnS%d0 = dot_product(dnVec1%d0(:),dnVec2%d0(:))
        ELSE IF (nderiv_loc == 1) THEN
           dnS%d0 = dot_product(dnVec1%d0(:),dnVec2%d0(:))
           DO id=1,dnVec1%nb_var_deriv
             dnS%d1(id) = dot_product(dnVec1%d1(:,id),dnVec2%d0(:)) +   &
                          dot_product(dnVec1%d0(:),dnVec2%d1(:,id))
           END DO
        ELSE IF (nderiv_loc == 2) THEN
           dnS%d0 = dot_product(dnVec1%d0(:),dnVec2%d0(:))
           DO id=1,dnVec1%nb_var_deriv
             dnS%d1(id) = dot_product(dnVec1%d1(:,id),dnVec2%d0(:)) +   &
                          dot_product(dnVec1%d0(:),dnVec2%d1(:,id))
           END DO
           DO id=1,dnVec1%nb_var_deriv
           DO jd=1,dnVec1%nb_var_deriv
             dnS%d2(id,jd) = dot_product(dnVec1%d2(:,id,jd),dnVec2%d0(:)) +&
                             dot_product(dnVec1%d1(:,id),dnVec2%d1(:,jd)) +&
                             dot_product(dnVec1%d1(:,jd),dnVec2%d1(:,id)) +&
                             dot_product(dnVec1%d0(:),dnVec2%d2(:,id,jd))
           END DO
           END DO
        ELSE IF (nderiv_loc == 3) THEN
           dnS%d0 = dot_product(dnVec1%d0(:),dnVec2%d0(:))
           DO id=1,dnVec1%nb_var_deriv
             dnS%d1(id) = dot_product(dnVec1%d1(:,id),dnVec2%d0(:)) +   &
                          dot_product(dnVec1%d0(:),dnVec2%d1(:,id))
           END DO
           DO id=1,dnVec1%nb_var_deriv
           DO jd=1,dnVec1%nb_var_deriv
             dnS%d2(id,jd) = dot_product(dnVec1%d2(:,id,jd),dnVec2%d0(:)) +&
                             dot_product(dnVec1%d1(:,id),dnVec2%d1(:,jd)) +&
                             dot_product(dnVec1%d1(:,jd),dnVec2%d1(:,id)) +&
                             dot_product(dnVec1%d0(:),dnVec2%d2(:,id,jd))
           END DO
           END DO
           DO id=1,dnVec1%nb_var_deriv
           DO jd=1,dnVec1%nb_var_deriv
           DO kd=1,dnVec1%nb_var_deriv
             dnS%d3(id,jd,kd) = dot_product(dnVec1%d3(:,id,jd,kd),dnVec2%d0(:)) +&
                                dot_product(dnVec1%d2(:,id,jd),dnVec2%d1(:,kd)) +&
                                dot_product(dnVec1%d2(:,id,kd),dnVec2%d1(:,jd)) +&
                                dot_product(dnVec1%d1(:,id),dnVec2%d2(:,jd,kd)) +&
                                dot_product(dnVec1%d2(:,jd,kd),dnVec2%d1(:,id)) +&
                                dot_product(dnVec1%d1(:,jd),dnVec2%d2(:,id,kd)) +&
                                dot_product(dnVec1%d1(:,kd),dnVec2%d2(:,id,jd)) +&
                                dot_product(dnVec1%d0(:),dnVec2%d3(:,id,jd,kd))
           END DO
           END DO
           END DO
        ELSE
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nderiv_loc > 4 is NOT possible',nderiv_loc
          write(out_unitp,*) 'It souhld never append! Check the source'
          STOP
        END IF
      END SUBROUTINE sub_dot_product_dnVec1_dnVec2_TO_dnS


!================================================================
!       cross product  v3 = v1 x v2
!================================================================
      SUBROUTINE Sub_crossproduct_dnVec1_dnVec2_TO_dnVec3(dnVec1,dnVec2,dnVec3,nderiv)
      USE mod_system
      IMPLICIT NONE

      TYPE (Type_dnVec) :: dnVec1,dnVec2,dnVec3
      integer           :: nderiv

      integer i,j,k


!     -----------------------------------------------------------------
!     logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      character (len=*), parameter :: name_sub= 'Sub_crossproduct_dnVec1_dnVec2_TO_dnVec3'
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Vec1'
        CALL Write_dnVec(dnVec1)
        write(out_unitp,*) 'Vec2'
        CALL Write_dnVec(dnVec2)
        write(out_unitp,*) 'Vec3'
        CALL Write_dnVec(dnVec3)
      END IF
!     -----------------------------------------------------------------

      IF (dnVec1%nb_var_vec /= 3 .OR. dnVec2%nb_var_vec /= 3 .OR. dnVec3%nb_var_vec /= 3) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' The size of the vectors must be 3',        &
                  dnVec1%nb_var_vec,dnVec2%nb_var_vec,dnVec3%nb_var_vec
        write(out_unitp,*) 'It should never append! Check the source'
        STOP
      END IF


!     -----------------------------------------------------------------
      IF (nderiv == 0) THEN
        dnVec3%d0(1) =  dnVec1%d0(2)*dnVec2%d0(3) -                     &
                        dnVec1%d0(3)*dnVec2%d0(2)
        dnVec3%d0(2) = -dnVec1%d0(1)*dnVec2%d0(3) +                     &
                        dnVec1%d0(3)*dnVec2%d0(1)
        dnVec3%d0(3) =  dnVec1%d0(1)*dnVec2%d0(2) -                     &
                        dnVec1%d0(2)*dnVec2%d0(1)
!     -----------------------------------------------------------------
      ELSE IF (nderiv == 1) THEN
        dnVec3%d0(1) =  dnVec1%d0(2)*dnVec2%d0(3) -                     &
                        dnVec1%d0(3)*dnVec2%d0(2)
        dnVec3%d0(2) = -dnVec1%d0(1)*dnVec2%d0(3) +                     &
                        dnVec1%d0(3)*dnVec2%d0(1)
        dnVec3%d0(3) =  dnVec1%d0(1)*dnVec2%d0(2) -                     &
                        dnVec1%d0(2)*dnVec2%d0(1)

        DO i=1,dnVec1%nb_var_deriv
          dnVec3%d1(1,i) =  dnVec1%d1(2,i)*dnVec2%d0(3) +               &
                            dnVec1%d0(2)*dnVec2%d1(3,i) -               &
                            dnVec1%d1(3,i)*dnVec2%d0(2) -               &
                            dnVec1%d0(3)*dnVec2%d1(2,i)
          dnVec3%d1(2,i) = -dnVec1%d1(1,i)*dnVec2%d0(3) -               &
                            dnVec1%d0(1)*dnVec2%d1(3,i) +               &
                            dnVec1%d1(3,i)*dnVec2%d0(1) +               &
                            dnVec1%d0(3)*dnVec2%d1(1,i)
          dnVec3%d1(3,i) =  dnVec1%d1(1,i)*dnVec2%d0(2) +               &
                            dnVec1%d0(1)*dnVec2%d1(2,i) -               &
                            dnVec1%d1(2,i)*dnVec2%d0(1) -               &
                            dnVec1%d0(2)*dnVec2%d1(1,i)
        END DO
!      -----------------------------------------------------------------
      ELSE IF (nderiv == 2) THEN
        dnVec3%d0(1) =  dnVec1%d0(2)*dnVec2%d0(3) -                     &
                        dnVec1%d0(3)*dnVec2%d0(2)
        dnVec3%d0(2) = -dnVec1%d0(1)*dnVec2%d0(3) +                     &
                        dnVec1%d0(3)*dnVec2%d0(1)
        dnVec3%d0(3) =  dnVec1%d0(1)*dnVec2%d0(2) -                     &
                        dnVec1%d0(2)*dnVec2%d0(1)

        DO i=1,dnVec1%nb_var_deriv
          dnVec3%d1(1,i) =  dnVec1%d1(2,i)*dnVec2%d0(3) +               &
                            dnVec1%d0(2)*dnVec2%d1(3,i) -               &
                            dnVec1%d1(3,i)*dnVec2%d0(2) -               &
                            dnVec1%d0(3)*dnVec2%d1(2,i)
          dnVec3%d1(2,i) = -dnVec1%d1(1,i)*dnVec2%d0(3) -               &
                            dnVec1%d0(1)*dnVec2%d1(3,i) +               &
                            dnVec1%d1(3,i)*dnVec2%d0(1) +               &
                            dnVec1%d0(3)*dnVec2%d1(1,i)
          dnVec3%d1(3,i) =  dnVec1%d1(1,i)*dnVec2%d0(2) +               &
                            dnVec1%d0(1)*dnVec2%d1(2,i) -               &
                            dnVec1%d1(2,i)*dnVec2%d0(1) -               &
                            dnVec1%d0(2)*dnVec2%d1(1,i)
        END DO


        DO i=1,dnVec1%nb_var_deriv
        DO j=1,dnVec1%nb_var_deriv

         dnVec3%d2(1,i,j) =  dnVec1%d2(2,i,j) * dnVec2%d0(3)     +      &
                             dnVec1%d1(2,i)   * dnVec2%d1(3,j)   +      &
                             dnVec1%d1(2,j)   * dnVec2%d1(3,i)   +      &
                             dnVec1%d0(2)     * dnVec2%d2(3,i,j) -      &
                             dnVec1%d2(3,i,j) * dnVec2%d0(2)     -      &
                             dnVec1%d1(3,i)   * dnVec2%d1(2,j)   -      &
                             dnVec1%d1(3,j)   * dnVec2%d1(2,i)   -      &
                             dnVec1%d0(3)     * dnVec2%d2(2,i,j)

         dnVec3%d2(2,i,j) = -dnVec1%d2(1,i,j) * dnVec2%d0(3)     -      &
                             dnVec1%d1(1,i)   * dnVec2%d1(3,j)   -      &
                             dnVec1%d1(1,j)   * dnVec2%d1(3,i)   -      &
                             dnVec1%d0(1)     * dnVec2%d2(3,i,j) +      &
                             dnVec1%d2(3,i,j) * dnVec2%d0(1)     +      &
                             dnVec1%d1(3,i)   * dnVec2%d1(1,j)   +      &
                             dnVec1%d1(3,j)   * dnVec2%d1(1,i)   +      &
                             dnVec1%d0(3)     * dnVec2%d2(1,i,j)

         dnVec3%d2(3,i,j) =  dnVec1%d2(1,i,j) * dnVec2%d0(2)     +      &
                             dnVec1%d1(1,i)   * dnVec2%d1(2,j)   +      &
                             dnVec1%d1(1,j)   * dnVec2%d1(2,i)   +      &
                             dnVec1%d0(1)     * dnVec2%d2(2,i,j) -      &
                             dnVec1%d2(2,i,j) * dnVec2%d0(1)     -      &
                             dnVec1%d1(2,i)   * dnVec2%d1(1,j)   -      &
                             dnVec1%d1(2,j)   * dnVec2%d1(1,i)   -      &
                             dnVec1%d0(2)     * dnVec2%d2(1,i,j)

        END DO
        END DO
!      -----------------------------------------------------------------
      ELSE IF (nderiv == 3) THEN
        dnVec3%d0(1) =  dnVec1%d0(2)*dnVec2%d0(3) -                     &
                        dnVec1%d0(3)*dnVec2%d0(2)
        dnVec3%d0(2) = -dnVec1%d0(1)*dnVec2%d0(3) +                     &
                        dnVec1%d0(3)*dnVec2%d0(1)
        dnVec3%d0(3) =  dnVec1%d0(1)*dnVec2%d0(2) -                     &
                        dnVec1%d0(2)*dnVec2%d0(1)

        DO i=1,dnVec1%nb_var_deriv
          dnVec3%d1(1,i) =  dnVec1%d1(2,i)*dnVec2%d0(3) +               &
                            dnVec1%d0(2)*dnVec2%d1(3,i) -               &
                            dnVec1%d1(3,i)*dnVec2%d0(2) -               &
                            dnVec1%d0(3)*dnVec2%d1(2,i)
          dnVec3%d1(2,i) = -dnVec1%d1(1,i)*dnVec2%d0(3) -               &
                            dnVec1%d0(1)*dnVec2%d1(3,i) +               &
                            dnVec1%d1(3,i)*dnVec2%d0(1) +               &
                            dnVec1%d0(3)*dnVec2%d1(1,i)
          dnVec3%d1(3,i) =  dnVec1%d1(1,i)*dnVec2%d0(2) +               &
                            dnVec1%d0(1)*dnVec2%d1(2,i) -               &
                            dnVec1%d1(2,i)*dnVec2%d0(1) -               &
                            dnVec1%d0(2)*dnVec2%d1(1,i)
        END DO

        DO i=1,dnVec1%nb_var_deriv
        DO j=1,dnVec1%nb_var_deriv

         dnVec3%d2(1,i,j) = dnVec1%d2(2,i,j) * dnVec2%d0(3)     +       &
                            dnVec1%d1(2,i)   * dnVec2%d1(3,j)   +       &
                            dnVec1%d1(2,j)   * dnVec2%d1(3,i)   +       &
                            dnVec1%d0(2)     * dnVec2%d2(3,i,j) -       &
                            dnVec1%d2(3,i,j) * dnVec2%d0(2)     -       &
                            dnVec1%d1(3,i)   * dnVec2%d1(2,j)   -       &
                            dnVec1%d1(3,j)   * dnVec2%d1(2,i)   -       &
                            dnVec1%d0(3)     * dnVec2%d2(2,i,j)

         dnVec3%d2(2,i,j) =-dnVec1%d2(1,i,j) * dnVec2%d0(3)     -       &
                            dnVec1%d1(1,i)   * dnVec2%d1(3,j)   -       &
                            dnVec1%d1(1,j)   * dnVec2%d1(3,i)   -       &
                            dnVec1%d0(1)     * dnVec2%d2(3,i,j) +       &
                            dnVec1%d2(3,i,j) * dnVec2%d0(1)     +       &
                            dnVec1%d1(3,i)   * dnVec2%d1(1,j)   +       &
                            dnVec1%d1(3,j)   * dnVec2%d1(1,i)   +       &
                            dnVec1%d0(3)     * dnVec2%d2(1,i,j)

         dnVec3%d2(3,i,j) = dnVec1%d2(1,i,j) * dnVec2%d0(2)     +       &
                            dnVec1%d1(1,i)   * dnVec2%d1(2,j)   +       &
                            dnVec1%d1(1,j)   * dnVec2%d1(2,i)   +       &
                            dnVec1%d0(1)     * dnVec2%d2(2,i,j) -       &
                            dnVec1%d2(2,i,j) * dnVec2%d0(1)     -       &
                            dnVec1%d1(2,i)   * dnVec2%d1(1,j)   -       &
                            dnVec1%d1(2,j)   * dnVec2%d1(1,i)   -       &
                            dnVec1%d0(2)     * dnVec2%d2(1,i,j)

        END DO
        END DO

        DO i=1,dnVec1%nb_var_deriv
        DO j=1,dnVec1%nb_var_deriv
        DO k=1,dnVec1%nb_var_deriv
         dnVec3%d3(1,i,j,k) = dnVec1%d3(2,i,j,k) * dnVec2%d0(3)       + &
                              dnVec1%d2(2,j,k)   * dnVec2%d1(3,i)     + &
                              dnVec1%d2(2,i,k)   * dnVec2%d1(3,j)     + &
                              dnVec1%d2(2,i,j)   * dnVec2%d1(3,k)     + &
                              dnVec1%d1(2,i)     * dnVec2%d2(3,j,k)   + &
                              dnVec1%d1(2,j)     * dnVec2%d2(3,i,k)   + &
                              dnVec1%d1(2,k)     * dnVec2%d2(3,i,j)   + &
                              dnVec1%d0(2)       * dnVec2%d3(3,i,j,k) - &
                              dnVec1%d3(3,i,j,k) * dnVec2%d0(2)       - &
                              dnVec1%d2(3,j,k)   * dnVec2%d1(2,i)     - &
                              dnVec1%d2(3,i,k)   * dnVec2%d1(2,j)     - &
                              dnVec1%d2(3,i,j)   * dnVec2%d1(2,k)     - &
                              dnVec1%d1(3,i)     * dnVec2%d2(2,j,k)   - &
                              dnVec1%d1(3,j)     * dnVec2%d2(2,i,k)   - &
                              dnVec1%d1(3,k)     * dnVec2%d2(2,i,j)   - &
                              dnVec1%d0(3)       * dnVec2%d3(2,i,j,k)

         dnVec3%d3(2,i,j,k) =-dnVec1%d3(1,i,j,k) * dnVec2%d0(3)       - &
                              dnVec1%d2(1,j,k)   * dnVec2%d1(3,i)     - &
                              dnVec1%d2(1,i,k)   * dnVec2%d1(3,j)     - &
                              dnVec1%d2(1,i,j)   * dnVec2%d1(3,k)     - &
                              dnVec1%d1(1,i)     * dnVec2%d2(3,j,k)   - &
                              dnVec1%d1(1,j)     * dnVec2%d2(3,i,k)   - &
                              dnVec1%d1(1,k)     * dnVec2%d2(3,i,j)   - &
                              dnVec1%d0(1)       * dnVec2%d3(3,i,j,k) + &
                              dnVec1%d3(3,i,j,k) * dnVec2%d0(1)       + &
                              dnVec1%d2(3,j,k)   * dnVec2%d1(1,i)     + &
                              dnVec1%d2(3,i,k)   * dnVec2%d1(1,j)     + &
                              dnVec1%d2(3,i,j)   * dnVec2%d1(1,k)     + &
                              dnVec1%d1(3,i)     * dnVec2%d2(1,j,k)   + &
                              dnVec1%d1(3,j)     * dnVec2%d2(1,i,k)   + &
                              dnVec1%d1(3,k)     * dnVec2%d2(1,i,j)   + &
                              dnVec1%d0(3)       * dnVec2%d3(1,i,j,k)


         dnVec3%d3(3,i,j,k) = dnVec1%d3(1,i,j,k) * dnVec2%d0(2)       + &
                              dnVec1%d2(1,j,k)   * dnVec2%d1(2,i)     + &
                              dnVec1%d2(1,i,k)   * dnVec2%d1(2,j)     + &
                              dnVec1%d2(1,i,j)   * dnVec2%d1(2,k)     + &
                              dnVec1%d1(1,i)     * dnVec2%d2(2,j,k)   + &
                              dnVec1%d1(1,j)     * dnVec2%d2(2,i,k)   + &
                              dnVec1%d1(1,k)     * dnVec2%d2(2,i,j)   + &
                              dnVec1%d0(1)       * dnVec2%d3(2,i,j,k) - &
                              dnVec1%d3(2,i,j,k) * dnVec2%d0(1)       - &
                              dnVec1%d2(2,j,k)   * dnVec2%d1(1,i)     - &
                              dnVec1%d2(2,i,k)   * dnVec2%d1(1,j)     - &
                              dnVec1%d2(2,i,j)   * dnVec2%d1(1,k)     - &
                              dnVec1%d1(2,i)     * dnVec2%d2(1,j,k)   - &
                              dnVec1%d1(2,j)     * dnVec2%d2(1,i,k)   - &
                              dnVec1%d1(2,k)     * dnVec2%d2(1,i,j)   - &
                              dnVec1%d0(2)       * dnVec2%d3(1,i,j,k)

        END DO
        END DO
        END DO
      END IF
!      -----------------------------------------------------------------

!      -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Vec3'
        CALL Write_dnVec(dnVec3)
        write(out_unitp,*) 'END ',name_sub
      END IF
!      -----------------------------------------------------------------

        END SUBROUTINE Sub_crossproduct_dnVec1_dnVec2_TO_dnVec3

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE sub_dnVec2_TO_dnVec1(dnVec1,dnVec2)
        CLASS (Type_dnVec), intent(inout) :: dnVec1
        TYPE (Type_dnVec), intent(in)     :: dnVec2

        character (len=*), parameter :: name_sub='sub_dnVec2_TO_dnVec1'

        CALL check_alloc_dnVec(dnVec2,'dnVec2',name_sub)

        CALL dealloc_dnVec(dnVec1)
        CALL alloc_dnVec(dnVec1,dnVec2%nb_var_vec,dnVec2%nb_var_deriv,dnVec2%nderiv)

        IF (dnVec2%nderiv == 0) THEN
           dnVec1%d0 = dnVec2%d0
        ELSE IF (dnVec2%nderiv == 1) THEN
           dnVec1%d0 = dnVec2%d0
           dnVec1%d1 = dnVec2%d1
        ELSE IF (dnVec2%nderiv == 2) THEN
           dnVec1%d0 = dnVec2%d0
           dnVec1%d1 = dnVec2%d1
           dnVec1%d2 = dnVec2%d2
        ELSE IF (dnVec2%nderiv == 3) THEN
           dnVec1%d0 = dnVec2%d0
           dnVec1%d1 = dnVec2%d1
           dnVec1%d2 = dnVec2%d2
           dnVec1%d3 = dnVec2%d3
        ELSE
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nderiv > 3 is NOT possible',dnVec2%nderiv
          write(out_unitp,*) 'It should never append! Check the source'
          STOP
        END IF

      END SUBROUTINE sub_dnVec2_TO_dnVec1

      SUBROUTINE sub_dnVec1_TO_dnVec2(dnVec1,dnVec2,nderiv)
        TYPE (Type_dnVec) :: dnVec1,dnVec2
        integer, optional :: nderiv

        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='sub_dnVec1_TO_dnVec2'

        CALL check_alloc_dnVec(dnVec1,'dnVec1',name_sub)

        IF (.NOT. dnVec2%alloc) THEN
          CALL alloc_dnVec(dnVec2,dnVec1%nb_var_vec,dnVec1%nb_var_deriv,&
                                                          dnVec1%nderiv)
        END IF

        nderiv_loc = min(dnVec1%nderiv,dnVec2%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (dnVec1%nb_var_deriv /= dnVec2%nb_var_deriv) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nb_var_deriv in dnVec1 and dnVec2 are different!',&
                    dnVec1%nb_var_deriv,dnVec2%nb_var_deriv
          STOP
        END IF
        IF (dnVec1%nb_var_Vec /= dnVec2%nb_var_Vec) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nb_var_Vec in dnVec1 and dnVec2 are different!', &
                    dnVec1%nb_var_Vec,dnVec2%nb_var_Vec
          STOP
        END IF

        IF (nderiv_loc == 0) THEN
           dnVec2%d0 = dnVec1%d0
        ELSE IF (nderiv_loc == 1) THEN
           dnVec2%d0 = dnVec1%d0
           dnVec2%d1 = dnVec1%d1
        ELSE IF (nderiv_loc == 2) THEN
           dnVec2%d0 = dnVec1%d0
           dnVec2%d1 = dnVec1%d1
           dnVec2%d2 = dnVec1%d2
        ELSE IF (nderiv_loc == 3) THEN
           dnVec2%d0 = dnVec1%d0
           dnVec2%d1 = dnVec1%d1
           dnVec2%d2 = dnVec1%d2
           dnVec2%d3 = dnVec1%d3
        ELSE
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nderiv > 4 is NOT possible',nderiv
          write(out_unitp,*) 'It should never append! Check the source'
          STOP
        END IF

      END SUBROUTINE sub_dnVec1_TO_dnVec2

      SUBROUTINE sub_PartdnVec1_TO_PartdnVec2(dnVec1,i1,dnVec2,i2,ndim,nderiv)
        TYPE (Type_dnVec) :: dnVec1,dnVec2
        integer           :: i1,i2,ndim
        integer, optional :: nderiv

        integer           :: f1,f2
        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='sub_PartdnVec1_TO_PartdnVec2'

        CALL check_alloc_dnVec(dnVec1,'dnVec1',name_sub)

        IF (.NOT. dnVec2%alloc) THEN
          CALL alloc_dnVec(dnVec2,dnVec1%nb_var_vec,dnVec1%nb_var_deriv,&
                                                          dnVec1%nderiv)
        END IF

        nderiv_loc = min(dnVec1%nderiv,dnVec2%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (dnVec1%nb_var_deriv /= dnVec2%nb_var_deriv) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nb_var_deriv in dnVec1 and dnVec2 are different!',&
                    dnVec1%nb_var_deriv,dnVec2%nb_var_deriv
          STOP
        END IF
        !IF (dnVec1%nb_var_Vec /= dnVec2%nb_var_Vec) THEN
        !  write(out_unitp,*) ' ERROR in ',name_sub
        !  write(out_unitp,*) ' nb_var_Vec in dnVec1 and dnVec2 are different!', &
        !            dnVec1%nb_var_Vec,dnVec2%nb_var_Vec
        !  STOP
        !END IF

        f1 = i1-1+ndim
        IF (i1 < 1 .OR. f1 > dnVec1%nb_var_vec) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' WRONG values of i1,f1',i1,f1
          write(out_unitp,*) ' ndim,dnVec1%nb_var_vec',ndim,dnVec1%nb_var_vec
          STOP
        END IF

        f2 = i2-1+ndim
        IF (i2 < 1 .OR. f2 > dnVec2%nb_var_vec) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' WRONG values of i2,f2',i2,f2
          write(out_unitp,*) ' ndim,dnVec2%nb_var_vec',ndim,dnVec2%nb_var_vec
          STOP
        END IF



        IF (nderiv_loc == 0) THEN
           dnVec2%d0(i2:f2) = dnVec1%d0(i1:f1)
        ELSE IF (nderiv_loc == 1) THEN
           dnVec2%d0(i2:f2) = dnVec1%d0(i1:f1)
           dnVec2%d1(i2:f2,:) = dnVec1%d1(i1:f1,:)
        ELSE IF (nderiv_loc == 2) THEN
           dnVec2%d0(i2:f2) = dnVec1%d0(i1:f1)
           dnVec2%d1(i2:f2,:) = dnVec1%d1(i1:f1,:)
           dnVec2%d2(i2:f2,:,:) = dnVec1%d2(i1:f1,:,:)
        ELSE IF (nderiv_loc == 3) THEN
           dnVec2%d0(i2:f2) = dnVec1%d0(i1:f1)
           dnVec2%d1(i2:f2,:) = dnVec1%d1(i1:f1,:)
           dnVec2%d2(i2:f2,:,:) = dnVec1%d2(i1:f1,:,:)
           dnVec2%d3(i2:f2,:,:,:) = dnVec1%d3(i1:f1,:,:,:)
        ELSE
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nderiv > 4 is NOT possible',nderiv
          write(out_unitp,*) 'It should never append! Check the source'
          STOP
        END IF

      END SUBROUTINE sub_PartdnVec1_TO_PartdnVec2

      SUBROUTINE sub_dnVec1_TO_dnVec2_WithIvec(dnVec1,dnVec2,iVec,nderiv)
        TYPE (Type_dnVec) :: dnVec1,dnVec2
        integer :: iVec

        integer, optional :: nderiv

        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='sub_dnVec1_TO_dnVec2_WithIvec'

        CALL check_alloc_dnVec(dnVec1,'dnVec1',name_sub)

        IF (.NOT. dnVec2%alloc) THEN
          CALL alloc_dnVec(dnVec2,dnVec1%nb_var_vec,dnVec1%nb_var_deriv,&
                                                          dnVec1%nderiv)
        END IF

        nderiv_loc = min(dnVec1%nderiv,dnVec2%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (dnVec1%nb_var_deriv /= dnVec2%nb_var_deriv) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nb_var_deriv in dnVec1 and dnVec2 are different!',&
                    dnVec1%nb_var_deriv,dnVec2%nb_var_deriv
          STOP
        END IF
        IF (dnVec1%nb_var_Vec /= dnVec2%nb_var_Vec) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nb_var_Vec in dnVec1 and dnVec2 are different!', &
                    dnVec1%nb_var_Vec,dnVec2%nb_var_Vec
          STOP
        END IF

        IF (iVec < 1 .OR. iVec > dnVec1%nb_var_Vec) THEN ! all

          IF (nderiv_loc == 0) THEN
             dnVec2%d0 = dnVec1%d0
          ELSE IF (nderiv_loc == 1) THEN
             dnVec2%d0 = dnVec1%d0
             dnVec2%d1 = dnVec1%d1
          ELSE IF (nderiv_loc == 2) THEN
             dnVec2%d0 = dnVec1%d0
             dnVec2%d1 = dnVec1%d1
             dnVec2%d2 = dnVec1%d2
          ELSE IF (nderiv_loc == 3) THEN
             dnVec2%d0 = dnVec1%d0
             dnVec2%d1 = dnVec1%d1
             dnVec2%d2 = dnVec1%d2
             dnVec2%d3 = dnVec1%d3
          ELSE
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' nderiv > 4 is NOT possible',nderiv
            write(out_unitp,*) 'It should never append! Check the source'
            STOP
          END IF
        ELSE

          IF (nderiv_loc == 0) THEN
             dnVec2%d0(iVec)       = dnVec1%d0(iVec)
          ELSE IF (nderiv_loc == 1) THEN
             dnVec2%d0(iVec)       = dnVec1%d0(iVec)
             dnVec2%d1(iVec,:)     = dnVec1%d1(iVec,:)
          ELSE IF (nderiv_loc == 2) THEN
             dnVec2%d0(iVec)       = dnVec1%d0(iVec)
             dnVec2%d1(iVec,:)     = dnVec1%d1(iVec,:)
             dnVec2%d2(iVec,:,:)   = dnVec1%d2(iVec,:,:)
          ELSE IF (nderiv_loc == 3) THEN
             dnVec2%d0(iVec)       = dnVec1%d0(iVec)
             dnVec2%d1(iVec,:)     = dnVec1%d1(iVec,:)
             dnVec2%d2(iVec,:,:)   = dnVec1%d2(iVec,:,:)
             dnVec2%d3(iVec,:,:,:) = dnVec1%d3(iVec,:,:,:)
          ELSE
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' nderiv > 4 is NOT possible',nderiv
            write(out_unitp,*) 'It should never append! Check the source'
            STOP
          END IF

        END IF
      END SUBROUTINE sub_dnVec1_TO_dnVec2_WithIvec

      SUBROUTINE sub_dnVec1_TO_dnVec2_partial(dnVec1,dnVec2,nderiv)
        TYPE (Type_dnVec) :: dnVec1,dnVec2
        integer, optional :: nderiv

        integer :: nderiv_loc,nd
        character (len=*), parameter :: name_sub='sub_dnVec1_TO_dnVec2_partial'

        CALL check_alloc_dnVec(dnVec1,'dnVec1',name_sub)

        IF (.NOT. dnVec2%alloc) THEN
          CALL alloc_dnVec(dnVec2,dnVec1%nb_var_vec,dnVec1%nb_var_deriv,&
                                                          dnVec1%nderiv)
        END IF

        nderiv_loc = min(dnVec1%nderiv,dnVec2%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        nd = min(dnVec1%nb_var_deriv,dnVec2%nb_var_deriv)

        IF (dnVec1%nb_var_Vec /= dnVec2%nb_var_Vec) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nb_var_Vec in dnVec1 and dnVec2 are different!', &
                    dnVec1%nb_var_Vec,dnVec2%nb_var_Vec
          STOP
        END IF

          IF (nderiv_loc == 0) THEN
             dnVec2%d0(:)                = dnVec1%d0(:)
          ELSE IF (nderiv_loc == 1) THEN
             dnVec2%d0(:)                = dnVec1%d0(:)
             dnVec2%d1(:,1:nd)           = dnVec1%d1(:,1:nd)
          ELSE IF (nderiv_loc == 2) THEN
             dnVec2%d0(:)                = dnVec1%d0(:)
             dnVec2%d1(:,1:nd)           = dnVec1%d1(:,1:nd)
             dnVec2%d2(:,1:nd,1:nd)      = dnVec1%d2(:,1:nd,1:nd)
          ELSE IF (nderiv_loc == 3) THEN
             dnVec2%d0(:)                = dnVec1%d0(:)
             dnVec2%d1(:,1:nd)           = dnVec1%d1(:,1:nd)
             dnVec2%d2(:,1:nd,1:nd)      = dnVec1%d2(:,1:nd,1:nd)
             dnVec2%d3(:,1:nd,1:nd,1:nd) = dnVec1%d3(:,1:nd,1:nd,1:nd)
          ELSE
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' nderiv > 4 is NOT possible',nderiv
            write(out_unitp,*) 'It should never append! Check the source'
            STOP
          END IF

      END SUBROUTINE sub_dnVec1_TO_dnVec2_partial
      SUBROUTINE sub_dnVec1_wTO_dnVec2(dnVec1,dnVec2,w1,nderiv)
        TYPE (Type_dnVec) :: dnVec1,dnVec2
        real (kind=Rkind) :: w1
        integer, optional :: nderiv

        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='sub_dnVec1_wTO_dnVec2'

        CALL check_alloc_dnVec(dnVec1,'dnVec1',name_sub)

!        CALL check_alloc_dnVec(dnVec2,'dnVec2',name_sub)
        IF (.NOT. dnVec2%alloc) THEN
          CALL alloc_dnVec(dnVec2,dnVec1%nb_var_vec,dnVec1%nb_var_deriv,&
                                                          dnVec1%nderiv)
        END IF

        nderiv_loc = min(dnVec1%nderiv,dnVec2%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)


        IF (dnVec1%nb_var_deriv /= dnVec2%nb_var_deriv) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nb_var_deriv in dnVec1 and dnVec2 are different!',&
                    dnVec1%nb_var_deriv,dnVec2%nb_var_deriv
          STOP
        END IF
        IF (dnVec1%nb_var_Vec /= dnVec2%nb_var_Vec) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nb_var_Vec in dnVec1 and dnVec2 are different!', &
                    dnVec1%nb_var_Vec,dnVec2%nb_var_Vec
          STOP
        END IF

          IF (nderiv_loc == 0) THEN
             dnVec2%d0 = w1*dnVec1%d0
          ELSE IF (nderiv_loc == 1) THEN
             dnVec2%d0 = w1*dnVec1%d0
             dnVec2%d1 = w1*dnVec1%d1
          ELSE IF (nderiv_loc == 2) THEN
             dnVec2%d0 = w1*dnVec1%d0
             dnVec2%d1 = w1*dnVec1%d1
             dnVec2%d2 = w1*dnVec1%d2
          ELSE IF (nderiv_loc == 3) THEN
             dnVec2%d0 = w1*dnVec1%d0
             dnVec2%d1 = w1*dnVec1%d1
             dnVec2%d2 = w1*dnVec1%d2
             dnVec2%d3 = w1*dnVec1%d3
          ELSE
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' nderiv > 4 is NOT possible',nderiv
            write(out_unitp,*) 'It should never append! Check the source'
            STOP
          END IF

      END SUBROUTINE sub_dnVec1_wTO_dnVec2

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE sub_dnVec1_wADDTO_dnVec2(dnVec1,iVec1,w1,dnVec2,iVec2,w2,&
                                          ndim,nderiv)
        TYPE (Type_dnVec)   :: dnVec1,dnVec2
        integer, intent(in) :: iVec1,iVec2,ndim
        integer, optional :: nderiv
        real (kind=Rkind) :: w1,w2

        integer :: nderiv_loc,i1,f1,i2,f2
        character (len=*), parameter :: name_sub='sub_dnVec1_wADDTO_dnVec2'

        CALL check_alloc_dnVec(dnVec1,'dnVec1',name_sub)
        CALL check_alloc_dnVec(dnVec2,'dnVec2',name_sub)

        nderiv_loc = min(dnVec1%nderiv,dnVec2%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)



        IF (dnVec1%nb_var_deriv /= dnVec2%nb_var_deriv) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nb_var_deriv in dnVec1 and dnVec2 are different!',&
                    dnVec1%nb_var_deriv,dnVec2%nb_var_deriv
          STOP
        END IF

        IF (iVec1 < 1 .OR.  iVec1+ndim-1 > dnVec1%nb_var_Vec) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' iVec1 < 1 .OR. iVec1+ndim-1 > nb_var_Vec', &
                    iVec1,ndim,dnVec1%nb_var_Vec
          write(out_unitp,*) ' CHECK the fortran source!!'
          STOP
        END IF

        IF (iVec2 < 1 .OR.  iVec2+ndim-1 > dnVec2%nb_var_Vec) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' iVec2 < 1 .OR. iVec2+ndim-1 > nb_var_Vec', &
                    iVec2,ndim,dnVec2%nb_var_Vec
          write(out_unitp,*) ' CHECK the fortran source!!'
          STOP
        END IF

        i1 = iVec1
        f1 = i1 + ndim-1
        i2 = iVec2
        f2 = i2 + ndim-1

        IF (nderiv_loc == 0) THEN
         dnVec2%d0(i2:f2) = w2*dnVec2%d0(i2:f2) + w1*dnVec1%d0(i1:f1)
        ELSE IF (nderiv_loc == 1) THEN
         dnVec2%d0(i2:f2) = w2*dnVec2%d0(i2:f2) + w1*dnVec1%d0(i1:f1)
         dnVec2%d1(i2:f2,:) = w2*dnVec2%d1(i2:f2,:) + w1*dnVec1%d1(i1:f1,:)
        ELSE IF (nderiv_loc == 2) THEN
         dnVec2%d0(i2:f2) = w2*dnVec2%d0(i2:f2) + w1*dnVec1%d0(i1:f1)
         dnVec2%d1(i2:f2,:) = w2*dnVec2%d1(i2:f2,:) + w1*dnVec1%d1(i1:f1,:)
         dnVec2%d2(i2:f2,:,:) = w2*dnVec2%d2(i2:f2,:,:) + w1*dnVec1%d2(i1:f1,:,:)
        ELSE IF (nderiv_loc == 3) THEN
         dnVec2%d0(i2:f2) = w2*dnVec2%d0(i2:f2) + w1*dnVec1%d0(i1:f1)
         dnVec2%d1(i2:f2,:) = w2*dnVec2%d1(i2:f2,:) + w1*dnVec1%d1(i1:f1,:)
         dnVec2%d2(i2:f2,:,:) = w2*dnVec2%d2(i2:f2,:,:) + w1*dnVec1%d2(i1:f1,:,:)
         dnVec2%d3(i2:f2,:,:,:) = w2*dnVec2%d3(i2:f2,:,:,:) + w1*dnVec1%d3(i1:f1,:,:,:)
        ELSE
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' nderiv > 4 is NOT possible',nderiv
            write(out_unitp,*) 'It souhld never append! Check the source'
            STOP
        END IF

      END SUBROUTINE sub_dnVec1_wADDTO_dnVec2

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE dnVec2_wPLUS_dnVec3_TO_dnVec1(dnVec1,iVec1,            &
                                               dnVec2,iVec2,w2,         &
                                               dnVec3,iVec3,w3,         &
                                               ndim,nderiv)
        TYPE (Type_dnVec)   :: dnVec1,dnVec2,dnVec3
        integer, intent(in) :: iVec1,iVec2,iVec3,ndim
        integer, optional :: nderiv
        real (kind=Rkind) :: w2,w3

        integer :: nderiv_loc,i1,f1,i2,f2,i3,f3
        character (len=*), parameter :: name_sub='dnVec2_wPLUS_dnVec3_TO_dnVec1'

        CALL check_alloc_dnVec(dnVec1,'dnVec1',name_sub)
        CALL check_alloc_dnVec(dnVec2,'dnVec2',name_sub)
        CALL check_alloc_dnVec(dnVec3,'dnVec3',name_sub)

        nderiv_loc = min(dnVec1%nderiv,dnVec2%nderiv,dnVec3%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)



        IF (dnVec1%nb_var_deriv /= dnVec2%nb_var_deriv .OR.             &
            dnVec1%nb_var_deriv /= dnVec2%nb_var_deriv) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nb_var_deriv in dnVec1, dnVec2, dnVec3 are differents!',&
                    dnVec1%nb_var_deriv,dnVec2%nb_var_deriv,dnVec3%nb_var_deriv
         write(out_unitp,*) ' CHECK the fortran source!!'

          STOP
        END IF

        IF (iVec1 < 1 .OR.  iVec1+ndim-1 > dnVec1%nb_var_Vec) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' iVec1 < 1 .OR. iVec1+ndim-1 > nb_var_Vec', &
                    iVec1,ndim,dnVec1%nb_var_Vec
          write(out_unitp,*) ' CHECK the fortran source!!'
          STOP
        END IF

        IF (iVec2 < 1 .OR.  iVec2+ndim-1 > dnVec2%nb_var_Vec) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' iVec2 < 1 .OR. iVec2+ndim-1 > nb_var_Vec', &
                    iVec2,ndim,dnVec2%nb_var_Vec
          write(out_unitp,*) ' CHECK the fortran source!!'
          STOP
        END IF

        IF (iVec3 < 1 .OR.  iVec3+ndim-1 > dnVec3%nb_var_Vec) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' iVec3 < 1 .OR. iVec3+ndim-1 > nb_var_Vec', &
                    iVec3,ndim,dnVec3%nb_var_Vec
          write(out_unitp,*) ' CHECK the fortran source!!'
          STOP
        END IF

        i1 = iVec1
        f1 = i1 + ndim-1
        i2 = iVec2
        f2 = i2 + ndim-1
        i3 = iVec3
        f3 = i3 + ndim-1

        IF (nderiv_loc == 0) THEN
         dnVec1%d0(i1:f1) = w2*dnVec2%d0(i2:f2) + w3*dnVec3%d0(i3:f3)
        ELSE IF (nderiv_loc == 1) THEN
         dnVec1%d0(i1:f1) = w2*dnVec2%d0(i2:f2) + w3*dnVec3%d0(i3:f3)
         dnVec1%d1(i1:f1,:) = w2*dnVec2%d1(i2:f2,:) + w3*dnVec3%d1(i3:f3,:)
        ELSE IF (nderiv_loc == 2) THEN
         dnVec1%d0(i1:f1) = w2*dnVec2%d0(i2:f2) + w3*dnVec3%d0(i3:f3)
         dnVec1%d1(i1:f1,:) = w2*dnVec2%d1(i2:f2,:) + w3*dnVec3%d1(i3:f3,:)
         dnVec1%d2(i1:f1,:,:) = w2*dnVec2%d2(i2:f2,:,:) + w3*dnVec3%d2(i3:f3,:,:)
        ELSE IF (nderiv_loc == 3) THEN
         dnVec1%d0(i1:f1) = w2*dnVec2%d0(i2:f2) + w3*dnVec3%d0(i3:f3)
         dnVec1%d1(i1:f1,:) = w2*dnVec2%d1(i2:f2,:) + w3*dnVec3%d1(i3:f3,:)
         dnVec1%d2(i1:f1,:,:) = w2*dnVec2%d2(i2:f2,:,:) + w3*dnVec3%d2(i3:f3,:,:)
         dnVec1%d3(i1:f1,:,:,:) = w2*dnVec2%d3(i2:f2,:,:,:) + w3*dnVec3%d3(i3:f3,:,:,:)
        ELSE
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' nderiv > 4 is NOT possible',nderiv
            write(out_unitp,*) 'It souhld never append! Check the source'
            STOP
        END IF

      END SUBROUTINE dnVec2_wPLUS_dnVec3_TO_dnVec1

!================================================================
!       v2 = v1+v2
!================================================================
      SUBROUTINE sub_dnVec1_PLUS_dnVec2_TO_dnVec3(dnVec1,dnVec2,dnVec3,nderiv)
      USE mod_system
      IMPLICIT NONE

      TYPE (Type_dnVec), intent(in)    :: dnVec1,dnVec2
      TYPE (Type_dnVec), intent(inout) :: dnVec3

      integer           :: nderiv


!     -----------------------------------------------------------------
!     logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      character (len=*), parameter ::                                   &
                    name_sub= 'sub_dnVec1_PLUS_dnVec2_TO_dnVec3'
!     -----------------------------------------------------------------
      IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'dnVec1'
         CALL Write_dnVec(dnVec1)
         write(out_unitp,*) 'dnVec2'
         CALL Write_dnVec(dnVec2)
      END IF
!     -----------------------------------------------------------------

!      -----------------------------------------------------------------
       IF (nderiv == 0) THEN
         dnVec3%d0(:)       =  dnVec1%d0(:)       + dnVec2%d0(:)
       ELSE IF (nderiv == 1) THEN
         dnVec3%d0(:)       =  dnVec1%d0(:)       + dnVec2%d0(:)
         dnVec3%d1(:,:)     =  dnVec1%d1(:,:)     + dnVec2%d1(:,:)
       ELSE IF (nderiv == 2) THEN
         dnVec3%d0(:)       =  dnVec1%d0(:)       + dnVec2%d0(:)
         dnVec3%d1(:,:)     =  dnVec1%d1(:,:)     + dnVec2%d1(:,:)
         dnVec3%d2(:,:,:)   =  dnVec1%d2(:,:,:)   + dnVec2%d2(:,:,:)
       ELSE IF (nderiv .GE. 3) THEN
         dnVec3%d0(:)       =  dnVec1%d0(:)       + dnVec2%d0(:)
         dnVec3%d1(:,:)     =  dnVec1%d1(:,:)     + dnVec2%d1(:,:)
         dnVec3%d2(:,:,:)   =  dnVec1%d2(:,:,:)   + dnVec2%d2(:,:,:)
         dnVec3%d3(:,:,:,:) =  dnVec1%d3(:,:,:,:) + dnVec2%d3(:,:,:,:)
       END IF
!      -----------------------------------------------------------------

!      -----------------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'dnVec3'
         CALL Write_dnVec(dnVec3)
         write(out_unitp,*) 'END ',name_sub
       END IF
!      -----------------------------------------------------------------

      END SUBROUTINE sub_dnVec1_PLUS_dnVec2_TO_dnVec3

      SUBROUTINE sub_ZERO_TO_dnVec(dnVec,nderiv)
        TYPE (Type_dnVec) :: dnVec
        integer, optional :: nderiv

        integer           :: nderiv_loc

        CALL check_alloc_dnVec(dnVec,'dnVec','sub_ZERO_TO_dnVec')

        nderiv_loc = dnVec%nderiv
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

!       write(out_unitp,*) 'BEGINNING sub_ZERO_TO_dnVec'
!       write(out_unitp,*) 'nderiv',dnVec%nderiv
!       write(out_unitp,*) 'nb_var_deriv',dnVec%nb_var_deriv
!       write(out_unitp,*) 'nb_var_vec',dnVec%nb_var_vec

          dnVec%d0(:) = ZERO
          IF (nderiv_loc == 1) THEN
            dnVec%d1(:,:) = ZERO
          ELSE IF (nderiv_loc == 2) THEN
            dnVec%d1(:,:) = ZERO
            dnVec%d2(:,:,:) = ZERO
          ELSE IF (nderiv_loc == 3) THEN
            dnVec%d1(:,:) = ZERO
            dnVec%d2(:,:,:) = ZERO
            dnVec%d3(:,:,:,:) = ZERO
          ELSE IF (nderiv_loc > 3) THEN
            write(out_unitp,*) ' ERROR in sub_ZERO_TO_dnVec'
            write(out_unitp,*) ' nderiv_loc MUST be < 4',nderiv_loc
            STOP
          END IF

!       write(out_unitp,*) 'END sub_ZERO_TO_dnVec'

      END SUBROUTINE sub_ZERO_TO_dnVec

      SUBROUTINE test_ZERO_OF_dnVec(dnVec,nderiv)
        TYPE (Type_dnVec) :: dnVec
        integer, optional :: nderiv

        integer           :: nderiv_loc

        CALL check_alloc_dnVec(dnVec,'dnVec','test_ZERO_OF_dnVec')

        nderiv_loc = dnVec%nderiv
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

!       write(out_unitp,*) 'BEGINNING test_ZERO_OF_dnVec'
!       write(out_unitp,*) 'nderiv',dnVec%nderiv
!       write(out_unitp,*) 'nb_var_deriv',dnVec%nb_var_deriv
!       write(out_unitp,*) 'nb_var_vec',dnVec%nb_var_vec

         write(out_unitp,*) 'max(abs(...)) for ider=0',maxval(abs(dnVec%d0))
         IF (nderiv_loc >= 1) THEN
            write(out_unitp,*) 'max(abs(...)) for ider=1',maxval(abs(dnVec%d1))
         END IF
         IF (nderiv_loc >= 2) THEN
            write(out_unitp,*) 'max(abs(...)) for ider=2',maxval(abs(dnVec%d2))
         END IF
         IF (nderiv_loc >= 3) THEN
            write(out_unitp,*) 'max(abs(...)) for ider=3',maxval(abs(dnVec%d3))
         END IF

!       write(out_unitp,*) 'END test_ZERO_OF_dnVec'

      END SUBROUTINE test_ZERO_OF_dnVec

      SUBROUTINE sub_dnVec1_PROD_dnS2_TO_dnVec3(dnVec1,dnS2,dnVec3)
      use mod_dnS, only: type_dns, check_alloc_dns, write_dns
      !USE mod_system
      IMPLICIT NONE

      TYPE (Type_dnVec), intent(in)    :: dnVec1
      TYPE (Type_dnS), intent(in)      :: dnS2

      TYPE (Type_dnVec), intent(inout) :: dnVec3


      integer           :: nderiv
      integer           :: i,j,k



!     -----------------------------------------------------------------
!      logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      character (len=*), parameter :: name_sub= 'sub_dnVec1_PROD_dnS2_TO_dnVec3'
!     -----------------------------------------------------------------
      IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         CALL Write_dnS(dnS2)
         CALL Write_dnVec(dnVec1)
      END IF
!     -----------------------------------------------------------------
      CALL check_alloc_dnVec(dnVec1,'dnVec1',name_sub)
      CALL check_alloc_dnS(dnS2,'dnS2',name_sub)


      nderiv = min(dnVec1%nderiv,dnS2%nderiv)
      IF (.NOT. dnVec3%alloc) THEN
        CALL alloc_dnVec(dnVec3,dnVec1%nb_var_vec,dnVec1%nb_var_deriv,nderiv)
      END IF
      nderiv = min(dnVec3%nderiv,nderiv)

!      -----------------------------------------------------------------
!      vector
!      -----------------------------------------------------------------
       IF (nderiv == 0) THEN
         dnVec3%d0(:) = dnVec1%d0(:) * dnS2%d0

       ELSE IF (nderiv == 1) THEN

         DO i=1,dnVec1%nb_var_deriv
          dnVec3%d1(:,i) = dnVec1%d1(:,i) * dnS2%d0 +                      &
                          dnVec1%d0(:)   * dnS2%d1(i)
         END DO
         dnVec3%d0(:) = dnVec1%d0(:) * dnS2%d0
!      -----------------------------------------------------------------
       ELSE IF (nderiv == 2) THEN

         DO i=1,dnVec1%nb_var_deriv
         DO j=1,dnVec1%nb_var_deriv
          dnVec3%d2(:,i,j) = dnVec1%d2(:,i,j) * dnS2%d0    +               &
                            dnVec1%d1(:,i)   * dnS2%d1(j) +               &
                            dnVec1%d1(:,j)   * dnS2%d1(i) +               &
                            dnVec1%d0(:)     * dnS2%d2(i,j)
         END DO
         END DO
         DO i=1,dnVec1%nb_var_deriv
          dnVec3%d1(:,i) = dnVec1%d1(:,i) * dnS2%d0 +                      &
                          dnVec1%d0(:)   * dnS2%d1(i)
         END DO
         dnVec3%d0(:) = dnVec1%d0(:) * dnS2%d0

!      -----------------------------------------------------------------
       ELSE IF (nderiv == 3) THEN

         DO i=1,dnVec1%nb_var_deriv
         DO j=1,dnVec1%nb_var_deriv
         DO k=1,dnVec1%nb_var_deriv
          dnVec3%d3(:,i,j,k) = dnVec1%d3(:,i,j,k) * dnS2%d0      +         &
                              dnVec1%d2(:,i,j)   * dnS2%d1(k)   +         &
                              dnVec1%d2(:,i,k)   * dnS2%d1(j)   +         &
                              dnVec1%d2(:,j,k)   * dnS2%d1(i)   +         &
                              dnVec1%d1(:,i)     * dnS2%d2(j,k) +         &
                              dnVec1%d1(:,j)     * dnS2%d2(i,k) +         &
                              dnVec1%d1(:,k)     * dnS2%d2(i,j) +         &
                              dnVec1%d0(:)       * dnS2%d3(i,j,k)
         END DO
         END DO
         END DO

         DO i=1,dnVec1%nb_var_deriv
         DO j=1,dnVec1%nb_var_deriv
          dnVec3%d2(:,i,j) = dnVec1%d2(:,i,j) * dnS2%d0    +               &
                            dnVec1%d1(:,i)   * dnS2%d1(j) +               &
                            dnVec1%d1(:,j)   * dnS2%d1(i) +               &
                            dnVec1%d0(:)     * dnS2%d2(i,j)
         END DO
         END DO
         DO i=1,dnVec1%nb_var_deriv
          dnVec3%d1(:,i) = dnVec1%d1(:,i) * dnS2%d0 +                      &
                          dnVec1%d0(:)   * dnS2%d1(i)
         END DO
         dnVec3%d0(:) = dnVec1%d0(:) * dnS2%d0

       ELSE
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nderiv MUST be <4',nderiv
         STOP
       END IF


!      -----------------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'new vector'
         CALL Write_dnVec(dnVec3)
         write(out_unitp,*) 'END ',name_sub
       END IF
!      -----------------------------------------------------------------

       END SUBROUTINE sub_dnVec1_PROD_dnS2_TO_dnVec3

!================================================================
!       norm of a vector d0v
!       then normalization of d0v
!================================================================
      SUBROUTINE sub_Normalize_dnVec(dnVec)
      use mod_dnS, only: type_dns, alloc_dns, check_alloc_dns, write_dns, sub_dns1_to_dntr2, dealloc_dns
      !USE mod_system
      IMPLICIT NONE

      TYPE (Type_dnVec), intent(inout) :: dnVec


      TYPE (Type_dnS)   :: dnNorm,dnNorm_inv
      integer           :: nderiv
      real (kind=Rkind) :: cte(20)

!     -----------------------------------------------------------------
      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      character (len=*), parameter :: name_sub= 'sub_Normalize_dnVec'
!      -----------------------------------------------------------------
      IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         CALL Write_dnVec(dnVec)
      END IF

      nderiv = dnVec%nderiv

      CALL alloc_dnS(dnNorm,    dnVec%nb_var_deriv,nderiv)
      CALL alloc_dnS(dnNorm_inv,dnVec%nb_var_deriv,nderiv)
      cte(:) = ZERO
      cte(1) = -HALF

      CALL sub_dot_product_dnVec1_dnVec2_TO_dnS(dnVec,dnVec,dnNorm)
      IF (dnNorm%d0 /= ZERO) THEN
        CALL sub_dnS1_TO_dntR2(dnNorm,dnNorm_inv,transfo_1D=99,cte=cte) ! 1/sqrt(Norm)
      ELSE
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' The norm of the vector is zero'
        CALL Write_dnS(dnNorm)
        STOP
      END IF

       !-----------------------------------------------------------------
       !then normalization of d0v
       CALL sub_dnVec1_PROD_dnS2_TO_dnVec3(dnVec,dnNorm_inv,dnVec)

       !-----------------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'dnNorm_inv ='
         CALL Write_dnS(dnNorm_inv)
         write(out_unitp,*) 'Normalized vector'
         CALL Write_dnVec(dnVec)
         write(out_unitp,*) 'END ',name_sub
       END IF

       CALL dealloc_dnS(dnNorm)
       CALL dealloc_dnS(dnNorm_inv)

       END SUBROUTINE sub_Normalize_dnVec

      ! composition dnVec3=dnVec2(dnVec1)
      ! dnVec3%d0(i23) = dnVec2%d0(i23)
      !   dnVec3%d0(i23,i31) =  ....
      !   d_dnVec3/d_Qi1der = Sum_(i2der or ) d_Qi2der/d_  * d_dnVec2/d_Qi2der
      SUBROUTINE dnVec2_O_dnVec1_TO_dnVec3(dnVec1,dnVec2,dnVec3,nderiv)
        TYPE (Type_dnVec), intent(in)    :: dnVec1,dnVec2
        TYPE (Type_dnVec), intent(inout) :: dnVec3
        integer, optional :: nderiv

        integer :: nderiv_loc,i3,i2,i1,j1,k1,i0,j0,k0

!     -----------------------------------------------------------------
      integer, parameter :: nderiv_debug=1
      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      character (len=*), parameter :: name_sub= 'dnVec2_O_dnVec1_TO_dnVec3'
!      -----------------------------------------------------------------
      IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'dnVec1, dnVec2'
         CALL Write_dnVec(dnVec1,nderiv_debug)
         CALL Write_dnVec(dnVec2,nderiv_debug)
      END IF

        CALL check_alloc_dnVec(dnVec1,'dnVec1',name_sub)
        CALL check_alloc_dnVec(dnVec2,'dnVec2',name_sub)
        CALL check_alloc_dnVec(dnVec3,'dnVec3',name_sub)

        nderiv_loc = min(dnVec1%nderiv,dnVec2%nderiv,dnVec3%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)
        IF (debug) write(out_unitp,*) 'nderiv',nderiv_loc


        IF (dnVec3%nb_var_vec   /= dnVec2%nb_var_vec   .OR.  &
            dnVec3%nb_var_deriv /= dnVec1%nb_var_deriv .OR.  &
            dnVec2%nb_var_deriv /= dnVec1%nb_var_vec  ) THEN

          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' WRONG nb_var_deriv nb_var_vec parameters'
          write(out_unitp,*) '  dnVec3%nb_var_vec   = dnVec2%nb_var_vec?  ', &
                                          dnVec3%nb_var_vec,dnVec2%nb_var_vec
          write(out_unitp,*) '  dnVec3%nb_var_deriv = dnVec1%nb_var_deriv?', &
                                      dnVec3%nb_var_deriv,dnVec1%nb_var_deriv
          write(out_unitp,*) '  dnVec2%nb_var_deriv = dnVec1%nb_var_vec?  ', &
                                        dnVec2%nb_var_deriv,dnVec1%nb_var_vec
          STOP
        END IF

        dnVec3%d0(:) = dnVec2%d0(:)

        IF (nderiv_loc >= 1) THEN
          DO i3=1,dnVec3%nb_var_vec
          DO i0=1,dnVec1%nb_var_deriv
            dnVec3%d1(i3,i0) = dot_product(dnVec1%d1(:,i0),dnVec2%d1(i3,:))
          END DO
          END DO
        END IF

        IF (nderiv_loc >= 2) THEN
          DO i3=1,dnVec3%nb_var_vec
          DO i0=1,dnVec1%nb_var_deriv
          DO j0=1,dnVec1%nb_var_deriv

            dnVec3%d2(i3,j0,i0) = dot_product(dnVec1%d2(:,j0,i0),dnVec2%d1(i3,:))

            DO i1=1,dnVec1%nb_var_vec
            DO j1=1,dnVec1%nb_var_vec
              dnVec3%d2(i3,j0,i0) = dnVec3%d2(i3,j0,i0) +               &
                  dnVec1%d1(i1,i0) * dnVec1%d1(j1,j0) * dnVec2%d2(i3,j1,i1)
            END DO
            END DO

          END DO
          END DO
          END DO
        END IF

        IF (nderiv_loc == 3) THEN

          DO i3=1,dnVec3%nb_var_vec
          DO i0=1,dnVec1%nb_var_deriv
          DO j0=1,dnVec1%nb_var_deriv
          DO k0=1,dnVec1%nb_var_deriv

            dnVec3%d3(i3,k0,j0,i0) = dot_product(dnVec1%d3(:,k0,j0,i0),dnVec2%d1(i3,:))


            DO i1=1,dnVec1%nb_var_vec
            DO j1=1,dnVec1%nb_var_vec
              dnVec3%d3(i3,k0,j0,i0) = dnVec3%d3(i3,k0,j0,i0) +         &
                  dnVec1%d2(i1,i0,j0) * dnVec1%d1(j1,k0) * dnVec2%d2(i3,j1,i1) + &
                  dnVec1%d2(i1,k0,i0) * dnVec1%d1(j1,j0) * dnVec2%d2(i3,j1,i1) + &
                  dnVec1%d2(i1,j0,k0) * dnVec1%d1(j1,i0) * dnVec2%d2(i3,j1,i1)
            END DO
            END DO

            DO i1=1,dnVec1%nb_var_vec
            DO j1=1,dnVec1%nb_var_vec
            DO k1=1,dnVec1%nb_var_vec
              dnVec3%d3(i3,k0,j0,i0) = dnVec3%d3(i3,k0,j0,i0) +               &
                  dnVec1%d1(i1,i0) * dnVec1%d1(j1,j0) * dnVec1%d1(k1,k0) * dnVec2%d3(i3,k1,j1,i1)
            END DO
            END DO
            END DO

          END DO
          END DO
          END DO
          END DO
        END IF


        IF (nderiv_loc > 3) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' nderiv > 4 is NOT possible',nderiv
            write(out_unitp,*) 'It should never append! Check the source'
            STOP
        END IF

      IF (debug) THEN
         write(out_unitp,*) 'dnVec3'
         CALL Write_dnVec(dnVec3,nderiv_debug)
         write(out_unitp,*) 'END ',name_sub
      END IF


      END SUBROUTINE dnVec2_O_dnVec1_TO_dnVec3

END MODULE mod_dnV

