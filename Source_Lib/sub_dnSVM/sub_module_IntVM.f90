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
MODULE mod_IntVM
      use mod_system, only: alloc_array, out_unitp, dealloc_array, &
                            write_error_not_null, sub_test_tab_ub, &
                            sub_test_tab_lb, error_memo_allo, write_error_null
      IMPLICIT NONE

      PRIVATE

      !==============================================
      !!@description: TODO
      !!@param: TODO
     TYPE Type_IntVec
          logical                     :: alloc=.FALSE.

          integer                     :: nb_var_vec   = 0

          integer, pointer  :: vec(:)     => null()

      END TYPE Type_IntVec
     TYPE Type_IntMat
          logical                     :: alloc=.FALSE.

          integer                     :: nb_var_Matl  = 0
          integer                     :: nb_var_Matc  = 0

          integer, pointer  :: mat(:,:)     => null()

      END TYPE Type_IntMat


      INTERFACE alloc_array
        MODULE PROCEDURE alloc_array_OF_IntVecdim1
        MODULE PROCEDURE alloc_array_OF_IntMatdim1
      END INTERFACE
      INTERFACE dealloc_array
        MODULE PROCEDURE dealloc_array_OF_IntVecdim1
        MODULE PROCEDURE dealloc_array_OF_IntMatdim1
      END INTERFACE

      PUBLIC :: Type_IntVec, alloc_IntVec, dealloc_IntVec, check_alloc_IntVec, &
                Write_IntVec, sub_IntVec1_TO_IntVec2, sub_ZERO_TO_IntVec
      PUBLIC :: Type_IntMat, alloc_IntMat, dealloc_IntMat, check_alloc_IntMat, &
                Write_IntMat, sub_IntMat1_TO_IntMat2
      PUBLIC :: alloc_array, dealloc_array

      CONTAINS
!
!================================================================
!
!     allocation
!
!================================================================
!
      !!@description: TODO
      !!@param: TODO
      SUBROUTINE alloc_IntVec(IntVec,nb_var_vec)
        TYPE (Type_IntVec) :: IntVec
        integer, optional :: nb_var_vec
        integer :: err_mem

        IF (present(nb_var_vec)) IntVec%nb_var_vec = nb_var_vec


        IF (IntVec%alloc) RETURN
        IntVec%alloc = .TRUE.

!        write(out_unitp,*) 'BEGINNING alloc_IntVec'
!        write(out_unitp,*) 'nb_var_vec',IntVec%nb_var_vec

        IF (IntVec%nb_var_vec > 0) THEN
          CALL alloc_array(IntVec%vec,(/ IntVec%nb_var_vec /),          &
                          'IntVec%vec','alloc_IntVec')
          IntVec%vec(:) = 0

        ELSE
          write(out_unitp,*) ' ERROR in alloc_IntVec'
          write(out_unitp,*) ' nb_var_vec MUST be > 0',IntVec%nb_var_vec
          STOP
        END IF

!        write(out_unitp,*) 'END alloc_IntVec'

      END SUBROUTINE alloc_IntVec

      SUBROUTINE alloc_IntMat(IntMat,nb_var_Matl,nb_var_Matc)
        TYPE (Type_IntMat) :: IntMat
        integer, optional :: nb_var_Matl,nb_var_Matc
        integer :: err_mem,memory

        IF (present(nb_var_Matl)) IntMat%nb_var_Matl = nb_var_Matl
        IF (present(nb_var_Matc)) IntMat%nb_var_Matc = nb_var_Matc


        IF (IntMat%alloc) RETURN
        IntMat%alloc = .TRUE.

        IF (IntMat%nb_var_Matc == 0) IntMat%nb_var_Matc = IntMat%nb_var_Matl

!        write(out_unitp,*) 'BEGINNING alloc_IntMat'
!        write(out_unitp,*) 'nb_var_Matl,nb_var_Matc',IntMat%nb_var_Matl,IntMat%nb_var_Matc

        IF (IntMat%nb_var_Matl > 0) THEN
          CALL alloc_array(IntMat%mat,                                  &
                           (/ IntMat%nb_var_Matl,IntMat%nb_var_Matc /), &
                          'IntMat%mat','alloc_IntMat')
          IntMat%mat(:,:) = 0
        ELSE
          write(out_unitp,*) ' ERROR in alloc_IntMat'
          write(out_unitp,*) ' nb_var_Matl and nb_var_Matc MUST be > 0',IntMat%nb_var_Matl,IntMat%nb_var_Matc
          STOP
        END IF

!        write(out_unitp,*) 'END alloc_IntMat'

      END SUBROUTINE alloc_IntMat

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE dealloc_IntVec(IntVec)
        TYPE (Type_IntVec) :: IntVec

        IF (associated(IntVec%vec)) THEN
          CALL dealloc_array(IntVec%vec,'IntVec%vec','dealloc_IntVec')
        END IF

        IntVec%alloc    = .FALSE.

        IntVec%nb_var_vec  = 0

      END SUBROUTINE dealloc_IntVec
      SUBROUTINE dealloc_IntMat(IntMat)
        TYPE (Type_IntMat) :: IntMat

!        write(out_unitp,*) 'BEGINNING alloc_IntMat'
!        write(out_unitp,*) 'nb_var_Matl,nb_var_Matc',IntMat%nb_var_Matl,IntMat%nb_var_Matc

        IF (associated(IntMat%mat)) THEN
          CALL dealloc_array(IntMat%mat,'IntMat%mat','dealloc_IntMat')
        END IF
        IntMat%nb_var_Matl = 0
        IntMat%nb_var_Matc = 0
        IntMat%alloc       = .FALSE.

!        write(out_unitp,*) 'END dealloc_IntMat'

      END SUBROUTINE dealloc_IntMat

      SUBROUTINE alloc_array_OF_IntVecdim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      TYPE (Type_IntVec), pointer, intent(inout) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=1
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_IntVecdim1'
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
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Type_IntVec')

      END SUBROUTINE alloc_array_OF_IntVecdim1
      SUBROUTINE dealloc_array_OF_IntVecdim1(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (Type_IntVec), pointer, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

      integer :: i
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_IntVecdim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN
       IF (.NOT. associated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       DO i=lbound(tab,dim=1),ubound(tab,dim=1)
         CALL dealloc_IntVec(tab(i))
       END DO

       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Type_IntVec')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_IntVecdim1
      SUBROUTINE alloc_array_OF_IntMatdim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      TYPE (Type_IntMat), pointer, intent(inout) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=1
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_IntMatdim1'
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
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Type_IntMat')

      END SUBROUTINE alloc_array_OF_IntMatdim1
      SUBROUTINE dealloc_array_OF_IntMatdim1(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (Type_IntMat), pointer, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

      integer :: i
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_IntMatdim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN
       IF (.NOT. associated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       DO i=lbound(tab,dim=1),ubound(tab,dim=1)
         CALL dealloc_IntMat(tab(i))
       END DO

       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Type_IntMat')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_IntMatdim1
!================================================================
!
!     check if alloc has been done
!
!================================================================
      !!@description: heck if alloc has been done
      !!@param: TODO
      SUBROUTINE check_alloc_IntVec(A,name_A,name_sub)
        TYPE (Type_IntVec), intent(in) :: A
        character (len=*), intent(in) :: name_A
        character (len=*), intent(in) :: name_sub

        IF ( .NOT. A%alloc) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) name_A,' has NOT been initiated with "alloc_IntVec"'
          write(out_unitp,*) ' CHECK the source!!!!!'
          STOP
        END IF
      END SUBROUTINE check_alloc_IntVec
      SUBROUTINE check_alloc_IntMat(A,name_A,name_sub)
        TYPE (Type_IntMat), intent(in) :: A
        character (len=*), intent(in) :: name_A
        character (len=*), intent(in) :: name_sub

        IF ( .NOT. A%alloc) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) name_A,' has NOT been initiated with "alloc_IntMat"'
          write(out_unitp,*) ' CHECK the source!!!!!'
          STOP
        END IF
      END SUBROUTINE check_alloc_IntMat


      SUBROUTINE Write_IntVec(IntVec)
        TYPE (Type_IntVec) :: IntVec

        CALL check_alloc_IntVec(IntVec,'IntVec','write_IntVec')


        write(out_unitp,*) 'BEGINNING Write IntVec'
        write(out_unitp,*) 'nb_var_vec',IntVec%nb_var_vec

        IF (associated(IntVec%vec)) write(out_unitp,*) IntVec%vec

        write(out_unitp,*) 'END Write IntVec'


      END SUBROUTINE Write_IntVec

      SUBROUTINE Write_IntMat(IntMat)
        TYPE (Type_IntMat) :: IntMat

        CALL check_alloc_IntMat(IntMat,'IntMat','write_IntMat')


        write(out_unitp,*) 'BEGINNING Write IntMat'
        write(out_unitp,*) 'nb_var_Matl,nb_var_Matc',                   &
                                   IntMat%nb_var_Matl,IntMat%nb_var_Matc

        IF (associated(IntMat%mat)) THEN
          CALL ecriture_intr4(IntMat%mat,out_unitp,IntMat%nb_var_Matl,  &
                                                   IntMat%nb_var_Matc,5)
        END IF

        write(out_unitp,*) 'END Write IntMat'


      END SUBROUTINE Write_IntMat

      SUBROUTINE sub_IntVec1_TO_IntVec2(IntVec1,IntVec2)
        TYPE (Type_IntVec), intent(in)    :: IntVec1
        TYPE (Type_IntVec), intent(inout) :: IntVec2

        character (len=*), parameter :: name_sub='sub_IntVec1_TO_IntVec2'

        IF (.NOT. IntVec1%alloc) RETURN

        IF (IntVec2%alloc .OR. IntVec1%nb_var_Vec /= IntVec2%nb_var_Vec) THEN
          CALL dealloc_IntVec(IntVec2)
        END IF
        CALL alloc_IntVec(IntVec2,IntVec1%nb_var_Vec)

        IntVec2%vec(:) = IntVec1%vec(:)
        IntVec2%alloc = .TRUE.

      END SUBROUTINE sub_IntVec1_TO_IntVec2

      SUBROUTINE sub_IntMat1_TO_IntMat2(IntMat1,IntMat2)
        TYPE (Type_IntMat), intent(in)    :: IntMat1
        TYPE (Type_IntMat), intent(inout) :: IntMat2


        character (len=*), parameter :: name_sub='sub_IntMat1_TO_IntMat2'

        IF (.NOT. IntMat1%alloc) RETURN

        IF (IntMat2%alloc) CALL dealloc_IntMat(IntMat2)
        CALL alloc_IntMat(IntMat2,IntMat1%nb_var_Matl,IntMat1%nb_var_Matc)


        IntMat2%Mat(:,:) = IntMat1%Mat(:,:)
        IntMat2%alloc = .TRUE.

      END SUBROUTINE sub_IntMat1_TO_IntMat2

      SUBROUTINE sub_ZERO_TO_IntVec(IntVec)
        TYPE (Type_IntVec) :: IntVec

        CALL check_alloc_IntVec(IntVec,'IntVec','sub_ZERO_TO_IntVec')

!       write(out_unitp,*) 'BEGINNING sub_ZERO_TO_IntVec'
!       write(out_unitp,*) 'nb_var_vec',IntVec%nb_var_vec

        IntVec%vec(:) = 0

!       write(out_unitp,*) 'END sub_ZERO_TO_IntVec'

      END SUBROUTINE sub_ZERO_TO_IntVec

END MODULE mod_IntVM

