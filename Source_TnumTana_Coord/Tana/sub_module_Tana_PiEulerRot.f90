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

 module mod_Tana_PiEulerRot
   !! @description: This module defines the data structures. It contains also
   !!               the some standard routine that initilize, 
   !!               allocate and delete the data structure
 USE mod_system
 USE mod_Tana_OpEl
 USE mod_Tana_Op1D
 USE mod_Tana_OpnD
 USE mod_Tana_sum_opnd
 USE mod_Tana_VecSumOpnD
 IMPLICIT NONE
! PRIVATE
 character (len=*), parameter, private :: mod_name = 'mod_Tana_PiEulerRot'

        !!@description: TODO
        !!@param: TODO
        TYPE Type_PiEulerRot
          type(vec_sum_opnd)           :: Pi                         ! conjugate  momentum
          type(vec_sum_opnd)           :: Pidag                      ! adjoint of Pi
          integer, pointer             :: Tab_num_Frame(:) => null() ! table of the frame numbers
          integer                      :: i_project=0                ! index of the projection
        END TYPE Type_PiEulerRot

      INTERFACE alloc_array
        MODULE PROCEDURE alloc_array_OF_PiEulerRotdim1
      END INTERFACE
      INTERFACE dealloc_array
        MODULE PROCEDURE dealloc_array_OF_PiEulerRotdim1
      END INTERFACE

   CONTAINS

      SUBROUTINE alloc_array_OF_PiEulerRotdim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      type(Type_PiEulerRot), pointer, intent(out) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=1
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc='alloc_array_OF_PiEulerRotdim1'
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
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'PiEulerRot')

      END SUBROUTINE alloc_array_OF_PiEulerRotdim1
      SUBROUTINE dealloc_array_OF_PiEulerRotdim1(tab,name_var,name_sub)
      IMPLICIT NONE

      type(Type_PiEulerRot), pointer, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

      integer :: i
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_PiEulerRotdim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN
       IF (.NOT. associated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       DO i=lbound(tab,dim=1),ubound(tab,dim=1)
        call delete_op(tab(i)%Pi)
        call delete_op(tab(i)%Pidag)
        IF (associated(tab(i)%Tab_num_Frame) )                          &
          CALL dealloc_array(tab(i)%Tab_num_Frame,'Tab_num_Frame',name_sub)
        tab(i)%i_project = 0                ! index of the projection
       END DO

       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'PiEulerRot')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_PiEulerRotdim1

      SUBROUTINE alloc_NParray_OF_PiEulerRotdim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      type(Type_PiEulerRot), allocatable, intent(out) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=1
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc='alloc_NParray_OF_PiEulerRotdim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       CALL dealloc_NParray_OF_PiEulerRotdim1(tab,name_var,name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)

         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1)),stat=err_mem)
       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'PiEulerRot')

      END SUBROUTINE alloc_NParray_OF_PiEulerRotdim1
      SUBROUTINE dealloc_NParray_OF_PiEulerRotdim1(tab,name_var,name_sub)
      IMPLICIT NONE

      type(Type_PiEulerRot), allocatable, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

      integer :: i
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_NParray_OF_PiEulerRotdim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       IF (.NOT. allocated(tab)) RETURN

       DO i=lbound(tab,dim=1),ubound(tab,dim=1)
        call delete_op(tab(i)%Pi)
        call delete_op(tab(i)%Pidag)

        IF (associated(tab(i)%Tab_num_Frame) )                          &
          CALL dealloc_array(tab(i)%Tab_num_Frame,'Tab_num_Frame',name_sub)

        tab(i)%i_project = 0                ! index of the projection

       END DO

       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'PiEulerRot')

      END SUBROUTINE dealloc_NParray_OF_PiEulerRotdim1

 !! @description: Deallocated a variable of type Type_PiEulerRot
 !! @param:    P_euler    The 1d operator (type: Type_PiEulerRot). 
 subroutine delete_P_euler(P_euler)

   type(Type_PiEulerRot),            intent(inout)    :: P_euler(:)

   character (len=*), parameter :: routine_name="delete_P_euler"

   integer                    :: i

    do i=1, size(P_euler)
      call delete_op(P_euler(i)%Pi)
      call delete_op(P_euler(i)%Pidag)
      CALL dealloc_array(P_euler(i)%Tab_num_Frame,'P_euler(i)%Tab_num_Frame',routine_name)
      P_euler(i)%i_project = 0                ! index of the projection
    end do
 end subroutine delete_P_euler



 end module mod_Tana_PiEulerRot
