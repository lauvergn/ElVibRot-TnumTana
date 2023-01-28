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
!      Tnum is written David Lauvergnat [1]
!      Tana is written by Mamadou Ndong [1] and David Lauvergnat [1]
!         with contributions
!          Emil Lund klinting (coupling with MidasCpp) [3]'
!
![1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France
![3]: Department of Chemistry, Aarhus University, DK-8000 Aarhus C, Denmark
!
!===========================================================================
!===========================================================================
 module mod_Tana_VecSumOpnD
 use mod_system
 USE mod_Tana_sum_opnd
 IMPLICIT NONE

 PRIVATE
      !-----------------------------------------------------------!
      !                        VEC_SUM_OPND                       !
      !-----------------------------------------------------------!
      !! @description: Definition of a type of a an array of sum of nd-operators.
      !!               This type is used for the analytical
      !!               computation of the KEO
      !! @param: vec_sum             Array of sum_opnd operators
      TYPE vec_sum_opnd
        type(sum_opnd), allocatable        :: vec_sum(:) ! Array of opnd
      END TYPE vec_sum_opnd

      INTERFACE alloc_array
        MODULE PROCEDURE alloc_array_OF_Vec_Sum_OpnDdim1
      END INTERFACE
      INTERFACE dealloc_array
        MODULE PROCEDURE dealloc_array_OF_Vec_Sum_OpnDdim1
      END INTERFACE

      INTERFACE alloc_NParray
        MODULE PROCEDURE alloc_NParray_OF_Vec_Sum_OpnDdim1
      END INTERFACE
      INTERFACE dealloc_NParray
        MODULE PROCEDURE dealloc_NParray_OF_Vec_Sum_OpnDdim1
      END INTERFACE

  !!@description: Generic routine that deletes an array of operator types
  interface delete_op
    module procedure delete_vec_sum_opnd
  end interface

  !!@description: Generic routine that allocate variables of operator type
  interface allocate_op
    module procedure allocate_vec_sum_opnd
  end interface

   INTERFACE write_op
     module procedure  write_vec_sum_opnd
   END INTERFACE

  !!@description: Generic routine that copy a operator F1 to another operator F2
  interface copy_F1_into_F2
    module procedure  copy_V1_sum_nd_into_V2_sum_nd
  end interface

  PUBLIC :: vec_sum_opnd, allocate_op, delete_op, write_op, copy_F1_into_F2
  PUBLIC :: alloc_array, dealloc_array, alloc_NParray, dealloc_NParray
  PUBLIC :: V1_scalar_V2_in_F_sum_nd
  PUBLIC :: V1_cross_V2_in_Vres
  PUBLIC :: V1_plus_V2_in_Vres,V1_PLUS_TO_Vres
  PUBLIC :: V1_MINUS_TO_Vres
  PUBLIC :: M1_times_M2_in_Mres, M_opnd_times_V_in_Vres, V_times_M_opnd_in_Vres
  PUBLIC :: F_sum_nd_times_V_in_Vres
  PUBLIC :: zero_TO_vec_sum_opnd

  contains

      SUBROUTINE alloc_array_OF_Vec_Sum_OpnDdim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      type(vec_sum_opnd), pointer, intent(inout) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=1
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_Vec_Sum_OpnDdim1'
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
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Vec_Sum_OpnD')

      END SUBROUTINE alloc_array_OF_Vec_Sum_OpnDdim1
      SUBROUTINE dealloc_array_OF_Vec_Sum_OpnDdim1(tab,name_var,name_sub)
      IMPLICIT NONE

      type(vec_sum_opnd), pointer, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

      integer :: i
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_Vec_Sum_OpnDdim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN
       IF (.NOT. associated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       DO i=lbound(tab,dim=1),ubound(tab,dim=1)
         CALL delete_op(tab(i))
       END DO

       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Vec_Sum_OpnD')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_Vec_Sum_OpnDdim1

      SUBROUTINE alloc_NParray_OF_Vec_Sum_OpnDdim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      type(vec_sum_opnd), allocatable, intent(inout) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=1
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_NParray_OF_Vec_Sum_OpnDdim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       CALL dealloc_NParray_OF_Vec_Sum_OpnDdim1(tab,name_var,name_sub)


       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)

         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1)),stat=err_mem)
       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Vec_Sum_OpnD')

      END SUBROUTINE alloc_NParray_OF_Vec_Sum_OpnDdim1
      SUBROUTINE dealloc_NParray_OF_Vec_Sum_OpnDdim1(tab,name_var,name_sub)
      IMPLICIT NONE

      type(vec_sum_opnd), allocatable, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

      integer :: i
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_NParray_OF_Vec_Sum_OpnDdim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       IF (.NOT. allocated(tab)) RETURN

       DO i=lbound(tab,dim=1),ubound(tab,dim=1)
         CALL delete_op(tab(i))
       END DO

       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Vec_Sum_OpnD')

      END SUBROUTINE dealloc_NParray_OF_Vec_Sum_OpnDdim1

 !! @description: Deallocated a nd operator
 !! @param:    F_vec    The nd operator (type: vec_sum_opnd).
 subroutine delete_vec_sum_opnd(F_vec)

   type(vec_sum_opnd),           intent(inout)    :: F_vec

   character (len=*), parameter :: routine_name="delete_vec_sum_opnd"

   CALL dealloc_NParray(F_vec%vec_sum,'F_vec%vec_sum',routine_name)

 end subroutine delete_vec_sum_opnd

 !! @description: Allocated vect_opnd operator
 !! @param:    F_vec_nd    The nd operator (type: vec_sum_opnd).
 !! @param:       ndim     Size of F_vec%vec_sum.
 subroutine allocate_vec_sum_opnd(F_vec, ndim)

   type(vec_sum_opnd),           intent(inout)    :: F_vec
   integer,                      intent(in)       :: ndim

   character (len=*), parameter :: routine_name="allocate_vec_sum_opnd"

   CALL dealloc_NParray(F_vec%vec_sum,'F_vec%vec_sum',routine_name)

   CALL alloc_NParray(F_vec%vec_sum,[ndim],'F_vec%vec_sum',routine_name)

 end subroutine allocate_vec_sum_opnd


 subroutine zero_TO_vec_sum_opnd(F_vec)

   type(vec_sum_opnd),           intent(inout)    :: F_vec

   integer :: i
   character (len=*), parameter :: routine_name="zero_TO_vec_sum_opnd"

   if(allocated(F_vec%vec_sum)) then
     DO i=1,size(F_vec%vec_sum)
       F_vec%vec_sum(i) = czero
     END DO
   end if

 end subroutine zero_TO_vec_sum_opnd
   !! @description: Write an array of sum of nd operators,
   !! @param:       V_sum_nd      The operator (type: vec_sum_opnd).
   !! @param:       header    If present write comment in the biginning of the
   !!                         of the file
   SUBROUTINE write_vec_sum_opnd(V_sum_nd, i_file, header, append, close_file)
     type(vec_sum_opnd),        intent(in)       :: V_sum_nd
     integer, optional,         intent(in)       :: i_file
     logical, optional,         intent(in)       :: header
     logical, optional,         intent(in)       :: append
     logical, optional,         intent(in)       :: close_file

     integer                        :: error
     integer                        :: i_open
     integer                        :: i, j
     logical                        :: header_loc

   character (len=*), parameter :: routine_name='write_vec_sum_opnd'

     IF ( .NOT. allocated(V_sum_nd%vec_sum) ) RETURN

     if(present(i_file)) then
       i_open = i_file
     else
       i_open = out_unitp
     end if

     if (present(header)) then
       header_loc = header
     else
       header_loc = .FALSE.
     end if

     if (header_loc) then
       write(i_open, '(A)')  "   idf          iqd        alfa       &
            &indexq          op_name                 qval               coeff&
            &           opval"
       else
       write(i_open,*) '@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@&
                  &@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@'
       write(i_open,*)
       write(i_open, *)
       write(i_open, *)
       write(i_open, '(40x, A)') '========== writing a vector of sum of nd operators ======== '
       write(i_open, *)
     end if
     do i = 1, size(V_sum_nd%vec_sum)
       write(i_open, '(40x, A, 1x, I2)') '============vector: component ============', i
       call write_op(V_sum_nd%vec_sum(i), i_open)
       write(i_open, *)
       write(i_open, *)
     end do
   END SUBROUTINE write_vec_sum_opnd

 !! @description: Copy a vector of operator sum_opnd into
 !!               another vector of operators sum_opnd
 !! @param:     V1_sum_nd    The operator which will be copied
 !! @param:    V2_sum_nd    The operator in which F1_sum_nd will be copied
 subroutine copy_V1_sum_nd_into_V2_sum_nd(V1_sum_nd, V2_sum_nd)

   type(vec_sum_opnd),           intent(in)    :: V1_sum_nd
   type(vec_sum_opnd),           intent(inout) :: V2_sum_nd

   integer                    :: i

   character (len=*), parameter :: routine_name="copy_V1_sum_nd_into_V2_sum_nd"

   call allocate_vec_sum_opnd(V2_sum_nd,size(V1_sum_nd%vec_sum))
   do i = 1, size(V1_sum_nd%vec_sum)
     call copy_F1_into_F2(V1_sum_nd%vec_sum(i),V2_sum_nd%vec_sum(i))
   end do
 end subroutine copy_V1_sum_nd_into_V2_sum_nd

   !! @description: Calculates the scalar product of V1 and V2
   !! @param:       V1         The first vector (type: vec_sum_opnd)
   !! @param:       V2         The second vector (type: vec_sum_opnd)
   !! @param:       F_sum_nd   The result (type: sum_opnd)
   SUBROUTINE V1_scalar_V2_in_F_sum_nd(V1, V2, F_sum_nd)
     type(vec_sum_opnd),     intent(in)        :: V1
     type(vec_sum_opnd),     intent(in)        :: V2
     type(sum_opnd),         intent(inout)     :: F_sum_nd

     type(sum_opnd)                  :: sum_tmp1
     type(sum_opnd)                  :: sum_tmp2
     integer                         :: i
     integer                         :: ndim1, ndim2

     character (len=*), parameter :: routine_name='V1_scalar_V2_in_F_sum_nd'

     if(.not.allocated(V1%vec_sum) .or. .not.allocated(V2%vec_sum)) then
       write(out_unitp,*) ' ERROR in ',routine_name
       write(out_unitp,*) "vectors V1 and V2 should be allocated"
       STOP
     end if
     ndim1 = size(V1%vec_sum)
     ndim2 = size(V2%vec_sum)
     if(ndim1 /= ndim2) then
       write(out_unitp,*) ' ERROR in ',routine_name
       write(out_unitp,*) "Illegal operation: vectors V1 and V2 should have the same size"
       STOP
     end if
     call init_to_opzero(F_sum_nd = F_sum_nd)
     do i = 1, ndim1
       call get_F1_times_F2_to_F_nd(F1_sum_nd = v1%vec_sum(i),&
                                   &F2_sum_nd = v2%vec_sum(i),&
                                   &Fres_sum_nd = sum_tmp2)
       call copy_F1_into_F2(F_sum_nd, sum_tmp1)
       call get_F1_plus_F2_to_F_sum_nd(F1_sum_nd = sum_tmp1,&
                                      & F2_sum_nd = sum_tmp2,&
                                      & Fres_sum_nd = F_sum_nd)
     end do
     call remove_opzero_in_F_sum_nd(F_sum_nd, 'Fres: from '//routine_name)
     call delete_op(sum_tmp1)
     call delete_op(sum_tmp2)
   END SUBROUTINE V1_scalar_V2_in_F_sum_nd

   !! @description: Calculates the cross product of V1 and V2
   !!               Only vectors of size = 3 are taken into account.
   !! @param:       V1    First vector (type: vec_sum_opnd)
   !! @param:       V2    Second vector (type: vec_sum_opnd)
   !! @param:       Vres  The result (type: vec_sum_opnd)
   SUBROUTINE V1_cross_V2_in_Vres(V1, V2, Vres)
     type(vec_sum_opnd),     intent(in)        :: V1
     type(vec_sum_opnd),     intent(in)        :: V2
     type(vec_sum_opnd),     intent(inout)     :: Vres

     type(sum_opnd)                  :: sum_tmp1
     type(sum_opnd)                  :: sum_tmp2
     integer                         :: i
     integer                         :: ndim1, ndim2
     logical                         :: minus

     !logical, parameter           :: debug=.TRUE.
     logical, parameter           :: debug=.FALSE.
     character (len=*), parameter :: routine_name='V1_cross_V2_in_Vres'

     IF (debug) THEN
       write(out_unitp,*) ' BEGINNING ',routine_name
       flush(out_unitp)
     END IF

     if(.not.allocated(V1%vec_sum) .or. .not.allocated(V2%vec_sum)) then
       write(out_unitp,*) ' ERROR in ',routine_name
       write(out_unitp,*) "vectors V1 and V2 should be allocated"
       STOP
     end if
     ndim1 = size(V1%vec_sum)
     ndim2 = size(V2%vec_sum)
     if(ndim1 /= ndim2  .and. ndim1 /= 3) then
       write(out_unitp,*) ' ERROR in ',routine_name
       write(out_unitp,*) "Illegal operation: vectors V1 and V2 should have the same size=3"
       STOP
     end if
     minus = .true.
     call allocate_op(Vres, ndim1)
     do i = 1, ndim1
       if(i ==1 ) then
         call get_F1_times_F2_to_F_nd(F1_sum_nd = v1%vec_sum(i+1),&
         & F2_sum_nd = v2%vec_sum(i+2),&
         & Fres_sum_nd = sum_tmp1)
         call get_F1_times_F2_to_F_nd(F1_sum_nd = v1%vec_sum(i+2),&
         & F2_sum_nd = v2%vec_sum(i+1),&
         & Fres_sum_nd = sum_tmp2)

       else if (i == 2) then
         call get_F1_times_F2_to_F_nd(F1_sum_nd = v1%vec_sum(i+1),&
         & F2_sum_nd = v2%vec_sum(i-1),&
         & Fres_sum_nd = sum_tmp1)
         call get_F1_times_F2_to_F_nd(F1_sum_nd = v1%vec_sum(i-1),&
         & F2_sum_nd = v2%vec_sum(i+1),&
         & Fres_sum_nd = sum_tmp2)

       else
         call get_F1_times_F2_to_F_nd(F1_sum_nd = v1%vec_sum(i-2),&
         & F2_sum_nd = v2%vec_sum(i-1),&
         & Fres_sum_nd = sum_tmp1)
         call get_F1_times_F2_to_F_nd(F1_sum_nd = v1%vec_sum(i-1),&
         & F2_sum_nd = v2%vec_sum(i-2),&
         & Fres_sum_nd = sum_tmp2)
       end if

       call get_F1_plus_F2_to_F_sum_nd(F1_sum_nd = sum_tmp1,&
       & F2_sum_nd = sum_tmp2,&
       & Fres_sum_nd = Vres%vec_sum(i),&
       & minus = minus)
       call remove_opzero_in_F_sum_nd(Vres%vec_sum(i), 'Fres: from '//routine_name)
     end do
     call delete_op(sum_tmp1)
     call delete_op(sum_tmp2)

     IF (debug) THEN
       write(out_unitp,*) ' END ',routine_name
       flush(out_unitp)
     END IF

   END SUBROUTINE V1_cross_V2_in_Vres


   !! @description: Calculates the sum of V1 and V2
   !! @param:       V1    First vector (type: vec_sum_opnd)
   !! @param:       V2    Second vector (type: vec_sum_opnd)
   !! @param:       Vres  The result (type: vec_sum_opnd)
   !! @param.in:    minus        Optional, if present, subtration operation
   !!                            if not present, addition operation
   SUBROUTINE V1_plus_V2_in_Vres(V1, V2, Vres, minus)
     type(vec_sum_opnd),     intent(in)        :: V1
     type(vec_sum_opnd),     intent(in)        :: V2
     type(vec_sum_opnd),     intent(inout)     :: Vres
     logical, optional,      intent(in)        :: minus

     integer                         :: i
     integer                         :: ndim1, ndim2
     logical                         :: minus_loc

     character (len=*), parameter :: routine_name='V1_plus_V2_in_Vres'

     if(.not.allocated(V1%vec_sum) .or. .not.allocated(V2%vec_sum)) then
       write(out_unitp,*) ' ERROR in ',routine_name
       write(out_unitp,*) "vectors V1 and V2 should be allocated"
       STOP
     end if
     ndim1 = size(V1%vec_sum)
     ndim2 = size(V2%vec_sum)
     if(ndim1 /= ndim2) then
       write(out_unitp,*) ' ERROR in ',routine_name
       write(out_unitp,*) "Illegal operation: vectors V1 and V2 should have the same size"
       STOP
     end if

     IF (present(minus)) THEN
       minus_loc = minus
     ELSE
       minus_loc = .FALSE.
     END IF

     call allocate_op(Vres, ndim1)
     do i = 1, ndim1
       if (minus_loc) then
         call get_F1_plus_F2_to_F_sum_nd(F1_sum_nd = V1%vec_sum(i),&
         & F2_sum_nd = V2%vec_sum(i),&
         & Fres_sum_nd = Vres%vec_sum(i),&
         & minus = .TRUE.)
       else
         call get_F1_plus_F2_to_F_sum_nd(F1_sum_nd = V1%vec_sum(i),&
         & F2_sum_nd = V2%vec_sum(i),&
         & Fres_sum_nd = Vres%vec_sum(i))
       end if
       call remove_opzero_in_F_sum_nd(Vres%vec_sum(i), 'Fres: from '//routine_name)
     end do
   END SUBROUTINE V1_plus_V2_in_Vres
   SUBROUTINE V1_PLUS_TO_Vres(V1,Vres)
     type(vec_sum_opnd),     intent(in)        :: V1
     type(vec_sum_opnd),     intent(inout)     :: Vres

     type(vec_sum_opnd)     :: Vres_loc


     character (len=*), parameter :: routine_name='V1_PLUS_TO_Vres'

     CALL copy_F1_into_F2(Vres, Vres_loc)

     CALL V1_plus_V2_in_Vres(Vres_loc, V1, Vres, minus=.FALSE.)

     CALL delete_op(Vres_loc)

   END SUBROUTINE V1_PLUS_TO_Vres
   SUBROUTINE V1_MINUS_TO_Vres(V1,Vres)
     type(vec_sum_opnd),     intent(in)        :: V1
     type(vec_sum_opnd),     intent(inout)     :: Vres

     type(vec_sum_opnd)     :: Vres_loc


     character (len=*), parameter :: routine_name='V1_MINUS_TO_Vres'

     CALL copy_F1_into_F2(Vres, Vres_loc)

     CALL V1_plus_V2_in_Vres(Vres_loc, V1, Vres, minus=.TRUE.)

     CALL delete_op(Vres_loc)

   END SUBROUTINE V1_MINUS_TO_Vres

   !! @description: Calculates the product of a F_sum_nd operator with
   !!               a vector V operator
   !! @param:       F_sum_nd The sum of nd operators (type: sum_opnd)
   !! @param:       V       The vector  (type: vec_sum_opnd)
   !! @param:       Vres  The result (type: vec_sum_opnd)
   SUBROUTINE F_sum_nd_times_V_in_Vres(F_sum_nd, V, Vres)
     type(sum_opnd),         intent(in)        :: F_sum_nd
     type(vec_sum_opnd),     intent(in)        :: V
     type(vec_sum_opnd),     intent(inout)     :: Vres

     integer                         :: i
     integer                         :: ndim

     character (len=*), parameter :: routine_name='F_sum_nd_times_V_in_Vres'

     if(.not.allocated(V%vec_sum)) then
       write(out_unitp,*) ' ERROR in ',routine_name
       write(out_unitp,*) "vector V should be allocated"
       STOP
     end if
     ndim = size(V%vec_sum)
     call allocate_op(Vres, ndim)
     do i = 1, ndim
       call get_F1_times_F2_to_F_nd(F1_sum_nd = F_sum_nd,&
                                   &F2_sum_nd = V%vec_sum(i), &
                                   &Fres_sum_nd = Vres%vec_sum(i))

     end do
     END SUBROUTINE F_sum_nd_times_V_in_Vres

   !! @description: Calculates the product of V1 with a F_sum_nd operator
   !! @param:       V       The vector  (type: vec_sum_opnd)
   !! @param:       F_sum_nd The sum of nd operators (type: sum_opnd)
   !! @param:       Vres  The result (type: vec_sum_opnd)
   SUBROUTINE V_times_F_sum_nd_in_Vres(V, F_sum_nd, Vres)
     type(vec_sum_opnd),     intent(in)        :: V
     type(sum_opnd),         intent(in)        :: F_sum_nd
     type(vec_sum_opnd),     intent(inout)     :: Vres

     integer                         :: i
     integer                         :: ndim

     character (len=*), parameter :: routine_name='V_times_F_sum_nd_in_Vres'

     if(.not.allocated(V%vec_sum)) then
       write(out_unitp,*) ' ERROR in ',routine_name
       write(out_unitp,*) "vector V should be allocated"
       STOP
     end if
     ndim = size(V%vec_sum)
     call allocate_op(Vres, ndim)
     do i = 1, ndim
       call get_F1_times_F2_to_F_nd(F1_sum_nd = V%vec_sum(i),&
                                   &F2_sum_nd = F_sum_nd, &
                                   &Fres_sum_nd = Vres%vec_sum(i))

     end do
   END SUBROUTINE V_times_F_sum_nd_in_Vres

   !! @description: Calculates the product of a matrix of type
   !!               sum_opnd with the vector of type vec_sum_opnd
   !! @param:       M_opnd    The matrix (type: sum_opnd)
   !! @param:       V         The vector (type: vec_sum_opnd)
   !! @param:       Vres      The result (type: vec_sum_opnd)
   SUBROUTINE M_opnd_times_V_in_Vres(M_opnd, V, Vres)
     type(sum_opnd),          intent(in)         :: M_opnd(:,:)
     type(vec_sum_opnd),      intent(in)         :: V
     type(vec_sum_opnd),      intent(inout)      :: Vres

     type(vec_sum_opnd)              :: Vtmp
     integer                         :: i, j
     integer                         :: ndim1, ndim2

     character (len=*), parameter :: routine_name='M_opnd_times_V_in_Vres'


     if(.not.allocated(V%vec_sum)) then
       write(out_unitp,*) ' ERROR in ',routine_name
       write(out_unitp,*) "vectors V1 and V2 should be allocated"
       STOP
     end if
     ndim1 = size(V%vec_sum)
     do i = 1, size(M_opnd(:,1))
       ndim2 = size(M_opnd(i,:))
       if(ndim1 /= ndim2) then
         write(out_unitp,*) ' ERROR in ',routine_name
         write(out_unitp,*) "Illegal operation: "
         write(out_unitp,*) "column number of the matrix should be equal "
         write(out_unitp,*) "to the size of V"
         STOP
       end if
     end do
     call allocate_op(Vres, size(M_opnd(:,1)))
     call allocate_op(Vtmp, ndim1)
     do i = 1, size(M_opnd(:,1))
       do j = 1, size(M_opnd(i,:))
         call copy_F1_into_F2(M_opnd(i,j), Vtmp%vec_sum(j))
       end do
       call V1_scalar_V2_in_F_sum_nd(Vtmp, V, Vres%vec_sum(i))
     end do
     call delete_op(Vtmp)
   END SUBROUTINE M_opnd_times_V_in_Vres

   !! @description: Calculates the product of a vector of type
   !!               vec_sum_opnd with a matrix of type sum_opnd
   !! @param:       M_opnd    The matrix (type: sum_opnd)
   !! @param:       V         The vector (type: vec_sum_opnd)
   !! @param:       Vres      The result (type: vec_sum_opnd)
   SUBROUTINE V_times_M_opnd_in_Vres(V, M_opnd, Vres)
     type(vec_sum_opnd),      intent(in)         :: V
     type(sum_opnd),          intent(in)         :: M_opnd(:,:)
     type(vec_sum_opnd),      intent(inout)      :: Vres

     type(vec_sum_opnd)              :: Vtmp
     integer                         :: i, j
     integer                         :: ndim1, ndim2

     character (len=*), parameter :: routine_name='V_times_M_opnd_in_Vres'

     if(.not.allocated(V%vec_sum)) then
       write(out_unitp,*) ' ERROR in ',routine_name
       write(out_unitp,*) "vector V should be allocated"
       STOP
     end if
     ndim1 = size(V%vec_sum)
     do i = 1, size(M_opnd(1,:))
       ndim2 = size(M_opnd(:,i))
       if(ndim1 /= ndim2) then
         write(out_unitp,*) ' ERROR in ',routine_name
         write(out_unitp,*) "Illegal operation: "
         write(out_unitp,*) "line number of the matrix should be equal "
         write(out_unitp,*) "to the size of V"
         STOP
       end if
     end do
     call allocate_op(Vres, size(M_opnd(1,:)))
     call allocate_op(Vtmp, ndim1)
     do i = 1, size(M_opnd(1,:))
       do j = 1, size(M_opnd(:,i))
         call copy_F1_into_F2(M_opnd(j,i), Vtmp%vec_sum(j))
       end do
       call V1_scalar_V2_in_F_sum_nd(V, Vtmp, Vres%vec_sum(i))
     end do
     call delete_op(Vtmp)
   END SUBROUTINE V_times_M_opnd_in_Vres


   !! @description: Calculates the product of two matrices
   !! @param:       M1    The first matrix (type: sum_opnd)
   !! @param:       M2    The second matrix (type: sum_opnd)
   !! @param:       M_res The matrix result (type: sum_opnd)
   SUBROUTINE M1_times_M2_in_Mres(M1, M2, Mres)
     type(sum_opnd),              intent(in)         :: M1(:,:)
     type(sum_opnd),              intent(in)         :: M2(:,:)
     type(sum_opnd),              intent(inout)      :: Mres(:,:)

     type(vec_sum_opnd)              :: V1_tmp
     type(vec_sum_opnd)              :: V2_tmp
     integer                         :: i, j, k
     integer                         :: ndim1, ndim2

     character (len=*), parameter :: routine_name='M1_times_M_2_in_Mres'

     ndim1 = size(M1(1,:))
     ndim2 = size(M2(:,1))
     if(ndim1 /= ndim2) then
       write(out_unitp,*) ' ERROR in ',routine_name
       write(out_unitp,*) "Illegal operation: "
       write(out_unitp,*) "the number of rows of M1 should be equal "
       write(out_unitp,*) "to the number of columns of M2 "
       STOP
     end if
     ndim1 = size(M1(:,1))
     ndim2 = size(M2(1,:))
     call allocate_op(V1_tmp, size(M2(:,1)))
     do j = 1, size(M2(1,:))
       do i = 1, size(M2(:,1))
         call copy_F1_into_F2(M2(i,j), V1_tmp%vec_sum(i))
       end do
       call M_opnd_times_V_in_Vres(M1, V1_tmp, V2_tmp)
       do k = 1, size(M1(:,1))
         call copy_F1_into_F2(V2_tmp%vec_sum(k), Mres(k,j))
       end do
       call delete_op(V2_tmp)
     end do
     call delete_op(V1_tmp)
   END SUBROUTINE M1_times_M2_in_Mres

   !! @description: Calculates the product of two matrices
   !! @param:       M1    The first matrix (type: sum_opnd)
   !! @param:       M2    The second matrix (type: sum_opnd)
   !! @param:       M_res The matrix result (type: sum_opnd)
   SUBROUTINE M1_times_M2_in_Mresbis(M1, M2, Mres)
     type(sum_opnd),              intent(in)         :: M1(:,:)
     type(sum_opnd),              intent(in)         :: M2(:,:)
     type(sum_opnd),              intent(inout)      :: Mres(:,:)

     type(vec_sum_opnd)              :: V1_tmp
     type(vec_sum_opnd)              :: V2_tmp
     integer                         :: i, j, k
     integer                         :: ndim1, ndim2

     character (len=*), parameter :: routine_name='M1_times_M_2_in_Mres'

     ndim1 = size(M1(1,:))
     ndim2 = size(M2(:,1))
     if(ndim1 /= ndim2) then
       write(out_unitp,*) ' ERROR in ',routine_name
       write(out_unitp,*) "Illegal operation: "
       write(out_unitp,*) "the number of columns of M1 should be equal "
       write(out_unitp,*) "to the number of rows of M2 "
       STOP
     end if
     ndim1 = size(M1(:,1))
     ndim2 = size(M2(1,:))
     call allocate_op(V1_tmp, size(M2(:,1)))
     do i = 1, size(M2(1,:))
       do j = 1, size(M2(:,1))
         call copy_F1_into_F2(M2(j,i), V1_tmp%vec_sum(j))
       end do
       call M_opnd_times_V_in_Vres(M1, V1_tmp, V2_tmp)
       do k = 1, size(M1(:,1))
         call copy_F1_into_F2(V2_tmp%vec_sum(k), Mres(k,i))
       end do
       call delete_op(V2_tmp)
     end do
     call delete_op(V1_tmp)
   END SUBROUTINE M1_times_M2_in_Mresbis

 end module mod_Tana_VecSumOpnD
