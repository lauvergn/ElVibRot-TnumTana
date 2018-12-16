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
module mod_Tana_Sum_OpnD
   USE mod_system
   USE mod_Tana_OpEl
   USE mod_Tana_Op1D
   USE mod_Tana_OpnD
   IMPLICIT NONE

   PRIVATE
      !-----------------------------------------------------------!
      !                        SUM_OPND                           !
      !-----------------------------------------------------------!
      !! @description: Definition of a type of a sum of nd-operators
      !!   like $ \hat{P}_{q_{k_1}} q_{k_2}^\alpha \cos^\alpha q_{k_2} +
      !!     \sin^\alpha q_{k_1}\hat{P}_{q_{k_3}} $.
      !!               This type is used for the analytical
      !!               computation of the KEO
      !! @param: sum_prod_op1d    Array of 1d operators
      !! @param: Cn               Array of complex numbers
      TYPE sum_opnd
        type(opnd), allocatable            :: sum_prod_op1d(:)
        complex(kind = Rkind), allocatable :: Cn(:)
      END TYPE sum_opnd

      INTERFACE alloc_array
        MODULE PROCEDURE alloc_array_OF_Sum_OpnDdim1,alloc_array_OF_Sum_OpnDdim2
      END INTERFACE
      INTERFACE dealloc_array
        MODULE PROCEDURE dealloc_array_OF_Sum_OpnDdim1,dealloc_array_OF_Sum_OpnDdim2
      END INTERFACE

      INTERFACE alloc_NParray
        MODULE PROCEDURE alloc_NParray_OF_Sum_OpnDdim1,alloc_NParray_OF_Sum_OpnDdim2
      END INTERFACE
      INTERFACE dealloc_NParray
        MODULE PROCEDURE dealloc_NParray_OF_Sum_OpnDdim1,dealloc_NParray_OF_Sum_OpnDdim2
      END INTERFACE

  !!@description: Generic routine that check the allocated array of operator types
  interface check_allocate_op
    module procedure check_allocate_sum_opnd
  end interface

  !!@description: Generic routine that deletes an array of operator types
  interface delete_op
    module procedure delete_sum_opnd
  end interface

  !!@description: Generic routine that initializes a variable of operator type to zero
  interface init_to_opzero
    module procedure init_opzero_sum_opnd
  end interface

  !!@description: Generic routine that allocate variables of operator type
  interface allocate_op
    module procedure allocate_sum_opnd
  end interface

   INTERFACE write_op
     module procedure  write_sum_opnd
   END INTERFACE

  !!@description: Generic routine that copy a operator F1 to another operator F2
  interface copy_F1_into_F2
    module procedure  copy_F1_sum_nd_into_F2_sum_nd, &
                      copy_F1_el_into_F2_sum_nd, copy_F1_1d_into_F2_sum_nd, &
                      copy_F1_nd_into_F2_sum_nd
  end interface

   interface get_F1_plus_F2_to_F_sum_nd
     module PROCEDURE get_F1_nd_plus_F2_nd_to_Fres_sum_nd, get_F1_nd_plus_F2_sum_nd_to_Fres_sum_nd, &
                      get_F1_sum_nd_plus_F2_nd_to_Fres_sum_nd, get_F1_sum_nd_plus_F2_sum_nd_to_Fres_sum_nd
   end interface

   INTERFACE get_F1_times_F2_to_F_nd
     module procedure get_F1_sum_nd_times_F2_sum_nd_to_Fres_sum_nd, &
                      get_F1_nd_times_F2_sum_nd_to_Fres_sum_nd, &
                      get_F1_sum_nd_times_F2_nd_to_Fres_sum_nd
   END INTERFACE

   INTERFACE operator (*)
      MODULE PROCEDURE F1_nd_times_F2_sum_nd_to_Fres_sum_nd
      MODULE PROCEDURE F1_sum_nd_times_F2_nd_to_Fres_sum_nd
      MODULE PROCEDURE F1_sum_nd_times_F2_sum_nd_to_Fres_sum_nd
   END INTERFACE


   INTERFACE assignment (=)
     MODULE PROCEDURE SumOpnD2_TO_SumOpnD1,OpnD2_TO_SumOpnD1,Op1D2_TO_SumOpnD1,OpEl2_TO_SumOpnD1
     MODULE PROCEDURE R_TO_SumOpnD1,C_TO_SumOpnD1
   END INTERFACE

   PUBLIC :: sum_opnd, allocate_op, delete_op, check_allocate_op, write_op, init_to_opzero
   PUBLIC :: Simplify_Sum_OpnD, Transpose_Mat_OF_sum_opnd

   PUBLIC :: alloc_array, dealloc_array, alloc_NParray, dealloc_NParray
   PUBLIC :: copy_F1_into_F2, get_F1_plus_F2_to_F_sum_nd, get_F1_times_F2_to_F_nd, operator (*), assignment (=)
   PUBLIC :: Der1_OF_OpnD_TO_Sum_OpnD, Der1_OF_Sum_OpnD_TO_Sum_OpnD
   PUBLIC :: Expand_Sum_OpnD_TO_Sum_OpnD, F1_sum_nd_PLUS_TO_Fres_sum_nd
   PUBLIC :: F1_nd_MINUS_TO_Fres_sum_nd, F1_sum_nd_MINUS_TO_Fres_sum_nd
   PUBLIC :: C_TO_Mat_OF_sum_opnd, remove_opzero_in_F_sum_nd


   CONTAINS 

      SUBROUTINE alloc_NParray_OF_Sum_OpnDdim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      type(sum_opnd), allocatable, intent(inout) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=1
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_NParray_OF_Sum_OpnDdim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       CALL dealloc_NParray_OF_Sum_OpnDdim1(tab,name_var,name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)

         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1)),stat=err_mem)
       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Sum_OpnD')

      END SUBROUTINE alloc_NParray_OF_Sum_OpnDdim1
      SUBROUTINE dealloc_NParray_OF_Sum_OpnDdim1(tab,name_var,name_sub)
      IMPLICIT NONE

      type(sum_opnd), allocatable, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

      integer :: i
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_NParray_OF_Sum_OpnDdim1'
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
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Sum_OpnD')

      END SUBROUTINE dealloc_NParray_OF_Sum_OpnDdim1

      SUBROUTINE alloc_array_OF_Sum_OpnDdim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      type(sum_opnd), pointer, intent(inout) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=1
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_Sum_OpnDdim1'
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
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Sum_OpnD')

      END SUBROUTINE alloc_array_OF_Sum_OpnDdim1
      SUBROUTINE dealloc_array_OF_Sum_OpnDdim1(tab,name_var,name_sub)
      IMPLICIT NONE

      type(sum_opnd), pointer, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

      integer :: i
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_Sum_OpnDdim1'
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
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Sum_OpnD')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_Sum_OpnDdim1
      SUBROUTINE alloc_array_OF_Sum_OpnDdim2(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      type(sum_opnd), pointer, intent(inout) :: tab(:,:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub


      integer, parameter :: ndim=2
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_Sum_OpnDdim2'
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
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'sum_opnd')

      END SUBROUTINE alloc_array_OF_Sum_OpnDdim2
      SUBROUTINE dealloc_array_OF_Sum_OpnDdim2(tab,name_var,name_sub)
      IMPLICIT NONE

      type(sum_opnd), pointer, intent(inout) :: tab(:,:)
      character (len=*), intent(in) :: name_var,name_sub

      integer :: i,j
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_Sum_OpnDdim2'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN

       IF (.NOT. associated(tab))                                       &
                 CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       DO i=lbound(tab,dim=2),ubound(tab,dim=2)
       DO j=lbound(tab,dim=1),ubound(tab,dim=1)
         CALL delete_op(tab(j,i))
       END DO
       END DO

       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'sum_opnd')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_Sum_OpnDdim2
      SUBROUTINE alloc_NParray_OF_Sum_OpnDdim2(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      type(sum_opnd), allocatable, intent(inout) :: tab(:,:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub


      integer, parameter :: ndim=2
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_NParray_OF_Sum_OpnDdim2'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       CALL dealloc_NParray_OF_Sum_OpnDdim2(tab,name_var,name_sub)

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
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'sum_opnd')

      END SUBROUTINE alloc_NParray_OF_Sum_OpnDdim2
      SUBROUTINE dealloc_NParray_OF_Sum_OpnDdim2(tab,name_var,name_sub)
      IMPLICIT NONE

      type(sum_opnd), allocatable, intent(inout) :: tab(:,:)
      character (len=*), intent(in) :: name_var,name_sub

      integer :: i,j
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_NParray_OF_Sum_OpnDdim2'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       IF (.NOT. allocated(tab)) RETURN

       DO i=lbound(tab,dim=2),ubound(tab,dim=2)
       DO j=lbound(tab,dim=1),ubound(tab,dim=1)
         CALL delete_op(tab(j,i))
       END DO
       END DO

       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'sum_opnd')

      END SUBROUTINE dealloc_NParray_OF_Sum_OpnDdim2

 subroutine check_allocate_sum_opnd(F_sum_nd)
   type(sum_opnd),           intent(in)    :: F_sum_nd


   character (len=*), parameter :: routine_name="check_allocate_sum_opnd"

   CALL check_NParray(F_sum_nd%sum_prod_op1d,'F_sum_nd%sum_prod_op1d',routine_name)

 end subroutine check_allocate_sum_opnd


 !! @description: Allocated a sum_opnd operator
 !! @param:    F_sum_nd    The nd operator (type: sum_opnd).
 !! @param:       ndim     Size of F_sum_nd%sum_prod_op1d.
 subroutine allocate_sum_opnd(F_sum_nd, ndim)

   type(sum_opnd),           intent(inout)    :: F_sum_nd
   integer,                  intent(in)       :: ndim

   integer                    :: i, error
   character (len=*), parameter :: routine_name="allocate_sum_opnd"

   CALL delete_sum_opnd(F_sum_nd)

   CALL alloc_NParray(F_sum_nd%sum_prod_op1d,(/ndim/),'F_sum_nd%sum_prod_op1d',routine_name)

   CALL alloc_NParray(F_sum_nd%Cn,(/ndim/),'F_sum_nd%Cn',routine_name)
   F_sum_nd%Cn = CONE

 end subroutine allocate_sum_opnd

 !! @description: Deallocated a nd operator
 !! @param:    F_sum_nd    The nd operator (type: sum_opnd).
 subroutine delete_sum_opnd(F_sum_nd)

   type(sum_opnd),           intent(inout)    :: F_sum_nd

   integer                    :: i
   character (len=*), parameter :: routine_name="delete_sum_opnd"

   if(allocated(F_sum_nd%sum_prod_op1d)) then
     CALL dealloc_NParray(F_sum_nd%sum_prod_op1d,'F_sum_nd%sum_prod_op1d',routine_name)
   end if
   if(allocated(F_sum_nd%Cn)) then
     CALL dealloc_NParray(F_sum_nd%Cn,'F_sum_nd%Cn',routine_name)
   end if

 end subroutine delete_sum_opnd

 !! @description: Initialized a sum_opnd operator to zero
 !! @param:       F_sum_nd    The nd operator (type: sum_opnd).
 subroutine init_opzero_sum_opnd(F_sum_nd)
   type(sum_opnd),           intent(inout)    :: F_sum_nd

   character (len=*), parameter :: routine_name="init_opzero_sum_nd"

   call allocate_sum_opnd(F_sum_nd,1)
   call init_to_opzero(F_sum_nd%sum_prod_op1d(1))

 end subroutine init_opzero_sum_opnd


   !! @description: Write an array of sum of nd operators,
   !! @param:       F_sum_nd      The operator (type: sum_opnd).
   !! @param:       header    If present write comment in the biginning of the
   !!                         of the file
   SUBROUTINE write_sum_opnd(F_sum_nd, i_file, header, append, close_file)
     type(sum_opnd),            intent(in)       :: F_sum_nd
     integer, optional,         intent(in)       :: i_file
     logical, optional,         intent(in)       :: header
     logical, optional,         intent(in)       :: append
     logical, optional,         intent(in)       :: close_file

     integer                        :: pq(2),J(2),L(2),nb_PJ
     character(len=1)               :: type_coef
     integer                        :: error
     integer                        :: i
     integer                        :: i_open
     logical                        :: header_loc
     character (len=*), parameter   :: routine_name='write_sum_opnd'

     if(present(i_file)) then
       i_open = i_file
     else
       i_open = out_unitp
     end if

     if(present(header)) then
       header_loc = header
     else
       header_loc = .FALSE.
     end if

     if (header_loc) then
       write(i_open, '(A)')  "   idf          iqd        alfa       &
              &indexq           op_name                  coeff     act"

     else
       write(i_open,*) '####################################&
                    &##############################################################################'
       write(i_open, *)
       write(i_open, *)
       write(i_open, '(40x, A)') '========== writing a sum of nd operators ======== '
       write(i_open, *)
     end if
     do i = 1, size(F_sum_nd%sum_prod_op1d)
       CALL get_pqJL_OF_OpnD(pq,J,L,F_sum_nd%sum_prod_op1d(i))
       nb_PJ = count(pq > 0)+count(J>0)

       IF (real(F_sum_nd%Cn(i),kind=Rkind) /= ZERO .AND. aimag(F_sum_nd%Cn(i)) /=ZERO) THEN
         type_coef = 'C'
         write(i_open, *) ' WARNING this term is complex!!'
       ELSE IF (real(F_sum_nd%Cn(i),kind=Rkind) == ZERO) THEN
         type_coef = 'I'
         IF (nb_PJ /= 1) write(i_open, *) ' WARNING this term MUST be real!!'
       ELSE
         type_coef = 'R'
         IF (nb_PJ == 1) write(i_open, *) ' WARNING this term MUST be imaginary!!'

       END IF
       write(i_open, "(A, 3x, I4, 3x, A, 1x, (E13.4,' Ix ',E13.4))") 'term', i, ', C_I=', F_sum_nd%Cn(i)
       write(i_open, *) 'term', i, ', C_I= ',type_coef, F_sum_nd%Cn(i),' nb_P+J',nb_PJ
       write(i_open, *)
       call write_op(F_sum_nd%sum_prod_op1d(i), i_open)
       write(i_open, *)
     end do
     CALL flush_perso(i_open)

   END SUBROUTINE write_sum_opnd


   !! @description: Write only the deformation part of T,
   !! @param:       F_sum_nd      The operator (type: sum_opnd).
   !! @param:       header    If present write comment in the biginning of the
   !!                         of the file
   SUBROUTINE write_sum_opnd_deform(F_sum_nd,  i_file, header, append, close_file)
     type(sum_opnd),            intent(in)       :: F_sum_nd
     integer, optional,         intent(in)       :: i_file
     logical, optional,         intent(in)       :: header
     logical, optional,         intent(in)       :: append
     logical, optional,         intent(in)       :: close_file

     integer                        :: error
     integer                        :: i, j, k, l
     integer                        :: n
     logical, pointer               :: l_J(:)

     integer                        :: i_open
     logical                        :: header_loc
     character (len=*), parameter :: routine_name='write_sum_opnd_deform'

     if(present(i_file)) then
       i_open = i_file
     else
       i_open = out_unitp
     end if

     if(present(header)) then
       header_loc = header
     else
       header_loc = .FALSE.
     end if

     nullify(l_J)
     CALL alloc_array(l_J,shape(F_sum_nd%sum_prod_op1d),'l_J',routine_name)
     l_J(:) = .true.

     if(i_open == out_unitp) then
       i_open = i_file
       do i = 1, size(F_sum_nd%sum_prod_op1d)
         do j = 1, size(F_sum_nd%sum_prod_op1d(i)%prod_op1d)
           do k = 1, size(F_sum_nd%sum_prod_op1d(i)%prod_op1d(j)%prod_opel)
             if(F_sum_nd%sum_prod_op1d(i)%prod_op1d(j)%prod_opel(k)%idf == 9 &
                &.or. F_sum_nd%sum_prod_op1d(i)%prod_op1d(j)%prod_opel(k)%idf == 10 &
                &.or.F_sum_nd%sum_prod_op1d(i)%prod_op1d(j)%prod_opel(k)%idf == 11) then
                L_J(i) = .false.
                exit
             end if
           end do
         end do
       end do
       n = 0
       do i = 1, size(F_sum_nd%sum_prod_op1d)
         if(l_J(i)) then
           n = n+1
           write(i_open, '(A, 3x, I4, 3x, A, 1x, E12.5)') 'term', n, ', C_I=', F_sum_nd%Cn(i)
           write(i_open, *)
           call write_op(F_sum_nd%sum_prod_op1d(i), i_open)
           write(i_open, *)
         end if
       end do
       IF (associated(l_J)) CALL dealloc_array(l_J,'l_J',routine_name)
       return
     end if

     if (header_loc) then
       write(i_open, '(A)')  "   idf          iqd        alfa       &
              &indexq          op_name                 qval               coeff&
              &           opval"
     else
       write(i_open,*) '####################################&
                    &##############################################################################'
       write(i_open, *)
       write(i_open, *)
       write(i_open, '(40x, A)') '========== writing a sum of nd operators ======== '
       write(i_open, *)
     end if
     do i = 1, size(F_sum_nd%sum_prod_op1d)
       do j = 1, size(F_sum_nd%sum_prod_op1d(i)%prod_op1d)
         do k = 1, size(F_sum_nd%sum_prod_op1d(i)%prod_op1d(j)%prod_opel)
           if(F_sum_nd%sum_prod_op1d(i)%prod_op1d(j)%prod_opel(k)%idf == 9 &
              &.or. F_sum_nd%sum_prod_op1d(i)%prod_op1d(j)%prod_opel(k)%idf == 10 &
              &.or.F_sum_nd%sum_prod_op1d(i)%prod_op1d(j)%prod_opel(k)%idf == 11) then
              L_J(i) = .false.
              exit
           end if
         end do
       end do
     end do
     n = 0
     do i = 1, size(F_sum_nd%sum_prod_op1d)
       if(l_J(i)) then
         n = n+1
         write(i_open, '(A, 3x, I4, 3x, A, 1x, E12.5)') 'term', n, ', C_I=', F_sum_nd%Cn(i)
         write(i_open, *)
         call write_op(F_sum_nd%sum_prod_op1d(i), i_open)
         write(i_open, *)
       end if
     end do

     IF (associated(l_J)) CALL dealloc_array(l_J,'l_J',routine_name)

   END SUBROUTINE write_sum_opnd_deform

   SUBROUTINE write_sum_Mat_OF_opnd(Mat_sum_nd, i_file, header, append, close_file)
     type(sum_opnd),            intent(in)       :: Mat_sum_nd(:,:)
     integer, optional,         intent(in)       :: i_file
     logical, optional,         intent(in)       :: header
     logical, optional,         intent(in)       :: append
     logical, optional,         intent(in)       :: close_file

     integer                        :: error
     integer                        :: i,j
     integer                        :: i_open
     logical                        :: header_loc
     character (len=*), parameter   :: routine_name='write_sum_opnd'

     if(present(i_file)) then
       i_open = i_file
     else
       i_open = out_unitp
     end if

     do i = 1, size(Mat_sum_nd(1,:))
     do j = 1, size(Mat_sum_nd(:,1))
       write(i_open,*)
       write(i_open,*) '========== writing a Mat sum of nd operators ======== '
       write(i_open,*) '========== i,j: ',i,j
       write(i_open,*)
       call write_op(Mat_sum_nd(i,j), i_open,header=.true.)
       write(i_open, *)
     end do
     end do
     CALL flush_perso(i_open)

   END SUBROUTINE write_sum_Mat_OF_opnd

 !! @description: Copy a nd operator F1_nd to
 !!               another a nd operator F2_sum_nd
 !! @param:     F1_nd        The operator which will be copied
 !! @param:    F2_sum_nd    The operator in which F1_el will be copied
 subroutine copy_F1_nd_into_F2_sum_nd(F1_nd, F2_sum_nd)

   type(opnd),           intent(in)    :: F1_nd
   type(sum_opnd),       intent(inout) :: F2_sum_nd


   character (len=*), parameter :: routine_name="copy_F1_nd_into_F2_sum_nd"

   call allocate_sum_opnd(F2_sum_nd,1)

   call copy_F1_into_F2(F1_nd, F2_sum_nd%sum_prod_op1d(1))
   F2_sum_nd%Cn(1) = cone

   CALL Simplify_Sum_OpnD(F2_sum_nd)

 end subroutine copy_F1_nd_into_F2_sum_nd

 subroutine SumOpnD2_TO_SumOpnD1(SumOpnD1,SumOpnD2)

   type(sum_opnd),       intent(in)    :: SumOpnD2
   type(sum_opnd),       intent(inout) :: SumOpnD1

   integer                    :: i

   character (len=*), parameter :: routine_name="SumOpnD2_TO_SumOpnD1"

   call allocate_op(SumOpnD1,size(SumOpnD2%sum_prod_op1d))

   do i = 1,size(SumOpnD2%sum_prod_op1d)
     SumOpnD1%sum_prod_op1d(i) = SumOpnD2%sum_prod_op1d(i)
     SumOpnD1%Cn(i)            = SumOpnD2%Cn(i)
   end do

   CALL simplify_sum_opnd(SumOpnD1)

 end subroutine SumOpnD2_TO_SumOpnD1
 subroutine OpnD2_TO_SumOpnD1(SumOpnD1,OpnD2)

   type(opnd),           intent(in)    :: OpnD2
   type(sum_opnd),       intent(inout) :: SumOpnD1

   integer                    :: i, j

   character (len=*), parameter :: routine_name="OpnD2_TO_SumOpnD1"

   call allocate_op(SumOpnD1,1)

   SumOpnD1%sum_prod_op1d(1) = OpnD2
   SumOpnD1%Cn(1) = cone

   CALL simplify_sum_opnd(SumOpnD1)

 end subroutine OpnD2_TO_SumOpnD1
 subroutine Op1D2_TO_SumOpnD1(SumOpnD1,Op1D2)

   type(op1d),           intent(in)    :: Op1D2
   type(sum_opnd),       intent(inout) :: SumOpnD1

   integer                    :: i, j

   character (len=*), parameter :: routine_name="Op1D2_TO_SumOpnD1"

   call allocate_op(SumOpnD1,1)

   SumOpnD1%sum_prod_op1d(1) = Op1D2
   SumOpnD1%Cn(1) = cone

   CALL simplify_sum_opnd(SumOpnD1)

 end subroutine Op1D2_TO_SumOpnD1
 subroutine OpEl2_TO_SumOpnD1(SumOpnD1,OpEl2)

   type(opel),           intent(in)    :: OpEl2
   type(sum_opnd),       intent(inout) :: SumOpnD1

   integer                    :: i, j

   character (len=*), parameter :: routine_name="OpEl2_TO_SumOpnD1"

   call allocate_op(SumOpnD1,1)

   SumOpnD1%sum_prod_op1d(1) = OpEl2
   SumOpnD1%Cn(1) = cone

   CALL simplify_sum_opnd(SumOpnD1)

 end subroutine OpEl2_TO_SumOpnD1
 subroutine R_TO_SumOpnD1(SumOpnD1,R)

   real (kind=Rkind),    intent(in)    :: R
   type(sum_opnd),       intent(inout) :: SumOpnD1

   character (len=*), parameter :: routine_name="R_TO_SumOpnD1"

   call allocate_op(SumOpnD1,1)

   SumOpnD1%sum_prod_op1d(1) = R
   SumOpnD1%Cn(1) = cone

   CALL simplify_sum_opnd(SumOpnD1)

 end subroutine R_TO_SumOpnD1
 subroutine C_TO_SumOpnD1(SumOpnD1,C)

   complex (kind=Rkind), intent(in)    :: C
   type(sum_opnd),       intent(inout) :: SumOpnD1

   character (len=*), parameter :: routine_name="C_TO_SumOpnD1"

   call allocate_op(SumOpnD1,1)

   SumOpnD1%sum_prod_op1d(1) = C
   SumOpnD1%Cn(1) = cone

   CALL simplify_sum_opnd(SumOpnD1)

 end subroutine C_TO_SumOpnD1
 !! @description: Copy a 1d operator F1_1d to
 !!               another a nd operator F2_sum_nd
 !! @param:     F1_1d    The operator which will be copied
 !! @param:    F2_sum_nd    The operator in which F1_el will be copied
 subroutine copy_F1_1d_into_F2_sum_nd(F1_1d, F2_sum_nd)
   type(op1d),           intent(in)    :: F1_1d
   type(sum_opnd),           intent(inout) :: F2_sum_nd

   integer                    :: i

   character (len=*), parameter :: routine_name="copy_F1_1d_into_F2_sum_nd"

   call allocate_sum_opnd(F2_sum_nd,1)

   call copy_F1_into_F2(F1_1d, F2_sum_nd%sum_prod_op1d(1))
   F2_sum_nd%Cn(1) = cone

   CALL simplify_sum_opnd(F2_sum_nd)

 end subroutine copy_F1_1d_into_F2_sum_nd


 !! @description: Copy an elementary operator F1_el to
 !!               another a nd operator F2_sum_nd
 !! @param:     F1_el    The operator which will be copied
 !! @param:    F2_sum_nd    The operator in which F1_el will be copied
 subroutine copy_F1_el_into_F2_sum_nd(F1_el, F2_sum_nd)

   type(opel),           intent(in)    :: F1_el
   type(sum_opnd),           intent(inout) :: F2_sum_nd

   character (len=*), parameter :: routine_name="copy_F1_el_into_F2_sum_nd"

   call allocate_sum_opnd(F2_sum_nd,1)

   call copy_F1_into_F2(F1_el, F2_sum_nd%sum_prod_op1d(1))
   F2_sum_nd%Cn(1) = cone

   CALL simplify_sum_opnd(F2_sum_nd)

 end subroutine copy_F1_el_into_F2_sum_nd


 !! @description: Copy a sum_opnd operator F1_sum_nd to
 !!               another sum_opnd operator F2_sum_nd
 !! @param:     F1_sum_nd    The operator which will be copied
 !! @param:    F2_sum_nd    The operator in which F1_sum_nd will be copied
 subroutine copy_F1_sum_nd_into_F2_sum_nd(F1_sum_nd, F2_sum_nd)

   type(sum_opnd),           intent(in)    :: F1_sum_nd
   type(sum_opnd),           intent(inout) :: F2_sum_nd

   integer                    :: i

   character (len=*), parameter :: routine_name="copy_F1_sum_nd_into_F2_sum_nd"

   call allocate_sum_opnd(F2_sum_nd,size(F1_sum_nd%sum_prod_op1d))

   do i = 1, size(F1_sum_nd%sum_prod_op1d)
     call copy_F1_into_F2(F1_sum_nd%sum_prod_op1d(i), &
                                F2_sum_nd%sum_prod_op1d(i))
     F2_sum_nd%Cn(i) = F1_sum_nd%Cn(i)
   end do

   CALL simplify_sum_opnd(F2_sum_nd)

 end subroutine copy_F1_sum_nd_into_F2_sum_nd


   !! @description: Does the sum of two nd operators.,
   !!               and save the result in Fres_sum_nd.
   !! @param:  F1_nd        First  operator (type: opnd). 
   !! @param:  F2_nd        second operator (type: opnd). 
   !! @param:  Fres_sum_nd  the result (type: sum_opnd).
   !! @param:  minus        Optional, if present, subtration operation
   !!                       if not present, addition operation. 
   subroutine get_F1_nd_plus_F2_nd_to_Fres_sum_nd(F1_nd, F2_nd, Fres_sum_nd, minus)
     type(opnd),       intent(in)       :: F1_nd
     type(opnd),       intent(in)       :: F2_nd
     type(sum_opnd),   intent(inout)    :: Fres_sum_nd
     logical, optional,intent(in)       :: minus

     logical                    :: minus_loc

     character (len=*), parameter :: routine_name='get_F1_1d_plus_F2_1d_to_Fres_sum_nd'

     if(present(minus)) then
       minus_loc = minus
     else
       minus_loc = .FALSE.
     end if

     call allocate_sum_opnd(Fres_sum_nd,2)

     Fres_sum_nd%sum_prod_op1d(1) = F1_nd
     Fres_sum_nd%Cn(1) = cone

     Fres_sum_nd%sum_prod_op1d(2) = F2_nd
     IF (minus_loc) THEN
       Fres_sum_nd%Cn(2) = -cone
     ELSE
       Fres_sum_nd%Cn(2) =  cone
     END IF

     CALL simplify_sum_opnd(Fres_sum_nd)

   end subroutine get_F1_nd_plus_F2_nd_to_Fres_sum_nd
   subroutine F1_nd_PLUS_TO_Fres_sum_nd(F1_nd, Fres_sum_nd)
     type(opnd),       intent(in)       :: F1_nd
     type(sum_opnd),   intent(inout)    :: Fres_sum_nd

     type(sum_opnd)                     :: Fres_sum_nd_loc

     character (len=*), parameter :: routine_name='F1_nd_PLUS_TO_Fres_sum_nd'

     Fres_sum_nd_loc = Fres_sum_nd

     CALL get_F1_plus_F2_to_F_sum_nd(F1_nd,Fres_sum_nd_loc,Fres_sum_nd)

     CALL delete_op(Fres_sum_nd_loc)

     CALL simplify_sum_opnd(Fres_sum_nd)

   end subroutine F1_nd_PLUS_TO_Fres_sum_nd
   subroutine F1_sum_nd_PLUS_TO_Fres_sum_nd(F1_sum_nd, Fres_sum_nd)
     type(sum_opnd),       intent(in)   :: F1_sum_nd
     type(sum_opnd),   intent(inout)    :: Fres_sum_nd

     type(sum_opnd)             :: Fres_sum_nd_loc

     character (len=*), parameter :: routine_name='F1_sum_nd_PLUS_TO_Fres_sum_nd'

     Fres_sum_nd_loc = Fres_sum_nd

     CALL get_F1_plus_F2_to_F_sum_nd (F1_sum_nd,Fres_sum_nd_loc,Fres_sum_nd)

     CALL delete_op(Fres_sum_nd_loc)

     CALL simplify_sum_opnd(Fres_sum_nd)

   end subroutine F1_sum_nd_PLUS_TO_Fres_sum_nd
   subroutine F1_nd_MINUS_TO_Fres_sum_nd(F1_nd, Fres_sum_nd)
     type(opnd),       intent(in)       :: F1_nd
     type(sum_opnd),   intent(inout)    :: Fres_sum_nd

     type(sum_opnd)                     :: Fres_sum_nd_loc


     character (len=*), parameter :: routine_name='F1_nd_MINUS_TO_Fres_sum_nd'

     Fres_sum_nd_loc = Fres_sum_nd

     CALL get_F1_plus_F2_to_F_sum_nd(Fres_sum_nd_loc,F1_nd,Fres_sum_nd,minus=.TRUE.)

     CALL delete_op(Fres_sum_nd_loc)

     CALL simplify_sum_opnd(Fres_sum_nd)

   end subroutine F1_nd_MINUS_TO_Fres_sum_nd
   subroutine F1_sum_nd_MINUS_TO_Fres_sum_nd(F1_sum_nd, Fres_sum_nd)
     type(sum_opnd),       intent(in)   :: F1_sum_nd
     type(sum_opnd),   intent(inout)    :: Fres_sum_nd

     type(sum_opnd)             :: Fres_sum_nd_loc

     character (len=*), parameter :: routine_name='F1_sum_nd_MINUS_TO_Fres_sum_nd'

     Fres_sum_nd_loc = Fres_sum_nd

     CALL get_F1_plus_F2_to_F_sum_nd (Fres_sum_nd_loc,F1_sum_nd,Fres_sum_nd,minus=.TRUE.)

     CALL delete_op(Fres_sum_nd_loc)

     CALL simplify_sum_opnd(Fres_sum_nd)

   end subroutine F1_sum_nd_MINUS_TO_Fres_sum_nd

   !! @description: Does the sum of an operator of type sum_opnd with 
   !!               another one of type opnd  
   !!               and save the result in Fres_sum_nd.
   !! @param:  F1_sum_nd    The first  operator (type: sum_opnd). 
   !! @param:  F2_nd        The econd  operator (type: opnd). 
   !! @param:  Fres_sum_nd  The result (type: sum_opnd).
   !! @param:  minus        Optional, if present and true, subtration operation
   !!                       if not present or false, addition operation.

   subroutine get_F1_sum_nd_plus_F2_nd_to_Fres_sum_nd(F1_sum_nd, F2_nd, Fres_sum_nd,&
                                                     & minus)

     type(sum_opnd),                intent(in)       :: F1_sum_nd
     type(opnd),                    intent(in)       :: F2_nd
     type(sum_opnd),                intent(inout)    :: Fres_sum_nd
     logical, optional,             intent(in)       :: minus
  
     integer                    :: i
     logical                    :: minus_loc

     character (len=*), parameter :: routine_name='get_F1_sum_nd_plus_F2_nd_to_Fres_sum_nd'

     if(present(minus)) then
       minus_loc = minus
     else
       minus_loc = .FALSE.
     end if

     call allocate_sum_opnd(Fres_sum_nd,size(F1_sum_nd%sum_prod_op1d)+1)

     DO i=1,size(F1_sum_nd%sum_prod_op1d)
       Fres_sum_nd%sum_prod_op1d(i) = F1_sum_nd%sum_prod_op1d(i)
       Fres_sum_nd%Cn(i) = F1_sum_nd%Cn(i)
     END DO

     Fres_sum_nd%sum_prod_op1d(size(F1_sum_nd%sum_prod_op1d)+1) = F2_nd
     IF (minus_loc) THEN
       Fres_sum_nd%Cn(size(F1_sum_nd%sum_prod_op1d)+1) = -cone
     ELSE
       Fres_sum_nd%Cn(size(F1_sum_nd%sum_prod_op1d)+1) =  cone
     END IF

     CALL simplify_sum_opnd(Fres_sum_nd)

   end subroutine get_F1_sum_nd_plus_F2_nd_to_Fres_sum_nd

   !! @description: Does the sum of an operator of type opnd with 
   !!               another one of type sum_opnd  
   !!               and save the result in Fres_sum_nd.
   !! @param:    F1_nd First  operator (type: opnd). 
   !! @param:    F2_sum_nd second  operator (type: sum_opnd). 
   !! @param:   Fres_sum_nd the result (type: sum_opnd).
   !! @param:    minus        Optional, if present and minus=t, subtraction operation
   !!                            if not present, addition operation. 
   subroutine get_F1_nd_plus_F2_sum_nd_to_Fres_sum_nd(F1_nd, F2_sum_nd, Fres_sum_nd,&
                                                     & minus)

     type(opnd),                    intent(in)       :: F1_nd
     type(sum_opnd),                intent(in)       :: F2_sum_nd
     type(sum_opnd),                intent(inout)    :: Fres_sum_nd
     logical, optional,             intent(in)       :: minus
  
     type(opel)                 :: op_zero
     integer                    :: i, j, k, m
     integer                    :: ndim1, ndim2
     integer                    :: ndim_op1d
     integer                    :: ndim_opel
     complex(kind = Rkind)         :: coeff, coeff2
     logical                    :: l_sum, l_permut
     logical                    :: minus_loc

     character (len=*), parameter :: routine_name='get_F1_nd_plus_F2_sum_nd_to_Fres_sum_nd'

     if(present(minus)) then
       minus_loc = minus
     else
       minus_loc = .FALSE.
     end if

     call allocate_sum_opnd(Fres_sum_nd,size(F2_sum_nd%sum_prod_op1d)+1)

     Fres_sum_nd%sum_prod_op1d(1) = F1_nd

     DO i=1,size(F2_sum_nd%sum_prod_op1d)
       Fres_sum_nd%sum_prod_op1d(i+1) = F2_sum_nd%sum_prod_op1d(i)
       IF (minus_loc) THEN
         Fres_sum_nd%Cn(i+1) = -F2_sum_nd%Cn(i)
       ELSE
         Fres_sum_nd%Cn(i+1) =  F2_sum_nd%Cn(i)
       END IF

     END DO

     CALL simplify_sum_opnd(Fres_sum_nd)

   end subroutine get_F1_nd_plus_F2_sum_nd_to_Fres_sum_nd


   !! @description: Does the sum of an operator of type sum_opnd with 
   !!               another one of type sum_opnd  
   !!               and save the result in Fres_sum_nd.
   !! @param:    F1_sum_nd First  operator (type: sum_opnd). 
   !! @param:    F2_sum_nd second  operator (type: sum_opnd). 
   !! @param:   Fres_sum_nd the result (type: sum_opnd).
   !! @param:    minus        Optional, if present, subtration operation
   !!                            if not present, addition operation. 
   subroutine get_F1_sum_nd_plus_F2_sum_nd_to_Fres_sum_nd(F1_sum_nd, F2_sum_nd, &
                                                         &Fres_sum_nd, minus)
     type(sum_opnd),       intent(in)       :: F1_sum_nd
     type(sum_opnd),       intent(in)       :: F2_sum_nd
     type(sum_opnd),       intent(inout)    :: Fres_sum_nd
     logical, optional,    intent(in)       :: minus
  
     integer                    :: i, n1
     logical                    :: minus_loc

     character (len=*), parameter :: routine_name='get_F1_sum_nd_plus_F2_sum_nd_to_Fres_sum_nd'

   CALL check_allocate_op(F1_sum_nd)
   CALL check_allocate_op(F2_sum_nd)

     if(present(minus)) then
       minus_loc = minus
     else
       minus_loc = .FALSE.
     end if

     n1 = size(F1_sum_nd%sum_prod_op1d)
     call allocate_sum_opnd(Fres_sum_nd,size(F1_sum_nd%sum_prod_op1d) + &
                                        size(F2_sum_nd%sum_prod_op1d))

     DO i=1,size(F1_sum_nd%sum_prod_op1d)
       Fres_sum_nd%sum_prod_op1d(i) = F1_sum_nd%sum_prod_op1d(i)
       Fres_sum_nd%Cn(i) = F1_sum_nd%Cn(i)
     END DO

     DO i=1,size(F2_sum_nd%sum_prod_op1d)
       Fres_sum_nd%sum_prod_op1d(i+n1) = F2_sum_nd%sum_prod_op1d(i)
       IF (minus_loc) THEN
         Fres_sum_nd%Cn(i+n1) = -F2_sum_nd%Cn(i)
       ELSE
         Fres_sum_nd%Cn(i+n1) =  F2_sum_nd%Cn(i)
       END IF
     END DO

     CALL simplify_sum_opnd(Fres_sum_nd)


   end subroutine get_F1_sum_nd_plus_F2_sum_nd_to_Fres_sum_nd

   !! @description: Does the product of a sum_opnd operator F1_nd with
   !!               a sum_opnd operator F2_nd
   !!               and save the result in Fres_nd.
   !! @param.in:    F1_nd the first nd operator (type: sum_opnd).
   !! @param.in:    F2_nd_the second nd operator (type: sum_opnd).
   !! @param.out:   Fres_nd  Array of nd operators in which
   !!               the result is saved (type: sum_opnd).
   SUBROUTINE get_F1_sum_nd_times_F2_sum_nd_to_Fres_sum_nd(F1_sum_nd,   &
                                                           F2_sum_nd,   &
                                                           Fres_sum_nd)
     type(sum_opnd),       intent(in)       :: F1_sum_nd
     type(sum_opnd),       intent(in)       :: F2_sum_nd
     type(sum_opnd),       intent(inout)    :: Fres_sum_nd

     integer                    :: i,j,ij,ndim1, ndim2


     !logical, parameter           :: debug=.TRUE.
     logical, parameter           :: debug=.FALSE.
     character (len=*), parameter :: routine_name='get_F1_sum_nd_times_F2_sum_nd_to_Fres_sum_nd'

     IF (debug) THEN
       write(out_unitp,*) ' BEGINNING ',routine_name
       CALL flush_perso(out_unitp)
     END IF

   CALL check_allocate_op(F1_sum_nd)
   CALL check_allocate_op(F2_sum_nd)

   ndim1 = size(F1_sum_nd%sum_prod_op1d)
   ndim2 = size(F2_sum_nd%sum_prod_op1d)
   call allocate_op(Fres_sum_nd,ndim1*ndim2)

   ij = 0
   DO i=1,ndim2
   DO j=1,ndim1
     ij = ij+1
     Fres_sum_nd%sum_prod_op1d(ij) = F1_sum_nd%sum_prod_op1d(j) * F2_sum_nd%sum_prod_op1d(i)
     Fres_sum_nd%Cn(ij)            = F1_sum_nd%Cn(j) * F2_sum_nd%Cn(i)
   END DO
   END DO

   CALL simplify_sum_opnd(Fres_sum_nd)


   IF (debug) THEN
     write(out_unitp,*) ' END ',routine_name
     CALL flush_perso(out_unitp)
   END IF

   END SUBROUTINE get_F1_sum_nd_times_F2_sum_nd_to_Fres_sum_nd

   !! @description: Does the product of a sum_opnd operator F1_sum_nd with
   !!               a opnd operator F2_nd
   !!               and save the result in Fres_nd.
   !! @param.in:    F1_sum_nd the first nd operator (type: sum_opnd).
   !! @param.in:    F2_nd_the second nd operator (type: opnd).
   !! @param.out:   Fres_nd  Array of nd operators in which
   !!               the result is saved (type: sum_opnd).
   SUBROUTINE get_F1_sum_nd_times_F2_nd_to_Fres_sum_nd(F1_sum_nd, &
                                                      & F2_nd, &
                                                      & Fres_sum_nd)
     type(sum_opnd),       intent(in)       :: F1_sum_nd
     type(opnd),           intent(in)       :: F2_nd
     type(sum_opnd),       intent(inout)    :: Fres_sum_nd

     integer                    :: i,ndim1, ndim2

     character (len=*), parameter :: routine_name='get_F1_sum_nd_times_F2_nd_to_Fres_sum_nd'

   CALL check_allocate_op(F1_sum_nd)
   CALL check_allocate_op(F2_nd)

   ndim1 = size(F1_sum_nd%sum_prod_op1d)
   ndim2 = 1
   call allocate_op(Fres_sum_nd,ndim1*ndim2)

   DO i=1,ndim1
     Fres_sum_nd%sum_prod_op1d(i) = F1_sum_nd%sum_prod_op1d(i) * F2_nd
     Fres_sum_nd%Cn(i)            = F1_sum_nd%Cn(i)
   END DO

   CALL simplify_sum_opnd(Fres_sum_nd)

   END SUBROUTINE get_F1_sum_nd_times_F2_nd_to_Fres_sum_nd

   !! @description: Does the product of a opnd operator F1_nd with
   !!               a sum_opnd operator F2_sum_nd
   !!               and save the result in Fres_nd.
   !! @param.in:    F1_nd the first nd operator (type: opnd).
   !! @param.in:    F2_sum_nd_the second nd operator (type: sum_opnd).
   !! @param.out:   Fres_nd  Array of nd operators in which
   !!               the result is saved (type: sum_opnd).
   SUBROUTINE get_F1_nd_times_F2_sum_nd_to_Fres_sum_nd(F1_nd, &
                                                                 & F2_sum_nd, &
                                                                 & Fres_sum_nd)
     type(opnd),           intent(in)       :: F1_nd
     type(sum_opnd),       intent(in)       :: F2_sum_nd
     type(sum_opnd),       intent(inout)    :: Fres_sum_nd

     integer                    :: i,ndim1, ndim2

     character (len=*), parameter :: routine_name='get_F1_nd_times_F2_sum_nd_to_Fres_sum_nd'

   CALL check_allocate_op(F2_sum_nd)
   CALL check_allocate_op(F1_nd)

   ndim1 = 1
   ndim2 = size(F2_sum_nd%sum_prod_op1d)
   call allocate_op(Fres_sum_nd,ndim1*ndim2)

   DO i=1,ndim2
     Fres_sum_nd%sum_prod_op1d(i) = F1_nd * F2_sum_nd%sum_prod_op1d(i)
     Fres_sum_nd%Cn(i)            = F2_sum_nd%Cn(i)
   END DO

   CALL simplify_sum_opnd(Fres_sum_nd)

   END SUBROUTINE get_F1_nd_times_F2_sum_nd_to_Fres_sum_nd

   FUNCTION F1_nd_times_F2_sum_nd_to_Fres_sum_nd(F1_nd,F2_sum_nd) RESULT(Fres_sum_nd)
     type(opnd),           intent(in)       :: F1_nd
     type(sum_opnd),       intent(in)       :: F2_sum_nd
     type(sum_opnd)   :: Fres_sum_nd

     character (len=*), parameter :: routine_name='F1_nd_times_F2_sum_nd_to_Fres_sum_nd'

     CALL get_F1_nd_times_F2_sum_nd_to_Fres_sum_nd(F1_nd,F2_sum_nd,Fres_sum_nd)

   END FUNCTION F1_nd_times_F2_sum_nd_to_Fres_sum_nd
   FUNCTION F1_sum_nd_times_F2_nd_to_Fres_sum_nd(F1_sum_nd,F2_nd) RESULT(Fres_sum_nd)
     type(opnd),           intent(in)       :: F2_nd
     type(sum_opnd),       intent(in)       :: F1_sum_nd
     type(sum_opnd)   :: Fres_sum_nd

     character (len=*), parameter :: routine_name='F1_sum_nd_times_F2_nd_to_Fres_sum_nd'

     CALL get_F1_sum_nd_times_F2_nd_to_Fres_sum_nd(F1_sum_nd,F2_nd,Fres_sum_nd)

   END FUNCTION F1_sum_nd_times_F2_nd_to_Fres_sum_nd
   FUNCTION F1_sum_nd_times_F2_sum_nd_to_Fres_sum_nd(F1_sum_nd,F2_sum_nd) RESULT(Fres_sum_nd)
     type(sum_opnd),       intent(in)       :: F2_sum_nd
     type(sum_opnd),       intent(in)       :: F1_sum_nd
     type(sum_opnd)   :: Fres_sum_nd

     character (len=*), parameter :: routine_name='F1_sum_nd_times_F2_sum_nd_to_Fres_sum_nd'

     CALL get_F1_sum_nd_times_F2_sum_nd_to_Fres_sum_nd(F1_sum_nd,F2_sum_nd,Fres_sum_nd)

   END FUNCTION F1_sum_nd_times_F2_sum_nd_to_Fres_sum_nd
 RECURSIVE SUBROUTINE Simplify_Sum_OpnD(SumOpnD,Expand_Sin2)

   type(sum_opnd),       intent(inout) :: SumOpnD
   logical, optional                   :: Expand_Sin2

   integer                    :: i,j
   logical                    :: Expand_Sin2_loc
   type(sum_opnd)             :: Expand_SumOpnD

   !logical, parameter           :: debug=.TRUE.
   logical, parameter           :: debug=.FALSE.
   character (len=*), parameter :: routine_name="Simplify_Sum_OpnD"

   IF (debug) THEN
       write(out_unitp,*) 'BEGINNING ',routine_name
       write(out_unitp,*) 'allocated: SumOpnD',allocated(SumOpnD%sum_prod_op1d)
       write(out_unitp,*) 'size: SumOpnD',size(SumOpnD%sum_prod_op1d)
       CALL flush_perso(out_unitp)
   END IF

   IF (present(Expand_Sin2)) THEN
     Expand_Sin2_loc = Expand_Sin2
   ELSE
     Expand_Sin2_loc = .FALSE.
   END IF

   DO i = 1, size(SumOpnD%sum_prod_op1d)
     SumOpnD%Cn(i) = SumOpnD%Cn(i) * get_coeff_OF_OpnD(SumOpnD%sum_prod_op1d(i))
     CALL Set_coeff_OF_OpnD_TO_ONE(SumOpnD%sum_prod_op1d(i))
   END DO

   CALL remove_opzero_in_F_sum_nd( SumOpnD,routine_name)

   DO i = 1,size(SumOpnD%sum_prod_op1d)
   DO j = i+1,size(SumOpnD%sum_prod_op1d)
     IF (compare_op(SumOpnD%sum_prod_op1d(i),SumOpnD%sum_prod_op1d(j))) THEN
        SumOpnD%Cn(i) = SumOpnD%Cn(i) + SumOpnD%Cn(j)
        SumOpnD%sum_prod_op1d(j) = czero
        SumOpnD%Cn(j)            = CZERO
     END IF
   END DO
   END DO

   CALL remove_opzero_in_F_sum_nd(SumOpnD,routine_name)

   IF (Expand_Sin2_loc) THEN
     CALL Expand_Sin2_IN_Sum_OpnD_TO_Sum_OpnD(SumOpnD,Expand_SumOpnD)
     SumOpnD = Expand_SumOpnD
     CALL delete_op(Expand_SumOpnD)
   END IF

   IF (debug) THEN
       write(out_unitp,*) 'allocated: SumOpnD',allocated(SumOpnD%sum_prod_op1d)
       write(out_unitp,*) 'size: SumOpnD',size(SumOpnD%sum_prod_op1d)
       write(out_unitp,*) 'END ',routine_name
       CALL flush_perso(out_unitp)
   END IF

 end subroutine Simplify_Sum_OpnD

 !! @description: Simplify a sum_op_nd operator by removind all zero operators
 !! @param:    F_sum_nd     The 1d operator (type: op1d).
 !! @param:    string      Message text. It just help to localize the problem.
 subroutine remove_opzero_in_F_sum_nd(F_sum_nd, string)

   type(sum_opnd),           intent(inout)       :: F_sum_nd
   character (len = *),      intent(in)          :: string

   type(opnd), allocatable    :: tab_sum_nd(:)
   integer                    :: i, ii
   integer                    :: j_opzero, k_opzero
   integer                    :: error, n_opzero
   logical, allocatable       :: l_zero(:)
   complex(kind=rkind), allocatable :: Cn(:)

   character (len=*), parameter :: routine_name='remove_opzero_in_F_sum_nd'


   CALL alloc_NParray(l_zero,shape(F_sum_nd%sum_prod_op1d),'l_zero',routine_name)
   l_zero(:) = .false.

   n_opzero = 0
   DO i = 1, size(F_sum_nd%sum_prod_op1d)
     CALL present_op_zero_in_F_nd(F_sum_nd%sum_prod_op1d(i), j_opzero, &
                                   & k_opzero, string, l_zero(i))
     IF (F_sum_nd%Cn(i) == czero) l_zero(i) = .TRUE.
     IF (l_zero(i)) n_opzero = n_opzero+1
   END DO

   IF (n_opzero == size(F_sum_nd%sum_prod_op1d)) THEN
     CALL allocate_op(F_sum_nd,1)
     F_sum_nd%sum_prod_op1d(1) = czero
     F_sum_nd%Cn(1) = czero

   ELSE IF (n_opzero /= 0) THEN

     ! we cannot put it in a sum_opnd temporary variable, because the copy will pass through
     ! this subroutine, and we don't want this append.
     ! therefore, we use two tables: tab_sum_nd(:) and Cn(:)
     CALL alloc_NParray(tab_sum_nd,shape(F_sum_nd%sum_prod_op1d),'tab_sum_nd',routine_name)
     CALL alloc_NParray(Cn,shape(F_sum_nd%sum_prod_op1d),'Cn',routine_name)

     DO i=1,size(tab_sum_nd)
       tab_sum_nd(i) = F_sum_nd%sum_prod_op1d(i)
       Cn(i) = F_sum_nd%Cn(i)
     END DO

     CALL allocate_op(F_sum_nd,size(tab_sum_nd)-n_opzero)

     ii=0
     DO i=1,size(tab_sum_nd)
       IF (.NOT. l_zero(i) ) THEN
         ii = ii+1
         F_sum_nd%sum_prod_op1d(ii) = tab_sum_nd(i)
         F_sum_nd%Cn(ii) = Cn(i)
       END IF
     END DO

     CALL dealloc_NParray(tab_sum_nd,'tab_sum_nd',routine_name)
     CALL dealloc_NParray(Cn,'Cn',routine_name)
   END IF

   CALL dealloc_NParray(l_zero,'l_zero',routine_name)

 end subroutine remove_opzero_in_F_sum_nd
subroutine Expand_Sin2_IN_OpnD_TO_SumOpnD(FOpnD,SumOpnD)

   type(sum_opnd),           intent(inout)    :: SumOpnD
   type(opnd),               intent(in)       :: FOpnD

   integer                    :: k,i

   type(Sum_OF_op1d)          :: temp_SumOp1D
   type(sum_opnd)             :: temp_SumOpnD



   !logical, parameter           :: debug=.TRUE.
   logical, parameter           :: debug=.FALSE.
   character (len =*), parameter :: routine_name='Expand_Sin2_IN_OpnD_TO_SumOpnD'

   IF (debug) THEN
     write(out_unitp,*) ' BEGINNING ',routine_name
     CALL write_op(FOpnD,header=.TRUE.)
     CALL flush_perso(out_unitp)
   END IF

   SumOpnD = CONE
   DO k=1,size(FOpnD%prod_op1d)
     CALL Expand_Sin2_IN_Op1D_TO_SumOp1D(FOpnD%prod_op1d(k),temp_SumOp1D)

     CALL allocate_op(temp_SumOpnD, size(temp_SumOp1D%Sum_op1D) )
     DO i=1,size(temp_SumOp1D%Sum_op1D)
       temp_SumOpnD%sum_prod_op1d(i) = temp_SumOp1D%Sum_op1D(i)
     END DO

     SumOpnD = SumOpnD * temp_SumOpnD

   END DO


  CALL delete_op(temp_SumOp1D)
  CALL delete_op(temp_SumOpnD)


  IF (debug) THEN
    write(out_unitp,*) ' Expand sin^2 in OpnD'
    CALL write_op(SumOpnD)
    write(out_unitp,*) ' END ',routine_name
    CALL flush_perso(out_unitp)
  END IF

 end subroutine Expand_Sin2_IN_OpnD_TO_SumOpnD
subroutine Expand_Sin2_IN_Sum_OpnD_TO_Sum_OpnD(F_Sum_nD,ExpandF_Sum_nD)

   type(sum_opnd),           intent(in)       :: F_Sum_nD
   type(sum_opnd),           intent(inout)    :: ExpandF_Sum_nD

   type(sum_opnd), allocatable                :: Temp_ExpandF_Sum_nD(:)

   integer                    :: i,j,ij,ndim


   !logical, parameter           :: debug=.TRUE.
   logical, parameter           :: debug=.FALSE.
   character (len =*), parameter :: routine_name='Expand_Sin2_IN_Sum_OpnD_TO_Sum_OpnD'

   IF (debug) THEN
     write(out_unitp,*) ' BEGINNING ',routine_name
     CALL write_op(F_sum_nd,header=.TRUE.)
     CALL flush_perso(out_unitp)
   END IF

   CALL alloc_NParray(Temp_ExpandF_Sum_nD,(/ size(F_sum_nd%sum_prod_op1d) /), &
                     'Temp_ExpandF_Sum_nD',routine_name)

   ndim = 0
   DO i=1,size(F_sum_nd%sum_prod_op1d)
     CALL Expand_Sin2_IN_OpnD_TO_SumOpnD(F_sum_nd%sum_prod_op1d(i),     &
                                         Temp_ExpandF_Sum_nD(i) )
     ndim = ndim + size(Temp_ExpandF_Sum_nD(i)%sum_prod_op1d)
   END DO

   CALL allocate_op(ExpandF_Sum_nD,ndim)

   ij = 0
   DO i=1,size(Temp_ExpandF_Sum_nD)
   DO j=1,size(Temp_ExpandF_Sum_nD(i)%sum_prod_op1d)
     ij = ij + 1
     ExpandF_Sum_nD%sum_prod_op1d(ij) = Temp_ExpandF_Sum_nD(i)%sum_prod_op1d(j)
     ExpandF_Sum_nD%Cn(ij) = F_sum_nd%Cn(i) * Temp_ExpandF_Sum_nD(i)%Cn(j)
   END DO
   END DO

   CALL dealloc_NParray(Temp_ExpandF_Sum_nD,                            &
                       'Temp_ExpandF_Sum_nD',routine_name)

   CALL Simplify_Sum_OpnD(ExpandF_Sum_nD,Expand_Sin2=.FALSE.)

  IF (debug) THEN
    write(out_unitp,*) ' ExpandF_sum_nd'
    CALL write_op(ExpandF_sum_nd,header=.TRUE.)
    write(out_unitp,*) ' END ',routine_name
    CALL flush_perso(out_unitp)
  END IF

 end subroutine Expand_Sin2_IN_Sum_OpnD_TO_Sum_OpnD
subroutine Expand_Sum_OpnD_TO_Sum_OpnD(F_Sum_nD,ExpandF_Sum_nD)

   type(sum_opnd),           intent(in)       :: F_Sum_nD
   type(sum_opnd),           intent(inout)    :: ExpandF_Sum_nD

   type(sum_opnd), allocatable                :: Temp_ExpandF_Sum_nD(:)

   integer                    :: i,j,ij,ndim,ndimi


   !logical, parameter           :: debug=.TRUE.
   logical, parameter           :: debug=.FALSE.
   character (len =*), parameter :: routine_name='Expand_Sum_OpnD_TO_Sum_OpnD'

   IF (debug) THEN
     write(out_unitp,*) ' BEGINNING ',routine_name
     CALL write_op(F_sum_nd,header=.TRUE.)
     CALL flush_perso(out_unitp)
   END IF

   CALL alloc_NParray(Temp_ExpandF_Sum_nD,(/ size(F_sum_nd%sum_prod_op1d) /), &
                     'Temp_ExpandF_Sum_nD',routine_name)

   ndim = 0
   DO i=1,size(F_sum_nd%sum_prod_op1d)
     !write(6,*)
     !write(6,*) 'i (sum)',i
     !CALL write_op(F_sum_nd%sum_prod_op1d(i))
     CALL Expand_OpnD_TO_SumOpnD(F_sum_nd%sum_prod_op1d(i),             &
                                 Temp_ExpandF_Sum_nD(i)%sum_prod_op1d)
     ndimi = size(Temp_ExpandF_Sum_nD(i)%sum_prod_op1d)
     CALL alloc_NParray(Temp_ExpandF_Sum_nD(i)%Cn,(/ ndimi /),'Cn',routine_name)
     Temp_ExpandF_Sum_nD(i)%Cn = CONE
     !write(6,*) 'Expansion',i
     !CALL write_op(Temp_ExpandF_Sum_nD(i))

     ndim = ndim + ndimi
   END DO

   CALL allocate_op(ExpandF_Sum_nD,ndim)

   ij = 0
   DO i=1,size(Temp_ExpandF_Sum_nD)
   DO j=1,size(Temp_ExpandF_Sum_nD(i)%sum_prod_op1d)
     ij = ij + 1
     ExpandF_Sum_nD%sum_prod_op1d(ij) = Temp_ExpandF_Sum_nD(i)%sum_prod_op1d(j)
     ExpandF_Sum_nD%Cn(ij) = F_sum_nd%Cn(i) * get_coeff_OF_OpnD(ExpandF_Sum_nD%sum_prod_op1d(ij))
     CALL Set_coeff_OF_OpnD_TO_ONE( ExpandF_Sum_nD%sum_prod_op1d(ij) )
   END DO
   END DO

   CALL dealloc_NParray(Temp_ExpandF_Sum_nD,                            &
                       'Temp_ExpandF_Sum_nD',routine_name)

   CALL Simplify_Sum_OpnD(ExpandF_Sum_nD,Expand_Sin2=.TRUE.)

  IF (debug) THEN
    write(out_unitp,*) ' ExpandF_sum_nd'
    CALL write_op(ExpandF_sum_nd,header=.TRUE.)
    write(out_unitp,*) ' END ',routine_name
    CALL flush_perso(out_unitp)
  END IF

 end subroutine Expand_Sum_OpnD_TO_Sum_OpnD

subroutine Der1_OF_OpnD_TO_Sum_OpnD(F_nD,ider,Der1_nD)

   type(opnd),           intent(in)       :: F_nD
   integer,              intent(in)       :: ider
   type(sum_opnd),       intent(inout)    :: Der1_nD

   type(Sum_OF_op1d)          :: d1Op1D ! it will contain a sum of Op1D
   type(opnd)                 :: Temp_nD,Temp2_nD

   integer                    :: i,indexQ


   !logical, parameter           :: debug=.TRUE.
   logical, parameter           :: debug=.FALSE.
   character (len =*), parameter :: routine_name='Der1_OF_OpnD_TO_Sum_OpnD'

   IF (debug) THEN
     write(out_unitp,*) ' BEGINNING ',routine_name
     write(out_unitp,*) ' ider: ',ider
     CALL write_op(F_nD,header=.TRUE.)
     CALL flush_perso(out_unitp)
     write(out_unitp,*) ' ---------------------------'
   END IF

   Temp_nD = CONE
   d1Op1D  = CZERO
   DO i=1,size(F_nD%prod_op1d)
     indexQ = get_indexQ_OF_Op1D(F_nD%prod_op1d(i))
     IF (indexQ == ider) THEN
       d1Op1D  = Der1_OF_d0Op1D(F_nD%prod_op1d(i))
     ELSE
       Temp_nD = Temp_nD * F_nD%prod_op1d(i)
     END IF
   END DO

   CALL allocate_op(Der1_nD,size(d1Op1D%Sum_Op1D))

   DO i=1,size(Der1_nD%sum_prod_Op1D)
     Der1_nD%sum_prod_Op1D(i) = Temp_nD * d1Op1D%sum_Op1D(i)
     Der1_nD%Cn(i) = CONE
   END DO

   CALL delete_op(d1Op1D)
   CALL delete_op(Temp_nD)
   CALL delete_op(Temp2_nD)

   CALL Simplify_Sum_OpnD(Der1_nD)

   IF (debug) THEN
     write(out_unitp,*) ' -------------------------'
     write(out_unitp,*) ' Der1_nD'
     CALL write_op(Der1_nD,header=.TRUE.)
     write(out_unitp,*) ' END ',routine_name
     CALL flush_perso(out_unitp)
   END IF

 end subroutine Der1_OF_OpnD_TO_Sum_OpnD
subroutine Der1_OF_Sum_OpnD_TO_Sum_OpnD(F_SumnD,ider,Der1_SumnD)

   type(sum_opnd),       intent(in)       :: F_SumnD
   integer,              intent(in)       :: ider
   type(sum_opnd),       intent(inout)    :: Der1_SumnD

   type(sum_opnd)                 :: Tmp_SumnD
   integer                        :: i


   !logical, parameter           :: debug=.TRUE.
   logical, parameter           :: debug=.FALSE.
   character (len =*), parameter :: routine_name='Der1_OF_Sum_OpnD_TO_Sum_OpnD'

   IF (debug) THEN
     write(out_unitp,*) ' BEGINNING ',routine_name
     write(out_unitp,*) ' ider: ',ider
     CALL write_op(F_SumnD,header=.TRUE.)
     CALL flush_perso(out_unitp)
   END IF

   Der1_SumnD = CZERO
   DO i=1,size(F_SumnD%sum_prod_Op1D)
     CALL Der1_OF_OpnD_TO_Sum_OpnD(F_SumnD%sum_prod_Op1D(i),ider,Tmp_SumnD)
     Tmp_SumnD%Cn(:) = Tmp_SumnD%Cn(:) * F_SumnD%Cn(i)
     !Der1_SumnD = Der1_SumnD + Tmp_SumnD
     CALL F1_sum_nd_PLUS_TO_Fres_sum_nd(Tmp_SumnD, Der1_SumnD)
   END DO

   CALL delete_op(Tmp_SumnD)

   CALL Simplify_Sum_OpnD(Der1_SumnD)

  IF (debug) THEN
    write(out_unitp,*) ' Der1_SumnD'
    CALL write_op(Der1_SumnD,header=.TRUE.)
    write(out_unitp,*) ' END ',routine_name
    CALL flush_perso(out_unitp)
  END IF

 end subroutine Der1_OF_Sum_OpnD_TO_Sum_OpnD
   !! @description: Calculates the transpose of M
   !! @param:       M    The first matrix (type: sum_opnd)
   !! @param:       Mt   Result, (type: sum_opnd)
   FUNCTION Transpose_Mat_OF_sum_opnd(M) result(Mt)
     type(sum_opnd), allocatable,         intent(inout)      :: M(:,:)
     type(sum_opnd), allocatable                             :: Mt(:,:)

     integer                         :: i, j
     integer                         :: ndim1, ndim2

     character (len=*), parameter :: routine_name='Transpose_Mat_OF_sum_opnd'


     IF (.NOT. allocated(M)) RETURN

     ndim1 = size(M(:,1))
     ndim2 = size(M(1,:))
     CALL alloc_NParray(Mt,(/ndim2,ndim1/),'M',routine_name)

     do i = 1,ndim1
     do j = 1,ndim2
        Mt(j,i) = M(i,j)
     end do
     end do

   END FUNCTION Transpose_Mat_OF_sum_opnd

   !! @description: Intiates a complex number, C, to M
   !! @param:       M    The first matrix (type: sum_opnd)
   SUBROUTINE C_TO_Mat_OF_sum_opnd(M,C)
     type(sum_opnd), allocatable,         intent(inout)      :: M(:,:)
     complex (kind=Rkind),                intent(in)         :: C

     integer                         :: i, j
     integer                         :: ndim1, ndim2

     character (len=*), parameter :: routine_name='C_TO_Mat_OF_sum_opnd'


     IF (.NOT. allocated(M)) RETURN


     do i = lbound(M,dim=1),ubound(M,dim=1)
     do j = lbound(M,dim=2),ubound(M,dim=2)
        M(i,j) = C
     end do
     end do

   END SUBROUTINE C_TO_Mat_OF_sum_opnd


end module mod_Tana_Sum_OpnD
