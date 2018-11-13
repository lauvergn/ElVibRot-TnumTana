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
MODULE mod_memory_NotPointer
      USE mod_NumParameters, only : Rkind, Ikind, ILkind
      use mod_memory, only: write_error_not_null, sub_test_tab_ub, sub_test_tab_lb, &
                            error_memo_allo, write_error_null, error_lmemo_allo,    &
                            sub_test_bigtab_ub, sub_test_bigtab_lb
      IMPLICIT NONE

      PRIVATE
      PUBLIC :: alloc_NParray,dealloc_NParray

      INTERFACE alloc_NParray
        MODULE PROCEDURE alloc_array_OF_Ldim1,alloc_array_OF_Ldim2

        MODULE PROCEDURE alloc_array_OF_I4dim1,alloc_array_OF_I8dim1
        MODULE PROCEDURE alloc_array_OF_Idim2
        MODULE PROCEDURE alloc_array_OF_Idim3,alloc_array_OF_Idim4
        MODULE PROCEDURE alloc_array_OF_Idim5

        MODULE PROCEDURE alloc_array_OF_Rdim1,alloc_array_OF_Rdim2
        MODULE PROCEDURE alloc_array_OF_Rdim3,alloc_array_OF_Rdim4
        MODULE PROCEDURE alloc_array_OF_Rdim5

        MODULE PROCEDURE alloc_array_OF_Cdim1,alloc_array_OF_Cdim2
        MODULE PROCEDURE alloc_array_OF_Cdim3,alloc_array_OF_Cdim4
        MODULE PROCEDURE alloc_array_OF_Cdim5

      END INTERFACE
      INTERFACE dealloc_NParray
        MODULE PROCEDURE dealloc_array_OF_Ldim1,dealloc_array_OF_Ldim2

        MODULE PROCEDURE dealloc_array_OF_I4dim1,dealloc_array_OF_I8dim1
        MODULE PROCEDURE dealloc_array_OF_Idim2
        MODULE PROCEDURE dealloc_array_OF_Idim3,dealloc_array_OF_Idim4
        MODULE PROCEDURE dealloc_array_OF_Idim5

        MODULE PROCEDURE dealloc_array_OF_Rdim1,dealloc_array_OF_Rdim2
        MODULE PROCEDURE dealloc_array_OF_Rdim3,dealloc_array_OF_Rdim4
        MODULE PROCEDURE dealloc_array_OF_Rdim5

        MODULE PROCEDURE dealloc_array_OF_Cdim1,dealloc_array_OF_Cdim2
        MODULE PROCEDURE dealloc_array_OF_Cdim3,dealloc_array_OF_Cdim4
        MODULE PROCEDURE dealloc_array_OF_Cdim5

      END INTERFACE

      CONTAINS

      !=================================================================
      ! Logical
      !=================================================================
      SUBROUTINE alloc_array_OF_Ldim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      logical, allocatable, intent(inout) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=1

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_Ldim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (allocated(tab))                                             &
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
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'logical')

      END SUBROUTINE alloc_array_OF_Ldim1
      SUBROUTINE dealloc_array_OF_Ldim1(tab,name_var,name_sub)
      IMPLICIT NONE

      logical, allocatable, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_Ldim1'
      integer :: err_mem
      integer (kind=ILkind) :: memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. allocated(tab)) RETURN
       IF (.NOT. allocated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab,kind=ILkind)
       deallocate(tab,stat=err_mem)
       CALL error_lmemo_allo(err_mem,-memory,name_var,name_sub,'logical')


      END SUBROUTINE dealloc_array_OF_Ldim1
      SUBROUTINE alloc_array_OF_Ldim2(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      logical, allocatable, intent(inout) :: tab(:,:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub


      integer, parameter :: ndim=2

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_Ldim2'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       IF (allocated(tab))                                             &
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
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'logical')

      END SUBROUTINE alloc_array_OF_Ldim2
      SUBROUTINE dealloc_array_OF_Ldim2(tab,name_var,name_sub)
      IMPLICIT NONE

      logical, allocatable, intent(inout) :: tab(:,:)
      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_Ldim2'
      integer :: err_mem
      integer (kind=ILkind) :: memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. allocated(tab)) RETURN

       IF (.NOT. allocated(tab))                                       &
                 CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab,kind=ILkind)
       deallocate(tab,stat=err_mem)
       CALL error_lmemo_allo(err_mem,-memory,name_var,name_sub,'logical')


      END SUBROUTINE dealloc_array_OF_Ldim2

      !=================================================================
      ! integer
      !=================================================================

      SUBROUTINE alloc_array_OF_I8dim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      integer (kind=ILkind), allocatable, intent(inout) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub


      integer, parameter :: ndim=1


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_I8dim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (allocated(tab))                                             &
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
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'integer')

      END SUBROUTINE alloc_array_OF_I8dim1
      SUBROUTINE dealloc_array_OF_I8dim1(tab,name_var,name_sub)
      IMPLICIT NONE

      integer (kind=ILkind), allocatable, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_I8dim1'
      integer :: err_mem
      integer (kind=ILkind) :: memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. allocated(tab)) RETURN

       IF (.NOT. allocated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab,kind=ILkind)
       deallocate(tab,stat=err_mem)
       CALL error_lmemo_allo(err_mem,-memory,name_var,name_sub,'integer')


      END SUBROUTINE dealloc_array_OF_I8dim1
      SUBROUTINE alloc_array_OF_I4dim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      integer (kind=Ikind), allocatable, intent(inout) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub


      integer, parameter :: ndim=1


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_I4dim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (allocated(tab))                                             &
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
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'integer')

      END SUBROUTINE alloc_array_OF_I4dim1
      SUBROUTINE dealloc_array_OF_I4dim1(tab,name_var,name_sub)
      IMPLICIT NONE

      integer (kind=Ikind), allocatable, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_I4dim1'
      integer :: err_mem
      integer (kind=ILkind) :: memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. allocated(tab)) RETURN

       IF (.NOT. allocated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab,kind=ILkind)
       deallocate(tab,stat=err_mem)
       CALL error_lmemo_allo(err_mem,-memory,name_var,name_sub,'integer')


      END SUBROUTINE dealloc_array_OF_I4dim1
      SUBROUTINE alloc_array_OF_Idim2(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      integer, allocatable, intent(inout) :: tab(:,:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub


      integer, parameter :: ndim=2


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_Idim2'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (allocated(tab))                                             &
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
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'integer')

      END SUBROUTINE alloc_array_OF_Idim2
      SUBROUTINE dealloc_array_OF_Idim2(tab,name_var,name_sub)
      IMPLICIT NONE

      integer, allocatable, intent(inout) :: tab(:,:)
      character (len=*), intent(in) :: name_var,name_sub



!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_Idim2'
      integer :: err_mem
      integer (kind=ILkind) :: memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. allocated(tab)) RETURN

       IF (.NOT. allocated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab,kind=ILkind)
       deallocate(tab,stat=err_mem)
       CALL error_lmemo_allo(err_mem,-memory,name_var,name_sub,'integer')


      END SUBROUTINE dealloc_array_OF_Idim2

      SUBROUTINE alloc_array_OF_Idim3(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      integer, allocatable, intent(inout) :: tab(:,:,:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub


      integer, parameter :: ndim=3

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_Idim3'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (allocated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)


         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1),                              &
                      tab_lb(2):tab_ub(2),                              &
                      tab_lb(3):tab_ub(3)),stat=err_mem)

       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1),tab_ub(2),tab_ub(3)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'integer')

      END SUBROUTINE alloc_array_OF_Idim3
      SUBROUTINE dealloc_array_OF_Idim3(tab,name_var,name_sub)
      IMPLICIT NONE

      integer, allocatable, intent(inout) :: tab(:,:,:)
      character (len=*), intent(in) :: name_var,name_sub



!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_Idim3'
      integer :: err_mem
      integer (kind=ILkind) :: memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. allocated(tab)) RETURN

       IF (.NOT. allocated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab,kind=ILkind)
       deallocate(tab,stat=err_mem)
       CALL error_lmemo_allo(err_mem,-memory,name_var,name_sub,'integer')


      END SUBROUTINE dealloc_array_OF_Idim3

      SUBROUTINE alloc_array_OF_Idim4(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      integer, allocatable, intent(inout) :: tab(:,:,:,:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub


      integer, parameter :: ndim=4

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_Idim4'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (allocated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)


         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1),                              &
                      tab_lb(2):tab_ub(2),                              &
                      tab_lb(3):tab_ub(3),                              &
                      tab_lb(4):tab_ub(4)),stat=err_mem)

       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1),tab_ub(2),tab_ub(3),tab_ub(4)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'integer')

      END SUBROUTINE alloc_array_OF_Idim4
      SUBROUTINE dealloc_array_OF_Idim4(tab,name_var,name_sub)
      IMPLICIT NONE

      integer, allocatable, intent(inout) :: tab(:,:,:,:)
      character (len=*), intent(in) :: name_var,name_sub



!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_Idim4'
      integer :: err_mem
      integer (kind=ILkind) :: memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. allocated(tab)) RETURN

       IF (.NOT. allocated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab,kind=ILkind)
       deallocate(tab,stat=err_mem)
       CALL error_lmemo_allo(err_mem,-memory,name_var,name_sub,'integer')


      END SUBROUTINE dealloc_array_OF_Idim4

      SUBROUTINE alloc_array_OF_Idim5(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      integer, allocatable, intent(inout) :: tab(:,:,:,:,:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub


      integer, parameter :: ndim=5

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_Idim5'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (allocated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)


         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1),                              &
                      tab_lb(2):tab_ub(2),                              &
                      tab_lb(3):tab_ub(3),                              &
                      tab_lb(4):tab_ub(4),                              &
                      tab_lb(5):tab_ub(5)),stat=err_mem)

       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1),tab_ub(2),tab_ub(3),tab_ub(4),tab_ub(5)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'integer')

      END SUBROUTINE alloc_array_OF_Idim5
      SUBROUTINE dealloc_array_OF_Idim5(tab,name_var,name_sub)
      IMPLICIT NONE

      integer, allocatable, intent(inout) :: tab(:,:,:,:,:)
      character (len=*), intent(in) :: name_var,name_sub



!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_Idim5'
      integer :: err_mem
      integer (kind=ILkind) :: memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. allocated(tab)) RETURN

       IF (.NOT. allocated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab,kind=ILkind)
       deallocate(tab,stat=err_mem)
       CALL error_lmemo_allo(err_mem,-memory,name_var,name_sub,'integer')


      END SUBROUTINE dealloc_array_OF_Idim5

      !=================================================================
      ! real (kind=Rkind)
      !=================================================================

      SUBROUTINE alloc_array_OF_Rdim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      real (kind=Rkind), allocatable, intent(inout) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub


      integer, parameter :: ndim=1

      !----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_Rdim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
      !----- for debuging --------------------------------------------------


       IF (allocated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         !write(6,*) 'tab_lb',tab_lb
         !write(6,*) 'tab_ub',tab_ub

         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)

         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1)),stat=err_mem)
       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'real8')

       !write(6,*) 'lbound',lbound(tab)
       !write(6,*) 'ubound',ubound(tab)

      END SUBROUTINE alloc_array_OF_Rdim1

      SUBROUTINE alloc_BigArray_OF_Rdim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      real (kind=Rkind), allocatable, intent(inout) :: tab(:)
      integer (kind=ILkind), intent(in) :: tab_ub(:)
      integer(kind=ILkind), intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub


      integer, parameter :: ndim=1

      !----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_BigArray_OF_Rdim1'
      integer :: err_mem
      integer (kind=ILkind) :: memory
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
      !----- for debuging --------------------------------------------------


       IF (allocated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       CALL sub_test_Bigtab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         !write(6,*) 'tab_lb',tab_lb
         !write(6,*) 'tab_ub',tab_ub

         CALL sub_test_Bigtab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)

         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1)),stat=err_mem)
       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1)),stat=err_mem)
       END IF
       CALL error_lmemo_allo(err_mem,memory,name_var,name_sub,'real8')

       !write(6,*) 'lbound',lbound(tab)
       !write(6,*) 'ubound',ubound(tab)

      END SUBROUTINE alloc_BigArray_OF_Rdim1
      SUBROUTINE dealloc_array_OF_Rdim1(tab,name_var,name_sub)
      IMPLICIT NONE

      real (kind=Rkind), allocatable, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_Rdim1'
      integer :: err_mem
      integer (kind=ILkind) :: memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. allocated(tab)) RETURN

       IF (.NOT. allocated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab,kind=ILkind)
       deallocate(tab,stat=err_mem)
       CALL error_lmemo_allo(err_mem,-memory,name_var,name_sub,'real8')


      END SUBROUTINE dealloc_array_OF_Rdim1
      SUBROUTINE alloc_array_OF_Rdim2(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      real (kind=Rkind), allocatable, intent(inout) :: tab(:,:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub


      integer, parameter :: ndim=2

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_Rdim2'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (allocated(tab))                                             &
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
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'real8')

      END SUBROUTINE alloc_array_OF_Rdim2
      SUBROUTINE dealloc_array_OF_Rdim2(tab,name_var,name_sub)
      IMPLICIT NONE

      real (kind=Rkind), allocatable, intent(inout) :: tab(:,:)
      character (len=*), intent(in) :: name_var,name_sub



!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_Rdim2'
      integer :: err_mem
      integer (kind=ILkind) :: memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. allocated(tab)) RETURN

       IF (.NOT. allocated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab,kind=ILkind)
       deallocate(tab,stat=err_mem)
       CALL error_lmemo_allo(err_mem,-memory,name_var,name_sub,'real8')


      END SUBROUTINE dealloc_array_OF_Rdim2

      SUBROUTINE alloc_array_OF_Rdim3(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      real (kind=Rkind), allocatable, intent(inout) :: tab(:,:,:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub


      integer, parameter :: ndim=3

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_Rdim3'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (allocated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)


         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1),                              &
                      tab_lb(2):tab_ub(2),                              &
                      tab_lb(3):tab_ub(3)),stat=err_mem)

       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1),tab_ub(2),tab_ub(3)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'real8')

      END SUBROUTINE alloc_array_OF_Rdim3
      SUBROUTINE dealloc_array_OF_Rdim3(tab,name_var,name_sub)
      IMPLICIT NONE

      real (kind=Rkind), allocatable, intent(inout) :: tab(:,:,:)
      character (len=*), intent(in) :: name_var,name_sub



!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_Rdim3'
      integer :: err_mem
      integer (kind=ILkind) :: memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. allocated(tab)) RETURN

       IF (.NOT. allocated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab,kind=ILkind)
       deallocate(tab,stat=err_mem)
       CALL error_lmemo_allo(err_mem,-memory,name_var,name_sub,'real8')


      END SUBROUTINE dealloc_array_OF_Rdim3

      SUBROUTINE alloc_array_OF_Rdim4(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      real (kind=Rkind), allocatable, intent(inout) :: tab(:,:,:,:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub


      integer, parameter :: ndim=4

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_Rdim4'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (allocated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)


         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1),                              &
                      tab_lb(2):tab_ub(2),                              &
                      tab_lb(3):tab_ub(3),                              &
                      tab_lb(4):tab_ub(4)),stat=err_mem)

       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1),tab_ub(2),tab_ub(3),tab_ub(4)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'real8')

      END SUBROUTINE alloc_array_OF_Rdim4
      SUBROUTINE dealloc_array_OF_Rdim4(tab,name_var,name_sub)
      IMPLICIT NONE

      real (kind=Rkind), allocatable, intent(inout) :: tab(:,:,:,:)
      character (len=*), intent(in) :: name_var,name_sub



!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_Rdim4'
      integer :: err_mem
      integer (kind=ILkind) :: memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. allocated(tab)) RETURN

       IF (.NOT. allocated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab,kind=ILkind)
       deallocate(tab,stat=err_mem)
       CALL error_lmemo_allo(err_mem,-memory,name_var,name_sub,'real8')


      END SUBROUTINE dealloc_array_OF_Rdim4

      SUBROUTINE alloc_array_OF_Rdim5(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      real (kind=Rkind), allocatable, intent(inout) :: tab(:,:,:,:,:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub


      integer, parameter :: ndim=5

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_Rdim5'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (allocated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)


         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1),                              &
                      tab_lb(2):tab_ub(2),                              &
                      tab_lb(3):tab_ub(3),                              &
                      tab_lb(4):tab_ub(4),                              &
                      tab_lb(5):tab_ub(5)),stat=err_mem)

       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1),tab_ub(2),tab_ub(3),tab_ub(4),tab_ub(5)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'real8')

      END SUBROUTINE alloc_array_OF_Rdim5
      SUBROUTINE dealloc_array_OF_Rdim5(tab,name_var,name_sub)
      IMPLICIT NONE

      real (kind=Rkind), allocatable, intent(inout) :: tab(:,:,:,:,:)
      character (len=*), intent(in) :: name_var,name_sub



!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_Rdim5'
      integer :: err_mem
      integer (kind=ILkind) :: memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. allocated(tab)) RETURN

       IF (.NOT. allocated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab,kind=ILkind)
       deallocate(tab,stat=err_mem)
       CALL error_lmemo_allo(err_mem,-memory,name_var,name_sub,'real8')


      END SUBROUTINE dealloc_array_OF_Rdim5

      !=================================================================
      ! complex (kind=Rkind)
      !=================================================================

      SUBROUTINE alloc_array_OF_Cdim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      complex (kind=Rkind), allocatable, intent(inout) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub


      integer, parameter :: ndim=1

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_Cdim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (allocated(tab))                                             &
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
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'complex8')

      END SUBROUTINE alloc_array_OF_Cdim1
      SUBROUTINE dealloc_array_OF_Cdim1(tab,name_var,name_sub)
      IMPLICIT NONE

      complex (kind=Rkind), allocatable, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_Cdim1'
      integer :: err_mem
      integer (kind=ILkind) :: memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. allocated(tab)) RETURN

       IF (.NOT. allocated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab,kind=ILkind)
       deallocate(tab,stat=err_mem)
       CALL error_lmemo_allo(err_mem,-memory,name_var,name_sub,'complex8')


      END SUBROUTINE dealloc_array_OF_Cdim1
      SUBROUTINE alloc_array_OF_Cdim2(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      complex (kind=Rkind), allocatable, intent(inout) :: tab(:,:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub


      integer, parameter :: ndim=2

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_Cdim2'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (allocated(tab))                                             &
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
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'complex8')

      END SUBROUTINE alloc_array_OF_Cdim2
      SUBROUTINE dealloc_array_OF_Cdim2(tab,name_var,name_sub)
      IMPLICIT NONE

      complex (kind=Rkind), allocatable, intent(inout) :: tab(:,:)
      character (len=*), intent(in) :: name_var,name_sub



!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_Cdim2'
      integer :: err_mem
      integer (kind=ILkind) :: memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. allocated(tab)) RETURN

       IF (.NOT. allocated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab,kind=ILkind)
       deallocate(tab,stat=err_mem)
       CALL error_lmemo_allo(err_mem,-memory,name_var,name_sub,'complex8')


      END SUBROUTINE dealloc_array_OF_Cdim2

      SUBROUTINE alloc_array_OF_Cdim3(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      complex (kind=Rkind), allocatable, intent(inout) :: tab(:,:,:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub


      integer, parameter :: ndim=3

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_Cdim3'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (allocated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)


         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1),                              &
                      tab_lb(2):tab_ub(2),                              &
                      tab_lb(3):tab_ub(3)),stat=err_mem)

       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1),tab_ub(2),tab_ub(3)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'complex8')

      END SUBROUTINE alloc_array_OF_Cdim3
      SUBROUTINE dealloc_array_OF_Cdim3(tab,name_var,name_sub)
      IMPLICIT NONE

      complex (kind=Rkind), allocatable, intent(inout) :: tab(:,:,:)
      character (len=*), intent(in) :: name_var,name_sub



!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_Cdim3'
      integer :: err_mem
      integer (kind=ILkind) :: memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. allocated(tab)) RETURN

       IF (.NOT. allocated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab,kind=ILkind)
       deallocate(tab,stat=err_mem)
       CALL error_lmemo_allo(err_mem,-memory,name_var,name_sub,'complex8')


      END SUBROUTINE dealloc_array_OF_Cdim3

      SUBROUTINE alloc_array_OF_Cdim4(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      complex (kind=Rkind), allocatable, intent(inout) :: tab(:,:,:,:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub


      integer, parameter :: ndim=4

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_Cdim4'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (allocated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)


         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1),                              &
                      tab_lb(2):tab_ub(2),                              &
                      tab_lb(3):tab_ub(3),                              &
                      tab_lb(4):tab_ub(4)),stat=err_mem)

       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1),tab_ub(2),tab_ub(3),tab_ub(4)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'complex8')

      END SUBROUTINE alloc_array_OF_Cdim4
      SUBROUTINE dealloc_array_OF_Cdim4(tab,name_var,name_sub)
      IMPLICIT NONE

      complex (kind=Rkind), allocatable, intent(inout) :: tab(:,:,:,:)
      character (len=*), intent(in) :: name_var,name_sub



!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_Cdim4'
      integer :: err_mem
      integer (kind=ILkind) :: memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. allocated(tab)) RETURN

       IF (.NOT. allocated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab,kind=ILkind)
       deallocate(tab,stat=err_mem)
       CALL error_lmemo_allo(err_mem,-memory,name_var,name_sub,'complex8')

      END SUBROUTINE dealloc_array_OF_Cdim4

      SUBROUTINE alloc_array_OF_Cdim5(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      complex (kind=Rkind), allocatable, intent(inout) :: tab(:,:,:,:,:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub


      integer, parameter :: ndim=5

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_Cdim5'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (allocated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)


         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1),                              &
                      tab_lb(2):tab_ub(2),                              &
                      tab_lb(3):tab_ub(3),                              &
                      tab_lb(4):tab_ub(4),                              &
                      tab_lb(5):tab_ub(5)),stat=err_mem)

       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1),tab_ub(2),tab_ub(3),tab_ub(4),tab_ub(5)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'complex8')

      END SUBROUTINE alloc_array_OF_Cdim5
      SUBROUTINE dealloc_array_OF_Cdim5(tab,name_var,name_sub)
      IMPLICIT NONE

      complex (kind=Rkind), allocatable, intent(inout) :: tab(:,:,:,:,:)
      character (len=*), intent(in) :: name_var,name_sub



!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_Cdim5'
      integer :: err_mem
      integer (kind=ILkind) :: memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. allocated(tab)) RETURN

       IF (.NOT. allocated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab,kind=ILkind)
       deallocate(tab,stat=err_mem)
       CALL error_lmemo_allo(err_mem,-memory,name_var,name_sub,'complex8')


      END SUBROUTINE dealloc_array_OF_Cdim5

END MODULE mod_memory_NotPointer

