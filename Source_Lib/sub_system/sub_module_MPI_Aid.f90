!=======================================================================================
!> module for MPI auxiliary functions 
!> include:
!>  - type multi_array
!>  - type multi_array4 integer4 restricted 
!>  - subroutine allocate_array   for allocate and initialize array
!>  - subroutine system_mem_usage for memmory monitoring
!>  - subroutine MPI_write        for MPI writing control
!>  -
!=======================================================================================
MODULE mod_MPI_Aid
#if(run_MPI)
  USE mod_system
  IMPLICIT NONE

  !-------------------------------------------------------------------------------------
  !> varilables for convenience, temprary here
  Integer                        :: i1_loop          !< indexs for loop
  Integer                        :: i2_loop
  Integer                        :: i3_loop
  Integer                        :: i4_loop
  Integer                        :: i5_loop
  Integer                        :: i6_loop
  Integer                        :: i1_length        !< boundary for looop
  Integer                        :: i2_length
  Integer                        :: i3_length
  Integer                        :: i4_length
  Integer                        :: i5_length
  Integer                        :: i6_length

  Complex(kind=Rkind)            :: temp_cplx        !< temparay Complex
  Complex(kind=Rkind)            :: temp_cplx1
  Complex(kind=Rkind)            :: temp_cplx2
  Real(kind=Rkind)               :: temp_real        !< temparay real
  Real(kind=Rkind)               :: temp_real1
  Real(kind=Rkind)               :: temp_real2
  Integer                        :: temp_int         !< temparay integer 
  Integer                        :: temp_int1
  Integer                        :: temp_int2
  Logical                        :: temp_logi
    
  !-------------------------------------------------------------------------------------
  INTEGER                       :: memory_RSS   !< memory usage
  
  !-------------------------------------------------------------------------------------
  TYPE multi_array
    Integer,allocatable :: array(:)
  END TYPE multi_array
  
  TYPE multi_array4
    Integer*4,allocatable :: array(:)
  END TYPE multi_array4  
  
  !> increase matrix size. the origin value is kept. 
  !! the increase way is decided by the variables presented
  INTERFACE increase_martix
    module procedure increase_martix_int
    module procedure increase_martix_real
  END INTERFACE
  
  !> allocate and initial array 
  INTERFACE allocate_array  
    module procedure allocate_array_int4_length4
    module procedure allocate_array_int4_length8 
    module procedure allocate_array_int8_length4
    module procedure allocate_array_int8_length8
    module procedure allocate_array_real_length4
    module procedure allocate_array_real_length8
    module procedure allocate_array_cplx_length4
    module procedure allocate_array_cplx_length8
    module procedure allocate_array_logic
    module procedure allocate_matrix_int4_length4
    module procedure allocate_matrix_int4_length8 
    module procedure allocate_matrix_int8_length4
    module procedure allocate_matrix_int8_length8
    module procedure allocate_matrix_real_length4
    module procedure allocate_matrix_real_length8
    module procedure allocate_matrix_cplx_length4
    module procedure allocate_matrix_cplx_length8
  END INTERFACE
  
  ! do not use for large matrix
  INTERFACE MPI_Bcast_matrix  
    module procedure MPI_Bcast_matrix_real 
    module procedure MPI_Bcast_matrix_int 
    module procedure MPI_Bcast_matrix_complex 
  END INTERFACE
  
  INTERFACE MPI_Reduce_sum_matrix  
    module procedure MPI_Reduce_sum_matrix_real 
    module procedure MPI_Reduce_sum_matrix_int 
    module procedure MPI_Reduce_sum_matrix_complex 
  END INTERFACE
  
  INTERFACE MPI_Reduce_sum_Bcast
    module procedure MPI_Reduce_sum_Bcast_real
    module procedure MPI_Reduce_sum_Bcast_int
    module procedure MPI_Reduce_sum_Bcast_complex
  END INTERFACE
  
  INTERFACE MPI_Send_matrix  
    module procedure MPI_Send_matrix_real 
    module procedure MPI_Send_matrix_int 
    module procedure MPI_Send_matrix_complex 
  END INTERFACE
  
  INTERFACE MPI_Recv_matrix  
    module procedure MPI_Recv_matrix_real 
    module procedure MPI_Recv_matrix_int 
    module procedure MPI_Recv_matrix_complex 
  END INTERFACE
  
  INTERFACE MPI_combine_array 
    module procedure MPI_combine_array_real
    module procedure MPI_combine_array_int
    module procedure MPI_combine_array_cplx
  END INTERFACE
  
  INTERFACE MPI_collect_info
    module procedure MPI_collect_info_int
    module procedure MPI_collect_info_real
  END INTERFACE

!---------------------------------------------------------------------------------------
  Contains
    !> check total memory used at certain point
!    SUBROUTINE system_mem_usage(memory_RSS,name)
!      USE mod_NumParameters
!      ! USE ifport ! if on intel compiler
!      IMPLICIT NONE
!      Integer, intent(out) :: memory_RSS
!      Character(len=200):: filename=' '
!      Character(len=80) :: line
!      Character(len=8)  :: pid_char=' '
!      Integer :: pid
!      Logical :: ifxst
!      Character (len=*), intent(in) :: name
!
!      memory_RSS=-1 ! return negative number if not found
!
!      !> get process ID
!      pid=getpid()
!      !write(*,*) 'pid=',pid
!      write(pid_char,'(I8)') pid
!      filename='/proc/'//trim(adjustl(pid_char))//'/status'
!
!      ! read system file
!      inquire (file=filename,exist=ifxst)
!      IF(.not.ifxst) THEN
!        !write (*,*) 'system file does not exist'
!      ELSE
!        OPEN(unit=100, file=filename, action='read')
!        DO
!          read(100,'(a)',end=120) line
!          IF(line(1:6).eq.'VmRSS:') THEN
!            read (line(7:),*) memory_RSS
!            EXIT
!          ENDIF
!        ENDDO
!120     CONTINUE
!        CLOSE(100)
!        write(out_unitp,121) name,memory_RSS,MPI_id
!121     format('memory check at ',a,': ',i4,' from ',i4)
!      ENDIF
!    ENDSUBROUTINE system_mem_usage

!---------------------------------------------------------------------------------------
!> write for MPI outpout
!---------------------------------------------------------------------------------------    
    SUBROUTINE MPI_write(out_channel,out_message)
      IMPLICIT NONE
      integer :: out_channel,MPIid
      Character (len=*), intent(in) :: out_message
      
      write(out_channel,*) out_message,' from ',MPI_id
    ENDSUBROUTINE MPI_write
    
    SUBROUTINE MPI0_write(out_channel,out_message)
      IMPLICIT NONE
      integer :: out_channel
      Character (len=*), intent(in) :: out_message
      
      IF(MPI_id==0) write(out_channel,*) out_message
    ENDSUBROUTINE MPI0_write
    
    SUBROUTINE MPI0_write_line(out_channel)
      IMPLICIT NONE
      integer :: out_channel,MPIid
            
      IF(MPI_id==0) write(out_channel,*) '------------------------------------------------------------'
    ENDSUBROUTINE MPI0_write_line
    
    SUBROUTINE MPI0_write_dline(out_channel)
      IMPLICIT NONE
      integer :: out_channel,MPIid
      
      IF(MPI_id==0) write(out_channel,*) '============================================================'
    ENDSUBROUTINE MPI0_write_dline
!---------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------
!> interface: increase_martix    
!---------------------------------------------------------------------------------------
    SUBROUTINE increase_martix_int(matrix,name_sub,ndim0,ndim,double_size)
      IMPLICIT NONE
      
      Integer,allocatable,intent(inout)            :: matrix(:,:)
      Integer,optional,intent(in)                  :: ndim0
      Integer,optional,intent(in)                  :: ndim
      Logical,optional,intent(in)                  :: double_size
      Character(len=*),intent(in)                  :: name_sub
      
      Integer,allocatable                          :: matrix_temp(:,:)
      Integer                                      :: new_dim

      IF(.NOT. allocated(matrix)) THEN
        CALL alloc_NParray(matrix,(/ ndim,ndim /),"increase matrix",name_sub)
      ELSE
        IF(present(ndim0)) THEN
          IF(present(ndim)) THEN
            new_dim=ndim
          ELSE IF(present(double_size)) THEN
            new_dim=ndim0*2
          ELSE
            STOP 'error in the variables for increase matrix size.'
          ENDIF
          
          CALL alloc_NParray(matrix_temp,(/ ndim,ndim /),"increase matrix",name_sub)
          matrix_temp(1:ndim0,1:ndim0)=matrix(1:ndim0,1:ndim0)
          CALL move_alloc(matrix_temp,matrix) ! matrix_temp is dellocated
        ELSE
          STOP 'error in the variables for increase matrix size.'
        ENDIF
      ENDIF
    ENDSUBROUTINE increase_martix_int
    
    SUBROUTINE increase_martix_real(matrix,name_sub,ndim0,ndim,double_size)
      IMPLICIT NONE
      
      Real(kind=Rkind),allocatable,intent(inout)   :: matrix(:,:)
      Integer,optional,            intent(in)      :: ndim0
      Integer,optional,            intent(in)      :: ndim
      Logical,optional,            intent(in)      :: double_size
      Character(len=*),            intent(in)      :: name_sub
      
      Real(kind=Rkind),allocatable                 :: matrix_temp(:,:)
      Integer                                      :: new_dim

      IF(.NOT. allocated(matrix)) THEN
        CALL alloc_NParray(matrix,(/ ndim,ndim /),"increase matrix",name_sub)
      ELSE
        IF(present(ndim0)) THEN
          IF(present(ndim)) THEN
            IF(ndim0>=ndim) write(*,*) 'error increase matrix, ndim<=ndim0'
            new_dim=ndim
          ELSE IF(present(double_size)) THEN
            new_dim=ndim0*2
          ELSE
            STOP 'error in the variables for increase matrix size.'
          ENDIF
          
          CALL alloc_NParray(matrix_temp,(/ ndim,ndim /),"increase matrix",name_sub)
          matrix_temp(1:ndim0,1:ndim0)=matrix(1:ndim0,1:ndim0)
          CALL move_alloc(matrix_temp,matrix) ! matrix_temp is dellocated
        ELSE
          STOP 'error in the variables for increase matrix size.'
        ENDIF
      ENDIF
    ENDSUBROUTINE increase_martix_real
!---------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------
!> interface: allocate_array    
!---------------------------------------------------------------------------------------
    SUBROUTINE allocate_array_int4_length4(array_in,d1_1,d1_2)
      IMPLICIT NONE
      Integer*4,allocatable,        intent(inout) :: array_in(:)
      Integer*4,                    intent(in)    :: d1_1
      Integer*4,                    intent(in)    :: d1_2
      
      IF(allocated(array_in)) deallocate(array_in)
      
      allocate(array_in(d1_1:d1_2))
      array_in=0
      
    ENDSUBROUTINE allocate_array_int4_length4
    
    SUBROUTINE allocate_matrix_int4_length4(array_in,d1_1,d1_2,d2_1,d2_2)
      IMPLICIT NONE
      Integer*4,allocatable,        intent(inout) :: array_in(:,:)
      Integer*4,                    intent(in)    :: d1_1
      Integer*4,                    intent(in)    :: d1_2
      Integer*4,optional,           intent(in)    :: d2_1
      Integer*4,optional,           intent(in)    :: d2_2
      
      IF(allocated(array_in)) deallocate(array_in)
      
      allocate(array_in(d1_1:d1_2,d2_1:d2_2))
      array_in=0
      
    ENDSUBROUTINE allocate_matrix_int4_length4
    
    SUBROUTINE allocate_array_int4_length8(array_in,d1_1,d1_2)
      IMPLICIT NONE
      Integer*4,allocatable,        intent(inout) :: array_in(:)
      Integer*8,                    intent(in)    :: d1_1
      Integer*8,                    intent(in)    :: d1_2

      IF(allocated(array_in)) deallocate(array_in)

      allocate(array_in(d1_1:d1_2))
      array_in=0

    ENDSUBROUTINE allocate_array_int4_length8
    
    SUBROUTINE allocate_matrix_int4_length8(array_in,d1_1,d1_2,d2_1,d2_2)
      IMPLICIT NONE
      Integer*4,allocatable,        intent(inout) :: array_in(:,:)
      Integer*8,                    intent(in)    :: d1_1
      Integer*8,                    intent(in)    :: d1_2
      Integer*8,optional,           intent(in)    :: d2_1
      Integer*8,optional,           intent(in)    :: d2_2

      IF(allocated(array_in)) deallocate(array_in)
      
      allocate(array_in(d1_1:d1_2,d2_1:d2_2))
      array_in=0
      
    ENDSUBROUTINE allocate_matrix_int4_length8
    
    SUBROUTINE allocate_array_int8_length4(array_in,d1_1,d1_2)
      IMPLICIT NONE
      Integer*8,allocatable,        intent(inout) :: array_in(:)
      Integer*4,                    intent(in)    :: d1_1
      Integer*4,                    intent(in)    :: d1_2

      IF(allocated(array_in)) deallocate(array_in)

      allocate(array_in(d1_1:d1_2))
      array_in=0

    ENDSUBROUTINE allocate_array_int8_length4
    
    SUBROUTINE allocate_matrix_int8_length4(array_in,d1_1,d1_2,d2_1,d2_2)
      IMPLICIT NONE
      Integer*8,allocatable,        intent(inout) :: array_in(:,:)
      Integer*4,                    intent(in)    :: d1_1
      Integer*4,                    intent(in)    :: d1_2
      Integer*4,optional,           intent(in)    :: d2_1
      Integer*4,optional,           intent(in)    :: d2_2

      IF(allocated(array_in)) deallocate(array_in)
      
      allocate(array_in(d1_1:d1_2,d2_1:d2_2))
      array_in=0
      
    ENDSUBROUTINE allocate_matrix_int8_length4
    
    SUBROUTINE allocate_array_int8_length8(array_in,d1_1,d1_2)
      IMPLICIT NONE
      Integer*8,allocatable,        intent(inout) :: array_in(:)
      Integer*8,                    intent(in)    :: d1_1
      Integer*8,                    intent(in)    :: d1_2

      IF(allocated(array_in)) deallocate(array_in)

      allocate(array_in(d1_1:d1_2))
      array_in=0
      
    ENDSUBROUTINE allocate_array_int8_length8
    
    SUBROUTINE allocate_matrix_int8_length8(array_in,d1_1,d1_2,d2_1,d2_2)
      IMPLICIT NONE
      Integer*8,allocatable,        intent(inout) :: array_in(:,:)
      Integer*8,                    intent(in)    :: d1_1
      Integer*8,                    intent(in)    :: d1_2
      Integer*8,optional,           intent(in)    :: d2_1
      Integer*8,optional,           intent(in)    :: d2_2

      IF(allocated(array_in)) deallocate(array_in)
      
      allocate(array_in(d1_1:d1_2,d2_1:d2_2))      
      array_in=0
      
    ENDSUBROUTINE allocate_matrix_int8_length8
    
    SUBROUTINE allocate_array_real_length4(array_in,d1_1,d1_2)
      IMPLICIT NONE
      Real(kind=Rkind),allocatable, intent(inout) :: array_in(:)
      Integer*4,                    intent(in)    :: d1_1
      Integer*4,                    intent(in)    :: d1_2

      IF(allocated(array_in)) deallocate(array_in)
      
      allocate(array_in(d1_1:d1_2))
      array_in=0.
      
    ENDSUBROUTINE allocate_array_real_length4
    
    SUBROUTINE allocate_matrix_real_length4(array_in,d1_1,d1_2,d2_1,d2_2)
      IMPLICIT NONE
      Real(kind=Rkind),allocatable, intent(inout) :: array_in(:,:)
      Integer*4,                    intent(in)    :: d1_1
      Integer*4,                    intent(in)    :: d1_2
      Integer*4,optional,           intent(in)    :: d2_1
      Integer*4,optional,           intent(in)    :: d2_2

      IF(allocated(array_in)) deallocate(array_in)
      
      allocate(array_in(d1_1:d1_2,d2_1:d2_2))
      array_in=0.

    ENDSUBROUTINE allocate_matrix_real_length4

    SUBROUTINE allocate_array_real_length8(array_in,d1_1,d1_2)
      IMPLICIT NONE
      Real(kind=Rkind),allocatable, intent(inout) :: array_in(:)
      Integer*8,                    intent(in)    :: d1_1
      Integer*8,                    intent(in)    :: d1_2

      IF(allocated(array_in)) deallocate(array_in)
      
      allocate(array_in(d1_1:d1_2))
      array_in=0.
      
    ENDSUBROUTINE allocate_array_real_length8
    
    SUBROUTINE allocate_matrix_real_length8(array_in,d1_1,d1_2,d2_1,d2_2)
      IMPLICIT NONE
      Real(kind=Rkind),allocatable, intent(inout) :: array_in(:,:)
      Integer*8,                    intent(in)    :: d1_1
      Integer*8,                    intent(in)    :: d1_2
      Integer*8,optional,           intent(in)    :: d2_1
      Integer*8,optional,           intent(in)    :: d2_2

      allocate(array_in(d1_1:d1_2,d2_1:d2_2))
      array_in=0.

    ENDSUBROUTINE allocate_matrix_real_length8
    
    SUBROUTINE allocate_array_cplx_length4(array_in,d1_1,d1_2)
      IMPLICIT NONE
      Complex(kind=Rkind),allocatable,intent(inout) :: array_in(:)
      Integer*4,                    intent(in)      :: d1_1
      Integer*4,                    intent(in)      :: d1_2

      IF(allocated(array_in)) deallocate(array_in)

      allocate(array_in(d1_1:d1_2))
      array_in=0.
      
    ENDSUBROUTINE allocate_array_cplx_length4
    
    SUBROUTINE allocate_matrix_cplx_length4(array_in,d1_1,d1_2,d2_1,d2_2)
      IMPLICIT NONE
      Complex(kind=Rkind),allocatable,intent(inout) :: array_in(:,:)
      Integer*4,                    intent(in)      :: d1_1
      Integer*4,                    intent(in)      :: d1_2
      Integer*4,optional,           intent(in)      :: d2_1
      Integer*4,optional,           intent(in)      :: d2_2

      IF(allocated(array_in)) deallocate(array_in)
      
      allocate(array_in(d1_1:d1_2,d2_1:d2_2))      
      array_in=0.
      
    ENDSUBROUTINE allocate_matrix_cplx_length4

    SUBROUTINE allocate_array_cplx_length8(array_in,d1_1,d1_2)
      IMPLICIT NONE
      Complex(kind=Rkind),allocatable,intent(inout) :: array_in(:)
      Integer*8,                    intent(in)      :: d1_1
      Integer*8,                    intent(in)      :: d1_2

      IF(allocated(array_in)) deallocate(array_in)
      
      allocate(array_in(d1_1:d1_2))
      array_in=0.
      
    ENDSUBROUTINE allocate_array_cplx_length8

    SUBROUTINE allocate_array_logic(array_in,d1_1,d1_2)
      IMPLICIT NONE
      Logical,allocatable,            intent(inout) :: array_in(:)
      Integer,                           intent(in) :: d1_1
      Integer,                           intent(in) :: d1_2

      IF(allocated(array_in)) deallocate(array_in)
      
      allocate(array_in(d1_1:d1_2))
      array_in=.FALSE.
      
    ENDSUBROUTINE allocate_array_logic

    SUBROUTINE allocate_matrix_cplx_length8(array_in,d1_1,d1_2,d2_1,d2_2)
      IMPLICIT NONE
      Complex(kind=Rkind),allocatable,intent(inout) :: array_in(:,:)
      Integer*8,                    intent(in)      :: d1_1
      Integer*8,                    intent(in)      :: d1_2
      Integer*8,optional,           intent(in)      :: d2_1
      Integer*8,optional,           intent(in)      :: d2_2

      IF(allocated(array_in)) deallocate(array_in)

      allocate(array_in(d1_1:d1_2,d2_1:d2_2))
      array_in=0.

    ENDSUBROUTINE allocate_matrix_cplx_length8
!---------------------------------------------------------------------------------------
!< interface: MPI_Bcast_matrix
!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_Bcast_matrix_real(matrix,d1_1,d1_2,d2_1,d2_2,source,shift1,shift2)
      IMPLICIT NONE
      
      Real(kind=Rkind),intent(inout)      :: matrix(:,:)
      Integer,                 intent(in) :: d1_1
      Integer,                 intent(in) :: d1_2
      Integer,                 intent(in) :: d2_1
      Integer,                 intent(in) :: d2_2
      Integer,                 intent(in) :: source
      Integer,Optional,        intent(in) :: shift1
      Integer,Optional,        intent(in) :: shift2
      
      Real(kind=Rkind),allocatable        :: array(:)
      Integer                             :: length
      Integer                             :: d1_l
      Integer                             :: d1_u
      Integer                             :: d2_l
      Integer                             :: d2_u
      Integer                             :: ii
      Integer                             :: jj
      Integer                             :: kk

      d1_l=d1_1
      d1_u=d1_2
      d2_l=d2_1
      d2_u=d2_2

      IF(present(shift1)) THEN
        d1_l=d1_l+shift1
        d1_u=d1_u+shift1
      ENDIF

      IF(present(shift2)) THEN
        d2_l=d2_l+shift2
        d2_u=d2_u+shift2
      ENDIF

      length=(d1_u-d1_l+1)*(d2_u-d2_l+1)
      allocate(array(length))
      IF(MPI_id==source) THEN
        kk=0
        DO ii=d2_l,d2_u
          DO jj=d1_l,d1_u
            kk=kk+1
            array(kk)=matrix(jj,ii)
          ENDDO 
        ENDDO
      ENDIF
      
      CALL MPI_Bcast(array,length,MPI_Real8,source,MPI_COMM_WORLD,MPI_err)
      
      IF(MPI_id/=source) THEN
        kk=0
        DO ii=d2_l,d2_u
          DO jj=d1_l,d1_u
            kk=kk+1
            matrix(jj,ii)=array(kk)
          ENDDO 
        ENDDO
      ENDIF
      
      deallocate(array)
      
    ENDSUBROUTINE MPI_Bcast_matrix_real

!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_Bcast_matrix_int(matrix,d1_1,d1_2,d2_1,d2_2,source,shift1,shift2)
      IMPLICIT NONE
      
      Integer,intent(inout)               :: matrix(:,:)
      Integer,                 intent(in) :: d1_1
      Integer,                 intent(in) :: d1_2
      Integer,                 intent(in) :: d2_1
      Integer,                 intent(in) :: d2_2
      Integer,                 intent(in) :: source
      Integer,Optional,        intent(in) :: shift1
      Integer,Optional,        intent(in) :: shift2
      
      Integer,allocatable                 :: array(:)
      Integer                             :: length
      Integer                             :: d1_l
      Integer                             :: d1_u
      Integer                             :: d2_l
      Integer                             :: d2_u
      Integer                             :: ii
      Integer                             :: jj
      Integer                             :: kk

      d1_l=d1_1
      d1_u=d1_2
      d2_l=d2_1
      d2_u=d2_2

      IF(present(shift1)) THEN
        d1_l=d1_l+shift1
        d1_u=d1_u+shift1
      ENDIF

      IF(present(shift2)) THEN
        d2_l=d2_l+shift2
        d2_u=d2_u+shift2
      ENDIF

      length=(d1_u-d1_l+1)*(d2_u-d2_l+1)
      allocate(array(length))
      IF(MPI_id==source) THEN
        kk=0
        DO ii=d2_l,d2_u
          DO jj=d1_l,d1_u
            kk=kk+1
            array(kk)=matrix(jj,ii)
          ENDDO 
        ENDDO
      ENDIF
      
      CALL MPI_Bcast(array,length,MPI_int_fortran,source,MPI_COMM_WORLD,MPI_err)
      
      IF(MPI_id/=source) THEN
        kk=0
        DO ii=d2_l,d2_u
          DO jj=d1_l,d1_u
            kk=kk+1
            matrix(jj,ii)=array(kk)
          ENDDO 
        ENDDO
      ENDIF
      
      deallocate(array)
      
    ENDSUBROUTINE MPI_Bcast_matrix_int
    
!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_Bcast_matrix_complex(matrix,d1_1,d1_2,d2_1,d2_2,source,shift1,shift2)
      IMPLICIT NONE
      
      Complex(kind=Rkind),  intent(inout) :: matrix(:,:)
      Integer,                 intent(in) :: d1_1
      Integer,                 intent(in) :: d1_2
      Integer,                 intent(in) :: d2_1
      Integer,                 intent(in) :: d2_2
      Integer,                 intent(in) :: source
      Integer,Optional,        intent(in) :: shift1
      Integer,Optional,        intent(in) :: shift2

      Complex(kind=Rkind),allocatable     :: array(:)
      Integer                             :: length
      Integer                             :: d1_l
      Integer                             :: d1_u
      Integer                             :: d2_l
      Integer                             :: d2_u
      Integer                             :: ii
      Integer                             :: jj
      Integer                             :: kk

      d1_l=d1_1
      d1_u=d1_2
      d2_l=d2_1
      d2_u=d2_2

      IF(present(shift1)) THEN
        d1_l=d1_l+shift1
        d1_u=d1_u+shift1
      ENDIF

      IF(present(shift2)) THEN
        d2_l=d2_l+shift2
        d2_u=d2_u+shift2
      ENDIF

      length=(d1_u-d1_l+1)*(d2_u-d2_l+1)
      allocate(array(length))
      IF(MPI_id==source) THEN
        kk=0
        DO ii=d2_l,d2_u
          DO jj=d1_l,d1_u
            kk=kk+1
            array(kk)=matrix(jj,ii)
          ENDDO
        ENDDO
      ENDIF
      
      CALL MPI_Bcast(array,length,MPI_Complex8,source,MPI_COMM_WORLD,MPI_err)
      
      IF(MPI_id/=source) THEN
        kk=0
        DO ii=d2_l,d2_u
          DO jj=d1_l,d1_u
            kk=kk+1
            matrix(jj,ii)=array(kk)
          ENDDO 
        ENDDO
      ENDIF
      
      deallocate(array)
      
    ENDSUBROUTINE MPI_Bcast_matrix_complex
!---------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------
! interface: MPI_Reduce_sum_matrix
!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_Reduce_sum_matrix_complex(matrix,d1_l,d1_u,d2_l,d2_u,destination)
      IMPLICIT NONE
      
      Complex(kind=Rkind),intent(inout)   :: matrix(:,:)
      Integer,intent(in)                  :: d1_l
      Integer,intent(in)                  :: d1_u
      Integer,intent(in)                  :: d2_l
      Integer,intent(in)                  :: d2_u
      Integer,intent(in)                  :: destination
      
      Complex(kind=Rkind),allocatable     :: array(:)
      Complex(kind=Rkind),allocatable     :: array_des(:)
      Integer                             :: length
      Integer                             :: ii
      Integer                             :: jj
      Integer                             :: kk

      length=(d1_u-d1_l+1)*(d2_u-d2_l+1)
      allocate(array(length))
      allocate(array_des(length))
      kk=0
      DO ii=d2_l,d2_u
        DO jj=d1_l,d1_u
          kk=kk+1
          array(kk)=matrix(jj,ii)
        ENDDO
      ENDDO
      
      CALL MPI_Reduce(array,array_des,length,MPI_Complex8,MPI_SUM,root_MPI,            &
                      MPI_COMM_WORLD,MPI_err)
      
      IF(MPI_id==destination) THEN
        kk=0
        DO ii=d2_l,d2_u
          DO jj=d1_l,d1_u
            kk=kk+1
            matrix(jj,ii)=array_des(kk)
          ENDDO 
        ENDDO
      ENDIF
      
      deallocate(array)
      deallocate(array_des)
      
    ENDSUBROUTINE MPI_Reduce_sum_matrix_complex
    
!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_Reduce_sum_matrix_real(matrix,d1_l,d1_u,d2_l,d2_u,destination)
      IMPLICIT NONE
      
      Real(kind=Rkind),intent(inout)      :: matrix(:,:)
      Integer,intent(in)                  :: d1_l
      Integer,intent(in)                  :: d1_u
      Integer,intent(in)                  :: d2_l
      Integer,intent(in)                  :: d2_u
      Integer,intent(in)                  :: destination
      
      Real(kind=Rkind),allocatable        :: array(:)
      Real(kind=Rkind),allocatable        :: array_des(:)
      Integer                             :: length
      Integer                             :: ii
      Integer                             :: jj
      Integer                             :: kk

      length=(d1_u-d1_l+1)*(d2_u-d2_l+1)
      allocate(array(length))
      allocate(array_des(length))
      kk=0
      DO ii=d2_l,d2_u
        DO jj=d1_l,d1_u
          kk=kk+1
          array(kk)=matrix(jj,ii)
        ENDDO
      ENDDO
      
      CALL MPI_Reduce(array,array_des,length,MPI_Real8,MPI_SUM,root_MPI,               &
                      MPI_COMM_WORLD,MPI_err)
      
      IF(MPI_id==destination) THEN
        kk=0
        DO ii=d2_l,d2_u
          DO jj=d1_l,d1_u
            kk=kk+1
            matrix(jj,ii)=array_des(kk)
          ENDDO 
        ENDDO
      ENDIF
      
      deallocate(array)
      deallocate(array_des)
      
    ENDSUBROUTINE MPI_Reduce_sum_matrix_real

!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_Reduce_sum_matrix_int(matrix,d1_l,d1_u,d2_l,d2_u,destination)
      IMPLICIT NONE
      
      Integer,intent(inout)               :: matrix(:,:)
      Integer,intent(in)                  :: d1_l
      Integer,intent(in)                  :: d1_u
      Integer,intent(in)                  :: d2_l
      Integer,intent(in)                  :: d2_u
      Integer,intent(in)                  :: destination
      
      Integer,allocatable                 :: array(:)
      Integer,allocatable                 :: array_des(:)
      Integer                             :: length
      Integer                             :: ii
      Integer                             :: jj
      Integer                             :: kk

      length=(d1_u-d1_l+1)*(d2_u-d2_l+1)
      allocate(array(length))
      allocate(array_des(length))
      kk=0
      DO ii=d2_l,d2_u
        DO jj=d1_l,d1_u
          kk=kk+1
          array(kk)=matrix(jj,ii)
        ENDDO
      ENDDO
      
      CALL MPI_Reduce(array,array_des,length,MPI_int_fortran,MPI_SUM,root_MPI,         &
                      MPI_COMM_WORLD,MPI_err)
      
      IF(MPI_id==destination) THEN
        kk=0
        DO ii=d2_l,d2_u
          DO jj=d1_l,d1_u
            kk=kk+1
            matrix(jj,ii)=array_des(kk)
          ENDDO 
        ENDDO
      ENDIF
      
      deallocate(array)
      deallocate(array_des)
      
    ENDSUBROUTINE MPI_Reduce_sum_matrix_int
!---------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------
!< interface: MPI_Reduce_sum_Bcast
!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_Reduce_sum_Bcast_real(value)
      IMPLICIT NONE
      Real(kind=Rkind),     intent(inout) :: value
      Real(kind=Rkind)                    :: value_temp
      
      CALL MPI_Reduce(value,value_temp,size1_MPI,MPI_real_fortran,MPI_SUM,root_MPI,    &
                      MPI_COMM_WORLD,MPI_err)
      IF(MPI_id==0) value=value_temp
      CALL MPI_Bcast(value,size1_MPI,MPI_real_fortran,root_MPI,MPI_COMM_WORLD,MPI_err)
      
    ENDSUBROUTINE MPI_Reduce_sum_Bcast_real
    
!---------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_Reduce_sum_Bcast_complex(value)
      IMPLICIT NONE
      Complex(kind=Rkind),  intent(inout) :: value
      Complex(kind=Rkind)                 :: value_temp
      
      CALL MPI_Reduce(value,value_temp,size1_MPI,MPI_cplx_fortran,MPI_SUM,root_MPI,    &
                      MPI_COMM_WORLD,MPI_err)
      IF(MPI_id==0) value=value_temp
      CALL MPI_Bcast(value,size1_MPI,MPI_cplx_fortran,root_MPI,MPI_COMM_WORLD,MPI_err)
      
    ENDSUBROUTINE MPI_Reduce_sum_Bcast_complex

!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_Reduce_sum_Bcast_int(value)
      IMPLICIT NONE
      Integer,              intent(inout) :: value
      Integer                             :: value_temp
      
      CALL MPI_Reduce(value,value_temp,size1_MPI,MPI_int_fortran,MPI_SUM,root_MPI,     &
                      MPI_COMM_WORLD,MPI_err)
      IF(MPI_id==0) value=value_temp
      CALL MPI_Bcast(value,size1_MPI,MPI_int_fortran,root_MPI,MPI_COMM_WORLD,MPI_err)
      
    ENDSUBROUTINE MPI_Reduce_sum_Bcast_int    
!---------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------
!< interface: MPI_Send_matrix
!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_Send_matrix_real(matrix,d1_l,d1_u,d2_l,d2_u,destination,tag)
      IMPLICIT NONE
      
      Real(kind=Rkind),intent(inout)      :: matrix(:,:)
      Integer,intent(in)                  :: d1_l
      Integer,intent(in)                  :: d1_u
      Integer,intent(in)                  :: d2_l
      Integer,intent(in)                  :: d2_u
      Integer,intent(in)                  :: destination
      Integer,intent(in)                  :: tag

      Real(kind=Rkind),allocatable        :: array(:)
      Integer                             :: length
      Integer                             :: ii
      Integer                             :: jj
      Integer                             :: kk

      length=(d1_u-d1_l+1)*(d2_u-d2_l+1)
      allocate(array(length))
      kk=0
      DO ii=d2_l,d2_u
        DO jj=d1_l,d1_u
          kk=kk+1
          array(kk)=matrix(jj,ii)
        ENDDO 
      ENDDO
      
      CALL MPI_Send(array,length,MPI_Real8,destination,tag,MPI_COMM_WORLD,MPI_err)
      
      deallocate(array)
      
    ENDSUBROUTINE MPI_Send_matrix_real

    !-----------------------------------------------------------------------------------
    SUBROUTINE MPI_Send_matrix_int(matrix,d1_l,d1_u,d2_l,d2_u,destination,tag)
      IMPLICIT NONE
      
      Integer,intent(inout)               :: matrix(:,:)
      Integer,intent(in)                  :: d1_l
      Integer,intent(in)                  :: d1_u
      Integer,intent(in)                  :: d2_l
      Integer,intent(in)                  :: d2_u
      Integer,intent(in)                  :: destination
      Integer,intent(in)                  :: tag

      Integer,allocatable                 :: array(:)
      Integer                             :: length
      Integer                             :: ii
      Integer                             :: jj
      Integer                             :: kk

      length=(d1_u-d1_l+1)*(d2_u-d2_l+1)
      allocate(array(length))
      kk=0
      DO ii=d2_l,d2_u
        DO jj=d1_l,d1_u
          kk=kk+1
          array(kk)=matrix(jj,ii)
        ENDDO 
      ENDDO
      
      CALL MPI_Send(array,length,MPI_int_fortran,destination,tag,MPI_COMM_WORLD,MPI_err)
      
      deallocate(array)
      
    ENDSUBROUTINE MPI_Send_matrix_int  

    !-----------------------------------------------------------------------------------
    SUBROUTINE MPI_Send_matrix_complex(matrix,d1_l,d1_u,d2_l,d2_u,destination,tag)
      IMPLICIT NONE
      
      Complex(kind=Rkind),intent(inout)   :: matrix(:,:)
      Integer,intent(in)                  :: d1_l
      Integer,intent(in)                  :: d1_u
      Integer,intent(in)                  :: d2_l
      Integer,intent(in)                  :: d2_u
      Integer,intent(in)                  :: destination
      Integer,intent(in)                  :: tag
            
      Complex(kind=Rkind),allocatable     :: array(:)
      Integer                             :: length
      Integer                             :: ii
      Integer                             :: jj
      Integer                             :: kk

      length=(d1_u-d1_l+1)*(d2_u-d2_l+1)
      allocate(array(length))
      kk=0
      DO ii=d2_l,d2_u
        DO jj=d1_l,d1_u
          kk=kk+1
          array(kk)=matrix(jj,ii)
        ENDDO 
      ENDDO
      
      CALL MPI_Send(array,length,MPI_Complex8,destination,tag,MPI_COMM_WORLD,MPI_err)
      
      deallocate(array)
      
    ENDSUBROUTINE MPI_Send_matrix_complex
    
    !-----------------------------------------------------------------------------------
    SUBROUTINE MPI_Recv_matrix_real(matrix,d1_l,d1_u,d2_l,d2_u,source,tag)
      IMPLICIT NONE
      
      Real(kind=Rkind),intent(inout)      :: matrix(:,:)
      Integer,intent(in)                  :: d1_l
      Integer,intent(in)                  :: d1_u
      Integer,intent(in)                  :: d2_l
      Integer,intent(in)                  :: d2_u
      Integer,intent(in)                  :: source
      Integer,intent(in)                  :: tag
      
      Real(kind=Rkind),allocatable        :: array(:)
      Integer                             :: length
      Integer                             :: ii
      Integer                             :: jj
      Integer                             :: kk

      length=(d1_u-d1_l+1)*(d2_u-d2_l+1)
      allocate(array(length))
      
      CALL MPI_Recv(array,length,MPI_Real8,source,tag,MPI_COMM_WORLD,MPI_stat,MPI_err)

      kk=0
      DO ii=d2_l,d2_u
        DO jj=d1_l,d1_u
          kk=kk+1
          matrix(jj,ii)=array(kk)
        ENDDO 
      ENDDO

      deallocate(array)
      
    ENDSUBROUTINE MPI_Recv_matrix_real
    
    !-----------------------------------------------------------------------------------
    SUBROUTINE MPI_Recv_matrix_int(matrix,d1_l,d1_u,d2_l,d2_u,source,tag)
      IMPLICIT NONE
      
      Integer,intent(inout)               :: matrix(:,:)
      Integer,intent(in)                  :: d1_l
      Integer,intent(in)                  :: d1_u
      Integer,intent(in)                  :: d2_l
      Integer,intent(in)                  :: d2_u
      Integer,intent(in)                  :: source
      Integer,intent(in)                  :: tag
      
      Integer,allocatable                 :: array(:)
      Integer                             :: length
      Integer                             :: ii
      Integer                             :: jj
      Integer                             :: kk

      length=(d1_u-d1_l+1)*(d2_u-d2_l+1)
      allocate(array(length))
      
      CALL MPI_Recv(array,length,MPI_int_fortran,source,tag,                           &
                    MPI_COMM_WORLD,MPI_stat,MPI_err)

      kk=0
      DO ii=d2_l,d2_u
        DO jj=d1_l,d1_u
          kk=kk+1
          matrix(jj,ii)=array(kk)
        ENDDO 
      ENDDO

      deallocate(array)
      
    ENDSUBROUTINE MPI_Recv_matrix_int
    
    !-----------------------------------------------------------------------------------
    SUBROUTINE MPI_Recv_matrix_complex(matrix,d1_l,d1_u,d2_l,d2_u,source,tag)
      IMPLICIT NONE
      
      Complex(kind=Rkind),intent(inout)   :: matrix(:,:)
      Integer,intent(in)                  :: d1_l
      Integer,intent(in)                  :: d1_u
      Integer,intent(in)                  :: d2_l
      Integer,intent(in)                  :: d2_u
      Integer,intent(in)                  :: source
      Integer,intent(in)                  :: tag
      
      Complex(kind=Rkind),allocatable     :: array(:)
      Integer                             :: length
      Integer                             :: ii
      Integer                             :: jj
      Integer                             :: kk

      length=(d1_u-d1_l+1)*(d2_u-d2_l+1)
      allocate(array(length))
      
      CALL MPI_Recv(array,length,MPI_Complex8,source,tag,MPI_COMM_WORLD,MPI_stat,MPI_err)

      kk=0
      DO ii=d2_l,d2_u
        DO jj=d1_l,d1_u
          kk=kk+1
          matrix(jj,ii)=array(kk)
        ENDDO 
      ENDDO

      deallocate(array)
      
    ENDSUBROUTINE MPI_Recv_matrix_complex

!---------------------------------------------------------------------------------------
!< interface: MPI_combine_array_real
!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_combine_array_real(array)
      IMPLICIT NONE
      
      Real(kind=Rkind),allocatable,   intent(inout) :: array(:)

      Integer                                       :: d1
      Integer                                       :: d2

      IF(MPI_id/=0) THEN
        d1=bounds_MPI(1,MPI_id)
        d2=bounds_MPI(2,MPI_id)
        CALL MPI_Send(array(d1:d2),d2-d1+1,MPI_real_fortran,root_MPI,MPI_id,           &
                     MPI_COMM_WORLD,MPI_err)
      ENDIF

      IF(MPI_id==0) THEN
        DO i_MPI=1,MPI_np-1
          d1=bounds_MPI(1,i_MPI)
          d2=bounds_MPI(2,i_MPI)
          CALL MPI_Recv(array(d1:d2),d2-d1+1,MPI_real_fortran,i_MPI,                   &
                        i_MPI,MPI_COMM_WORLD,MPI_stat,MPI_err)
        ENDDO
      ENDIF
    ENDSUBROUTINE MPI_combine_array_real

    !-----------------------------------------------------------------------------------
    SUBROUTINE MPI_combine_array_int(array)
      IMPLICIT NONE
      
      Integer,allocatable,intent(inout)             :: array(:)

      Integer                                       :: d1
      Integer                                       :: d2

      IF(MPI_id/=0) THEN
        d1=bounds_MPI(1,MPI_id)
        d2=bounds_MPI(2,MPI_id)
        CALL MPI_Send(array(d1:d2),d2-d1+1,MPI_int_fortran,root_MPI,MPI_id,            &
                     MPI_COMM_WORLD,MPI_err)
      ENDIF

      IF(MPI_id==0) THEN
        DO i_MPI=1,MPI_np-1
          d1=bounds_MPI(1,i_MPI)
          d2=bounds_MPI(2,i_MPI)
          CALL MPI_Recv(array(d1:d2),d2-d1+1,MPI_int_fortran,i_MPI,                    &
                        i_MPI,MPI_COMM_WORLD,MPI_stat,MPI_err)
        ENDDO
      ENDIF
    ENDSUBROUTINE MPI_combine_array_int
    
    !-----------------------------------------------------------------------------------
    SUBROUTINE MPI_combine_array_cplx(array)
      IMPLICIT NONE
      
      Complex(kind=Rkind),allocatable,intent(inout) :: array(:)

      Integer                                       :: d1
      Integer                                       :: d2

      IF(MPI_id/=0) THEN
        d1=bounds_MPI(1,MPI_id)
        d2=bounds_MPI(2,MPI_id)
        CALL MPI_Send(array(d1:d2),d2-d1+1,MPI_cplx_fortran,root_MPI,MPI_id,           &
                     MPI_COMM_WORLD,MPI_err)
      ENDIF

      IF(MPI_id==0) THEN
        DO i_MPI=1,MPI_np-1
          d1=bounds_MPI(1,i_MPI)
          d2=bounds_MPI(2,i_MPI)
          CALL MPI_Recv(array(d1:d2),d2-d1+1,MPI_cplx_fortran,i_MPI,                   &
                        i_MPI,MPI_COMM_WORLD,MPI_stat,MPI_err)
        ENDDO
      ENDIF
    ENDSUBROUTINE MPI_combine_array_cplx
!---------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------
!< interface: MPI_collect_info_int
!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_collect_info_int(array)
      IMPLICIT NONE

      Integer,                        intent(inout) :: array(0:MPI_np-1)

      IF(MPI_id/=0) THEN
        CALL MPI_Send(array(MPI_id),size1_MPI,MPI_int_fortran,root_MPI,MPI_id,         &
                      MPI_COMM_WORLD,MPI_err)
      ENDIF

      IF(MPI_id==0) THEN
        DO i_MPI=1,MPI_np-1
          CALL MPI_Recv(array(i_MPI),size1_MPI,MPI_int_fortran,i_MPI,                  &
                        i_MPI,MPI_COMM_WORLD,MPI_stat,MPI_err)
        ENDDO
      ENDIF

!      CALL MPI_Bcast(array(0:MPI_np-1),MPI_np,MPI_int_fortran,root_MPI,                &
!                     MPI_COMM_WORLD,MPI_err)

    ENDSUBROUTINE MPI_collect_info_int
    
    !-----------------------------------------------------------------------------------
    SUBROUTINE MPI_collect_info_real(array)
      IMPLICIT NONE

      Real(kind=Rkind),               intent(inout) :: array(0:MPI_np-1)

      IF(MPI_id/=0) THEN
        CALL MPI_Send(array(MPI_id),size1_MPI,MPI_real_fortran,root_MPI,MPI_id,        &
                      MPI_COMM_WORLD,MPI_err)
      ENDIF

      IF(MPI_id==0) THEN
        DO i_MPI=1,MPI_np-1
          CALL MPI_Recv(array(i_MPI),size1_MPI,MPI_real_fortran,i_MPI,                 &
                        i_MPI,MPI_COMM_WORLD,MPI_stat,MPI_err)
        ENDDO
      ENDIF

!      CALL MPI_Bcast(array(0:MPI_np-1),MPI_np,MPI_real_fortran,root_MPI,               &
!                     MPI_COMM_WORLD,MPI_err)

    ENDSUBROUTINE MPI_collect_info_real

#endif

END MODULE mod_MPI_Aid


