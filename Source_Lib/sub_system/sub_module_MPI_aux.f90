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
MODULE mod_MPI_aux
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
  
!  INTERFACE MPI__Send
!    module procedure MPI__Send_int4
!    module procedure MPI__Send_int8
!    module procedure MPI__Send_real
!    module procedure MPI__Send_cplx
!  END INTERFACE
  
  
  INTERFACE MPI_Bcast_
    module procedure MPI_Bcast_array_int4
    module procedure MPI_Bcast_array_int8
    module procedure MPI_Bcast_array_real
    module procedure MPI_Bcast_array_cplx
    module procedure MPI_Bcast_array_logical
    module procedure MPI_Bcast_variable_int4
    module procedure MPI_Bcast_variable_int8
    module procedure MPI_Bcast_variable_real
    module procedure MPI_Bcast_variable_cplx
    module procedure MPI_Bcast_variable_logical
  END INTERFACE

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
    module procedure MPI_Reduce_sum_Bcast_array_int
    module procedure MPI_Reduce_sum_Bcast_array_real
  END INTERFACE

  INTERFACE MPI_Reduce_max_Bcast
    module procedure MPI_Reduce_max_Bcast_real
    module procedure MPI_Reduce_max_Bcast_int
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
    module procedure MPI_combine_array_int4
    module procedure MPI_combine_array_int8
    module procedure MPI_combine_array_cplx
    module procedure MPI_combine_array_general_cplx
  END INTERFACE
  
  INTERFACE MPI_collect_info
    module procedure MPI_collect_info_int
    module procedure MPI_collect_info_real
  END INTERFACE

!---------------------------------------------------------------------------------------
  Contains

!---------------------------------------------------------------------------------------
! we need a better realtime memory check subroutine.
!---------------------------------------------------------------------------------------
    !> check total memory used at certain point
    SUBROUTINE system_mem_usage(memory_RSS,name)
      USE mod_NumParameters
#if(run_MPI_ifort)
      USE ifport ! if on intel compiler
#endif
      IMPLICIT NONE
      Integer, intent(out) :: memory_RSS
      Character(len=200):: filename=' '
      Character(len=80) :: line
      Character(len=8)  :: pid_char=' '
      Integer :: pid
      Logical :: ifxst
      Character (len=*), intent(in) :: name

#if(run_MPI)

      memory_RSS=-1 ! return negative number if not found

      !> get process ID
      pid=getpid()
      !write(*,*) 'pid=',pid
      write(pid_char,'(I8)') pid
      filename='/proc/'//trim(adjustl(pid_char))//'/status'

      ! read system file
      inquire (file=filename,exist=ifxst)
      IF(.not.ifxst) THEN
        !write (*,*) 'system file does not exist'
      ELSE
        OPEN(unit=100, file=filename, action='read')
        DO
          read(100,'(a)',end=120) line
          IF(line(1:6).eq.'VmRSS:') THEN
            read (line(7:),*) memory_RSS
            EXIT
          ENDIF
        ENDDO
120     CONTINUE
        CLOSE(100)
        write(out_unitp,121) name,memory_RSS,MPI_id    ! kB
121     format('memory check at ',a,': ',i4,' from ',i4)
      ENDIF

#endif
    ENDSUBROUTINE system_mem_usage

!---------------------------------------------------------------------------------------
!> Modified MPI_send
!---------------------------------------------------------------------------------------
!    SUBROUTINE MPI__Send_int4(array,length,destination,tag)
!      Integer(kind=Ikind),              intent(in) :: array(:)
!      Integer,                          intent(in) :: length
!      Integer(kind=MPI_INTEGER_KIND),   intent(in) :: destination
!      Integer(kind=MPI_INTEGER_KIND),   intent(in) :: tag
!
!      IF(length<=huge(0_4)) THEN
!        CALL MPI_Send(array,length,MPI_Integer4,destination,tag,MPI_COMM_WORLD,MPI_err)
!      ELSE
!        
!      ENDIF
!    ENDSUBROUTINE MPI__Send_int4

!---------------------------------------------------------------------------------------
!> interface: MPI_Bcast_
!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_Bcast_array_int4(array,length,source)
      IMPLICIT NONE

      Integer(kind=Ikind),              intent(inout) :: array(:)
      Integer(kind=MPI_INTEGER_KIND),   intent(in)    :: length
      Integer(kind=MPI_INTEGER_KIND),   intent(in)    :: source

#if(run_MPI)

      CALL MPI_Bcast(array,length,MPI_Integer4,source,MPI_COMM_WORLD,MPI_err)

#endif
    ENDSUBROUTINE MPI_Bcast_array_int4
!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_Bcast_array_int8(array,length,source)
      IMPLICIT NONE

      Integer(kind=ILkind),             intent(inout) :: array(:)
      Integer(kind=MPI_INTEGER_KIND),   intent(in)    :: length
      Integer(kind=MPI_INTEGER_KIND),   intent(in)    :: source

#if(run_MPI)

      CALL MPI_Bcast(array,length,MPI_Integer8,source,MPI_COMM_WORLD,MPI_err)

#endif
    ENDSUBROUTINE MPI_Bcast_array_int8
!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_Bcast_array_real(array,length,source)
      IMPLICIT NONE

      Real(kind=Rkind),                 intent(inout) :: array(:)
      Integer(kind=MPI_INTEGER_KIND),   intent(in)    :: length
      Integer(kind=MPI_INTEGER_KIND),   intent(in)    :: source

#if(run_MPI)

      CALL MPI_Bcast(array,length,Real_MPI,source,MPI_COMM_WORLD,MPI_err)

#endif
    ENDSUBROUTINE MPI_Bcast_array_real
!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_Bcast_array_cplx(array,length,source)
      IMPLICIT NONE

      Complex(kind=Rkind),              intent(inout) :: array(:)
      Integer(kind=MPI_INTEGER_KIND),   intent(in)    :: length
      Integer(kind=MPI_INTEGER_KIND),   intent(in)    :: source

#if(run_MPI)

      CALL MPI_Bcast(array,length,cplx_MPI,source,MPI_COMM_WORLD,MPI_err)

#endif
    ENDSUBROUTINE MPI_Bcast_array_cplx
!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_Bcast_array_logical(array,length,source)
      IMPLICIT NONE

      Logical,                          intent(inout) :: array(:)
      Integer(kind=MPI_INTEGER_KIND),   intent(in)    :: length
      Integer(kind=MPI_INTEGER_KIND),   intent(in)    :: source

#if(run_MPI)

      CALL MPI_Bcast(array,length,MPI_logical,source,MPI_COMM_WORLD,MPI_err)

#endif
    ENDSUBROUTINE MPI_Bcast_array_logical
!---------------------------------------------------------------------------------------

    SUBROUTINE MPI_Bcast_variable_int4(variable,length,source)
      IMPLICIT NONE

      Integer(kind=Ikind),              intent(inout) :: variable
      Integer(kind=MPI_INTEGER_KIND),   intent(in)    :: length
      Integer(kind=MPI_INTEGER_KIND),   intent(in)    :: source

#if(run_MPI)

      CALL MPI_Bcast(variable,size1_MPI,MPI_Integer4,source,MPI_COMM_WORLD,MPI_err)

#endif
    ENDSUBROUTINE MPI_Bcast_variable_int4
!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_Bcast_variable_int8(variable,length,source)
      IMPLICIT NONE

      Integer(kind=ILkind),             intent(inout) :: variable
      Integer(kind=MPI_INTEGER_KIND),   intent(in)    :: length
      Integer(kind=MPI_INTEGER_KIND),   intent(in)    :: source

#if(run_MPI)

      CALL MPI_Bcast(variable,length,MPI_Integer8,source,MPI_COMM_WORLD,MPI_err)

#endif
    ENDSUBROUTINE MPI_Bcast_variable_int8
!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_Bcast_variable_real(variable,length,source)
      IMPLICIT NONE

      Real(kind=Rkind),                 intent(inout) :: variable
      Integer(kind=MPI_INTEGER_KIND),   intent(in)    :: length
      Integer(kind=MPI_INTEGER_KIND),   intent(in)    :: source

#if(run_MPI)

      CALL MPI_Bcast(variable,length,Real_MPI,source,MPI_COMM_WORLD,MPI_err)

#endif
    ENDSUBROUTINE MPI_Bcast_variable_real
!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_Bcast_variable_cplx(variable,length,source)
      IMPLICIT NONE

      Complex(kind=Rkind),              intent(inout) :: variable
      Integer(kind=MPI_INTEGER_KIND),   intent(in)    :: length
      Integer(kind=MPI_INTEGER_KIND),   intent(in)    :: source

#if(run_MPI)

      CALL MPI_Bcast(variable,length,Cplx_MPI,source,MPI_COMM_WORLD,MPI_err)

#endif
    ENDSUBROUTINE MPI_Bcast_variable_cplx
!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_Bcast_variable_logical(variable,length,source)
      IMPLICIT NONE

      Logical,                          intent(inout) :: variable
      Integer(kind=MPI_INTEGER_KIND),   intent(in)    :: length
      Integer(kind=MPI_INTEGER_KIND),   intent(in)    :: source

#if(run_MPI)

      CALL MPI_Bcast(variable,length,MPI_Logical,source,MPI_COMM_WORLD,MPI_err)

#endif
    ENDSUBROUTINE MPI_Bcast_variable_logical
!---------------------------------------------------------------------------------------


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
            IF(ndim0>=ndim) write(out_unitp,*) 'error increase matrix, ndim<=ndim0'
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
      Integer(kind=Ikind),allocatable, intent(inout) :: array_in(:)
      Integer(kind=Ikind),                intent(in) :: d1_1
      Integer(kind=Ikind),                intent(in) :: d1_2
      
      IF(allocated(array_in)) deallocate(array_in)
      
      allocate(array_in(d1_1:d1_2))
      array_in=0
      
    ENDSUBROUTINE allocate_array_int4_length4
    
    SUBROUTINE allocate_matrix_int4_length4(array_in,d1_1,d1_2,d2_1,d2_2)
      IMPLICIT NONE
      Integer(kind=Ikind),allocatable, intent(inout) :: array_in(:,:)
      Integer(kind=Ikind),                intent(in) :: d1_1
      Integer(kind=Ikind),                intent(in) :: d1_2
      Integer(kind=Ikind),                intent(in) :: d2_1
      Integer(kind=Ikind),                intent(in) :: d2_2
      
      IF(allocated(array_in)) deallocate(array_in)
      
      allocate(array_in(d1_1:d1_2,d2_1:d2_2))
      array_in=0
      
    ENDSUBROUTINE allocate_matrix_int4_length4
    
    SUBROUTINE allocate_array_int4_length8(array_in,d1_1,d1_2)
      IMPLICIT NONE
      Integer(kind=Ikind),allocatable, intent(inout) :: array_in(:)
      Integer(kind=ILkind),               intent(in) :: d1_1
      Integer(kind=ILkind),               intent(in) :: d1_2

      IF(allocated(array_in)) deallocate(array_in)

      allocate(array_in(d1_1:d1_2))
      array_in=0

    ENDSUBROUTINE allocate_array_int4_length8
    
    SUBROUTINE allocate_matrix_int4_length8(array_in,d1_1,d1_2,d2_1,d2_2)
      IMPLICIT NONE
      Integer(kind=Ikind),allocatable, intent(inout) :: array_in(:,:)
      Integer(kind=ILkind),               intent(in) :: d1_1
      Integer(kind=ILkind),               intent(in) :: d1_2
      Integer(kind=ILkind),               intent(in) :: d2_1
      Integer(kind=ILkind),               intent(in) :: d2_2

      IF(allocated(array_in)) deallocate(array_in)
      
      allocate(array_in(d1_1:d1_2,d2_1:d2_2))
      array_in=0
      
    ENDSUBROUTINE allocate_matrix_int4_length8
    
    SUBROUTINE allocate_array_int8_length4(array_in,d1_1,d1_2)
      IMPLICIT NONE
      Integer(kind=ILkind),allocatable,intent(inout) :: array_in(:)
      Integer(kind=Ikind),                intent(in) :: d1_1
      Integer(kind=Ikind),                intent(in) :: d1_2

      IF(allocated(array_in)) deallocate(array_in)

      allocate(array_in(d1_1:d1_2))
      array_in=0

    ENDSUBROUTINE allocate_array_int8_length4
    
    SUBROUTINE allocate_matrix_int8_length4(array_in,d1_1,d1_2,d2_1,d2_2)
      IMPLICIT NONE
      Integer(kind=ILkind),allocatable,intent(inout) :: array_in(:,:)
      Integer(kind=Ikind),                intent(in) :: d1_1
      Integer(kind=Ikind),                intent(in) :: d1_2
      Integer(kind=Ikind),                intent(in) :: d2_1
      Integer(kind=Ikind),                intent(in) :: d2_2

      IF(allocated(array_in)) deallocate(array_in)
      
      allocate(array_in(d1_1:d1_2,d2_1:d2_2))
      array_in=0
      
    ENDSUBROUTINE allocate_matrix_int8_length4
    
    SUBROUTINE allocate_array_int8_length8(array_in,d1_1,d1_2)
      IMPLICIT NONE
      Integer(kind=ILkind),allocatable,intent(inout) :: array_in(:)
      Integer(kind=ILkind),               intent(in) :: d1_1
      Integer(kind=ILkind),               intent(in) :: d1_2

      IF(allocated(array_in)) deallocate(array_in)

      allocate(array_in(d1_1:d1_2))
      array_in=0
      
    ENDSUBROUTINE allocate_array_int8_length8
    
    SUBROUTINE allocate_matrix_int8_length8(array_in,d1_1,d1_2,d2_1,d2_2)
      IMPLICIT NONE
      Integer(kind=ILkind),allocatable,intent(inout) :: array_in(:,:)
      Integer(kind=ILkind),               intent(in) :: d1_1
      Integer(kind=ILkind),               intent(in) :: d1_2
      Integer(kind=ILkind),               intent(in) :: d2_1
      Integer(kind=ILkind),               intent(in) :: d2_2

      IF(allocated(array_in)) deallocate(array_in)
      
      allocate(array_in(d1_1:d1_2,d2_1:d2_2))      
      array_in=0
      
    ENDSUBROUTINE allocate_matrix_int8_length8
    
    SUBROUTINE allocate_array_real_length4(array_in,d1_1,d1_2)
      IMPLICIT NONE
      Real(kind=Rkind),allocatable,    intent(inout) :: array_in(:)
      Integer(kind=Ikind),                intent(in) :: d1_1
      Integer(kind=Ikind),                intent(in) :: d1_2

      IF(allocated(array_in)) deallocate(array_in)
      
      allocate(array_in(d1_1:d1_2))
      array_in=0.
      
    ENDSUBROUTINE allocate_array_real_length4
    
    SUBROUTINE allocate_matrix_real_length4(array_in,d1_1,d1_2,d2_1,d2_2)
      IMPLICIT NONE
      Real(kind=Rkind),allocatable,    intent(inout) :: array_in(:,:)
      Integer(kind=Ikind),                intent(in) :: d1_1
      Integer(kind=Ikind),                intent(in) :: d1_2
      Integer(kind=Ikind),                intent(in) :: d2_1
      Integer(kind=Ikind),                intent(in) :: d2_2

      IF(allocated(array_in)) deallocate(array_in)
      
      allocate(array_in(d1_1:d1_2,d2_1:d2_2))
      array_in=0.

    ENDSUBROUTINE allocate_matrix_real_length4

    SUBROUTINE allocate_array_real_length8(array_in,d1_1,d1_2)
      IMPLICIT NONE
      Real(kind=Rkind),allocatable,    intent(inout) :: array_in(:)
      Integer(kind=ILkind),               intent(in) :: d1_1
      Integer(kind=ILkind),               intent(in) :: d1_2

      IF(allocated(array_in)) deallocate(array_in)
      
      allocate(array_in(d1_1:d1_2))
      array_in=0.
      
    ENDSUBROUTINE allocate_array_real_length8
    
    SUBROUTINE allocate_matrix_real_length8(array_in,d1_1,d1_2,d2_1,d2_2)
      IMPLICIT NONE
      Real(kind=Rkind),allocatable,    intent(inout) :: array_in(:,:)
      Integer(kind=ILkind),               intent(in) :: d1_1
      Integer(kind=ILkind),               intent(in) :: d1_2
      Integer(kind=ILkind),               intent(in) :: d2_1
      Integer(kind=ILkind),               intent(in) :: d2_2

      IF(allocated(array_in)) deallocate(array_in)

      allocate(array_in(d1_1:d1_2,d2_1:d2_2))
      array_in=0.

    ENDSUBROUTINE allocate_matrix_real_length8
    
    SUBROUTINE allocate_array_cplx_length4(array_in,d1_1,d1_2)
      IMPLICIT NONE
      Complex(kind=Rkind),allocatable, intent(inout) :: array_in(:)
      Integer(kind=Ikind),                intent(in) :: d1_1
      Integer(kind=Ikind),                intent(in) :: d1_2

      IF(allocated(array_in)) deallocate(array_in)

      allocate(array_in(d1_1:d1_2))
      array_in=0.
      
    ENDSUBROUTINE allocate_array_cplx_length4
    
    SUBROUTINE allocate_matrix_cplx_length4(array_in,d1_1,d1_2,d2_1,d2_2)
      IMPLICIT NONE
      Complex(kind=Rkind),allocatable, intent(inout) :: array_in(:,:)
      Integer(kind=Ikind),                intent(in) :: d1_1
      Integer(kind=Ikind),                intent(in) :: d1_2
      Integer(kind=Ikind),                intent(in) :: d2_1
      Integer(kind=Ikind),                intent(in) :: d2_2

      IF(allocated(array_in)) deallocate(array_in)
      
      allocate(array_in(d1_1:d1_2,d2_1:d2_2))      
      array_in=0.
      
    ENDSUBROUTINE allocate_matrix_cplx_length4

    SUBROUTINE allocate_array_cplx_length8(array_in,d1_1,d1_2)
      IMPLICIT NONE
      Complex(kind=Rkind),allocatable, intent(inout) :: array_in(:)
      Integer(kind=ILkind),               intent(in) :: d1_1
      Integer(kind=ILkind),               intent(in) :: d1_2

      IF(allocated(array_in)) deallocate(array_in)
      
      allocate(array_in(d1_1:d1_2))
      array_in=0.
      
    ENDSUBROUTINE allocate_array_cplx_length8

    SUBROUTINE allocate_array_logic(array_in,d1_1,d1_2)
      IMPLICIT NONE
      Logical,allocatable,             intent(inout) :: array_in(:)
      Integer,                            intent(in) :: d1_1
      Integer,                            intent(in) :: d1_2

      IF(allocated(array_in)) deallocate(array_in)
      
      allocate(array_in(d1_1:d1_2))
      array_in=.FALSE.
      
    ENDSUBROUTINE allocate_array_logic

    SUBROUTINE allocate_matrix_cplx_length8(array_in,d1_1,d1_2,d2_1,d2_2)
      IMPLICIT NONE
      Complex(kind=Rkind),allocatable, intent(inout) :: array_in(:,:)
      Integer*8,                          intent(in) :: d1_1
      Integer*8,                          intent(in) :: d1_2
      Integer*8,                          intent(in) :: d2_1
      Integer*8,                          intent(in) :: d2_2

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
      
#if(run_MPI)

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

      CALL MPI_Bcast(array,Int(length,kind=MPI_INTEGER_KIND),Real_MPI,                 &
                           Int(source,kind=MPI_INTEGER_KIND),MPI_COMM_WORLD,MPI_err)

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

#endif
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

#if(run_MPI)

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

      CALL MPI_Bcast(array,Int(length,kind=MPI_INTEGER_KIND),Int_MPI,                  &
                           Int(source,kind=MPI_INTEGER_KIND),MPI_COMM_WORLD,MPI_err)

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

#endif
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

#if(run_MPI)

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

      CALL MPI_Bcast(array,Int(length,kind=MPI_INTEGER_KIND),Cplx_MPI,                 &
                           Int(source,kind=MPI_INTEGER_KIND),MPI_COMM_WORLD,MPI_err)

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

#endif
    ENDSUBROUTINE MPI_Bcast_matrix_complex
!---------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------
! interface: MPI_Reduce_sum_matrix
!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_Reduce_sum_matrix_complex(matrix,d1_l,d1_u,d2_l,d2_u,destination)
      IMPLICIT NONE
      
      Complex(kind=Rkind),intent(inout)         :: matrix(:,:)
      Integer,intent(in)                        :: d1_l
      Integer,intent(in)                        :: d1_u
      Integer,intent(in)                        :: d2_l
      Integer,intent(in)                        :: d2_u
      Integer(kind=MPI_INTEGER_KIND),intent(in) :: destination
      
      Complex(kind=Rkind),allocatable           :: array(:)
      Complex(kind=Rkind),allocatable           :: array_des(:)
      Integer                                   :: length
      Integer                                   :: ii
      Integer                                   :: jj
      Integer                                   :: kk

#if(run_MPI)

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

      IF(length>huge(0_4)) Stop 'MPI_Reduce_sum_matrix_real: length > 32-bit integer'
      
      CALL MPI_Reduce(array,array_des,Int(length,kind=MPI_INTEGER_KIND),               &
                      Cplx_MPI,MPI_SUM,destination,MPI_COMM_WORLD,MPI_err)

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

#endif
    ENDSUBROUTINE MPI_Reduce_sum_matrix_complex
    
!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_Reduce_sum_matrix_real(matrix,d1_l,d1_u,d2_l,d2_u,destination)
      IMPLICIT NONE
      
      Real(kind=Rkind),intent(inout)            :: matrix(:,:)
      Integer,intent(in)                        :: d1_l
      Integer,intent(in)                        :: d1_u
      Integer,intent(in)                        :: d2_l
      Integer,intent(in)                        :: d2_u
      Integer(kind=MPI_INTEGER_KIND),intent(in) :: destination
      
      Real(kind=Rkind),allocatable              :: array(:)
      Real(kind=Rkind),allocatable              :: array_des(:)
      Integer                                   :: length
      Integer                                   :: ii
      Integer                                   :: jj
      Integer                                   :: kk

#if(run_MPI)

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

      IF(length>huge(0_4)) Stop 'MPI_Reduce_sum_matrix_real: length > 32-bit integer'
      
      CALL MPI_Reduce(array,array_des,Int(length,kind=MPI_INTEGER_KIND),               &
                      Real_MPI,MPI_SUM,destination,MPI_COMM_WORLD,MPI_err)

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

#endif
    ENDSUBROUTINE MPI_Reduce_sum_matrix_real

!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_Reduce_sum_matrix_int(matrix,d1_l,d1_u,d2_l,d2_u,destination)
      IMPLICIT NONE
      
      Integer,intent(inout)                     :: matrix(:,:)
      Integer,intent(in)                        :: d1_l
      Integer,intent(in)                        :: d1_u
      Integer,intent(in)                        :: d2_l
      Integer,intent(in)                        :: d2_u
      Integer(kind=MPI_INTEGER_KIND),intent(in) :: destination
      
      Integer,allocatable                       :: array(:)
      Integer,allocatable                       :: array_des(:)
      Integer                                   :: length
      Integer                                   :: ii
      Integer                                   :: jj
      Integer                                   :: kk

#if(run_MPI)

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

      IF(length>huge(0_4)) Stop 'MPI_Reduce_sum_matrix_real: length > 32-bit integer'

      CALL MPI_Reduce(array,array_des,Int(length,kind=MPI_INTEGER_KIND),               &
                      Int_MPI,MPI_SUM,destination,MPI_COMM_WORLD,MPI_err)
      
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

#endif
    ENDSUBROUTINE MPI_Reduce_sum_matrix_int
!---------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------
!< interface: MPI_Reduce_sum_Bcast
!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_Reduce_sum_Bcast_real(value,MS)
      IMPLICIT NONE

      Real(kind=Rkind),            intent(inout) :: value
      Integer,optional,               intent(in) :: MS
      Real(kind=Rkind)                           :: value_temp

#if(run_MPI)

      MPI_COMM_current=MPI_COMM_WORLD
      IF(present(MS)) THEN
        IF(MS==3) MPI_COMM_current=MPI_NODE_0_COMM
      ENDIF

      CALL MPI_Reduce(value,value_temp,size1_MPI,Real_MPI,MPI_SUM,root_MPI,            &
                      MPI_COMM_current,MPI_err)

      IF(MPI_id==0) value=value_temp
      CALL MPI_Bcast(value,size1_MPI,Real_MPI,root_MPI,MPI_COMM_current,MPI_err)

#endif
    ENDSUBROUTINE MPI_Reduce_sum_Bcast_real

!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_Reduce_sum_Bcast_complex(value,MS)
      IMPLICIT NONE

      Complex(kind=Rkind),         intent(inout) :: value
      Integer,optional,               intent(in) :: MS
      Complex(kind=Rkind)                        :: value_temp

#if(run_MPI)

      MPI_COMM_current=MPI_COMM_WORLD
      IF(present(MS)) THEN
        IF(MS==3) MPI_COMM_current=MPI_NODE_0_COMM
      ENDIF
      
      CALL MPI_Reduce(value,value_temp,size1_MPI,Cplx_MPI,MPI_SUM,root_MPI,            &
                      MPI_COMM_current,MPI_err)
      IF(MPI_id==0) value=value_temp
      CALL MPI_Bcast(value,size1_MPI,Cplx_MPI,root_MPI,MPI_COMM_current,MPI_err)

#endif
    ENDSUBROUTINE MPI_Reduce_sum_Bcast_complex

!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_Reduce_sum_Bcast_int(value,MS)
      IMPLICIT NONE

      Integer,                     intent(inout) :: value
      Integer,optional,               intent(in) :: MS
      Integer                                    :: value_temp

#if(run_MPI)

      MPI_COMM_current=MPI_COMM_WORLD
      IF(present(MS)) THEN
        IF(MS==3) MPI_COMM_current=MPI_NODE_0_COMM
      ENDIF

      CALL MPI_Reduce(value,value_temp,size1_MPI,Int_MPI,MPI_SUM,root_MPI,           &
                      MPI_COMM_current,MPI_err)
      IF(MPI_id==0) value=value_temp
      CALL MPI_Bcast(value,size1_MPI,Int_MPI,root_MPI,MPI_COMM_current,MPI_err)

#endif
    ENDSUBROUTINE MPI_Reduce_sum_Bcast_int

!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_Reduce_sum_Bcast_array_int(array,length,MS)
      IMPLICIT NONE

      Integer,                     intent(inout) :: array(length)
      Integer,                        intent(in) :: length 
      Integer,optional,               intent(in) :: MS
      Integer                                    :: array_temp(length)

#if(run_MPI)

      MPI_COMM_current=MPI_COMM_WORLD
      IF(present(MS)) THEN
        IF(MS==3) MPI_COMM_current=MPI_NODE_0_COMM
      ENDIF

      array_temp=0

      CALL MPI_Reduce(array,array_temp,length,Int_MPI,MPI_SUM,root_MPI,                &
                      MPI_COMM_current,MPI_err)
      IF(MPI_id==0) array=array_temp

      CALL MPI_Bcast(array,length,Int_MPI,root_MPI,MPI_COMM_current,MPI_err)

#endif
    ENDSUBROUTINE MPI_Reduce_sum_Bcast_array_int

!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_Reduce_sum_Bcast_array_real(array,length,MS)
      IMPLICIT NONE

      Real(kind=Rkind),            intent(inout) :: array(length)
      Integer,                        intent(in) :: length 
      Integer,optional,               intent(in) :: MS
      Real(kind=Rkind)                           :: array_temp(length)

#if(run_MPI)

      MPI_COMM_current=MPI_COMM_WORLD
      IF(present(MS)) THEN
        IF(MS==3) MPI_COMM_current=MPI_NODE_0_COMM
      ENDIF

      array_temp=0

      CALL MPI_Reduce(array,array_temp,length,Real_MPI,MPI_SUM,root_MPI,               &
                      MPI_COMM_current,MPI_err)
      IF(MPI_id==0) array=array_temp

      CALL MPI_Bcast(array,length,Real_MPI,root_MPI,MPI_COMM_current,MPI_err)

#endif
    ENDSUBROUTINE MPI_Reduce_sum_Bcast_array_real
!---------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------
!> interface: MPI_Reduce_max_Bcast
! combine with MPI_Reduce_sum_Bcast later 
!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_Reduce_max_Bcast_real(value)
      IMPLICIT NONE

      Real(kind=Rkind),     intent(inout) :: value
      Real(kind=Rkind)                    :: value_temp

#if(run_MPI)

      CALL MPI_Reduce(value,value_temp,size1_MPI,Real_MPI,MPI_MAX,root_MPI,            &
                      MPI_COMM_WORLD,MPI_err)
      IF(MPI_id==0) value=value_temp
      CALL MPI_Bcast(value,size1_MPI,Real_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)

#endif
    ENDSUBROUTINE MPI_Reduce_max_Bcast_real

!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_Reduce_max_Bcast_int(value)
      IMPLICIT NONE

      Integer,              intent(inout) :: value
      Integer                             :: value_temp

#if(run_MPI)

      CALL MPI_Reduce(value,value_temp,size1_MPI,Int_MPI,MPI_MAX,root_MPI,             &
                      MPI_COMM_WORLD,MPI_err)
      IF(MPI_id==0) value=value_temp
      CALL MPI_Bcast(value,size1_MPI,Int_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)

#endif
    ENDSUBROUTINE MPI_Reduce_max_Bcast_int
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

#if(run_MPI)

      length=(d1_u-d1_l+1)*(d2_u-d2_l+1)
      allocate(array(length))
      kk=0
      DO ii=d2_l,d2_u
        DO jj=d1_l,d1_u
          kk=kk+1
          array(kk)=matrix(jj,ii)
        ENDDO 
      ENDDO

      CALL MPI_Send(array,Int(length,kind=MPI_INTEGER_KIND),Real_MPI,                  &
                          Int(destination,kind=MPI_INTEGER_KIND),                      &
                          Int(tag,kind=MPI_INTEGER_KIND),MPI_COMM_WORLD,MPI_err)

      deallocate(array)

#endif
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

#if(run_MPI)

      length=(d1_u-d1_l+1)*(d2_u-d2_l+1)
      allocate(array(length))
      kk=0
      DO ii=d2_l,d2_u
        DO jj=d1_l,d1_u
          kk=kk+1
          array(kk)=matrix(jj,ii)
        ENDDO 
      ENDDO

      CALL MPI_Send(array,Int(length,kind=MPI_INTEGER_KIND),Int_MPI,                   &
                          Int(destination,kind=MPI_INTEGER_KIND),                      &
                          Int(tag,kind=MPI_INTEGER_KIND),MPI_COMM_WORLD,MPI_err)

      deallocate(array)

#endif
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

#if(run_MPI)

      length=(d1_u-d1_l+1)*(d2_u-d2_l+1)
      allocate(array(length))
      kk=0
      DO ii=d2_l,d2_u
        DO jj=d1_l,d1_u
          kk=kk+1
          array(kk)=matrix(jj,ii)
        ENDDO 
      ENDDO
      
      CALL MPI_Send(array,Int(length,kind=MPI_INTEGER_KIND),Cplx_MPI,                  &
                          Int(destination,kind=MPI_INTEGER_KIND),                      &
                          Int(tag,kind=MPI_INTEGER_KIND),MPI_COMM_WORLD,MPI_err)
      
      deallocate(array)

#endif
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

#if(run_MPI)

      length=(d1_u-d1_l+1)*(d2_u-d2_l+1)
      allocate(array(length))
      
      CALL MPI_Recv(array,Int(length,kind=MPI_INTEGER_KIND),Real_MPI,                  &
                    Int(source,kind=MPI_INTEGER_KIND),Int(tag,kind=MPI_INTEGER_KIND),  &
                    MPI_COMM_WORLD,MPI_stat,MPI_err)

      kk=0
      DO ii=d2_l,d2_u
        DO jj=d1_l,d1_u
          kk=kk+1
          matrix(jj,ii)=array(kk)
        ENDDO 
      ENDDO

      deallocate(array)

#endif
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

#if(run_MPI)

      length=(d1_u-d1_l+1)*(d2_u-d2_l+1)
      allocate(array(length))

      CALL MPI_Recv(array,Int(length,kind=MPI_INTEGER_KIND),Int_MPI,                   & 
                    Int(source,kind=MPI_INTEGER_KIND),Int(tag,kind=MPI_INTEGER_KIND),  &
                    MPI_COMM_WORLD,MPI_stat,MPI_err)

      kk=0
      DO ii=d2_l,d2_u
        DO jj=d1_l,d1_u
          kk=kk+1
          matrix(jj,ii)=array(kk)
        ENDDO 
      ENDDO

      deallocate(array)

#endif
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

#if(run_MPI)

      length=(d1_u-d1_l+1)*(d2_u-d2_l+1)
      allocate(array(length))
      
      CALL MPI_Recv(array,Int(length,kind=MPI_INTEGER_KIND),Cplx_MPI,                  &
                    Int(source,kind=MPI_INTEGER_KIND),Int(tag,kind=MPI_INTEGER_KIND),  &
                    MPI_COMM_WORLD,MPI_stat,MPI_err)

      kk=0
      DO ii=d2_l,d2_u
        DO jj=d1_l,d1_u
          kk=kk+1
          matrix(jj,ii)=array(kk)
        ENDDO 
      ENDDO

      deallocate(array)

#endif
    ENDSUBROUTINE MPI_Recv_matrix_complex

!---------------------------------------------------------------------------------------
!< interface: MPI_combine_array
! test effeiciency if use MPI_Gather
!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_combine_array_real(array,MS)
      IMPLICIT NONE
      
      Real(kind=Rkind),allocatable,intent(inout) :: array(:)
      Integer,optional,            intent(in)    :: MS

      Integer                                    :: d1
      Integer                                    :: d2
      Integer                                    :: ii

#if(run_MPI)

      IF(MPI_id/=0) THEN
        d1=bounds_MPI(1,MPI_id)
        d2=bounds_MPI(2,MPI_id)
        CALL MPI_Send(array(d1:d2),Int(d2-d1+1,kind=MPI_INTEGER_KIND),                 &
                      Real_MPI,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_err)
      ENDIF

      IF(MPI_id==0) THEN
        DO i_MPI=1,MPI_np-1
          d1=bounds_MPI(1,i_MPI)
          d2=bounds_MPI(2,i_MPI)
          CALL MPI_Recv(array(d1:d2),Int(d2-d1+1,kind=MPI_INTEGER_KIND),               &
                        Real_MPI,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_stat,MPI_err)
        ENDDO
      ENDIF

      IF(present(MS)) THEN
        d1=bounds_MPI(1,0)
        d2=bounds_MPI(2,MPI_np-1)
        IF(MS==1) THEN
          CALL MPI_Bcast(array(d1:d2),d2-d1+1,Real_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
        ELSE IF(MS==3) THEN
          IF(MPI_id==0) THEN
            DO ii=1,MPI_nodes_num-1
              CALL MPI_send(array(d1:d2),d2-d1+1,Real_MPI,MPI_nodes_p00(ii),           &
                            MPI_nodes_p00(ii),MPI_COMM_WORLD,MPI_err)
            ENDDO
          ELSEIF(MPI_nodes_p0) THEN
            CALL MPI_Recv(array(d1:d2),d2-d1+1,Real_MPI,root_MPI,MPI_id,               &
                          MPI_COMM_WORLD,MPI_stat,MPI_err)
          ENDIF
        ELSE
          ! nothing
        ENDIF
      ENDIF

#endif
    ENDSUBROUTINE MPI_combine_array_real

    !-----------------------------------------------------------------------------------
    SUBROUTINE MPI_combine_array_int4(array,MS)
      IMPLICIT NONE

      Integer(kind=Ikind),allocatable,intent(inout) :: array(:)
      Integer,optional,            intent(in)    :: MS

      Integer                                    :: d1
      Integer                                    :: d2
      Integer                                    :: ii

#if(run_MPI)

      integer(kind=Ikind) :: tempint

      IF(MPI_id/=0) THEN
        d1=bounds_MPI(1,MPI_id)
        d2=bounds_MPI(2,MPI_id)
        CALL MPI_Send(array(d1:d2),Int(d2-d1+1,kind=MPI_INTEGER_KIND),MPI_integer4,    &
                            root_MPI,MPI_id,MPI_COMM_WORLD,MPI_err)
      ENDIF

      IF(MPI_id==0) THEN
        DO i_MPI=1,MPI_np-1
          d1=bounds_MPI(1,i_MPI)
          d2=bounds_MPI(2,i_MPI)
          CALL MPI_Recv(array(d1:d2),Int(d2-d1+1,kind=MPI_INTEGER_KIND),MPI_integer4,  &
                        i_MPI,i_MPI,MPI_COMM_WORLD,MPI_stat,MPI_err)
        ENDDO
      ENDIF

      ! IF(present(bcast) .AND. bcast) THEN
      !   d1=bounds_MPI(1,0)
      !   d2=bounds_MPI(2,MPI_np-1)
      !   CALL MPI_Bcast(array(d1:d2),d2-d1+1,MPI_integer4,root_MPI,MPI_COMM_WORLD,MPI_err)
      ! ENDIF

      IF(present(MS)) THEN
        d1=bounds_MPI(1,0)
        d2=bounds_MPI(2,MPI_np-1)
        IF(MS==1) THEN
          CALL MPI_Bcast(array(d1:d2),d2-d1+1,MPI_integer4,root_MPI,                   &
                         MPI_COMM_WORLD,MPI_err)
        ELSE IF(MS==3) THEN
          IF(MPI_id==0) THEN
            DO ii=1,MPI_nodes_num-1
              CALL MPI_send(array(d1:d2),d2-d1+1,MPI_integer4,MPI_nodes_p00(ii),       &
                            MPI_nodes_p00(ii),MPI_COMM_WORLD,MPI_err)
            ENDDO
          ELSEIF(MPI_nodes_p0) THEN
            CALL MPI_Recv(array(d1:d2),d2-d1+1,MPI_integer4,root_MPI,MPI_id,           &
                          MPI_COMM_WORLD,MPI_stat,MPI_err)
          ENDIF
        ELSE
          ! nothing
        ENDIF
      ENDIF

#endif
    ENDSUBROUTINE MPI_combine_array_int4

    !-----------------------------------------------------------------------------------
    SUBROUTINE MPI_combine_array_int8(array,MS)
      IMPLICIT NONE
      
      Integer(kind=ILkind),allocatable,intent(inout) :: array(:)
      Integer,optional,            intent(in)    :: MS

      Integer                                    :: d1
      Integer                                    :: d2
      Integer                                    :: ii

#if(run_MPI)

      IF(MPI_id/=0) THEN
        d1=bounds_MPI(1,MPI_id)
        d2=bounds_MPI(2,MPI_id)
        CALL MPI_Send(array(d1:d2),Int(d2-d1+1,kind=MPI_INTEGER_KIND),                 &
                      MPI_integer8,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_err)
      ENDIF

      IF(MPI_id==0) THEN
        DO i_MPI=1,MPI_np-1
          d1=bounds_MPI(1,i_MPI)
          d2=bounds_MPI(2,i_MPI)
          CALL MPI_Recv(array(d1:d2),Int(d2-d1+1,kind=MPI_INTEGER_KIND),               &
                        MPI_integer8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_stat,MPI_err)
        ENDDO
      ENDIF

      ! IF(present(bcast) .AND. bcast) THEN
      !   d1=bounds_MPI(1,0)
      !   d2=bounds_MPI(2,MPI_np-1)
      !   CALL MPI_Bcast(array(d1:d2),d2-d1+1,MPI_integer8,root_MPI,MPI_COMM_WORLD,MPI_err)
      ! ENDIF

      IF(present(MS)) THEN
        d1=bounds_MPI(1,0)
        d2=bounds_MPI(2,MPI_np-1)
        IF(MS==1) THEN
          CALL MPI_Bcast(array(d1:d2),d2-d1+1,MPI_integer8,root_MPI,                   &
                         MPI_COMM_WORLD,MPI_err)
        ELSE IF(MS==3) THEN
          IF(MPI_id==0) THEN
            DO ii=1,MPI_nodes_num-1
              CALL MPI_send(array(d1:d2),d2-d1+1,MPI_integer8,MPI_nodes_p00(ii),       &
                            MPI_nodes_p00(ii),MPI_COMM_WORLD,MPI_err)
            ENDDO
          ELSEIF(MPI_nodes_p0) THEN
            CALL MPI_Recv(array(d1:d2),d2-d1+1,MPI_integer8,root_MPI,MPI_id,           &
                          MPI_COMM_WORLD,MPI_stat,MPI_err)
          ENDIF
        ELSE
          ! nothing
        ENDIF
      ENDIF

#endif
    ENDSUBROUTINE MPI_combine_array_int8

    !-----------------------------------------------------------------------------------
    SUBROUTINE MPI_combine_array_cplx(array,MS)
      IMPLICIT NONE

      Complex(kind=Rkind),allocatable,intent(inout)  :: array(:)
      Integer,optional,            intent(in)    :: MS

      Integer                                    :: d1
      Integer                                    :: d2
      Integer                                    :: ii

#if(run_MPI)

      IF(MPI_id/=0) THEN
        d1=bounds_MPI(1,MPI_id)
        d2=bounds_MPI(2,MPI_id)
        CALL MPI_Send(array(d1:d2),Int(d2-d1+1,kind=MPI_INTEGER_KIND),                 &
                      Cplx_MPI,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_err)
      ENDIF

      IF(MPI_id==0) THEN
        DO i_MPI=1,MPI_np-1
          d1=bounds_MPI(1,i_MPI)
          d2=bounds_MPI(2,i_MPI)
          CALL MPI_Recv(array(d1:d2),Int(d2-d1+1,kind=MPI_INTEGER_KIND),               &
                        Cplx_MPI,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_stat,MPI_err)
        ENDDO
      ENDIF

      ! IF(present(bcast) .AND. bcast) THEN
      !   d1=bounds_MPI(1,0)
      !   d2=bounds_MPI(2,MPI_np-1)
      !   CALL MPI_Bcast(array(d1:d2),d2-d1+1,Cplx_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
      ! ENDIF

      IF(present(MS)) THEN
        d1=bounds_MPI(1,0)
        d2=bounds_MPI(2,MPI_np-1)
        IF(MS==1) THEN
          CALL MPI_Bcast(array(d1:d2),d2-d1+1,Cplx_MPI,root_MPI,                       &
                         MPI_COMM_WORLD,MPI_err)
        ELSE IF(MS==3) THEN
          IF(MPI_id==0) THEN
            DO ii=1,MPI_nodes_num-1
              CALL MPI_send(array(d1:d2),d2-d1+1,Cplx_MPI,MPI_nodes_p00(ii),           &
                            MPI_nodes_p00(ii),MPI_COMM_WORLD,MPI_err)
            ENDDO
          ELSEIF(MPI_nodes_p0) THEN
            CALL MPI_Recv(array(d1:d2),d2-d1+1,Cplx_MPI,root_MPI,MPI_id,               &
                          MPI_COMM_WORLD,MPI_stat,MPI_err)
          ENDIF
        ELSE
          ! nothing
        ENDIF
      ENDIF

#endif
    ENDSUBROUTINE MPI_combine_array_cplx

    !-----------------------------------------------------------------------------------
    ! dependents on the length array comtain the length on each threads 
    ! consider MPI_Gether
    ! combine with previous ones later
    ! add for the other scheme later
    SUBROUTINE MPI_combine_array_general_cplx(array_all,array,lengths,MS)
      IMPLICIT NONE

      Complex(kind=Rkind),         intent(inout) :: array_all(:)
      Complex(kind=Rkind),         intent(inout) :: array(:)
      Integer,                        intent(in) :: lengths(0:MPI_np-1)
      Integer,optional,               intent(in) :: MS

      Integer                                    :: d1
      Integer                                    :: d2

#if(run_MPI)

      IF(MPI_id/=0) THEN
        IF(present(MS)) THEN
          IF(MS/=3 .OR. (MS==3 .AND. MPI_nodes_p0)) THEN
            d1=1
            d2=lengths(MPI_id)
            CALL MPI_Send(array(d1:d2),Int(d2-d1+1,kind=MPI_INTEGER_KIND),             &
                          Cplx_MPI,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_err)
          ENDIF
        ELSE
          d1=1
          d2=lengths(MPI_id)
          CALL MPI_Send(array(d1:d2),Int(d2-d1+1,kind=MPI_INTEGER_KIND),               &
                        Cplx_MPI,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_err)
        ENDIF
      ENDIF

      IF(MPI_id==0) THEN
        array_all(1:lengths(0))=array(1:lengths(0))

        d1=lengths(0)+1
        DO i_MPI=1,MPI_np-1
          IF(present(MS)) THEN
            IF(MS/=3 .OR. (MS==3 .AND. MPI_nodes_00(i_MPI))) THEN
              d2=d1+lengths(i_MPI)-1
              CALL MPI_Recv(array_all(d1:d2),Int(d2-d1+1,kind=MPI_INTEGER_KIND),       &
                            Cplx_MPI,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_stat,MPI_err)
              d1=d2+1
            ENDIF
          ELSE
            d2=d1+lengths(i_MPI)-1
            CALL MPI_Recv(array_all(d1:d2),Int(d2-d1+1,kind=MPI_INTEGER_KIND),         &
                          Cplx_MPI,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_stat,MPI_err)
            d1=d2+1
          ENDIF
        ENDDO
      ENDIF

#endif
    ENDSUBROUTINE MPI_combine_array_general_cplx
!---------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------
!< interface: MPI_collect_info_int
!---------------------------------------------------------------------------------------
    SUBROUTINE MPI_collect_info_int(array,bcast,MS)
      IMPLICIT NONE

      Integer,                     intent(inout) :: array(0:MPI_np-1)
      Logical,optional,               intent(in) :: bcast
      Integer,optional,               intent(in) :: MS

#if(run_MPI)

      IF(MPI_id/=0) THEN
        IF(present(MS)) THEN
          IF(MS/=3 .OR. (MS==3 .AND. MPI_nodes_p0)) THEN
            CALL MPI_Send(array(MPI_id),size1_MPI,Int_MPI,root_MPI,MPI_id,             &
                          MPI_COMM_WORLD,MPI_err)
          ENDIF
        ELSE
          CALL MPI_Send(array(MPI_id),size1_MPI,Int_MPI,root_MPI,MPI_id,               &
                        MPI_COMM_WORLD,MPI_err)
        ENDIF
      ENDIF

      IF(MPI_id==0) THEN
        DO i_MPI=1,MPI_np-1
          IF(present(MS)) THEN
            IF(MS/=3 .OR. (MS==3 .AND. MPI_nodes_00(i_MPI))) THEN
              CALL MPI_Recv(array(i_MPI),size1_MPI,Int_MPI,i_MPI,i_MPI,                &
                            MPI_COMM_WORLD,MPI_stat,MPI_err)
            ENDIF
          ELSE
            CALL MPI_Recv(array(i_MPI),size1_MPI,Int_MPI,i_MPI,i_MPI,                  &
                          MPI_COMM_WORLD,MPI_stat,MPI_err)
          ENDIF
        ENDDO
      ENDIF

      MPI_COMM_current=MPI_COMM_WORLD
      IF(present(MS)) THEN
        IF(MS==3) MPI_COMM_current=MPI_NODE_0_COMM
      ENDIF

      IF(present(bcast) .AND. bcast) THEN
        CALL MPI_Bcast(array(0:MPI_np-1),MPI_np,Int_MPI,root_MPI,MPI_COMM_current,MPI_err)
      ENDIF

#endif
    ENDSUBROUTINE MPI_collect_info_int
    
    !-----------------------------------------------------------------------------------
    SUBROUTINE MPI_collect_info_real(array,bcast,MS)
      IMPLICIT NONE

      Real(kind=Rkind),            intent(inout) :: array(0:MPI_np-1)
      Logical,optional,               intent(in) :: bcast
      Integer,optional,               intent(in) :: MS

#if(run_MPI)

      IF(MPI_id/=0) THEN
        IF(present(MS)) THEN
          IF(MS/=3 .OR. (MS==3 .AND. MPI_nodes_p0)) THEN
            CALL MPI_Send(array(MPI_id),size1_MPI,Real_MPI,root_MPI,MPI_id,            &
                          MPI_COMM_WORLD,MPI_err)
          ENDIF
        ELSE
          CALL MPI_Send(array(MPI_id),size1_MPI,Real_MPI,root_MPI,MPI_id,              &
                        MPI_COMM_WORLD,MPI_err)
        ENDIF
      ENDIF

      IF(MPI_id==0) THEN
        DO i_MPI=1,MPI_np-1
          IF(present(MS)) THEN
            IF(MS/=3 .OR. (MS==3 .AND. MPI_nodes_00(i_MPI))) THEN
              CALL MPI_Recv(array(i_MPI),size1_MPI,Real_MPI,i_MPI,i_MPI,               &
                            MPI_COMM_WORLD,MPI_stat,MPI_err)
            ENDIF
          ELSE
            CALL MPI_Recv(array(i_MPI),size1_MPI,Real_MPI,i_MPI,i_MPI,                 &
                          MPI_COMM_WORLD,MPI_stat,MPI_err)
          ENDIF
        ENDDO
      ENDIF

      MPI_COMM_current=MPI_COMM_WORLD
      IF(present(MS)) THEN
        IF(MS==3) MPI_COMM_current=MPI_NODE_0_COMM
      ENDIF

      IF(present(bcast) .AND. bcast) THEN
        CALL MPI_Bcast(array(0:MPI_np-1),MPI_np,Real_MPI,root_MPI,MPI_COMM_current,MPI_err)
      ENDIF

#endif
    ENDSUBROUTINE MPI_collect_info_real
END MODULE mod_MPI_aux


