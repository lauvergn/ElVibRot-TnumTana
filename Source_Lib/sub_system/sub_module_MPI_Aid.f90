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
  USE mod_MPI
  IMPLICIT NONE
  
  INTEGER                       :: memory_RSS   !< memory usage
  
  TYPE multi_array
    Integer,allocatable :: array(:)
  END TYPE multi_array
  
  TYPE multi_array4
    Integer*4,allocatable :: array(:)
  END TYPE multi_array4  
  
  INTERFACE allocate_array  
    module procedure allocate_array_int4_length4
    module procedure allocate_array_int4_length8 
    module procedure allocate_array_int8_length4
    module procedure allocate_array_int8_length8
    module procedure allocate_array_real_length4
    module procedure allocate_array_real_length8  
  END INTERFACE
!---------------------------------------------------------------------------------------
  Contains
    !> check total memory used at certain point
    SUBROUTINE system_mem_usage(memory_RSS,name)
      USE mod_NumParameters
      ! USE ifport ! if on intel compiler
      IMPLICIT NONE
      Integer, intent(out) :: memory_RSS
      Character(len=200):: filename=' '
      Character(len=80) :: line
      Character(len=8)  :: pid_char=' '
      Integer :: pid
      Logical :: ifxst
      Character (len=*), intent(in) :: name

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
        write(out_unitp,121) name,memory_RSS,MPI_id
121     format('memory check at ',a,': ',i4,' from ',i4)
      ENDIF
    ENDSUBROUTINE system_mem_usage

!---------------------------------------------------------------------------------------
!> write for MPI outpout
!---------------------------------------------------------------------------------------    
    SUBROUTINE MPI_write(out_channel,out_message)
      USE mod_system
      IMPLICIT NONE
      integer :: out_channel,MPIid
      Character (len=*), intent(in) :: out_message
      
      write(out_channel,*) out_message,' from ',MPI_id
    ENDSUBROUTINE
    
    SUBROUTINE MPI0_write(out_channel,out_message)
      USE mod_system
      IMPLICIT NONE
      integer :: out_channel
      Character (len=*), intent(in) :: out_message
      
      IF(MPI_id==0) write(out_channel,*) out_message
    ENDSUBROUTINE
    
    SUBROUTINE MPI0_write_line(out_channel)
      USE mod_system
      IMPLICIT NONE
      integer :: out_channel,MPIid
            
      IF(MPI_id==0) write(out_channel,*) '------------------------------------------------------------'
    ENDSUBROUTINE
    
    SUBROUTINE MPI0_write_dline(out_channel)
      USE mod_system
      IMPLICIT NONE
      integer :: out_channel,MPIid
      
      IF(MPI_id==0) write(out_channel,*) '============================================================'
    ENDSUBROUTINE
!---------------------------------------------------------------------------------------

!---------------------------------------------------------------------------------------
!> interface: allocate_array    
!---------------------------------------------------------------------------------------
    SUBROUTINE allocate_array_int4_length4(array_in,length)
      USE mod_system
      IMPLICIT NONE
      Integer*4,allocatable, intent(inout) :: array_in(:)
      Integer*4, intent(in)                :: length
      
      IF(allocated(array_in)) deallocate(array_in)
      allocate(array_in(length))
      array_in=0
    END SUBROUTINE
    
    SUBROUTINE allocate_array_int4_length8(array_in,length)
      USE mod_system
      IMPLICIT NONE
      Integer*4,allocatable, intent(inout) :: array_in(:)
      Integer*8, intent(in)                  :: length
      
      IF(allocated(array_in)) deallocate(array_in)
      allocate(array_in(length))
      array_in=0
    END SUBROUTINE
    
    SUBROUTINE allocate_array_int8_length4(array_in,length)
      USE mod_system
      IMPLICIT NONE
      Integer*8,allocatable, intent(inout) :: array_in(:)
      Integer*4, intent(in)                  :: length
      
      IF(allocated(array_in)) deallocate(array_in)
      allocate(array_in(length))
      array_in=0
    END SUBROUTINE
    
    SUBROUTINE allocate_array_int8_length8(array_in,length)
      USE mod_system
      IMPLICIT NONE
      Integer*8,allocatable, intent(inout) :: array_in(:)
      Integer*8, intent(in)                  :: length
      
      IF(allocated(array_in)) deallocate(array_in)
      allocate(array_in(length))
      array_in=0
    END SUBROUTINE
    
    SUBROUTINE allocate_array_real_length4(array_in,length)
      USE mod_system
      IMPLICIT NONE
      Real(kind=Rkind),allocatable, intent(inout) :: array_in(:)
      Integer*4, intent(in)                         :: length

      IF(allocated(array_in)) deallocate(array_in)
      allocate(array_in(length))
      array_in=0.
    END SUBROUTINE

    SUBROUTINE allocate_array_real_length8(array_in,length)
      USE mod_system
      IMPLICIT NONE
      Real(kind=Rkind),allocatable, intent(inout) :: array_in(:)
      Integer*8, intent(in)                       :: length

      IF(allocated(array_in)) deallocate(array_in)
      allocate(array_in(length))
      array_in=0.
    END SUBROUTINE
!---------------------------------------------------------------------------------------

END MODULE mod_MPI_Aid


