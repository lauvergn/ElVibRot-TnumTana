!=======================================================================================
!> module for implementation MPI
!> to use long integer, use 64-bit compiled MPI
!> include:
!>  - basic MPI variables
!>  - asisstant MPI variables
!>  - varilables for recording time usage 
!>  - varilables for process control (to be added to certain type later)
!>  - 
!=======================================================================================
MODULE mod_MPI
#if(run_MPI)
  USE MPI_F08
  IMPLICIT NONE
  
  !-------------------------------------------------------------------------------------
  !> use save to keep MPI id cross the project.
  !> other options:
  !>  -call MPI_Comm_rank & MPI_Comm_size in each subroutine
  !>  -pass as variables. 
  !>  -use as commmon varilables. (not recomm.) 
  !
  !> variables for MPI
  !>  - MPI_err: error flag for MPI
  !>  - MPI_id: process ID, 0~MPI_np-1  
  !>  - MPI_np: total MPI
  !
  !> @note 'offical' MPI variables are started with MPI_...
  !!  while the other assistant variables are end with ..._MPI
  !-------------------------------------------------------------------------------------
  Integer(kind=MPI_INTEGER_KIND) :: MPI_err      !< error flag for MPI
  Integer(kind=MPI_INTEGER_KIND) :: MPI_id=0     !< rocess ID, 0~MPI_np-1
  Integer(kind=MPI_INTEGER_KIND) :: MPI_np       !< total number of MPI threads
  
  !> note the difference on MPI_stat in "USE MPI" and "USE MPI_F08"
  ! Integer(kind=MPI_INTEGER_KIND) :: MPI_stat(MPI_STATUS_SIZE) !< status of MPI process
  TYPE(MPI_Status)               :: MPI_stat          !< status of MPI process
  TYPE(MPI_Datatype)             :: MPI_int           !< integer type of default MPI
  TYPE(MPI_Datatype)             :: MPI_int_fortran   !< integer type of default fortran
  TYPE(MPI_Datatype)             :: MPI_real_fortran  !< real type for fortran (Rkind)
  TYPE(MPI_Datatype)             :: MPI_clpx_fortran  !< complex type for fortran (Rkind)
  Integer(kind=MPI_INTEGER_KIND) :: MPI_rec_source    !< for MPI_RECV, not used
  Integer(kind=MPI_INTEGER_KIND) :: MPI_tag1     !< tag for MPI send and receive
  Integer(kind=MPI_INTEGER_KIND) :: MPI_tag2     !< tag for MPI send and receive
  Integer                        :: integer_MPI  !< get integer type for MPI

  Integer(kind=MPI_INTEGER_KIND) :: root_MPI=0   !< used for MPI functions, master ID 
  Integer(kind=MPI_INTEGER_KIND) :: size1_MPI=1  !< used for MPI functions, size 1
  Integer(kind=MPI_INTEGER_KIND) :: i_MPI        !< fake MPI thread id
  Integer                        :: bound1_MPI   !< up boundary of works for each thread
  Integer                        :: bound2_MPI   !< dn boundary of works for each thread
  Integer                        :: iG1_MPI      !< up boundary of iGs in action 
  Integer                        :: iG2_MPI      !< dn boundary of iGs in action 
  Integer                        :: nb_per_MPI   !< number of distribed works per thread
  Integer                        :: nb_rem_MPI   !< remainder of distribed works
  Integer,allocatable            :: iGs_MPI(:,:) !< iG boundary of each theard in action
  
  Logical                        :: Srep_MPI=.FALSE. !< MPI working on full Smolyak Rep. 
  Integer                        :: action_MPI=2     !< MPI working on SR scheme #. 


  !Common /group_MPI_world/     MPI_err, MPI_id, MPI_np, MPI_status
  !Common /group_MPI_tag/       MPI_tag0,MPI_tag1,MPI_tag2,MPI_tag3
  !-------------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------------
  !> varilables for recording time usage
  !Real*8                        :: time_MPI_action !< total time used in action 
  !Real*8                        :: time_comm       !< total time used in comm. in action
  !Real*8                        :: time_point1     !< tags for save current time
  !Real*8                        :: time_point2
  !Real*8                        :: time_temp1
  !Real*8                        :: time_temp2
  Integer                        :: time_MPI_action !< total time used in action 
  Integer                        :: time_comm       !< total time used in comm. in action
  Integer                        :: time_point1     !< tags for save current time
  Integer                        :: time_point2
  Integer                        :: time_temp1
  Integer                        :: time_temp2
  Integer                        :: time_rate       !< for function system_clock()
  Integer                        :: time_max        !< for function system_clock()
  ! time_rate in kind=4: COUNT in system_clock represents milliseconds
  ! time_rate in kind>4: COUNT in system_clock represents micro- or nanoseconds
  
  !-------------------------------------------------------------------------------------
  ! varilables for process control, add to certain type later
  Logical                        :: if_propa
  Logical                        :: Grid_allco
  
  !Common /group_MPI_time/     time_MPI_action,time_point1,time_point2,time_action
  !Common /group_MPI_control/  if_propa,Grid_allco

!=======================================================================================
  Contains  
  !-------------------------------------------------------------------------------------
  !> MPI initialization 
  SUBROUTINE MPI_initialization(Rkind)
    IMPLICIT NONE
    Integer,intent(in)    :: Rkind

    CALL MPI_Init(MPI_err)
    CALL MPI_Comm_rank(MPI_COMM_WORLD, MPI_id, MPI_err)
    CALL MPI_Comm_size(MPI_COMM_WORLD, MPI_np, MPI_err)
    
    time_MPI_action=0
    time_comm=0
    Grid_allco=.True.
    
    !> get Fortran default bit
    IF(sizeof(integer_MPI)==4) THEN
      MPI_int_fortran=MPI_integer4
    ELSEIF(sizeof(integer_MPI)==8) THEN
      MPI_int_fortran=MPI_integer8
    ELSE
      STOP 'integer neither 4 or 8 in default fortran'
    ENDIF
    
    !> define Fortran Real type according to Rkind
    IF(Rkind==8) THEN
      MPI_real_fortran=MPI_Real8
      MPI_clpx_fortran=MPI_Complex8
    ELSEIF(Rkind==16) THEN
      MPI_real_fortran=MPI_Real16
      MPI_clpx_fortran=MPI_Complex16
    ELSE
      STOP 'Rkind neither 64 or 128 bit, define MPI_real_fortran in sub_module_MPI.f90'
    ENDIF
    
    !> get MPI default bit
    IF(MPI_INTEGER_KIND==4) THEN
      MPI_int=MPI_integer4
    ELSEIF(MPI_INTEGER_KIND==8) THEN
      MPI_int=MPI_integer8
    ELSE
      STOP 'integer neither 4 or 8 in default MPI'
    ENDIF
    
  END SUBROUTINE MPI_initialization
  
#else
  Integer                        :: MPI_id=0     !< fake MPI_id, for convenience
#endif

END MODULE mod_MPI


