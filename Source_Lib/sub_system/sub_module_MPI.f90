!=======================================================================================
!> @brief 
!> module for implementating MPI
!> to use long integer, use 64-bit compiled MPI
!>
!> include:
!> basic MPI variables,
!> asisstant MPI variables,
!> varilables for recording time usage,
!> varilables for process control (to be added to certain type later),
!>  
!> @note there is a bug in openMPI for working on 64-bit fortran USE MPI_F08
!>        which will possiblely be solved in future release MPI_v4.x or MPI_v5.x.
!>        The current strategy apply USE MPI or 
!>        keep using 32-bit with extra subroutine for large array, need further test
!
!> @note further implemention: using 32-bit MPI in 64-bit Fortran
!=======================================================================================
MODULE mod_MPI
#if(run_MPI)
  USE MPI
  !USE MPI_F08
  IMPLICIT NONE
  
  !-------------------------------------------------------------------------------------
  !> @note 'offical' MPI variables are started with MPI_...
  !!        while the other assistant variables are end with ..._MPI
  !-------------------------------------------------------------------------------------

  !> note the difference in "USE MPI" and "USE MPI_F08"
  !-USE MPI-----------------------------------------------------------------------------
  Integer(kind=MPI_INTEGER_KIND) :: MPI_stat(MPI_STATUS_SIZE) !< status of MPI process
  Integer(kind=MPI_INTEGER_KIND) :: Int_MPI         !< integer type of default MPI
  Integer(kind=MPI_INTEGER_KIND) :: Real_MPI        !< real type for fortran (Rkind)
  Integer(kind=MPI_INTEGER_KIND) :: Cplx_MPI        !< complex type for fortran (Rkind)
  Integer(kind=MPI_INTEGER_KIND) :: Int_fortran     !< integer type of default fortran
  !-USE MPI_F08-------------------------------------------------------------------------
!  TYPE(MPI_Status)               :: MPI_stat        !< status of MPI process
!  TYPE(MPI_Datatype)             :: Int_MPI         !< integer type of default MPI
!  TYPE(MPI_Datatype)             :: Real_MPI        !< real type for fortran (Rkind)
!  TYPE(MPI_Datatype)             :: Cplx_MPI        !< complex type for fortran (Rkind)
!  TYPE(MPI_Datatype)             :: Int_fortran     !< integer type of default fortran
  !-------------------------------------------------------------------------------------
  Integer(kind=MPI_INTEGER_KIND) :: MPI_err         !< error flag for MPI
  Integer(kind=MPI_INTEGER_KIND) :: MPI_id          !< rocess ID, 0~MPI_np-1
  Integer(kind=MPI_INTEGER_KIND) :: MPI_np          !< total number of MPI threads
  Integer(kind=MPI_INTEGER_KIND) :: MPI_tag1        !< tag for MPI send and receive
  Integer(kind=MPI_INTEGER_KIND) :: MPI_tag2        !< tag for MPI send and receive
  
  Integer(kind=MPI_INTEGER_KIND) :: source_MPI      !< for MPI_RECV, not used
  Integer                        :: integer_MPI     !< get integer type for MPI

  Integer(kind=MPI_INTEGER_KIND) :: root_MPI=0      !< used for MPI functions, master ID 
  Integer(kind=MPI_INTEGER_KIND) :: size1_MPI=1     !< used for MPI functions, size 1
  Integer                        :: ONE_1=1
  Integer(kind=MPI_INTEGER_KIND) :: i_MPI           !< fake MPI thread id
  Integer                        :: bound1_MPI      !< up boundary of works for threads
  Integer                        :: bound2_MPI      !< dn boundary of works for threads
  Integer                        :: nb_per_MPI      !< number of works per thread
  Integer                        :: nb_rem_MPI      !< remainder of distribed works
  Integer,allocatable            :: iGs_MPI(:,:)    !< iG boundary of theards in action
  Integer,allocatable            :: iGs_MPI0(:,:)
  Integer,allocatable            :: iGs_MPI_mc(:,:,:) !< iG boundary of theards in action for memory control
  Integer,allocatable            :: iGs_MPI_mc0(:,:,:)
  Integer,allocatable            :: bounds_MPI(:,:) !< boundary of each theards
  Integer                        :: Num_L2=6        !<1/Num_L2 threads used for level2 distributor
  Integer                        :: n_level2        !< number of level2 distributor 
  Logical                        :: keep_MPI=.FALSE.!< if keep valables on current thread 
  Logical                        :: Srep_MPI=.FALSE.!< MPI working on full Smolyak Rep. 
  Integer                        :: MPI_scheme=0    !< MPI working on MPI scheme #. 
                                                    !! 0 auto
  Integer                        :: MPI_nb_WP=1     ! nb of levels to converge
  Integer                        :: MPI_mc=1        !< memory control for acrion
                                                    !< 1/MPI_mc of smolyak terms are 
                                                    !< processed each time, 
                                                    !< for reducing memory
  Logical                        :: iGs_auto=.True. !< if auto-adjust the distribution 
                                                    !< of Smolyak terms
  !-------------------------------------------------------------------------------------
  !> varilables for recording time usage
  Integer                        :: time_MPI_action=0  !< total time used in action 
  Integer                        :: time_comm=0        !< total time used in comm. in action
  Integer,allocatable            :: time_MPI_calcu(:)  !< time used in each action for sharable claculation
  Integer,allocatable            :: time_MPI_local(:)  !< local time used in each action
  Integer,allocatable            :: time_MPI_act_all(:)!< time used in each action
  Integer                        :: time_point1        !< tags for save current time
  Integer                        :: time_point2
  Integer                        :: time_temp1
  Integer                        :: time_temp2
  Integer                        :: time_rate          !< for function system_clock()
  Integer                        :: time_max           !< for function system_clock()
  ! time_rate in kind=4: COUNT in system_clock represents milliseconds
  ! time_rate in kind>4: COUNT in system_clock represents micro- or nanoseconds
  
  !-------------------------------------------------------------------------------------
  ! varilables for process control, add to certain type later
!  Logical                        :: Grid_allco=.True.

!=======================================================================================
  Contains  

  !-------------------------------------------------------------------------------------
  !> MPI initialization 
  SUBROUTINE ini_MPI()
    USE mod_NumParameters
    IMPLICIT NONE

    CALL MPI_Init(MPI_err)
    CALL MPI_Comm_rank(MPI_COMM_WORLD,MPI_id,MPI_err)
    CALL MPI_Comm_size(MPI_COMM_WORLD,MPI_np,MPI_err)

    !> get Fortran default bit
    IF(sizeof(integer_MPI)==4) THEN
      Int_fortran=MPI_integer4
    ELSEIF(sizeof(integer_MPI)==8) THEN
      Int_fortran=MPI_integer8
    ELSE
      STOP 'integer neither 4 or 8 in default fortran'
    ENDIF
    
    !> define Fortran Real type according to Rkind
    IF(Rkind==8) THEN
      Real_MPI=MPI_Real8
      Cplx_MPI=MPI_Complex8
    ELSEIF(Rkind==16) THEN
      Real_MPI=MPI_Real16
      Cplx_MPI=MPI_Complex16
    ELSE
      STOP 'Rkind neither 64 or 128 bit, define Real_MPI in sub_module_MPI.f90'
    ENDIF
    
    !> get MPI default bit
    IF(MPI_INTEGER_KIND==4) THEN
      Int_MPI=MPI_integer4
    ELSEIF(MPI_INTEGER_KIND==8) THEN
      Int_MPI=MPI_integer8
    ELSE
      STOP 'integer neither 4 or 8 in default MPI'
    ENDIF

    write(out_unitp,*) 'Initiaize MPI with ', MPI_np, 'cores.'
    write(out_unitp,*) 'NOTE: MPI in progress. If get memory error, check if           &
                                      the variables are just allocated on root threads.'
    write(out_unitp,*)
    write(out_unitp,*) 'Integer type of default Fortran Compiler:',                    &
                                         sizeof(integer_MPI),', MPI: ',MPI_INTEGER_KIND
    IF(sizeof(integer_MPI)/=MPI_INTEGER_KIND)                                          &
                                 STOP 'Please use same integer type for Fortran and MPI'

  END SUBROUTINE ini_MPI

  !-------------------------------------------------------------------------------------
  !> MPI_Finalize
  SUBROUTINE end_MPI()
    USE mod_NumParameters
    IMPLICIT NONE

    IF(MPI_id==0) THEN
      write(out_unitp,*) 'time used for action: ',                                     &
             real(time_MPI_action,kind=Rkind)/real(time_rate,kind=Rkind),' from ',MPI_id
      write(out_unitp,*) 'time used for MPI communication: ',                          &
                  real(time_comm,kind=Rkind)/real(time_rate,kind=Rkind),' from ', MPI_id
    ENDIF

    CALL MPI_Finalize(MPI_err)
  END SUBROUTINE end_MPI

!=======================================================================================

#else
  !< fake variables, for convenience
  Integer                        :: MPI_err=0
  Integer                        :: MPI_id=0
  Integer                        :: MPI_np=1
  
  Integer,parameter              :: MPI_INTEGER_KIND=4
  Integer,parameter              :: MPI_ADDRESS_KIND=4
  Integer                        :: MPI_scheme=0
  Integer                        :: MPI_nb_WP=2
  Integer                        :: iG1_MPI      
  Integer                        :: iG2_MPI  
  Integer(kind=MPI_INTEGER_KIND) :: root_MPI=0      
  Integer(kind=MPI_INTEGER_KIND) :: size1_MPI=1
  Integer(kind=MPI_INTEGER_KIND) :: i_MPI    
  Logical                        :: Srep_MPI=.FALSE.  
  Integer,allocatable            :: iGs_MPI(:,:)  
  Integer,allocatable            :: bounds_MPI(:,:) 
  Integer                        :: Num_L2
  Integer                        :: n_level2
  Integer                        :: MPI_mc=1
  Logical                        :: iGs_auto=.True.
  Logical                        :: keep_MPI=.True.

  Integer                        :: time_rate       !< for function system_clock()
  Integer                        :: time_max        !< for function system_clock()

  Contains
  SUBROUTINE ini_MPI()
  ENDSUBROUTINE ini_MPI

  SUBROUTINE end_MPI()
  ENDSUBROUTINE end_MPI
#endif

END MODULE mod_MPI


