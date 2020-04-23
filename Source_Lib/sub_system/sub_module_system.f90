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
MODULE mod_system
      USE mod_NumParameters
      USE mod_string
      USE mod_file
      USE mod_RW_MatVec
      USE mod_Frac
      USE mod_memory
      USE mod_memory_Pointer
      USE mod_memory_NotPointer
      USE mod_MPI
      IMPLICIT NONE

      INTERFACE compare_tab
        MODULE PROCEDURE compare_la, compare_tab_int, compare_tab_real, &
                         compare_tab_cmplx
      END INTERFACE
      PRIVATE :: compare_la, compare_tab_int, compare_tab_real, compare_tab_cmplx

      INTERFACE inferior_tab
        MODULE PROCEDURE inferior_tab_real, inferior_tab_int
      END INTERFACE
      PRIVATE :: inferior_tab_real, inferior_tab_int

#if defined(__TNUM_VER)
      character (len=Name_len) :: Tnum_version = __TNUM_VER
#else
      character (len=Name_len) :: Tnum_version = "unknown: -D__TNUM_VER=?"
#endif

#if defined(__TANA_VER)
      character (len=Name_len) :: Tana_version = __TANA_VER
#else
      character (len=Name_len) :: Tana_version = "unknown: -D__TANA_VER=?"
#endif

#if defined(__EVR_VER)
      character (len=Name_len) :: EVR_version = __EVR_VER
#else
      character (len=Name_len) :: EVR_version = "unknown: -D__EVR_VER=?"
#endif

#if defined(__EVRTPATH)
      character (len=Line_len) :: EVRT_path   =                         &
       __EVRTPATH
#else
      character (len=Line_len) :: EVRT_path   = '~/ElVibRot'
#endif

#if defined(__COMPILE_DATE)
      character (len=Line_len) :: compile_date = __COMPILE_DATE
#else
      character (len=Line_len) :: compile_date = "unknown: -D__COMPILE_DATE=?"
#endif

#if defined(__COMPILE_HOST)
      character (len=Line_len) :: compile_host = __COMPILE_HOST
#else
      character (len=Line_len) :: compile_host = "unknown: -D__COMPILE_HOST=?"
#endif
#if defined(__COMPILER)
      character (len=Line_len) :: compiler = __COMPILER
#else
      character (len=Line_len) :: compiler = "unknown: -D__COMPILER=?"
#endif
#if defined(__COMPILER_VER)
      character (len=Line_len) :: compiler_ver = __COMPILER_VER
#else
      character (len=Line_len) :: compiler_ver = "unknown: -D__COMPILER_VER=?"
#endif
#if defined(__COMPILER_OPT)
      character (len=Line_len) :: compiler_opt = &
      __COMPILER_OPT
#else
      character (len=Line_len) :: compiler_opt = "unknown: -D__COMPILER_OPT=?"
#endif
#if defined(__COMPILER_LIBS)
      character (len=Line_len) :: compiler_libs = __COMPILER_LIBS
#else
      character (len=Line_len) :: compiler_libs = "unknown: -D__COMPILER_LIBS=?"
#endif


      logical :: openmp = .FALSE.
      logical :: openmpi= .FALSE. 
      integer :: MatOp_omp,OpPsi_omp,BasisTOGrid_omp,Grid_omp,SG4_omp
      integer :: MatOp_maxth,OpPsi_maxth,BasisTOGrid_maxth,Grid_maxth,SG4_maxth

      integer (kind=ILkind) :: nb_mult_BTOG  = 0
      integer (kind=ILkind) :: nb_mult_GTOB  = 0
      integer (kind=ILkind) :: nb_mult_OpPsi = 0

      integer, parameter :: max_HADA = 5000
      integer, parameter :: max_nb_G_FOR_print = 2000
      !integer, parameter :: max_nb_G_FOR_print = 20000

      integer :: SGtype               = -1
      integer :: FilePsiVersion       = 0
      logical :: NewBasisEl           = .FALSE.
      logical :: print_CoordType_done = .FALSE.! if T, the CoordType has been already print


      TYPE param_FOR_optimization

        integer                        :: nb_OptParam    = 0
        integer                        :: i_OptParam     = 0

        real (kind=Rkind), allocatable :: Val_RVec(:)
        integer, allocatable           :: opt_RVec(:)

        character (len=Name_len) :: Optimization_param  = 'geometry'

      END TYPE param_FOR_optimization

      TYPE param_EVRT_calc
        integer :: optimization     = 0
        logical :: EVR              = .TRUE.   ! ElVibRot (default)
        logical :: analysis_only    = .FALSE.
        logical :: intensity_only   = .FALSE.
        logical :: cart             = .FALSE.
        logical :: GridTOBasis_test = .FALSE.
        logical :: OpPsi_test       = .FALSE.
        logical :: nDfit            = .FALSE.
        logical :: nDGrid           = .FALSE.
        logical :: main_test        = .FALSE.

      END TYPE param_EVRT_calc

      TYPE param_time
        integer :: count_old,count_ini
        real    :: t_cpu_old,t_cpu_ini
        logical :: begin = .TRUE.
      END TYPE param_time

      TYPE (param_FOR_optimization), save :: para_FOR_optimization
      TYPE (param_EVRT_calc),        save :: para_EVRT_calc

      CONTAINS

!===============================================================================
      !!@description: TODO
      !!@param: TODO
      SUBROUTINE time_perso(name)
      USE mod_MPI
      IMPLICIT NONE

        character (len=*), intent(in) :: name

        !local variables
        integer           :: tab_time(8) = 0
        real (kind=Rkind) :: dt_real,t_real
        real              :: dt_cpu,t_cpu
        integer           :: seconds,minutes,hours,days


        CALL date_and_time(values=tab_time)
        IF(MPI_id==0) write(out_unitp,21) name,tab_time(5:8),tab_time(3:1:-1)
 21     format('     Time and date in ',a,' : ',i2,'h:',                &
               i2,'m:',i2,'.',i3,'s, the ',i2,'/',i2,'/',i4)

        CALL DeltaTime(dt_real,t_real,dt_cpu,t_cpu)

        !============================================
        !real and cpu delta times in the subroutine: "name"
        seconds = int(dt_real)
        minutes = seconds/60
        seconds = mod(seconds,60)
        hours   = minutes/60
        minutes = mod(minutes,60)
        days    = hours/24
        hours   = mod(hours,24)

#if(run_MPI)
        write(out_unitp,31) dt_real,name,MPI_id
31      format('        real (s): ',f18.3,' in ',a, ' from MPI id ',i4)
#else
        write(out_unitp,31) dt_real,name
31      format('        real (s): ',f18.3,' in ',a)
#endif
        write(out_unitp,32) days,hours,minutes,seconds,name
32      format('        real    : ',i3,'d ',i2,'h ',i2,'m ',i2,'s in ',a)

         write(out_unitp,33) dt_cpu,name
33      format('        cpu (s): ',f18.3,' in ',a)


        !============================================
        !real and cpu total time
        seconds = int(t_real)
        minutes = seconds/60
        seconds = mod(seconds,60)
        hours   = minutes/60
        minutes = mod(minutes,60)
        days    = hours/24
        hours   = mod(hours,24)

#if(run_MPI)
        write(out_unitp,41) t_real,MPI_id
41      format('  Total real (s): ',f18.3,' from MPI id ',i4)
#else
        write(out_unitp,41) t_real
41      format('  Total real (s): ',f18.3)
#endif
        write(out_unitp,42) days,hours,minutes,seconds
42      format('  Total real    : ',i3,'d ',i2,'h ',i2,'m ',i2,'s')
        write(out_unitp,43) t_cpu
43      format('  Total cpu (s): ',f18.3)

        CALL flush_perso(out_unitp)
        !============================================


      END SUBROUTINE time_perso
!===============================================================================

      SUBROUTINE time_perso_v0(name)
      IMPLICIT NONE

        character (len=*) :: name


        integer :: tab_time(8) = 0
        real (kind=Rkind) :: t_real
        integer       :: count,count_work,freq,count_max
        real          :: t_cpu
        integer, save :: count_old,count_ini
        real, save    :: t_cpu_old,t_cpu_ini
        integer       :: seconds,minutes,hours,days
        logical, save :: begin = .TRUE.



        CALL date_and_time(values=tab_time)
        write(out_unitp,21) name,tab_time(5:8),tab_time(3:1:-1)
 21     format('     Time and date in ',a,' : ',i2,'h:',                &
               i2,'m:',i2,'.',i3,'s, the ',i2,'/',i2,'/',i4)

        CALL system_clock(count=count,count_rate=freq,count_max=count_max)
        call cpu_time(t_cpu)

        IF (begin) THEN
          begin = .FALSE.
          count_old = count
          count_ini = count
          t_cpu_old = t_cpu
          t_cpu_ini = t_cpu
        END IF


!       ============================================
!       cpu time in the subroutine: "name"

        !count_work = count-count_old
        count_work=merge(count-count_old,count-count_old+count_max,count>=count_old)
        seconds = count_work/freq

        minutes = seconds/60
        seconds = mod(seconds,60)
        hours   = minutes/60
        minutes = mod(minutes,60)
        days    = hours/24
        hours   = mod(hours,24)


        t_real = real(count_work,kind=Rkind)/real(freq,kind=Rkind)
        write(out_unitp,31) t_real,name
 31     format('        real (s): ',f18.3,' in ',a)
        write(out_unitp,32) days,hours,minutes,seconds,name
 32     format('        real    : ',i3,'d ',i2,'h ',i2,'m ',i2,'s in ',a)

         write(out_unitp,33) t_cpu-t_cpu_old,name
 33     format('        cpu (s): ',f18.3,' in ',a)


!       ============================================
!       Total cpu time

        !count_work = count-count_ini
        count_work=merge(count-count_ini,count-count_ini+count_max,count>=count_ini)
        seconds = count_work/freq

        minutes = seconds/60
        seconds = mod(seconds,60)
        hours   = minutes/60
        minutes = mod(minutes,60)
        days    = hours/24
        hours   = mod(hours,24)

        t_real = real(count_work,kind=Rkind)/real(freq,kind=Rkind)
#if(run_MPI)
        write(out_unitp,41) t_real,MPI_id
41      format('  Total real (s): ',f18.3,' from MPI id ',i4)
#else
        write(out_unitp,41) t_real
41      format('  Total real (s): ',f18.3)
#endif
        write(out_unitp,42) days,hours,minutes,seconds
42      format('  Total real    : ',i3,'d ',i2,'h ',i2,'m ',i2,'s')
        write(out_unitp,43) t_cpu-t_cpu_ini
43      format('  Total cpu (s): ',f18.3)

        write(out_unitp,51) '  Total memory: ',para_mem%mem_tot,' in ',name
51      format(a,i10,a,a)


        CALL flush_perso(out_unitp)
!       ============================================

        count_old = count
        t_cpu_old = t_cpu


      END SUBROUTINE time_perso_v0
      SUBROUTINE DeltaTime(dt_real,t_real,dt_cpu,t_cpu,LocalTime)
      IMPLICIT NONE

        real (kind=Rkind), intent(inout)           :: dt_real,t_real
        real,              intent(inout)           :: dt_cpu,t_cpu
        TYPE (param_time), intent(inout), optional :: LocalTime


        integer       :: count,count_work,freq

        TYPE (param_time), save :: MainTime

        IF (present(LocalTime)) THEN
          CALL DeltaTime_withParam_time(dt_real,t_real,dt_cpu,t_cpu,LocalTime)
        ELSE
          CALL DeltaTime_withParam_time(dt_real,t_real,dt_cpu,t_cpu,MainTime)
        END IF

      END SUBROUTINE DeltaTime
      
      SUBROUTINE DeltaTime_withParam_time(dt_real,t_real,dt_cpu,t_cpu,LocalTime)
      IMPLICIT NONE

        real (kind=Rkind), intent(inout) :: dt_real,t_real
        real,              intent(inout) :: dt_cpu,t_cpu
        TYPE (param_time), intent(inout) :: LocalTime


        integer       :: count,count_work,freq,count_max



        CALL system_clock(count=count,count_rate=freq,count_max=count_max)
        call cpu_time(t_cpu)

        IF (LocalTime%begin) THEN
          LocalTime%begin     = .FALSE.
          LocalTime%count_old = count
          LocalTime%count_ini = count
          LocalTime%t_cpu_old = t_cpu
          LocalTime%t_cpu_ini = t_cpu
        END IF


        ! real time
        !count_work = count-LocalTime%count_old
        count_work=merge(count-LocalTime%count_old,count-LocalTime%count_old+count_max,&
                         count>=LocalTime%count_old)
        dt_real    = real(count_work,kind=Rkind)/real(freq,kind=Rkind)
        !count_work = count-LocalTime%count_ini
        count_work=merge(count-LocalTime%count_ini,count-LocalTime%count_ini+count_max,&
                         count>=LocalTime%count_ini)
        t_real     = real(count_work,kind=Rkind)/real(freq,kind=Rkind)

        ! cpu time
        dt_cpu  = t_cpu-LocalTime%t_cpu_old
        t_cpu   = t_cpu-LocalTime%t_cpu_ini

        ! change the save variable
        LocalTime%count_old = count
        LocalTime%t_cpu_old = t_cpu

      END SUBROUTINE DeltaTime_withParam_time
      FUNCTION Delta_RealTime(LocalTime)
      IMPLICIT NONE

        TYPE (param_time), intent(inout), optional :: LocalTime

        real (kind=Rkind) :: Delta_RealTime

        real (kind=Rkind) :: dt_real,t_real
        real              :: dt_cpu,t_cpu

        IF (present(LocalTime)) THEN
          CALL DeltaTime(dt_real,t_real,dt_cpu,t_cpu,LocalTime)
        ELSE
          CALL DeltaTime(dt_real,t_real,dt_cpu,t_cpu)
        END IF

        Delta_RealTime = dt_real

      END FUNCTION Delta_RealTime

      !! @description: Compare two arrays of complex numbers
      !!               L1 and L2 of equal size
      !! @param: L1 First  array
      !! @param: L2 Second array
      logical FUNCTION compare_tab_cmplx(L1, L2)

       complex(kind=Rkind), intent(in) :: L1(:), L2(:)

       integer :: i

       if (size(L1) /= size(L2)) then
         compare_tab_cmplx = .false.
         return
       end if

       compare_tab_cmplx = .true.
       do i=1, size(L1)
         if (abs(L1(i)-L2(i)) > ONETENTH**13) then
           compare_tab_cmplx = .false.
            return
         end if
       end do

      END FUNCTION compare_tab_cmplx

      !! @description: Compare two arrays of real L1 and L2 of equal size
      !! @param: L1 First  array
      !! @param: L2 Second array
      logical FUNCTION compare_tab_real(L1, L2)

       real(kind=Rkind), intent(in) :: L1(:), L2(:)

       integer :: i

       if (size(L1) /= size(L2)) then
         compare_tab_real = .false.
         return
       end if

       compare_tab_real = .true.
       do i=1, size(L1)
         if (abs(L1(i)-L2(i)) > ONETENTH**13) then
           compare_tab_real = .false.
            return
         end if
       end do

      END FUNCTION compare_tab_real

      !! @description: Compare two arrays of integer L1 and L2 of equal size
      !! @param: L1 First  array
      !! @param: L2 Second array
      logical FUNCTION compare_tab_int(L1, L2)

       integer, intent(in) :: L1(:), L2(:)

       integer :: i

       if (size(L1) /= size(L2)) then
         compare_tab_int = .false.
         return
       end if

       compare_tab_int = .true.
       do i=1, size(L1)
         if (abs(L1(i)-L2(i)) /= 0) then
           compare_tab_int = .false.
            return
         end if
       end do

      END FUNCTION compare_tab_int

      !! @description: Compare two logical arrays a and b of equal size
      !! @param: L1 First logical array
      !! @param: L2 Second logical array
      logical FUNCTION compare_la(L1, L2)

       logical, intent(in) :: L1(:), L2(:)

       integer :: i

       if (size(L1) /= size(L2)) then
         compare_la = .false.
         return
       end if

       compare_la = .true.
       do i=1, size(L1)
         if (L1(i) .neqv. L2(i)) then
           compare_la = .false.
            return
         end if
       end do

      END FUNCTION compare_la

      logical FUNCTION inferior_tab_real(x1,x2)
      IMPLICIT NONE


      logical :: inf_loc
      integer       :: i
      real (kind=Rkind), intent(in) :: x1(:),x2(:)


       IF (size(x1) /= size(x2)) then
         write(out_unitp,*) 'the size of the tab are different !!'
         write(out_unitp,*) 'x1(:)',x1(:)
         write(out_unitp,*) 'x2(:)',x2(:)
         write(out_unitp,*) 'Check the fortran'
         STOP
       END IF

      inf_loc = .FALSE.

      DO i=1,size(x1)
        inf_loc = (x1(i) < x2(i))
        IF (x1(i) == x2(i)) CYCLE
        EXIT
      END DO

!     write(out_unitp,*) 'x1,x2,inf_loc',x1,x2,inf_loc

      inferior_tab_real = inf_loc

      END FUNCTION inferior_tab_real
      logical FUNCTION inferior_tab_int(x1,x2)
      IMPLICIT NONE


      logical :: inf_loc
      integer       :: i
      integer, intent(in) :: x1(:),x2(:)


       IF (size(x1) /= size(x2)) then
         write(out_unitp,*) 'the size of the tab are different !!'
         write(out_unitp,*) 'x1(:)',x1(:)
         write(out_unitp,*) 'x2(:)',x2(:)
         write(out_unitp,*) 'Check the fortran'
         STOP
       END IF

      inf_loc = .FALSE.

      DO i=1,size(x1)
        inf_loc = (x1(i) < x2(i))
        IF (x1(i) == x2(i)) CYCLE
        EXIT
      END DO

!     write(out_unitp,*) 'x1,x2,inf_loc',x1,x2,inf_loc

      inferior_tab_int = inf_loc

      END FUNCTION inferior_tab_int

      SUBROUTINE dihedral_range(angle,itype_dihedral)

        real (kind=Rkind), intent(inout) :: angle
        integer, optional :: itype_dihedral

        integer :: itype_dihedral_loc

        itype_dihedral_loc = 0
        IF (present(itype_dihedral)) itype_dihedral_loc = itype_dihedral

        SELECT CASE (itype_dihedral_loc)
        CASE (1) ! [-pi:pi]
          angle = modulo(angle,TWO*pi)
          IF (angle > pi) angle = angle - TWO*pi
        CASE (2) ! [0:2pi]
          angle = modulo(angle,TWO*pi)
        CASE Default
          ! nothing
        END SELECT

      END SUBROUTINE dihedral_range
END MODULE mod_system

