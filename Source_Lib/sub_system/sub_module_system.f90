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
      USE QDUtil_m, out_unitp => out_unit, in_unitp => in_unit
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
#if defined(__GIT)
      character (len=Line_len) :: git_branch = __GIT
#else
      character (len=Line_len) :: git_branch = "unknown: -D__GIT=?"
#endif


      logical :: openmp = .FALSE.
      logical :: openmpi= .FALSE.
      integer :: MatOp_omp,OpPsi_omp,BasisTOGrid_omp,Grid_omp,SG4_omp
      integer :: MatOp_maxth,OpPsi_maxth,BasisTOGrid_maxth,Grid_maxth,SG4_maxth

      integer :: CRP_omp,CRP_maxth

      logical :: Tune_SG4_omp  = .FALSE.
      logical :: Tune_Grid_omp = .FALSE.

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
        logical :: Grid_only        = .FALSE.
        logical :: cart             = .FALSE.
        logical :: GridTOBasis_test = .FALSE.
        logical :: OpPsi_test       = .FALSE.
        logical :: nDfit            = .FALSE.
        logical :: nDGrid           = .FALSE.
        logical :: main_test        = .FALSE.

      END TYPE param_EVRT_calc

      TYPE (param_FOR_optimization), save :: para_FOR_optimization
      TYPE (param_EVRT_calc),        save :: para_EVRT_calc

      CONTAINS

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

!=======================================================================================
!> @brief subroutine for recording time
!> @param time_sum should be initialized before calling this function
!=======================================================================================
      SUBROUTINE time_record(time_sum,time1,time2,point)
        USE mod_MPI
        IMPLICIT NONE

        Integer,                        intent(inout) :: time_sum
        Integer,                        intent(inout) :: time1
        Integer,                        intent(inout) :: time2
        Integer,                           intent(in) :: point

        IF(point==1) THEN
          CALL system_clock(time1,time_rate,time_max)
        ELSEIF(point==2) THEN
          CALL system_clock(time2,time_rate,time_max)
          time_sum=time_sum+merge(time2-time1,time2-time1+time_max,time2>=time1)
        ELSE
          STOP 'error when calling time_record'
        ENDIF
      ENDSUBROUTINE
  SUBROUTINE versionEVRT(write_version)
    IMPLICIT NONE
  
        logical :: write_version
  
        character (len=*), parameter :: EVR_name='ElVibRot'
        character (len=*), parameter :: Tnum_name='Tnum'
        character (len=*), parameter :: Tana_name='Tana'
  
  
  
        IF (write_version .AND. MPI_id==0) THEN
          write(out_unitp,*) '==============================================='
          write(out_unitp,*) '==============================================='
          write(out_unitp,*) 'Working with ',                             &
                     EVR_name,trim(adjustl(EVR_version)),'-',             &
                     Tnum_name,trim(adjustl(Tnum_version)),'-',           &
                     Tana_name,trim(adjustl(Tana_version))
  
          write(out_unitp,*) 'Compiled on "',trim(compile_host), '" the ',trim(compile_date)
          write(out_unitp,*) 'Compiler version: ',trim(compiler_ver)
          write(out_unitp,*) 'Compiler options: ',trim(compiler_opt)
          write(out_unitp,*) 'Compiler libs: ',trim(compiler_libs)
  
          write(out_unitp,*) 'EVRT_path: ',trim(EVRT_path)
          write(out_unitp,*) 'git ',trim(git_branch)
  
          write(out_unitp,*) '-----------------------------------------------'
  
          write(out_unitp,*) EVR_name,' is written by David Lauvergnat [1] '
          write(out_unitp,*) '  with contributions of'
          write(out_unitp,*) '     Josep Maria Luis (optimization) [2]'
          write(out_unitp,*) '     Ahai Chen (MPI) [1,4]'
          write(out_unitp,*) '     Lucien Dupuy (CRP) [5]'
  
          write(out_unitp,*) EVR_name,' is under GNU LGPL3 license.'
          write(out_unitp,*)
  
          write(out_unitp,*) Tnum_name,' is written David Lauvergnat [1]'
          write(out_unitp,*) Tana_name,' is written by Mamadou Ndong [1] and David Lauvergnat [1]'
          write(out_unitp,*) '  with contributions'
          write(out_unitp,*) '      Emil Lund klinting (coupling with MidasCpp) [3]'
  
          write(out_unitp,*) Tnum_name,' and ',Tana_name,' are under GNU LGPL3 license.'
          write(out_unitp,*)
          write(out_unitp,*) '[1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France'
          write(out_unitp,*) '[2]: Institut de Química Computacional and Departament de Química',&
                                     ' Universitat de Girona, Catalonia, Spain'
          write(out_unitp,*) '[3]: Department of Chemistry, Aarhus University, DK-8000 Aarhus C, Denmark'
          write(out_unitp,*) '[4]: Maison de la Simulation USR 3441, CEA Saclay, France'
          write(out_unitp,*) '[5]: Laboratoire Univers et Particule de Montpellier, UMR 5299,', &
                                     ' Université de Montpellier, France'
          write(out_unitp,*) '==============================================='
          write(out_unitp,*) '==============================================='
        END IF
  END SUBROUTINE versionEVRT

END MODULE mod_system
