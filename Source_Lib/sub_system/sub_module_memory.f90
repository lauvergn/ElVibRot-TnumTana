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
      MODULE mod_memory
      USE QDUtil_NumParameters_m, only : Rkind, ILkind, out_unitp => out_unit
      IMPLICIT NONE

      PRIVATE

      TYPE param_memory
        integer (kind=ILkind)           :: max_mem      =  4000000000_ILkind/Rkind   ! (8GO) max_memory

        integer (kind=ILkind)           :: max_mem_used =  0   ! the maximal memory used
        integer (kind=ILkind)           :: mem_tot      =  0   ! memory used
        integer (kind=ILkind)           :: memory       =  0   ! asked memory

        integer                         :: nb_alloc     =  0   ! nb of allocations
        integer                        :: nb_dealloc   =  0   ! nb of deallocations

        logical                         :: mem_debug    = .FALSE.
        logical                         :: mem_print    = .FALSE.

        integer                         :: mem_unit     = -1
        character (len=:), allocatable  :: mem_file

      END TYPE param_memory

      TYPE (param_memory), save, public :: para_mem

      PUBLIC :: param_memory
      PUBLIC :: Write_error_NOT_null, Write_error_null, Write_mem_tot,Write_mem_file
      PUBLIC :: sub_test_tab_ub, sub_test_tab_lb, sub_test_Bigtab_ub, sub_test_Bigtab_lb
      PUBLIC :: Check_mem, UnCheck_mem
      PUBLIC :: error_memo_allo, error_lmemo_allo
      PUBLIC :: convertMem

      CONTAINS

      SUBROUTINE Write_error_NOT_null(name_sub_alloc,name_var,name_sub)
      IMPLICIT NONE


      character (len=*), intent(in) :: name_var,name_sub,name_sub_alloc

      write(out_unitp,*) ' ERROR in ',name_sub_alloc
      write(out_unitp,*) '   CALLED from "',name_sub,'" for the array "',name_var,'"'
      write(out_unitp,*) '   The array IS ALREADY allocated'
      write(out_unitp,*) '               OR'
      write(out_unitp,*) '   The pointer (array) IS ALREADY associated'
      write(out_unitp,*) '   => it cannot be allocated!'
      write(out_unitp,*) ' CHECK the fortran! '
      STOP

      END SUBROUTINE Write_error_NOT_null
      SUBROUTINE Write_error_null(name_sub_alloc,name_var,name_sub)
      IMPLICIT NONE


      character (len=*), intent(in) :: name_var,name_sub,name_sub_alloc

      write(out_unitp,*) ' ERROR in ',name_sub_alloc
      write(out_unitp,*) '   CALLED from "',name_sub,'" for the array "',name_var,'"'
      write(out_unitp,*) '   The array IS NOT allocated'
      write(out_unitp,*) '               OR'
      write(out_unitp,*) '   The pointer (array) IS NOT associated'
      write(out_unitp,*) '   => it cannot be used!'
      write(out_unitp,*) '   or it cannot be deallocated!'
      write(out_unitp,*) ' CHECK the fortran! '
      STOP

      END SUBROUTINE Write_error_null
      SUBROUTINE sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)
      IMPLICIT NONE

      integer, intent(in) :: ndim
      integer, intent(in) :: tab_ub(:)
      character (len=*), intent(in) :: name_var,name_sub,name_sub_alloc

       IF (sum(shape(tab_ub)) /= ndim) THEN
         write(out_unitp,*) ' ERROR in ',name_sub_alloc
         write(out_unitp,*) '   CALLED from "',name_sub,'" for the array "',name_var,'"'
         write(out_unitp,*) '   The number of element(s) of tab_ub(:) MUST be ',ndim
         write(out_unitp,*) '   The actual number IS ',shape(tab_ub)
         write(out_unitp,*) '      tab_ub(:) ',tab_ub(:)

         write(out_unitp,*) ' CHECK the fortran! '

         STOP
       END IF

      END SUBROUTINE sub_test_tab_ub
      SUBROUTINE sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)
      IMPLICIT NONE

      integer, intent(in) :: ndim
      integer, intent(in) :: tab_lb(:)
      character (len=*), intent(in) :: name_var,name_sub,name_sub_alloc

       IF (sum(shape(tab_lb)) /= ndim) THEN
         write(out_unitp,*) ' ERROR in ',name_sub_alloc
          write(out_unitp,*) '   CALLED from "',name_sub,'" for the array "',name_var,'"'
         write(out_unitp,*) '   The number of element(s) of tab_lb(:) MUST be ',ndim
         write(out_unitp,*) '   The actual number IS ',shape(tab_lb)
         write(out_unitp,*) '      tab_lb(:) ',tab_lb(:)

         write(out_unitp,*) ' CHECK the fortran! '

         STOP
       END IF

      END SUBROUTINE sub_test_tab_lb
      SUBROUTINE sub_test_Bigtab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)
      IMPLICIT NONE

      integer,               intent(in) :: ndim
      integer (kind=ILkind), intent(in) :: tab_ub(:)
      character (len=*),     intent(in) :: name_var,name_sub,name_sub_alloc

       IF (sum(shape(tab_ub)) /= ndim) THEN
         write(out_unitp,*) ' ERROR in ',name_sub_alloc
         write(out_unitp,*) '   CALLED from "',name_sub,'" for the array "',name_var,'"'
         write(out_unitp,*) '   The number of element(s) of tab_ub(:) MUST be ',ndim
         write(out_unitp,*) '   The actual number IS ',shape(tab_ub)
         write(out_unitp,*) '      tab_ub(:) ',tab_ub(:)

         write(out_unitp,*) ' CHECK the fortran! '

         STOP
       END IF

      END SUBROUTINE sub_test_Bigtab_ub
      SUBROUTINE sub_test_Bigtab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)
      IMPLICIT NONE

      integer,               intent(in) :: ndim
      integer (kind=ILkind), intent(in) :: tab_lb(:)
      character (len=*),     intent(in) :: name_var,name_sub,name_sub_alloc

       IF (sum(shape(tab_lb)) /= ndim) THEN
         write(out_unitp,*) ' ERROR in ',name_sub_alloc
          write(out_unitp,*) '   CALLED from "',name_sub,'" for the array "',name_var,'"'
         write(out_unitp,*) '   The number of element(s) of tab_lb(:) MUST be ',ndim
         write(out_unitp,*) '   The actual number IS ',shape(tab_lb)
         write(out_unitp,*) '      tab_lb(:) ',tab_lb(:)

         write(out_unitp,*) ' CHECK the fortran! '

         STOP
       END IF

      END SUBROUTINE sub_test_Bigtab_lb
      SUBROUTINE Check_mem()
      IMPLICIT NONE



      para_mem%mem_print = .FALSE.
      para_mem%mem_debug = .TRUE.

      write(out_unitp,*) '=============================================='
      write(out_unitp,*) '========= CHECK MEMORY ======================='
      write(out_unitp,*) '=============================================='
      write(out_unitp,*) 'mem_tot,max_mem_used',para_mem%mem_tot,para_mem%max_mem_used
      write(out_unitp,*) 'nb_alloc,nb_dealloc',para_mem%nb_alloc,para_mem%nb_dealloc
      write(out_unitp,*) '=============================================='

      END SUBROUTINE Check_mem
      SUBROUTINE UnCheck_mem()
      IMPLICIT NONE

      para_mem%mem_print = .TRUE.
      para_mem%mem_debug = .FALSE.

      write(out_unitp,*) '=============================================='
      write(out_unitp,*) '========= UNCHECK MEMORY ====================='
      write(out_unitp,*) '=============================================='
      write(out_unitp,*) 'mem_tot,max_mem_used',para_mem%mem_tot,para_mem%max_mem_used
      write(out_unitp,*) 'nb_alloc,nb_dealloc',para_mem%nb_alloc,para_mem%nb_dealloc
      write(out_unitp,*) '=============================================='

      END SUBROUTINE UnCheck_mem
      SUBROUTINE error_memo_allo(err,memory,name_var,name_sub,var_type)
      IMPLICIT NONE

      integer, intent(in) :: err
      integer, intent(in) :: memory
      character (len=*), intent(in) :: name_var,name_sub
      character (len=*), intent(in), optional :: var_type


      logical :: memory_test

!----- for debuging --------------------------------------------------
      logical,parameter :: debug=.FALSE.
      ! logical,parameter :: debug=.TRUE.

!      para_mem%mem_tot = para_mem%mem_tot + int(memory,kind=ILkind)
!      IF (para_mem%mem_tot > para_mem%max_mem_used)                     &
!                               para_mem%max_mem_used = para_mem%mem_tot
!      IF (memory > 0) THEN
!        para_mem%nb_alloc = para_mem%nb_alloc + 1
!      ELSE IF (memory < 0) THEN
!        para_mem%nb_dealloc = para_mem%nb_dealloc + 1
!      END IF

!      IF (abs(memory) > 10**6-1) THEN
!!$OMP CRITICAL (error_memo_allo_CRIT)
!      IF (present(var_type)) THEN
!        write(999,*) para_mem%mem_tot,para_mem%mem_tot+int(memory,kind=ILkind),&
!                memory,' var_type=',var_type,' name_var=',name_var,' ',name_sub
!        !write(out_unitp,*) para_mem%mem_tot,memory,' var_type=',var_type,&
!        !      ' name_var=',name_var,' ',name_sub
!      ELSE
!        write(999,*) para_mem%mem_tot,para_mem%mem_tot+int(memory,kind=ILkind),&
!                         memory,' no_var_type name_var=',name_var,' ',name_sub
!        !write(out_unitp,*) para_mem%mem_tot,memory,' no_var_type name_var=',  &
!        !        name_var,' ',name_sub
!      END IF
!      flush(999)
!!$OMP END CRITICAL (error_memo_allo_CRIT)
!      END IF

      IF (err /= 0) THEN
        write(out_unitp,*) ' ERROR in  the routine "',name_sub,'"'
        write(out_unitp,*) ' cannot allocate or deallocate the ',       &
                                             'variable: "',name_var,'"'
        IF (present(var_type) ) write(out_unitp,*) 'Variable TYPE: ',var_type

        write(out_unitp,*) 'err_mem',err
        write(out_unitp,*) 'mem_tot,memory',para_mem%mem_tot,memory
        STOP
      END IF

      IF (.NOT. debug .AND. .NOT. para_mem%mem_debug) RETURN

!$OMP CRITICAL (error_memo_allo_CRIT)
      CALL open_mem_file()

      IF (present(var_type)) THEN
        write(para_mem%mem_unit,*) para_mem%mem_tot,para_mem%mem_tot+int(memory,kind=ILkind),&
                memory,' var_type=',var_type,' name_var=',name_var,' ',name_sub
      ELSE
        write(para_mem%mem_unit,*) para_mem%mem_tot,para_mem%mem_tot+int(memory,kind=ILkind),&
                         memory,' no_var_type name_var=',name_var,' ',name_sub
      END IF
      flush(para_mem%mem_unit)

      IF (memory > 0) THEN
        memory_test = ( para_mem%mem_tot > huge(1_ILkind)-int(memory,kind=ILkind) )
        IF (memory_test) THEN
          write(out_unitp,*) ' ERROR in error_memo_allo'
          write(out_unitp,*) ' Variable and subroutine: ',                      &
                                name_var,' in ',name_sub
          IF (present(var_type) ) write(out_unitp,*) 'Variable TYPE: ',var_type

          write(out_unitp,*) ' mem_tot WILL be larger than the largest integer!'
          write(out_unitp,*) ' mem_tot,memory,huge(int)',para_mem%mem_tot,      &
                                                         memory,huge(1_ILkind)
          write(out_unitp,*) ' huge(int)-memory',huge(1_ILkind)-int(memory,kind=ILkind)
          write(out_unitp,*) ' => calculation TOO large or memory leak'
          STOP
        END IF
      END IF

      para_mem%mem_tot = para_mem%mem_tot + int(memory,kind=ILkind)
      IF (para_mem%mem_tot > para_mem%max_mem_used)                     &
                               para_mem%max_mem_used = para_mem%mem_tot
      IF (memory > 0) THEN
        para_mem%nb_alloc = para_mem%nb_alloc + 1
      ELSE IF (memory < 0) THEN
        para_mem%nb_dealloc = para_mem%nb_dealloc + 1
      END IF

      IF (para_mem%mem_print) THEN
        IF (memory == 0 .AND. err == 0) THEN
          write(out_unitp,*) ' WARNING memory = 0'
          write(out_unitp,*) 'old_mem_tot,new_mem_tot,memory',          &
                              para_mem%mem_tot-int(memory,kind=ILkind),para_mem%mem_tot, &
                              memory,'nothing: ',name_var,' in ',name_sub
          IF (present(var_type) ) write(out_unitp,*) 'Variable TYPE: ',var_type
        END IF

        IF ((debug .OR. para_mem%mem_debug) .AND. memory < 0 .AND. err == 0) THEN
          write(out_unitp,*) 'old_mem_tot,new_mem_tot,memory',          &
                              para_mem%mem_tot-int(memory,kind=ILkind),para_mem%mem_tot, &
                              memory,' dealloc of var name "',name_var, &
                              '" in routine "',name_sub,'"'
          IF (present(var_type) ) write(out_unitp,*) 'Variable TYPE: ',var_type
        END IF
        IF ((debug .OR. para_mem%mem_debug) .AND. memory > 0 .AND. err == 0) THEN
          write(out_unitp,*) 'old_mem_tot,new_mem_tot,memory',          &
                              para_mem%mem_tot-int(memory,kind=ILkind),para_mem%mem_tot, &
                              memory,' alloc of var name "',name_var,   &
                              '" in routine "',name_sub,'"'
          IF (present(var_type) ) write(out_unitp,*) 'Variable TYPE: ',var_type
        END IF

        IF ((debug .OR. para_mem%mem_debug) .AND. err == 0) THEN
          write(out_unitp,*) 'max_mem_used,mem_tot,nb_alloc,nb_dealloc',&
                                para_mem%max_mem_used,para_mem%mem_tot, &
                                 para_mem%nb_alloc,para_mem%nb_dealloc, &
                 ' of var name "',name_var,'" in routine "',name_sub,'"'
        END IF

      END IF

!$OMP END CRITICAL (error_memo_allo_CRIT)
      END SUBROUTINE error_memo_allo

      SUBROUTINE error_lmemo_allo(err,memory,name_var,name_sub,var_type)
      IMPLICIT NONE

      integer,               intent(in)            :: err
      integer (kind=ILkind), intent(in)            :: memory
      character (len=*),     intent(in)            :: name_var,name_sub
      character (len=*),     intent(in), optional :: var_type


      logical :: memory_test

!----- for debuging --------------------------------------------------
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.

      IF (err /= 0) THEN
        write(out_unitp,*) ' ERROR in  the routine "',name_sub,'"'
        write(out_unitp,*) ' cannot allocate or deallocate the ',       &
                                             'variable: "',name_var,'"'
        IF (present(var_type) ) write(out_unitp,*) 'Variable TYPE: ',var_type

        write(out_unitp,*) 'err_mem',err
        write(out_unitp,*) 'mem_tot,memory',para_mem%mem_tot,memory
        STOP
      END IF

      IF (.NOT. debug .AND. .NOT. para_mem%mem_debug) RETURN

!$OMP CRITICAL (error_lmemo_allo_CRIT)

      CALL open_mem_file()

      IF (present(var_type)) THEN
        write(para_mem%mem_unit,*) para_mem%mem_tot,para_mem%mem_tot+memory,memory,&
                        ' var_type=',var_type,' name_var=',name_var,' ',name_sub
      ELSE
        write(para_mem%mem_unit,*) para_mem%mem_tot,para_mem%mem_tot+memory,memory,&
                        ' no_var_type name_var=',name_var,' ',name_sub
      END IF
      flush(para_mem%mem_unit)

      IF (memory > 0) THEN
        memory_test = ( para_mem%mem_tot > huge(1_ILkind)-memory )
        IF (memory_test) THEN
          write(out_unitp,*) ' ERROR in error_memo_allo'
          write(out_unitp,*) ' Variable and subroutine: ',                      &
                                name_var,' in ',name_sub
          IF (present(var_type) ) write(out_unitp,*) 'Variable TYPE: ',var_type

          write(out_unitp,*) ' mem_tot WILL be larger than the largest integer!'
          write(out_unitp,*) ' mem_tot,memory,huge(int)',para_mem%mem_tot,      &
                                                         memory,huge(1_ILkind)
          write(out_unitp,*) ' huge(int)-memory',huge(1_ILkind)-memory
          write(out_unitp,*) ' => calculation TOO large or memory leak'
          STOP
        END IF
      END IF

      para_mem%mem_tot = para_mem%mem_tot + memory
      IF (para_mem%mem_tot > para_mem%max_mem_used)                     &
                               para_mem%max_mem_used = para_mem%mem_tot
      IF (memory > 0) THEN
        para_mem%nb_alloc = para_mem%nb_alloc + 1
      ELSE IF (memory < 0) THEN
        para_mem%nb_dealloc = para_mem%nb_dealloc + 1
      END IF

      IF (para_mem%mem_print) THEN
        IF (memory == 0 .AND. err == 0) THEN
          write(out_unitp,*) ' WARNING memory = 0'
          write(out_unitp,*) 'old_mem_tot,new_mem_tot,memory',          &
                              para_mem%mem_tot-memory,para_mem%mem_tot, &
                              memory,'nothing: ',name_var,' in ',name_sub
          IF (present(var_type) ) write(out_unitp,*) 'Variable TYPE: ',var_type
        END IF

        IF ((debug .OR. para_mem%mem_debug) .AND. memory < 0 .AND. err == 0) THEN
          write(out_unitp,*) 'old_mem_tot,new_mem_tot,memory',          &
                              para_mem%mem_tot-memory,para_mem%mem_tot, &
                              memory,' dealloc of var name "',name_var, &
                              '" in routine "',name_sub,'"'
          IF (present(var_type) ) write(out_unitp,*) 'Variable TYPE: ',var_type
        END IF
        IF ((debug .OR. para_mem%mem_debug) .AND. memory > 0 .AND. err == 0) THEN
          write(out_unitp,*) 'old_mem_tot,new_mem_tot,memory',          &
                              para_mem%mem_tot-memory,para_mem%mem_tot, &
                              memory,' alloc of var name "',name_var,   &
                              '" in routine "',name_sub,'"'
          IF (present(var_type) ) write(out_unitp,*) 'Variable TYPE: ',var_type
        END IF

        IF ((debug .OR. para_mem%mem_debug) .AND. err == 0) THEN
          write(out_unitp,*) 'max_mem_used,mem_tot,nb_alloc,nb_dealloc',&
                                para_mem%max_mem_used,para_mem%mem_tot, &
                                 para_mem%nb_alloc,para_mem%nb_dealloc, &
                 ' of var name "',name_var,'" in routine "',name_sub,'"'
        END IF

      END IF

!$OMP END CRITICAL (error_lmemo_allo_CRIT)
      END SUBROUTINE error_lmemo_allo

 SUBROUTINE Write_mem_tot(info)
 IMPLICIT NONE
 character (len=*), optional, intent(in) :: info


!$OMP CRITICAL (Write_mem_tot_CRIT)
   IF (present(info)) THEN
     write(out_unitp,*) 'max_mem_used,mem_tot,nb_alloc,nb_dealloc',     &
                                para_mem%max_mem_used,para_mem%mem_tot, &
                                 para_mem%nb_alloc,para_mem%nb_dealloc, &
                                 ' info: ',info
   ELSE
     write(out_unitp,*) 'max_mem_used,mem_tot,nb_alloc,nb_dealloc',     &
                                para_mem%max_mem_used,para_mem%mem_tot, &
                                 para_mem%nb_alloc,para_mem%nb_dealloc
   END IF
!$OMP END CRITICAL (Write_mem_tot_CRIT)
 END SUBROUTINE Write_mem_tot

 SUBROUTINE Write_mem_file(info)
 IMPLICIT NONE

 character (len=*), intent(in) :: info

 IF (.NOT. para_mem%mem_debug) RETURN

 !$OMP CRITICAL (write_mem_file_CRIT)
 CALL open_mem_file()

 write(para_mem%mem_unit,*) "--------------------------------------------------"
 write(para_mem%mem_unit,*) "--------------------------------------------------"
 write(para_mem%mem_unit,*)
 write(para_mem%mem_unit,*) "Memory analysis: ",info
 write(para_mem%mem_unit,*)
 write(para_mem%mem_unit,*) 'max_mem_used,mem_tot,nb_alloc,nb_dealloc',         &
                            para_mem%max_mem_used,para_mem%mem_tot,             &
                             para_mem%nb_alloc,para_mem%nb_dealloc
 write(para_mem%mem_unit,*) "--------------------------------------------------"
 write(para_mem%mem_unit,*) "--------------------------------------------------"

 flush(para_mem%mem_unit)

 !$OMP END CRITICAL (write_mem_file_CRIT)

END SUBROUTINE Write_mem_file

 SUBROUTINE open_mem_file()
 IMPLICIT NONE

 IF (para_mem%mem_unit == -1) THEN
   open(newunit=para_mem%mem_unit,file='EVRT_mem.log')
   write(out_unitp,*) 'Memory file: "EVRT_mem.log", unit:',para_mem%mem_unit
   IF (para_mem%mem_unit == -1) STOP 'ERROR cannot open the file: "EVRT_mem.log"'
 END IF

 END SUBROUTINE open_mem_file

 SUBROUTINE convertMem(mem,MemUnit)
 real(kind=Rkind),   intent(inout) :: mem
 character (len=2),  intent(inout) :: MemUnit

   IF (mem < 1024_Rkind) THEN
     MemUnit='O'
   ELSE IF (mem < 1024_Rkind**2) THEN
     mem = mem / 1024_Rkind
     MemUnit='kO'
   ELSE IF (mem < 1024_Rkind**3) THEN
     mem = mem / 1024_Rkind**2
     MemUnit='MO'
   ELSE
     mem = mem / 1024_Rkind**3
     MemUnit='GO'
   END IF
 END SUBROUTINE convertMem
END MODULE mod_memory
