!===========================================================================
!===========================================================================
!This file is part of ElVibRot.
!
!    ElVibRot is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ElVibRot is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ElVibRot.  If not, see <http://www.gnu.org/licenses/>.
!
!    Copyright 2015  David Lauvergnat
!      with contributions of Mamadou Ndong, Josep Maria Luis
!
!    ElVibRot includes:
!        - Tnum-Tana under the GNU LGPL3 license
!        - Somme subroutines of John Burkardt under GNU LGPL license
!             http://people.sc.fsu.edu/~jburkardt/
!        - Somme subroutines of SHTOOLS written by Mark A. Wieczorek under BSD license
!             http://shtools.ipgp.fr
!===========================================================================
!===========================================================================

MODULE mod_file
      use mod_NumParameters, only: line_len, name_len, Name_longlen, out_unitp
      use mod_string
      !$ USE omp_lib, only : OMP_GET_THREAD_NUM
      IMPLICIT NONE

      PRIVATE

      character (len=Line_len), public :: File_path = ''

      !!@description: TODO
      !!@param: TODO
      TYPE param_file
        character (len=Line_len) :: name      = " "     ! name of the file
        integer                  :: unit      = 0       ! unit of the file
        logical                  :: formatted = .TRUE.
        logical                  :: append    =.FALSE.
        logical                  :: old       =.FALSE.
        logical                  :: seq       = .TRUE.
        integer                  :: frecl     = 0
        logical                  :: init      = .FALSE.


        ! to store/use data of several files using threads
        integer                               :: nb_thread      = 0
        character (len=Line_len), allocatable :: tab_name_th(:)
        integer, allocatable                  :: tab_unit(:)
      CONTAINS
        PROCEDURE, PRIVATE, PASS(file1) :: file2TOfile1
        GENERIC,   PUBLIC  :: assignment(=) => file2TOfile1
      END TYPE param_file

      PUBLIC :: param_file,file_GetUnit, file_open, file_open2
      PUBLIC :: file_close, file_delete, file_dealloc, file_write, make_FileName
      PUBLIC :: err_file_name,check_file_exist_WITH_file_name
      PUBLIC :: flush_perso,join_path
      PUBLIC :: exit_Davidson_external

      CONTAINS

      SUBROUTINE file2TOfile1(file1,file2)
      CLASS(param_file), intent(inout)  :: file1
      TYPE(param_file),  intent(in)     :: file2

      integer :: err_mem,memory
      !write(out_unitp,*) 'BEGINNING file_GetUnit'

        ! to store/use data of several files using threads
        integer                               :: nb_thread      = 0
        character (len=Line_len), allocatable :: tab_name_th(:)
        integer, allocatable                  :: tab_unit(:)
      !IF (.NOT. file2%init) RETURN

      file1%name      = file2%name
      file1%unit      = file2%unit
      file1%formatted = file2%formatted
      file1%append    = file2%append
      file1%old       = file2%old
      file1%seq       = file2%seq
      file1%frecl     = file2%frecl
      file1%init      = file2%init

      IF (allocated(file2%tab_unit) .AND. allocated(file2%tab_name_th)) THEN
        memory = file1%nb_thread
        allocate(file1%tab_unit(0:file1%nb_thread-1),stat=err_mem) ! change alloc done

        memory = file1%nb_thread
        allocate(file1%tab_name_th(0:file1%nb_thread-1),stat=err_mem) ! change alloc done

        file1%tab_unit       = file2%tab_unit
        file1%tab_name_th    = file2%tab_name_th

      END IF


      END SUBROUTINE file2TOfile1

      FUNCTION err_file_name(file_name,name_sub)
      integer                                 :: err_file_name
      character (len=*), intent(in)           :: file_name
      character (len=*), intent(in), optional :: name_sub

        IF (string_IS_empty(file_name) ) THEN
          IF (present(name_sub)) THEN
            write(out_unitp,*) ' ERROR in err_file_name, called from: ',name_sub
          ELSE
            write(out_unitp,*) ' ERROR in err_file_name'
          END IF
          write(out_unitp,*) '   The file name is empty'
          err_file_name = 1
        ELSE
          err_file_name = 0
        END IF

      END FUNCTION err_file_name

      FUNCTION check_file_exist_WITH_file_name(file_name,err_file,name_sub)
      logical                                    :: check_file_exist_WITH_file_name
      character (len=*), intent(in)              :: file_name
      character (len=*), intent(in)   , optional :: name_sub
      integer,           intent(inout), optional :: err_file

      logical                          :: file_exist


        IF (string_IS_empty(file_name) ) THEN
          IF (present(name_sub)) THEN
            write(out_unitp,*) ' ERROR in check_file_exist_WITH_file_name, called from: ',name_sub
          ELSE
            write(out_unitp,*) ' ERROR in check_file_exist_WITH_file_name'
          END IF
          write(out_unitp,*) '   The file name is empty'
          IF (present(err_file)) err_file = 1
          file_exist = .FALSE.
        ELSE
          IF (present(err_file)) err_file = 0
          INQUIRE(FILE=trim(file_name), EXIST=file_exist)
        END IF
        check_file_exist_WITH_file_name = file_exist

      END FUNCTION check_file_exist_WITH_file_name

      FUNCTION file_GetUnit(ffile,err_file)

      integer           :: file_GetUnit
      TYPE(param_file)  :: ffile
      integer, optional :: err_file


      logical                  :: unit_opened
      integer                  :: ithread,err_file_loc

      !write(out_unitp,*) 'BEGINNING file_GetUnit'

      IF (.NOT. ffile%init) THEN
        file_GetUnit = 0 ! the file is not opened
        RETURN
      END IF


!     - check if the file is already open ------------------
      IF (ffile%nb_thread > 1) THEN
         ithread      = 0
!$       ithread      = OMP_GET_THREAD_NUM()

        err_file_loc = err_file_name(ffile%tab_name_th(ithread),'file_GetUnit')
        IF (.NOT. present(err_file) .AND. err_file_loc /= 0) STOP ' ERROR, the file name is empty!'
        IF (present(err_file)) err_file = err_file_loc


        inquire(FILE=ffile%tab_name_th(ithread),OPENED=unit_opened)
        IF (.NOT. unit_opened) THEN ! the file is not open
          file_GetUnit = 0 ! the file is not opened
        ELSE
          file_GetUnit = ffile%tab_unit(ithread)
        END IF

      ELSE

        err_file_loc = err_file_name(ffile%name,'file_GetUnit')
        IF (.NOT. present(err_file) .AND. err_file_loc /= 0) STOP ' ERROR, the file name is empty!'
        IF (present(err_file)) err_file = err_file_loc

        inquire(FILE=ffile%name,OPENED=unit_opened)
        IF (.NOT. unit_opened) THEN ! the file is not open
          file_GetUnit = 0 ! the file is not opened
        ELSE
          file_GetUnit = ffile%unit
        END IF
      END IF

      END FUNCTION file_GetUnit

  FUNCTION GetUnit_NewFile(file_name,err_file)

      integer                                    :: GetUnit_NewFile
      character (len=*), intent(in)              :: file_name
      integer,           intent(inout), optional :: err_file


      logical                  :: unit_opened
      integer                  :: iunit,err_file_loc

      !write(out_unitp,*) 'BEGINNING GetUnit_NewFile'


      err_file_loc = err_file_name(file_name,'GetUnit_NewFile')
      IF (.NOT. present(err_file) .AND. err_file_loc /= 0) STOP ' ERROR, the file name is empty!'
      IF (present(err_file)) err_file = err_file_loc


      !- check if the file is already open ------------------
      inquire(FILE=file_name,NUMBER=iunit,OPENED=unit_opened)
      IF (.NOT. unit_opened) THEN ! the file is not open

        !- the file is not open, find an unused UNIT ---------
        iunit = 66
        DO
          iunit = iunit + 1
          inquire(UNIT=iunit,OPENED=unit_opened)
          IF (.NOT. unit_opened) EXIT
        END DO
        GetUnit_NewFile = iunit

      ELSE
        GetUnit_NewFile = 0
      END IF

  END FUNCTION GetUnit_NewFile


      !!@description: TODO
      !!@param: TODO
      SUBROUTINE file_open(ffile,iunit,lformatted,append,old,seq,lrecl,err_file)
      USE mod_string, ONLY : int_TO_char

      TYPE(param_file)  :: ffile
      integer           :: iunit
      integer, optional :: lrecl
      logical, optional :: lformatted,append,old,seq
      integer, optional :: err_file


      character (len=Name_len) :: fform,fstatus,fposition,faccess

      logical                  :: unit_opened
      integer                  :: ith,err_file_loc

      integer :: err_mem,memory
      !write(out_unitp,*) 'BEGINNING file_open'

!     - test if optional arguments are present ---------
      IF (.NOT. ffile%init) THEN  ! IF init=.T., those parameters are already set-up


        IF (present(lformatted)) THEN
          ffile%formatted = lformatted
        ELSE
          ffile%formatted = .TRUE.
        END IF

        IF (present(old)) THEN
          ffile%old = old
        ELSE
          ffile%old = .FALSE.
        END IF

        IF (present(append)) THEN
          ffile%append = append
        ELSE
          ffile%append = .FALSE.
        END IF


        IF (present(seq)) THEN
          IF (seq) THEN
            ffile%seq = .TRUE.
          ELSE
            ffile%seq = .FALSE.
            IF (present(lrecl)) THEN
              ffile%frecl = lrecl
            ELSE
              write(out_unitp,*) 'ERROR in file_open'
              write(out_unitp,*) 'The file access is DIRECT but lrecl is not present!'
              STOP
            END IF
          END IF
        ELSE
          ffile%seq = .TRUE.
        END IF



        ffile%init  = .TRUE.
      END IF
      !-------------------------------------


      !-------------------------------------
      IF (ffile%formatted) THEN
        fform = 'formatted'
      ELSE
        fform = 'unformatted'
      END IF

      IF (ffile%old) THEN
        fstatus = 'old'
      ELSE
        fstatus = 'unknown'
      END IF

      IF (ffile%append) THEN
        fposition = 'append'
      ELSE
        fposition = 'asis'
      END IF

      IF (ffile%seq) THEN
        faccess = 'sequential'
      ELSE
        faccess = 'direct'
      END IF
      !-------------------------------------
      !CALL file_Write(ffile)


      IF (.NOT. ffile%seq) ffile%nb_thread = 0

      err_file_loc = err_file_name(ffile%name,'file_open')
      IF (.NOT. present(err_file) .AND. err_file_loc /= 0) STOP ' ERROR, the file name is empty!'
      IF (present(err_file)) err_file = err_file_loc

      !write(out_unitp,*) 'ffile%name,iunit ',ffile%name,iunit ; flush(out_unitp)


      !- check if the file is already open ------------------
      inquire(FILE=ffile%name,NUMBER=iunit,OPENED=unit_opened)
      IF (.NOT. unit_opened) THEN ! the file is not open

        !- the file is not open, find an unused UNIT ---------
        iunit = 9
        DO
          iunit = iunit + 1
          inquire(UNIT=iunit,OPENED=unit_opened)
          IF (.NOT. unit_opened) EXIT
        END DO


        !-- open the file
        IF (ffile%seq) THEN
          IF (present(err_file)) THEN
            open(UNIT=iunit,FILE=trim(ffile%name),FORM=trim(fform),STATUS=trim(fstatus),  &
                 POSITION=trim(fposition),ACCESS='SEQUENTIAL',IOSTAT=err_file)
            IF (err_file /= 0) RETURN
          ELSE
            open(UNIT=iunit,FILE=trim(ffile%name),FORM=trim(fform),STATUS=trim(fstatus),  &
                 POSITION=trim(fposition),ACCESS='SEQUENTIAL')
          END IF
        ELSE
          IF (present(err_file)) THEN
            open(UNIT=iunit,FILE=ffile%name,FORM=fform,STATUS=fstatus,  &
                       ACCESS='DIRECT',RECL=ffile%frecl,IOSTAT=err_file)
            IF (err_file /= 0) RETURN
          ELSE
            open(UNIT=iunit,FILE=ffile%name,FORM=fform,STATUS=fstatus,  &
                 ACCESS='DIRECT',RECL=ffile%frecl)
          END IF
        END IF

        ffile%unit = iunit
      ELSE
        ffile%unit = iunit
      END IF

      !write(out_unitp,*) 'open ',iunit,ffile%name

      IF (ffile%nb_thread > 1) THEN
        IF (allocated(ffile%tab_unit)) deallocate(ffile%tab_unit,stat=err_mem)
        memory = ffile%nb_thread
        allocate(ffile%tab_unit(0:ffile%nb_thread-1),stat=err_mem) ! change alloc done
!        CALL error_memo_allo(err_mem,memory,"ffile%tab_unit",           &
!                                                            "file_open")
        IF (allocated(ffile%tab_name_th)) deallocate(ffile%tab_name_th,stat=err_mem)
        memory = ffile%nb_thread
        allocate(ffile%tab_name_th(0:ffile%nb_thread-1),stat=err_mem) ! change alloc done
!        CALL error_memo_allo(err_mem,memory,"ffile%tab_name_th",        &
!                                                            "file_open")

        DO ith=0,ffile%nb_thread-1

        err_file_loc = err_file_name(ffile%tab_name_th(ith),'file_open')
        IF (.NOT. present(err_file) .AND. err_file_loc /= 0) STOP ' ERROR, the file name is empty!'
        IF (present(err_file)) err_file = err_file_loc

          ffile%tab_name_th(ith) = trim(adjustl(ffile%name)) //         &
                                     "." // int_TO_char(ith)

          inquire(FILE=ffile%tab_name_th(ith),NUMBER=iunit,OPENED=unit_opened)
          IF (.NOT. unit_opened) THEN ! the file is not open

            !- find an unused UNIT ---------
            iunit = 9
            DO
              iunit = iunit + 1
              inquire(UNIT=iunit,OPENED=unit_opened)
              !write(out_unitp,*) 'name,iunit,unit_opened ',ffile%name,iunit,unit_opened
              IF (.NOT. unit_opened) EXIT
            END DO

            !-- open the file
            IF (ffile%seq) THEN
              open(UNIT=iunit,FILE=ffile%tab_name_th(ith),FORM=fform,   &
                STATUS=fstatus,POSITION=fposition,ACCESS='SEQUENTIAL')
            ELSE
              open(UNIT=iunit,FILE=ffile%tab_name_th(ith),FORM=fform,   &
                STATUS=fstatus,ACCESS='DIRECT',RECL=ffile%frecl)
            END IF

            ffile%tab_unit(ith) = iunit
            !write(out_unitp,*) 'open ',ffile%name,iunit
          ELSE
            ffile%tab_unit(ith) = iunit
          END IF
          !write(out_unitp,*) 'open ',ffile%tab_unit(ith),ffile%tab_name_th(ith)

        END DO
      END IF

      iunit = ffile%unit
      !write(out_unitp,*) 'END file_open'

      END SUBROUTINE file_open

      SUBROUTINE file_close(ffile)

      TYPE(param_file)  :: ffile

      integer                   :: ith
      logical                   :: op

      inquire(unit=ffile%unit,OPENED=op)
      IF (op) close(ffile%unit)

      IF (ffile%nb_thread > 1) THEN
        DO ith=0,ffile%nb_thread-1
          inquire(unit=ffile%tab_unit(ith),OPENED=op)
          IF (op)  close(ffile%tab_unit(ith))
        END DO
      END IF
      ffile%init = .FALSE.

      END SUBROUTINE file_close

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE file_open2(name_file,iunit,lformatted,append,old,err_file)

      character (len=*) :: name_file
      integer           :: iunit
      logical, optional :: lformatted,append,old
      integer, optional :: err_file

      character (len=Name_len) :: fform,fstatus,fposition

      logical           :: unit_opened
      integer           :: err_file_loc

!     - default for the open ---------------------------

!     - test if optional arguments are present ---------
      IF (present(lformatted)) THEN
        IF (.NOT. lformatted) THEN
          fform = 'unformatted'
        ELSE
          fform = 'formatted'
        END IF
      ELSE
        fform = 'formatted'
      END IF

      IF (present(append)) THEN
        IF (append) THEN
          fposition = 'append'
        ELSE
          fposition = 'asis'
        END IF
      ELSE
        fposition = 'asis'
      END IF

      IF (present(old)) THEN
        IF (old) THEN
          fstatus = 'old'
        ELSE
          fstatus = 'unknown'
        END IF
      ELSE
          fstatus = 'unknown'
      END IF

      err_file_loc = err_file_name(name_file,'file_open2')
      IF (.NOT. present(err_file) .AND. err_file_loc /= 0) STOP ' ERROR, the file name is empty!'
      IF (present(err_file)) err_file = err_file_loc

!     - check if the file is already open ------------------
      !write(out_unitp,*) 'name_file,iunit ',name_file,iunit ; flush(out_unitp)

      inquire(FILE=name_file,NUMBER=iunit,OPENED=unit_opened)
!     write(out_unitp,*) 'name,unit,unit_opened ',name_file,unit,unit_opened


!     - the file is not open, find an unused UNIT ---------
      IF (unit_opened) RETURN ! the file is already open

      iunit = 9
      DO
        iunit = iunit + 1
        inquire(UNIT=iunit,OPENED=unit_opened)
!       write(out_unitp,*) 'name,iunit,unit_opened ',name_file,iunit,unit_opened
        IF (.NOT. unit_opened) exit
      END DO


!     -- open the file
      IF (present(err_file)) THEN
        open(UNIT=iunit,FILE=name_file,FORM=fform,STATUS=fstatus,       &
             POSITION=fposition,ACCESS='SEQUENTIAL',IOSTAT=err_file)
      ELSE
        open(UNIT=iunit,FILE=name_file,FORM=fform,STATUS=fstatus,       &
             POSITION=fposition,ACCESS='SEQUENTIAL')
      END IF

!     write(out_unitp,*) 'open ',name_file,iunit

      END SUBROUTINE file_open2

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE file_delete(ffile)

      TYPE(param_file)  :: ffile
      integer           :: unit

      integer                  :: ithread,nio
      !integer :: err_mem,memory

      !write(out_unitp,*) 'BEGINNING file_delete'

      IF (len_trim(ffile%name) == 0) RETURN
      CALL file_open(ffile,unit)

      close(unit,status='delete')
      !write(out_unitp,*) 'delete file: ',unit,file%name

      IF (ffile%nb_thread > 1) THEN
        DO ithread=0,ffile%nb_thread-1
          nio = ffile%tab_unit(ithread)
          close(nio,status='delete')
          !write(out_unitp,*) 'delete file: ',nio,ffile%tab_name_th(ithread)
        END DO
      END IF
      ffile%init = .FALSE.
!      memory = ffile%nb_thread
!      deallocate(ffile%tab_unit,stat=err_mem) ! change dealloc done
!      CALL error_memo_allo(err_mem,-memory,"ffile%tab_unit","file_delete")
!
!      memory = ffile%nb_thread
!      deallocate(ffile%tab_name_th,stat=err_mem) ! change dealloc done
!      CALL error_memo_allo(err_mem,-memory,"ffile%tab_name_th","file_delete")

      !write(out_unitp,*) 'END file_delete'


      END SUBROUTINE file_delete

      SUBROUTINE file_dealloc(ffile)

      TYPE(param_file)  :: ffile

      !write(out_unitp,*) 'BEGINNING file_dealloc'

      ! first close the file
      CALL file_close(ffile)

      ffile%init = .FALSE.

      ffile%nb_thread = 0
      IF (allocated(ffile%tab_unit))    deallocate(ffile%tab_unit)
      IF (allocated(ffile%tab_name_th)) deallocate(ffile%tab_name_th)

      !write(out_unitp,*) 'END file_dealloc'

      END SUBROUTINE file_dealloc

      SUBROUTINE file_Write(ffile)

      TYPE(param_file)  :: ffile

      integer :: ith

      write(out_unitp,*) 'BEGINNING file_Write'
      write(out_unitp,*) 'name:    ',trim(adjustl(ffile%name))
      write(out_unitp,*) 'unit     ',ffile%unit
      write(out_unitp,*) 'formatted',ffile%formatted
      write(out_unitp,*) 'old      ',ffile%old
      write(out_unitp,*) 'seq      ',ffile%seq
      write(out_unitp,*) 'frecl    ',ffile%frecl
      write(out_unitp,*) 'init     ',ffile%init
      write(out_unitp,*) 'nb_thread',ffile%nb_thread

      IF (allocated(ffile%tab_name_th) .AND. allocated(ffile%tab_unit)) THEN
        write(out_unitp,*) 'tab_name_th,tab_unit'
        DO ith=0,ffile%nb_thread-1
          write(out_unitp,*) 'name,unit: ',                             &
                             trim(adjustl(ffile%tab_name_th(ith))),' ', &
                             ffile%tab_unit(ith)
        END DO
      END IF

      write(out_unitp,*) 'END file_Write'


      END SUBROUTINE file_Write


      FUNCTION make_FileName(FileName)
        character(len=*), intent(in)    :: FileName

        character (len=:), allocatable  :: make_FileName
        integer :: ilast_char

        ilast_char = len_trim(File_path)

        IF (FileName(1:1) == "/" .OR. FileName(1:1) == "" .OR. ilast_char == 0) THEN
          make_FileName = trim(adjustl(FileName))
        ELSE
          IF (File_path(ilast_char:ilast_char) == "/") THEN
            make_FileName = trim(adjustl(File_path)) // trim(adjustl(FileName))
          ELSE
            make_FileName = trim(adjustl(File_path)) // '/' // trim(adjustl(FileName))
          END IF
        END IF
!write(out_unitp66,*) 'make_FileName: ',make_FileName
!stop
      END FUNCTION make_FileName

  !!@description: TODO
  !!@param: TODO
  SUBROUTINE flush_perso(nio)

  integer, intent(in) :: nio

    flush(nio)

  END  SUBROUTINE flush_perso

      !! @description: Join two path1 with path2. Return "path1/path2"
      !!               If path2 is absolute (starting with /), return path2
      !! @param: path1 First path
      !! @param: path2 Second path
      character(len=Name_longlen) function join_path(path1, path2)

        character(len=*), intent(in) :: path1
        character(len=*), intent(in) :: path2

        character (len=*), parameter :: routine_name = 'join_path'

        if (path2(1:1) == '/') then
          join_path = path2
          return
        end if
        if (path1(len_trim(path1):len_trim(path1)) == '/') then
          join_path = trim(path1)//trim(path2)
        else
          join_path = trim(path1)//"/"//trim(path2)
        end if

      end function join_path

!=======================================================================================
!@brief set save_WP=.true. with external control
!
! create a file with name "Davidson_exit", the program wil exit
! and change the file name to "done_Davidson_exit".
!=======================================================================================
  SUBROUTINE exit_Davidson_external(exit_Davidson,save_WP,it)
    IMPLICIT NONE

    Logical,                       intent(inout) :: exit_Davidson
    Logical,                       intent(inout) :: save_WP
    Integer,                       intent(in)    :: it
    Logical                                      :: exist
    Integer                                      :: stat

    IF(it>2) THEN
      INQUIRE(FILE='Davidson_exit',EXIST=exist)
      IF(exist) THEN
        save_WP=.TRUE.
        exit_Davidson=.TRUE.
        !CALL RENAME('Davidson_exit','done_Davidson_exit')
        open(unit=111, FILE='Davidson_exit')
        close(111, status='delete')
        open(unit=112, FILE='done_Davidson_exit')
        close(112)
      ENDIF
    ENDIF
  END SUBROUTINE exit_Davidson_external
!=======================================================================================
END MODULE mod_file
