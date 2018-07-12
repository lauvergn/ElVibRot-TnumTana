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
MODULE mod_string
  USE mod_NumParameters
  USE mod_memory
  IMPLICIT NONE

  INTERFACE alloc_array
    MODULE PROCEDURE alloc_array_OF_ChLendim1
  END INTERFACE
  INTERFACE dealloc_array
    MODULE PROCEDURE dealloc_array_OF_ChLendim1
  END INTERFACE


  CONTAINS

  FUNCTION String_TO_String(string,ltrim)
  character(len=:), allocatable     :: String_TO_String
  character(len=*), intent(in)      :: string
  logical, optional,intent(in)      :: ltrim

  logical :: ltrim_loc

  IF (allocated(String_TO_String)) deallocate(String_TO_String)

  IF (present(ltrim)) THEN
    ltrim_loc = ltrim
  ELSE
    ltrim_loc = .TRUE.
  END IF

  IF (ltrim_loc) THEN
    allocate(character(len=len_trim(string)) :: String_TO_String)
    String_TO_String = trim(string)
  ELSE
    allocate(character(len=len(string)) :: String_TO_String)
    String_TO_String = string
  END IF

  END FUNCTION String_TO_String


  !! @description: Write an interger type in the character typ
  !! @param:       i       The interger varaiable, input
  !! @param:       ci      The character variable, output
  SUBROUTINE Write_int_IN_char(i, ci)
  integer,   intent(in)                :: i
  character(len=*), intent(inout)      :: ci

  character (len=:), allocatable  :: temp

    temp = int_TO_char(i)

    IF (len(ci) < len(temp)) THEN
      write(out_unitp,*) ' ERROR in Write_int_IN_char'
      write(out_unitp,*) ' len of ci is too small !!'
      write(out_unitp,*) '   Check the fortran !!'
      write(out_unitp,*) '   USE int_TO_char(i) instead of Write_int_IN_char'
      STOP
    END IF

    ci = temp

    deallocate(temp)

  END SUBROUTINE Write_int_IN_char

  SUBROUTINE Write_real_IN_char(R, cR)
  real (kind=Rkind),   intent(in)      :: R
  character(len=*), intent(inout)      :: cR

    ! modif DML 9/11/2012
    write(cR,'(' // RMatIO_format // ')') R
    cR = trim(adjustl(cR))

  END SUBROUTINE Write_real_IN_char

  !!@description: TODO
  !!@param: TODO
  SUBROUTINE make_nameQ(nameQ,baseQ,iQ,it)
    character(len=Name_len), intent(inout) :: nameQ
    character(len=*), intent(in)  :: baseQ
    integer, intent(in), optional :: it,iq

    character(len=Name_len) :: baseW
    integer :: ic

    baseW = trim(adjustl(baseQ))
    DO ic=1,len_trim(baseQ)
      IF (baseW(ic:ic) == " ") baseW(ic:ic) = "_"
    END DO


    nameQ = trim(adjustl(baseW))

    IF (present(it)) THEN
      nameQ = trim(adjustl(nameQ)) // int_TO_char(it)
    END IF

    IF (present(iq)) THEN
      nameQ = trim(adjustl(nameQ)) // "_" // int_TO_char(iq)
    END IF


    !write(out_unitp,*) 'nameQ...: ',nameQ,nameW1,nameW2
  END SUBROUTINE make_nameQ

  FUNCTION logical_TO_char(l)
    character (len=1)  :: logical_TO_char
    logical, intent(in) :: l

    IF (l) THEN
      logical_TO_char = 'T'
    ELSE
      logical_TO_char = 'F'
    END IF

  END FUNCTION logical_TO_char

  FUNCTION int_TO_char(i)
    character (len=:), allocatable  :: int_TO_char
    integer, intent(in)             :: i

    character (len=:), allocatable  :: name_int
    integer :: clen

    IF (allocated(int_TO_char)) deallocate(int_TO_char)

    ! first approximated size of name_int
    IF (i == 0) THEN
      clen = 1
    ELSE IF (i < 0) THEN
      clen = int(log10(abs(real(i,kind=Rkind))))+2
    ELSE
      clen = int(log10(real(i,kind=Rkind)))+1
    END IF

    ! allocate name_int
    allocate(character(len=clen) :: name_int)

    ! write i in name_int
    write(name_int,'(i0)') i

    ! transfert name_int in int_TO_char
    int_TO_char = String_TO_String(name_int)


    ! deallocate name_int
    deallocate(name_int)

  END FUNCTION int_TO_char
  FUNCTION real_TO_char(r,Rformat) RESULT(string)
    character (len=:), allocatable           :: string
    real (kind=Rkind), intent(in)            :: r
    character (len=*), intent(in), optional  :: Rformat

    character(len=Line_len) :: name_real
    integer :: clen,i

    IF (allocated(string)) deallocate(string)


    IF (present(Rformat)) THEN
      write(name_real,'(' // Rformat // ')') r
    ELSE
      write(name_real,*) r
    END IF

    clen = len_trim(adjustl(name_real))
    allocate(character(len=clen) :: string)

    string = trim(adjustl(name_real))

    DO i=len(string),2,-1
      IF (string(i:i) == '0') THEN
        string(i:i) = ' '
      ELSE
        EXIT
      END IF
    END DO

    string = String_TO_String(string)

  END FUNCTION real_TO_char

  FUNCTION nom_i(nom1,i1)
  IMPLICIT NONE

  character (len=14) :: nom_i
  character (len=10) :: nom1
  character (len=14) :: nom2
  integer            :: j,i1

  write(out_unitp,*) nom1,i1
  IF (i1 .GT. 100 ) STOP ' in nom_i: i1 too big'

  write(nom2,'(a10,i2)') nom1,i1
  DO j=1,12  ! it has to be 12 and not 14
    IF (nom2(j:j) .EQ. ' ') nom2(j:j)='_'
  END DO
  nom_i=nom2

  END FUNCTION nom_i

  FUNCTION nom_ii(nom1,i1,i2)
  IMPLICIT NONE

  character (len=14) :: nom_ii

  character (len=10) :: nom1
  character (len=14) :: nom2
  integer            :: j,i1,i2

  !write(out_unitp,*) nom1,i1,i2
  IF (i1 .GT. 100 .OR. i2 .GT. 100) STOP ' in nom_ii: i1 or i2 too big'

  write(nom2,'(a10,2i2)') nom1,i1,i2
  DO j=1,14
    IF (nom2(j:j) .EQ. ' ') nom2(j:j)='_'
  END DO
  nom_ii=nom2

  END FUNCTION nom_ii

  SUBROUTINE read_name_advNo(nio,Read_name,err_io)
    character(len=*), intent(inout) :: Read_name
    integer, intent(inout) :: err_io
    integer, intent(in) :: nio

    character(len=1) :: chara
    logical          :: first
    integer :: ic

    Read_name = ''
    first     = .TRUE.
    ic        = 0
    err_io    = 0
    DO
      read(nio,'(a1)',IOSTAT=err_io,advance='no') chara

      IF (err_io /= 0)   EXIT
      !write(6,*) 'ic,chara',ic,'"',chara,'"'
      IF (chara == ' ' .AND. .NOT. first) EXIT

      IF (chara == ' ' .AND. first) CYCLE

      ic = ic + 1
      Read_name(ic:ic) = chara
      first = .FALSE.

    END DO
    !write(6,*) 'Read_name: ',trim(Read_name)

  END SUBROUTINE read_name_advNo

  !!@description: Change the case of a string (default lowercase)
  !!@param: name_string the string
  !!@param: lower If the variable is present and its value is F,
  !!              the string will be converted into a uppercase string, otherwise,
  !!              it will be convert into a lowercase string.
  SUBROUTINE string_uppercase_TO_lowercase(name_string,lower)

   character (len=*), intent(inout)  :: name_string
   logical, optional  :: lower

   logical  :: lower_loc
   integer  :: i,ascii_char

   IF (present(lower)) THEN
     lower_loc = lower
   ELSE
     lower_loc = .TRUE.
   END IF

   !write(out_unitp,*) 'name_string: ',name_string
   IF (lower_loc) THEN ! uppercase => lowercase
     DO i=1,len_trim(name_string)
       ascii_char = iachar(name_string(i:i))
       IF (ascii_char >= 65 .AND. ascii_char <= 90)                 &
                           name_string(i:i) = achar(ascii_char+32)
     END DO
   ELSE ! lowercase => uppercase
     DO i=1,len_trim(name_string)
       ascii_char = iachar(name_string(i:i))
       IF (ascii_char >= 97 .AND. ascii_char <= 122)                 &
                            name_string(i:i) = achar(ascii_char-32)

     END DO
   END IF
   !write(out_unitp,*) 'name_string: ',name_string


  END SUBROUTINE string_uppercase_TO_lowercase

  SUBROUTINE alloc_array_OF_ChLendim1(tab,tab_ub,ChLen,name_var,name_sub,tab_lb)
  IMPLICIT NONE

  integer, intent(in) :: ChLen
  character (len=ChLen), pointer, intent(out) :: tab(:)
  integer, intent(in) :: tab_ub(:)
  integer, intent(in), optional :: tab_lb(:)

  character (len=*), intent(in) :: name_var,name_sub

  integer, parameter :: ndim=1

  !----- for debuging --------------------------------------------------
  character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_ChLendim1'
  integer :: err_mem,memory
  logical,parameter :: debug=.FALSE.
  !logical,parameter :: debug=.TRUE.
  !----- for debuging --------------------------------------------------


   IF (associated(tab))                                             &
         CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

   CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

   IF (present(tab_lb)) THEN
     CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)

     memory = ChLen * product(tab_ub(:)-tab_lb(:)+1)
     allocate(tab(tab_lb(1):tab_ub(1)),stat=err_mem)
   ELSE
     memory = ChLen * product(tab_ub(:))
     allocate(tab(tab_ub(1)),stat=err_mem)
   END IF
   CALL error_memo_allo(err_mem,memory,name_var,name_sub,'character')

  END SUBROUTINE alloc_array_OF_ChLendim1
  SUBROUTINE dealloc_array_OF_ChLendim1(tab,name_var,name_sub)
  IMPLICIT NONE

  character (len=*), pointer, intent(out) :: tab(:)
  character (len=*), intent(in) :: name_var,name_sub

  !----- for debuging --------------------------------------------------
  character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_ChLendim1'
  integer :: err_mem,memory
  logical,parameter :: debug=.FALSE.
  !logical,parameter :: debug=.TRUE.
  !----- for debuging --------------------------------------------------

   !IF (.NOT. associated(tab)) RETURN
   IF (.NOT. associated(tab))                                       &
         CALL Write_error_null(name_sub_alloc,name_var,name_sub)

   memory = size(tab) * len(tab(lbound(tab,dim=1)))
   deallocate(tab,stat=err_mem)
   CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'character')
   nullify(tab)

  END SUBROUTINE dealloc_array_OF_ChLendim1

END MODULE mod_string

