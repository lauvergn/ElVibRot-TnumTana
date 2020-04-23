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
  MODULE mod_RealWithUnit
  USE mod_NumParameters
  USE mod_string
  IMPLICIT NONE

  PRIVATE

    TYPE REAL_WU ! real with unit

      real (kind=Rkind)        :: val      = ZERO
      character (len=Name_len) :: unit     = ''
      character (len=Name_len) :: quantity = ''
    CONTAINS
      PROCEDURE, PRIVATE, PASS(RWU2) :: RWU2_TO_R1
      PROCEDURE, PRIVATE, PASS(RWU1) :: RWU2_TO_RWU1
      GENERIC,   PUBLIC  :: assignment(=) => RWU2_TO_R1,RWU2_TO_RWU1
    END TYPE REAL_WU

    TYPE Type_TabConvRWU ! real with unit
      character (len=Name_len)   :: quantity = ''
      TYPE(REAL_WU), allocatable :: conv(:)
      TYPE(REAL_WU)              :: Work_unit
      TYPE(REAL_WU)              :: Write_unit
    CONTAINS
      PROCEDURE, PRIVATE, PASS(TabConvRWU1) :: TabConvRWU2_TO_TabConvRWU1
      GENERIC,   PUBLIC  :: assignment(=) => TabConvRWU2_TO_TabConvRWU1
    END TYPE Type_TabConvRWU

    TYPE(Type_TabConvRWU), allocatable, public :: Tab_conv_FOR_quantity(:) ! conversion factor to the working unit


    PUBLIC :: REAL_WU, Type_TabConvRWU
    PUBLIC :: dealloc_TabConvRWU, dealloc_TabConvRWU_dim1, Write_TabConvRWU, Write_TabConvRWU_dim1
    PUBLIC :: ADD_RWU_TO_TabConvRWU, ADD_RWU_TO_Tab_conv_FOR_quantity
    PUBLIC :: convRWU_TO_R,convRWU_TO_R_WITH_WorkingUnit,convRWU_TO_R_WITH_WritingUnit,convRWU_TO_RWU
    PUBLIC :: RWU_Write,RWU_WriteUnit
    PUBLIC :: get_Conv_au_TO_WriteUnit,get_Conv_au_TO_Unit,get_val_FROM_RWU
    PUBLIC :: Test_RWU


  CONTAINS
  ELEMENTAL SUBROUTINE RWU2_TO_R1(R1,RWU2)
    CLASS(REAL_WU),    intent(in)      :: RWU2
    real (kind=Rkind), intent(inout)   :: R1

    R1 = RWU2%val

  END SUBROUTINE RWU2_TO_R1
  ELEMENTAL FUNCTION get_val_FROM_RWU(RWU) RESULT(val)
    TYPE(REAL_WU),     intent(in)      :: RWU
    real (kind=Rkind)                  :: val

    val = RWU%val

  END FUNCTION get_val_FROM_RWU
  ELEMENTAL SUBROUTINE RWU2_TO_RWU1(RWU1,RWU2)
    TYPE(REAL_WU),  intent(in)      :: RWU2
    CLASS(REAL_WU), intent(inout)   :: RWU1

    RWU1%val      = RWU2%val
    RWU1%unit     = RWU2%unit
    RWU1%quantity = RWU2%quantity

  END SUBROUTINE RWU2_TO_RWU1

  ELEMENTAL SUBROUTINE dealloc_TabConvRWU(TabConvRWU)
    TYPE(Type_TabConvRWU), intent(inout)   :: TabConvRWU

    TabConvRWU%quantity   = ''
    TabConvRWU%Work_unit  = REAL_WU(ZERO,'','')
    TabConvRWU%Write_unit = REAL_WU(ZERO,'','')

    IF (allocated(TabConvRWU%conv)) THEN
      deallocate(TabConvRWU%conv)
    END IF

  END SUBROUTINE dealloc_TabConvRWU

  ELEMENTAL SUBROUTINE TabConvRWU2_TO_TabConvRWU1(TabConvRWU1,TabConvRWU2)
    TYPE(Type_TabConvRWU),  intent(in)      :: TabConvRWU2
    CLASS(Type_TabConvRWU), intent(inout)   :: TabConvRWU1

    integer :: i

    TabConvRWU1%quantity   = TabConvRWU2%quantity
    TabConvRWU1%Work_unit  = TabConvRWU2%Work_unit
    TabConvRWU1%Write_unit = TabConvRWU2%Write_unit

    IF (allocated(TabConvRWU2%conv)) THEN
      allocate(TabConvRWU1%conv(size(TabConvRWU2%conv)))
      DO i=1,size(TabConvRWU2%conv)
        TabConvRWU1%conv(i) = TabConvRWU2%conv(i)
      END DO
    END IF

  END SUBROUTINE TabConvRWU2_TO_TabConvRWU1

  SUBROUTINE dealloc_TabConvRWU_dim1(TabConvRWU)
    TYPE(Type_TabConvRWU), intent(inout), allocatable   :: TabConvRWU(:)

    integer :: i

    IF (allocated(TabConvRWU)) THEN
      DO i=1,size(TabConvRWU)
        CALL dealloc_TabConvRWU(TabConvRWU(i))
      END DO
      deallocate(TabConvRWU)
    END IF

  END SUBROUTINE dealloc_TabConvRWU_dim1

  SUBROUTINE TabConvRWU2_TO_TabConvRWU1_dim1(TabConvRWU1,TabConvRWU2)
    TYPE(Type_TabConvRWU), intent(in),    allocatable :: TabConvRWU2(:)
    TYPE(Type_TabConvRWU), intent(inout), allocatable :: TabConvRWU1(:)

    integer :: i

    IF (allocated(TabConvRWU1)) CALL dealloc_TabConvRWU_dim1(TabConvRWU1)

    IF (allocated(TabConvRWU2)) THEN
      allocate(TabConvRWU1(size(TabConvRWU2)))
      DO i=1,size(TabConvRWU2)
        TabConvRWU1(i) = TabConvRWU2(i)
      END DO
    END IF

  END SUBROUTINE TabConvRWU2_TO_TabConvRWU1_dim1


  SUBROUTINE Write_TabConvRWU(TabConvRWU)
    TYPE(Type_TabConvRWU), intent(in)      :: TabConvRWU

    integer :: i

    write(out_unitp,*) "======================================"
    write(out_unitp,*) "quantity:   ",TabConvRWU%quantity
    write(out_unitp,*) "Work_unit:  ",TabConvRWU%Work_unit
    write(out_unitp,*) "Write_unit: ",TabConvRWU%Write_unit
    IF (allocated(TabConvRWU%conv)) THEN
      DO i=1,size(TabConvRWU%conv)
        write(out_unitp,*) "conv(i):  ",i,TabConvRWU%conv(i)
      END DO
    END IF
    write(out_unitp,*) "======================================"

  END SUBROUTINE Write_TabConvRWU

  SUBROUTINE Write_TabConvRWU_dim1(TabConvRWU)
    TYPE(Type_TabConvRWU), intent(in), allocatable  :: TabConvRWU(:)

    integer :: i

    IF (allocated(TabConvRWU)) THEN
      DO i=1,size(TabConvRWU)
        CALL Write_TabConvRWU(TabConvRWU(i))
      END DO
    END IF

  END SUBROUTINE Write_TabConvRWU_dim1

  SUBROUTINE ADD_RWU_TO_TabConvRWU(TabConvRWU,RWU,Work_unit,Write_unit)
    USE mod_MPI

    TYPE(REAL_WU), intent(in)           :: RWU
    TYPE(Type_TabConvRWU), intent(inout)   :: TabConvRWU
    logical, optional :: Work_unit,Write_unit

    integer               :: i
    TYPE(Type_TabConvRWU) :: TabConvRWU_loc
    logical               :: unit_present

    character (len=:), allocatable  :: name_quantity,name_RWUquantity
    character (len=:), allocatable  :: name_unit,name_RWUunit


    IF (TabConvRWU%quantity == '') TabConvRWU%quantity = RWU%quantity

    allocate(character(len=len(RWU%quantity)) :: name_RWUquantity)
    name_RWUquantity = RWU%quantity
    CALL string_uppercase_TO_lowercase(name_RWUquantity)

    allocate(character(len=len(TabConvRWU%quantity)) :: name_quantity)
    name_quantity = TabConvRWU%quantity
    CALL string_uppercase_TO_lowercase(name_quantity)

    IF (name_quantity /= name_RWUquantity) THEN
      write(out_unitp,*) ' ERROR in ADD_RWU_TO_TabConvRWU'
      write(out_unitp,*) ' quantities of RWU and TabConvRWU are different!'
      write(out_unitp,*) ' RWU%quantity:        ',RWU%quantity
      write(out_unitp,*) ' TabConvRWU%quantity: ',TabConvRWU%quantity
      STOP
    END IF

    IF (allocated(name_quantity))    deallocate(name_quantity)
    IF (allocated(name_RWUquantity)) deallocate(name_RWUquantity)

    ! check if the unit is already present
    unit_present = .FALSE.
    IF (allocated(TabConvRWU%conv)) THEN

      allocate(character(len=len(RWU%unit)) :: name_RWUunit)
      name_RWUunit = RWU%unit
      CALL string_uppercase_TO_lowercase(name_RWUunit)

      DO i=1,size(TabConvRWU%conv)

        allocate(character(len=len(TabConvRWU%conv(i)%unit)) :: name_unit)
        name_unit = TabConvRWU%conv(i)%unit
        CALL string_uppercase_TO_lowercase(name_unit)

        unit_present = (name_unit == name_RWUunit)

        deallocate(name_unit)

        IF (unit_present) THEN
          IF(MPI_id==0) write(out_unitp,*) 'The unit "',trim(adjustl(RWU%unit)),      &
                                           '" is already present'
          EXIT
        END IF

      END DO

      deallocate(name_RWUunit)

    ELSE
      unit_present = .FALSE.
    END IF

    IF (.NOT. unit_present) THEN
      TabConvRWU_loc = TabConvRWU

      IF (allocated(TabConvRWU%conv)) deallocate(TabConvRWU%conv)

      IF (allocated(TabConvRWU_loc%conv)) THEN

        allocate(TabConvRWU%conv(size(TabConvRWU_loc%conv)+1))

        DO i=1,size(TabConvRWU_loc%conv)
          TabConvRWU%conv(i) = TabConvRWU_loc%conv(i)
        END DO

        TabConvRWU%conv(size(TabConvRWU_loc%conv)+1) = RWU

        deallocate(TabConvRWU_loc%conv)
      ELSE
        allocate(TabConvRWU%conv(1))
        TabConvRWU%conv(1) = RWU
      END IF
    END IF

    IF (present(Write_unit)) THEN
    IF (Write_unit) THEN
      TabConvRWU%Write_unit = RWU
    END IF
    END IF

    IF (present(Work_unit)) THEN
    IF (Work_unit) THEN
      TabConvRWU%Work_unit = RWU
    END IF
    END IF

  END SUBROUTINE ADD_RWU_TO_TabConvRWU

  SUBROUTINE ADD_RWU_TO_Tab_conv_FOR_quantity(RWU,Work_unit,Write_unit)
   TYPE(REAL_WU), intent(in)           :: RWU
   logical, optional :: Work_unit,Write_unit

   integer            :: i,iq
   TYPE(Type_TabConvRWU), allocatable :: Tab_conv_FOR_quantity_loc(:)
   logical :: Work_unit_loc,Write_unit_loc

   character (len=:), allocatable  :: name_quantity,name_RWUquantity


   IF (present(Write_unit)) THEN
     Write_unit_loc = Write_unit
   ELSE
     Write_unit_loc = .FALSE.
   END IF

   IF (present(Work_unit)) THEN
     Work_unit_loc = Work_unit
   ELSE
     Work_unit_loc = .FALSE.
   END IF

   allocate(character(len=len(RWU%quantity)) :: name_RWUquantity)
   name_RWUquantity = RWU%quantity
   CALL string_uppercase_TO_lowercase(name_RWUquantity)

   IF (allocated(Tab_conv_FOR_quantity)) THEN
     ! first find the index of RWU%quantity
     DO i=1,size(Tab_conv_FOR_quantity)
       allocate(character(len=len(Tab_conv_FOR_quantity(i)%quantity)) :: name_quantity)
       name_quantity = Tab_conv_FOR_quantity(i)%quantity
       CALL string_uppercase_TO_lowercase(name_quantity)
       IF (name_quantity == name_RWUquantity) EXIT
       deallocate(name_quantity)
     END DO
     iq = i
     IF (iq > size(Tab_conv_FOR_quantity)) THEN
       !!!--------------------------------------------------
       !!!! Compiler ERROR with ifort (ifort version 16.0.3) with the "= assignment"
       ! Tab_conv_FOR_quantity_loc=Tab_conv_FOR_quantity
       !!!  We have to use the corresponding subroutine
       CALL TabConvRWU2_TO_TabConvRWU1_dim1(Tab_conv_FOR_quantity_loc,&
                                            Tab_conv_FOR_quantity)
       !!!--------------------------------------------------
       CALL dealloc_TabConvRWU_dim1(Tab_conv_FOR_quantity)
       allocate(Tab_conv_FOR_quantity(size(Tab_conv_FOR_quantity_loc)+1))
       DO i=1,size(Tab_conv_FOR_quantity_loc)
         Tab_conv_FOR_quantity(i) = Tab_conv_FOR_quantity_loc(i)
       END DO
       i = size(Tab_conv_FOR_quantity_loc)+1
       CALL ADD_RWU_TO_TabConvRWU(Tab_conv_FOR_quantity(i),RWU,      &
                   Work_unit=Work_unit_loc,Write_unit=Write_unit_loc)
     ELSE
       CALL ADD_RWU_TO_TabConvRWU(Tab_conv_FOR_quantity(iq),RWU,     &
                   Work_unit=Work_unit_loc,Write_unit=Write_unit_loc)
     END IF
   ELSE
     allocate(Tab_conv_FOR_quantity(1))
     CALL ADD_RWU_TO_TabConvRWU(Tab_conv_FOR_quantity(1),RWU,        &
                   Work_unit=Work_unit_loc,Write_unit=Write_unit_loc)
   END IF

   IF (allocated(name_quantity))    deallocate(name_quantity)
   IF (allocated(name_RWUquantity)) deallocate(name_RWUquantity)

  END SUBROUTINE ADD_RWU_TO_Tab_conv_FOR_quantity

  FUNCTION convRWU_TO_R(RWU,WorkingUnit)
    real (kind=Rkind)                  :: convRWU_TO_R
    TYPE(REAL_WU),          intent(in) :: RWU
    logical,      optional, intent(in) :: WorkingUnit

    TYPE(REAL_WU)                 :: RWU_loc


    IF (present(WorkingUnit)) THEN
      RWU_loc = convRWU_TO_RWU(RWU,WorkingUnit=WorkingUnit)
    ELSE
      RWU_loc = convRWU_TO_RWU(RWU,WorkingUnit=.TRUE.)
    END IF

    convRWU_TO_R = RWU_loc%val

  END FUNCTION convRWU_TO_R

  FUNCTION convRWU_TO_R_WITH_WorkingUnit(RWU)
    real (kind=Rkind)             :: convRWU_TO_R_WITH_WorkingUnit
    TYPE(REAL_WU), intent(in)     :: RWU

    TYPE(REAL_WU)                 :: RWU_loc

    RWU_loc = convRWU_TO_RWU(RWU,WorkingUnit=.TRUE.)

    convRWU_TO_R_WITH_WorkingUnit = RWU_loc%val

  END FUNCTION convRWU_TO_R_WITH_WorkingUnit
  FUNCTION convRWU_TO_R_WITH_WritingUnit(RWU)
    real (kind=Rkind)             :: convRWU_TO_R_WITH_WritingUnit
    TYPE(REAL_WU), intent(in)     :: RWU

    TYPE(REAL_WU)                 :: RWU_loc

    RWU_loc = convRWU_TO_RWU(RWU,WorkingUnit=.FALSE.)

    convRWU_TO_R_WITH_WritingUnit = RWU_loc%val

    END FUNCTION convRWU_TO_R_WITH_WritingUnit
    FUNCTION convRWU_TO_RWU(RWU,WorkingUnit)
    TYPE(REAL_WU)                       :: convRWU_TO_RWU
    TYPE(REAL_WU), intent(in)           :: RWU
    logical,       intent(in), optional :: WorkingUnit


    logical                             :: skip_conv,WorkingUnit_loc

    real (kind=Rkind)                   :: conv
    integer :: i,iq
    character (len=Name_len)  :: name_quantity,name_RWUquantity
    character (len=Name_len)  :: name_unit,name_RWUunit


    IF (present(WorkingUnit)) THEN
      WorkingUnit_loc = WorkingUnit
    ELSE
      WorkingUnit_loc = .TRUE.
    END IF

    convRWU_TO_RWU     = RWU

    skip_conv = (RWU%val == huge(ONE) .OR. RWU%val == -huge(ONE))

    !write(out_unitp,*) 'RWU',RWU,skip_conv

    ! first find the quantity
    iq = get_Index_OF_Quantity(RWU%quantity)

    name_RWUunit = RWU%unit
    CALL string_uppercase_TO_lowercase(name_RWUunit)
    !write(out_unitp,*) 'name_RWUunit',name_RWUunit

    IF (iq <= size(Tab_conv_FOR_quantity)) THEN
      ! modify quantity to have the correct case
      convRWU_TO_RWU%quantity = Tab_conv_FOR_quantity(iq)%quantity

      conv = ONE
      DO i=1,size(Tab_conv_FOR_quantity(iq)%conv)

        name_unit = Tab_conv_FOR_quantity(iq)%conv(i)%unit
        CALL string_uppercase_TO_lowercase(name_unit)

        !write(out_unitp,*) 'i,name_unit',i,name_unit,name_RWUunit,(name_RWUunit == name_unit)

        IF (name_RWUunit == name_unit) THEN
          IF (WorkingUnit_loc) THEN ! working unit
            conv = Tab_conv_FOR_quantity(iq)%conv(i)%val ! working unit
            convRWU_TO_RWU%unit = Tab_conv_FOR_quantity(iq)%Work_unit%unit
          ELSE ! writing unit
            conv = Tab_conv_FOR_quantity(iq)%conv(i)%val /          &
                            Tab_conv_FOR_quantity(iq)%Write_unit%val ! printing unit
            convRWU_TO_RWU%unit = Tab_conv_FOR_quantity(iq)%Write_unit%unit
          END IF
          EXIT
        END IF
      END DO
      skip_conv = skip_conv .OR. (i > size(Tab_conv_FOR_quantity(iq)%conv))

      IF (.NOT. skip_conv) convRWU_TO_RWU%val  = RWU%val * conv
    END IF

    !write(out_unitp,*) 'conv',conv,skip_conv

  END FUNCTION convRWU_TO_RWU

  FUNCTION RWU_Write(RWU,WithUnit,WorkingUnit)
    character (len=:), allocatable  :: RWU_Write
    TYPE(REAL_WU), intent(in)       :: RWU
    logical, intent(in)             :: WithUnit
    logical, intent(in)             :: WorkingUnit


    TYPE(REAL_WU)                 :: RWU_loc
    character (len=Name_longlen)  :: RWU_Write_loc
    integer :: clen

    RWU_loc = convRWU_TO_RWU(RWU,WorkingUnit)

    IF (abs(RWU_loc%val) < TEN**5) THEN
      write(RWU_Write_loc,'(f14.6)') RWU_loc%val
    ELSE
      write(RWU_Write_loc,'(e14.6)') RWU_loc%val
    END IF

    IF (WithUnit) THEN
      RWU_Write_loc = trim(adjustl(RWU_Write_loc)) // ' ' // RWU_loc%unit
    ELSE
      RWU_Write_loc = trim(adjustl(RWU_Write_loc))
    END IF
    clen = len_trim(RWU_Write_loc)
    allocate(character(len=clen) :: RWU_Write)
    RWU_Write = trim(adjustl(RWU_Write_loc))

  END function RWU_Write

  FUNCTION RWU_WriteUnit(quantity,WorkingUnit)
    character (len=:), allocatable  :: RWU_WriteUnit

    logical,           intent(in) :: WorkingUnit
    character (len=*), intent(in) :: quantity

    integer :: clen,i

    ! first find the quantity
    i = get_Index_OF_Quantity(quantity)
    IF (i > size(Tab_conv_FOR_quantity)) THEN
      clen = len_trim('quantity "' // trim(quantity) // '" not found!')
      allocate(character(len=clen) :: RWU_WriteUnit)
      RWU_WriteUnit = 'quantity "' // trim(quantity) // '" not found!'
    ELSE

      IF (WorkingUnit) THEN
        clen = len_trim(Tab_conv_FOR_quantity(i)%Work_unit%unit)
        allocate(character(len=clen) :: RWU_WriteUnit)
        RWU_WriteUnit = trim(adjustl(Tab_conv_FOR_quantity(i)%Work_unit%unit))
      ELSE
        clen = len_trim(Tab_conv_FOR_quantity(i)%Write_unit%unit)
        allocate(character(len=clen) :: RWU_WriteUnit)
        RWU_WriteUnit = trim(adjustl(Tab_conv_FOR_quantity(i)%Write_unit%unit))
      END IF
    END IF

  END function RWU_WriteUnit

  FUNCTION get_Conv_au_TO_WriteUnit(quantity,Unit)
    real (kind=Rkind)  :: get_Conv_au_TO_WriteUnit

    character (len=*), intent(in)              :: quantity
    character (len=*), intent(inout), optional :: unit

    integer :: i

    i = get_Index_OF_Quantity(quantity)
    get_Conv_au_TO_WriteUnit = ONE/Tab_conv_FOR_quantity(i)%Write_unit%val

    IF (present(unit)) THEN
      unit = Tab_conv_FOR_quantity(i)%Write_unit%unit
    END IF

  END function get_Conv_au_TO_WriteUnit

  FUNCTION get_Conv_au_TO_unit(quantity,Unit,WorkingUnit,err_unit)
    real (kind=Rkind)  :: get_Conv_au_TO_unit

    character (len=*), intent(in)              :: quantity
    character (len=*), intent(in),    optional :: unit
    logical,           intent(in),    optional :: WorkingUnit
    integer,           intent(inout), optional :: err_unit

    integer :: iq,iu,err_unit_loc

    IF (present(err_unit)) err_unit = 0

    iq = get_Index_OF_Quantity(quantity)
    IF (present(WorkingUnit)) THEN
      ! first find the quantity
      IF (WorkingUnit) THEN
        get_Conv_au_TO_unit = ONE/Tab_conv_FOR_quantity(iq)%Work_unit%val
      ELSE
        get_Conv_au_TO_unit = ONE/Tab_conv_FOR_quantity(iq)%Write_unit%val
      END IF
    ELSE IF (present(unit)) THEN
      iu = get_Index_OF_Unit(unit,iq,err_unit=err_unit_loc)
      IF (err_unit_loc /= 0) THEN
        get_Conv_au_TO_unit = ONE
        IF (present(err_unit)) err_unit = err_unit_loc
      ELSE
        get_Conv_au_TO_unit = ONE/Tab_conv_FOR_quantity(iq)%conv(iu)%val
      END IF

    ELSE ! problem, unit and WorkingUnit are not present (conversion factor = ONE)
      get_Conv_au_TO_unit = ONE
      IF (present(err_unit)) err_unit = -2
    END IF

  END function get_Conv_au_TO_unit

  FUNCTION get_Index_OF_Quantity(quantity)
    integer  :: get_Index_OF_Quantity
    character (len=*), intent(in) :: quantity

    character (len=Name_len)            :: quantity_loc,name_quantity
    integer :: i

    ! first find the quantity
    quantity_loc = quantity
    CALL string_uppercase_TO_lowercase(quantity_loc)

    DO i=1,size(Tab_conv_FOR_quantity)

      name_quantity = Tab_conv_FOR_quantity(i)%quantity
      CALL string_uppercase_TO_lowercase(name_quantity)

      IF (name_quantity == quantity_loc) THEN
        EXIT
      END IF
    END DO

    get_Index_OF_Quantity = i
  END FUNCTION get_Index_OF_Quantity

  FUNCTION get_Index_OF_Unit(unit,iq,err_unit)
    integer  :: get_Index_OF_Unit
    character (len=*), intent(in)              :: unit
    integer,           intent(in)              :: iq
    integer,           intent(inout), optional :: err_unit


    character (len=Name_len)            :: unit_loc,name_unit
    integer :: iu


    IF (present(err_unit)) err_unit = 0
    iu = 0
    ! first find the unit
    unit_loc = unit
    CALL string_uppercase_TO_lowercase(unit_loc)

    IF (iq <= size(Tab_conv_FOR_quantity)) THEN
      DO iu=1,size(Tab_conv_FOR_quantity(iq)%conv)

        name_unit = Tab_conv_FOR_quantity(iq)%conv(iu)%unit
        CALL string_uppercase_TO_lowercase(name_unit)

        !write(out_unitp,*) 'i,name_unit',i,name_unit,name_RWUunit,(name_RWUunit == name_unit)

        IF (unit_loc == name_unit) EXIT
      END DO

      IF (iu > size(Tab_conv_FOR_quantity(iq)%conv)) THEN
        IF (present(err_unit)) THEN
          err_unit = -2
        ELSE
          STOP 'iu is out-of-range (unity cannot be found!!)'
        END IF
      END IF
    ELSE
      IF (present(err_unit)) THEN
        err_unit = -1
      ELSE
        STOP 'iq is out-of-range (quantity cannot be found!!)'
      END IF
    END IF

    get_Index_OF_Unit = iu

  END FUNCTION get_Index_OF_Unit
  SUBROUTINE Test_RWU()

    TYPE(REAL_WU)     :: RWU1,RWU2 ! test the real with unit convertion
    real (kind=Rkind) :: Rval

    integer :: i,iq

    write(out_unitp,*) "======================================"
    write(out_unitp,*) "======================================"
    write(out_unitp,*) "======================================"
    write(out_unitp,*) "======= Set up conversion factor ====="

    CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(ONE,'au','l'),Work_unit=.TRUE.)
    CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(ONE,'AU','l'))
    CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(ONE/0.5291772575_Rkind,'Angs','L'),Write_unit=.TRUE.)

    CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(ONE,'au','E'),Work_unit=.TRUE.)
    CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(ONE,'hartree','E'))
    CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(ONE/27.211396_Rkind,'eV','E'),Write_unit=.TRUE.)

    CALL Write_TabConvRWU_dim1(Tab_conv_FOR_quantity)


    DO iq=1,size(Tab_conv_FOR_quantity)
      write(out_unitp,*) "======================================"
      write(out_unitp,*) "======================================"
      write(out_unitp,*) "======================================"
      write(out_unitp,*) "======= ",Tab_conv_FOR_quantity(iq)%quantity," ========="
      flush(out_unitp)

      DO i=1,size(Tab_conv_FOR_quantity(iq)%conv)

        RWU1 = REAL_WU(ONE,Tab_conv_FOR_quantity(iq)%conv(i)%unit,Tab_conv_FOR_quantity(iq)%quantity)

        write(out_unitp,*) 'test RWU (without conv): ',RWU1
        write(out_unitp,*) 'test RWU (working unit): ',RWU_Write(RWU1,WithUnit=.FALSE.,WorkingUnit=.TRUE.)
        write(out_unitp,*) 'test RWU (working unit): ',RWU_Write(RWU1,WithUnit=.TRUE.,WorkingUnit=.TRUE.)
        write(out_unitp,*) 'test RWU (writing unit): ',RWU_Write(RWU1,WithUnit=.FALSE.,WorkingUnit=.FALSE.)
        write(out_unitp,*) 'test RWU (writing unit): ',RWU_Write(RWU1,WithUnit=.TRUE.,WorkingUnit=.FALSE.)
        RWU2 = convRWU_TO_RWU(RWU1)
        write(out_unitp,*) 'test RWU (after conv, working unit default) :',RWU2
        RWU2 = convRWU_TO_RWU(RWU1,WorkingUnit=.TRUE.)
        write(out_unitp,*) 'test RWU (after conv, working unit)         :',RWU2
        RWU2 = convRWU_TO_RWU(RWU1,WorkingUnit=.FALSE.)
        write(out_unitp,*) 'test RWU (after conv, writing unit)         :',RWU2

        Rval = convRWU_TO_R(RWU1)
        write(out_unitp,*) 'test Rval (after conv, working unit default):',Rval
        Rval = RWU1
        write(out_unitp,*) 'test Rval (after conv, = default):           ',Rval
        Rval = convRWU_TO_R(RWU1,WorkingUnit=.TRUE.)
        write(out_unitp,*) 'test Rval (after conv, working unit)        :',Rval
        Rval = convRWU_TO_R(RWU1,WorkingUnit=.FALSE.)
        write(out_unitp,*) 'test Rval (after conv, writing unit)        :',Rval

        write(out_unitp,*) "======================================"
        flush(out_unitp)
      END DO
      write(out_unitp,*) "======================================"
      write(out_unitp,*) "======================================"
      write(out_unitp,*) "======================================"
      flush(out_unitp)
    END DO

    write(out_unitp,*) "======================================"
    write(out_unitp,*) "======== test with a wrong unit ======"

    RWU1 = REAL_WU(ONE,'xx','E')

    write(out_unitp,*) 'test RWU (without conv): ',RWU1
    write(out_unitp,*) 'test RWU (working unit): ',RWU_Write(RWU1,WithUnit=.FALSE.,WorkingUnit=.TRUE.)
    write(out_unitp,*) 'test RWU (working unit): ',RWU_Write(RWU1,WithUnit=.TRUE.,WorkingUnit=.TRUE.)
    write(out_unitp,*) 'test RWU (writing unit): ',RWU_Write(RWU1,WithUnit=.FALSE.,WorkingUnit=.FALSE.)
    write(out_unitp,*) 'test RWU (writing unit): ',RWU_Write(RWU1,WithUnit=.TRUE.,WorkingUnit=.FALSE.)
    RWU2 = convRWU_TO_RWU(RWU1)
    write(out_unitp,*) 'test RWU (after conv, working unit default) :',RWU2
    RWU2 = convRWU_TO_RWU(RWU1,WorkingUnit=.TRUE.)
    write(out_unitp,*) 'test RWU (after conv, working unit)         :',RWU2
    RWU2 = convRWU_TO_RWU(RWU1,WorkingUnit=.FALSE.)
    write(out_unitp,*) 'test RWU (after conv, writing unit)         :',RWU2

    Rval = convRWU_TO_R(RWU1)
    write(out_unitp,*) 'test Rval (after conv, working unit default):',Rval
    Rval = RWU1
    write(out_unitp,*) 'test Rval (after conv, = default):           ',Rval
    Rval = convRWU_TO_R(RWU1,WorkingUnit=.TRUE.)
    write(out_unitp,*) 'test Rval (after conv, working unit)        :',Rval
    Rval = convRWU_TO_R(RWU1,WorkingUnit=.FALSE.)
    write(out_unitp,*) 'test Rval (after conv, writing unit)        :',Rval

    write(out_unitp,*) "======================================"
    flush(out_unitp)

    write(out_unitp,*) "======================================"
    write(out_unitp,*) "======================================"
    write(out_unitp,*) "======================================"

    write(out_unitp,*) "======================================"
    write(out_unitp,*) "======== test with a wrong quantity =="

    RWU1 = REAL_WU(ONE,'au','x')

    write(out_unitp,*) 'test RWU (without conv): ',RWU1
    write(out_unitp,*) 'test RWU (working unit): ',RWU_Write(RWU1,WithUnit=.FALSE.,WorkingUnit=.TRUE.)
    write(out_unitp,*) 'test RWU (working unit): ',RWU_Write(RWU1,WithUnit=.TRUE.,WorkingUnit=.TRUE.)
    write(out_unitp,*) 'test RWU (writing unit): ',RWU_Write(RWU1,WithUnit=.FALSE.,WorkingUnit=.FALSE.)
    write(out_unitp,*) 'test RWU (writing unit): ',RWU_Write(RWU1,WithUnit=.TRUE.,WorkingUnit=.FALSE.)
    RWU2 = convRWU_TO_RWU(RWU1)
    write(out_unitp,*) 'test RWU (after conv, working unit default) :',RWU2
    RWU2 = convRWU_TO_RWU(RWU1,WorkingUnit=.TRUE.)
    write(out_unitp,*) 'test RWU (after conv, working unit)         :',RWU2
    RWU2 = convRWU_TO_RWU(RWU1,WorkingUnit=.FALSE.)
    write(out_unitp,*) 'test RWU (after conv, writing unit)         :',RWU2

    Rval = convRWU_TO_R(RWU1)
    write(out_unitp,*) 'test Rval (after conv, working unit default):',Rval
    Rval = RWU1
    write(out_unitp,*) 'test Rval (after conv, = default):           ',Rval
    Rval = convRWU_TO_R(RWU1,WorkingUnit=.TRUE.)
    write(out_unitp,*) 'test Rval (after conv, working unit)        :',Rval
    Rval = convRWU_TO_R(RWU1,WorkingUnit=.FALSE.)
    write(out_unitp,*) 'test Rval (after conv, writing unit)        :',Rval

    write(out_unitp,*) "======================================"
    flush(out_unitp)


    write(out_unitp,*) "======================================"
    write(out_unitp,*) "== test with a quantity with different case =="

    RWU1 = REAL_WU(ONE,'au','e')

    write(out_unitp,*) 'test RWU (without conv): ',RWU1
    write(out_unitp,*) 'test RWU (working unit): ',RWU_Write(RWU1,WithUnit=.FALSE.,WorkingUnit=.TRUE.)
    write(out_unitp,*) 'test RWU (working unit): ',RWU_Write(RWU1,WithUnit=.TRUE.,WorkingUnit=.TRUE.)
    write(out_unitp,*) 'test RWU (writing unit): ',RWU_Write(RWU1,WithUnit=.FALSE.,WorkingUnit=.FALSE.)
    write(out_unitp,*) 'test RWU (writing unit): ',RWU_Write(RWU1,WithUnit=.TRUE.,WorkingUnit=.FALSE.)
    RWU2 = convRWU_TO_RWU(RWU1)
    write(out_unitp,*) 'test RWU (after conv, working unit default) :',RWU2
    RWU2 = convRWU_TO_RWU(RWU1,WorkingUnit=.TRUE.)
    write(out_unitp,*) 'test RWU (after conv, working unit)         :',RWU2
    RWU2 = convRWU_TO_RWU(RWU1,WorkingUnit=.FALSE.)
    write(out_unitp,*) 'test RWU (after conv, writing unit)         :',RWU2

    Rval = convRWU_TO_R(RWU1)
    write(out_unitp,*) 'test Rval (after conv, working unit default):',Rval
    Rval = RWU1
    write(out_unitp,*) 'test Rval (after conv, = default):           ',Rval
    Rval = convRWU_TO_R(RWU1,WorkingUnit=.TRUE.)
    write(out_unitp,*) 'test Rval (after conv, working unit)        :',Rval
    Rval = convRWU_TO_R(RWU1,WorkingUnit=.FALSE.)
    write(out_unitp,*) 'test Rval (after conv, writing unit)        :',Rval

    write(out_unitp,*) "======================================"
    flush(out_unitp)

    write(out_unitp,*) "======================================"
    write(out_unitp,*) "===== test with with different case =="

    RWU1 = REAL_WU(ONE,'ANGS','l')

    write(out_unitp,*) 'test RWU (without conv): ',RWU1
    write(out_unitp,*) 'test RWU (working unit): ',RWU_Write(RWU1,WithUnit=.FALSE.,WorkingUnit=.TRUE.)
    write(out_unitp,*) 'test RWU (working unit): ',RWU_Write(RWU1,WithUnit=.TRUE.,WorkingUnit=.TRUE.)
    write(out_unitp,*) 'test RWU (writing unit): ',RWU_Write(RWU1,WithUnit=.FALSE.,WorkingUnit=.FALSE.)
    write(out_unitp,*) 'test RWU (writing unit): ',RWU_Write(RWU1,WithUnit=.TRUE.,WorkingUnit=.FALSE.)
    RWU2 = convRWU_TO_RWU(RWU1)
    write(out_unitp,*) 'test RWU (after conv, working unit default) :',RWU2
    RWU2 = convRWU_TO_RWU(RWU1,WorkingUnit=.TRUE.)
    write(out_unitp,*) 'test RWU (after conv, working unit)         :',RWU2
    RWU2 = convRWU_TO_RWU(RWU1,WorkingUnit=.FALSE.)
    write(out_unitp,*) 'test RWU (after conv, writing unit)         :',RWU2

    Rval = convRWU_TO_R(RWU1)
    write(out_unitp,*) 'test Rval (after conv, working unit default):',Rval
    Rval = RWU1
    write(out_unitp,*) 'test Rval (after conv, = default):           ',Rval
    Rval = convRWU_TO_R(RWU1,WorkingUnit=.TRUE.)
    write(out_unitp,*) 'test Rval (after conv, working unit)        :',Rval
    Rval = convRWU_TO_R(RWU1,WorkingUnit=.FALSE.)
    write(out_unitp,*) 'test Rval (after conv, writing unit)        :',Rval

    write(out_unitp,*) "======================================"
    flush(out_unitp)

    CALL dealloc_TabConvRWU_dim1(Tab_conv_FOR_quantity)

  END SUBROUTINE Test_RWU


  END MODULE mod_RealWithUnit

