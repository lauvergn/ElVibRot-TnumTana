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
 MODULE mod_Frac
 USE mod_NumParameters
 IMPLICIT NONE

 ! fraction num/den
 TYPE Frac_t
   integer  :: num      =  0
   integer  :: den      =  1
 CONTAINS
   PROCEDURE, PRIVATE, PASS(frac)   :: PLUS_frac
   PROCEDURE, PRIVATE, PASS(frac1)  :: frac1_PLUS_frac2
   PROCEDURE, PRIVATE, PASS(frac1)  :: frac1_PLUS_Int2
   PROCEDURE, PRIVATE, PASS(frac2)  :: Int1_PLUS_frac2
   GENERIC,   PUBLIC  :: operator(+) => PLUS_frac,frac1_PLUS_frac2,     &
                                        frac1_PLUS_Int2,Int1_PLUS_frac2

   PROCEDURE, PRIVATE, PASS(frac)   :: MINUS_frac
   PROCEDURE, PRIVATE, PASS(frac1)  :: frac1_MINUS_frac2
   PROCEDURE, PRIVATE, PASS(frac1)  :: frac1_MINUS_Int2
   PROCEDURE, PRIVATE, PASS(frac2)  :: Int1_MINUS_frac2
   GENERIC,   PUBLIC  :: operator(-) => MINUS_frac,frac1_MINUS_frac2,   &
                                        frac1_MINUS_Int2,Int1_MINUS_frac2

   PROCEDURE, PRIVATE, PASS(frac1)  :: frac1_TIME_frac2
   PROCEDURE, PRIVATE, PASS(frac1)  :: frac1_TIME_Int2
   PROCEDURE, PRIVATE, PASS(frac2)  :: Int1_TIME_frac2
   GENERIC,   PUBLIC  :: operator(*) => frac1_TIME_frac2,frac1_TIME_Int2,Int1_TIME_frac2

   PROCEDURE, PRIVATE, PASS(frac1)  :: frac1_DIVIDEBY_frac2
   PROCEDURE, PRIVATE, PASS(frac1)  :: frac1_DIVIDEBY_Int2
   PROCEDURE, PRIVATE, PASS(frac2)  :: Int1_DIVIDEBY_frac2
   GENERIC,   PUBLIC  :: operator(/) => frac1_DIVIDEBY_frac2,frac1_DIVIDEBY_Int2,Int1_DIVIDEBY_frac2

   PROCEDURE, PRIVATE, PASS(frac)   :: Int_TO_frac
   PROCEDURE, PRIVATE, PASS(frac)   :: string_TO_frac
   GENERIC,   PUBLIC  :: assignment(=) => Int_TO_frac,string_TO_frac

   PROCEDURE, PRIVATE, PASS(frac1)  :: frac1_EQ_frac2
   PROCEDURE, PRIVATE, PASS(frac1)  :: frac1_EQ_Int2
   PROCEDURE, PRIVATE, PASS(frac2)  :: Int1_EQ_frac2
   GENERIC,   PUBLIC  :: operator(==) => frac1_EQ_frac2,frac1_EQ_Int2,Int1_EQ_frac2

   PROCEDURE, PRIVATE, PASS(frac1)  :: frac1_NEQ_frac2
   PROCEDURE, PRIVATE, PASS(frac1)  :: frac1_NEQ_Int2
   PROCEDURE, PRIVATE, PASS(frac2)  :: Int1_NEQ_frac2
   GENERIC,   PUBLIC  :: operator(/=) => frac1_NEQ_frac2,frac1_NEQ_Int2,Int1_NEQ_frac2

   PROCEDURE, PRIVATE, PASS(frac1)  :: frac1_GT_frac2
   PROCEDURE, PRIVATE, PASS(frac1)  :: frac1_GT_Int2
   PROCEDURE, PRIVATE, PASS(frac2)  :: Int1_GT_frac2
   GENERIC,   PUBLIC  :: operator(>) => frac1_GT_frac2,frac1_GT_Int2,Int1_GT_frac2

   PROCEDURE, PRIVATE, PASS(frac1)  :: frac1_LT_frac2
   PROCEDURE, PRIVATE, PASS(frac1)  :: frac1_LT_Int2
   PROCEDURE, PRIVATE, PASS(frac2)  :: Int1_LT_frac2
   GENERIC,   PUBLIC  :: operator(<) => frac1_LT_frac2,frac1_LT_Int2,Int1_LT_frac2

   PROCEDURE, PRIVATE, PASS(frac1)  :: frac1_GE_frac2
   PROCEDURE, PRIVATE, PASS(frac1)  :: frac1_GE_Int2
   PROCEDURE, PRIVATE, PASS(frac2)  :: Int1_GE_frac2
   GENERIC,   PUBLIC  :: operator(>=) => frac1_GE_frac2,frac1_GE_Int2,Int1_GE_frac2

   PROCEDURE, PRIVATE, PASS(frac1)  :: frac1_LE_frac2
   PROCEDURE, PRIVATE, PASS(frac1)  :: frac1_LE_Int2
   PROCEDURE, PRIVATE, PASS(frac2)  :: Int1_LE_frac2
   GENERIC,   PUBLIC  :: operator(<=) => frac1_LE_frac2,frac1_LE_Int2,Int1_LE_frac2
 END TYPE Frac_t

 INTERFACE Frac_t
   MODULE PROCEDURE construct_Frac
 END INTERFACE Frac_t

PRIVATE
PUBLIC :: Frac_t,frac_IS_integer
PUBLIC :: frac_TO_string,frac_TO_real
PUBLIC :: frac_simplification

 CONTAINS

  ELEMENTAL FUNCTION frac_IS_integer(frac) RESULT(test)
   TYPE(Frac_t), intent(in) :: frac
   logical                       :: test

   test = (frac%den == 1) .OR. (frac%num == 0)

 END FUNCTION frac_IS_integer

 ELEMENTAL SUBROUTINE Int_TO_frac(frac,i)
   CLASS(Frac_t), intent(inout) :: frac
   integer,           intent(in)    :: i

   frac%num = i
   frac%den = 1

 END SUBROUTINE Int_TO_frac


 ELEMENTAL FUNCTION construct_Frac(num,den) RESULT(frac)
   TYPE(Frac_t)                :: frac
   integer,           intent(in)    :: num,den

   frac%num = num
   frac%den = den

   CALL frac_simplification(frac%num,frac%den)

 END FUNCTION construct_Frac

 SUBROUTINE string_TO_frac(frac,String)
   character (len=*),  intent(in)    :: String
   CLASS(Frac_t), intent(inout) :: frac

   integer :: islash

   IF (len_trim(String) == 0) THEN
     frac = 0
   ELSE
     islash = index(String,'/')
     IF (islash == 0 .OR. islash == 1 .OR. islash == len(String)) THEN
       write(out_unitp,*) ' ERROR in string_TO_frac'
       write(out_unitp,*) ' The string does not contain a "/" or its postion is wrong'
       write(out_unitp,*) ' String: ',String
       STOP ' ERROR in string_TO_frac'
     END IF
     read(String(1:islash-1),*) frac%num
     read(String(islash+1:len(String)),*) frac%den
   END IF

   CALL frac_simplification(frac%num,frac%den)

 END SUBROUTINE string_TO_frac
 FUNCTION frac_TO_string(frac) RESULT(String)
   USE mod_string, only : String_TO_String,int_TO_char
   character (len=:), allocatable  :: String
   TYPE(Frac_t), intent(in)   :: frac

   IF (frac%den == 1) THEN
     String = int_TO_char(frac%num)
   ELSE
     String = String_TO_String( int_TO_char(frac%num) // '/' // int_TO_char(frac%den) )
   END IF

 END FUNCTION frac_TO_string
 ELEMENTAL FUNCTION frac_TO_real(frac) RESULT(R)
 USE mod_NumParameters, only : Rkind
   real(kind=Rkind)              :: R
   TYPE(Frac_t), intent(in) :: frac

   R = real(frac%num,kind=Rkind) / real(frac%den,kind=Rkind)

 END FUNCTION frac_TO_real

 ELEMENTAL FUNCTION frac1_EQ_frac2(frac1,frac2) RESULT(leq)
   logical                       :: leq
   CLASS(Frac_t), intent(in) :: frac1,frac2

   TYPE(Frac_t) :: frac

   frac = frac1-frac2
   leq  = (frac%num == 0)

 END FUNCTION frac1_EQ_frac2
 ELEMENTAL FUNCTION frac1_EQ_Int2(frac1,Int2) RESULT(leq)
   logical                       :: leq
   CLASS(Frac_t), intent(in) :: frac1
   Integer,           intent(in) :: Int2

   TYPE(Frac_t) :: frac

   frac = frac1-Frac_t(Int2,1)
   leq  = (frac%num == 0)

 END FUNCTION frac1_EQ_Int2
 ELEMENTAL FUNCTION Int1_EQ_frac2(Int1,frac2) RESULT(leq)
   logical                        :: leq
   Integer,            intent(in) :: Int1
   CLASS(Frac_t), intent(in) :: frac2

   TYPE(Frac_t) :: frac

   frac = Frac_t(Int1,1) - frac2
   leq  = (frac%num == 0)

 END FUNCTION Int1_EQ_frac2
 ELEMENTAL FUNCTION frac1_NEQ_frac2(frac1,frac2) RESULT(leq)
   logical                       :: leq
   CLASS(Frac_t), intent(in) :: frac1,frac2

   TYPE(Frac_t) :: frac

   frac = frac1-frac2
   leq  = (frac%num /= 0)

 END FUNCTION frac1_NEQ_frac2
 ELEMENTAL FUNCTION frac1_NEQ_Int2(frac1,Int2) RESULT(leq)
   logical                       :: leq
   CLASS(Frac_t), intent(in) :: frac1
   Integer,           intent(in) :: Int2

   TYPE(Frac_t) :: frac

   frac = frac1-Frac_t(Int2,1)
   leq  = (frac%num /= 0)

 END FUNCTION frac1_NEQ_Int2
 ELEMENTAL FUNCTION Int1_NEQ_frac2(Int1,frac2) RESULT(leq)
   logical                        :: leq
   Integer,            intent(in) :: Int1
   CLASS(Frac_t), intent(in) :: frac2

   TYPE(Frac_t) :: frac

   frac = Frac_t(Int1,1) - frac2
   leq  = (frac%num /= 0)

 END FUNCTION Int1_NEQ_frac2
 ELEMENTAL FUNCTION frac1_GT_frac2(frac1,frac2) RESULT(leq)
   logical                       :: leq
   CLASS(Frac_t), intent(in) :: frac1,frac2

   TYPE(Frac_t) :: frac

   frac = frac1-frac2
   leq  = (frac%num > 0)

 END FUNCTION frac1_GT_frac2
 ELEMENTAL FUNCTION frac1_GT_Int2(frac1,Int2) RESULT(leq)
   logical                       :: leq
   CLASS(Frac_t), intent(in) :: frac1
   Integer,           intent(in) :: Int2

   TYPE(Frac_t) :: frac

   frac = frac1-Frac_t(Int2,1)
   leq  = (frac%num > 0)

 END FUNCTION frac1_GT_Int2
 ELEMENTAL FUNCTION Int1_GT_frac2(Int1,frac2) RESULT(leq)
   logical                        :: leq
   Integer,            intent(in) :: Int1
   CLASS(Frac_t), intent(in) :: frac2

   TYPE(Frac_t) :: frac

   frac = Frac_t(Int1,1) - frac2
   leq  = (frac%num > 0)

 END FUNCTION Int1_GT_frac2
 ELEMENTAL FUNCTION frac1_LT_frac2(frac1,frac2) RESULT(leq)
   logical                       :: leq
   CLASS(Frac_t), intent(in) :: frac1,frac2

   TYPE(Frac_t) :: frac

   frac = frac1-frac2
   leq  = (frac%num < 0)

 END FUNCTION frac1_LT_frac2
 ELEMENTAL FUNCTION frac1_LT_Int2(frac1,Int2) RESULT(leq)
   logical                       :: leq
   CLASS(Frac_t), intent(in) :: frac1
   Integer,           intent(in) :: Int2

   TYPE(Frac_t) :: frac

   frac = frac1-Frac_t(Int2,1)
   leq  = (frac%num < 0)

 END FUNCTION frac1_LT_Int2
 ELEMENTAL FUNCTION Int1_LT_frac2(Int1,frac2) RESULT(leq)
   logical                        :: leq
   Integer,            intent(in) :: Int1
   CLASS(Frac_t), intent(in) :: frac2

   TYPE(Frac_t) :: frac

   frac = Frac_t(Int1,1) - frac2
   leq  = (frac%num < 0)

 END FUNCTION Int1_LT_frac2
 ELEMENTAL FUNCTION frac1_GE_frac2(frac1,frac2) RESULT(leq)
   logical                       :: leq
   CLASS(Frac_t), intent(in) :: frac1,frac2

   TYPE(Frac_t) :: frac

   frac = frac1-frac2
   leq  = (frac%num >= 0)

 END FUNCTION frac1_GE_frac2
 ELEMENTAL FUNCTION frac1_GE_Int2(frac1,Int2) RESULT(leq)
   logical                       :: leq
   CLASS(Frac_t), intent(in) :: frac1
   Integer,           intent(in) :: Int2

   TYPE(Frac_t) :: frac

   frac = frac1-Frac_t(Int2,1)
   leq  = (frac%num >= 0)

 END FUNCTION frac1_GE_Int2
 ELEMENTAL FUNCTION Int1_GE_frac2(Int1,frac2) RESULT(leq)
   logical                        :: leq
   Integer,            intent(in) :: Int1
   CLASS(Frac_t), intent(in) :: frac2

   TYPE(Frac_t) :: frac

   frac = Frac_t(Int1,1) - frac2
   leq  = (frac%num >= 0)

 END FUNCTION Int1_GE_frac2
 ELEMENTAL FUNCTION frac1_LE_frac2(frac1,frac2) RESULT(leq)
   logical                       :: leq
   CLASS(Frac_t), intent(in) :: frac1,frac2

   TYPE(Frac_t) :: frac

   frac = frac1-frac2
   leq  = (frac%num <= 0)

 END FUNCTION frac1_LE_frac2
 ELEMENTAL FUNCTION frac1_LE_Int2(frac1,Int2) RESULT(leq)
   logical                       :: leq
   CLASS(Frac_t), intent(in) :: frac1
   Integer,           intent(in) :: Int2

   TYPE(Frac_t) :: frac

   frac = frac1-Frac_t(Int2,1)
   leq  = (frac%num <= 0)

 END FUNCTION frac1_LE_Int2
 ELEMENTAL FUNCTION Int1_LE_frac2(Int1,frac2) RESULT(leq)
   logical                        :: leq
   Integer,            intent(in) :: Int1
   CLASS(Frac_t), intent(in) :: frac2

   TYPE(Frac_t) :: frac

   frac = Frac_t(Int1,1) - frac2
   leq  = (frac%num <= 0)

 END FUNCTION Int1_LE_frac2
 ELEMENTAL FUNCTION PLUS_frac(frac) RESULT(Resfrac)
   TYPE(Frac_t)              :: Resfrac
   CLASS(Frac_t), intent(in) :: frac

   ResFrac = frac

 END FUNCTION PLUS_frac
 ELEMENTAL FUNCTION frac1_PLUS_frac2(frac1,frac2) RESULT(frac)
   TYPE(Frac_t)             :: frac
   CLASS(Frac_t), intent(in) :: frac1,frac2

   Frac%num = frac1%num * frac2%den + frac1%den * frac2%num
   Frac%den = frac1%den * frac2%den

   CALL frac_simplification(Frac%num, Frac%den)

 END FUNCTION frac1_PLUS_frac2
 ELEMENTAL FUNCTION frac1_PLUS_Int2(frac1,Int2) RESULT(frac)
   TYPE(Frac_t)             :: frac
   CLASS(Frac_t), intent(in) :: frac1
   Integer,           intent(in) :: Int2

   Frac%num = frac1%num + frac1%den * Int2
   Frac%den = frac1%den

   CALL frac_simplification(Frac%num, Frac%den)

 END FUNCTION frac1_PLUS_Int2
 ELEMENTAL FUNCTION Int1_PLUS_frac2(Int1,frac2) RESULT(frac)
   TYPE(Frac_t)             :: frac
   CLASS(Frac_t), intent(in) :: frac2
   Integer,           intent(in) :: Int1

   Frac%num = frac2%num + frac2%den * Int1
   Frac%den = frac2%den

   CALL frac_simplification(Frac%num, Frac%den)

 END FUNCTION Int1_PLUS_frac2
 ELEMENTAL FUNCTION MINUS_frac(frac) RESULT(Resfrac)
   TYPE(Frac_t)              :: Resfrac
   CLASS(Frac_t), intent(in) :: frac

   ResFrac%num = -frac%num
   ResFrac%den = frac%den

 END FUNCTION MINUS_frac
 ELEMENTAL FUNCTION frac1_MINUS_frac2(frac1,frac2) RESULT(frac)
   TYPE(Frac_t)             :: frac
   CLASS(Frac_t), intent(in) :: frac1,frac2

   Frac%num = frac1%num * frac2%den - frac1%den * frac2%num
   Frac%den = frac1%den * frac2%den

   CALL frac_simplification(Frac%num, Frac%den)

 END FUNCTION frac1_MINUS_frac2
 ELEMENTAL FUNCTION frac1_MINUS_Int2(frac1,Int2) RESULT(frac)
   TYPE(Frac_t)             :: frac
   CLASS(Frac_t), intent(in) :: frac1
   Integer,           intent(in) :: Int2

   Frac%num = frac1%num - frac1%den * Int2
   Frac%den = frac1%den

   CALL frac_simplification(Frac%num, Frac%den)

 END FUNCTION frac1_MINUS_Int2
 ELEMENTAL FUNCTION Int1_MINUS_frac2(Int1,frac2) RESULT(frac)
   TYPE(Frac_t)             :: frac
   CLASS(Frac_t), intent(in) :: frac2
   Integer,           intent(in) :: Int1

   Frac%num = frac2%den * Int1 - frac2%num
   Frac%den = frac2%den

   CALL frac_simplification(Frac%num, Frac%den)

 END FUNCTION Int1_MINUS_frac2

 ELEMENTAL FUNCTION frac1_TIME_frac2(frac1,frac2) RESULT(frac)
   TYPE(Frac_t)             :: frac
   CLASS(Frac_t), intent(in) :: frac1,frac2

   Frac%num = frac1%num * frac2%num
   Frac%den = frac1%den * frac2%den

   CALL frac_simplification(Frac%num, Frac%den)

 END FUNCTION frac1_TIME_frac2
 ELEMENTAL FUNCTION Int1_TIME_frac2(Int1,frac2) RESULT(frac)
   TYPE(Frac_t)             :: frac
   CLASS(Frac_t), intent(in) :: frac2
   Integer,           intent(in) :: Int1

   Frac%num = Int1 * frac2%num
   Frac%den = frac2%den

   CALL frac_simplification(Frac%num, Frac%den)

 END FUNCTION Int1_TIME_frac2
 ELEMENTAL FUNCTION frac1_TIME_Int2(frac1,Int2) RESULT(frac)
   TYPE(Frac_t)             :: frac
   CLASS(Frac_t), intent(in) :: frac1
   Integer,           intent(in) :: Int2

   Frac%num = frac1%num * Int2
   Frac%den = frac1%den

   CALL frac_simplification(Frac%num, Frac%den)

 END FUNCTION frac1_TIME_Int2
 ELEMENTAL FUNCTION frac1_DIVIDEBY_frac2(frac1,frac2) RESULT(frac)
   TYPE(Frac_t)             :: frac
   CLASS(Frac_t), intent(in) :: frac1,frac2

   Frac%num = frac1%num * frac2%den
   Frac%den = frac1%den * frac2%num

   CALL frac_simplification(Frac%num, Frac%den)

 END FUNCTION frac1_DIVIDEBY_frac2
 ELEMENTAL FUNCTION Int1_DIVIDEBY_frac2(Int1,frac2) RESULT(frac)
   TYPE(Frac_t)             :: frac
   CLASS(Frac_t), intent(in) :: frac2
   Integer,           intent(in) :: Int1

   Frac%num = Int1 * frac2%den
   Frac%den = frac2%num

   CALL frac_simplification(Frac%num, Frac%den)

 END FUNCTION Int1_DIVIDEBY_frac2
 ELEMENTAL FUNCTION frac1_DIVIDEBY_Int2(frac1,Int2) RESULT(frac)
   TYPE(Frac_t)             :: frac
   CLASS(Frac_t), intent(in) :: frac1
   Integer,           intent(in) :: Int2

   Frac%num = frac1%num
   Frac%den = frac1%den * Int2

   CALL frac_simplification(Frac%num, Frac%den)

 END FUNCTION frac1_DIVIDEBY_Int2
!================================================================
!    greatest common divisor
!================================================================
    ELEMENTAL FUNCTION gcd(a, b)
    integer              :: gcd
    integer, intent (in) :: a,b

    integer :: aa,bb,t

    aa = a
    bb = b
    DO
      IF (bb == 0) EXIT
      t = bb
      bb = mod(aa,bb)
      aa = t
    END DO
    gcd = aa

    END FUNCTION gcd
!================================================================
!    fraction simplification a/b
!================================================================
    ELEMENTAL SUBROUTINE frac_simplification(a, b)
    integer, intent (inout) :: a,b

    integer :: aa,bb,t

    aa = a
    bb = b
    t = gcd(aa,bb)

    a = aa/t
    b = bb/t

    IF (b < 0) THEN
      b = -b
      a = -a
    END IF

    END SUBROUTINE frac_simplification

 END MODULE mod_Frac

