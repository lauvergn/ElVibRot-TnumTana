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
 MODULE mod_FracInteger
 IMPLICIT NONE

 ! fraction num/den
 TYPE FracInteger
   integer  :: num      =  0
   integer  :: den      =  1
 END TYPE FracInteger

 INTERFACE operator (+)
  MODULE PROCEDURE frac1_PLUS_frac2,frac1_PLUS_Int2
 END INTERFACE
 INTERFACE operator (-)
  MODULE PROCEDURE frac1_MINUS_frac2,frac1_MINUS_Int2
 END INTERFACE
 INTERFACE operator (*)
  MODULE PROCEDURE frac1_TIME_frac2
 END INTERFACE
 INTERFACE operator (/)
  MODULE PROCEDURE frac1_DIVIDEBY_frac2
 END INTERFACE

 INTERFACE operator (==)
  MODULE PROCEDURE frac1_EQ_frac2,frac1_EQ_Int2
 END INTERFACE
 INTERFACE operator (/=)
  MODULE PROCEDURE frac1_NEQ_frac2,frac1_NEQ_Int2
 END INTERFACE

 INTERFACE operator (>)
  MODULE PROCEDURE frac1_GT_frac2,frac1_GT_Int2
 END INTERFACE
 INTERFACE operator (<)
  MODULE PROCEDURE frac1_LT_frac2,frac1_LT_Int2
 END INTERFACE
 INTERFACE operator (>=)
  MODULE PROCEDURE frac1_GE_frac2,frac1_GE_Int2
 END INTERFACE
 INTERFACE operator (<=)
  MODULE PROCEDURE frac1_LE_frac2,frac1_LE_Int2
 END INTERFACE

 INTERFACE assignment (=)
  MODULE PROCEDURE Int_TO_frac
 END INTERFACE

PRIVATE
PUBLIC :: FracInteger, test_FracInteger
PUBLIC :: operator (+),operator (-),operator (*),operator (/)
PUBLIC :: operator (==),operator (/=),operator (>),operator (<),operator (>=),operator (<=)
PUBLIC :: assignment (=),frac_TO_string,frac_TO_real
PUBLIC :: frac_simplification

 CONTAINS

  SUBROUTINE test_FracInteger()
   TYPE(FracInteger) :: frac1,frac2


   frac1 = FracInteger(1,2)
   frac2 = 2

   write(6,*) 'frac1: ',frac_TO_string(frac1)
   write(6,*) 'frac2: ',frac_TO_string(frac2)
   write(6,*) 'frac1: ',frac_TO_real(frac1)
   write(6,*)
   write(6,*) 'frac1+frac2: ',frac_TO_string(frac1+frac2)
   write(6,*) 'frac1-frac2: ',frac_TO_string(frac1-frac2)
   write(6,*) 'frac1*frac2: ',frac_TO_string(frac1*frac2)
   write(6,*) 'frac1/frac2: ',frac_TO_string(frac1/frac2)
   write(6,*)
   write(6,*) 'frac1 == frac2 ?: ',(frac1 == frac2)
   write(6,*) 'frac1 /= frac2 ?: ',(frac1 /= frac2)
   write(6,*) 'frac1 == frac1 ?: ',(frac1 == frac1)
   write(6,*) 'frac1 /= frac1 ?: ',(frac1 /= frac1)
   write(6,*)
   write(6,*) 'frac1 < frac2 ?: ',(frac1 < frac2)
   write(6,*) 'frac1 > frac2 ?: ',(frac1 > frac2)
   write(6,*) 'frac1 >= frac1 ?: ',(frac1 >= frac1)
   write(6,*) 'frac1 <= frac1 ?: ',(frac1 <= frac1)

   write(6,*)
   write(6,*)
   frac1 = FracInteger(1,6)
   frac2 = FracInteger(18,12)
   write(6,*) 'frac1: ',frac_TO_string(frac1)
   write(6,*) 'frac2: ',frac_TO_string(frac2)
   write(6,*) 'frac1+frac2: ',frac_TO_string(frac1+frac2)

   frac1 = 6
   write(6,*) 'frac1: ',frac_TO_string(frac1)


 END SUBROUTINE test_FracInteger

  SUBROUTINE Int_TO_frac(frac,i)
   TYPE(FracInteger), intent(inout) :: frac
   integer,           intent(in)    :: i

   frac = FracInteger(i,1)

 END SUBROUTINE Int_TO_frac

  FUNCTION frac_TO_string(frac) RESULT(String)
   USE mod_string, only : String_TO_String,int_TO_char
   character (len=:), allocatable  :: String
   TYPE(FracInteger), intent(in) :: frac

   IF (frac%den == 1) THEN
     String = int_TO_char(frac%num)
   ELSE
     String = String_TO_String( int_TO_char(frac%num) // '/' // int_TO_char(frac%den) )
   END IF

 END FUNCTION frac_TO_string
 FUNCTION frac_TO_real(frac) RESULT(R)
 USE mod_NumParameters, only : Rkind
   real(kind=Rkind)              :: R
   TYPE(FracInteger), intent(in) :: frac

   R = real(frac%num,kind=Rkind) / real(frac%den,kind=Rkind)

 END FUNCTION frac_TO_real

 FUNCTION frac1_EQ_frac2(frac1,frac2) RESULT(leq)
   logical                       :: leq
   TYPE(FracInteger), intent(in) :: frac1,frac2

   TYPE(FracInteger) :: frac

   frac = frac1-frac2
   leq  = (frac%num == 0)

 END FUNCTION frac1_EQ_frac2
 FUNCTION frac1_EQ_Int2(frac1,Int2) RESULT(leq)
   logical                       :: leq
   TYPE(FracInteger), intent(in) :: frac1
   Integer,           intent(in) :: Int2

   TYPE(FracInteger) :: frac

   frac = frac1-FracInteger(Int2,1)
   leq  = (frac%num == 0)

 END FUNCTION frac1_EQ_Int2
 FUNCTION frac1_NEQ_frac2(frac1,frac2) RESULT(leq)
   logical                       :: leq
   TYPE(FracInteger), intent(in) :: frac1,frac2

   TYPE(FracInteger) :: frac

   frac = frac1-frac2
   leq  = (frac%num /= 0)

 END FUNCTION frac1_NEQ_frac2
 FUNCTION frac1_NEQ_Int2(frac1,Int2) RESULT(leq)
   logical                       :: leq
   TYPE(FracInteger), intent(in) :: frac1
   Integer,           intent(in) :: Int2

   TYPE(FracInteger) :: frac

   frac = frac1-FracInteger(Int2,1)
   leq  = (frac%num /= 0)

 END FUNCTION frac1_NEQ_Int2
 FUNCTION frac1_GT_frac2(frac1,frac2) RESULT(leq)
   logical                       :: leq
   TYPE(FracInteger), intent(in) :: frac1,frac2

   TYPE(FracInteger) :: frac

   frac = frac1-frac2
   leq  = (frac%num > 0)

 END FUNCTION frac1_GT_frac2
 FUNCTION frac1_GT_Int2(frac1,Int2) RESULT(leq)
   logical                       :: leq
   TYPE(FracInteger), intent(in) :: frac1
   Integer,           intent(in) :: Int2

   TYPE(FracInteger) :: frac

   frac = frac1-FracInteger(Int2,1)
   leq  = (frac%num > 0)

 END FUNCTION frac1_GT_Int2
 FUNCTION frac1_LT_frac2(frac1,frac2) RESULT(leq)
   logical                       :: leq
   TYPE(FracInteger), intent(in) :: frac1,frac2

   TYPE(FracInteger) :: frac

   frac = frac1-frac2
   leq  = (frac%num < 0)

 END FUNCTION frac1_LT_frac2
 FUNCTION frac1_LT_Int2(frac1,Int2) RESULT(leq)
   logical                       :: leq
   TYPE(FracInteger), intent(in) :: frac1
   Integer,           intent(in) :: Int2

   TYPE(FracInteger) :: frac

   frac = frac1-FracInteger(Int2,1)
   leq  = (frac%num < 0)

 END FUNCTION frac1_LT_Int2
 FUNCTION frac1_GE_frac2(frac1,frac2) RESULT(leq)
   logical                       :: leq
   TYPE(FracInteger), intent(in) :: frac1,frac2

   TYPE(FracInteger) :: frac

   frac = frac1-frac2
   leq  = (frac%num >= 0)

 END FUNCTION frac1_GE_frac2
 FUNCTION frac1_GE_Int2(frac1,Int2) RESULT(leq)
   logical                       :: leq
   TYPE(FracInteger), intent(in) :: frac1
   Integer,           intent(in) :: Int2

   TYPE(FracInteger) :: frac

   frac = frac1-FracInteger(Int2,1)
   leq  = (frac%num >= 0)

 END FUNCTION frac1_GE_Int2
 FUNCTION frac1_LE_frac2(frac1,frac2) RESULT(leq)
   logical                       :: leq
   TYPE(FracInteger), intent(in) :: frac1,frac2

   TYPE(FracInteger) :: frac

   frac = frac1-frac2
   leq  = (frac%num <= 0)

 END FUNCTION frac1_LE_frac2
 FUNCTION frac1_LE_Int2(frac1,Int2) RESULT(leq)
   logical                       :: leq
   TYPE(FracInteger), intent(in) :: frac1
   Integer,           intent(in) :: Int2

   TYPE(FracInteger) :: frac

   frac = frac1-FracInteger(Int2,1)
   leq  = (frac%num <= 0)

 END FUNCTION frac1_LE_Int2

 FUNCTION frac1_PLUS_frac2(frac1,frac2) RESULT(frac)
   TYPE(FracInteger)             :: frac
   TYPE(FracInteger), intent(in) :: frac1,frac2

   Frac%num = frac1%num * frac2%den + frac1%den * frac2%num
   Frac%den = frac1%den * frac2%den

   CALL frac_simplification(Frac%num, Frac%den)

 END FUNCTION frac1_PLUS_frac2
 FUNCTION frac1_PLUS_Int2(frac1,Int2) RESULT(frac)
   TYPE(FracInteger)             :: frac
   TYPE(FracInteger), intent(in) :: frac1
   Integer,           intent(in) :: Int2

   Frac%num = frac1%num + frac1%den * Int2
   Frac%den = frac1%den

   CALL frac_simplification(Frac%num, Frac%den)

 END FUNCTION frac1_PLUS_Int2
 FUNCTION frac1_MINUS_frac2(frac1,frac2) RESULT(frac)
   TYPE(FracInteger)             :: frac
   TYPE(FracInteger), intent(in) :: frac1,frac2

   Frac%num = frac1%num * frac2%den - frac1%den * frac2%num
   Frac%den = frac1%den * frac2%den

   CALL frac_simplification(Frac%num, Frac%den)

 END FUNCTION frac1_MINUS_frac2
 FUNCTION frac1_MINUS_Int2(frac1,Int2) RESULT(frac)
   TYPE(FracInteger)             :: frac
   TYPE(FracInteger), intent(in) :: frac1
   Integer,           intent(in) :: Int2

   Frac%num = frac1%num - frac1%den * Int2
   Frac%den = frac1%den

   CALL frac_simplification(Frac%num, Frac%den)

 END FUNCTION frac1_MINUS_Int2
 FUNCTION frac1_TIME_frac2(frac1,frac2) RESULT(frac)
   TYPE(FracInteger)             :: frac
   TYPE(FracInteger), intent(in) :: frac1,frac2

   Frac%num = frac1%num * frac2%num
   Frac%den = frac1%den * frac2%den

   CALL frac_simplification(Frac%num, Frac%den)

 END FUNCTION frac1_TIME_frac2
 FUNCTION frac1_DIVIDEBY_frac2(frac1,frac2) RESULT(frac)
   TYPE(FracInteger)             :: frac
   TYPE(FracInteger), intent(in) :: frac1,frac2

   Frac%num = frac1%num * frac2%den
   Frac%den = frac1%den * frac2%num

   CALL frac_simplification(Frac%num, Frac%den)

 END FUNCTION frac1_DIVIDEBY_frac2

!================================================================
!    greatest common divisor
!================================================================
    FUNCTION gcd(a, b)
    integer :: gcd
    integer :: a,b

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
    SUBROUTINE frac_simplification(a, b)
    !integer :: gcd
    integer :: a,b

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

 END MODULE mod_FracInteger

