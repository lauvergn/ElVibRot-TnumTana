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
MODULE mod_NumParameters
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64,real128,int32,int64
  IMPLICIT NONE

  !integer, parameter :: Rkind        = real128
  integer, parameter :: Rkind        = real64

  integer, parameter :: Ikind        = int32
  integer, parameter :: ILkind       = int64
  integer, parameter :: Name_len     = 20
  integer, parameter :: Name_longlen = 50
  integer, parameter :: Line_len     = 255
  integer, parameter :: error_l      = 80


  real (kind=Rkind), parameter :: ZERO    = 0._Rkind
  real (kind=Rkind), parameter :: ONE     = 1._Rkind
  real (kind=Rkind), parameter :: TWO     = 2._Rkind
  real (kind=Rkind), parameter :: THREE   = 3._Rkind
  real (kind=Rkind), parameter :: FOUR    = 4._Rkind
  real (kind=Rkind), parameter :: FIVE    = 5._Rkind
  real (kind=Rkind), parameter :: SIX     = 6._Rkind
  real (kind=Rkind), parameter :: SEVEN   = 7._Rkind
  real (kind=Rkind), parameter :: EIGHT   = 8._Rkind
  real (kind=Rkind), parameter :: NINE    = 9._Rkind
  real (kind=Rkind), parameter :: TEN     = 10._Rkind
  real (kind=Rkind), parameter :: ELEVEN  = 11._Rkind
  real (kind=Rkind), parameter :: TWELVE  = 12._Rkind
  real (kind=Rkind), parameter :: HUNDRED = 100._Rkind

  real (kind=Rkind), parameter :: HALF      = ONE/TWO
  real (kind=Rkind), parameter :: THIRD     = ONE/THREE
  real (kind=Rkind), parameter :: FOURTH    = ONE/FOUR
  real (kind=Rkind), parameter :: QUARTER   = ONE/FOUR
  real (kind=Rkind), parameter :: FIFTH     = ONE/FIVE
  real (kind=Rkind), parameter :: SIXTH     = ONE/SIX
  real (kind=Rkind), parameter :: ONETENTH  = ONE/TEN
  real (kind=Rkind), parameter :: TWOTENTHS = TWO/TEN

  real (kind=Rkind), parameter ::                                   &
   pi = 3.14159265358979323846264338327950288419716939937511_Rkind

  complex (kind=Rkind), parameter :: EYE      = (0._Rkind,1._Rkind)
  complex (kind=Rkind), parameter :: CZERO    = (0._Rkind,0._Rkind)
  complex (kind=Rkind), parameter :: CONE     = (1._Rkind,0._Rkind)
  complex (kind=Rkind), parameter :: CTWO     = (2._Rkind,0._Rkind)
  complex (kind=Rkind), parameter :: CHALF    = (0.5_Rkind,0._Rkind)

  integer :: print_level = 0        ! 0 minimal, 1 default, 2 large, -1 nothing

  character (len=Name_longlen) :: EneIO_format = "f20.5"
  character (len=Name_longlen) :: RMatIO_format = "f18.10"
  character (len=Name_longlen) :: CMatIO_format = "'(',f15.7,' +i',f15.7,')'"

  !integer :: in_unitp  = 5 ! Unit for input and the ouptput files
  !integer :: out_unitp = 6 ! Unit for input and the ouptput files
  integer :: in_unitp  = INPUT_UNIT  ! Unit for the ouptput files, with the ISO_FORTRAN_ENV
  integer :: out_unitp = OUTPUT_UNIT ! Unit for the input, with the ISO_FORTRAN_ENV

END MODULE mod_NumParameters
!-------------------------------------------------------------------------------
!   This is a Fortran translation of the 64-bit version of
!   the Mersenne Twister pseudorandom number generator
!
!   Before using, initialize the state by using
!       call init_genrand64(seed)
!   or
!       call init_by_array64(init_key)
!
!   Translated from C-program for MT19937-64 (2004/9/29 version)
!   originally coded by Takuji Nishimura and Makoto Matsumoto
!   http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt64.html
!
!   Fortran translation by Rémi Piatek
!   The University of Copenhagen
!   Department of Economics
!   email: {first}.{last}@econ.ku.dk
!
!-------------------------------------------------------------------------------
!   A C-program for MT19937-64 (2004/9/29 version).
!   Coded by Takuji Nishimura and Makoto Matsumoto.
!
!   This is a 64-bit version of Mersenne Twister pseudorandom number
!   generator.
!
!   Before using, initialize the state by using init_genrand64(seed)
!   or init_by_array64(init_key, key_length).
!
!   Copyright (C) 2004, Makoto Matsumoto and Takuji Nishimura,
!   All rights reserved.
!
!   Redistribution and use in source and binary forms, with or without
!   modification, are permitted provided that the following conditions
!   are met:
!
!     1. Redistributions of source code must retain the above copyright
!        notice, this list of conditions and the following disclaimer.
!
!     2. Redistributions in binary form must reproduce the above copyright
!        notice, this list of conditions and the following disclaimer in the
!        documentation and/or other materials provided with the distribution.
!
!     3. The names of its contributors may not be used to endorse or promote
!        products derived from this software without specific prior written
!        permission.
!
!   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
!   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
!   OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
!   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
!   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
!   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
!   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
!   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
!   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!   References:
!   T. Nishimura, ``Tables of 64-bit Mersenne Twisters''
!     ACM Transactions on Modeling and
!     Computer Simulation 10. (2000) 348--357.
!   M. Matsumoto and T. Nishimura,
!     ``Mersenne Twister: a 623-dimensionally equidistributed
!       uniform pseudorandom number generator''
!     ACM Transactions on Modeling and
!     Computer Simulation 8. (Jan. 1998) 3--30.
!
!   Any feedback is very welcome.
!   http://www.math.hiroshima-u.ac.jp/~m-mat/MT/emt.html
!   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove spaces)
!-------------------------------------------------------------------------------

module mt19937_64
  use, intrinsic :: iso_fortran_env
  implicit none

  private
  public :: init_genrand64
  public :: init_by_array64
  public :: genrand64_real1
  public :: genrand64_real2
  public :: genrand64_real3

! NOTE: genrand64_int64 is kept private, as it generates different numbers
!       compared to the original C code. This is because the original C code
!       uses unsigned integers, while Fortran relies on signed integers.
!       This, however, has no impact on the generation of real numbers
!       (they are identical to those produced by the original C code).
! public :: genrand64_int64

  integer, parameter :: r8 = real64
  integer, parameter :: i8 = int64

  integer(i8), parameter :: nn       = 312_i8
  integer(i8), parameter :: mm       = 156_i8
  integer(i8), parameter :: seed_def = 5489_i8
  integer(i8), parameter :: matrix_a = -5403634167711393303_i8
  integer(i8), parameter :: um       = -2147483648_i8 ! most significant 33 bits
  integer(i8), parameter :: lm       =  2147483647_i8 ! least significant 31 bits

  real(r8),    parameter :: pi253_1  = 1._r8/(2._r8**53 - 1._r8)
  real(r8),    parameter :: pi253    = 1._r8/(2._r8**53)
  real(r8),    parameter :: pi252    = 1._r8/(2._r8**52)

  integer(i8) :: mt(nn)       ! array for the state vector
  integer     :: mti = nn+1   ! mti==nn+1 means mt(nn) is not initialized


contains


  !-----------------------------------------------------------------------------
  ! Initializes mt(nn) with a seed
  subroutine init_genrand64(seed)
    implicit none
    integer(i8), intent(in) :: seed
    integer :: i

    mt(1) = seed
    do i = 1, nn-1
      mt(i+1) = 6364136223846793005_i8 * ieor(mt(i), ishft(mt(i), -62)) + i
    end do

    mti = nn

  end subroutine init_genrand64


  !-----------------------------------------------------------------------------
  ! Initializes by an array with array-length
  !   init_key is the array for initializing keys

  subroutine init_by_array64(init_key)
    implicit none
    integer(i8), intent(in) :: init_key(:)
    integer(i8), parameter  :: c1 = 3935559000370003845_i8
    integer(i8), parameter  :: c2 = 2862933555777941757_i8
    integer(i8) :: i, j, k, kk, key_length

    call init_genrand64(19650218_i8)
    key_length = size(init_key)
    i = 1_i8; j = 0_i8
    k = max(nn, key_length)

    do kk = 1, k
      mt(i+1) = ieor(mt(i+1), c1 * ieor(mt(i), ishft(mt(i), -62))) &
                  + init_key(j+1) + j
      i = i+1; j = j+1
      if(i >= nn) then
        mt(1) = mt(nn)
        i = 1
      end if
      if(j >= key_length) j = 0
    end do

    do kk = 1, nn-1
      mt(i+1) = ieor(mt(i+1), c2 * ieor(mt(i), ishft(mt(i), -62))) - i
      i = i+1
      if(i >= nn) then
        mt(1) = mt(nn)
        i = 1
      end if
    end do

    mt(1) = ishft(1_i8, 63)  ! MSB is 1; assuring non-zero initial array

  end subroutine init_by_array64


  !-----------------------------------------------------------------------------
  ! Generates a random number on [-2^63, 2^63-1]-interval

  integer(r8) function genrand64_int64()
    implicit none
    integer(i8) :: mag01(0:1) = (/0_i8, matrix_a/)
    integer(i8) :: x
    integer     :: i

    if(mti >= nn) then ! generate nn words at one time

      ! if init_genrand64() has not been called, a default initial seed is used
      if(mti == nn+1) call init_genrand64(seed_def)

      do i = 1, nn-mm
        x = ior(iand(mt(i),um), iand(mt(i+1), lm))
        mt(i) = ieor(ieor(mt(i+mm), ishft(x, -1)), mag01(iand(x, 1_i8)))
      end do

      do i = nn-mm+1, nn-1
        x = ior(iand(mt(i), um), iand(mt(i+1), lm))
        mt(i) = ieor(ieor(mt(i+mm-nn), ishft(x, -1)), mag01(iand(x, 1_i8)))
      end do

      x = ior(iand(mt(nn), um), iand(mt(1), lm))
      mt(nn) = ieor(ieor(mt(mm), ishft(x, -1)), mag01(iand(x, 1_i8)))

      mti = 0

    end if

    mti = mti + 1
    x = mt(mti)

    x = ieor(x, iand(ishft(x,-29), 6148914691236517205_i8))
    x = ieor(x, iand(ishft(x, 17), 8202884508482404352_i8))
    x = ieor(x, iand(ishft(x, 37),   -2270628950310912_i8))
    x = ieor(x, ishft(x, -43))

    genrand64_int64 = x

  end function genrand64_int64


  !-----------------------------------------------------------------------------
  ! Generates a random number on [0,1]-real-interval

  real(r8) function genrand64_real1()
    implicit none

    genrand64_real1 = real(ishft(genrand64_int64(), -11), kind=r8) * pi253_1

  end function genrand64_real1


  !-----------------------------------------------------------------------------
  ! Generates a random number on [0,1)-real-interval

  real(r8) function genrand64_real2()
    implicit none

    genrand64_real2 = real(ishft(genrand64_int64(), -11), kind=r8) * pi253

  end function genrand64_real2


  !-----------------------------------------------------------------------------
  ! Generates a random number on (0,1)-real-interval

  real(r8) function genrand64_real3()
    implicit none

    genrand64_real3 = real(ishft(genrand64_int64(), -12), kind=r8)
    genrand64_real3 = (genrand64_real3 + 0.5_r8) * pi252

  end function genrand64_real3


end module mt19937_64
