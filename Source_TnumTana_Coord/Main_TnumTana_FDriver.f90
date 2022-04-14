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
!      Tnum is written David Lauvergnat [1]
!      Tana is written by Mamadou Ndong [1] and David Lauvergnat [1]
!         with contributions
!          Emil Lund klinting (coupling with MidasCpp) [3]'
!
![1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France
![3]: Department of Chemistry, Aarhus University, DK-8000 Aarhus C, Denmark
!
!===========================================================================
!===========================================================================
PROGRAM Main_TnumTana_FDriver
  USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT,real64,real128,int32,int64
  IMPLICIT NONE

  integer, parameter :: Rk  = real64 ! 8

  integer                         :: nb_act,nb_cart,init_sub
  real (kind=Rk), allocatable     :: Qact(:),Qcart(:)
  real (kind=Rk)                  :: mass
  real (kind=Rk)                  :: EckartRot(3,3)
  integer                         :: Z,A
  character (len=10)              :: Atomic_Symbol
  integer                         :: InputUnit,ResUnit
  character (len=:), allocatable  :: InputName,OutputName
  logical                         :: X2Q
  integer                         :: nb_eval
  character (len=10)              :: time0,time1
  INTEGER                         :: Time_val0(8),Time_val1(8)


  integer :: i
  character (len=*), parameter :: name_sub='Main_TnumTana'

  CALL ReadArguments(InputName,OutputName,X2Q,nb_eval)

  open(newunit=ResUnit,file=OutputName)
  open(newunit=InputUnit,file=InputName)

  CALL Init_InputUnit_Driver(InputUnit)
  CALL Init_OutputUnit_Driver(ResUnit)
  write(OUTPUT_UNIT,*) 'InputUnit',InputUnit
  write(OUTPUT_UNIT,*) 'ResUnit',ResUnit

  CALL Init_TnumTana_FOR_Driver(nb_act,nb_cart,init_sub)
  write(OUTPUT_UNIT,*) 'nb_act,nb_cart,init_sub',nb_act,nb_cart,init_sub

  !=================================
  ! to get isotopic masses
  Z = -1
  A = -1
  Atomic_Symbol = 'C'
  CALL Tnum_get_mass(mass,Z,A,Atomic_Symbol)
  write(OUTPUT_UNIT,*) 'Z,A,Atomic_Symbol,mass',Z,A,' ',trim(Atomic_Symbol),' ',mass
  flush(OUTPUT_UNIT)

  Z = -1
  A = -1
  Atomic_Symbol = '6_13'
  CALL Tnum_get_mass(mass,Z,A,Atomic_Symbol)
  write(OUTPUT_UNIT,*) 'Z,A,Atomic_Symbol,mass',Z,A,' ',trim(Atomic_Symbol),' ',mass
  flush(OUTPUT_UNIT)

  Z = 8
  A = 17
  Atomic_Symbol = ''
  CALL Tnum_get_mass(mass,Z,A,Atomic_Symbol)
  write(OUTPUT_UNIT,*) 'Z,A,Atomic_Symbol,mass',Z,A,' ',trim(Atomic_Symbol),' ',mass
  flush(OUTPUT_UNIT)

  Z = 8
  A = -1
  Atomic_Symbol = ''
  CALL Tnum_get_mass(mass,Z,A,Atomic_Symbol)
  write(OUTPUT_UNIT,*) 'Z,A,Atomic_Symbol,mass',Z,A,' ',trim(Atomic_Symbol),' ',mass
  flush(OUTPUT_UNIT)
  !=================================



  allocate(Qact(nb_act))
  allocate(Qcart(nb_cart))

  IF (X2Q) THEN
    Qact(:) = 0.5_Rk
    write(OUTPUT_UNIT,*) 'Qact (initial values)',Qact
    CALL Qact_TO_cart(Qact,size(Qact),Qcart,size(Qcart))
    CALL cart_TO_Qact(Qact,size(Qact),Qcart,size(Qcart))
    write(OUTPUT_UNIT,*) 'Qact (from cart_TO_Qact)',Qact
  END IF

  Qact(:) = 0.5_Rk
  Qact(1) = 0.5_Rk + nb_eval*0.001_Rk
  write(OUTPUT_UNIT,*) 'Qact (final values)',Qact
  CALL Qact_TO_cart(Qact,size(Qact),Qcart,size(Qcart))
  write(OUTPUT_UNIT,*) 'Qcart (not recenter / COM)'
  DO i=1,nb_cart,3
    write(OUTPUT_UNIT,*) (i-1)/3+1,Qcart(i:i+2)
  END DO

  CALL Qact_TO_cartCOM(Qact,size(Qact),Qcart,size(Qcart))
  write(OUTPUT_UNIT,*) 'Qcart (recenter / COM)'
  DO i=1,nb_cart,3
    write(OUTPUT_UNIT,*) (i-1)/3+1,Qcart(i:i+2)
  END DO

  CALL Tnum_get_EckartRot(Qact,size(Qact),EckartRot)
  write(OUTPUT_UNIT,*) 'Eckart Rotation matrix'
  write(OUTPUT_UNIT,*) '1 ',EckartRot(:,1)
  write(OUTPUT_UNIT,*) '2 ',EckartRot(:,2)
  write(OUTPUT_UNIT,*) '3 ',EckartRot(:,3)

  deallocate(Qact)
  deallocate(Qcart)

  CALL date_and_time(TIME=time0,VALUES=Time_val0)
  write(OUTPUT_UNIT,*) 'Beginning time:',time0
  write(OUTPUT_UNIT,*) 'Beginning loop:',nb_eval
 !$OMP   PARALLEL &
 !$OMP   DEFAULT(NONE) &
 !$OMP   SHARED (nb_act,nb_cart,nb_eval) &
 !$OMP   PRIVATE(i,Qact,Qcart,EckartRot)

  allocate(Qact(nb_act))
  Qact(:) = 0.5_Rk
  allocate(Qcart(nb_cart))

 !$OMP   DO SCHEDULE(STATIC)
  DO i=1,nb_eval
    IF (mod(i,100) == 0) write(OUTPUT_UNIT,'(".")',advance='no')
    Qact(1) = 0.5_Rk + i*0.5_Rk/nb_eval
    CALL Qact_TO_cart(Qact,nb_act,Qcart,nb_cart)
    CALL Tnum_get_EckartRot(Qact,nb_act,EckartRot)
  END DO
 !$OMP   END DO

  deallocate(Qact)
  deallocate(Qcart)

 !$OMP   END PARALLEL
  write(OUTPUT_UNIT,*)
  write(OUTPUT_UNIT,*) 'END loop'
  CALL date_and_time(TIME=time1,VALUES=Time_val1)

  write(OUTPUT_UNIT,*) 'End time:',time1
  Time_val1 = Time_val1 - Time_val0
  !write(OUTPUT_UNIT,*) 'Delta time:',Time_val1
  write(OUTPUT_UNIT,*) 'Delta time:',Time_val1(5)*3600._Rk+Time_val1(6)*60._Rk+Time_val1(7)+Time_val1(7)/1000._Rk

CONTAINS
  SUBROUTINE ReadArguments(InputName,OutputName,X2Q,nb_eval)
    USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : INPUT_UNIT,OUTPUT_UNIT
    IMPLICIT NONE

    character (len=:), allocatable, intent(inout)  :: InputName,OutputName
    logical,                        intent(inout)  :: X2Q
    integer,                        intent(inout)  :: nb_eval

    character(len=:), allocatable :: arg,arg2
    integer :: long , i,n

    IF (mod(COMMAND_ARGUMENT_COUNT(),2) == 1) THEN
      write(OUTPUT_UNIT,*) 'Wrong number of arguments: ',COMMAND_ARGUMENT_COUNT()
      write(OUTPUT_UNIT,*) ' It MUST be an even number.'

      STOP 'ERROR in ReadArguments: Wrong number of arguments'
    END IF

    X2Q     = .TRUE.
    nb_eval = 10**4

    DO i=1,COMMAND_ARGUMENT_COUNT(),2
      CALL GET_COMMAND_ARGUMENT( NUMBER=i, LENGTH=long )
      allocate( character(len=long) :: arg )
      CALL GET_COMMAND_ARGUMENT( NUMBER=i, VALUE=arg )

      CALL GET_COMMAND_ARGUMENT( NUMBER=i+1, LENGTH=long )
      allocate( character(len=long) :: arg2 )
      CALL GET_COMMAND_ARGUMENT( NUMBER=i+1, VALUE=arg2 )

      SELECT CASE(arg)
      CASE("-o","-output")
        OutputName = arg2
      CASE("-i","-input")
        InputName = arg2
      CASE("-x2q","-X2Q")
        read(arg2,*) X2Q
      CASE("-nb_eval")
        read(arg2,*) nb_eval
      CASE Default
        write(OUTPUT_UNIT,*) 'Number of argument(s): ',COMMAND_ARGUMENT_COUNT()
        STOP 'ERROR in ReadArguments: no default argument'
      END SELECT

      deallocate(arg)
      deallocate(arg2)

    END DO
    IF (.NOT. allocated(OutputName)) OutputName = 'res_driver'
    IF (.NOT. allocated(InputName))  InputName  = 'dat_driver'

    write(OUTPUT_UNIT,*) 'OutputName:       ',OutputName
    write(OUTPUT_UNIT,*) 'InputName:        ',InputName
    write(OUTPUT_UNIT,*) 'X2Q (Cart_TO_Q):  ',X2Q
    write(OUTPUT_UNIT,*) 'nb_eval:          ',nb_eval

  END SUBROUTINE ReadArguments
END PROGRAM Main_TnumTana_FDriver
