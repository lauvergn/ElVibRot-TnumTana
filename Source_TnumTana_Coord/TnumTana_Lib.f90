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
SUBROUTINE Qact_TO_cart(Qact,nb_act,Qcart,nb_cart)
  USE mod_system
  USE Module_ForTnumTana_Driver
  IMPLICIT NONE

  integer,           intent(in)     :: nb_act,nb_cart

  real (kind=Rkind), intent(in)     :: Qact(nb_act)
  real (kind=Rkind), intent(inout)  :: Qcart(nb_cart)


  character (len=*), parameter :: name_sub='Qact_TO_cart'

!===========================================================
!===========================================================
  !$OMP    CRITICAL (Qact_TO_cart_CRIT)
  IF (Init == 0) THEN
    Init = 1
    CALL versionEVRT(.TRUE.)
    print_level=-1
    !-----------------------------------------------------------------
    !     - read the coordinate transformations :
    !     -   zmatrix, polysperical, bunch...
    !     ------------------------------------------------------------
    CALL Read_CoordType(mole,para_Tnum,const_phys)
    !     ------------------------------------------------------------
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    !     - read coordinate values -----------------------------------
    !     ------------------------------------------------------------
    CALL read_RefGeom(mole,para_Tnum)
    !     ------------------------------------------------------------
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    !     ---- TO finalize the coordinates (NM) and the KEO ----------
    !     ------------------------------------------------------------
    para_Tnum%Tana=.FALSE.
    CALL Finalize_TnumTana_Coord_PrimOp(para_Tnum,mole,para_PES)

    IF (nb_act /= mole%nb_act .OR. nb_cart /= mole%ncart_act) THEN
       write(out_unitp,*) ' ERROR in ', name_sub
       write(out_unitp,*) ' nb_act is different from the Tnum one ',nb_act,mole%nb_act
       write(out_unitp,*) '    or '
       write(out_unitp,*) ' nb_cart is different from the Tnum one ',nb_cart,mole%ncart_act
       STOP
    END IF

  END IF
  !$OMP   END CRITICAL (Qact_TO_cart_CRIT)

!===========================================================
!===========================================================

  CALL sub_QactTOd0x(Qcart,Qact,mole,Gcenter=.FALSE.)


END SUBROUTINE Qact_TO_cart
SUBROUTINE Init_InputUnit_Driver(InputUnit)
  USE mod_system
  IMPLICIT NONE

  integer,           intent(in)     :: InputUnit

  character (len=*), parameter :: name_sub='Init_InputUnit_Driver'


  !$OMP    CRITICAL (Init_InputUnit_Driver_CRIT)
  in_unitp  = InputUnit
  !$OMP   END CRITICAL (Init_InputUnit_Driver_CRIT)

END SUBROUTINE Init_InputUnit_Driver
SUBROUTINE Init_OutputUnit_Driver(OutputUnit)
  USE mod_system
  IMPLICIT NONE

  integer,           intent(in)     :: OutputUnit

  character (len=*), parameter :: name_sub='Init_OutputUnit_Driver'


  !$OMP    CRITICAL (Init_OutputUnit_Driver_CRIT)
  out_unitp = OutputUnit
  !$OMP   END CRITICAL (Init_OutputUnit_Driver_CRIT)

END SUBROUTINE Init_OutputUnit_Driver
SUBROUTINE Init_TnumTana_FOR_Driver(nb_act,nb_cart,init_sub)
  USE Module_ForTnumTana_Driver
  IMPLICIT NONE

  integer,           intent(inout)     :: nb_act,nb_cart

  integer,           intent(inout)     :: init_sub

  character (len=*), parameter :: name_sub='Init_TnumTana_FOR_Driver'


  !$OMP    CRITICAL (Init_TnumTana_FOR_Driver_CRIT)
  IF (Init == 0 .OR. init_sub == 0) THEN
    init     = 1
    init_sub = 1

    CALL versionEVRT(.TRUE.)
    print_level=-1
    !-----------------------------------------------------------------
    !     - read the coordinate transformations :
    !     -   zmatrix, polysperical, bunch...
    !     ------------------------------------------------------------
    CALL Read_CoordType(mole,para_Tnum,const_phys)
    !     ------------------------------------------------------------
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    !     - read coordinate values -----------------------------------
    !     ------------------------------------------------------------
    CALL read_RefGeom(mole,para_Tnum)
    !     ------------------------------------------------------------
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    !     ---- TO finalize the coordinates (NM) and the KEO ----------
    !     ------------------------------------------------------------
    CALL Finalize_TnumTana_Coord_PrimOp(para_Tnum,mole,para_PES)

  END IF

  nb_act  = mole%nb_act
  nb_cart = mole%ncart_act
  !$OMP   END CRITICAL (Init_TnumTana_FOR_Driver_CRIT)


END SUBROUTINE Init_TnumTana_FOR_Driver
SUBROUTINE Init_TnumTana_FOR_Driver_FOR_c(nb_act,nb_cart,init_sub)  BIND(C, name="Init_TnumTana_FOR_Driver_FOR_c")
  USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT
  IMPLICIT NONE

  integer (kind=C_INT), intent(inout)     :: nb_act,nb_cart
  integer (kind=C_INT), intent(inout)     :: init_sub



  integer               :: nb_act_loc,nb_cart_loc
  integer               :: init_sub_loc

  character (len=*), parameter :: name_sub='Init_TnumTana_FOR_Driver_FOR_c'



  CALL Init_TnumTana_FOR_Driver(nb_act_loc,nb_cart_loc,init_sub_loc)


  nb_act   = nb_act_loc
  nb_cart  = nb_cart_loc
  init_sub = init_sub_loc


END SUBROUTINE Init_TnumTana_FOR_Driver_FOR_c
SUBROUTINE Qact_TO_Qcart_TnumTanaDriver_FOR_c(Qact,nb_act,Qcart,nb_cart) BIND(C, name="Qact_TO_Qcart_TnumTanaDriver_FOR_c")
  USE, INTRINSIC :: ISO_C_BINDING,             ONLY : C_INT,C_DOUBLE
  USE            :: mod_system
  USE            :: Module_ForTnumTana_Driver, ONLY : mole,init,sub_QactTOd0x
  IMPLICIT NONE

  integer (kind=C_INT), intent(in)     :: nb_act,nb_cart

  real (kind=C_DOUBLE), intent(in)     :: Qact(nb_act)
  real (kind=C_DOUBLE), intent(inout)  :: Qcart(nb_cart)


  !- local parameters for para_Tnum -----------------------
  real (kind=Rkind)      :: Qact_loc(nb_act)
  real (kind=Rkind)      :: Qcart_loc(nb_cart)


  character (len=*), parameter :: name_sub='Qact_TO_Qcart_TnumTanaDriver_FOR_c'


  IF (Init == 0) THEN
    write(out_unitp,*) ' ERROR in ', name_sub
    write(out_unitp,*) '   The intialization IS NOT done!'
    write(out_unitp,*) ' First, you MUST call Init_TnumTana_FOR_Driver_FOR_c'
    STOP
  END IF
  IF (nb_act /= mole%nb_act .OR. nb_cart /= mole%ncart_act) THEN
     write(out_unitp,*) ' ERROR in ', name_sub
     write(out_unitp,*) ' nb_act is different from the Tnum one ',nb_act,mole%nb_act
     write(out_unitp,*) '    or '
     write(out_unitp,*) ' nb_cart is different from the Tnum one ',nb_cart,mole%ncart_act
     STOP
  END IF

  Qact_loc(:)  = Qact
  CALL sub_QactTOd0x(Qcart_loc,Qact_loc,mole,Gcenter=.FALSE.)
  Qcart(:)     = Qcart_loc

END SUBROUTINE Qact_TO_Qcart_TnumTanaDriver_FOR_c
SUBROUTINE Qcart_TO_Qact_TnumTanaDriver_FOR_c(Qact,nb_act,Qcart,nb_cart) BIND(C, name="Qcart_TO_Qact_TnumTanaDriver_FOR_c")
  USE, INTRINSIC :: ISO_C_BINDING,             ONLY : C_INT,C_DOUBLE
  USE            :: mod_system
  USE            :: Module_ForTnumTana_Driver, ONLY : mole,Init,sub_d0xTOQact
  IMPLICIT NONE

  integer (kind=C_INT), intent(in)     :: nb_act,nb_cart

  real (kind=C_DOUBLE), intent(inout)  :: Qact(nb_act)
  real (kind=C_DOUBLE), intent(in)     :: Qcart(nb_cart)


  !- local parameters for para_Tnum -----------------------
  real (kind=Rkind)      :: Qact_loc(nb_act)
  real (kind=Rkind)      :: Qcart_loc(nb_cart)


  character (len=*), parameter :: name_sub='Qcart_TO_Qact_TnumTanaDriver_FOR_c'


  IF (Init == 0) THEN
    write(out_unitp,*) ' ERROR in ', name_sub
    write(out_unitp,*) '   The intialization IS NOT done!'
    write(out_unitp,*) ' First, you MUST call Init_TnumTana_FOR_Driver_FOR_c'
    STOP
  END IF
  IF (nb_act /= mole%nb_act .OR. nb_cart /= mole%ncart_act) THEN
     write(out_unitp,*) ' ERROR in ', name_sub
     write(out_unitp,*) ' nb_act is different from the Tnum one ',nb_act,mole%nb_act
     write(out_unitp,*) '    or '
     write(out_unitp,*) ' nb_cart is different from the Tnum one ',nb_cart,mole%ncart_act
     STOP
  END IF

  Qcart_loc(:) = Qcart(:)
  CALL sub_d0xTOQact(Qcart_loc,Qact_loc,mole)
  Qact(:)      = Qact_loc(:)


END SUBROUTINE Qcart_TO_Qact_TnumTanaDriver_FOR_c
