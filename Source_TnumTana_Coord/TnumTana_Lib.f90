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
  USE mod_ActiveTransfo
  IMPLICIT NONE

  integer,           intent(in)     :: nb_act,nb_cart

  real (kind=Rkind), intent(in)     :: Qact(nb_act)
  real (kind=Rkind), intent(inout)  :: Qcart(nb_cart)

  real (kind=Rkind), allocatable     :: Qact_loc(:)


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
    CALL Finalize_TnumTana_Coord_PrimOp(para_Tnum,mole,PrimOp)

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

  IF (nb_act == mole%nb_var) THEN
    CALL sub_QactTOd0x(Qcart,Qact,mole,Gcenter=.FALSE.)
  ELSE IF (nb_act < mole%nb_var) THEN
    allocate(Qact_loc(mole%nb_var))
    CALL get_Qact0(Qact_loc,mole%ActiveTransfo)
    Qact_loc(1:nb_act) = Qact
    CALL sub_QactTOd0x(Qcart,Qact_loc,mole,Gcenter=.FALSE.)
    deallocate(Qact_loc)
  ELSE
    write(out_unitp,*) ' ERROR in ', name_sub
    write(out_unitp,*) ' nb_act is larger than mole%nb_var'
    write(out_unitp,*) ' nb_act     ',nb_act
    write(out_unitp,*) ' mole%nb_var',mole%nb_var
    STOP 'ERROR in Qact_TO_cart: nb_act is larger than mole%nb_var'
  END IF


END SUBROUTINE Qact_TO_cart
SUBROUTINE cart_TO_Qact(Qact,nb_act,Qcart,nb_cart)
  USE mod_system
  USE Module_ForTnumTana_Driver
  IMPLICIT NONE

  integer,           intent(in)     :: nb_act,nb_cart

  real (kind=Rkind), intent(in)     :: Qact(nb_act)
  real (kind=Rkind), intent(inout)  :: Qcart(nb_cart)


  character (len=*), parameter :: name_sub='cart_TO_Qact'

!===========================================================
!===========================================================
  !$OMP    CRITICAL (cart_TO_Qact_CRIT)
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
    CALL Finalize_TnumTana_Coord_PrimOp(para_Tnum,mole,PrimOp)

    IF (nb_act /= mole%nb_act .OR. nb_cart /= mole%ncart_act) THEN
       write(out_unitp,*) ' ERROR in ', name_sub
       write(out_unitp,*) ' nb_act is different from the Tnum one ',nb_act,mole%nb_act
       write(out_unitp,*) '    or '
       write(out_unitp,*) ' nb_cart is different from the Tnum one ',nb_cart,mole%ncart_act
       STOP
    END IF

  END IF
  !$OMP   END CRITICAL (cart_TO_Qact_CRIT)
!===========================================================
!===========================================================



  CALL sub_d0xTOQact(Qcart,Qact,mole)


END SUBROUTINE cart_TO_Qact
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
SUBROUTINE Tnum_Check_NMTransfo(Check)
  USE mod_system
  USE Module_ForTnumTana_Driver
  IMPLICIT NONE

  logical,           intent(inout)     :: Check

  character (len=*), parameter :: name_sub='Tnum_Check_NMTransfo'


  !$OMP    CRITICAL (Tnum_Check_NMTransfo_CRIT)
  Check = associated(mole%NMTransfo)
  !$OMP   END CRITICAL (Tnum_Check_NMTransfo_CRIT)

END SUBROUTINE Tnum_Check_NMTransfo
SUBROUTINE Tnum_Set_skip_NM(skip_NM_sub)
  USE mod_system
  USE Module_ForTnumTana_Driver
  IMPLICIT NONE

  integer,           intent(in)     :: skip_NM_sub

  character (len=*), parameter :: name_sub='Tnum_Set_skip_NM'


  !$OMP    CRITICAL (Tnum_Set_skip_NM_CRIT)
  IF (skip_NM_sub == 0 .OR. skip_NM_sub == 1) THEN
    skip_NM = skip_NM_sub
    IF (associated(mole%NMTransfo)) THEN
      mole%tab_Qtransfo(mole%itNM)%skip_transfo = (skip_NM == 1)
    END IF

  ELSE
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) ' Wrong skip_NM value. Possible values: 0 or 1'
    write(out_unitp,*) ' skip_NM',skip_NM_sub
    STOP 'ERROR in Tnum_Set_skip_NM: wrong skip_NM value!'
  END IF
  !$OMP   END CRITICAL (Tnum_Set_skip_NM_CRIT)

END SUBROUTINE Tnum_Set_skip_NM
SUBROUTINE Tnum_Set_k_Half(k_Half_sub)
  USE mod_system
  USE Module_ForTnumTana_Driver
  IMPLICIT NONE

  integer,           intent(in)     :: k_Half_sub

  character (len=*), parameter :: name_sub='Tnum_Set_k_Half'


  !$OMP    CRITICAL (Tnum_Set_k_Half_CRIT)
  IF (k_Half_sub == 0 .OR. k_Half_sub == 1) THEN
    k_Half = k_Half_sub
  ELSE
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) ' Wrong k_Half value. Possible values: 0 or 1'
    write(out_unitp,*) ' k_Half',k_Half_sub
    STOP 'ERROR in Tnum_Set_skip_NM: wrong k_Half value!'
  END IF
  !$OMP   END CRITICAL (Tnum_Set_k_Half_CRIT)

END SUBROUTINE Tnum_Set_k_Half
SUBROUTINE Tnum_Set_active_masses(active_masses,ncart_act)
  USE mod_system
  USE Module_ForTnumTana_Driver
  IMPLICIT NONE

  integer,           intent(in)        :: ncart_act
  integer,           intent(in)        :: active_masses(ncart_act)

  character (len=*), parameter :: name_sub='Tnum_Set_active_masses'

  CALL Check_TnumInit(name_sub)

  IF (ncart_act == mole%ncart_act) THEN
    mole%active_masses(:) = active_masses
  ELSE
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) ' Wrong ncart_act value. It MUST be: ',mole%ncart_act
    write(out_unitp,*) ' ncart_act',ncart_act
    STOP 'ERROR in Tnum_Set_active_masses: Wrong ncart_act value!'
  END IF

END SUBROUTINE Tnum_Set_active_masses
SUBROUTINE Tnum_get_ndimG(ndimG)
  USE mod_system
  USE Module_ForTnumTana_Driver
  IMPLICIT NONE

  integer,           intent(inout)     :: ndimG

  character (len=*), parameter :: name_sub='Tnum_get_ndimG'

  CALL Check_TnumInit(name_sub)

  ndimG = mole%ndimG

END SUBROUTINE Tnum_get_ndimG
SUBROUTINE Tnum_get_Z(Z,nb_at)
  USE mod_system
  USE Module_ForTnumTana_Driver
  IMPLICIT NONE

  integer,           intent(in)        :: nb_at
  integer,           intent(inout)     :: Z(nb_at)

  character (len=*), parameter :: name_sub='Tnum_get_Z'

  CALL Check_TnumInit(name_sub)

  IF (nb_at <= mole%nat .AND. nb_at >= 1) THEN
    Z(:) = mole%Z(1:nb_at)
  ELSE
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) ' Wrong nb_at value. It has to be in the range: [1:',mole%nat,']'
    write(out_unitp,*) ' nb_at',nb_at
    STOP 'ERROR in Tnum_get_Z: Wrong nb_at value!'
  END IF

END SUBROUTINE Tnum_get_Z
SUBROUTINE Tnum_get_Qact0(Qact0,nb_act)
  USE mod_system
  USE Module_ForTnumTana_Driver
  USE mod_ActiveTransfo
  IMPLICIT NONE

  integer,           intent(in)        :: nb_act
  real (kind=Rkind), intent(inout)     :: Qact0(nb_act)

  character (len=*), parameter :: name_sub='Tnum_get_Qact0'

  CALL Check_TnumInit(name_sub)

  IF (nb_act == mole%nb_act) THEN
    CALL get_Qact0(Qact0,mole%ActiveTransfo)
  ELSE
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) ' Wrong nb_act value. It must be equal to: ',mole%nb_act
    write(out_unitp,*) ' nb_act',nb_act
    STOP 'ERROR in Tnum_get_Qact0:  Wrong nb_act value!'
  END IF

END SUBROUTINE Tnum_get_Qact0
SUBROUTINE Tnum_get_mass(mass,Z,A,name)
  USE mod_system
  USE Module_ForTnumTana_Driver
  IMPLICIT NONE

  integer,            intent(inout)     :: Z,A
  real (kind=Rkind),  intent(inout)     :: mass
  character (len=10), intent(in)        :: name


  integer     :: err_mass


  character (len=*), parameter :: name_sub='Tnum_get_mass'

  CALL Check_TnumInit(name_sub)

  IF (len_trim(adjustl(name)) == 0) THEN
    mass = get_mass_Tnum(const_phys%mendeleev,Z,A,err_mass=err_mass)
  ELSE
    mass = get_mass_Tnum(const_phys%mendeleev,Z,A,name,err_mass)
  END IF

  IF (err_mass /= 0) THEN
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) ' The mass cannot be find!'
    write(out_unitp,*) ' Z,A,name',Z,A,trim(name)
    STOP 'ERROR in Tnum_get_mass: The mass cannot be find!'
  END IF

END SUBROUTINE Tnum_get_mass
SUBROUTINE Tnum_get_EigNM_ForPVSCF(EigNM,nb_NM)
  USE mod_system
  USE mod_Constant
  USE Module_ForTnumTana_Driver
  IMPLICIT NONE

  integer,           intent(in)        :: nb_NM
  real (kind=Rkind), intent(inout)     :: EigNM(nb_NM)

  character (len=*), parameter :: name_sub='Tnum_get_EigNM_ForPVSCF'

  CALL Check_TnumInit(name_sub)

  IF (nb_NM /= mole%nb_act) THEN
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) ' Wrong nb_NM value. It must be equal to: ',mole%nb_act
    write(out_unitp,*) ' nb_NM',nb_NM
    STOP 'ERROR in Tnum_get_EigNM_ForPVSCF:  Wrong nb_NM value!'
  END IF
  IF (.NOT. associated(mole%NMTransfo%d0eh)) THEN
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) ' mole%NMTransfo%d0eh(:) is NOT associated!'
    write(out_unitp,*) ' => You MUST call InitTnum3_NM_TO_LinearTransfo'
    STOP 'ERROR in Tnum_get_EigNM_ForPVSCF: d0eh is not associated!'
  END IF

  EigNM(:) = mole%NMTransfo%d0eh(:)**2 * const_phys%inv_Name/TEN**3

END SUBROUTINE Tnum_get_EigNM_ForPVSCF
SUBROUTINE Tnum_get_GG(Qact,nb_act,GG,ndimG,def)
  USE mod_system
  USE Module_ForTnumTana_Driver
  USE mod_ActiveTransfo
  USE mod_dnGG_dng
  IMPLICIT NONE

  integer,           intent(in)        :: ndimG,nb_act
  logical,           intent(in)        :: def
  real (kind=Rkind), intent(in)        :: Qact(nb_act)
  real (kind=Rkind), intent(inout)     :: GG(ndimG,ndimG)

  real (kind=Rkind), allocatable       :: Qact_loc(:)

  character (len=*), parameter :: name_sub='Tnum_get_GG'

  CALL Check_TnumInit(name_sub)

  IF (nb_act /= mole%nb_act) THEN
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) ' Wrong nb_act value. It must be equal to: ',mole%nb_act
    write(out_unitp,*) ' nb_act',nb_act
    STOP 'ERROR in Tnum_get_GG:  Wrong nb_act value!'
  END IF

  IF (ndimG /= mole%ndimG) THEN
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) ' Wrong ndimG value. It must be equal to: ',mole%ndimG
    write(out_unitp,*) ' ndimG',ndimG
    STOP 'ERROR in Tnum_get_GG:  Wrong ndimG value!'
  END IF

  IF (nb_act <= mole%nb_var) THEN
    allocate(Qact_loc(mole%nb_var))
    CALL get_Qact0(Qact_loc,mole%ActiveTransfo)
    Qact_loc(1:nb_act) = Qact
  ELSE
    write(out_unitp,*) ' ERROR in ', name_sub
    write(out_unitp,*) ' nb_act is larger than mole%nb_var'
    write(out_unitp,*) ' nb_act     ',nb_act
    write(out_unitp,*) ' mole%nb_var',mole%nb_var
    STOP 'ERROR in Tnum_get_GG: nb_act is larger than mole%nb_var'
  END IF

  CALL get_d0GG(Qact_loc,para_Tnum,mole,GG,def=def)

  deallocate(Qact_loc)

END SUBROUTINE Tnum_get_GG
SUBROUTINE InitTnum3_NM_TO_LinearTransfo(Qact,nb_act,Hess,nb_cart)
  USE Module_ForTnumTana_Driver
  USE mod_ActiveTransfo
  USE mod_PrimOp
  IMPLICIT NONE

  integer,           intent(in)     :: nb_act,nb_cart

  real (kind=Rkind)                 :: Qact(nb_act),Hess(nb_cart,nb_cart)

  real (kind=Rkind), allocatable    :: Qact_loc(:)

  character (len=*), parameter :: name_sub='InitTnum3_NM_TO_LinearTransfo'


    CALL Check_TnumInit(name_sub)

    CALL Tnum_Set_skip_NM(0) ! To be able to initialize the Normal Modes
    print_level=0

    IF (nb_act <= mole%nb_var) THEN
      allocate(Qact_loc(mole%nb_var))
      CALL get_Qact0(Qact_loc,mole%ActiveTransfo)
      Qact_loc(1:nb_act) = Qact
    ELSE
      write(out_unitp,*) ' ERROR in ', name_sub
      write(out_unitp,*) ' nb_act is larger than mole%nb_var'
      write(out_unitp,*) ' nb_act     ',nb_act
      write(out_unitp,*) ' mole%nb_var',mole%nb_var
      STOP 'ERROR in InitTnum3_NM_TO_LinearTransfo: nb_act is larger than mole%nb_var'
    END IF

    CALL calc3_NM_TO_sym(Qact_loc,mole,para_Tnum,PrimOp,Hess,.TRUE.)

    deallocate(Qact_loc)

END SUBROUTINE InitTnum3_NM_TO_LinearTransfo

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
    IF (associated(mole%NMTransfo)) THEN
      mole%tab_Qtransfo(mole%itNM)%skip_transfo = (skip_NM == 1) ! the NM are not calculated
      IF (k_Half == 0 .OR. k_Half == 1) mole%NMTransfo%k_Half = (k_Half == 1)
    END IF
    CALL Finalize_TnumTana_Coord_PrimOp(para_Tnum,mole,PrimOp)

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
