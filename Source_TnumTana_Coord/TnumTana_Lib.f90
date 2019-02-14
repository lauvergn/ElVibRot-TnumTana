      SUBROUTINE Qact_TO_cart(Qact,nb_act,Qcart,nb_cart)
      USE mod_system
      USE mod_dnSVM
      USE mod_Constant
      USE mod_Coord_KEO
      USE mod_PrimOp
      IMPLICIT NONE

      integer,           intent(in)     :: nb_act,nb_cart

      real (kind=Rkind), intent(in)     :: Qact(nb_act)
      real (kind=Rkind), intent(inout)  :: Qcart(nb_cart)


      !- parameters for para_Tnum -----------------------
      TYPE (constant),  save :: const_phys
      TYPE (zmatrix),   save :: mole
      TYPE (Tnum),      save :: para_Tnum
      TYPE (param_PES), save :: para_PES

      logical,          save :: begin=.TRUE.

      character (len=*), parameter :: name_sub='Qact_TO_cart'

!===========================================================
!===========================================================
  !$OMP    CRITICAL (Qact_TO_cart_CRIT)
  IF (begin) THEN
    begin = .FALSE.
    CALL versionEVRT(.TRUE.)
    print_level=-1
    !-----------------------------------------------------------------
    !     - read the coordinate transformations :
    !     -   zmatrix, polysperical, bunch...
    !     ------------------------------------------------------------
    CALL Read_mole(mole,para_Tnum,const_phys)
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
    CALL Finalyze_TnumTana_Coord_PrimOp(para_Tnum,mole,para_PES)

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
