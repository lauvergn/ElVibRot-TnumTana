!
!=======================================================================
!    calcN_op: calculation of the potential and scalar operator matrices
!     - mat_V(nb_be,nb_be):   potential (real part)
!     - mat_imV(nb_be,nb_be): potential (imaginary part)
!     - mat_ScalOp(nb_be,nb_be,nb_ScalOp): scalar operators (dipole ...)
!    nb_be:     nb of diabatic electronic states
!    nb_ScalOp: nb of scalar operators
!
!    Qpot are the coordinates (active, dynamic Cartesian ... ones)
!    nb_Qpot: nb of coordinates
!
!     mole: coordinate definition
!=======================================================================
      SUBROUTINE calcN_op(mat_V,mat_imV,mat_ScalOp,nb_be,nb_ScalOp,     &
                          Qpot,nb_Qpot,mole,calc_ScalOp,pot_cplx)
      USE mod_system
      USE mod_Tnum
      IMPLICIT NONE

      !----- for the coordinate definition------------------------------
      TYPE (CoordType)    :: mole

      integer           :: nb_be,nb_ScalOp,nb_Qpot
      logical           :: calc_ScalOp,pot_cplx
      real (kind=Rkind) :: mat_V(nb_be,nb_be)
      real (kind=Rkind) :: mat_imV(nb_be,nb_be)
      real (kind=Rkind) :: mat_ScalOp(nb_be,nb_be,nb_ScalOp)
      real (kind=Rkind) :: Qpot(nb_Qpot)


      !-----------------------------------------------------------------
      mat_V(:,:) = ZERO
      IF (pot_cplx) THEN
        mat_imV(:,:) = ZERO
      END IF
      IF (calc_ScalOp) THEN
          mat_ScalOp(:,:,:) = ZERO
      END IF
      !-----------------------------------------------------------------

      STOP 'STOP in calcN_op: you have to set-up the potential!!'

      END SUBROUTINE calcN_op

!=======================================================================
!    fonction pot_rest(x)
!=======================================================================
      FUNCTION pot_rest(Qact,Delta_Qact,nb_inact2n)
      USE mod_system
      IMPLICIT NONE
      real (kind=Rkind) :: pot_rest


      real (kind=Rkind) :: Qact(1)
      integer           :: nb_inact2n
      real (kind=Rkind) :: Delta_Qact(nb_inact2n)

      !-----------------------------------------------------------------
      pot_rest = ZERO
      !-----------------------------------------------------------------


      !-----------------------------------------------------------------
      STOP 'The function pot_rest MUST be make'
      !-----------------------------------------------------------------


      END FUNCTION pot_rest

!=======================================================================
!    sub hessian
!=======================================================================
      SUBROUTINE sub_hessian(h)
      USE mod_system
      IMPLICIT NONE

       real (kind=Rkind) :: h

      !-----------------------------------------------------------------
       h = ZERO
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      STOP 'The subroutine sub_hessian MUST be make'
      !-----------------------------------------------------------------


      END SUBROUTINE sub_hessian
!=======================================================================
!     analytical gradient along the path
!=======================================================================
      SUBROUTINE d0d1d2_g(d0g,d1g,d2g,Qdyn,mole,deriv,num,step)
      USE mod_system
      USE mod_Tnum
      IMPLICIT NONE

      !----- for the CoordType and Tnum ---------------------------------
      TYPE (CoordType)    :: mole

      real (kind=Rkind) :: d0g(mole%nb_inact2n)
      real (kind=Rkind) :: d1g(mole%nb_inact2n,mole%nb_act1)
      real (kind=Rkind) :: d2g(mole%nb_inact2n,mole%nb_act1,mole%nb_act1)

      real (kind=Rkind) :: Qdyn(mole%nb_var)
      real (kind=Rkind) :: step
      logical           :: deriv,num

      real (kind=Rkind) :: Qact1(mole%nb_act1)

      !----- for debuging ----------------------------------------------
      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      character (len=*), parameter :: name_sub='d0d1d2_g'
      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_var',mole%nb_var
        write(out_unitp,*) 'nb_act1',mole%nb_act1
        write(out_unitp,*) 'nb_inact22,nb_inact21',mole%nb_inact22,mole%nb_inact21
        write(out_unitp,*) 'nb_inact2n',mole%nb_inact2n
        write(out_unitp,*) 'deriv',deriv
      END IF
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      Qact1(:) = Qdyn(mole%liste_QactTOQdyn(1:mole%nb_act1))
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      d0g(:)     = ZERO
      d1g(:,:)   = ZERO
      d2g(:,:,:) = ZERO
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      STOP 'The subroutine d0d1d2_g MUST be make'
      !-----------------------------------------------------------------


      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'd0g at Qact:',Qact1
        write(out_unitp,*) d0g(:)
        write(out_unitp,*) 'END ',name_sub
      END IF
      !-----------------------------------------------------------------

      END SUBROUTINE d0d1d2_g
!=======================================================================
!     analytical hessian along the path (only d0h is used!!)
!=======================================================================
      SUBROUTINE d0d1d2_h(d0h,d1h,d2h,Qdyn,mole,deriv,num,step)
      USE mod_system
      USE mod_Tnum
      IMPLICIT NONE

      !----- for the CoordType and Tnum ----------------------------------
      TYPE (CoordType)    :: mole

      real (kind=Rkind) :: Qdyn(mole%nb_var)


      real (kind=Rkind) :: d0h(mole%nb_inact2n,mole%nb_inact2n)
      real (kind=Rkind) :: d1h(mole%nb_inact2n,mole%nb_inact2n)
      real (kind=Rkind) :: d2h(mole%nb_inact2n,mole%nb_inact2n)

      real (kind=Rkind) :: step
      logical           :: deriv,num


      real (kind=Rkind) :: Qact1(mole%nb_act1)

      !----- for debuging ----------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='d0d1d2_h'
      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_var',mole%nb_var
        write(out_unitp,*) 'nb_act1',mole%nb_act1
        write(out_unitp,*) 'nb_inact22,nb_inact21',mole%nb_inact22,mole%nb_inact21
        write(out_unitp,*) 'nb_inact2n',mole%nb_inact2n
      END IF
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      Qact1(:) = Qdyn(mole%liste_QactTOQdyn(1:mole%nb_act1))
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      d0h(:,:) = ZERO
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      STOP 'The subroutine d0d1d2_h MUST be make'
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Qact1',Qact1
        write(out_unitp,*) 'd0h at Qact1'
        CALL Write_Mat(d0h,6,4)
        write(out_unitp,*) 'END ',name_sub
      END IF
      !-----------------------------------------------------------------

      END SUBROUTINE d0d1d2_h
!=======================================================================
!     analytical derivative (dnQflex : Qflex Qflex' Qflex" Qflex'") calculation
!     for the variable iq
!=======================================================================
      SUBROUTINE calc_dnQflex(iq,dnQflex,Qact,nb_act,nderiv,it)
      USE mod_system
      USE mod_dnSVM
      IMPLICIT NONE

      integer           :: iq,nb_act
      real (kind=Rkind) :: Qact(nb_act)
      integer           :: nderiv,it
      TYPE (Type_dnS)   :: dnQflex


      ! for debuging ---------------------------------------------------
      character (len=*), parameter :: name_sub='dnQflex'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      ! for debuging ---------------------------------------------------


      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_act',nb_act
        write(out_unitp,*) 'iq',iq
      END IF
      !-----------------------------------------------------------------


      !-----------------------------------------------------------------
      CALL sub_ZERO_TO_dnS(dnQflex)

      ! Zero order dervivative
      dnQflex%d0 = ZERO
      ! First order dervivatives
      dnQflex%d1(:) = ZERO
      ! Second order dervivatives
      dnQflex%d2(:,:) = ZERO
      ! Third order dervivatives
      dnQflex%d3(:,:,:) = ZERO
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      STOP 'The subroutine calc_dnQflex MUST be make'
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'dnQflex : ',Qact
        CALL write_dnS(dnQflex,nderiv)
        write(out_unitp,*) 'END ',name_sub
      END IF
      !-----------------------------------------------------------------

       END SUBROUTINE calc_dnQflex
!================================================================
!     analytical derivative (Qopt Qopt' Qopt" Qopt'") calculation
!     for the variable i_Qdyn
!     derivative with respect to Qact1(:) coordinates
!================================================================
      SUBROUTINE d0d1d2d3_Qeq(i_Qdyn,d0Qopt,d1Qopt,d2Qopt,d3Qopt,Qdyn,mole,nderiv)
      USE mod_system
      USE mod_Tnum
      IMPLICIT NONE

      !----- for the CoordType and Tnum ----------------------------------
      TYPE (CoordType)    :: mole

      integer           :: i_Qdyn

      real (kind=Rkind) :: Qdyn(mole%nb_var)

      integer           :: nderiv

      real (kind=Rkind) :: d0Qopt
      real (kind=Rkind) :: d1Qopt(mole%nb_act1)
      real (kind=Rkind) :: d2Qopt(mole%nb_act1,mole%nb_act1)
      real (kind=Rkind) :: d3Qopt(mole%nb_act1,mole%nb_act1,mole%nb_act1)


      !local variables
      real (kind=Rkind) :: Qact1(mole%nb_act1)


      !----- for debuging ----------------------------------------------
      character (len=*), parameter :: name_sub='d0d1d2d3_Qeq'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      !----- for debuging ----------------------------------------------

      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_inact20,nb_act',mole%nb_inact20,mole%nb_act
        write(out_unitp,*) 'nb_var',mole%nb_var
        write(out_unitp,*) 'i_Qdyn',i_Qdyn
      END IF
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      Qact1(:) = Qdyn(mole%liste_QactTOQdyn(1:mole%nb_act1))
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      d0Qopt        = ZERO
      d1Qopt(:)     = ZERO
      d2Qopt(:,:)   = ZERO
      d3Qopt(:,:,:) = ZERO
      !-----------------------------------------------------------------

      STOP 'The subroutine d0d1d2d3_Qeq MUST be make'


      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Qact1 : ',Qact1
        write(out_unitp,*) 'd0Qopt : ',d0Qopt
        write(out_unitp,*) 'd1Qopt : ',d1Qopt
        write(out_unitp,*) 'd2Qopt : ',d2Qopt
        write(out_unitp,*) 'd3Qopt : ',d3Qopt
        write(out_unitp,*) 'END ',name_sub
      END IF
      !-----------------------------------------------------------------

      END SUBROUTINE d0d1d2d3_Qeq
