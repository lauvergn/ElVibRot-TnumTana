!c
!C================================================================
!C    calc_Op : calculation of the potential and dipolar matrices
!c    mat_V(nb_be,nb_be) and mat_dip(nb_be,nb_be,3)
!c    nb_be : nb of elctronic surface
!c    Q are the coordinates in active order or syl order
!c    dipolar calculation if calc_dip = T
!C================================================================
      SUBROUTINE calcN_op(mat_V,mat_imV,mat_ScalOp,nb_be,nb_ScalOp, &
                        Qact,nb_var,mole,calc_ScalOp,pot_cplx)

      USE mod_system
      USE mod_Tnum
      IMPLICIT NONE

!c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

      integer           :: nb_be,nb_ScalOp,nb_var
      logical           :: calc_ScalOp,pot_cplx
      real (kind=Rkind) :: mat_V(nb_be,nb_be),mat_imV(nb_be,nb_be)
      real (kind=Rkind) :: mat_ScalOp(nb_be,nb_be,nb_ScalOp)
      real (kind=Rkind) :: Qact(nb_var)

      real (kind=Rkind) :: im_pot0
      real (kind=Rkind) :: ScalOp(nb_ScalOp)

      
      !real (kind=Rkind), parameter :: b     = ZERO
      real (kind=Rkind), parameter :: b     = ONE
      real (kind=Rkind), parameter :: alpha = 3._Rkind
      !real (kind=Rkind), parameter :: v0    = 0.001822534_Rkind
      real (kind=Rkind), parameter :: v0    = 0.0018_Rkind
      real (kind=Rkind), parameter :: zmu   = 2000._Rkind


      IF (nb_be == 1 ) THEN
         mat_V(1,1) = ZERO
        !IF (Qact(1) > ZERO) mat_V(1,1) = TEN
         CALL pot_bottleneck_nd (Qact(1:mole%nb_act1), mole%nb_act1, mat_V(1,1), alpha, V0, b, zmu)
        !mat_V(1,1) = half * zmu * V0**2 * sum(Qact(2:mole%nb_act1)**2)

        !write(6,*) 'Qact,mat_V',Qact(1:mole%nb_act1),mat_V(1,1)
        IF (pot_cplx) mat_imV(1,1) = im_pot0(Qact(mole%nb_act1))
        IF (calc_ScalOp) THEN
          CALL sub_ScalarOp(ScalOp,nb_ScalOp,Qact(1:mole%nb_act1),mole)
          mat_ScalOp(1,1,:) = ScalOp(:)
        END IF
      ELSE
        STOP 'nb_be > 1'
      END IF

      END
!C================================================================
!C    fonction pot_rest(x)
!C================================================================
      FUNCTION pot_rest(Qact,Delta_Qact,nb_inact2n)
      USE mod_system
      IMPLICIT NONE
      real(kind=Rkind) :: pot_rest

       real(kind=Rkind) Qact(1)
       integer nb_inact2n
       real(kind=Rkind) Delta_Qact(nb_inact2n)

       pot_rest = ZERO

       END
!C================================================================
!C    fonction im_pot0(x)
!C================================================================
      FUNCTION im_pot0(Q)
      USE mod_system
      IMPLICIT NONE
      real(kind=Rkind) :: im_pot0

      real(kind=Rkind) Q(3)

       real(kind=Rkind), parameter :: a=-0.01,Q0=2.5

       IF (Q(3) > Q0) THEN
         im_pot0 = a*(Q(3) - Q0)
       ELSE
         im_pot0 = 0.d0
       END IF

       RETURN
       END
!C================================================================
!C    sub hessian
!C================================================================
      SUBROUTINE sub_hessian (h)
      USE mod_system
      IMPLICIT NONE

       real(kind=Rkind) h

       h = ZERO


      END
!C
!C================================================================
!C    fonction pot0(x) 1 D (avec x=cos(theta))
!C    pour une tri atomique en jacobie
!C================================================================
      SUBROUTINE sub_ScalarOp(ScalOp,nb_ScalOp,Q,mole)
      USE mod_system
      USE mod_Tnum
      IMPLICIT NONE

!C----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

       integer nb_Q,nb_ScalOp
       real(kind=Rkind) Q(mole%nb_act1)
       real(kind=Rkind) ScalOp(nb_ScalOp)
       integer :: i

       
       !real(kind=Rkind), parameter :: x0 = 20._Rkind
       !real(kind=Rkind), parameter :: x  = 25._Rkind
       !real(kind=Rkind), parameter :: A   = 0.006_Rkind

       real(kind=Rkind), parameter :: x0 = 9._Rkind
       real(kind=Rkind), parameter :: x  = 17._Rkind
       real(kind=Rkind), parameter :: A   = 0.003_Rkind

       IF (nb_ScalOp /= 2) STOP 'WRONG nb_ScalOp value. It MUST be 2'

       IF (Q(1) < -x0) THEN
         ScalOp(1) = A * ((-x0-Q(1))/(x-x0))**4
       ELSE
         ScalOp(1) = ZERO
       END IF
       IF (Q(1) > x0) THEN
         ScalOp(2) = A * ((x0-Q(1))/(x-x0))**4
       ELSE
         ScalOp(2) = ZERO
       END IF

       !write(6,*) 'Q,ScalOp',Q,ScalOp

       END
!C================================================================
!C    analytical gradient along the path
!C================================================================
      SUBROUTINE d0d1d2_g(d0g,d1g,d2g,Qdyn,mole,deriv,num,step)
      USE mod_system
      USE mod_Tnum
      IMPLICIT NONE

!C----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

      real (kind=Rkind) :: d0g(mole%nb_inact2n)
      real (kind=Rkind) :: d1g(mole%nb_inact2n,mole%nb_act1)
      real (kind=Rkind) :: d2g(mole%nb_inact2n,mole%nb_act1,mole%nb_act1)

      real (kind=Rkind) :: Qdyn(mole%nb_var)
      real (kind=Rkind) :: step
      logical       :: deriv,num

      real (kind=Rkind) :: Qact(mole%nb_act1)

!C----- for debuging ----------------------------------
!C     logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
!C---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING d0d1d2_g'
        write(out_unitp,*) 'nb_var',mole%nb_var
        write(out_unitp,*) 'nb_act1',mole%nb_act1
        write(out_unitp,*) 'nb_inact22,nb_inact21',mole%nb_inact22,mole%nb_inact21
        write(out_unitp,*) 'nb_inact2n',mole%nb_inact2n
        write(out_unitp,*) 'deriv',deriv
      END IF

!C---------------------------------------------------------------------
      Qact(1) = Qdyn(mole%liste_QactTOQsym(1))

      d0g(:)     = ZERO
      d1g(:,:)   = ZERO
      d2g(:,:,:) = ZERO

!C---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'd0g at Qact:',Qact
        write(out_unitp,*) d0g(:)
        write(out_unitp,*) 'END d0d1d2_g'
      END IF
!C---------------------------------------------------------------------

      END SUBROUTINE d0d1d2_g
!C================================================================
!C    analytical hessian along the path
!C================================================================
     SUBROUTINE d0d1d2_h(d0h,d1h,d2h,Qdyn,mole,deriv,num,step)
      USE mod_system
      USE mod_Tnum
      IMPLICIT NONE

!c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole


       real (kind=Rkind) :: Qdyn(mole%nb_var)
       real (kind=Rkind) :: Qact


       real (kind=Rkind) :: d0h(mole%nb_inact2n,mole%nb_inact2n)
       real (kind=Rkind) :: d1h(mole%nb_inact2n,mole%nb_inact2n)
       real (kind=Rkind) :: d2h(mole%nb_inact2n,mole%nb_inact2n)

       real (kind=Rkind) :: hess(mole%nb_act,mole%nb_act)

       logical           :: deriv,num
       real (kind=Rkind) :: step
       integer :: i,j,iQdyn,jQdyn

       real (kind=Rkind) :: k0

       !real (kind=Rkind), parameter :: b     = ZERO
       real (kind=Rkind), parameter :: b     = ONE
       real (kind=Rkind), parameter :: alpha = 3._Rkind
       !real (kind=Rkind), parameter :: v0    = 0.001822534_Rkind
       real (kind=Rkind), parameter :: v0    = 0.0018_Rkind
       real (kind=Rkind), parameter :: zmu   = 2000._Rkind

!c----- for debuging ----------------------------------
      logical, parameter :: debug = .FALSE.
!c     logical, parameter :: debug = .TRUE.
!c---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING d0d1d2_h'
        write(out_unitp,*) 'nb_var',mole%nb_var
        write(out_unitp,*) 'nb_act1',mole%nb_act1
        write(out_unitp,*) 'nb_inact22,nb_inact21',mole%nb_inact22,mole%nb_inact21
        write(out_unitp,*) 'nb_inact2n',mole%nb_inact2n
      END IF
!c---------------------------------------------------------------------

      k0 = zmu * V0**2


      Qact = Qdyn(mole%liste_QactTOQsym(1))

      IF (deriv) THEN
        write(out_unitp,*) 'ERROR in d0d1d2_h'
        write(out_unitp,*) '  deriv CANNOT be true!!'
        write(out_unitp,*) ' check the fortran source'
        STOP
      END IF

      d0h(:,:) = ZERO
      DO i=1,mole%nb_inact2n
       d0h(i,i) = k0 * ( ONE + b / COSH(alpha * Qact)**2 )
      END DO

!c---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Qact1',Qact
        CALL Write_RMat(d0h,6,4)
        write(out_unitp,*) 'END d0d1d2_h'
      END IF
!c---------------------------------------------------------------------

      END SUBROUTINE d0d1d2_h
!C================================================================
!C    analytical derivative (dnQflex : Qflex Qflex' Qflex" Qflex'") calculation
!C    for the variable iq
!C================================================================
      SUBROUTINE calc_dnQflex(iq,dnQflex,Qact,nb_act,nderiv,it)
      USE mod_system
      USE mod_dnSVM      
      IMPLICIT NONE

       integer           :: iq,nb_act
       real (kind=Rkind) :: Qact(nb_act)
       integer           :: nderiv,it
       TYPE (Type_dnS)   :: dnQflex




!C----- for debuging ----------------------------------
      character (len=*), parameter :: name_sub='dnQflex'
      logical, parameter :: debug=.FALSE.
!C     logical, parameter :: debug=.TRUE.
!C----- for debuging ----------------------------------


!C---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_act',nb_act
        write(out_unitp,*) 'iq',iq
      END IF
!C---------------------------------------------------------------------

!C---------------------------------------------------------------------
       CALL sub_ZERO_TO_dnS(dnQflex)

!C---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'dnQflex : ',Qact
        CALL write_dnS(dnQflex,nderiv)
        write(out_unitp,*) 'END ',name_sub
      END IF
!C---------------------------------------------------------------------

       END SUBROUTINE calc_dnQflex
!C================================================================
!C    analytical derivative (Qeq Qeq' Qeq" Qeq'") calculation
!C    for the variable i_Qdyn
!C================================================================
      SUBROUTINE d0d1d2d3_Qeq(i_Qdyn,d0req,d1req,d2req,d3req,Qdyn,mole,nderiv)
      USE mod_system
      USE mod_Tnum
      IMPLICIT NONE

!C----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

       integer :: i_Qdyn
       real (kind=Rkind) ::  Qdyn(mole%nb_var)

       integer :: nderiv

       real (kind=Rkind) :: d0req
       real (kind=Rkind) :: d1req(mole%nb_act1)
       real (kind=Rkind) :: d2req(mole%nb_act1,mole%nb_act1)
       real (kind=Rkind) :: d3req(mole%nb_act1,mole%nb_act1,mole%nb_act1)

      real (kind=Rkind) :: Qact,Qtot(mole%nb_act)


!c----- for debuging ----------------------------------
      character (len=*), parameter :: name_sub='d0d1d2d3_Qeq'
      logical, parameter :: debug=.FALSE.
!c     logical, parameter :: debug=.TRUE.
!c----- for debuging ----------------------------------

!c---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_inact20,nb_act',mole%nb_inact20,mole%nb_act
        write(out_unitp,*) 'nb_var',mole%nb_var
        write(out_unitp,*) 'i_Qdyn',i_Qdyn
      END IF
!c---------------------------------------------------------------------

!c---------------------------------------------------------------------
!c      Qact value. Rq: only ONE active variable is possible
!c---------------------------------------------------------------------
       IF (mole%nb_act1 /= 1) THEN
         write(out_unitp,*) ' ERROR : d0d1d2d3_Qeq'
         write(out_unitp,*) ' the number of Active variable'
         write(out_unitp,*) ' should be 1. But nb_act1 =',mole%nb_act1
         STOP
       END IF
!c---------------------------------------------------------------------
!c---------------------------------------------------------------------

!c---------------------------------------------------------------------

      d0req  = zero
      d1req  = zero
      d2req  = zero
      d3req  = zero
!c---------------------------------------------------------------------

!c---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'd0req : ',Qact,d0req
        IF (nderiv > 0) write(out_unitp,*) 'd1req : ',Qact,d1req
        IF (nderiv > 1) write(out_unitp,*) 'd2req : ',Qact,d2req
        IF (nderiv > 2) write(out_unitp,*) 'd3req : ',Qact,d3req
        write(out_unitp,*) 'END d0d1d2d3_Qeq'
      END IF
!c---------------------------------------------------------------------
      END SUBROUTINE d0d1d2d3_Qeq

SUBROUTINE pot_bottleneck_nd (r, n, Vpot, alpha, V0, b, zmu)

  ! - Time-stamp: <pot_bottleneck_nd.f90 - Apr 05, 2012 16:08:31 - G. Parlant>
  !
  ! - Language: F90/95
  !
  ! ////// description //////
  !
  !     (n+1)D bottleneck potential (Eckart + n harmonic oscillators):
  !
  !             V(x,y_i) = V0 * Sech(alpha * x)^2 + k(x)/2 * Sum_i |y_i|^2
  !   
  !             with    k(x) = k0 * [1 + b * Sech(alpha * x)^2]
  !   
  !                     k0 = m * Omega^2
  !                     hbar*Omega = V0
  !   
  !     Refs:
  !             Maddox & Poirier, Multidim. Quantum Mechanics with Trajs, CCP (2009))
  !
  !             Trahan, Wyatt and Poirier, J Chem Phys 122, 164104 (2005)
  !             Multidimensional quantum trajectories: Applications of the
  !             derivative propagation method.
  !
  ! ////////////////////////////////////////////////////////////////////////////

  USE mod_system

  IMPLICIT NONE 

  ! Subroutine scalars - intent(in)

  REAL(Rkind), INTENT(in) :: alpha, V0               !! Eckart constants
  REAL(Rkind), INTENT(in) :: b                       !! coupling parameter
  REAL(Rkind), INTENT(in) :: zmu                     !! mass
  INTEGER    , INTENT(in) :: n                       !! size of r (number of coordinates)

  ! Subroutine arrays - intent(in)

  REAL(Rkind), INTENT(in) :: r(n)                  !! internuclear distance

  ! Subroutine arrays - intent(out)

  REAL(Rkind), INTENT(out) :: Vpot                !! potential

  ! Local scalars

  INTEGER :: icoord                               !! do-loop index
  REAL(Rkind) :: k0

  ! Local arrays

  REAL(Rkind) :: k

  ! ============================   end of header ===============================

  k0 = zmu * V0**2

  k = k0 * ( one + b / COSH(alpha * r(1))**2 )

  ! Eckart barrier
  Vpot = V0 / COSH(alpha * r(1))**2

  ! add classical part
  Vpot = Vpot + half * k * sum(r(2:n)**2)
  
END SUBROUTINE pot_bottleneck_nd
