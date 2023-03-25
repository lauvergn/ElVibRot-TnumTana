!c
!C================================================================
!C    calc_Op : calculation of the potential and dipolar matrices
!c    mat_V(nb_be,nb_be) and mat_dip(nb_be,nb_be,3)
!c    nb_be : nb of elctronic surface
!c    Q are the coordinates in active order or syl order
!c    dipolar calculation if calc_dip = T
!C================================================================
      SUBROUTINE calcN_op(mat_V,mat_imV,mat_ScalOp,nb_be,nb_ScalOp, &
                          Q,nb_Q,mole,calc_ScalOp,pot_cplx)

      USE mod_system
      USE mod_Tnum
      IMPLICIT NONE

!c----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole

      integer           :: nb_be,nb_ScalOp,nb_Q
      logical           :: calc_ScalOp,pot_cplx
      real (kind=Rkind) :: mat_V(nb_be,nb_be),mat_imV(nb_be,nb_be)
      real (kind=Rkind) :: mat_ScalOp(nb_be,nb_be,nb_ScalOp)
      real (kind=Rkind) :: Q(nb_Q)

      real (kind=Rkind) :: im_pot0
      real (kind=Rkind) :: ScalOp(nb_ScalOp)

      
      logical,           save :: begin = .TRUE.
      real (kind=Rkind) :: Qdist(3),R1(3),R2(3),R3(3),GOO(3),R2Jac(3),R1Jac(3)
      integer :: ndim

      IF (nb_be == 1 ) THEN

!$OMP    CRITICAL (calcN_op_CRIT)
        IF (begin) THEN
          ndim=3
          !    Coordinates for H3 SLTH potential with the 3 distances
          CALL sub_Init_Qmodel(ndim,nb_be,'H3_LSTH',.FALSE.,0)
          write(6,*) 'end sub_Init_Qmodel in calcN_op'
          flush(6)
          begin = .FALSE.
        END IF
!$OMP   END CRITICAL (calcN_op_CRIT)

        ! we assume, Q(:) are the Cartesian coordinates
        R1(:) = Q(4:6)-Q(7:9)      ! H3->H2 also R1Jac
        R2(:) = Q(1:3)-Q(4:6)      ! H2->H1
        R3(:) = Q(1:3)-Q(7:9)      ! H3->H1
        Qdist(1) = sqrt(dot_product(R1,R1))
        Qdist(2) = sqrt(dot_product(R2,R2))
        Qdist(3) = sqrt(dot_product(R3,R3))
        CALL sub_Qmodel_V(mat_V,Qdist)

        !write(6,*) 'Q pot    ',Q(:),mat_V
         write(6,*) 'Qdist pot',Qdist(:),mat_V

        IF (pot_cplx) THEN
          mat_imV(:,:) = im_pot0(Qdist)
        END IF
        IF (calc_ScalOp) THEN
          !GOO(:) = (Q(4:6)+Q(7:9))*HALF
          !R1Jac(:) = R1(:)
          !R2Jac(:) = Q(1:3)-GOO(:)   ! GHH->H1
          !Qdist(1) = sqrt(dot_product(R1Jac,R1Jac))
          !Qdist(2) = sqrt(dot_product(R2Jac,R2Jac))
          CALL sub_ScalarOp(ScalOp,nb_ScalOp,Qdist,mole)
          mat_ScalOp(1,1,:) = ScalOp(:)
          !write(6,*) 'R1,R2',Qdist(1:2),'CAP',ScalOp(:)
        END IF
      ELSE
        write(6,*) ' ERROR in calc_op'
        write(6,*) ' It has to be used with ONE electronic surfaces'
        write(6,*) ' Rq: nb_be',nb_be
        STOP
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
      USE mod_CAP
      IMPLICIT NONE
      real(kind=Rkind) :: im_pot0
      real(kind=Rkind) Q(3)

      TYPE (CAP_t), parameter :: CAP1=CAP_t(type_CAP=1,n_exp=4,A=0.03_Rkind,   &
                                            B=FIVE/TWO,Q0=10.7,LQ=5._Rkind,    &
                                            ind_Q=1)
      TYPE (CAP_t), parameter :: CAP2=CAP_t(type_CAP=1,n_exp=4,A=0.03_Rkind,   &
                                            B=FIVE/TWO,Q0=10.7,LQ=5._Rkind,    &
                                            ind_Q=2)

       im_pot0 = calc_CAP(CAP1,Q) + calc_CAP(CAP2,Q)

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
! CAP for H+O2 (in Jacobi)
!
          !-----------------------------------
          !
          !
          !
          !                     H1
          !                     ^
          !                    /
          !                   / R2Jac=Q(2)=RR
          !                  /
          !                 /
          !                /
          !               / )th
          !      H2------/-------->H3 - - -> zBF
          !               R1Jac=Q(1)=r
          !     
          !
          !
          ! the coordinates are [R1,R2,cth] with cth=cos(th)
          !-----------------------------------
!C================================================================
      SUBROUTINE sub_ScalarOp(ScalOp,nb_ScalOp,Q,mole)
      USE mod_system
      USE mod_Tnum
      IMPLICIT NONE

!C----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole

       integer nb_Q,nb_ScalOp
       real(kind=Rkind) Q(3)
       real(kind=Rkind) ScalOp(nb_ScalOp)
       integer :: i

       real(kind=Rkind) :: r,RR
       ! Q(1): O2->O3:  r
       ! Q(2): GOO->H1: RR
       
       real(kind=Rkind), parameter :: RRcut = 5._Rkind
       real(kind=Rkind), parameter :: LRR   = 1._Rkind
       real(kind=Rkind), parameter :: ARR   = 1._Rkind
       real(kind=Rkind), parameter :: rcut  = 5._Rkind
       real(kind=Rkind), parameter :: Lr    = 1._Rkind
       real(kind=Rkind), parameter :: Ar    = 1._Rkind

       IF (nb_ScalOp /= 2) STOP 'WRONG nb_ScalOp value. It MUST be 2'

       r  = Q(1) ! O-O
       !RR = Q(2)/sqrt(THREE)/TWO ! H-O2
       RR = Q(2) ! H-H2

       ! defined the reactant CAP (along RR: H+H2)
       IF (RR > RRcut) THEN
         ScalOp(1) = 13.22_Rkind*ARR*exp(-TWO*LRR/(RR-RRcut))
       ELSE
         ScalOp(1) = ZERO
       END IF

       ! defined the product CAP (along r: O+O)
       IF (r > rcut) THEN
         ScalOp(2) = 13.22_Rkind*Ar*exp(-TWO*Lr/(r-rcut))
       ELSE
         ScalOp(2) = ZERO
       END IF


       END
!C================================================================
!C    analytical gradient along the path
!C================================================================
      SUBROUTINE d0d1d2_g(d0g,d1g,d2g,Qdyn,mole,deriv,num,step)
      USE mod_system
      USE mod_Tnum
      IMPLICIT NONE

!C----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole

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

!c----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole


       real (kind=Rkind) :: Qdyn(mole%nb_var)
       real (kind=Rkind) :: Qact


       real (kind=Rkind) :: d0h(mole%nb_inact2n,mole%nb_inact2n)
       real (kind=Rkind) :: d1h(mole%nb_inact2n,mole%nb_inact2n)
       real (kind=Rkind) :: d2h(mole%nb_inact2n,mole%nb_inact2n)

       real (kind=Rkind) :: hess(mole%nb_act,mole%nb_act)

       logical           :: deriv,num
       real (kind=Rkind) :: step
       integer :: i,j,iQdyn,jQdyn


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

      Qact = Qdyn(mole%liste_QactTOQsym(1))

      IF (deriv) THEN
        write(out_unitp,*) 'ERROR in d0d1d2_h'
        write(out_unitp,*) '  deriv CANNOT be true!!'
        write(out_unitp,*) ' check the fortran source'
        STOP
      END IF

      CALL get_CurviRPH( (/ Qact /),mole%CurviRPH,Hess=hess)

      DO i=1,mole%nb_inact2n
      DO j=1,mole%nb_inact2n
       iQdyn = mole%liste_QactTOQsym(i+1)
       jQdyn = mole%liste_QactTOQsym(j+1)
       d0h(j,i) = hess(jQdyn,iQdyn)
      END DO
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



       real (kind=Rkind), parameter :: x0  = 1.0_Rkind
       real (kind=Rkind), parameter :: RH2 = 1.401036_Rkind
       real (kind=Rkind), parameter :: RTS = 1.7570_Rkind
       real (kind=Rkind), parameter :: RR0 = RTS-RH2-x0
       real (kind=Rkind), parameter :: b = TWO
       real (kind=Rkind) :: Rm
       real (kind=Rkind) :: d0,d1,d2,d3
       real (kind=Rkind) :: r0,r1,r2,r3

!C----- for debuging ----------------------------------
      character (len=*), parameter :: name_sub='dnQflex'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!C----- for debuging ----------------------------------


!C---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_act',nb_act
        write(out_unitp,*) 'iq',iq
        write(out_unitp,*) 'Qact',Qact(:)
      END IF
!C---------------------------------------------------------------------

!C---------------------------------------------------------------------
       CALL sub_ZERO_TO_dnS(dnQflex)
       Rm = Qact(1) ! R-

       IF (iQ == 1) THEN
         d0 = sqrt(x0**2 + Rm**2)
         r0 = RR0*exp(-b*Rm**2)
         dnQflex%d0 = RH2 + r0 + d0
         IF (nderiv > 0) THEN 
           d1     = Rm/d0
           r1     = -TWO*b*Rm * r0
           dnQflex%d1(1)   = d1 + r1
         END IF
         IF (nderiv > 1) THEN 
           d2     = (ONE-d1**2)/d0
           r2     = TWO*b*(-ONE+TWO*b*Rm**2) * r0
           dnQflex%d2(1,1)   = d2 + r2
         END IF
         IF (nderiv > 2) THEN 
           d3     = -THREE*d2*d1/d0
           r3     = FOUR*b**2*Rm*(THREE-TWO*b*Rm**2) * r0
           dnQflex%d3(1,1,1) = d3 + r3
         END IF

       END IF

!C---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'dnQflex : ',Qact,dnQflex%d0
        !ALL write_dnS(dnQflex,nderiv)
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

!C----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole

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
