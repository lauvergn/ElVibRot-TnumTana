C**********************************************************************
C
C     Copyright (C) 1992  Roland W. Freund and Noel M. Nachtigal
C     All rights reserved.
C
C     This code is part of a copyrighted package.  For details, see the
C     file `cpyrit.doc' in the top-level directory.
C
C     *****************************************************************
C     ANY USE OF  THIS CODE CONSTITUTES ACCEPTANCE OF  THE TERMS OF THE
C     COPYRIGHT NOTICE
C     *****************************************************************
C
C**********************************************************************
C
C     This  file  contains  the  routines  for  the  QMR algorithm  for
C     symmetric matrices, using the coupled two-term recurrence variant
C     of the Lanczos algorithm without look-ahead.
C
C**********************************************************************
C
      SUBROUTINE ZSCPX (NDIM,NLEN,NLIM,VECS,TOL,INFO)
      use mod_system, only : Rkind

C     Purpose:
C     This subroutine uses the QMR algorithm based  on the coupled two-
C     term variant of the Lanczos  process without look-ahead  to solve
C     linear systems.  It runs the algorithm to convergence  or until a
C     user-specified limit on the number of  iterations is reached.  It
C     is  set up  to solve  symmetric  systems  starting with identical
C     starting vectors.
C
C     The  code is  set up  to solve  the system  A x = b  with initial
C     guess x_0 = 0.  Here  A x = b  denotes the preconditioned system,
C     and  it  is  connected with the  original system as follows.  Let
C     B y = c be the original unpreconditioned system to be solved, and
C     let y_0 be an arbitrary initial guess for its solution.  Then:
C     A x = b, where  A = M_1^{-1} B M_2^{-1},
C     x = M_2 (y - y_0), b = M_1^{-1} (c - B y_0).
C     Here M = M_1 M_2 is the preconditioner.
C
C     To recover the final iterate y_n for  the original system B y = c
C     from the final iterate x_n for the preconditioned system A x = b,
C     set
C     y_n = y_0 + M_2^{-1} x_n.
C
C     The algorithm  was first described  in the RIACS Technical Report
C     92.15, "An Implementation of the QMR Method Based on Coupled Two-
C     Term Recurrences",  June 1992.  This implementation does not have
C     look-ahead, so it is less robust than the full version.
C
C     Parameters:
C     For a description of  the parameters, see the file `zscpx.doc' in
C     the current directory.
C
C     External routines used:
C     double precision dlamch(ch)
C     LAPACK routine, computes machine-related constants.
C     double precision dznrm2(n,x,incx)
C     BLAS-1 routine, computes the 2-norm of x.
C     subroutine zaxpby(n,z,a,x,b,y)
C     Library routine, computes z = a * x + b * y.
C     double precision zdotu(n,x,incx,y,incy)
C     BLAS-1 routine, computes y' * x.
C     subroutine zrandn(n,x,seed)
C     Library routine, fills x with random numbers.
C     subroutine zrotg(a,b,cos,sin)
C     BLAS-1 routine, computes the Givens rotation which rotates the
C     vector [a; b] into [ sqrt(a**2 + b**2); 0 ].
C     double precision zscpxo(n)
C     User-supplied routine, specifies the QMR scaling factors.
C
C     Noel M. Nachtigal
C     March 30, 1993
C
C**********************************************************************
C
      !INTRINSIC CDABS, DABS, DBLE, DCMPLX, DCONJG, DMAX1, DSQRT, MAX0
      EXTERNAL DLAMCH, DZNRM2, ZAXPBY, ZDOTU, ZRANDN, ZROTG, ZSCPXO
      DOUBLE COMPLEX ZDOTU
      DOUBLE PRECISION DLAMCH, DZNRM2, ZSCPXO
C
      INTEGER INFO(4), NDIM, NLEN, NLIM
      DOUBLE COMPLEX VECS(NDIM,6)
      DOUBLE PRECISION TOL
C
C     Miscellaneous parameters.
C
      DOUBLE COMPLEX ZONE, ZZERO
      PARAMETER (ZONE = (1.0D0,0.0D0),ZZERO = (0.0D0,0.0D0))
      DOUBLE PRECISION DHUN, DONE, DTEN, DZERO
      PARAMETER (DHUN = 1.0D2,DONE = 1.0D0,DTEN = 1.0D1,DZERO = 0.0D0)
C
C     Local variables, permanent.
C
      INTEGER IERR, N, RETLBL, TF, TRES, VF
      SAVE    IERR, N, RETLBL, TF, TRES, VF
      DOUBLE COMPLEX DNN, ENN, SCS, SINN, RHSN
      SAVE           DNN, ENN, SCS, SINN, RHSN
      DOUBLE PRECISION COSN, GAMN, LNP1N, MAXOMG, OMG, R0, SCPN, SCQN
      SAVE             COSN, GAMN, LNP1N, MAXOMG, OMG, R0, SCPN, SCQN
      DOUBLE PRECISION SCV, RESN, TMAX, TMIN, TNRM, UCHK, UNRM
      SAVE             SCV, RESN, TMAX, TMIN, TNRM, UCHK, UNRM
C
C     Local variables, transient.
C
      INTEGER INIT, REVCOM
      DOUBLE COMPLEX LNN, RHN, RHNM1, RHNP1, RHSNP1, UNM1N, ZTMP
      DOUBLE PRECISION GAMNM1, SCW
C
C     Initialize some of the permanent variables.
C
      DATA RETLBL /0/
C
C     Check the reverse communication flag to see where to branch.
C     REVCOM   RETLBL      Comment
C     0        0    first call, go to label 10
C     1       30    returning from AXB, go to label 30
C     1       50    returning from AXB, go to label 50
C     2       40    returning from ATXB, go to label 40
C
      REVCOM  = INFO(2)
      INFO(2) = 0
      IF (REVCOM.EQ.0) THEN
         N = 0
         IF (RETLBL.EQ.0) GO TO 10
      ELSE IF (REVCOM.EQ.1) THEN
         IF (RETLBL.EQ.30) THEN
            GO TO 30
         ELSE IF (RETLBL.EQ.50) THEN
            GO TO 50
         END IF
      ELSE IF (REVCOM.EQ.2) THEN
         IF (RETLBL.EQ.40) GO TO 40
      END IF
      IERR = 1
      GO TO 70
C
C     Check whether the inputs are valid.
C
 10   IERR = 0
      IF (NDIM.LT.1)        IERR = 2
      IF (NLEN.LT.1)        IERR = 2
      IF (NLIM.LT.1)        IERR = 2
      IF (NLEN.GT.NDIM)     IERR = 2
      IF (IERR.NE.0) GO TO 70
C
C     Extract from INFO the output units TF and VF, the true residual
C     flag TRES, and the left starting vector flag INIT.
C
      VF   = MAX0(INFO(1),0)
      INIT = VF / 100000
      VF   = VF - INIT * 100000
      TRES = VF / 10000
      VF   = VF - TRES * 10000
      TF   = VF / 100
      VF   = VF - TF * 100
C
C     Extract and check the various tolerances.
C
      TNRM = DLAMCH('E') * DTEN
      TMIN = DSQRT(DSQRT(DLAMCH('S')))
      TMAX = DONE / TMIN
      IF (TOL.LE.DZERO) TOL = DSQRT(DLAMCH('E'))
C
C     Start the trace messages and convergence history.
C
      IF (VF.NE.0) WRITE (VF,'(I8,2E11.4)') 0, DONE, DONE
      IF (TF.NE.0) WRITE (TF,'(I8,2E11.4)') 0, DONE, DONE
C
C     Set x_0 = 0 and compute the norm of the initial residual.
C
      CALL ZAXPBY (NLEN,VECS(1,3),ZONE,VECS(1,2),ZZERO,VECS(1,3))
      CALL ZAXPBY (NLEN,VECS(1,1),ZZERO,VECS(1,1),ZZERO,VECS(1,1))
      R0 = DZNRM2(NLEN,VECS(1,3),1)
      IF ((TOL.GE.DONE).OR.(R0.EQ.DZERO)) GO TO 70
C
C     Check whether the auxiliary vector must be supplied.
C
C     SYM  IF (INIT.EQ.0) CALL ZRANDN (NLEN,VECS(1,3),1)
C
C     Initialize the variables.
C
      N      = 1
      SCV    = R0
      ENN    = ZONE
      COSN   = DONE
      GAMN   = DONE
      RESN   = DONE
      SCPN   = DONE
      SCQN   = DONE
      SCS    = ZZERO
      SINN   = ZZERO
      LNP1N  = DZERO
      OMG    = ZSCPXO(N)
      RHSN   = OMG * R0
      MAXOMG = DONE / OMG
C
C     This is one step of the coupled two-term Lanczos algorithm.
C     Check whether E_n is nonsingular.
C
 20   IF (abs(ENN).EQ.DZERO) THEN
         IERR = 8
         GO TO 70
      END IF
C
C     Compute scale factor for the vector w_{n}.
C     Check for invariant subspaces, and scale the vectors if needed.
C
      IERR = 0
C     SYM  SCW  = DZNRM2(NLEN,VECS(1,3),1)
      SCW  = SCV
      IF (SCPN*SCV.LT.TNRM) IERR = IERR + 16
      IF (SCQN*SCW.LT.TNRM) IERR = IERR + 32
      IF (IERR.NE.0) GO TO 70
      GAMNM1 = GAMN
      GAMN   = GAMN * SCPN / SCQN * SCV / SCW
      DNN    = ZDOTU(NLEN,VECS(1,3),1,VECS(1,3),1) / ( SCV * SCW )
      IF ((SCV.GE.TMAX).OR.(SCV.LE.TMIN)) THEN
         ZTMP = cmplx(DONE / SCV,DZERO,kind=Rkind)
         CALL ZAXPBY (NLEN,VECS(1,3),ZTMP,VECS(1,3),ZZERO,VECS(1,3))
         SCV = DONE
      END IF
      IF ((SCW.GE.TMAX).OR.(SCW.LE.TMIN)) THEN
         ZTMP = cmplx(DONE / SCW,DZERO,kind=Rkind)
C     SYM     CALL ZAXPBY (NLEN,VECS(1,3),ZTMP,VECS(1,3),ZZERO,VECS(1,3))
         SCW = DONE
      END IF
      SCV = DONE / SCV
      SCW = DONE / SCW
C
C     Build the vectors p_n and q_n.
C
      UNM1N = DNN * LNP1N * GAMNM1 / ( GAMN * ENN )
      ZTMP  = UNM1N * SCPN / SCV
      CALL ZAXPBY (NLEN,VECS(1,4),ZONE,VECS(1,3),-ZTMP,VECS(1,4))
      ZTMP  = UNM1N * SCQN / SCW * GAMN / GAMNM1
C     SYM  CALL ZAXPBY (NLEN,VECS(1,4),ZONE,VECS(1,3),-ZTMP,VECS(1,4))
      SCPN  = SCV
      SCQN  = SCW
C
C     Check whether D_n is nonsingular.
C
      IF (abs(DNN).EQ.DZERO) THEN
         IERR = 8
         GO TO 70
      END IF
C
C     Have the caller carry out AXB, then return here.
C     CALL AXB (VECS(1,4),VECS(1,6))
C
      INFO(2) = 1
      INFO(3) = 4
      INFO(4) = 6
      RETLBL  = 30
      RETURN
C
C     Compute q_n^T A p_n.
C
 30   ENN = SCPN * SCQN * ZDOTU(NLEN,VECS(1,4),1,VECS(1,6),1)
C
C     Build the vector v_{n+1}.
C
      LNN = ENN / DNN
      CALL ZAXPBY (NLEN,VECS(1,3),ZONE,VECS(1,6),-LNN,VECS(1,3))
C
C     Have the caller carry out ATXB, then return here.
C     CALL ATXB (VECS(1,4),VECS(1,6))
C
      INFO(2) = 2
      INFO(3) = 8
      INFO(4) = 6
      RETLBL  = 40
C     SYM  RETURN
C
C     Build the vector w_{n+1}.
C
 40   LNN = ENN / DNN
C     SYM  CALL ZAXPBY (NLEN,VECS(1,3),ZONE,VECS(1,6),-LNN,VECS(1,3))
C
C     Compute scale factor for the vector v_{n+1}.
C
      SCV    = DZNRM2(NLEN,VECS(1,3),1)
      LNP1N  = SCPN * SCV
C
C     The QMR code starts here.
C     Multiply the new column by the previous omega's.
C     Get the next scaling factor omega(i) and update MAXOMG.
C
      RHN    = OMG * LNN
      OMG    = ZSCPXO(N+1)
      RHNP1  = OMG * LNP1N
      MAXOMG = DMAX1(MAXOMG,DONE/OMG)
C
C     Apply the previous rotation.
C
      RHNM1 = SINN * RHN
      RHN   = COSN * RHN
C
C     Compute the rotation for the last element (this also applies it).
C
      CALL ZROTG (RHN,RHNP1,COSN,SINN)
C
C     Apply the new rotation to the right-hand side vector.
C
      RHSNP1 = -conjg(SINN) * RHSN
      RHSN   =  COSN * RHSN
C
C     Compute the next search direction s_i.
C
      ZTMP = RHNM1 * SCS / SCPN
      CALL ZAXPBY (NLEN,VECS(1,5),ZONE,VECS(1,4),-ZTMP,VECS(1,5))
C
C     Compute the new QMR iterate, then scale the search direction.
C
      SCS  = SCPN / RHN
      ZTMP = SCS * RHSN
      CALL ZAXPBY (NLEN,VECS(1,1),ZONE,VECS(1,1),ZTMP,VECS(1,5))
      IF ((abs(SCS).GE.TMAX).OR.(abs(SCS).LE.TMIN)) THEN
         CALL ZAXPBY (NLEN,VECS(1,5),SCS,VECS(1,5),ZZERO,VECS(1,5))
         SCS = ZONE
      END IF
C
C     Compute the residual norm upper bound.
C     If the scaled upper bound is within one order of magnitude of the
C     target convergence norm, compute the true residual norm.
C
      RHSN = RHSNP1
      UNRM = DSQRT(DBLE(N+1)) * MAXOMG * abs(RHSNP1) / R0
      UCHK = UNRM
      IF ((TRES.EQ.0).AND.(UNRM/TOL.GT.DTEN).AND.(N.LT.NLIM)) GO TO 60
C
C     Have the caller carry out AXB, then return here.
C     CALL AXB (VECS(1,1),VECS(1,6))
C
      INFO(2) = 1
      INFO(3) = 1
      INFO(4) = 6
      RETLBL  = 50
      RETURN
 50   CALL ZAXPBY (NLEN,VECS(1,6),ZONE,VECS(1,2),-ZONE,VECS(1,6))
      RESN = DZNRM2(NLEN,VECS(1,6),1) / R0
      UCHK = RESN
C
C     Output the convergence history.
C
 60   IF (VF.NE.0) WRITE (VF,'(I8,2E11.4)') N, UNRM, RESN
      IF (TF.NE.0) WRITE (TF,'(I8,2E11.4)') N, UNRM, RESN
C
C     Check for convergence or termination.  Stop if:
C     1. algorithm converged;
C     2. there is an error condition;
C     3. the residual norm upper bound is smaller than the computed
C     residual norm by a factor of at least 100;
C     4. algorithm exceeded the iterations limit.
C
      IF (RESN.LE.TOL) THEN
         IERR = 0
         GO TO 70
      ELSE IF (IERR.NE.0) THEN
         GO TO 70
      ELSE IF (UNRM.LT.UCHK/DHUN) THEN
         IERR = 4
         GO TO 70
      ELSE IF (N.GE.NLIM) THEN
         IERR = 4
         GO TO 70
      END IF
C
C     Update the running counter.
C
      N = N + 1
      GO TO 20
C
C     That's all.
C
 70   NLIM    = N
      RETLBL  = 0
      INFO(1) = IERR
C
      RETURN
      END
C
C**********************************************************************
C
      DOUBLE PRECISION FUNCTION ZSCPXO (I)
C
C     Purpose:
C     Returns the scaling parameter OMEGA(I).
C
C     Parameters:
C     I = the index of the parameter OMEGA (input).
C
C     Noel M. Nachtigal
C     March 30, 1993
C
C**********************************************************************
C
      INTEGER I
C
      ZSCPXO = 1.0D0
C
      RETURN
      END
C
C**********************************************************************
      DOUBLE PRECISION FUNCTION DLAMCH( CMACH )
*
*     -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          CMACH
*     ..
*
*     Purpose
*     =======
*
*     DLAMCH determines double precision machine parameters.
*
*     Arguments
*     =========
*
*     CMACH   (input) CHARACTER*1
*     Specifies the value to be returned by DLAMCH:
*     = 'E' or 'e',   DLAMCH := eps
*     = 'S' or 's ,   DLAMCH := sfmin
*     = 'B' or 'b',   DLAMCH := base
*     = 'P' or 'p',   DLAMCH := eps*base
*     = 'N' or 'n',   DLAMCH := t
*     = 'R' or 'r',   DLAMCH := rnd
*     = 'M' or 'm',   DLAMCH := emin
*     = 'U' or 'u',   DLAMCH := rmin
*     = 'L' or 'l',   DLAMCH := emax
*     = 'O' or 'o',   DLAMCH := rmax
*
*     where
*
*     eps   = relative machine precision
*     sfmin = safe minimum, such that 1/sfmin does not overflow
*     base  = base of the machine
*     prec  = eps*base
*     t     = number of (base) digits in the mantissa
*     rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*     emin  = minimum exponent before (gradual) underflow
*     rmin  = underflow threshold - base**(emin-1)
*     emax  = largest exponent before overflow
*     rmax  = overflow threshold  - (base**emax)*(1-eps)
*
*     =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRST, LRND
      INTEGER            BETA, IMAX, IMIN, IT
      DOUBLE PRECISION   BASE, EMAX, EMIN, EPS, PREC, RMACH, RMAX, RMIN,
     $   RND, SFMIN, SMALL, T
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAMC2
*     ..
*     .. Save statement ..
      SAVE               FIRST, EPS, SFMIN, BASE, T, RND, EMIN, RMIN,
     $   EMAX, RMAX, PREC
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         CALL DLAMC2( BETA, IT, LRND, EPS, IMIN, RMIN, IMAX, RMAX )
         BASE = BETA
         T = IT
         IF( LRND ) THEN
            RND = ONE
            EPS = ( BASE**( 1-IT ) ) / 2
         ELSE
            RND = ZERO
            EPS = BASE**( 1-IT )
         END IF
         PREC = EPS*BASE
         EMIN = IMIN
         EMAX = IMAX
         SFMIN = RMIN
         SMALL = ONE / RMAX
         IF( SMALL.GE.SFMIN ) THEN
*
*     Use SMALL plus a bit, to avoid the possibility of rounding
*     causing overflow when computing  1/sfmin.
*
            SFMIN = SMALL*( ONE+EPS )
         END IF
      END IF
*
      IF( LSAME( CMACH, 'E' ) ) THEN
         RMACH = EPS
      ELSE IF( LSAME( CMACH, 'S' ) ) THEN
         RMACH = SFMIN
      ELSE IF( LSAME( CMACH, 'B' ) ) THEN
         RMACH = BASE
      ELSE IF( LSAME( CMACH, 'P' ) ) THEN
         RMACH = PREC
      ELSE IF( LSAME( CMACH, 'N' ) ) THEN
         RMACH = T
      ELSE IF( LSAME( CMACH, 'R' ) ) THEN
         RMACH = RND
      ELSE IF( LSAME( CMACH, 'M' ) ) THEN
         RMACH = EMIN
      ELSE IF( LSAME( CMACH, 'U' ) ) THEN
         RMACH = RMIN
      ELSE IF( LSAME( CMACH, 'L' ) ) THEN
         RMACH = EMAX
      ELSE IF( LSAME( CMACH, 'O' ) ) THEN
         RMACH = RMAX
      END IF
*
      DLAMCH = RMACH
      RETURN
*
*     End of DLAMCH
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC1( BETA, T, RND, IEEE1 )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE1, RND
      INTEGER            BETA, T
*     ..
*
*  Purpose
*  =======
*
*  DLAMC1 determines the machine parameters given by BETA, T, RND, and
*  IEEE1.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  IEEE1   (output) LOGICAL
*          Specifies whether rounding appears to be done in the IEEE
*          'round to nearest' style.
*
*  Further Details
*  ===============
*
*  The routine is based on the routine  ENVRON  by Malcolm and
*  incorporates suggestions by Gentleman and Marovich. See
*
*     Malcolm M. A. (1972) Algorithms to reveal properties of
*        floating-point arithmetic. Comms. of the ACM, 15, 949-951.
*
*     Gentleman W. M. and Marovich S. B. (1974) More on algorithms
*        that reveal properties of floating point arithmetic units.
*        Comms. of the ACM, 17, 276-277.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, LIEEE1, LRND
      INTEGER            LBETA, LT
      DOUBLE PRECISION   A, B, C, F, ONE, QTR, SAVEC, T1, T2
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Save statement ..
      SAVE               FIRST, LIEEE1, LBETA, LRND, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ONE = 1
*
*        LBETA,  LIEEE1,  LT and  LRND  are the  local values  of  BETA,
*        IEEE1, T and RND.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are  stored and not held in registers,  or
*        are not affected by optimizers.
*
*        Compute  a = 2.0**m  with the  smallest positive integer m such
*        that
*
*           fl( a + 1.0 ) = a.
*
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   10    CONTINUE
         IF( C.EQ.ONE ) THEN
            A = 2*A
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 10
         END IF
*+       END WHILE
*
*        Now compute  b = 2.0**m  with the smallest positive integer m
*        such that
*
*           fl( a + b ) .gt. a.
*
         B = 1
         C = DLAMC3( A, B )
*
*+       WHILE( C.EQ.A )LOOP
   20    CONTINUE
         IF( C.EQ.A ) THEN
            B = 2*B
            C = DLAMC3( A, B )
            GO TO 20
         END IF
*+       END WHILE
*
*        Now compute the base.  a and c  are neighbouring floating point
*        numbers  in the  interval  ( beta**t, beta**( t + 1 ) )  and so
*        their difference is beta. Adding 0.25 to c is to ensure that it
*        is truncated to beta and not ( beta - 1 ).
*
         QTR = ONE / 4
         SAVEC = C
         C = DLAMC3( C, -A )
         LBETA = C + QTR
*
*        Now determine whether rounding or chopping occurs,  by adding a
*        bit  less  than  beta/2  and a  bit  more  than  beta/2  to  a.
*
         B = LBETA
         F = DLAMC3( B / 2, -B / 100 )
         C = DLAMC3( F, A )
         IF( C.EQ.A ) THEN
            LRND = .TRUE.
         ELSE
            LRND = .FALSE.
         END IF
         F = DLAMC3( B / 2, B / 100 )
         C = DLAMC3( F, A )
         IF( ( LRND ) .AND. ( C.EQ.A ) )
     $      LRND = .FALSE.
*
*        Try and decide whether rounding is done in the  IEEE  'round to
*        nearest' style. B/2 is half a unit in the last place of the two
*        numbers A and SAVEC. Furthermore, A is even, i.e. has last  bit
*        zero, and SAVEC is odd. Thus adding B/2 to A should not  change
*        A, but adding B/2 to SAVEC should change SAVEC.
*
         T1 = DLAMC3( B / 2, A )
         T2 = DLAMC3( B / 2, SAVEC )
         LIEEE1 = ( T1.EQ.A ) .AND. ( T2.GT.SAVEC ) .AND. LRND
*
*        Now find  the  mantissa, t.  It should  be the  integer part of
*        log to the base beta of a,  however it is safer to determine  t
*        by powering.  So we find t as the smallest positive integer for
*        which
*
*           fl( beta**t + 1.0 ) = 1.0.
*
         LT = 0
         A = 1
         C = 1
*
*+       WHILE( C.EQ.ONE )LOOP
   30    CONTINUE
         IF( C.EQ.ONE ) THEN
            LT = LT + 1
            A = A*LBETA
            C = DLAMC3( A, ONE )
            C = DLAMC3( C, -A )
            GO TO 30
         END IF
*+       END WHILE
*
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      IEEE1 = LIEEE1
      RETURN
*
*     End of DLAMC1
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC2( BETA, T, RND, EPS, EMIN, RMIN, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            RND
      INTEGER            BETA, EMAX, EMIN, T
      DOUBLE PRECISION   EPS, RMAX, RMIN
*     ..
*
*  Purpose
*  =======
*
*  DLAMC2 determines the machine parameters specified in its argument
*  list.
*
*  Arguments
*  =========
*
*  BETA    (output) INTEGER
*          The base of the machine.
*
*  T       (output) INTEGER
*          The number of ( BETA ) digits in the mantissa.
*
*  RND     (output) LOGICAL
*          Specifies whether proper rounding  ( RND = .TRUE. )  or
*          chopping  ( RND = .FALSE. )  occurs in addition. This may not
*          be a reliable guide to the way in which the machine performs
*          its arithmetic.
*
*  EPS     (output) DOUBLE PRECISION
*          The smallest positive number such that
*
*             fl( 1.0 - EPS ) .LT. 1.0,
*
*          where fl denotes the computed value.
*
*  EMIN    (output) INTEGER
*          The minimum exponent before (gradual) underflow occurs.
*
*  RMIN    (output) DOUBLE PRECISION
*          The smallest normalized number for the machine, given by
*          BASE**( EMIN - 1 ), where  BASE  is the floating point value
*          of BETA.
*
*  EMAX    (output) INTEGER
*          The maximum exponent before overflow occurs.
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest positive number for the machine, given by
*          BASE**EMAX * ( 1 - EPS ), where  BASE  is the floating point
*          value of BETA.
*
*  Further Details
*  ===============
*
*  The computation of  EPS  is based on a routine PARANOIA by
*  W. Kahan of the University of California at Berkeley.
*
* =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            FIRST, IEEE, IWARN, LIEEE1, LRND
      INTEGER            GNMIN, GPMIN, I, LBETA, LEMAX, LEMIN, LT,
     $                   NGNMIN, NGPMIN
      DOUBLE PRECISION   A, B, C, HALF, LEPS, LRMAX, LRMIN, ONE, RBASE,
     $                   SIXTH, SMALL, THIRD, TWO, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAMC1, DLAMC4, DLAMC5
*     ..
*     .. Intrinsic Functions ..
      !INTRINSIC          ABS, MAX, MIN
*     ..
*     .. Save statement ..
      SAVE               FIRST, IWARN, LBETA, LEMAX, LEMIN, LEPS, LRMAX,
     $                   LRMIN, LT
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. / , IWARN / .FALSE. /
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         ZERO = 0
         ONE = 1
         TWO = 2
*
*        LBETA, LT, LRND, LEPS, LEMIN and LRMIN  are the local values of
*        BETA, T, RND, EPS, EMIN and RMIN.
*
*        Throughout this routine  we use the function  DLAMC3  to ensure
*        that relevant values are stored  and not held in registers,  or
*        are not affected by optimizers.
*
*        DLAMC1 returns the parameters  LBETA, LT, LRND and LIEEE1.
*
         CALL DLAMC1( LBETA, LT, LRND, LIEEE1 )
*
*        Start to find EPS.
*
         B = LBETA
         A = B**( -LT )
         LEPS = A
*
*        Try some tricks to see whether or not this is the correct  EPS.
*
         B = TWO / 3
         HALF = ONE / 2
         SIXTH = DLAMC3( B, -HALF )
         THIRD = DLAMC3( SIXTH, SIXTH )
         B = DLAMC3( THIRD, -HALF )
         B = DLAMC3( B, SIXTH )
         B = ABS( B )
         IF( B.LT.LEPS )
     $      B = LEPS
*
         LEPS = 1
*
*+       WHILE( ( LEPS.GT.B ).AND.( B.GT.ZERO ) )LOOP
   10    CONTINUE
         IF( ( LEPS.GT.B ) .AND. ( B.GT.ZERO ) ) THEN
            LEPS = B
            C = DLAMC3( HALF*LEPS, ( TWO**5 )*( LEPS**2 ) )
            C = DLAMC3( HALF, -C )
            B = DLAMC3( HALF, C )
            C = DLAMC3( HALF, -B )
            B = DLAMC3( HALF, C )
            GO TO 10
         END IF
*+       END WHILE
*
         IF( A.LT.LEPS )
     $      LEPS = A
*
*        Computation of EPS complete.
*
*        Now find  EMIN.  Let A = + or - 1, and + or - (1 + BASE**(-3)).
*        Keep dividing  A by BETA until (gradual) underflow occurs. This
*        is detected when we cannot recover the previous A.
*
         RBASE = ONE / LBETA
         SMALL = ONE
         DO 20 I = 1, 3
            SMALL = DLAMC3( SMALL*RBASE, ZERO )
   20    CONTINUE
         A = DLAMC3( ONE, SMALL )
         CALL DLAMC4( NGPMIN, ONE, LBETA )
         CALL DLAMC4( NGNMIN, -ONE, LBETA )
         CALL DLAMC4( GPMIN, A, LBETA )
         CALL DLAMC4( GNMIN, -A, LBETA )
         IEEE = .FALSE.
*
         IF( ( NGPMIN.EQ.NGNMIN ) .AND. ( GPMIN.EQ.GNMIN ) ) THEN
            IF( NGPMIN.EQ.GPMIN ) THEN
               LEMIN = NGPMIN
*            ( Non twos-complement machines, no gradual underflow;
*              e.g.,  VAX )
            ELSE IF( ( GPMIN-NGPMIN ).EQ.3 ) THEN
               LEMIN = NGPMIN - 1 + LT
               IEEE = .TRUE.
*            ( Non twos-complement machines, with gradual underflow;
*              e.g., IEEE standard followers )
            ELSE
               LEMIN = MIN( NGPMIN, GPMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( NGPMIN.EQ.GPMIN ) .AND. ( NGNMIN.EQ.GNMIN ) ) THEN
            IF( ABS( NGPMIN-NGNMIN ).EQ.1 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN )
*            ( Twos-complement machines, no gradual underflow;
*              e.g., CYBER 205 )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE IF( ( ABS( NGPMIN-NGNMIN ).EQ.1 ) .AND.
     $            ( GPMIN.EQ.GNMIN ) ) THEN
            IF( ( GPMIN-MIN( NGPMIN, NGNMIN ) ).EQ.3 ) THEN
               LEMIN = MAX( NGPMIN, NGNMIN ) - 1 + LT
*            ( Twos-complement machines with gradual underflow;
*              no known machine )
            ELSE
               LEMIN = MIN( NGPMIN, NGNMIN )
*            ( A guess; no known machine )
               IWARN = .TRUE.
            END IF
*
         ELSE
            LEMIN = MIN( NGPMIN, NGNMIN, GPMIN, GNMIN )
*         ( A guess; no known machine )
            IWARN = .TRUE.
         END IF
***
* Comment out this if block if EMIN is ok
         IF( IWARN ) THEN
            FIRST = .TRUE.
            WRITE( 6, FMT = 9999 )LEMIN
         END IF
***
*
*        Assume IEEE arithmetic if we found denormalised  numbers above,
*        or if arithmetic seems to round in the  IEEE style,  determined
*        in routine DLAMC1. A true IEEE machine should have both  things
*        true; however, faulty machines may have one or the other.
*
         IEEE = IEEE .OR. LIEEE1
*
*        Compute  RMIN by successive division by  BETA. We could compute
*        RMIN as BASE**( EMIN - 1 ),  but some machines underflow during
*        this computation.
*
         LRMIN = 1
         DO 30 I = 1, 1 - LEMIN
            LRMIN = DLAMC3( LRMIN*RBASE, ZERO )
   30    CONTINUE
*
*        Finally, call DLAMC5 to compute EMAX and RMAX.
*
         CALL DLAMC5( LBETA, LT, LEMIN, IEEE, LEMAX, LRMAX )
      END IF
*
      BETA = LBETA
      T = LT
      RND = LRND
      EPS = LEPS
      EMIN = LEMIN
      RMIN = LRMIN
      EMAX = LEMAX
      RMAX = LRMAX
*
      RETURN
*
 9999 FORMAT( / / ' WARNING. The value EMIN may be incorrect:-',
     $      '  EMIN = ', I8, /
     $      ' If, after inspection, the value EMIN looks',
     $      ' acceptable please comment out ',
     $      / ' the IF block as marked within the code of routine',
     $      ' DLAMC2,', / ' otherwise supply EMIN explicitly.', / )
*
*     End of DLAMC2
*
      END
*
************************************************************************
*
      DOUBLE PRECISION FUNCTION DLAMC3( A, B )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B
*     ..
*
*  Purpose
*  =======
*
*  DLAMC3  is intended to force  A  and  B  to be stored prior to doing
*  the addition of  A  and  B ,  for use in situations where optimizers
*  might hold one of these in a register.
*
*  Arguments
*  =========
*
*  A, B    (input) DOUBLE PRECISION
*          The values A and B.
*
* =====================================================================
*
*     .. Executable Statements ..
*
      DLAMC3 = A + B
*
      RETURN
*
*     End of DLAMC3
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC4( EMIN, START, BASE )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      INTEGER            BASE, EMIN
      DOUBLE PRECISION   START
*     ..
*
*  Purpose
*  =======
*
*  DLAMC4 is a service routine for DLAMC2.
*
*  Arguments
*  =========
*
*  EMIN    (output) EMIN
*          The minimum exponent before (gradual) underflow, computed by
*          setting A = START and dividing by BASE until the previous A
*          can not be recovered.
*
*  START   (input) DOUBLE PRECISION
*          The starting point for determining EMIN.
*
*  BASE    (input) INTEGER
*          The base of the machine.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I
      DOUBLE PRECISION   A, B1, B2, C1, C2, D1, D2, ONE, RBASE, ZERO
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Executable Statements ..
*
      A = START
      ONE = 1
      RBASE = ONE / BASE
      ZERO = 0
      EMIN = 1
      B1 = DLAMC3( A*RBASE, ZERO )
      C1 = A
      C2 = A
      D1 = A
      D2 = A
*+    WHILE( ( C1.EQ.A ).AND.( C2.EQ.A ).AND.
*    $       ( D1.EQ.A ).AND.( D2.EQ.A )      )LOOP
   10 CONTINUE
      IF( ( C1.EQ.A ) .AND. ( C2.EQ.A ) .AND. ( D1.EQ.A ) .AND.
     $    ( D2.EQ.A ) ) THEN
         EMIN = EMIN - 1
         A = B1
         B1 = DLAMC3( A / BASE, ZERO )
         C1 = DLAMC3( B1*BASE, ZERO )
         D1 = ZERO
         DO 20 I = 1, BASE
            D1 = D1 + B1
   20    CONTINUE
         B2 = DLAMC3( A*RBASE, ZERO )
         C2 = DLAMC3( B2 / RBASE, ZERO )
         D2 = ZERO
         DO 30 I = 1, BASE
            D2 = D2 + B2
   30    CONTINUE
         GO TO 10
      END IF
*+    END WHILE
*
      RETURN
*
*     End of DLAMC4
*
      END
*
************************************************************************
*
      SUBROUTINE DLAMC5( BETA, P, EMIN, IEEE, EMAX, RMAX )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      LOGICAL            IEEE
      INTEGER            BETA, EMAX, EMIN, P
      DOUBLE PRECISION   RMAX
*     ..
*
*  Purpose
*  =======
*
*  DLAMC5 attempts to compute RMAX, the largest machine floating-point
*  number, without overflow.  It assumes that EMAX + abs(EMIN) sum
*  approximately to a power of 2.  It will fail on machines where this
*  assumption does not hold, for example, the Cyber 205 (EMIN = -28625,
*  EMAX = 28718).  It will also fail if the value supplied for EMIN is
*  too large (i.e. too close to zero), probably with overflow.
*
*  Arguments
*  =========
*
*  BETA    (input) INTEGER
*          The base of floating-point arithmetic.
*
*  P       (input) INTEGER
*          The number of base BETA digits in the mantissa of a
*          floating-point value.
*
*  EMIN    (input) INTEGER
*          The minimum exponent before (gradual) underflow.
*
*  IEEE    (input) LOGICAL
*          A logical flag specifying whether or not the arithmetic
*          system is thought to comply with the IEEE standard.
*
*  EMAX    (output) INTEGER
*          The largest exponent before overflow
*
*  RMAX    (output) DOUBLE PRECISION
*          The largest machine floating-point number.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      INTEGER            EXBITS, EXPSUM, I, LEXP, NBITS, TRY, UEXP
      DOUBLE PRECISION   OLDY, RECBAS, Y, Z
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMC3
      EXTERNAL           DLAMC3
*     ..
*     .. Intrinsic Functions ..
      !INTRINSIC          MOD
*     ..
*     .. Executable Statements ..
*
*     First compute LEXP and UEXP, two powers of 2 that bound
*     abs(EMIN). We then assume that EMAX + abs(EMIN) will sum
*     approximately to the bound that is closest to abs(EMIN).
*     (EMAX is the exponent of the required number RMAX).
*
      LEXP = 1
      EXBITS = 1
   10 CONTINUE
      TRY = LEXP*2
      IF( TRY.LE.( -EMIN ) ) THEN
         LEXP = TRY
         EXBITS = EXBITS + 1
         GO TO 10
      END IF
      IF( LEXP.EQ.-EMIN ) THEN
         UEXP = LEXP
      ELSE
         UEXP = TRY
         EXBITS = EXBITS + 1
      END IF
*
*     Now -LEXP is less than or equal to EMIN, and -UEXP is greater
*     than or equal to EMIN. EXBITS is the number of bits needed to
*     store the exponent.
*
      IF( ( UEXP+EMIN ).GT.( -LEXP-EMIN ) ) THEN
         EXPSUM = 2*LEXP
      ELSE
         EXPSUM = 2*UEXP
      END IF
*
*     EXPSUM is the exponent range, approximately equal to
*     EMAX - EMIN + 1 .
*
      EMAX = EXPSUM + EMIN - 1
      NBITS = 1 + EXBITS + P
*
*     NBITS is the total number of bits needed to store a
*     floating-point number.
*
      IF( ( MOD( NBITS, 2 ).EQ.1 ) .AND. ( BETA.EQ.2 ) ) THEN
*
*        Either there are an odd number of bits used to store a
*        floating-point number, which is unlikely, or some bits are
*        not used in the representation of numbers, which is possible,
*        (e.g. Cray machines) or the mantissa has an implicit bit,
*        (e.g. IEEE machines, Dec Vax machines), which is perhaps the
*        most likely. We have to assume the last alternative.
*        If this is true, then we need to reduce EMAX by one because
*        there must be some way of representing zero in an implicit-bit
*        system. On machines like Cray, we are reducing EMAX by one
*        unnecessarily.
*
         EMAX = EMAX - 1
      END IF
*
      IF( IEEE ) THEN
*
*        Assume we are on an IEEE machine which reserves one exponent
*        for infinity and NaN.
*
         EMAX = EMAX - 1
      END IF
*
*     Now create RMAX, the largest machine number, which should
*     be equal to (1.0 - BETA**(-P)) * BETA**EMAX .
*
*     First compute 1.0 - BETA**(-P), being careful that the
*     result is less than 1.0 .
*
      RECBAS = ONE / BETA
      Z = BETA - ONE
      Y = ZERO
      DO 20 I = 1, P
         Z = Z*RECBAS
         IF( Y.LT.ONE )
     $      OLDY = Y
         Y = DLAMC3( Y, Z )
   20 CONTINUE
      IF( Y.GE.ONE )
     $   Y = OLDY
*
*     Now multiply by BETA**EMAX to get RMAX.
*
      DO 30 I = 1, EMAX
         Y = DLAMC3( Y*BETA, ZERO )
   30 CONTINUE
*
      RMAX = Y
      RETURN
*
*     End of DLAMC5
*
      END

      double precision function dznrm2( n, zx, incx)
      logical imag, scale
      integer i, incx, ix, n, next
      double precision cutlo, cuthi, hitest, sum, xmax, absx, zero, one
      double complex      zx(n)
      double precision dreal,dimag
      double complex zdumr,zdumi
      dreal(zdumr) = zdumr
      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
      data         zero, one /0.0d0, 1.0d0/
c
c     unitary norm of the complex n-vector stored in zx() with storage
c     increment incx .
c     if    n .le. 0 return with result = 0.
c     if n .ge. 1 then incx must be .ge. 1
c
c           c.l.lawson , 1978 jan 08
c     modified to correct problem with negative increment, 8/21/90.
c
c     four phase method     using two built-in constants that are
c     hopefully applicable to all machines.
c         cutlo = maximum of  sqrt(u/eps)  over all known machines.
c         cuthi = minimum of  sqrt(v)      over all known machines.
c     where
c         eps = smallest no. such that eps + 1. .gt. 1.
c         u   = smallest positive no.   (underflow limit)
c         v   = largest  no.            (overflow  limit)
c
c     brief outline of algorithm..
c
c     phase 1    scans zero components.
c     move to phase 2 when a component is nonzero and .le. cutlo
c     move to phase 3 when a component is .gt. cutlo
c     move to phase 4 when a component is .ge. cuthi/m
c     where m = n for x() real and m = 2*n for complex.
c
c     values for cutlo and cuthi..
c     from the environmental parameters listed in the imsl converter
c     document the limiting values are as follows..
c     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
c                   univac and dec at 2**(-103)
c                   thus cutlo = 2**(-51) = 4.44089e-16
c     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
c                   thus cuthi = 2**(63.5) = 1.30438e19
c     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
c                   thus cutlo = 2**(-33.5) = 8.23181d-11
c     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
c     data cutlo, cuthi / 8.232d-11,  1.304d19 /
c     data cutlo, cuthi / 4.441e-16,  1.304e19 /
      data cutlo, cuthi / 8.232d-11,  1.304d19 /
c
      if(n .gt. 0) go to 10
         dznrm2  = zero
         go to 300
c
   10 assign 30 to next
      sum = zero
      i = 1
      if( incx .lt. 0 )i = (-n+1)*incx + 1
c                                                 begin main loop
      do 220 ix = 1,n
         absx = dabs(dreal(zx(i)))
         imag = .false.
         go to next,(30, 50, 70, 90, 110)
   30 if( absx .gt. cutlo) go to 85
      assign 50 to next
      scale = .false.
c
c                        phase 1.  sum is zero
c
   50 if( absx .eq. zero) go to 200
      if( absx .gt. cutlo) go to 85
c
c                                prepare for phase 2.
      assign 70 to next
      go to 105
c
c                                prepare for phase 4.
c
  100 assign 110 to next
      sum = (sum / absx) / absx
  105 scale = .true.
      xmax = absx
      go to 115
c
c                   phase 2.  sum is small.
c                             scale to avoid destructive underflow.
c
   70 if( absx .gt. cutlo ) go to 75
c
c                     common code for phases 2 and 4.
c                     in phase 4 sum is large.  scale to avoid overflow.
c
  110 if( absx .le. xmax ) go to 115
         sum = one + sum * (xmax / absx)**2
         xmax = absx
         go to 200
c
  115 sum = sum + (absx/xmax)**2
      go to 200
c
c
c                  prepare for phase 3.
c
   75 sum = (sum * xmax) * xmax
c
   85 assign 90 to next
      scale = .false.
c
c     for real or d.p. set hitest = cuthi/n
c     for complex      set hitest = cuthi/(2*n)
c
      hitest = cuthi/float( 2*n )
c
c                   phase 3.  sum is mid-range.  no scaling.
c
   90 if(absx .ge. hitest) go to 100
         sum = sum + absx**2
  200 continue
c                  control selection of real and imaginary parts.
c
      if(imag) go to 210
         absx = dabs(dimag(zx(i)))
         imag = .true.
      go to next,(  50, 70, 90, 110 )
c
  210 continue
      i = i + incx
  220 continue
c
c              end of main loop.
c              compute square root and adjust for scaling.
c
      dznrm2 = dsqrt(sum)
      if(scale) dznrm2 = dznrm2 * xmax
  300 continue
      return
      end

      LOGICAL          FUNCTION LSAME( CA, CB )
*
*  -- LAPACK auxiliary routine (version 1.1) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          CA, CB
*     ..
*
*  Purpose
*  =======
*
*  LSAME returns .TRUE. if CA is the same letter as CB regardless of
*  case.
*
*  Arguments
*  =========
*
*  CA      (input) CHARACTER*1
*  CB      (input) CHARACTER*1
*          CA and CB specify the single characters to be compared.
*
*     .. Intrinsic Functions ..
      !INTRINSIC          ICHAR
*     ..
*     .. Local Scalars ..
      INTEGER            INTA, INTB, ZCODE
*     ..
*     .. Executable Statements ..
*
*     Test if the characters are equal
*
      LSAME = CA.EQ.CB
      IF( LSAME )
     $   RETURN
*
*     Now test for equivalence if both characters are alphabetic.
*
      ZCODE = ICHAR( 'Z' )
*
*     Use 'Z' rather than 'A' so that ASCII can be detected on Prime
*     machines, on which ICHAR returns a value with bit 8 set.
*     ICHAR('A') on Prime machines returns 193 which is the same as
*     ICHAR('A') on an EBCDIC machine.
*
      INTA = ICHAR( CA )
      INTB = ICHAR( CB )
*
      IF( ZCODE.EQ.90 .OR. ZCODE.EQ.122 ) THEN
*
*        ASCII is assumed - ZCODE is the ASCII code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.97 .AND. INTA.LE.122 ) INTA = INTA - 32
         IF( INTB.GE.97 .AND. INTB.LE.122 ) INTB = INTB - 32
*
      ELSE IF( ZCODE.EQ.233 .OR. ZCODE.EQ.169 ) THEN
*
*        EBCDIC is assumed - ZCODE is the EBCDIC code of either lower or
*        upper case 'Z'.
*
         IF( INTA.GE.129 .AND. INTA.LE.137 .OR.
     $       INTA.GE.145 .AND. INTA.LE.153 .OR.
     $       INTA.GE.162 .AND. INTA.LE.169 ) INTA = INTA + 64
         IF( INTB.GE.129 .AND. INTB.LE.137 .OR.
     $       INTB.GE.145 .AND. INTB.LE.153 .OR.
     $       INTB.GE.162 .AND. INTB.LE.169 ) INTB = INTB + 64
*
      ELSE IF( ZCODE.EQ.218 .OR. ZCODE.EQ.250 ) THEN
*
*        ASCII is assumed, on Prime machines - ZCODE is the ASCII code
*        plus 128 of either lower or upper case 'Z'.
*
         IF( INTA.GE.225 .AND. INTA.LE.250 ) INTA = INTA - 32
         IF( INTB.GE.225 .AND. INTB.LE.250 ) INTB = INTB - 32
      END IF
      LSAME = INTA.EQ.INTB
*
*     RETURN
*
*     End of LSAME
*
      END
**********************************************************************
C
C     Copyright (C) 1992  Roland W. Freund and Noel M. Nachtigal
C     All rights reserved.
C
C     This code is part of a copyrighted package.  For details, see the
C     file `cpyrit.doc' in the top-level directory.
C
C     *****************************************************************
C     ANY USE OF  THIS CODE CONSTITUTES ACCEPTANCE OF  THE TERMS OF THE
C     COPYRIGHT NOTICE
C     *****************************************************************
C
C**********************************************************************
C
      SUBROUTINE ZAXPBY (N,ZZ,ZA,ZX,ZB,ZY)
      use mod_system, only : Rkind

C
C     Purpose:
C     This subroutine computes ZZ = ZA * ZX + ZB * ZY.  Several special
C     cases are handled separately:
C     ZA =  0.0, ZB =  0.0 => ZZ = 0.0
C     ZA =  0.0, ZB =  1.0 => ZZ = ZY  (this is COPY)
C     ZA =  0.0, ZB = -1.0 => ZZ = -ZY
C     ZA =  0.0, ZB =   ZB => ZZ = ZB * ZY  (this is SCAL)
C     ZA =  1.0, ZB =  0.0 => ZZ = ZX  (this is COPY)
C     ZA =  1.0, ZB =  1.0 => ZZ = ZX + ZY
C     ZA =  1.0, ZB = -1.0 => ZZ = ZX - ZY
C     ZA =  1.0, ZB =   ZB => ZZ = ZX + ZB * ZY (this is AXPY)
C     ZA = -1.0, ZB =  0.0 => ZZ = -ZX
C     ZA = -1.0, ZB =  1.0 => ZZ = -ZX + ZY
C     ZA = -1.0, ZB = -1.0 => ZZ = -ZX - ZY
C     ZA = -1.0, ZB =   ZB => ZZ = -ZX + ZB * ZY
C     ZA =   ZA, ZB =  0.0 => ZZ = ZA * ZX  (this is SCAL)
C     ZA =   ZA, ZB =  1.0 => ZZ = ZA * ZX + ZY  (this is AXPY)
C     ZA =   ZA, ZB = -1.0 => ZZ = ZA * ZX - ZY
C     ZA =   ZA, ZB =   ZB => ZZ = ZA * ZX + ZB * ZY
C     ZZ may be the same as ZX or ZY.
C
C     Parameters:
C     N  = the dimension of the vectors (input).
C     ZZ = the vector result (output).
C     ZA = scalar multiplier for ZX (input).
C     ZX = one of the vectors (input).
C     ZB = scalar multiplier for ZY (input).
C     ZY = the other vector (input).
C
C     Noel M. Nachtigal
C     March 23, 1993
C
C**********************************************************************
C
      !INTRINSIC DIMAG, DREAL
C
      INTEGER N
      DOUBLE COMPLEX ZA, ZB, ZX(N), ZY(N), ZZ(N)
C
C     Local variables.
C
      INTEGER I
      DOUBLE PRECISION DAI, DAR, DBI, DBR
C
      IF (N.LE.0) RETURN
C
      DAI = aimag(ZA)
      DAR = real(ZA,kind=Rkind)
      DBI = aimag(ZB)
      DBR = real(ZB,kind=Rkind)
      IF ((DAR.EQ.0.0D0).AND.(DAI.EQ.0.0D0)) THEN
         IF ((DBR.EQ.0.0D0).AND.(DBI.EQ.0.0D0)) THEN
C     ZA = 0.0, ZB = 0.0 => ZZ = 0.0.
            DO 10 I = 1, N
               ZZ(I) = (0.0D0,0.0D0)
 10         CONTINUE
         ELSE IF ((DBR.EQ.1.0D0).AND.(DBI.EQ.0.0D0)) THEN
C     ZA = 0.0, ZB = 1.0 => ZZ = ZY (this is COPY).
            DO 20 I = 1, N
               ZZ(I) = ZY(I)
 20         CONTINUE
         ELSE IF ((DBR.EQ.-1.0D0).AND.(DBI.EQ.0.0D0)) THEN
C     ZA = 0.0, ZB = -1.0 => ZZ = -ZY.
            DO 30 I = 1, N
               ZZ(I) = -ZY(I)
 30         CONTINUE
         ELSE
C     ZA = 0.0, ZB = ZB => ZZ = ZB * ZY (this is SCAL).
            DO 40 I = 1, N
               ZZ(I) = ZB * ZY(I)
 40         CONTINUE
         END IF
      ELSE IF ((DAR.EQ.1.0D0).AND.(DAI.EQ.0.0D0)) THEN
         IF ((DBR.EQ.0.0D0).AND.(DBI.EQ.0.0D0)) THEN
C     ZA = 1.0, ZB = 0.0 => ZZ = ZX (this is COPY).
            DO 50 I = 1, N
               ZZ(I) = ZX(I)
 50         CONTINUE
         ELSE IF ((DBR.EQ.1.0D0).AND.(DBI.EQ.0.0D0)) THEN
C     ZA = 1.0, ZB = 1.0 => ZZ = ZX + ZY.
            DO 60 I = 1, N
               ZZ(I) = ZX(I) + ZY(I)
 60         CONTINUE
         ELSE IF ((DBR.EQ.-1.0D0).AND.(DBI.EQ.0.0D0)) THEN
C     ZA = 1.0, ZB = -1.0 => ZZ = ZX - ZY.
            DO 70 I = 1, N
               ZZ(I) = ZX(I) - ZY(I)
 70         CONTINUE
         ELSE
C     ZA = 1.0, ZB = ZB => ZZ = ZX + ZB * ZY (this is AXPY).
            DO 80 I = 1, N
               ZZ(I) = ZX(I) + ZB * ZY(I)
 80         CONTINUE
         END IF
      ELSE IF ((DAR.EQ.-1.0D0).AND.(DAI.EQ.0.0D0)) THEN
         IF ((DBR.EQ.0.0D0).AND.(DBI.EQ.0.0D0)) THEN
C     ZA = -1.0, ZB = 0.0 => ZZ = -ZX
            DO 90 I = 1, N
               ZZ(I) = -ZX(I)
 90         CONTINUE
         ELSE IF ((DBR.EQ.1.0D0).AND.(DBI.EQ.0.0D0)) THEN
C     ZA = -1.0, ZB = 1.0 => ZZ = -ZX + ZY
            DO 100 I = 1, N
               ZZ(I) = -ZX(I) + ZY(I)
 100        CONTINUE
         ELSE IF ((DBR.EQ.-1.0D0).AND.(DBI.EQ.0.0D0)) THEN
C     ZA = -1.0, ZB = -1.0 => ZZ = -ZX - ZY.
            DO 110 I = 1, N
               ZZ(I) = -ZX(I) - ZY(I)
 110        CONTINUE
         ELSE
C     ZA = -1.0, ZB = ZB => ZZ = -ZX + ZB * ZY
            DO 120 I = 1, N
               ZZ(I) = -ZX(I) + ZB * ZY(I)
 120        CONTINUE
         END IF
      ELSE
         IF ((DBR.EQ.0.0D0).AND.(DBI.EQ.0.0D0)) THEN
C     ZA = ZA, ZB = 0.0 => ZZ = ZA * ZX (this is SCAL).
            DO 130 I = 1, N
               ZZ(I) = ZA * ZX(I)
 130        CONTINUE
         ELSE IF ((DBR.EQ.1.0D0).AND.(DBI.EQ.0.0D0)) THEN
C     ZA = ZA, ZB = 1.0 => ZZ = ZA * ZX + ZY (this is AXPY)
            DO 140 I = 1, N
               ZZ(I) = ZA * ZX(I) + ZY(I)
 140        CONTINUE
         ELSE IF ((DBR.EQ.-1.0D0).AND.(DBI.EQ.0.0D0)) THEN
C     ZA = ZA, ZB = -1.0 => ZZ = ZA * ZX - ZY.
            DO 150 I = 1, N
               ZZ(I) = ZA * ZX(I) - ZY(I)
 150        CONTINUE
         ELSE
C     ZA = ZA, ZB = ZB => ZZ = ZA * ZX + ZB * ZY.
            DO 160 I = 1, N
               ZZ(I) = ZA * ZX(I) + ZB * ZY(I)
 160        CONTINUE
         END IF
      END IF
C
      RETURN
      END
C
C**********************************************************************
      double complex function zdotu(n,zx,incx,zy,incy)
c
c     forms the dot product of a vector.
c     jack dongarra, 3/11/78.
c
      double complex zx(n),zy(n),ztemp
      integer i,ix,iy,n,incx,incy
      ztemp = (0.0d0,0.0d0)
      zdotu = (0.0d0,0.0d0)
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ztemp = ztemp + zx(ix)*zy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      zdotu = ztemp
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        ztemp = ztemp + zx(i)*zy(i)
   30 continue
      zdotu = ztemp
      return
      end
C**********************************************************************
C
      SUBROUTINE ZRANDN (N,ZX,SEED)
      use mod_system, only : Rkind

C
C     Purpose:
C     Fills the vector ZX with random numbers  between 0 and 1.  If the
C     SEED is given, it should be odd and positive.  The generator is a
C     fairly unsophisticated one, from Pearson's  "Numerical methods in
C     engineering and science" book.
C
C     Parameters:
C     N    = the dimension of the vector (input).
C     ZX   = the vector to fill with random numbers (output).
C     SEED = the seed for the generator (input).
C
C     Noel M. Nachtigal
C     April 23, 1993
C
C**********************************************************************
C
      !INTRINSIC DBLE, DCMPLX, IABS, MOD
C
      INTEGER N, SEED
      DOUBLE COMPLEX ZX(N)
C
C     Local variables.
C
      INTEGER I, J
      DOUBLE PRECISION IMAGX, REALX
C
C     Local variables that are saved from one call to the next.
C
      DOUBLE PRECISION DMAX
      INTEGER IM, IMAX, IS
      SAVE DMAX, IM, IMAX, IS
      DATA IM/0/
C
C     Initialize the generator data.
C
      IF (IM.EQ.0) THEN
         J  = 0
         IM = 1
         DO 10 I = 1, 31
            J = J + 1
            IF (IM*2.LE.IM) GO TO 20
            IM = IM * 2
 10      CONTINUE
 20      IMAX = (IM-1) * 2 + 1
         DMAX = DBLE(IMAX)
         DO 30 I = 1, MOD(J,3)
            J = J - 1
            IM = IM / 2
 30      CONTINUE
         IM = IM + 5
         IS = IABS(MOD(IM*30107,IMAX))
      END IF
C
C     Check whether there is a new seed.
C
      IF (SEED.GT.0) IS = (SEED / 2) * 2 + 1
C
C     Here goes the rest.
C
      DO 40 I = 1, N
         REALX = DBLE(IS) / DMAX
         IS    = IABS(MOD(IM*IS,IMAX))
         IMAGX = DBLE(IS) / DMAX
         IS    = IABS(MOD(IM*IS,IMAX))
         ZX(I) = cmplx(REALX,IMAGX,kind=Rkind)
 40   CONTINUE
C
      RETURN
      END
C
C**********************************************************************
      subroutine zrotg(ca,cb,c,s)
      use mod_system, only : Rkind
      double complex ca,cb,s
      double precision c
      double precision norm,scale
      double complex alpha
      if (abs(ca) .ne. 0.0d0) go to 10
         c = 0.0d0
         s = (1.0d0,0.0d0)
         ca = cb
         go to 20
   10 continue
         scale = abs(ca) + abs(cb)
         norm  = scale*sqrt((abs(ca/cmplx(scale,0.0d0,kind=Rkind)))**2 +
     *                      (abs(cb/cmplx(scale,0.0d0,kind=Rkind)))**2)
         alpha = ca /abs(ca)
         c = abs(ca) / norm
         s = alpha * conjg(cb) / norm
         ca = alpha * norm
   20 continue
      return
      end
C**********************************************************************
C
C     Copyright (C) 1992  Roland W. Freund and Noel M. Nachtigal
C     All rights reserved.
C
C     This code is part of a copyrighted package.  For details, see the
C     file `cpyrit.doc' in the top-level directory.
C
C     *****************************************************************
C     ANY USE OF  THIS CODE CONSTITUTES ACCEPTANCE OF  THE TERMS OF THE
C     COPYRIGHT NOTICE
C     *****************************************************************
C
C**********************************************************************
C
C     This  file  contains  the  routine  for  the  QMR  algorithm  for
C     symmetric matrices,  using the  three-term recurrence  variant of
C     the Lanczos algorithm without look-ahead.
C
C**********************************************************************
C
      SUBROUTINE ZSQMX (NDIM,NLEN,NLIM,VECS,TOL,INFO)
      use mod_system, only : Rkind

C
C     Purpose:
C     This subroutine uses  the QMR algorithm to solve  linear systems.
C     It runs  the algorithm  to convergence or until  a user-specified
C     limit on the number of  iterations is reached.  It is set  up  to
C     solve symmetric systems starting with identical starting vectors.
C
C     The  code is  set up  to solve  the system  A x = b  with initial
C     guess x_0 = 0.  Here  A x = b  denotes the preconditioned system,
C     and  it  is  connected with the  original system as follows.  Let
C     B y = c be the original unpreconditioned system to be solved, and
C     let y_0 be an arbitrary initial guess for its solution.  Then:
C     A x = b, where  A = M_1^{-1} B M_2^{-1},
C     x = M_2 (y - y_0), b = M_1^{-1} (c - B y_0).
C     Here M = M_1 M_2 is the preconditioner.
C
C     To recover the final iterate y_n for  the original system B y = c
C     from the final iterate x_n for the preconditioned system A x = b,
C     set
C     y_n = y_0 + M_2^{-1} x_n.
C
C     The implementation does not have look-ahead, so it is less robust
C     than the full version.
C
C     Parameters:
C     For a description of  the parameters, see the file `zsqmx.doc' in
C     the current directory.
C
C     External routines used:
C     double precision dlamch(ch)
C     LAPACK routine, computes machine-related constants.
C     double precision dznrm2(n,x,incx)
C     BLAS-1 routine, computes the 2-norm of x.
C     subroutine zaxpby(n,z,a,x,b,y)
C     Library routine, computes z = a * x + b * y.
C     double precision zdotu(n,x,incx,y,incy)
C     BLAS-1 routine, computes y' * x.
C     double precision zqmxom(n)
C     User-supplied routine, specifies the QMR scaling factors.
C     subroutine zrandn(n,x,seed)
C     Library routine, fills x with random numbers.
C     subroutine zrotg(a,b,cos,sin)
C     BLAS-1 routine, computes the Givens rotation which rotates the
C     vector [a; b] into [ sqrt(a**2 + b**2); 0 ].
C
C     Noel M. Nachtigal
C     May 25, 1993
C
C**********************************************************************
C
      !INTRINSIC CDABS, DABS, DBLE, DCMPLX, DCONJG, DMAX1, DSQRT, MAX0
      !INTRINSIC MOD
      EXTERNAL DLAMCH, DZNRM2, ZAXPBY, ZDOTU, ZRANDN, ZROTG, ZSQMXO
      DOUBLE COMPLEX ZDOTU
      DOUBLE PRECISION DLAMCH, DZNRM2, ZSQMXO
C
      INTEGER INFO(4), NDIM, NLEN, NLIM
      DOUBLE COMPLEX VECS(NDIM,7)
      DOUBLE PRECISION TOL
C
C     Miscellaneous parameters.
C
      DOUBLE COMPLEX ZONE, ZZERO
      PARAMETER (ZONE = (1.0D0,0.0D0),ZZERO = (0.0D0,0.0D0))
      DOUBLE PRECISION DHUN, DONE, DTEN, DZERO
      PARAMETER (DHUN = 1.0D2,DONE = 1.0D0,DTEN = 1.0D1,DZERO = 0.0D0)
C
C     Local variables, permanent.
C
      INTEGER IERR, ISN, ISNM1, ISNM2, IVN, IVNM1, IVNP1, IWN, IWNM1
      SAVE    IERR, ISN, ISNM1, ISNM2, IVN, IVNM1, IVNP1, IWN, IWNM1
      INTEGER IWNP1, N, RETLBL, TF, TRES, VF
      SAVE    IWNP1, N, RETLBL, TF, TRES, VF
      DOUBLE COMPLEX DN, DNM1, DNP1, RHSN, SCSN, SCSNM1, SINN, SINNM1
      SAVE           DN, DNM1, DNP1, RHSN, SCSN, SCSNM1, SINN, SINNM1
      DOUBLE PRECISION COSN, COSNM1, CSIN, CSINP1, MAXOMG, OMGN
      SAVE             COSN, COSNM1, CSIN, CSINP1, MAXOMG, OMGN
      DOUBLE PRECISION OMGNP1, R0, RESN, RHON, RHONP1, SCVN, SCVNP1
      SAVE             OMGNP1, R0, RESN, RHON, RHONP1, SCVN, SCVNP1
      DOUBLE PRECISION SCWN, SCWNP1, TMAX, TMIN, TNRM, UCHK, UNRM
      SAVE             SCWN, SCWNP1, TMAX, TMIN, TNRM, UCHK, UNRM
C
C     Local variables, transient.
C
      INTEGER INIT, REVCOM
      DOUBLE COMPLEX RHN, RHNM1, RHNM2, RHNP1, RHSNP1, ZTMP
      DOUBLE PRECISION DTMP, SCVNM1, SCWNM1
C
C     Initialize some of the permanent variables.
C
      DATA RETLBL /0/
C
C     Check the reverse communication flag to see where to branch.
C     REVCOM   RETLBL      Comment
C     0        0    first call, go to label 10
C     1       30    returning from AXB, go to label 30
C     1       60    returning from AXB, go to label 60
C     2       40    returning from ATXB, go to label 40
C
      REVCOM  = INFO(2)
      INFO(2) = 0
      IF (REVCOM.EQ.0) THEN
         N = 0
         IF (RETLBL.EQ.0) GO TO 10
      ELSE IF (REVCOM.EQ.1) THEN
         IF (RETLBL.EQ.30) THEN
            GO TO 30
         ELSE IF (RETLBL.EQ.60) THEN
            GO TO 60
         END IF
      ELSE IF (REVCOM.EQ.2) THEN
         IF (RETLBL.EQ.40) GO TO 40
      END IF
      IERR = 1
      GO TO 80
C
C     Check whether the inputs are valid.
C
 10   IERR = 0
      IF (NDIM.LT.1)        IERR = 2
      IF (NLEN.LT.1)        IERR = 2
      IF (NLIM.LT.1)        IERR = 2
      IF (NLEN.GT.NDIM)     IERR = 2
      IF (IERR.NE.0) GO TO 80
C
C     Extract from INFO the output units TF and VF, the true residual
C     flag TRES, and the left starting vector flag INIT.
C
      VF   = MAX0(INFO(1),0)
      INIT = VF / 100000
      VF   = VF - INIT * 100000
      TRES = VF / 10000
      VF   = VF - TRES * 10000
      TF   = VF / 100
      VF   = VF - TF * 100
C
C     Extract and check the various tolerances.
C
      TNRM = DLAMCH('E') * DTEN
      TMIN = DSQRT(DSQRT(DLAMCH('S')))
      TMAX = DONE / TMIN
      IF (TOL.LE.DZERO) TOL = DSQRT(DLAMCH('E'))
C
C     Start the trace messages and convergence history.
C
      IF (VF.NE.0) WRITE (VF,'(I8,2E11.4)') 0, DONE, DONE
      IF (TF.NE.0) WRITE (TF,'(I8,2E11.4)') 0, DONE, DONE
C
C     Initialize the wrapped indices.
C
      ISNM1 = 5
      ISN   = ISNM1
      IVN   = 3
      IVNP1 = IVN
      IWN   = 3
      IWNP1 = IWN
C
C     Set x_0 = 0 and compute the norm of the initial residual.
C
      CALL ZAXPBY (NLEN,VECS(1,IVN),ZONE,VECS(1,2),ZZERO,VECS(1,IVN))
      CALL ZAXPBY (NLEN,VECS(1,1),ZZERO,VECS(1,1),ZZERO,VECS(1,1))
      R0 = DZNRM2(NLEN,VECS(1,IVN),1)
      IF ((TOL.GE.DONE).OR.(R0.EQ.DZERO)) GO TO 80
C
C     Check whether the auxiliary vector must be supplied.
C
C     SYM  IF (INIT.EQ.0) CALL ZRANDN (NLEN,VECS(1,IWN),1)
C
C     Compute scale factors and check for invariant subspaces.
C
      SCVNP1 = R0
C     SYM  SCWNP1 = DZNRM2(NLEN,VECS(1,IWN),1)
      SCWNP1 = SCVNP1
      IF (SCVNP1.LT.TNRM) IERR = IERR + 16
      IF (SCWNP1.LT.TNRM) IERR = IERR + 32
      IF (IERR.NE.0) GO TO 80
      DNP1 = ZDOTU(NLEN,VECS(1,IWN),1,VECS(1,IVN),1) / ( SCVNP1 * SCWNP1
     $   )
      IF ((SCVNP1.GE.TMAX).OR.(SCVNP1.LE.TMIN)) THEN
         ZTMP = cmplx(DONE / SCVNP1,DZERO,kind=Rkind)
         CALL ZAXPBY (NLEN,VECS(1,IVN),ZTMP,VECS(1,IVN),ZZERO,VECS(1,IVN
     $      ))
         SCVNP1 = DONE
      END IF
      IF ((SCWNP1.GE.TMAX).OR.(SCWNP1.LE.TMIN)) THEN
         ZTMP = cmplx(DONE / SCWNP1,DZERO,kind=Rkind)
C     SYM     CALL ZAXPBY (NLEN,VECS(1,IWN),ZTMP,VECS(1,IWN),ZZERO,VECS(1,IWN))
         SCWNP1 = DONE
      END IF
      RHONP1 = SCVNP1
      CSINP1 = SCWNP1
      SCVNP1 = DONE / SCVNP1
      SCWNP1 = DONE / SCWNP1
C
C     Initialize the variables.
C
      N      = 1
      DN     = ZONE
      COSN   = DONE
      RESN   = DONE
      COSNM1 = DONE
      OMGN   = DZERO
      SCSN   = DZERO
      SCVN   = DZERO
      SCWN   = DZERO
      SCSNM1 = DZERO
      SINN   = ZZERO
      SINNM1 = ZZERO
      OMGNP1 = ZSQMXO(N)
      RHSN   = OMGNP1 * R0
      MAXOMG = DONE / OMGNP1
C
C     This is one step of the classical Lanczos algorithm.
C
 20   IVNM1 = IVN
      IVN   = IVNP1
      IVNP1 = MOD(N,2) + 3
      IWNM1 = IWN
      IWN   = IWNP1
      IWNP1 = MOD(N,2) + 3
C
C     Check whether D_n is nonsingular.
C
      DNM1 = DN
      DN   = DNP1
      IF (abs(DN).EQ.DZERO) THEN
         IERR = 8
         GO TO 80
      END IF
C
C     Have the caller carry out AXB, then return here.
C     CALL AXB (VECS(1,IVN),VECS(1,7))
C
      INFO(2) = 1
      INFO(3) = IVN
      INFO(4) = 7
      RETLBL  = 30
      RETURN
 30   RETLBL = 0
C
C     Compute H_{n-1,n} and build part of the vector v_{n+1}.
C
      SCVNM1 = SCVN
      CSIN   = CSINP1
      SCVN   = SCVNP1
      ZTMP   = CSIN * DN / DNM1 * SCVNM1 / SCVN
      CALL ZAXPBY (NLEN,VECS(1,IVNP1),ZONE,VECS(1,7),-ZTMP,VECS(1,IVNM1)
     $   )
C
C     Have the caller carry out ATXB, then return here.
C     CALL ATXB (VECS(1,IWN),VECS(1,7))
C
      INFO(2) = 2
      INFO(3) = IWN
      INFO(4) = 7
      RETLBL  = 40
C     SYM  RETURN
 40   RETLBL = 0
C
C     Build part of the vector w_{n+1}.
C
      SCWNM1 = SCWN
      RHON   = RHONP1
      SCWN   = SCWNP1
      ZTMP   = RHON * DN / DNM1 * SCWNM1 / SCWN
C     SYM  CALL ZAXPBY (NLEN,VECS(1,IWNP1),ZONE,VECS(1,7),-ZTMP,VECS(1,IWNM1))
C
C     Compute H_{nn} and finish the new vectors.
C
      RHN = SCVN * SCWN * ZDOTU(NLEN,VECS(1,IWN),1,VECS(1,IVNP1),1) / DN
      CALL ZAXPBY (NLEN,VECS(1,IVNP1),ZONE,VECS(1,IVNP1),-RHN,VECS(1,IVN
     $   ))
C     SYM CALL ZAXPBY (NLEN,VECS(1,IWNP1),ZONE,VECS(1,IWNP1),-RHN,VECS(1,IWN))
C
C     Compute scale factors and check for invariant subspaces.
C
      IERR   = 0
      SCVNP1 = DZNRM2(NLEN,VECS(1,IVNP1),1)
C     SYM  SCWNP1 = DZNRM2(NLEN,VECS(1,IWNP1),1)
      SCWNP1 = SCVNP1
      RHONP1 = SCVN * SCVNP1
      CSINP1 = SCWN * SCWNP1
      RHNP1  = RHONP1
      IF (SCVNP1.LT.TNRM) IERR = IERR + 16
      IF (SCWNP1.LT.TNRM) IERR = IERR + 32
      IF (IERR.NE.0) GO TO 50
      DNP1   = ZDOTU(NLEN,VECS(1,IWNP1),1,VECS(1,IVNP1),1) / ( SCVNP1 *
     $   SCWNP1 )
      IF ((SCVNP1.GE.TMAX).OR.(SCVNP1.LE.TMIN)) THEN
         ZTMP = cmplx(DONE / SCVNP1,DZERO,kind=Rkind)
         CALL ZAXPBY (NLEN,VECS(1,IVNP1),ZTMP,VECS(1,IVNP1),ZZERO,VECS(1
     $      ,IVNP1))
         SCVNP1 = DONE
      END IF
      IF ((SCWNP1.GE.TMAX).OR.(SCWNP1.LE.TMIN)) THEN
         ZTMP = cmplx(DONE / SCWNP1,DZERO,kind=Rkind)
C
         SCWNP1 = DONE
      END IF
      SCVNP1 = DONE / SCVNP1
      SCWNP1 = DONE / SCWNP1
C
C     The QMR code starts here.
C     Multiply the new column by the previous omega's.
C     Get the next scaling factor omega(i) and update MAXOMG.
C
 50   RHNM1  = CSIN * DN * OMGN / DNM1
      OMGN   = OMGNP1
      RHN    = OMGN * RHN
      OMGNP1 = ZSQMXO(N+1)
      RHNP1  = OMGNP1 * RHNP1
      MAXOMG = DMAX1(MAXOMG,DONE/OMGN)
C
C     Apply the previous rotations.
C
      RHNM2  = SINNM1 * RHNM1
      RHNM1  = COSNM1 * RHNM1
      COSNM1 = COSN
      SINNM1 = SINN
      ZTMP   = RHNM1
      RHNM1  =  COSNM1 * ZTMP + SINNM1 * RHN
      RHN    = -conjg(SINNM1) * ZTMP + COSNM1 * RHN
C
C     Compute the rotation for the last element (this also applies it).
C
      CALL ZROTG (RHN,RHNP1,COSN,SINN)
C
C     Apply the new rotation to the right-hand side vector.
C
      RHSNP1 = -conjg(SINN) * RHSN
      RHSN   =  COSN * RHSN
C
C     Compute the next search direction s_n.
C
      ISNM2  = ISNM1
      ISNM1  = ISN
      ISN    = MOD(N-1,2) + 5
      ZTMP   = SCSNM1 * RHNM2 / SCVN
      CALL ZAXPBY (NLEN,VECS(1,ISN),ZONE,VECS(1,IVN),-ZTMP,VECS(1,ISNM2)
     $   )
      SCSNM1 = SCSN
      ZTMP   = SCSNM1 * RHNM1 / SCVN
      CALL ZAXPBY (NLEN,VECS(1,ISN),ZONE,VECS(1,ISN),-ZTMP,VECS(1,ISNM1)
     $   )
      SCSN   = SCVN / RHN
C
C     Compute the new QMR iterate, then scale the search direction.
C
      ZTMP = SCSN * RHSN
      CALL ZAXPBY (NLEN,VECS(1,1),ZONE,VECS(1,1),ZTMP,VECS(1,ISN))
      DTMP = abs(SCSN)
      IF ((DTMP.GE.TMAX).OR.(DTMP.LE.TMIN)) THEN
         CALL ZAXPBY (NLEN,VECS(1,ISN),SCSN,VECS(1,ISN),ZZERO,VECS(1,ISN
     $      ))
         SCSN = ZONE
      END IF
C
C     Compute the residual norm upper bound.
C     If the scaled upper bound is within one order of magnitude of the
C     target convergence norm, compute the true residual norm.
C
      RHSN = RHSNP1
      UNRM = DSQRT(DBLE(N+1)) * MAXOMG * abs(RHSNP1) / R0
      UCHK = UNRM
      IF ((TRES.EQ.0).AND.(UNRM/TOL.GT.DTEN).AND.(N.LT.NLIM)) GO TO 70
C
C     Have the caller carry out AXB, then return here.
C     CALL AXB (VECS(1,1),VECS(1,7))
C
      INFO(2) = 1
      INFO(3) = 1
      INFO(4) = 7
      RETLBL  = 60
      RETURN
 60   RETLBL = 0
      CALL ZAXPBY (NLEN,VECS(1,7),ZONE,VECS(1,2),-ZONE,VECS(1,7))
      RESN = DZNRM2(NLEN,VECS(1,7),1) / R0
      UCHK = RESN
C
C     Output the trace messages and convergence history.
C
 70   IF (VF.NE.0) WRITE (VF,'(I8,2E11.4)') N, UNRM, RESN
      IF (TF.NE.0) WRITE (TF,'(I8,2E11.4)') N, UNRM, RESN
C
C     Check for convergence or termination.  Stop if:
C     1. algorithm converged;
C     2. there is an error condition;
C     3. the residual norm upper bound is smaller than the computed
C     residual norm by a factor of at least 100;
C     4. algorithm exceeded the iterations limit.
C
      IF (RESN.LE.TOL) THEN
         IERR = 0
         GO TO 80
      ELSE IF (IERR.NE.0) THEN
         GO TO 80
      ELSE IF (UNRM.LT.UCHK/DHUN) THEN
         IERR = 4
         GO TO 80
      ELSE IF (N.GE.NLIM) THEN
         IERR = 4
         GO TO 80
      END IF
C
C     Update the running counter.
C
      N = N + 1
      GO TO 20
C
C     That's all.
C
 80   NLIM    = N
      RETLBL  = 0
      INFO(1) = IERR
C
      RETURN
      END
C
C**********************************************************************
C
      DOUBLE PRECISION FUNCTION ZSQMXO (I)
C
C     Purpose:
C     Returns the scaling parameter OMEGA(I).
C
C     Parameters:
C     I = the index of the parameter OMEGA (input).
C
C     Noel M. Nachtigal
C     March 30, 1993
C
C**********************************************************************
C
      INTEGER I
C
      ZSQMXO = 1.0D0
C
      RETURN
      END
C
C**********************************************************************




C**********************************************************************
C
C     Copyright (C) 1992  Roland W. Freund and Noel M. Nachtigal
C     All rights reserved.
C
C     This code is part of a copyrighted package.  For details, see the
C     file `cpyrit.doc' in the top-level directory.
C
C     *****************************************************************
C     ANY USE OF  THIS CODE CONSTITUTES ACCEPTANCE OF  THE TERMS OF THE
C                             COPYRIGHT NOTICE
C     *****************************************************************
C
C**********************************************************************
C
C     This  file  contains  the  routine  for  the  QMR  algorithm  for
C     unsymmetric matrices, using the  three-term recurrence variant of
C     the Lanczos algorithm without look-ahead.
C
C**********************************************************************
C
      SUBROUTINE ZUQMX (NDIM,NLEN,NLIM,VECS,TOL,INFO)
      use mod_system, only : Rkind
C
C     Purpose:
C     This subroutine uses  the QMR algorithm to solve  linear systems.
C     It runs  the algorithm  to convergence or until  a user-specified
C     limit on the number of  iterations is reached.
C
C     The  code is  set up  to solve  the system  A x = b  with initial
C     guess x_0 = 0.  Here  A x = b  denotes the preconditioned system,
C     and  it  is  connected with the  original system as follows.  Let
C     B y = c be the original unpreconditioned system to be solved, and
C     let y_0 be an arbitrary initial guess for its solution.  Then:
C          A x = b, where  A = M_1^{-1} B M_2^{-1},
C                          x = M_2 (y - y_0), b = M_1^{-1} (c - B y_0).
C     Here M = M_1 M_2 is the preconditioner.
C
C     To recover the final iterate y_n for  the original system B y = c
C     from the final iterate x_n for the preconditioned system A x = b,
C     set
C               y_n = y_0 + M_2^{-1} x_n.
C
C     The implementation does not have look-ahead, so it is less robust
C     than the full version.
C
C     Parameters:
C     For a description of  the parameters, see the file `zuqmx.doc' in
C     the current directory.
C
C     External routines used:
C     double precision dlamch(ch)
C        LAPACK routine, computes machine-related constants.
C     double precision dznrm2(n,x,incx)
C        BLAS-1 routine, computes the 2-norm of x.
C     subroutine zaxpby(n,z,a,x,b,y)
C        Library routine, computes z = a * x + b * y.
C     double precision zdotu(n,x,incx,y,incy)
C        BLAS-1 routine, computes y' * x.
C     double precision zqmxom(n)
C        User-supplied routine, specifies the QMR scaling factors.
C     subroutine zrandn(n,x,seed)
C        Library routine, fills x with random numbers.
C     subroutine zrotg(a,b,cos,sin)
C        BLAS-1 routine, computes the Givens rotation which rotates the
C        vector [a; b] into [ sqrt(a**2 + b**2); 0 ].
C
C     Noel M. Nachtigal
C     May 25, 1993
C
C**********************************************************************
C
      !INTRINSIC abs, DABS, DBLE, DCMPLX, DCONJG, DMAX1, DSQRT, MAX0
      !INTRINSIC MOD
      EXTERNAL DLAMCH, DZNRM2, ZAXPBY, ZDOTU, ZRANDN, ZROTG, ZUQMXO
      DOUBLE COMPLEX ZDOTU
      DOUBLE PRECISION DLAMCH, DZNRM2, ZUQMXO
C
      INTEGER INFO(4), NDIM, NLEN, NLIM
      DOUBLE COMPLEX VECS(NDIM,9)
      DOUBLE PRECISION TOL
C
C     Miscellaneous parameters.
C
      DOUBLE COMPLEX ZONE, ZZERO
      PARAMETER (ZONE = (1.0D0,0.0D0),ZZERO = (0.0D0,0.0D0))
      DOUBLE PRECISION DHUN, DONE, DTEN, DZERO
      PARAMETER (DHUN = 1.0D2,DONE = 1.0D0,DTEN = 1.0D1,DZERO = 0.0D0)
C
C     Local variables, permanent.
C
      INTEGER IERR, ISN, ISNM1, ISNM2, IVN, IVNM1, IVNP1, IWN, IWNM1
      SAVE    IERR, ISN, ISNM1, ISNM2, IVN, IVNM1, IVNP1, IWN, IWNM1
      INTEGER IWNP1, N, RETLBL, TF, TRES, VF
      SAVE    IWNP1, N, RETLBL, TF, TRES, VF
      DOUBLE COMPLEX DN, DNM1, DNP1, RHSN, SCSN, SCSNM1, SINN, SINNM1
      SAVE           DN, DNM1, DNP1, RHSN, SCSN, SCSNM1, SINN, SINNM1
      DOUBLE PRECISION COSN, COSNM1, CSIN, CSINP1, MAXOMG, OMGN
      SAVE             COSN, COSNM1, CSIN, CSINP1, MAXOMG, OMGN
      DOUBLE PRECISION OMGNP1, R0, RESN, RHON, RHONP1, SCVN, SCVNP1
      SAVE             OMGNP1, R0, RESN, RHON, RHONP1, SCVN, SCVNP1
      DOUBLE PRECISION SCWN, SCWNP1, TMAX, TMIN, TNRM, UCHK, UNRM
      SAVE             SCWN, SCWNP1, TMAX, TMIN, TNRM, UCHK, UNRM
C
C     Local variables, transient.
C
      INTEGER INIT, REVCOM
      DOUBLE COMPLEX RHN, RHNM1, RHNM2, RHNP1, RHSNP1, ZTMP
      DOUBLE PRECISION DTMP, SCVNM1, SCWNM1
C
C     Initialize some of the permanent variables.
C
      DATA RETLBL /0/
C
C     Check the reverse communication flag to see where to branch.
C        REVCOM   RETLBL      Comment
C           0        0    first call, go to label 10
C           1       30    returning from AXB, go to label 30
C           1       60    returning from AXB, go to label 60
C           2       40    returning from ATXB, go to label 40
C
      REVCOM  = INFO(2)
      INFO(2) = 0
      IF (REVCOM.EQ.0) THEN
         N = 0
         IF (RETLBL.EQ.0) GO TO 10
      ELSE IF (REVCOM.EQ.1) THEN
         IF (RETLBL.EQ.30) THEN
            GO TO 30
         ELSE IF (RETLBL.EQ.60) THEN
            GO TO 60
         END IF
      ELSE IF (REVCOM.EQ.2) THEN
         IF (RETLBL.EQ.40) GO TO 40
      END IF
      IERR = 1
      GO TO 80
C
C     Check whether the inputs are valid.
C
 10   IERR = 0
      IF (NDIM.LT.1)        IERR = 2
      IF (NLEN.LT.1)        IERR = 2
      IF (NLIM.LT.1)        IERR = 2
      IF (NLEN.GT.NDIM)     IERR = 2
      IF (IERR.NE.0) GO TO 80
C
C     Extract from INFO the output units TF and VF, the true residual
C     flag TRES, and the left starting vector flag INIT.
C
      VF   = MAX0(INFO(1),0)
      INIT = VF / 100000
      VF   = VF - INIT * 100000
      TRES = VF / 10000
      VF   = VF - TRES * 10000
      TF   = VF / 100
      VF   = VF - TF * 100
C
C     Extract and check the various tolerances.
C
      TNRM = DLAMCH('E') * DTEN
      TMIN = DSQRT(DSQRT(DLAMCH('S')))
      TMAX = DONE / TMIN
      IF (TOL.LE.DZERO) TOL = DSQRT(DLAMCH('E'))
C
C     Start the trace messages and convergence history.
C
      IF (VF.NE.0) WRITE (VF,'(I8,2E11.4)') 0, DONE, DONE
      IF (TF.NE.0) WRITE (TF,'(I8,2E11.4)') 0, DONE, DONE
C
C     Initialize the wrapped indices.
C
      ISNM1 = 5
      ISN   = ISNM1
      IVN   = 3
      IVNP1 = IVN
      IWN   = 8
      IWNP1 = IWN
C
C     Set x_0 = 0 and compute the norm of the initial residual.
C
      CALL ZAXPBY (NLEN,VECS(1,IVN),ZONE,VECS(1,2),ZZERO,VECS(1,IVN))
      CALL ZAXPBY (NLEN,VECS(1,1),ZZERO,VECS(1,1),ZZERO,VECS(1,1))
      R0 = DZNRM2(NLEN,VECS(1,IVN),1)
      IF ((TOL.GE.DONE).OR.(R0.EQ.DZERO)) GO TO 80
C
C     Check whether the auxiliary vector must be supplied.
C
      IF (INIT.EQ.0) CALL ZRANDN (NLEN,VECS(1,IWN),1)
C
C     Compute scale factors and check for invariant subspaces.
C
      SCVNP1 = R0
      SCWNP1 = DZNRM2(NLEN,VECS(1,IWN),1)
      IF (SCVNP1.LT.TNRM) IERR = IERR + 16
      IF (SCWNP1.LT.TNRM) IERR = IERR + 32
      IF (IERR.NE.0) GO TO 80
      DNP1 = ZDOTU(NLEN,VECS(1,IWN),1,VECS(1,IVN),1) / ( SCVNP1 * SCWNP1
     $ )
      IF ((SCVNP1.GE.TMAX).OR.(SCVNP1.LE.TMIN)) THEN
         ZTMP = cmplx(DONE / SCVNP1,DZERO,kind=Rkind)
         CALL ZAXPBY (NLEN,VECS(1,IVN),ZTMP,VECS(1,IVN),ZZERO,VECS(1,IVN
     $))
         SCVNP1 = DONE
      END IF
      IF ((SCWNP1.GE.TMAX).OR.(SCWNP1.LE.TMIN)) THEN
         ZTMP = cmplx(DONE / SCWNP1,DZERO,kind=Rkind)
         CALL ZAXPBY (NLEN,VECS(1,IWN),ZTMP,VECS(1,IWN),ZZERO,VECS(1,IWN
     $))
         SCWNP1 = DONE
      END IF
      RHONP1 = SCVNP1
      CSINP1 = SCWNP1
      SCVNP1 = DONE / SCVNP1
      SCWNP1 = DONE / SCWNP1
C
C     Initialize the variables.
C
      N      = 1
      DN     = ZONE
      COSN   = DONE
      RESN   = DONE
      COSNM1 = DONE
      OMGN   = DZERO
      SCSN   = DZERO
      SCVN   = DZERO
      SCWN   = DZERO
      SCSNM1 = DZERO
      SINN   = ZZERO
      SINNM1 = ZZERO
      OMGNP1 = ZUQMXO(N)
      RHSN   = OMGNP1 * R0
      MAXOMG = DONE / OMGNP1
C
C     This is one step of the classical Lanczos algorithm.
C
 20   IVNM1 = IVN
      IVN   = IVNP1
      IVNP1 = MOD(N,2) + 3
      IWNM1 = IWN
      IWN   = IWNP1
      IWNP1 = MOD(N,2) + 8
C
C     Check whether D_n is nonsingular.
C
      DNM1 = DN
      DN   = DNP1
      IF (abs(DN).EQ.DZERO) THEN
         IERR = 8
         GO TO 80
      END IF
C
C     Have the caller carry out AXB, then return here.
C        CALL AXB (VECS(1,IVN),VECS(1,7))
C
      INFO(2) = 1
      INFO(3) = IVN
      INFO(4) = 7
      RETLBL  = 30
      RETURN
 30   RETLBL = 0
C
C     Compute H_{n-1,n} and build part of the vector v_{n+1}.
C
      SCVNM1 = SCVN
      CSIN   = CSINP1
      SCVN   = SCVNP1
      ZTMP   = CSIN * DN / DNM1 * SCVNM1 / SCVN
      CALL ZAXPBY (NLEN,VECS(1,IVNP1),ZONE,VECS(1,7),-ZTMP,VECS(1,IVNM1)
     $)
C
C     Have the caller carry out ATXB, then return here.
C        CALL ATXB (VECS(1,IWN),VECS(1,7))
C
      INFO(2) = 2
      INFO(3) = IWN
      INFO(4) = 7
      RETLBL  = 40
      RETURN
 40   RETLBL = 0
C
C     Build part of the vector w_{n+1}.
C
      SCWNM1 = SCWN
      RHON   = RHONP1
      SCWN   = SCWNP1
      ZTMP   = RHON * DN / DNM1 * SCWNM1 / SCWN
      CALL ZAXPBY (NLEN,VECS(1,IWNP1),ZONE,VECS(1,7),-ZTMP,VECS(1,IWNM1)
     $)
C
C     Compute H_{nn} and finish the new vectors.
C
      RHN = SCVN * SCWN * ZDOTU(NLEN,VECS(1,IWN),1,VECS(1,IVNP1),1) / DN
      CALL ZAXPBY (NLEN,VECS(1,IVNP1),ZONE,VECS(1,IVNP1),-RHN,VECS(1,IVN
     $))
      CALL ZAXPBY (NLEN,VECS(1,IWNP1),ZONE,VECS(1,IWNP1),-RHN,VECS(1,IWN
     $))
C
C     Compute scale factors and check for invariant subspaces.
C
      IERR   = 0
      SCVNP1 = DZNRM2(NLEN,VECS(1,IVNP1),1)
      SCWNP1 = DZNRM2(NLEN,VECS(1,IWNP1),1)
      RHONP1 = SCVN * SCVNP1
      CSINP1 = SCWN * SCWNP1
      RHNP1  = RHONP1
      IF (SCVNP1.LT.TNRM) IERR = IERR + 16
      IF (SCWNP1.LT.TNRM) IERR = IERR + 32
      IF (IERR.NE.0) GO TO 50
      DNP1   = ZDOTU(NLEN,VECS(1,IWNP1),1,VECS(1,IVNP1),1) / ( SCVNP1 *
     $SCWNP1 )
      IF ((SCVNP1.GE.TMAX).OR.(SCVNP1.LE.TMIN)) THEN
         ZTMP = cmplx(DONE / SCVNP1,DZERO,kind=Rkind)
         CALL ZAXPBY (NLEN,VECS(1,IVNP1),ZTMP,VECS(1,IVNP1),ZZERO,VECS(1
     $,IVNP1))
         SCVNP1 = DONE
      END IF
      IF ((SCWNP1.GE.TMAX).OR.(SCWNP1.LE.TMIN)) THEN
         ZTMP = cmplx(DONE / SCWNP1,DZERO,kind=Rkind)
         CALL ZAXPBY (NLEN,VECS(1,IWNP1),ZTMP,VECS(1,IWNP1),ZZERO,VECS(1
     $,IWNP1))
         SCWNP1 = DONE
      END IF
      SCVNP1 = DONE / SCVNP1
      SCWNP1 = DONE / SCWNP1
C
C     The QMR code starts here.
C     Multiply the new column by the previous omega's.
C     Get the next scaling factor omega(i) and update MAXOMG.
C
 50   RHNM1  = CSIN * DN * OMGN / DNM1
      OMGN   = OMGNP1
      RHN    = OMGN * RHN
      OMGNP1 = ZUQMXO(N+1)
      RHNP1  = OMGNP1 * RHNP1
      MAXOMG = DMAX1(MAXOMG,DONE/OMGN)
C
C     Apply the previous rotations.
C
      RHNM2  = SINNM1 * RHNM1
      RHNM1  = COSNM1 * RHNM1
      COSNM1 = COSN
      SINNM1 = SINN
      ZTMP   = RHNM1
      RHNM1  =  COSNM1 * ZTMP + SINNM1 * RHN
      RHN    = -conjg(SINNM1) * ZTMP + COSNM1 * RHN
C
C     Compute the rotation for the last element (this also applies it).
C
      CALL ZROTG (RHN,RHNP1,COSN,SINN)
C
C     Apply the new rotation to the right-hand side vector.
C
      RHSNP1 = -conjg(SINN) * RHSN
      RHSN   =  COSN * RHSN
C
C     Compute the next search direction s_n.
C
      ISNM2  = ISNM1
      ISNM1  = ISN
      ISN    = MOD(N-1,2) + 5
      ZTMP   = SCSNM1 * RHNM2 / SCVN
      CALL ZAXPBY (NLEN,VECS(1,ISN),ZONE,VECS(1,IVN),-ZTMP,VECS(1,ISNM2)
     $)
      SCSNM1 = SCSN
      ZTMP   = SCSNM1 * RHNM1 / SCVN
      CALL ZAXPBY (NLEN,VECS(1,ISN),ZONE,VECS(1,ISN),-ZTMP,VECS(1,ISNM1)
     $)
      SCSN   = SCVN / RHN
C
C     Compute the new QMR iterate, then scale the search direction.
C
      ZTMP = SCSN * RHSN
      CALL ZAXPBY (NLEN,VECS(1,1),ZONE,VECS(1,1),ZTMP,VECS(1,ISN))
      DTMP = abs(SCSN)
      IF ((DTMP.GE.TMAX).OR.(DTMP.LE.TMIN)) THEN
         CALL ZAXPBY (NLEN,VECS(1,ISN),SCSN,VECS(1,ISN),ZZERO,VECS(1,ISN
     $))
         SCSN = ZONE
      END IF
C
C     Compute the residual norm upper bound.
C     If the scaled upper bound is within one order of magnitude of the
C     target convergence norm, compute the true residual norm.
C
      RHSN = RHSNP1
      UNRM = DSQRT(DBLE(N+1)) * MAXOMG * abs(RHSNP1) / R0
      UCHK = UNRM
      IF ((TRES.EQ.0).AND.(UNRM/TOL.GT.DTEN).AND.(N.LT.NLIM)) GO TO 70
C
C     Have the caller carry out AXB, then return here.
C        CALL AXB (VECS(1,1),VECS(1,7))
C
      INFO(2) = 1
      INFO(3) = 1
      INFO(4) = 7
      RETLBL  = 60
      RETURN
 60   RETLBL = 0
      CALL ZAXPBY (NLEN,VECS(1,7),ZONE,VECS(1,2),-ZONE,VECS(1,7))
      RESN = DZNRM2(NLEN,VECS(1,7),1) / R0
      UCHK = RESN
C
C     Output the trace messages and convergence history.
C
 70   IF (VF.NE.0) WRITE (VF,'(I8,2E11.4)') N, UNRM, RESN
      IF (TF.NE.0) WRITE (TF,'(I8,2E11.4)') N, UNRM, RESN
C
C     Check for convergence or termination.  Stop if:
C         1. algorithm converged;
C         2. there is an error condition;
C         3. the residual norm upper bound is smaller than the computed
C     residual norm by a factor of at least 100;
C         4. algorithm exceeded the iterations limit.
C
      IF (RESN.LE.TOL) THEN
         IERR = 0
         GO TO 80
      ELSE IF (IERR.NE.0) THEN
         GO TO 80
      ELSE IF (UNRM.LT.UCHK/DHUN) THEN
         IERR = 4
         GO TO 80
      ELSE IF (N.GE.NLIM) THEN
         IERR = 4
         GO TO 80
      END IF
C
C     Update the running counter.
C
      N = N + 1
      GO TO 20
C
C     That's all.
C
 80   NLIM    = N
      RETLBL  = 0
      INFO(1) = IERR
C
      RETURN
      END
C
C**********************************************************************
C
      DOUBLE PRECISION FUNCTION ZUQMXO (I)
C
C     Purpose:
C     Returns the scaling parameter OMEGA(I).
C
C     Parameters:
C     I = the index of the parameter OMEGA (input).
C
C     Noel M. Nachtigal
C     March 30, 1993
C
C**********************************************************************
C
      INTEGER I
C
      ZUQMXO = 1.0D0
C
      RETURN
      END
C
C**********************************************************************
      subroutine p_multiplyQMR(Vin,Vut,tab_Op,nb_Op,Ene,N,M,eps_in,
     *                         iOp_CAP_Reactif,iOp_CAP_Product)
      use mod_system
      USE mod_Op
      USE mod_CRP
      implicit none

      integer i,N
      real(kind=Rkind) eps_in
      complex(kind=Rkind) Vin(N),Vut(N),b(N)
      complex(kind=Rkind) x(N),y(N), M(N)
      character (len=*) name_sub
      parameter ( name_sub ='p_multiplyQMR' )
      logical precon_flag
!----- Operator variables ----------------------------------------------
      integer           :: nb_Op,iOp_CAP_Reactif,iOp_CAP_Product
      TYPE (param_Op)   :: tab_Op(nb_Op)
      real (kind=Rkind) :: Ene

c     ------------------------------------------
      !M = ONE ! DML
      precon_flag=.TRUE.


c     |b>=e_r*|0>

      b(:)=Vin(:)
      call OpOnVec(b,tab_Op(iOp_CAP_Product),'NOC')

c     |x>=1/(H-E-ie)*|b>

      write(6,*) '# here before QMR 1 '

c      if ( precon_flag ) then
      call qm_precon(N,b,x,M,tab_Op,nb_Op,Ene,eps_in,'CJG')
c      else
c      call qm(N,b,x,tab_Op,eps_in)
c      end if
c      x(:)=b(:)


c     |b>=e_p*|x>


      b(:)=x(:)
      call OpOnVec(b,tab_Op(iOp_CAP_Reactif),'NOC')



c     |x>=1/(H-E+ie)*|b>

      write(6,*) '# here before QMR 2 '

c      if ( precon_flag ) then
      call qm_precon(N,b,x,M,tab_Op,nb_Op,Ene,eps_in,'NOC')

c      else
c         call qm(N,b,x,tab_Op,eps_in)
c      end if

      Vut(:)=x(:)

      End
      subroutine qm(N,x,y,eps_in)
      use mod_system
      implicit NONE
c
!      include 'headers/sizes.h'
c
      integer i,N,nlim,info(4),revcom,colx,colb
      real*8  eps_in, tol
      complex*16 vecs(N,7),x(N)
      complex*16 y(N)
c
c     Input with random start vector
c
c      call zrandn(N,y,-1)       ! fill y with random numbers
c      call g_inv_multiply(N,y,vecs(1,2))    ! vecs(:,2)=A*y
c
c     x is given right-hand side (= b)
      vecs(:,2)=(x(:)-vecs(:,2)) ! vecs(:,2)=M1*(x-A*y)

c     -------------------------------------------------------
      nlim=10*N
      tol=eps_in
      info(1)=10000
      info(2)=0
      info(3)=0
      info(4)=0
c
 10   call zscpx(N,N,nlim,vecs,tol,info)
c
      revcom = info(2)
      colx   = info(3)
      colb   = info(4)
c
      if(revcom == 1) then        !Multiply y=H*x
c       Right preconditioner
        x(:)=vecs(:,colx)

c       Matrix-vector mult.
c       call g_inv_multiply(N,vecs(1,colx),vecs(1,colb))
c
c       Left preconditioner
        vecs(:,colb)=vecs(:,colb)
         go to 10
      end if
c
      if(info(1) == 0) then
         write(6,*)'Iter',nlim
      else if(info(i).ne.0) then
         go to 999
      end if
c     -------------------------------------------------
c     Output
      y(:) = y(:)+vecs(:,1)

      return
c
c     -------------------- E R R O R ----------------------
 999  continue
      if(info(1).eq.1) write(6,*)'Invalid reverse comm. call'
      if(info(1).eq.2) write(6,*)'Invalid inputs'
      if(info(1).eq.4) write(6,*)'No convergence after Nlim steps',nlim
      if(info(1).eq.8) write(6,*)'Exact B R E A K D O W N'
      if(info(1).eq.16) write(6,*) 'A invariant subspace found'
      if(info(1).eq.32) write(6,*) 'A^T invariant subspace found'
      if(info(1).eq.48) write(6,*) ' Both subspaces found'
      stop
c     -------------------- E R R O R ----------------------
c
      end


      subroutine qm_precon(N,x,y,M,tab_Op,nb_Op,Ene,eps_in,cjg)
      use mod_system
      USE mod_Op
      USE mod_CRP
      implicit NONE

      integer i,N,nlim,info(4),revcom,colx,colb,it
      real(kind=Rkind)  eps_in,tol
      complex(kind=Rkind) vecs(N,7),x(N),x2(N),y(N)
      character (len=3) cjg
!----- Operator variables ----------------------------------------------
      integer           :: nb_Op
      TYPE (param_Op)   :: tab_Op(nb_Op)


! Matrices de préconditionnement ( diagonale donc stockable comme un vecteur)
      complex(kind=Rkind)  :: M(N)

      real (kind=Rkind) :: Ene


c     x is given right-hand side (= b)
      vecs(:,2)=M(:)*(x(:))  ! vecs(:,2)=M*(x-A*y)


      it = 0

c     -------------------------------------------------------
      nlim=10*N
      tol=eps_in
      info(1)=00000
      info(2)=0
      info(3)=0
      info(4)=0
      write(out_unitp,'(a)',advance='no') 'QMR it:'

c
c      write(*,*) 'nlim = ', nlim, ' info(1) = ', info(1), ' N =', N
c      write(*,*) 'tol = ', tol, ' maxDVRpoints = ', maxDVRpoints
 10   call zscpx(N,N,nlim,vecs,tol,info)
c      write(*,*) 'nlim = ', nlim, ' info(1) = ', info(1), ' N =',N
c      write(*,*) 'tol = ', tol, ' maxDVRpoints = ', maxDVRpoints
       it = it + 1
       !write(*,*) 'QMR it',it
       IF (mod(it,10) == 0 .AND. mod(it,100) /= 0) THEN
          write(out_unitp,'(a)',advance='no') '.'
       END IF
       IF (mod(it,100) == 0) write(out_unitp,'(i0)',advance='no') it


      revcom = info(2)
      colx   = info(3)
      colb   = info(4)
c
      if(revcom == 1) then        !Multiply y=H*x
c       Right preconditioner
        x(:)=M(:)*vecs(:,colx)

c       Matrix-vector mult.
        call Gpsi(x,tab_Op,nb_Op,Ene,cjg)

        vecs(:,colb)=x(:)

c       Left preconditioner
        vecs(:,colb)=M(:)*vecs(:,colb)
        go to 10
      end if
c
      if(info(1) == 0) then
         write(6,*)'# Iter',nlim
      else if(info(1).ne.0) then
         go to 999
      end if
c     -------------------------------------------------
c     Output
      y(:)=M(:)*vecs(:,1)
      return
c
c     -------------------- E R R O R ----------------------
 999  continue
      if(info(1).eq.1) write(6,*)'Invalid reverse comm. call'
      if(info(1).eq.2) write(6,*)'Invalid inputs'
      if(info(1).eq.4) write(6,*)'No convergence after Nlim steps',nlim
      if(info(1).eq.8) write(6,*)'Exact B R E A K D O W N'
      if(info(1).eq.16) write(6,*) 'A invariant subspace found'
      if(info(1).eq.32) write(6,*) 'A^T invariant subspace found'
      if(info(1).eq.48) write(6,*) ' Both subspaces found'
      stop
c     -------------------- E R R O R ----------------------
c
      end
