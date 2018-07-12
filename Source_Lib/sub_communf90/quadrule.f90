!  quadrule.f  25 September 1998
!
      subroutine bashset ( norder, xtab, weight )
!
!***********************************************************************
!
!! BASHSET sets abscissas and weights for Adams-Bashforth quadrature.
!
!
!  Definition:
!
!    Adams-Bashforth quadrature formulas are normally used in solving
!    ordinary differential equations, and are not really suitable for
!    general quadrature computations.  However, an Adams-Bashforth formula
!    is equivalent to approximating the integral of F(Y(X)) between X(M)
!    and X(M+1), using an explicit formula that relies only on known values
!    of F(Y(X)) at X(M-N+1) through X(M).  For this reason, the formulas
!    have been included here.
!
!    Suppose the unknown function is denoted by Y(X), with derivative
!    F(Y(X)), and that approximate values of the function are known at a
!    series of X values, which we write as X(1), X(2), ..., X(M).  We write
!    the value Y(X(1)) as Y(1) and so on.
!
!    Then the solution of the ODE Y'=F(X,Y) at the next point X(M+1) is
!    computed by:
!
!      Y(M+1) = Y(M) + Integral ( X(M) < X < X(M+1) ) F(Y(X)) dX
!             = Y(M) + H * Sum ( I = 1 to N ) W(I) * F(Y(M+1-I)) approximately.
!
!    In the documentation that follows, we replace F(Y(X)) by F(X).
!
!  Integration interval:
!
!    [ 0, 1 ].
!
!  Weight function:
!
!    1.
!
!  Integral to approximate:
!
!    INTEGRAL ( 0 <= X <= 1 ) F(X) dX.
!
!  Approximate integral:
!
!    SUM ( I = 1 to NORDER ) WEIGHT(I) * F ( 1 - I ),
!
!  Note:
!
!    The Adams-Bashforth formulas require equally spaced data.
!
!    Here is how the formula is applied in the general case with non-unit spacing:
!
!      INTEGRAL ( A <= X <= A+H ) F(X) dX =
!      H * SUM ( I = 1 to NORDER ) WEIGHT(I) * F ( A - (I-1)*H ), approximately.
!
!    The reference lists the second coefficient of the order 8 Adams-Bashforth
!    formula as
!      weight(2) =  -1162169.0 / 120960.0
!    but this should be
!      weight(2) =  -1152169.0 / 120960.0
!
!  Reference:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964.
!
!    Jean Lapidus and John Seinfeld,
!    Numerical Solution of Ordinary Differential Equations,
!    Academic Press, 1971.
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Modified:
!
!    15 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.  NORDER should be
!    between 1 and 8.
!
!    Output, real(kind=Rkind) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real(kind=Rkind) WEIGHT(NORDER), the weights of the rule.
!    WEIGHT(1) is the weight at X = 0, WEIGHT(2) the weight at X = -1,
!    and so on.  The weights are rational, and should sum to 1.  Some
!    weights may be negative.
!
      integer norder
!
      real(kind=Rkind) d
      integer i
      real(kind=Rkind) weight(norder)
      real(kind=Rkind) xtab(norder)
!
      if ( norder .eq. 1 ) then

        weight(1) = 1.0 D+00

      else if ( norder .eq. 2 ) then

        d = 2.0 D+00

        weight(1) =   3.0 D+00 / d
        weight(2) = - 1.0 D+00 / d

      else if ( norder .eq. 3 ) then

        d = 12.0 D+00

        weight(1) =   23.0 D+00 / d
        weight(2) = - 16.0 D+00 / d
        weight(3) =    5.0 D+00 / d

      else if ( norder .eq. 4 ) then

        d = 24.0 D+00

        weight(1) =   55.0 D+00 / d
        weight(2) = - 59.0 D+00 / d
        weight(3) =   37.0 D+00 / d
        weight(4) =  - 9.0 D+00 / d

      else if ( norder .eq. 5 ) then

        d = 720.0 D+00

        weight(1) =   1901.0 D+00 / d
        weight(2) = - 2774.0 D+00 / d
        weight(3) =   2616.0 D+00 / d
        weight(4) = - 1274.0 D+00 / d
        weight(5) =    251.0 D+00 / d

      else if ( norder .eq. 6 ) then

        d = 1440.0 D+00

        weight(1) =   4277.0 D+00 / d
        weight(2) = - 7923.0 D+00 / d
        weight(3) =   9982.0 D+00 / d
        weight(4) = - 7298.0 D+00 / d
        weight(5) =   2877.0 D+00 / d
        weight(6) =  - 475.0 D+00 / d

      else if ( norder .eq. 7 ) then

        d = 60480.0 D+00

        weight(1) =    198721.0 D+00 / d
        weight(2) =  - 447288.0 D+00 / d
        weight(3) =    705549.0 D+00 / d
        weight(4) =  - 688256.0 D+00 / d
        weight(5) =    407139.0 D+00 / d
        weight(6) =  - 134472.0 D+00 / d
        weight(7) =     19087.0 D+00 / d

      else if ( norder .eq. 8 ) then

        d = 120960.0 D+00

        weight(1) =     434241.0 D+00 / d
        weight(2) =  - 1152169.0 D+00 / d
        weight(3) =    2183877.0 D+00 / d
        weight(4) =  - 2664477.0 D+00 / d
        weight(5) =    2102243.0 D+00 / d
        weight(6) =  - 1041723.0 D+00 / d
        weight(7) =     295767.0 D+00 / d
        weight(8) =    - 36799.0 D+00 / d

      else

        write ( *, * ) ' '
        write ( *, * ) 'BASHSET - Fatal error!'
        write ( *, * ) '  Illegal value of NORDER = ', norder
        write ( *, * ) '  Legal values are 1 through 8.'
        stop

      end if

      do i = 1, norder
        xtab(i) = dble ( 1 - i )
      end do

      return
      end subroutine bashset
      subroutine bdfcset ( norder, weight, xtab )
!
!***********************************************************************
!
!! BDFCSET sets weights for backward differentiation corrector quadrature.
!
!
!  Definition:
!
!    A backward differentiation corrector formula is defined for a set
!    of evenly spaced abscissas X(I) with X(1) = 1 and X(2) = 0.  Assuming
!    that the values of the function to be integrated are known at the
!    abscissas, the formula is written in terms of the function value at
!    X(1), and the backward differences at X(1) that approximate the
!    derivatives there.  A backward differentiation corrector formula
!    is equivalent to an Adams-Moulton formula of the same order.
!
!  Integration interval:
!
!    [ 0, 1 ].
!
!  Weight function:
!
!    1.
!
!  Integral to approximate:
!
!    INTEGRAL ( 0 <= X <= 1 ) F(X) dX.
!
!  Approximate integral:
!
!    SUM ( I = 1 to NORDER ) WEIGHT(I) * BD**(I-1) F ( 1 ).
!
!    Here, "BD**(I-1) F ( 0 )" denotes the (I-1)st backward difference
!    of F at X = 0, using a spacing of 1.  In particular,
!
!    BD**0 F(1) = F(1)
!    BD**1 F(1) = F(1) - F(0)
!    BD**2 F(1) = F(1) - 2 * F(0) + F(-1 )
!
!  Reference:
!
!    Simeon Fatunla,
!    Numerical Methods for Initial Value Problems in Ordinary Differential
!      Equations,
!    Academic Press, 1988.
!
!  Modified:
!
!    15 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule, which can be
!    any value from 1 to 19.
!
!    Output, real(kind=Rkind) WEIGHT(NORDER), the weights of the rule.
!
!    Output, real(kind=Rkind) XTAB(NORDER), the abscissas of the rule.
!
      integer maxord
      parameter ( maxord = 19 )
!
      integer norder
!
      integer i
      real(kind=Rkind) w(maxord)
      real(kind=Rkind) weight(norder)
      real(kind=Rkind) xtab(norder)
!
      w(1) =                 1.0 D+00
      w(2) =               - 1.0 D+00 /               2.0 D+00
      w(3) =               - 1.0 D+00 /              12.0 D+00
      w(4) =               - 1.0 D+00 /              24.0 D+00
      w(5) =              - 19.0 D+00 /             720.0 D+00
      w(6) =               - 3.0 D+00 /             160.0 D+00
      w(7) =             - 863.0 D+00 /           60480.0 D+00
      w(8) =             - 275.0 D+00 /           24792.0 D+00
      w(9) =           - 33953.0 D+00 /         3628800.0 D+00
      w(10) =           - 8183.0 D+00 /         1036800.0 D+00
      w(11) =        - 3250433.0 D+00 /       479001600.0 D+00
      w(12) =           - 4671.0 D+00 /          788480.0 D+00
      w(13) =    - 13695779093.0 D+00 /   2615348736000.0 D+00
      w(14) =     - 2224234463.0 D+00 /    475517952000.0 D+00
      w(15) =   - 132282840127.0 D+00 /  31384184832000.0 D+00
      w(16) =     - 2639651053.0 D+00 /    689762304000.0 D+00
      w(17) =  111956703448001.0 D+00 /   3201186852864.0 D+00
      w(18) =         50188465.0 D+00 /     15613165568.0 D+00
      w(19) = 2334028946344463.0 D+00 / 786014494949376.0 D+00

      do i = 1, min ( norder, maxord )
        weight(i) = w(i)
      end do

      do i = 1, norder
        xtab(i) = dble ( 1 - i )
      end do

      return
      end subroutine bdfcset
      subroutine bdfpset ( norder, weight, xtab )
!
!***********************************************************************
!
!! BDFPSET sets weights for backward differentiation predictor quadrature.
!
!
!  Definition:
!
!    A backward differentiation predictor formula is defined for a set
!    of evenly spaced abscissas X(I) with X(1) = 1 and X(2) = 0.  Assuming
!    that the values of the function to be integrated are known at the
!    abscissas, the formula is written in terms of the function value at
!    X(2), and the backward differences at X(2) that approximate the
!    derivatives there.  A backward differentiation predictor formula
!    is equivalent to an Adams-Bashforth formula of the same order.
!
!  Integration interval:
!
!    [ 0, 1 ].
!
!  Weight function:
!
!    1.
!
!  Integral to approximate:
!
!    INTEGRAL ( 0 <= X <= 1 ) F(X) dX.
!
!  Approximate integral:
!
!    SUM ( I = 1 to NORDER ) WEIGHT(I) * BD**(I-1) F ( 0 ),
!
!    Here, "BD**(I-1) F ( 0 )" denotes the (I-1)st backward difference
!    of F at X = 0, using a spacing of 1.  In particular,
!
!    BD**0 F(0) = F(0)
!    BD**1 F(0) = F(0) - F(-1)
!    BD**2 F(0) = F(0) - 2 * F(-1) + F(-2 )
!
!  Note:
!
!    The relationship between a backward difference predictor and the
!    corresponding Adams-Bashforth formula may be illustrated for the
!    BDF predictor of order 3:
!
!      BD**0 F(0) + 0.5 * BD**1 F(0) + 5/12 * BD**2 F(0)
!      = F(0) + 0.5 * ( F(0) - F(1) + 5/12 * ( F(0) - 2 * F(-1) + F(-2) )
!      = 23/12 F(0) - 16/12 F(-1) + 5/12 F(-2)
!
!    which is the Adams-Bashforth formula of order 3.
!
!  Reference:
!
!    Simeon Fatunla,
!    Numerical Methods for Initial Value Problems in Ordinary Differential
!      Equations,
!    Academic Press, 1988.
!
!  Modified:
!
!    15 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule, which can be
!    any value from 1 to 19.
!
!    Output, real(kind=Rkind) WEIGHT(NORDER), the weight of the rule.
!
!    Output, real(kind=Rkind) XTAB(NORDER), the abscissas of the rule.
!
      integer maxord
      parameter ( maxord = 19 )
!
      integer norder
!
      integer i
      real(kind=Rkind) w(maxord)
      real(kind=Rkind) weight(norder)
      real(kind=Rkind) xtab(norder)
!
      w(1) =                       1.0 D+00
      w(2) =                       1.0 D+00 /                2.0 D+00
      w(3) =                       5.0 D+00 /               12.0 D+00
      w(4) =                       3.0 D+00 /                8.0 D+00
      w(5) =                     251.0 D+00 /              720.0 D+00
      w(6) =                      95.0 D+00 /              288.0 D+00
      w(7) =                   19087.0 D+00 /            60480.0 D+00
      w(8) =                    5257.0 D+00 /            17280.0 D+00
      w(9) =                 1070017.0 D+00 /          3628800.0 D+00
      w(10) =                  25713.0 D+00 /            89600.0 D+00
      w(11) =               26842253.0 D+00 /         95800320.0 D+00
      w(12) =                4777223.0 D+00 /         17418240.0 D+00
      w(13) =           703604254357.0 D+00 /    2615348736000.0 D+00
      w(14) =           106364763817.0 D+00 /     402361344000.0 D+00
      w(15) =          1166309819657.0 D+00 /    4483454976000.0 D+00
      w(16) =               25221445.0 D+00 /         98402304.0 D+00
      w(17) =       8092989203533249.0 D+00 /    3201186852864.0 D+00
      w(18) =         85455477715379.0 D+00 /      34237292544.0 D+00
      w(19) =   12600467236042756559.0 D+00 / 5109094217170944.0 D+00

      do i = 1, min ( norder, maxord )
        weight(i) = w(i)
      end do

      do i = 1, norder
        xtab(i) = dble ( 2 - i )
      end do

      return
      end subroutine bdfpset
      subroutine bdfsum ( func, norder, weight, xtab, diftab, result )
!
!***********************************************************************
!
!! BDFSUM carries out an explicit backward difference quadrature rule for [0,1].
!
!
!  Integral to approximate:
!
!    INTEGRAL ( 0 <= X <= 1 ) F(X) dX.
!
!  Formula:
!
!    RESULT = SUM ( I = 1 to NORDER ) WEIGHT(I) * BDF**(I-1) FUNC ( 0 ).
!
!  Note:
!
!    The integral from 0 to 1 is approximated using data at X = 0,
!    -1, -2, ..., -NORDER+1.  This is a form of extrapolation, and
!    the approximation can become poor as NORDER increases.
!
!  Modified:
!
!    15 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FUNC, the name of the FORTRAN function which evaluates
!    the integrand.  The function must have the form
!      real(kind=Rkind) func ( x ).
!
!    Input, integer NORDER, the order of the rule.
!
!    Input, real(kind=Rkind) WEIGHT(NORDER), the weights of the rule.
!
!    Input, real(kind=Rkind) XTAB(NORDER), the abscissas of the rule.
!
!    Workspace, real(kind=Rkind) DIFTAB(NORDER).
!
!    Output, real(kind=Rkind) RESULT, the approximate value of the integral.
!
      integer norder
!
      real(kind=Rkind) diftab(norder)
      real(kind=Rkind) func
      integer i
      integer j
      real(kind=Rkind) result
      real(kind=Rkind) weight(norder)
      real(kind=Rkind) xtab(norder)
!
      external func
!
      do i = 1, norder
        diftab(i) = func ( xtab(i) )
      end do

      do i = 2, norder
        do j = i, norder
          diftab(norder+i-j) =                                          &
            ( diftab(norder+i-j-1) - diftab(norder+i-j) )
        end do
      end do

      result = 0.0 D+00
      do i = 1, norder
        result = result + weight(i) * diftab(i)
      end do

      return
      end subroutine bdfsum
      subroutine chebset ( norder, xtab, weight )
!
!***********************************************************************
!
!! CHEBSET sets abscissas and weights for Chebyshev quadrature.
!
!
!  Integration interval:
!
!    [ -1, 1 ].
!
!  Weight function:
!
!    1.
!
!  Integral to approximate:
!
!    INTEGRAL ( -1 <= X <= 1 ) F(X) dX.
!
!  Approximate integral:
!
!    SUM ( I = 1 to NORDER ) WEIGHT(I) * F ( XTAB(I) ).
!
!  Note:
!
!    The Chebyshev rule is distinguished by using equal weights.
!
!  Reference:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964.
!
!    Zdenek Kopal,
!    Numerical Analysis,
!    John Wiley, 1955.
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Modified:
!
!    15 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!    NORDER may only have the values 1, 2, 3, 4, 5, 6, 7 or 9.
!    There are NO other Chebyshev rules with real abscissas.
!
!    Output, real(kind=Rkind) XTAB(NORDER), the abscissas of the rule,
!    which are symmetric in [-1,1].
!
!    Output, real(kind=Rkind) WEIGHT(NORDER), the weights of the rule,
!    which should each equal 2 / NORDER.
!
      integer norder
!
      integer i
      real(kind=Rkind) weight(norder)
      real(kind=Rkind) xtab(norder)
!
      if ( norder .eq. 1 ) then

        xtab(1) = 0.0 D+00

      else if ( norder .eq. 2 ) then

        xtab(1) = - 1.0 D+00 / sqrt ( 3.0 D+00 )
        xtab(2) =   1.0 D+00 / sqrt ( 3.0 D+00 )

      else if ( norder .eq. 3 ) then

        xtab(1) = - 1.0 D+00 / sqrt ( 2.0 D+00 )
        xtab(2) =   0.0 D+00
        xtab(3) =   1.0 D+00 / sqrt ( 2.0 D+00 )

      else if ( norder .eq. 4 ) then

        xtab(1) = - 0.7946544723 D+00
        xtab(2) = - 0.1875924741 D+00
        xtab(3) =   0.1875924741 D+00
        xtab(4) =   0.7946544723 D+00

      else if ( norder .eq. 5 ) then

        xtab(1) = - 0.8324974870 D+00
        xtab(2) = - 0.3745414096 D+00
        xtab(3) =   0.0 D+00
        xtab(4) =   0.3745414096 D+00
        xtab(5) =   0.8324974870 D+00

      else if ( norder .eq. 6 ) then

        xtab(1) = - 0.8662468181 D+00
        xtab(2) = - 0.4225186538 D+00
        xtab(3) = - 0.2666354015 D+00
        xtab(4) =   0.2666354015 D+00
        xtab(5) =   0.4225186538 D+00
        xtab(6) =   0.8662468181 D+00

      else if ( norder .eq. 7 ) then

        xtab(1) = - 0.8838617008 D+00
        xtab(2) = - 0.5296567753 D+00
        xtab(3) = - 0.3239118105 D+00
        xtab(4) =   0.0 D+00
        xtab(5) =   0.3239118105 D+00
        xtab(6) =   0.5296567753 D+00
        xtab(7) =   0.8838617008 D+00

      else if ( norder .eq. 9 ) then

        xtab(1) = - 0.9115893077 D+00
        xtab(2) = - 0.6010186554 D+00
        xtab(3) = - 0.5287617831 D+00
        xtab(4) = - 0.1679061842 D+00
        xtab(5) =   0.0 D+00
        xtab(6) =   0.1679061842 D+00
        xtab(7) =   0.5287617831 D+00
        xtab(8) =   0.6010186554 D+00
        xtab(9) =   0.9115893077 D+00

      else

        write ( *, * ) ' '
        write ( *, * ) 'CHEBSET - Fatal error!'
        write ( *, * ) '  Illegal value of NORDER = ', norder
        write ( *, * ) '  Legal values are 1 through 7, and 9.'
        stop

      end if

      do i = 1, norder
        weight(i) = 2.0 D+00 / dble ( norder )
      end do

      return
      end subroutine chebset
      subroutine chebtoset ( norder, xtab, weight )
!
!***********************************************************************
!
!! CHEBTOSET sets up open Gauss-Chebyshev (first kind) quadrature.
!
!
!  Integration interval:
!
!    [ -1, 1 ]
!
!  Weight function:
!
!    1 / SQRT ( 1 - X**2 ).
!
!  Integral to approximate:
!
!    INTEGRAL ( -1 <= X <= 1 ) F(X) / SQRT ( 1 - X**2 ) dX.
!
!  Approximate integral:
!
!    SUM ( I = 1 to NORDER ) WEIGHT(I) * F ( XTAB(I) ).
!
!  Precision:
!
!    If NORDER points are used, then Gauss-Chebyshev quadrature
!    will compute the integral exactly, whenever F(X) is a polynomial
!    of degree 2*NORDER-1 or less.
!
!  Note:
!
!    The abscissas of the rule are zeroes of the Chebyshev polynomials
!    of the first kind, T(NORDER)(X).
!
!  Reference:
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Modified:
!
!    15 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!
!    Output, real(kind=Rkind) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real(kind=Rkind) WEIGHT(NORDER), the weights of the rule,
!    which are all equal to PI / NORDER.
!
      real(kind=Rkind) pi
      parameter ( pi = 3.141592653589793238462643383279D+00 )
!
      integer norder
!
      real(kind=Rkind) angle
      integer i
      real(kind=Rkind) weight(norder)
      real(kind=Rkind) xtab(norder)
!
      do i = 1, norder
        angle = dble ( 2 * i - 1 ) * pi / dble ( 2 * norder )
        xtab(i) = cos ( angle )
        weight(i) = pi / dble ( norder )
      end do

      return
      end subroutine chebtoset
      subroutine chebtcset ( norder, xtab, weight )
!
!***********************************************************************
!
!! CHEBTCSET sets up closed Gauss-Chebyshev (first kind) quadrature.
!
!
!  Integration interval:
!
!    [ -1, 1 ]
!
!  Weight function:
!
!    1 / SQRT ( 1 - X**2 ).
!
!  Integral to approximate:
!
!    INTEGRAL ( -1 <= X <= 1 ) F(X) / SQRT ( 1 - X**2 ) dX.
!
!  Approximate integral:
!
!    SUM ( I = 1 to NORDER ) WEIGHT(I) * F ( XTAB(I) ).
!
!  Precision:
!
!    If NORDER points are used, then Gauss-Chebyshev quadrature
!    will compute the integral exactly, whenever F(X) is a polynomial
!    of degree 2*NORDER-3 or less.
!
!  Note:
!
!    The abscissas include -1 and 1.
!
!    If the order is doubled, the abscissas of the new rule include
!    all the points of the old rule.  This fact can be used to
!    efficiently implement error estimation.
!
!  Reference:
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Modified:
!
!    15 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule, which must be
!    at least 2.
!
!    Output, real(kind=Rkind) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real(kind=Rkind) WEIGHT(NORDER), the weights of the rule.
!    The first and last weights are 0.5 * PI / ( NORDER - 1),
!    and all other weights are PI / ( NORDER - 1 ).
!
      real(kind=Rkind) pi
      parameter ( pi = 3.141592653589793238462643383279D+00 )
!
      integer norder
!
      real(kind=Rkind) angle
      integer i
      real(kind=Rkind) weight(norder)
      real(kind=Rkind) xtab(norder)
!
      if ( norder .lt. 2 ) then
        write ( *, * ) ' '
        write ( *, * ) 'CHEBTCSET - Fatal error!'
        write ( *, * ) '  NORDER must be at least 2.'
        write ( *, * ) '  The input value was NORDER = ', norder
        stop
      end if

      do i = 1, norder

        angle = dble ( i - 1 ) * pi / dble ( norder - 1 )
        xtab(i) = cos ( angle )

        if ( i .eq. 1 .or. i .eq. norder ) then
          weight(i) = pi / dble ( 2 * ( norder - 1 ) )
        else
          weight(i) = pi / dble ( norder - 1 )
        end if

      end do

      return
      end subroutine chebtcset
      subroutine chebuset ( norder, xtab, weight )
!
!***********************************************************************
!
!! CHEBUSET sets abscissas and weights for Gauss-Chebyshev quadrature.
!
!
!  Integration interval:
!
!    [ -1, 1 ]
!
!  Weight function:
!
!    SQRT ( 1 - X**2 ).
!
!  Integral to approximate:
!
!    INTEGRAL ( -1 <= X <= 1 ) SQRT ( 1 - X**2 ) * F(X) dX.
!
!  Approximate integral:
!
!    SUM ( I = 1 to NORDER ) WEIGHT(I) * F ( XTAB(I) ).
!
!  Precision:
!
!    If NORDER points are used, then Gauss-Chebyshev quadrature
!    will compute the integral exactly, whenever F(X) is a polynomial
!    of degree 2*NORDER-1 or less.
!
!  Note:
!
!    The abscissas are zeroes of the Chebyshev polynomials
!    of the second kind, U(NORDER)(X).
!
!  Reference:
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Modified:
!
!    15 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!
!    Output, real(kind=Rkind) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real(kind=Rkind) WEIGHT(NORDER), the weights of the rule,
!    which are all equal to PI / NORDER.
!
      real(kind=Rkind) pi
      parameter ( pi = 3.141592653589793238462643383279D+00 )
!
      integer norder
!
      real(kind=Rkind) angle
      integer i
      real(kind=Rkind) weight(norder)
      real(kind=Rkind) xtab(norder)
!
      do i = 1, norder
        angle = dble ( i ) * pi / dble ( norder + 1 )
        xtab(i) = cos ( angle )
        weight(i) = pi * ( sin ( angle ) )**2 / dble ( norder + 1 )
      end do

      return
      end subroutine chebuset
      function gamma ( x )
!
!***********************************************************************
!
!! GAMMA computes the gamma function using Hastings's approximation.
!
!
!  Reference:
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!  Modified:
!
!    19 September 1998
!
!  Parameters:
!
!    Input, real(kind=Rkind) X, the argument at which the gamma function
!    is to be evaluated.  X must be greater than 0, and less than 70.
!
!    Output, real(kind=Rkind) GAMMA, the gamma function at X.
!
      real(kind=Rkind) gam
      real(kind=Rkind) gamma
      real(kind=Rkind) x
      real(kind=Rkind) y
      real(kind=Rkind) z
      real(kind=Rkind) za
!
      gam ( y ) = (((((((                                               &
          0.035868343 D+00   * y                                        &
        - 0.193527818 D+00 ) * y                                        &
        + 0.482199394 D+00 ) * y                                        &
        - 0.756704078 D+00 ) * y                                        &
        + 0.918206857 D+00 ) * y                                        &
        - 0.897056937 D+00 ) * y                                        &
        + 0.988205891 D+00 ) * y                                        &
        - 0.577191652 D+00 ) * y + 1.0 D+00
!
      if ( x .le. 0.0 ) then
        gamma = 0.0
        write ( *, * ) ' '
        write ( *, * ) 'GAMMA - Fatal error!'
        write ( *, * ) '  Input argument X <= 0.'
        stop
      end if

      if ( x .ge. 70.0 ) then
        write ( *, * ) ' '
        write ( *, * ) 'GAMMA - Fatal error!'
        write ( *, * ) '  Input argument X >= 70.'
        stop
      end if

      if ( x .eq. 1.0 ) then
        gamma = 1.0 D+00
        return
      end if

      if ( x .le. 1.0 ) then
        gamma = gam ( x ) / x
        return
      end if

      z = x

      za = 1.0

   10 continue

      z = z - 1.0

      if ( z .lt. 1.0 ) then
        gamma = za * gam ( z )
        return
      else if ( z .eq. 1.0 ) then
        gamma = za
        return
      else if ( z .gt. 1.0 ) then
        za = za * z
        go to 10
      end if

      return
      end function gamma
      subroutine hercom ( norder, xtab, weight )
!
!***********************************************************************
!
!! HERCOM computes the abscissa and weights for Gauss-Hermite quadrature.
!
!
!  The abscissas are the zeros of the N-th order Hermite polynomial.
!
!  Integration interval:
!
!    ( -Infinity, +Infinity )
!
!  Weight function:
!
!    EXP ( - X**2 ).
!
!  Integral to approximate:
!
!    INTEGRAL ( -INFINITY < X < +INFINITY ) EXP ( - X**2 ) * F(X) dX.
!
!  Approximate integral:
!
!    SUM ( I = 1 to NORDER ) WEIGHT(I) * F ( XTAB(I) ).
!
!  Reference:
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!  Modified:
!
!    19 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NORDER, the order of the formula to be computed.
!
!    Output, real(kind=Rkind) XTAB(NORDER), the Gauss-Hermite abscissas.
!
!    Output, real(kind=Rkind) WEIGHT(NORDER), the Gauss-Hermite weights.
!
      integer norder
!
      real(kind=Rkind) cc
      real(kind=Rkind) dp2
      real(kind=Rkind) gamma
      integer i
      real(kind=Rkind) p1
      real(kind=Rkind) s
      real(kind=Rkind) temp
      real(kind=Rkind) weight(norder)
      real(kind=Rkind) x
      real(kind=Rkind) xtab(norder)
!
      cc = 1.7724538509 * gamma ( dble ( norder ) )                     &
        / ( 2.0**(norder-1) )

      s = ( 2.0 * dble ( norder ) + 1.0 )**( 1.0 / 6.0 )

      do i = 1, ( norder + 1 ) / 2

        if ( i .eq. 1 ) then

          x = s**3 - 1.85575 / s

        else if ( i .eq. 2 ) then

          x = x - 1.14 * ( ( dble ( norder ) )**0.426 ) / x

        else if ( i .eq. 3 ) then

          x = 1.86 * x - 0.86 * xtab(1)

        else if ( i .eq. 4 ) then

          x = 1.91 * x - 0.91 * xtab(2)

        else

          x = 2.0 * x - xtab(i-2)

        end if

        call herroot ( x, norder, dp2, p1 )

        xtab(i) = x
        weight(i) = cc / dp2 / p1

        xtab(norder-i+1) = - x
        weight(norder-i+1) = weight(i)

      end do
!
!  Reverse the order of the XTAB values.
!
      do i = 1, norder/2
        temp = xtab(i)
        xtab(i) = xtab(norder+1-i)
        xtab(norder+1-i) = temp
      end do

      return
      end subroutine hercom
      subroutine herrec ( p2, dp2, p1, x, norder )
!
!***********************************************************************
!
!! HERREC finds the value and derivative of the NORDER-th Hermite polynomial.
!
!
!  Reference:
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!  Modified:
!
!    19 September 1998
!
!  Parameters:
!
!    Output, real(kind=Rkind) P2, the value of H(NORDER)(X).
!
!    Output, real(kind=Rkind) DP2, the value of H'(NORDER)(X).
!
!    Output, real(kind=Rkind) P1, the value of H(NORDER-1)(X).
!
!    Input, real(kind=Rkind) X, the point at which polynomials are evaluated.
!
!    Input, integer NORDER, the order of the polynomial to be computed.
!
      integer i
      real(kind=Rkind) dp0
      real(kind=Rkind) dp1
      real(kind=Rkind) dp2
      integer norder
      real(kind=Rkind) p0
      real(kind=Rkind) p1
      real(kind=Rkind) p2
      real(kind=Rkind) x
!
      p1 = 1.0
      dp1 = 0.0

      p2 = x
      dp2 = 1.0

      do i = 2, norder

        p0 = p1
        dp0 = dp1

        p1 = p2
        dp1 = dp2

        p2  = x * p1 - 0.5 * ( dble ( i ) - 1.0 ) * p0
        dp2 = x * dp1 + p1 - 0.5 * ( dble ( i ) - 1.0 ) * dp0

      end do

      return
      end subroutine herrec
      subroutine herroot ( x, norder, dp2, p1 )
!
!***********************************************************************
!
!! HERROOT improves an approximate root of a Hermite polynomial.
!
!
!  Reference:
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!  Modified:
!
!    19 September 1998
!
!  Parameters:
!
!    Input/output, real(kind=Rkind) X, the approximate root, which
!    should be improved on output.
!
!    Input, integer NORDER, the order of the Hermite polynomial.
!
!    Output, real(kind=Rkind) DP2, the value of H'(NORDER)(X).
!
!    Output, real(kind=Rkind) P1, the value of H(NORDER-1)(X).
!
      real(kind=Rkind) eps
      parameter ( eps = 1.0D-12 )
!
      real(kind=Rkind) d
      real(kind=Rkind) dp2
      integer i
      integer norder
      real(kind=Rkind) p1
      real(kind=Rkind) p2
      real(kind=Rkind) x
!
      do i = 1, 10

        call herrec ( p2, dp2, p1, x, norder )

        d = p2 / dp2
        x = x - d

        if ( abs ( d ) .le. eps * ( abs ( x ) + 1.0 ) ) then
          return
        end if

      end do

      return
      end subroutine herroot
      subroutine herset ( norder, xtab, weight )
!
!***********************************************************************
!
!! HERSET sets abscissas and weights for Hermite quadrature.
!
!
!  Integration interval:
!
!    ( -Infinity, +Infinity )
!
!  Weight function:
!
!    EXP ( - X**2 ).
!
!  Integral to approximate:
!
!    INTEGRAL ( -INFINITY < X < +INFINITY ) EXP ( - X**2 ) * F(X) dX.
!
!  Approximate integral:
!
!    SUM ( I = 1 to NORDER ) WEIGHT(I) * F ( XTAB(I) ).
!
!  Reference:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964.
!
!    Vladimir Krylov,
!    Approximate Calculation of Integrals,
!    MacMillan, 1962.
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Modified:
!
!    16 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!    NORDER must be between 1 and 20.
!
!    Output, real(kind=Rkind) XTAB(NORDER), the abscissas of the rule,
!    which are symmetrically placed around 0.
!
!    Output, real(kind=Rkind) WEIGHT(NORDER), the weights of the rule.
!    The weights are positive and symmetric, and should sum
!    to SQRT(PI).
!
      real(kind=Rkind) pi
      parameter ( pi = 3.141592653589793238462643383279D+00 )
!
      integer norder
!
      real(kind=Rkind) xtab(norder)
      real(kind=Rkind) weight(norder)
!
      if ( norder .eq. 1 ) then

        xtab(1) = 0.0 D+00

        weight(1) = sqrt ( pi )

      else if ( norder .eq. 2 )then

        xtab(1) = - 0.707106781186547524400844362105D+00
        xtab(2) =   0.707106781186547524400844362105D+00

        weight(1) = 0.886226925452758013649083741671D+00
        weight(2) = 0.886226925452758013649083741671D+00

      else if ( norder .eq. 3 ) then

        xtab(1) = - 0.122474487139158904909864203735D+01
        xtab(2) =   0.0 D+00
        xtab(3) =   0.122474487139158904909864203735D+01

        weight(1) = 0.295408975150919337883027913890D+00
        weight(2) = 0.118163590060367735153211165556D+01
        weight(3) = 0.295408975150919337883027913890D+00

      else if ( norder .eq. 4 ) then

        xtab(1) = - 0.165068012388578455588334111112D+01
        xtab(2) = - 0.524647623275290317884060253835D+00
        xtab(3) =   0.524647623275290317884060253835D+00
        xtab(4) =   0.165068012388578455588334111112D+01

        weight(1) = 0.813128354472451771430345571899D-01
        weight(2) = 0.804914090005512836506049184481D+00
        weight(3) = 0.804914090005512836506049184481D+00
        weight(4) = 0.813128354472451771430345571899D-01

      else if ( norder .eq. 5 ) then

        xtab(1) = - 0.202018287045608563292872408814 D+01
        xtab(2) = - 0.958572464613818507112770593893 D+00
        xtab(3) =   0.0 D+00
        xtab(4) =   0.958572464613818507112770593893 D+00
        xtab(5) =   0.202018287045608563292872408814 D+01

        weight(1) = 0.199532420590459132077434585942 D-01
        weight(2) = 0.393619323152241159828495620852 D+00
        weight(3) = 0.945308720482941881225689324449 D+00
        weight(4) = 0.393619323152241159828495620852 D+00
        weight(5) = 0.199532420590459132077434585942 D-01

      else if ( norder .eq. 6 ) then

        xtab(1) = - 0.235060497367449222283392198706 D+01
        xtab(2) = - 0.133584907401369694971489528297 D+01
        xtab(3) = - 0.436077411927616508679215948251 D+00
        xtab(4) =   0.436077411927616508679215948251 D+00
        xtab(5) =   0.133584907401369694971489528297 D+01
        xtab(6) =   0.235060497367449222283392198706 D+01

        weight(1) = 0.453000990550884564085747256463 D-02
        weight(2) = 0.157067320322856643916311563508 D+00
        weight(3) = 0.724629595224392524091914705598 D+00
        weight(4) = 0.724629595224392524091914705598 D+00
        weight(5) = 0.157067320322856643916311563508 D+00
        weight(6) = 0.453000990550884564085747256463 D-02

      else if ( norder .eq. 7 ) then

        xtab(1) = - 0.265196135683523349244708200652 D+01
        xtab(2) = - 0.167355162876747144503180139830 D+01
        xtab(3) = - 0.816287882858964663038710959027 D+00
        xtab(4) =   0.0 D+00
        xtab(5) =   0.816287882858964663038710959027 D+00
        xtab(6) =   0.167355162876747144503180139830 D+01
        xtab(7) =   0.265196135683523349244708200652 D+01

        weight(1) = 0.971781245099519154149424255939 D-03
        weight(2) = 0.545155828191270305921785688417 D-01
        weight(3) = 0.425607252610127800520317466666 D+00
        weight(4) = 0.810264617556807326764876563813 D+00
        weight(5) = 0.425607252610127800520317466666 D+00
        weight(6) = 0.545155828191270305921785688417 D-01
        weight(7) = 0.971781245099519154149424255939 D-03

      else if ( norder .eq. 8 ) then

        xtab(1) = - 0.293063742025724401922350270524 D+01
        xtab(2) = - 0.198165675669584292585463063977 D+01
        xtab(3) = - 0.115719371244678019472076577906 D+01
        xtab(4) = - 0.381186990207322116854718885584 D+00
        xtab(5) =   0.381186990207322116854718885584 D+00
        xtab(6) =   0.115719371244678019472076577906 D+01
        xtab(7) =   0.198165675669584292585463063977 D+01
        xtab(8) =   0.293063742025724401922350270524 D+01

        weight(1) = 0.199604072211367619206090452544 D-03
        weight(2) = 0.170779830074134754562030564364 D-01
        weight(3) = 0.207802325814891879543258620286 D+00
        weight(4) = 0.661147012558241291030415974496 D+00
        weight(5) = 0.661147012558241291030415974496 D+00
        weight(6) = 0.207802325814891879543258620286 D+00
        weight(7) = 0.170779830074134754562030564364 D-01
        weight(8) = 0.199604072211367619206090452544 D-03

      else if ( norder .eq. 9 ) then

        xtab(1) = - 0.319099320178152760723004779538 D+01
        xtab(2) = - 0.226658058453184311180209693284 D+01
        xtab(3) = - 0.146855328921666793166701573925 D+01
        xtab(4) = - 0.723551018752837573322639864579 D+00
        xtab(5) =   0.0 D+00
        xtab(6) =   0.723551018752837573322639864579 D+00
        xtab(7) =   0.146855328921666793166701573925 D+01
        xtab(8) =   0.226658058453184311180209693284 D+01
        xtab(9) =   0.319099320178152760723004779538 D+01

        weight(1) = 0.396069772632643819045862946425 D-04
        weight(2) = 0.494362427553694721722456597763 D-02
        weight(3) = 0.884745273943765732879751147476 D-01
        weight(4) = 0.432651559002555750199812112956 D+00
        weight(5) = 0.720235215606050957124334723389 D+00
        weight(6) = 0.432651559002555750199812112956 D+00
        weight(7) = 0.884745273943765732879751147476 D-01
        weight(8) = 0.494362427553694721722456597763 D-02
        weight(9) = 0.396069772632643819045862946425 D-04

      else if ( norder .eq. 10 ) then

        xtab(1) =  - 0.343615911883773760332672549432 D+01
        xtab(2) =  - 0.253273167423278979640896079775 D+01
        xtab(3) =  - 0.175668364929988177345140122011 D+01
        xtab(4) =  - 0.103661082978951365417749191676 D+01
        xtab(5) =  - 0.342901327223704608789165025557 D+00
        xtab(6) =    0.342901327223704608789165025557 D+00
        xtab(7) =    0.103661082978951365417749191676 D+01
        xtab(8) =    0.175668364929988177345140122011 D+01
        xtab(9) =    0.253273167423278979640896079775 D+01
        xtab(10) =   0.343615911883773760332672549432 D+01

        weight(1) =  0.764043285523262062915936785960 D-05
        weight(2) =  0.134364574678123269220156558585 D-02
        weight(3) =  0.338743944554810631361647312776 D-01
        weight(4) =  0.240138611082314686416523295006 D+00
        weight(5) =  0.610862633735325798783564990433 D+00
        weight(6) =  0.610862633735325798783564990433 D+00
        weight(7) =  0.240138611082314686416523295006 D+00
        weight(8) =  0.338743944554810631361647312776 D-01
        weight(9) =  0.134364574678123269220156558585 D-02
        weight(10) = 0.764043285523262062915936785960 D-05

      else if ( norder .eq. 11 ) then

        xtab(1) =  - 0.366847084655958251845837146485 D+01
        xtab(2) =  - 0.278329009978165177083671870152 D+01
        xtab(3) =  - 0.202594801582575533516591283121 D+01
        xtab(4) =  - 0.132655708449493285594973473558 D+01
        xtab(5) =  - 0.656809566882099765024611575383 D+00
        xtab(6) =    0.0 D+00
        xtab(7) =    0.656809566882099765024611575383 D+00
        xtab(8) =    0.132655708449493285594973473558 D+01
        xtab(9) =    0.202594801582575533516591283121 D+01
        xtab(10) =   0.278329009978165177083671870152 D+01
        xtab(11) =   0.366847084655958251845837146485 D+01

        weight(1) =  0.143956039371425822033088366032 D-05
        weight(2) =  0.346819466323345510643413772940 D-03
        weight(3) =  0.119113954449115324503874202916 D-01
        weight(4) =  0.117227875167708503381788649308 D+00
        weight(5) =  0.429359752356125028446073598601 D+00
        weight(6) =  0.654759286914591779203940657627 D+00
        weight(7) =  0.429359752356125028446073598601 D+00
        weight(8) =  0.117227875167708503381788649308 D+00
        weight(9) =  0.119113954449115324503874202916 D-01
        weight(10) = 0.346819466323345510643413772940 D-03
        weight(11) = 0.143956039371425822033088366032 D-05

      else if ( norder .eq. 12 ) then

        xtab(1) =  - 0.388972489786978191927164274724 D+01
        xtab(2) =  - 0.302063702512088977171067937518 D+01
        xtab(3) =  - 0.227950708050105990018772856942 D+01
        xtab(4) =  - 0.159768263515260479670966277090 D+01
        xtab(5) =  - 0.947788391240163743704578131060 D+00
        xtab(6) =  - 0.314240376254359111276611634095 D+00
        xtab(7) =    0.314240376254359111276611634095 D+00
        xtab(8) =    0.947788391240163743704578131060 D+00
        xtab(9) =    0.159768263515260479670966277090 D+01
        xtab(10) =   0.227950708050105990018772856942 D+01
        xtab(11) =   0.302063702512088977171067937518 D+01
        xtab(12) =   0.388972489786978191927164274724 D+01

        weight(1) =  0.265855168435630160602311400877 D-06
        weight(2) =  0.857368704358785865456906323153 D-04
        weight(3) =  0.390539058462906185999438432620 D-02
        weight(4) =  0.516079856158839299918734423606 D-01
        weight(5) =  0.260492310264161129233396139765 D+00
        weight(6) =  0.570135236262479578347113482275 D+00
        weight(7) =  0.570135236262479578347113482275 D+00
        weight(8) =  0.260492310264161129233396139765 D+00
        weight(9) =  0.516079856158839299918734423606 D-01
        weight(10) = 0.390539058462906185999438432620 D-02
        weight(11) = 0.857368704358785865456906323153 D-04
        weight(12) = 0.265855168435630160602311400877 D-06

      else if ( norder .eq. 13 ) then

        xtab(1) =  - 0.410133759617863964117891508007 D+01
        xtab(2) =  - 0.324660897837240998812205115236 D+01
        xtab(3) =  - 0.251973568567823788343040913628 D+01
        xtab(4) =  - 0.185310765160151214200350644316 D+01
        xtab(5) =  - 0.122005503659074842622205526637 D+01
        xtab(6) =  - 0.605763879171060113080537108602 D+00
        xtab(7) =    0.0 D+00
        xtab(8) =    0.605763879171060113080537108602 D+00
        xtab(9) =    0.122005503659074842622205526637 D+01
        xtab(10) =   0.185310765160151214200350644316 D+01
        xtab(11) =   0.251973568567823788343040913628 D+01
        xtab(12) =   0.324660897837240998812205115236 D+01
        xtab(13) =   0.410133759617863964117891508007 D+01

        weight(1) =  0.482573185007313108834997332342 D-07
        weight(2) =  0.204303604027070731248669432937 D-04
        weight(3) =  0.120745999271938594730924899224 D-02
        weight(4) =  0.208627752961699392166033805050 D-01
        weight(5) =  0.140323320687023437762792268873 D+00
        weight(6) =  0.421616296898543221746893558568 D+00
        weight(7) =  0.604393187921161642342099068579 D+00
        weight(8) =  0.421616296898543221746893558568 D+00
        weight(9) =  0.140323320687023437762792268873 D+00
        weight(10) = 0.208627752961699392166033805050 D-01
        weight(11) = 0.120745999271938594730924899224 D-02
        weight(12) = 0.204303604027070731248669432937 D-04
        weight(13) = 0.482573185007313108834997332342 D-07

      else if ( norder .eq. 14 ) then

        xtab(1) =  - 0.430444857047363181262129810037 D+01
        xtab(2) =  - 0.346265693360227055020891736115 D+01
        xtab(3) =  - 0.274847072498540256862499852415 D+01
        xtab(4) =  - 0.209518325850771681573497272630 D+01
        xtab(5) =  - 0.147668273114114087058350654421 D+01
        xtab(6) =  - 0.878713787329399416114679311861 D+00
        xtab(7) =  - 0.291745510672562078446113075799 D+00
        xtab(8) =    0.291745510672562078446113075799 D+00
        xtab(9) =    0.878713787329399416114679311861 D+00
        xtab(10) =   0.147668273114114087058350654421 D+01
        xtab(11) =   0.209518325850771681573497272630 D+01
        xtab(12) =   0.274847072498540256862499852415 D+01
        xtab(13) =   0.346265693360227055020891736115 D+01
        xtab(14) =   0.430444857047363181262129810037 D+01

        weight(1) =  0.862859116812515794532041783429 D-08
        weight(2) =  0.471648435501891674887688950105 D-05
        weight(3) =  0.355092613551923610483661076691 D-03
        weight(4) =  0.785005472645794431048644334608 D-02
        weight(5) =  0.685055342234652055387163312367 D-01
        weight(6) =  0.273105609064246603352569187026 D+00
        weight(7) =  0.536405909712090149794921296776 D+00
        weight(8) =  0.536405909712090149794921296776 D+00
        weight(9) =  0.273105609064246603352569187026 D+00
        weight(10) = 0.685055342234652055387163312367 D-01
        weight(11) = 0.785005472645794431048644334608 D-02
        weight(12) = 0.355092613551923610483661076691 D-03
        weight(13) = 0.471648435501891674887688950105 D-05
        weight(14) = 0.862859116812515794532041783429 D-08

      else if ( norder .eq. 15 ) then

        xtab(1) =  - 0.449999070730939155366438053053 D+01
        xtab(2) =  - 0.366995037340445253472922383312 D+01
        xtab(3) =  - 0.296716692790560324848896036355 D+01
        xtab(4) =  - 0.232573248617385774545404479449 D+01
        xtab(5) =  - 0.171999257518648893241583152515 D+01
        xtab(6) =  - 0.113611558521092066631913490556 D+01
        xtab(7) =  - 0.565069583255575748526020337198 D+00
        xtab(8) =    0.0 D+00
        xtab(9) =    0.565069583255575748526020337198 D+00
        xtab(10) =   0.113611558521092066631913490556 D+01
        xtab(11) =   0.171999257518648893241583152515 D+01
        xtab(12) =   0.232573248617385774545404479449 D+01
        xtab(13) =   0.296716692790560324848896036355 D+01
        xtab(14) =   0.366995037340445253472922383312 D+01
        xtab(15) =   0.449999070730939155366438053053 D+01

        weight(1) =  0.152247580425351702016062666965 D-08
        weight(2) =  0.105911554771106663577520791055 D-05
        weight(3) =  0.100004441232499868127296736177 D-03
        weight(4) =  0.277806884291277589607887049229 D-02
        weight(5) =  0.307800338725460822286814158758 D-01
        weight(6) =  0.158488915795935746883839384960 D+00
        weight(7) =  0.412028687498898627025891079568 D+00
        weight(8) =  0.564100308726417532852625797340 D+00
        weight(9) =  0.412028687498898627025891079568 D+00
        weight(10) = 0.158488915795935746883839384960 D+00
        weight(11) = 0.307800338725460822286814158758 D-01
        weight(12) = 0.277806884291277589607887049229 D-02
        weight(13) = 0.100004441232499868127296736177 D-03
        weight(14) = 0.105911554771106663577520791055 D-05
        weight(15) = 0.152247580425351702016062666965 D-08

      else if ( norder .eq. 16 ) then

        xtab(1) =  - 0.468873893930581836468849864875 D+01
        xtab(2) =  - 0.386944790486012269871942409801 D+01
        xtab(3) =  - 0.317699916197995602681399455926 D+01
        xtab(4) =  - 0.254620215784748136215932870545 D+01
        xtab(5) =  - 0.195178799091625397743465541496 D+01
        xtab(6) =  - 0.138025853919888079637208966969 D+01
        xtab(7) =  - 0.822951449144655892582454496734 D+00
        xtab(8) =  - 0.273481046138152452158280401965 D+00
        xtab(9) =    0.273481046138152452158280401965 D+00
        xtab(10) =   0.822951449144655892582454496734 D+00
        xtab(11) =   0.138025853919888079637208966969 D+01
        xtab(12) =   0.195178799091625397743465541496 D+01
        xtab(13) =   0.254620215784748136215932870545 D+01
        xtab(14) =   0.317699916197995602681399455926 D+01
        xtab(15) =   0.386944790486012269871942409801 D+01
        xtab(16) =   0.468873893930581836468849864875 D+01

        weight(1) =  0.265480747401118224470926366050 D-09
        weight(2) =  0.232098084486521065338749423185 D-06
        weight(3) =  0.271186009253788151201891432244 D-04
        weight(4) =  0.932284008624180529914277305537 D-03
        weight(5) =  0.128803115355099736834642999312 D-01
        weight(6) =  0.838100413989858294154207349001 D-01
        weight(7) =  0.280647458528533675369463335380 D+00
        weight(8) =  0.507929479016613741913517341791 D+00
        weight(9) =  0.507929479016613741913517341791 D+00
        weight(10) = 0.280647458528533675369463335380 D+00
        weight(11) = 0.838100413989858294154207349001 D-01
        weight(12) = 0.128803115355099736834642999312 D-01
        weight(13) = 0.932284008624180529914277305537 D-03
        weight(14) = 0.271186009253788151201891432244 D-04
        weight(15) = 0.232098084486521065338749423185 D-06
        weight(16) = 0.265480747401118224470926366050 D-09

      else if ( norder .eq. 17 ) then

        xtab(1) =  - 0.487134519367440308834927655662 D+01
        xtab(2) =  - 0.406194667587547430689245559698 D+01
        xtab(3) =  - 0.337893209114149408338327069289 D+01
        xtab(4) =  - 0.275776291570388873092640349574 D+01
        xtab(5) =  - 0.217350282666662081927537907149 D+01
        xtab(6) =  - 0.161292431422123133311288254454 D+01
        xtab(7) =  - 0.106764872574345055363045773799 D+01
        xtab(8) =  - 0.531633001342654731349086553718 D+00
        xtab(9) =    0.0 D+00
        xtab(10) =   0.531633001342654731349086553718 D+00
        xtab(11) =   0.106764872574345055363045773799 D+01
        xtab(12) =   0.161292431422123133311288254454 D+01
        xtab(13) =   0.217350282666662081927537907149 D+01
        xtab(14) =   0.275776291570388873092640349574 D+01
        xtab(15) =   0.337893209114149408338327069289 D+01
        xtab(16) =   0.406194667587547430689245559698 D+01
        xtab(17) =   0.487134519367440308834927655662 D+01

        weight(1) =  0.458057893079863330580889281222 D-10
        weight(2) =  0.497707898163079405227863353715 D-07
        weight(3) =  0.711228914002130958353327376218 D-05
        weight(4) =  0.298643286697753041151336643059 D-03
        weight(5) =  0.506734995762753791170069495879 D-02
        weight(6) =  0.409200341495762798094994877854 D-01
        weight(7) =  0.172648297670097079217645196219 D+00
        weight(8) =  0.401826469470411956577635085257 D+00
        weight(9) =  0.530917937624863560331883103379 D+00
        weight(10) = 0.401826469470411956577635085257 D+00
        weight(11) = 0.172648297670097079217645196219 D+00
        weight(12) = 0.409200341495762798094994877854 D-01
        weight(13) = 0.506734995762753791170069495879 D-02
        weight(14) = 0.298643286697753041151336643059 D-03
        weight(15) = 0.711228914002130958353327376218 D-05
        weight(16) = 0.497707898163079405227863353715 D-07
        weight(17) = 0.458057893079863330580889281222 D-10

      else if ( norder .eq. 18 ) then

        xtab(1) =  - 0.504836400887446676837203757885 D+01
        xtab(2) =  - 0.424811787356812646302342016090 D+01
        xtab(3) =  - 0.357376906848626607950067599377 D+01
        xtab(4) =  - 0.296137750553160684477863254906 D+01
        xtab(5) =  - 0.238629908916668600026459301424 D+01
        xtab(6) =  - 0.183553160426162889225383944409 D+01
        xtab(7) =  - 0.130092085838961736566626555439 D+01
        xtab(8) =  - 0.776682919267411661316659462284 D+00
        xtab(9) =  - 0.258267750519096759258116098711 D+00
        xtab(10) =   0.258267750519096759258116098711 D+00
        xtab(11) =   0.776682919267411661316659462284 D+00
        xtab(12) =   0.130092085838961736566626555439 D+01
        xtab(13) =   0.183553160426162889225383944409 D+01
        xtab(14) =   0.238629908916668600026459301424 D+01
        xtab(15) =   0.296137750553160684477863254906 D+01
        xtab(16) =   0.357376906848626607950067599377 D+01
        xtab(17) =   0.424811787356812646302342016090 D+01
        xtab(18) =   0.504836400887446676837203757885 D+01

        weight(1) =  0.782819977211589102925147471012 D-11
        weight(2) =  0.104672057957920824443559608435 D-07
        weight(3) =  0.181065448109343040959702385911 D-05
        weight(4) =  0.918112686792940352914675407371 D-04
        weight(5) =  0.188852263026841789438175325426 D-02
        weight(6) =  0.186400423875446519219315221973 D-01
        weight(7) =  0.973017476413154293308537234155 D-01
        weight(8) =  0.284807285669979578595606820713 D+00
        weight(9) =  0.483495694725455552876410522141 D+00
        weight(10) = 0.483495694725455552876410522141 D+00
        weight(11) = 0.284807285669979578595606820713 D+00
        weight(12) = 0.973017476413154293308537234155 D-01
        weight(13) = 0.186400423875446519219315221973 D-01
        weight(14) = 0.188852263026841789438175325426 D-02
        weight(15) = 0.918112686792940352914675407371 D-04
        weight(16) = 0.181065448109343040959702385911 D-05
        weight(17) = 0.104672057957920824443559608435 D-07
        weight(18) = 0.782819977211589102925147471012 D-11

      else if ( norder .eq. 19 ) then

        xtab(1) =  - 0.522027169053748216460967142500 D+01
        xtab(2) =  - 0.442853280660377943723498532226 D+01
        xtab(3) =  - 0.376218735196402009751489394104 D+01
        xtab(4) =  - 0.315784881834760228184318034120 D+01
        xtab(5) =  - 0.259113378979454256492128084112 D+01
        xtab(6) =  - 0.204923170985061937575050838669 D+01
        xtab(7) =  - 0.152417061939353303183354859367 D+01
        xtab(8) =  - 0.101036838713431135136859873726 D+01
        xtab(9) =  - 0.503520163423888209373811765050 D+00
        xtab(10) =   0.0 D+00
        xtab(11) =   0.503520163423888209373811765050 D+00
        xtab(12) =   0.101036838713431135136859873726 D+01
        xtab(13) =   0.152417061939353303183354859367 D+01
        xtab(14) =   0.204923170985061937575050838669 D+01
        xtab(15) =   0.259113378979454256492128084112 D+01
        xtab(16) =   0.315784881834760228184318034120 D+01
        xtab(17) =   0.376218735196402009751489394104 D+01
        xtab(18) =   0.442853280660377943723498532226 D+01
        xtab(19) =   0.522027169053748216460967142500 D+01

        weight(1) =  0.132629709449851575185289154385 D-11
        weight(2) =  0.216305100986355475019693077221 D-08
        weight(3) =  0.448824314722312295179447915594 D-06
        weight(4) =  0.272091977631616257711941025214 D-04
        weight(5) =  0.670877521407181106194696282100 D-03
        weight(6) =  0.798886677772299020922211491861 D-02
        weight(7) =  0.508103869090520673569908110358 D-01
        weight(8) =  0.183632701306997074156148485766 D+00
        weight(9) =  0.391608988613030244504042313621 D+00
        weight(10) = 0.502974888276186530840731361096 D+00
        weight(11) = 0.391608988613030244504042313621 D+00
        weight(12) = 0.183632701306997074156148485766 D+00
        weight(13) = 0.508103869090520673569908110358 D-01
        weight(14) = 0.798886677772299020922211491861 D-02
        weight(15) = 0.670877521407181106194696282100 D-03
        weight(16) = 0.272091977631616257711941025214 D-04
        weight(17) = 0.448824314722312295179447915594 D-06
        weight(18) = 0.216305100986355475019693077221 D-08
        weight(19) = 0.132629709449851575185289154385 D-11

      else if ( norder .eq. 20 ) then

        xtab(1) =  - 0.538748089001123286201690041068 D+01
        xtab(2) =  - 0.460368244955074427307767524898 D+01
        xtab(3) =  - 0.394476404011562521037562880052 D+01
        xtab(4) =  - 0.334785456738321632691492452300 D+01
        xtab(5) =  - 0.278880605842813048052503375640 D+01
        xtab(6) =  - 0.225497400208927552308233334473 D+01
        xtab(7) =  - 0.173853771211658620678086566214 D+01
        xtab(8) =  - 0.123407621539532300788581834696 D+01
        xtab(9) =  - 0.737473728545394358705605144252 D+00
        xtab(10) = - 0.245340708300901249903836530634 D+00
        xtab(11) =   0.245340708300901249903836530634 D+00
        xtab(12) =   0.737473728545394358705605144252 D+00
        xtab(13) =   0.123407621539532300788581834696 D+01
        xtab(14) =   0.173853771211658620678086566214 D+01
        xtab(15) =   0.225497400208927552308233334473 D+01
        xtab(16) =   0.278880605842813048052503375640 D+01
        xtab(17) =   0.334785456738321632691492452300 D+01
        xtab(18) =   0.394476404011562521037562880052 D+01
        xtab(19) =   0.460368244955074427307767524898 D+01
        xtab(20) =   0.538748089001123286201690041068 D+01

        weight(1) =  0.222939364553415129252250061603 D-12
        weight(2) =  0.439934099227318055362885145547 D-09
        weight(3) =  0.108606937076928169399952456345 D-06
        weight(4) =  0.780255647853206369414599199965 D-05
        weight(5) =  0.228338636016353967257145917963 D-03
        weight(6) =  0.324377334223786183218324713235 D-02
        weight(7) =  0.248105208874636108821649525589 D-01
        weight(8) =  0.109017206020023320013755033535 D+00
        weight(9) =  0.286675505362834129719659706228 D+00
        weight(10) = 0.462243669600610089650328639861 D+00
        weight(11) = 0.462243669600610089650328639861 D+00
        weight(12) = 0.286675505362834129719659706228 D+00
        weight(13) = 0.109017206020023320013755033535 D+00
        weight(14) = 0.248105208874636108821649525589 D-01
        weight(15) = 0.324377334223786183218324713235 D-02
        weight(16) = 0.228338636016353967257145917963 D-03
        weight(17) = 0.780255647853206369414599199965 D-05
        weight(18) = 0.108606937076928169399952456345 D-06
        weight(19) = 0.439934099227318055362885145547 D-09
        weight(20) = 0.222939364553415129252250061603 D-12

      else

        write ( *, * ) ' '
        write ( *, * ) 'HERSET - Fatal error!'
        write ( *, * ) '  Illegal value of NORDER = ', norder
        write ( *, * ) '  Legal values are 1 to 20.'
        stop

      end if

      return
      end subroutine herset
      subroutine jaccom ( norder, xtab, weight, alpha, beta, b, c,      &
        csx, csa, tsx, tsa )
!
!***********************************************************************
!
!! JACCOM computes the abscissa and weights for Gauss-Jacobi quadrature.
!
!
!  Integration interval:
!
!    [ -1, 1 ]
!
!  Weight function:
!
!    1.
!
!  Integral to approximate:
!
!    INTEGRAL ( -1 <= X <= 1 ) (1+X)**ALPHA * (1-X)**BETA * F(X) dX.
!
!  Approximate integral:
!
!    SUM ( I = 1 to NORDER ) WEIGHT(I) * F ( XTAB(I) ).
!
!  Reference:
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!  Modified:
!
!    19 September 1998
!
!  Parameters:
!
!    Input, integer NORDER, the order of the quadrature rule to be computed.
!
!    Output, real(kind=Rkind) XTAB(NORDER), the Gauss-Jacobi abscissas.
!
!    Output, real(kind=Rkind) WEIGHT(NORDER), the Gauss-Jacobi weights.
!
!    Input, real(kind=Rkind) ALPHA, BETA, the exponents of (1+X) and
!    (1-X) in the quadrature rule.  For simple Gauss-Legendre quadrature,
!    set ALPHA = BETA = 0.0.
!
!    Workspace, real(kind=Rkind) B(NORDER), C(NORDER), the recursion
!    coefficients.
!
!    Output, real(kind=Rkind) CSX, the sum of the computed abscissas.
!
!    Output, real(kind=Rkind) CSA, the sum of the computed weights.
!
!    Output, real(kind=Rkind) TSX, the true sum of the abscissas.
!
!    Output, real(kind=Rkind) TSA, the true sum of the weights.
!
      integer norder
!
      real(kind=Rkind) alpha
      real(kind=Rkind) an
      real(kind=Rkind) b(norder)
      real(kind=Rkind) beta
      real(kind=Rkind) bn
      real(kind=Rkind) c(norder)
      real(kind=Rkind) cc
      real(kind=Rkind) csa
      real(kind=Rkind) csx
      real(kind=Rkind) delta
      real(kind=Rkind) dp2
      real(kind=Rkind) log_gamma
      integer i
      real(kind=Rkind) p1
      real(kind=Rkind) r1
      real(kind=Rkind) r2
      real(kind=Rkind) r3
      real(kind=Rkind) temp
      real(kind=Rkind) tsa
      real(kind=Rkind) tsx
      real(kind=Rkind) weight(norder)
      real(kind=Rkind) x
      real(kind=Rkind) xtab(norder)
!
!  Set the recursion coefficients.
!
      do i = 1, norder

        b(i) = ( alpha + beta ) * ( beta - alpha ) /                    &
          ( ( alpha + beta + dble ( 2 * i ) )                           &
          * ( alpha + beta + dble ( 2 * i - 2 ) ) )

        c(i) = 4.0D+00 * dble ( i - 1 ) * ( alpha + dble ( i - 1 ) )    &
          * ( beta + dble ( i - 1 ) )                                   &
          * ( alpha + beta + dble ( i - 1 ) ) /                         &
          ( ( alpha + beta + dble ( 2 * i - 1 ) )                       &
          * ( alpha + beta + dble ( 2 * i - 2 ) )**2                    &
          * ( alpha + beta + dble ( 2 * i - 3 ) ) )

      end do

      csx = 0.0
      csa = 0.0

      delta = exp (                                                     &
          log_gamma ( alpha + 1.0 )                                     &
        + log_gamma ( beta + 1.0 )                                      &
        + log_gamma ( alpha + beta + 2.0 ) )

      cc = delta * 2.0**( alpha + beta + 1.0 )

      tsx = dble ( norder ) * ( beta - alpha ) /                        &
        ( alpha + beta + 2.0 * dble ( norder ) )

      tsa = cc

      do i = 2, norder
        cc  = cc * c(i)
      end do

      do i = 1, norder

        if ( i .eq. 1 ) then

          an = alpha / dble ( norder )
          bn = beta / dble ( norder )

          r1 = ( 1.0 + alpha ) * (                                      &
            2.78 / ( 4.0 + dble ( norder**2 ) )                         &
            + 0.768 * an / dble ( norder ) )

          r2 = 1.0 + 1.48 * an + 0.96 * bn + 0.452 * an**2              &
            + 0.83 * an * bn

          x = ( r2 - r1 ) / r2

        else if ( i .eq. 2 ) then

          r1 = ( 4.1 + alpha ) / ( ( 1.0 + alpha ) *                    &
            ( 1.0 + 0.156 * alpha ) )

          r2 = 1.0 + 0.06 * ( dble ( norder ) - 8.0 ) *                 &
            ( 1.0 + 0.12 * alpha ) / dble ( norder )

          r3 = 1.0 + 0.012 * beta *                                     &
            ( 1.0 + 0.25 * abs ( alpha ) ) / dble ( norder )

          x = x - r1 * r2 * r3 * ( 1.0 - x )

        else if ( i .eq. 3 ) then

          r1 = ( 1.67 + 0.28 * alpha ) / ( 1.0 + 0.37 * alpha )

          r2 = 1.0 + 0.22 * ( dble ( norder ) - 8.0 ) / dble ( norder )

          r3 = 1.0 + 8.0 * beta /                                       &
            ( ( 6.28 + beta ) * dble ( norder**2 ) )

          x = x - r1 * r2 * r3 * ( xtab(1) - x )

        else if ( i .lt. norder - 1 ) then

          x = 3.0 * xtab(i-1) - 3.0 * xtab(i-2) + xtab(i-3)

        else if ( i .eq. norder - 1 ) then

          r1 = ( 1.0 + 0.235 * beta ) / ( 0.766 + 0.119 * beta )

          r2 = 1.0 / ( 1.0 + 0.639 * ( dble ( norder ) - 4.0 )          &
            / ( 1.0 + 0.71 * ( dble ( norder ) - 4.0 ) ) )

          r3 = 1.0 / ( 1.0 + 20.0 * alpha                               &
            / ( ( 7.5 + alpha ) * dble ( norder**2 ) ) )

          x = x + r1 * r2 * r3 * ( x - xtab(i-2) )

        else if ( i .eq. norder ) then

          r1 = ( 1.0 + 0.37 * beta ) / ( 1.67 + 0.28 * beta )

          r2 = 1.0 / ( 1.0 + 0.22 * ( dble ( norder ) - 8.0 )           &
            / dble ( norder ) )

          r3 = 1.0 / ( 1.0 + 8.0 * alpha /                              &
            ( ( 6.28 + alpha ) * dble ( norder**2 ) ) )

          x = x + r1 * r2 * r3 * ( x - xtab(i-2) )

        end if

        call jacroot ( x, norder, alpha, beta, dp2, p1, b, c )

        xtab(i) = x
        weight(i) = cc / ( dp2 * p1 )

        csx = csx + xtab(i)
        csa = csa + weight(i)

      end do
!
!  Reverse the order of the XTAB values.
!
      do i = 1, norder/2
        temp = xtab(i)
        xtab(i) = xtab(norder+1-i)
        xtab(norder+1-i) = temp
      end do

      return
      end subroutine jaccom
      subroutine jacrec ( p2, dp2, p1, x, norder, alpha, beta, b, c )
!
!***********************************************************************
!
!! JACREC finds the value and derivative of the NORDER-th Jacobi polynomial.
!
!
!  Reference:
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!  Modified:
!
!    19 September 1998
!
!  Parameters:
!
!    Output, real(kind=Rkind) P2, the value of J(NORDER)(X).
!
!    Output, real(kind=Rkind) DP2, the value of J'(NORDER)(X).
!
!    Output, real(kind=Rkind) P1, the value of J(NORDER-1)(X).
!
!    Input, real(kind=Rkind) X, the point at which polynomials are evaluated.
!
!    Input, integer NORDER, the order of the polynomial to be computed.
!
!    Input, real(kind=Rkind) ALPHA, BETA, the exponents of (1+X) and
!    (1-X) in the quadrature rule.
!
!    Input, real(kind=Rkind) B(NORDER), C(NORDER), the recursion
!    coefficients.
!
      integer norder
!
      real(kind=Rkind) alpha
      real(kind=Rkind) b(norder)
      real(kind=Rkind) beta
      real(kind=Rkind) c(norder)
      real(kind=Rkind) dp0
      real(kind=Rkind) dp1
      real(kind=Rkind) dp2
      integer i
      real(kind=Rkind) p0
      real(kind=Rkind) p1
      real(kind=Rkind) p2
      real(kind=Rkind) x
!
      p1 = 1.0
      dp1 = 0.0

      p2 = x + ( alpha - beta ) / ( alpha + beta + 2.0 )
      dp2 = 1.0

      do i = 2, norder

        p0 = p1
        dp0 = dp1

        p1 = p2
        dp1 = dp2

        p2 = ( x - b(i) ) * p1 - c(i) * p0
        dp2 = ( x - b(i) ) * dp1 + p1 - c(i) * dp0

      end do

      return
      end subroutine jacrec
      subroutine jacroot ( x, norder, alpha, beta, dp2, p1, b, c )
!
!***********************************************************************
!
!! JACROOT improves an approximate root of a Jacobi polynomial.
!
!
!  Reference:
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!  Modified:
!
!    19 September 1998
!
!  Parameters:
!
!    Input/output, real(kind=Rkind) X, the approximate root, which
!    should be improved on output.
!
!    Input, integer NORDER, the order of the polynomial to be computed.
!
!    Input, real(kind=Rkind) ALPHA, BETA, the exponents of (1+X) and
!    (1-X) in the quadrature rule.
!
!    Output, real(kind=Rkind) DP2, the value of J'(NORDER)(X).
!
!    Output, real(kind=Rkind) P1, the value of J(NORDER-1)(X).
!
!    Input, real(kind=Rkind) B(NORDER), C(NORDER), the recursion
!    coefficients.
!
      real(kind=Rkind) eps
      parameter ( eps = 1.0D-12 )
!
      integer norder
!
      real(kind=Rkind) alpha
      real(kind=Rkind) b(norder)
      real(kind=Rkind) beta
      real(kind=Rkind) c(norder)
      real(kind=Rkind) d
      real(kind=Rkind) dp2
      integer i
      real(kind=Rkind) p1
      real(kind=Rkind) p2
      real(kind=Rkind) x
!
      do i = 1, 10

        call jacrec ( p2, dp2, p1, x, norder, alpha, beta, b, c )

        d = p2 / dp2
        x = x - d

        if ( abs ( d ) .le. eps * ( abs ( x ) + 1.0 ) ) then
          return
        end if

      end do

      return
      end subroutine jacroot
      subroutine kronset ( norder, xtab, weight )
!
!***********************************************************************
!
!! KRONSET sets abscissas and weights for Gauss-Kronrod quadrature.
!
!
!  Integration interval:
!
!    [ -1, 1 ].
!
!  Weight function:
!
!    1.
!
!  Integral to approximate:
!
!    INTEGRAL ( -1 <= X <= 1 ) F(X) dX.
!
!  Approximate integral:
!
!    SUM ( I = 1 to NORDER ) WEIGHT(I) * F ( XTAB(I) ).
!
!  Note:
!
!    A Kronrod rule is used in conjunction with a lower order
!    Gauss rule, and provides an efficient error estimation.
!
!    The error may be estimated as the difference in the two integral
!    approximations.
!
!    The efficiency comes about because the Kronrod uses the abscissas
!    of the Gauss rule, thus saving on the number of function evaluations
!    necessary.  If the Kronrod rule were replaced by a Gauss rule of
!    the same order, a higher precision integral estimate would be
!    made, but the function would have to be evaluated at many more
!    points.
!
!    The Gauss Kronrod pair of rules involves an ( NORDER + 1 ) / 2
!    point Gauss-Legendre rule and an NORDER point Kronrod rule.
!    Thus, the 15 point Kronrod rule should be paired with the
!    Gauss-Legendre 7 point rule.
!
!  Reference:
!
!    R Piessens, E de Doncker-Kapenger, C Ueberhuber, D Kahaner,
!    QUADPACK, A Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983.
!
!  Modified:
!
!    16 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule, which may be
!    15, 21, 31 or 41, corresponding to Gauss-Legendre rules of
!    order 7, 10, 15 or 20.
!
!    Output, real(kind=Rkind) XTAB(NORDER), the abscissas of the rule, which
!    are symmetrically places in [-1,1].
!
!    Output, real(kind=Rkind) WEIGHT(NORDER), the weights of the rule.
!    The weights are positive, symmetric, and should sum to 2.
!
      integer norder
!
      real(kind=Rkind) weight(norder)
      real(kind=Rkind) xtab(norder)
!
      if ( norder .eq. 15 ) then

        xtab(1) =  - 0.9914553711208126 D+00
        xtab(2) =  - 0.9491079123427585 D+00
        xtab(3) =  - 0.8648644233597691 D+00
        xtab(4) =  - 0.7415311855993944 D+00
        xtab(5) =  - 0.5860872354676911 D+00
        xtab(6) =  - 0.4058451513773972 D+00
        xtab(7) =  - 0.2077849550789850 D+00
        xtab(8) =    0.0 D+00
        xtab(9) =    0.2077849550789850 D+00
        xtab(10) =   0.4058451513773972 D+00
        xtab(11) =   0.5860872354676911 D+00
        xtab(12) =   0.7415311855993944 D+00
        xtab(13) =   0.8648644233597691 D+00
        xtab(14) =   0.9491079123427585 D+00
        xtab(15) =   0.9914553711208126 D+00

        weight(1) =  0.2293532201052922 D-01
        weight(2) =  0.6309209262997855 D-01
        weight(3) =  0.1047900103222502 D+00
        weight(4) =  0.1406532597155259 D+00
        weight(5) =  0.1690047266392679 D+00
        weight(6) =  0.1903505780647854 D+00
        weight(7) =  0.2044329400752989 D+00
        weight(8) =  0.2094821410847278 D+00
        weight(9) =  0.2044329400752989 D+00
        weight(10) = 0.1903505780647854 D+00
        weight(11) = 0.1690047266392679 D+00
        weight(12) = 0.1406532597155259 D+00
        weight(13) = 0.1047900103222502 D+00
        weight(14) = 0.6309209262997855 D-01
        weight(15) = 0.2293532201052922 D-01

      else if ( norder .eq. 21 ) then

        xtab(1) =  - 0.9956571630258081 D+00
        xtab(2) =  - 0.9739065285171717 D+00
        xtab(3) =  - 0.9301574913557082 D+00
        xtab(4) =  - 0.8650633666889845 D+00
        xtab(5) =  - 0.7808177265864169 D+00
        xtab(6) =  - 0.6794095682990244 D+00
        xtab(7) =  - 0.5627571346686047 D+00
        xtab(8) =  - 0.4333953941292472 D+00
        xtab(9) =  - 0.2943928627014602 D+00
        xtab(10) = - 0.1488743389816312 D+00
        xtab(11) =   0.0 D+00
        xtab(12) =   0.1488743389816312 D+00
        xtab(13) =   0.2943928627014602 D+00
        xtab(14) =   0.4333953941292472 D+00
        xtab(15) =   0.5627571346686047 D+00
        xtab(16) =   0.6794095682990244 D+00
        xtab(17) =   0.7808177265864169 D+00
        xtab(18) =   0.8650633666889845 D+00
        xtab(19) =   0.9301574913557082 D+00
        xtab(20) =   0.9739065285171717 D+00
        xtab(21) =   0.9956571630258081 D+00

        weight(1) =  0.1169463886737187 D-01
        weight(2) =  0.3255816230796473 D-01
        weight(3) =  0.5475589657435200 D-01
        weight(4) =  0.7503967481091995 D-01
        weight(5) =  0.9312545458369761 D-01
        weight(6) =  0.1093871588022976 D+00
        weight(7) =  0.1234919762620659 D+00
        weight(8) =  0.1347092173114733 D+00
        weight(9) =  0.1427759385770601 D+00
        weight(10) = 0.1477391049013385 D+00
        weight(11) = 0.1494455540029169 D+00
        weight(12) = 0.1477391049013385 D+00
        weight(13) = 0.1427759385770601 D+00
        weight(14) = 0.1347092173114733 D+00
        weight(15) = 0.1234919762620659 D+00
        weight(16) = 0.1093871588022976 D+00
        weight(17) = 0.9312545458369761 D-01
        weight(18) = 0.7503967481091995 D-01
        weight(19) = 0.5475589657435200 D-01
        weight(20) = 0.3255816230796473 D-01
        weight(21) = 0.1169463886737187 D-01

      else if ( norder .eq. 31 ) then

        xtab(1) =  - 0.9980022986933971 D+00
        xtab(2) =  - 0.9879925180204854 D+00
        xtab(3) =  - 0.9677390756791391 D+00
        xtab(4) =  - 0.9372733924007059 D+00
        xtab(5) =  - 0.8972645323440819 D+00
        xtab(6) =  - 0.8482065834104272 D+00
        xtab(7) =  - 0.7904185014424659 D+00
        xtab(8) =  - 0.7244177313601700 D+00
        xtab(9) =  - 0.6509967412974170 D+00
        xtab(10) = - 0.5709721726085388 D+00
        xtab(11) = - 0.4850818636402397 D+00
        xtab(12) = - 0.3941513470775634 D+00
        xtab(13) = - 0.2991800071531688 D+00
        xtab(14) = - 0.2011940939974345 D+00
        xtab(15) = - 0.1011420669187175 D+00
        xtab(16) =   0.0 D+00
        xtab(17) =   0.1011420669187175 D+00
        xtab(18) =   0.2011940939974345 D+00
        xtab(19) =   0.2991800071531688 D+00
        xtab(20) =   0.3941513470775634 D+00
        xtab(21) =   0.4850818636402397 D+00
        xtab(22) =   0.5709721726085388 D+00
        xtab(23) =   0.6509967412974170 D+00
        xtab(24) =   0.7244177313601700 D+00
        xtab(25) =   0.7904185014424659 D+00
        xtab(26) =   0.8482065834104272 D+00
        xtab(27) =   0.8972645323440819 D+00
        xtab(28) =   0.9372733924007059 D+00
        xtab(29) =   0.9677390756791391 D+00
        xtab(30) =   0.9879925180204854 D+00
        xtab(31) =   0.9980022986933971 D+00

        weight(1) =  0.5377479872923349 D-02
        weight(2) =  0.1500794732931612 D-01
        weight(3) =  0.2546084732671532 D-01
        weight(4) =  0.3534636079137585 D-01
        weight(5) =  0.4458975132476488 D-01
        weight(6) =  0.5348152469092809 D-01
        weight(7) =  0.6200956780067064 D-01
        weight(8) =  0.6985412131872826 D-01
        weight(9) =  0.7684968075772038 D-01
        weight(10) = 0.8308050282313302 D-01
        weight(11) = 0.8856444305621177 D-01
        weight(12) = 0.9312659817082532 D-01
        weight(13) = 0.9664272698362368 D-01
        weight(14) = 0.9917359872179196 D-01
        weight(15) = 0.1007698455238756 D+00
        weight(16) = 0.1013300070147915 D+00
        weight(17) = 0.1007698455238756 D+00
        weight(18) = 0.9917359872179196 D-01
        weight(19) = 0.9664272698362368 D-01
        weight(20) = 0.9312659817082532 D-01
        weight(21) = 0.8856444305621177 D-01
        weight(22) = 0.8308050282313302 D-01
        weight(23) = 0.7684968075772038 D-01
        weight(24) = 0.6985412131872826 D-01
        weight(25) = 0.6200956780067064 D-01
        weight(26) = 0.5348152469092809 D-01
        weight(27) = 0.4458975132476488 D-01
        weight(28) = 0.3534636079137585 D-01
        weight(29) = 0.2546084732671532 D-01
        weight(30) = 0.1500794732931612 D-01
        weight(31) = 0.5377479872923349 D-02

      else if ( norder .eq. 41 ) then

        xtab(1) =  - 0.9988590315882777 D+00
        xtab(2) =  - 0.9931285991850949 D+00
        xtab(3) =  - 0.9815078774502503 D+00
        xtab(4) =  - 0.9639719272779138 D+00
        xtab(5) =  - 0.9408226338317548 D+00
        xtab(6) =  - 0.9122344282513259 D+00
        xtab(7) =  - 0.8782768112522820 D+00
        xtab(8) =  - 0.8391169718222188 D+00
        xtab(9) =  - 0.7950414288375512 D+00
        xtab(10) = - 0.7463319064601508 D+00
        xtab(11) = - 0.6932376563347514 D+00
        xtab(12) = - 0.6360536807265150 D+00
        xtab(13) = - 0.5751404468197103 D+00
        xtab(14) = - 0.5108670019508271 D+00
        xtab(15) = - 0.4435931752387251 D+00
        xtab(16) = - 0.3737060887154196 D+00
        xtab(17) = - 0.3016278681149130 D+00
        xtab(18) = - 0.2277858511416451 D+00
        xtab(19) = - 0.1526054652409227 D+00
        xtab(20) = - 0.7652652113349733 D-01
        xtab(21) =   0.0 D+00
        xtab(22) =   0.7652652113349733 D-01
        xtab(23) =   0.1526054652409227 D+00
        xtab(24) =   0.2277858511416451 D+00
        xtab(25) =   0.3016278681149130 D+00
        xtab(26) =   0.3737060887154196 D+00
        xtab(27) =   0.4435931752387251 D+00
        xtab(28) =   0.5108670019508271 D+00
        xtab(29) =   0.5751404468197103 D+00
        xtab(30) =   0.6360536807265150 D+00
        xtab(31) =   0.6932376563347514 D+00
        xtab(32) =   0.7463319064601508 D+00
        xtab(33) =   0.7950414288375512 D+00
        xtab(34) =   0.8391169718222188 D+00
        xtab(35) =   0.8782768112522820 D+00
        xtab(36) =   0.9122344282513259 D+00
        xtab(37) =   0.9408226338317548 D+00
        xtab(38) =   0.9639719272779138 D+00
        xtab(39) =   0.9815078774502503 D+00
        xtab(40) =   0.9931285991850949 D+00
        xtab(41) =   0.9988590315882777 D+00

        weight(1) =  0.3073583718520532 D-02
        weight(2) =  0.8600269855642942 D-02
        weight(3) =  0.1462616925697125 D-01
        weight(4) =  0.2038837346126652 D-01
        weight(5) =  0.2588213360495116 D-01
        weight(6) =  0.3128730677703280 D-01
        weight(7) =  0.3660016975820080 D-01
        weight(8) =  0.4166887332797369 D-01
        weight(9) =  0.4643482186749767 D-01
        weight(10) = 0.5094457392372869 D-01
        weight(11) = 0.5519510534828599 D-01
        weight(12) = 0.5911140088063957 D-01
        weight(13) = 0.6265323755478117 D-01
        weight(14) = 0.6583459713361842 D-01
        weight(15) = 0.6864867292852162 D-01
        weight(16) = 0.7105442355344407 D-01
        weight(17) = 0.7303069033278667 D-01
        weight(18) = 0.7458287540049919 D-01
        weight(19) = 0.7570449768455667 D-01
        weight(20) = 0.7637786767208074 D-01
        weight(21) = 0.7660071191799966 D-01
        weight(22) = 0.7637786767208074 D-01
        weight(23) = 0.7570449768455667 D-01
        weight(24) = 0.7458287540049919 D-01
        weight(25) = 0.7303069033278667 D-01
        weight(26) = 0.7105442355344407 D-01
        weight(27) = 0.6864867292852162 D-01
        weight(28) = 0.6583459713361842 D-01
        weight(29) = 0.6265323755478117 D-01
        weight(30) = 0.5911140088063957 D-01
        weight(31) = 0.5519510534828599 D-01
        weight(32) = 0.5094457392372869 D-01
        weight(33) = 0.4643482186749767 D-01
        weight(34) = 0.4166887332797369 D-01
        weight(35) = 0.3660016975820080 D-01
        weight(36) = 0.3128730677703280 D-01
        weight(37) = 0.2588213360495116 D-01
        weight(38) = 0.2038837346126652 D-01
        weight(39) = 0.1462616925697125 D-01
        weight(40) = 0.8600269855642942 D-02
        weight(41) = 0.3073583718520532 D-02

      else

        write ( *, * ) ' '
        write ( *, * ) 'KRONSET - Fatal error!'
        write ( *, * ) '  Illegal value of NORDER = ', norder
        write ( *, * ) '  Legal values are 15, 21, 31 or 41.'
        stop

      end if

      return
      end subroutine kronset
      subroutine lagcom ( norder, xtab, weight, alpha, b, c, csx,       &
        csa, tsx, tsa )
!
!***********************************************************************
!
!! LAGCOM computes the abscissa and weights for Gauss-Laguerre quadrature.
!
!
!  Integration interval:
!
!    [ 0, +Infinity )
!
!  Weight function:
!
!    EXP ( - X ).
!
!  Integral to approximate:
!
!    INTEGRAL ( 0 <= X < +INFINITY ) EXP ( - X ) * X**ALPHA * F(X) dX.
!
!  Approximate integral:
!
!    SUM ( I = 1 to NORDER ) WEIGHT(I) * F ( XTAB(I) ).
!
!  Reference:
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!  Modified:
!
!    19 September 1998
!
!  Parameters:
!
!    Input, integer NORDER, the order of the quadrature rule to be computed.
!
!    Output, real(kind=Rkind) XTAB(NORDER), the Gauss-Laguerre abscissas.
!
!    Output, real(kind=Rkind) WEIGHT(NORDER), the Gauss-Laguerre weights.
!
!    Input, real(kind=Rkind) ALPHA, the exponent of the X factor in the
!    integrand.  Set ALPHA = 0.0 for the simplest rule.
!
!    Workspace, real(kind=Rkind) B(NORDER), C(NORDER), used to hold
!    the recursion coefficients.
!
!    Output, real(kind=Rkind) CSX, the sum of the computed abscissas.
!
!    Output, real(kind=Rkind) CSA, the sum of the computed weights.
!
!    Output, real(kind=Rkind) TSX, the true sum of the abscissas.
!
!    Output, real(kind=Rkind) TSA, the true sum of the weights.
!
      integer norder
!
      real(kind=Rkind) alpha
      real(kind=Rkind) b(norder)
      real(kind=Rkind) c(norder)
      real(kind=Rkind) cc
      real(kind=Rkind) csa
      real(kind=Rkind) csx
      real(kind=Rkind) dp2
      real(kind=Rkind) gamma
      integer i
      real(kind=Rkind) p1
      real(kind=Rkind) r1
      real(kind=Rkind) r2
      real(kind=Rkind) ratio
      real(kind=Rkind) tsa
      real(kind=Rkind) tsx
      real(kind=Rkind) weight(norder)
      real(kind=Rkind) x
      real(kind=Rkind) xtab(norder)
!
!  Set the recursion coefficients.
!
      do i = 1, norder
        b(i) = ( alpha + dble ( 2 * i - 1 ) )
        c(i) = dble ( i - 1 ) * ( alpha + dble ( i - 1 ) )
      end do

      csx = 0.0
      csa = 0.0
      cc = gamma ( alpha + 1.0 )
      tsx = dble ( norder ) * ( dble ( norder ) + alpha )
      tsa = cc

      do i = 2, norder
        cc = cc * c(i)
      end do

      do i = 1, norder

        if ( i .eq. 1 ) then

          x = ( 1.0 + alpha ) * ( 3.0 + 0.92 * alpha ) /                &
            ( 1.0 + 2.4 * dble ( norder ) + 1.8 * alpha )

        else if ( i .eq. 2 ) then

          x = x + ( 15.0 + 6.25 * alpha ) /                             &
            ( 1.0 + 0.9 * alpha + 2.5 * dble ( norder ) )

        else

          r1 = ( 1.0 + 2.55 * dble ( i - 2 ) ) /                        &
            ( 1.9 * dble ( i - 2 ) )

          r2 = 1.26 * dble ( i - 2 ) * alpha /                          &
            ( 1.0 + 3.5 * dble ( i - 2 ) )

          ratio = ( r1 + r2 ) / ( 1.0 + 0.3 * alpha )

          x = x + ratio * ( x - xtab(i-2) )

        end if

        call lagroot ( x, norder, alpha, dp2, p1, b, c )

        xtab(i) = x
        weight(i) = cc / dp2 / p1

        csx = csx + xtab(i)
        csa = csa + weight(i)

      end do

      return
      end subroutine lagcom
      subroutine lagrec ( p2, dp2, p1, x, norder, alpha, b, c )
!
!***********************************************************************
!
!! LAGREC finds the value and derivative of the NORDER-th Laguerre polynomial.
!
!
!  Reference:
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!  Modified:
!
!    19 September 1998
!
!  Parameters:
!
!    Output, real(kind=Rkind) P2, the value of L(NORDER)(X).
!
!    Output, real(kind=Rkind) DP2, the value of L'(NORDER)(X).
!
!    Output, real(kind=Rkind) P1, the value of L(NORDER-1)(X).
!
!    Input, real(kind=Rkind) X, the point at which polynomials are evaluated.
!
!    Input, integer NORDER, the order of the polynomial to be computed.
!
!    Input, real(kind=Rkind) ALPHA, the exponent of the X factor in the
!    integrand.
!
!    Input, real(kind=Rkind) B(NORDER), C(NORDER), the recursion
!    coefficients.
!
      integer norder
!
      real(kind=Rkind) alpha
      real(kind=Rkind) b(norder)
      real(kind=Rkind) c(norder)
      real(kind=Rkind) dp0
      real(kind=Rkind) dp1
      real(kind=Rkind) dp2
      integer i
      real(kind=Rkind) p0
      real(kind=Rkind) p1
      real(kind=Rkind) p2
      real(kind=Rkind) x
!
      p1 = 1.0
      dp1 = 0.0

      p2 = x - alpha - 1.0
      dp2 = 1.0

      do i = 2, norder

        p0 = p1
        dp0 = dp1

        p1 = p2
        dp1 = dp2

        p2 = ( x - b(i) ) * p1 - c(i) * p0
        dp2 = ( x - b(i) ) * dp1 + p1 - c(i) * dp0

      end do

      return
      end subroutine lagrec
      subroutine lagroot ( x, norder, alpha, dp2, p1, b, c )
!
!***********************************************************************
!
!! LAGROOT improves an approximate root of a Laguerre polynomial.
!
!
!  Reference:
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!  Modified:
!
!    19 September 1998
!
!  Parameters:
!
!    Input/output, real(kind=Rkind) X, the approximate root, which
!    should be improved on output.
!
!    Input, integer NORDER, the order of the polynomial to be computed.
!
!    Input, real(kind=Rkind) ALPHA, the exponent of the X factor in the
!    integrand.
!
!    Output, real(kind=Rkind) DP2, the value of L'(NORDER)(X).
!
!    Output, real(kind=Rkind) P1, the value of L(NORDER-1)(X).
!
!    Input, real(kind=Rkind) B(NORDER), C(NORDER), the recursion
!    coefficients.
!
      real(kind=Rkind) eps
      parameter ( eps = 1.0D-12 )
!
      integer norder
!
      real(kind=Rkind) alpha
      real(kind=Rkind) b(norder)
      real(kind=Rkind) c(norder)
      real(kind=Rkind) d
      real(kind=Rkind) dp2
      integer i
      real(kind=Rkind) p1
      real(kind=Rkind) p2
      real(kind=Rkind) x
!
      do i = 1, 10

        call lagrec ( p2, dp2, p1, x, norder, alpha, b, c )

        d = p2 / dp2
        x = x - d

        if ( abs ( d ) .le. eps * ( abs ( x ) + 1.0 ) ) then
          return
        end if

      end do

      return
      end subroutine lagroot
      subroutine lagset ( norder, xtab, weight )
!
!***********************************************************************
!
!! LAGSET sets abscissas and weights for Laguerre quadrature.
!
!
!  Integration interval:
!
!    [ 0, +Infinity )
!
!  Weight function:
!
!    EXP ( - X ).
!
!  Integral to approximate:
!
!    INTEGRAL ( 0 <= X < +INFINITY ) EXP ( - X ) * F(X) dX.
!
!  Approximate integral:
!
!    SUM ( I = 1 to NORDER ) WEIGHT(I) * F ( XTAB(I) ).
!
!  Note:
!
!    The abscissas are the zeroes of the Laguerre polynomial L(NORDER)(X).
!
!  Reference:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964.
!
!    Vladimir Krylov,
!    Approximate Calculation of Integrals,
!    MacMillan, 1962.
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Modified:
!
!    17 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!    NORDER must be between 1 and 20.
!
!    Output, real(kind=Rkind) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real(kind=Rkind) WEIGHT(NORDER), the weights of the rule.
!    The weights are positive, and should add to 1.
!
      integer norder
!
      real(kind=Rkind) weight(norder)
      real(kind=Rkind) xtab(norder)
!
      if ( norder .eq. 1 ) then

        xtab(1) =    1.0 D+00
        weight(1) =  1.0 D+00

      else if ( norder .eq. 2 ) then

        xtab(1) =    0.585786437626904951198311275790 D+00
        xtab(2) =    0.341421356237309504880168872421 D+01

        weight(1) =  0.853553390593273762200422181052 D+00
        weight(2) =  0.146446609406726237799577818948 D+00

      else if ( norder .eq. 3 ) then

        xtab(1) =    0.415774556783479083311533873128 D+00
        xtab(2) =    0.229428036027904171982205036136 D+01
        xtab(3) =    0.628994508293747919686641576551 D+01

        weight(1) =  0.711093009929173015449590191143 D+00
        weight(2) =  0.278517733569240848801444888457 D+00
        weight(3) =  0.103892565015861357489649204007 D-01

      else if ( norder .eq. 4 ) then

        xtab(1) =    0.322547689619392311800361943361 D+00
        xtab(2) =    0.174576110115834657568681671252 D+01
        xtab(3) =    0.453662029692112798327928538496 D+01
        xtab(4) =    0.939507091230113312923353644342 D+01

        weight(1) =  0.603154104341633601635966023818 D+00
        weight(2) =  0.357418692437799686641492017458 D+00
        weight(3) =  0.388879085150053842724381681562 D-01
        weight(4) =  0.539294705561327450103790567621 D-03

      else if ( norder .eq. 5 ) then

        xtab(1) =    0.263560319718140910203061943361 D+00
        xtab(2) =    0.141340305910651679221840798019 D+01
        xtab(3) =    0.359642577104072208122318658878 D+01
        xtab(4) =    0.708581000585883755692212418111 D+01
        xtab(5) =    0.126408008442757826594332193066 D+02

        weight(1) =  0.521755610582808652475860928792 D+00
        weight(2) =  0.398666811083175927454133348144 D+00
        weight(3) =  0.759424496817075953876533114055 D-01
        weight(4) =  0.361175867992204845446126257304 D-02
        weight(5) =  0.233699723857762278911490845516 D-04

      else if ( norder .eq. 6 ) then

        xtab(1) =    0.222846604179260689464354826787 D+00
        xtab(2) =    0.118893210167262303074315092194 D+01
        xtab(3) =    0.299273632605931407769132528451 D+01
        xtab(4) =    0.577514356910451050183983036943 D+01
        xtab(5) =    0.983746741838258991771554702994 D+01
        xtab(6) =    0.159828739806017017825457915674 D+02

        weight(1) =  0.458964673949963593568284877709 D+00
        weight(2) =  0.417000830772120994113377566193 D+00
        weight(3) =  0.113373382074044975738706185098 D+00
        weight(4) =  0.103991974531490748989133028469 D-01
        weight(5) =  0.261017202814932059479242860001 D-03
        weight(6) =  0.898547906429621238825292052825 D-06

      else if ( norder .eq. 7 ) then

        xtab(1) =    0.193043676560362413838247885004 D+00
        xtab(2) =    0.102666489533919195034519944317 D+01
        xtab(3) =    0.256787674495074620690778622666 D+01
        xtab(4) =    0.490035308452648456810171437810 D+01
        xtab(5) =    0.818215344456286079108182755123 D+01
        xtab(6) =    0.127341802917978137580126424582 D+02
        xtab(7) =    0.193957278622625403117125820576 D+02

        weight(1) =  0.409318951701273902130432880018 D+00
        weight(2) =  0.421831277861719779929281005417 D+00
        weight(3) =  0.147126348657505278395374184637 D+00
        weight(4) =  0.206335144687169398657056149642 D-01
        weight(5) =  0.107401014328074552213195962843 D-02
        weight(6) =  0.158654643485642012687326223234 D-04
        weight(7) =  0.317031547899558056227132215385 D-07

      else if ( norder .eq. 8 ) then

        xtab(1) =    0.170279632305100999788861856608 D+00
        xtab(2) =    0.903701776799379912186020223555 D+00
        xtab(3) =    0.225108662986613068930711836697 D+01
        xtab(4) =    0.426670017028765879364942182690 D+01
        xtab(5) =    0.704590540239346569727932548212 D+01
        xtab(6) =    0.107585160101809952240599567880 D+02
        xtab(7) =    0.157406786412780045780287611584 D+02
        xtab(8) =    0.228631317368892641057005342974 D+02

        weight(1) =  0.369188589341637529920582839376 D+00
        weight(2) =  0.418786780814342956076978581333 D+00
        weight(3) =  0.175794986637171805699659866777 D+00
        weight(4) =  0.333434922612156515221325349344 D-01
        weight(5) =  0.279453623522567252493892414793 D-02
        weight(6) =  0.907650877335821310423850149336 D-04
        weight(7) =  0.848574671627253154486801830893 D-06
        weight(8) =  0.104800117487151038161508853552 D-08

      else if ( norder .eq. 9 ) then

        xtab(1) =    0.152322227731808247428107073127 D+00
        xtab(2) =    0.807220022742255847741419210952 D+00
        xtab(3) =    0.200513515561934712298303324701 D+01
        xtab(4) =    0.378347397333123299167540609364 D+01
        xtab(5) =    0.620495677787661260697353521006 D+01
        xtab(6) =    0.937298525168757620180971073215 D+01
        xtab(7) =    0.134662369110920935710978818397 D+02
        xtab(8) =    0.188335977889916966141498992996 D+02
        xtab(9) =    0.263740718909273767961410072937 D+02

        weight(1) =  0.336126421797962519673467717606 D+00
        weight(2) =  0.411213980423984387309146942793 D+00
        weight(3) =  0.199287525370885580860575607212 D+00
        weight(4) =  0.474605627656515992621163600479 D-01
        weight(5) =  0.559962661079458317700419900556 D-02
        weight(6) =  0.305249767093210566305412824291 D-03
        weight(7) =  0.659212302607535239225572284875 D-05
        weight(8) =  0.411076933034954844290241040330 D-07
        weight(9) =  0.329087403035070757646681380323 D-10

      else if ( norder .eq. 10 ) then

        xtab(1) =    0.137793470540492430830772505653 D+00
        xtab(2) =    0.729454549503170498160373121676 D+00
        xtab(3) =    0.180834290174031604823292007575 D+01
        xtab(4) =    0.340143369785489951448253222141 D+01
        xtab(5) =    0.555249614006380363241755848687 D+01
        xtab(6) =    0.833015274676449670023876719727 D+01
        xtab(7) =    0.118437858379000655649185389191 D+02
        xtab(8) =    0.162792578313781020995326539358 D+02
        xtab(9) =    0.219965858119807619512770901956 D+02
        xtab(10) =   0.299206970122738915599087933408 D+02

        weight(1) =  0.308441115765020141547470834678 D+00
        weight(2) =  0.401119929155273551515780309913 D+00
        weight(3) =  0.218068287611809421588648523475 D+00
        weight(4) =  0.620874560986777473929021293135 D-01
        weight(5) =  0.950151697518110055383907219417 D-02
        weight(6) =  0.753008388587538775455964353676 D-03
        weight(7) =  0.282592334959956556742256382685 D-04
        weight(8) =  0.424931398496268637258657665975 D-06
        weight(9) =  0.183956482397963078092153522436 D-08
        weight(10) = 0.991182721960900855837754728324 D-12

      else if ( norder .eq. 11 ) then

        xtab(1) =    0.125796442187967522675794577516 D+00
        xtab(2) =    0.665418255839227841678127839420 D+00
        xtab(3) =    0.164715054587216930958700321365 D+01
        xtab(4) =    0.309113814303525495330195934259 D+01
        xtab(5) =    0.502928440157983321236999508366 D+01
        xtab(6) =    0.750988786380661681941099714450 D+01
        xtab(7) =    0.106059509995469677805559216457 D+02
        xtab(8) =    0.144316137580641855353200450349 D+02
        xtab(9) =    0.191788574032146786478174853989 D+02
        xtab(10) =   0.252177093396775611040909447797 D+02
        xtab(11) =   0.334971928471755372731917259395 D+02

        weight(1) =  0.284933212894200605056051024724 D+00
        weight(2) =  0.389720889527849377937553508048 D+00
        weight(3) =  0.232781831848991333940223795543 D+00
        weight(4) =  0.765644535461966864008541790132 D-01
        weight(5) =  0.143932827673506950918639187409 D-01
        weight(6) =  0.151888084648487306984777640042 D-02
        weight(7) =  0.851312243547192259720424170600 D-04
        weight(8) =  0.229240387957450407857683270709 D-05
        weight(9) =  0.248635370276779587373391491114 D-07
        weight(10) = 0.771262693369132047028152590222 D-10
        weight(11) = 0.288377586832362386159777761217 D-13

      else if ( norder .eq. 12 ) then

        xtab(1) =    0.115722117358020675267196428240 D+00
        xtab(2) =    0.611757484515130665391630053042 D+00
        xtab(3) =    0.151261026977641878678173792687 D+01
        xtab(4) =    0.283375133774350722862747177657 D+01
        xtab(5) =    0.459922763941834848460572922485 D+01
        xtab(6) =    0.684452545311517734775433041849 D+01
        xtab(7) =    0.962131684245686704391238234923 D+01
        xtab(8) =    0.130060549933063477203460524294 D+02
        xtab(9) =    0.171168551874622557281840528008 D+02
        xtab(10) =   0.221510903793970056699218950837 D+02
        xtab(11) =   0.284879672509840003125686072325 D+02
        xtab(12) =   0.370991210444669203366389142764 D+02

        weight(1) =  0.264731371055443190349738892056 D+00
        weight(2) =  0.377759275873137982024490556707 D+00
        weight(3) =  0.244082011319877564254870818274 D+00
        weight(4) =  0.904492222116809307275054934667 D-01
        weight(5) =  0.201023811546340965226612867827 D-01
        weight(6) =  0.266397354186531588105415760678 D-02
        weight(7) =  0.203231592662999392121432860438 D-03
        weight(8) =  0.836505585681979874533632766396 D-05
        weight(9) =  0.166849387654091026116989532619 D-06
        weight(10) = 0.134239103051500414552392025055 D-08
        weight(11) = 0.306160163503502078142407718971 D-11
        weight(12) = 0.814807746742624168247311868103 D-15

      else if ( norder .eq. 13 ) then

        xtab(1) =    0.107142388472252310648493376977 D+00
        xtab(2) =    0.566131899040401853406036347177 D+00
        xtab(3) =    0.139856433645101971792750259921 D+01
        xtab(4) =    0.261659710840641129808364008472 D+01
        xtab(5) =    0.423884592901703327937303389926 D+01
        xtab(6) =    0.629225627114007378039376523025 D+01
        xtab(7) =    0.881500194118697804733348868036 D+01
        xtab(8) =    0.118614035888112425762212021880 D+02
        xtab(9) =    0.155107620377037527818478532958 D+02
        xtab(10) =   0.198846356638802283332036594634 D+02
        xtab(11) =   0.251852638646777580842970297823 D+02
        xtab(12) =   0.318003863019472683713663283526 D+02
        xtab(13) =   0.407230086692655795658979667001 D+02

        weight(1) =  0.247188708429962621346249185964 D+00
        weight(2) =  0.365688822900521945306717530893 D+00
        weight(3) =  0.252562420057658502356824288815 D+00
        weight(4) =  0.103470758024183705114218631672 D+00
        weight(5) =  0.264327544155616157781587735702 D-01
        weight(6) =  0.422039604025475276555209292644 D-02
        weight(7) =  0.411881770472734774892472527082 D-03
        weight(8) =  0.235154739815532386882897300772 D-04
        weight(9) =  0.731731162024909910401047197761 D-06
        weight(10) = 0.110884162570398067979150974759 D-07
        weight(11) = 0.677082669220589884064621459082 D-10
        weight(12) = 0.115997995990507606094507145382 D-12
        weight(13) = 0.224509320389275841599187226865 D-16

      else if ( norder .eq. 14 ) then

        xtab(1) =    0.997475070325975745736829452514 D-01
        xtab(2) =    0.526857648851902896404583451502 D+00
        xtab(3) =    0.130062912125149648170842022116 D+01
        xtab(4) =    0.243080107873084463616999751038 D+01
        xtab(5) =    0.393210282229321888213134366778 D+01
        xtab(6) =    0.582553621830170841933899983898 D+01
        xtab(7) =    0.814024014156514503005978046052 D+01
        xtab(8) =    0.109164995073660188408130510904 D+02
        xtab(9) =    0.142108050111612886831059780825 D+02
        xtab(10) =   0.181048922202180984125546272083 D+02
        xtab(11) =   0.227233816282696248232280886985 D+02
        xtab(12) =   0.282729817232482056954158923218 D+02
        xtab(13) =   0.351494436605924265828643121364 D+02
        xtab(14) =   0.443660817111174230416312423666 D+02

        weight(1) =  0.231815577144864977840774861104 D+00
        weight(2) =  0.353784691597543151802331301273 D+00
        weight(3) =  0.258734610245428085987320561144 D+00
        weight(4) =  0.115482893556923210087304988673 D+00
        weight(5) =  0.331920921593373600387499587137 D-01
        weight(6) =  0.619286943700661021678785967675 D-02
        weight(7) =  0.739890377867385942425890907080 D-03
        weight(8) =  0.549071946684169837857331777667 D-04
        weight(9) =  0.240958576408537749675775256553 D-05
        weight(10) = 0.580154398167649518088619303904 D-07
        weight(11) = 0.681931469248497411961562387084 D-09
        weight(12) = 0.322120775189484793980885399656 D-11
        weight(13) = 0.422135244051658735159797335643 D-14
        weight(14) = 0.605237502228918880839870806281 D-18

      else if ( norder .eq. 15 ) then

        xtab(1) =    0.933078120172818047629030383672 D-01
        xtab(2) =    0.492691740301883908960101791412 D+00
        xtab(3) =    0.121559541207094946372992716488 D+01
        xtab(4) =    0.226994952620374320247421741375 D+01
        xtab(5) =    0.366762272175143727724905959436 D+01
        xtab(6) =    0.542533662741355316534358132596 D+01
        xtab(7) =    0.756591622661306786049739555812 D+01
        xtab(8) =    0.101202285680191127347927394568 D+02
        xtab(9) =    0.131302824821757235640991204176 D+02
        xtab(10) =   0.166544077083299578225202408430 D+02
        xtab(11) =   0.207764788994487667729157175676 D+02
        xtab(12) =   0.256238942267287801445868285977 D+02
        xtab(13) =   0.314075191697539385152432196202 D+02
        xtab(14) =   0.385306833064860094162515167595 D+02
        xtab(15) =   0.480260855726857943465734308508 D+02

        weight(1) =  0.218234885940086889856413236448 D+00
        weight(2) =  0.342210177922883329638948956807 D+00
        weight(3) =  0.263027577941680097414812275022 D+00
        weight(4) =  0.126425818105930535843030549378 D+00
        weight(5) =  0.402068649210009148415854789871 D-01
        weight(6) =  0.856387780361183836391575987649 D-02
        weight(7) =  0.121243614721425207621920522467 D-02
        weight(8) =  0.111674392344251941992578595518 D-03
        weight(9) =  0.645992676202290092465319025312 D-05
        weight(10) = 0.222631690709627263033182809179 D-06
        weight(11) = 0.422743038497936500735127949331 D-08
        weight(12) = 0.392189726704108929038460981949 D-10
        weight(13) = 0.145651526407312640633273963455 D-12
        weight(14) = 0.148302705111330133546164737187 D-15
        weight(15) = 0.160059490621113323104997812370 D-19

      else if ( norder .eq. 16 ) then

        xtab(1) =    0.876494104789278403601980973401 D-01
        xtab(2) =    0.462696328915080831880838260664 D+00
        xtab(3) =    0.114105777483122685687794501811 D+01
        xtab(4) =    0.212928364509838061632615907066 D+01
        xtab(5) =    0.343708663389320664523510701675 D+01
        xtab(6) =    0.507801861454976791292305830814 D+01
        xtab(7) =    0.707033853504823413039598947080 D+01
        xtab(8) =    0.943831433639193878394724672911 D+01
        xtab(9) =    0.122142233688661587369391246088 D+02
        xtab(10) =   0.154415273687816170767647741622 D+02
        xtab(11) =   0.191801568567531348546631409497 D+02
        xtab(12) =   0.235159056939919085318231872752 D+02
        xtab(13) =   0.285787297428821403675206137099 D+02
        xtab(14) =   0.345833987022866258145276871778 D+02
        xtab(15) =   0.419404526476883326354722330252 D+02
        xtab(16) =   0.517011603395433183643426971197 D+02

        weight(1) =  0.206151714957800994334273636741 D+00
        weight(2) =  0.331057854950884165992983098710 D+00
        weight(3) =  0.265795777644214152599502020650 D+00
        weight(4) =  0.136296934296377539975547513526 D+00
        weight(5) =  0.473289286941252189780623392781 D-01
        weight(6) =  0.112999000803394532312490459701 D-01
        weight(7) =  0.184907094352631086429176783252 D-02
        weight(8) =  0.204271915308278460126018338421 D-03
        weight(9) =  0.148445868739812987713515067551 D-04
        weight(10) = 0.682831933087119956439559590327 D-06
        weight(11) = 0.188102484107967321388159920418 D-07
        weight(12) = 0.286235024297388161963062629156 D-09
        weight(13) = 0.212707903322410296739033610978 D-11
        weight(14) = 0.629796700251786778717446214552 D-14
        weight(15) = 0.505047370003551282040213233303 D-17
        weight(16) = 0.416146237037285519042648356116 D-21

      else if ( norder .eq. 17 ) then

        xtab(1) =    0.826382147089476690543986151980 D-01
        xtab(2) =    0.436150323558710436375959029847 D+00
        xtab(3) =    0.107517657751142857732980316755 D+01
        xtab(4) =    0.200519353164923224070293371933 D+01
        xtab(5) =    0.323425612404744376157380120696 D+01
        xtab(6) =    0.477351351370019726480932076262 D+01
        xtab(7) =    0.663782920536495266541643929703 D+01
        xtab(8) =    0.884668551116980005369470571184 D+01
        xtab(9) =    0.114255293193733525869726151469 D+02
        xtab(10) =   0.144078230374813180021982874959 D+02
        xtab(11) =   0.178382847307011409290658752412 D+02
        xtab(12) =   0.217782682577222653261749080522 D+02
        xtab(13) =   0.263153178112487997766149598369 D+02
        xtab(14) =   0.315817716804567331343908517497 D+02
        xtab(15) =   0.377960938374771007286092846663 D+02
        xtab(16) =   0.453757165339889661829258363215 D+02
        xtab(17) =   0.553897517898396106640900199790 D+02

        weight(1) =  0.195332205251770832145927297697 D+00
        weight(2) =  0.320375357274540281336625631970 D+00
        weight(3) =  0.267329726357171097238809604160 D+00
        weight(4) =  0.145129854358758625407426447473 D+00
        weight(5) =  0.544369432453384577793805803066 D-01
        weight(6) =  0.143572977660618672917767247431 D-01
        weight(7) =  0.266282473557277256843236250006 D-02
        weight(8) =  0.343679727156299920611775097985 D-03
        weight(9) =  0.302755178378287010943703641131 D-04
        weight(10) = 0.176851505323167689538081156159 D-05
        weight(11) = 0.657627288681043332199222748162 D-07
        weight(12) = 0.146973093215954679034375821888 D-08
        weight(13) = 0.181691036255544979555476861323 D-10
        weight(14) = 0.109540138892868740297645078918 D-12
        weight(15) = 0.261737388222337042155132062413 D-15
        weight(16) = 0.167293569314615469085022374652 D-18
        weight(17) = 0.106562631627404278815253271162 D-22

      else if ( norder .eq. 18 ) then

        xtab(1) =    0.781691666697054712986747615334 D-01
        xtab(2) =    0.412490085259129291039101536536 D+00
        xtab(3) =    0.101652017962353968919093686187 D+01
        xtab(4) =    0.189488850996976091426727831954 D+01
        xtab(5) =    0.305435311320265975115241130719 D+01
        xtab(6) =    0.450420553888989282633795571455 D+01
        xtab(7) =    0.625672507394911145274209116326 D+01
        xtab(8) =    0.832782515660563002170470261564 D+01
        xtab(9) =    0.107379900477576093352179033397 D+02
        xtab(10) =   0.135136562075550898190863812108 D+02
        xtab(11) =   0.166893062819301059378183984163 D+02
        xtab(12) =   0.203107676262677428561313764553 D+02
        xtab(13) =   0.244406813592837027656442257980 D+02
        xtab(14) =   0.291682086625796161312980677805 D+02
        xtab(15) =   0.346279270656601721454012429438 D+02
        xtab(16) =   0.410418167728087581392948614284 D+02
        xtab(17) =   0.488339227160865227486586093290 D+02
        xtab(18) =   0.590905464359012507037157810181 D+02

        weight(1) =  0.185588603146918805623337752284 D+00
        weight(2) =  0.310181766370225293649597595713 D+00
        weight(3) =  0.267866567148536354820854394783 D+00
        weight(4) =  0.152979747468074906553843082053 D+00
        weight(5) =  0.614349178609616527076780103487 D-01
        weight(6) =  0.176872130807729312772600233761 D-01
        weight(7) =  0.366017976775991779802657207890 D-02
        weight(8) =  0.540622787007735323128416319257 D-03
        weight(9) =  0.561696505121423113817929049294 D-04
        weight(10) = 0.401530788370115755858883625279 D-05
        weight(11) = 0.191466985667567497969210011321 D-06
        weight(12) = 0.583609526863159412918086289717 D-08
        weight(13) = 0.107171126695539012772851317562 D-09
        weight(14) = 0.108909871388883385562011298291 D-11
        weight(15) = 0.538666474837830887608094323164 D-14
        weight(16) = 0.104986597803570340877859934846 D-16
        weight(17) = 0.540539845163105364356554467358 D-20
        weight(18) = 0.269165326920102862708377715980 D-24

      else if ( norder .eq. 19 ) then

        xtab(1) =    0.741587837572050877131369916024 D-01
        xtab(2) =    0.391268613319994607337648350299 D+00
        xtab(3) =    0.963957343997958058624879377130 D+00
        xtab(4) =    0.179617558206832812557725825252 D+01
        xtab(5) =    0.289365138187378399116494713237 D+01
        xtab(6) =    0.426421553962776647436040018167 D+01
        xtab(7) =    0.591814156164404855815360191408 D+01
        xtab(8) =    0.786861891533473373105668358176 D+01
        xtab(9) =    0.101324237168152659251627415800 D+02
        xtab(10) =   0.127308814638423980045092979656 D+02
        xtab(11) =   0.156912783398358885454136069861 D+02
        xtab(12) =   0.190489932098235501532136429732 D+02
        xtab(13) =   0.228508497608294829323930586693 D+02
        xtab(14) =   0.271606693274114488789963947149 D+02
        xtab(15) =   0.320691222518622423224362865906 D+02
        xtab(16) =   0.377129058012196494770647508283 D+02
        xtab(17) =   0.443173627958314961196067736013 D+02
        xtab(18) =   0.523129024574043831658644222420 D+02
        xtab(19) =   0.628024231535003758413504690673 D+02

        weight(1) =  0.176768474915912502251035479815 D+00
        weight(2) =  0.300478143607254379482156807712 D+00
        weight(3) =  0.267599547038175030772695440648 D+00
        weight(4) =  0.159913372135580216785512147895 D+00
        weight(5) =  0.682493799761491134552355368344 D-01
        weight(6) =  0.212393076065443249244062193091 D-01
        weight(7) =  0.484162735114839596725013121019 D-02
        weight(8) =  0.804912747381366766594647138204 D-03
        weight(9) =  0.965247209315350170843161738801 D-04
        weight(10) = 0.820730525805103054408982992869 D-05
        weight(11) = 0.483056672473077253944806671560 D-06
        weight(12) = 0.190499136112328569993615674552 D-07
        weight(13) = 0.481668463092806155766936380273 D-09
        weight(14) = 0.734825883955114437684376840171 D-11
        weight(15) = 0.620227538757261639893719012423 D-13
        weight(16) = 0.254143084301542272371866857954 D-15
        weight(17) = 0.407886129682571235007187465134 D-18
        weight(18) = 0.170775018759383706100412325084 D-21
        weight(19) = 0.671506464990818995998969111749 D-26

      else if ( norder .eq. 20 ) then

        xtab(1) =    0.705398896919887533666890045842 D-01
        xtab(2) =    0.372126818001611443794241388761 D+00
        xtab(3) =    0.916582102483273564667716277074 D+00
        xtab(4) =    0.170730653102834388068768966741 D+01
        xtab(5) =    0.274919925530943212964503046049 D+01
        xtab(6) =    0.404892531385088692237495336913 D+01
        xtab(7) =    0.561517497086161651410453988565 D+01
        xtab(8) =    0.745901745367106330976886021837 D+01
        xtab(9) =    0.959439286958109677247367273428 D+01
        xtab(10) =   0.120388025469643163096234092989 D+02
        xtab(11) =   0.148142934426307399785126797100 D+02
        xtab(12) =   0.179488955205193760173657909926 D+02
        xtab(13) =   0.214787882402850109757351703696 D+02
        xtab(14) =   0.254517027931869055035186774846 D+02
        xtab(15) =   0.299325546317006120067136561352 D+02
        xtab(16) =   0.350134342404790000062849359067 D+02
        xtab(17) =   0.408330570567285710620295677078 D+02
        xtab(18) =   0.476199940473465021399416271529 D+02
        xtab(19) =   0.558107957500638988907507734445 D+02
        xtab(20) =   0.665244165256157538186403187915 D+02

        weight(1) =  0.168746801851113862149223899689 D+00
        weight(2) =  0.291254362006068281716795323812 D+00
        weight(3) =  0.266686102867001288549520868998 D+00
        weight(4) =  0.166002453269506840031469127816 D+00
        weight(5) =  0.748260646687923705400624639615 D-01
        weight(6) =  0.249644173092832210728227383234 D-01
        weight(7) =  0.620255084457223684744754785395 D-02
        weight(8) =  0.114496238647690824203955356969 D-02
        weight(9) =  0.155741773027811974779809513214 D-03
        weight(10) = 0.154014408652249156893806714048 D-04
        weight(11) = 0.108648636651798235147970004439 D-05
        weight(12) = 0.533012090955671475092780244305 D-07
        weight(13) = 0.175798117905058200357787637840 D-08
        weight(14) = 0.372550240251232087262924585338 D-10
        weight(15) = 0.476752925157819052449488071613 D-12
        weight(16) = 0.337284424336243841236506064991 D-14
        weight(17) = 0.115501433950039883096396247181 D-16
        weight(18) = 0.153952214058234355346383319667 D-19
        weight(19) = 0.528644272556915782880273587683 D-23
        weight(20) = 0.165645661249902329590781908529 D-27

      else

        write ( *, * ) ' '
        write ( *, * ) 'LAGSET - Fatal error!'
        write ( *, * ) '  Illegal value of NORDER = ', norder
        write ( *, * ) '  Legal values are 1 to 20.'
        stop

      end if

      return
      end subroutine lagset
      subroutine lagsum ( func, a, norder, xtab, weight, result )
!
!***********************************************************************
!
!! LAGSUM carries out Laguerre quadrature over [ A, +Infinity ).
!
!
!  Integration interval:
!
!    [ A, +Infinity ).
!
!  Weight function:
!
!    EXP ( - X ).
!
!  Integral to approximate:
!
!    INTEGRAL ( A <= X <= +Infinity ) EXP ( -X ) * F(X) dX.
!
!  Approximate integral:
!
!    EXP ( - A ) * SUM ( I = 1 to NORDER ) WEIGHT(I) * F ( XTAB(I) + A ).
!
!  Reference:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964.
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Modified:
!
!    15 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FUNC, the name of the FORTRAN function which
!    evaluates the integrand.  The function must have the form
!      real(kind=Rkind) func ( x ).
!
!    Input, real(kind=Rkind) A, the beginning of the integration interval.
!
!    Input, integer NORDER, the order of the rule.
!
!    Input, real(kind=Rkind) XTAB(NORDER), the abscissas of the rule.
!
!    Input, real(kind=Rkind) WEIGHT(NORDER), the weights of the rule.
!
!    Output, real(kind=Rkind) RESULT, the approximate value of the integral.
!
      integer norder
!
      real(kind=Rkind) a
      real(kind=Rkind) func
      integer i
      real(kind=Rkind) result
      real(kind=Rkind) xtab(norder)
      real(kind=Rkind) weight(norder)
!
      external func
!
      if ( norder .lt. 1 ) then
        write ( *, * ) ' '
        write ( *, * ) 'LAGSUM - Fatal error!'
        write ( *, * ) '  Nonpositive NORDER = ', norder
        stop
      end if

      result = 0.0 D+00
      do i = 1, norder
        result = result + weight(i) * func ( xtab(i) + a )
      end do
      result = exp ( - a ) * result

      return
      end subroutine lagsum
      subroutine legcom ( norder, xtab, weight )
!
!***********************************************************************
!
!! LEGCOM computes abscissas and weights for Gauss-Legendre quadrature.
!
!
!  Integration interval:
!
!    [ -1, 1 ]
!
!  Weight function:
!
!    1.
!
!  Integral to approximate:
!
!    INTEGRAL ( -1 <= X <= 1 ) F(X) dX.
!
!  Approximate integral:
!
!    SUM ( I = 1 to NORDER ) WEIGHT(I) * F ( XTAB(I) ).
!
!  Note:
!
!    LEGCOM computes the values using Newton iteration and real(kind=Rkind).
!
!  Modified:
!
!    16 September 1998
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!    NORDER must be greater than 0.
!
!    Output, real(kind=Rkind) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real(kind=Rkind) WEIGHT(NORDER), the weights of the rule.
!    The weights are positive, symmetric, and should sum to 2.
!
      real(kind=Rkind) pi
      parameter ( pi = 3.141592653589793238462643383279D+00 )
!
      integer norder
!
      real(kind=Rkind) d1
      real(kind=Rkind) d2pn
      real(kind=Rkind) d3pn
      real(kind=Rkind) d4pn
      real(kind=Rkind) dp
      real(kind=Rkind) dpn
      real(kind=Rkind) e1
      real(kind=Rkind) fx
      real(kind=Rkind) h
      integer i
      integer iback
      integer k
      integer m
      integer mp1mi
      integer ncopy
      integer nmove
      real(kind=Rkind) p
      real(kind=Rkind) pk
      real(kind=Rkind) pkm1
      real(kind=Rkind) pkp1
      real(kind=Rkind) t
      real(kind=Rkind) u
      real(kind=Rkind) v
      real(kind=Rkind) x0
      real(kind=Rkind) xtab(norder)
      real(kind=Rkind) xtemp
      real(kind=Rkind) weight(norder)
!
      if ( norder .lt. 1 ) then
        write ( *, * ) ' '
        write ( *, * ) 'LEGCOM - Fatal error!'
        write ( *, * ) '  Illegal value of NORDER = ', norder
        stop
      end if

      e1 = dble ( norder * ( norder + 1 ) )

      m = ( norder + 1 ) / 2

      do i = 1, ( norder + 1 ) / 2

        mp1mi = m + 1 - i
        t = dble ( 4 * i - 1 ) * pi / dble ( 4 * norder + 2 )
        x0 = cos(t) * ( 1.0D+00 - ( 1.0D+00 - 1.0D+00 /                 &
          dble ( norder ) ) / dble ( 8 * norder**2) )

        pkm1 = 1.0D+00
        pk = x0

        do k = 2, norder
          pkp1 = 2.0D+00 * x0 * pk - pkm1 - ( x0 * pk - pkm1 )          &
            / dble ( k )
          pkm1 = pk
          pk = pkp1
        end do

        d1 = dble ( norder ) * ( pkm1 - x0 * pk )

        dpn = d1 / ( 1.0D+00 - x0**2 )

        d2pn = ( 2.0D+00 * x0 * dpn - e1 * pk ) / ( 1.0D+00 - x0**2 )

        d3pn = ( 4.0D+00 * x0 * d2pn + ( 2.0D+00 - e1 ) * dpn ) /       &
          ( 1.0D+00 - x0**2 )

        d4pn = ( 6.0D+00 * x0 * d3pn + ( 6.0D+00 - e1 ) * d2pn ) /      &
          ( 1.0D+00 - x0**2 )

        u = pk / dpn
        v = d2pn / dpn
!
!  Initial approximation H:
!
        h = - u * ( 1.0D+00 + 0.5D+00 * u * ( v + u * ( v**2 - d3pn     &
          / ( 3.0D+00 * dpn ) ) ) )
!
!  Refine H using one step of Newton's method:
!
        p = pk + h * ( dpn + 0.5D+00 * h * ( d2pn + h / 3.0D+00         &
          * ( d3pn + 0.25D+00 * h * d4pn ) ) )

        dp = dpn + h * ( d2pn + 0.5D+00 * h *                           &
          ( d3pn + h * d4pn / 3.0D+00 ) )

        h = h - p / dp

        xtemp = x0 + h

        xtab(mp1mi) = xtemp

        fx = d1 - h * e1 * ( pk + 0.5D+00 * h * ( dpn + h / 3.0D+00     &
          * ( d2pn + 0.25D+00 * h * ( d3pn + 0.2D+00 * h * d4pn ) ) ) )

        weight(mp1mi) = 2.0D+00 * ( 1.0D+00 - xtemp**2 ) / fx**2

      end do

      if ( mod ( norder, 2 ) .eq. 1 ) then
        xtab(1) = 0.0D+00
      end if
!
!  Shift the data up.
!
      nmove = ( norder + 1 ) / 2
      ncopy = norder - nmove

      do i = 1, nmove
        iback = norder + 1 - i
        xtab(iback) = xtab(iback-ncopy)
        weight(iback) = weight(iback-ncopy)
      end do
!
!  Reflect values for the negative abscissas.
!
      do i = 1, norder - nmove
        xtab(i) = - xtab(norder+1-i)
        weight(i) = weight(norder+1-i)
      end do

      return
      end subroutine legcom
      subroutine legrec ( p2, dp2, p1, x, norder )
!
!***********************************************************************
!
!! LEGREC finds the value and derivative of the NORDER-th Legendre polynomial.
!
!
!  Reference:
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!  Modified:
!
!    19 September 1998
!
!  Parameters:
!
!    Output, real(kind=Rkind) P2, the value of P(NORDER)(X).
!
!    Output, real(kind=Rkind) DP2, the value of P'(NORDER)(X).
!
!    Output, real(kind=Rkind) P1, the value of P(NORDER-1)(X).
!
!    Input, real(kind=Rkind) X, the point at which polynomials are evaluated.
!
!    Input, integer NORDER, the order of the polynomial to be computed.
!
      integer norder
!
      real(kind=Rkind) dp0
      real(kind=Rkind) dp1
      real(kind=Rkind) dp2
      integer i
      real(kind=Rkind) p0
      real(kind=Rkind) p1
      real(kind=Rkind) p2
      real(kind=Rkind) x
!
      p1 = 1.0
      dp1 = 0.0

      p2 = x
      dp2 = 1.0

      do i = 2, norder

        p0 = p1
        dp0 = dp1

        p1 = p2
        dp1 = dp2

        p2 = ( dble ( 2 * i - 1 ) * x * p1                              &
          - dble ( i - 1 ) * p0 ) / dble ( i )

        dp2 = ( dble ( 2 * i - 1 ) * ( p1 + x * dp1 )                   &
          - dble ( i - 1 ) * dp0 ) / dble ( i )

      end do

      return
      end subroutine legrec
      subroutine legset ( norder, xtab, weight )
!
!***********************************************************************
!
!! LEGSET sets abscissas and weights for Gauss-Legendre quadrature.
!
!
!  Integration interval:
!
!    [ -1, 1 ]
!
!  Weight function:
!
!    1.
!
!  Integral to approximate:
!
!    INTEGRAL ( -1 <=  X <= 1 ) F(X) dX.
!
!  Approximate integral:
!
!    SUM ( I = 1 to NORDER ) WEIGHT(I) * F ( XTAB(I) ).
!
!  Precision:
!
!    The quadrature rule will integrate exactly all polynomials up to
!    X**(2*NORDER-1).
!
!  Note:
!
!    The abscissas of the rule are the zeroes of the Legendre polynomial
!    P(NORDER)(X).
!
!    The integral produced by a Gauss-Legendre rule is equal to the
!    integral of the unique polynomial of degree NORDER-1 which
!    agrees with the function at the NORDER abscissas of the rule.
!
!  Reference:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964.
!
!    Vladimir Krylov,
!    Approximate Calculation of Integrals,
!    MacMillan, 1962.
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Modified:
!
!    20 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!    NORDER must be between 1 and 20, 32 or 64.
!
!    Output, real(kind=Rkind) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real(kind=Rkind) WEIGHT(NORDER), the weights of the rule.
!    The weights are positive, symmetric and should sum to 2.
!
      integer norder
!
      real(kind=Rkind) xtab(norder)
      real(kind=Rkind) weight(norder)
!
      if ( norder .eq. 1 ) then

        xtab(1) =   0.0 D+00

        weight(1) = 2.0 D+00

      else if ( norder .eq. 2 ) then

        xtab(1) = - 0.577350269189625764509148780502 D+00
        xtab(2) =   0.577350269189625764509148780502 D+00

        weight(1) = 1.0 D+00
        weight(2) = 1.0 D+00

      else if ( norder .eq. 3 ) then

        xtab(1) = - 0.774596669241483377035853079956 D+00
        xtab(2) =   0.0 D+00
        xtab(3) =   0.774596669241483377035853079956 D+00

        weight(1) = 5.0 D+00 / 9.0 D+00
        weight(2) = 8.0 D+00 / 9.0 D+00
        weight(3) = 5.0 D+00 / 9.0 D+00

      else if ( norder .eq. 4 ) then

        xtab(1) = - 0.861136311594052575223946488893 D+00
        xtab(2) = - 0.339981043584856264802665759103 D+00
        xtab(3) =   0.339981043584856264802665759103 D+00
        xtab(4) =   0.861136311594052575223946488893 D+00

        weight(1) = 0.347854845137453857373063949222 D+00
        weight(2) = 0.652145154862546142626936050778 D+00
        weight(3) = 0.652145154862546142626936050778 D+00
        weight(4) = 0.347854845137453857373063949222 D+00

      else if ( norder .eq. 5 ) then

        xtab(1) = - 0.906179845938663992797626878299 D+00
        xtab(2) = - 0.538469310105683091036314420700 D+00
        xtab(3) =   0.0 D+00
        xtab(4) =   0.538469310105683091036314420700 D+00
        xtab(5) =   0.906179845938663992797626878299 D+00

        weight(1) = 0.236926885056189087514264040720 D+00
        weight(2) = 0.478628670499366468041291514836 D+00
        weight(3) = 0.568888888888888888888888888889 D+00
        weight(4) = 0.478628670499366468041291514836 D+00
        weight(5) = 0.236926885056189087514264040720 D+00

      else if ( norder .eq. 6 ) then

        xtab(1) = - 0.932469514203152027812301554494 D+00
        xtab(2) = - 0.661209386466264513661399595020 D+00
        xtab(3) = - 0.238619186083196908630501721681 D+00
        xtab(4) =   0.238619186083196908630501721681 D+00
        xtab(5) =   0.661209386466264513661399595020 D+00
        xtab(6) =   0.932469514203152027812301554494 D+00

        weight(1) = 0.171324492379170345040296142173 D+00
        weight(2) = 0.360761573048138607569833513838 D+00
        weight(3) = 0.467913934572691047389870343990 D+00
        weight(4) = 0.467913934572691047389870343990 D+00
        weight(5) = 0.360761573048138607569833513838 D+00
        weight(6) = 0.171324492379170345040296142173 D+00

      else if ( norder .eq. 7 ) then

        xtab(1) = - 0.949107912342758524526189684048 D+00
        xtab(2) = - 0.741531185599394439863864773281 D+00
        xtab(3) = - 0.405845151377397166906606412077 D+00
        xtab(4) =   0.0 D+00
        xtab(5) =   0.405845151377397166906606412077 D+00
        xtab(6) =   0.741531185599394439863864773281 D+00
        xtab(7) =   0.949107912342758524526189684048 D+00

        weight(1) = 0.129484966168869693270611432679 D+00
        weight(2) = 0.279705391489276667901467771424 D+00
        weight(3) = 0.381830050505118944950369775489 D+00
        weight(4) = 0.417959183673469387755102040816 D+00
        weight(5) = 0.381830050505118944950369775489 D+00
        weight(6) = 0.279705391489276667901467771424 D+00
        weight(7) = 0.129484966168869693270611432679 D+00

      else if ( norder .eq. 8 ) then

        xtab(1) = - 0.960289856497536231683560868569 D+00
        xtab(2) = - 0.796666477413626739591553936476 D+00
        xtab(3) = - 0.525532409916328985817739049189 D+00
        xtab(4) = - 0.183434642495649804939476142360 D+00
        xtab(5) =   0.183434642495649804939476142360 D+00
        xtab(6) =   0.525532409916328985817739049189 D+00
        xtab(7) =   0.796666477413626739591553936476 D+00
        xtab(8) =   0.960289856497536231683560868569 D+00

        weight(1) = 0.101228536290376259152531354310 D+00
        weight(2) = 0.222381034453374470544355994426 D+00
        weight(3) = 0.313706645877887287337962201987 D+00
        weight(4) = 0.362683783378361982965150449277 D+00
        weight(5) = 0.362683783378361982965150449277 D+00
        weight(6) = 0.313706645877887287337962201987 D+00
        weight(7) = 0.222381034453374470544355994426 D+00
        weight(8) = 0.101228536290376259152531354310 D+00

      else if ( norder .eq. 9 ) then

        xtab(1) = - 0.968160239507626089835576202904 D+00
        xtab(2) = - 0.836031107326635794299429788070 D+00
        xtab(3) = - 0.613371432700590397308702039341 D+00
        xtab(4) = - 0.324253423403808929038538014643 D+00
        xtab(5) =   0.0 D+00
        xtab(6) =   0.324253423403808929038538014643 D+00
        xtab(7) =   0.613371432700590397308702039341 D+00
        xtab(8) =   0.836031107326635794299429788070 D+00
        xtab(9) =   0.968160239507626089835576202904 D+00

        weight(1) = 0.812743883615744119718921581105 D-01
        weight(2) = 0.180648160694857404058472031243 D+00
        weight(3) = 0.260610696402935462318742869419 D+00
        weight(4) = 0.312347077040002840068630406584 D+00
        weight(5) = 0.330239355001259763164525069287 D+00
        weight(6) = 0.312347077040002840068630406584 D+00
        weight(7) = 0.260610696402935462318742869419 D+00
        weight(8) = 0.180648160694857404058472031243 D+00
        weight(9) = 0.812743883615744119718921581105 D-01

      else if ( norder .eq. 10 ) then

        xtab(1) =  - 0.973906528517171720077964012084 D+00
        xtab(2) =  - 0.865063366688984510732096688423 D+00
        xtab(3) =  - 0.679409568299024406234327365115 D+00
        xtab(4) =  - 0.433395394129247290799265943166 D+00
        xtab(5) =  - 0.148874338981631210884826001130 D+00
        xtab(6) =    0.148874338981631210884826001130 D+00
        xtab(7) =    0.433395394129247290799265943166 D+00
        xtab(8) =    0.679409568299024406234327365115 D+00
        xtab(9) =    0.865063366688984510732096688423 D+00
        xtab(10) =   0.973906528517171720077964012084 D+00

        weight(1) =  0.666713443086881375935688098933 D-01
        weight(2) =  0.149451349150580593145776339658 D+00
        weight(3) =  0.219086362515982043995534934228 D+00
        weight(4) =  0.269266719309996355091226921569 D+00
        weight(5) =  0.295524224714752870173892994651 D+00
        weight(6) =  0.295524224714752870173892994651 D+00
        weight(7) =  0.269266719309996355091226921569 D+00
        weight(8) =  0.219086362515982043995534934228 D+00
        weight(9) =  0.149451349150580593145776339658 D+00
        weight(10) = 0.666713443086881375935688098933 D-01

      else if ( norder .eq. 11 ) then

        xtab(1) =  - 0.978228658146056992803938001123 D+00
        xtab(2) =  - 0.887062599768095299075157759304 D+00
        xtab(3) =  - 0.730152005574049324093416252031 D+00
        xtab(4) =  - 0.519096129206811815925725669459 D+00
        xtab(5) =  - 0.269543155952344972331531985401 D+00
        xtab(6) =    0.0 D+00
        xtab(7) =    0.269543155952344972331531985401 D+00
        xtab(8) =    0.519096129206811815925725669459 D+00
        xtab(9) =    0.730152005574049324093416252031 D+00
        xtab(10) =   0.887062599768095299075157759304 D+00
        xtab(11) =   0.978228658146056992803938001123 D+00

        weight(1) =  0.556685671161736664827537204425 D-01
        weight(2) =  0.125580369464904624634694299224 D+00
        weight(3) =  0.186290210927734251426097641432 D+00
        weight(4) =  0.233193764591990479918523704843 D+00
        weight(5) =  0.262804544510246662180688869891 D+00
        weight(6) =  0.272925086777900630714483528336 D+00
        weight(7) =  0.262804544510246662180688869891 D+00
        weight(8) =  0.233193764591990479918523704843 D+00
        weight(9) =  0.186290210927734251426097641432 D+00
        weight(10) = 0.125580369464904624634694299224 D+00
        weight(11) = 0.556685671161736664827537204425 D-01

      else if ( norder .eq. 12 ) then

        xtab(1) =  - 0.981560634246719250690549090149 D+00
        xtab(2) =  - 0.904117256370474856678465866119 D+00
        xtab(3) =  - 0.769902674194304687036893833213 D+00
        xtab(4) =  - 0.587317954286617447296702418941 D+00
        xtab(5) =  - 0.367831498998180193752691536644 D+00
        xtab(6) =  - 0.125233408511468915472441369464 D+00
        xtab(7) =    0.125233408511468915472441369464 D+00
        xtab(8) =    0.367831498998180193752691536644 D+00
        xtab(9) =    0.587317954286617447296702418941 D+00
        xtab(10) =   0.769902674194304687036893833213 D+00
        xtab(11) =   0.904117256370474856678465866119 D+00
        xtab(12) =   0.981560634246719250690549090149 D+00

        weight(1) =  0.471753363865118271946159614850 D-01
        weight(2) =  0.106939325995318430960254718194 D+00
        weight(3) =  0.160078328543346226334652529543 D+00
        weight(4) =  0.203167426723065921749064455810 D+00
        weight(5) =  0.233492536538354808760849898925 D+00
        weight(6) =  0.249147045813402785000562436043 D+00
        weight(7) =  0.249147045813402785000562436043 D+00
        weight(8) =  0.233492536538354808760849898925 D+00
        weight(9) =  0.203167426723065921749064455810 D+00
        weight(10) = 0.160078328543346226334652529543 D+00
        weight(11) = 0.106939325995318430960254718194 D+00
        weight(12) = 0.471753363865118271946159614850 D-01

      else if ( norder .eq. 13 ) then

        xtab(1) =  - 0.984183054718588149472829448807 D+00
        xtab(2) =  - 0.917598399222977965206547836501 D+00
        xtab(3) =  - 0.801578090733309912794206489583 D+00
        xtab(4) =  - 0.642349339440340220643984606996 D+00
        xtab(5) =  - 0.448492751036446852877912852128 D+00
        xtab(6) =  - 0.230458315955134794065528121098 D+00
        xtab(7) =    0.0 D+00
        xtab(8) =    0.230458315955134794065528121098 D+00
        xtab(9) =    0.448492751036446852877912852128 D+00
        xtab(10) =   0.642349339440340220643984606996 D+00
        xtab(11) =   0.801578090733309912794206489583 D+00
        xtab(12) =   0.917598399222977965206547836501 D+00
        xtab(13) =   0.984183054718588149472829448807 D+00

        weight(1) =  0.404840047653158795200215922010 D-01
        weight(2) =  0.921214998377284479144217759538 D-01
        weight(3) =  0.138873510219787238463601776869 D+00
        weight(4) =  0.178145980761945738280046691996 D+00
        weight(5) =  0.207816047536888502312523219306 D+00
        weight(6) =  0.226283180262897238412090186040 D+00
        weight(7) =  0.232551553230873910194589515269 D+00
        weight(8) =  0.226283180262897238412090186040 D+00
        weight(9) =  0.207816047536888502312523219306 D+00
        weight(10) = 0.178145980761945738280046691996 D+00
        weight(11) = 0.138873510219787238463601776869 D+00
        weight(12) = 0.921214998377284479144217759538 D-01
        weight(13) = 0.404840047653158795200215922010 D-01

      else if ( norder .eq. 14 ) then

        xtab(1) =  - 0.986283808696812338841597266704 D+00
        xtab(2) =  - 0.928434883663573517336391139378 D+00
        xtab(3) =  - 0.827201315069764993189794742650 D+00
        xtab(4) =  - 0.687292904811685470148019803019 D+00
        xtab(5) =  - 0.515248636358154091965290718551 D+00
        xtab(6) =  - 0.319112368927889760435671824168 D+00
        xtab(7) =  - 0.108054948707343662066244650220 D+00
        xtab(8) =    0.108054948707343662066244650220 D+00
        xtab(9) =    0.319112368927889760435671824168 D+00
        xtab(10) =   0.515248636358154091965290718551 D+00
        xtab(11) =   0.687292904811685470148019803019 D+00
        xtab(12) =   0.827201315069764993189794742650 D+00
        xtab(13) =   0.928434883663573517336391139378 D+00
        xtab(14) =   0.986283808696812338841597266704 D+00

        weight(1) =  0.351194603317518630318328761382 D-01
        weight(2) =  0.801580871597602098056332770629 D-01
        weight(3) =  0.121518570687903184689414809072 D+00
        weight(4) =  0.157203167158193534569601938624 D+00
        weight(5) =  0.185538397477937813741716590125 D+00
        weight(6) =  0.205198463721295603965924065661 D+00
        weight(7) =  0.215263853463157790195876443316 D+00
        weight(8) =  0.215263853463157790195876443316 D+00
        weight(9) =  0.205198463721295603965924065661 D+00
        weight(10) = 0.185538397477937813741716590125 D+00
        weight(11) = 0.157203167158193534569601938624 D+00
        weight(12) = 0.121518570687903184689414809072 D+00
        weight(13) = 0.801580871597602098056332770629 D-01
        weight(14) = 0.351194603317518630318328761382 D-01

      else if ( norder .eq. 15 ) then

        xtab(1) =  - 0.987992518020485428489565718587 D+00
        xtab(2) =  - 0.937273392400705904307758947710 D+00
        xtab(3) =  - 0.848206583410427216200648320774 D+00
        xtab(4) =  - 0.724417731360170047416186054614 D+00
        xtab(5) =  - 0.570972172608538847537226737254 D+00
        xtab(6) =  - 0.394151347077563369897207370981 D+00
        xtab(7) =  - 0.201194093997434522300628303395 D+00
        xtab(8) =    0.0 D+00
        xtab(9) =    0.201194093997434522300628303395 D+00
        xtab(10) =   0.394151347077563369897207370981 D+00
        xtab(11) =   0.570972172608538847537226737254 D+00
        xtab(12) =   0.724417731360170047416186054614 D+00
        xtab(13) =   0.848206583410427216200648320774 D+00
        xtab(14) =   0.937273392400705904307758947710 D+00
        xtab(15) =   0.987992518020485428489565718587 D+00

        weight(1) =  0.307532419961172683546283935772 D-01
        weight(2) =  0.703660474881081247092674164507 D-01
        weight(3) =  0.107159220467171935011869546686 D+00
        weight(4) =  0.139570677926154314447804794511 D+00
        weight(5) =  0.166269205816993933553200860481 D+00
        weight(6) =  0.186161000015562211026800561866 D+00
        weight(7) =  0.198431485327111576456118326444 D+00
        weight(8) =  0.202578241925561272880620199968 D+00
        weight(9) =  0.198431485327111576456118326444 D+00
        weight(10) = 0.186161000015562211026800561866 D+00
        weight(11) = 0.166269205816993933553200860481 D+00
        weight(12) = 0.139570677926154314447804794511 D+00
        weight(13) = 0.107159220467171935011869546686 D+00
        weight(14) = 0.703660474881081247092674164507 D-01
        weight(15) = 0.307532419961172683546283935772 D-01

      else if ( norder .eq. 16 ) then

        xtab(1) =  - 0.9894009349916499325961541734 D+00
        xtab(2) =  - 0.9445750230732325760779884155 D+00
        xtab(3) =  - 0.8656312023878317438804678977 D+00
        xtab(4) =  - 0.7554044083550030338951011948 D+00
        xtab(5) =  - 0.6178762444026437484466717640 D+00
        xtab(6) =  - 0.4580167776572273863424194429 D+00
        xtab(7) =  - 0.2816035507792589132304605014 D+00
        xtab(8) =  - 0.9501250983763744018531933542 D-01
        xtab(9) =    0.9501250983763744018531933542 D-01
        xtab(10) =   0.2816035507792589132304605014 D+00
        xtab(11) =   0.4580167776572273863424194429 D+00
        xtab(12) =   0.6178762444026437484466717640 D+00
        xtab(13) =   0.7554044083550030338951011948 D+00
        xtab(14) =   0.8656312023878317438804678977 D+00
        xtab(15) =   0.9445750230732325760779884155 D+00
        xtab(16) =   0.9894009349916499325961541734 D+00

        weight(1) =  0.2715245941175409485178057245 D-01
        weight(2) =  0.6225352393864789286284383699 D-01
        weight(3) =  0.9515851168249278480992510760 D-01
        weight(4) =  0.1246289712555338720524762821 D+00
        weight(5) =  0.1495959888165767320815017305 D+00
        weight(6) =  0.1691565193950025381893120790 D+00
        weight(7) =  0.1826034150449235888667636679 D+00
        weight(8) =  0.1894506104550684962853957232 D+00
        weight(9) =  0.1894506104550684962853957232 D+00
        weight(10) = 0.1826034150449235888667636679 D+00
        weight(11) = 0.1691565193950025381893120790 D+00
        weight(12) = 0.1495959888165767320815017305 D+00
        weight(13) = 0.1246289712555338720524762821 D+00
        weight(14) = 0.9515851168249278480992510760 D-01
        weight(15) = 0.6225352393864789286284383699 D-01
        weight(16) = 0.2715245941175409485178057245 D-01

      else if ( norder .eq. 17 ) then

        xtab(1) =  - 0.990575475314417335675434019941 D+00
        xtab(2) =  - 0.950675521768767761222716957896 D+00
        xtab(3) =  - 0.880239153726985902122955694488 D+00
        xtab(4) =  - 0.781514003896801406925230055520 D+00
        xtab(5) =  - 0.657671159216690765850302216643 D+00
        xtab(6) =  - 0.512690537086476967886246568630 D+00
        xtab(7) =  - 0.351231763453876315297185517095 D+00
        xtab(8) =  - 0.178484181495847855850677493654 D+00
        xtab(9) =    0.0 D+00
        xtab(10) =   0.178484181495847855850677493654 D+00
        xtab(11) =   0.351231763453876315297185517095 D+00
        xtab(12) =   0.512690537086476967886246568630 D+00
        xtab(13) =   0.657671159216690765850302216643 D+00
        xtab(14) =   0.781514003896801406925230055520 D+00
        xtab(15) =   0.880239153726985902122955694488 D+00
        xtab(16) =   0.950675521768767761222716957896 D+00
        xtab(17) =   0.990575475314417335675434019941 D+00

        weight(1) =  0.241483028685479319601100262876 D-01
        weight(2) =  0.554595293739872011294401653582 D-01
        weight(3) =  0.850361483171791808835353701911 D-01
        weight(4) =  0.111883847193403971094788385626 D+00
        weight(5) =  0.135136368468525473286319981702 D+00
        weight(6) =  0.154045761076810288081431594802 D+00
        weight(7) =  0.168004102156450044509970663788 D+00
        weight(8) =  0.176562705366992646325270990113 D+00
        weight(9) =  0.179446470356206525458265644262 D+00
        weight(10) = 0.176562705366992646325270990113 D+00
        weight(11) = 0.168004102156450044509970663788 D+00
        weight(12) = 0.154045761076810288081431594802 D+00
        weight(13) = 0.135136368468525473286319981702 D+00
        weight(14) = 0.111883847193403971094788385626 D+00
        weight(15) = 0.850361483171791808835353701911 D-01
        weight(16) = 0.554595293739872011294401653582 D-01
        weight(17) = 0.241483028685479319601100262876 D-01

      else if ( norder .eq. 18 ) then

        xtab(1) =  - 0.991565168420930946730016004706 D+00
        xtab(2) =  - 0.955823949571397755181195892930 D+00
        xtab(3) =  - 0.892602466497555739206060591127 D+00
        xtab(4) =  - 0.803704958972523115682417455015 D+00
        xtab(5) =  - 0.691687043060353207874891081289 D+00
        xtab(6) =  - 0.559770831073947534607871548525 D+00
        xtab(7) =  - 0.411751161462842646035931793833 D+00
        xtab(8) =  - 0.251886225691505509588972854878 D+00
        xtab(9) =  - 0.847750130417353012422618529358 D-01
        xtab(10) =   0.847750130417353012422618529358 D-01
        xtab(11) =   0.251886225691505509588972854878 D+00
        xtab(12) =   0.411751161462842646035931793833 D+00
        xtab(13) =   0.559770831073947534607871548525 D+00
        xtab(14) =   0.691687043060353207874891081289 D+00
        xtab(15) =   0.803704958972523115682417455015 D+00
        xtab(16) =   0.892602466497555739206060591127 D+00
        xtab(17) =   0.955823949571397755181195892930 D+00
        xtab(18) =   0.991565168420930946730016004706 D+00

        weight(1) =  0.216160135264833103133427102665 D-01
        weight(2) =  0.497145488949697964533349462026 D-01
        weight(3) =  0.764257302548890565291296776166 D-01
        weight(4) =  0.100942044106287165562813984925 D+00
        weight(5) =  0.122555206711478460184519126800 D+00
        weight(6) =  0.140642914670650651204731303752 D+00
        weight(7) =  0.154684675126265244925418003836 D+00
        weight(8) =  0.164276483745832722986053776466 D+00
        weight(9) =  0.169142382963143591840656470135 D+00
        weight(10) = 0.169142382963143591840656470135 D+00
        weight(11) = 0.164276483745832722986053776466 D+00
        weight(12) = 0.154684675126265244925418003836 D+00
        weight(13) = 0.140642914670650651204731303752 D+00
        weight(14) = 0.122555206711478460184519126800 D+00
        weight(15) = 0.100942044106287165562813984925 D+00
        weight(16) = 0.764257302548890565291296776166 D-01
        weight(17) = 0.497145488949697964533349462026 D-01
        weight(18) = 0.216160135264833103133427102665 D-01

      else if ( norder .eq. 19 ) then

        xtab(1) =  - 0.992406843843584403189017670253 D+00
        xtab(2) =  - 0.960208152134830030852778840688 D+00
        xtab(3) =  - 0.903155903614817901642660928532 D+00
        xtab(4) =  - 0.822714656537142824978922486713 D+00
        xtab(5) =  - 0.720966177335229378617095860824 D+00
        xtab(6) =  - 0.600545304661681023469638164946 D+00
        xtab(7) =  - 0.464570741375960945717267148104 D+00
        xtab(8) =  - 0.316564099963629831990117328850 D+00
        xtab(9) =  - 0.160358645640225375868096115741 D+00
        xtab(10) =   0.0 D+00
        xtab(11) =   0.160358645640225375868096115741 D+00
        xtab(12) =   0.316564099963629831990117328850 D+00
        xtab(13) =   0.464570741375960945717267148104 D+00
        xtab(14) =   0.600545304661681023469638164946 D+00
        xtab(15) =   0.720966177335229378617095860824 D+00
        xtab(16) =   0.822714656537142824978922486713 D+00
        xtab(17) =   0.903155903614817901642660928532 D+00
        xtab(18) =   0.960208152134830030852778840688 D+00
        xtab(19) =   0.992406843843584403189017670253 D+00

        weight(1) =  0.194617882297264770363120414644 D-01
        weight(2) =  0.448142267656996003328381574020 D-01
        weight(3) =  0.690445427376412265807082580060 D-01
        weight(4) =  0.914900216224499994644620941238 D-01
        weight(5) =  0.111566645547333994716023901682 D+00
        weight(6) =  0.128753962539336227675515784857 D+00
        weight(7) =  0.142606702173606611775746109442 D+00
        weight(8) =  0.152766042065859666778855400898 D+00
        weight(9) =  0.158968843393954347649956439465 D+00
        weight(10) = 0.161054449848783695979163625321 D+00
        weight(11) = 0.158968843393954347649956439465 D+00
        weight(12) = 0.152766042065859666778855400898 D+00
        weight(13) = 0.142606702173606611775746109442 D+00
        weight(14) = 0.128753962539336227675515784857 D+00
        weight(15) = 0.111566645547333994716023901682 D+00
        weight(16) = 0.914900216224499994644620941238 D-01
        weight(17) = 0.690445427376412265807082580060 D-01
        weight(18) = 0.448142267656996003328381574020 D-01
        weight(19) = 0.194617882297264770363120414644 D-01

      else if ( norder .eq. 20 ) then

        xtab(1) =  - 0.993128599185094924786122388471 D+00
        xtab(2) =  - 0.963971927277913791267666131197 D+00
        xtab(3) =  - 0.912234428251325905867752441203 D+00
        xtab(4) =  - 0.839116971822218823394529061702 D+00
        xtab(5) =  - 0.746331906460150792614305070356 D+00
        xtab(6) =  - 0.636053680726515025452836696226 D+00
        xtab(7) =  - 0.510867001950827098004364050955 D+00
        xtab(8) =  - 0.373706088715419560672548177025 D+00
        xtab(9) =  - 0.227785851141645078080496195369 D+00
        xtab(10) = - 0.765265211334973337546404093988 D-01
        xtab(11) =   0.765265211334973337546404093988 D-01
        xtab(12) =   0.227785851141645078080496195369 D+00
        xtab(13) =   0.373706088715419560672548177025 D+00
        xtab(14) =   0.510867001950827098004364050955 D+00
        xtab(15) =   0.636053680726515025452836696226 D+00
        xtab(16) =   0.746331906460150792614305070356 D+00
        xtab(17) =   0.839116971822218823394529061702 D+00
        xtab(18) =   0.912234428251325905867752441203 D+00
        xtab(19) =   0.963971927277913791267666131197 D+00
        xtab(20) =   0.993128599185094924786122388471 D+00

        weight(1) =  0.176140071391521183118619623519 D-01
        weight(2) =  0.406014298003869413310399522749 D-01
        weight(3) =  0.626720483341090635695065351870 D-01
        weight(4) =  0.832767415767047487247581432220 D-01
        weight(5) =  0.101930119817240435036750135480 D+00
        weight(6) =  0.118194531961518417312377377711 D+00
        weight(7) =  0.131688638449176626898494499748 D+00
        weight(8) =  0.142096109318382051329298325067 D+00
        weight(9) =  0.149172986472603746787828737002 D+00
        weight(10) = 0.152753387130725850698084331955 D+00
        weight(11) = 0.152753387130725850698084331955 D+00
        weight(12) = 0.149172986472603746787828737002 D+00
        weight(13) = 0.142096109318382051329298325067 D+00
        weight(14) = 0.131688638449176626898494499748 D+00
        weight(15) = 0.118194531961518417312377377711 D+00
        weight(16) = 0.101930119817240435036750135480 D+00
        weight(17) = 0.832767415767047487247581432220 D-01
        weight(18) = 0.626720483341090635695065351870 D-01
        weight(19) = 0.406014298003869413310399522749 D-01
        weight(20) = 0.176140071391521183118619623519 D-01

      else if ( norder .eq. 32 ) then

        xtab(1) =  - 0.997263861849481563544981128665 D+00
        xtab(2) =  - 0.985611511545268335400175044631 D+00
        xtab(3) =  - 0.964762255587506430773811928118 D+00
        xtab(4) =  - 0.934906075937739689170919134835 D+00
        xtab(5) =  - 0.896321155766052123965307243719 D+00
        xtab(6) =  - 0.849367613732569970133693004968 D+00
        xtab(7) =  - 0.794483795967942406963097298970 D+00
        xtab(8) =  - 0.732182118740289680387426665091 D+00
        xtab(9) =  - 0.663044266930215200975115168663 D+00
        xtab(10) = - 0.587715757240762329040745476402 D+00
        xtab(11) = - 0.506899908932229390023747474378 D+00
        xtab(12) = - 0.421351276130635345364119436172 D+00
        xtab(13) = - 0.331868602282127649779916805730 D+00
        xtab(14) = - 0.239287362252137074544603209166 D+00
        xtab(15) = - 0.144471961582796493485186373599 D+00
        xtab(16) = - 0.483076656877383162348125704405 D-01
        xtab(17) =   0.483076656877383162348125704405 D-01
        xtab(18) =   0.144471961582796493485186373599 D+00
        xtab(19) =   0.239287362252137074544603209166 D+00
        xtab(20) =   0.331868602282127649779916805730 D+00
        xtab(21) =   0.421351276130635345364119436172 D+00
        xtab(22) =   0.506899908932229390023747474378 D+00
        xtab(23) =   0.587715757240762329040745476402 D+00
        xtab(24) =   0.663044266930215200975115168663 D+00
        xtab(25) =   0.732182118740289680387426665091 D+00
        xtab(26) =   0.794483795967942406963097298970 D+00
        xtab(27) =   0.849367613732569970133693004968 D+00
        xtab(28) =   0.896321155766052123965307243719 D+00
        xtab(29) =   0.934906075937739689170919134835 D+00
        xtab(30) =   0.964762255587506430773811928118 D+00
        xtab(31) =   0.985611511545268335400175044631 D+00
        xtab(32) =   0.997263861849481563544981128665 D+00

        weight(1) =  0.701861000947009660040706373885 D-02
        weight(2) =  0.162743947309056706051705622064 D-01
        weight(3) =  0.253920653092620594557525897892 D-01
        weight(4) =  0.342738629130214331026877322524 D-01
        weight(5) =  0.428358980222266806568786466061 D-01
        weight(6) =  0.509980592623761761961632446895 D-01
        weight(7) =  0.586840934785355471452836373002 D-01
        weight(8) =  0.658222227763618468376500637069 D-01
        weight(9) =  0.723457941088485062253993564785 D-01
        weight(10) = 0.781938957870703064717409188283 D-01
        weight(11) = 0.833119242269467552221990746043 D-01
        weight(12) = 0.876520930044038111427714627518 D-01
        weight(13) = 0.911738786957638847128685771116 D-01
        weight(14) = 0.938443990808045656391802376681 D-01
        weight(15) = 0.956387200792748594190820022041 D-01
        weight(16) = 0.965400885147278005667648300636 D-01
        weight(17) = 0.965400885147278005667648300636 D-01
        weight(18) = 0.956387200792748594190820022041 D-01
        weight(19) = 0.938443990808045656391802376681 D-01
        weight(20) = 0.911738786957638847128685771116 D-01
        weight(21) = 0.876520930044038111427714627518 D-01
        weight(22) = 0.833119242269467552221990746043 D-01
        weight(23) = 0.781938957870703064717409188283 D-01
        weight(24) = 0.723457941088485062253993564785 D-01
        weight(25) = 0.658222227763618468376500637069 D-01
        weight(26) = 0.586840934785355471452836373002 D-01
        weight(27) = 0.509980592623761761961632446895 D-01
        weight(28) = 0.428358980222266806568786466061 D-01
        weight(29) = 0.342738629130214331026877322524 D-01
        weight(30) = 0.253920653092620594557525897892 D-01
        weight(31) = 0.162743947309056706051705622064 D-01
        weight(32) = 0.701861000947009660040706373885 D-02

      else if ( norder .eq. 64 ) then

        xtab(1) =  - 0.999305041735772139456905624346 D+00
        xtab(2) =  - 0.996340116771955279346924500676 D+00
        xtab(3) =  - 0.991013371476744320739382383443 D+00
        xtab(4) =  - 0.983336253884625956931299302157 D+00
        xtab(5) =  - 0.973326827789910963741853507352 D+00
        xtab(6) =  - 0.961008799652053718918614121897 D+00
        xtab(7) =  - 0.946411374858402816062481491347 D+00
        xtab(8) =  - 0.929569172131939575821490154559 D+00
        xtab(9) =  - 0.910522137078502805756380668008 D+00
        xtab(10) = - 0.889315445995114105853404038273 D+00
        xtab(11) = - 0.865999398154092819760783385070 D+00
        xtab(12) = - 0.840629296252580362751691544696 D+00
        xtab(13) = - 0.813265315122797559741923338086 D+00
        xtab(14) = - 0.783972358943341407610220525214 D+00
        xtab(15) = - 0.752819907260531896611863774886 D+00
        xtab(16) = - 0.719881850171610826848940217832 D+00
        xtab(17) = - 0.685236313054233242563558371031 D+00
        xtab(18) = - 0.648965471254657339857761231993 D+00
        xtab(19) = - 0.611155355172393250248852971019 D+00
        xtab(20) = - 0.571895646202634034283878116659 D+00
        xtab(21) = - 0.531279464019894545658013903544 D+00
        xtab(22) = - 0.489403145707052957478526307022 D+00
        xtab(23) = - 0.446366017253464087984947714759 D+00
        xtab(24) = - 0.402270157963991603695766771260 D+00
        xtab(25) = - 0.357220158337668115950442615046 D+00
        xtab(26) = - 0.311322871990210956157512698560 D+00
        xtab(27) = - 0.264687162208767416373964172510 D+00
        xtab(28) = - 0.217423643740007084149648748989 D+00
        xtab(29) = - 0.169644420423992818037313629748 D+00
        xtab(30) = - 0.121462819296120554470376463492 D+00
        xtab(31) = - 0.729931217877990394495429419403 D-01
        xtab(32) = - 0.243502926634244325089558428537 D-01
        xtab(33) =   0.243502926634244325089558428537 D-01
        xtab(34) =   0.729931217877990394495429419403 D-01
        xtab(35) =   0.121462819296120554470376463492 D+00
        xtab(36) =   0.169644420423992818037313629748 D+00
        xtab(37) =   0.217423643740007084149648748989 D+00
        xtab(38) =   0.264687162208767416373964172510 D+00
        xtab(39) =   0.311322871990210956157512698560 D+00
        xtab(40) =   0.357220158337668115950442615046 D+00
        xtab(41) =   0.402270157963991603695766771260 D+00
        xtab(42) =   0.446366017253464087984947714759 D+00
        xtab(43) =   0.489403145707052957478526307022 D+00
        xtab(44) =   0.531279464019894545658013903544 D+00
        xtab(45) =   0.571895646202634034283878116659 D+00
        xtab(46) =   0.611155355172393250248852971019 D+00
        xtab(47) =   0.648965471254657339857761231993 D+00
        xtab(48) =   0.685236313054233242563558371031 D+00
        xtab(49) =   0.719881850171610826848940217832 D+00
        xtab(50) =   0.752819907260531896611863774886 D+00
        xtab(51) =   0.783972358943341407610220525214 D+00
        xtab(52) =   0.813265315122797559741923338086 D+00
        xtab(53) =   0.840629296252580362751691544696 D+00
        xtab(54) =   0.865999398154092819760783385070 D+00
        xtab(55) =   0.889315445995114105853404038273 D+00
        xtab(56) =   0.910522137078502805756380668008 D+00
        xtab(57) =   0.929569172131939575821490154559 D+00
        xtab(58) =   0.946411374858402816062481491347 D+00
        xtab(59) =   0.961008799652053718918614121897 D+00
        xtab(60) =   0.973326827789910963741853507352 D+00
        xtab(61) =   0.983336253884625956931299302157 D+00
        xtab(62) =   0.991013371476744320739382383443 D+00
        xtab(63) =   0.996340116771955279346924500676 D+00
        xtab(64) =   0.999305041735772139456905624346 D+00

        weight(1) =  0.178328072169643294729607914497 D-02
        weight(2) =  0.414703326056246763528753572855 D-02
        weight(3) =  0.650445796897836285611736039998 D-02
        weight(4) =  0.884675982636394772303091465973 D-02
        weight(5) =  0.111681394601311288185904930192 D-01
        weight(6) =  0.134630478967186425980607666860 D-01
        weight(7) =  0.157260304760247193219659952975 D-01
        weight(8) =  0.179517157756973430850453020011 D-01
        weight(9) =  0.201348231535302093723403167285 D-01
        weight(10) = 0.222701738083832541592983303842 D-01
        weight(11) = 0.243527025687108733381775504091 D-01
        weight(12) = 0.263774697150546586716917926252 D-01
        weight(13) = 0.283396726142594832275113052002 D-01
        weight(14) = 0.302346570724024788679740598195 D-01
        weight(15) = 0.320579283548515535854675043479 D-01
        weight(16) = 0.338051618371416093915654821107 D-01
        weight(17) = 0.354722132568823838106931467152 D-01
        weight(18) = 0.370551285402400460404151018096 D-01
        weight(19) = 0.385501531786156291289624969468 D-01
        weight(20) = 0.399537411327203413866569261283 D-01
        weight(21) = 0.412625632426235286101562974736 D-01
        weight(22) = 0.424735151236535890073397679088 D-01
        weight(23) = 0.435837245293234533768278609737 D-01
        weight(24) = 0.445905581637565630601347100309 D-01
        weight(25) = 0.454916279274181444797709969713 D-01
        weight(26) = 0.462847965813144172959532492323 D-01
        weight(27) = 0.469681828162100173253262857546 D-01
        weight(28) = 0.475401657148303086622822069442 D-01
        weight(29) = 0.479993885964583077281261798713 D-01
        weight(30) = 0.483447622348029571697695271580 D-01
        weight(31) = 0.485754674415034269347990667840 D-01
        weight(32) = 0.486909570091397203833653907347 D-01
        weight(33) = 0.486909570091397203833653907347 D-01
        weight(34) = 0.485754674415034269347990667840 D-01
        weight(35) = 0.483447622348029571697695271580 D-01
        weight(36) = 0.479993885964583077281261798713 D-01
        weight(37) = 0.475401657148303086622822069442 D-01
        weight(38) = 0.469681828162100173253262857546 D-01
        weight(39) = 0.462847965813144172959532492323 D-01
        weight(40) = 0.454916279274181444797709969713 D-01
        weight(41) = 0.445905581637565630601347100309 D-01
        weight(42) = 0.435837245293234533768278609737 D-01
        weight(43) = 0.424735151236535890073397679088 D-01
        weight(44) = 0.412625632426235286101562974736 D-01
        weight(45) = 0.399537411327203413866569261283 D-01
        weight(46) = 0.385501531786156291289624969468 D-01
        weight(47) = 0.370551285402400460404151018096 D-01
        weight(48) = 0.354722132568823838106931467152 D-01
        weight(49) = 0.338051618371416093915654821107 D-01
        weight(50) = 0.320579283548515535854675043479 D-01
        weight(51) = 0.302346570724024788679740598195 D-01
        weight(52) = 0.283396726142594832275113052002 D-01
        weight(53) = 0.263774697150546586716917926252 D-01
        weight(54) = 0.243527025687108733381775504091 D-01
        weight(55) = 0.222701738083832541592983303842 D-01
        weight(56) = 0.201348231535302093723403167285 D-01
        weight(57) = 0.179517157756973430850453020011 D-01
        weight(58) = 0.157260304760247193219659952975 D-01
        weight(59) = 0.134630478967186425980607666860 D-01
        weight(60) = 0.111681394601311288185904930192 D-01
        weight(61) = 0.884675982636394772303091465973 D-02
        weight(62) = 0.650445796897836285611736039998 D-02
        weight(63) = 0.414703326056246763528753572855 D-02
        weight(64) = 0.178328072169643294729607914497 D-02

      else

        write ( *, * ) ' '
        write ( *, * ) 'LEGSET - Fatal error!'
        write ( *, * ) '  Illegal value of NORDER = ', norder
        write ( *, * ) '  Legal values are 1 to 20, 32 or 64.'
        stop

      end if

      return
      end subroutine legset
      subroutine lobset ( norder, xtab, weight )
!
!***********************************************************************
!
!! LOBSET sets abscissas and weights for Lobatto quadrature.
!
!
!  Integration interval:
!
!    [ -1, 1 ].
!
!  Weight function:
!
!    1.
!
!  Integral to approximate:
!
!    INTEGRAL ( -1 <= X <= 1 ) F(X) dX.
!
!  Approximate integral:
!
!    SUM ( I = 1 to NORDER ) WEIGHT(I) * F ( XTAB(I) ).
!
!  Precision:
!
!    The quadrature rule will integrate exactly all polynomials up to
!    X**(2*NORDER-3).
!
!  Note:
!
!    The Lobatto rule is distinguished by the fact that both endpoints
!    (-1 and 1) are always abscissas of the rule.
!
!  Reference:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964.
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Modified:
!
!    20 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!    NORDER must be between 2 and 20.
!
!    Output, real(kind=Rkind) XTAB(NORDER), the abscissas for the rule.
!
!    Output, real(kind=Rkind) WEIGHT(NORDER), the weights of the rule.
!    The weights are positive, symmetric and should sum to 2.
!
      integer norder
!
      real(kind=Rkind) xtab(norder)
      real(kind=Rkind) weight(norder)
!
      if ( norder .eq. 2 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =    1.0 D+00

        weight(1) =  1.0 D+00
        weight(2) =  1.0 D+00

      else if ( norder .eq. 3 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =    0.0 D+00
        xtab(3) =    1.0 D+00

        weight(1) =  1.0 D+00 / 3.0 D+00
        weight(2) =  4.0 D+00 / 3.0 D+00
        weight(3) =  1.0 D+00 / 3.0 D+00

      else if ( norder .eq. 4 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =  - 0.447213595499957939281834733746 D+00
        xtab(3) =    0.447213595499957939281834733746 D+00
        xtab(4) =    1.0 D+00

        weight(1) =  1.0 D+00 / 6.0 D+00
        weight(2) =  5.0 D+00 / 6.0 D+00
        weight(3) =  5.0 D+00 / 6.0 D+00
        weight(4) =  1.0 D+00 / 6.0 D+00

      else if ( norder .eq. 5 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =  - 0.654653670707977143798292456247 D+00
        xtab(3) =    0.0 D+00
        xtab(4) =    0.654653670707977143798292456247 D+00
        xtab(5) =    1.0 D+00

        weight(1) =  9.0 D+00 / 90.0 D+00
        weight(2) = 49.0 D+00 / 90.0 D+00
        weight(3) = 64.0 D+00 / 90.0 D+00
        weight(4) = 49.0 D+00 / 90.0 D+00
        weight(5) =  9.0 D+00 / 90.0 D+00

      else if ( norder .eq. 6 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =  - 0.765055323929464692851002973959 D+00
        xtab(3) =  - 0.285231516480645096314150994041 D+00
        xtab(4) =    0.285231516480645096314150994041 D+00
        xtab(5) =    0.765055323929464692851002973959 D+00
        xtab(6) =    1.0 D+00

        weight(1) =  0.066666666666666666666666666667 D+00
        weight(2) =  0.378474956297846980316612808212 D+00
        weight(3) =  0.554858377035486353016720525121 D+00
        weight(4) =  0.554858377035486353016720525121 D+00
        weight(5) =  0.378474956297846980316612808212 D+00
        weight(6) =  0.066666666666666666666666666667 D+00

      else if ( norder .eq. 7 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =  - 0.830223896278566929872032213967 D+00
        xtab(3) =  - 0.468848793470714213803771881909 D+00
        xtab(4) =    0.0 D+00
        xtab(5) =    0.468848793470714213803771881909 D+00
        xtab(6) =    0.830223896278566929872032213967 D+00
        xtab(7) =    1.0 D+00

        weight(1) =  0.476190476190476190476190476190 D-01
        weight(2) =  0.276826047361565948010700406290 D+00
        weight(3) =  0.431745381209862623417871022281 D+00
        weight(4) =  0.487619047619047619047619047619 D+00
        weight(5) =  0.431745381209862623417871022281 D+00
        weight(6) =  0.276826047361565948010700406290 D+00
        weight(7) =  0.476190476190476190476190476190 D-01

      else if ( norder .eq. 8 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =  - 0.871740148509606615337445761221 D+00
        xtab(3) =  - 0.591700181433142302144510731398 D+00
        xtab(4) =  - 0.209299217902478868768657260345 D+00
        xtab(5) =    0.209299217902478868768657260345 D+00
        xtab(6) =    0.591700181433142302144510731398 D+00
        xtab(7) =    0.871740148509606615337445761221 D+00
        xtab(8) =    1.0 D+00

        weight(1) =  0.357142857142857142857142857143 D-01
        weight(2) =  0.210704227143506039382991065776 D+00
        weight(3) =  0.341122692483504364764240677108 D+00
        weight(4) =  0.412458794658703881567052971402 D+00
        weight(5) =  0.412458794658703881567052971402 D+00
        weight(6) =  0.341122692483504364764240677108 D+00
        weight(7) =  0.210704227143506039382991065776 D+00
        weight(8) =  0.357142857142857142857142857143 D-01

      else if ( norder .eq. 9 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =  - 0.899757995411460157312345244418 D+00
        xtab(3) =  - 0.677186279510737753445885427091 D+00
        xtab(4) =  - 0.363117463826178158710752068709 D+00
        xtab(5) =    0.0 D+00
        xtab(6) =    0.363117463826178158710752068709 D+00
        xtab(7) =    0.677186279510737753445885427091 D+00
        xtab(8) =    0.899757995411460157312345244418 D+00
        xtab(9) =    1.0 D+00

        weight(1) =  0.277777777777777777777777777778 D-01
        weight(2) =  0.165495361560805525046339720029 D+00
        weight(3) =  0.274538712500161735280705618579 D+00
        weight(4) =  0.346428510973046345115131532140 D+00
        weight(5) =  0.371519274376417233560090702948 D+00
        weight(6) =  0.346428510973046345115131532140 D+00
        weight(7) =  0.274538712500161735280705618579 D+00
        weight(8) =  0.165495361560805525046339720029 D+00
        weight(9) =  0.277777777777777777777777777778 D-01

      else if ( norder .eq. 10 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =  - 0.919533908166458813828932660822 D+00
        xtab(3) =  - 0.738773865105505075003106174860 D+00
        xtab(4) =  - 0.477924949810444495661175092731 D+00
        xtab(5) =  - 0.165278957666387024626219765958 D+00
        xtab(6) =    0.165278957666387024626219765958 D+00
        xtab(7) =    0.477924949810444495661175092731 D+00
        xtab(8) =    0.738773865105505075003106174860 D+00
        xtab(9) =    0.919533908166458813828932660822 D+00
        xtab(10) =   1.0 D+00

        weight(1) =  0.222222222222222222222222222222 D-01
        weight(2) =  0.133305990851070111126227170755 D+00
        weight(3) =  0.224889342063126452119457821731 D+00
        weight(4) =  0.292042683679683757875582257374 D+00
        weight(5) =  0.327539761183897456656510527917 D+00
        weight(6) =  0.327539761183897456656510527917 D+00
        weight(7) =  0.292042683679683757875582257374 D+00
        weight(8) =  0.224889342063126452119457821731 D+00
        weight(9) =  0.133305990851070111126227170755 D+00
        weight(10) = 0.222222222222222222222222222222 D-01

      else if ( norder .eq. 11 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =  - 0.934001430408059134332274136099 D+00
        xtab(3) =  - 0.784483473663144418622417816108 D+00
        xtab(4) =  - 0.565235326996205006470963969478 D+00
        xtab(5) =  - 0.295758135586939391431911515559 D+00
        xtab(6) =    0.0 D+00
        xtab(7) =    0.295758135586939391431911515559 D+00
        xtab(8) =    0.565235326996205006470963969478 D+00
        xtab(9) =    0.784483473663144418622417816108 D+00
        xtab(10) =   0.934001430408059134332274136099 D+00
        xtab(11) =   1.0 D+00

        weight(1) =  0.181818181818181818181818181818 D-01
        weight(2) =  0.109612273266994864461403449580 D+00
        weight(3) =  0.187169881780305204108141521899 D+00
        weight(4) =  0.248048104264028314040084866422 D+00
        weight(5) =  0.286879124779008088679222403332 D+00
        weight(6) =  0.300217595455690693785931881170 D+00
        weight(7) =  0.286879124779008088679222403332 D+00
        weight(8) =  0.248048104264028314040084866422 D+00
        weight(9) =  0.187169881780305204108141521899 D+00
        weight(10) = 0.109612273266994864461403449580 D+00
        weight(11) = 0.181818181818181818181818181818 D-01

      else if ( norder .eq. 12 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =  - 0.944899272222882223407580138303 D+00
        xtab(3) =  - 0.819279321644006678348641581717 D+00
        xtab(4) =  - 0.632876153031869677662404854444 D+00
        xtab(5) =  - 0.399530940965348932264349791567 D+00
        xtab(6) =  - 0.136552932854927554864061855740 D+00
        xtab(7) =    0.136552932854927554864061855740 D+00
        xtab(8) =    0.399530940965348932264349791567 D+00
        xtab(9) =    0.632876153031869677662404854444 D+00
        xtab(10) =   0.819279321644006678348641581717 D+00
        xtab(11) =   0.944899272222882223407580138303 D+00
        xtab(12) =   1.0 D+00

        weight(1) =  0.151515151515151515151515151515 D-01
        weight(2) =  0.916845174131961306683425941341 D-01
        weight(3) =  0.157974705564370115164671062700 D+00
        weight(4) =  0.212508417761021145358302077367 D+00
        weight(5) =  0.251275603199201280293244412148 D+00
        weight(6) =  0.271405240910696177000288338500 D+00
        weight(7) =  0.271405240910696177000288338500 D+00
        weight(8) =  0.251275603199201280293244412148 D+00
        weight(9) =  0.212508417761021145358302077367 D+00
        weight(10) = 0.157974705564370115164671062700 D+00
        weight(11) = 0.916845174131961306683425941341 D-01
        weight(12) = 0.151515151515151515151515151515 D-01

      else if ( norder .eq. 13 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =  - 0.953309846642163911896905464755 D+00
        xtab(3) =  - 0.846347564651872316865925607099 D+00
        xtab(4) =  - 0.686188469081757426072759039566 D+00
        xtab(5) =  - 0.482909821091336201746937233637 D+00
        xtab(6) =  - 0.249286930106239992568673700374 D+00
        xtab(7) =    0.0 D+00
        xtab(8) =    0.249286930106239992568673700374 D+00
        xtab(9) =    0.482909821091336201746937233637 D+00
        xtab(10) =   0.686188469081757426072759039566 D+00
        xtab(11) =   0.846347564651872316865925607099 D+00
        xtab(12) =   0.953309846642163911896905464755 D+00
        xtab(13) =   1.0 D+00

        weight(1) =  0.128205128205128205128205128205 D-01
        weight(2) =  0.778016867468189277935889883331 D-01
        weight(3) =  0.134981926689608349119914762589 D+00
        weight(4) =  0.183646865203550092007494258747 D+00
        weight(5) =  0.220767793566110086085534008379 D+00
        weight(6) =  0.244015790306676356458578148360 D+00
        weight(7) =  0.251930849333446736044138641541 D+00
        weight(8) =  0.244015790306676356458578148360 D+00
        weight(9) =  0.220767793566110086085534008379 D+00
        weight(10) = 0.183646865203550092007494258747 D+00
        weight(11) = 0.134981926689608349119914762589 D+00
        weight(12) = 0.778016867468189277935889883331 D-01
        weight(13) = 0.128205128205128205128205128205 D-01

      else if ( norder .eq. 14 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =  - 0.959935045267260901355100162015 D+00
        xtab(3) =  - 0.867801053830347251000220202908 D+00
        xtab(4) =  - 0.728868599091326140584672400521 D+00
        xtab(5) =  - 0.550639402928647055316622705859 D+00
        xtab(6) =  - 0.342724013342712845043903403642 D+00
        xtab(7) =  - 0.116331868883703867658776709736 D+00
        xtab(8) =    0.116331868883703867658776709736 D+00
        xtab(9) =    0.342724013342712845043903403642 D+00
        xtab(10) =   0.550639402928647055316622705859 D+00
        xtab(11) =   0.728868599091326140584672400521 D+00
        xtab(12) =   0.867801053830347251000220202908 D+00
        xtab(13) =   0.959935045267260901355100162015 D+00
        xtab(14) =   1.0 D+00

        weight(1) =  0.109890109890109890109890109890 D-01
        weight(2) =  0.668372844976812846340706607461 D-01
        weight(3) =  0.116586655898711651540996670655 D+00
        weight(4) =  0.160021851762952142412820997988 D+00
        weight(5) =  0.194826149373416118640331778376 D+00
        weight(6) =  0.219126253009770754871162523954 D+00
        weight(7) =  0.231612794468457058889628357293 D+00
        weight(8) =  0.231612794468457058889628357293 D+00
        weight(9) =  0.219126253009770754871162523954 D+00
        weight(10) = 0.194826149373416118640331778376 D+00
        weight(11) = 0.160021851762952142412820997988 D+00
        weight(12) = 0.116586655898711651540996670655 D+00
        weight(13) = 0.668372844976812846340706607461 D-01
        weight(14) = 0.109890109890109890109890109890 D-01

      else if ( norder .eq. 15 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =  - 0.965245926503838572795851392070 D+00
        xtab(3) =  - 0.885082044222976298825401631482 D+00
        xtab(4) =  - 0.763519689951815200704118475976 D+00
        xtab(5) =  - 0.606253205469845711123529938637 D+00
        xtab(6) =  - 0.420638054713672480921896938739 D+00
        xtab(7) =  - 0.215353955363794238225679446273 D+00
        xtab(8) =    0.0 D+00
        xtab(9) =    0.215353955363794238225679446273 D+00
        xtab(10) =   0.420638054713672480921896938739 D+00
        xtab(11) =   0.606253205469845711123529938637 D+00
        xtab(12) =   0.763519689951815200704118475976 D+00
        xtab(13) =   0.885082044222976298825401631482 D+00
        xtab(14) =   0.965245926503838572795851392070 D+00
        xtab(15) =   1.0 D+00

        weight(1) =  0.952380952380952380952380952381 D-02
        weight(2) =  0.580298930286012490968805840253 D-01
        weight(3) =  0.101660070325718067603666170789 D+00
        weight(4) =  0.140511699802428109460446805644 D+00
        weight(5) =  0.172789647253600949052077099408 D+00
        weight(6) =  0.196987235964613356092500346507 D+00
        weight(7) =  0.211973585926820920127430076977 D+00
        weight(8) =  0.217048116348815649514950214251 D+00
        weight(9) =  0.211973585926820920127430076977 D+00
        weight(10) = 0.196987235964613356092500346507 D+00
        weight(11) = 0.172789647253600949052077099408 D+00
        weight(12) = 0.140511699802428109460446805644 D+00
        weight(13) = 0.101660070325718067603666170789 D+00
        weight(14) = 0.580298930286012490968805840253 D-01
        weight(15) = 0.952380952380952380952380952381 D-02

      else if ( norder .eq. 16 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =  - 0.969568046270217932952242738367 D+00
        xtab(3) =  - 0.899200533093472092994628261520 D+00
        xtab(4) =  - 0.792008291861815063931088270963 D+00
        xtab(5) =  - 0.652388702882493089467883219641 D+00
        xtab(6) =  - 0.486059421887137611781890785847 D+00
        xtab(7) =  - 0.299830468900763208098353454722 D+00
        xtab(8) =  - 0.101326273521949447843033005046 D+00
        xtab(9) =    0.101326273521949447843033005046 D+00
        xtab(10) =   0.299830468900763208098353454722 D+00
        xtab(11) =   0.486059421887137611781890785847 D+00
        xtab(12) =   0.652388702882493089467883219641 D+00
        xtab(13) =   0.792008291861815063931088270963 D+00
        xtab(14) =   0.899200533093472092994628261520 D+00
        xtab(15) =   0.969568046270217932952242738367 D+00
        xtab(16) =   1.0 D+00

        weight(1) =  0.833333333333333333333333333333 D-02
        weight(2) =  0.508503610059199054032449195655 D-01
        weight(3) =  0.893936973259308009910520801661 D-01
        weight(4) =  0.124255382132514098349536332657 D+00
        weight(5) =  0.154026980807164280815644940485 D+00
        weight(6) =  0.177491913391704125301075669528 D+00
        weight(7) =  0.193690023825203584316913598854 D+00
        weight(8) =  0.201958308178229871489199125411 D+00
        weight(9) =  0.201958308178229871489199125411 D+00
        weight(10) = 0.193690023825203584316913598854 D+00
        weight(11) = 0.177491913391704125301075669528 D+00
        weight(12) = 0.154026980807164280815644940485 D+00
        weight(13) = 0.124255382132514098349536332657 D+00
        weight(14) = 0.893936973259308009910520801661 D-01
        weight(15) = 0.508503610059199054032449195655 D-01
        weight(16) = 0.833333333333333333333333333333 D-02

      else if ( norder .eq. 17 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =  - 0.973132176631418314156979501874 D+00
        xtab(3) =  - 0.910879995915573595623802506398 D+00
        xtab(4) =  - 0.815696251221770307106750553238 D+00
        xtab(5) =  - 0.691028980627684705394919357372 D+00
        xtab(6) =  - 0.541385399330101539123733407504 D+00
        xtab(7) =  - 0.372174433565477041907234680735 D+00
        xtab(8) =  - 0.189511973518317388304263014753 D+00
        xtab(9) =    0.0 D+00
        xtab(10) =   0.189511973518317388304263014753 D+00
        xtab(11) =   0.372174433565477041907234680735 D+00
        xtab(12) =   0.541385399330101539123733407504 D+00
        xtab(13) =   0.691028980627684705394919357372 D+00
        xtab(14) =   0.815696251221770307106750553238 D+00
        xtab(15) =   0.910879995915573595623802506398 D+00
        xtab(16) =   0.973132176631418314156979501874 D+00
        xtab(17) =   1.0 D+00

        weight(1) =  0.735294117647058823529411764706 D-02
        weight(2) =  0.449219405432542096474009546232 D-01
        weight(3) =  0.791982705036871191902644299528 D-01
        weight(4) =  0.110592909007028161375772705220 D+00
        weight(5) =  0.137987746201926559056201574954 D+00
        weight(6) =  0.160394661997621539516328365865 D+00
        weight(7) =  0.177004253515657870436945745363 D+00
        weight(8) =  0.187216339677619235892088482861 D+00
        weight(9) =  0.190661874753469433299407247028 D+00
        weight(10) = 0.187216339677619235892088482861 D+00
        weight(11) = 0.177004253515657870436945745363 D+00
        weight(12) = 0.160394661997621539516328365865 D+00
        weight(13) = 0.137987746201926559056201574954 D+00
        weight(14) = 0.110592909007028161375772705220 D+00
        weight(15) = 0.791982705036871191902644299528 D-01
        weight(16) = 0.449219405432542096474009546232 D-01
        weight(17) = 0.735294117647058823529411764706 D-02

      else if ( norder .eq. 18 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =  - 0.976105557412198542864518924342 D+00
        xtab(3) =  - 0.920649185347533873837854625431 D+00
        xtab(4) =  - 0.835593535218090213713646362328 D+00
        xtab(5) =  - 0.723679329283242681306210365302 D+00
        xtab(6) =  - 0.588504834318661761173535893194 D+00
        xtab(7) =  - 0.434415036912123975342287136741 D+00
        xtab(8) =  - 0.266362652878280984167665332026 D+00
        xtab(9) =  - 0.897490934846521110226450100886 D-01
        xtab(10) =   0.897490934846521110226450100886 D-01
        xtab(11) =   0.266362652878280984167665332026 D+00
        xtab(12) =   0.434415036912123975342287136741 D+00
        xtab(13) =   0.588504834318661761173535893194 D+00
        xtab(14) =   0.723679329283242681306210365302 D+00
        xtab(15) =   0.835593535218090213713646362328 D+00
        xtab(16) =   0.920649185347533873837854625431 D+00
        xtab(17) =   0.976105557412198542864518924342 D+00
        xtab(18) =   1.0 D+00

        weight(1) =  0.653594771241830065359477124183 D-02
        weight(2) =  0.399706288109140661375991764101 D-01
        weight(3) =  0.706371668856336649992229601678 D-01
        weight(4) =  0.990162717175028023944236053187 D-01
        weight(5) =  0.124210533132967100263396358897 D+00
        weight(6) =  0.145411961573802267983003210494 D+00
        weight(7) =  0.161939517237602489264326706700 D+00
        weight(8) =  0.173262109489456226010614403827 D+00
        weight(9) =  0.179015863439703082293818806944 D+00
        weight(10) = 0.179015863439703082293818806944 D+00
        weight(11) = 0.173262109489456226010614403827 D+00
        weight(12) = 0.161939517237602489264326706700 D+00
        weight(13) = 0.145411961573802267983003210494 D+00
        weight(14) = 0.124210533132967100263396358897 D+00
        weight(15) = 0.990162717175028023944236053187 D-01
        weight(16) = 0.706371668856336649992229601678 D-01
        weight(17) = 0.399706288109140661375991764101 D-01
        weight(18) = 0.653594771241830065359477124183 D-02

      else if ( norder .eq. 19 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =  - 0.978611766222080095152634063110 D+00
        xtab(3) =  - 0.928901528152586243717940258797 D+00
        xtab(4) =  - 0.852460577796646093085955970041 D+00
        xtab(5) =  - 0.751494202552613014163637489634 D+00
        xtab(6) =  - 0.628908137265220497766832306229 D+00
        xtab(7) =  - 0.488229285680713502777909637625 D+00
        xtab(8) =  - 0.333504847824498610298500103845 D+00
        xtab(9) =  - 0.169186023409281571375154153445 D+00
        xtab(10) =   0.0 D+00
        xtab(11) =   0.169186023409281571375154153445 D+00
        xtab(12) =   0.333504847824498610298500103845 D+00
        xtab(13) =   0.488229285680713502777909637625 D+00
        xtab(14) =   0.628908137265220497766832306229 D+00
        xtab(15) =   0.751494202552613014163637489634 D+00
        xtab(16) =   0.852460577796646093085955970041 D+00
        xtab(17) =   0.928901528152586243717940258797 D+00
        xtab(18) =   0.978611766222080095152634063110 D+00
        xtab(19) =   1.0 D+00

        weight(1) =  0.584795321637426900584795321637 D-02
        weight(2) =  0.357933651861764771154255690351 D-01
        weight(3) =  0.633818917626297368516956904183 D-01
        weight(4) =  0.891317570992070844480087905562 D-01
        weight(5) =  0.112315341477305044070910015464 D+00
        weight(6) =  0.132267280448750776926046733910 D+00
        weight(7) =  0.148413942595938885009680643668 D+00
        weight(8) =  0.160290924044061241979910968184 D+00
        weight(9) =  0.167556584527142867270137277740 D+00
        weight(10) = 0.170001919284827234644672715617 D+00
        weight(11) = 0.167556584527142867270137277740 D+00
        weight(12) = 0.160290924044061241979910968184 D+00
        weight(13) = 0.148413942595938885009680643668 D+00
        weight(14) = 0.132267280448750776926046733910 D+00
        weight(15) = 0.112315341477305044070910015464 D+00
        weight(16) = 0.891317570992070844480087905562 D-01
        weight(17) = 0.633818917626297368516956904183 D-01
        weight(18) = 0.357933651861764771154255690351 D-01
        weight(19) = 0.584795321637426900584795321637 D-02

      else if ( norder .eq. 20 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =  - 0.980743704893914171925446438584 D+00
        xtab(3) =  - 0.935934498812665435716181584931 D+00
        xtab(4) =  - 0.866877978089950141309847214616 D+00
        xtab(5) =  - 0.775368260952055870414317527595 D+00
        xtab(6) =  - 0.663776402290311289846403322971 D+00
        xtab(7) =  - 0.534992864031886261648135961829 D+00
        xtab(8) =  - 0.392353183713909299386474703816 D+00
        xtab(9) =  - 0.239551705922986495182401356927 D+00
        xtab(10) = - 0.805459372388218379759445181596 D-01
        xtab(11) =   0.805459372388218379759445181596 D-01
        xtab(12) =   0.239551705922986495182401356927 D+00
        xtab(13) =   0.392353183713909299386474703816 D+00
        xtab(14) =   0.534992864031886261648135961829 D+00
        xtab(15) =   0.663776402290311289846403322971 D+00
        xtab(16) =   0.775368260952055870414317527595 D+00
        xtab(17) =   0.866877978089950141309847214616 D+00
        xtab(18) =   0.935934498812665435716181584931 D+00
        xtab(19) =   0.980743704893914171925446438584 D+00
        xtab(20) =   1.0 D+00

        weight(1) =  0.526315789473684210526315789474 D-02
        weight(2) =  0.322371231884889414916050281173 D-01
        weight(3) =  0.571818021275668260047536271732 D-01
        weight(4) =  0.806317639961196031447768461137 D-01
        weight(5) =  0.101991499699450815683781205733 D+00
        weight(6) =  0.120709227628674725099429705002 D+00
        weight(7) =  0.136300482358724184489780792989 D+00
        weight(8) =  0.148361554070916825814713013734 D+00
        weight(9) =  0.156580102647475487158169896794 D+00
        weight(10) = 0.160743286387845749007726726449 D+00
        weight(11) = 0.160743286387845749007726726449 D+00
        weight(12) = 0.156580102647475487158169896794 D+00
        weight(13) = 0.148361554070916825814713013734 D+00
        weight(14) = 0.136300482358724184489780792989 D+00
        weight(15) = 0.120709227628674725099429705002 D+00
        weight(16) = 0.101991499699450815683781205733 D+00
        weight(17) = 0.806317639961196031447768461137 D-01
        weight(18) = 0.571818021275668260047536271732 D-01
        weight(19) = 0.322371231884889414916050281173 D-01
        weight(20) = 0.526315789473684210526315789474 D-02

      else

        write ( *, * ) ' '
        write ( *, * ) 'LOBSET - Fatal error!'
        write ( *, * ) '  Illegal value of NORDER = ', norder
        write ( *, * ) '  Legal values are between 2 and 20.'
        stop

      end if

      return
      end subroutine lobset
      function log_gamma ( x )
!
!***********************************************************************
!
!! LOG_GAMMA calculates the natural logarithm of GAMMA(X).
!
!
!  The method uses Stirling's approximation, and is accurate to about
!  12 decimal places.
!
!  Reference:
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!  Modified:
!
!    19 September 1998
!
!  Parameters:
!
!    Input, real(kind=Rkind) X, the evaluation point.  The routine
!    will fail if GAMMA(X) is not positive.  X should be greater than 0.
!
!    Output, real(kind=Rkind) LOG_GAMMA, the natural logarithm of the
!    gamma function of X.
!
      real(kind=Rkind) pi
      parameter ( pi = 3.141592653589793238462643383279D+00 )
!
      integer i
      integer k
      real(kind=Rkind) log_gamma
      integer m
      real(kind=Rkind) p
      real(kind=Rkind) x
      real(kind=Rkind) x2
      real(kind=Rkind) y
      real(kind=Rkind) z
!
      if ( x .lt. 0.5 ) then

        m = 1
        x2 = 1.0 - x

      else

        m = 0
        x2 = x

      end if

      k = - 1

!0    continue

      k = k + 1

      if ( x2 + dble ( k ) .le. 6.0 ) then
        go to 10
      end if

      z = x2 + dble ( k )

      y = ( z - 0.5 ) * log ( z ) - z + 0.9189385332047 +               &
        ( ( ( ( (                                                       &
        - 4146.0 / z**2                                                 &
        + 1820.0 ) / z**2                                               &
        - 1287.0 ) / z**2                                               &
        + 1716.0 ) / z**2                                               &
        - 6006.0 ) / z**2                                               &
        + 180180.0 ) / z / 2162160.0

      if ( k .gt. 0 ) then

        do i = 1, k
          y = y - log ( x2 + dble ( k - i ) )
        end do

      end if

      if ( m .ne. 0 ) then

        p = pi / sin ( pi * ( 1.0 - x2 ) )

        if ( p .le. 0.0 ) then

          write ( *, * ) ' '
          write ( *, * ) 'LOG_GAMMA - fatal error!'
          stop

        else

          y = log ( p ) - y

        end if

      end if

      log_gamma = y

      return
      end function log_gamma
      subroutine moulset ( norder, xtab, weight )
!
!***********************************************************************
!
!! MOULSET sets weights for Adams-Moulton quadrature.
!
!
!  Definition:
!
!    Adams-Moulton quadrature formulas are normally used in solving
!    ordinary differential equations, and are not suitable for general
!    quadrature computations.  However, an Adams-Moulton formula is
!    equivalent to approximating the integral of F(Y(X)) between X(M)
!    and X(M+1), using an implicit formula that relies on known values
!    of F(Y(X)) at X(M-N+1) through X(M), plus the unknown value at X(M+1).
!
!    Suppose the unknown function is denoted by Y(X), with derivative F(Y(X)),
!    and that approximate values of the function are known at a series of
!    X values, which we write as X(1), X(2), ..., X(M).  We write the value
!    Y(X(1)) as Y(1) and so on.
!
!    Then the solution of the ODE Y' = F(X,Y) at the next point X(M+1) is
!    computed by:
!
!      Y(M+1) = Y(M) + Integral ( X(M) < X < X(M+1) ) F(Y(X)) dX
!             = Y(M) + H * Sum ( I = 1 to N ) W(I) * F(Y(M+2-I)) approximately.
!
!    Note that this formula is implicit, since the unknown value Y(M+1)
!    appears on the right hand side.  Hence, in ODE applications, this
!    equation must be solved via a nonlinear equation solver.  For
!    quadrature problems, where the function to be integrated is known
!    beforehand, this is not a problem, and the calculation is explicit.
!
!    In the documentation that follows, we replace F(Y(X)) by F(X).
!
!  Integration interval:
!
!    [ 0, 1 ].
!
!  Weight function:
!
!    1.
!
!  Integral to approximate:
!
!    INTEGRAL ( 0 <= X <= 1 ) F(X) dX.
!
!  Approximate integral:
!
!    SUM ( I = 1 to NORDER ) WEIGHT(I) * F ( 2 - I ),
!
!  Note:
!
!    The Adams-Moulton formulas require equally spaced data.
!
!    Here is how the formula is applied in the general case with non-unit spacing:
!
!      INTEGRAL ( A <= X <= A+H ) F(X) dX =
!      H * SUM ( I = 1 to NORDER ) WEIGHT(I) * F ( A - (I-2)*H ), approximately.
!
!  Reference:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964,
!    page 915 ("Lagrangian Integration Coefficients").
!
!    Jean Lapidus and John Seinfeld,
!    Numerical Solution of Ordinary Differential Equations,
!    Academic Press, 1971.
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Modified:
!
!    16 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.  NORDER must be
!    between 1 and 10.
!
!    Output, real(kind=Rkind) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real(kind=Rkind) WEIGHT(NORDER), the weights of the rule.  WEIGHT(1) is
!    the weight at X = 1, WEIGHT(2) the weight at X = 0, and so on.
!    The weights are rational.  The weights are not symmetric, and
!    some weights may be negative.  They should sum to 1.
!
      integer norder
!
      real(kind=Rkind) d
      integer i
      real(kind=Rkind) weight(norder)
      real(kind=Rkind) xtab(norder)
!
      if ( norder .eq. 1 ) then

        weight(1) =  1.0 D+00

      else if ( norder .eq. 2 ) then

        d = 2.0 D+00

        weight(1) =  1.0 D+00 / d
        weight(2) =  1.0 D+00 / d

      else if ( norder .eq. 3 ) then

        d = 12.0 D+00

        weight(1) =    5.0 D+00 / d
        weight(2) =    8.0 D+00 / d
        weight(3) =  - 1.0 D+00 / d

      else if ( norder .eq. 4 ) then

        d = 24.0 D+00

        weight(1) =    9.0 D+00 / d
        weight(2) =   19.0 D+00 / d
        weight(3) =  - 5.0 D+00 / d
        weight(4) =    1.0 D+00 / d

      else if ( norder .eq. 5 ) then

        d = 720.0 D+00

        weight(1) =    251.0 D+00 / d
        weight(2) =    646.0 D+00 / d
        weight(3) =  - 264.0 D+00 / d
        weight(4) =    106.0 D+00 / d
        weight(5) =   - 19.0 D+00 / d

      else if ( norder .eq. 6 ) then

        d = 1440.0 D+00

        weight(1) =    475.0 D+00 / d
        weight(2) =   1427.0 D+00 / d
        weight(3) =  - 798.0 D+00 / d
        weight(4) =    482.0 D+00 / d
        weight(5) =  - 173.0 D+00 / d
        weight(6) =     27.0 D+00 / d

      else if ( norder .eq. 7 ) then

        d = 60480.0 D+00

        weight(1) =    19087.0 D+00 / d
        weight(2) =    65112.0 D+00 / d
        weight(3) =  - 46461.0 D+00 / d
        weight(4) =    37504.0 D+00 / d
        weight(5) =  - 20211.0 D+00 / d
        weight(6) =     6312.0 D+00 / d
        weight(7) =    - 863.0 D+00 / d

      else if ( norder .eq. 8 ) then

        d = 120960.0 D+00

        weight(1) =    36799.0 D+00 / d
        weight(2) =   139849.0 D+00 / d
        weight(3) = - 121797.0 D+00 / d
        weight(4) =   123133.0 D+00 / d
        weight(5) =  - 88547.0 D+00 / d
        weight(6) =    41499.0 D+00 / d
        weight(7) =  - 11351.0 D+00 / d
        weight(8) =     1375.0 D+00 / d

      else if ( norder .eq. 9 ) then

        d = 3628800.0 D+00

        weight(1) =   1070017.0 D+00 / d
        weight(2) =   4467094.0 D+00 / d
        weight(3) = - 4604594.0 D+00 / d
        weight(4) =   5595358.0 D+00 / d
        weight(5) = - 5033120.0 D+00 / d
        weight(6) =   3146338.0 D+00 / d
        weight(7) = - 1291214.0 D+00 / d
        weight(8) =    312874.0 D+00 / d
        weight(9) =   - 33953.0 D+00 / d

      else if ( norder .eq. 10 ) then

        d = 7257600.0 D+00

        weight(1) =    2082753.0 D+00 / d
        weight(2) =    9449717.0 D+00 / d
        weight(3) = - 11271304.0 D+00 / d
        weight(4) =   16002320.0 D+00 / d
        weight(5) = - 17283646.0 D+00 / d
        weight(6) =   13510082.0 D+00 / d
        weight(7) =  - 7394032.0 D+00 / d
        weight(8) =    2687864.0 D+00 / d
        weight(9) =   - 583435.0 D+00 / d
        weight(10) =     57281.0 D+00 / d

      else

        write ( *, * ) ' '
        write ( *, * ) 'MOULSET - Fatal error!'
        write ( *, * ) '  Illegal value of NORDER = ', norder
        write ( *, * ) '  Legal values are 1 through 10.'
        stop

      end if

      do i = 1, norder
        xtab(i) = dble ( 2 - i )
      end do

      return
      end subroutine moulset
      subroutine ncccom ( norder, xtab, weight )
!
!***********************************************************************
!
!! NCCCOM computes the coefficients of a Newton-Cotes closed quadrature rule.
!
!
!  Definition:
!
!    For the interval [-1,1], the Newton-Cotes open quadrature rule
!    estimates
!
!      Integral ( -1 <= X <= 1 ) F(X) dX
!
!    using NORDER equally spaced abscissas XTAB(I) and a weight vector
!    WEIGHT(I):
!
!      Sum ( I = 1 to N ) WEIGHT(I) * F ( XTAB(I) ).
!
!    For the CLOSED rule, the abscissas include A and B.
!
!  Modified:
!
!    25 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule, which should be
!    at least 2.
!
!    Output, real(kind=Rkind) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real(kind=Rkind) WEIGHT(NORDER), the weights of the rule.
!
      integer norder
!
      real(kind=Rkind) a
      real(kind=Rkind) b
      integer i
      real(kind=Rkind) weight(norder)
      real(kind=Rkind) xtab(norder)
!
!  Compute a closed quadrature rule.
!
      a = -1.0D+00
      b =  1.0D+00

      if ( norder .eq. 1 ) then

        xtab(1) = 0.0D+00
        weight(1) = 2.0D+00
        return

      else

        do i = 1, norder
          xtab(i) = ( dble ( norder - i ) * a + dble ( i - 1 ) * b )    &
            / dble ( norder - 1 )
        end do

      end if

      call nccom ( norder, a, b, xtab, weight )

      return
      end subroutine ncccom
      subroutine nccom ( norder, a, b, xtab, weight )
!
!***********************************************************************
!
!! NCCOM computes the coefficients of a Newton-Cotes quadrature rule.
!
!
!  Definition:
!
!    For the interval [A,B], the Newton-Cotes quadrature rule estimates
!
!      Integral ( A <= X <= B ) F(X) dX
!
!    using NORDER equally spaced abscissas XTAB(I) and a weight vector
!    WEIGHT(I):
!
!      Sum ( I = 1 to N ) WEIGHT(I) * F ( XTAB(I) ).
!
!    For the CLOSED rule, the abscissas include the points A and B.
!    For the OPEN rule, the abscissas do not include A and B.
!
!  Modified:
!
!    25 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!
!    Input, real(kind=Rkind) A, B, the left and right endpoints of the interval
!    over which the quadrature rule is to be applied.
!
!    Input, real(kind=Rkind) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real(kind=Rkind) WEIGHT(NORDER), the weights of the rule.
!
      integer maxorder
      parameter ( maxorder = 100 )
!
      integer norder
!
      real(kind=Rkind) a
      real(kind=Rkind) b
      real(kind=Rkind) diftab(maxorder)
      integer i
      integer j
      integer k
      real(kind=Rkind) weight(norder)
      real(kind=Rkind) xtab(norder)
      real(kind=Rkind) yvala
      real(kind=Rkind) yvalb
!
      if ( norder .gt. maxorder ) then
        write ( *, * ) ' '
        write ( *, * ) 'NCCOM - Fatal error!'
        write ( *, * ) '  Internal workspace vector is too small.'
        stop
      end if

      do i = 1, norder
!
!  Compute the Lagrange basis polynomial which is 1 at XTAB(I),
!  and zero at the other nodes.
!
        do j = 1, norder
          if ( j .eq. i ) then
            diftab(j) = 1.0 D+00
          else
            diftab(j) = 0.0 D+00
          end if
        end do

        do j = 2, norder
          do k = j, norder
            diftab(norder+j-k) =                                        &
              ( diftab(norder+j-k-1) - diftab(norder+j-k) )             &
              / ( xtab(norder+1-k) - xtab(norder+j-k) )
          end do
        end do

        do j = 1, norder-1
          do k = 1, norder-j
            diftab(norder-k) = diftab(norder-k)                         &
              - xtab(norder-k-j+1) * diftab(norder-k+1)
          end do
        end do
!
!  Evaluate the antiderivative of the polynomial at the left and
!  right endpoints.
!
        yvala = diftab(norder) / dble ( norder )
        do j = norder-1, 1, -1
          yvala = yvala * a + diftab(j) / dble ( j )
        end do
        yvala = yvala * a

        yvalb = diftab(norder) / dble ( norder )
        do j = norder-1, 1, -1
          yvalb = yvalb * b + diftab(j) / dble ( j )
        end do
        yvalb = yvalb * b

        weight(i) = yvalb - yvala

      end do

      return
      end subroutine nccom

      subroutine nccset ( norder, xtab, weight )
!
!***********************************************************************
!
!! NCCSET sets abscissas and weights for closed Newton-Cotes quadrature.
!
!
!  Integration interval:
!
!    [ -1, 1 ].
!
!  Weight function:
!
!    1.
!
!  Integral to approximate:
!
!    INTEGRAL ( -1 <= X <= 1 ) F(X) dX.
!
!  Approximate integral:
!
!    SUM ( I = 1 to NORDER ) WEIGHT(I) * F ( XTAB(I) ).
!
!  Note:
!
!    The closed Newton-Cotes rules use equally spaced abscissas, and
!    hence may be used with tabulated function data.
!
!    The rules are called "closed" because they include the endpoints.
!
!    The higher order rules involve negative weights.  These can produce
!    loss of accuracy due to the subtraction of large, nearly equal quantities.
!
!    NORDER = 2 is the trapezoidal rule.
!    NORDER = 3 is Simpson's rule.
!    NORDER = 4 is Simpson's 3/8 rule.
!    NORDER = 5 is Bode's rule.
!
!    The Kopal reference for NORDER = 12 lists
!      WEIGHT(6) = 15494566.0 / 43545600.0
!    but this results in a set of coeffients that don't add up to 2.
!    The correct value is
!      WEIGHT(6) = 15493566.0 / 43545600.0.
!
!  Reference:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964.
!
!    Johnson,
!    Quarterly Journal of Mathematics,
!    Volume 46, Number 52, 1915.
!
!    Zdenek Kopal,
!    Numerical Analysis,
!    John Wiley, 1955.
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Modified:
!
!    16 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!    NORDER must be between 2 and 20.
!
!    Output, real(kind=Rkind) XTAB(NORDER), the abscissas of the rule.
!    The abscissas are uniformly spaced in the interval, and include
!    -1 and 1.
!
!    Output, real(kind=Rkind) WEIGHT(NORDER), the weights of the rule.
!    The weights are symmetric, rational, and should sum to 2.0.
!    Some weights may be negative.
!
      integer norder
!
      real(kind=Rkind) d
      integer i
      real(kind=Rkind) weight(norder)
      real(kind=Rkind) xtab(norder)
!
      if ( norder .eq. 2 ) then

        weight(1) = 1.0 D+00
        weight(2) = 1.0 D+00

      else if ( norder .eq. 3 ) then

        d = 3.0 D+00

        weight(1) = 1.0 D+00 / d
        weight(2) = 4.0 D+00 / d
        weight(3) = 1.0 D+00 / d

      else if ( norder .eq. 4 ) then

        d = 4.0 D+00

        weight(1) = 1.0 D+00 / d
        weight(2) = 3.0 D+00 / d
        weight(3) = 3.0 D+00 / d
        weight(4) = 1.0 D+00 / d

      else if ( norder .eq. 5 ) then

        d = 45.0 D+00

        weight(1) =  7.0 D+00 / d
        weight(2) = 32.0 D+00 / d
        weight(3) = 12.0 D+00 / d
        weight(4) = 32.0 D+00 / d
        weight(5) =  7.0 D+00 / d

      else if ( norder .eq. 6 ) then

        d = 144.0 D+00

        weight(1) = 19.0 D+00 / d
        weight(2) = 75.0 D+00 / d
        weight(3) = 50.0 D+00 / d
        weight(4) = 50.0 D+00 / d
        weight(5) = 75.0 D+00 / d
        weight(6) = 19.0 D+00 / d

      else if ( norder .eq. 7 ) then

        d = 420.0 D+00

        weight(1) =  41.0 D+00 / d
        weight(2) = 216.0 D+00 / d
        weight(3) =  27.0 D+00 / d
        weight(4) = 272.0 D+00 / d
        weight(5) =  27.0 D+00 / d
        weight(6) = 216.0 D+00 / d
        weight(7) =  41.0 D+00 / d

      else if ( norder .eq. 8 ) then

        d = 8640.0 D+00

        weight(1) =  751.0 D+00 / d
        weight(2) = 3577.0 D+00 / d
        weight(3) = 1323.0 D+00 / d
        weight(4) = 2989.0 D+00 / d
        weight(5) = 2989.0 D+00 / d
        weight(6) = 1323.0 D+00 / d
        weight(7) = 3577.0 D+00 / d
        weight(8) =  751.0 D+00 / d

      else if ( norder .eq. 9 ) then

        d = 14175.0 D+00

        weight(1) =    989.0 D+00 / d
        weight(2) =   5888.0 D+00 / d
        weight(3) =  - 928.0 D+00 / d
        weight(4) =  10496.0 D+00 / d
        weight(5) = - 4540.0 D+00 / d
        weight(6) =  10496.0 D+00 / d
        weight(7) =  - 928.0 D+00 / d
        weight(8) =   5888.0 D+00 / d
        weight(9) =    989.0 D+00 / d

      else if ( norder .eq. 10 ) then

        d = 44800.0 D+00

        weight(1) =   2857.0 D+00 / d
        weight(2) =  15741.0 D+00 / d
        weight(3) =   1080.0 D+00 / d
        weight(4) =  19344.0 D+00 / d
        weight(5) =   5778.0 D+00 / d
        weight(6) =   5778.0 D+00 / d
        weight(7) =  19344.0 D+00 / d
        weight(8) =   1080.0 D+00 / d
        weight(9) =  15741.0 D+00 / d
        weight(10) =  2857.0 D+00 / d

      else if ( norder .eq. 11 ) then

        d = 299376.0 D+00

        weight(1) =     16067.0 D+00 / d
        weight(2) =    106300.0 D+00 / d
        weight(3) =   - 48525.0 D+00 / d
        weight(4) =    272400.0 D+00 / d
        weight(5) =  - 260550.0 D+00 / d
        weight(6) =    427368.0 D+00 / d
        weight(7) =  - 260550.0 D+00 / d
        weight(8) =    272400.0 D+00 / d
        weight(9) =   - 48525.0 D+00 / d
        weight(10) =   106300.0 D+00 / d
        weight(11) =    16067.0 D+00 / d

      else if ( norder .eq. 12 ) then

        d = 43545600.0 D+00

        weight(1) =     2171465.0 D+00 / d
        weight(2) =    13486539.0 D+00 / d
        weight(3) =   - 3237113.0 D+00 / d
        weight(4) =    25226685.0 D+00 / d
        weight(5) =   - 9595542.0 D+00 / d
        weight(6) =    15493566.0 D+00 / d
        weight(7) =    15493566.0 D+00 / d
        weight(8) =   - 9595542.0 D+00 / d
        weight(9) =    25226685.0 D+00 / d
        weight(10) =  - 3237113.0 D+00 / d
        weight(11) =   13486539.0 D+00 / d
        weight(12) =    2171465.0 D+00 / d

      else if ( norder .eq. 13 ) then

        d = 31531500.0 D+00

        weight(1) =      1364651.0 D+00 / d
        weight(2) =      9903168.0 D+00 / d
        weight(3) =    - 7587864.0 D+00 / d
        weight(4) =     35725120.0 D+00 / d
        weight(5) =   - 51491295.0 D+00 / d
        weight(6) =     87516288.0 D+00 / d
        weight(7) =   - 87797136.0 D+00 / d
        weight(8) =     87516288.0 D+00 / d
        weight(9) =   - 51491295.0 D+00 / d
        weight(10) =    35725120.0 D+00 / d
        weight(11) =   - 7587864.0 D+00 / d
        weight(12) =     9903168.0 D+00 / d
        weight(13) =     1364651.0 D+00 / d

      else if ( norder .eq. 14 ) then

        d = 150885504000.0 D+00

        weight(1) =      6137698213.0 D+00 / d
        weight(2) =     42194238652.0 D+00 / d
        weight(3) =   - 23361540993.0 D+00 / d
        weight(4) =    116778274403.0 D+00 / d
        weight(5) =  - 113219777650.0 D+00 / d
        weight(6) =    154424590209.0 D+00 / d
        weight(7) =   - 32067978834.0 D+00 / d
        weight(8) =   - 32067978834.0 D+00 / d
        weight(9) =    154424590209.0 D+00 / d
        weight(10) = - 113219777650.0 D+00 / d
        weight(11) =   116778274403.0 D+00 / d
        weight(12) =  - 23361540993.0 D+00 / d
        weight(13) =    42194238652.0 D+00 / d
        weight(14) =     6137698213.0 D+00 / d

      else if ( norder .eq. 15 ) then

        d = 2501928000.0 D+00

        weight(1) =       90241897.0 D+00 / d
        weight(2) =      710986864.0 D+00 / d
        weight(3) =    - 770720657.0 D+00 / d
        weight(4) =     3501442784.0 D+00 / d
        weight(5) =   - 6625093363.0 D+00 / d
        weight(6) =    12630121616.0 D+00 / d
        weight(7) =  - 16802270373.0 D+00 / d
        weight(8) =    19534438464.0 D+00 / d
        weight(9) =  - 16802270373.0 D+00 / d
        weight(10) =   12630121616.0 D+00 / d
        weight(11) =  - 6625093363.0 D+00 / d
        weight(12) =    3501442784.0 D+00 / d
        weight(13) =   - 770720657.0 D+00 / d
        weight(14) =     710986864.0 D+00 / d
        weight(15) =      90241897.0 D+00 / d

      else if ( norder .eq. 16 ) then

        d = 3099672576.0 D+00

        weight(1) =     105930069.0 D+00 / d
        weight(2) =     796661595.0 D+00 / d
        weight(3) =   - 698808195.0 D+00 / d
        weight(4) =    3143332755.0 D+00 / d
        weight(5) =  - 4688522055.0 D+00 / d
        weight(6) =    7385654007.0 D+00 / d
        weight(7) =  - 6000998415.0 D+00 / d
        weight(8) =    3056422815.0 D+00 / d
        weight(9) =    3056422815.0 D+00 / d
        weight(10) = - 6000998415.0 D+00 / d
        weight(11) =   7385654007.0 D+00 / d
        weight(12) = - 4688522055.0 D+00 / d
        weight(13) =   3143332755.0 D+00 / d
        weight(14) =  - 698808195.0 D+00 / d
        weight(15) =    796661595.0 D+00 / d
        weight(16) =    105930069.0 D+00 / d

      else if ( norder .eq. 17 ) then

        d = 488462349375.0 D+00

        weight(1) =       15043611773.0 D+00 / d
        weight(2) =      127626606592.0 D+00 / d
        weight(3) =    - 179731134720.0 D+00 / d
        weight(4) =      832211855360.0 D+00 / d
        weight(5) =   - 1929498607520.0 D+00 / d
        weight(6) =     4177588893696.0 D+00 / d
        weight(7) =   - 6806534407936.0 D+00 / d
        weight(8) =     9368875018240.0 D+00 / d
        weight(9) =  - 10234238972220.0 D+00 / d
        weight(10) =    9368875018240.0 D+00 / d
        weight(11) =  - 6806534407936.0 D+00 / d
        weight(12) =    4177588893696.0 D+00 / d
        weight(13) =  - 1929498607520.0 D+00 / d
        weight(14) =     832211855360.0 D+00 / d
        weight(15) =   - 179731134720.0 D+00 / d
        weight(16) =     127626606592.0 D+00 / d
        weight(17) =      15043611773.0 D+00 / d

      else if ( norder .eq. 18 ) then

        d = 1883051089920000.0 D+00

        weight(1) =       55294720874657.0 D+00 / d
        weight(2) =      450185515446285.0 D+00 / d
        weight(3) =    - 542023437008852.0 D+00 / d
        weight(4) =     2428636525764260.0 D+00 / d
        weight(5) =   - 4768916800123440.0 D+00 / d
        weight(6) =     8855416648684984.0 D+00 / d
        weight(7) =  - 10905371859796660.0 D+00 / d
        weight(8) =    10069615750132836.0 D+00 / d
        weight(9) =   - 3759785974054070.0 D+00 / d
        weight(10) =  - 3759785974054070.0 D+00 / d
        weight(11) =   10069615750132836.0 D+00 / d
        weight(12) = - 10905371859796660.0 D+00 / d
        weight(13) =    8855416648684984.0 D+00 / d
        weight(14) =  - 4768916800123440.0 D+00 / d
        weight(15) =    2428636525764260.0 D+00 / d
        weight(16) =   - 542023437008852.0 D+00 / d
        weight(17) =     450185515446285.0 D+00 / d
        weight(18) =      55294720874657.0 D+00 / d

      else if ( norder .eq. 19 ) then

        d = 7604556960000.0 D+00

        weight(1) =       203732352169.0 D+00 / d
        weight(2) =      1848730221900.0 D+00 / d
        weight(3) =    - 3212744374395.0 D+00 / d
        weight(4) =     15529830312096.0 D+00 / d
        weight(5) =   - 42368630685840.0 D+00 / d
        weight(6) =    103680563465808.0 D+00 / d
        weight(7) =  - 198648429867720.0 D+00 / d
        weight(8) =    319035784479840.0 D+00 / d
        weight(9) =  - 419127951114198.0 D+00 / d
        weight(10) =   461327344340680.0 D+00 / d
        weight(11) = - 419127951114198.0 D+00 / d
        weight(12) =   319035784479840.0 D+00 / d
        weight(13) = - 198648429867720.0 D+00 / d
        weight(14) =   103680563465808.0 D+00 / d
        weight(15) =  - 42368630685840.0 D+00 / d
        weight(16) =    15529830312096.0 D+00 / d
        weight(17) =   - 3212744374395.0 D+00 / d
        weight(18) =     1848730221900.0 D+00 / d
        weight(19) =      203732352169.0 D+00 / d

      else if ( norder .eq. 20 ) then

        d = 2688996956405760000.0 D+00

        weight(1) =       69028763155644023.0 D+00 / d
        weight(2) =      603652082270808125.0 D+00 / d
        weight(3) =    - 926840515700222955.0 D+00 / d
        weight(4) =     4301581538450500095.0 D+00 / d
        weight(5) =  - 10343692234243192788.0 D+00 / d
        weight(6) =    22336420328479961316.0 D+00 / d
        weight(7) =  - 35331888421114781580.0 D+00 / d
        weight(8) =    43920768370565135580.0 D+00 / d
        weight(9) =  - 37088370261379851390.0 D+00 / d
        weight(10) =   15148337305921759574.0 D+00 / d
        weight(11) =   15148337305921759574.0 D+00 / d
        weight(12) = - 37088370261379851390.0 D+00 / d
        weight(13) =   43920768370565135580.0 D+00 / d
        weight(14) = - 35331888421114781580.0 D+00 / d
        weight(15) =   22336420328479961316.0 D+00 / d
        weight(16) = - 10343692234243192788.0 D+00 / d
        weight(17) =    4301581538450500095.0 D+00 / d
        weight(18) =   - 926840515700222955.0 D+00 / d
        weight(19) =     603652082270808125.0 D+00 / d
        weight(20) =      69028763155644023.0 D+00 / d

      else

        write ( *, * ) ' '
        write ( *, * ) 'NCCSET - Fatal error!'
        write ( *, * ) '  Illegal value of NORDER = ', norder
        write ( *, * ) '  Legal values are 2 through 20.'
        stop

      end if
!
!  The abscissas are uniformly spaced.
!
      do i = 1, norder
        xtab(i) = dble ( 2 * i - 1 - norder ) / dble ( norder - 1 )
      end do

      return
      end subroutine nccset
      subroutine ncocom ( norder, xtab, weight )
!
!***********************************************************************
!
!! NCOCOM computes the coefficients of a Newton-Cotes open quadrature rule.
!
!
!  Definition:
!
!    For the interval [-1,1], the Newton-Cotes open quadrature rule
!    estimates
!
!      Integral ( -1 <= X <= 1 ) F(X) dX
!
!    using NORDER equally spaced abscissas XTAB(I) and a weight vector
!    WEIGHT(I):
!
!      Sum ( I = 1 to N ) WEIGHT(I) * F ( XTAB(I) ).
!
!    For the OPEN rule, the abscissas do not include A and B.
!
!  Modified:
!
!    25 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!
!    Output, real(kind=Rkind) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real(kind=Rkind) WEIGHT(NORDER), the weights of the  rule.
!
      integer norder
!
      real(kind=Rkind) a
      real(kind=Rkind) b
      integer i
      real(kind=Rkind) weight(norder)
      real(kind=Rkind) xtab(norder)
!
      a = -1.0D+00
      b =  1.0D+00

      do i = 1, norder
        xtab(i) = ( dble ( norder + 1 - i ) * a + dble ( i ) * b )      &
          / dble ( norder + 1 )
      end do

      call nccom ( norder, a, b, xtab, weight )

      return
      end subroutine ncocom
      subroutine ncoset ( norder, xtab, weight )
!
!***********************************************************************
!
!! NCOSET sets abscissas and weights for open Newton-Cotes quadrature.
!
!
!  Integration interval:
!
!    [ -1, 1 ].
!
!  Weight function:
!
!    1.
!
!  Integral to approximate:
!
!    INTEGRAL ( -1 <= X <= 1 ) F(X) dX.
!
!  Approximate integral:
!
!    SUM ( I = 1 to NORDER ) WEIGHT(I) * F ( XTAB(I) ).
!
!  Note:
!
!    The open Newton-Cotes rules use equally spaced abscissas, and
!    hence may be used with equally spaced data.
!
!    The rules are called "open" because they do not include the interval
!    endpoints.
!
!    Most of the rules involve negative weights.  These can produce loss
!    of accuracy due to the subtraction of large, nearly equal quantities.
!
!  Reference:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964.
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Modified:
!
!    15 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!    NORDER must be between 1 and 7, and 9.
!
!    Output, real(kind=Rkind) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real(kind=Rkind) WEIGHT(NORDER), the weights of the rule.
!    The weights are rational, symmetric, and should sum to 2.
!    Some weights may be negative.
!
      integer norder
!
      real(kind=Rkind) d
      integer i
      real(kind=Rkind) weight(norder)
      real(kind=Rkind) xtab(norder)
!
      if ( norder .eq. 1 ) then

        weight(1) = 2.0 D+00

      else if ( norder .eq. 2 ) then

        weight(1) = 1.0 D+00
        weight(2) = 1.0 D+00

      else if ( norder .eq. 3 ) then

        d = 3.0 D+00

        weight(1) =   4.0 D+00 / d
        weight(2) = - 2.0 D+00 / d
        weight(3) =   4.0 D+00 / d

      else if ( norder .eq. 4 ) then

        d = 12.0 D+00

        weight(1) = 11.0 D+00 / d
        weight(2) =  1.0 D+00 / d
        weight(3) =  1.0 D+00 / d
        weight(4) = 11.0 D+00 / d

      else if ( norder .eq. 5 ) then

        d = 10.0 D+00

        weight(1) =   11.0 D+00 / d
        weight(2) = - 14.0 D+00 / d
        weight(3) =   26.0 D+00 / d
        weight(4) = - 14.0 D+00 / d
        weight(5) =   11.0 D+00 / d

      else if ( norder .eq. 6 ) then

        d = 1440.0 D+00

        weight(1) =  1222.0 D+00 / d
        weight(2) = - 906.0 D+00 / d
        weight(3) =  1124.0 D+00 / d
        weight(4) =  1124.0 D+00 / d
        weight(5) = - 906.0 D+00 / d
        weight(6) =  1222.0 D+00 / d

      else if ( norder .eq. 7 ) then

        d = 945.0 D+00

        weight(1) =    920.0 D+00 / d
        weight(2) = - 1908.0 D+00 / d
        weight(3) =   4392.0 D+00 / d
        weight(4) = - 4918.0 D+00 / d
        weight(5) =   4392.0 D+00 / d
        weight(6) = - 1908.0 D+00 / d
        weight(7) =    920.0 D+00 / d

      else if ( norder .eq. 9 ) then

        d = 4536.0 D+00

        weight(1) =    4045.0 D+00 / d
        weight(2) = - 11690.0 D+00 / d
        weight(3) =   33340.0 D+00 / d
        weight(4) = - 55070.0 D+00 / d
        weight(5) =   67822.0 D+00 / d
        weight(6) = - 55070.0 D+00 / d
        weight(7) =   33340.0 D+00 / d
        weight(8) = - 11690.0 D+00 / d
        weight(9) =    4045.0 D+00 / d

      else

        write ( *, * ) ' '
        write ( *, * ) 'NCOSET - Fatal error!'
        write ( *, * ) '  Illegal value of NORDER = ', norder
        write ( *, * ) '  Legal values are 1 to 7, and 9.'
        stop

      end if
!
!  Set the abscissas.
!
      do i = 1, norder
        xtab(i) = dble ( 2 * i - norder - 1 ) / dble ( norder + 1 )
      end do

      return
      end subroutine ncoset
      subroutine radset ( norder, xtab, weight )
!
!***********************************************************************
!
!! RADSET sets abscissas and weights for Radau quadrature.
!
!
!  Integration interval:
!
!    [ -1, 1 ].
!
!  Weight function:
!
!    1.
!
!  Integral to approximate:
!
!    INTEGRAL ( -1 <= X <= 1 ) F(X) dX.
!
!  Approximate integral:
!
!    SUM ( I = 1 to NORDER ) WEIGHT(I) * F ( XTAB(I) ).
!
!  Precision:
!
!    The quadrature rule will integrate exactly all polynomials up to
!    X**(2*NORDER-2).
!
!  Note:
!
!    The Radau rule is distinguished by the fact that the left endpoint
!    (-1) is always an abscissa.
!
!  Reference:
!
!    Abramowitz and Stegun,
!    Handbook of Mathematical Functions,
!    National Bureau of Standards, 1964.
!
!    Arthur Stroud and Don Secrest,
!    Gaussian Quadrature Formulas,
!    Prentice Hall, 1966.
!
!    Daniel Zwillinger, editor,
!    Standard Mathematical Tables and Formulae,
!    30th Edition,
!    CRC Press, 1996.
!
!  Modified:
!
!    16 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer NORDER, the order of the rule.
!    NORDER must be between 1 and 15.
!
!    Output, real(kind=Rkind) XTAB(NORDER), the abscissas of the rule.
!
!    Output, real(kind=Rkind) WEIGHT(NORDER), the weights of the rule.
!    The weights are positive.  The weights are not symmetric.
!    The weights should sum to 2.  WEIGHT(1) should equal 2 / NORDER**2.
!
      integer norder
!
      real(kind=Rkind) xtab(norder)
      real(kind=Rkind) weight(norder)
!
      if ( norder .eq. 1 ) then

        xtab(1) =   - 1.0 D+00
        weight(1) =   2.0 D+00

      else if ( norder .eq. 2 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =    1.0 D+00 / 3.0 D+00

        weight(1) =  0.5 D+00
        weight(2) =  1.5 D+00

      else if ( norder .eq. 3 ) then

        xtab(1) =   - 1.0 D+00
        xtab(2) =   - 0.289897948556635619639456814941 D+00
        xtab(3) =     0.689897948556635619639456814941 D+00

        weight(1) =  0.222222222222222222222222222222 D+00
        weight(2) =  0.102497165237684322767762689304 D+01
        weight(3) =  0.752806125400934550100150884739 D+00

      else if ( norder .eq. 4 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =  - 0.575318923521694112050483779752 D+00
        xtab(3) =    0.181066271118530578270147495862 D+00
        xtab(4) =    0.822824080974592105208907712461 D+00

        weight(1) =  0.125 D+00
        weight(2) =  0.657688639960119487888578442146 D+00
        weight(3) =  0.776386937686343761560464613780 D+00
        weight(4) =  0.440924422353536750550956944074 D+00

      else if ( norder .eq. 5 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =  - 0.720480271312438895695825837750 D+00
        xtab(3) =  - 0.167180864737833640113395337326 D+00
        xtab(4) =    0.446313972723752344639908004629 D+00
        xtab(5) =    0.885791607770964635613757614892 D+00

        weight(1) =  0.08 D+00
        weight(2) =  0.446207802167141488805120436457 D+00
        weight(3) =  0.623653045951482508163709823153 D+00
        weight(4) =  0.562712030298924120384345300681 D+00
        weight(5) =  0.287427121582451882646824439708 D+00

      else if ( norder .eq. 6 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =  - 0.802929828402347147753002204224 D+00
        xtab(3) =  - 0.390928546707272189029229647442 D+00
        xtab(4) =    0.124050379505227711989974959990 D+00
        xtab(5) =    0.603973164252783654928415726409 D+00
        xtab(6) =    0.920380285897062515318386619813 D+00

        weight(1) =  0.555555555555555555555555555556 D-01
        weight(2) =  0.319640753220510966545779983796 D+00
        weight(3) =  0.485387188468969916159827915587 D+00
        weight(4) =  0.520926783189574982570229406570 D+00
        weight(5) =  0.416901334311907738959406382743 D+00
        weight(6) =  0.201588385253480840209200755749 D+00

      else if ( norder .eq. 7 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =  - 0.853891342639482229703747931639 D+00
        xtab(3) =  - 0.538467724060109001833766720231 D+00
        xtab(4) =  - 0.117343037543100264162786683611 D+00
        xtab(5) =    0.326030619437691401805894055838 D+00
        xtab(6) =    0.703842800663031416300046295008 D+00
        xtab(7) =    0.941367145680430216055899446174 D+00

        weight(1) =  0.408163265306122448979591836735 D-01
        weight(2) =  0.239227489225312405787077480770 D+00
        weight(3) =  0.380949873644231153805938347876 D+00
        weight(4) =  0.447109829014566469499348953642 D+00
        weight(5) =  0.424703779005955608398308039150 D+00
        weight(6) =  0.318204231467301481744870434470 D+00
        weight(7) =  0.148988471112020635866497560418 D+00

      else if ( norder .eq. 8 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =  - 0.887474878926155707068695617935 D+00
        xtab(3) =  - 0.639518616526215270024840114382 D+00
        xtab(4) =  - 0.294750565773660725252184459658 D+00
        xtab(5) =    0.943072526611107660028971153047 D-01
        xtab(6) =    0.468420354430821063046421216613 D+00
        xtab(7) =    0.770641893678191536180719525865 D+00
        xtab(8) =    0.955041227122575003782349000858 D+00

        weight(1) =  0.03125 D+00
        weight(2) =  0.185358154802979278540728972699 D+00
        weight(3) =  0.304130620646785128975743291400 D+00
        weight(4) =  0.376517545389118556572129261442 D+00
        weight(5) =  0.391572167452493593082499534004 D+00
        weight(6) =  0.347014795634501280228675918422 D+00
        weight(7) =  0.249647901329864963257869293513 D+00
        weight(8) =  0.114508814744257199342353728520 D+00

      else if ( norder .eq. 9 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =  - 0.910732089420060298533757956283 D+00
        xtab(3) =  - 0.711267485915708857029562959544 D+00
        xtab(4) =  - 0.426350485711138962102627520502 D+00
        xtab(5) =  - 0.903733696068532980645444599064 D-01
        xtab(6) =    0.256135670833455395138292079035 D+00
        xtab(7) =    0.571383041208738483284917464837 D+00
        xtab(8) =    0.817352784200412087992517083851 D+00
        xtab(9) =    0.964440169705273096373589797925 D+00

        weight(1) =  0.246913580246913580246913580247 D-01
        weight(2) =  0.147654019046315385819588499802 D+00
        weight(3) =  0.247189378204593052361239794969 D+00
        weight(4) =  0.316843775670437978338000849642 D+00
        weight(5) =  0.348273002772966594071991031186 D+00
        weight(6) =  0.337693966975929585803724239792 D+00
        weight(7) =  0.286386696357231171146705637752 D+00
        weight(8) =  0.200553298024551957421165090417 D+00
        weight(9) =  0.907145049232829170128934984159 D-01

      else if ( norder .eq. 10 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =  - 0.927484374233581078117671398464 D+00
        xtab(3) =  - 0.763842042420002599615429776011 D+00
        xtab(4) =  - 0.525646030370079229365386614293 D+00
        xtab(5) =  - 0.236234469390588049278459503207 D+00
        xtab(6) =    0.760591978379781302337137826389 D-01
        xtab(7) =    0.380664840144724365880759065541 D+00
        xtab(8) =    0.647766687674009436273648507855 D+00
        xtab(9) =    0.851225220581607910728163628088 D+00
        xtab(10) =   0.971175180702246902734346518378 D+00

        weight(1) =  0.02 D+00
        weight(2) =  0.120296670557481631517310522702 D+00
        weight(3) =  0.204270131879000675555788672223 D+00
        weight(4) =  0.268194837841178696058554475262 D+00
        weight(5) =  0.305859287724422621016275475401 D+00
        weight(6) =  0.313582457226938376695902847302 D+00
        weight(7) =  0.290610164832918311146863077963 D+00
        weight(8) =  0.239193431714379713376571966160 D+00
        weight(9) =  0.164376012736921475701681668908 D+00
        weight(10) = 0.736170054867584989310512940790 D-01

      else if ( norder .eq. 11 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =  - 0.939941935677027005913871284731 D+00
        xtab(3) =  - 0.803421975580293540697597956820 D+00
        xtab(4) =  - 0.601957842073797690275892603234 D+00
        xtab(5) =  - 0.351888923353330214714301017870 D+00
        xtab(6) =  - 0.734775314313212657461903554238 D-01
        xtab(7) =    0.210720306228426314076095789845 D+00
        xtab(8) =    0.477680647983087519467896683890 D+00
        xtab(9) =    0.705777100713859519144801128840 D+00
        xtab(10) =   0.876535856245703748954741265611 D+00
        xtab(11) =   0.976164773135168806180508826082 D+00

        weight(1) =  0.165289256198347107438016528926 D-01
        weight(2) =  0.998460819079680638957534695802 D-01
        weight(3) =  0.171317619206659836486712649042 D+00
        weight(4) =  0.228866123848976624401683231126 D+00
        weight(5) =  0.267867086189684177806638163355 D+00
        weight(6) =  0.285165563941007337460004408915 D+00
        weight(7) =  0.279361333103383045188962195720 D+00
        weight(8) =  0.250925377697128394649140267633 D+00
        weight(9) =  0.202163108540024418349931754266 D+00
        weight(10) = 0.137033682133202256310153880580 D+00
        weight(11) = 0.609250978121311347072183268883 D-01

      else if ( norder .eq. 12 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =  - 0.949452759204959300493337627077 D+00
        xtab(3) =  - 0.833916773105189706586269254036 D+00
        xtab(4) =  - 0.661649799245637148061133087811 D+00
        xtab(5) =  - 0.444406569781935851126642615609 D+00
        xtab(6) =  - 0.196994559534278366455441427346 D+00
        xtab(7) =    0.637247738208319158337792384845 D-01
        xtab(8) =    0.319983684170669623532789532206 D+00
        xtab(9) =    0.554318785912324288984337093085 D+00
        xtab(10) =   0.750761549711113852529400825472 D+00
        xtab(11) =   0.895929097745638894832914608454 D+00
        xtab(12) =   0.979963439076639188313950540264 D+00

        weight(1) =  0.138888888888888888888888888888 D-01
        weight(2) =  0.841721349386809762415796536813 D-01
        weight(3) =  0.145563668853995128522547654706 D+00
        weight(4) =  0.196998534826089634656049637969 D+00
        weight(5) =  0.235003115144985839348633985940 D+00
        weight(6) =  0.256991338152707776127974253598 D+00
        weight(7) =  0.261465660552133103438074715743 D+00
        weight(8) =  0.248121560804009959403073107079 D+00
        weight(9) =  0.217868879026192438848747482023 D+00
        weight(10) = 0.172770639313308564306065766966 D+00
        weight(11) = 0.115907480291738392750341908272 D+00
        weight(12) = 0.512480992072692974680229451351 D-01

      else if ( norder .eq. 13 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =  - 0.956875873668299278183813833834 D+00
        xtab(3) =  - 0.857884202528822035697620310269 D+00
        xtab(4) =  - 0.709105087529871761580423832811 D+00
        xtab(5) =  - 0.519197779050454107485205148087 D+00
        xtab(6) =  - 0.299201300554509985532583446686 D+00
        xtab(7) =  - 0.619016986256353412578604857936 D-01
        xtab(8) =    0.178909837597084635021931298881 D+00
        xtab(9) =    0.409238231474839556754166331248 D+00
        xtab(10) =   0.615697890940291918017885487543 D+00
        xtab(11) =   0.786291018233046684731786459135 D+00
        xtab(12) =   0.911107073689184553949066402429 D+00
        xtab(13) =   0.982921890023145161262671078244 D+00

        weight(1) =  0.118343195266272189349112426036 D-01
        weight(2) =  0.719024162924955289397537405641 D-01
        weight(3) =  0.125103834331152358133769287976 D+00
        weight(4) =  0.171003460470616642463758674512 D+00
        weight(5) =  0.206960611455877074631132560829 D+00
        weight(6) =  0.230888862886995434012203758668 D+00
        weight(7) =  0.241398342287691148630866924129 D+00
        weight(8) =  0.237878547660712031342685189180 D+00
        weight(9) =  0.220534229288451464691077164199 D+00
        weight(10) = 0.190373715559631732254759820746 D+00
        weight(11) = 0.149150950090000205151491864242 D+00
        weight(12) = 0.992678068818470859847363877478 D-01
        weight(13) = 0.437029032679020748288533846051 D-01

      else if ( norder .eq. 14 ) then

        xtab(1) =  - 1.0 D+00
        xtab(2) =  - 0.962779269978024297120561244319 D+00
        xtab(3) =  - 0.877048918201462024795266773531 D+00
        xtab(4) =  - 0.747389642613378838735429134263 D+00
        xtab(5) =  - 0.580314056546874971105726664999 D+00
        xtab(6) =  - 0.384202003439203313794083903375 D+00
        xtab(7) =  - 0.168887928042680911008441695622 D+00
        xtab(8) =    0.548312279917645496498107146428 D-01
        xtab(9) =    0.275737205435522399182637403545 D+00
        xtab(10) =   0.482752918588474966820418534355 D+00
        xtab(11) =   0.665497977216884537008955042481 D+00
        xtab(12) =   0.814809550601994729434217249123 D+00
        xtab(13) =   0.923203722520643299246334950272 D+00
        xtab(14) =   0.985270697947821356698617003172 D+00

        weight(1) =  0.102040816326530612244897959184 D-01
        weight(2) =  0.621220169077714601661329164668 D-01
        weight(3) =  0.108607722744362826826720935229 D+00
        weight(4) =  0.149620539353121355950520836946 D+00
        weight(5) =  0.183127002125729654123867302103 D+00
        weight(6) =  0.207449763335175672668082886489 D+00
        weight(7) =  0.221369811499570948931671683021 D+00
        weight(8) =  0.224189348002707794238414632220 D+00
        weight(9) =  0.215767100604618851381187446115 D+00
        weight(10) = 0.196525518452982430324613091930 D+00
        weight(11) = 0.167429727891086278990102277038 D+00
        weight(12) = 0.129939668737342347807425737146 D+00
        weight(13) = 0.859405354429804030893077310866 D-01
        weight(14) = 0.377071632698969142774627282919 D-01

      else if ( norder .eq. 15 ) then

        xtab(1) =  - 1.0
        xtab(2) =  - 0.967550468197200476562456018282 D+00
        xtab(3) =  - 0.892605400120550767066811886849 D+00
        xtab(4) =  - 0.778685617639031079381743321893 D+00
        xtab(5) =  - 0.630779478886949283946148437224 D+00
        xtab(6) =  - 0.455352905778529370872053455981 D+00
        xtab(7) =  - 0.260073376740807915768961188263 D+00
        xtab(8) =  - 0.534757226797460641074538896258 D-01
        xtab(9) =    0.155410685384859484319182024964 D+00
        xtab(10) =   0.357456512022127651195319205174 D+00
        xtab(11) =   0.543831458701484016930711802760 D+00
        xtab(12) =   0.706390264637572540152679669478 D+00
        xtab(13) =   0.838029000636089631215097384520 D+00
        xtab(14) =   0.932997190935973719928072142859 D+00
        xtab(15) =   0.987166478414363086378359071811 D+00

        weight(1) =  0.888888888888888888888888888889 D-02
        weight(2) =  0.542027800486444943382142368018 D-01
        weight(3) =  0.951295994604808992038477266346 D-01
        weight(4) =  0.131875462504951632186262157944 D+00
        weight(5) =  0.162854477303832629448732245828 D+00
        weight(6) =  0.186715145839450908083795103799 D+00
        weight(7) =  0.202415187030618429872703310435 D+00
        weight(8) =  0.209268608147694581430889790306 D+00
        weight(9) =  0.206975960249553755479027321787 D+00
        weight(10) = 0.195637503045116116473556617575 D+00
        weight(11) = 0.175748872642447685670310440476 D+00
        weight(12) = 0.148179527003467253924682058743 D+00
        weight(13) = 0.114135203489752753013075582569 D+00
        weight(14) = 0.751083927605064397329716653914 D-01
        weight(15) = 0.328643915845935322530428528231 D-01

      else

        write ( *, * ) ' '
        write ( *, * ) 'RADSET - Fatal error!'
        write ( *, * ) '  Illegal value of NORDER = ', norder
        write ( *, * ) '  Legal values are 1 to 15.'
        stop

      end if

      return
      end subroutine radset
      subroutine summer ( func, norder, xtab, weight, result )
!
!***********************************************************************
!
!! SUMMER carries out a quadrature rule over a single interval.
!
!
!  Formula:
!
!    RESULT = SUM ( I = 1 to NORDER ) WEIGHT(I) * FUNC ( XTAB(I) ).
!
!  Modified:
!
!    16 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FUNC, the name of the FORTRAN function which
!    evaluates the integrand.  The function must have the form
!      real(kind=Rkind) func ( x ).
!
!    Input, integer NORDER, the order of the rule.
!
!    Input, real(kind=Rkind) XTAB(NORDER), the abscissas of the rule.
!
!    Input, real(kind=Rkind) WEIGHT(NORDER), the weights of the rule.
!
!    Output, real(kind=Rkind) RESULT, the approximate value of the integral.
!
      integer norder
!
      real(kind=Rkind) func
      integer i
      real(kind=Rkind) result
      real(kind=Rkind) weight(norder)
      real(kind=Rkind) xtab(norder)
!
      external func
!
      if ( norder .lt. 1 ) then
        write ( *, * ) ' '
        write ( *, * ) 'SUMMER - Fatal error!'
        write ( *, * ) '  NORDER must be at least 1.'
        write ( *, * ) '  The input value was NORDER = ', norder
        stop
      end if

      result = 0.0 D+00
      do i = 1, norder
        result = result + weight(i) * func ( xtab(i) )
      end do

      return
      end subroutine summer
      subroutine summergk ( func, norderg, weightg, resultg, norderk,   &
        xtabk, weightk, resultk )
!
!***********************************************************************
!
!! SUMMERGK carries out Gauss-Kronrod quadrature over a single interval.
!
!
!  Note:
!
!    The abscissas for the Gauss-Legendre rule of order NORDERG are
!    not required, since they are assumed to be the even-indexed
!    entries of the corresponding Kronrod rule.
!
!  Modified:
!
!    16 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FUNC, the name of the FORTRAN function which
!    evaluates the integrand.  The function must have the form
!      real(kind=Rkind) func ( x ).
!
!    Input, integer NORDERG, the order of the Gauss-Legendre rule.
!
!    Input, real(kind=Rkind) WEIGHTG(NORDERG), the weights of the
!    Gauss-Legendre rule.
!
!    Output, real(kind=Rkind) RESULTG, the approximate value of the
!    integral, based on the Gauss-Legendre rule.
!
!    Input, integer NORDERK, the order of the Kronrod rule.  NORDERK
!    must equal 2 * NORDERG + 1.
!
!    Input, real(kind=Rkind) XTABK(NORDERK), the abscissas of the Kronrod rule.
!
!    Input, real(kind=Rkind) WEIGHTK(NORDERK), the weights of the Kronrod rule.
!
!    Output, real(kind=Rkind) RESULTK, the approximate value of the integral,
!    based on the Kronrod rule.
!
      integer norderg
      integer norderk
!
      real(kind=Rkind) fk
      real(kind=Rkind) func
      integer i
      real(kind=Rkind) resultg
      real(kind=Rkind) resultk
      real(kind=Rkind) weightg(norderg)
      real(kind=Rkind) weightk(norderk)
      real(kind=Rkind) xtabk(norderk)
!
      external func
!
      if ( norderk .ne. 2 * norderg + 1 ) then
        write ( *, * ) ' '
        write ( *, * ) 'SUMMERGK - Fatal error!'
        write ( *, * ) '  NORDERK must equal 2 * NORDERG + 1.'
        write ( *, * ) '  The input value was NORDERG = ', norderg
        write ( *, * ) '  The input value was NORDERK = ', norderk
        stop
      end if

      resultg = 0.0 D+00
      resultk = 0.0 D+00

      do i = 1, norderk

        fk = func ( xtabk(i) )

        resultk = resultk + weightk(i) * fk

        if ( mod ( i, 2 ) .eq. 0 )then
          resultg = resultg + weightg(i/2) * fk
        end if

      end do

      return
      end subroutine summergk
      subroutine sumsub ( func, a, b, nsub, norder, xtab, weight,       &
        result )
!
!***********************************************************************
!
!! SUMSUB carries out a composite quadrature rule.
!
!
!  Integration interval:
!
!    [ A, B ]
!
!  Integral to approximate:
!
!    INTEGRAL ( A <= X <= B ) F(X) dX.
!
!  Approximate integral:
!
!    H = ( B - A ) / NSUB;
!    XMID(J) = A + 0.5 * H * REAL ( 2 * J - 1 );
!
!    SUM ( J = 1 to NSUB )
!      SUM ( I = 1 to NORDER )
!        WEIGHT(I) * F ( XMID(J) + 0.5 * H * XTAB(I) ).
!
!  Note:
!
!    The quadrature rule weights and abscissas are assumed to be
!    defined for the interval [-1,1].  These values may be computed
!    by any of the routines CHEBSET, CLOSET, LEGCOM, LEGSET, LOBSET,
!    OPNSET, RADSET.
!
!  Modified:
!
!    16 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FUNC, the name of the FORTRAN function which
!    evaluates the integrand.  The function must have the form
!      real(kind=Rkind) func ( x ).
!
!    Input, real(kind=Rkind) A, B, the lower and upper limits of integration.
!
!    Input, integer NSUB, the number of equal subintervals into
!    which the finite interval (A,B) is to be subdivided for
!    higher accuracy.  NSUB must be at least 1.
!
!    Input, integer NORDER, the order of the rule.
!    NORDER must be at least 1.
!
!    Input, real(kind=Rkind) XTAB(NORDER), the abscissas of the rule.
!
!    Input, real(kind=Rkind) WEIGHT(NORDER), the weights of the rule.
!
!    Output, real(kind=Rkind) RESULT, the approximate value of the integral.
!
      integer norder
!
      real(kind=Rkind) a
      real(kind=Rkind) b
      real(kind=Rkind) func
      real(kind=Rkind) h
      integer i
      integer j
      integer nsub
      real(kind=Rkind) result
      real(kind=Rkind) x
      real(kind=Rkind) xmid
      real(kind=Rkind) xtab(norder)
      real(kind=Rkind) weight(norder)
!
      external func
!
      result = 0.0 D+00

      if ( a .eq. b ) then
        return
      end if

      if ( norder .lt. 1 ) then
        write ( *, * ) ' '
        write ( *, * ) 'SUMSUB - Fatal error!'
        write ( *, * ) '  Nonpositive value of NORDER = ', norder
        stop
      end if

      if ( nsub .lt. 1 ) then
        write ( *, * ) ' '
        write ( *, * ) 'SUMSUB - Fatal error!'
        write ( *, * ) '  Nonpositive value of NSUB = ', nsub
        stop
      end if

      h = ( b - a ) / dble ( nsub )

      do j = 1, nsub
        xmid = a + 0.5D+00 * h * dble ( 2 * j - 1 )
        do i = 1, norder
          x = xmid + 0.5D+00 * h * xtab(i)
          result = result + weight(i) * func ( x )
        end do
      end do

      result = 0.5D+00 * h * result

      return
      end subroutine sumsub
      subroutine sumsubgk ( func, a, b, nsub, norderg, weightg,         &
        resultg, norderk, xtabk, weightk, resultk, error )
!
!***********************************************************************
!
!! SUMSUBGK carries out a composite Gauss-Kronrod rule.
!
!
!  Integration interval:
!
!    [ A, B ]
!
!  Integral to approximate:
!
!    INTEGRAL ( A <= X <= B ) F(X) dX.
!
!  Approximate integral:
!
!    H = ( B - A ) / NSUB;
!    XMID(J) = A + 0.5 * H * REAL ( 2 * J - 1 );
!
!    SUM ( J = 1 to NSUB )
!      SUM ( I = 1 to NORDERK )
!        WEIGHTK(I) * F ( XMID(J) + 0.5 * H * XTABK(I) ).
!
!  Note:
!
!    The Gauss-Legendre weights should be computed by LEGCOM or LEGSET.
!    The Kronrod abscissas and weights should be computed by KRONSET.
!
!    The orders of the Gauss-Legendre and Kronrod rules must satisfy
!    NORDERK = 2 * NORDERG + 1.
!
!    The Kronrod rule uses the abscissas of the Gauss-Legendre rule,
!    plus more points, resulting in an efficient and higher order estimate.
!
!    The difference between the Gauss-Legendre and Kronrod estimates
!    is taken as an estimate of the error in the approximation to the
!    integral.
!
!  Modified:
!
!    15 September 1998
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, external FUNC, the name of the FORTRAN function which
!    evaluates the integrand.  The function must have the form
!      real(kind=Rkind) func ( x ).
!
!    Input, real(kind=Rkind) A, B, the lower and upper limits of integration.
!
!    Input, integer NSUB, the number of equal subintervals into
!    which the finite interval (A,B) is to be subdivided for
!    higher accuracy.  NSUB must be at least 1.
!
!    Input, integer NORDERG, the order of the Gauss-Legendre rule.
!    NORDERG must be at least 1.
!
!    Input, real(kind=Rkind) WEIGHTG(NORDERG), the weights of the
!    Gauss-Legendre rule.
!
!    Output, real(kind=Rkind) RESULTG, the approximate value of the
!    integral based on the Gauss-Legendre rule.
!
!    Input, integer NORDERK, the order of the Kronrod rule.
!    NORDERK must be at least 1.
!
!    Input, real(kind=Rkind) XTABK(NORDERK), the abscissas of the
!    Kronrod rule.
!
!    Input, real(kind=Rkind) WEIGHTK(NORDERK), the weights of the
!    Kronrod rule.
!
!    Output, real(kind=Rkind) RESULTK, the approximate value of the
!    integral based on the Kronrod rule.
!
!    Output, real(kind=Rkind) ERROR, an estimate of the approximation
!    error.  This is computed by taking the sum of the absolute values of
!    the differences between the Gauss-Legendre and Kronrod rules
!    over each subinterval.  This is usually a good estimate of
!    the error in the value RESULTG.  The error in the Kronrod
!    estimate RESULTK is usually much smaller.
!
      integer norderg
      integer norderk
!
      real(kind=Rkind) a
      real(kind=Rkind) b
      real(kind=Rkind) error
      real(kind=Rkind) fk
      real(kind=Rkind) func
      real(kind=Rkind) h
      integer i
      integer j
      integer nsub
      real(kind=Rkind) partg
      real(kind=Rkind) partk
      real(kind=Rkind) resultg
      real(kind=Rkind) resultk
      real(kind=Rkind) x
      real(kind=Rkind) xmid
      real(kind=Rkind) xtabk(norderk)
      real(kind=Rkind) weightg(norderg)
      real(kind=Rkind) weightk(norderk)
!
      external func
!
      resultg = 0.0 D+00
      resultk = 0.0 D+00
      error = 0

      if ( a .eq. b ) then
        return
      end if

      if ( norderk .ne. 2 * norderg + 1 ) then
        write ( *, * ) ' '
        write ( *, * ) 'SUMSUBGK - Fatal error!'
        write ( *, * ) '  NORDERK must equal 2 * NORDERG + 1.'
        write ( *, * ) '  The input value was NORDERG = ', norderg
        write ( *, * ) '  The input value was NORDERK = ', norderk
        stop
      end if

      h = ( b - a ) / dble ( nsub )

      do j = 1, nsub

        xmid = a + 0.5D0 * h * dble ( 2 * j - 1 )

        partg = 0.0 D+00
        partk = 0.0 D+00

        do i = 1, norderk

          x = xmid + 0.5D0 * h * xtabk(i)
          fk = func ( x )
          partk = partk + 0.5D0 * h * weightk(i) * fk

          if ( mod ( i, 2 ) .eq. 0 ) then
            partg = partg + 0.5D0 * h * weightg(i/2) * fk
          end if

        end do

        resultg = resultg + partg
        resultk = resultk + partk
        error = error + abs ( partk - partg )

      end do

      return
      end subroutine sumsubgk

