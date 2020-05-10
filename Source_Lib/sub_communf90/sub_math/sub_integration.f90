!=============================================================
!
!  ++  integration de -Pn1(x)*pot(x)*Pn2(x) de 1 a -1
!      en utilisant une quadrature de Gauss-Legendre
!
!=============================================================
      FUNCTION integration_pot(n1,n2,x,w,nb_gauss,poly_g,pot_g,nb_niv)
      USE mod_system
      IMPLICIT NONE
      real(kind=Rkind) integration_pot

       integer n1,n2,nb_gauss,nb_niv
       real(kind=Rkind) val
!      variables pour la subroutine gauleg
       real(kind=Rkind) x(nb_gauss),w(nb_gauss)
       real(kind=Rkind) poly_g(nb_niv,nb_gauss)
       real(kind=Rkind) pot_g(nb_gauss)

       integer i

       val = ZERO
       DO i=1,nb_gauss
         val = val + w(i)*poly_g(n1,i)*poly_g(n2,i)*pot_g(i)
       END DO

!      write(out_unitp,*) val,n1,n2
       integration_pot =  val

       end function integration_pot
!=============================================================
!
! ++   integration de -Pn1(x)*T_g1(x)*P'n2(x) de 1 a -1
!      en utilisant une quadrature de Gauss-Legendre
!
!=============================================================
      FUNCTION integration_T_g1(n1,n2,x,w,nb_gauss,poly_g,T_g1,nb_niv)
      USE mod_system
      IMPLICIT NONE
      real(kind=Rkind) :: integration_T_g1

       integer n1,n2,nb_gauss,nb_niv
       real(kind=Rkind) val
!      variables pour la subroutine gauleg
       real(kind=Rkind) x(nb_gauss),w(nb_gauss)
       real(kind=Rkind) poly_g(nb_niv,nb_gauss)
       real(kind=Rkind) T_g1(nb_gauss)
       real(kind=Rkind) d1poly_legendre

       integer i

       val = ZERO
       DO i=1,nb_gauss
!        write(out_unitp,*) ' integration_T_g1 :',i,x(i),
!    *                 d1poly_legendre(x(i),n2),T_g1(i)
         val = val + w(i)*poly_g(n1,i)*d1poly_legendre(x(i),n2)*T_g1(i)
       END DO

!      write(out_unitp,*) val,n1,n2
       integration_T_g1 =  val

       end function integration_T_g1
!=============================================================
!
! ++   generation des poids (w) et des abscisses (x) pour la quadrature
!      d'une particule dans une boite  [0,pi]
!
!=============================================================
      SUBROUTINE gauss_box(x,w,n)
      USE mod_system
      IMPLICIT NONE

       integer n
       real(kind=Rkind) w(n),x(n)

       integer i,nn

       real(kind=Rkind) pas,x0,wi


       pas = pi/real(n,kind=Rkind)
       wi  = pi/real(n,kind=Rkind)

       x0 = -pas*HALF

       DO i=1,n

         x0 = x0 + pas
         w(i)   = wi
         x(i)   = x0

!        write(out_unitp,*) i,nn+i,nn+1-i,wi,x0

       END DO

       end subroutine gauss_box
!=============================================================
!
! ++   generation des poids (w) et des abscisses (x) pour la quadrature
!      d'une particule dans une boite  [0,pi]
!
!=============================================================
      SUBROUTINE gauss_box_nosym(x,w,n)
      USE mod_system
      IMPLICIT NONE

       integer n
       real(kind=Rkind) w(n),x(n)

       integer i,nn

       real(kind=Rkind) pas,x0,wi


       pas = pi/real(n,kind=Rkind)
       wi  = pi/real(n,kind=Rkind)

       x0 = ZERO

       DO i=1,n

         x0 = x0 + pas
         w(i)   = wi
         x(i)   = x0

!        write(out_unitp,*) i,nn+i,nn+1-i,wi,x0

       END DO

       end subroutine gauss_box_nosym
!=============================================================
!
! ++   generation des poids (w) et des abscisses (x) pour la quadrature
!      de gauss-fourier
!
!=============================================================
      SUBROUTINE gauss_fourier(x,w,n)
      USE mod_system
      IMPLICIT NONE

       integer n
       real(kind=Rkind) w(n),x(n)

       integer i,nn

       real(kind=Rkind) pas,x0,wi


       pas = TWO*pi/real(n,kind=Rkind)
       wi  = pas
       x0 = -pi - HALF*pas

       DO i=1,n

!        x0 = x0 + pas
         x0 = -pi + (real(i,kind=Rkind)-HALF )*pas

         w(i)   = wi

         x(i)   = x0

!        write(out_unitp,*) i,wi,x0

       END DO

       end subroutine gauss_fourier
!=============================================================
!
! ++   gauss-chebychev quadrature
!
!=============================================================
      SUBROUTINE gauss_cheby(x,w,n)
      USE mod_system
      IMPLICIT NONE

       integer :: n
       real(kind=Rkind) w(n),x(n)

       integer :: i
       real(kind=Rkind) :: step,x0


       step = pi/real(n,kind=Rkind)
       x0 = - HALF*step

       DO i=1,n

         x(i) = cos(x0 + real(i,kind=Rkind)*step)

         !w(i) = step/sqrt(ONE-x(i)*x(i))
         w(i) = step

       END DO

       end subroutine gauss_cheby
      SUBROUTINE gauss_chebyWeight(x,w,n)
      USE mod_system
      IMPLICIT NONE

       integer :: n
       real(kind=Rkind) w(n),x(n)

       integer :: i
       real(kind=Rkind) :: step,x0


       step = pi/real(n,kind=Rkind)
       x0 = - HALF*step

       DO i=1,n

         x(i) = cos(x0 + real(i,kind=Rkind)*step)

         w(i) = step*sqrt(ONE-x(i)*x(i))
         !w(i) = step

       END DO

       end subroutine gauss_chebyWeight
!=============================================================

! ++   generation des poids (w) et des abscisses (x) pour la quadrature
!      de gauss-legendre
!      Numerical Recipes pp125 (1er ed?)
!
!=============================================================
      SUBROUTINE gauleg128(x1,x2,x,w,n)
      USE, intrinsic :: ISO_FORTRAN_ENV, ONLY : real128
      USE mod_system
      IMPLICIT NONE

       integer            :: n
       real(kind=Rkind)   :: x1,x2,x(n),w(n)

       integer :: m,i,j,it
       integer, parameter :: max_it = 100
       real (kind=real128) :: xm,xl,z,z1,p1,p2,p3,pp
       real (kind=real128), parameter ::                                   &
         pi128 = 3.14159265358979323846264338327950288419716939937511_real128

       m = (n+1)/2

       xm = 0.5_real128*(x2+x1)
       xl = 0.5_real128*(x2-x1)

       DO i=1,m
         z = cos(pi128*(real(i,kind=real128)-0.25_real128)/(0.5_real128+real(n,kind=real128)))

         DO it=1,max_it
           p1 = 1.0_real128
           p2 = 0.0_real128

           DO j=1,n
             p3 = p2
             p2 = p1
             p1 = (real(2*j-1,kind=real128)*z*p2 - real(j-1,kind=real128)*p3)/real(j,kind=real128)
           END DO

           pp = real(n,kind=real128)*(z*p1-p2)/(z*z-1.0_real128)
           z1 = z
           z  = z1-p1/pp

           IF (abs(z-z1) <= spacing(z)) EXIT
         END DO

         x(i)     = xm-xl*z
         x(n+1-i) = xm+xl*z
         w(i)     = 2.0_real128*xl/((1.0_real128-z*z)*pp*pp)
         w(n+1-i) = w(i)

       END DO

       RETURN
       end subroutine gauleg128
      SUBROUTINE gauleg(x1,x2,x,w,n)
      USE mod_system
      IMPLICIT NONE

       integer n,m,i,j,it
       !real(kind=Rkind), parameter :: eps = ONETENTH**14
       real(kind=Rkind), parameter :: eps = TINY(ONE)
       integer, parameter :: max_it = 50

       real(kind=Rkind) x1,x2,x(n),w(n)

       real(kind=Rkind) xm,xl,z,z1,p1,p2,p3,pp

       m = (n+1)/2

       xm = HALF*(x2+x1)
       xl = HALF*(x2-x1)

       DO i=1,m
         z = cos(pi*(real(i,kind=Rkind)-QUARTER)/(HALF+real(n,kind=Rkind)))

         DO it=1,max_it
           p1 = ONE
           p2 = ZERO

           DO j=1,n
             p3 = p2
             p2 = p1
             p1 = (real(2*j-1,kind=Rkind)*z*p2 - real(j-1,kind=Rkind)*p3)/real(j,kind=Rkind)
           END DO

           pp = real(n,kind=Rkind)*(z*p1-p2)/(z*z-ONE)
           z1 = z
           z  = z1-p1/pp

           IF (abs(z-z1) <= eps) EXIT
           !IF (abs(z-z1) <= TWO*spacing(z)) EXIT
         END DO

         x(i)     = xm-xl*z
         x(n+1-i) = xm+xl*z
         w(i)     = TWO*xl/((ONE-z*z)*pp*pp)
         w(n+1-i) = w(i)

       END DO

!      DO i=1,n
!       write(out_unitp,*) x(i),w(i)
!      END DO
!
       RETURN
       end subroutine gauleg
!=============================================================
!
! ++   generation des poids (w) et des abscisses (x) pour la quadrature
!      de Gauss-Hermite
!
!=============================================================

      SUBROUTINE hercom ( norder, xtab, weight )
      USE mod_system
      IMPLICIT NONE
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
      integer i
      real(kind=Rkind) p1
      real(kind=Rkind) s
      real(kind=Rkind) temp
      real(kind=Rkind) weight(norder)
      real(kind=Rkind) x
      real(kind=Rkind) xtab(norder)

!----- function -----------------------------------
      real(kind=Rkind) :: gamma,gamma_perso
!----- function -----------------------------------
!


      cc = 1.7724538509_Rkind * gamma_perso ( norder )                      &
        / ( TWO**(norder-1) )

      s = ( TWO * real(norder,kind=Rkind) + ONE )**(ONE/SIX)

      do i = 1, ( norder + 1 ) / 2

        if ( i .eq. 1 ) then

          x = s**3 - 1.85575_Rkind / s

        else if ( i .eq. 2 ) then

          x = x - 1.14_Rkind * ( real(norder,kind=Rkind)**0.426_Rkind ) / x

        else if ( i .eq. 3 ) then

          x = 1.86_Rkind * x - 0.86_Rkind * xtab(1)

        else if ( i .eq. 4 ) then

          x = 1.91_Rkind * x - 0.91_Rkind * xtab(2)

        else

          x = TWO * x - xtab(i-2)

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
!
!  Suprime exp(-x*x) de la varaiable weigth()
!
      do i = 1, norder
        weight(i) = weight(i)*exp(xtab(i)*xtab(i))
      end do



      return
      end subroutine hercom
      SUBROUTINE herrec ( p2, dp2, p1, x, norder )
      USE mod_system
      IMPLICIT NONE
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
      p1  = ONE
      dp1 = ZERO

      p2 = x
      dp2 = ONE

      do i = 2, norder

        p0 = p1
        dp0 = dp1

        p1 = p2
        dp1 = dp2

        p2  = x * p1       - HALF * real(i-1,kind=Rkind) * p0
        dp2 = x * dp1 + p1 - HALF * real(i-1,kind=Rkind) * dp0

      end do

      end subroutine herrec
      SUBROUTINE herroot ( x, norder, dp2, p1 )
      USE mod_system
      IMPLICIT NONE
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
      parameter ( eps = 1.0_Rkind / TEN**12 )
      !parameter ( eps = TINY(1.0_Rkind) )
      integer :: max_it = 30

!
      real(kind=Rkind) d
      real(kind=Rkind) dp2
      integer i
      integer norder
      real(kind=Rkind) p1
      real(kind=Rkind) p2
      real(kind=Rkind) x
!
      do i = 1, max_it

        call herrec ( p2, dp2, p1, x, norder )

        d = p2 / dp2
        x = x - d

        if ( abs ( d ) .le. eps * ( abs ( x ) + ONE ) ) then
          return
        end if

      end do

      return
      end subroutine herroot
      SUBROUTINE herset ( norder, xtab, weight )
      USE mod_system
      IMPLICIT NONE
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
!
      integer norder
!
      real(kind=Rkind) xtab(norder)
      real(kind=Rkind) weight(norder)
!
      if ( norder .eq. 1 ) then

        xtab(1) = ZERO

        weight(1) = sqrt ( pi )

      else if ( norder .eq. 2 )then

        xtab(1) = - 0.707106781186547524400844362105_Rkind
        xtab(2) =   0.707106781186547524400844362105_Rkind

        weight(1) = 0.886226925452758013649083741671_Rkind
        weight(2) = 0.886226925452758013649083741671_Rkind

      else if ( norder .eq. 3 ) then

        xtab(1) = - 0.122474487139158904909864203735_Rkind * TEN
        xtab(2) =   ZERO
        xtab(3) =   0.122474487139158904909864203735_Rkind * TEN

        weight(1) = 0.295408975150919337883027913890_Rkind
        weight(2) = 0.118163590060367735153211165556_Rkind * TEN
        weight(3) = 0.295408975150919337883027913890_Rkind

      else if ( norder .eq. 4 ) then

        xtab(1) = - 0.165068012388578455588334111112_Rkind * TEN
        xtab(2) = - 0.524647623275290317884060253835_Rkind
        xtab(3) =   0.524647623275290317884060253835_Rkind
        xtab(4) =   0.165068012388578455588334111112_Rkind * TEN

        weight(1) = 0.813128354472451771430345571899_Rkind / TEN
        weight(2) = 0.804914090005512836506049184481_Rkind
        weight(3) = 0.804914090005512836506049184481_Rkind
        weight(4) = 0.813128354472451771430345571899_Rkind / TEN

      else if ( norder .eq. 5 ) then

        xtab(1) = - 0.202018287045608563292872408814_Rkind * TEN
        xtab(2) = - 0.958572464613818507112770593893_Rkind
        xtab(3) =   ZERO
        xtab(4) =   0.958572464613818507112770593893_Rkind
        xtab(5) =   0.202018287045608563292872408814_Rkind * TEN

        weight(1) = 0.199532420590459132077434585942_Rkind / TEN
        weight(2) = 0.393619323152241159828495620852_Rkind
        weight(3) = 0.945308720482941881225689324449_Rkind
        weight(4) = 0.393619323152241159828495620852_Rkind
        weight(5) = 0.199532420590459132077434585942_Rkind / TEN

      else if ( norder .eq. 6 ) then

        xtab(1) = - 0.235060497367449222283392198706_Rkind * TEN
        xtab(2) = - 0.133584907401369694971489528297_Rkind * TEN
        xtab(3) = - 0.436077411927616508679215948251_Rkind
        xtab(4) =   0.436077411927616508679215948251_Rkind
        xtab(5) =   0.133584907401369694971489528297_Rkind * TEN
        xtab(6) =   0.235060497367449222283392198706_Rkind * TEN

        weight(1) = 0.453000990550884564085747256463_Rkind / TEN**2
        weight(2) = 0.157067320322856643916311563508_Rkind
        weight(3) = 0.724629595224392524091914705598_Rkind
        weight(4) = 0.724629595224392524091914705598_Rkind
        weight(5) = 0.157067320322856643916311563508_Rkind
        weight(6) = 0.453000990550884564085747256463_Rkind / TEN**2

      else if ( norder .eq. 7 ) then

        xtab(1) = - 0.265196135683523349244708200652_Rkind * TEN
        xtab(2) = - 0.167355162876747144503180139830_Rkind * TEN
        xtab(3) = - 0.816287882858964663038710959027_Rkind
        xtab(4) =   ZERO
        xtab(5) =   0.816287882858964663038710959027_Rkind
        xtab(6) =   0.167355162876747144503180139830_Rkind * TEN
        xtab(7) =   0.265196135683523349244708200652_Rkind * TEN

        weight(1) = 0.971781245099519154149424255939_Rkind / TEN**3
        weight(2) = 0.545155828191270305921785688417_Rkind / TEN
        weight(3) = 0.425607252610127800520317466666_Rkind
        weight(4) = 0.810264617556807326764876563813_Rkind
        weight(5) = 0.425607252610127800520317466666_Rkind
        weight(6) = 0.545155828191270305921785688417_Rkind / TEN
        weight(7) = 0.971781245099519154149424255939_Rkind / TEN**3

      else if ( norder .eq. 8 ) then

        xtab(1) = - 0.293063742025724401922350270524_Rkind * TEN
        xtab(2) = - 0.198165675669584292585463063977_Rkind * TEN
        xtab(3) = - 0.115719371244678019472076577906_Rkind * TEN
        xtab(4) = - 0.381186990207322116854718885584_Rkind
        xtab(5) =   0.381186990207322116854718885584_Rkind
        xtab(6) =   0.115719371244678019472076577906_Rkind * TEN
        xtab(7) =   0.198165675669584292585463063977_Rkind * TEN
        xtab(8) =   0.293063742025724401922350270524_Rkind * TEN

        weight(1) = 0.199604072211367619206090452544_Rkind / TEN**3
        weight(2) = 0.170779830074134754562030564364_Rkind / TEN
        weight(3) = 0.207802325814891879543258620286_Rkind
        weight(4) = 0.661147012558241291030415974496_Rkind
        weight(5) = 0.661147012558241291030415974496_Rkind
        weight(6) = 0.207802325814891879543258620286_Rkind
        weight(7) = 0.170779830074134754562030564364_Rkind / TEN
        weight(8) = 0.199604072211367619206090452544_Rkind / TEN**3

      else if ( norder .eq. 9 ) then

        xtab(1) = - 0.319099320178152760723004779538_Rkind * TEN
        xtab(2) = - 0.226658058453184311180209693284_Rkind * TEN
        xtab(3) = - 0.146855328921666793166701573925_Rkind * TEN
        xtab(4) = - 0.723551018752837573322639864579_Rkind
        xtab(5) =   ZERO
        xtab(6) =   0.723551018752837573322639864579_Rkind
        xtab(7) =   0.146855328921666793166701573925_Rkind * TEN
        xtab(8) =   0.226658058453184311180209693284_Rkind * TEN
        xtab(9) =   0.319099320178152760723004779538_Rkind * TEN

        weight(1) = 0.396069772632643819045862946425_Rkind / TEN**4
        weight(2) = 0.494362427553694721722456597763_Rkind / TEN**2
        weight(3) = 0.884745273943765732879751147476_Rkind / TEN
        weight(4) = 0.432651559002555750199812112956_Rkind
        weight(5) = 0.720235215606050957124334723389_Rkind
        weight(6) = 0.432651559002555750199812112956_Rkind
        weight(7) = 0.884745273943765732879751147476_Rkind / TEN
        weight(8) = 0.494362427553694721722456597763_Rkind / TEN**2
        weight(9) = 0.396069772632643819045862946425_Rkind / TEN**4

      else if ( norder .eq. 10 ) then

        xtab(1) =  - 0.343615911883773760332672549432_Rkind * TEN
        xtab(2) =  - 0.253273167423278979640896079775_Rkind * TEN
        xtab(3) =  - 0.175668364929988177345140122011_Rkind * TEN
        xtab(4) =  - 0.103661082978951365417749191676_Rkind * TEN
        xtab(5) =  - 0.342901327223704608789165025557_Rkind
        xtab(6) =    0.342901327223704608789165025557_Rkind
        xtab(7) =    0.103661082978951365417749191676_Rkind * TEN
        xtab(8) =    0.175668364929988177345140122011_Rkind * TEN
        xtab(9) =    0.253273167423278979640896079775_Rkind * TEN
        xtab(10) =   0.343615911883773760332672549432_Rkind * TEN

        weight(1) =  0.764043285523262062915936785960_Rkind / TEN**5
        weight(2) =  0.134364574678123269220156558585_Rkind / TEN**2
        weight(3) =  0.338743944554810631361647312776_Rkind / TEN
        weight(4) =  0.240138611082314686416523295006_Rkind
        weight(5) =  0.610862633735325798783564990433_Rkind
        weight(6) =  0.610862633735325798783564990433_Rkind
        weight(7) =  0.240138611082314686416523295006_Rkind
        weight(8) =  0.338743944554810631361647312776_Rkind / TEN
        weight(9) =  0.134364574678123269220156558585_Rkind / TEN**2
        weight(10) = 0.764043285523262062915936785960_Rkind / TEN**5

      else if ( norder .eq. 11 ) then

        xtab(1) =  - 0.366847084655958251845837146485_Rkind * TEN
        xtab(2) =  - 0.278329009978165177083671870152_Rkind * TEN
        xtab(3) =  - 0.202594801582575533516591283121_Rkind * TEN
        xtab(4) =  - 0.132655708449493285594973473558_Rkind * TEN
        xtab(5) =  - 0.656809566882099765024611575383_Rkind
        xtab(6) =    ZERO
        xtab(7) =    0.656809566882099765024611575383_Rkind
        xtab(8) =    0.132655708449493285594973473558_Rkind * TEN
        xtab(9) =    0.202594801582575533516591283121_Rkind * TEN
        xtab(10) =   0.278329009978165177083671870152_Rkind * TEN
        xtab(11) =   0.366847084655958251845837146485_Rkind * TEN

        weight(1) =  0.143956039371425822033088366032_Rkind / TEN**5
        weight(2) =  0.346819466323345510643413772940_Rkind / TEN**3
        weight(3) =  0.119113954449115324503874202916_Rkind / TEN
        weight(4) =  0.117227875167708503381788649308_Rkind
        weight(5) =  0.429359752356125028446073598601_Rkind
        weight(6) =  0.654759286914591779203940657627_Rkind
        weight(7) =  0.429359752356125028446073598601_Rkind
        weight(8) =  0.117227875167708503381788649308_Rkind
        weight(9) =  0.119113954449115324503874202916_Rkind / TEN
        weight(10) = 0.346819466323345510643413772940_Rkind / TEN**3
        weight(11) = 0.143956039371425822033088366032_Rkind / TEN**5

      else if ( norder .eq. 12 ) then

        xtab(1) =  - 0.388972489786978191927164274724_Rkind * TEN
        xtab(2) =  - 0.302063702512088977171067937518_Rkind * TEN
        xtab(3) =  - 0.227950708050105990018772856942_Rkind * TEN
        xtab(4) =  - 0.159768263515260479670966277090_Rkind * TEN
        xtab(5) =  - 0.947788391240163743704578131060_Rkind
        xtab(6) =  - 0.314240376254359111276611634095_Rkind
        xtab(7) =    0.314240376254359111276611634095_Rkind
        xtab(8) =    0.947788391240163743704578131060_Rkind
        xtab(9) =    0.159768263515260479670966277090_Rkind * TEN
        xtab(10) =   0.227950708050105990018772856942_Rkind * TEN
        xtab(11) =   0.302063702512088977171067937518_Rkind * TEN
        xtab(12) =   0.388972489786978191927164274724_Rkind * TEN

        weight(1) =  0.265855168435630160602311400877_Rkind / TEN**6
        weight(2) =  0.857368704358785865456906323153_Rkind / TEN**4
        weight(3) =  0.390539058462906185999438432620_Rkind / TEN**2
        weight(4) =  0.516079856158839299918734423606_Rkind / TEN
        weight(5) =  0.260492310264161129233396139765_Rkind
        weight(6) =  0.570135236262479578347113482275_Rkind
        weight(7) =  0.570135236262479578347113482275_Rkind
        weight(8) =  0.260492310264161129233396139765_Rkind
        weight(9) =  0.516079856158839299918734423606_Rkind / TEN
        weight(10) = 0.390539058462906185999438432620_Rkind / TEN**2
        weight(11) = 0.857368704358785865456906323153_Rkind / TEN**4
        weight(12) = 0.265855168435630160602311400877_Rkind / TEN**6

      else if ( norder .eq. 13 ) then

        xtab(1) =  - 0.410133759617863964117891508007_Rkind * TEN
        xtab(2) =  - 0.324660897837240998812205115236_Rkind * TEN
        xtab(3) =  - 0.251973568567823788343040913628_Rkind * TEN
        xtab(4) =  - 0.185310765160151214200350644316_Rkind * TEN
        xtab(5) =  - 0.122005503659074842622205526637_Rkind * TEN
        xtab(6) =  - 0.605763879171060113080537108602_Rkind
        xtab(7) =    ZERO
        xtab(8) =    0.605763879171060113080537108602_Rkind
        xtab(9) =    0.122005503659074842622205526637_Rkind * TEN
        xtab(10) =   0.185310765160151214200350644316_Rkind * TEN
        xtab(11) =   0.251973568567823788343040913628_Rkind * TEN
        xtab(12) =   0.324660897837240998812205115236_Rkind * TEN
        xtab(13) =   0.410133759617863964117891508007_Rkind * TEN

        weight(1) =  0.482573185007313108834997332342_Rkind / TEN**7
        weight(2) =  0.204303604027070731248669432937_Rkind / TEN**4
        weight(3) =  0.120745999271938594730924899224_Rkind / TEN**2
        weight(4) =  0.208627752961699392166033805050_Rkind / TEN
        weight(5) =  0.140323320687023437762792268873_Rkind
        weight(6) =  0.421616296898543221746893558568_Rkind
        weight(7) =  0.604393187921161642342099068579_Rkind
        weight(8) =  0.421616296898543221746893558568_Rkind
        weight(9) =  0.140323320687023437762792268873_Rkind
        weight(10) = 0.208627752961699392166033805050_Rkind / TEN
        weight(11) = 0.120745999271938594730924899224_Rkind / TEN**2
        weight(12) = 0.204303604027070731248669432937_Rkind / TEN**4
        weight(13) = 0.482573185007313108834997332342_Rkind / TEN**7

      else if ( norder .eq. 14 ) then

        xtab(1) =  - 0.430444857047363181262129810037_Rkind * TEN
        xtab(2) =  - 0.346265693360227055020891736115_Rkind * TEN
        xtab(3) =  - 0.274847072498540256862499852415_Rkind * TEN
        xtab(4) =  - 0.209518325850771681573497272630_Rkind * TEN
        xtab(5) =  - 0.147668273114114087058350654421_Rkind * TEN
        xtab(6) =  - 0.878713787329399416114679311861_Rkind
        xtab(7) =  - 0.291745510672562078446113075799_Rkind
        xtab(8) =    0.291745510672562078446113075799_Rkind
        xtab(9) =    0.878713787329399416114679311861_Rkind
        xtab(10) =   0.147668273114114087058350654421_Rkind * TEN
        xtab(11) =   0.209518325850771681573497272630_Rkind * TEN
        xtab(12) =   0.274847072498540256862499852415_Rkind * TEN
        xtab(13) =   0.346265693360227055020891736115_Rkind * TEN
        xtab(14) =   0.430444857047363181262129810037_Rkind * TEN

        weight(1) =  0.862859116812515794532041783429_Rkind / TEN**8
        weight(2) =  0.471648435501891674887688950105_Rkind / TEN**5
        weight(3) =  0.355092613551923610483661076691_Rkind / TEN**3
        weight(4) =  0.785005472645794431048644334608_Rkind / TEN**2
        weight(5) =  0.685055342234652055387163312367_Rkind / TEN
        weight(6) =  0.273105609064246603352569187026_Rkind
        weight(7) =  0.536405909712090149794921296776_Rkind
        weight(8) =  0.536405909712090149794921296776_Rkind
        weight(9) =  0.273105609064246603352569187026_Rkind
        weight(10) = 0.685055342234652055387163312367_Rkind / TEN
        weight(11) = 0.785005472645794431048644334608_Rkind / TEN**2
        weight(12) = 0.355092613551923610483661076691_Rkind / TEN**3
        weight(13) = 0.471648435501891674887688950105_Rkind / TEN**5
        weight(14) = 0.862859116812515794532041783429_Rkind / TEN**8

      else if ( norder .eq. 15 ) then

        xtab(1) =  - 0.449999070730939155366438053053_Rkind * TEN
        xtab(2) =  - 0.366995037340445253472922383312_Rkind * TEN
        xtab(3) =  - 0.296716692790560324848896036355_Rkind * TEN
        xtab(4) =  - 0.232573248617385774545404479449_Rkind * TEN
        xtab(5) =  - 0.171999257518648893241583152515_Rkind * TEN
        xtab(6) =  - 0.113611558521092066631913490556_Rkind * TEN
        xtab(7) =  - 0.565069583255575748526020337198_Rkind
        xtab(8) =    ZERO
        xtab(9) =    0.565069583255575748526020337198_Rkind
        xtab(10) =   0.113611558521092066631913490556_Rkind * TEN
        xtab(11) =   0.171999257518648893241583152515_Rkind * TEN
        xtab(12) =   0.232573248617385774545404479449_Rkind * TEN
        xtab(13) =   0.296716692790560324848896036355_Rkind * TEN
        xtab(14) =   0.366995037340445253472922383312_Rkind * TEN
        xtab(15) =   0.449999070730939155366438053053_Rkind * TEN

        weight(1) =  0.152247580425351702016062666965_Rkind / TEN**8
        weight(2) =  0.105911554771106663577520791055_Rkind / TEN**5
        weight(3) =  0.100004441232499868127296736177_Rkind / TEN**3
        weight(4) =  0.277806884291277589607887049229_Rkind / TEN**2
        weight(5) =  0.307800338725460822286814158758_Rkind / TEN
        weight(6) =  0.158488915795935746883839384960_Rkind
        weight(7) =  0.412028687498898627025891079568_Rkind
        weight(8) =  0.564100308726417532852625797340_Rkind
        weight(9) =  0.412028687498898627025891079568_Rkind
        weight(10) = 0.158488915795935746883839384960_Rkind
        weight(11) = 0.307800338725460822286814158758_Rkind / TEN
        weight(12) = 0.277806884291277589607887049229_Rkind / TEN**2
        weight(13) = 0.100004441232499868127296736177_Rkind / TEN**3
        weight(14) = 0.105911554771106663577520791055_Rkind / TEN**5
        weight(15) = 0.152247580425351702016062666965_Rkind / TEN**8

      else if ( norder .eq. 16 ) then

        xtab(1) =  - 0.468873893930581836468849864875_Rkind * TEN
        xtab(2) =  - 0.386944790486012269871942409801_Rkind * TEN
        xtab(3) =  - 0.317699916197995602681399455926_Rkind * TEN
        xtab(4) =  - 0.254620215784748136215932870545_Rkind * TEN
        xtab(5) =  - 0.195178799091625397743465541496_Rkind * TEN
        xtab(6) =  - 0.138025853919888079637208966969_Rkind * TEN
        xtab(7) =  - 0.822951449144655892582454496734_Rkind
        xtab(8) =  - 0.273481046138152452158280401965_Rkind
        xtab(9) =    0.273481046138152452158280401965_Rkind
        xtab(10) =   0.822951449144655892582454496734_Rkind
        xtab(11) =   0.138025853919888079637208966969_Rkind * TEN
        xtab(12) =   0.195178799091625397743465541496_Rkind * TEN
        xtab(13) =   0.254620215784748136215932870545_Rkind * TEN
        xtab(14) =   0.317699916197995602681399455926_Rkind * TEN
        xtab(15) =   0.386944790486012269871942409801_Rkind * TEN
        xtab(16) =   0.468873893930581836468849864875_Rkind * TEN

        weight(1) =  0.265480747401118224470926366050_Rkind / TEN**9
        weight(2) =  0.232098084486521065338749423185_Rkind / TEN**6
        weight(3) =  0.271186009253788151201891432244_Rkind / TEN**4
        weight(4) =  0.932284008624180529914277305537_Rkind / TEN**3
        weight(5) =  0.128803115355099736834642999312_Rkind / TEN
        weight(6) =  0.838100413989858294154207349001_Rkind / TEN
        weight(7) =  0.280647458528533675369463335380_Rkind
        weight(8) =  0.507929479016613741913517341791_Rkind
        weight(9) =  0.507929479016613741913517341791_Rkind
        weight(10) = 0.280647458528533675369463335380_Rkind
        weight(11) = 0.838100413989858294154207349001_Rkind / TEN
        weight(12) = 0.128803115355099736834642999312_Rkind / TEN
        weight(13) = 0.932284008624180529914277305537_Rkind / TEN**3
        weight(14) = 0.271186009253788151201891432244_Rkind / TEN**4
        weight(15) = 0.232098084486521065338749423185_Rkind / TEN**6
        weight(16) = 0.265480747401118224470926366050_Rkind / TEN**9

      else if ( norder .eq. 17 ) then

        xtab(1) =  - 0.487134519367440308834927655662_Rkind * TEN
        xtab(2) =  - 0.406194667587547430689245559698_Rkind * TEN
        xtab(3) =  - 0.337893209114149408338327069289_Rkind * TEN
        xtab(4) =  - 0.275776291570388873092640349574_Rkind * TEN
        xtab(5) =  - 0.217350282666662081927537907149_Rkind * TEN
        xtab(6) =  - 0.161292431422123133311288254454_Rkind * TEN
        xtab(7) =  - 0.106764872574345055363045773799_Rkind * TEN
        xtab(8) =  - 0.531633001342654731349086553718_Rkind
        xtab(9) =    ZERO
        xtab(10) =   0.531633001342654731349086553718_Rkind
        xtab(11) =   0.106764872574345055363045773799_Rkind * TEN
        xtab(12) =   0.161292431422123133311288254454_Rkind * TEN
        xtab(13) =   0.217350282666662081927537907149_Rkind * TEN
        xtab(14) =   0.275776291570388873092640349574_Rkind * TEN
        xtab(15) =   0.337893209114149408338327069289_Rkind * TEN
        xtab(16) =   0.406194667587547430689245559698_Rkind * TEN
        xtab(17) =   0.487134519367440308834927655662_Rkind * TEN

        weight(1) =  0.458057893079863330580889281222_Rkind / TEN**10
        weight(2) =  0.497707898163079405227863353715_Rkind / TEN**7
        weight(3) =  0.711228914002130958353327376218_Rkind / TEN**5
        weight(4) =  0.298643286697753041151336643059_Rkind / TEN**3
        weight(5) =  0.506734995762753791170069495879_Rkind / TEN**2
        weight(6) =  0.409200341495762798094994877854_Rkind / TEN
        weight(7) =  0.172648297670097079217645196219_Rkind
        weight(8) =  0.401826469470411956577635085257_Rkind
        weight(9) =  0.530917937624863560331883103379_Rkind
        weight(10) = 0.401826469470411956577635085257_Rkind
        weight(11) = 0.172648297670097079217645196219_Rkind
        weight(12) = 0.409200341495762798094994877854_Rkind / TEN
        weight(13) = 0.506734995762753791170069495879_Rkind / TEN**2
        weight(14) = 0.298643286697753041151336643059_Rkind / TEN**3
        weight(15) = 0.711228914002130958353327376218_Rkind / TEN**5
        weight(16) = 0.497707898163079405227863353715_Rkind / TEN**7
        weight(17) = 0.458057893079863330580889281222_Rkind / TEN**10

      else if ( norder .eq. 18 ) then

        xtab(1) =  - 0.504836400887446676837203757885_Rkind * TEN
        xtab(2) =  - 0.424811787356812646302342016090_Rkind * TEN
        xtab(3) =  - 0.357376906848626607950067599377_Rkind * TEN
        xtab(4) =  - 0.296137750553160684477863254906_Rkind * TEN
        xtab(5) =  - 0.238629908916668600026459301424_Rkind * TEN
        xtab(6) =  - 0.183553160426162889225383944409_Rkind * TEN
        xtab(7) =  - 0.130092085838961736566626555439_Rkind * TEN
        xtab(8) =  - 0.776682919267411661316659462284_Rkind
        xtab(9) =  - 0.258267750519096759258116098711_Rkind
        xtab(10) =   0.258267750519096759258116098711_Rkind
        xtab(11) =   0.776682919267411661316659462284_Rkind
        xtab(12) =   0.130092085838961736566626555439_Rkind * TEN
        xtab(13) =   0.183553160426162889225383944409_Rkind * TEN
        xtab(14) =   0.238629908916668600026459301424_Rkind * TEN
        xtab(15) =   0.296137750553160684477863254906_Rkind * TEN
        xtab(16) =   0.357376906848626607950067599377_Rkind * TEN
        xtab(17) =   0.424811787356812646302342016090_Rkind * TEN
        xtab(18) =   0.504836400887446676837203757885_Rkind * TEN

        weight(1) =  0.782819977211589102925147471012_Rkind / TEN**11
        weight(2) =  0.104672057957920824443559608435_Rkind / TEN**7
        weight(3) =  0.181065448109343040959702385911_Rkind / TEN**5
        weight(4) =  0.918112686792940352914675407371_Rkind / TEN**4
        weight(5) =  0.188852263026841789438175325426_Rkind / TEN**2
        weight(6) =  0.186400423875446519219315221973_Rkind / TEN
        weight(7) =  0.973017476413154293308537234155_Rkind / TEN
        weight(8) =  0.284807285669979578595606820713_Rkind
        weight(9) =  0.483495694725455552876410522141_Rkind
        weight(10) = 0.483495694725455552876410522141_Rkind
        weight(11) = 0.284807285669979578595606820713_Rkind
        weight(12) = 0.973017476413154293308537234155_Rkind / TEN
        weight(13) = 0.186400423875446519219315221973_Rkind / TEN
        weight(14) = 0.188852263026841789438175325426_Rkind / TEN**2
        weight(15) = 0.918112686792940352914675407371_Rkind / TEN**4
        weight(16) = 0.181065448109343040959702385911_Rkind / TEN**5
        weight(17) = 0.104672057957920824443559608435_Rkind / TEN**7
        weight(18) = 0.782819977211589102925147471012_Rkind / TEN**11

      else if ( norder .eq. 19 ) then

        xtab(1) =  - 0.522027169053748216460967142500_Rkind * TEN
        xtab(2) =  - 0.442853280660377943723498532226_Rkind * TEN
        xtab(3) =  - 0.376218735196402009751489394104_Rkind * TEN
        xtab(4) =  - 0.315784881834760228184318034120_Rkind * TEN
        xtab(5) =  - 0.259113378979454256492128084112_Rkind * TEN
        xtab(6) =  - 0.204923170985061937575050838669_Rkind * TEN
        xtab(7) =  - 0.152417061939353303183354859367_Rkind * TEN
        xtab(8) =  - 0.101036838713431135136859873726_Rkind * TEN
        xtab(9) =  - 0.503520163423888209373811765050_Rkind
        xtab(10) =   ZERO
        xtab(11) =   0.503520163423888209373811765050_Rkind
        xtab(12) =   0.101036838713431135136859873726_Rkind * TEN
        xtab(13) =   0.152417061939353303183354859367_Rkind * TEN
        xtab(14) =   0.204923170985061937575050838669_Rkind * TEN
        xtab(15) =   0.259113378979454256492128084112_Rkind * TEN
        xtab(16) =   0.315784881834760228184318034120_Rkind * TEN
        xtab(17) =   0.376218735196402009751489394104_Rkind * TEN
        xtab(18) =   0.442853280660377943723498532226_Rkind * TEN
        xtab(19) =   0.522027169053748216460967142500_Rkind * TEN

        weight(1) =  0.132629709449851575185289154385_Rkind / TEN**11
        weight(2) =  0.216305100986355475019693077221_Rkind / TEN**8
        weight(3) =  0.448824314722312295179447915594_Rkind / TEN**6
        weight(4) =  0.272091977631616257711941025214_Rkind / TEN**4
        weight(5) =  0.670877521407181106194696282100_Rkind / TEN**3
        weight(6) =  0.798886677772299020922211491861_Rkind / TEN**2
        weight(7) =  0.508103869090520673569908110358_Rkind / TEN
        weight(8) =  0.183632701306997074156148485766_Rkind
        weight(9) =  0.391608988613030244504042313621_Rkind
        weight(10) = 0.502974888276186530840731361096_Rkind
        weight(11) = 0.391608988613030244504042313621_Rkind
        weight(12) = 0.183632701306997074156148485766_Rkind
        weight(13) = 0.508103869090520673569908110358_Rkind / TEN
        weight(14) = 0.798886677772299020922211491861_Rkind / TEN**2
        weight(15) = 0.670877521407181106194696282100_Rkind / TEN**3
        weight(16) = 0.272091977631616257711941025214_Rkind / TEN**4
        weight(17) = 0.448824314722312295179447915594_Rkind / TEN**6
        weight(18) = 0.216305100986355475019693077221_Rkind / TEN**8
        weight(19) = 0.132629709449851575185289154385_Rkind / TEN**11

      else if ( norder .eq. 20 ) then

        xtab(1) =  - 0.538748089001123286201690041068_Rkind * TEN
        xtab(2) =  - 0.460368244955074427307767524898_Rkind * TEN
        xtab(3) =  - 0.394476404011562521037562880052_Rkind * TEN
        xtab(4) =  - 0.334785456738321632691492452300_Rkind * TEN
        xtab(5) =  - 0.278880605842813048052503375640_Rkind * TEN
        xtab(6) =  - 0.225497400208927552308233334473_Rkind * TEN
        xtab(7) =  - 0.173853771211658620678086566214_Rkind * TEN
        xtab(8) =  - 0.123407621539532300788581834696_Rkind * TEN
        xtab(9) =  - 0.737473728545394358705605144252_Rkind
        xtab(10) = - 0.245340708300901249903836530634_Rkind
        xtab(11) =   0.245340708300901249903836530634_Rkind
        xtab(12) =   0.737473728545394358705605144252_Rkind
        xtab(13) =   0.123407621539532300788581834696_Rkind * TEN
        xtab(14) =   0.173853771211658620678086566214_Rkind * TEN
        xtab(15) =   0.225497400208927552308233334473_Rkind * TEN
        xtab(16) =   0.278880605842813048052503375640_Rkind * TEN
        xtab(17) =   0.334785456738321632691492452300_Rkind * TEN
        xtab(18) =   0.394476404011562521037562880052_Rkind * TEN
        xtab(19) =   0.460368244955074427307767524898_Rkind * TEN
        xtab(20) =   0.538748089001123286201690041068_Rkind * TEN

        weight(1) =  0.222939364553415129252250061603_Rkind / TEN**12
        weight(2) =  0.439934099227318055362885145547_Rkind / TEN**9
        weight(3) =  0.108606937076928169399952456345_Rkind / TEN**6
        weight(4) =  0.780255647853206369414599199965_Rkind / TEN**5
        weight(5) =  0.228338636016353967257145917963_Rkind / TEN**3
        weight(6) =  0.324377334223786183218324713235_Rkind / TEN**2
        weight(7) =  0.248105208874636108821649525589_Rkind / TEN
        weight(8) =  0.109017206020023320013755033535_Rkind
        weight(9) =  0.286675505362834129719659706228_Rkind
        weight(10) = 0.462243669600610089650328639861_Rkind
        weight(11) = 0.462243669600610089650328639861_Rkind
        weight(12) = 0.286675505362834129719659706228_Rkind
        weight(13) = 0.109017206020023320013755033535_Rkind
        weight(14) = 0.248105208874636108821649525589_Rkind / TEN
        weight(15) = 0.324377334223786183218324713235_Rkind / TEN**2
        weight(16) = 0.228338636016353967257145917963_Rkind / TEN**3
        weight(17) = 0.780255647853206369414599199965_Rkind / TEN**5
        weight(18) = 0.108606937076928169399952456345_Rkind / TEN**6
        weight(19) = 0.439934099227318055362885145547_Rkind / TEN**9
        weight(20) = 0.222939364553415129252250061603_Rkind / TEN**12

      else

        write ( *, * ) ' '
        write ( *, * ) 'HERSET - Fatal error!'
        write ( *, * ) '  Illegal value of NORDER = ', norder
        write ( *, * ) '  Legal values are 1 to 20.'
        stop

      end if

      return
      end subroutine herset

! cubature from Burkardt
subroutine en_r2_01_1 ( n, o, x, w )

!*****************************************************************************80
!
!! EN_R2_01_1 implements the Stroud rule 1.1 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 1.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind=Rkind ) X(N,O), the abscissas.
!
!    Output, real ( kind=Rkind ) W(O), the weights.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  integer ( kind = 4 ) k
  real ( kind=Rkind ) volume
  real ( kind=Rkind ) w(o)
  real ( kind=Rkind ) x(n,o)

  volume = sqrt ( pi**n )

  x(1:n,1:o) = ZERO

  k = 0
!
!  1 point.
!
  k = k + 1
  w(k) = volume

  return
end
subroutine en_r2_01_1_size ( n, o )

!*****************************************************************************80
!
!! EN_R2_01_1_SIZE sizes the Stroud rule 1.1 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 1.
!
!    The rule has precision P = 1.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 1

  return
end
subroutine en_r2_02_xiu ( n, o, x, w )

!*****************************************************************************80
!
!! EN_R2_02_XIU implements the Xiu precision 2 rule for region EN_R2.
!
!  Discussion:
!
!    The rule has order
!
!      O = N + 1.
!
!    The rule has precision P = 2.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 March 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind=Rkind ) X(N,O), the abscissas.
!
!    Output, real ( kind=Rkind ) W(O), the weights.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind=Rkind ) arg
  real ( kind=Rkind ) c1
  real ( kind=Rkind ) delta0
  real ( kind=Rkind ) gamma0
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) r
  real ( kind=Rkind ) r8_mop
  real ( kind=Rkind ) volume
  real ( kind=Rkind ) volume_1d
  real ( kind=Rkind ) w(o)
  real ( kind=Rkind ) x(n,o)

  do j = 1, o

    i = 0
    do r = 1, n / 2

      arg = real ( 2 * r * ( j - 1 ), kind=Rkind ) * pi / real ( n + 1, kind=Rkind )

      i = i + 1
      x(i,j) = sqrt ( TWO ) * cos ( arg )
      i = i + 1
      x(i,j) = sqrt ( TWO ) * sin ( arg )

    end do

    if ( i < n ) then
      i = i + 1
      x(i,j) = (-ONE) ** ( j - 1 )
    end if

  end do

  gamma0 = TWO
  delta0 = ZERO
  c1 = ONE

  x(1:n,1:o) = ( sqrt ( gamma0 * c1 ) * x(1:n,1:o) - delta0 ) / gamma0

  volume_1d = sqrt ( pi )
  volume = volume_1d ** n

  w(1:o) = volume / real ( o, kind=Rkind )

  return
end
subroutine en_r2_02_xiu_size ( n, o )

!*****************************************************************************80
!
!! EN_R2_02_XIU_SIZE sizes the Xiu rule for region EN_R2.
!
!  Discussion:
!
!    The rule has order
!
!      O = N + 1.
!
!    The rule has precision P = 2.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = n + 1

  return
end
subroutine en_r2_03_1 ( n, o, x, w )

!*****************************************************************************80
!
!! EN_R2_03_1 implements the Stroud rule 3.1 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 2 * N.
!
!    The rule has precision P = 3.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = exp ( - x1^2 - x2^2 ... - xn^2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind=Rkind ) X(N,O), the abscissas.
!
!    Output, real ( kind=Rkind ) W(O), the weights.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind=Rkind ) a
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  real ( kind=Rkind ) r
  real ( kind=Rkind ) volume
  real ( kind=Rkind ) w(o)
  real ( kind=Rkind ) x(n,o)

  volume = sqrt ( pi**n )

  a = volume / real ( o, kind=Rkind )
  r = sqrt ( real ( n, kind=Rkind ) / TWO )

  x(1:n,1:o) = ZERO

  k = 0
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - r
    w(k) = a
    k = k + 1
    x(i,k) = + r
    w(k) = a
  end do

  return
end
subroutine en_r2_03_1_size ( n, o )

!*****************************************************************************80
!
!! EN_R2_03_1_SIZE sizes the Stroud rule 3.1 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 2 * N.
!
!    The rule has precision P = 3.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 2 * n

  return
end
subroutine en_r2_03_2 ( n, o, x, w )

!*****************************************************************************80
!
!! EN_R2_03_2 implements the Stroud rule 3.2 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 2^N.
!
!    The rule has precision P = 3.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    06 June 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind=Rkind ) X(N,O), the abscissas.
!
!    Output, real ( kind=Rkind ) W(O), the weights.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind=Rkind ) a
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  logical more
  real ( kind=Rkind ) r
  real ( kind=Rkind ) volume
  real ( kind=Rkind ) w(o)
  real ( kind=Rkind ) x(n,o)

  volume = sqrt ( pi ** n )

  a = volume / real ( o, kind=Rkind )
  r = sqrt ( HALF );

  x(1:n,1:o) = ZERO

  k = 0
!
!  2^N points.
!
  k = k + 1
  x(1:n,k) = - r
  w(k) = a
  more = .true.

  do while ( more )
    more = .false.
    do i = n, 1, -1
      if ( x(i,k) < ZERO ) then
        k = k + 1
        x(1:i-1,k) = x(1:i-1,k-1)
        x(i,k)     = + r
        x(i+1:n,k) = - r
        w(k) = a
        more = .true.
        exit
      end if
    end do
  end do

  return
end
subroutine en_r2_03_2_size ( n, o )

!*****************************************************************************80
!
!! EN_R2_03_2_SIZE sizes the Stroud rule 3.2 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 2^N.
!
!    The rule has precision P = 3.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 2 ** n;

  return
end
subroutine en_r2_03_xiu ( n, o, x, w )

!*****************************************************************************80
!
!! EN_R2_03_XIU implements the Xiu precision 3 rule for region EN_R2.
!
!  Discussion:
!
!    The rule has order
!
!      O = 2 * N.
!
!    The rule has precision P = 3.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind=Rkind ) X(N,O), the abscissas.
!
!    Output, real ( kind=Rkind ) W(O), the weights.
!
  USE mod_system
  implicit none

  real ( kind=Rkind ) arg
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) n
  integer ( kind = 4 ) o
  integer ( kind = 4 ) r
  real ( kind=Rkind ) r8_mop
  real ( kind=Rkind ) volume
  real ( kind=Rkind ) :: w(o)
  real ( kind=Rkind ) :: x(n,o)

  volume = sqrt ( pi ** n )

  do j = 1, o

    i = 0
    do r = 1, n / 2
      arg = real ( ( 2 * r - 1 ) * j, kind=Rkind ) * pi / real ( n, kind=Rkind )
      i = i + 1
      x(i,j) = cos ( arg )
      i = i + 1
      x(i,j) = sin ( arg )
    end do

    if ( i < n ) then
      i = i + 1
      x(i,j) = (-ONE) ** ( j )
      if ( n == 1 ) then
        x(i,j) = x(i,j) / sqrt ( TWO )
      end if
    end if

  end do

  w(1:o) = volume / real ( o, kind=Rkind )

  return
end
subroutine en_r2_03_xiu_size ( n, o )

!*****************************************************************************80
!
!! EN_R2_03_XIU_SIZE sizes the Xiu precision 3 rule for region EN_R2.
!
!  Discussion:
!
!    The rule has order
!
!      O = 2 * N.
!
!    The rule has precision P = 3.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    29 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Dongbin Xiu,
!    Numerical integration formulas of degree two,
!    Applied Numerical Mathematics,
!    Volume 58, 2008, pages 1515-1520.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 2 * n

  return
end
subroutine en_r2_05_1 ( n, option, o, x, w )

!*****************************************************************************80
!
!! EN_R2_05_1 implements the Stroud rule 5.1 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = N^2 + N + 2.
!
!    The rule has precision P = 5.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!    For N = 3, 5 and 6, there are two versions of the rule, chosen by setting
!    the OPTION variable to 1 or 2.
!
!    Versions of this rule are only available for N = 2 through 7.
!
!    There is a typographical error in the reference.
!    For the second version of the rule for N = 2, the line
!      gamma =    0.313300683022281D+00
!    should read
!      gamma =    0.312200683022281D+00
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    2 <= N <= 7.
!
!    Input, integer ( kind = 4 ) OPTION, selects option 1 or 2.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind=Rkind ) X(N,O), the abscissas.
!
!    Output, real ( kind=Rkind ) W(O), the weights.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind=Rkind ) a
  real ( kind=Rkind ) b
  real ( kind=Rkind ) c
  real ( kind=Rkind ) eta
  real ( kind=Rkind ) gamma
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind=Rkind ) lambda
  real ( kind=Rkind ) mu
  integer ( kind = 4 ) option
  real ( kind=Rkind ) volume
  real ( kind=Rkind ) w(o)
  real ( kind=Rkind ) x(n,o)
  real ( kind=Rkind ) xsi

  if ( n < 2 .or. 7 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_05_1 - Fatal error!'
    write ( *, '(a)' ) '  2 <= N <= 7 required.'
    stop
  end if

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_05_1 - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    stop
  end if

  if ( option == 2 ) then
    if ( n /= 3 .and. n /= 5 .and. n /= 6 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EN_R2_05_1 - Fatal error!'
      write ( *, '(a)' ) '  OPTION = 2 requires N = 3, 5 or 6.'
      stop
    end if
  end if

  volume = sqrt ( pi ** n )

  if ( n == 2 ) then
    eta =      0.446103183094540_Rkind
    lambda =   0.136602540378444_Rkind * TEN
    xsi =    - 0.366025403784439_Rkind
    mu =       0.198167882945871_Rkind * TEN
    gamma =    0.000000000000000_Rkind
    a =        0.328774019778636_Rkind * volume
    b =        0.833333333333333_Rkind / TEN * volume
    c =        0.455931355469736_Rkind / TEN**2 * volume
  else if ( n == 3 .and. option == 1 ) then
    eta =      0.476731294622796_Rkind
    lambda =   0.935429018879534_Rkind
    xsi =    - 0.731237647787132_Rkind
    mu =       0.433155309477649_Rkind
    gamma =    0.266922328697744_Rkind * TEN
    a =        0.242000000000000_Rkind * volume
    b =        0.810000000000000_Rkind / TEN * volume
    c =        0.500000000000000_Rkind / TEN**2 * volume
!
!  The value of gamma that follows corrects an error in the reference.
!
  else if ( n == 3 .and. option == 2 ) then
    eta =      0.476731294622796_Rkind
    lambda =   0.128679320334269_Rkind * TEN
    xsi =    - 0.379873463323979_Rkind
    mu =     - 0.192386729447751_Rkind * TEN
    gamma =    0.312200683022281_Rkind
    a =        0.242000000000000_Rkind * volume
    b =        0.810000000000000_Rkind / TEN * volume
    c =        0.500000000000000_Rkind / TEN**2 * volume
  else if ( n == 4 ) then
    eta =      0.523945658287507_Rkind
    lambda =   0.119433782552719_Rkind * TEN
    xsi =    - 0.398112608509063_Rkind
    mu =     - 0.318569372920112_Rkind
    gamma =    0.185675837424096_Rkind * TEN
    a =        0.155502116982037_Rkind * volume
    b =        0.777510584910183_Rkind / TEN * volume
    c =        0.558227484231506_Rkind / TEN**2 * volume
  else if ( n == 5 .and. option == 1 ) then
    eta =      0.214972564378798_Rkind * TEN
    lambda =   0.464252986016289_Rkind * TEN
    xsi =    - 0.623201054093728_Rkind
    mu =     - 0.447108700673434_Rkind
    gamma =    0.812171426076311_Rkind
    a =        0.487749259189752_Rkind / TEN**3 * volume
    b =        0.487749259189752_Rkind / TEN**3 * volume
    c =        0.497073504444862_Rkind / TEN * volume
  else if ( n == 5 .and. option == 2 ) then
    eta =      0.615369528365158_Rkind
    lambda =   0.132894698387445_Rkind * TEN
    xsi =    - 0.178394363877324_Rkind
    mu =     - 0.745963266507289_Rkind
    gamma =    0.135503972310817_Rkind * TEN
    a =        0.726415024414905_Rkind / TEN * volume
    b =        0.726415024414905_Rkind / TEN * volume
    c =        0.641509853510569_Rkind / TEN**2 * volume
  else if ( n == 6 .and. option == 1 ) then
    eta =      0.100000000000000_Rkind * TEN
    lambda =   0.141421356237309_Rkind * TEN
    xsi =      0.000000000000000_Rkind
    mu =     - 0.100000000000000_Rkind * TEN
    gamma =    0.100000000000000_Rkind * TEN
    a =        0.781250000000000_Rkind / TEN**2 * volume
    b =        0.625000000000000_Rkind / TEN * volume
    c =        0.781250000000000_Rkind / TEN**2 * volume
  else if ( n == 6 .and. option == 2 ) then
    eta =      ONE
    lambda =   0.942809041582063_Rkind
    xsi =    - 0.471404520791032_Rkind
    mu =     - 0.166666666666667_Rkind * TEN
    gamma =    0.333333333333333_Rkind
    a =        0.781250000000000_Rkind / TEN**2 * volume
    b =        0.625000000000000_Rkind / TEN * volume
    c =        0.781250000000000_Rkind / TEN**2 * volume
  else if ( n == 7 ) then
    eta =      ZERO
    lambda =   0.959724318748357_Rkind
    xsi =    - 0.772326488820521_Rkind
    mu =     - 0.141214270131942_Rkind * TEN
    gamma =    0.319908106249452_Rkind
    a =        0.111111111111111_Rkind * volume
    b =        0.138888888888889_Rkind / TEN * volume
    c =        0.138888888888889_Rkind / TEN * volume
  end if

  x(1:n,1:o) = ZERO

  k = 0
!
!  2 points.
!
  k = k + 1
  x(1:n,k) = - eta
  w(k) = a
  k = k + 1
  x(1:n,k) = + eta
  w(k) = a
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(1:n,k) = - xsi
    x(i,k) = - lambda
    w(k) = b
    k = k + 1
    x(1:n,k) = + xsi
    x(i,k) = + lambda
    w(k) = b
  end do
!
!  2 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(1:n,k) = - gamma
      x(i,k) = - mu
      x(j,k) = - mu
      w(k) = c
      k = k + 1
      x(1:n,k) = + gamma
      x(i,k) = + mu
      x(j,k) = + mu
      w(k) = c
    end do
  end do

  return
end
subroutine en_r2_05_1_size ( n, option, o)

!*****************************************************************************80
!
!! EN_R2_05_1_SIZE sizes the Stroud rule 5.1 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = N^2 + N + 2.
!
!    The rule has precision P = 5.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!    For N = 3, 5 and 6, there are two versions of the rule, chosen by setting
!    the OPTION variable to 1 or 2.
!
!    Versions of this rule are only available for N = 2 through 7.
!
!    There is a typographical error in the reference.
!    For the second version of the rule for N = 2, the line
!      gamma =    0.313300683022281D+00
!    should read
!      gamma =    0.312200683022281D+00
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    2 <= N <= 7.
!
!    Input, integer ( kind = 4 ) OPTION, selects option 1 or 2.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o
  integer ( kind = 4 ) option
  integer ( kind = 4 ) err

  o = -1

  if ( n < 2 .or. 7 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_05_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  2 <= N <= 7 required.'
    RETURN
  end if

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_05_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    RETURN
  end if

  if ( option == 2 ) then
    if ( n /= 3 .and. n /= 5 .and. n /= 6 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EN_R2_05_1_SIZE - Fatal error!'
      write ( *, '(a)' ) '  OPTION = 2 requires N = 3, 5 or 6.'
      RETURN
    end if
  end if

  o = n * n + n + 2

  return
end
subroutine en_r2_05_2 ( n, o, x, w )

!*****************************************************************************80
!
!! EN_R2_05_2 implements the Stroud rule 5.2 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 2 * N^2 + 1.
!
!    The rule has precision P = 5.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind=Rkind ) X(N,O), the abscissas.
!
!    Output, real ( kind=Rkind ) W(O), the weights.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind=Rkind ) a
  real ( kind=Rkind ) b
  real ( kind=Rkind ) c
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  real ( kind=Rkind ) r
  real ( kind=Rkind ) s
  real ( kind=Rkind ) volume
  real ( kind=Rkind ) w(o)
  real ( kind=Rkind ) x(n,o)

  volume = sqrt ( pi ** n )

  a = TWO * volume / real ( n + 2, kind=Rkind )
  b = real ( 4 - n, kind=Rkind ) * volume / TWO &
    / real ( ( n + 2 ) * ( n + 2 ), kind=Rkind )
  c = volume / real ( ( n + 2 ) * ( n + 2 ), kind=Rkind )

  r = sqrt ( real ( n + 2, kind=Rkind ) / TWO )
  s = sqrt ( real ( n + 2, kind=Rkind ) / FOUR )

  x(1:n,1:o) = ZERO

  k = 0
!
!  1 point.
!
  k = k + 1
! x(1:n,k) = ZERO
  w(k) = a
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - r
    w(k) = b
    k = k + 1
    x(i,k) = + r
    w(k) = b
  end do
!
!  4 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - s
      x(j,k) = - s
      w(k) = c
      k = k + 1
      x(i,k) = - s
      x(j,k) = + s
      w(k) = c
      k = k + 1
      x(i,k) = + s
      x(j,k) = - s
      w(k) = c
      k = k + 1
      x(i,k) = + s
      x(j,k) = + s
      w(k) = c
    end do
  end do

  return
end
subroutine en_r2_05_2_size ( n, o )

!*****************************************************************************80
!
!! EN_R2_05_2_SIZE sizes the Stroud rule 5.2 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 2 * N^2 + 1.
!
!    The rule has precision P = 5.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 2 * n * n + 1

  return
end
subroutine en_r2_05_3 ( n, o, x, w )

!*****************************************************************************80
!
!! EN_R2_05_3 implements the Stroud rule 5.3 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 2^N + 2 * N.
!
!    The rule has precision P = 5.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!    The rule requires 3 <= N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind=Rkind ) X(N,O), the abscissas.
!
!    Output, real ( kind=Rkind ) W(O), the weights.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind=Rkind ) a
  real ( kind=Rkind ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) k
  logical more
  real ( kind=Rkind ) r
  real ( kind=Rkind ) s
  real ( kind=Rkind ) volume
  real ( kind=Rkind ) w(o)
  real ( kind=Rkind ) x(n,o)

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_05_3 - Fatal error!'
    write ( *, '(a)' ) '  3 <= N is required.'
    stop
  end if

  volume = sqrt ( pi ** n )

  a = FOUR * volume / real ( ( n + 2 ) * ( n + 2 ), kind=Rkind )
  b = real ( ( n - 2 ) * ( n - 2 ), kind=Rkind ) * volume &
    / real ( 2**n, kind=Rkind ) / real ( ( n + 2 ) * ( n + 2 ), kind=Rkind )
  r = sqrt ( real ( n + 2, kind=Rkind ) / FOUR )
  s = sqrt ( real ( n + 2, kind=Rkind ) / TWO / real ( n - 2, kind=Rkind ) )

  x(1:n,1:o) = ZERO

  k = 0
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - r
    w(k) = a
    k = k + 1
    x(i,k) = + r
    w(k) = a
  end do
!
!  2^N points.
!
  k = k + 1
  x(1:n,k) = - s
  w(k) = b
  more = .true.
  do while ( more )
    more = .false.
    do i = n, 1, -1
      if ( x(i,k) < ZERO ) then
        k = k + 1
        x(1:i-1,k) = x(1:i-1,k-1)
        x(i,k)     = + s
        x(i+1:n,k) = - s
        w(k) = b
        more = .true.
        exit
      end if
    end do
  end do

  return
end
subroutine en_r2_05_3_size ( n, o )

!*****************************************************************************80
!
!! EN_R2_05_3_SIZE sizes the Stroud rule 5.3 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 2^N + 2 * N.
!
!    The rule has precision P = 5.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!    The rule requires 3 <= N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = -1

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_05_3_SIZE - Fatal error!'
    write ( *, '(a)' ) '  3 <= N is required.'
    RETURN
  end if

  o = 2**n + 2 * n

  return
end
subroutine en_r2_05_4 ( n, o, x, w )

!*****************************************************************************80
!
!! EN_R2_05_4 implements the Stroud rule 5.4 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 2^(N+1) - 1.
!
!    The rule has precision P = 5.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind=Rkind ) X(N,O), the abscissas.
!
!    Output, real ( kind=Rkind ) W(O), the weights.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind=Rkind ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical more
  real ( kind=Rkind ) r
  real ( kind=Rkind ) s
  real ( kind=Rkind ) volume
  real ( kind=Rkind ) w(o)
  real ( kind=Rkind ) x(n,o)

  volume = sqrt ( pi ** n )

  s = sqrt ( HALF )

  x(1:n,1:o) = ZERO

  k = 0
!
!  2^N + 2^(N-1) + 2^(N-2) + ... + 1 = 2^(N+1)-1 points.
!  but do the last point separately.
!
  do i = 1, n

    r = sqrt ( real ( i + 2, kind=Rkind ) / TWO )
    b = TWO ** ( i - n ) * volume / real ( i + 1, kind=Rkind ) &
      / real ( i + 2, kind=Rkind )

    k = k + 1
    x(i,k) = - r
    x(i+1:n,k) = - s
    w(k) = b
    more = .true.
    do while ( more )
      more = .false.
      do j = n, i, -1
        if ( x(j,k) < ZERO ) then
          k = k + 1
          x(1:j-1,k) = x(1:j-1,k-1)
          x(j,k)     =   abs ( x(j,k) )
          x(j+1:n,k) = - abs ( x(j+1:n,k) )
          w(k) = b
          more = .true.
          exit
        end if
      end do
    end do

  end do
!
!  Last point.
!
  k = k + 1
! x(1:n,k) = ZERO
  w(k) = TWO * volume / real ( n + 2, kind=Rkind )

  return
end
subroutine en_r2_05_4_size ( n, o )

!*****************************************************************************80
!
!! EN_R2_05_4_SIZE sizes the Stroud rule 5.4 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 2^(N+1) - 1.
!
!    The rule has precision P = 5.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = 2 ** ( n + 1 ) - 1

  return
end
subroutine en_r2_05_5 ( n, o, x, w )

!*****************************************************************************80
!
!! EN_R2_05_5 implements the Stroud rule 5.5 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = N * 2^N + 1.
!
!    The rule has precision P = 5.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!    There is a second version of this rule however it results in
!    complex abscissas, and so it has been disabled.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    21 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind=Rkind ) X(N,O), the abscissas.
!
!    Output, real ( kind=Rkind ) W(O), the weights.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind=Rkind ) a
  real ( kind=Rkind ) b
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical more
  real ( kind=Rkind ) n_r8
  integer ( kind = 4 ) option
  real ( kind=Rkind ) r
  real ( kind=Rkind ) s
  real ( kind=Rkind ) volume
  real ( kind=Rkind ) w(o)
  real ( kind=Rkind ) x(n,o)

  volume = sqrt ( pi ** n )

  n_r8 = real ( n, kind=Rkind )

  a = TWO * volume / ( n_r8 + TWO )
  b =           volume / ( n_r8 + TWO ) / ( TWO ** n )

  option = 1

  if ( option == 1 ) then
    r = sqrt ( ( n_r8 + TWO &
      + ( n_r8 - ONE ) * sqrt ( TWO * ( n_r8 + TWO ) ) ) &
      / TWO / n_r8 )
    s = sqrt ( ( n_r8 + TWO &
      -                      sqrt ( TWO * ( n_r8 + TWO ) ) ) &
      / TWO / n_r8 )
  else if ( option == 2 ) then
    r = sqrt ( ( n_r8 + TWO &
      - ( n_r8 - ONE ) * sqrt ( TWO * ( n_r8 + TWO ) ) ) &
      / TWO / n_r8 )
    s = sqrt ( ( n_r8 + TWO &
      +                      sqrt ( TWO * ( n_r8 + TWO ) ) ) &
      / TWO / n_r8 )
  end if

  x(1:n,1:o) = ZERO

  k = 0
!
!  1 point.
!
  k = k + 1
! x(1:n,k) = ZERO
  w(k) = a
!
!  N * 2^N points:
!  N choices for location of R, 2^N choices of sign pattern.
!
  do i = 1, n

    k = k + 1
    x(1:n,k) = - s
    x(i,k)   = - r
    w(k) = b

    more = .true.

    do while ( more )
      more = .false.
      do j = n, 1, -1
        if ( x(j,k) < ZERO ) then
          k = k + 1
          x(1:j-1,k) = x(1:j-1,k-1)
          x(j,k)     =   abs ( x(j,k) )
          x(j+1:n,k) = - abs ( x(j+1:n,k) )
          w(k) = b
          more = .true.
          exit
        end if
      end do
    end do

  end do

  return
end
subroutine en_r2_05_5_size ( n, o )

!*****************************************************************************80
!
!! EN_R2_05_5_SIZE sizes the Stroud rule 5.5 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = N * 2^N + 1.
!
!    The rule has precision P = 5.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!    There is a second version of this rule however it results in
!    complex abscissas, and so it has been disabled.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = n * 2 ** n + 1

  return
end
subroutine en_r2_05_6 ( n, o, x, w )

!*****************************************************************************80
!
!! EN_R2_05_6 implements the Stroud rule 5.6 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = ( N + 1 ) * 2^N.
!
!    The rule has precision P = 5.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!    The rule requires 5 <= N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    22 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    5 <= N.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind=Rkind ) X(N,O), the abscissas.
!
!    Output, real ( kind=Rkind ) W(O), the weights.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind=Rkind ) a
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical more
  real ( kind=Rkind ) n_r8
  integer ( kind = 4 ) option
  real ( kind=Rkind ) r
  real ( kind=Rkind ) s
  real ( kind=Rkind ) t
  real ( kind=Rkind ) volume
  real ( kind=Rkind ) w(o)
  real ( kind=Rkind ) x(n,o)

  if ( n < 5 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_05_6 - Fatal error!'
    write ( *, '(a)' ) '  5 <= N is required.'
    stop
  end if

  volume = sqrt ( pi ** n )

  n_r8 = real ( n, kind=Rkind )

  a = volume / ( TWO ** n ) / ( n_r8 + ONE )

  r = sqrt ( ( n_r8 - sqrt ( TWO ) &
    + ( n_r8 - ONE ) * sqrt ( TWO * ( n_r8 + ONE ) ) ) &
    / TWO / n_r8 )
  s = sqrt ( ( n_r8 - sqrt ( TWO ) &
    -                      sqrt ( TWO * ( n_r8 + ONE ) ) ) &
    / TWO / n_r8 )
  t = sqrt ( ( ONE + sqrt ( TWO ) ) / TWO )

  x(1:n,1:o) = ZERO

  k = 0
!
!  N * 2^N points.
!
  do i = 1, n

    k = k + 1
    x(1:n,k) = - s
    x(i,k)   = - r
    w(k) = a

    more = .true.

    do while ( more )
      more = .false.
      do j = n, 1, -1
        if ( x(j,k) < ZERO ) then
          k = k + 1
          x(1:j-1,k) = x(1:j-1,k-1)
          x(j,k)     =   abs ( x(j,k) )
          x(j+1:n,k) = - abs ( x(j+1:n,k) )
          w(k) = a
          more = .true.
          exit
        end if
      end do
    end do

  end do
!
!  2^N points.
!
  k = k + 1
  x(1:n,k) = - t
  w(k) = a
  more = .true.
  do while ( more )
    more = .false.
    do j = n, 1, -1
      if ( x(j,k) < ZERO ) then
        k = k + 1
        x(1:j-1,k) = x(1:j-1,k-1)
        x(j,k)     =   abs ( x(j,k) )
        x(j+1:n,k) = - abs ( x(j+1:n,k) )
        w(k) = a
        more = .true.
        exit
      end if
    end do
  end do

  return
end
subroutine en_r2_05_6_size ( n, o )

!*****************************************************************************80
!
!! EN_R2_05_6_SIZE sizes the Stroud rule 5.6 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = ( N + 1 ) * 2^N.
!
!    The rule has precision P = 5.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!    The rule requires 5 <= N.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    5 <= N.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = -1

  if ( n < 5 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_05_6_SIZE - Fatal error!'
    write ( *, '(a)' ) '  5 <= N is required.'
    RETURN
  end if

  o = ( 2 ** n ) * ( n + 1 )

  return
end
subroutine en_r2_07_1 ( n, option, o, x, w )

!*****************************************************************************80
!
!! EN_R2_07_1 implements the Stroud rule 7.1 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 2^N + 2 * N^2 + 1.
!
!    The rule has precision P = 7.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!    There are two versions of the rule, chosen by setting the
!    OPTION variable to 1 or 2.
!
!    Option 1 is only valid for N = 3, 4, 6 or 7.
!    Option 2 is only valid for N = 3 or 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    N = 3, 4, 6 or 7.
!
!    Input, integer ( kind = 4 ) OPTION, chooses rule option 1 or 2.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind=Rkind ) X(N,O), the abscissas.
!
!    Output, real ( kind=Rkind ) W(O), the weights.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind=Rkind ) a
  real ( kind=Rkind ) b
  real ( kind=Rkind ) c
  real ( kind=Rkind ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical more
  real ( kind=Rkind ) n_r8
  integer ( kind = 4 ) option
  real ( kind=Rkind ) r
  real ( kind=Rkind ) s
  real ( kind=Rkind ) t
  real ( kind=Rkind ) volume
  real ( kind=Rkind ) w(o)
  real ( kind=Rkind ) x(n,o)

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_07_1 - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    stop
  end if

  if ( option == 1 ) then
    if ( n /= 3 .and. n /= 4 .and. n /= 6 .and. n /= 7 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EN_R2_07_1 - Fatal error!'
      write ( *, '(a)' ) '  OPTION 1 requires N =  3, 4, 6 or 7.'
      stop
    end if
  end if

  if ( option == 2 ) then
    if ( n /= 3 .and. n /= 4 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EN_R2_07_1 - Fatal error!'
      write ( *, '(a)' ) '  OPTION 2 requires N =  3 or 4.'
      stop
    end if
  end if

  volume = sqrt ( pi ** n )

  n_r8 = real ( n, kind=Rkind )

  if ( option == 1 ) then
    r = sqrt ( ( THREE * ( EIGHT - n_r8 ) - ( n_r8 - TWO ) &
      * sqrt ( THREE * ( EIGHT - n_r8 ) ) ) / TWO / ( FIVE - n_r8 ) )
    s = sqrt ( ( THREE *             n_r8   -          TWO   &
      * sqrt ( THREE * ( EIGHT - n_r8 ) ) ) / TWO &
      / ( THREE * n_r8 - EIGHT ) )
    t = sqrt ( ( SIX + sqrt ( THREE * ( EIGHT - n_r8 ) ) ) / TWO )
  else if ( option == 2 ) then
    r = sqrt ( ( THREE * ( EIGHT - n_r8 ) + ( n_r8 - TWO ) &
      * sqrt ( THREE * ( EIGHT - n_r8 ) ) ) / TWO / ( FIVE - n_r8 ) )
    s = sqrt ( ( THREE *             n_r8   +          TWO   &
      * sqrt ( THREE * ( EIGHT - n_r8 ) ) ) / TWO &
      / ( THREE * n_r8 - EIGHT ) )
    t = sqrt ( ( SIX - sqrt ( THREE * ( EIGHT - n_r8 ) ) ) / TWO )
  end if

  b = ( EIGHT - n_r8 ) * volume / EIGHT / r ** 6
  c = volume / TWO ** ( n + 3 ) / s ** 6
  d = volume / 16.0_Rkind / t ** 6
  a = volume - TWO * n_r8 * b - TWO ** n * c - TWO * n_r8 &
    * ( n_r8 - ONE ) * d

  x(1:n,1:o) = ZERO

  k = 0
!
!  1 point.
!
  k = k + 1
! x(1:n,k) = ZERO
  w(k) = a
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - r
    w(k) = b
    k = k + 1
    x(i,k) = + r
    w(k) = b
  end do
!
!  2^N points.
!
  k = k + 1
  x(1:n,k) = - s
  w(k) = c
  more = .true.
  do while ( more )
    more = .false.
    do i = n, 1, -1
      if ( x(i,k) < ZERO ) then
        k = k + 1
        x(1:i-1,k) = x(1:i-1,k-1)
        x(i,k)     =   abs ( x(i,k) )
        x(i+1:n,k) = - abs ( x(i+1:n,k) )
        w(k) = c
        more = .true.
        exit
      end if
    end do
  end do
!
!  2 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - t
      x(j,k) = - t
      w(k) = d
      k = k + 1
      x(i,k) = - t
      x(j,k) = + t
      w(k) = d
      k = k + 1
      x(i,k) = + t
      x(j,k) = - t
      w(k) = d
      k = k + 1
      x(i,k) = + t
      x(j,k) = + t
      w(k) = d
    end do
  end do

  return
end
subroutine en_r2_07_1_size ( n, option, o )

!*****************************************************************************80
!
!! EN_R2_07_1_SIZE sizes the Stroud rule 7.1 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 2^N + 2 * N^2 + 1.
!
!    The rule has precision P = 7.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!    There are two versions of the rule, chosen by setting the
!    OPTION variable to 1 or 2.
!
!    Option 1 is only valid for N = 3, 4, 6 or 7.
!    Option 2 is only valid for N = 3 or 4.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    N = 3, 4, 6 or 7.
!
!    Input, integer ( kind = 4 ) OPTION, chooses rule option 1 or 2.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o
  integer ( kind = 4 ) option

  o = -1

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_07_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    RETURN
  end if

  if ( option == 1 ) then
    if ( n /= 3 .and. n /= 4 .and. n /= 6 .and. n /= 7 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EN_R2_07_1_SIZE - Fatal error!'
      write ( *, '(a)' ) '  OPTION 1 requires N =  3, 4, 6 or 7.'
      RETURN
    end if
  end if

  if ( option == 2 ) then
    if ( n /= 3 .and. n /= 4 ) then
      write ( *, '(a)' ) ' '
      write ( *, '(a)' ) 'EN_R2_07_1_SIZE - Fatal error!'
      write ( *, '(a)' ) '  OPTION 2 requires N =  3 or 4.'
      RETURN
    end if
  end if

  o = 2 ** n + 2 * n ** 2 + 1

  return
end
subroutine en_r2_07_2 ( n, o, x, w )

!*****************************************************************************80
!
!! EN_R2_07_2 implements the Stroud rule 7.2 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 2^(N+1) + 4 * N^2.
!
!    The rule has precision P = 7.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!    The rule requires 3 <= N.
!
!    The reference has a typographical error in the description of this rule.
!    The formula:
!
!      (t,t,t,...,t)FS
!
!    should read
!
!      (t,t,0,...,0)FS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    07 January 2012
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind=Rkind ) X(N,O), the abscissas.
!
!    Output, real ( kind=Rkind ) W(O), the weights.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind=Rkind ) a1
  real ( kind=Rkind ) a2
  real ( kind=Rkind ) b
  real ( kind=Rkind ) c
  real ( kind=Rkind ) d
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  logical more
  real ( kind=Rkind ) n_r8
  integer ( kind = 4 ) option
  real ( kind=Rkind ) r
  real ( kind=Rkind ) rho1
  real ( kind=Rkind ) rho2
  real ( kind=Rkind ) s
  real ( kind=Rkind ) t
  real ( kind=Rkind ) volume
  real ( kind=Rkind ) w(o)
  real ( kind=Rkind ) x(n,o)

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_07_2 - Fatal error!'
    write ( *, '(a)' ) '  3 <= N is required.'
    stop
  end if

  volume = sqrt ( pi ** n )

  n_r8 = real ( n, kind=Rkind )

  rho1 = sqrt ( ( n_r8 + TWO - sqrt ( TWO * ( n_r8 + TWO ) ) ) &
    / TWO )
  rho2 = sqrt ( ( n_r8 + TWO + sqrt ( TWO * ( n_r8 + TWO ) ) ) &
    / TWO )
  a1 = ( n_r8 + TWO + sqrt ( TWO * ( n_r8 + TWO ) ) ) / TWO &
    / ( n_r8 + TWO )
  a2 = ( n_r8 + TWO - sqrt ( TWO * ( n_r8 + TWO ) ) ) / TWO &
    / ( n_r8 + TWO )

  r = ONE
  s = sqrt ( ONE / n_r8 )
  t = sqrt ( HALF )
  b = ( EIGHT - n_r8 ) * volume / n_r8 / ( n_r8 + TWO ) &
    / ( n_r8 + FOUR )
  c = n_r8 ** 3 * volume / TWO ** n / n_r8 / ( n_r8 + TWO ) &
    / ( n_r8 + FOUR )
  d = FOUR * volume / n_r8 / ( n_r8 + TWO ) / ( n_r8 + FOUR )

  x(1:n,1:o) = ZERO

  k = 0
!
!  2 * 2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - rho1 * r
    w(k) = a1 * b
    k = k + 1
    x(i,k) = - rho2 * r
    w(k) = a2 * b
    k = k + 1
    x(i,k) = + rho1 * r
    w(k) = a1 * b
    k = k + 1
    x(i,k) = + rho2 * r
    w(k) = a2 * b
  end do
!
!  2 * 2^N points.
!
  k = k + 1
  x(1:n,k) = - rho1 * s
  w(k) = a1 * c
  k = k + 1
  x(1:n,k) = - rho2 * s
  w(k) = a2 * c
  more = .true.
  do while ( more )
    more = .false.
    do i = n, 1, -1
      if ( x(i,k) < ZERO ) then
        k = k + 1
        x(1:i-1,k) =     x(1:i-1,k-2)
        x(i,k)     =   abs ( x(i,k) )
        x(i+1:n,k) = - abs ( x(i+1:n,k) )
        w(k) = a1 * c
        k = k + 1
        x(1:i-1,k) =     x(1:i-1,k-2)
        x(i,k)     =   abs ( x(i,k) )
        x(i+1:n,k) = - abs ( x(i+1:n,k) )
        w(k) = a2 * c
        more = .true.
        exit
      end if
    end do
  end do
!
!  2 * 4 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - rho1 * t
      x(j,k) = - rho1 * t
      w(k) = a1 * d
      k = k + 1
      x(i,k) = - rho1 * t
      x(j,k) = + rho1 * t
      w(k) = a1 * d
      k = k + 1
      x(i,k) = + rho1 * t
      x(j,k) = - rho1 * t
      w(k) = a1 * d
      k = k + 1
      x(i,k) = + rho1 * t
      x(j,k) = + rho1 * t
      w(k) = a1 * d
      k = k + 1
      x(i,k) = - rho2 * t
      x(j,k) = - rho2 * t
      w(k) = a2 * d
      k = k + 1
      x(i,k) = - rho2 * t
      x(j,k) = + rho2 * t
      w(k) = a2 * d
      k = k + 1
      x(i,k) = + rho2 * t
      x(j,k) = - rho2 * t
      w(k) = a2 * d
      k = k + 1
      x(i,k) = + rho2 * t
      x(j,k) = + rho2 * t
      w(k) = a2 * d
    end do
  end do

  return
end
subroutine en_r2_07_2_size ( n, o )

!*****************************************************************************80
!
!! EN_R2_07_2_SIZE sizes the Stroud rule 7.2 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = 2^(N+1) + 4 * N^2.
!
!    The rule has precision P = 7.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!    The rule requires 3 <= N.
!
!    The reference has a typographical error in the description of this rule.
!    The formula:
!
!      (t,t,t,...,t)FS
!
!    should read
!
!      (t,t,0,...,0)FS.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  o = -1

  if ( n < 3 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_07_2_SIZE - Fatal error!'
    write ( *, '(a)' ) '  3 <= N is required.'
    RETURN
  end if

  o = 2 ** ( n + 1 ) + 4 * n * n

  return
end
subroutine en_r2_07_3 ( n, option, o, x, w )

!*****************************************************************************80
!
!! EN_R2_07_3 implements the Stroud rule 7.3 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = ( 4 * N^3 + 8 * N + 3 ) / 3.
!
!    The rule has precision P = 7.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!    There are two versions of each rule, chosen by setting the
!    OPTION variable to 1 or 2.
!
!    The rule as tabulated by Stenger is available for N = 2 through 20.
!    This function accepts N = 3 through 6.
!
!     N    O
!    __  ___
!     3   45
!     4   97
!     5  181
!     6  305
!
!    The reference has a typographical error for N = 5, OPTION 1, B4:
!
!      -(1)0.736330882774831
!
!    should read
!
!      (-1)0.736330882774831
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N <= 6.
!
!    Input, integer ( kind = 4 ) OPTION, chooses rule option 1 or 2.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind=Rkind ) X(N,O), the abscissas.
!
!    Output, real ( kind=Rkind ) W(O), the weights.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind=Rkind ) b0
  real ( kind=Rkind ) b1
  real ( kind=Rkind ) b2
  real ( kind=Rkind ) b3
  real ( kind=Rkind ) b4
  real ( kind=Rkind ) b5
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  logical more
  integer ( kind = 4 ) option
  real ( kind=Rkind ) u
  real ( kind=Rkind ) v
  real ( kind=Rkind ) volume
  real ( kind=Rkind ) w(o)
  real ( kind=Rkind ) x(n,o)

  if ( n < 3 .or. 6 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_07_3 - Fatal error!'
    write ( *, '(a)' ) '  3 <= N <= 6 required.'
    stop
  end if

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_07_3 - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    stop
  end if

  volume = sqrt ( pi ** n )

  if ( n == 3 .and. option == 1 ) then
    u =    0.524647623275290_Rkind
    v =    0.165068012388578_Rkind * TEN
    b0 = - 0.166705761599566_Rkind * TEN**2
    b1 =   0.100296981655678_Rkind * TEN**2
    b2 =   0.161699246687754_Rkind
    b3 = - 0.604719151221535_Rkind * TEN
    b4 =   0.234381399489666_Rkind / TEN
    b5 =   0.417194501880647_Rkind * TEN
  else if ( n == 3 .and. option == 2 ) then
    u =    0.165068012388578_Rkind * TEN
    v =    0.524647623275290_Rkind
    b0 =   0.166705761599566_Rkind * TEN**2
    b1 =   0.178903161957074_Rkind
    b2 = - 0.665808190965810_Rkind * TEN
    b3 =   0.148361823143070_Rkind / TEN
    b4 =   0.229669852539758_Rkind * TEN
    b5 =   0.430097881732984_Rkind / TEN**2
  else if ( n == 4 .and. option == 1 ) then
    u  =   0.524647623275290_Rkind
    v  =   0.165068012388578_Rkind * TEN
    b0 = - 0.167539329651562_Rkind * TEN**3
    b1 =   0.687922329603575_Rkind * TEN**2
    b2 =   0.203518409659014_Rkind
    b3 = - 0.255075279116885_Rkind * TEN**2
    b4 =   0.415430214106084_Rkind / TEN
    b5 =   0.739458001434961_Rkind * TEN
  else if ( n == 4 .and. option == 2 ) then
    u =    0.165068012388578_Rkind * TEN
    v =    0.524647623275290_Rkind
    b0 =   0.688432856406677_Rkind * TEN**2
    b1 =   0.294997847268286_Rkind
    b2 = - 0.199427272118378_Rkind * TEN**2
    b3 =   0.110498755408511_Rkind / TEN
    b4 =   0.407079214570997_Rkind * TEN
    b5 =   0.762328646743931_Rkind / TEN**2
  else if ( n == 5 .and. option == 1 ) then
    u  =   0.524647623275290_Rkind
    v  =   0.165068012388578_Rkind * TEN
    b0 = - 0.826940846964452_Rkind * TEN**3
    b1 =   0.264779097660331_Rkind * TEN**3
    b2 =   0.213460812375320_Rkind
    b3 = - 0.714240197186780_Rkind * TEN**2
    b4 =   0.736330882774831_Rkind / TEN
    b5 =   0.131065518222629_Rkind * TEN**2
  else if ( n == 5 .and. option == 2 ) then
    u =    0.165068012388578_Rkind * TEN
    v =    0.524647623275290_Rkind
    b0 =   0.220502344940121_Rkind * TEN**3
    b1 =   0.537746975313769_Rkind
    b2 = - 0.497781460739792_Rkind * TEN**2
    b3 = - 0.743845245712926_Rkind / TEN**2
    b4 =   0.721529121489956_Rkind * TEN
    b5 =   0.135119234557687_Rkind / TEN
  else if ( n == 6 .and. option == 1 ) then
    u  =   0.524647623275290_Rkind
    v  =   0.165068012388578_Rkind * TEN
    b0 = - 0.309679578630802_Rkind * TEN**4
    b1 =   0.815423321880237_Rkind * TEN**3
    b2 =   0.117326937169073_Rkind
    b3 = - 0.173057295296448_Rkind * TEN**3
    b4 =   0.130511250871491_Rkind
    b5 =   0.232307582494626_Rkind * TEN**2
  else if ( n == 6 .and. option == 2 ) then
    u =    0.165068012388578_Rkind * TEN
    v =    0.524647623275290_Rkind
    b0 =   0.616293651884027_Rkind * TEN**3
    b1 =   0.107529736766179_Rkind * TEN
    b2 = - 0.113807008098269_Rkind * TEN**3
    b3 = - 0.610828352270520_Rkind / TEN
    b4 =   0.127887706992535_Rkind * TEN**2
    b5 =   0.239492607623178_Rkind / TEN
  end if

  x(1:n,1:o) = ZERO

  k = 0
!
!  1 point.
!
  k = k + 1
! x(1:n,k) = ZERO
  w(k) = b0
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - u
    w(k) = b1
    k = k + 1
    x(i,k) = + u
    w(k) = b1
  end do
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - v
    w(k) = b2
    k = k + 1
    x(i,k) = + v
    w(k) = b2
  end do
!
!  4 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - u
      x(j,k) = - u
      w(k) = b3
      k = k + 1
      x(i,k) = - u
      x(j,k) = + u
      w(k) = b3
      k = k + 1
      x(i,k) = + u
      x(j,k) = - u
      w(k) = b3
      k = k + 1
      x(i,k) = + u
      x(j,k) = + u
      w(k) = b3
    end do
  end do
!
!  4 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - v
      x(j,k) = - v
      w(k) = b4
      k = k + 1
      x(i,k) = - v
      x(j,k) = + v
      w(k) = b4
      k = k + 1
      x(i,k) = + v
      x(j,k) = - v
      w(k) = b4
      k = k + 1
      x(i,k) = + v
      x(j,k) = + v
      w(k) = b4
    end do
  end do
!
!  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
!
  do i = 1, n - 2
    do j = i + 1, n - 1
      do l = j + 1, n
        k = k + 1
        x(i,k) = - u
        x(j,k) = - u
        x(l,k) = - u
        w(k) = b5
        k = k + 1
        x(i,k) = - u
        x(j,k) = - u
        x(l,k) = + u
        w(k) = b5
        k = k + 1
        x(i,k) = - u
        x(j,k) = + u
        x(l,k) = - u
        w(k) = b5
        k = k + 1
        x(i,k) = - u
        x(j,k) = + u
        x(l,k) = + u
        w(k) = b5
        k = k + 1
        x(i,k) = + u
        x(j,k) = - u
        x(l,k) = - u
        w(k) = b5
        k = k + 1
        x(i,k) = + u
        x(j,k) = - u
        x(l,k) = + u
        w(k) = b5
        k = k + 1
        x(i,k) = + u
        x(j,k) = + u
        x(l,k) = - u
        w(k) = b5
        k = k + 1
        x(i,k) = + u
        x(j,k) = + u
        x(l,k) = + u
        w(k) = b5
      end do
    end do
  end do

  return
end
subroutine en_r2_07_3_size ( n, option, o )

!*****************************************************************************80
!
!! EN_R2_07_3_SIZE sizes the Stroud rule 7.3 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = ( 4 * N^3 + 8 * N + 3 ) / 3.
!
!    The rule has precision P = 7.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!    There are two versions of each rule, chosen by setting the
!    OPTION variable to 1 or 2.
!
!    The rule as tabulated by Stenger is available for N = 2 through 20.
!    This function accepts N = 3 through 6.
!
!     N    O
!    __  ___
!     3   45
!     4   97
!     5  181
!     6  305
!
!    The reference has a typographical error for N = 5, OPTION 1, B4:
!
!      -(1)0.736330882774831
!
!    should read
!
!      (-1)0.736330882774831
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N <= 6.
!
!    Input, integer ( kind = 4 ) OPTION, chooses rule option 1 or 2.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o
  integer ( kind = 4 ) option

  o = -1

  if ( n < 3 .or. 6 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_07_3_SIZE - Fatal error!'
    write ( *, '(a)' ) '  3 <= N <= 6 required.'
    RETURN
  end if

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_07_3_SIZE - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    RETURN
  end if

  o = ( 4 * n ** 3 + 8 * n + 3 ) / 3

  return
end
subroutine en_r2_09_1 ( n, option, o, x, w )

!*****************************************************************************80
!
!! EN_R2_09_1 implements the Stroud rule 9.1 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = ( 2 * N^4 - 4 * N^3 + 22 * N^2 - 8 * N + 3 ) / 3.
!
!    The rule has precision P = 9.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!    There are two versions of each rule, chosen by setting the
!    OPTION variable to 1 or 2.
!
!    The rule as tabulated by Stenger is available for N = 2 through 20.
!    This function accepts N = 3 through 6.
!
!     N    O
!    __  ___
!     3   77
!     4  193
!     5  421
!     6  825
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N <= 6.
!
!    Input, integer ( kind = 4 ) OPTION, chooses rule option 1 or 2.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind=Rkind ) X(N,O), the abscissas.
!
!    Output, real ( kind=Rkind ) W(O), the weights.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind=Rkind ) b0
  real ( kind=Rkind ) b1
  real ( kind=Rkind ) b2
  real ( kind=Rkind ) b3
  real ( kind=Rkind ) b4
  real ( kind=Rkind ) b5
  real ( kind=Rkind ) b6
  real ( kind=Rkind ) b7
  real ( kind=Rkind ) b8
  integer ( kind = 4 ) i
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  logical more
  integer ( kind = 4 ) option
  real ( kind=Rkind ) u
  real ( kind=Rkind ) v
  real ( kind=Rkind ) volume
  real ( kind=Rkind ) w(o)
  real ( kind=Rkind ) x(n,o)

  if ( n < 3 .or. 6 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_09_1 - Fatal error!'
    write ( *, '(a)' ) '  3 <= N <= 6 required.'
    stop
  end if

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_09_1 - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    stop
  end if

  volume = sqrt ( pi ** n )

  if ( n == 3 ) then
    u =    0.202018287045609_Rkind * TEN
    v =    0.958572464613819_Rkind
    b0 =   0.676448734429924_Rkind
    b1 =   0.511989106291551_Rkind / TEN**2
    b2 =   0.448595723493744_Rkind
    b3 =   0.235223454595606_Rkind / TEN**3
    b4 =   0.915390713080005_Rkind / TEN
    b5 =   0.139208199920793_Rkind / TEN
    b6 =   0.235223454595606_Rkind / TEN**3
    b7 =   0.915390713080008_Rkind / TEN
    b8 =   ZERO
  else if ( n == 4 .and. option == 1 ) then
    u =    0.202018287045609_Rkind * TEN
    v =    0.958572464613819_Rkind
    b0 = - 0.860452945007048_Rkind
    b1 = - 0.405511998533795_Rkind / TEN
    b2 =   0.107026475449715_Rkind * TEN
    b3 =   0.138974239307092_Rkind / TEN**3
    b4 = - 0.162248779448181_Rkind
    b5 =   0.246740110027234_Rkind / TEN
    b6 =   0.138974239307094_Rkind / TEN**3
    b7 =   0.162248779448181_Rkind
    b8 =   0.138974239307094_Rkind / TEN**3
  else if ( n == 4 .and. option == 2 ) then
    u =    0.958572464613819_Rkind
    v =    0.202018287045609_Rkind * TEN
    b0 =   0.265029088766810_Rkind / TEN**2
    b1 =   0.637601342635332_Rkind
    b2 = - 0.394394059389228_Rkind / TEN
    b3 =   0.540829264827264_Rkind / TEN
    b4 = - 0.416922717921281_Rkind / TEN**3
    b5 =   0.246740110027234_Rkind / TEN
    b6 =   0.540829264827270_Rkind / TEN
    b7 =   0.416922717921281_Rkind / TEN**3
    b8 =   0.540829264827269_Rkind / TEN
  else if ( n == 5 .and. option == 1 ) then
    u =    0.202018287045609_Rkind * TEN
    v =    0.958572464613819_Rkind
    b0 = - 0.827347006200826_Rkind * TEN
    b1 = - 0.160820174530905_Rkind
    b2 =   0.353499863758467_Rkind * TEN
    b3 =   0.738976276909564_Rkind / TEN**3
    b4 = - 0.862735421812943_Rkind
    b5 =   0.437335458190621_Rkind / TEN
    b6 = - 0.246325425636523_Rkind / TEN**3
    b7 =   0.287578473937648_Rkind
    b8 =   0.246325425636523_Rkind / TEN**3
  else if ( n == 5 .and. option == 2 ) then
    u =    0.958572464613819_Rkind
    v =    0.202018287045609_Rkind * TEN
    b0 = - 0.624416791055272_Rkind
    b1 =   0.467494915583104_Rkind
    b2 = - 0.152937760910536_Rkind
    b3 =   0.287578473937646_Rkind
    b4 = - 0.221692883072871_Rkind / TEN**2
    b5 =   0.437335458190621_Rkind / TEN
    b6 = - 0.958594913125490_Rkind / TEN
    b7 =   0.738976276909568_Rkind / TEN**3
    b8 =   0.958594913125492_Rkind / TEN
  else if ( n == 6 .and. option == 1 ) then
    u =    0.202018287045609_Rkind * TEN
    v =    0.958572464613819_Rkind
    b0 = - 0.361840434143098_Rkind * TEN**2
    b1 = - 0.447936529138517_Rkind
    b2 =   0.112077863004144_Rkind * TEN**2
    b3 =   0.392940404320855_Rkind / TEN**2
    b4 = - 0.254859786784158_Rkind * TEN
    b5 =   0.775156917007496_Rkind / TEN
    b6 = - 0.130980134773619_Rkind / TEN**2
    b7 =   0.509719573568315_Rkind
    b8 =   0.436600449245395_Rkind / TEN**3
  else if ( n == 6 .and. option == 2 ) then
    u =    0.958572464613819_Rkind
    v =    0.202018287045609_Rkind * TEN
    b0 =   0.448873836333650_Rkind * TEN
    b1 = - 0.238473566140736_Rkind * TEN
    b2 = - 0.413008493198885_Rkind
    b3 =   0.152915872070494_Rkind * TEN
    b4 = - 0.654900673868093_Rkind / TEN**2
    b5 =   0.775156917007496_Rkind / TEN
    b6 = - 0.509719573568314_Rkind
    b7 =   0.130980134773618_Rkind / TEN**2
    b8 =   0.169906524522772_Rkind
  end if

  x(1:n,1:o) = ZERO

  k = 0
!
!  1 point.
!
  k = k + 1
! x(1:n,k) = ZERO
  w(k) = b0
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - u
    w(k) = b1
    k = k + 1
    x(i,k) = + u
    w(k) = b1
  end do
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - v
    w(k) = b2
    k = k + 1
    x(i,k) = + v
    w(k) = b2
  end do
!
!  4 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - u
      x(j,k) = - u
      w(k) = b3
      k = k + 1
      x(i,k) = - u
      x(j,k) = + u
      w(k) = b3
      k = k + 1
      x(i,k) = + u
      x(j,k) = - u
      w(k) = b3
      k = k + 1
      x(i,k) = + u
      x(j,k) = + u
      w(k) = b3
    end do
  end do
!
!  4 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - v
      x(j,k) = - v
      w(k) = b4
      k = k + 1
      x(i,k) = - v
      x(j,k) = + v
      w(k) = b4
      k = k + 1
      x(i,k) = + v
      x(j,k) = - v
      w(k) = b4
      k = k + 1
      x(i,k) = + v
      x(j,k) = + v
      w(k) = b4
    end do
  end do
!
!  4 * ( N * ( N - 1 ) ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - u
      x(j,k) = - v
      w(k) = b5
      k = k + 1
      x(i,k) = - u
      x(j,k) = + v
      w(k) = b5
      k = k + 1
      x(i,k) = + u
      x(j,k) = - v
      w(k) = b5
      k = k + 1
      x(i,k) = + u
      x(j,k) = + v
      w(k) = b5
      k = k + 1
      x(i,k) = - v
      x(j,k) = - u
      w(k) = b5
      k = k + 1
      x(i,k) = - v
      x(j,k) = + u
      w(k) = b5
      k = k + 1
      x(i,k) = + v
      x(j,k) = - u
      w(k) = b5
      k = k + 1
      x(i,k) = + v
      x(j,k) = + u
      w(k) = b5
    end do
  end do
!
!  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
!
  do i = 1, n - 2
    do j = i + 1, n - 1
      do l = j + 1, n
        k = k + 1
        x(i,k) = - u
        x(j,k) = - u
        x(l,k) = - u
        w(k) = b6
        k = k + 1
        x(i,k) = - u
        x(j,k) = - u
        x(l,k) = + u
        w(k) = b6
        k = k + 1
        x(i,k) = - u
        x(j,k) = + u
        x(l,k) = - u
        w(k) = b6
        k = k + 1
        x(i,k) = - u
        x(j,k) = + u
        x(l,k) = + u
        w(k) = b6
        k = k + 1
        x(i,k) = + u
        x(j,k) = - u
        x(l,k) = - u
        w(k) = b6
        k = k + 1
        x(i,k) = + u
        x(j,k) = - u
        x(l,k) = + u
        w(k) = b6
        k = k + 1
        x(i,k) = + u
        x(j,k) = + u
        x(l,k) = - u
        w(k) = b6
        k = k + 1
        x(i,k) = + u
        x(j,k) = + u
        x(l,k) = + u
        w(k) = b6
      end do
    end do
  end do
!
!  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
!
  do i = 1, n - 2
    do j = i + 1, n - 1
      do l = j + 1, n
        k = k + 1
        x(i,k) = - v
        x(j,k) = - v
        x(l,k) = - v
        w(k) = b7
        k = k + 1
        x(i,k) = - v
        x(j,k) = - v
        x(l,k) = + v
        w(k) = b7
        k = k + 1
        x(i,k) = - v
        x(j,k) = + v
        x(l,k) = - v
        w(k) = b7
        k = k + 1
        x(i,k) = - v
        x(j,k) = + v
        x(l,k) = + v
        w(k) = b7
        k = k + 1
        x(i,k) = + v
        x(j,k) = - v
        x(l,k) = - v
        w(k) = b7
        k = k + 1
        x(i,k) = + v
        x(j,k) = - v
        x(l,k) = + v
        w(k) = b7
        k = k + 1
        x(i,k) = + v
        x(j,k) = + v
        x(l,k) = - v
        w(k) = b7
        k = k + 1
        x(i,k) = + v
        x(j,k) = + v
        x(l,k) = + v
        w(k) = b7
      end do
    end do
  end do
!
!  16 * ( N * ( N - 1 ) * ( N - 2 ) * ( N - 3 ) / 24 ) points.
!
  do i = 1, n - 3
    do j = i + 1, n - 2
      do l = j + 1, n - 1
        do m = l + 1, n
          k = k + 1
          x(i,k) = - u
          x(j,k) = - u
          x(l,k) = - u
          x(m,k) = - u
          w(k) = b8
          k = k + 1
          x(i,k) = - u
          x(j,k) = - u
          x(l,k) = - u
          x(m,k) = + u
          w(k) = b8
          k = k + 1
          x(i,k) = - u
          x(j,k) = - u
          x(l,k) = + u
          x(m,k) = - u
          w(k) = b8
          k = k + 1
          x(i,k) = - u
          x(j,k) = - u
          x(l,k) = + u
          x(m,k) = + u
          w(k) = b8
          k = k + 1
          x(i,k) = - u
          x(j,k) = + u
          x(l,k) = - u
          x(m,k) = - u
          w(k) = b8
          k = k + 1
          x(i,k) = - u
          x(j,k) = + u
          x(l,k) = - u
          x(m,k) = + u
          w(k) = b8
          k = k + 1
          x(i,k) = - u
          x(j,k) = + u
          x(l,k) = + u
          x(m,k) = - u
          w(k) = b8
          k = k + 1
          x(i,k) = - u
          x(j,k) = + u
          x(l,k) = + u
          x(m,k) = + u
          w(k) = b8
          k = k + 1
          x(i,k) = + u
          x(j,k) = - u
          x(l,k) = - u
          x(m,k) = - u
          w(k) = b8
          k = k + 1
          x(i,k) = + u
          x(j,k) = - u
          x(l,k) = - u
          x(m,k) = + u
          w(k) = b8
          k = k + 1
          x(i,k) = + u
          x(j,k) = - u
          x(l,k) = + u
          x(m,k) = - u
          w(k) = b8
          k = k + 1
          x(i,k) = + u
          x(j,k) = - u
          x(l,k) = + u
          x(m,k) = + u
          w(k) = b8
          k = k + 1
          x(i,k) = + u
          x(j,k) = + u
          x(l,k) = - u
          x(m,k) = - u
          w(k) = b8
          k = k + 1
          x(i,k) = + u
          x(j,k) = + u
          x(l,k) = - u
          x(m,k) = + u
          w(k) = b8
          k = k + 1
          x(i,k) = + u
          x(j,k) = + u
          x(l,k) = + u
          x(m,k) = - u
          w(k) = b8
          k = k + 1
          x(i,k) = + u
          x(j,k) = + u
          x(l,k) = + u
          x(m,k) = + u
          w(k) = b8
        end do
      end do
    end do
  end do

  return
end
subroutine en_r2_09_1_size ( n, option, o )

!*****************************************************************************80
!
!! EN_R2_09_1_SIZE sizes the Stroud rule 9.1 for region EN_R2.
!
!  Discussion:
!
!    The rule has order O = ( 2 * N^4 - 4 * N^3 + 22 * N^2 - 8 * N + 3 ) / 3.
!
!    The rule has precision P = 9.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!    There are two versions of each rule, chosen by setting the
!    OPTION variable to 1 or 2.
!
!    The rule as tabulated by Stenger is available for N = 2 through 20.
!    This function accepts N = 3 through 6.
!
!     N    O
!    __  ___
!     3   77
!     4  193
!     5  421
!     6  825
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N <= 6.
!
!    Input, integer ( kind = 4 ) OPTION, chooses rule option 1 or 2.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o
  integer ( kind = 4 ) option

  o = -1

  if ( n < 3 .or. 6 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_09_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  3 <= N <= 6 required.'
    RETURN
  end if

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_09_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    RETURN
  end if

  o = ( 2 * n ** 4 - 4 * n ** 3 + 22 * n ** 2 - 8 * n + 3 ) / 3

  return
end
subroutine en_r2_11_1 ( n, option, o, x, w )

!*****************************************************************************80
!
!! EN_R2_11_1 implements the Stroud rule 11.1 for region EN_R2.
!
!  Discussion:
!
!    The rule has order
!
!      O = ( 4 * N^5 - 20 * N^4 + 140 * N^3 - 130 * N^2 + 96 * N + 15 ) / 15.
!
!    The rule has precision P = 11.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!    There are two versions of each rule, chosen by setting the
!    OPTION variable to 1 or 2.
!
!    The rule as tabulated by Stenger is available for N = 2 through 20.
!    This function accepts N = 3 through 5.
!
!     N    O
!    __  ___
!     3  151
!     4  417
!     5  983
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    23 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N <= 5.
!
!    Input, integer ( kind = 4 ) OPTION, chooses rule option 1 or 2.
!
!    Input, integer ( kind = 4 ) O, the order.
!
!    Output, real ( kind=Rkind ) X(N,O), the abscissas.
!
!    Output, real ( kind=Rkind ) W(O), the weights.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o

  real ( kind=Rkind ) b0
  real ( kind=Rkind ) b1
  real ( kind=Rkind ) b2
  real ( kind=Rkind ) b3
  real ( kind=Rkind ) b4
  real ( kind=Rkind ) b5
  real ( kind=Rkind ) b6
  real ( kind=Rkind ) b7
  real ( kind=Rkind ) b8
  real ( kind=Rkind ) b9
  real ( kind=Rkind ) b10
  real ( kind=Rkind ) b11
  real ( kind=Rkind ) b12
  real ( kind=Rkind ) b13
  real ( kind=Rkind ) b14
  real ( kind=Rkind ) b15
  integer ( kind = 4 ) i
  integer ( kind = 4 ) i1
  integer ( kind = 4 ) i2
  integer ( kind = 4 ) i3
  integer ( kind = 4 ) i4
  integer ( kind = 4 ) i5
  integer ( kind = 4 ) j
  integer ( kind = 4 ) k
  integer ( kind = 4 ) l
  integer ( kind = 4 ) m
  logical more
  integer ( kind = 4 ) option
  real ( kind=Rkind ) u
  real ( kind=Rkind ) v
  real ( kind=Rkind ) volume
  real ( kind=Rkind ) w2
  real ( kind=Rkind ) w(o)
  real ( kind=Rkind ) x(n,o)

  if ( n < 3 .or. 5 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_11_1 - Fatal error!'
    write ( *, '(a)' ) '  3 <= N <= 5 required.'
    stop
  end if

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_11_1 - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    stop
  end if

  volume = sqrt ( pi ** n )

  if ( n == 3 .and. option == 1 ) then
    u =     0.235060497367449_Rkind * TEN
    v =     0.436077411927617_Rkind
    w2 =    0.133584907401370_Rkind * TEN
    b0 =  - 0.881591029957858_Rkind * TEN
    b1 =  - 0.751996143360650_Rkind / TEN
    b2 =    0.621743189471515_Rkind * TEN
    b3 =    0.241426451456494_Rkind
    b4 =  - 0.120709739276065_Rkind / TEN**2
    b5 =  - 0.427751221210138_Rkind * TEN
    b6 =    0.550169924840163_Rkind / TEN
    b7 =    0.237084999634707_Rkind / TEN
    b8 =  - 0.169791992887741_Rkind / TEN**2
    b9 =  - 0.252266276123350_Rkind / TEN**4
    b10 =   0.326777873717691_Rkind * TEN
    b11 =   0.968469949206802_Rkind / TEN**2
    b12 =   0.789754514877422_Rkind / TEN**3
    b13 =   ZERO
    b14 =   ZERO
    b15 =   ZERO
  else if ( n == 3 .and. option == 2 ) then
    u =     0.235060497367449_Rkind * TEN
    v =     0.133584907401370_Rkind * TEN
    w2 =    0.436077411927617_Rkind
    b0 =  - 0.141214037032900_Rkind * TEN**2
    b1 =  - 0.803730274707282_Rkind / TEN
    b2 =    0.235546545595906_Rkind
    b3 =    0.888123191556611_Rkind * TEN
    b4 =    0.142467131155533_Rkind / TEN**3
    b5 =    0.582993124006494_Rkind / TEN
    b6 =  - 0.561099173155661_Rkind * TEN
    b7 =  - 0.204028691521686_Rkind / TEN**2
    b8 =    0.252880089932256_Rkind / TEN
    b9 =  - 0.814378678627283_Rkind / TEN**4
    b10 =   0.804353953375146_Rkind / TEN**2
    b11 =   0.393451849690453_Rkind * TEN
    b12 =   0.171183493169724_Rkind / TEN**3
    b13 =   ZERO
    b14 =   ZERO
    b15 =   ZERO
  else if ( n == 4 .and. option == 1 ) then
    u =     0.235060497367449_Rkind * TEN
    v =     0.436077411927617_Rkind
    w2 =    0.133584907401370_Rkind * TEN
    b0 =    0.241502736147339_Rkind * TEN**3
    b1 =  - 0.196095938531478_Rkind
    b2 =  - 0.128675737999280_Rkind * TEN**3
    b3 =    0.307568784278696_Rkind
    b4 =  - 0.480908422319460_Rkind / TEN**2
    b5 =    0.698087019367085_Rkind * TEN**2
    b6 =    0.631837143743771_Rkind / TEN
    b7 =    0.392226151971179_Rkind / TEN
    b8 =  - 0.300948471646799_Rkind / TEN**2
    b9 =  - 0.650235306755170_Rkind / TEN**4
    b10 = - 0.386951974646715_Rkind * TEN**2
    b11 =   0.171656829095787_Rkind / TEN
    b12 =   0.139980343116450_Rkind / TEN**2
    b13 =   0.101552487093372_Rkind / TEN**4
    b14 =   0.222435922356439_Rkind * TEN**2
    b15 =   ZERO
  else if ( n == 4 .and. option == 2 ) then
    u =     0.235060497367449_Rkind * TEN
    v =     0.133584907401370_Rkind * TEN
    w2 =    0.436077411927617_Rkind
    b0 =  - 0.151944464736584_Rkind * TEN**3
    b1 =  - 0.223498438689039_Rkind
    b2 =    0.243574919068010_Rkind
    b3 =    0.634373877008693_Rkind * TEN**2
    b4 =  - 0.782065187814018_Rkind / TEN**4
    b5 =    0.911833754536616_Rkind / TEN
    b6 =  - 0.238927288245914_Rkind * TEN**2
    b7 =  - 0.422314408318853_Rkind / TEN**2
    b8 =    0.448218289217760_Rkind / TEN
    b9 =  - 0.138053374667391_Rkind / TEN**3
    b10 =   0.607473265800655_Rkind / TEN**2
    b11 =   0.697375246129742_Rkind * TEN
    b12 =   0.303414841680135_Rkind / TEN**3
    b13 = - 0.314574391771792_Rkind / TEN**5
    b14 =   0.409103498175100_Rkind / TEN**2
    b15 =   ZERO
  else if ( n == 5 .and. option == 1 ) then
    u =     0.235060497367449_Rkind * TEN
    v =     0.436077411927617_Rkind
    w2 =    0.133584907401370_Rkind * TEN
    b0 =    0.255885269311763_Rkind * TEN**4
    b1 =  - 0.439598677491526_Rkind
    b2 =  - 0.106541406144610_Rkind * TEN**4
    b3 =    0.453540909054264_Rkind
    b4 =  - 0.132100905623778_Rkind / TEN
    b5 =    0.418606568954203_Rkind * TEN**3
    b6 =    0.511394563043680_Rkind / TEN
    b7 =    0.645581013845604_Rkind / TEN
    b8 =  - 0.533417277494500_Rkind / TEN**2
    b9 =  - 0.137981626254496_Rkind / TEN**3
    b10 = - 0.147436933189884_Rkind * TEN**3
    b11 =   0.304253807765057_Rkind / TEN
    b12 =   0.248108698207828_Rkind / TEN**2
    b13 =   0.113652094546015_Rkind / TEN**4
    b14 =   0.394257407160391_Rkind * TEN**2
    b15 =   0.331725011358320_Rkind / TEN**5
  else if ( n == 5 .and. option == 2 ) then
    u =     0.235060497367449_Rkind * TEN
    v =     0.133584907401370_Rkind * TEN
    w2 =    0.436077411927617_Rkind
    b0 =  - 0.761305347548192_Rkind * TEN**3
    b1 =  - 0.536360805019297_Rkind
    b2 =    0.110669832078736_Rkind
    b3 =    0.246421088923968_Rkind * TEN**3
    b4 =  - 0.773649327968607_Rkind / TEN**3
    b5 =    0.169088641205970_Rkind
    b6 =  - 0.670700680243651_Rkind * TEN**2
    b7 =  - 0.856090560229205_Rkind / TEN**2
    b8 =    0.794446232770302_Rkind / TEN
    b9 =  - 0.220272863263544_Rkind / TEN**3
    b10 = - 0.373515812228225_Rkind / TEN**2
    b11 =   0.123606544052884_Rkind * TEN**2
    b12 =   0.537788804557843_Rkind / TEN**3
    b13 = - 0.122101861480881_Rkind / TEN**4
    b14 =   0.725117070759373_Rkind / TEN**2
    b15 =   0.331725011358320_Rkind / TEN**5
  end if

  x(1:n,1:o) = ZERO

  k = 0
!
!  1 point.
!
  k = k + 1
! x(1:n,k) = ZERO
  w(k) = b0
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - u
    w(k) = b1
    k = k + 1
    x(i,k) = + u
    w(k) = b1
  end do
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - v
    w(k) = b2
    k = k + 1
    x(i,k) = + v
    w(k) = b2
  end do
!
!  2 * N points.
!
  do i = 1, n
    k = k + 1
    x(i,k) = - w2
    w(k) = b3
    k = k + 1
    x(i,k) = + w2
    w(k) = b3
  end do
!
!  4 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - u
      x(j,k) = - u
      w(k) = b4
      k = k + 1
      x(i,k) = - u
      x(j,k) = + u
      w(k) = b4
      k = k + 1
      x(i,k) = + u
      x(j,k) = - u
      w(k) = b4
      k = k + 1
      x(i,k) = + u
      x(j,k) = + u
      w(k) = b4
    end do
  end do
!
!  4 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - v
      x(j,k) = - v
      w(k) = b5
      k = k + 1
      x(i,k) = - v
      x(j,k) = + v
      w(k) = b5
      k = k + 1
      x(i,k) = + v
      x(j,k) = - v
      w(k) = b5
      k = k + 1
      x(i,k) = + v
      x(j,k) = + v
      w(k) = b5
    end do
  end do
!
!  4 * ( N * ( N - 1 ) / 2 ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - w2
      x(j,k) = - w2
      w(k) = b6
      k = k + 1
      x(i,k) = - w2
      x(j,k) = + w2
      w(k) = b6
      k = k + 1
      x(i,k) = + w2
      x(j,k) = - w2
      w(k) = b6
      k = k + 1
      x(i,k) = + w2
      x(j,k) = + w2
      w(k) = b6
    end do
  end do
!
!  4 * ( N * ( N - 1 ) ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - u
      x(j,k) = - v
      w(k) = b7
      k = k + 1
      x(i,k) = - u
      x(j,k) = + v
      w(k) = b7
      k = k + 1
      x(i,k) = + u
      x(j,k) = - v
      w(k) = b7
      k = k + 1
      x(i,k) = + u
      x(j,k) = + v
      w(k) = b7
      k = k + 1
      x(i,k) = - v
      x(j,k) = - u
      w(k) = b7
      k = k + 1
      x(i,k) = - v
      x(j,k) = + u
      w(k) = b7
      k = k + 1
      x(i,k) = + v
      x(j,k) = - u
      w(k) = b7
      k = k + 1
      x(i,k) = + v
      x(j,k) = + u
      w(k) = b7
    end do
  end do
!
!  4 * ( N * ( N - 1 ) ) points.
!
  do i = 1, n - 1
    do j = i + 1, n
      k = k + 1
      x(i,k) = - u
      x(j,k) = - w2
      w(k) = b8
      k = k + 1
      x(i,k) = - u
      x(j,k) = + w2
      w(k) = b8
      k = k + 1
      x(i,k) = + u
      x(j,k) = - w2
      w(k) = b8
      k = k + 1
      x(i,k) = + u
      x(j,k) = + w2
      w(k) = b8
      k = k + 1
      x(i,k) = - w2
      x(j,k) = - u
      w(k) = b8
      k = k + 1
      x(i,k) = - w2
      x(j,k) = + u
      w(k) = b8
      k = k + 1
      x(i,k) = + w2
      x(j,k) = - u
      w(k) = b8
      k = k + 1
      x(i,k) = + w2
      x(j,k) = + u
      w(k) = b8
    end do
  end do
!
!  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
!
  do i = 1, n - 2
    do j = i + 1, n - 1
      do l = j + 1, n
        k = k + 1
        x(i,k) = - u
        x(j,k) = - u
        x(l,k) = - u
        w(k) = b9
        k = k + 1
        x(i,k) = - u
        x(j,k) = - u
        x(l,k) = + u
        w(k) = b9
        k = k + 1
        x(i,k) = - u
        x(j,k) = + u
        x(l,k) = - u
        w(k) = b9
        k = k + 1
        x(i,k) = - u
        x(j,k) = + u
        x(l,k) = + u
        w(k) = b9
        k = k + 1
        x(i,k) = + u
        x(j,k) = - u
        x(l,k) = - u
        w(k) = b9
        k = k + 1
        x(i,k) = + u
        x(j,k) = - u
        x(l,k) = + u
        w(k) = b9
        k = k + 1
        x(i,k) = + u
        x(j,k) = + u
        x(l,k) = - u
        w(k) = b9
        k = k + 1
        x(i,k) = + u
        x(j,k) = + u
        x(l,k) = + u
        w(k) = b9
      end do
    end do
  end do
!
!  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
!
  do i = 1, n - 2
    do j = i + 1, n - 1
      do l = j + 1, n
        k = k + 1
        x(i,k) = - v
        x(j,k) = - v
        x(l,k) = - v
        w(k) = b10
        k = k + 1
        x(i,k) = - v
        x(j,k) = - v
        x(l,k) = + v
        w(k) = b10
        k = k + 1
        x(i,k) = - v
        x(j,k) = + v
        x(l,k) = - v
        w(k) = b10
        k = k + 1
        x(i,k) = - v
        x(j,k) = + v
        x(l,k) = + v
        w(k) = b10
        k = k + 1
        x(i,k) = + v
        x(j,k) = - v
        x(l,k) = - v
        w(k) = b10
        k = k + 1
        x(i,k) = + v
        x(j,k) = - v
        x(l,k) = + v
        w(k) = b10
        k = k + 1
        x(i,k) = + v
        x(j,k) = + v
        x(l,k) = - v
        w(k) = b10
        k = k + 1
        x(i,k) = + v
        x(j,k) = + v
        x(l,k) = + v
        w(k) = b10
      end do
    end do
  end do
!
!  8 * ( N * ( N - 1 ) * ( N - 2 ) / 6 ) points.
!
  do i = 1, n - 2
    do j = i + 1, n - 1
      do l = j + 1, n
        k = k + 1
        x(i,k) = - w2
        x(j,k) = - w2
        x(l,k) = - w2
        w(k) = b11
        k = k + 1
        x(i,k) = - w2
        x(j,k) = - w2
        x(l,k) = + w2
        w(k) = b11
        k = k + 1
        x(i,k) = - w2
        x(j,k) = + w2
        x(l,k) = - w2
        w(k) = b11
        k = k + 1
        x(i,k) = - w2
        x(j,k) = + w2
        x(l,k) = + w2
        w(k) = b11
        k = k + 1
        x(i,k) = + w2
        x(j,k) = - w2
        x(l,k) = - w2
        w(k) = b11
        k = k + 1
        x(i,k) = + w2
        x(j,k) = - w2
        x(l,k) = + w2
        w(k) = b11
        k = k + 1
        x(i,k) = + w2
        x(j,k) = + w2
        x(l,k) = - w2
        w(k) = b11
        k = k + 1
        x(i,k) = + w2
        x(j,k) = + w2
        x(l,k) = + w2
        w(k) = b11
      end do
    end do
  end do
!
!  8 * ( N * ( N - 1 ) * ( N - 2 ) / 2 ) points.
!
  do i = 1, n - 2
    do j = i + 1, n - 1
      do l = j + 1, n
        k = k + 1
        x(i,k) = - u
        x(j,k) = - u
        x(l,k) = - v
        w(k) = b12
        k = k + 1
        x(i,k) = - u
        x(j,k) = - u
        x(l,k) = + v
        w(k) = b12
        k = k + 1
        x(i,k) = - u
        x(j,k) = + u
        x(l,k) = - v
        w(k) = b12
        k = k + 1
        x(i,k) = - u
        x(j,k) = + u
        x(l,k) = + v
        w(k) = b12
        k = k + 1
        x(i,k) = + u
        x(j,k) = - u
        x(l,k) = - v
        w(k) = b12
        k = k + 1
        x(i,k) = + u
        x(j,k) = - u
        x(l,k) = + v
        w(k) = b12
        k = k + 1
        x(i,k) = + u
        x(j,k) = + u
        x(l,k) = - v
        w(k) = b12
        k = k + 1
        x(i,k) = + u
        x(j,k) = + u
        x(l,k) = + v
        w(k) = b12
        k = k + 1
        x(i,k) = - u
        x(j,k) = - v
        x(l,k) = - u
        w(k) = b12
        k = k + 1
        x(i,k) = - u
        x(j,k) = - v
        x(l,k) = + u
        w(k) = b12
        k = k + 1
        x(i,k) = - u
        x(j,k) = + v
        x(l,k) = - u
        w(k) = b12
        k = k + 1
        x(i,k) = - u
        x(j,k) = + v
        x(l,k) = + u
        w(k) = b12
        k = k + 1
        x(i,k) = + u
        x(j,k) = - v
        x(l,k) = - u
        w(k) = b12
        k = k + 1
        x(i,k) = + u
        x(j,k) = - v
        x(l,k) = + u
        w(k) = b12
        k = k + 1
        x(i,k) = + u
        x(j,k) = + v
        x(l,k) = - u
        w(k) = b12
        k = k + 1
        x(i,k) = + u
        x(j,k) = + v
        x(l,k) = + u
        w(k) = b12
        k = k + 1
        x(i,k) = - v
        x(j,k) = - u
        x(l,k) = - u
        w(k) = b12
        k = k + 1
        x(i,k) = - v
        x(j,k) = - u
        x(l,k) = + u
        w(k) = b12
        k = k + 1
        x(i,k) = - v
        x(j,k) = + u
        x(l,k) = - u
        w(k) = b12
        k = k + 1
        x(i,k) = - v
        x(j,k) = + u
        x(l,k) = + u
        w(k) = b12
        k = k + 1
        x(i,k) = + v
        x(j,k) = - u
        x(l,k) = - u
        w(k) = b12
        k = k + 1
        x(i,k) = + v
        x(j,k) = - u
        x(l,k) = + u
        w(k) = b12
        k = k + 1
        x(i,k) = + v
        x(j,k) = + u
        x(l,k) = - u
        w(k) = b12
        k = k + 1
        x(i,k) = + v
        x(j,k) = + u
        x(l,k) = + u
        w(k) = b12
      end do
    end do
  end do
!
!  16 * ( N * ( N - 1 ) * ( N - 2 ) * ( N - 3 ) / 24 ) points.
!
  do i = 1, n - 3
    do j = i + 1, n - 2
      do l = j + 1, n - 1
        do m = l + 1, n
          k = k + 1
          x(i,k) = - u
          x(j,k) = - u
          x(l,k) = - u
          x(m,k) = - u
          w(k) = b13
          k = k + 1
          x(i,k) = - u
          x(j,k) = - u
          x(l,k) = - u
          x(m,k) = + u
          w(k) = b13
          k = k + 1
          x(i,k) = - u
          x(j,k) = - u
          x(l,k) = + u
          x(m,k) = - u
          w(k) = b13
          k = k + 1
          x(i,k) = - u
          x(j,k) = - u
          x(l,k) = + u
          x(m,k) = + u
          w(k) = b13
          k = k + 1
          x(i,k) = - u
          x(j,k) = + u
          x(l,k) = - u
          x(m,k) = - u
          w(k) = b13
          k = k + 1
          x(i,k) = - u
          x(j,k) = + u
          x(l,k) = - u
          x(m,k) = + u
          w(k) = b13
          k = k + 1
          x(i,k) = - u
          x(j,k) = + u
          x(l,k) = + u
          x(m,k) = - u
          w(k) = b13
          k = k + 1
          x(i,k) = - u
          x(j,k) = + u
          x(l,k) = + u
          x(m,k) = + u
          w(k) = b13
          k = k + 1
          x(i,k) = + u
          x(j,k) = - u
          x(l,k) = - u
          x(m,k) = - u
          w(k) = b13
          k = k + 1
          x(i,k) = + u
          x(j,k) = - u
          x(l,k) = - u
          x(m,k) = + u
          w(k) = b13
          k = k + 1
          x(i,k) = + u
          x(j,k) = - u
          x(l,k) = + u
          x(m,k) = - u
          w(k) = b13
          k = k + 1
          x(i,k) = + u
          x(j,k) = - u
          x(l,k) = + u
          x(m,k) = + u
          w(k) = b13
          k = k + 1
          x(i,k) = + u
          x(j,k) = + u
          x(l,k) = - u
          x(m,k) = - u
          w(k) = b13
          k = k + 1
          x(i,k) = + u
          x(j,k) = + u
          x(l,k) = - u
          x(m,k) = + u
          w(k) = b13
          k = k + 1
          x(i,k) = + u
          x(j,k) = + u
          x(l,k) = + u
          x(m,k) = - u
          w(k) = b13
          k = k + 1
          x(i,k) = + u
          x(j,k) = + u
          x(l,k) = + u
          x(m,k) = + u
          w(k) = b13
        end do
      end do
    end do
  end do
!
!  16 * ( N * ( N - 1 ) * ( N - 2 ) * ( N - 3 ) / 24 ) points.
!
  do i = 1, n - 3
    do j = i + 1, n - 2
      do l = j + 1, n - 1
        do m = l + 1, n
          k = k + 1
          x(i,k) = - v
          x(j,k) = - v
          x(l,k) = - v
          x(m,k) = - v
          w(k) = b14
          k = k + 1
          x(i,k) = - v
          x(j,k) = - v
          x(l,k) = - v
          x(m,k) = + v
          w(k) = b14
          k = k + 1
          x(i,k) = - v
          x(j,k) = - v
          x(l,k) = + v
          x(m,k) = - v
          w(k) = b14
          k = k + 1
          x(i,k) = - v
          x(j,k) = - v
          x(l,k) = + v
          x(m,k) = + v
          w(k) = b14
          k = k + 1
          x(i,k) = - v
          x(j,k) = + v
          x(l,k) = - v
          x(m,k) = - v
          w(k) = b14
          k = k + 1
          x(i,k) = - v
          x(j,k) = + v
          x(l,k) = - v
          x(m,k) = + v
          w(k) = b14
          k = k + 1
          x(i,k) = - v
          x(j,k) = + v
          x(l,k) = + v
          x(m,k) = - v
          w(k) = b14
          k = k + 1
          x(i,k) = - v
          x(j,k) = + v
          x(l,k) = + v
          x(m,k) = + v
          w(k) = b14
          k = k + 1
          x(i,k) = + v
          x(j,k) = - v
          x(l,k) = - v
          x(m,k) = - v
          w(k) = b14
          k = k + 1
          x(i,k) = + v
          x(j,k) = - v
          x(l,k) = - v
          x(m,k) = + v
          w(k) = b14
          k = k + 1
          x(i,k) = + v
          x(j,k) = - v
          x(l,k) = + v
          x(m,k) = - v
          w(k) = b14
          k = k + 1
          x(i,k) = + v
          x(j,k) = - v
          x(l,k) = + v
          x(m,k) = + v
          w(k) = b14
          k = k + 1
          x(i,k) = + v
          x(j,k) = + v
          x(l,k) = - v
          x(m,k) = - v
          w(k) = b14
          k = k + 1
          x(i,k) = + v
          x(j,k) = + v
          x(l,k) = - v
          x(m,k) = + v
          w(k) = b14
          k = k + 1
          x(i,k) = + v
          x(j,k) = + v
          x(l,k) = + v
          x(m,k) = - v
          w(k) = b14
          k = k + 1
          x(i,k) = + v
          x(j,k) = + v
          x(l,k) = + v
          x(m,k) = + v
          w(k) = b14
        end do
      end do
    end do
  end do
!
!  All quintuples UUUUU with 32 sign combinations.
!
  do i1 = 1, n - 4
    do i2 = i1 + 1, n - 3
      do i3 = i2 + 1, n - 2
        do i4 = i3 + 1, n - 1
          do i5 = i4 + 1, n
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = - u
            x(i3,k) = - u
            x(i4,k) = - u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = - u
            x(i3,k) = - u
            x(i4,k) = - u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = - u
            x(i3,k) = - u
            x(i4,k) = + u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = - u
            x(i3,k) = - u
            x(i4,k) = + u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = - u
            x(i3,k) = + u
            x(i4,k) = - u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = - u
            x(i3,k) = + u
            x(i4,k) = - u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = - u
            x(i3,k) = + u
            x(i4,k) = + u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = - u
            x(i3,k) = + u
            x(i4,k) = + u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = + u
            x(i3,k) = - u
            x(i4,k) = - u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = + u
            x(i3,k) = - u
            x(i4,k) = - u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = + u
            x(i3,k) = - u
            x(i4,k) = + u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = + u
            x(i3,k) = - u
            x(i4,k) = + u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = + u
            x(i3,k) = + u
            x(i4,k) = - u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = + u
            x(i3,k) = + u
            x(i4,k) = - u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = + u
            x(i3,k) = + u
            x(i4,k) = + u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = - u
            x(i2,k) = + u
            x(i3,k) = + u
            x(i4,k) = + u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = - u
            x(i3,k) = - u
            x(i4,k) = - u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = - u
            x(i3,k) = - u
            x(i4,k) = - u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = - u
            x(i3,k) = - u
            x(i4,k) = + u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = - u
            x(i3,k) = - u
            x(i4,k) = + u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = - u
            x(i3,k) = + u
            x(i4,k) = - u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = - u
            x(i3,k) = + u
            x(i4,k) = - u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = - u
            x(i3,k) = + u
            x(i4,k) = + u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = - u
            x(i3,k) = + u
            x(i4,k) = + u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = + u
            x(i3,k) = - u
            x(i4,k) = - u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = + u
            x(i3,k) = - u
            x(i4,k) = - u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = + u
            x(i3,k) = - u
            x(i4,k) = + u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = + u
            x(i3,k) = - u
            x(i4,k) = + u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = + u
            x(i3,k) = + u
            x(i4,k) = - u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = + u
            x(i3,k) = + u
            x(i4,k) = - u
            x(i5,k) = + u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = + u
            x(i3,k) = + u
            x(i4,k) = + u
            x(i5,k) = - u
            w(k) = b15
            k = k + 1
            x(i1,k) = + u
            x(i2,k) = + u
            x(i3,k) = + u
            x(i4,k) = + u
            x(i5,k) = + u
            w(k) = b15
          end do
        end do
      end do
    end do
  end do

  return
end
subroutine en_r2_11_1_size ( n, option, o )

!*****************************************************************************80
!
!! EN_R2_11_1_SIZE sizes the Stroud rule 11.1 for region EN_R2.
!
!  Discussion:
!
!    The rule has order
!
!      O = ( 4 * N^5 - 20 * N^4 + 140 * N^3 - 130 * N^2 + 96 * N + 15 ) / 15.
!
!    The rule has precision P = 11.
!
!    EN_R2 is the entire N-dimensional space with weight function
!
!      w(x) = product ( 1 <= i <= n ) ( exp ( - x(i)^2 )
!
!    There are two versions of each rule, chosen by setting the
!    OPTION variable to 1 or 2.
!
!    The rule as tabulated by Stenger is available for N = 2 through 20.
!    This function accepts N = 3 through 5.
!
!     N    O
!    __  ___
!     3  151
!     4  417
!     5  983
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    27 January 2010
!
!  Author:
!
!    John Burkardt
!
!  Reference:
!
!    Arthur Stroud,
!    Approximate Calculation of Multiple Integrals,
!    Prentice Hall, 1971,
!    ISBN: 0130438936,
!    LC: QA311.S85.
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) N, the spatial dimension.
!    3 <= N <= 5.
!
!    Input, integer ( kind = 4 ) OPTION, chooses rule option 1 or 2.
!
!    Output, integer ( kind = 4 ) O, the order.
!
  USE mod_system
  implicit none

  integer ( kind = 4 ) n
  integer ( kind = 4 ) o
  integer ( kind = 4 ) option

  o = -1

  if ( n < 3 .or. 5 < n ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_11_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  3 <= N <= 5 required.'
    RETURN
  end if

  if ( option < 1 .or. 2 < option ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'EN_R2_11_1_SIZE - Fatal error!'
    write ( *, '(a)' ) '  1 <= OPTION <= 2 required.'
    RETURN
  end if

  o = ( 4 * n ** 5 - 20 * n ** 4 + 140 * n ** 3 - 130 * n ** 2 &
    + 96 * n + 15 ) / 15


  return
end
