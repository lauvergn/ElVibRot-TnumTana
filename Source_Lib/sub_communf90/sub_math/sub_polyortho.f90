      FUNCTION Funct_1D(x,i,ntyp,first_i)
      USE mod_system
      IMPLICIT NONE
      real(kind=Rkind)  :: Funct_1D

       real(kind=Rkind), intent(in) :: x
       integer, intent(in)          :: i,first_i,ntyp

       real(kind=Rkind) :: fourier,poly_legendre,poly_Hermite,poly_Hermite_exp

       SELECT CASE (ntyp)
       CASE (15) ! poly
         Funct_1D = x**(i-first_i)
       CASE (6) ! Huffaker
         STOP 'funct v6'
       CASE (7) ! BO
         STOP 'funct v7'

       CASE (8) ! (x-Req)/(c(1)*x+c(2)*Req)
         !c(1)=0 c(2)=1 : Dunham
         !c(1)=1 c(2)=0 : SPF
         !c(1)=0 c(2)=1 : Ogilvie
         STOP 'funct v8'
       CASE (12) ! legendre
         Funct_1D = poly_legendre(x,i+1-first_i,0)
       CASE (26) !fourier en x
         Funct_1D = fourier(x,i+1-first_i)
       CASE (32) !Hermite*Exp
         Funct_1D = poly_Hermite_exp(x,i-first_i)
       CASE (132) !Hermite
         Funct_1D = poly_Hermite(x,i-first_i)

       CASE default ! ERROR: wrong function !
         write(out_unitp,*) ' ERROR wrong function, ntyp',ntyp
         STOP
       END SELECT


       END FUNCTION Funct_1D

!================================================================
!    v15 x^(i-1)
!
!    rq: i=[1,2,3,....]
!================================================================
      SUBROUTINE d0d1d2d3poly(x,d0,d1,d2,d3,i)
      USE mod_system
      IMPLICIT NONE

       real(kind=Rkind) d0,d1,d2,d3
       real(kind=Rkind) x
       integer i

       integer ii

       ii = i-1
       IF (ii .EQ. 0) THEN
         d0 = ONE
         d1 = ZERO
         d2 = ZERO
         d3 = ZERO
       ELSE IF (ii .EQ. 1) THEN
         d0 = x
         d1 = ONE
         d2 = ZERO
         d3 = ZERO
       ELSE IF (ii .EQ. 2) THEN
         d0 = x*x
         d1 = x+x
         d2 = TWO
         d3 = ZERO
       ELSE IF (ii .EQ. 3) THEN
         d2 = x
         d1 = x*d2
         d0 = x*d1

         d1 = THREE * d1
         d2 = SIX * d2
         d3 = SIX
       ELSE
         d3 = x**(ii-3)
         d2 = x * d3
         d1 = x * d2
         d0 = x * d1

         d1 = real(ii,                  kind=Rkind) * d1
         d2 = real(ii * (ii-1),         kind=Rkind) * d2
         d3 = real(ii * (ii-1) * (ii-2),kind=Rkind) * d3
       END IF


       end subroutine d0d1d2d3poly
!================================================================
!    v26 serie de fourier
!================================================================
      SUBROUTINE d0d1d2d3fourier(x,d0,d1,d2,d3,i)
      USE mod_system
      IMPLICIT NONE

       real(kind=Rkind) d0,d1,d2,d3
       real(kind=Rkind) x,xx,c,s
       integer i

       integer ii

!---------------------------------------------------------------------
      real(kind=Rkind) sq2pi,sqpi
!---------------------------------------------------------------------

       sqpi = ONE/sqrt(pi)
       sq2pi = ONE/sqrt(pi+pi)


       ii = i/2
       xx = mod(x*real(ii,kind=Rkind),pi+pi)

       IF (ii .EQ. 0) THEN
         d0 = sq2pi
         d1 = ZERO
         d2 = ZERO
         d3 = ZERO
       ELSE
           s = sin(xx) * sqpi
           c = cos(xx) * sqpi
         IF (mod(i,2) .EQ. 0) THEN
           d0 = s
           d1 =  real(ii,kind=Rkind)      * c
           d2 = -real(ii * ii,kind=Rkind) * d0
           d3 = -real(ii * ii,kind=Rkind) * d1
         ELSE
           d0 =  c
           d1 = -real(ii,kind=Rkind)      * s
           d2 = -real(ii * ii,kind=Rkind) * d0
           d3 = -real(ii * ii,kind=Rkind) * d1
         END IF
       END IF


       end subroutine d0d1d2d3fourier
!================================================================
!    eigenfunctions of a particle in a box [0,pi]
!    sin(k*x) with x E [0,pi]
!================================================================
      SUBROUTINE d0d1d2d3box(x,d0,d1,d2,d3,i)
      USE mod_system
      IMPLICIT NONE

       real(kind=Rkind) d0,d1,d2,d3
       real(kind=Rkind) x,xx,c,s
       integer i


!---------------------------------------------------------------------
      real(kind=Rkind) norme
!---------------------------------------------------------------------

       norme = ONE/sqrt(pi*HALF)


       xx = mod(x*real(i,kind=Rkind),pi+pi)
       s = sin(xx) * norme
       c = cos(xx) * norme

       d0 = s
       d1 =  real(i,kind=Rkind)     * c
       d2 = -real(i * i,kind=Rkind) * d0
       d3 = -real(i * i,kind=Rkind) * d1

       end subroutine d0d1d2d3box
!================================================================
!    v26 serie de fourier
!================================================================
      FUNCTION fourier(x,n1)
      USE mod_system
      IMPLICIT NONE
      real(kind=Rkind) :: fourier

       real(kind=Rkind) x,xx
       integer n1

       integer ii

!---------------------------------------------------------------------
      real(kind=Rkind) sq2pi,sqpi
!---------------------------------------------------------------------

       sqpi = ONE/sqrt(pi)
       sq2pi = ONE/sqrt(pi+pi)


       ii = n1/2
       xx = mod(x*real(ii,kind=Rkind),pi+pi)

       IF (ii == 0) THEN
         fourier = sq2pi
         !write(out_unitp,*) 'n1,ii: cte',n1,ii
       ELSE
         IF (mod(n1,2) == 0) THEN
           fourier = sin(xx) * sqpi
           !write(out_unitp,*) 'n1,ii: sin',n1,ii
         ELSE
           fourier = cos(xx) * sqpi
           !write(out_unitp,*) 'n1,ii: cos',n1,ii
         END IF
       END IF

       end function fourier
!===================================================
!
!   Orthonormal chebychev polynomial, Tn(x)
!    x ( -1 =< x =< 1 )
!
!   Rq: nn=n+1
!
!===================================================
      FUNCTION poly_cheby(x,nn)
      USE mod_system
      IMPLICIT NONE
      real(kind=Rkind) :: poly_cheby

       real(kind=Rkind), intent(in) :: x
       integer, intent(in)          :: nn

       real(kind=Rkind) :: p0,p1,p2
       integer :: i,n


       n=nn-1


       IF (n < 0 .OR. abs(x) > ONE) THEN
         write(out_unitp,*) ' ERROR in  poly_cheby :'
         write(out_unitp,*) ' n : ',n,' et x = ',x
         STOP
       END IF

       IF (n == 0) THEN
         poly_cheby = ONE/sqrt(pi)
       ELSE IF (n == 1) THEN
         poly_cheby = x / sqrt(pi/TWO)
       ELSE
         p0 = ONE
         p1 = x
         DO i=2,n
           p2 = TWO*x* p1 - p0
           p0 = p1
           p1 = p2
         END DO
         poly_cheby = p2 / sqrt(pi/TWO)
       END IF

       end function poly_cheby
!===================================================
!
!   calcule pour les polynomes d'hermite 0 a nb_herm
!   la valeur et les dirivees sur les points de la grille
!
!===================================================
      SUBROUTINE d0poly_cheby_grille(x,d0p,nb,nq)
      USE mod_system
      IMPLICIT NONE

      integer, intent(in) :: nq,nb
      real(kind=Rkind), intent(in) :: x(nq)

      real(kind=Rkind) :: d0p(nq,nb)

      real(kind=Rkind) :: p0(nq),p1(nq),p2(nq)

      integer :: i,k

       IF (nb >= 1) THEN ! first poly (n=0)
         d0p(:,1) = ONE/sqrt(pi)
       END IF

       IF (nb >= 2) THEN ! second poly (n=1)
         d0p(:,2) = x(:)/sqrt(pi)
       END IF

       IF (nb >= 3) THEN ! then the other poly (n>1)
         p0 = ONE
         p1 = x
         DO k=3,nb
           p2 = TWO*x* p1 - p0
           d0p(:,k) = p2 / sqrt(pi/TWO)
           p0 = p1
           p1 = p2
         END DO
       END IF


      end subroutine d0poly_cheby_grille
!===================================================
!
!   calcule pour les polynomes d'hermite 0 a nb_herm
!   la valeur et les dirivees sur les points de la grille
!
!===================================================
      SUBROUTINE d0poly_chebyWeight_grid(x,d0p,nb,nq)
      USE mod_system
      IMPLICIT NONE

      integer, intent(in) :: nq,nb
      real(kind=Rkind), intent(in) :: x(nq)

      real(kind=Rkind) :: d0p(nq,nb)

      real(kind=Rkind) :: p0(nq),p1(nq),p2(nq),sqw(nq)

      integer :: i,k

       sqw(:) = ONE/sqrt(sqrt(1-x*x))
       IF (nb >= 1) THEN ! first poly (n=0)
         d0p(:,1) = ONE/sqrt(pi) * sqw(:)
       END IF

       IF (nb >= 2) THEN ! second poly (n=1)
         d0p(:,2) = x(:)/sqrt(pi/TWO) * sqw(:)
       END IF

       IF (nb >= 3) THEN ! then the other poly (n>1)
         p0 = ONE
         p1 = x
         DO k=3,nb
           p2 = TWO*x* p1 - p0
           d0p(:,k) = p2 / sqrt(pi/TWO) * sqw(:)
           p0 = p1
           p1 = p2
         END DO
       END IF


      end subroutine d0poly_chebyWeight_grid
!===================================================
!
!   calcule la valeur d'un polynome de Legendre l,m
!   pour un x ( -1 =< x =< 1 )
!
!   Rq: lll=l+1
!
!===================================================
      FUNCTION poly_legendre(xx,lll,mmm)
      USE mod_system
      IMPLICIT NONE
      real(kind=Rkind) :: poly_legendre

       real(kind=Rkind) x,xx

       real(kind=Rkind) pmm,somx2,fact,pmmp1,pll,poly,norme2

       integer l,ll,lll,m,mmm,i

       l=lll-1
       m=mmm

       x = xx

       IF (x .GT. ONE) x = TWO-x
       IF (x .LT. -ONE) x = -TWO-x


       IF (m < 0 .OR. l < 0 .OR. abs(x) > ONE) THEN
         write(out_unitp,*) 'mauvais arguments dans poly_legendre :'
         write(out_unitp,*) ' m l : ',m,l,' et x = ',x
         STOP
       END IF

       IF (m > l) THEN
         poly_legendre = ZERO
         RETURN
       END IF

       pmm = ONE

       IF (m .GT. 0) THEN
         somx2 = sqrt(ONE - x*x)
         fact = ONE
         DO i=1,m
           pmm = -pmm*fact*somx2
           fact = fact+TWO
         END DO
       END IF


       IF ( m .EQ. l) THEN
         poly = pmm
       ELSE
         pmmp1 = x*real(2*m+1,kind=Rkind)*pmm
         IF (l .EQ. m+1) THEN
           poly = pmmp1
         ELSE
           DO ll=m+2,l

             pll = (x*real(2*ll-1,kind=Rkind)*pmmp1-real(ll+m-1,kind=Rkind)*pmm) / &
                    real(ll-m,kind=Rkind)
             pmm=pmmp1
             pmmp1=pll
           END DO
           poly = pll
         END IF
       END IF

!      nomalisation de Pl,m(x)
       norme2 = TWO/real(2*l+1,kind=Rkind)
       DO i=l-m+1,l+m
         norme2 = norme2 * real(i,kind=Rkind)
       END DO

!      write(out_unitp,*) l,m,norme2
       poly_legendre = poly/sqrt(norme2)

       RETURN
       end function poly_legendre
!===================================================
!
!   calcule la derivee d'un polynome de Legendre l,m
!   pour un x ( -1 =< x =< 1 )
!
!   (x*x-1)*P'l(x) = l*(x*Pl(x) - Pl-1(x))  et il y la normalisation
! => P'l(x) = l/(x*x-1) * ( x*Pl(x)- sqrt(2l+1/2l-1)*Pl-1(x) )
!
!   Rq: lll=l+1
!
!===================================================
      FUNCTION d1poly_legendre(x,lll)
      USE mod_system
      IMPLICIT NONE
      real(kind=Rkind) :: d1poly_legendre

       real(kind=Rkind) x
       real(kind=Rkind) poly_legendre
       integer l,lll

       l=lll-1

       IF (l .EQ. 0) THEN
         d1poly_legendre = ZERO
       ELSE
         d1poly_legendre = real(l,kind=Rkind)/(x*x-ONE) * (             &
                    x*poly_legendre(x,lll,0) -                          &
                   sqrt(real(2*l+2,kind=Rkind)/real(2*l-1,kind=Rkind))* &
                   poly_legendre(x,lll-1,0) )

!      write(out_unitp,*) 'derive :',lll,x,d1poly_legendre,d1

       END IF

       RETURN
       end function d1poly_legendre
!===================================================
!
!   calcule la derivee d'un polynome de Legendre l,m
!   pour un x ( -1 =< x =< 1 )
!
!   (x*x-1)*P'l(x) = l*(x*Pl(x) - Pl-1(x))  et il y la normalisation
! => P'l(x) = l/(x*x-1) * ( x*Pl(x)- sqrt(2l+1/2l-1)*Pl-1(x) )
!
!   Rq: lll=l+1
!
!===================================================
      SUBROUTINE d0d1d2poly_legendre(x,lll,d0,d1,d2,num,step)
      USE mod_system
      IMPLICIT NONE

       real(kind=Rkind) x
       real(kind=Rkind) poly_legendre

       real(kind=Rkind) ep,em,step,d0,d1,d2,d01
       logical num

       integer l,lll


       l=lll-1

       d0 = poly_legendre(x,lll,0)

       IF (l .EQ. 0) THEN
         d1 = ZERO
         d2 = ZERO
       ELSE
         IF (num) THEN
           ep = poly_legendre(x+step,lll,0)
           em = poly_legendre(x-step,lll,0)
           CALL d1d2_b(d0,ep,em,d1,d2,step)
         ELSE
           d01 = poly_legendre(x,lll-1,0)
           d1 = real(l,kind=Rkind)/(x*x-ONE) * (x*d0 -                  &
           sqrt(real(2*l+1,kind=Rkind)/real(2*l-1,kind=Rkind))*d01 )

           d2 = (-d0 * real(l*(l+1),kind=Rkind) + TWO*x*d1)/(ONE-x*x)
         END IF
       END IF
!      write(out_unitp,*) 'derive :',lll,x,d0,d1,d2

       RETURN
       end subroutine d0d1d2poly_legendre
!===================================================
!
!   calcule la derivee d'un polynome de Legendre l,m
!   pour un x ( -1 =< x =< 1 )
!
!   (x*x-1)*P'l(x) = l*(x*Pl(x) - Pl-1(x))  et il y la normalisation
! => P'l(x) = l/(x*x-1) * ( x*Pl(x)- sqrt(2l+1/2l-1)*Pl-1(x) )
!
!   Rq: lll=l+1
!
!===================================================
      SUBROUTINE d0d1d2d3poly_legendre(                                 &
                                        x,lll,                          &
                                        d0,d1,d2,d3,                    &
                                        nderiv)
      USE mod_system
      IMPLICIT NONE

       real(kind=Rkind) x
       real(kind=Rkind) poly_legendre

       real(kind=Rkind) d0,d1,d2,d01,d3

       integer l,lll,nderiv

       l=lll-1

       d0 = poly_legendre(x,lll,0)

       IF (l .EQ. 0) THEN
         d1 = ZERO
         d2 = ZERO
         d3 = ZERO
       ELSE
         d01 = poly_legendre(x,lll-1,0)
         d1 = real(l,kind=Rkind)/(x*x-ONE) * (x*d0 -                    &
           sqrt(real(2*l+1,kind=Rkind)/real(2*l-1,kind=Rkind))*d01 )

         d2 = (-d0 * real(l*(l+1),kind=Rkind) + TWO*x*d1)/(ONE-x*x)
         d3 = (d1 * real(2-l*(l+1),kind=Rkind) + d2 * FOUR*x ) / (ONE-x*x)

       END IF
!      write(out_unitp,*) 'derive :',lll,x,d0,d1,d2,d3

       RETURN
       end subroutine d0d1d2d3poly_legendre
!===================================================
!
!   calcule la derivee d'un polynome de Legendre l,m
!   pour un c=cos(theta)
!   rq : s=sin(theta)
!   on veut les derivees en theta
!
!   (x*x-1)*P'l(x) = l*(x*Pl(x) - Pl-1(x))  et il y la normalisation
! => P'l(x) = l/(x*x-1) * ( x*Pl(x)- sqrt(2l+1/2l-1)*Pl-1(x) )
!
!   Rq: lll=l+1
!
!===================================================
      SUBROUTINE d0d1d2d3poly_legendre_theta(                           &
                                              c,s,lll,                  &
                                              dc0,dc1,dc2,dc3,          &
                                              nderiv)
      USE mod_system
      IMPLICIT NONE

       real(kind=Rkind) s,c

       real(kind=Rkind) d0,d1,d2,d3
       real(kind=Rkind) dc0,dc1,dc2,dc3

       integer l,lll,nderiv

       l=lll-1


       CALL d0d1d2d3poly_legendre(c,lll,d0,d1,d2,d3,nderiv)
!      write(out_unitp,*) 'derive :',lll,c,d0,d1,d2,d3

!      transfo des derivees en theta

       dc0 = d0
       dc1 = -s * d1
       dc2 = -c * d1 +       s*s * d2
       dc3 =  s * d1 + THREE*s*c * d2 - s*s*s * d3

!      write(out_unitp,*) 'derive :',lll,c,dc0,dc1,dc2,dc3

       RETURN
       end subroutine d0d1d2d3poly_legendre_theta
!===================================================
!
!   calcule pour les polynomes d'hermite 0 a nb_herm
!   la valeur et les dirivees sur les points de la grille
!
!===================================================
      SUBROUTINE d0d1d2poly_legendre_grille(xl,                         &
        d0l,d1l,d2l,nb_legendre,nb_quadra,deriv,num,step)
      USE mod_system
      IMPLICIT NONE

      integer nb_legendre,nb_quadra
      real(kind=Rkind) step
      logical num,deriv

      real(kind=Rkind) xl(nb_quadra)
      real(kind=Rkind) d0l(nb_quadra,nb_legendre)
      real(kind=Rkind) d1l(nb_quadra,nb_legendre)
      real(kind=Rkind) d2l(nb_quadra,nb_legendre)

      integer i,k

!     write(out_unitp,*) 'd0d1d2poly_legendre_grille : deriv',deriv
!     write(out_unitp,*) 'd0d1d2poly_legendre_grille : num',num
!     write(out_unitp,*) 'd0d1d2poly_legendre_grille : step',step
      DO k=1,nb_quadra
        DO i=1,nb_legendre
          CALL d0d1d2poly_legendre(xl(k),i,                             &
                    d0l(k,i),d1l(k,i),d2l(k,i),num,step)
        END DO
      END DO
!     CALL ecrit_dib(d0l,d2l,d2l,nb_legendre,nb_quadra)

      RETURN
      end subroutine d0d1d2poly_legendre_grille
!===================================================
!
!   calcule pour les polynomes d'hermite 0 a nb_herm
!   la valeur et les dirivees sur les points de la grille
!
!===================================================
      SUBROUTINE d0d1d2Plm_grid(xl,                                     &
        d0l,d1l,d2l,nb_legendre,nb_quadra,deriv,num,step)
      USE mod_system
      IMPLICIT NONE

      integer nb_legendre,nb_quadra
      real(kind=Rkind) step
      logical num,deriv

      real(kind=Rkind) xl(nb_quadra)
      real(kind=Rkind) d0l(nb_quadra,nb_legendre)
      real(kind=Rkind) d1l(nb_quadra,nb_legendre)
      real(kind=Rkind) d2l(nb_quadra,nb_legendre)

      integer i,k

!     write(out_unitp,*) 'd0d1d2poly_legendre_grille : deriv',deriv
!     write(out_unitp,*) 'd0d1d2poly_legendre_grille : num',num
!     write(out_unitp,*) 'd0d1d2poly_legendre_grille : step',step
      DO k=1,nb_quadra
        DO i=1,nb_legendre
          CALL d0d1d2poly_legendre(xl(k),i,                             &
                    d0l(k,i),d1l(k,i),d2l(k,i),num,step)
        END DO
      END DO

      RETURN
      end subroutine d0d1d2Plm_grid
!===================================================
!
!   calcule pour les polynomes d'hermite 0 a nb_herm
!   la valeur et les dirivees sur les points de la grille
!
!   polynome paire
!
!===================================================
      SUBROUTINE d0d1d2Plm_0_grid(xl,                                   &
        d0l,d1l,d2l,nb_legendre,nb_quadra,deriv,num,step)
      USE mod_system
      IMPLICIT NONE

      integer nb_legendre,nb_quadra
      real(kind=Rkind) step
      logical num,deriv

      real(kind=Rkind) xl(nb_quadra)
      real(kind=Rkind) d0l(nb_quadra,nb_legendre)
      real(kind=Rkind) d1l(nb_quadra,nb_legendre)
      real(kind=Rkind) d2l(nb_quadra,nb_legendre)

      integer i,ii,k

!     write(out_unitp,*) 'd0d1d2poly_legendre_grille : deriv',deriv
!     write(out_unitp,*) 'd0d1d2poly_legendre_grille : num',num
!     write(out_unitp,*) 'd0d1d2poly_legendre_grille : step',step
      DO k=1,nb_quadra
        ii = 1
        DO i=1,nb_legendre
          CALL d0d1d2poly_legendre(xl(k),ii,                            &
                    d0l(k,i),d1l(k,i),d2l(k,i),num,step)
          ii = ii + 2
        END DO
      END DO

      RETURN
      end subroutine d0d1d2Plm_0_grid
!===================================================
!
!   calcule pour les polynomes d'hermite 0 a nb_herm
!   la valeur et les dirivees sur les points de la grille
!
!   polynome impaire
!
!===================================================
      SUBROUTINE d0d1d2Plm_1_grid(xl,                                   &
        d0l,d1l,d2l,nb_legendre,nb_quadra,deriv,num,step)
      USE mod_system
      IMPLICIT NONE

      integer nb_legendre,nb_quadra
      real(kind=Rkind) step
      logical num,deriv

      real(kind=Rkind) xl(nb_quadra)
      real(kind=Rkind) d0l(nb_quadra,nb_legendre)
      real(kind=Rkind) d1l(nb_quadra,nb_legendre)
      real(kind=Rkind) d2l(nb_quadra,nb_legendre)

      integer i,ii,k

!     write(out_unitp,*) 'd0d1d2poly_legendre_grille : deriv',deriv
!     write(out_unitp,*) 'd0d1d2poly_legendre_grille : num',num
!     write(out_unitp,*) 'd0d1d2poly_legendre_grille : step',step
      DO k=1,nb_quadra
        ii = 2
        DO i=1,nb_legendre
          CALL d0d1d2poly_legendre(xl(k),ii,                            &
                    d0l(k,i),d1l(k,i),d2l(k,i),num,step)
          ii = ii + 2
        END DO
      END DO

      RETURN
      end subroutine d0d1d2Plm_1_grid
!===================================================
!
!   calcule la valeur d'un polynome de Hermite l
!   on rajoute la partie exp
!   pour un x ( -inf =< x =< inf )
!
!===================================================
      FUNCTION poly_Hermite_exp(x,l)
      USE mod_system
      IMPLICIT NONE
      real(kind=Rkind) :: poly_Hermite_exp

      real(kind=Rkind) x
      real(kind=Rkind) poly_Hermite

      integer l

      poly_Hermite_exp = poly_Hermite(x,l) * exp(-x*x*HALF)

!      write(out_unitp,*) x,poly_Hermite_exp
       RETURN
       end function poly_Hermite_exp
!===================================================
!
!   calcule la valeur d'un polynome de Hermite l
!   pour un x ( -inf =< x =< inf )
!
!===================================================
      FUNCTION poly_Hermite(x,l)
      USE mod_system
      IMPLICIT NONE
      real(kind=Rkind) :: poly_Hermite


       real(kind=Rkind) x

!      valeur en x du polynome en l, l-1 et l-2
       real(kind=Rkind) pl0,pl1,pl2,norme


      integer i,l



       IF ( l .LT. 0 ) THEN
         write(out_unitp,*) 'mauvais arguments dans poly_hermite :'
         write(out_unitp,*) ' l < 0 : ',l
         STOP
       END IF

       norme=sqrt(pi)

       IF (l .EQ. 0) THEN
         poly_Hermite = ONE/sqrt(norme)
       ELSE IF (l .EQ. 1) THEN
         norme = norme*TWO
         poly_Hermite = TWO*x/sqrt(norme)
       ELSE

         pl2 = ONE
         pl1 = TWO*x
         norme = norme*TWO

         DO i=2,l
           norme = norme*TWO*real(i,kind=Rkind)
           pl0 = TWO*( x*pl1 - real(i-1,kind=Rkind)*pl2 )
           pl2 = pl1
           pl1 = pl0
         END DO
         poly_Hermite = pl0/sqrt(norme)
       END IF


!      write(out_unitp,*) x,poly_Hermite
       RETURN
       end function poly_Hermite
!===================================================
!
!   calcule la derivee d'un polynome de Hermite (l>0)
!   pour un x ( -inf =< x =< +inf )
!
!    avec la norme sqrt(2**l * l! *sqrt(pi) )
!    Pl'(x) = sqrt(2*l) * Pl-1(x)
!
!===================================================
      SUBROUTINE d0d1d2poly_Hermite(x,l,d0,d1,d2,deriv,num,step)
      USE mod_system
      IMPLICIT NONE

       real(kind=Rkind) x,d0,d1,d2

       real(kind=Rkind) xp,xm,ep,em,step
       logical num,deriv

      integer l

!----- function -----------------------------------
       real(kind=Rkind) :: poly_Hermite
!----- function -----------------------------------

       IF (deriv) THEN

         IF (num) THEN
           xp = x+step
           xm = x-step
           d0 = poly_Hermite(x ,l)
           ep = poly_Hermite(xp,l)
           em = poly_Hermite(xm,l)
           CALL d1d2_b(d0,ep,em,d1,d2,step)
         ELSE
           d0 = poly_Hermite(x,l)
           IF (l .EQ. 0) THEN
             d1 = ZERO
             d2 = ZERO
           ELSE IF (l .EQ. 1) THEN
              d1 = sqrt(TWO)*poly_Hermite(x,0)
              d2 = ZERO
           ELSE
             d1 = sqrt(real(2*l,kind=Rkind)) * poly_Hermite(x,l-1)
             d2 = TWO*(x*d1-d0*real(l,kind=Rkind))
           END IF


         END IF

!        on rajoute la partie exponentielle
!        pexp = exp(-.5d0*x*x)
!        uniquement sur d0
!        pour d1 et d2 : on donne +/- d1/d0 et d2/d0


         d2 = d2/d0
         d1 = d1/d0

         d2 = d2 - d1*d1 -ONE
         d1 = d1 - x
         d0 = d0 * exp(-HALF*x*x)


       ELSE
         d0 = poly_Hermite(x ,l)*exp(-HALF*x*x)
         d1 = ZERO
         d2 = ZERO
       END IF

!      write(out_unitp,*) 'l x d0 d1 d2 hermite :',l,x,d0,d1,d2


       RETURN
       end subroutine d0d1d2poly_Hermite
!===================================================
!
!   calcule la derivee d'un polynome de Hermite (l>0)
!   pour un x ( -inf =< x =< +inf )
!
!    avec la norme sqrt(2**l * l! *sqrt(pi) )
!    Pl'(x) = sqrt(2*l) * Pl-1(x)
!
!===================================================
      SUBROUTINE d0d1d2d3poly_Hermite_exp(x,l,d0,d1,d2,d3,deriv)
      USE mod_system
      IMPLICIT NONE

       real(kind=Rkind) x,d0,d1,d2,d3,pexp

       logical deriv

      integer l

!----- function -----------------------------------
       real(kind=Rkind) poly_Hermite
!----- function -----------------------------------

       IF (deriv) THEN

         d0 = poly_Hermite(x,l)
         IF (l .EQ. 0) THEN
           d1 = ZERO
           d2 = ZERO
           d3 = ZERO
         ELSE IF (l .EQ. 1) THEN
           d1 = sqrt(TWO)*poly_Hermite(x,0)
           d2 = ZERO
           d3 = ZERO
         ELSE IF (l .EQ. 2) THEN
           d1 = sqrt(real(2*l,kind=Rkind)) * poly_Hermite(x,l-1)
           d2 = TWO*(x*d1-d0*real(l,kind=Rkind))
           d3 = ZERO
         ELSE
           d1 = sqrt(real(2*l,kind=Rkind)) * poly_Hermite(x,l-1)
           d2 = TWO*(x*d1-d0*real(l,kind=Rkind))
           d3 = TWO*(x*d2-d1*real(l-1,kind=Rkind))
         END IF

!        on rajoute la partie exponentielle
!        pour d0 d1 d2, d3
         pexp = exp(-HALF*x*x)
         d3 = ( d3-THREE*x*d2-THREE*(ONE-x*x)*d1+                       &
               THREE*x*(ONE-x*x)*d0 ) * pexp
         d2 = (d2-TWO*x*d1+(x*x-ONE)*d0)*pexp
         d1 = (d1-x*d0)*pexp
         d0 = d0*pexp

       ELSE
         d0 = poly_Hermite(x ,l)*exp(-HALF*x*x)
         d1 = ZERO
         d2 = ZERO
       END IF

!      write(out_unitp,*) 'l x d0 d1 d2 hermite :',l,x,d0,d1,d2


       RETURN
       end subroutine d0d1d2d3poly_Hermite_exp
!===================================================
!
!   calcule pour les polynomes d'hermite 0 a nb_herm
!   la valeur et les dirivees sur les points de la grille
!
!===================================================
      SUBROUTINE d0d1d2poly_Hermite_grille(xh,                          &
        d0h,d1h,d2h,nb_herm,nb_gauss_h,deriv,num,step)
      USE mod_system
      IMPLICIT NONE

       real(kind=Rkind) step
       logical num,deriv

       integer nb_herm,nb_gauss_h

       real(kind=Rkind) xh(nb_gauss_h)
       real(kind=Rkind) d0h(nb_gauss_h,0:nb_herm)
       real(kind=Rkind) d1h(nb_gauss_h,0:nb_herm)
       real(kind=Rkind) d2h(nb_gauss_h,0:nb_herm)

       integer i,k

       DO k=1,nb_gauss_h
         DO i=0,nb_herm
           CALL d0d1d2poly_Hermite(xh(k),i,                             &
                     d0h(k,i),d1h(k,i),d2h(k,i),deriv,num,step)

!          write(out_unitp,*) i,k,d0h(k,i)
!          write(out_unitp,*) i,k,d1h(k,i)
!          write(out_unitp,*) i,k,d2h(k,i)
         END DO
       END DO

       end subroutine d0d1d2poly_Hermite_grille
!===================================================
!
!   calcule la derivee d'un polynome de Hermite (l>0)
!   pour un x ( -inf =< x =< +inf )
!
!    avec la norme sqrt(2**l * l! *sqrt(pi) )
!    Pl'(x) = sqrt(2*l) * Pl-1(x)
!
!===================================================
      SUBROUTINE d0d1d2poly_Hermite_exp(x,l,d0,d1,d2,deriv,num,step)
      USE mod_system
      IMPLICIT NONE

       real(kind=Rkind) x,d0,d1,d2,pexp

       real(kind=Rkind) xp,xm,ep,em,step
       logical num,deriv

      integer l

!----- function -----------------------------------
       real(kind=Rkind) poly_Hermite
!----- function -----------------------------------


!      write(out_unitp,*) 'num',num

       IF (deriv) THEN

         IF (num) THEN
           xp = x+step
           xm = x-step
           d0 = poly_Hermite(x ,l)
           ep = poly_Hermite(xp,l)
           em = poly_Hermite(xm,l)
           CALL d1d2_b(d0,ep,em,d1,d2,step)
         ELSE
           d0 = poly_Hermite(x,l)
           IF (l .EQ. 0) THEN
             d1 = ZERO
             d2 = ZERO
           ELSE IF (l .EQ. 1) THEN
              d1 = sqrt(TWO)*poly_Hermite(x,0)
              d2 = ZERO
           ELSE
             d1 = sqrt(TWO*real(l,kind=Rkind)) * poly_Hermite(x,l-1)
             d2 = TWO*(x*d1-d0*real(l,kind=Rkind))
!            write(out_unitp,*) 'l,x,d0,d1,d2',l,x,d0,d1,d2
           END IF


         END IF

!      write(out_unitp,*) 'l x d0 d1 d2 hermite :',l,x,d0,d1,d2
!        on rajoute la partie exponentielle
!        pour d0 d1 d2
         pexp = exp(-HALF*x*x)
         d2 = (d2-TWO*x*d1+(x*x-ONE)*d0)*pexp
         d1 = (d1-x*d0)*pexp
         d0 = d0*pexp

       ELSE
         d0 = poly_Hermite(x ,l)*exp(-HALF*x*x)
         d1 = ZERO
         d2 = ZERO
       END IF

!      write(out_unitp,*) 'l x d0 d1 d2 hermite :',l,x,d0,d1,d2


       RETURN
       end subroutine d0d1d2poly_Hermite_exp
!===================================================
!
!   calcule pour les polynomes d'hermite 0 a nb_herm
!   la valeur et les dirivees sur les points de la grille
!
!   tout
!   0 : paire
!   1 : impaire
!
!===================================================
      SUBROUTINE d0d1d2poly_Hermite_exp_grille(xh,                      &
        d0h,d1h,d2h,nb_herm,nb_gauss_h,deriv,num,step)
      USE mod_system
      IMPLICIT NONE

       real(kind=Rkind) step
       logical num,deriv

       integer nb_herm,nb_gauss_h

       real(kind=Rkind) xh(nb_gauss_h)
       real(kind=Rkind) d0h(nb_gauss_h,nb_herm)
       real(kind=Rkind) d1h(nb_gauss_h,nb_herm)
       real(kind=Rkind) d2h(nb_gauss_h,nb_herm)

       integer i,k,id

!      write(out_unitp,*) 'num',num
       DO k=1,nb_gauss_h
         id = 0
         DO i=1,nb_herm
           CALL d0d1d2poly_Hermite_exp(xh(k),id,                        &
                     d0h(k,i),d1h(k,i),d2h(k,i),deriv,num,step)

!          write(out_unitp,*) i,k,d0h(k,i)
!          write(out_unitp,*) i,k,d1h(k,i)
!          write(out_unitp,*) i,k,d2h(k,i)
           id = id+1
         END DO
       END DO

       end subroutine d0d1d2poly_Hermite_exp_grille
      SUBROUTINE d0d1d2poly_Hermite_0_exp_grille(xh,                    &
        d0h,d1h,d2h,nb_herm,nb_gauss_h,deriv,num,step)
      USE mod_system
      IMPLICIT NONE

       real(kind=Rkind) step
       logical num,deriv

       integer nb_herm,nb_gauss_h

       real(kind=Rkind) xh(nb_gauss_h)
       real(kind=Rkind) d0h(nb_gauss_h,nb_herm)
       real(kind=Rkind) d1h(nb_gauss_h,nb_herm)
       real(kind=Rkind) d2h(nb_gauss_h,nb_herm)

       integer i,k,id

!      write(out_unitp,*) 'num',num
       DO k=1,nb_gauss_h
         id = 0
         DO i=1,nb_herm
           CALL d0d1d2poly_Hermite_exp(xh(k),id,                        &
                     d0h(k,i),d1h(k,i),d2h(k,i),deriv,num,step)

!          write(out_unitp,*) i,k,d0h(k,i)
!          write(out_unitp,*) i,k,d1h(k,i)
!          write(out_unitp,*) i,k,d2h(k,i)
           id = id+2
         END DO
       END DO

       end subroutine d0d1d2poly_Hermite_0_exp_grille
      SUBROUTINE d0d1d2poly_Hermite_1_exp_grille(xh,                    &
        d0h,d1h,d2h,nb_herm,nb_gauss_h,deriv,num,step)
      USE mod_system
      IMPLICIT NONE

       real(kind=Rkind) step
       logical num,deriv

       integer nb_herm,nb_gauss_h

       real(kind=Rkind) xh(nb_gauss_h)
       real(kind=Rkind) d0h(nb_gauss_h,nb_herm)
       real(kind=Rkind) d1h(nb_gauss_h,nb_herm)
       real(kind=Rkind) d2h(nb_gauss_h,nb_herm)

       integer i,k,id

!      write(out_unitp,*) 'num',num
       DO k=1,nb_gauss_h
         id = 1
         DO i=1,nb_herm
           CALL d0d1d2poly_Hermite_exp(xh(k),id,                        &
                     d0h(k,i),d1h(k,i),d2h(k,i),deriv,num,step)

!          write(out_unitp,*) i,k,d0h(k,i)
!          write(out_unitp,*) i,k,d1h(k,i)
!          write(out_unitp,*) i,k,d2h(k,i)
           id = id+2
         END DO
       END DO

       end subroutine d0d1d2poly_Hermite_1_exp_grille
!===================================================
!
!   calcule la derivee d'un polynome de Hermite (l>0)
!   pour un x ( -inf =< x =< +inf )
!
!    avec la norme sqrt(2**l * l! *sqrt(pi) )
!    Pl'(x) = sqrt(2*l) * Pl-1(x)
!
!===================================================
      SUBROUTINE d0d1d2poly_Hermite_exp_noexp(                          &
                                            x,l,d0,d1,d2,deriv,num,step)
      USE mod_system
      IMPLICIT NONE

       real(kind=Rkind) x,d0,d1,d2,pexp

       real(kind=Rkind) xp,xm,ep,em,step
       logical num,deriv

      integer l

!----- function -----------------------------------
       real(kind=Rkind) poly_Hermite
!----- function -----------------------------------


!      write(out_unitp,*) 'num',num

       IF (deriv) THEN

         IF (num) THEN
           xp = x+step
           xm = x-step
           d0 = poly_Hermite(x ,l)
           ep = poly_Hermite(xp,l)
           em = poly_Hermite(xm,l)
           CALL d1d2_b(d0,ep,em,d1,d2,step)
         ELSE
           d0 = poly_Hermite(x,l)
           IF (l .EQ. 0) THEN
             d1 = ZERO
             d2 = ZERO
           ELSE IF (l .EQ. 1) THEN
              d1 = sqrt(TWO)*poly_Hermite(x,0)
              d2 = ZERO
           ELSE
             d1 = sqrt(TWO*real(l,kind=Rkind)) * poly_Hermite(x,l-1)
             d2 = TWO*(x*d1-d0*real(l,kind=Rkind))
!            write(out_unitp,*) 'l,x,d0,d1,d2',l,x,d0,d1,d2
           END IF


         END IF

!      write(out_unitp,*) 'l x d0 d1 d2 hermite :',l,x,d0,d1,d2
!        On tient compte de exp() pour les derivees mais il n'est pas dans les fonctions
!        pour d0 d1 d2
         d2 = (d2-TWO*x*d1+(x*x-ONE)*d0)
         d1 = (d1-x*d0)
!        d0 = d0

       ELSE
         d0 = poly_Hermite(x ,l)
         d1 = ZERO
         d2 = ZERO
       END IF

!      write(out_unitp,*) 'l x d0 d1 d2 hermite :',l,x,d0,d1,d2


       RETURN
       end subroutine d0d1d2poly_Hermite_exp_noexp
!===================================================
!
!   calcule pour les polynomes d'hermite 0 a nb_herm
!   la valeur et les dirivees sur les points de la grille
!
!===================================================
      SUBROUTINE d0d1d2poly_Hermite_exp_noexp_G(xh,                &
        d0h,d1h,d2h,nb_herm,nb_gauss_h,deriv,num,step)
      USE mod_system
      IMPLICIT NONE

       real(kind=Rkind) step
       logical num,deriv

       integer nb_herm,nb_gauss_h

       real(kind=Rkind) xh(nb_gauss_h)
       real(kind=Rkind) d0h(nb_gauss_h,nb_herm)
       real(kind=Rkind) d1h(nb_gauss_h,nb_herm)
       real(kind=Rkind) d2h(nb_gauss_h,nb_herm)

       integer i,k,id

!      write(out_unitp,*) 'num',num
       DO k=1,nb_gauss_h
         id = 0
         DO i=1,nb_herm
           CALL d0d1d2poly_Hermite_exp_noexp(xh(k),id,                  &
                     d0h(k,i),d1h(k,i),d2h(k,i),deriv,num,step)

!          write(out_unitp,*) i,k,d0h(k,i)
!          write(out_unitp,*) i,k,d1h(k,i)
!          write(out_unitp,*) i,k,d2h(k,i)
           id = id+1
         END DO
       END DO

       end subroutine d0d1d2poly_Hermite_exp_noexp_G


      FUNCTION poly_laguerre(xx,nnn,alpha)
      USE mod_system
      IMPLICIT NONE
      real(kind=Rkind) :: poly_laguerre

       integer, intent (in) :: nnn,alpha

       real(kind=Rkind),intent(in) :: xx

       real(kind=Rkind) :: x,pmm,somx2,fact,pmmp1,pll,poly,norme2

       integer :: n

       n = nnn
       x = xx

       IF (x < ZERO) THEN
         STOP 'x < 0 in poly_laguerre'
       END IF

STOP 'laguerre'

       END FUNCTION poly_laguerre


!=============================================================
!
!      determination des tous les Ln(xi)=serie_fourier(n,i)
!      + les derivees 1er et 2d.
!
!=============================================================
      SUBROUTINE d0d1d2poly_fourier_grille(xf,d0f,d1f,d2f,              &
                              nb_fourier,nb_quadra,deriv,num,step)
      USE mod_system
      IMPLICIT NONE

!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      integer nb_quadra,nb_fourier

      logical num,deriv
      real(kind=Rkind)  step

      real(kind=Rkind)  xf(nb_quadra)
      real(kind=Rkind)  d0f(nb_quadra,nb_fourier)
      real(kind=Rkind)  d1f(nb_quadra,nb_fourier)
      real(kind=Rkind)  d2f(nb_quadra,nb_fourier)

      integer ii,k

      integer is,ic
      real(kind=Rkind)  xiik,c,s
!---------------------------------------------------------------------
      real(kind=Rkind) sq2pi,sqpi
!---------------------------------------------------------------------

      sqpi = ONE/sqrt(pi)
      sq2pi = ONE/sqrt(pi+pi)

!     write(out_unitp,*) 'd0d1d2poly_fourier_grille : deriv',deriv
!     write(out_unitp,*) 'd0d1d2poly_fourier_grille : num',num
!     write(out_unitp,*) 'd0d1d2poly_fourier_grille : step',step

!     pour les 1/2pi   ic=1
      ic=1
      DO k=1,nb_quadra
        d0f(k,ic) = sq2pi
        d1f(k,ic) = ZERO
        d2f(k,ic) = ZERO
      END DO

      ic = 1
      is = nb_fourier/2 + 1
      DO ii=1,nb_fourier/2
        ic = ic + 1
        is = is + 1
!       write(out_unitp,*) 'ii,is,ic',ii,is,ic
        DO k=1,nb_quadra
!         xiik = mod(xf(k)*ii,pi2)
          xiik = xf(k)*real(ii,kind=Rkind)
          s = sin(xiik)*sqpi
          c = cos(xiik)*sqpi

          d0f(k,ic) = c
          d0f(k,is) = s

          d1f(k,ic) = -real(ii,kind=Rkind) * s
          d1f(k,is) =  real(ii,kind=Rkind) * c

          d2f(k,ic) = -real(ii * ii,kind=Rkind) * c
          d2f(k,is) = -real(ii * ii,kind=Rkind) * s
        END DO

      END DO

      end subroutine d0d1d2poly_fourier_grille
!=============================================================
!
!      bessel bj(x) in x for j=0,to nmax
!
!=============================================================
      SUBROUTINE bessel (x, nmax, bj)
      USE mod_system
      IMPLICIT NONE

! Computes Bessel functions Jn(x), from order 0 to order nmax.
! Taken from http://iris-lee3.ece.uiuc.edu/~jjin/routines/mjynb.for

      real(kind=Rkind) x           !  in
      INTEGER nmax                 !

      real(kind=Rkind) bj(0:nmax)  !  out

                                   !  local
      real(kind=Rkind) bs, f2, f1, f, s0, a(4), b(4), a1(4), b1(4),     &
                       t1, p0, q0, cu, bj0, t2, p1, q1, bj1, bjk
      INTEGER k, nm, m
!     real(kind=Rkind) pi
!     DATA pi /3.141592653589793/

      INTEGER msta1, msta2         !  called functions

      bj(:) = ZERO

      IF (x.LT. ONETENTH**100) THEN
         bj(0)= ONE
         RETURN
      ENDIF

      nm= nmax
      IF (x.LE.300_Rkind .OR. nmax.GT.INT(0.9_Rkind*x)) THEN
         IF (nmax.EQ.0) nm= 1
         m= msta1(x,200)
         IF (m.LT.nm) THEN
            nm= m
         ELSE
            m= msta2(x,nm,15)
         ENDIF

         bs= ZERO
         f2= ZERO
         f1= tiny(ONE)
         DO k=m, 0, -1
            f= TWO*real(k+1,kind=Rkind)/x*f1-f2
            IF (k.LE.nm) bj(k)= f
            IF (k.EQ.2*INT(k/2) .AND. k.NE.0) bs= bs+TWO*f
            f2= f1
            f1= f
         ENDDO

         s0= bs+f
         DO k=0, nm
            bj(k)= bj(k)/s0
         ENDDO

      ELSE
         DATA a  /-.0703125000000000_Rkind,  .1121520996093750_Rkind,   &
                  -.5725014209747314_Rkind,  6.074042001273483_Rkind /
         DATA b  / .0732421875000000_Rkind, -.2271080017089844_Rkind,   &
                   1.727727502584457_Rkind, -24.38052969955606_Rkind /
         DATA a1 / .1171875000000000_Rkind, -.1441955566406250_Rkind,   &
                   .6765925884246826_Rkind, -6.883914268109947_Rkind /
         DATA b1 /-.1025390625000000_Rkind,  .2775764465332031_Rkind,   &
                  -1.993531733751297_Rkind,  27.24882731126854_Rkind /

         t1= x-QUARTER*pi
         p0= ONE
         q0= -0.125_Rkind/x
         DO k=1, 4
            p0= p0+a(k)*x**(-2*k)
            q0= q0+b(k)*x**(-2*k-1)
         ENDDO

         cu= sqrt(.63661977236758_Rkind/x)
         bj0= cu*(p0*cos(t1)-q0*sin(t1))
         bj(0)= bj0

         t2= x-0.75_Rkind*pi
         p1= ONE
         q1= 0.375_Rkind/x
         DO k=1, 4
            p1= p1+a1(k)*x**(-2*k)
            q1= q1+b1(k)*x**(-2*k-1)
         ENDDO

         bj1= cu*(p1*cos(t2)-q1*sin(t2))
         bj(1)= bj1

         DO k=2, nm
            bjk= real(2*(k-1),kind=Rkind)/x*bj1-bj0
            bj(k)= bjk
            bj0= bj1
            bj1= bjk
         ENDDO
      ENDIF
      end subroutine bessel

!-----------------------------------------------------------------------

      FUNCTION msta1 (x, mp)
      USE mod_system
      IMPLICIT NONE
      integer :: msta1

! Determine the starting point for backward recurrence such that
! the magnitude of Jn(x) at that point is about 10^(-mp)

      real(kind=Rkind) x              !  in
      INTEGER mp                      !

      real(kind=Rkind) a0, f0, f1, f  !  local
      INTEGER n, n0, n1, it, nn       !

      real(kind=Rkind) envj           !  intrinsic
      envj(n,x)= HALF*log10(6.28_Rkind*real(n,kind=Rkind)) -            &
              real(n,kind=Rkind)*log10(1.36_Rkind*x/real(n,kind=Rkind))

      a0= abs(x)
      n0= INT(1.1_Rkind * a0)+1
      f0= envj(n0,a0)-real(mp,kind=Rkind)
      n1= n0+5
      f1= envj(n1,a0)-real(mp,kind=Rkind)

      DO it=1, 20
         nn= n1-int(real(n1-n0,kind=Rkind)/(ONE-f0/f1))
         f= envj(nn,a0)-real(mp,kind=Rkind)
         IF (ABS(nn-n1).LT.1) EXIT
         n0= n1
         f0= f1
         n1= nn
         f1= f
      ENDDO

      msta1= nn

      end function msta1

!-----------------------------------------------------------------------

      FUNCTION msta2 (x, n, mp)
      USE mod_system
      IMPLICIT NONE
      INTEGER :: msta2

! Determine the starting point for backward recurrence such that
! all Jn(x) has mp significant digits

      real(kind=Rkind) x                             !  in
      INTEGER n, mp                                  !

      real(kind=Rkind) a0, hmp, ejn, obj, f0, f1, f  !  local
      INTEGER n0, n1, it, nn                         !

      real(kind=Rkind) envj                          !  intrinsic
      envj(n,x)= HALF*log10(6.28_Rkind*real(n,kind=Rkind)) -            &
              real(n,kind=Rkind)*log10(1.36_Rkind*x/real(n,kind=Rkind))

      a0= abs(x)
      hmp= HALF*real(mp,kind=Rkind)
      ejn= envj(n,a0)
      IF (ejn.LE.hmp) THEN
         obj= real(mp,kind=Rkind)
         n0= INT(1.1_Rkind*a0)
      ELSE
         obj= hmp+ejn
         n0= n
      ENDIF

      f0= envj(n0,a0)-obj
      n1= n0+5
      f1= envj(n1,a0)-obj
      DO it=1, 20
         nn= n1-int(real(n1-n0,kind=Rkind)/(ONE-f0/f1))
         f= envj(nn,a0)-obj
         IF (ABS(nn-n1).LT.1) EXIT
         n0= n1
         f0= f1
         n1= nn
         f1= f
      ENDDO

      msta2= nn+10

      end function msta2
!###############################################################################
!   subroutine / function  mmbsjn
!
!   imsl routine name   - mmbsjn
!
!-----------------------------------------------------------------------
!
!   computer            - cray/single
!
!   latest revision     - november 1, 1984
!
!   purpose             - bessel function of the first kind of
!                           nonnegative integer order for
!                           real arguments
!
!   usage               - call mmbsjn (arg,n,b,ier)
!
!   arguments    arg    - input argument. the absolute value of arg must
!                           be less than or equal to 100000. arg must be
!                           typed appropriately in the calling program.
!                           (see the precision/hardware section.)
!                n      - input parameter specifying the number of
!                           function values to be computed.
!                b      - output vector of length n containing the
!                           computed function values. b must be typed
!                           appropriately in the calling program.
!                           (see the precision/hardware section.)
!                           b(1) will contain the computed value for
!                           order zero, b(2) will contain the computed
!                           value for order 1, b(3) for order 2, etc.
!                ier    - error parameter. (output)
!                         terminal error
!                           ier = 129 indicates that either arg or n is
!                             out of range. b(i), (i=1,n) is set to
!                             machine infinity.
!                           ier = 129 + j indicates that b(i), (i=1,j)
!                             are computed to machine precision, but
!                             precision is lost for b(i), (i=j+1,n.)
!                             see the programming notes.
!
!   precision/hardware  - double/h32,h36
!                       - single/h48,h60
!
!   reqd. imsl routines - uertst,ugetio
!
!   notation            - information on special notation and
!                           conventions is available in the manual
!                           introduction or through imsl routine uhelp
!
!   copyright           - 1982 by imsl, inc. all rights reserved.
!
!   warranty            - imsl warrants only that imsl testing has been
!                           applied to this code. no other warranty,
!                           expressed or implied, is applicable.
!
!-----------------------------------------------------------------------
!
      SUBROUTINE mmbsjn (arg,n,b,ier)
      USE mod_system
      IMPLICIT NONE
!                                  specifications for arguments
      integer            n,ier
      real(kind=Rkind)             arg,b(*)
!                                  specifications for local variables
      integer            l,largex,magx,ncalc,nn,nbmx,m,nstart,nend
      integer            itemp
      real(kind=Rkind)  test,tempa,tempb,tempc,p
      real(kind=Rkind)  rsign,sum,tover,plast,pold,psave,psavel
      real(kind=Rkind)  dsig,rten,tmpa4,smallx,xinf
      data               dsig/SEVEN/
      data               rten/44._Rkind/
      data               xinf/1.7e+38_Rkind/
      data               largex /100000/
!                                  first executable statement
      ier = 0
      tempa = abs(arg)
      magx =  int(tempa)
      if(n.gt.0 .and. magx.le.largex) go to 10
!                                  error return -- arg,n is out of range
      ier = 129
      b(1) = xinf
      if(n.lt.2) go to 9000
      do 5 l=2,n
         b(l) = xinf
    5 continue
      go to 9000
   10 rsign = ONE
      ncalc = n
!                                  use 2-term ascending series for
!                                    small arg
      tmpa4 = tempa**FOUR
      smallx = ONETENTH**dsig
      if(tmpa4.ge.smallx) go to 20
!                                  two-term ascending series for
!                                    small arg
      tempa = ONE
      tempb = -FOURTH*arg*arg*rsign
      b(1) = ONE+tempb
      if(n.eq.1) go to 9005
      do 15 nn=2,n
         tempa = tempa*arg/real(2*nn-2,kind=Rkind)
         b(nn) = tempa*(ONE+tempb/real(nn,kind=Rkind))
   15 continue
      go to 9005
!                                  initialize the calculation of p*s
   20 nbmx = n-magx
      nn = magx+1
      plast = ONE
      p = real(2*nn,kind=Rkind)/tempa
!                                  calculate general significance test
      test = TWO * TEN**dsig
      m = 0
      if(nbmx.lt.3) go to 30
!                                  calculate p*s until nn=n-1.
!                                    check for possible overflow.
      tover = TEN**(rten-dsig)
      nstart = magx+2
      nend = n-1
      do 25 nn=nstart,nend
         pold = plast
         plast = p
         p = real(2*nn,kind=Rkind)*plast/tempa-rsign*pold
         if(p-tover) 25, 25, 35
   25 continue
      nn = nend
!                                  calculate special significance test
!                                    for nbmx.gt.2.
!
      test = max(test,sqrt(plast*TEN**dsig)*sqrt(TWO*p))
!
!                                  calculate p*s until significance
!                                    test passes
   30 nn = nn+1
      pold = plast
      plast = p
      p = real(2*nn,kind=Rkind)*plast/tempa-rsign*pold
      if(p.lt.test) go to 30
      if(m.eq.1) go to 55
!                                  for j*s, a strong variant of the test
!                                    is necessary. calculate it, and
!                                    calculate p*s until this test is
!                                    passed.
      m = 1
      tempb = p/plast
      tempc = real(nn+1,kind=Rkind)/tempa
      if(tempb+ONE/tempb .gt. TWO*tempc) tempb=tempc+sqrt(tempc**2-ONE)
      test = test/sqrt(tempb-ONE/tempb)
      if(p-test) 30, 55, 55
!                                  to avoid overflow, divide p*s by
!                                    tover.  calculate p*s until
!                                    abs(p).gt.1.
   35 tover = TEN**rten
      p = p/tover
      plast = plast/tover
      psave = p
      psavel = plast
      nstart = nn+1
   40 nn = nn+1
      pold = plast
      plast = p
      p = real(2*nn,kind=Rkind)*plast/tempa-rsign*pold
      if(p .le. ONE) go to 40
      tempb = real(2*nn,kind=Rkind)/tempa
      tempc = HALF*tempb
      tempb = plast/pold
      if(tempb+ONE/tempb .gt. TWO*tempc) tempb=tempc+sqrt(tempc**2-ONE)
!
!                                  calculate backward test, and find
!                                    ncalc, the highest nn such that the
!                                    test is passed.
      test = HALF*pold*plast*(ONE-ONE/tempb**2)/TEN**dsig
      p = plast*tover
      nn = nn-1
      nend = min(n,nn)
      do 45 ncalc=nstart,nend
         pold = psavel
         psavel = psave
         psave = real(2*nn,kind=Rkind)*psavel/tempa-rsign*pold
         if(psave*psavel-test) 45, 45, 50
   45 continue
      ncalc = nend+1
   50 ncalc = ncalc-1
!                                  the sum b(1)+2b(3)+2b(5)... is used
!                                    to normalize. m, the coefficient of
!                                    b(nn), is initialized to 2 or 0.
   55 nn = nn+1
      m = 2*nn-4*(nn/2)
!                                  initialize the backward recursion and
!                                    the normalization sum
      tempb = ZERO
      tempa = ONE/p
      sum = real(m,kind=Rkind)*tempa
      nend = nn-n
      if(nend) 80, 70, 60
!                                  recur backward via difference
!                                    equation, calculating (but not
!                                    storing) b(nn), until nn=n.
   60 do 65 l=1,nend
         nn = nn-1
         tempc = tempb
         tempb = tempa
         tempa = real(2*nn,kind=Rkind)*tempb/arg-rsign*tempc
         m = 2-m
         sum = sum+real(m,kind=Rkind)*tempa
   65 continue
!                                  store b(nn)
   70 b(nn) = tempa
      if(n.gt.1) go to 75
!                                  n=1.  since 2*tempa is added to the
!                                    sum, tempa must be subtracted
      sum = sum-tempa
      go to 110
!                                  calculate and store b(nn-1)
   75 nn = nn-1
      b(nn) = real(2*nn,kind=Rkind)*tempa/arg-rsign*tempb
      if(nn.eq.1) go to 105
      m = 2-m
      sum = sum+real(m,kind=Rkind)*b(nn)
      go to 90
!                                  nn.lt.n, so store b(nn) and set
!                                  higher orders to zero
   80 b(nn) = tempa
      nend = -nend
      do 85 l=1,nend
         itemp = nn+l
         b(itemp) = ZERO
   85 continue
   90 nend = nn-2
      if(nend.eq.0) go to 100
!                                  calculate via difference equation and
!                                    store b(nn), until nn=2
      do 95 l=1,nend
         nn = nn-1
         b(nn) = real(2*nn,kind=Rkind)*b(nn+1)/arg-rsign*b(nn+2)
         m = 2-m
         sum = sum+real(m,kind=Rkind)*b(nn)
   95 continue
!                                  calculate b(1)
  100 b(1) = TWO*b(2)/arg-rsign*b(3)
  105 sum = sum+b(1)
!                                  normalize--if ize=1, divide sum by
!                                    cosh(arg). divide all b(nn) by sum.
  110 continue
      do 115 nn=1,n
  115 b(nn) = b(nn)/sum
      if(ncalc.eq.n) go to 9005
      ier = 129+ncalc
 9000 continue
!      call uertst(ier,'mmbsjn')
 9005 return
      end subroutine mmbsjn
!================================================================
!    fonction d0ylm
!================================================================

      SUBROUTINE d0Ylm(d0,x,i,ndim)
      USE mod_system
      IMPLICIT NONE

       real(kind=Rkind) d0,x(2)
       integer i,ndim

       integer l,m,lll,mmm
       real(kind=Rkind) th,phi

!      ---------------------------------------------------------------
       real(kind=Rkind) Ylm,poly_legendre,fourier
!      ---------------------------------------------------------------

       IF (ndim .NE. 2) THEN
         write(out_unitp,*) ' ERROR in d0Ylm'
         write(out_unitp,*) ' ndim MUST set to 2 (ndim=',ndim,')'
         STOP
       END IF

       th=x(1)
       phi=x(2)


       l = int(sqrt(real(i,kind=Rkind)-HALF))
       lll = l+1
       mmm = i-l*l
       m = mmm/2



       !write(out_unitp,*) 'ylm',i,lll,mmm

       d0 = poly_legendre(cos(th),lll,m)*fourier(phi,mmm)


       end subroutine d0Ylm
!================================================================
!    fonction d0ylm
!================================================================

      SUBROUTINE d0d1d2Ylm(d0,d1,d2,x,i,num,step)
      USE mod_system
      IMPLICIT NONE

       real(kind=Rkind) d0,d1(2),d2(2,2),x(2)
       integer i
       logical num
       real(kind=Rkind) step

       integer l,m,lll,mmm
       real(kind=Rkind) th,phi,cth,sth
       real(kind=Rkind) d0fm,d1fm,d2fm,d3fm
       real(kind=Rkind) d0plm,d1plm,d2plm



!      ---------------------------------------------------------------
       real(kind=Rkind) poly_legendre
!      ---------------------------------------------------------------

!      write(out_unitp,*) i,num,step,x,d0,d1,d2

       th=x(1)
       phi=x(2)
       cth = cos(th)
       sth = sin(th)


       l = int(sqrt(real(i,kind=Rkind)-HALF))
       lll = l+1
       mmm = i-l*l
       m = mmm/2




       CALL d0d1d2d3fourier(phi,d0fm,d1fm,d2fm,d3fm,mmm)

       d0plm  = poly_legendre(cth,lll,m)
       IF (num) THEN
         d1plm = poly_legendre(cos(th+step),lll,m) -                    &
                  poly_legendre(cos(th-step),lll,m)
         d1plm = d1plm /(TWO*step)
         d2plm = poly_legendre(cos(th+step),lll,m) +                    &
                  poly_legendre(cos(th-step),lll,m)
         d2plm = (d2plm-d0plm-d0plm)/(step*step)
       ELSE
         d1plm  = real(m,kind=Rkind)*cth/sth*d0plm +                    &
               sqrt(real(l+m+1,kind=Rkind)*real(l-m,kind=Rkind)) *      &
                 poly_legendre(cth,lll,m+1)
         d2plm  = (real(m*m,kind=Rkind)/(sth*sth) -                     &
                   real(l*l+l,kind=Rkind))*d0plm -                      &
                  cth/sth*d1plm
       END IF



!      write(out_unitp,*) i,l,m,d0plm,d1plm,d2plm


       d0      = d0plm * d0fm
       d1(1)   = d1plm * d0fm
       d1(2)   = d0plm * d1fm
       d2(1,1) = d2plm * d0fm
       d2(1,2) = d1plm * d1fm
       d2(2,1) = d2(1,2)
       d2(2,2) = d0plm * d2fm

     END SUBROUTINE d0d1d2Ylm
     SUBROUTINE d0d1d2d3Ylm(d0,d1,d2,d3,x,i)
      USE mod_system
      IMPLICIT NONE

       real(kind=Rkind) d0,d1(2),d2(2,2),d3(2,2,2),x(2)
       integer i
       logical num
       real(kind=Rkind) step

       integer l,m,lll,mmm
       real(kind=Rkind) th,phi,cth,sth
       real(kind=Rkind) d0fm,d1fm,d2fm,d3fm
       real(kind=Rkind) d0plm,d1plm,d2plm,d3plm



       th  = x(1)
       phi = x(2)
       cth = cos(th)
       sth = sin(th)


       l = int(sqrt(real(i,kind=Rkind)-HALF))
       lll = l+1
       mmm = i-l*l
       m = mmm/2




       CALL d0d1d2d3fourier(phi,d0fm,d1fm,d2fm,d3fm,mmm)
       CALL d0d1d2d3Plm(d0plm,d1plm,d2plm,d3plm,th,l,m)

!      write(out_unitp,*) i,l,m,d0plm,d1plm,d2plm,d3plm


       d0        = d0plm * d0fm

       d1(1)     = d1plm * d0fm
       d1(2)     = d0plm * d1fm

       d2(1,1)   = d2plm * d0fm
       d2(1,2)   = d1plm * d1fm
       d2(2,1)   = d2(1,2)
       d2(2,2)   = d0plm * d2fm

       d3(1,1,1) = d3plm * d0fm
       d3(2,1,1) = d2plm * d1fm
       d3(1,2,1) = d3(2,1,1)
       d3(1,1,2) = d3(2,1,1)
       d3(1,2,2) = d1plm * d2fm
       d3(2,1,2) = d3(1,2,2)
       d3(2,2,1) = d3(1,2,2)
       d3(2,2,2) = d0plm * d3fm

      END SUBROUTINE d0d1d2d3Ylm
!===================================================
!
!   Ylm(the,phi)
!    m >= 0  plm(th,l, m)*cos( m phi) *norme
!    m <  0  plm(th,l,-m)*sin(-m phi) *norme
!
!    l=lll m = mmm/2
!    mmm dans l'ordre de la serie de fourier v26
!===================================================
      FUNCTION Ylm(th,phi,lll,mmm)
      USE mod_system
      IMPLICIT NONE
      real(kind=Rkind) :: Ylm

       real(kind=Rkind) th,phi
       integer lll,mmm
       integer l,m

!     ----------------------------------------------------------------
      real(kind=Rkind) sq2pi,sqpi
!     ----------------------------------------------------------------

!     - function -----------------------------------------------------
      real(kind=Rkind) poly_legendre,v26,v27,v28,fourier
!     ----------------------------------------------------------------


       l=lll
       m = mmm/2

       IF ( m > l .OR. l < 0 .OR. th >= pi .OR. th <= ZERO) THEN
         write(out_unitp,*) 'mauvais arguments dans Ylm :'
         write(out_unitp,*) ' m l : ',m,l,' et th = ',th
         write(out_unitp,*) ' mmm lll : ',mmm,lll
         STOP
       END IF


       Ylm = poly_legendre(cos(th),lll,m)*fourier(phi,mmm)


       end function Ylm

      SUBROUTINE d0d1d2Plm(d0plm,d1plm,d2plm,th,l,m)
      USE mod_system
      IMPLICIT NONE

       real(kind=Rkind) :: th,d0plm,d1plm,d2plm
       integer          :: l,m


       integer          :: lll,mmm
       real(kind=Rkind) :: cth,sth



!      ---------------------------------------------------------------
       real(kind=Rkind) :: poly_legendre
!      ---------------------------------------------------------------

       cth = cos(th)
       sth = sin(th)
       lll = l+1
       mmm = abs(m)


       d0plm  = poly_legendre(cth,lll,mmm)
       d1plm  = real(mmm,kind=Rkind)*cth/sth*d0plm +                    &
                sqrt(real(l+mmm+1,kind=Rkind)*real(l-mmm,kind=Rkind)) * &
                poly_legendre(cth,lll,mmm+1)
       d2plm  = (real(mmm*mmm,kind=Rkind)/(sth*sth) -                   &
                 real(l*l+l,kind=Rkind))*d0plm - cth/sth*d1plm


!      write(out_unitp,*) lll,m,d0plm,d1plm,d2plm

      END SUBROUTINE d0d1d2Plm
      SUBROUTINE d0d1d2d3Plm(d0plm,d1plm,d2plm,d3plm,th,l,m)
      USE mod_system
      IMPLICIT NONE

       real(kind=Rkind) :: th,d0plm,d1plm,d2plm,d3plm
       integer          :: l,m


       integer          :: lll,mmm
       real(kind=Rkind) :: cth,sth,Rm2,Rllp1



!      ---------------------------------------------------------------
       real(kind=Rkind) :: poly_legendre
!      ---------------------------------------------------------------

       cth = cos(th)
       sth = sin(th)
       lll = l+1
       mmm = abs(m)
       Rm2 = real(mmm*mmm,kind=Rkind)
       Rllp1 = real(l*l+l,kind=Rkind)



       d0plm  = poly_legendre(cth,lll,mmm)

       d1plm  = real(mmm,kind=Rkind)*cth/sth*d0plm +                    &
                sqrt(real(l+mmm+1,kind=Rkind)*real(l-mmm,kind=Rkind)) * &
                poly_legendre(cth,lll,mmm+1)


       d2plm  = (Rm2/sth**2 - Rllp1) * d0plm - cth/sth * d1plm


       d3plm  = Rm2*(-TWO*cth/sth**3)  * d0plm +                        &
               (Rm2/sth**2 - Rllp1) * d1plm - &
               (-ONE/sth**2) * d1plm - cth/sth * d2plm

!      write(out_unitp,*) lll,m,d0plm,d1plm,d2plm

       END SUBROUTINE d0d1d2d3Plm
