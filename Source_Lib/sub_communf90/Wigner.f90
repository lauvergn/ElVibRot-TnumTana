MODULE wigner_m
  USE mod_system
  Implicit none

  PRIVATE
  PUBLIC :: Wigner
  CONTAINS

  function lnfac(i)
    real(kind=Rkind)            :: lnfac
    integer,        intent(in)  :: i

    IF (i == 0 .OR. i == 1) THEN 
      lnfac = ZERO
    ELSE IF (mod(i,2) == 0) THEN 
      lnfac = log(real(i,kind=Rkind))
    ELSE
     lnfac = log(real(i,kind=Rkind)) + log(real(i-1,kind=Rkind))
    END IF

  end function lnfac

  ! complex(kind=Rkind) :: Djmk,Wigner
      ! Integer :: j,m,k
      ! real(kind=Rkind) :: phi,theta,chi

      ! j=3
      ! k=2
      ! m=1
      ! phi=0d0
      ! theta=1d0
      ! chi=2d00
      ! Djmk=Wigner(j,m,k,phi,theta,chi)
      ! Write (6,*) 'Djmk=',Djmk
      ! End
!     ===+=========+=========+=========+=========+=========+=========+==
      Function Wigner(j,m,k,phi,theta,chi) Result(Djmk)
!     ===+=========+=========+=========+=========+=========+=========+==
      Implicit none
      complex(kind=Rkind)             :: Djmk
      Integer,             intent(in) :: j,k,m
      real(kind=Rkind),    intent(in) :: phi,theta,chi

      !Local variables
      complex(kind=Rkind) :: cmplx
      Integer             :: iabs,itemp,keff,weff
      real(kind=Rkind)    :: dJ_MK,cos,sin,norm,cosmFi,coskQi,sinmFi,sinkQi
      Logical             :: first
      Integer, parameter  :: nbFac=100
      real(kind=Rkind)    :: fact(0:nbFac)

      first = .true.
!     norm=sqrt(j+half)
      norm = ONE
      keff = k
      weff = m

      if(weff+keff < 0) Then
         norm = norm*(-1)**abs(m-k)
         weff = -weff
         keff = -keff
      Endif

      if(weff-keff < 0) Then
         norm = norm*(-1)**abs(m-k)
         itemp = weff
         weff = keff
         keff = itemp
      Endif

      Djmk = norm * dJ_MK(j,weff,keff,cos(theta),fact,first) *                  &
                   cmplx(cos(m*phi),-sin(m*phi)) * cmplx(cos(k*chi),-sin(k*chi))

      End Function Wigner
!     ===+=========+=========+=========+=========+=========+=========+==
      Function dJ_MK(J,M,K,u,fact,first)
!     ===+=========+=========+=========+=========+=========+=========+==
!     Nikiforov et al.,Eq.(5.1.26)
!     ===+=========+=========+=========+=========+=========+=========+==
      real(kind=Rkind)                   :: dJ_MK
      Integer,          INTENT(IN)       :: J,M,K
      real(kind=Rkind), INTENT(IN)       :: u

      Integer, parameter  :: nbFac=100
      real(kind=Rkind), INTENT(INOUT)    :: fact(0:nbFac)
      Logical,          INTENT(INOUT)    :: first

      real(kind=Rkind)    :: Jacobi ! function
      real(kind=Rkind)    :: alf,bet,sign

      if(first) then
         if(2*J.gt.nbFac) stop 'Increase nbFac in <dJ_MK>'
         fact(0)=zero
         fact(1)=zero
         do i=2,2*J
            fact(i)=alog(float(i))
         enddo
         do i=3,2*J
            fact(i)=fact(i)+fact(i-1)
         enddo
         first=.false.
      endif

      alf=M-K
      bet=M+K
      sign=ONE
      if(mod(int(alf),2).ne.0) sign=-ONE

      dJ_MK=sign/TWO**M                                                         &
           *exp(HALF*(fact(J+M)+fact(J-M)-fact(J+K)-fact(J-K)))                 &
           *(1-u)**(alf/2)*(1+u)**(bet/2)*Jacobi(J-M,alf,bet,u)

      End Function dJ_MK
!     ===+=========+=========+=========+=========+=========+=========+==
      Function Jacobi(J,alf,bet,x)
!     ================================================================
!     Abramowitz & Stegun , p782 + (22.4.1) for n=1 value
!     ===+=========+=========+=========+=========+=========+=========+==

      real(kind=Rkind)             ::  Jacobi
      real(kind=Rkind), INTENT(IN) :: alf,bet,x
      integer,          INTENT(IN) :: J


      real(kind=Rkind) ::  Jnp1,Jn,Jn_1,a1n,a2n,a3n,a4n,b3
      integer          :: n

      if(J.eq.0) Then
         Jnp1 = one
      Else
         Jnp1 = HALF*(alf-bet+(alf+bet+2)*x)
         Jn   = one
         Jn_1 = zero
      Endif
      do n=1,J-1
         Jn_1 = Jn
         Jn   = Jnp1
         a1n  = 2*(n+1)*(n+alf+bet+1)*(2*n+alf+bet)
         a2n  = (2*n+alf+bet+1)*(alf**2-bet**2)
         b3   = (2*n+alf+bet)
         a3n  = b3*(b3+1)*(b3+2)
         a4n  = 2*(n+alf)*(n+bet)*(2*n+alf+bet+2)
         Jnp1 = (a2n+x*a3n)*Jn-a4n*Jn_1
         Jnp1 = Jnp1/a1n
      enddo

      Jacobi = Jnp1

    End Function Jacobi
!     ===+=========+=========+=========+=========+=========+=========+==
      Function F3J(j1,j2,j3,m1,m2,m3,fact,nbFac,init)
!     ===+=========+=========+=========+=========+=========+=========+==
!     Calculates 3j coefficients from Racah formula
!     (Messiah: vol.2, p 910; formula 21) .
!     ------------------------------------------------------------------
!     has been tested for j up to 200.
!     j.m.l. (1975)
!     ===+=========+=========+=========+=========+=========+=========+==
      Implicit real(kind=Rkind)(a-h,o-z)
      real(kind=Rkind)                   :: f3j
      Integer,          INTENT(IN)       :: j1,j2,j3,m1,m2,m3
      Integer,          INTENT(IN)       :: nbFac
      real(kind=Rkind), INTENT(INOUT)    :: fact(0:nbFac)
      Logical,          INTENT(INOUT)    :: init

      real(kind=Rkind)  :: phase
      Integer           :: t,tmin,tmax
      Integer           :: k1,k2,k3,k4,k5
      Integer           :: n1,n2,n3,n4,n5,n6,n7,n8,n9,n10,Nmax

      if(init) Then
         fact(0)=0d0
         fact(1)=0d0
         do 20 i=2,nbFac
   20    fact(i)=log(float(i))
         do 30 i=3,nbFac
   30    fact(i)=fact(i)+fact(i-1)
         init=.false.
      Endif

      f3j = ZERO

      if (j3.gt.j1+j2) goto 100
      if (abs(j1-j2).gt.j3) goto 100
      if (abs(m1+m2+m3).gt.0) goto 100
      if (abs(m1).gt.j1) goto 100
      if (abs(m2).gt.j2) goto 100
      if (abs(m3).gt.j3) goto 100

      k1=j3-j2+m1
      k2=j3-j1-m2
      k3=j1-m1
      k4=j2+m2
      k5=j1+j2-j3
      tmin = 0
      if (k1+tmin .lt. 0) tmin = -k1
      if (k2+tmin .lt. 0) tmin = -k2
      tmax = k3
      if (k4-tmax .lt. 0) tmax = k4
      if (k5-tmax .lt. 0) tmax = k5
      n1 = j1+j2-j3
      n2 = j2+j3-j1
      n3 = j3+j1-j2
      n4 = j1+m1
      n5 = j2+m2
      n6 = j3+m3
      n7 = j1-m1
      n8 = j2-m2
      n9 = j3-m3
      n10 = j1+j2+j3+1
      Nmax=Max(n1,n2,n3,n4,n5,n6,n7,n8,n9,n10)
      Nmax=Max(Nmax,Tmax,Tmax+k1,Tmax+2,k3-Tmin,k4-Tmin,k5-Tmin)

      if(Nmax.gt.nbFac) then
          Write (6,1000) Nmax
 1000     format[/' $$$$$ In F3j Nmax(',i3,') > nbFac $$$$$')
          stop
      Endif

      x = fact(n1)+fact(n2)+fact(n3)+fact(n4)+fact(n5)+fact(n6)               &
         +fact(n7)+fact(n8)+fact(n9)-fact(n10)
      x = 0.5*x
      do  t = tmin,tmax
         phase = one
         if (mod(t,2) .ne. 0) phase = -one
         f3j=f3j+phase*exp(-fact(t) - fact(k1+t) - fact(k2+t)                 &
                           -fact(k3-t) - fact(k4-t) - fact(k5-t)+x)
      end do

      if (mod(abs(j1-j2-m3),2) .gt. 0) f3j=-f3j

      100  return
      End Function F3J
END MODULE wigner_m
