!================================================================
!    fonction v_inter(x)
!    v1  2D Ylm
!    v2  2D Ylm m>=0 : cos(m phi)
!    v3  2D Ylm m<0 : sin(m phi) not yet
!
!    v5  poly nD
!
!    v6  1D Huffaker
!    v7  1D BO
!    v8  1D (x-Req)/(c(1)*x+c(2)*Req)
!     c(1)=0 c(2)=1 : Dunham
!     c(1)=1 c(2)=0 : SPF
!     c(1)=0 c(2)=1 : Ogilvie
!    v9  1D x**i*(x-1)**k
!
!    v12 1D poly legendre          en cos(x)
!    v13 1D poly legendre paire    en cos(x)
!    v14 1D poly legendre impaire  en cos(x)
!
!    v15 1D x^n
!    v16 1D x^n paire
!    v17 1D x^n impaire
!    v18 1D x^n paire sans n=0
!    v19 1D x^n sans n=0
!
!    v22 1D poly legendre          en x
!    v23 1D poly legendre paire    en x
!    v24 1D poly legendre impaire  en x
!
!    v26 1D serie de fourier en x
!    v27 1D serie de fourier paire en x
!    v28 1D serie de fourier impaire en x
!    v29 1D serie de fourier en 2x
!    v30 1D serie de fourier en 3x
!
!
!    v32 1D poly Hermite_exp          en x
!    v33 1D poly Hermite_exp paire    en x
!    v34 1D poly Hermite_exp impaire  en x
!
!    Fourier suite:
!    v41 1D serie de fourier en cos(nx)
!    v42 1D serie de fourier en cos(2nx)
!    v43 1D serie de fourier en cos(3nx)
!
!    v51 1D serie de fourier en sin(nx)
!    v52 1D serie de fourier en sin(2nx)
!    v53 1D serie de fourier en sin(3nx)
!================================================================

       FUNCTION v(x,ndim,kl,n,ntyp)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v ! function
       integer ntyp,ndim

       real (kind=Rkind) :: v1,v2,v3,v4,v_DirProd,v_poly,v6,v7,v8,v9
       real (kind=Rkind) :: v10,v11,v12,v13,v14,v15,v16,v17,v18,v19
       real (kind=Rkind) :: v20,v21,v22,v23,v24,v25,v26,v27,v28,v29
       real (kind=Rkind) :: v30,v31,v32,v33,v34,v35,v36,v37,v38,v39
       real (kind=Rkind) :: v40,v41,v42,v43,v44,v45,v46,v47,v48,v49
       real (kind=Rkind) :: v50,v51,v52,v53,v54,v55,v56,v57,v58,v59
       real (kind=Rkind) :: x(ndim)
       integer kl,n(0:ndim)


        SELECT CASE (ntyp)
        CASE (1)
          v=v1(x,ndim,kl,n)
        CASE (2)
          v=v2(x,ndim,kl,n)
        CASE (3)
          v=v3(x,ndim,kl,n)
        CASE (4)
          !v=v4(x,ndim,kl,n)
          v=v_DirProd(x,ndim,kl,n)
        CASE (5)
          !v=v5(x,ndim,kl,n)
          v=v_poly(x,ndim,kl,n)
        CASE (6)
          v=v6(x,ndim,kl,n)
        CASE (7)
          v=v7(x,ndim,kl,n)
        CASE (8)
          v=v8(x,ndim,kl,n)
        CASE (9)
          v=v9(x,ndim,kl,n)
        CASE (10)
          v=v10(x,ndim,kl,n)
        CASE (11)
          v=v11(x,ndim,kl,n)
        CASE (12)
          v=v12(x,ndim,kl,n)
        CASE (13)
          v=v13(x,ndim,kl,n)
        CASE (14)
          v=v14(x,ndim,kl,n)
        CASE (15)
          v=v15(x,ndim,kl,n)
        CASE (16)
          v=v16(x,ndim,kl,n)
        CASE (17)
          v=v17(x,ndim,kl,n)
        CASE (18)
          v=v18(x,ndim,kl,n)
        CASE (19)
          v=v19(x,ndim,kl,n)

        CASE (22)
          v=v22(x,ndim,kl,n)
        CASE (23)
          v=v23(x,ndim,kl,n)
        CASE (24)
          v=v24(x,ndim,kl,n)

        CASE (26)
          v=v26(x,ndim,kl,n)
        CASE (27)
          v=v27(x,ndim,kl,n)
        CASE (28)
          v=v28(x,ndim,kl,n)
        CASE (29)
          v=v29(x,ndim,kl,n)
        CASE (30)
          v=v30(x,ndim,kl,n)

        CASE (41)
          v=v41(x,ndim,kl,n)
        CASE (42)
          v=v42(x,ndim,kl,n)
        CASE (43)
          v=v43(x,ndim,kl,n)
        CASE (51)
          v=v51(x,ndim,kl,n)
        CASE (52)
          v=v52(x,ndim,kl,n)
        CASE (53)
          v=v53(x,ndim,kl,n)

        CASE (32)
          v=v32(x,ndim,kl,n)
        CASE (33)
          v=v33(x,ndim,kl,n)
        CASE (34)
          v=v34(x,ndim,kl,n)

        CASE default ! ERROR: wrong function !
          write(out_unitp,*) ' ERROR wrong function, ntyp',ntyp
          STOP
        END SELECT

       end function v

!================================================================
!    fonction v1(x,i)  ylm
!================================================================

       FUNCTION v1(x,ndim,i,n)
       USE mod_system
       IMPLICIT NONE
       real (kind=Rkind) :: v1 ! function

       integer ndim
       real (kind=Rkind) :: x(ndim)
       integer n(0:ndim)

       integer i,ii,l,m,lll,mmm
       real (kind=Rkind) :: th,phi

!      ---------------------------------------------------------------
       real (kind=Rkind) :: Ylm
!      ---------------------------------------------------------------

       th=x(1)
       phi=x(2)

       l = int(sqrt(real(i,kind=Rkind)-HALF))
       lll = l+1
       mmm = i-l*l


!      write(out_unitp,*) 'ylm',i,lll,mmm

       v1 = Ylm(th,phi,lll,mmm)



       end function v1
!================================================================
!    fonction v2(x,i)  ylm m>=0
!================================================================

       FUNCTION v2(x,ndim,i,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v2 ! function

       integer ndim
       real (kind=Rkind) :: x(ndim)
       integer n(0:ndim)

       integer i,ii,l,m,lll,mmm
       real (kind=Rkind) :: xi,xl,th,phi(1)

!      ---------------------------------------------------------------
       real (kind=Rkind) :: Ylm,poly_legendre,v27
!      ---------------------------------------------------------------

       th=x(1)
       phi=x(2)


       xi = real(i,kind=Rkind)-HALF

       xl = (sqrt(EIGHT*xi+ONE)-THREE)*HALF

       l = int(xl)
       IF (xl .LT. ZERO) l=-1
       lll = l+2

       mmm = i-(l+1)*(l+2)/2
       IF (l .EQ. -1) mmm=1
       m = mmm-1


!      write(out_unitp,*) 'ylm i,lll,mmm',i,lll,mmm
!      write(out_unitp,*) 'ylm i,l,m',i,l,m

       v2 = poly_legendre(cos(th),lll,m)*v27(phi,1,mmm,n)



       end function v2
!================================================================
!    fonction v2(x,ndim,ij,n)
!    x = cos(x(1) et R=x(2)
!    Pi0(x) * R^j * exp(-R)
!================================================================

       FUNCTION v2_old(x,ndim,ij,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v2_old ! function

       integer ndim
       real (kind=Rkind) :: x(ndim)
       integer i,j,ii,jj,ij,n(0:ndim)



       real (kind=Rkind) :: c,ep,em,th,xx
       real (kind=Rkind) :: poly_legendre

       double precision Req,beta,beta2
       NAMELIST /molec/Req,beta,beta2
       logical begin
       data begin/.true./
       SAVE

       IF (begin) THEN
         Req=ZERO
         beta=ONE
         beta2=ONE
         read(5,molec)
         write(out_unitp,molec)
         begin=.FALSE.
       END IF
! fin     initialisation la premiere fois

!      -----------------------------------------------
!      determine les indices ii et jj en fonction de
!      l'indice double ij et du tableau n(.)

       j=ij/n(2)
       i=ij-j*n(2)
       IF (i .EQ. 0) THEN
         i=n(2)
         j=j-1
       END IF

       ii=i
       jj=j
!      -----------------------------------------------

       c  = cos(x(1))


!      beta2 = beta * Req*Req
       xx = x(2)
       th = exp(-beta*xx-beta2/xx)

       v2_old = poly_legendre(c,ii,0) * x(2)**jj * th

!      write(out_unitp,*) ij,ii,jj,x,v2

       RETURN
       end function v2_old
!================================================================
!    fonction v3(x,ndim,ij,n)
!    x = cos(x(1) et R=x(2)
!    Pi0(x) * v7
!================================================================

       FUNCTION v3(x,ndim,ij,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v3 ! function

       integer ndim
       real (kind=Rkind) :: x(ndim)
       integer i,j,ii,jj,ij,n(0:ndim)

       real (kind=Rkind) :: c,th,ep,em

       real (kind=Rkind) :: poly_legendre,v7

       double precision Req,beta,beta2
       NAMELIST /molec/Req,beta,beta2
       logical begin
       data begin/.true./
       SAVE

       IF (begin) THEN
         Req=ZERO
         beta=ONE
         beta2=ONE
         read(5,molec)
         write(out_unitp,molec)
         begin=.FALSE.
       END IF
! fin     initialisation la premiere fois


       j=ij/n(2)
       i=ij-j*n(2)
       IF (i .EQ. 0) THEN
         i=n(2)
         j=j-1
       END IF

       ii=i
       jj=j

       c=cos(x(1))

!       ep = exp(beta*(x(2)-Req))
!       em = exp(-beta*(x(2)-Req))
!      th = ( (ep-em)/(ep+em) + ONE)*HALF * exp(-beta2*x(2))
       th = exp(-beta2*x(2))

       v3 = poly_legendre(c,ii,0) * x(2)**jj * th
!      write(out_unitp,*) ij,ii,jj,x(1),x(2),v3


       RETURN
       end function v3
!================================================================
!    fonction v_poly(x,i) id v5
!================================================================

       FUNCTION v_DirProd(x,ndim,ii,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v_DirProd ! function


       integer max_fit,max_dim
       parameter (max_fit=2000,max_dim=6)
       integer ind_n(max_dim)
       integer ijk(max_dim,max_fit)

       integer ndim
       real (kind=Rkind) :: z,x(ndim)
       integer type_v(max_dim),n_v(0:1)
       integer n(0:ndim),ii,kk
       integer i,max_p,i_func

       real (kind=Rkind) :: v ! recursive function

       NAMELIST /dirprod_typ/type_v

       logical begin
       data begin/.true./
       SAVE ijk,begin,type_v


       IF (begin) THEN


!        STOP
         type_v = (/ 15,15,15,15,15,15 /)
         read(5,dirprod_typ)
         write(out_unitp,dirprod_typ)

!        determine un tableau ijk(i,i_func)
!        qui donne un indice pour la dimension i en fonction de l'ordre de la fonction i_func
!        initialisation des ind_n(i)
         n(0)=1
         ind_n(:)=0
         max_p = 0
         DO i=1,ndim
           ind_n(i)=0
           n(0)=n(0)*n(i)
           IF (max_p .LT. n(i)) max_p = n(i)
         END DO
         ind_n(ndim)=-1
         write(out_unitp,*) 'nn,max_p',n(0),max_p


         kk=0

         write(out_unitp,*) 'max_p',max_p
         DO i_func=1,n(0)

!          determine les indices (ind_n) en fonction de i_func
           ind_n(ndim)=ind_n(ndim)+1
           DO i=ndim,1,-1
             IF (ind_n(i) .EQ. n(i)) THEN
               ind_n(i)=0
               IF (i .NE. 0) ind_n(i-1)=ind_n(i-1)+1
             END IF
           END DO

           IF (sum(ind_n) < max_p) THEN
             kk = kk+1

             DO i=1,ndim
               ijk(i,kk)=ind_n(i)
             END DO

             write(out_unitp,11) kk,(ijk(i,kk),i=1,ndim)
 11          format(10i3)

           END IF
         END DO

         n(0) = kk

         IF (n(0) .GT. max_fit) THEN
           write(out_unitp,*) ' ERROR : n(0) > max_fit',n(0),max_fit
           STOP
         END IF

         begin=.FALSE.
       END IF

       IF (ii .GT. n(0)) THEN
         write(out_unitp,*) ' ERROR in v_poly or v5'
         write(out_unitp,*) ' number of function :',n(0)
         write(out_unitp,*) ' You want the function :',ii
         STOP
       END IF




       z=ONE
       DO i=1,ndim
         z = z * v(x(i),1,ijk(i,ii)+1,n_v,type_v(i))
       END DO
       v_DirProd=z

       RETURN
       end function v_DirProd
!================================================================
!    fonction v4(x,i)
!================================================================

       FUNCTION v4(x,ndim,i,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v4 ! function

       integer ndim
       real (kind=Rkind) :: x(ndim)
       integer n(0:ndim),i
       integer :: max_fit,max_dim
       parameter (max_fit=2000,max_dim=6)
       integer ind_n(max_dim)
       integer ijk(max_dim,max_fit)

       integer k,i1,i2,ii1
       real (kind=Rkind) :: phi1(1),phi2(1)

       character (len=8) :: name

!      ---------------------------------------------------------------
       real (kind=Rkind) :: v26,v22,v15,v12
!      real (kind=Rkind) :: v30
!      real (kind=Rkind) :: v27
!      ---------------------------------------------------------------

       logical begin,all
       data begin/.true./
       SAVE begin,ijk


       IF (begin) THEN
         IF (n(0) == 0 ) n(0) = n(1)*n(2)
         all = .FALSE.
         write(out_unitp,*) 'n',n
         k = 0
         i1 = 0

         IF (.NOT. all) THEN
         DO
           IF (k == n(0) ) EXIT

!          write(out_unitp,*) 'i1,mod(i1,12)',i1,mod(i1,12)
           IF (mod(i1,12) == 0 ) THEN
!            cos(n*3 ) * PLeg(paire)
             name='cos'
             DO i2=0,n(2)-1,2
               k = k + 1
               ijk(1,k)=i1+1
               ijk(2,k)=i2+1
               write(out_unitp,*) k,name,(i1+1)/2,' * PLeg',i2
               IF (k == n(0) ) EXIT
             END DO
           ELSE IF (mod(i1,12) == 6) THEN
!            sin(n*6 ) * PLeg(impaire)
             name='cos'
             DO i2=1,n(2)-1,2
               k = k + 1
               ijk(1,k)=i1+1
               ijk(2,k)=i2+1
               write(out_unitp,*) k,name,(i1+1)/2,' * PLeg',i2
               IF (k == n(0) ) EXIT
             END DO

           END IF

           i1=i1+1

         END DO

         ELSE

         DO k=1,n(0)
           i1=(k-1)/n(2)
           i2=k-n(2)*i1-1
           ii1=i1
           IF (mod(ii1,2) == 1) name='sin'
           IF (mod(ii1,2) == 0) name='cos'
           ijk(1,k)=ii1+1
           ijk(2,k)=i2+1
           write(out_unitp,*) k,name,(ii1+1)/2,' * Leg',i2

         END DO
         END IF

         begin=.FALSE.
       END IF

       phi1=x(1)
       phi2=x(2)
       phi2=cos(x(2))



       i1=(i-1)/n(2)+1
       i2=i-n(2)*(i1-1)


       v4 = v26(phi1,1,ijk(1,i),n)*v12(phi2,1,ijk(2,i),n)
!      v4 = v26(phi1,1,i1,1)*v12(phi2,1,i2,1)




       end function v4

!================================================================
!    fonction v_poly(x,i) id v5
!================================================================

       FUNCTION v_poly(x,ndim,ii,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v_poly ! function


       integer max_fit,max_dim
       parameter (max_fit=2000,max_dim=6)
       integer ind_n(max_dim)
       integer ijk(max_dim,max_fit)

       integer ndim
       real (kind=Rkind) :: z,x(ndim)
       integer n(0:ndim),ii,kk
       integer i,max_p,i_func

       logical begin
       data begin/.true./
       SAVE ijk,begin


       IF (begin) THEN

!        determine un tableau ijk(i,i_func)
!        qui donne un indice pour la dimension i en fonction de l'ordre de la fonction i_func
!        initialisation des ind_n(i)
         n(0)=1
         ind_n(:)=0
         max_p = 0
         DO i=1,ndim
           ind_n(i)=0
           n(0)=n(0)*n(i)
           IF (max_p .LT. n(i)) max_p = n(i)
         END DO
         ind_n(ndim)=-1
         write(out_unitp,*) 'nn,max_p',n(0),max_p


         kk=0

         write(out_unitp,*) 'max_p',max_p
         DO i_func=1,n(0)

!          determine les indices (ind_n) en fonction de i_func
           ind_n(ndim)=ind_n(ndim)+1
           DO i=ndim,1,-1
             IF (ind_n(i) .EQ. n(i)) THEN
               ind_n(i)=0
               IF (i .NE. 0) ind_n(i-1)=ind_n(i-1)+1
             END IF
           END DO

           IF (sum(ind_n) < max_p) THEN
             kk = kk+1

             DO i=1,ndim
               ijk(i,kk)=ind_n(i)
             END DO

             write(out_unitp,11) kk,(ijk(i,kk),i=1,ndim)
 11          format(10i3)

           END IF
         END DO

         n(0) = kk

         IF (n(0) .GT. max_fit) THEN
           write(out_unitp,*) ' ERROR : n(0) > max_fit',n(0),max_fit
           STOP
         END IF

         begin=.FALSE.
       END IF

       IF (ii .GT. n(0)) THEN
         write(out_unitp,*) ' ERROR in v_poly or v5'
         write(out_unitp,*) ' number of function :',n(0)
         write(out_unitp,*) ' You want the function :',ii
         STOP
       END IF




       z=ONE
       DO i=1,ndim

         z = z * x(i)**ijk(i,ii)

       END DO
       v_poly=z

       RETURN
       end function v_poly
!================================================================
!    fonction v6(x,i) 1 D   Huffaker
!================================================================

       FUNCTION v6(x,ndim,i,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v6 ! function

       integer ndim
      real (kind=Rkind) :: x(ndim),m
      integer n(0:ndim)
       integer i

      double precision Req,beta
      NAMELIST /molec/Req,beta
      logical begin
      data begin/.true./
      SAVE

       IF (begin) THEN
         Req=1.40104_Rkind
         beta=ONE
         read(5,molec)
         write(out_unitp,molec)
         begin=.FALSE.
       END IF
! fin     initialisation la premiere fois

       m = ONE - exp(-beta*(x(1)-Req))
       v6= m**(i+1)


       RETURN
       end function v6
!================================================================
!    fonction v7(x,i) 1 BO
!================================================================

       FUNCTION v7(x,ndim,i,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v7 ! function

       integer ndim
      real (kind=Rkind) :: x(ndim),m
      integer n(0:ndim)
       integer i

      double precision Req,beta,beta2
      NAMELIST /molec/Req,beta,beta2
      logical begin
      data begin/.true./
      SAVE

       IF (begin) THEN
         Req=1.40104_Rkind
         beta=ONE
         beta2=ZERO
         read(5,molec)
         write(out_unitp,molec)
         begin=.FALSE.
       END IF
! fin     initialisation la premiere fois


!      m = exp(-beta*(x(1)-Req)-beta2/(x(1)-Req))
       m = exp(-beta*(x(1)-Req))

!      write(out_unitp,*) 'x(1),Req,beta,m',x(1),Req,beta,m

       v7= m**(i-1)

!      write(out_unitp,*) i,x(1),m,v7


       RETURN
       end function v7
!================================================================
!    fonction v8(x,i) 1 D (x-Req)/(c(1)*x+c(2)*Req)
!     c(1)=0 c(2)=1 : Dunham
!     c(1)=1 c(2)=0 : SPF
!     c(1)=0 c(2)=1 : Ogilvie
!================================================================

       FUNCTION v8(x,ndim,i,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v8 ! function

       integer ndim
      real (kind=Rkind) :: x(ndim),m
      integer n(0:ndim)
       integer i

      double precision Req,c(2)
      NAMELIST /molec/Req,c
      logical begin
      data begin/.true./
      SAVE

       IF (begin) THEN
         Req=1.40104_Rkind
         C(1)=ZERO
         C(2)=ONE
         read(5,molec)
         write(out_unitp,molec)
         begin=.FALSE.
       END IF
! fin     initialisation la premiere fois

       m = (x(1)-Req)/(c(1)*x(1)+c(2)*Req)
       v8= m**(i+1)


       RETURN
       end function v8

!================================================================
!    fonction v9(x,i) 1 D x**i*(x-1)**k
!================================================================

       FUNCTION v9(x,ndim,i,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v9 ! function

       integer ndim
      real (kind=Rkind) :: x(ndim)
      integer n(0:ndim)
       integer i,k

      NAMELIST /molec/k
      logical begin
      data begin/.true./
      SAVE

       IF (begin) THEN
         read(5,molec)
         write(out_unitp,molec)
         begin=.FALSE.
       END IF
! fin     initialisation la premiere fois

       v9= x(1)**i * (x(1)-ONE)**k


       RETURN
       end function v9
!================================================================
!    fonction v10(x,i) 1 D 1/x**(k*i)
!================================================================

       FUNCTION v10(x,ndim,i,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v10 ! function

       integer ndim
      real (kind=Rkind) :: x(ndim)
      integer n(0:ndim)
       integer i,k

      NAMELIST /molec/k
      logical begin
      data begin/.true./
      SAVE

       IF (begin) THEN
         read(5,molec)
         write(out_unitp,molec)
         begin=.FALSE.
       END IF
! fin     initialisation la premiere fois

       v10 = 1/x(1)**((i-1)*k)

!      write(out_unitp,*) 'v10=',v10


       RETURN
       end function v10
!================================================================
!    fonction v11(x,i) x^(i-1)*exp(-beta(x-xeq))
!================================================================

       FUNCTION v11(x,ndim,i,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v11 ! function

       integer ndim
      real (kind=Rkind) :: x(ndim),m
      integer n(0:ndim)
       integer i

      double precision Req,beta
      NAMELIST /molec/Req,beta
      logical begin
      data begin/.true./
      SAVE

       IF (begin) THEN
         Req=1.40104_Rkind
         beta=ONE
         read(5,molec)
         write(out_unitp,molec)
         begin=.FALSE.
       END IF
! fin     initialisation la premiere fois


!      IF (i .EQ. 1) THEN
!        v11 = ONE
!      ELSE
!        m = exp(-beta*(x(1)-Req))
!        v11= x(1)**(i-1) * m
!      END IF
       m = exp(-beta*(x(1)-Req))
       v11= x(1)**(i-1) * m

!      write(out_unitp,*) i,x(1),m,v7


       RETURN
       end function v11



!================================================================
!   calcule la valeur d'un polynome de Legendre (n1-1),0
!   pour un xx ( -1 =< x =< 1 )
!    v12 poly legendre
!    v13 poly legendre paire
!    v14 poly legendre impaire
!================================================================
       FUNCTION v12(x,ndim,n1,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v12 ! function

       integer ndim
       real (kind=Rkind) :: x(ndim),xx
       integer n1,n(0:ndim)
       real (kind=Rkind) :: poly_legendre
       integer l,m

       xx=x(1)
       l = n1
       m=0

       v12 = poly_legendre(xx,l,m)

       end function v12
!---------------------------------------------------
       FUNCTION v13(x,ndim,n1,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v13 ! function

       integer ndim
       real (kind=Rkind) :: x(ndim),xx
       integer n1,n(0:ndim)
       real (kind=Rkind) :: poly_legendre
       integer l,m

       xx=x(1)
       l = n1+n1-1
       m=0

!      write(out_unitp,*) xx,l,m

       v13 = poly_legendre(xx,l,m)

       end function v13
!---------------------------------------------------
       FUNCTION v14(x,ndim,n1,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v14 ! function

       integer ndim
       real (kind=Rkind) :: x(ndim),xx
       integer n1,n(0:ndim)
       real (kind=Rkind) :: poly_legendre
       integer l,m

       xx=x(1)
       l = n1+n1
       m=0

       v14 = poly_legendre(xx,l,m)

       end function v14
!================================================================
!   calcule la valeur d'un polynome x^i
!    v15 x^i
!    v16 x^i paire
!    v17 x^i impaire
!    v18 x^i paire sans i=0
!    v19 x^i sans i=0
!================================================================
       FUNCTION v15(x,ndim,n1,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v15 ! function

       integer ndim
       real (kind=Rkind) :: x(ndim),xx
       integer n(0:ndim)
       integer i,n1

       xx=x(1)
       i=n1-1

       v15 = xx**i

!      write(out_unitp,*) 'v15:',xx,i

       end function v15
!      ---------------------------------------------------------
       FUNCTION v16(x,ndim,n1,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v16 ! function

       integer ndim
       real (kind=Rkind) :: x(ndim),xx
       integer n(0:ndim)
       integer i,n1

       xx=x(1)
       i=n1-1

       v16 = xx**(i+i)

       end function v16
!      ---------------------------------------------------------
       FUNCTION v17(x,ndim,n1,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v17 ! function

       integer ndim
       real (kind=Rkind) :: x(ndim),xx
       integer n(0:ndim)
       integer i,n1

       xx=x(1)
       i=n1-1

       v17 = xx**(i+i+1)
       end function v17
!      ---------------------------------------------------------
       FUNCTION v18(x,ndim,n1,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v18 ! function

       integer ndim
       real (kind=Rkind) :: x(ndim),xx
       integer n(0:ndim)
       integer i,n1

       xx=x(1)
       i=n1-1

       v18 = xx**(i+i+2)

       end function v18
!      ---------------------------------------------------------
       FUNCTION v19(x,ndim,n1,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v19 ! function

       integer ndim
       real (kind=Rkind) :: x(ndim),xx
       integer n(0:ndim)
       integer i,n1

       xx=x(1)
       i=n1

       v19 = xx**i

       end function v19
!      ---------------------------------------------------------
!================================================================
!   calcule la valeur d'un polynome de Legendre (n1-1),0
!   pour un x ( 0 =< x =< pi )
!    v22 poly legendre
!    v23 poly legendre paire
!    v24 poly legendre impaire
!================================================================
       FUNCTION v22(x,ndim,n1,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v22 ! function

       integer ndim
       real (kind=Rkind) :: x(ndim),xx
       integer n1,n(0:ndim)
       real (kind=Rkind) :: poly_legendre
       integer l,m

       xx=cos(x(1))
       l = n1
       m=0

       v22 = poly_legendre(xx,l,m)

       end function v22
!---------------------------------------------------
       FUNCTION v23(x,ndim,n1,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v23 ! function

       integer ndim
       real (kind=Rkind) :: x(ndim),xx
       integer n1,n(0:ndim)
       real (kind=Rkind) :: poly_legendre
       integer l,m

       xx=cos(x(1))
       l = n1+n1-1
       m=0

!      write(out_unitp,*) xx,l,m

       v23 = poly_legendre(xx,l,m)

       end function v23
!---------------------------------------------------
       FUNCTION v24(x,ndim,n1,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v24 ! function

       integer ndim
       real (kind=Rkind) :: x(ndim),xx
       integer n1,n(0:ndim)
       real (kind=Rkind) :: poly_legendre
       integer l,m

       xx=cos(x(1))
       l = n1+n1
       m=0

       v24 = poly_legendre(xx,l,m)

       end function v24
!================================================================
!    v26 serie de fourier
!    v27 serie de fourier paire
!    v28 serie de fourier impaire
!    v29 serie de fourier en 2x
!    v30 serie de fourier en 3x
!================================================================
       FUNCTION v26(x,ndim,n1,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v26 ! function

       integer ndim
       real (kind=Rkind) :: x(ndim),xx
       integer n1,n(0:ndim)

       integer ii

!---------------------------------------------------------------------
      real (kind=Rkind), parameter :: pi2=pi+pi

      real (kind=Rkind) :: sq2pi,sqpi
!---------------------------------------------------------------------

       sqpi = ONE/sqrt(pi)
       sq2pi = ONE/sqrt(pi2)

       xx=x(1)
       ii = n1/2
       xx = mod(xx*real(ii,kind=Rkind),pi2)

       IF (ii .EQ. 0) THEN
         v26 = sq2pi
       ELSE
         IF (mod(n1,2) .EQ. 0) THEN
           v26 = sin(xx) * sqpi
         ELSE
           v26 = cos(xx) * sqpi
         END IF
       END IF


       end function v26
       FUNCTION v27(x,ndim,n1,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v27 ! function

       integer ndim
       real (kind=Rkind) :: x(ndim),xx
       integer n1,n(0:ndim)

       integer ii

!---------------------------------------------------------------------
      real (kind=Rkind), parameter :: pi2=pi+pi

      real (kind=Rkind) :: sq2pi,sqpi
!---------------------------------------------------------------------

       sqpi = ONE/sqrt(pi)
       sq2pi = ONE/sqrt(pi+pi)

       xx=x(1)
       ii = n1 - 1
       xx = mod(xx*real(ii,kind=Rkind),pi2)

       IF (ii .EQ. 0) THEN
         v27 = sq2pi
       ELSE
         v27 = cos(xx) * sqpi
       END IF


       end function v27
       FUNCTION v28(x,ndim,n1,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v28 ! function

       integer ndim
       real (kind=Rkind) :: x(ndim),xx
       integer n1,n(0:ndim)

       integer ii

!---------------------------------------------------------------------
      real (kind=Rkind), parameter :: pi2=pi+pi

      real (kind=Rkind) :: sq2pi,sqpi
!---------------------------------------------------------------------

       sqpi = ONE/sqrt(pi)
       sq2pi = ONE/sqrt(pi+pi)

       xx=x(1)
       ii = n1
       xx = mod(xx*real(ii,kind=Rkind),pi2)

       v28 = sin(xx) * sqpi


       end function v28
       FUNCTION v29(x,ndim,n1,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v29 ! function

       integer ndim
       real (kind=Rkind) :: x(ndim),xx
       integer n1,n(0:ndim)

       integer ii

!     ----------------------------------------------------------------
      real (kind=Rkind), parameter :: pi2=pi+pi

      real (kind=Rkind) :: sq2pi,sqpi
!     ----------------------------------------------------------------

       sqpi = ONE/sqrt(pi)
       sq2pi = ONE/sqrt(pi+pi)

       xx=x(1)+x(1)
       ii = n1/2
       xx = mod(xx*real(ii,kind=Rkind),pi2)

       IF (ii .EQ. 0) THEN
         v29 = sq2pi
       ELSE
         IF (mod(n1,2) .EQ. 0) THEN
           v29 = sin(xx) * sqpi
         ELSE
           v29 = cos(xx) * sqpi
         END IF
       END IF
       end function v29
       FUNCTION v30(x,ndim,n1,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v30 ! function

       integer ndim
       real (kind=Rkind) :: x(ndim),xx
       integer n1,n(0:ndim)

       integer ii

!     ----------------------------------------------------------------
      real (kind=Rkind), parameter :: pi2=pi+pi

      real (kind=Rkind) :: sq2pi,sqpi
!     ----------------------------------------------------------------

       sqpi = ONE/sqrt(pi)
       sq2pi = ONE/sqrt(pi+pi)

       xx=x(1)+x(1)+x(1)
       ii = n1/2
       xx = mod(xx*real(ii,kind=Rkind),pi2)

       IF (ii .EQ. 0) THEN
         v30 = sq2pi
       ELSE
         IF (mod(n1,2) .EQ. 0) THEN
           v30 = sin(xx) * sqpi
         ELSE
           v30 = cos(xx) * sqpi
         END IF
       END IF
       end function v30

       FUNCTION v41(x,ndim,n1,n)
       USE mod_system
       IMPLICIT NONE
       real (kind=Rkind) :: v41 ! function

       integer ndim
       real (kind=Rkind) :: x(ndim),xx
       integer n1,n(0:ndim)

       integer ii

!---------------------------------------------------------------------
      real (kind=Rkind), parameter :: pi2=pi+pi

      real (kind=Rkind) :: sq2pi,sqpi
!---------------------------------------------------------------------

       sqpi = ONE/sqrt(pi)
       sq2pi = ONE/sqrt(pi+pi)

       xx=x(1)
       ii = n1 - 1
       xx = mod(xx*real(ii,kind=Rkind),pi2)

       IF (ii .EQ. 0) THEN
         v41 = sq2pi
       ELSE
         v41 = cos(xx) * sqpi
       END IF


       END
       FUNCTION v51(x,ndim,n1,n)
       USE mod_system
       IMPLICIT NONE
       real (kind=Rkind) :: v51 ! function

       integer ndim
       real (kind=Rkind) :: x(ndim),xx
       integer n1,n(0:ndim)

       integer ii

!---------------------------------------------------------------------
      real (kind=Rkind), parameter :: pi2=pi+pi

      real (kind=Rkind) :: sq2pi,sqpi
!---------------------------------------------------------------------

       sqpi = ONE/sqrt(pi)
       sq2pi = ONE/sqrt(pi+pi)

       xx=x(1)
       ii = n1
       xx = mod(xx*real(ii,kind=Rkind),pi2)

       v51= sin(xx) * sqpi


       END
       FUNCTION v42(x,ndim,n1,n)
       USE mod_system
       IMPLICIT NONE
       real (kind=Rkind) :: v42 ! function

       integer ndim
       real (kind=Rkind) :: x(ndim),xx
       integer n1,n(0:ndim)

       integer ii

!---------------------------------------------------------------------
      real (kind=Rkind), parameter :: pi2=pi+pi

      real (kind=Rkind) :: sq2pi,sqpi
!---------------------------------------------------------------------

       sqpi = ONE/sqrt(pi)
       sq2pi = ONE/sqrt(pi+pi)

       xx=x(1)
       ii = 2*(n1 - 1)
       xx = mod(xx*real(ii,kind=Rkind),pi2)

       IF (ii .EQ. 0) THEN
         v42 = sq2pi
       ELSE
         v42 = cos(xx) * sqpi
       END IF


       END
       FUNCTION v52(x,ndim,n1,n)
       USE mod_system
       IMPLICIT NONE
       real (kind=Rkind) :: v52 ! function

       integer ndim
       real (kind=Rkind) :: x(ndim),xx
       integer n1,n(0:ndim)

       integer ii

!---------------------------------------------------------------------
      real (kind=Rkind), parameter :: pi2=pi+pi

      real (kind=Rkind) :: sq2pi,sqpi
!---------------------------------------------------------------------

       sqpi = ONE/sqrt(pi)
       sq2pi = ONE/sqrt(pi+pi)

       xx=x(1)
       ii = 2*n1
       xx = mod(xx*real(ii,kind=Rkind),pi2)

       v52= sin(xx) * sqpi


       END
       FUNCTION v43(x,ndim,n1,n)
       USE mod_system
       IMPLICIT NONE
       real (kind=Rkind) :: v43 ! function

       integer ndim
       real (kind=Rkind) :: x(ndim),xx
       integer n1,n(0:ndim)

       integer ii

!---------------------------------------------------------------------
      real (kind=Rkind), parameter :: pi2=pi+pi

      real (kind=Rkind) :: sq2pi,sqpi
!---------------------------------------------------------------------

       sqpi = ONE/sqrt(pi)
       sq2pi = ONE/sqrt(pi+pi)

       xx=x(1)
       ii = 3*(n1 - 1)
       xx = mod(xx*real(ii,kind=Rkind),pi2)

       IF (ii .EQ. 0) THEN
         v43 = sq2pi
       ELSE
         v43 = cos(xx) * sqpi
       END IF


       END
       FUNCTION v53(x,ndim,n1,n)
       USE mod_system
       IMPLICIT NONE
       real (kind=Rkind) :: v53 ! function

       integer ndim
       real (kind=Rkind) :: x(ndim),xx
       integer n1,n(0:ndim)

       integer ii

!---------------------------------------------------------------------
      real (kind=Rkind), parameter :: pi2=pi+pi

      real (kind=Rkind) :: sq2pi,sqpi
!---------------------------------------------------------------------

       sqpi = ONE/sqrt(pi)
       sq2pi = ONE/sqrt(pi+pi)

       xx=x(1)
       ii = 3*n1
       xx = mod(xx*real(ii,kind=Rkind),pi2)

       v53= sin(xx) * sqpi


       END
!================================================================
!   calcule la valeur d'un polynome d'Hermite*exp
!   pour un x ( 0 =< x =< pi )
!    v32 poly legendre
!    v33 poly legendre paire
!    v34 poly legendre impaire
!================================================================
       FUNCTION v32(x,ndim,n1,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v32 ! function

       integer ndim
       real (kind=Rkind) :: x(ndim),xx
       integer n1,n(0:ndim)
       real (kind=Rkind) :: poly_Hermite
       integer l

       xx=x(1)
       l = n1-1

       v32 = poly_Hermite(xx,l) * exp(-xx*xx*HALF)

       end function v32
!---------------------------------------------------
       FUNCTION v33(x,ndim,n1,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v33 ! function

       integer ndim
       real (kind=Rkind) :: x(ndim),xx
       integer n1,n(0:ndim)
       real (kind=Rkind) :: poly_Hermite
       integer l

       xx=x(1)
       l = n1+n1-2

       v33 = poly_Hermite(xx,l) * exp(-xx*xx*HALF)

       end function v33
!---------------------------------------------------
       FUNCTION v34(x,ndim,n1,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: v34 ! function

       integer ndim
       integer l
       real (kind=Rkind) :: x(ndim),xx
       integer n1,n(0:ndim)
       real (kind=Rkind) :: poly_Hermite

       xx=x(1)
       l = n1+n1-1

       v34 = poly_Hermite(xx,l) * exp(-xx*xx*HALF)

       end function v34
!================================================================
!    fonction v_inter(x)
!
!    v6  1D Huffaker
!    v7  1D BO
!    v8  1D (x-Req)/(c(1)*x+c(2)*Req)
!     c(1)=0 c(2)=1 : Dunham
!     c(1)=1 c(2)=0 : SPF
!     c(1)=0 c(2)=1 : Ogilvie
!    v9  1D x**i*(x-1)**k
!
!    v12 1D poly legendre          en cos(x)
!    v13 1D poly legendre paire    en cos(x)
!    v14 1D poly legendre impaire  en cos(x)
!
!    v15 1D x^n
!    v16 1D x^n paire
!    v17 1D x^n impaire
!    v18 1D x^n paire sans n=0
!    v19 1D x^n sans n=0
!
!    v22 1D poly legendre          en x
!    v23 1D poly legendre paire    en x
!    v24 1D poly legendre impaire  en x
!
!    v26 1D serie de fourier en x
!================================================================


!================================================================
!    fonction vgene_inter(x,ndim) 1 D
!================================================================
       FUNCTION vgene_inter(x,ndim,v_typ,F,nn,n)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: vgene_inter ! function

       integer   ndim
       real (kind=Rkind) ::    x(ndim)
       integer   nn,n(0:ndim)
       real (kind=Rkind) ::    z
       real (kind=Rkind) ::    F(nn)
       integer   kl

!      - function ----------------------------------------------
       real (kind=Rkind) :: v_typ
!      ---------------------------------------------------------

!      write(out_unitp,*) 'BEGINING vgene_inter',ndim,nn,n


       z=ZERO
       DO kl=1,nn
         z = z + F(kl) * v_typ(x,ndim,kl,n)
!        write(out_unitp,*) z,F(kl),ndim,nn,n
       END DO

       vgene_inter = z

!      write(out_unitp,*) 'END vgene_inter',x,z,ndim,nn,n

       end function vgene_inter
!================================================================
!    fonction vgene_inter(x,ndim) 1 D
!================================================================
       FUNCTION vgene_inter2(x,ndim,F,nn,n,ntyp)
       USE mod_system
       IMPLICIT NONE

       real (kind=Rkind) :: vgene_inter2 ! function

       integer   ndim
       real (kind=Rkind) ::    x(ndim)
       integer   ntyp,nn,n(0:ndim)
       real (kind=Rkind) ::    F(nn)

       integer   kl
       real (kind=Rkind) ::    z

!      - function ----------------------------------------------
       real (kind=Rkind) :: v
!      ---------------------------------------------------------

!      write(out_unitp,*) 'BEGINING vgene_inter2',ndim,nn,n
!      write(out_unitp,*) 'F(:)',F
!      write(out_unitp,*) 'x',x


       z=ZERO
       DO kl=1,nn
         z = z + F(kl) * v(x,ndim,kl,n,ntyp)
         !write(out_unitp,*) z,F(kl),ndim,nn,n
         !CALL flush_perso(6)
       END DO

       vgene_inter2 = z

!      write(out_unitp,*) 'END vgene_inter2',x,z,ndim,nn,n

       end function vgene_inter2
!
!================================================================
!
! subroutine splin1d
!
!
!================================================================
!
!
      SUBROUTINE SPLINE_1D(X,Y,N,Y2)
      USE mod_system
      IMPLICIT NONE


      integer    N
      real(kind=Rkind)     X(N),Y(N),Y2(N),U(N)
      real(kind=Rkind)     YP1,YPN

      integer    I,K
      real(kind=Rkind)     SIG,P,QN,UN


      YP1=(Y(1)-Y(2))/(X(1)-X(2))
      YPN=(Y(N)-Y(N-1))/(X(N)-X(N-1))

      IF (YP1.GT. TEN**30) THEN
        Y2(1)=ZERO
        U(1)=ZERO
      ELSE
        Y2(1)=-HALF
        U(1)=(THREE/(X(2)-X(1)))*((Y(2)-Y(1))/(X(2)-X(1))-YP1)
      ENDIF
      DO 11 I=2,N-1
        SIG=(X(I)-X(I-1))/(X(I+1)-X(I-1))
        P=SIG*Y2(I-1)+TWO
        Y2(I)=(SIG-ONE)/P
        U(I)=(SIX*((Y(I+1)-Y(I))/(X(I+1)-X(I))-(Y(I)-Y(I-1))            &
            /(X(I)-X(I-1)))/(X(I+1)-X(I-1))-SIG*U(I-1))/P
11    CONTINUE
      IF (YPN.GT.TEN**30) THEN
        QN=ZERO
        UN=ZERO
      ELSE
        QN=HALF
        UN=(THREE/(X(N)-X(N-1)))*(YPN-(Y(N)-Y(N-1))/(X(N)-X(N-1)))
      ENDIF
      Y2(N)=(UN-QN*U(N-1))/(QN*Y2(N-1)+ONE)
      DO 12 K=N-1,1,-1
        Y2(K)=Y2(K)*Y2(K+1)+U(K)
12    CONTINUE

      RETURN
      end subroutine SPLINE_1D

      FUNCTION SPLINT_1D(X,XA,YA,Y2A,N)
      USE mod_system
      IMPLICIT NONE

      real(kind=Rkind)  :: SPLINT_1D
      integer    N
      real(kind=Rkind)     XA(N),YA(N),Y2A(N)
      real(kind=Rkind)     X,Y

      integer    KLO,KHI,I,K
      real(kind=Rkind)     H,A,B

      KLO=1
      KHI=N
1     IF (KHI-KLO.GT.1) THEN
        K=(KHI+KLO)/2
        IF(XA(K).GT.X)THEN
          KHI=K
        ELSE
          KLO=K
        ENDIF
      GOTO 1
      ENDIF
      H=XA(KHI)-XA(KLO)
      IF (H.EQ.ZERO) STOP 'Bad XA input.'
      A=(XA(KHI)-X)/H
      B=(X-XA(KLO))/H
      Y=A*YA(KLO)+B*YA(KHI)+                                            &
            ((A**3-A)*Y2A(KLO)+(B**3-B)*Y2A(KHI))*(H**2)/SIX


      SPLINT_1D=Y
      RETURN
      end function SPLINT_1D

