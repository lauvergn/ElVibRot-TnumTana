c
C================================================================
C    calc_Op : calculation of the potential and dipolar matrices
c    mat_V(nb_be,nb_be) and mat_dip(nb_be,nb_be,3)
c    nb_be : nb of elctronic surface
c    Q are the coordinates in active order or syl order
c    dipolar calculation if calc_dip = T
C================================================================
      SUBROUTINE calcN_op(mat_V,mat_imV,mat_dip,nb_be,nb_Dip,
     *                    Qcart,nb_cart,mole,calc_dip,pot_cplx)

      use mod_system
      USE mod_Tnum
      USE mod_Lib_QTransfo
      implicit none

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole


      integer       :: nb_be,nb_cart,nb_Dip
      logical       :: calc_dip,pot_cplx
      real (kind=8) :: mat_V(nb_be,nb_be),mat_imV(nb_be,nb_be)
      real (kind=8) :: mat_dip(nb_be,nb_be,nb_Dip)
      real (kind=8) :: Qcart(nb_cart)

      real (kind=8) :: pot0,im_pot0
      real (kind=8) :: v1(3),r1,v2(3),r2,v3(3),r3
      real (kind=8) :: v1xv2(3),rv1xv2,v2xv3(3),rv2xv3
      real (kind=8) :: th2,th3,phi3,Qdyn(6)
      real (kind=8) :: dip(3)
      !# Third, calculate the energy from harmonic oscillator potential
      real (kind=8), parameter :: r1_eq = 2.2791110577179996_Rkind
      real (kind=8), parameter :: r2_eq = 2.0824950811500056_Rkind
      real (kind=8), parameter :: theta2_eq = 2.1237280405061139_Rkind
      real (kind=8), parameter :: r3_eq = 2.0824950811500056_Rkind
      real (kind=8), parameter :: theta3_eq = 2.1237280405061139_Rkind
      real (kind=8), parameter :: phi3_eq = 3.1415926535897931_Rkind
      real (kind=8), parameter :: Qeq(6) =
     *  (/ r1_eq,r2_eq,theta2_eq,r3_eq,theta3_eq,phi3_eq /)


      real (kind=8), parameter :: k_r1 =
     *    TWO*0.44369846935865775439112_Rkind
      real (kind=8), parameter :: k_r2 =
     *    TWO*0.16461456237399485491579_Rkind
      real (kind=8), parameter :: k_r3 =
     *    TWO*0.16461456237399485491579_Rkind
      real (kind=8), parameter :: k_theta2 =
     *     TWO*0.15598412935297364945164_Rkind
      real (kind=8), parameter :: k_theta3 =
     *     TWO*0.15598412935297364945164_Rkind
      real (kind=8), parameter :: k_phi3 =
     *     TWO*.036413510801518840509505_Rkind
      real (kind=8), parameter :: k(6) =
     *   (/k_r1,k_r2,k_theta2,k_r2,k_theta3,k_phi3 /)


      v1(:) = Qcart(4:6)-Qcart(1:3)
      v2(:) = Qcart(4:6)-Qcart(7:9)
      v3(:) = Qcart(4:6)-Qcart(10:12)
      r1    = sqrt(dot_product(v1,v1))
      r2    = sqrt(dot_product(v2,v2))
      r3    = sqrt(dot_product(v3,v3))
      th2   = acos(dot_product(v1,v2)/(r1 * r2))
      th3   = acos(dot_product(v1,v3)/(r1 * r3))

      CALL calc_cross_product(v1,r1,v2,r2,v1xv2,rv1xv2)
      CALL calc_cross_product(v2,r2,v3,r3,v2xv3,rv2xv3)


      phi3 = acos(-dot_product(v1xv2,v2xv3) / (rv1xv2*rv2xv3))

      Qdyn(:) = (/ r1,r2,th2,r3,th3,phi3 /)

      mat_V(1,1) = HALF*dot_product(k,(Qdyn-Qeq)**2)


      !write(6,*) 'Qdyn',Qdyn,mat_V(1,1)
      !STOP


      RETURN
      END
      SUBROUTINE pot_ori(V,Qdyn)
      USE mod_system
      USE mod_file
      implicit none

       integer,parameter :: ndim = 6
       real (kind=Rkind) :: Qdyn(ndim)
       real (kind=Rkind) :: Qref(ndim)
       real (kind=Rkind) :: DQ(ndim)
       real (kind=Rkind) :: V,vFnd

       integer,parameter :: max_points = 175
       real (kind=Rkind) :: F(max_points)
       integer :: tab_func(ndim,max_points)

       integer       ::  i,j,k,nb_columns,nb_func,nio

       logical :: begin
       data begin/.true./
       SAVE begin,nb_func,F,tab_func,Qref


c      GOTO 10
c---------------------------------------------------------------
c      initialisation la premiere fois
c$OMP  CRITICAL (pot0_CRIT)
       IF (begin) THEN
         CALL file_open2(name_file='pot_h2co_tot.txt',iunit=nio)
         read(nio,*) nb_func
         IF (max_points < nb_func) THEN
           write(6,*) ' ERROR in pot0'
           write(6,*) ' nb_func > max_points',nb_func,max_points
           STOP
         END IF

         DO j=1,nb_func
           read(nio,*) tab_func(1:ndim,j),F(j)
         END DO
         close(nio)
         begin=.FALSE.
         ! R1opt=R2Opt=1.099 Angs, RCO=1.203 Angs)
         ! (theta1-121.7*Pi/180.), (theta2-121.7*Pi/180.), (phi-Pi)
         Qref(:) = (/2.07680591407806069035d0,2.07680591407806069035d0,
     *             2.27333713797625751637d0,
     *             2.12406569967709909509d0,2.12406569967709909509d0,
     *             3.14159265358979323844d0 /)

c        write(6,*) 'Qref', Qref(:)
c        DO j=1,nb_func
c          write(6,*) j,': ',tab_func(1:ndim,j),F(j)
c        END DO
       END IF
c$OMP  END CRITICAL (pot0_CRIT)
c fin     initialisation la premiere fois
c---------------------------------------------------------------


c     FQ1 = (1 - Exp[-a1*(Q1 - 1.099/0.529178)]);
c     FQ2 = (1 - Exp[-a2*(Q2 - 1.099/0.529178)]);
c     FQ3 = (1 - Exp[-a3*(Q3 - 1.203/0.529178)]);
c     FQ4 = (Q4 - 121.7*N[Pi]/180.);
c     FQ5 = (Q5 - 121.7*N[Pi]/180.);
c     FQ6 = (Q6 - N[Pi]);

      DQ(:) = Qdyn(:) - Qref(:)

      DQ(1) = ONE - exp(-0.872313d0*DQ(1))
      DQ(2) = ONE - exp(-0.872313d0*DQ(2))
      DQ(3) = ONE - exp(-0.910102d0*DQ(3))


c     write(6,*) 'Qdyn',Qdyn
c     write(6,*) 'DQ',DQ

      V = ZERO
      DO j=1,nb_func
        vFnd =  F(j)
        DO i=1,ndim
        vFnd = vFnd * DQ(i)**tab_func(i,j)
        END DO
        V = V + vFnd
      END DO
c10   CONTINUE

c     write(66,*) Qdyn(:),V

       END
C================================================================
C    fonction im_pot0(x) imaginary part of pot0
C================================================================
       real*8 FUNCTION im_pot0(Qsym0)


       real*8 Qsym0(1)

       im_pot0 = 0.0

       END
C================================================================
C    fonction pot_rest(x) rest of the DL : pot0 + v2 + pot_rest
C================================================================
       real*8 FUNCTION pot_rest(Qact,Delta_Qact,nb_inact2n)

       real*8 Qact(1)
       integer nb_inact2n
       real*8 Delta_Qact(nb_inact2n)

       pot_rest = 0.0

      END
C================================================================
C    subroutine calculant le gradient
C================================================================
       SUBROUTINE d0d1d2_g(d0g,d1g,d2g,Qsym0,mole,deriv,num,step)

      USE mod_Tnum
      implicit none

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

      real (kind=8) :: d0g(mole%nb_inact2n)
      real (kind=8) :: d1g(mole%nb_inact2n,mole%nb_act1)
      real (kind=8) :: d2g(mole%nb_inact2n,mole%nb_act1,mole%nb_act1)

      real (kind=8) :: Qsym0(mole%nb_var)
      real (kind=8) :: step
      logical       :: deriv,num




      real (kind=8) :: Qact(mole%nb_act1)
      real (kind=8) :: d0f
      real (kind=8) :: d1f(mole%nb_act1)
      real (kind=8) :: d2f(mole%nb_act1,mole%nb_act1)
      integer       :: i,k,vi,ibo

      real (kind=8) :: d0gzmt(mole%nb_var)

c----- for debuging ----------------------------------
c     logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
c---------------------------------------------------------------------
      IF (debug) THEN
      write(6,*)
      write(6,*) 'BEGINNING d0d1d2_g'
      write(6,*) 'nb_var',mole%nb_var
      write(6,*) 'nb_act1',mole%nb_act1
      write(6,*) 'nb_inact22,nb_inact21',mole%nb_inact22,mole%nb_inact21
      write(6,*) 'nb_inact2n',mole%nb_inact2n
      write(6,*) 'deriv',deriv
      END IF

c---------------------------------------------------------------------
       Qact(1) = Qsym0(mole%liste_QactTOQsym(1))
c---------------------------------------------------------------------

      d0g(:) = 0.d0

c---------------------------------------------------------------------
       IF (debug) THEN
         write(6,*) 'd0g at Qact:',Qact
         write(6,*) d0g(:)
         write(6,*) 'END d0d1d2_g'
       END IF
c---------------------------------------------------------------------

       RETURN
       END
C================================================================
C    subroutine calculant la matrice hessienne
C    en fonction de x=cos(theta)
C================================================================
       SUBROUTINE sub_hessian (h)

       real*8 h

       h = 0.d0


       RETURN
       END
       SUBROUTINE H0_sym(h)

       real*8 h

       RETURN
       END
       SUBROUTINE d0d1d2_h(d0h,d1h,d2h,
     *                     Qsym0,mole,deriv,num,step)

      USE mod_Tnum
      implicit none

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

       real*8  Qact1,Qsym0(mole%nb_var)



       real*8 d0h(mole%nb_inact2n,mole%nb_inact2n)
       real*8 d1h(mole%nb_inact2n,mole%nb_inact2n)
       real*8 d2h(mole%nb_inact2n,mole%nb_inact2n)
       real*8 step
      logical       :: deriv,num


       RETURN
       END
C================================================================
C    analytical derivative (Qeq Qeq' Qeq" Qeq'") calculation
c    for the variable i_qsym
C================================================================
      SUBROUTINE calc_dnQflex(i_qsym)
      integer :: i_qsym
      STOP
      END 
      SUBROUTINE d0d1d2d3_Qeq(i_qsym,
     *                        d0req,d1req,d2req,d3req,
     *                        Qsym0,mole,nderiv)


      USE mod_Tnum
      implicit none

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

       integer i_qact,i_qsym
       real*8  Qact1,Qsym0(mole%nb_var)

       integer nderiv

       real*8  d0req
       real*8  d1req(mole%nb_act1)
       real*8  d2req(mole%nb_act1,mole%nb_act1)
       real*8  d3req(mole%nb_act1,mole%nb_act1,mole%nb_act1)
       integer vi


       real*8 dc0,dc1,dc2,dc3
       real*8 c_act,s_act,ss,sc,sss

       character*14 nom_i,nom
       logical exist

       integer   max_points,nb_inactb
       parameter (max_points=200)
       parameter (nb_inactb=10)
       real*8 F(max_points,nb_inactb)
       integer nn(nb_inactb)

       integer k_act,i_act,j_act,k,kl
       integer i_sym_act,i_type_act,i_qsym_inact

c----- function --------------------------------------
       real*8 poly_legendre
c----- function --------------------------------------

c----- for debuging ----------------------------------
      logical debug
      parameter (debug=.FALSE.)
c     parameter (debug=.TRUE.)
c----- for debuging ----------------------------------

       logical begin
       data begin/.true./
       SAVE begin,F,nn,i_type_act

c---------------------------------------------------------------------
      IF (debug) THEN
        write(6,*) 'BEGINNING d0d1d2d3_Qeq'
        write(6,*) 'nb_inact20,nb_act',mole%nb_inact20,mole%nb_act
        write(6,*) 'nb_var',mole%nb_var
        write(6,*) 'i_qsym',i_qsym
      END IF
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c      Qact value. Rq: only ONE active variable is possible
c---------------------------------------------------------------------
       IF (mole%nb_act1 .NE. 1) THEN
         write(6,*) ' ERROR : d0d1d2d3_Qeq'
         write(6,*) ' the number of Active variable'
         write(6,*) ' should be 1. But nb_act1 =',mole%nb_act1
         STOP
       END IF

c---------------------------------------------------------------------
c---------------------------------------------------------------------


c---------------------------------------------------------------------
       d0req = 0.d0
       d1req = 0.d0
       d2req = 0.d0
       d3req = 0.d0
c---------------------------------------------------------------------

       END
C================================================================
C    The tri component of the dipole moment.
C================================================================
       SUBROUTINE sub_dipole(dip,Q)

       real*8 Q(3)
       real*8 dip(3)


       dip(:) = 0.d0


       END
