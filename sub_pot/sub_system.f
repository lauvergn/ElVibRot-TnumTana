c
c================================================================
c    calc_Op : calculation of the potential and dipolar matrices
c    mat_V(nb_be,nb_be) and mat_dip(nb_be,nb_be,3)
c    nb_be : nb of elctronic surface
c    Q are the coordinates in active order or syl order
c    dipolar calculation if calc_dip = T
c================================================================
      SUBROUTINE calcN_op(mat_V,mat_imV,mat_ScalOp,nb_be,nb_ScalOp,
     *                    Q,nb_var,mole,
     *                    calc_ScalOp,pot_cplx)

      USE mod_Tnum
      USE mod_system
      IMPLICIT NONE

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole


      integer           :: nb_be,nb_ScalOp,nb_var
      logical           :: calc_ScalOp,pot_cplx
      real (kind=Rkind) :: mat_V(nb_be,nb_be),mat_imV(nb_be,nb_be)
      real (kind=Rkind) :: mat_ScalOp(nb_be,nb_be,nb_ScalOp)
      real (kind=Rkind) :: Q(nb_var)

      real (kind=Rkind) :: pot0,im_pot0


      IF (nb_be == 1 ) THEN
        CALL POTS(mat_V(1,1),Q(1),Q(2),Q(3))
        !mat_V(1,1) = pot0(Q)
        IF (pot_cplx) mat_imV(1,1) = im_pot0(Q)
        IF (calc_ScalOp) THEN
          CALL sub_dipole(mat_ScalOp(1,1,:),Q,mole)
        END IF
      ELSE
        write(6,*) ' ERROR in calc_op'
        write(6,*) ' more than ONE electronic surface is impossible'
        write(6,*) ' Rq: nb_be',nb_be
        STOP
      END IF

      RETURN
      END
c
C================================================================
C    fonction pot0(x) 3+9 D pour h2o en cartesiennes (calcul direct)
C================================================================
      FUNCTION pot0(Q)
      USE mod_system
      IMPLICIT NONE

      integer, parameter :: ndim=3
      real (kind=Rkind) :: pot0
      real (kind=Rkind) :: Q(ndim)
      real (kind=Rkind) :: DQ(ndim)
      real (kind=Rkind) :: G(ndim)

      real (kind=Rkind), parameter :: Qref(ndim) = 
     *  (/ 1.869713_Rkind,1.869713_Rkind,1.745803_Rkind /)

      real (kind=Rkind) :: hess(ndim,ndim)

c1          0.65408812100000002232        -0.03064317947323269911         0.03770002118955283199
c2         -0.03064317947323269911         0.65408822251663933933         0.03770005404468122767
c3          0.03770002118955283199         0.03770005404468122767         0.29688447831465542004

      hess(:,1) = (/
     1 0.65408812100000002232d0,-0.03064317947323269911d0,
     1    0.03770002118955283199d0/)
      hess(:,2) = (/
     2-0.03064317947323269911d0, 0.65408822251663933933d0,
     2      0.03770005404468122767d0/)
      hess(:,3) = (/
     3 0.03770002118955283199d0, 0.03770005404468122767d0,
     3        0.29688447831465542004d0 /)

      !write(6,*) hess
      DQ(:) = Q(:)-Qref(:)
      write(6,*) 'DQ',DQ
      pot0 = -74.9659012171d0 + HALF * dot_product(DQ,matmul(hess,DQ))

      END
C================================================================
C    subroutine calculant le gradient
C================================================================
      SUBROUTINE d0d1d2_g(d0g,d1g,d2g,Qsym0,mole,deriv,num,step)
      USE mod_Tnum
      USE mod_system
      IMPLICIT NONE

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

      real (kind=Rkind) :: d0g(mole%nb_inact2n)
      real (kind=Rkind) :: d1g(mole%nb_inact2n,mole%nb_act1)
      real (kind=Rkind) :: 
     *                d2g(mole%nb_inact2n,mole%nb_act1,mole%nb_act1)

      real (kind=Rkind) :: Qsym0(mole%nb_var)
      real (kind=Rkind) :: step
      logical       :: deriv,num


      d0g = 0.d0


      END
C================================================================
C    subroutine calculant la matrice hessienne en coordonnees cartesiennes
C    !!! il faut changer le paramtre n  (3*nb_at)
C    et le nom de file_FChk%name
C================================================================
      SUBROUTINE sub_hessian(hh)
      USE mod_file
      USE mod_system
      USE mod_OTF
      IMPLICIT NONE

      integer, parameter :: n = 9
      real (kind=Rkind) :: h(n,n),hh(n,n)
      integer  ::  err
      TYPE (param_file) :: file_pun
      integer  :: nio
      logical  :: located

      integer  ::  i,j,k,idum,jdum,nbligne,nbreste


      hh(:,:) = ZERO

      END
      SUBROUTINE H0_sym(h,n)
      USE mod_system
      IMPLICIT NONE
        integer       :: n
        real (kind=Rkind) :: h(n,n)
        integer       :: n1 = 9
        integer       :: n2 = 18


        real (kind=Rkind) :: d


        RETURN


        d = h(n1,n1)
        h(:,n1) = 0.d0
        h(n1,:) = 0.d0
        h(n1,n1) = d
        

        d = h(n2,n2)
        h(:,n2) = 0.d0
        h(n2,:) = 0.d0
        h(n2,n2) = d
        write(6,*) 'hessian sym'
        CALL ecriture(h,n,n,5,.TRUE.,n)
        
      END
C================================================================
C    fonction pot_rest(x)
C================================================================
      FUNCTION pot_rest(Qact,Delta_Qact,nb_inact2n)
      USE mod_system
      IMPLICIT NONE

       real (kind=Rkind) :: pot_rest
       real (kind=Rkind) :: Qact(1)
       integer       :: nb_inact2n
       real (kind=Rkind) :: Delta_Qact(nb_inact2n)

       pot_rest = 0.d0

       END
C================================================================
C    fonction im_pot0(x)
C================================================================
      FUNCTION im_pot0(Qsym0)
      USE mod_system
      IMPLICIT NONE

       real (kind=Rkind) :: im_pot0
       real (kind=Rkind) :: Qsym0(1)
       real (kind=Rkind) :: z

       z = 0.d0

       im_pot0 = z

       RETURN
       END
C================================================================
C    subroutine calculant la matrice hessienne
C    en fonction de x=cos(theta)
C================================================================
       SUBROUTINE d0d1d2_h(d0h,d1h,d2h,
     *                     Qsym0,mole,deriv,num,step)

      USE mod_Tnum
      USE mod_system
      IMPLICIT NONE

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

      
      real (kind=Rkind) ::  Qsym0(mole%nb_var)


      real (kind=Rkind) :: step
      logical deriv,num

      real (kind=Rkind) :: d0h
      real (kind=Rkind) :: d1h
      real (kind=Rkind) :: d2h

c----- for debuging ----------------------------------
      logical debug
      parameter (debug=.FALSE.)
c     parameter (debug=.TRUE.)
c---------------------------------------------------------------------
      IF (debug) THEN
      write(6,*)
      write(6,*) 'BEGINNING d0d1d2_h'
      END IF
c---------------------------------------------------------------------


      STOP 'd0d1d2_h'

      END
C================================================================
C    analytical derivative (Qeq Qeq' Qeq" Qeq'") calculation
c    for the variable i_qsym
C================================================================
      SUBROUTINE d0d1d2d3_Qeq(i_qsym,
     *                        d0req,d1req,d2req,d3req,
     *                        Qsym0,mole,nderiv)

      USE mod_Tnum
      USE mod_system
      IMPLICIT NONE

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

       integer i_qsym
       real (kind=Rkind) ::  Qsym0(mole%nb_var)

       integer nderiv

       real (kind=Rkind) ::  d0req
       real (kind=Rkind) ::  d1req
       real (kind=Rkind) ::  d2req
       real (kind=Rkind) ::  d3req


c----- for debuging ----------------------------------
      logical debug
      parameter (debug=.FALSE.)
c     parameter (debug=.TRUE.)
c---------------------------------------------------------------------
      IF (debug) THEN
        write(6,*) 'BEGINNING d0d1d2d3_Qeq'
        write(6,*) 'i_qsym',i_qsym
      END IF
c---------------------------------------------------------------------

      STOP 'd0d1d2d3_Qeq'

      RETURN
      END

C================================================================
C    analytical derivative (dnQflex : Qflex Qflex' Qflex" Qflex'") calculation
c    for the variable iq
C================================================================
      SUBROUTINE calc_dnQflex(iq,dnQflex,Qact,nb_act,nderiv,it)
      USE mod_system
      USE mod_dnSVM
      IMPLICIT NONE

       integer :: iq,nb_act
       real (kind=Rkind) ::  Qact(nb_act)
       integer :: nderiv,it
       TYPE (Type_dnS)   :: dnQflex
       STOP 'dnQflex'
       END
C================================================================
C    analytical derivative (dnQgene : Qgene Qgene' Qgene" Qgene'") calculation
C    for the variable iq_gene
C================================================================
      SUBROUTINE calc_dnQgene(iq_gene,dnQgene,Qgene,nb_Qgene,nderiv,it,
     *                        inTOout)
      USE mod_system
      USE mod_dnSVM
      IMPLICIT NONE

       integer           :: iq_gene,nb_Qgene
       real (kind=Rkind) :: Qgene(nb_Qgene)
       integer           :: nderiv,it
       TYPE (Type_dnS)   :: dnQgene
       logical           :: inTOout

       STOP 'calc_dnQgene'
       END
c
C================================================================
c    dipole read
C================================================================
      SUBROUTINE sub_dipole(dip,Q,mole)
      USE mod_Tnum
      USE mod_system
      IMPLICIT NONE

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

      real (kind=Rkind) :: Q(mole%nb_var)
      real (kind=Rkind) :: DQ(mole%nb_var)
      real (kind=Rkind) :: dip(3)

      real (kind=Rkind), parameter :: Qref(3) = 
     *  (/ 1.869713_Rkind,1.869713_Rkind,1.745803_Rkind /)

      real (kind=Rkind), parameter :: MuX = 0.51523633071858288d0
      real (kind=Rkind), parameter :: MuY = ZERO
      real (kind=Rkind), parameter :: MuZ = 0.43212653883837798d0
      real (kind=Rkind), parameter :: dMuX(3) = (/
     x   -2.89121043244283428d-002,-0.19662504130644246d0,
     x    0.11731395143422289d0/)
      real (kind=Rkind), parameter :: dMuY(3) = ZERO
      real (kind=Rkind), parameter :: dMuZ(3) = (/
     z -0.19456328715229212d0,    5.40468620138319714d-003,
     z -0.34043381267231121d0/)

      DQ(:) = Q(:)-Qref(:)


      dip(1) = MuX + dot_product(DQ,dMuX)
      dip(2) = MuY + dot_product(DQ,dMuY)
      dip(3) = MuZ + dot_product(DQ,dMuZ)

      !dip = dip * 0.5291772084d0

      END
      SUBROUTINE POTS(V,Q1,Q2,THETA)
 
C     Potential PJT2 due Polyansky, Jensen and Tennyson, 
C     J. Chem. Phys., 105, 6490-6497 (1996)
C     Update of Polyansky, Jensen and Tennyson, J Chem Phys 101, 7651 (1994))
C     Units: Hartree and Bohr
 
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
 
C     RZ = OH equilibrium value
C     RHO = equilibrium value of pi - bond angle(THETA)
 
      DATA TOANG/0.5291772/, CMTOAU/219474.624/
      DATA X1/1.0/
      DATA RHO1    /    75.50035308/
      DATA FA1     /     .00000000/
      DATA FA2     /18902.44193433/
      DATA FA3     /  1893.99788146/
      DATA FA4     /  4096.73443772/
      DATA FA5     /-1959.60113289/
      DATA FA6     /  4484.15893388/
      DATA FA7     /  4044.55388819/
      DATA FA8     / -4771.45043545/
      DATA FA9     /     0.00000000/
      DATA FA10    /     0.00000000/
      DATA RZ    /     .95792059/
      DATA A     /    2.22600000/
      DATA F1A1    /  -6152.40141181/
      DATA F2A1    / -2902.13912267/
      DATA F3A1    / -5732.68460689/
      DATA F4A1    /  953.88760833/
      DATA F11     / 42909.88869093/
      DATA F1A11   /  -2767.19197173/
      DATA F2A11   /  -3394.24705517/
      DATA F3A11   /     .00000000/
      DATA F13     /  -1031.93055205/
      DATA F1A13   /  6023.83435258/
      DATA F2A13   /     .00000000/
      DATA F3A13   /     .00000000/
      DATA F111    /     .00000000/
      DATA F1A111  /   124.23529382/
      DATA F2A111  /  -1282.50661226/
      DATA F113    /  -1146.49109522/
      DATA F1A113  /  9884.41685141/
      DATA F2A113  /  3040.34021836/ 
      DATA F1111   /  2040.96745268/
      DATA FA1111  /  .00000000/
      DATA F1113   /  -422.03394198/
      DATA FA1113  /-7238.09979404/
      DATA F1133   /     .00000000/
      DATA FA1133  /     .00000000/
      DATA F11111  / -4969.24544932/
      DATA f111111/  8108.49652354/
      DATA F71   /  90.00000000/
 
 
 
      data c1/50.0/,c2/10.0/,beta1/22.0/,beta2/13.5/,gammas/0.05/,
     *     gammaa/0.10/,delta/0.85/,rhh0/1.40/
                 RHO=RHO1*3.141592654/180.000000000
      fa11=0.0
      f1a3=f1a1
      f2a3=f2a1
      f3a3=f3a1
      f4a3=f4a1
      f33=f11
      f1a33=f1a11
      f2a33=f2a11
      f333=f111
      f1a333=f1a111
      f2a333=f2a111
      f133=f113
      f1a133=f1a113
      f2a133=f2a113
      f3333=f1111
      fa3333=fa1111
      f1333=f1113
      fa1333=fa1113
      f33333=f11111
      f333333 =f111111
      f73     =f71
 
C     Find value for DR and DS
      DR = TOANG*Q1 - RZ
      DS = TOANG*Q2 - RZ
 
C     Transform to Morse coordinates
      Y1 = X1 - EXP(-A * DR)
      Y3 = X1 - EXP(-A * DS)
 
C     transform to Jensens angular coordinate
      CORO = DCOS(THETA) + DCOS(RHO)
 
C     Now for the potential
      V0=(FA2+FA3*CORO+FA4*CORO**2+FA6*CORO**4+FA7*CORO**5)*CORO**2
      V0=V0+(FA8*CORO**6+FA5*CORO**3+FA9*CORO**7+FA10*CORO**8 )*CORO**2
      V0=V0+(                                    FA11*CORO**9 )*CORO**2
      FE1= F1A1*CORO+F2A1*CORO**2+F3A1*CORO**3+F4A1*CORO**4
      FE3= F1A3*CORO+F2A3*CORO**2+F3A3*CORO**3+F4A3*CORO**4
      FE11= F11+F1A11*CORO+F2A11*CORO**2
      FE33= F33+F1A33*CORO+F2A33*CORO**2
      FE13= F13+F1A13*CORO
      FE111= F111+F1A111*CORO+F2A111*CORO**2
      FE333= F333+F1A333*CORO+F2A333*CORO**2
      FE113= F113+F1A113*CORO+F2A113*CORO**2
      FE133= F133+F1A133*CORO+F2A133*CORO**2
      FE1111= F1111+FA1111*CORO
      FE3333= F3333+FA3333*CORO
      FE1113= F1113+FA1113*CORO
      FE1333= F1333+FA1333*CORO
      FE1133=       FA1133*CORO
      FE11111=F11111
      FE33333=F33333
      FE111111=F111111
      FE333333=F333333
      FE71    =F71
      FE73    =F73
      V   = V0 +  FE1*Y1+FE3*Y3
     1         +FE11*Y1**2+FE33*Y3**2+FE13*Y1*Y3
     2         +FE111*Y1**3+FE333*Y3**3+FE113*Y1**2*Y3
     3         +FE133*Y1*Y3**2
     4         +FE1111*Y1**4+FE3333*Y3**4+FE1113*Y1**3*Y3
     5         +FE1333*Y1*Y3**3+FE1133*Y1**2*Y3**2
     6         +FE11111*Y1**5+FE33333*Y3**5
     7         +FE111111*Y1**6+FE333333*Y3**6
     8         +FE71    *Y1**7+FE73    *Y3**7
C     modification by Choi & Light, J. Chem. Phys., 97, 7031 (1992).
      sqrt2=sqrt(2.0)
      xmup1=sqrt2/3.0+0.5
      xmum1=xmup1-x1
      term=2.0*xmum1*xmup1*q1*q2*cos(theta)
      r1=toang*sqrt((xmup1*q1)**2+(xmum1*q2)**2-term)
      r2=toang*sqrt((xmum1*q1)**2+(xmup1*q2)**2-term)
      rhh=sqrt(q1**2+q2**2-2.0*q1*q2*cos(theta))
      rbig=(r1+r2)/sqrt2
      rlit=(r1-r2)/sqrt2
 
      alpha=(x1-tanh(gammas*rbig**2))*(x1-tanh(gammaa*rlit**2))
      alpha1=beta1*alpha
      alpha2=beta2*alpha
      drhh=toang*(rhh-delta*rhh0)
      DOLEG=     (1.4500-THETA)
C     IF (THETA.LE.0.64  ) V=0.1E17
C     IF((DR.LE.-0.4).AND.(THETA.LE.1.1)) V=0.1E17
C     IF((DS.LE.-0.4).AND.(THETA.LE.1.1)) V=0.1E17
C     IF (DS.LE. 0.0  ) V=0.1E17
      v = v + c1*exp(-alpha1*drhh) + c2*exp(-alpha2*drhh)
 
C     Convert to Hartree
      V=V/CMTOAU
      RETURN
      END
