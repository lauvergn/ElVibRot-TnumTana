c
C================================================================
C    calc_Op : calculation of the potential and dipolar matrices
c    mat_V(nb_be,nb_be) and mat_dip(nb_be,nb_be,3)
c    nb_be : nb of elctronic surface
c    Q are the coordinates in active order or syl order
c    dipolar calculation if calc_dip = T
C================================================================
      SUBROUTINE calcN_op(mat_V,mat_imV,mat_ScalOp,nb_be,nb_ScalOp,
     *                   Q,nb_var,mole,
     *                   calc_ScalOp,pot_cplx)

      USE mod_system
      USE mod_Tnum
      IMPLICIT NONE

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

      integer           :: nb_be,nb_ScalOp,nb_var
      logical           :: calc_ScalOp,pot_cplx
      real (kind=Rkind) :: mat_V(nb_be,nb_be),mat_imV(nb_be,nb_be)
      real (kind=Rkind) :: mat_ScalOp(nb_be,nb_be,nb_ScalOp)
      real (kind=Rkind) :: Q(nb_var)

      real (kind=Rkind) :: murrell,im_pot0
      real (kind=Rkind) :: ScalOp(nb_ScalOp)



      IF (nb_be == 1 ) THEN
        mat_V(1,1) = murrell(Q)
        !write(6,*) 'Q,mat_V',Q,mat_V(1,1)
        IF (pot_cplx) mat_imV(1,1) = im_pot0(Q)
        IF (calc_ScalOp) THEN
          CALL sub_ScalarOp(ScalOp,nb_ScalOp,Q,mole)
          mat_ScalOp(1,1,:) = ScalOp(:)
        END IF
      ELSE
        STOP 'nb_be > 1'
      END IF

      END
c=================================================================
c
c     potentiel de Murrell pour HCN en coordonnees de Jacobi
c     Ref: J. N. Murrell, S. Carter and L. O. Halonene, i
c          J. Mol. Spectrosc. vo93 p307 1982
c          
c                       TT  C
c                   H_____ / gR=Q(3)
c                  pR=Q(2)/
c                        N
c    Q(1) = T = cos(TT)
c
c    le vrai potentiel est en fonction des 3 distances
c    R1 = R_CH
c    R2 = R_CN
c    R3 = R_HN
c=================================================================
      FUNCTION murrell(Q)
      USE mod_system
      IMPLICIT NONE
      real(kind=Rkind) :: murrell

      real(kind=Rkind) autoA,autoeV
c     parameter (autoA=.529178d0)
c     parameter (autoeV=27.21183d0)
      parameter (autoA=0.52917715_Rkind)
      parameter (autoeV=27.21183_Rkind)


      real(kind=Rkind) Q(3)
      real(kind=Rkind) gR,pR

      real(kind=Rkind) mB,mC,mH
      real(kind=Rkind) T,RN,RC,r1,r2,r3

      real(kind=Rkind) Z1,Z12,V1
      real(kind=Rkind) Z2,Z22,V2
      real(kind=Rkind) Z3,Z32,V3

      real(kind=Rkind) S1,S2,S3,S12,S22,S32,POLY
      real(kind=Rkind) E1,E2,E3,HY1,HY2,HY3,SWITCH,V123

      T  = Q(1)
      pR = Q(2)
      gR = Q(3)
      !write(out_unitp,*) Q

      mB=12._Rkind
      mC=14.003074_Rkind
      mH=1.007825_Rkind

c     T=cos(TT)
      RN=mC/(mB+mC)
      RC=ONE-RN
      r1=pR**2+(RN*gR)**2-TWO*pR*gR*RN*T
      R1=sqrt(r1)
      R2=gR
      R3=sqrt(pR**2+(RC*gR)**2+TWO*pR*gR*RC*T)
      !write(out_unitp,*) r1,r2,r3
      !flush(out_unitp)
      R1=R1*autoA
      R2=R2*autoA
      R3=R3*autoA

C.....CH
      Z1=R1-1.0823_Rkind
      Z12=Z1*Z1
      V1=-2.8521_Rkind*(ONE+5.5297_Rkind*Z1+
     *    8.7166_Rkind*Z12+5.3082_Rkind*Z1*Z12)*
     *    exp(-5.5297_Rkind*Z1)
C.....CN
      Z2=R2-1.1718_Rkind
      Z22=Z2*Z2 
      V2=-7.9282_Rkind*(ONE+5.2448_Rkind*Z2+
     *       7.3416_Rkind*Z22+4.9785_Rkind*Z2*Z22)*
     *    exp(-5.2448_Rkind*Z2)
C.....NH
      Z3=R3-1.0370_Rkind
      Z32=Z3*Z3
      V3=-3.9938_Rkind*(ONE+3.0704_Rkind*Z3)*exp(-3.0704_Rkind*Z3)

c     write(out_unitp,*) v1,v2,v3

C.....THREE BODY TERMS
      Z1=R1-1.9607_Rkind
      Z2=R2-2.2794_Rkind
      Z3=R3-1.8687_Rkind
      S1=0.4436_Rkind*Z1+0.6091_Rkind*Z2+0.6575_Rkind*Z3
      S2=-.8941_Rkind*Z1+0.2498_Rkind*Z2+0.3718_Rkind*Z3
      S3=0.0622_Rkind*Z1-0.7527_Rkind*Z2+0.6554_Rkind*Z3

      S12=S1*S1
      S22=S2*S2
      S32=S3*S3
      POLY=-3.0578_Rkind*(ONE+1.9076_Rkind*S1-0.5008_Rkind*S2-
     *      0.0149_Rkind*S3+0.6695_Rkind*S12-
     +      1.3535_Rkind*S22-1.0501_Rkind*S32+
     *      0.2698_Rkind*S1*S2-1.1120_Rkind*S1*S3+
     +      1.9310_Rkind*S2*S3-0.0877_Rkind*S1*S12+
     *      0.0044_Rkind*S2*S22+0.0700_Rkind*S3*S32+
     +      0.0898_Rkind*S12*S2-1.0186_Rkind*S1*S22-
     *      0.0911_Rkind*S12*S3+
     +      0.0017_Rkind*S1*S32+0.4567_Rkind*S22*S3-
     *      0.8840_Rkind*S2*S32+
     +      0.3333_Rkind*S1*S2*S3-0.0367_Rkind*S12*S12+
     *      0.4821_Rkind*S22*S22+
     +      0.2564_Rkind*S32*S32-0.0017_Rkind*S12*S1*S2-
     *      0.2278_Rkind*S12*S22-
     +      0.1287_Rkind*S1*S2*S22+0.1759_Rkind*S1*S12*S3-
     *      0.0399_Rkind*S12*S32-
     +      0.1447_Rkind*S1*S3*S32-0.3147_Rkind*S2*S22*S3+
     *      0.1233_Rkind*S22*S32+
     +      0.3161_Rkind*S2*S3*S32+0.0919_Rkind*S12*S2*S3-
     *      0.0954_Rkind*S1*S22*S3+
     +      0.1778_Rkind*S1*S2*S32-0.1892_Rkind*S22*S22*S2)


      E1=exp(3.9742_Rkind*Z1/TWO)
      E2=exp(4.3688_Rkind*Z2/TWO)
      E3=exp(1.5176_Rkind*Z3/TWO)
      HY1=(E1-ONE/E1)/(E1+ONE/E1)
      HY2=(E2-ONE/E2)/(E2+ONE/E2)
      HY3=(E3-ONE/E3)/(E3+ONE/E3)
      SWITCH=(ONE-HY1)*(ONE-HY2)*(ONE-HY3)
      V123=SWITCH*POLY

c     write(out_unitp,*) v123,V1+V2+V3+V123

      murrell=(V1+V2+V3+V123)/autoeV


      RETURN
      END
C================================================================
C    fonction pot_rest(x)
C================================================================
      FUNCTION pot_rest(Qact,Delta_Qact,nb_inact2n)
      USE mod_system
      IMPLICIT NONE
      real(kind=Rkind) :: pot_rest

       real(kind=Rkind) Qact(1)
       integer nb_inact2n
       real(kind=Rkind) Delta_Qact(nb_inact2n)

       pot_rest = ZERO

       END
C================================================================
C    fonction im_pot0(x)
C================================================================
      FUNCTION im_pot0(Q)
      USE mod_system
      IMPLICIT NONE
      real(kind=Rkind) :: im_pot0

      real(kind=Rkind) Q(3)

       real(kind=Rkind), parameter :: a=-0.01_Rkind,Q0=2.5_Rkind

       IF (Q(3) > Q0) THEN
         im_pot0 = a*(Q(3) - Q0)
       ELSE
         im_pot0 = ZERO
       END IF

       RETURN
       END
C================================================================
C    sub hessian
C================================================================
      SUBROUTINE sub_hessian (h)
      USE mod_system
      IMPLICIT NONE

       real(kind=Rkind) h

       h = ZERO


      END
c
C================================================================
C    fonction pot0(x) 1 D (avec x=cos(theta))
c    pour une tri atomique en jacobie
C================================================================
      SUBROUTINE sub_ScalarOp(ScalOp,nb_ScalOp,Q,mole)
      USE mod_system
      USE mod_Tnum
      IMPLICIT NONE

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

       integer nb_Q,nb_ScalOp
       parameter (nb_Q=3)
       real(kind=Rkind) Q(nb_Q)
       real(kind=Rkind) ScalOp(nb_ScalOp)
       integer :: i


       ScalOp(:) = ZERO
       DO i=1,nb_ScalOp
         ScalOp(i) = real(i,kind=Rkind) + 0.01_Rkind * Q(1) + 
     *                                   0.02_Rkind*(Q(2)+Q(3))
       END DO
       ScalOp(2) = ZERO
       ScalOp(3) = ZERO


       END
C================================================================
C    analytical gradient along the path
C================================================================
      SUBROUTINE d0d1d2_g(d0g,d1g,d2g,Qdyn,mole,deriv,num,step)
      USE mod_system
      USE mod_Tnum
      IMPLICIT NONE

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

      real (kind=Rkind) :: d0g(mole%nb_inact2n)
      real (kind=Rkind) :: d1g(mole%nb_inact2n,mole%nb_act1)
      real (kind=Rkind) :: 
     *              d2g(mole%nb_inact2n,mole%nb_act1,mole%nb_act1)

      real (kind=Rkind) :: Qdyn(mole%nb_var)
      real (kind=Rkind) :: step
      logical       :: deriv,num

      real (kind=Rkind) :: Qact(mole%nb_act1)

c----- for debuging ----------------------------------
c     logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
c---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING d0d1d2_g'
        write(out_unitp,*) 'nb_var',mole%nb_var
        write(out_unitp,*) 'nb_act1',mole%nb_act1
        write(out_unitp,*) 'nb_inact22,nb_inact21',
     *                   mole%nb_inact22,mole%nb_inact21
        write(out_unitp,*) 'nb_inact2n',mole%nb_inact2n
        write(out_unitp,*) 'deriv',deriv
      END IF

c---------------------------------------------------------------------
      Qact(1) = Qdyn(mole%liste_QactTOQsym(1))

      d0g(:)     = ZERO
      d1g(:,:)   = ZERO
      d2g(:,:,:) = ZERO

c---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'd0g at Qact:',Qact
        write(out_unitp,*) d0g(:)
        write(out_unitp,*) 'END d0d1d2_g'
      END IF
c---------------------------------------------------------------------

      END SUBROUTINE d0d1d2_g
C================================================================
C    analytical hessian along the path
C================================================================
      SUBROUTINE d0d1d2_h(d0h,d1h,d2h,
     *                     Qdyn,mole,deriv,num,step)
      USE mod_system
      USE mod_Tnum
      IMPLICIT NONE

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

       real (kind=Rkind) :: Qdyn(mole%nb_var)
       real (kind=Rkind) :: c_act


       real (kind=Rkind) :: d0h(mole%nb_inact2n,mole%nb_inact2n)
       real (kind=Rkind) :: d1h(mole%nb_inact2n,mole%nb_inact2n)
       real (kind=Rkind) :: d2h(mole%nb_inact2n,mole%nb_inact2n)

       real (kind=Rkind)  :: poly_legendre ! function
       character (len=14) :: nom
       logical :: exist
       integer :: iv,jv,i,j,kl
       logical :: deriv,num
       real (kind=Rkind) :: step

       integer, parameter      ::  max_points=200
       integer, parameter      ::  nb_inactb = 5
       real (kind=Rkind), save :: F(max_points,nb_inactb,nb_inactb)
       integer, save           :: nn(nb_inactb,nb_inactb)
       logical, save           :: begin = .TRUE.

c----- for debuging ----------------------------------
      logical, parameter :: debug = .FALSE.
c     logical, parameter :: debug = .TRUE.
c---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING d0d1d2_h'
        write(out_unitp,*) 'nb_var',mole%nb_var
        write(out_unitp,*) 'nb_act1',mole%nb_act1
        write(out_unitp,*) 'nb_inact22,nb_inact21',
     *            mole%nb_inact22,mole%nb_inact21
        write(out_unitp,*) 'nb_inact2n',mole%nb_inact2n
      END IF
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c---------------------------------------------------------------------
c     initialization (only once)
c$OMP CRITICAL (d0d1d2_h_CRIT)
      IF (begin) THEN

        IF (nb_inactb < mole%nb_inact2n ) THEN
          write(out_unitp,*) 'ERROR : nb_inactb is TO small',nb_inactb
          write(out_unitp,*) 'it should at least equal to',
     *                       mole%nb_inact2n
          STOP
        END IF

        DO i=1,mole%nb_inact2n
        DO j=i,mole%nb_inact2n

          iv = mole%liste_QactTOQsym(mole%nb_act1+i)
          jv = mole%liste_QactTOQsym(mole%nb_act1+j)
          IF (iv > jv) THEN
            nom=nom_ii('inter12h__',jv,iv)
          ELSE
            nom=nom_ii('inter12h__',iv,jv)
          END IF

          CALL read_para0d(F(1,i,j),nn(i,j),max_points,nom,exist)
          IF ( .NOT. exist ) THEN
            write(out_unitp,*) 'F(1,i,i) tq d0h = 1 ...'
            write(out_unitp,*) '... and F(1,i,j) tq d0h = 0'
            nn(i,j) = 1
            IF ( i .EQ. j ) THEN
              F(1,i,j) = ONE/poly_legendre(ONE,1,0)
            ELSE
              F(1,j,i) = ZERO
              F(1,i,j) = ZERO
            END IF
          END IF

        END DO
        END DO

        begin=.FALSE.
      END IF
c$OMP END CRITICAL (d0d1d2_h_CRIT)
c     END initialization
c---------------------------------------------------------------------
c---------------------------------------------------------------------

      c_act = Qdyn(mole%liste_QactTOQsym(1))
 
      IF (deriv) THEN
        write(out_unitp,*) 'ERROR in d0d1d2_h'
        write(out_unitp,*) '  deriv CANNOT be true!!'
        write(out_unitp,*) ' check the fortran source'
        STOP
      END IF

      DO i=1,mole%nb_inact2n
      DO j=i,mole%nb_inact2n
        d0h(i,j)=ZERO
        DO kl=1,nn(i,j)
           d0h(i,j) = d0h(i,j) + F(kl,i,j)*poly_legendre(c_act,kl,0)
        END DO
        d0h(j,i) = d0h(i,j)
c       write(out_unitp,*) 'd0h(i,j)',i,j,d0h(i,j)
      END DO
      END DO

c---------------------------------------------------------------------
      IF (debug) THEN       
       !write(out_unitp,*) 'Qact1',c_act
       !DO i=1,mole%nb_inact2n
       !DO j=i,mole%nb_inact2n
       !  write(out_unitp,*) 'F(.,i,j)',i,j,nn(i,j)
       !  write(out_unitp,*) (F(k,i,j),k=1,nn(i,j))
       !END DO
       !END DO
        write(out_unitp,*) 'd0h at c_act:',c_act
        CALL Write_Mat(d0h,6,4)
        write(out_unitp,*) 'END d0d1d2_h'
      END IF
c---------------------------------------------------------------------

      END SUBROUTINE d0d1d2_h
C================================================================
C    analytical derivative (dnQflex) calculation
c    for the variable iq
C================================================================
      SUBROUTINE calc_dnQflex(iq,dnQflex,Qact,nb_act,nderiv,it)
      USE mod_system
      USE mod_dnSVM      
      IMPLICIT NONE

       integer           :: iq,nb_act
       real (kind=Rkind) :: Qact(nb_act)
       integer           :: nderiv,it
       TYPE (Type_dnS)   :: dnQflex




       real (kind=Rkind) :: dc0,dc1,dc2,dc3
       real (kind=Rkind) :: c_act

       character (len=14) :: nom
       logical :: exist

       integer :: vi,kl

       integer, parameter      ::  max_points=200
       integer, parameter      ::  nb_inactb = 5
       real (kind=Rkind), save :: F(max_points,nb_inactb)
       integer, save           :: nn(nb_inactb)
       logical, save           :: begin = .TRUE.


c----- for debuging ----------------------------------
      character (len=*), parameter :: name_sub='dnQflex'
      logical, parameter :: debug=.FALSE.
c     logical, parameter :: debug=.TRUE.
c----- for debuging ----------------------------------


c---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_act',nb_act
        write(out_unitp,*) 'iq',iq
      END IF
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c      Qact value. Rq: only ONE active variable is possible
c---------------------------------------------------------------------
       IF (nb_act /= 1) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' the number of Active variable'
         write(out_unitp,*) ' should be 1. But nb_act =',nb_act
         STOP
       END IF

c---------------------------------------------------------------------
c---------------------------------------------------------------------


c---------------------------------------------------------------------
c---------------------------------------------------------------------
c      initialization (only once)
c$OMP    CRITICAL (dnQflex_CRIT)
       IF (begin) THEN
         write(out_unitp,*) ' INITIALIZATION of ',name_sub
         begin=.FALSE.
         nn(:) = 0
         F(:,:) = ZERO
         DO vi=2,3
           nom=nom_i('inter12___',vi)
           write(out_unitp,*) 'read file :',nom,vi

           CALL read_para0d(F(1,vi),nn(vi),max_points,nom,exist)
           IF ( .NOT. exist ) STOP

c          write(out_unitp,*) vi,(F(k,vi),k=1,nn(vi))
         END DO
         write(out_unitp,*) ' END INITIALIZATION of ',name_sub

       END IF
c$OMP    END CRITICAL (dnQflex_CRIT)
c     end  initialization
c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
       CALL sub_ZERO_TO_dnS(dnQflex)
       c_act = Qact(1)

       IF (iQ == 1) THEN
         dnQflex%d0 = c_act
         IF (nderiv >= 1) dnQflex%d1(1)     = ONE
       ELSE

       DO kl=1,nn(iQ)
         CALL d0d1d2d3poly_legendre(c_act,kl,dc0,dc1,dc2,dc3,nderiv)

         dnQflex%d0 = dnQflex%d0 + F(kl,iQ) * dc0

         IF (nderiv >= 1) 
     *     dnQflex%d1(1)     = dnQflex%d1(1)     + F(kl,iQ)*dc1

         IF (nderiv >= 2) 
     *     dnQflex%d2(1,1)   = dnQflex%d2(1,1)   + F(kl,iQ)*dc2

         IF (nderiv >= 3) 
     *     dnQflex%d3(1,1,1) = dnQflex%d3(1,1,1) + F(kl,iQ)*dc3

      END DO
      END IF

c---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'dnQflex : ',Qact
        CALL write_dnS(dnQflex,nderiv)
        write(out_unitp,*) 'END ',name_sub
      END IF
c---------------------------------------------------------------------

       END SUBROUTINE calc_dnQflex
C================================================================
C    analytical derivative (Qeq) calculation
c    for the variable i_Qdyn
C================================================================
      SUBROUTINE d0d1d2d3_Qeq(i_Qdyn,
     *                        d0req,d1req,d2req,d3req,
     *                        Qdyn,mole,nderiv)
      USE mod_system
      USE mod_Tnum
      IMPLICIT NONE

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

       integer :: i_Qdyn
       real (kind=Rkind) ::  Qdyn(mole%nb_var)

       integer :: nderiv

       real (kind=Rkind) ::  d0req
       real (kind=Rkind) ::  d1req(mole%nb_act1)
       real (kind=Rkind) ::  d2req(mole%nb_act1,mole%nb_act1)
       real (kind=Rkind) :: 
     *             d3req(mole%nb_act1,mole%nb_act1,mole%nb_act1)



       real (kind=Rkind) :: dc0,dc1,dc2,dc3
       real (kind=Rkind) :: c_act

       character (len=14) :: nom
       logical :: exist

       integer :: vi,kl,i_act,j_act,k_act

       integer, parameter      ::  max_points=200
       integer, parameter      ::  nb_inactb = 5
       real (kind=Rkind), save :: F(max_points,nb_inactb)
       integer, save           :: nn(nb_inactb)
       logical, save           :: begin = .TRUE.


c----- for debuging ----------------------------------
      character (len=*), parameter :: name_sub='d0d1d2d3_Qeq'
      logical, parameter :: debug=.FALSE.
c     logical, parameter :: debug=.TRUE.
c----- for debuging ----------------------------------

c---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_inact20,nb_act',
     *                     mole%nb_inact20,mole%nb_act
        write(out_unitp,*) 'nb_var',mole%nb_var
        write(out_unitp,*) 'i_Qdyn',i_Qdyn
      END IF
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c      Qact value. Rq: only ONE active variable is possible
c---------------------------------------------------------------------
       IF (mole%nb_act1 /= 1) THEN
         write(out_unitp,*) ' ERROR : d0d1d2d3_Qeq'
         write(out_unitp,*) ' the number of Active variable'
         write(out_unitp,*) ' should be 1. But nb_act1 =',mole%nb_act1
         STOP
       END IF

c---------------------------------------------------------------------
c---------------------------------------------------------------------


c---------------------------------------------------------------------
c---------------------------------------------------------------------
c      initialization the first time
c$OMP  CRITICAL (d0d1d2d3_Qeq_CRIT)
       IF (begin) THEN
         begin=.FALSE.

c        -------------------------------------------------------------
c         nb_inact (=nb_inact20+nb_inact21+nb_inact22) > nb_inactb ??
c        -------------------------------------------------------------
         IF (mole%nb_inact > nb_inactb) THEN
           write(out_unitp,*) ' ERROR : in ',name_sub
           write(out_unitp,*) 'nb_inact(',mole%nb_inact,
     *                    ')>nb_inactb(',nb_inactb,')'
           STOP
         END IF

         nn(:) = 0
         F(:,:) = ZERO
         DO vi=2,3
           nom=nom_i('inter12___',vi)
           write(out_unitp,*) 'read file :',nom,vi

           CALL read_para0d(F(1,vi),nn(vi),max_points,nom,exist)
           IF ( .NOT. exist ) STOP

c          write(out_unitp,*) vi,(F(k,vi),k=1,nn(vi))
         END DO
       END IF
c$OMP  END CRITICAL (d0d1d2d3_Qeq_CRIT)
c      end  initialisation
c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
       c_act = Qdyn(mole%liste_QactTOQsym(1))
c---------------------------------------------------------------------


c---------------------------------------------------------------------
       d0req = ZERO
       d1req(:) = ZERO
       d2req(:,:) = ZERO
       d3req(:,:,:) = ZERO
c---------------------------------------------------------------------

       i_act = mole%nb_act1
       j_act = mole%nb_act1
       k_act = mole%nb_act1


         DO kl=1,nn(i_Qdyn)
           CALL d0d1d2d3poly_legendre(c_act,kl,dc0,dc1,dc2,dc3,nderiv)

           d0req = d0req + F(kl,i_Qdyn) * dc0

           IF (nderiv .GE. 1)  d1req(i_act) =
     *                         d1req(i_act) + F(kl,i_Qdyn) * dc1

           IF (nderiv .GE. 2)  d2req(i_act,j_act) =
     *                         d2req(i_act,j_act) + F(kl,i_Qdyn)*dc2

           IF (nderiv .GE. 3)  d3req(i_act,j_act,k_act) =
     *                      d3req(i_act,j_act,k_act) + F(kl,i_Qdyn)*dc3

         END DO

c---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'd0req : ',c_act,d0req
        IF (nderiv > 0) write(out_unitp,*) 'd1req : ',c_act,d1req
        IF (nderiv > 1) write(out_unitp,*) 'd2req : ',c_act,d2req
        IF (nderiv > 2) write(out_unitp,*) 'd3req : ',c_act,d3req
        write(out_unitp,*) 'END d0d1d2d3_Qeq'
      END IF
c---------------------------------------------------------------------

      END SUBROUTINE d0d1d2d3_Qeq
