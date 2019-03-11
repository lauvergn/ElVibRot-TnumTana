c
C================================================================
C    calc_Op : calculation of the potential and scalar operator matrices
c    mat_V(nb_be,nb_be) and Mat_Scal(nb_be,nb_be,nb_ScalOp)
c    nb_be : nb of elctronic surfaces
c    Qop are the coordinates in active order or dyn order
C================================================================
      SUBROUTINE calcN_op(mat_V,mat_imV,mat_ScalOp,nb_be,nb_ScalOp,
     *                    Qop,nb_QOp,mole,calc_ScalOp,pot_cplx)

      USE mod_system
      USE mod_Tnum
      IMPLICIT NONE

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

      integer           :: nb_be,nb_ScalOp,nb_QOp
      logical           :: calc_ScalOp,pot_cplx
      real (kind=Rkind) :: mat_V(nb_be,nb_be),mat_imV(nb_be,nb_be)
      real (kind=Rkind) :: mat_ScalOp(nb_be,nb_be,nb_ScalOp)
      real (kind=Rkind) :: Qop(nb_QOp)

      integer           :: nb_QOp_loc=3


      IF (nb_be == 1 ) THEN
        write(out_unitp,*) ' ERROR sub_system for 2 PES'
        STOP
      END IF

      CALL Mat_pot0(mat_V,Qop(1:nb_QOp_loc),nb_be,nb_QOp_loc) 
      mat_V(1,1) = min(mat_V(1,1),ONE)
      mat_V(2,2) = min(mat_V(2,2),ONE)
      !mat_V(1,2) = ZERO
      !mat_V(2,1) = ZERO

      IF (pot_cplx) THEN
        CALL Mat_im_pot0(mat_imV,Qop(1:nb_QOp_loc),nb_be,nb_QOp_loc)
      END IF

      IF (calc_ScalOp) THEN
        CALL Mat_ScalarOp(mat_ScalOp,Qop(1:nb_QOp_loc),mole,
     *                    nb_be,nb_ScalOp,nb_QOp_loc)
      END IF

      END SUBROUTINE calcN_op
C================================================================
C     pot0(x) 1 D 2 surfaces
C================================================================
      SUBROUTINE Mat_pot0(mat_V,Qop,nb_be,nb_QOp)
      USE mod_system
      USE mod_Constant
      IMPLICIT NONE

      integer           :: nb_be,nb_QOp
      real (kind=Rkind) :: Qop(nb_QOp)
      real (kind=Rkind) :: mat_V(nb_be,nb_be)
!declaration of the functions (2 diabatic states and 1 coupling term)
      real (kind=Rkind) :: Hdir2D
      real (kind=Rkind) :: Hct2D
      real (kind=Rkind) :: Hdircorr
      real (kind=Rkind) :: Hctcorr
      real (kind=Rkind) :: Hdir
      real (kind=Rkind) :: Hct
      real (kind=Rkind) :: Hcp
      real (kind=Rkind) :: HOOP, Tors, BLA
!declaration of the parameters: 
      real (kind=Rkind), parameter :: d1 = 2571.28_Rkind
      real (kind=Rkind), parameter :: d2 = 50.5629_Rkind
      real (kind=Rkind), parameter :: d3 = 3.7302_Rkind
      real (kind=Rkind), parameter :: d4 = 594.966_Rkind
      real (kind=Rkind), parameter :: c1 = 437.068_Rkind
      real (kind=Rkind), parameter :: c2 = 16.7317_Rkind
      real (kind=Rkind), parameter :: c3 = 7.35468_Rkind
      real (kind=Rkind), parameter :: c4 = 88.517_Rkind
      real (kind=Rkind), parameter :: c5 = 5.95221_Rkind
      real (kind=Rkind), parameter :: hd1 = 126.216_Rkind
      real (kind=Rkind), parameter :: hd2 = 75.0032_Rkind
      real (kind=Rkind), parameter :: hc1 = 86.4639_Rkind
      real (kind=Rkind), parameter :: hc2 = 43.4614_Rkind
      real (kind=Rkind), parameter :: cp1 = 13.6991_Rkind
      real (kind=Rkind), parameter :: cp2 = 1.13469_Rkind
       
        BLA = Qop(1) * get_Conv_au_TO_unit('L','Angs')
        Tors = Qop(2)
        HOOP = Qop(3) + TWO * Pi
      
       Hdir2D = SIN(Tors)**2 * (d1 * BLA**2 + d2) +
     *     d3 * COS(Tors / 2_Rkind)**2 + d4 * (BLA - 0.091_Rkind)**2
       Hct2D = (1_Rkind + c5 * SIN(Tors)**2) *
     *     (c1 * BLA**2 + c2 * BLA + c3) + c4 * COS(Tors)**2
       Hdircorr = hd1 * SIN(HOOP / 4_Rkind)**2  -
     *     hd2 * SIN(HOOP / 4_Rkind) * SIN(Tors * 2_Rkind)
       Hctcorr = hc1 * SIN(HOOP / 4_Rkind)**2  +
     *     hc2 * SIN(HOOP / 4_Rkind) * SIN(Tors * 2_Rkind)

        mat_V(1,1) = Hdir2D + Hdircorr 
        mat_V(2,2) = Hct2D + Hctcorr
        mat_V(1,2) = (1_Rkind + cp2 * SIN(Tors)**2) *
     *     cp1 * SIN((-HOOP/2_Rkind + Tors) * 2_Rkind)
        mat_V(2,1) = mat_V(1,2)

       
        mat_V(:,:) = mat_V(:,:) / 627.509_Rkind

        !write(6,*) 'mat_V',BLA,Tors,HOOP,mat_V(1,1)

      END SUBROUTINE Mat_pot0
      
C================================================================
C    fonction im_pot0(x) imaginary part of pot0
C================================================================
      SUBROUTINE Mat_im_pot0(mat_imV,Qop,nb_be,nb_QOp)
      USE mod_system
      IMPLICIT NONE

      integer           :: nb_be,nb_QOp
      real (kind=Rkind) :: mat_imV(nb_be,nb_be)
      real (kind=Rkind) :: Qop(nb_QOp)

      integer           :: i
      real (kind=Rkind) :: im_pot0 ! imaginary function

      mat_imV = ZERO

      DO i=1,nb_be
        mat_imV(i,i) = im_pot0(Qop)
      END DO

      END SUBROUTINE Mat_im_pot0
      FUNCTION im_pot0(Qop)
      USE mod_system
      IMPLICIT NONE
      real (kind=Rkind) :: im_pot0


       real (kind=Rkind) :: Qop(1)
       real (kind=Rkind) :: z
       real (kind=Rkind), parameter :: a=ONETENTH, Q0=0.9_Rkind

       z = ZERO
       IF (Qop(1) > Q0) z = -a * (Qop(1) - Q0)**2

       im_pot0 = z

       END FUNCTION im_pot0
C================================================================
C    fonction pot_rest(x)
C================================================================
      FUNCTION pot_rest(Qact,Delta_Qact,nb_inact2n)
      USE mod_system
      IMPLICIT NONE
      real (kind=Rkind) :: pot_rest


       real (kind=Rkind) :: Qact(1)
       integer :: nb_inact2n
       real (kind=Rkind) :: Delta_Qact(nb_inact2n)

       pot_rest = ZERO

       END FUNCTION pot_rest
C================================================================
C    sub hessian
C================================================================
      SUBROUTINE sub_hessian (h)
      USE mod_system
      IMPLICIT NONE

       real (kind=Rkind) :: h

       h = ZERO

       END SUBROUTINE sub_hessian
c
C================================================================
C    fonction pot0(x) 1 D (avec x=cos(theta))
c    pour une tri atomique en jacobie
C================================================================
      SUBROUTINE Mat_ScalarOp(mat_ScalOp,Qop,mole,
     *                        nb_be,nb_ScalOp,nb_QOp)
      USE mod_system
      USE mod_Tnum
      IMPLICIT NONE

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole


      integer           :: nb_be,nb_ScalOp,nb_QOp
      real (kind=Rkind) :: mat_ScalOp(nb_be,nb_be,nb_ScalOp)
      real (kind=Rkind) :: Qop(nb_QOp)


       mat_ScalOp(:,:,:) = ZERO

       mat_ScalOp(1,2,1) = ONE
       mat_ScalOp(2,1,1) = ONE

       END SUBROUTINE Mat_ScalarOp
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
       integer :: iv,jv,i,j,kl,k
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
        write(out_unitp,*) 'Qact1',c_act
        DO i=1,mole%nb_inact2n
        DO j=i,mole%nb_inact2n
          write(out_unitp,*) 'F(.,i,j)',i,j,nn(i,j)
          write(out_unitp,*) (F(k,i,j),k=1,nn(i,j))
        END DO
        END DO
        write(out_unitp,*) 'd0h at c_act:',c_act
        CALL Write_Mat(d0h,6,4)
        write(out_unitp,*) 'END d0d1d2_h'
      END IF
c---------------------------------------------------------------------

      END SUBROUTINE d0d1d2_h
C================================================================
C    analytical derivative (dnQflex : Qflex Qflex' Qflex" Qflex'") calculation
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

       integer :: vi,kl,k

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

c---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'dnQflex : ',Qact
        CALL write_dnS(dnQflex,nderiv)
        write(out_unitp,*) 'END ',name_sub
      END IF
c---------------------------------------------------------------------

       END SUBROUTINE calc_dnQflex
C================================================================
C    analytical derivative (Qeq Qeq' Qeq" Qeq'") calculation
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

       integer :: vi,kl,k,i_act,j_act,k_act

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
        write(out_unitp,*) 'd1req : ',c_act,d1req
        write(out_unitp,*) 'd2req : ',c_act,d2req
        write(out_unitp,*) 'd3req : ',c_act,d3req
        write(out_unitp,*) 'END d0d1d2d3_Qeq'
      END IF
c---------------------------------------------------------------------

      END SUBROUTINE d0d1d2d3_Qeq
