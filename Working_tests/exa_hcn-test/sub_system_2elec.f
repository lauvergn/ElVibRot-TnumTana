c
C================================================================
C    calc_Op : calculation of the potential and scalar operator matrices
c    mat_V(nb_be,nb_be) and mat_dip(nb_be,nb_be,nb_ScalOp)
c    nb_be : nb of elctronic surfaces
c    Q are the coordinates in active order or dyn order
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

      real (kind=Rkind) :: pot0,im_pot0
      real (kind=Rkind) :: dip(3)

      IF (nb_be == 1 ) THEN
        write(out_unitp,*) ' ERROR sub_system for 2 PES'
        STOP
      ELSE
        CALL Mat_pot0(mat_V,Q) 
        IF (pot_cplx) CALL Mat_im_pot0(mat_imV,Q)
        IF (calc_ScalOp) THEN
          IF (nb_ScalOp /= 3) THEN
            write(out_unitp,*) ' ERROR in calcN_op'
            write(out_unitp,*) ' nb_ScalOp /= 3',nb_ScalOp
            STOP
          END IF
          CALL sub_dipole(mat_ScalOp,Q,mole)
        END IF
      END IF

      END
C================================================================
C     pot0(x) 1 D 2 surfaces
C================================================================
      SUBROUTINE Mat_pot0(mat_V,Q)
      USE mod_system
      IMPLICIT NONE

       real (kind=Rkind) :: Q(3)
       real (kind=Rkind) :: mat_V(2,2)

       real (kind=Rkind), parameter :: k11  = 0.1d0
       real (kind=Rkind), parameter :: k22  = 0.5d0
       real (kind=Rkind), parameter :: k33  = 0.5d0
       real (kind=Rkind), parameter :: q1e1   = 0.80d0
       real (kind=Rkind), parameter :: q1e2   = 0.80d0
       real (kind=Rkind), parameter :: q2e1   = 3.20d0
       real (kind=Rkind), parameter :: q2e2   = 3.00d0
       real (kind=Rkind), parameter :: q3e1   = 2.20d0
       real (kind=Rkind), parameter :: q3e2   = 2.00d0
       real (kind=Rkind), parameter :: e11  = 0.00d0
       real (kind=Rkind), parameter :: e22  = 0.03d0
       real (kind=Rkind), parameter :: h12  = 5.0d-3

       mat_V(1,1) = e11 + 0.5d0 * 
     *   (k11*(Q(1)-q1e1)**2 + k22*(Q(2)-q2e1)**2 + k33*(Q(3)-q3e1)**2)
       mat_V(2,2) = e22 + 0.5d0 * 
     *   (k11*(Q(1)-q1e2)**2 + k22*(Q(2)-q2e2)**2 + k33*(Q(3)-q3e2)**2)
       mat_V(1,2) = h12
       mat_V(2,1) = h12

       RETURN
       END
C================================================================
C    fonction im_pot0(x) imaginary part of pot0
C================================================================
      SUBROUTINE Mat_im_pot0(mat_imV,Q)
      USE mod_system
      IMPLICIT NONE

      real(kind=Rkind) :: Q(3)
      real(kind=Rkind) :: mat_imV(2,2)
      real(kind=Rkind) :: im_pot0 ! imaginary function

      mat_imV = 0.0
      mat_imV(1,1) = im_pot0(Q)
      mat_imV(2,2) = im_pot0(Q)

       END
C================================================================
C    fonction pot0(x) 1 D (avec x=cos(theta))
c    pour une tri atomique en jacobie
C================================================================
      FUNCTION pot0(Qsym0)
      USE mod_system
      IMPLICIT NONE
      real (kind=Rkind) :: pot0

       real (kind=Rkind) :: Qsym0(1)
       integer :: nn
       real (kind=Rkind) :: z
       real (kind=Rkind) :: poly_legendre

       character*14 nom
       logical :: exist

       integer :: max_points
       parameter (max_points=200)
       real (kind=Rkind) :: F(max_points)

       integer :: kl,i_type_act,i_sym_act
       real (kind=Rkind) ::  c_act,Qact1

       logical :: begin
       data begin/.true./
       SAVE begin,F,nn,i_type_act

c---------------------------------------------------------------
c      initialisation la premiere fois
       IF (begin) THEN
         nom='inter12-ene'
         CALL read_para0d(F,nn,max_points,nom,exist)
         IF ( .NOT. exist) STOP

         begin=.FALSE.
       END IF
c fin     initialisation la premiere fois
c---------------------------------------------------------------

       c_act = Qsym0(1)
       z = ZERO
       DO kl=1,nn
         z = z + F(kl) * poly_legendre(c_act,kl,0)
       END DO

       pot0 = z

c      write(out_unitp,*) Qsym0(1),c_act,z

       END
C================================================================
C    impaginary part of the potential (CAP)
C================================================================
      FUNCTION im_pot0(Q)
      IMPLICIT NONE
      real(kind=8) :: im_pot0

      real(kind=8) Q(3) ! coordinates in the Qsym order (cos(a), R, r) since pot_act=f

c      parameter of cap
       real(kind=8), parameter :: a=-0.71,Q0=2.0

       IF (Q(3) > Q0) THEN
         im_pot0 = a*(Q(3) - Q0)**2
       ELSE
         im_pot0 = 0.d0
       END IF
       
       RETURN
       END
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

       END
C================================================================
C    subroutine calculant le gradient
C================================================================
      SUBROUTINE d0d1d2_g(d0g,d1g,d2g,Qsym0,mole,deriv,num,step)
      USE mod_system
      USE mod_Tnum
      IMPLICIT NONE

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

      real (kind=Rkind) :: d0g(mole%nb_inact2n)
      real (kind=Rkind) :: d1g(mole%nb_inact2n,mole%nb_act1)
      real (kind=Rkind) :: 
     *              d2g(mole%nb_inact2n,mole%nb_act1,mole%nb_act1)

      real (kind=Rkind) :: Qsym0(mole%nb_var)
      real (kind=Rkind) :: step
      logical       :: deriv,num




      real (kind=Rkind) :: Qact(2),c2,s2
      real (kind=Rkind) :: d0f,d1f,d2f,d3f
      integer       :: ix_inact2n,iz_inact2n

c----- for debuging ----------------------------------
c     logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
c---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unitp,*)
      write(out_unitp,*) 'BEGINNING d0d1d2_g'
      write(out_unitp,*) 'nb_var',mole%nb_var
      write(out_unitp,*) 'nb_act1',mole%nb_act1
      write(out_unitp,*) 'nb_inact2n',mole%nb_inact2n
      write(out_unitp,*) 'deriv',deriv
      END IF

c---------------------------------------------------------------------
       Qact(1) = Qsym0(mole%liste_QactTOQsym(1))

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

       RETURN
       END
C================================================================
C    sub hessian
C================================================================
      SUBROUTINE sub_hessian (h)
      USE mod_system
      IMPLICIT NONE

       real (kind=Rkind) :: h

       h = ZERO


       END
      SUBROUTINE h0_sym(h)
      USE mod_system
      IMPLICIT NONE

       real (kind=Rkind) :: h

       RETURN
       END
C================================================================
C    subroutine calculant la matrice hessienne
C    en fonction de x=cos(theta)
C================================================================
       SUBROUTINE d0d1d2_h(d0h,d1h,d2h,
     *                     Qsym0,mole,deriv,num,step)
      USE mod_system
      USE mod_Tnum
      IMPLICIT NONE

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

       integer :: i_qact
       real (kind=Rkind) ::  Qact1,Qsym0(mole%nb_var)

       real (kind=Rkind) :: c_act,s_act,ss


       real (kind=Rkind) :: d0h(mole%nb_inact2n,mole%nb_inact2n)
       real (kind=Rkind) :: d1h(mole%nb_inact2n,mole%nb_inact2n)
       real (kind=Rkind) :: d2h(mole%nb_inact2n,mole%nb_inact2n)
       real (kind=Rkind) :: d0,d1,d2
       real (kind=Rkind) :: step

       real (kind=Rkind) :: poly_legendre
       character*14 nom_ii,nom
       logical :: exist

       integer ::  max_points,nb_inactb
       parameter (max_points=200)
       parameter (nb_inactb=5)
       real (kind=Rkind) :: F(max_points,nb_inactb,nb_inactb)
       integer :: nn(nb_inactb,nb_inactb)


       integer :: jv,i_sym_act,kl,i,j,k,iv,i_type_act

       logical :: deriv,num,begin
       data begin/.true./
       SAVE begin,F,nn,i_type_act

c----- for debuging ----------------------------------
      logical debug
      parameter (debug=.FALSE.)
c     parameter (debug=.TRUE.)
c---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unitp,*)
      write(out_unitp,*) 'BEGINNING d0d1d2_h'
      write(out_unitp,*) 'nb_var',mole%nb_var
      write(out_unitp,*) 'nb_act1',mole%nb_act1
      write(out_unitp,*) 'nb_inact2n',mole%nb_inact2n
      END IF
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c---------------------------------------------------------------------
c      initialisation la premiere fois sauf si direct=.true.
       IF (begin) THEN

         IF (nb_inactb .LT. mole%nb_inact2n ) THEN
           write(out_unitp,*) 'ERROR : nb_inactb is TO small',nb_inactb
           write(out_unitp,*) 'it should at least equal to',
     *              mole%nb_inact2n
           STOP
         END IF

         DO i=1,mole%nb_inact2n
         DO j=i,mole%nb_inact2n

           iv = mole%liste_QactTOQsym(mole%nb_act1+i)
           jv = mole%liste_QactTOQsym(mole%nb_act1+j)
           IF (iv .GT. jv) THEN
             nom=nom_ii('inter12h__',jv,iv)
           ELSE
             nom=nom_ii('inter12h__',iv,jv)
           END IF

           CALL read_para0d(F(1,i,j),nn(i,j),max_points,nom,exist)
           IF ( .NOT. exist ) THEN
             write(out_unitp,*) 'F(1,i,i) tq d0h = 1 '
             write(out_unitp,*) ' and F(1,i,j) tq d0h = 0'
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
c fin     initialisation
c---------------------------------------------------------------------
c---------------------------------------------------------------------

       c_act = Qsym0(1)
       IF (deriv) THEN

         DO i=1,mole%nb_inact2n
           DO j=i,mole%nb_inact2n
             d0h(i,j)=ZERO
             d1h(i,j)=ZERO
             d2h(i,j)=ZERO
             DO kl=1,nn(i,j)
               CALL d0d1d2poly_legendre(c_act,kl,d0,d1,d2,num,step)
               d0h(i,j) = d0h(i,j) + F(kl,i,j)*d0
               d1h(i,j) = d1h(i,j) + F(kl,i,j)*d1
               d2h(i,j) = d2h(i,j) + F(kl,i,j)*d2
             END DO
             d0h(j,i) = d0h(i,j)
             d1h(j,i) = d1h(i,j)
             d2h(j,i) = d2h(i,j)
c            write(out_unitp,*) 'd0h(i,j)',i,j,d0h(i,j)
           END DO
         END DO
       ELSE

         DO i=1,mole%nb_inact2n
           DO j=i,mole%nb_inact2n
             d0h(i,j)=ZERO
             DO kl=1,nn(i,j)
               d0 = poly_legendre(c_act,kl,0)
               d0h(i,j) = d0h(i,j) + F(kl,i,j)*d0
             END DO
             d0h(j,i) = d0h(i,j)
c            write(out_unitp,*) 'd0h(i,j)',i,j,d0h(i,j)
           END DO
         END DO
       END IF

c---------------------------------------------------------------------
       IF (debug) THEN       
         write(out_unitp,*) 'i_qact i_sym_act Qact1',
     *                 i_qact,i_sym_act,Qact1
         DO i=1,mole%nb_inact2n
         DO j=i,mole%nb_inact2n
           write(out_unitp,*) 'F(.,i,j)',i,j,nn(i,j)
           write(out_unitp,*) (F(k,i,j),k=1,nn(i,j))
         END DO
         END DO
         write(out_unitp,*) 'd0h at c_act:',c_act
         CALL ecriture(d0h,mole%nb_inact2n,
     *                 mole%nb_inact2n,4,.TRUE.,mole%nb_inact2n)
         write(out_unitp,*) 'END d0d1d2_h'
       END IF
c---------------------------------------------------------------------

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

       integer :: it,iq,nb_act
       real (kind=Rkind) ::  Qact(nb_act)

       integer :: nderiv

       TYPE (Type_dnS)   :: dnQflex


       integer :: vi


       real (kind=Rkind) :: dc0,dc1,dc2,dc3
       real (kind=Rkind) :: c_act,s_act,ss,sc,sss

       character*14 nom_i,nom
       logical :: exist

       integer, parameter ::  nb_inactb=10,max_points=200
       real (kind=Rkind), save :: F(max_points,nb_inactb)
       integer, save :: nn(nb_inactb)
       logical, save :: begin = .true.

       integer :: k_act,i_act,j_act,k,kl
       integer :: i_sym_act,i_type_act,i_qsym_inact

c----- function --------------------------------------
       real (kind=Rkind) :: poly_legendre
c----- function --------------------------------------


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
c      initialisation the first time
       IF (begin) THEN
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

       END IF
c end  initialisation the first time
c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
       dnQflex%d0 = ZERO
       IF (nderiv >= 1)  dnQflex%d1(:) = ZERO
       IF (nderiv >= 2)  dnQflex%d2(:,:) = ZERO
       IF (nderiv >= 3)  dnQflex%d3(:,:,:) = ZERO
c---------------------------------------------------------------------

         c_act = Qact(1)
         DO kl=1,nn(iQ)
           CALL d0d1d2d3poly_legendre(c_act,kl,dc0,dc1,dc2,dc3,nderiv)

           dnQflex%d0 = dnQflex%d0 + F(kl,iQ) * dc0

       IF (nderiv >= 1) dnQflex%d1(1) = dnQflex%d1(1) + F(kl,iQ) * dc1
       IF (nderiv >= 2) dnQflex%d2(1,1) = dnQflex%d2(1,1) + F(kl,iQ)*dc2
       IF (nderiv >= 3) dnQflex%d3(1,1,1)=dnQflex%d3(1,1,1)+F(kl,iQ)*dc3

         END DO

c---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) ' F(.,iQ) iQ',iQ
        write(out_unitp,*) (F(k,iQ),k=1,nn(iQ))
        write(out_unitp,*) 'dnQflex : ',Qact
        CALL write_dnS(dnQflex,nderiv)
        write(out_unitp,*) 'END ',name_sub
      END IF
c---------------------------------------------------------------------

       RETURN
       END
C================================================================
C    analytical derivative (Qeq Qeq' Qeq" Qeq'") calculation
c    for the variable i_qsym
C================================================================
      SUBROUTINE d0d1d2d3_Qeq(i_qsym,
     *                        d0req,d1req,d2req,d3req,
     *                        Qsym0,mole,nderiv)
      USE mod_system
      USE mod_Tnum
      IMPLICIT NONE

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

       integer :: i_qact,i_qsym
       real (kind=Rkind) ::  Qact1,Qsym0(mole%nb_var)

       integer :: nderiv

       real (kind=Rkind) ::  d0req
       real (kind=Rkind) ::  d1req(mole%nb_act1)
       real (kind=Rkind) ::  d2req(mole%nb_act1,mole%nb_act1)
       real (kind=Rkind) :: 
     *             d3req(mole%nb_act1,mole%nb_act1,mole%nb_act1)


       integer :: vi


       real (kind=Rkind) :: dc0,dc1,dc2,dc3
       real (kind=Rkind) :: c_act,s_act,ss,sc,sss

       character*14 nom_i,nom
       logical :: exist

       integer  ::  max_points,nb_inactb
       parameter (max_points=200)
       parameter (nb_inactb=10)
       real (kind=Rkind) :: F(max_points,nb_inactb)
       integer :: nn(nb_inactb)

       integer :: k_act,i_act,j_act,k,kl
       integer :: i_sym_act,i_type_act,i_qsym_inact

c----- function --------------------------------------
       real (kind=Rkind) :: poly_legendre
c----- function --------------------------------------

c----- for debuging ----------------------------------
      logical :: debug
      parameter (debug=.FALSE.)
c     parameter (debug=.TRUE.)
c----- for debuging ----------------------------------

       logical :: begin
       data begin/.true./
       SAVE begin,F,nn,i_type_act

c---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING d0d1d2d3_Qeq'
        write(out_unitp,*) 'nb_inact20,nb_act',
     *             mole%nb_inact20,mole%nb_act
        write(out_unitp,*) 'nb_var',mole%nb_var
        write(out_unitp,*) 'i_qsym',i_qsym
      END IF
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c      Qact value. Rq: only ONE active variable is possible
c---------------------------------------------------------------------
       IF (mole%nb_act1 .NE. 1) THEN
         write(out_unitp,*) ' ERROR : d0d1d2d3_Qeq'
         write(out_unitp,*) ' the number of Active variable'
         write(out_unitp,*) ' should be 1. But nb_act1 =',mole%nb_act1
         STOP
       END IF

c---------------------------------------------------------------------
c---------------------------------------------------------------------


c---------------------------------------------------------------------
c---------------------------------------------------------------------
c      initialisation the first time
       IF (begin) THEN
         begin=.FALSE.

c        -------------------------------------------------------------
c         nb_inact (=nb_inact20+nb_inact21+nb_inact22) > nb_inactb ??
c        -------------------------------------------------------------
         IF (mole%nb_inact .GT. nb_inactb) THEN
           write(out_unitp,*) ' ERROR : in d0d1d2d3_Qeq '
           write(out_unitp,*) 'nb_inact(',mole%nb_inact,
     *                    ')>nb_inactb(',nb_inactb,')'
           STOP
         END IF

         write(out_unitp,*) 'liste_QactTOQsym',mole%liste_QactTOQsym
         DO vi=1,mole%nb_inact

           i_qsym_inact = mole%liste_QactTOQsym(mole%nb_act1+vi)
           nom=nom_i('inter12___',i_qsym_inact)
           write(out_unitp,*) 'read file :',nom,i_qsym_inact

           CALL read_para0d(F(1,vi),nn(vi),max_points,nom,exist)
           IF ( .NOT. exist ) THEN
             write(out_unitp,*) 'F(1) tq d0req =',Qsym0(i_qsym_inact)
             nn(vi) = 1
             F(1,vi) = Qsym0(i_qsym_inact)/poly_legendre(ONE,1,0)
           END IF

c          write(out_unitp,*) vi,(F(k,vi),k=1,nn(vi))
         END DO

       END IF
c end  initialisation the first time
c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
       i_sym_act = mole%liste_QactTOQsym(1)
       Qact1  = Qsym0(i_sym_act)

c      write(out_unitp,*) 'i_qact i_sym_act Qact1',i_qact,i_sym_act,Qact1

       vi =  mole%liste_QsymTOQact(i_qsym)-mole%nb_act1
c      write(out_unitp,*) 'vi,i_qsym',vi,i_qsym
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


         c_act = Qact1
         DO kl=1,nn(vi)
           CALL d0d1d2d3poly_legendre(c_act,kl,dc0,dc1,dc2,dc3,nderiv)

           d0req = d0req + F(kl,vi) * dc0

           IF (nderiv .GE. 1)  d1req(i_act) =
     *                         d1req(i_act) + F(kl,vi) * dc1

           IF (nderiv .GE. 2)  d2req(i_act,j_act) =
     *                         d2req(i_act,j_act) + F(kl,vi)*dc2

           IF (nderiv .GE. 3)  d3req(i_act,j_act,k_act) =
     *                         d3req(i_act,j_act,k_act) + F(kl,vi)*dc3

         END DO

c---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) ' F(.,vi) vi',vi
        write(out_unitp,*) (F(k,vi),k=1,nn(vi))
        write(out_unitp,*) 'd0req : ',Qact1,d0req
        write(out_unitp,*) 'd1req : ',Qact1,d1req
        write(out_unitp,*) 'd2req : ',Qact1,d2req
        write(out_unitp,*) 'd3req : ',Qact1,d3req
        write(out_unitp,*) 'END d0d1d2d3_Qeq'
      END IF
c---------------------------------------------------------------------

       RETURN
       END
c
C================================================================
C    fonction pot0(x) 1 D (avec x=cos(theta))
c    pour une tri atomique en jacobie
C================================================================
      SUBROUTINE sub_dipole(Mat_dip,Qsym0,mole)
      USE mod_system
      USE mod_Tnum
      IMPLICIT NONE

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

       integer, parameter :: nb_Q=3
       real (kind=Rkind) :: Qsym0(nb_Q)
       real (kind=Rkind) :: Mat_dip(2,2,3)

       Mat_dip(:,:,:) = 0.d0

       Mat_dip(1,2,1) = 1.d0
       Mat_dip(2,1,1) = 1.d0


       END
