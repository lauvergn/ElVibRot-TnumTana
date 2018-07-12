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
      real (kind=Rkind) :: ScalOp(nb_ScalOp)


      IF (nb_be == 1 ) THEN
        mat_V(1,1) = pot0(Q)
        IF (pot_cplx) mat_imV(1,1) = im_pot0(Q)
        IF (calc_ScalOp) THEN
          CALL sub_ScalarOp(ScalOp,nb_ScalOp,Q,mole)
          mat_ScalOp(1,1,:) = ScalOp(:)
        END IF
      END IF

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

       real (kind=Rkind) :: z ! pot0
       real (kind=Rkind) :: poly_legendre ! function

       character*14 nom
       logical :: exist
       integer :: kl
       real (kind=Rkind) ::  c_act


       integer, save           :: nn
       integer, parameter      :: max_points = 200
       real (kind=Rkind), save :: F(max_points)
       logical, save           :: begin = .TRUE.

c---------------------------------------------------------------
c      initialization (only once)
c$OMP    CRITICAL (pot0_CRIT)
       IF (begin) THEN
c$       write(out_unitp,*) "F def thread",omp_get_thread_num()
         nom='inter12-ene'
         CALL read_para0d(F,nn,max_points,nom,exist)
         IF ( .NOT. exist) STOP

         begin=.FALSE.
       END IF
c$OMP    END CRITICAL (pot0_CRIT)
c      END initialization
c---------------------------------------------------------------

       c_act = Qsym0(1)
       z = ZERO
       DO kl=1,nn
         z = z + F(kl) * poly_legendre(c_act,kl,0)
       END DO

       pot0 = z

c      write(out_unitp,*) 'pot0',Qsym0(1),c_act,z

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
C    fonction im_pot0(x)
C================================================================
      FUNCTION im_pot0(Qsym0)
      USE mod_system
      IMPLICIT NONE
      real (kind=Rkind) :: im_pot0


       real (kind=Rkind) :: Qsym0(1)
       real (kind=Rkind) :: z
       real (kind=Rkind), parameter :: a=ONETENTH, Q0=0.9_Rkind

       z = ZERO
       IF (Qsym0(1) > Q0) z = -a * (Qsym0(1) - Q0)**2

       im_pot0 = z

       RETURN
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




      real (kind=Rkind) :: Qact(2)

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

       real (kind=Rkind) ::  Qsym0(mole%nb_var)
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
c      initialization (only once)
c$OMP    CRITICAL (d0d1d2_h_CRIT)
       IF (begin) THEN
c$       write(out_unitp,*) "F def thread",omp_get_thread_num()

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
c$OMP  END CRITICAL (d0d1d2_h_CRIT)
c      END initialization
c---------------------------------------------------------------------
c---------------------------------------------------------------------

       c_act = Qsym0(1)
       IF (deriv) THEN
         write(out_unitp,*) 'ERROR in d0d1d2_h'
         write(out_unitp,*) '  deriv CANNOT be true!!'
         write(out_unitp,*) ' check the fortran source'
         STOP
       ELSE

         DO i=1,mole%nb_inact2n
           DO j=i,mole%nb_inact2n
             d0h(i,j)=ZERO
             DO kl=1,nn(i,j)
               d0h(i,j) = d0h(i,j) + F(kl,i,j)*poly_legendre(c_act,kl,0)
             END DO
             d0h(j,i) = d0h(i,j)
c            write(out_unitp,*) 'd0h(i,j)',i,j,d0h(i,j)
           END DO
         END DO
       END IF

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
         CALL Write_RMat(d0h,6,4)
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

       integer :: iq,nb_act
       real (kind=Rkind) ::  Qact(nb_act)
       integer :: nderiv,it
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
c$       write(out_unitp,*) "F thread",omp_get_thread_num()
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
C    analytical derivative (dnQgene : Qgene Qgene' Qgene" Qgene'") calculation
c    for the variable iq
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

       integer :: i_qsym
       real (kind=Rkind) ::  Qsym0(mole%nb_var)

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
        write(out_unitp,*) 'i_qsym',i_qsym
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
c$OMP    CRITICAL (d0d1d2d3_Qeq_CRIT)
       IF (begin) THEN
c$       write(out_unitp,*) "F def thread",omp_get_thread_num()
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
c$OMP    END CRITICAL (d0d1d2d3_Qeq_CRIT)
c     end  initialisation
c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
       c_act  = Qsym0(1)
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


         DO kl=1,nn(i_qsym)
           CALL d0d1d2d3poly_legendre(c_act,kl,dc0,dc1,dc2,dc3,nderiv)

           d0req = d0req + F(kl,i_qsym) * dc0

           IF (nderiv .GE. 1)  d1req(i_act) =
     *                         d1req(i_act) + F(kl,i_qsym) * dc1

           IF (nderiv .GE. 2)  d2req(i_act,j_act) =
     *                         d2req(i_act,j_act) + F(kl,i_qsym)*dc2

           IF (nderiv .GE. 3)  d3req(i_act,j_act,k_act) =
     *                      d3req(i_act,j_act,k_act) + F(kl,i_qsym)*dc3

         END DO

c---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) ' F(.,i_qsym) i_qsym',i_qsym
        write(out_unitp,*) (F(k,i_qsym),k=1,nn(i_qsym))
        write(out_unitp,*) 'd0req : ',c_act,d0req
        write(out_unitp,*) 'd1req : ',c_act,d1req
        write(out_unitp,*) 'd2req : ',c_act,d2req
        write(out_unitp,*) 'd3req : ',c_act,d3req
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
       ScalOp(1) = Q(1)
       ScalOp(2) = Q(2)
       ScalOp(3) = Q(3)


       END
