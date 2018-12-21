c
C================================================================
C    calc_Op : calculation of the potential and dipolar matrices
c    mat_V(nb_be,nb_be) and mat_dip(nb_be,nb_be,3)
c    nb_be : nb of elctronic surface
c    Q are the coordinates in active order or syl order
c    dipolar calculation if calc_dip = T
C================================================================
      SUBROUTINE calcN_op(mat_V,mat_imV,mat_dip,nb_be,nb_Dip,
     *                    Qzmt,nb_var,mole,calc_dip,pot_cplx)

      USE mod_Coord_KEO
      implicit none

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole


      integer       :: nb_be,nb_var,nb_Dip
      logical       :: calc_dip,pot_cplx
      real (kind=8) :: mat_V(nb_be,nb_be),mat_imV(nb_be,nb_be)
      real (kind=8) :: mat_dip(nb_be,nb_be,nb_Dip)
      real (kind=8) :: Qzmt(nb_var)
      real (kind=8) :: Qloc(nb_var)

      real (kind=8) :: im_pot0
      real (kind=8) :: dip(3)
      integer, save :: iq = -1


      IF (nb_be == 1 ) THEN
        !write(6,*) 'Qzmt',Qzmt
        Qloc(1) = Qzmt(2) ! NH1
        Qloc(2) = Qzmt(4) ! NH2
        Qloc(3) = Qzmt(1) ! SiN
        Qloc(4) = Qzmt(3) ! H1NSi
        Qloc(5) = Qzmt(5) ! H2NSi
        Qloc(6) = Qzmt(6) - 3.141592653589793d0
        CALL pot_ori(mat_V(1,1),Qloc)
        !write(6,*) 'Qloc',Qloc,mat_V(1,1)
        iq = iq + 1
        !write(99,*) 'grid',iq,Qzmt,mat_V(1,1)
        !flush(99)
        IF (pot_cplx) mat_imV(1,1) = im_pot0(Qloc)
        IF (calc_dip ) THEN
          CALL sub_dipole(dip,Qloc)
          mat_dip(1,1,1) = dip(1)
          mat_dip(1,1,2) = dip(2)
          mat_dip(1,1,3) = dip(3)
        END IF
      END IF
      !write(6,*) 'Qloc',Qloc,mat_V(1,1)

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

       integer,parameter :: max_points = 210
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
         CALL file_open2(name_file='h2sinf12a.pot',iunit=nio)
         read(nio,*) nb_func
         IF (max_points < nb_func) THEN
           write(6,*) ' ERROR in pot0'
           write(6,*) ' nb_func > max_points',nb_func,max_points
           STOP
         END IF

         k = 0
         DO
          nb_columns = min(6,nb_func-k)
c         write(6,*) k+1,k+nb_columns,nb_columns
          IF (nb_columns == 0) exit
           read(nio,11) (tab_func(1:ndim,j),F(j),j=k+1,k+nb_columns)
c11        format(6(2x,6i1,f14.8))
 11        format(6i1,f15.8,5(2x,6i1,f15.8))
           k = k + nb_columns
         END DO
         read(nio,*) Qref(:)
         close(nio)
         begin=.FALSE.
c        write(6,*) 'Qref', Qref(:)
c        DO j=1,nb_func
c          write(6,*) j,': ',tab_func(1:ndim,j),F(j)
c        END DO
       END IF
c$OMP  END CRITICAL (pot0_CRIT)
c fin     initialisation la premiere fois
c---------------------------------------------------------------

      DQ(:) = Qdyn(:) - Qref(:)
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
c     V = dot_product(Qdyn,Qdyn)

c     write(66,*) Qdyn(1:6),V

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

      USE mod_Coord_KEO
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

      USE mod_system
      USE mod_Coord_KEO
      implicit none

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

       integer i_qact
       real*8  Qact1,Qsym0(mole%nb_var)

       real*8 c_act,s_act,ss


       real*8 d0h(mole%nb_inact2n,mole%nb_inact2n)
       real*8 d1h(mole%nb_inact2n,mole%nb_inact2n)
       real*8 d2h(mole%nb_inact2n,mole%nb_inact2n)
       real*8 step

       real*8 poly_legendre
c      character*14 nom_ii
       character*14 nom
       logical exist

       integer   max_points,nb_inactb
       parameter (max_points=200)
       parameter (nb_inactb=5)
       real*8 F(max_points,nb_inactb,nb_inactb)
       integer nn(nb_inactb,nb_inactb)


       integer jv,i_sym_act,i,j,iv,i_type_act

       logical deriv,num,begin
       data begin/.true./
       SAVE begin,F,nn,i_type_act

c----- for debuging ----------------------------------
      logical debug
      parameter (debug=.FALSE.)
c     parameter (debug=.TRUE.)
c---------------------------------------------------------------------
      IF (debug) THEN
      write(6,*)
      write(6,*) 'BEGINNING d0d1d2_h'
      write(6,*) 'nb_var',mole%nb_var
      write(6,*) 'nb_act1',mole%nb_act1
      write(6,*) 'nb_inact22,nb_inact21',mole%nb_inact22,mole%nb_inact21
      write(6,*) 'nb_inact2n',mole%nb_inact2n
      END IF
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c---------------------------------------------------------------------
c      initialisation la premiere fois sauf si direct=.true.
       IF (begin) THEN

         IF (nb_inactb .LT. mole%nb_inact2n ) THEN
           write(6,*) 'ERROR : nb_inactb is TO small',nb_inactb
           write(6,*) 'it should at least equal to',mole%nb_inact2n
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
             write(6,*) 'F(1,i,i) tq d0h = 1 and F(1,i,j) tq d0h = 0'
             nn(i,j) = 1
             IF ( i .EQ. j ) THEN
               F(1,i,j) = 1.d0/poly_legendre(1.d0,1,0)
             ELSE
               F(1,j,i) = 0.d0
               F(1,i,j) = 0.d0
             END IF
           END IF

         END DO
         END DO

c        -------------------------------------------------------------
         write(6,*) ' WARNING : I suppose the ONLY active variable '
c        write(6,*) '           is an angle'
c        i_type_act =  3
         write(6,*) '           is cos(angle)'
         i_type_act = -3
c        -------------------------------------------------------------

         begin=.FALSE.
       END IF
c fin     initialisation
c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
       i_qact = mole%nb_act1
       i_sym_act = mole%liste_QactTOQsym(1)
       Qact1  = Qsym0(i_sym_act)
c      write(6,*) 'test d0h',deriv
c      write(6,*) 'i_qact i_sym_act Qact1',i_qact,i_sym_act,Qact1

       IF (i_type_act .GT. 0) THEN
         c_act = dcos(Qact1)
         s_act = dsin(Qact1)

         ss    = s_act * s_act
       ELSE
         c_act = Qact1
       END IF
c---------------------------------------------------------------------

       d0h(:,:)=0.d0
c---------------------------------------------------------------------

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


      USE mod_Coord_KEO
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


       real*8 s_act,ss,sc,sss

       character*14 nom
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

       real*8 Q(21)
       real*8 dip(3)

c     ----------------------------------------------------------------
      real*8 pi,pi2
      parameter (pi = 3.141592653589793238462643383279
     *                          50288419716939937511d0)
      parameter (pi2 =pi+pi)

      real*8 sq2pi,sqpi
c     ----------------------------------------------------------------

       real*8 x
       integer kl


       character*14 nom
       logical exist

       integer max_points,nn
       parameter (max_points=30)
       real*8 Fx(max_points)
       real*8 Fy(max_points)
       real*8 Fz(max_points)

       logical begin
       data begin/.true./
       SAVE begin,Fx,Fy,Fz,nn

       dip(:) = 0.d0
       RETURN
c---------------------------------------------------------------
c      initialisation la premiere fois
       IF (begin) THEN
         nom='inter27-mux'
         CALL read_para0d(Fx,nn,max_points,nom,exist)
         IF ( .NOT. exist) STOP
         nom='inter28-muy'
         CALL read_para0d(Fy,nn,max_points,nom,exist)
         IF ( .NOT. exist) STOP
         nom='inter27-muz'
         CALL read_para0d(Fz,nn,max_points,nom,exist)
         IF ( .NOT. exist) STOP

         begin=.FALSE.
       END IF
c fin     initialisation la premiere fois
c---------------------------------------------------------------
       sqpi = 1.d0/sqrt(pi)
       sq2pi = 1.d0/sqrt(pi+pi)
       x = Q(1)


       dip(1) = Fx(1)*sq2pi
       DO kl=2,nn
         dip(1) = dip(1) + Fx(kl) * dcos((kl-1)*x)*sqpi
       END DO

       dip(2) = 0.d0
       DO kl=1,nn
         dip(2) = dip(2) + Fy(kl) * dsin((kl)*x)*sqpi
       END DO

       dip(3) = Fz(1)*sq2pi
       DO kl=2,nn
         dip(3) = dip(3) + Fz(kl) * dcos((kl-1)*x)*sqpi
       END DO


c      dip(1) = 0.d0
c      dip(2) = 0.d0
c      dip(3) = 0.d0


       END
