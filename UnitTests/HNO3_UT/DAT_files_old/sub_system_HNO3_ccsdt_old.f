c
C================================================================
C    calc_Op : calculation of the potential and dipolar matrices
c    mat_V(nb_be,nb_be) and mat_dip(nb_be,nb_be,3)
c    nb_be : nb of elctronic surface
c    Q are the coordinates in active order or syl order
c    dipolar calculation if calc_dip = T
C================================================================
      SUBROUTINE calcN_op(mat_V,mat_imV,mat_ScalOp,nb_be,nb_ScalOp,
     *                   Qdyn,nb_var,mole,
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
      real (kind=Rkind) :: Qdyn(nb_var)
      real (kind=Rkind) :: Qsym(nb_var)

      real (kind=Rkind) :: pot0,im_pot0,pot0_9DQflex
      real (kind=Rkind) :: dip(3)

      IF (nb_be == 1 ) THEN
        !write(6,*) 'mole%nb_Qtransfo',mole%nb_Qtransfo
        IF (mole%nb_Qtransfo == 4) THEN ! sym + NM
          !write(6,*) 'name_transfo',mole%tab_Qtransfo(3)%name_transfo
          Qsym(:) = matmul(mole%tab_Qtransfo(3)%LinearTransfo%mat,Qdyn)
        ELSE ! only sym
          Qsym(:) = Qdyn(:)
        END IF
!       write(6,*) 'Qsym for pot',Qsym(:)
        mat_V(1,1) = pot0(Qsym)
!       mat_V(1,1) = pot0_9DQflex(Qsym,mole)
        IF (pot_cplx) mat_imV(1,1) = im_pot0(Qsym)
        IF (calc_ScalOp) THEN
          CALL sub_dipole(dip,Qsym)
          mat_ScalOp(1,1,1) = dip(1)
          mat_ScalOp(1,1,2) = dip(2)
          mat_ScalOp(1,1,3) = dip(3)
        END IF
      END IF

      RETURN
      END
C================================================================
C================================================================
      FUNCTION pot0_9DQflex(Qsym,mole)
      USE mod_system
      USE mod_Tnum
      implicit none

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

      real (kind=Rkind) :: pot0_9DQflex,V


      integer, parameter :: ndim=12
      integer, parameter :: ndimh=8

      real (kind=Rkind) :: Qsym(ndim),Qact1(1)
      real (kind=Rkind) :: DQ(ndim)
      real (kind=Rkind) :: d0v,d1v,d2v,d3v,vh
      integer       :: i1,i2

      !IF (mole%nb_act1 /= 9) STOP ' CANNOT be use with constraints'
      Qact1  = Qsym(9)
      CALL d0d1d2d3vfour(Qact1,d0v,d1v,d2v,d3v,0,0)

      V = d0v
      DO i1=1,ndimh
        CALL d0d1d2d3vfour(Qact1,d0v,d1v,d2v,d3v,i1,0)
        DQ(i1) = Qsym(i1)-d0v
      END DO

      vh = ZERO
      DO i1=1,ndimh
      DO i2=1,ndimh
        CALL d0d1d2d3vfour(Qact1,d0v,d1v,d2v,d3v,i1,i2)
        vh = vh + d0v * DQ(i1)*DQ(i2)
      END DO
      END DO
      V = V + vh * HALF
      !write(6,*) 'vh',Qact1,Qsym(8),vh*HALF,V

c---------------------------------------------------------------------

      pot0_9DQflex = V
c     write(66,*) Qact1,V

      END

C================================================================
C================================================================
      FUNCTION pot0(Q)
      USE mod_system
      IMPLICIT NONE
      real (kind=Rkind) :: pot0,V
      integer, parameter :: ndim=12

      real (kind=Rkind) :: Q(ndim)
      real (kind=Rkind) :: Qact1
      real (kind=Rkind) :: d0v,d1v,d2v,d3v

      CALL d0d1d2d3vfour(Q(9),d0v,d1v,d2v,d3v,0,0)

      pot0 = d0v
c     write(66,*) Qact1,d0v

      END
C================================================================
c      Si jq=0, on regarde seulement iq
c      avec iq (ou jq) E [0....8]
c      iq=0 : energy
c      iq=1,...8 : r1,r+,a+,d2,r-,a-,rh,ah
c   
C================================================================
      SUBROUTINE d0d1d2d3vfour(Qact1,d0v,d1v,d2v,d3v,iq,jq)
      USE mod_system

c     ----------------------------------------------------------------
      real (kind=Rkind) Qact1(1)
      real (kind=Rkind) d0v,d1v,d2v,d3v
      
c     ----------------------------------------------------------------
      real (kind=Rkind) :: d0,d1,d2,d3
      integer kl
      integer i,iq,jq


      integer, parameter :: ndim=1
      integer, parameter :: max_fit=8
      integer, parameter :: max_nn=10
      real (kind=Rkind)      :: F(max_nn,0:max_fit,0:max_fit)
      integer            :: nn(0:ndim,0:max_fit,0:max_fit)
      integer            :: nt(0:max_fit,0:max_fit)
      character*50 nom
      logical exist

       logical begin
       data begin/.true./
       SAVE begin,F,nn,nt

c---------------------------------------------------------------
c      initialisation la premiere fois
!$OMP  CRITICAL (d0d1d2d3vfour_CRIT)
       IF (begin) THEN
         F(:,:,:) = 0.d0
         nn(:,:,:) = 0
         nt(:,:) = 0

c        energy + Qeq(phi)
         jj=0
         DO ii=0,max_fit
           write(nom,'("HESS-fit/inter_",i1)') ii
           flush(6)
           CALL read_para4d(F(:,ii,jj),nn(:,ii,jj),ndim,nt(ii,jj),
     *                      max_nn,nom,exist)
           IF ( .NOT. exist) STOP
         END DO


c        hess(i,j)
         DO jj=1,max_fit
         DO ii=1,max_fit
           write(nom,'("HESS-fit/inter_",i1,"_",i1)') ii,jj
           flush(6)
           CALL read_para4d(F(:,ii,jj),nn(:,ii,jj),ndim,nt(ii,jj),
     *                      max_nn,nom,exist)
           IF ( .NOT. exist) STOP
         END DO
         END DO

         begin=.FALSE.
       END IF
!$OMP  END CRITICAL (d0d1d2d3vfour_CRIT)
c fin     initialisation la premiere fois
c---------------------------------------------------------------

      IF (iq > max_fit .OR. iq < 0 .OR. jq > max_fit .OR. jq < 0) THEN
        write(6,*) ' ERROR in d0d1d2d3vfour'
        write(6,*) ' wrong value for iq or jq E [0...8]',iq,jq
        STOP
      END IF
 
      d0v = 0
      d1v = 0
      d2v = 0
      d3v = 0
      DO i=1,nn(0,iq,jq)
        IF (nt(iq,jq) == 41) kl = 2*i-1      ! cos(kl.x)
        IF (nt(iq,jq) == 42) kl = 4*i-3      ! cos(2.kl.x)
        IF (nt(iq,jq) == 51) kl = 2*i-0      ! sin(kl.x)
        IF (nt(iq,jq) == 52) kl = 4*i-0    ! sin(2.kl.x)

        CALL d0d1d2d3fourier(Qact1(1),d0,d1,d2,d3,kl)
c       write(6,*) 'kl,nn,F,d0',kl,nn(0,n_fit),F(i,n_fit),d0
        d0v = d0v + F(i,iq,jq) * d0
        d1v = d1v + F(i,iq,jq) * d1
        d2v = d2v + F(i,iq,jq) * d2
        d3v = d3v + F(i,iq,jq) * d3
      END DO

      END
C================================================================
C    fonction im_pot0(x) imaginary part of pot0
C================================================================
       real (kind=Rkind) FUNCTION im_pot0(Qsym0)
      USE mod_system


       real (kind=Rkind) Qsym0(1)

       im_pot0 = 0.0

       END
C================================================================
C    fonction pot_rest(x) rest of the DL : pot0 + v2 + pot_rest
C================================================================
       real (kind=Rkind) FUNCTION pot_rest(Qact,Delta_Qact,nb_inact2n)
      USE mod_system

       real (kind=Rkind) Qact(1)
       integer nb_inact2n
       real (kind=Rkind) Delta_Qact(nb_inact2n)

       pot_rest = 0.0

      END
C================================================================
C    subroutine calculant le gradient
C================================================================
       SUBROUTINE d0d1d2_g(d0g,d1g,d2g,Qsym0,mole,deriv,num,step)
      USE mod_system

      USE mod_Tnum
      implicit none

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

      real (kind=Rkind) :: d0g(mole%nb_inact2n)
      real (kind=Rkind) :: d1g(mole%nb_inact2n,mole%nb_act1)
      real (kind=Rkind) :: 
     *        d2g(mole%nb_inact2n,mole%nb_act1,mole%nb_act1)

      real (kind=Rkind) :: Qsym0(mole%nb_var)
      real (kind=Rkind) :: step
      logical       :: deriv,num




      real (kind=Rkind) :: Qact(mole%nb_act1)
      real (kind=Rkind) :: d0f
      real (kind=Rkind) :: d1f(mole%nb_act1)
      real (kind=Rkind) :: d2f(mole%nb_act1,mole%nb_act1)
      integer       :: i,k,vi,ibo

      real (kind=Rkind) :: d0gzmt(mole%nb_var)

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
      USE mod_system

       real (kind=Rkind) h

       h = 0.d0


       RETURN
       END
       SUBROUTINE H0_sym(h)
      USE mod_system

       real (kind=Rkind) h

       RETURN
       END
       SUBROUTINE d0d1d2_h(d0h,d1h,d2h,
     *                     Qsym0,mole,deriv,num,step)
      USE mod_system

      USE mod_Tnum
      implicit none

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

       integer i_qact
       real (kind=Rkind)  Qact1(1),Qsym0(mole%nb_var)
       real (kind=Rkind)  step
       logical :: deriv,num

       integer :: i1,i2,i_qsym1,i_qsym2,i_sym_act,n_fit
       real (kind=Rkind) d0hzmat(9,9)

       real (kind=Rkind) d0h(mole%nb_inact2n,mole%nb_inact2n)
       real (kind=Rkind) d1h(mole%nb_inact2n,mole%nb_inact2n)
       real (kind=Rkind) d2h(mole%nb_inact2n,mole%nb_inact2n)

      real (kind=Rkind) d0v,d1v,d2v,d3v

c----- for debuging ----------------------------------
c     logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
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
       i_qact = mole%nb_act1
       i_sym_act = mole%liste_QactTOQsym(1)
       Qact1(1)  = Qsym0(i_sym_act)


      d0h(:,:)=0.d0

      DO i1=1,mole%nb_inact21
      DO i2=1,mole%nb_inact21
        i_qsym1 = mole%liste_QactTOQsym(i1+1)
        i_qsym2 = mole%liste_QactTOQsym(i2+1)
        
        CALL d0d1d2d3vfour(Qact1,d0v,d1v,d2v,d3v,i_qsym1,i_qsym2)
        d0h(i1,i2) = d0v
      END DO
      END DO

c---------------------------------------------------------------------

      IF (debug) THEN
      write(6,*) 'd0h'
      CALL ecriture(d0h,mole%nb_inact21,mole%nb_inact21,
     *              5,.TRUE.,mole%nb_inact21)
      write(6,*) 'END d0d1d2_h'
      END IF

      END
C================================================================
C    analytical derivative (Qeq Qeq' Qeq" Qeq'") calculation
c    for the variable i_qsym
C================================================================
      SUBROUTINE calc_dnQflex(iq,dnQflex,Qact,nb_act,nderiv,it)
      USE mod_system
      USE mod_dnSVM
      IMPLICIT NONE

       integer :: iq,nb_act
       real (kind=Rkind) ::  Qact(nb_act)

       integer :: nderiv,it ! it number of the transformation

       TYPE (Type_dnS)   :: dnQflex

       real (kind=Rkind) :: d0v,d1v,d2v,d3v

c----- for debuging ----------------------------------
      character (len=*), parameter :: name_sub='dnQflex'
      logical, parameter :: debug=.FALSE.
c     logical, parameter :: debug=.TRUE.
c----- for debuging ----------------------------------


c---------------------------------------------------------------------
      IF (debug) THEN
        write(6,*) 'BEGINNING ',name_sub
        write(6,*) 'nb_act',nb_act
        write(6,*) 'iq',iq
        flush(6)
      END IF
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c      Qact value. Rq: only ONE active variable is possible
c---------------------------------------------------------------------
       IF (nb_act /= 1) THEN
         write(6,*) ' ERROR in ',name_sub
         write(6,*) ' the number of Active variable'
         write(6,*) ' should be 1. But nb_act =',nb_act
         STOP
       END IF

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
       dnQflex%d0                          = ZERO
       IF (nderiv >= 1)  dnQflex%d1(:)     = ZERO
       IF (nderiv >= 2)  dnQflex%d2(:,:)   = ZERO
       IF (nderiv >= 3)  dnQflex%d3(:,:,:) = ZERO
c---------------------------------------------------------------------

       CALL d0d1d2d3vfour(Qact,d0v,d1v,d2v,d3v,iq,0)

       IF (nderiv >= 0) dnQflex%d0 = d0v
       IF (nderiv >= 1) dnQflex%d1 = d1v
       IF (nderiv >= 2) dnQflex%d2 = d2v
       IF (nderiv == 3) dnQflex%d3 = d3v


       !write(6,*) 'Qflex',iq,Qact(1),dnQflex%d0

c---------------------------------------------------------------------
      IF (debug) THEN
        write(6,*) 'dnQflex : ',Qact
        CALL write_dnS(dnQflex,nderiv)
        write(6,*) 'END ',name_sub
        flush(6)
      END IF
c---------------------------------------------------------------------

      END
      SUBROUTINE d0d1d2d3_Qeq(i_qsym,
     *                        d0req,d1req,d2req,d3req,
     *                        Qsym,mole,nderiv)
      USE mod_system


      USE mod_Tnum
      implicit none

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

      integer i_qsym
      real (kind=Rkind)  Qact1(1),Qsym(mole%nb_var)

      integer nderiv

      real (kind=Rkind) d0req
      real (kind=Rkind) d1req(mole%nb_act1)
      real (kind=Rkind) d2req(mole%nb_act1,mole%nb_act1)
      real (kind=Rkind) d3req(mole%nb_act1,mole%nb_act1,mole%nb_act1)


      real (kind=Rkind) d0v,d1v,d2v,d3v
      integer nio


c----- for debuging ----------------------------------
c     logical, parameter ::  debug = .TRUE.
      logical, parameter ::  debug = .FALSE.
c---------------------------------------------------------------------
      IF (debug) THEN
        write(6,*) 'BEGINNING d0d1d2d3_Qeq'
        write(6,*) 'nb_inact20,nb_act',mole%nb_inact20,mole%nb_act
        write(6,*) 'nb_var',mole%nb_var
        write(6,*) 'i_qsym',i_qsym
        write(6,*) 'nderiv',nderiv
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


      Qact1(1) = Qsym(9)
      CALL d0d1d2d3vfour(Qact1,d0v,d1v,d2v,d3v,i_qsym,0)

c---------------------------------------------------------------------
       d0req = 0.d0
       d1req = 0.d0
       d2req = 0.d0
       d3req = 0.d0

       IF (nderiv >= 0) d0req = d0v
       IF (nderiv >= 1) d1req = d1v
       IF (nderiv >= 2) d2req = d2v
       IF (nderiv == 3) d3req = d3v
c---------------------------------------------------------------------

c     nio = 660+i_qsym
c     write(nio,*) 'd0req',i_qsym,Qact1,d0v

      IF (debug) THEN
        write(6,*) 'd0req',d0req
        write(6,*) 'END d0d1d2d3_Qeq'
      END IF

       END
C================================================================
C    The tri component of the dipole moment.
C================================================================
       SUBROUTINE sub_dipole(dip,Q)
      USE mod_system

       real (kind=Rkind) Q(21)
       real (kind=Rkind) dip(3)

c     ----------------------------------------------------------------
      real (kind=Rkind) pi2
      parameter (pi2 =pi+pi)

      real (kind=Rkind) sq2pi,sqpi
c     ----------------------------------------------------------------

       real (kind=Rkind) x,z
       integer kl


       character*14 nom
       logical exist

       integer max_points,nn
       parameter (max_points=30)
       real (kind=Rkind) Fx(max_points)
       real (kind=Rkind) Fy(max_points)
       real (kind=Rkind) Fz(max_points)

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
