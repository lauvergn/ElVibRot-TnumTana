c
c================================================================
c    calc_Op : calculation of the potential and dipolar matrices
c    mat_V(nb_be,nb_be) and mat_dip(nb_be,nb_be,3)
c    nb_be : nb of elctronic surface
c    Q are the coordinates in active order or syl order
c    dipolar calculation if calc_dip = T
c================================================================
      SUBROUTINE calcN_op(mat_V,mat_imV,mat_ScalOp,nb_be,nb_ScalOp,
     *                    Qact,nb_var,mole,
     *                    calc_ScalOp,pot_cplx)

      USE mod_Tnum
      USE mod_system
      USE mod_Constant
      IMPLICIT NONE

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole


      integer           :: nb_be,nb_ScalOp,nb_var
      logical           :: calc_ScalOp,pot_cplx
      real (kind=Rkind) :: mat_V(nb_be,nb_be),mat_imV(nb_be,nb_be)
      real (kind=Rkind) :: mat_ScalOp(nb_be,nb_be,nb_ScalOp)
      real (kind=Rkind) :: Qact(nb_var)

      integer,parameter :: ndim=24
      real (kind=Rkind) :: im_pot0
      real (kind=Rkind) :: Q(ndim)
      real (kind=Rkind) :: Q2(ndim)

      real (kind=Rkind), parameter :: w(ndim) = (/
     *  0.09357_Rkind, 0.0740_Rkind, 0.1273_Rkind, 0.1568_Rkind,
     *  0.0400_Rkind, 0.0589_Rkind, 0.0778_Rkind, 0.0968_Rkind,
     *  0.1157_Rkind, 0.1347_Rkind, 0.1536_Rkind, 0.1726_Rkind,
     *  0.1915_Rkind, 0.2105_Rkind, 0.2294_Rkind, 0.2484_Rkind,
     *  0.2673_Rkind, 0.2863_Rkind, 0.3052_Rkind, 0.3242_Rkind,
     *  0.3431_Rkind, 0.3621_Rkind, 0.3810_Rkind, 0.4000_Rkind /)

      real (kind=Rkind), parameter :: k1(ndim) = (/
     *  ZERO, -0.0964_Rkind, 0.0470_Rkind, 0.1594_Rkind,
     *  0.0069_Rkind, 0.0112_Rkind, 0.0102_Rkind, 0.0188_Rkind,
     *  0.0261_Rkind, 0.0308_Rkind, 0.0210_Rkind, 0.0265_Rkind,
     *  0.0196_Rkind, 0.0281_Rkind, 0.0284_Rkind, 0.0361_Rkind,
     *  0.0560_Rkind, 0.0433_Rkind, 0.0625_Rkind, 0.0717_Rkind,
     *  0.0782_Rkind, 0.0780_Rkind, 0.0269_Rkind, 0.0306_Rkind/)
      real (kind=Rkind) :: k2(ndim) = (/
     *  ZERO,0.1194_Rkind, 0.2012_Rkind, 0.0484_Rkind,
     *  -0.0069_Rkind, -0.0112_Rkind, -0.0102_Rkind, -0.0188_Rkind,
     *  -0.0261_Rkind, -0.0308_Rkind, -0.0210_Rkind, -0.0265_Rkind,
     *  -0.0196_Rkind, -0.0281_Rkind, -0.0284_Rkind, -0.0361_Rkind,
     *  -0.0560_Rkind, -0.0433_Rkind, -0.0625_Rkind, -0.0717_Rkind,
     *  -0.0782_Rkind, -0.0780_Rkind, -0.0269_Rkind, -0.0306_Rkind /)
      real (kind=Rkind), parameter :: delta=0.46165_Rkind
      real (kind=Rkind), parameter :: lambda=0.1825_Rkind

      real (kind=Rkind) :: eVTOau

      
      eVTOau = ONE/get_Conv_au_TO_unit(quantity='E',Unit='eV')
     
      !write(6,*) 'w (au)',w10a,w6a,w1,w9a

      Q  = Qact(4:ndim+3)
      Q2 = Q*Q

      IF (nb_be == 2 ) THEN
        mat_V(:,:) = ZERO
        mat_V(1,1) = HALF* dot_product(w(:),Q2(:))

        mat_V(2,2) = mat_V(1,1) + delta + dot_product(k2(:),Q(:))
        mat_V(1,1) = mat_V(1,1) - delta + dot_product(k1(:),Q(:))
        mat_V(1,2) = lambda*Q(1)
        mat_V(2,1) = lambda*Q(1)

        mat_V = mat_V * eVTOau

c       write(6,*) 'Q,V',Q,mat_V
        IF (pot_cplx) mat_imV(1,1) = im_pot0(Q,ndim)
c       write(6,*) 'Q,imV',Q,mat_imV
        IF (calc_ScalOp) THEN
          CALL sub_scalar(mat_ScalOp(1,1,:),nb_ScalOp,Q,ndim,
     *                    mole)
        END IF
      ELSE
        write(6,*) ' ERROR in calc_op'
        write(6,*) ' It needs two PES'
        write(6,*) ' Rq: nb_be',nb_be
        STOP
      END IF

      RETURN
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
      IMPLICIT NONE

      integer, parameter :: n = 9
      real (kind=Rkind) :: h(n,n),hh(n,n)
      integer  ::  err

      hh = zero
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
      FUNCTION im_pot0(Q,n)
      USE mod_system
      IMPLICIT NONE

       real (kind=Rkind) :: im_pot0
       integer :: i,n
       real (kind=Rkind) :: Q(n)
       real (kind=Rkind) :: Q0=8.d0,z

       z = 0.d0
       DO i=1,n
        IF (abs(Q(i)) > Q0) THEN
           z = z -(abs(Q(i)) - Q0)**3
        END IF
       END DO
       !write(6,*) Q,z

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
c
C================================================================
c    dipole read
C================================================================
      SUBROUTINE sub_scalar(scalar,nb_scalar,Q,n,mole)
      USE mod_Tnum
      USE mod_paramQ
      USE mod_system
      IMPLICIT NONE

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

      integer           :: n,nb_scalar
      real (kind=Rkind) :: Q(n)
      real (kind=Rkind) :: scalar(nb_scalar)



      scalar(:)   = ZERO

      END
