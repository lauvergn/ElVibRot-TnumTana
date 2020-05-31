c
c================================================================
c    calc_Op : calculation of the potential and dipolar matrices
c    mat_V(nb_be,nb_be) and mat_dip(nb_be,nb_be,3)
c    nb_be : nb of electronic surface(s)
c    Qact are the coordinates in active order
c    Scalar operator (dipole moment) calculation if calc_dip = T
c================================================================
      SUBROUTINE calcN_op(mat_V,mat_imV,mat_ScalOp,nb_be,nb_ScalOp,
     *                    Qact,nb_var,mole,
     *                    calc_ScalOp,pot_cplx)

      USE mod_Tnum
      USE mod_system
      IMPLICIT NONE

c----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole


      integer           :: nb_be,nb_ScalOp,nb_var
      logical           :: calc_ScalOp,pot_cplx
      real (kind=Rkind) :: mat_V(nb_be,nb_be),mat_imV(nb_be,nb_be)
      real (kind=Rkind) :: mat_ScalOp(nb_be,nb_be,nb_ScalOp)
      real (kind=Rkind) :: Qact(nb_var)

      real (kind=Rkind) :: im_pot0
      real (kind=Rkind) :: x,y
      integer           :: i,ndim,nsurf
      real (kind=Rkind), parameter :: E0=HALF
      !real (kind=Rkind), parameter :: V12=ZERO
      real (kind=Rkind), parameter :: V12=HALF
     
        x= Qact(1)
        y= Qact(2)

        mat_V(1,1) = HALF*(x**2+y**2) - E0
        IF (nb_be > 1) THEN
          mat_V(2,2) = HALF*(x**2+y**2) + E0
          mat_V(1,2) = V12
          mat_V(2,1) = V12
        END IF


        !write(66,*) 'pot',Qact(1),mat_V

        IF (pot_cplx) THEN
          mat_imV(:,:) = ZERO
          DO i=1,nb_be
            mat_imV(i,i) = im_pot0(Qact,ndim)
          END DO
        END IF

      END
C================================================================
C    subroutine calculant le gradient
C================================================================
      SUBROUTINE d0d1d2_g(d0g,d1g,d2g,Qsym0,mole,deriv,num,step)
      USE mod_Tnum
      USE mod_system
      IMPLICIT NONE

c----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole

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

c----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole

      
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

c----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole

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

c----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole

      integer           :: n,nb_scalar
      real (kind=Rkind) :: Q(n)
      real (kind=Rkind) :: scalar(nb_scalar)



      scalar(:)   = ZERO

      END
