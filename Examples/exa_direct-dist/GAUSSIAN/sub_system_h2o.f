c
C================================================================
C    calc_Op : calculation of the potential and dipolar matrices
c    mat_V(nb_be,nb_be) and mat_dip(nb_be,nb_be,3)
c    nb_be : nb of elctronic surface
c    Q are the coordinates in active order or syl order
c    dipolar calculation if calc_dip = T
C================================================================
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

      STOP 'calcN_op not with direct!!'

      IF (nb_be == 1 ) THEN
        mat_V(1,1) = pot0(Q)
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

      integer, parameter :: ndim=6
      real (kind=Rkind) :: pot0
      real (kind=Rkind) :: Q(ndim)


       pot0 = ZERO

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
      USE mod_OTF
      IMPLICIT NONE

      integer, parameter :: n = 9
      real (kind=Rkind) :: hh(n,n)
      integer  ::  err
      TYPE (param_file) :: file_FChk
      integer  :: nio
      logical  :: located

      integer  ::  i,j


      file_FChk%name = 'h2o_freq.fchk'
      file_FChk%old = .TRUE.
c     - read the hessain matrix
      CALL file_open(file_FChk,nio)


      CALL Find_Label(nio,'Cartesian Force Constants',located)
      IF (located) THEN
        read(nio,*,iostat=err)
        read(nio,*,iostat=err) ((hh(i,j),i=1,j),j=1,n)
      END IF
      IF (.NOT. located .OR. err /=0) THEN
        write(out_unitp,*) 'ERROR in sub_hessian'
        write(out_unitp,*) 'I cannot find the hessian in :',
     *                        file_FChk%name
        write(out_unitp,*) 'located,err',located,err
        STOP
      END IF
!     CALL Write_Mat(hh,out_unitp,5)
      DO j=1,n
      DO i=1,j-1
        hh(j,i) = hh(i,j)
      END DO
      END DO
      CALL Write_Mat(hh,out_unitp,5)

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
      real (kind=Rkind) :: dip(3)

       logical :: begin = .TRUE.
       SAVE begin

       IF (begin) THEN
         begin = .FALSE.
         open(3,file='dip3D',form='formatted')
       END IF

      read(3,*) dip(1)
      read(3,*) dip(2)
      read(3,*) dip(3)

      END
