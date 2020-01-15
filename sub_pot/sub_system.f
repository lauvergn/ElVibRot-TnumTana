c
c================================================================
c    calc_Op : calculation of the potential and dipolar matrices
c    mat_V(nb_be,nb_be) and mat_dip(nb_be,nb_be,3)
c    nb_be : nb of elctronic surface
c    Q are the coordinates in active order or syl order
c    dipolar calculation if calc_dip = T
c================================================================
      SUBROUTINE calcN_op(mat_V,mat_imV,mat_ScalOp,nb_be,nb_ScalOp,
     *                    Qcart,nb_Qcart,mole,
     *                    calc_ScalOp,pot_cplx)
      USE mod_Tnum
      USE mod_system
      USE mod_constant, only : get_Conv_au_TO_unit
      IMPLICIT NONE

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole


      integer           :: nb_be,nb_ScalOp,nb_Qcart
      logical           :: calc_ScalOp,pot_cplx
      real (kind=Rkind) :: mat_V(nb_be,nb_be),mat_imV(nb_be,nb_be)
      real (kind=Rkind) :: mat_ScalOp(nb_be,nb_be,nb_ScalOp)
      real (kind=Rkind) :: Qcart(nb_Qcart)

      real (kind=Rkind) :: pot0,im_pot0,x(3),x1(3),x2(3)
      real (kind=Rkind) :: pot2_H2atWn_WITH_cage_FlexSPCFW_SPCE
      integer           :: i,option
      integer           :: nat = 122
      real (kind=Rkind) :: XH2cage(3,122)
      integer, parameter :: n2dshell=60
      character (len=10) :: name_dum

      real (kind=Rkind), save :: X2dShell(3,60)
      real (kind=Rkind), save :: auTOcm=219474.631443_Rkind ! from Tnum
      real (kind=Rkind), save :: a0
      logical,           save :: begin=.TRUE.

      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.

      !mat_V(1,1) = ZERO
      !RETURN

      IF (nb_be == 1 ) THEN
c---------------------------------------------------------------
c$OMP    CRITICAL (calcN_op_CRIT)
        IF (begin) THEN
          a0 = get_Conv_au_TO_unit(quantity='L',Unit='Angs')
          auTOcm = get_Conv_au_TO_unit(quantity='E',Unit='cm-1')
          open(unit=99,file='2dShell.xyz')
          read(99,*)
          read(99,*)
          DO i=1,n2dshell
            read(99,*) name_dum,X2dShell(:,i)
            write(6,'(a,i0,3(x,F12.6))') 'xyz,i',i,X2dShell(:,i)
          END DO
          X2dShell = X2dShell / a0
          close(unit=99)
          begin=.FALSE.
        END IF
c$OMP    END CRITICAL (calcN_op_CRIT)
c---------------------------------------------------------------
        IF (debug) write(6,*) 'Qcart',Qcart
        !write(6,*) 'Qcart',Qcart(1:21)

        XH2cage(:,1:62)   = reshape(Qcart(1:186),shape=(/3,62/) )
        XH2cage(:,63:122) = X2dShell(:,:)
        
        IF (debug) THEN
          a0 = get_Conv_au_TO_unit(quantity='L',Unit='Angs')
          write(6,*) '==========================================='
          write(6,'(a,3(x,F12.6))') 'H ',XH2cage(:,1)*a0
          write(6,'(a,3(x,F12.6))') 'H ',XH2cage(:,2)*a0
          DO i=3,size(XH2cage(1,:)),3
            write(6,'(a,3(x,F12.6))') 'O ',XH2cage(:,i)*a0
            write(6,'(a,3(x,F12.6))') 'H ',XH2cage(:,i+1)*a0
            write(6,'(a,3(x,F12.6))') 'H ',XH2cage(:,i+2)*a0
          END DO
        END IF

c       CALL pot_H2_Cage(reshape(XH2cage,shape=[366]),mat_V(1,1))

        mat_V(1,1) = pot2_H2atWn_WITH_cage_FlexSPCFW_SPCE(XH2cage,
     *                                                    122,0,'./')

        IF (debug) write(6,*) 'nat,mat_V(1,1)',nat,mat_V(1,1)
        !write(6,*) 'nat,mat_V(1,1)',nat,mat_V(1,1)
        mat_V(1,1) = mat_V(1,1) / auTOcm
        !STOP
        
        IF (debug)  write(6,*) 'nat,mat_V(1,1)',nat,mat_V(1,1)

        IF (pot_cplx) mat_imV(1,1) = im_pot0(Qcart)
        IF (calc_ScalOp) THEN
          CALL sub_dipole(mat_ScalOp(1,1,:),Qcart,mole)
        END IF
      ELSE
        write(6,*) ' ERROR in calc_op'
        write(6,*) ' more than ONE electronic surface is impossible'
        write(6,*) ' Rq: nb_be',nb_be
        STOP
      END IF

      END
c
C================================================================
C    fonction pot0(x) 3+9 D pour h2o en cartesiennes (calcul direct)
C================================================================
      FUNCTION pot0(Q)
      USE mod_system
      USE mod_file
      IMPLICIT NONE

      integer, parameter :: ndim=6
      real (kind=Rkind) :: pot0
      real (kind=Rkind) :: Q(ndim)


      real (kind=Rkind) :: DQ(ndim)
      real (kind=Rkind) :: G(ndim)

      integer :: nio,i,err,idum
      real (kind=Rkind) :: V
   
      real (kind=Rkind), save :: Qref(ndim)
      real (kind=Rkind), save :: hess(ndim,ndim)

      logical, save     :: begin = .TRUE.

c---------------------------------------------------------------
c      initialisation la premiere fois
c$OMP  CRITICAL (pot0_CRIT)
       IF (begin) THEN
         CALL file_open2(name_file='hessien6D.inp',iunit=nio)

         read(nio,*) 

         DO i=1,ndim
           read(nio,*) idum,Qref(i)
         END DO
         write(6,*) 'Qref',Qref(:)

         read(nio,*) 
         CALL Read_Mat(hess,nio,5,err)
         IF (err /= 0) STOP 'error read hessien'

         close(nio)
         begin = .FALSE.
       END IF
c$OMP  END CRITICAL (pot0_CRIT)
c fin     initialisation la premiere fois
c---------------------------------------------------------------

      DQ(:) = Q(:)-Qref(:)
      V = HALF * dot_product(DQ,matmul(hess,DQ))
      pot0 = V
c     write(6,*) 'DQ,pot0',DQ,V

c     STOP 

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
      real (kind=Rkind) :: hh(n,n)

      hh = ZERO
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
C    analytical derivative (Qeq + derivatives) calculation
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
C    analytical derivative (dnQflex : Qflex + derivatives) calculation
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
      SUBROUTINE sub_dipole(dip,XH2,mole)
      USE mod_Tnum
      USE mod_system
      IMPLICIT NONE

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

      real (kind=Rkind) :: XH2(3,2)
      real (kind=Rkind) :: dip(3)


      real (kind=Rkind) :: G(3),R(3)

      G=HALF*(XH2(:,1)+XH2(:,2))
      R=(XH2(:,1)-XH2(:,2))

      dip(1) = G(1) + R(1)
      dip(2) = G(2) + R(2)
      dip(3) = G(3) + R(3)

      !dip(2) = R(3)/sqrt(dot_product(R,R))

      END
