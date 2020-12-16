!======================================================================
!
!      Calculation of Tdef Tcor Trot at Q
!      Tdef = Tdef2 * d2./dQ1dQ2 + Tdef1 * d./dQ1 + vep
!      Tcor = Tcor2 * d./dQ1*Jx  + Tcor1 * Jx
!      Trot = Trot  * Jx*Jy
!
!      Calculation of rho at Q
!      dT = rho*dQ
!
!
!======================================================================
      SUBROUTINE calc_f2_f1Q_ana(Qsym0,                                 &
                                 Tdef2,Tdef1,vep,rho,                   &
                                 Tcor2,Tcor1,Trot,                      &
                                 para_Tnum,mole)
      USE mod_system
      USE mod_Tnum
      USE mod_Constant
      IMPLICIT NONE

!----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole
      TYPE (Tnum)    :: para_Tnum

      real (kind=Rkind) ::  Qsym0(mole%nb_var)

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!--- variables for Tnum --------------------------------------------------
!
!      Tdef = Tdef2(1,2) * d2./dQ1dQ2 + Tdef1(1) * d./dQ1 + vep
!      Tcor = Tcor2(1,x) * d./dQ1*Jx  + Tcor1(x) * Jx
!      Trot = Trot(x,y)  * Jx*Jy
!
!      Calculation of rho at Q
!      dT = rho*dQ
!
!
!      nrho: type of normalization
!            1 => Wilson (rho = 1.)
!           10 => Wilson (rho = 1.) (without vep (vep=0))
!            2 => clever choice ( Q=R =>1.; Q=val => sin(val); Q=cos(val) =>1. Q=dih =>1.)
!           20 => id 2 but without vep (vep=0)
!            0 => Euclidian (rho = Jac)
!
!     JJ: Total angular momentum
!
!
!     stepT: displacement for numerical calculation of Tnum
!            see num_H,num_A,num_x
!
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
      real (kind=Rkind) :: Tdef2(mole%nb_act,mole%nb_act)
      real (kind=Rkind) :: Tdef1(mole%nb_act)
      real (kind=Rkind) :: vep
      real (kind=Rkind) :: Tcor2(mole%nb_act,3),Tcor1(3)
      real (kind=Rkind) :: Trot(3,3)

      real (kind=Rkind) :: rho

      real (kind=Rkind) :: Tdef2_tot(mole%nb_var,mole%nb_var)
      real (kind=Rkind) :: Tdef1_tot(mole%nb_var)

      integer           :: i

      integer, parameter :: ndim=12
      real (kind=Rkind), parameter :: w(ndim) = (/               &
        0.09357_Rkind, 0.0740_Rkind, 0.1273_Rkind, 0.1568_Rkind, &
        0.1347_Rkind, 0.3431_Rkind, 0.1157_Rkind, 0.3242_Rkind,  &
        0.3621_Rkind, 0.2673_Rkind, 0.3052_Rkind, 0.0968_Rkind/)

      real (kind=Rkind) :: eVTOau
!-------------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      !-----------------------------------------------------------
       IF (debug .OR. para_Tnum%WriteT) THEN
         write(out_unitp,*) 'BEGINNING calc_f2_f1Q_ana'
         write(out_unitp,*) 'ndimG',mole%ndimG
         write(out_unitp,*) 'WriteCC',mole%WriteCC
         write(out_unitp,*) 'Qsym0',Qsym0
         IF (debug) THEN
           write(out_unitp,*)
           CALL Write_mole(mole)
           write(out_unitp,*)
         END IF
         write(out_unitp,*) 'JJ',para_Tnum%JJ
         write(out_unitp,*)
       END IF
!-----------------------------------------------------------

      eVTOau = ONE/get_Conv_au_TO_unit(quantity='E',Unit='eV')

      Tdef2(:,:) = ZERO
      Tdef1(:)   = ZERO
      vep        = ZERO
      rho        = ZERO
      Tcor2(:,:) = ZERO
      Tcor1(:)   = ZERO
      Trot(:,:)  = ZERO

      DO i=1,ndim
        Tdef2(i,i) = -HALF*w(i)*eVTOau
      END DO

      !write(out_unitp,*) 'Tdef',Tdef2

!-----------------------------------------------------------
      IF (debug .OR. para_Tnum%WriteT) THEN
        write(out_unitp,*) 'END calc_f2_f1Q_ana'
      END IF
!-----------------------------------------------------------


      RETURN
      end subroutine calc_f2_f1Q_ana

      SUBROUTINE q2x_RisOO(qq,x)
      USE mod_system
      IMPLICIT NONE
      real (kind=Rkind) :: qq(0:14)
      real (kind=Rkind) :: x(0:20)
      integer :: i

      CALL q2x_RisOO_Oriol_OK(qq,x)
!     CALL q2x_RisOO_Oriol(qq,x)

      !write(out_unitp,*) ' Oriol subroutine: W2H+'
      !DO i=0,20,3
      !  write(out_unitp,*) x(i+0),x(i+1),x(i+2)
      !END DO

      end subroutine q2x_RisOO

      SUBROUTINE q2x_RisOO_Oriol_OK(qq,x)
      USE mod_system
      IMPLICIT NONE
      real (kind=Rkind) :: q(0:14),qq(0:14)
      real (kind=Rkind) :: dum,q_OriolOrder(0:14)
      real (kind=Rkind) :: x(0:20)
      real (kind=Rkind), dimension (0:2) :: xHp,xOa,xH1a,                   &
           xH2a,xCMa,xOb,xH1b,xH2b,xCMb
      real (kind=Rkind) :: R(0:2,0:2)

      logical :: z_prim = .FALSE. ! IF true => z'=z/(R-3)


      q_OriolOrder(:) = qq(:)
      q_OriolOrder(4) = qq(2) !a
      q_OriolOrder(2) = qq(4) !z
!     write(out_unitp,*) 'a,z',q_OriolOrder(4),q_OriolOrder(2)



!   Coordinate transformation, from internal to cartesian

!        Internals are given in the following order:

!        x:        q(0)
!        y:        q(1)
!        z:        q(2)

!        R:        q(3)
!        alpha:    q(4)

!        R1_a:     q(5)
!        R2_a:     q(6)
!        theta_a:  q(7)

!        beta_a:   q(8)
!        lambda_a: q(9)

!        R1_b:     q(10)
!        R2_b:     q(11)
!        theta_b:  q(12)

!        beta_b:   q(13)
!        lambda_b: q(14)
!   # Transform angles given as cosines back to radians
       CALL qcos2q(q_OriolOrder,q)
       q(0) = q(0) !x
       q(1) = q(1) !y
       q(2) = q(2) !z
       q(3) = q(3) !R
       IF (z_prim) q(2) = q(2) * (q(3)-THREE)
       q(4) = q(4) !a

       q(5) = q(5) !r1a
       q(6) = q(6) !r2a
       q(7) = q(7) !tha
       q(8) = q(8) !ba
       q(9) = q(9) !la

       q(10) = q(10) !r1b
       q(11) = q(11) !r2b
       q(12) = q(12) !thb
       q(13) = q(13) !bb
       q(14) = q(14) !lb
!   # Initialize empty array for the final cartesian coordinates
       x(:) = ZERO
!   # Initialize xyz arrays for the cartesians of each atom
       xHp = ZERO
       xOa = ZERO
       xH1a= ZERO
       xH2a= ZERO
       xCMa= ZERO
       xOb = ZERO
       xH1b= ZERO
       xH2b= ZERO
       xCMb= ZERO
!   # Central proton in BFA
      xHp(0) = q(0)
      xHp(1) = q(1)
      xHp(2) = q(2)
!   # Wat A
!   # Set watA on the XZ_BFA plane, R2 parallel to Z_BFA
!   # and calculate the position of the oxygen in relation
!   # to the R2 vector
      xH1a(2) = HALF*q(6)
      xH2a(2) =-HALF*q(6)
      xOa(0) = sin(q(7))*q(5)
      xOa(2) = cos(q(7))*q(5)
!   # Shift atoms so that the oxygen lies at (0,0,0)_BFA
      xCMa = xOa
      xOa = xOa - xCMa
      xH1a = xH1a - xCMa
      xH2a = xH2a - xCMa
!   # Apply Euler rotations
      CALL eulerMatrix(ZERO,q(8),q(9),R,.FALSE.)
      xOa = matmul(R,xOa)
      xH1a = matmul(R,xH1a)
      xH2a = matmul(R,xH2a)
!   # Displace watA by 0.5*R to the right (Wrong F Gatti)
!   # Displace watA by 0.5*R to the left
      xOa(2)  = xOa(2)  - HALF*q(3)
      xH1a(2) = xH1a(2) - HALF*q(3)
      xH2a(2) = xH2a(2) - HALF*q(3)
!   # Wat B
!   # Set watB on the XZ_BFB plane, R2 parallel to Z_BFB
!   # and calculate the position of the oxygen in relation
!   # to the R2 vector
      xH1b(2) = HALF*q(11)
      xH2b(2) =-HALF*q(11)
      xOb(0) = sin(q(12))*q(10)
      xOb(2) = cos(q(12))*q(10)
!   # Shift atoms so that the oxygen lies at (0,0,0)_BFB
      xCMb = xOb
      xOb = xOb - xCMb
      xH1b = xH1b - xCMb
      xH2b = xH2b - xCMb
!   # Apply Euler rotations
      CALL eulerMatrix(q(4),q(13),q(14),R,.FALSE.)
      xOb = matmul(R,xOb)
      xH1b = matmul(R,xH1b)
      xH2b = matmul(R,xH2b)
!   # Displace watB by 0.5*R to the left (Wrong)
!   # Displace watB by 0.5*R to the right
      xOb(2)  = xOb(2)  + HALF*q(3)
      xH1b(2) = xH1b(2) + HALF*q(3)
      xH2b(2) = xH2b(2) + HALF*q(3)
!   # Set the total x vector
      x(0)  = xOa(0)
      x(1)  = xOa(1)
      x(2)  = xOa(2)
      x(3)  = xOb(0)
      x(4)  = xOb(1)
      x(5)  = xOb(2)
      x(6)  = xHp(0)
      x(7)  = xHp(1)
      x(8)  = xHp(2)
      x(9)  = xH1b(0)
      x(10) = xH1b(1)
      x(11) = xH1b(2)
      x(12) = xH2b(0)
      x(13) = xH2b(1)
      x(14) = xH2b(2)
      x(15) = xH1a(0)
      x(16) = xH1a(1)
      x(17) = xH1a(2)
      x(18) = xH2a(0)
      x(19) = xH2a(1)
      x(20) = xH2a(2)
      end subroutine q2x_RisOO_Oriol_OK

      SUBROUTINE q2x_RisOO_Oriol(qq,x)
      USE mod_system
      IMPLICIT NONE
      real (kind=Rkind) :: q(0:14),qq(0:14)
      real (kind=Rkind) :: dum,q_OriolOrder(0:14)
      real (kind=Rkind) :: x(0:20)
      real (kind=Rkind), dimension (0:2) :: xHp,xOa,xH1a,                   &
           xH2a,xCMa,xOb,xH1b,xH2b,xCMb
      real (kind=Rkind) :: R(0:2,0:2)


      q_OriolOrder(:) = qq(:)
      q_OriolOrder(4) = qq(2) !a
      q_OriolOrder(2) = qq(4) !z
!     write(out_unitp,*) 'a,z',q_OriolOrder(4),q_OriolOrder(2)



!   Coordinate transformation, from internal to cartesian

!        Internals are given in the following order:

!        x:        q(0)
!        y:        q(1)
!        z:        q(2)

!        R:        q(3)
!        alpha:    q(4)

!        R1_a:     q(5)
!        R2_a:     q(6)
!        theta_a:  q(7)

!        beta_a:   q(8)
!        lambda_a: q(9)

!        R1_b:     q(10)
!        R2_b:     q(11)
!        theta_b:  q(12)

!        beta_b:   q(13)
!        lambda_b: q(14)
!   # Transform angles given as cosines back to radians
       CALL qcos2q(q_OriolOrder,q)
       q(0) =  q(0) !x
       q(1) =  q(1) !y
       q(2) = -q(2) !z
       q(3) = q(3) !R
       q(4) = q(4) !a

       q(5) = q(5) !r1a
       q(6) = q(6) !r2a
       q(7) = q(7) !tha
       q(8) = q(8) !ba
       q(9) = pi-q(9) !la

       q(10) = q(10) !r1b
       q(11) = q(11) !r2b
       q(12) = q(12) !thb
       q(13) = q(13) !bb
       q(14) = pi-q(14) !lb
!   # Initialize empty array for the final cartesian coordinates
       x(:) = ZERO
!   # Initialize xyz arrays for the cartesians of each atom
       xHp = ZERO
       xOa = ZERO
       xH1a= ZERO
       xH2a= ZERO
       xCMa= ZERO
       xOb = ZERO
       xH1b= ZERO
       xH2b= ZERO
       xCMb= ZERO
!   # Central proton in BFA
      xHp(0) = q(0)
      xHp(1) = q(1)
      xHp(2) = q(2)
!   # Wat A
!   # Set watA on the XZ_BFA plane, R2 parallel to Z_BFA
!   # and calculate the position of the oxygen in relation
!   # to the R2 vector
      xH1a(2) = HALF*q(6)
      xH2a(2) =-HALF*q(6)
      xOa(0) = sin(q(7))*q(5)
      xOa(2) = cos(q(7))*q(5)
!   # Shift atoms so that the oxygen lies at (0,0,0)_BFA
      xCMa = xOa
      xOa = xOa - xCMa
      xH1a = xH1a - xCMa
      xH2a = xH2a - xCMa
!   # Apply Euler rotations
      CALL eulerMatrix(ZERO,q(8),q(9),R,.FALSE.)
      xOa = matmul(R,xOa)
      xH1a = matmul(R,xH1a)
      xH2a = matmul(R,xH2a)
!   # Displace watA by 0.5*R to the right
      xOa(2) = xOa(2) + HALF*q(3)
      xH1a(2) = xH1a(2) + HALF*q(3)
      xH2a(2) = xH2a(2) + HALF*q(3)
!   # Wat B
!   # Set watB on the XZ_BFB plane, R2 parallel to Z_BFB
!   # and calculate the position of the oxygen in relation
!   # to the R2 vector
      xH1b(2) = HALF*q(11)
      xH2b(2) =-HALF*q(11)
      xOb(0) = sin(q(12))*q(10)
      xOb(2) = cos(q(12))*q(10)
!   # Shift atoms so that the oxygen lies at (0,0,0)_BFB
      xCMb = xOb
      xOb = xOb - xCMb
      xH1b = xH1b - xCMb
      xH2b = xH2b - xCMb
!   # Apply Euler rotations
      CALL eulerMatrix(q(4),q(13),q(14),R,.FALSE.)
      xOb = matmul(R,xOb)
      xH1b = matmul(R,xH1b)
      xH2b = matmul(R,xH2b)
!   # Displace watB by 0.5*R to the left
      xOb(2) = xOb(2) - HALF*q(3)
      xH1b(2) = xH1b(2) - HALF*q(3)
      xH2b(2) = xH2b(2) - HALF*q(3)
!   # Set the total x vector
      x(0)  = xOa(0)
      x(1)  = xOa(1)
      x(2)  = xOa(2)
      x(3)  = xOb(0)
      x(4)  = xOb(1)
      x(5)  = xOb(2)
      x(6)  = xHp(0)
      x(7)  = xHp(1)
      x(8)  = xHp(2)
      x(9)  = xH1b(0)
      x(10) = xH1b(1)
      x(11) = xH1b(2)
      x(12) = xH2b(0)
      x(13) = xH2b(1)
      x(14) = xH2b(2)
      x(15) = xH1a(0)
      x(16) = xH1a(1)
      x(17) = xH1a(2)
      x(18) = xH2a(0)
      x(19) = xH2a(1)
      x(20) = xH2a(2)
      end subroutine q2x_RisOO_Oriol

      SUBROUTINE qcos2q(q,qout)
      USE mod_system
      IMPLICIT NONE
      real (kind=Rkind) :: q(0:14),qout(0:14)
!     For the angles defined as cosines transform them to radians
      qout(:) = q(:)
      qout(7) = acos(qout(7))
      qout(8) = acos(qout(8))
      qout(12) = acos(qout(12))
      qout(13) = acos(qout(13))
      end subroutine qcos2q


      SUBROUTINE eulerMatrix(a,b,c,r,inv)
      USE mod_system
      IMPLICIT NONE
      real (kind=Rkind) :: a,b,c
      logical :: inv
      real (kind=Rkind) :: r(0:2,0:2)
      real (kind=Rkind) :: Am(0:2,0:2)
      real (kind=Rkind) :: Bm(0:2,0:2)
      real (kind=Rkind) :: Cm(0:2,0:2)

!def eulerMatrix(a,b,c,inv=False)
!   """Return a 3x3 rotation matrix

!   xyz is the body fixed frame
!   XYZ is the space fixed frame

!   The normal convention (default) is the rotation XYZ -> xyz

!   - First rotate around Z by an ammount of xhi,
!   - then around Y by an ammount of theta
!   - and finally around Z by an ammount of phi.

!   The theta and phi angles can also be seen as the two spherical angles
!   of a vector with respect to the Z vector.
!        phi       azimutal spherical angle (0,2pi)
!        theta     polar spherical angle (0,pi)

!   phi   (a) -- rotation around Z axis
!   theta (b) -- rotation around Y axis
!   xhi   (c) -- rotation around Z axis

!   The inverse convention corresponds to xyz -> XYZ
!   The vector is rotated by the specified Euler angles that align it
!   with the space fixed frame.

!   Normal and Inverse conventions are inverse operations with respect to
!   each other
!   """
!   # Initilization of the rotation matrices:
      Am = ZERO
      Bm = ZERO
      Cm = ZERO

      Am(0,0) = cos(a)
      Am(1,1) = cos(a)
      Am(2,2) = ONE
      Bm(0,0) = cos(b)
      Bm(1,1) = ONE
      Bm(2,2) = cos(b)
      Cm(0,0) = cos(c)
      Cm(1,1) = cos(c)
      Cm(2,2) = ONE
      IF (.not. inv) THEN
        Am(0,1) = -sin(a)
        Am(1,0) = sin(a)
        Bm(0,2) = sin(b)
        Bm(2,0) = -sin(b)
        Cm(0,1) = -sin(c)
        Cm(1,0) = sin(c)
        r=matmul(matmul(Am,Bm),Cm)
      ELSE
        Am(0,1) = sin(a)
        Am(1,0) = -sin(a)
        Bm(0,2) = -sin(b)
        Bm(2,0) = sin(b)
        Cm(0,1) = sin(c)
        Cm(1,0) = -sin(c)
        r=matmul(matmul(Cm,Bm),Am)
      END IF
      end subroutine eulerMatrix

