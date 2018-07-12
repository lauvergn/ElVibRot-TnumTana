SUBROUTINE Q_TO_X_ana(Q,nb_Q,X,nb_X,inTOout)
  USE mod_system
  IMPLICIT NONE
  integer, intent(in)              :: nb_Q,nb_X
  real (kind=Rkind), intent(in)    :: Q(nb_Q)
  real (kind=Rkind), intent(inout) :: X(nb_X)
  logical, intent(in)              :: inTOout

  IF (nb_Q /= 15 .OR. nb_X < 21) THEN
    write(out_unitp,*) ' ERROR in Q_TO_X_ana'
    write(out_unitp,*) ' For protonated water dimer subroutine "q2x_RisOO".'
    write(out_unitp,*) ' nb_Q /= 15',nb_Q
    write(out_unitp,*) ' nb_X < 21',nb_X
    STOP
  END IF

  IF (inTOout) THEN
    CALL q2x_RisOO(Q,X(1:21))
  ELSE
    write(out_unitp,*) ' ERROR in Q_TO_X_ana'
    write(out_unitp,*) ' For protonated water dimer subroutine "q2x_RisOO".'
    write(out_unitp,*) ' inTOout=f is not possible'
    STOP
  END IF

END SUBROUTINE Q_TO_X_ana
! C O O R D I N A T E   T R A N S F O R M A T I O N S
!      real (kind=8) :: q(0:14)
!      real (kind=8) :: x(0:20)
!      character (len=10) :: name_v
!
!      DO i=0,14
!        read(in_unitp,*) name_v,q(i)
!      END DO
!      CALL q2x_RisOO(q,x)
!      DO i=0,20,3
!        write(out_unitp,*) x(i+0),x(i+1),x(i+2)
!      END DO
!      END
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

