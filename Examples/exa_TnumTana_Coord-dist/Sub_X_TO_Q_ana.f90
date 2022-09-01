SUBROUTINE Q_TO_X_ana(Q,nb_Q,X,nb_X,inTOout)
  USE mod_system
  IMPLICIT NONE
  integer, intent(in)              :: nb_Q,nb_X
  real (kind=Rkind), intent(inout) :: Q(nb_Q)
  real (kind=Rkind), intent(inout) :: X(nb_X)
  logical, intent(in)              :: inTOout

  real (kind=Rkind) :: R1(3),R2(3)


  IF (inTOout) THEN

    X(:) = ZERO
    X(4:6) = [ZERO                  , ZERO, Q(1)     ]
    X(7:9) = [sqrt(ONE-Q(3)**2)*Q(2), ZERO, Q(3)*Q(2)]

  ELSE

    R1(:) = X(4:6)-X(1:3)
    R2(:) = X(7:9)-X(1:3)

    Q(1) = sqrt(dot_product(R1,R1))
    Q(2) = sqrt(dot_product(R2,R2))
    Q(3) = dot_product(R1,R2)/(Q(1)*Q(2))  ! because Q(3) is the cosine

  END IF

END SUBROUTINE Q_TO_X_ana
