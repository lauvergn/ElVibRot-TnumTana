SUBROUTINE Q_TO_X_ana(Q,nb_Q,X,nb_X,inTOout)
  USE mod_system
  IMPLICIT NONE
  integer, intent(in)              :: nb_Q,nb_X
  real (kind=Rkind), intent(inout) :: Q(nb_Q)
  real (kind=Rkind), intent(inout) :: X(nb_X)
  logical, intent(in)              :: inTOout


  IF (inTOout) THEN  ! X(Q)

    X(:) = ZERO

  ELSE ! Q(X)

    Q(:) = ZERO

  END IF

  STOP 'The subroutine Q_TO_X_ana MUST be make.'


END SUBROUTINE Q_TO_X_ana
