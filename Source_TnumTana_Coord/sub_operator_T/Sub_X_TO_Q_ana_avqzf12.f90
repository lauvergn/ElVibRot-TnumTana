SUBROUTINE Q_TO_X_ana(Q,nb_Q,X,nb_X,inTOout)
  USE mod_system
  IMPLICIT NONE
  integer, intent(in)              :: nb_Q,nb_X
  real (kind=Rkind), intent(inout) :: Q(nb_Q)
  real (kind=Rkind), intent(inout) :: X(nb_X)
  logical, intent(in)              :: inTOout


       integer                    :: nb_act
       integer                    :: nat,ncart
       real (kind=8)              :: XG(3),vH2(3),R,th,phi


  IF (inTOout) THEN

    CALL cartesian(Q(1),Q(3),X)

  ELSE
    STOP
  END IF

END SUBROUTINE Q_TO_X_ana
        subroutine cartesian(al1,al3,xp)
        implicit none
        real*8 :: al1,al3
        real*8 x(12),xp(12),dx1(12),dx3(12)

        data x/
     C       0.000000000  , 0.000000000 , -1.138906322,
     C       0.000000000  , 0.000000000 ,  1.138906322,
     C       0.000000000  , 0.000000000 ,  3.147958336,
     C       0.000000000  , 0.000000000 , -3.147958336 /
        data dx1/
     C   2.25586426495d-18, 7.37208318208d-18, 8.84446685715d-02,
     C  -7.79760546727d-18,-1.02834010480d-17,-8.84446685715d-02,
     C  -1.54322165131d-16, 9.26059125202d-17, 6.34804844577d-01,
     C  -2.53198333823d-17, 4.75653781134d-18,-6.34804844577d-01 /
        data dx3/
     C   1.05582172833d-17,-7.72471469781d-18,-5.68178243690d-02,
     C  -9.31279490677d-18,-3.68652493552d-17,-5.68178243690d-02,
     C   3.41147965384d-18, 1.31308278774d-16, 6.76520122470d-01,
     C   8.79495911308d-17,-3.49347059094d-17, 6.76520122470d-01 /

        xp(:) = x + al1*dx1 + al2*dx2
 
        end
