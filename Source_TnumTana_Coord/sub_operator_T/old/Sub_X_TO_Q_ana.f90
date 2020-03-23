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

        data x/ &
             0.000000000  , 0.000000000 , -1.139410848 ,&
             0.000000000  , 0.000000000 ,  1.139410848 ,&
             0.000000000  , 0.000000000 ,  3.148893846 ,&
             0.000000000  , 0.000000000 , -3.148893846 /
        data dx1/&
        -2.010519459090D-17,-9.730796852890D-18, 8.835574708410D-02,&
         6.174747737470D-17, 2.398192948080D-17,-8.835574708410D-02,&
         1.591357846140D-16,-1.866174769650D-17, 6.349522676520D-01,&
        -2.213969178550D-18,-1.228283829360D-17,-6.349522676520D-01 /
        data dx3/&
        1.427997226810D-17,-1.377773418730D-17,-5.681782436900D-02,&
        8.256861611810D-17, 6.231783581810D-18,-5.681782436900D-02,&
        3.016169137950D-16,-3.133896886330D-17, 6.765201224700D-01,&
       -8.518403579320D-17,-1.179905255980D-17, 6.765201224700D-01 /


        xp(:) = x + al1*dx1 + al3*dx3
 
        end
