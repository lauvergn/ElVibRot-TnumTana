 PROGRAM Main_TnumTana_FDriver
 IMPLICIT NONE


  integer, parameter :: nat=13
  real (kind=8) :: Qact(3*nat-6),Qcart(3*nat)


  integer :: i
  character (len=*), parameter :: name_sub='Main_TnumTana'



  Qact(:) = 0.5d0
  CALL Qact_TO_cart(Qact,size(Qact),Qcart,size(Qcart))

 !$OMP   PARALLEL &
 !$OMP   DEFAULT(NONE) &
 !$OMP   PRIVATE(i,Qact,Qcart)

 !$OMP   DO SCHEDULE(STATIC)
  DO i=1,10**6
    IF (mod(i,10) == 0) write(6,'(".")',advance='no')
    Qact(1) = 0.5d0 + real(i,kind=8)*0.001d0
    CALL Qact_TO_cart(Qact,size(Qact),Qcart,size(Qcart))
  END DO
 !$OMP   END DO
 !$OMP   END PARALLEL
 write(6,*)

  DO i=1,3*nat,3
    write(6,*) (i-1)/3+1,Qcart(i:i+2)
  END DO

 END PROGRAM Main_TnumTana_FDriver
