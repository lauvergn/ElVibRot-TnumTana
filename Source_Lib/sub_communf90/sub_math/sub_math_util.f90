FUNCTION gamma_perso(n)
USE mod_system
IMPLICIT NONE
  real(kind=Rkind) :: gamma_perso
  real(kind=Rkind) a
  integer i,n
  IF (n < 0) THEN
   write(out_unitp,*) 'ERROR: gamma( n<=0)',n
   STOP
  END IF
  a = ONE
  DO i=1,n-1
   a = a * i
  END DO
  gamma_perso = a
end function gamma_perso
      FUNCTION factor(n)
      USE mod_system
      IMPLICIT NONE
      real(kind=Rkind) :: factor
         real(kind=Rkind) a
         integer i,n
         IF (n .LT. 0) THEN
           write(out_unitp,*) 'ERROR: factor( n<=0)',n
           STOP
         END IF
         a = ONE
         DO i=1,n
           a = a * real(i,kind=Rkind)
         END DO
         factor = a
      end function factor
      FUNCTION binomial(n,i)
      USE mod_system
      IMPLICIT NONE
      real(kind=Rkind) :: binomial
         real(kind=Rkind) a
         integer i,k,n
         IF (n .LT. 0 .OR. i .GT. n .OR. i .LT. 0) THEN
           write(out_unitp,*) 'ERROR: binomial( n<=0 i<0 i>n)',n,i
           STOP
         END IF
         a = ONE
         DO k=1,n
           a = a * real(k,kind=Rkind)
         END DO
         DO k=1,n-i
           a = a / real(k,kind=Rkind)
         END DO
         DO k=1,i
           a = a / real(k,kind=Rkind)
         END DO
         binomial = a

!        write(out_unitp,*) 'binomial',n,i,a
      end function binomial
      FUNCTION combi(n,i)
      USE mod_system
      IMPLICIT NONE
      integer i,k,n
      real(kind=Rkind) :: combi(n,i)
      real(kind=Rkind) a
         IF (n .LT. 0 .OR. i .GT. n .OR. i .LT. 0) THEN
           write(out_unitp,*) 'ERROR: combi( n<0 i<0 i>n)',n,i
           STOP
         END IF
         a = ONE
         DO k=1,n
           a = a * real(k,kind=Rkind)
         END DO
         DO k=1,n-i
           a = a / real(k,kind=Rkind)
         END DO
         DO k=1,i
           a = a / real(k,kind=Rkind)
         END DO
         combi = a

!        write(out_unitp,*) 'combi',n,i,a
      end function combi
      FUNCTION icombi(n,i)
      USE mod_system
      IMPLICIT NONE
         integer :: icombi
         integer :: a
         integer :: i,n
         integer :: f1,f2,k
         IF (n .LE. 0 .OR. i .GT. n .OR. i .LT. 0) THEN
           write(out_unitp,*) 'ERROR: icombi( n<=0 i<0 i>n)',n,i
           STOP
         END IF

         IF (i > n-i) THEN
           f1 = i
           f2 = n-i
         ELSE
           f2 = i
           f1 = n-i
         END IF
         a = 1
         DO k=f1+1,n
           a = a * k
         END DO
         DO k=1,f2
           a = a / k
         END DO
         icombi = a

!        write(out_unitp,*) 'combi',n,i,a
!        STOP
      end function icombi
