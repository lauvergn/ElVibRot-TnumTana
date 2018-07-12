!
!=====================================================================
!
! ++   correction au 1er ordre de l energie avec v
!
!=====================================================================
!
      SUBROUTINE ene_1er(d1w,v,nb_act)
      USE mod_system
      IMPLICIT NONE

!      nombre de variable active
       integer i,nb_act
       real(kind=Rkind) v(nb_act,nb_act)
       real(kind=Rkind) d1w(nb_act)

!      correction au 1er ordre avec d0k' => derivee 1er de w
       DO i=1,nb_act
         d1w(i) = v(i,i)
       END DO

       end subroutine ene_1er
!
!
!=====================================================================
!
! ++   matrice de perturbation pour le 2d ordre
!
!=====================================================================
!
      SUBROUTINE pert1(v1,v0,w,nb_act)
      USE mod_system
      IMPLICIT NONE

!      nombre de variable active
       integer i,j,nb_act
       real(kind=Rkind) v0(nb_act,nb_act)
       real(kind=Rkind) v1(nb_act,nb_act)
       real(kind=Rkind) w(nb_act)


       DO i=1,nb_act
         v1(i,i) = ZERO
         DO j=1,i-1
           v1(i,j) = v0(i,j)/(w(i)-w(j))
           v1(j,i) = -v1(i,j)
         END DO
!        DO j=i+1,nb_act
!          v1(i,j) = v0(i,j)/(w(i)-w(j))
!        END DO
       END DO

       end subroutine pert1
!
!=====================================================================
!
! ++   matrice de perturbation pour le 2d ordre
!
!=====================================================================
!
      SUBROUTINE pert2(v2,v1,v0,w,nb_act)
      USE mod_system
      IMPLICIT NONE

!      nombre de variable active
       integer i,j,k,nb_act
       real(kind=Rkind) v2(nb_act,nb_act)
       real(kind=Rkind) v1(nb_act,nb_act)
       real(kind=Rkind) v0(nb_act,nb_act)
       real(kind=Rkind) w(nb_act)

!      correction au 2d  ordre avec d0k'  => derivee 2d de c

       DO i=1,nb_act

         DO j=1,nb_act
           IF (j .NE. i) THEN
!            correction d2c(i,j)
             v2(i,j) = -v0(i,i)*v1(i,j)
             DO k=1,nb_act
               IF (k .NE. i)                                            &
                 v2(i,j) = v2(i,j) + v0(j,k)*v1(i,k)
             END DO
             v2(i,j) = v2(i,j)/(w(i)-w(j))
           ELSE
!            correction d2c(i,i)
             v2(i,i) = ZERO
             DO k=1,nb_act
               IF (k .NE. i)                                            &
                 v2(i,i) = v2(i,i) - v1(i,k)*v1(i,k)
             END DO
             v2(i,i) = v2(i,i)*HALF
           END IF
         END DO
       END DO
       end subroutine pert2
!
!
!=====================================================================
!
!  ++  correction au 2er ordre avec v
!
!=====================================================================
!
      SUBROUTINE ene_2dbis(d2w,d0w,v1,v2,nb_act)
      USE mod_system
      IMPLICIT NONE

!      nombre de variable active
       integer i,j,nb_act
       real(kind=Rkind) v1(nb_act,nb_act)
       real(kind=Rkind) v2(nb_act,nb_act)
       real(kind=Rkind) d2w(nb_act)
       real(kind=Rkind) d0w(nb_act)

!      correction au 2d  ordre avec v
       DO i=1,nb_act
         d2w(i) = ZERO
         DO j=1,i-1
           d2w(i) = d2w(i) + v1(i,j)*v2(i,j)/(d0w(i)-d0w(j))
         END DO
         DO j=i+1,nb_act
           d2w(i) = d2w(i) + v1(i,j)*v2(i,j)/(d0w(i)-d0w(j))
         END DO
       END DO

       end subroutine ene_2dbis
!
!
!=====================================================================
!
! ++   correction au 2er ordre avec v0 et v1(i,j) = v0(i,j)/(d0w(i)-d0w(j))
!
!=====================================================================
!
      SUBROUTINE ene_2d(d2w,v1,v0,nb_act)
      USE mod_system
      IMPLICIT NONE

!      nombre de variable active
       integer i,j,nb_act
       real(kind=Rkind) v0(nb_act,nb_act)
       real(kind=Rkind) v1(nb_act,nb_act)
       real(kind=Rkind) d2w(nb_act)

!      correction au 2d  ordre avec v
       DO i=1,nb_act
         d2w(i) = ZERO
         DO j=1,nb_act
           d2w(i) = d2w(i) + v0(i,j)*v1(i,j)
         END DO
       END DO

       end subroutine ene_2d
!
!
!=====================================================================
!
! ++   correction au 1er ordre avec vij(i,j) = v(i,j)/(d0w(i)-d0w(j))
!
!=====================================================================
!
      SUBROUTINE vec_1er(d1c,d0c,vij,nb_act)
      USE mod_system
      IMPLICIT NONE

!      nombre de variable active
       integer i,j,k,nb_act
       real(kind=Rkind) vij(nb_act,nb_act)
       real(kind=Rkind) d1c(nb_act,nb_act)
       real(kind=Rkind) d0c(nb_act,nb_act)

!      correction au 1er ordre avec d0k'  => derivee 1er de c
       DO i=1,nb_act
         DO k=1,nb_act
           d1c(k,i) = ZERO
         END DO
       END DO
       DO i=1,nb_act
         DO j=1,nb_act
           DO k=1,nb_act
             d1c(k,i) = d1c(k,i) + vij(i,j)*d0c(k,j)
           END DO
         END DO
       END DO

      end subroutine vec_1er

