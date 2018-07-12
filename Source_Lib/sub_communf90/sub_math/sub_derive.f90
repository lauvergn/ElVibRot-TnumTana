!
!=====================================================================
!
!  ++   derivee 1ere et seconde en x
!
!=====================================================================
!
      SUBROUTINE d1d2(d0,d1,d2,step)
      USE mod_system
      IMPLICIT NONE

       real(kind=Rkind) d0,ep,em,d1,d2,step

       ep=d1
       em=d2
       d1 = (ep-em)/(step+step)
       d2 = (ep+em-d0-d0)/(step*step)

       end subroutine d1d2
!
!=====================================================================
!
!  ++   derive 1er et second en x
!
!=====================================================================
!
      SUBROUTINE d1d2_b(d0,ep,em,d1,d2,step)
      USE mod_system
      IMPLICIT NONE

       real(kind=Rkind) d0,ep,em,d1,d2,step

       d1 = (ep-em)/(step+step)
       d2 = (ep+em-d0-d0)/(step*step)

       end subroutine d1d2_b

