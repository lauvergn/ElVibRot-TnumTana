!================================================================
!    analytical derivative (dnQflex : Qflex Qflex' Qflex" Qflex'") calculation
!    for the variable iq
!================================================================
  SUBROUTINE calc_Tab_dnQflex(Tab_dnQflex,nb_var,Qact,nb_act,nderiv,it)
  USE mod_system
  USE mod_dnSVM      
  IMPLICIT NONE

       integer           :: nb_var,nb_act
       real (kind=Rkind) :: Qact(nb_act)
       integer           :: nderiv,it
       TYPE (Type_dnS)   :: Tab_dnQflex(nb_var)




       real (kind=Rkind) :: dc0,dc1,dc2,dc3
       real (kind=Rkind) :: c_act

       character (len=14) :: nom
       logical :: exist

       integer :: vi,kl,k,iQ

       integer, parameter      ::  max_points=200
       integer, parameter      ::  nb_inactb = 5
       real (kind=Rkind), save :: F(max_points,nb_inactb)
       integer, save           :: nn(nb_inactb)
       logical, save           :: begin = .TRUE.


!----- for debuging ----------------------------------
      character (len=*), parameter :: name_sub='calc_Tab_dnQflex'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!----- for debuging ----------------------------------


!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_act',nb_act
      END IF
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!      Qact value. Rq: only ONE active variable is possible
!---------------------------------------------------------------------
       IF (nb_act /= 1) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' the number of Active variable'
         write(out_unitp,*) ' should be 1. But nb_act =',nb_act
         STOP
       END IF

!---------------------------------------------------------------------
!---------------------------------------------------------------------


!---------------------------------------------------------------------
!---------------------------------------------------------------------
!      initialization (only once)
!$OMP    CRITICAL (dnQflex_CRIT)
       IF (begin) THEN
         write(out_unitp,*) ' INITIALIZATION of ',name_sub
         begin=.FALSE.
         nn(:) = 0
         F(:,:) = ZERO
         DO vi=2,3
           nom=nom_i('inter12___',vi)
           write(out_unitp,*) 'read file :',nom,vi

           CALL read_para0d(F(1,vi),nn(vi),max_points,nom,exist)
           IF ( .NOT. exist ) STOP

           !write(out_unitp,*) vi,(F(k,vi),k=1,nn(vi))
         END DO
         write(out_unitp,*) ' END INITIALIZATION of ',name_sub

       END IF
!$OMP    END CRITICAL (dnQflex_CRIT)
!     end  initialization
!---------------------------------------------------------------------
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      c_act = Qact(1)

      DO iQ=1,nb_var
       IF (iQ < 2 .OR. iQ > 3) CYCLE

       CALL sub_ZERO_TO_dnS(tab_dnQflex(iQ))

       DO kl=1,nn(iQ)
         CALL d0d1d2d3poly_legendre(c_act,kl,dc0,dc1,dc2,dc3,nderiv)

         tab_dnQflex(iQ)%d0 = tab_dnQflex(iQ)%d0 + F(kl,iQ) * dc0

         IF (nderiv >= 1) THEN
            tab_dnQflex(iQ)%d1(1)     = tab_dnQflex(iQ)%d1(1)     + F(kl,iQ)*dc1
         END IF

         IF (nderiv >= 2) THEN
           tab_dnQflex(iQ)%d2(1,1)   = tab_dnQflex(iQ)%d2(1,1)   + F(kl,iQ)*dc2
         END IF

         IF (nderiv >= 3) THEN
           tab_dnQflex(iQ)%d3(1,1,1) = tab_dnQflex(iQ)%d3(1,1,1) + F(kl,iQ)*dc3
         END IF

      END DO
    END DO

!---------------------------------------------------------------------
    IF (debug) THEN
      DO iQ=1,nb_var
        write(out_unitp,*) 'tab_dnQflex : ',iQ,Qact
        CALL write_dnS(tab_dnQflex(iQ),nderiv)
      END DO
      write(out_unitp,*) 'END ',name_sub
    END IF
!---------------------------------------------------------------------

  END SUBROUTINE calc_Tab_dnQflex
