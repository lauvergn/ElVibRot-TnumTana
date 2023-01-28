!================================================================
!    analytical derivative (Tab_dnQflex : Qflex(:) Qflex' Qflex" Qflex'") calculation
!    at Qact, for all nb_var coordinates
!================================================================
  SUBROUTINE calc_Tab_dnQflex(Tab_dnQflex,nb_var,Qact,nb_act,nderiv,it)
  USE mod_system
  USE mod_dnSVM
  IMPLICIT NONE

    integer           :: nb_var,nb_act
    real (kind=Rkind) :: Qact(nb_act)
    integer           :: nderiv,it
    TYPE (Type_dnS)   :: Tab_dnQflex(nb_var)

    integer :: iQ

    !----- for debuging ----------------------------------
    character (len=*), parameter :: name_sub='calc_Tab_dnQflex'
    logical, parameter :: debug=.FALSE.
    !logical, parameter :: debug=.TRUE.
    !----- for debuging ----------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) 'nb_act',nb_act
      write(out_unitp,*) 'nb_var',nb_var
      flush(out_unitp)
    END IF
    !---------------------------------------------------------------------

    write(out_unitp,*) 'The subroutine ',name_sub
    write(out_unitp,*) ' must be adapted to the sytem'
    flush(out_unitp)
    STOP ' ERROR in calc_Tab_dnQflex: the subroutine must be adapted.'

    !---------------------------------------------------------------------
    IF (debug) THEN
      DO iQ=1,nb_var
        write(out_unitp,*) 'tab_dnQflex : ',iQ,Qact
        CALL write_dnS(tab_dnQflex(iQ),nderiv)
      END DO
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF
    !---------------------------------------------------------------------

  END SUBROUTINE calc_Tab_dnQflex
