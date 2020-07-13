!===========================================================================
!===========================================================================
!This file is part of Tnum-Tana.
!
!    Tnum-Tana is a free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    Tnum-Tana is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ElVibRot.  If not, see <http://www.gnu.org/licenses/>.
!
!    Copyright 2015  David Lauvergnat
!      Tnum is written David Lauvergnat [1]
!      Tana is written by Mamadou Ndong [1] and David Lauvergnat [1]
!         with contributions
!          Emil Lund klinting (coupling with MidasCpp) [3]'
!
![1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France
![3]: Department of Chemistry, Aarhus University, DK-8000 Aarhus C, Denmark
!
!===========================================================================
!===========================================================================
      MODULE mod_OneDTransfo
      use mod_system
      USE mod_dnSVM
      IMPLICIT NONE

      PRIVATE

      !!@description: TODO
      !!@param: TODO
      TYPE Type_oneDTransfo
        integer                    :: iQin        = 0
        integer                    :: type_oneD   = 0     ! identity
        logical                    :: inTOout     =.TRUE. ! T => no inversion (inTOout), F => inversion (outTOin)
        character (len=Name_len)   :: name_oneD   = "identity"
        real (kind=Rkind), pointer :: cte(:)      => null()
        integer, pointer           :: opt_cte(:)  => null()
      END TYPE Type_oneDTransfo

      INTERFACE alloc_array
        MODULE PROCEDURE alloc_array_OF_OneDTransfodim1
      END INTERFACE
      INTERFACE dealloc_array
        MODULE PROCEDURE dealloc_array_OF_OneDTransfodim1
      END INTERFACE

      PUBLIC :: Type_oneDTransfo, alloc_oneDTransfo, dealloc_oneDTransfo, &
                Read_oneDTransfo, Write_oneDTransfo, calc_oneDTransfo,    &
                oneDTransfo1TOoneDTransfo2, alloc_array, dealloc_array

      CONTAINS

!=======================================================================
!     oneD transfo
!=======================================================================
      SUBROUTINE alloc_oneDTransfo(oneDTransfo,nb_transfo)

      TYPE (Type_oneDTransfo), pointer, intent(inout) :: oneDTransfo(:)
      integer, intent(in) :: nb_transfo

      integer :: it
      character (len=*), parameter :: name_sub='alloc_oneDTransfo'

      IF (associated(oneDTransfo)) THEN
        CALL dealloc_oneDTransfo(oneDTransfo)
      END IF
      IF (nb_transfo < 1) RETURN

      CALL alloc_array(oneDTransfo,(/nb_transfo/),"oneDTransfo",name_sub)

      DO it=1,nb_transfo
        CALL alloc_array(oneDTransfo(it)%cte,(/20/),                    &
                        "oneDTransfo(it)%cte",name_sub)

        CALL alloc_array(oneDTransfo(it)%opt_cte,(/20/),                &
                        "oneDTransfo(it)%opt_cte",name_sub)
      END DO

      END SUBROUTINE alloc_oneDTransfo
      SUBROUTINE dealloc_oneDTransfo(oneDTransfo)

      TYPE (Type_oneDTransfo), pointer, intent(inout) :: oneDTransfo(:)

      integer :: it
      character (len=*), parameter :: name_sub='dealloc_oneDTransfo'

      IF (.NOT. associated(oneDTransfo)) RETURN

      DO it=1,size(oneDTransfo)
        CALL dealloc_array(oneDTransfo(it)%cte,                         &
                          "oneDTransfo(it)%cte",name_sub)
        CALL dealloc_array(oneDTransfo(it)%opt_cte,                     &
                          "oneDTransfo(it)%opt_cte",name_sub)
      END DO
      CALL dealloc_array(oneDTransfo,"oneDTransfo",name_sub)

      END SUBROUTINE dealloc_oneDTransfo

      SUBROUTINE alloc_array_OF_OneDTransfodim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      TYPE (Type_oneDTransfo), pointer, intent(inout) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=1
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_OneDTransfodim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (associated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)

         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1)),stat=err_mem)
       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Type_oneDTransfo')

      END SUBROUTINE alloc_array_OF_OneDTransfodim1
      SUBROUTINE dealloc_array_OF_OneDTransfodim1(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (Type_oneDTransfo), pointer, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_OneDTransfodim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN
       IF (.NOT. associated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Type_oneDTransfo')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_OneDTransfodim1

      SUBROUTINE oneDTransfo1TOoneDTransfo2(oneDTransfo1,oneDTransfo2)

      !-- oneDTransfo --------------------------------------
      TYPE (Type_oneDTransfo), pointer, intent(in) :: oneDTransfo1(:)
      TYPE (Type_oneDTransfo), pointer, intent(inout) :: oneDTransfo2(:)

      integer :: it
      character (len=*), parameter :: name_sub = 'oneDTransfo1TOoneDTransfo2'

      CALL dealloc_oneDTransfo(oneDTransfo2)
      IF (.NOT. associated(oneDTransfo1)) RETURN

      CALL alloc_oneDTransfo(oneDTransfo2,size(oneDTransfo1))

      DO it=1,size(oneDTransfo2)

        oneDTransfo2(it)%iQin         = oneDTransfo1(it)%iQin
        oneDTransfo2(it)%type_oneD    = oneDTransfo1(it)%type_oneD
        oneDTransfo2(it)%inTOout      = oneDTransfo1(it)%inTOout
        oneDTransfo2(it)%name_oneD    = oneDTransfo1(it)%name_oneD
        oneDTransfo2(it)%cte          = oneDTransfo1(it)%cte
        oneDTransfo2(it)%opt_cte      = oneDTransfo1(it)%opt_cte

      END DO

      END SUBROUTINE oneDTransfo1TOoneDTransfo2


      SUBROUTINE Read_oneDTransfo(oneDTransfo,nb_transfo,nb_Qin)

      TYPE (Type_oneDTransfo), pointer, intent(inout) :: oneDTransfo(:)
      integer, intent(in) :: nb_Qin,nb_transfo

      integer :: i,it,nb_flex_act,err,nbcol

      logical                    :: inTOout
      integer                    :: iQin,type_oneD
      character (len=Name_len)   :: name_oneD
      real (kind=Rkind)          :: cte(20)
      integer                    :: opt_cte(20) = 0

       NAMELIST /oneD / iQin,inTOout,name_oneD,cte,opt_cte

      character (len=*), parameter :: name_sub='Read_oneDTransfo'

      CALL alloc_oneDTransfo(oneDTransfo,nb_transfo)

      DO i=1,nb_transfo
        cte(:)    = ZERO
        cte(1)    = ONE
        inTOout   = .TRUE. ! T => no inversion (inTOout), F => inversion (outTOin)
        name_oneD = "identity"
        type_oneD = 0 ! identity
        iQin      = 0
        opt_cte   = 0

        read(in_unitp,oneD,IOSTAT=err)
        IF (err /= 0) THEN
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) '  while reading the "oneD" namelist'
           write(out_unitp,*) ' end of file or end of record'
           write(out_unitp,*) ' Check your data !!'
           STOP
        END IF

        write(out_unitp,oneD)

        IF (iQin < 1 .OR. iQin > nb_Qin) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  iQin is out of range',iQin
          write(out_unitp,*) '  range: [1:',nb_Qin,']'
          write(out_unitp,*) ' Check your data !!'
          STOP
        END IF
        oneDTransfo(i)%iQin      = iQin
        oneDTransfo(i)%cte       = cte
        oneDTransfo(i)%opt_cte   = opt_cte
        oneDTransfo(i)%inTOout   = inTOout
        oneDTransfo(i)%name_oneD = name_oneD

        !special test when the name is not defined (just the number):
        read(name_oneD,*,IOSTAT=err) type_oneD
        IF (err == 0) THEN
          write(out_unitp,*) ' name_oneD is a number',name_oneD
          oneDTransfo(i)%type_oneD = type_oneD
        ELSE


        SELECT CASE (name_oneD)
        CASE ('identity')
          oneDTransfo(i)%type_oneD = 0
        CASE ('affine')
          IF ( abs(cte(1)) < ONETENTH**4 .OR. abs(cte(1)) > TEN**4) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) '  cte(1) is too small or too large:',cte(1)
            write(out_unitp,*) ' for an affine transformation:'
            write(out_unitp,*) ' Qout = cte(1) * Qin + cte(2) or'
            write(out_unitp,*) ' Qin  = ( Qold - cte(2) ) / cte(1)'
            write(out_unitp,*) ' 10**-4 < abs(cte(1)) < 10**4'
            write(out_unitp,*) ' Check your data !!'
            STOP
          END IF
          IF (inTOout) THEN
            oneDTransfo(i)%type_oneD = 100
          ELSE
            oneDTransfo(i)%type_oneD = -100
          END IF
        CASE ('cos')
          IF (inTOout) THEN
            oneDTransfo(i)%type_oneD = 2
          ELSE
            oneDTransfo(i)%type_oneD = 5
          END IF
        CASE ('acos')
          IF (inTOout) THEN
            oneDTransfo(i)%type_oneD = 5
          ELSE
            oneDTransfo(i)%type_oneD = 2
          END IF
        CASE ('thetaTOx')
          IF ( cte(1) > ONE .OR. cte(1) < ZERO) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) '  For this transformation: ',trim(name_oneD)
            write(out_unitp,*) '    wrong value of cte(1):',cte(1)
            write(out_unitp,*) '    It has to be:   0 < cte(1) <= 1'
            write(out_unitp,*) ' Check your data !!'
            STOP
          END IF

          IF (inTOout) THEN
            ! t(x) = tan((x-Pi/2)/c1)  x E ]0,Pi[
            oneDTransfo(i)%type_oneD = -71
          ELSE
            ! t(x) = Pi/2 + c1*Atan(x) x E ]-inf,inf[
            oneDTransfo(i)%type_oneD = 71
          END IF
        CASE ('xTOtheta')
          IF ( cte(1) > ONE .OR. cte(1) < ZERO) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) '  For this transformation: ',trim(name_oneD)
            write(out_unitp,*) '    wrong value of cte(1):',cte(1)
            write(out_unitp,*) '    It has to be:   0 < cte(1) <= 1'
            write(out_unitp,*) ' Check your data !!'
            STOP
          END IF

          IF (inTOout) THEN
            oneDTransfo(i)%type_oneD = 71   ! x => theta
          ELSE
            oneDTransfo(i)%type_oneD = -71  ! theta => x
          END IF
        CASE ('xTOR')
         !  transfo R ]0,inf[ => x ]-inf,inf[
         ! 111      =>    (-a + x^2)/x = x-a/x x E ]0,inf[
         !-111      =>    1/2(x+sqrt(4a+x^2)) x E ]-inf,inf[
         ! a = cte(1)^2
          IF ( cte(1) == ZERO) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) '  For this transformation: ',trim(name_oneD)
            write(out_unitp,*) '    wrong value of cte(1):',cte(1)
            write(out_unitp,*) '    It has to be:  cte(1) /= 0'
            write(out_unitp,*) ' Check your data !!'
            STOP
          END IF

          IF (inTOout) THEN
            oneDTransfo(i)%type_oneD = 111   ! x => R
          ELSE
            oneDTransfo(i)%type_oneD = -111  ! R => x
          END IF

        CASE ('xTOconstX','xTOu')
         ! invers of R0.tanh(x/R0) x E ]-inf,inf[  (invers)
         ! t(x) = R0 atanh(x/R0) R0=cte(1)
          IF ( cte(1) == ZERO) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) '  For this transformation: ',trim(name_oneD)
            write(out_unitp,*) '    wrong value of cte(1):',cte(1)
            write(out_unitp,*) '    It has to be:  cte(1) /= 0'
            write(out_unitp,*) ' Check your data !!'
            STOP
          END IF

          IF (inTOout) THEN
            oneDTransfo(i)%type_oneD = 74   ! x => R
          ELSE
            oneDTransfo(i)%type_oneD = -74  ! R => x
          END IF

        CASE default ! ERROR: wrong transformation !
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' The oneD transformation is UNKNOWN: ',     &
                                                    trim(name_oneD)
          write(out_unitp,*) ' Check your data !!'
          STOP
        END SELECT
        END IF
      END DO

      END SUBROUTINE Read_oneDTransfo

      SUBROUTINE Write_oneDTransfo(oneDTransfo)

      TYPE (Type_oneDTransfo), pointer, intent(in) :: oneDTransfo(:)

      integer :: it
      character (len=*), parameter :: name_sub='Write_oneDTransfo'

      IF (.NOT. associated(oneDTransfo)) RETURN

      write(out_unitp,*) 'BEGINNING Write_oneDTransfo ',size(oneDTransfo)
      DO it=1,size(oneDTransfo)
        write(out_unitp,*) 'it,iQin       ',it,oneDTransfo(it)%iQin
        write(out_unitp,*) 'it,type_oneD  ',it,oneDTransfo(it)%type_oneD
        write(out_unitp,*) 'it,inTOout    ',it,oneDTransfo(it)%inTOout
        write(out_unitp,*) 'it,name_oneD  ',it,trim(adjustl(oneDTransfo(it)%name_oneD))
        write(out_unitp,*) 'it,cte(:)     ',it,oneDTransfo(it)%cte(:)
        write(out_unitp,*) 'it,opt_cte(:) ',it,oneDTransfo(it)%opt_cte(:)
      END DO
      write(out_unitp,*) 'END Write_oneDTransfo '

      END SUBROUTINE Write_oneDTransfo


      !!@description: TODO
      !!@param: TODO
      SUBROUTINE calc_oneDTransfo(dnQin,dnQout,oneDTransfo,nderiv,inTOout)

        TYPE (Type_dnVec), intent(inout)    :: dnQin,dnQout
        TYPE (Type_oneDTransfo), intent(in) :: oneDTransfo(:)
        integer, intent(in)                 :: nderiv
        logical                             :: inTOout


        integer                    :: iQin,type_oneD
        character (len=Name_len)   :: name_oneD
        real (kind=Rkind)          :: cte(20)

        TYPE (Type_dnS)   :: dnR,dntR

        integer :: nb_act

        integer :: i,j,k

        integer :: iQ,it=0

!----- for debuging ----------------------------------
       character (len=*),parameter :: name_sub='calc_oneDTransfo'
       logical, parameter :: debug=.FALSE.
!       logical, parameter :: debug=.TRUE.
!----- for debuging ----------------------------------


!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
      END IF
!---------------------------------------------------------------------

      CALL check_alloc_dnVec(dnQin,'dnQin',name_sub)
      CALL check_alloc_dnVec(dnQout,'dnQout',name_sub)

      nb_act = dnQin%nb_var_deriv
      CALL alloc_dnSVM(dnR,nb_act,nderiv)
      CALL alloc_dnSVM(dntR,nb_act,nderiv)

      IF (inTOout) THEN   ! => Qout=oneT(Qin)
        CALL sub_dnVec1_TO_dnVec2(dnQin,dnQout,nderiv=nderiv)

        DO i=1,size(oneDTransfo)
          iQin   = oneDTransfo(i)%iQin
          name_oneD = oneDTransfo(i)%name_oneD
          type_oneD = oneDTransfo(i)%type_oneD

          CALL sub_dnVec_TO_dnS(dnQin,dnR,iQin,nderiv)

          CALL sub_dnS1_TO_dntR2(dnR,dntR,type_oneD,nderiv,             &
                                 oneDTransfo(i)%cte)
          IF (debug) THEN
            write(out_unitp,*) 'i,iQin,type_oneD',i,iQin,type_oneD
            write(out_unitp,*) 'dnR'
            CALL Write_dnS(dnR)
            write(out_unitp,*) 'dntR'
            CALL Write_dnS(dntR)
          END IF

          CALL sub_dnS_TO_dnVec(dntR,dnQout,iQin,nderiv)

        END DO
      ELSE  ! => Qin=oneT^-1(Qout)
        CALL sub_dnVec1_TO_dnVec2(dnQout,dnQin,nderiv=nderiv)
        DO i=1,size(oneDTransfo)
          iQin   = oneDTransfo(i)%iQin
          name_oneD = oneDTransfo(i)%name_oneD
          type_oneD = -oneDTransfo(i)%type_oneD ! the invers

          CALL sub_dnVec_TO_dnS(dnQout,dnR,iQin,nderiv)

          CALL sub_dnS1_TO_dntR2(dnR,dntR,type_oneD,nderiv,             &
                                 oneDTransfo(i)%cte)
          IF (debug) THEN
            write(out_unitp,*) 'i,iQin,type_oneD',i,iQin,type_oneD
            write(out_unitp,*) 'dnR'
            CALL Write_dnS(dnR)
            write(out_unitp,*) 'dntR'
            CALL Write_dnS(dntR)
          END IF

          CALL sub_dnS_TO_dnVec(dntR,dnQin,iQin,nderiv)

        END DO
      END IF
      CALL dealloc_dnSVM(dnR)
      CALL dealloc_dnSVM(dntR)
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------
      END SUBROUTINE calc_oneDTransfo

      END MODULE mod_OneDTransfo
