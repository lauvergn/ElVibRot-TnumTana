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
      MODULE mod_TwoDTransfo
      use mod_system
      USE mod_dnSVM
      IMPLICIT NONE

      PRIVATE

      TYPE Type_TwoDTransfo
        integer                  :: Type_2D ! 0: cart: identity; 1: polar; 2: spherical, 3: R, x,z
        character (len=Name_len) :: name_Transfo_2D = ''
        integer                  :: list_TwoD_coord(2) = [0,0]
        real (kind=Rkind)        :: Twod0 ! for the zundel 2d transformation
        real (kind=Rkind)        :: theta,phi ! for the spherical transformation

      END TYPE Type_TwoDTransfo

      INTERFACE alloc_array
        MODULE PROCEDURE alloc_array_OF_TwoDTransfodim1
      END INTERFACE
      INTERFACE dealloc_array
        MODULE PROCEDURE dealloc_array_OF_TwoDTransfodim1
      END INTERFACE

      PUBLIC :: Type_TwoDTransfo, alloc_TwoDTransfo, dealloc_TwoDTransfo
      PUBLIC :: Read_TwoDTransfo, Write_TwoDTransfo, calc_TwoDTransfo
      PUBLIC :: TwoDTransfo1TOTwoDTransfo2

      CONTAINS

      SUBROUTINE alloc_TwoDTransfo(TwoDTransfo,nb_transfo)
      TYPE (Type_TwoDTransfo), pointer, intent(inout) :: TwoDTransfo(:)
      integer,                          intent(in)    :: nb_transfo

      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='alloc_TwoDTransfo'

      IF (associated(TwoDTransfo)) THEN
        CALL dealloc_TwoDTransfo(TwoDTransfo)
      END IF
      IF (nb_transfo < 1) RETURN

      CALL alloc_array(TwoDTransfo,[nb_transfo],"TwoDTransfo",name_sub)

      END SUBROUTINE alloc_TwoDTransfo

      SUBROUTINE dealloc_TwoDTransfo(TwoDTransfo)

      TYPE (Type_TwoDTransfo), pointer, intent(inout) :: TwoDTransfo(:)

      character (len=*), parameter :: name_sub='dealloc_TwoDTransfo'

      IF (.NOT. associated(TwoDTransfo)) RETURN


      CALL dealloc_array(TwoDTransfo,"TwoDTransfo",name_sub)

      END SUBROUTINE dealloc_TwoDTransfo

      SUBROUTINE alloc_array_OF_TwoDTransfodim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      TYPE (Type_TwoDTransfo), pointer, intent(inout) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=1
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_TwoDTransfodim1'
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
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Type_TwoDTransfo')

      END SUBROUTINE alloc_array_OF_TwoDTransfodim1
      SUBROUTINE dealloc_array_OF_TwoDTransfodim1(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (Type_TwoDTransfo), pointer, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_TwoDTransfodim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN
       IF (.NOT. associated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Type_TwoDTransfo')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_TwoDTransfodim1

      SUBROUTINE Read_TwoDTransfo(TwoDTransfo,nb_transfo,nb_Qin)

      TYPE (Type_TwoDTransfo), pointer, intent(inout) :: TwoDTransfo(:)
      integer,                          intent(in)    :: nb_transfo,nb_Qin

      integer           :: nb_coord,Type_2D,list_TwoD_coord(2)
      real (kind=Rkind) :: d0,theta,phi
      logical           :: multiple

       NAMELIST / TwoD / Type_2D,list_TwoD_coord,d0,theta,phi


      integer :: i,err_io,err_mem,memory
      character (len=*), parameter :: name_sub='Read_TwoDTransfo'

      CALL alloc_TwoDTransfo(TwoDTransfo,nb_transfo)

      DO i=1,nb_transfo

        Type_2D            = 0
        list_TwoD_coord(:) = 0
        d0                 = 1.5d0 ! in bohr
        ! spherical angle
        theta              = PI/TWO
        phi                = 0
        read(in_unitp,TwoD,IOSTAT=err_io)
        IF (err_io /= 0) THEN
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) '  while reading "TwoD" namelist'
           write(out_unitp,*) ' end of file or end of record'
           write(out_unitp,*) ' Check your data !!'
           STOP
        END IF
        TwoDTransfo(i)%Type_2D         = Type_2D
        TwoDTransfo(i)%list_TwoD_coord = list_TwoD_coord
        TwoDTransfo(i)%Twod0           = d0+d0
        TwoDTransfo(i)%theta           = theta
        TwoDTransfo(i)%phi             = phi

        nb_coord = count(list_TwoD_coord /= 0)

        IF (nb_coord /= 2) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' The number of coordinates is different from 2'
          write(out_unitp,*) ' Check your data !!'
          write(out_unitp,*) 'list_TwoD_coord: ',list_TwoD_coord(:)
          STOP
        END IF

        multiple = (count(list_TwoD_coord == list_TwoD_coord(1)) /= 1) .OR.       &
                   (count(list_TwoD_coord == list_TwoD_coord(2)) /= 1)
        IF (multiple) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' Some values are identical in the list_TwoD_coord'
          write(out_unitp,*) 'list_TwoD_coord: ',list_TwoD_coord(:)
          write(out_unitp,*) ' Check your data !!'
          STOP
        END IF

        SELECT CASE (Type_2D)
        CASE (0) ! identity
          TwoDTransfo(i)%name_Transfo_2D = 'identity: x,y'
        CASE (1) ! polar
          TwoDTransfo(i)%name_Transfo_2D = 'Polar: R,theta'
        CASE (2) ! special Zundel
          TwoDTransfo(i)%name_Transfo_2D = 'Zundel: z,R'
        CASE (3) ! spherical
          TwoDTransfo(i)%name_Transfo_2D = 'Spherical'
        CASE (-3) ! spherical
          TwoDTransfo(i)%name_Transfo_2D = 'Spherical'
        CASE default ! ERROR: wrong transformation !
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' The type of TwoD transformation is UNKNOWN: ',Type_2D
          write(out_unitp,*) ' The possible values are:'
          write(out_unitp,*) '    0: identity: x,y'
          write(out_unitp,*) '    1: Polar: R,theta'
          write(out_unitp,*) '    2: Zundel: z,R'
          write(out_unitp,*) '    3: Spherical: rotation of the spherical angles (theta,phi)'
          write(out_unitp,*) '   -3: Spherical: rotation of the spherical angles (u,phi)'

          write(out_unitp,*) ' Check your data !!'
          STOP

        END SELECT

      END DO
      CALL Write_TwoDTransfo(TwoDTransfo)

      END SUBROUTINE Read_TwoDTransfo

      SUBROUTINE Write_TwoDTransfo(TwoDTransfo)

      TYPE (Type_TwoDTransfo), pointer, intent(in) :: TwoDTransfo(:)

      integer :: i,err_mem,memory
      character (len=*), parameter :: name_sub='Write_TwoDTransfo'

      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) 'asso TwoDTransfo: ',associated(TwoDTransfo)
      write(out_unitp,*) 'nb_transfo:       ',size(TwoDTransfo)

      IF (associated(TwoDTransfo)) THEN
      DO i=1,size(TwoDTransfo)
        write(out_unitp,*) 'Type_2D:           ',TwoDTransfo(i)%Type_2D
        write(out_unitp,*) 'name_Transfo_2D:   ',TwoDTransfo(i)%name_Transfo_2D
        write(out_unitp,*) 'list_TwoD_coord:   ',TwoDTransfo(i)%list_TwoD_coord
        write(out_unitp,*) 'Twod0:             ',TwoDTransfo(i)%Twod0
        write(out_unitp,*) ' => d0:            ',TwoDTransfo(i)%Twod0*HALF
        write(out_unitp,*) 'theta,phi:         ',TwoDTransfo(i)%theta,TwoDTransfo(i)%phi

      END DO
      END IF
      write(out_unitp,*) 'END ',name_sub

      END SUBROUTINE Write_TwoDTransfo

      SUBROUTINE calc_TwoDTransfo(dnQin,dnQout,TwoDTransfo,nderiv,inTOout)
        USE mod_QML_dnS
        TYPE (Type_dnVec), intent(inout)              :: dnQin,dnQout
        TYPE (Type_TwoDTransfo),pointer, intent(in)   :: TwoDTransfo(:)
        integer, intent(in)                           :: nderiv
        logical, intent(in)                           :: inTOout


        TYPE (Type_dnS) :: dnR,dnZ,dnZp,dntR
        TYPE (dnS_t)    :: dntho,dnctho,dnphio,  dnthn,dncthn,dnphin
        TYPE (dnS_t)    :: dnXo,dnYo,dnZo,dnXn,dnYn,dnZn

        real(kind=Rkind) :: cte(20)

        integer :: i,i1,i2
        integer :: dnErr

!----- for debuging ----------------------------------
       character (len=*),parameter :: name_sub='calc_TwoDTransfo'
       logical, parameter :: debug=.FALSE.
!       logical, parameter :: debug=.TRUE.
!----- for debuging ----------------------------------


!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'dnQin'
        CALL Write_dnSVM(dnQin,nderiv)
      END IF
!---------------------------------------------------------------------

      IF (.NOT. associated(TwoDTransfo)) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' TwoDTransfo is NOT associated'
        write(out_unitp,*) ' Check source !!'
        STOP
      END IF

      CALL check_alloc_dnVec(dnQin,'dnQin',name_sub)
      CALL check_alloc_dnVec(dnQout,'dnQout',name_sub)

      IF (inTOout) THEN
        CALL sub_dnVec1_TO_dnVec2(dnQin,dnQout,nderiv)

        DO i=1,size(TwoDTransfo)
          SELECT CASE (TwoDTransfo(i)%Type_2D)
          CASE (0) ! identity
            ! nothing
          CASE (1) ! polar
            STOP 'polar not yet'
          CASE (2) ! Zundel: z,R => z'= z/(R-2d0) and R'=R
                  ! {z',R'} => {z,R}: z=z'*(R'-2d0) and R=R'

            i1 = TwoDTransfo(i)%list_TwoD_coord(1)
            i2 = TwoDTransfo(i)%list_TwoD_coord(2)

            !write(out_unitp,*) 'i1,i2',i1,i2
            CALL sub_dnVec_TO_dnS(dnQin,dnZp,i1,nderiv)
            CALL sub_dnVec_TO_dnS(dnQin,dnR,i2,nderiv)

            ! 100 (affine) =>    cte(1) * x + cte(2): => tR=R-2d0
            cte(:) = ZERO
            cte(1:2) = [ONE,-TwoDTransfo(i)%Twod0]
            CALL sub_dnS1_TO_dntR2(dnR,dntR,100,nderiv,cte,dnErr)

            CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnZp,dntR,dnZ,nderiv)

            ! transfert dnZp in dnQout
            CALL sub_dnS_TO_dnVec(dnZ,dnQout,i1,nderiv)

            CALL dealloc_dnSVM(dnR)
            CALL dealloc_dnSVM(dntR)
            CALL dealloc_dnSVM(dnZ)
            CALL dealloc_dnSVM(dnZp)

          CASE (3) ! spherical
            i1 = TwoDTransfo(i)%list_TwoD_coord(1)
            i2 = TwoDTransfo(i)%list_TwoD_coord(2)

            !write(out_unitp,*) 'i1,i2',i1,i2
            CALL sub_dnVec_TO_dnSt(dnQin,dnthn, i1)
            CALL sub_dnVec_TO_dnSt(dnQin,dnphin,i2)

            dnXn = cos(dnphin) * sin(dnthn)
            dnYn = sin(dnphin) * sin(dnthn)
            dnZn = cos(dnthn)

            dnXo = cos(TwoDTransfo(i)%theta) * dnXn - sin(TwoDTransfo(i)%theta) * dnZn
            dnYo = dnYn
            dnZo = sin(TwoDTransfo(i)%theta) * dnXn + cos(TwoDTransfo(i)%theta) * dnZn

            dntho  = acos(dnZo)
            dnphio = atan2(dnYo,dnXo)

           ! transfert in dnQout
           CALL sub_dnSt_TO_dnVec(dntho,dnQout,i1)
           CALL sub_dnSt_TO_dnVec(dnphio,dnQout,i2)

           CALL QML_dealloc_dnS(dnthn)
           CALL QML_dealloc_dnS(dnphin)
           CALL QML_dealloc_dnS(dnXn)
           CALL QML_dealloc_dnS(dnYn)
           CALL QML_dealloc_dnS(dnZn)
           CALL QML_dealloc_dnS(dntho)
           CALL QML_dealloc_dnS(dnphio)
           CALL QML_dealloc_dnS(dnXo)
           CALL QML_dealloc_dnS(dnYo)
           CALL QML_dealloc_dnS(dnZo)
         CASE (-3) ! spherical
            i1 = TwoDTransfo(i)%list_TwoD_coord(1)
            i2 = TwoDTransfo(i)%list_TwoD_coord(2)

            !write(out_unitp,*) 'i1,i2',i1,i2
            CALL sub_dnVec_TO_dnSt(dnQin,dncthn, i1)
            CALL sub_dnVec_TO_dnSt(dnQin,dnphin,i2)

            dnXn = cos(dnphin) * sqrt(ONE-dncthn*dncthn)
            dnYn = sin(dnphin) * sqrt(ONE-dncthn*dncthn)
            dnZn = dncthn

            dnXo = cos(TwoDTransfo(i)%theta) * dnXn - sin(TwoDTransfo(i)%theta) * dnZn
            dnYo = dnYn
            dnZo = sin(TwoDTransfo(i)%theta) * dnXn + cos(TwoDTransfo(i)%theta) * dnZn

            dntho  = dnZo
            dnphio = atan2(dnYo,dnXo)

           ! transfert in dnQout
           CALL sub_dnSt_TO_dnVec(dntho,dnQout,i1)
           CALL sub_dnSt_TO_dnVec(dnphio,dnQout,i2)

           CALL QML_dealloc_dnS(dncthn)
           CALL QML_dealloc_dnS(dnphin)
           CALL QML_dealloc_dnS(dnXn)
           CALL QML_dealloc_dnS(dnYn)
           CALL QML_dealloc_dnS(dnZn)
           CALL QML_dealloc_dnS(dnctho)
           CALL QML_dealloc_dnS(dnphio)
           CALL QML_dealloc_dnS(dnXo)
           CALL QML_dealloc_dnS(dnYo)
           CALL QML_dealloc_dnS(dnZo)

          CASE default ! ERROR: wrong transformation !
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' The type of TwoD transformation is UNKNOWN: ',TwoDTransfo%Type_2D
            write(out_unitp,*) ' The possible values are:'
            write(out_unitp,*) '    0: identity: x,y'
            write(out_unitp,*) '    1: Polar: R,theta'
            write(out_unitp,*) '    2: Zundel: z,R'
            write(out_unitp,*) '    3: Spherical: rotation of the spherical angles'
            write(out_unitp,*) ' Check your data !!'
            STOP

          END SELECT
        END DO
      ELSE
        CALL sub_dnVec1_TO_dnVec2(dnQout,dnQin,nderiv)
        DO i=1,size(TwoDTransfo)

          SELECT CASE (TwoDTransfo(i)%Type_2D)
          CASE (0) ! identity
            ! nothing
          CASE (1) ! polar
            STOP 'polar not yet'
          CASE (2) ! Zundel: z,R => z'= z/(R-2d0) and R'=R

            i1 = TwoDTransfo(i)%list_TwoD_coord(1)
            i2 = TwoDTransfo(i)%list_TwoD_coord(2)

            !write(out_unitp,*) 'i1,i2',i1,i2
            CALL sub_dnVec_TO_dnS(dnQout,dnz,i1,nderiv)
            CALL sub_dnVec_TO_dnS(dnQout,dnR,i2,nderiv)

            ! 100 (affine) =>    cte(1) * x + cte(2): => tR=R-2d0
            cte(:) = ZERO
            cte(1:2) = [ONE,-TwoDTransfo(i)%Twod0]
            CALL sub_dnS1_TO_dntR2(dnR,dntR,100,nderiv,cte,dnErr)

            ! 90 (-90) =>    1/x: => R = 1/tR = 1/(R-2d0)
            CALL sub_dnS1_TO_dntR2(dntR,dnR,90,nderiv,cte,dnErr)

            CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnZ,dnR,dnZp,nderiv)

            ! transfert dnZp in dnQout
            CALL sub_dnS_TO_dnVec(dnZp,dnQin,i1,nderiv)

            CALL dealloc_dnSVM(dnR)
            CALL dealloc_dnSVM(dntR)
            CALL dealloc_dnSVM(dnZ)
            CALL dealloc_dnSVM(dnZp)

          CASE (3) ! spherical
            i1 = TwoDTransfo(i)%list_TwoD_coord(1)
            i2 = TwoDTransfo(i)%list_TwoD_coord(2)

            !write(out_unitp,*) 'i1,i2',i1,i2
            CALL sub_dnVec_TO_dnSt(dnQout,dntho, i1)
            CALL sub_dnVec_TO_dnSt(dnQout,dnphio,i2)

            dnXo = cos(dnphio) * sin(dntho)
            dnYo = sin(dnphio) * sin(dntho)
            dnZo = cos(dntho)

            dnXn = cos(-TwoDTransfo(i)%theta) * dnXo - sin(-TwoDTransfo(i)%theta) * dnZo
            dnYn = dnYo
            dnZn = sin(-TwoDTransfo(i)%theta) * dnXo + cos(-TwoDTransfo(i)%theta) * dnZo

            dnthn  = acos(dnZn)
            dnphin = atan2(dnYn,dnXn)

           ! transfert in dnQout
           CALL sub_dnSt_TO_dnVec(dnthn,dnQin,i1)
           CALL sub_dnSt_TO_dnVec(dnphin,dnQin,i2)

           CALL QML_dealloc_dnS(dnthn)
           CALL QML_dealloc_dnS(dnphin)
           CALL QML_dealloc_dnS(dnXn)
           CALL QML_dealloc_dnS(dnYn)
           CALL QML_dealloc_dnS(dnZn)
           CALL QML_dealloc_dnS(dntho)
           CALL QML_dealloc_dnS(dnphio)
           CALL QML_dealloc_dnS(dnXo)
           CALL QML_dealloc_dnS(dnYo)
           CALL QML_dealloc_dnS(dnZo)
         CASE (-3) ! spherical
            i1 = TwoDTransfo(i)%list_TwoD_coord(1)
            i2 = TwoDTransfo(i)%list_TwoD_coord(2)

            !write(out_unitp,*) 'i1,i2',i1,i2
            CALL sub_dnVec_TO_dnSt(dnQout,dnctho, i1)
            CALL sub_dnVec_TO_dnSt(dnQout,dnphio,i2)

            dnXo = cos(dnphio) * sqrt(ONE-dnctho*dnctho)
            dnYo = sin(dnphio) * sqrt(ONE-dnctho*dnctho)
            dnZo = dnctho

            dnXn = cos(-TwoDTransfo(i)%theta) * dnXo - sin(-TwoDTransfo(i)%theta) * dnZo
            dnYn = dnYo
            dnZn = sin(-TwoDTransfo(i)%theta) * dnXo + cos(-TwoDTransfo(i)%theta) * dnZo

            dnthn  = dnZn
            dnphin = atan2(dnYn,dnXn)

           ! transfert in dnQout
           CALL sub_dnSt_TO_dnVec(dnthn,dnQin,i1)
           CALL sub_dnSt_TO_dnVec(dnphin,dnQin,i2)

           CALL QML_dealloc_dnS(dncthn)
           CALL QML_dealloc_dnS(dnphin)
           CALL QML_dealloc_dnS(dnXn)
           CALL QML_dealloc_dnS(dnYn)
           CALL QML_dealloc_dnS(dnZn)
           CALL QML_dealloc_dnS(dnctho)
           CALL QML_dealloc_dnS(dnphio)
           CALL QML_dealloc_dnS(dnXo)
           CALL QML_dealloc_dnS(dnYo)
           CALL QML_dealloc_dnS(dnZo)
          CASE default ! ERROR: wrong transformation !
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' The type of TwoD transformation is UNKNOWN: ',TwoDTransfo%Type_2D
            write(out_unitp,*) ' The possible values are:'
            write(out_unitp,*) '    0: identity: x,y'
            write(out_unitp,*) '    1: Polar: R,theta'
            write(out_unitp,*) '    2: Zundel: z,R'
            write(out_unitp,*) '    3: Spherical: rotation of the spherical angles'
            write(out_unitp,*) ' Check your data !!'
            STOP

          END SELECT
        END DO
      END IF


!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------
      END SUBROUTINE calc_TwoDTransfo

      SUBROUTINE TwoDTransfo1TOTwoDTransfo2(TwoDTransfo1,TwoDTransfo2)
        TYPE (Type_TwoDTransfo),pointer, intent(in)    :: TwoDTransfo1(:)
        TYPE (Type_TwoDTransfo),pointer, intent(inout) :: TwoDTransfo2(:)

        integer :: i

        IF (.NOT. associated(TwoDTransfo1)) RETURN

        CALL alloc_TwoDTransfo(TwoDTransfo2,nb_transfo=size(TwoDTransfo1))
        DO i=1,size(TwoDTransfo1)
          TwoDTransfo2(i)%Type_2D            = TwoDTransfo1(i)%Type_2D
          TwoDTransfo2(i)%name_Transfo_2D    = TwoDTransfo1(i)%name_Transfo_2D
          TwoDTransfo2(i)%list_TwoD_coord(:) = TwoDTransfo1(i)%list_TwoD_coord(:)
          TwoDTransfo2(i)%Twod0              = TwoDTransfo1(i)%Twod0
          TwoDTransfo2(i)%theta              = TwoDTransfo1(i)%theta
          TwoDTransfo2(i)%phi                = TwoDTransfo1(i)%phi

        END DO

      END SUBROUTINE TwoDTransfo1TOTwoDTransfo2

      END MODULE mod_TwoDTransfo
