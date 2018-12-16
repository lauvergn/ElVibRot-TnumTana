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
!      with contributions of Mamadou Ndong
!
!===========================================================================
!===========================================================================
MODULE mod_CartesianTransfo
      use mod_system
      use mod_dnSVM
      use mod_Lib_QTransfo, only: write_dnx, calc_cross_product
      IMPLICIT NONE

      PRIVATE

      !!@description: enables to change the BF frame orientation (Eckart or ....)
      !!@param: TODO
      TYPE Type_CartesianTransfo

        logical           :: New_Orient          = .FALSE. ! (F) T => Can use different orientation for the z-matrix
        real (kind=Rkind) :: vAt1(3)             = ZERO
        real (kind=Rkind) :: vAt2(3)             = ZERO
        real (kind=Rkind) :: vAt3(3)             = ZERO

        logical           :: P_Axis_Ref          = .FALSE.
        real (kind=Rkind) :: InertiaT(3,3)       = ZERO ! inertia tenssor
        real (kind=Rkind) :: ABC(3)              = ZERO ! rotational constants

        real (kind=Rkind) :: Rot_initial(3,3)    = ZERO      ! Initial Rotational matrix
                                                             ! The other ones (Rot_New_Orient and Rot_PA)
                                                             !    are not used anymore

        logical                    :: Eckart           = .FALSE.

        logical                    :: MultiRefEckart     = .FALSE.
        integer                    :: nb_Qact            = 0
        integer, pointer           :: list_Qact(:)       => null() ! list of active coordinates
        real (kind=Rkind), pointer :: tab_sc(:)          => null()
        real (kind=Rkind), pointer :: tab_phi0(:)        => null()
        integer, pointer           :: tab_expo(:)        => null()
        integer, pointer           :: tab_switch_type(:) => null()
        integer, pointer           :: tab_nb_Ref(:)      => null()

        integer                    :: MultiRef_type      = 0 ! 0: old periodic
                                                             ! 1: along coordinates (high barrier)
                                                             ! 2: along coordinates (low barrier)
                                                             ! 3: using distances (high barrier)
                                                             ! 4: using distances + Eckart (high barrier)

        logical           :: P_Axis_Always    = .FALSE.

        logical           :: ReadRefGeometry  = .FALSE.
        integer           :: nb_RefGeometry   = 0            ! nb of cart. coordinates
        integer           :: ncart_act        = 0            ! nb of cart. coordinates
        integer           :: nat_act          = 0            ! nb of cart. coordinates


        real (kind=Rkind), pointer :: Qxyz(:,:,:) => null()  ! Qxyz(3,nat_act,nb_RefGeometry) cart. coordinates
        real (kind=Rkind), pointer :: d0sm(:) => null()      ! sqrt of the masses

        integer           :: type_diago  = 4 ! enables to select diagonalization type (MatOFdnS)
        integer           :: type_cs     = 0 ! enables to select the way to calculate dnCos and dnSin
        logical           :: check_dnT   = .FALSE. ! if true, we check det(dnT) == 1
        real (kind=Rkind) :: dnTErr(0:3) = ZERO ! Error for the det(dnT)

      END TYPE Type_CartesianTransfo

      PUBLIC :: Type_CartesianTransfo, alloc_CartesianTransfo, dealloc_CartesianTransfo
      PUBLIC :: Read_CartesianTransfo, Write_CartesianTransfo
      PUBLIC :: CartesianTransfo1TOCartesianTransfo2, calc_CartesianTransfo_new
      PUBLIC :: P_Axis_CartesianTransfo
      PUBLIC :: calc_dnTxdnXin_TO_dnXout, calc_EckartRot, calc_dnTEckart, dnX_MultiRef
      PUBLIC :: centre_masse, sub3_dncentre_masse, sub3_NOdncentre_masse
      PUBLIC :: sub_dnxMassWeight, sub_dnxNOMassWeight

      CONTAINS

      SUBROUTINE alloc_CartesianTransfo(CartesianTransfo,ncart_act,     &
                                                 nb_RefGeometry,nb_Qact)

      TYPE (Type_CartesianTransfo), intent(inout) :: CartesianTransfo
      integer, intent(in) :: ncart_act
      integer, intent(in), optional :: nb_RefGeometry
      integer, intent(in), optional :: nb_Qact

      integer :: nb_RefGeometry_loc
      character (len=*), parameter :: name_sub='alloc_CartesianTransfo'

      IF (present(nb_RefGeometry)) THEN
        nb_RefGeometry_loc = nb_RefGeometry
      ELSE
        nb_RefGeometry_loc = 1
      END IF

      IF (associated(CartesianTransfo%Qxyz))                            &
        CALL dealloc_array(CartesianTransfo%Qxyz,                       &
                                       'CartesianTransfo%Qxyz',name_sub)

      IF (ncart_act > 0 .AND. nb_RefGeometry_loc > 0) THEN
        CartesianTransfo%ncart_act      = ncart_act
        CartesianTransfo%nat_act        = ncart_act/3
        CartesianTransfo%nb_RefGeometry = nb_RefGeometry_loc


        CALL alloc_array(CartesianTransfo%Qxyz,                         &
                                (/ 3,ncart_act/3,nb_RefGeometry_loc /), &
                                       'CartesianTransfo%Qxyz',name_sub)
        CartesianTransfo%Qxyz(:,:,:) = ZERO


        IF (.NOT. associated(CartesianTransfo%d0sm)) THEN
          CALL alloc_array(CartesianTransfo%d0sm,(/ ncart_act /),       &
                                       'CartesianTransfo%d0sm',name_sub)
          CartesianTransfo%d0sm(:) = ZERO
        END IF
      ELSE
        CartesianTransfo%ncart_act      = 0
        CartesianTransfo%nat_act        = 0
        CartesianTransfo%nb_RefGeometry = 0

        nullify(CartesianTransfo%Qxyz)
      END IF

      IF (present(nb_Qact)) THEN
         CartesianTransfo%nb_Qact = nb_Qact

         IF (associated(CartesianTransfo%list_Qact))                    &
         CALL dealloc_array(CartesianTransfo%list_Qact,                 &
                           'CartesianTransfo%list_Qact',name_sub)

         IF (associated(CartesianTransfo%tab_sc))                       &
         CALL dealloc_array(CartesianTransfo%tab_sc,                    &
                           'CartesianTransfo%tab_sc',name_sub)

         IF (associated(CartesianTransfo%tab_expo))                     &
         CALL dealloc_array(CartesianTransfo%tab_expo,                  &
                           'CartesianTransfo%tab_expo',name_sub)

         IF (associated(CartesianTransfo%tab_phi0))                     &
         CALL dealloc_array(CartesianTransfo%tab_phi0,                  &
                           'CartesianTransfo%tab_phi0',name_sub)

         IF (associated(CartesianTransfo%tab_switch_type))              &
         CALL dealloc_array(CartesianTransfo%tab_switch_type,           &
                           'CartesianTransfo%tab_switch_type',name_sub)

         IF (associated(CartesianTransfo%tab_nb_Ref))                   &
         CALL dealloc_array(CartesianTransfo%tab_nb_Ref,                &
                           'CartesianTransfo%tab_nb_Ref',name_sub)

        IF (nb_Qact > 0) THEN
          CALL alloc_array(CartesianTransfo%list_Qact,(/ nb_Qact /),    &
                          'CartesianTransfo%list_Qact',name_sub)
          CartesianTransfo%list_Qact(:) = 0

          CALL alloc_array(CartesianTransfo%tab_sc,(/ nb_Qact /),       &
                          'CartesianTransfo%tab_sc',name_sub)
          CartesianTransfo%tab_sc(:) = SIX

          CALL alloc_array(CartesianTransfo%tab_expo,(/ nb_Qact /),     &
                          'CartesianTransfo%tab_expo',name_sub)
          CartesianTransfo%tab_expo(:) = 0

          CALL alloc_array(CartesianTransfo%tab_phi0,(/ nb_Qact /),     &
                          'CartesianTransfo%tab_phi0',name_sub)
          CartesianTransfo%tab_phi0(:) = ZERO

          CALL alloc_array(CartesianTransfo%tab_switch_type,(/ nb_Qact /),&
                          'CartesianTransfo%tab_switch_type',name_sub)
          CartesianTransfo%tab_switch_type(:) = 0

          CALL alloc_array(CartesianTransfo%tab_nb_Ref,(/ nb_Qact /),   &
                          'CartesianTransfo%tab_nb_Ref',name_sub)
          CartesianTransfo%tab_nb_Ref(:) = 2

        END IF
      END IF

      END SUBROUTINE alloc_CartesianTransfo
      SUBROUTINE dealloc_CartesianTransfo(CartesianTransfo)

      TYPE (Type_CartesianTransfo), intent(inout) :: CartesianTransfo

      character (len=*), parameter :: name_sub='dealloc_CartesianTransfo'


      CartesianTransfo%New_Orient          = .FALSE.
      CartesianTransfo%vAt1(:)             = ZERO
      CartesianTransfo%vAt2(:)             = ZERO
      CartesianTransfo%vAt3(:)             = ZERO

      CartesianTransfo%P_Axis_Ref       = .FALSE.
      CartesianTransfo%P_Axis_Always    = .FALSE.
      CartesianTransfo%InertiaT(3,3)    = ZERO ! inertia tenssor
      CartesianTransfo%ABC(3)           = ZERO ! rotational constants

      CartesianTransfo%Rot_initial(3,3) = ZERO ! Initial Rotational matrix


      CartesianTransfo%Eckart           = .FALSE.
      CartesianTransfo%MultiRefEckart   = .FALSE.
      CartesianTransfo%MultiRef_type    = 0
      CartesianTransfo%nb_Qact          = 0

      CartesianTransfo%ReadRefGeometry  = .FALSE.
      CartesianTransfo%ncart_act      = 0
      CartesianTransfo%nat_act        = 0

      CartesianTransfo%nb_RefGeometry = 0


      IF (associated(CartesianTransfo%Qxyz))                            &
      CALL dealloc_array(CartesianTransfo%Qxyz,                         &
                        'CartesianTransfo%Qxyz',name_sub)

      IF (associated(CartesianTransfo%d0sm))                            &
      CALL dealloc_array(CartesianTransfo%d0sm,                         &
                        'CartesianTransfo%d0sm',name_sub)

      IF (associated(CartesianTransfo%list_Qact))                       &
      CALL dealloc_array(CartesianTransfo%list_Qact,                    &
                        'CartesianTransfo%list_Qact',name_sub)

      IF (associated(CartesianTransfo%tab_sc))                          &
      CALL dealloc_array(CartesianTransfo%tab_sc,                       &
                        'CartesianTransfo%tab_sc',name_sub)

         IF (associated(CartesianTransfo%tab_expo))                     &
         CALL dealloc_array(CartesianTransfo%tab_expo,                  &
                           'CartesianTransfo%tab_expo',name_sub)

      IF (associated(CartesianTransfo%tab_phi0))                        &
      CALL dealloc_array(CartesianTransfo%tab_phi0,                     &
                        'CartesianTransfo%tab_phi0',name_sub)

      IF (associated(CartesianTransfo%tab_switch_type))                 &
      CALL dealloc_array(CartesianTransfo%tab_switch_type,              &
                        'CartesianTransfo%tab_switch_type',name_sub)

      IF (associated(CartesianTransfo%tab_nb_Ref))                      &
      CALL dealloc_array(CartesianTransfo%tab_nb_Ref,                   &
                        'CartesianTransfo%tab_nb_Ref',name_sub)

      CartesianTransfo%type_diago      = 4
      CartesianTransfo%type_cs         = 0
      CartesianTransfo%check_dnT       = .FALSE.
      CartesianTransfo%dnTErr(:)       = ZERO


      END SUBROUTINE dealloc_CartesianTransfo
      SUBROUTINE Read_CartesianTransfo(CartesianTransfo)

      TYPE (Type_CartesianTransfo), intent(inout) :: CartesianTransfo


      logical           :: New_Orient
      real (kind=Rkind) :: vAt1(3)
      real (kind=Rkind) :: vAt2(3)
      real (kind=Rkind) :: vAt3(3)
      logical           :: Eckart,MultiRefEckart
      logical           :: P_Axis_Ref
      logical           :: P_Axis_Always
      logical           :: ReadRefGeometry
      integer           :: nat,nb_RefGeometry,MultiRef_type
      real (kind=Rkind) :: Mtot_inv
      real (kind=Rkind), pointer :: d0x(:)

      character (len=Name_len) :: QnameRead,name_int
      character (len=Name_len) :: unit
      character (len=1) :: chara

      integer :: i,j,nb_Qact
      integer, allocatable :: list_Qact(:)
      real (kind=Rkind) :: sc
      integer :: type_diago,type_cs
      logical :: check_dnT



      integer :: err_io
      character (len=*), parameter :: name_sub='Read_CartesianTransfo'

      NAMELIST /Cartesian/ New_Orient,Eckart,P_Axis_Ref,P_Axis_Always,  &
                           ReadRefGeometry,nat,unit,nb_RefGeometry,     &
                           vAt1,vAt2,vAt3,                              &
                           MultiRefEckart,MultiRef_type,nb_Qact,sc,     &
                           type_diago,type_cs,check_dnT

      write(out_unitp,*) 'BEGINNING ',name_sub

       nullify(d0x)

      New_Orient       = .FALSE. ! (F) T => Can use different orientation for the z-matrix
      vAt1(:)          = ZERO
      vAt2(:)          = ZERO
      vAt3(:)          = ZERO

      Eckart           = .FALSE.
      MultiRefEckart   = .FALSE.
      MultiRef_type    = 0
      nb_Qact          = 0
      P_Axis_Ref       = .FALSE.
      P_Axis_Always    = .FALSE.

      ReadRefGeometry  = .FALSE.
      nb_RefGeometry   = 1
      nat              = 0
      unit             = 'au'
      sc               = SIX
      type_diago       = 4
      type_cs          = 0
      check_dnT        = .FALSE.

      read(in_unitp,Cartesian,IOSTAT=err_io)
      IF (err_io < 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) '  while reading the namelist "Cartesian"'
        write(out_unitp,*) ' end of file or end of record'
        write(out_unitp,*) ' Probably you have forgotten the "Cartesian transfo"'
        write(out_unitp,*) ' Check your data !!'
        STOP
      END IF
      IF (err_io > 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) '  while reading the namelist "Cartesian"'
        write(out_unitp,*) ' Probably, some arguments of namelist are wrong.'
        write(out_unitp,*) ' Check your data !!'
        STOP
      END IF
      write(out_unitp,Cartesian)
      CALL flush_perso(out_unitp)

      CartesianTransfo%New_Orient      = New_Orient
      CartesianTransfo%vAt1            = vAt1
      CartesianTransfo%vAt2            = vAt2
      CartesianTransfo%vAt3            = vAt3
      CartesianTransfo%Eckart          = Eckart
      CartesianTransfo%MultiRefEckart  = MultiRefEckart
      CartesianTransfo%MultiRef_type   = MultiRef_type

      CartesianTransfo%nb_RefGeometry  = nb_RefGeometry
      CartesianTransfo%nb_Qact         = nb_Qact

      CartesianTransfo%P_Axis_Ref      = P_Axis_Ref
      CartesianTransfo%P_Axis_Always   = P_Axis_Always

      CartesianTransfo%ReadRefGeometry = ReadRefGeometry
      CartesianTransfo%nat_act         = nat
      CartesianTransfo%ncart_act       = 3*nat

      CartesianTransfo%type_diago      = type_diago
      CartesianTransfo%type_cs         = type_cs
      CartesianTransfo%check_dnT       = check_dnT


      IF (ReadRefGeometry .AND. nat < 1) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' ReadRefGeometry=t and nat < 1'
        write(out_unitp,*) ' Parameters incompatible, check the Fortran!'
        STOP
      END IF


      IF (MultiRefEckart) THEN

        SELECT CASE (MultiRef_type)
        CASE (0) ! original periodic (do not use)
          CALL Read_MultiRef_type0_OF_CartesianTransfo(CartesianTransfo)
        CASE (1) ! along coordinates (high barrier)
          CALL Read_MultiRef_type1_OF_CartesianTransfo(CartesianTransfo)
        CASE (2) ! along coordinates (low barrier)
          !CALL Read_MultiRef_type2_OF_CartesianTransfo(CartesianTransfo)
          STOP 'case2 not yet'
        CASE (3) ! using distances (high barrier)
          nb_Qact = 1
          CartesianTransfo%nb_Qact         = nb_Qact
          CALL alloc_CartesianTransfo(CartesianTransfo,3*nat,           &
                                                 nb_RefGeometry,nb_Qact)
          CartesianTransfo%tab_sc(:)       = sc
        CASE Default ! using distances (high barrier)
          nb_Qact = 1
          CartesianTransfo%nb_Qact         = nb_Qact
          CALL alloc_CartesianTransfo(CartesianTransfo,3*nat,           &
                                                 nb_RefGeometry,nb_Qact)
          CartesianTransfo%tab_sc(:)       = sc
        END SELECT
      ELSE
        CALL alloc_CartesianTransfo(CartesianTransfo,                   &
                                            CartesianTransfo%ncart_act, &
                                        CartesianTransfo%nb_RefGeometry)
      END IF
      !CALL Write_CartesianTransfo(CartesianTransfo)


      CALL mat_id(CartesianTransfo%Rot_initial,3,3)


      IF (ReadRefGeometry) THEN
        CALL alloc_array(d0x,(/ 3*nat /),'d0x','Read_CartesianTransfo')

        Mtot_inv = THREE / sum(CartesianTransfo%d0sm(:)**2)

        DO j=1,CartesianTransfo%nb_RefGeometry
          DO i=1,CartesianTransfo%nat_act
            read(in_unitp,*,IOSTAT=err_io) QnameRead,CartesianTransfo%Qxyz(:,i,j)
            IF (err_io /= 0) THEN
              write(out_unitp,*) ' ERROR in ',name_sub
              write(out_unitp,*) '  while reading a Cartessian reference geometry ...'
              write(out_unitp,*) '   ... just after the namelist "Cartesian".'
              write(out_unitp,'(a,i0,a,i0,a,i0,a)') '  Trying to read the atom:',i, &
                                      ' among ',CartesianTransfo%nat_act,&
                                      ' of the reference geometry:',j,'.'
              write(out_unitp,*) ' Check your data !!'
              STOP
            END IF

          END DO
          read(in_unitp,*,IOSTAT=err_io)
          IF (err_io /= 0) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) '  while reading an empty line ...'
            write(out_unitp,*) '   ... just after Cartessian reference ',&
                             'geometry of the namelist "Cartesian".'
            write(out_unitp,*) ' Probably, end-of-file.'
            write(out_unitp,*) ' Check your data !!'
            STOP
          END IF


          d0x(:) = reshape(CartesianTransfo%Qxyz(:,:,j),(/CartesianTransfo%ncart_act/) )

          CALL centre_masse(CartesianTransfo%ncart_act,                 &
                            CartesianTransfo%ncart_act,                 &
                            d0x,CartesianTransfo%d0sm**2,Mtot_inv)

          CartesianTransfo%Qxyz(:,:,j) = reshape(d0x,(/3,nat/) )

        END DO

        CALL string_uppercase_TO_lowercase(unit)
        IF (unit == 'angs' ) THEN
          CartesianTransfo%Qxyz = CartesianTransfo%Qxyz / 0.52917720835354106_Rkind
        END IF

        DO j=1,CartesianTransfo%nb_RefGeometry
          DO i=1,CartesianTransfo%nat_act
            CartesianTransfo%Qxyz(:,i,j) = CartesianTransfo%Qxyz(:,i,j) * &
                                            CartesianTransfo%d0sm(3*i-2)
          END DO

          IF (CartesianTransfo%P_Axis_Ref) THEN
            CALL P_Axis_CartesianTransfo(CartesianTransfo,i_ref=j)
          END IF

        END DO

        CALL dealloc_array(d0x,'d0x','Read_CartesianTransfo')

      END IF

      !CALL Write_CartesianTransfo(CartesianTransfo)
      write(out_unitp,*) 'END ',name_sub
      END SUBROUTINE Read_CartesianTransfo

      SUBROUTINE Read_MultiRef_type0_OF_CartesianTransfo(CartesianTransfo)

      TYPE (Type_CartesianTransfo), intent(inout) :: CartesianTransfo


      character (len=Name_len) :: name_int
      integer :: i,nb_Qact
      integer, pointer :: list_Qact(:)



      integer :: err_io
      character (len=*), parameter :: name_sub='Read_MultiRef_type0_OF_CartesianTransfo'

      write(out_unitp,*) 'BEGINNING ',name_sub

      nullify(list_Qact)
      CALL alloc_array(list_Qact,(/CartesianTransfo%ncart_act/),'list_Qact',name_sub)
      list_Qact(:) = 0

      DO i=1,CartesianTransfo%ncart_act
        CALL read_name_advNo(in_unitp,name_int,err_io)

        IF (len_trim(name_int) == 0) EXIT
        !write(out_unitp,*) 'i,err_io',i,err_io
        !write(out_unitp,*) 'i,name_int',i,name_int
        read(name_int,*) list_Qact(i)
        IF (err_io /= 0) EXIT ! end of the liste

      END DO
      nb_Qact = count(list_Qact(:) > 0)
      CartesianTransfo%nb_Qact = nb_Qact
      !write(out_unitp,*) 'nb_Qact',nb_Qact
      !write(out_unitp,*) 'list_Qact',list_Qact(:)

      IF (nb_Qact == 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' nb_Qact=0 and MultiRefEckart=t'
        write(out_unitp,*) ' list_Qact',list_Qact(:)
        write(out_unitp,*) ' Check your data!'
        STOP
      END IF

      CALL alloc_CartesianTransfo(CartesianTransfo,CartesianTransfo%ncart_act, &
                                       CartesianTransfo%nb_RefGeometry,nb_Qact)
      CartesianTransfo%list_Qact(:)       = list_Qact(1:nb_Qact)
      CartesianTransfo%tab_sc(:)          = SIX
      CartesianTransfo%tab_phi0(:)        = ZERO
      CartesianTransfo%tab_switch_type(:) = 0
      CartesianTransfo%tab_nb_Ref(:)      = 2

      write(out_unitp,*) 'list_Qact',CartesianTransfo%list_Qact(:)

      IF (product(CartesianTransfo%tab_nb_Ref) /= CartesianTransfo%nb_RefGeometry) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' The product of nb_Ref should be equal to nb_RefGeometry'
        write(out_unitp,*) ' tab_nb_Ref(:)',CartesianTransfo%tab_nb_Ref(:)
        write(out_unitp,*) ' nb_RefGeometry',CartesianTransfo%nb_RefGeometry
        write(out_unitp,*) ' Check your data!'
        STOP
      END IF
      CALL dealloc_array(list_Qact,'list_Qact',name_sub)
      write(out_unitp,*) 'END ',name_sub
      END SUBROUTINE Read_MultiRef_type0_OF_CartesianTransfo

      SUBROUTINE Read_MultiRef_type1_OF_CartesianTransfo(CartesianTransfo)

      TYPE (Type_CartesianTransfo), intent(inout) :: CartesianTransfo


      integer :: i,iQact,nb_Ref,Switch_type,nb_Qact,expo
      real (kind=Rkind) :: sc,phi0

      NAMELIST /MultiRef_type1/ iQact,sc,phi0,nb_Ref,Switch_type,expo

      integer :: err_io
      character (len=*), parameter :: name_sub='Read_MultiRef_type1_OF_CartesianTransfo'

      write(out_unitp,*) 'BEGINNING ',name_sub

      nb_Qact = CartesianTransfo%nb_Qact

      IF (nb_Qact == 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' nb_Qact=0 and MultiRefEckart=t'
        write(out_unitp,*) ' Check your data!'
        STOP
      END IF

      CALL alloc_CartesianTransfo(CartesianTransfo,CartesianTransfo%ncart_act, &
                                       CartesianTransfo%nb_RefGeometry,nb_Qact)

      CartesianTransfo%list_Qact(:)       = 0
      CartesianTransfo%tab_sc(:)          = SIX
      CartesianTransfo%tab_expo(:)        = 0
      CartesianTransfo%tab_phi0(:)        = ZERO
      CartesianTransfo%tab_switch_type(:) = 0
      CartesianTransfo%tab_nb_Ref(:)      = 2

      DO i=1,nb_Qact

        iQact        = 0
        sc           = SIX
        phi0         = ZERO
        nb_Ref       = 2
        Switch_type  = 0
        expo         = 0
        read(in_unitp,MultiRef_type1,IOSTAT=err_io)
        IF (err_io < 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  while reading the namelist "MultiRef_type1" ...'
          write(out_unitp,*) ' end of file or end of record just after the namelist "Coord_transfo".'
          write(out_unitp,*) ' Probably, you have forgotten the namelist'
          write(out_unitp,*) ' Check your data !!'
          STOP
        END IF
        IF (err_io > 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  while reading the namelist "MultiRef_type1"'
          write(out_unitp,*) ' Probably, some arguments of namelist are wrong.'
          write(out_unitp,*) ' Check your data !!'
          STOP
        END IF
        write(out_unitp,MultiRef_type1)

        IF (iQact == 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' You have to specify iQact (iQact=0)'
          write(out_unitp,*) ' Check your data!'
          STOP
        END IF

        IF (Switch_type == 0) THEN ! automatic selection
          IF (expo == 0) THEN
            Switch_type = 1
          ELSE
            Switch_type = 2
          END IF
        END IF

        IF (Switch_type /= 1 .AND. Switch_type /= 2) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' Wrong value of Switch_type',Switch_type
          write(out_unitp,*) '   The possible values are 1 or 2'
          write(out_unitp,*) ' Check your data!'
          STOP
        END IF

        IF (Switch_type == 2 .AND. (expo /= 2 .AND. expo /= 4 .AND. expo /= 6)) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' The exponent has to be 2 or 4 or 6 for Switch_type=2'
          write(out_unitp,*) ' Your value, expo:',expo
          write(out_unitp,*) ' Check your data!'
          STOP
        END IF

        IF (Switch_type == 2 .AND. CartesianTransfo%tab_nb_Ref(i) < expo) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' The exponent has to be larger than or equal to the number of references'
          write(out_unitp,*) ' Your value, expo:  ',expo
          write(out_unitp,*) ' Your value, nb_ref:',CartesianTransfo%tab_nb_Ref(i)
          write(out_unitp,*) ' Check your data!'
          STOP
        END IF

        CartesianTransfo%list_Qact(i)       = iQact
        CartesianTransfo%tab_sc(i)          = sc
        CartesianTransfo%tab_expo(i)        = expo
        CartesianTransfo%tab_phi0(i)        = phi0
        CartesianTransfo%tab_switch_type(i) = switch_type
        CartesianTransfo%tab_nb_Ref(i)      = nb_Ref



      END DO

      IF (product(CartesianTransfo%tab_nb_Ref) /= CartesianTransfo%nb_RefGeometry) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' The product of nb_Ref should be equal to nb_RefGeometry'
        write(out_unitp,*) ' tab_nb_Ref(:)',CartesianTransfo%tab_nb_Ref(:)
        write(out_unitp,*) ' nb_RefGeometry',CartesianTransfo%nb_RefGeometry
        write(out_unitp,*) ' Check your data!'
        STOP
      END IF

      write(out_unitp,*) 'END ',name_sub
      END SUBROUTINE Read_MultiRef_type1_OF_CartesianTransfo

      SUBROUTINE Write_CartesianTransfo(CartesianTransfo)

      TYPE (Type_CartesianTransfo) :: CartesianTransfo


      integer :: i,j

      character (len=*), parameter :: name_sub='Write_CartesianTransfo'

      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) 'New_Orient',CartesianTransfo%New_Orient
      write(out_unitp,*) 'vAt1',CartesianTransfo%vAt1
      write(out_unitp,*) 'vAt2',CartesianTransfo%vAt2
      write(out_unitp,*) 'vAt3',CartesianTransfo%vAt3

      write(out_unitp,*) 'P_Axis_Ref',CartesianTransfo%P_Axis_Ref
      write(out_unitp,*) 'Inertia Tensor'
      CALL Write_VecMat(CartesianTransfo%InertiaT,nio=out_unitp,nbcol1=3)
      write(out_unitp,*) 'Rotational constants'
      CALL Write_VecMat(CartesianTransfo%ABC,nio=out_unitp,nbcol1=3)

      write(out_unitp,*) 'Initial Rotational matrix :'
      CALL Write_VecMat(CartesianTransfo%Rot_initial,nio=out_unitp,nbcol1=3)

      write(out_unitp,*) 'Eckart',CartesianTransfo%Eckart
      write(out_unitp,*) 'P_Axis_Always',CartesianTransfo%P_Axis_Always


      write(out_unitp,*) 'ReadRefGeometry',CartesianTransfo%ReadRefGeometry
      write(out_unitp,*) 'ncart_act',CartesianTransfo%ncart_act

      write(out_unitp,*) 'MultiRefEckart',CartesianTransfo%MultiRefEckart

      IF (CartesianTransfo%MultiRefEckart) THEN
        write(out_unitp,*) ' MultiRef_type',CartesianTransfo%MultiRef_type
        write(out_unitp,*) ' nb_RefGeometry',CartesianTransfo%nb_RefGeometry

        write(out_unitp,*) ' nb_Qact',CartesianTransfo%nb_Qact
        DO i=1,CartesianTransfo%nb_Qact
        write(out_unitp,*) ' i,Qact,sc,expo,phi0,switch_type,nb_Ref',i,' :', &
                                   CartesianTransfo%list_Qact(i),       &
                                   CartesianTransfo%tab_sc(i),          &
                                   CartesianTransfo%tab_expo(i),        &
                                   CartesianTransfo%tab_phi0(i),        &
                                   CartesianTransfo%tab_switch_type(i), &
                                   CartesianTransfo%tab_nb_Ref(i)
        END DO
      END IF

      write(out_unitp,*) 'type_diago',CartesianTransfo%type_diago
      write(out_unitp,*) 'type_cs',CartesianTransfo%type_cs
      write(out_unitp,*) 'check_dnT',CartesianTransfo%check_dnT
      write(out_unitp,*) 'dnTErr(:)',CartesianTransfo%dnTErr(:)


      IF (associated(CartesianTransfo%Qxyz)) THEN
        DO j=1,CartesianTransfo%nb_RefGeometry
        write(out_unitp,*) 'Qxyz (not mass-weigthed)',j
        DO i=1,CartesianTransfo%nat_act
          write(out_unitp,*) CartesianTransfo%Qxyz(:,i,j) /   &
                                   CartesianTransfo%d0sm(3*i-2)
        END DO
        END DO

        write(out_unitp,*) 'd0sm'
        DO i=1,size(CartesianTransfo%d0sm),3
          write(out_unitp,*) CartesianTransfo%d0sm(i:i+2)
        END DO
      END IF
      write(out_unitp,*) 'END ',name_sub

      END SUBROUTINE Write_CartesianTransfo

      SUBROUTINE P_Axis_CartesianTransfo(CartesianTransfo,i_ref)

      TYPE (Type_CartesianTransfo), intent(inout) :: CartesianTransfo
      integer, optional :: i_ref

      real (kind=Rkind) :: Ixx,Iyy,Izz,Ixy,Ixz,Iyz,Rot_PA(3,3)
      integer :: i,k,kx,ky,kz,i_ref_loc

      character (len=*), parameter :: name_sub='P_Axis_CartesianTransfo'

      i_ref_loc = 1
      IF (present(i_ref)) i_ref_loc = i_ref
      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) 'Reference geometry (in Bohr), iref:',i_ref_loc
      DO i=1,CartesianTransfo%nat_act
        write(out_unitp,*) CartesianTransfo%Qxyz(:,i,i_ref_loc) /       &
                                            CartesianTransfo%d0sm(3*i-2)
      END DO


      Ixx = ZERO
      Iyy = ZERO
      Izz = ZERO
      Ixy = ZERO
      Ixz = ZERO
      Iyz = ZERO

      DO k=1,CartesianTransfo%nat_act
         Ixx = Ixx +  CartesianTransfo%Qxyz(1,k,i_ref_loc)*CartesianTransfo%Qxyz(1,k,i_ref_loc)
         Iyy = Iyy +  CartesianTransfo%Qxyz(2,k,i_ref_loc)*CartesianTransfo%Qxyz(2,k,i_ref_loc)
         Izz = Izz +  CartesianTransfo%Qxyz(3,k,i_ref_loc)*CartesianTransfo%Qxyz(3,k,i_ref_loc)
         Ixy = Ixy +  CartesianTransfo%Qxyz(1,k,i_ref_loc)*CartesianTransfo%Qxyz(2,k,i_ref_loc)
         Ixz = Ixz +  CartesianTransfo%Qxyz(1,k,i_ref_loc)*CartesianTransfo%Qxyz(3,k,i_ref_loc)
         Iyz = Iyz +  CartesianTransfo%Qxyz(2,k,i_ref_loc)*CartesianTransfo%Qxyz(3,k,i_ref_loc)
      END DO

      CartesianTransfo%InertiaT(1,1) =  Iyy+Izz   ! Ixx
      CartesianTransfo%InertiaT(2,2) =  Ixx+Izz   ! Iyy
      CartesianTransfo%InertiaT(3,3) =  Ixx+Iyy   ! Izz

      CartesianTransfo%InertiaT(1,2) = -Ixy
      CartesianTransfo%InertiaT(1,3) = -Ixz
      CartesianTransfo%InertiaT(2,3) = -Iyz

      CartesianTransfo%InertiaT(2,1) = -Ixy
      CartesianTransfo%InertiaT(3,1) = -Ixz
      CartesianTransfo%InertiaT(3,2) = -Iyz

      write(out_unitp,*) 'Inertia Tensor, i_ref',i_ref_loc
      CALL Write_VecMat(CartesianTransfo%InertiaT,nio=out_unitp,nbcol1=3)
      CALL diagonalization(CartesianTransfo%InertiaT,                   &
                           CartesianTransfo%ABC,Rot_PA,3,1,1,.FALSE.)
      CartesianTransfo%Rot_initial(:,:) = transpose(Rot_PA)
      CartesianTransfo%ABC(:) = HALF / CartesianTransfo%ABC(:)
      CartesianTransfo%ABC(:) = CartesianTransfo%ABC(:) * 219474.63144319772_Rkind

      write(out_unitp,*) 'Rotational Matrix (for principal axis)'
      CALL Write_VecMat(CartesianTransfo%Rot_initial,nio=out_unitp,nbcol1=3)
      write(out_unitp,*) 'Rotational constants'
      CALL Write_VecMat(CartesianTransfo%ABC,nio=out_unitp,nbcol1=3)

      DO i=1,CartesianTransfo%nat_act
       CartesianTransfo%Qxyz(:,i,i_ref_loc) =                           &
                                   matmul(CartesianTransfo%Rot_initial, &
                                    CartesianTransfo%Qxyz(:,i,i_ref_loc))
      END DO

      write(out_unitp,*) 'Reference geometry in the PA frame (in Bohr), iref:',i_ref_loc
      DO i=1,CartesianTransfo%nat_act
        write(out_unitp,*) CartesianTransfo%Qxyz(:,i,i_ref_loc) /       &
                                            CartesianTransfo%d0sm(3*i-2)
      END DO

      write(out_unitp,*) 'END ',name_sub

      END SUBROUTINE P_Axis_CartesianTransfo

      SUBROUTINE CartesianTransfo1TOCartesianTransfo2(CartesianTransfo1,CartesianTransfo2)

      TYPE (Type_CartesianTransfo), intent(in)    :: CartesianTransfo1
      TYPE (Type_CartesianTransfo), intent(inout) :: CartesianTransfo2

      integer :: ncart_act,nb_RefGeometry,nb_Qact

        CartesianTransfo2%New_Orient      = CartesianTransfo1%New_Orient
        CartesianTransfo2%vAt1            = CartesianTransfo1%vAt1
        CartesianTransfo2%vAt2            = CartesianTransfo1%vAt2
        CartesianTransfo2%vAt3            = CartesianTransfo1%vAt3


        CartesianTransfo2%Eckart          = CartesianTransfo1%Eckart

        CartesianTransfo2%P_Axis_Ref      = CartesianTransfo1%P_Axis_Ref
        CartesianTransfo2%P_Axis_Always   = CartesianTransfo1%P_Axis_Always
        CartesianTransfo2%InertiaT        = CartesianTransfo1%InertiaT
        CartesianTransfo2%ABC             = CartesianTransfo1%ABC

        CartesianTransfo2%ReadRefGeometry = CartesianTransfo1%ReadRefGeometry
        CartesianTransfo2%ncart_act       = CartesianTransfo1%ncart_act
        CartesianTransfo2%nat_act         = CartesianTransfo1%nat_act
        CartesianTransfo2%Rot_initial     = CartesianTransfo1%Rot_initial

        CartesianTransfo2%nb_RefGeometry  = CartesianTransfo1%nb_RefGeometry
        CartesianTransfo2%MultiRefEckart  = CartesianTransfo1%MultiRefEckart
        CartesianTransfo2%nb_Qact         = CartesianTransfo1%nb_Qact
        CartesianTransfo2%MultiRef_type   = CartesianTransfo1%MultiRef_type


        ncart_act      = CartesianTransfo1%ncart_act
        nb_RefGeometry = CartesianTransfo1%nb_RefGeometry
        nb_Qact        = CartesianTransfo1%nb_Qact

        CALL alloc_CartesianTransfo(CartesianTransfo2,                  &
                                       ncart_act,nb_RefGeometry,nb_Qact)

        CartesianTransfo2%Qxyz       = CartesianTransfo1%Qxyz
        CartesianTransfo2%d0sm       = CartesianTransfo1%d0sm

        IF (nb_Qact > 0) THEN
          CartesianTransfo2%list_Qact       = CartesianTransfo1%list_Qact
          CartesianTransfo2%tab_sc          = CartesianTransfo1%tab_sc
          CartesianTransfo2%tab_expo        = CartesianTransfo1%tab_expo
          CartesianTransfo2%tab_phi0        = CartesianTransfo1%tab_phi0
          CartesianTransfo2%tab_switch_type = CartesianTransfo1%tab_switch_type
          CartesianTransfo2%tab_nb_Ref      = CartesianTransfo1%tab_nb_Ref
        END IF

        CartesianTransfo2%type_diago      = CartesianTransfo1%type_diago
        CartesianTransfo2%type_cs         = CartesianTransfo1%type_cs
        CartesianTransfo2%check_dnT       = CartesianTransfo1%check_dnT
        CartesianTransfo2%dnTErr(:)       = CartesianTransfo1%dnTErr(:)

      END SUBROUTINE CartesianTransfo1TOCartesianTransfo2


      SUBROUTINE calc_CartesianTransfo_new(dnQin,dnQout,CartesianTransfo,&
                                           Qact,nderiv,inTOout)

        TYPE (Type_dnVec), intent(inout)         :: dnQin,dnQout
        TYPE (Type_CartesianTransfo), intent(in) :: CartesianTransfo
        real (kind=Rkind), intent(in)            :: Qact(:)
        integer, intent(in)                      :: nderiv
        logical, intent(in)                      :: inTOout

        ! from Dymarsky and Kudin ref JCP v122, p124103, 2005
        TYPE(Type_dnS) :: dnT(3,3)
        TYPE(Type_dnS) :: dnXref(3,CartesianTransfo%nat_act)

        character (len=Name_longlen) :: RMatIO_format_save


!----- for debuging --------------------------------------------------
        integer :: nderiv_debug = 1
        character (len=*), parameter :: name_sub='calc_CartesianTransfo_new'
        logical, parameter :: debug=.FALSE.
        !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
          write(out_unitp,*) 'ncart_act ',CartesianTransfo%ncart_act
          write(out_unitp,*) 'nat_act ',CartesianTransfo%nat_act
          write(out_unitp,*) 'Qxyz (ref) ',CartesianTransfo%Qxyz
          write(out_unitp,*) 'dnQin%d0 '
          CALL write_dnx(1,CartesianTransfo%ncart_act,dnQin,nderiv_debug)
          write(out_unitp,*) 'Qact ',Qact
          CALL Write_CartesianTransfo(CartesianTransfo)
          CALL flush_perso(out_unitp)
        END IF
!-----------------------------------------------------------
        RMatIO_format_save = RMatIO_format
        RMatIO_format = "f30.10"
        CALL check_alloc_dnVec(dnQin,'dnQin',name_sub)
        CALL check_alloc_dnVec(dnQout,'dnQout',name_sub)


        IF (inTOout) THEN ! for Q TO x

          !==============================================================
          ! rotation of the CC with an initial and constant rotational matrix
          CALL alloc_MatOFdnS(dnT,dnQin%nb_var_deriv,nderiv)
          CALL sub_ZERO_TO_MatOFdnS(dnT)
          dnT(:,:)%d0 = CartesianTransfo%Rot_initial

          CALL calc_dnTxdnXin_TO_dnXout(dnQin,dnT,dnQout,CartesianTransfo,nderiv)
          CALL sub_dnVec1_TO_dnVec2(dnQout,dnQin)
          !==============================================================

          !==============================================================
          ! from Dymarsky and Kudin, JCP v122, p124103, 2005
          IF (CartesianTransfo%Eckart) THEN

            CALL alloc_MatOFdnS(dnXref,dnQin%nb_var_deriv,nderiv)

            CALL dnX_MultiRef(dnXref,CartesianTransfo,Qact(1:CartesianTransfo%nb_Qact),dnQin)
            CALL calc_dnTEckart(dnQin,dnT,dnXref,CartesianTransfo,nderiv)
            CALL calc_dnTxdnXin_TO_dnXout(dnQin,dnT,dnQout,CartesianTransfo,nderiv)

            IF (debug)                                                  &
              CALL calc_Analysis_dnXout(dnQout,dnT,dnXref,CartesianTransfo,Qact,nderiv)

            CALL dealloc_MatOFdnS(dnXref)

          END IF

          CALL dealloc_MatOFdnS(dnT)


        ELSE
          CALL sub_dnVec1_TO_dnVec2(dnQout,dnQin)
        END IF

        IF (debug) THEN
          write(out_unitp,*) 'dnQout%d0 '
          CALL write_dnx(1,CartesianTransfo%ncart_act,dnQout,nderiv_debug)

          write(out_unitp,*) 'END ',name_sub
          CALL flush_perso(out_unitp)
        END IF
        RMatIO_format = RMatIO_format_save

      END SUBROUTINE calc_CartesianTransfo_new

      SUBROUTINE calc_CartesianTransfo_old(dnQin,dnQout,CartesianTransfo,&
                                           Qact,nderiv,inTOout)

        TYPE (Type_dnVec), intent(inout)         :: dnQin,dnQout
        TYPE (Type_CartesianTransfo), intent(in) :: CartesianTransfo
        real (kind=Rkind), intent(in)            :: Qact(:)


        integer, intent(in)               :: nderiv
        logical, intent(in)               :: inTOout

        real (kind=Rkind) :: Rot_Eckart(3)
        integer           :: i,j,ic,ix,iy,iz,ixyz,jxyz,iref,nb_ref

        ! from Dymarsky and Kudin ref JCP v122, p124103, 2005
        TYPE(Type_dnS) :: dnA(3,3),dntA(3,3),dnA1(3,3),dnA2(3,3),dnT(3,3)
        TYPE(Type_dnS) :: dnEig1(3),dnEig2(3),dnVec1(3,3),dnVec2(3,3)
        TYPE(Type_dnS) :: dnWork,dnWork2,dnVec1Work(3),dnVec2Work(3)


        TYPE(Type_dnS) :: dnXref(3,CartesianTransfo%ncart_act/3)


        real (kind=Rkind) :: norm,dp(3),dp_save(3)
        integer           :: isort_dp(3),idp(1)

        real(kind=Rkind) :: max_val
        integer          :: max_j

        character (len=Name_longlen) :: RMatIO_format_save


!----- for debuging --------------------------------------------------
        integer :: nderiv_debug = 1
        character (len=*), parameter :: name_sub='calc_CartesianTransfo_old'
        logical, parameter :: debug=.FALSE.
        !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
          write(out_unitp,*) 'ncart_act ',CartesianTransfo%ncart_act
          write(out_unitp,*) 'nat_act ',CartesianTransfo%nat_act
          write(out_unitp,*) 'Qxyz (ref) ',CartesianTransfo%Qxyz
          CALL write_dnx(1,CartesianTransfo%ncart_act,dnQin,nderiv_debug)
          !write(out_unitp,*) 'dnQin%d0 ',dnQin%d0(1:CartesianTransfo%ncart_act)
          write(out_unitp,*) 'Qact ',Qact
          CALL Write_CartesianTransfo(CartesianTransfo)
          CALL flush_perso(out_unitp)
        END IF
!-----------------------------------------------------------
        RMatIO_format_save = RMatIO_format
        RMatIO_format = "f30.10"
        CALL check_alloc_dnVec(dnQin,'dnQin',name_sub)
        CALL check_alloc_dnVec(dnQout,'dnQout',name_sub)


        IF (inTOout) THEN ! for Q TO x

          !==============================================================
          ! allocation
          CALL alloc_MatOFdnS(dnA,dnQin%nb_var_deriv,nderiv)
          CALL alloc_MatOFdnS(dntA,dnQin%nb_var_deriv,nderiv)
          CALL alloc_MatOFdnS(dnA1,dnQin%nb_var_deriv,nderiv)
          CALL alloc_MatOFdnS(dnA2,dnQin%nb_var_deriv,nderiv)
          CALL alloc_MatOFdnS(dnT,dnQin%nb_var_deriv,nderiv)
          CALL alloc_MatOFdnS(dnVec1,dnQin%nb_var_deriv,nderiv)
          CALL alloc_MatOFdnS(dnVec2,dnQin%nb_var_deriv,nderiv)

          CALL alloc_VecOFdnS(dnEig1,dnQin%nb_var_deriv,nderiv)
          CALL alloc_VecOFdnS(dnEig2,dnQin%nb_var_deriv,nderiv)
          CALL alloc_VecOFdnS(dnVec1Work,dnQin%nb_var_deriv,nderiv)
          CALL alloc_VecOFdnS(dnVec2Work,dnQin%nb_var_deriv,nderiv)

          CALL alloc_dnS(dnWork,dnQin%nb_var_deriv,nderiv)
          CALL alloc_dnS(dnWork2,dnQin%nb_var_deriv,nderiv)
          !==============================================================

          !==============================================================
          ! rotation of the CC with an initial and constant rotational matrix
          dnT(:,:)%d0 = CartesianTransfo%Rot_initial
          DO ic=1,CartesianTransfo%ncart_act,3
            CALL sub_dnVec_TO_dnS(dnQin,dnVec1Work(1),ic+0,nderiv)
            CALL sub_dnVec_TO_dnS(dnQin,dnVec1Work(2),ic+1,nderiv)
            CALL sub_dnVec_TO_dnS(dnQin,dnVec1Work(3),ic+2,nderiv)

            CALL Mat1OFdnS_MUL_Vec2OFdnS_TO_Vec3OFdnS(dnT,dnVec1Work,dnVec2Work,nderiv)

            CALL sub_dnS_TO_dnVec(dnVec2Work(1),dnQin,ic+0,nderiv)
            CALL sub_dnS_TO_dnVec(dnVec2Work(2),dnQin,ic+1,nderiv)
            CALL sub_dnS_TO_dnVec(dnVec2Work(3),dnQin,ic+2,nderiv)
          END DO
          !==============================================================

          !==============================================================
          ! from Dymarsky and Kudin, JCP v122, p124103, 2005
          !write(out_unitp,*) 'Eckart',CartesianTransfo%Eckart

          IF (CartesianTransfo%Eckart) THEN

            CALL sub_ZERO_TO_MatOFdnS(dnA)

            CALL alloc_MatOFdnS(dnXref,dnQin%nb_var_deriv,nderiv)

            CALL dnX_MultiRef(dnXref,CartesianTransfo,Qact,dnQin)

            CALL sub_ZERO_TO_MatOFdnS(dnA,nderiv)
            DO ic=1,CartesianTransfo%ncart_act,3 ! loop on atoms
            DO i=1,3
              CALL sub_dnVec_TO_dnS(dnQin,dnWork,ic-1+i,nderiv)

              DO j=1,3
                CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnWork,dnXref(j,(ic+2)/3),dnWork2,nderiv)
                CALL sub_dnS1_PLUS_dnS2_TO_dnS2(dnWork2,dnA(i,j),nderiv)
              END DO
            END DO
            END DO
            CALL dealloc_dnS(dnWork2)
            IF (debug) THEN
              write(out_unitp,*) 'dnA'
              CALL Write_MatOFdnS(dnA,nderiv=0)
            END IF

            CALL TRANS_Mat1OFdnS_TO_Mat2OFdnS(dnA,dntA,nderiv)

            CALL Mat1OFdnS_MUL_Mat2OFdnS_TO_Mat3OFdnS(dntA,dnA,dnA1,nderiv)
            CALL Mat1OFdnS_MUL_Mat2OFdnS_TO_Mat3OFdnS(dnA,dntA,dnA2,nderiv)

            IF (debug) THEN
              write(out_unitp,*) 'dnA1'
              CALL Write_MatOFdnS(dnA1,nderiv=0)

              write(out_unitp,*) 'dnA2'
              CALL Write_MatOFdnS(dnA2,nderiv=0)

            END IF

            CALL DIAG_MatOFdnS(dnA1,dnVec1,nderiv,sort=1,type_diago=4)
            CALL DIAG_MatOFdnS(dnA2,dnVec2,nderiv,sort=1,type_diago=4)

            IF (debug) THEN
              write(out_unitp,*) 'dnVec1 before change sign (in column)'
              CALL Write_MatOFdnS(dnVec1,nderiv=0)
              write(out_unitp,*) 'dnVec2 before change sign (in column)'
              CALL Write_MatOFdnS(dnVec2,nderiv=0)
            END IF


            ! change the sign of Vec2(:,i)
            DO i=1,3
              dp(i) = dot_product(dnVec1(:,i)%d0,dnVec2(:,i)%d0)
            END DO
            dp_save(:) = dp(:)
            !write(out_unitp,*) 'dot_pro',dp(:)
            idp(:) = maxloc(abs(dp))
            isort_dp(1) = idp(1)
            dp(isort_dp(1)) = ZERO
            idp(:) = maxloc(abs(dp))
            isort_dp(2) = idp(1)
            dp(isort_dp(2)) = ZERO
            idp(:) = maxloc(abs(dp))
            isort_dp(3) = idp(1)
            !write(out_unitp,*) 'sort order of dot_pro',isort_dp(:)

            DO i=1,2
              IF (dp_save(isort_dp(i)) < ZERO) THEN ! change sign of dnVec2
                !write(6,*) 'change sign of vec: ',isort_dp(i)
                DO j=1,3
                  CALL sub_dnS1_PROD_w_TO_dnS2(dnVec2(j,isort_dp(i)),-ONE,&
                                               dnVec2(j,isort_dp(i)),nderiv)
                END DO
              END IF
            END DO

            ! Cross product for the third vector
            CALL Vec1OFdnS_CROSSPRODUCT_Vec2OFdnS_TO_Vec3OFdnS(         &
                                           dnVec1(:,isort_dp(1)),       &
                                           dnVec1(:,isort_dp(2)),       &
                                           dnVec1(:,isort_dp(3)),nderiv)
            CALL Vec1OFdnS_CROSSPRODUCT_Vec2OFdnS_TO_Vec3OFdnS(         &
                                           dnVec2(:,isort_dp(1)),       &
                                           dnVec2(:,isort_dp(2)),       &
                                           dnVec2(:,isort_dp(3)),nderiv)

            IF (debug) THEN
              write(out_unitp,*) 'dnVec1 (in column)'
              CALL Write_MatOFdnS(dnVec1,nderiv=0)
              write(out_unitp,*) 'dnVec2 (in column)'
              CALL Write_MatOFdnS(dnVec2,nderiv=0)
            END IF

            DO i=1,3
            DO j=1,3
              CALL Vec1OFdnS_DOTPRODUCT_Vec2OFdnS_TO_dnS3(dnVec1(i,:),  &
                                            dnVec2(j,:),dnT(i,j),nderiv)
            END DO
            END DO

            IF (debug) THEN
              write(out_unitp,*) 'eig1 ',(dnA1(i,i)%d0,i=1,3)
              write(out_unitp,*) 'eig2 ',(dnA2(i,i)%d0,i=1,3)
              CALL flush_perso(out_unitp)

              ! check the determniant of dnT (and its derivatives)
              !first line
              CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnT(2,2),dnT(3,3),        &
                                                          dnVec1Work(1))
              CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnT(2,3),dnT(3,2),        &
                                                          dnVec1Work(2))
              CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnVec1Work(1),ONE,       &
                                       dnVec1Work(2),-ONE,dnVec1Work(3))
              CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnT(1,1),dnVec1Work(3),   &
                                                                 dnWork)

              !second line
              CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnT(1,2),dnT(3,3),        &
                                                          dnVec1Work(1))
              CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnT(1,3),dnT(3,2),        &
                                                          dnVec1Work(2))
              CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnVec1Work(1),ONE,       &
                                       dnVec1Work(2),-ONE,dnVec1Work(3))
              CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnT(2,1),dnVec1Work(3),   &
                                                                 dnWork2)
              CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnWork,ONE,              &
                                                    dnWork2,-ONE,dnWork)

              !second line
              CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnT(1,2),dnT(2,3),        &
                                                          dnVec1Work(1))
              CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnT(1,3),dnT(2,2),        &
                                                          dnVec1Work(2))
              CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnVec1Work(1),ONE,       &
                                       dnVec1Work(2),-ONE,dnVec1Work(3))
              CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnT(3,1),dnVec1Work(3),   &
                                                                 dnWork2)
              CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnWork,ONE,              &
                                                    dnWork2, ONE,dnWork)

              write(out_unitp,*) 'det(dnT)'
              CALL Write_dnS(dnWork)

              write(out_unitp,*) 'Eckart rotational matrix, T + det(T)',dnWork%d0
              CALL Write_MatOFdnS(dnT,nderiv=0)
              IF (dnWork%d0 < ZERO) STOP 'det(T)=-1'

            END IF

            ! rotation of the CC
            DO ic=1,CartesianTransfo%ncart_act,3
              CALL sub_dnVec_TO_dnS(dnQin,dnVec1Work(1),ic+0,nderiv)
              CALL sub_dnVec_TO_dnS(dnQin,dnVec1Work(2),ic+1,nderiv)
              CALL sub_dnVec_TO_dnS(dnQin,dnVec1Work(3),ic+2,nderiv)

              CALL Mat1OFdnS_MUL_Vec2OFdnS_TO_Vec3OFdnS(dnT,dnVec1Work,dnVec2Work,nderiv)

              CALL sub_dnS_TO_dnVec(dnVec2Work(1),dnQin,ic+0,nderiv)
              CALL sub_dnS_TO_dnVec(dnVec2Work(2),dnQin,ic+1,nderiv)
              CALL sub_dnS_TO_dnVec(dnVec2Work(3),dnQin,ic+2,nderiv)
            END DO

            ! Check the rotational Eckart condition
            IF (debug) THEN

              DO iref=1,CartesianTransfo%nb_RefGeometry
                Rot_Eckart(:) = ZERO
                DO ic=1,CartesianTransfo%ncart_act,3

                  ix = ic+0
                  iy = ic+1
                  iz = ic+2

                  Rot_Eckart(1) = Rot_Eckart(1)  +                      &
                  dnQin%d0(iy)*CartesianTransfo%Qxyz(3,(ic+2)/3,iref) - &
                  dnQin%d0(iz)*CartesianTransfo%Qxyz(2,(ic+2)/3,iref)

                  Rot_Eckart(2) = Rot_Eckart(2) -                       &
                  dnQin%d0(ix)*CartesianTransfo%Qxyz(3,(ic+2)/3,iref) + &
                  dnQin%d0(iz)*CartesianTransfo%Qxyz(1,(ic+2)/3,iref)

                  Rot_Eckart(3) = Rot_Eckart(3) +                       &
                  dnQin%d0(ix)*CartesianTransfo%Qxyz(2,(ic+2)/3,iref) - &
                  dnQin%d0(iy)*CartesianTransfo%Qxyz(1,(ic+2)/3,iref)

                END DO
                norm = sqrt(dot_product(Rot_Eckart(:),Rot_Eckart(:)))
                IF (CartesianTransfo%MultiRefEckart) THEN
                   write(out_unitp,*) 'Norm Rot_Eckart',iref,            &
                               Qact(CartesianTransfo%list_Qact(1)),norm
                ELSE
                   write(out_unitp,*) 'Norm Rot_Eckart',norm
                END IF
              END DO
            END IF

            IF (debug) THEN
              Rot_Eckart(:) = ZERO
              DO ic=1,CartesianTransfo%ncart_act,3
                ix = ic+0
                iy = ic+1
                iz = ic+2

                Rot_Eckart(1) = Rot_Eckart(1) +                         &
                                   dnQin%d0(iy)*dnXref(3,(ic+2)/3)%d0 - &
                                   dnQin%d0(iz)*dnXref(2,(ic+2)/3)%d0

                Rot_Eckart(2) = Rot_Eckart(2) -                         &
                                   dnQin%d0(ix)*dnXref(3,(ic+2)/3)%d0 + &
                                   dnQin%d0(iz)*dnXref(1,(ic+2)/3)%d0

                Rot_Eckart(3) = Rot_Eckart(3) +                         &
                                   dnQin%d0(ix)*dnXref(2,(ic+2)/3)%d0 - &
                                   dnQin%d0(iy)*dnXref(1,(ic+2)/3)%d0

              END DO

              norm = sqrt(dot_product(Rot_Eckart(:),Rot_Eckart(:)))
              write(out_unitp,*) 'Norm Rot_Eckart (av ref)',norm
            END IF

            IF (debug) THEN
              write(out_unitp,*) 'Actual geometry (mass weighted): '
              DO ic=1,CartesianTransfo%ncart_act,3
                write(out_unitp,*) (ic+2)/3,dnQin%d0(ic+0:ic+2)
              END DO
              write(out_unitp,*) 'Reference geometry (mass weighted): '
              DO ic=1,CartesianTransfo%nat_act
                write(out_unitp,*) ic,dnXref(:,ic)%d0
              END DO

              write(out_unitp,*) 'Actual geometry (not mass weighted): '
              DO ic=1,CartesianTransfo%ncart_act,3
                write(out_unitp,*) (ic+2)/3,dnQin%d0(ic+0:ic+2) / CartesianTransfo%d0sm(ic)
              END DO
              write(out_unitp,*) 'Reference geometry (not mass weighted): '
              DO ic=1,CartesianTransfo%nat_act
                write(out_unitp,*) ic,dnXref(:,ic)%d0 / CartesianTransfo%d0sm(3*ic-2)
              END DO

            END IF

            CALL dealloc_MatOFdnS(dnXref)


            ! Check the "RMS" between the reference and the actual geometry
            norm = ZERO
            i = 0
            DO ic=1,CartesianTransfo%nat_act
            DO j=1,3
              i = i + 1
              norm = norm + ( dnQin%d0(i) - CartesianTransfo%Qxyz(j,ic,1) )**2
            END DO
            END DO

            norm = sqrt(norm/real(CartesianTransfo%ncart_act,kind=Rkind))
            IF (debug .AND. norm > TEN) write(out_unitp,*) 'large RMS',norm
!            write(out_unitp,*) 'RMS',norm
          END IF ! End Eckart
          !==============================================================


          !==============================================================
          ! deallocation
          CALL dealloc_MatOFdnS(dnA)
          CALL dealloc_MatOFdnS(dntA)
          CALL dealloc_MatOFdnS(dnA1)
          CALL dealloc_MatOFdnS(dnA2)
          CALL dealloc_MatOFdnS(dnT)
          CALL dealloc_MatOFdnS(dnVec1)
          CALL dealloc_MatOFdnS(dnVec2)

          CALL dealloc_VecOFdnS(dnEig1)
          CALL dealloc_VecOFdnS(dnEig2)
          CALL dealloc_VecOFdnS(dnVec1Work)
          CALL dealloc_VecOFdnS(dnVec2Work)

          CALL dealloc_dnS(dnWork)
          CALL dealloc_dnS(dnWork2)
          !==============================================================

          CALL sub_dnVec1_TO_dnVec2(dnQin,dnQout)
        ELSE
          CALL sub_dnVec1_TO_dnVec2(dnQout,dnQin)
        END IF

        IF (debug) THEN
          !write(out_unitp,*) 'dnQout%d0 ',dnQout%d0(1:CartesianTransfo%ncart_act)
          CALL write_dnx(1,CartesianTransfo%ncart_act,dnQout,nderiv_debug)

          write(out_unitp,*) 'END ',name_sub
          CALL flush_perso(out_unitp)
        END IF
        RMatIO_format = RMatIO_format_save

      END SUBROUTINE calc_CartesianTransfo_old

      SUBROUTINE calc_dnTxdnXin_TO_dnXout(dnXin,dnT,dnXout,CartesianTransfo,nderiv)

        TYPE (Type_dnVec), intent(in)            :: dnXin
        TYPE(Type_dnS), intent(in)               :: dnT(3,3)
        TYPE (Type_dnVec), intent(inout)         :: dnXout
        TYPE (Type_CartesianTransfo), intent(in) :: CartesianTransfo

        integer, intent(in)                      :: nderiv

        integer           :: ic

        TYPE(Type_dnS) :: dnVec1Work(3),dnVec2Work(3)

!----- for debuging --------------------------------------------------
        integer :: nderiv_debug = 1
        character (len=*), parameter :: name_sub='calc_dnTxdnXin_TO_dnXout'
        logical, parameter :: debug=.FALSE.
        !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
          CALL flush_perso(out_unitp)
        END IF
!-----------------------------------------------------------
        CALL check_alloc_dnVec(dnXin,'dnXin',name_sub)
        CALL check_alloc_MatOFdnS(dnT,'dnT',name_sub)


        !==============================================================
        ! allocation
        CALL alloc_VecOFdnS(dnVec1Work,dnXin%nb_var_deriv,nderiv)
        CALL alloc_VecOFdnS(dnVec2Work,dnXin%nb_var_deriv,nderiv)
        !==============================================================

        !==============================================================
        ! rotation of the CC with an initial and constant rotational matrix
        DO ic=1,CartesianTransfo%ncart_act,3
          CALL sub_dnVec_TO_dnS(dnXin,dnVec1Work(1),ic+0,nderiv)
          CALL sub_dnVec_TO_dnS(dnXin,dnVec1Work(2),ic+1,nderiv)
          CALL sub_dnVec_TO_dnS(dnXin,dnVec1Work(3),ic+2,nderiv)

          CALL Mat1OFdnS_MUL_Vec2OFdnS_TO_Vec3OFdnS(dnT,dnVec1Work,dnVec2Work,nderiv)

          CALL sub_dnS_TO_dnVec(dnVec2Work(1),dnXout,ic+0,nderiv)
          CALL sub_dnS_TO_dnVec(dnVec2Work(2),dnXout,ic+1,nderiv)
          CALL sub_dnS_TO_dnVec(dnVec2Work(3),dnXout,ic+2,nderiv)
        END DO
        !==============================================================

        !==============================================================
        ! deallocation
        CALL dealloc_VecOFdnS(dnVec1Work)
        CALL dealloc_VecOFdnS(dnVec2Work)
        !==============================================================


        IF (debug) THEN
          CALL write_dnx(1,CartesianTransfo%ncart_act,dnXout,nderiv_debug)

          write(out_unitp,*) 'END ',name_sub
          CALL flush_perso(out_unitp)
        END IF

      END SUBROUTINE calc_dnTxdnXin_TO_dnXout

      SUBROUTINE calc_RMS_d0Xout(sign_Vec,Xin,Vec1,Vec2,Xref,nat)

      integer, intent(in)              :: nat
      real (kind=Rkind), intent(in)    :: Xin(3,nat),Xref(3,nat)
      real (kind=Rkind), intent(in)    :: Vec1(3,3),Vec2(3,3)
      real (kind=Rkind), intent(inout) :: sign_Vec(3)


      real (kind=Rkind) :: RMS(4),RMS_min,det,normEC,Xout(3,nat)
      real (kind=Rkind) :: T(3,3),Rot_Eckart(3)
      integer :: iat,i,j,k
      integer :: irot,irot_min
      real (kind=Rkind), parameter :: sig(3,4) = reshape( (/            &
        ONE, ONE, ONE, -ONE, ONE,-ONE, ONE,-ONE,-ONE, -ONE,-ONE, ONE /),&
                                                              (/ 3,4 /))

!----- for debuging --------------------------------------------------
        character (len=*), parameter :: name_sub='calc_RMS_d0Xout'
        logical, parameter :: debug=.FALSE.
        !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------

      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL flush_perso(out_unitp)
      END IF

      RMS_min  = huge(ONE)
      irot_min = -1

      DO irot=1,4
        DO i=1,3
        DO j=1,3
          T(i,j) = ZERO
          DO k=1,3
            T(i,j) = T(i,j) + Vec1(i,k)*Vec2(j,k)*sig(k,irot)
          END DO
        END DO
        END DO

        DO iat=1,nat
          Xout(:,iat) = matmul(T,Xin(:,iat))
        END DO

        IF (debug) THEN
          write(out_unitp,*) 'T'
          CALL Write_Mat(T,out_unitp,3)
          CALL Det_OF_m1(T,det,3)
          write(out_unitp,*) 'det T',irot,det
          write(out_unitp,*) 'Xin',Xin
          write(out_unitp,*) 'Xref',Xref
          write(out_unitp,*) 'Xout',Xout
        END IF


        ! Check the "RMS" between the reference and the actual geometry
        RMS(irot) = sum((Xout-Xref)**2)
        RMS(irot) = sqrt(RMS(irot)/real(3*nat,kind=Rkind))
        IF (RMS(irot) < RMS_min) THEN
          irot_min = irot
          RMS_min  = RMS(irot)
          sign_Vec(:) = sig(:,irot)
        END IF

        IF (debug) write(out_unitp,*) 'RMS',irot,RMS(irot)

        Rot_Eckart(:) = ZERO
        DO iat=1,nat

          Rot_Eckart(1) = Rot_Eckart(1)  +                              &
                       Xout(2,iat)*Xref(3,iat) - Xout(3,iat)*Xref(2,iat)

          Rot_Eckart(2) = Rot_Eckart(2) -                               &
                       Xout(1,iat)*Xref(3,iat) + Xout(3,iat)*Xref(1,iat)

          Rot_Eckart(3) = Rot_Eckart(3) +                               &
                       Xout(1,iat)*Xref(2,iat) - Xout(2,iat)*Xref(1,iat)

        END DO
        normEC = sqrt(dot_product(Rot_Eckart(:),Rot_Eckart(:)))
        IF (debug) write(out_unitp,*) 'normEC',irot,normEC


      END DO

      !write(out_unitp,*) 'RMS_min',irot_min,RMS_min,sign_Vec(:)

      IF (debug) THEN
        write(out_unitp,*) 'RMS_min',irot_min,RMS_min,sign_Vec(:)
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF

      END SUBROUTINE calc_RMS_d0Xout

      SUBROUTINE calc_Analysis_dnXout(dnXout,dnT,dnXref,CartesianTransfo,Qact,nderiv)

        TYPE (Type_dnVec), intent(in)            :: dnXout
        TYPE(Type_dnS), intent(in)               :: dnT(3,3)
        TYPE(Type_dnS), intent(in)               :: dnXref(:,:)
        TYPE (Type_CartesianTransfo), intent(in) :: CartesianTransfo
        real (kind=Rkind), intent(in)            :: Qact(:)

        integer, intent(in)                      :: nderiv

        integer           :: ic,iref,ix,iy,iz
        real (kind=Rkind) :: norm,Rot_Eckart(3),RMS2

!----- for debuging --------------------------------------------------
        integer :: nderiv_debug = 1
        character (len=*), parameter :: name_sub='calc_Analysis_dnXout'
        logical, parameter :: debug=.FALSE.
        !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
!$OMP  CRITICAL (calc_Analysis_dnXout_CRIT)

        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL flush_perso(out_unitp)
!-----------------------------------------------------------
        CALL check_alloc_dnVec(dnXout,'dnXout',name_sub)
        CALL check_alloc_MatOFdnS(dnT,'dnT',name_sub)


        !Check the rotational Eckart condition
        RMS2 = sum((dnXout%d0(:)-reshape(dnXref(:,:)%d0, (/ CartesianTransfo%ncart_act /)) )**2)
        RMS2 = RMS2/real(3*CartesianTransfo%nat_act,kind=Rkind)
        write(out_unitp,*) 'RMS^2',Qact,RMS2

        !Check the rotational Eckart condition
        DO iref=1,CartesianTransfo%nb_RefGeometry
          Rot_Eckart(:) = ZERO
          DO ic=1,CartesianTransfo%ncart_act,3

            ix = ic+0
            iy = ic+1
            iz = ic+2

            Rot_Eckart(1) = Rot_Eckart(1)  +                            &
                dnXout%d0(iy)*CartesianTransfo%Qxyz(3,(ic+2)/3,iref) -  &
                dnXout%d0(iz)*CartesianTransfo%Qxyz(2,(ic+2)/3,iref)

            Rot_Eckart(2) = Rot_Eckart(2) -                             &
                 dnXout%d0(ix)*CartesianTransfo%Qxyz(3,(ic+2)/3,iref) + &
                 dnXout%d0(iz)*CartesianTransfo%Qxyz(1,(ic+2)/3,iref)

            Rot_Eckart(3) = Rot_Eckart(3) +                             &
                 dnXout%d0(ix)*CartesianTransfo%Qxyz(2,(ic+2)/3,iref) - &
                 dnXout%d0(iy)*CartesianTransfo%Qxyz(1,(ic+2)/3,iref)

          END DO
          norm = sqrt(dot_product(Rot_Eckart(:),Rot_Eckart(:)))
          IF (CartesianTransfo%MultiRefEckart) THEN
             write(out_unitp,*) 'Norm Rot_Eckart',iref,Qact(:),norm
          ELSE
             write(out_unitp,*) 'Norm Rot_Eckart',norm
          END IF
        END DO

        IF (CartesianTransfo%nb_RefGeometry > 1) THEN
          Rot_Eckart(:) = ZERO
          DO ic=1,CartesianTransfo%ncart_act,3
            ix = ic+0
            iy = ic+1
            iz = ic+2

            Rot_Eckart(1) = Rot_Eckart(1) +                             &
                                  dnXout%d0(iy)*dnXref(3,(ic+2)/3)%d0 - &
                                  dnXout%d0(iz)*dnXref(2,(ic+2)/3)%d0

            Rot_Eckart(2) = Rot_Eckart(2) -                             &
                                  dnXout%d0(ix)*dnXref(3,(ic+2)/3)%d0 + &
                                  dnXout%d0(iz)*dnXref(1,(ic+2)/3)%d0

            Rot_Eckart(3) = Rot_Eckart(3) +                             &
                                  dnXout%d0(ix)*dnXref(2,(ic+2)/3)%d0 - &
                                  dnXout%d0(iy)*dnXref(1,(ic+2)/3)%d0

          END DO

          norm = sqrt(dot_product(Rot_Eckart(:),Rot_Eckart(:)))
          write(out_unitp,*) 'Norm Rot_Eckart (av ref)',norm
        END IF

        write(out_unitp,*) 'Actual geometry (mass weighted): '
        DO ic=1,CartesianTransfo%ncart_act,3
          write(out_unitp,*) (ic+2)/3,dnXout%d0(ic+0:ic+2)
        END DO
        write(out_unitp,*) 'Reference geometry (mass weighted): '
        DO ic=1,CartesianTransfo%nat_act
          write(out_unitp,*) ic,dnXref(:,ic)%d0
        END DO

        write(out_unitp,*) 'Actual geometry (not mass weighted): '
        DO ic=1,CartesianTransfo%ncart_act,3
          write(out_unitp,*) (ic+2)/3,dnXout%d0(ic+0:ic+2) / CartesianTransfo%d0sm(ic)
        END DO
        write(out_unitp,*) 'Reference geometry (not mass weighted): '
        DO ic=1,CartesianTransfo%nat_act
          write(out_unitp,*) ic,dnXref(:,ic)%d0 / CartesianTransfo%d0sm(3*ic-2)
        END DO

        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
!$OMP  END CRITICAL (calc_Analysis_dnXout_CRIT)


      END SUBROUTINE calc_Analysis_dnXout

      SUBROUTINE calc_EckartRot(dnx,T,CartesianTransfo,Qact)

        TYPE (Type_dnVec), intent(inout)  :: dnx
        TYPE (Type_CartesianTransfo), intent(in) :: CartesianTransfo
        real (kind=Rkind), intent(in) :: Qact(:)

        real (kind=Rkind) :: T(3,3),d0x(3,CartesianTransfo%ncart_act/3)
        real (kind=Rkind) :: Rot_Eckart(3)

        integer           :: i,j,ic,ix,iy,iz,iref

        ! from Dymarsky and Kudin ref JCP v122, p124103, 2005
        real (kind=Rkind) :: A(3,3),A1(3,3),A2(3,3)
        real (kind=Rkind) :: vec1(3,3),vec2(3,3),eig1(3),eig2(3)
        real (kind=Rkind) :: normx,normy,normz,norm

        TYPE(Type_dnS) :: dnXref(3,CartesianTransfo%ncart_act/3)


!------ for debuging --------------------------------------------------
        character (len=*), parameter :: name_sub='calc_EckartRot'
        logical, parameter :: debug=.FALSE.
        !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
          write(out_unitp,*) 'ncart_act ',CartesianTransfo%ncart_act
          write(out_unitp,*) 'Qxyz (ref) ',CartesianTransfo%Qxyz
          write(out_unitp,*) 'dnx%d0 ',dnx%d0(1:CartesianTransfo%ncart_act)
        END IF
!-----------------------------------------------------------

        CALL check_alloc_dnVec(dnx,'dnx',name_sub)

        !==============================================================
        ! rotation of the CC with an initial and constant matrix
        DO ic=1,CartesianTransfo%ncart_act,3
          ix = ic+0
          iz = ic+2
          dnx%d0(ix:iz) = dnx%d0(ix:iz)/CartesianTransfo%d0sm(ix)
          dnx%d0(ix:iz) = matmul(CartesianTransfo%Rot_initial,dnx%d0(ix:iz))
          dnx%d0(ix:iz) = dnx%d0(ix:iz)*CartesianTransfo%d0sm(ix)
        END DO

        CALL alloc_MatOFdnS(dnXref,dnx%nb_var_deriv,0)

        CALL dnX_MultiRef(dnXref,CartesianTransfo,Qact,dnx)

        ! from Dymarsky and Kudin ref JCP v122, p124103, 2005
        A(:,:) = ZERO
        DO ic=1,CartesianTransfo%ncart_act,3 ! loop on atoms
          ix = ic+0
          iz = ic+2
          A(:,1) = A(:,1) + dnx%d0(ix:iz)*dnXref(1,(ic+2)/3)%d0
          A(:,2) = A(:,2) + dnx%d0(ix:iz)*dnXref(2,(ic+2)/3)%d0
          A(:,3) = A(:,3) + dnx%d0(ix:iz)*dnXref(3,(ic+2)/3)%d0
        END DO
        CALL dealloc_MatOFdnS(dnXref)

        A1 = matmul(transpose(A),A)
        A2 = matmul(A,transpose(A))

        IF (debug) THEN
          write(out_unitp,*) 'A'
          CALL Write_Mat(A,out_unitp,3)
          write(out_unitp,*) 'A1'
          CALL Write_Mat(A1,out_unitp,3)
          write(out_unitp,*) 'A2'
          CALL Write_Mat(A2,out_unitp,3)
        END IF

        CALL diagonalization(A1,eig1,Vec1,3,1,1,.FALSE.) ! jacobi + sort
        CALL diagonalization(A2,eig2,Vec2,3,1,1,.FALSE.) ! jacobi + sort

        ! change the sign of eig2(:,i)
        DO i=1,2
          IF (dot_product(Vec1(:,i),Vec2(:,i)) < ZERO) Vec2(:,i) = -Vec2(:,i)
        END DO
        CALL calc_cross_product(Vec1(:,1),normx,Vec1(:,2),normy,Vec1(:,3),normz)
        CALL calc_cross_product(Vec2(:,1),normx,Vec2(:,2),normy,Vec2(:,3),normz)

        IF (debug) THEN
          write(out_unitp,*) 'Vec1'
          CALL Write_Mat(Vec1,out_unitp,3)
          write(out_unitp,*) 'Vec2'
          CALL Write_Mat(Vec2,out_unitp,3)
        END IF

        T(1,1) = sum(Vec1(1,:)*Vec2(1,:))
        T(1,2) = sum(Vec1(1,:)*Vec2(2,:))
        T(1,3) = sum(Vec1(1,:)*Vec2(3,:))
        T(2,1) = sum(Vec1(2,:)*Vec2(1,:))
        T(2,2) = sum(Vec1(2,:)*Vec2(2,:))
        T(2,3) = sum(Vec1(2,:)*Vec2(3,:))
        T(3,1) = sum(Vec1(3,:)*Vec2(1,:))
        T(3,2) = sum(Vec1(3,:)*Vec2(2,:))
        T(3,3) = sum(Vec1(3,:)*Vec2(3,:))


        ! Check the rotational Eckart condition
        IF (debug) THEN
          d0x(:,:) = matmul(T,reshape(dnx%d0(1:CartesianTransfo%ncart_act), shape(d0x) ) )
          DO iref=1,CartesianTransfo%nb_RefGeometry
            Rot_Eckart(:) = ZERO
            DO i=1,CartesianTransfo%ncart_act/3

              Rot_Eckart(1) = Rot_Eckart(1)  +                          &
                             d0x(2,i)*CartesianTransfo%Qxyz(3,i,iref) - &
                             d0x(3,i)*CartesianTransfo%Qxyz(2,i,iref)

              Rot_Eckart(2) = Rot_Eckart(2) -                           &
                             d0x(1,i)*CartesianTransfo%Qxyz(3,i,iref) + &
                             d0x(3,i)*CartesianTransfo%Qxyz(1,i,iref)

              Rot_Eckart(3) = Rot_Eckart(3) +                           &
                             d0x(1,i)*CartesianTransfo%Qxyz(2,i,iref) - &
                             d0x(2,i)*CartesianTransfo%Qxyz(1,i,iref)
            END DO
            norm = sqrt(dot_product(Rot_Eckart(:),Rot_Eckart(:)))
            write(out_unitp,*) 'Norm Rot_Eckart',iref,Qact(1),norm
          END DO
        END IF

        IF (debug) write(out_unitp,*) 'Rot EC diff ID?',sum(abs(T(:,1)))+sum(abs(T(:,2)))+sum(abs(T(:,3)))-THREE

        IF (debug) THEN
          write(out_unitp,*) 'eig1 ',eig1
          write(out_unitp,*) 'eig2 ',eig2
          write(out_unitp,*) 'Eckart rotational matrix, T'
          CALL Write_Mat(T,out_unitp,3)
          write(out_unitp,*) 'END ',name_sub
        END IF

      END SUBROUTINE calc_EckartRot

      SUBROUTINE calc_dnTEckart(dnXin,dnT,dnXref,CartesianTransfo,nderiv)

        TYPE (Type_dnVec), intent(in)              :: dnXin
        TYPE(Type_dnS), intent(in)                 :: dnT(3,3)
        TYPE(Type_dnS), intent(in)                 :: dnXref(:,:)

        TYPE (Type_CartesianTransfo)   :: CartesianTransfo
        integer, intent(in)                        :: nderiv


        integer           :: i,j,ic

        ! from Dymarsky and Kudin ref JCP v122, p124103, 2005
        TYPE(Type_dnS) :: dnA(3,3),dntA(3,3),dnA1(3,3),dnA2(3,3)
        TYPE(Type_dnS) :: dnEig1(3),dnEig2(3),dnVec1(3,3),dnVec2(3,3)
        TYPE(Type_dnS) :: dnWork,dnWork2

        real (kind=Rkind) :: T(3,3),Xin(3,CartesianTransfo%nat_act),Vec1(3,3),Vec2(3,3)

        real (kind=Rkind) :: norm,dp(3),dp_save(3),sign_Vec(3)
        integer           :: isort_dp(3),idp(1)

        real(kind=Rkind) :: max_val
        integer          :: max_j
        logical          :: DymarskyKundin_only = .FALSE.  ! when false, the new way


!----- for debuging --------------------------------------------------
        integer :: nderiv_debug = 1
        character (len=*), parameter :: name_sub='calc_dnTEckart'
        logical, parameter :: debug=.FALSE.
        !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
          write(out_unitp,*) 'dnXin'
          CALL write_dnx(1,CartesianTransfo%ncart_act,dnXin,nderiv_debug)
          CALL Write_CartesianTransfo(CartesianTransfo)
          CALL flush_perso(out_unitp)
        END IF
!-----------------------------------------------------------
        CALL check_alloc_dnVec(dnXin,'dnXin',name_sub)
        CALL check_alloc_MatOFdnS(dnT,'dnT',name_sub)

        !==============================================================
        ! allocation
        CALL alloc_MatOFdnS(dnA,dnXin%nb_var_deriv,nderiv)
        CALL alloc_MatOFdnS(dntA,dnXin%nb_var_deriv,nderiv)
        CALL alloc_MatOFdnS(dnA1,dnXin%nb_var_deriv,nderiv)
        CALL alloc_MatOFdnS(dnA2,dnXin%nb_var_deriv,nderiv)
        CALL alloc_MatOFdnS(dnVec1,dnXin%nb_var_deriv,nderiv)
        CALL alloc_MatOFdnS(dnVec2,dnXin%nb_var_deriv,nderiv)

        CALL alloc_VecOFdnS(dnEig1,dnXin%nb_var_deriv,nderiv)
        CALL alloc_VecOFdnS(dnEig2,dnXin%nb_var_deriv,nderiv)

        CALL alloc_dnS(dnWork,dnXin%nb_var_deriv,nderiv)
        CALL alloc_dnS(dnWork2,dnXin%nb_var_deriv,nderiv)
        !==============================================================


        !==============================================================
        ! from Dymarsky and Kudin, JCP v122, p124103, 2005
        !write(out_unitp,*) 'Eckart',CartesianTransfo%Eckart

        CALL sub_ZERO_TO_MatOFdnS(dnA)
        DO ic=1,CartesianTransfo%ncart_act,3 ! loop on atoms
        DO i=1,3
          CALL sub_dnVec_TO_dnS(dnXin,dnWork,ic-1+i,nderiv)
          Xin(i,(ic+2)/3) = dnWork%d0

          DO j=1,3
            CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnWork,dnXref(j,(ic+2)/3),dnWork2,nderiv)
            CALL sub_dnS1_PLUS_dnS2_TO_dnS2(dnWork2,dnA(i,j),nderiv)
          END DO
        END DO
        END DO
        IF (debug) THEN
          write(out_unitp,*) 'dnA'
          CALL Write_MatOFdnS(dnA,nderiv=0)
        END IF

        CALL TRANS_Mat1OFdnS_TO_Mat2OFdnS(dnA,dntA,nderiv)

        CALL Mat1OFdnS_MUL_Mat2OFdnS_TO_Mat3OFdnS(dntA,dnA,dnA1,nderiv)
        CALL Mat1OFdnS_MUL_Mat2OFdnS_TO_Mat3OFdnS(dnA,dntA,dnA2,nderiv)

        IF (debug) THEN
          write(out_unitp,*) 'dnA1'
          CALL Write_MatOFdnS(dnA1,nderiv=0)

          write(out_unitp,*) 'dnA2'
          CALL Write_MatOFdnS(dnA2,nderiv=0)
        END IF

        CALL DIAG_MatOFdnS(dnA1,dnVec1,nderiv,sort=1,phase=.TRUE.,      &
                                type_diago=CartesianTransfo%type_diago, &
                                type_cs=CartesianTransfo%type_cs)
        CALL DIAG_MatOFdnS(dnA2,dnVec2,nderiv,sort=1,phase=.TRUE.,      &
                                type_diago=CartesianTransfo%type_diago, &
                                type_cs=CartesianTransfo%type_cs)

        IF (debug) THEN
          write(out_unitp,*) 'eig1 ',(dnA1(i,i)%d0,i=1,3)
          write(out_unitp,*) 'eig2 ',(dnA2(i,i)%d0,i=1,3)
          CALL flush_perso(out_unitp)
          write(out_unitp,*) 'dnVec1 before change sign (in column)'
          CALL Write_MatOFdnS(dnVec1,nderiv=0)
          write(out_unitp,*) 'dnVec2 before change sign (in column)'
          CALL Write_MatOFdnS(dnVec2,nderiv=0)
        END IF

        DO i=1,2
          dp(i) = dot_product(dnVec1(:,i)%d0,dnVec2(:,i)%d0)
          IF (dp(i) < ZERO) THEN ! change sign of dnVec2
            !write(6,*) 'change sign of vec: ',isort_dp(i)
            DO j=1,3
              CALL sub_dnS1_PROD_w_TO_dnS2(dnVec2(j,i),-ONE,dnVec2(j,i),nderiv)
            END DO
          END IF
        END DO

        ! Cross product for the third vector
        CALL Vec1OFdnS_CROSSPRODUCT_Vec2OFdnS_TO_Vec3OFdnS(         &
                         dnVec1(:,1),dnVec1(:,2),dnVec1(:,3),nderiv)
        CALL Vec1OFdnS_CROSSPRODUCT_Vec2OFdnS_TO_Vec3OFdnS(         &
                         dnVec2(:,1),dnVec2(:,2),dnVec2(:,3),nderiv)

        IF (.NOT. DymarskyKundin_only) THEN
          sign_Vec(:) = ONE
          CALL calc_RMS_d0Xout(sign_Vec,Xin,dnVec1(:,:)%d0,dnVec2(:,:)%d0,&
                                  dnXref(:,:)%d0,CartesianTransfo%nat_act)
          !write(out_unitp,*) 'sign_Vec',sign_Vec(:)

          DO i=1,3
            IF (sign_Vec(i) == -ONE) THEN
              DO j=1,3
                CALL sub_dnS1_PROD_w_TO_dnS2(dnVec2(j,i),sign_Vec(i),dnVec2(j,i),nderiv)
              END DO
            END IF
          END DO
        END IF


        IF (debug) THEN
          write(out_unitp,*) 'dnVec1 (in column)'
          CALL Write_MatOFdnS(dnVec1,nderiv=0)
          write(out_unitp,*) 'dnVec2 (in column)'
          CALL Write_MatOFdnS(dnVec2,nderiv=0)
        END IF
        !==============================================================

        DO i=1,3
        DO j=1,3
          CALL Vec1OFdnS_DOTPRODUCT_Vec2OFdnS_TO_dnS3(dnVec1(i,:),  &
                                        dnVec2(j,:),dnT(i,j),nderiv)
        END DO
        END DO
        !write(out_unitp,*) 'T',dnT(:,:)%d0

        IF (debug .OR. CartesianTransfo%check_dnT) THEN
          CALL DET_Mat3x3OFdnS_TO_dnS(dnT,dnWork,nderiv)
          IF (debug) THEN
            write(out_unitp,*) 'det(dnT)'
            CALL Write_dnS(dnWork)

            IF (abs(dnWork%d0-ONE) > ONETENTH**6) THEN
              write(out_unitp,*) 'ERROR in ',name_sub
              write(out_unitp,*) ' The determinant of T is different from ONE'
              write(out_unitp,*) ' It should never append!!'
              write(out_unitp,*) ' Check the fortran'
              STOP
            END IF
          END IF

!$OMP CRITICAL (calc_dnTEckart_CRIT)
          CartesianTransfo%dnTErr(0) =                                  &
                      max(CartesianTransfo%dnTErr(0),abs(dnWork%d0-ONE))
          IF (dnWork%nderiv > 0) THEN
             CartesianTransfo%dnTErr(1) =                               &
                  max(CartesianTransfo%dnTErr(1),maxval(abs(dnWork%d1)))
          END IF
          IF (dnWork%nderiv > 1) THEN
             CartesianTransfo%dnTErr(2) =                               &
                  max(CartesianTransfo%dnTErr(2),maxval(abs(dnWork%d2)))
          END IF
          IF (dnWork%nderiv > 2) THEN
             CartesianTransfo%dnTErr(3) =                               &
                  max(CartesianTransfo%dnTErr(3),maxval(abs(dnWork%d3)))
          END IF
!$OMP END CRITICAL (calc_dnTEckart_CRIT)


          IF (debug) THEN
            write(out_unitp,*) 'Eckart rotational matrix, T + det(T)',dnWork%d0
            CALL Write_MatOFdnS(dnT,nderiv=nderiv_debug)

            write(out_unitp,*) 'END ',name_sub
            CALL flush_perso(out_unitp)
          END IF
        END IF


        !==============================================================
        ! deallocation
        CALL dealloc_MatOFdnS(dnA)
        CALL dealloc_MatOFdnS(dntA)
        CALL dealloc_MatOFdnS(dnA1)
        CALL dealloc_MatOFdnS(dnA2)
        CALL dealloc_MatOFdnS(dnVec1)
        CALL dealloc_MatOFdnS(dnVec2)

        CALL dealloc_VecOFdnS(dnEig1)
        CALL dealloc_VecOFdnS(dnEig2)

        CALL dealloc_dnS(dnWork)
        CALL dealloc_dnS(dnWork2)
        !==============================================================


      END SUBROUTINE calc_dnTEckart

      RECURSIVE SUBROUTINE dnX_MultiRef(dnXref,CartesianTransfo,Qact,dnx,iref)

        TYPE (Type_dnS), intent(inout)  :: dnXref(:,:)

        TYPE (Type_CartesianTransfo), intent(in) :: CartesianTransfo
        TYPE (Type_dnVec), intent(in)            :: dnx
        real (kind=Rkind), intent(in)            :: Qact(:)
        integer, optional                        :: iref


        TYPE (Type_dnS)                 :: dnSwitch(CartesianTransfo%nb_RefGeometry)
        integer                         :: iref_loc,iat,i,iQact



!----- for debuging --------------------------------------------------
        character (len=*), parameter :: name_sub='dnX_MultiRef'
        logical, parameter :: debug=.FALSE.
!        logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
          DO iref_loc=1,CartesianTransfo%nb_RefGeometry
            write(out_unitp,*) 'Qxyz ',CartesianTransfo%Qxyz(:,:,iref_loc)
          END DO
          CALL flush_perso(out_unitp)
        END IF
!-----------------------------------------------------------
        iref_loc = 1
        IF (present(iref)) iref_loc = iref

        CALL sub_ZERO_TO_MatOFdnS(dnXref)


        IF (CartesianTransfo%nb_RefGeometry == 1 .OR. present(iref) ) THEN

          IF (iref_loc > CartesianTransfo%nb_RefGeometry) THEN
            write(out_unitp,*) 'ERROR in ',name_sub
            write(out_unitp,*) ' iref is larger than nb_RefGeometry',   &
                               iref_loc,CartesianTransfo%nb_RefGeometry
            write(out_unitp,*) ' It should never append!!'
            write(out_unitp,*) ' Check the fortran'
            STOP
          END IF

          DO iat=1,CartesianTransfo%nat_act
          DO i=1,3
            dnXref(i,iat)%d0 = CartesianTransfo%Qxyz(i,iat,iref_loc)
          END DO
          END DO

        ELSE

          CALL alloc_VecOFdnS(dnSwitch,dnXref(1,1)%nb_var_deriv,dnXref(1,1)%nderiv)

          SELECT CASE (CartesianTransfo%MultiRef_type)
          CASE (0) ! original periodic (do not use)
            iQact = CartesianTransfo%list_Qact(1)
            CALL Switch_type0(dnSwitch,CartesianTransfo,Qact(iQact),iQact)
          CASE (1) ! original periodic
            iQact = CartesianTransfo%list_Qact(1)
            CALL Switch_type1(dnSwitch,CartesianTransfo,Qact(iQact),iQact)
          CASE (2) ! along coordinates (low barrier)
            STOP 'case2 not yet'
          CASE (3) ! using distances (high barrier)
            CALL Switch_type3(dnSwitch,CartesianTransfo,dnx)
          CASE (4) ! using distances (high barrier)
            CALL Switch_type4(dnSwitch,CartesianTransfo,dnx)
          CASE Default ! using distances (high barrier)
            CALL Switch_type3(dnSwitch,CartesianTransfo,dnx)
          END SELECT

          DO iref_loc=1,CartesianTransfo%nb_RefGeometry
          DO iat=1,CartesianTransfo%nat_act
          DO i=1,3
            CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnSwitch(iref_loc),        &
                                CartesianTransfo%Qxyz(i,iat,iref_loc),  &
                                       dnXref(i,iat),ONE,dnXref(i,iat))
          END DO
          END DO
          END DO
          !write(99,*) Qact,dnSwitch(:)%d0,sum(dnSwitch(:)%d0)


          CALL dealloc_VecOFdnS(dnSwitch)


        END IF

!-----------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'dnXref'
          CALL Write_MatOFdnS(dnXref)
          write(out_unitp,*) 'END ',name_sub
          CALL flush_perso(out_unitp)
        END IF

      END SUBROUTINE dnX_MultiRef

      SUBROUTINE Switch_type0(dnSwitch,CartesianTransfo,Qact,iQact)

        TYPE (Type_dnS), intent(inout)           :: dnSwitch(:)

        TYPE (Type_CartesianTransfo), intent(in) :: CartesianTransfo
        real (kind=Rkind), intent(in)            :: Qact
        integer                                  :: iQact


        TYPE (Type_dnS)          :: dnQact
        TYPE (Type_dnS)          :: dnW1,dnW2

        integer                  :: iref,iat,i
        real (kind=Rkind)        :: cte(20)


!----- for debuging --------------------------------------------------
        character (len=*), parameter :: name_sub='Switch_type0'
        logical, parameter :: debug=.FALSE.
!        logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
          write(out_unitp,*) 'iQact ',iQact
          CALL flush_perso(out_unitp)
        END IF
!-----------------------------------------------------------

        CALL alloc_dnS(dnQact,dnSwitch(1)%nb_var_deriv,dnSwitch(1)%nderiv)
        CALL alloc_dnS(dnW1,  dnSwitch(1)%nb_var_deriv,dnSwitch(1)%nderiv)
        CALL alloc_dnS(dnW2,  dnSwitch(1)%nb_var_deriv,dnSwitch(1)%nderiv)


        CALL sub_ZERO_TO_dnS(dnQact)
        dnQact%d0 = Qact
        IF (dnQact%nderiv > 0) dnQact%d1(iQact) = ONE

        cte(:) =  ZERO
        cte(1) =  ONE
        cte(2) = -CartesianTransfo%tab_phi0(1)
        CALL sub_dnS1_TO_dntR2(dnQact,dnW1,100,cte=cte) ! x=>x-phi0

        CALL sub_dnS1_TO_dntR2(dnW1,dnW2,2) ! x=>cos(x-phi0)


        cte(:) =  ZERO
        cte(1) =  CartesianTransfo%tab_sc(1)
        CALL sub_dnS1_TO_dntR2(dnW2,dnW1,100,cte=cte) ! x=>sc*cos(x-phi0)

        CALL sub_dnS1_TO_dntR2(dnW1,dnSwitch(1),72) ! x=>(1+tanh(x))/2 : x= (1+tanh(sc*cos(x-phi0)))/2


        cte(:) = ZERO ; cte(1) = -ONE ; cte(2) = ONE
        CALL sub_dnS1_TO_dntR2(dnSwitch(1),dnSwitch(2),100,cte=cte) ! x=>1-x : : x= 1-(1+tanh(sc*cos(x-phi0)))/2

        CALL dealloc_dnS(dnQact)
        CALL dealloc_dnS(dnW1)
        CALL dealloc_dnS(dnW2)

!-----------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'dnSwitch'
          CALL Write_VecOFdnS(dnSwitch)
          write(out_unitp,*) 'END ',name_sub
          CALL flush_perso(out_unitp)
        END IF

      END SUBROUTINE Switch_type0

      SUBROUTINE Switch_type1(dnSwitch,CartesianTransfo,Qact,iQact)

        TYPE (Type_dnS), intent(inout) :: dnSwitch(:)

        TYPE (Type_CartesianTransfo), intent(in) :: CartesianTransfo
        real (kind=Rkind), intent(in)            :: Qact
        integer, intent(in)                      :: iQact


        TYPE (Type_dnS)          :: dnQact

        TYPE (Type_dnS)          :: dnW1,dnW2

        integer              :: iref,iat,i
        real (kind=Rkind)    :: cte(20),phase,a,b,coef


!----- for debuging --------------------------------------------------
        character (len=*), parameter :: name_sub='Switch_type1'
        logical, parameter :: debug=.FALSE.
!        logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
          write(out_unitp,*) 'iQact ',iQact
          CALL flush_perso(out_unitp)
        END IF
!-----------------------------------------------------------

        CALL alloc_dnS(dnQact,dnSwitch(1)%nb_var_deriv,dnSwitch(1)%nderiv)

        CALL sub_ZERO_TO_dnS(dnQact)
        dnQact%d0 = Qact
        IF (dnQact%nderiv > 0) dnQact%d1(iQact) = ONE


        SELECT CASE (CartesianTransfo%tab_switch_type(1))
        CASE (1) ! High barier
          CALL alloc_dnS(dnW1,  dnSwitch(1)%nb_var_deriv,dnSwitch(1)%nderiv)
          CALL alloc_dnS(dnW2,  dnSwitch(1)%nb_var_deriv,dnSwitch(1)%nderiv)

          phase = TWO*pi/real(CartesianTransfo%nb_RefGeometry,kind=Rkind)
          ! a*(cos(x)-b): b=cos(phase/2)
          !   => a*(cos(x)-b)=0 if x=phase/2 (in the middle of two reference geometries, x=0 and x=phase)
          b  = cos(phase*HALF)
          a  = ONE/(ONE-b)
          !write(6,*) 'a,b',a,b

          DO iref=1,CartesianTransfo%nb_RefGeometry
            cte(:) =  ZERO
            cte(1) =  ONE
            cte(2) = -CartesianTransfo%tab_phi0(1) + phase*real(iref-1,kind=Rkind)
            CALL sub_dnS1_TO_dntR2(dnQact,dnW1,100,cte=cte) ! x=>x-phi0

            CALL sub_dnS1_TO_dntR2(dnW1,dnW2,2) ! x=>cos(x-phi0)


            cte(:) =  ZERO
            cte(1) =  CartesianTransfo%tab_sc(1)*a
            cte(2) = -cte(1)*b
            CALL sub_dnS1_TO_dntR2(dnW2,dnW1,100,cte=cte) ! x=>sc*[a * (cos(x-phi0) - b)]

            CALL sub_dnS1_TO_dntR2(dnW1,dnSwitch(iref),72) ! x=>(1+tanh(x))/2 : x= (1+tanh(sc*cos(x-phi0)))/2
          END DO

          CALL dealloc_dnS(dnW1)
          CALL dealloc_dnS(dnW2)

        CASE (2)  ! Low barier
          CALL alloc_dnS(dnW1,  dnSwitch(1)%nb_var_deriv,dnSwitch(1)%nderiv)
          CALL alloc_dnS(dnW2,  dnSwitch(1)%nb_var_deriv,dnSwitch(1)%nderiv)

          !  2/n  {sum_i=1,n sin((phi-phi_ref1)/2+i·pi/n)^2 · RefGeom_i} (for n>=2)
          !  8/3n {sum_i=1,n sin((phi-phi_ref1)/2+i·pi/n)^4 · RefGeom_i} (for n>=3)

          phase = TWO*pi/real(CartesianTransfo%nb_RefGeometry,kind=Rkind)

          IF (CartesianTransfo%nb_RefGeometry == 2)                     &
             coef = TWO/real(CartesianTransfo%nb_RefGeometry,kind=Rkind)
          IF (CartesianTransfo%nb_RefGeometry == 4)                     &
             coef = EIGHT/THREE / real(CartesianTransfo%nb_RefGeometry,kind=Rkind)
          IF (CartesianTransfo%nb_RefGeometry == 6)                     &
             coef = TWO**4/FIVE / real(CartesianTransfo%nb_RefGeometry,kind=Rkind)

          DO iref=1,CartesianTransfo%nb_RefGeometry

            cte(:) =  ZERO
            cte(1) =  TWO
            cte(2) = CartesianTransfo%tab_phi0(1) - phase*real(iref,kind=Rkind)
            CALL sub_dnS1_TO_dntR2(dnQact,dnW1,-100,cte=cte) ! x=>(x-cte2)/cte1

            CALL sub_dnS1_TO_dntR2(dnW1,dnW2,3) ! x=>sin(x-phi0)

            cte(:) =  ZERO
            cte(1) =  real(CartesianTransfo%tab_expo(1),kind=Rkind)
            CALL sub_dnS1_TO_dntR2(dnW2,dnW1,99,cte=cte) ! x=>x^tab_expo(1)

            cte(:) =  ZERO
            cte(1) =  coef
            CALL sub_dnS1_TO_dntR2(dnW1,dnSwitch(iref),100,cte=cte) ! x=> coef * x

          END DO

          CALL dealloc_dnS(dnW1)
          CALL dealloc_dnS(dnW2)

        CASE DEFAULT
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' Wrong value of CartesianTransfo%tab_switch_type',&
                                   CartesianTransfo%tab_switch_type(:)
          write(out_unitp,*) ' check the Fortran!'
          STOP
        END SELECT

        CALL dealloc_dnS(dnQact)


!-----------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'dnSwitch'
          CALL Write_VecOFdnS(dnSwitch)
          write(out_unitp,*) 'END ',name_sub
          CALL flush_perso(out_unitp)
        END IF

      END SUBROUTINE Switch_type1

      SUBROUTINE Switch_type3(dnSwitch,CartesianTransfo,dnx)

        TYPE (Type_dnS), intent(inout) :: dnSwitch(:)

        TYPE (Type_dnVec), intent(in)            :: dnx
        TYPE (Type_CartesianTransfo), intent(in) :: CartesianTransfo


        TYPE (Type_dnS), pointer :: dnDist2(:)

        TYPE (Type_dnS)          :: dnW1,dnW2,dnSumExp

        integer              :: iref,kref,iat,ic,ixyz
        real (kind=Rkind)    :: cte(20),sc


!----- for debuging --------------------------------------------------
        character (len=*), parameter :: name_sub='Switch_type3'
        logical, parameter :: debug=.FALSE.
!        logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
          DO iref=1,CartesianTransfo%nb_RefGeometry
            write(out_unitp,*) 'Qxyz ',CartesianTransfo%Qxyz(:,:,iref)
          END DO
          write(out_unitp,*) 'dnx%d0 ',dnx%d0(:)
          CALL flush_perso(out_unitp)
        END IF
!-----------------------------------------------------------

        !---------------------------------------------------------------
        ! allocation
        nullify(dnDist2)
        CALL alloc_array(dnDist2,(/CartesianTransfo%nb_RefGeometry/),    &
                        "dnDist2",name_sub)
        CALL alloc_VecOFdnS(dnDist2,dnSwitch(1)%nb_var_deriv,dnSwitch(1)%nderiv)

        CALL alloc_dnS(dnW1,    dnSwitch(1)%nb_var_deriv,dnSwitch(1)%nderiv)
        CALL alloc_dnS(dnW2,    dnSwitch(1)%nb_var_deriv,dnSwitch(1)%nderiv)
        CALL alloc_dnS(dnSumExp,dnSwitch(1)%nb_var_deriv,dnSwitch(1)%nderiv)

        !---------------------------------------------------------------

        sc = CartesianTransfo%tab_sc(1) ! all of them must be identical

        DO iref=1,CartesianTransfo%nb_RefGeometry

          CALL sub_ZERO_TO_dnS(dnDist2(iref))

          DO ic=1,CartesianTransfo%ncart_act
            ! first extrac cartesian CC, ic, of atom iat from dnxEC
            CALL sub_dnVec_TO_dnS(dnx,dnW1,ic)
            iat = (ic+2)/3
            ixyz = mod(ic-1,3)+1
            dnW1%d0 = dnW1%d0 - CartesianTransfo%Qxyz(ixyz,iat,iref) ! x-xref

            CALL sub_dnS1_PROD_w_TO_dnS2(dnW1,ONE/CartesianTransfo%d0sm(ic),dnW1) ! remove the mass-weighted

            CALL sub_dnS1_TO_dntR2(dnW1,dnW2,-91) ! (x-xref)^2

            CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnW2,ONE,dnDist2(iref),ONE,&
                                             dnDist2(iref))             ! add (x-xref)^2 to dnDist2
          END DO
          CALL sub_dnS1_PROD_w_TO_dnS2(dnDist2(iref),                   &
            ONE/real(CartesianTransfo%nat_act,kind=Rkind),dnDist2(iref))  ! divide by nat_act

        END DO
        IF (debug) write(out_unitp,*) 'dnDist2',dnDist2(:)%d0
        !write(98,*) dnDist2(:)%d0

        DO iref=1,CartesianTransfo%nb_RefGeometry

          CALL sub_ZERO_TO_dnS(dnSumExp) ! the sum of the exp
          dnSumExp%d0 = ONE ! because the exp with kref = iref is not the next loop

          DO kref=1,CartesianTransfo%nb_RefGeometry
            IF (iref == kref) CYCLE
            CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnDist2(iref), sc,         &
                                             dnDist2(kref),-sc,dnW1)
            IF (debug) write(out_unitp,*) 'iref,kref,DeltaDist2',iref,kref,dnW1%d0

            cte(:) = ZERO ; cte(1) = ONE
            CALL sub_dnS1_TO_dntR2(dnW1,dnW2,80,cte=cte) ! exp(sc*(dist2_i-dist2_k))
            IF (debug) write(out_unitp,*) 'iref,kref,dnExp',iref,kref,dnW2%d0
            CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnW2,ONE,dnSumExp,ONE,     &
                                             dnSumExp)             ! sum of the exp
          END DO
          CALL sub_dnS1_TO_dntR2(dnSumExp,dnSwitch(iref),90) ! 1/sum(exp ....)

        END DO
        IF (debug) write(out_unitp,*) 'dnSwitch',dnSwitch(:)%d0

        !---------------------------------------------------------------
        ! deallocation
        CALL dealloc_dnS(dnW1)
        CALL dealloc_dnS(dnW2)
        CALL dealloc_dnS(dnSumExp)

        CALL dealloc_VecOFdnS(dnDist2)
        CALL dealloc_array(dnDist2,"dnDist2",name_sub)
        !---------------------------------------------------------------
!stop
!-----------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'dnSwitch'
          CALL Write_VecOFdnS(dnSwitch)
          write(out_unitp,*) 'END ',name_sub
          CALL flush_perso(out_unitp)
        END IF

      END SUBROUTINE Switch_type3

      SUBROUTINE Switch_type4(dnSwitch,CartesianTransfo,dnx)

        TYPE (Type_dnS), intent(inout)           :: dnSwitch(:)
        TYPE (Type_dnVec), intent(in)            :: dnx
        TYPE (Type_CartesianTransfo), intent(in) :: CartesianTransfo


        TYPE (Type_dnS), pointer :: dnDist2(:)

        TYPE (Type_dnS)          :: dnW1,dnW2,dnSumExp
        TYPE (Type_dnS)          :: dnXref(3,CartesianTransfo%nat_act)
        TYPE(Type_dnS)           :: dnT(3,3)
        TYPE (Type_dnVec)        :: dnxEC

        integer              :: iref,kref,iat,ic,ixyz
        real (kind=Rkind)    :: cte(20),sc
        real (kind=Rkind)    :: Qact(CartesianTransfo%nb_Qact)


!----- for debuging --------------------------------------------------
        character (len=*), parameter :: name_sub='Switch_type4'
        logical, parameter :: debug=.FALSE.
        !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
          DO iref=1,CartesianTransfo%nb_RefGeometry
            write(out_unitp,*) 'Qxyz ',CartesianTransfo%Qxyz(:,:,iref)
          END DO
          write(out_unitp,*) 'dnx%d0 ',dnx%d0(:)
          CALL flush_perso(out_unitp)
        END IF
!-----------------------------------------------------------

        !---------------------------------------------------------------
        ! allocation
        nullify(dnDist2)
        CALL alloc_array(dnDist2,(/CartesianTransfo%nb_RefGeometry/),    &
                        "dnDist2",name_sub)
        CALL alloc_VecOFdnS(dnDist2,dnSwitch(1)%nb_var_deriv,dnSwitch(1)%nderiv)

        CALL alloc_dnS(dnW1,    dnSwitch(1)%nb_var_deriv,dnSwitch(1)%nderiv)
        CALL alloc_dnS(dnW2,    dnSwitch(1)%nb_var_deriv,dnSwitch(1)%nderiv)
        CALL alloc_dnS(dnSumExp,dnSwitch(1)%nb_var_deriv,dnSwitch(1)%nderiv)

        CALL alloc_MatOFdnS(dnXref,dnSwitch(1)%nb_var_deriv,dnSwitch(1)%nderiv)
        CALL alloc_MatOFdnS(dnT,dnSwitch(1)%nb_var_deriv,dnSwitch(1)%nderiv)
        CALL alloc_dnSVM(dnxEC,dnx%nb_var_vec,dnx%nb_var_deriv,dnx%nderiv)

        !---------------------------------------------------------------

        sc = CartesianTransfo%tab_sc(1) ! all of them must be identical

        DO iref=1,CartesianTransfo%nb_RefGeometry

          Qact(:) = ZERO
          CALL dnX_MultiRef(dnXref,CartesianTransfo,Qact,dnx,iref)
          CALL calc_dnTEckart(dnx,dnT,dnXref,CartesianTransfo,dnSwitch(1)%nderiv)
          CALL calc_dnTxdnXin_TO_dnXout(dnx,dnT,dnxEC,CartesianTransfo,dnSwitch(1)%nderiv)

          CALL sub_ZERO_TO_dnS(dnDist2(iref))

          DO ic=1,CartesianTransfo%ncart_act
            ! first extrac cartesian CC, ic, of atom iat from dnxEC
            CALL sub_dnVec_TO_dnS(dnxEC,dnW1,ic)
            iat = (ic+2)/3
            ixyz = mod(ic-1,3)+1

            dnW1%d0 = dnW1%d0 - CartesianTransfo%Qxyz(ixyz,iat,iref) ! x-xref

            CALL sub_dnS1_PROD_w_TO_dnS2(dnW1,ONE/CartesianTransfo%d0sm(ic),dnW1) ! remove the mass-weighted

            CALL sub_dnS1_TO_dntR2(dnW1,dnW2,-91) ! (x-xref)^2

            CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnW2,ONE,dnDist2(iref),ONE,&
                                             dnDist2(iref))             ! add (x-xref)^2 to dnDist2
          END DO
          CALL sub_dnS1_PROD_w_TO_dnS2(dnDist2(iref),                   &
            ONE/real(CartesianTransfo%nat_act,kind=Rkind),dnDist2(iref))  ! divide by nat_act
        END DO
        IF (debug) write(out_unitp,*) 'dnDist2',dnDist2(:)%d0
        !write(98,*) dnDist2(:)%d0

        DO iref=1,CartesianTransfo%nb_RefGeometry

          CALL sub_ZERO_TO_dnS(dnSumExp) ! the sum of the exp
          dnSumExp%d0 = ONE ! because the exp with kref = iref is not the next loop

          DO kref=1,CartesianTransfo%nb_RefGeometry
            IF (iref == kref) CYCLE
            CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnDist2(iref), sc,         &
                                             dnDist2(kref),-sc,dnW1)
            IF (debug) write(out_unitp,*) 'iref,kref,DeltaDist2',iref,kref,dnW1%d0

            cte(:) = ZERO ; cte(1) = ONE
            CALL sub_dnS1_TO_dntR2(dnW1,dnW2,80,cte=cte) ! exp(sc*(dist2_i-dist2_k))
            IF (debug) write(out_unitp,*) 'iref,kref,dnExp',iref,kref,dnW2%d0
            CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnW2,ONE,dnSumExp,ONE,     &
                                             dnSumExp)             ! sum of the exp
          END DO
          CALL sub_dnS1_TO_dntR2(dnSumExp,dnSwitch(iref),90) ! 1/sum(exp ....)

        END DO
        IF (debug) write(out_unitp,*) 'dnSwitch',dnSwitch(:)%d0

        !---------------------------------------------------------------
        ! deallocation
        CALL dealloc_MatOFdnS(dnXref)
        CALL dealloc_MatOFdnS(dnT)
        CALL dealloc_dnSVM(dnxEC)

        CALL dealloc_dnS(dnW1)
        CALL dealloc_dnS(dnW2)
        CALL dealloc_dnS(dnSumExp)

        CALL dealloc_VecOFdnS(dnDist2)
        CALL dealloc_array(dnDist2,"dnDist2",name_sub)
        !---------------------------------------------------------------

!-----------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'dnSwitch'
          CALL Write_VecOFdnS(dnSwitch)
          write(out_unitp,*) 'END ',name_sub
          CALL flush_perso(out_unitp)
        END IF

      END SUBROUTINE Switch_type4

!================================================================
!       Calcul le centre de masse G
!================================================================
      SUBROUTINE sub3_dncentre_masse(ncart_act,nb_act,ncart,            &
                                     dnx,                               &
                                     masses,d0sm,Mtot_inv,icG,          &
                                     nderiv)
      USE mod_system
      USE mod_dnSVM
      IMPLICIT NONE

        integer :: ncart_act,nb_act,ncart
        integer :: icG,nderiv
        real (kind=Rkind) ::  Mtot_inv,masses(ncart),d0sm(ncart)

        TYPE (Type_dnVec) :: dnx


        integer :: i,j,k


!======================================================================
!         G : Center of mass calculation
!         => molecule centers on G
!======================================================================

          CALL centre_masse(ncart_act,ncart,dnx%d0,            &
                            masses,Mtot_inv,icG)

          IF (nderiv >= 1) THEN
            DO i=1,nb_act
!           first derivatives
            CALL centre_masse(ncart_act,ncart,dnx%d1(1:ncart,i),        &
                              masses,Mtot_inv,icG)
            END DO
          END IF

          IF (nderiv >= 2) THEN
          DO i=1,nb_act
          DO j=1,nb_act
!         2d derivatives
          CALL centre_masse(ncart_act,ncart,dnx%d2(1:ncart,i,j),        &
                            masses,Mtot_inv,icG)
          END DO
          END DO
          END IF

          IF (nderiv >= 3) THEN
          DO i=1,nb_act
          DO j=1,nb_act
          DO k=1,nb_act
!         3d derivatives
          CALL centre_masse(ncart_act,ncart,dnx%d3(1:ncart,i,j,k),      &
                            masses,Mtot_inv,icG)
          END DO
          END DO
          END DO
          END IF
!======================================================================

      end subroutine sub3_dncentre_masse
!================================================================
!         G : masse center calculation
!         => molecule centers on G
!       can be use for d1G d2G ...
!================================================================
      SUBROUTINE centre_masse(ncart_act,ncart,d0x,                      &
                              masses,Mtot_inv,icG)
      USE mod_system
      IMPLICIT NONE

        integer :: ncart_act,ncart
        integer, optional :: icG
        real (kind=Rkind) :: Mtot_inv,masses(ncart)
        real (kind=Rkind) :: d0x(ncart)

        integer :: i,ic
        real (kind=Rkind) ::  xG(3)

        !write(6,*) 'ncart_act,ncart',ncart_act,ncart
        !write(6,*) 'Mtot_inv',Mtot_inv
        !write(6,*) 'masses',masses(:)
        !write(6,*) 'd0x',d0x

        DO i=1,3
          xG(i) = ZERO
          DO ic=i,ncart_act,3
            xG(i) = xG(i) + masses(ic)*d0x(ic)
          END DO
          xG(i) = xG(i) * Mtot_inv
          DO ic=i,ncart_act,3
             d0x(ic) = (d0x(ic) - xG(i))
          END DO
        END DO

        IF (present(icG)) d0x(icG+0:icG+2) = xG(:)

        !write(out_unitp,*) 'G =',xG(:)
        !write(6,*) 'd0x',d0x

      END SUBROUTINE centre_masse

      SUBROUTINE sub3_NOdncentre_masse(ncart_act,nb_act,ncart,dnx,icG,nderiv)
      USE mod_system
      USE mod_dnSVM
      IMPLICIT NONE

        integer :: ncart_act,nb_act,ncart
        integer :: icG,nderiv

        TYPE (Type_dnVec) :: dnx


        integer :: i,j,k


!======================================================================
!         G : Center of mass calculation
!         => molecule centers on G
!======================================================================

          CALL NOcentre_masse(ncart_act,ncart,dnx%d0,icG)

          IF (nderiv >= 1) THEN
            DO i=1,nb_act
!           first derivatives
            CALL NOcentre_masse(ncart_act,ncart,dnx%d1(1:ncart,i),icG)
            END DO
          END IF

          IF (nderiv >= 2) THEN
          DO i=1,nb_act
          DO j=1,nb_act
!         2d derivatives
          CALL NOcentre_masse(ncart_act,ncart,dnx%d2(1:ncart,i,j),icG)
          END DO
          END DO
          END IF

          IF (nderiv >= 3) THEN
          DO i=1,nb_act
          DO j=1,nb_act
          DO k=1,nb_act
!         3d derivatives
          CALL NOcentre_masse(ncart_act,ncart,dnx%d3(1:ncart,i,j,k),icG)
          END DO
          END DO
          END DO
          END IF
!======================================================================


      end subroutine sub3_NOdncentre_masse

!================================================================
!         G : masse center calculation
!         => molecule centers on G
!       can be use for d1G d2G ...
!================================================================
      SUBROUTINE NOcentre_masse(ncart_act,ncart,d0x,icG)
      USE mod_system
      IMPLICIT NONE

        integer :: ncart_act,ncart
        integer :: icG
        real (kind=Rkind) :: d0x(ncart)

        integer :: i,ic
        real (kind=Rkind) ::  xG(3)

        !write(6,*) 'ncart_act,ncart',ncart_act,ncart
        !write(6,*) 'd0x',d0x

        xG(:) = d0x(icG+0:icG+2)


        DO i=1,3
          DO ic=i,ncart_act,3
             d0x(ic) = d0x(ic) + xG(i)
          END DO
        END DO

        !write(out_unitp,*) 'G =',xG(:)
        !write(6,*) 'd0x',d0x

      END SUBROUTINE NOcentre_masse
!================================================================
!       Calcul le centre de masse G
!================================================================
      SUBROUTINE sub_dnxMassWeight(dnx,d0sm,ncart,ncart_act,nderiv)
      USE mod_system
      USE mod_dnSVM
      IMPLICIT NONE

      integer :: ncart_act,ncart
      integer :: nderiv
      real (kind=Rkind) :: d0sm(ncart)
      TYPE (Type_dnVec) :: dnx


        integer :: i,j,k

        dnx%d0(1:ncart_act) = dnx%d0(1:ncart_act) * d0sm(1:ncart_act)

        IF (nderiv >= 1) THEN
          DO i=1,dnx%nb_var_deriv
            ! first derivatives
            dnx%d1(1:ncart_act,i) = dnx%d1(1:ncart_act,i) * d0sm(1:ncart_act)
          END DO
        END IF

        IF (nderiv >= 2) THEN
          DO i=1,dnx%nb_var_deriv
          DO j=1,dnx%nb_var_deriv
            ! 2d derivatives
            dnx%d2(1:ncart_act,j,i) = dnx%d2(1:ncart_act,j,i) * d0sm(1:ncart_act)
          END DO
          END DO
        END IF

        IF (nderiv >= 3) THEN
          DO i=1,dnx%nb_var_deriv
          DO j=1,dnx%nb_var_deriv
          DO k=1,dnx%nb_var_deriv
            ! 2d derivatives
            dnx%d3(1:ncart_act,k,j,i) = dnx%d3(1:ncart_act,k,j,i) * d0sm(1:ncart_act)
          END DO
          END DO
          END DO
        END IF

      END SUBROUTINE sub_dnxMassWeight
      SUBROUTINE sub_dnxNOMassWeight(dnx,d0sm,ncart,ncart_act,nderiv)
      USE mod_system
      USE mod_dnSVM
      IMPLICIT NONE

      integer :: ncart_act,ncart
      integer :: nderiv
      real (kind=Rkind) :: d0sm(ncart)
      TYPE (Type_dnVec) :: dnx


        integer :: i,j,k

        dnx%d0(1:ncart_act) = dnx%d0(1:ncart_act) / d0sm(1:ncart_act)

        IF (nderiv >= 1) THEN
          DO i=1,dnx%nb_var_deriv
            ! first derivatives
            dnx%d1(1:ncart_act,i) = dnx%d1(1:ncart_act,i) / d0sm(1:ncart_act)
          END DO
        END IF

        IF (nderiv >= 2) THEN
          DO i=1,dnx%nb_var_deriv
          DO j=1,dnx%nb_var_deriv
            ! 2d derivatives
            dnx%d2(1:ncart_act,j,i) = dnx%d2(1:ncart_act,j,i) / d0sm(1:ncart_act)
          END DO
          END DO
        END IF

        IF (nderiv >= 3) THEN
          DO i=1,dnx%nb_var_deriv
          DO j=1,dnx%nb_var_deriv
          DO k=1,dnx%nb_var_deriv
            ! 2d derivatives
            dnx%d3(1:ncart_act,k,j,i) = dnx%d3(1:ncart_act,k,j,i) / d0sm(1:ncart_act)
          END DO
          END DO
          END DO
        END IF

      END SUBROUTINE sub_dnxNOMassWeight



      END MODULE mod_CartesianTransfo

