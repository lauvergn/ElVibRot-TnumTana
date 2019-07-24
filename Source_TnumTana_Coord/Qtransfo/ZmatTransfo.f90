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
      MODULE mod_ZmatTransfo
      use mod_system
      USE mod_dnSVM
      use mod_constant,     only: table_atom, get_mass_tnum
      use mod_Lib_QTransfo ! only all
      IMPLICIT NONE

      PRIVATE

      !!@description: TODO
      !!@param: TODO
      TYPE Type_ZmatTransfo

        integer           :: ncart=0,ncart_act=0
        integer           :: nat0=0,nat=0,nat_act=0
        integer           :: nb_var=0
        integer, pointer  :: ind2_zmat(:,:) => null()
        integer, pointer  :: ind_zmat(:,:) => null()
        logical           :: New_Orient = .FALSE. ! (F) T => Can use different orientation for the z-matrix
        real (kind=Rkind) :: vAt1(3)=ZERO,vAt2(3)=ZERO,vAt3(3)=ZERO
        logical           :: cos_th  = .FALSE. ! T => coordinate (valence angle) => cos(th)
                                                 ! F => coordinate (valence angle) => th

        ! just for read the input data
        real (kind=Rkind),        pointer :: masses(:)     => null()
        integer, pointer                  :: Z(:)          => null()
        character (len=Name_len),pointer  :: symbole(:)    => null()

        integer, pointer                  :: type_Qin(:)   => null() ! TRUE pointer
        character (len=Name_len), pointer :: name_Qin(:)   => null() ! TRUE pointer

      END TYPE Type_ZmatTransfo

      PUBLIC :: Type_ZmatTransfo, alloc_ZmatTransfo, dealloc_ZmatTransfo
      PUBLIC :: read_ZmatTransfo, Write_ZmatTransfo, calc_ZmatTransfo, calc_ZmatTransfo_outTOin
      PUBLIC :: ZmatTransfo1TOZmatTransfo2

      CONTAINS

!================================================================
!       Read Zmat Transfo
!
!       -analysis of the z-matrix or cart --------------
!       => type_Qin(i):  (before type_zmat)
!                       1 => cartesian
!                       2 => distance
!                       -3 => cos(valence angle)
!                       3 => valence angle
!                       4 => dihedral angle
!================================================================
      !!@description: Read Zmat Transfo
      !!
      !!       -analysis of the z-matrix or cart --------------
      !!       => type_Qin(i):  (before type_zmat)
      !!                       1 => cartesian
      !!                       2 => distance
      !!                       -3 => cos(valence angle)
      !!                       3 => valence angle
      !!                       4 => dihedral angle
      !!@param: TODO
      SUBROUTINE alloc_ZmatTransfo(ZmatTransfo)
      TYPE (Type_ZmatTransfo), intent(inout) :: ZmatTransfo

!      write(out_unitp,*) 'BEGINNING alloc_ZmatTransfo'
!      write(out_unitp,*) 'nat',ZmatTransfo%nat

       IF (ZmatTransfo%nat < 3) THEN
         write(out_unitp,*) ' ERROR in alloc_ZmatTransfo'
         write(out_unitp,*) ' wrong value of nat',ZmatTransfo%nat
         write(out_unitp,*) ' CHECK the source !!'
         STOP
       END IF

       IF (associated(ZmatTransfo%ind2_zmat))  THEN
         CALL dealloc_array(ZmatTransfo%ind2_zmat,                      &
                           "ZmatTransfo%ind2_zmat","alloc_ZmatTransfo")
       END IF
       CALL alloc_array(ZmatTransfo%ind2_zmat,(/5,ZmatTransfo%nat/),    &
                       "ZmatTransfo%ind2_zmat","alloc_ZmatTransfo")
       ZmatTransfo%ind2_zmat(:,:) = 0

       IF (associated(ZmatTransfo%ind_zmat))   THEN
         CALL dealloc_array(ZmatTransfo%ind_zmat,                       &
                           "ZmatTransfo%ind_zmat","alloc_ZmatTransfo")
       END IF
       CALL alloc_array(ZmatTransfo%ind_zmat,(/5,ZmatTransfo%nat/),     &
                       "ZmatTransfo%ind_zmat","alloc_ZmatTransfo")
       ZmatTransfo%ind_zmat(:,:) = 0


       IF (associated(ZmatTransfo%Z))  THEN
         CALL dealloc_array(ZmatTransfo%Z,                              &
                           "ZmatTransfo%Z","alloc_ZmatTransfo")
       END IF
       CALL alloc_array(ZmatTransfo%Z,(/ZmatTransfo%nat/),              &
                       "ZmatTransfo%Z","alloc_ZmatTransfo")
       ZmatTransfo%Z(:) = 0

       IF (associated(ZmatTransfo%masses))  THEN
         CALL dealloc_array(ZmatTransfo%masses,                         &
                           "ZmatTransfo%masses","alloc_ZmatTransfo")
       END IF
       CALL alloc_array(ZmatTransfo%masses,(/ZmatTransfo%ncart/), &
                       "ZmatTransfo%masses","alloc_ZmatTransfo")
       ZmatTransfo%masses(:) = ZERO

       IF (associated(ZmatTransfo%symbole))  THEN
         CALL dealloc_array(ZmatTransfo%symbole,                        &
                           "ZmatTransfo%symbole","alloc_ZmatTransfo")
       END IF
       CALL alloc_array(ZmatTransfo%symbole,(/ZmatTransfo%nat/),  &
              Name_len,"ZmatTransfo%symbole","alloc_ZmatTransfo")
       ZmatTransfo%symbole(:) = ""


!      write(out_unitp,*) 'END alloc_ZmatTransfo'

      END SUBROUTINE alloc_ZmatTransfo

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE dealloc_ZmatTransfo(ZmatTransfo)

       TYPE (Type_ZmatTransfo), intent(inout) :: ZmatTransfo

       !write(out_unitp,*) 'BEGINNING dealloc_ZmatTransfo'; call flush_perso(out_unitp)

       IF (associated(ZmatTransfo%ind2_zmat))  THEN
         CALL dealloc_array(ZmatTransfo%ind2_zmat,                      &
                           "ZmatTransfo%ind2_zmat","dealloc_ZmatTransfo")
       END IF

       IF (associated(ZmatTransfo%ind_zmat))   THEN
         CALL dealloc_array(ZmatTransfo%ind_zmat,                       &
                           "ZmatTransfo%ind_zmat","dealloc_ZmatTransfo")
       END IF

       IF (associated(ZmatTransfo%Z))  THEN
         CALL dealloc_array(ZmatTransfo%Z,                              &
                           "ZmatTransfo%Z","dealloc_ZmatTransfo")
       END IF

       IF (associated(ZmatTransfo%masses))  THEN
         CALL dealloc_array(ZmatTransfo%masses,                         &
                           "ZmatTransfo%masses","dealloc_ZmatTransfo")
       END IF

       IF (associated(ZmatTransfo%symbole))  THEN
         CALL dealloc_array(ZmatTransfo%symbole,                        &
                           "ZmatTransfo%symbole","dealloc_ZmatTransfo")
       END IF

        ZmatTransfo%ncart     = 0
        ZmatTransfo%ncart_act = 0
        ZmatTransfo%nat0      = 0
        ZmatTransfo%nat       = 0
        ZmatTransfo%nat_act   = 0
        ZmatTransfo%nb_var    = 0

        ZmatTransfo%New_Orient  = .FALSE.
        ZmatTransfo%vAt1(:)     = ZERO
        ZmatTransfo%vAt2(:)     = ZERO
        ZmatTransfo%vAt3(:)     = ZERO

        ZmatTransfo%cos_th      = .FALSE.


       nullify(ZmatTransfo%type_Qin)
       nullify(ZmatTransfo%name_Qin)

       !write(out_unitp,*) 'END dealloc_ZmatTransfo'; call flush_perso(out_unitp)

      END SUBROUTINE dealloc_ZmatTransfo

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE Read_ZmatTransfo(ZmatTransfo,mendeleev)


       TYPE (Type_ZmatTransfo),intent(inout) :: ZmatTransfo
       TYPE (table_atom), intent(in)         :: mendeleev

      integer                  :: n1,n2,n3
      real (kind=Rkind)        :: at
      character (len=Name_len) :: name_at

      real (kind=Rkind)        :: masses(3*ZmatTransfo%nat)
      integer                  :: Z(ZmatTransfo%nat)
      character (len=Name_len) :: symbole(ZmatTransfo%nat)


        integer                :: ic1,ic2,ic3,icf
        integer                :: nat_dum
        integer                :: i,j
        integer                :: ZZ,iz,it
        real (kind=Rkind)      :: d1

       integer :: err_mem,memory,err_io
       logical, parameter :: debug=.FALSE.
       !logical, parameter :: debug=.TRUE.
       character (len=*), parameter :: name_sub = 'Read_ZmatTransfo'


!-----------------------------------------------------------------------
       IF (print_level > 1) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'nat0,nat',ZmatTransfo%nat0,ZmatTransfo%nat
         write(out_unitp,*) 'nb_var',ZmatTransfo%nb_var
         ZmatTransfo%ncart = 3 * ZmatTransfo%nat
         write(out_unitp,*) 'ncart',ZmatTransfo%ncart
         write(out_unitp,*) 'cos_th',ZmatTransfo%cos_th
       END IF

       it = 0

       ! allocation of the variables:
       CALL alloc_ZmatTransfo(ZmatTransfo)


        ZmatTransfo%nat_act = 0
        nat_dum = ZmatTransfo%nat


        IF (ZmatTransfo%nat0 >= 1) THEN
          iz  = 0
          i   = 1
          IF (print_level > 1) write(out_unitp,*) "==================",i
          read(in_unitp,*,IOSTAT=err_io) name_at
          IF (err_io /= 0) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) '  while reading the first line ',       &
                                         'of the "Zmat" transformation.'
            write(out_unitp,*) ' Check your data !!'
            STOP
          END IF
          ZZ = -1
          at = get_mass_Tnum(mendeleev,Z=ZZ,name=name_at)
          IF (print_level > 1) write(out_unitp,*) i,ZZ,at

          ZmatTransfo%ind2_zmat(1,i) = i
          ZmatTransfo%ind2_zmat(2,i) = 0
          ZmatTransfo%ind2_zmat(3,i) = 0
          ZmatTransfo%ind2_zmat(4,i) = 0
          ZmatTransfo%ind2_zmat(5,i) = 0

          IF (at > ZERO) THEN
             ZmatTransfo%nat_act          = ZmatTransfo%nat_act + 1
             symbole(ZmatTransfo%nat_act) = name_at
             Z(ZmatTransfo%nat_act)       = ZZ
             icf                          = func_ic(ZmatTransfo%nat_act)
          ELSE
             nat_dum          = nat_dum - 1
             symbole(nat_dum) = name_at
             Z(nat_dum)       = ZZ
             icf              = func_ic(nat_dum)
          END IF
          masses(icf+0:icf+2)       = at
          ZmatTransfo%ind_zmat(1,i) = icf
          ZmatTransfo%ind_zmat(2,i) = 0
          ZmatTransfo%ind_zmat(3,i) = 0
          ZmatTransfo%ind_zmat(4,i) = 0
          ZmatTransfo%ind_zmat(5,i) = 0


          IF (ZmatTransfo%nat0 >= 2) THEN

            i   = 2
            IF (print_level > 1) write(out_unitp,*) "==================",i
            read(in_unitp,*,IOSTAT=err_io) name_at,n1
            IF (err_io /= 0) THEN
              write(out_unitp,*) ' ERROR in ',name_sub
              write(out_unitp,*) '  while reading the second line ',    &
                                         'of the "Zmat" transformation.'
              write(out_unitp,*) ' Check your data !!'
              STOP
            END IF
            ZZ = -1
            at = get_mass_Tnum(mendeleev,Z=ZZ,name=name_at)
            IF (print_level > 1) write(out_unitp,*) i,ZZ,at,n1

            iz = iz+1
            ZmatTransfo%type_Qin(iz) = 2
            CALL make_nameQ(ZmatTransfo%name_Qin(iz),"Qzmat_d",iz,it)

            ZmatTransfo%ind2_zmat(1,i) = i
            ZmatTransfo%ind2_zmat(2,i) = n1
            ZmatTransfo%ind2_zmat(3,i) = 0
            ZmatTransfo%ind2_zmat(4,i) = 0
            ZmatTransfo%ind2_zmat(5,i) = 0


           IF (n1 == 0 ) THEN
              write(out_unitp,*) 'ERROR in ',name_sub
              write(out_unitp,*) 'The second atom can NOT be in cartesian'
              STOP
            END IF

            IF (at > ZERO) THEN
              ZmatTransfo%nat_act          = ZmatTransfo%nat_act + 1
              symbole(ZmatTransfo%nat_act) = name_at
              Z(ZmatTransfo%nat_act)       = ZZ
              icf                          = func_ic(ZmatTransfo%nat_act)
            ELSE
              nat_dum          = nat_dum - 1
              symbole(nat_dum) = name_at
              Z(nat_dum)       = ZZ
              icf              = func_ic(nat_dum)
            END IF
            ic1 = ZmatTransfo%ind_zmat(1,n1)
            masses(icf+0:icf+2)       = at
            ZmatTransfo%ind_zmat(1,i) = icf
            ZmatTransfo%ind_zmat(2,i) = ic1
            ZmatTransfo%ind_zmat(3,i) = 0
            ZmatTransfo%ind_zmat(4,i) = 0
            ZmatTransfo%ind_zmat(5,i) = 0

            IF (ZmatTransfo%nat0 >= 3) THEN

              i   = 3
              IF (print_level > 1) write(out_unitp,*) "==================",i
              read(in_unitp,*,IOSTAT=err_io) name_at,n1,n2
              IF (err_io /= 0) THEN
                write(out_unitp,*) ' ERROR in ',name_sub
                write(out_unitp,*) '  while reading the third line ',   &
                                         'of the "Zmat" transformation.'
                write(out_unitp,*) ' Check your data !!'
                STOP
              END IF

              ZZ = -1
              at = get_mass_Tnum(mendeleev,Z=ZZ,name=name_at)
              IF (print_level > 1) write(out_unitp,*) i,ZZ,at,n1,n2

              ZmatTransfo%ind2_zmat(1,i) = i
              ZmatTransfo%ind2_zmat(2,i) = n1
              ZmatTransfo%ind2_zmat(3,i) = n2
              ZmatTransfo%ind2_zmat(4,i) = 0
              ZmatTransfo%ind2_zmat(5,i) = 0

              IF (n1 == 0) THEN
                write(out_unitp,*) 'ERROR in ',name_sub
                write(out_unitp,*) 'The third atom can NOT be in cartesian'
                STOP
              END IF

              iz = iz+1
              ZmatTransfo%type_Qin(iz) = 2
              CALL make_nameQ(ZmatTransfo%name_Qin(iz),"Qzmat_d",iz,it)

              IF (at > ZERO) THEN
                ZmatTransfo%nat_act          = ZmatTransfo%nat_act + 1
                symbole(ZmatTransfo%nat_act) = name_at
                Z(ZmatTransfo%nat_act)       = ZZ
                icf                          = func_ic(ZmatTransfo%nat_act)
              ELSE
                nat_dum          = nat_dum - 1
                symbole(nat_dum) = name_at
                Z(nat_dum)       = ZZ
                icf              = func_ic(nat_dum)
              END IF
              ic1 = ZmatTransfo%ind_zmat(1,n1)
              iz = iz+1
              IF (n2 == 0) THEN
                ic2 = 0
                IF (ZmatTransfo%cos_th) THEN
                  ZmatTransfo%type_Qin(iz) = -3 ! cos(angle)
                  CALL make_nameQ(ZmatTransfo%name_Qin(iz),"Qzmat_Costh",iz,it)
                  IF (print_level > 1) write(out_unitp,*) at,n1,'polyspherical with cos(th)'
                ELSE
                  ZmatTransfo%type_Qin(iz) = 3  ! angle
                  CALL make_nameQ(ZmatTransfo%name_Qin(iz),"Qzmat_th",iz,it)
                  IF (print_level > 1) write(out_unitp,*) at,n1,'polyspherical with th'
                END IF
              ELSE IF (n2 > 0) THEN
                ic2 = ZmatTransfo%ind_zmat(1,n2)
                ZmatTransfo%type_Qin(iz) = 3 ! valence angle
                CALL make_nameQ(ZmatTransfo%name_Qin(iz),"Qzmat_th",iz,it)
              ELSE
                ic2 = ZmatTransfo%ind_zmat(1,-n2)
                ZmatTransfo%type_Qin(iz) = -3 ! cos(angle)
                CALL make_nameQ(ZmatTransfo%name_Qin(iz),"Qzmat_Costh",iz,it)
              END IF
              masses(icf+0:icf+2)       = at
              ZmatTransfo%ind_zmat(1,i) = icf
              ZmatTransfo%ind_zmat(2,i) = ic1
              ZmatTransfo%ind_zmat(3,i) = ic2
              ZmatTransfo%ind_zmat(4,i) = 0
              ZmatTransfo%ind_zmat(5,i) = 0

              DO i=4,ZmatTransfo%nat0

                IF (print_level > 1) write(out_unitp,*) "==================",i
                read(in_unitp,*,IOSTAT=err_io) name_at,n1,n2,n3
                IF (err_io /= 0) THEN
                  write(out_unitp,*) ' ERROR in ',name_sub
                  write(out_unitp,'(a,i0,a)') '  while reading the ',i, &
                                 'th line of the "Zmat" transformation.'
                  write(out_unitp,*) ' Check your data !!'
                  STOP
                END IF
                ZZ = -1
                at = get_mass_Tnum(mendeleev,Z=ZZ,name=name_at)
                IF (print_level > 1) write(out_unitp,*) i,ZZ,at,n1,n2,n3
                ZmatTransfo%ind2_zmat(1,i)=i
                ZmatTransfo%ind2_zmat(2,i)=n1
                ZmatTransfo%ind2_zmat(3,i)=n2
                ZmatTransfo%ind2_zmat(4,i)=n3
                ZmatTransfo%ind2_zmat(5,i)=0

                IF (n1 == 0) THEN
!                 l'atome est defini en coordonnees cartesiennes
                  IF (print_level > 1) write(out_unitp,*) at,'cart'
                  IF (at > ZERO) THEN
                   ZmatTransfo%nat_act          = ZmatTransfo%nat_act + 1
                   symbole(ZmatTransfo%nat_act) = name_at
                   Z(ZmatTransfo%nat_act)       = ZZ
                   icf                          = func_ic(ZmatTransfo%nat_act)
                  ELSE
                   nat_dum          = nat_dum - 1
                   symbole(nat_dum) = name_at
                   Z(nat_dum)       = ZZ
                   icf              = func_ic(nat_dum)
                  END IF
                  ic1 = 0
                  ic2 = 0
                  ic3 = 0
                  iz = iz+1
                  ZmatTransfo%type_Qin(iz) = 1  ! cartesian
                  CALL make_nameQ(ZmatTransfo%name_Qin(iz),"Qzmat_x",iz,it)
                  iz = iz+1
                  ZmatTransfo%type_Qin(iz) = 1 ! cartesian
                  CALL make_nameQ(ZmatTransfo%name_Qin(iz),"Qzmat_y",iz,it)
                  iz = iz+1
                  ZmatTransfo%type_Qin(iz) = 1 ! cartesian
                  CALL make_nameQ(ZmatTransfo%name_Qin(iz),"Qzmat_z",iz,it)

                ELSE
!                 at en coord internes
                  IF (print_level > 1) write(out_unitp,*) at,n1,n2,n3
                  IF (at > ZERO) THEN
                    ZmatTransfo%nat_act          = ZmatTransfo%nat_act + 1
                    symbole(ZmatTransfo%nat_act) = name_at
                    Z(ZmatTransfo%nat_act)       = ZZ
                    icf                          = func_ic(ZmatTransfo%nat_act)
                  ELSE
                    nat_dum          = nat_dum - 1
                    symbole(nat_dum) = name_at
                    Z(nat_dum)       = ZZ
                    icf              = func_ic(nat_dum)
                  END IF

                  iz = iz+1
                  ZmatTransfo%type_Qin(iz) = 2 ! distance
                  CALL make_nameQ(ZmatTransfo%name_Qin(iz),"Qzmat_d",iz,it)

                  ic1 = ZmatTransfo%ind_zmat(1,n1)
                  iz = iz+1
                  IF (n2 == 0) THEN
                    ic2 = 0
                    ic3 = 0
                    IF (ZmatTransfo%cos_th) THEN
                      ZmatTransfo%type_Qin(iz) = -3 ! cos(angle)
                      CALL make_nameQ(ZmatTransfo%name_Qin(iz),"Qzmat_Costh",iz,it)
                      IF (print_level > 1) write(out_unitp,*) at,n1,'polyspherical with cos(th)'
                    ELSE
                      ZmatTransfo%type_Qin(iz) = 3 ! angle
                      CALL make_nameQ(ZmatTransfo%name_Qin(iz),"Qzmat_th",iz,it)
                      IF (print_level > 1) write(out_unitp,*) at,n1,'polyspherical with th'
                    END IF
                  ELSE IF (n2 > 0) THEN
                    ic2 = ZmatTransfo%ind_zmat(1,n2)
                    ic3 = ZmatTransfo%ind_zmat(1,n3)
                    ZmatTransfo%type_Qin(iz) = 3 ! valence angle
                    CALL make_nameQ(ZmatTransfo%name_Qin(iz),"Qzmat_th",iz,it)
                  ELSE
                    ic2 = ZmatTransfo%ind_zmat(1,-n2)
                    ic3 = ZmatTransfo%ind_zmat(1,n3)
                    ZmatTransfo%type_Qin(iz) = -3 ! cos(angle)
                    CALL make_nameQ(ZmatTransfo%name_Qin(iz),"Qzmat_Costh",iz,it)
                  END IF

                  iz = iz+1
                  ZmatTransfo%type_Qin(iz) = 4 ! diedral angle
                  CALL make_nameQ(ZmatTransfo%name_Qin(iz),"Qzmat_phi",iz,it)
                ENDIF
                masses(icf+0:icf+2)       = at
                ZmatTransfo%ind_zmat(1,i) = icf
                ZmatTransfo%ind_zmat(2,i) = ic1
                ZmatTransfo%ind_zmat(3,i) = ic2
                ZmatTransfo%ind_zmat(4,i) = ic3
                ZmatTransfo%ind_zmat(5,i) = 0
              END DO
            END IF
          END IF
        ELSE
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' There is no atoms !!'
          STOP
        END IF

!       ncart_act number of active cartesian coordinates (without dummy atom and G)
        ZmatTransfo%ncart_act = 3 * ZmatTransfo%nat_act

!--------------------------------------------------------------
!     -------------------------------------------------------
!     Mtot_inv and sqrt(masses(i)) calculations
!     -------------------------------------------------------

!     -------------------------------------------------------
!--------------------------------------------------------------
       ZmatTransfo%Z(:)       = Z(:)
       ZmatTransfo%symbole(:) = symbole(:)
       ZmatTransfo%masses(:)  = masses(:)

      IF (print_level > 1) write(out_unitp,*) 'END ',name_sub
      END SUBROUTINE Read_ZmatTransfo


!=================================================================
!
!       ZmatTransfo
!
!=================================================================
      SUBROUTINE calc_ZmatTransfo(dnQzmat,dnx,ZmatTransfo,nderiv)
      IMPLICIT NONE

      TYPE (Type_dnVec), intent(inout)   :: dnQzmat,dnx
      TYPE (Type_ZmatTransfo),intent(in) :: ZmatTransfo
      integer, intent(in)                :: nderiv

      TYPE (Type_dnS)    :: dnd,dnQval,dnCval,dnSval,dnQdih,dnCdih,dnSdih

      TYPE (Type_dnS)    :: dnf1,dnf2,dnf3
      TYPE (Type_dnVec)  :: dnv1,dnv2,dnv3

       TYPE (Type_dnVec) :: dnEz2,dnEz3,dnEy3,dnEx3,dnAt1
       real (kind=Rkind) :: d1,s12,nEx3,nEy3,nEz3
!      ---------------------------------------------


       logical :: case1

       integer :: ic,ic1,ic2,ic3,icf,icG
       integer :: i_q
       integer :: i
       integer :: nb_act

       logical :: check

!      -----------------------------------------------------------------
      integer :: nderiv_debug = 0
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='calc_ZmatTransfo'
!      -----------------------------------------------------------------
       IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nderiv',nderiv
        write(out_unitp,*)
        CALL Write_ZmatTransfo(ZmatTransfo)
        write(out_unitp,*) 'dnQzmat'
        CALL Write_dnSVM(dnQzmat,nderiv)
       END IF
      !-----------------------------------------------------------------
      nb_act = dnQzmat%nb_var_deriv

        CALL alloc_dnSVM(dnAt1,3,nb_act,nderiv)
        CALL alloc_dnSVM(dnEz2,3,nb_act,nderiv)
        CALL alloc_dnSVM(dnEz3,3,nb_act,nderiv)
        CALL alloc_dnSVM(dnEx3,3,nb_act,nderiv)
        CALL alloc_dnSVM(dnEy3,3,nb_act,nderiv)

        CALL alloc_dnSVM(dnd,nb_act,nderiv)
        CALL alloc_dnSVM(dnQval,nb_act,nderiv)
        CALL alloc_dnSVM(dnCval,nb_act,nderiv)
        CALL alloc_dnSVM(dnSval,nb_act,nderiv)
        CALL alloc_dnSVM(dnQdih,nb_act,nderiv)
        CALL alloc_dnSVM(dnCdih,nb_act,nderiv)
        CALL alloc_dnSVM(dnSdih,nb_act,nderiv)

        CALL alloc_dnSVM(dnf1,nb_act,nderiv)
        CALL alloc_dnSVM(dnf2,nb_act,nderiv)
        CALL alloc_dnSVM(dnf3,nb_act,nderiv)

        CALL alloc_dnSVM(dnv1,3,nb_act,nderiv)
        CALL alloc_dnSVM(dnv2,3,nb_act,nderiv)
        CALL alloc_dnSVM(dnv3,3,nb_act,nderiv)

        !=================================================
        ! initialization: useless
        !=================================================
        CALL Set_ZERO_TO_dnSVM(dnx)

        !=================================================
        ! first atom
        !=================================================
        i   = 1
        icf = ZmatTransfo%ind_zmat(1,i)

        CALL Set_ZERO_TO_dnSVM(dnAt1)
        IF (ZmatTransfo%New_Orient) THEN
          dnAt1%d0(:) = ZmatTransfo%vAt1(:)
        END IF

        CALL sub3_dnVec_TOxf(dnx,icf,dnAt1,nderiv)
        !-----------------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*)
          write(out_unitp,*) '-------------------------------------------------'
          write(out_unitp,*) 'atom :',i
          CALL write_dnx(1,dnx%nb_var_vec,dnx,nderiv_debug)
        END IF
        !-----------------------------------------------------------------


        !-IF check=t => check distances ----------------------------------
        check = .FALSE.
        IF (ZmatTransfo%nat0 == 2) check = .FALSE.
        !check = .NOT. ZmatTransfo%nat0 == 2

        i_q = 1
        IF (ZmatTransfo%nat0 >= 2) THEN

          !=================================================
          ! 2d    atom

          CALL sub_dnVec_TO_dnS(dnQzmat,dnd,i_q)
!         CALL Write_dnSVM(dnd,nderiv)
          CALL check_Valence(i_q,dnd%d0,ZmatTransfo%type_Qin(i_q))

          i   = 2
          icf = ZmatTransfo%ind_zmat(1,i)
          ic1 = ZmatTransfo%ind_zmat(2,i)
!         write(out_unitp,*) 'icf,ic1',icf,ic1

          CALL Set_ZERO_TO_dnSVM(dnEz2)
          IF (ZmatTransfo%New_Orient) THEN
            dnEz2%d0(:) = ZmatTransfo%vAt2(:)-                 &
                          ZmatTransfo%vAt1(:)
            d1 = sqrt(dot_product(dnEz2%d0,dnEz2%d0))
            dnEz2%d0(:) = dnEz2%d0(:)/d1
          ELSE
            !--- Z2 axis along z_BF ------
            dnEz2%d0(3) = ONE
          END IF

          !write(6,*) 'dnEz2',dnEz2%d0

          CALL sub3_dnx_AT2_new(dnx,icf,ic1,dnd,dnEz2,nderiv,check)

          !-----------------------------------------------------------------
          IF (debug) THEN
            write(out_unitp,*)
            write(out_unitp,*) '-------------------------------------------------'
            write(out_unitp,*) 'atom :',i
            CALL write_dnx(1,dnx%nb_var_vec,dnx,nderiv_debug)
          END IF
          !-----------------------------------------------------------------
          !=================================================

          IF (ZmatTransfo%nat0 >= 3) THEN

            !=================================================
            !3d    atom

            i   = 3

            i_q = i_q + 1
            CALL sub_dnVec_TO_dnS(dnQzmat,dnd,i_q)
!           CALL Write_dnSVM(dnd,nderiv)
            CALL check_Valence(i_q,dnd%d0,ZmatTransfo%type_Qin(i_q))

            i_q = i_q + 1
            CALL sub_dnVec_TO_dnS(dnQzmat,dnQval,i_q)
!           CALL Write_dnSVM(dnQval,nderiv)

            IF (ZmatTransfo%type_Qin(i_q) == 3) THEN
              CALL sub_dnS1_TO_dntR2(dnQval,dnCval,2,nderiv)
              CALL sub_dnS1_TO_dntR2(dnQval,dnSval,3,nderiv)
            ELSE ! type_Qin(i_q) == -3 (using Q=cos(val) as coordinate)
              CALL sub_dnS1_TO_dnS2(dnQval,dnCval,nderiv)
              CALL sub_dnS1_TO_dntR2(dnQval,dnSval,4,nderiv)
            END IF
!           CALL Write_dnSVM(dnCval,nderiv)
!           CALL Write_dnSVM(dnSval,nderiv)

            icf = ZmatTransfo%ind_zmat(1,i)
            ic1 = ZmatTransfo%ind_zmat(2,i)
            ic2 = ZmatTransfo%ind_zmat(3,i)

            CALL Set_ZERO_TO_dnSVM(dnEx3)
            CALL Set_ZERO_TO_dnSVM(dnEz3)
            IF (ic2 == 0) THEN    ! polyspherical
              dnEx3%d0(1) = ONE
              dnEz3%d0(3) = ONE
            ELSE                  ! true zmat

              case1 = (ZmatTransfo%ind_zmat(2,i) ==            &
                       ZmatTransfo%ind_zmat(1,1) )
!             write(out_unitp,*) 'icf,ic1,ic2,case1',icf,ic1,ic2,case1
              CALL check_Valence(i_q,dnQval%d0,ZmatTransfo%type_Qin(i_q))

              IF (ZmatTransfo%New_Orient) THEN
                IF (case1) THEN
                  dnEz3%d0(:) = ZmatTransfo%vAt2(:)-           &
                                ZmatTransfo%vAt1(:)
                  dnEx3%d0(:) = ZmatTransfo%vAt3(:)-           &
                                ZmatTransfo%vAt1(:)
                ELSE
                  dnEz3%d0(:) = ZmatTransfo%vAt1(:)-           &
                                ZmatTransfo%vAt2(:)
                  dnEx3%d0(:) = ZmatTransfo%vAt3(:)-           &
                                ZmatTransfo%vAt2(:)
                END IF
                d1 = sqrt(dot_product(dnEz3%d0,dnEz3%d0))
                dnEz3%d0(:) = dnEz3%d0(:)/d1
                s12 = dot_product(dnEz3%d0,dnEx3%d0)

                dnEx3%d0(:) = dnEx3%d0(:) - dnEz3%d0(:) * s12
                dnEx3%d0(:) = dnEx3%d0(:) /                             &
                    sqrt(dot_product(dnEx3%d0,dnEx3%d0))
              ELSE
                !--- Z3 axis along z_BF and x3 along x_BF ------
                IF (case1) THEN
                  dnEx3%d0(1) = ONE
                  dnEz3%d0(3) = ONE
                ELSE
                  dnEx3%d0(1) = ONE
                  dnEz3%d0(3) = -ONE
                END IF
              END IF
            END IF
            !write(6,*) 'New_Orient',ZmatTransfo%New_Orient
            !write(6,*) 'dnEx3',dnEx3%d0
            !write(6,*) 'dnEz3',dnEz3%d0

            CALL sub3_dnx_AT3_new(dnx,icf,ic1,check,                    &
                                  dnd,dnCval,dnSval,                    &
                                  dnEz3,dnEx3,dnf1,nderiv)

            !-----------------------------------------------------------------
            IF (debug) THEN
              write(out_unitp,*)
              write(out_unitp,*) '-------------------------------------------------'
              write(out_unitp,*) 'atom :',i
              CALL write_dnx(1,dnx%nb_var_vec,dnx,nderiv_debug)
            END IF
            !-----------------------------------------------------------------
            !=================================================
            !we used nat0(=nat-1), because the atom "nat" is the center of masse
            DO i=4,ZmatTransfo%nat0

              !=================================================
              ! 4th ... atom

              icf = ZmatTransfo%ind_zmat(1,i)
              ic1 = ZmatTransfo%ind_zmat(2,i)
              ic2 = ZmatTransfo%ind_zmat(3,i)
              ic3 = ZmatTransfo%ind_zmat(4,i)

              IF (ic1 == 0) THEN !  atome en cartesiennes

                IF (ZmatTransfo%New_Orient) THEN
                  IF (case1) THEN
                    dnEz3%d0(:) = ZmatTransfo%vAt2(:)-           &
                                  ZmatTransfo%vAt1(:)
                    dnEx3%d0(:) = ZmatTransfo%vAt3(:)-           &
                                  ZmatTransfo%vAt1(:)
                  ELSE
                    dnEz3%d0(:) = ZmatTransfo%vAt1(:)-           &
                                  ZmatTransfo%vAt2(:)
                    dnEx3%d0(:) = ZmatTransfo%vAt3(:)-           &
                                  ZmatTransfo%vAt2(:)
                  END IF
                  d1 = sqrt(dot_product(dnEz3%d0,dnEz3%d0))
                  dnEz3%d0(:) = dnEz3%d0(:)/d1
                  s12 = dot_product(dnEz3%d0,dnEx3%d0)

                  dnEx3%d0(:) = dnEx3%d0(:) - dnEz3%d0(:) * s12
                  dnEx3%d0(:) = dnEx3%d0(:) /                           &
                      sqrt(dot_product(dnEx3%d0,dnEx3%d0))
                ELSE
                  !--- Z3 axis along z_BF and x3 along x_BF ------
                  IF (case1) THEN
                    dnEx3%d0(1) = ONE
                    dnEz3%d0(3) = ONE
                  ELSE
                    dnEx3%d0(1) = ONE
                    dnEz3%d0(3) = -ONE
                  END IF
                END IF
                CALL calc_cross_product(dnEz3%d0,nEz3,dnEx3%d0,nEx3,dnEy3%d0,nEy3)
                !write(6,*) 'dnEx3',dnEx3%d0
                !write(6,*) 'dnEy3',dnEy3%d0
                !write(6,*) 'dnEz3',dnEz3%d0

                i_q = i_q + 1
                CALL sub_dnVec_TO_dnS(dnQzmat,dnd,i_q)

                i_q = i_q + 1
                CALL sub_dnVec_TO_dnS(dnQzmat,dnQval,i_q)

                i_q = i_q + 1
                CALL sub_dnVec_TO_dnS(dnQzmat,dnQdih,i_q)


                CALL sub3_dnx_AT4_cart_new(dnx,icf,dnd,dnQval,dnQdih,dnEx3,dnEy3,dnEz3,nderiv)

                !CALL sub3_dnx_AT4_cart(dnx,icf,dnd,dnQval,dnQdih,nderiv)
                IF (ZmatTransfo%New_Orient) THEN
                  icf = ZmatTransfo%ind_zmat(1,i)
                  dnx%d0(icf:icf+2) = dnx%d0(icf:icf+2) + ZmatTransfo%vAt1(:)
                END IF

              ELSE IF (ic2 == 0) THEN    ! polyspherical

                i_q = i_q + 1
                CALL sub_dnVec_TO_dnS(dnQzmat,dnd,i_q)
                CALL check_Valence(i_q,dnd%d0,ZmatTransfo%type_Qin(i_q))

                i_q = i_q + 1
                CALL sub_dnVec_TO_dnS(dnQzmat,dnQval,i_q)
                IF (ZmatTransfo%type_Qin(i_q) == 3) THEN
                  CALL sub_dnS1_TO_dntR2(dnQval,dnCval,2,nderiv)
                  CALL sub_dnS1_TO_dntR2(dnQval,dnSval,3,nderiv)
                ELSE ! type_Qin(i_q) == -3 (using Q=cos(val) as coordinate)
                  CALL sub_dnS1_TO_dnS2(dnQval,dnCval,nderiv)
                  CALL sub_dnS1_TO_dntR2(dnQval,dnSval,4,nderiv)
                END IF

                IF (ZmatTransfo%ind_zmat(2,i) /= 0) THEN
                  CALL check_Valence(i_q,dnQval%d0,ZmatTransfo%type_Qin(i_q))
                END IF

                i_q = i_q + 1
                CALL sub_dnVec_TO_dnS(dnQzmat,dnQdih,i_q)
                CALL sub_dnS1_TO_dntR2(dnQdih,dnCdih,2,nderiv)
                CALL sub_dnS1_TO_dntR2(dnQdih,dnSdih,3,nderiv)

                CALL sub3_dnx_AT4_poly(dnx,icf,ic1,ic2,ic3,check,       &
                                       dnd,dncval,dnsval,dncdih,dnsdih, &
                                       dnv1,dnf1,dnf2,dnf3,nderiv)
              ELSE                  ! true zmat

                i_q = i_q + 1
                CALL sub_dnVec_TO_dnS(dnQzmat,dnd,i_q)
                CALL check_Valence(i_q,dnd%d0,ZmatTransfo%type_Qin(i_q))

                i_q = i_q + 1
                CALL sub_dnVec_TO_dnS(dnQzmat,dnQval,i_q)
                IF (ZmatTransfo%type_Qin(i_q) == 3) THEN
                  CALL sub_dnS1_TO_dntR2(dnQval,dnCval,2,nderiv)
                  CALL sub_dnS1_TO_dntR2(dnQval,dnSval,3,nderiv)
                ELSE ! type_Qin(i_q) == -3 (using Q=cos(val) as coordinate)
                  CALL sub_dnS1_TO_dnS2(dnQval,dnCval,nderiv)
                  CALL sub_dnS1_TO_dntR2(dnQval,dnSval,4,nderiv)
                END IF

                IF (ZmatTransfo%ind_zmat(2,i) /= 0) THEN
                  CALL check_Valence(i_q,dnQval%d0,ZmatTransfo%type_Qin(i_q))
                END IF

                i_q = i_q + 1
                CALL sub_dnVec_TO_dnS(dnQzmat,dnQdih,i_q)
                CALL sub_dnS1_TO_dntR2(dnQdih,dnCdih,2,nderiv)
                CALL sub_dnS1_TO_dntR2(dnQdih,dnSdih,3,nderiv)

                CALL sub3_dnx_AT4(dnx,icf,ic1,ic2,ic3,check,            &
                                  dnd,dncval,dnsval,dncdih,dnsdih,      &
                                  dnv1,dnv2,dnv3,dnf1,dnf2,dnf3,nderiv)
              END IF

              !-----------------------------------------------------------------
              IF (debug) THEN
                write(out_unitp,*)
                write(out_unitp,*) '-------------------------------------------------'
                write(out_unitp,*) 'atom :',i
                CALL write_dnx(1,dnx%nb_var_vec,dnx,nderiv_debug)
              END IF
              !-----------------------------------------------------------------
              !=================================================
            END DO
          END IF
        ELSE
          write(out_unitp,*) ' STOP in ',name_sub
          write(out_unitp,*) ' ERROR : there is no atoms !!'
          STOP
        END IF
        !=================================================


        !=================================================
        !-----------------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'Final Cartesian coordinates:'
          CALL write_dnx(1,dnx%nb_var_vec,dnx,nderiv_debug)
          write(out_unitp,*) 'END ',name_sub
          write(out_unitp,*)
        END IF
        !-----------------------------------------------------------------
        !=================================================

        CALL dealloc_dnSVM(dnd)
        CALL dealloc_dnSVM(dnQval)
        CALL dealloc_dnSVM(dnCval)
        CALL dealloc_dnSVM(dnSval)
        CALL dealloc_dnSVM(dnQdih)
        CALL dealloc_dnSVM(dnCdih)
        CALL dealloc_dnSVM(dnSdih)

        CALL dealloc_dnSVM(dnf1)
        CALL dealloc_dnSVM(dnf2)
        CALL dealloc_dnSVM(dnf3)

        CALL dealloc_dnSVM(dnEz2)
        CALL dealloc_dnSVM(dnEz3)
        CALL dealloc_dnSVM(dnEx3)
        CALL dealloc_dnSVM(dnAt1)

        CALL dealloc_dnSVM(dnv1)
        CALL dealloc_dnSVM(dnv2)
        CALL dealloc_dnSVM(dnv3)

      END SUBROUTINE calc_ZmatTransfo


      SUBROUTINE calc_ZmatTransfo_outTOin(dnQzmat,dnx,ZmatTransfo,nderiv)
      IMPLICIT NONE

      TYPE (Type_dnVec), intent(inout) :: dnQzmat,dnx
      TYPE (Type_ZmatTransfo),intent(in) :: ZmatTransfo
      integer, intent(in)              :: nderiv

       integer           :: ncart0
       real (kind=Rkind) :: v1(3),norm1,v2(3),norm2,v3(3),norm3
       real (kind=Rkind) :: v4(3),norm4,v5(3),norm5
       real (kind=Rkind) :: angle_v,angle_d
       integer           :: nc0,nc1,nc2,nc3,nc4,idum,iqz,i

       real (kind=Rkind) :: ex(3),nx,ey(3),ny,ez(3),nz


!      -----------------------------------------------------------------
      integer :: nderiv_debug = 0
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='calc_ZmatTransfo_outTOin'
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nderiv',nderiv
        write(out_unitp,*)
        CALL Write_ZmatTransfo(ZmatTransfo)

        write(out_unitp,*) 'Cartesian coordinates:'
        CALL write_dnx(1,dnx%nb_var_vec,dnx,nderiv_debug)
      END IF
      !-----------------------------------------------------------------
      ncart0 = 3*ZmatTransfo%nat0

      IF (nderiv > 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' nderiv > 0 is not possible'
        write(out_unitp,*) ' It should NOT append ! Check the Fortran source'
        STOP
      END IF

      iqz = 0
      DO i=2,ZmatTransfo%nat0
        nc1 = ZmatTransfo%ind_zmat(1,i)
        nc2 = ZmatTransfo%ind_zmat(2,i)
        nc3 = abs(ZmatTransfo%ind_zmat(3,i))
        nc4 = ZmatTransfo%ind_zmat(4,i)
        IF (debug) write(out_unitp,*) '-------------------',i
        IF (debug) write(out_unitp,*) 'nc1,nc2,nc3,nc4',nc1,nc2,nc3,nc4


        IF (i == 2) THEN

          IF (nc2==0) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' Your zmatrix are using:'
            write(out_unitp,*) ' cartesian coordinates for the 2d atom'
            write(out_unitp,*) 'zmat at:',i,ZmatTransfo%ind2_zmat(:,i)
            STOP
          END IF

          CALL calc_vector2(v1,norm1,nc2,nc1,dnx%d0,ncart0)
          iqz = iqz + 1
          dnQzmat%d0(iqz) = norm1
          ez(:) = v1(:) / norm1
          IF (debug) write(out_unitp,*) ' nc1,nc2,v1',nc1,nc2,v1(:)
          IF (debug) write(out_unitp,*) ' nc1,nc2,d',nc1,nc2,norm1

        ELSE IF (i == 3) THEN

          IF (nc2==0) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' Your zmatrix are using:'
            write(out_unitp,*) ' cartesian coordinates for the 3d atom'
            write(out_unitp,*) 'zmat at:',i,ZmatTransfo%ind2_zmat(:,i)
            STOP
          END IF

          CALL calc_vector2(v1,norm1,nc2,nc1,dnx%d0,ncart0)
          CALL calc_vector2(v2,norm2,nc2,nc3,dnx%d0,ncart0)


          iqz = iqz + 1
          dnQzmat%d0(iqz) = norm1

          IF (debug) write(out_unitp,*) nc1,nc2,v1,norm1
          IF (debug) write(out_unitp,*) nc3,nc2,v2,norm2

          CALL calc_angle(angle_v,v1,norm1,v2,norm2)
          IF ( abs(sin(angle_v)) < ONETENTH**5 ) THEN
            write(out_unitp,*) ' WARNNING in ',name_sub
            write(out_unitp,*) ' 3 atoms are aligned!',nc1,nc2,nc3
            write(out_unitp,*) ' I cannot calculate the valence angle'
            write(out_unitp,*) ' angle,cos(angle)',angle_v,cos(angle_v)
            write(out_unitp,*) ' Check your data !'
            DO idum=1,ZmatTransfo%nat0
              write(out_unitp,*) idum,dnx%d0(3*idum-2:3*idum)
            END DO
          END IF
          iqz = iqz + 1
          IF (ZmatTransfo%ind2_zmat(3,i) < 0) THEN
              dnQzmat%d0(iqz) = cos(angle_v)
          ELSE
              dnQzmat%d0(iqz) = angle_v
          END IF

          ex(:) = v1(:)-ez(:)*dot_product(v1,ez)
          ex(:) = ex(:)/sqrt(dot_product(ex,ex))
          CALL calc_cross_product(ez,nz,ex,nx,ey,ny)

          IF (debug) write(out_unitp,*) ' ex(:)',ex(:)
          IF (debug) write(out_unitp,*) ' ey(:)',ey(:)
          IF (debug) write(out_unitp,*) ' ez(:)',ez(:)


          IF (debug) write(out_unitp,*) ' nc1,nc2,nc3,d,angle',         &
                                       nc1,nc2,nc3,norm1,dnQzmat%d0(iqz)

        ELSE ! i>3

          IF (nc2==0 .AND. nc3==0 .AND. nc4==0) THEN
            nc0 = ZmatTransfo%ind_zmat(1,1) ! first atom (can be dummy)
            CALL calc_vector2(v1,norm1,nc0,nc1,dnx%d0,ncart0)
            IF (debug) write(out_unitp,*) ' nc1,nc0,v1,norm1',nc1,nc0,v1,norm1

            iqz = iqz + 1
            dnQzmat%d0(iqz) = dot_product(v1,ex) ! x
            iqz = iqz + 1
            dnQzmat%d0(iqz) = dot_product(v1,ey) ! y
            iqz = iqz + 1
            dnQzmat%d0(iqz) = dot_product(v1,ez) ! z
            IF (debug) write(out_unitp,*) ' nc1,nc2,nc3,x,y,z',         &
                                        nc1,nc2,nc3,dnQzmat%d0(iqz-2:iqz),sqrt(dnQzmat%d0(iqz-2:iqz)**2)
          ELSE IF (nc2/=0 .AND. nc3==0 .AND. nc4==0) THEN
            CALL calc_vector2(v1,norm1,nc1,nc2,dnx%d0,ncart0)
            iqz = iqz + 1
            dnQzmat%d0(iqz) = norm1

            CALL calc_angle(angle_v,ez,norm2,v1,norm1)
            IF ( abs(sin(angle_v)) < ONETENTH**5 ) THEN
              write(out_unitp,*) ' WARNNING in ',name_sub
              write(out_unitp,*) ' 3 atoms are aligned!',nc1,nc2,nc3
              write(out_unitp,*) ' I cannot calculate the valence angle'
              write(out_unitp,*) ' angle,cos(angle)',angle_v,cos(angle_v)
              write(out_unitp,*) ' Check your data !'
              DO idum=1,ZmatTransfo%nat0
                write(out_unitp,*) idum,dnx%d0(3*idum-2:3*idum)
              END DO
            END IF
            iqz = iqz + 1
            IF (ZmatTransfo%cos_th) THEN
              dnQzmat%d0(iqz) = cos(angle_v)
            ELSE
              dnQzmat%d0(iqz) = angle_v
            END IF

            angle_d = atan2(dot_product(v1,ey),dot_product(v1,ex))
            CALL dihedral_range(angle_d,2) ! [0:2pi]

            iqz = iqz + 1
            dnQzmat%d0(iqz) = angle_d
            IF (debug) write(out_unitp,*) ' nc1,nc2,nc3,R,th,phi',      &
                                        nc1,nc2,nc3,dnQzmat%d0(iqz-2:iqz)
          ELSE

            CALL calc_vector2(v1,norm1,nc2,nc1,dnx%d0,ncart0)

            CALL calc_vector2(v2,norm2,nc2,nc3,dnx%d0,ncart0)
            CALL calc_vector2(v3,norm3,nc3,nc4,dnx%d0,ncart0)

            iqz = iqz + 1
            dnQzmat%d0(iqz) = norm1

            CALL calc_angle(angle_v,v2,norm2,v1,norm1)
            IF ( abs(sin(angle_v)) < ONETENTH**5 ) THEN
              write(out_unitp,*) ' WARNNING in ',name_sub
              write(out_unitp,*) ' 3 atoms are aligned!',nc1,nc2,nc3
              write(out_unitp,*) ' I cannot calculate the valence angle'
              write(out_unitp,*) ' angle,cos(angle)',angle_v,cos(angle_v)
              write(out_unitp,*) ' Check your data !'
              DO idum=1,ZmatTransfo%nat0
                write(out_unitp,*) idum,dnx%d0(3*idum-2:3*idum)
              END DO
            END IF
            iqz = iqz + 1
            IF (ZmatTransfo%ind2_zmat(3,i) < 0) THEN
              dnQzmat%d0(iqz) = cos(angle_v)
            ELSE
              dnQzmat%d0(iqz) = angle_v
            END IF

            CALL calc_cross_product(v1,norm1,v2,norm2,v4,norm4)
            CALL calc_cross_product(v3,norm3,v2,norm2,v5,norm5)
            CALL calc_angle_d(angle_d,v4,norm4,v5,norm5,v2,norm2)
            CALL dihedral_range(angle_d,2) ! [0:2pi]

            iqz = iqz + 1
            dnQzmat%d0(iqz) = angle_d

            IF (debug) write(out_unitp,*) 'nc1,nc2,nc3,nc4,d,angle,angle_d',&
                           nc1,nc2,nc3,nc4,dnQzmat%d0(iqz-2:iqz)
          END IF

        END IF
      END DO



        !=================================================
        !-----------------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'dnQzmat'
          CALL Write_dnSVM(dnQzmat,nderiv)
          write(out_unitp,*) 'END ',name_sub
          write(out_unitp,*)
        END IF
        !-----------------------------------------------------------------
        !=================================================

      END SUBROUTINE calc_ZmatTransfo_outTOin

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE ZmatTransfo1TOZmatTransfo2(ZmatTransfo1,ZmatTransfo2)

!      for the zmatrix and Tnum --------------------------------------
      TYPE (Type_ZmatTransfo), intent(in)    :: ZmatTransfo1
      TYPE (Type_ZmatTransfo), intent(inout) :: ZmatTransfo2

      integer :: it
      character (len=*), parameter ::                                   &
                                name_sub = 'ZmatTransfo1TOZmatTransfo2'

      CALL dealloc_ZmatTransfo(ZmatTransfo2)

      ZmatTransfo2%ncart        = ZmatTransfo1%ncart
      ZmatTransfo2%ncart_act    = ZmatTransfo1%ncart_act
      ZmatTransfo2%nat          = ZmatTransfo1%nat
      ZmatTransfo2%nat0         = ZmatTransfo1%nat0
      ZmatTransfo2%nat_act      = ZmatTransfo1%nat_act
      ZmatTransfo2%nb_var       = ZmatTransfo1%nb_var

      IF (.NOT. associated(ZmatTransfo1%ind2_zmat) .OR.                 &
          .NOT. associated(ZmatTransfo1%ind_zmat) ) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' ZmatTransfo1 is NOT allocated !!'
        write(out_unitp,*) ' Check the source'
        STOP
      END IF

      CALL alloc_ZmatTransfo(ZmatTransfo2)

      ZmatTransfo2%ind2_zmat       = ZmatTransfo1%ind2_zmat
      ZmatTransfo2%ind_zmat        = ZmatTransfo1%ind_zmat

      ZmatTransfo2%masses          = ZmatTransfo1%masses
      ZmatTransfo2%Z               = ZmatTransfo1%Z
      ZmatTransfo2%symbole         = ZmatTransfo1%symbole


      ZmatTransfo2%New_Orient  = ZmatTransfo1%New_Orient
      ZmatTransfo2%vAt1(:)     = ZmatTransfo1%vAt1(:)
      ZmatTransfo2%vAt2(:)     = ZmatTransfo1%vAt2(:)
      ZmatTransfo2%vAt3(:)     = ZmatTransfo1%vAt3(:)

      ZmatTransfo2%cos_th       = ZmatTransfo1%cos_th

!     write(out_unitp,*) 'END ZmatTransfo1TOZmatTransfo2'

      END SUBROUTINE ZmatTransfo1TOZmatTransfo2

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE Write_ZmatTransfo(ZmatTransfo)
      TYPE (Type_ZmatTransfo), intent(in) :: ZmatTransfo

      integer :: i
      character (len=*), parameter :: name_sub='Write_ZmatTransfo'


      write(out_unitp,*) 'BEGINNING ',name_sub

      write(out_unitp,*) 'ncart_act,ncart',                             &
                  ZmatTransfo%ncart_act,ZmatTransfo%ncart

      write(out_unitp,*) 'nat_act,nat0,nat,',                           &
                  ZmatTransfo%nat_act,ZmatTransfo%nat0,ZmatTransfo%nat

      write(out_unitp,*) 'nb_var',ZmatTransfo%nb_var

      write(out_unitp,*) 'ind2_zmat'
      IF (associated(ZmatTransfo%ind2_zmat)) THEN
        DO i=1,ZmatTransfo%nat
          write(out_unitp,*) i,ZmatTransfo%ind2_zmat(:,i)
        END DO
      END IF
      write(out_unitp,*) 'ind_zmat,cos_th',ZmatTransfo%cos_th
      IF (associated(ZmatTransfo%ind_zmat)) THEN
        DO i=1,ZmatTransfo%nat
          write(out_unitp,*) i,ZmatTransfo%ind_zmat(:,i)
        END DO
      END IF

      write(out_unitp,*) 'New_Orient',ZmatTransfo%New_Orient
      write(out_unitp,*) 'vAt1',ZmatTransfo%vAt1(:)
      write(out_unitp,*) 'vAt2',ZmatTransfo%vAt2(:)
      write(out_unitp,*) 'vAt3',ZmatTransfo%vAt3(:)

      write(out_unitp,*) 'END ',name_sub

      END SUBROUTINE Write_ZmatTransfo

      END MODULE mod_ZmatTransfo

