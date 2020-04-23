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
      MODULE mod_ThreeDTransfo
      use mod_system
      USE mod_dnSVM
      IMPLICIT NONE

      PRIVATE

      TYPE Type_ThreeDTransfo
        integer                  :: Type_3D ! 0: cart: identity; 1: polar; 2: spherical, 3: R, x,z
        character (len=Name_len) :: name_Transfo_3D = ''
        integer                  :: list_ThreeD_coord(3) = (/ 0,0,0 /)
      END TYPE Type_ThreeDTransfo

      PUBLIC :: Type_ThreeDTransfo, alloc_ThreeDTransfo, dealloc_ThreeDTransfo
      PUBLIC :: Read_ThreeDTransfo, Write_ThreeDTransfo, calc_ThreeDTransfo
      PUBLIC :: ThreeDTransfo1TOThreeDTransfo2

      CONTAINS

      SUBROUTINE alloc_ThreeDTransfo(ThreeDTransfo)
      TYPE (Type_ThreeDTransfo), pointer, intent(inout) :: ThreeDTransfo
      integer :: err_mem,memory

      IF (.NOT. associated(ThreeDTransfo)) THEN
        allocate(ThreeDTransfo,stat=err_mem)
        memory = 1
        CALL error_memo_allo(err_mem,memory,'ThreeDTransfo','alloc_ThreeDTransfo','Type_ThreeDTransfo')
      END IF

      END SUBROUTINE alloc_ThreeDTransfo

      SUBROUTINE dealloc_ThreeDTransfo(ThreeDTransfo)
      TYPE (Type_ThreeDTransfo), pointer, intent(inout) :: ThreeDTransfo
      integer :: err_mem,memory

      IF (.NOT. associated(ThreeDTransfo)) RETURN
      deallocate(ThreeDTransfo,stat=err_mem)
      memory = 1
      CALL error_memo_allo(err_mem,-memory,'ThreeDTransfo','dealloc_ThreeDTransfo','Type_ThreeDTransfo')
      nullify(ThreeDTransfo)

      END SUBROUTINE dealloc_ThreeDTransfo

      SUBROUTINE Read_ThreeDTransfo(ThreeDTransfo,nb_Qin)

      TYPE (Type_ThreeDTransfo), pointer, intent(inout) :: ThreeDTransfo
      integer, intent(in) :: nb_Qin

      integer :: nb_coord,Type_3D,list_ThreeD_coord(3)
      logical :: multiple

       NAMELIST / ThreeD / Type_3D,list_ThreeD_coord


      integer :: err_io,err_mem,memory
      character (len=*), parameter :: name_sub='Read_ThreeDTransfo'

      CALL alloc_ThreeDTransfo(ThreeDTransfo)

      Type_3D              = 0
      list_ThreeD_coord(:) = 0
      read(in_unitp,threeD,IOSTAT=err_io)
      IF (err_io /= 0) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) '  while reading "threeD" namelist'
         write(out_unitp,*) ' end of file or end of record'
         write(out_unitp,*) ' Check your data !!'
         STOP
      END IF
      ThreeDTransfo%Type_3D           = Type_3D
      ThreeDTransfo%list_ThreeD_coord = list_ThreeD_coord

      nb_coord = count(list_ThreeD_coord /= 0)

      IF (nb_coord /= 3) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' The number of coordinates is different from 3'
        write(out_unitp,*) ' Check your data !!'
        write(out_unitp,*) 'list_ThreeD_coord: ',list_ThreeD_coord(:)
        STOP
      END IF

      multiple = (count(list_ThreeD_coord == list_ThreeD_coord(1)) /= 1)
      multiple = multiple .OR. (count(list_ThreeD_coord == list_ThreeD_coord(2)) /= 1)
      multiple = multiple .OR. (count(list_ThreeD_coord == list_ThreeD_coord(3)) /= 1)
      IF (multiple) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' Some values are identical in the list_ThreeD_coord'
        write(out_unitp,*) 'list_ThreeD_coord: ',list_ThreeD_coord(:)
        write(out_unitp,*) ' Check your data !!'
        STOP
      END IF

      SELECT CASE (Type_3D)
      CASE (0) ! identity
        ThreeDTransfo%name_Transfo_3D = 'identity: x,y,z'
      CASE (1) ! polar
        ThreeDTransfo%name_Transfo_3D = 'Polar: R,theta,z'
      CASE (2)
        ThreeDTransfo%name_Transfo_3D = 'Spherical: R,theta,phi'
      CASE (3)
        ThreeDTransfo%name_Transfo_3D = 'x,y,R'
      CASE default ! ERROR: wrong transformation !
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' The type of threeD transformation is UNKNOWN: ',Type_3D
        write(out_unitp,*) ' The possible values are:'
        write(out_unitp,*) '    0: identity: x,y,z'
        write(out_unitp,*) '    1: Polar: R,theta,z'
        write(out_unitp,*) '    2: Spherical: R,theta,phi'
        write(out_unitp,*) '    3: x,y,R'
        write(out_unitp,*) ' Check your data !!'
        STOP

      END SELECT

      CALL Write_ThreeDTransfo(ThreeDTransfo)

      END SUBROUTINE Read_ThreeDTransfo

      SUBROUTINE Write_ThreeDTransfo(ThreeDTransfo)

      TYPE (Type_ThreeDTransfo), pointer, intent(in) :: ThreeDTransfo

      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Write_ThreeDTransfo'

      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) 'asso ThreeDTransfo:',associated(ThreeDTransfo)
      IF (associated(ThreeDTransfo)) THEN
        write(out_unitp,*) 'Type_3D:           ',ThreeDTransfo%Type_3D
        write(out_unitp,*) 'name_Transfo_3D:   ',ThreeDTransfo%name_Transfo_3D
        write(out_unitp,*) 'list_ThreeD_coord: ',ThreeDTransfo%list_ThreeD_coord
      END IF
      write(out_unitp,*) 'END ',name_sub

      END SUBROUTINE Write_ThreeDTransfo

      SUBROUTINE calc_ThreeDTransfo(dnQin,dnQout,ThreeDTransfo,   &
                                    nderiv,inTOout)

        TYPE (Type_dnVec), intent(inout)              :: dnQin,dnQout
        TYPE (Type_ThreeDTransfo),pointer, intent(in) :: ThreeDTransfo
        integer, intent(in)                           :: nderiv
        logical, intent(in)                           :: inTOout


        TYPE (Type_dnS) :: dnR,dnTheta,dnCosTheta,dnSinTheta
        TYPE (Type_dnS) :: dnX,dnY,dnZ
        TYPE (Type_dnS) :: dnX2,dnY2,dnZ2,dnR2

        TYPE (Type_dnS) :: dnQ1,dnQ2,dnQ3

        integer :: i1,i2,i3
!----- for debuging ----------------------------------
       character (len=*),parameter :: name_sub='calc_ThreeDTransfo'
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

      IF (.NOT. associated(ThreeDTransfo)) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' ThreeDTransfo is NOT associated'
        write(out_unitp,*) ' Check source !!'
        STOP
      END IF

      CALL check_alloc_dnVec(dnQin,'dnQin',name_sub)
      CALL check_alloc_dnVec(dnQout,'dnQout',name_sub)

      IF (inTOout) THEN
        CALL sub_dnVec1_TO_dnVec2(dnQin,dnQout,nderiv)

        SELECT CASE (ThreeDTransfo%Type_3D)
        CASE (0) ! identity
          ! nothing
        CASE (1) ! polar
          STOP 'polar not yet'
        CASE (2)
          STOP 'spherical not yet'
        CASE (3)
          i1 = ThreeDTransfo%list_ThreeD_coord(1)
          i2 = ThreeDTransfo%list_ThreeD_coord(2)
          i3 = ThreeDTransfo%list_ThreeD_coord(3)

          !write(out_unitp,*) 'i1,i2,i3',i1,i2,3
          CALL sub_dnVec_TO_dnS(dnQin,dnX,i1,nderiv)
          CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnX,dnX,dnX2,nderiv)

          CALL sub_dnVec_TO_dnS(dnQin,dnY,i2,nderiv)
          CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnY,dnY,dnY2,nderiv)

          CALL sub_dnVec_TO_dnS(dnQin,dnR,i3,nderiv)
          CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnR,dnR,dnR2,nderiv)


          ! Z2 = R2 - X2 - Y2
          CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnR2,ONE,dnX2,-ONE,dnZ2,nderiv)
          CALL sub_dnS1_wPLUS_dnS2_TO_dnS2(dnY2,-ONE,dnZ2,ONE,nderiv)

          CALL sub_dnS1_TO_dntR2(dnZ2,dnZ,91,nderiv) ! dnZ=sqrt(dnZ2)

          ! Change the sign of dnY as function of the sign of dnR%d0
          IF (dnR%d0 < 0) CALL sub_Weight_dnS(dnZ,-ONE,nderiv)

          ! transfert dnY in dnQout
          CALL sub_dnS_TO_dnVec(dnZ,dnQout,i3,nderiv)

          CALL dealloc_dnSVM(dnR)
          CALL dealloc_dnSVM(dnX)
          CALL dealloc_dnSVM(dnY)
          CALL dealloc_dnSVM(dnY)
          CALL dealloc_dnSVM(dnR2)
          CALL dealloc_dnSVM(dnX2)
          CALL dealloc_dnSVM(dnY2)
          CALL dealloc_dnSVM(dnY2)

        CASE default ! ERROR: wrong transformation !
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' The type of threeD transformation is UNKNOWN: ',ThreeDTransfo%Type_3D
          write(out_unitp,*) ' The possible values are:'
          write(out_unitp,*) '    0: identity: x,y,z'
          write(out_unitp,*) '    1: Polar: R,theta,z'
          write(out_unitp,*) '    2: Spherical: R,theta,phi'
          write(out_unitp,*) '    3: x,R,z'
          write(out_unitp,*) ' Check your data !!'
          STOP

        END SELECT

      ELSE
        STOP 'inTOout=f not yet'
      END IF


!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------
      END SUBROUTINE calc_ThreeDTransfo

      SUBROUTINE ThreeDTransfo1TOThreeDTransfo2(ThreeDTransfo1,ThreeDTransfo2)
        TYPE (Type_ThreeDTransfo),pointer, intent(in)    :: ThreeDTransfo1
        TYPE (Type_ThreeDTransfo),pointer, intent(inout) :: ThreeDTransfo2

        IF (.NOT. associated(ThreeDTransfo1)) RETURN

        CALL alloc_ThreeDTransfo(ThreeDTransfo2)

        ThreeDTransfo2%Type_3D              = ThreeDTransfo1%Type_3D
        ThreeDTransfo2%name_Transfo_3D      = ThreeDTransfo1%name_Transfo_3D
        ThreeDTransfo2%list_ThreeD_coord(:) = ThreeDTransfo1%list_ThreeD_coord(:)

      END SUBROUTINE ThreeDTransfo1TOThreeDTransfo2


      END MODULE mod_ThreeDTransfo
