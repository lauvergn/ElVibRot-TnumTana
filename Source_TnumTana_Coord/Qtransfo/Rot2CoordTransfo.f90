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
      MODULE mod_Rot2CoordTransfo
      use mod_system
      USE mod_dnSVM
      IMPLICIT NONE

      PRIVATE

      TYPE Type_Rot2CoordTransfo
        integer                  :: num_Rot
        integer                  :: list_2Coord(2) = (/ 0,0 /)
      END TYPE Type_Rot2CoordTransfo

      PUBLIC :: Type_Rot2CoordTransfo, alloc_Rot2CoordTransfo, dealloc_Rot2CoordTransfo
      PUBLIC :: Read_Rot2CoordTransfo, Write_Rot2CoordTransfo, calc_Rot2CoordTransfo
      PUBLIC :: Rot2CoordTransfo1TORot2CoordTransfo2

      INTERFACE calc_Rot2CoordTransfo
       !MODULE PROCEDURE calc_Rot2CoordTransfo_new
       MODULE PROCEDURE calc_Rot2CoordTransfo_old
      END INTERFACE

      CONTAINS

      SUBROUTINE alloc_Rot2CoordTransfo(Rot2CoordTransfo,nb_transfo)
      TYPE (Type_Rot2CoordTransfo), pointer, intent(inout) :: Rot2CoordTransfo(:)
      integer,                               intent(in)    :: nb_transfo

      integer :: err_mem,memory

      IF (associated(Rot2CoordTransfo)) THEN
        CALL dealloc_Rot2CoordTransfo(Rot2CoordTransfo)
      END IF
      IF (nb_transfo < 1) RETURN

      allocate(Rot2CoordTransfo(nb_transfo),stat=err_mem)
      memory = nb_transfo
      CALL error_memo_allo(err_mem,memory,'Rot2CoordTransfo','alloc_Rot2CoordTransfo','Type_Rot2CoordTransfo')

      END SUBROUTINE alloc_Rot2CoordTransfo

      SUBROUTINE dealloc_Rot2CoordTransfo(Rot2CoordTransfo)
      TYPE (Type_Rot2CoordTransfo), pointer, intent(inout) :: Rot2CoordTransfo(:)
      integer :: err_mem,memory

      IF (.NOT. associated(Rot2CoordTransfo)) RETURN
      memory = size(Rot2CoordTransfo)

      deallocate(Rot2CoordTransfo,stat=err_mem)
      CALL error_memo_allo(err_mem,-memory,'Rot2CoordTransfo','dealloc_Rot2CoordTransfo','Type_Rot2CoordTransfo')
      nullify(Rot2CoordTransfo)

      END SUBROUTINE dealloc_Rot2CoordTransfo

      SUBROUTINE Read_Rot2CoordTransfo(Rot2CoordTransfo,nb_transfo,nb_Qin)

      TYPE (Type_Rot2CoordTransfo), pointer, intent(inout) :: Rot2CoordTransfo(:)
      integer,                               intent(in)    :: nb_transfo,nb_Qin

      integer :: num_Rot,list_2Coord(2)
      integer :: i

       NAMELIST / Rot2Coord / num_Rot,list_2Coord


      integer :: err_io,err_mem,memory
      character (len=*), parameter :: name_sub='Read_Rot2CoordTransfo'

      CALL alloc_Rot2CoordTransfo(Rot2CoordTransfo,nb_transfo)

      DO i=1,nb_transfo
        num_Rot         = 0
        list_2Coord(:)  = 0
        read(in_unitp,Rot2Coord,IOSTAT=err_io)
        IF (err_io /= 0) THEN
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) '  while reading the ',int_TO_char(i),'th "Rot2Coord" namelist'
           write(out_unitp,*) ' end of file or end of record'
           write(out_unitp,*) ' Check your data !!'
           STOP
        END IF
        Rot2CoordTransfo(i)%num_Rot        = num_Rot
        Rot2CoordTransfo(i)%list_2Coord(:) = list_2Coord(:)

        !test if num_Rot are indentical. Most of the time it should be identical
        IF (Rot2CoordTransfo(i)%num_Rot /= Rot2CoordTransfo(1)%num_Rot) THEN
           write(out_unitp,*) ' WARNING in ',name_sub
           write(out_unitp,*) '  The "num_Rot" are not identical!!'
           write(out_unitp,*) ' Rot2Coord transformation:',i
           write(out_unitp,*) ' First   "num_Rot"',Rot2CoordTransfo(1)%num_Rot
           write(out_unitp,*) ' Current "num_Rot"',Rot2CoordTransfo(i)%num_Rot
        END IF

      END DO

      CALL Write_Rot2CoordTransfo(Rot2CoordTransfo)


      END SUBROUTINE Read_Rot2CoordTransfo

      SUBROUTINE Write_Rot2CoordTransfo(Rot2CoordTransfo)

      TYPE (Type_Rot2CoordTransfo), pointer, intent(in) :: Rot2CoordTransfo(:)

      integer :: i
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Write_Rot2CoordTransfo'

      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) 'asso Rot2CoordTransfo:',associated(Rot2CoordTransfo)
      write(out_unitp,*) 'size Rot2CoordTransfo:',size(Rot2CoordTransfo)
      IF (associated(Rot2CoordTransfo)) THEN
        DO i=1,size(Rot2CoordTransfo)
          write(out_unitp,*) 'Rot2CoordTransfo: ',i
          write(out_unitp,*) 'num_Rot:          ',Rot2CoordTransfo(i)%num_Rot
          write(out_unitp,*) 'list_2Coord:      ',Rot2CoordTransfo(i)%list_2Coord(:)
        END DO
      END IF
      write(out_unitp,*) 'END ',name_sub

      END SUBROUTINE Write_Rot2CoordTransfo

!      SUBROUTINE calc_Rot2CoordTransfo_new(dnQin,dnQout,Rot2CoordTransfo,   &
!                                       nderiv,inTOout)
!      USE mod_QML_dnS
!
!        TYPE (Type_dnVec), intent(inout)              :: dnQin,dnQout
!        TYPE (Type_Rot2CoordTransfo),pointer, intent(in) :: Rot2CoordTransfo(:)
!        integer, intent(in)                           :: nderiv
!        logical, intent(in)                           :: inTOout
!
!        TYPE (dnS_t) :: dnTheta,dnCosTheta,dnSinTheta
!
!        TYPE (dnS_t) :: dnQ1old,dnQ2old
!        TYPE (dns_t) :: dnQ1new,dnQ2new
!
!        integer :: i
!
!!----- for debuging ----------------------------------
!       character (len=*),parameter :: name_sub='calc_Rot2CoordTransfo_new'
!       logical, parameter :: debug=.FALSE.
!       !logical, parameter :: debug=.TRUE.
!!----- for debuging ----------------------------------
!
!
!!---------------------------------------------------------------------
!      IF (debug) THEN
!        write(out_unitp,*) 'BEGINNING ',name_sub
!        IF (inTOout) THEN
!          write(out_unitp,*) 'dnQin'
!          CALL Write_dnSVM(dnQin,nderiv)
!        ELSE
!          write(out_unitp,*) 'dnQout'
!          CALL Write_dnSVM(dnQout,nderiv)
!        END IF
!      END IF
!!---------------------------------------------------------------------
!
!      IF (.NOT. associated(Rot2CoordTransfo)) THEN
!        write(out_unitp,*) ' ERROR in ',name_sub
!        write(out_unitp,*) ' Rot2CoordTransfo is NOT associated'
!        write(out_unitp,*) ' Check source !!'
!        STOP
!      END IF
!
!      CALL check_alloc_dnVec(dnQin,'dnQin',name_sub)
!      CALL check_alloc_dnVec(dnQout,'dnQout',name_sub)
!
!      IF (inTOout) THEN
!        CALL sub_dnVec1_TO_dnVec2(dnQin,dnQout,nderiv)
!
!        DO i=1,size(Rot2CoordTransfo)
!          ! get from dnQin
!          CALL sub_dnVec_TO_dnSt(dnQin,dnTheta,Rot2CoordTransfo(i)%num_Rot)
!
!          CALL sub_dnVec_TO_dnSt(dnQin,dnQ1old,Rot2CoordTransfo(i)%list_2Coord(1))
!          CALL sub_dnVec_TO_dnSt(dnQin,dnQ2old,Rot2CoordTransfo(i)%list_2Coord(2))
!
!
!          dnCosTheta = cos(dnTheta)
!          dnSinTheta = sin(dnTheta)
!
!          dnQ1new = dnCosTheta * dnQ1old - dnSinTheta * dnQ2old
!          dnQ2new = dnSinTheta * dnQ1old + dnCosTheta * dnQ2old
!
!          ! transfert in dnQout
!          CALL sub_dnSt_TO_dnVec(dnQ1new,dnQout,Rot2CoordTransfo(i)%list_2Coord(1))
!          CALL sub_dnSt_TO_dnVec(dnQ2new,dnQout,Rot2CoordTransfo(i)%list_2Coord(2))
!
!        END DO
!
!      ELSE
!        CALL sub_dnVec1_TO_dnVec2(dnQout,dnQin,nderiv)
!
!        DO i=1,size(Rot2CoordTransfo)
!          ! get from dnQout
!          CALL sub_dnVec_TO_dnSt(dnQout,dnTheta,Rot2CoordTransfo(i)%num_Rot)
!
!          CALL sub_dnVec_TO_dnSt(dnQout,dnQ1old,Rot2CoordTransfo(i)%list_2Coord(1))
!          CALL sub_dnVec_TO_dnSt(dnQout,dnQ2old,Rot2CoordTransfo(i)%list_2Coord(2))
!
!
!          dnCosTheta = cos(dnTheta)
!          dnSinTheta = sin(dnTheta)
!
!          dnQ1new =  dnCosTheta * dnQ1old + dnSinTheta * dnQ2old
!          dnQ2new = -dnSinTheta * dnQ1old + dnCosTheta * dnQ2old
!
!          ! transfert in dnQin
!          CALL sub_dnSt_TO_dnVec(dnQ1new,dnQin,Rot2CoordTransfo(i)%list_2Coord(1))
!          CALL sub_dnSt_TO_dnVec(dnQ2new,dnQin,Rot2CoordTransfo(i)%list_2Coord(2))
!
!        END DO
!
!
!
!      END IF
!
!      CALL QML_dealloc_dnS(dnTheta)
!      CALL QML_dealloc_dnS(dnCosTheta)
!      CALL QML_dealloc_dnS(dnSinTheta)
!      CALL QML_dealloc_dnS(dnQ1old)
!      CALL QML_dealloc_dnS(dnQ2old)
!      CALL QML_dealloc_dnS(dnQ1new)
!      CALL QML_dealloc_dnS(dnQ2new)
!
!!---------------------------------------------------------------------
!      IF (debug) THEN
!        IF (inTOout) THEN
!          write(out_unitp,*)
!          write(out_unitp,*) 'dnQout'
!          CALL Write_dnSVM(dnQout,nderiv)
!        ELSE
!          write(out_unitp,*) 'dnQin'
!          CALL Write_dnSVM(dnQin,nderiv)
!        END IF
!        write(out_unitp,*) 'END ',name_sub
!      END IF
!!---------------------------------------------------------------------
!      END SUBROUTINE calc_Rot2CoordTransfo_new

      SUBROUTINE calc_Rot2CoordTransfo_old(dnQin,dnQout,Rot2CoordTransfo,   &
                                    nderiv,inTOout)

        TYPE (Type_dnVec), intent(inout)              :: dnQin,dnQout
        TYPE (Type_Rot2CoordTransfo),pointer, intent(in) :: Rot2CoordTransfo(:)
        integer, intent(in)                           :: nderiv
        logical, intent(in)                           :: inTOout


        TYPE (Type_dnS) :: dnTheta,dnCosTheta,dnSinTheta
        TYPE (Type_dnS) :: dnTemp

        TYPE (Type_dnS) :: dnQ1old,dnQ2old
        TYPE (Type_dnS) :: dnQ1new,dnQ2new

        integer :: i

!----- for debuging ----------------------------------
       character (len=*),parameter :: name_sub='calc_Rot2CoordTransfo_old'
       logical, parameter :: debug=.FALSE.
       !logical, parameter :: debug=.TRUE.
!----- for debuging ----------------------------------


!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'dnQin'
        CALL Write_dnSVM(dnQin,nderiv)
      END IF
!---------------------------------------------------------------------

      IF (.NOT. associated(Rot2CoordTransfo)) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' Rot2CoordTransfo is NOT associated'
        write(out_unitp,*) ' Check source !!'
        STOP
      END IF

      CALL check_alloc_dnVec(dnQin,'dnQin',name_sub)
      CALL check_alloc_dnVec(dnQout,'dnQout',name_sub)

      IF (inTOout) THEN
        CALL sub_dnVec1_TO_dnVec2(dnQin,dnQout,nderiv)

        DO i=1,size(Rot2CoordTransfo)
          CALL sub_dnVec_TO_dnS(dnQin,dnTheta,Rot2CoordTransfo(i)%num_Rot,nderiv)

          CALL sub_dnS1_TO_dntR2(dnTheta,dnCosTheta,2,nderiv)
          CALL sub_dnS1_TO_dntR2(dnTheta,dnSinTheta,3,nderiv)


          CALL sub_dnVec_TO_dnS(dnQin,dnQ1old,Rot2CoordTransfo(i)%list_2Coord(1),nderiv)
          CALL sub_dnVec_TO_dnS(dnQin,dnQ2old,Rot2CoordTransfo(i)%list_2Coord(2),nderiv)


          CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnQ1old,dnCosTheta,dnQ1new,nderiv)
          CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnQ2old,dnSinTheta,dnTemp,nderiv)
          CALL sub_dnS1_wPLUS_dnS2_TO_dnS2(dnTemp,-ONE,dnQ1new,ONE,nderiv)

          CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnQ2old,dnCosTheta,dnQ2new,nderiv)
          CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnQ1old,dnSinTheta,dnTemp,nderiv)
          CALL sub_dnS1_PLUS_dnS2_TO_dnS2(dnTemp,dnQ2new,nderiv)

          ! transfert dnY in dnQout
          CALL sub_dnS_TO_dnVec(dnQ1new,dnQout,Rot2CoordTransfo(i)%list_2Coord(1),nderiv)
          CALL sub_dnS_TO_dnVec(dnQ2new,dnQout,Rot2CoordTransfo(i)%list_2Coord(2),nderiv)

        END DO

      ELSE
        CALL sub_dnVec1_TO_dnVec2(dnQout,dnQin,nderiv)

        DO i=1,size(Rot2CoordTransfo)

          CALL sub_dnVec_TO_dnS(dnQout,dnTheta,Rot2CoordTransfo(i)%num_Rot,nderiv)

          CALL sub_dnS1_TO_dntR2(dnTheta,dnCosTheta,2,nderiv)
          CALL sub_dnS1_TO_dntR2(dnTheta,dnSinTheta,3,nderiv)


          CALL sub_dnVec_TO_dnS(dnQout,dnQ1old,Rot2CoordTransfo(i)%list_2Coord(1),nderiv)
          CALL sub_dnVec_TO_dnS(dnQout,dnQ2old,Rot2CoordTransfo(i)%list_2Coord(2),nderiv)


          CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnQ1old,dnCosTheta,dnQ1new,nderiv)
          CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnQ2old,dnSinTheta,dnTemp,nderiv)
          CALL sub_dnS1_PLUS_dnS2_TO_dnS2(dnTemp,dnQ1new,nderiv)

          CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnQ2old,dnCosTheta,dnQ2new,nderiv)
          CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnQ1old,dnSinTheta,dnTemp,nderiv)
          CALL sub_dnS1_wPLUS_dnS2_TO_dnS2(dnTemp,-ONE,dnQ2new,ONE,nderiv)


          ! transfert dnY in dnQout
          CALL sub_dnS_TO_dnVec(dnQ1new,dnQin,Rot2CoordTransfo(i)%list_2Coord(1),nderiv)
          CALL sub_dnS_TO_dnVec(dnQ2new,dnQin,Rot2CoordTransfo(i)%list_2Coord(2),nderiv)

        END DO



      END IF

      CALL dealloc_dnSVM(dnTemp)
      CALL dealloc_dnSVM(dnQ1new)
      CALL dealloc_dnSVM(dnQ1old)
      CALL dealloc_dnSVM(dnQ2old)
      CALL dealloc_dnSVM(dnQ2new)
      CALL dealloc_dnSVM(dnTheta)
      CALL dealloc_dnSVM(dnCosTheta)
      CALL dealloc_dnSVM(dnSinTheta)

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------
      END SUBROUTINE calc_Rot2CoordTransfo_old

      SUBROUTINE Rot2CoordTransfo1TORot2CoordTransfo2(Rot2CoordTransfo1,Rot2CoordTransfo2)
        TYPE (Type_Rot2CoordTransfo),pointer, intent(in)    :: Rot2CoordTransfo1(:)
        TYPE (Type_Rot2CoordTransfo),pointer, intent(inout) :: Rot2CoordTransfo2(:)

        integer :: i

        IF (.NOT. associated(Rot2CoordTransfo1)) RETURN

        CALL alloc_Rot2CoordTransfo(Rot2CoordTransfo2,size(Rot2CoordTransfo1))

        DO i=1,size(Rot2CoordTransfo1)
          Rot2CoordTransfo2(i)%num_Rot        = Rot2CoordTransfo1(i)%num_Rot
          Rot2CoordTransfo2(i)%list_2Coord(:) = Rot2CoordTransfo1(i)%list_2Coord(:)
        END DO

      END SUBROUTINE Rot2CoordTransfo1TORot2CoordTransfo2

      END MODULE mod_Rot2CoordTransfo
