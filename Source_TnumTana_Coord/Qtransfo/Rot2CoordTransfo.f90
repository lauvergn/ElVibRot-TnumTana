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
      MODULE mod_Rot2CoordTransfo
      use mod_system, only: error_memo_allo, in_unitp, out_unitp, one
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

      CONTAINS

      SUBROUTINE alloc_Rot2CoordTransfo(Rot2CoordTransfo)
      TYPE (Type_Rot2CoordTransfo), pointer, intent(inout) :: Rot2CoordTransfo
      integer :: err_mem,memory

      IF (.NOT. associated(Rot2CoordTransfo)) THEN
        allocate(Rot2CoordTransfo,stat=err_mem)
        memory = 1
        CALL error_memo_allo(err_mem,memory,'Rot2CoordTransfo','alloc_Rot2CoordTransfo','Type_Rot2CoordTransfo')
      END IF

      END SUBROUTINE alloc_Rot2CoordTransfo

      SUBROUTINE dealloc_Rot2CoordTransfo(Rot2CoordTransfo)
      TYPE (Type_Rot2CoordTransfo), pointer, intent(inout) :: Rot2CoordTransfo
      integer :: err_mem,memory

      IF (.NOT. associated(Rot2CoordTransfo)) RETURN
      deallocate(Rot2CoordTransfo,stat=err_mem)
      memory = 1
      CALL error_memo_allo(err_mem,-memory,'Rot2CoordTransfo','dealloc_Rot2CoordTransfo','Type_Rot2CoordTransfo')
      nullify(Rot2CoordTransfo)

      END SUBROUTINE dealloc_Rot2CoordTransfo

      SUBROUTINE Read_Rot2CoordTransfo(Rot2CoordTransfo,nb_Qin)

      TYPE (Type_Rot2CoordTransfo), pointer, intent(inout) :: Rot2CoordTransfo
      integer, intent(in) :: nb_Qin

      integer :: num_Rot,list_2Coord(2)

       NAMELIST / Rot2Coord / num_Rot,list_2Coord


      integer :: err_io,err_mem,memory
      character (len=*), parameter :: name_sub='Read_Rot2CoordTransfo'

      CALL alloc_Rot2CoordTransfo(Rot2CoordTransfo)

      num_Rot         = 0
      list_2Coord(:)  = 0
      read(in_unitp,Rot2Coord,IOSTAT=err_io)
      IF (err_io /= 0) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) '  while reading "Rot2Coord" namelist'
         write(out_unitp,*) ' end of file or end of record'
         write(out_unitp,*) ' Check your data !!'
         STOP
      END IF
      Rot2CoordTransfo%num_Rot        = num_Rot
      Rot2CoordTransfo%list_2Coord(:) = list_2Coord(:)


      CALL Write_Rot2CoordTransfo(Rot2CoordTransfo)

      END SUBROUTINE Read_Rot2CoordTransfo

      SUBROUTINE Write_Rot2CoordTransfo(Rot2CoordTransfo)

      TYPE (Type_Rot2CoordTransfo), pointer, intent(in) :: Rot2CoordTransfo

      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Write_Rot2CoordTransfo'

      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) 'asso Rot2CoordTransfo:',associated(Rot2CoordTransfo)
      IF (associated(Rot2CoordTransfo)) THEN
        write(out_unitp,*) 'num_Rot:     ',Rot2CoordTransfo%num_Rot
        write(out_unitp,*) 'list_2Coord: ',Rot2CoordTransfo%list_2Coord(:)
      END IF
      write(out_unitp,*) 'END ',name_sub

      END SUBROUTINE Write_Rot2CoordTransfo

      SUBROUTINE calc_Rot2CoordTransfo(dnQin,dnQout,Rot2CoordTransfo,   &
                                    nderiv,inTOout)

        TYPE (Type_dnVec), intent(inout)              :: dnQin,dnQout
        TYPE (Type_Rot2CoordTransfo),pointer, intent(in) :: Rot2CoordTransfo
        integer, intent(in)                           :: nderiv
        logical, intent(in)                           :: inTOout


        TYPE (Type_dnS) :: dnTheta,dnCosTheta,dnSinTheta
        TYPE (Type_dnS) :: dnTemp

        TYPE (Type_dnS) :: dnQ1old,dnQ2old
        TYPE (Type_dnS) :: dnQ1new,dnQ2new

!----- for debuging ----------------------------------
       character (len=*),parameter :: name_sub='calc_Rot2CoordTransfo'
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

        CALL sub_dnVec_TO_dnS(dnQin,dnTheta,Rot2CoordTransfo%num_Rot,nderiv)

        CALL sub_dnS1_TO_dntR2(dnTheta,dnCosTheta,2,nderiv)
        CALL sub_dnS1_TO_dntR2(dnTheta,dnSinTheta,3,nderiv)


        CALL sub_dnVec_TO_dnS(dnQin,dnQ1old,Rot2CoordTransfo%list_2Coord(1),nderiv)
        CALL sub_dnVec_TO_dnS(dnQin,dnQ2old,Rot2CoordTransfo%list_2Coord(2),nderiv)


        CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnQ1old,dnCosTheta,dnQ1new,nderiv)
        CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnQ2old,dnSinTheta,dnTemp,nderiv)
        CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnQ1new,ONE,dnTemp,-ONE,dnQ1new,nderiv)

        CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnQ2old,dnCosTheta,dnQ2new,nderiv)
        CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnQ1old,dnSinTheta,dnTemp,nderiv)
        CALL sub_dnS1_PLUS_dnS2_TO_dnS2(dnTemp,dnQ2new,nderiv)



        ! transfert dnY in dnQout
        CALL sub_dnS_TO_dnVec(dnQ1new,dnQout,Rot2CoordTransfo%list_2Coord(1),nderiv)
        CALL sub_dnS_TO_dnVec(dnQ2new,dnQout,Rot2CoordTransfo%list_2Coord(2),nderiv)

      ELSE
        CALL sub_dnVec1_TO_dnVec2(dnQout,dnQin,nderiv)

        CALL sub_dnVec_TO_dnS(dnQout,dnTheta,Rot2CoordTransfo%num_Rot,nderiv)

        CALL sub_dnS1_TO_dntR2(dnTheta,dnCosTheta,2,nderiv)
        CALL sub_dnS1_TO_dntR2(dnTheta,dnSinTheta,3,nderiv)


        CALL sub_dnVec_TO_dnS(dnQout,dnQ1old,Rot2CoordTransfo%list_2Coord(1),nderiv)
        CALL sub_dnVec_TO_dnS(dnQout,dnQ2old,Rot2CoordTransfo%list_2Coord(2),nderiv)


        CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnQ1old,dnCosTheta,dnQ1new,nderiv)
        CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnQ2old,dnSinTheta,dnTemp,nderiv)
        CALL sub_dnS1_PLUS_dnS2_TO_dnS2(dnTemp,dnQ1new,nderiv)

        CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnQ2old,dnCosTheta,dnQ2new,nderiv)
        CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnQ1old,dnSinTheta,dnTemp,nderiv)
        CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnQ2new,ONE,dnTemp,-ONE,dnQ2new,nderiv)


        ! transfert dnY in dnQout
        CALL sub_dnS_TO_dnVec(dnQ1new,dnQin,Rot2CoordTransfo%list_2Coord(1),nderiv)
        CALL sub_dnS_TO_dnVec(dnQ2new,dnQin,Rot2CoordTransfo%list_2Coord(2),nderiv)



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
      END SUBROUTINE calc_Rot2CoordTransfo

      SUBROUTINE Rot2CoordTransfo1TORot2CoordTransfo2(Rot2CoordTransfo1,Rot2CoordTransfo2)
        TYPE (Type_Rot2CoordTransfo),pointer, intent(in)    :: Rot2CoordTransfo1
        TYPE (Type_Rot2CoordTransfo),pointer, intent(inout) :: Rot2CoordTransfo2

        IF (.NOT. associated(Rot2CoordTransfo1)) RETURN

        CALL alloc_Rot2CoordTransfo(Rot2CoordTransfo2)

        Rot2CoordTransfo2%num_Rot        = Rot2CoordTransfo1%num_Rot
        Rot2CoordTransfo2%list_2Coord(:) = Rot2CoordTransfo1%list_2Coord(:)

      END SUBROUTINE Rot2CoordTransfo1TORot2CoordTransfo2


      END MODULE mod_Rot2CoordTransfo
