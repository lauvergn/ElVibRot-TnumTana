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
      MODULE mod_FlexibleTransfo
      use mod_system
      use mod_dnSVM, only: type_dnvec, type_dns, check_alloc_dnvec,    &
                           alloc_dnsvm, sub_dnvec1_to_dnvec2_withivec, &
                           dealloc_dnsvm, write_dnvec
      IMPLICIT NONE

      PRIVATE

      !!@description: TODO
      !!@param: TODO
      TYPE Type_FlexibleTransfo
        integer              :: nb_flex_act  = 0
        integer, allocatable :: list_flex(:)
        integer, allocatable :: list_act(:)

        CONTAINS
            PROCEDURE :: Read_FlexibleTransfo
            GENERIC :: QtransfoRead => Read_FlexibleTransfo
      END TYPE Type_FlexibleTransfo

      PUBLIC :: Type_FlexibleTransfo, Read_FlexibleTransfo
      PUBLIC :: alloc_FlexibleTransfo, dealloc_FlexibleTransfo
      PUBLIC :: FlexibleTransfo1TOFlexibleTransfo2
      PUBLIC :: calc_FlexibleTransfo

      CONTAINS

!=======================================================================
!     Felxible transfo
!=======================================================================
      SUBROUTINE alloc_FlexibleTransfo(FlexibleTransfo,nb_Qin)
      TYPE (Type_FlexibleTransfo), intent(inout) :: FlexibleTransfo
      integer                    , intent(in)    :: nb_Qin
      integer :: err_mem,memory

      CALL dealloc_FlexibleTransfo(FlexibleTransfo)


      CALL alloc_NParray(FlexibleTransfo%list_act,(/nb_Qin/),           &
                        "FlexibleTransfo%list_act","alloc_FlexibleTransfo")
      FlexibleTransfo%list_act(:) = 0
      CALL alloc_NParray(FlexibleTransfo%list_flex,(/nb_Qin/),          &
                        "FlexibleTransfo%list_flex","alloc_FlexibleTransfo")
      FlexibleTransfo%list_flex(:) = 0

      END SUBROUTINE alloc_FlexibleTransfo

      SUBROUTINE dealloc_FlexibleTransfo(FlexibleTransfo)
      TYPE (Type_FlexibleTransfo), intent(inout) :: FlexibleTransfo
      integer :: err_mem,memory

      FlexibleTransfo%nb_flex_act = 0

      IF (allocated(FlexibleTransfo%list_act) ) THEN
        CALL dealloc_NParray(FlexibleTransfo%list_act,                  &
                            "FlexibleTransfo%list_act","dealloc_FlexibleTransfo")
      END IF
      IF (allocated(FlexibleTransfo%list_flex) ) THEN
        CALL dealloc_NParray(FlexibleTransfo%list_flex,                 &
                            "FlexibleTransfo%list_flex","dealloc_FlexibleTransfo")
      END IF

      END SUBROUTINE dealloc_FlexibleTransfo

      SUBROUTINE FlexibleTransfo1TOFlexibleTransfo2(FlexibleTransfo1,FlexibleTransfo2)

      !-- FlexibleTransfo --------------------------------------
      TYPE (Type_FlexibleTransfo), intent(in)    :: FlexibleTransfo1
      TYPE (Type_FlexibleTransfo), intent(inout) :: FlexibleTransfo2

      integer :: nb_Qin
      character (len=*), parameter :: name_sub = 'FlexibleTransfo1TOFlexibleTransfo2'

      nb_Qin = size(FlexibleTransfo1%list_flex)

      CALL alloc_FlexibleTransfo(FlexibleTransfo2,nb_Qin)

      FlexibleTransfo2%nb_flex_act = FlexibleTransfo1%nb_flex_act
      FlexibleTransfo2%list_flex   = FlexibleTransfo1%list_flex
      FlexibleTransfo2%list_act    = FlexibleTransfo1%list_act

      END SUBROUTINE FlexibleTransfo1TOFlexibleTransfo2

      SUBROUTINE Read_FlexibleTransfo(FlexibleTransfo,nb_Qin,list_flex)

      !TYPE (Type_FlexibleTransfo), intent(inout) :: FlexibleTransfo
      CLASS (Type_FlexibleTransfo), intent(inout) :: FlexibleTransfo

      integer, intent(in) :: nb_Qin
      integer, intent(in), optional :: list_flex(:)


      integer :: i,it,nb_flex_act,err,nbcol

      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Read_FlexibleTransfo'


      CALL alloc_FlexibleTransfo(FlexibleTransfo,nb_Qin)

      IF (present(list_flex)) THEN
        IF (size(list_flex) /= nb_Qin) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  The "list_flex" has the wrong size'
          write(out_unitp,*) '  size(list_flex)',size(list_flex)
          write(out_unitp,*) '           nb_Qin',nb_Qin
          write(out_unitp,*) ' Check the FORTRAN source !!'
          STOP
        END IF
        FlexibleTransfo%list_flex(:) = list_flex(:)
      ELSE
        read(in_unitp,*,IOSTAT=err) FlexibleTransfo%list_flex(:)
        IF (err /= 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  while reading "list_flex"'
          write(out_unitp,*) '  end of file or end of record'
          write(out_unitp,*) ' Check your data !!'
          STOP
        END IF
      END IF
      nb_flex_act = count(FlexibleTransfo%list_flex == 1)
      FlexibleTransfo%nb_flex_act = nb_flex_act
      nb_flex_act = 0
      DO i=1,nb_Qin
        IF (FlexibleTransfo%list_flex(i) == 1) THEN
          nb_flex_act = nb_flex_act + 1
          FlexibleTransfo%list_act(nb_flex_act) = i
        END IF
      END DO
      IF (nb_flex_act < 1) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' nb_flex_act is smaller than 1 !!'
        write(out_unitp,*) ' Check your data !!'
        write(out_unitp,*) 'list_flex',FlexibleTransfo%list_flex(:)
        STOP
      END IF

      write(out_unitp,*) 'nb_flex_act',nb_flex_act,':',                 &
                 FlexibleTransfo%list_act(1:nb_flex_act)
      write(out_unitp,*) 'list_flex: ',FlexibleTransfo%list_flex(:)

      END SUBROUTINE Read_FlexibleTransfo

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE calc_FlexibleTransfo(dnQin,dnQout,FlexibleTransfo,nderiv,inTOout)

        TYPE (Type_dnVec), intent(inout)        :: dnQin,dnQout
        TYPE (Type_flexibleTransfo), intent(in) :: FlexibleTransfo
        integer, intent(in)                     :: nderiv
        logical                                 :: inTOout

        TYPE (Type_dnS)   :: dnQeq

        integer :: nb_flex_act

        integer :: list_act(FlexibleTransfo%nb_flex_act)
        integer :: i,j,k,iact,jact,kact,id,jd,kd
        real (kind=Rkind) :: Qact_flex(FlexibleTransfo%nb_flex_act)

        integer :: iQ,it=0

!----- for debuging ----------------------------------
       character (len=*),parameter :: name_sub='calc_FlexibleTransfo'
       logical, parameter :: debug=.FALSE.
!        logical, parameter :: debug=.TRUE.
!----- for debuging ----------------------------------

!---------------------------------------------------------------------


      nb_flex_act = FlexibleTransfo%nb_flex_act
      list_act(:) = FlexibleTransfo%list_act(1:nb_flex_act)

      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_flex_act',nb_flex_act
        write(out_unitp,*) 'list_act ',FlexibleTransfo%list_act
        write(out_unitp,*) 'list_flex',FlexibleTransfo%list_flex

        CALL flush_perso(out_unitp)
      END IF
!---------------------------------------------------------------------

       CALL check_alloc_dnVec(dnQin,'dnQin',name_sub)
       CALL check_alloc_dnVec(dnQout,'dnQout',name_sub)

       CALL alloc_dnSVM(dnQeq,nb_flex_act,nderiv)

       Qact_flex(:) = dnQin%d0(list_act)

       IF (inTOout) THEN
         !write(out_unitp,*) 'list_act,Qact_flex',list_act,Qact_flex
         DO iQ=1,dnQin%nb_var_vec
           CALL sub_dnVec1_TO_dnVec2_WithIvec(dnQin,dnQout,iQ,nderiv)

           IF (FlexibleTransfo%list_flex(iQ) == 20) THEN

             CALL calc_dnQflex(iQ,dnQeq,Qact_flex,nb_flex_act,nderiv,it)
             !write(out_unitp,*) 'dnQout%d0(iQ)',iQ,dnQout%d0(iQ)
             !write(out_unitp,*) 'dnQeq,iQ',iQ
             !CALL Write_dnSVM(dnQeq,nderiv)

             dnQout%d0(iQ) = dnQin%d0(iQ) + dnQeq%d0

             IF (nderiv > 0) THEN
               DO i=1,nb_flex_act
                 iact = list_act(i)
                 dnQout%d1(iQ,:) = dnQout%d1(iQ,:) +          &
                                          dnQeq%d1(i) * dnQin%d1(iact,:)
               END DO
               !write(out_unitp,*) 'first der all',iQ,dnQout%d1(iQ,:)
             END IF

             IF (nderiv > 1) THEN
               DO i=1,nb_flex_act
                 iact = list_act(i)
                 dnQout%d2(iQ,:,:) = dnQout%d2(iQ,:,:) +      &
                                       dnQeq%d1(i) * dnQin%d2(iact,:,:)
               END DO
               !write(out_unitp,*) 'first der all',iQ,dnQout%d1(iQ,:)

               DO i=1,nb_flex_act
               DO j=1,nb_flex_act
                 iact = list_act(i)
                 jact = list_act(j)

                 DO id=1,dnQout%nb_var_deriv
                 DO jd=1,dnQout%nb_var_deriv
                   dnQout%d2(iQ,id,jd) =                            &
                                             dnQout%d2(iQ,id,jd) +&
                     dnQeq%d2(i,j) * dnQin%d1(iact,id) * dnQin%d1(jact,jd)
                 END DO
                 END DO

               END DO
               END DO
             END IF
             IF (nderiv > 2) THEN

               DO i=1,nb_flex_act
                 iact = list_act(i)
                 dnQout%d3(iQ,:,:,:) = dnQout%d3(iQ,:,:,:) +  &
                                     dnQeq%d1(i) * dnQin%d3(iact,:,:,:)
               END DO

               DO i=1,nb_flex_act
               DO j=1,nb_flex_act
                 iact = list_act(i)
                 jact = list_act(j)

                 DO id=1,dnQout%nb_var_deriv
                 DO jd=1,dnQout%nb_var_deriv
                 DO kd=1,dnQout%nb_var_deriv

                   dnQout%d3(iQ,id,jd,kd) =                         &
                       dnQout%d3(iQ,id,jd,kd) + dnQeq%d2(i,j) * ( &
                          dnQin%d1(iact,id) * dnQin%d2(jact,jd,kd) +    &
                          dnQin%d1(iact,jd) * dnQin%d2(jact,kd,id) +    &
                          dnQin%d1(iact,kd) * dnQin%d2(jact,id,jd) )

                 END DO
                 END DO
                 END DO

               END DO
               END DO

               DO i=1,nb_flex_act
               DO j=1,nb_flex_act
               DO k=1,nb_flex_act
                 iact = list_act(i)
                 jact = list_act(j)
                 kact = list_act(k)

                 DO id=1,dnQout%nb_var_deriv
                 DO jd=1,dnQout%nb_var_deriv
                 DO kd=1,dnQout%nb_var_deriv

                     dnQout%d3(iQ,id,jd,kd) =                         &
                       dnQout%d3(iQ,id,jd,kd) + dnQeq%d3(i,j,k) * &
                            dnQin%d1(iact,id) *                         &
                            dnQin%d1(jact,jd) *                         &
                            dnQin%d1(kact,kd)
                 END DO
                 END DO
                 END DO

               END DO
               END DO
               END DO
             END IF
           END IF
         END DO
       ELSE
         !write(out_unitp,*) 'list_act,Qact_flex',list_act,Qact_flex
         DO iQ=1,dnQin%nb_var_vec
           CALL sub_dnVec1_TO_dnVec2_WithIvec(dnQout,dnQin,iQ,nderiv)

           IF (FlexibleTransfo%list_flex(iQ) == 20) THEN

             CALL calc_dnQflex(iQ,dnQeq,Qact_flex,nb_flex_act,nderiv,it)
             !write(out_unitp,*) 'dnQeq,iQ',iQ
             !CALL Write_dnSVM(dnQeq,nderiv)

             dnQin%d0(iQ) = dnQout%d0(iQ) - dnQeq%d0
           END IF

         END DO
       END IF
       CALL dealloc_dnSVM(dnQeq)

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'dnQout'
        CALL Write_dnVec(dnQout)
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
      !stop
!---------------------------------------------------------------------
      END SUBROUTINE calc_FlexibleTransfo

      END MODULE mod_FlexibleTransfo
