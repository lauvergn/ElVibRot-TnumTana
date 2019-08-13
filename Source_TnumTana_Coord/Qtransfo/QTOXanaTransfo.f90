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
      MODULE mod_QTOXanaTransfo
      use mod_system
      USE mod_dnSVM
      use mod_constant,  only: table_atom, get_mass_tnum
      IMPLICIT NONE

      PRIVATE

      !!@description: TODO
      !!@param: TODO
      TYPE Type_QTOXanaTransfo

        integer           :: ncart     = 0
        integer           :: ncart_act = 0

        integer           :: nat0      = 0
        integer           :: nat       = 0
        integer           :: nat_act   = 0

        integer           :: nb_var    = 0

        real (kind=Rkind),        pointer :: masses(:)     => null()
        integer, pointer                  :: Z(:)          => null()
        character (len=Name_len),pointer  :: symbole(:)    => null()
        integer, pointer                  :: type_Qin(:)   => null() ! TRUE pointer

      END TYPE Type_QTOXanaTransfo

      PUBLIC :: Type_QTOXanaTransfo, alloc_QTOXanaTransfo, dealloc_QTOXanaTransfo, &
                Read_QTOXanaTransfo, Write_QTOXanaTransfo, QTOXanaTransfo1TOQTOXanaTransfo2

      CONTAINS

!================================================================
!       Read QTOXana Transfo
!================================================================
      SUBROUTINE alloc_QTOXanaTransfo(QTOXanaTransfo)
      TYPE (Type_QTOXanaTransfo), intent(inout) :: QTOXanaTransfo

       character (len=*), parameter :: name_sub = 'alloc_QTOXanaTransfo'

!      write(out_unitp,*) 'BEGINNING ',name_sub
!      write(out_unitp,*) 'nat',QTOXanaTransfo%nat

       IF (QTOXanaTransfo%nat < 3) THEN
         write(out_unitp,*) ' ERROR in alloc_QTOXanaTransfo'
         write(out_unitp,*) ' wrong value of nat',QTOXanaTransfo%nat
         write(out_unitp,*) ' CHECK the source !!'
         STOP
       END IF

       CALL alloc_array(QTOXanaTransfo%Z,(/QTOXanaTransfo%nat/),        &
                       "QTOXanaTransfo%Z",name_sub)
       QTOXanaTransfo%Z(:) = 0

       CALL alloc_array(QTOXanaTransfo%symbole,(/QTOXanaTransfo%nat/),  &
              Name_len,"QTOXanaTransfo%symbole",name_sub)
       QTOXanaTransfo%symbole(:) = ""

       CALL alloc_array(QTOXanaTransfo%masses,(/QTOXanaTransfo%ncart/), &
                       "QTOXanaTransfo%masses",name_sub)
       QTOXanaTransfo%masses(:) = ZERO

!      write(out_unitp,*) 'END ',name_sub

      END SUBROUTINE alloc_QTOXanaTransfo

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE dealloc_QTOXanaTransfo(QTOXanaTransfo)

       TYPE (Type_QTOXanaTransfo), intent(inout) :: QTOXanaTransfo

       !write(out_unitp,*) 'BEGINNING dealloc_QTOXanaTransfo'; call flush_perso(out_unitp)

       IF (associated(QTOXanaTransfo%Z))  THEN
         CALL dealloc_array(QTOXanaTransfo%Z,                           &
                           "QTOXanaTransfo%Z","dealloc_QTOXanaTransfo")
       END IF

       IF (associated(QTOXanaTransfo%masses))  THEN
         CALL dealloc_array(QTOXanaTransfo%masses,                      &
                           "QTOXanaTransfo%masses","dealloc_QTOXanaTransfo")
       END IF

       IF (associated(QTOXanaTransfo%symbole))  THEN
         CALL dealloc_array(QTOXanaTransfo%symbole,                     &
                           "QTOXanaTransfo%symbole","dealloc_QTOXanaTransfo")
       END IF


        QTOXanaTransfo%ncart     = 0
        QTOXanaTransfo%ncart_act = 0
        QTOXanaTransfo%nat0      = 0
        QTOXanaTransfo%nat       = 0
        QTOXanaTransfo%nat_act   = 0
        QTOXanaTransfo%nb_var    = 0

       nullify(QTOXanaTransfo%type_Qin)

       !write(out_unitp,*) 'END dealloc_QTOXanaTransfo'; call flush_perso(out_unitp)

      END SUBROUTINE dealloc_QTOXanaTransfo

      SUBROUTINE Read_QTOXanaTransfo(QTOXanaTransfo,mendeleev)

       TYPE (Type_QTOXanaTransfo),intent(inout) :: QTOXanaTransfo
       TYPE (table_atom), intent(in)            :: mendeleev


       real (kind=Rkind)        :: at
       integer :: i
        character (len=Name_len), pointer :: name_at(:)

       !-----------------------------------------------------------------------
       integer :: err_mem,memory,err_io
       logical, parameter :: debug=.FALSE.
       !logical, parameter :: debug=.TRUE.
       character (len=*), parameter :: name_sub = 'Read_QTOXanaTransfo'
       !-----------------------------------------------------------------------
       IF (print_level > 1) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'nat0,nat ',QTOXanaTransfo%nat0,QTOXanaTransfo%nat
         write(out_unitp,*) 'nb_var   ',QTOXanaTransfo%nb_var
       END IF


          QTOXanaTransfo%type_Qin(:) = 0
          read(in_unitp,*,IOSTAT=err_io) QTOXanaTransfo%type_Qin(:)
          IF (err_io /= 0) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) '  while reading the type of the coordinates.'
            write(out_unitp,*) ' QTOXanaTransfo%type_Qin',QTOXanaTransfo%type_Qin(:)
            write(out_unitp,*) ' end of file or end of record'
            write(out_unitp,*) ' Check your data !!'
            STOP
          END IF

          CALL alloc_QTOXanaTransfo(QTOXanaTransfo)

          nullify(name_at)
          CALL alloc_array(name_at,(/QTOXanaTransfo%nat0/),Name_len,    &
                          "name_at",name_sub)
          read(in_unitp,*,IOSTAT=err_io) (name_at(i),i=1,QTOXanaTransfo%nat0)
          IF (err_io /= 0) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) '  while reading a mass.'
            write(out_unitp,*) ' end of file or end of record'
            write(out_unitp,*) ' Check your data !!'
            STOP
          END IF

          DO i=1,QTOXanaTransfo%nat0
            QTOXanaTransfo%Z(i) = -1
            QTOXanaTransfo%symbole(i) = name_at(i)
            at = get_mass_Tnum(mendeleev,Z=QTOXanaTransfo%Z(i),name=name_at(i))
            IF (print_level > 0) write(out_unitp,*) i,QTOXanaTransfo%Z(i),at

!            IF (at == ZERO) THEN
!              write(out_unitp,*) ' ERROR in ',name_sub
!              write(out_unitp,*) '  One mass is ZERO.'
!              write(out_unitp,*) '  It is not possible for this transformation'
!              write(out_unitp,*) ' Check your data'
!              STOP
!            END IF

            QTOXanaTransfo%masses( 3*i -2: 3*i -0) = at
          END DO
          CALL dealloc_array(name_at,"name_at",name_sub)

      IF (print_level > 1) write(out_unitp,*) 'END ',name_sub
      END SUBROUTINE Read_QTOXanaTransfo


      SUBROUTINE QTOXanaTransfo1TOQTOXanaTransfo2(QTOXanaTransfo1,QTOXanaTransfo2)

!      for the QTOXanarix and Tnum --------------------------------------
      TYPE (Type_QTOXanaTransfo), intent(in)    :: QTOXanaTransfo1
      TYPE (Type_QTOXanaTransfo), intent(inout) :: QTOXanaTransfo2

      character (len=*), parameter ::                                   &
                                name_sub = 'QTOXanaTransfo1TOQTOXanaTransfo2'

      CALL dealloc_QTOXanaTransfo(QTOXanaTransfo2)

      QTOXanaTransfo2%ncart        = QTOXanaTransfo1%ncart
      QTOXanaTransfo2%ncart_act    = QTOXanaTransfo1%ncart_act
      QTOXanaTransfo2%nat          = QTOXanaTransfo1%nat
      QTOXanaTransfo2%nat0         = QTOXanaTransfo1%nat0
      QTOXanaTransfo2%nat_act      = QTOXanaTransfo1%nat_act
      QTOXanaTransfo2%nb_var       = QTOXanaTransfo1%nb_var

      CALL alloc_QTOXanaTransfo(QTOXanaTransfo2)


      QTOXanaTransfo2%masses       = QTOXanaTransfo1%masses
      QTOXanaTransfo2%Z(:)         = QTOXanaTransfo1%Z(:)
      QTOXanaTransfo2%symbole(:)   = QTOXanaTransfo1%symbole(:)

      !write(out_unitp,*) 'END QTOXanaTransfo1TOQTOXanaTransfo2'

      END SUBROUTINE QTOXanaTransfo1TOQTOXanaTransfo2

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE Write_QTOXanaTransfo(QTOXanaTransfo)
      TYPE (Type_QTOXanaTransfo), intent(in) :: QTOXanaTransfo

      integer :: i
      character (len=*), parameter :: name_sub='Write_QTOXanaTransfo'


      write(out_unitp,*) 'BEGINNING ',name_sub

      write(out_unitp,*) 'ncart_act,ncart',                             &
                  QTOXanaTransfo%ncart_act,QTOXanaTransfo%ncart

      write(out_unitp,*) 'nat_act,nat0,nat,',                           &
                  QTOXanaTransfo%nat_act,QTOXanaTransfo%nat0,QTOXanaTransfo%nat

      write(out_unitp,*) 'nb_var',QTOXanaTransfo%nb_var

      write(out_unitp,*) 'masses : ',QTOXanaTransfo%masses(:)
      write(out_unitp,*) 'Z      : ',QTOXanaTransfo%Z(:)
      write(out_unitp,*) 'symbole: ',QTOXanaTransfo%symbole(:)

      write(out_unitp,*) 'END ',name_sub

      END SUBROUTINE Write_QTOXanaTransfo

      END MODULE mod_QTOXanaTransfo

