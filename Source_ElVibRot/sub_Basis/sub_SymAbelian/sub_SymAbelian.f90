!===========================================================================
!===========================================================================
!This file is part of ElVibRot.
!
!    ElVibRot is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ElVibRot is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ElVibRot.  If not, see <http://www.gnu.org/licenses/>.
!
!    Copyright 2015  David Lauvergnat
!      with contributions of Mamadou Ndong, Josep Maria Luis
!
!    ElVibRot includes:
!        - Tnum-Tana under the GNU LGPL3 license
!        - Somme subroutines of John Burkardt under GNU LGPL license
!             http://people.sc.fsu.edu/~jburkardt/
!        - Somme subroutines of SHTOOLS written by Mark A. Wieczorek under BSD license
!             http://shtools.ipgp.fr
!===========================================================================
!===========================================================================
      MODULE mod_SymAbelian
      USE mod_system
      IMPLICIT NONE

        PRIVATE

        TYPE Type_SymAbelian
          PRIVATE
          integer :: Read_symab = -1
          integer :: nb_PER_symab(-1:7) = 0
          integer :: nb = 0                            ! size of tab_symab
          integer, allocatable :: tab_symab(:)         ! tab_symab(nb). We use 3 bits of each integer
                                                       ! 0: => (0,0,0)
                                                       ! 1: => (0,0,1)
                                                       ! 2: => (0,1,0)
                                                       ! ....
        END TYPE Type_SymAbelian

      INTERFACE alloc_array
        MODULE PROCEDURE alloc_array_OF_SymAbeliandim0
      END INTERFACE

      INTERFACE dealloc_array
        MODULE PROCEDURE dealloc_array_OF_SymAbeliandim0
      END INTERFACE

      PUBLIC  Type_SymAbelian, alloc_SymAbelian, dealloc_SymAbelian,     &
              SymAbelian1_TO_SymAbelian2, Write_SymAbelian,              &
              SymAbelian_IS_initialized, Get_Read_symabOFSymAbelian,     &
              Get_symabOFSymAbelian_AT_ib, Set_symabOFSymAbelian_AT_ib,  &
              Set_tab_symabOFSymAbelian_WITH_tab, Set_tab_SymAbelian,    &
              Set_nbPERsym_FROM_SymAbelian, Get_nbPERsym_FROM_SymAbelian,&
              Set_ReadsymabOFSymAbelian, Calc_symab1_EOR_symab2,         &
              Write_symab, WriteTOstring_symab

      PUBLIC  alloc_array, dealloc_array


      CONTAINS

       SUBROUTINE alloc_SymAbelian(P_SymAbelian,nb)
         TYPE (Type_SymAbelian), pointer, intent(inout) :: P_SymAbelian
         integer, intent(in)                   :: nb

         integer :: Read_symab

         character (len=*), parameter :: name_sub='alloc_SymAbelian'

         IF (nb < 1) THEN
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) '  WRONG paramter values: nb',nb
           write(out_unitp,*) '  CHECK the fortran !!'
           STOP
         END IF

         IF (associated(P_SymAbelian)) THEN
           Read_symab = P_SymAbelian%Read_symab
         ELSE
           Read_symab = -1
         END IF
         CALL dealloc_SymAbelian(P_SymAbelian)
         IF (.NOT. associated(P_SymAbelian)) THEN
           CALL alloc_array(P_SymAbelian,'P_SymAbelian',name_sub)
         END IF

         P_SymAbelian%nb          = nb
         P_SymAbelian%Read_symab  = Read_symab


         IF (nb > 0) THEN
           CALL alloc_NParray(P_SymAbelian%tab_symab,(/ nb /),          &
                             'P_SymAbelian%tab_symab',name_sub)
           P_SymAbelian%tab_symab(:) = 0
         END IF

       END SUBROUTINE alloc_SymAbelian

       SUBROUTINE dealloc_SymAbelian(P_SymAbelian)
         TYPE (Type_SymAbelian), pointer, intent(inout) :: P_SymAbelian

         character (len=*), parameter :: name_sub='dealloc_SymAbelian'

         IF (.NOT. associated(P_SymAbelian)) RETURN

         P_SymAbelian%nb          = 0

         IF (allocated(P_SymAbelian%tab_symab)) THEN
           CALL dealloc_NParray(P_SymAbelian%tab_symab,                 &
                               'P_SymAbelian%tab_symab',name_sub)
         END IF
         CALL dealloc_array(P_SymAbelian,'P_SymAbelian',name_sub)

       END SUBROUTINE dealloc_SymAbelian

     SUBROUTINE alloc_array_OF_SymAbeliandim0(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (Type_SymAbelian), pointer, intent(inout) :: tab

      character (len=*), intent(in) :: name_var,name_sub


      integer, parameter :: ndim=0
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_SymAbeliandim0'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       IF (associated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       memory = 1
       allocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Type_SymAbelian')

      END SUBROUTINE alloc_array_OF_SymAbeliandim0
      SUBROUTINE dealloc_array_OF_SymAbeliandim0(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (Type_SymAbelian), pointer, intent(inout) :: tab
      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_SymAbeliandim0'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN

       IF (.NOT. associated(tab))                                       &
                 CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = 1
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Type_SymAbelian')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_SymAbeliandim0

       SUBROUTINE SymAbelian1_TO_SymAbelian2(P_SymAbelian1,P_SymAbelian2)
         IMPLICIT NONE
         TYPE (Type_SymAbelian), pointer, intent(in)    :: P_SymAbelian1
         TYPE (Type_SymAbelian), pointer, intent(inout) :: P_SymAbelian2

         character (len=*), parameter :: name_sub='SymAbelian1_TO_SymAbelian2'

         !write(out_unitp,*) 'BEGINNING ',name_sub

         IF (.NOT. associated(P_SymAbelian1)) THEN
           IF (print_level > 1) THEN
             write(out_unitp,*) 'WARNING in ',name_sub
             write(out_unitp,*) '    SymAbelian1 is not associated!'
           END IF
         ELSE
           IF (.NOT. associated(P_SymAbelian2)) THEN
             CALL alloc_array(P_SymAbelian2,'P_SymAbelian2',name_sub)
           END IF

           IF (SymAbelian_IS_initialized(P_SymAbelian1)) THEN
             CALL alloc_SymAbelian(P_SymAbelian2,P_SymAbelian1%nb)
             P_SymAbelian2%tab_symab(:)   = P_SymAbelian1%tab_symab(:)
           END IF
           P_SymAbelian2%Read_symab     = P_SymAbelian1%Read_symab
           P_SymAbelian2%nb_PER_symab(:) = P_SymAbelian1%nb_PER_symab(:)

         END IF

         !write(out_unitp,*) 'END ',name_sub


       END SUBROUTINE SymAbelian1_TO_SymAbelian2

       SUBROUTINE Write_SymAbelian(P_SymAbelian)
         IMPLICIT NONE
         TYPE (Type_SymAbelian), pointer, intent(in)    :: P_SymAbelian

         integer :: ib

         character (len=*), parameter :: name_sub='Write_SymAbelian'


         IF (.NOT. associated(P_SymAbelian)) THEN
           write(out_unitp,*) 'BEGINNING ',name_sub
           write(out_unitp,*) 'WARNING: "SymAbelian" is not associated!'
           write(out_unitp,*) 'END ',name_sub
           CALL flush_perso(out_unitp)
           RETURN
         END IF

         write(out_unitp,*) 'BEGINNING ',name_sub
         CALL flush_perso(out_unitp)


         write(out_unitp,*) 'nb',P_SymAbelian%nb
         write(out_unitp,*) 'Read_symab',P_SymAbelian%Read_symab
         write(out_unitp,*) 'nb_PER_symab(:)',P_SymAbelian%nb_PER_symab(:)

         write(out_unitp,*) 'alloc tab_symab',allocated(P_SymAbelian%tab_symab)
         IF (allocated(P_SymAbelian%tab_symab)) THEN
           DO ib=1,P_SymAbelian%nb
             write(out_unitp,*) 'ib,tab_symab,bits(tab_symab)',ib,      &
                           WriteTOstring_symab(P_SymAbelian%tab_symab(ib))
           END DO
         END IF

         write(out_unitp,*) 'END ',name_sub
         CALL flush_perso(out_unitp)

       END SUBROUTINE Write_SymAbelian

       FUNCTION SymAbelian_IS_initialized(P_SymAbelian)
         IMPLICIT NONE
         TYPE (Type_SymAbelian), pointer, intent(in)    :: P_SymAbelian
         logical :: SymAbelian_IS_initialized

         IF (.NOT. associated(P_SymAbelian)) THEN
           SymAbelian_IS_initialized = .FALSE.
         ELSE
           SymAbelian_IS_initialized = allocated(P_SymAbelian%tab_symab)
         END IF

       END FUNCTION SymAbelian_IS_initialized

       FUNCTION Get_Read_symabOFSymAbelian(P_SymAbelian)
         IMPLICIT NONE
         TYPE (Type_SymAbelian), pointer, intent(in)    :: P_SymAbelian
         integer :: Get_Read_symabOFSymAbelian


         character (len=*), parameter :: name_sub='Get_Read_symabOFSymAbelian'

         IF (.NOT. associated(P_SymAbelian)) THEN
           Get_Read_symabOFSymAbelian = -2
         ELSE
           Get_Read_symabOFSymAbelian = P_SymAbelian%Read_symab
         END IF

       END FUNCTION Get_Read_symabOFSymAbelian

       FUNCTION Get_symabOFSymAbelian_AT_ib(P_SymAbelian,ib)
         IMPLICIT NONE
         TYPE (Type_SymAbelian), pointer, intent(in)    :: P_SymAbelian
         integer, intent(in) :: ib
         integer :: Get_symabOFSymAbelian_AT_ib


         character (len=*), parameter :: name_sub='Get_symabOFSymAbelian_AT_ib'

         IF (.NOT. associated(P_SymAbelian)) THEN
           Get_symabOFSymAbelian_AT_ib = -1
         ELSE IF (.NOT. SymAbelian_IS_initialized(P_SymAbelian) .OR.    &
                  ib < 1 .OR. ib > P_SymAbelian%nb) THEN
           Get_symabOFSymAbelian_AT_ib = -1
         ELSE
           Get_symabOFSymAbelian_AT_ib = P_SymAbelian%tab_symab(ib)
         END IF


       END FUNCTION Get_symabOFSymAbelian_AT_ib

       SUBROUTINE Set_symabOFSymAbelian_AT_ib(P_SymAbelian,ib,symab)
         IMPLICIT NONE
         TYPE (Type_SymAbelian), pointer, intent(inout)    :: P_SymAbelian
         integer, intent(in) :: ib,symab


         character (len=*), parameter :: name_sub='Set_symabOFSymAbelian_AT_ib'

         IF (.NOT. associated(P_SymAbelian)) THEN
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) '  SymAbelian is not associated !!'
           write(out_unitp,*) '  CHECK the fortran !!'
           STOP
          !ELSE IF (.NOT.  SymAbelian_IS_initialized(P_SymAbelian) .OR.     &
          !        ib < 1 .OR. ib > P_SymAbelian%nb) THEN
         ELSE IF (.NOT.  SymAbelian_IS_initialized(P_SymAbelian) .OR. ib < 1) THEN
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) '  WRONG parameter values'
           write(out_unitp,*) ' ib',ib
           CALL Write_SymAbelian(P_SymAbelian)
           write(out_unitp,*) '  CHECK the fortran !!'
           STOP
         END IF
         IF (ib > P_SymAbelian%nb) RETURN

         IF (P_SymAbelian%Read_symab == -1) THEN
           P_SymAbelian%tab_symab(ib) = -1
         ELSE
           P_SymAbelian%tab_symab(ib) = symab
         END IF

       END SUBROUTINE Set_symabOFSymAbelian_AT_ib

       SUBROUTINE Set_tab_symabOFSymAbelian_WITH_tab(P_SymAbelian,tab)
         IMPLICIT NONE
         TYPE (Type_SymAbelian), pointer, intent(inout)    :: P_SymAbelian
         integer, intent(in) :: tab(:)


         character (len=*), parameter :: name_sub='Set_tab_symabOFSymAbelian_WITH_tab'


         IF ( size(tab) < 1 ) THEN
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) '  the size of "tab" is < 1'
           write(out_unitp,*) '  probably, it is not allocated !'
           write(out_unitp,*) '  CHECK the fortran !!'
           STOP
         END IF

         CALL alloc_SymAbelian(P_SymAbelian,size(tab))

         P_SymAbelian%tab_symab(:) = tab(:)

         CALL Set_nbPERsym_FROM_SymAbelian(P_SymAbelian)

       END SUBROUTINE Set_tab_symabOFSymAbelian_WITH_tab

       SUBROUTINE Set_tab_SymAbelian(P_SymAbelian,nb,symab)
         IMPLICIT NONE
         TYPE (Type_SymAbelian), pointer, intent(inout)    :: P_SymAbelian
         integer, intent(in) :: nb

         integer, intent(in), optional :: symab     ! symab=-1 (nosym)
                                                    ! symab=0 (000), 1 (001), 2 (010), 4 (100)


         integer :: ib,symab_loc

         character (len=*), parameter :: name_sub='Set_tab_SymAbelian'

         IF (.NOT. associated(P_SymAbelian)) THEN
           CALL alloc_array(P_SymAbelian,'P_SymAbelian',name_sub)
         END IF

         IF (present(symab)) THEN
           symab_loc = symab
         ELSE
           symab_loc = P_SymAbelian%Read_symab
         END IF

         IF (nb > 0) THEN
           IF (SymAbelian_IS_initialized(P_SymAbelian)) THEN
             CALL dealloc_SymAbelian(P_SymAbelian)
           END IF
           CALL alloc_SymAbelian(P_SymAbelian,nb)

           CALL Set_ReadsymabOFSymAbelian(P_SymAbelian,symab_loc)


           SELECT CASE (P_SymAbelian%Read_symab)
           CASE (-1)
             P_SymAbelian%tab_symab(:) = -1
           CASE (0,1,2,3,4,5,6,7)
             P_SymAbelian%tab_symab(:) = 0
             DO ib=2,nb,2  ! tab_symab = [0 s 0 s 0 s ....] with s=symab
               P_SymAbelian%tab_symab(ib) = P_SymAbelian%Read_symab
             END DO
           CASE DEFAULT
             write(out_unitp,*) ' ERROR in ',name_sub
             write(out_unitp,*) '  it should never append. The error should come from'
             write(out_unitp,*) ' "Set_ReadsymabOFSymAbelian" subroutine'
             write(out_unitp,*) ' CHECK the fortran!!'
             STOP
           END SELECT
         END IF
         CALL Set_nbPERsym_FROM_SymAbelian(P_SymAbelian)

       END SUBROUTINE Set_tab_SymAbelian

       SUBROUTINE Set_nbPERsym_FROM_SymAbelian(P_SymAbelian)
         IMPLICIT NONE
         TYPE (Type_SymAbelian), pointer, intent(inout)    :: P_SymAbelian

         integer :: isym

         character (len=*), parameter :: name_sub='Set_nbPERsym_FROM_SymAbelian'

         IF (.NOT.  SymAbelian_IS_initialized(P_SymAbelian)) RETURN


         DO isym=-1,ubound(P_SymAbelian%nb_PER_symab,dim=1)
           P_SymAbelian%nb_PER_symab(isym) = count(P_SymAbelian%tab_symab == isym)
         END DO
         IF (P_SymAbelian%nb_PER_symab(-1) == 0) THEN
           P_SymAbelian%nb_PER_symab(-1) = sum(P_SymAbelian%nb_PER_symab)
         END IF

       END SUBROUTINE Set_nbPERsym_FROM_SymAbelian
       integer FUNCTION Get_nbPERsym_FROM_SymAbelian(P_SymAbelian,symab)
         IMPLICIT NONE
         TYPE (Type_SymAbelian), pointer, intent(inout)    :: P_SymAbelian

         integer :: symab

         character (len=*), parameter :: name_sub='Get_nbPERsym_FROM_symab'

         IF (.NOT.  SymAbelian_IS_initialized(P_SymAbelian)) THEN
           Get_nbPERsym_FROM_SymAbelian = 0
         ELSE
           Get_nbPERsym_FROM_SymAbelian = P_SymAbelian%nb_PER_symab(symab)
         END IF

       END FUNCTION Get_nbPERsym_FROM_SymAbelian
       SUBROUTINE Set_ReadsymabOFSymAbelian(P_SymAbelian,Read_symab)
         IMPLICIT NONE
         TYPE (Type_SymAbelian), pointer, intent(inout)    :: P_SymAbelian

         integer, intent(in) :: Read_symab          ! symab=-1 (nosym)
                                                    ! symab=0 (000), 1 (001), 2 (010), 4 (100)


         character (len=*), parameter :: name_sub='Set_ReadsymabOFSymAbelian'

         IF (.NOT. associated(P_SymAbelian)) THEN
           CALL alloc_array(P_SymAbelian,'P_SymAbelian',name_sub)
         END IF

         IF (Read_symab >= -1 .AND.  Read_symab <= 8) THEN

           P_SymAbelian%Read_symab = Read_symab

         ELSE
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) '  WRONG Read_symab',Read_symab
           write(out_unitp,*) '  symab = -1     : no symmetry'
           write(out_unitp,*) '  symab = 0 or 1 : one   symmetry element , Cs, C2...'
           write(out_unitp,*) '  symab = 0-3    : two   symmetry elements, C2v'
           write(out_unitp,*) '  symab = 0-7    : three symmetry elements, D4h'
           write(out_unitp,*) '  CHECK your data'
           STOP
         END IF

       END SUBROUTINE Set_ReadsymabOFSymAbelian

       SUBROUTINE Write_symab(symab)
         IMPLICIT NONE
         integer, intent(in)    :: symab

         character (len=*), parameter :: name_sub='Write_symab'

         write(out_unitp,*) 'symab,bits(symab)',WriteTOstring_symab(symab)

       END SUBROUTINE Write_symab
       FUNCTION WriteTOstring_symab(symab)
         IMPLICIT NONE
         integer, intent(in)  :: symab
         character (len=10)   :: WriteTOstring_symab

         character (len=*), parameter :: name_sub='WriteTOstring_symab'

         IF (symab == -2) THEN
           write(WriteTOstring_symab,11) -2,-2,-2,-2
         ELSE IF (symab < 0 .OR. symab > 7) THEN
           write(WriteTOstring_symab,11) -1,-1,-1,-1
         ELSE
           write(WriteTOstring_symab,11) symab,                         &
                           iand(symab,4)/4,iand(symab,2)/2,iand(symab,1)
 11        format(i3,':',3i2)
         END IF

       END FUNCTION WriteTOstring_symab

       integer FUNCTION Calc_symab1_EOR_symab2(symab1,symab2)
         IMPLICIT NONE
         integer, intent(in)    :: symab1,symab2

         character (len=*), parameter :: name_sub='Calc_symab1_EOR_symab2'

         IF (symab1 == -2) THEN
           Calc_symab1_EOR_symab2 = symab2
         ELSE IF (symab2 == -2) THEN
           Calc_symab1_EOR_symab2 = symab1
         ELSE
           IF (symab1 < 0 .OR. symab1 > 7) THEN
             Calc_symab1_EOR_symab2 = -1
           ELSE
             IF (symab2 < 0 .OR. symab2 > 7) THEN
               Calc_symab1_EOR_symab2 = -1
             ELSE
               Calc_symab1_EOR_symab2 = ieor(symab1,symab2)
             END IF
           END IF
         END IF

       END FUNCTION Calc_symab1_EOR_symab2

      END MODULE mod_SymAbelian
