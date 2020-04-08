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

      MODULE mod_cart
      USE mod_system
      USE mod_dnSVM
      IMPLICIT NONE

      !!@description: TODO
      !!@param: TODO
        TYPE Type_cart

          integer                    :: nb_at     = 0
          integer                    :: nb_vect   = 0

          TYPE (Type_dnVec), pointer :: dnAt(:)   => null() ! table of vectors (ndim=3)
                                                            ! for cartesian coordinates of the atoms
          TYPE (Type_dnVec), pointer :: dnVect(:) => null() ! table of vectors (ndim=3)
                                                            ! for cartesian coordinates of the vectors
          real(kind=Rkind), pointer :: masses(:)  => null() ! masses(nb_at)
          real(kind=Rkind)          :: Mtot       = ZERO

        END TYPE Type_cart
      CONTAINS

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE alloc_Type_cart(para_cart,nb_at,nb_vect)
        TYPE (Type_cart), intent(inout) :: para_cart
        integer, optional :: nb_at,nb_vect

        integer :: i
        integer :: err_mem,memory

        IF (present(nb_at)) THEN
          IF (nb_at < 1) THEN
             write(out_unitp,*) ' ERROR in alloc_Type_cart'
             write(out_unitp,*) ' nb_at is present and < 1 !!'
             STOP
          END IF
          para_cart%nb_at = nb_at
          CALL alloc_array(para_cart%masses,(/nb_at/),                  &
                          "para_cart%masses","alloc_Type_cart")
          para_cart%masses(:) = ZERO

          CALL alloc_array(para_cart%dnAt,(/nb_at/),                    &
                          "para_cart%dnAt","alloc_Type_cart")
          DO i=1,nb_at
            CALL alloc_dnSVM(para_cart%dnAt(i),nb_var_vec=3,nderiv=0)
          END DO
        END IF

        IF (present(nb_vect)) THEN
          IF (nb_vect < 1) THEN
             write(out_unitp,*) ' ERROR in alloc_Type_cart'
             write(out_unitp,*) ' nb_vect is present and < 1 !!'
             STOP
          END IF
          para_cart%nb_vect = nb_vect

          CALL alloc_array(para_cart%dnVect,(/nb_vect/),                &
                          "para_cart%dnVect","alloc_Type_cart")
          DO i=1,nb_vect
            CALL alloc_dnSVM(para_cart%dnVect(i),nb_var_vec=3,nderiv=0)
          END DO
        END IF


      END SUBROUTINE alloc_Type_cart

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE dealloc_Type_cart(para_cart)
        TYPE (Type_cart), intent(inout) :: para_cart

        integer :: i
        integer :: err_mem,memory


        IF (associated(para_cart%masses))  THEN
          CALL dealloc_array(para_cart%masses,                          &
                            "para_cart%masses","dealloc_Type_cart")
        END IF

        IF (associated(para_cart%dnAt)) THEN
          DO i=1,para_cart%nb_at
            CALL dealloc_dnSVM(para_cart%dnAt(i))
          END DO
          CALL dealloc_array(para_cart%dnAt,                            &
                            "para_cart%dnAt","dealloc_Type_cart")
        END IF

        IF (associated(para_cart%dnVect)) THEN
          DO i=1,para_cart%nb_vect
            CALL dealloc_dnSVM(para_cart%dnVect(i))
          END DO
          CALL dealloc_array(para_cart%dnVect,                          &
                            "para_cart%dnVect","dealloc_Type_cart")
        END IF

        para_cart%nb_at     = 0
        para_cart%nb_vect   = 0
        para_cart%Mtot      = ZERO


      END SUBROUTINE dealloc_Type_cart

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE write_Type_cart(para_cart)
      TYPE (Type_cart), intent(in) :: para_cart

      integer           :: i

      write(out_unitp,*) 'WRITE Type_cart'
      write(out_unitp,*)
      IF (para_cart%nb_at > 0) THEN
        write(out_unitp,*) 'Cartesian coordinates of the atoms:'
        DO i=1,para_cart%nb_at
          write(out_unitp,*) 'Atom,',i,para_cart%masses(i),para_cart%dnAt(i)%d0(:)
        END DO
        write(out_unitp,*) 'Mtot: ',para_cart%Mtot
      END IF

      IF (para_cart%nb_vect > 0) THEN
        write(out_unitp,*) 'Cartesian coordinates of the vectors:'
        DO i=1,para_cart%nb_vect
          write(out_unitp,*) 'Vector,',i,para_cart%dnVect(i)%d0(:)
        END DO
      END IF

      write(out_unitp,*)
      write(out_unitp,*) 'END WRITE Type_cart'


      END SUBROUTINE write_Type_cart

      END MODULE mod_cart

