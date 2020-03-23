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
      MODULE mod_RectilinearNM_Transfo
      use mod_system
      USE mod_dnSVM
      use mod_constant,     only: table_atom, get_mass_tnum
      IMPLICIT NONE

      PRIVATE

      !!@description: TODO
      !!@param: TODO
      TYPE Type_RectilinearNM_Transfo

        integer           :: ncart=0,ncart_act=0
        integer           :: nat0=0,nat=0,nat_act=0
        integer           :: nb_var=0

        ! just for read the input data
        real (kind=Rkind),        pointer :: masses(:)     => null() ! TRUE pointer
        integer, pointer                  :: type_Qin(:)   => null() ! TRUE pointer
        character (len=Name_len), pointer :: name_Qin(:)   => null() ! TRUE pointer
        integer, pointer                  :: Z(:)          => null() ! TRUE pointer
        character (len=Name_len),pointer  :: symbole(:)    => null() ! TRUE pointer


        real (kind=Rkind), pointer  :: mat(:,:)=>null()
        real (kind=Rkind), pointer  :: mat_inv(:,:)=>null()
        logical :: inv=.FALSE.

      END TYPE Type_RectilinearNM_Transfo

      PUBLIC :: Type_RectilinearNM_Transfo, alloc_RectilinearNM_Transfo, dealloc_RectilinearNM_Transfo
      PUBLIC :: Read_RectilinearNM_Transfo, Write_RectilinearNM_Transfo, calc_RectilinearNM_Transfo
      PUBLIC :: RectilinearNM_Transfo1TORectilinearNM_Transfo2

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

      SUBROUTINE alloc_RectilinearNM_Transfo(RectilinearNM_Transfo,nb_Qin)
      TYPE (Type_RectilinearNM_Transfo), intent(inout) :: RectilinearNM_Transfo
      integer, intent(in) :: nb_Qin

      integer :: err
      character (len=*), parameter :: name_sub='alloc_RectilinearNM_Transfo'

!      write(out_unitp,*) 'BEGINNING ',name_sub
!      write(out_unitp,*) 'nat',RectilinearNM_Transfo%nat

       IF (RectilinearNM_Transfo%nat < 3) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' wrong value of nat',RectilinearNM_Transfo%nat
         write(out_unitp,*) ' CHECK the source !!'
         STOP
       END IF

      IF (associated(RectilinearNM_Transfo%mat))  THEN
        CALL dealloc_array(RectilinearNM_Transfo%mat,                   &
                          "RectilinearNM_Transfo%mat",name_sub)
      END IF
      IF (associated(RectilinearNM_Transfo%mat_inv))  THEN
        CALL dealloc_array(RectilinearNM_Transfo%mat_inv,               &
                          "RectilinearNM_Transfo%mat_inv",name_sub)
      END IF

      CALL alloc_array(RectilinearNM_Transfo%mat,(/nb_Qin,nb_Qin/),     &
                      "RectilinearNM_Transfo%mat",name_sub)
      RectilinearNM_Transfo%mat(:,:) = ZERO

      CALL alloc_array(RectilinearNM_Transfo%mat_inv,(/nb_Qin,nb_Qin/), &
                      "RectilinearNM_Transfo%mat_inv",name_sub)
      RectilinearNM_Transfo%mat_inv(:,:) = ZERO


!      write(out_unitp,*) 'END ',name_sub

      END SUBROUTINE alloc_RectilinearNM_Transfo

      SUBROUTINE dealloc_RectilinearNM_Transfo(RectilinearNM_Transfo)

       TYPE (Type_RectilinearNM_Transfo), intent(inout) :: RectilinearNM_Transfo

      character (len=*), parameter :: name_sub='dealloc_RectilinearNM_Transfo'
       !write(out_unitp,*) 'BEGINNING ',name_sub; call flush_perso(out_unitp)

      RectilinearNM_Transfo%ncart     = 0
      RectilinearNM_Transfo%ncart_act = 0
      RectilinearNM_Transfo%nat0      = 0
      RectilinearNM_Transfo%nat       = 0
      RectilinearNM_Transfo%nat_act   = 0
      RectilinearNM_Transfo%nb_var    = 0

      nullify(RectilinearNM_Transfo%masses)
      nullify(RectilinearNM_Transfo%type_Qin)
      nullify(RectilinearNM_Transfo%name_Qin)
      nullify(RectilinearNM_Transfo%Z)
      nullify(RectilinearNM_Transfo%symbole)

      IF (associated(RectilinearNM_Transfo%mat))  THEN
        CALL dealloc_array(RectilinearNM_Transfo%mat,                   &
                          "RectilinearNM_Transfo%mat",name_sub)
      END IF
      IF (associated(RectilinearNM_Transfo%mat_inv))  THEN
        CALL dealloc_array(RectilinearNM_Transfo%mat_inv,               &
                          "RectilinearNM_Transfo%mat_inv",name_sub)
      END IF

      RectilinearNM_Transfo%inv                 = .FALSE.

       !write(out_unitp,*) 'END ',name_sub; call flush_perso(out_unitp)

      END SUBROUTINE dealloc_RectilinearNM_Transfo

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE Read_RectilinearNM_Transfo(RectilinearNM_Transfo,mendeleev)


       TYPE (Type_RectilinearNM_Transfo),intent(inout) :: RectilinearNM_Transfo
       TYPE (table_atom), intent(in)         :: mendeleev


       real (kind=Rkind)        :: at,Mtot
       integer                  :: i,i_at
       character (len=Name_len), pointer :: name_at(:)


       integer :: err_mem,memory,err_io
       character (len=*), parameter :: name_sub = 'Read_RectilinearNM_Transfo'


!-----------------------------------------------------------------------
        RectilinearNM_Transfo%ncart = 3 * RectilinearNM_Transfo%nat
        RectilinearNM_Transfo%nat_act = RectilinearNM_Transfo%nat
        RectilinearNM_Transfo%ncart_act = 3 * RectilinearNM_Transfo%nat_act

        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nat',RectilinearNM_Transfo%nat
        write(out_unitp,*) 'nb_var',RectilinearNM_Transfo%nb_var
        write(out_unitp,*) 'ncart',RectilinearNM_Transfo%ncart


        ! allocation of the variables:
        CALL alloc_RectilinearNM_Transfo(RectilinearNM_Transfo,RectilinearNM_Transfo%ncart)



        CALL alloc_array(name_at,(/RectilinearNM_Transfo%nat/),Name_len,&
                        "name_at",name_sub)

        read(in_unitp,*,IOSTAT=err_io) (name_at(i),i=1,RectilinearNM_Transfo%nat)
        IF (err_io /= 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  while reading a mass in the "RectilinearNM" transformation'
          write(out_unitp,*) ' end of file or end of record'
          write(out_unitp,*) ' Check your data !!'
          STOP
        END IF

        Mtot = ZERO
        DO i_at=1,RectilinearNM_Transfo%nat
          RectilinearNM_Transfo%Z(i_at) = -1
          RectilinearNM_Transfo%symbole(i_at) = name_at(i_at)
          at =get_mass_Tnum(mendeleev,Z=RectilinearNM_Transfo%Z(i_at),name=name_at(i_at))

          RectilinearNM_Transfo%masses((i_at-1)*3+1:(i_at-1)*3+3) = at
          !write(out_unitp,*) 'atom:',i_at,Z(i_at),at
          Mtot = Mtot + at
          IF (at == ZERO) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' The read atom cannot be dummy !'
            write(out_unitp,*) 'atom:',i_at,RectilinearNM_Transfo%Z(i_at),at
            STOP
          END IF

        END DO
        CALL dealloc_array(name_at,"name_at",name_sub)
        !write(out_unitp,*) 'Masses: ',masses(:)


      write(out_unitp,*) 'END ',name_sub
      END SUBROUTINE Read_RectilinearNM_Transfo


!=================================================================
!
!       Rectilinear Normal modes
!
!=================================================================
      SUBROUTINE calc_RectilinearNM_Transfo(dnQin,dnQout,RectilinearNM_Transfo,nderiv,inTOout)

        TYPE (Type_dnVec), intent(inout)      :: dnQin,dnQout
        TYPE (Type_RectilinearNM_Transfo), intent(in) :: RectilinearNM_Transfo

        integer, intent(in)                   :: nderiv
        logical, intent(in)                   :: inTOout


        integer :: i,j,k
        character (len=*), parameter :: name_sub='calc_RectilinearNM_Transfo'



        CALL check_alloc_dnVec(dnQin,'dnQin',name_sub)
        CALL check_alloc_dnVec(dnQout,'dnQout',name_sub)

        IF (inTOout) THEN
          IF (nderiv == 0) THEN
            dnQout%d0 = matmul(RectilinearNM_Transfo%mat,dnQin%d0)
          ELSE IF (nderiv == 1) THEN
            dnQout%d0 = matmul(RectilinearNM_Transfo%mat,dnQin%d0)
            DO i=1,dnQin%nb_var_deriv ! mole%nb_act
              dnQout%d1(:,i) = matmul(RectilinearNM_Transfo%mat,dnQin%d1(:,i))
            END DO
          ELSE IF (nderiv == 2) THEN
            dnQout%d0 = matmul(RectilinearNM_Transfo%mat,dnQin%d0)
            DO i=1,dnQin%nb_var_deriv ! mole%nb_act
              dnQout%d1(:,i) = matmul(RectilinearNM_Transfo%mat,dnQin%d1(:,i))
            END DO
            DO i=1,dnQin%nb_var_deriv ! mole%nb_act
            DO j=1,dnQin%nb_var_deriv ! mole%nb_act
              dnQout%d2(:,i,j) = matmul(RectilinearNM_Transfo%mat,dnQin%d2(:,i,j))
            END DO
            END DO
          ELSE IF (nderiv == 3) THEN
            dnQout%d0 = matmul(RectilinearNM_Transfo%mat,dnQin%d0)
            DO i=1,dnQin%nb_var_deriv ! mole%nb_act
              dnQout%d1(:,i) = matmul(RectilinearNM_Transfo%mat,dnQin%d1(:,i))
            END DO
            DO i=1,dnQin%nb_var_deriv ! mole%nb_act
            DO j=1,dnQin%nb_var_deriv ! mole%nb_act
              dnQout%d2(:,i,j) = matmul(RectilinearNM_Transfo%mat,dnQin%d2(:,i,j))
            END DO
            END DO
            DO i=1,dnQin%nb_var_deriv ! mole%nb_act
            DO j=1,dnQin%nb_var_deriv ! mole%nb_act
            DO k=1,dnQin%nb_var_deriv ! mole%nb_act
              dnQout%d3(:,i,j,k) = matmul(RectilinearNM_Transfo%mat,dnQin%d3(:,i,j,k))
            END DO
            END DO
            END DO
          ELSE
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' nderiv > 4 is NOT possible',nderiv
            write(out_unitp,*) 'It should never append! Check the source'
            STOP
          END IF
        ELSE
          IF (nderiv == 0) THEN
            dnQin%d0 = matmul(RectilinearNM_Transfo%mat_inv,dnQout%d0)
          ELSE IF (nderiv == 1) THEN
            dnQin%d0 = matmul(RectilinearNM_Transfo%mat_inv,dnQout%d0)
            DO i=1,dnQout%nb_var_deriv ! mole%nb_act
              dnQin%d1(:,i) = matmul(RectilinearNM_Transfo%mat_inv,dnQout%d1(:,i))
            END DO
          ELSE IF (nderiv == 2) THEN
            dnQin%d0 = matmul(RectilinearNM_Transfo%mat_inv,dnQout%d0)
            DO i=1,dnQout%nb_var_deriv ! mole%nb_act
              dnQin%d1(:,i) = matmul(RectilinearNM_Transfo%mat_inv,dnQout%d1(:,i))
            END DO
            DO i=1,dnQout%nb_var_deriv ! mole%nb_act
            DO j=1,dnQout%nb_var_deriv ! mole%nb_act
              dnQin%d2(:,i,j) = matmul(RectilinearNM_Transfo%mat_inv,dnQout%d2(:,i,j))
            END DO
            END DO
          ELSE IF (nderiv == 3) THEN
            dnQin%d0 = matmul(RectilinearNM_Transfo%mat_inv,dnQout%d0)
            DO i=1,dnQout%nb_var_deriv ! mole%nb_act
              dnQin%d1(:,i) = matmul(RectilinearNM_Transfo%mat_inv,dnQout%d1(:,i))
            END DO
            DO i=1,dnQout%nb_var_deriv ! mole%nb_act
            DO j=1,dnQout%nb_var_deriv ! mole%nb_act
              dnQin%d2(:,i,j) = matmul(RectilinearNM_Transfo%mat_inv,dnQout%d2(:,i,j))
            END DO
            END DO
            DO i=1,dnQout%nb_var_deriv ! mole%nb_act
            DO j=1,dnQout%nb_var_deriv ! mole%nb_act
            DO k=1,dnQout%nb_var_deriv ! mole%nb_act
              dnQin%d3(:,i,j,k) = matmul(RectilinearNM_Transfo%mat_inv,dnQout%d3(:,i,j,k))
            END DO
            END DO
            END DO
          ELSE
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' nderiv > 4 is NOT possible',nderiv
            write(out_unitp,*) 'It should never append! Check the source'
            STOP
          END IF
        END IF


      END SUBROUTINE calc_RectilinearNM_Transfo



      !!@description: TODO
      !!@param: TODO
      SUBROUTINE RectilinearNM_Transfo1TORectilinearNM_Transfo2(RectilinearNM_Transfo1,RectilinearNM_Transfo2)

      TYPE (Type_RectilinearNM_Transfo), intent(in)    :: RectilinearNM_Transfo1
      TYPE (Type_RectilinearNM_Transfo), intent(inout) :: RectilinearNM_Transfo2

      integer :: it
      character (len=*), parameter ::                                   &
                                name_sub = 'RectilinearNM_Transfo1TORectilinearNM_Transfo2'

      CALL dealloc_RectilinearNM_Transfo(RectilinearNM_Transfo2)

      RectilinearNM_Transfo2%ncart        = RectilinearNM_Transfo1%ncart
      RectilinearNM_Transfo2%ncart_act    = RectilinearNM_Transfo1%ncart_act
      RectilinearNM_Transfo2%nat          = RectilinearNM_Transfo1%nat
      RectilinearNM_Transfo2%nat0         = RectilinearNM_Transfo1%nat0
      RectilinearNM_Transfo2%nat_act      = RectilinearNM_Transfo1%nat_act
      RectilinearNM_Transfo2%nb_var       = RectilinearNM_Transfo1%nb_var


      CALL alloc_RectilinearNM_Transfo(RectilinearNM_Transfo2,RectilinearNM_Transfo2%ncart)

      RectilinearNM_Transfo2%mat       = RectilinearNM_Transfo1%mat
      RectilinearNM_Transfo2%mat_inv   = RectilinearNM_Transfo1%mat_inv

      RectilinearNM_Transfo2%inv       = RectilinearNM_Transfo1%inv

!     write(out_unitp,*) 'END RectilinearNM_Transfo1TORectilinearNM_Transfo2'

      END SUBROUTINE RectilinearNM_Transfo1TORectilinearNM_Transfo2

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE Write_RectilinearNM_Transfo(RectilinearNM_Transfo)
      TYPE (Type_RectilinearNM_Transfo), intent(in) :: RectilinearNM_Transfo

      integer :: i
      character (len=*), parameter :: name_sub='Write_RectilinearNM_Transfo'


      write(out_unitp,*) 'BEGINNING ',name_sub

      write(out_unitp,*) 'ncart_act,ncart',                             &
                  RectilinearNM_Transfo%ncart_act,RectilinearNM_Transfo%ncart

      write(out_unitp,*) 'nat_act,nat0,nat,',                           &
                  RectilinearNM_Transfo%nat_act,RectilinearNM_Transfo%nat0,RectilinearNM_Transfo%nat

      write(out_unitp,*) 'nb_var',RectilinearNM_Transfo%nb_var

      write(out_unitp,*) 'inv',RectilinearNM_Transfo%inv
      write(out_unitp,*)  'Mat of RectilinearNM_Transfo: '
      CALL Write_Mat(RectilinearNM_Transfo%mat,out_unitp,4)

      write(out_unitp,*)  'Mat_inv of RectilinearNM_Transfo: '
      CALL Write_Mat(RectilinearNM_Transfo%mat_inv,out_unitp,4)

      write(out_unitp,*) 'END ',name_sub

      END SUBROUTINE Write_RectilinearNM_Transfo

      END MODULE mod_RectilinearNM_Transfo

