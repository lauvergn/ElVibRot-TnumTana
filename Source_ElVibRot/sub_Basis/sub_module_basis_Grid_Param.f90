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
      MODULE mod_Basis_Grid_Param
      USE mod_system
      USE mod_nDindex
      IMPLICIT NONE

        TYPE Basis_Grid_Param
          integer :: LGrid_max = -1
          integer :: nq        =  0
          integer :: nq_init   =  0


        END TYPE Basis_Grid_Param

        INTERFACE assignment (=)
          MODULE PROCEDURE Basis_Grid_Param2TOBasis_Grid_Param1
        END INTERFACE

      CONTAINS

      SUBROUTINE Write_Basis_Grid_Param(Basis_Grid_Para,Rec_line)

      TYPE (Basis_Grid_Param), intent(in) :: Basis_Grid_Para
      character (len=*), optional         :: Rec_line

      integer :: L
!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Write_Basis_Grid_Param'
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------

      IF (present(Rec_line)) THEN
       write(out_unitp,*) trim(Rec_line),'LGridmax',Basis_Grid_Para%LGrid_max
       write(out_unitp,*) trim(Rec_line),'nq',Basis_Grid_Para%nq
       write(out_unitp,*) trim(Rec_line),'nq_ini',Basis_Grid_Para%nq_init

      ELSE
       write(out_unitp,*) 'LGridmax',Basis_Grid_Para%LGrid_max
       write(out_unitp,*) 'nq',Basis_Grid_Para%nq
       write(out_unitp,*) 'nq_ini',Basis_Grid_Para%nq_init
      END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------

      END SUBROUTINE Write_Basis_Grid_Param

      SUBROUTINE Basis_Grid_Param2TOBasis_Grid_Param1(Basis_Grid_Para1, &
                                                      Basis_Grid_Para2)


      TYPE (Basis_Grid_Param), intent(inout) :: Basis_Grid_Para1
      TYPE (Basis_Grid_Param), intent(in) :: Basis_Grid_Para2


!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Basis_Grid_Param2TOBasis_Grid_Param1'
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------


      Basis_Grid_Para1%LGrid_max = Basis_Grid_Para2%LGrid_max
      Basis_Grid_Para1%nq        = Basis_Grid_Para2%nq
      Basis_Grid_Para1%nq_init   = Basis_Grid_Para2%nq_init

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        CALL Write_Basis_Grid_Param(Basis_Grid_Para1)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE Basis_Grid_Param2TOBasis_Grid_Param1

      SUBROUTINE Basis_Grid_ParamTOBasis_Grid_Param_init(Basis_Grid_Para)


      TYPE (Basis_Grid_Param), intent(inout) :: Basis_Grid_Para


!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Basis_Grid_ParamTOBasis_Grid_Param_init'
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------


      Basis_Grid_Para%nq_init   = Basis_Grid_Para%nq

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        CALL Write_Basis_Grid_Param(Basis_Grid_Para)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE Basis_Grid_ParamTOBasis_Grid_Param_init

      SUBROUTINE Basis_Grid_Param_initTOBasis_Grid_Param(Basis_Grid_Para)


      TYPE (Basis_Grid_Param), intent(inout) :: Basis_Grid_Para


!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Basis_Grid_Param_initTOBasis_Grid_Param'
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------


      Basis_Grid_Para%nq   = Basis_Grid_Para%nq_init

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        CALL Write_Basis_Grid_Param(Basis_Grid_Para)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE Basis_Grid_Param_initTOBasis_Grid_Param

      END MODULE mod_Basis_Grid_Param
