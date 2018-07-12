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
      MODULE mod_ReadOp

      USE mod_system
      USE mod_OpGrid
      IMPLICIT NONE

!       ====================================================================
!
!       ====================================================================
        TYPE param_ReadOp ! used for transfert info from read_active to para_H

           logical               :: OpPsi_WithGrid = .FALSE.

           logical               :: pack_Op        = .FALSE.
           real (kind=Rkind)     :: tol_pack       = ONETENTH**7
           real (kind=Rkind)     :: tol_nopack     = NINE*ONETENTH
           logical               :: read_Op        = .FALSE.
           logical               :: make_Mat       = .FALSE.
           logical               :: spectral       = .FALSE.       ! IF T, spectral represention
           integer               :: spectral_Op    = 0             ! IF sepctral, we use the H (nOp=0)

           logical               :: pot_only       = .FALSE.
           logical               :: T_only         = .FALSE.
           integer               :: nb_bRot        = 0

           TYPE (param_FileGrid) :: para_FileGrid               ! parameters to tranfer to OpGrid%...

           logical               :: comput_S    = .FALSE.       ! calculation of the active overlap matrix

           logical           :: Op_Transfo = .FALSE.   ! true => we are using Transfo(Op) instead of Op
           real (kind=Rkind) :: E0_Transfo = ZERO      ! ScaledOp = Op - E0_transfo * I
           integer           :: degree_Transfo = -1     ! degree of the transformation
           real (kind=Rkind), allocatable :: Poly_Transfo(:) ! Poly_transfo(0:degree_Transfo)
                                               ! Transfo(Op) = Sum_k Poly_transfo(k) * ScaledOp^k



        END TYPE param_ReadOp

      INTERFACE assignment (=)
        MODULE PROCEDURE ReadOp2_TO_ReadOp1
      END INTERFACE

      CONTAINS

      SUBROUTINE init_ReadOp(para_ReadOp)
      TYPE (param_ReadOp) :: para_ReadOp

      para_ReadOp%pack_Op     = .FALSE.
      para_ReadOp%tol_pack    = ONETENTH**7
      para_ReadOp%tol_nopack  = NINE*ONETENTH
      para_ReadOp%read_Op     = .FALSE.
      para_ReadOp%make_Mat    = .FALSE.
      para_ReadOp%spectral    = .FALSE.
      para_ReadOp%spectral_Op = 0
      para_ReadOp%nb_bRot     = 0


      para_ReadOp%pot_only    = .FALSE.
      para_ReadOp%T_only      = .FALSE.
      para_ReadOp%comput_S    = .FALSE.   ! calculation of the active overlap matrix

      CALL init_FileGrid(para_ReadOp%para_FileGrid)

      para_ReadOp%Op_Transfo = .FALSE.   ! true => we are using Transfo(Op) instead of Op
      para_ReadOp%E0_Transfo = ZERO      ! ScaledOp = Op - E0_transfo * I
      para_ReadOp%degree_Transfo = -1     ! degree of the transformation

      END SUBROUTINE init_ReadOp

      SUBROUTINE ReadOp2_TO_ReadOp1(para_ReadOp1,para_ReadOp2)
      TYPE (param_ReadOp), intent(inout) :: para_ReadOp1
      TYPE (param_ReadOp), intent(in)    :: para_ReadOp2


      para_ReadOp1%OpPsi_WithGrid  = para_ReadOp2%OpPsi_WithGrid

      para_ReadOp1%pack_Op         = para_ReadOp2%pack_Op
      para_ReadOp1%tol_pack        = para_ReadOp2%tol_pack
      para_ReadOp1%tol_nopack      = para_ReadOp2%tol_nopack
      para_ReadOp1%read_Op         = para_ReadOp2%read_Op
      para_ReadOp1%make_Mat        = para_ReadOp2%make_Mat
      para_ReadOp1%spectral        = para_ReadOp2%spectral
      para_ReadOp1%spectral_Op     = para_ReadOp2%spectral_Op
      para_ReadOp1%nb_bRot         = para_ReadOp2%nb_bRot

      para_ReadOp1%pot_only        = para_ReadOp2%pot_only
      para_ReadOp1%T_only          = para_ReadOp2%T_only
      para_ReadOp1%comput_S        = para_ReadOp2%comput_S

      para_ReadOp1%para_FileGrid   = para_ReadOp2%para_FileGrid

      para_ReadOp1%Op_Transfo     = para_ReadOp2%Op_Transfo
      para_ReadOp1%E0_Transfo     = para_ReadOp2%E0_Transfo
      para_ReadOp1%degree_Transfo = para_ReadOp2%degree_Transfo
      IF (allocated(para_ReadOp2%Poly_Transfo)) THEN
        para_ReadOp1%Poly_Transfo = para_ReadOp2%Poly_Transfo
      END IF

      END SUBROUTINE ReadOp2_TO_ReadOp1


      END MODULE mod_ReadOp

