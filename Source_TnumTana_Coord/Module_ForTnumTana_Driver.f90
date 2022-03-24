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
MODULE Module_ForTnumTana_Driver
  USE mod_system
  USE mod_Constant
  USE mod_Coord_KEO,             ONLY: CoordType,Tnum,Read_CoordType,           &
                                       read_RefGeom,sub_QactTOd0x,sub_d0xTOQact,&
                                       sub_QactTOdnx
  USE mod_PrimOp,                ONLY: PrimOp_t,Finalize_TnumTana_Coord_PrimOp
  IMPLICIT NONE

  TYPE (constant)  :: const_phys
  TYPE (CoordType) :: mole
  TYPE (Tnum)      :: para_Tnum
  TYPE (PrimOp_t)  :: PrimOp

  integer          :: Init    =  0  ! Initialization is not done
  integer          :: skip_NM =  0  ! if 0 => calc the NM, if 1 no NM calculation
  integer          :: k_Half  = -1  ! if -1, k_Half is not modified.
                                    ! 1 => k_Half is set to T, 0 => k_Half is set to F

CONTAINS
SUBROUTINE Check_TnumInit(name_sub)
  IMPLICIT NONE

  character(len=*), intent(in) :: name_sub


  IF (Init == 0) THEN
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) ' Tnum is not initialized!'
    STOP 'ERROR:  Tnum is not initialized!'
  END IF
END SUBROUTINE Check_TnumInit
END MODULE Module_ForTnumTana_Driver
