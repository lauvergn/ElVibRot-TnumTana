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
  USE mod_system,                ONLY : Rkind,out_unitp,print_level
  USE mod_Constant,              ONLY : constant
  USE mod_Coord_KEO,             ONLY : CoordType,Tnum,Read_CoordType,  &
                                        read_RefGeom,sub_QactTOd0x,sub_d0xTOQact
  USE mod_PrimOp,                ONLY : param_PES,Finalyze_TnumTana_Coord_PrimOp
  IMPLICIT NONE

  TYPE (constant)  :: const_phys
  TYPE (CoordType) :: mole
  TYPE (Tnum)      :: para_Tnum
  TYPE (param_PES) :: para_PES

  integer          :: Init = 0 ! Initialization is not done

END MODULE Module_ForTnumTana_Driver
