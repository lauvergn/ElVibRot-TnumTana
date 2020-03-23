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
MODULE mod_Coord_KEO


  USE mod_Constant

  USE mod_Lib_QTransfo,    ONLY : Write_dnx
  USE mod_freq,            ONLY : gaussian_width
  USE mod_ActiveTransfo
  USE mod_RPHTransfo
  USE mod_CartesianTransfo
  USE mod_export_KEO  ! no only
  USE mod_Tnum
  USE mod_paramQ
  USE mod_dnRho
  USE mod_dnGG_dng
  USE mod_dnDetGG_dnDetg
  USE mod_f2f2Vep

  USE mod_Tana_keo    ! no only
  USE mod_Tana_Tnum   ! no only
  USE mod_Tana_Sum_OpnD
  IMPLICIT NONE

END MODULE mod_Coord_KEO
