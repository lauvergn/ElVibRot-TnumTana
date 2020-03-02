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
MODULE mod_EVR
 USE mod_system
!$ USE omp_lib, only : omp_get_max_threads
 USE mod_Constant
 USE mod_Coord_KEO
 USE mod_PrimOp
 USE mod_basis

 USE mod_psi_set_alloc
 USE mod_psi_Op
 USE mod_ana_psi
 USE mod_psi_SimpleOp

 USE mod_propa
 USE mod_Op
 USE mod_analysis
 USE mod_MPI
 IMPLICIT NONE

 TYPE param_EVRT

   !----- physical and mathematical constants ---------------------------
   TYPE (constant) :: const_phys

   !----- On the fly parameters (at this time for gaussian) -------------
   TYPE (param_OTF) :: para_OTF

   !----- for the zmatrix and Tnum --------------------------------------
   TYPE (zmatrix) :: mole
   TYPE (Tnum)    :: para_Tnum

   !----- for the basis set ----------------------------------------------
   TYPE (param_AllBasis) :: para_AllBasis

   !----- variables pour la namelist minimum ----------------------------
   TYPE (param_PES) :: para_PES

   !----- variables for the construction of H ----------------------------
   TYPE (param_ComOp)  :: ComOp
   TYPE (param_AllOp)  :: para_AllOp


   !----- variables pour la namelist analyse ----------------------------
   TYPE (param_ana)           :: para_ana
   TYPE (param_intensity)     :: para_intensity
   logical                    :: intensity_only
   integer                    :: nio_res_int

   !----- variables for the WP propagation ----------------------------
   TYPE (param_propa) :: para_propa
   TYPE (param_psi)   :: WP0
 END TYPE param_EVRT

 TYPE (param_EVRT)              :: para_EVRT
 TYPE (param_EVRT), allocatable :: tab_EVRT(:) ! for the openmp


END MODULE mod_EVR

