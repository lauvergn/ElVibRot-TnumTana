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
MODULE mod_dnSVM
      USE mod_dnS
      USE mod_VecOFdnS
      USE mod_MatOFdnS
      USE mod_dnV
      USE mod_dnM
      USE mod_IntVM

      IMPLICIT NONE

      INTERFACE Write_dnSVM
        MODULE PROCEDURE Write_dnS,Write_dnVec,Write_dnMat,Write_IntVec,Write_dnCplxMat
      END INTERFACE
      INTERFACE alloc_dnSVM
        MODULE PROCEDURE alloc_dnS,alloc_dnVec,alloc_dnMat,alloc_IntVec,alloc_dnCplxMat
      END INTERFACE
      INTERFACE dealloc_dnSVM
        MODULE PROCEDURE dealloc_dnS,dealloc_dnVec,dealloc_dnMat,       &
                         dealloc_IntVec,dealloc_dnCplxMat
      END INTERFACE
      INTERFACE Set_ZERO_TO_dnSVM
        MODULE PROCEDURE sub_ZERO_TO_dnS,sub_ZERO_TO_dnVec,             &
                         sub_ZERO_TO_dnMat,sub_ZERO_TO_IntVec,sub_ZERO_TO_dnCplxMat
      END INTERFACE

END MODULE mod_dnSVM

