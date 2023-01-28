!===========================================================================
!===========================================================================
!This file is part of ElVibRot.
!
!    ElVibRot is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ElVibRot is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ElVibRot.  If not, see <http://www.gnu.org/licenses/>.
!
!    Copyright 2015 David Lauvergnat [1]
!      with contributions of
!        Josep Maria Luis (optimization) [2]
!        Ahai Chen (MPI) [1,4]
!        Lucien Dupuy (CRP) [5]
!
![1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France
![2]: Institut de Química Computacional and Departament de Química,
!        Universitat de Girona, Catalonia, Spain
![3]: Department of Chemistry, Aarhus University, DK-8000 Aarhus C, Denmark
![4]: Maison de la Simulation USR 3441, CEA Saclay, France
![5]: Laboratoire Univers et Particule de Montpellier, UMR 5299,
!         Université de Montpellier, France
!
!    ElVibRot includes:
!        - Tnum-Tana under the GNU LGPL3 license
!        - Somme subroutines of John Burkardt under GNU LGPL license
!             http://people.sc.fsu.edu/~jburkardt/
!        - Somme subroutines of SHTOOLS written by Mark A. Wieczorek under BSD license
!             http://shtools.ipgp.fr
!        - Some subroutine of QMRPack (see cpyrit.doc) Roland W. Freund and Noel M. Nachtigal:
!             https://www.netlib.org/linalg/qmr/
!
!===========================================================================
!===========================================================================
      MODULE mod_ReadOp
      USE mod_system
      USE mod_OpGrid
      USE mod_PrimOp
      IMPLICIT NONE

      PRIVATE

        TYPE, EXTENDS(PrimOp_t) :: param_ReadOp ! used for transfert info from read_active to para_Op

           logical                  :: OpPsi_WithGrid      = .FALSE.

           logical                  :: pack_Op             = .FALSE.
           real (kind=Rkind)        :: tol_pack            = ONETENTH**7
           real (kind=Rkind)        :: tol_nopack          = NINE*ONETENTH
           logical                  :: read_Op             = .FALSE.

           logical                  :: make_Mat            = .FALSE.
           logical                  :: save_MatOp          = .FALSE.
           logical                  :: restart_MatOp       = .FALSE.
           integer                  :: Partial_MatOp_i     = 0
           integer                  :: Partial_MatOp_f     = huge(1)
           logical                  :: formatted_Mat       = .TRUE.
           character (len=Line_len) :: name_Mat            = 'MatOp'

           logical               :: spectral            = .FALSE.       ! IF T, spectral represention
           integer               :: spectral_Op         = 0             ! IF sepctral, we use the H (nOp=0)
           logical               :: Op_WithContracRVec  = .FALSE.
           logical               :: pot_only            = .FALSE.
           logical               :: T_only              = .FALSE.
           integer               :: nb_bRot             = 0

           TYPE (File_tGrid) :: para_FileGrid                       ! parameters to tranfer to OpGrid%...
           TYPE (File_t)     :: FileMat                             ! file Operator Matrix

           logical               :: comput_S            = .FALSE.       ! calculation of the active overlap matrix

           logical               :: Op_Transfo          = .FALSE.       ! true => we are using Transfo(Op) instead of Op
           real (kind=Rkind)     :: E0_Transfo          = ZERO          ! ScaledOp = Op - E0_transfo * I
           integer               :: degree_Transfo      = -1            ! degree of the transformation
           real (kind=Rkind), allocatable :: Poly_Transfo(:) ! Poly_transfo(0:degree_Transfo)
                                               ! Transfo(Op) = Sum_k Poly_transfo(k) * ScaledOp^k

        CONTAINS
          PROCEDURE, PRIVATE, PASS(para_ReadOp1) :: ReadOp2_TO_ReadOp1
          GENERIC,   PUBLIC  :: assignment(=) => ReadOp2_TO_ReadOp1
        END TYPE param_ReadOp

        PUBLIC :: param_ReadOp, init_ReadOp, dealloc_ReadOp, ReadOp2_TO_ReadOp1_FOR_AutoBasis

      CONTAINS

      SUBROUTINE init_ReadOp(para_ReadOp)
      TYPE (param_ReadOp) :: para_ReadOp

      para_ReadOp%pack_Op             = .FALSE.
      para_ReadOp%tol_pack            = ONETENTH**7
      para_ReadOp%tol_nopack          = NINE*ONETENTH
      para_ReadOp%read_Op             = .FALSE.
      para_ReadOp%make_Mat            = .FALSE.
      para_ReadOp%save_MatOp          = .FALSE.
      para_ReadOp%restart_MatOp       = .FALSE.
      para_ReadOp%Partial_MatOp_i     = 0
      para_ReadOp%Partial_MatOp_f     = huge(1)
      para_ReadOp%formatted_Mat       = .TRUE.
      para_ReadOp%name_Mat            = 'MatOp'

      para_ReadOp%spectral            = .FALSE.
      para_ReadOp%Op_WithContracRVec  = .FALSE.

      para_ReadOp%spectral_Op         = 0
      para_ReadOp%nb_bRot             = 0

      para_ReadOp%pot_only            = .FALSE.
      para_ReadOp%T_only              = .FALSE.
      para_ReadOp%comput_S            = .FALSE.   ! calculation of the active overlap matrix

      CALL init_FileGrid(para_ReadOp%para_FileGrid)
      CALL file_dealloc(para_ReadOp%FileMat)

      para_ReadOp%Op_Transfo          = .FALSE.   ! true => we are using Transfo(Op) instead of Op
      para_ReadOp%E0_Transfo          = ZERO      ! ScaledOp = Op - E0_transfo * I
      para_ReadOp%degree_Transfo      = -1     ! degree of the transformation

      END SUBROUTINE init_ReadOp

      SUBROUTINE ReadOp2_TO_ReadOp1(para_ReadOp1,para_ReadOp2)
      CLASS (param_ReadOp), intent(inout) :: para_ReadOp1
      TYPE (param_ReadOp),  intent(in)    :: para_ReadOp2

      para_ReadOp1%PrimOp_t           = para_ReadOp2%PrimOp_t

      para_ReadOp1%OpPsi_WithGrid     = para_ReadOp2%OpPsi_WithGrid

      para_ReadOp1%pack_Op            = para_ReadOp2%pack_Op
      para_ReadOp1%tol_pack           = para_ReadOp2%tol_pack
      para_ReadOp1%tol_nopack         = para_ReadOp2%tol_nopack
      para_ReadOp1%read_Op            = para_ReadOp2%read_Op
      para_ReadOp1%make_Mat           = para_ReadOp2%make_Mat
      para_ReadOp1%save_MatOp         = para_ReadOp2%save_MatOp
      para_ReadOp1%restart_MatOp      = para_ReadOp2%restart_MatOp
      para_ReadOp1%formatted_Mat      = para_ReadOp2%formatted_Mat
      para_ReadOp1%Partial_MatOp_i    = para_ReadOp2%Partial_MatOp_i
      para_ReadOp1%Partial_MatOp_f    = para_ReadOp2%Partial_MatOp_f
      para_ReadOp1%name_Mat           = para_ReadOp2%name_Mat
      para_ReadOp1%FileMat            = para_ReadOp2%FileMat

      para_ReadOp1%spectral           = para_ReadOp2%spectral
      para_ReadOp1%spectral_Op        = para_ReadOp2%spectral_Op
      para_ReadOp1%Op_WithContracRVec = para_ReadOp2%Op_WithContracRVec

      para_ReadOp1%nb_bRot            = para_ReadOp2%nb_bRot

      para_ReadOp1%pot_only           = para_ReadOp2%pot_only
      para_ReadOp1%T_only             = para_ReadOp2%T_only
      para_ReadOp1%comput_S           = para_ReadOp2%comput_S

      para_ReadOp1%para_FileGrid      = para_ReadOp2%para_FileGrid

      para_ReadOp1%Op_Transfo         = para_ReadOp2%Op_Transfo
      para_ReadOp1%E0_Transfo         = para_ReadOp2%E0_Transfo
      para_ReadOp1%degree_Transfo     = para_ReadOp2%degree_Transfo
      IF (allocated(para_ReadOp2%Poly_Transfo)) THEN
        para_ReadOp1%Poly_Transfo     = para_ReadOp2%Poly_Transfo
      END IF

      END SUBROUTINE ReadOp2_TO_ReadOp1

      SUBROUTINE ReadOp2_TO_ReadOp1_FOR_AutoBasis(para_ReadOp1,para_ReadOp2,nq)
      TYPE (param_ReadOp),  intent(inout) :: para_ReadOp1
      TYPE (param_ReadOp),  intent(in)    :: para_ReadOp2
      integer,              intent(in)    :: nq

      para_ReadOp1%PrimOp_t           = para_ReadOp2%PrimOp_t
      para_ReadOp1%nb_scalar_Op       = 0
      para_ReadOp1%nb_CAP             = 0
      para_ReadOp1%nb_FluxOp          = 0
      para_ReadOp1%calc_scalar_Op     = .FALSE.
      para_ReadOp1%type_HamilOp       = 1
      para_ReadOp1%direct_KEO         = .FALSE.


      para_ReadOp1%OpPsi_WithGrid     = .FALSE.

      para_ReadOp1%pack_Op            = .FALSE.
      para_ReadOp1%read_Op            = .FALSE.
      para_ReadOp1%make_Mat           = .TRUE.
      para_ReadOp1%save_MatOp         = .FALSE.
      para_ReadOp1%restart_MatOp      = .FALSE.
      para_ReadOp1%formatted_Mat      = .TRUE.
      para_ReadOp1%name_Mat           = 'MatOp'
      para_ReadOp1%Partial_MatOp_i    = 0
      para_ReadOp1%Partial_MatOp_f    = huge(1)
      CALL file_dealloc(para_ReadOp1%FileMat)

      para_ReadOp1%spectral           = .FALSE.
      para_ReadOp1%Op_WithContracRVec = .FALSE.

      para_ReadOp1%nb_bRot            = 1

      para_ReadOp1%pot_only           = .FALSE.
      para_ReadOp1%T_only             = .FALSE.
      para_ReadOp1%comput_S           = .FALSE.

      para_ReadOp1%para_FileGrid      = para_ReadOp2%para_FileGrid
      para_ReadOp1%para_FileGrid%Save_FileGrid   = .FALSE.
      para_ReadOp1%para_FileGrid%First_GridPoint = 1
      para_ReadOp1%para_FileGrid%Last_GridPoint  = nq
      para_ReadOp1%para_FileGrid%Restart_Grid    = .FALSE.
      para_ReadOp1%para_FileGrid%Test_Grid       = .FALSE.
      para_ReadOp1%para_FileGrid%Read_FileGrid   = .FALSE.
      para_ReadOp1%para_FileGrid%Type_FileGrid   = 0

      para_ReadOp1%Op_Transfo         = .FALSE.
      para_ReadOp1%E0_Transfo         = ZERO
      para_ReadOp1%degree_Transfo     = -1

    END SUBROUTINE ReadOp2_TO_ReadOp1_FOR_AutoBasis

      SUBROUTINE dealloc_ReadOp(para_ReadOp)
      CLASS (param_ReadOp), intent(inout) :: para_ReadOp


      CALL dealloc_PrimOp(para_ReadOp%PrimOp_t)

      para_ReadOp%OpPsi_WithGrid      =  .FALSE.
      para_ReadOp%pack_Op             =  .FALSE.
      para_ReadOp%tol_pack            = ONETENTH**7
      para_ReadOp%tol_nopack          = 0.9_Rkind
      para_ReadOp%read_Op             = .FALSE.
      para_ReadOp%make_Mat            = .FALSE.
      para_ReadOp%save_MatOp          = .FALSE.
      para_ReadOp%restart_MatOp       = .FALSE.
      para_ReadOp%formatted_Mat       = .TRUE.
      para_ReadOp%Partial_MatOp_i     = 0
      para_ReadOp%Partial_MatOp_f     = huge(1)
      para_ReadOp%name_Mat            = 'MatOp'
      CALL file_dealloc(para_ReadOp%FileMat)
      para_ReadOp%spectral            = .FALSE.
      para_ReadOp%spectral_Op         = 0
      para_ReadOp%Op_WithContracRVec  = .FALSE.

      para_ReadOp%nb_bRot             = 0

      para_ReadOp%pot_only            = .FALSE.
      para_ReadOp%T_only              = .FALSE.
      para_ReadOp%comput_S            = .FALSE.

      CALL dealloc_FileGrid(para_ReadOp%para_FileGrid)

      para_ReadOp%Op_Transfo          = .FALSE.
      para_ReadOp%E0_Transfo          = ZERO
      para_ReadOp%degree_Transfo      = -1
      IF (allocated(para_ReadOp%Poly_Transfo)) deallocate(para_ReadOp%Poly_Transfo)

      END SUBROUTINE dealloc_ReadOp

      END MODULE mod_ReadOp
