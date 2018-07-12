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
MODULE mod_param_SGType2
USE mod_system
USE mod_dnSVM
USE mod_nDindex
IMPLICIT NONE

  TYPE param_SGType2
    integer                      :: L1_SparseGrid = huge(1)
    integer                      :: L2_SparseGrid = huge(1)
    integer                      :: Num_OF_Lmax   = 0 ! use normal L_SparseGrid

    integer                      :: nb_SG = 0
    TYPE (Type_nDindex)          :: nDind_SmolyakGrids  ! multidimensional index smolyak grids

    TYPE (Type_nDindex), allocatable :: nDind_DPG(:)    ! multidimensional DP index (nb_SG)
    TYPE (Type_nDindex), allocatable :: nDind_DPB(:)    ! multidimensional DP index (nb_SG)

    !for SGtype=4
    integer, allocatable :: tab_iB_OF_SRep_TO_iB(:)   ! size (SRep)

    integer, allocatable :: tab_Sum_nq_OF_SRep(:)     ! size (SRep)
    integer, allocatable :: tab_nq_OF_SRep(:)         ! size (SRep)

    integer, allocatable :: tab_Sum_nb_OF_SRep(:)     ! size (SRep)
    integer, allocatable :: tab_nb_OF_SRep(:)         ! size (SRep)

  END TYPE param_SGType2

  TYPE OldParam
    integer                      :: i_SG = 0

    integer                      :: iq   = 0
    integer                      :: iq_SG = 0

    integer                      :: ib   = 0
    integer                      :: ib_SG = 0

    integer, allocatable         :: tab_l_AT_SG(:) ! associated to i_SG
  END TYPE OldParam

INTERFACE assignment (=)
  MODULE PROCEDURE SGType2_2TOSGType2_1
END INTERFACE

CONTAINS
SUBROUTINE Write_OldParam(OldPara)

TYPE (OldParam), intent(out) :: OldPara

character (len=*), parameter :: name_sub='Write_OldParam'

  write(out_unitp,*) 'BEGINNING ',name_sub
  write(out_unitp,*) 'iq,ib          ',OldPara%iq,OldPara%ib
  write(out_unitp,*) 'i_SG           ',OldPara%i_SG
  write(out_unitp,*) 'iq_SG,ib_SG    ',OldPara%iq_SG,OldPara%ib_SG

  IF (allocated(OldPara%tab_l_AT_SG)) &
  write(out_unitp,*) 'tab_l_AT_SG(:) ',OldPara%tab_l_AT_SG(:)

  write(out_unitp,*) 'END ',name_sub
  CALL flush_perso(out_unitp)

END SUBROUTINE Write_OldParam


SUBROUTINE dealloc_SGType2(SGType2)

TYPE (param_SGType2), intent(inout) :: SGType2

character (len=*), parameter :: name_sub='dealloc_SGType2'

SGType2%L1_SparseGrid = huge(1)
SGType2%L2_SparseGrid = huge(1)
SGType2%Num_OF_Lmax   = 0 ! use normal L_SparseGrid

SGType2%nb_SG = 0

CALL dealloc_nDindex(SGType2%nDind_SmolyakGrids)

IF (allocated(SGType2%nDind_DPG)) THEN
  CALL dealloc_NParray(SGType2%nDind_DPG,'SGType2%nDind_DPG',name_sub)
END IF

IF (allocated(SGType2%nDind_DPB)) THEN
  CALL dealloc_NParray(SGType2%nDind_DPB,'SGType2%nDind_DPB',name_sub)
END IF

IF (allocated(SGType2%tab_iB_OF_SRep_TO_iB)) THEN
  CALL dealloc_NParray(SGType2%tab_iB_OF_SRep_TO_iB,        &
                      'SGType2%tab_iB_OF_SRep_TO_iB',name_sub)
END IF

IF (allocated(SGType2%tab_Sum_nq_OF_SRep)) THEN
  CALL dealloc_NParray(SGType2%tab_Sum_nq_OF_SRep,        &
                      'SGType2%tab_Sum_nq_OF_SRep',name_sub)
END IF
IF (allocated(SGType2%tab_nq_OF_SRep)) THEN
  CALL dealloc_NParray(SGType2%tab_nq_OF_SRep,             &
                      'SGType2%tab_nq_OF_SRep',name_sub)
END IF

IF (allocated(SGType2%tab_Sum_nb_OF_SRep)) THEN
  CALL dealloc_NParray(SGType2%tab_Sum_nb_OF_SRep,        &
                      'SGType2%tab_Sum_nb_OF_SRep',name_sub)
END IF
IF (allocated(SGType2%tab_nb_OF_SRep)) THEN
  CALL dealloc_NParray(SGType2%tab_nb_OF_SRep,        &
                      'SGType2%tab_nb_OF_SRep',name_sub)
END IF
END SUBROUTINE dealloc_SGType2

SUBROUTINE SGType2_2TOSGType2_1(SGType2_1,SGType2_2)

TYPE (param_SGType2), intent(inout) :: SGType2_1
TYPE (param_SGType2), intent(in)    :: SGType2_2

integer :: i

character (len=*), parameter :: name_sub='SGType2_2TOSGType2_1'


SGType2_1%L1_SparseGrid = SGType2_2%L1_SparseGrid
SGType2_1%L2_SparseGrid = SGType2_2%L2_SparseGrid
SGType2_1%Num_OF_Lmax   = SGType2_2%Num_OF_Lmax

SGType2_1%nDind_SmolyakGrids = SGType2_2%nDind_SmolyakGrids
SGType2_1%nb_SG              = SGType2_2%nb_SG


IF (allocated(SGType2_2%nDind_DPG)) THEN
  CALL alloc_NParray(SGType2_1%nDind_DPG,(/ SGType2_1%nb_SG /),            &
                    'SGType2_1%nDind_DPG',name_sub)
  DO i=1,SGType2_1%nb_SG
    SGType2_1%nDind_DPG(i) = SGType2_2%nDind_DPG(i)
  END DO
END IF

IF (allocated(SGType2_2%nDind_DPB)) THEN
  CALL alloc_NParray(SGType2_1%nDind_DPB,(/ SGType2_1%nb_SG /),            &
                    'SGType2_1%nDind_DPB',name_sub)
  DO i=1,SGType2_1%nb_SG
    SGType2_1%nDind_DPB(i) = SGType2_2%nDind_DPB(i)
  END DO
END IF

IF (allocated(SGType2_2%tab_iB_OF_SRep_TO_iB)) THEN
  CALL alloc_NParray(SGType2_1%tab_iB_OF_SRep_TO_iB,                    &
                          shape(SGType2_2%tab_iB_OF_SRep_TO_iB),        &
                    'SGType2_1%tab_iB_OF_SRep_TO_iB',name_sub)
  SGType2_1%tab_iB_OF_SRep_TO_iB(:) = SGType2_2%tab_iB_OF_SRep_TO_iB
END IF

IF (allocated(SGType2_2%tab_Sum_nq_OF_SRep)) THEN
  CALL alloc_NParray(SGType2_1%tab_Sum_nq_OF_SRep,                      &
                          shape(SGType2_2%tab_Sum_nq_OF_SRep),          &
                    'SGType2_1%tab_Sum_nq_OF_SRep',name_sub)
  SGType2_1%tab_Sum_nq_OF_SRep(:) = SGType2_2%tab_Sum_nq_OF_SRep
END IF
IF (allocated(SGType2_2%tab_nq_OF_SRep)) THEN
  CALL alloc_NParray(SGType2_1%tab_nq_OF_SRep,                      &
                          shape(SGType2_2%tab_nq_OF_SRep),          &
                    'SGType2_1%tab_nq_OF_SRep',name_sub)
  SGType2_1%tab_nq_OF_SRep(:) = SGType2_2%tab_nq_OF_SRep
END IF


IF (allocated(SGType2_2%tab_Sum_nb_OF_SRep)) THEN
  CALL alloc_NParray(SGType2_1%tab_Sum_nb_OF_SRep,                      &
                          shape(SGType2_2%tab_Sum_nb_OF_SRep),          &
                    'SGType2_1%tab_Sum_nb_OF_SRep',name_sub)
  SGType2_1%tab_Sum_nb_OF_SRep(:) = SGType2_2%tab_Sum_nb_OF_SRep
END IF
IF (allocated(SGType2_2%tab_nb_OF_SRep)) THEN
  CALL alloc_NParray(SGType2_1%tab_nb_OF_SRep,                      &
                          shape(SGType2_2%tab_nb_OF_SRep),          &
                    'SGType2_1%tab_nb_OF_SRep',name_sub)
  SGType2_1%tab_nb_OF_SRep(:) = SGType2_2%tab_nb_OF_SRep
END IF

END SUBROUTINE SGType2_2TOSGType2_1

! from an index iq (global index of the multidimentional Smolyak grid) get:
!  - iSG:  the numero of the direct-product grid of the Smolyak grid
!  - iqSG: the numero of the grid point from the "iSG" direct-product grid
SUBROUTINE get_iqSG_iSG_FROM_iq(iSG,iqSG,iq,SGType2,OldPara,err_sub)

integer, intent(inout) :: iSG,iqSG
integer, intent(inout) :: err_sub

 TYPE (OldParam), intent(inout), optional :: OldPara

integer, intent(in) :: iq
TYPE (param_SGType2), intent(in) :: SGType2

integer :: nq,iSG_loc
!----- for debuging --------------------------------------------------
logical, parameter :: debug=.FALSE.
!logical, parameter :: debug=.TRUE.
character (len=*), parameter :: name_sub = 'get_iqSG_iSG_FROM_iq'
!-----------------------------------------------------------
err_sub = 0
IF (debug) THEN
  write(out_unitp,*) 'BEGINNING ',name_sub
  CALL flush_perso(out_unitp)
END IF

IF (present(OldPara)) THEN
  !write(6,*) 'OldPara ',name_sub,OldPara
  iSG = OldPara%i_SG
END IF


!write(6,*) 'alloc tab_Sum_nq_OF_SRep',allocated(SGType2%tab_Sum_nq_OF_SRep)
!flush(6)
IF (iSG > 1 .AND. iSG <= size(SGType2%tab_Sum_nq_OF_SRep)) THEN
  iqSG = iq - SGType2%tab_Sum_nq_OF_SRep(iSG-1)

  IF (iqSG < 1) THEN
    iqSG = iq
    iSG  = 1
  END IF
ELSE
  iqSG  = iq
  iSG   = 1
END IF


DO iSG_loc=iSG,SGType2%nb_SG
  nq = SGType2%tab_nq_OF_SRep(iSG_loc)
  IF (iqSG <= nq) EXIT
  iqSG = iqSG - nq
END DO

iSG       = iSG_loc
IF (present(OldPara)) THEN
  OldPara%i_SG = iSG
  OldPara%iq   = iq
END IF

IF (iqSG < 1) STOP 'iqSG < 1'


IF (debug) THEN
  write(out_unitp,*) 'iq,iSG,iqSG',iq,iSG,iqSG
  write(out_unitp,*) 'END ',name_sub
END IF

END SUBROUTINE get_iqSG_iSG_FROM_iq

! from an index iq (global index of the multidimentional Smolyak grid) get:
!  - i_SG:  the numero of the direct-product grid of the Smolyak grid
!  - iq_SG: the numero of the grid point from the "i_SG" direct-product grid
!  - Tabil(:): the "l" indexes corresponding to "i_SG"
!  - Tabiq(:): the "iq" indexes corresponding to "iq_SG"
SUBROUTINE get_Tabiq_Tabil_FROM_iq(Tabiq,Tabil,i_SG,iq_SG,iq,SGType2,OldPara,err_sub)

integer, intent(inout) :: Tabiq(:)
integer, intent(inout) :: Tabil(:)
integer, intent(inout) :: i_SG,iq_SG
integer, intent(inout) :: err_sub

 TYPE (OldParam), intent(inout), optional :: OldPara

integer, intent(in) :: iq
TYPE (param_SGType2), intent(in) :: SGType2

integer :: nq,i_SG_loc
logical :: old_l
!----- for debuging --------------------------------------------------
logical, parameter :: debug=.FALSE.
!logical, parameter :: debug=.TRUE.
character (len=*), parameter :: name_sub = 'get_Tabiq_Tabil_FROM_iq'
!-----------------------------------------------------------
err_sub = 0
IF (debug) THEN
  write(out_unitp,*) 'BEGINNING ',name_sub
  CALL flush_perso(out_unitp)
END IF

!first calculation of i_SG and iq_SG from iq and OldPara (if available)
IF (present(OldPara)) THEN
  !write(6,*) 'OldPara ',name_sub,OldPara
  i_SG = OldPara%i_SG
END IF

IF (i_SG > 1 .AND. i_SG <= size(SGType2%tab_Sum_nq_OF_SRep)) THEN
  iq_SG = iq - SGType2%tab_Sum_nq_OF_SRep(i_SG-1)

  IF (iq_SG < 1) THEN
    iq_SG = iq
    i_SG  = 1
  END IF
ELSE
  iq_SG  = iq
  i_SG   = 1
END IF


DO i_SG_loc=i_SG,SGType2%nb_SG
  nq = SGType2%tab_nq_OF_SRep(i_SG_loc)
  IF (iq_SG <= nq) EXIT
  iq_SG = iq_SG - nq
END DO

i_SG       = i_SG_loc
IF (iq_SG < 1) STOP 'iq_SG < 1'

  IF (debug) write(out_unitp,*) ' i_SG,iq_SG,iq',i_SG,iq_SG,iq
  CALL flush_perso(out_unitp)


!2d calculation of Tabil from i_SG and OldPara (if available)
IF (present(OldPara)) THEN

  old_l = (i_SG >= OldPara%i_SG) .AND. allocated(OldPara%tab_l_AT_SG)

  IF (old_l) THEN

    Tabil(:) = OldPara%tab_l_AT_SG
    DO i_SG_loc=OldPara%i_SG+1,i_SG
      CALL ADD_ONE_TO_nDindex(SGType2%nDind_SmolyakGrids,Tabil,iG=i_SG_loc,err_sub=err_sub)
    END DO

  ELSE
    CALL calc_nDindex(SGType2%nDind_SmolyakGrids,i_SG,Tabil,err_sub)
  END IF

ELSE
  CALL calc_nDindex(SGType2%nDind_SmolyakGrids,i_SG,Tabil,err_sub)
END IF

IF (err_sub /= 0) THEN
  write(out_unitp,*) ' SGType2%nDind_SmolyakGrids'
  write(out_unitp,*) ' ERROR in ',name_sub
  write(out_unitp,*) ' i_SG,iq_SG,iq',i_SG,iq_SG,iq
  write(out_unitp,*) ' Tabil',Tabil
  write(out_unitp,*) '  from SGType2%nDind_SmolyakGrids',i_SG
  err_sub = 1
  RETURN
    !STOP 'calc_nDindex'
END IF
  IF (debug) write(out_unitp,*) ' Tabil',i_SG,' : ',Tabil
  CALL flush_perso(out_unitp)


! Save the parameters in OldPara if OldPara is present
IF (present(OldPara)) THEN
  OldPara%i_SG        = i_SG
  OldPara%iq_SG       = iq_SG
  OldPara%iq          = iq
  OldPara%tab_l_AT_SG = Tabil
END IF

!3d calculation of Tabiq from i_SG and iq_SG
  Tabiq(:) = 0
  !CALL calc_nDval_m1(Tabiq,SGType2%tab_nq_OF_SRep(i_SG),nDsize,size(Tabil))
  !CALL calc_nDindex(SGType2%nDind_DPG(i_SG),iq_SG,Tabiq,err_sub)
  IF (err_sub /= 0) THEN
    write(out_unitp,*) ' SGType2%nDind_DPG(i_SG)'
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) ' i_SG,iq_SG,iq',i_SG,iq_SG,iq
    write(out_unitp,*) ' Tabiq',Tabiq
    write(out_unitp,*) ' Tabil',Tabil
    write(out_unitp,*) '  from SGType2%nDind_DPG(i_SG)',i_SG
    err_sub = 2
    RETURN
    !STOP 'calc_nDindex'
   END IF

  !IF (debug) write(out_unitp,*) ' Tabiq',Tabiq
  !CALL flush_perso(out_unitp)

IF (debug) THEN
  write(out_unitp,*) 'iq,i_SG,iq_SG',iq,i_SG,iq_SG
  write(out_unitp,*) 'iq,i_SG,Tabil',iq,i_SG,':',Tabil
  write(out_unitp,*) 'iq,i_SG,Tabiq',iq,i_SG,':',Tabiq
  write(out_unitp,*) 'END ',name_sub
  CALL flush_perso(out_unitp)
END IF

END SUBROUTINE get_Tabiq_Tabil_FROM_iq



SUBROUTINE get_Tabiq_Tabil_FROM_iq_old(Tabiq,Tabil,i_SG,iq_SG,iq,SGType2)

integer, intent(inout) :: Tabiq(:)
integer, intent(inout) :: Tabil(:)
integer, intent(inout) :: i_SG,iq_SG

integer, intent(in) :: iq
TYPE (param_SGType2), intent(in) :: SGType2


integer :: nq
integer :: err_sub

!----- for debuging --------------------------------------------------
 logical, parameter :: debug=.FALSE.
!logical, parameter :: debug=.TRUE.
character (len=*), parameter :: name_sub = 'get_Tabiq_Tabil_FROM_iq_old'
!-----------------------------------------------------------
IF (debug) THEN
  write(out_unitp,*) 'BEGINNING ',name_sub
END IF

iq_SG           = iq

DO i_SG=1,SGType2%nb_SG
  nq = SGType2%nDind_DPG(i_SG)%Max_nDI
  IF (iq_SG <= nq) EXIT
  iq_SG = iq_SG - nq
END DO

  CALL calc_nDindex(SGType2%nDind_SmolyakGrids,i_SG,Tabil,err_sub)
  IF (err_sub /= 0) THEN
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) '  from SGType2%nDind_SmolyakGrids',i_SG
    STOP 'calc_nDindex'
  END IF

  CALL calc_nDindex(SGType2%nDind_DPG(i_SG),iq_SG,Tabiq,err_sub)
  IF (err_sub /= 0) THEN
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) '  from SGType2%nDind_DPG(i_SG)',i_SG
    STOP 'calc_nDindex'
  END IF

IF (debug) THEN
  write(out_unitp,*) 'iq,i_SG,iq_SG',iq,i_SG,iq_SG
  write(out_unitp,*) 'iq,i_SG,Tabil',iq,i_SG,':',Tabil
  write(out_unitp,*) 'iq,i_SG,Tabiq',iq,i_SG,':',Tabiq
  write(out_unitp,*) 'END ',name_sub
END IF

END SUBROUTINE get_Tabiq_Tabil_FROM_iq_old

END MODULE mod_param_SGType2
