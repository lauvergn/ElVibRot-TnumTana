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
MODULE mod_basis_BtoG_GtoB_SGType4
USE mod_system
!$ USE omp_lib, only : OMP_GET_THREAD_NUM
USE mod_nDindex
USE mod_basis_set_alloc
USE mod_param_SGType2
USE mod_basis_RCVec_SGType4, only: typervec, typecvec, &
                                   alloc_typervec, alloc_typecvec, &
                                   dealloc_typervec, dealloc_typecvec
IMPLICIT NONE

PRIVATE

TYPE Type_SmolyakRep
  integer                      :: nb0      = 0         ! to deal with several electronic PES, rot basis, or channels (HAC)
  logical                      :: Grid     = .FALSE.
  logical                      :: Delta    = .FALSE.
  TYPE (TypeRVec), allocatable :: SmolyakRep(:)
CONTAINS
  PROCEDURE, PRIVATE, PASS(SRep1) :: SmolyakRep2_TO_SmolyakRep1
  PROCEDURE, PRIVATE, PASS(SRep2) :: SmolyakRep2_TO_tabR1
  PROCEDURE, PRIVATE, PASS(SRep1) :: R2_TO_SmolyakRep1
  PROCEDURE, PRIVATE, PASS(SRep1) :: tabR2_TO_SmolyakRep1
  GENERIC,   PUBLIC  :: assignment(=) => SmolyakRep2_TO_SmolyakRep1,    &
                              SmolyakRep2_TO_tabR1,tabR2_TO_SmolyakRep1,&
                               R2_TO_SmolyakRep1

END TYPE Type_SmolyakRep
TYPE Type_SmolyakRepC
  integer                      :: nb0      = 0         ! to deal with several electronic PES, rot basis, or channels (HAC)
  logical                      :: Grid     = .FALSE.
  logical                      :: Delta    = .FALSE.
  TYPE (TypeCVec), allocatable :: SmolyakRep(:)
END TYPE Type_SmolyakRepC


INTERFACE operator(*)
  module procedure SmolyakRep1_TIME_SmolyakRep2,SmolyakRepC1_TIME_SmolyakRepC2
  module procedure SmolyakRep1_TIME_R2,R1_TIME_SmolyakRep2
END INTERFACE
INTERFACE operator(+)
  module procedure SmolyakRep1_PLUS_SmolyakRep2,SmolyakRepC1_PLUS_SmolyakRepC2
END INTERFACE
INTERFACE operator(-)
  module procedure SmolyakRep1_MINUS_SmolyakRep2,SmolyakRepC1_MINUS_SmolyakRepC2
END INTERFACE

PUBLIC  getbis_tab_nq, getbis_tab_nb
PUBLIC  Type_SmolyakRep,  alloc2_SmolyakRep,  dealloc_SmolyakRep,  DerivOp_TO3_GSmolyakRep
PUBLIC  Type_SmolyakRepC, alloc2_SmolyakRepC, dealloc_SmolyakRepC, DerivOp_TO3_GSmolyakRepC
PUBLIC  alloc_SmolyakRep_only
PUBLIC  Write_SmolyakRep, alloc_SmolyakRep
PUBLIC  tabR_AT_iG_TO_tabPackedBasis, tabPackedBasis_TO_tabR_AT_iG
PUBLIC  tabR2bis_TO_SmolyakRep1
PUBLIC  tabR2grid_TO_tabR1_AT_iG, tabR2gridbis_TO_tabR1_AT_iG
PUBLIC  BDP_TO_GDP_OF_SmolyakRep, GDP_TO_BDP_OF_SmolyakRep
PUBLIC  DerivOp_TO_RDP_OF_SmolaykRep
PUBLIC  Set_weight_TO_SmolyakRep, Get_weight_FROM_OneDP, &
        dot_product_SmolyakRep_Grid, dot_product_SmolyakRep_Basis
PUBLIC  GSmolyakRep_TO3_BSmolyakRep, BSmolyakRep_TO3_GSmolyakRep
PUBLIC  GSmolyakRep_TO_BSmolyakRep,  BSmolyakRep_TO_GSmolyakRep
PUBLIC  Set_tables_FOR_SmolyakRepBasis_TO_tabPackedBasis
PUBLIC  SmolyakRep2_TO_tabR1bis, SmolyakRepC2_TO_tabC1bis, tabC2bis_TO_SmolyakRepC1
PUBLIC  SmolyakRepBasis_TO_tabPackedBasis, tabPackedBasis_TO_SmolyakRepBasis

PUBLIC  typeRvec, alloc_typeRvec, dealloc_typeRvec
PUBLIC  typeCvec, alloc_typeCvec, dealloc_typeCvec
PUBLIC  operator(*), operator(+), operator(-)

CONTAINS

!=======================================================================================
! added for SmolyakRep
!---------------------------------------------------------------------------------------
SUBROUTINE alloc_SmolyakRep_only(SRep,nb_SG,delta,grid,nb0)
  USE mod_system
  USE mod_basis_set_alloc
  IMPLICIT NONE

  TYPE(Type_SmolyakRep),           intent(inout)         :: SRep
  Integer,                         intent(in)            :: nb_SG
  logical,                         intent(in),  optional :: delta,grid
  integer,                         intent(in),  optional :: nb0

  integer               :: iG,nb_B
  integer, allocatable  :: tab_n(:)


  !write(out_unitp,*) 'in alloc_SmolyakRep, shape tab_ind',shape(tab_ind) ; flush(out_unitp)
  CALL dealloc_SmolyakRep(SRep)

  IF (present(nb0)) THEN
    SRep%nb0 = nb0
  ELSE
    SRep%nb0 = 1
  END IF
  IF (SRep%nb0 < 1) SRep%nb0 = 1

  IF (present(delta)) THEN
    SRep%delta = delta
  ELSE
    SRep%delta = .FALSE.
  END IF

  IF (present(grid)) THEN
    SRep%grid = grid
  ELSE
    SRep%grid = .FALSE.
  END IF

  !write(out_unitp,*) 'Alloc Smolyak Rep' ; flush(out_unitp)

  allocate(SRep%SmolyakRep( nb_SG ))

  !write(out_unitp,*) 'Size Smolyak Rep:',nb_B ; flush(out_unitp)

END SUBROUTINE alloc_SmolyakRep_only
!=======================================================================================

SUBROUTINE alloc_SmolyakRep(SRep,tab_ind,tab_ba,delta,grid,nb0)
USE mod_system
USE mod_basis_set_alloc
IMPLICIT NONE

TYPE(Type_SmolyakRep),           intent(inout)         :: SRep
integer,                         intent(in)            :: tab_ind(:,:) ! tab_ind(D,MaxnD)
TYPE(basis),                     intent(in)            :: tab_ba(0:,:) ! tab_ba(0:L,D)
logical,                         intent(in),  optional :: delta,grid
integer,                         intent(in),  optional :: nb0

integer               :: iG,nb_B
integer, allocatable  :: tab_n(:)


!write(out_unitp,*) 'in alloc_SmolyakRep, shape tab_ind',shape(tab_ind) ; flush(out_unitp)
CALL dealloc_SmolyakRep(SRep)

IF (present(nb0)) THEN
  SRep%nb0 = nb0
ELSE
  SRep%nb0 = 1
END IF
IF (SRep%nb0 < 1) SRep%nb0 = 1

IF (present(delta)) THEN
  SRep%delta = delta
ELSE
  SRep%delta = .FALSE.
END IF

IF (present(grid)) THEN
  SRep%grid = grid
ELSE
  SRep%grid = .FALSE.
END IF

!write(out_unitp,*) 'Alloc Smolyak Rep' ; flush(out_unitp)

allocate(SRep%SmolyakRep( size(tab_ind(1,:)) ))
CALL alloc_NParray(tab_n,shape(tab_ind(:,1)),'tab_n','alloc_SmolyakRep')

DO iG=1,size(tab_ind(1,:))

  IF (SRep%grid) THEN
    tab_n = getbis_tab_nq(tab_ind(:,iG),tab_ba)
  ELSE
    tab_n = getbis_tab_nb(tab_ind(:,iG),tab_ba)
  END IF
  !write(out_unitp,*) iG,'tab_n',tab_n ; flush(out_unitp)
  CALL alloc_TypeRVec(SRep%SmolyakRep(iG),nvec=product(tab_n)*SRep%nb0)
END DO

IF (allocated(tab_n)) CALL dealloc_NParray(tab_n,'tab_n','alloc_SmolyakRep')

nb_B = Size_SmolyakRep(SRep)
!write(out_unitp,*) 'Size Smolyak Rep:',nb_B ; flush(out_unitp)

END SUBROUTINE alloc_SmolyakRep
SUBROUTINE alloc2_SmolyakRep(SRep,nDind_SmolyakRep,tab_ba,delta,grid,nb0)
USE mod_system
USE mod_basis_set_alloc
IMPLICIT NONE

TYPE(Type_SmolyakRep),           intent(inout)         :: SRep
TYPE(basis),                     intent(in)            :: tab_ba(0:,:) ! tab_ba(0:L,D)
logical,                         intent(in),  optional :: delta,grid
TYPE (Type_nDindex),             intent(in)            :: nDind_SmolyakRep
integer,                         intent(in),  optional :: nb0

integer  :: iG,nb_B,err_sub
integer, allocatable  :: tab_n(:),tab_l(:)


CALL dealloc_SmolyakRep(SRep)

IF (present(nb0)) THEN
  SRep%nb0 = nb0
ELSE
  SRep%nb0 = 1
END IF
IF (SRep%nb0 < 1) SRep%nb0 = 1

IF (present(delta)) THEN
  SRep%delta = delta
ELSE
  SRep%delta = .FALSE.
END IF

IF (present(grid)) THEN
  SRep%grid = grid
ELSE
  SRep%grid = .FALSE.
END IF

!write(out_unitp,*) 'Alloc Smolyak Rep'

allocate(SRep%SmolyakRep(nDind_SmolyakRep%Max_nDI))
CALL alloc_NParray(tab_n,[nDind_SmolyakRep%ndim],'tab_n','alloc2_SmolyakRep')

IF (allocated(nDind_SmolyakRep%Tab_nDval)) THEN

  DO iG=1,nDind_SmolyakRep%Max_nDI

    IF (SRep%grid) THEN
      tab_n = getbis_tab_nq(nDind_SmolyakRep%Tab_nDval(:,iG),tab_ba)
    ELSE
      tab_n = getbis_tab_nb(nDind_SmolyakRep%Tab_nDval(:,iG),tab_ba)
    END IF
    CALL alloc_TypeRVec(SRep%SmolyakRep(iG),nvec=product(tab_n)*SRep%nb0)
  END DO
ELSE
  CALL alloc_NParray(tab_l,[nDind_SmolyakRep%ndim],'tab_l','alloc2_SmolyakRep')

  CALL init_nDval_OF_nDindex(nDind_SmolyakRep,tab_l,err_sub)
  IF (err_sub /= 0) STOP 'init_nDval_OF_nDindex'

  DO iG=1,nDind_SmolyakRep%Max_nDI
    CALL ADD_ONE_TO_nDindex(nDind_SmolyakRep,tab_l,iG=iG,err_sub=err_sub)
    IF (err_sub /= 0) STOP 'ADD_ONE_TO_nDindex'

    IF (SRep%grid) THEN
      tab_n = getbis_tab_nq(tab_l,tab_ba)
    ELSE
      tab_n = getbis_tab_nb(tab_l,tab_ba)
    END IF
    CALL alloc_TypeRVec(SRep%SmolyakRep(iG),nvec=product(tab_n)*SRep%nb0)
  END DO
  CALL dealloc_NParray(tab_l,'tab_l','alloc2_SmolyakRep')

END IF
CALL dealloc_NParray(tab_n,'tab_n','alloc2_SmolyakRep')

!nb_B = Size_SmolyakRep(SRep)
!write(out_unitp,*) 'Size Smolyak Rep:',nb_B

END SUBROUTINE alloc2_SmolyakRep

SUBROUTINE alloc_SmolyakRepC(SRep,tab_ind,tab_ba,delta,grid,nb0)
USE mod_system
USE mod_basis_set_alloc
IMPLICIT NONE

TYPE(Type_SmolyakRepC),          intent(inout)         :: SRep
integer,                         intent(in)            :: tab_ind(:,:) ! tab_ind(D,MaxnD)
TYPE(basis),                     intent(in)            :: tab_ba(0:,:) ! tab_ba(0:L,D)
logical,                         intent(in),  optional :: delta,grid
integer,                         intent(in),  optional :: nb0

integer               :: iG,nb_B
integer, allocatable  :: tab_n(:)


!write(out_unitp,*) 'in alloc_SmolyakRep, shape tab_ind',shape(tab_ind)
CALL dealloc_SmolyakRepC(SRep)

IF (present(nb0)) THEN
  SRep%nb0 = nb0
ELSE
  SRep%nb0 = 1
END IF
IF (SRep%nb0 < 1) SRep%nb0 = 1

IF (present(delta)) THEN
  SRep%delta = delta
ELSE
  SRep%delta = .FALSE.
END IF

IF (present(grid)) THEN
  SRep%grid = grid
ELSE
  SRep%grid = .FALSE.
END IF

!write(out_unitp,*) 'Alloc Smolyak Rep'

allocate(SRep%SmolyakRep( size(tab_ind(1,:)) ))
CALL alloc_NParray(tab_n,shape(tab_ind(:,1)),'tab_n','alloc_SmolyakRepC')

DO iG=1,size(tab_ind(1,:))

  IF (SRep%grid) THEN
    tab_n = getbis_tab_nq(tab_ind(:,iG),tab_ba)
  ELSE
    tab_n = getbis_tab_nb(tab_ind(:,iG),tab_ba)
  END IF
  !write(out_unitp,*) iG,'tab_n',tab_n
  CALL alloc_TypeCVec(SRep%SmolyakRep(iG),nvec=product(tab_n)*SRep%nb0)
END DO

IF (allocated(tab_n)) CALL dealloc_NParray(tab_n,'tab_n','alloc_SmolyakRepC')

!nb_B = Size_SmolyakRep(SRep)
!write(out_unitp,*) 'Size Smolyak Rep:',nb_B

END SUBROUTINE alloc_SmolyakRepC
SUBROUTINE alloc2_SmolyakRepC(SRep,nDind_SmolyakRep,tab_ba,delta,grid,nb0)
USE mod_system
USE mod_basis_set_alloc
IMPLICIT NONE

TYPE(Type_SmolyakRepC),          intent(inout)         :: SRep
TYPE(basis),                     intent(in)            :: tab_ba(0:,:) ! tab_ba(0:L,D)
logical,                         intent(in),  optional :: delta,grid
TYPE (Type_nDindex),             intent(in)            :: nDind_SmolyakRep
integer,                         intent(in),  optional :: nb0

integer  :: iG,nb_B,err_sub
integer, allocatable  :: tab_n(:),tab_l(:)


CALL dealloc_SmolyakRepC(SRep)

IF (present(nb0)) THEN
  SRep%nb0 = nb0
ELSE
  SRep%nb0 = 1
END IF
IF (SRep%nb0 < 1) SRep%nb0 = 1

IF (present(delta)) THEN
  SRep%delta = delta
ELSE
  SRep%delta = .FALSE.
END IF

IF (present(grid)) THEN
  SRep%grid = grid
ELSE
  SRep%grid = .FALSE.
END IF

!write(out_unitp,*) 'Alloc Smolyak Rep'

allocate(SRep%SmolyakRep(nDind_SmolyakRep%Max_nDI))
CALL alloc_NParray(tab_n,[nDind_SmolyakRep%ndim],'tab_n','alloc2_SmolyakRepC')

IF (allocated(nDind_SmolyakRep%Tab_nDval)) THEN

  DO iG=1,nDind_SmolyakRep%Max_nDI

    IF (SRep%grid) THEN
      tab_n = getbis_tab_nq(nDind_SmolyakRep%Tab_nDval(:,iG),tab_ba)
    ELSE
      tab_n = getbis_tab_nb(nDind_SmolyakRep%Tab_nDval(:,iG),tab_ba)
    END IF
    CALL alloc_TypeCVec(SRep%SmolyakRep(iG),nvec=product(tab_n)*SRep%nb0)
  END DO
ELSE
  CALL alloc_NParray(tab_l,[nDind_SmolyakRep%ndim],'tab_l','alloc2_SmolyakRepC')

  CALL init_nDval_OF_nDindex(nDind_SmolyakRep,tab_l,err_sub)
  IF (err_sub /= 0) STOP 'init_nDval_OF_nDindex'

  DO iG=1,nDind_SmolyakRep%Max_nDI
    CALL ADD_ONE_TO_nDindex(nDind_SmolyakRep,tab_l,iG=iG,err_sub=err_sub)
    IF (err_sub /= 0) STOP 'ADD_ONE_TO_nDindex'

    IF (SRep%grid) THEN
      tab_n = getbis_tab_nq(tab_l,tab_ba)
    ELSE
      tab_n = getbis_tab_nb(tab_l,tab_ba)
    END IF
    CALL alloc_TypeCVec(SRep%SmolyakRep(iG),nvec=product(tab_n)*SRep%nb0)
  END DO
  CALL dealloc_NParray(tab_l,'tab_l','alloc2_SmolyakRepC')

END IF
CALL dealloc_NParray(tab_n,'tab_n','alloc2_SmolyakRepC')

!nb_B = Size_SmolyakRep(SRep)
!write(out_unitp,*) 'Size Smolyak Rep:',nb_B

END SUBROUTINE alloc2_SmolyakRepC

SUBROUTINE dealloc_SmolyakRep(SRep)
  USE mod_system
  USE mod_MPI_aux
  IMPLICIT NONE

  TYPE(Type_SmolyakRep), intent(inout)     :: SRep

  integer               :: iG
  integer               :: d1,d2

  IF (allocated(SRep%SmolyakRep)) THEN

    d1=lbound(SRep%SmolyakRep,dim=1)
    d2=ubound(SRep%SmolyakRep,dim=1)
    ! IF(openmpi .AND. MPI_scheme==1) THEN
    !   d1=iGs_MPI(1,MPI_id)
    !   d2=iGs_MPI(2,MPI_id)
    ! ENDIF
    IF(openmpi) THEN
      IF(MPI_scheme==1) THEN
        d1=iGs_MPI(1,MPI_id)
        d2=iGs_MPI(2,MPI_id)
      ELSEIF(MPI_scheme==3) THEN
        IF(MPI_nodes_p0) THEN
          d1=iGs_MPI(1,MPI_id)
          d2=iGs_MPI(2,MPI_sub_id(2))
        ELSE
          d1=0
          d2=-1
        ENDIF
      ENDIF
    ENDIF

    !DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)
    DO iG=d1,d2
      CALL dealloc_TypeRVec(SRep%SmolyakRep(iG))
    END DO
    deallocate(SRep%SmolyakRep)
  END IF

  SRep%Grid     = .FALSE.
  SRep%Delta    = .FALSE.
  SRep%nb0      = 0

END SUBROUTINE dealloc_SmolyakRep

SUBROUTINE dealloc_SmolyakRepC(SRep)
  USE mod_system
  IMPLICIT NONE

  TYPE(Type_SmolyakRepC), intent(inout)     :: SRep

  integer               :: iG

  IF (allocated(SRep%SmolyakRep)) THEN
    DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)
      CALL dealloc_TypeCVec(SRep%SmolyakRep(iG))
    END DO
    deallocate(SRep%SmolyakRep)
  END IF

  SRep%Grid     = .FALSE.
  SRep%Delta    = .FALSE.
  SRep%nb0      = 0

END SUBROUTINE dealloc_SmolyakRepC

SUBROUTINE Write_SmolyakRep(SRep)
USE mod_system
IMPLICIT NONE

TYPE(Type_SmolyakRep), intent(in)     :: SRep

integer               :: iG

  write(out_unitp,*) 'Grid',SRep%Grid
  write(out_unitp,*) 'Delta',SRep%Delta
  write(out_unitp,*) 'nb0',SRep%nb0
  write(out_unitp,*) 'alloc?',allocated(SRep%SmolyakRep)

IF (allocated(SRep%SmolyakRep)) THEN
  write(out_unitp,*) '======== Smolyak Rep ============================'
  DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)
    write(out_unitp,*) iG,size(SRep%SmolyakRep(iG)%V),SRep%SmolyakRep(iG)%V
  END DO
END IF
CALL flush_perso(out_unitp)
END SUBROUTINE Write_SmolyakRep

SUBROUTINE Write_SmolyakRepC(SRep)
USE mod_system
IMPLICIT NONE

TYPE(Type_SmolyakRepC), intent(in)     :: SRep

integer               :: iG

  write(out_unitp,*) 'Grid',SRep%Grid
  write(out_unitp,*) 'Delta',SRep%Delta
  write(out_unitp,*) 'nb0',SRep%nb0
  write(out_unitp,*) 'alloc?',allocated(SRep%SmolyakRep)

IF (allocated(SRep%SmolyakRep)) THEN
  write(out_unitp,*) '======== Smolyak Rep ============================'
  DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)
    write(out_unitp,*) iG,size(SRep%SmolyakRep(iG)%V),SRep%SmolyakRep(iG)%V
  END DO
END IF
CALL flush_perso(out_unitp)
END SUBROUTINE Write_SmolyakRepC


FUNCTION Size_SmolyakRep(SRep) RESULT(nb_BG)
USE mod_system
IMPLICIT NONE

TYPE(Type_SmolyakRep), intent(in)     :: SRep
integer                               :: nb_BG

integer               :: iG,nG


nb_BG = 0
IF (allocated(SRep%SmolyakRep)) THEN
  DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)
    nb_BG = nb_BG + size(SRep%SmolyakRep(iG)%V)
  END DO
END IF

END FUNCTION Size_SmolyakRep
FUNCTION MaxVal_SmolyakRep(SRep) RESULT(maxSRep)
USE mod_system
IMPLICIT NONE

TYPE(Type_SmolyakRep), intent(in)     :: SRep
real (kind=Rkind)                     :: maxSRep

integer               :: iG

maxSRep = ZERO
IF (allocated(SRep%SmolyakRep)) THEN
  DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)
    maxSRep = max(maxSRep,maxval(abs(SRep%SmolyakRep(iG)%V)))
  END DO
END IF

END FUNCTION MaxVal_SmolyakRep

FUNCTION Size_SmolyakRepC(SRep) RESULT(nb_BG)
USE mod_system
IMPLICIT NONE

TYPE(Type_SmolyakRepC), intent(in)    :: SRep
integer                               :: nb_BG

integer               :: iG,nG


nb_BG = 0
IF (allocated(SRep%SmolyakRep)) THEN
  DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)
    nb_BG = nb_BG + size(SRep%SmolyakRep(iG)%V)
  END DO
END IF

END FUNCTION Size_SmolyakRepC
FUNCTION MaxVal_SmolyakRepC(SRep) RESULT(maxSRep)
USE mod_system
IMPLICIT NONE

TYPE(Type_SmolyakRepC), intent(in)    :: SRep
real (kind=Rkind)                     :: maxSRep

integer               :: iG

maxSRep = ZERO
IF (allocated(SRep%SmolyakRep)) THEN
  DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)
    maxSRep = max(maxSRep,maxval(abs(SRep%SmolyakRep(iG)%V)))
  END DO
END IF

END FUNCTION MaxVal_SmolyakRepC

!=======================================================================================
! mapping table
!---------------------------------------------------------------------------------------
SUBROUTINE Set_tables_FOR_SmolyakRepBasis_TO_tabPackedBasis(basis_SG)
  USE mod_system
  USE mod_basis_set_alloc
  USE mod_param_SGType2
  USE mod_nDindex
  USE mod_basis_BtoG_GtoB_SGType4_MPI,ONLY:Set_scheme_MPI,Mapping_table_MPI,           &
                                           Mapping_table_allocate_MPI
  USE mod_MPI_aux
  IMPLICIT NONE

  TYPE (basis),     intent(inout)        :: basis_SG

  integer               :: L,i,k,lk,iG,iBDP,iBSRep,nb_BG,nR,itabR,nDI,Max_Srep

  integer               :: tab_n(basis_SG%nDindB%ndim)

  integer, allocatable  :: tab_l(:)
  integer, allocatable  :: tab_ib(:),tab_nb(:)


  integer               :: MaxnD_with_id_and_L(basis_SG%nDindB%ndim,0:basis_SG%nDindB%Lmax)
  integer               :: max_NBB,max_nb,id,ib,iib,iBB,iVal,LL
  integer               :: Tab_inD_nDindB(basis_SG%nDindB%Max_nDI)
  integer               :: nDval(basis_SG%nDindB%ndim)

  TYPE (Type_nDindex)   :: nDind_DPB

  integer               :: ith,nb_thread
  integer               :: d1,d2

  integer (kind=ILkind) :: lMax_Srep

  INTEGER ::      nb_ticks_initial  ! initial value of the clock tick counter
  INTEGER ::      nb_ticks_final    ! final value of the clock tick counter
  INTEGER ::      nb_ticks_max      ! maximum value of the clock counter
  INTEGER ::      nb_ticks_sec      ! number of clock ticks per second
  INTEGER ::      nb_ticks          ! number of clock ticks of the code
  REAL    ::      elapsed_time      ! real time in seconds
  integer ::      err_sub

  !logical, parameter :: debug=.TRUE.
  logical, parameter :: debug=.FALSE.
  character (len=*), parameter :: name_sub='Set_tables_FOR_SmolyakRepBasis_TO_tabPackedBasis'

  write(out_unitp,*) 'BEGINNING ',name_sub
  IF (debug) THEN
    write(out_unitp,*) 'Write nDindB'
    !CALL write_nDindex(basis_SG%nDindB)
  !  IF (associated(basis_SG%nDindB%Tab_i_TO_l)) THEN
  !    write(out_unitp,*) 'tab_i_TO_l'
  !    DO id=1,basis_SG%nDindB%ndim
  !      CALL Write_IntVec(basis_SG%nDindB%tab_i_TO_l(id))
  !    END DO
  !  END IF
    write(out_unitp,*) 'END Write nDindB'
    CALL flush_perso(out_unitp)
  END IF

  CALL SYSTEM_CLOCK(COUNT_RATE=nb_ticks_sec, COUNT_MAX=nb_ticks_max)
  CALL SYSTEM_CLOCK(COUNT=nb_ticks_initial)

  ! Initialisations
  IF (debug) write(out_unitp,*) '---------------------------------------------------------------'
  MaxnD_with_id_and_L(:,:) = 0

  max_NBB = basis_SG%nDindB%Max_nDI
  Tab_inD_nDindB(:) = 0

  DO id=1,basis_SG%nDindB%ndim
    IF (debug) THEN
      write(out_unitp,*) '------------------------------------ id',id
      write(out_unitp,*) 'id,tab_i_TO_l',id,':',basis_SG%nDindB%tab_i_TO_l(id)%vec
      CALL flush_perso(out_unitp)
    END IF
    ib    = 1
    iVal  = 1

    ith = 0
    CALL init_nDval_OF_nDindex(basis_SG%nDindB,nDval,err_sub)
    IF (err_sub /= 0) STOP 'init_nDval_OF_nDindex in Mapping'
    DO iBB=1,max_NBB
      CALL ADD_ONE_TO_nDindex(basis_SG%nDindB,nDval,iG=iBB)
      !nDval(:) = basis_SG%nDindB%Tab_nDval(:,iBB)

      IF (nDval(id) /= ib) THEN
        ib   = nDval(id)
        iVal = 1
      END IF
      IF (nDval(id) == ib) THEN
        l = basis_SG%nDindB%tab_i_TO_l(id)%vec(ib)
        MaxnD_with_id_and_L(id,l) = iVal
      END IF
      iVal = iVal+1
    END DO
    max_NBB = MaxnD_with_id_and_L(id,0)
  END DO

  IF (debug) THEN
    write(out_unitp,*) '---------------------------------------------------------------'
    DO l=0,basis_SG%nDindB%Lmax
      write(out_unitp,*) ' MaxnD_with_id_and_l',l,':',MaxnD_with_id_and_L(:,l)
    END DO
    write(out_unitp,*) '---------------------------------------------------------------'
    CALL flush_perso(out_unitp)
  END IF

  IF (allocated(basis_SG%para_SGType2%tab_iB_OF_SRep_TO_iB)) THEN
    CALL dealloc_NParray(basis_SG%para_SGType2%tab_iB_OF_SRep_TO_iB,      &
                        'basis_SG%para_SGType2%tab_iB_OF_SRep_TO_iB',name_sub)
  END IF
  ! we calculate Max_Srep in short and long integer to check if Max_Srep is larger than huge(1)
  ! It is useless when the program is compliled with kind=8 for the integer (Equivalent to ILkind=8)
  lMax_Srep = sum(int(basis_SG%para_SGType2%tab_nb_OF_SRep(:),kind=ILkind))
  Max_Srep  = sum(basis_SG%para_SGType2%tab_nb_OF_SRep(:))

  write(out_unitp,*) 'nb of terms (grids)',size(basis_SG%para_SGType2%tab_nb_OF_SRep)
  write(out_unitp,*) 'nb_ba',basis_SG%nDindB%Max_nDI
  write(out_unitp,*) 'lMax_Srep',lMax_Srep
  write(out_unitp,*) 'Max_Srep',Max_Srep

  IF (lMax_Srep /= int(Max_Srep,kind=ILkind)) STOP 'ERROR Max_Srep is too large!!'
  CALL flush_perso(out_unitp)

  IF(openmpi) THEN
    CALL Set_scheme_MPI(basis_SG,lMax_Srep)
    CALL allocate_array(bounds_MPI,1,2,0,MPI_np-1)
    DO i_MPI=0,MPI_np-1
      bounds_MPI(1,i_MPI)=basis_SG%para_SGType2%tab_Sum_nb_OF_SRep(iGs_MPI(1,i_MPI))   &
                         -basis_SG%para_SGType2%tab_nb_OF_SRep(iGs_MPI(1,i_MPI))+1
      bounds_MPI(2,i_MPI)=basis_SG%para_SGType2%tab_Sum_nb_OF_SRep(iGs_MPI(2,i_MPI))
    ENDDO
  ENDIF

  IF(.NOT. openmpi) THEN
    CALL alloc_NParray(basis_SG%para_SGType2%tab_iB_OF_SRep_TO_iB,[Max_Srep], &
                      'basis_SG%para_SGType2%tab_iB_OF_SRep_TO_iB',name_sub)
    basis_SG%para_SGType2%tab_iB_OF_SRep_TO_iB(:) = 0
  ELSE
    ! allocate tab_iB_OF_SRep_TO_iB according to different MPI_scheme
    CALL Mapping_table_allocate_MPI(basis_SG,Max_Srep)
  ENDIF

  IF (.NOT. allocated(basis_SG%para_SGType2%nDind_SmolyakRep%Tab_nDval) ) THEN
    CALL alloc_NParray(tab_l, [basis_SG%nDindB%ndim],'tab_l',name_sub)
    CALL init_nDval_OF_nDindex(basis_SG%para_SGType2%nDind_SmolyakRep,tab_l)
    CALL dealloc_NParray(tab_l,'tab_l',name_sub)
  END IF

  !to be sure to have the correct number of threads
  nb_thread = basis_SG%para_SGType2%nb_tasks
  !nb_thread = 1

  IF(.NOT. debug) THEN
    write(out_unitp,'(a)')              'Tab(:) (%): [--0-10-20-30-40-50-60-70-80-90-100]'
    write(out_unitp,'(a)',ADVANCE='no') 'Tab(:) (%): ['
    CALL flush_perso(out_unitp)
  END IF

  !$OMP parallel                                                &
  !$OMP default(none)                                           &
  !$OMP shared(basis_SG,nb_thread,Tab_inD_nDindB)               &
  !$OMP shared(out_unitp,MaxnD_with_id_and_L)                   &
  !$OMP shared(iGs_MPI,MPI_id,MPI_scheme)                       &
  !$OMP shared(MPI_np,n_level2,openmpi)                         &
  !$OMP private(iG,tab_l,ith,tab_nb,tab_ib,nDI,iBDP,iBSRep)     &
  !$OMP private(LL,k,l,ib,iib,err_sub,d1,d2)                    &
  !$OMP num_threads(nb_thread)

  CALL alloc_NParray(tab_l, [basis_SG%nDindB%ndim],'tab_l',name_sub)
  CALL alloc_NParray(tab_ib,[basis_SG%nDindB%ndim],'tab_ib',name_sub)
  CALL alloc_NParray(tab_nb,[basis_SG%nDindB%ndim],'tab_nb',name_sub)

  !--------------------------------------------------------------
  !-- For the initialization of tab_l(:) and the use of ADD_ONE_TO_nDindex in the parallel loop
  ith = 0
  !$ ith = OMP_GET_THREAD_NUM()
  tab_l(:) = basis_SG%para_SGType2%nDval_init(:,ith+1)
  !--------------------------------------------------------------

  d1=basis_SG%para_SGType2%iG_th(ith+1)
  d2=basis_SG%para_SGType2%fG_th(ith+1)
  IF(openmpi) THEN
    d1=iGs_MPI(1,MPI_id)
    d2=iGs_MPI(2,MPI_id)
  ENDIF

  DO iG=d1,d2
!  ! we are not using the parallel do, to be able to use the correct initialized tab_l with nDval_init
!  DO iG=basis_SG%para_SGType2%iG_th(ith+1),basis_SG%para_SGType2%fG_th(ith+1)

    IF (debug) THEN
      write(out_unitp,*) '============================== iG,nb_G',iG,size(basis_SG%para_SGType2%tab_nb_OF_SRep)
      write(out_unitp,*) '======================= size(DPB)',basis_SG%para_SGType2%tab_nb_OF_SRep(iG)
      CALL flush_perso(out_unitp)
    END IF

  !  IF (max(1,size(basis_SG%para_SGType2%tab_nb_OF_SRep)/100) == 0 .OR. &
  !      size(basis_SG%para_SGType2%tab_nb_OF_SRep) < 100) THEN
  !    IF(MPI_id==0) write(out_unitp,*) 'iG,nb_G',iG,size(basis_SG%para_SGType2%tab_nb_OF_SRep)
  !    CALL flush_perso(out_unitp)
  !  END IF

    CALL ADD_ONE_TO_nDindex(basis_SG%para_SGType2%nDind_SmolyakRep,tab_l,iG=iG)

    tab_nb(:) = getbis_tab_nb(tab_l,basis_SG%tab_basisPrimSG)

    IF (iG == 1) THEN
      iBSRep = 0
    ELSE
      iBSRep = basis_SG%para_SGType2%tab_Sum_nb_OF_SRep(iG-1)
    END IF

    nDI = 0 ! old value
    tab_ib(:) = 1 ; tab_ib(1) = 0
    DO iBDP=1,basis_SG%para_SGType2%tab_nb_OF_SRep(iG)

      iBSRep = iBSRep + 1

      CALL ADD_ONE_TO_nDval_m1(tab_ib,tab_nb)

      LL = calc_L_OF_nDval(tab_ib,basis_SG%nDindB)

      IF (debug) THEN
        write(out_unitp,*) 'iG,iBDP,',iG,iBDP,'tab',tab_ib,'LL',LL
        CALL flush_perso(out_unitp)
      END IF
      IF (LL > basis_SG%nDindB%Lmax) CYCLE

      ! first estimation of nDI. It works most of the time !!!
      nDI    = 1
      LL     = 0
      DO k=1,basis_SG%nDindB%ndim

        ib = tab_ib(k)
        DO iib=1,ib-1
          l   = LL  + basis_SG%nDindB%tab_i_TO_l(k)%vec(iib)
          nDI = nDI + MaxnD_with_id_and_L(k,l)
        END DO
        LL = LL + basis_SG%nDindB%tab_i_TO_l(k)%vec(ib)
      END DO

      ! Then calculation of nDI with calc_nDI with the estimated nDI
      CALL calc_nDI(nDI,tab_ib,basis_SG%nDindB,err_sub)

      IF (err_sub == 0) THEN
        IF (debug) THEN
          write(out_unitp,*) '     nDI',nDI,nDI,'tab',tab_ib(:)
          CALL flush_perso(out_unitp)
        END IF
        IF (allocated(basis_SG%nDindB%Tab_nDval)) THEN
        IF ( .NOT. all(tab_ib == basis_SG%nDindB%Tab_nDval(:,nDI)) ) THEN
           write(out_unitp,*) 'ERROR in ',name_sub
           write(out_unitp,*) 'The Tab_ib(:) from both representation are different!'
           write(out_unitp,*) 'Tab_ib from Smolyak Rep',tab_ib
           write(out_unitp,*) 'Tab_ib from nDindB     ',basis_SG%nDindB%Tab_nDval(:,nDI)
           STOP
        END IF
        END IF
      END IF

      !write(out_unitp,*) 'nDI,Max_nDI', nDI,basis_SG%nDindB%Max_nDI ; flush(out_unitp)

      IF (nDI > 0 .AND. nDI <= basis_SG%nDindB%Max_nDI) THEN
        basis_SG%para_SGType2%tab_iB_OF_SRep_TO_iB(iBSRep) = nDI
        !$OMP ATOMIC
        Tab_inD_nDindB(nDI) = Tab_inD_nDindB(nDI) + 1
      END IF
    END DO

    IF (mod(iG,max(1,int(basis_SG%nb_SG/10))) == 0 .AND. .NOT. debug) THEN
      write(out_unitp,'(a)',ADVANCE='no') '---'
      CALL flush_perso(out_unitp)
    END IF

  END DO
  CALL dealloc_NParray(tab_l, 'tab_l', name_sub)
  CALL dealloc_NParray(tab_ib,'tab_ib',name_sub)
  CALL dealloc_NParray(tab_nb,'tab_nb',name_sub)
  !$OMP   END PARALLEL

  IF(openmpi) CALL Mapping_table_MPI(basis_SG,Max_Srep)

  IF (.NOT. debug) THEN
    write(out_unitp,'(a)',ADVANCE='yes') '----]'
    CALL flush_perso(out_unitp)
  END IF

  IF(MPI_id==0) THEN
    write(out_unitp,*) 'count 0',count(basis_SG%para_SGType2%tab_iB_OF_SRep_TO_iB == 0)
    IF (count(basis_SG%para_SGType2%tab_iB_OF_SRep_TO_iB == 0) > 0) THEN
      write(out_unitp,*) 'WARNING in ',name_sub
      write(out_unitp,*) 'The Smolyak Basis has more basis function than the nD-Basis'
      write(out_unitp,*) ' Probably LB < LG'
    END IF
  ENDIF

  IF (count(Tab_inD_nDindB == 0) > 0 .AND. MPI_np<2) THEN
    write(out_unitp,*) 'ERROR in ',name_sub
    write(out_unitp,*) ' Propblem with the mapping!'
    DO nDI=1,basis_SG%nDindB%Max_nDI
      write(out_unitp,*) nDI,'Tab_ib from nDindB     ',basis_SG%nDindB%Tab_nDval(:,nDI),':',Tab_inD_nDindB(nDI)
    END DO
    STOP
  END IF

  IF (debug) THEN
    write(out_unitp,*) ' Mapping nDindB%Tab_nDval'
    DO nDI=1,basis_SG%nDindB%Max_nDI
      write(out_unitp,*) nDI,'Tab_inD_nDindB',Tab_inD_nDindB(nDI)
    END DO
    write(out_unitp,*) ' Mapping tab_iB_OF_SRep_TO_iB'
    DO iBSRep=1,size(basis_SG%para_SGType2%tab_iB_OF_SRep_TO_iB)
     write(out_unitp,*) iBSRep,'tab_iB_OF_SRep_TO_iB',basis_SG%para_SGType2%tab_iB_OF_SRep_TO_iB(iBSRep)
   END DO
  END IF

  CALL SYSTEM_CLOCK(COUNT=nb_ticks_final)
  nb_ticks = nb_ticks_final - nb_ticks_initial
  IF (nb_ticks_final < nb_ticks_initial) nb_ticks = nb_ticks + nb_ticks_max
  write(out_unitp,*) 'real time:',REAL(nb_ticks) / real(nb_ticks_sec,kind=Rkind)

  write(out_unitp,*) 'END ',name_sub
  CALL flush_perso(out_unitp)

END SUBROUTINE Set_tables_FOR_SmolyakRepBasis_TO_tabPackedBasis

SUBROUTINE SmolyakRepBasis_TO_tabPackedBasis(SRep,tabR,tab_ba,nDindB,SGType2,WeightSG)
USE mod_system
USE mod_param_SGType2
USE mod_nDindex
USE mod_MPI_aux
IMPLICIT NONE

TYPE(Type_SmolyakRep),           intent(in)        :: SRep
real(kind=Rkind),                intent(inout)     :: tabR(:)
TYPE (Type_nDindex),             intent(in)        :: nDindB
TYPE(basis),                     intent(in)        :: tab_ba(0:,:) ! tab_ba(0:L,D)
TYPE (param_SGType2),            intent(in)        :: SGType2
real(kind=Rkind),                intent(in)        :: WeightSG(:)

integer :: iBSRep,iG,iB,nb_BG,nR,itabR,nDI
integer :: ib0,nb_AT_iG,iB_ib0,nDI_ib0
integer :: d1,d2

!write(out_unitp,*) 'SmolyakRepBasis_TO_tabPackedBasis'
!CALL write_SmolyakRep(SRep)
flush(out_unitp)

IF (size(tabR) == 0) STOP ' ERROR in SmolyakRepBasis_TO_tabPackedBasis. size(tabR)=0'
!write(out_unitp,*) 'nb_size',Size_SmolyakRep(SRep) ; flush(out_unitp)

tabR(:) = ZERO

d1=lbound(SRep%SmolyakRep,dim=1)
d2=ubound(SRep%SmolyakRep,dim=1)
! IF(openmpi .AND. MPI_scheme==1) THEN
!   d1=iGs_MPI(1,MPI_id)
!   d2=iGs_MPI(2,MPI_id)
! ENDIF
IF(openmpi) THEN
  IF(MPI_scheme==1) THEN
    d1=iGs_MPI(1,MPI_id)
    d2=iGs_MPI(2,MPI_id)
  ELSEIF(MPI_scheme==3) THEN
    IF(MPI_nodes_p0) THEN
      d1=iGs_MPI(1,MPI_id)
      d2=iGs_MPI(2,MPI_sub_id(2))
    ELSE
      d1=0
      d2=-1
    ENDIF
  ENDIF
ENDIF

DO ib0=1,SRep%nb0
  iBSRep = 0
  !DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)
  DO iG=d1,d2

    IF (abs(WeightSG(iG)) < ONETENTH**6) CYCLE

    iBSRep = SGType2%tab_Sum_nb_OF_SRep(iG)-SGType2%tab_nb_OF_SRep(iG)
    nb_AT_iG = SGType2%tab_nb_OF_SRep(iG)

    DO iB_ib0=1,nb_AT_iG
      iBSRep = iBSRep + 1

      nDI_ib0 = SGType2%tab_iB_OF_SRep_TO_iB(iBSRep)

      IF (nDI_ib0 > 0 .AND. nDI_ib0 <= nDindB%Max_nDI) THEN
        nDI = (ib0-1)*nDindB%Max_nDI + nDI_ib0
        iB  = (ib0-1)*nb_AT_iG       + iB_ib0
        tabR(nDI) = tabR(nDI) + WeightSG(iG)*SRep%SmolyakRep(iG)%V(iB)
      END IF
    END DO

  END DO
END DO

!write(out_unitp,*) 'tabR',tabR
!write(out_unitp,*) 'END SmolyakRepBasis_TO_tabPackedBasis'


END SUBROUTINE SmolyakRepBasis_TO_tabPackedBasis

! warning, for MPI scheme 1, the "SRep" are MPI calculated on different processor
! thus not shared by all processor. Be careful with this.
SUBROUTINE tabPackedBasis_TO_SmolyakRepBasis(SRep,tabR,tab_ba,nDindB,SGType2)
USE mod_system
USE mod_basis_set_alloc
USE mod_param_SGType2
USE mod_nDindex
USE mod_MPI_aux
IMPLICIT NONE

TYPE(Type_SmolyakRep),           intent(inout)     :: SRep
real(kind=Rkind),                intent(in)        :: tabR(:)
TYPE(basis),                     intent(in)        :: tab_ba(0:,:) ! tab_ba(0:L,D)
TYPE (Type_nDindex),             intent(in)        :: nDindB
TYPE (param_SGType2),            intent(in)        :: SGType2

integer :: iBSRep,iG,iB,nb_BG,nR,itabR,nDI
integer :: ib0,nb_AT_iG,iB_ib0,nDI_ib0
integer :: d1,d2

!write(out_unitp,*) 'BEGINNING tabPackedBasis_TO_SmolyakRepBasis' ; flush(out_unitp)

CALL alloc2_SmolyakRep(SRep,SGType2%nDind_SmolyakRep,tab_ba,grid=.FALSE.,nb0=SGType2%nb0)
SRep = ZERO
!write(out_unitp,*) 'nb_size',Size_SmolyakRep(SRep)

d1=lbound(SRep%SmolyakRep,dim=1)
d2=ubound(SRep%SmolyakRep,dim=1)
! IF(openmpi .AND. MPI_scheme==1) THEN
!   d1=iGs_MPI(1,MPI_id)
!   d2=iGs_MPI(2,MPI_id)
! ENDIF

IF(openmpi) THEN
  IF(MPI_scheme==1) THEN
    d1=iGs_MPI(1,MPI_id)
    d2=iGs_MPI(2,MPI_id)
  ELSEIF(MPI_scheme==3) THEN
    IF(MPI_nodes_p0) THEN
      d1=iGs_MPI(1,MPI_id)
      d2=iGs_MPI(2,MPI_sub_id(2))
    ELSE
      d1=0
      d2=-1
    ENDIF
  ENDIF
ENDIF

DO ib0=1,SRep%nb0
  iBSRep = 0
  !DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)
  DO iG=d1,d2

    iBSRep = SGType2%tab_Sum_nb_OF_SRep(iG)-SGType2%tab_nb_OF_SRep(iG)
    nb_AT_iG = SGType2%tab_nb_OF_SRep(iG)

    DO iB_ib0=1,nb_AT_iG
      iBSRep = iBSRep + 1

      nDI_ib0 = SGType2%tab_iB_OF_SRep_TO_iB(iBSRep)

      !write(out_unitp,*) 'iG,iB,nDI',iG,iB,nDI ; flush(out_unitp)
      IF (nDI_ib0 > 0 .AND. nDI_ib0 <= nDindB%Max_nDI) THEN

        nDI = (ib0-1)*nDindB%Max_nDI + nDI_ib0
        iB  = (ib0-1)*nb_AT_iG       + iB_ib0

        SRep%SmolyakRep(iG)%V(iB) = tabR(nDI)
      END IF
    END DO

  END DO
END DO
!write(out_unitp,*) 'END tabPackedBasis_TO_SmolyakRepBasis' ; flush(out_unitp)


END SUBROUTINE tabPackedBasis_TO_SmolyakRepBasis
SUBROUTINE tabPackedBasis_TO2_tabR_AT_iG(tabR_iG,tabR,iG,tab_l,tab_ba,SGType2)
USE mod_system
USE mod_basis_set_alloc
USE mod_param_SGType2
USE mod_nDindex
IMPLICIT NONE

real(kind=Rkind),     allocatable,  intent(inout)     :: tabR_iG(:)
real(kind=Rkind),                   intent(in)        :: tabR(:)
integer,                            intent(in)        :: iG,tab_l(:)
TYPE (param_SGType2),               intent(in)        :: SGType2
TYPE(basis),                        intent(in)        :: tab_ba(0:,:) ! tab_ba(0:L,D)

integer               :: iBSRep,iB,nDI,i
integer               :: tab_ldelta(size(tab_ba,dim=2))
integer               :: tab_ib(size(tab_ba,dim=2))
integer               :: tab_nb(size(tab_ba,dim=2))
integer :: ib0,nb_AT_iG,iB_ib0,nDI_ib0,Max_nDI_ib0

  IF (allocated(tabR_iG)) THEN
    CALL dealloc_NParray(tabR_iG,'tabR_iG','tabPackedBasis_TO_tabR_AT_iG')
  END IF

CALL alloc_NParray(tabR_iG,[SGType2%tab_nb_OF_SRep(iG)*SGType2%nb0], &
                    'tabR_iG','tabPackedBasis_TO_tabR_AT_iG')

nb_AT_iG = SGType2%tab_nb_OF_SRep(iG)
Max_nDI_ib0 = size(tabR)/SGType2%nb0

tabR_iG(:) = ZERO
tab_nb(:)  = getbis_tab_nb(tab_l,tab_ba)


DO ib0=1,SGType2%nb0
  iBSRep = SGType2%tab_Sum_nb_OF_SRep(iG)-SGType2%tab_nb_OF_SRep(iG)

  tab_ib(:) = 1 ; tab_ib(1) = 0
  DO iB_ib0=1,SGType2%tab_nb_OF_SRep(iG)
    CALL ADD_ONE_TO_nDval_m1(tab_ib,tab_nb)

    !????????????????????????????????????????????, to be removed !!!
    DO i=1,size(tab_ba,dim=2)
      tab_ldelta(i) = tab_ba(tab_l(i),i)%L_TO_nb%Tab_n_TO_L(tab_ib(i))
    END DO
    ! in tab_ldelta(:), we have the l(:) for the DP in DeltaB
    !????????????????????????????????????????????,


    iBSRep     = iBSRep + 1
    nDI_ib0    = SGType2%tab_iB_OF_SRep_TO_iB(iBSRep)
    !write(out_unitp,*) 'iG,iB,nDI',iG,iB,nDI ; flush(out_unitp)
    IF (nDI_ib0 > 0 .AND. nDI_ib0 <= Max_nDI_ib0) THEN
      nDI = (ib0-1)*Max_nDI_ib0                + nDI_ib0
      iB  = (ib0-1)*SGType2%tab_nb_OF_SRep(iG) + iB_ib0
      tabR_iG(iB) = tabR(nDI)
    END IF
  END DO
END DO

END SUBROUTINE tabPackedBasis_TO2_tabR_AT_iG

!=======================================================================================
!> in sub_TabOpPsi_FOR_SGtype4:
!>  CALL tabPackedBasis_TO_tabR_AT_iG(PsiR(itab)%V,psi(itab)%RvecB,iG,BasisnD%para_SGType2)
!>  PsiR(itab)%V     --> tabR_iG
!>  psi(itab)%RvecB  --> tabR
!
!> tabR_iG is assigned in the subroutine according to iG using the informtation in tabR
!=======================================================================================
SUBROUTINE tabPackedBasis_TO_tabR_AT_iG(tabR_iG,tabR,iG,SGType2)
USE mod_system
USE mod_basis_set_alloc
USE mod_param_SGType2
USE mod_nDindex
IMPLICIT NONE

real(kind=Rkind),     allocatable,  intent(inout)     :: tabR_iG(:)
real(kind=Rkind),                   intent(in)        :: tabR(:)
integer,                            intent(in)        :: iG
TYPE (param_SGType2),               intent(in)        :: SGType2
!TYPE(basis),                        intent(in)        :: tab_ba(0:,:) ! tab_ba(0:L,D)

integer :: iBSRep,iB,nDI
integer :: ib0,nb_AT_iG,iB_ib0,nDI_ib0,Max_nDI_ib0

 !----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='tabPackedBasis_TO_tabR_AT_iG'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'tabR',tabR(:)
        write(out_unitp,*) 'nb0',SGType2%nb0
        write(out_unitp,*) 'iG',iG
        write(out_unitp,*) 'nb at iG',SGType2%tab_nb_OF_SRep(iG)
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------


IF (allocated(tabR_iG)) THEN
  CALL dealloc_NParray(tabR_iG,'tabR_iG',name_sub)
END IF

CALL alloc_NParray(tabR_iG,[SGType2%tab_nb_OF_SRep(iG)*SGType2%nb0],&
                  'tabR_iG',name_sub)

tabR_iG(:)  = ZERO
Max_nDI_ib0 = size(tabR)/SGType2%nb0
!write(out_unitp,*) 'iG,nb0,nb,Max_nDI_ib0',iG,SGType2%nb0,SGType2%tab_nb_OF_SRep(iG),Max_nDI_ib0 ; flush(out_unitp)

DO ib0=1,SGType2%nb0
  iBSRep = SGType2%tab_Sum_nb_OF_SRep(iG)-SGType2%tab_nb_OF_SRep(iG)

  DO iB_ib0=1,SGType2%tab_nb_OF_SRep(iG)
    iBSRep  = iBSRep + 1
    nDI_ib0 = SGType2%tab_iB_OF_SRep_TO_iB(iBSRep)
    !write(out_unitp,*) 'ib0,iB_ib0',ib0,iB_ib0 ; flush(out_unitp)
    !write(out_unitp,*) 'iBSRep,nDI_ib0',iBSRep,nDI_ib0 ; flush(out_unitp)

    IF (nDI_ib0 > 0 .AND. nDI_ib0 <= Max_nDI_ib0 ) THEN
      nDI = (ib0-1)*Max_nDI_ib0                + nDI_ib0
      iB  = (ib0-1)*SGType2%tab_nb_OF_SRep(iG) + iB_ib0
      !write(out_unitp,*) 'iB,nDI',iB,nDI ; flush(out_unitp)

      tabR_iG(iB) = tabR(nDI)
    END IF

  END DO
END DO
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'tabR_iG: ',tabR_iG(:)
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------
END SUBROUTINE tabPackedBasis_TO_tabR_AT_iG
!=======================================================================================

!=======================================================================================
SUBROUTINE tabR_AT_iG_TO_tabPackedBasis(tabR,tabR_iG,iG,SGType2,WeightiG)
USE mod_system
USE mod_basis_set_alloc
USE mod_param_SGType2
USE mod_nDindex
IMPLICIT NONE

real(kind=Rkind),                   intent(in)        :: tabR_iG(:)
real(kind=Rkind),                   intent(inout)     :: tabR(:)
integer,                            intent(in)        :: iG
!TYPE(basis),                        intent(in)        :: tab_ba(0:,:) ! tab_ba(0:L,D)
TYPE (param_SGType2),               intent(in)        :: SGType2
real(kind=Rkind),                   intent(in)        :: WeightiG


integer :: iBSRep,iB,nDI
integer :: ib0,nb_AT_iG,iB_ib0,nDI_ib0,Max_nDI_ib0


IF (size(tabR) == 0) STOP ' ERROR in SmolyakRepBasis_TO_tabPackedBasis. size(tabR)=0'

Max_nDI_ib0 = size(tabR)/SGType2%nb0

DO ib0=1,SGType2%nb0
  iBSRep = SGType2%tab_Sum_nb_OF_SRep(iG)-SGType2%tab_nb_OF_SRep(iG)

  DO iB_ib0=1,SGType2%tab_nb_OF_SRep(iG)
    iBSRep = iBSRep + 1
    nDI_ib0    = SGType2%tab_iB_OF_SRep_TO_iB(iBSRep)
    !write(out_unitp,*) 'iG,iB,nDI',iG,iB,nDI ; flush(out_unitp)
    IF (nDI_ib0 > 0 .AND. nDI_ib0 <= Max_nDI_ib0 ) THEN
      nDI = (ib0-1)*Max_nDI_ib0                + nDI_ib0
      iB  = (ib0-1)*SGType2%tab_nb_OF_SRep(iG) + iB_ib0
      !$OMP ATOMIC
      tabR(nDI) = tabR(nDI) + WeightiG*tabR_iG(iB)
    END IF
  END DO
END DO

END SUBROUTINE tabR_AT_iG_TO_tabPackedBasis
!=======================================================================================

SUBROUTINE SmolyakRep2_TO_tabR1(tabR1,SRep2)
USE mod_system
IMPLICIT NONE

real(kind=Rkind), allocatable,   intent(inout)  :: tabR1(:)
CLASS (Type_SmolyakRep),         intent(in)     :: SRep2

real(kind=Rkind), allocatable  :: tab(:,:)
real(kind=Rkind), allocatable  :: V(:,:)

integer               :: iG,nb_BG,itabR
integer               :: ib0,nb,nb_AT_iG


IF (allocated(tabR1)) CALL dealloc_NParray(tabR1,'tabR1','SmolyakRep2_TO_tabR1')

nb_BG = Size_SmolyakRep(SRep2)

nb = nb_BG/SRep2%nb0
CALL alloc_NParray(tab,[nb,SRep2%nb0],                              &
                  'tab','SmolyakRep2_TO_tabR1')

IF (nb_BG > 0) THEN

  itabR = 0
  DO iG=lbound(SRep2%SmolyakRep,dim=1),ubound(SRep2%SmolyakRep,dim=1)
    nb_AT_iG = size(SRep2%SmolyakRep(iG)%V) / SRep2%nb0

    V = reshape(SRep2%SmolyakRep(iG)%V,shape=[nb_AT_iG,SRep2%nb0] )

    tab(itabR+1:itabR+nb_AT_iG,:) = V(:,:)

    itabR = itabR + nb_AT_iG

  END DO
  CALL alloc_NParray(tabR1,[nb_BG],'tabR1','SmolyakRep2_TO_tabR1')
  tabR1(:) = reshape(tab,shape=[nb_BG] )

  CALL dealloc_NParray(tab,'tab','SmolyakRep2_TO_tabR1bis')
  deallocate(V)

END IF

END SUBROUTINE SmolyakRep2_TO_tabR1
SUBROUTINE SmolyakRep2_TO_tabR1bis(tabR1,SRep2)
USE mod_system
IMPLICIT NONE

real(kind=Rkind),                intent(inout)  :: tabR1(:)
TYPE(Type_SmolyakRep),           intent(in)     :: SRep2


real(kind=Rkind), allocatable  :: tab(:,:)
real(kind=Rkind), allocatable  :: V(:,:)

integer               :: iG,nb_BG,itabR
integer               :: ib0,nb,nb_AT_iG

nb_BG = Size_SmolyakRep(SRep2)
IF (nb_BG /= size(tabR1)) THEN
  write(out_unitp,*) ' ERROR in SmolyakRep2_TO_tabR1bis'
  write(out_unitp,*) ' size of tabR1 and SRep2 are different',size(tabR1),nb_BG
  STOP
END IF

nb = nb_BG/SRep2%nb0
CALL alloc_NParray(tab,[nb,SRep2%nb0],                              &
                  'tab','SmolyakRep2_TO_tabR1bis')

IF (nb_BG > 0) THEN

  itabR = 0
  DO iG=lbound(SRep2%SmolyakRep,dim=1),ubound(SRep2%SmolyakRep,dim=1)
    nb_AT_iG = size(SRep2%SmolyakRep(iG)%V) / SRep2%nb0

    V = reshape(SRep2%SmolyakRep(iG)%V,shape=[nb_AT_iG,SRep2%nb0] )

    tab(itabR+1:itabR+nb_AT_iG,:) = V(:,:)

    itabR = itabR + nb_AT_iG

  END DO
  tabR1 = reshape(tab,shape=[nb_BG] )

  CALL dealloc_NParray(tab,'tab','SmolyakRep2_TO_tabR1bis')
  deallocate(V)
END IF

END SUBROUTINE SmolyakRep2_TO_tabR1bis
SUBROUTINE SmolyakRepC2_TO_tabC1bis(tabC1,SRep2)
USE mod_system
IMPLICIT NONE

complex(kind=Rkind),              intent(inout)  :: tabC1(:)
TYPE(Type_SmolyakRepC),           intent(in)     :: SRep2

integer               :: iG,nb_BG,nR,itabR

nb_BG = Size_SmolyakRepC(SRep2)
IF (nb_BG /= size(tabC1)) THEN
  write(out_unitp,*) ' ERROR in SmolyakRepC2_TO_tabC1bis'
  write(out_unitp,*) ' size of tabR1 and SRep2 are different',size(tabC1),nb_BG
  STOP
END IF

IF (nb_BG > 0) THEN
  itabR = 0
  DO iG=lbound(SRep2%SmolyakRep,dim=1),ubound(SRep2%SmolyakRep,dim=1)
    nR = size(SRep2%SmolyakRep(iG)%V)
    tabC1(itabR+1:itabR+nR) = SRep2%SmolyakRep(iG)%V
    itabR = itabR + nR
  END DO
END IF

END SUBROUTINE SmolyakRepC2_TO_tabC1bis
SUBROUTINE tabR2_TO_SmolyakRep1(SRep1,tabR2)
USE mod_system
IMPLICIT NONE

real(kind=Rkind),                intent(in)     :: tabR2(:)
CLASS (Type_SmolyakRep),         intent(inout)  :: SRep1

real(kind=Rkind), allocatable  :: tab(:,:)
real(kind=Rkind), allocatable  :: V(:,:)

integer               :: iG,nb_BG,itabR
integer               :: ib0,nb,nb_AT_iG

nb_BG = Size_SmolyakRep(SRep1)
IF (size(tabR2) /= nb_BG) THEN
  write(out_unitp,*) ' ERROR in tabR2_TO_SmolyakRep1'
  write(out_unitp,*) ' sizes are different!!'
  write(out_unitp,*) ' sizes of tabR2 and SRep1',size(tabR2),nb_BG
  STOP
END IF

IF (nb_BG > 0) THEN

  nb = nb_BG/SRep1%nb0

  tab = reshape(tabR2,shape=[nb,SRep1%nb0] )

  itabR = 0
  DO iG=lbound(SRep1%SmolyakRep,dim=1),ubound(SRep1%SmolyakRep,dim=1)
    nb_AT_iG = size(SRep1%SmolyakRep(iG)%V) / SRep1%nb0

    V = tab(itabR+1:itabR+nb_AT_iG,:)

    SRep1%SmolyakRep(iG)%V(:) = reshape(V,shape(SRep1%SmolyakRep(iG)%V))

    itabR = itabR + nb_AT_iG

  END DO

  deallocate(tab)
  deallocate(V)

END IF

END SUBROUTINE tabR2_TO_SmolyakRep1
SUBROUTINE tabR2bis_TO_SmolyakRep1(SRep1,tabR2)
USE mod_system
IMPLICIT NONE

real(kind=Rkind),                intent(in)     :: tabR2(:)
TYPE(Type_SmolyakRep),           intent(inout)  :: SRep1

real(kind=Rkind), allocatable  :: tab(:,:)
real(kind=Rkind), allocatable  :: V(:,:)

integer               :: iG,nb_BG,itabR
integer               :: ib0,nb,nb_AT_iG

nb_BG = Size_SmolyakRep(SRep1)
IF (size(tabR2) /= nb_BG) THEN
  write(out_unitp,*) ' ERROR in tabR2bis_TO_SmolyakRep1'
  write(out_unitp,*) ' sizes are different!!'
  write(out_unitp,*) ' sizes of tabR2 and SRep1',size(tabR2),nb_BG
  STOP
END IF

IF (nb_BG > 0) THEN

  nb = nb_BG/SRep1%nb0

  tab = reshape(tabR2,shape=[nb,SRep1%nb0] )

  itabR = 0
  DO iG=lbound(SRep1%SmolyakRep,dim=1),ubound(SRep1%SmolyakRep,dim=1)
    nb_AT_iG = size(SRep1%SmolyakRep(iG)%V) / SRep1%nb0

    V = tab(itabR+1:itabR+nb_AT_iG,:)

    SRep1%SmolyakRep(iG)%V(:) = reshape(V,shape(SRep1%SmolyakRep(iG)%V))

    itabR = itabR + nb_AT_iG

  END DO

  deallocate(tab)
  deallocate(V)

END IF

END SUBROUTINE tabR2bis_TO_SmolyakRep1
SUBROUTINE tabC2bis_TO_SmolyakRepC1(SRep1,tabC2)
USE mod_system
IMPLICIT NONE

complex(kind=Rkind),              intent(in)     :: tabC2(:)
TYPE(Type_SmolyakRepC),           intent(inout)  :: SRep1

complex(kind=Rkind), allocatable  :: tab(:,:)
complex(kind=Rkind), allocatable  :: V(:,:)

integer               :: iG,nb_BG,itabR,nR
integer               :: ib0,nb,nb_AT_iG

nb_BG = Size_SmolyakRepC(SRep1)
IF (size(tabC2) /= nb_BG) THEN
  write(out_unitp,*) ' ERROR in tabC2bis_TO_SmolyakRepC1'
  write(out_unitp,*) ' sizes are different!!'
  write(out_unitp,*) ' sizes of tabR2 and SRep1',size(tabC2),nb_BG
  STOP
END IF

IF (nb_BG > 0) THEN

  itabR = 0
  DO iG=lbound(SRep1%SmolyakRep,dim=1),ubound(SRep1%SmolyakRep,dim=1)
    nR = size(SRep1%SmolyakRep(iG)%V)
    SRep1%SmolyakRep(iG)%V(:) = tabC2(itabR+1:itabR+nR)
    itabR = itabR + nR
  END DO

END IF

END SUBROUTINE tabC2bis_TO_SmolyakRepC1
SUBROUTINE tabR2grid_TO_tabR1_AT_iG(tabR1,tabR2,iG,SGType2)
USE mod_system
USE mod_param_SGType2
IMPLICIT NONE

real(kind=Rkind), allocatable,   intent(inout)  :: tabR1(:)
real(kind=Rkind),                intent(in)     :: tabR2(:)
integer,                         intent(in)     :: iG
TYPE (param_SGType2),            intent(in)     :: SGType2

integer               :: nR,itabR


IF (size(tabR2) > 0) THEN

  itabR = SGType2%tab_Sum_nq_OF_SRep(iG)
  nR    = SGType2%tab_nq_OF_SRep(iG)

  IF (allocated(tabR1)) CALL dealloc_NParray(tabR1,'tabR1','tabR2grid_TO_tabR1_AT_iG')
  CALL alloc_NParray(tabR1,[nR],'tabR1','tabR2grid_TO_tabR1_AT_iG')


  tabR1(:) = tabR2(itabR-nR+1:itabR)


END IF

END SUBROUTINE tabR2grid_TO_tabR1_AT_iG

!=======================================================================================
! tabR2 ==> tabR1
! tabR2: grid  (:,j,i)
! tabR1: GridOp(:,j,i,iterm)
!---------------------------------------------------------------------------------------
SUBROUTINE tabR2gridbis_TO_tabR1_AT_iG(tabR1,tabR2,iG,SGType2)
USE mod_system
USE mod_param_SGType2
IMPLICIT NONE

real(kind=Rkind),                intent(inout)  :: tabR1(:)
real(kind=Rkind),                intent(in)     :: tabR2(:)
integer,                         intent(in)     :: iG
TYPE (param_SGType2),            intent(in)     :: SGType2

integer               :: nR,itabR


IF (size(tabR2) > 0) THEN

  itabR = SGType2%tab_Sum_nq_OF_SRep(iG)
  nR    = SGType2%tab_nq_OF_SRep(iG)
  tabR1(:) = tabR2(itabR-nR+1:itabR)
END IF

END SUBROUTINE tabR2gridbis_TO_tabR1_AT_iG
!=======================================================================================

SUBROUTINE tabR1_AT_iG_TO_tabR2grid(tabR1,tabR2,iG,SGType2)
USE mod_system
USE mod_param_SGType2
IMPLICIT NONE

real(kind=Rkind),                intent(in)     :: tabR1(:)
real(kind=Rkind), allocatable,   intent(inout)  :: tabR2(:)
integer,                         intent(in)     :: iG
TYPE (param_SGType2),            intent(in)     :: SGType2

integer               :: nR,itabR


IF (size(tabR1) > 0) THEN

  itabR = SGType2%tab_Sum_nq_OF_SRep(iG)
  nR    = SGType2%tab_nq_OF_SRep(iG)

  IF (.NOT. allocated(tabR2)) THEN
    CALL alloc_NParray(tabR2,[SGType2%tab_nq_OF_SRep(SGType2%nb_SG)],&
                      'tabR2','tabR1_AT_iG_TO_tabR2grid')
  END IF


  tabR2(itabR-nR+1:itabR) = tabR1(:)


END IF

END SUBROUTINE tabR1_AT_iG_TO_tabR2grid

SUBROUTINE R2_TO_SmolyakRep1(SRep1,R2)
USE mod_system
IMPLICIT NONE

real(kind=Rkind),                intent(in)     :: R2
CLASS (Type_SmolyakRep),         intent(inout)  :: SRep1

integer               :: iG,nb_BG


nb_BG = Size_SmolyakRep(SRep1)

IF (nb_BG > 0) THEN

  DO iG=lbound(SRep1%SmolyakRep,dim=1),ubound(SRep1%SmolyakRep,dim=1)
    SRep1%SmolyakRep(iG)%V(:) = R2
  END DO

END IF

END SUBROUTINE R2_TO_SmolyakRep1

SUBROUTINE R2_TO_SmolyakRep1_with_tab_i(SRep1,R2,tab_i,tab_ind,tab_ba)
USE mod_system
USE mod_basis_set_alloc
IMPLICIT NONE

real(kind=Rkind),                intent(in)             :: R2
TYPE(Type_SmolyakRep),           intent(inout)          :: SRep1
integer,         allocatable,    intent(in)             :: tab_ind(:,:) ! tab_ind(D,MaxnD)
integer,                         intent(in)             :: tab_i(:)
TYPE(basis),                     intent(in)             :: tab_ba(0:,:) ! tab_ba(0:L,D)

integer               :: i,II,iG,nb_BG,D
integer, allocatable  :: tab_n(:)


nb_BG = Size_SmolyakRep(SRep1)

IF (SRep1%Grid) STOP 'Grid is not possible in R2_TO_SmolyakRep1_with_tab_i'

IF (nb_BG > 0) THEN
IF (SRep1%delta) THEN
  STOP 'SRep1%delta=t impossible yet'
ELSE
  D = size(tab_ind(:,1))
  DO iG=lbound(SRep1%SmolyakRep,dim=1),ubound(SRep1%SmolyakRep,dim=1)

    tab_n = get_tab_nb(tab_ind(:,iG),tab_ba)

    IF (minval(tab_n - tab_i) < 0) CYCLE

    II = tab_i(D)
    DO i=D-1,1,-1
      II = (II-1)*tab_n(i) + tab_i(i)
    END DO

    IF (II > size(SRep1%SmolyakRep(iG)%V)) STOP 'ERROR in R2_TO_SmolyakRep1_with_tab_i'

    SRep1%SmolyakRep(iG)%V(II) = R2
  END DO
  IF (allocated(tab_n)) CALL dealloc_NParray(tab_n,'tab_n','R2_TO_SmolyakRep1_with_tab_i')
END IF
END IF

END SUBROUTINE R2_TO_SmolyakRep1_with_tab_i

FUNCTION dot_product_SmolyakRep(SRep1,SRep2,WSRep) RESULT(R)
USE mod_system
IMPLICIT NONE

real(kind=Rkind)  :: R
TYPE(Type_SmolyakRep),           intent(in)     :: SRep1,SRep2
real(kind=Rkind),                intent(in)     :: WSRep(:)


integer               :: iG,nb_BG

nb_BG = Size_SmolyakRep(SRep1)
IF (nb_BG /= Size_SmolyakRep(SRep2) .AND. nb_BG /= size(WSRep)) THEN
  write(out_unitp,*) 'ERROR in dot_product_SmolyakRep'
  write(out_unitp,*) 'sizes are different',Size_SmolyakRep(SRep1),Size_SmolyakRep(SRep2),size(WSRep)
  STOP
END IF

R = ZERO
IF (nb_BG > 0) THEN

  DO iG=lbound(SRep1%SmolyakRep,dim=1),ubound(SRep1%SmolyakRep,dim=1)
    R = R + WSRep(iG) * dot_product(SRep1%SmolyakRep(iG)%V,SRep2%SmolyakRep(iG)%V)
  END DO

END IF

END FUNCTION dot_product_SmolyakRep

FUNCTION dot_product_SmolyakRep_Grid(SRep1,SRep2,SRep_w,WSRep,ib0) RESULT(R)
USE mod_system
IMPLICIT NONE

real(kind=Rkind)  :: R
TYPE(Type_SmolyakRep),           intent(in)            :: SRep1,SRep2
TYPE(Type_SmolyakRep),           intent(in)            :: SRep_w
real(kind=Rkind),                intent(in)            :: WSRep(:)
integer,                         intent(in), optional  :: ib0


integer               :: iG,nb_BG,nq_AT_iG,i1,i2

!write(out_unitp,*) 'dot_product_SmolyakRep_Grid: size',Size_SmolyakRep(SRep1),Size_SmolyakRep(SRep2),Size_SmolyakRep(SRep_w)
!write(out_unitp,*) 'nb DP terms',size(SRep1%SmolyakRep),size(SRep2%SmolyakRep),size(WSRep),size(WSRep)

!flush(out_unitp)

nb_BG = Size_SmolyakRep(SRep1)
IF (nb_BG /= Size_SmolyakRep(SRep2) .OR. size(SRep1%SmolyakRep) /= size(WSRep) .OR. size(SRep2%SmolyakRep) /= size(WSRep)) THEN
  write(out_unitp,*) 'ERROR in dot_product_SmolyakRep'
  write(out_unitp,*) 'sizes are different',Size_SmolyakRep(SRep1),Size_SmolyakRep(SRep2)
  write(out_unitp,*) 'nb DP terms',size(SRep1%SmolyakRep),size(SRep2%SmolyakRep),size(WSRep)
  STOP
END IF

R = ZERO
IF (nb_BG > 0) THEN
  IF (present(ib0)) THEN
    DO iG=lbound(SRep1%SmolyakRep,dim=1),ubound(SRep1%SmolyakRep,dim=1)
      nq_AT_iG = size(SRep1%SmolyakRep(iG)%V) / SRep1%nb0
      i1 = 1+(ib0-1)*nq_AT_iG
      i2 = ib0*nq_AT_iG
      R = R + WSRep(iG) * dot_product(SRep1%SmolyakRep(iG)%V(i1:i2),    &
                  SRep_w%SmolyakRep(iG)%V*SRep2%SmolyakRep(iG)%V(i1:i2))
    END DO
  ELSE
    DO iG=lbound(SRep1%SmolyakRep,dim=1),ubound(SRep1%SmolyakRep,dim=1)
      R = R + WSRep(iG) * dot_product(SRep1%SmolyakRep(iG)%V,           &
                         SRep_w%SmolyakRep(iG)%V*SRep2%SmolyakRep(iG)%V)
    END DO
  END IF

END IF

END FUNCTION dot_product_SmolyakRep_Grid

FUNCTION dot_product_SmolyakRep_Basis(SRep1,SRep2,WSRep,ib0) RESULT(R)
USE mod_system
USE mod_MPI_aux
IMPLICIT NONE

real(kind=Rkind)  :: R
TYPE(Type_SmolyakRep),           intent(in)            :: SRep1,SRep2
real(kind=Rkind),                intent(in)            :: WSRep(:)
integer,                         intent(in), optional  :: ib0


integer               :: iG,nb_BG,nb_AT_iG,i1,i2
integer               :: d1,d2

nb_BG = Size_SmolyakRep(SRep1)
IF (nb_BG /= Size_SmolyakRep(SRep2) .OR. size(SRep1%SmolyakRep) /= size(WSRep) .OR. size(SRep2%SmolyakRep) /= size(WSRep)) THEN
  write(out_unitp,*) 'ERROR in dot_product_SmolyakRep_Basis'
  write(out_unitp,*) 'sizes are different',Size_SmolyakRep(SRep1),Size_SmolyakRep(SRep2)
  write(out_unitp,*) 'nb DP terms are different',size(SRep1%SmolyakRep),size(SRep2%SmolyakRep),size(WSRep)
  STOP
END IF

!write(out_unitp,*) 'dot_product_SmolyakRep_Basis: size',Size_SmolyakRep(SRep1),Size_SmolyakRep(SRep2),size(WSRep)
!flush(out_unitp)

d1=lbound(SRep1%SmolyakRep,dim=1)
d2=ubound(SRep1%SmolyakRep,dim=1)
! IF(openmpi .AND. MPI_scheme==1) THEN
!   d1=iGs_MPI(1,MPI_id)
!   d2=iGs_MPI(2,MPI_id)
! ENDIF
IF(openmpi) THEN
  IF(MPI_scheme==1) THEN
    d1=iGs_MPI(1,MPI_id)
    d2=iGs_MPI(2,MPI_id)
  ELSEIF(MPI_scheme==3) THEN
    IF(MPI_nodes_p0) THEN
      d1=iGs_MPI(1,MPI_id)
      d2=iGs_MPI(2,MPI_sub_id(2))
    ELSE
      d1=0
      d2=-1
    ENDIF
  ENDIF
ENDIF

R = ZERO
IF (nb_BG > 0) THEN
  IF (present(ib0)) THEN
    !DO iG=lbound(SRep1%SmolyakRep,dim=1),ubound(SRep1%SmolyakRep,dim=1)
    DO iG=d1,d2
      nb_AT_iG = size(SRep1%SmolyakRep(iG)%V) / SRep1%nb0
      i1 = 1+(ib0-1)*nb_AT_iG
      i2 = ib0*nb_AT_iG
      R = R + WSRep(iG) * dot_product(SRep1%SmolyakRep(iG)%V(i1:i2),SRep2%SmolyakRep(iG)%V(i1:i2))
    END DO
  ELSE
    !DO iG=lbound(SRep1%SmolyakRep,dim=1),ubound(SRep1%SmolyakRep,dim=1)
    DO iG=d1,d2
      R = R + WSRep(iG) * dot_product(SRep1%SmolyakRep(iG)%V,SRep2%SmolyakRep(iG)%V)
    END DO
  END IF

END IF

IF(openmpi) THEN
  IF(MPI_scheme==1) CALL MPI_Reduce_sum_Bcast(R)
  IF(MPI_scheme==3) CALL MPI_Reduce_sum_Bcast(R,3)
ENDIF

END FUNCTION dot_product_SmolyakRep_Basis

SUBROUTINE SmolyakRep2_TO_SmolyakRep1(SRep1,SRep2)
USE mod_system
IMPLICIT NONE

TYPE(Type_SmolyakRep),           intent(in)     :: SRep2
CLASS(Type_SmolyakRep),          intent(inout)  :: SRep1

integer               :: iG

CALL dealloc_SmolyakRep(SRep1)

IF (allocated(SRep2%SmolyakRep)) THEN
IF (size(SRep2%SmolyakRep) > 0) THEN
  SRep1%Grid  = SRep2%Grid
  SRep1%Delta = SRep2%Delta
  SRep1%nb0   = SRep2%nb0

  allocate(SRep1%SmolyakRep(size(SRep2%SmolyakRep)) )

  DO iG=1,size(SRep2%SmolyakRep)
    SRep1%SmolyakRep(iG) = SRep2%SmolyakRep(iG)
  END DO

END IF
END IF

END SUBROUTINE SmolyakRep2_TO_SmolyakRep1
SUBROUTINE SmolyakRepC2_TO_SmolyakRepC1(SRep1,SRep2)
USE mod_system
IMPLICIT NONE

TYPE(Type_SmolyakRepC),           intent(in)     :: SRep2
TYPE(Type_SmolyakRepC),           intent(inout)  :: SRep1

integer               :: iG

CALL dealloc_SmolyakRepC(SRep1)

IF (size(SRep2%SmolyakRep) > 0) THEN
  SRep1%Grid  = SRep2%Grid
  SRep1%Delta = SRep2%Delta
  SRep1%nb0   = SRep2%nb0

  allocate(SRep1%SmolyakRep(size(SRep2%SmolyakRep)) )

  DO iG=1,size(SRep2%SmolyakRep)
    SRep1%SmolyakRep(iG) = SRep2%SmolyakRep(iG)
  END DO

END IF

END SUBROUTINE SmolyakRepC2_TO_SmolyakRepC1
FUNCTION SmolyakRep1_TIME_SmolyakRep2(SRep1,SRep2) RESULT(SRep)
USE mod_system
IMPLICIT NONE

TYPE(Type_SmolyakRep),           intent(in)     :: SRep1,SRep2
TYPE(Type_SmolyakRep)                           :: SRep


integer               :: iG,nb_BG

nb_BG = Size_SmolyakRep(SRep1)
IF (nb_BG /= Size_SmolyakRep(SRep2)) THEN
  write(out_unitp,*) 'ERROR in SmolyakRep1_TIME_SmolyakRep2'
  write(out_unitp,*) 'sizes are different',Size_SmolyakRep(SRep1),Size_SmolyakRep(SRep2)
  STOP
END IF

SRep = SRep1 !  for the allocation

IF (nb_BG > 0) THEN

  DO iG=lbound(SRep1%SmolyakRep,dim=1),ubound(SRep1%SmolyakRep,dim=1)
    SRep%SmolyakRep(iG)%V(:) = SRep1%SmolyakRep(iG)%V(:) * SRep2%SmolyakRep(iG)%V(:)
  END DO

END IF

END FUNCTION SmolyakRep1_TIME_SmolyakRep2
FUNCTION R1_TIME_SmolyakRep2(R1,SRep2) RESULT(SRep)
USE mod_system
IMPLICIT NONE

TYPE(Type_SmolyakRep),           intent(in)     :: SRep2
real (kind=Rkind),               intent(in)     :: R1
TYPE(Type_SmolyakRep)                           :: SRep


integer               :: iG,nb_BG

nb_BG = Size_SmolyakRep(SRep2)

SRep = SRep2 !  for the allocation

IF (nb_BG > 0) THEN

  DO iG=lbound(SRep2%SmolyakRep,dim=1),ubound(SRep2%SmolyakRep,dim=1)
    SRep%SmolyakRep(iG)%V(:) = SRep2%SmolyakRep(iG)%V(:) * R1
  END DO

END IF

END FUNCTION R1_TIME_SmolyakRep2
FUNCTION SmolyakRep1_TIME_R2(SRep1,R2) RESULT(SRep)
USE mod_system
IMPLICIT NONE

TYPE(Type_SmolyakRep),           intent(in)     :: SRep1
real (kind=Rkind),               intent(in)     :: R2
TYPE(Type_SmolyakRep)                           :: SRep


integer               :: iG,nb_BG

nb_BG = Size_SmolyakRep(SRep1)

SRep = SRep1 !  for the allocation

IF (nb_BG > 0) THEN

  DO iG=lbound(SRep1%SmolyakRep,dim=1),ubound(SRep1%SmolyakRep,dim=1)
    SRep%SmolyakRep(iG)%V(:) = SRep1%SmolyakRep(iG)%V(:) * R2
  END DO

END IF

END FUNCTION SmolyakRep1_TIME_R2
FUNCTION SmolyakRepC1_TIME_SmolyakRepC2(SRep1,SRep2) RESULT(SRep)
USE mod_system
IMPLICIT NONE

TYPE(Type_SmolyakRepC),           intent(in)     :: SRep1,SRep2
TYPE(Type_SmolyakRepC)                           :: SRep


integer               :: iG,nb_BG

nb_BG = Size_SmolyakRepC(SRep1)
IF (nb_BG /= Size_SmolyakRepC(SRep2)) THEN
  write(out_unitp,*) 'ERROR in SmolyakRep1_TIME_SmolyakRep2'
  write(out_unitp,*) 'sizes are different',Size_SmolyakRepC(SRep1),Size_SmolyakRepC(SRep2)
  STOP
END IF

SRep = SRep1 !  for the allocation

IF (nb_BG > 0) THEN

  DO iG=lbound(SRep1%SmolyakRep,dim=1),ubound(SRep1%SmolyakRep,dim=1)
    SRep%SmolyakRep(iG)%V(:) = SRep1%SmolyakRep(iG)%V(:) * SRep2%SmolyakRep(iG)%V(:)
  END DO

END IF

END FUNCTION SmolyakRepC1_TIME_SmolyakRepC2
FUNCTION SmolyakRep1_PLUS_SmolyakRep2(SRep1,SRep2) RESULT(SRep)
USE mod_system
IMPLICIT NONE

TYPE(Type_SmolyakRep),           intent(in)     :: SRep1,SRep2
TYPE(Type_SmolyakRep)                           :: SRep


integer               :: iG,nb_BG

nb_BG = Size_SmolyakRep(SRep1)
IF (nb_BG /= Size_SmolyakRep(SRep2)) THEN
  write(out_unitp,*) 'ERROR in SmolyakRep1_PLUS_SmolyakRep2'
  write(out_unitp,*) 'sizes are different',Size_SmolyakRep(SRep1),Size_SmolyakRep(SRep2)
  STOP
END IF

SRep = SRep1 !  for the allocation

IF (nb_BG > 0) THEN

  DO iG=lbound(SRep1%SmolyakRep,dim=1),ubound(SRep1%SmolyakRep,dim=1)
    SRep%SmolyakRep(iG)%V(:) = SRep1%SmolyakRep(iG)%V(:) + SRep2%SmolyakRep(iG)%V(:)
  END DO

END IF

END FUNCTION SmolyakRep1_PLUS_SmolyakRep2
FUNCTION SmolyakRepC1_PLUS_SmolyakRepC2(SRep1,SRep2) RESULT(SRep)
USE mod_system
IMPLICIT NONE

TYPE(Type_SmolyakRepC),           intent(in)     :: SRep1,SRep2
TYPE(Type_SmolyakRepC)                           :: SRep


integer               :: iG,nb_BG

nb_BG = Size_SmolyakRepC(SRep1)
IF (nb_BG /= Size_SmolyakRepC(SRep2)) THEN
  write(out_unitp,*) 'ERROR in SmolyakRep1_PLUS_SmolyakRep2'
  write(out_unitp,*) 'sizes are different',Size_SmolyakRepC(SRep1),Size_SmolyakRepC(SRep2)
  STOP
END IF

SRep = SRep1 !  for the allocation

IF (nb_BG > 0) THEN

  DO iG=lbound(SRep1%SmolyakRep,dim=1),ubound(SRep1%SmolyakRep,dim=1)
    SRep%SmolyakRep(iG)%V(:) = SRep1%SmolyakRep(iG)%V(:) + SRep2%SmolyakRep(iG)%V(:)
  END DO

END IF

END FUNCTION SmolyakRepC1_PLUS_SmolyakRepC2
FUNCTION SmolyakRep1_MINUS_SmolyakRep2(SRep1,SRep2) RESULT(SRep)
USE mod_system
IMPLICIT NONE

TYPE(Type_SmolyakRep),           intent(in)     :: SRep1,SRep2
TYPE(Type_SmolyakRep)                           :: SRep


integer               :: iG,nb_BG

nb_BG = Size_SmolyakRep(SRep1)
IF (nb_BG /= Size_SmolyakRep(SRep2)) THEN
  write(out_unitp,*) 'ERROR in SmolyakRep1_MINUS_SmolyakRep2'
  write(out_unitp,*) 'sizes are different',Size_SmolyakRep(SRep1),Size_SmolyakRep(SRep2)
  STOP
END IF

SRep = SRep1 !  for the allocation

IF (nb_BG > 0) THEN

  DO iG=lbound(SRep1%SmolyakRep,dim=1),ubound(SRep1%SmolyakRep,dim=1)
    SRep%SmolyakRep(iG)%V(:) = SRep1%SmolyakRep(iG)%V(:) - SRep2%SmolyakRep(iG)%V(:)
  END DO

END IF

END FUNCTION SmolyakRep1_MINUS_SmolyakRep2
FUNCTION SmolyakRepC1_MINUS_SmolyakRepC2(SRep1,SRep2) RESULT(SRep)
USE mod_system
IMPLICIT NONE

TYPE(Type_SmolyakRepC),           intent(in)     :: SRep1,SRep2
TYPE(Type_SmolyakRepC)                           :: SRep


integer               :: iG,nb_BG

nb_BG = Size_SmolyakRepC(SRep1)
IF (nb_BG /= Size_SmolyakRepC(SRep2)) THEN
  write(out_unitp,*) 'ERROR in SmolyakRep1_MINUS_SmolyakRep2'
  write(out_unitp,*) 'sizes are different',Size_SmolyakRepC(SRep1),Size_SmolyakRepC(SRep2)
  STOP
END IF

SRep = SRep1 !  for the allocation

IF (nb_BG > 0) THEN

  DO iG=lbound(SRep1%SmolyakRep,dim=1),ubound(SRep1%SmolyakRep,dim=1)
    SRep%SmolyakRep(iG)%V(:) = SRep1%SmolyakRep(iG)%V(:) - SRep2%SmolyakRep(iG)%V(:)
  END DO

END IF

END FUNCTION SmolyakRepC1_MINUS_SmolyakRepC2
SUBROUTINE GSmolyakRep_TO_BSmolyakRep(SRep,tab_ind,tab_ba,nb0)
USE mod_system
USE mod_basis_set_alloc
IMPLICIT NONE

TYPE(Type_SmolyakRep),           intent(inout)          :: SRep
integer,                         intent(in)             :: tab_ind(:,:) ! tab_ind(D,MaxnD)
TYPE(basis),                     intent(in)             :: tab_ba(0:,:) ! tab_ba(L,D)
integer,                         intent(in)             :: nb0

integer               :: D,iG
integer, allocatable  :: tab_nb(:),tab_nq(:)


IF (.NOT. SRep%Grid) STOP 'Basis is not possible in GSmolyakRep_TO_BSmolyakRep'

!!$OMP   SHARED(nb_mult_GTOB,D,SRep,tab_ind,tab_ba) &
!!$OMP   PRIVATE(iG,tab_nb,tab_nq,i,ib,iq,nnb,nnq,nb2,nq2,RTempG,RTempB) &

D = size(tab_ba(0,:))
nb_mult_GTOB = 0
!$OMP   PARALLEL DEFAULT(NONE) &
!$OMP   SHARED(SRep,tab_ind,tab_ba,nb0,D) &
!$OMP   PRIVATE(iG,tab_nb,tab_nq) &
!$OMP   NUM_THREADS(SG4_maxth)

 allocate(tab_nb(D))
 allocate(tab_nq(D))

!$OMP   DO SCHEDULE(STATIC)

DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)

  tab_nq = getbis_tab_nq(tab_ind(:,iG),tab_ba)
  tab_nb = getbis_tab_nb(tab_ind(:,iG),tab_ba)

  CALL GDP_TO_BDP_OF_SmolyakRep(SRep%SmolyakRep(iG)%V,tab_ba,           &
                                tab_ind(:,iG),tab_nq,tab_nb,nb0)

END DO
!$OMP   END DO

IF (allocated(tab_nb)) deallocate(tab_nb)
IF (allocated(tab_nq)) deallocate(tab_nq)


!$OMP   END PARALLEL

SRep%Grid = .FALSE.

END SUBROUTINE GSmolyakRep_TO_BSmolyakRep

SUBROUTINE GSmolyakRep_TO3_BSmolyakRep(SRep,SGType2,tab_ba)
USE mod_system
USE mod_basis_set_alloc
IMPLICIT NONE

TYPE(Type_SmolyakRep),           intent(inout)          :: SRep
TYPE(basis),                     intent(in)             :: tab_ba(0:,:) ! tab_ba(L,D)
TYPE (param_SGType2),            intent(in)             :: SGType2

integer               :: D,iG,ith,nb_thread,err_sub
integer, allocatable  :: tab_nb(:),tab_nq(:),tab_l(:)


IF (.NOT. SRep%Grid) STOP 'Basis is not possible in GSmolyakRep_TO3_BSmolyakRep'

D = SGType2%nDind_SmolyakRep%ndim
nb_mult_GTOB = 0

!to be sure to have the correct number of threads
nb_thread = SGType2%nb_threads

!!$OMP   SHARED(nb_mult_GTOB,D,SRep,SGType2,tab_ba,nb_thread) &
!!$OMP   PRIVATE(iG,tab_nb,tab_nq,i,ib,iq,nnb,nnq,nb2,nq2,RTempG,RTempB) &
!!$OMP   PRIVATE(ith,tab_l,err_sub) &

!$OMP   PARALLEL DEFAULT(NONE) &
!$OMP   SHARED(D,SRep,SGType2,tab_ba,nb_thread) &
!$OMP   PRIVATE(iG,tab_nb,tab_nq,ith,tab_l,err_sub) &
!$OMP   NUM_THREADS(nb_thread)

allocate(tab_l(D))
allocate(tab_nb(D))
allocate(tab_nq(D))

!--------------------------------------------------------------
!-- For the initialization of tab_l(:) and the use of ADD_ONE_TO_nDindex in the parallel loop
ith = 0
!$ ith = OMP_GET_THREAD_NUM()
tab_l(:) = SGType2%nDval_init(:,ith+1)
!--------------------------------------------------------------

! we are not using the parallel do, to be able to use the correct initialized tab_l with nDval_init
DO iG=SGType2%iG_th(ith+1),SGType2%fG_th(ith+1)
  CALL ADD_ONE_TO_nDindex(SGType2%nDind_SmolyakRep,tab_l,iG=iG,err_sub=err_sub)

  tab_nq = getbis_tab_nq(tab_l,tab_ba)
  tab_nb = getbis_tab_nb(tab_l,tab_ba)

  CALL GDP_TO_BDP_OF_SmolyakRep(SRep%SmolyakRep(iG)%V,tab_ba,           &
                                tab_l,tab_nq,tab_nb,nb0=SGType2%nb0)

END DO

IF (allocated(tab_l)) deallocate(tab_l)
IF (allocated(tab_nb)) deallocate(tab_nb)
IF (allocated(tab_nq)) deallocate(tab_nq)

!$OMP   END PARALLEL

SRep%Grid = .FALSE.

END SUBROUTINE GSmolyakRep_TO3_BSmolyakRep

SUBROUTINE BSmolyakRep_TO_GSmolyakRep(SRep,tab_ind,tab_ba,nb0)
USE mod_system
USE mod_basis_set_alloc
IMPLICIT NONE

real(kind=Rkind)  :: R
TYPE(Type_SmolyakRep),           intent(inout)          :: SRep
integer,                         intent(in)             :: tab_ind(:,:) ! tab_ind(D,MaxnD)
TYPE(basis),                     intent(in)             :: tab_ba(0:,:) ! tab_ba(0:L,D)
integer,                         intent(in)             :: nb0


integer                            :: iG
integer, allocatable               :: tab_nb(:),tab_nq(:)


integer                            :: i,D,nb_BG
!integer                            :: nnb,nnq,nb2,nq2,ib,iq
!real(kind=Rkind), allocatable      :: RTempG(:,:,:),RTempB(:,:,:)
integer                            :: d1,d2

IF (SRep%Grid) STOP 'Grid is not possible in BSmolyakRep_TO_GSmolyakRep'

!write(6,*) 'BSmolyakRep_TO_GSmolyakRep: shape(tab_ind) ',shape(tab_ind)
!write(6,*) 'BSmolyakRep_TO_GSmolyakRep: shape(tab_ba) ',shape(tab_ba)
!write(6,*) 'BSmolyakRep_TO_GSmolyakRep: nb0 ',nb0
!flush(6)

D = size(tab_ba(0,:))
!nb_mult_BTOG = 0

d1=lbound(SRep%SmolyakRep,dim=1)
d2=ubound(SRep%SmolyakRep,dim=1)
! IF(openmpi .AND. MPI_scheme==1) THEN
!   d1=iGs_MPI(1,MPI_id)
!   d2=iGs_MPI(2,MPI_id)
! ENDIF
IF(openmpi) THEN
  IF(MPI_scheme==1) THEN
    d1=iGs_MPI(1,MPI_id)
    d2=iGs_MPI(2,MPI_id)
  ELSEIF(MPI_scheme==3) THEN
    IF(MPI_nodes_p0) THEN
      d1=iGs_MPI(1,MPI_id)
      d2=iGs_MPI(2,MPI_sub_id(2))
    ELSE
      d1=0
      d2=-1
    ENDIF
  ENDIF
ENDIF

!!!$OMP   SHARED(nb_mult_BTOG,D,SRep,tab_ind,tab_ba) &
!!!$OMP   PRIVATE(iG,tab_nb,tab_nq,i,ib,iq,nnb,nnq,nb2,nq2,RTempG,RTempB) &

!$OMP   PARALLEL DEFAULT(NONE) &
!$OMP   SHARED(SRep,tab_ind,tab_ba,nb0,D) &
!$OMP   SHARED(openmpi,MPI_scheme,iGs_MPI,MPI_id,d1,d2) &
!$OMP   PRIVATE(iG,tab_nb,tab_nq) &
!$OMP   NUM_THREADS(SG4_maxth)

allocate(tab_nb(D))
allocate(tab_nq(D))

!$OMP   DO SCHEDULE(STATIC)

!DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)
DO iG=d1,d2

  !write(6,*) iG,'tab_ind(:,iG)',tab_ind(:,iG)
  tab_nq = getbis_tab_nq(tab_ind(:,iG),tab_ba)
  tab_nb = getbis_tab_nb(tab_ind(:,iG),tab_ba)

  CALL BDP_TO_GDP_OF_SmolyakRep(SRep%SmolyakRep(iG)%V,tab_ba,           &
                                tab_ind(:,iG),tab_nq,tab_nb,nb0=nb0)

END DO
!$OMP   END DO

IF (allocated(tab_nb)) deallocate(tab_nb)
IF (allocated(tab_nq)) deallocate(tab_nq)


!$OMP   END PARALLEL

SRep%Grid = .TRUE.
!write(6,*) 'END BSmolyakRep_TO_GSmolyakRep'
!flush(6)

END SUBROUTINE BSmolyakRep_TO_GSmolyakRep

SUBROUTINE BSmolyakRep_TO3_GSmolyakRep(SRep,SGType2,tab_ba)
USE mod_system
USE mod_basis_set_alloc
IMPLICIT NONE

TYPE(Type_SmolyakRep),           intent(inout)          :: SRep
TYPE(basis),                     intent(in)             :: tab_ba(0:,:) ! tab_ba(0:L,D)
TYPE (param_SGType2),            intent(in)             :: SGType2


integer                            :: iG,D
integer, allocatable               :: tab_nb(:),tab_nq(:),tab_l(:)


integer               :: i,nb_BG,ith,nb_thread,err_sub

!integer                            :: nnb,nnq,nb2,nq2,ib,iq
!real(kind=Rkind), allocatable      :: RTempG(:,:,:),RTempB(:,:,:)
integer                            :: d1,d2

IF (SRep%Grid) STOP 'Grid is not possible in BSmolyakRep_TO3_GSmolyakRep'


D = SGType2%nDind_SmolyakRep%ndim
nb_mult_BTOG = 0

!to be sure to have the correct number of threads
nb_thread = SGType2%nb_threads

!$OMP   PARALLEL DEFAULT(NONE) &
!$OMP   SHARED(D,SRep,tab_ba,SGType2,nb_thread) &
!$OMP   SHARED(openmpi,MPI_scheme,iGs_MPI,MPI_id) &
!$OMP   PRIVATE(iG,tab_nb,tab_nq,tab_l,ith,err_sub,d1,d2) &
!$OMP   NUM_THREADS(nb_thread)

allocate(tab_l(D))
allocate(tab_nb(D))
allocate(tab_nq(D))

!--------------------------------------------------------------
!-- For the initialization of tab_l(:) and the use of ADD_ONE_TO_nDindex in the parallel loop
ith = 0
!$ ith = OMP_GET_THREAD_NUM()
tab_l(:) = SGType2%nDval_init(:,ith+1)
!--------------------------------------------------------------

d1=SGType2%iG_th(ith+1)
d2=SGType2%fG_th(ith+1)
IF(openmpi .AND. MPI_scheme==1) THEN
  d1=iGs_MPI(1,MPI_id)
  d2=iGs_MPI(2,MPI_id)
ENDIF

! we are not using the parallel do, to be able to use the correct initialized tab_l with nDval_init
!DO iG=SGType2%iG_th(ith+1),SGType2%fG_th(ith+1)
DO iG=d1,d2

  CALL ADD_ONE_TO_nDindex(SGType2%nDind_SmolyakRep,tab_l,iG=iG,err_sub=err_sub)

  tab_nq = getbis_tab_nq(tab_l,tab_ba)
  tab_nb = getbis_tab_nb(tab_l,tab_ba)

  CALL BDP_TO_GDP_OF_SmolyakRep(SRep%SmolyakRep(iG)%V,tab_ba,           &
                                tab_l,tab_nq,tab_nb,nb0=SGType2%nb0)

END DO

IF (allocated(tab_l)) deallocate(tab_l)
IF (allocated(tab_nb)) deallocate(tab_nb)
IF (allocated(tab_nq)) deallocate(tab_nq)


!$OMP   END PARALLEL

SRep%Grid = .TRUE.

END SUBROUTINE BSmolyakRep_TO3_GSmolyakRep

SUBROUTINE BDP_TO_GDP_OF_SmolyakRep(R,tab_ba,tab_l,tab_nq,tab_nb,nb0)
USE mod_system
USE mod_basis_set_alloc
IMPLICIT NONE

real(kind=Rkind), allocatable,   intent(inout)          :: R(:) ! from SRep%SmolyakRep(iG)%V

TYPE(basis),                     intent(in)             :: tab_ba(0:,:) ! tab_ba(0:L,D)
integer,                         intent(in)             :: tab_nb(:),tab_nq(:),tab_l(:)
integer,                         intent(in), optional   :: nb0

integer                            :: i,ib0,nb0_loc
integer                            :: nnb,nnq,nb2,nq2,ib,iq
real(kind=Rkind), allocatable      :: RTempG(:,:,:),RTempB(:,:,:)
real(kind=Rkind), allocatable      :: RG(:,:),RB(:,:)


  !----- for debuging ----------------------------------------------
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub='BDP_TO_GDP_OF_SmolyakRep'
  logical, parameter :: debug = .FALSE.
  !logical, parameter :: debug = .TRUE.
  !-----------------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*) 'tab_l(:) ',tab_l(:)
    write(out_unitp,*) 'tab_nb(:)',tab_nb(:)
    write(out_unitp,*) 'tab_nq(:)',tab_nq(:)

    IF (present(nb0)) write(out_unitp,*) 'nb0',nb0
    write(out_unitp,*) 'size(R)',size(R)
    CALL flush_perso(out_unitp)
  END IF
  !-----------------------------------------------------------------
 nb0_loc = 1
 IF (present(nb0)) nb0_loc = nb0

 nnb = product(tab_nb)
 CALL alloc_NParray(RB,[nnb,nb0_loc],'RB',name_sub)
 RB(:,:)  = reshape(R,shape=[nnb,nb0_loc])

 nnq = product(tab_nq)
 CALL alloc_NParray(RG,[nnq,nb0_loc],'RG',name_sub)
 RG(:,:)  = ZERO

 CALL dealloc_NParray(R,'R',name_sub)

 DO ib0=1,nb0_loc

  nnb = product(tab_nb)
  nnq = 1
  nb2 = 1
  nq2 = 1

  CALL alloc_NParray(RTempG,[nnq,nq2,nnb],'RTempG',name_sub)
  RTempG(:,:,:) = reshape(RB(:,ib0),shape=[nnq,nq2,nnb])


  ! B order : b1 * b2 * b3 * ... bD
  ! G order : g1 * g2 * g3 * ... gD

  DO i=1,size(tab_nb) ! D=size(tab_nb)
    nb2 = tab_nb(i)
    nq2 = tab_nq(i)
    nnb = nnb / nb2

    CALL alloc_NParray(RTempB,[nnq,nb2,nnb],'RTempB',name_sub)
    RTempB(:,:,:) = reshape(RTempG,shape=[nnq,nb2,nnb])

    CALL dealloc_NParray(RTempG,'RTempG',name_sub)
    CALL alloc_NParray(RTempG,[nnq,nq2,nnb],'RTempG',name_sub)

    DO ib=1,nnb
    DO iq=1,nnq
       RTempG(iq,:,ib) = matmul(tab_ba(tab_l(i),i)%dnRGB%d0,RTempB(iq,:,ib))
       !dgemv (TRANS, M, N, ALPHA, A, LDA, X, INCX, BETA, Y, INCY)
       ! Y = ALPHA * A.X + BETA * Y (if TRANS='n')
       ! Y = ALPHA * A^t.X + BETA * Y (if TRANS='t')

       !CALL dgemv('n',nq2,nb2,ONE,tab_ba(tab_l(i),i)%dnRGB%d0, &
       !           nq2,RTempB(iq,:,ib),1,ZERO,RTempG(iq,:,ib),1)

      !$OMP ATOMIC
      nb_mult_BTOG = nb_mult_BTOG + int(nq2,kind=ILkind)*int(nb2,kind=ILkind)

    END DO
    END DO

    nnq = nnq * nq2
    CALL dealloc_NParray(RTempB,'RTempB',name_sub)

  END DO

  RG(:,ib0) = reshape(RTempG, shape=[nnq] )
  CALL dealloc_NParray(RTempG,'RTempG',name_sub)

 END DO

  nnq = product(tab_nq)
  CALL alloc_NParray(R,[nnq*nb0_loc],'R',name_sub)
  R(:) = reshape(RG, shape=[nnq*nb0_loc] )

  CALL dealloc_NParray(RB,'RB',name_sub)
  CALL dealloc_NParray(RG,'RG',name_sub)

  IF (debug) THEN
    write(out_unitp,*) 'size(R)',size(R)
    write(out_unitp,*) 'END ',name_sub
    CALL flush_perso(out_unitp)
  END IF

END SUBROUTINE BDP_TO_GDP_OF_SmolyakRep
SUBROUTINE GDP_TO_BDP_OF_SmolyakRep(R,tab_ba,tab_l,tab_nq,tab_nb,nb0)
USE mod_system
USE mod_basis_set_alloc
IMPLICIT NONE

real(kind=Rkind), allocatable,   intent(inout)          :: R(:) ! from SRep%SmolyakRep(iG)%V

TYPE(basis),                     intent(in)             :: tab_ba(0:,:) ! tab_ba(0:L,D)
integer,                         intent(in)             :: tab_nb(:),tab_nq(:),tab_l(:)
integer,                         intent(in), optional   :: nb0

integer                            :: i,ib0,nb0_loc
integer                            :: nnb,nnq,nb2,nq2,ib,iq
real(kind=Rkind), allocatable      :: RTempG(:,:,:),RTempB(:,:,:)
real(kind=Rkind), allocatable      :: RG(:,:),RB(:,:)

 character (len=*), parameter :: name_sub='GDP_TO_BDP_OF_SmolyakRep'

 nb0_loc = 1
 IF (present(nb0)) nb0_loc = nb0

 nnq = product(tab_nq)
 CALL alloc_NParray(RG,[nnq,nb0_loc],'RG',name_sub)
 RG(:,:)  = reshape(R,shape=[nnq,nb0_loc])
 CALL dealloc_NParray(R,'R',name_sub)
 nnb = product(tab_nb)
 CALL alloc_NParray(RB,[nnb,nb0_loc],'RB',name_sub)
 RB(:,:) = ZERO

 DO ib0=1,nb0_loc
  nnq = product(tab_nq)
  nnb = 1
  nb2 = 1
  nq2 = 1

  CALL alloc_NParray(RTempB,[1,1,nnq],'RTempB',name_sub)
  RTempB(:,:,:) = reshape(RG(:,ib0),shape=[1,1,nnq])

! B order : b1 * b2 * b3 * ... bD
! G order : g1 * g2 * g3 * ... gD

  DO i=1,size(tab_nb) ! D=size(tab_nb)
    nb2 = tab_nb(i)
    nq2 = tab_nq(i)
    nnq = nnq / nq2

    CALL alloc_NParray(RTempG,[nnb,nq2,nnq],'RTempG',name_sub)
    RTempG(:,:,:) = reshape(RTempB,shape=[nnb,nq2,nnq])

    CALL dealloc_NParray(RTempB,'RTempB',name_sub)
    CALL alloc_NParray(RTempB,[nnb,nb2,nnq],'RTempB',name_sub)

    DO iq=1,nnq
    DO ib=1,nnb
!#if __LAPACK == 1
!       CALL dgemv('n',nb2,nq2,ONE,tab_ba(tab_l(i),i)%dnRBGwrho%d0, &
!                  nb2,RTempG(ib,:,iq),1,ZERO,RTempB(ib,:,iq),1)
!#else
       RTempB(ib,:,iq) = matmul(tab_ba(tab_l(i),i)%dnRBGwrho%d0 , RTempG(ib,:,iq))
!#endif
      !$OMP ATOMIC
      nb_mult_GTOB = nb_mult_GTOB + int(nq2,kind=ILkind)*int(nb2,kind=ILkind)

    END DO
    END DO

    nnb = nnb * nb2
    CALL dealloc_NParray(RTempG,'RTempG',name_sub)

  END DO

  RB(:,ib0) = reshape(RTempB, shape=[nnb] )

  CALL dealloc_NParray(RTempB,'RTempB',name_sub)
 END DO


  CALL alloc_NParray(R,[nnb*nb0_loc],'R',name_sub)
  R(:) = reshape(RB, shape=[nnb*nb0_loc] )

  CALL dealloc_NParray(RB,'RB',name_sub)
  CALL dealloc_NParray(RG,'RG',name_sub)


 END SUBROUTINE GDP_TO_BDP_OF_SmolyakRep

SUBROUTINE DerivOp_TO3_GSmolyakRep(SRep,SGType2,tab_ba,tab_der)
USE mod_system
USE mod_basis_set_alloc
IMPLICIT NONE

TYPE(Type_SmolyakRep),           intent(inout)          :: SRep
TYPE(basis),                     intent(in)             :: tab_ba(0:,:) ! tab_ba(L,D)
TYPE (param_SGType2),            intent(in)             :: SGType2
integer,                         intent(in)             :: tab_der(2)

integer               :: D,iG,ith,nb_thread,err_sub
integer, allocatable  :: tab_nq(:),tab_l(:)


IF (.NOT. SRep%Grid) STOP 'Basis is not possible in DerivOp_TO3_GSmolyakRep'

D = SGType2%nDind_SmolyakRep%ndim

!to be sure to have the correct number of threads
nb_thread = SGType2%nb_threads

!$OMP   PARALLEL DEFAULT(NONE) &
!$OMP   SHARED(D,SRep,SGType2,tab_ba,tab_der) &
!$OMP   PRIVATE(iG,ith,tab_l,tab_nq,err_sub) &
!$OMP   NUM_THREADS(SGType2%nb_threads)

allocate(tab_l(D))
allocate(tab_nq(D))

!--------------------------------------------------------------
!-- For the initialization of tab_l(:) and the use of ADD_ONE_TO_nDindex in the parallel loop
ith = 0
!$ ith = OMP_GET_THREAD_NUM()
tab_l(:) = SGType2%nDval_init(:,ith+1)
!--------------------------------------------------------------

! we are not using the parallel do, to be able to use the correct initialized tab_l with nDval_init
DO iG=SGType2%iG_th(ith+1),SGType2%fG_th(ith+1)
  CALL ADD_ONE_TO_nDindex(SGType2%nDind_SmolyakRep,tab_l,iG=iG,err_sub=err_sub)

  tab_nq = getbis_tab_nq(tab_l,tab_ba)

  CALL DerivOp_TO_RDP_OF_SmolaykRep(SRep%SmolyakRep(iG)%V,tab_ba,tab_l,tab_nq,tab_der)

END DO

IF (allocated(tab_l)) deallocate(tab_l)
IF (allocated(tab_nq)) deallocate(tab_nq)

!$OMP   END PARALLEL

END SUBROUTINE DerivOp_TO3_GSmolyakRep

SUBROUTINE DerivOp_TO3_GSmolyakRepC(SRep,SGType2,tab_ba,tab_der)
USE mod_system
USE mod_basis_set_alloc
IMPLICIT NONE

TYPE(Type_SmolyakRepC),          intent(inout)          :: SRep
TYPE(basis),                     intent(in)             :: tab_ba(0:,:) ! tab_ba(L,D)
TYPE (param_SGType2),            intent(in)             :: SGType2
integer,                         intent(in)             :: tab_der(2)

integer               :: D,iG,ith,nb_thread,err_sub
integer, allocatable  :: tab_nq(:),tab_l(:)


IF (.NOT. SRep%Grid) STOP 'Basis is not possible in DerivOp_TO3_GSmolyakRepC'

D = SGType2%nDind_SmolyakRep%ndim

!to be sure to have the correct number of threads
nb_thread = SGType2%nb_threads

!$OMP   PARALLEL DEFAULT(NONE) &
!$OMP   SHARED(D,SRep,SGType2,tab_ba,tab_der) &
!$OMP   PRIVATE(iG,ith,tab_l,tab_nq,err_sub) &
!$OMP   NUM_THREADS(SGType2%nb_threads)

allocate(tab_l(D))
allocate(tab_nq(D))

!--------------------------------------------------------------
!-- For the initialization of tab_l(:) and the use of ADD_ONE_TO_nDindex in the parallel loop
ith = 0
!$ ith = OMP_GET_THREAD_NUM()
tab_l(:) = SGType2%nDval_init(:,ith+1)
!--------------------------------------------------------------

! we are not using the parallel do, to be able to use the correct initialized tab_l with nDval_init
DO iG=SGType2%iG_th(ith+1),SGType2%fG_th(ith+1)
  CALL ADD_ONE_TO_nDindex(SGType2%nDind_SmolyakRep,tab_l,iG=iG,err_sub=err_sub)

  tab_nq = getbis_tab_nq(tab_l,tab_ba)

  CALL DerivOp_TO_RDP_OF_SmolaykRepC(SRep%SmolyakRep(iG)%V,tab_ba,tab_l,tab_nq,tab_der)

END DO

IF (allocated(tab_l)) deallocate(tab_l)
IF (allocated(tab_nq)) deallocate(tab_nq)

!$OMP   END PARALLEL

END SUBROUTINE DerivOp_TO3_GSmolyakRepC


SUBROUTINE DerivOp_TO_RDP_OF_SmolaykRep(R,tab_ba,tab_l,tab_nq,tab_der)
  USE mod_basis_set_alloc
  IMPLICIT NONE
  TYPE(basis),                     intent(in)             :: tab_ba(0:,:) ! tab_ba(0:L,D)
  integer,                         intent(in)             :: tab_nq(:),tab_l(:)
  real (kind=Rkind),               intent(inout)          :: R(:)
  integer, optional,               intent(in)             :: tab_der(2)

  integer :: tab_der_loc(2),dnba_ind(2),iQact_ba,k
  real (kind=Rkind), allocatable       :: RG1(:,:,:)
  real (kind=Rkind), allocatable       :: RG2(:,:,:)
  real (kind=Rkind), allocatable       :: BGG(:,:)

  integer                          :: ibasis,nq

  integer                          :: nnq1,nnq3,nq2,iq1,iq2,iq3

  !----- for debuging --------------------------------------------------
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub='DerivOp_TO_RDP_OF_SmolaykRep'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
  !---------------------------------------------------------------------

  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*) 'R(:)',R(:)
    CALL flush_perso(out_unitp)
  END IF

  IF (present(tab_der)) THEN
    tab_der_loc(:) = tab_der(:)
  ELSE
    tab_der_loc(:) = 0
  END IF
  IF (debug) write(out_unitp,*) 'tab_der_loc',tab_der_loc
  WHERE (tab_der_loc < 0) tab_der_loc = 0

  nq = size(R)

  IF (nq < 1 .OR. any(tab_nq < 1) .OR. any(tab_l < 0)) THEN
    write(out_unitp,*) ' Wrong arguments values'
    write(out_unitp,*) ' nq size(R): ',nq
    write(out_unitp,*) ' tab_l(:)  : ',tab_l
    write(out_unitp,*) ' tab_nq(:) : ',tab_nq
    STOP 'ERROR in DerivOp_TO_RDP_OF_SmolaykRep'
  END IF

  IF (any(tab_der_loc /= 0)) THEN

     nnq3 = nq
     nnq1 = 1
     nq2  = 1

     CALL alloc_NParray(RG1,[nnq1,nq2,nnq3],"RG1",name_sub)
     RG1(:,:,:) = reshape(R,shape=[nnq1,nq2,nnq3])

     DO ibasis=1,size(tab_l)

       nq2  = tab_nq(ibasis)
       nnq3 = nnq3 / nq2

       dnba_ind(:) = tab_ba(tab_l(ibasis),ibasis)%Tabder_Qdyn_TO_Qbasis(tab_der_loc(:))

       IF (dnba_ind(1) /= 0 .OR. dnba_ind(2) /= 0) THEN
         CALL alloc_NParray(RG2,[nnq1,nq2,nnq3],"RG2",name_sub)
         RG2(:,:,:) = reshape(RG1,shape=[nnq1,nq2,nnq3])

         IF (tab_ba(tab_l(ibasis),ibasis)%packed) THEN

           CALL alloc_NParray(BGG,[nq2,nq2],"BGG",name_sub)

           CALL Get_MatdnRGG(tab_ba(tab_l(ibasis),ibasis),BGG,dnba_ind)

           DO iq3=1,nnq3
           DO iq1=1,nnq1
              RG2(iq1,:,iq3) = matmul(BGG,RG2(iq1,:,iq3))
           END DO
           END DO

           CALL dealloc_NParray(BGG,"BGG",name_sub)

         ELSE
           STOP 'not packed'
         END IF

         CALL dealloc_NParray(RG1,"RG1",name_sub)
         CALL alloc_NParray(RG1,[nnq1,nq2,nnq3],"RG1",name_sub)
         RG1(:,:,:) = RG2
         CALL dealloc_NParray(RG2,"RG2",name_sub)
       END IF

       nnq1 = nnq1 * nq2

     END DO

     R(:) = reshape(RG1, shape=[nq] )
     CALL dealloc_NParray(RG1,"RG1",name_sub)
   END IF

  IF (debug) THEN
    write(out_unitp,*) 'R(:)',R(:)
    write(out_unitp,*) 'END ',name_sub
  END IF

END SUBROUTINE DerivOp_TO_RDP_OF_SmolaykRep

SUBROUTINE DerivOp_TO_RDP_OF_SmolaykRepC(R,tab_ba,tab_l,tab_nq,tab_der)
  USE mod_basis_set_alloc
  IMPLICIT NONE
  TYPE(basis),                     intent(in)             :: tab_ba(0:,:) ! tab_ba(0:L,D)
  integer,                         intent(in)             :: tab_nq(:),tab_l(:)
  complex (kind=Rkind),            intent(inout)          :: R(:)
  integer, optional,               intent(in)             :: tab_der(2)

  integer :: tab_der_loc(2),dnba_ind(2),iQact_ba,k
  complex (kind=Rkind), allocatable       :: RG1(:,:,:)
  complex (kind=Rkind), allocatable       :: RG2(:,:,:)
  real (kind=Rkind), allocatable          :: BGG(:,:)

  integer                          :: ibasis,nq

  integer                          :: nnq1,nnq3,nq2,iq1,iq2,iq3

  !----- for debuging --------------------------------------------------
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub='DerivOp_TO_RDP_OF_SmolaykRep'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
  !---------------------------------------------------------------------

  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*) 'R(:)',R(:)
    CALL flush_perso(out_unitp)
  END IF

  IF (present(tab_der)) THEN
    tab_der_loc(:) = tab_der(:)
  ELSE
    tab_der_loc(:) = 0
  END IF
  IF (debug) write(out_unitp,*) 'tab_der_loc',tab_der_loc
  WHERE (tab_der_loc < 0) tab_der_loc = 0

  nq = size(R)

  IF (nq < 1 .OR. any(tab_nq < 1) .OR. any(tab_l < 0)) THEN
    write(out_unitp,*) ' Wrong arguments values'
    write(out_unitp,*) ' nq size(R): ',nq
    write(out_unitp,*) ' tab_l(:)  : ',tab_l
    write(out_unitp,*) ' tab_nq(:) : ',tab_nq
    STOP 'ERROR in DerivOp_TO_RDP_OF_SmolaykRep'
  END IF

  IF (any(tab_der_loc /= 0)) THEN

     nnq3 = nq
     nnq1 = 1
     nq2  = 1

     CALL alloc_NParray(RG1,[nnq1,nq2,nnq3],"RG1",name_sub)
     RG1(:,:,:) = reshape(R,shape=[nnq1,nq2,nnq3])

     DO ibasis=1,size(tab_l)

       nq2  = tab_nq(ibasis)
       nnq3 = nnq3 / nq2

       dnba_ind(:) = tab_ba(tab_l(ibasis),ibasis)%Tabder_Qdyn_TO_Qbasis(tab_der_loc(:))

       IF (any(dnba_ind /= 0 )) THEN
         CALL alloc_NParray(RG2,[nnq1,nq2,nnq3],"RG2",name_sub)
         RG2(:,:,:) = reshape(RG1,shape=[nnq1,nq2,nnq3])

         IF (tab_ba(tab_l(ibasis),ibasis)%packed) THEN

           CALL alloc_NParray(BGG,[nq2,nq2],"BGG",name_sub)

           CALL Get_MatdnRGG(tab_ba(tab_l(ibasis),ibasis),BGG,dnba_ind)

           DO iq3=1,nnq3
           DO iq1=1,nnq1
              RG2(iq1,:,iq3) = matmul(BGG,RG2(iq1,:,iq3))
           END DO
           END DO

           CALL dealloc_NParray(BGG,"BGG",name_sub)

         ELSE
           STOP 'not packed'
         END IF

         CALL dealloc_NParray(RG1,"RG1",name_sub)
         CALL alloc_NParray(RG1,[nnq1,nq2,nnq3],"RG1",name_sub)
         RG1(:,:,:) = RG2
         CALL dealloc_NParray(RG2,"RG2",name_sub)
       END IF

       nnq1 = nnq1 * nq2

     END DO

     R(:) = reshape(RG1, shape=[nq] )
     CALL dealloc_NParray(RG1,"RG1",name_sub)
   END IF

  IF (debug) THEN
    write(out_unitp,*) 'R(:)',R(:)
    write(out_unitp,*) 'END ',name_sub
  END IF

END SUBROUTINE DerivOp_TO_RDP_OF_SmolaykRepC

SUBROUTINE Set_weight_TO_SmolyakRep(SRep,tab_ind,tab_ba)
USE mod_system
USE mod_basis_set_alloc
IMPLICIT NONE

TYPE(Type_SmolyakRep),           intent(inout)     :: SRep
integer,         allocatable,    intent(in)        :: tab_ind(:,:) ! tab_ind(D,MaxnD)
TYPE(basis),                     intent(in)        :: tab_ba(0:,:) ! tab_ba(0:L,D)

integer               :: i,D,iG,nb_BG

integer                            :: nnq,nnq1,nq2,nnq3,iq1,iq3
integer,          allocatable      :: tab_n(:)
real(kind=Rkind), allocatable      :: RTempG(:,:,:)

character (len=*), parameter :: name_sub='Set_weight_TO_SmolyakRep'


!write(6,*) 'Set_weight_TO_SmolyakRep: shape tab_ba',shape(tab_ba) ; flush(6)

  CALL alloc_SmolyakRep(SRep,tab_ind,tab_ba,grid=.TRUE.)

  !CALL Write_SmolyakRep(Srep)
  CALL alloc_NParray(tab_n,shape(tab_ind(:,1)),'tab_n','alloc_SmolyakRep')


  DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)

    tab_n = getbis_tab_nq(tab_ind(:,iG),tab_ba)

    nnq = product(tab_n)

    nnq1 = nnq
    nq2  = 1
    nnq3 = 1

    SRep%SmolyakRep(iG)%V(:) = ONE

    RTempG = reshape(SRep%SmolyakRep(iG)%V,shape=[nnq,1,1])

    DO i=size(tab_ind(:,iG)),1,-1
      nq2  = tab_n(i)
      nnq1 = nnq1/nq2

      RTempG = reshape(RTempG,shape=[nnq1,nq2,nnq3])

      DO iq3=1,nnq3
      DO iq1=1,nnq1
        RTempG(iq1,:,iq3) = tab_ba(tab_ind(i,iG),i)%wrho(1:nq2) * RTempG(iq1,1:nq2,iq3)
      END DO
      END DO

      nnq3 = nnq3 * nq2


    END DO

    SRep%SmolyakRep(iG)%V(:) = reshape(RTempG, shape=[nnq] )

    deallocate(RTempG)

  END DO

  IF (allocated(tab_n)) CALL dealloc_NParray(tab_n,'tab_n',name_sub)

  !CALL Write_SmolyakRep(Srep)

END SUBROUTINE Set_weight_TO_SmolyakRep

SUBROUTINE Get_weight_FROM_OneDP(wrho,iG,tab_ind,tab_ba)
USE mod_system
USE mod_basis_set_alloc
IMPLICIT NONE

real(kind=Rkind), allocatable,   intent(inout)     :: wrho(:)
integer,                         intent(in)        :: iG
integer,          allocatable,   intent(in)        :: tab_ind(:,:) ! tab_ind(D,MaxnD)
TYPE(basis),                     intent(in)        :: tab_ba(0:,:) ! tab_ba(0:L,D)

integer               :: i,D,nb_BG

integer                            :: nnq,nnq1,nq2,nnq3,iq1,iq3
integer, allocatable               :: tab_n(:)
real(kind=Rkind), allocatable      :: RTempG(:,:,:)

  character (len=*), parameter :: name_sub='Get_weight_FROM_OneDP'


!write(6,*) 'Get_weight_FROM_OneDP: shape tab_ba',shape(tab_ba) ; flush(6)

  CALL alloc_NParray(tab_n,shape(tab_ind(:,1)),'tab_n','alloc_SmolyakRep')

  tab_n = getbis_tab_nq(tab_ind(:,iG),tab_ba)

  nnq = product(tab_n)

  nnq1 = nnq
  nq2  = 1
  nnq3 = 1

  IF (allocated(wrho)) CALL dealloc_NParray(wrho,'wrho',name_sub)
  CALL alloc_NParray(wrho,[nnq],'wrho',name_sub)
  wrho(:) = ONE


  RTempG = reshape(wrho,shape=[nnq,1,1])

  DO i=size(tab_ind(:,iG)),1,-1
    nq2  = tab_n(i)
    nnq1 = nnq1/nq2

    RTempG = reshape(RTempG,shape=[nnq1,nq2,nnq3])

    DO iq3=1,nnq3
    DO iq1=1,nnq1
      RTempG(iq1,:,iq3) = tab_ba(tab_ind(i,iG),i)%wrho(1:nq2) * RTempG(iq1,1:nq2,iq3)
    END DO
    END DO

    nnq3 = nnq3 * nq2

  END DO

  wrho(:) = reshape(RTempG, shape=[nnq] )

  deallocate(RTempG)

  IF (allocated(tab_n)) CALL dealloc_NParray(tab_n,'tab_n',name_sub)

END SUBROUTINE Get_weight_FROM_OneDP

FUNCTION get_tab_nq(tab_l,tab_ba) RESULT (tab_nq)
USE mod_basis_set_alloc
IMPLICIT NONE

  integer,         allocatable                           :: tab_nq(:) ! result
  integer ,                        intent(in)            :: tab_l(:)
  TYPE(basis),                     intent(in)            :: tab_ba(0:,:) ! tab_ba(0:L,D)

  integer :: i

  IF (allocated(tab_nq)) CALL dealloc_NParray(tab_nq,'tab_nq','get_tab_nq')
  CALL alloc_NParray(tab_nq,shape(tab_l),'tab_nq','get_tab_nq')

  DO i=1,size(tab_l)
    tab_nq(i) = get_nq_FROM_basis(tab_ba(tab_l(i),i))
  END DO

END FUNCTION get_tab_nq
FUNCTION getbis_tab_nq(tab_l,tab_ba) RESULT (tab_nq)
USE mod_basis_set_alloc
IMPLICIT NONE

  integer ,                        intent(in)            :: tab_l(:)
  TYPE(basis),                     intent(in)            :: tab_ba(0:,:) ! tab_ba(0:L,D)
  integer                                                :: tab_nq(size(tab_l)) ! result

  integer :: i

  DO i=1,size(tab_l)
    tab_nq(i) = get_nq_FROM_basis(tab_ba(tab_l(i),i))
  END DO

END FUNCTION getbis_tab_nq
FUNCTION get_tab_nb(tab_l,tab_ba) RESULT (tab_nb)
USE mod_basis_set_alloc
IMPLICIT NONE

  integer,         allocatable                           :: tab_nb(:) ! result
  integer ,                        intent(in)            :: tab_l(:)
  TYPE(basis),                     intent(in)            :: tab_ba(0:,:) ! tab_ba(0:L,D)

  integer :: i

  IF (allocated(tab_nb)) CALL dealloc_NParray(tab_nb,'tab_nb','get_tab_nb')
  CALL alloc_NParray(tab_nb,shape(tab_l),'tab_nb','get_tab_nb')

  DO i=1,size(tab_l)
    tab_nb(i) = get_nb_FROM_basis(tab_ba(tab_l(i),i))
  END DO

END FUNCTION get_tab_nb
FUNCTION getbis_tab_nb(tab_l,tab_ba) RESULT (tab_nb)
USE mod_basis_set_alloc
IMPLICIT NONE

  integer ,                        intent(in)            :: tab_l(:)
  TYPE(basis),                     intent(in)            :: tab_ba(0:,:) ! tab_ba(0:L,D)
  integer                                                :: tab_nb(size(tab_l)) ! result

  integer :: i

  DO i=1,size(tab_l)
    tab_nb(i) = get_nb_FROM_basis(tab_ba(tab_l(i),i))
  END DO

END FUNCTION getbis_tab_nb

END MODULE mod_basis_BtoG_GtoB_SGType4
