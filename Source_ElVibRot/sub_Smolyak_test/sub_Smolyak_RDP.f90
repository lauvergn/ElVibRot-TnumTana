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
MODULE mod_Smolyak_RDP
use mod_system
IMPLICIT NONE


TYPE TypeRDP
  integer :: n1 = 0
  integer :: n2 = 0
  integer :: n3 = 0
  real(kind=Rkind), allocatable :: RDP(:,:,:)
END TYPE TypeRDP

TYPE TypeRVec
  real(kind=Rkind), allocatable :: R(:)
END TYPE TypeRVec

TYPE Type_SmolyakRep
  logical :: Grid     = .FALSE.
  logical :: Delta    = .FALSE.
  integer :: k        = -1
  logical :: k_type_b = .TRUE.
  TYPE (TypeRVec), allocatable :: SmolyakRep(:)
END TYPE Type_SmolyakRep

INTERFACE assignment(=)
  module procedure TypeRVec2_TO_TypeRVec1,tabR2_TO_TypeRVec1
  module procedure SmolyakRep2_TO_tabR1,tabR2_TO_SmolyakRep1,R2_TO_SmolyakRep1
END INTERFACE

INTERFACE operator(*)
  module procedure SmolyakRep1_TIME_SmolyakRe2
END INTERFACE
INTERFACE operator(+)
  module procedure SmolyakRep1_PLUS_SmolyakRe2
END INTERFACE
INTERFACE operator(-)
  module procedure SmolyakRep1_MINUS_SmolyakRe2
END INTERFACE
CONTAINS

SUBROUTINE alloc_TypeRVec(Rvec,nvec)
IMPLICIT NONE

  TYPE (TypeRVec), intent(inout) :: Rvec
  integer,         intent(in)    :: nvec


  CALL dealloc_TypeRVec(Rvec)

  IF (nvec < 1) THEN
    write(out_unitp,*) ' ERROR in alloc_TypeRVec'
    write(out_unitp,*) ' nvec < 1',nvec
    STOP
  END IF

  allocate(Rvec%R(nvec))

END SUBROUTINE alloc_TypeRVec
SUBROUTINE dealloc_TypeRVec(Rvec)
!USE mod_system
IMPLICIT NONE

  TYPE (TypeRVec), intent(inout) :: Rvec

  IF (allocated(Rvec%R)) THEN
    deallocate(Rvec%R)
  END IF

END SUBROUTINE dealloc_TypeRVec

SUBROUTINE Write_TypeRVec(Rvec)
!USE mod_system
IMPLICIT NONE

  TYPE (TypeRVec), intent(in) :: Rvec

  IF (allocated(Rvec%R)) THEN
    write(out_unitp,*) 'R:',Rvec%R
  END IF

END SUBROUTINE Write_TypeRVec

SUBROUTINE TypeRVec2_TO_TypeRVec1(Rvec1,Rvec2)
IMPLICIT NONE

  TYPE (TypeRVec), intent(inout) :: Rvec1
  TYPE (TypeRVec), intent(in)    :: Rvec2

  CALL dealloc_TypeRVec(Rvec1)

  IF (allocated(Rvec2%R)) Rvec1%R = Rvec2%R

END SUBROUTINE TypeRVec2_TO_TypeRVec1
SUBROUTINE tabR2_TO_TypeRVec1(Rvec1,tabR2)
IMPLICIT NONE

  TYPE (TypeRVec),                intent(inout) :: Rvec1
  real(kind=Rkind), allocatable,  intent(in)    :: tabR2(:)

  CALL dealloc_TypeRVec(Rvec1)

  IF (allocated(tabR2)) Rvec1%R = tabR2

END SUBROUTINE tabR2_TO_TypeRVec1


SUBROUTINE dealloc_TabRDP(TabRDP)
IMPLICIT NONE

TYPE(TypeRDP), allocatable :: TabRDP(:)


integer               :: iG


IF (allocated(TabRDP)) THEN
  DO iG=lbound(TabRDP,dim=1),ubound(TabRDP,dim=1)

    TabRDP(iG)%n1 = 0
    TabRDP(iG)%n2 = 0
    TabRDP(iG)%n3 = 0

    deallocate(TabRDP(iG)%RDP)

  END DO

  deallocate(TabRDP)

END IF


END SUBROUTINE dealloc_TabRDP

SUBROUTINE alloc_SmolyakRep(SRep,tab_ind,tab_ba,delta,grid)
use mod_Smolyak_ba, only: typeba, get_tab_nq, get_tab_nb
IMPLICIT NONE

TYPE(Type_SmolyakRep),           intent(inout)         :: SRep
integer,         allocatable,    intent(in)            :: tab_ind(:,:) ! tab_ind(D,MaxnD)
TYPE(TypeBa),    allocatable,    intent(in)            :: tab_ba(:,:) ! tab_ba(L,D)
logical,                         intent(in),  optional :: delta,grid

integer               :: i,iG,nb_B,D,MaxnD
integer, allocatable  :: tab_n(:)

CALL dealloc_SmolyakRep(SRep)

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

D = size(tab_ind(:,1))
MaxnD = size(tab_ind(1,:))
!write(out_unitp,*) 'Alloc Smolyak Rep'


allocate(SRep%SmolyakRep(MaxnD))

DO iG=1,MaxnD
  IF (SRep%grid) THEN
    tab_n = get_tab_nq(tab_ind(:,iG),tab_ba)
  ELSE
    tab_n = get_tab_nb(tab_ind(:,iG),tab_ba)
  END IF
  !write(out_unitp,*) iG,'tab_n',tab_n
  allocate(SRep%SmolyakRep(iG)%R(product(tab_n)))
  !allocate(SRep%SmolyakRep(iG)%R(tab_n(1),tab_n(2)))

END DO

IF (allocated(tab_n)) deallocate(tab_n)

!nb_B = Size_SmolyakRep(SRep)
!write(out_unitp,*) 'Size Smolyak Rep:',nb_B

END SUBROUTINE alloc_SmolyakRep
SUBROUTINE dealloc_SmolyakRep(SRep)
IMPLICIT NONE

TYPE(Type_SmolyakRep), intent(inout)     :: SRep

integer               :: iG

IF (allocated(SRep%SmolyakRep)) THEN
  DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)
    CALL dealloc_TypeRVec(SRep%SmolyakRep(iG))
  END DO
  deallocate(SRep%SmolyakRep)
END IF

  SRep%k        = -1
  SRep%k_type_b = .TRUE.

  SRep%Grid     = .FALSE.
  SRep%Delta    = .FALSE.

END SUBROUTINE dealloc_SmolyakRep
SUBROUTINE Write_TabRDP(TabRDP)
IMPLICIT NONE

TYPE(TypeRDP), allocatable :: TabRDP(:)


integer               :: iG,i1,i2,i3


IF (allocated(TabRDP)) THEN
  write(out_unitp,*) '======== TabRDP ============================',shape(TabRDP)
  DO iG=lbound(TabRDP,dim=1),ubound(TabRDP,dim=1)
    write(out_unitp,*) '    iG',iG
    DO i3=1,TabRDP(iG)%n3
    DO i2=1,TabRDP(iG)%n2
    DO i1=1,TabRDP(iG)%n1
      write(out_unitp,*) 'i1,i2,i3',i1,i2,i3,TabRDP(iG)%RDP(i1,i2,i3)
    END DO
    END DO
    END DO

  END DO

END IF


END SUBROUTINE Write_TabRDP
SUBROUTINE Write_TabRDP_pack(TabRDP)
IMPLICIT NONE

TYPE(TypeRDP), allocatable :: TabRDP(:)


integer               :: igg,iG,i1,i2,i3

IF (allocated(TabRDP)) THEN
  igg = 0
  DO iG=lbound(TabRDP,dim=1),ubound(TabRDP,dim=1)
    DO i3=1,TabRDP(iG)%n3
    DO i2=1,TabRDP(iG)%n2
    DO i1=1,TabRDP(iG)%n1
      igg = igg+1
      write(out_unitp,*) igg,TabRDP(iG)%RDP(i1,i2,i3)
    END DO
    END DO
    END DO

  END DO

END IF


END SUBROUTINE Write_TabRDP_pack


SUBROUTINE Write_SmolyakRep(SRep)
IMPLICIT NONE

TYPE(Type_SmolyakRep), intent(in)     :: SRep

integer               :: iG

  write(out_unitp,*) 'Grid',SRep%Grid
  write(out_unitp,*) 'Delta',SRep%Delta
  write(out_unitp,*) 'k',SRep%k
  write(out_unitp,*) 'k_type_b',SRep%k_type_b
  write(out_unitp,*) 'alloc?',allocated(SRep%SmolyakRep)

IF (allocated(SRep%SmolyakRep)) THEN
  write(out_unitp,*) '======== Smolyak Rep ============================'
  DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)
    write(out_unitp,*) iG,size(SRep%SmolyakRep(iG)%R),SRep%SmolyakRep(iG)%R
  END DO
END IF

END SUBROUTINE Write_SmolyakRep
SUBROUTINE Write_SmolyakRep_pack(SRep)
IMPLICIT NONE

TYPE(Type_SmolyakRep), intent(in)     :: SRep

integer               :: i,iG,igg

real(kind=Rkind), allocatable :: R(:)


IF (allocated(SRep%SmolyakRep)) THEN
  igg = 0
  DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)
    R = reshape(SRep%SmolyakRep(iG)%R,(/ size(SRep%SmolyakRep(iG)%R) /) )
    DO i=1,size(SRep%SmolyakRep(iG)%R)
      igg = igg+1
      !write(out_unitp,*) igg,SRep%SmolyakRep(iG)%R(i)
      write(out_unitp,*) igg,R(i)
    END DO
  END DO
END IF

IF (allocated(R)) deallocate(R)

END SUBROUTINE Write_SmolyakRep_pack
SUBROUTINE SumSq_TabRDP(TabRDP)
IMPLICIT NONE

TYPE(TypeRDP), allocatable :: TabRDP(:)

integer               :: iG,i1,i2,i3
real (kind=Rkind)     :: SS,S


SS = ZERO
IF (allocated(TabRDP)) THEN
  DO iG=lbound(TabRDP,dim=1),ubound(TabRDP,dim=1)
    S  = sum(TabRDP(iG)%RDP**2)
    SS = SS + S

!    DO i3=1,TabRDP(iG)%n3
!    DO i2=1,TabRDP(iG)%n2
!    DO i1=1,TabRDP(iG)%n1
!      IF (abs(TabRDP(iG)%RDP(i1,i2,i3))>1.d-5) write(out_unitp,*) 'i1,i2,i3',i1,i2,i3,TabRDP(iG)%RDP(i1,i2,i3)
!    END DO
!    END DO
!    END DO

  END DO
  write(out_unitp,*) 'SumSq TabRDP',SS

END IF


END SUBROUTINE SumSq_TabRDP

SUBROUTINE Size_TabRDP(TabRDP,nb_BG)
IMPLICIT NONE

TYPE(TypeRDP), allocatable :: TabRDP(:)
integer  :: nb_BG

integer               :: iG


nb_BG = 0
IF (allocated(TabRDP)) THEN
  DO iG=lbound(TabRDP,dim=1),ubound(TabRDP,dim=1)
    nb_BG = nb_BG + size(TabRDP(iG)%RDP)
  END DO
END IF


END SUBROUTINE Size_TabRDP

FUNCTION Size_SmolyakRep(SRep) RESULT(nb_BG)
IMPLICIT NONE

TYPE(Type_SmolyakRep), intent(in)     :: SRep
integer                               :: nb_BG

integer               :: iG


nb_BG = 0
IF (allocated(SRep%SmolyakRep)) THEN
  DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)
    nb_BG = nb_BG + size(SRep%SmolyakRep(iG)%R)
  END DO
END IF

END FUNCTION Size_SmolyakRep
FUNCTION MaxVal_SmolyakRep(SRep) RESULT(maxSRep)
IMPLICIT NONE

TYPE(Type_SmolyakRep), intent(in)     :: SRep
real (kind=Rkind)                     :: maxSRep

integer               :: iG

maxSRep = ZERO
IF (allocated(SRep%SmolyakRep)) THEN
  DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)
    maxSRep = max(maxSRep,maxval(abs(SRep%SmolyakRep(iG)%R)))
  END DO
END IF

END FUNCTION MaxVal_SmolyakRep
SUBROUTINE SmolyakRep2_TO_tabR1(tabR1,SRep2)
IMPLICIT NONE

real(kind=Rkind), allocatable,   intent(inout)  :: tabR1(:)
TYPE(Type_SmolyakRep),           intent(in)     :: SRep2

integer               :: iG,nb_BG,nR,itabR

IF (allocated(tabR1)) deallocate(tabR1)

nb_BG = Size_SmolyakRep(SRep2)

IF (nb_BG > 0) THEN
  allocate(tabR1(nb_BG))

  itabR = 0
  DO iG=lbound(SRep2%SmolyakRep,dim=1),ubound(SRep2%SmolyakRep,dim=1)
    nR = size(SRep2%SmolyakRep(iG)%R)
    tabR1(itabR+1:itabR+nR) = reshape(SRep2%SmolyakRep(iG)%R,shape=(/nR/))
    itabR = itabR + nR
  END DO

END IF

END SUBROUTINE SmolyakRep2_TO_tabR1
SUBROUTINE tabR2_TO_SmolyakRep1(SRep1,tabR2)
IMPLICIT NONE

real(kind=Rkind), allocatable,   intent(in)     :: tabR2(:)
TYPE(Type_SmolyakRep),           intent(inout)  :: SRep1

integer               :: iG,nb_BG,nR,itabR


nb_BG = Size_SmolyakRep(SRep1)
IF (size(tabR2) /= nb_BG) THEN
  write(out_unitp,*) ' ERROR in tabR2_TO_SmolyakRep1'
  write(out_unitp,*) ' sizes are different!!'
  write(out_unitp,*) ' sizes of tabR2 and SRep1',size(tabR2),nb_BG
  STOP
END IF

IF (nb_BG > 0) THEN

  itabR = 0
  DO iG=lbound(SRep1%SmolyakRep,dim=1),ubound(SRep1%SmolyakRep,dim=1)
    nR = size(SRep1%SmolyakRep(iG)%R)
    !SRep1%SmolyakRep(iG)%R = tabR2(itabR+1:itabR+nR)
    SRep1%SmolyakRep(iG)%R = reshape(tabR2(itabR+1:itabR+nR),shape=shape(SRep1%SmolyakRep(iG)%R))
    itabR = itabR + nR
  END DO

END IF

END SUBROUTINE tabR2_TO_SmolyakRep1
SUBROUTINE R2_TO_SmolyakRep1(SRep1,R2)
IMPLICIT NONE

real(kind=Rkind),                intent(in)     :: R2
TYPE(Type_SmolyakRep),           intent(inout)  :: SRep1

integer               :: iG,nb_BG


nb_BG = Size_SmolyakRep(SRep1)

IF (nb_BG > 0) THEN

  DO iG=lbound(SRep1%SmolyakRep,dim=1),ubound(SRep1%SmolyakRep,dim=1)
    SRep1%SmolyakRep(iG)%R = R2
  END DO

END IF

END SUBROUTINE R2_TO_SmolyakRep1

SUBROUTINE R2_TO_SmolyakRep1_with_tab_i(SRep1,R2,tab_i,tab_ba,tab_ind)
USE mod_Smolyak_DInd, ONLY : l_TO_n
USE mod_Smolyak_ba, only: typeba, get_tab_nq, get_tab_nb
IMPLICIT NONE

real(kind=Rkind),                intent(in)             :: R2
TYPE(Type_SmolyakRep),           intent(inout)          :: SRep1
integer,                         intent(in)             :: tab_i(:)
TYPE(TypeBa),    allocatable,    intent(in)             :: tab_ba(:,:) ! tab_ba(L,D)
integer,         allocatable,    intent(in)             :: tab_ind(:,:) ! tab_ind(D,MaxnD)

integer               :: i,II,iG,nb_BG,D
integer, allocatable  :: tab_n(:),tab_n1(:),tab_n2(:)


nb_BG = Size_SmolyakRep(SRep1)
D     = size(tab_ind(:,1))

IF (SRep1%Grid) STOP 'Grid is not possible in R2_TO_SmolyakRep1_with_tab_i'

IF (nb_BG > 0) THEN
IF (SRep1%delta) THEN
  allocate(tab_n1(D))
  allocate(tab_n2(D))
  DO iG=lbound(SRep1%SmolyakRep,dim=1),ubound(SRep1%SmolyakRep,dim=1)
    tab_n = tab_ind(:,iG)
    DO i=1,D
      tab_n2(i) = l_TO_n(tab_n(i),1)
      tab_n1(i) = l_TO_n(tab_n(i)-1,1)
    END DO

    !write(out_unitp,*) 'iG',iG
    !write(out_unitp,*) 'tab_i',tab_i
    !write(out_unitp,*) 'tab_l',tab_n
    !write(out_unitp,*) 'tab_n1',tab_n1
    !write(out_unitp,*) 'tab_n2',tab_n2

    IF (minval(tab_n2 - tab_i) < 0) CYCLE
    IF (minval(tab_i - tab_n1) < 1) CYCLE

    II = tab_i(D)-tab_n1(D)
    DO i=D-1,1,-1
      II = (II-1)*(tab_n2(i)-tab_n1(i)) + tab_i(i)-tab_n1(i)
    END DO

    !write(out_unitp,*) 'II size R',II,size(SRep1%SmolyakRep(iG)%R)

    IF (II > size(SRep1%SmolyakRep(iG)%R) .OR. II < 1) THEN
      write(out_unitp,*) 'iG',iG
      write(out_unitp,*) 'II size R',II,size(SRep1%SmolyakRep(iG)%R)
      write(out_unitp,*) 'tab_i',tab_i
      write(out_unitp,*) 'tab_l',tab_n
      write(out_unitp,*) 'tab_n1',tab_n1
      write(out_unitp,*) 'tab_n2',tab_n2

      STOP 'ERROR in R2_TO_SmolyakRep1_with_tab_i'
    END IF

    SRep1%SmolyakRep(iG)%R(II) = R2

  END DO
  IF (allocated(tab_n)) deallocate(tab_n)
  IF (allocated(tab_n1)) deallocate(tab_n1)
  IF (allocated(tab_n2)) deallocate(tab_n2)

ELSE
  DO iG=lbound(SRep1%SmolyakRep,dim=1),ubound(SRep1%SmolyakRep,dim=1)

    tab_n = get_tab_nb(tab_ind(:,iG),tab_ba)

    IF (minval(tab_n - tab_i) < 0) CYCLE

    II = tab_i(D)
    DO i=D-1,1,-1
      II = (II-1)*tab_n(i) + tab_i(i)
    END DO

    IF (II > size(SRep1%SmolyakRep(iG)%R)) THEN
      write(out_unitp,*) 'tab_n',tab_n
      write(out_unitp,*) 'tab_i',tab_i
      write(out_unitp,*) 'II',II
      STOP 'ERROR in R2_TO_SmolyakRep1_with_tab_i'
    END IF

    SRep1%SmolyakRep(iG)%R( II ) = R2

  END DO
  IF (allocated(tab_n)) deallocate(tab_n)
END IF
END IF

END SUBROUTINE R2_TO_SmolyakRep1_with_tab_i

FUNCTION dot_product_SmolyakRep(SRep1,SRep2,WSRep) RESULT(R)
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
    R = R + WSRep(iG) * dot_product(SRep1%SmolyakRep(iG)%R,SRep2%SmolyakRep(iG)%R)
  END DO

END IF

END FUNCTION dot_product_SmolyakRep
FUNCTION SmolyakRep1_TIME_SmolyakRe2(SRep1,SRep2) RESULT(SRep)
IMPLICIT NONE

TYPE(Type_SmolyakRep),           intent(in)     :: SRep1,SRep2
TYPE(Type_SmolyakRep)                           :: SRep


integer               :: iG,nb_BG

nb_BG = Size_SmolyakRep(SRep1)
IF (nb_BG /= Size_SmolyakRep(SRep2)) THEN
  write(out_unitp,*) 'ERROR in SmolyakRep1_TIME_SmolyakRe2'
  write(out_unitp,*) 'sizes are different',Size_SmolyakRep(SRep1),Size_SmolyakRep(SRep2)
  STOP
END IF

SRep = SRep1 !  for the allocation

IF (nb_BG > 0) THEN

  DO iG=lbound(SRep1%SmolyakRep,dim=1),ubound(SRep1%SmolyakRep,dim=1)
    SRep%SmolyakRep(iG)%R = SRep1%SmolyakRep(iG)%R * SRep2%SmolyakRep(iG)%R
  END DO

END IF

END FUNCTION SmolyakRep1_TIME_SmolyakRe2
FUNCTION SmolyakRep1_PLUS_SmolyakRe2(SRep1,SRep2) RESULT(SRep)
IMPLICIT NONE

TYPE(Type_SmolyakRep),           intent(in)     :: SRep1,SRep2
TYPE(Type_SmolyakRep)                           :: SRep


integer               :: iG,nb_BG

nb_BG = Size_SmolyakRep(SRep1)
IF (nb_BG /= Size_SmolyakRep(SRep2)) THEN
  write(out_unitp,*) 'ERROR in SmolyakRep1_PLUS_SmolyakRe2'
  write(out_unitp,*) 'sizes are different',Size_SmolyakRep(SRep1),Size_SmolyakRep(SRep2)
  STOP
END IF

SRep = SRep1 !  for the allocation

IF (nb_BG > 0) THEN

  DO iG=lbound(SRep1%SmolyakRep,dim=1),ubound(SRep1%SmolyakRep,dim=1)
    SRep%SmolyakRep(iG)%R = SRep1%SmolyakRep(iG)%R + SRep2%SmolyakRep(iG)%R
  END DO

END IF

END FUNCTION SmolyakRep1_PLUS_SmolyakRe2
FUNCTION SmolyakRep1_MINUS_SmolyakRe2(SRep1,SRep2) RESULT(SRep)
IMPLICIT NONE

TYPE(Type_SmolyakRep),           intent(in)     :: SRep1,SRep2
TYPE(Type_SmolyakRep)                           :: SRep


integer               :: iG,nb_BG

nb_BG = Size_SmolyakRep(SRep1)
IF (nb_BG /= Size_SmolyakRep(SRep2)) THEN
  write(out_unitp,*) 'ERROR in SmolyakRep1_MINUS_SmolyakRe2'
  write(out_unitp,*) 'sizes are different',Size_SmolyakRep(SRep1),Size_SmolyakRep(SRep2)
  STOP
END IF

SRep = SRep1 !  for the allocation

IF (nb_BG > 0) THEN

  DO iG=lbound(SRep1%SmolyakRep,dim=1),ubound(SRep1%SmolyakRep,dim=1)
    SRep%SmolyakRep(iG)%R = SRep1%SmolyakRep(iG)%R - SRep2%SmolyakRep(iG)%R
  END DO

END IF

END FUNCTION SmolyakRep1_MINUS_SmolyakRe2
SUBROUTINE GSmolyakRep_TO_BSmolyakRep(SRep,tab_ind,tab_ba)
USE mod_Smolyak_ba, only: typeba, get_tab_nq, get_tab_nb
IMPLICIT NONE

real(kind=Rkind)  :: R
TYPE(Type_SmolyakRep),           intent(inout)          :: SRep
integer,         allocatable,    intent(in)             :: tab_ind(:,:) ! tab_ind(D,MaxnD)
TYPE(TypeBa),    allocatable,    intent(in)             :: tab_ba(:,:) ! tab_ba(L,D)

integer               :: i,D,iG,nb_BG

integer                            :: nnb,nnq,nb2,nq2,ib,iq
integer, allocatable               :: tab_nb(:),tab_nq(:)
real(kind=Rkind), allocatable      :: RTempG(:,:,:),RTempB(:,:,:)

IF (.NOT. SRep%Grid) STOP 'Grid is not possible in GSmolyakRep_TO_BSmolyakRep'

nb_mult_GTOB = 0

D = size(tab_ind(:,1))

!$OMP   PARALLEL DEFAULT(NONE) &
!$OMP   SHARED(D,SRep,tab_ind,tab_ba,nb_mult_GTOB) &
!$OMP   PRIVATE(iG,tab_nb,tab_nq,i,ib,iq,nnb,nnq,nb2,nq2,RTempG,RTempB) &
!$OMP   NUM_THREADS(BasisTOGrid_maxth)

!$OMP   DO SCHEDULE(STATIC)

DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)

  tab_nb = get_tab_nb(tab_ind(:,iG),tab_ba)
  tab_nq = get_tab_nq(tab_ind(:,iG),tab_ba)

  nnq = product(tab_nq)
  nnb = 1
  nb2 = 1
  nq2 = 1

  RTempB = reshape(SRep%SmolyakRep(iG)%R,shape=(/ 1,1,nnq /))

! B order : b1 * b2 * b3 * ... bD ???
! G order : gD ... * g3 * g2 * g1 ???

  DO i=D,1,-1
    nb2 = tab_nb(i)
    nq2 = tab_nq(i)
    nnq = nnq / nq2

    RTempG = reshape(RTempB,shape=(/ nnb,nq2,nnq /))

    deallocate(RTempB)
    allocate(RTempB(nnb,nb2,nnq))

    !write(out_unitp,*) 'i',i
    !write(out_unitp,*) 'nnq,nq2,nb2,nnb',nnq,nq2,nb2,nnb
    !write(out_unitp,*) 'shape RTempB',shape(RTempB)
    !write(out_unitp,*) 'shape RTempG',shape(RTempG)

    DO iq=1,nnq
    DO ib=1,nnb
       !RTempB(ib,:,iq) = matmul(tab_ba(tab_l(i),i)%w*RTempG(ib,:,iq) , tab_ba(tab_l(i),i)%d0b)
       RTempB(ib,:,iq) = matmul(tab_ba(tab_ind(i,iG),i)%twd0b , RTempG(ib,:,iq))
       !$OMP ATOMIC
       nb_mult_GTOB = nb_mult_GTOB + nb2*nq2
    END DO
    END DO

    nnb = nnb * nb2
    deallocate(RTempG)

  END DO

  SRep%SmolyakRep(iG)%R = reshape(RTempB, shape=(/ nnb /) )
  !SRep%SmolyakRep(iG)%R = reshape(RTempB, (/ tab_nb(1),tab_nb(2) /) )

  deallocate(RTempB)
END DO
!$OMP   END DO

IF (allocated(tab_nb)) deallocate(tab_nb)
IF (allocated(tab_nq)) deallocate(tab_nq)


!$OMP   END PARALLEL

SRep%Grid = .FALSE.

END SUBROUTINE GSmolyakRep_TO_BSmolyakRep
SUBROUTINE GSmolyakRep_TO_BSmolyakRep_01(SRep,tab_ind,tab_ba)
USE mod_Smolyak_ba, only: typeba, get_tab_nq, get_tab_nb
IMPLICIT NONE

real(kind=Rkind)  :: R
TYPE(Type_SmolyakRep),           intent(inout)          :: SRep
integer,         allocatable,    intent(in)             :: tab_ind(:,:) ! tab_ind(D,MaxnD)
TYPE(TypeBa),    allocatable,    intent(in)             :: tab_ba(:,:) ! tab_ba(L,D)

integer               :: i,D,iG,nb_BG

integer                            :: nnb,nnq,nb2,nq2,ib,iq
integer, allocatable               :: tab_nb(:),tab_nq(:)
real(kind=Rkind), allocatable      :: RTempG(:,:,:),RTempB(:,:,:)

IF (.NOT. SRep%Grid) STOP 'Grid is not possible in GSmolyakRep_TO_BSmolyakRep_01'

nb_mult_GTOB = 0

D = size(tab_ind(:,1))

!$OMP   PARALLEL DEFAULT(NONE) &
!$OMP   SHARED(D,SRep,tab_ind,tab_ba,nb_mult_GTOB) &
!$OMP   PRIVATE(iG,tab_nb,tab_nq,i,ib,iq,nnb,nnq,nb2,nq2,RTempG,RTempB) &
!$OMP   NUM_THREADS(BasisTOGrid_maxth)

!$OMP   DO SCHEDULE(STATIC)

DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)

  tab_nb = get_tab_nb(tab_ind(:,iG),tab_ba)
  tab_nq = get_tab_nq(tab_ind(:,iG),tab_ba)

  nnq = product(tab_nq)
  nnb = 1
  nb2 = 1
  nq2 = 1

  RTempB = reshape(SRep%SmolyakRep(iG)%R,shape=(/ 1,1,nnq /))

! B order : b1 * b2 * b3 * ... bD ???
! G order : gD ... * g3 * g2 * g1 ???

  DO i=1,D
    nb2 = tab_nb(i)
    nq2 = tab_nq(i)
    nnq = nnq / nq2

    RTempG = reshape(RTempB,shape=(/ nnb,nq2,nnq /))

    deallocate(RTempB)
    allocate(RTempB(nnb,nb2,nnq))

    !write(out_unitp,*) 'i',i
    !write(out_unitp,*) 'nnq,nq2,nb2,nnb',nnq,nq2,nb2,nnb
    !write(out_unitp,*) 'shape RTempB',shape(RTempB)
    !write(out_unitp,*) 'shape RTempG',shape(RTempG)

    DO iq=1,nnq
    DO ib=1,nnb
       RTempB(ib,:,iq) = matmul(tab_ba(tab_ind(i,iG),i)%twd0b , RTempG(ib,:,iq))
       !$OMP ATOMIC
       nb_mult_GTOB = nb_mult_GTOB + nb2*nq2
    END DO
    END DO

    nnb = nnb * nb2
    deallocate(RTempG)

  END DO

  SRep%SmolyakRep(iG)%R = reshape(RTempB, shape=(/ nnb /) )
  !SRep%SmolyakRep(iG)%R = reshape(RTempB, (/ tab_nb(1),tab_nb(2) /) )

  deallocate(RTempB)
END DO
!$OMP   END DO

IF (allocated(tab_nb)) deallocate(tab_nb)
IF (allocated(tab_nq)) deallocate(tab_nq)


!$OMP   END PARALLEL

SRep%Grid = .FALSE.

END SUBROUTINE GSmolyakRep_TO_BSmolyakRep_01


SUBROUTINE GSmolyakRep_TO_BSmolyakRep_old(SRep,tab_ind,tab_ba)
USE mod_Smolyak_ba, only: typeba, get_tab_nq, get_tab_nb
IMPLICIT NONE

real(kind=Rkind)  :: R
TYPE(Type_SmolyakRep),           intent(inout)          :: SRep
integer,         allocatable,    intent(in)             :: tab_ind(:,:) ! tab_ind(D,MaxnD)
TYPE(TypeBa),    allocatable,    intent(in)             :: tab_ba(:,:) ! tab_ba(L,D)

integer               :: i,D,iG,nb_BG

integer                            :: nnb,nnq,nb2,nq2,ib,iq
integer, allocatable               :: tab_nb(:),tab_nq(:)
real(kind=Rkind), allocatable      :: RTempG(:,:,:),RTempB(:,:,:)

IF (.NOT. SRep%Grid) STOP 'Grid is not possible in GSmolyakRep_TO_BSmolyakRep_old'

nb_mult_GTOB = 0

D = size(tab_ind(:,1))

!$OMP   PARALLEL DEFAULT(NONE) &
!$OMP   SHARED(D,SRep,tab_ind,tab_ba,nb_mult_GTOB) &
!$OMP   PRIVATE(iG,tab_nb,tab_nq,i,ib,iq,nnb,nnq,nb2,nq2,RTempG,RTempB) &
!$OMP   NUM_THREADS(BasisTOGrid_maxth)

!$OMP   DO SCHEDULE(STATIC)

DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)

  tab_nb = get_tab_nb(tab_ind(:,iG),tab_ba)
  tab_nq = get_tab_nq(tab_ind(:,iG),tab_ba)

  nnq = product(tab_nq)
  nnb = 1
  nb2 = 1
  nq2 = 1

  RTempB = reshape(SRep%SmolyakRep(iG)%R,shape=(/ nnq,1,1 /))

  DO i=1,D
    nb2 = tab_nb(i)
    nq2 = tab_nq(i)
    nnq = nnq / nq2

    RTempG = reshape(RTempB,shape=(/ nnq,nq2,nnb /))

    deallocate(RTempB)
    allocate(RTempB(nnq,nb2,nnb))

    DO ib=1,nnb
    DO iq=1,nnq
       !RTempB(iq,:,ib) = matmul(tab_ba(tab_l(i),i)%w*RTempG(iq,:,ib) , tab_ba(tab_l(i),i)%d0b)
       RTempB(iq,:,ib) = matmul(tab_ba(tab_ind(i,iG),i)%twd0b , RTempG(iq,:,ib))
       !$OMP ATOMIC
       nb_mult_GTOB = nb_mult_GTOB + nb2*nq2
    END DO
    END DO

    nnb = nnb * nb2
    deallocate(RTempG)

  END DO

  SRep%SmolyakRep(iG)%R = reshape(RTempB, shape=(/ nnb /) )
  !SRep%SmolyakRep(iG)%R = reshape(RTempB, (/ tab_nb(1),tab_nb(2) /) )

  deallocate(RTempB)
END DO
!$OMP   END DO

IF (allocated(tab_nb)) deallocate(tab_nb)
IF (allocated(tab_nq)) deallocate(tab_nq)


!$OMP   END PARALLEL

SRep%Grid = .FALSE.

END SUBROUTINE GSmolyakRep_TO_BSmolyakRep_old

SUBROUTINE BSmolyakRep_TO_GSmolyakRep(SRep,tab_ind,tab_ba)
USE mod_Smolyak_ba, only: typeba, get_tab_nq, get_tab_nb
IMPLICIT NONE

real(kind=Rkind)  :: R
TYPE(Type_SmolyakRep),           intent(inout)          :: SRep
integer,         allocatable,    intent(in)             :: tab_ind(:,:) ! tab_ind(D,MaxnD)
TYPE(TypeBa),    allocatable,    intent(in)             :: tab_ba(:,:) ! tab_ba(L,D)

integer               :: i,D,iG,nb_BG

integer                            :: nnb,nnq,nb2,nq2,ib,iq
integer, allocatable               :: tab_nb(:),tab_nq(:)
real(kind=Rkind), allocatable      :: RTempG(:,:,:),RTempB(:,:,:)

IF (SRep%Grid) STOP 'Grid is not possible in BSmolyakRep_TO_GSmolyakRep'

nb_mult_BTOG = 0

D = size(tab_ind(:,1))

!$OMP   PARALLEL DEFAULT(NONE) &
!$OMP   SHARED(D,SRep,tab_ind,tab_ba,nb_mult_BTOG) &
!$OMP   PRIVATE(iG,tab_nb,tab_nq,i,ib,iq,nnb,nnq,nb2,nq2,RTempG,RTempB) &
!$OMP   NUM_THREADS(BasisTOGrid_maxth)

!$OMP   DO SCHEDULE(STATIC)

DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)

  tab_nb = get_tab_nb(tab_ind(:,iG),tab_ba)
  tab_nq = get_tab_nq(tab_ind(:,iG),tab_ba)

  nnb = product(tab_nb)
  nnq = 1
  nb2 = 1
  nq2 = 1

  RTempG = reshape(SRep%SmolyakRep(iG)%R,shape=(/ nnq,nq2,nnb /))
! B order : b1 * b2 * b3 * ... bD ???
! G order : gD ... * g3 * g2 * g1 ???

  DO i=D,1,-1
    nb2 = tab_nb(i)
    nq2 = tab_nq(i)
    nnb = nnb / nb2



    RTempB = reshape(RTempG,shape=(/ nnq,nb2,nnb /))

    deallocate(RTempG)
    allocate(RTempG(nnq,nq2,nnb))

    !write(out_unitp,*) 'i',i
    !write(out_unitp,*) 'nnq,nq2,nb2,nnb',nnq,nq2,nb2,nnb
    !write(out_unitp,*) 'shape RTempB',shape(RTempB)
    !write(out_unitp,*) 'shape RTempG',shape(RTempG)


    DO ib=1,nnb
    DO iq=1,nnq
       RTempG(iq,:,ib) = matmul(tab_ba(tab_ind(i,iG),i)%d0b,RTempB(iq,:,ib))
       !$OMP ATOMIC
       nb_mult_BTOG = nb_mult_BTOG + nb2*nq2
    END DO
    END DO

    nnq = nnq * nq2
    deallocate(RTempB)

  END DO

  SRep%SmolyakRep(iG)%R = reshape(RTempG, shape=(/ nnq*nnb /) )
  !SRep%SmolyakRep(iG)%R = reshape(RTempG, (/ tab_nq(1),tab_nq(2) /) )

  deallocate(RTempG)
END DO
!$OMP   END DO

IF (allocated(tab_nb)) deallocate(tab_nb)
IF (allocated(tab_nq)) deallocate(tab_nq)


!$OMP   END PARALLEL

SRep%Grid = .TRUE.


END SUBROUTINE BSmolyakRep_TO_GSmolyakRep
SUBROUTINE BSmolyakRep_TO_GSmolyakRep_01(SRep,tab_ind,tab_ba)
USE mod_Smolyak_ba, only: typeba, get_tab_nq, get_tab_nb
IMPLICIT NONE

real(kind=Rkind)  :: R
TYPE(Type_SmolyakRep),           intent(inout)          :: SRep
integer,         allocatable,    intent(in)             :: tab_ind(:,:) ! tab_ind(D,MaxnD)
TYPE(TypeBa),    allocatable,    intent(in)             :: tab_ba(:,:) ! tab_ba(L,D)

integer               :: i,D,iG,nb_BG

integer                            :: nnb,nnq,nb2,nq2,ib,iq
integer, allocatable               :: tab_nb(:),tab_nq(:)
real(kind=Rkind), allocatable      :: RTempG(:,:,:),RTempB(:,:,:)

IF (SRep%Grid) STOP 'Grid is not possible in BSmolyakRep_TO_GSmolyakRep_01'

nb_mult_BTOG = 0

D = size(tab_ind(:,1))

!$OMP   PARALLEL DEFAULT(NONE) &
!$OMP   SHARED(D,SRep,tab_ind,tab_ba,nb_mult_BTOG) &
!$OMP   PRIVATE(iG,tab_nb,tab_nq,i,ib,iq,nnb,nnq,nb2,nq2,RTempG,RTempB) &
!$OMP   NUM_THREADS(BasisTOGrid_maxth)

!$OMP   DO SCHEDULE(STATIC)

DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)

  tab_nb = get_tab_nb(tab_ind(:,iG),tab_ba)
  tab_nq = get_tab_nq(tab_ind(:,iG),tab_ba)

  nnb = product(tab_nb)
  nnq = 1
  nb2 = 1
  nq2 = 1

  RTempG = reshape(SRep%SmolyakRep(iG)%R,shape=(/ nnq,nq2,nnb /))
! B order : b1 * b2 * b3 * ... bD ???
! G order : gD ... * g3 * g2 * g1 ???

  DO i=1,D
    nb2 = tab_nb(i)
    nq2 = tab_nq(i)
    nnb = nnb / nb2



    RTempB = reshape(RTempG,shape=(/ nnq,nb2,nnb /))

    deallocate(RTempG)
    allocate(RTempG(nnq,nq2,nnb))

    !write(out_unitp,*) 'i',i
    !write(out_unitp,*) 'nnq,nq2,nb2,nnb',nnq,nq2,nb2,nnb
    !write(out_unitp,*) 'shape RTempB',shape(RTempB)
    !write(out_unitp,*) 'shape RTempG',shape(RTempG)


    DO ib=1,nnb
    DO iq=1,nnq
       !write(out_unitp,*) 'i,iG,ibb,iqq',i,iG,ib,iq,'b',RTempB(iq,:,ib)
       RTempG(iq,:,ib) = matmul(tab_ba(tab_ind(i,iG),i)%d0b,RTempB(iq,:,ib))
       !write(out_unitp,*) 'i,iG,ibb,iqq',i,iG,ib,iq,'g',RTempG(iq,:,ib)
       !$OMP ATOMIC
       nb_mult_BTOG = nb_mult_BTOG + nb2*nq2
    END DO
    END DO

    nnq = nnq * nq2
    deallocate(RTempB)

  END DO

  SRep%SmolyakRep(iG)%R = reshape(RTempG, shape=(/ nnq*nnb /) )
  !SRep%SmolyakRep(iG)%R = reshape(RTempG, (/ tab_nq(1),tab_nq(2) /) )

  deallocate(RTempG)
END DO
!$OMP   END DO

IF (allocated(tab_nb)) deallocate(tab_nb)
IF (allocated(tab_nq)) deallocate(tab_nq)


!$OMP   END PARALLEL

SRep%Grid = .TRUE.

END SUBROUTINE BSmolyakRep_TO_GSmolyakRep_01
SUBROUTINE BSmolyakRep_TO_GSmolyakRep_01_v2(SRep,tab_ind,tab_ba)
USE mod_Smolyak_ba, only: typeba, get_tab_nq, get_tab_nb
IMPLICIT NONE

real(kind=Rkind)  :: R
TYPE(Type_SmolyakRep),           intent(inout)          :: SRep
integer,         allocatable,    intent(in)             :: tab_ind(:,:) ! tab_ind(D,MaxnD)
TYPE(TypeBa),    allocatable,    intent(in)             :: tab_ba(:,:) ! tab_ba(L,D)

integer               :: i,D,iG,nb_BG

integer                            :: nnb,nnq,nb2,nq2,ib,iq
integer, allocatable               :: tab_nb(:),tab_nq(:)
real(kind=Rkind), allocatable      :: RTempG(:,:,:),RTempB(:,:,:)

IF (SRep%Grid) STOP 'Grid is not possible in BSmolyakRep_TO_GSmolyakRep_01'

nb_mult_BTOG = 0

D = size(tab_ind(:,1))

!$OMP   PARALLEL DEFAULT(NONE) &
!$OMP   SHARED(D,SRep,tab_ind,tab_ba,nb_mult_BTOG) &
!$OMP   PRIVATE(iG,tab_nb,tab_nq,i,ib,iq,nnb,nnq,nb2,nq2,RTempG,RTempB) &
!$OMP   NUM_THREADS(BasisTOGrid_maxth)

!$OMP   DO SCHEDULE(STATIC)

DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)

  tab_nb = get_tab_nb(tab_ind(:,iG),tab_ba)
  tab_nq = get_tab_nq(tab_ind(:,iG),tab_ba)

  nnb = product(tab_nb)
  nnq = 1
  nb2 = 1
  nq2 = 1

  RTempG = reshape(SRep%SmolyakRep(iG)%R,shape=(/ nnq,nq2,nnb /))
! B order : b1 * b2 * b3 * ... bD ???
! G order : gD ... * g3 * g2 * g1 ???

  DO i=1,D
    nb2 = tab_nb(i)
    nq2 = tab_nq(i)
    nnb = nnb / nb2



    RTempB = reshape(RTempG,shape=(/ nnq,nb2,nnb /))

    deallocate(RTempG)
    allocate(RTempG(nnq,nq2,nnb))

    !write(out_unitp,*) 'i',i
    !write(out_unitp,*) 'nnq,nq2,nb2,nnb',nnq,nq2,nb2,nnb
    !write(out_unitp,*) 'shape RTempB',shape(RTempB)
    !write(out_unitp,*) 'shape RTempG',shape(RTempG)


    DO ib=1,nnb
    DO iq=1,nnq
       !write(out_unitp,*) 'i,iG,ibb,iqq',i,iG,ib,iq,'b',RTempB(iq,:,ib)
       !RTempG(iq,:,ib) = matmul(tab_ba(tab_ind(i,iG),i)%d0b,RTempB(iq,:,ib))
       RTempG(iq,:,ib) = matmul(RTempB(iq,:,ib),tab_ba(tab_ind(i,iG),i)%td0b)

       !write(out_unitp,*) 'i,iG,ibb,iqq',i,iG,ib,iq,'g',RTempG(iq,:,ib)
       !$OMP ATOMIC
       nb_mult_BTOG = nb_mult_BTOG + nb2*nq2
    END DO
    END DO

    nnq = nnq * nq2
    deallocate(RTempB)

  END DO

  SRep%SmolyakRep(iG)%R = reshape(RTempG, shape=(/ nnq*nnb /) )
  !SRep%SmolyakRep(iG)%R = reshape(RTempG, (/ tab_nq(1),tab_nq(2) /) )

  deallocate(RTempG)
END DO
!$OMP   END DO

IF (allocated(tab_nb)) deallocate(tab_nb)
IF (allocated(tab_nq)) deallocate(tab_nq)


!$OMP   END PARALLEL

SRep%Grid = .TRUE.

END SUBROUTINE BSmolyakRep_TO_GSmolyakRep_01_v2
SUBROUTINE BSmolyakRep_TO_GSmolyakRep_withsqrtW(SRep,tab_ind,tab_ba)
USE mod_Smolyak_ba, only: typeba, get_tab_nq, get_tab_nb
IMPLICIT NONE

real(kind=Rkind)  :: R
TYPE(Type_SmolyakRep),           intent(inout)          :: SRep
integer,         allocatable,    intent(in)             :: tab_ind(:,:) ! tab_ind(D,MaxnD)
TYPE(TypeBa),    allocatable,    intent(in)             :: tab_ba(:,:) ! tab_ba(L,D)

integer               :: i,D,iG,nb_BG

integer                            :: nnb,nnq,nb2,nq2,ib,iq
integer, allocatable               :: tab_nb(:),tab_nq(:)
real(kind=Rkind), allocatable      :: RTempG(:,:,:),RTempB(:,:,:)

IF (SRep%Grid) STOP 'Grid is not possible in BSmolyakRep_TO_GSmolyakRep_withsqrtW'
IF (SRep%delta) STOP 'delta=t is not possible in BSmolyakRep_TO_GSmolyakRep_withsqrtW'

nb_mult_BTOG = 0

D = size(tab_ind(:,1))

!$OMP   PARALLEL DEFAULT(NONE) &
!$OMP   SHARED(D,SRep,tab_ind,tab_ba,nb_mult_BTOG) &
!$OMP   PRIVATE(iG,tab_nb,tab_nq,i,ib,iq,nnb,nnq,nb2,nq2,RTempG,RTempB) &
!$OMP   NUM_THREADS(BasisTOGrid_maxth)

!$OMP   DO SCHEDULE(STATIC)

DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)

  tab_nb = get_tab_nb(tab_ind(:,iG),tab_ba)
  tab_nq = get_tab_nq(tab_ind(:,iG),tab_ba)

  nnb = product(tab_nb)
  nnq = 1
  nb2 = 1
  nq2 = 1

  RTempG = reshape(SRep%SmolyakRep(iG)%R,shape=(/ nnq,nq2,nnb /))

  DO i=D,1,-1
    nb2 = tab_nb(i)
    nq2 = tab_nq(i)
    nnb = nnb / nb2

    RTempB = reshape(RTempG,shape=(/ nnq,nb2,nnb /))

    deallocate(RTempG)
    allocate(RTempG(nnq,nq2,nnb))

    DO ib=1,nnb
    DO iq=1,nnq
       RTempG(iq,:,ib) = sqrt(tab_ba(tab_ind(i,iG),i)%w(1:nq2)) * &
             matmul( tab_ba(tab_ind(i,iG),i)%d0b(1:nq2,1:nb2) , RTempB(iq,:,ib) )

       !$OMP ATOMIC
       nb_mult_BTOG = nb_mult_BTOG + nb2*nq2

    END DO
    END DO

    nnq = nnq * nq2
    deallocate(RTempB)

  END DO

  SRep%SmolyakRep(iG)%R = reshape(RTempG, shape=(/ nnq*nnb /) )
  !SRep%SmolyakRep(iG)%R = reshape(RTempG, (/ tab_nq(1),tab_nq(2) /) )

  deallocate(RTempG)
END DO
!$OMP   END DO

IF (allocated(tab_nb)) deallocate(tab_nb)
IF (allocated(tab_nq)) deallocate(tab_nq)


!$OMP   END PARALLEL

SRep%Grid = .TRUE.

END SUBROUTINE BSmolyakRep_TO_GSmolyakRep_withsqrtW
FUNCTION Set_weight_TO_SmolyakRep(tab_ind,tab_ba) RESULT (SRep)
USE mod_Smolyak_ba, only: typeba, get_tab_nq, get_tab_nb
IMPLICIT NONE

TYPE(Type_SmolyakRep)                              :: SRep
integer,         allocatable,    intent(in)        :: tab_ind(:,:) ! tab_ind(D,MaxnD)
TYPE(TypeBa),    allocatable,    intent(in)        :: tab_ba(:,:) ! tab_ba(L,D)

integer               :: i,D,iG,nb_BG

integer                            :: nnq,nnq1,nq2,nnq3,iq1,iq3
integer, allocatable               :: tab_nq(:)
real(kind=Rkind), allocatable      :: RTempG(:,:,:)


CALL alloc_SmolyakRep(SRep,tab_ind,tab_ba,grid=.TRUE.)
!CALL Write_SmolyakRep(Srep)

D = size(tab_ind(:,1))

DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)

  tab_nq = get_tab_nq(tab_ind(:,iG),tab_ba)

  nnq = product(tab_nq)


  nnq1 = nnq
  nq2  = 1
  nnq3 = 1

  SRep%SmolyakRep(iG)%R = ONE
  RTempG = reshape(SRep%SmolyakRep(iG)%R,shape=(/ nnq,1,1 /))

  DO i=1,D
    nq2  = tab_nq(i)
    nnq1 = nnq1/nq2

    RTempG = reshape(RTempG,shape=(/ nnq1,nq2,nnq3 /))

    DO iq3=1,nnq3
    DO iq1=1,nnq1
       RTempG(iq1,:,iq3) = tab_ba(tab_ind(i,iG),i)%w(1:nq2) * RTempG(iq1,1:nq2,iq3)
    END DO
    END DO

    nnq3 = nnq3 * nq2


  END DO

  SRep%SmolyakRep(iG)%R = reshape(RTempG, shape=(/ nnq /) )
  !SRep%SmolyakRep(iG)%R = reshape(RTempG, (/ tab_nq(1),tab_nq(2) /) )

  deallocate(RTempG)

END DO

IF (allocated(tab_nq)) deallocate(tab_nq)

!CALL Write_SmolyakRep(Srep)

END FUNCTION Set_weight_TO_SmolyakRep
FUNCTION Set_V_TO_SmolyakRep(tab_ind,tab_ba) RESULT (SRep)
USE mod_Smolyak_DInd, only : InD_TO_tabi
USE mod_Smolyak_ba, only: typeba, get_tab_nq, get_tab_nb
IMPLICIT NONE

TYPE(Type_SmolyakRep)                              :: SRep         ! potential
integer,         allocatable,    intent(in)        :: tab_ind(:,:) ! tab_ind(D,MaxnD)
TYPE(TypeBa),    allocatable,    intent(in)        :: tab_ba(:,:)  ! tab_ba(L,D)

integer               :: i,D,iG,nb_BG

integer                            :: iq,nnq
integer, allocatable               :: tab_nq(:),tab_i(:)
real (kind=Rkind), allocatable     :: tab_x(:)


CALL alloc_SmolyakRep(SRep,tab_ind,tab_ba,grid=.TRUE.)
!CALL Write_SmolyakRep(Srep)

D = size(tab_ind(:,1))
allocate(tab_x(D))
allocate(tab_i(D))

DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)

  tab_nq = get_tab_nq(tab_ind(:,iG),tab_ba)

  !write(out_unitp,*) 'iG,tabnq',iG,':',tab_nq
  DO iq=1,product(tab_nq)
    CALL InD_TO_tabi(iq,D,tab_nq,tab_i)
    DO i=1,D
      tab_x(i) = tab_ba(tab_ind(i,iG),i)%x(tab_i(i))
    END DO
    !write(out_unitp,*) 'iq,tab_i',iq,':',tab_i,'tab_nq',tab_nq,'x:',tab_x

    SRep%SmolyakRep(iG)%R(iq) = HALF * dot_product(tab_x,tab_x)
    !SRep%SmolyakRep(iG)%R(tab_i(1),tab_i(2)) = HALF * dot_product(tab_x,tab_x)

    !write(out_unitp,*) 'iG,iq,tab_i',iG,iq,':',tab_i,'x:',tab_x,'V:',SRep%SmolyakRep(iG)%R(iq)


  END DO
END DO

IF (allocated(tab_nq)) deallocate(tab_nq)
IF (allocated(tab_i))  deallocate(tab_i)
IF (allocated(tab_x))  deallocate(tab_x)

!CALL Write_SmolyakRep(Srep)

END FUNCTION Set_V_TO_SmolyakRep
SUBROUTINE Norm_OFF_Diff_TabRDP(TabRDP1,TabRDP2)
IMPLICIT NONE

TYPE(TypeRDP)      :: TabRDP1(:),TabRDP2(:)


integer            :: iG
real(kind=Rkind)   :: Norm


IF (TabRDP1(1)%n1 /= TabRDP2(1)%n1 .OR.                                 &
    TabRDP1(1)%n2 /= TabRDP2(1)%n2 .OR.                                 &
    TabRDP1(1)%n3 /= TabRDP2(1)%n3 .OR.                                 &
    size(TabRDP1) /= size(TabRDP2) ) THEN
  write(out_unitp,*) ' ERROR in Norm_OFF_Diff_TabRDP'
  write(out_unitp,*) ' incompatible TabRDP1 and TabRDP2'
  STOP
END IF

Norm = ZERO
DO iG=1,ubound(TabRDP1,dim=1)
  Norm = Norm + sum( (TabRDP1(iG)%RDP(:,:,:) - TabRDP2(iG)%RDP(:,:,:))**2 )
END DO
write(out_unitp,*) 'Norm_OFF_Diff_TabRDP',Norm


END SUBROUTINE Norm_OFF_Diff_TabRDP
SUBROUTINE TabRDP2_TO_TabRDP1(TabRDP1,TabRDP2)
IMPLICIT NONE

TYPE(TypeRDP), allocatable :: TabRDP1(:),TabRDP2(:)


integer            :: iG
real(kind=Rkind)   :: Norm

IF (allocated(TabRDP1)) CALL dealloc_TabRDP(TabRDP1)

IF (allocated(TabRDP2)) THEN
  allocate(TabRDP1( lbound(TabRDP2,dim=1) : ubound(TabRDP2,dim=1) ))

  DO iG=lbound(TabRDP1,dim=1),ubound(TabRDP1,dim=1)

    TabRDP1(iG)%n1 = TabRDP2(iG)%n1
    TabRDP1(iG)%n2 = TabRDP2(iG)%n2
    TabRDP1(iG)%n3 = TabRDP2(iG)%n3

    allocate(TabRDP1(iG)%RDP(TabRDP1(iG)%n1,TabRDP1(iG)%n2,TabRDP1(iG)%n3))
    TabRDP1(iG)%RDP(:,:,:) = TabRDP2(iG)%RDP(:,:,:)

  END DO
END IF

END SUBROUTINE TabRDP2_TO_TabRDP1


END MODULE mod_Smolyak_RDP
