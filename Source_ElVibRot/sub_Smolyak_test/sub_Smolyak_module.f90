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
MODULE mod_Smolyak_test
USE mod_Smolyak_DInd
USE mod_Smolyak_RDP
USE mod_Smolyak_ba
USE mod_system
IMPLICIT NONE

logical :: nested = .FALSE.

TYPE Type_AllParam

  integer :: D  = 0
  integer :: LB = 0
  integer :: LG = 0

  TYPE(TypeDInd),  allocatable :: ind_Grid(:)
  TYPE(TypeDInd),  allocatable :: ind_Basis(:)
  TYPE(TypeBa),    allocatable :: tab_ba(:,:) ! tab_ba(L,D)

  real(kind=Rkind), allocatable :: WSG(:) ! weight of the Smolyak grids

END TYPE Type_AllParam

CONTAINS


SUBROUTINE Set_SmolyakWeight(WSG,indGrid,D,LG)
USE mod_system
IMPLICIT NONE

integer          :: D,LG
TYPE(TypeDInd)   :: indGrid
real(kind=Rkind), allocatable :: WSG(:)

real(kind=Rkind) :: binomial ! function
integer          :: iG,L

IF (allocated(WSG)) deallocate(WSG)

allocate(WSG(indGrid%MaxnD))

!write(out_unitp,*) 'D,LG',D,LG
!CALL Write_TypeDInd(indGrid)


!weight of the Smolyak grids
IF (nested) THEN
  WSG(:) = ONE
ELSE
  DO iG=1,indGrid%MaxnD
    L = sum(indGrid%tab_ind(:,iG))
    !write(out_unitp,*) 'iG,L',iG,L,'ind:',indGrid%tab_ind(:,iG)
    !write(out_unitp,*) 'LG-L,D-1',LG-L,D-1
    IF (LG-L > D-1) THEN
      WSG(iG) =  ZERO
    ELSE IF (mod(LG-L,2) == 0) THEN
      WSG(iG) =  binomial(D-1,LG-L)
    ELSE
      WSG(iG) = -binomial(D-1,LG-L)
    END IF
    !write(out_unitp,*) 'iG,l(:)',iG,indGrid%tab_ind(:,iG),WSG(iG)
  END DO
END IF

END SUBROUTINE Set_SmolyakWeight

SUBROUTINE Size_OF_WP0_ON_GRID(nq_WP,indGrid,tab_ba,D,LG)
USE mod_system
IMPLICIT NONE


integer          :: nq_WP,D,LG


TYPE(TypeDInd) :: indGrid
TYPE(TypeBa), allocatable   :: tab_ba(:,:)

integer          :: iG,nqq,i,li

!-------------------------------------------
!-------------------------------------------
! Size of the WP
nq_WP = 0

DO iG=1,indGrid%MaxnD
  nqq = 1
  DO i=1,D
    li = indGrid%tab_ind(i,iG)
    !write(out_unitp,*) 'li,i',li,i ; flush(out_unitp)
    nqq = nqq * tab_ba(li,i)%nq
  END DO
  nq_WP = nq_WP + nqq
END DO

write(out_unitp,*) 'nq_WP',nq_WP
!-------------------------------------------
!-------------------------------------------
END SUBROUTINE Size_OF_WP0_ON_GRID

SUBROUTINE WP0_ON_GRID(RWPG,tablb0,WSG,indGrid,tab_ba,D,LG)
USE mod_system
IMPLICIT NONE


integer          :: nq_WP,D,LG
TYPE(TypeDInd)   :: indGrid
TYPE(TypeBa), allocatable   :: tab_ba(:,:)

real(kind=Rkind), allocatable :: RWPG(:)
integer          :: tablb0(:)

real(kind=Rkind), allocatable :: WSG(:)


integer               :: i_WP,iG,i,nqq,iqq
real(kind=Rkind)      :: x,w,g,Norm
integer, allocatable  :: tabiq(:),tabnq(:),tabl(:)


real(kind=Rkind) :: poly_Hermite_exp ! function


IF (allocated(RWPG)) deallocate(RWPG)


CALL Size_OF_WP0_ON_GRID(nq_WP,indGrid,tab_ba,D,LG)

allocate(RWPG(nq_WP))
!-------------------------------------------
!-------------------------------------------
! Allocation of WP on the grid (old way)
allocate(tabl(D))
allocate(tabnq(D))
allocate(tabiq(D))

! then the WP
i_WP = 0
DO iG=1,indGrid%MaxnD
  DO i=1,D
    tabl(i)  = indGrid%tab_ind(i,iG)
    tabnq(i) = tab_ba(tabl(i),i)%nq
  END DO
  nqq = product(tabnq)

  !write(out_unitp,*) 'iG,tabnq',iG,':',tabnq
  DO iqq=1,nqq
    i_WP = i_WP + 1
    CALL InD_TO_tabi(iqq,D,tabnq,tabiq)
    !write(out_unitp,*) 'iqq,tabiq',iqq,':',tabiq
    g = ONE
    DO i=1,D
      x = tab_ba(tabl(i),i)%x(tabiq(i))
      g = g * poly_Hermite_exp(x,tablb0(i))
    END DO
    RWPG(i_WP) = g
  END DO
END DO



i_WP = 0
Norm = ZERO
DO iG=1,indGrid%MaxnD
  DO i=1,D
    tabl(i)  = indGrid%tab_ind(i,iG)
    tabnq(i) = tab_ba(tabl(i),i)%nq
  END DO
  nqq = product(tabnq)

  DO iqq=1,nqq
    i_WP = i_WP + 1

    ! the weight with Smolyak coeficient
    CALL InD_TO_tabi(iqq,D,tabnq,tabiq)
    w = WSG(iG)
    DO i=1,D
      w = w * tab_ba(tabl(i),i)%w(tabiq(i))
    END DO
    Norm = Norm + RWPG(i_WP)**2 * w
  END DO
END DO

write(out_unitp,*) 'tablb0',tablb0

write(out_unitp,*) 'Norm of WP0',Norm

deallocate(tabl)
deallocate(tabnq)
deallocate(tabiq)

!-------------------------------------------
!-------------------------------------------
END SUBROUTINE WP0_ON_GRID

SUBROUTINE Set_BgG_FOR_id(BgG,ind_Grid,ind_Basis,tab_ba,D,LG,id)
USE mod_system
IMPLICIT NONE

TYPE(TypeRDP), allocatable    :: BgG(:)

integer                       :: D,LG,id
TYPE(TypeDInd), allocatable   :: ind_Grid(:)
TYPE(TypeDInd), allocatable   :: ind_Basis(:)
TYPE(TypeBa),   allocatable   :: tab_ba(:,:)

integer          :: iG,i,li,nb,nbb,nq,nqq,iq,ll
logical, parameter :: debug = .FALSE.
!logical, parameter :: debug = .TRUE.

CALL dealloc_TabRDP(BgG)

IF (debug) THEN
  write(out_unitp,*) 'BEGINNING in Set_BgG_FOR_id',id
  write(out_unitp,*) 'tab_q, for B',ind_Basis(id+1)%tab_q
  write(out_unitp,*) 'tab_q, for G',ind_Grid(id)%tab_q
  write(out_unitp,*) 'tab_q, for g',0
END IF

allocate(BgG(ind_Grid(id)%MaxnD))
DO iG=1,ind_Grid(id)%MaxnD

  nqq = 1
  !write(out_unitp,*) 'iG,l(:)',ind_Grid(id)%tab_ind(:,iG)
  DO i=1,ind_Grid(id)%ndim
    li = ind_Grid(id)%tab_ind(i,iG)
    iq = ind_Grid(id)%tab_q(i)
    IF (iq < 1 .OR. iq > D) CYCLE
    nqq = nqq * tab_ba(li,iq)%nq
  END DO

  nq = 1

  nbb = ind_Basis(id+1)%MaxnD


  BgG(iG)%n1 = nq
  BgG(iG)%n2 = nqq
  BgG(iG)%n3 = nbb

  allocate(BgG(iG)%RDP(nq,nqq,nbb))
  BgG(iG)%RDP(:,:,:) = ZERO

  IF (debug) write(out_unitp,*) 'BgG, at id and iG',id,iG,shape(BgG(iG)%RDP)

END DO

IF (debug) THEN
  CALL Write_TabRDP(BgG)
  write(out_unitp,*) 'END in Set_BgG_FOR_id',id
END IF

END SUBROUTINE Set_BgG_FOR_id

SUBROUTINE Set_BbG_FOR_id(BbG,ind_Grid,ind_Basis,tab_ba,D,LG,id)
USE mod_system
IMPLICIT NONE

TYPE(TypeRDP), allocatable    :: BbG(:)

integer          :: D,LG,id
TYPE(TypeDInd), allocatable   :: ind_Grid(:)
TYPE(TypeDInd), allocatable   :: ind_Basis(:)
TYPE(TypeBa), allocatable     :: tab_ba(:,:)



integer               :: iG,i,li,nb,nbb,nq,nqq,iq
integer               :: idba,idg
logical, parameter    :: debug = .FALSE.
!logical, parameter    :: debug = .TRUE.

CALL dealloc_TabRDP(BbG)

IF (debug) THEN
  write(out_unitp,*) 'BEGINNING in Set_BbG_FOR_id',id
  write(out_unitp,*)
  write(out_unitp,*) 'tab_q, for B',ind_Basis(id)%tab_q
  write(out_unitp,*) 'tab_q, for G',ind_Grid(id)%tab_q
  iq = id
  IF (iq>=1 .AND. iq <=D) THEN
    write(out_unitp,*) 'tab_q, for b',id
  ELSE
    write(out_unitp,*) 'tab_q, for b',0
  END IF
END IF

allocate(BbG(ind_Grid(id)%MaxnD)) ! smaller number of Somlyak grid

DO iG=1,ind_Grid(id)%MaxnD
  nqq = 1
  DO i=1,ind_Grid(id)%ndim
    li = ind_Grid(id)%tab_ind(i,iG)
    iq = ind_Grid(id)%tab_q(i)
    IF (iq < 1 .OR. iq > D) CYCLE
    nqq = nqq * tab_ba(li,iq)%nq
  END DO


  iq = id
  IF (iq>=1 .AND. iq <=D) THEN
    nb = tab_ba(LG,id)%nb  ! always the largest value ??
  ELSE
    nb = 1
  END IF

  nbb = ind_Basis(id)%MaxnD

  BbG(iG)%n1 = nb
  BbG(iG)%n2 = nqq
  BbG(iG)%n3 = nbb

  allocate(BbG(iG)%RDP(nb,nqq,nbb))
  BbG(iG)%RDP(:,:,:) = ZERO

  IF (debug) write(out_unitp,*) 'BbG, at id and iG',id,iG,':',shape(BbG(iG)%RDP)

END DO

IF (debug) THEN
  CALL Write_TabRDP(BbG)
  write(out_unitp,*) 'END in Set_BbG_FOR_id',id
END IF

END SUBROUTINE Set_BbG_FOR_id
SUBROUTINE Transfer_WP0_TO_BgG(WPG,BgG)
USE mod_system
IMPLICIT NONE

TYPE(TypeRDP)      :: BgG(:)
real(kind=Rkind)   :: WPG(:)


integer               :: i_WP,iG,nqq


IF (BgG(1)%n1 /= 1 .OR. BgG(1)%n3 /= 1) THEN
  write(out_unitp,*) ' ERROR in Transfer_WP0_TO_BgG'
  write(out_unitp,*) ' BgG is not completely on the G grid'
  STOP
END IF
IF (size(WPG) /= sum(BgG(:)%n2)) THEN
  write(out_unitp,*) ' ERROR in Transfer_WP0_TO_BgG'
  write(out_unitp,*) ' The sizes of WPG and BgG are not compatible'
  write(out_unitp,*) 'size(WPG)',size(WPG)
  write(out_unitp,*) 'sum(BgG(:)%n2)',sum(BgG(:)%n2)
  STOP
END IF
write(out_unitp,*) 'sum(sq...)',sum(WPG**2)


i_WP = 0
DO iG=1,ubound(BgG,dim=1)
  nqq = BgG(iG)%n2
  BgG(iG)%RDP(1,:,1) = WPG(i_WP+1:i_WP+nqq)
  i_WP = i_WP+nqq
END DO
!write(out_unitp,*) 'i_WP',i_WP


END SUBROUTINE Transfer_WP0_TO_BgG

SUBROUTINE Norm_OFF_Diff_WP0_BgG(WPG,BgG)
USE mod_system
IMPLICIT NONE

TYPE(TypeRDP)      :: BgG(:)
real(kind=Rkind)   :: WPG(:)


integer            :: i_WP,iG,nqq
real(kind=Rkind)   :: Norm


IF (BgG(1)%n1 /= 1 .OR. BgG(1)%n3 /= 1) THEN
  write(out_unitp,*) ' ERROR in Norm_OFF_Diff_WP0_BgG'
  write(out_unitp,*) ' BgG is not completely on the G grid'
  STOP
END IF
IF (size(WPG) /= sum(BgG(:)%n2)) THEN
  write(out_unitp,*) ' ERROR in Norm_OFF_Diff_WP0_BgG'
  write(out_unitp,*) ' The sizes of WPG and BgG are not compatible'
  write(out_unitp,*) 'size(WPG)',size(WPG)
  write(out_unitp,*) 'sum(BgG(:)%n2)',sum(BgG(:)%n2)
  STOP
END IF
!write(out_unitp,*) 'sum(sq...)',sum(WPG**2)


i_WP = 0
Norm = ZERO
DO iG=1,ubound(BgG,dim=1)
  nqq = BgG(iG)%n2
  Norm = Norm + sum( (BgG(iG)%RDP(1,:,1) - WPG(i_WP+1:i_WP+nqq))**2 )
  i_WP = i_WP+nqq
END DO
write(out_unitp,*) 'Norm_OFF_Diff_WP0_BgG',Norm


END SUBROUTINE Norm_OFF_Diff_WP0_BgG

SUBROUTINE Compare_WP0_BbG_ON_basis(tablb0,BgG,AllPara)
USE mod_system
IMPLICIT NONE

TYPE(TypeRDP)        :: BgG(:)
TYPE (Type_AllParam) :: AllPara
integer              :: tablb0(:)


integer            :: ibb,iG,nqq,nb_not_zero
real(kind=Rkind)   :: Norm,coef


IF (BgG(1)%n1 /= 1 .OR. BgG(1)%n2 /= 1 .OR. size(BgG) /= 1) THEN
  write(out_unitp,*) ' ERROR in Compare_WP0_BbG_ON_basis'
  write(out_unitp,*) ' BgG is not completely on the basis'
  STOP
END IF


write(out_unitp,*) 'shape BgG(1)%RDP',shape(BgG(1)%RDP)
Norm        = ZERO
nb_not_zero = 0

DO ibb=1,BgG(1)%n3 ! nbb (the others are 1)
  coef = BgG(1)%RDP(1,1,ibb)
  Norm = Norm + coef**2
  IF (abs(coef) > 1.d-10) THEN
    nb_not_zero = nb_not_zero + 1
    write(out_unitp,*) 'WP',AllPara%ind_Basis(AllPara%D+1)%tab_ind(:,ibb),coef
    IF (sum(abs(tablb0-AllPara%ind_Basis(AllPara%D+1)%tab_ind(:,ibb))) == 0) THEN
      write(out_unitp,*) 'WP0 OK',ibb,'/',shape(BgG(1)%RDP)
    ELSE
      write(out_unitp,*) 'WP0 NOT OK',ibb,'/',shape(BgG(1)%RDP)
    END IF
  END IF
END DO
IF (nb_not_zero == 0) THEN
  write(out_unitp,*) 'WP0 NOT OK (0)',ibb,'/',shape(BgG(1)%RDP)
END IF
write(out_unitp,*) 'Norm',Norm



END SUBROUTINE Compare_WP0_BbG_ON_basis

SUBROUTINE Analysis_BbG_ON_basis(BgG,AllPara)
USE mod_system
IMPLICIT NONE

TYPE(TypeRDP)        :: BgG(:)
TYPE (Type_AllParam) :: AllPara


integer            :: ibb,iG,i,ib(AllPara%D)
real(kind=Rkind)   :: Norm,coef,harmo


IF (BgG(1)%n1 /= 1 .OR. BgG(1)%n2 /= 1 .OR. size(BgG) /= 1) THEN
  write(out_unitp,*) ' ERROR in Compare_WP0_BbG_ON_basis'
  write(out_unitp,*) ' BgG is not completely on the basis'
  STOP
END IF


write(out_unitp,*) 'shape BgG(1)%RDP',shape(BgG(1)%RDP)

Norm        = ZERO
DO ibb=1,BgG(1)%n3 ! nbb (the others are 1)
  coef = BgG(1)%RDP(1,1,ibb)
  Norm = Norm + coef**2
  IF (abs(coef) > ONETENTH**8) THEN
    ib(:) = AllPara%ind_Basis(AllPara%D+1)%tab_ind(:,ibb)
    harmo = ZERO
    DO i=1,AllPara%D
      harmo = harmo + HALF + ONE*real(ib(i)-1,kind=Rkind)
    END DO
    write(out_unitp,*) 'WP, ib(:)-1,coef,diff: ',ib(:)-1,coef,coef-harmo
    IF (abs(coef-harmo) > ONETENTH**8) write(out_unitp,*) 'WARNING large difference'
  END IF
END DO
write(out_unitp,*) 'Norm',Norm



END SUBROUTINE Analysis_BbG_ON_basis

SUBROUTINE BgG_TO_BbG(BgG,BbG,WSG,ind_Grid,ind_Basis,tab_ba,D,LG,LB,id,nb_mult)

USE mod_system
IMPLICIT NONE

TYPE(TypeRDP), allocatable    :: BgG(:)
TYPE(TypeRDP), allocatable    :: BbG(:)
real(kind=Rkind)   :: WSG(:)

integer            :: D,LG,LB,id
TYPE(TypeDInd), allocatable   :: ind_Grid(:)
TYPE(TypeDInd), allocatable   :: ind_Basis(:)
TYPE(TypeBa),   allocatable   :: tab_ba(:,:)



integer               :: ibbqq,i_WP,iG,iq,iqq,nq,nqq,nGp1,ibb,nbb,iqqi,iqqf,ib,nb,li,i,iGm1,lbb,lqq
integer               :: nb_mult
integer               :: tabl(D)
integer               :: tabnq(D)
real(kind=Rkind), allocatable  :: gwc(:)
real(kind=Rkind), allocatable  :: b(:)
logical, parameter :: debug = .FALSE.
!logical, parameter :: debug = .TRUE.

IF (debug) THEN
  write(out_unitp,*)
  write(out_unitp,*) 'BEGINNING BgG_TO_BbG',id
END IF

CALL Set_BbG_FOR_id(BbG,ind_Grid,ind_Basis,tab_ba,D,LG,id)

IF (debug) THEN
  write(out_unitp,*) 'size(BgG),size(BbG)',size(BgG),size(BbG)
  flush(out_unitp)
END IF

! id index of "b" of BbG
IF (size(BgG) /= ind_Grid(id-1)%MaxnD) THEN
  write(out_unitp,*) ' ERROR in BgG_TO_BbG'
  write(out_unitp,*) ' size(BgG) /= ind_Grid(id-1)%MaxnD',size(BgG),ind_Grid(id-1)%MaxnD
  STOP
END IF

IF (size(BbG) /= ind_Grid(id)%MaxnD) THEN
  write(out_unitp,*) ' ERROR in BgG_TO_BbG'
  write(out_unitp,*) ' size(BbG) /= ind_Grid(id)%MaxnD',size(BbG),ind_Grid(id)%MaxnD
  STOP
END IF

nGp1 = maxval(ind_Grid(id-1)%indD_OF_Dm1)
IF (nGp1 /= size(BbG)) THEN
  write(out_unitp,*) ' ERROR in BgG_TO_BbG'
  write(out_unitp,*) ' size(BbG) /= nGp1',size(BbG),nGp1
  STOP
END IF


IF (BgG(1)%n3 /= BbG(1)%n3) THEN
  write(out_unitp,*) ' ERROR in BgG_TO_BbG'
  write(out_unitp,*) ' BgG%n3 /= BbG%n3',BgG%n3,BbG%n3
  STOP
END IF

!write(out_unitp,*) 'max_val nqi',maxval(tab_ba(:,:)%nq)
allocate(gwc(maxval(tab_ba(:,:)%nq)))
allocate(b(maxval(tab_ba(:,:)%nb)))

nb_mult = 0
! !$ write(out_unitp,*) 'nb of threads',BasisTOGrid_maxth

DO iG=1,size(BgG)

  DO i=1,ind_Grid(id-1)%ndim
    tabl(i) = ind_Grid(id-1)%tab_ind(i,iG)
    iq = ind_Grid(id-1)%tab_q(i)
    IF (iq < 1 .OR. iq > D) CYCLE
    tabnq(i) = tab_ba(tabl(i),iq)%nq
  END DO
  !write(out_unitp,*) 'iG,tabl ',iG,tabl(1:ind_Grid(id-1)%ndim)
  !write(out_unitp,*) 'iG,tabnq',iG,tabnq(1:ind_Grid(id-1)%ndim)

  nbb  = BgG(iG)%n3
  nqq  = product(tabnq(2:ind_Grid(id-1)%ndim))
  iGm1 = ind_Grid(id-1)%indD_OF_Dm1(iG)

  !write(out_unitp,*) '---------------------' ; flush(out_unitp)
  !write(out_unitp,*) 'iG,nbb,nqq',iG,nbb,nqq ; flush(out_unitp)

  !$OMP   PARALLEL &
  !$OMP   DEFAULT(NONE) &
  !$OMP   SHARED(LB,id,iG,nbb,nqq,ind_Basis,ind_Grid,tab_ba,tabl,tabnq,BbG,BgG,WSG,nb_mult) &
  !$OMP   SHARED(iGm1) &
  !$OMP   PRIVATE(ibbqq,ibb,lbb,iqq,li,iq,ib,nb,nq,gwc,b) &
  !$OMP   NUM_THREADS(BasisTOGrid_maxth)

  allocate(gwc(maxval(tabnq)))
  allocate(b(maxval(tab_ba(:,:)%nb)))

  !$OMP   DO SCHEDULE(dynamic)
  DO ibbqq=1,nbb*nqq
    ibb = (ibbqq-1)/nqq + 1
    iqq = mod((ibbqq-1),nqq)+1
    !write(out_unitp,*) 'nbb,nqq',nbb,nqq
    !write(out_unitp,*) 'ibbqq,ibb,iqq',ibbqq,ibb,iqq

    lbb = ind_Basis(id)%i_TO_l(ibb)

      li = min(tabl(1),LB-lbb)

      iq = ind_Grid(id-1)%tab_q(1)
      nb = tab_ba(li,iq)%nb

      li = tabl(1)
      nq = tabnq(1)

      gwc(1:nq) = BgG(iG)%RDP(1,(iqq-1)*nq+1:iqq*nq,ibb)*tab_ba(li,iq)%w(:)
      IF (id == 1) THEN
        b(1:nb)   = matmul(gwc(1:nq),tab_ba(li,iq)%d0b(:,1:nb)) * WSG(iG)
      ELSE
        b(1:nb)   = matmul(gwc(1:nq),tab_ba(li,iq)%d0b(:,1:nb))
      END IF

!     !$OMP CRITICAL (BgG_TO_BbG_CRIT)
!     BbG(iGm1)%RDP(1:nb,iqq,ibb) = BbG(iGm1)%RDP(1:nb,iqq,ibb) + b(1:nb)
!     !$OMP END CRITICAL (BgG_TO_BbG_CRIT)

      DO ib=1,nb
        !$OMP ATOMIC
        BbG(iGm1)%RDP(ib,iqq,ibb) = BbG(iGm1)%RDP(ib,iqq,ibb) + b(ib)
      END DO

      !$OMP ATOMIC
      nb_mult = nb_mult + nq*nb

  END DO
  !$OMP   END DO
  deallocate(gwc)
  deallocate(b)
  !$OMP   END PARALLEL

END DO


CALL dealloc_TabRDP(BgG)

IF (debug) THEN
  !CALL Write_TabRDP(BbG)
  !CALL SumSq_TabRDP(BbG)
  write(out_unitp,*) 'END BgG_TO_BbG',id
  flush(out_unitp)
END IF

END SUBROUTINE BgG_TO_BbG

SUBROUTINE BbG_TO_BgG(BbG,BgG,ind_Grid,ind_Basis,tab_ba,D,LG,LB,id,nb_mult)

USE mod_system
IMPLICIT NONE

TYPE(TypeRDP), allocatable    :: BgG(:)
TYPE(TypeRDP), allocatable    :: BbG(:)

integer            :: D,LG,LB,id
TYPE(TypeDInd), allocatable   :: ind_Grid(:)
TYPE(TypeDInd), allocatable   :: ind_Basis(:)
TYPE(TypeBa),   allocatable   :: tab_ba(:,:)



integer               :: i_WP,iG,iq,iqq,nq,nqq,ibb,nbb,iqqi,iqqf,ib,nb,li,i,iGm1,lbb,lqq
integer               :: nb_mult
integer               :: tabl(D)
integer               :: tabnq(D)
real(kind=Rkind), allocatable  :: b(:)
logical, parameter :: debug = .FALSE.
!logical, parameter :: debug = .TRUE.

IF (debug) THEN
  write(out_unitp,*)
  write(out_unitp,*) 'BEGINNING BbG_TO_BgG',id
END IF

CALL Set_BgG_FOR_id(BgG,ind_Grid,ind_Basis,tab_ba,D,LG,id-1)

IF (debug) THEN
  write(out_unitp,*) 'size(BgG),size(BbG)',size(BgG),size(BbG)
  write(out_unitp,*) 'BbG'
  CALL Write_TabRDP(BbG)
  flush(out_unitp)
END IF

! id index of "b" of BbG
IF (size(BgG) /= ind_Grid(id-1)%MaxnD) THEN
  write(out_unitp,*) ' ERROR in BbG_TO_BgG'
  write(out_unitp,*) ' size(BgG) /= ind_Grid(id-1)%MaxnD',size(BgG),ind_Grid(id-1)%MaxnD
  STOP
END IF

IF (size(BbG) /= ind_Grid(id)%MaxnD) THEN
  write(out_unitp,*) ' ERROR in BbG_TO_BgG'
  write(out_unitp,*) ' size(BbG) /= ind_Grid(id)%MaxnD',size(BbG),ind_Grid(id)%MaxnD
  STOP
END IF

iGm1 = maxval(ind_Grid(id-1)%indD_OF_Dm1)
IF (iGm1 /= size(BbG)) THEN
  write(out_unitp,*) ' ERROR in BbG_TO_BgG'
  write(out_unitp,*) ' size(BbG) /= iGm1',size(BbG),iGm1
  STOP
END IF

IF (BgG(1)%n3 /= BbG(1)%n3) THEN
  write(out_unitp,*) ' ERROR in BbG_TO_BgG'
  write(out_unitp,*) ' BgG%n3 /= BbG%n3',BgG(1)%n3,BbG(1)%n3
  STOP
END IF

!write(out_unitp,*) 'max_val nqi',maxval(tab_ba(:,:)%nq)
allocate(b(maxval(tab_ba(:,:)%nb)))

nb_mult = 0

DO iG=1,size(BgG)

  !write(out_unitp,*) 'B.G=>BbG'
  DO i=1,ind_Grid(id-1)%ndim
    tabl(i) = ind_Grid(id-1)%tab_ind(i,iG)
    iq = ind_Grid(id-1)%tab_q(i)
    IF (iq < 1 .OR. iq > D) CYCLE
    tabnq(i) = tab_ba(tabl(i),iq)%nq
  END DO
  !write(out_unitp,*) 'iG,tabl ',iG,tabl(1:ind_Grid(id-1)%ndim)
  !write(out_unitp,*) 'iG,tabnq',iG,tabnq(1:ind_Grid(id-1)%ndim)

  nbb  = BgG(iG)%n3
  nqq  = product(tabnq(2:ind_Grid(id-1)%ndim))
  iGm1 = ind_Grid(id-1)%indD_OF_Dm1(iG)

  !write(out_unitp,*) '---------------------' ; flush(out_unitp)
  !write(out_unitp,*) 'iG,nbb,nqq',iG,nbb,nqq ; flush(out_unitp)

  DO ibb=1,nbb

    iqqi = 0
    iqqf = 0
    DO iqq=1,nqq

      li = min(tabl(1),LB)

      iq = ind_Grid(id-1)%tab_q(1)
      nb = tab_ba(li,iq)%nb

      li = tabl(1)
      nq = tabnq(1)
      iqqf = iqqi + nq

      b(1:nb)   = BbG(iGm1)%RDP(1:nb,iqq,ibb)
      !write(out_unitp,*) 'iG,ibb,iqq',iG,ibb,iqq,'b',b(1:nb)

      BgG(iG)%RDP(1,iqqi+1:iqqf,ibb) = matmul(tab_ba(li,iq)%d0b(:,1:nb),b(1:nb))

      !write(out_unitp,*) 'iG,ibb,iqq',iG,ibb,iqq,'g',BgG(iG)%RDP(1,iqqi+1:iqqf,ibb)

      nb_mult = nb_mult + nq*nb

      iqqi = iqqf

    END DO

  END DO

END DO

deallocate(b)

CALL dealloc_TabRDP(BbG)

IF (debug) THEN
  write(out_unitp,*) 'BgG'
  CALL Write_TabRDP(BgG)
  CALL SumSq_TabRDP(BgG)
  write(out_unitp,*) 'END BbG_TO_BgG',id
  flush(out_unitp)
END IF

END SUBROUTINE BbG_TO_BgG
SUBROUTINE Transfer_BbG_TO_BgG(BbG,BgG,ind_Grid,ind_Basis,tab_ba,D,LG,LB,id)

USE mod_system
IMPLICIT NONE

TYPE(TypeRDP), allocatable    :: BgG(:)
TYPE(TypeRDP), allocatable    :: BbG(:)

integer            :: D,LG,LB,id
TYPE(TypeDInd), allocatable   :: ind_Grid(:)
TYPE(TypeDInd), allocatable   :: ind_Basis(:)
TYPE(TypeBa),   allocatable   :: tab_ba(:,:)



integer            :: iG,iqq,nq,nqq,ibb,nbb,ib,nb,lbb,l,ibbNew
logical, parameter :: debug = .FALSE.


IF (debug) THEN
  write(out_unitp,*)
  write(out_unitp,*) 'BEGINNING Transfer_BbG_TO_BgG',id
  write(out_unitp,*) 'LG,LB',LG,LB
END IF

CALL Set_BgG_FOR_id(BgG,ind_Grid,ind_Basis,tab_ba,D,LG,id)

IF (debug) THEN
  write(out_unitp,*) 'size(BgG),size(BbG)',size(BgG),size(BbG)
  flush(out_unitp)
END IF


! id index of "b" of BbG
IF (size(BgG) /= size(BbG)) THEN
  write(out_unitp,*) ' ERROR in Transfer_BbG_TO_BgG'
  write(out_unitp,*) ' size(BgG) /= size(BbG)',size(BgG),size(BbG)
  STOP
END IF


IF (BgG(1)%n2 /= BbG(1)%n2) THEN ! nqq
  write(out_unitp,*) ' ERROR in Transfer_BbG_TO_BgG'
  write(out_unitp,*) ' BgG(1)%n2 /= BbG(1)%n2',BgG(1)%n2,BbG(1)%n2
  STOP
END IF

IF (BgG(1)%n1 /= 1) THEN ! nqq
  write(out_unitp,*) ' ERROR in Transfer_BbG_TO_BgG'
  write(out_unitp,*) ' BgG(1)%n1 /= 1',BgG(1)%n1
  STOP
END IF

DO iG=1,size(BgG)

  nbb  = BbG(iG)%n3
  nqq  = BbG(iG)%n2
  !write(out_unitp,*) 'iG,nbb,nqq, shape tab_ind',iG,nbb,nqq,shape(ind_Basis(id)%tab_ind(:,:))
  DO iqq=1,nqq
    ibbNew = 0
    DO ibb=1,nbb
      lbb = ind_Basis(id)%i_TO_l(ibb)
      !write(out_unitp,*) 'iG,iqq,ibb,lbb,LB-lbb',iG,iqq,ibb,':',lbb,LB-lbb
      nb = tab_ba(LB-lbb,id)%nb
      BgG(iG)%RDP(1,iqq,ibbNew+1:ibbNew+nb) = BbG(iG)%RDP(1:nb,iqq,ibb)
      ibbNew = ibbNew+nb
    END DO
  END DO
END DO

CALL dealloc_TabRDP(BbG)

IF (debug) THEN
  !CALL Write_TabRDP(BgG)
  !CALL SumSq_TabRDP(BgG)
  write(out_unitp,*) 'END Transfer_BbG_TO_BgG',id
  flush(out_unitp)
END IF


END SUBROUTINE Transfer_BbG_TO_BgG

SUBROUTINE Transfer_BgG_TO_BbG(BgG,BbG,ind_Grid,ind_Basis,tab_ba,D,LG,LB,id)

USE mod_system
IMPLICIT NONE

TYPE(TypeRDP), allocatable    :: BgG(:)
TYPE(TypeRDP), allocatable    :: BbG(:)

integer            :: D,LG,LB,id
TYPE(TypeDInd), allocatable   :: ind_Grid(:)
TYPE(TypeDInd), allocatable   :: ind_Basis(:)
TYPE(TypeBa),   allocatable   :: tab_ba(:,:)



integer               :: iG,iqq,nq,nqq,ibb,nbb,ib,nb,lbb,l,ibbNew
logical, parameter :: debug = .FALSE.
!logical, parameter :: debug = .TRUE.

IF (debug) THEN
  write(out_unitp,*)
  write(out_unitp,*) 'BEGINNING Transfer_BgG_TO_BbG',id
  write(out_unitp,*) 'LG,LB',LG,LB
END IF

CALL Set_BbG_FOR_id(BbG,ind_Grid,ind_Basis,tab_ba,D,LG,id)

IF (debug) THEN
  write(out_unitp,*) 'size(BgG),size(BbG)',size(BgG),size(BbG)
  flush(out_unitp)
END IF


! id index of "b" of BbG
IF (size(BgG) /= size(BbG)) THEN
  write(out_unitp,*) ' ERROR in Transfer_BgG_TO_BbG'
  write(out_unitp,*) ' size(BgG) /= size(BbG)',size(BgG),size(BbG)
  STOP
END IF


IF (BgG(1)%n2 /= BbG(1)%n2) THEN ! nqq
  write(out_unitp,*) ' ERROR in Transfer_BgG_TO_BbG'
  write(out_unitp,*) ' BgG(1)%n2 /= BbG(1)%n2',BgG(1)%n2,BbG(1)%n2
  STOP
END IF

IF (BgG(1)%n1 /= 1) THEN ! nqq
  write(out_unitp,*) ' ERROR in Transfer_BgG_TO_BbG'
  write(out_unitp,*) ' BgG(1)%n1 /= 1',BgG(1)%n1
  STOP
END IF

DO iG=1,size(BgG)
  !write(out_unitp,*) 'shape BgG  %RDP',shape(BgG(iG)%RDP)
  !write(out_unitp,*) 'shape BbG  %RDP',shape(BbG(iG)%RDP)

  nbb  = BbG(iG)%n3
  nqq  = BbG(iG)%n2
  !write(out_unitp,*) 'iG,nbb,nqq, shape tab_ind',iG,nbb,nqq,shape(ind_Basis(id)%tab_ind(:,:))
  DO iqq=1,nqq
    ibbNew = 0
    DO ibb=1,nbb
      lbb = ind_Basis(id)%i_TO_l(ibb)
      !write(out_unitp,*) 'iG,iqq,ibb,lbb,LB-lbb',iG,iqq,ibb,':',lbb,LB-lbb
      nb = tab_ba(LB-lbb,id)%nb
      BbG(iG)%RDP(1:nb,iqq,ibb) = BgG(iG)%RDP(1,iqq,ibbNew+1:ibbNew+nb)
      ibbNew = ibbNew+nb
    END DO
  END DO
END DO

CALL dealloc_TabRDP(BgG)

IF (debug) THEN
  CALL Write_TabRDP(BbG)
  CALL SumSq_TabRDP(BbG)
  write(out_unitp,*) 'END Transfer_BgG_TO_BbG',id
  flush(out_unitp)
END IF


END SUBROUTINE Transfer_BgG_TO_BbG

SUBROUTINE Norm_OF_BgG(BgG,WSG,indGrid,tab_ba,D,LG)
USE mod_system
IMPLICIT NONE


TYPE(TypeRDP), allocatable    :: BgG(:)
real(kind=Rkind), allocatable :: WSG(:)

integer          :: D,LG
TYPE(TypeDInd)   :: indGrid
TYPE(TypeBa),     allocatable :: tab_ba(:,:)



integer               :: iG,i,nqq,iqq
real(kind=Rkind)      :: w,Norm
integer  :: tabiq(D),tabnq(D),tabl(D)

IF (BgG(1)%n1 /= 1 .OR. BgG(1)%n3 /= 1) THEN
  write(out_unitp,*) ' ERROR in Transfer_WP0_TO_BgG'
  write(out_unitp,*) ' BgG is not completely on the G grid'
  STOP
END IF

IF (size(BgG) /= size(WSG) ) THEN
  write(out_unitp,*) ' ERROR in Transfer_WP0_TO_BgG'
  write(out_unitp,*) ' The number of SG are different'
  write(out_unitp,*) ' size(BgG) /= size(WSG)',size(BgG),size(WSG)
  STOP
END IF

Norm = ZERO
DO iG=1,indGrid%MaxnD
  DO i=1,D
    tabl(i)  = indGrid%tab_ind(i,iG)
    tabnq(i) = tab_ba(tabl(i),i)%nq
  END DO
  nqq = product(tabnq)

  DO iqq=1,nqq

    ! the weight with Smolyak coeficient
    CALL InD_TO_tabi(iqq,D,tabnq,tabiq)
    w = WSG(iG)
    DO i=1,D
      w = w * tab_ba(tabl(i),i)%w(tabiq(i))
    END DO
    Norm = Norm + BgG(iG)%RDP(1,iqq,1)**2 * w
  END DO
END DO

write(out_unitp,*) 'Norm of BgB',Norm

!-------------------------------------------
!-------------------------------------------
END SUBROUTINE Norm_OF_BgG

SUBROUTINE sub_G_TO_B(WPG,WPB,AllPara)
USE mod_system
IMPLICIT NONE

TYPE(TypeRDP), allocatable    :: WPG(:)
TYPE(TypeRDP), allocatable    :: WPB(:)

TYPE (Type_AllParam) :: AllPara

integer :: id

TYPE(TypeRDP), allocatable    :: WPG_temp(:)

integer :: nb_BG,nb_mult_id

logical, parameter :: debug=.FALSE.

IF (debug) THEN
  write(out_unitp,*) '====================================='
  write(out_unitp,*) '====================================='
  write(out_unitp,*) '====================================='
  CALL time_perso('sub_G_TO_B')
  write(out_unitp,*) '====================================='
  flush(out_unitp)
END IF

CALL TabRDP2_TO_TabRDP1(WPG_temp,WPG)

DO id=1,AllPara%D
  IF (debug) THEN
    CALL Size_TabRDP(WPG_temp,nb_BG)
    write(out_unitp,*) 'id, size WPG',id,int(real(nb_BG,kind=Rkind)*EIGHT/(1024.d0**2)),' MB'
    flush(out_unitp)
  END IF

  CALL BgG_TO_BbG(WPG_temp,WPB,AllPara%WSG,                             &
                   AllPara%ind_Grid,AllPara%ind_Basis,AllPara%tab_ba,   &
                   AllPara%D,AllPara%LG,AllPara%LB,id,nb_mult_id)

  IF (debug) THEN
    CALL Size_TabRDP(WPB,nb_BG)
    write(out_unitp,*) 'id, size WPB',id,int(real(nb_BG,kind=Rkind)*EIGHT/(1024.d0**2)),' MB'
    write(out_unitp,*) 'id, nb_mult_id ',id,nb_mult_id
    flush(out_unitp)
  END IF

  CALL Transfer_BbG_TO_BgG(WPB,WPG_temp,                                &
                   AllPara%ind_Grid,AllPara%ind_Basis,AllPara%tab_ba,   &
                   AllPara%D,AllPara%LG,AllPara%LB,id)
  !write(out_unitp,*) 'sub_G_TO_B',id,' done' ; flush(out_unitp)
END DO

CALL TabRDP2_TO_TabRDP1(WPB,WPG_temp)
CALL dealloc_TabRDP(WPG_temp)

IF (debug) THEN
  write(out_unitp,*) '====================================='
  CALL time_perso('sub_G_TO_B')
  write(out_unitp,*) '====================================='
  write(out_unitp,*) '====================================='
  write(out_unitp,*) '====================================='
  flush(out_unitp)
END IF


END SUBROUTINE sub_G_TO_B
SUBROUTINE sub_B_TO_G(WPB,WPG,AllPara)
USE mod_system
IMPLICIT NONE

TYPE(TypeRDP), allocatable    :: WPG(:)
TYPE(TypeRDP), allocatable    :: WPB(:)

TYPE (Type_AllParam) :: AllPara

integer :: id
TYPE(TypeRDP), allocatable    :: WPB_temp(:)


integer :: nb_BG,nb_mult_id

logical, parameter :: debug=.FALSE.

IF (debug) THEN
  write(out_unitp,*) '====================================='
  write(out_unitp,*) '====================================='
  write(out_unitp,*) '====================================='
  CALL time_perso('sub_B_TO_G')
  write(out_unitp,*) '====================================='
  flush(out_unitp)
END IF

CALL TabRDP2_TO_TabRDP1(WPG,WPB)


DO id=AllPara%D,1,-1
  IF (debug) THEN
    CALL Size_TabRDP(WPG,nb_BG)
    write(out_unitp,*) 'id, size WPG',id,int(real(nb_BG,kind=Rkind)*EIGHT/(1024.d0**2)),' MB'
    flush(out_unitp)
  END IF

  CALL Transfer_BgG_TO_BbG(WPG,WPB_temp,                                &
                    AllPara%ind_Grid,AllPara%ind_Basis,AllPara%tab_ba,  &
                    AllPara%D,AllPara%LG,AllPara%LB,id)

  IF (debug) THEN
    CALL Size_TabRDP(WPB_temp,nb_BG)
    write(out_unitp,*) 'id, size WPB_temp',id,int(real(nb_BG,kind=Rkind)*EIGHT/(1024.d0**2)),' MB'
    write(out_unitp,*) 'id, nb_mult_id ',id,nb_mult_id
    flush(out_unitp)
  END IF

  CALL BbG_TO_BgG(WPB_temp,WPG,                                         &
                    AllPara%ind_Grid,AllPara%ind_Basis,AllPara%tab_ba,  &
                    AllPara%D,AllPara%LG,AllPara%LB,id,nb_mult_id)
  !write(out_unitp,*) 'sub_B_TO_G',id,' done' ; flush(out_unitp)

END DO



IF (debug) THEN
  write(out_unitp,*) '====================================='
  CALL time_perso('sub_B_TO_G')
  write(out_unitp,*) '====================================='
  write(out_unitp,*) '====================================='
  write(out_unitp,*) '====================================='
  flush(out_unitp)
END  IF


END SUBROUTINE sub_B_TO_G

SUBROUTINE V_ON_GRID(V_ON_G,AllPara)
USE mod_system
IMPLICIT NONE

TYPE(TypeRDP), allocatable    :: V_ON_G(:)
TYPE (Type_AllParam)          :: AllPara


integer            :: iG,i,nqq,iqq
real(kind=Rkind)   :: V,lambda
integer            :: tabiq(AllPara%D),tabnq(AllPara%D),tabl(AllPara%D)
real(kind=Rkind)   :: x(AllPara%D)

lambda = 0.111803_Rkind
!lambda = ZERO
CALL Set_BgG_FOR_id(V_ON_G,AllPara%ind_Grid,AllPara%ind_Basis,AllPara%tab_ba,AllPara%D,AllPara%LG,0)


DO iG=1,AllPara%ind_Grid(0)%MaxnD
  DO i=1,AllPara%D
    tabl(i)  = AllPara%ind_Grid(0)%tab_ind(i,iG)
    tabnq(i) = AllPara%tab_ba(tabl(i),i)%nq
  END DO
  nqq = product(tabnq)

  !write(out_unitp,*) 'iG,tabnq',iG,':',tabnq
  DO iqq=1,nqq
    CALL InD_TO_tabi(iqq,AllPara%D,tabnq,tabiq)
    DO i=1,AllPara%D
      x(i) = AllPara%tab_ba(tabl(i),i)%x(tabiq(i))
    END DO
    !write(out_unitp,*) 'iqq,tabiq',iqq,':',tabiq,'x:',x
    V = ZERO
    DO i=1,AllPara%D
      V = V + HALF * x(i)**2
    END DO

    !DO i=1,AllPara%D-1
    !  V = V + lambda * (x(i)**2 * x(i+1) - x(i+1)*3)
    !END DO

    V_ON_G(iG)%RDP(1,iqq,1) = V
  END DO
END DO

END SUBROUTINE V_ON_GRID
SUBROUTINE T1Psi_ON_GRID(T2WPG,WPG,iq_d1,AllPara)
USE mod_system
IMPLICIT NONE

TYPE(TypeRDP), allocatable    :: WPG(:),T2WPG(:)
TYPE (Type_AllParam)          :: AllPara
integer            :: iq_d1
real(kind=Rkind), allocatable :: RG(:,:,:)


integer            :: iG,i,nqq,li,iq1,nq1,iq,nq_d1,iq3,nq3,l_d1
integer            :: tabiq(AllPara%D),tabnq(AllPara%D),tabl(AllPara%D)


CALL Set_BgG_FOR_id(T2WPG,AllPara%ind_Grid,AllPara%ind_Basis,AllPara%tab_ba,AllPara%D,AllPara%LG,0)


DO iG=1,size(WPG)
  DO i=1,AllPara%D
    li  = AllPara%ind_Grid(0)%tab_ind(i,iG)
    tabnq(i) = AllPara%tab_ba(li,i)%nq
  END DO
  nqq = product(tabnq)
  IF (iq_d1 == 1) THEN
    nq1 = 1
  ELSE
    nq1 = product(tabnq(1:iq_d1-1))
  END IF
  nq_d1 = tabnq(iq_d1)
  l_d1  = AllPara%ind_Grid(0)%tab_ind(iq_d1,iG)
  IF (iq_d1 == AllPara%D) THEN
    nq3 = 1
  ELSE
    nq3 = product(tabnq(iq_d1+1:AllPara%D))
  END IF

  IF (nq1*nq_d1*nq3 /= nqq) THEN
    STOP 'ERROR with nqi'
  END IF

  allocate(RG(nq1,nq_d1,nq3))
  RG(:,:,:) = reshape(WPG(iG)%RDP,(/ nq1,nq_d1,nq3 /))

  DO iq1=1,nq1
  DO iq3=1,nq3

    RG(iq1,:,iq3) = matmul(AllPara%tab_ba(l_d1,iq_d1)%d1bGG,RG(iq1,:,iq3))

  END DO
  END DO

  T2WPG(iG)%RDP = reshape(RG,(/ 1,nqq,1 /))
  deallocate(RG)

END DO

END SUBROUTINE T1Psi_ON_GRID
SUBROUTINE T11Psi_ON_GRID(T2WPG,WPG,iq1_d1,iq2_d1,AllPara)
USE mod_system
IMPLICIT NONE

TYPE(TypeRDP), allocatable    :: WPG(:),T2WPG(:)
TYPE (Type_AllParam)          :: AllPara
integer            :: iq1_d1,iq2_d1
real(kind=Rkind), allocatable :: RG(:,:,:)


integer            :: iG,i,nqq,li,iq1,nq1,iq,nq_d1,iq3,nq3,l_d1,iq_d1
integer            :: tabiq(AllPara%D),tabnq(AllPara%D),tabl(AllPara%D)


CALL Set_BgG_FOR_id(T2WPG,AllPara%ind_Grid,AllPara%ind_Basis,AllPara%tab_ba,AllPara%D,AllPara%LG,0)


DO iG=1,size(WPG)
  DO i=1,AllPara%D
    li       = AllPara%ind_Grid(0)%tab_ind(i,iG)
    tabnq(i) = AllPara%tab_ba(li,i)%nq
  END DO
  nqq = product(tabnq)

  ! first derivative: iq1_d1
  iq_d1 = iq1_d1

  IF (iq_d1 == 1) THEN
    nq1 = 1
  ELSE
    nq1 = product(tabnq(1:iq_d1-1))
  END IF
  nq_d1 = tabnq(iq_d1)
  l_d1  = AllPara%ind_Grid(0)%tab_ind(iq_d1,iG)
  IF (iq_d1 == AllPara%D) THEN
    nq3 = 1
  ELSE
    nq3 = product(tabnq(iq_d1+1:AllPara%D))
  END IF

  IF (nq1*nq_d1*nq3 /= nqq) THEN
    STOP 'ERROR with nqi'
  END IF

  allocate(RG(nq1,nq_d1,nq3))
  RG(:,:,:) = reshape(WPG(iG)%RDP,(/ nq1,nq_d1,nq3 /))

  DO iq1=1,nq1
  DO iq3=1,nq3

    RG(iq1,:,iq3) = matmul(AllPara%tab_ba(l_d1,iq_d1)%d1bGG,RG(iq1,:,iq3))

  END DO
  END DO

  T2WPG(iG)%RDP = reshape(RG,(/ 1,nqq,1 /))
  deallocate(RG)

  ! second derivative: iq2_d1
  iq_d1 = iq2_d1

  IF (iq_d1 == 1) THEN
    nq1 = 1
  ELSE
    nq1 = product(tabnq(1:iq_d1-1))
  END IF
  nq_d1 = tabnq(iq_d1)
  l_d1  = AllPara%ind_Grid(0)%tab_ind(iq_d1,iG)
  IF (iq_d1 == AllPara%D) THEN
    nq3 = 1
  ELSE
    nq3 = product(tabnq(iq_d1+1:AllPara%D))
  END IF

  IF (nq1*nq_d1*nq3 /= nqq) THEN
    STOP 'ERROR with nqi'
  END IF

  allocate(RG(nq1,nq_d1,nq3))
  RG(:,:,:) = reshape(WPG(iG)%RDP,(/ nq1,nq_d1,nq3 /))

  DO iq1=1,nq1
  DO iq3=1,nq3

    RG(iq1,:,iq3) = matmul(AllPara%tab_ba(l_d1,iq_d1)%d1bGG,RG(iq1,:,iq3))

  END DO
  END DO

  T2WPG(iG)%RDP = reshape(RG,(/ 1,nqq,1 /))
  deallocate(RG)

END DO

END SUBROUTINE T11Psi_ON_GRID

SUBROUTINE T2Psi_ON_GRID(T2WPG,WPG,iq_d2,AllPara)
USE mod_system
IMPLICIT NONE

TYPE(TypeRDP), allocatable    :: WPG(:),T2WPG(:)
TYPE (Type_AllParam)          :: AllPara
integer            :: iq_d2
real(kind=Rkind), allocatable :: RG(:,:,:)


integer            :: iG,i,nqq,li,iq1,nq1,iq,nq_d2,iq3,nq3,l_d2
integer            :: tabiq(AllPara%D),tabnq(AllPara%D),tabl(AllPara%D)


CALL Set_BgG_FOR_id(T2WPG,AllPara%ind_Grid,AllPara%ind_Basis,AllPara%tab_ba,AllPara%D,AllPara%LG,0)


DO iG=1,size(WPG)
  DO i=1,AllPara%D
    li  = AllPara%ind_Grid(0)%tab_ind(i,iG)
    tabnq(i) = AllPara%tab_ba(li,i)%nq
  END DO
  nqq = product(tabnq)
  IF (iq_d2 == 1) THEN
    nq1 = 1
  ELSE
    nq1 = product(tabnq(1:iq_d2-1))
  END IF
  nq_d2 = tabnq(iq_d2)
  l_d2  = AllPara%ind_Grid(0)%tab_ind(iq_d2,iG)
  IF (iq_d2 == AllPara%D) THEN
    nq3 = 1
  ELSE
    nq3 = product(tabnq(iq_d2+1:AllPara%D))
  END IF

  IF (nq1*nq_d2*nq3 /= nqq) THEN
    STOP 'ERROR with nqi'
  END IF

  allocate(RG(nq1,nq_d2,nq3))
  RG(:,:,:) = reshape(WPG(iG)%RDP,(/ nq1,nq_d2,nq3 /))

  DO iq1=1,nq1
  DO iq3=1,nq3

    RG(iq1,:,iq3) = matmul(AllPara%tab_ba(l_d2,iq_d2)%d2bGG,RG(iq1,:,iq3))

  END DO
  END DO

  T2WPG(iG)%RDP = reshape(RG,(/ 1,nqq,1 /))
  deallocate(RG)

END DO

END SUBROUTINE T2Psi_ON_GRID

END MODULE mod_Smolyak_test
