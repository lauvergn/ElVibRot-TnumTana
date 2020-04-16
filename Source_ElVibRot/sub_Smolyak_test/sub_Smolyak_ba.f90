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
MODULE mod_Smolyak_ba
use mod_Smolyak_DInd, only: l_to_n
use mod_system
IMPLICIT NONE

integer :: Bmin=1

TYPE TypeBa
  integer :: nq = 0
  integer :: nb = 0
  real(kind=Rkind), allocatable :: x(:)
  real(kind=Rkind), allocatable :: w(:)
  real(kind=Rkind), allocatable :: d0b(:,:)
  real(kind=Rkind), allocatable :: twd0b(:,:)
  real(kind=Rkind), allocatable :: td0b(:,:)



  real(kind=Rkind), allocatable :: d1b(:,:)
  real(kind=Rkind), allocatable :: d2b(:,:)

  real(kind=Rkind), allocatable :: d1bGG(:,:)
  real(kind=Rkind), allocatable :: d2bGG(:,:)

  integer,          allocatable :: ib_TO_lb(:)
  integer                       :: Nested = 0 ! default not nested
  integer                       :: nq_max_Nested = 0

END TYPE TypeBa

PRIVATE
PUBLIC :: TypeBa, alloc_TypeBa, dealloc_TypeBa, Set_ba, Set_Delatba, Set_Delatba_nested1, &
          Set_tab_Ba,Set_tab_DelatBa, get_tab_nq, get_tab_nb

CONTAINS

SUBROUTINE alloc_TypeBa(Ba,nb,nq)
!USE mod_system
IMPLICIT NONE

integer      :: nb,nq
TYPE(TypeBa) :: Ba

CALL dealloc_TypeBa(Ba)

Ba%nb = nb
Ba%nq = nq

allocate(Ba%x(nq))
allocate(Ba%w(nq))

allocate(Ba%d0b(nq,nb))
allocate(Ba%d1b(nq,nb))
allocate(Ba%d2b(nq,nb))

allocate(Ba%twd0b(nb,nq))
allocate(Ba%td0b(nb,nq))


allocate(Ba%d1bGG(nq,nq))
allocate(Ba%d2bGG(nq,nq))

allocate(Ba%ib_TO_lb(nb))


END SUBROUTINE alloc_TypeBa
SUBROUTINE dealloc_TypeBa(Ba)
!USE mod_system
IMPLICIT NONE

TYPE(TypeBa) :: Ba

IF (allocated(Ba%x))        deallocate(Ba%x)
IF (allocated(Ba%w))        deallocate(Ba%w)

IF (allocated(Ba%d0b))      deallocate(Ba%d0b)
IF (allocated(Ba%d1b))      deallocate(Ba%d1b)
IF (allocated(Ba%d2b))      deallocate(Ba%d2b)

IF (allocated(Ba%twd0b))    deallocate(Ba%twd0b)


IF (allocated(Ba%d1bGG))    deallocate(Ba%d1bGG)
IF (allocated(Ba%d2bGG))    deallocate(Ba%d2bGG)

IF (allocated(Ba%ib_TO_lb)) deallocate(Ba%ib_TO_lb)

Ba%nb     = 0
Ba%nq     = 0

END SUBROUTINE dealloc_TypeBa


SUBROUTINE dealloc_tab_ba(tab_ba)
!USE mod_system
IMPLICIT NONE

TYPE(TypeBa), allocatable :: tab_ba(:,:)

integer :: id,l

IF (allocated(tab_ba)) THEN
  DO id=lbound(tab_ba,dim=2),ubound(tab_ba,dim=2)
  DO l =lbound(tab_ba,dim=1),ubound(tab_ba,dim=1)
    CALL dealloc_TypeBa(tab_ba(l,id))
  END DO
  END DO
  deallocate(tab_ba)
END IF

END SUBROUTINE dealloc_tab_ba

SUBROUTINE Set_ba(ba,nb,nq,l,LB,B)
!USE mod_system
IMPLICIT NONE

integer       :: nb,nq
TYPE(TypeBa)  :: ba


integer          :: id,l,LB,nql,nqlm1,iq,ib,jb,ilb,B
real(kind=Rkind) :: max_S,S
real(kind=Rkind) :: poly_Hermite_exp ! function

real(kind=Rkind) :: step,d0,d1,d2
logical          :: num,deriv

real(kind=Rkind), allocatable :: d0bxd0bT(:,:),d0bxd0bT_inv(:,:),d0b_pseudoInv(:,:)
real(kind=Rkind), allocatable :: Check_bGB(:,:)     ! (nb,nq)
real(kind=Rkind), allocatable :: d0b(:,:)
!logical                       :: new = .FALSE.
logical                       :: new = .TRUE.
logical                       :: debug = .FALSE.


num   =.TRUE.
deriv =.TRUE.
step  = ONETENTH**4

! the 1D-basis
  !write(out_unitp,*) 'B',B

  CALL alloc_TypeBa(ba,nb,nq)

  DO ilb=min(l,LB),0,-1
    ib = l_TO_n(ilb,1,B=B)
    ba%ib_TO_lb(1:ib) = ilb
  END DO

  IF (debug) THEN
    write(out_unitp,*) '==== ba =========== nb,nq',nb,nq
    write(out_unitp,*) '  ib_TO_lb',ba%ib_TO_lb(:)
  END IF

  CALL hercom(nq,ba%x,ba%w)
  !write(out_unitp,*) 'l,x',ba%x
  !write(out_unitp,*) 'l,w',ba%w
  DO ib=1,nb
  DO iq=1,nq
    CALL d0d1d2poly_Hermite_exp(ba%x(iq),ib-1,d0,d1,d2,deriv,num,step)
    ba%d0b(iq,ib) = d0
    ba%d1b(iq,ib) = d1
    ba%d2b(iq,ib) = d2

    ba%twd0b(ib,iq) = d0 * ba%w(iq)
    ba%td0b(ib,iq) = d0

  END DO
  END DO

  IF (debug) write(out_unitp,*) '===================== ortho ?'
  max_S = ZERO
  DO ib=1,nb
  DO jb=1,nb
    S = dot_product(ba%d0b(:,ib),ba%twd0b(:,jb))
    IF (ib == jb) S=S-ONE
    IF (abs(S) > max_S) max_S = S
    !write(out_unitp,*) 'ib,jb',ib,jb,S
  END DO
  END DO
  IF (debug) write(out_unitp,*) '===================== max_S:',max_S

  IF (debug) write(out_unitp,*) '===================== Set d1bGG and d2bGG'
  allocate(Check_bGB(nb,nq))

  IF (new) THEN
    allocate(d0b(nq,nq))
    allocate(d0b_pseudoInv(nq,nq))

    d0b = ba%d0b
    CALL inv_m1_TO_m2(d0b,d0b_pseudoInv,nq,0,ZERO) ! not SVD
    !write(out_unitp,*) 'd0b_pseudoInv'
    !CALL Write_VecMat(d0b_pseudoInv,out_unitp,5)
    !flush(out_unitp)
    deallocate(d0b)

  ELSE
    allocate(d0bxd0bT(nq,nq))
    allocate(d0bxd0bT_inv(nq,nq))
    allocate(d0b_pseudoInv(nb,nq))

    d0bxd0bT = matmul(ba%d0b,transpose(ba%d0b))
    !write(out_unitp,*) 'd0bxd0bT'
    !CALL Write_VecMat(d0bxd0bT,out_unitp,5)

    CALL inv_m1_TO_m2(d0bxd0bT,d0bxd0bT_inv,nq,1,ONETENTH**10) ! SVD
    !write(out_unitp,*) 'd0bxd0bT_inv'
    !CALL Write_VecMat(d0bxd0bT_inv,out_unitp,5)

    d0b_pseudoInv =  matmul(transpose(ba%d0b),d0bxd0bT_inv)
    !write(out_unitp,*) 'd0b_pseudoInv'
    !CALL Write_VecMat(d0b_pseudoInv,out_unitp,5)

    deallocate(d0bxd0bT)
    deallocate(d0bxd0bT_inv)
  END IF

  ba%d1bGG =  matmul(ba%d1b,d0b_pseudoInv(1:nb,:))
  ba%d2bGG =  matmul(ba%d2b,d0b_pseudoInv(1:nb,:))

  Check_bGB = ba%d1b-matmul(ba%d1bGG,ba%d0b)
  IF (debug) write(out_unitp,*) 'WARNING Check_bGB%d1',maxval(abs(Check_bGB))
  IF (maxval(abs(Check_bGB)) > ONETENTH**10) write(out_unitp,*) 'WARNING Check_bGB%d1',maxval(abs(Check_bGB))

  Check_bGB = ba%d2b-matmul(ba%d2bGG,ba%d0b)
  IF (debug) write(out_unitp,*) 'WARNING Check_bGB%d2',maxval(abs(Check_bGB))
  IF (maxval(abs(Check_bGB)) > ONETENTH**10) write(out_unitp,*) 'WARNING Check_bGB%d2',maxval(abs(Check_bGB))

  deallocate(Check_bGB)
  deallocate(d0b_pseudoInv)

END SUBROUTINE Set_ba

SUBROUTINE Set_DelatBa(DelatBa,nb1,nq1,nb2,nq2,l,D,LB,LG)
!USE mod_system
IMPLICIT NONE

integer       :: D,LB,LG
TYPE(TypeBa)  :: DelatBa


integer          :: id,l,nb,nq,nql,nqlm1,iq,ib,jb,ilb,B
integer          :: nb1,nq1,nb2,nq2

real(kind=Rkind) :: max_S,S
real(kind=Rkind) :: poly_Hermite_exp ! function

real(kind=Rkind) :: step,d0,d1,d2
logical          :: num,deriv


!logical                       :: debug = .FALSE.
logical                       :: debug = .TRUE.


num   =.TRUE.
deriv =.TRUE.
step  = ONETENTH**4

! the 1D-basis

  nb = nb2-nb1
  nq = nq1+nq2

  CALL alloc_TypeBa(DelatBa,nb,nq)

  DelatBa%ib_TO_lb(:) = l


  IF (debug) THEN
    write(out_unitp,*) '==== tab_DelatBa =========== nb,nq',nb,nq
    write(out_unitp,*) '  ib_TO_lb',DelatBa%ib_TO_lb(:)
  END IF

  IF (nq1 > 0) THEN
    CALL hercom(nq1,DelatBa%x(1:nq1),DelatBa%w(1:nq1))
    DelatBa%w(1:nq1) = -DelatBa%w(1:nq1)
  END IF
  CALL hercom(nq2,DelatBa%x(nq1+1:nq),DelatBa%w(nq1+1:nq))
  !write(out_unitp,*) 'l,x',DelatBa%x
  !write(out_unitp,*) 'l,w',DelatBa%w
  DO ib=1,nb
  DO iq=1,nq
    CALL d0d1d2poly_Hermite_exp(DelatBa%x(iq),nb1+ib-1,d0,d1,d2,deriv,num,step)
    DelatBa%d0b(iq,ib) = d0
    DelatBa%d1b(iq,ib) = d1
    DelatBa%d2b(iq,ib) = d2

    DelatBa%twd0b(ib,iq) = d0 * DelatBa%w(iq)

  END DO
  END DO

  !IF (debug) write(out_unitp,*) '===================== d0b'
  !IF (debug) CALL Write_Mat(DelatBa%d0b,6,5)

  !IF (debug) write(out_unitp,*) '===================== ortho ?'
  max_S = ZERO
  DO ib=1,nb
  DO jb=1,nb
    S = dot_product(DelatBa%d0b(:,ib),DelatBa%twd0b(jb,:))
    IF (ib == jb) S=S-ONE
    IF (abs(S) > max_S) max_S = S
    !write(out_unitp,*) 'ib,jb',ib,jb,S
  END DO
  END DO
  IF (debug) write(out_unitp,*) '===================== max_S:',max_S


END SUBROUTINE Set_Delatba

SUBROUTINE Set_Delatba_nested1(DelatBa,nb1,nq1,nb2,nq2,l,D,LB,LG)
!USE mod_system
IMPLICIT NONE

integer       :: D,LB,LG
TYPE(TypeBa)  :: DelatBa


integer          :: id,l,nb,nq,nql,nqlm1,iq,ib,jb,ilb,B
integer          :: nb1,nq1,nb2,nq2,dnq,nq0

real(kind=Rkind) :: max_S,S
real(kind=Rkind) :: poly_Hermite_exp ! function

real(kind=Rkind) :: x1(nq1),w1(nq1)

real(kind=Rkind) :: step,d0,d1,d2
logical          :: num,deriv

logical                       :: debug = .TRUE.


num   =.TRUE.
deriv =.TRUE.
step  = ONETENTH**4

! the 1D-basis

  nb = nb2-nb1
  nq = max(nq1,nq2)
  IF (l > 1 .AND. 2*nb /= (nq2-nq1) ) THEN

    write(out_unitp,*) 'l,D,LB,LG',l,D,LB,LG

    write(out_unitp,*) 'nb1,nb2,nb',nb1,nb2,nb
    write(out_unitp,*) 'nq1,nq2,nq',nq1,nq2,nq

    STOP 'ERROR in Set_Delatba_nested1'
  END IF

  CALL alloc_TypeBa(DelatBa,nb,nq)

  DelatBa%nq_max_Nested = 2*LG+1

  CALL grid_HermiteNested1(DelatBa%x,DelatBa%w,nq2,DelatBa%nq_max_Nested)
  IF (nq1 > 0) THEN
    CALL grid_HermiteNested1(x1,w1,nq1,DelatBa%nq_max_Nested)
    write(out_unitp,*) 'max diff x1 x2',maxval(DelatBa%x(1+nb:nq2-nb)-x1)
    DelatBa%w(1+nb:nq2-nb) = DelatBa%w(1+nb:nq2-nb) - w1
  END IF


  DelatBa%ib_TO_lb(:) = l


  IF (debug) THEN
    write(out_unitp,*) '==== tab_DelatBa =========== nb,nq',nb,nq
    write(out_unitp,*) '  ib_TO_lb',DelatBa%ib_TO_lb(:)
  END IF

  !write(out_unitp,*) 'l,x',DelatBa%x
  !write(out_unitp,*) 'l,w',DelatBa%w
  DO ib=1,nb
  DO iq=1,nq
    CALL d0d1d2poly_Hermite_exp(DelatBa%x(iq),nb1+ib-1,d0,d1,d2,deriv,num,step)
    DelatBa%d0b(iq,ib) = d0
    DelatBa%d1b(iq,ib) = d1
    DelatBa%d2b(iq,ib) = d2

    DelatBa%twd0b(ib,iq) = d0 * DelatBa%w(iq)

  END DO
  END DO

  !IF (debug) write(out_unitp,*) '===================== d0b'
  !IF (debug) CALL Write_Mat(DelatBa%d0b,6,5)

  !IF (debug) write(out_unitp,*) '===================== ortho ?'
  max_S = ZERO
  DO ib=1,nb
  DO jb=1,nb
    S = dot_product(DelatBa%d0b(:,ib),DelatBa%twd0b(jb,:))
    IF (ib == jb) S=S-ONE
    IF (abs(S) > max_S) max_S = S
    !write(out_unitp,*) 'ib,jb',ib,jb,S
  END DO
  END DO
  IF (debug) write(out_unitp,*) '===================== max_S:',max_S


END SUBROUTINE Set_Delatba_nested1

SUBROUTINE Set_tab_Ba(tab_ba,D,LB,LG)
!USE mod_system
IMPLICIT NONE

integer       :: D,LB,LG
TYPE(TypeBa), allocatable :: tab_ba(:,:)


integer          :: id,l,nb,nq,B


!-------------------------------------------
!-------------------------------------------
IF (.NOT. allocated(tab_ba)) allocate(tab_ba(0:LG,D))
!-------------------------------------------
!-------------------------------------------

! the 1D-basis
DO id=1,D
  write(out_unitp,*) 'id,B',id,min(id,Bmin)

DO l=0,ubound(tab_ba,dim=1)
  B = min(id,Bmin)
  nq =  l_TO_n(l,1,B=B)
  nb =  l_TO_n(min(LB,l),1,B=B)

  CALL Set_ba(tab_ba(l,id),nb,nq,l,LB,B)

END DO
END DO

END SUBROUTINE Set_tab_Ba

SUBROUTINE Set_tab_DelatBa(tab_DelatBa,D,LB,LG)
!USE mod_system
IMPLICIT NONE

integer       :: D,LB,LG
TYPE(TypeBa), allocatable :: tab_DelatBa(:,:)


integer          :: id,l,nb,nq,nql,nqlm1,iq,ib,jb,ilb,B
integer          :: nb1,nq1,nb2,nq2

real(kind=Rkind) :: max_S,S
real(kind=Rkind) :: poly_Hermite_exp ! function

real(kind=Rkind) :: step,d0,d1,d2
logical          :: num,deriv


logical                       :: nested = .TRUE.
logical                       :: debug = .FALSE.

!-------------------------------------------
!-------------------------------------------
IF (.NOT. allocated(tab_DelatBa)) allocate(tab_DelatBa(0:LG,D))
!-------------------------------------------
!-------------------------------------------
num   =.TRUE.
deriv =.TRUE.
step  = ONETENTH**4

! the 1D-basis
DO id=1,D
DO l=0,ubound(tab_DelatBa,dim=1)
  IF (nested) THEN

    IF (l>0) THEN
      nq1 =  l_TO_n(l-1,1,B=2)
      nb1 =  l_TO_n(min(LB,l-1),1,B=1)
    ELSE
      nq1 = 0
      nb1 = 0
    END IF

    nq2 =  l_TO_n(l,1,B=2)
    nb2 =  l_TO_n(min(LB,l),1,B=1)

    CALL Set_DelatBa_nested1(tab_DelatBa(l,id),nb1,nq1,nb2,nq2,l,D,LB,LG)

  ELSE

    B = min(id,Bmin)

    IF (l>0) THEN
      nq1 =  l_TO_n(l-1,1,B=B)
      nb1 =  l_TO_n(min(LB,l-1),1,B=B)
    ELSE
      nq1 = 0
      nb1 = 0
    END IF

    nq2 =  l_TO_n(l,1,B=B)
    nb2 =  l_TO_n(min(LB,l),1,B=B)

    CALL Set_DelatBa(tab_DelatBa(l,id),nb1,nq1,nb2,nq2,l,D,LB,LG)
  END IF

END DO
END DO

END SUBROUTINE Set_tab_Delatba

FUNCTION get_tab_nq(tab_l,tab_ba) RESULT (tab_nq)
  integer,         allocatable                           :: tab_nq(:) ! result
  integer ,                        intent(in)            :: tab_l(:)
  TYPE(TypeBa),    allocatable,    intent(in)            :: tab_ba(:,:) ! tab_ba(L,D)

  integer :: i

  tab_nq = tab_l
  DO i=1,size(tab_l)
    tab_nq(i) = tab_ba(tab_l(i),i)%nq
  END DO

END FUNCTION get_tab_nq
FUNCTION get_tab_nb(tab_l,tab_ba) RESULT (tab_nb)
  integer,         allocatable                           :: tab_nb(:) ! result
  integer ,                        intent(in)            :: tab_l(:)
  TYPE(TypeBa),    allocatable,    intent(in)            :: tab_ba(:,:) ! tab_ba(L,D)

  integer :: i

  tab_nb = tab_l
  DO i=1,size(tab_l)
    tab_nb(i) = tab_ba(tab_l(i),i)%nb
    !tab_nb(size(tab_l)+1-i) = tab_ba(tab_l(i),i)%nb
  END DO

END FUNCTION get_tab_nb
END MODULE mod_Smolyak_ba
