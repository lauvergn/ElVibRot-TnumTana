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
MODULE mod_module_DInd
use mod_system
use mod_dnSVM, only: type_intvec, write_intvec
IMPLICIT NONE

PRIVATE

TYPE TypeDInd
  integer              :: ndim  = 0
  integer              :: MaxnD = -1
  integer, allocatable :: tab_ind(:,:)
  integer, allocatable :: indD_OF_Dm1(:)
  integer, allocatable :: i_TO_l(:)         ! give the l value for i (usefull when i /= l+1)
  integer, allocatable :: lmax_TO_nb(:)     ! give the nb value for l (number of terms with l<=lmax)

  integer, allocatable :: tab_q(:)          ! size: ndim
CONTAINS
  PROCEDURE, PRIVATE, PASS(DInd1) :: TypeDInd2TOTypeDInd1
  GENERIC,   PUBLIC  :: assignment(=) => TypeDInd2TOTypeDInd1
END TYPE TypeDInd

PUBLIC :: TypeDInd, alloc_TypeDInd, dealloc_TypeDInd, Write_TypeDInd
PUBLIC :: Write_Tab_nDInd,dealloc_nDInd,nDInd2TOnDInd1
PUBLIC :: Set_nDInd_01order,Set_nDInd_10order, Set_nDInd_01order_L,Set_nDInd_10order_L
PUBLIC :: InD_TO_tabi,tabi_TO_InD

CONTAINS

ELEMENTAL SUBROUTINE alloc_TypeDInd(DInd,ndim,MaxnD)
IMPLICIT NONE

integer,        intent(in)    :: ndim,MaxnD
TYPE(TypeDInd), intent(inout) :: DInd


CALL dealloc_TypeDInd(DInd)

DInd%ndim  = ndim
DInd%MaxnD = MaxnD
allocate(DInd%tab_ind(ndim,MaxnD))
allocate(DInd%i_TO_l(MaxnD))  ! give the l value for i (usefull when i /= l+1)
DInd%i_TO_l(:) = 0
allocate(DInd%indD_OF_Dm1(MaxnD))
allocate(DInd%tab_q(ndim))

! lmax_TO_nb will be allocated later

END SUBROUTINE alloc_TypeDInd
ELEMENTAL SUBROUTINE dealloc_TypeDInd(DInd)
IMPLICIT NONE

TYPE(TypeDInd), intent(inout) :: DInd


DInd%ndim  = 0
DInd%MaxnD = -1
IF (allocated(DInd%tab_ind))          deallocate(DInd%tab_ind)
IF (allocated(DInd%i_TO_l))           deallocate(DInd%i_TO_l)
IF (allocated(DInd%lmax_TO_nb))       deallocate(DInd%lmax_TO_nb)
IF (allocated(DInd%indD_OF_Dm1))      deallocate(DInd%indD_OF_Dm1)
IF (allocated(DInd%tab_q))            deallocate(DInd%tab_q)


END SUBROUTINE dealloc_TypeDInd
ELEMENTAL SUBROUTINE TypeDInd2TOTypeDInd1(DInd1,DInd2)
IMPLICIT NONE

CLASS(TypeDInd), intent(inout) :: DInd1
TYPE(TypeDInd),  intent(in)     :: DInd2

integer :: lmax

CALL alloc_TypeDInd(DInd1,DInd2%ndim,DInd2%MaxnD)
DInd1%tab_ind     = DInd2%tab_ind
DInd1%i_TO_l      = DInd2%i_TO_l
DInd1%indD_OF_Dm1 = DInd2%indD_OF_Dm1
DInd1%tab_q       = DInd2%tab_q

IF (allocated(DInd2%lmax_TO_nb)) THEN
  lmax = ubound(DInd2%lmax_TO_nb,dim=1)
  allocate(DInd1%lmax_TO_nb(0:lmax))
  DInd1%lmax_TO_nb = DInd2%lmax_TO_nb
END IF

END SUBROUTINE TypeDInd2TOTypeDInd1
SUBROUTINE Write_TypeDInd(DInd)
IMPLICIT NONE

TYPE(TypeDInd), intent(in) :: DInd

integer :: I

write(out_unitp,*) 'BEGINNING Write_TypeDind'
write(out_unitp,*) 'ndim',DInd%ndim
write(out_unitp,*) 'MaxnD',DInd%MaxnD
write(out_unitp,*) 'tab_q(:)',DInd%tab_q(:)

IF (allocated(DInd%lmax_TO_nb)) write(out_unitp,*) 'lmax_TO_nb(:)',DInd%lmax_TO_nb(:)

write(out_unitp,*) '         I,        L,   indD_OF_Dm1,   ind(:)'
DO I=1,DInd%MaxnD
  write(out_unitp,*) I,DInd%i_TO_l(I),DInd%indD_OF_Dm1(I),':',DInd%tab_ind(:,I)
END DO
write(out_unitp,*) 'END Write_TypeDind'
flush(out_unitp)

END SUBROUTINE Write_TypeDInd

SUBROUTINE Set_nDInd_01order(nDind,D,Lmin,Lmax,tab_i_TO_l)
IMPLICIT NONE

integer                         :: D,Lmin,Lmax
TYPE (TypeDInd), allocatable    :: nDind(:)
TYPE (Type_IntVec), pointer     :: tab_i_TO_l(:)

integer :: i,id,iGm1,iG,nG,l,ll,lll,ndimGm1,n
integer, allocatable :: i_TO_l(:)
logical :: test

logical, parameter :: debug=.FALSE.
!logical, parameter :: debug=.TRUE.

IF (.NOT. allocated(nDind)) allocate(nDind(0:D+1))

IF (.NOT. associated(tab_i_TO_l)) STOP 'ERROR in Set_nDInd_01order: tab_i_TO_l is not associated'


IF (associated(tab_i_TO_l) .AND. debug) THEN
  write(out_unitp,*) 'Lmin,Lmax',Lmin,Lmax
  write(out_unitp,*) 'tab_i_TO_l'
  DO id=1,D
    CALL Write_IntVec(tab_i_TO_l(id))
  END DO
END IF
!-------------------------------------------
!-------------------------------------------
!-------------------------------------------
!-------------------------------------------
id=0
CALL alloc_TypeDInd(nDind(id),ndim=1,MaxnD=1)
nDind(id)%tab_ind(1,1)   = 0
nDind(id)%i_TO_l(1)      = 0
nDind(id)%indD_OF_Dm1(1) = 1
nDind(id)%tab_q(1)       = 0

id=1
CALL alloc_TypeDInd(nDind(id),ndim=1,MaxnD=1)
nDind(id)%tab_ind(1,1)   = 0
nDind(id)%i_TO_l(1)      = 0
nDind(id)%indD_OF_Dm1(1) = 1
nDind(id)%tab_q(1)       = 0


DO id=2,D+1
  ! first the number of points
  n = tab_i_TO_l(id-1)%nb_var_vec
  allocate(i_TO_l(n))
  i_TO_l(:) = tab_i_TO_l(id-1)%vec(:)

  nG = 0
  DO iGm1=1,nDind(id-1)%MaxnD
    ll = nDind(id-1)%i_TO_l(iGm1)
    DO i=1,n
      l  = i_TO_l(i)
      lll=ll+l
      test = lll <= Lmax
      IF (id==D+1) test = test .AND. lll >= Lmin
      IF (test) nG = nG + 1
      !write(out_unitp,*) 'iGm1,i,ll,l,test,nG',iGm1,i,ll,l,test,nG
    END DO
  END DO
  CALL alloc_TypeDInd(nDind(id),ndim=id-1,MaxnD=nG)
  nDind(id)%tab_q(:) = (/ (i,i=1,id-1) /)


  iG = 0
  ndimGm1 = nDind(id-1)%ndim

  DO iGm1=1,nDind(id-1)%MaxnD
    ll = nDind(id-1)%i_TO_l(iGm1)
    DO i=1,n
      l  = i_TO_l(i)
      lll=ll+l
      test = lll <= Lmax
      IF (id==D+1) test = test .AND. lll >= Lmin
      IF (test) THEN
        iG = iG + 1
        nDind(id)%tab_ind(1:ndimGm1,iG)      = nDind(id-1)%tab_ind(:,iGm1)
        nDind(id)%tab_ind(nDind(id)%ndim,iG) = i
        nDind(id)%i_TO_l(iG)                 = lll
        nDind(id)%indD_OF_Dm1(iG)            = iGm1
        !write(out_unitp,*) 'id,iG,l(:)',id,iG,nDind(id)%tab_ind(:,iG),nDind(id)%indD_OF_Dm1(iG) ; flush(out_unitp)
      END IF
    END DO
  END DO
  IF (debug) THEN
    write(out_unitp,*) '======================================='
    write(out_unitp,*) 'id,tab_q ',id,':',nDind(id)%tab_q
    write(out_unitp,*) 'id,MaxnD ',id,':',nDind(id)%MaxnD
    write(out_unitp,*) 'id,i_TO_l',id,':',i_TO_l(:)
    CALL Write_TypeDInd(nDind(id))
    flush(out_unitp)
  END IF
  deallocate(i_TO_l)


END DO

END SUBROUTINE Set_nDInd_01order
SUBROUTINE Set_nDInd_10order(nDind,D,Lmin,Lmax,tab_i_TO_l)
IMPLICIT NONE

integer                         :: D,Lmin,Lmax
TYPE (TypeDInd), allocatable    :: nDind(:)
TYPE (Type_IntVec), pointer     :: tab_i_TO_l(:)

integer :: i,id,iGp1,iG,nG,l,ll,lll,ndimGp1,n
integer, allocatable :: i_TO_l(:)
logical :: test

logical, parameter :: debug=.FALSE.
!logical, parameter :: debug=.TRUE.

IF (.NOT. allocated(nDind)) allocate(nDind(0:D+1))

IF (.NOT. associated(tab_i_TO_l)) STOP 'ERROR in Set_nDInd_10order: tab_i_TO_l is not associated'

IF (associated(tab_i_TO_l) .AND. debug) THEN
  write(out_unitp,*) 'tab_i_TO_l'
  DO id=1,D
    CALL Write_IntVec(tab_i_TO_l(id))
  END DO
END IF
!-------------------------------------------
!-------------------------------------------
!-------------------------------------------
!-------------------------------------------
id=D+1
CALL alloc_TypeDInd(nDind(id),ndim=1,MaxnD=1)
nDind(id)%tab_ind(1,1)   = 0
nDind(id)%i_TO_l(1)      = 0
nDind(id)%indD_OF_Dm1(1) = 1
nDind(id)%tab_q(1)       = 0

id=D
nDind(id)%ndim  = 1
nDind(id)%MaxnD = 1
CALL alloc_TypeDInd(nDind(id),ndim=1,MaxnD=1)
nDind(id)%tab_ind(1,1)   = 0
nDind(id)%i_TO_l(1)      = 0
nDind(id)%indD_OF_Dm1(1) = 1
nDind(id)%tab_q(1)       = 0


DO id=D-1,0,-1
  ! first the number of points
  n = tab_i_TO_l(id+1)%nb_var_vec
  allocate(i_TO_l(n))
  i_TO_l(:) = tab_i_TO_l(id+1)%vec(:)

  nG = 0
  DO iGp1=1,nDind(id+1)%MaxnD
    ll = nDind(id+1)%i_TO_l(iGp1)
    DO i=1,n
      l  = i_TO_l(i)
      lll=ll+l
      test = lll <= Lmax
      IF (id==0) test = test .AND. lll >= Lmin
      IF (test) nG = nG + 1
    END DO
  END DO
  CALL alloc_TypeDInd(nDind(id),ndim=D-id,MaxnD=nG)
  nDind(id)%tab_q(:) = (/ (i,i=id+1,D) /)


  iG      = 0
  ndimGp1 = nDind(id+1)%ndim

  DO iGp1=1,nDind(id+1)%MaxnD
    !ll = sum(nDind(id+1)%tab_ind(:,iGp1))
    ll = nDind(id+1)%i_TO_l(iGp1)
    DO i=1,n
      l  = i_TO_l(i)
      lll=ll+l
      test = lll <= Lmax
      IF (id==0) test = test .AND. lll >= Lmin
      IF (test) THEN
        iG = iG + 1
        IF (id < D-1) nDind(id)%tab_ind(2:1+ndimGp1,iG) = nDind(id+1)%tab_ind(:,iGp1)
        nDind(id)%tab_ind(1,iG)              = i
        nDind(id)%i_TO_l(iG)                 = lll
        nDind(id)%indD_OF_Dm1(iG)            = iGp1
        !write(out_unitp,*) 'id,iG,l(:)',id,iG,nDind(id)%tab_ind(:,iG),nDind(id)%indD_OF_Dm1(iG) ; flush(out_unitp)

      END IF
    END DO
  END DO
  IF (debug) THEN
    write(out_unitp,*) '======================================='
    write(out_unitp,*) 'id,tab_q',id,':',nDind(id)%tab_q
    write(out_unitp,*) 'id,MaxnD',id,':',nDind(id)%MaxnD
    write(out_unitp,*) 'id,i_TO_l',id,':',i_TO_l(:)
    CALL Write_TypeDInd(nDind(id))
    flush(out_unitp)
  END IF
  deallocate(i_TO_l)

END DO

END SUBROUTINE Set_nDInd_10order

SUBROUTINE Set_nDInd_01order_L(nDind,D,Lmin,Lmax)
IMPLICIT NONE

integer                         :: D,Lmin,Lmax
TYPE (TypeDInd), allocatable    :: nDind(:)

integer :: i,id,iGm1,iG,nG,l,ll,lll,ndimGm1
logical :: test

logical, parameter :: debug=.FALSE.
!logical, parameter :: debug=.TRUE.


IF (.NOT. allocated(nDind)) allocate(nDind(0:D+1))

!-------------------------------------------
!-------------------------------------------
!-------------------------------------------
!-------------------------------------------
id=0
CALL alloc_TypeDInd(nDind(id),ndim=1,MaxnD=1)
nDind(id)%tab_ind(1,1)   = 0
nDind(id)%indD_OF_Dm1(1) = 1
nDind(id)%tab_q(1)       = 0

id=1
CALL alloc_TypeDInd(nDind(id),ndim=1,MaxnD=1)
nDind(id)%tab_ind(1,1)   = 0
nDind(id)%indD_OF_Dm1(1) = 1
nDind(id)%tab_q(1)       = 0


DO id=2,D+1
  ! first the number of points

  nG = 0
  DO iGm1=1,nDind(id-1)%MaxnD
    ll = nDind(id-1)%i_TO_l(iGm1)
    DO l=0,Lmax
      lll=ll+l
      test = lll <= Lmax
      IF (id==D+1) test = test .AND. lll >= Lmin
      IF (test) nG = nG + 1
    END DO
  END DO
  CALL alloc_TypeDInd(nDind(id),ndim=id-1,MaxnD=nG)
  nDind(id)%tab_q(:) = (/ (i,i=1,id-1) /)


  iG = 0
  ndimGm1 = nDind(id-1)%ndim

  DO iGm1=1,nDind(id-1)%MaxnD
    ll = nDind(id-1)%i_TO_l(iGm1)
    DO l=0,Lmax
      lll=ll+l
      test = lll <= Lmax
      IF (id==D+1) test = test .AND. lll >= Lmin
      IF (test) THEN
        iG = iG + 1
        nDind(id)%tab_ind(1:ndimGm1,iG)      = nDind(id-1)%tab_ind(:,iGm1)
        nDind(id)%tab_ind(nDind(id)%ndim,iG) = l
        nDind(id)%indD_OF_Dm1(iG)            = iGm1
      END IF
    END DO
  END DO
  IF (debug) THEN
    write(out_unitp,*) '======================================='
    write(out_unitp,*) 'id,tab_q ',id,':',nDind(id)%tab_q
    write(out_unitp,*) 'id,MaxnD ',id,':',nDind(id)%MaxnD
    CALL Write_TypeDInd(nDind(id))
    flush(out_unitp)
  END IF


END DO

END SUBROUTINE Set_nDInd_01order_L
SUBROUTINE Set_nDInd_10order_L(nDind,D,Lmin,Lmax)
IMPLICIT NONE

integer        :: D,Lmin,Lmax
TYPE(TypeDInd), allocatable :: nDind(:)

integer :: i,id,iGp1,iG,nG,l,ll,lll,ndimGp1
logical :: test

logical, parameter :: debug=.FALSE.
!logical, parameter :: debug=.TRUE.

IF (.NOT. allocated(nDind)) allocate(nDind(0:D+1))
!-------------------------------------------
!-------------------------------------------
!-------------------------------------------
!-------------------------------------------
id=D+1
CALL alloc_TypeDInd(nDind(id),ndim=1,MaxnD=1)
nDind(id)%tab_ind(1,1)   = 0
nDind(id)%indD_OF_Dm1(1) = 1
nDind(id)%tab_q(1)       = 0

id=D
nDind(id)%ndim  = 1
nDind(id)%MaxnD = 1
CALL alloc_TypeDInd(nDind(id),ndim=1,MaxnD=1)
nDind(id)%tab_ind(1,1)   = 0
nDind(id)%indD_OF_Dm1(1) = 1
nDind(id)%tab_q(1)       = 0


DO id=D-1,0,-1

  ! first the number of points
  nG = 0
  DO iGp1=1,nDind(id+1)%MaxnD
    ll = sum(nDind(id+1)%tab_ind(:,iGp1))
    DO l=0,Lmax
      lll=ll+l
      test = lll <= Lmax
      IF (id==0) test = test .AND. lll >= Lmin
      IF (test) nG = nG + 1
    END DO
  END DO
  CALL alloc_TypeDInd(nDind(id),ndim=D-id,MaxnD=nG)
  nDind(id)%tab_q(:) = (/ (i,i=id+1,D) /)


  iG      = 0
  ndimGp1 = nDind(id+1)%ndim
  DO iGp1=1,nDind(id+1)%MaxnD
    ll = sum(nDind(id+1)%tab_ind(:,iGp1))

    DO l=0,Lmax
      lll=ll+l
      test = lll <= Lmax
      IF (id==0) test = test .AND. lll >= Lmin
      IF (test) THEN
        iG = iG + 1
        IF (id < D-1) nDind(id)%tab_ind(2:1+ndimGp1,iG) = nDind(id+1)%tab_ind(:,iGp1)
        nDind(id)%tab_ind(1,iG) = l

        nDind(id)%indD_OF_Dm1(iG)    = iGp1
        !write(out_unitp,*) 'id,iG,l(:)',id,iG,nDind(id)%tab_ind(:,iG),nDind(id)%indD_OF_Dm1(iG) ; flush(out_unitp)

      END IF
    END DO
  END DO
  IF (debug) THEN
    write(out_unitp,*) '======================================='
    write(out_unitp,*) 'id,tab_q',id,':',nDind(id)%tab_q
    write(out_unitp,*) 'id,MaxnD',id,':',nDind(id)%MaxnD
    CALL Write_TypeDInd(nDind(id))
    flush(out_unitp)
  END IF

END DO

END SUBROUTINE Set_nDInd_10order_L

SUBROUTINE dealloc_nDInd(nDind)
IMPLICIT NONE

TYPE(TypeDInd), allocatable :: nDind(:)

integer :: id


IF (allocated(nDind)) THEN

  DO id=lbound(nDind,dim=1),ubound(nDind,dim=1)
    CALL dealloc_TypeDInd(nDind(id))
  END DO

  deallocate(nDind)
END IF


END SUBROUTINE dealloc_nDInd

SUBROUTINE nDInd2TOnDInd1(nDInd1,nDInd2)
IMPLICIT NONE

TYPE(TypeDInd), allocatable, intent(inout) :: nDInd1(:)
TYPE(TypeDInd), allocatable, intent(in)    :: nDInd2(:)

integer :: i,li,ui

CALL dealloc_nDInd(nDind1)
IF (allocated(nDInd2)) THEN
  li = lbound(nDInd2,dim=1)
  ui = ubound(nDInd2,dim=1)
  allocate(nDind1(li:ui))
  DO i=li,ui
    nDind1(i) = nDind2(i)
  END DO
END IF

END SUBROUTINE nDInd2TOnDInd1

SUBROUTINE Write_Tab_nDInd(Tab_nDInd)
IMPLICIT NONE

TYPE(TypeDInd), allocatable :: Tab_nDInd(:)

integer :: i

write(out_unitp,*) 'BEGINNING Write_Tab_nDInd'

DO i=lbound(Tab_nDInd,dim=1),ubound(Tab_nDInd,dim=1)
  write(out_unitp,*) 'index:',i
  CALL Write_TypeDInd(Tab_nDInd(i))
END DO
write(out_unitp,*) 'END Write_Tab_nDInd'

END SUBROUTINE Write_Tab_nDInd

SUBROUTINE InD_TO_tabi(InD,D,tabn,tabi)
IMPLICIT NONE

integer          :: D,InD
integer          :: tabn(D),tabi(D)

integer          :: II,id,NN


II = InD-1
!DO id=D,1,-1
DO id=1,D
  tabi(id) = mod(II,tabn(id))
  II       = II/tabn(id)
END DO
tabi(:) = tabi(:) + 1

CALL tabi_TO_InD(II,D,tabn,tabi)


IF (II /= InD) STOP 'II /= InD'
!write(out_unitp,*) 'InD,tabn',InD,tabn
!write(out_unitp,*) 'InD,tabi',II,tabi

END SUBROUTINE InD_TO_tabi
SUBROUTINE tabi_TO_InD(InD,D,tabn,tabi)
IMPLICIT NONE

integer          :: D,InD
integer          :: tabn(D),tabi(D)

integer          :: II,id,NN


InD = 1
!DO id=1,D
DO id=D,1,-1
  InD = tabi(id) + tabn(id)*(InD-1)
END DO

!write(out_unitp,*) 'InD,tabn',InD,tabn
!write(out_unitp,*) 'InD,tabi',InD,tabi

END SUBROUTINE tabi_TO_InD

END MODULE mod_module_DInd
