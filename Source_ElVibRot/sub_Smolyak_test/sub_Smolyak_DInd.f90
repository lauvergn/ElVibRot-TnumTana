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
MODULE mod_Smolyak_DInd
USE mod_system
IMPLICIT NONE

TYPE TypeDInd
  integer :: ndim  = 0
  integer :: MaxnD = -1
  integer, allocatable :: tab_ind(:,:)
  integer, allocatable :: indD_OF_Dm1(:)
  integer, allocatable :: i_TO_l(:)  ! give the l value for i (usefull when i /= l+1)
  integer, allocatable :: tab_q(:)   ! size: ndim
END TYPE TypeDInd

INTERFACE assignment(=)
  module procedure TypeDInd2_TO_TypeDInd1
END INTERFACE

CONTAINS

SUBROUTINE alloc_TypeDInd(DInd,ndim,MaxnD)
USE mod_system
IMPLICIT NONE

integer        :: ndim,MaxnD
TYPE(TypeDInd) :: DInd


CALL dealloc_TypeDInd(DInd)

DInd%ndim  = ndim
DInd%MaxnD = MaxnD
allocate(DInd%tab_ind(ndim,MaxnD))
allocate(DInd%i_TO_l(MaxnD))  ! give the l value for i (usefull when i /= l+1)
DInd%i_TO_l(:) = 0
allocate(DInd%indD_OF_Dm1(MaxnD))
allocate(DInd%tab_q(ndim))

END SUBROUTINE alloc_TypeDInd
SUBROUTINE dealloc_TypeDInd(DInd)
USE mod_system
IMPLICIT NONE

TYPE(TypeDInd) :: DInd


DInd%ndim  = 0
DInd%MaxnD = 0
IF (allocated(DInd%tab_ind))     deallocate(DInd%tab_ind)
IF (allocated(DInd%i_TO_l))      deallocate(DInd%i_TO_l)
IF (allocated(DInd%indD_OF_Dm1)) deallocate(DInd%indD_OF_Dm1)
IF (allocated(DInd%tab_q))       deallocate(DInd%tab_q)

END SUBROUTINE dealloc_TypeDInd
SUBROUTINE TypeDInd2_TO_TypeDInd1(DInd1,DInd2)
USE mod_system
IMPLICIT NONE

TYPE(TypeDInd),    intent(in)    :: DInd2
TYPE(TypeDInd),    intent(inout) :: DInd1


CALL dealloc_TypeDInd(DInd1)

DInd1%ndim        = DInd1%ndim
DInd1%MaxnD       = DInd2%MaxnD

DInd1%tab_ind     = DInd2%tab_ind
DInd1%indD_OF_Dm1 = DInd2%indD_OF_Dm1

DInd1%i_TO_l      = DInd2%i_TO_l

DInd1%tab_q       = DInd2%tab_q

END SUBROUTINE TypeDInd2_TO_TypeDInd1
SUBROUTINE Write_TypeDInd(DInd)
USE mod_system
IMPLICIT NONE

TYPE(TypeDInd) :: DInd

integer :: I

write(6,*) 'BEGINNING Write_TypeDind'
write(6,*) 'ndim',DInd%ndim
write(6,*) 'MaxnD',DInd%MaxnD
write(6,*) 'tab_q(:)',DInd%tab_q(:)
write(6,*) 'I,L,indD_OF_Dm1,ind(:)'
DO I=1,DInd%MaxnD
  write(6,*) I,DInd%i_TO_l(I),DInd%indD_OF_Dm1(I),':',DInd%tab_ind(:,I)
END DO
write(6,*) 'END Write_TypeDind'
flush(6)

END SUBROUTINE Write_TypeDInd
SUBROUTINE Set_Smolyak_nDInd(SnDind,D,Lmin,Lmax)
USE mod_system
IMPLICIT NONE

integer        :: D,Lmin,Lmax
TYPE(TypeDInd) :: SnDind


TYPE(TypeDInd), allocatable :: nDind(:)


CALL Set_nDInd_10order(nDind,D,Lmin,Lmax)

SnDind = nDind(0)

SnDind%ndim = D

CALL dealloc_nDInd(nDind)

END SUBROUTINE Set_Smolyak_nDInd
SUBROUTINE Set_nDInd(nDind,D,Lmin,Lmax)
USE mod_system
IMPLICIT NONE

integer        :: D,Lmin,Lmax
TYPE(TypeDInd), allocatable :: nDind(:)

CALL Set_nDInd_01order(nDind,D,Lmin,Lmax,0)

END SUBROUTINE Set_nDInd
SUBROUTINE Set_nDInd_01order(nDind,D,Lmin,Lmax,type_l_TO_n)
USE mod_system
IMPLICIT NONE

integer        :: D,Lmin,Lmax,type_l_TO_n
TYPE(TypeDInd), allocatable :: nDind(:)

integer :: i,id,iGm1,iG,nG,l,ll,lll,ndimGm1,n
integer, allocatable :: i_TO_l(:)
logical :: test


IF (.NOT. allocated(nDind)) allocate(nDind(0:D+1))
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
  n = l_TO_n(Lmax,type_l_TO_n)
  allocate(i_TO_l(n))
  DO l=Lmax,0,-1
    i = l_TO_n(l,type_l_TO_n)
    i_TO_l(1:i) = l
  END DO

  nG = 0
  DO iGm1=1,nDind(id-1)%MaxnD
    ll = nDind(id-1)%i_TO_l(iGm1)
    DO i=1,n
      l  = i_TO_l(i)
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
        !write(6,*) 'id,iG,l(:)',id,iG,nDind(id)%tab_ind(:,iG),nDind(id)%indD_OF_Dm1(iG) ; flush(6)
      END IF
    END DO
  END DO
  !write(6,*) '======================================='
  !write(6,*) 'id,tab_q ',id,':',nDind(id)%tab_q
  !write(6,*) 'id,MaxnD ',id,':',nDind(id)%MaxnD
  !write(6,*) 'id,i_TO_l',id,':',i_TO_l(:)
  deallocate(i_TO_l)
  !CALL Write_TypeDInd(nDind(id))
  !flush(6)

END DO


END SUBROUTINE Set_nDInd_01order
SUBROUTINE Set_nDInd_10order(nDind,D,Lmin,Lmax)
USE mod_system
IMPLICIT NONE

integer        :: D,Lmin,Lmax
TYPE(TypeDInd), allocatable :: nDind(:)

integer :: i,id,iGp1,iG,nG,l,ll,lll,ndimGp1
logical :: test


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
        !write(6,*) 'id,iG,l(:)',id,iG,nDind(id)%tab_ind(:,iG),nDind(id)%indD_OF_Dm1(iG) ; flush(6)

      END IF
    END DO
  END DO
  !write(6,*) '======================================='
  !write(6,*) 'id,tab_q',id,':',nDind(id)%tab_q
  !write(6,*) 'id,MaxnD',id,':',nDind(id)%MaxnD
  !CALL Write_TypeDInd(nDind(id))
  flush(6)

END DO

END SUBROUTINE Set_nDInd_10order

SUBROUTINE dealloc_nDInd(nDind)
USE mod_system
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
SUBROUTINE Write_nDInd(nDind)
USE mod_system
IMPLICIT NONE

TYPE(TypeDInd), allocatable :: nDind(:)

integer :: i

write(6,*) 'BEGINNING Write_nDInd'

DO i=lbound(nDind,dim=1),ubound(nDind,dim=1)
  write(6,*) 'index:',i
  CALL Write_TypeDInd(nDind(i))
END DO
write(6,*) 'END Write_nDInd'

END SUBROUTINE Write_nDInd

SUBROUTINE InD_TO_tabi(InD,D,tabn,tabi)
USE mod_system
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
!write(6,*) 'InD,tabn',InD,tabn
!write(6,*) 'InD,tabi',II,tabi

END SUBROUTINE InD_TO_tabi
SUBROUTINE tabi_TO_InD(InD,D,tabn,tabi)
USE mod_system
IMPLICIT NONE

integer          :: D,InD
integer          :: tabn(D),tabi(D)

integer          :: II,id,NN


InD = 1
!DO id=1,D
DO id=D,1,-1
  InD = tabi(id) + tabn(id)*(InD-1)
END DO

!write(6,*) 'InD,tabn',InD,tabn
!write(6,*) 'InD,tabi',InD,tabi

END SUBROUTINE tabi_TO_InD

FUNCTION l_TO_n(l,typ,B)
USE mod_system
IMPLICIT NONE
integer :: l_TO_n

integer, intent(in), optional   :: B

integer, intent(in)             :: l,typ

SELECT CASE (typ)
CASE (0)
  l_TO_n = l+1
CASE (1)
  IF (present(B)) THEN
    l_TO_n = B*l+1
    !l_TO_n = 2+B
  ELSE
    l_TO_n = l+1
  END IF
CASE (2)
  l_TO_n = 1*l+1
CASE DEFAULT
  l_TO_n = l+1
END SELECT

IF (l < 0) l_TO_n = 0

END FUNCTION l_TO_n


END MODULE mod_Smolyak_DInd
