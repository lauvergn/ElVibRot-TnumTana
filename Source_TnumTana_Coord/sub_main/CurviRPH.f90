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
!      Tnum is written David Lauvergnat [1]
!      Tana is written by Mamadou Ndong [1] and David Lauvergnat [1]
!         with contributions
!          Emil Lund klinting (coupling with MidasCpp) [3]'
!
![1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France
![3]: Department of Chemistry, Aarhus University, DK-8000 Aarhus C, Denmark
!
!===========================================================================
!===========================================================================
PROGRAM CurviRPH
use mod_system
use mod_dnSVM
use mod_Constant
! in the use mod_Coord_KEO, we have to use "only", because "calc_freq" is
!   a subroutine in mod_Coord_KEO and also a variable in the namelist.
use mod_Coord_KEO,  ONLY: CoordType,Tnum,Read_CoordType,              &
                          read_RefGeom,get_Qact0,sub_QactTOdnx,       &
                          get_d0g_d0GG, &
                          dealloc_CoordType
use mod_PrimOp
use ADdnSVM_m, only : dns_t,set_dnS,Write_dnS,get_d0, &
                      dnMat_t,Write_dnMat
implicit NONE

! - parameters for para_Tnum -----------------------
  TYPE (constant)  :: const_phys
  TYPE (CoordType) :: mole
  TYPE (Tnum)      :: para_Tnum
  TYPE (PrimOp_t)  :: PrimOp
  real (kind=Rkind), allocatable :: Qact(:)
  real (kind=Rkind), allocatable :: d0GG(:,:),d0g(:,:)
  real (kind=Rkind), allocatable :: d0GGNew(:,:),d0gNew(:,:)

  real (kind=Rkind), allocatable :: betaO(:),alphaON(:,:),alphaON0(:,:)

  real (kind=Rkind), allocatable :: d0Qop(:),d1Qop(:),d2Qop(:),hessNew(:,:)

  real (kind=Rkind), allocatable :: JacON(:,:)


  integer :: i,j,is
  real (kind=Rkind) :: norm2,s,x,Qact_QML(1)
  ! for QML
  integer :: ndim,nsurf,nb_Func,ndimFunc
  real(kind=Rkind), allocatable  :: d0Func(:)
  real(kind=Rkind), allocatable  :: d1Func(:,:)
  real(kind=Rkind), allocatable  :: d2Func(:,:,:)
  real(kind=Rkind), allocatable  :: d3Func(:,:,:,:)
  real(kind=Rkind), allocatable  :: V(:,:),G(:,:,:),H(:,:,:,:),Q_QML(:)

  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.

  !===========================================================================
  !===========================================================================
  ! Tnum data read (without flexible transformation)

  CALL versionEVRT(.TRUE.)
  print_level=2

  !-----------------------------------------------------------------
  !     - read the coordinate transformations :
  !     -   zmatrix, polysperical, bunch...
  !     ------------------------------------------------------------
  CALL Read_CoordType(mole,para_Tnum,const_phys)
  !     ------------------------------------------------------------
  !-----------------------------------------------------------------

  !-----------------------------------------------------------------
  !     - read coordinate values -----------------------------------
  !     ------------------------------------------------------------
  CALL read_RefGeom(mole,para_Tnum)
  !     ------------------------------------------------------------
  !-----------------------------------------------------------------
  !-----------------------------------------------------------------
  !     ---- TO finalize the coordinates (NM) and the KEO ----------
  !     ------------------------------------------------------------
  CALL Finalize_TnumTana_Coord_PrimOp(para_Tnum,mole,PrimOp)
  !===========================================================================

  !===========================================================================
  ! QML initialization
  ndim  = 0
  nsurf = 0
  !CALL sub_Init_Qmodel_Cart(ndim,nsurf,'H3_LSTH',.FALSE.,0) ! initialization
  CALL sub_Init_Qmodel(ndim,nsurf,'H3_LSTH',.FALSE.,0) ! initialization
  CALL get_Qmodel_nb_Func_ndimFunc(nb_Func,ndimFunc)
  write(out_unitp,*) 'ndimFunc,nb_Func',ndimFunc,nb_Func
  flush(out_unitp)

  allocate(d0Func(nb_Func))
  allocate(d1Func(ndimFunc,nb_Func))
  allocate(d2Func(ndimFunc,ndimFunc,nb_Func))
  allocate(d3Func(ndimFunc,ndimFunc,ndimFunc,nb_Func))

  ! QML Cartesian coord, potential, gradient and hessian
  allocate(Q_QML(ndim))
  allocate(V(nsurf,nsurf))
  allocate(G(nsurf,nsurf,ndim))
  allocate(H(nsurf,nsurf,ndim,ndim))
  !===========================================================================

  allocate(Qact(mole%nb_var))
  allocate(betaO(mole%nb_act))
  allocate(alphaON(mole%nb_act,mole%nb_act))
  allocate(alphaON0(mole%nb_act,mole%nb_act))
  allocate(d0GG(mole%ndimG,mole%ndimG))
  allocate(d0g(mole%ndimG,mole%ndimG))
  allocate(d0GGNew(mole%ndimG,mole%ndimG))
  allocate(d0gNew(mole%ndimG,mole%ndimG))

  allocate(d0Qop(mole%nb_act))
  allocate(d1Qop(mole%nb_act))
  allocate(d2Qop(mole%nb_act))
  allocate(hessNew(mole%nb_act,mole%nb_act))
  allocate(JacON(mole%ndimG,mole%ndimG))


  alphaON0(:,:) = ZERO

  s = -100._Rkind
  DO is=-40,40
    s = ONETENTH**1 * is

  !===========================================================================
  ! get Qop(s) with QML
  Qact_QML(1) = s
  d0Qop(:) = ZERO
  d1Qop(:) = ZERO

  CALL get_Qmodel_d0d1d2d3Func(d0Func,d1Func,d2Func,d3Func,Qact_QML,nb_Func,ndimFunc)
  d0Qop(1:2) = [d0Func(2:3)]
  d1Qop(1:2) = [d1Func(1,2:3)]
  d2Qop(1:2) = [d2Func(1,1,2:3)]

  IF (debug) write(out_unitp,*) 's,d0Qop,d1Qop',s,d0Qop,d1Qop
  flush(out_unitp)


  !Q_QML(:) = [ZERO,ZERO,-d0Func(2),  ZERO,ZERO,ZERO,  ZERO,ZERO,d0Func(3)] ! in Cartesian
  Q_QML(:) = [d0Qop(1:2),sum(d0Qop(1:2))]
  CALL sub_Qmodel_VGH(V,G,H,Q_QML)

  IF (debug)  CALL Write_Vec(G(1,1,:),out_unitp,5,name_info='grad')
  IF (debug)  CALL Write_Mat(H(1,1,:,:),out_unitp,5,name_info='hess')

  Qact(:) = [d0Qop(1:2),-0.999999_Rkind] ! in 3D Valence coordinates (to be changed)
  IF (debug)  write(out_unitp,*) 's,Qact',s,Qact(:)
  !===========================================================================

  !===========================================================================
  ! get the metric tensors
  CALL get_d0g_d0GG(Qact,para_Tnum,mole,d0g=d0g,d0GG=d0GG,def=.FALSE.)
  IF (debug) CALL Write_Mat(d0g,out_unitp,5,name_info='d0g')
  IF (debug) CALL Write_Mat(d0GG,out_unitp,5,name_info='d0GG')
  !===========================================================================

  !===========================================================================
  ! get the betaO(i) = sum_i' d0g(i,i')d1Qop(i')
  ! and the alphaON(:,:)
  CALL make_betaO(betaO,d0g,d1Qop)

  CALL make_alphaON(alphaON,betaO)
  !CALL make_alphaON_fit(s,alphaON)
  DO i=1,size(alphaON,dim=2)
    norm2 = dot_product(alphaON(:,i),alphaON(:,i))
    IF (sqrt(abs(norm2)) < ONETENTH**10) CYCLE

    norm2 = dot_product(alphaON0(:,i),alphaON0(:,i))
    IF (norm2 > ONETENTH**10 .AND. dot_product(alphaON(:,i),alphaON0(:,i)) < ZERO) THEN
      alphaON(:,i) = -alphaON(:,i)
    END IF
    alphaON0(:,i) = alphaON(:,i)
    !write(out_unitp,*) 's,alphaON(:,i)',s,i,alphaON(:,i)
  END DO



  IF (debug) CALL check_gnew(d0g,alphaON,betaO)

  CALL make_JacON(JacON,alphaON,d1Qop)

  !===========================================================================

  !===========================================================================
  ! get the metric tensors
  d0gNew = matmul(transpose(JacON),matmul(d0g,JacON))
  IF (debug) CALL Write_Mat(d0gNew,out_unitp,5,name_info='d0g new')

  CALL inv_m1_TO_m2(d0gNew,d0GGNew,mole%ndimG,0,ZERO)
  IF (debug) CALL Write_Mat(d0GGNew(1:2,1:2),out_unitp,5,name_info='d0GG new')

  ! get the new hessian (without the gradient contribution)
  IF (debug) write(out_unitp,*) 'gradNew',matmul(transpose(JacON(1:2,1:2)),G(1,1,1:2))
  hessNew = matmul(transpose(JacON(1:2,1:2)),matmul(H(1,1,1:2,1:2),JacON(1:2,1:2)))
  IF (debug) CALL Write_Mat(hessNew,out_unitp,5,name_info='hessNew')

  write(out_unitp,*) 's,d0Qop,G,V,freq',s,d0Qop,d0GGNew(1,1),d0GGNew(2,2),V,    &
                      hessNew(2,2),sqrt(hessNew(2,2)*d0GGNew(2,2))*219475._Rkind



  ! modification of alphaON, to take into account the new hessian and metric tensor
  DO i=1,size(alphaON,dim=2)
    norm2 = dot_product(alphaON(:,i),alphaON(:,i))
    IF (sqrt(abs(norm2)) < ONETENTH**10) CYCLE
    alphaON(:,i) = alphaON(:,i) / sqrt(sqrt(hessNew(2,2)/d0GGNew(2,2)))
    write(out_unitp,*) 's,alphaON(:,i)',s,i,alphaON(:,i)
  END DO

  IF (debug)  CALL check_gnew(d0g,alphaON,betaO)
  CALL make_JacON(JacON,alphaON,d1Qop)
  CALL Write_Mat(JacON,out_unitp,5,name_info='JacON')

  !===========================================================================
  ! get the metric tensors
  d0gNew = matmul(transpose(JacON),matmul(d0g,JacON))
  CALL Write_Mat(d0gNew,out_unitp,5,name_info='d0g new')

  CALL inv_m1_TO_m2(d0gNew,d0GGNew,mole%ndimG,0,ZERO)
  CALL Write_Mat(d0GGNew(1:2,1:2),out_unitp,5,name_info='d0GG new')

  ! get the new hessian (without the gradient contribution)
  IF (debug) write(out_unitp,*) 'gradNew',matmul(transpose(JacON(1:2,1:2)),G(1,1,1:2))
  hessNew = matmul(transpose(JacON(1:2,1:2)),matmul(H(1,1,1:2,1:2),JacON(1:2,1:2)))
  IF (debug) CALL Write_Mat(hessNew,out_unitp,5,name_info='hessNew')

  write(out_unitp,*) 's,G(2,2),hess(2,2),freq',s,d0GGNew(2,2),hessNew(2,2),     &
                      sqrt(hessNew(2,2)*d0GGNew(2,2))

  END DO
CONTAINS
  SUBROUTINE make_JacON(JacON,alphaON,d1Qop)
    use mod_system
    implicit NONE

    real(kind=Rkind), intent(inout) :: JacON(:,:)
    real(kind=Rkind), intent(in)    :: alphaON(:,:)
    real(kind=Rkind), intent(in)    :: d1Qop(:)

    integer :: i,j,n,ndimG
    logical, parameter :: debug=.FALSE.

    ndimG = size(JacON,dim=1)
    n = size(d1Qop)
    if (size(alphaON,dim=1) /= n .OR. size(alphaON,dim=2) /= n)                       &
                                    STOP 'ERROR in make1_alpha: inconsistent size'

    CALL mat_id(JacON,ndimG,ndimG)

    JacON(1:n,1) = d1Qop ! along s
    j = 2
    DO i=1,n ! new
      IF (sqrt(abs(dot_product(alphaON(:,i),alphaON(:,i)))) < ONETENTH**10) CYCLE
      JacON(1:n,j) = alphaON(:,i)
      j = j + 1
    END DO
    IF (debug) CALL Write_Mat(JacON,out_unitp,5,name_info='JacON')

  END SUBROUTINE make_JacON

  SUBROUTINE check_gnew(d0g,alphaON,betaO)
    use mod_system
    implicit NONE

    real(kind=Rkind), intent(in) :: d0g(:,:)
    real(kind=Rkind), intent(in) :: alphaON(:,:)
    real(kind=Rkind), intent(in) :: betaO(:)

    integer :: i,j,n

    n = size(betaO)
    if (size(alphaON,dim=1) /= n .OR. size(alphaON,dim=2) /= n)                       &
                                    STOP 'ERROR in check_gnew: inconsistent size'

    DO j=1,n
      write(6,*) 'g1j new',j,dot_product(alphaON(:,j),betaO)
    END DO

  END SUBROUTINE check_gnew

  SUBROUTINE make_alphaON(alphaON,betaO)
    use mod_system
    implicit NONE

    real(kind=Rkind), intent(inout) :: alphaON(:,:)
    real(kind=Rkind), intent(in)    :: betaO(:)

    integer :: i,j,n
    real(kind=Rkind) :: x,norm2
    logical, parameter :: debug=.FALSE.

    n = size(betaO)
    if (size(alphaON,dim=1) /= n .OR. size(alphaON,dim=2) /= n)                 &
                                STOP 'ERROR in make_alphaON: inconsistent size'

    ! Schmidt ortho
    CALL mat_id(alphaON,n,n)
    !CALL Write_Mat(alphaON,out_unitp,5,name_info='alphaON')
    DO i=1,n ! new coord
      ! first ortho against betaO
      x = dot_product(alphaON(:,i),betaO)
      IF (debug) write(out_unitp,*) i,x
      alphaON(:,i) = alphaON(:,i) - x*betaO
      x = dot_product(alphaON(:,i),betaO)
      IF (debug) write(out_unitp,*) i,x

      ! then ortho against the previous vectors
      DO j=1,i-1
        norm2 = dot_product(alphaON(:,j),alphaON(:,j))
        IF (sqrt(abs(norm2)) < ONETENTH**10) CYCLE

        x = dot_product(alphaON(:,i),alphaON(:,j))/norm2
        IF (debug) write(out_unitp,*) i,j,x
        alphaON(:,i) = alphaON(:,i) - x*alphaON(:,j)
        x = dot_product(alphaON(:,i),alphaON(:,j))/norm2
        IF (debug) write(out_unitp,*) i,j,x
      END DO

      norm2 = dot_product(alphaON(:,i),alphaON(:,i))

      IF (sqrt(abs(norm2)) < ONETENTH**10) THEN
        alphaON(:,i) = ZERO
      ELSE
        alphaON(:,i) = alphaON(:,i)/sqrt(norm2)
      END IF

      IF (debug) write(out_unitp,*) i,alphaON(:,i)
    END DO

    IF (debug) THEN
      DO i=1,n
        write(out_unitp,*) 'alphaON(:,i)',i,alphaON(:,i)
      END DO
    END IF

  END SUBROUTINE make_alphaON

  SUBROUTINE make_alphaON_fit(s,alphaON)
    use mod_system
    implicit NONE

    real(kind=Rkind), intent(inout) :: alphaON(:,:)
    real(kind=Rkind), intent(in)    :: s

    real(kind=Rkind) :: UNpoly_legendre ! function


    real(kind=Rkind) :: tx,th

    tx = tanh(0.8_Rkind*s)
    th = &
    (0.7853980685664477_Rkind)      * UNpoly_legendre(tx, 1,0) +  &
    (1.6552831655323903_Rkind)      * UNpoly_legendre(tx, 2,0) +  &
    (3.9331177710684384e-07_Rkind)  * UNpoly_legendre(tx, 3,0) +  &
    (-0.49565178041423397_Rkind)    * UNpoly_legendre(tx, 4,0) +  &
    (-3.190325029540911e-07_Rkind)  * UNpoly_legendre(tx, 5,0) +  &
    (0.10656909615661768_Rkind)     * UNpoly_legendre(tx, 6,0) +  &
    (1.8690480109475923e-07_Rkind)  * UNpoly_legendre(tx, 7,0) +  &
    (-0.0032496111357101162_Rkind)  * UNpoly_legendre(tx, 8,0) +  &
    (-9.3180214164917e-08_Rkind)    * UNpoly_legendre(tx, 9,0) +  &
    (-0.02826693610529469_Rkind)    * UNpoly_legendre(tx,10,0) +  &
    (3.377042785834256e-08_Rkind)   * UNpoly_legendre(tx,11,0) +  &
    (0.025134143523729502_Rkind)    * UNpoly_legendre(tx,12,0) +  &
    (-4.9310394603674515e-09_Rkind) * UNpoly_legendre(tx,13,0) +  &
    (-0.015786871028081358_Rkind)   * UNpoly_legendre(tx,14,0) +  &
    (-5.808070293114866e-09_Rkind)  * UNpoly_legendre(tx,15,0) +  &
    (0.00725533736619045_Rkind)     * UNpoly_legendre(tx,16,0) +  &
    (7.018569477713918e-09_Rkind)   * UNpoly_legendre(tx,17,0) +  &
    (-0.0029716961471153826_Rkind)  * UNpoly_legendre(tx,18,0) +  &
    (-3.2559298856508953e-09_Rkind) * UNpoly_legendre(tx,19,0) +  &
    (0.0007270727526154396_Rkind)   * UNpoly_legendre(tx,20,0)

    alphaON(:,1) = [cos(th),sin(th)]
    alphaON(:,2) = ZERO

  END SUBROUTINE make_alphaON_fit

  SUBROUTINE make_betaO(betaO,d0g,d1Qop)
    use mod_system
    implicit NONE

    real(kind=Rkind), intent(in)    :: d0g(:,:),d1Qop(:)
    real(kind=Rkind), intent(inout) :: betaO(:)

    integer :: i,j,n
    real(kind=Rkind) :: x
    logical, parameter :: debug=.FALSE.

    n = size(betaO)
    if (size(d1Qop) /= n)  STOP 'ERROR in make_betaO: inconsistent size'


    DO i=1,n
      betaO(i) = dot_product(d0g(1:n,i),d1Qop)
    END DO

    IF (debug) write(out_unitp,*) 'betaO',betaO(:)
    betaO = betaO / sqrt(dot_product(betaO,betaO))
    IF (debug) write(out_unitp,*) 'Normalized betaO',betaO(:)


  END SUBROUTINE make_betaO

END PROGRAM CurviRPH
