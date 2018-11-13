PROGRAM CurviRPH
USE mod_system
USE CurviRPH_mod
implicit NONE


  TYPE (CurviRPH_type) :: CurviRPH_param

  real(kind=Rkind), allocatable :: Qpath(:),Q21(:),Grad(:),Hess(:,:)

  integer :: i,iq,jq,nb_Qpath,nb_Q21

  nb_Qpath = 1
  nb_Q21   = 2

  !CALL init_CurviRPH(CurviRPH_param,nb_Qpath,nb_Q21)

  allocate(Qpath(nb_Qpath))
  allocate(Q21(nb_Q21))
  allocate(Grad(nb_Q21))
  allocate(Hess(nb_Q21,nb_Q21))

  Qpath = 0.5_Rkind
  CALL get_CurviRPH(Qpath,CurviRPH_param,Q21,Grad,Hess)

  write(6,*) 'points:',i,Qpath(:)
  write(6,*) 'Q21',Q21(:)
  IF (CurviRPH_param%gradient) write(6,*) 'Grad',Grad(:)

  write(6,*) 'Hess'
  CALL Write_Mat(Hess,6,5)

END PROGRAM CurviRPH

!PROGRAM CurviRPH
!USE mod_system
!implicit NONE
!
!integer :: nb_pts,ndim,nb_dev
!character (len=256) :: name_dum
!real(kind=Rkind), allocatable :: Q(:,:),g(:,:),hess(:,:,:)
!real(kind=Rkind), allocatable :: aQ(:,:),ag(:,:),ahess(:,:,:)
!
!real(kind=Rkind), allocatable :: Qpath(:),fQ(:,:),fQ_inv(:,:),fQpath(:)
!logical :: gradient = .FALSE.
!
!
!integer :: i,j,iq,jq,IOerr
!
!read(5,*) nb_pts,ndim
!nb_dev = nb_pts
!write(6,*) 'nb_pts,nb_dev,ndim',nb_pts,nb_dev,ndim
!flush(6)
!
!allocate(Q(ndim,nb_pts))
!allocate(g(ndim,nb_pts))
!allocate(hess(ndim,ndim,nb_pts))
!
!allocate(aQ(nb_pts,ndim))
!allocate(ag(nb_pts,ndim))
!allocate(ahess(nb_pts,ndim,ndim))
!
!allocate(Qpath(nb_pts))
!allocate(fQ(nb_dev,nb_pts))
!allocate(fQ_inv(nb_pts,nb_dev))
!allocate(fQpath(nb_dev))
!
!
!DO i=1,nb_pts
!  read(5,*) Qpath(i)
!  write(6,*) 'Qpath',Qpath(i)
!  DO j=1,nb_dev
!    fQ(j,i) = Qpath(i)**(j-1)
!  END DO
!  !read geometry
!  read(5,*) Q(:,i)
!
!  !DO iq=1,ndim
!  !  read(5,*) name_dum,Q(iq,i)
!  !  write(6,*) 'Q',i,Q(iq,i)
!  !END DO
!  !read gradient
!  IF (gradient) THEN
!    DO iq=1,ndim
!      read(5,*) name_dum,g(iq,i)
!    END DO
!  END IF
!  !read hessian
!  CALL Read_Mat(hess(:,:,i),5,5,IOerr)
!  write(6,*) 'IOerr',IOerr
!  CALL Write_Mat(hess(:,:,i),6,5)
!
!END DO
!!CALL Write_Mat(fQ,6,5)
!
!CALL inv_m1_TO_m2(fQ,fQ_inv,nb_pts,0,ZERO)
!
!!CALL Write_Mat(fQ_inv,6,5)
!
!!for the fit of Q
!DO iq=1,ndim
!  aq(:,iq) = matmul(Q(iq,:),fQ_inv(:,:))
!  write(6,*) 'a(:)',iq,aq(:,iq)
!END DO
!
!!for the fit of g
!IF (gradient) THEN
!  DO iq=1,ndim
!    ag(:,iq) = matmul(g(iq,:),fQ_inv(:,:))
!    write(6,*) 'a(:)',iq,ag(:,iq)
!  END DO
!END IF
!
!!for the fit of hess
!DO iq=1,ndim
!DO jq=1,ndim
!  ahess(:,jq,iq) = matmul(hess(jq,iq,:),fQ_inv(:,:))
!  write(6,*) 'a(:)',iq,jq,ahess(:,jq,iq)
!END DO
!END DO
!
!
!DO i=1,nb_pts
!  write(6,*) 'points:',i,Qpath(i)
!  fQpath(:) = fQ(:,i)
!
!  DO iq=1,ndim
!    write(6,*) 'Q',iq,Q(iq,i),Q(iq,i)-dot_product(fQpath,aq(:,iq))
!  END DO
!
!  DO iq=1,ndim
!  DO jq=1,ndim
!    write(6,*) 'Q',jq,iq,hess(jq,iq,i),hess(jq,iq,i)-dot_product(fQpath,ahess(:,jq,iq))
!  END DO
!  END DO
!
!END DO
!
!END PROGRAM CurviRPH

