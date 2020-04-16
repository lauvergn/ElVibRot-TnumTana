Program test
IMPLICIT NONE

integer, parameter :: Rkind=8

integer, parameter :: nb_Q0=6
real(kind=Rkind)   :: th,Q0(nb_Q0)
integer            :: i,nb_pts

integer                        :: nb,nq           ! numbers of basis functions and grid points
integer                        :: nb_vec          ! number of eigenvectors/eigenvalues
real (kind=Rkind), allocatable :: EigenVal(:)     ! eigenvalues.               EigenVal(i)
real (kind=Rkind), allocatable :: EigenVecB(:,:)  ! eigenvectors on the basis. EigenVecB(:,i)
real (kind=Rkind), allocatable :: EigenVecG(:,:)  ! eigenvectors on the grid.  EigenVecG(:,i)
real (kind=Rkind), allocatable :: RhoWeight(:)    ! rho(Q).Weight(Q), on the grid points.


Q0(:) = [0.900000_8,3.187000_8, 2.179000_8,0.000000_8,3.141593_8,0.000000_8]

CALL init_EVR()

CALL get_nb_nq(nb,nq)
allocate(EigenVal(nb))
allocate(EigenVecB(nb,nb))
allocate(EigenVecG(nq,nb))
allocate(RhoWeight(nq))
write(out_unitp,*) 'END init_EVR'

CALL levels_EVR(EigenVal,EigenVecB,EigenVecG,RhoWeight,nb,nq,nb_vec)
write(out_unitp,*) 'nb_vec',nb_vec
write(out_unitp,*) 'EigenVal(:)',EigenVal(1:nb_vec)

write(out_unitp,*) 'modify Q0'
nb_pts = 10
DO i=1,2*nb_pts-1
  th = -1.d0+real(i,kind=8)/real(nb_pts,kind=8)
  Q0(:) = [th,3.187_8, 2.179_8,0._8,3.141593_8,0._8]
  CALL Modify_TnumRefGeom_Q0(Q0,nb_Q0,.TRUE.)
  CALL levels_EVR(EigenVal,EigenVecB,EigenVecG,RhoWeight,nb,nq,nb_vec)
  write(out_unitp,*) 'EigenVal(:)',EigenVal(1:nb_vec)
END DO
CALL finalyze_EVR()
END
