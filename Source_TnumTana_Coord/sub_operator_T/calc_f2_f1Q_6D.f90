!======================================================================
!
!      Calculation of Tdef Tcor Trot at Q
!      Tdef = Tdef2 * d2./dQ1dQ2 + Tdef1 * d./dQ1 + vep
!      Tcor = Tcor2 * d./dQ1*Jx  + Tcor1 * Jx
!      Trot = Trot  * Jx*Jy
!
!      Calculation of rho at Q
!      dT = rho*dQ
!
!
!======================================================================
      SUBROUTINE calc_f2_f1Q_ana(Qdyn,                                 &
                                 Tdef2,Tdef1,vep,rho,                   &
                                 Tcor2,Tcor1,Trot,                      &
                                 para_Tnum,mole)
      USE mod_system
      USE mod_Tnum
      USE mod_file
      IMPLICIT NONE

!----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole
      TYPE (Tnum)    :: para_Tnum

      real (kind=Rkind) ::  Qdyn(mole%nb_var)

!-------------------------------------------------------------------------
      TYPE (param_file) :: file_Tnum
      integer :: ni ! unit of the file
      integer :: n0,n1,n2 ! number of term in the Taylor expansion
      integer, parameter :: max_var = 6
      real (kind=Rkind) ::  Qact_ref(max_var) ! G is expand around Qsym_ref
      real (kind=Rkind) ::  d0G(max_var,max_var)
      real (kind=Rkind) ::  d1G(max_var,max_var,max_var)
      real (kind=Rkind) ::  d2G(max_var,max_var,max_var,max_var)

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
!--- variables for Tnum --------------------------------------------------
!
!      Tdef = Tdef2(1,2) * d2./dQ1dQ2 + Tdef1(1) * d./dQ1 + vep
!      Tcor = Tcor2(1,x) * d./dQ1*Jx  + Tcor1(x) * Jx
!      Trot = Trot(x,y)  * Jx*Jy
!
!      Calculation of rho at Q
!      dT = rho*dQ
!
!
!      nrho: type of normalization
!            1 => Wilson (rho = 1.)
!           10 => Wilson (rho = 1.) (without vep (vep=0))
!            2 => clever choice ( Q=R =>1.; Q=val => sin(val); Q=cos(val) =>1. Q=dih =>1.)
!           20 => id 2 but without vep (vep=0)
!            0 => Euclidian (rho = Jac)
!
!     JJ: Total angular momentum
!
!
!     stepT: displacement for numerical calculation of Tnum
!            see num_GG,num_g,num_x
!
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
      real (kind=Rkind) :: Tdef2(mole%nb_act,mole%nb_act)
      real (kind=Rkind) :: Tdef1(mole%nb_act)
      real (kind=Rkind) :: vep
      real (kind=Rkind) :: Tcor2(mole%nb_act,3),Tcor1(3)
      real (kind=Rkind) :: Trot(3,3)

      real (kind=Rkind) :: rho

!-------------------------------------------------------------------------

      real (kind=Rkind) :: GG
      integer       :: i,j,k,l,ij
      real (kind=Rkind) ::  DQact(mole%nb_var)

      logical :: begin
      data begin/.true./
      SAVE begin,Qact_ref,n0,n1,n2,d0G,d1G,d2G


!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       IF (debug .OR. para_Tnum%WriteT) THEN
         write(out_unitp,*) 'BEGINNING calc_f2_f1Q_ana'
         write(out_unitp,*) 'ndimG',mole%ndimG
         write(out_unitp,*) 'WriteCC',mole%WriteCC
         write(out_unitp,*) 'Qdyn',Qdyn
         write(out_unitp,*)
         write(out_unitp,*) 'num_GG,num_g',para_Tnum%num_GG,para_Tnum%num_g
         write(out_unitp,*) 'num_x,nrho',para_Tnum%num_x,para_Tnum%nrho
         write(out_unitp,*) 'JJ',para_Tnum%JJ
         write(out_unitp,*)
       END IF
!-----------------------------------------------------------

      IF (max_var /= mole%nb_act) THEN
        write(out_unitp,*) ' ERROR in calc_f2_f1Q_ana'
        write(out_unitp,*) ' max_var /= mole%nb_act',max_var,mole%nb_act
        write(out_unitp,*) ' Change max_var in calc_f2_f1Q_ana and recompile'
        STOP
      END IF

      Tdef2(:,:) = ZERO
      Tdef1(:)   = ZERO
      vep        = ZERO
      rho        = ZERO
      Tcor2(:,:) = ZERO
      Tcor1(:)   = ZERO
      Trot(:,:)  = ZERO


!---------------------------------------------------------------
!     initialisation la premiere fois
      IF (begin) THEN
         begin = .FALSE.
         d0G(:,:)     = ZERO
         d1G(:,:,:)   = ZERO
         d2G(:,:,:,:) = ZERO
         file_Tnum%name='Tnum.op'
         CALL file_open(file_Tnum,ni)
         IF (mole%nb_act > max_var) THEN
           write(out_unitp,*) ' ERROR in calc_f2_f1Q_ana'
           write(out_unitp,*) ' mole%nb_act > max_var',mole%nb_act,max_var
           STOP
         END IF
         read(ni,*) Qact_ref(1:mole%nb_act)
         read(ni,*) n0
         write(out_unitp,*) 'Tnum.op,n0',n0
         DO ij=1,n0
           read(ni,*) i,j,GG
           d0G(i,j) = GG
           d0G(j,i) = GG
         END DO
         read(ni,*) n1
         write(out_unitp,*) 'n1',n1
         DO ij=1,n1
           read(ni,*) i,j,k,GG
           d1G(i,j,k) = GG
           d1G(j,i,k) = GG
         END DO
         read(ni,*) n2
         write(out_unitp,*) 'n2',n2
         DO ij=1,n2
           read(ni,*) i,j,k,l,GG
           d2G(i,j,k,l) = GG
           d2G(j,i,k,l) = GG
           d2G(i,j,l,k) = GG
           d2G(j,i,l,k) = GG
         END DO
         close(ni)
      END IF
!---------------------------------------------------------------

      rho = ONE

      DO i=1,mole%nb_act
        DQact(i) = Qdyn(mole%ActiveTransfo%list_QactTOQdyn(i)) - Qact_ref(i)
      END DO


!     Tdef2 : second order
      DO i=1,mole%nb_act
      DO j=1,mole%nb_act
        Tdef2(i,j) = d0G(i,j)
        DO k=1,mole%nb_act
          Tdef2(i,j) = Tdef2(i,j) + d1G(i,j,k)*DQact(k)
        END DO
        DO k=1,mole%nb_act
        DO l=1,mole%nb_act
          Tdef2(i,j) = Tdef2(i,j) + HALF*d2G(i,j,k,l)*DQact(k)*DQact(l)
        END DO
        END DO
      END DO
      END DO

      DO i=1,mole%nb_act
        Tdef2(i,i) = HALF * Tdef2(i,i)
      END DO
      Tdef2(:,:) = -Tdef2(:,:)

!     Tdef1 : second order
      DO i=1,mole%nb_act
        Tdef1(i) = ZERO
        DO j=1,mole%nb_act
          Tdef1(i) = Tdef1(i) + d1G(i,j,j)
        END DO
        DO j=1,mole%nb_act
        DO k=1,mole%nb_act
          Tdef1(i) = Tdef1(i) + d2G(i,j,j,k)*DQact(k)
        END DO
        END DO
      END DO
      Tdef1(:) = -HALF * Tdef1(:)


      DO i=1,3
        Trot(i,i) = -HALF
      END DO

!-----------------------------------------------------------
      IF (debug .OR. para_Tnum%WriteT) THEN


        CALL Write_f2f1vep(Tdef2,Tdef1,vep,rho,mole%nb_act)
        IF (para_Tnum%JJ .GT. 0) CALL Write_TcorTrot(Tcor2,Tcor1,Trot,   &
                                           mole%nb_act)
        write(out_unitp,*) 'END calc_f2_f1Q_ana'
      END IF
!-----------------------------------------------------------


      end subroutine calc_f2_f1Q_ana
