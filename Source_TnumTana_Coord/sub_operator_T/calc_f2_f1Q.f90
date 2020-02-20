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
      SUBROUTINE calc_f2_f1Q_ana(Qdyn0,                                 &
                                 Tdef2,Tdef1,vep,rho,                   &
                                 Tcor2,Tcor1,Trot,                      &
                                 para_Tnum,mole)
      USE mod_system
      USE mod_Tnum
      IMPLICIT NONE

!----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole
      TYPE (Tnum)    :: para_Tnum

      real (kind=Rkind) ::  Qdyn0(mole%nb_var)

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

      real (kind=Rkind), parameter :: G11 = 0.00007981_Rkind
      real (kind=Rkind), parameter :: G22 = 0.00002599_Rkind
      real (kind=Rkind), parameter :: G33 = 0.00040375_Rkind
      real (kind=Rkind), parameter :: G23 = 0.00000000_Rkind

      integer :: i


!-------------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      !-----------------------------------------------------------
       IF (debug .OR. para_Tnum%WriteT) THEN
         write(out_unitp,*) 'BEGINNING calc_f2_f1Q_ana'
         write(out_unitp,*) 'ndimG',mole%ndimG
         write(out_unitp,*) 'WriteCC',mole%WriteCC
         write(out_unitp,*) 'Qdyn0',Qdyn0
         IF (debug) THEN
           write(out_unitp,*)
           CALL Write_mole(mole)
           write(out_unitp,*)
         END IF
         write(out_unitp,*)
        write(out_unitp,*) 'num_GG,num_g',para_Tnum%num_GG,para_Tnum%num_g
         write(out_unitp,*) 'num_x,nrho',para_Tnum%num_x,para_Tnum%nrho
         write(out_unitp,*) 'JJ',para_Tnum%JJ
         write(out_unitp,*)
       END IF
!-----------------------------------------------------------
      Tdef2(:,:) = ZERO
      Tdef1(:)   = ZERO
      vep        = ZERO
      rho        = ONE
      Tcor2(:,:) = ZERO
      Tcor1(:)   = ZERO
      Trot(:,:)  = ZERO

      Tdef2(1,1) = -HALF * G11
      Tdef2(2,2) = -HALF * G22
      Tdef2(3,3) = -HALF * G33
      Tdef2(2,3) = -G23
      Tdef2(3,2) = -G23


!-----------------------------------------------------------
      IF (debug .OR. para_Tnum%WriteT) THEN

        CALL Write_f2f1vep(Tdef2,Tdef1,vep,rho,mole%nb_act)
        IF (para_Tnum%JJ .GT. 0) CALL Write_TcorTrot(Tcor2,Tcor1,Trot,   &
                                           mole%nb_act)
        write(out_unitp,*) 'END calc_f2_f1Q_ana'
      END IF
!-----------------------------------------------------------

      end subroutine calc_f2_f1Q_ana
