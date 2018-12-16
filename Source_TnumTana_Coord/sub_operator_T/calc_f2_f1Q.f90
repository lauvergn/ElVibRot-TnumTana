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
      SUBROUTINE calc_f2_f1Q_ana(Qsym0,                                 &
                                 Tdef2,Tdef1,vep,rho,                   &
                                 Tcor2,Tcor1,Trot,                      &
                                 para_Tnum,mole)
      USE mod_system
      USE mod_Tnum
      IMPLICIT NONE

!----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole
      TYPE (Tnum)    :: para_Tnum

      real (kind=Rkind) ::  Qsym0(mole%nb_var)

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
!            see num_H,num_A,num_x
!
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
      real (kind=Rkind) :: Tdef2(mole%nb_act,mole%nb_act)
      real (kind=Rkind) :: Tdef1(mole%nb_act)
      real (kind=Rkind) :: vep
      real (kind=Rkind) :: Tcor2(mole%nb_act,3),Tcor1(3)
      real (kind=Rkind) :: Trot(3,3)

      real (kind=Rkind) :: rho

      real (kind=Rkind) :: Tdef2_tot(mole%nb_var,mole%nb_var)
      real (kind=Rkind) :: Tdef1_tot(mole%nb_var)


!-------------------------------------------------------------------------

      real (kind=Rkind), parameter :: auTOcm_inv = 219474.63144319772_Rkind
      real (kind=Rkind), parameter :: inv_Name   = 1822.888485541_Rkind

      !real (kind=Rkind), parameter :: mH   = 1837.1526464003414_Rkind ! mH ! Tnum
       real (kind=Rkind), parameter :: mD   = 3671.4829394591770_Rkind ! mD ! Tnum
      real (kind=Rkind), parameter :: mH   = 1.00800_Rkind * inv_Name ! mH ! Bacic

      real (kind=Rkind), parameter :: mH2  = TWO*mD ! mH2

      integer       :: i,iQdyn

      real (kind=Rkind) :: BH2

      real (kind=Rkind) :: R,th,phi,c,s

!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      !-----------------------------------------------------------
       IF (debug .OR. para_Tnum%WriteT) THEN
         write(out_unitp,*) 'BEGINNING calc_f2_f1Q_ana'
         write(out_unitp,*) 'ndimG',mole%ndimG
         write(out_unitp,*) 'WriteCC',mole%WriteCC
         write(out_unitp,*) 'Qsym0',Qsym0
         IF (debug) THEN
           write(out_unitp,*)
           CALL Write_mole(mole)
           write(out_unitp,*)
         END IF
       END IF
!-----------------------------------------------------------

      Tdef2(:,:) = ZERO
      Tdef1(:)   = ZERO
      vep        = ZERO
      rho        = ZERO
      Tcor2(:,:) = ZERO
      Tcor1(:)   = ZERO
      Trot(:,:)  = ZERO

      Tdef2_tot(:,:) = ZERO
      Tdef1_tot(:)   = ZERO


      DO i=1,mole%nb_var
        Tdef2_tot(i,i) = -HALF
      END DO


      R   = Qsym0(1)

      th  = Qsym0(9)
      phi = Qsym0(10)

      BH2 = ONE/(mD*R**2)

      s = sin(th)
      c = cos(th)

      rho = s

      Tdef2_tot(1,1)   = -HALF/(mD/TWO)

      Tdef2_tot(6,6)   = -HALF/(mD*TWO)
      Tdef2_tot(7,7)   = -HALF/(mD*TWO)
      Tdef2_tot(8,8)   = -HALF/(mD*TWO)

      Tdef2_tot(9,9)   = -BH2
      Tdef2_tot(10,10) = -BH2/(s*s)

      Tdef1_tot(9)     = -BH2 * c/s

      DO i=1,mole%nb_act
        iQdyn = mole%liste_QactTOQsym(i)
        Tdef2(i,i) = Tdef2_tot(iQdyn,iQdyn)
        Tdef1(i)   = Tdef1_tot(iQdyn)
      END DO


      DO i=1,3
        Trot(i,i) = -HALF
      END DO

!-----------------------------------------------------------
      IF (debug .OR. para_Tnum%WriteT) THEN

        CALL Write_f2f1vep(Tdef2,Tdef1,vep,rho,mole%nb_act)
        IF (para_Tnum%JJ .GT. 0) CALL Write_TcorTrot(Tcor2,Tcor1,Trot, &
                                           mole%nb_act)
        write(out_unitp,*) 'END calc_f2_f1Q_ana'
      END IF
!-----------------------------------------------------------


      end subroutine calc_f2_f1Q_ana
