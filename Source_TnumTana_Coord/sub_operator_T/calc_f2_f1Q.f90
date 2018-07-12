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

      real (kind=Rkind) :: Tdef2_tot(mole%nb_var,mole%nb_var)
      real (kind=Rkind) :: Tdef1_tot(mole%nb_var)


!-------------------------------------------------------------------------

      real (kind=Rkind), parameter :: auTOcm_inv = 219474.63144319772_Rkind
      real (kind=Rkind), parameter :: inv_Name   = 1822.888485541_Rkind

      !real (kind=Rkind), parameter :: mH   = 1837.1526464003414_Rkind ! mH ! Tnum
      real (kind=Rkind), parameter :: mH   = 1.00800_Rkind * inv_Name ! mH ! Bacic

      real (kind=Rkind), parameter :: mH2  = TWO*mH ! mH2

      integer       :: i,iQdyn

      !real (kind=Rkind), parameter :: BH2 = 59.322_Rkind / auTOcm_inv ! from Valdez, PCCP, 2011, 13, 2935 and Bacic paper
      !real (kind=Rkind), parameter :: R=1.41897242230202301767_Rkind !  to recover the BH2 of Bacic
      real (kind=Rkind) :: BH2,R
      

      !real (kind=Rkind), parameter :: BH2 = 59.32266_Rkind / auTOcm_inv ! from Bacic program

      real (kind=Rkind) :: th,phi,c,s

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
      rho        = ZERO
      Tcor2(:,:) = ZERO
      Tcor1(:)   = ZERO
      Trot(:,:)  = ZERO

      Tdef2_tot(:,:) = ZERO
      Tdef1_tot(:)   = ZERO


      DO i=1,mole%nb_var
        Tdef2_tot(i,i) = -HALF/mH2
      END DO


      th  = Qdyn0(8)
      phi = Qdyn0(9)
      R   = Qdyn0(7)+Qdyn0(10)

      s = sin(th)
      c = cos(th)

      rho = s
      !Tdef2_tot(9,9) = -HALF/(Iner*s*s)
      !Tdef2_tot(8,8) = -HALF/Iner
      !Tdef1_tot(8)   = -HALF/Iner * c/s


      BH2 = ONE/(mH*R**2)
      !write(6,*) 'R/2',Qdyn0(7),Qdyn0(10)
      !write(6,*) 'BH2',BH2*auTOcm_inv
      !STOP
      Tdef2_tot(9,9) = -BH2/(s*s)
      Tdef2_tot(8,8) = -BH2
      Tdef1_tot(8)   = -BH2 * c/s

      DO i=1,mole%nb_act
        iQdyn = mole%liste_QactTOQsym(i)
        Tdef2(i,i) = Tdef2_tot(iQdyn,iQdyn)
        Tdef1(i)   = Tdef1_tot(iQdyn)
      END DO


      !IF (mole%nb_inact2n > 0) THEN
      !  DO i=3,2+mole%nb_inact2n
      !    Tdef2(i,i) = -HALF
      !  END DO
      !END IF

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
