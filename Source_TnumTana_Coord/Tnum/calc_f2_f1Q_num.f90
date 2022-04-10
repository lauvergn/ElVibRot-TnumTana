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
MODULE mod_f2f2Vep
  USE mod_system
  USE mod_dnSVM
  USE mod_Tnum
  USE mod_paramQ
  USE mod_dnGG_dng
  USE mod_dnDetGG_dnDetg
  USE mod_dnRho
  IMPLICIT NONE

  INTERFACE calc3_f2_f1Q_num
    MODULE PROCEDURE calc3_f2_f1Q_num_CoordType,calc3_f2_f1Q_num_zmatrix
  END INTERFACE

  PRIVATE
  PUBLIC :: calc3_f2_f1Q_num,calc3_f2_f1Q_numTay0Qinact2n

  CONTAINS

!======================================================================
!
!      Calculation of Tdef Tcor Trot at Q
!      Tdef = Tdef2 * d2./dQ1dQ2 + Tdef1 * d./dQ1 + vep
!      Tcor = i * (Tcor2 * d./dQ1*Jx  + Tcor1 * Jx)
!      Trot = Trot  * Jx*Jy
!
!      Calculation of rho at Q
!      dT = rho*dQ
!
!
!======================================================================
      SUBROUTINE calc3_f2_f1Q_num_zmatrix(Qact,Tdef2,Tdef1,vep,rho,      &
                                          Tcor2,Tcor1,Trot,para_Tnum,mole)
      USE mod_system
      USE mod_dnSVM
      use mod_ActiveTransfo,    only: qact_to_qdyn_from_activetransfo
      USE mod_Tnum
      USE mod_paramQ
      USE mod_dnGG_dng
      USE mod_dnRho ! all
      USE mod_Tana_NumKEO
      IMPLICIT NONE

!----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole
      TYPE (Tnum)    :: para_Tnum

      real (kind=Rkind), intent(inout) :: Qact(mole%nb_var)
      real (kind=Rkind) :: Tdef2(mole%nb_act,mole%nb_act),Tdef1(mole%nb_act)
      real (kind=Rkind) :: vep,rho
      real (kind=Rkind) :: Tcor2(mole%nb_act,3),Tcor1(3),Trot(3,3)



        CALL calc3_f2_f1Q_num_CoordType(Qact,Tdef2,Tdef1,vep,rho,       &
                              Tcor2,Tcor1,Trot,para_Tnum,mole%CoordType)


      END SUBROUTINE calc3_f2_f1Q_num_zmatrix
      SUBROUTINE calc3_f2_f1Q_num_CoordType(Qact,Tdef2,Tdef1,vep,rho,   &
                                        Tcor2,Tcor1,Trot,para_Tnum,mole)
      USE mod_system
      USE mod_dnSVM
      use mod_ActiveTransfo,    only: qact_to_qdyn_from_activetransfo
      USE mod_Tnum
      USE mod_paramQ
      USE mod_dnGG_dng
      USE mod_dnRho ! all
      USE mod_Tana_NumKEO
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole
      TYPE (Tnum)      :: para_Tnum

      real (kind=Rkind), intent(inout) :: Qact(mole%nb_var)

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
!                             or adapted to the basis set)
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
      real (kind=Rkind) ::  rho

!-------------------------------------------------------------------------




!     - for memory ---------------------------------------------
      real (kind=Rkind) :: Qdyn(mole%nb_var)
      TYPE(Type_dnS)    :: dnJac,dnrho

      TYPE(Type_dnMat) :: dng,dnGG
      integer          :: calc_g,calc_GG
      !TYPE (zmatrix)   :: mole_zmatrix

!     ----------------------------------------------------------

      integer :: nderiv

      integer :: err,memory

!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'calc3_f2_f1Q_num'
!-----------------------------------------------------------
       IF (debug .OR. para_Tnum%WriteT) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'ndimG',mole%ndimG
         write(out_unitp,*) 'WriteCC',mole%WriteCC
         write(out_unitp,*) 'Qact',Qact
         !IF (debug) THEN
         !  write(out_unitp,*)
         !  CALL Write_CoordType(mole)
         !  write(out_unitp,*)
         !END IF
         write(out_unitp,*)
         write(out_unitp,*) 'num_GG,num_g',para_Tnum%num_GG,para_Tnum%num_g
         write(out_unitp,*) 'num_x,nrho',para_Tnum%num_x,para_Tnum%nrho
         write(out_unitp,*) 'JJ',para_Tnum%JJ
         write(out_unitp,*)
       END IF
!-----------------------------------------------------------
      IF (para_Tnum%Tana) THEN

        CALL get_Numf2f1vep_WITH_AnaKEO(para_Tnum%ExpandTWOxKEO,Qact,   &
                                        mole,Tdef2,Tdef1,vep,rho)

        IF (debug .OR. para_Tnum%WriteT) THEN
          write(out_unitp,*) ' f2,f1,vep with Tana'
          CALL Write_f2f1vep(Tdef2,Tdef1,vep,rho,mole%nb_act)
          !IF (para_Tnum%JJ > 0 .AND. .NOT. mole%Without_Rot)              &
          !      CALL Write_TcorTrot(Tcor2,Tcor1,Trot,mole%nb_act)

          write(out_unitp,*) 'END ',name_sub
        END IF
        RETURN
      END IF


      IF (mole%nb_Qtransfo == -1 .OR. para_Tnum%f2f1_ana) THEN
        CALL Qact_TO_Qdyn_FROM_ActiveTransfo(Qact,Qdyn,mole%ActiveTransfo)

        CALL calc_f2_f1Q_ana(Qdyn,                                      &
                             Tdef2,Tdef1,vep,rho,                       &
                             Tcor2,Tcor1,Trot,                          &
                             para_Tnum,mole)
        RETURN
      END IF


      IF (para_Tnum%vep_type == 0) THEN
        nderiv = 1
      ELSE
        nderiv = 2
      END IF

      CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,nderiv)
      CALL alloc_dnSVM(dng,mole%ndimG,mole%ndimG,mole%nb_act,nderiv)

      CALL get_dng_dnGG(Qact,para_Tnum,mole,dng,dnGG,vep=vep,nderiv=nderiv)
      ! write(out_unitp,*) 'dng'
      ! CALL write_dnSVM(dng)
      ! write(out_unitp,*) 'dnGG'
      ! CALL write_dnSVM(dnGG)
      ! write(out_unitp,*) 'vep',vep

      !----- For dnrho -------------------------------------------
      nderiv = 1
      ! for dnrho (because, we need dnrho%d1 in sub_H1def and sub_Tcor1
      CALL alloc_dnSVM(dnrho,dng%nb_var_deriv,nderiv)
      CALL alloc_dnSVM(dnJac,dng%nb_var_deriv,nderiv)

      IF (para_Tnum%nrho == 0 .AND. .NOT. para_Tnum%Gcte) THEN
        ! dnJac
        CALL sub3_dndetGG(dnJac,dnGG,nderiv,                            &
                          mole%masses,mole%Mtot_inv,mole%ncart)
        !write(out_unitp,*) 'dnJac'
        !CALL write_dnS(dnJac)
      ELSE
        CALL sub_ZERO_TO_dnS(dnJac)
      END IF
      CALL sub3_dnrho(dnrho,dnJac,Qact,mole,                            &
                      nderiv,para_Tnum%num_x,para_Tnum%stepT,           &
                      para_Tnum%nrho)
      !write(out_unitp,*) 'dnrho'
      !CALL write_dnS(dnrho)
!-----------------------------------------------------------

      rho = dnrho%d0
      !write(out_unitp,*) 'rho :',rho
!-----------------------------------------------------------

!==================================================================
!==================================================================
!
!      Hamiltonien
!
!==================================================================
!==================================================================

!     ------------------------------------------------------------
!     -- 2d derivative terms -------------------------------------
      CALL sub_H2def(Tdef2,dnGG%d0,mole%ndimG,mole%nb_act)

!     ------------------------------------------------------------

!     ------------------------------------------------------------
!     -- 1st derivative terms ------------------------------------
      IF (.NOT. para_Tnum%Gcte) THEN
        CALL sub_H1def(Tdef1,dnGG%d0,dnGG%d1,dnrho%d1,mole%ndimG,mole%nb_act)
      ELSE
        Tdef1(:) = ZERO
      END IF
!     ------------------------------------------------------------

!     ------------------------------------------------------------

      IF ( .NOT. mole%Without_Rot) THEN
        CALL sub_Tcor2(Tcor2,dnGG%d0,mole%ndimG,mole%nb_act)
        CALL sub_Trot(Trot,dnGG%d0,mole%ndimG,mole%nb_act)
        CALL sub_Tcor1(Tcor1,dnGG%d0,dnGG%d1,dnrho%d1,                  &
                       mole%ndimG,mole%nb_act)
      END IF

      CALL dealloc_dnSVM(dnJac)
      CALL dealloc_dnSVM(dnrho)
      CALL dealloc_dnSVM(dnGG)
      CALL dealloc_dnSVM(dng)

!-----------------------------------------------------------
      IF (debug .OR. para_Tnum%WriteT) THEN

        CALL Write_f2f1vep(Tdef2,Tdef1,vep,rho,mole%nb_act)
        IF (para_Tnum%JJ > 0 .AND. .NOT. mole%Without_Rot)              &
                CALL Write_TcorTrot(Tcor2,Tcor1,Trot,mole%nb_act)

        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE calc3_f2_f1Q_num_CoordType
      SUBROUTINE calc3_f2_f1Q_numTay0Qinact2n(Qact,dnQinact2n,          &
                                              Tdef2,Tdef1,vep,rho,      &
                                              Tcor2,Tcor1,Trot,         &
                                              para_Tnum,mole)
      USE mod_system
      USE mod_dnSVM
      USE mod_Tnum
      USE mod_paramQ
      USE mod_dnGG_dng
      USE mod_dnRho ! all
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType)  :: mole
      TYPE (Tnum)       :: para_Tnum

      real (kind=Rkind) :: Tdef2(mole%nb_act,mole%nb_act)
      real (kind=Rkind) :: Tdef1(mole%nb_act)
      real (kind=Rkind) :: vep
      real (kind=Rkind) :: Tcor2(mole%nb_act,3),Tcor1(3)
      real (kind=Rkind) :: Trot(3,3)
      real (kind=Rkind) ::  rho

      real (kind=Rkind), intent(inout) :: Qact(mole%nb_var)
      TYPE(Type_dnVec)  :: dnQinact2n


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

!-------------------------------------------------------------------------




!     - for memory ---------------------------------------------
      TYPE(Type_dnS)    :: dnJac,dnrho
      TYPE(Type_dnMat)  :: dng,dnGG
      integer           :: calc_g,calc_GG
!     ----------------------------------------------------------

      integer :: nderiv
      integer :: i,j

      integer :: err,memory

!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
!      logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'calc3_f2_f1Q_numTay0Qinact2n'
!-----------------------------------------------------------
       IF (debug .OR. para_Tnum%WriteT) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'ndimG',mole%ndimG
         write(out_unitp,*) 'WriteCC',mole%WriteCC
         write(out_unitp,*) 'Qact',Qact
         IF (debug) THEN
           write(out_unitp,*)
           CALL Write_CoordType(mole)
           write(out_unitp,*)
         END IF
         write(out_unitp,*)
         write(out_unitp,*) 'num_GG,num_g',para_Tnum%num_GG,para_Tnum%num_g
         write(out_unitp,*) 'num_x,nrho',para_Tnum%num_x,para_Tnum%nrho
         write(out_unitp,*) 'JJ',para_Tnum%JJ
         write(out_unitp,*)
       END IF
!-----------------------------------------------------------
      IF (mole%nb_Qtransfo == -1) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'You cannot use Taylor expansion and "calc_f2_f1Q_ana"'
        STOP
      END IF


      IF (para_Tnum%vep_type == 0) THEN
        nderiv = 1
      ELSE
        nderiv = 2
      END IF

      CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,nderiv)
      CALL alloc_dnSVM(dng,mole%ndimG,mole%ndimG,mole%nb_act,nderiv)

      CALL get_dng_dnGG(Qact,para_Tnum,mole,dng,dnGG,nderiv=nderiv)

      ! transformation of dnGG%d1
      DO i=1,mole%nb_act1 ! first the derivatives with respect to the type1 coord.
      DO j=1,mole%nb_inact2n
         dnGG%d1(:,:,i) = dnGG%d1(:,:,i) +                              &
                        dnGG%d1(:,:,j+mole%nb_act1) * dnQinact2n%d1(j,i)
      END DO
      END DO
      DO j=1,mole%nb_inact2n ! first the derivatives with respect to the type2n coord.
         dnGG%d1(:,:,j+mole%nb_act1) = ZERO
      END DO

      CALL alloc_dnSVM(dnJac,dng%nb_var_deriv,nderiv)
      CALL alloc_dnSVM(dnrho,dng%nb_var_deriv,nderiv)

      IF ( (para_Tnum%vep_type /= 0 .OR. para_Tnum%nrho == 0) .AND. .NOT. para_Tnum%num_GG) THEN
        write(out_unitp,*) ' ERROR in name_sub'
        write(out_unitp,*) '  You have to use nrho=10 or 20 or vep_type=0 to use the taylor expansion'
        write(out_unitp,*) '   => without extrapotential term (vep)'
        write(out_unitp,*) ' Check your data !'
        STOP

        CALL sub3_dndetA(dnJac,dng,nderiv,                              &
                         mole%masses,mole%Mtot_inv,mole%ncart)

      ELSE ! nrho = 20 or 10 (jac0 is not used => jac0=1)
        CALL Set_ZERO_TO_dnSVM(dnJac)
        dnJac%d0 = ONE
      END IF
      !write(out_unitp,*) 'jac0 :'
      !CALL Write_dnS(dnJac)

!-----------------------------------------------------------

!-----------------------------------------------------------
!     f0,fi,Fij calculation
!-----------------------------------------------------------
      CALL sub3_dnrho(dnrho,dnJac,Qact,mole,                            &
                      nderiv,para_Tnum%num_x,para_Tnum%stepT,           &
                      para_Tnum%nrho)
!-----------------------------------------------------------
      rho = dnrho%d0
      !write(out_unitp,*) 'rho :'
      !CALL Write_dnS(dnrho)

!==================================================================
!==================================================================
!
!      Hamiltonien
!
!==================================================================
!==================================================================

!     ------------------------------------------------------------
!     -- 2d derivative terms -------------------------------------
      CALL sub_H2def(Tdef2,dnGG%d0,mole%ndimG,mole%nb_act)
!     ------------------------------------------------------------

!     ------------------------------------------------------------
!     -- 1st derivative terms ------------------------------------
      CALL sub_H1def(Tdef1,dnGG%d0,dnGG%d1,dnrho%d1,mole%ndimG,mole%nb_act)
!     ------------------------------------------------------------

!     ------------------------------------------------------------
!     -- extra potential term (is always zero for the taylor expansion)-
      vep = ZERO
!     ------------------------------------------------------------

      IF (.NOT. mole%Without_Rot) THEN
        CALL sub_Tcor2(Tcor2,dnGG%d0,mole%ndimG,mole%nb_act)
        CALL sub_Trot(Trot,dnGG%d0,mole%ndimG,mole%nb_act)
        CALL sub_Tcor1(Tcor1,dnGG%d0,dnGG%d1,dnrho%d1,mole%ndimG,mole%nb_act)
      END IF

      CALL dealloc_dnSVM(dnJac)
      CALL dealloc_dnSVM(dnrho)
      CALL dealloc_dnSVM(dnGG)
      CALL dealloc_dnSVM(dng)

!-----------------------------------------------------------
      IF (debug .OR. para_Tnum%WriteT) THEN

        CALL Write_f2f1vep(Tdef2,Tdef1,vep,rho,mole%nb_act)
        IF (para_Tnum%JJ > 0 .AND. .NOT. mole%Without_Rot)              &
                        CALL Write_TcorTrot(Tcor2,Tcor1,Trot,mole%nb_act)

        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE calc3_f2_f1Q_numTay0Qinact2n


      SUBROUTINE sub_H2def(H2ij,d0invA,ndimA,nb_act)
      USE mod_system
      IMPLICIT NONE

      integer :: ndimA,nb_act

      real (kind=Rkind) ::  d0invA(ndimA,ndimA)
      real (kind=Rkind) ::  H2ij(nb_act,nb_act)

      integer i,j

!     -- 2d derivative terms -------------------------------------
      DO i=1,nb_act
        H2ij(i,i)=-HALF*d0invA(i,i)
        DO j=i+1,nb_act
          H2ij(i,j)=-d0invA(i,j)
          H2ij(j,i)=H2ij(i,j)
        END DO
      END DO

      end subroutine sub_H2def
      SUBROUTINE sub_Tcor2(Tcor2,d0invA,ndimA,nb_act)
      USE mod_system
      IMPLICIT NONE

      integer ndimA,nb_act

      real (kind=Rkind) ::  d0invA(ndimA,ndimA)
      real (kind=Rkind) ::  Tcor2(nb_act,3)

      integer i,j

      DO i=1,nb_act
!       - terms in i * d./dQi*Jx ---------------------------------
        Tcor2(i,1) = - d0invA(i,nb_act+1)
!       - terms in i * d./dQi*Jy ---------------------------------
        Tcor2(i,2) = - d0invA(i,nb_act+2)
!       - terms in i * d./dQi*Jz ---------------------------------
        Tcor2(i,3) = - d0invA(i,nb_act+3)
      END DO

      end subroutine sub_Tcor2
      SUBROUTINE sub_Tcor1(Tcor1,d0invA,d1invA,fi,ndimA,nb_act)
      USE mod_system
      IMPLICIT NONE

       integer ndimA,nb_act
       real (kind=Rkind) ::  d0invA(ndimA,ndimA)
       real (kind=Rkind) ::  d1invA(ndimA,ndimA,nb_act)
       real (kind=Rkind) ::  fi(nb_act)
       real (kind=Rkind) ::  Tcor1(3)

       integer i

!     -- Jx Jy Jz terms ------------------------------------------
      !write(out_unitp,*) 'fi:',fi(:)
      Tcor1(:) = ZERO
      DO i=1,nb_act
        Tcor1(1) = Tcor1(1) - d1invA(i,nb_act+1,i) - fi(i)*d0invA(i,nb_act+1)
        Tcor1(2) = Tcor1(2) - d1invA(i,nb_act+2,i) - fi(i)*d0invA(i,nb_act+2)
        Tcor1(3) = Tcor1(3) - d1invA(i,nb_act+3,i) - fi(i)*d0invA(i,nb_act+3)
      END DO
      Tcor1(:) = HALF * Tcor1(:)
!     ------------------------------------------------------------


      end subroutine sub_Tcor1
      SUBROUTINE sub_Tcor1_2(Tcor1,d0invA,d1invA,fi,ndimA,nb_act)
      USE mod_system
      IMPLICIT NONE

       integer ndimA,nb_act
       real (kind=Rkind) ::  d0invA(ndimA,ndimA)
       real (kind=Rkind) ::  d1invA(ndimA,ndimA)
       real (kind=Rkind) ::  fi(nb_act)
       real (kind=Rkind) ::  Tcor1(3)

       integer i,j

!     -- Jx Jy Jz terms ------------------------------------------
      Tcor1(:) = ZERO
      DO i=1,nb_act
        Tcor1(1) = Tcor1(1) - d1invA(i,nb_act+1) - fi(i)*d0invA(i,nb_act+1)
        Tcor1(2) = Tcor1(2) - d1invA(i,nb_act+2) - fi(i)*d0invA(i,nb_act+2)
        Tcor1(3) = Tcor1(3) - d1invA(i,nb_act+3) - fi(i)*d0invA(i,nb_act+3)
      END DO
      Tcor1(:) = HALF * Tcor1(:)
!     ------------------------------------------------------------

      end subroutine sub_Tcor1_2



      SUBROUTINE sub_Trot(Trot,d0invA,ndimA,nb_act)
      USE mod_system
      IMPLICIT NONE

      integer ndimA,nb_act

      real (kind=Rkind) ::  d0invA(ndimA,ndimA)
      real (kind=Rkind) ::  Trot(3,3)

      integer i


      Trot(:,:) = HALF * d0invA(nb_act+1:nb_act+3,nb_act+1:nb_act+3)

      END SUBROUTINE sub_Trot



      SUBROUTINE sub_H1def_2(H1i,d0invA,d1invA,fi,ndimA,nb_act)
      USE mod_system
      IMPLICIT NONE

       integer ndimA,nb_act
       real (kind=Rkind) ::  d0invA(ndimA,ndimA)
       real (kind=Rkind) ::  d1invA(ndimA,ndimA)
       real (kind=Rkind) ::  fi(nb_act)
       real (kind=Rkind) ::  H1i(nb_act)

       integer i,j

!     -- 1st derivative terms -------------------------------------
      DO i=1,nb_act
        H1i(i) = ZERO
        DO j=1,nb_act
          H1i(i) = H1i(i) + d1invA(j,i) + d0invA(j,i)*fi(j)
        END DO
        H1i(i) = -HALF*H1i(i)

!       write(out_unitp,*) ' H1i',i,H1i(i)
!       write(out_unitp,*) 'fi',fi
!       write(out_unitp,*) 'd1invA(.,i)',(d1invA(j,i),j=1,nb_act)
!       write(out_unitp,*) 'd0invA(.,i)',(d0invA(j,i),j=1,nb_act)
      END DO
!     ------------------------------------------------------------

      end subroutine sub_H1def_2
      SUBROUTINE sub_H1def(H1i,d0invA,d1invA,fi,ndimA,nb_act)
      USE mod_system
      IMPLICIT NONE

       integer :: ndimA,nb_act
       real (kind=Rkind) ::  d0invA(ndimA,ndimA)
       real (kind=Rkind) ::  d1invA(ndimA,ndimA,nb_act)
       real (kind=Rkind) ::  fi(nb_act)
       real (kind=Rkind) ::  H1i(nb_act)

       integer :: i,j

!     -- 1st derivative terms -------------------------------------
        !write(out_unitp,*) 'fi:',fi(:)
        DO i=1,nb_act
          H1i(i) = ZERO
          DO j=1,nb_act
            H1i(i) = H1i(i) + d1invA(j,i,j) + d0invA(j,i)*fi(j)
          END DO
          H1i(i) = -HALF*H1i(i)
        END DO
!     ------------------------------------------------------------

      end subroutine sub_H1def


END MODULE mod_f2f2Vep
