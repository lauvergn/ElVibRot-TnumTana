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
MODULE mod_dnGG_dng
  use mod_system
  use mod_dnSVM
  use mod_paramQ,  only: sub_QactTOdnMWx, Write_dnx, analyze_dnx
  use mod_Tnum,    only: tnum, CoordType,Write_CoordType,dealloc_CoordType
  USE mod_dnRho ! all

  IMPLICIT NONE

  PRIVATE

  logical, parameter :: new_vep = .FALSE.

  INTERFACE get_dng_dnGG
    MODULE PROCEDURE get_dng_dnGG_CoordType
  END INTERFACE
  INTERFACE get_d0GG
    MODULE PROCEDURE get_d0GG_CoordType
  END INTERFACE

  PUBLIC :: new_vep,sub_vep_new,sub_vep,Calc_vep_rho_from_dnGG
  PUBLIC :: get_d0g_d0GG,get_d0GG,get_dnGG_vep,get_dng_dnGG
  PUBLIC :: Set_dnVepTaylor

  CONTAINS
!======================================================================
!
!      Calculation of d0g,d1g,d2g ...
!      Also vep if needed
!
!      Calculation of rho and vep at Qdyn
!      dT = rho*dQ
!
!
!======================================================================
  SUBROUTINE get_dnGG_vep(Qact,para_Tnum,mole,dnGG,vep,rho,nderiv)
  USE mod_dnDetGG_dnDetg, only : sub3_dnDetGG
  IMPLICIT NONE

    !-----------------------------------------------------------------
    real (kind=Rkind), intent(in)               :: Qact(:)
    TYPE(Type_dnMat),  intent(inout)            :: dnGG
    real (kind=Rkind), intent(inout)            :: vep,rho
    integer,           intent(in)               :: nderiv
    !----- for the CoordType and Tnum --------------------------------
    TYPE (CoordType),  intent(in)               :: mole
    TYPE (Tnum),       intent(in)               :: para_Tnum
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! local variables
    integer           :: nderiv_loc
    !-----------------------------------------------------------------


    !--- for debuging --------------------------------------------------
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    character (len=*), parameter :: name_sub = 'get_dnGG_vep'
    !---------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*)
      write(out_unitp,*) 'ndimG',mole%ndimG
      write(out_unitp,*)
      write(out_unitp,*) 'Qact',Qact
      write(out_unitp,*)
      CALL Write_CoordType(mole)
      write(out_unitp,*)
      write(out_unitp,*)
      write(out_unitp,*) 'num_GG,num_g',para_Tnum%num_GG,para_Tnum%num_g
      write(out_unitp,*) 'num_x,nrho',para_Tnum%num_x,para_Tnum%nrho
      write(out_unitp,*) 'JJ',para_Tnum%JJ
      write(out_unitp,*)
    END IF
    !---------------------------------------------------------

    nderiv_loc = nderiv
    IF (nderiv_loc < 0) nderiv_loc = 0
    IF (nderiv_loc > 2) nderiv_loc = 2

    CALL get_dng_dnGG_CoordType(Qact,para_Tnum,mole,dnGG=dnGG,vep=vep,nderiv=nderiv_loc)

    !---------------------------------------------------------
    IF (debug .OR. para_Tnum%WriteT) THEN
      write(out_unitp,*) 'vep',vep
      write(out_unitp,*)
      write(out_unitp,*) 'dnGG'
      CALL Write_dnSVM(dnGG)
      write(out_unitp,*) 'END ',name_sub
    END IF
    !-----------------------------------------------------------


  END SUBROUTINE get_dnGG_vep
!======================================================================
!
!      Calculation of d0g,d1g,d2g ...
!
!      at Qdyn
!
!======================================================================
  SUBROUTINE get_dng_dnGG_CoordType(Qact,para_Tnum,mole,dng,dnGG,vep,nderiv)
  IMPLICIT NONE

    !-----------------------------------------------------------------
    real (kind=Rkind), intent(in)               :: Qact(:)
    TYPE(Type_dnMat),  intent(inout), optional  :: dng,dnGG
    real (kind=Rkind), intent(inout), optional  :: vep
    integer,           intent(in)               :: nderiv
    !----- for the CoordType and Tnum --------------------------------
    TYPE (CoordType),  intent(in)               :: mole
    TYPE (Tnum),       intent(in)               :: para_Tnum
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! local variables
    integer           :: i,j
    real (kind=Rkind) :: vep_loc,rho
    logical           :: vep_done
    !-----------------------------------------------------------------


    !----- for debuging --------------------------------------------------
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    character (len=*), parameter :: name_sub = 'get_dng_dnGG_CoordType'
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*)
      write(out_unitp,*) 'ndimG',mole%ndimG
      write(out_unitp,*)
      write(out_unitp,*) 'Qact',Qact
      !write(out_unitp,*)
      !CALL Write_CoordType(mole)
      write(out_unitp,*)
      write(out_unitp,*)
      write(out_unitp,*) 'num_GG,num_g',para_Tnum%num_GG,para_Tnum%num_g
      write(out_unitp,*) 'num_x,nrho',para_Tnum%num_x,para_Tnum%nrho
      write(out_unitp,*) 'JJ',para_Tnum%JJ
      write(out_unitp,*) 'GTaylor_Order',para_Tnum%GTaylor_Order
      write(out_unitp,*)
      write(out_unitp,*) 'present dnGG',present(dnGG)
      write(out_unitp,*) 'present dng ',present(dng)
      write(out_unitp,*) 'present vep ',present(vep)
      write(out_unitp,*)

      CALL flush_perso(out_unitp)
    END IF
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    IF (.NOT. present(dnGG) .AND. .NOT. present(dng)) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) ' dng and dnGG are not present !!'
        write(out_unitp,*) '  check the fortran!!'
        STOP
    END IF
    !-----------------------------------------------------------------

    vep_done = .FALSE.
    vep_loc  = ZERO

    IF (para_Tnum%GTaylor_Order > -1) THEN
      !-----------------------------------------------------------------
      ! with constant metric tensor
      IF (present(dnGG)) THEN
        IF (present(dng)) THEN
          CALL get_dng_dnGG_WITH_GTaylor(Qact,para_Tnum,mole,dng=dng,dnGG=dnGG) ! it doesn't work yet with dng
        ELSE
          CALL get_dng_dnGG_WITH_GTaylor(Qact,para_Tnum,mole,dnGG=dnGG)
        END IF
      ELSE
        IF (present(dng)) THEN
          CALL get_dng_dnGG_WITH_GTaylor(Qact,para_Tnum,mole,dng=dng)  ! it doesn't work yet with dng
        ELSE
          ! nothing to do !!
        END IF
      END IF
      !-----------------------------------------------------------------
    ELSE IF (para_Tnum%Gcte .AND. associated(para_Tnum%Gref)) THEN
      !-----------------------------------------------------------------
      ! with constant metric tensor
      IF (present(dnGG)) THEN
        IF (present(dng)) THEN
          CALL get_dng_dnGG_WITH_Gcte(Qact,para_Tnum,mole,dng=dng,dnGG=dnGG)
        ELSE
          CALL get_dng_dnGG_WITH_Gcte(Qact,para_Tnum,mole,dnGG=dnGG)
        END IF
      ELSE
        IF (present(dng)) THEN
          CALL get_dng_dnGG_WITH_Gcte(Qact,para_Tnum,mole,dng=dng)
        ELSE
          ! nothing to do !!
        END IF
      END IF
      !-----------------------------------------------------------------

   ELSE IF (mole%nb_Qtransfo == -1 .OR. para_Tnum%f2f1_ana) THEN
      !-----------------------------------------------------------------
      ! with anlytical expression of f2, f1
      IF (present(dnGG)) THEN
        IF (present(dng)) THEN
          CALL get_dng_dnGG_WITH_f2f1_ana(Qact,para_Tnum,mole,          &
                                          dng=dng,dnGG=dnGG,vep=vep_loc,&
                                          nderiv=nderiv)
        ELSE
          CALL get_dng_dnGG_WITH_f2f1_ana(Qact,para_Tnum,mole,          &
                                    dnGG=dnGG,vep=vep_loc,nderiv=nderiv)
        END IF
      ELSE
        IF (present(dng)) THEN
          CALL get_dng_dnGG_WITH_f2f1_ana(Qact,para_Tnum,mole,          &
                                      dng=dng,vep=vep_loc,nderiv=nderiv)
        ELSE
          ! nothing to do !!
        END IF
      END IF
      !-----------------------------------------------------------------
      vep_done = .TRUE.
    ELSE IF (mole%nb_rigid100 > 0) THEN
      !-----------------------------------------------------------------
      ! with type100 coordinates
      IF (present(dnGG)) THEN
        IF (present(dng)) THEN
          IF (para_Tnum%vep_type == -100) THEN
            ! This version has a bug. We keep it to recover published results (HNO3)
            CALL get_dng_dnGG_WITH_type100_bug_dng(Qact,para_Tnum,mole, &
                                         dng=dng,dnGG=dnGG,             &
                                         vep=vep_loc,vep_done=vep_done, &
                                         nderiv=nderiv)
            write(out_unitp,*) 'WARNING vep_type=-100. vep',vep_loc
          ELSE
            CALL get_dng_dnGG_WITH_type100(Qact,para_Tnum,mole,         &
                                         dng=dng,dnGG=dnGG,             &
                                         vep=vep_loc,vep_done=vep_done, &
                                         nderiv=nderiv)
          END IF
        ELSE
          IF (para_Tnum%vep_type == -100) THEN
            ! This version has a bug. We keep it to recover published results (HNO3)
            CALL get_dng_dnGG_WITH_type100_bug_dng(Qact,para_Tnum,mole, &
                                         dnGG=dnGG,                     &
                                         vep=vep_loc,vep_done=vep_done, &
                                         nderiv=nderiv)
            write(out_unitp,*) 'WARNING vep_type=-100. vep',vep_loc
          ELSE
            CALL get_dng_dnGG_WITH_type100(Qact,para_Tnum,mole,         &
                                         dnGG=dnGG,                     &
                                         vep=vep_loc,vep_done=vep_done, &
                                         nderiv=nderiv)
          END IF
        END IF
      ELSE
        IF (present(dng)) THEN
          CALL get_dng_dnGG_WITH_type100(Qact,para_Tnum,mole,           &
                                         dng=dng,                       &
                                         vep=vep_loc,vep_done=vep_done, &
                                         nderiv=nderiv)
        ELSE
          ! nothing to do !!
        END IF
      END IF
      !-----------------------------------------------------------------
      !write(6,*) 'vep_type,vep_done',para_Tnum%vep_type,vep_done
    ELSE
      !-----------------------------------------------------------------
      ! Normal computation ...
      IF (present(dng) .AND. present(dnGG)) THEN
        IF (para_Tnum%num_GG) THEN
          ! here we need to:
          !   (1) calculate dnGG numerically with get_dnGG.
          !   (2) inversion to get dng
          !  This is important for type200
          CALL get_dnGG(Qact,dnGG,nderiv,para_Tnum,mole)
          CALL INV_dnMat1_TO_dnMat2(dnGG,dng,nderiv)
        ELSE
          CALL get_dng(Qact,dng,nderiv,para_Tnum,mole)
          CALL INV_dnMat1_TO_dnMat2(dng,dnGG,nderiv)
        END IF
      ELSE IF (present(dnGG)) THEN
        CALL get_dnGG(Qact,dnGG,nderiv,para_Tnum,mole)
      ELSE IF (present(dng)) THEN
        CALL get_dng(Qact,dng,nderiv,para_Tnum,mole)
      END IF
      !-----------------------------------------------------------------

    END IF

    !-----------------------------------------------------------------
    ! To remove off-diagonal components
    IF (para_Tnum%Gdiago) THEN
      DO i=1,mole%nb_act
      DO j=1,mole%nb_act
        IF (i == j) CYCLE
        IF (present(dnGG)) THEN
          dnGG%d0(i,j) = ZERO
          IF (dnGG%nderiv > 0) dnGG%d1(i,j,:)   = ZERO
          IF (dnGG%nderiv > 1) dnGG%d2(i,j,:,:) = ZERO
        END IF
        IF (present(dng)) THEN
          dng%d0(i,j) = ZERO
          IF (dnGG%nderiv > 0) dng%d1(i,j,:)   = ZERO
          IF (dnGG%nderiv > 1) dng%d2(i,j,:,:) = ZERO
        END IF
      END DO
      END DO
    END IF
    !-----------------------------------------------------------------
    IF (present(vep)) THEN
      IF (vep_done) THEN
        vep = vep_loc
      ELSE
        IF (present(dnGG)) THEN
          vep = get_vep_rho(Qact,mole,para_Tnum,dnGG,rho)
        ELSE
          vep = get_vep_rho(Qact,mole,para_Tnum,rho=rho)
        END IF
        vep_done = .TRUE.
      END IF
    END IF

!-----------------------------------------------------------
      IF (debug) THEN
        IF (present(dnGG)) THEN
          write(out_unitp,*) 'dnGG'
          CALL Write_dnSVM(dnGG)
        END IF

        IF (present(dng)) THEN
          write(out_unitp,*) 'dng'
          CALL Write_dnSVM(dng)
        END IF
        IF (present(vep) .AND. vep_done) write(out_unitp,*) 'vep',vep

        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

  END SUBROUTINE get_dng_dnGG_CoordType

  SUBROUTINE get_dng_dnGG_WITH_Gcte(Qact,para_Tnum,mole,dng,dnGG)
      IMPLICIT NONE

      real (kind=Rkind), intent(in)               :: Qact(:)
      TYPE(Type_dnMat),  intent(inout), optional  :: dng,dnGG

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType),  intent(in)               :: mole
      TYPE (Tnum),       intent(in)               :: para_Tnum



!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'get_dng_dnGG_WITH_Gcte'
!-----------------------------------------------------------

      IF (associated(para_Tnum%Gref)) THEN
        IF (present(dnGG)) THEN
          CALL sub_ZERO_TO_dnMat(dnGG,nderiv=0)
          dnGG%d0(:,:) = para_Tnum%Gref(:,:)
        END IF
        IF (present(dng)) THEN
          CALL sub_ZERO_TO_dnMat(dng,nderiv=0)
          CALL inv_m1_TO_m2(para_Tnum%Gref,dng%d0,mole%ndimG,0,ZERO)
        END IF
      ELSE
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) ' You cannot call this subroutine when para_Tnum%Gref in not associated.'
        write(out_unitp,*) '  check the fortran!!'
        STOP 'ERROR in get_dng_dnGG_WITH_Gcte, para_Tnum%Gref is not associated.'
      END IF

  END SUBROUTINE get_dng_dnGG_WITH_Gcte

  SUBROUTINE get_dng_dnGG_WITH_GTaylor(Qact,para_Tnum,mole,dng,dnGG)
    IMPLICIT NONE

    real (kind=Rkind), intent(in)               :: Qact(:)
    TYPE(Type_dnMat),  intent(inout), optional  :: dng,dnGG

!----- for the CoordType and Tnum --------------------------------------
    TYPE (CoordType),  intent(in)               :: mole
    TYPE (Tnum),       intent(in)               :: para_Tnum


    real (kind=Rkind), allocatable    :: DQ(:)
    integer                           :: i,j

!----- for debuging --------------------------------------------------
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    character (len=*), parameter :: name_sub = 'get_dng_dnGG_WITH_GTaylor'
!-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*)
      write(out_unitp,*) 'ndimG',mole%ndimG
      write(out_unitp,*)
      write(out_unitp,*) 'Qact',Qact
      write(out_unitp,*) 'Qact0',mole%ActiveTransfo%Qact0
      write(out_unitp,*)
      write(out_unitp,*) 'GTaylor_Order',para_Tnum%GTaylor_Order
      write(out_unitp,*)
      write(out_unitp,*) 'present dnGG',present(dnGG)
      write(out_unitp,*) 'present dng ',present(dng)
      write(out_unitp,*)
      CALL flush_perso(out_unitp)
    END IF

    IF (.NOT. present(dnGG) .OR. present(dng)) THEN
      write(out_unitp,*) 'ERROR in ',name_sub
      write(out_unitp,*) ' This subroutine works ONLY with dnGG, and ...'
      write(out_unitp,*) ' dnGG is not present or dng is present'
      write(out_unitp,*) ' present(dnGG)',present(dnGG)
      write(out_unitp,*) ' present(dng) ',present(dng)
      write(out_unitp,*) '  check the fortran!!'
      STOP 'ERROR in get_dng_dnGG_WITH_GTaylor, Problem with dnGG or dng.'
    END IF
    IF (.NOT. present(dnGG) .OR. present(dng)) THEN
      write(out_unitp,*) 'ERROR in ',name_sub
      write(out_unitp,*) ' This subroutine works ONLY with dnGG, and ...'
      write(out_unitp,*) ' dnGG is not present or dng is present'
      write(out_unitp,*) ' present(dnGG)',present(dnGG)
      write(out_unitp,*) ' present(dng) ',present(dng)
      write(out_unitp,*) '  check the fortran!!'
      STOP 'ERROR in get_dng_dnGG_WITH_GTaylor, Problem with dnGG or dng.'
    END IF

    CALL check_alloc_dnMat(dnGG,'dnGG',name_sub)
    CALL check_alloc_dnMat(para_Tnum%dnGGref,'para_Tnum%dnGGref',name_sub)

    DQ = Qact - mole%ActiveTransfo%Qact0
    IF (debug) THEN 
      write(out_unitp,*) 'DQ',DQ(:)
      write(out_unitp,*) 'para_Tnum%dnGGref'
      CALL Write_dnSVM(para_Tnum%dnGGref)
    END IF

    dnGG%d0 = para_Tnum%dnGGref%d0

    IF (para_Tnum%GTaylor_Order > 0) THEN
      DO i=1,mole%nb_act
        dnGG%d0 = dnGG%d0 + para_Tnum%dnGGref%d1(:,:,i)*DQ(i)
      END DO
    END IF

    IF (para_Tnum%GTaylor_Order > 1) THEN
      DO i=1,mole%nb_act
      DO j=1,mole%nb_act
        dnGG%d0 = dnGG%d0 + HALF*para_Tnum%dnGGref%d2(:,:,j,i)*DQ(i)*DQ(j)
      END DO
      END DO
    END IF

    IF (dnGG%nderiv > 0) THEN
      dnGG%d1 = ZERO

      IF (para_Tnum%GTaylor_Order > 0) THEN
        dnGG%d1 = para_Tnum%dnGGref%d1
      END IF

      IF (para_Tnum%GTaylor_Order > 1) THEN
        dnGG%d1 = para_Tnum%dnGGref%d1

        DO i=1,mole%nb_act
        DO j=1,mole%nb_act
          dnGG%d1(:,:,i) = dnGG%d1(:,:,i) + para_Tnum%dnGGref%d2(:,:,j,i)*DQ(j)
        END DO
        END DO

      END IF

    END IF

    IF (dnGG%nderiv > 1) THEN
      IF (para_Tnum%GTaylor_Order > 1) THEN
        dnGG%d2 = para_Tnum%dnGGref%d2
      ELSE
        dnGG%d2 = ZERO
      END IF
    END IF

    deallocate(DQ)

    IF (debug) THEN
      write(out_unitp,*) 'dnGG'
      CALL Write_dnSVM(dnGG)
      write(out_unitp,*) 'END ',name_sub
      CALL flush_perso(out_unitp)
    END IF

  END SUBROUTINE get_dng_dnGG_WITH_GTaylor


  SUBROUTINE get_dng_dnGG_WITH_f2f1_ana(Qact,para_Tnum,mole,dng,dnGG,vep,nderiv)
  USE mod_ActiveTransfo,    only: qact_to_qdyn_from_activetransfo
  IMPLICIT NONE

    real (kind=Rkind), intent(in)               :: Qact(:)
    TYPE(Type_dnMat),  intent(inout), optional  :: dng,dnGG
    real (kind=Rkind), intent(inout), optional  :: vep
    integer,           intent(in)               :: nderiv

    !----- for the CoordType and Tnum --------------------------------
    TYPE (CoordType),  intent(in)               :: mole
    TYPE (Tnum),       intent(in)               :: para_Tnum

    !-----------------------------------------------------------------
    ! local variables
    real (kind=Rkind) :: Gref(mole%ndimG,mole%ndimG)
    real (kind=Rkind) :: Tdef2(mole%nb_act,mole%nb_act)
    real (kind=Rkind) :: Tdef1(mole%nb_act)
    real (kind=Rkind) :: Tcor2(mole%nb_act,3),Tcor1(3)
    real (kind=Rkind) :: Trot(3,3)
    real (kind=Rkind) :: rho,vep_loc
    real (kind=Rkind) :: Qdyn(mole%nb_var)
    integer           :: i
    !-----------------------------------------------------------------


    !----- for debuging --------------------------------------------------
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    character (len=*), parameter :: name_sub = 'get_dng_dnGG_WITH_f2f1_ana'
    !---------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      CALL flush_perso(out_unitp)
    END IF

    !-----------------------------------------------------------
    !-----------------------------------------------------------
    !-----------------------------------------------------------
    !-----------------------------------------------------------
    !     with anlytical expression of f2, f1
    !-----------------------------------------------------------
    !-----------------------------------------------------------
    !-----------------------------------------------------------
    CALL Qact_TO_Qdyn_FROM_ActiveTransfo(Qact,Qdyn,mole%ActiveTransfo)

    CALL calc_f2_f1Q_ana(Qdyn,                                      &
                         Tdef2,Tdef1,vep_loc,rho,                   &
                         Tcor2,Tcor1,Trot,                          &
                         para_Tnum,mole)

    IF (present(vep)) vep = vep_loc

    CALL mat_id(Gref,mole%ndimG,mole%ndimG)
    Gref(1:mole%nb_act,1:mole%nb_act) = -Tdef2
    DO i=1,mole%nb_act
        Gref(i,i) = TWO*Gref(i,i)
    END DO
    !!dnGG%alloc not present
    IF (present(dnGG)) THEN
      CALL sub_ZERO_TO_dnMat(dnGG,nderiv)
      dnGG%d0(:,:) = Gref
    END IF

    IF (present(dng)) THEN
      CALL sub_ZERO_TO_dnMat(dng,nderiv)
      CALL inv_m1_TO_m2(Gref,dng%d0,mole%ndimG,0,ZERO)
    END IF

    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'END ',name_sub
      CALL flush_perso(out_unitp)
    END IF
    !-----------------------------------------------------------

  END SUBROUTINE get_dng_dnGG_WITH_f2f1_ana
  SUBROUTINE get_dng_dnGG_WITH_type100(Qact,para_Tnum,mole,dng,dnGG,vep,vep_done,nderiv)
  IMPLICIT NONE

    real (kind=Rkind), intent(in)               :: Qact(:)
    TYPE(Type_dnMat),  intent(inout), optional  :: dng,dnGG
    real (kind=Rkind), intent(inout)            :: vep
    logical,           intent(inout)            :: vep_done
    integer,           intent(in)               :: nderiv

    !----- for the CoordType and Tnum --------------------------------
    TYPE (CoordType),  intent(in)               :: mole
    TYPE (Tnum),       intent(in)               :: para_Tnum

    !-----------------------------------------------------------------
    ! local variables
    TYPE (Type_dnMat) :: dnGG100,dnGG_temp
    TYPE (CoordType)  :: mole100
    integer           :: i,j
    real (kind=Rkind) :: rho
    !-----------------------------------------------------------------


    !----- for debuging --------------------------------------------------
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    character (len=*), parameter :: name_sub = 'get_dng_dnGG_WITH_type100'
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNIG ',name_sub
      CALL flush_perso(out_unitp)
    END IF
    vep      = ZERO
    vep_done = .FALSE.

    !---------------------------------------------------------------
    IF (mole%nb_rigid100 > 0) THEN

      !-----------------------------------------------------------------
      ! first copy mole to mole100 and compute dnGG100
      !-----------------------------------------------------------------
      IF (present(dnGG) .OR. present(dng)) THEN ! with type100 we need to calculate dnGG100 first
                                                ! then we can calculate dng form the inversion of dnGG
        mole100 = mole

        DO i=1,mole100%nb_var
          IF (mole100%ActiveTransfo%list_act_OF_Qdyn(i) == 100)         &
                           mole100%ActiveTransfo%list_act_OF_Qdyn(i) = 1
        END DO
        mole100%nb_act      = mole%nb_act  + mole%nb_rigid100
        mole100%nb_act1     = mole%nb_act1 + mole%nb_rigid100
        mole100%ndimG       = mole%ndimG   + mole%nb_rigid100
        mole100%nb_rigid100 = 0
        mole100%nb_rigid    = mole100%nb_rigid - mole%nb_rigid100

        mole100%tab_Qtransfo(:)%nb_act = mole100%nb_act

        IF (debug) THEN
          write(out_unitp,*) 'mole%nb_rigid100 > 0',mole%nb_rigid100
          write(out_unitp,*)
          CALL Write_CoordType(mole100)
        END IF

        CALL alloc_dnSVM(dnGG100,mole100%ndimG,mole100%ndimG,mole100%nb_act,nderiv)
        CALL get_dnGG(Qact,dnGG100,nderiv,para_Tnum,mole100)

        IF (para_Tnum%vep_type == 100 .AND. dnGG100%nderiv == 2) THEN
          CALL Calc_vep_rho_from_dnGG(vep,rho,Qact,dnGG100,mole,para_Tnum)
          vep_done = .TRUE.
        END IF

        IF (debug) THEN
          write(out_unitp,*) 'dnGG100'
          CALL Write_dnSVM(dnGG100)
        END IF
      END IF
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! transfert dnGG100 in dnGG and dng if they are present
      IF (present(dnGG)) THEN
        !- transfo G100 => G
        CALL dngG100_TO_dngG(dnGG100,dnGG,mole100,mole)
      END IF

      IF (present(dng)) THEN
        IF (present(dnGG)) THEN
          CALL INV_dnMat1_TO_dnMat2(dnGG,dng,nderiv)
        ELSE
          !- transfo G100 => G
          CALL alloc_dnSVM(dnGG_temp,mole%ndimG,mole%ndimG,mole%nb_act,nderiv)
          CALL dngG100_TO_dngG(dnGG100,dnGG_temp,mole100,mole)
          CALL INV_dnMat1_TO_dnMat2(dnGG_temp,dng)
          CALL dealloc_dnSVM(dnGG_temp)
        END IF
      END IF
      !-----------------------------------------------------------------

      CALL dealloc_dnSVM(dnGG100)
      CALL dealloc_CoordType(mole100)

    END IF

    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'END ',name_sub
      CALL flush_perso(out_unitp)
    END IF
    !-----------------------------------------------------------

  END SUBROUTINE get_dng_dnGG_WITH_type100
  SUBROUTINE get_dng_dnGG_WITH_type100_bug_dng(Qact,para_Tnum,mole,     &
                                               dng,dnGG,vep,vep_done,nderiv)
  IMPLICIT NONE

    real (kind=Rkind), intent(in)               :: Qact(:)
    TYPE(Type_dnMat),  intent(inout), optional  :: dng,dnGG
    real (kind=Rkind), intent(inout)            :: vep
    logical,           intent(inout)            :: vep_done
    integer,           intent(in)               :: nderiv

    !----- for the CoordType and Tnum --------------------------------
    TYPE (CoordType),  intent(in)               :: mole
    TYPE (Tnum),       intent(in)               :: para_Tnum

    !-----------------------------------------------------------------
    ! local variables
    TYPE (Type_dnMat) :: dnGG100,dng100,dng_loc
    TYPE (CoordType)  :: mole100
    integer           :: i,j
    real (kind=Rkind) :: rho
    !-----------------------------------------------------------------


    !----- for debuging --------------------------------------------------
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    character (len=*), parameter :: name_sub = 'get_dng_dnGG_WITH_type100_bug_dng'
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNIG ',name_sub
      CALL flush_perso(out_unitp)
    END IF
    vep      = ZERO
    vep_done = .FALSE.
    !---------------------------------------------------------------
    IF (mole%nb_rigid100 > 0) THEN

      !-----------------------------------------------------------------
      ! first copy mole to mole100 and compute dnGG100
      !-----------------------------------------------------------------
      IF (present(dnGG) .OR. present(dng)) THEN ! with type100 we need to calculate dnGG100 first
                                                ! then we can calculate dng form the inversion of dnGG
        mole100 = mole
        DO i=1,mole100%nb_var
          IF (mole100%ActiveTransfo%list_act_OF_Qdyn(i) == 100)         &
                           mole100%ActiveTransfo%list_act_OF_Qdyn(i) = 1
        END DO
        mole100%nb_act      = mole%nb_act  + mole%nb_rigid100
        mole100%nb_act1     = mole%nb_act1 + mole%nb_rigid100
        mole100%ndimG       = mole%ndimG   + mole%nb_rigid100
        mole100%nb_rigid100 = 0
        mole100%nb_rigid    = mole100%nb_rigid - mole%nb_rigid100

        mole100%tab_Qtransfo(:)%nb_act = mole100%nb_act

        IF (debug) THEN
          write(out_unitp,*) 'mole%nb_rigid100 > 0',mole%nb_rigid100
          write(out_unitp,*)
          CALL Write_CoordType(mole100)
          write(out_unitp,*)
        END IF

        CALL alloc_dnSVM(dnGG100,                                     &
                         mole100%ndimG,mole100%ndimG,mole100%nb_act,  &
                         nderiv)
        CALL alloc_dnSVM(dng100,                                      &
                         mole100%ndimG,mole100%ndimG,mole100%nb_act,  &
                         nderiv)

        CALL get_dng(Qact,dng100,nderiv,para_Tnum,mole100)
        CALL INV_dnMat1_TO_dnMat2(dng100,dnGG100,nderiv)

        ! Up to here, it is correct

        IF (debug) THEN
          write(out_unitp,*) 'dng100'
          CALL Write_dnSVM(dng100)
          write(out_unitp,*) 'dnGG100'
          CALL Write_dnSVM(dnGG100)
        END IF
      END IF
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! transfert dnGG100 in dnGG and compute dng if they are present
      IF (present(dnGG)) THEN
        !- transfo G100 => G: here, it is correct as well.
        CALL dngG100_TO_dngG(dnGG100,dnGG,mole100,mole)
      END IF
      IF (present(dng)) THEN ! here it is wrong. dng MUST be dnGG^-1
        CALL dngG100_TO_dngG(dng100,dng,mole100,mole)
      END IF

      IF (present(dng) .AND. present(dnGG)) THEN ! here it is wrong. because dng is wrong
        CALL Calc_vep_rho_from_dng_AND_dnGG(vep,rho,Qact,dng,dnGG,mole,para_Tnum)
        vep_done = .TRUE.
      ELSE IF (present(dnGG)) THEN ! here it is wrong. because dng is wrong
        CALL alloc_dnSVM(dng_loc,mole%ndimG,mole%ndimG,mole%nb_act,nderiv)

        CALL dngG100_TO_dngG(dng100,dng_loc,mole100,mole) ! here it is wrong. dng MUST be dnGG^-1

        CALL Calc_vep_rho_from_dng_AND_dnGG(vep,rho,Qact,dng_loc,dnGG,mole,para_Tnum)

        CALL dealloc_dnSVM(dng_loc)
        vep_done = .TRUE.

      END IF
      !-----------------------------------------------------------------

      CALL dealloc_dnSVM(dng100)
      CALL dealloc_dnSVM(dnGG100)
      CALL dealloc_CoordType(mole100)


    END IF

    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'END ',name_sub
      CALL flush_perso(out_unitp)
    END IF
    !-----------------------------------------------------------

  END SUBROUTINE get_dng_dnGG_WITH_type100_bug_dng

  SUBROUTINE get_d0GG_CoordType(Qact,para_Tnum,mole,d0GG,def,Jac,Rho)
  USE mod_dnDetGG_dnDetg, only : sub3_dnDetGG
  IMPLICIT NONE

    !-----------------------------------------------------------------
    real (kind=Rkind),           intent(in)     :: Qact(:)
    real (kind=Rkind),           intent(inout)  :: d0GG(:,:)
    logical,                     intent(in)     :: def
    real (kind=Rkind), optional, intent(inout)  :: Jac,Rho
    !----- for the CoordType and Tnum --------------------------------
    TYPE (CoordType),            intent(in)     :: mole
    TYPE (Tnum),                 intent(in)     :: para_Tnum
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! local variables
    TYPE(Type_dnMat)  :: dnGG
    TYPE(Type_dnS)    :: dnS
    !-----------------------------------------------------------------


    !----- for debuging --------------------------------------------------
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    character (len=*), parameter :: name_sub = 'get_d0GG'
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*)
      write(out_unitp,*) 'ndimG',mole%ndimG
      write(out_unitp,*)
      write(out_unitp,*) 'Qact',Qact
      write(out_unitp,*)
      !CALL Write_CoordType(mole)
      write(out_unitp,*)
      write(out_unitp,*)
      write(out_unitp,*) 'num_GG,num_g',para_Tnum%num_GG,para_Tnum%num_g
      write(out_unitp,*) 'num_x,nrho',para_Tnum%num_x,para_Tnum%nrho
      write(out_unitp,*) 'JJ',para_Tnum%JJ
      write(out_unitp,*)
      CALL flush_perso(out_unitp)
    END IF


    CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,0)

    CALL get_dng_dnGG(Qact,para_Tnum,mole,dnGG=dnGG,nderiv=0)

    IF (def) THEN
      d0GG(:,:) = dnGG%d0(1:mole%nb_act,1:mole%nb_act)
    ELSE
      d0GG(:,:) = dnGG%d0(:,:)
    END IF

    ! for Jac and Rho if required
    IF (present(Jac) .OR. present(Rho)) CALL alloc_dnSVM(dnS,mole%nb_act,nderiv=0)

    IF (present(Jac)) THEN
      ! here dnS contains dnJac
      CALL sub3_dndetGG(dnS,dnGG,0,mole%masses,mole%Mtot_inv,mole%ncart)
      Jac = dnS%d0
    END IF

    IF (present(Rho)) THEN
      ! now the true dnrho
      CALL sub3_dnrho_ana(dnS,Qact,mole,0)
      Rho = dnS%d0 ! here we should have the good rho for the basis set
    END IF

    IF (present(Jac) .OR. present(Rho)) CALL dealloc_dnSVM(dnS)
    CALL dealloc_dnSVM(dnGG)

    !-----------------------------------------------------------
    IF (debug) THEN
      CALL Write_Mat(d0GG,out_unitp,5,info='d0GG')
      write(out_unitp,*) 'END ',name_sub
    END IF
    !-----------------------------------------------------------

  END SUBROUTINE get_d0GG_CoordType

  SUBROUTINE get_d0g_d0GG(Qact,para_Tnum,mole,d0g,d0GG,def)
  IMPLICIT NONE

    !-----------------------------------------------------------------
    real (kind=Rkind),           intent(in)     :: Qact(:)
    real (kind=Rkind), optional, intent(inout)  :: d0g(:,:),d0GG(:,:)
    logical,                     intent(in)     :: def
    !----- for the CoordType and Tnum --------------------------------
    TYPE (CoordType),            intent(in)     :: mole
    TYPE (Tnum),                 intent(in)     :: para_Tnum
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! local variables
    TYPE(Type_dnMat)  :: dng,dnGG
    !-----------------------------------------------------------------


    !----- for debuging --------------------------------------------------
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    character (len=*), parameter :: name_sub = 'get_d0g_d0GG'
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*)
      write(out_unitp,*) 'ndimG',mole%ndimG
      write(out_unitp,*)
      write(out_unitp,*) 'Qact',Qact
      write(out_unitp,*)
      CALL Write_CoordType(mole)
      write(out_unitp,*)
      write(out_unitp,*)
      write(out_unitp,*) 'num_GG,num_g',para_Tnum%num_GG,para_Tnum%num_g
      write(out_unitp,*) 'num_x,nrho',para_Tnum%num_x,para_Tnum%nrho
      write(out_unitp,*) 'JJ',para_Tnum%JJ
      write(out_unitp,*)
      CALL flush_perso(out_unitp)
    END IF

    IF ( present(d0g) .AND. present(d0GG) ) THEN

      CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,0)
      CALL alloc_dnSVM(dng, mole%ndimG,mole%ndimG,mole%nb_act,0)

      CALL get_dng_dnGG(Qact,para_Tnum,mole,dng=dng,dnGG=dnGG,nderiv=0)

    ELSE IF ( .NOT. present(d0g) .AND. present(d0GG) ) THEN

      CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,0)

      CALL get_dng_dnGG(Qact,para_Tnum,mole,dnGG=dnGG,nderiv=0)

    ELSE IF ( present(d0g) .AND. .NOT. present(d0GG) ) THEN

      CALL alloc_dnSVM(dng, mole%ndimG,mole%ndimG,mole%nb_act,0)

      CALL get_dng_dnGG(Qact,para_Tnum,mole,dng=dng,nderiv=0)

    ELSE
      write(out_unitp,*) 'ERROR in ',name_sub
      write(out_unitp,*) ' d0g and d0GG are not present !!'
      write(out_unitp,*) '  check the fortran!!'
      STOP
    END IF

    IF (present(d0g)) THEN
      IF (def) THEN
        d0g(:,:) = dng%d0(1:mole%nb_act,1:mole%nb_act)
      ELSE
        d0g(:,:) = dng%d0(:,:)
      END IF
    END IF

    IF (present(d0GG)) THEN
      IF (def) THEN
        d0GG(:,:) = dnGG%d0(1:mole%nb_act,1:mole%nb_act)
      ELSE
        d0GG(:,:) = dnGG%d0(:,:)
      END IF
    END IF

    CALL dealloc_dnSVM(dnGG)
    CALL dealloc_dnSVM(dng)

    !-----------------------------------------------------------
    IF (debug) THEN
      IF (present(d0GG)) THEN
        CALL Write_Mat(d0GG,out_unitp,5,info='d0GG')
      END IF
      IF (present(d0g)) THEN
        CALL Write_Mat(d0g,out_unitp,5,info='d0g')
      END IF
      write(out_unitp,*) 'END ',name_sub
    END IF
    !-----------------------------------------------------------

  END SUBROUTINE get_d0g_d0GG

  ! this subroutine does NOT deal with type100. You have to use get_dng_dnGG
  RECURSIVE SUBROUTINE get_dnGG(Qact,dnGG,nderiv,para_Tnum,mole)
  IMPLICIT NONE

      !-----------------------------------------------------------------
      real (kind=Rkind), intent(in)               :: Qact(:)
      TYPE(Type_dnMat),  intent(inout)            :: dnGG
      integer,           intent(in)               :: nderiv
      !----- for the CoordType and Tnum --------------------------------
      TYPE (CoordType),  intent(in)               :: mole
      TYPE (Tnum),       intent(in)               :: para_Tnum
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! local variables
      real (kind=Rkind) :: Qact_loc(mole%nb_var)
      TYPE(Type_dnMat)  :: dng,dnGG_loc
      integer           :: i,j
      real (kind=Rkind) :: step2,stepp,step24
      !-----------------------------------------------------------------

      !----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'get_dnGG'
      !-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nderiv',nderiv
        write(out_unitp,*) 'ndimG',mole%ndimG
        write(out_unitp,*) 'Qact',Qact
        write(out_unitp,*) 'num_GG,num_g',para_Tnum%num_GG,para_Tnum%num_g
        write(out_unitp,*) 'num_x,nrho',para_Tnum%num_x,para_Tnum%nrho
        write(out_unitp,*) 'step',para_Tnum%stepT
        write(out_unitp,*) 'JJ',para_Tnum%JJ
        write(out_unitp,*)
        CALL flush_perso(out_unitp)
      END IF

      !---------------------------------------------------------------
      IF (para_Tnum%num_GG .AND. nderiv > 0) THEN

        IF (para_Tnum%stepT == ZERO) THEN
          write(out_unitp,*) ' ERROR : stepT is zero'
          STOP
        END IF
        step2 = ONE/(para_Tnum%stepT*para_Tnum%stepT)
        step24 = step2/FOUR
        stepp = ONE/(para_Tnum%stepT+para_Tnum%stepT)

        Qact_loc(:) = Qact(:)

        ! Calculation at Qact
        IF (nderiv >=0) THEN
          CALL get_dnGG(Qact_loc,dnGG,0,para_Tnum,mole)
        END IF

        ! Calculation at Qacti+/-step
        CALL alloc_dnSVM(dnGG_loc,mole%ndimG,mole%ndimG,mole%nb_act,0)
        IF (nderiv >=1) THEN


          DO i=1,mole%nb_act

            Qact_loc(i) = Qact(i) + para_Tnum%stepT
            CALL get_dnGG(Qact_loc,dnGG_loc,0,para_Tnum,mole)

            dnGG%d1(:,:,i) = dnGG_loc%d0
            IF (nderiv == 2) dnGG%d2(:,:,i,i) = dnGG_loc%d0

            Qact_loc(i) = Qact(i) - para_Tnum%stepT
            CALL get_dnGG(Qact_loc,dnGG_loc,0,para_Tnum,mole)

            dnGG%d1(:,:,i) = (dnGG%d1(:,:,i)-dnGG_loc%d0) * stepp
            IF (nderiv == 2) THEN
              dnGG%d2(:,:,i,i) = ( dnGG%d2(:,:,i,i) + dnGG_loc%d0 -     &
                  TWO*dnGG%d0 ) * step2
            END IF

          END DO
        END IF

        ! Calculation at Qacti+/-step and Qactj+/-step
        IF (nderiv == 2) THEN
          DO i=1,mole%nb_act
          DO j=i+1,mole%nb_act

            Qact_loc(i) = Qact(i) + para_Tnum%stepT
            Qact_loc(j) = Qact(j) + para_Tnum%stepT
            CALL get_dnGG(Qact_loc,dnGG_loc,0,para_Tnum,mole)

            dnGG%d2(:,:,i,j) = dnGG_loc%d0

            Qact_loc(i) = Qact(i) - para_Tnum%stepT
            Qact_loc(j) = Qact(j) - para_Tnum%stepT
            CALL get_dnGG(Qact_loc,dnGG_loc,0,para_Tnum,mole)

            dnGG%d2(:,:,i,j) = dnGG%d2(:,:,i,j) + dnGG_loc%d0

            Qact_loc(i) = Qact(i) - para_Tnum%stepT
            Qact_loc(j) = Qact(j) + para_Tnum%stepT
            CALL get_dnGG(Qact_loc,dnGG_loc,0,para_Tnum,mole)

            dnGG%d2(:,:,i,j) = dnGG%d2(:,:,i,j) - dnGG_loc%d0

            Qact_loc(i) = Qact(i) + para_Tnum%stepT
            Qact_loc(j) = Qact(j) - para_Tnum%stepT
            CALL get_dnGG(Qact_loc,dnGG_loc,0,para_Tnum,mole)

            dnGG%d2(:,:,i,j) = dnGG%d2(:,:,i,j) - dnGG_loc%d0


            !-- d2A/dQidQj -----------------------------------------
            dnGG%d2(:,:,i,j) = dnGG%d2(:,:,i,j) * step24
            dnGG%d2(:,:,j,i) = dnGG%d2(:,:,i,j)

          END DO
          END DO
        END IF

        CALL dealloc_dnSVM(dnGG_loc)


      ELSE
        !------------------------------------------------------------
        !------------------------------------------------------------
        CALL alloc_dnSVM(dng,mole%ndimG,mole%ndimG,mole%nb_act,nderiv)
        CALL get_dng(Qact,dng,nderiv,para_Tnum,mole)
        !write(out_unitp,*) 'dng'
        !CALL Write_dnSVM(dng)

        CALL INV_dnMat1_TO_dnMat2(dng,dnGG,nderiv)

        CALL dealloc_dnSVM(dng)

      END IF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'dnGG'
        CALL Write_dnSVM(dnGG)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------


  END SUBROUTINE get_dnGG

  ! this subroutine does NOT deal with type100. You have to use get_dng_dnGG
  RECURSIVE SUBROUTINE get_dng(Qact,dng,nderiv,para_Tnum,mole)
  IMPLICIT NONE

    !-----------------------------------------------------------------
    real (kind=Rkind), intent(in)               :: Qact(:)
    TYPE(Type_dnMat),  intent(inout)            :: dng
    integer,           intent(in)               :: nderiv
    !----- for the CoordType and Tnum --------------------------------
    TYPE (CoordType),  intent(in)               :: mole
    TYPE (Tnum),       intent(in)               :: para_Tnum
    !-----------------------------------------------------------------

    !-----------------------------------------------------------------
    ! local variables
    TYPE (Type_dnVec) :: dnMWx
    !-----------------------------------------------------------------

    !----- for debuging --------------------------------------------------
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    character (len=*), parameter :: name_sub = 'get_dng'
    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) 'nderiv',nderiv
      write(out_unitp,*) 'ndimG',mole%ndimG
      write(out_unitp,*) 'Qact',Qact
      write(out_unitp,*) 'num_x,nrho',para_Tnum%num_x,para_Tnum%nrho
      write(out_unitp,*) 'step',para_Tnum%stepT
      write(out_unitp,*) 'JJ',para_Tnum%JJ
      write(out_unitp,*)
      CALL flush_perso(out_unitp)
    END IF
    !---------------------------------------------------------------


    CALL alloc_dnSVM(dnMWx,mole%ncart,mole%nb_act,nderiv+1)

    IF (para_Tnum%num_g) THEN
      CALL sub3_dnA_num(Qact,dng,dnMWx,mole,para_Tnum,nderiv)
    ELSE
      CALL sub3_dnA_ana(Qact,dng,dnMWx,mole,para_Tnum,nderiv)
    END IF

    IF (debug) CALL analyze_dnx(dnMWx,Qact,mole)

    CALL dealloc_dnSVM(dnMWx)

    !-----------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'dng'
      CALL Write_dnSVM(dng)
      write(out_unitp,*) 'END ',name_sub
    END IF
    !-----------------------------------------------------------

  END SUBROUTINE get_dng
!
!=====================================================================
!
! ++   analytical derivative of d0A => d1A, d2A
!
!=====================================================================
!
  SUBROUTINE sub3_dnA_ana(Qact,dnA,dnMWx,mole,para_Tnum,nderivA)
  IMPLICIT NONE

      real (kind=Rkind), intent(in)               :: Qact(:)
      TYPE(Type_dnMat),  intent(inout)            :: dnA
      TYPE(Type_dnVec),  intent(inout)            :: dnMWx
      integer,           intent(in)               :: nderivA

      !----- for the CoordType and Tnum --------------------------------
      TYPE (CoordType),  intent(in)               :: mole
      TYPE (Tnum),       intent(in)               :: para_Tnum

      !-----------------------------------------------------------------
      ! local variables
      integer :: nderivX
      logical :: Gcenter

!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'sub3_dnA_ana'
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'ndimA',dnA%nb_var_Matl
         write(out_unitp,*) 'ndimA',dnA%nb_var_Matc
         write(out_unitp,*) 'nderivA',nderivA
         write(out_unitp,*) 'Qact',Qact
         write(out_unitp,*)
         CALL Write_CoordType(mole)
         write(out_unitp,*)
       END IF
!-----------------------------------------------------------

!===================================================================
!
!       Calcul en Qact
!
!===================================================================
      Gcenter = .TRUE.
      nderivX = nderivA + 1

      CALL sub_QactTOdnMWx(Qact,dnMWx,mole,nderivX,Gcenter,             &
                               Cart_Transfo=para_Tnum%With_Cart_Transfo)

      IF (nderivA >= 0) THEN
        CALL sub_d0A(dnA%d0,dnMWx%d0,dnMWx%d1,                          &
                     mole%nb_act,mole%nat_act,mole%ncart,               &
                     mole%ncart_act,dnA%nb_var_Matl,                    &
                     mole%Without_Rot,mole%With_VecCOM,mole%Mtot)

!test linear dep:
!       CALL Linear_Dep(dnA%d0(1:mole%nb_act,1:mole%nb_act))

      END IF

      IF (nderivA >= 1) THEN
        CALL sub_d1A(dnA%d1,dnMWx%d0,dnMWx%d1,dnMWx%d2,                 &
                     mole%nb_act,mole%nat_act,mole%ncart,               &
                     mole%ncart_act,dnA%nb_var_Matl,                    &
                     mole%Without_Rot)
      END IF

      IF (nderivA == 2) THEN
        CALL sub_d2A(dnA%d2,dnMWx%d0,dnMWx%d1,dnMWx%d2,dnMWx%d3,         &
                     mole%nb_act,mole%nat_act,mole%ncart,               &
                     mole%ncart_act,dnA%nb_var_Matl,                    &
                     mole%Without_Rot)
      END IF

!-----------------------------------------------------------
      IF (debug) THEN

        write(out_unitp,*) 'dnMWx'
        CALL write_dnx(1,mole%ncart,dnMWx,nderivX)

        write(out_unitp,*) 'dnA'
        CALL Write_dnSVM(dnA)

        write(out_unitp,*) 'END ',name_sub
       END IF
!-----------------------------------------------------------
  END SUBROUTINE sub3_dnA_ana
!
!=====================================================================
!
! ++   numerical derivative of d0A => d1A, d2A
!
!=====================================================================
!
  SUBROUTINE sub3_dnA_num(Qact,dnA,dnMWx,mole,para_Tnum,nderivA)
  IMPLICIT NONE

      real (kind=Rkind), intent(in)               :: Qact(:)
      TYPE(Type_dnMat),  intent(inout)            :: dnA
      TYPE(Type_dnVec),  intent(inout)            :: dnMWx
      integer,           intent(in)               :: nderivA

      !----- for the CoordType and Tnum --------------------------------
      TYPE (CoordType),  intent(in)               :: mole
      TYPE (Tnum),       intent(in)               :: para_Tnum

      !-----------------------------------------------------------------
      ! local variables
      real (kind=Rkind) :: Qact_loc(mole%nb_var)
      logical           :: Gcenter
      real (kind=Rkind) :: step2,step24,stepp
      integer           :: i,j
      real (kind=Rkind), allocatable :: At(:,:)

!----- for debuging --------------------------------------------------

      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'sub3_dnA_num'
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'ndimA',dnA%nb_var_Matl
         write(out_unitp,*) 'ndimA',dnA%nb_var_Matc
         write(out_unitp,*) 'nderivA',nderivA
         write(out_unitp,*) 'num_x,nrho',para_Tnum%num_x,para_Tnum%nrho
         write(out_unitp,*) 'step',para_Tnum%stepT
         write(out_unitp,*) 'Qact',Qact
         write(out_unitp,*)
         CALL Write_CoordType(mole)
         write(out_unitp,*)
       END IF
!-----------------------------------------------------------


!----- some step ----------------------------------------------------
      Gcenter = .TRUE.
      IF (para_Tnum%stepT == ZERO) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' para_Tnum%stepT is zero'
        STOP ' ERROR : para_Tnum%stepT is zero'
      END IF
      step2 = ONE/(para_Tnum%stepT*para_Tnum%stepT)
      step24 = step2/FOUR
      stepp = ONE/(para_Tnum%stepT+para_Tnum%stepT)
      Qact_loc(:) = Qact(:)

      CALL alloc_NParray(At,[dnA%nb_var_Matl,dnA%nb_var_Matc],"At",name_sub)

!===================================================================
!
!       Calcul en Qact
!
!===================================================================
      IF (nderivA >=0) THEN

        CALL sub_QactTOdnMWx(Qact_loc,dnMWx,mole,1,Gcenter,             &
                               Cart_Transfo=para_Tnum%With_Cart_Transfo)

        CALL sub_d0A(dnA%d0,dnMWx%d0,dnMWx%d1,                          &
                     mole%nb_act,mole%nat_act,mole%ncart,               &
                     mole%ncart_act,dnA%nb_var_Matl,                    &
                     mole%Without_Rot,mole%With_VecCOM,mole%Mtot)
      END IF

!===================================================================
!
!       Calcul en Qact(i)+step
!           et en Qact(i)-step
!
!===================================================================

      IF (nderivA >=1) THEN
        DO i=1,mole%nb_act
          Qact_loc(:) = Qact(:)

          Qact_loc(i) = Qact(i) + para_Tnum%stepT
          CALL sub_QactTOdnMWx(Qact_loc,dnMWx,mole,1,Gcenter,            &
                               Cart_Transfo=para_Tnum%With_Cart_Transfo)

          CALL sub_d0A(At,dnMWx%d0,dnMWx%d1,                            &
                       mole%nb_act,mole%nat_act,mole%ncart,             &
                       mole%ncart_act,dnA%nb_var_Matl,                  &
                       mole%Without_Rot,mole%With_VecCOM,mole%Mtot)
          dnA%d1(:,:,i) = At(:,:)
          IF (nderivA == 2) dnA%d2(:,:,i,i) = At(:,:)

          Qact_loc(i) = Qact(i) - para_Tnum%stepT
          CALL sub_QactTOdnMWx(Qact_loc,dnMWx,mole,1,Gcenter,           &
                               Cart_Transfo=para_Tnum%With_Cart_Transfo)
          CALL sub_d0A(At,dnMWx%d0,dnMWx%d1,                            &
                       mole%nb_act,mole%nat_act,mole%ncart,             &
                       mole%ncart_act,dnA%nb_var_Matl,                  &
                       mole%Without_Rot,mole%With_VecCOM,mole%Mtot)
          dnA%d1(:,:,i) = ( dnA%d1(:,:,i) - At(:,:) ) * stepp

          IF (nderivA == 2) dnA%d2(:,:,i,i) = ( dnA%d2(:,:,i,i) +       &
                                  At(:,:) - TWO*dnA%d0(:,:) ) * step2

        END DO
      END IF

!===================================================================
!
!       Calcul en Qact(i)+/-step
!           et en Qact(j)+/-step
!
!===================================================================

      IF (nderivA ==2) THEN
        DO i=1,mole%nb_act
        DO j=i+1,mole%nb_act
          Qact_loc(:) = Qact(:)

          Qact_loc(i) = Qact(i) + para_Tnum%stepT
          Qact_loc(j) = Qact(j) + para_Tnum%stepT
          CALL sub_QactTOdnMWx(Qact_loc,dnMWx,mole,1,Gcenter,           &
                               Cart_Transfo=para_Tnum%With_Cart_Transfo)
          CALL sub_d0A(At,dnMWx%d0,dnMWx%d1,                            &
                       mole%nb_act,mole%nat_act,mole%ncart,             &
                       mole%ncart_act,dnA%nb_var_Matl,                  &
                       mole%Without_Rot,mole%With_VecCOM,mole%Mtot)
          dnA%d2(:,:,i,j) = At(:,:)

          Qact_loc(i) = Qact(i) - para_Tnum%stepT
          Qact_loc(j) = Qact(j) - para_Tnum%stepT
          CALL sub_QactTOdnMWx(Qact_loc,dnMWx,mole,1,Gcenter,           &
                               Cart_Transfo=para_Tnum%With_Cart_Transfo)
          CALL sub_d0A(At,dnMWx%d0,dnMWx%d1,                            &
                       mole%nb_act,mole%nat_act,mole%ncart,             &
                       mole%ncart_act,dnA%nb_var_Matl,                  &
                       mole%Without_Rot,mole%With_VecCOM,mole%Mtot)
          dnA%d2(:,:,i,j) = dnA%d2(:,:,i,j) + At(:,:)

          Qact_loc(i) = Qact(i) - para_Tnum%stepT
          Qact_loc(j) = Qact(j) + para_Tnum%stepT
          CALL sub_QactTOdnMWx(Qact_loc,dnMWx,mole,1,Gcenter,           &
                               Cart_Transfo=para_Tnum%With_Cart_Transfo)
          CALL sub_d0A(At,dnMWx%d0,dnMWx%d1,                            &
                       mole%nb_act,mole%nat_act,mole%ncart,             &
                       mole%ncart_act,dnA%nb_var_Matl,                  &
                       mole%Without_Rot,mole%With_VecCOM,mole%Mtot)
          dnA%d2(:,:,i,j) = dnA%d2(:,:,i,j) - At(:,:)

          Qact_loc(i) = Qact(i) + para_Tnum%stepT
          Qact_loc(j) = Qact(j) - para_Tnum%stepT
          CALL sub_QactTOdnMWx(Qact_loc,dnMWx,mole,1,Gcenter,           &
                               Cart_Transfo=para_Tnum%With_Cart_Transfo)
          CALL sub_d0A(At,dnMWx%d0,dnMWx%d1,                            &
                       mole%nb_act,mole%nat_act,mole%ncart,             &
                       mole%ncart_act,dnA%nb_var_Matl,                  &
                       mole%Without_Rot,mole%With_VecCOM,mole%Mtot)
          dnA%d2(:,:,i,j) = dnA%d2(:,:,i,j) - At(:,:)

!        -- d2A/dQidQj -----------------------------------------

          dnA%d2(:,:,i,j) = dnA%d2(:,:,i,j) * step24
          dnA%d2(:,:,j,i) = dnA%d2(:,:,i,j)

        END DO
        END DO
      END IF


      CALL dealloc_NParray(At,"At",name_sub)

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'nderivA',nderivA
        write(out_unitp,*) 'dnA'
        CALL Write_dnSVM(dnA)

        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

  END SUBROUTINE sub3_dnA_num
!
!================================================================
!
!    Calcul de la matrice A
!
!================================================================

      SUBROUTINE sub_d0A(A,d0x,d1x,                                     &
                         nb_act,nat_act,ncart,ncart_act,ndimA,          &
                         Without_Rot,With_VecCOM,Mtot)
      USE mod_paramQ, ONLY : Write_d0Q
      IMPLICIT NONE


       integer :: nb_act,nat_act,ncart,ncart_act,ndimA
       real (kind=Rkind) :: A(ndimA,ndimA),Mtot
       real (kind=Rkind) :: d0x(ncart)
       real (kind=Rkind) :: d1x(ncart,nb_act)


       real (kind=Rkind) :: Ixx,Iyy,Izz,Ixy,Ixz,Iyz

       real (kind=Rkind) :: Cxq,Cyq,Czq

       integer :: i,j,k,kx,ky,kz
       logical :: Without_Rot,With_VecCOM

!----- for debuging --------------------------------------------------
       logical, parameter :: debug = .FALSE.
       !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING sub_d0A'
         write(out_unitp,*) 'nb_act,nat_act,ncart,ncart_act,ndimA',             &
                     nb_act,nat_act,ncart,ncart_act,ndimA
         write(out_unitp,*) 'Without_Rot',Without_Rot
         CALL Write_d0Q(0,'d0x',d0x,3)
         write(out_unitp,*)
       END IF
!-----------------------------------------------------------

!================================================================
!    Calcul de la matrice S : 1er bloc de A
!================================================================
       DO i=1,nb_act
         DO j=i,nb_act
           A(i,j) = ZERO
           DO k=1,ncart_act
             A(i,j) = A(i,j) + d1x(k,i)*d1x(k,j)
           END DO
           A(j,i) = A(i,j)
         END DO
       END DO
!================================================================


       IF (Without_Rot) RETURN

!================================================================
!    Calcul de la matrice I : 2d bloc sur la diagonale de A
!================================================================
       Ixx = ZERO
       Iyy = ZERO
       Izz = ZERO
       Ixy = ZERO
       Ixz = ZERO
       Iyz = ZERO

       kx = 1
       ky = 2
       kz = 3
       DO k=1,nat_act
         Ixx = Ixx +  d0x(kx)*d0x(kx)
         Iyy = Iyy +  d0x(ky)*d0x(ky)
         Izz = Izz +  d0x(kz)*d0x(kz)
         Ixy = Ixy +  d0x(kx)*d0x(ky)
         Ixz = Ixz +  d0x(kx)*d0x(kz)
         Iyz = Iyz +  d0x(ky)*d0x(kz)

         kx = kx + 3
         ky = ky + 3
         kz = kz + 3
       END DO


       A(nb_act+1,nb_act+1) =  Iyy+Izz
       A(nb_act+2,nb_act+2) =  Ixx+Izz
       A(nb_act+3,nb_act+3) =  Ixx+Iyy

       A(nb_act+1,nb_act+2) = -Ixy
       A(nb_act+1,nb_act+3) = -Ixz
       A(nb_act+2,nb_act+3) = -Iyz

       A(nb_act+2,nb_act+1) = -Ixy
       A(nb_act+3,nb_act+1) = -Ixz
       A(nb_act+3,nb_act+2) = -Iyz
!================================================================


!================================================================
!    Calcul de la matrice C : bloc hors de la diagonale de A
!================================================================
       DO i=1,nb_act

         kx = 1
         ky = 2
         kz = 3

         Cxq = ZERO
         Cyq = ZERO
         Czq = ZERO

         DO k=1,nat_act

           Cxq = Cxq +                                                  &
                 (d0x(ky)*d1x(kz,i)-d0x(kz)*d1x(ky,i))
           Cyq = Cyq +                                                  &
                 (d0x(kz)*d1x(kx,i)-d0x(kx)*d1x(kz,i))
           Czq = Czq +                                                  &
                 (d0x(kx)*d1x(ky,i)-d0x(ky)*d1x(kx,i))

           kx = kx + 3
           ky = ky + 3
           kz = kz + 3
         END DO
         A(i,nb_act+1) = Cxq
         A(i,nb_act+2) = Cyq
         A(i,nb_act+3) = Czq
         A(nb_act+1,i) = Cxq
         A(nb_act+2,i) = Cyq
         A(nb_act+3,i) = Czq
       END DO
!================================================================



!================================================================
!  spaecial case when nat=2
!================================================================
       IF ( nat_act == 2 ) THEN
          A(ndimA,ndimA) = ONETENTH**9
       END IF


!================================================================

      IF (With_VecCOM) THEN
         DO i=ndimA-3,ndimA
           A(i,i) = Mtot
         END DO
      END IF


!-----------------------------------------------------------
       IF (debug) THEN
         CALL Write_Mat(A,out_unitp,4,Rformat='e30.20',info='d0A')
         write(out_unitp,*) 'END sub_d0A'
       END IF
!-----------------------------------------------------------

       END SUBROUTINE sub_d0A
!
!================================================================
!
!    Calcul de la matrice d1A
!
!================================================================

      SUBROUTINE sub_d1A(d1A,d0x,d1x,d2x,                               &
                         nb_act,nat_act,ncart,ncart_act,ndimA,          &
                         Without_Rot)
      IMPLICIT NONE

       integer :: nb_act,nat_act,ncart,ncart_act,ndimA
       real (kind=Rkind) :: d1A(ndimA,ndimA,nb_act)
       real (kind=Rkind) :: d0x(ncart)
       real (kind=Rkind) :: d1x(ncart,nb_act)
       real (kind=Rkind) :: d2x(ncart,nb_act,nb_act)


       real (kind=Rkind) :: Ixx,Iyy,Izz,Ixy,Ixz,Iyz

       real (kind=Rkind) :: Cxq,Cyq,Czq

       integer :: i,j,k,ii,kx,ky,kz
       logical :: Without_Rot

!----- for debuging --------------------------------------------------
       logical, parameter :: debug = .FALSE.
       !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING sub_d1A'
         write(out_unitp,*) 'nb_act,nat_act,ncart,ncart_act,ndimA',             &
                     nb_act,nat_act,ncart,ncart_act,ndimA
       END IF
!-----------------------------------------------------------

!================================================================
!      derive par rapport a Qii
       DO ii=1,nb_act
!      ==========================================================
!      Calcul de la matrice S : 1er bloc de A
!      ==========================================================
       DO i=1,nb_act
         DO j=i,nb_act
           d1A(i,j,ii) = ZERO
           DO k=1,ncart_act
             d1A(i,j,ii) = d1A(i,j,ii) +                                &
                                      d2x(k,i,ii) * d1x(k,j) +          &
                                      d1x(k,i)    * d2x(k,j,ii)
           END DO
           d1A(j,i,ii) = d1A(i,j,ii)
         END DO
       END DO
!      ==========================================================



       IF (Without_Rot) CYCLE

!      ==========================================================
!      Calcul de la matrice I : 2d bloc sur la diagonale de A
!      ==========================================================
       Ixx = ZERO
       Iyy = ZERO
       Izz = ZERO
       Ixy = ZERO
       Ixz = ZERO
       Iyz = ZERO

       kx = 1
       ky = 2
       kz = 3
       DO k=1,nat_act
         Ixx = Ixx +  TWO * d1x(kx,ii) * d0x(kx)
         Iyy = Iyy +  TWO * d1x(ky,ii) * d0x(ky)
         Izz = Izz +  TWO * d1x(kz,ii) * d0x(kz)
         Ixy = Ixy +  d1x(kx,ii)*d0x(ky) + d0x(kx)*d1x(ky,ii)
         Ixz = Ixz +  d1x(kx,ii)*d0x(kz) + d0x(kx)*d1x(kz,ii)
         Iyz = Iyz +  d1x(ky,ii)*d0x(kz) + d0x(ky)*d1x(kz,ii)

         kx = kx + 3
         ky = ky + 3
         kz = kz + 3
       END DO


       d1A(nb_act+1,nb_act+1,ii) =  Iyy+Izz
       d1A(nb_act+2,nb_act+2,ii) =  Ixx+Izz
       d1A(nb_act+3,nb_act+3,ii) =  Ixx+Iyy

       d1A(nb_act+1,nb_act+2,ii) = -Ixy
       d1A(nb_act+1,nb_act+3,ii) = -Ixz
       d1A(nb_act+2,nb_act+3,ii) = -Iyz

       d1A(nb_act+2,nb_act+1,ii) = -Ixy
       d1A(nb_act+3,nb_act+1,ii) = -Ixz
       d1A(nb_act+3,nb_act+2,ii) = -Iyz
!      ==========================================================


!      ==========================================================
!      Calcul de la matrice C : bloc hors de la diagonale de A
!      ==========================================================
       DO i=1,nb_act

         kx = 1
         ky = 2
         kz = 3

         Cxq = ZERO
         Cyq = ZERO
         Czq = ZERO

         DO k=1,nat_act

           Cxq = Cxq +                                                  &
                 (d1x(ky,ii)*d1x(kz,i)-d1x(kz,ii)*d1x(ky,i)) +          &
                 (d0x(ky)*d2x(kz,i,ii)-d0x(kz)*d2x(ky,i,ii))
           Cyq = Cyq +                                                  &
                 (d1x(kz,ii)*d1x(kx,i)-d1x(kx,ii)*d1x(kz,i)) +          &
                 (d0x(kz)*d2x(kx,i,ii)-d0x(kx)*d2x(kz,i,ii))
           Czq = Czq +                                                  &
                 (d1x(kx,ii)*d1x(ky,i)-d1x(ky,ii)*d1x(kx,i)) +          &
                 (d0x(kx)*d2x(ky,i,ii)-d0x(ky)*d2x(kx,i,ii))

           kx = kx + 3
           ky = ky + 3
           kz = kz + 3
         END DO
         d1A(i,nb_act+1,ii) = Cxq
         d1A(i,nb_act+2,ii) = Cyq
         d1A(i,nb_act+3,ii) = Czq
         d1A(nb_act+1,i,ii) = Cxq
         d1A(nb_act+2,i,ii) = Cyq
         d1A(nb_act+3,i,ii) = Czq
       END DO
       END DO

      IF (debug) THEN
        DO ii=1,nb_act
          CALL Write_Mat(d1A(:,:,ii),out_unitp,4,Rformat='e30.20',info='d1A')
        END DO
        write(out_unitp,*) 'END sub_d1A'
      END IF

END SUBROUTINE sub_d1A
!
!================================================================
!
!    Calcul de la matrice d2A
!
!================================================================

      SUBROUTINE sub_d2A(d2A,d0x,d1x,d2x,d3x,                           &
                         nb_act,nat_act,ncart,ncart_act,ndimA,          &
                         Without_Rot)
      IMPLICIT NONE

       integer :: nb_act,nat_act,ncart,ncart_act,ndimA
       real (kind=Rkind) :: d2A(ndimA,ndimA,nb_act,nb_act)
       real (kind=Rkind) :: d0x(ncart)
       real (kind=Rkind) :: d1x(ncart,nb_act)
       real (kind=Rkind) :: d2x(ncart,nb_act,nb_act)
       real (kind=Rkind) :: d3x(ncart,nb_act,nb_act,nb_act)

       real (kind=Rkind) :: Ixx,Iyy,Izz,Ixy,Ixz,Iyz

       real (kind=Rkind) :: Cxq,Cyq,Czq

       integer :: i,j,k,ii,jj,kx,ky,kz
       logical :: Without_Rot

!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING sub_d2A'
         write(out_unitp,*) 'nb_act,nat_act,ncart,ncart_act,ndimA',             &
                     nb_act,nat_act,ncart,ncart_act,ndimA
       END IF
!-----------------------------------------------------------

!================================================================
!      derive par rapport a Qii Qjj
       DO ii=1,nb_act
       DO jj=ii,nb_act
!      ==========================================================
!      Calcul de la matrice S : 1er bloc de A
!      ==========================================================
       DO i=1,nb_act
         DO j=i,nb_act
           d2A(i,j,ii,jj) = ZERO
           DO k=1,ncart_act
             d2A(i,j,ii,jj) = d2A(i,j,ii,jj) +                          &
                                      d3x(k,i,ii,jj) * d1x(k,j) +       &
                                      d2x(k,i,ii) * d2x(k,j,jj) +       &
                                      d2x(k,i,jj) * d2x(k,j,ii) +       &
                                      d1x(k,i)    * d3x(k,j,ii,jj)
           END DO
           d2A(j,i,ii,jj) = d2A(i,j,ii,jj)
           d2A(j,i,jj,ii) = d2A(i,j,ii,jj)
           d2A(i,j,jj,ii) = d2A(i,j,ii,jj)
         END DO
       END DO
!      ==========================================================


       IF (Without_Rot) CYCLE

!      ==========================================================
!      Calcul de la matrice I : 2d bloc sur la diagonale de A
!      ==========================================================
       Ixx = ZERO
       Iyy = ZERO
       Izz = ZERO
       Ixy = ZERO
       Ixz = ZERO
       Iyz = ZERO

       kx = 1
       ky = 2
       kz = 3
       DO k=1,nat_act
         Ixx = Ixx +  TWO * (                                           &
                              d2x(kx,ii,jj) * d0x(kx)    +              &
                              d1x(kx,ii)    * d1x(kx,jj) )
         Iyy = Iyy +  TWO * (                                           &
                              d2x(ky,ii,jj) * d0x(ky)    +              &
                              d1x(ky,ii)    * d1x(ky,jj) )
         Izz = Izz +  TWO * (                                           &
                              d2x(kz,ii,jj) * d0x(kz)    +              &
                              d1x(kz,ii)    * d1x(kz,jj) )

         Ixy = Ixy +  d2x(kx,ii,jj)*d0x(ky) + d0x(kx)*d2x(ky,ii,jj) +   &
                      d1x(kx,ii)*d1x(ky,jj) + d1x(kx,jj)*d1x(ky,ii)

         Ixz = Ixz +  d2x(kx,ii,jj)*d0x(kz) + d0x(kx)*d2x(kz,ii,jj) +   &
                      d1x(kx,ii)*d1x(kz,jj) + d1x(kx,jj)*d1x(kz,ii)

         Iyz = Iyz +  d2x(ky,ii,jj)*d0x(kz) + d0x(ky)*d2x(kz,ii,jj) +   &
                      d1x(ky,ii)*d1x(kz,jj) + d1x(ky,jj)*d1x(kz,ii)

         kx = kx + 3
         ky = ky + 3
         kz = kz + 3
       END DO


       d2A(nb_act+1,nb_act+1,ii,jj) =  Iyy+Izz
       d2A(nb_act+2,nb_act+2,ii,jj) =  Ixx+Izz
       d2A(nb_act+3,nb_act+3,ii,jj) =  Ixx+Iyy

       d2A(nb_act+1,nb_act+2,ii,jj) = -Ixy
       d2A(nb_act+1,nb_act+3,ii,jj) = -Ixz
       d2A(nb_act+2,nb_act+3,ii,jj) = -Iyz

       d2A(nb_act+2,nb_act+1,ii,jj) = -Ixy
       d2A(nb_act+3,nb_act+1,ii,jj) = -Ixz
       d2A(nb_act+3,nb_act+2,ii,jj) = -Iyz

       d2A(nb_act+1,nb_act+1,jj,ii) = d2A(nb_act+1,nb_act+1,ii,jj)
       d2A(nb_act+2,nb_act+2,jj,ii) = d2A(nb_act+2,nb_act+2,ii,jj)
       d2A(nb_act+3,nb_act+3,jj,ii) = d2A(nb_act+3,nb_act+3,ii,jj)

       d2A(nb_act+1,nb_act+2,jj,ii) = d2A(nb_act+1,nb_act+2,ii,jj)
       d2A(nb_act+1,nb_act+3,jj,ii) = d2A(nb_act+1,nb_act+3,ii,jj)
       d2A(nb_act+2,nb_act+3,jj,ii) = d2A(nb_act+2,nb_act+3,ii,jj)

       d2A(nb_act+2,nb_act+1,jj,ii) = d2A(nb_act+2,nb_act+1,ii,jj)
       d2A(nb_act+3,nb_act+1,jj,ii) = d2A(nb_act+3,nb_act+1,ii,jj)
       d2A(nb_act+3,nb_act+2,jj,ii) = d2A(nb_act+3,nb_act+2,ii,jj)
!      ==========================================================


!      ==========================================================
!      Calcul de la matrice C : bloc hors de la diagonale de A
!      ==========================================================
       DO i=1,nb_act

         kx = 1
         ky = 2
         kz = 3

         Cxq = ZERO
         Cyq = ZERO
         Czq = ZERO

         DO k=1,nat_act

           Cxq = Cxq +                                                  &
               ( d2x(ky,ii,jj)*d1x(kz,i) + d1x(ky,ii)*d2x(kz,i,jj) ) -  &
               ( d2x(kz,ii,jj)*d1x(ky,i) + d1x(kz,ii)*d2x(ky,i,jj) ) +  &
               ( d1x(ky,jj)*d2x(kz,i,ii) + d0x(ky)*d3x(kz,i,ii,jj) ) -  &
               ( d1x(kz,jj)*d2x(ky,i,ii) + d0x(kz)*d3x(ky,i,ii,jj) )
           Cyq = Cyq +                                                  &
               ( d2x(kz,ii,jj)*d1x(kx,i) + d1x(kz,ii)*d2x(kx,i,jj) ) -  &
               ( d2x(kx,ii,jj)*d1x(kz,i) + d1x(kx,ii)*d2x(kz,i,jj) ) +  &
               ( d1x(kz,jj)*d2x(kx,i,ii) + d0x(kz)*d3x(kx,i,ii,jj) ) -  &
               ( d1x(kx,jj)*d2x(kz,i,ii) + d0x(kx)*d3x(kz,i,ii,jj) )
           Czq = Czq +                                                  &
               ( d2x(kx,ii,jj)*d1x(ky,i) + d1x(kx,ii)*d2x(ky,i,jj) ) -  &
               ( d2x(ky,ii,jj)*d1x(kx,i) + d1x(ky,ii)*d2x(kx,i,jj) ) +  &
               ( d1x(kx,jj)*d2x(ky,i,ii) + d0x(kx)*d3x(ky,i,ii,jj) ) -  &
               ( d1x(ky,jj)*d2x(kx,i,ii) + d0x(ky)*d3x(kx,i,ii,jj) )

           kx = kx + 3
           ky = ky + 3
           kz = kz + 3
         END DO
         d2A(i,nb_act+1,ii,jj) = Cxq
         d2A(i,nb_act+2,ii,jj) = Cyq
         d2A(i,nb_act+3,ii,jj) = Czq
         d2A(nb_act+1,i,ii,jj) = Cxq
         d2A(nb_act+2,i,ii,jj) = Cyq
         d2A(nb_act+3,i,ii,jj) = Czq

         d2A(i,nb_act+1,jj,ii) = Cxq
         d2A(i,nb_act+2,jj,ii) = Cyq
         d2A(i,nb_act+3,jj,ii) = Czq
         d2A(nb_act+1,i,jj,ii) = Cxq
         d2A(nb_act+2,i,jj,ii) = Cyq
         d2A(nb_act+3,i,jj,ii) = Czq
       END DO
!      ==========================================================
       END DO
       END DO
!================================================================

  IF (debug) THEN
    DO ii=1,nb_act
    DO jj=1,nb_act
      CALL Write_Mat(d2A(:,:,ii,jj),out_unitp,4,Rformat='e30.20',info='d2A')
    END DO
    END DO
    write(out_unitp,*) 'END sub_d2A'
  END IF
END SUBROUTINE sub_d2A
FUNCTION get_vep_rho(Qact,mole,para_Tnum,dnGG,rho,WithTaylor) RESULT (vep)
  USE mod_dnDetGG_dnDetg, only : sub3_dnDetGG
  IMPLICIT NONE

    real (kind=Rkind)                             :: vep

    !-----------------------------------------------------------------
    real (kind=Rkind), intent(in)               :: Qact(:)
    TYPE(Type_dnMat),  intent(in),    optional  :: dnGG
    real (kind=Rkind), intent(inout), optional  :: rho
    logical,           intent(in),    optional  :: WithTaylor

    !----- for the CoordType and Tnum --------------------------------
    TYPE (CoordType),  intent(in)               :: mole
    TYPE (Tnum),       intent(in)               :: para_Tnum
    !-----------------------------------------------------------------


    !-----------------------------------------------------------------
    ! local variables
    TYPE(Type_dnS)    :: dnJac,dnrho
    TYPE(Type_dnMat)  :: dnGG_loc

    integer           :: nderiv_loc
    logical           :: dnGG_OK,WithTaylor_loc,vep_done,rho_done
    !-----------------------------------------------------------------


    !--- for debuging --------------------------------------------------
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    character (len=*), parameter :: name_sub = 'get_vep_rho'
    !---------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) 'nrho       ',para_Tnum%nrho
      write(out_unitp,*) 'vep_type   ',para_Tnum%vep_type
      write(out_unitp,*)
      write(out_unitp,*) 'Qact',Qact

      IF (present(dnGG)) THEN
        write(out_unitp,*) 'dnGG'
        write(out_unitp,*) 'ndimG      ',mole%ndimG
        write(out_unitp,*) 'dnGG%nderiv',dnGG%nderiv
        CALL write_dnMat(dnGG)
      END IF
      !write(out_unitp,*)
      !CALL Write_CoordType(mole)
      !write(out_unitp,*)
      flush(out_unitp)
    END IF
    !---------------------------------------------------------

    IF (present(WithTaylor)) THEN 
      WithTaylor_loc = (para_Tnum%vepTaylor_Order > -1 .AND. WithTaylor)
    ELSE
      WithTaylor_loc = (para_Tnum%vepTaylor_Order > -1)
    END IF
    IF (debug) THEN
      write(out_unitp,*) 'WithTaylor_loc',WithTaylor_loc
      flush(out_unitp)
    END IF

    vep_done = .FALSE.
    rho_done = .FALSE.

    !------------------------------------------------------------
    !-- extra potential term ------------------------------------
    IF ( para_Tnum%vep_type == 1 .OR. para_Tnum%vep_type == 100) THEN

      IF (WithTaylor_loc) THEN
        vep      = get_vepTaylor(Qact,mole,para_Tnum)
        vep_done = .TRUE.
      END IF

      IF (present(dnGG)) THEN 
        dnGG_OK = (dnGG%nderiv == 2)
      ELSE
        dnGG_OK = .FALSE.
      END IF

      IF (debug) THEN
        write(out_unitp,*) 'dnGG_OK,vep_done',dnGG_OK,vep_done
        flush(out_unitp)
      END IF

      IF (.NOT. vep_done) THEN

        !----------------------------------------------------
        !-- jac0,jaci,JACij calculation
        !----------------------------------------------------
        nderiv_loc = 2
        CALL alloc_dnSVM(dnJac,mole%nb_act,nderiv=nderiv_loc)
        CALL alloc_dnSVM(dnrho,mole%nb_act,nderiv=nderiv_loc)

        IF (dnGG_OK) THEN
          CALL sub3_dnDetGG(dnJac,dnGG,nderiv_loc,mole%masses,mole%Mtot_inv,mole%ncart)
        ELSE
          CALL alloc_dnSVM(dnGG_loc,mole%ndimG,mole%ndimG,mole%nb_act,nderiv=nderiv_loc)

          CALL get_dng_dnGG(Qact,para_Tnum,mole,dnGG=dnGG_loc,nderiv=nderiv_loc)
          CALL sub3_dnDetGG(dnJac,dnGG_loc,nderiv_loc,mole%masses,mole%Mtot_inv,mole%ncart)
        END IF

        IF (debug) THEN
          write(out_unitp,*) 'dnJac'
          CALL write_dnS(dnJac)
          flush(out_unitp)
        END IF

        !----------------------------------------------------
        !-- f0,fi,Fij calculation
        !----------------------------------------------------
        CALL sub3_dnrho(dnrho,dnJac,Qact,mole,                          &
                        nderiv_loc,para_Tnum%num_x,para_Tnum%stepT,     &
                        para_Tnum%nrho)

        IF (debug) THEN
          write(out_unitp,*) 'dnrho'
          CALL write_dnS(dnrho)
          flush(out_unitp)
        END IF

        !----------------------------------------------------
        !-- vep calculation
        !----------------------------------------------------
        IF (dnGG_OK) THEN
          IF (new_vep) THEN !new vep
            CALL sub_vep_new(vep,dnGG%d0,dnGG%d1,dnrho%d1,dnrho%d2,       &
                             dnJac%d1,dnJac%d2,mole%ndimG,mole%nb_act)
          ELSE !old vep
            CALL sub_vep(vep,dnGG%d0,dnGG%d1,dnrho%d1,dnrho%d2,           &
                             dnJac%d1,dnJac%d2,mole%ndimG,mole%nb_act)
          END IF
        ELSE
          IF (new_vep) THEN !new vep
            CALL sub_vep_new(vep,dnGG_loc%d0,dnGG_loc%d1,dnrho%d1,dnrho%d2, &
                             dnJac%d1,dnJac%d2,mole%ndimG,mole%nb_act)
          ELSE !old vep
            CALL sub_vep(vep,dnGG_loc%d0,dnGG_loc%d1,dnrho%d1,dnrho%d2,      &
                             dnJac%d1,dnJac%d2,mole%ndimG,mole%nb_act)
          END IF
          CALL dealloc_dnSVM(dnGG_loc)
        END IF
        vep_done = .TRUE.

        IF (present(rho)) THEN 
          rho = dnrho%d0
          !write(out_unitp,*) 'rho :',rho
          rho_done = .TRUE.
        END IF
  
        CALL dealloc_dnSVM(dnJac)
        CALL dealloc_dnSVM(dnrho)
      END IF
    ELSE
      vep = ZERO
      vep_done = .TRUE.
    END IF

    IF (present(rho) .AND. .NOT. rho_done) THEN ! rho needs to calculate (nderiv_loc = 0). Vep calculation is done.

      !----------------------------------------------------
      !-- jac0,jaci,JACij calculation
      !----------------------------------------------------
      nderiv_loc = 0
      CALL alloc_dnSVM(dnJac,mole%nb_act,nderiv=nderiv_loc)
      CALL alloc_dnSVM(dnrho,mole%nb_act,nderiv=nderiv_loc)
      IF (present(dnGG)) THEN
        CALL sub3_dnDetGG(dnJac,dnGG,nderiv_loc,mole%masses,mole%Mtot_inv,mole%ncart)
      ELSE
        CALL alloc_dnSVM(dnGG_loc,mole%ndimG,mole%ndimG,mole%nb_act,nderiv=nderiv_loc)

        CALL get_dng_dnGG(Qact,para_Tnum,mole,dnGG=dnGG_loc,nderiv=nderiv_loc)
        CALL sub3_dnDetGG(dnJac,dnGG_loc,nderiv_loc,mole%masses,mole%Mtot_inv,mole%ncart)
        CALL dealloc_dnSVM(dnGG_loc)
      END IF
      IF (debug) write(out_unitp,*) 'dnJac'
      IF (debug) CALL write_dnS(dnJac)

      !----------------------------------------------------
      !-- f0,fi,Fij calculation
      !----------------------------------------------------
      CALL sub3_dnrho(dnrho,dnJac,Qact,mole,                          &
                      nderiv_loc,para_Tnum%num_x,para_Tnum%stepT,     &
                      para_Tnum%nrho)

      IF (debug) write(out_unitp,*) 'dnrho'
      IF (debug) CALL write_dnS(dnrho)
      rho = dnrho%d0
      rho_done = .TRUE.

      CALL dealloc_dnSVM(dnJac)
      CALL dealloc_dnSVM(dnrho)
    END IF

    !---------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'vep',vep
      write(out_unitp,*) 'END ',name_sub
      CALL flush_perso(out_unitp)
    END IF
    !-----------------------------------------------------------

END FUNCTION get_vep_rho
FUNCTION get_vepTaylor(Qact,mole,para_Tnum) RESULT (vep)
  IMPLICIT NONE

  real (kind=Rkind)                             :: vep
    !-----------------------------------------------------------------
    real (kind=Rkind), intent(in)               :: Qact(:)
    !----- for the CoordType and Tnum --------------------------------
    TYPE (CoordType),  intent(in)               :: mole
    TYPE (Tnum),       intent(in)               :: para_Tnum
    !-----------------------------------------------------------------


    !-----------------------------------------------------------------
    ! local variables
    integer                           :: nderiv_loc,i,j
    real (kind=Rkind), allocatable    :: DQ(:)

    !-----------------------------------------------------------------


    !--- for debuging --------------------------------------------------
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    character (len=*), parameter :: name_sub = 'get_vepTaylor'
    !---------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) 'vep_type   ',para_Tnum%vep_type
      write(out_unitp,*)
      write(out_unitp,*) 'Qact',Qact
      !write(out_unitp,*)
      !CALL Write_CoordType(mole)
      !write(out_unitp,*)
    END IF
    !---------------------------------------------------------

    CALL check_alloc_dnS(para_Tnum%dnVepref,'para_Tnum%dnVepref',name_sub)

    DQ = Qact - mole%ActiveTransfo%Qact0
    IF (debug) THEN 
      write(out_unitp,*) 'DQ',DQ(:)
      write(out_unitp,*) 'para_Tnum%dnVepref'
      CALL Write_dnSVM(para_Tnum%dnVepref)
    END IF

    vep = para_Tnum%dnVepref%d0

    IF (para_Tnum%vepTaylor_Order > 0) THEN
      DO i=1,mole%nb_act
        vep = vep + para_Tnum%dnVepref%d1(i)*DQ(i)
      END DO
    END IF

    IF (para_Tnum%vepTaylor_Order > 1) THEN
      DO i=1,mole%nb_act
      DO j=1,mole%nb_act
        vep = vep + HALF*para_Tnum%dnVepref%d2(j,i)*DQ(i)*DQ(j)
      END DO
      END DO
    END IF

    deallocate(DQ)

    !---------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'vep',vep
      write(out_unitp,*) 'END ',name_sub
      CALL flush_perso(out_unitp)
    END IF
    !-----------------------------------------------------------

END FUNCTION get_vepTaylor
SUBROUTINE Calc_vep_rho_from_dnGG(vep,rho,Qact,dnGG,mole,para_Tnum)
  USE mod_dnDetGG_dnDetg, only : sub3_dnDetGG
  IMPLICIT NONE


    !-----------------------------------------------------------------
    real (kind=Rkind), intent(inout)            :: vep,rho
    real (kind=Rkind), intent(in)               :: Qact(:)
    TYPE(Type_dnMat),  intent(in)               :: dnGG
    !----- for the CoordType and Tnum --------------------------------
    TYPE (CoordType),  intent(in)               :: mole
    TYPE (Tnum),       intent(in)               :: para_Tnum
    !-----------------------------------------------------------------


    !-----------------------------------------------------------------
    ! local variables
    TYPE(Type_dnS)    :: dnJac,dnrho
    integer           :: nderiv_loc
    !-----------------------------------------------------------------


    !--- for debuging --------------------------------------------------
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    character (len=*), parameter :: name_sub = 'Calc_vep_rho_from_dnGG'
    !---------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) 'nrho       ',para_Tnum%nrho
      write(out_unitp,*) 'vep_type   ',para_Tnum%vep_type
      write(out_unitp,*) 'dnGG%nderiv',dnGG%nderiv
      write(out_unitp,*) 'ndimG      ',mole%ndimG
      write(out_unitp,*)
      write(out_unitp,*) 'Qact',Qact

      write(out_unitp,*) 'dnGG'
      CALL write_dnMat(dnGG)
      !write(out_unitp,*)
      !CALL Write_CoordType(mole)
      !write(out_unitp,*)
      write(out_unitp,*)
      write(out_unitp,*) 'num_GG',para_Tnum%num_GG
      write(out_unitp,*) 'num_g ',para_Tnum%num_g
      write(out_unitp,*) 'num_x ',para_Tnum%num_x
      write(out_unitp,*) 'JJ',para_Tnum%JJ
      write(out_unitp,*)
    END IF
    !---------------------------------------------------------

    !------------------------------------------------------------
    !-- extra potential term ------------------------------------
    IF ( (para_Tnum%vep_type == 1 .OR. para_Tnum%vep_type == 100) .AND. &
                                                dnGG%nderiv == 2) THEN

      nderiv_loc = 2

      !----------------------------------------------------
      !-- jac0,jaci,JACij calculation
      !----------------------------------------------------
      CALL alloc_dnSVM(dnJac,dnGG%nb_var_deriv,nderiv=nderiv_loc)
      CALL alloc_dnSVM(dnrho,dnGG%nb_var_deriv,nderiv=nderiv_loc)

      CALL sub3_dnDetGG(dnJac,dnGG,nderiv_loc,                          &
                                 mole%masses,mole%Mtot_inv,mole%ncart)
      IF (debug) write(out_unitp,*) 'dnJac'
      IF (debug) CALL write_dnS(dnJac)
      !----------------------------------------------------
      !-- f0,fi,Fij calculation
      !----------------------------------------------------
      CALL sub3_dnrho(dnrho,dnJac,Qact,mole,                          &
                      nderiv_loc,para_Tnum%num_x,para_Tnum%stepT,     &
                      para_Tnum%nrho)
      rho = dnrho%d0
      !write(out_unitp,*) 'rho :',rho
      IF (debug) write(out_unitp,*) 'dnrho'
      IF (debug) CALL write_dnS(dnrho)
      !----------------------------------------------------
      !-- vep calculation
      !----------------------------------------------------
      IF (new_vep) THEN !new vep
        CALL sub_vep_new(vep,dnGG%d0,dnGG%d1,dnrho%d1,dnrho%d2,       &
                         dnJac%d1,dnJac%d2,mole%ndimG,mole%nb_act)
      ELSE !old vep
        CALL sub_vep(vep,dnGG%d0,dnGG%d1,dnrho%d1,dnrho%d2,           &
                         dnJac%d1,dnJac%d2,mole%ndimG,mole%nb_act)
      END IF

      CALL dealloc_dnSVM(dnJac)
      CALL dealloc_dnSVM(dnrho)
    ELSE ! here, the vpe=0, but rho needs to calculate rho (nderiv_loc = 0)
      nderiv_loc = 0

      !----------------------------------------------------
      !-- jac0,jaci,JACij calculation
      !----------------------------------------------------
      CALL alloc_dnSVM(dnJac,dnGG%nb_var_deriv,nderiv=nderiv_loc)
      CALL alloc_dnSVM(dnrho,dnGG%nb_var_deriv,nderiv=nderiv_loc)

      CALL sub3_dnDetGG(dnJac,dnGG,nderiv_loc,                          &
                                 mole%masses,mole%Mtot_inv,mole%ncart)
      !write(out_unitp,*) 'dnJac'
      !CALL write_dnS(dnJac)
      !----------------------------------------------------
      !-- f0,fi,Fij calculation
      !----------------------------------------------------
      CALL sub3_dnrho(dnrho,dnJac,Qact,mole,                          &
                      nderiv_loc,para_Tnum%num_x,para_Tnum%stepT,     &
                      para_Tnum%nrho)
      rho = dnrho%d0

      vep = ZERO

      CALL dealloc_dnSVM(dnJac)
      CALL dealloc_dnSVM(dnrho)

    END IF


    !---------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'vep',vep
      write(out_unitp,*) 'END ',name_sub
      CALL flush_perso(out_unitp)
    END IF
    !-----------------------------------------------------------

  END SUBROUTINE Calc_vep_rho_from_dnGG
  SUBROUTINE Calc_vepTalylor_rho_from_dnGG(vep,rho,Qact,dnGG,mole,para_Tnum)
    USE mod_dnDetGG_dnDetg, only : sub3_dnDetGG
    IMPLICIT NONE
  
  
      !-----------------------------------------------------------------
      real (kind=Rkind), intent(inout)            :: vep,rho
      real (kind=Rkind), intent(in)               :: Qact(:)
      TYPE(Type_dnMat),  intent(in)               :: dnGG
      !----- for the CoordType and Tnum --------------------------------
      TYPE (CoordType),  intent(in)               :: mole
      TYPE (Tnum),       intent(in)               :: para_Tnum
      !-----------------------------------------------------------------
  
  
      !-----------------------------------------------------------------
      ! local variables
      TYPE(Type_dnS)                    :: dnJac,dnrho
      integer                           :: nderiv_loc,i,j
      real (kind=Rkind), allocatable    :: DQ(:)

      !-----------------------------------------------------------------
  
  
      !--- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'Calc_vepTalylor_rho_from_dnGG'
      !---------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nrho       ',para_Tnum%nrho
        write(out_unitp,*) 'vep_type   ',para_Tnum%vep_type
        write(out_unitp,*) 'dnGG%nderiv',dnGG%nderiv
        write(out_unitp,*) 'ndimG      ',mole%ndimG
        write(out_unitp,*)
        write(out_unitp,*) 'Qact',Qact
  
        write(out_unitp,*) 'dnGG'
        CALL write_dnMat(dnGG)
        !write(out_unitp,*)
        !CALL Write_CoordType(mole)
        !write(out_unitp,*)
        write(out_unitp,*)
        write(out_unitp,*) 'num_GG',para_Tnum%num_GG
        write(out_unitp,*) 'num_g ',para_Tnum%num_g
        write(out_unitp,*) 'num_x ',para_Tnum%num_x
        write(out_unitp,*) 'JJ',para_Tnum%JJ
        write(out_unitp,*)
      END IF
      !---------------------------------------------------------

      CALL check_alloc_dnMat(dnGG,'dnGG',name_sub)
      CALL check_alloc_dnS(para_Tnum%dnVepref,'para_Tnum%dnVepref',name_sub)
  
      DQ = Qact - mole%ActiveTransfo%Qact0
      IF (debug) THEN 
        write(out_unitp,*) 'DQ',DQ(:)
        write(out_unitp,*) 'para_Tnum%dnVepref'
        CALL Write_dnSVM(para_Tnum%dnVepref)
      END IF
  
      vep = para_Tnum%dnVepref%d0
  
      IF (para_Tnum%vepTaylor_Order > 0) THEN
        DO i=1,mole%nb_act
          vep = vep + para_Tnum%dnVepref%d1(i)*DQ(i)
        END DO
      END IF
  
      IF (para_Tnum%vepTaylor_Order > 1) THEN
        DO i=1,mole%nb_act
        DO j=1,mole%nb_act
          vep = vep + HALF*para_Tnum%dnVepref%d2(j,i)*DQ(i)*DQ(j)
        END DO
        END DO
      END IF
  
  
      deallocate(DQ)

      !------------------------------------------------------------
      !-- rho only ------------------------------------
      nderiv_loc = 0
  
      !----------------------------------------------------
      !-- jac0,jaci,JACij calculation
      !----------------------------------------------------
      CALL alloc_dnSVM(dnJac,dnGG%nb_var_deriv,nderiv=nderiv_loc)
      CALL alloc_dnSVM(dnrho,dnGG%nb_var_deriv,nderiv=nderiv_loc)
  
      CALL sub3_dnDetGG(dnJac,dnGG,nderiv_loc,                          &
                                 mole%masses,mole%Mtot_inv,mole%ncart)
      !write(out_unitp,*) 'dnJac'
      !CALL write_dnS(dnJac)
      !----------------------------------------------------
      !-- f0,fi,Fij calculation
      !----------------------------------------------------
      CALL sub3_dnrho(dnrho,dnJac,Qact,mole,                          &
                      nderiv_loc,para_Tnum%num_x,para_Tnum%stepT,     &
                      para_Tnum%nrho)
      rho = dnrho%d0

      CALL dealloc_dnSVM(dnJac)
      CALL dealloc_dnSVM(dnrho)

      !---------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'vep',vep
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
      !-----------------------------------------------------------
  
  END SUBROUTINE Calc_vepTalylor_rho_from_dnGG

  SUBROUTINE Set_dnVepTaylor(dnVep,Qact,mole,para_Tnum,TaylorOrder)
    IMPLICIT NONE

    !-----------------------------------------------------------------
    TYPE(Type_dnS),    intent(inout)            :: dnVep
    real (kind=Rkind), intent(in)               :: Qact(:)
    integer,           intent(in)               :: TaylorOrder

    !----- for the CoordType and Tnum --------------------------------
    TYPE (CoordType),  intent(in)               :: mole
    TYPE (Tnum),       intent(inout)            :: para_Tnum
    !-----------------------------------------------------------------
  
  
      !-----------------------------------------------------------------
      ! local variables
      integer                           :: i,j
      real (kind=Rkind), allocatable    :: Qact_loc(:)
      real (kind=Rkind)                 :: vep,step2,step24,stepp
      !-----------------------------------------------------------------
  
  
      !--- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'Set_dnVepTaylor'
      !---------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'vepTaylor_Order   ',para_Tnum%vepTaylor_Order
        write(out_unitp,*)
        write(out_unitp,*) 'Qact',Qact
        flush(out_unitp)
      END IF
      !---------------------------------------------------------

      IF (TaylorOrder < 0) RETURN
      IF (para_Tnum%stepT == ZERO) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' para_Tnum%stepT is zero'
        STOP ' ERROR : para_Tnum%stepT is zero'
      END IF
      step2  = ONE/(para_Tnum%stepT*para_Tnum%stepT)
      step24 = step2/FOUR
      stepp  = ONE/(para_Tnum%stepT+para_Tnum%stepT)

      IF (.NOT. dnVep%alloc) CALL alloc_dnSVM(dnVep,mole%nb_act,TaylorOrder)

!===================================================================
!
!       Computation in Qact
!
!===================================================================
      Qact_loc = Qact(:)

      dnVep%d0 = get_vep_rho(Qact_loc,mole,para_Tnum,WithTaylor=.FALSE.)
 
!===================================================================
!
!       Computation in Qact(i)+step
!               and in Qact(i)-step
!
!===================================================================

      IF (TaylorOrder >=1) THEN
        DO i=1,mole%nb_act
          Qact_loc(:) = Qact(:)

          Qact_loc(i) = Qact(i) + para_Tnum%stepT
          vep = get_vep_rho(Qact_loc,mole,para_Tnum,WithTaylor=.FALSE.)

          dnVep%d1(i) = vep
          IF (TaylorOrder == 2) dnVep%d2(i,i) = vep

          Qact_loc(i) = Qact(i) - para_Tnum%stepT
          vep = get_vep_rho(Qact_loc,mole,para_Tnum,WithTaylor=.FALSE.)

          dnVep%d1(i) = ( dnVep%d1(i) - vep ) * stepp

          IF (TaylorOrder == 2) THEN 
            dnVep%d2(i,i) = ( dnVep%d2(i,i) + vep - TWO*dnVep%d0 ) * step2
          END IF

        END DO
      END IF

!===================================================================
!
!       Computation in Qact(i)+/-step
!               and in Qact(j)+/-step
!
!===================================================================

      IF (TaylorOrder ==2) THEN
        DO i=1,mole%nb_act
        DO j=i+1,mole%nb_act
          Qact_loc(:) = Qact(:)

          Qact_loc(i) = Qact(i) + para_Tnum%stepT
          Qact_loc(j) = Qact(j) + para_Tnum%stepT
          vep = get_vep_rho(Qact_loc,mole,para_Tnum,WithTaylor=.FALSE.)

          dnVep%d2(i,j) = vep

          Qact_loc(i) = Qact(i) - para_Tnum%stepT
          Qact_loc(j) = Qact(j) - para_Tnum%stepT
          vep = get_vep_rho(Qact_loc,mole,para_Tnum,WithTaylor=.FALSE.)

          dnVep%d2(i,j) = dnVep%d2(i,j) + vep

          Qact_loc(i) = Qact(i) - para_Tnum%stepT
          Qact_loc(j) = Qact(j) + para_Tnum%stepT
          vep = get_vep_rho(Qact_loc,mole,para_Tnum,WithTaylor=.FALSE.)

          dnVep%d2(i,j) = dnVep%d2(i,j) - vep

          Qact_loc(i) = Qact(i) + para_Tnum%stepT
          Qact_loc(j) = Qact(j) - para_Tnum%stepT
          vep = get_vep_rho(Qact_loc,mole,para_Tnum,WithTaylor=.FALSE.)

          dnVep%d2(i,j) = dnVep%d2(i,j) - vep

          !-- d2A/dQidQj -----------------------------------------
          dnVep%d2(i,j) = dnVep%d2(i,j) * step24
          dnVep%d2(j,i) = dnVep%d2(i,j)

        END DO
        END DO
      END IF

      !---------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'dnVep'
        CALL write_dnS(dnVep)
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
      !-----------------------------------------------------------
  
  END SUBROUTINE Set_dnVepTaylor


  SUBROUTINE Calc_vep_rho_from_dng_AND_dnGG(vep,rho,Qact,dng,dnGG,mole,para_Tnum)
  USE mod_dnDetGG_dnDetg, only : sub3_dndetA
  IMPLICIT NONE


    !-----------------------------------------------------------------
    real (kind=Rkind), intent(inout)            :: vep,rho
    real (kind=Rkind), intent(in)               :: Qact(:)
    TYPE(Type_dnMat),  intent(in)               :: dng,dnGG
    !----- for the CoordType and Tnum --------------------------------
    TYPE (CoordType),  intent(in)               :: mole
    TYPE (Tnum),       intent(in)               :: para_Tnum
    !-----------------------------------------------------------------


    !-----------------------------------------------------------------
    ! local variables
    TYPE(Type_dnS)    :: dnJac,dnrho
    integer           :: nderiv_loc
    !-----------------------------------------------------------------


    !--- for debuging --------------------------------------------------
    logical, parameter :: debug = .FALSE.
    !logical, parameter :: debug = .TRUE.
    character (len=*), parameter :: name_sub = 'Calc_vep_rho_from_dng_AND_dnGG'
    !---------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) 'nrho       ',para_Tnum%nrho
      write(out_unitp,*) 'vep_type   ',para_Tnum%vep_type
      write(out_unitp,*) 'dng%nderiv',dng%nderiv
      write(out_unitp,*) 'ndimG      ',mole%ndimG
      write(out_unitp,*)
      write(out_unitp,*) 'Qact',Qact
      write(out_unitp,*) 'dng'
      CALL write_dnMat(dng)
      !write(out_unitp,*)
      !CALL Write_CoordType(mole)
      !write(out_unitp,*)
      write(out_unitp,*)
      write(out_unitp,*) 'num_GG',para_Tnum%num_GG
      write(out_unitp,*) 'num_g ',para_Tnum%num_g
      write(out_unitp,*) 'num_x ',para_Tnum%num_x
      write(out_unitp,*) 'JJ',para_Tnum%JJ
      write(out_unitp,*)
    END IF
    !---------------------------------------------------------

    !------------------------------------------------------------
    !-- extra potential term ------------------------------------
    IF ( (para_Tnum%vep_type == 1 .OR. abs(para_Tnum%vep_type) == 100) .AND. &
                                                dng%nderiv == 2) THEN

      nderiv_loc = 2

      !----------------------------------------------------
      !-- jac0,jaci,JACij calculation
      !----------------------------------------------------
      CALL alloc_dnSVM(dnJac,dng%nb_var_deriv,nderiv=nderiv_loc)
      CALL alloc_dnSVM(dnrho,dng%nb_var_deriv,nderiv=nderiv_loc)

      CALL sub3_dndetA(dnJac,dng,nderiv_loc,                            &
                       mole%masses,mole%Mtot_inv,mole%ncart)

      IF (debug) write(out_unitp,*) 'dnJac'
      IF (debug) CALL write_dnS(dnJac)
      !----------------------------------------------------
      !-- f0,fi,Fij calculation
      !----------------------------------------------------
      CALL sub3_dnrho(dnrho,dnJac,Qact,mole,                          &
                      nderiv_loc,para_Tnum%num_x,para_Tnum%stepT,     &
                      para_Tnum%nrho)
      rho = dnrho%d0
      !write(out_unitp,*) 'rho :',rho
      IF (debug) write(out_unitp,*) 'dnrho'
      IF (debug) CALL write_dnS(dnrho)
      !----------------------------------------------------
      !-- vep calculation
      !----------------------------------------------------
      IF (new_vep) THEN !new vep
        CALL sub_vep_new(vep,dnGG%d0,dnGG%d1,dnrho%d1,dnrho%d2,       &
                         dnJac%d1,dnJac%d2,mole%ndimG,mole%nb_act)
      ELSE !old vep
        CALL sub_vep(vep,dnGG%d0,dnGG%d1,dnrho%d1,dnrho%d2,           &
                         dnJac%d1,dnJac%d2,mole%ndimG,mole%nb_act)
      END IF

      CALL dealloc_dnSVM(dnJac)
      CALL dealloc_dnSVM(dnrho)
    ELSE ! here, the vpe=0, but rho needs to be calculate (nderiv_loc = 0)
      nderiv_loc = 0

      !----------------------------------------------------
      !-- jac0,jaci,JACij calculation
      !----------------------------------------------------
      CALL alloc_dnSVM(dnJac,dng%nb_var_deriv,nderiv=nderiv_loc)
      CALL alloc_dnSVM(dnrho,dng%nb_var_deriv,nderiv=nderiv_loc)

      CALL sub3_dndetA(dnJac,dng,nderiv_loc,                            &
                       mole%masses,mole%Mtot_inv,mole%ncart)

      !write(out_unitp,*) 'dnJac'
      !CALL write_dnS(dnJac)
      !----------------------------------------------------
      !-- f0,fi,Fij calculation
      !----------------------------------------------------
      CALL sub3_dnrho(dnrho,dnJac,Qact,mole,                          &
                      nderiv_loc,para_Tnum%num_x,para_Tnum%stepT,     &
                      para_Tnum%nrho)
      rho = dnrho%d0

      vep = ZERO

      CALL dealloc_dnSVM(dnJac)
      CALL dealloc_dnSVM(dnrho)

    END IF


    !---------------------------------------------------------
    IF (debug) THEN
      write(out_unitp,*) 'vep',vep
      write(out_unitp,*) 'END ',name_sub
      CALL flush_perso(out_unitp)
    END IF
    !-----------------------------------------------------------

  END SUBROUTINE Calc_vep_rho_from_dng_AND_dnGG

      SUBROUTINE sub_vep_new(vep,d0invA,d1invA,                         &
                             fi,Fij,jaci,JACij,ndimA,nb_act)
      IMPLICIT NONE

       integer ndimA,nb_act
       real (kind=Rkind) ::  d0invA(ndimA,ndimA)
       real (kind=Rkind) ::  d1invA(ndimA,ndimA,nb_act)
       real (kind=Rkind) ::  fi(nb_act)
       real (kind=Rkind) ::  Fij(nb_act,nb_act)
       real (kind=Rkind) ::  jaci(nb_act)
       real (kind=Rkind) ::  JACij(nb_act,nb_act)
       real (kind=Rkind) ::  vep

       integer i,j
!     -- extra potential term ------------------------------------
      vep = ZERO
      DO i=1,nb_act
        DO j=1,nb_act
          vep = vep + d1invA(i,j,i)/FOUR * ( jaci(j) -fi(j) )        +  &
                  d0invA(i,j)/FOUR   * ( JACij(i,j)-Fij(i,j)         +  &
                                HALF * ( fi(i)*fi(j)-jaci(i)*jaci(j) +  &
                                         fi(i)*jaci(j)-jaci(i)*fi(j) ) )
      END DO
      END DO

!     write(out_unitp,*) 'd0invA',d0invA
!     write(out_unitp,*) 'd1invA',d1invA
!     write(out_unitp,*) 'fi,jaci',fi,jaci
!     write(out_unitp,*) 'Fij,JACij',Fij,JACij
!     write(out_unitp,*) 'vep',vep
!     ------------------------------------------------------------
      end subroutine sub_vep_new


      SUBROUTINE sub_vep_2(vep,d0invA,d1invA,                           &
                           fi,Fij,jaci,JACij,ndimA,nb_act)
      IMPLICIT NONE

       integer ndimA,nb_act
       real (kind=Rkind) ::  d0invA(ndimA,ndimA)
       real (kind=Rkind) ::  d1invA(ndimA,ndimA)
       real (kind=Rkind) ::  fi(nb_act)
       real (kind=Rkind) ::  Fij(nb_act,nb_act)
       real (kind=Rkind) ::  jaci(nb_act)
       real (kind=Rkind) ::  JACij(nb_act,nb_act)
       real (kind=Rkind) ::  vep

       integer i,j

!     -- extra potential term ------------------------------------
      vep = ZERO
      DO i=1,nb_act
        DO j=1,nb_act
        vep = vep + d1invA(i,j)*  (fi(j)-jaci(j))        +              &
                    d0invA(i,j)*( (Fij(i,j)-JACij(i,j))  -              &
                          HALF*(fi(i)*fi(j)-jaci(i)*jaci(j)) )
      END DO
      END DO
      vep = -vep/FOUR
!     ------------------------------------------------------------
      end subroutine sub_vep_2

      SUBROUTINE sub_vep(vep,d0invA,d1invA,                             &
                         fi,Fij,jaci,JACij,ndimA,nb_act)
      IMPLICIT NONE

       integer ndimA,nb_act
       real (kind=Rkind) ::  d0invA(ndimA,ndimA)
       real (kind=Rkind) ::  d1invA(ndimA,ndimA,nb_act)
       real (kind=Rkind) ::  fi(nb_act)
       real (kind=Rkind) ::  Fij(nb_act,nb_act)
       real (kind=Rkind) ::  jaci(nb_act)
       real (kind=Rkind) ::  JACij(nb_act,nb_act)
       real (kind=Rkind) ::  vep

       integer i,j
!     -- extra potential term ------------------------------------

      vep = ZERO
      DO i=1,nb_act
        DO j=1,nb_act
          vep = vep + d1invA(i,j,i)*( fi(j)-jaci(j))        +           &
                      d0invA(i,j)*  ( Fij(i,j)-JACij(i,j)   -           &
                             HALF*(fi(i)*fi(j)-jaci(i)*jaci(j)) )

      END DO
      END DO
      vep = -vep/FOUR

!     write(out_unitp,*) 'd0invA',d0invA
!     write(out_unitp,*) 'd1invA',d1invA
!     write(out_unitp,*) 'fi,jaci',fi,jaci
!     write(out_unitp,*) 'Fij,JACij',Fij,JACij
!     write(out_unitp,*) 'vep',vep
!     ------------------------------------------------------------
      end subroutine sub_vep


!======================================================================
!
!      Transfert G100 (g or GG) to G (g or GG)
!       G100 includes coordinates of type 100
!
!
!======================================================================
      SUBROUTINE dngG100_TO_dngG(dnG100,dnG,mole100,mole)
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType), intent(in)    :: mole,mole100
      TYPE(Type_dnMat), intent(in)    :: dnG100
      TYPE(Type_dnMat), intent(inout) :: dnG

      integer :: i,j


      IF (mole100%ndimG /= mole%ndimG+mole%nb_rigid100) THEN
        write(out_unitp,*) ' ERROR in dngG100_TO_dngG'
        write(out_unitp,*) ' mole100%ndimG NOT = mole%ndimG+mole%nb_rigid100'
        write(out_unitp,*) mole100%ndimG,mole%ndimG,mole%nb_rigid100
        STOP
      END IF
      IF (dnG100%nderiv /= dnG%nderiv) THEN
        write(out_unitp,*) ' ERROR in dngG100_TO_dngG'
        write(out_unitp,*) ' dnG100%nderiv NOT = dnG%nderiv'
        write(out_unitp,*) dnG100%nderiv,dnG%nderiv
        STOP
      END IF

!     - transfo G100 => G
      CALL gG100TOgG(dnG100%d0,dnG%d0,mole100,mole)

      IF (dnG%nderiv > 0) THEN
        DO i=1,mole%nb_act
          CALL gG100TOgG(dnG100%d1(:,:,i),dnG%d1(:,:,i),mole100,mole)
        END DO
      END IF

      IF (dnG%nderiv > 1) THEN
        DO i=1,mole%nb_act
        DO j=1,mole%nb_act
         CALL gG100TOgG(dnG100%d2(:,:,i,j),dnG%d2(:,:,i,j),mole100,mole)
        END DO
        END DO

      END IF

      END SUBROUTINE dngG100_TO_dngG
      SUBROUTINE gG100TOgG(d0G100,d0G,mole100,mole)
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType),  intent(in)    :: mole,mole100
      real (kind=Rkind), intent(inout) :: d0G(mole%ndimG,mole%ndimG)
      real (kind=Rkind), intent(in)    :: d0G100(mole100%ndimG,mole100%ndimG)

      integer :: i,j,i100,j100

      IF (mole100%ndimG /= mole%ndimG+mole%nb_rigid100) THEN
        write(out_unitp,*) ' ERROR in gG100TOgG'
        write(out_unitp,*) ' mole100%ndimG NOT = mole%ndimG+mole%nb_rigid100'
        write(out_unitp,*) mole100%ndimG,mole%ndimG,mole%nb_rigid100
        STOP
      END IF


      DO i=1,mole%ndimG
      DO j=1,mole%ndimG
        i100 = i
        j100 = j
        IF (i > mole%nb_act) i100 = i + mole%nb_rigid100
        IF (j > mole%nb_act) j100 = j + mole%nb_rigid100
        !IF (j==1) write(6,*) 'gG100TOgG',i,i100
        d0G(i,j)  = d0G100(i100,j100)
      END DO
      END DO

      END SUBROUTINE gG100TOgG
      SUBROUTINE Linear_Dep(g)
      IMPLICIT NONE
      real (kind=Rkind), intent(in) :: g(:,:)

      integer :: n
      real (kind=Rkind), allocatable :: Vec(:,:),vp(:)
      real (kind=Rkind), allocatable :: gij(:,:),gsi(:),d1f(:),gnsi(:),gn(:,:)
      integer :: i,j

RETURN
      n = size(g,dim=1)
      allocate(Vec(n,n))
      allocate(vp(n))

      CALL Write_Mat(g,out_unitp,5,info='g')
      CALL diagonalization(g,vp,Vec,n,2,1,.FALSE.)
      CALL Write_Mat(Vec,out_unitp,5,info='Vec')
      CALL Write_Vec(vp,out_unitp,5,info='vp')
      write(out_unitp,*)
      CALL flush_perso(out_unitp)

      gij = g(1:n-1,1:n-1)
      gsi = g(1:n-1,n)
      CALL Write_Vec(gsi,out_unitp,5,info='gsi')
      write(out_unitp,*)

      allocate(d1f(n-1))

      CALL Linear_Sys(gij,gsi,d1f,n-1)
      gnsi = gsi(:) - matmul(gij,d1f)

      CALL Write_Vec(gnsi,out_unitp,5,info='gnsi')
      write(out_unitp,*)
      CALL flush_perso(out_unitp)


      gn = g
      gn(n,n) = gn(n,n)-TWO*dot_product(gsi,d1f)
      DO i=1,n-1
      DO j=1,n-1
        gn(n,n) = g(n,n) + gij(i,j)*d1f(i)*d1f(j)
      END DO
      END DO
      gn(n,1:n-1) = gsi(:) - matmul(gij,d1f) ! zero
      gn(1:n-1,n) = gn(n,1:n-1)               ! zero
      CALL Write_Mat(gn,out_unitp,5,info='gn')

      CALL diagonalization(gn,vp,Vec,n,2,1,.FALSE.)
      CALL Write_Mat(Vec,out_unitp,5,info='Vec')
      CALL Write_Vec(vp,out_unitp,5,info='vp')
      CALL flush_perso(out_unitp)


      deallocate(Vec)
      deallocate(vp)
      STOP
      END SUBROUTINE Linear_Dep
END MODULE mod_dnGG_dng
