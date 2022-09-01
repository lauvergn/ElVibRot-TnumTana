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
      MODULE mod_GWP
      USE mod_system
      USE mod_poly
      IMPLICIT NONE

        TYPE para_GWP
        !----------------------------------------------------------------
        !- enable to check if CALL init0_GWP() has been done -------------
          logical                     :: init0,notinit0

          logical          :: cplx
          integer          :: ndim  ! number of degrees of freedom
          TYPE (para_poly) :: poly
          real (kind=Rkind),    pointer :: Qmean(:)
          real (kind=Rkind),    pointer :: Pmean(:)
          real (kind=Rkind),    pointer :: RAmean(:,:)
          real (kind=Rkind)             :: Rphase
          complex (kind=Rkind), pointer :: CAmean(:,:)
          complex (kind=Rkind)          :: Cphase

          logical :: linearization
          logical :: trajectory
          real (kind=Rkind),    pointer :: RZ(:,:),RPZ(:,:)
          complex (kind=Rkind), pointer :: CZ(:,:),CPZ(:,:)

          logical :: alloc_GWP
        END TYPE para_GWP

        TYPE para_LHA
        !----------------------------------------------------------------
        !- enable to check if CALL init0_LHA() has been done -------------
          logical                     :: init0,notinit0
          logical :: alloc
          integer :: ndim  ! number of degrees of freedom

          real (kind=Rkind)             :: d0Ene,pot0
          real (kind=Rkind),    pointer :: d1Ene(:)
          real (kind=Rkind),    pointer :: d2Ene(:,:)

          real (kind=Rkind),    pointer :: d0G(:,:)
          real (kind=Rkind),    pointer :: d1G(:,:,:)
          real (kind=Rkind),    pointer :: d2G(:,:,:,:)

          real (kind=Rkind)             :: d0mu(3)
          real (kind=Rkind),    pointer :: d1mu(:,:)

        END TYPE para_LHA
        CONTAINS

!       ==============================================================
!         init0_GWP    : initalization
!         alloc_GWP    : allocation
!         dealloc_GWP  : deallocation
!         GWP1_TO_GWP2 : GWP2 = GWP1
!       ==============================================================

        SUBROUTINE init0_GWP(GWP)
          TYPE (para_GWP), intent(inout) :: GWP
          GWP%init0      = .TRUE.
          GWP%notinit0   = .FALSE.

          GWP%cplx           = .FALSE.
          GWP%linearization  = .FALSE.
          GWP%trajectory     = .FALSE.
          GWP%ndim           = 0

          GWP%alloc_GWP  = .FALSE.

          GWP%Cphase = CZERO
          GWP%Rphase =  ZERO
          nullify(GWP%Qmean)
          nullify(GWP%Pmean)
          nullify(GWP%RAmean)
          nullify(GWP%CAmean)
          nullify(GWP%RZ)
          nullify(GWP%RPZ)
          nullify(GWP%CZ)
          nullify(GWP%CPZ)

        END SUBROUTINE init0_GWP
!================================================================
!
!     check if init0 has been done
!
!================================================================
      SUBROUTINE check_init0_GWP(A,name_A,name_sub)
        TYPE (para_GWP), intent(in) :: A
        character (len=*), intent(in) :: name_A
        character (len=*), intent(in) :: name_sub

        IF ( (A%init0 .EQV. A%notinit0) .OR.                            &
             (A%notinit0 .AND. .NOT. A%init0) ) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) name_A,' has NOT been initiated with "init0_"'
          write(out_unitp,*) ' CHECK the source!!!!!'
          STOP
        END IF
      END SUBROUTINE check_init0_GWP

        SUBROUTINE alloc_GWP(GWP,ndim,cplx,linearization,trajectory)
           TYPE (para_GWP), intent(inout)      :: GWP
           integer, intent(in)                 :: ndim
           logical, intent(in), optional       :: cplx,linearization
           logical, intent(in), optional       :: trajectory

           integer :: err_mem,memory

           CALL check_init0_GWP(GWP,'GWP','alloc_GWP')

           IF (ndim <=0) THEN
             write(out_unitp,*) ' ERROR in alloc_GWP'
             write(out_unitp,*) ' ndim MUST be > 0',ndim
             STOP
           END IF
           GWP%ndim = ndim

           IF (present(linearization) ) THEN
             GWP%linearization = linearization
           ELSE
             GWP%linearization = .FALSE.
           END IF

           IF (present(cplx) ) THEN
             GWP%cplx = cplx
           ELSE
             GWP%cplx = .FALSE.
           END IF

           IF (present(trajectory) ) THEN
             GWP%trajectory = trajectory
           ELSE
             GWP%trajectory = .FALSE.
           END IF

           CALL alloc_array(GWP%Qmean,[ndim],"GWP%Qmean","alloc_GWP")
           GWP%Qmean(:) = ZERO
           CALL alloc_array(GWP%Pmean,[ndim],"GWP%Pmean","alloc_GWP")
           GWP%Pmean(:) = ZERO

           IF (GWP%cplx .AND. .NOT. GWP%trajectory) THEN
             CALL alloc_array(GWP%CAmean,[ndim,ndim],"GWP%CAmean","alloc_GWP")
             GWP%CAmean(:,:) = CZERO
             IF (GWP%linearization) THEN
               CALL alloc_array(GWP%CZ,[ndim,ndim],"GWP%CZ","alloc_GWP")
               GWP%CZ(:,:) = CZERO
               CALL alloc_array(GWP%CPZ,[ndim,ndim],"GWP%CPZ","alloc_GWP")
               GWP%CPZ(:,:) = CZERO
             END IF
           ELSE IF (.NOT. GWP%cplx .AND. .NOT. GWP%trajectory) THEN
             CALL alloc_array(GWP%RAmean,[ndim,ndim],"GWP%RAmean","alloc_GWP")
             GWP%RAmean(:,:) = ZERO
             IF (GWP%linearization) THEN
               CALL alloc_array(GWP%RZ,[ndim,ndim],"GWP%RZ","alloc_GWP")
               GWP%RZ(:,:) = ZERO
               CALL alloc_array(GWP%RPZ,[ndim,ndim],"GWP%RPZ","alloc_GWP")
               GWP%RPZ(:,:) = ZERO
             END IF
           END IF

           GWP%alloc_GWP = .TRUE.

        END SUBROUTINE alloc_GWP

        SUBROUTINE dealloc_GWP(GWP)
           TYPE (para_GWP), intent(inout)      :: GWP
           integer :: err_mem,memory

           CALL check_init0_GWP(GWP,'GWP','dealloc_GWP')

           IF (associated(GWP%CAmean))  THEN
             CALL dealloc_array(GWP%CAmean,"GWP%CAmean","dealloc_GWP")
           END IF
           IF (associated(GWP%CZ))      THEN
             CALL dealloc_array(GWP%CZ,"GWP%CZ","dealloc_GWP")
           END IF
           IF (associated(GWP%CPZ))     THEN
             CALL dealloc_array(GWP%CPZ,"GWP%CPZ","dealloc_GWP")
           END IF

           IF (associated(GWP%RAmean))  THEN
             CALL dealloc_array(GWP%RAmean,"GWP%RAmean","dealloc_GWP")
           END IF
           IF (associated(GWP%RZ))      THEN
             CALL dealloc_array(GWP%RZ,"GWP%RZ","dealloc_GWP")
           END IF
           IF (associated(GWP%RPZ))     THEN
             CALL dealloc_array(GWP%RPZ,"GWP%RPZ","dealloc_GWP")
           END IF

           IF (associated(GWP%Qmean))   THEN
             CALL dealloc_array(GWP%Qmean,"GWP%Qmean","dealloc_GWP")
           END IF
           IF (associated(GWP%Pmean))   THEN
             CALL dealloc_array(GWP%Pmean,"GWP%Pmean","dealloc_GWP")
           END IF
           GWP%alloc_GWP = .FALSE.

           CALL init0_GWP(GWP)

        END SUBROUTINE dealloc_GWP

        SUBROUTINE write_GWP(GWP)
           TYPE (para_GWP), intent(in)      :: GWP
           integer :: ndim


           CALL check_init0_GWP(GWP,'GWP','write_GWP')

           IF (GWP%alloc_GWP) THEN
             write(out_unitp,*) '---------------------'
             write(out_unitp,*) '-- write_GWP --------'
             write(out_unitp,*) 'init0,alloc_GWP',GWP%init0,GWP%alloc_GWP
             write(out_unitp,*) 'cplx,linearization',GWP%cplx,GWP%linearization
             write(out_unitp,*) 'trajectory',GWP%trajectory
             write(out_unitp,*) 'ndim',GWP%ndim
             ndim = GWP%ndim
             write(out_unitp,*) 'Qmean'
             CALL Write_Vec(GWP%Qmean,out_unitp,5)
             write(out_unitp,*) 'Pmean'
             CALL Write_Vec(GWP%Pmean,out_unitp,5)
             IF (GWP%cplx .AND. .NOT. GWP%trajectory) THEN
               write(out_unitp,*) 'phase',GWP%Cphase
               write(out_unitp,*) 'Amean'
               CALL Write_Mat(GWP%CAmean,out_unitp,5)
               IF (GWP%linearization) THEN
                 write(out_unitp,*) 'Z'
                 CALL Write_Mat(GWP%CZ,out_unitp,5)
                 write(out_unitp,*) 'PZ'
                 CALL Write_Mat(GWP%CPZ,out_unitp,5)
               END IF
             ELSE IF (.NOT. GWP%cplx .AND. .NOT. GWP%trajectory) THEN
               write(out_unitp,*) 'phase',GWP%Rphase
               write(out_unitp,*) 'Amean'
               CALL Write_Mat(GWP%RAmean,out_unitp,5)
               IF (GWP%linearization) THEN
                 write(out_unitp,*) 'Z'
                 CALL Write_Mat(GWP%RZ,out_unitp,5)
                 write(out_unitp,*) 'PZ'
                 CALL Write_Mat(GWP%RPZ,out_unitp,5)
               END IF
             END IF
           ELSE
             write(out_unitp,*) '---------------------'
             write(out_unitp,*) ' The GWP is NOT allocated : no write!'
             write(out_unitp,*) 'init0,alloc_GWP',GWP%init0,GWP%alloc_GWP
             write(out_unitp,*) '---------------------'
           END IF

        END SUBROUTINE write_GWP

        SUBROUTINE GWP1_TO_GWP2(GWP1,GWP2)
           TYPE (para_GWP), intent(in)      :: GWP1
           TYPE (para_GWP), intent(inout)   :: GWP2


           CALL check_init0_GWP(GWP1,'GWP1','GWP1_TO_GWP2')
           CALL check_init0_GWP(GWP2,'GWP2','GWP1_TO_GWP2')

           IF (GWP2%alloc_GWP) CALL dealloc_GWP(GWP2)

           CALL alloc_GWP(GWP2,ndim=GWP1%ndim,                          &
                          cplx=GWP1%cplx,                               &
                          linearization=GWP1%linearization,             &
                          trajectory=GWP1%trajectory)

           GWP2%Qmean = GWP1%Qmean
           GWP2%Pmean = GWP1%Pmean

           IF (GWP2%cplx .AND. .NOT. GWP2%trajectory) THEN
             GWP2%Cphase = GWP1%Cphase
             GWP2%CAmean = GWP1%CAmean
             IF (GWP2%linearization) THEN
               GWP2%CZ = GWP1%CZ
               GWP2%CPZ = GWP1%CPZ
             END IF
           ELSE IF (.NOT. GWP2%cplx .AND. .NOT. GWP2%trajectory) THEN
             GWP2%Rphase = GWP1%Rphase
             GWP2%RAmean = GWP1%RAmean
             IF (GWP2%linearization) THEN
               GWP2%RZ = GWP1%RZ
               GWP2%RPZ = GWP1%RPZ
             END IF
           END IF
        END SUBROUTINE GWP1_TO_GWP2

        SUBROUTINE init0_LHA(param_LHA)
          TYPE (para_LHA), intent(inout) :: param_LHA

          param_LHA%init0   = .TRUE.
          param_LHA%notinit0= .FALSE.

          param_LHA%ndim    = 0
          param_LHA%init0   = .TRUE.
          param_LHA%alloc   = .FALSE.

          param_LHA%pot0    = -huge(ONE)
          param_LHA%d0Ene   = ZERO
          nullify(param_LHA%d1Ene)
          nullify(param_LHA%d2Ene)

          nullify(param_LHA%d0G)
          nullify(param_LHA%d1G)
          nullify(param_LHA%d2G)

          param_LHA%d0mu(:) = ZERO
          nullify(param_LHA%d1mu)

        END SUBROUTINE init0_LHA
!================================================================
!
!     check if init0 has been done
!
!================================================================
      SUBROUTINE check_init0_LHA(A,name_A,name_sub)
        TYPE (para_LHA), intent(inout) :: A
        character (len=*), intent(in) :: name_A
        character (len=*), intent(in) :: name_sub

        IF ( (A%init0 .EQV. A%notinit0) .OR.                            &
             (A%notinit0 .AND. .NOT. A%init0) ) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) name_A,' has NOT been initiated with "init0_"'
          write(out_unitp,*) ' CHECK the source!!!!!'
          STOP
        END IF
      END SUBROUTINE check_init0_LHA

        SUBROUTINE alloc_LHA(param_LHA,ndim)
          TYPE (para_LHA), intent(inout) :: param_LHA
          integer, intent(in)            :: ndim
          integer :: err_mem,memory

          CALL check_init0_LHA(param_LHA,'param_LHA','alloc_LHA')

          IF (ndim <=0) THEN
            write(out_unitp,*) ' ERROR in alloc_param_LHA'
            write(out_unitp,*) ' ndim MUST be > 0',ndim
            STOP
          END IF
          param_LHA%ndim = ndim

          param_LHA%pot0  = -huge(ONE)
          param_LHA%d0Ene = ZERO
          CALL alloc_array(param_LHA%d1Ene,[ndim],                    &
                          "param_LHA%d1Ene","alloc_LHA")
          param_LHA%d1Ene(:) = ZERO
          CALL alloc_array(param_LHA%d2Ene,[ndim,ndim],               &
                          "param_LHA%d2Ene","alloc_LHA")
          param_LHA%d2Ene(:,:) = ZERO

          CALL alloc_array(param_LHA%d0G,[ndim,ndim],                 &
                          "param_LHA%d0G","alloc_LHA")
          param_LHA%d0G(:,:) = ZERO
          CALL alloc_array(param_LHA%d1G,[ndim,ndim,ndim],            &
                          "param_LHA%d1G","alloc_LHA")
          param_LHA%d1G(:,:,:) = ZERO
          CALL alloc_array(param_LHA%d2G,[ndim,ndim,ndim,ndim],       &
                          "param_LHA%d2G","alloc_LHA")
          param_LHA%d2G(:,:,:,:) = ZERO

          param_LHA%d0mu(:) = ZERO
          CALL alloc_array(param_LHA%d1mu,[3,ndim],                   &
                          "param_LHA%d1mu","alloc_LHA")
          param_LHA%d1mu(:,:) = ZERO

          param_LHA%alloc = .TRUE.

        END SUBROUTINE alloc_LHA

        SUBROUTINE dealloc_LHA(param_LHA)
          TYPE (para_LHA), intent(inout)      :: param_LHA
          integer :: err_mem,memory

          CALL check_init0_LHA(param_LHA,'param_LHA','dealloc_LHA')

          IF (param_LHA%alloc) THEN
            CALL dealloc_array(param_LHA%d1Ene,"param_LHA%d1Ene","dealloc_LHA")
            CALL dealloc_array(param_LHA%d2Ene,"param_LHA%d2Ene","dealloc_LHA")
            CALL dealloc_array(param_LHA%d0G,"param_LHA%d0G","dealloc_LHA")
            CALL dealloc_array(param_LHA%d1G,"param_LHA%d1G","dealloc_LHA")
            CALL dealloc_array(param_LHA%d2G,"param_LHA%d2G","dealloc_LHA")
            CALL dealloc_array(param_LHA%d1mu,"param_LHA%d1mu","dealloc_LHA")
          END IF

          CALL init0_LHA(param_LHA)
        END SUBROUTINE dealloc_LHA

      END MODULE mod_GWP

