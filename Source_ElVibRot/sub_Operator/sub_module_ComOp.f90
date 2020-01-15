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

      MODULE mod_ComOp

      USE mod_system
      USE mod_basis_BtoG_GtoB_SGType4, only : Type_SmolyakRep, dealloc_SmolyakRep

      IMPLICIT NONE

        TYPE param_ComOp

          integer                       :: nb_act1              = 0
          integer                       :: nb_ba                = 0
          integer                       :: nb_bi                = 0
          integer                       :: nb_be                = 0
          integer                       :: nb_bie               = 0

          integer                       :: MaxNormIndexBasisRep = huge(1)
          integer, pointer              :: NormIndexBasisRep(:) => null()  ! NormIndexBasisRep(nb_ba)


          logical                       :: ADA              = .FALSE. ! (t) => true adiabatic channels (not used yet)
          logical                       :: contrac_ba_ON_HAC= .FALSE. ! (t) => the active basis are contracted on each adiabatic channel
          integer                       :: max_nb_ba_ON_HAC = huge(1) ! largest dimension of d0Cba_ON_HAC for all adiabatic channels
          real (kind=Rkind)             :: max_ene_ON_HAC   = huge(ONE)  ! energy maximal

          real (kind=Rkind), pointer    :: d0Cba_ON_HAC(:,:,:) => null() ! d0Cba_ON_HAC(nb_ba,nb_ba,nb_bie)
          real (kind=Rkind), pointer    :: Eneba_ON_HAC(:,:)   => null() ! Eneba_ON_HAC(nb_ba,nb_bie)
          integer, pointer              :: nb_ba_ON_HAC(:)     => null() ! dimension of nb_ba_ON_HAC for a given i_h (nb_ba_ON_HAC(nb_bie))
          real (kind=Rkind), pointer    :: d0C_bhe_ini(:,:)    => null() ! d0C_bhe_ini(nb_bie,nb_bie) for ADA method

          !file name HADA
          TYPE (param_file)             :: file_HADA
          logical                       :: calc_grid_HADA   = .FALSE. ! (F). If T the grid has been calculated
          logical                       :: save_file_HADA   = .FALSE. ! save the HADA file when direct=0,1,3

          real (kind=Rkind), allocatable:: sqRhoOVERJac(:)            ! sqrt(rho/jac), nb_qa pts
          real (kind=Rkind), allocatable:: Jac(:)                     ! jac, nb_qa pts
          TYPE(Type_SmolyakRep)         :: SRep_sqRhoOVERJac          ! sqrt(rho/jac) in Smolyap Rep.
          TYPE(Type_SmolyakRep)         :: SRep_Jac                   ! jac in Smolyap Rep.


          !- the vectors are never allocated, but they point to other vectors
          integer                       :: nb_vp_spec       = huge(1) ! number of eigenvectors for the spectral represention
          integer, pointer              :: liste_spec(:)    => null() ! liste of eigenvectors (nb_vp_spec)
                                                                      ! default : 1,2,3,..... nb_vp_spec
          real (kind=Rkind), pointer    :: Rvp_spec(:,:)    => null() ! eigenvectors for the spectral represention
          complex (kind=Rkind), pointer :: Cvp_spec(:,:)    => null() ! Both are TRUE pointer

          real (kind=Rkind)             :: ZPE              = huge(ONE)
          logical                       :: Set_ZPE          = .FALSE.

        END TYPE param_ComOp
        !====================================================================

      CONTAINS

      SUBROUTINE dealloc_ComOp(ComOp,keep_init)

      TYPE (param_ComOp),          intent(inout) :: ComOp
      logical,           optional, intent(in)    :: keep_init

      logical :: keep_init_loc

      character (len=*), parameter :: name_sub='dealloc_ComOp'

      keep_init_loc = .FALSE.
      IF (present(keep_init)) keep_init_loc = keep_init

      IF (.NOT. keep_init_loc) THEN
        ComOp%nb_act1     = 0
        ComOp%nb_ba       = 0
        ComOp%nb_bi       = 0
        ComOp%nb_be       = 0
        ComOp%nb_bie      = 0
      END IF

      IF (associated(ComOp%NormIndexBasisRep))  THEN
        CALL dealloc_array(ComOp%NormIndexBasisRep,                     &
                          "ComOp%NormIndexBasisRep",name_sub)
      END IF
      ComOp%MaxNormIndexBasisRep     = huge(1)

      IF (.NOT. keep_init_loc) THEN
        ComOp%contrac_ba_ON_HAC = .FALSE.
        ComOp%max_ene_ON_HAC    = huge(ONE)  ! maximal energy
      END IF

      IF (associated(ComOp%d0Cba_ON_HAC))  THEN
        CALL dealloc_array(ComOp%d0Cba_ON_HAC,                          &
                          "ComOp%d0Cba_ON_HAC",name_sub)
      END IF
      IF (associated(ComOp%Eneba_ON_HAC))  THEN
        CALL dealloc_array(ComOp%Eneba_ON_HAC,                          &
                          "ComOp%Eneba_ON_HAC",name_sub)
      END IF
      IF (.NOT. keep_init_loc .AND. .NOT. ComOp%contrac_ba_ON_HAC) THEN
        IF (associated(ComOp%nb_ba_ON_HAC))  THEN
          CALL dealloc_array(ComOp%nb_ba_ON_HAC,                        &
                            "ComOp%nb_ba_ON_HAC",name_sub)
        END IF
        ComOp%max_nb_ba_ON_HAC = huge(1)
      END IF

      IF (associated(ComOp%d0C_bhe_ini))  THEN
        CALL dealloc_array(ComOp%d0C_bhe_ini,                           &
                          "ComOp%d0C_bhe_ini",name_sub)
      END IF

      IF (.NOT. keep_init_loc) THEN
        ComOp%ADA                 = .FALSE.

        ComOp%file_HADA%name      = make_FileName('SH_HADA')
        ComOp%file_HADA%unit      = 0
        ComOp%file_HADA%formatted = .TRUE.
        ComOp%file_HADA%append    = .FALSE.
        ComOp%file_HADA%old       = .TRUE.
        ComOp%calc_grid_HADA      = .FALSE.
        ComOp%save_file_HADA      = .FALSE.
      END IF

      nullify(ComOp%Rvp_spec)
      nullify(ComOp%Cvp_spec)

      IF (associated(ComOp%liste_spec))  THEN
        CALL dealloc_array(ComOp%liste_spec,                            &
                          "ComOp%liste_spec",name_sub)
      END IF
      ComOp%nb_vp_spec  = huge(1)


      IF (allocated(ComOp%sqRhoOVERJac)) THEN
        CALL dealloc_NParray(ComOp%sqRhoOVERJac,"ComOp%sqRhoOVERJac",name_sub)
      END IF
      IF (allocated(ComOp%Jac)) THEN
        CALL dealloc_NParray(ComOp%Jac,"ComOp%Jac",name_sub)
      END IF
      CALL dealloc_SmolyakRep(ComOp%SRep_sqRhoOVERJac)
      CALL dealloc_SmolyakRep(ComOp%SRep_Jac)

      END SUBROUTINE dealloc_ComOp

!================================================================
! ++    write the type param_ComOp
!================================================================
!
      SUBROUTINE write_param_ComOp(ComOp)
      USE mod_system
      IMPLICIT NONE

      TYPE (param_ComOp)   :: ComOp

!----- for debuging ----------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

      write(out_unitp,*) 'WRITE param_ComOp'
      write(out_unitp,*)

      write(out_unitp,*) 'nb_act1',ComOp%nb_act1
      write(out_unitp,*) 'nb_ba,nb_bi,nb_be',ComOp%nb_ba,ComOp%nb_bi,ComOp%nb_be
      write(out_unitp,*) 'nb_bie',ComOp%nb_bie

      write(out_unitp,*) 'MaxNormIndexBasisRep',ComOp%MaxNormIndexBasisRep
      write(out_unitp,*) 'NormIndexBasisRep',associated(ComOp%NormIndexBasisRep)
      IF ( associated(ComOp%NormIndexBasisRep))                              &
                                 write(out_unitp,'(10i6)') ComOp%NormIndexBasisRep(:)

      write(out_unitp,*)
      write(out_unitp,*) 'ADA',ComOp%ADA
      write(out_unitp,*) 'contrac_ba_ON_HAC',ComOp%contrac_ba_ON_HAC
      write(out_unitp,*) 'max_nb_ba_ON_HAC',ComOp%max_nb_ba_ON_HAC
      write(out_unitp,*) 'max_ene_ON_HAC',ComOp%max_ene_ON_HAC
      write(out_unitp,*) 'asso d0Cba_ON_HAC',associated(ComOp%d0Cba_ON_HAC)
      IF (associated(ComOp%d0Cba_ON_HAC))                                   &
               write(out_unitp,*) shape(ComOp%d0Cba_ON_HAC)
      write(out_unitp,*) 'asso Eneba_ON_HAC',associated(ComOp%Eneba_ON_HAC)
      IF (associated(ComOp%Eneba_ON_HAC))                                   &
               write(out_unitp,*) shape(ComOp%Eneba_ON_HAC)
      write(out_unitp,*) 'asso nb_ba_ON_HAC',associated(ComOp%nb_ba_ON_HAC)
      IF (associated(ComOp%nb_ba_ON_HAC))                                &
               write(out_unitp,*) shape(ComOp%nb_ba_ON_HAC)


      write(out_unitp,*)
      write(out_unitp,*) 'HADA file name :',ComOp%file_HADA%name
      write(out_unitp,*)


      write(out_unitp,*) 'nb_vp_spec',ComOp%nb_vp_spec
      IF (associated(ComOp%Rvp_spec))                                   &
               write(out_unitp,*) shape(ComOp%Rvp_spec)
      IF (associated(ComOp%Cvp_spec))                                   &
               write(out_unitp,*) shape(ComOp%Cvp_spec)
      IF (associated(ComOp%liste_spec))                                 &
               write(out_unitp,*) shape(ComOp%liste_spec)


      IF (allocated(ComOp%sqRhoOVERJac)) THEN
        write(out_unitp,*) 'shape sqRhoOVERJac',shape(ComOp%sqRhoOVERJac)
      END IF
      IF (allocated(ComOp%Jac)) THEN
        write(out_unitp,*) 'shape Jac ',shape(ComOp%Jac)
      END IF

      write(out_unitp,*) 'END WRITE param_ComOp'


      END SUBROUTINE write_param_ComOp

      !================================================
      ! Ene     : a table of energy (not necessary sorted)
      ! ZPE     : the ZPE value
      ! Ene_min : if this value is given, it means that all physical energies are larger than Ene_min
      !            => Therefore, ZPE is larger than Ene_min (to avoid holes)
      ! forced  : To force to set-up the ZPE even if ComOp%Set_ZPE=.TRUE.
      !================================================
      SUBROUTINE Set_ZPE_OF_ComOp(ComOp,Ene,ZPE,Ene_min,forced)
      USE mod_system
      IMPLICIT NONE

      TYPE (param_ComOp), intent(inout) :: ComOp

      real (kind=Rkind), intent(in), optional :: ZPE,Ene_min
      real (kind=Rkind), intent(in), optional :: Ene(:)
      logical, intent(in), optional :: forced


      logical :: forced_loc
      real (kind=Rkind) :: ZPE_loc,Ene_min_loc
      integer :: i


      IF (present(forced)) THEN
        forced_loc = forced
      ELSE
        forced_loc = .FALSE.
      END IF

      IF (present(ZPE)) THEN
        ZPE_loc = ZPE
      ELSE
        ZPE_loc = huge(ONE)
      END IF
      IF (present(Ene_min)) THEN
        Ene_min_loc = Ene_min
      ELSE
        Ene_min_loc = -huge(ONE)
      END IF
      ZPE_loc = max(ZPE_loc,Ene_min_loc)

      IF (.NOT. ComOp%Set_ZPE .OR. forced_loc) THEN


        IF (present(Ene)) THEN
          DO i=1,size(Ene)
            IF (Ene(i) >= Ene_min_loc .AND. Ene(i) < ZPE_loc) ZPE_loc = Ene(i)
          END DO
        END IF

        ComOp%ZPE     = ZPE_loc
        ComOp%Set_ZPE = .TRUE.
      END IF

      END SUBROUTINE Set_ZPE_OF_ComOp


      FUNCTION Get_ZPE(Ene,min_Ene)
      USE mod_system
      IMPLICIT NONE

      real (kind=Rkind) :: Get_ZPE

      real (kind=Rkind), intent(in), optional :: min_Ene

      real (kind=Rkind), intent(in) :: Ene(:)

      real (kind=Rkind) :: ZPE

      IF (present(min_Ene)) THEN
        ZPE = min_Ene
      ELSE
        ZPE = Huge(ONE)
      END IF

      IF (size(Ene) > 0) ZPE = min(ZPE,minval(Ene))

      Get_ZPE = ZPE


      END FUNCTION Get_ZPE

      END MODULE mod_ComOp

