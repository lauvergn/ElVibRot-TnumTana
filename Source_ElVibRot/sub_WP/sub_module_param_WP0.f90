!===========================================================================
!===========================================================================
!This file is part of ElVibRot.
!
!    ElVibRot is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ElVibRot is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ElVibRot.  If not, see <http://www.gnu.org/licenses/>.
!
!    Copyright 2015 David Lauvergnat [1]
!      with contributions of
!        Josep Maria Luis (optimization) [2]
!        Ahai Chen (MPI) [1,4]
!        Lucien Dupuy (CRP) [5]
!
![1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France
![2]: Institut de Química Computacional and Departament de Química,
!        Universitat de Girona, Catalonia, Spain
![3]: Department of Chemistry, Aarhus University, DK-8000 Aarhus C, Denmark
![4]: Maison de la Simulation USR 3441, CEA Saclay, France
![5]: Laboratoire Univers et Particule de Montpellier, UMR 5299,
!         Université de Montpellier, France
!
!    ElVibRot includes:
!        - Tnum-Tana under the GNU LGPL3 license
!        - Somme subroutines of John Burkardt under GNU LGPL license
!             http://people.sc.fsu.edu/~jburkardt/
!        - Somme subroutines of SHTOOLS written by Mark A. Wieczorek under BSD license
!             http://shtools.ipgp.fr
!        - Some subroutine of QMRPack (see cpyrit.doc) Roland W. Freund and Noel M. Nachtigal:
!             https://www.netlib.org/linalg/qmr/
!
!===========================================================================
!===========================================================================
      MODULE mod_param_WP0
      USE mod_system
      USE mod_file
      IMPLICIT NONE

        TYPE param_WP0

        integer         :: nb_WP0             ! default: 1
        logical         :: read_file          ! default: F
        logical         :: read_listWP0       ! default: F
        logical         :: New_Read_WP0       ! default: F. if T, use a new subroutine to read the WP0
        logical         :: WP0restart         ! restart
        real (kind=Rkind)   :: Trestart       ! time of the restart (id T0)
        TYPE (param_file)   :: file_WP0
        logical         :: WP0cplx            ! default = F
        logical         :: lect_WP0GridRep        ! read WP0 on the grid
        logical         :: lect_WP0BasisRep        ! read WP0 on the basis
        logical         :: lect_WP0BasisRepall     ! read WP0 on the basis
                                              ! (for some given basis functions index)
                                              ! Usefull for spectral paropagation or test
        integer           :: WP0n_h ,WP0nb_elec ! indices of the harmonic and electronic channel for WP0
        integer           :: WP0_DIP            ! =0 no dipole moment =1,2,3 => mhux,mhuy,mhuy
        real (kind=Rkind) :: th_WP0 = ZERO      ! mixing angle between the initial WP0 and Dip.WP0
                                                ! so that NewWP0=cos(th_WP0)*Dip.WP0 + sin(th_WP0).WP0

        integer         :: WP0nrho            ! WP0nrho : define how to normalize WP0
        logical         :: WP0BasisRep             ! calculation of WP0BasisRep with the following parameters
        integer         :: WP0_nb_CleanChannel   ! remove some adiabatic or electronic channels
        integer, allocatable:: WP0_CleanChannellist(:) !list of adiabatic or electronic channels to be remove

        integer         ::       nb_act1      ! nb active dimension
!       for each variable Qi : exp[-((Q-Qeq)/sigma)2+i*imp_k*(Q-Qeq)]
        real (kind=Rkind), allocatable :: WP0sigma(:)  ! WP0sigma(nb_act1) : sigma for WP0
        real (kind=Rkind), allocatable :: WP0Qeq(:)    ! WP0Qeq(nb_act1)   : position of WP0
        real (kind=Rkind), allocatable :: WP0imp_k(:)  ! WP0imp_k(nb_act1) : impultion for WP0
        real (kind=Rkind), allocatable :: WP0phase(:)  ! WP0imp_k(nb_act1) : impultion for WP0

        END TYPE param_WP0

        CONTAINS

!================================================================
!
!    alloc / dealloc param_WP0
!
!================================================================
      !!@description: alloc param_WP0
      !!@param: para_poly
      SUBROUTINE alloc_param_WP0(para_WP0,WP0Grid_Gaussian,WP0_CleanChannel)
      IMPLICIT NONE
      TYPE (param_WP0), intent(inout) :: para_WP0
      logical, intent(in) :: WP0Grid_Gaussian,WP0_CleanChannel

        IF (para_WP0%nb_act1 > 0 .AND. WP0Grid_Gaussian) THEN
          IF (.NOT. allocated(para_WP0%WP0sigma)) THEN
            CALL alloc_NParray(para_WP0%WP0sigma,(/para_WP0%nb_act1/),    &
                              "para_WP0%WP0sigma","alloc_param_WP0")
          END IF

          IF (.NOT. allocated(para_WP0%WP0Qeq)) THEN
            CALL alloc_NParray(para_WP0%WP0Qeq,(/para_WP0%nb_act1/),      &
                              "para_WP0%WP0Qeq","alloc_param_WP0")
          END IF

          IF (.NOT. allocated(para_WP0%WP0imp_k)) THEN
            CALL alloc_NParray(para_WP0%WP0imp_k,(/para_WP0%nb_act1/),    &
                              "para_WP0%WP0imp_k","alloc_param_WP0")
          END IF

          IF (.NOT. allocated(para_WP0%WP0phase)) THEN
            CALL alloc_NParray(para_WP0%WP0phase,(/para_WP0%nb_act1/),    &
                              "para_WP0%WP0phase","alloc_param_WP0")
          END IF

        END IF

        IF (WP0_CleanChannel .AND.                                      &
            .NOT. allocated(para_WP0%WP0_CleanChannellist)) THEN
          CALL alloc_NParray(para_WP0%WP0_CleanChannellist,               &
                                     (/para_WP0%WP0_nb_CleanChannel/),  &
                          "para_WP0%WP0_CleanChannellist","alloc_param_WP0")
        END IF

      END SUBROUTINE alloc_param_WP0

      !!@description: dealloc param_WP0
      !!@param: para_poly
      SUBROUTINE dealloc_param_WP0(para_WP0)
      IMPLICIT NONE
      TYPE (param_WP0), intent(inout) :: para_WP0

        IF (allocated(para_WP0%WP0sigma)) THEN
          CALL dealloc_NParray(para_WP0%WP0sigma,                         &
                              "para_WP0%WP0sigma","dealloc_param_WP0")
        END IF

        IF (allocated(para_WP0%WP0Qeq)) THEN
          CALL dealloc_NParray(para_WP0%WP0Qeq,                           &
                              "para_WP0%WP0Qeq","dealloc_param_WP0")
        END IF

        IF (allocated(para_WP0%WP0imp_k)) THEN
          CALL dealloc_NParray(para_WP0%WP0imp_k,                         &
                              "para_WP0%WP0imp_k","dealloc_param_WP0")
        END IF

        IF (allocated(para_WP0%WP0phase)) THEN
          CALL dealloc_NParray(para_WP0%WP0phase,                         &
                              "para_WP0%WP0phase","dealloc_param_WP0")
        END IF

        IF (allocated(para_WP0%WP0_CleanChannellist)) THEN
          CALL dealloc_NParray(para_WP0%WP0_CleanChannellist,             &
                            "para_WP0%WP0_CleanChannellist","dealloc_param_WP0")
        END IF
      END SUBROUTINE dealloc_param_WP0

      END MODULE mod_param_WP0

