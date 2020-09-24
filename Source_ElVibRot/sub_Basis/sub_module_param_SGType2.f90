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
MODULE mod_param_SGType2
USE mod_system
use mod_nDindex, only: type_nDindex, dealloc_nDindex,Write_nDindex,     &
                       alloc_nparray, init_ndval_of_nDindex,            &
                       add_one_to_nDindex, calc_ndi, calc_nDindex,      &
                       dealloc_nparray
USE mod_MPI_aux
IMPLICIT NONE

  PRIVATE

  TYPE param_SGType2
    integer                      :: L1_SparseGrid = huge(1)
    integer                      :: L2_SparseGrid = huge(1)

    integer                      :: L1_SparseBasis = huge(1)
    integer                      :: L2_SparseBasis = huge(1)

    integer                      :: Num_OF_Lmax   = 0 ! use normal L_SparseGrid

    integer                      :: nb0   = 0 ! to deal with several electronic PES, rotational basis, or channels (HAC)
    integer                      :: nb_SG = 0 ! numer of terms
    TYPE (Type_nDindex)          :: nDind_SmolyakRep  ! multidimensional index smolyak grids

    TYPE (Type_nDindex), allocatable :: nDind_DPG(:)    ! multidimensional DP index (nb_SG)
    TYPE (Type_nDindex), allocatable :: nDind_DPB(:)    ! multidimensional DP index (nb_SG)

    !for SGtype=4
    integer (kind=Ikind), allocatable :: tab_iB_OF_SRep_TO_iB(:)   ! size (SRep)

    integer, allocatable :: tab_Sum_nq_OF_SRep(:)     ! size (SRep)
    integer, allocatable :: tab_nq_OF_SRep(:)         ! size (SRep)

    integer, allocatable :: tab_Sum_nb_OF_SRep(:)     ! size (SRep)
    integer, allocatable :: tab_nb_OF_SRep(:)         ! size (SRep)

    ! To intialize nDval(:) for each thread and when Tab_nDval is not allocated
    integer                       :: nb_threads = 0
    integer                       :: nb_tasks   = 0
    integer, allocatable          :: nDval_init(:,:)   ! nDval_init(ndim,nb_threads) table the individual indexes
    integer, allocatable          :: iG_th(:),fG_th(:) ! iG indexes associated to the OpenMP threads

    TYPE (multi_array4),allocatable :: nDI_index_master(:) !< for MPI in action, scheme3 
    Integer(kind=MPI_INTEGER_KIND),allocatable :: reduce_Vlength_master(:) 
    Integer(kind=MPI_INTEGER_KIND),allocatable :: size_ST(:) 
    Integer(kind=MPI_INTEGER_KIND),allocatable :: size_ST_mc(:,:)
    Integer(kind=MPI_INTEGER_KIND)             :: size_psi
    Integer(kind=MPI_INTEGER_KIND)             :: reduce_Vlength !< reduced size of V 
    Integer                                    :: Max_nDI_ib0
    Integer(kind=Ikind),allocatable            :: nDI_index(:)
    Integer(kind=Ikind),allocatable            :: nDI_index_list(:) 
    Integer                                    :: num_nDI_index
    Integer                                    :: V_allcount
    Integer                                    :: V_allcount2
    Logical                                    :: once_action=.TRUE.

  CONTAINS
    PROCEDURE, PRIVATE, PASS(SGType2_1) :: SGType2_2TOSGType2_1
    GENERIC,   PUBLIC  :: assignment(=) => SGType2_2TOSGType2_1
  END TYPE param_SGType2

  TYPE OldParam
    integer                      :: i_SG = 0

    integer                      :: iq    = 0
    integer                      :: iq_SG = 0

    integer                      :: ib    = 0
    integer                      :: ib_SG = 0

    integer, allocatable         :: tab_l_AT_SG(:) ! associated to i_SG
  END TYPE OldParam

 PUBLIC :: param_SGType2, dealloc_SGType2, Set_nDval_init_FOR_SG4, Write_SGType2
 PUBLIC :: OldParam, Write_OldParam
 PUBLIC :: get_iqSG_iSG_FROM_iq, get_Tabiq_Tabil_FROM_iq, get_Tabiq_Tabil_FROM_iq_old
 PUBLIC :: calc_Weight_OF_SRep,dealloc_OldParam

CONTAINS
SUBROUTINE Write_OldParam(OldPara)

TYPE (OldParam), intent(in) :: OldPara

character (len=*), parameter :: name_sub='Write_OldParam'

  write(out_unitp,*) 'BEGINNING ',name_sub
  write(out_unitp,*) 'iq,ib          ',OldPara%iq,OldPara%ib
  write(out_unitp,*) 'i_SG           ',OldPara%i_SG
  write(out_unitp,*) 'iq_SG,ib_SG    ',OldPara%iq_SG,OldPara%ib_SG

  IF (allocated(OldPara%tab_l_AT_SG)) &
  write(out_unitp,*) 'tab_l_AT_SG(:) ',OldPara%tab_l_AT_SG(:)

  write(out_unitp,*) 'END ',name_sub
  CALL flush_perso(out_unitp)

END SUBROUTINE Write_OldParam
SUBROUTINE dealloc_OldParam(OldPara)

TYPE (OldParam), intent(inout) :: OldPara

character (len=*), parameter :: name_sub='dealloc_OldParam'

    OldPara%i_SG = 0

    OldPara%iq    = 0
    OldPara%iq_SG = 0

    OldPara%ib    = 0
    OldPara%ib_SG = 0

    IF (allocated(OldPara%tab_l_AT_SG)) deallocate(OldPara%tab_l_AT_SG)

END SUBROUTINE dealloc_OldParam

SUBROUTINE Write_SGType2(SGType2)
  TYPE (param_SGType2), intent(in) :: SGType2

character (len=*), parameter :: name_sub='Write_SGType2'

  write(out_unitp,*) 'BEGINNING ',name_sub

  write(out_unitp,*) 'L1_SparseGrid ',SGType2%L1_SparseGrid
  write(out_unitp,*) 'L2_SparseGrid ',SGType2%L2_SparseGrid
  write(out_unitp,*) 'L1_SparseBasis ',SGType2%L1_SparseBasis
  write(out_unitp,*) 'L2_SparseBasis ',SGType2%L2_SparseBasis
  write(out_unitp,*) 'Num_OF_Lmax ',SGType2%Num_OF_Lmax
  write(out_unitp,*)

  write(out_unitp,*) 'nb0 ',SGType2%nb0
  write(out_unitp,*) 'nb_SG ',SGType2%nb_SG
  write(out_unitp,*)

  CALL Write_nDindex(SGType2%nDind_SmolyakRep,'nDind_SmolyakRep')
  write(out_unitp,*) 'alloc nDind_DPG',allocated(SGType2%nDind_DPG)
  write(out_unitp,*) 'alloc nDind_DPB',allocated(SGType2%nDind_DPB)
  write(out_unitp,*)

  write(out_unitp,*) '  FOR SG4:'
  write(out_unitp,*) 'alloc tab_iB_OF_SRep_TO_iB',allocated(SGType2%tab_iB_OF_SRep_TO_iB)
  write(out_unitp,*) 'alloc tab_Sum_nq_OF_SRep',allocated(SGType2%tab_Sum_nq_OF_SRep)
  write(out_unitp,*) 'alloc tab_nq_OF_SRep',allocated(SGType2%tab_nq_OF_SRep)
  write(out_unitp,*) 'alloc tab_Sum_nb_OF_SRep',allocated(SGType2%tab_Sum_nb_OF_SRep)
  write(out_unitp,*) 'alloc tab_nb_OF_SRep',allocated(SGType2%tab_nb_OF_SRep)
  write(out_unitp,*)

  write(out_unitp,*) 'nb_threads ',SGType2%nb_threads
  write(out_unitp,*) 'nb_tasks ',SGType2%nb_tasks
  write(out_unitp,*) 'alloc nDval_init',allocated(SGType2%nDval_init)
  write(out_unitp,*) 'alloc iG_th',allocated(SGType2%iG_th)
  write(out_unitp,*) 'alloc fG_th',allocated(SGType2%fG_th)
  write(out_unitp,*)

  write(out_unitp,*) 'END ',name_sub
  CALL flush_perso(out_unitp)

END SUBROUTINE Write_SGType2

SUBROUTINE dealloc_SGType2(SGType2)

TYPE (param_SGType2), intent(inout) :: SGType2

character (len=*), parameter :: name_sub='dealloc_SGType2'

SGType2%L1_SparseGrid = huge(1)
SGType2%L2_SparseGrid = huge(1)
SGType2%Num_OF_Lmax   = 0 ! use normal L_SparseGrid

SGType2%nb_SG = 0
SGType2%nb0   = 0

CALL dealloc_nDindex(SGType2%nDind_SmolyakRep)

IF (allocated(SGType2%nDind_DPG)) THEN
  CALL dealloc_NParray(SGType2%nDind_DPG,'SGType2%nDind_DPG',name_sub)
END IF

IF (allocated(SGType2%nDind_DPB)) THEN
  CALL dealloc_NParray(SGType2%nDind_DPB,'SGType2%nDind_DPB',name_sub)
END IF

IF (allocated(SGType2%tab_iB_OF_SRep_TO_iB)) THEN
  CALL dealloc_NParray(SGType2%tab_iB_OF_SRep_TO_iB,        &
                      'SGType2%tab_iB_OF_SRep_TO_iB',name_sub)
END IF

IF (allocated(SGType2%tab_Sum_nq_OF_SRep)) THEN
  CALL dealloc_NParray(SGType2%tab_Sum_nq_OF_SRep,        &
                      'SGType2%tab_Sum_nq_OF_SRep',name_sub)
END IF
IF (allocated(SGType2%tab_nq_OF_SRep)) THEN
  CALL dealloc_NParray(SGType2%tab_nq_OF_SRep,             &
                      'SGType2%tab_nq_OF_SRep',name_sub)
END IF

IF (allocated(SGType2%tab_Sum_nb_OF_SRep)) THEN
  CALL dealloc_NParray(SGType2%tab_Sum_nb_OF_SRep,        &
                      'SGType2%tab_Sum_nb_OF_SRep',name_sub)
END IF
IF (allocated(SGType2%tab_nb_OF_SRep)) THEN
  CALL dealloc_NParray(SGType2%tab_nb_OF_SRep,        &
                      'SGType2%tab_nb_OF_SRep',name_sub)
END IF


SGType2%nb_threads = 0

IF (allocated(SGType2%nDval_init)) THEN
  CALL dealloc_NParray(SGType2%nDval_init,'SGType2%nDval_init',name_sub)
END IF

IF (allocated(SGType2%iG_th)) THEN
  CALL dealloc_NParray(SGType2%iG_th,'SGType2%iG_th',name_sub)
END IF

IF (allocated(SGType2%fG_th)) THEN
  CALL dealloc_NParray(SGType2%fG_th,'SGType2%fG_th',name_sub)
END IF

END SUBROUTINE dealloc_SGType2

SUBROUTINE SGType2_2TOSGType2_1(SGType2_1,SGType2_2)

CLASS (param_SGType2), intent(inout) :: SGType2_1
TYPE (param_SGType2),  intent(in)    :: SGType2_2

integer :: i

character (len=*), parameter :: name_sub='SGType2_2TOSGType2_1'


SGType2_1%L1_SparseGrid = SGType2_2%L1_SparseGrid
SGType2_1%L2_SparseGrid = SGType2_2%L2_SparseGrid

SGType2_1%L1_SparseBasis = SGType2_2%L1_SparseBasis
SGType2_1%L2_SparseBasis = SGType2_2%L2_SparseBasis

SGType2_1%Num_OF_Lmax   = SGType2_2%Num_OF_Lmax

SGType2_1%nDind_SmolyakRep = SGType2_2%nDind_SmolyakRep
SGType2_1%nb_SG            = SGType2_2%nb_SG
SGType2_1%nb0              = SGType2_2%nb0

IF (allocated(SGType2_2%nDind_DPG)) THEN
  CALL alloc_NParray(SGType2_1%nDind_DPG,(/ SGType2_1%nb_SG /),            &
                    'SGType2_1%nDind_DPG',name_sub)
  DO i=1,SGType2_1%nb_SG
    SGType2_1%nDind_DPG(i) = SGType2_2%nDind_DPG(i)
  END DO
END IF

IF (allocated(SGType2_2%nDind_DPB)) THEN
  CALL alloc_NParray(SGType2_1%nDind_DPB,(/ SGType2_1%nb_SG /),            &
                    'SGType2_1%nDind_DPB',name_sub)
  DO i=1,SGType2_1%nb_SG
    SGType2_1%nDind_DPB(i) = SGType2_2%nDind_DPB(i)
  END DO
END IF

IF (allocated(SGType2_2%tab_iB_OF_SRep_TO_iB)) THEN
  CALL alloc_NParray(SGType2_1%tab_iB_OF_SRep_TO_iB,                    &
                          shape(SGType2_2%tab_iB_OF_SRep_TO_iB),        &
                    'SGType2_1%tab_iB_OF_SRep_TO_iB',name_sub)
  SGType2_1%tab_iB_OF_SRep_TO_iB(:) = SGType2_2%tab_iB_OF_SRep_TO_iB
END IF

IF (allocated(SGType2_2%tab_Sum_nq_OF_SRep)) THEN
  CALL alloc_NParray(SGType2_1%tab_Sum_nq_OF_SRep,                      &
                          shape(SGType2_2%tab_Sum_nq_OF_SRep),          &
                    'SGType2_1%tab_Sum_nq_OF_SRep',name_sub)
  SGType2_1%tab_Sum_nq_OF_SRep(:) = SGType2_2%tab_Sum_nq_OF_SRep
END IF
IF (allocated(SGType2_2%tab_nq_OF_SRep)) THEN
  CALL alloc_NParray(SGType2_1%tab_nq_OF_SRep,                      &
                          shape(SGType2_2%tab_nq_OF_SRep),          &
                    'SGType2_1%tab_nq_OF_SRep',name_sub)
  SGType2_1%tab_nq_OF_SRep(:) = SGType2_2%tab_nq_OF_SRep
END IF


IF (allocated(SGType2_2%tab_Sum_nb_OF_SRep)) THEN
  CALL alloc_NParray(SGType2_1%tab_Sum_nb_OF_SRep,                      &
                          shape(SGType2_2%tab_Sum_nb_OF_SRep),          &
                    'SGType2_1%tab_Sum_nb_OF_SRep',name_sub)
  SGType2_1%tab_Sum_nb_OF_SRep(:) = SGType2_2%tab_Sum_nb_OF_SRep
END IF
IF (allocated(SGType2_2%tab_nb_OF_SRep)) THEN
  CALL alloc_NParray(SGType2_1%tab_nb_OF_SRep,                      &
                          shape(SGType2_2%tab_nb_OF_SRep),          &
                    'SGType2_1%tab_nb_OF_SRep',name_sub)
  SGType2_1%tab_nb_OF_SRep(:) = SGType2_2%tab_nb_OF_SRep
END IF


SGType2_1%nb_threads = SGType2_2%nb_threads
SGType2_1%nb_tasks   = SGType2_2%nb_tasks

IF (allocated(SGType2_2%nDval_init)) THEN
  CALL alloc_NParray(SGType2_1%nDval_init,shape(SGType2_2%nDval_init),  &
                    'SGType2_1%nDval_init',name_sub)
  SGType2_1%nDval_init(:,:) = SGType2_2%nDval_init
END IF

IF (allocated(SGType2_2%iG_th)) THEN
  CALL alloc_NParray(SGType2_1%iG_th,shape(SGType2_2%iG_th),  &
                    'SGType2_1%iG_th',name_sub)
  SGType2_1%iG_th(:) = SGType2_2%iG_th
END IF
IF (allocated(SGType2_2%fG_th)) THEN
  CALL alloc_NParray(SGType2_1%fG_th,shape(SGType2_2%fG_th),  &
                    'SGType2_1%fG_th',name_sub)
  SGType2_1%fG_th(:) = SGType2_2%fG_th
END IF

END SUBROUTINE SGType2_2TOSGType2_1

      SUBROUTINE Set_nDval_init_FOR_SG4(SGType2,version)
      USE mod_system
      !$ USE omp_lib, only : omp_get_max_threads
      IMPLICIT NONE

!----- for the basis set ----------------------------------------------
      TYPE (param_SGType2), intent(inout) :: SGType2
      integer,              intent(in)    :: version


      integer             :: ith,err_sub,nb_threads,ndim


      character (len=:), allocatable :: fformat

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Set_nDval_init_FOR_SG4'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'ndim (nb_basis)',SGType2%nDind_SmolyakRep%ndim
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

  !this line is fine for both openMP and MPI
  IF (SG4_omp == 0) THEN
    nb_threads = 1
  ELSE
    nb_threads = SG4_maxth
  END IF

      ndim        = SGType2%nDind_SmolyakRep%ndim

      ! version 1 currently
      SELECT CASE (version)

      CASE (0)
        CALL Set_nDval_init_FOR_SG4_v0(SGType2,nb_threads,err_sub)

      CASE (1)

        DO
          write(out_unitp,*) ' nb_threads',nb_threads
          ! nDval_init setup
          CALL Set_nDval_init_FOR_SG4_v1(SGType2,nb_threads,err_sub)

          IF (err_sub /= 0) THEN
            nb_threads = nb_threads/2
            IF (nb_threads < 1) EXIT
          ELSE
            EXIT
          END IF

        END DO

      CASE (2)
        STOP 'do not use this version: Set_nDval_init_FOR_SG4_v2'
        CALL Set_nDval_init_FOR_SG4_v2(SGType2,nb_threads,err_sub)

      CASE default
        CALL Set_nDval_init_FOR_SG4_v0(SGType2,nb_threads,err_sub)

      END SELECT

      IF (err_sub /= 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        STOP
      END IF

!-----------------------------------------------------------
      IF (debug) THEN
        fformat = '(i0,a,i0,x,i0,a,' // int_TO_char(ndim) // '(1x,i0))'

        DO ith=1,SGType2%nb_tasks
          write(out_unitp,fformat) ith-1,' iG_th,fG_th ',               &
                                 SGType2%iG_th(ith),SGType2%fG_th(ith), &
                    ' nDval_init: ',SGType2%nDval_init(:,ith)
        END DO
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

      END SUBROUTINE Set_nDval_init_FOR_SG4
      SUBROUTINE Set_nDval_init_FOR_SG4_v0(SGType2,nb_threads,err_sub)
      USE mod_system
      IMPLICIT NONE

!----- for the basis set ----------------------------------------------
      TYPE (param_SGType2), intent(inout) :: SGType2
      integer,              intent(inout) :: nb_threads,err_sub


      integer             :: ith,nqq,nqq_Th,ndim,i_SG,iiG,diG
      integer             :: tab_l(SGType2%nDind_SmolyakRep%ndim)
      integer             :: tab_l0(SGType2%nDind_SmolyakRep%ndim)

      character (len=:), allocatable :: fformat

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Set_nDval_init_FOR_SG4_v0'
      !logical,parameter :: debug=.FALSE.
      logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'ndim (nb_basis)',SGType2%nDind_SmolyakRep%ndim
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

      ndim        = SGType2%nDind_SmolyakRep%ndim

      SGType2%nb_threads = nb_threads
      SGType2%nb_tasks   = nb_threads

      IF (SGType2%nb_tasks > SGType2%nDind_SmolyakRep%max_nDI) THEN
        SGType2%nb_tasks = SGType2%nDind_SmolyakRep%max_nDI
      END IF




      CALL alloc_NParray(SGType2%nDval_init,(/ ndim,SGType2%nb_tasks /),&
                        'SGType2%nDval_init',name_sub)
      SGType2%nDval_init(:,:) = 0

      CALL alloc_NParray(SGType2%iG_th,(/ SGType2%nb_tasks /),  &
                        'SGType2%iG_th',name_sub)
      SGType2%iG_th(:) = 0

      CALL alloc_NParray(SGType2%fG_th,(/ SGType2%nb_tasks /),  &
                        'SGType2%fG_th',name_sub)
      SGType2%fG_th(:) = 0


        DO ith=0,SGType2%nb_tasks-1

          SGType2%iG_th(ith+1) =                          &
              (ith*SGType2%nDind_SmolyakRep%max_nDI)/SGType2%nb_tasks+1
          SGType2%fG_th(ith+1) =                          &
            ((ith+1)*SGType2%nDind_SmolyakRep%max_nDI)/SGType2%nb_tasks


          i_SG = SGType2%iG_th(ith+1)
          CALL init_nDval_OF_nDindex(SGType2%nDind_SmolyakRep,tab_l)

          DO iiG=1,i_SG-1
            CALL ADD_ONE_TO_nDindex(SGType2%nDind_SmolyakRep,tab_l)
          END DO
          SGType2%nDval_init(:,ith+1) = tab_l(:)

        END DO

        err_sub = 0
        IF (SGType2%fG_th(SGType2%nb_tasks) /= SGType2%nDind_SmolyakRep%max_nDI) THEN
          err_sub = 1
        END IF


!-----------------------------------------------------------
      IF (debug) THEN
        fformat = '(i0,a,i0,x,i0,a,' // int_TO_char(ndim) // '(1x,i0))'

        DO ith=1,SGType2%nb_tasks
          write(out_unitp,fformat) ith-1,' iG_th,fG_th ',               &
                                 SGType2%iG_th(ith),SGType2%fG_th(ith), &
                    ' nDval_init: ',SGType2%nDval_init(:,ith)
        END DO
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------
      END SUBROUTINE Set_nDval_init_FOR_SG4_v0
      SUBROUTINE Set_nDval_init_FOR_SG4_v1(SGType2,nb_threads,err_sub)
      USE mod_system
      IMPLICIT NONE

!----- for the basis set ----------------------------------------------
      TYPE (param_SGType2), intent(inout) :: SGType2
      integer,              intent(inout) :: nb_threads,err_sub


      integer             :: ith,nqq,nqq_Th,ndim,i_SG,iiG
      integer             :: tab_l(SGType2%nDind_SmolyakRep%ndim)
      integer             :: tab_l0(SGType2%nDind_SmolyakRep%ndim)

      character (len=:), allocatable :: fformat

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Set_nDval_init_FOR_SG4_v1'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'ndim (nb_basis)',SGType2%nDind_SmolyakRep%ndim
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

      ndim        = SGType2%nDind_SmolyakRep%ndim

      SGType2%nb_threads = nb_threads
      SGType2%nb_tasks   = nb_threads


      CALL alloc_NParray(SGType2%nDval_init,(/ ndim,SGType2%nb_tasks /),&
                        'SGType2%nDval_init',name_sub)
      SGType2%nDval_init(:,:) = 0

      CALL alloc_NParray(SGType2%iG_th,(/ SGType2%nb_tasks /),  &
                        'SGType2%iG_th',name_sub)
      SGType2%iG_th(:) = 0

      CALL alloc_NParray(SGType2%fG_th,(/ SGType2%nb_tasks /),  &
                        'SGType2%fG_th',name_sub)
      SGType2%fG_th(:) = 0


        nqq_Th      = sum(SGType2%tab_nq_OF_SRep(:)) / SGType2%nb_tasks
        nqq         = 0

        ith = 1
        CALL init_nDval_OF_nDindex(SGType2%nDind_SmolyakRep,tab_l)

        SGType2%nDval_init(:,ith) = tab_l(:)
        SGType2%iG_th(ith)        = 1

        DO i_SG=1,SGType2%nb_SG
          CALL ADD_ONE_TO_nDindex(SGType2%nDind_SmolyakRep,tab_l,iG=i_SG)

          nqq         = nqq         + SGType2%tab_nq_OF_SRep(i_SG)

          IF (nqq > nqq_Th) THEN
            ! end for the current thread
            SGType2%fG_th(ith) = i_SG-1

            ! for the new thread
            ith = ith+1
            SGType2%iG_th(ith)        = i_SG
            SGType2%nDval_init(:,ith) = tab_l0(:) ! tab_l before ADD_ONE

            nqq = 0
          END IF
          tab_l0(:) = tab_l(:)

        END DO
        ! end of the last thread
        SGType2%fG_th(ith) = SGType2%nb_SG

        err_sub = 0
        IF (count(SGType2%iG_th == 0) > 0) err_sub = 1
        IF (count(SGType2%fG_th == 0) > 0) err_sub = 1

      IF (debug) THEN
        fformat = '(i0,a,i0,x,i0,a,' // int_TO_char(ndim) // '(1x,i0))'

        DO ith=1,SGType2%nb_tasks
          write(out_unitp,fformat) ith-1,' iG_th,fG_th ',               &
                                 SGType2%iG_th(ith),SGType2%fG_th(ith), &
                    ' nDval_init: ',SGType2%nDval_init(:,ith)
        END DO
        CALL flush_perso(out_unitp)
      END IF

      !err /= 0 means nb_threads is too large => table are deallocated
      IF (err_sub /= 0) THEN
        CALL dealloc_NParray(SGType2%nDval_init,'SGType2%nDval_init',name_sub)
        CALL dealloc_NParray(SGType2%iG_th,'SGType2%iG_th',name_sub)
        CALL dealloc_NParray(SGType2%fG_th,'SGType2%fG_th',name_sub)
      END IF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

      END SUBROUTINE Set_nDval_init_FOR_SG4_v1

      SUBROUTINE Set_nDval_init_FOR_SG4_v2(SGType2,nb_threads,err_sub)
      USE mod_system
      IMPLICIT NONE

!----- for the basis set ----------------------------------------------
      TYPE (param_SGType2), intent(inout) :: SGType2
      integer,              intent(inout) :: nb_threads,err_sub


      integer             :: ith,nqq,nqq_Th,ndim,i_SG,iiG,diG
      integer             :: tab_l(SGType2%nDind_SmolyakRep%ndim)
      integer             :: tab_l0(SGType2%nDind_SmolyakRep%ndim)

      character (len=:), allocatable :: fformat

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Set_nDval_init_FOR_SG4_v2'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'ndim (nb_basis): ',SGType2%nDind_SmolyakRep%ndim
        write(out_unitp,*) 'nb_threads:      ',nb_threads

        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

      ndim        = SGType2%nDind_SmolyakRep%ndim

      SGType2%nb_threads = nb_threads
      SGType2%nb_tasks   = nb_threads*10

      IF (SGType2%nb_threads > SGType2%nDind_SmolyakRep%max_nDI) THEN
        SGType2%nb_tasks   = SGType2%nDind_SmolyakRep%max_nDI
        SGType2%nb_threads = SGType2%nDind_SmolyakRep%max_nDI
      END IF

      IF (SGType2%nb_tasks > SGType2%nDind_SmolyakRep%max_nDI) THEN
        SGType2%nb_tasks   = SGType2%nDind_SmolyakRep%max_nDI
      END IF

      IF (debug) THEN
        write(out_unitp,*) 'SGType2%nb_tasks',SGType2%nb_tasks
        CALL flush_perso(out_unitp)
      END IF


      CALL alloc_NParray(SGType2%nDval_init,(/ ndim,SGType2%nb_tasks /),&
                        'SGType2%nDval_init',name_sub)
      SGType2%nDval_init(:,:) = 0

      CALL alloc_NParray(SGType2%iG_th,(/ SGType2%nb_tasks /),  &
                        'SGType2%iG_th',name_sub)
      SGType2%iG_th(:) = 0

      CALL alloc_NParray(SGType2%fG_th,(/ SGType2%nb_tasks /),  &
                        'SGType2%fG_th',name_sub)
      SGType2%fG_th(:) = 0


        DO ith=0,SGType2%nb_tasks-1

          SGType2%iG_th(ith+1) =                          &
              (ith*SGType2%nDind_SmolyakRep%max_nDI)/SGType2%nb_tasks+1
          SGType2%fG_th(ith+1) =                          &
            ((ith+1)*SGType2%nDind_SmolyakRep%max_nDI)/SGType2%nb_tasks


          i_SG = SGType2%iG_th(ith+1)
          CALL init_nDval_OF_nDindex(SGType2%nDind_SmolyakRep,tab_l)

          DO iiG=1,i_SG-1
            CALL ADD_ONE_TO_nDindex(SGType2%nDind_SmolyakRep,tab_l)
          END DO
          SGType2%nDval_init(:,ith+1) = tab_l(:)

        END DO

        err_sub = 0
        IF (SGType2%fG_th(SGType2%nb_tasks) /= SGType2%nDind_SmolyakRep%max_nDI) THEN
          err_sub = 1
        END IF


!-----------------------------------------------------------
      IF (debug) THEN
        fformat = '(i0,a,i0,x,i0,a,' // int_TO_char(ndim) // '(1x,i0))'

        DO ith=1,SGType2%nb_tasks
          write(out_unitp,fformat) ith-1,' iG_th,fG_th ',               &
                                 SGType2%iG_th(ith),SGType2%fG_th(ith), &
                              ' nDval_init: ',SGType2%nDval_init(:,ith)
        END DO
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------
      END SUBROUTINE Set_nDval_init_FOR_SG4_v2

      RECURSIVE SUBROUTINE calc_Weight_OF_SRep(WeightSG,nDind_SmolyakRep)
      USE mod_system
      IMPLICIT NONE

      TYPE (Type_nDindex),             intent(in)    :: nDind_SmolyakRep
      real (kind=Rkind),               intent(inout) :: WeightSG(nDind_SmolyakRep%Max_nDI)

!---------------------------------------------------------------------
      real (kind=Rkind) :: binomial ! function
!---------------------------------------------------------------------

      integer             :: i,i_SG,i_SGm,DeltaL,max_print
      integer             :: tab_l(nDind_SmolyakRep%ndim)
      integer             :: tab_lm(nDind_SmolyakRep%ndim)

      logical             :: old = .TRUE.
      real (kind=Rkind)   :: WeightSG_tmp(nDind_SmolyakRep%Max_nDI)


!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='calc_Weight_OF_SRep'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
      END IF
!-----------------------------------------------------------

    IF (count(nDind_SmolyakRep%nDNum_OF_Lmax == 0) == nDind_SmolyakRep%ndim .AND. &
        nDind_SmolyakRep%MaxCoupling >= nDind_SmolyakRep%ndim) THEN
      ! it works only when L1max or L2max are not used and
      !     when the max number of coupling terms is >= than ndim
      CALL init_nDval_OF_nDindex(nDind_SmolyakRep,tab_l)
      DO i_SG=1,nDind_SmolyakRep%Max_nDI
        CALL ADD_ONE_TO_nDindex(nDind_SmolyakRep,tab_l,iG=i_SG)
        DeltaL = nDind_SmolyakRep%Lmax - sum(tab_l)
        !IF (DeltaL < 0) STOP 'DeltaL < 0'
        !IF (DeltaL > nDind_SmolyakRep%ndim -1) STOP 'DeltaL > ndim-1'
        IF (DeltaL < 0 .OR. DeltaL > nDind_SmolyakRep%ndim -1) THEN
          WeightSG(i_SG) = ZERO
        ELSE
          IF (mod(DeltaL,2) == 0) THEN
            WeightSG(i_SG) =  binomial(nDind_SmolyakRep%ndim-1,deltaL)
          ELSE
            WeightSG(i_SG) = -binomial(nDind_SmolyakRep%ndim-1,deltaL)
          END IF
        END IF

        IF (debug) write(out_unitp,*) 'i_SG,nDval,coef',i_SG,tab_l(:),WeightSG(i_SG)
      END DO
    ELSE ! here the Smolyak rep in Delta_S is transformed in S to get the correct WeightSG
      WeightSG(:) = ONE
      DO i=1,nDind_SmolyakRep%ndim
        WeightSG_tmp(:) = ZERO

        CALL init_nDval_OF_nDindex(nDind_SmolyakRep,tab_l)
        DO i_SG=1,nDind_SmolyakRep%Max_nDI
          CALL ADD_ONE_TO_nDindex(nDind_SmolyakRep,tab_l,iG=i_SG)
          ! DeltaS_(li) = S_(li) - S_(li-1)

          ! S_(li) contribution
          WeightSG_tmp(i_SG) = WeightSG_tmp(i_SG) + WeightSG(i_SG)

          ! -S_(li-1) contribution
          IF (tab_l(i) > 0) THEN
            tab_lm(:) = tab_l(:)
            tab_lm(i) = tab_l(i) -1
            CALL calc_nDI(i_SGm,tab_lm,nDind_SmolyakRep)
            WeightSG_tmp(i_SGm) = WeightSG_tmp(i_SGm) - WeightSG(i_SG)
          END IF

        END DO
        WeightSG(:) = WeightSG_tmp(:)
      END DO
      !STOP 'not yet'
    END IF

    IF (debug) write(out_unitp,*) 'count zero weight: ',count(abs(WeightSG) <= ONETENTH**6)
!-----------------------------------------------------------
    IF (debug .OR. print_level > 1) THEN
      max_print = nDind_SmolyakRep%Max_nDI
      IF (.NOT. debug) max_print = min(100,max_print)

      CALL init_nDval_OF_nDindex(nDind_SmolyakRep,tab_l)
      DO i_SG=1,max_print
        CALL ADD_ONE_TO_nDindex(nDind_SmolyakRep,tab_l,iG=i_SG)
        write(out_unitp,*) 'i_SG,nDval,coef',i_SG,tab_l(:),WeightSG(i_SG)
      END DO
      IF (max_print < nDind_SmolyakRep%Max_nDI) THEN
         write(out_unitp,*) 'i_SG,nDval,coef ....'
      END IF
    END IF

    IF (debug) THEN
      write(out_unitp,*) 'END ',name_sub
    END IF
!-----------------------------------------------------------
  END SUBROUTINE calc_Weight_OF_SRep


! from an index iq (global index of the multidimentional Smolyak grid) get:
!  - iSG:  the numero of the direct-product grid of the Smolyak grid
!  - iqSG: the numero of the grid point from the "iSG" direct-product grid
SUBROUTINE get_iqSG_iSG_FROM_iq(iSG,iqSG,iq,SGType2,OldPara,err_sub)

integer, intent(inout) :: iSG,iqSG
integer, intent(inout) :: err_sub

 TYPE (OldParam), intent(inout), optional :: OldPara

integer, intent(in) :: iq
TYPE (param_SGType2), intent(in) :: SGType2

integer :: nq,iSG_loc
!----- for debuging --------------------------------------------------
logical, parameter :: debug=.FALSE.
!logical, parameter :: debug=.TRUE.
character (len=*), parameter :: name_sub = 'get_iqSG_iSG_FROM_iq'
!-----------------------------------------------------------
err_sub = 0
IF (debug) THEN
  write(out_unitp,*) 'BEGINNING ',name_sub
  CALL flush_perso(out_unitp)
END IF

IF (present(OldPara)) THEN
  !write(out_unitp,*) 'OldPara ',name_sub,OldPara
  iSG = OldPara%i_SG
END IF


!write(out_unitp,*) 'alloc tab_Sum_nq_OF_SRep',allocated(SGType2%tab_Sum_nq_OF_SRep)
!flush(out_unitp)
IF (iSG > 1 .AND. iSG <= size(SGType2%tab_Sum_nq_OF_SRep)) THEN
  iqSG = iq - SGType2%tab_Sum_nq_OF_SRep(iSG-1)

  IF (iqSG < 1) THEN
    iqSG = iq
    iSG  = 1
  END IF
ELSE
  iqSG  = iq
  iSG   = 1
END IF


DO iSG_loc=iSG,SGType2%nb_SG
  nq = SGType2%tab_nq_OF_SRep(iSG_loc)
  IF (iqSG <= nq) EXIT
  iqSG = iqSG - nq
END DO

iSG       = iSG_loc
IF (present(OldPara)) THEN
  OldPara%i_SG = iSG
  OldPara%iq   = iq
END IF

IF (iqSG < 1) STOP 'iqSG < 1'


IF (debug) THEN
  write(out_unitp,*) 'iq,iSG,iqSG',iq,iSG,iqSG
  write(out_unitp,*) 'END ',name_sub
END IF

END SUBROUTINE get_iqSG_iSG_FROM_iq

! from an index iq (global index of the multidimentional Smolyak grid) get:
!  - i_SG:  the numero of the direct-product grid of the Smolyak grid
!  - iq_SG: the numero of the grid point from the "i_SG" direct-product grid
!  - Tabil(:): the "l" indexes corresponding to "i_SG"
!  - Tabiq(:): the "iq" indexes corresponding to "iq_SG"
SUBROUTINE get_Tabiq_Tabil_FROM_iq(Tabiq,Tabil,i_SG,iq_SG,iq,SGType2,OldPara,err_sub)
!$    USE omp_lib, only : omp_get_thread_num
IMPLICIT NONE

integer, intent(inout) :: Tabiq(:)
integer, intent(inout) :: Tabil(:)
integer, intent(inout) :: i_SG,iq_SG
integer, intent(inout) :: err_sub

 TYPE (OldParam), intent(inout), optional :: OldPara

integer, intent(in) :: iq
TYPE (param_SGType2), intent(in) :: SGType2

integer :: nq,i_SG_loc
logical :: old_l
!----- for debuging --------------------------------------------------
logical, parameter :: debug=.FALSE.
!logical, parameter :: debug=.TRUE.
character (len=*), parameter :: name_sub = 'get_Tabiq_Tabil_FROM_iq'
!-----------------------------------------------------------

err_sub = 0
IF (debug) THEN
  write(out_unitp,*) 'BEGINNING ',name_sub
  !$ write(out_unitp,*) 'thread_num:',omp_get_thread_num()
  write(out_unitp,*) ' iq',iq
  write(out_unitp,*) ' i_SG,iq_SG',i_SG,iq_SG
  write(out_unitp,*) ' Tabiq',Tabiq
  write(out_unitp,*) ' Tabil',Tabil
  IF (present(OldPara)) CALL Write_OldParam(OldPara)
  CALL flush_perso(out_unitp)
END IF

!first calculation of i_SG and iq_SG from iq and OldPara (if available)
IF (present(OldPara)) THEN
  !write(out_unitp,*) 'OldPara ',name_sub,OldPara
  i_SG = OldPara%i_SG
ELSE
  i_SG = 0
END IF

IF (i_SG > 1 .AND. i_SG <= size(SGType2%tab_Sum_nq_OF_SRep)) THEN
  iq_SG = iq - SGType2%tab_Sum_nq_OF_SRep(i_SG-1)

  IF (iq_SG < 1) THEN
    iq_SG = iq
    i_SG  = 1
  END IF
ELSE
  iq_SG  = iq
  i_SG   = 1
END IF


DO i_SG_loc=i_SG,SGType2%nb_SG
  nq = SGType2%tab_nq_OF_SRep(i_SG_loc)
  IF (iq_SG <= nq) EXIT
  iq_SG = iq_SG - nq
END DO

i_SG       = i_SG_loc

IF (i_SG > SGType2%nDind_SmolyakRep%Max_nDI) THEN
  write(out_unitp,*) 'ERROR in ',name_sub
  write(out_unitp,*) ' iq_SG',iq_SG
  write(out_unitp,*) ' i_SG',i_SG
  write(out_unitp,*) ' nDind_SmolyakRep%Max_nDI',SGType2%nDind_SmolyakRep%Max_nDI
  IF (present(OldPara)) CALL Write_OldParam(OldPara)
  CALL flush_perso(out_unitp)
  STOP 'i_SG too large'
END IF
IF (iq_SG < 1) THEN
  write(out_unitp,*) 'ERROR in ',name_sub
  write(out_unitp,*) ' iq_SG',iq_SG
  write(out_unitp,*) ' iq_SG < 1'
  IF (present(OldPara)) CALL Write_OldParam(OldPara)
  CALL flush_perso(out_unitp)
  STOP 'iq_SG < 1'
END IF

  IF (debug) write(out_unitp,*) ' i_SG,iq_SG,iq',i_SG,iq_SG,iq
  CALL flush_perso(out_unitp)


!2d calculation of Tabil from i_SG and OldPara (if available)
IF (present(OldPara)) THEN

  old_l = (i_SG >= OldPara%i_SG) .AND. allocated(OldPara%tab_l_AT_SG)

  IF (old_l) THEN

    Tabil(:) = OldPara%tab_l_AT_SG
    DO i_SG_loc=OldPara%i_SG+1,i_SG
      CALL ADD_ONE_TO_nDindex(SGType2%nDind_SmolyakRep,Tabil,iG=i_SG_loc,err_sub=err_sub)
    END DO

  ELSE
    CALL calc_nDindex(SGType2%nDind_SmolyakRep,i_SG,Tabil,err_sub)
  END IF

ELSE
  CALL calc_nDindex(SGType2%nDind_SmolyakRep,i_SG,Tabil,err_sub)
END IF

IF (err_sub /= 0) THEN
  write(out_unitp,*) ' SGType2%nDind_SmolyakRep'
  write(out_unitp,*) ' ERROR in ',name_sub
  write(out_unitp,*) ' i_SG,iq_SG,iq',i_SG,iq_SG,iq
  write(out_unitp,*) ' Tabil',Tabil
  write(out_unitp,*) '  from SGType2%nDind_SmolyakRep',i_SG
  IF (present(OldPara)) CALL Write_OldParam(OldPara)
  err_sub = 1
  RETURN
  !STOP 'calc_nDindex'
END IF
  IF (debug) write(out_unitp,*) ' Tabil',i_SG,' : ',Tabil
  CALL flush_perso(out_unitp)


! Save the parameters in OldPara if OldPara is present
IF (present(OldPara)) THEN
  OldPara%i_SG        = i_SG
  OldPara%iq_SG       = iq_SG
  OldPara%iq          = iq
  OldPara%tab_l_AT_SG = Tabil
END IF

!3d calculation of Tabiq from i_SG and iq_SG
  Tabiq(:) = 0
  !CALL calc_nDval_m1(Tabiq,SGType2%tab_nq_OF_SRep(i_SG),nDsize,size(Tabil))
  !CALL calc_nDindex(SGType2%nDind_DPG(i_SG),iq_SG,Tabiq,err_sub)
  IF (err_sub /= 0) THEN
    write(out_unitp,*) ' SGType2%nDind_DPG(i_SG)'
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) ' i_SG,iq_SG,iq',i_SG,iq_SG,iq
    write(out_unitp,*) ' Tabiq',Tabiq
    write(out_unitp,*) ' Tabil',Tabil
    write(out_unitp,*) '  from SGType2%nDind_DPG(i_SG)',i_SG
    err_sub = 2
    RETURN
    !STOP 'calc_nDindex'
   END IF

  !IF (debug) write(out_unitp,*) ' Tabiq',Tabiq
  !CALL flush_perso(out_unitp)

IF (debug) THEN
  write(out_unitp,*) 'iq,i_SG,iq_SG',iq,i_SG,iq_SG
  write(out_unitp,*) 'iq,i_SG,Tabil',iq,i_SG,':',Tabil
  write(out_unitp,*) 'iq,i_SG,Tabiq',iq,i_SG,':',Tabiq
  write(out_unitp,*) 'END ',name_sub
  CALL flush_perso(out_unitp)
END IF

END SUBROUTINE get_Tabiq_Tabil_FROM_iq



SUBROUTINE get_Tabiq_Tabil_FROM_iq_old(Tabiq,Tabil,i_SG,iq_SG,iq,SGType2)

integer, intent(inout) :: Tabiq(:)
integer, intent(inout) :: Tabil(:)
integer, intent(inout) :: i_SG,iq_SG

integer, intent(in) :: iq
TYPE (param_SGType2), intent(in) :: SGType2


integer :: nq
integer :: err_sub

!----- for debuging --------------------------------------------------
 logical, parameter :: debug=.FALSE.
!logical, parameter :: debug=.TRUE.
character (len=*), parameter :: name_sub = 'get_Tabiq_Tabil_FROM_iq_old'
!-----------------------------------------------------------
IF (debug) THEN
  write(out_unitp,*) 'BEGINNING ',name_sub
END IF

iq_SG           = iq

DO i_SG=1,SGType2%nb_SG
  nq = SGType2%nDind_DPG(i_SG)%Max_nDI
  IF (iq_SG <= nq) EXIT
  iq_SG = iq_SG - nq
END DO

  CALL calc_nDindex(SGType2%nDind_SmolyakRep,i_SG,Tabil,err_sub)
  IF (err_sub /= 0) THEN
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) '  from SGType2%nDind_SmolyakRep',i_SG
    STOP 'calc_nDindex'
  END IF

  CALL calc_nDindex(SGType2%nDind_DPG(i_SG),iq_SG,Tabiq,err_sub)
  IF (err_sub /= 0) THEN
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) '  from SGType2%nDind_DPG(i_SG)',i_SG
    STOP 'calc_nDindex'
  END IF

IF (debug) THEN
  write(out_unitp,*) 'iq,i_SG,iq_SG',iq,i_SG,iq_SG
  write(out_unitp,*) 'iq,i_SG,Tabil',iq,i_SG,':',Tabil
  write(out_unitp,*) 'iq,i_SG,Tabiq',iq,i_SG,':',Tabiq
  write(out_unitp,*) 'END ',name_sub
END IF

END SUBROUTINE get_Tabiq_Tabil_FROM_iq_old

END MODULE mod_param_SGType2
