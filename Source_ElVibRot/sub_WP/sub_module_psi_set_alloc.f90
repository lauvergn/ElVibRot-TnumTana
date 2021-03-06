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
 MODULE mod_psi_set_alloc
 USE mod_system
 USE mod_basis
 USE mod_type_ana_psi
 IMPLICIT NONE

 PRIVATE

        TYPE param_psi
        !-- for the basis set ----------------------------------------------
        TYPE (param_AllBasis), pointer       :: para_AllBasis => null()   ! .TRUE. pointer
        TYPE (Basis),          pointer       :: BasisnD       => null()   ! .TRUE. pointer
        TYPE (Basis),          pointer       :: Basis2n       => null()   ! .TRUE. pointer
        integer                              :: symab         = -1

        logical       :: init             = .FALSE.   ! (F) all constantes are initialized
        logical       :: cplx             = .FALSE.   ! (T) if T => psi is complex (use of Cvec)
                                                      ! otherwise psi is real (use of Rvec)
        logical       :: BasisRep         = .TRUE.    ! (T) BasisRep
        logical       :: GridRep          = .FALSE.   ! (F) GridRep psi
        logical       :: BasisRep_contrac = .FALSE.    ! (T) BasisRep

        Logical       :: SRG_MPI          = .FALSE.    ! (F) Smolyak rep. on grids with MPI
        Logical       :: SRB_MPI          = .FALSE.    ! (F) Smolyak rep. on Basis with MPI

        integer       :: nb_ba            =0
        integer       :: nb_bi            =0
        integer       :: nb_be            =0
        integer       :: nb_bRot          =0
        integer       :: nb_baie          =0
        integer       :: nb_tot           =0
        integer       :: nb_tot_contrac   =0
        integer       :: nb_tot_uncontrac =0

        integer       :: nb_qa            =0
        integer       :: nb_qaie          =0
        integer       :: nb_paie          =0

        integer       :: nb_act1          =0
        integer       :: nb_act           =0

        integer       :: nb_basis_act1    =0
        integer       :: nb_basis         =0
        integer       :: max_dim          =0

        real (kind=Rkind),    allocatable :: RvecB(:) ! RvecB(nb_tot)
        complex (kind=Rkind), allocatable :: CvecB(:) ! CvecB(nb_tot)

        real (kind=Rkind),    allocatable :: TDParam(:) ! Non linear parameters (basis-set)
        integer                           :: nb_TDParam =0    ! number of extra Parameters


        real (kind=Rkind),    allocatable :: RvecG(:) ! RvecG(nb_qaie)
        complex (kind=Rkind), allocatable :: CvecG(:) ! CvecG(nb_qaie)

        Real(kind=Rkind),allocatable      :: SR_B(:,:) ! Smolyak rep. on basis
        Real(kind=Rkind),allocatable      :: SR_G(:,:) ! Smolyak rep. on grid
        Integer,allocatable               :: SR_B_index(:)    !< index for iG in SR_B
        Integer,allocatable               :: SR_B_length(:)   !< length of SR_B
        Integer,allocatable               :: SR_G_index(:)    !< index for iG in SR_G
        Integer,allocatable               :: SR_G_length(:)   !< length of SR_G

        complex (kind=Rkind) :: CAvOp    = (ZERO,ZERO) ! average value for an operator (usualy H)
        integer              :: IndAvOp  = -1          ! operator type  (usualy H, IndAvOp=0)
        logical              :: convAvOp = .FALSE.     ! TRUE, if the average value is calculated with a converged wavefunction

        real (kind=Rkind)    :: norm2    = -ONE       ! norm^2 of psi

        CONTAINS
          PROCEDURE, PRIVATE, PASS(psi1)  :: psi1_pluseq_psi2
          GENERIC,   PUBLIC  :: pluseq => psi1_pluseq_psi2

          PROCEDURE, PRIVATE, PASS(psi)  :: CplxPsi_TO_RCpsi
          PROCEDURE, PRIVATE, PASS(psi)  :: RCPsi_TO_CplxPsi
          PROCEDURE, PRIVATE, PASS(psi1) :: psi2TOpsi1
          PROCEDURE, PRIVATE, PASS(psi)  :: R_TOpsi
          PROCEDURE, PRIVATE, PASS(psi)  :: C_TOpsi
          GENERIC,   PUBLIC  :: assignment(=) => psi2TOpsi1,            &
                                    CplxPsi_TO_RCpsi, RCPsi_TO_CplxPsi, &
                                    R_TOpsi, C_TOpsi

          PROCEDURE, PRIVATE             :: psi1_plus_psi2
          PROCEDURE, PRIVATE, PASS(psi)  :: R_plus_psi
          PROCEDURE, PRIVATE, PASS(psi)  :: psi_plus_R
          PROCEDURE, PRIVATE, PASS(psi)  :: C_plus_psi
          PROCEDURE, PRIVATE, PASS(psi)  :: psi_plus_C
          PROCEDURE, PRIVATE, PASS(psi)  :: plus_psi
          GENERIC,   PUBLIC  :: operator(+) => psi1_plus_psi2,          &
                                               R_plus_psi, psi_plus_R,  &
                                               C_plus_psi, psi_plus_C,  &
                                               plus_psi

          PROCEDURE, PRIVATE             :: psi1_minus_psi2
          PROCEDURE, PRIVATE, PASS(psi)  :: R_minus_psi
          PROCEDURE, PRIVATE, PASS(psi)  :: psi_minus_R
          PROCEDURE, PRIVATE, PASS(psi)  :: C_minus_psi
          PROCEDURE, PRIVATE, PASS(psi)  :: psi_minus_C
          PROCEDURE, PRIVATE, PASS(psi)  :: minus_psi
          GENERIC,   PUBLIC  :: operator(-) => psi1_minus_psi2,         &
                                              R_minus_psi, psi_minus_R, &
                                              C_minus_psi, psi_minus_C, &
                                              minus_psi

          PROCEDURE, PRIVATE, PASS(psi)  :: R_time_psi
          PROCEDURE, PRIVATE, PASS(psi)  :: psi_time_R
          PROCEDURE, PRIVATE, PASS(psi)  :: C_time_psi
          PROCEDURE, PRIVATE, PASS(psi)  :: psi_time_C
          GENERIC,   PUBLIC  :: operator(*) => R_time_psi, psi_time_R,  &
                                               C_time_psi, psi_time_C

          PROCEDURE, PRIVATE, PASS(psi)  :: psi_over_R
          PROCEDURE, PRIVATE, PASS(psi)  :: psi_over_C
          GENERIC,   PUBLIC  :: operator(/) => psi_over_R, psi_over_C


          FINAL              :: dealloc_all_psi
        END TYPE param_psi

 PUBLIC :: param_psi,Set_psi_With_index,Set_Random_psi,copy_psi2TOpsi1
 PUBLIC :: ecri_init_psi,ecri_psi
 PUBLIC :: get_nb_be_FROM_psi,get_nb_bi_FROM_psi
 PUBLIC :: alloc_array,dealloc_array,alloc_NParray,dealloc_NParray
 PUBLIC :: alloc_psi,dealloc_psi

 PUBLIC :: get_RVec_OF_psi_AT_ind_a,get_CVec_OF_psi_AT_ind_a
 PUBLIC :: set_RVec_OF_psi_AT_ind_a,set_CVec_OF_psi_AT_ind_a
 PUBLIC :: get_RMat_OF_psi_AT_ind_a,get_CMat_OF_psi_AT_ind_a
 PUBLIC :: set_RMat_OF_psi_AT_ind_a,set_CMat_OF_psi_AT_ind_a
 !> subroutine for working on full Smolyak rep.
 PUBLIC :: psi_plus_SR_MPI,psi_minus_SR_MPI,psi_times_SR_MPI

      INTERFACE alloc_array
        MODULE PROCEDURE alloc_array_OF_Psidim1
      END INTERFACE

      INTERFACE dealloc_array
        MODULE PROCEDURE dealloc_array_OF_Psidim1
      END INTERFACE

      INTERFACE alloc_NParray
        MODULE PROCEDURE alloc_NParray_OF_Psidim1
      END INTERFACE

      INTERFACE dealloc_NParray
        MODULE PROCEDURE dealloc_NParray_OF_Psidim1
      END INTERFACE

      !!@description: TODO
      INTERFACE Set_psi_With_index
        MODULE PROCEDURE Set_psi_With_index_R
        MODULE PROCEDURE Set_psi_With_index_C
      END INTERFACE

      !> operators for working in Smolyak rep.
      INTERFACE psi_times_SR_MPI
        MODULE PROCEDURE psi_times_R_SR_MPI
        MODULE PROCEDURE psi_times_C_SR_MPI
        MODULE PROCEDURE psi_times_psi_SR_MPI
      END INTERFACE

      INTERFACE psi_plus_SR_MPI
        MODULE PROCEDURE psi_plus_R_SR_MPI
        MODULE PROCEDURE psi_plus_C_SR_MPI
        MODULE PROCEDURE psi_plus_psi_SR_MPI
      END INTERFACE

      INTERFACE psi_minus_SR_MPI
        MODULE PROCEDURE psi_minus_R_SR_MPI
        MODULE PROCEDURE psi_minus_C_SR_MPI
        MODULE PROCEDURE psi_minus_psi_SR_MPI
      END INTERFACE

 CONTAINS

!=======================================================================================
!     allocation of psi
!     deallocation of psi
!=======================================================================================
      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      SUBROUTINE alloc_psi(psi,BasisRep,GridRep,BasisRep_contrac)
      USE mod_MPI_aux

      TYPE (param_psi),  intent(inout) :: psi
      logical, optional, intent(in)    :: BasisRep,GridRep,BasisRep_contrac

      integer :: i,ib

!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING alloc_psi'
        write(out_unitp,*) ' 1 test_alloc_psi'
        write(out_unitp,*) 'BasisRep,GridRep',psi%BasisRep,psi%GridRep
      END IF
!-----------------------------------------------------------
      IF (present(BasisRep))          psi%BasisRep          = BasisRep
      IF (present(GridRep))           psi%GridRep           = GridRep


      IF (present(BasisRep_contrac)) THEN
        psi%BasisRep_contrac  = BasisRep_contrac
        IF (BasisRep_contrac) THEN
          psi%nb_tot = psi%nb_tot_contrac
        ELSE
          psi%nb_tot = psi%nb_tot_uncontrac
        END IF
      END IF

      IF (psi%BasisRep .AND. psi%nb_tot < 1) THEN
        write(out_unitp,*) ' ERROR in alloc_psi'
        write(out_unitp,*) ' nb_tot <1 ',psi%nb_tot
        STOP
      END IF
      IF (psi%GridRep .AND. psi%nb_qaie < 1) THEN
        write(out_unitp,*) ' ERROR in alloc_psi'
        write(out_unitp,*) ' nb_qaie <1 ',psi%nb_qaie
        STOP
      END IF

      IF (psi%BasisRep) THEN ! allocate psi%BasisRep
        IF (psi%cplx) THEN
          IF ( .NOT. allocated(psi%CvecB) ) THEN
            IF(keep_MPI) CALL alloc_NParray(psi%CvecB,(/psi%nb_tot/),'psi%CvecB','alloc_psi')
            IF (debug) write(out_unitp,*) 'alloc: CvecB'
            IF(keep_MPI) psi%CvecB(:) = CZERO
          END IF
          IF ( allocated(psi%RvecB) ) THEN
            CALL dealloc_NParray(psi%RvecB,'psi%RvecB','alloc_psi')
            IF (debug) write(out_unitp,*) 'dealloc: RvecB'
          END IF
        ELSE
          IF ( .NOT. allocated(psi%RvecB) ) THEN
            IF(keep_MPI) CALL alloc_NParray(psi%RvecB,(/psi%nb_tot/),'psi%RvecB','alloc_psi')
            IF (debug) write(out_unitp,*) 'alloc: RvecB'
            IF(keep_MPI) psi%RvecB(:) = ZERO
          END IF
          IF ( allocated(psi%CvecB) ) THEN
            CALL dealloc_NParray(psi%CvecB,'psi%CvecB','alloc_psi')
            IF (debug) write(out_unitp,*) 'dealloc: CvecB'
          END IF
        END IF
      ELSE ! deallocate psi%BasisRep
        IF ( allocated(psi%RvecB) ) THEN
          CALL dealloc_NParray(psi%RvecB,'psi%RvecB','alloc_psi')
          IF (debug) write(out_unitp,*) 'dealloc: RvecB'
        END IF
        IF ( allocated(psi%CvecB) ) THEN
          CALL dealloc_NParray(psi%CvecB,'psi%CvecB','alloc_psi')
          IF (debug) write(out_unitp,*) 'dealloc: CvecB'
        END IF
      END IF


      IF (psi%GridRep) THEN ! allocate psi%GridRep
        IF (psi%cplx) THEN
          IF ( .NOT. allocated(psi%CvecG) ) THEN
            IF(keep_MPI) CALL alloc_NParray(psi%CvecG,(/psi%nb_qaie/),'psi%CvecG','alloc_psi')
            IF (debug) write(out_unitp,*) 'alloc: CvecG'
            IF(keep_MPI) psi%CvecG(:) = CZERO
          END IF
          IF ( allocated(psi%RvecG) ) THEN
            CALL dealloc_NParray(psi%RvecG,'psi%RvecG','alloc_psi')
            IF (debug) write(out_unitp,*) 'dealloc: RvecG'
          END IF
        ELSE
          IF ( .NOT. allocated(psi%RvecG) ) THEN
            IF(keep_MPI) CALL alloc_NParray(psi%RvecG,(/psi%nb_qaie/),'psi%RvecG','alloc_psi')
            IF (debug) write(out_unitp,*) 'alloc: RvecG'
            IF(keep_MPI) psi%RvecG(:) = ZERO
          END IF
          IF ( allocated(psi%CvecG) ) THEN
            CALL dealloc_NParray(psi%CvecG,'psi%CvecG','alloc_psi')
            IF (debug) write(out_unitp,*) 'dealloc: CvecG'
          END IF
        END IF
      ELSE ! deallocate psi%GridRep
        IF ( allocated(psi%RvecG) ) THEN
          IF(keep_MPI) CALL dealloc_NParray(psi%RvecG,'psi%RvecG','alloc_psi')
          IF (debug) write(out_unitp,*) 'dealloc: RvecG'
        END IF
        IF ( allocated(psi%CvecG) ) THEN
          IF(keep_MPI) CALL dealloc_NParray(psi%CvecG,'psi%CvecG','alloc_psi')
          IF (debug) write(out_unitp,*) 'dealloc: CvecG'
        END IF
      END IF

      if ( psi%nb_TDParam > 0 .AND. .NOT. allocated(psi%TDParam)) then
        CALL alloc_NParray(psi%TDParam,[psi%nb_TDParam],'psi%TDParam','alloc_psi')
      end if

      IF(psi%SRG_MPI) THEN
        IF(.NOT. allocated(psi%SR_G_length)) THEN
          CALL allocate_array(psi%SR_G_length,0,MPI_np-1)
        ENDIF

        IF(.NOT. allocated(psi%SR_G_index)) THEN
          CALL allocate_array(psi%SR_G_index,iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)+1)
        ENDIF

        IF (.NOT. allocated(psi%SR_G)) THEN
          IF(psi%cplx) THEN
            CALL alloc_NParray(psi%SR_G,(/psi%SR_G_length(MPI_id),2/),'psi%SR_G','alloc_psi')
            !CALL allocate_array(psi%SR_G,1,psi%SR_G_length(MPI_id),1,2)
          ELSE
            CALL alloc_NParray(psi%SR_G,(/psi%SR_G_length(MPI_id),1/),'psi%SR_G','alloc_psi')
            !CALL allocate_array(psi%SR_G,1,psi%SR_G_length(MPI_id),1,1)
          ENDIF
        ENDIF
        ! consider deallocate RvecB,RvecG ...
      ENDIF

      IF(psi%SRB_MPI) THEN
        IF(.NOT. allocated(psi%SR_B_length)) THEN
          CALL allocate_array(psi%SR_B_length,0,MPI_np-1)
        ENDIF

        IF(.NOT. allocated(psi%SR_B_index)) THEN
          CALL allocate_array(psi%SR_B_index,iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)+1)
        ENDIF

        IF (.NOT. allocated(psi%SR_B)) THEN
          IF(psi%cplx) THEN
            CALL alloc_NParray(psi%SR_B,(/psi%SR_B_length(MPI_id),2/),'psi%SR_B','alloc_psi')
          ELSE
            CALL alloc_NParray(psi%SR_B,(/psi%SR_B_length(MPI_id),1/),'psi%SR_B','alloc_psi')
          ENDIF
        ENDIF
        ! consider deallocate RvecB,RvecG ...
      ENDIF

      IF (debug) write(out_unitp,*) ' 2 test_alloc_psi'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END alloc_psi'
      END IF
!-----------------------------------------------------------

      END SUBROUTINE alloc_psi
!=======================================================================================

 FUNCTION print_alloc_psi(psi)
 USE mod_system
 IMPLICIT NONE

  character(len=:), allocatable     :: print_alloc_psi

  TYPE (param_psi), intent(in) :: psi


  print_alloc_psi = String_TO_String(                                           &
          'RvecB:' // logical_TO_char(allocated(psi%RvecB)) //                  &
         ',RvecG:' // logical_TO_char(allocated(psi%RvecG)) //                  &
         ',CvecB:' // logical_TO_char(allocated(psi%CvecB)) //                  &
         ',CvecG:' // logical_TO_char(allocated(psi%CvecG)) //                  &
         ',TDParam:' // logical_TO_char(allocated(psi%TDParam))   )

 END FUNCTION print_alloc_psi

      SUBROUTINE alloc_vecB_OF_psi(psi)

      TYPE (param_psi)   :: psi
      integer :: i,ib

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='alloc_vecB_OF_psi'
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'BasisRep,GridRep',psi%BasisRep,psi%GridRep
      END IF
!-----------------------------------------------------------

      IF (psi%BasisRep .AND. psi%nb_tot < 1) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' nb_tot <1 ',psi%nb_tot
        STOP
      END IF

      IF (psi%cplx) THEN
        IF ( .NOT. allocated(psi%CvecB) ) THEN
          CALL alloc_NParray(psi%CvecB,(/psi%nb_tot/),'psi%CvecB',name_sub)
          IF (debug) write(out_unitp,*) 'alloc: CvecB'
          psi%CvecB(:) = CZERO
        END IF
        IF ( allocated(psi%RvecB) ) THEN
          CALL dealloc_NParray(psi%RvecB,'psi%RvecB','name_sub')
          IF (debug) write(out_unitp,*) 'dealloc: RvecB'
        END IF
      ELSE
        IF ( .NOT. allocated(psi%RvecB) ) THEN
          CALL alloc_NParray(psi%RvecB,(/psi%nb_tot/),'psi%RvecB',name_sub)
          IF (debug) write(out_unitp,*) 'alloc: RvecB'
          psi%RvecB(:) = ZERO
        END IF
        IF ( allocated(psi%CvecB) ) THEN
          CALL dealloc_NParray(psi%CvecB,'psi%CvecB','name_sub')
          IF (debug) write(out_unitp,*) 'dealloc: CvecB'
        END IF
      END IF


!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------


      END SUBROUTINE alloc_vecB_OF_psi


      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
  SUBROUTINE dealloc_all_psi(psi)
      TYPE (param_psi)   :: psi

      !----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='dealloc_all_psi'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      !-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL flush_perso(out_unitp)
      END IF
      !-----------------------------------------------------------

      CALL dealloc_psi(psi,delete_all=.TRUE.)
      psi%BasisRep        = .TRUE.
      psi%GridRep         = .FALSE.

      nullify(psi%BasisnD)
      nullify(psi%Basis2n)

      !-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
      !-----------------------------------------------------------

  END SUBROUTINE dealloc_all_psi
      SUBROUTINE dealloc_psi(psi,delete_all)


!----- variables for the WP propagation ----------------------------
      TYPE (param_psi)   :: psi
      logical, optional  :: delete_all

      integer            :: i
      integer,save       :: nb = 0
      logical            :: delete_all_loc


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='dealloc_psi'
      logical, parameter :: debug=.FALSE.
!      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) ' 1 test_',name_sub
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

      delete_all_loc = .FALSE.
      IF (present(delete_all)) delete_all_loc = delete_all
      nb = nb+1
      IF (debug) THEN
        write(out_unitp,*) 'dealloc_psi:',nb
        CALL flush_perso(out_unitp)
      END IF

      IF (allocated(psi%CvecB)) THEN
        CALL dealloc_NParray(psi%CvecB,'psi%CvecB',name_sub)
        IF (debug) THEN
          write(out_unitp,*) 'dealloc: CvecB'
          CALL flush_perso(out_unitp)
        END IF
      END IF

      IF (allocated(psi%RvecB)) THEN
        CALL dealloc_NParray(psi%RvecB,'psi%RvecB',name_sub)
        IF (debug) THEN
          write(out_unitp,*) 'dealloc: RvecB'
          CALL flush_perso(out_unitp)
        END IF
      END IF

      IF (allocated(psi%CvecG)) THEN
        CALL dealloc_NParray(psi%CvecG,'psi%CvecG',name_sub)
        IF (debug) THEN
          write(out_unitp,*) 'dealloc: CvecG'
          CALL flush_perso(out_unitp)
        END IF
      END IF

      IF (allocated(psi%RvecG)) THEN
        CALL dealloc_NParray(psi%RvecG,'psi%RvecG',name_sub)
        IF (debug) THEN
          write(out_unitp,*) 'dealloc: RvecG'
          CALL flush_perso(out_unitp)
        END IF
      END IF

      if (allocated(psi%TDParam)) then
        CALL dealloc_NParray(psi%TDParam,'psi%TDParam',name_sub)
      end if

      IF(allocated(psi%SR_G)) THEN
        CALL dealloc_NParray(psi%SR_G,'psi%SR_G',name_sub)
        IF (debug) THEN
          write(out_unitp,*) 'dealloc: SR_G'
          CALL flush_perso(out_unitp)
        END IF
      ENDIF

      IF(allocated(psi%SR_B)) THEN
        CALL dealloc_NParray(psi%SR_B,'psi%SR_B',name_sub)
        IF (debug) THEN
          write(out_unitp,*) 'dealloc: SR_B'
          CALL flush_perso(out_unitp)
        END IF
      ENDIF

      psi%symab = -1


      IF (delete_all_loc) THEN
        IF (debug) write(out_unitp,*) 'dealloc: init param'

        nullify(psi%Basis2n)
        nullify(psi%BasisnD)
        nullify(psi%para_AllBasis)

        psi%init            = .FALSE.
        psi%cplx            = .FALSE.
        !psi%BasisRep        = .TRUE.
        !psi%GridRep         = .FALSE.

        psi%BasisRep_contrac = .FALSE.

        psi%nb_tot           = 0
        psi%nb_tot_contrac   = 0
        psi%nb_tot_uncontrac = 0

        psi%nb_paie          = 0
        psi%nb_baie          = 0
        psi%nb_ba            = 0
        psi%nb_bi            = 0
        psi%nb_be            = 0
        psi%nb_bRot          = 0
        psi%nb_qa            = 0
        psi%nb_qaie          = 0

        psi%nb_act1          = 0
        psi%nb_act           = 0
        psi%nb_basis_act1    = 0
        psi%nb_basis         = 0
        psi%max_dim          = 0

        psi%IndAvOp          = -1
        psi%CAvOp            = (ZERO,ZERO)
        psi%convAvOp         = .FALSE.

        IF (debug) write(out_unitp,*) 'end dealloc: init param'
        CALL flush_perso(out_unitp)


      END IF


!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) ' 2 test_',name_sub
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

      END SUBROUTINE dealloc_psi

      SUBROUTINE alloc_array_OF_Psidim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      TYPE (param_psi), pointer, intent(inout) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=1
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_Psidim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (associated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)

         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1)),stat=err_mem)
       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'param_psi')

      END SUBROUTINE alloc_array_OF_Psidim1
      SUBROUTINE dealloc_array_OF_Psidim1(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (param_psi),  pointer, intent(inout) :: tab(:)
      character (len=*),          intent(in)    :: name_var,name_sub

      integer :: i
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_Psidim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN
       IF (.NOT. associated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       DO i=lbound(tab,dim=1),ubound(tab,dim=1)
         CALL dealloc_psi(tab(i))
       END DO

       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'param_psi')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_Psidim1

      SUBROUTINE alloc_NParray_OF_Psidim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      TYPE (param_psi), allocatable, intent(inout)        :: tab(:)
      integer,                       intent(in)           :: tab_ub(:)
      integer,                       intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=1
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_NParray_OF_Psidim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (allocated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)

         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1)),stat=err_mem)
       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'param_psi')

      END SUBROUTINE alloc_NParray_OF_Psidim1
      SUBROUTINE dealloc_NParray_OF_Psidim1(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (param_psi),  allocatable, intent(inout) :: tab(:)
      character (len=*),              intent(in)    :: name_var,name_sub

      integer :: i
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_NParray_OF_Psidim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. allocated(tab)) RETURN
       IF (.NOT. allocated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       DO i=lbound(tab,dim=1),ubound(tab,dim=1)
         CALL dealloc_psi(tab(i))
       END DO

       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'param_psi')

      END SUBROUTINE dealloc_NParray_OF_Psidim1

      FUNCTION get_nb_be_FROM_psi(psi)
      USE mod_system
      IMPLICIT NONE

      integer                      :: get_nb_be_FROM_psi
      TYPE (param_psi), intent(in) :: psi

      integer                      :: nb_be_psi,nb_be_NewBasisEl

!----- for debuging --------------------------------------------------
      !logical, parameter :: debug = .TRUE.
       logical, parameter :: debug = .FALSE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING get_nb_be_FROM_psi'
      END IF

     IF (NewBasisEl) THEN
       nb_be_NewBasisEl = get_nb_be_FROM_basis(psi%BasisnD)
       IF (nb_be_NewBasisEl == -1) THEN
         write(out_unitp,*) ' ERROR in get_nb_be_FROM_psi'
         write(out_unitp,*) ' NewBasisEl = t and no El basis set!!'
         write(out_unitp,*) ' CHECK the source'
         STOP
       END IF
     ELSE
       nb_be_NewBasisEl = -1
     END IF

     nb_be_psi = psi%nb_be


     IF (NewBasisEl .AND. nb_be_psi /= nb_be_NewBasisEl ) THEN
       write(out_unitp,*) ' ERROR in get_nb_be_FROM_psi'
       write(out_unitp,*) ' inconsistent value in psi and NewBasisEl'
       write(out_unitp,*) ' nb_be in psi        ',nb_be_psi
       write(out_unitp,*) ' nb_be in NewBasisEl ',nb_be_NewBasisEl
       write(out_unitp,*) ' CHECK the source'
       STOP
     END IF
     get_nb_be_FROM_psi = nb_be_psi

      IF (debug) THEN
        write(out_unitp,*) 'END get_nb_be_FROM_psi'
      END IF
      END FUNCTION get_nb_be_FROM_psi
      FUNCTION get_nb_bi_FROM_psi(psi)
      USE mod_system
      IMPLICIT NONE

      integer                      :: get_nb_bi_FROM_psi
      TYPE (param_psi), intent(in) :: psi

      integer                      :: nb_bi_Basis2n,nb_bi_psi

!----- for debuging --------------------------------------------------
      !logical, parameter :: debug = .TRUE.
       logical, parameter :: debug = .FALSE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING get_nb_bi_FROM_psi'
      END IF

     nb_bi_psi = psi%nb_bi

!     nb_bi_Basis2n = get_nb_FROM_basis(psi%Basis2n)

!     IF (nb_bi_psi /= nb_bi_Basis2n) THEN
!       write(out_unitp,*) ' ERROR in get_nb_bi_FROM_psi'
!       write(out_unitp,*) ' inconsistent value in psi and psi%Basis2n'
!       write(out_unitp,*) ' nb_bi in psi        ',nb_bi_psi
!       write(out_unitp,*) ' nb_bi in psi%Basis2n',nb_bi_Basis2n
!       write(out_unitp,*) ' CHECK the source'
!       STOP
!     END IF

     get_nb_bi_FROM_psi = nb_bi_psi

      IF (debug) THEN
        write(out_unitp,*) 'END get_nb_bi_FROM_psi'
      END IF
      END FUNCTION get_nb_bi_FROM_psi
!================================================================
!
!     copy the initialization parameters
!
!================================================================
      !!@description: y the initialization parameters
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      SUBROUTINE psi2TOpsi1(psi1,psi2)

!----- variables for the WP propagation ----------------------------
      CLASS (param_psi), intent(inout) :: psi1
      TYPE (param_psi),  intent(in)    :: psi2


      integer  :: i,n1

      !write(out_unitp,*) 'BEGINNING psi2TOpsi1'

      CALL dealloc_psi(psi1)

      IF (associated(psi2%para_AllBasis)) THEN
         psi1%para_AllBasis => psi2%para_AllBasis
      ELSE
         write(out_unitp,*) ' ERROR in psi2TOpsi1'
         write(out_unitp,*) ' para_AllBasis CANNOT be associated'
         write(out_unitp,*) ' asso psi2%para_AllBasis',associated(psi2%para_AllBasis)
         write(out_unitp,*) ' CHECK the source'
         STOP
      END IF

      IF (associated(psi1%para_AllBasis%BasisnD)) THEN
         psi1%BasisnD => psi1%para_AllBasis%BasisnD
      ELSE
         write(out_unitp,*) ' ERROR in psi2TOpsi1'
         write(out_unitp,*) ' BasisnD CANNOT be associated'
         write(out_unitp,*) ' asso psi1%para_AllBasis%BasisnD',associated(psi1%para_AllBasis%BasisnD)
         write(out_unitp,*) ' CHECK the source'
         STOP
      END IF

      IF (associated(psi1%para_AllBasis%Basis2n)) THEN
         psi1%Basis2n => psi1%para_AllBasis%Basis2n
      ELSE
         write(out_unitp,*) ' ERROR in psi2TOpsi1'
         write(out_unitp,*) ' BasisnD CANNOT be associated'
         write(out_unitp,*) ' asso psi1%para_AllBasis%Basis2n',associated(psi1%para_AllBasis%Basis2n)
         write(out_unitp,*) ' CHECK the source'
         STOP
      END IF

      psi1%init             = psi2%init
      psi1%cplx             = psi2%cplx

      psi1%BasisRep         = psi2%BasisRep
      psi1%GridRep          = psi2%GridRep
      psi1%BasisRep_contrac = psi2%BasisRep_contrac

      psi1%SRG_MPI          = psi2%SRG_MPI      !< if working on full smolyak rep. on grid
      psi1%SRB_MPI          = psi2%SRB_MPI     !< if working on full smolyak rep. on Basis

      psi1%symab            = psi2%symab

      psi1%nb_tot           = psi2%nb_tot
      psi1%nb_tot_contrac   = psi2%nb_tot_contrac
      psi1%nb_tot_uncontrac = psi2%nb_tot_uncontrac

      psi1%nb_baie          = psi2%nb_baie
      psi1%nb_ba            = psi2%nb_ba
      psi1%nb_bi            = psi2%nb_bi
      psi1%nb_be            = psi2%nb_be
      psi1%nb_bRot          = psi2%nb_bRot

      psi1%nb_qa            = psi2%nb_qa
      psi1%nb_qaie          = psi2%nb_qaie

      psi1%nb_act1          = psi2%nb_act1
      psi1%nb_act           = psi2%nb_act
      psi1%nb_basis_act1    = psi2%nb_basis_act1
      psi1%nb_basis         = psi2%nb_basis

      psi1%max_dim          = psi2%max_dim

      psi1%nb_TDParam       = psi2%nb_TDParam


      IF (psi1%nb_tot <= 0) THEN
        write(out_unitp,*) ' ERROR in psi2TOpsi1'
        write(out_unitp,*) 'nb_tot = 0',psi1%nb_tot
        STOP
      END IF
      IF (psi1%nb_baie <= 0 .OR. psi1%nb_qaie <= 0) THEN
        write(out_unitp,*) ' ERROR in psi2TOpsi1'
        write(out_unitp,*) 'nb_baie OR nb_qaie <= 0',psi1%nb_baie,psi1%nb_qaie
        STOP
      END IF

      CALL alloc_psi(psi1)

      IF (psi1%BasisRep) THEN
        IF (allocated(psi2%CvecB)) psi1%CvecB(:) = psi2%CvecB(:)
        IF (allocated(psi2%RvecB)) psi1%RvecB(:) = psi2%RvecB(:)
      END IF
      IF (psi1%GridRep) THEN
        IF (allocated(psi2%CvecG)) psi1%CvecG(:) = psi2%CvecG(:)
        IF (allocated(psi2%RvecG)) psi1%RvecG(:) = psi2%RvecG(:)
      END IF

      IF (allocated(psi2%TDParam)) psi1%TDParam(:) = psi2%TDParam(:)

      IF(psi1%SRG_MPI) THEN
        IF(allocated(psi2%SR_G))        psi1%SR_G=psi2%SR_G
        IF(allocated(psi2%SR_G_length)) psi1%SR_G_length=psi2%SR_G_length
        IF(allocated(psi2%SR_G_index))  psi1%SR_G_index=psi2%SR_G_index
      ENDIF
      IF(psi1%SRB_MPI) THEN
        IF(allocated(psi2%SR_B))        psi1%SR_B=psi2%SR_B
        IF(allocated(psi2%SR_B_length)) psi1%SR_B_length=psi2%SR_B_length
        IF(allocated(psi2%SR_B_index))  psi1%SR_B_index=psi2%SR_B_index
      ENDIF


     psi1%IndAvOp         = psi2%IndAvOp
     psi1%CAvOp           = psi2%CAvOp
     psi1%convAvOp        = psi2%convAvOp

      psi1%norm2          = psi2%norm2

      END SUBROUTINE psi2TOpsi1

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
!=======================================================================================
      SUBROUTINE copy_psi2TOpsi1(psi1,psi2,BasisRep,GridRep,alloc)

!----- variables for the WP propagation ----------------------------
      TYPE (param_psi),   intent(inout)         :: psi1
      TYPE (param_psi),   intent(in)            :: psi2
      logical,            intent(in), optional  :: BasisRep
      logical,            intent(in), optional  :: GridRep
      logical,            intent(in), optional  :: alloc


      logical  :: loc_alloc
      integer  :: i,n1

!     write(out_unitp,*) 'BEGINNING copy_psi2TOpsi1'

      CALL dealloc_psi(psi1)

      IF (associated(psi2%para_AllBasis)) THEN
         psi1%para_AllBasis => psi2%para_AllBasis
      ELSE
         write(out_unitp,*) ' ERROR in copy_psi2TOpsi1'
         write(out_unitp,*) ' para_AllBasis CANNOT be associated'
         write(out_unitp,*) ' asso psi2%para_AllBasis',associated(psi2%para_AllBasis)
         write(out_unitp,*) ' CHECK the source'
         STOP
      END IF

      IF (associated(psi1%para_AllBasis%BasisnD)) THEN
         psi1%BasisnD => psi1%para_AllBasis%BasisnD
      ELSE
         write(out_unitp,*) ' ERROR in copy_psi2TOpsi1'
         write(out_unitp,*) ' BasisnD CANNOT be associated'
         write(out_unitp,*) ' asso psi1%para_AllBasis%BasisnD',associated(psi1%para_AllBasis%BasisnD)
         write(out_unitp,*) ' CHECK the source'
         STOP
      END IF

      IF (associated(psi1%para_AllBasis%Basis2n)) THEN
         psi1%Basis2n => psi1%para_AllBasis%Basis2n
      ELSE
         write(out_unitp,*) ' ERROR in copy_psi2TOpsi1'
         write(out_unitp,*) ' BasisnD CANNOT be associated'
         write(out_unitp,*) ' asso psi1%para_AllBasis%Basis2n',associated(psi1%para_AllBasis%Basis2n)
         write(out_unitp,*) ' CHECK the source'
         STOP
      END IF


      psi1%init = psi2%init
      psi1%cplx = psi2%cplx

      IF (present(alloc)) THEN
        loc_alloc  = alloc
      ELSE
        loc_alloc  = .TRUE.
      END IF

      IF (present(BasisRep)) THEN
        psi1%BasisRep  = BasisRep
      ELSE
        psi1%BasisRep  = psi2%BasisRep
      END IF
      IF (present(GridRep)) THEN
        psi1%GridRep  = GridRep
      ELSE
        psi1%GridRep  = psi2%GridRep
      END IF

      psi1%BasisRep_contrac = psi2%BasisRep_contrac


      psi1%SRG_MPI=psi2%SRG_MPI
      psi1%SRB_MPI=psi2%SRB_MPI

      psi1%symab = psi2%symab


      psi1%nb_tot           = psi2%nb_tot
      psi1%nb_tot_contrac   = psi2%nb_tot_contrac
      psi1%nb_tot_uncontrac = psi2%nb_tot_uncontrac

      psi1%nb_baie = psi2%nb_baie
      psi1%nb_ba   = psi2%nb_ba
      psi1%nb_bi   = psi2%nb_bi
      psi1%nb_be   = psi2%nb_be
      psi1%nb_bRot = psi2%nb_bRot
      psi1%nb_qa   = psi2%nb_qa
      psi1%nb_qaie = psi2%nb_qaie

      psi1%nb_act1       = psi2%nb_act1
      psi1%nb_act        = psi2%nb_act
      psi1%nb_basis_act1 = psi2%nb_basis_act1
      psi1%nb_basis      = psi2%nb_basis

      psi1%max_dim       = psi2%max_dim

      psi1%nb_TDParam    = psi2%nb_TDParam


      IF (psi1%nb_tot <= 0) THEN
        write(out_unitp,*) ' ERROR in copy_psi2TOpsi1'
        write(out_unitp,*) 'nb_tot = 0',psi1%nb_tot
        STOP
      END IF
      IF (psi1%nb_baie <= 0 .OR. psi1%nb_qaie <= 0) THEN
        write(out_unitp,*) ' ERROR in copy_psi2TOpsi1'
        write(out_unitp,*) 'nb_baie OR nb_qaie = 0',psi1%nb_baie,psi1%nb_qaie
        STOP
      END IF

      IF (loc_alloc) THEN
        CALL alloc_psi(psi1)

        IF (allocated(psi2%TDParam)) psi1%TDParam(:) = psi2%TDParam(:)

        IF (psi1%BasisRep) THEN

          IF (allocated(psi2%RvecB)) THEN
            psi1%RvecB(:) = psi2%RvecB(:)
          END IF
          IF (allocated(psi2%CvecB)) THEN
            psi1%CvecB(:) = psi2%CvecB(:)
          END IF
        END IF

        IF (psi1%GridRep) THEN
          IF (allocated(psi2%RvecG)) THEN
            psi1%RvecG(:) = psi2%RvecG(:)
          END IF
          IF (allocated(psi2%CvecG)) THEN
            psi1%CvecG(:) = psi2%CvecG(:)
          END IF
        END IF

        IF(psi1%SRG_MPI) THEN
          psi1%SR_G_length=psi2%SR_G_length
          IF(allocated(psi2%SR_G)) psi1%SR_G=psi2%SR_G
        ENDIF

        IF(psi1%SRB_MPI) THEN
          psi1%SR_B_length=psi2%SR_B_length
          IF(allocated(psi2%SR_B)) psi1%SR_B=psi2%SR_B
        ENDIF

        psi1%IndAvOp         = psi2%IndAvOp
        psi1%CAvOp           = psi2%CAvOp
        psi1%convAvOp        = psi2%convAvOp

        psi1%norm2           = psi2%norm2

      END IF

!     write(out_unitp,*) 'END copy_psi2TOpsi1'

      END SUBROUTINE copy_psi2TOpsi1
!=======================================================================================

      SUBROUTINE CplxPsi_TO_RCpsi(RCPsi,Psi)
      USE mod_MPI

!----- variables for the WP propagation ----------------------------
      TYPE (param_psi),  intent(inout) :: RCPsi(2)
      CLASS (param_psi), intent(in)    :: psi

!     write(out_unitp,*) 'BEGINNING CplxPsi_TO_RCpsi'
      IF (.NOT. Psi%cplx) STOP 'ERROR psi MUST be complex'

      CALL copy_psi2TOpsi1(RCPsi(1),psi,alloc=.FALSE.)
      RCPsi(1)%cplx = .FALSE.
      CALL alloc_psi(RCPsi(1),psi%BasisRep,psi%GridRep)

      CALL copy_psi2TOpsi1(RCPsi(2),psi,alloc=.FALSE.)
      RCPsi(2)%cplx = .FALSE.
      CALL alloc_psi(RCPsi(2),psi%BasisRep,psi%GridRep)

      IF (psi%BasisRep .AND. allocated(psi%CvecB)) THEN
         RCPsi(1)%RvecB(:) = real(psi%CvecB(:),kind=Rkind)
         RCPsi(2)%RvecB(:) = aimag(psi%CvecB(:))
      END IF

      IF (psi%GridRep .AND. allocated(psi%CvecG)) THEN
         RCPsi(1)%RvecG(:) = real(psi%CvecG(:),kind=Rkind)
         RCPsi(2)%RvecG(:) = aimag(psi%CvecG(:))
      END IF



!     write(out_unitp,*) 'END CplxPsi_TO_RCpsi'

      END SUBROUTINE CplxPsi_TO_RCpsi
      SUBROUTINE RCPsi_TO_CplxPsi(Psi,RCPsi)

!----- variables for the WP propagation ----------------------------
      TYPE (param_psi),  intent(in)    :: RCPsi(2)
      CLASS (param_psi), intent(inout) :: psi

!     write(out_unitp,*) 'BEGINNING RCPsi_TO_CplxPsi'

      CALL copy_psi2TOpsi1(psi,RCPsi(1),alloc=.FALSE.)
      psi%cplx = .TRUE.
      CALL alloc_psi(psi,RCPsi(1)%BasisRep,RCPsi(1)%GridRep)

      IF (RCPsi(1)%BasisRep .AND. allocated(RCPsi(1)%RvecB) .AND.       &
          RCPsi(2)%BasisRep .AND. allocated(RCPsi(2)%RvecB)) THEN
         psi%CvecB(:) = cmplx(RCPsi(1)%RvecB,RCPsi(2)%RvecB,kind=Rkind)
      END IF

      IF (RCPsi(1)%GridRep .AND. allocated(RCPsi(1)%RvecG) .AND.        &
          RCPsi(2)%GridRep .AND. allocated(RCPsi(2)%RvecG)) THEN
         psi%CvecG(:) = cmplx(RCPsi(1)%RvecG,RCPsi(2)%RvecG,kind=Rkind)
      END IF

!     write(out_unitp,*) 'END RCPsi_TO_CplxPsi'

      END SUBROUTINE RCPsi_TO_CplxPsi

  SUBROUTINE get_RVec_OF_psi_AT_ind_a(RVec_AT_iq,psi,iq,ib,OldPara)

!----- variables for the WP propagation ----------------------------
      TYPE (param_psi),intent(in)              :: psi
      real(kind=Rkind),intent(inout)           :: RVec_AT_iq(:)
      integer,         intent(in),    optional :: iq,ib
      TYPE (OldParam), intent(inout), optional :: OldPara

      integer :: i_qaie,i_baie,i_bi,i_be,i_bie,iSG,iqSG,nqSG,nb0
      integer :: err_sub
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='get_RVec_OF_psi_AT_ind_a'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) write(out_unitp,*) 'BEGINNING ',name_sub

      IF (      present(iq) .AND.       present(ib) .OR.                       &
          .NOT. present(iq) .AND. .NOT. present(ib))     THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) ' present(iq),present(ib)',present(iq),present(ib)
        write(out_unitp,*) ' One and only one ib or iq must be present'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF

      IF (present(iq) .AND. .NOT. allocated(psi%RvecG)) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) '  psi%RvecG is not allocated !!'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF

      IF (present(ib) .AND. .NOT. allocated(psi%RvecB)) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) '  psi%RvecB is not allocated !!'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF


      IF (any(shape(RVec_AT_iq) /= [psi%nb_bi*psi%nb_be])) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) '  Inconsitent shape of RVec_AT_iq'
        write(out_unitp,*) '  shape(RVec_AT_iq) ',shape(RVec_AT_iq)
        write(out_unitp,*) '  psi%nb_bi*psi%nb_be',psi%nb_bi*psi%nb_be
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF

    RVec_AT_iq(:) = ZERO

    IF (psi%BasisnD%SparseGrid_type == 4 .AND. present(iq)) THEN
      IF (present(OldPara)) THEN
        CALL get_iqSG_iSG_FROM_iq(iSG,iqSG,iq,psi%BasisnD%para_SGType2,OldPara,err_sub)
      ELSE
        CALL get_iqSG_iSG_FROM_iq(iSG,iqSG,iq,psi%BasisnD%para_SGType2,err_sub=err_sub)
      END if
      IF (err_sub /= 0) STOP 'Error in get_iqSG_iSG_FROM_iq'
      nqSG = psi%BasisnD%para_SGType2%tab_nq_OF_SRep(iSG)
      nb0  = psi%BasisnD%para_SGType2%nb0 ! this MUST be equal to psi%nb_bi*psi%nb_be
      IF (iSG == 1) THEN
         i_qaie = iqSG
      ELSE
         i_qaie = iqSG + psi%BasisnD%para_SGType2%tab_Sum_nq_OF_SRep(iSG-1) * &
                                                 psi%BasisnD%para_SGType2%nb0
      END IF
      RVec_AT_iq(:) = psi%RvecG(i_qaie:i_qaie-1+nqSG*nb0:nqSG)

    ELSE

      IF (present(iq)) THEN
        DO i_be=1,psi%nb_be
        DO i_bi=1,psi%nb_bi
          i_bie  = (i_bi-1)+ (i_be-1)*psi%nb_bi
          i_qaie = iq + i_bie * psi%nb_qa

          RVec_AT_iq(i_bie+1) = psi%RvecG(i_qaie)

        END DO
        END DO
      ELSE ! ib is present
        DO i_be=1,psi%nb_be
        DO i_bi=1,psi%nb_bi
          i_bie  = (i_bi-1)+ (i_be-1)*psi%nb_bi
          i_baie = ib + i_bie * psi%nb_ba

          RVec_AT_iq(i_bie+1) = psi%RvecB(i_baie)

        END DO
        END DO
      END IF
    END IF

      IF (debug) THEN
        write(out_unitp,*) 'RVec_AT_iq ',RVec_AT_iq
        write(out_unitp,*) 'END ',name_sub
        flush(out_unitp)
      END IF

  END SUBROUTINE get_RVec_OF_psi_AT_ind_a
  SUBROUTINE get_CVec_OF_psi_AT_ind_a(CVec_AT_iq,psi,iq,ib,OldPara)

!----- variables for the WP propagation ----------------------------
      TYPE (param_psi),   intent(in)             :: psi
      complex(kind=Rkind),intent(inout)          :: CVec_AT_iq(:)
      integer,            intent(in), optional   :: iq,ib
      TYPE (OldParam), intent(inout), optional   :: OldPara

      integer :: i_qaie,i_baie,i_bi,i_be,i_bie,iSG,iqSG,nqSG,nb0
      integer :: err_sub
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='get_CVec_OF_psi_AT_ind_a'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) write(out_unitp,*) 'BEGINNING ',name_sub

      IF (      present(iq) .AND.       present(ib) .OR.                       &
          .NOT. present(iq) .AND. .NOT. present(ib))     THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) ' present(iq),present(ib)',present(iq),present(ib)
        write(out_unitp,*) ' One and only one ib or iq must be present'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF

      IF (present(iq) .AND. .NOT. allocated(psi%CvecG)) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) '  psi%CvecG is not allocated !!'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF

      IF (present(ib) .AND. .NOT. allocated(psi%CvecB)) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) '  psi%CvecB is not allocated !!'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF


      IF (any(shape(CVec_AT_iq) /= [psi%nb_bi*psi%nb_be])) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) '  Inconsitent shape of CVec_AT_iq'
        write(out_unitp,*) '  shape(CVec_AT_iq) ',shape(CVec_AT_iq)
        write(out_unitp,*) '  psi%nb_bi*psi%nb_be',psi%nb_bi*psi%nb_be
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF

    CVec_AT_iq(:) = ZERO
    IF (psi%BasisnD%SparseGrid_type == 4 .AND. present(iq)) THEN
      ! not with ib since psi is not in Smolyak rep, we don't need that
     IF (present(OldPara)) THEN
       CALL get_iqSG_iSG_FROM_iq(iSG,iqSG,iq,psi%BasisnD%para_SGType2,OldPara,err_sub)
     ELSE
       CALL get_iqSG_iSG_FROM_iq(iSG,iqSG,iq,psi%BasisnD%para_SGType2,err_sub=err_sub)
     END if
     IF (err_sub /= 0) STOP 'Error in get_iqSG_iSG_FROM_iq'
     nqSG = psi%BasisnD%para_SGType2%tab_nq_OF_SRep(iSG)
     nb0  = psi%BasisnD%para_SGType2%nb0
     IF (iSG == 1) THEN
        i_qaie = iqSG
     ELSE
        i_qaie = iqSG + psi%BasisnD%para_SGType2%tab_Sum_nq_OF_SRep(iSG-1) * &
                                                psi%BasisnD%para_SGType2%nb0
     END IF

     ! DO i_be=1,psi%nb_be
     ! DO i_bi=1,psi%nb_bi
     !   i_bie  = (i_bi-1)+ (i_be-1)*psi%nb_bi
     !   CVec_AT_iq(i_bie) = psi%CvecG(i_qaie)
     !   i_qaie = i_qaie + nqSG
     ! END DO
     ! END DO
     CVec_AT_iq(:) = psi%CvecG(i_qaie:i_qaie-1+nqSG*nb0:nqSG)

    ELSE
      IF (present(iq)) THEN

        DO i_be=1,psi%nb_be
        DO i_bi=1,psi%nb_bi
          i_bie  = (i_bi-1)+ (i_be-1)*psi%nb_bi
          i_qaie = iq + i_bie * psi%nb_qa

          CVec_AT_iq(i_bie+1) = psi%CvecG(i_qaie)

        END DO
        END DO
      ELSE ! ib is present
        DO i_be=1,psi%nb_be
        DO i_bi=1,psi%nb_bi
          i_bie  = (i_bi-1)+ (i_be-1)*psi%nb_bi

          i_baie = ib + i_bie * psi%nb_ba

          CVec_AT_iq(i_bie+1) = psi%CvecB(i_baie)

        END DO
        END DO
      END IF
    END IF

    IF (debug) THEN
      write(out_unitp,*) 'CVec_AT_iq ',CVec_AT_iq
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END SUBROUTINE get_CVec_OF_psi_AT_ind_a

  SUBROUTINE set_RVec_OF_psi_AT_ind_a(RVec_AT_iq,psi,iq,ib,OldPara)

!----- variables for the WP propagation ----------------------------
    TYPE (param_psi),intent(inout)           :: psi
    real(kind=Rkind),intent(in)              :: RVec_AT_iq(:)
    integer,         intent(in), optional    :: iq,ib
    TYPE (OldParam), intent(inout), optional :: OldPara

    integer :: i_qaie,i_baie,i_bi,i_be,i_bie,iSG,iqSG,nqSG,nb0
    INTEGER :: err_sub
!---- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='set_RVec_OF_psi_AT_ind_a'
    logical, parameter :: debug=.FALSE.
    !logical, parameter :: debug=.TRUE.
!----------------------------------------------------------
    IF (debug) write(out_unitp,*) 'BEGINNING ',name_sub
    IF (debug) write(out_unitp,*) 'RVec_AT_iq ',RVec_AT_iq


    IF (      present(iq) .AND.       present(ib) .OR.                       &
        .NOT. present(iq) .AND. .NOT. present(ib))     THEN
      write(out_unitp,*) ' ERROR ',name_sub
      write(out_unitp,*) ' present(iq),present(ib)',present(iq),present(ib)
      write(out_unitp,*) ' One and only one ib or iq must be present'
      write(out_unitp,*) ' CHECK the fortran !!'
      STOP
    END IF

      IF (present(iq) .AND. .NOT. allocated(psi%RvecG)) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) '  psi%RvecG is not allocated !!'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF

      IF (present(ib) .AND. .NOT. allocated(psi%RvecB)) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) '  psi%RvecB is not allocated !!'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF


    IF (any(shape(RVec_AT_iq) /= [psi%nb_bi*psi%nb_be])) THEN
      write(out_unitp,*) ' ERROR ',name_sub
      write(out_unitp,*) '  Inconsitent shape of RVec_AT_iq'
      write(out_unitp,*) '  shape(RVec_AT_iq) ',shape(RVec_AT_iq)
      write(out_unitp,*) '  psi%nb_bi*psi%nb_be',psi%nb_bi*psi%nb_be
      write(out_unitp,*) ' CHECK the fortran !!'
      STOP
    END IF

    IF (psi%BasisnD%SparseGrid_type == 4 .AND. present(iq)) THEN
      ! not with ib since psi is not in Smolyak rep, we don't need that
     IF (present(OldPara)) THEN
       CALL get_iqSG_iSG_FROM_iq(iSG,iqSG,iq,psi%BasisnD%para_SGType2,OldPara,err_sub)
     ELSE
       CALL get_iqSG_iSG_FROM_iq(iSG,iqSG,iq,psi%BasisnD%para_SGType2,err_sub=err_sub)
     END if
     IF (err_sub /= 0) STOP 'Error in get_iqSG_iSG_FROM_iq'
     nqSG = psi%BasisnD%para_SGType2%tab_nq_OF_SRep(iSG)
     nb0  = psi%BasisnD%para_SGType2%nb0
     IF (iSG == 1) THEN
        i_qaie = iqSG
     ELSE
        i_qaie = iqSG + psi%BasisnD%para_SGType2%tab_Sum_nq_OF_SRep(iSG-1) * &
                                                psi%BasisnD%para_SGType2%nb0
     END IF

     psi%RvecG(i_qaie:i_qaie-1+nqSG*nb0:nqSG) = RVec_AT_iq

    ELSE
     IF (present(iq)) THEN

       DO i_be=1,psi%nb_be
       DO i_bi=1,psi%nb_bi
         i_bie  = (i_bi-1)+ (i_be-1)*psi%nb_bi
         i_qaie = iq + i_bie * psi%nb_qa

         psi%RvecG(i_qaie) = RVec_AT_iq(i_bie+1)

       END DO
       END DO
     ELSE ! ib is present
       DO i_be=1,psi%nb_be
       DO i_bi=1,psi%nb_bi
         i_bie  = (i_bi-1)+ (i_be-1)*psi%nb_bi
         i_baie = ib + i_bie * psi%nb_ba

         psi%RvecB(i_baie) = RVec_AT_iq(i_bie+1)

       END DO
       END DO
     END IF
    END IF

      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
        flush(out_unitp)
      END IF

  END SUBROUTINE set_RVec_OF_psi_AT_ind_a
  SUBROUTINE set_CVec_OF_psi_AT_ind_a(CVec_AT_iq,psi,iq,ib,OldPara)

!----- variables for the WP propagation ----------------------------
    TYPE (param_psi),   intent(inout)          :: psi
    complex(kind=Rkind),intent(in)             :: CVec_AT_iq(:)
    integer,            intent(in), optional   :: iq,ib
    TYPE (OldParam), intent(inout), optional   :: OldPara

    integer :: i_qaie,i_baie,i_bi,i_be,i_bie,iSG,iqSG,nqSG,nb0
    INTEGER :: err_sub
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='set_CVec_OF_psi_AT_ind_a'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) write(out_unitp,*) 'BEGINNING ',name_sub
      IF (debug) write(out_unitp,*) 'CVec_AT_iq ',CVec_AT_iq

      IF (      present(iq) .AND.       present(ib) .OR.                       &
          .NOT. present(iq) .AND. .NOT. present(ib))     THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) ' present(iq),present(ib)',present(iq),present(ib)
        write(out_unitp,*) ' One and only one ib or iq must be present'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF

      IF (present(iq) .AND. .NOT. allocated(psi%CvecG)) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) '  psi%CvecG is not allocated !!'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF

      IF (present(ib) .AND. .NOT. allocated(psi%CvecB)) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) '  psi%CvecB is not allocated !!'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF


      IF (any(shape(CVec_AT_iq) /= [psi%nb_bi*psi%nb_be])) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) '  Inconsitent shape of CVec_AT_iq'
        write(out_unitp,*) '  shape(CVec_AT_iq) ',shape(CVec_AT_iq)
        write(out_unitp,*) '  psi%nb_bi*psi%nb_be',psi%nb_bi*psi%nb_be
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF

    IF (psi%BasisnD%SparseGrid_type == 4 .AND. present(iq)) THEN
      ! not with ib since psi is not in Smolyak rep, we don't need that
     IF (present(OldPara)) THEN
       CALL get_iqSG_iSG_FROM_iq(iSG,iqSG,iq,psi%BasisnD%para_SGType2,OldPara,err_sub)
     ELSE
       CALL get_iqSG_iSG_FROM_iq(iSG,iqSG,iq,psi%BasisnD%para_SGType2,err_sub=err_sub)
     END if
     IF (err_sub /= 0) STOP 'Error in get_iqSG_iSG_FROM_iq'
     nqSG = psi%BasisnD%para_SGType2%tab_nq_OF_SRep(iSG)
     nb0  = psi%BasisnD%para_SGType2%nb0
     IF (iSG == 1) THEN
        i_qaie = iqSG
     ELSE
        i_qaie = iqSG + psi%BasisnD%para_SGType2%tab_Sum_nq_OF_SRep(iSG-1) * &
                                                psi%BasisnD%para_SGType2%nb0
     END IF

     psi%CvecG(i_qaie:i_qaie-1+nqSG*nb0:nqSG) = CVec_AT_iq

    ELSE
      IF (present(iq)) THEN

        DO i_be=1,psi%nb_be
        DO i_bi=1,psi%nb_bi
          i_bie  = (i_bi-1)+ (i_be-1)*psi%nb_bi

          i_qaie = iq + i_bie * psi%nb_qa

          psi%CvecG(i_qaie) = CVec_AT_iq(i_bie+1)

        END DO
        END DO
      ELSE ! ib is present
        DO i_be=1,psi%nb_be
        DO i_bi=1,psi%nb_bi
          i_bie  = (i_bi-1)+ (i_be-1)*psi%nb_bi
          i_baie = ib + i_bie * psi%nb_ba

          psi%CvecB(i_baie) = CVec_AT_iq(i_bie+1)

        END DO
        END DO
      END IF
    END IF

      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
        flush(out_unitp)
      END IF

    END SUBROUTINE set_CVec_OF_psi_AT_ind_a

  SUBROUTINE get_RMat_OF_psi_AT_ind_a(RMat_AT_iq,psi,iq,ib,OldPara)

!----- variables for the WP propagation ----------------------------
      TYPE (param_psi),intent(in)              :: psi
      real(kind=Rkind),intent(inout)           :: RMat_AT_iq(:,:)
      integer,         intent(in),    optional :: iq,ib
      TYPE (OldParam), intent(inout), optional :: OldPara

      integer :: i_qaie,i_baie,i_bi,i_be,iSG,iqSG,nqSG,nb0
      integer :: err_sub
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='get_RMat_OF_psi_AT_ind_a'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) write(out_unitp,*) 'BEGINNING ',name_sub

      IF (      present(iq) .AND.       present(ib) .OR.                       &
          .NOT. present(iq) .AND. .NOT. present(ib))     THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) ' present(iq),present(ib)',present(iq),present(ib)
        write(out_unitp,*) ' One and only one ib or iq must be present'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF

      IF (present(iq) .AND. .NOT. allocated(psi%RvecG)) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) '  psi%RvecG is not allocated !!'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF

      IF (present(ib) .AND. .NOT. allocated(psi%RvecB)) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) '  psi%RvecB is not allocated !!'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF


      IF (any(shape(RMat_AT_iq) /= [psi%nb_bi,psi%nb_be])) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) '  Inconsitent shape of RMat_AT_iq'
        write(out_unitp,*) '  shape(RMat_AT_iq) ',shape(RMat_AT_iq)
        write(out_unitp,*) '  psi%nb_bi,psi%nb_be',psi%nb_bi,psi%nb_be
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF

    RMat_AT_iq(:,:) = ZERO

    IF (psi%BasisnD%SparseGrid_type == 4 .AND. present(iq)) THEN
      IF (present(OldPara)) THEN
        CALL get_iqSG_iSG_FROM_iq(iSG,iqSG,iq,psi%BasisnD%para_SGType2,OldPara,err_sub)
      ELSE
        CALL get_iqSG_iSG_FROM_iq(iSG,iqSG,iq,psi%BasisnD%para_SGType2,err_sub=err_sub)
      END if
      IF (err_sub /= 0) STOP 'Error in get_iqSG_iSG_FROM_iq'
      nqSG = psi%BasisnD%para_SGType2%tab_nq_OF_SRep(iSG)
      nb0  = psi%BasisnD%para_SGType2%nb0
      IF (iSG == 1) THEN
         i_qaie = iqSG
      ELSE
         i_qaie = iqSG + psi%BasisnD%para_SGType2%tab_Sum_nq_OF_SRep(iSG-1) * &
                                                 psi%BasisnD%para_SGType2%nb0
      END IF
      RMat_AT_iq(:,:) = reshape(psi%RvecG(i_qaie:i_qaie-1+nqSG*nb0:nqSG),     &
                                shape=[psi%nb_bi,psi%nb_be])

    ELSE

      IF (present(iq)) THEN

        DO i_be=1,psi%nb_be
        DO i_bi=1,psi%nb_bi

          i_qaie = iq + ( (i_bi-1)+ (i_be-1)*psi%nb_bi ) * psi%nb_qa

          RMat_AT_iq(i_bi,i_be) = psi%RvecG(i_qaie)

        END DO
        END DO
      ELSE ! ib is present
        DO i_be=1,psi%nb_be
        DO i_bi=1,psi%nb_bi

          i_baie = ib + ( (i_bi-1)+ (i_be-1)*psi%nb_bi ) * psi%nb_ba

          RMat_AT_iq(i_bi,i_be) = psi%RvecB(i_baie)

        END DO
        END DO
      END IF
    END IF

      IF (debug) THEN
        write(out_unitp,*) 'RMat_AT_iq ',RMat_AT_iq
        write(out_unitp,*) 'END ',name_sub
        flush(out_unitp)
      END IF

  END SUBROUTINE get_RMat_OF_psi_AT_ind_a
  SUBROUTINE get_CMat_OF_psi_AT_ind_a(CMat_AT_iq,psi,iq,ib,OldPara)

!----- variables for the WP propagation ----------------------------
      TYPE (param_psi),   intent(in)             :: psi
      complex(kind=Rkind),intent(inout)          :: CMat_AT_iq(:,:)
      integer,            intent(in), optional   :: iq,ib
      TYPE (OldParam), intent(inout), optional   :: OldPara

      integer :: i_qaie,i_baie,i_bi,i_be,iSG,iqSG,nqSG,nb0
      integer :: err_sub
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='get_CMat_OF_psi_AT_ind_a'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) write(out_unitp,*) 'BEGINNING ',name_sub

      IF (      present(iq) .AND.       present(ib) .OR.                       &
          .NOT. present(iq) .AND. .NOT. present(ib))     THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) ' present(iq),present(ib)',present(iq),present(ib)
        write(out_unitp,*) ' One and only one ib or iq must be present'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF

      IF (present(iq) .AND. .NOT. allocated(psi%CvecG)) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) '  psi%CvecG is not allocated !!'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF

      IF (present(ib) .AND. .NOT. allocated(psi%CvecB)) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) '  psi%CvecB is not allocated !!'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF


      IF (any(shape(CMat_AT_iq) /= [psi%nb_bi,psi%nb_be])) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) '  Inconsitent shape of CMat_AT_iq'
        write(out_unitp,*) '  shape(CMat_AT_iq) ',shape(CMat_AT_iq)
        write(out_unitp,*) '  psi%nb_bi,psi%nb_be',psi%nb_bi,psi%nb_be
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF

    CMat_AT_iq(:,:) = ZERO
    IF (psi%BasisnD%SparseGrid_type == 4 .AND. present(iq)) THEN
      ! not with ib since psi is not in Smolyak rep, we don't need that
     IF (present(OldPara)) THEN
       CALL get_iqSG_iSG_FROM_iq(iSG,iqSG,iq,psi%BasisnD%para_SGType2,OldPara,err_sub)
     ELSE
       CALL get_iqSG_iSG_FROM_iq(iSG,iqSG,iq,psi%BasisnD%para_SGType2,err_sub=err_sub)
     END if
     IF (err_sub /= 0) STOP 'Error in get_iqSG_iSG_FROM_iq'
     nqSG = psi%BasisnD%para_SGType2%tab_nq_OF_SRep(iSG)
     nb0  = psi%BasisnD%para_SGType2%nb0
     IF (iSG == 1) THEN
        i_qaie = iqSG
     ELSE
        i_qaie = iqSG + psi%BasisnD%para_SGType2%tab_Sum_nq_OF_SRep(iSG-1) * &
                                                psi%BasisnD%para_SGType2%nb0
     END IF

     ! DO i_be=1,psi%nb_be
     ! DO i_bi=1,psi%nb_bi
     !
     !   CMat_AT_iq(i_bi,i_be) = psi%CvecG(i_qaie)
     !   i_qaie = i_qaie + nqSG
     !
     ! END DO
     ! END DO
     CMat_AT_iq(:,:) = reshape(psi%CvecG(i_qaie:i_qaie-1+nqSG*nb0:nqSG),     &
                               shape=[psi%nb_bi,psi%nb_be])

    ELSE
      IF (present(iq)) THEN

        DO i_be=1,psi%nb_be
        DO i_bi=1,psi%nb_bi

          i_qaie = iq + ( (i_bi-1)+ (i_be-1)*psi%nb_bi ) * psi%nb_qa

          CMat_AT_iq(i_bi,i_be) = psi%CvecG(i_qaie)

        END DO
        END DO
      ELSE ! ib is present
        DO i_be=1,psi%nb_be
        DO i_bi=1,psi%nb_bi

          i_baie = ib + ( (i_bi-1)+ (i_be-1)*psi%nb_bi ) * psi%nb_ba

          CMat_AT_iq(i_bi,i_be) = psi%CvecB(i_baie)

        END DO
        END DO
      END IF
    END IF

    IF (debug) THEN
      write(out_unitp,*) 'CMat_AT_iq ',CMat_AT_iq
      write(out_unitp,*) 'END ',name_sub
      flush(out_unitp)
    END IF

  END SUBROUTINE get_CMat_OF_psi_AT_ind_a

  SUBROUTINE set_RMat_OF_psi_AT_ind_a(RMat_AT_iq,psi,iq,ib,OldPara)

!----- variables for the WP propagation ----------------------------
    TYPE (param_psi),intent(inout)           :: psi
    real(kind=Rkind),intent(in)              :: RMat_AT_iq(:,:)
    integer,         intent(in), optional    :: iq,ib
    TYPE (OldParam), intent(inout), optional :: OldPara

    integer :: i_qaie,i_baie,i_bi,i_be,iSG,iqSG,nqSG,nb0
    INTEGER :: err_sub
!---- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='set_RMat_OF_psi_AT_ind_a'
    logical, parameter :: debug=.FALSE.
    !logical, parameter :: debug=.TRUE.
!----------------------------------------------------------
    IF (debug) write(out_unitp,*) 'BEGINNING ',name_sub
    IF (debug) write(out_unitp,*) 'RMat_AT_iq ',RMat_AT_iq


    IF (      present(iq) .AND.       present(ib) .OR.                       &
        .NOT. present(iq) .AND. .NOT. present(ib))     THEN
      write(out_unitp,*) ' ERROR ',name_sub
      write(out_unitp,*) ' present(iq),present(ib)',present(iq),present(ib)
      write(out_unitp,*) ' One and only one ib or iq must be present'
      write(out_unitp,*) ' CHECK the fortran !!'
      STOP
    END IF

      IF (present(iq) .AND. .NOT. allocated(psi%RvecG)) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) '  psi%RvecG is not allocated !!'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF

      IF (present(ib) .AND. .NOT. allocated(psi%RvecB)) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) '  psi%RvecB is not allocated !!'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF


    IF (any(shape(RMat_AT_iq) /= [psi%nb_bi,psi%nb_be])) THEN
      write(out_unitp,*) ' ERROR ',name_sub
      write(out_unitp,*) '  Inconsitent shape of RMat_AT_iq'
      write(out_unitp,*) '  shape(RMat_AT_iq) ',shape(RMat_AT_iq)
      write(out_unitp,*) '  psi%nb_bi,psi%nb_be',psi%nb_bi,psi%nb_be
      write(out_unitp,*) ' CHECK the fortran !!'
      STOP
    END IF

    IF (psi%BasisnD%SparseGrid_type == 4 .AND. present(iq)) THEN
      ! not with ib since psi is not in Smolyak rep, we don't need that
     IF (present(OldPara)) THEN
       CALL get_iqSG_iSG_FROM_iq(iSG,iqSG,iq,psi%BasisnD%para_SGType2,OldPara,err_sub)
     ELSE
       CALL get_iqSG_iSG_FROM_iq(iSG,iqSG,iq,psi%BasisnD%para_SGType2,err_sub=err_sub)
     END if
     IF (err_sub /= 0) STOP 'Error in get_iqSG_iSG_FROM_iq'
     nqSG = psi%BasisnD%para_SGType2%tab_nq_OF_SRep(iSG)
     nb0  = psi%BasisnD%para_SGType2%nb0
     IF (iSG == 1) THEN
        i_qaie = iqSG
     ELSE
        i_qaie = iqSG + psi%BasisnD%para_SGType2%tab_Sum_nq_OF_SRep(iSG-1) * &
                                                psi%BasisnD%para_SGType2%nb0
     END IF

     psi%RvecG(i_qaie:i_qaie-1+nqSG*nb0:nqSG) = reshape(RMat_AT_iq,shape=[nb0])

    ELSE
     IF (present(iq)) THEN

       DO i_be=1,psi%nb_be
       DO i_bi=1,psi%nb_bi

         i_qaie = iq + ( (i_bi-1)+ (i_be-1)*psi%nb_bi ) * psi%nb_qa

         psi%RvecG(i_qaie) = RMat_AT_iq(i_bi,i_be)

       END DO
       END DO
     ELSE ! ib is present
       DO i_be=1,psi%nb_be
       DO i_bi=1,psi%nb_bi

         i_baie = ib + ( (i_bi-1)+ (i_be-1)*psi%nb_bi ) * psi%nb_ba

         psi%RvecB(i_baie) = RMat_AT_iq(i_bi,i_be)

       END DO
       END DO
     END IF
    END IF

      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
        flush(out_unitp)
      END IF

  END SUBROUTINE set_RMat_OF_psi_AT_ind_a
  SUBROUTINE set_CMat_OF_psi_AT_ind_a(CMat_AT_iq,psi,iq,ib,OldPara)

!----- variables for the WP propagation ----------------------------
    TYPE (param_psi),   intent(inout)          :: psi
    complex(kind=Rkind),intent(in)             :: CMat_AT_iq(:,:)
    integer,            intent(in), optional   :: iq,ib
    TYPE (OldParam), intent(inout), optional   :: OldPara

    integer :: i_qaie,i_baie,i_bi,i_be,iSG,iqSG,nqSG,nb0
    INTEGER :: err_sub
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='set_CMat_OF_psi_AT_ind_a'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) write(out_unitp,*) 'BEGINNING ',name_sub
      IF (debug) write(out_unitp,*) 'CMat_AT_iq ',CMat_AT_iq

      IF (      present(iq) .AND.       present(ib) .OR.                       &
          .NOT. present(iq) .AND. .NOT. present(ib))     THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) ' present(iq),present(ib)',present(iq),present(ib)
        write(out_unitp,*) ' One and only one ib or iq must be present'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF

      IF (present(iq) .AND. .NOT. allocated(psi%CvecG)) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) '  psi%CvecG is not allocated !!'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF

      IF (present(ib) .AND. .NOT. allocated(psi%CvecB)) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) '  psi%CvecB is not allocated !!'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF


      IF (any(shape(CMat_AT_iq) /= [psi%nb_bi,psi%nb_be])) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) '  Inconsitent shape of CMat_AT_iq'
        write(out_unitp,*) '  shape(CMat_AT_iq) ',shape(CMat_AT_iq)
        write(out_unitp,*) '  psi%nb_bi,psi%nb_be',psi%nb_bi,psi%nb_be
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF

    IF (psi%BasisnD%SparseGrid_type == 4 .AND. present(iq)) THEN
      ! not with ib since psi is not in Smolyak rep, we don't need that
     IF (present(OldPara)) THEN
       CALL get_iqSG_iSG_FROM_iq(iSG,iqSG,iq,psi%BasisnD%para_SGType2,OldPara,err_sub)
     ELSE
       CALL get_iqSG_iSG_FROM_iq(iSG,iqSG,iq,psi%BasisnD%para_SGType2,err_sub=err_sub)
     END if
     IF (err_sub /= 0) STOP 'Error in get_iqSG_iSG_FROM_iq'
     nqSG = psi%BasisnD%para_SGType2%tab_nq_OF_SRep(iSG)
     nb0  = psi%BasisnD%para_SGType2%nb0
     IF (iSG == 1) THEN
        i_qaie = iqSG
     ELSE
        i_qaie = iqSG + psi%BasisnD%para_SGType2%tab_Sum_nq_OF_SRep(iSG-1) * &
                                                psi%BasisnD%para_SGType2%nb0
     END IF

     psi%CvecG(i_qaie:i_qaie-1+nqSG*nb0:nqSG) = reshape(CMat_AT_iq,shape=[nb0])

    ELSE
      IF (present(iq)) THEN

        DO i_be=1,psi%nb_be
        DO i_bi=1,psi%nb_bi

          i_qaie = iq + ( (i_bi-1)+ (i_be-1)*psi%nb_bi ) * psi%nb_qa

          psi%CvecG(i_qaie) = CMat_AT_iq(i_bi,i_be)

        END DO
        END DO
      ELSE ! ib is present
        DO i_be=1,psi%nb_be
        DO i_bi=1,psi%nb_bi

          i_baie = ib + ( (i_bi-1)+ (i_be-1)*psi%nb_bi ) * psi%nb_ba

          psi%CvecB(i_baie) = CMat_AT_iq(i_bi,i_be)

        END DO
        END DO
      END IF
    END IF

      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
        flush(out_unitp)
      END IF

    END SUBROUTINE set_CMat_OF_psi_AT_ind_a


!================================================================
!
!     psi = R
!     psi = C
!
!================================================================
      SUBROUTINE R_TOpsi(psi,R)

!----- variables for the WP propagation ----------------------------
      CLASS (param_psi), intent(inout) :: psi
      real (kind=Rkind), intent(in)    :: R

      integer :: i

!     write(out_unitp,*) 'BEGINNING R_TOpsi'

       CALL alloc_psi(psi)

      IF (psi%BasisRep) THEN
        IF (allocated(psi%RvecB)) THEN
          psi%RvecB(:) = R
        END IF
        IF (allocated(psi%CvecB)) THEN
          psi%CvecB(:) = cmplx(R,ZERO,kind=Rkind)
        END IF
      END IF

      IF (psi%GridRep) THEN
        IF (allocated(psi%RvecG)) THEN
          psi%RvecG(:) = R
        END IF
        IF (allocated(psi%CvecG)) THEN
          psi%CvecG(:) = cmplx(R,ZERO,kind=Rkind)
        END IF
      END IF

      IF(psi%SRG_MPI) THEN
        IF(allocated(psi%SR_G)) THEN
          psi%SR_G(:,1)=R
          IF(psi%cplx) psi%SR_G(:,2)=ZERO
        ENDIF
      ENDIF

      IF(psi%SRB_MPI) THEN
        IF(allocated(psi%SR_B)) THEN
          psi%SR_B(:,1)=R
          IF(psi%cplx) psi%SR_B(:,2)=ZERO
        ENDIF
      ENDIF

      psi%norm2 = ZERO

      IF (R == ZERO) THEN
        psi%symab = -2
      ELSE
        psi%symab = -1
      END IF

!     write(out_unitp,*) 'END R_TOpsi'

      END SUBROUTINE R_TOpsi

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      SUBROUTINE C_TOpsi(psi,C)

!----- variables for the WP propagation ----------------------------
      CLASS (param_psi),    intent(inout) :: psi
      complex (kind=Rkind), intent(in)    :: C

      integer :: i

!     write(out_unitp,*) 'BEGINNING C_TOpsi'

       CALL alloc_psi(psi)

      IF (psi%BasisRep) THEN
        IF (allocated(psi%RvecB)) THEN
          psi%RvecB(:) = real(C,kind=Rkind)
        END IF
        IF (allocated(psi%CvecB)) THEN
          psi%CvecB(:) = C
        END IF
      END IF

      IF (psi%GridRep) THEN
        IF (allocated(psi%RvecG)) THEN
          psi%RvecG(:) = real(C,kind=Rkind)
        END IF
        IF (allocated(psi%CvecG)) THEN
          psi%CvecG(:) = C
        END IF
      END IF

      IF(psi%SRG_MPI) THEN
        IF(allocated(psi%SR_G)) THEN
          psi%SR_G(:,1)=real(C,kind=Rkind)
          IF(psi%cplx) psi%SR_G(:,2)=aimag(C)
        ENDIF
      ENDIF

      IF(psi%SRB_MPI) THEN
        IF(allocated(psi%SR_B)) THEN
          psi%SR_B(:,1)=real(C,kind=Rkind)
          IF(psi%cplx) psi%SR_B(:,2)=aimag(C)
        ENDIF
      ENDIF

      psi%norm2 = ZERO


      IF (abs(C) == ZERO) THEN
        psi%symab = -2
      ELSE
        psi%symab = -1
      END IF

!     write(out_unitp,*) 'END C_TOpsi'

      END SUBROUTINE C_TOpsi

      SUBROUTINE Set_Random_psi(psi,option)

!----- variables for the WP propagation ----------------------------
      TYPE (param_psi),intent(inout)          :: psi
      integer,         intent(in),   optional :: option

      integer          :: i
      real(kind=Rkind) :: a,b


!     write(out_unitp,*) 'BEGINNING Set_Random_psi'

       CALL alloc_psi(psi)

      IF (psi%BasisRep) THEN
        IF (allocated(psi%RvecB)) THEN
          CALL random_number(psi%RvecB(:))
        END IF
        IF (allocated(psi%CvecB)) THEN
          DO i=1,size(psi%CvecB)
            CALL random_number(a)
            CALL random_number(b)
            psi%CvecB(i) = cmplx(a,b,kind=Rkind)
          END DO
        END IF
      END IF

      IF (psi%GridRep) THEN
        IF (allocated(psi%RvecG)) THEN
          CALL random_number(psi%RvecG(:))
        END IF
        IF (allocated(psi%CvecG)) THEN
          DO i=1,size(psi%CvecG)
            CALL random_number(a)
            CALL random_number(b)
            psi%CvecG(i) = cmplx(a,b,kind=Rkind)
          END DO
        END IF
      END IF

      psi%norm2 = ZERO
      psi%symab = -1

!     write(out_unitp,*) 'END Set_Random_psi'

      END SUBROUTINE Set_Random_psi

      SUBROUTINE Set_psi_With_index_R(psi,R,ind_a,ind_i,ind_e,ind_aie)

!----- variables for the WP propagation ----------------------------
      TYPE (param_psi),intent(inout)          :: psi
      real(kind=Rkind),intent(in)             :: R
      integer,         intent(in),   optional :: ind_a,ind_i,ind_e,ind_aie

      integer          :: ind_a_loc,ind_i_loc,ind_e_loc,ind_aie_loc

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Set_psi_With_index_R'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) write(out_unitp,*) 'BEGINNING ',name_sub

      IF (present(ind_aie) .AND. present(ind_a)) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) ' Both ind_aie and ind_a are present !!'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF
      IF (.NOT. present(ind_aie) .AND. .NOT. present(ind_a)) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) ' Both ind_aie and ind_a are absent !!'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF

      IF (.NOT. allocated(psi%RvecB) .AND. .NOT. allocated(psi%CvecB)) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) ' Both psi%RvecB and psi%CvecB are not allocated !!'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF

      IF (present(ind_a)) THEN
        ind_a_loc = ind_a

        IF (present(ind_i)) THEN
          ind_i_loc = ind_i-1
        ELSE
          ind_i_loc = 0
        END IF

        IF (present(ind_e)) THEN
          ind_e_loc = ind_e-1
        ELSE
          ind_e_loc = 0
        END IF

        ind_aie_loc = ind_a_loc + (ind_i_loc+ind_e_loc*psi%nb_bi) * psi%nb_ba

      ELSE ! it means ind_aie is present
        ind_aie_loc = ind_aie
      END IF

      IF (ind_aie_loc < 1 .OR. ind_aie_loc > psi%nb_tot) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) ' ind_aie_loc is out of range. ind_aie_loc:',ind_aie_loc
        write(out_unitp,*) '    ind_aie_loc:',ind_aie_loc
        write(out_unitp,*) '    range: [1   ...',psi%nb_tot,']'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF


      IF (allocated(psi%RvecB)) psi%RvecB(ind_aie_loc) = R
      IF (allocated(psi%CvecB)) psi%CvecB(ind_aie_loc) = R


      psi%norm2 = ZERO
      psi%symab = -1

      IF (debug) THEN
        write(out_unitp,*) 'ind_aie_loc ',ind_aie_loc
        write(out_unitp,*) 'END ',name_sub
        flush(out_unitp)
      END IF

      END SUBROUTINE Set_psi_With_index_R

      SUBROUTINE Set_psi_With_index_C(psi,C,ind_a,ind_i,ind_e,ind_aie)

!----- variables for the WP propagation ----------------------------
      TYPE (param_psi),intent(inout)          :: psi
      complex(kind=Rkind),intent(in)          :: C
      integer,         intent(in),   optional :: ind_a,ind_i,ind_e,ind_aie

      integer          :: ind_a_loc,ind_i_loc,ind_e_loc,ind_aie_loc

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Set_psi_With_index_C'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) write(out_unitp,*) 'BEGINNING ',name_sub

      IF (present(ind_aie) .AND. present(ind_a)) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) ' Both ind_aie and ind_a are present !!'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF
      IF (.NOT. present(ind_aie) .AND. .NOT. present(ind_a)) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) ' Both ind_aie and ind_a are absent !!'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF

      IF (.NOT. allocated(psi%RvecB) .AND. .NOT. allocated(psi%CvecB)) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) ' Both psi%RvecB and psi%CvecB are not allocated !!'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF

      IF (present(ind_a)) THEN
        ind_a_loc = ind_a

        IF (present(ind_i)) THEN
          ind_i_loc = ind_i-1
        ELSE
          ind_i_loc = 0
        END IF

        IF (present(ind_e)) THEN
          ind_e_loc = ind_e-1
        ELSE
          ind_e_loc = 0
        END IF

        ind_aie_loc = ind_a_loc + (ind_i_loc+ind_e_loc*psi%nb_bi) * psi%nb_ba

      ELSE ! it means ind_aie is present
        ind_aie_loc = ind_aie
      END IF


      IF (ind_aie_loc < 1 .OR. ind_aie_loc > psi%nb_tot) THEN
        write(out_unitp,*) ' ERROR ',name_sub
        write(out_unitp,*) ' ind_aie_loc is out of range. ind_aie_loc:',ind_aie_loc
        write(out_unitp,*) '    ind_aie_loc:',ind_aie_loc
        write(out_unitp,*) '    range: [1   ...',psi%nb_tot,']'
        write(out_unitp,*) ' CHECK the fortran !!'
        STOP
      END IF


      IF (allocated(psi%RvecB)) psi%RvecB(ind_aie_loc) = C
      IF (allocated(psi%CvecB)) psi%CvecB(ind_aie_loc) = C


      psi%norm2 = ZERO
      psi%symab = -1

      IF (debug) write(out_unitp,*) 'END ',name_sub

      END SUBROUTINE Set_psi_With_index_C
!================================================================
!
!     psi+psi, psi+R, R+psi....
!
!================================================================
      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
          FUNCTION psi1_plus_psi2(psi1,psi2)
            TYPE (param_psi)              :: psi1_plus_psi2
            CLASS (param_psi), intent(in) :: psi1,psi2

            integer           :: err,i

            !write(out_unitp,*) 'BEGINNING psi1_plus_psi2'


            !- define and allocate psi1_plus_psi2 ----
            IF (.NOT. psi1_plus_psi2%init) THEN
              !write(out_unitp,*) 'psi1_plus_psi2 is not initialized yet'
              CALL copy_psi2TOpsi1(psi1_plus_psi2,psi1)
            ELSE
              !write(out_unitp,*) 'psi1_plus_psi2 is initialized already'
            END IF


            IF (psi1%GridRep .AND. psi2%GridRep) THEN
              IF (psi1%cplx .AND. psi2%cplx) THEN
                psi1_plus_psi2%CvecG = psi1%CvecG + psi2%CvecG
              ELSE IF (.NOT. psi1%cplx .AND. .NOT. psi2%cplx) THEN
                psi1_plus_psi2%RvecG = psi1%RvecG + psi2%RvecG
              ELSE
                write(out_unitp,*) ' ERROR : I CANNOT mix real and complex psi !!'
                write(out_unitp,*) ' psi1%cplx,psi2%cplx',psi1%cplx,psi2%cplx
                STOP
              END IF
            END IF

            IF (psi1%BasisRep .AND. psi2%BasisRep) THEN
              IF (psi1%cplx .AND. psi2%cplx) THEN
                psi1_plus_psi2%CvecB = psi1%CvecB + psi2%CvecB
              ELSE IF (.NOT. psi1%cplx .AND. .NOT. psi2%cplx) THEN
                psi1_plus_psi2%RvecB = psi1%RvecB + psi2%RvecB
              ELSE
                write(out_unitp,*) ' ERROR : I CANNOT mix real and complex psi !!'
                write(out_unitp,*) ' psi1%cplx,psi2%cplx',psi1%cplx,psi2%cplx
                STOP
              END IF
            END IF

            IF (psi1%symab == psi2%symab) THEN
               psi1_plus_psi2%symab = psi1%symab
            ELSE IF (psi1%symab == -2) THEN
               psi1_plus_psi2%symab = psi2%symab
            ELSE IF (psi2%symab == -2) THEN
               psi1_plus_psi2%symab = psi1%symab
            ELSE
              psi1_plus_psi2%symab = -1
            END IF

            !write(out_unitp,*) 'END psi1_plus_psi2'

          END FUNCTION psi1_plus_psi2
          SUBROUTINE psi1_pluseq_psi2(psi1,psi2)
            CLASS (param_psi), intent(inout)   :: psi1
            TYPE  (param_psi), intent(in)      :: psi2

            integer           :: err,i

            !write(out_unitp,*) 'BEGINNING psi1_pluseq_psi2'


            !- define and allocate psi1_plus_psi2 ----
            IF (.NOT. psi1%init) THEN
              write(out_unitp,*) 'psi1 is not initialized'
              STOP 'ERROR in psi1_pluseq_psi2, psi1 is not initalized'
            ELSE
              write(out_unitp,*) 'psi1 is initialized'
            END IF


            IF (psi1%GridRep .AND. psi2%GridRep) THEN
              IF (psi1%cplx .AND. psi2%cplx) THEN
                psi1%CvecG = psi1%CvecG + psi2%CvecG
              ELSE IF (.NOT. psi1%cplx .AND. .NOT. psi2%cplx) THEN
                psi1%RvecG = psi1%RvecG + psi2%RvecG
              ELSE
                write(out_unitp,*) ' ERROR : I CANNOT mix real and complex psi !!'
                write(out_unitp,*) ' psi1%cplx,psi2%cplx',psi1%cplx,psi2%cplx
                STOP
              END IF
            END IF

            IF (psi1%BasisRep .AND. psi2%BasisRep) THEN
              IF (psi1%cplx .AND. psi2%cplx) THEN
                psi1%CvecB = psi1%CvecB + psi2%CvecB
              ELSE IF (.NOT. psi1%cplx .AND. .NOT. psi2%cplx) THEN
                psi1%RvecB = psi1%RvecB + psi2%RvecB
              ELSE
                write(out_unitp,*) ' ERROR : I CANNOT mix real and complex psi !!'
                write(out_unitp,*) ' psi1%cplx,psi2%cplx',psi1%cplx,psi2%cplx
                STOP
              END IF
            END IF

            IF (psi1%symab == psi2%symab) THEN
               CONTINUE ! nothing
            ELSE IF (psi1%symab == -2) THEN
               psi1%symab = psi2%symab
            ELSE IF (psi2%symab == -2) THEN
               CONTINUE ! nothing
            ELSE
              psi1%symab = -1
            END IF

            !write(out_unitp,*) 'END psi1_pluseq_psi2'

          END SUBROUTINE psi1_pluseq_psi2

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
          FUNCTION R_plus_psi(R,psi)
            CLASS (param_psi), intent (in) :: psi
            real (kind=Rkind),    intent (in) :: R
            TYPE (param_psi) :: R_plus_psi
            integer           :: err,i


!           write(out_unitp,*) 'BEGINNING R_plus_psi'

!           - define and allocate R_plus_psi ----
            CALL copy_psi2TOpsi1(R_plus_psi,psi)
!           -----------------------------------------


            IF (.NOT. psi%BasisRep) THEN
              write(out_unitp,*) ' ERROR in R_plus_psi: I CANNOT use psi in GridRep !!'
              write(out_unitp,*) ' psi%BasisRep',psi%BasisRep
              STOP
            END IF


            IF (psi%cplx) THEN
              R_plus_psi%CvecB = psi%CvecB + cmplx(R,kind=Rkind)
            ELSE
              R_plus_psi%RvecB = psi%RvecB + R
            END IF

            IF (R == ZERO) THEN
              R_plus_psi%symab = psi%symab
            ELSE
              R_plus_psi%symab = -1
            END IF

!           write(out_unitp,*) 'END R_plus_psi'

          END FUNCTION R_plus_psi
      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
          FUNCTION psi_plus_R(psi,R)
            CLASS (param_psi), intent (in) :: psi
            real (kind=Rkind),    intent (in) :: R
            TYPE (param_psi)  :: psi_plus_R
            integer           :: err,i

!           write(out_unitp,*) 'BEGINNING psi_plus_R'

!           - define and allocate psi_plus_R ----
            CALL copy_psi2TOpsi1(psi_plus_R,psi)
!           -----------------------------------------


            IF (.NOT. psi%BasisRep) THEN
              write(out_unitp,*) ' ERROR in psi_plus_R: I CANNOT use psi in GridRep !!'
              write(out_unitp,*) ' psi%BasisRep',psi%BasisRep
              STOP
            END IF


            IF (psi%cplx) THEN
              psi_plus_R%CvecB = psi%CvecB + cmplx(R,kind=Rkind)
            ELSE
              psi_plus_R%RvecB = psi%RvecB + R
            END IF

            IF (R == ZERO) THEN
              psi_plus_R%symab = psi%symab
            ELSE
              psi_plus_R%symab = -1
            END IF

!           write(out_unitp,*) 'END psi_plus_R'

          END FUNCTION psi_plus_R

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
          FUNCTION C_plus_psi(C,psi)
            CLASS (param_psi), intent (in) :: psi
            complex (kind=Rkind), intent (in) :: C
            TYPE (param_psi)  :: C_plus_psi
            integer           :: err,i

!           write(out_unitp,*) 'BEGINNING C_plus_psi'

!           - define and allocate C_plus_psi ----
            CALL copy_psi2TOpsi1(C_plus_psi,psi)
!           -----------------------------------------


            IF (.NOT. psi%BasisRep) THEN
              write(out_unitp,*) ' ERROR C_plus_psi: I CANNOT use psi in GridRep !!'
              write(out_unitp,*) ' psi%BasisRep',psi%BasisRep
              STOP
            END IF


            IF (psi%cplx) THEN
              C_plus_psi%CvecB = psi%CvecB + C
            ELSE
              write(out_unitp,*) ' ERROR : in C_plus_psi'
              write(out_unitp,*) ' I cannot sum a real psi and a complex'
              write(out_unitp,*) 'psi%cplx,C',psi%cplx,C
              STOP
            END IF

            IF (abs(C) == ZERO) THEN
              C_plus_psi%symab = psi%symab
            ELSE
              C_plus_psi%symab = -1
            END IF

!           write(out_unitp,*) 'END C_plus_psi'

          END FUNCTION C_plus_psi

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
          FUNCTION psi_plus_C(psi,C)
            CLASS (param_psi), intent (in) :: psi
            complex (kind=Rkind), intent (in) :: C
            TYPE (param_psi)  :: psi_plus_C
            integer           :: err,i

!           write(out_unitp,*) 'BEGINNING psi_plus_C'

!           - define and allocate psi_plus_C ----
            CALL copy_psi2TOpsi1(psi_plus_C,psi)
!           -----------------------------------------


            IF (.NOT. psi%BasisRep) THEN
              write(out_unitp,*) ' ERROR psi_plus_C: I CANNOT use psi in GridRep !!'
              write(out_unitp,*) ' psi%BasisRep',psi%BasisRep
              STOP
            END IF


            IF (psi%cplx) THEN
              psi_plus_C%CvecB = psi%CvecB + C
            ELSE
              write(out_unitp,*) ' ERROR : in psi_plus_C'
              write(out_unitp,*) ' I cannot sum a real psi and a complex'
              write(out_unitp,*) 'psi%cplx,C',psi%cplx,C
              STOP
            END IF

            IF (abs(C) == ZERO) THEN
              psi_plus_C%symab = psi%symab
            ELSE
              psi_plus_C%symab = -1
            END IF

!           write(out_unitp,*) 'END psi_plus_C'

          END FUNCTION psi_plus_C
          FUNCTION plus_psi(psi)
            CLASS (param_psi), intent (in) :: psi
            TYPE (param_psi) :: plus_psi


!           write(out_unitp,*) 'BEGINNING plus_psi'

            CALL copy_psi2TOpsi1(plus_psi,psi)

!           write(out_unitp,*) 'END plus_psi'

          END FUNCTION plus_psi
!================================================================
!
!     psi-psi, psi-R, R-psi....
!
!================================================================
      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
          FUNCTION psi1_minus_psi2(psi1,psi2)
            CLASS (param_psi), intent(in) :: psi1,psi2
            TYPE (param_psi) :: psi1_minus_psi2
            integer          :: err,i


!           write(out_unitp,*) 'BEGINNING psi1_minus_psi2'

!           - define and allocate psi1_minus_psi2 ----
            CALL copy_psi2TOpsi1(psi1_minus_psi2,psi1)
!           -----------------------------------------

            IF (psi1%GridRep .AND. psi2%GridRep) THEN
              IF (psi1%cplx .AND. psi2%cplx) THEN
                psi1_minus_psi2%CvecG = psi1%CvecG - psi2%CvecG
              ELSE IF (.NOT. psi1%cplx .AND. .NOT. psi2%cplx) THEN
                psi1_minus_psi2%RvecG = psi1%RvecG - psi2%RvecG
              ELSE
                write(out_unitp,*) ' ERROR : I CANNOT mix real and complex psi !!'
                write(out_unitp,*) ' psi1%cplx,psi2%cplx',psi1%cplx,psi2%cplx
                STOP
              END IF
            END IF

            IF (psi1%BasisRep .AND. psi2%BasisRep) THEN
              IF (psi1%cplx .AND. psi2%cplx) THEN
                psi1_minus_psi2%CvecB = psi1%CvecB - psi2%CvecB
              ELSE IF (.NOT. psi1%cplx .AND. .NOT. psi2%cplx) THEN
                psi1_minus_psi2%RvecB = psi1%RvecB - psi2%RvecB
              ELSE
                write(out_unitp,*) ' ERROR : I CANNOT mix real and complex psi !!'
                write(out_unitp,*) ' psi1%cplx,psi2%cplx',psi1%cplx,psi2%cplx
                STOP
              END IF
            END IF

            IF (psi1%symab == psi2%symab) THEN
               psi1_minus_psi2%symab = psi1%symab
            ELSE IF (psi1%symab == -2) THEN
               psi1_minus_psi2%symab = psi2%symab
            ELSE IF (psi2%symab == -2) THEN
               psi1_minus_psi2%symab = psi1%symab
            ELSE
              psi1_minus_psi2%symab = -1
            END IF

!           write(out_unitp,*) 'END psi1_minus_psi2'

          END FUNCTION psi1_minus_psi2

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
          FUNCTION R_minus_psi(R,psi)
            CLASS (param_psi),    intent (in) :: psi
            real (kind=Rkind),    intent (in) :: R
            TYPE (param_psi) :: R_minus_psi


!           write(out_unitp,*) 'BEGINNING R_minus_psi'

!           - define and allocate R_minus_psi ----
            CALL copy_psi2TOpsi1(R_minus_psi,psi)
!           -----------------------------------------


            IF (.NOT. psi%BasisRep) THEN
              write(out_unitp,*) ' ERROR R_minus_psi: I CANNOT use psi in GridRep !!'
              write(out_unitp,*) ' psi%BasisRep',psi%BasisRep
              STOP
            END IF


            IF (psi%cplx) THEN
              R_minus_psi%CvecB = cmplx(R,kind=Rkind) - psi%CvecB
            ELSE
              R_minus_psi%RvecB = R - psi%RvecB
            END IF

            IF (R == ZERO) THEN
              R_minus_psi%symab = psi%symab
            ELSE
              R_minus_psi%symab = -1
            END IF


!           write(out_unitp,*) 'END R_minus_psi'

          END FUNCTION R_minus_psi

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
          FUNCTION psi_minus_R(psi,R)
            CLASS (param_psi), intent (in) :: psi
            real (kind=Rkind),    intent (in) :: R
            TYPE (param_psi)  :: psi_minus_R
            integer           :: err,i

!           write(out_unitp,*) 'BEGINNING psi_minus_R'

!           - define and allocate psi_minus_R ----
            CALL copy_psi2TOpsi1(psi_minus_R,psi)
!           -----------------------------------------


            IF (.NOT. psi%BasisRep) THEN
              write(out_unitp,*) ' ERROR psi_minus_R: I CANNOT use psi in GridRep !!'
              write(out_unitp,*) ' psi%BasisRep',psi%BasisRep
              STOP
            END IF


            IF (psi%cplx) THEN
              psi_minus_R%CvecB = psi%CvecB - cmplx(R,kind=Rkind)
            ELSE
              psi_minus_R%RvecB = psi%RvecB - R
            END IF

            IF (R == ZERO) THEN
              psi_minus_R%symab = psi%symab
            ELSE
              psi_minus_R%symab = -1
            END IF

!           write(out_unitp,*) 'END psi_minus_R'

          END FUNCTION psi_minus_R

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
          FUNCTION C_minus_psi(C,psi)
            CLASS (param_psi), intent (in) :: psi
            complex (kind=Rkind), intent (in) :: C
            TYPE (param_psi)  :: C_minus_psi
            integer           :: err,i

!           write(out_unitp,*) 'BEGINNING C_minus_psi'

!           - define and allocate C_minus_psi ----
            CALL copy_psi2TOpsi1(C_minus_psi,psi)
!           -----------------------------------------


            IF (.NOT. psi%BasisRep) THEN
              write(out_unitp,*) ' ERROR C_minus_psi: I CANNOT use psi in GridRep !!'
              write(out_unitp,*) ' psi%BasisRep',psi%BasisRep
              STOP
            END IF


            IF (psi%cplx) THEN
              C_minus_psi%CvecB = C - psi%CvecB
            ELSE
              write(out_unitp,*) ' ERROR : in C_minus_psi'
              write(out_unitp,*) ' I cannot sum a real psi and a complex'
              write(out_unitp,*) 'psi%cplx,C',psi%cplx,C
              STOP
            END IF

            IF (abs(C) == ZERO) THEN
              C_minus_psi%symab = psi%symab
            ELSE
              C_minus_psi%symab = -1
            END IF

!           write(out_unitp,*) 'END C_minus_psi'

          END FUNCTION C_minus_psi

          !!@description: TODO
          FUNCTION psi_minus_C(psi,C)
            CLASS (param_psi), intent (in) :: psi
            complex (kind=Rkind), intent (in) :: C
            TYPE (param_psi)  :: psi_minus_C
            integer           :: err,i

!           write(out_unitp,*) 'BEGINNING psi_minus_C'

!           - define and allocate psi_minus_C ----
            CALL copy_psi2TOpsi1(psi_minus_C,psi)
!           -----------------------------------------


            IF (.NOT. psi%BasisRep) THEN
              write(out_unitp,*) ' ERROR psi_minus_C: I CANNOT use psi in GridRep !!'
              write(out_unitp,*) ' psi%BasisRep',psi%BasisRep
              STOP
            END IF


            IF (psi%cplx) THEN
              psi_minus_C%CvecB = psi%CvecB - C
            ELSE
              write(out_unitp,*) ' ERROR : in psi_minus_C'
              write(out_unitp,*) ' I cannot sum a real psi and a complex'
              write(out_unitp,*) 'psi%cplx,C',psi%cplx,C
              STOP
            END IF

            IF (abs(C) == ZERO) THEN
              psi_minus_C%symab = psi%symab
            ELSE
              psi_minus_C%symab = -1
            END IF

!           write(out_unitp,*) 'END psi_minus_C'

          END FUNCTION psi_minus_C

          FUNCTION minus_psi(psi)
            CLASS (param_psi),    intent (in) :: psi
            TYPE (param_psi) :: minus_psi


!           write(out_unitp,*) 'BEGINNING minus_psi'

!           - define and allocate minus_psi ----
            CALL copy_psi2TOpsi1(minus_psi,psi)
!           -----------------------------------------

          IF (psi%BasisRep) THEN
            IF (psi%cplx) THEN
              minus_psi%CvecB(:) =  - psi%CvecB
            ELSE
              minus_psi%RvecB(:) =  - psi%RvecB
            END IF
          END IF

          IF (psi%GridRep) THEN
            IF (psi%cplx) THEN
              minus_psi%CvecG(:) =  - psi%CvecG
            ELSE
              minus_psi%RvecG(:) =  - psi%RvecG
            END IF
          END IF

!           write(out_unitp,*) 'END minus_psi'

          END FUNCTION minus_psi
!================================================================
!
!         psi*R ....
!
!================================================================
          !!@description: TODO
          FUNCTION R_time_psi(R,psi)
            CLASS (param_psi), intent (in) :: psi
            real (kind=Rkind),    intent (in) :: R
            TYPE (param_psi) :: R_time_psi
            integer           :: err,i

!           write(out_unitp,*) 'BEGINNING R_time_psi'

!           - define and allocate R_time_psi ----
            CALL copy_psi2TOpsi1(R_time_psi,psi)
!           -----------------------------------------


            IF (psi%GridRep) THEN
              IF (psi%cplx) THEN
                R_time_psi%CvecG = psi%CvecG * cmplx(R,kind=Rkind)
              ELSE
                R_time_psi%RvecG = psi%RvecG * R
              END IF
            END IF

            IF (psi%BasisRep) THEN
              IF (psi%cplx) THEN
                R_time_psi%CvecB = psi%CvecB * cmplx(R,kind=Rkind)
              ELSE
                R_time_psi%RvecB = psi%RvecB * R
              END IF
            END IF

            IF(psi%SRG_MPI) THEN
              R_time_psi%SR_G(:,1)=psi%SR_G(:,1)*R
              IF(psi%cplx) R_time_psi%SR_G(:,2)=psi%SR_G(:,2)*R
            ENDIF

            IF(psi%SRB_MPI) THEN
              R_time_psi%SR_B(:,1)=psi%SR_B(:,1)*R
              IF(psi%cplx) R_time_psi%SR_B(:,2)=psi%SR_B(:,2)*R
            ENDIF

            R_time_psi%symab = psi%symab

!           write(out_unitp,*) 'END R_time_psi'

          END FUNCTION R_time_psi

          FUNCTION psi_time_R(psi,R)
            USE mod_MPI

            CLASS (param_psi), intent (in) :: psi
            real (kind=Rkind),    intent (in) :: R
            TYPE (param_psi)  :: psi_time_R
            integer           :: err,i

!           write(out_unitp,*) 'BEGINNING psi_time_R'

!           - define and allocate psi_time_R ----
            CALL copy_psi2TOpsi1(psi_time_R,psi)
!           -----------------------------------------


            IF (psi%GridRep) THEN
              IF (psi%cplx) THEN
                IF(keep_MPI) psi_time_R%CvecG = psi%CvecG * cmplx(R,kind=Rkind)
              ELSE
                IF(keep_MPI) psi_time_R%RvecG = psi%RvecG * R
              END IF
            END IF

            IF (psi%BasisRep) THEN
              IF (psi%cplx) THEN
                IF(keep_MPI) psi_time_R%CvecB = psi%CvecB * cmplx(R,kind=Rkind)
              ELSE
                IF(keep_MPI) psi_time_R%RvecB = psi%RvecB * R
              END IF
            END IF

            IF(psi%SRG_MPI) THEN
              psi_time_R%SR_G(:,1)=psi%SR_G(:,1)*R
              IF(psi%cplx) psi_time_R%SR_G(:,2)=psi%SR_G(:,2)*R
            ENDIF

            IF(psi%SRB_MPI) THEN
              psi_time_R%SR_B(:,1)=psi%SR_B(:,1)*R
              IF(psi%cplx) psi_time_R%SR_B(:,2)=psi%SR_B(:,2)*R
            ENDIF

            psi_time_R%symab = psi%symab

!           write(out_unitp,*) 'END psi_time_R'

          END FUNCTION psi_time_R

          !!@description: TODO
          FUNCTION C_time_psi(C,psi)
            CLASS (param_psi), intent (in) :: psi
            complex (kind=Rkind), intent (in) :: C
            TYPE (param_psi)  :: C_time_psi
            integer           :: err,i

!           write(out_unitp,*) 'BEGINNING C_time_psi'

!           - define and allocate C_time_psi ----
            CALL copy_psi2TOpsi1(C_time_psi,psi)
!           -----------------------------------------


            IF (psi%BasisRep) THEN
              IF (psi%cplx) THEN
                C_time_psi%CvecB = psi%CvecB * C
              ELSE
                write(out_unitp,*) ' ERROR : in C_time_psi'
                write(out_unitp,*) ' I cannot multiply a real psi and a complex'
                write(out_unitp,*) 'psi%cplx,C',psi%cplx,C
                STOP
              END IF
            END IF

            IF (psi%GridRep) THEN
              IF (psi%cplx) THEN
                C_time_psi%CvecG = psi%CvecG * C
              ELSE
                write(out_unitp,*) ' ERROR : in C_time_psi'
                write(out_unitp,*) ' I cannot multiply a real psi and a complex'
                write(out_unitp,*) 'psi%cplx,C',psi%cplx,C
                STOP
              END IF
            END IF

            C_time_psi%symab = psi%symab

!           write(out_unitp,*) 'END C_time_psi'

          END FUNCTION C_time_psi

          !!@description: TODO
!=======================================================================================
          FUNCTION psi_time_C(psi,C)
            TYPE (param_psi)                  :: psi_time_C

            CLASS (param_psi),     intent (in) :: psi
            complex (kind=Rkind),  intent (in) :: C

            integer           :: err,i

            !write(out_unitp,*) 'BEGINNING psi_time_C'
            !write(out_unitp,*) 'psi%CvecB',psi%CvecB
            !write(out_unitp,*) 'C',C
            !CALL flush_perso(out_unitp)

!           - define and allocate psi_time_C ----
            CALL copy_psi2TOpsi1(psi_time_C,psi)
!           -----------------------------------------

            IF (psi%BasisRep) THEN
              IF (psi%cplx) THEN
                IF(keep_MPI) psi_time_C%CvecB(:) = psi%CvecB * C
              ELSE
                write(out_unitp,*) ' ERROR : in psi_time_C'
                write(out_unitp,*) ' I cannot multiply a real psi and a complex'
                write(out_unitp,*) 'psi%cplx,C',psi%cplx,C
                STOP
              END IF
            END IF

            IF (psi%GridRep) THEN
              IF (psi%cplx) THEN
                IF(keep_MPI) psi_time_C%CvecG(:) = psi%CvecG * C
              ELSE
                write(out_unitp,*) ' ERROR : in psi_time_C'
                write(out_unitp,*) ' I cannot multiply a real psi and a complex'
                write(out_unitp,*) 'psi%cplx,C',psi%cplx,C
                STOP
              END IF
            END IF

            psi_time_C%symab = psi%symab

            !write(out_unitp,*) 'END psi_time_C'
            !CALL flush_perso(out_unitp)

          END FUNCTION psi_time_C


          FUNCTION psi_over_R(psi,R)
            USE mod_MPI

            CLASS (param_psi),    intent (in) :: psi
            real (kind=Rkind),    intent (in) :: R
            TYPE (param_psi)  :: psi_over_R

!           write(out_unitp,*) 'BEGINNING psi_over_R'

!           - define and allocate psi_over_R ----
            CALL copy_psi2TOpsi1(psi_over_R,psi)
!           -----------------------------------------


            IF (psi%GridRep) THEN
              IF (psi%cplx) THEN
                IF(keep_MPI) psi_over_R%CvecG = psi%CvecG / cmplx(R,kind=Rkind)
              ELSE
                IF(keep_MPI) psi_over_R%RvecG = psi%RvecG / R
              END IF
            END IF

            IF (psi%BasisRep) THEN
              IF (psi%cplx) THEN
                IF(keep_MPI) psi_over_R%CvecB = psi%CvecB / cmplx(R,kind=Rkind)
              ELSE
                IF(keep_MPI) psi_over_R%RvecB = psi%RvecB / R
              END IF
            END IF

            psi_over_R%symab = psi%symab

!           write(out_unitp,*) 'END psi_over_R'

          END FUNCTION psi_over_R
          FUNCTION psi_over_C(psi,C)
            TYPE (param_psi)                   :: psi_over_C

            CLASS (param_psi),     intent (in) :: psi
            complex (kind=Rkind),  intent (in) :: C

            !write(out_unitp,*) 'BEGINNING psi_over_C'
            !write(out_unitp,*) 'psi%CvecB',psi%CvecB
            !write(out_unitp,*) 'C',C
            !CALL flush_perso(out_unitp)

!           - define and allocate psi_over_C ----
            CALL copy_psi2TOpsi1(psi_over_C,psi)
!           -----------------------------------------

            IF (psi%BasisRep) THEN
              IF (psi%cplx) THEN
                IF(keep_MPI) psi_over_C%CvecB(:) = psi%CvecB / C
              ELSE
                write(out_unitp,*) ' ERROR : in psi_over_C'
                write(out_unitp,*) ' I cannot divide a real psi by a complex'
                write(out_unitp,*) 'psi%cplx,C',psi%cplx,C
                STOP
              END IF
            END IF

            IF (psi%GridRep) THEN
              IF (psi%cplx) THEN
                IF(keep_MPI) psi_over_C%CvecG(:) = psi%CvecG / C
              ELSE
                write(out_unitp,*) ' ERROR : in psi_over_C'
                write(out_unitp,*) ' I cannot divide a real psi by a complex'
                write(out_unitp,*) 'psi%cplx,C',psi%cplx,C
                STOP
              END IF
            END IF

            psi_over_C%symab = psi%symab

            !write(out_unitp,*) 'END psi_over_C'
            !CALL flush_perso(out_unitp)

          END FUNCTION psi_over_C

      !---------------------------------------------------------------------------------
      !> INTERFACE: times_psi_SR_MPI
      !---------------------------------------------------------------------------------
      SUBROUTINE psi_times_R_SR_MPI(psi,R_const)
        USE mod_system
        IMPLICIT NONE

        TYPE(param_psi),                  intent(inout) :: psi
        Real(kind=Rkind),                 intent(in)    :: R_const

        IF(psi%SRG_MPI) THEN
          psi%SR_G=psi%SR_G*R_const
        ELSE IF(psi%SRB_MPI) THEN
          psi%SR_B=psi%SR_B*R_const
        ELSE
          STOP 'error in R_times_psi_SR_MPI, psi%SR*_MPI=.FALSE.'
        ENDIF

      ENDSUBROUTINE psi_times_R_SR_MPI

      !---------------------------------------------------------------------------------
      ! note: psi_new should be already well allocated here
      SUBROUTINE psi_times_C_SR_MPI(psi,C_const,psi_new)
        USE mod_system
        IMPLICIT NONE

        TYPE(param_psi),                  intent(inout) :: psi
        TYPE(param_psi),optional,         intent(inout) :: psi_new
        Complex(kind=Rkind),              intent(in)    :: C_const

        Real(kind=Rkind),allocatable                    :: SR_temp(:)
        Real(kind=Rkind)                                :: C
        Real(kind=Rkind)                                :: R


        R=Real(C_const,kind=Rkind)
        C=aimag(C_const)

        IF(psi%SRG_MPI) THEN
          IF(psi%cplx) THEN
            IF(present(psi_new)) THEN
              IF(.NOT. psi_new%SRG_MPI) STOP 'psi_new error in psi_times_C_SR_MPI'
              psi_new%SR_G(:,1)=psi%SR_G(:,1)*R-psi%SR_G(:,2)*C
              psi_new%SR_G(:,2)=psi%SR_G(:,1)*C+psi%SR_G(:,2)*R
            ELSE
              allocate(SR_temp(size(psi%SR_G,1)))
              SR_temp(:)=psi%SR_G(:,1)
              psi%SR_G(:,1)=psi%SR_G(:,1)*R-psi%SR_G(:,2)*C
              psi%SR_G(:,2)=SR_temp      *C+psi%SR_G(:,2)*R
              deallocate(SR_temp)
            ENDIF
          ELSE
            psi%SR_G(:,1)=psi%SR_G(:,1)*R
          ENDIF
        ELSE IF(psi%SRB_MPI) THEN
          IF(psi%cplx) THEN
            IF(present(psi_new)) THEN
              IF(.NOT. psi_new%SRB_MPI) STOP 'psi_new error in psi_times_C_SR_MPI'
              psi_new%SR_B(:,1)=psi%SR_B(:,1)*R-psi%SR_B(:,2)*C
              psi_new%SR_B(:,2)=psi%SR_B(:,1)*C+psi%SR_B(:,2)*R
            ELSE
              allocate(SR_temp(size(psi%SR_B,1)))
              SR_temp(:)=psi%SR_B(:,1)
              psi%SR_B(:,1)=psi%SR_B(:,1)*R-psi%SR_B(:,2)*C
              psi%SR_B(:,2)=SR_temp      *C+psi%SR_B(:,2)*R
              deallocate(SR_temp)
            ENDIF
          ELSE
            psi%SR_B(:,1)=psi%SR_B(:,1)*R
          ENDIF
        ELSE
          STOP 'error in R_times_psi_SR_MPI, psi%SR*_MPI=.FALSE.'
        ENDIF

      ENDSUBROUTINE psi_times_C_SR_MPI

      !---------------------------------------------------------------------------------
      SUBROUTINE psi_times_psi_SR_MPI(psi,psi1)
        USE mod_system
        IMPLICIT NONE

        TYPE(param_psi),                  intent(inout) :: psi
        TYPE(param_psi),                  intent(in)    :: psi1

        Real(kind=Rkind),allocatable                    :: SR_temp(:)

        IF(psi%SRG_MPI) THEN
          IF(psi%cplx) THEN
            allocate(SR_temp(size(psi%SR_G,1)))
            SR_temp=psi%SR_G(:,1)
            psi%SR_G(:,1)=psi%SR_G(:,1)*psi1%SR_G(:,1)-psi%SR_G(:,2)*psi1%SR_G(:,2)
            psi%SR_G(:,2)=SR_temp      *psi1%SR_G(:,2)+psi%SR_G(:,2)*SR_temp
            deallocate(SR_temp)
          ELSE
            psi%SR_G=psi1%SR_G*psi%SR_G
          ENDIF
        ELSE IF(psi%SRB_MPI) THEN
          IF(psi%cplx) THEN
            allocate(SR_temp(size(psi%SR_B,1)))
            SR_temp=psi%SR_B(:,1)
            psi%SR_B(:,1)=psi%SR_B(:,1)*psi1%SR_B(:,1)-psi%SR_B(:,2)*psi1%SR_B(:,2)
            psi%SR_B(:,2)=SR_temp      *psi1%SR_B(:,2)+psi%SR_B(:,2)*SR_temp
            deallocate(SR_temp)
          ELSE
            psi%SR_B=psi1%SR_B*psi%SR_B
          ENDIF
        ELSE
          STOP 'error in R_times_psi_SR_MPI, psi%SR*_MPI=.FALSE.'
        ENDIF

      ENDSUBROUTINE psi_times_psi_SR_MPI

      !---------------------------------------------------------------------------------
      !> INTERFACE: plus_psi_SR_MPI
      !---------------------------------------------------------------------------------
      SUBROUTINE psi_plus_R_SR_MPI(psi,R_const)
        USE mod_system
        IMPLICIT NONE

        TYPE(param_psi),                  intent(inout) :: psi
        Real(kind=Rkind),                 intent(in)    :: R_const

        IF(psi%SRG_MPI) THEN
          psi%SR_G(:,1)=psi%SR_G(:,1)+R_const
        ELSE IF(psi%SRB_MPI) THEN
          psi%SR_B(:,1)=psi%SR_B(:,1)+R_const
        ELSE
          STOP 'error in R_times_psi_SR_MPI, psi%SR*_MPI=.FALSE.'
        ENDIF

        IF(R_const/=ZERO) psi%symab=-1

      ENDSUBROUTINE psi_plus_R_SR_MPI

      !---------------------------------------------------------------------------------
      SUBROUTINE psi_plus_C_SR_MPI(psi,C_const)
        USE mod_system
        IMPLICIT NONE

        TYPE(param_psi),                  intent(inout) :: psi
        Complex(kind=Rkind),              intent(in)    :: C_const

        IF(psi%SRG_MPI) THEN
          psi%SR_G(:,1)=psi%SR_G(:,1)+Real(C_const,kind=Rkind)
          IF(psi%cplx) psi%SR_G(:,2)=psi%SR_G(:,2)+aimag(C_const)
        ELSE IF(psi%SRB_MPI) THEN
          psi%SR_B(:,1)=psi%SR_B(:,1)+Real(C_const,kind=Rkind)
          IF(psi%cplx) psi%SR_B(:,2)=psi%SR_B(:,2)+aimag(C_const)
        ELSE
          STOP 'error in R_times_psi_SR_MPI, psi%SR*_MPI=.FALSE.'
        ENDIF

        IF(abs(C_const)/=ZERO) psi%symab=-1

      ENDSUBROUTINE psi_plus_C_SR_MPI

      !---------------------------------------------------------------------------------
      SUBROUTINE psi_plus_psi_SR_MPI(psi,psi1)
        USE mod_system
        IMPLICIT NONE

        TYPE(param_psi),                  intent(inout) :: psi
        TYPE(param_psi),                  intent(in)    :: psi1

        IF(psi%SRG_MPI) THEN
          psi%SR_G(:,:)=psi%SR_G(:,:)+psi1%SR_G(:,:)
        ELSE IF(psi%SRB_MPI) THEN
          psi%SR_B(:,:)=psi%SR_B(:,:)+psi1%SR_B(:,:)
        ELSE
          STOP 'error in R_times_psi_SR_MPI, psi%SR*_MPI=.FALSE.'
        ENDIF

        IF(psi%symab/=psi%symab) THEN
          IF(psi%symab==-2) THEN
            psi%symab=psi1%symab
          ELSE
            psi%symab=-1
          ENDIF
        ENDIF

      ENDSUBROUTINE psi_plus_psi_SR_MPI

      !---------------------------------------------------------------------------------
      !> INTERFACE: minus_psi_SR_MPI
      !---------------------------------------------------------------------------------
      SUBROUTINE psi_minus_R_SR_MPI(psi,R_const)
        USE mod_system
        IMPLICIT NONE

        TYPE(param_psi),                  intent(inout) :: psi
        Real(kind=Rkind),                 intent(in)    :: R_const

        IF(psi%SRG_MPI) THEN
          psi%SR_G(:,1)=psi%SR_G(:,1)-R_const
        ELSE IF(psi%SRB_MPI) THEN
          psi%SR_B(:,1)=psi%SR_B(:,1)-R_const
        ELSE
          STOP 'error in R_times_psi_SR_MPI, psi%SRG*_MPI=.FALSE.'
        ENDIF

        IF(R_const/=ZERO) psi%symab=-1

      ENDSUBROUTINE psi_minus_R_SR_MPI

      !---------------------------------------------------------------------------------
      SUBROUTINE psi_minus_C_SR_MPI(psi,C_const)
        USE mod_system
        IMPLICIT NONE

        TYPE(param_psi),                  intent(inout) :: psi
        Complex(kind=Rkind),              intent(in)    :: C_const

        IF(psi%SRG_MPI) THEN
          psi%SR_G(:,1)=psi%SR_G(:,1)-Real(C_const,kind=Rkind)
          IF(psi%cplx) psi%SR_G(:,2)=psi%SR_G(:,2)-aimag(C_const)
        ELSE IF(psi%SRB_MPI) THEN
          psi%SR_B(:,1)=psi%SR_B(:,1)-Real(C_const)
          IF(psi%cplx) psi%SR_B(:,2)=psi%SR_B(:,2)-aimag(C_const)
        ELSE
          STOP 'error in R_times_psi_SR_MPI, psi%SR*_MPI=.FALSE.'
        ENDIF

        IF(abs(C_const)/=ZERO) psi%symab=-1

      ENDSUBROUTINE psi_minus_C_SR_MPI

      !---------------------------------------------------------------------------------
      SUBROUTINE psi_minus_psi_SR_MPI(psi,psi1)
        USE mod_system
        IMPLICIT NONE

        TYPE(param_psi),                  intent(inout) :: psi
        TYPE(param_psi),                  intent(in)    :: psi1

        IF(psi%SRG_MPI) THEN
          psi%SR_G(:,:)=psi%SR_G(:,:)-psi1%SR_G(:,:)
        ELSE IF(psi%SRB_MPI) THEN
          psi%SR_B(:,:)=psi%SR_B(:,:)-psi1%SR_B(:,:)
        ELSE
          STOP 'error in R_times_psi_SR_MPI, psi%SR*_MPI=.FALSE.'
        ENDIF

        IF(psi%symab/=psi%symab) THEN
          IF(psi%symab==-2) THEN
            psi%symab=psi1%symab
          ELSE
            psi%symab=-1
          ENDIF
        ENDIF

      ENDSUBROUTINE psi_minus_psi_SR_MPI

!=======================================================================================
!================================================================
!
!     write the initialization parameters
!
!================================================================
      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      SUBROUTINE ecri_init_psi(psi)

      USE mod_system
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_psi)   :: psi


      write(out_unitp,*) ' BEGINNING Write_init_psi'

      write(out_unitp,*) 'para_AllBasis is linked?',associated(psi%para_AllBasis)
      write(out_unitp,*) 'BasisnD       is linked?',associated(psi%BasisnD)
      write(out_unitp,*) 'Basis2n       is linked?',associated(psi%Basis2n)
      write(out_unitp,*) 'symab, bits(symab)',WriteTOstring_symab(psi%symab)

      write(out_unitp,*) 'init,cplx',psi%init,psi%cplx

      write(out_unitp,*)
      write(out_unitp,*) 'nb_tot          ',psi%nb_tot
      write(out_unitp,*) 'nb_tot_contrac  ',psi%nb_tot_contrac
      write(out_unitp,*) 'nb_tot_uncontrac',psi%nb_tot_uncontrac
      write(out_unitp,*) 'BasisRep_contrac',psi%BasisRep_contrac

      write(out_unitp,*)
      write(out_unitp,*) 'nb_ba,nb_bi,nb_be,nb_bRot',psi%nb_ba,psi%nb_bi,psi%nb_be,psi%nb_bRot
      write(out_unitp,*) 'nb_qa',psi%nb_qa
      write(out_unitp,*) 'nb_baie,nb_qaie,nb_paie',psi%nb_baie,psi%nb_qaie,psi%nb_paie
      write(out_unitp,*)
      write(out_unitp,*) 'nb_act1,nb_act',psi%nb_act1,psi%nb_act
      write(out_unitp,*) 'nb_basis_act1,nb_basis',psi%nb_basis_act1,psi%nb_basis
      write(out_unitp,*) 'max_dim',psi%max_dim
      write(out_unitp,*)

      write(out_unitp,*) 'nb_TDParam',psi%nb_TDParam
      write(out_unitp,*)


      write(out_unitp,*) 'RvecB_alloc',allocated(psi%RvecB)
      write(out_unitp,*) 'CvecB_alloc',allocated(psi%CvecB)
      write(out_unitp,*) 'RvecG_alloc',allocated(psi%RvecG)
      write(out_unitp,*) 'CvecG_alloc',allocated(psi%CvecG)
      write(out_unitp,*) 'TDParam_alloc',allocated(psi%TDParam)

      write(out_unitp,*) 'IndAvOp,CAvOp,convAvOp',psi%IndAvOp,psi%CAvOp,psi%convAvOp

      write(out_unitp,*)
      write(out_unitp,*) 'norm^2',psi%norm2


      write(out_unitp,*) ' END ecri_init_psi'

      END SUBROUTINE ecri_init_psi
!===============================================================================
!     writing psi
!===============================================================================
      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      SUBROUTINE ecri_psi(T,psi,nioWP,ecri_GridRep,ecri_BasisRep,       &
                          channel_all,channel_pack,ecri_nume,ecri_numi, &
                          ecri_psi2)

!----- variables for the WP ----------------------------------------
      TYPE (param_psi)            :: psi
      real (kind=Rkind), optional :: T ! time
      integer,           optional :: nioWP

      logical,           optional :: ecri_GridRep,ecri_BasisRep
      logical,           optional :: channel_all,channel_pack
      logical,           optional :: ecri_psi2
      integer,           optional :: ecri_nume,ecri_numi



!----- working variables ------------------------------------------
      real (kind=Rkind) :: locT ! time
      integer           :: loc_nioWP

      logical           :: loc_ecri_GridRep,loc_ecri_BasisRep
      logical           :: loc_channel_all,loc_channel_pack
      logical           :: loc_ecri_psi2
      logical           :: new = .TRUE.
      integer           :: loc_ecri_nume,loc_ecri_numi
      real (kind=Rkind) :: x(psi%BasisnD%ndim),x2

      real (kind=Rkind),    allocatable :: psiQchannel(:)
      real (kind=Rkind),    allocatable :: RVec(:)
      complex (kind=Rkind), allocatable :: CVec(:)

      real (kind=Rkind) :: psi2

      integer       :: i_d,i_wp,nb_ei
      integer       :: i_qa,i_qaie
      integer       :: i_ba,i_bi,i_be,i_baie
      integer       :: ii_bi,ii_be,if_bi,if_be

      character(len=:), allocatable     :: fformat
      character (len=*),parameter       :: Rformat='f20.10'
      integer            :: ni

      logical           :: Write_Grid_done,Write_Basis_done
      TYPE (OldParam)   :: OldPara


!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical,parameter :: debug = .FALSE.
!     logical,parameter :: debug = .TRUE.
!-----------------------------------------------------------

      IF (present(T)) THEN
        locT = T
      ELSE
        locT = ZERO
      END IF

      IF (present(nioWP)) THEN
        loc_nioWP = nioWP
      ELSE
        loc_nioWP = out_unitp
      END IF

      IF (present(ecri_GridRep)) THEN
        loc_ecri_GridRep = ecri_GridRep
      ELSE
        loc_ecri_GridRep = psi%GridRep
      END IF

      IF (present(ecri_BasisRep)) THEN
        loc_ecri_BasisRep = ecri_BasisRep
      ELSE
        loc_ecri_BasisRep = psi%BasisRep
      END IF

      IF (present(channel_all)) THEN
        loc_channel_all = channel_all
      ELSE
        loc_channel_all = .TRUE.
      END IF

      IF (present(channel_pack)) THEN
        loc_channel_pack = channel_pack
      ELSE
        loc_channel_pack = .FALSE.
      END IF

      IF (present(ecri_psi2)) THEN
        loc_ecri_psi2 = ecri_psi2
      ELSE
        loc_ecri_psi2 = .FALSE.
      END IF

      IF (present(ecri_nume)) THEN
        loc_ecri_nume = ecri_nume
      ELSE
        loc_ecri_nume = 1
      END IF

      IF (present(ecri_numi)) THEN
        loc_ecri_numi = ecri_numi
      ELSE
        loc_ecri_numi = 1
      END IF

      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ecri_psi'
        write(out_unitp,*) 'T',locT
        write(out_unitp,*) 'ecri_GridRep,ecri_BasisRep',loc_ecri_GridRep,loc_ecri_BasisRep
        write(out_unitp,*) 'channel_all,channel_pack',                          &
                    loc_channel_all,loc_channel_pack
        write(out_unitp,*)
        write(out_unitp,*) 'nioWP',loc_nioWP
        write(out_unitp,*) 'ecri_psi2',loc_ecri_psi2
        write(out_unitp,*) 'ecri_nume,ecri_numi',loc_ecri_nume,loc_ecri_numi
        write(out_unitp,*)
      END IF
!-----------------------------------------------------------


      IF (loc_ecri_numi <= 0) loc_ecri_numi = 1
      IF (loc_ecri_nume <= 0) loc_ecri_nume = 1

      IF (loc_channel_all) THEN
        ii_bi=1
        ii_be=1
        if_bi=psi%nb_bi
        if_be=psi%nb_be
      ELSE
        ii_bi=ecri_numi
        ii_be=ecri_nume
        if_bi=ecri_numi
        if_be=ecri_nume
      END IF

      Write_Grid_done  = .FALSE.
      Write_Basis_done = .FALSE.

      IF (psi%nb_baie == psi%nb_tot .AND. loc_ecri_GridRep .AND.            &
          (allocated(psi%CvecG) .OR. allocated(psi%RvecG)) ) THEN

        IF (loc_channel_pack) THEN
          ni=2+psi%nb_act1

          DO i_qa=1,psi%nb_qa

            psi2 = ZERO
            DO i_be=1,psi%nb_be
            DO i_bi=1,psi%nb_bi

              i_qaie = i_qa + ( (i_bi-1)+ (i_be-1)*psi%nb_bi ) * psi%nb_qa
              IF (psi%cplx) THEN
                psi2 = psi2 + abs(psi%CvecG(i_qaie))**2
              ELSE
                psi2 = psi2 + psi%RvecG(i_qaie)**2
              END IF

            END DO
            END DO

!           -x(psi%BasisnD%ndim) calculation -----------------------------------
            CALL Rec_x(x,psi%BasisnD,i_qa)

!           - write psi2 ----------------------------------------
            IF (i_qa == psi%nb_qa) THEN
              fformat = String_TO_String( '(' // int_TO_char(ni) // Rformat // '/)' )
            ELSE
              fformat = String_TO_String( '(' // int_TO_char(ni) // Rformat //  ')' )
            END IF

            IF (psi%BasisnD%ndim > 1) THEN
              IF (i_qa == 1) THEN
                x2 = x(psi%BasisnD%ndim-1)
              END IF
              IF (x2 /=  x(psi%BasisnD%ndim-1)) THEN
                x2 = x(psi%nb_act1-1)
                write(loc_nioWP,*)
              END IF
            END IF
            write(loc_nioWP,fformat) locT,x,psi2

            deallocate(fformat)
          END DO

        ELSE IF (new) THEN
          nb_ei = (if_bi-ii_bi+1)*(if_be-ii_be+1)
          !write(out_unitp,*) 'nb_ei',nb_ei,ii_be,if_be,ii_bi,if_bi
          ni = 1+psi%nb_act1
          IF (loc_ecri_psi2) THEN
            ni = ni + nb_ei
            CALL alloc_NParray(psiQchannel,(/nb_ei/),"psiQchannel","ecri_psi")
          ELSE
            ni = ni + 2*nb_ei
            CALL alloc_NParray(psiQchannel,(/2*nb_ei/),"psiQchannel","ecri_psi")
          END IF

          IF (psi%cplx) THEN
            CALL alloc_NParray(CVec,(/psi%nb_bi*psi%nb_be/),"CVec","ecri_psi")
          ELSE
            CALL alloc_NParray(RVec,(/psi%nb_bi*psi%nb_be/),"RVec","ecri_psi")
          END IF

          DO i_qa=1,psi%nb_qa

            IF (psi%cplx) THEN
              CALL get_CVec_OF_psi_AT_ind_a(CVec,psi,i_qa,OldPara=OldPara)
              IF (loc_ecri_psi2) THEN
                psiQchannel(:) = abs(CVec)**2
              ELSE
                psiQchannel(1:2*nb_ei:2) = real(CVec,kind=Rkind)
                psiQchannel(2:2*nb_ei:2) = aimag(CVec)
              END IF
            ELSE
              CALL get_RVec_OF_psi_AT_ind_a(RVec,psi,i_qa,OldPara=OldPara)
              IF (loc_ecri_psi2) THEN
                psiQchannel(:) = RVec**2
              ELSE
                psiQchannel(:) = ZERO
                psiQchannel(1:2*nb_ei:2) = RVec
              END IF
            END IF

!-x(psi%BasisnD%ndim) calculation -----------------------------------
            CALL Rec_x(x,psi%BasisnD,i_qa)

!- write psi2 ----------------------------------------
            ! IF (i_qa == psi%nb_qa) THEN
            !   fformat = '(' // int_TO_char(ni) // Rformat // '/)'
            ! ELSE
            !   fformat = '(' // int_TO_char(ni) // Rformat //  ')'
            ! END IF
            IF (i_qa == psi%nb_qa) THEN
              fformat = '(' // int_TO_char(ni) // 'f8.3/)'
            ELSE
              fformat = '(' // int_TO_char(ni) // 'f8.3)'
            END IF

            IF (psi%BasisnD%ndim > 1) THEN
              IF (i_qa == 1) THEN
                x2 = x(psi%BasisnD%ndim-1)
              END IF
              IF (x2 /=  x(psi%BasisnD%ndim-1)) THEN
                x2 = x(psi%nb_act1-1)
                write(loc_nioWP,*)
              END IF
            END IF
            write(loc_nioWP,fformat) locT,x,psiQchannel(:)

            deallocate(fformat)

          END DO
          CALL dealloc_NParray(psiQchannel,"psiQchannel","ecri_psi")
          IF (allocated(CVec)) CALL dealloc_NParray(CVec,"CVec","ecri_psi")
          IF (allocated(RVec)) CALL dealloc_NParray(RVec,"RVec","ecri_psi")

        ELSE
          ni=2+psi%nb_act1
          IF (.NOT. loc_ecri_psi2 .AND. psi%cplx) ni=3+psi%nb_act1

          DO i_be=ii_be,if_be
          DO i_bi=ii_bi,if_bi

            DO i_qa=1,psi%nb_qa
              i_qaie = i_qa + ( (i_bi-1)+ (i_be-1)*psi%nb_bi ) *        &
                                                            psi%nb_qa

              !-x(psi%BasisnD%ndim) calculation ------------------------
              CALL Rec_x(x,psi%BasisnD,i_qa)

              IF (i_qa == psi%nb_qa) THEN
                fformat = String_TO_String( '(' // Rformat // ',3i8,' // int_TO_char(ni) // Rformat // '/)' )
              ELSE
                fformat = String_TO_String( '(' // Rformat // ',3i8,' // int_TO_char(ni) // Rformat //  ')' )
              END IF

              IF (loc_ecri_psi2) THEN
                IF (psi%cplx) THEN
                  write(loc_nioWP,fformat) locT,i_be,i_bi,i_qa,x,abs(psi%CvecG(i_qaie))**2
                ELSE
                  write(loc_nioWP,fformat) locT,i_be,i_bi,i_qa,x,psi%RvecG(i_qaie)**2
                END IF
              ELSE
                IF (psi%cplx) THEN
                  write(loc_nioWP,fformat) locT,i_be,i_bi,i_qa,x,psi%CvecG(i_qaie)
                ELSE
                  write(loc_nioWP,fformat) locT,i_be,i_bi,i_qa,x,psi%RvecG(i_qaie)
                END IF
              END IF

              deallocate(fformat)

            END DO
          END DO
          END DO
        END IF
        Write_Grid_done  = .TRUE.

      END IF

      IF (psi%nb_baie == psi%nb_tot .AND. loc_ecri_BasisRep .AND.       &
          (allocated(psi%CvecB) .OR. allocated(psi%RvecB)) ) THEN

        ni=1
        IF (.NOT. loc_ecri_psi2 .AND. psi%cplx) ni=2

        DO i_be=ii_be,if_be
        DO i_bi=ii_bi,if_bi
        DO i_ba=1,psi%nb_ba
          i_baie = i_ba + ( (i_bi-1)+ (i_be-1)*psi%nb_bi ) * psi%nb_ba

          IF (i_ba == psi%nb_ba) THEN
            fformat = String_TO_String( '(' // Rformat // ',3i8,' // int_TO_char(ni) // Rformat // '/)' )
          ELSE
            fformat = String_TO_String( '(' // Rformat // ',3i8,' // int_TO_char(ni) // Rformat //  ')' )
          END IF

          IF (loc_ecri_psi2) THEN
            IF (psi%cplx) THEN
              write(loc_nioWP,fformat) locT,i_be,i_bi,i_ba,abs(psi%CvecB(i_baie))**2
            ELSE
              write(loc_nioWP,fformat) locT,i_be,i_bi,i_ba,psi%RvecB(i_baie)**2
            END IF
          ELSE
            IF (psi%cplx) THEN
              write(loc_nioWP,fformat) locT,i_be,i_bi,i_ba,psi%CvecB(i_baie)
            ELSE
              write(loc_nioWP,fformat) locT,i_be,i_bi,i_ba,psi%RvecB(i_baie)
            END IF
          END IF

          deallocate(fformat)


        END DO
        END DO
        END DO
      ELSE IF (psi%nb_baie > psi%nb_tot .AND. loc_ecri_BasisRep .AND.        &
                         (allocated(psi%CvecB) .OR. allocated(psi%RvecB)) ) THEN
        ni=1
        IF (.NOT. loc_ecri_psi2 .AND. psi%cplx) ni=2

        i_be =1
        i_bi =1
        DO i_baie=1,psi%nb_tot

          IF (i_baie == psi%nb_tot) THEN
            fformat = String_TO_String( '(' // Rformat // ',3i8,' // int_TO_char(ni) // Rformat // '/)' )
          ELSE
            fformat = String_TO_String( '(' // Rformat // ',3i8,' // int_TO_char(ni) // Rformat //  ')' )
          END IF

          IF (loc_ecri_psi2) THEN
            IF (psi%cplx) THEN
              write(loc_nioWP,fformat) locT,i_be,i_bi,i_baie,abs(psi%CvecB(i_baie))**2
            ELSE
              write(loc_nioWP,fformat) locT,i_be,i_bi,i_baie,psi%RvecB(i_baie)**2
            END IF
          ELSE
            IF (psi%cplx) THEN
            write(loc_nioWP,fformat) locT,i_be,i_bi,i_baie,psi%CvecB(i_baie)
            ELSE
            write(loc_nioWP,fformat) locT,i_be,i_bi,i_baie,psi%RvecB(i_baie)
            END IF
          END IF

          deallocate(fformat)

        END DO
        Write_Basis_done = .TRUE.
      END IF

      IF (.NOT. Write_Basis_done .AND. .NOT. Write_Grid_done) THEN
        IF (loc_nioWP == 6) THEN
          write(loc_nioWP,*) ' WARNING in ecri_psi'
          IF (loc_ecri_GridRep)                                             &
            write(loc_nioWP,*) ' impossible to write CvecG or RvecG'
          IF (loc_ecri_BasisRep)                                             &
            write(loc_nioWP,*) ' impossible to write CvecB or RvecB'
        END IF
      END IF


!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'END ecri_psi'
       END IF
!-----------------------------------------------------------


      END SUBROUTINE ecri_psi
!===============================================================================

  END MODULE mod_psi_set_alloc
