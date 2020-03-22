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
 MODULE mod_psi_set_alloc
 USE mod_system
 USE mod_basis
 USE mod_ComOp
 USE mod_type_ana_psi
 IMPLICIT NONE

        TYPE param_psi
        !-- for the basis set ----------------------------------------------
        TYPE (Basis), pointer       :: BasisnD => null()   ! .TRUE. pointer
        TYPE (Basis), pointer       :: Basis2n => null()   ! .TRUE. pointer
        integer                     :: symab   = -1
        !-- for ComOp ------------------------------------------------------
        TYPE (param_ComOp), pointer :: ComOp => null()     ! .TRUE. pointer

        logical       :: builtINsub = .FALSE.            ! (F) for the use of memory collector

        logical       :: init       = .FALSE.            ! (F) all constantes are initialized
        logical       :: cplx       = .FALSE.            ! (T) if T => psi is complex (use of Cvec)
                                                         ! otherwise psi is real (use of Rvec)
        logical       :: BasisRep   = .TRUE.             ! (T) BasisRep
        logical       :: GridRep    = .FALSE.            ! (F) GridRep psi

        integer       :: nb_ba      =0
        integer       :: nb_bi      =0
        integer       :: nb_be      =0
        integer       :: nb_bRot    =0
        integer       :: nb_baie    =0
        integer       :: nb_tot     =0

        integer       :: nb_qa      =0
        integer       :: nb_qaie    =0
        integer       :: nb_paie    =0

        integer       :: nb_act1    =0
        integer       :: nb_act     =0

        integer       :: nb_basis_act1 =0
        integer       :: nb_basis      =0
        integer       :: max_dim       =0

        real (kind=Rkind),    allocatable :: RvecB(:) ! RvecB(nb_tot)
        complex (kind=Rkind), allocatable :: CvecB(:) ! CvecB(nb_tot)

        real (kind=Rkind),    allocatable :: RvecG(:) ! RvecG(nb_qaie)
        complex (kind=Rkind), allocatable :: CvecG(:) ! CvecG(nb_qaie)

        complex (kind=Rkind) :: CAvOp    = (ZERO,ZERO) ! average value for an operator (usualy H)
        integer              :: IndAvOp  = -1          ! operator type  (usualy H, IndAvOp=0)
        logical              :: convAvOp = .FALSE.     ! TRUE, if the average value is calculated with a converged wavefunction

        real (kind=Rkind)              :: norm2             = -ONE       ! norm^2 of psi

        END TYPE param_psi

        !!@description: TODO
        INTERFACE assignment (=)
          MODULE PROCEDURE psi2TOpsi1, CplxPsi_TO_RCpsi, RCPsi_TO_CplxPsi
        END INTERFACE

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

 CONTAINS

!=======================================================================================
!     allocation of psi
!     deallocation of psi
!=======================================================================================
      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      SUBROUTINE alloc_psi(psi,BasisRep,GridRep)
      USE mod_MPI

      TYPE (param_psi), intent(inout) :: psi
      logical, optional, intent(in)   :: BasisRep,GridRep

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
      IF (present(BasisRep)) psi%BasisRep = BasisRep
      IF (present(GridRep))  psi%GridRep  = GridRep

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
            IF(MPI_id==0) CALL alloc_NParray(psi%CvecB,(/psi%nb_tot/),'psi%CvecB','alloc_psi')
            IF (debug) write(out_unitp,*) 'alloc: CvecB'
            IF(MPI_id==0) psi%CvecB(:) = CZERO
          END IF
          IF ( allocated(psi%RvecB) ) THEN
            CALL dealloc_NParray(psi%RvecB,'psi%RvecB','alloc_psi')
            IF (debug) write(out_unitp,*) 'dealloc: RvecB'
          END IF
        ELSE
          IF ( .NOT. allocated(psi%RvecB) ) THEN
            IF(MPI_id==0) CALL alloc_NParray(psi%RvecB,(/psi%nb_tot/),'psi%RvecB','alloc_psi')
            IF (debug) write(out_unitp,*) 'alloc: RvecB'
            IF(MPI_id==0) psi%RvecB(:) = ZERO
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
            IF(MPI_id==0) CALL alloc_NParray(psi%CvecG,(/psi%nb_qaie/),'psi%CvecG','alloc_psi')
            IF (debug) write(out_unitp,*) 'alloc: CvecG'
            IF(MPI_id==0) psi%CvecG(:) = CZERO
          END IF
          IF ( allocated(psi%RvecG) ) THEN
            CALL dealloc_NParray(psi%RvecG,'psi%RvecG','alloc_psi')
            IF (debug) write(out_unitp,*) 'dealloc: RvecG'
          END IF
        ELSE
          IF ( .NOT. allocated(psi%RvecG) ) THEN
            IF(MPI_id==0) CALL alloc_NParray(psi%RvecG,(/psi%nb_qaie/),'psi%RvecG','alloc_psi')
            IF (debug) write(out_unitp,*) 'alloc: RvecG'
            IF(MPI_id==0) psi%RvecG(:) = ZERO
          END IF
          IF ( allocated(psi%CvecG) ) THEN
            CALL dealloc_NParray(psi%CvecG,'psi%CvecG','alloc_psi')
            IF (debug) write(out_unitp,*) 'dealloc: CvecG'
          END IF
        END IF
      ELSE ! deallocate psi%GridRep
        IF ( allocated(psi%RvecG) ) THEN
          CALL dealloc_NParray(psi%RvecG,'psi%RvecG','alloc_psi')
          IF (debug) write(out_unitp,*) 'dealloc: RvecG'
        END IF
        IF ( allocated(psi%CvecG) ) THEN
          CALL dealloc_NParray(psi%CvecG,'psi%CvecG','alloc_psi')
          IF (debug) write(out_unitp,*) 'dealloc: CvecG'
        END IF
      END IF

      IF (debug) write(out_unitp,*) ' 2 test_alloc_psi'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END alloc_psi'
      END IF
!-----------------------------------------------------------

      END SUBROUTINE alloc_psi
!=======================================================================================

  FUNCTION print_alloc_psi(psi)
  USE mod_string
  character(len=:), allocatable     :: print_alloc_psi

  TYPE (param_psi), intent(in) :: psi


  print_alloc_psi = String_TO_String(                                   &
          'RvecB:' // logical_TO_char(allocated(psi%RvecB)) //          &
         ',RvecG:' // logical_TO_char(allocated(psi%RvecG)) //          &
         ',CvecB:' // logical_TO_char(allocated(psi%CvecB)) //          &
         ',CvecG:' // logical_TO_char(allocated(psi%CvecG)) )

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


      psi%symab = -1


      IF (delete_all_loc) THEN
        IF (debug) write(out_unitp,*) 'dealloc: init param'

        nullify(psi%ComOp)

        psi%builtINsub      = .FALSE.

        psi%init            = .FALSE.
        psi%cplx            = .FALSE.
        !psi%BasisRep        = .TRUE.
        !psi%GridRep         = .FALSE.

        psi%nb_tot          = 0
        psi%nb_paie         = 0
        psi%nb_baie         = 0
        psi%nb_ba           = 0
        psi%nb_bi           = 0
        psi%nb_be           = 0
        psi%nb_bRot         = 0
        psi%nb_qa           = 0
        psi%nb_qaie         = 0

        psi%nb_act1         = 0
        psi%nb_act          = 0
        psi%nb_basis_act1   = 0
        psi%nb_basis        = 0
        psi%max_dim         = 0

        psi%IndAvOp         = -1
        psi%CAvOp           = (ZERO,ZERO)
        psi%convAvOp        = .FALSE.

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

      TYPE (param_psi), pointer, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

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

      TYPE (param_psi), allocatable, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

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

!----- for debuging --------------------------------------------------
      !logical, parameter :: debug = .TRUE.
       logical, parameter :: debug = .FALSE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING get_nb_be_FROM_psi'
      END IF

     IF (NewBasisEl) THEN
       get_nb_be_FROM_psi = get_nb_be_FROM_basis(psi%BasisnD)
       IF (get_nb_be_FROM_psi == -1) THEN
         write(out_unitp,*) ' ERROR in get_nb_be_FROM_psi'
         write(out_unitp,*) ' NewBasisEl = t and no El basis set!!'
         write(out_unitp,*) ' CHECK the source'
         STOP
       END IF
     ELSE IF (associated(psi%ComOp)) THEN
       get_nb_be_FROM_psi = psi%ComOp%nb_be
     ELSE
       write(out_unitp,*) ' ERROR in get_nb_be_FROM_psi'
       write(out_unitp,*) ' psi%ComOp in not associated'
       write(out_unitp,*) ' CHECK the source'
       STOP
     END IF

      IF (debug) THEN
        write(out_unitp,*) 'END get_nb_be_FROM_psi'
      END IF
      END FUNCTION get_nb_be_FROM_psi
      FUNCTION get_nb_bi_FROM_psi(psi)
      USE mod_system
      IMPLICIT NONE

      integer                      :: get_nb_bi_FROM_psi
      TYPE (param_psi), intent(in) :: psi

!----- for debuging --------------------------------------------------
      !logical, parameter :: debug = .TRUE.
       logical, parameter :: debug = .FALSE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING get_nb_bi_FROM_psi'
      END IF

     IF (associated(psi%ComOp)) THEN
       get_nb_bi_FROM_psi = psi%ComOp%nb_bi
     ELSE
       write(out_unitp,*) ' ERROR in get_nb_bi_FROM_psi'
       write(out_unitp,*) ' psi%ComOp in not associated'
       write(out_unitp,*) ' CHECK the source'
       STOP
     END IF

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
      TYPE (param_psi), intent(inout) :: psi1
      TYPE (param_psi), intent(in)    :: psi2


      integer  :: i,n1

!     write(out_unitp,*) 'BEGINNING psi2TOpsi1'

      CALL dealloc_psi(psi1)

      IF (associated(psi2%BasisnD)) THEN
         psi1%BasisnD => psi2%BasisnD
      ELSE
         write(out_unitp,*) ' ERROR in psi2TOpsi1'
         write(out_unitp,*) ' BasisnD CANNOT be associated'
         write(out_unitp,*) ' asso psi2%BasisnD',associated(psi2%BasisnD)
         write(out_unitp,*) ' CHECK the source'
         STOP
      END IF

      IF (associated(psi2%Basis2n)) THEN
         psi1%Basis2n => psi2%Basis2n
      ELSE
         write(out_unitp,*) ' ERROR in psi2TOpsi1'
         write(out_unitp,*) ' BasisnD CANNOT be associated'
         write(out_unitp,*) ' asso psi2%Basis2n',associated(psi2%Basis2n)
         write(out_unitp,*) ' CHECK the source'
         STOP
      END IF

      IF (associated(psi2%ComOp)) THEN
         psi1%ComOp => psi2%ComOp
      ELSE
         write(out_unitp,*) ' ERROR in psi2TOpsi1'
         write(out_unitp,*) ' ComOp CANNOT be associated'
         write(out_unitp,*) ' asso psi2%ComOp',associated(psi2%ComOp)
         write(out_unitp,*) ' CHECK the source'
         STOP
      END IF

      psi1%init = psi2%init
      psi1%cplx = psi2%cplx

      psi1%BasisRep  = psi2%BasisRep
      psi1%GridRep   = psi2%GridRep

      psi1%symab = psi2%symab

      psi1%nb_tot  = psi2%nb_tot
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


     psi1%IndAvOp         = psi2%IndAvOp
     psi1%CAvOp           = psi2%CAvOp
     psi1%convAvOp        = psi2%convAvOp

      psi1%norm2          = psi2%norm2

      IF (psi2%builtINsub) CALL dealloc_psi(psi2)

      END SUBROUTINE psi2TOpsi1

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
!=======================================================================================
      SUBROUTINE copy_psi2TOpsi1(psi1,psi2,BasisRep,GridRep,alloc)

!----- variables for the WP propagation ----------------------------
      TYPE (param_psi),intent(inout) :: psi1
      TYPE (param_psi), intent(in)   :: psi2
      logical, intent(in), optional  :: BasisRep
      logical, intent(in), optional  :: GridRep
      logical, intent(in), optional  :: alloc


      logical  :: loc_alloc
      integer  :: i,n1

!     write(out_unitp,*) 'BEGINNING copy_psi2TOpsi1'

      CALL dealloc_psi(psi1)

     IF (associated(psi2%BasisnD)) THEN
         psi1%BasisnD => psi2%BasisnD
      ELSE
         write(out_unitp,*) ' ERROR in copy_psi2TOpsi1'
         write(out_unitp,*) ' BasisnD CANNOT be associated'
         write(out_unitp,*) ' asso psi2%BasisnD',associated(psi2%BasisnD)
         write(out_unitp,*) ' CHECK the source'
         STOP
      END IF

      IF (associated(psi2%Basis2n)) THEN
         psi1%Basis2n => psi2%Basis2n
      ELSE
         write(out_unitp,*) ' ERROR in psi2TOpsi1'
         write(out_unitp,*) ' BasisnD CANNOT be associated'
         write(out_unitp,*) ' asso psi2%Basis2n',associated(psi2%Basis2n)
         write(out_unitp,*) ' CHECK the source'
         STOP
      END IF

      IF (associated(psi2%ComOp)) THEN
         psi1%ComOp => psi2%ComOp
      ELSE
         write(out_unitp,*) ' ERROR in copy_psi2TOpsi1'
         write(out_unitp,*) ' ComOp CANNOT be associated'
         write(out_unitp,*) ' asso psi2%ComOp',associated(psi2%ComOp)
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

      psi1%symab = psi2%symab


      psi1%nb_tot  = psi2%nb_tot
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


        psi1%IndAvOp         = psi2%IndAvOp
        psi1%CAvOp           = psi2%CAvOp
        psi1%convAvOp        = psi2%convAvOp

        psi1%norm2           = psi2%norm2

      END IF

      IF (psi2%builtINsub) CALL dealloc_psi(psi2)

!     write(out_unitp,*) 'END copy_psi2TOpsi1'

      END SUBROUTINE copy_psi2TOpsi1
!=======================================================================================

      SUBROUTINE CplxPsi_TO_RCpsi(RCPsi,Psi)
      USE mod_MPI

!----- variables for the WP propagation ----------------------------
      TYPE (param_psi),intent(inout) :: RCPsi(2)
      TYPE (param_psi), intent(in)   :: psi

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
      TYPE (param_psi), intent(in)    :: RCPsi(2)
      TYPE (param_psi), intent(inout) :: psi

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
      integer            :: i


      write(out_unitp,*) ' BEGINNING Write_init_psi'

      write(out_unitp,*) 'BasisnD, Basis2n,ComOp are linked?',          &
                  associated(psi%BasisnD),associated(psi%Basis2n),      &
                  associated(psi%ComOp)
      write(out_unitp,*) 'symab, bits(symab)',WriteTOstring_symab(psi%symab)

      write(out_unitp,*) 'builtINsub',psi%builtINsub

      write(out_unitp,*) 'init,cplx',psi%init,psi%cplx

      write(out_unitp,*)
      write(out_unitp,*) 'nb_tot',psi%nb_tot
      write(out_unitp,*)
      write(out_unitp,*) 'nb_ba,nb_bi,nb_be,nb_bRot',psi%nb_ba,psi%nb_bi,psi%nb_be,psi%nb_bRot
      write(out_unitp,*) 'nb_qa',psi%nb_qa
      write(out_unitp,*) 'nb_baie,nb_qaie,nb_paie',psi%nb_baie,psi%nb_qaie,psi%nb_paie
      write(out_unitp,*)
      write(out_unitp,*) 'nb_act1,nb_act',psi%nb_act1,psi%nb_act
      write(out_unitp,*) 'nb_basis_act1,nb_basis',psi%nb_basis_act1,psi%nb_basis
      write(out_unitp,*) 'max_dim',psi%max_dim
      write(out_unitp,*)


      write(out_unitp,*) 'RvecB_alloc',allocated(psi%RvecB)
!      IF ( allocated(psi%RvecB) .AND. psi%nb_tot > 0)                &
!       CALL Write_Vec(psi%RvecB,out_unitp,5,Rformat='e30.23')
      write(out_unitp,*) 'CvecB_alloc',allocated(psi%CvecB)
!      IF ( allocated(psi%CvecB) .AND. psi%nb_tot > 0)                &
!       CALL Write_Vec(psi%CvecB,out_unitp,5,Rformat='e30.23')
!      write(out_unitp,*)
      write(out_unitp,*) 'RvecG_alloc',allocated(psi%RvecG)
!      IF ( allocated(psi%RvecG) .AND. psi%nb_qaie .GT. 0)            &
!       CALL Write_Vec(psi%RvecG,out_unitp,5,Rformat='e30.23')
      write(out_unitp,*) 'CvecG_alloc',allocated(psi%CvecG)
!      IF ( allocated(psi%CvecG) .AND. psi%nb_qaie .GT. 0)            &
!       CALL Write_Vec(psi%CvecG,out_unitp,5,Rformat='e30.23')
!      write(out_unitp,*)

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
      integer, optional           :: nioWP

      logical, optional           :: ecri_GridRep,ecri_BasisRep
      logical, optional           :: channel_all,channel_pack
      logical, optional           :: ecri_psi2
      integer, optional           :: ecri_nume,ecri_numi



!----- working variables ------------------------------------------
      real (kind=Rkind) :: locT ! time
      integer           :: loc_nioWP

      logical           :: loc_ecri_GridRep,loc_ecri_BasisRep
      logical           :: loc_channel_all,loc_channel_pack
      logical           :: loc_ecri_psi2
      logical           :: new = .TRUE.
      integer           :: loc_ecri_nume,loc_ecri_numi
      real (kind=Rkind) :: x(psi%BasisnD%ndim),x2

      real (kind=Rkind), allocatable :: psiQchannel(:)
      real (kind=Rkind) :: psi2

      integer       :: i_d,i_wp,nb_ei
      integer       :: i_qa,i_qaie
      integer       :: i_ba,i_bi,i_be,i_baie
      integer       :: ii_bi,ii_be,if_bi,if_be

      character(len=:), allocatable     :: fformat
      character (len=*),parameter       :: Rformat='f20.10'
      integer            :: ni

      logical           :: Write_Grid_done,Write_Basis_done


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

          DO i_qa=1,psi%nb_qa

            i_d = 0
            i_wp = -1

            DO i_be=1,psi%nb_be
            DO i_bi=1,psi%nb_bi
              i_qaie = i_qa + ( (i_bi-1)+ (i_be-1)*psi%nb_bi ) * psi%nb_qa
              i_d    = i_d+1
              i_wp   = i_wp+2

              IF (psi%cplx) THEN
                IF (loc_ecri_psi2) THEN
                  psiQchannel(i_d) = abs(psi%CvecG(i_qaie))**2
                ELSE
                  psiQchannel(i_wp+0) = real(psi%CvecG(i_qaie),kind=Rkind)
                  psiQchannel(i_wp+1) = aimag(psi%CvecG(i_qaie))
                END IF
              ELSE
                IF (loc_ecri_psi2) THEN
                  psiQchannel(i_d) = psi%RvecG(i_qaie)**2
                ELSE
                  psiQchannel(i_wp+0) = psi%RvecG(i_qaie)
                  psiQchannel(i_wp+1) = ZERO
                END IF
              END IF

!             write(out_unitp,*) 'psiQchannel',i_qaie,i_wp,psiQchannel(i_wp+0),
!    *                                          psiQchannel(i_wp+1)

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
            write(loc_nioWP,fformat) locT,x,psiQchannel(:)

            deallocate(fformat)

          END DO
          CALL dealloc_NParray(psiQchannel,"psiQchannel","ecri_psi")

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
