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
MODULE mod_param_RD
USE mod_system
use mod_nDindex
IMPLICIT NONE

  PRIVATE

  TYPE param_RD
    logical                        :: RD_analysis    = .FALSE.
    integer                        :: basis_index    = 0  ! index of the basis set
    integer                        :: nb             = 0  ! size of the basis set

    TYPE (Type_nDindex)            :: nDindex_ComplBasis  ! multidimensional index for the complementary basis set
                                                          ! (without the basis set with this basis_index)
    integer                        :: nbb_ComplBasis = 0

    integer,           allocatable :: tab_OF_iBComplBasis_AND_ib_TO_iB(:,:)   ! size (nbb_ComplBasis,nb)
    real (kind=rkind), allocatable :: cbb(:,:)         ! coefficient of the contracted basis set
  CONTAINS
    PROCEDURE, PRIVATE, PASS(para_RD1) :: RD2_TO_RD1
    GENERIC,   PUBLIC  :: assignment(=) => RD2_TO_RD1
  END TYPE param_RD

PUBLIC :: param_RD,dealloc_RD,init_RD,calc_RD,calc_CRD,dealloc_tab_RD

CONTAINS

SUBROUTINE dealloc_RD(para_RD)

TYPE (param_RD), intent(inout) :: para_RD

character (len=*), parameter :: name_sub='dealloc_RD'

    para_RD%RD_analysis    = .FALSE.
    para_RD%basis_index    = 0  ! index of the basis set
    para_RD%nb             = 0  ! size of the basis set

    CALL dealloc_nDindex(para_RD%nDindex_ComplBasis)
    para_RD%nbb_ComplBasis = 0

    IF (allocated(para_RD%tab_OF_iBComplBasis_AND_ib_TO_iB)) THEN
      CALL dealloc_NParray(para_RD%tab_OF_iBComplBasis_AND_ib_TO_iB,'para_RD%tab_OF_iBComplBasis_AND_ib_TO_iB',name_sub)
    END IF

    IF (allocated(para_RD%cbb)) THEN
      CALL dealloc_NParray(para_RD%cbb,'para_RD%cbb',name_sub)
    END IF

END SUBROUTINE dealloc_RD

SUBROUTINE RD2_TO_RD1(para_RD1,para_RD2)

CLASS (param_RD), intent(inout) :: para_RD1
TYPE (param_RD),  intent(in)    :: para_RD2

integer :: i

character (len=*), parameter :: name_sub='RD2_TO_RD1'


    para_RD1%RD_analysis    = para_RD2%RD_analysis
    para_RD1%basis_index    = para_RD2%basis_index
    para_RD1%nb             = para_RD1%nb

    para_RD1%nDindex_ComplBasis = para_RD1%nDindex_ComplBasis

    para_RD1%nbb_ComplBasis = para_RD2%nbb_ComplBasis

    IF (allocated(para_RD2%tab_OF_iBComplBasis_AND_ib_TO_iB)) THEN
      CALL alloc_NParray(para_RD1%tab_OF_iBComplBasis_AND_ib_TO_iB,     &
                   shape(para_RD2%tab_OF_iBComplBasis_AND_ib_TO_iB),    &
                        'para_RD1%tab_OF_iBComplBasis_AND_ib_TO_iB',name_sub)
      para_RD1%tab_OF_iBComplBasis_AND_ib_TO_iB(:,:) =                  &
                               para_RD2%tab_OF_iBComplBasis_AND_ib_TO_iB
    END IF

    IF (allocated(para_RD2%cbb)) THEN
      CALL alloc_NParray(para_RD1%cbb,shape(para_RD2%cbb),              &
                        'para_RD1%cbb',name_sub)
      para_RD1%cbb(:,:) = para_RD2%cbb
    END IF

END SUBROUTINE RD2_TO_RD1

SUBROUTINE init_RD(para_RD,nDindB,Rvec)

TYPE (param_RD),                  intent(inout)        :: para_RD
TYPE (Type_nDindex),              intent(in)           :: nDindB    ! multidimensional index for the full basis set
real (kind=Rkind),   allocatable, intent(in), optional :: Rvec(:,:) ! real eigenvectors for the contraction

integer :: i,ibasis,IBb,IB,nbc
integer :: nDval(nDindB%ndim)

logical,parameter :: debug=.FALSE.
!logical,parameter :: debug=.TRUE.
character (len=*), parameter :: name_sub='init_RD'

  IF (.NOT. para_RD%RD_analysis) RETURN

  !-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*) 'nDindB: '
    CALL Write_nDindex(nDindB)
  END IF
  !-----------------------------------------------------------

  ibasis     = para_RD%basis_index
  !para_RD%nb = maxval(nDindB%Tab_nDval(ibasis,:))
  para_RD%nb = nDindB%nDsize(ibasis)

  IF (debug) write(out_unitp,*) 'para_RD%basis_index: ',para_RD%basis_index
  IF (debug) write(out_unitp,*) 'para_RD%nb:          ',para_RD%nb

  para_RD%nDindex_ComplBasis = nDindB

  para_RD%nDindex_ComplBasis%nDsize(ibasis) = 1
  para_RD%nDindex_ComplBasis%nDend(ibasis)  = 1

  IF (para_RD%nDindex_ComplBasis%packed) THEN
    CALL dealloc_NParray(para_RD%nDindex_ComplBasis%Tab_nDval,          &
                        "para_RD%nDindex_ComplBasis%Tab_nDval",name_sub)
    CALL dealloc_NParray(para_RD%nDindex_ComplBasis%Tab_Norm,           &
                        "para_RD%nDindex_ComplBasis%Tab_Norm",name_sub)
    CALL dealloc_NParray(para_RD%nDindex_ComplBasis%Tab_L,              &
                        "para_RD%nDindex_ComplBasis%Tab_L",name_sub)
  END IF

  CALL init_nDindexPrim(para_RD%nDindex_ComplBasis,                     &
     para_RD%nDindex_ComplBasis%ndim,para_RD%nDindex_ComplBasis%nDsize, &
     With_init=.FALSE.,With_nDindex=.TRUE.)

  IF (debug) THEN
    write(out_unitp,*) 'nDindex_ComplBasis: '
    CALL Write_nDindex(para_RD%nDindex_ComplBasis)
  END IF

  para_RD%nbb_ComplBasis = para_RD%nDindex_ComplBasis%Max_nDI

  CALL alloc_NParray(para_RD%tab_OF_iBComplBasis_AND_ib_TO_iB,          &
                    [para_RD%nb,para_RD%nbb_ComplBasis],                &
                    'para_RD%tab_OF_iBComplBasis_AND_ib_TO_iB',name_sub)

  para_RD%tab_OF_iBComplBasis_AND_ib_TO_iB(:,:) = 0

  DO iBb=1,nDindB%Max_nDI
    CALL calc_nDindex(nDindB,iBb,nDval)
    !nDval(:) = nDindB%Tab_nDval(:,iBb)
    i = nDval(ibasis)
    nDval(ibasis) = 1
    CALL calc_nDI(iB,nDval,para_RD%nDindex_ComplBasis)

    para_RD%tab_OF_iBComplBasis_AND_ib_TO_iB(i,iB) = iBb

  END DO
  IF (nDindB%Max_nDI /= count(para_RD%tab_OF_iBComplBasis_AND_ib_TO_iB > 0) ) THEN
    STOP 'Wrong init_RD'
  END IF

  IF (present(Rvec)) THEN
  IF (allocated(Rvec)) THEN

    !write(out_unitp,*) 'shape(Rvec)',shape(Rvec)

    CALL alloc_NParray(para_RD%cbb,shape(Rvec),'para_RD%cbb',name_sub)
    para_RD%cbb(:,:) = Rvec

    IF (debug) THEN
      write(out_unitp,*) 'para_RD%cbb(:,:)'
      CALL Write_Mat(para_RD%cbb,out_unitp,5)
    END IF

  END IF
  END IF

  !-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'tab_OF_iBComplBasis_AND_ib_TO_iB'
    DO iB=1,para_RD%nDindex_ComplBasis%Max_nDI
      CALL calc_nDindex(nDindB,iB,nDval)
      write(out_unitp,*) 'iB:',nDval(:),':',para_RD%tab_OF_iBComplBasis_AND_ib_TO_iB(:,iB)
    END DO
    write(out_unitp,*) 'END ',name_sub
  END IF
  !-----------------------------------------------------------

END SUBROUTINE init_RD
SUBROUTINE calc_RD(para_RD,RvecB,printRD,RD,RDcontrac)

TYPE (param_RD),                  intent(in)              :: para_RD
real (kind=Rkind),                intent(in)              :: RvecB(:) ! RvecB(nb_tot)
logical,                          intent(in),    optional :: printRD
real (kind=Rkind),   allocatable, intent(inout), optional :: RD(:,:) ! reduced density matrix with the contracted basis set (nbc,nbc)
real (kind=Rkind),   allocatable, intent(inout), optional :: RDcontrac(:,:) ! reduced density matrix with the contracted basis set (nbc,nbc)

integer :: i,j,ibasis,IBb,JBb,IB,nbc

real (kind=Rkind),  allocatable    :: RDcontrac_loc(:,:) ! reduced density matrix with the contracted basis set (nbc,nbc)
real (kind=Rkind),  allocatable    :: RD_loc(:,:) ! reduced density matrix with the contracted basis set (nbc,nbc)

logical :: printRD_loc


logical,parameter :: debug=.FALSE.
!logical,parameter :: debug=.TRUE.
character (len=*), parameter :: name_sub='calc_RD'

  IF (.NOT. para_RD%RD_analysis) RETURN

  !-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
  END IF
  !-----------------------------------------------------------

  IF (size(RvecB) /= maxval(para_RD%tab_OF_iBComplBasis_AND_ib_TO_iB)) THEN
    write(out_unitp,*) 'ERROR in ',name_sub
    write(out_unitp,*) 'the size of RvecB is not consistent with tab_OF_iBComplBasis_AND_ib_TO_iB'
    write(out_unitp,*) 'size(RvecB):',size(RvecB)
    write(out_unitp,*) 'maxval(tab_OF_iBComplBasis_AND_ib_TO_iB):',     &
                        maxval(para_RD%tab_OF_iBComplBasis_AND_ib_TO_iB)
    write(out_unitp,*) ' Check the fortran.'
    STOP ' ERROR calc_RD: inconsistent parameters'
  END IF

  printRD_loc = .FALSE.
  IF (present(printRD)) printRD_loc = printRD
  printRD_loc = printRD_loc .OR. debug

  CALL alloc_NParray(RD_loc,[para_RD%nb,para_RD%nb],'RD',name_sub)
  RD_loc(:,:) = ZERO
  DO iB=1,para_RD%nbb_ComplBasis
    DO i=1,para_RD%nb
    DO j=1,para_RD%nb
      iBb = para_RD%tab_OF_iBComplBasis_AND_ib_TO_iB(i,iB)
      jBb = para_RD%tab_OF_iBComplBasis_AND_ib_TO_iB(j,iB)
      IF (iBb > 0 .AND. jBb > 0 ) RD_loc(i,j) = RD_loc(i,j) + RvecB(iBb)*RvecB(jBb)

    END DO
    END DO
  END DO

  IF (debug) THEN
    write(out_unitp,*) 'RD_loc(:,:)'
    CALL Write_Mat(RD_loc,out_unitp,5)
  END IF

  IF (printRD_loc) write(out_unitp,*) 'Diag RD           ',para_RD%basis_index,(RD_loc(i,i),i=1,para_RD%nb)

  IF (allocated(para_RD%cbb)) THEN

    nbc = size(para_RD%cbb,dim=2)
    CALL alloc_NParray(RDcontrac_loc,[nbc,nbc],'RDcontrac',name_sub)
    RDcontrac_loc(:,:) = matmul(transpose(para_RD%cbb),matmul(RD_loc,para_RD%cbb))

    IF (printRD_loc) write(out_unitp,*) 'Diag RD contracted',para_RD%basis_index,(RDcontrac_loc(i,i),i=1,nbc)

    IF (debug) THEN
      write(out_unitp,*) 'RDcontrac_loc(:,:)'
      CALL Write_Mat(RDcontrac_loc,out_unitp,5)
    END IF

    IF (present(RDcontrac)) CALL MOVE_ALLOC( TO=RDcontrac, FROM=RDcontrac_loc )


  END IF

  IF (present(RD)) CALL MOVE_ALLOC( TO=RD, FROM=RD_loc )

  IF (allocated(RD_loc))        CALL dealloc_NParray(RD_loc,       'RD',       name_sub)
  IF (allocated(RDcontrac_loc)) CALL dealloc_NParray(RDcontrac_loc,'RDcontrac',name_sub)

  !-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'END ',name_sub
  END IF
  !-----------------------------------------------------------

END SUBROUTINE calc_RD
SUBROUTINE calc_CRD(para_RD,CvecB,printRD,CRD,CRDcontrac)

TYPE (param_RD),                   intent(in)              :: para_RD
complex (kind=Rkind),              intent(in)              :: CvecB(:) ! CvecB(nb_tot)
logical,                           intent(in),    optional :: printRD
complex (kind=Rkind), allocatable, intent(inout), optional :: CRD(:,:) ! reduced density matrix with the contracted basis set (nbc,nbc)
complex (kind=Rkind), allocatable, intent(inout), optional :: CRDcontrac(:,:) ! reduced density matrix with the contracted basis set (nbc,nbc)

integer :: i,j,ibasis,IBb,JBb,IB,nbc

complex (kind=Rkind),  allocatable    :: CRDcontrac_loc(:,:) ! reduced density matrix with the contracted basis set (nbc,nbc)
complex (kind=Rkind),  allocatable    :: CRD_loc(:,:) ! reduced density matrix with the contracted basis set (nbc,nbc)

logical :: printRD_loc


logical,parameter :: debug=.FALSE.
!logical,parameter :: debug=.TRUE.
character (len=*), parameter :: name_sub='calc_CRD'

  IF (.NOT. para_RD%RD_analysis) RETURN

  !-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
  END IF
  !-----------------------------------------------------------

  IF (size(CvecB) /= maxval(para_RD%tab_OF_iBComplBasis_AND_ib_TO_iB)) THEN
    write(out_unitp,*) 'ERROR in ',name_sub
    write(out_unitp,*) 'the size of CvecB is not consistent with tab_OF_iBComplBasis_AND_ib_TO_iB'
    write(out_unitp,*) 'size(CvecB):',size(CvecB)
    write(out_unitp,*) 'maxval(tab_OF_iBComplBasis_AND_ib_TO_iB):',     &
                        maxval(para_RD%tab_OF_iBComplBasis_AND_ib_TO_iB)
    write(out_unitp,*) ' Check the fortran.'
    STOP ' ERROR calc_CRD: inconsistent parameters'
  END IF

  printRD_loc = .FALSE.
  IF (present(printRD)) printRD_loc = printRD
  printRD_loc = printRD_loc .OR. debug

  CALL alloc_NParray(CRD_loc,[para_RD%nb,para_RD%nb],'CRD',name_sub)
  CRD_loc(:,:) = CZERO
  DO iB=1,para_RD%nbb_ComplBasis
    DO i=1,para_RD%nb
    DO j=1,para_RD%nb
      iBb = para_RD%tab_OF_iBComplBasis_AND_ib_TO_iB(i,iB)
      jBb = para_RD%tab_OF_iBComplBasis_AND_ib_TO_iB(j,iB)
      IF (iBb > 0 .AND. jBb > 0 ) CRD_loc(i,j) = CRD_loc(i,j) + conjg(CvecB(iBb))*CvecB(jBb)

    END DO
    END DO
  END DO

  IF (debug) THEN
    write(out_unitp,*) 'CRD_loc(:,:)'
    CALL Write_Mat(CRD_loc,out_unitp,5)
  END IF

  IF (printRD_loc) write(out_unitp,*) 'Diag CRD           ',para_RD%basis_index,(CRD_loc(i,i),i=1,para_RD%nb)

  IF (allocated(para_RD%cbb)) THEN

    nbc = size(para_RD%cbb,dim=2)
    CALL alloc_NParray(CRDcontrac_loc,[nbc,nbc],'CRDcontrac',name_sub)
    CRDcontrac_loc(:,:) = matmul(transpose(para_RD%cbb),matmul(CRD_loc,para_RD%cbb))

    IF (printRD_loc) write(out_unitp,*) 'Diag CRD contracted',para_RD%basis_index,(CRDcontrac_loc(i,i),i=1,nbc)

    IF (debug) THEN
      write(out_unitp,*) 'CRDcontrac_loc(:,:)'
      CALL Write_Mat(CRDcontrac_loc,out_unitp,5)
    END IF

    IF (present(CRDcontrac)) CALL MOVE_ALLOC( TO=CRDcontrac, FROM=CRDcontrac_loc )


  END IF

  IF (present(CRD)) CALL MOVE_ALLOC( TO=CRD, FROM=CRD_loc )

  IF (allocated(CRD_loc))        CALL dealloc_NParray(CRD_loc,       'CRD',       name_sub)
  IF (allocated(CRDcontrac_loc)) CALL dealloc_NParray(CRDcontrac_loc,'CRDcontrac',name_sub)

  !-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'END ',name_sub
  END IF
  !-----------------------------------------------------------

END SUBROUTINE calc_CRD
SUBROUTINE dealloc_tab_RD(para_RD)

  TYPE (param_RD),     intent(inout), allocatable :: para_RD(:)

  integer :: i

  character (len=*), parameter :: name_sub='dealloc_tab_RD'

  IF (allocated(para_RD)) THEN
    DO i=1,size(para_RD)
       CALL dealloc_RD(para_RD(i))
    END DO
    deallocate(para_RD)
  END IF


END SUBROUTINE dealloc_tab_RD
SUBROUTINE tab_RD2_TO_RD1(para_RD1,para_RD2)

TYPE (param_RD), intent(inout), allocatable :: para_RD1(:)
TYPE (param_RD), intent(in)   , allocatable :: para_RD2(:)

integer :: i

character (len=*), parameter :: name_sub='tab_RD2_TO_RD1'

    CALL dealloc_tab_RD(para_RD1)

    IF (allocated(para_RD2)) THEN
      allocate(para_RD1(size(para_RD2)))
      DO i=1,size(para_RD2)
        para_RD1(i) = para_RD2(i)
      END DO
    END IF

END SUBROUTINE tab_RD2_TO_RD1

END MODULE mod_param_RD
