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
MODULE mod_basis_RCVec_SGType4
USE mod_system
IMPLICIT NONE

PRIVATE

TYPE TypeRVec
  real(kind=Rkind), allocatable :: V(:)
END TYPE TypeRVec
TYPE TypeCVec
  complex(kind=Rkind), allocatable :: V(:)
END TYPE TypeCVec

INTERFACE assignment(=)
  module procedure TypeRVec2_TO_TypeRVec1,tabR2_TO_TypeRVec1
END INTERFACE

PUBLIC  TypeRVec, alloc_TypeRVec, dealloc_TypeRVec, Write_TypeRVec, &
        TypeRVec2_TO_TypeRVec1, tabR2_TO_TypeRVec1,                 &
        sub_ReadRVec, sub_WriteRVec

PUBLIC  TypeCVec, alloc_TypeCVec, dealloc_TypeCVec, Write_TypeCVec, &
        TypeCVec2_TO_TypeCVec1, tabC2_TO_TypeCVec1

CONTAINS

SUBROUTINE alloc_TypeRVec(Rvec,nvec)
USE mod_system
IMPLICIT NONE

  TYPE (TypeRVec),    intent(inout) :: Rvec
  integer,            intent(in)    :: nvec


  CALL dealloc_TypeRVec(Rvec)

  IF (nvec < 1) THEN
    write(6,*) ' ERROR in alloc_TypeRVec'
    write(6,*) ' nvec < 1',nvec
    STOP
  END IF


  CALL alloc_NParray(Rvec%V,(/nvec/),'Rvec%V','alloc_TypeRVec')

END SUBROUTINE alloc_TypeRVec
SUBROUTINE dealloc_TypeRVec(Rvec)
USE mod_system
IMPLICIT NONE

  TYPE (TypeRVec), intent(inout) :: Rvec

  IF (allocated(Rvec%V)) THEN
    CALL dealloc_NParray(Rvec%V,'Rvec%V','dealloc_TypeRVec')
  END IF

END SUBROUTINE dealloc_TypeRVec

SUBROUTINE Write_TypeRVec(Rvec)
USE mod_system
IMPLICIT NONE

  TYPE (TypeRVec), intent(in) :: Rvec

  IF (allocated(Rvec%V)) THEN
    CALL Write_VecMat(Rvec%V,out_unitp,5,name_info='R:')
  END IF

END SUBROUTINE Write_TypeRVec

SUBROUTINE TypeRVec2_TO_TypeRVec1(Rvec1,Rvec2)
USE mod_system
IMPLICIT NONE

  TYPE (TypeRVec), intent(inout) :: Rvec1
  TYPE (TypeRVec), intent(in)    :: Rvec2

  IF (allocated(Rvec2%V)) THEN
    CALL alloc_TypeRVec(Rvec1,nvec=size(Rvec2%V))
    Rvec1%V(:) = Rvec2%V
  ELSE
    CALL dealloc_TypeRVec(Rvec1)
  END IF

END SUBROUTINE TypeRVec2_TO_TypeRVec1
SUBROUTINE tabR2_TO_TypeRVec1(Rvec1,tabR2)
USE mod_system
IMPLICIT NONE

  TYPE (TypeRVec),                intent(inout) :: Rvec1
  real(kind=Rkind), allocatable,  intent(in)    :: tabR2(:)

  IF (allocated(tabR2)) THEN
    CALL alloc_TypeRVec(Rvec1,nvec=size(tabR2))
    Rvec1%V(:) = tabR2(:)
  ELSE
    CALL dealloc_TypeRVec(Rvec1)
  END IF

END SUBROUTINE tabR2_TO_TypeRVec1

SUBROUTINE sub_ReadRVec(RVec,FileName_RVec,err_sub)
  USE mod_system
  IMPLICIT NONE

  TYPE (TypeRVec),                 intent(inout)           :: RVec
  character (len=Line_len),        intent(in)              :: FileName_RVec
  integer,                         intent(inout), optional :: err_sub

  character (len=*),               parameter              :: Name_sub='sub_ReadRVec'


  integer :: nio,error,err_file,nvec

  error = 0
  nvec  = 0

  CALL file_open2(FileName_RVec,nio,lformatted=.FALSE.,old=.TRUE.,err_file=err_file)
  IF (err_file /= 0) THEN
    write(out_unitp,*) ' ERROR in ',Name_sub
    write(out_unitp,*) '   Problem with the file associated to RVec'
    write(out_unitp,*) '   err_file: ',err_file
    error = 2
  ELSE
    read(nio,iostat=err_file) nvec
    IF (err_file /= 0 .OR. nvec < 1) THEN
      write(out_unitp,*) ' ERROR in ',Name_sub
      write(out_unitp,*) '   Error while reading nvec',nvec
      error = 3
      nvec = 0
    END IF
  END IF

  IF (error == 0) THEN
    CALL alloc_TypeRVec(Rvec,nvec)
    read(nio,iostat=err_file) Rvec%V
    IF (err_file /= 0) THEN
      write(out_unitp,*) ' ERROR in ',Name_sub
      write(out_unitp,*) '   Error while reading Rvec'
      error = 4
    END IF
  END IF

  close(nio,iostat=err_file)

  IF (present(err_sub)) THEN
    err_sub = error
  ELSE
    STOP ' in sub_ReadRVec'
  END IF

END SUBROUTINE sub_ReadRVec
SUBROUTINE sub_WriteRVec(RVec,FileName_RVec,err_sub)
  USE mod_system
  IMPLICIT NONE

  TYPE (TypeRVec),                 intent(in)              :: RVec
  character (len=Line_len),        intent(in)              :: FileName_RVec
  integer,                         intent(inout), optional :: err_sub

  character (len=*),               parameter              :: Name_sub='sub_WriteRVec'


  integer :: nio,error,err_file,nvec

  error = 0
  nvec  = size(Rvec%V)

  IF (nvec < 1) THEN
    error = 1
  ELSE
    CALL file_open2(FileName_RVec,nio,lformatted=.FALSE.,old=.FALSE.,err_file=err_file)
    IF (err_file /= 0) THEN
      write(out_unitp,*) ' ERROR in ',Name_sub
      write(out_unitp,*) '   Problem with the file associated to RVec'
      write(out_unitp,*) '   err_file: ',err_file
      error = 2
    ELSE
      write(nio,iostat=err_file) nvec
      IF (err_file /= 0 .OR. nvec < 1) THEN
        write(out_unitp,*) ' ERROR in ',Name_sub
        write(out_unitp,*) '   Error while writing nvec',nvec
        error = 3
        nvec = 0
      ELSE
        write(nio,iostat=err_file) Rvec%V
        IF (err_file /= 0) THEN
          write(out_unitp,*) ' ERROR in ',Name_sub
          write(out_unitp,*) '   Error while reading Rvec'
          error = 4
        END IF
      END IF
    END IF
  END IF

  close(nio,iostat=err_file)

  IF (present(err_sub)) THEN
    err_sub = error
  ELSE
    STOP ' in sub_WriteRVec'
  END IF

END SUBROUTINE sub_WriteRVec
SUBROUTINE alloc_TypeCVec(vec,nvec)
USE mod_system
IMPLICIT NONE

  TYPE (TypeCVec), intent(inout) :: vec
  integer,         intent(in)    :: nvec


  CALL dealloc_TypeCVec(vec)

  IF (nvec < 1) THEN
    write(6,*) ' ERROR in alloc_TypeCVec'
    write(6,*) ' nvec < 1',nvec
    STOP
  END IF

  CALL alloc_NParray(vec%V,(/nvec/),'vec%V','alloc_TypeCVec')

END SUBROUTINE alloc_TypeCVec
SUBROUTINE dealloc_TypeCVec(vec)
USE mod_system
IMPLICIT NONE

  TYPE (TypeCVec), intent(inout) :: vec

  IF (allocated(vec%V)) THEN
    CALL dealloc_NParray(vec%V,'vec%V','dealloc_TypeCVec')
  END IF

END SUBROUTINE dealloc_TypeCVec

SUBROUTINE Write_TypeCVec(vec)
USE mod_system
IMPLICIT NONE

  TYPE (TypeCVec), intent(in) :: vec

  IF (allocated(vec%V)) THEN
    CALL Write_VecMat(vec%V,out_unitp,5,name_info='R:')
  END IF

END SUBROUTINE Write_TypeCVec

SUBROUTINE TypeCVec2_TO_TypeCVec1(vec1,vec2)
USE mod_system
IMPLICIT NONE

  TYPE (TypeCVec), intent(inout) :: vec1
  TYPE (TypeCVec), intent(in)    :: vec2

  IF (allocated(vec2%V)) THEN
    CALL alloc_TypeCVec(vec1,nvec=size(vec2%V))
    vec1%V(:) = vec2%V
  ELSE
    CALL dealloc_TypeCVec(vec1)
  END IF

END SUBROUTINE TypeCVec2_TO_TypeCVec1
SUBROUTINE tabC2_TO_TypeCVec1(vec1,tab2)
USE mod_system
IMPLICIT NONE

  TYPE (TypeCVec),                   intent(inout) :: vec1
  complex(kind=Rkind), allocatable,  intent(in)    :: tab2(:)

  IF (allocated(tab2)) THEN
    CALL alloc_TypeCVec(vec1,nvec=size(tab2))
    vec1%V(:) = tab2(:)
  ELSE
    CALL dealloc_TypeCVec(vec1)
  END IF

END SUBROUTINE tabC2_TO_TypeCVec1

END MODULE mod_basis_RCVec_SGType4
