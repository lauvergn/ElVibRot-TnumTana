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
      MODULE mod_ComOp

      USE mod_system
      USE mod_basis_BtoG_GtoB_SGType4, only : Type_SmolyakRep, dealloc_SmolyakRep

      IMPLICIT NONE

      PRIVATE

      TYPE param_ComOp

          real (kind=Rkind), allocatable:: sqRhoOVERJac(:)            ! sqrt(rho/jac), nb_qa pts
          real (kind=Rkind), allocatable:: Jac(:)                     ! jac, nb_qa pts
          TYPE(Type_SmolyakRep)         :: SRep_sqRhoOVERJac          ! sqrt(rho/jac) in Smolyap Rep.
          TYPE(Type_SmolyakRep)         :: SRep_Jac                   ! jac in Smolyap Rep.

      CONTAINS
          PROCEDURE, PRIVATE, PASS(ComOp1) :: ComOp2_TO_ComOp1
          GENERIC,   PUBLIC  :: assignment(=) => ComOp2_TO_ComOp1
      END TYPE param_ComOp
      !====================================================================

      PUBLIC :: param_ComOp,dealloc_ComOp,write_param_ComOp

      CONTAINS

      SUBROUTINE dealloc_ComOp(ComOp,keep_init)

      TYPE (param_ComOp),          intent(inout) :: ComOp
      logical,           optional, intent(in)    :: keep_init

      logical :: keep_init_loc

      character (len=*), parameter :: name_sub='dealloc_ComOp'

      keep_init_loc = .FALSE.
      IF (present(keep_init)) keep_init_loc = keep_init

      IF (allocated(ComOp%sqRhoOVERJac)) THEN
        CALL dealloc_NParray(ComOp%sqRhoOVERJac,"ComOp%sqRhoOVERJac",name_sub)
      END IF
      IF (allocated(ComOp%Jac)) THEN
        CALL dealloc_NParray(ComOp%Jac,"ComOp%Jac",name_sub)
      END IF
      CALL dealloc_SmolyakRep(ComOp%SRep_sqRhoOVERJac)
      CALL dealloc_SmolyakRep(ComOp%SRep_Jac)

      END SUBROUTINE dealloc_ComOp

      SUBROUTINE ComOp2_TO_ComOp1(ComOp1,ComOp2)

      TYPE (param_ComOp),           intent(in)    :: ComOp2
      CLASS (param_ComOp),          intent(inout) :: ComOp1

      character (len=*), parameter :: name_sub='ComOp2_TO_ComOp1'

      IF (allocated(ComOp2%sqRhoOVERJac)) THEN
        CALL alloc_NParray(ComOp1%sqRhoOVERJac,shape(ComOp2%sqRhoOVERJac),&
                          "ComOp1%sqRhoOVERJac",name_sub)
        ComOp1%sqRhoOVERJac(:) = ComOp2%sqRhoOVERJac
      END IF
      IF (allocated(ComOp2%Jac)) THEN
        CALL alloc_NParray(ComOp1%Jac,shape(ComOp2%Jac),                &
                          "ComOp1%Jac",name_sub)
        ComOp1%Jac(:) = ComOp2%Jac
      END IF

      ComOp1%SRep_sqRhoOVERJac = ComOp2%SRep_sqRhoOVERJac
      ComOp1%SRep_Jac          = ComOp2%SRep_Jac

      END SUBROUTINE ComOp2_TO_ComOp1

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

      IF (allocated(ComOp%sqRhoOVERJac)) THEN
        write(out_unitp,*) 'shape sqRhoOVERJac',shape(ComOp%sqRhoOVERJac)
      END IF
      IF (allocated(ComOp%Jac)) THEN
        write(out_unitp,*) 'shape Jac ',shape(ComOp%Jac)
      END IF

      write(out_unitp,*) 'END WRITE param_ComOp'


      END SUBROUTINE write_param_ComOp

      END MODULE mod_ComOp

