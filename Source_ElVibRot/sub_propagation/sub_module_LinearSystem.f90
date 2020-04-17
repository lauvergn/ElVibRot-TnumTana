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
MODULE mod_LinearSystem
USE mod_Constant
USE mod_MPI 
IMPLICIT NONE

PRIVATE
PUBLIC :: sub_LinearSystem

CONTAINS
!===============================================================================
!
!  Solve with an iterative method: Op|TabPsi> = |TabOpPsi>
!     => find TabPsi from TabOpPsi
!
!===============================================================================
SUBROUTINE sub_LinearSystem(TabOpPsi,TabPsi,Op,para_propa)
      USE mod_system
      USE mod_psi,    ONLY : param_psi
      USE mod_propa,  ONLY : param_propa,param_Davidson
      USE mod_Op
      USE mod_MPI
      IMPLICIT NONE

      !----- psi ---------------------------------------------
      TYPE (param_psi),        intent(in)    :: OpTabPsi(:)
      TYPE (param_psi),        intent(inout) :: TabPsi(:)

      !----- Operator: Hamiltonian ----------------------------
      TYPE (param_Op),         intent(in)    :: Op

      !----- WP, energy ... -----------------------------------
      TYPE (param_propa),      intent(in)    :: para_propa


      !----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_LinearSystem'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      !-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
        CALL flush_perso(out_unitp)
      END IF
      !-----------------------------------------------------------

      CALL sub_LinearSystem_Jacobi(TabOpPsi,TabPsi,Op,para_propa)

      !----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
      !----------------------------------------------------------

 END SUBROUTINE sub_LinearSystem

 ! solve iteratively Op.psi = OpPsi with:
 !   psi(k) = psi(k-1) + D0^-1(OpPsi - Op.
 SUBROUTINE sub_LinearSystem_Jacobi(TabOpPsi,TabPsi,para_H,Ene,para_propa)
      USE mod_system
      USE mod_psi,   ONLY : param_psi,alloc_psi,alloc_NParray,norm2_psi,&
                            sub_LCpsi_TO_psi,sub_save_psi
      USE mod_Op
      USE mod_propa, ONLY : param_propa,param_Davidson
      USE mod_MPI
      IMPLICIT NONE

      !----- psi ---------------------------------------------
      TYPE (param_psi), intent(in)    :: OpTabPsi(:)
      TYPE (param_psi), intent(inout) :: TabPsi(:)

      !----- Operator: Hamiltonian ----------------------------
      TYPE (param_Op),  intent(in)    :: para_H

      !----- WP, energy ... -----------------------------------
      TYPE (param_propa), intent(in)  :: para_propa




      !------ working parameters --------------------------------
      TYPE(param_file)  :: Log_file
      integer           :: iunit

      TYPE (param_psi) :: HDiagInv


      !----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_LinearSystem'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      !-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
        CALL flush_perso(out_unitp)
      END IF


      !-----------------------------------------------------------
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')
      auTOene    = get_Conv_au_TO_WriteUnit('E',WriteUnit)

      !------ initialization -------------------------------------
      Log_file%name='LinearSystem.log'
      CALL file_open(Log_file,iunit)


      ! first the diagonal part of Ene-H
      CALL init_psi(HDiagInv,para_H,Op%cplx)
      IF (allocated(para_H%BasisnD%EneH0)) THEN
        write(out_unitp,*) 'precon /= 1. DML'
        IF (OpInvDiag_part%cplx) THEN
          HDiagInv%CvecB(:) = ONE/(Ene-para_H%BasisnD%EneH0(:)) ! approximation of 1/(Ene-H(i,i))
        ELSE
          HDiagInv%RvecB(:) = ONE/(Ene-para_H%BasisnD%EneH0(:)) ! approximation of 1/(Ene-H(i,i))
        END IF
      ELSE
        write(out_unitp,*) 'precon = 1. DML'
        STOP 'ERROR in sub_LinearSystem: Op%BasisnD%EneH0 is not allocated!!'
      END IF

      DO it=1,max_it
      END DO

STOP 'not finished!!'
      !----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
      !----------------------------------------------------------

 END SUBROUTINE sub_LinearSystem_Jacobi

END MODULE mod_LinearSystem
