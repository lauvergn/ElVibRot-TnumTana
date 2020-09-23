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

MODULE mod_OpPsi_MPI
  PUBLIC :: sub_scaledOpPsi_SR_MPI

  CONTAINS
!=======================================================================================
!> OpPsi = (OpPsi - E0*Psi) / Esc, for working on full Smolyak rep.
!=======================================================================================
    SUBROUTINE sub_scaledOpPsi_SR_MPI(Psi,OpPsi,E0,Esc)
      USE mod_system
      USE mod_psi,ONLY:param_psi,ecri_psi
      IMPLICIT NONE

      TYPE(param_psi)                          :: OpPsi
      TYPE(param_psi)                          :: Psi
      Real(kind=Rkind)                         :: E0
      Real(kind=Rkind)                         :: Esc

      IF(Psi%SRG_MPI .AND. OpPsi%SRG_MPI) THEN
        OpPsi%SR_G(:,:)=(OpPsi%SR_G(:,:)-E0*Psi%SR_G(:,:))/Esc
      ELSE IF(Psi%SRB_MPI .AND. OpPsi%SRB_MPI) THEN
        OpPsi%SR_B(:,:)=(OpPsi%SR_B(:,:)-E0*Psi%SR_B(:,:))/Esc
      ELSE
        STOP 'error in sub_scaledOpPsi_SR_MPI'
      ENDIF
      
      IF(OpPsi%symab/=Psi%symab) THEN
        IF(OpPsi%symab==-2) THEN
          OpPsi%symab=Psi%symab
        ELSE
          OpPsi%symab=-1
        ENDIF
      ENDIF

    ENDSUBROUTINE sub_scaledOpPsi_SR_MPI
!=======================================================================================

END MODULE mod_OpPsi_MPI
