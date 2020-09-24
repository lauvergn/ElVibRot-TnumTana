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
Module mod_Hmax_MPI

Private 
Public get_Hmin_MPI

CONTAINS
  SUBROUTINE get_Hmin_MPI(para_H)
    USE mod_system
    USE mod_SetOp
    USE mod_MPI_aux
    IMPLICIT NONE

    TYPE(param_Op),                   intent(inout) :: para_H

#if(run_MPI) 

    para_H%Hmin=1.0e10
    DO i1_loop=iGs_MPI(1,MPI_id),iGs_MPI(2,MPI_id)
      ! Minimal value of Veff
      temp_real=minval(para_H%OpGrid(1)%SRep%SmolyakRep(i1_loop)%V)
      IF(temp_real<para_H%Hmin) para_H%Hmin=temp_real
    ENDDO

    IF(MPI_np>1) THEN
      CALL MPI_Reduce(para_H%Hmin,temp_real,size1_MPI,Real_MPI,MPI_MIN,                &
                      root_MPI,MPI_COMM_WORLD,MPI_err)
      IF(MPI_id==0) para_H%Hmin=temp_real
      CALL MPI_Bcast(para_H%Hmin,size1_MPI,Real_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
    ENDIF

#endif
  ENDSUBROUTINE

ENDMODULE  mod_Hmax_MPI
