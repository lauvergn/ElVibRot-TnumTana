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
!      with contributions of:
!          Mamadou Ndong:       (Tana)
!          Josep Maria Luis:    geometry optimization (ElVibRot)
!          Ahai Chen:           MPI (ElVibRot)
!          Emil Lund klinting:  coupling with MidasCpp (Tana)
!          Lucien Dupuy:        CRP (ElVibRot)
!
!    ElVibRot includes:
!        - Tnum-Tana under the GNU LGPL3 license
!        - Somme subroutines of John Burkardt under GNU LGPL license
!             http://people.sc.fsu.edu/~jburkardt/
!        - Somme subroutines of SHTOOLS written by Mark A. Wieczorek under BSD license
!             http://shtools.ipgp.fr
!        - Some subroutine of QMRPack (see cpyrit.doc) Roland W. Freund and Noel M. Nachtigal:
!             https://www.netlib.org/linalg/qmr/
!===========================================================================
!===========================================================================
      SUBROUTINE versionEVRT(write_version)
      USE mod_system
      USE mod_MPI
      IMPLICIT NONE

      logical :: write_version

      character (len=*), parameter :: EVR_name='ElVibRot'
      character (len=*), parameter :: Tnum_name='Tnum'
      character (len=*), parameter :: Tana_name='Tana'



      IF (write_version .AND. MPI_id==0) THEN
        write(out_unitp,*) '==============================================='
        write(out_unitp,*) '==============================================='
        write(out_unitp,*) 'Working with ',                             &
                   EVR_name,trim(adjustl(EVR_version)),'-',             &
                   Tnum_name,trim(adjustl(Tnum_version)),'-',           &
                   Tana_name,trim(adjustl(Tana_version))

        write(out_unitp,*) 'Compiled on "',trim(compile_host), '" the ',trim(compile_date)
        write(out_unitp,*) 'Compiler version: ',trim(compiler_ver)
        write(out_unitp,*) 'Compiler options: ',trim(compiler_opt)
        write(out_unitp,*) 'Compiler libs: ',trim(compiler_libs)

        write(out_unitp,*) 'EVRT_path: ',trim(EVRT_path)

        write(out_unitp,*) '-----------------------------------------------'

        write(out_unitp,*) EVR_name,' is written by David Lauvergnat [1] '
        write(out_unitp,*) '  with contributions of'
        write(out_unitp,*) '     Josep Maria Luis (optimization) [2]'
        write(out_unitp,*) '     Ahai Chen (MPI) [1,4]'
        write(out_unitp,*) '     Lucien Dupuy (CRP) [5]'

        write(out_unitp,*) EVR_name,' is under GNU LGPL3 license.'
        write(out_unitp,*)

        write(out_unitp,*) Tnum_name,' is written David Lauvergnat [1]'
        write(out_unitp,*) Tana_name,' is written by Mamadou Ndong [1] and David Lauvergnat [1]'
        write(out_unitp,*) '  with contributions'
        write(out_unitp,*) '      Emil Lund klinting (coupling with MidasCpp) [3]'

        write(out_unitp,*) Tnum_name,' and ',Tana_name,' are under GNU LGPL3 license.'
        write(out_unitp,*)
        write(out_unitp,*) '[1]: Laboratoire de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France'
        write(out_unitp,*) '[2]: Institut de Química Computacional and Departament de Química',&
                                   ' Universitat de Girona, Catalonia, Spain'
        write(out_unitp,*) '[3]: Department of Chemistry, Aarhus University, DK-8000 Aarhus C, Denmark'
        write(out_unitp,*) '[4]: Maison de la Simulation USR 3441, CEA Saclay, France'
        write(out_unitp,*) '[5]: Laboratoire Univers et Particule de Montpellier, UMR 5299,', &
                                   ' Université de Montpellier, France'
        write(out_unitp,*) '==============================================='
        write(out_unitp,*) '==============================================='
      END IF
      END SUBROUTINE versionEVRT

