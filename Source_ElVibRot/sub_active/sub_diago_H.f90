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
      SUBROUTINE sub_diago_H(H,E,Vec,n,sym)
      USE mod_system
      USE mod_Constant, ONLY: get_Conv_au_TO_unit
      IMPLICIT NONE

!------ active Matrix H Vec E ------------------------------------

      integer           :: n
      real (kind=Rkind) :: H(n,n)
      real (kind=Rkind) :: Vec(n,n)
      real (kind=Rkind) :: E(n)
      logical           :: sym

!----- divers --------------------------------------------------------
      real (kind=Rkind) :: auTOcm_inv
      real (kind=Rkind) :: A,B

      real (kind=Rkind), allocatable :: Hsave(:,:),r(:)
      logical           :: residual = .FALSE.

      integer :: i


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_diago_H'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'n',n
         flush(out_unitp)
       END IF
!---------------------------------------------------------------------------------------
!       H matrix diagonalisation
!---------------------------------------------------------------------------------------
  auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')

  IF(keep_MPI) THEN
    IF (residual) THEN
      CALL alloc_NParray(Hsave,[n,n],'Hsave',name_sub)
      CALL alloc_NParray(r,[n],'r',name_sub)
      Hsave(:,:) = H(:,:)
    END IF
  ENDIF

  IF (.NOT. sym) THEN
    IF(keep_MPI) CALL diagonalization(H,E,Vec,n,4,1,.TRUE.)
  ELSE
    IF(keep_MPI) CALL diagonalization(H,E,Vec,n,3,1,.TRUE.)
  END IF

  IF (debug) THEN
    write(out_unitp,*) 'OpMin,OpMax (ua)  : ',minval(E),maxval(E)
    write(out_unitp,*) 'OpMin,OpMax (cm-1): ',minval(E)*auTOcm_inv,maxval(E)*auTOcm_inv
  END IF

  IF(keep_MPI) THEN
    DO i=1,n
      A = maxval(Vec(:,i))
      B = minval(Vec(:,i))
      IF (abs(A) < abs(B)) Vec(:,i) = -Vec(:,i)
    END DO
  ENDIF
  
  IF(keep_MPI) THEN
    IF (residual) THEN
      DO i=1,n
        r(:) = matmul(Hsave,Vec(:,i))-E(i)*Vec(:,i)
        write(out_unitp,*) 'residual (cm-1)',i,sqrt(dot_product(r,r))*&
                                                         auTOcm_inv
        write(out_unitp,*) 'residual (au)',i,sqrt(dot_product(r,r))
        write(out_unitp,*) 'err (cm-1)',i,dot_product(Vec(:,i),r)*    &
                                                         auTOcm_inv
        write(out_unitp,*) 'err (au)',i,dot_product(Vec(:,i),r)
      END DO
      CALL dealloc_NParray(Hsave,'Hsave',name_sub)
      CALL dealloc_NParray(r,'r',name_sub)
    END IF
  END IF !for keep_MPI!

!------------------------------------------------------
  IF (debug) THEN

    write(out_unitp,*) ' level energy :'
    DO i=1,n
      write(out_unitp,*) i,E(i)*auTOcm_inv,(E(i)-E(1))*auTOcm_inv
    END DO

    write(out_unitp,*) ' Vec:'
    CALL Write_Mat(Vec,out_unitp,5)

    write(out_unitp,*) 'END ',name_sub
    flush(out_unitp)
  END IF
!------------------------------------------------------

END SUBROUTINE sub_diago_H
!=======================================================================================

!
!=====================================================================
!
!     diagonalisation complexe
!
!=====================================================================
      SUBROUTINE sub_diago_CH(CH,CE,CVec,n)
      USE mod_system
      USE mod_Constant, ONLY: get_Conv_au_TO_unit
      IMPLICIT NONE
      !
      !------ active Matrix H Vec E ------------------------------------
      integer              :: n
      complex (kind=Rkind) :: CH(n,n)
      complex (kind=Rkind) :: CVec(n,n)
      complex (kind=Rkind) :: CE(n)
      complex (kind=Rkind) :: trav(n)

      real (kind=Rkind)    :: auTOcm_inv


!----- divers --------------------------------------------------------
      integer              :: i
      integer              :: ierr


!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING sub_diago_CH'
      END IF
!-----------------------------------------------------------
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')

!=====================================================================
!
!     H matrix diagonalisation
!
!=====================================================================

      CALL cTred2(n,n,CH,CE,trav,CVec)
      CALL cTql2(n,n,CE,trav,CVec,ierr)
      write(out_unitp,*)'ierr=',ierr

!=====================================================================
!
!
!=====================================================================

!-----------------------------------------------------------
      IF (debug) THEN

        write(out_unitp,*) ' level energy :'
        DO i=1,n
          write(out_unitp,*) i,CE(i)*cmplx(auTOcm_inv,kind=Rkind),     &
                             (CE(i)-CE(1))*cmplx(auTOcm_inv,kind=Rkind)
        END DO

        write(out_unitp,*) ' Vec:'
        CALL Write_Mat(CVec,out_unitp,5)

        write(out_unitp,*) 'END sub_diago_CH'
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_diago_CH
