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
      MODULE mod_psi_SimpleOp
      USE mod_system
      USE mod_psi_set_alloc
      IMPLICIT NONE

        !> @todo operator reload for psi operation 
        INTERFACE Set_psi_With_index
          MODULE PROCEDURE Set_psi_With_index_R
          MODULE PROCEDURE Set_psi_With_index_C
        END INTERFACE

        !!@description: TODO
        INTERFACE operator (+)
          MODULE PROCEDURE psi1_plus_psi2
          MODULE PROCEDURE R_plus_psi
          MODULE PROCEDURE psi_plus_R
          MODULE PROCEDURE C_plus_psi
          MODULE PROCEDURE psi_plus_C
        END INTERFACE

        !!@description: TODO
        INTERFACE operator (-)
          MODULE PROCEDURE psi1_minus_psi2
          MODULE PROCEDURE R_minus_psi
          MODULE PROCEDURE psi_minus_R
          MODULE PROCEDURE C_minus_psi
          MODULE PROCEDURE psi_minus_C
        END INTERFACE

        !!@description: TODO
        INTERFACE operator (*)
          MODULE PROCEDURE R_time_psi
          MODULE PROCEDURE psi_time_R
          MODULE PROCEDURE C_time_psi
          MODULE PROCEDURE psi_time_C
        END INTERFACE

        !!@description: TODO
        INTERFACE assignment (=)
          MODULE PROCEDURE R_TOpsi
          MODULE PROCEDURE C_TOpsi
        END INTERFACE

        CONTAINS

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
!     psi = R
!     psi = C
!
!================================================================
      SUBROUTINE R_TOpsi(psi,R)

!----- variables for the WP propagation ----------------------------
      TYPE (param_psi),intent(inout) :: psi
      real (kind=Rkind),intent(in)       :: R

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
      TYPE (param_psi),intent(inout) :: psi
      complex (kind=Rkind),intent(in)    :: C

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


      psi%norm2 = ZERO


      IF (abs(C) == ZERO) THEN
        psi%symab = -2
      ELSE
        psi%symab = -1
      END IF

!     write(out_unitp,*) 'END C_TOpsi'

      END SUBROUTINE C_TOpsi
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
            TYPE (param_psi), intent (in) :: psi1,psi2
            TYPE (param_psi) :: psi1_plus_psi2
            integer           :: err,i

!           write(out_unitp,*) 'BEGINNING psi1_plus_psi2'


!           - define and allocate psi1_plus_psi2 ----
            CALL copy_psi2TOpsi1(psi1_plus_psi2,psi1)
            psi1_plus_psi2%builtINsub = .TRUE.
!           -----------------------------------------


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
                IF(Srep_MPI) psi1_plus_psi2%RS_G=psi1%RS_G+psi2%RS_G
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

            IF (psi1%builtINsub) CALL dealloc_psi(psi1)
            IF (psi2%builtINsub) CALL dealloc_psi(psi2)

!           write(out_unitp,*) 'END psi1_plus_psi2'

          END FUNCTION psi1_plus_psi2


      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
          FUNCTION R_plus_psi(R,psi)
            TYPE (param_psi), intent (in) :: psi
            real (kind=Rkind),    intent (in) :: R
            TYPE (param_psi) :: R_plus_psi
            integer           :: err,i


!           write(out_unitp,*) 'BEGINNING R_plus_psi'

!           - define and allocate R_plus_psi ----
            CALL copy_psi2TOpsi1(R_plus_psi,psi)
            R_plus_psi%builtINsub = .TRUE.
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

            IF (psi%builtINsub) CALL dealloc_psi(psi)
!           write(out_unitp,*) 'END R_plus_psi'

          END FUNCTION R_plus_psi

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
          FUNCTION psi_plus_R(psi,R)
            TYPE (param_psi), intent (in) :: psi
            real (kind=Rkind),    intent (in) :: R
            TYPE (param_psi)  :: psi_plus_R
            integer           :: err,i

!           write(out_unitp,*) 'BEGINNING psi_plus_R'

!           - define and allocate psi_plus_R ----
            CALL copy_psi2TOpsi1(psi_plus_R,psi)
            psi_plus_R%builtINsub = .TRUE.
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

            IF (psi%builtINsub) CALL dealloc_psi(psi)
!           write(out_unitp,*) 'END psi_plus_R'

          END FUNCTION psi_plus_R

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
          FUNCTION C_plus_psi(C,psi)
            TYPE (param_psi), intent (in) :: psi
            complex (kind=Rkind), intent (in) :: C
            TYPE (param_psi)  :: C_plus_psi
            integer           :: err,i

!           write(out_unitp,*) 'BEGINNING C_plus_psi'

!           - define and allocate C_plus_psi ----
            CALL copy_psi2TOpsi1(C_plus_psi,psi)
            C_plus_psi%builtINsub = .TRUE.
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

            IF (psi%builtINsub) CALL dealloc_psi(psi)
!           write(out_unitp,*) 'END C_plus_psi'

          END FUNCTION C_plus_psi

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
          FUNCTION psi_plus_C(psi,C)
            TYPE (param_psi), intent (in) :: psi
            complex (kind=Rkind), intent (in) :: C
            TYPE (param_psi)  :: psi_plus_C
            integer           :: err,i

!           write(out_unitp,*) 'BEGINNING psi_plus_C'

!           - define and allocate psi_plus_C ----
            CALL copy_psi2TOpsi1(psi_plus_C,psi)
            psi_plus_C%builtINsub = .TRUE.
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

            IF (psi%builtINsub) CALL dealloc_psi(psi)
!           write(out_unitp,*) 'END psi_plus_C'

          END FUNCTION psi_plus_C
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
            TYPE (param_psi), intent (in) :: psi1,psi2
            TYPE (param_psi) :: psi1_minus_psi2
            integer          :: err,i


!           write(out_unitp,*) 'BEGINNING psi1_minus_psi2'

!           - define and allocate psi1_minus_psi2 ----
            CALL copy_psi2TOpsi1(psi1_minus_psi2,psi1)
            psi1_minus_psi2%builtINsub = .TRUE.
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
                IF(Srep_MPI) psi1_plus_psi2%RS_G=psi1%RS_G-psi2%RS_G
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


            IF (psi1%builtINsub) CALL dealloc_psi(psi1)
            IF (psi2%builtINsub) CALL dealloc_psi(psi2)

!           write(out_unitp,*) 'END psi1_minus_psi2'

          END FUNCTION psi1_minus_psi2

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
          FUNCTION R_minus_psi(R,psi)
            TYPE (param_psi), intent (in) :: psi
            real (kind=Rkind),    intent (in) :: R
            TYPE (param_psi) :: R_minus_psi
            integer           :: err,i


!           write(out_unitp,*) 'BEGINNING R_minus_psi'

!           - define and allocate R_minus_psi ----
            CALL copy_psi2TOpsi1(R_minus_psi,psi)
            R_minus_psi%builtINsub = .TRUE.
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


            IF (psi%builtINsub) CALL dealloc_psi(psi)
!           write(out_unitp,*) 'END R_minus_psi'

          END FUNCTION R_minus_psi

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
          FUNCTION psi_minus_R(psi,R)
            TYPE (param_psi), intent (in) :: psi
            real (kind=Rkind),    intent (in) :: R
            TYPE (param_psi)  :: psi_minus_R
            integer           :: err,i

!           write(out_unitp,*) 'BEGINNING psi_minus_R'

!           - define and allocate psi_minus_R ----
            CALL copy_psi2TOpsi1(psi_minus_R,psi)
            psi_minus_R%builtINsub = .TRUE.
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


            IF (psi%builtINsub) CALL dealloc_psi(psi)
!           write(out_unitp,*) 'END psi_minus_R'

          END FUNCTION psi_minus_R

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
          FUNCTION C_minus_psi(C,psi)
            TYPE (param_psi), intent (in) :: psi
            complex (kind=Rkind), intent (in) :: C
            TYPE (param_psi)  :: C_minus_psi
            integer           :: err,i

!           write(out_unitp,*) 'BEGINNING C_minus_psi'

!           - define and allocate C_minus_psi ----
            CALL copy_psi2TOpsi1(C_minus_psi,psi)
            C_minus_psi%builtINsub = .TRUE.
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

            IF (psi%builtINsub) CALL dealloc_psi(psi)
!           write(out_unitp,*) 'END C_minus_psi'

          END FUNCTION C_minus_psi

          !!@description: TODO
          FUNCTION psi_minus_C(psi,C)
            TYPE (param_psi), intent (in) :: psi
            complex (kind=Rkind), intent (in) :: C
            TYPE (param_psi)  :: psi_minus_C
            integer           :: err,i

!           write(out_unitp,*) 'BEGINNING psi_minus_C'

!           - define and allocate psi_minus_C ----
            CALL copy_psi2TOpsi1(psi_minus_C,psi)
            psi_minus_C%builtINsub = .TRUE.
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

            IF (psi%builtINsub) CALL dealloc_psi(psi)
!           write(out_unitp,*) 'END psi_minus_C'

          END FUNCTION psi_minus_C
!================================================================
!
!         psi*R ....
!
!================================================================
          !!@description: TODO
          FUNCTION R_time_psi(R,psi)
            TYPE (param_psi), intent (in) :: psi
            real (kind=Rkind),    intent (in) :: R
            TYPE (param_psi) :: R_time_psi
            integer           :: err,i

!           write(out_unitp,*) 'BEGINNING R_time_psi'

!           - define and allocate R_time_psi ----
            CALL copy_psi2TOpsi1(R_time_psi,psi)
            R_time_psi%builtINsub = .TRUE.
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

            R_time_psi%symab = psi%symab

            IF (psi%builtINsub) CALL dealloc_psi(psi)
!           write(out_unitp,*) 'END R_time_psi'

          END FUNCTION R_time_psi

          !!@description: TODO
          FUNCTION psi_time_R(psi,R)
            USE mod_MPI

            TYPE (param_psi), intent (in) :: psi
            real (kind=Rkind),    intent (in) :: R
            TYPE (param_psi)  :: psi_time_R
            integer           :: err,i

!           write(out_unitp,*) 'BEGINNING psi_time_R'

!           - define and allocate psi_time_R ----
            CALL copy_psi2TOpsi1(psi_time_R,psi)
            psi_time_R%builtINsub = .TRUE.
!           -----------------------------------------


            IF (psi%GridRep) THEN
              IF (psi%cplx) THEN
                IF(MPI_id==0) psi_time_R%CvecG = psi%CvecG * cmplx(R,kind=Rkind)
              ELSE
                IF(MPI_id==0) psi_time_R%RvecG = psi%RvecG * R
              END IF
            END IF

            IF (psi%BasisRep) THEN
              IF (psi%cplx) THEN
                IF(MPI_id==0) psi_time_R%CvecB = psi%CvecB * cmplx(R,kind=Rkind)
              ELSE
                IF(MPI_id==0) psi_time_R%RvecB = psi%RvecB * R
              END IF
            END IF

            psi_time_R%symab = psi%symab

            IF (psi%builtINsub) CALL dealloc_psi(psi)
!           write(out_unitp,*) 'END psi_time_R'

          END FUNCTION psi_time_R

          !!@description: TODO
          FUNCTION C_time_psi(C,psi)
            TYPE (param_psi), intent (in) :: psi
            complex (kind=Rkind), intent (in) :: C
            TYPE (param_psi)  :: C_time_psi
            integer           :: err,i

!           write(out_unitp,*) 'BEGINNING C_time_psi'

!           - define and allocate C_time_psi ----
            CALL copy_psi2TOpsi1(C_time_psi,psi)
            C_time_psi%builtINsub = .TRUE.
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

            IF (psi%builtINsub) CALL dealloc_psi(psi)
!           write(out_unitp,*) 'END C_time_psi'

          END FUNCTION C_time_psi

          !!@description: TODO
!=======================================================================================
          FUNCTION psi_time_C(psi,C)
            TYPE (param_psi), intent (in) :: psi
            complex (kind=Rkind), intent (in) :: C
            TYPE (param_psi)  :: psi_time_C
            integer           :: err,i

            !write(out_unitp,*) 'BEGINNING psi_time_C'
            !CALL flush_perso(out_unitp)

!           - define and allocate psi_time_C ----
            CALL copy_psi2TOpsi1(psi_time_C,psi)
            psi_time_C%builtINsub = .TRUE.
!           -----------------------------------------

            IF (psi%BasisRep) THEN
              IF (psi%cplx) THEN
                IF(MPI_id==0) psi_time_C%CvecB = psi%CvecB * C
              ELSE
                write(out_unitp,*) ' ERROR : in psi_time_C'
                write(out_unitp,*) ' I cannot multiply a real psi and a complex'
                write(out_unitp,*) 'psi%cplx,C',psi%cplx,C
                STOP
              END IF
            END IF

            IF (psi%GridRep) THEN
              IF (psi%cplx) THEN
                IF(MPI_id==0) psi_time_C%CvecG = psi%CvecG * C
              ELSE
                write(out_unitp,*) ' ERROR : in psi_time_C'
                write(out_unitp,*) ' I cannot multiply a real psi and a complex'
                write(out_unitp,*) 'psi%cplx,C',psi%cplx,C
                STOP
              END IF
            END IF

            psi_time_C%symab = psi%symab

            IF (psi%builtINsub) CALL dealloc_psi(psi)

            !write(out_unitp,*) 'END psi_time_C'
            !CALL flush_perso(out_unitp)

          END FUNCTION psi_time_C
!=======================================================================================

      END MODULE mod_psi_SimpleOp

