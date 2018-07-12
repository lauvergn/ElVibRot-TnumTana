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

      MODULE mod_psi_SimpleOp
      USE mod_system
      USE mod_psi_set_alloc
      IMPLICIT NONE

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


      psi%norme = ZERO

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


      psi%norme = ZERO


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
                psi_time_R%CvecG = psi%CvecG * cmplx(R,kind=Rkind)
              ELSE
                psi_time_R%RvecG = psi%RvecG * R
              END IF
            END IF

            IF (psi%BasisRep) THEN
              IF (psi%cplx) THEN
                psi_time_R%CvecB = psi%CvecB * cmplx(R,kind=Rkind)
              ELSE
                psi_time_R%RvecB = psi%RvecB * R
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
          FUNCTION psi_time_C(psi,C)
            TYPE (param_psi), intent (in) :: psi
            complex (kind=Rkind), intent (in) :: C
            TYPE (param_psi)  :: psi_time_C
            integer           :: err,i

!           write(out_unitp,*) 'BEGINNING psi_time_C'

!           - define and allocate psi_time_C ----
            CALL copy_psi2TOpsi1(psi_time_C,psi)
            psi_time_C%builtINsub = .TRUE.
!           -----------------------------------------

            IF (psi%BasisRep) THEN
              IF (psi%cplx) THEN
                psi_time_C%CvecB = psi%CvecB * C
              ELSE
                write(out_unitp,*) ' ERROR : in psi_time_C'
                write(out_unitp,*) ' I cannot multiply a real psi and a complex'
                write(out_unitp,*) 'psi%cplx,C',psi%cplx,C
                STOP
              END IF
            END IF

            IF (psi%GridRep) THEN
              IF (psi%cplx) THEN
                psi_time_C%CvecG = psi%CvecG * C
              ELSE
                write(out_unitp,*) ' ERROR : in psi_time_C'
                write(out_unitp,*) ' I cannot multiply a real psi and a complex'
                write(out_unitp,*) 'psi%cplx,C',psi%cplx,C
                STOP
              END IF
            END IF

            psi_time_C%symab = psi%symab

            IF (psi%builtINsub) CALL dealloc_psi(psi)
!           write(out_unitp,*) 'END psi_time_C'

          END FUNCTION psi_time_C

      END MODULE mod_psi_SimpleOp

