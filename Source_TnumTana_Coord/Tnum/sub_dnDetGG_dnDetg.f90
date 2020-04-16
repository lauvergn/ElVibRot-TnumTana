!===========================================================================
!===========================================================================
!This file is part of Tnum-Tana.
!
!    Tnum-Tana is a free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    Tnum-Tana is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ElVibRot.  If not, see <http://www.gnu.org/licenses/>.
!
!    Copyright 2015  David Lauvergnat
!      Tnum is written David Lauvergnat [1]
!      Tana is written by Mamadou Ndong [1] and David Lauvergnat [1]
!         with contributions
!          Emil Lund klinting (coupling with MidasCpp) [3]'
!
![1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France
![3]: Department of Chemistry, Aarhus University, DK-8000 Aarhus C, Denmark
!
!===========================================================================
!===========================================================================
MODULE mod_dnDetGG_dnDetg
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: sub3_dndetGG,sub3_dndetA

  CONTAINS

!=====================================================================
!
! ++   calculation of  the Jacobian J=sqrt(det(g)) and its log derivatives
!      using pivot algorithmes
!      J               = sqrt(det(g))
!      (dJ/dQi) / J    = 1/2 ( ddet(g)/dQi ) / det(g)
!      (dJ/dQidQj) / J = 1/2 ( ddet(g)/dQidQj ) / det(g) - (dJ/dQi) / J * (dJ/dQj) / J
!=====================================================================
!
      SUBROUTINE sub3_dndetA(dndetA,dnA,nderiv,                         &
                             masses,Mtot_inv,ncart)
      USE mod_system
      USE mod_dnSVM
      IMPLICIT NONE


      integer :: nderiv

      TYPE(Type_dnS)    :: dndetA

      TYPE(Type_dnMat)  ::  dnA
      real (kind=Rkind) ::  d0Aii_inv(dnA%nb_var_Matl)
      real (kind=Rkind) ::  d0Aii_sqinv(dnA%nb_var_Matl)

      integer :: ncart
      real (kind=Rkind) :: mass,masses(ncart),Mtot_inv

      TYPE(Type_dnMat)  :: dnA_save
      integer :: i
      integer :: ider,jder

      integer :: err,memory
      logical :: keep = .FALSE.
!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
!      logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'sub3_dndetA'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nderiv',nderiv

        write(out_unitp,*) 'dnA'
        CALL Write_dnSVM(dnA,nderiv)
      END IF
!-----------------------------------------------------------

!-----------------------------------------------------------
      IF (keep) THEN
        CALL alloc_dnSVM(dnA_save,dnA%nb_var_Matl,dnA%nb_var_Matc,      &
                         dnA%nb_var_deriv,nderiv)
        CALL sub_dnMat1_TO_dnMat2(dnA,dnA_save,nderiv)
      END IF

      dndetA%d0 = ONE
!-----------------------------------------------------------


!-----------------------------------------------------------
!-----------------------------------------------------------
      DO i=1,dnA%nb_var_Matl

        CALL sub3_pivot(i,dndetA%d0,dnA,nderiv)
        CALL sub3_newline(i,dnA,nderiv)

      END DO
!-----------------------------------------------------------
!-----------------------------------------------------------

!-----------------------------------------------------------
!     determinant and calculation of 1./d0A(i,i) 1./d0A(i,i)**2
      DO i=1,dnA%nb_var_Matl
        mass = masses(i)
        IF (mass == ZERO) mass = ONE
!       d0d = d0d*d0A(i,i)/(mass*mass)
        dndetA%d0 = dndetA%d0 * dnA%d0(i,i)/mass
        d0Aii_inv(i) = ONE / dnA%d0(i,i)
        d0Aii_sqinv(i) = d0Aii_inv(i)*d0Aii_inv(i)
      END DO
!     write(out_unitp,*) ' det(A)=det(g) :',d0d

      DO i=dnA%nb_var_Matl+1,ncart
        mass = masses(i)
        IF (mass == ZERO) mass = ONE
        dndetA%d0 = dndetA%d0/mass
      END DO
      dndetA%d0 = dndetA%d0/(Mtot_inv**3)
!     write(out_unitp,*) ' det(A)=det(g) :',d0d
      dndetA%d0 = sqrt(dndetA%d0)
!----------------------------------------------------------

!-----------------------------------------------------------
      IF (nderiv .GE. 1) THEN
      DO ider=1,dnA%nb_var_deriv
        dndetA%d1(ider) = ZERO
        DO i=1,dnA%nb_var_Matl
          dndetA%d1(ider) = dndetA%d1(ider) + dnA%d1(i,i,ider)*d0Aii_inv(i)
        END DO
        dndetA%d1(ider) = HALF*dndetA%d1(ider)
      END DO
      END IF
!-----------------------------------------------------------


!-----------------------------------------------------------
      IF (nderiv >= 2) THEN
      DO ider=1,dnA%nb_var_deriv
      DO jder=ider,dnA%nb_var_deriv
        dndetA%d2(ider,jder) = ZERO
        DO i=1,dnA%nb_var_Matl
          dndetA%d2(ider,jder) = dndetA%d2(ider,jder) +                 &
                               dnA%d2(i,i,ider,jder)*d0Aii_inv(i) -     &
                  dnA%d1(i,i,ider)*dnA%d1(i,i,jder)*d0Aii_sqinv(i)
        END DO
        dndetA%d2(ider,jder) = HALF * dndetA%d2(ider,jder) +            &
                                         dndetA%d1(ider)*dndetA%d1(jder)
        dndetA%d2(jder,ider) = dndetA%d2(ider,jder)
      END DO
      END DO
      END IF
!-----------------------------------------------------------

!-----------------------------------------------------------
!     - restore d0A d1A d2A ----
      IF (keep) THEN
        CALL sub_dnMat1_TO_dnMat2(dnA_save,dnA,nderiv)
        CALL dealloc_dnSVM(dnA_save)
      END IF
!-----------------------------------------------------------

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'dn lnJ (from g)'
         CALL Write_dnSVM(dndetA)
         write(out_unitp,*) 'END ',name_sub
       END IF
!-----------------------------------------------------------

      end subroutine sub3_dndetA

!=====================================================================
!
! ++   calculation of  the Jacobian J=sqrt(det(g)) and its log derivatives
!      using pivot algorithmes
!      J               = sqrt(det(g))
!      (dJ/dQi) / J    = 1/2 ( ddet(g)/dQi ) / det(g)
!      (dJ/dQidQj) / J = 1/2 ( ddet(g)/dQidQj ) / det(g) - (dJ/dQi) / J * (dJ/dQj) / J
!=====================================================================
!
      SUBROUTINE sub3_dnDetGG(dndetA,dnGG,nderiv,                      &
                              masses,Mtot_inv,ncart)
      USE mod_system
      USE mod_dnSVM
      IMPLICIT NONE


      integer :: nderiv

      TYPE(Type_dnS)    :: dndetA

      TYPE(Type_dnMat)  :: dnGG
      real (kind=Rkind) :: d0Aii_inv(dnGG%nb_var_Matl)
      real (kind=Rkind) :: d0Aii_sqinv(dnGG%nb_var_Matl)

      integer :: ncart
      real (kind=Rkind) :: mass,masses(ncart),Mtot_inv,det

      TYPE(Type_dnMat)  :: dnGG_save
      integer :: i
      integer :: ider,jder

      integer :: err,memory
      logical :: keep = .TRUE.
      !logical :: new  = .TRUE.
      logical :: new  = .FALSE.

!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'sub3_dnDetGG'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nderiv',nderiv

        write(out_unitp,*) 'dnGG'
        CALL Write_dnSVM(dnGg,nderiv)
      END IF
!-----------------------------------------------------------

!-----------------------------------------------------------
      IF (keep) THEN
        CALL alloc_dnSVM(dnGG_save,dnGG%nb_var_Matl,dnGG%nb_var_Matc,   &
                         dnGG%nb_var_deriv,nderiv)
        CALL sub_dnMat1_TO_dnMat2(dnGG,dnGG_save,nderiv)
      END IF


      IF (new) THEN
        CALL Det_OF_m1(dnGG%d0,det,dnGG%nb_var_Matl)
        write(out_unitp,*) 'nGG%d0 and det,jac',det,ONE/sqrt(det)
        CALL Write_Mat(dnGG%d0,6,6)
        CALL Det_OF_dnMat_TO_dnS(dnGG,dndetA,nderiv)
        mass = product(masses,mask=(masses > ONETENTH**5))
        write(out_unitp,*) 'masses',masses

        write(out_unitp,*) 'mass',mass


        IF (nderiv >= 1) THEN ! d_i det / det
          dndetA%d1(:) = -HALF * dndetA%d1(:)/dndetA%d0
        END IF

        IF (nderiv >= 2) THEN ! d_ij det / det
          DO ider=1,dnGG%nb_var_deriv
          DO jder=ider,dnGG%nb_var_deriv
            dndetA%d2(ider,jder) = -HALF * dndetA%d2(ider,jder)/dndetA%d0 + &
                                     THREE * dndetA%d1(ider)*dndetA%d1(jder)
            dndetA%d2(jder,ider) = dndetA%d2(ider,jder)
          END DO
          END DO
        END IF

        dndetA%d0 = ONE / sqrt(dndetA%d0 * mass)
        !dndetA%d0 = sqrt(dndetA%d0 * mass)


      ELSE



!-----------------------------------------------------------

      dndetA%d0 = ONE

!-----------------------------------------------------------
!-----------------------------------------------------------
      DO i=1,dnGG%nb_var_Matl

        CALL sub3_pivot(i,dndetA%d0,dnGG,nderiv)
        CALL sub3_newline(i,dnGG,nderiv)

      END DO
!-----------------------------------------------------------
!-----------------------------------------------------------

!-----------------------------------------------------------
!     determinant and calculation of 1./d0A(i,i) 1./d0A(i,i)**2
      DO i=1,dnGG%nb_var_Matl
        mass = masses(i)
        IF (mass == ZERO) mass = ONE
!       d0d = d0d*d0A(i,i)/(mass*mass)
        dndetA%d0 = dndetA%d0 * dnGG%d0(i,i)*mass
        d0Aii_inv(i) = ONE / dnGG%d0(i,i)
        d0Aii_sqinv(i) = d0Aii_inv(i)*d0Aii_inv(i)
      END DO
      !write(out_unitp,*) ' det(GG) :',dndetA%d0

!      DO i=dnGG%nb_var_Matl+1,ncart
!        mass = masses(i)
!        IF (mass == ZERO) mass = ONE
!        dndetA%d0 = dndetA%d0 * mass
!      END DO
!      dndetA%d0 = dndetA%d0 * (Mtot_inv**3)
      !write(out_unitp,*) ' det(GG) :',dndetA%d0
      dndetA%d0 = ONE/sqrt(dndetA%d0)
!----------------------------------------------------------

!-----------------------------------------------------------
      IF (nderiv >= 1) THEN
      DO ider=1,dnGG%nb_var_deriv
        dndetA%d1(ider) = ZERO
        DO i=1,dnGG%nb_var_Matl
          dndetA%d1(ider) = dndetA%d1(ider) +                           &
                              dnGG%d1(i,i,ider)*d0Aii_inv(i)
        END DO
        dndetA%d1(ider) = -HALF*dndetA%d1(ider)
      END DO
      END IF
!-----------------------------------------------------------


!-----------------------------------------------------------
      IF (nderiv >= 2) THEN
      DO ider=1,dnGG%nb_var_deriv
      DO jder=ider,dnGG%nb_var_deriv
        dndetA%d2(ider,jder) = ZERO
        DO i=1,dnGG%nb_var_Matl
          dndetA%d2(ider,jder) = dndetA%d2(ider,jder) +                 &
                               dnGG%d2(i,i,ider,jder)*d0Aii_inv(i) -    &
                  dnGG%d1(i,i,ider)*dnGG%d1(i,i,jder)*d0Aii_sqinv(i)
        END DO
        dndetA%d2(ider,jder) = -HALF * dndetA%d2(ider,jder) +           &
                                    dndetA%d1(ider)*dndetA%d1(jder)
        dndetA%d2(jder,ider) = dndetA%d2(ider,jder)
      END DO
      END DO
      END IF
!-----------------------------------------------------------

      END IF


!-----------------------------------------------------------
!     - restore d0A d1A d2A ----
      IF (keep) THEN
        CALL sub_dnMat1_TO_dnMat2(dnGG_save,dnGG,nderiv)
        CALL dealloc_dnSVM(dnGG_save)
      END IF
!-----------------------------------------------------------

!-----------------------------------------------------------

       IF (debug) THEN
         write(out_unitp,*) 'dn lnJ (from G)'
         CALL Write_dnSVM(dndetA)
         write(out_unitp,*) 'END ',name_sub
       END IF
!-----------------------------------------------------------

      end subroutine sub3_dnDetGG
!=====================================================================
!
! ++   calculation of determinant of a matrix, Mat, LndetMat = det(Mat) and
!          its log derivatives using pivot algorithmes
!      LndetMat         = det(Mat)
!      dLndetMat/dQi    = ddet(Mat)/dQi / det(Mat)
!      dLndetMat/dQidQj = ddet(Mat)/dQidQj / det(Mat) - dLndetMat/dQi * dLndetMat/dQj
!=====================================================================
!
      SUBROUTINE sub3_dndetMat(dnLndetMat,dnMat,nderiv,                 &
                               masses,Mtot_inv,ncart)
      USE mod_system
      USE mod_dnSVM
      IMPLICIT NONE


      TYPE(Type_dnS)    :: dnLndetMat
      TYPE(Type_dnMat)  :: dnMat
      integer           :: nderiv
      integer           :: ncart
      real (kind=Rkind) :: masses(ncart),Mtot_inv

      TYPE(Type_dnMat)  :: dnMat_save
      real (kind=Rkind) :: d0Matii_inv(dnMat%nb_var_Matl)
      real (kind=Rkind) :: d0Matii_sqinv(dnMat%nb_var_Matl)
      real (kind=Rkind) :: mass

      integer           :: i
      integer           :: ider,jder

      logical           :: keep = .TRUE.
!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'sub3_dndetGG'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nderiv',nderiv

        write(out_unitp,*) 'dnMat'
        CALL Write_dnSVM(dnMat,nderiv)
      END IF
!-----------------------------------------------------------

!-----------------------------------------------------------
      IF (keep) THEN
        CALL alloc_dnSVM(dnMat_save,dnMat%nb_var_Matl,dnMat%nb_var_Matc,   &
                         dnMat%nb_var_deriv,nderiv)
        CALL sub_dnMat1_TO_dnMat2(dnMat,dnMat_save,nderiv)
      END IF

      dnLndetMat%d0 = ONE
!-----------------------------------------------------------


!-----------------------------------------------------------
!-----------------------------------------------------------
      DO i=1,dnMat%nb_var_Matl

        CALL sub3_pivot(i,dnLndetMat%d0,dnMat,nderiv)
        CALL sub3_newline(i,dnMat,nderiv)

      END DO
!-----------------------------------------------------------
!-----------------------------------------------------------

!-----------------------------------------------------------
!     determinant and calculation of 1./d0A(i,i) 1./d0A(i,i)**2
      DO i=1,dnMat%nb_var_Matl
        mass = masses(i)
        IF (mass == ZERO) mass = ONE
!       d0d = d0d*d0A(i,i)/(mass*mass)
        dnLndetMat%d0 = dnLndetMat%d0 * dnMat%d0(i,i)*mass
        d0Matii_inv(i) = ONE / dnMat%d0(i,i)
        d0Matii_sqinv(i) = d0Matii_inv(i)*d0Matii_inv(i)
      END DO
!     write(out_unitp,*) ' det(A)=det(g) :',d0d

      DO i=dnMat%nb_var_Matl+1,ncart
        mass = masses(i)
        IF (mass == ZERO) mass = ONE
        dnLndetMat%d0 = dnLndetMat%d0 * mass
      END DO
      dnLndetMat%d0 = dnLndetMat%d0 * (Mtot_inv**3)
!     write(out_unitp,*) ' det(A)=det(g) :',d0d
!----------------------------------------------------------

!-----------------------------------------------------------
      IF (nderiv >= 1) THEN
      DO ider=1,dnMat%nb_var_deriv
        dnLndetMat%d1(ider) = ZERO
        DO i=1,dnMat%nb_var_Matl
          dnLndetMat%d1(ider) = dnLndetMat%d1(ider) +                   &
                                dnMat%d1(i,i,ider)*d0Matii_inv(i)
        END DO
      END DO
      END IF
!-----------------------------------------------------------


!-----------------------------------------------------------
      IF (nderiv >= 2) THEN
      DO ider=1,dnMat%nb_var_deriv
      DO jder=ider,dnMat%nb_var_deriv
        dnLndetMat%d2(ider,jder) = ZERO
        DO i=1,dnMat%nb_var_Matl
          dnLndetMat%d2(ider,jder) = dnLndetMat%d2(ider,jder) +         &
                             dnMat%d2(i,i,ider,jder) * d0Matii_inv(i) - &
                dnMat%d1(i,i,ider)*dnMat%d1(i,i,jder) * d0Matii_sqinv(i)
        END DO
        dnLndetMat%d2(jder,ider) = dnLndetMat%d2(ider,jder)
      END DO
      END DO
      END IF
!-----------------------------------------------------------

!-----------------------------------------------------------
!     - restore d0A d1A d2A ----
      IF (keep) THEN
        CALL sub_dnMat1_TO_dnMat2(dnMat_save,dnMat,nderiv)
        CALL dealloc_dnSVM(dnMat_save)
      END IF
!-----------------------------------------------------------


!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'dn ln(det(Mat))'
         CALL Write_dnSVM(dnLndetMat)
         write(out_unitp,*) 'END ',name_sub
       END IF
!-----------------------------------------------------------

      end subroutine sub3_dndetMat
!================================================================
!       seach of the pivot
!       and permutation
!================================================================
      SUBROUTINE sub3_pivot(i,d0d,dnA,nderiv)
      USE mod_system
      USE mod_dnSVM
      IMPLICIT NONE

      integer nderiv
      real (kind=Rkind) ::  d0d

      TYPE(Type_dnMat)  :: dnA
      real (kind=Rkind) ::  piv_i,piv_j
      integer i,j,line
      integer ider,jder

!----------------------------------------------------------------
!     seach for the larger d0A(j,i)
!----------------------------------------------------------------
      piv_i = abs(dnA%d0(i,i))
      line = i
      DO j=i+1,dnA%nb_var_Matl
        piv_j = abs(dnA%d0(j,i))
        IF (piv_j > piv_i ) THEN
          piv_i = piv_j
          line = j
        END IF
      END DO
!----------------------------------------------------------------

!----------------------------------------------------------------
!     permutation of the line : i and line
!----------------------------------------------------------------
      CALL sub3_permutation(i,line,dnA%d0,dnA%nb_var_Matl)
      IF (nderiv .GE. 1) THEN
        DO ider=1,dnA%nb_var_deriv
          CALL sub3_permutation(i,line,dnA%d1(:,:,ider),dnA%nb_var_Matl)
        END DO
      END IF
      IF (nderiv .GE. 2) THEN
        DO ider=1,dnA%nb_var_deriv
        DO jder=ider,dnA%nb_var_deriv
         CALL sub3_permutation(i,line,                                  &
                 dnA%d2(:,:,ider,jder),dnA%nb_var_Matl)
        END DO
        END DO
      END IF
!----------------------------------------------------------------

!----------------------------------------------------------------
!     sign of the permutation
!----------------------------------------------------------------
      IF ( line .NE. i) THEN
        d0d = -d0d
      END IF
!----------------------------------------------------------------

!     write(out_unitp,*) 'perm',i,line,d0d,d0A(i,i)

      end subroutine sub3_pivot
!================================================================
!     determine the new line with the pivot
!================================================================
      SUBROUTINE sub3_newline(i,dnA,nderiv)
      USE mod_system
      USE mod_dnSVM
      IMPLICIT NONE


      integer           :: nderiv
      TYPE(Type_dnMat)  :: dnA
      TYPE(Type_dnS)    :: dnpj,dnai
      real (kind=Rkind) :: d0ai2,d0ai3
      integer :: i,j,k
      integer :: ider,jder


      CALL alloc_dnSVM(dnai,dnA%nb_var_deriv,nderiv)
      CALL alloc_dnSVM(dnpj,dnA%nb_var_deriv,nderiv)


      dnai%d0 = ONE / dnA%d0(i,i)
      d0ai2   = dnai%d0 * dnai%d0
      d0ai3   = d0ai2   * dnai%d0

      DO j=i+1,dnA%nb_var_Matl
!       ---------------------------------------------------------
        dnpj%d0 = dnA%d0(j,i)*dnai%d0
!       ---------------------------------------------------------

!       ---------------------------------------------------------
        IF (nderiv .GE. 1) THEN
          DO ider=1,dnA%nb_var_deriv
            dnai%d1(ider) =-dnA%d1(i,i,ider) * d0ai2
            dnpj%d1(ider) = dnA%d1(j,i,ider) * dnai%d0      +           &
                            dnA%d0(j,i)      * dnai%d1(ider)
          END DO
        END IF
!       ---------------------------------------------------------

!       ---------------------------------------------------------
        IF (nderiv .GE. 2) THEN
          DO ider=1,dnA%nb_var_deriv
          DO jder=ider,dnA%nb_var_deriv
            dnai%d2(ider,jder) = -dnA%d2(i,i,ider,jder)*d0ai2 +         &
                           TWO*dnA%d1(i,i,ider)*dnA%d1(i,i,jder)*d0ai3
            dnpj%d2(ider,jder) = dnA%d2(j,i,ider,jder) * dnai%d0       +&
                                dnA%d1(j,i,ider)       * dnai%d1(jder) +&
                                dnA%d1(j,i,jder)       * dnai%d1(ider) +&
                                dnA%d0(j,i)         * dnai%d2(ider,jder)
          END DO
          END DO
        END IF
!       ---------------------------------------------------------

        IF (nderiv .GE. 2) THEN
          DO ider=1,dnA%nb_var_deriv
          DO jder=ider,dnA%nb_var_deriv
            DO k=i+1,dnA%nb_var_Matl
             dnA%d2(j,k,ider,jder) = dnA%d2(j,k,ider,jder)           -  &
                                dnpj%d2(ider,jder) * dnA%d0(i,k)      - &
                                dnpj%d1(ider)      * dnA%d1(i,k,jder) - &
                                dnpj%d1(jder)      * dnA%d1(i,k,ider) - &
                                dnpj%d0          * dnA%d2(i,k,ider,jder)
           END DO
          END DO
          END DO
        END IF

        IF (nderiv .GE. 1) THEN
          DO ider=1,dnA%nb_var_deriv
            DO k=i+1,dnA%nb_var_Matl
             dnA%d1(j,k,ider) = dnA%d1(j,k,ider) -                      &
                                dnpj%d1(ider) * dnA%d0(i,k) -           &
                                dnpj%d0       * dnA%d1(i,k,ider)
           END DO
          END DO
        END IF

        DO k=i+1,dnA%nb_var_Matl
          dnA%d0(j,k) = dnA%d0(j,k) - dnpj%d0*dnA%d0(i,k)
        END DO



      END DO

      CALL dealloc_dnSVM(dnai)
      CALL dealloc_dnSVM(dnpj)

      end subroutine sub3_newline
!================================================================
!       permutation of lines i and j
!================================================================
      SUBROUTINE sub3_permutation(i,j,A,ndimA)
      USE mod_system
      IMPLICIT NONE

      integer ndimA
      real (kind=Rkind) ::  A(ndimA,ndimA)
      real (kind=Rkind) ::  piv
      integer i,j,k

!----------------------------------------------------------------
!     permutation of the line i and j
!----------------------------------------------------------------
      IF (j .NE. i) THEN
        DO k=i,ndimA
          piv = A(i,k)
          A(i,k) = A(j,k)
          A(j,k) = piv
        END DO
      END IF
!     write(out_unitp,*) 'perm',i,j,d,A(i,i)
!----------------------------------------------------------------

      end subroutine sub3_permutation
!=====================================================================
!
! ++   calcul du determinant de A
!
!=====================================================================
!
      SUBROUTINE sub_detA(d,A,trav,index,ndimA,masses,Mtot_inv,ncart)
      USE mod_system
      IMPLICIT NONE

       integer ndimA
       real (kind=Rkind) ::  A(ndimA,ndimA)
       real (kind=Rkind) :: trav(ndimA)
       integer index(ndimA)
       real (kind=Rkind) ::  d
       integer ncart
       real (kind=Rkind) :: mass,masses(ncart),Mtot_inv

       integer i


!      determinant de A

       CALL ludcmp(A,ndimA,trav,index,d)

!      determinant
       DO i=1,ndimA
        mass = masses(i)
        IF (mass == ZERO) mass = ONE
        d = d*A(i,i)/(mass*mass)
       END DO

       DO i=ndimA+1,ncart
        mass = masses(i)
        IF (mass == ZERO) mass = ONE
        d = d/(mass*mass)
       END DO

       d = d/(Mtot_inv**3)
!      write(out_unitp,*) d

       end subroutine sub_detA
END MODULE mod_dnDetGG_dnDetg
