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
MODULE mod_freq
      use mod_system
      USE mod_Constant, ONLY: get_Conv_au_TO_unit
      IMPLICIT NONE

      PRIVATE

      PUBLIC :: calc_freq, calc_freq_new, calc_freq_block, calc_freq_WITH_d0c,  &
                calc_freqNM, calc_freq_width
      PUBLIC :: H0_symmetrization, sort_with_Tab, gaussian_width
      PUBLIC :: degenerate_freq_t,Write_degenerate_freq,Read_degenerate_freq,   &
                Init_degenerate_freq,dealloc_degenerate_freq

      TYPE degenerate_freq_t
        integer                      :: option = 0 ! 0 old one, 1 new one
        logical, allocatable         :: list_k_check(:)
        real (kind=Rkind)            :: epsi = ONETENTH**8
      END TYPE degenerate_freq_t


      CONTAINS

      SUBROUTINE Write_degenerate_freq(degenerate_freq)
      IMPLICIT NONE
      TYPE (degenerate_freq_t), INTENT(IN) :: degenerate_freq

        write(out_unitp,*) 'BEGINNING Write_degenerate_freq'

        write(out_unitp,*) 'option',degenerate_freq%option
        write(out_unitp,*) 'epsi',degenerate_freq%epsi
        IF (allocated(degenerate_freq%list_k_check)) THEN
          write(out_unitp,*) 'list_k_check(:)',degenerate_freq%list_k_check
        END IF

        write(out_unitp,*) 'END Write_degenerate_freq'
        flush(out_unitp)

      END SUBROUTINE Write_degenerate_freq
      SUBROUTINE Read_degenerate_freq(degene_freq,n)
      IMPLICIT NONE
      TYPE (degenerate_freq_t), INTENT(inout) :: degene_freq
      integer,                  INTENT(in)    :: n

      integer                      :: option,err_read
      real (kind=Rkind)            :: epsi

      NAMELIST /degenerate_freq/ option,epsi
      character(len=*), PARAMETER :: name_sub = 'Read_degenerate_freq'

        write(out_unitp,*) 'BEGINNING ',name_sub

        CALL alloc_NParray(degene_freq%list_k_check,[n],'list_k_check',name_sub)

        epsi   = ONETENTH**8
        option = 0
       read(in_unitp,degenerate_freq,IOSTAT=err_read)
       IF (err_read /= 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  while reading the "degenerate_freq" namelist'
          write(out_unitp,*) ' end of file or end of record'
          write(out_unitp,*) ' Check your data !!'
          STOP
       END IF
       write(out_unitp,degenerate_freq)

       degene_freq%option = option
       degene_freq%epsi   = epsi
       IF (option == 1) THEN
         read(in_unitp,*,IOSTAT=err_read) degene_freq%list_k_check(:)
         !write(out_unitp,*) 'list_k_check',degene_freq%list_k_check
         IF (err_read /= 0) THEN
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) '  while reading "list_k_check"'
           write(out_unitp,*) ' end of file or end of record'
           write(out_unitp,*) ' Check your data !!'
           STOP
         END IF
       END IF
       flush(out_unitp)

       CALL Write_degenerate_freq(degene_freq)

        write(out_unitp,*) 'END ',name_sub
        flush(out_unitp)

      END SUBROUTINE Read_degenerate_freq
      SUBROUTINE Init_degenerate_freq(degene_freq,n)
      IMPLICIT NONE
      TYPE (degenerate_freq_t), INTENT(inout) :: degene_freq
      integer,                  INTENT(in)    :: n


        write(out_unitp,*) 'BEGINNING Init_degenerate_freq'

        CALL alloc_NParray(degene_freq%list_k_check,[n],'list_k_check','Init_degenerate_freq')
        degene_freq%list_k_check(:) = .TRUE.
        degene_freq%epsi   = ONETENTH**8
        degene_freq%option = 0

        write(out_unitp,*) 'END Init_degenerate_freq'
        flush(out_unitp)

      END SUBROUTINE Init_degenerate_freq
      SUBROUTINE dealloc_degenerate_freq(degene_freq)
      IMPLICIT NONE
      TYPE (degenerate_freq_t), INTENT(inout) :: degene_freq

        IF (allocated(degene_freq%list_k_check)) THEN
          CALL dealloc_NParray(degene_freq%list_k_check,'list_k_check','dealloc_degenerate_freq')
        END IF
        degene_freq%option = 0
        degene_freq%epsi   = ONETENTH**8

      END SUBROUTINE dealloc_degenerate_freq
!=====================================================================
!
! ++   calcule les frequences et les modes normaux
!
!=====================================================================
      SUBROUTINE calc_freq_new(nb_var,d0h,d0k,d0eh,                         &
                           d0c,d0c_inv,norme,d0c_ini,diab_freq,degenerate_freq)
      IMPLICIT NONE

      integer :: nb_var
      logical :: diab_freq

      real (kind=Rkind) :: d0h(nb_var,nb_var)
      real (kind=Rkind) :: d0k(nb_var,nb_var)

      real (kind=Rkind) :: d0eh(nb_var)
      real (kind=Rkind) :: d0ek(nb_var)

      real (kind=Rkind) :: d0ch(nb_var,nb_var)

      real (kind=Rkind) :: d0ck(nb_var,nb_var)

      real (kind=Rkind) :: d0c(nb_var,nb_var)
      real (kind=Rkind) :: d0c_inv(nb_var,nb_var)
      real (kind=Rkind) :: d0c_ini(nb_var,nb_var)
      TYPE (degenerate_freq_t), INTENT(IN) :: degenerate_freq

      real (kind=Rkind) :: mat1(nb_var,nb_var)
      real (kind=Rkind) :: mat2(nb_var,nb_var)

      real (kind=Rkind) :: val,val1
      !real (kind=Rkind),parameter :: epsi_freq = ONETENTH**7
      real (kind=Rkind),parameter :: epsi_freq = ONETENTH**10

      !----- pour le determinant de d0c
      real (kind=Rkind) ::    d,norme,max_err
      real (kind=Rkind) ::    dh,dk


      integer :: i,j,k,l
      integer :: ierr

      !-----------------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      !-----------------------------------------------------------
      IF (debug) THEN
         write(out_unitp,*) 'BEGINNING calc_freq_new'
         write(out_unitp,*) 'd0h',nb_var
         CALL Write_Mat(d0h,out_unitp,5)
         write(out_unitp,*) 'd0k',nb_var
         CALL Write_Mat(d0k,out_unitp,5)
         write(out_unitp,*) 'd0c_ini',nb_var
         CALL Write_Mat(d0c_ini,out_unitp,5)
         flush(out_unitp)
      END IF

      !-----------------------------------------------------------

!      pour le cas 3 du calcul de det(d0c)
!      - cas 3 --------------------------------
!      trav2 est utilise pour index ...
!      CALL copy_mat(d0ck,d0k,nb_var,nb_var)
!      CALL ludcmp(d0ck,nb_var,trav1,trav2,dk)
!      CALL copy_mat(d0ch,d0h,nb_var,nb_var)
!      CALL ludcmp(d0ch,nb_var,trav1,trav2,dh)
!      d = dh/dk
!      DO i=1,nb_var
!       d = d * d0ch(i,i)/d0ck(i,i)
!      END DO
!      write(out_unitp,*) 'det d0c',sqrt(sqrt(d))
!      ----------------------------------------


!-----------------------------------------------------------
!----- orthonormalisation de Lowdin ------------------------
!-----------------------------------------------------------
!      de la partie cinetique d0k
!      passage de la base b0 (coordonnees initiales)
!      a b1 (coordonnees tq d0k soit diagonale)
       CALL diagonalization(d0k,d0ek,d0ck,nb_var,2,1,.TRUE.)
       CALL rota_denerated(d0ek,d0ck,nb_var)
       CALL mat_epsiTOzero(d0ck,nb_var,epsi_freq,nb_var)

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'vp de d0k'
         write(out_unitp,*) (d0ek(i),i=1,nb_var)
         flush(out_unitp)
       END IF
!-----------------------------------------------------------
       DO i=1,nb_var
         IF (d0ek(i) >= ZERO) THEN
           d0ek(i) =  sqrt(d0ek(i))
         ELSE
           d0ek(i) = -sqrt(-d0ek(i))
           write(out_unitp,*) 'ERROR: d0k has one negative eigenvalue !'
           !STOP
         END IF
       END DO

!----- vecteurs normalises par sqrt(d0ek(i)) ---------------
       DO i=1,nb_var
         d0ck(:,i) = d0ck(:,i)*d0ek(i)
       END DO

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*)
         write(out_unitp,*) 'd0ck(:,i)*sqrt(d0ek(i))'
         CALL Write_Mat(d0ck,out_unitp,5)
         flush(out_unitp)
       END IF
!-----------------------------------------------------------



!-----------------------------------------------------------
!----- calcul du hessien dans la nouvelle base -------------
!      et diagonalisation obtention des modes normaux (base b2)
!      exprimes dans la base b1
!      -- problem with gfortran ---
!      d0k = matmul( transpose(d0ck) , matmul(d0h,d0ck) )
!      the line is split ---
       mat1 = matmul(d0h,d0ck)
       mat2 = transpose(d0ck)
       d0k = matmul(mat2,mat1)
!      -- problem with gfortran ---
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*)
         write(out_unitp,*) 'd0h dans la nouvelle base'
         CALL Write_Mat(d0k,out_unitp,5)
         flush(out_unitp)
       END IF
!-----------------------------------------------------------

       CALL diagonalization(d0k,d0eh,d0ch,nb_var,2,1,.TRUE.)
       CALL rota_denerated(d0eh,d0ch,nb_var)
       CALL mat_epsiTOzero(d0ch,nb_var,epsi_freq,nb_var)


       DO i=1,nb_var
         IF (d0eh(i) >= ZERO) THEN
           d0eh(i) =  sqrt(d0eh(i))
         ELSE
           d0eh(i) =  sqrt(-d0eh(i))
           write(out_unitp,*) ' ERROR : one imaginary frequency',               &
                                    d0eh(i)*get_Conv_au_TO_unit('E','cm-1')
!          STOP
         END IF

       END DO

      !CALL order_ini4(d0ch,d0c_inv,d0eh,d0c_ini,nb_var,diab_freq)


!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'ZPE (cm-1): ',HALF*sum(d0eh(:))*get_Conv_au_TO_unit('E','cm-1')
         !write(out_unitp,*) 'ZPE   (eV): ',HALF*sum(d0eh(:))*get_Conv_au_TO_unit('E','eV')
         write(out_unitp,*) 'ZPE   (au): ',HALF*sum(d0eh(:))

         write(out_unitp,*) 'frequencies (cm-1): ',d0eh(:)*get_Conv_au_TO_unit('E','cm-1')
         write(out_unitp,*)
         write(out_unitp,*) 'modes normaux: d0ch'
         CALL Write_Mat(d0ch,out_unitp,5)
         flush(out_unitp)
       END IF
!-----------------------------------------------------------


!-----------------------------------------------------------
!----- calcul de la matrice qui passe des modes normaux ----
!      exprimes en fonction des coordonnees
!      passage de la base b2 a b0
!-----------------------------------------------------------



! ca marche si d0k n est pas diagonale
       mat1 = matmul(d0ck,d0ch)
       d0c_inv = transpose( mat1 )


       DO k=1,nb_var
         val = ONE/(d0ek(k)*d0ek(k))
         d0ck(:,k) = d0ck(:,k) * val
       END DO
       d0c = matmul(d0ck,d0ch)

       DO i=1,nb_var
         d0c(:,i)     = d0c(:,i)     * sqrt(d0eh(i))
         d0c_inv(i,:) = d0c_inv(i,:) / sqrt(d0eh(i))
       END DO

       !write(out_unitp,*) 'rota_denerated on d0c and d0c_inv'

       CALL rota_degenerate_opt1(d0eh,d0c,nb_var,degenerate_freq)
       !write(out_unitp,*) 'd0c'
       !CALL Write_Mat(d0c,out_unitp,5)

       CALL inv_m1_TO_m2(d0c,d0c_inv,nb_var,0,ZERO)
       !write(out_unitp,*) 'transpose(d0c_inv)'
       !CALL Write_Mat(transpose(d0c_inv),out_unitp,5)

!     on reordonne d0c, d0c_inv et d0eh
      !CALL order_ini4(d0c,d0c_inv,d0eh,d0c_ini,nb_var,diab_freq)
      CALL order_ini5(d0c,d0c_inv,d0eh,d0c_ini,nb_var,diab_freq)


! ca marche si d0k n est pas diagonale
!      DO i=1,nb_var
!        val = sqrt(d0eh(i))
!        DO l=1,nb_var
!          d0c(l,i) = ONE
!          d0c_inv(i,l) = ZERO
!          DO k=1,nb_var
!            d0c(l,i) = d0c(l,i) +
!    *          d0ch(k,i)*d0ck(l,k)*val/(d0ek(k)*d0ek(k))
!            d0c_inv(i,l) = d0c_inv(i,l) + d0ch(k,i)*d0ck(l,k)/val
!          END DO
!        END DO
!      END DO

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'matrice de passage de b2 a b0'
         CALL Write_Mat(d0c,out_unitp,5)
         write(out_unitp,*) 'matrice de passage de b0 a b2'
         CALL Write_Mat(d0c_inv,out_unitp,5)

!        test inversion ------------------------------------
         mat1 = matmul(d0c,d0c_inv)
         write(out_unitp,*) 'test inversion de d0c'
         CALL Write_Mat(mat1,out_unitp,5)
         max_err = ZERO
         DO i=1,nb_var
         DO j=1,nb_var
            IF (i /= j) max_err = max(max_err,abs(mat1(i,j)))
            IF (i == j) max_err = max(max_err,abs(mat1(i,j)-ONE))
         END DO
         END DO
         write(out_unitp,*) ' max_err: ',max_err
         flush(out_unitp)
       END IF
!-----------------------------------------------------------

!-----------------------------------------------------------
!-----------------------------------------------------------
!      calcul le determinant de d0c :
!
!      3 facons :
!
!      1)  on calcule directement d=det(d0c)
!      2)  on utilise les valeurs propres d0eh et d0ek
!      3)  d = [det(d0h)/det(d0k)]^1/4
!-----------------------------------------------------------
!-----------------------------------------------------------

!
!      - cas 1 --------------------------------
!      CALL copy_mat(d0ch,d0c,nb_var,nb_var)
!      trav2 est utilise pour index
!      CALL ludcmp(d0ch,nb_var,trav1,trav2,d)
!      DO i=1,nb_var
!       d = d * d0ch(i,i)
!      END DO
!      write(out_unitp,*) 'det d0c',d
!      ----------------------------------------


!      - cas 2 --------------------------------
!      utilise le fait que det(d0ch) et det(d0k) =1
       d = ONE
       DO i=1,nb_var
        d = d * sqrt(d0eh(i)) / d0ek(i)
       END DO
!      write(out_unitp,*) 'det d0c',d
!      ----------------------------------------
!
!      - cas 3 --------------------------------
!      il est fait au debut
!      ----------------------------------------


!-----------------------------------------------------------
!-----------------------------------------------------------

!      norme = ONE/sqrt(d)
       norme = sqrt(d)

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'determinants :',norme
         write(out_unitp,*) 'END calc_freq_new'
         flush(out_unitp)
       END IF
!-----------------------------------------------------------
!stop

end subroutine calc_freq_new

      SUBROUTINE calc_freq(nb_var,d0h,d0k,d0eh,                         &
                           d0c,d0c_inv,norme,d0c_ini,diab_freq)
      IMPLICIT NONE

      integer :: nb_var
      logical :: diab_freq

      real (kind=Rkind) :: d0h(nb_var,nb_var)
      real (kind=Rkind) :: d0k(nb_var,nb_var)

      real (kind=Rkind) :: d0eh(nb_var)
      real (kind=Rkind) :: d0ek(nb_var)

      real (kind=Rkind) :: d0ch(nb_var,nb_var)

      real (kind=Rkind) :: d0ck(nb_var,nb_var)

      real (kind=Rkind) :: d0c(nb_var,nb_var)
      real (kind=Rkind) :: d0c_inv(nb_var,nb_var)
      real (kind=Rkind) :: d0c_ini(nb_var,nb_var)
      logical           :: list_k_check(nb_var)

      real (kind=Rkind) :: mat1(nb_var,nb_var)
      real (kind=Rkind) :: mat2(nb_var,nb_var)

      real (kind=Rkind) :: val,val1
      !real (kind=Rkind),parameter :: epsi_freq = ONETENTH**7
      real (kind=Rkind),parameter :: epsi_freq = ONETENTH**10

      !----- pour le determinant de d0c
      real (kind=Rkind) ::    d,norme,max_err
      real (kind=Rkind) ::    dh,dk


      integer :: i,j,k,l
      integer :: ierr

      !-----------------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      !-----------------------------------------------------------
      IF (debug) THEN
         write(out_unitp,*) 'BEGINNING calc_freq'
         write(out_unitp,*) 'd0h',nb_var
         CALL Write_Mat(d0h,out_unitp,5)
         write(out_unitp,*) 'd0k',nb_var
         CALL Write_Mat(d0k,out_unitp,5)
         write(out_unitp,*) 'd0c_ini',nb_var
         CALL Write_Mat(d0c_ini,out_unitp,5)
         flush(out_unitp)
      END IF
      !-----------------------------------------------------------

!      pour le cas 3 du calcul de det(d0c)
!      - cas 3 --------------------------------
!      trav2 est utilise pour index ...
!      CALL copy_mat(d0ck,d0k,nb_var,nb_var)
!      CALL ludcmp(d0ck,nb_var,trav1,trav2,dk)
!      CALL copy_mat(d0ch,d0h,nb_var,nb_var)
!      CALL ludcmp(d0ch,nb_var,trav1,trav2,dh)
!      d = dh/dk
!      DO i=1,nb_var
!       d = d * d0ch(i,i)/d0ck(i,i)
!      END DO
!      write(out_unitp,*) 'det d0c',sqrt(sqrt(d))
!      ----------------------------------------


!-----------------------------------------------------------
!----- orthonormalisation de Lowdin ------------------------
!-----------------------------------------------------------
!      de la partie cinetique d0k
!      passage de la base b0 (coordonnees initiales)
!      a b1 (coordonnees tq d0k soit diagonale)
       CALL diagonalization(d0k,d0ek,d0ck,nb_var,2,1,.TRUE.)
       CALL rota_denerated(d0ek,d0ck,nb_var)
       CALL mat_epsiTOzero(d0ck,nb_var,epsi_freq,nb_var)

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'vp de d0k'
         write(out_unitp,*) (d0ek(i),i=1,nb_var)
         flush(out_unitp)
       END IF
!-----------------------------------------------------------
       DO i=1,nb_var
         IF (d0ek(i) >= ZERO) THEN
           d0ek(i) =  sqrt(d0ek(i))
         ELSE
           d0ek(i) = -sqrt(-d0ek(i))
           write(out_unitp,*) 'ERROR: d0k has one negative eigenvalue !'
           !STOP
         END IF
       END DO

!----- vecteurs normalises par sqrt(d0ek(i)) ---------------
       DO i=1,nb_var
         d0ck(:,i) = d0ck(:,i)*d0ek(i)
       END DO

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*)
         write(out_unitp,*) 'd0ck(:,i)*sqrt(d0ek(i))'
         CALL Write_Mat(d0ck,out_unitp,5)
         flush(out_unitp)
       END IF
!-----------------------------------------------------------



!-----------------------------------------------------------
!----- calcul du hessien dans la nouvelle base -------------
!      et diagonalisation obtention des modes normaux (base b2)
!      exprimes dans la base b1
!      -- problem with gfortran ---
!      d0k = matmul( transpose(d0ck) , matmul(d0h,d0ck) )
!      the line is split ---
       mat1 = matmul(d0h,d0ck)
       mat2 = transpose(d0ck)
       d0k = matmul(mat2,mat1)
!      -- problem with gfortran ---
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*)
         write(out_unitp,*) 'd0h dans la nouvelle base'
         CALL Write_Mat(d0k,out_unitp,5)
         flush(out_unitp)
       END IF
!-----------------------------------------------------------

       CALL diagonalization(d0k,d0eh,d0ch,nb_var,2,1,.TRUE.)
       CALL rota_denerated(d0eh,d0ch,nb_var)
       CALL mat_epsiTOzero(d0ch,nb_var,epsi_freq,nb_var)


       DO i=1,nb_var
         IF (d0eh(i) >= ZERO) THEN
           d0eh(i) =  sqrt(d0eh(i))
         ELSE
           d0eh(i) =  sqrt(-d0eh(i))
           write(out_unitp,*) ' ERROR : one imaginary frequency',               &
                                    d0eh(i)*get_Conv_au_TO_unit('E','cm-1')
!          STOP
         END IF

       END DO

      !CALL order_ini4(d0ch,d0c_inv,d0eh,d0c_ini,nb_var,diab_freq)


!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'ZPE (cm-1): ',HALF*sum(d0eh(:))*get_Conv_au_TO_unit('E','cm-1')
         !write(out_unitp,*) 'ZPE   (eV): ',HALF*sum(d0eh(:))*get_Conv_au_TO_unit('E','eV')
         write(out_unitp,*) 'ZPE   (au): ',HALF*sum(d0eh(:))

         write(out_unitp,*) 'frequencies (cm-1): ',d0eh(:)*get_Conv_au_TO_unit('E','cm-1')
         write(out_unitp,*)
         write(out_unitp,*) 'modes normaux: d0ch'
         CALL Write_Mat(d0ch,out_unitp,5)
         flush(out_unitp)
       END IF
!-----------------------------------------------------------


!-----------------------------------------------------------
!----- calcul de la matrice qui passe des modes normaux ----
!      exprimes en fonction des coordonnees
!      passage de la base b2 a b0
!-----------------------------------------------------------



! ca marche si d0k n est pas diagonale
       mat1 = matmul(d0ck,d0ch)
       d0c_inv = transpose( mat1 )


       DO k=1,nb_var
         val = ONE/(d0ek(k)*d0ek(k))
         d0ck(:,k) = d0ck(:,k) * val
       END DO
       d0c = matmul(d0ck,d0ch)

       DO i=1,nb_var
         d0c(:,i)     = d0c(:,i)     * sqrt(d0eh(i))
         d0c_inv(i,:) = d0c_inv(i,:) / sqrt(d0eh(i))
       END DO

!     on reordonne d0c, d0c_inv et d0eh
      !CALL order_ini4(d0c,d0c_inv,d0eh,d0c_ini,nb_var,diab_freq)
      CALL order_ini5(d0c,d0c_inv,d0eh,d0c_ini,nb_var,diab_freq)


! ca marche si d0k n est pas diagonale
!      DO i=1,nb_var
!        val = sqrt(d0eh(i))
!        DO l=1,nb_var
!          d0c(l,i) = ONE
!          d0c_inv(i,l) = ZERO
!          DO k=1,nb_var
!            d0c(l,i) = d0c(l,i) +
!    *          d0ch(k,i)*d0ck(l,k)*val/(d0ek(k)*d0ek(k))
!            d0c_inv(i,l) = d0c_inv(i,l) + d0ch(k,i)*d0ck(l,k)/val
!          END DO
!        END DO
!      END DO

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'matrice de passage de b2 a b0'
         CALL Write_Mat(d0c,out_unitp,5)
         write(out_unitp,*) 'matrice de passage de b0 a b2'
         CALL Write_Mat(d0c_inv,out_unitp,5)

!        test inversion ------------------------------------
         mat1 = matmul(d0c,d0c_inv)
         write(out_unitp,*) 'test inversion de d0c'
         CALL Write_Mat(mat1,out_unitp,5)
         max_err = ZERO
         DO i=1,nb_var
         DO j=1,nb_var
            IF (i /= j) max_err = max(max_err,abs(mat1(i,j)))
            IF (i == j) max_err = max(max_err,abs(mat1(i,j)-ONE))
         END DO
         END DO
         write(out_unitp,*) ' max_err: ',max_err
         flush(out_unitp)
       END IF
!-----------------------------------------------------------

!-----------------------------------------------------------
!-----------------------------------------------------------
!      calcul le determinant de d0c :
!
!      3 facons :
!
!      1)  on calcule directement d=det(d0c)
!      2)  on utilise les valeurs propres d0eh et d0ek
!      3)  d = [det(d0h)/det(d0k)]^1/4
!-----------------------------------------------------------
!-----------------------------------------------------------

!
!      - cas 1 --------------------------------
!      CALL copy_mat(d0ch,d0c,nb_var,nb_var)
!      trav2 est utilise pour index
!      CALL ludcmp(d0ch,nb_var,trav1,trav2,d)
!      DO i=1,nb_var
!       d = d * d0ch(i,i)
!      END DO
!      write(out_unitp,*) 'det d0c',d
!      ----------------------------------------


!      - cas 2 --------------------------------
!      utilise le fait que det(d0ch) et det(d0k) =1
       d = ONE
       DO i=1,nb_var
        d = d * sqrt(d0eh(i)) / d0ek(i)
       END DO
!      write(out_unitp,*) 'det d0c',d
!      ----------------------------------------
!
!      - cas 3 --------------------------------
!      il est fait au debut
!      ----------------------------------------


!-----------------------------------------------------------
!-----------------------------------------------------------

!      norme = ONE/sqrt(d)
       norme = sqrt(d)

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'determinants :',norme
         write(out_unitp,*) 'END calc_freq'
         flush(out_unitp)
       END IF
!-----------------------------------------------------------
!stop

      end subroutine calc_freq
      SUBROUTINE calc_freq_block(nb_var,d0h,d0k,d0eh,                   &
                                 d0c,d0c_inv,norme,d0c_ini,diab_freq,sym)
      IMPLICIT NONE


      integer :: nb_var
      real (kind=Rkind) :: norme
      logical :: diab_freq
      integer :: sym(nb_var)
      real (kind=Rkind) :: d0k(nb_var,nb_var),d0h(nb_var,nb_var),d0eh(nb_var)
      real (kind=Rkind) :: d0c(nb_var,nb_var),d0c_inv(nb_var,nb_var),d0c_ini(nb_var,nb_var)



      real (kind=Rkind) :: norm_PerBlock
      real (kind=Rkind), allocatable :: d0k_PerBlock(:,:),d0h_PerBlock(:,:),d0eh_PerBlock(:)
      real (kind=Rkind), allocatable :: d0c_PerBlock(:,:),d0c_inv_PerBlock(:,:),d0c_ini_PerBlock(:,:)


      integer :: ierr

      integer :: nb_NM
      integer, allocatable :: Ind_Coord_PerBlock(:)
      integer, allocatable :: nb_PerBlock(:),Ind_Coord_AtBlock(:)
      integer :: i_Block,nb_Block
      integer ::i,j,i2,ib,jb,iNM

!      -----------------------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'calc_freq_block'
!      -----------------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*)
         write(out_unitp,*) '========================================='
         write(out_unitp,*) '========== calc_freq_block =============='
         write(out_unitp,*) '========================================='
         write(out_unitp,*)
         write(out_unitp,*) '=========, sym: ',sym(:)
       END IF

      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      ! analysis the number of block ... => nb_NM
      CALL alloc_NParray(Ind_Coord_PerBlock,[nb_var],          &
                        "Ind_Coord_PerBlock",name_sub)

      Ind_Coord_PerBlock(:) = sym(:)

      ! first count the blocks
      nb_Block = 0
      DO
        i = maxval(Ind_Coord_PerBlock)
        IF ( i /= -Huge(1)) THEN
           nb_Block = nb_Block + 1
           WHERE (Ind_Coord_PerBlock == i) Ind_Coord_PerBlock = -Huge(1)
        ELSE
           EXIT
        END IF
      END DO
      Ind_Coord_PerBlock(:) = sym(:)

      ! then, count the number of coordinates per block
      CALL alloc_NParray(nb_PerBlock,      [nb_Block],"nb_PerBlock",      name_sub)
      CALL alloc_NParray(Ind_Coord_AtBlock,[nb_Block],"Ind_Coord_AtBlock",name_sub)
      Ind_Coord_AtBlock(:) = 0
      nb_PerBlock(:)       = 0

      DO i_Block=1,nb_Block

        i = maxval(Ind_Coord_PerBlock)
        Ind_Coord_AtBlock(i_Block) = i

        IF (i /= 0) nb_PerBlock(i_Block) = count(Ind_Coord_PerBlock == i)

        WHERE (Ind_Coord_PerBlock == i) Ind_Coord_PerBlock = -Huge(1)

      END DO
      Ind_Coord_PerBlock(:) = sym(:)

      nb_NM = sum(nb_PerBlock)

      IF (debug) THEN
        write(out_unitp,*) '  nb_Block            ',nb_Block
        write(out_unitp,*) '  nb_PerBlock(:)      ',nb_PerBlock(:)
        write(out_unitp,*) '  Ind_Coord_AtBlock(:)',Ind_Coord_AtBlock(:)
        write(out_unitp,*) '  nb_NM tot',nb_NM
      END IF

      !-----------------------------------------------------------------
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      iNM          = 0
      norme        = ONE
      d0c(:,:)     = ZERO
      d0C_inv(:,:) = ZERO
      DO i_Block=1,nb_Block
        IF (Ind_Coord_AtBlock(i_Block) == 0) CYCLE

        nb_NM = nb_PerBlock(i_Block)
        IF (debug) THEN
          write(out_unitp,*) '========================================='
          write(out_unitp,*) '=========             Block: ',i_Block
          write(out_unitp,*) '=========       nb_PerBlock: ',nb_PerBlock(i_Block)
          write(out_unitp,*) '========= Ind_Coord_AtBlock: ',Ind_Coord_AtBlock(i_Block)
          write(out_unitp,*) '========================================='
          write(out_unitp,*)
        END IF

        !d0k_PerBlock, d0h_PerBlock
        CALL alloc_NParray(d0k_PerBlock,      [nb_NM,nb_NM],"d0k_PerBlock",      name_sub)
        CALL alloc_NParray(d0h_PerBlock,      [nb_NM,nb_NM],"d0h_PerBlock",      name_sub)
        CALL alloc_NParray(d0c_PerBlock,      [nb_NM,nb_NM],"d0c_PerBlock",      name_sub)
        CALL alloc_NParray(d0c_inv_PerBlock,  [nb_NM,nb_NM],"d0c_inv_PerBlock",  name_sub)
        CALL alloc_NParray(d0c_ini_PerBlock,  [nb_NM,nb_NM],"d0c_ini_PerBlock",  name_sub)

        CALL alloc_NParray(d0eh_PerBlock,     [nb_NM],      "d0eh_PerBlock",     name_sub)

        !d0k => d0k_PerBlock, d0h => d0h_PerBlock
        ib = 0
        DO i=1,nb_var
          IF (sym(i) /= Ind_Coord_AtBlock(i_Block)) CYCLE
          ib = ib + 1
          jb = 0
          DO j=1,nb_var
            IF (sym(j) /= Ind_Coord_AtBlock(i_Block)) CYCLE
            jb = jb + 1
            d0h_PerBlock(ib,jb) = d0h(i,j)
            d0k_PerBlock(ib,jb) = d0k(i,j)
          END DO
        END DO

        !d0c_ini => d0c_ini_PerBlock
        ib = 0
        DO i=1,nb_var
          IF (sym(i) /= Ind_Coord_AtBlock(i_Block)) CYCLE
          ib = ib + 1

          d0c_ini_PerBlock(ib,:) = d0c_ini(i,iNM+1:iNM+nb_NM)

        END DO

        !freq
        CALL calc_freq(nb_NM,d0h_PerBlock,d0k_PerBlock,d0eh_PerBlock,   &
                       d0c_PerBlock,d0c_inv_PerBlock,norm_PerBlock,     &
                       d0c_ini_PerBlock,diab_freq)

        norme = norme * norm_PerBlock

        !d0c_ini_PerBlock => d0c_ini
        !d0c_PerBlock => d0c
        !d0c_inv_PerBlock => d0c_inv
        !d0eh_PerBlock => d0eh
        ib = 0
        DO i=1,nb_var
          IF (sym(i) /= Ind_Coord_AtBlock(i_Block)) CYCLE
          ib = ib + 1

          d0c_ini(i,iNM+1:iNM+nb_NM) = d0c_ini_PerBlock(ib,:)
          d0c(    i,iNM+1:iNM+nb_NM) = d0c_PerBlock(    ib,:)
          d0c_inv(iNM+1:iNM+nb_NM,i) = d0c_inv_PerBlock(:,ib)
        END DO
        d0eh(iNM+1:iNM+nb_NM) = d0eh_PerBlock(:)


        CALL dealloc_NParray(d0k_PerBlock,      "d0k_PerBlock",      name_sub)
        CALL dealloc_NParray(d0h_PerBlock,      "d0h_PerBlock",      name_sub)
        CALL dealloc_NParray(d0c_PerBlock,      "d0c_PerBlock",      name_sub)
        CALL dealloc_NParray(d0c_inv_PerBlock,  "d0c_inv_PerBlock",  name_sub)
        CALL dealloc_NParray(d0c_ini_PerBlock,  "d0c_ini_PerBlock",  name_sub)

        CALL dealloc_NParray(d0eh_PerBlock,     "d0eh_PerBlock",     name_sub)

        IF (debug) THEN
          write(out_unitp,*) '========================================='
          write(out_unitp,*) '===== END Block: ',i_Block
          write(out_unitp,*) '========================================='
          flush(out_unitp)
        END IF


        iNM = iNM + nb_NM
      END DO


      CALL dealloc_NParray(nb_PerBlock,"nb_PerBlock",name_sub)
      CALL dealloc_NParray(Ind_Coord_AtBlock,"Ind_Coord_AtBlock",name_sub)
      CALL dealloc_NParray(Ind_Coord_PerBlock,"Ind_Coord_PerBlock",name_sub)

      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'determinants :',norme

        DO i=1,nb_var,3
          i2 = min(i+2,nb_var)
          write(out_unitp,'("frequencies (cm-1): ",i0,"-",i0,3(1x,f0.4))') &
                            i,i2,d0eh(i:i2)*get_Conv_au_TO_unit('E','cm-1')
        END DO

        write(out_unitp,*)
        write(out_unitp,*)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF

!     -----------------------------------------------------------------

      END SUBROUTINE calc_freq_block



      SUBROUTINE calc_freq_WITH_d0c(nb_var,d0h,d0k,d0eh,                &
                                    d0c,d0c_inv,norme)
      IMPLICIT NONE

      integer :: nb_var
      logical :: diab_freq

      real (kind=Rkind) :: d0h(nb_var,nb_var)
      real (kind=Rkind) :: d0k(nb_var,nb_var)

      real (kind=Rkind) :: d0eh(nb_var)


      real (kind=Rkind) :: d0c(nb_var,nb_var)
      real (kind=Rkind) :: d0c_inv(nb_var,nb_var)

      real (kind=Rkind) :: mat1(nb_var,nb_var)
      real (kind=Rkind) :: mat2(nb_var,nb_var)

      real (kind=Rkind) :: val,val1
!     real (kind=Rkind),parameter :: epsi_freq = ONETENTH**7
      real (kind=Rkind),parameter :: epsi_freq = ONETENTH**10

!----- pour le determinant de d0c
      real (kind=Rkind) ::    d,norme,max_err
      real (kind=Rkind) ::    dh,dk


      integer :: i,j,k,l
      integer :: ierr

!-----------------------------------------------------------
       logical, parameter :: debug = .FALSE.
       !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING calc_freq_WITH_d0c'
         write(out_unitp,*) 'd0h',nb_var
         CALL Write_Mat(d0h,out_unitp,5)
         write(out_unitp,*) 'd0k',nb_var
         CALL Write_Mat(d0k,out_unitp,5)
         write(out_unitp,*) 'd0c',nb_var
         CALL Write_Mat(d0c,out_unitp,5)
         flush(out_unitp)
       END IF
!-----------------------------------------------------------

!-----------------------------------------------------------
!-----------------------------------------------------------
!      calcul le determinant de d0c :
!
!      3 facons :
!
!      1)  on calcule directement d=det(d0c)
!      2)  on utilise les valeurs propres d0eh et d0ek
!      3)  d = [det(d0h)/det(d0k)]^1/4


!
!      - cas 1 --------------------------------
       CALL Det_OF_m1(d0c,d,nb_var)
       write(out_unitp,*) 'det d0c',d
       norme = sqrt(abs(d))
!      ----------------------------------------
!-----------------------------------------------------------
!-----------------------------------------------------------


!-----------------------------------------------------------
!-----------------------------------------------------------
       ! d0c_inv directly from d0c
       CALL inv_m1_TO_m2(d0c,d0c_inv,nb_var,0,ZERO) ! not SVD
!-----------------------------------------------------------
!-----------------------------------------------------------


      mat1 = matmul(d0c_inv,matmul(d0h,transpose(d0c_inv)))
      !write(out_unitp,*) 'diago hess ?'
      !CALL Write_Mat(mat1,out_unitp,5)

       DO i=1,nb_var
         d0eh(i) = mat1(i,i)
         IF (d0eh(i) < ZERO) THEN
           write(out_unitp,*) ' ERROR : one imaginary frequency',               &
                                    d0eh(i)*get_Conv_au_TO_unit('E','cm-1')
         END IF
       END DO

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'ZPE (cm-1): ',HALF*sum(d0eh(:))*get_Conv_au_TO_unit('E','cm-1')
         write(out_unitp,*) 'ZPE   (au): ',HALF*sum(d0eh(:))

         write(out_unitp,*) 'frequencies (cm-1): ',d0eh(:)*get_Conv_au_TO_unit('E','cm-1')
         flush(out_unitp)
       END IF
!-----------------------------------------------------------


!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'determinants :',norme
         write(out_unitp,*) 'END calc_freq_WITH_d0c'
         flush(out_unitp)
       END IF
!-----------------------------------------------------------


      end subroutine calc_freq_WITH_d0c
!================================================================
!    SUBROUTINE H0 symmetrization
!================================================================
      SUBROUTINE H0_symmetrization(h,n,sym,dim_equi,tab_equi)
      IMPLICIT NONE

       integer       :: n
       real (kind=Rkind) :: h(n,n)
       real (kind=Rkind) :: h2(n,n)
       integer       :: sym(n)
       integer       :: dim_equi(n),tab_equi(n,n)

       integer       :: i,j,k,l,leq,keq

!-----------------------------------------------------------
       logical, parameter :: debug = .FALSE.
       !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING H0_symmetrization'
         write(out_unitp,*) 'h: ',n
         CALL Write_Mat(h,out_unitp,5)
         write(out_unitp,*)
         write(out_unitp,*) 'sym',sym
         write(out_unitp,*) 'eq'
         DO i=1,n
          write(out_unitp,*) 'tab_equi:',i,dim_equi(i),':',                     &
                       (tab_equi(i,k),k=1,dim_equi(i))
        END DO
       END IF
!-----------------------------------------------------------


!---- h becomes block-diagonal ----------------
      DO i=1,n
      DO j=i+1,n
        IF (sym(i) /= sym(j)) THEN
          h(i,j) = ZERO
          h(j,i) = ZERO
        END IF
      END DO
      END DO

!----- if sym=-1 => h(i,i) = 1 and h(i,j) = 0
      DO i=1,n
        IF (sym(i) == -1) THEN
          h(i,:) = ZERO
          h(:,i) = ZERO
          h(i,i) = ONE
        END IF
      END DO

!----- if sym=-2 => h(i,i) = -1 and h(i,j) = 0
      DO i=1,n
        IF (sym(i) == -2) THEN
          h(i,:) = ZERO
          h(:,i) = ZERO
          h(i,i) = -ONE
        END IF
      END DO

      h2(:,:) = h(:,:)

!---- for the equivalence --------------------------
      h(:,:) = ZERO
      DO i=1,n
        DO k=1,dim_equi(i)
          keq = tab_equi(i,k)
          h(i,i) = h(i,i) + h2(keq,keq)
        END DO
        h(i,i) = h(i,i) / real(dim_equi(i),kind=Rkind)
      END DO

      DO i=1,n
      DO j=i+1,n
        IF (sym(i) == sym(j)) THEN
          DO k=1,dim_equi(i)
          DO l=1,dim_equi(j)
            keq = tab_equi(i,k)
            leq = tab_equi(j,l)
            IF (sym(keq) == sym(leq)) h(i,j) = h(i,j) + h2(keq,leq)
          END DO
          END DO
          h(i,j) = h(i,j) / real(dim_equi(i),kind=Rkind)
          h(j,i) = h(i,j)
        END IF
      END DO
      END DO


!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'h: ',n
         CALL Write_Mat(h,out_unitp,5)
         write(out_unitp,*) 'END H0_symmetrization'
       END IF
!-----------------------------------------------------------

       RETURN
       end subroutine H0_symmetrization
!=====================================================================
!
! ++   calcule les frequences et les modes normaux
!
!=====================================================================
!
      SUBROUTINE calc_freqNM(nb_var,d0h,d0sm,d0eh,d0c,d0c_inv,print)
      IMPLICIT NONE

      integer :: nb_var

      real (kind=Rkind) :: d0h(nb_var,nb_var)
      real (kind=Rkind) :: d0sm(nb_var)
      real (kind=Rkind) :: d0eh(nb_var)
      real (kind=Rkind) :: d0c(nb_var,nb_var)
      real (kind=Rkind) :: d0c_inv(nb_var,nb_var)

      real (kind=Rkind) :: d0ek(nb_var)
      real (kind=Rkind) :: d0k(nb_var,nb_var)
      real (kind=Rkind) :: d0ch(nb_var,nb_var)
      real (kind=Rkind) :: mat1(nb_var,nb_var)
      real (kind=Rkind) :: mat2(nb_var,nb_var)

      real (kind=Rkind) :: val,val1

!----- pour le determinant de d0c
      real (kind=Rkind) ::    d
      real (kind=Rkind) ::    dh,dk


      logical :: print


      integer :: i,j,k,l

!-----------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING calc_freqNM'
         write(out_unitp,*) 'd0h',nb_var
         CALL Write_Mat(d0h,out_unitp,5)
         write(out_unitp,*) 'sqrt(masses)',d0sm
       END IF
!-----------------------------------------------------------

!-----------------------------------------------------------
!----- mass-weighted hessian -------------------------------
       DO i=1,nb_var
       DO j=1,nb_var
         d0k(i,j) = d0h(i,j)/(d0sm(i)*d0sm(j))
       END DO
       END DO

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*)
         write(out_unitp,*) 'mass-weighted hessian:'
         CALL Write_Mat(d0k,out_unitp,5)
       END IF
!-----------------------------------------------------------

       CALL diagonalization(d0k,d0eh,d0ch,nb_var,2,2,.TRUE.)

       DO i=1,nb_var
         IF (d0eh(i) >= ZERO) THEN
           d0eh(i) =  sqrt(d0eh(i))
         ELSE
           d0eh(i) = -sqrt(-d0eh(i))
!          d0eh(i) =  sqrt(-d0eh(i))
!          write(out_unitp,*) ' ERROR : one imaginary frequency',d0eh(i),'au'
         END IF

       END DO


!-----------------------------------------------------------
       IF (debug .OR. print) THEN
         write(out_unitp,*) 'frequencies : ',d0eh(:)*get_Conv_au_TO_unit('E','cm-1')
         write(out_unitp,*)
         write(out_unitp,*) 'modes normaux'
         CALL Write_Mat(d0ch,out_unitp,5)
       END IF
!-----------------------------------------------------------
!      mat1 = transpose(d0ch)
!      mat2 = matmul(mat1,d0ch)
!      write(out_unitp,*) 'test inversion de d0ch'
!      CALL Write_Mat(mat2,out_unitp,5)





!-----------------------------------------------------------
!----- calcul de la matrice qui passe des modes normaux ----
!      exprimes en fonction des coordonnees
!      passage de la base b2 a b0
!-----------------------------------------------------------

       DO i=1,nb_var
       DO k=1,nb_var
         d0c(k,i) = d0ch(k,i) * d0sm(k)
         d0c_inv(i,k) = d0ch(k,i) / d0sm(k)
       END DO
       END DO



!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'matrice de passage de b2 a b0'
         CALL Write_Mat(d0c,out_unitp,5)
         write(out_unitp,*) 'matrice de passage de b0 a b2'
         CALL Write_Mat(d0c_inv,out_unitp,5)

!        test inversion ------------------------------------
         mat1 = matmul(d0c,d0c_inv)
         write(out_unitp,*) 'test inversion de d0c'
         CALL Write_Mat(mat1,out_unitp,5)
       END IF
!-----------------------------------------------------------

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'END calc_freqNM'
       END IF
!-----------------------------------------------------------

      RETURN
      end subroutine calc_freqNM
!=====================================================================
!
! ++   calculation gaussian_width and freq with curvilinear coordinates
!
!=====================================================================
!
      SUBROUTINE calc_freq_width(nb_var,A,d0c,d0eh,d0h,d0k)
      IMPLICIT NONE

      integer :: nb_var
      real (kind=Rkind) :: A(nb_var,nb_var)
      real (kind=Rkind) :: d0h(nb_var,nb_var),d0k(nb_var,nb_var)
      real (kind=Rkind) :: d0c(nb_var,nb_var)
      real (kind=Rkind) :: d0eh(nb_var)


      real (kind=Rkind), allocatable :: d0c_inv(:,:),d0c_ini(:,:)
      real (kind=Rkind), allocatable :: d0k_save(:,:)
      real (kind=Rkind) :: norme
      integer :: err_mem,memory

      CALL alloc_NParray(d0c_inv,[nb_var,nb_var],"d0c_inv","calc_freq_width")
      CALL alloc_NParray(d0c_ini,[nb_var,nb_var],"d0c_ini","calc_freq_width")
      CALL alloc_NParray(d0k_save,[nb_var,nb_var],"d0k_save","calc_freq_width")
      d0c_ini(:,:)  = ZERO
      d0k_save(:,:) = d0k(:,:)

      CALL calc_freq(nb_var,d0h,d0k_save,d0eh,                          &
                     d0c,d0c_inv,norme,d0c_ini,.FALSE.)

      CALL dealloc_NParray(d0c_inv,"d0c_inv","calc_freq_width")
      CALL dealloc_NParray(d0c_ini,"d0c_ini","calc_freq_width")
      CALL dealloc_NParray(d0k_save,"d0k_save","calc_freq_width")

      CALL gaussian_width(nb_var,A,d0c)


      end subroutine calc_freq_width


!=====================================================================
!
! ++   calculation gaussian_width
!
!=====================================================================
!
      SUBROUTINE gaussian_width(nb_var,A,d0c)
      IMPLICIT NONE


      integer :: nb_var
      real (kind=Rkind) :: A(nb_var,nb_var)
      real (kind=Rkind) :: d0c(nb_var,nb_var)
      real (kind=Rkind) :: trav1(nb_var),trav2(nb_var)

      real (kind=Rkind) :: val,val1

!----- pour le determinant de d0c
      real (kind=Rkind) ::    d,norme
      real (kind=Rkind) ::    dh,dk


      integer :: i,j,k

!-----------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING gaussian_width'
         write(out_unitp,*) 'd0c',nb_var
         CALL Write_Mat(d0c,out_unitp,5)
       END IF

       DO i=1,nb_var
       DO j=1,nb_var
         A(i,j) = ZERO
         DO k=1,nb_var
           A(i,j) = A(i,j) -HALF*d0c(i,k)*d0c(j,k)
         END DO
       END DO
       END DO


       IF (debug) THEN
         write(out_unitp,*) 'A',nb_var
         CALL Write_Mat(A,out_unitp,5)
         write(out_unitp,*) 'END gaussian_width'
       END IF
       end subroutine gaussian_width

!=====================================================================
!
! ++   if v(i) and v(j), rotation of vector c(.,i) and c(.,j)
!      such that ???
!
!      it is working only if 2 vectors are degenerated !!!!
!
!=====================================================================
!
      SUBROUTINE rota_denerated(v,c,n)
      IMPLICIT NONE

      integer       :: n
      real (kind=Rkind) :: v(n)
      real (kind=Rkind) :: c(n,n)

      integer       :: i,j,k,ind_maxi,ind_maxj
      real (kind=Rkind) :: max_ci,max_cj,norme,ai,aj,cc,ss

      real (kind=Rkind), parameter :: epsi = ONETENTH**8

!---------------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unitp,*) 'BEGINNING rota_denerated'
      write(out_unitp,*) 'v',v
      write(out_unitp,*) 'c'
      CALL Write_Mat(c,out_unitp,5)
      END IF
!---------------------------------------------------------------------

      DO i=1,n
!       swap sign of c(.,i) such max_ci is positive
        max_ci = ZERO
        DO k=1,n
          IF ( abs(c(k,i)) > max_ci ) THEN
            max_ci = abs(c(k,i))
            ind_maxi = k
          END IF
        END DO
        IF (c(ind_maxi,i) < ZERO) c(:,i) = -c(:,i)
      END DO

      DO i=1,n-1

        IF ( abs(v(i)-v(i+1)) < epsi) THEN

          j = i+1
!         determination of max_ci and max_cj
          max_ci = ZERO
          DO k=1,n
            IF ( abs(c(k,i)) > max_ci ) THEN
              max_ci = abs(c(k,i))
              ind_maxi = k
            END IF
          END DO
          max_cj = ZERO
          DO k=1,n
            IF ( abs(c(k,j)) > max_cj .AND. k /= ind_maxi) THEN
              max_cj = abs(c(k,j))
              ind_maxj = k
            END IF
          END DO

!         rotation of c(.i) and c(.,j)
          norme = sqrt(max_ci*max_ci + c(ind_maxj,i)*c(ind_maxj,i))
          cc =  max_ci/norme
          ss = -c(ind_maxj,i)/norme

          DO k=1,n
           ai = c(k,i)
           aj = c(k,j)

           c(k,i) =  cc * ai + ss * aj
           c(k,j) = -ss * ai + cc * aj

          END DO

        END IF
      END DO



      DO i=1,n
!       swap sign of c(.,i) such max_ci is positive
        max_ci = ZERO
        DO k=1,n
          IF ( abs(c(k,i)) > max_ci ) THEN
            max_ci = abs(c(k,i))
            ind_maxi = k
          END IF
        END DO
        IF (c(ind_maxi,i) < ZERO) c(:,i) = -c(:,i)


      END DO


!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unitp,*) 'i et i+1 dege',j-1
      write(out_unitp,*) 'ind_maxi ind_maxj',ind_maxi,ind_maxj
      write(out_unitp,*) 'cos et sin',cc,ss,norme
      write(out_unitp,*) 'new c'
      CALL Write_Mat(c,out_unitp,5)
      write(out_unitp,*) 'END rota_denerated'
      END IF
!---------------------------------------------------------------------

      RETURN
      end subroutine rota_denerated
      SUBROUTINE rota_degenerate_opt1(v,c,n,degenerate_freq)
      IMPLICIT NONE

      integer           :: n
      real (kind=Rkind) :: v(n)
      real (kind=Rkind) :: c(n,n)
      TYPE (degenerate_freq_t), INTENT(IN) :: degenerate_freq

      integer           :: i,j,k,ind_maxi,ind_maxj
      real (kind=Rkind) :: max_ci,max_cj,norme,ai,aj,cc,ss


!---------------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unitp,*) 'BEGINNING rota_degenerate_opt1'
      write(out_unitp,*) 'v',v
      write(out_unitp,*) 'c'
      CALL Write_Mat(c,out_unitp,5)
      END IF
!---------------------------------------------------------------------

      DO i=1,n
        !swap sign of c(.,i) such max_ci is positive
        max_ci = ZERO
        DO k=1,n
          IF ( abs(c(k,i)) > max_ci ) THEN
            max_ci = abs(c(k,i))
            ind_maxi = k
          END IF
        END DO
        IF (c(ind_maxi,i) < ZERO) c(:,i) = -c(:,i)
      END DO

      DO i=1,n-1

        IF ( abs(v(i)-v(i+1)) < degenerate_freq%epsi) THEN

          j = i+1
          !determination of max_ci and max_cj
          max_ci   = ZERO
          ind_maxi = 0
          DO k=1,n
            IF (.NOT. degenerate_freq%list_k_check(k)) CYCLE
            IF ( abs(c(k,i)) > max_ci ) THEN
              max_ci   = abs(c(k,i))
              ind_maxi = k
            END IF
          END DO
          max_cj   = ZERO
          ind_maxj = 0
          DO k=1,n
            IF (.NOT. degenerate_freq%list_k_check(k)) CYCLE
            IF ( abs(c(k,j)) > max_cj .AND. k /= ind_maxi) THEN
              max_cj   = abs(c(k,j))
              ind_maxj = k
            END IF
          END DO

          IF (ind_maxj > 0 .AND. ind_maxi > 0) THEN
            !rotation of c(.i) and c(.,j)
            norme = sqrt(max_ci*max_ci + c(ind_maxj,i)*c(ind_maxj,i))
            cc =  max_ci/norme
            ss = -c(ind_maxj,i)/norme
            IF (debug) then
              write(out_unitp,*) 'i et i+1 dege',i,i+1
              write(out_unitp,*) 'ind_maxi ind_maxj',ind_maxi,ind_maxj
              write(out_unitp,*) 'cos et sin',cc,ss,norme
            END IF

            DO k=1,n
              ai = c(k,i)
              aj = c(k,j)

              c(k,i) =  cc * ai + ss * aj
              c(k,j) = -ss * ai + cc * aj

            END DO
          END IF

        END IF
      END DO



      DO i=1,n
        !swap sign of c(.,i) such max_ci is positive
        max_ci = ZERO
        DO k=1,n
          IF ( abs(c(k,i)) > max_ci ) THEN
            max_ci = abs(c(k,i))
            ind_maxi = k
          END IF
        END DO
        IF (c(ind_maxi,i) < ZERO) c(:,i) = -c(:,i)


      END DO


!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unitp,*) 'new c'
      CALL Write_Mat(c,out_unitp,5)
      write(out_unitp,*) 'END rota_degenerate_opt1'
      END IF
!---------------------------------------------------------------------

    end subroutine rota_degenerate_opt1
!=====================================================================
!
! ++   order the vectors cf(.,i) such that cf(.,i) is "closed" to ci(.,i)
!      order also vi(i) and vf(i)
!
!=====================================================================
!
      SUBROUTINE order_ini4(c,c_inv,e,c_ini,n,dia_freq)
      IMPLICIT NONE

      logical       :: dia_freq
      integer       :: n
      real (kind=Rkind) :: c(n,n),c_ini(n,n),c_inv(n,n)
      real (kind=Rkind) :: e(n)

      integer       :: i,j,k,ind_maxi
      real (kind=Rkind) :: de,a,Sii,Sjj,Sij,max_Si,tab_Sji(n)

      !---------------------------------------------------------------------
      character (len=*), parameter :: name_sub='order_ini4'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      !---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) '  dia_freq',dia_freq
        write(out_unitp,*) '  e',e
       write(out_unitp,*) '  c'
       CALL Write_Mat(c,out_unitp,5)
       write(out_unitp,*) '  c_inv'
       CALL Write_Mat(c_inv,out_unitp,5)
        write(out_unitp,*) '  c_ini'
        CALL Write_Mat(c_ini,out_unitp,5)
      END IF
!---------------------------------------------------------------------

      IF ( sum(abs(c_ini)) /= ZERO ) THEN

        IF (dia_freq) THEN

          DO i=1,n

            max_Si = ZERO
            DO j=i,n
              !  calculation of Sij
              Sij = dot_product(c_ini(:,i) , c(:,j))
              Sii = dot_product(c_ini(:,i) , c_ini(:,i))
              Sjj = dot_product(c(:,j)     , c(:,j))

              Sij = Sij/sqrt(Sii*Sjj)
              IF ( abs(Sij) > max_Si ) THEN
                max_Si = abs(Sij)
                ind_maxi = j
                de = abs(  sqrt(abs(e(i))) -sqrt(abs(e(j))) )
              END IF
            END DO
            !write(out_unitp,*) 'max_Si',i,ind_maxi,max_Si

            !IF (ind_maxi /= i .AND. de < TWO*ONETENTH**4) THEN
            IF (ind_maxi /= i) THEN
              IF (debug) write(out_unitp,*) 'permutation: max_Si',i,ind_maxi,max_Si

              ! permutation of the vectors c(.,i)     and c(.,ind_maxi)
              ! permutation of the vectors c_inv(i,.) and c_inv(ind_maxi,.)
              ! and e(i) and e(ind_maxi)
              a           = e(i)
              e(i)        = e(ind_maxi)
              e(ind_maxi) = a
              DO k=1,n
                a             = c(k,i)
                c(k,i)        = c(k,ind_maxi)
                c(k,ind_maxi) = a

                a                 = c_inv(i,k)
                c_inv(i,k)        = c_inv(ind_maxi,k)
                c_inv(ind_maxi,k) = a
              END DO
            END IF
          END DO
        END IF

        !-- phase changement -------------------------
        DO i=1,n
          ! calculation of Sii
          Sii = dot_product(c_ini(:,i) , c(:,i))
          IF (debug) write(out_unitp,*) 'Sii',i,Sii

          ! if Sii<0 => c(.,i)=-c(.,i) and c_inv(i,.) = -c_inv(i,.)
          IF (Sii < ZERO ) THEN
            c(:,i)     = -c(:,i)
            c_inv(i,:) = -c_inv(i,:)
          END IF
        END DO

        ! check the overlap with c_ini
        IF (dia_freq) THEN

          IF (debug) write(out_unitp,*) ' Overlapp between c(:,i) and c_ini(:,i)'
          DO i=1,n


            max_Si = ZERO
            DO j=1,n
              !  calculation of Sij
              Sij = dot_product(c_ini(:,i) , c(:,j))
              Sii = dot_product(c_ini(:,i) , c_ini(:,i))  ! should be one
              Sjj = dot_product(c(:,j)     , c(:,j))      ! should be one

              tab_Sji(j) = Sij/sqrt(Sii*Sjj)

            END DO
            !write(out_unitp,'(a,i0,50(xf5.2))') 'Error on Sji: ',i,tab_Sji
            tab_Sji(i) = tab_Sji(i) - ONE

            IF (debug) &
              write(out_unitp,'(a,i0,50(1x,f5.2))') 'Sji and Error on Sji: ',i,tab_Sji(i)+ONE,maxval(abs(tab_Sji))

          END DO
        END IF



        !IF (Grid_maxth == 1) c_ini(:,:) = c(:,:)   ! omp_get_max_threads
        c_ini(:,:) = c(:,:)
      ELSE
        IF (debug) write(out_unitp,*) 'Initialization of C_ini in ',name_sub
        c_ini(:,:) = c(:,:)
      END IF

      !---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) '  e',e
        write(out_unitp,*) '  new c'
        CALL Write_Mat(c,out_unitp,5)
        write(out_unitp,*) '  new c_inv'
        CALL Write_Mat(c_inv,out_unitp,5)
        write(out_unitp,*) 'END ',name_sub
      END IF
      !---------------------------------------------------------------------

      END SUBROUTINE order_ini4
      SUBROUTINE order_ini5(c,c_inv,e,c_ini,n,dia_freq)
      IMPLICIT NONE

      logical       :: dia_freq
      integer       :: n
      real (kind=Rkind) :: c(n,n),c_ini(n,n),c_inv(n,n)
      real (kind=Rkind) :: e(n)

      integer       :: i,j,k,ind_maxi
      real (kind=Rkind) :: Nini,Ni,de,a,Sii,Sjj,Sij,max_Si,tab_Sji(n)

      !---------------------------------------------------------------------
      character (len=*), parameter :: name_sub='order_ini5'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      !---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) '  dia_freq',dia_freq
        write(out_unitp,*) '  e',e
       write(out_unitp,*) '  c'
       CALL Write_Mat(c,out_unitp,5)
       write(out_unitp,*) '  c_inv'
       CALL Write_Mat(c_inv,out_unitp,5)
        write(out_unitp,*) '  c_ini'
        CALL Write_Mat(c_ini,out_unitp,5)
      END IF
!---------------------------------------------------------------------

      IF ( sum(abs(c_ini)) /= ZERO ) THEN

        IF (dia_freq) THEN

          DO i=1,n

            max_Si = ZERO
            DO j=i,n
              !  calculation of Sij
              Sij = dot_product(c_ini(:,i) , c(:,j))
              Sii = dot_product(c_ini(:,i) , c_ini(:,i))
              Sjj = dot_product(c(:,j)     , c(:,j))

              Sij = Sij/sqrt(Sii*Sjj)
              IF ( abs(Sij) > max_Si ) THEN
                max_Si = abs(Sij)
                ind_maxi = j
                de = abs(  sqrt(abs(e(i))) -sqrt(abs(e(j))) )
              END IF
            END DO
            !write(out_unitp,*) 'max_Si',i,ind_maxi,max_Si

            !IF (ind_maxi /= i .AND. de < TWO*ONETENTH**4) THEN
            IF (ind_maxi /= i) THEN
              IF (debug) write(out_unitp,*) 'permutation: max_Si',i,ind_maxi,max_Si

              ! permutation of the vectors c(.,i)     and c(.,ind_maxi)
              ! permutation of the vectors c_inv(i,.) and c_inv(ind_maxi,.)
              ! and e(i) and e(ind_maxi)
              a           = e(i)
              e(i)        = e(ind_maxi)
              e(ind_maxi) = a
              DO k=1,n
                a             = c(k,i)
                c(k,i)        = c(k,ind_maxi)
                c(k,ind_maxi) = a

                a                 = c_inv(i,k)
                c_inv(i,k)        = c_inv(ind_maxi,k)
                c_inv(ind_maxi,k) = a
              END DO
            END IF
          END DO
        END IF

        !-- phase changement -------------------------
        DO i=1,n
          ! calculation of Sii
          Nini = sqrt(dot_product(c_ini(:,i) , c_ini(:,i)))
          Ni = sqrt(dot_product(c(:,i) , c(:,i)))
          Sii = dot_product(c_ini(:,i) , c(:,i))/(Nini*Ni)
          IF (debug) write(out_unitp,*) 'Sii',i,Sii

          ! if Sii<0 => c(.,i)=-c(.,i) and c_inv(i,.) = -c_inv(i,.)
          IF (Sii < ZERO ) THEN
            c(:,i)     = -c(:,i)
            c_inv(i,:) = -c_inv(i,:)
          END IF
        END DO

        ! check the overlap with c_ini
        IF (dia_freq) THEN

          IF (debug) write(out_unitp,*) ' Overlapp between c(:,i) and c_ini(:,i)'
          DO i=1,n
            DO j=1,n
              !  calculation of Sij
              Sij = dot_product(c_ini(:,i) , c(:,j))
              Sii = dot_product(c_ini(:,i) , c_ini(:,i))  ! should be one
              Sjj = dot_product(c(:,j)     , c(:,j))      ! should be one

              tab_Sji(j) = Sij/sqrt(Sii*Sjj)

            END DO
            !write(out_unitp,'(a,i0,50(xf5.2))') 'Error on Sji: ',i,tab_Sji
            tab_Sji(i) = tab_Sji(i) - ONE
            IF (debug) write(out_unitp,'(3a,i0,2(1x,f5.2))') 'In ',name_sub,               &
               'Sji(i) and Error on Sji: ',i,tab_Sji(i)+ONE,maxval(abs(tab_Sji))

            IF (debug) THEN
              tab_Sji(i) = tab_Sji(i) + ONE
              write(out_unitp,'(a,i0,50(1x,f5.2))') 'Sji(:): ',i,tab_Sji(:)
            END IF
          END DO
        END IF



        !IF (Grid_maxth == 1) c_ini(:,:) = c(:,:)   ! omp_get_max_threads
        c_ini(:,:) = c(:,:)
      ELSE
        IF (debug) write(out_unitp,*) 'Initialization of C_ini in ',name_sub
        c_ini(:,:) = c(:,:)
      END IF

      !---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) '  e',e
        write(out_unitp,*) '  new c'
        CALL Write_Mat(c,out_unitp,5)
        write(out_unitp,*) '  new c_inv'
        CALL Write_Mat(c_inv,out_unitp,5)
        write(out_unitp,*) 'END ',name_sub
      END IF
      !---------------------------------------------------------------------

      END SUBROUTINE order_ini5
!=====================================================================
!
! ++   order the vectors cf(.,i) such that cf(.,i) is "closed" to ci(.,i)
!      order also vi(i) and vf(i)
!
!=====================================================================
!
      SUBROUTINE sort_with_Tab(c,c_inv,e,tab_order,n)
      IMPLICIT NONE

      integer           :: n
      real (kind=Rkind) :: c(n,n),c_inv(n,n)
      real (kind=Rkind) :: e(n)
      real (kind=Rkind) :: tab_order(n)


      real (kind=Rkind) :: tab_order_loc(n),Max_val
      integer           :: i,j,k
      real (kind=Rkind) :: a


!---------------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unitp,*) 'BEGINNING sort_with_Tab'
      write(out_unitp,*) 'e',e
      write(out_unitp,*) 'c, c_inv'
      CALL Write_Mat(c,out_unitp,5)
      CALL Write_Mat(c_inv,out_unitp,5)
      END IF
!---------------------------------------------------------------------



      ! first initialization of tab_order_loc
      tab_order_loc(:) = tab_order(:)
      Max_val = maxval(tab_order_loc)
      IF (minval(tab_order) == -ONE) THEN
        DO i=1,n
          IF (tab_order(i) == -ONE) tab_order_loc(i) = Max_val + ONE
        END DO
      END IF

      ! sort with tab_order_loc
      DO i=1,n
      DO j=i+1,n
       IF (tab_order_loc(i) > tab_order_loc(j)) THEN
         ! permutation
         a                = tab_order_loc(i)
         tab_order_loc(i) = tab_order_loc(j)
         tab_order_loc(j) = a

         a           = e(i)
         e(i)        = e(j)
         e(j)        = a
         DO k=1,n
           a             = c(k,i)
           c(k,i)        = c(k,j)
           c(k,j)        = a

           a                 = c_inv(i,k)
           c_inv(i,k)        = c_inv(j,k)
           c_inv(j,k)        = a
          END DO

        END IF
      END DO
      END DO

!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unitp,*) 'e',e
      write(out_unitp,*) 'new c and c_inv'
      CALL Write_Mat(c,out_unitp,5)
      CALL Write_Mat(c_inv,out_unitp,5)
      write(out_unitp,*) 'END sort_with_Tab'
      END IF
!---------------------------------------------------------------------

      end subroutine sort_with_Tab
END MODULE mod_freq
