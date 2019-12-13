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
!      with contributions of Mamadou Ndong
!
!===========================================================================
!===========================================================================
MODULE mod_dnRho
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: sub3_dnrho, sub3_dnrho_ana ,Write_rho

  CONTAINS
!=====================================================================
!
! ++   fi Fij f0=rho calculation
!      fi(i)    = ( drho/dQi) / rho
!      Fij(i,j) = ( d2rho/dQidQj) / rho
!
! nrho = 0 : Normal condition : rho = sqrt(jac)
! nrho = 1 : Wilson condition : rho = 1
! nrho = 10: Wilson condition : rho = 1 (but without vep)
! nrho = 2 : ana              : analitical derivation
! nrho = 20: ana              : analitical derivation (but without vep)
! nrho = 3 : num              : numerical derivation
!
!=====================================================================
!
      SUBROUTINE sub3_dnrho(dnrho,dnjac,                                &
                            Qact,mole,nderiv,num,step,nrho)
      USE mod_system
      USE mod_dnSVM
      USE mod_Tnum
      IMPLICIT NONE

      TYPE(Type_dnS)    :: dnJac,dnrho

      TYPE (zmatrix), intent(in)    :: mole
      real (kind=Rkind), intent(inout) :: Qact(mole%nb_var)

      logical           :: num
      real (kind=Rkind) :: step
      integer           :: nderiv
      integer           :: nrho


      integer :: i



!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'sub3_dnrho'
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'step',step
         write(out_unitp,*)
         CALL Write_mole(mole)
         write(out_unitp,*)
       END IF
!-----------------------------------------------------------

       IF (nrho .EQ. 0) THEN
!      euclidien

         CALL sub_dnS1_TO_dnS2(dnJac,dnrho,nderiv)

       ELSE IF (nrho .EQ. 1 .OR. nrho .EQ. 10) THEN
!      wilson

         CALL Set_ZERO_TO_dnSVM(dnrho)
         dnrho%d0 = ONE

       ELSE IF (nrho .EQ. 3) THEN
!      numerical
         CALL    sub3_dnrho_num(dnrho,Qact,mole,nderiv,step)
       ELSE IF (nrho .EQ. 2) THEN
!      analitical
         CALL    sub3_dnrho_ana(dnrho,Qact,mole,2)
       ELSE IF (nrho .EQ. 20) THEN
!      analitical (Without vep)
         CALL    sub3_dnrho_ana(dnrho,Qact,mole,1)
       ELSE
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nrho =',nrho,' is not defined'
          STOP
       END IF

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'rho'
         CALL Write_dnSVM(dnrho)
         write(out_unitp,*) 'END ',name_sub
       END IF
!-----------------------------------------------------------

       end subroutine sub3_dnrho
!
!=====================================================================
!
! ++   fi Fij f0=rho calculation
!      fi(i)    = ( drho/dQi) / rho
!      Fij(i,j) = ( d2rho/dQidQj) / rho
!
!      num              : numerical derivation
!
!=====================================================================
!
      SUBROUTINE sub3_dnrho_num(dnrho,Qact,mole,nderiv,step)
      USE mod_system
      USE mod_dnSVM
      USE mod_Tnum
      !USE mod_paramQ
      IMPLICIT NONE

      TYPE (zmatrix), intent(in)    :: mole
      real (kind=Rkind), intent(inout) :: Qact(mole%nb_var)

      TYPE(Type_dnS)    :: dnrho
      TYPE(Type_dnS)    :: dnrho1,dnrho2

      logical           :: num
      real (kind=Rkind) :: step
      integer           :: nderiv
      integer           :: nrho


      real (kind=Rkind) :: Qacti,Qactj
      integer           :: i,j



!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'sub3_dnrho_num'
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         CALL Write_mole(mole)
         write(out_unitp,*)
       END IF
!-----------------------------------------------------------

!===================================================================
!
!       Calcul en Qact
!
!===================================================================

       CALL sub3_dnrho_ana(dnrho,Qact,mole,0)
!      -----------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'step',step
         write(out_unitp,*) 'Qact',Qact
         write(out_unitp,*) 'f0',dnrho%d0
       END IF
!      -----------------------------------------------------

!===================================================================
!
!       Calcul en Qact(i)+step
!           et en Qact(i)-step
!
!===================================================================

      CALL alloc_dnSVM(dnrho1,dnrho%nb_var_deriv,0)
      CALL alloc_dnSVM(dnrho2,dnrho%nb_var_deriv,0)

      IF (nderiv >= 1) THEN
       DO i=1,mole%nb_act

         Qacti = Qact(i)
         Qact(i) = Qacti + step

         CALL sub3_dnrho_ana(dnrho1,Qact,mole,0)
!        -----------------------------------------------------
         IF (debug) THEN
           write(out_unitp,*) 'Qact',Qact
           write(out_unitp,*) 'f1',dnrho1%d0
         END IF
!        -----------------------------------------------------

         Qact(i) = Qacti - step
         CALL sub3_dnrho_ana(dnrho2,Qact,mole,0)
!        -----------------------------------------------------
         IF (debug) THEN
           write(out_unitp,*) 'Qact',Qact
           write(out_unitp,*) 'f2',dnrho2%d0
         END IF
!        -----------------------------------------------------


         CALL d1d2(dnrho%d0,dnrho1%d0,dnrho2%d0,step)
         dnrho%d1(i)    = dnrho1%d0/dnrho%d0
         IF (nderiv == 2) dnrho%d2(i,i) = dnrho2%d0/dnrho%d0

         Qact(i) = Qacti

       END DO
      END IF


!===================================================================
!
!       Calcul en Qact(i)+/-step
!           et en Qact(j)+/-step
!
!===================================================================

      IF (nderiv == 2) THEN
       DO i=1,mole%nb_act
       DO j=i+1,mole%nb_act

         Qacti = Qact(i)
         Qactj = Qact(j)

         Qact(i) = Qacti + step
         Qact(j) = Qactj + step

         CALL sub3_dnrho_ana(dnrho1,Qact,mole,0)
         dnrho%d2(i,j) = dnrho1%d0
!        -----------------------------------------------------
         IF (debug) THEN
           write(out_unitp,*) 'Qact',Qact
           write(out_unitp,*) 'f1',dnrho1%d0,dnrho%d2(i,j)
         END IF
!        -----------------------------------------------------


         Qact(i) = Qacti - step
         Qact(j) = Qactj - step
         CALL sub3_dnrho_ana(dnrho1,Qact,mole,0)
         dnrho%d2(i,j) = dnrho%d2(i,j) + dnrho1%d0
!        -----------------------------------------------------
         IF (debug) THEN
           write(out_unitp,*) 'Qact',Qact
           write(out_unitp,*) 'f1',dnrho1%d0,dnrho%d2(i,j)
         END IF
!        -----------------------------------------------------


         Qact(i) = Qacti - step
         Qact(j) = Qactj + step
         CALL sub3_dnrho_ana(dnrho1,Qact,mole,0)
         dnrho%d2(i,j) = dnrho%d2(i,j) - dnrho1%d0
!        -----------------------------------------------------
         IF (debug) THEN
           write(out_unitp,*) 'Qact',Qact
           write(out_unitp,*) 'f1',dnrho1%d0,dnrho%d2(i,j)
         END IF
!        -----------------------------------------------------



         Qact(i) = Qacti + step
         Qact(j) = Qactj - step
         CALL sub3_dnrho_ana(dnrho1,Qact,mole,0)
         dnrho%d2(i,j) = dnrho%d2(i,j) - dnrho1%d0
!        -----------------------------------------------------
         IF (debug) THEN
           write(out_unitp,*) 'Qact',Qact
           write(out_unitp,*) 'f1',dnrho1%d0,dnrho%d2(i,j)
         END IF
!        -----------------------------------------------------



         dnrho%d2(i,j) = dnrho%d2(i,j)/(FOUR*dnrho%d0*step*step)
         dnrho%d2(j,j) = dnrho%d2(i,j)


         Qact(i) = Qacti
         Qact(j) = Qactj

       END DO
       END DO
      END IF


!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'rho'
         CALL Write_dnSVM(dnrho)
       write(out_unitp,*) 'END ',name_sub
       END IF
!-----------------------------------------------------------
      CALL dealloc_dnSVM(dnrho1)
      CALL dealloc_dnSVM(dnrho2)

      end subroutine sub3_dnrho_num
!
!=====================================================================
!
! ++   f0=rho= rho(Qact_1) * rho(Qact_2) ... rho(Qact_nb_act)
!
!      si nderiv = 0 on ne calcule pas les derivees
!      Rq : l'initialisation des fi et Fij doit ce faire APRES le test de nderiv
!
!=====================================================================
!
      SUBROUTINE sub3_dnrho_ana(dnrho,Qact,mole,nderiv)
      USE mod_system
      USE mod_dnSVM
      USE mod_Tnum
      !USE mod_paramQ
      IMPLICIT NONE

      TYPE (zmatrix), intent(in)    :: mole
      real (kind=Rkind), intent(inout) :: Qact(mole%nb_var)
      TYPE(Type_dnS), intent(inout) :: dnrho
      integer, intent(in)           :: nderiv

      integer                       :: iQact,jQact,i,dnErr
      integer                       :: type_act,iQact_transfo
      real (kind=Rkind)             :: d0tf,d1tf,d2tf,d3tf
      TYPE (Type_dnS)               :: dntf

      character (len=Name_longlen)  :: nom

!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
!      logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'sub3_dnrho_ana'
!-----------------------------------------------------------
      IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*)
         CALL Write_mole(mole)
         write(out_unitp,*)
      END IF
!-----------------------------------------------------------

       !write(6,*) 'mole%nrho_OF_Qact',mole%nrho_OF_Qact
!      - initialisation --------------------------------------------------
       CALL alloc_dnSVM(dntf,1,3)

       CALL Set_ZERO_TO_dnSVM(dnrho)
       dnrho%d0 = ONE

!      - rho calculation -------------------------------------------------
!      - diagonal terms ----
       DO iQact=1,mole%nb_act1

         iQact_transfo = mole%nrho_OF_Qact(iQact)

         IF (iQact_transfo == 0) THEN
           type_act = mole%tab_Qtransfo(mole%nb_Qtransfo)%type_Qin(iQact)

           IF ( type_act == 3 ) THEN
             iQact_transfo = 2
           ELSE
             iQact_transfo = 1
           END IF
         END IF
         !write(out_unitp,*) 'iQact_transfo,type_act',type_act,iQact_transfo

         CALL sub_dntf(iQact_transfo,dntf,Qact(iQact),(/ (ZERO,i=1,20) /), dnErr )
         IF (dnErr /= 0) THEN
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) '   ERROR in the sub_dntf call for the coordinates, iQact:',iQact
           STOP 'ERROR in sub_dntf called from sub3_dnrho_ana'
         END IF
         IF (iQact_transfo == 2) THEN
           dntf%d0 = -dntf%d0
           dntf%d1 = -dntf%d1
           dntf%d2 = -dntf%d2
           dntf%d3 = -dntf%d3
         END IF

         IF (nderiv == 0 .OR. iQact_transfo == 1) THEN
           dnrho%d0              = dnrho%d0       * dntf%d1(1)
         ELSE IF (nderiv == 1) THEN
           dnrho%d0              = dnrho%d0       * dntf%d1(1)
           dnrho%d1(iQact)       = dntf%d2(1,1)   / dntf%d1(1)
         ELSE IF (nderiv == 2) THEN
           dnrho%d0              = dnrho%d0       * dntf%d1(1)
           dnrho%d1(iQact)       = dntf%d2(1,1)   / dntf%d1(1)
           dnrho%d2(iQact,iQact) = dntf%d3(1,1,1) / dntf%d1(1)
         END IF
       END DO
!      - non diagonal terms ----
       IF (nderiv == 2) THEN
         DO iQact=1,mole%nb_act1
         DO jQact=iQact+1,mole%nb_act1
             dnrho%d2(iQact,jQact) = dnrho%d1(iQact) * dnrho%d1(jQact)
             dnrho%d2(jQact,iQact) = dnrho%d2(iQact,jQact)
         END DO
         END DO
       END IF
!      -------------------------------------------------------------------

       CALL dealloc_dnSVM(dntf)

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'rho'
         CALL Write_dnSVM(dnrho)
         write(out_unitp,*) 'END ',name_sub
       END IF
!-----------------------------------------------------------

       end subroutine sub3_dnrho_ana
       
      SUBROUTINE Write_rho(mole)
      USE mod_system
      USE mod_Tnum
      USE mod_MPI
      IMPLICIT NONE

      TYPE (zmatrix), intent(in)    :: mole

      integer                       :: iQact,type_act,iQact_transfo
      character (len=Name_longlen)  :: name_rho


!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
!      logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'Write_rho'
!-----------------------------------------------------------
      IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*)
         CALL Write_mole(mole)
         write(out_unitp,*)
      END IF
!-----------------------------------------------------------

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!      analysis of rho as a function of the zmatrix definition
       IF (print_level > -1 .AND. MPI_id==0) THEN

         write(out_unitp,*)
         write(out_unitp,*) '----------------------------------------------'
         write(out_unitp,*) ' Definition of the volume element dV=rho.dQact1.dQact2...'
         write(out_unitp,*) '    rho=rho(Qact1)*rho(Qact2)...'
         write(out_unitp,*) ' Remark: only for variables of type 1'
         write(out_unitp,'(a)') ' iQact iQact_transfo : rho(iQact)'
!
!
!        ------------------------------------------------------------------
!        - loop only on nb_act1 and not on nb_inact21 and nb_inact22  -----
!          because rho_i = 1 for those variables
!        ------------------------------------------------------------------
         DO iQact=1,mole%nb_act1

           iQact_transfo = mole%nrho_OF_Qact(iQact)

           IF (iQact_transfo == 0) THEN
             type_act = mole%tab_Qtransfo(mole%nb_Qtransfo)%type_Qin(iQact)

             IF ( type_act == 3 ) THEN
               iQact_transfo = 2
             ELSE
               iQact_transfo = 1
             END IF
           END IF


           IF (iQact_transfo == 2) THEN
             name_rho = 'sin(Qact)'
           ELSE
             name_rho = '1.'
           END IF
           IF(MPI_id==0) write(out_unitp,'(i6,8x,i6,a,a)') iQact,iQact_transfo,' : ',name_rho

         END DO
         IF(MPI_id==0) write(out_unitp,*) '----------------------------------------------'
         IF(MPI_id==0) write(out_unitp,*)

       END IF
      end subroutine Write_rho
END MODULE mod_dnRho
