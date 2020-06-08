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
      SUBROUTINE sub_mat6_HST(PrimOp,                                   &
                              d0MatHADAOp,nb_Op,nb_harm,                &
                              d0f_harm,                                 &
                              rho,                                      &
                              nb_inact2n,ind_quadra,                    &
                              d1xa,d2xaa,                               &
                              d0c,d1c,d2c,                              &
                              d0cd0c,                                   &
                              d1lnN,d2lnN,                              &
                              f2Qaa,f2Qii,f2Qai,                        &
                              f1Qa,f1Qi,nb_act1,JJ,                     &
                              tcor2a,tcor2i,tcor1a,trota,               &
                              Basis2n,                                  &
                              wherm,Vinact,ScalOp)
      USE mod_system
      USE mod_nDindex
      use mod_PrimOp, only: PrimOp_t, param_d0matop
      USE mod_basis
      IMPLICIT NONE

!----- variables for the construction of H ---------------------------
      TYPE (PrimOp_t)   :: PrimOp

      integer                  :: nb_Op
      TYPE (param_d0MatOp)     :: d0MatHADAOp(nb_Op)


      integer :: nb_inact2n,nb_act1
!------ f2Q and f1Q in programm order ------------------
      real (kind=Rkind) :: f1Qa(nb_act1),f1Qi(nb_inact2n)
      real (kind=Rkind) :: f2Qaa(nb_act1,nb_act1)
      real (kind=Rkind) :: f2Qii(nb_inact2n,nb_inact2n)
      real (kind=Rkind) :: f2Qai(nb_act1,nb_inact2n)
      real (kind=Rkind) :: tcor2a(nb_act1,3),tcor2i(nb_inact2n,3)
      real (kind=Rkind) :: tcor1a(3),trota(3,3)

      real (kind=Rkind) :: wherm

      real (kind=Rkind) :: d0herm_ij(nb_inact2n)
      real (kind=Rkind) :: d1herm_ij(nb_inact2n)
      real (kind=Rkind) :: d2herm_ij(nb_inact2n)

      real (kind=Rkind) :: gi(nb_inact2n),ga(nb_act1)

      real (kind=Rkind) :: d1xa(nb_inact2n,nb_act1)
      real (kind=Rkind) :: d2xaa(nb_inact2n,nb_act1,nb_act1)

      real (kind=Rkind) :: d0c(nb_inact2n,nb_inact2n)
      real (kind=Rkind) :: d1c(nb_inact2n,nb_inact2n,nb_act1)
      real (kind=Rkind) :: d2c(nb_inact2n,nb_inact2n,nb_act1,nb_act1)

      real (kind=Rkind) :: d0cd0c(nb_inact2n,nb_inact2n,nb_inact2n)
      real (kind=Rkind) :: d1xad1xa(nb_inact2n,nb_act1,nb_act1)
      real (kind=Rkind) :: d1xad0c(nb_inact2n,nb_inact2n,nb_act1)

      real (kind=Rkind) :: d1lnN(nb_act1),d2lnN(nb_act1,nb_act1)

      integer :: ind_quadra(nb_inact2n)
      integer :: ind_basis(nb_inact2n)

!----- for the inactive basis sets ----------------------------------
      TYPE (basis)  :: Basis2n

!     variables pour la diagonalisation sur les fonctions harmoniques
      integer           :: nb_harm


      real (kind=Rkind) :: d0f_harm(Basis2n%nb,1)


      real (kind=Rkind) :: d0fW(1,Basis2n%nb)
      real (kind=Rkind) :: Veffd0fW(1,Basis2n%nb)
      real (kind=Rkind) :: T1d0fW(1,Basis2n%nb,nb_act1)
      real (kind=Rkind) :: T2d0fW(1,Basis2n%nb,nb_act1,nb_act1)
      real (kind=Rkind) :: ScalOpd0fW(1,Basis2n%nb,PrimOp%nb_scalar_Op)
      real (kind=Rkind) :: Tcor2d0fW(1,Basis2n%nb,nb_act1,3)
      real (kind=Rkind) :: Tcor1d0fW(1,Basis2n%nb,3)
      real (kind=Rkind) :: Trotd0fW(1,Basis2n%nb,3,3)

      real (kind=Rkind) :: valS
      real (kind=Rkind) :: Tinact,Vinact,Veff

      real (kind=Rkind) :: valS_ij,TVinact_ij,Veff_ij
      real (kind=Rkind) :: valScalOp_ij(PrimOp%nb_scalar_Op)
      real (kind=Rkind) :: valScalOp(PrimOp%nb_scalar_Op)
      real (kind=Rkind) :: rho,T1(nb_act1),T2(nb_act1,nb_act1)
      real (kind=Rkind) :: Tcor2(nb_act1,3),Tcor1(3),Trot(3,3)
      real (kind=Rkind) :: ScalOp(PrimOp%nb_scalar_Op)

      integer :: JJ

!     - local variables -----------
      integer :: i,j,k,ij,kl,nx
      integer :: iterm,i_term,k_term,iOp,err_sub


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_mat6_HST'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Basis2n%nb or nb_harm',Basis2n%nb
        write(out_unitp,*)
        write(out_unitp,*) 'd0f_harm',d0f_harm
        !---------------------------------------------------
        write(out_unitp,*) 'Matrices of Heff',Basis2n%nb
        write(out_unitp,*)

        DO iOp=1,nb_Op
          write(out_unitp,*) ' iOp:',iOp
          DO k_term=1,d0MatHADAOp(iOp)%nb_term
            write(out_unitp,*) ' deriv_term:',                          &
                             d0MatHADAOp(iOp)%derive_termQact(:,k_term)
            CALL Write_Mat(d0MatHADAOp(iOp)%ReVal(:,:,k_term),out_unitp,5)
          END DO
          IF (d0MatHADAOp(iOp)%cplx) THEN
            write(out_unitp,*) ' cplx Op:'
            CALL Write_Mat(d0MatHADAOp(iOp)%ImVal(:,:),out_unitp,5)
          END IF
        END DO
        !---------------------------------------------------

         write(out_unitp,*) 'nb_inact2n',nb_inact2n
         write(out_unitp,*)
         write(out_unitp,*) 'rho,Vinact',rho,Vinact
         write(out_unitp,*) 'T1,T2',T1,T2
         write(out_unitp,*)
         write(out_unitp,*) 'd1xa',d1xa
         write(out_unitp,*) 'd2xaa',d2xaa
         write(out_unitp,*) 'd0c',d0c
         write(out_unitp,*) 'd1c',d1c
         write(out_unitp,*) 'd2c',d2c
         write(out_unitp,*) 'd1lnN,d2lnN',d1lnN,d2lnN
         write(out_unitp,*) 'f2Qaa,f2Qii,f2Qai',f2Qaa,f2Qii,f2Qai
         write(out_unitp,*) 'tcor2a,tcor2i',tcor2a,tcor2i
         write(out_unitp,*) 'tcor1a,trota,',tcor1a,trota

         write(out_unitp,*)
         write(out_unitp,*) 'ind_quadra',ind_quadra
         write(out_unitp,*) 'wherm,Vinact',wherm,Vinact

       END IF
!-----------------------------------------------------------

!     -------------------------------------------------
!     - d0c * d0c
!     done in sub_HST_harm
!     -------------------------------------------------

!     -------------------------------------------------
!     - d1xa * d1xa
!     -------------------------------------------------
      DO i=1,nb_act1
      DO j=i,nb_act1
      DO k=1,nb_inact2n
        d1xad1xa(k,j,i) = d1xa(k,i) * d1xa(k,j)
      END DO
      END DO
      END DO
!     -------------------------------------------------
!     -------------------------------------------------

!     -------------------------------------------------
!     - d1xa * d0c
!     -------------------------------------------------
      DO i=1,nb_act1
      DO j=1,nb_inact2n
      DO k=1,nb_inact2n
        d1xad0c(k,j,i) = d1xa(k,i) * d0c(j,k)
      END DO
      END DO
      END DO
!     -------------------------------------------------
!     -------------------------------------------------

!     write(out_unitp,*) 'ind_quadra',ind_quadra
!     --------------------------------------------------------
!     construction des matrices
      DO ij=1,Basis2n%nb
        CALL calc_nDindex(Basis2n%nDindB,ij,ind_basis(:),err_sub)

        IF (err_sub /= 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) '  from Basis2n%nDindB'
          STOP 'calc_nDindex'
        END IF

        DO i=1,Basis2n%nb_basis
          !CALL Rec_d0d1d2bnD(d0herm_ij(i),                              &
          !                       reshape((/ d1herm_ij(i) /),(/1,1/)),   &
          !                       reshape((/ d2herm_ij(i) /),(/1,1,1/)), &
          !      Basis2n%tab_Pbasis(i)%Pbasis,ind_quadra(i),ind_basis(i))
          d0herm_ij(i) = Basis2n%tab_Pbasis(i)%Pbasis%dnRGB%d0(          &
                                             ind_quadra(i),ind_basis(i))
          d1herm_ij(i) = Basis2n%tab_Pbasis(i)%Pbasis%dnRGB%d1(          &
                                             ind_quadra(i),ind_basis(i),1)
          d2herm_ij(i) = Basis2n%tab_Pbasis(i)%Pbasis%dnRGB%d2(          &
                                           ind_quadra(i),ind_basis(i),1,1)
        END DO

!       write(out_unitp,*) ' tab_Pbasis(i)%nq is odd'
!       here I should multiply by rho, but since rho MUST not be a function
!       of the inactive variables, I'll do that after
        CALL calc_Tinact_new(Tinact,Veff,T1,T2,                         &
                             d0c,d1c,d1xa,d2xaa,                        &
                             d0cd0c,d1xad0c,d1xad1xa,                   &
                             d1lnN,d2lnN,                               &
                             f2Qaa,f2Qii,f2Qai,                         &
                             f1Qa,f1Qi,                                 &
                             d0herm_ij,d1herm_ij,d2herm_ij,             &
                             nb_inact2n,nb_act1)

        d0fW(1,ij)       = d0f_harm(ij,1)      * wherm
        Veffd0fW(1,ij)   = Veff*wherm   + Vinact*d0fW(1,ij)
        IF (nb_act1 > 0) THEN
          T1d0fW(1,ij,:)   = T1(:)               * wherm
          T2d0fW(1,ij,:,:) = T2(:,:)             * wherm
        END IF


        IF (JJ > 0) THEN
          CALL calc_Tcorrot(Tcor2,Tcor1,Trot,                           &
                            d0c,d1c,d1xa,d2xaa,                         &
                            d0cd0c,d1xad0c,d1xad1xa,                    &
                            d1lnN,d2lnN,                                &
                            tcor2a,tcor2i,tcor1a,trota,                 &
                            d0herm_ij,d1herm_ij,d2herm_ij,              &
                            gi,ga,                                      &
                            nb_inact2n,nb_act1)


          Tcor1d0fW(1,ij,:)   = Tcor1(:)   * d0fW(1,ij)
          Tcor2d0fW(1,ij,:,:) = Tcor2(:,:) * d0fW(1,ij)
          Trotd0fW(1,ij,:,:)  = Trot(:,:)  * d0fW(1,ij)
        END IF

        IF (PrimOp%calc_scalar_Op) THEN
          ScalOpd0fW(1,ij,:) = ScalOp(:) * d0fW(1,ij)
        END IF

      END DO


      !---------------------------------------------------------
      !-- matrices construction --------------------------------
      iOp=1 ! H
      !-- Veff -------------------------------------------------
      iterm = d0MatHADAOp(iOp)%derive_term_TO_iterm(0,0)
      d0MatHADAOp(iOp)%ReVal(:,:,iterm) =                               &
                                   d0MatHADAOp(iOp)%ReVal(:,:,iterm) +  &
                        matmul(d0f_harm,Veffd0fW(:,:))
      !-- T1 ---------------------------------------------------
      DO i=1,nb_act1
        iterm = d0MatHADAOp(iOp)%derive_term_TO_iterm(i,0)
        d0MatHADAOp(iOp)%ReVal(:,:,iterm) =                             &
                                   d0MatHADAOp(iOp)%ReVal(:,:,iterm) +  &
                        matmul(d0f_harm,T1d0fW(:,:,i))
      END DO
      !-- T2 ---------------------------------------------------
      DO i=1,nb_act1
      DO j=i,nb_act1
        iterm = d0MatHADAOp(iOp)%derive_term_TO_iterm(i,j)
        d0MatHADAOp(iOp)%ReVal(:,:,iterm) =                             &
                                   d0MatHADAOp(iOp)%ReVal(:,:,iterm) +  &
                        matmul(d0f_harm(:,:),T2d0fW(:,:,i,j))
      END DO
      END DO

      IF (JJ > 0) THEN
        DO i=-3,-1
        DO j=i,-1
          iterm = d0MatHADAOp(iOp)%derive_term_TO_iterm(i,j)
          d0MatHADAOp(iOp)%ReVal(:,:,iterm) =                           &
                                   d0MatHADAOp(iOp)%ReVal(:,:,iterm) +  &
                        matmul(d0f_harm,Trotd0fW(:,:,-j,-i))
        END DO
        END DO

        DO i=-3,-1
          iterm = d0MatHADAOp(iOp)%derive_term_TO_iterm(i,0)
          d0MatHADAOp(iOp)%ReVal(:,:,iterm) =                           &
                                   d0MatHADAOp(iOp)%ReVal(:,:,iterm) +  &
                        matmul(d0f_harm,Tcor1d0fW(:,:,-i))
          DO j=1,nb_act1
            iterm = d0MatHADAOp(iOp)%derive_term_TO_iterm(i,j)
            d0MatHADAOp(iOp)%ReVal(:,:,iterm) =                         &
                                   d0MatHADAOp(iOp)%ReVal(:,:,iterm) +  &
                                     matmul(d0f_harm,Tcor2d0fW(:,:,j,-i))
          END DO
        END DO
      END IF

      iOp = 2 ! S
      iterm = d0MatHADAOp(iOp)%derive_term_TO_iterm(0,0)
      d0MatHADAOp(iOp)%ReVal(:,:,iterm) =                               &
                                   d0MatHADAOp(iOp)%ReVal(:,:,iterm) +  &
                        matmul(d0f_harm,d0fW(:,:))

      DO iOp=3,nb_Op
        iterm = d0MatHADAOp(iOp)%derive_term_TO_iterm(0,0)
        d0MatHADAOp(iOp)%ReVal(:,:,iterm) =                             &
                                   d0MatHADAOp(iOp)%ReVal(:,:,iterm) +  &
                        matmul(d0f_harm,ScalOpd0fW(:,:,iOp-2))
      END DO
      !---------------------------------------------------------
      !---------------------------------------------------------


      !---------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Matrices of Heff',Basis2n%nb
        DO iOp=1,nb_Op
          write(out_unitp,*) ' iOp:',iOp
          DO k_term=1,d0MatHADAOp(iOp)%nb_term
            write(out_unitp,*) ' deriv_term:',                          &
                             d0MatHADAOp(iOp)%derive_termQact(:,k_term)
            CALL Write_Mat(d0MatHADAOp(iOp)%ReVal(:,:,k_term),out_unitp,5)
          END DO
          IF (d0MatHADAOp(iOp)%cplx) THEN
            write(out_unitp,*) ' cplx Op:'
            CALL Write_Mat(d0MatHADAOp(iOp)%ImVal(:,:),out_unitp,5)
          END IF
        END DO
        write(out_unitp,*) 'END ',name_sub
      END IF
      !---------------------------------------------------------

      END SUBROUTINE sub_mat6_HST

