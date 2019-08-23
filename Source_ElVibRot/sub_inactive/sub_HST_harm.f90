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
      SUBROUTINE sub_HSOp_inact(iq,freq_only,para_AllOp,                &
                                max_Sii,max_Sij,test,OldPara)

      USE mod_system
      USE mod_dnSVM
      USE mod_Coord_KEO, only : zmatrix, Tnum, get_Qact, qact_to_qdyn_from_activetransfo
      USE mod_basis
      USE mod_Op
      USE mod_PrimOp
      IMPLICIT NONE

      integer       :: iq ! grid point number
      logical, intent(in) :: freq_only
      TYPE (OldParam), intent(inout) :: OldPara


!----- variables for the construction of H ---------------------------
      TYPE (param_AllOp) :: para_AllOp
      real (kind=Rkind) :: max_Sii,max_Sij

!---- variable pour le test -----------------------------------------
      logical :: test



!----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix), pointer     :: mole      ! true pointer
      TYPE (Tnum),    pointer     :: para_Tnum ! true pointer

!----- working variables ----------------------------------------
      logical                           :: pot,keo,KEO_bis
      real (kind=Rkind)                 :: WrhonD,rho
      real (kind=Rkind)                 :: max1_Sii,max1_Sij

      real (kind=Rkind) , allocatable   :: d0ehess(:),d0Qeq(:)

!----- variables to store the operators of a single grid point ----------------
      TYPE (param_d0MatOp), allocatable :: d0MatOp(:)

      integer                           :: iOp,k_term
      logical                           :: JacSave,sqRhoOVERJacSave

      real (kind=Rkind), allocatable    :: Qact(:),Qdyn(:)
      real (kind=Rkind), save           :: Q1

      integer                           :: i_SG,iq_SG,err_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub = 'sub_HSOp_inact'
      integer :: err_mem,memory
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      mole       => para_AllOp%tab_Op(1)%mole
      para_Tnum  => para_AllOp%tab_Op(1)%para_Tnum

      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_Op',para_AllOp%nb_Op,shape(para_AllOp%tab_Op)

        !-------------------------------------------------------
        write(out_unitp,*) 'nb_bie',para_AllOp%tab_Op(1)%nb_bie
        write(out_unitp,*) 'nrho',para_AllOp%tab_Op(1)%para_Tnum%nrho
        write(out_unitp,*) 'JJ',para_AllOp%tab_Op(1)%para_Tnum%JJ

        !-------------------------------------------------------
        write(out_unitp,*) ' Max Overlap',max_Sii,max_Sij

        DO iOp=1,para_AllOp%nb_Op
            write(out_unitp,*) ' iOp, n_Op, name_Op:',iOp,              &
                               para_AllOp%tab_Op(iOp)%n_Op,             &
                               para_AllOp%tab_Op(iOp)%name_Op
            DO k_term=1,para_AllOp%tab_Op(iOp)%nb_term
              write(out_unitp,*) ' name_Op, deriv_term:',               &
                             para_AllOp%tab_Op(iOp)%name_Op,            &
                 para_AllOp%tab_Op(iOp)%derive_termQact(:,k_term)
            END DO
            IF (para_AllOp%tab_Op(iOp)%cplx) THEN
              write(out_unitp,*) ' cplx name_Op:',para_AllOp%tab_Op(iOp)%name_Op
            END IF
        END DO
        !-----------------------------------------------------

      END IF
      !-----------------------------------------------------------

      !-----------------------------------------------------------
      !----- New nD-Grid points ----------------------------------
      !-----------------------------------------------------------
      CALL alloc_NParray(Qact,(/mole%nb_var/),'Qact',name_sub)

      CALL get_Qact(Qact,mole%ActiveTransfo) ! rigid, flexible coordinates

      IF (iq > 0 .OR. .NOT. test) THEN
        !-----------------------------------------------------
        !calculation of Qact
        CALL Rec_Qact(Qact,                                             &
                      para_AllOp%tab_Op(1)%para_AllBasis%BasisnD,iq,    &
                      mole,OldPara)
      END IF
      !-----------------------------------------------------------
      !-----------------------------------------------------------

      !-----------------------------------------------------------
      !------ special case if nb_inact2n=0 -----------------------
      !-----------------------------------------------------------
      IF (mole%nb_inact2n == 0) THEN
        pot = .TRUE.
        KEO = .TRUE.
        iOp=1 ! H
        IF (para_AllOp%tab_Op(iOp)%pot_only) KEO = .FALSE.
        IF (para_AllOp%tab_Op(iOp)%T_only)   pot = .FALSE.

        IF (test) THEN
          para_Tnum%WriteT = .TRUE.
          mole%WriteCC     = .TRUE.
        END IF

        allocate(d0MatOp(para_AllOp%tab_Op(1)%para_PES%nb_scalar_Op+2))
        DO iOp=1,size(d0MatOp)
          CALL Init_d0MatOp(d0MatOp(iOp),para_AllOp%tab_Op(iOp)%param_TypeOp,&
                            para_AllOp%tab_Op(iOp)%para_PES%nb_elec)
        END DO

        CALL get_d0MatOp_AT_Qact(Qact,d0MatOp,mole,                     &
                                 para_Tnum,para_AllOp%tab_Op(1)%para_PES)

        IF (.NOT. pot) THEN ! remove the potential part
          d0MatOp(1)%ReVal(:,:,1) = ZERO
        END IF

        JacSave = (para_Tnum%nrho == 0) .AND.                           &
                             (para_AllOp%tab_Op(1)%type_Op == 10) .AND. &
                                   .NOT. para_AllOp%tab_Op(1)%direct_KEO
        JacSave = JacSave .OR. (para_Tnum%nrho == 0) .AND.              &
                             (para_AllOp%tab_Op(1)%type_Op == 10) .AND. &
                                  para_AllOp%tab_Op(1)%direct_KEO .AND. &
                       para_AllOp%tab_Op(1)%BasisnD%SparseGrid_type /= 4
!JacSave = .TRUE.
        sqRhoOVERJacSave = JacSave .OR. (para_Tnum%nrho == 0) .AND.     &
                                     (para_AllOp%tab_Op(1)%type_Op == 1)

        KEO_bis = (JacSave .OR. sqRhoOVERJacSave) .OR.                  &
             (d0MatOp(1)%type_Op /= 10) .OR. .NOT. d0MatOp(1)%direct_KEO

        IF (.NOT. para_AllOp%tab_Op(1)%para_Tnum%Gcte .AND. KEO .AND. KEO_bis) THEN
          CALL TnumKEO_TO_tab_d0H(Qact,d0MatOp(1),mole,para_Tnum)
        END IF

        IF (sqRhoOVERJacSave) THEN
            !$OMP  CRITICAL (sub_HSOp_inact1_CRIT)
            IF (allocated(para_AllOp%tab_Op(1)%ComOp%sqRhoOVERJac)) THEN
              IF (size(para_AllOp%tab_Op(1)%ComOp%sqRhoOVERJac) /= para_AllOp%tab_Op(1)%nb_qa) THEN
                CALL dealloc_NParray(para_AllOp%tab_Op(1)%ComOp%sqRhoOVERJac, &
                                    'para_AllOp%tab_Op(1)%ComOp%sqRhoOVERJac',name_sub)
              END IF
            END IF

            IF (.NOT. allocated(para_AllOp%tab_Op(1)%ComOp%sqRhoOVERJac)) THEN
                CALL alloc_NParray(para_AllOp%tab_Op(1)%ComOp%sqRhoOVERJac, &
                                        (/ para_AllOp%tab_Op(1)%nb_qa /),   &
                                  'para_AllOp%tab_Op(1)%ComOp%sqRhoOVERJac',name_sub)
            END IF
            !$OMP  END CRITICAL (sub_HSOp_inact1_CRIT)

             para_AllOp%tab_Op(1)%ComOp%sqRhoOVERJac(iq) =              &
                                     sqrt(d0MatOp(1)%rho/d0MatOp(1)%Jac)
        END IF

        IF (JacSave) THEN
            !$OMP  CRITICAL (sub_HSOp_inact2_CRIT)
            IF (allocated(para_AllOp%tab_Op(1)%ComOp%Jac)) THEN
              IF (size(para_AllOp%tab_Op(1)%ComOp%Jac) /= para_AllOp%tab_Op(1)%nb_qa) THEN
                CALL dealloc_NParray(para_AllOp%tab_Op(1)%ComOp%Jac, &
                                    'para_AllOp%tab_Op(1)%ComOp%Jac',name_sub)
              END IF
            END IF

            IF (.NOT. allocated(para_AllOp%tab_Op(1)%ComOp%Jac)) THEN
                CALL alloc_NParray(para_AllOp%tab_Op(1)%ComOp%Jac, &
                                        (/ para_AllOp%tab_Op(1)%nb_qa /),   &
                                  'para_AllOp%tab_Op(1)%ComOp%Jac',name_sub)
            END IF
            !$OMP  END CRITICAL (sub_HSOp_inact2_CRIT)

            para_AllOp%tab_Op(1)%ComOp%Jac(iq) = d0MatOp(1)%Jac
        END IF

        IF (iq > 0 .OR. .NOT. test) THEN
          !calculation of WrhonD
          WrhonD = Rec_WrhonD(para_AllOp%tab_Op(1)%para_AllBasis%BasisnD,iq,OldPara)
        END IF

      !-----------------------------------------------------------
      !------ end special case if nb_inact2n=0 -------------------
      !-----------------------------------------------------------
      ELSE
      !-----------------------------------------------------------
      !------ case nb_inact2n>0 ----------------------------------
      !-----------------------------------------------------------

        IF (iq > 0 .OR. .NOT. test) THEN
          !calculation of WrhonD
          WrhonD = Rec_WrhonD(para_AllOp%tab_Op(1)%para_AllBasis%BasisnD,iq,OldPara)
        END IF

        IF (para_AllOp%tab_Op(1)%para_PES%Type_HamilOp /= 1) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '    Type_HamilOp',para_AllOp%tab_Op(1)%para_PES%Type_HamilOp
          write(out_unitp,*) '    Type_HamilOp MUST be equal to 1 for HADA or cHAC'
          write(out_unitp,*) '    CHECK your data!!'
          STOP
        END IF

        allocate(d0MatOp(para_AllOp%tab_Op(1)%para_PES%nb_scalar_Op+2))

        DO iOp=1,size(d0MatOp)
          CALL Init_d0MatOp(d0MatOp(iOp),para_AllOp%tab_Op(iOp)%param_TypeOp,&
                        para_AllOp%tab_Op(1)%nb_bi*para_AllOp%tab_Op(1)%para_PES%nb_elec)
        END DO


        CALL alloc_NParray(d0ehess,(/ mole%nb_inact2n /),   "d0ehess",name_sub)
        CALL alloc_NParray(d0Qeq,  (/ mole%nb_inact2n /),   "d0Qeq",  name_sub)

        CALL sub_HST7_bhe(Qact,d0Qeq,d0ehess,                           &
                          d0MatOp,para_AllOp%nb_Op,rho,para_AllOp,      &
                          para_AllOp%tab_Op(1)%para_AllBasis%Basis2n,   &
                          freq_only,test)

        !------------------------------------------------------------
        !------------------------------------------------------------
        IF (.NOT. freq_only) THEN
          !-------------------------------------------------------
          ! Matrix analysis of S_bhe
          ! S
          iOp = 2
          CALL sub_ana_S(d0MatOp(iOp)%ReVal(:,:,1),para_AllOp%tab_Op(1)%nb_bie,&
                         max1_Sii,max1_Sij,.FALSE.)
          IF (max1_Sii > max_Sii) max_Sii = max1_Sii
          IF (max1_Sij > max_Sij) max_Sij = max1_Sij
          !-------------------------------------------------------
          !-------------------------------------------------------
        END IF
        !------------------------------------------------------------
        !------------------------------------------------------------

        CALL dealloc_NParray(d0ehess,"d0ehess",name_sub)
        CALL dealloc_NParray(d0Qeq,  "d0Qeq",  name_sub)

      !-----------------------------------------------------------
      !------ end case nb_inact2n>0 ------------------------------
      !-----------------------------------------------------------
      END IF
      !-----------------------------------------------------------
      !-----------------------------------------------------------
      !-----------------------------------------------------------



      !-----------------------------------------------------------
      !-----------------------------------------------------------
      IF (print_level > 0) THEN

        IF (para_AllOp%tab_Op(1)%nb_qa <= max_nb_G_FOR_print) THEN
          IF (freq_only) write(out_unitp,*) 'freq_only,iq',iq

          IF (JacSave) THEN
            write(out_unitp,111) iq,Qact(1:mole%nb_act1),d0MatOp(1)%ReVal(1,1,:),d0MatOp(1)%Jac,d0MatOp(1)%rho
          ELSE
            write(out_unitp,111) iq,Qact(1:mole%nb_act1),d0MatOp(1)%ReVal(1,1,:)
          END IF
 111      format('Grid: ',i6,50(1x,f18.10))

        ELSE
          IF (mod(iq,max(1,int(para_AllOp%tab_Op(1)%nb_qa/10))) == 0)   &
            write(out_unitp,'(a)',ADVANCE='no') '---'
        END IF
        CALL flush_perso(out_unitp)
      END IF
      !-----------------------------------------------------------
      !-----------------------------------------------------------

      !-----------------------------------------------------------
      !-----------------------------------------------------------
      IF ( .NOT. freq_only .AND. iq > 0 .AND. .NOT. test) THEN
        CALL alloc_NParray(Qdyn,(/mole%nb_var/),'Qdyn',name_sub)

        CALL Qact_TO_Qdyn_FROM_ActiveTransfo(Qact,Qdyn,mole%ActiveTransfo)

        CALL sub_Save_GridFile_AllOp(iq,d0MatOp,para_AllOp%nb_Op,       &
                                     Qdyn,mole%nb_var,Qact,mole%nb_act1,&
                                     para_AllOp,WrhonD)
        CALL sub_Save_GridMem_AllOp(iq,d0MatOp,para_AllOp%nb_Op,para_AllOp)

        CALL dealloc_NParray(Qdyn,'Qdyn',name_sub)
      END IF
      !-----------------------------------------------------------
      !-----------------------------------------------------------

      !-----------------------------------------------------------
      !-----------------------------------------------------------
      IF (test .OR. debug) THEN

        write(out_unitp,*) 'JJ',para_Tnum%JJ
        write(out_unitp,*) 'Qact',Qact

        IF (mole%nb_inact2n > 0) THEN
          write(out_unitp,*) ' Max Overlap',max1_Sii,max1_Sij
        END IF

        DO iOp=1,para_AllOp%nb_Op
          write(out_unitp,*) ' iOp, name_Op:',iOp,                      &
                               para_AllOp%tab_Op(iOp)%name_Op
          DO k_term=1,d0MatOp(iOp)%nb_term
            IF (para_Tnum%JJ == 0 .AND.                                 &
                count(d0MatOp(iOp)%derive_termQact(:,k_term) <0) > 0) CYCLE
            write(out_unitp,*) ' name_Op, deriv_term:',                 &
                             para_AllOp%tab_Op(iOp)%name_Op,            &
                             d0MatOp(iOp)%derive_termQact(:,k_term)
            CALL Write_Mat(d0MatOp(iOp)%ReVal(:,:,k_term),out_unitp,5)
          END DO
          IF (d0MatOp(iOp)%cplx) THEN
            write(out_unitp,*) ' cplx name_Op:',para_AllOp%tab_Op(iOp)%name_Op
            CALL Write_Mat(d0MatOp(iOp)%Imval(:,:),out_unitp,5)
          END IF
        END DO

      END IF

      !-----------------------------------------------------------
      ! deallocations ....
      CALL dealloc_NParray(Qact,'Qact',name_sub)

      DO iOp=1,size(d0MatOp)
        CALL dealloc_d0MatOp(d0MatOp(iOp))
      END DO
      deallocate(d0MatOp)

      !-----------------------------------------------------------
      !-----------------------------------------------------------


      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF

      END SUBROUTINE sub_HSOp_inact

!=============================================================
!
!     Hamiltonian harmonic nD matrix calculation:
!                                      H_bhe(nb_bie,nb_bie)
!     Overlap harmonic matrix calculation:
!                                      S_bhe(nb_bie,nb_bie)
!     Kinetic matrix calculation (non adiabatic coupling):
!                             T1_bhe(nb_bie,nb_bie,nb_act1)
!                     T2_bhe(nb_bie,nb_bie,nb_act1,nb_act1)
!     Effectivie potential matrix (non adiabatic coupling):
!                                   Veff_bhe(nb_bie,nb_bie)
!
!=============================================================
     SUBROUTINE sub_HST7_bhe(Qact,d0Qeq,d0ehess,d0MatHADAOp,nb_Op,    &
                             rho,para_AllOp,Basis2n,freq_only,test)

      USE mod_system
      USE mod_dnSVM
      USE mod_nDindex
      USE mod_Coord_KEO, only : zmatrix, Tnum, Qinact2n_TO_Qact_FROM_ActiveTransfo
      USE mod_basis
      USE mod_Op
      USE mod_PrimOp
      IMPLICIT NONE

!----- variables for the construction of H ---------------------------
      TYPE (param_AllOp)       :: para_AllOp

      integer                  :: nb_Op
      TYPE (param_d0MatOp)     :: d0MatHADAOp(nb_Op)
      TYPE (param_d0MatOp)     :: d0MatOp(nb_Op)

!----- variables for the active and inactive namelists ----------------
      integer :: ind_quadra(para_AllOp%tab_Op(1)%mole%nb_inact2n)
      integer :: ind_basis(para_AllOp%tab_Op(1)%mole%nb_inact2n)

!----- for the inactive basis-sets ------------------------------------
      TYPE (basis)            :: Basis2n
      TYPE (P_basis), pointer :: tab_Pbasis2n(:)

!----- Coordinates ------------------------------------
      real (kind=Rkind)       :: Qact(para_AllOp%tab_Op(1)%mole%nb_var)

!----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix),pointer     :: mole      ! true pointer
      TYPE (Tnum),pointer        :: para_Tnum ! true pointer

!------ harmonic matrix : H_bhe, S_bhe  ---------------------
      integer  :: nb_bie

!------ pour les frequences -------------------------------
      real (kind=Rkind) :: d0norme,                                     &
                           d0Qeq(para_AllOp%tab_Op(1)%mole%nb_inact2n), &
                           d0ehess(para_AllOp%tab_Op(1)%mole%nb_inact2n)

!------ for the potential ---------------------------------
!       Vinact   : harmonic part (always real)
!       V0       : active part (real)
!       HarD     : if .TRUE. => Harmonic Domain (path)
!       pot_cplx : if .TRUE. => complex potential (for nb_act1 variables)
      real (kind=Rkind) ::  Vinact,V0


!---- variable pour le test -----------------------------------------
      logical           :: test,freq_only
      real (kind=Rkind) :: d0psi
      integer           :: ih


!----- working variables ----------------------------------------
      integer :: iOp
      integer :: i_point,i_modif,i,j,k,l,ib,iq
      real (kind=Rkind) :: ScalOp(para_AllOp%tab_Op(1)%para_PES%nb_scalar_Op)
      real (kind=Rkind) :: vep,rho,wnDh
      real (kind=Rkind) :: pot0_corgrad
      logical :: pot,KEO

!------ for the frequencies -------------------------------
        integer :: nderiv
        TYPE (Type_dnMat)     :: dnC,dnC_inv      ! derivative with respect to Qact1
        TYPE (Type_dnVec)     :: dnQeq            ! derivative with respect to Qact1
        TYPE (Type_dnVec)     :: dnEHess          ! derivative with respect to Qact1
        TYPE (Type_dnVec)     :: dnGrad           ! derivative with respect to Qact1
        TYPE (Type_dnMat)     :: dnHess           ! derivative with respect to Qact1
        TYPE (Type_dnS)       :: dnLnN            ! derivative with respect to Qact1

      real (kind=Rkind) , allocatable ::                                &
       f1Q(:),f2Q(:,:),Tcor2(:,:),                                      &
       f1Qa(:),f1Qi(:),f2Qaa(:,:),f2Qii(:,:),f2Qai(:,:),                &
       tcor2a(:,:),tcor2i(:,:),tcor1a(:),trota(:,:)


      real (kind=Rkind) , allocatable ::                                &
       d0x(:),d1x(:,:),d2x(:,:,:),                                      &
       d0cd0c(:,:,:),                                                   &
       Qinact(:),deltaQ(:),deltaQ2(:)

      real (kind=Rkind) , allocatable :: d0f_bhe(:)


      integer :: nb_act1,nb_act,nb_inact2n
      integer :: i_term,k_term
      integer :: ihe
      integer :: err_sub
      integer :: nq

      integer, parameter :: nq_write_HADA = 10000

!----- FUNCTION --------------------------------------------------
      real (kind=Rkind) ::      pot_rest
!---- FUNCTION ---------------------------------------------------
!-----------------------------------------------------------------
!-----------------------------------------------------------------


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_HST7_bhe'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------

      mole         => para_AllOp%tab_Op(1)%mole
      para_Tnum    => para_AllOp%tab_Op(1)%para_Tnum
      tab_Pbasis2n => Basis2n%tab_Pbasis

      nb_act      = mole%nb_act
      nb_act1     = mole%nb_act1
      nb_inact2n  = mole%nb_inact2n
      nb_bie      = para_AllOp%tab_Op(1)%nb_bie



!-----------------------------------------------------------------
!-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_Op',para_AllOp%nb_Op,shape(para_AllOp%tab_Op)

        !-------------------------------------------------------
        write(out_unitp,*) 'Qact',Qact
        write(out_unitp,*) 'nb_bie',nb_bie
        write(out_unitp,*) 'nrho',para_Tnum%nrho
        write(out_unitp,*) 'JJ',para_Tnum%JJ
        !-------------------------------------------------------
        DO iOp=1,nb_Op
          write(out_unitp,*) ' iOp, name_Op:',iOp,                      &
                               para_AllOp%tab_Op(iOp)%name_Op
          DO k_term=1,d0MatHADAOp(iOp)%nb_term
            IF (para_AllOp%tab_Op(1)%para_Tnum%JJ == 0 .AND.            &
                count(d0MatHADAOp(iOp)%derive_termQact(:,k_term) <0) > 0) CYCLE
            write(out_unitp,*) ' name_Op, deriv_term:',                 &
                             para_AllOp%tab_Op(iOp)%name_Op,            &
                             d0MatHADAOp(iOp)%derive_termQact(:,k_term)
            CALL Write_Mat(d0MatHADAOp(iOp)%ReVal(:,:,k_term),out_unitp,5)
          END DO
          IF (d0MatHADAOp(iOp)%cplx) THEN
            write(out_unitp,*) ' cplx name_Op:',para_AllOp%tab_Op(iOp)%name_Op
            CALL Write_Mat(d0MatHADAOp(iOp)%ImVal(:,:),out_unitp,5)
          END IF
        END DO
        CALL flush_perso(out_unitp)
        !-------------------------------------------------------

!       -----------------------------------------------------
        write(out_unitp,*)
        CALL Write_mole(mole)
        write(out_unitp,*)
!       -----------------------------------------------------

        CALL RecWrite_basis(Basis2n)
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

!-----------------------------------------------------------------
!-----------------------------------------------------------------
!     -- Matrix initialisation -----------------------------------
      !write(6,*) 'nb_Op',nb_Op
      DO iOp=1,nb_Op
        d0MatHADAOp(iOp)%ReVal(:,:,:) = ZERO
        IF (d0MatHADAOp(iOp)%cplx) d0MatHADAOp(iOp)%ImVal(:,:) = ZERO
      END DO


!     ------ this memory is free at the subroutine end -----------
      DO iOp=1,size(d0MatOp)
        CALL Init_d0MatOp(d0MatOp(iOp),d0MatHADAOp(iOp)%param_TypeOp,   &
                                para_AllOp%tab_Op(iOp)%para_PES%nb_elec)
      END DO

      CALL alloc_NParray(f1Q,    (/nb_act/),                    "f1Q",   name_sub)
      CALL alloc_NParray(f2Q,    (/nb_act, nb_act/),            "f2Q",   name_sub)
      CALL alloc_NParray(Tcor2,  (/nb_act,3/),                  "Tcor2", name_sub)
      CALL alloc_NParray(f1Qa,   (/nb_act1/),                   "f1Qa",  name_sub)
      CALL alloc_NParray(f1Qi,   (/nb_inact2n/),                "f1Qi",  name_sub)
      CALL alloc_NParray(f2Qaa,  (/nb_act1, nb_act1/),          "f2Qaa", name_sub)
      CALL alloc_NParray(f2Qii,  (/nb_inact2n, nb_inact2n/),    "f2Qii", name_sub)
      CALL alloc_NParray(f2Qai,  (/nb_act1, nb_inact2n/),       "f2Qai", name_sub)
      CALL alloc_NParray(tcor2a, (/nb_act1,3/),                 "tcor2a",name_sub)
      CALL alloc_NParray(tcor2i, (/nb_inact2n,3/),              "tcor2i",name_sub)
      CALL alloc_NParray(tcor1a, (/3/),                         "tcor1a",name_sub)
      CALL alloc_NParray(trota,  (/3,3/),                       "trota", name_sub)

      CALL alloc_NParray(d0x,    (/nb_inact2n/),                        "d0x",    name_sub)
      CALL alloc_NParray(d1x,    (/nb_inact2n, nb_act1/),               "d1x",    name_sub)
      CALL alloc_NParray(d2x,    (/nb_inact2n, nb_act1, nb_act1/),      "d2x",    name_sub)
      CALL alloc_NParray(d0cd0c, (/nb_inact2n, nb_inact2n, nb_inact2n/),"d0cd0c", name_sub)
      CALL alloc_NParray(Qinact, (/nb_inact2n/),                        "Qinact", name_sub)
      CALL alloc_NParray(deltaQ2,(/nb_inact2n/),                        "deltaQ2",name_sub)
      CALL alloc_NParray(deltaQ, (/nb_inact2n/),                        "deltaQ", name_sub)

!-----------------------------------------------------------------
!-----------------------------------------------------------------
!     --- frequencies and normal modes calculation
!     at Qact(nb_var)
!     equilibrium inactive variables and hessian matrix

      nderiv = 2
      CALL alloc_dnSVM(dnC,    nb_inact2n,nb_inact2n,nb_act1,nderiv)
      CALL alloc_dnSVM(dnC_inv,nb_inact2n,nb_inact2n,nb_act1,nderiv)
      CALL alloc_dnSVM(dnQeq,  nb_inact2n,           nb_act1,nderiv)
      CALL alloc_dnSVM(dnEHess,nb_inact2n,           nb_act1,nderiv)
      CALL alloc_dnSVM(dnHess, nb_inact2n,nb_inact2n,nb_act1,nderiv)
      CALL alloc_dnSVM(dnGrad, nb_inact2n,           nb_act1,nderiv)
      CALL alloc_dnSVM(dnLnN,                        nb_act1,nderiv)


      CALL sub_dnfreq_4p(dnQeq,dnC,dnLnN,dnEHess,dnHess,dnGrad,dnC_inv, &
                         pot0_corgrad,    &
                         Qact,para_Tnum,mole,mole%RPHTransfo_inact2n,   &
                         nderiv,test)

      d0Qeq      = dnQeq%d0
      d0ehess    = dnEHess%d0

      IF (print_level > 0) write(out_unitp,12) Qact(1:nb_act1),dnQeq%d0
 12   format('Qeq',20(' ',f10.6))
      IF (print_level > 0) write(out_unitp,11) Qact(1:nb_act1),pot0_corgrad
 11   format('pot0_corgrad',10(' ',f10.6))
      CALL flush_perso(out_unitp)

!-----------------------------------------------------------------
!-----------------------------------------------------------------

      IF (.NOT. freq_only) THEN

!------------------------------------------------------------
!------------------------------------------------------------
!     --- Gauss-Hermite quadrature loop ---------------------
!     -------------------------------------------------------
!     -------------------------------------------------------
!     -------------------------------------------------------

      pot = .TRUE.
      KEO = .TRUE.
      IF (para_AllOp%tab_Op(1)%pot_only) KEO = .FALSE.
      IF (para_AllOp%tab_Op(1)%T_only)   pot = .FALSE.


      CALL alloc_NParray(d0f_bhe,(/nb_bie/),"d0f_bhe",name_sub)

      CALL calc_d0cd0c(d0cd0c,dnC%d0,nb_inact2n)


      nq = get_nq_FROM_basis(Basis2n)
      IF (print_level > 0 .AND. nq > nq_write_HADA) &
                write(out_unitp,'(a)',ADVANCE='no') 'Grid_HADA (%): 0'
      CALL flush_perso(out_unitp)


      DO i_point=1,nq
        IF (print_level > 0 .AND. nq > nq_write_HADA .AND.              &
                                  mod(i_point,max(1,(nq/10))) == 0) THEN

          write(out_unitp,'(a,i3)',ADVANCE='no') ' -',                  &
              int(real(i_point,kind=Rkind)*HUNDRED/real(nq,kind=Rkind))
          CALL flush_perso(out_unitp)
        END IF


!       -----------------------------------------------------
!       IF the weight of the sparse grid is zero we skip the iteration
!       the weight in nD
        IF (Basis2n%SparseGrid_type == 3) THEN
          IF (Basis2n%w(i_point) == ZERO) CYCLE
          wnDh = Basis2n%w(i_point)
        ELSE
          wnDh = Rec_WrhonD(Basis2n,i_point)
        END IF

!       -----------------------------------------------------
!       d0f_bhe in nD
        DO ib=1,Basis2n%nb
          d0f_bhe(ib) = Rec_d0bnD(Basis2n,i_point,ib)
        END DO
!       write(out_unitp,*) 'd0f_bhe',d0f_bhe
!       -----------------------------------------------------

        CALL calc_nDindex(Basis2n%nDindG,i_point,ind_quadra(:),err_sub)

        IF (err_sub /= 0) THEN
           write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) '  from Basis2n%nDindG'
           STOP 'calc_nDindex'
         END IF

!       -----------------------------------------------------
!       determine les d0x en fonction de i_point
        CALL Rec_x(d0x,Basis2n,i_point)
!       -----------------------------------------------------

!       -----------------------------------------------------
!       determine Qinact et deltaQ en fonction de d0Qeq d0x
        CALL calc_d0xTOQ(Qinact,dnQeq%d0,deltaQ,d0x,dnC_inv%d0,nb_inact2n)

!       -----------------------------------------------------
!       determine les d1x d2x en fonction de i_point
        CALL calc_d1d2x(d1x,d2x,nb_inact2n,                             &
                        deltaQ,dnC%d0,dnQeq%d1,dnC%d1,dnQeq%d2,dnC%d2,nb_act1)
!       -----------------------------------------------------

!       -- merge Qact(nb_var) (active and rigid) and Qinact(nb_inact2n)
!       ---here only nb_inact2n variables have been modified --------
        CALL Qinact2n_TO_Qact_FROM_ActiveTransfo(Qinact,Qact,mole%ActiveTransfo)

        CALL get_d0MatOp_AT_Qact(Qact,d0MatOp,mole,                     &
                                 para_Tnum,para_AllOp%tab_Op(1)%para_PES)

        IF (.NOT. pot) THEN ! remove the potential part
          d0MatOp(1)%ReVal(:,:,1) = ZERO
        END IF


        IF (KEO) THEN
          IF (para_Tnum%KEO_TalyorOFQinact2n == -1) THEN
            CALL calc3_f2_f1Q_num(Qact,                                 &
                                  f2Q,f1Q,vep,rho,Tcor2,tcor1a,trota,   &
                                  para_Tnum,mole)
          ELSE
            CALL calc3_f2_f1Q_numTay0Qinact2n(Qact,dnQeq,               &
                                    f2Q,f1Q,vep,rho,Tcor2,tcor1a,trota, &
                                    para_Tnum,mole)
          END IF
        ELSE
          vep        = ZERO
          f1Q(:)     = ZERO
          f2Q(:,:)   = ZERO
          Tcor2(:,:) = ZERO
          Tcor1a(:)  = ZERO
          Trota(:,:) = ZERO
        END IF


!       ------------------------------------------------------------
!       -- split f2Q in f2Qaa f2Qii F2Qai and f1Q in f1Qa f1Qi -----
        f2Qaa(:,:) = f2Q(1:nb_act1,1:nb_act1)
        f2Qii(:,:) = f2Q(nb_act1+1:nb_act1+nb_inact2n,                  &
                        nb_act1+1:nb_act1+nb_inact2n)
        f2Qai(:,:) = f2Q(1:nb_act1,nb_act1+1:nb_act1+nb_inact2n)

        f1Qa(:) = f1Q(1:nb_act1)
        f1Qi(:) = f1Q(nb_act1+1:nb_act1+nb_inact2n)

!       -- split Tcor2 in Tcor2a Tcor2i
        tcor2a(:,:)=Tcor2(1:nb_act1,:)
        tcor2i(:,:)=Tcor2(nb_act1+1:nb_act1+nb_inact2n,:)
!       ------------------------------------------------------------


!        -----------------------------------------------------
         Vinact = vep
         IF (pot) THEN
           Vinact = Vinact + d0MatOp(1)%ReVal(1,1,1)
           IF (para_AllOp%tab_Op(1)%para_PES%HarD) THEN
             Vinact = Vinact + pot2(dnHess%d0,deltaQ,nb_inact2n)
             Vinact = Vinact + pot_rest(Qact,deltaQ,nb_inact2n)
           END IF

         END IF
!        Scalar operators (Dip)
         DO iOp=3,nb_Op
           ScalOp(iOp-2) = d0MatOp(iOp)%ReVal(1,1,1)
         END DO

!        -----------------------------------------------------

        !write(out_unitp,*) 'Qact,Ta',Qact(1:nb_act1),f2Qaa,f1Qa,vep
        !write(out_unitp,*) 'Qact,Tii',Qact(1:nb_act1),f2Qii,f1Qi
        !write(out_unitp,*) 'Qact,Tai',Qact(1:nb_act1),f2Qai
        !write(out_unitp,*) 'Qact,V',Qact(1:nb_act1),Vinact


         CALL sub_mat6_HST(para_AllOp%tab_Op(1)%para_PES,               &
                           d0MatHADAOp,nb_Op,nb_bie,                    &
                           d0f_bhe,                                     &
                           rho,                                         &
                           nb_inact2n,ind_quadra,                       &
                           d1x,d2x,                                     &
                           dnC%d0,dnC%d1,dnC%d2,                        &
                           d0cd0c,                                      &
                           dnLnN%d1,dnLnN%d2,                           &
                           f2Qaa,f2Qii,f2Qai,                           &
                           f1Qa,f1Qi,nb_act1,                           &
                           para_AllOp%tab_Op(1)%para_Tnum%JJ,           &
                           tcor2a,tcor2i,tcor1a,trota,                  &
                           Basis2n,wnDh,Vinact,ScalOp)

!       -----------------------------------------------------

      END DO
      IF (print_level > 0 .AND. nq > nq_write_HADA)                     &
                           write(out_unitp,'(a)',ADVANCE='yes') ' - 100'
      CALL flush_perso(out_unitp)

!     --- END Gauss-Hermite quadrature loop ---------------------
!     -------------------------------------------------------
!     -------------------------------------------------------
!     -------------------------------------------------------
!------------------------------------------------------------
!------------------------------------------------------------


! since the imaginary part is only on the active coordinates (Qact1)
      IF (d0MatHADAOp(1)%cplx) THEN
         DO ihe=1,nb_bie
           d0MatHADAOp(1)%ImVal(ihe,ihe) = d0MatOp(1)%ImVal(1,1)
         END DO
      END IF


!------------------------------------------------------------
!------------------------------------------------------------
!     -------------------------------------------------------
      IF (test .OR. debug) THEN
        DO iOp=1,nb_Op
          write(out_unitp,*) ' iOp, name_Op:',iOp,                      &
                                          para_AllOp%tab_Op(iOp)%name_Op
          DO k_term=1,d0MatHADAOp(iOp)%nb_term
            IF (para_AllOp%tab_Op(1)%para_Tnum%JJ == 0 .AND.            &
                count(d0MatHADAOp(iOp)%derive_termQact(:,k_term) <0) > 0) CYCLE
            write(out_unitp,*) ' name_Op, deriv_term:',                 &
                             para_AllOp%tab_Op(iOp)%name_Op,            &
                             d0MatHADAOp(iOp)%derive_termQact(:,k_term)
            CALL Write_Mat(d0MatHADAOp(iOp)%ReVal(:,:,k_term),out_unitp,5)
          END DO
          IF (d0MatHADAOp(iOp)%cplx) THEN
            write(out_unitp,*) ' cplx name_Op:',para_AllOp%tab_Op(iOp)%name_Op
            CALL Write_Mat(d0MatHADAOp(iOp)%ImVal(:,:),out_unitp,5)
          END IF
        END DO
        CALL flush_perso(out_unitp)
      END IF

!     ------ free memory -----------------------------------------
      CALL dealloc_NParray(d0f_bhe,"d0f_bhe",name_sub)


      END IF

!     ------ free memory -----------------------------------------

      DO iOp=1,size(d0MatOp)
        CALL dealloc_d0MatOp(d0MatOp(iOp)) ! Scalar Operator
      END DO

      CALL dealloc_dnSVM(dnC)
      CALL dealloc_dnSVM(dnC_inv)
      CALL dealloc_dnSVM(dnQeq)
      CALL dealloc_dnSVM(dnEHess)
      CALL dealloc_dnSVM(dnHess)
      CALL dealloc_dnSVM(dnGrad)
      CALL dealloc_dnSVM(dnLnN)

      CALL dealloc_NParray(f1Q,  "f1Q",   name_sub)
      CALL dealloc_NParray(f2Q,  "f2Q",   name_sub)
      CALL dealloc_NParray(Tcor2,"Tcor2", name_sub)
      CALL dealloc_NParray(f1Qa, "f1Qa",  name_sub)
      CALL dealloc_NParray(f1Qi, "f1Qi",  name_sub)
      CALL dealloc_NParray(f2Qaa,"f2Qaa", name_sub)
      CALL dealloc_NParray(f2Qii,"f2Qii", name_sub)
      CALL dealloc_NParray(f2Qai,"f2Qai", name_sub)
      CALL dealloc_NParray(tcor2a,"tcor2a",name_sub)
      CALL dealloc_NParray(tcor2i, "tcor2i", name_sub)
      CALL dealloc_NParray(tcor1a, "tcor1a", name_sub)
      CALL dealloc_NParray(trota,  "trota",  name_sub)

      CALL dealloc_NParray(d0x,    "d0x",    name_sub)
      CALL dealloc_NParray(d1x,    "d1x",    name_sub)
      CALL dealloc_NParray(d2x,    "d2x",    name_sub)

      CALL dealloc_NParray(d0cd0c, "d0cd0c", name_sub)

      CALL dealloc_NParray(Qinact, "Qinact", name_sub)
      CALL dealloc_NParray(deltaQ2,"deltaQ2",name_sub)
      CALL dealloc_NParray(deltaQ, "deltaQ", name_sub)
!     ------------------------------------------------------------

!     -------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -------------------------------------------------------


     END SUBROUTINE sub_HST7_bhe

