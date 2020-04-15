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
MODULE Mod_MatOp
 IMPLICIT NONE
 PRIVATE
 PUBLIC :: sub_MatOp,sub_build_MatOp
CONTAINS

!=======================================================================================
      SUBROUTINE sub_MatOp(para_Op,print_Op)
      USE mod_system
      USE mod_Constant
      USE mod_SetOp
      USE mod_MPI   
      IMPLICIT NONE

!----- Operator variables --------------------------------------------
      TYPE (param_Op)  :: para_Op


!------ for the H matrix analysis -------------------------
      logical           :: print_Op
      real (kind=Rkind) :: non_hermitic,SS,max_Sii,max_Sij,auTOcm_inv

!----- divers ----------------------------------------------------
      integer  :: n,ib,jb,iq
      !logical  :: test = .TRUE.
      logical  :: test = .FALSE.

!----- for debuging --------------------------------------------------
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (para_Op%mat_done) RETURN
      IF (.NOT. para_Op%alloc_mat)                                      &
             CALL alloc_para_Op(para_Op,Mat=.TRUE.,Grid=.FALSE.)

      n = para_Op%nb_tot
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING sub_MatOp',para_Op%n_Op
        write(out_unitp,*) 'nb_act1,nb_var',para_Op%mole%nb_act1,       &
                                   para_Op%mole%nb_var
        write(out_unitp,*) 'read MatOp',para_Op%n_Op,para_Op%read_Op
        write(out_unitp,*) 'nb_tot',para_Op%nb_tot
        write(out_unitp,*)
        CALL write_param_Op(para_Op)
!        IF (associated(para_Op%Cmat)) CALL Write_Mat(para_Op%Cmat,out_unitp,3)
!        IF (associated(para_Op%Rmat)) CALL Write_Mat(para_Op%Rmat,out_unitp,5)
      END IF
!-----------------------------------------------------------

      para_Op%Make_mat = .FALSE.

      IF (para_Op%read_Op) THEN
        CALL sub_Read_MatOp(para_Op)
      ELSE
        IF (para_Op%para_ReadOp%para_FileGrid%Type_FileGrid == 0 .AND.  &
               .NOT. para_Op%para_ReadOp%para_FileGrid%Save_MemGrid) THEN
          IF (para_Op%Para_Tnum%JJ == 0) THEN
            CALL sub_MatOp_WITH_FileGrid_type0(para_Op)
            !CALL sub_MatOpVibRot_WITH_FileGrid_type0(para_Op)
          ELSE
            CALL sub_MatOpVibRot_WITH_FileGrid_type0(para_Op)
          END IF
        ELSE
          IF (test) THEN
            !CALL time_perso('sub_MatOp_OpExact_SG4')
            !CALL sub_MatOp_OpExact_SG4(para_Op)
            !CALL time_perso('sub_MatOp_OpExact_SG4')

            CALL time_perso('sub_MatOp_direct1_Overlap')
            CALL sub_MatOp_direct1_Overlap(para_Op)
            CALL sub_ana_S(para_Op%Rmat,para_Op%nb_tot,max_Sii,max_Sij,.TRUE.)

            CALL time_perso('sub_MatOp_direct1_Overlap')

            !CALL sub_MatOp_Overlap_SG4(para_Op)
            !CALL sub_MatOp_V_SG4(para_Op)
            STOP
          END IF
#if(run_MPI)          
          IF(openmpi) THEN
            ! add MPI for sub_MatOp_direct1 later
            CALL sub_MatOp_direct1(para_Op)
          ELSE
#endif
            IF (MatOp_omp == 2) THEN
              CALL sub_MatOp_direct2(para_Op) ! for openmp
            ELSE IF (MatOp_omp == 1) THEN
              CALL sub_MatOp_direct1_old(para_Op)  ! for openmp but more memory
            ELSE ! no openmp (nb_thread=1)
              CALL sub_MatOp_direct1(para_Op)
            END IF
#if(run_MPI)             
          ENDIF
#endif
        END IF
      END IF
      para_Op%Make_mat = .TRUE.
      para_Op%mat_done = .TRUE.

      IF (debug) THEN
        write(out_unitp,*) para_Op%name_Op,' non symmetrized'
        IF (para_Op%cplx) THEN
          CALL Write_Mat(para_Op%Cmat,out_unitp,3)
        ELSE
          CALL Write_Mat(para_Op%Rmat,out_unitp,5)
        END IF
      END IF

!     - For the Hamiltonian (n_Op=0) -------------------------
      IF (para_Op%n_Op == 0 ) THEN
        !IF (print_level>-1 .AND. MPI_id==0)                                            &
        !                    write(out_unitp,*) 'Hmin and Hmax',para_Op%Hmin,para_Op%Hmax
        IF (para_Op%cplx) THEN
          CALL sub_hermitic_cplxH(para_Op%Cmat,para_Op%nb_tot,          &
                                  non_hermitic,para_Op%sym_Hamil)
        ELSE
          CALL sub_hermitic_H(para_Op%Rmat,para_Op%nb_tot,non_hermitic,para_Op%sym_Hamil)
        END IF

        auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')
        IF (non_hermitic >= FOUR/TEN**4) THEN
          If(MPI_id==0) write(out_unitp,*) 'WARNING: non_hermitic is BIG'
          If(MPI_id==0) write(out_unitp,31) non_hermitic
 31       format(' Hamiltonien: ',f16.12,' au')
        ELSE
          IF (print_level>-1) write(out_unitp,21) non_hermitic*auTOcm_inv
 21       format(' Hamiltonien: ',f16.12,' cm-1')
        END IF

      END IF

      IF ((print_Op .OR. debug) .AND. MPI_id==0) THEN
        write(out_unitp,*) para_Op%name_Op,' symmetrized'
        IF (para_Op%cplx) THEN
          CALL Write_Mat(para_Op%Cmat,out_unitp,3)
        ELSE
          CALL Write_Mat(para_Op%Rmat,out_unitp,5)
        END IF
      END IF
!     --------------------------------------------------------
!     - for sthe spectral representation -------------------
      IF (para_Op%spectral) CALL sub_Spectral_Op(para_Op)
!     --------------------------------------------------------

!     - pack the operator ------------------------------------
      IF (para_Op%pack_Op) CALL pack_MatOp(para_Op)
!     --------------------------------------------------------


      IF (para_Op%spectral .AND. (print_Op .OR. debug)) THEN
        IF(MPI_id==0) write(out_unitp,*) para_Op%name_Op,' spectral'
        IF (para_Op%cplx) THEN
          CALL Write_Mat(para_Op%Cmat,out_unitp,3)
        ELSE
          CALL Write_Mat(para_Op%Rmat,out_unitp,5)
        END IF
      END IF

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*)
         write(out_unitp,*) 'END sub_MatOp'
       END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_MatOp
!=======================================================================================


!=======================================================================================
!     Construct an Operator on a set of WP
!=======================================================================================
      SUBROUTINE sub_build_MatOp(WP,nb_WP,para_Op,hermitic,print_mat)
      USE mod_system
      USE mod_Constant
      USE mod_psi_set_alloc
      USE mod_psi_Op
      USE mod_SetOp
      USE mod_OpPsi
      IMPLICIT NONE


!----- variables for the WP propagation ----------------------------
      integer               :: nb_WP
      TYPE (param_psi)      :: WP(nb_WP)
      TYPE (param_Op)       :: para_Op
      logical               :: print_mat,hermitic

!------ working parameters --------------------------------
      TYPE (param_psi), allocatable     :: OpWP(:)
      real (kind=Rkind), allocatable    :: RMatOp(:,:)
      complex (kind=Rkind), allocatable :: CMatOp(:,:)

      complex (kind=Rkind)  :: Overlap
      real (kind=Rkind)     :: non_hermitic,auTOcm_inv
      integer       :: i,j
      integer       :: err
      character (len=Name_len) :: info
      logical                  :: spectral_save

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter ::name_sub='sub_build_MatOp'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) ' nb_diago',nb_WP
        write(out_unitp,*)
        CALL write_param_Op(para_Op)
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------
      spectral_save = para_Op%spectral
      para_Op%spectral = .FALSE.
!     - init and allocation of OpWP --
      CALL alloc_NParray(OpWP,(/nb_WP/),"OpWP",name_sub)

      IF (para_Op%cplx) THEN
        CALL alloc_NParray(CMatOp,(/nb_WP,nb_WP/),'CmatOp',name_sub)
      ELSE
        CALL alloc_NParray(RMatOp,(/nb_WP,nb_WP/),'RmatOp',name_sub)
      END IF
      DO i=1,nb_WP
        CALL copy_psi2TOpsi1(OpWP(i),WP(i),alloc=.FALSE.)
        !@chen was CALL init_psi(OpWP(i),para_Op,para_Op%cplx)
      END DO

!     - calculation of the matrix ----
      CALL sub_TabOpPsi(WP,OpWP,para_Op)
      DO i=1,nb_WP
        CALL sub_scaledOpPsi(WP(i),OpWP(i),para_Op%E0,ONE)
        DO j=1,nb_WP
          CALL Overlap_psi1_psi2(Overlap,WP(j),OpWP(i))
          IF (para_Op%cplx) THEN
            CMatOp(i,j) = Overlap
          ELSE
            RMatOp(i,j) = real(Overlap,kind=Rkind)
          END IF
        END DO
      END DO

      IF (hermitic .AND. para_Op%n_Op == 0) THEN
        IF (para_Op%cplx) THEN
          CALL sub_hermitic_cplxH(CMatOp,nb_WP,non_hermitic,para_Op%sym_Hamil)
          IF (.NOT. associated(para_Op%Cdiag)) THEN
            CALL alloc_array(para_Op%Cdiag,(/nb_WP/),"para_Op%Cdiag",name_sub)
          END IF
          IF (.NOT. associated(para_Op%Cvp)) THEN
            CALL alloc_array(para_Op%Cvp,(/nb_WP,nb_WP/),"para_Op%Cvp",name_sub)
          END IF
          para_Op%Cvp(:,:) = CZERO
          DO i=1,nb_WP
            para_Op%Cdiag(i) = CMatOp(i,i)
            para_Op%Cvp(i,i) = CONE
          END DO
        ELSE
          CALL sub_hermitic_H(RMatOp,nb_WP,non_hermitic,para_Op%sym_Hamil)
          IF (.NOT. associated(para_Op%Rdiag)) THEN
            CALL alloc_array(para_Op%Rdiag,(/nb_WP/),"para_Op%Rdiag",name_sub)
          END IF
          IF (.NOT. associated(para_Op%Rvp)) THEN
            CALL alloc_array(para_Op%Rvp,(/nb_WP,nb_WP/),"para_Op%Rvp",name_sub)
          END IF
          para_Op%Rvp(:,:) = ZERO
          DO i=1,nb_WP
            para_Op%Rdiag(i) = RMatOp(i,i)
            para_Op%Rvp(i,i) = ONE
          END DO
        END IF
        IF (non_hermitic >= FOUR/TEN**4) THEN
          If(MPI_id==0) write(out_unitp,*) 'WARNING: non_hermitic is BIG'
          If(MPI_id==0) write(out_unitp,31) non_hermitic
 31       format(' non-hermitic Hamiltonien: ',f16.12,' au')
        ELSE
          auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')
          If(MPI_id==0) write(out_unitp,21) non_hermitic*auTOcm_inv
 21       format(' non-hermitic Hamiltonien: ',f16.12,' cm-1')
        END IF
        CALL flush_perso(out_unitp)
      END IF

!     - Write the matrix ----
      IF (print_mat .OR. debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) '==== Write Op: ',para_Op%name_Op
        IF (para_Op%cplx) THEN
          CALL Write_Mat(CMatOp,out_unitp,5)
        ELSE
          CALL Write_Mat(RMatOp,out_unitp,5)
        END IF
        CALL flush_perso(out_unitp)
      END IF


!     - Save the matrix for a spectral representation
      para_Op%spectral    = .TRUE.
      para_Op%spectral_Op = 0
      para_Op%nb_tot      = nb_WP
      para_Op%mat_done    = .TRUE.
      write(out_unitp,*) 'nb_tot_ini',para_Op%nb_tot_ini
      CALL flush_perso(out_unitp)

      IF (para_Op%cplx) THEN
        IF (associated(para_Op%Cmat))  THEN
          CALL dealloc_array(para_Op%Cmat,"para_Op%Cmat",name_sub)
        END IF
        CALL alloc_array(para_Op%Cmat,(/nb_WP,nb_WP/),'para_Op%Cmat',name_sub)
        para_Op%Cmat(:,:) = CMatOp(:,:)
      ELSE
        IF (associated(para_Op%Rmat))  THEN
          CALL dealloc_array(para_Op%Rmat,"para_Op%Rmat",name_sub)
        END IF
        CALL alloc_array(para_Op%Rmat,(/nb_WP,nb_WP/),'para_Op%Rmat',name_sub)
        para_Op%Rmat(:,:) = RMatOp(:,:)
      END IF

!     - free memories -------
      DO i=1,nb_WP
        CALL dealloc_psi(OpWP(i))
      END DO
      CALL dealloc_NParray(OpWP,"OpWP",name_sub)

      IF (allocated(CMatOp)) THEN
        CALL dealloc_NParray(CMatOp,"CMatOp",name_sub)
      END IF
      IF (allocated(RMatOp)) THEN
        CALL dealloc_NParray(RMatOp,"RMatOp",name_sub)
      END IF

      para_Op%spectral = spectral_save

!----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'END ',name_sub
         CALL flush_perso(out_unitp)
       END IF
!----------------------------------------------------------

      END SUBROUTINE sub_build_MatOp
!=======================================================================================


!=====================================================================
!
!  analysis of a matrix M(n,n) => dim_M(i) and ind_M(:,i)
!
!=====================================================================
      SUBROUTINE pack_MatOp(para_Op)
      USE mod_system
      USE mod_SetOp
      IMPLICIT NONE



!----- Operator variables --------------------------------------------
      TYPE (param_Op)  :: para_Op


      real (kind=Rkind) :: rate
      integer :: i,j,jj

!----- for debuging ----------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING pack_MatOp'
        write(out_unitp,*) 'nb_tot',para_Op%nb_tot
        write(out_unitp,*) 'pack_Op,tol_pack',para_Op%pack_Op,para_Op%tol_pack
      END IF
!-----------------------------------------------------------

      IF (.NOT. para_Op%mat_done) para_Op%pack_Op = .FALSE.

      IF (para_Op%pack_Op) THEN

        CALL alloc_para_Op(para_Op)

        IF (para_Op%cplx) THEN
          DO i=1,para_Op%nb_tot
            para_Op%dim_Op(i) =                                         &
                          count(abs(para_Op%Cmat(:,i))>para_Op%tol_pack)
            jj = 0
            DO j=1,para_Op%nb_tot
              IF ( abs(para_Op%Cmat(j,i))>para_Op%tol_pack) THEN
                jj = jj + 1
                para_Op%ind_Op(jj,i) = j
              END IF
            END DO
          END DO

        ELSE

          DO i=1,para_Op%nb_tot
            para_Op%dim_Op(i) =                                         &
                   count(abs(para_Op%Rmat(:,i))>para_Op%tol_pack)
            jj = 0
            DO j=1,para_Op%nb_tot
              IF ( abs(para_Op%Rmat(j,i))>para_Op%tol_pack) THEN
                jj = jj + 1
                para_Op%ind_Op(jj,i) = j
              END IF
            END DO
          END DO


        END IF

      END IF

      para_Op%ratio_pack = real(sum(para_Op%dim_Op(:)),kind=Rkind) /    &
                                      real(para_Op%nb_tot,kind=Rkind)**2
      write(out_unitp,*) 'pack_mat: Sum(dim)',sum(para_Op%dim_Op(:))
      write(out_unitp,*) 'pack_mat: Sum(dim)/n^2',para_Op%ratio_pack

      IF (para_Op%ratio_pack >= para_Op%tol_nopack) THEN
        write(out_unitp,*) 'pack_mat: ratio >= tol_nopack',             &
                     para_Op%ratio_pack,para_Op%tol_nopack
        write(out_unitp,*) 'pack_mat: rate too large => the Op is NOT packed'
        CALL dealloc_array(para_Op%dim_Op,"para_Op%dim_Op","pack_MatOp")
        CALL dealloc_array(para_Op%ind_Op,"para_Op%ind_Op","pack_MatOp")
      END IF


!     -------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        DO i=1,para_Op%nb_tot
          write(out_unitp,*) i,para_Op%dim_Op(i),': ',                          &
                 para_Op%ind_Op(1:para_Op%dim_Op(i),i)
        END DO
        write(out_unitp,*) 'END pack_MatOp'
      END IF
!     -------------------------------------------------------

      END SUBROUTINE pack_MatOp
!================================================================
!
!     The Operator are transformed into the spectral representation
!
!================================================================
      SUBROUTINE sub_Spectral_Op(para_Op)
      USE mod_system
      USE mod_SetOp
      IMPLICIT NONE

!----- Operator variables --------------------------------------------
      TYPE (param_Op)  :: para_Op


!----- for the gestion of the memory ---------------------------------
      real (kind=Rkind), allocatable    :: Rmat1(:,:)
      real (kind=Rkind), allocatable    :: Rmat2(:,:)
      complex (kind=Rkind), allocatable :: Cmat1(:,:)
      complex (kind=Rkind), allocatable :: Cmat2(:,:)

!----- divers ----------------------------------------------------
      integer  :: n
      integer  :: nio
      integer  :: i,j,k,l

!----- for debuging --------------------------------------------------
      integer :: err,err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='sub_Spectral_Op'
!-----------------------------------------------------------

      IF (.NOT. para_Op%alloc) CALL alloc_para_Op(para_Op)
      IF (.NOT. para_Op%mat_done) RETURN

      n = para_Op%nb_tot
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub,para_Op%n_Op
        write(out_unitp,*) 'nb_act1,nb_var',para_Op%mole%nb_act1,               &
                                   para_Op%mole%nb_var
        write(out_unitp,*) 'nb_tot',para_Op%nb_tot
        write(out_unitp,*)
        CALL write_param_Op(para_Op)
!       IF (associated(para_Op%Cmat)) CALL Write_Mat(para_Op%Cmat,out_unitp,3)
!       IF (associated(para_Op%Rmat)) CALL Write_Mat(para_Op%Rmat,out_unitp,5)
      END IF
!-----------------------------------------------------------

      IF (.NOT. para_Op%mat_done) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' the Matrix of the operator is NOT calculated !'
        write(out_unitp,*) ' mat_done',para_Op%mat_done
        CALL write_param_Op(para_Op)
        STOP
      END IF


!     - spectral representation ------------------------------
      write(out_unitp,*) 'spectral',para_Op%spectral,para_Op%spectral_Op
      IF ( para_Op%spectral .AND. para_Op%n_Op == para_Op%spectral_Op ) THEN
        para_Op%diago = .TRUE.
        CALL alloc_para_Op(para_Op)

        IF (para_Op%cplx) THEN
          CALL sub_diago_CH(para_Op%Cmat,para_Op%Cdiag,para_Op%Cvp,     &
                          para_Op%nb_baie)
          para_Op%ComOp%nb_vp_spec = min(para_Op%ComOp%nb_vp_spec,      &
                                         para_Op%nb_baie)

          para_Op%nb_tot = para_Op%ComOp%nb_vp_spec
          IF (associated(para_Op%Cmat)) THEN
            CALL dealloc_array(para_Op%Cmat,"para_Op%Cmat",name_sub)
          END IF
          CALL alloc_array(para_Op%Cmat,(/para_Op%nb_tot,para_Op%nb_tot/),&
                          "para_Op%Cmat",name_sub)
          para_Op%Cmat(:,:) = CZERO

          IF (associated(para_Op%ComOp%liste_spec)) THEN
            DO i=1,para_Op%nb_tot
              para_Op%Cmat(i,i) =                                       &
                              para_Op%Cdiag(para_Op%ComOp%liste_spec(i))
            END DO
          ELSE
            DO i=1,para_Op%nb_tot
              para_Op%Cmat(i,i) = para_Op%Cdiag(i)
            END DO
          END IF
          para_Op%ComOp%Cvp_spec => para_Op%Cvp
        ELSE

          CALL sub_diago_H(para_Op%Rmat,para_Op%Rdiag,para_Op%Rvp,      &
                           para_Op%nb_baie,para_Op%sym_Hamil)

          para_Op%ComOp%nb_vp_spec = min( para_Op%ComOp%nb_vp_spec,     &
                                            para_Op%nb_baie)

          para_Op%nb_tot = para_Op%ComOp%nb_vp_spec
          IF (associated(para_Op%Rmat)) THEN
            CALL dealloc_array(para_Op%Rmat,"para_Op%Rmat",name_sub)
          END IF
          CALL alloc_array(para_Op%Rmat,(/para_Op%nb_tot,para_Op%nb_tot/),&
                          "para_Op%Rmat",name_sub)
          para_Op%Rmat(:,:) = ZERO


          IF (associated(para_Op%ComOp%liste_spec)) THEN
            DO i=1,para_Op%nb_tot
              para_Op%Rmat(i,i) =                                       &
                      para_Op%Rdiag(para_Op%ComOp%liste_spec(i))
            END DO
          ELSE
            DO i=1,para_Op%nb_tot
              para_Op%Rmat(i,i) = para_Op%Rdiag(i)
            END DO
          END IF
          para_Op%ComOp%Rvp_spec => para_Op%Rvp
        END IF
        IF (para_Op%pack_Op)  THEN
          CALL dealloc_array(para_Op%ind_Op,"para_Op%ind_Op",name_sub)
          CALL dealloc_array(para_Op%dim_Op,"para_Op%dim_Op",name_sub)
        END IF

      ELSE IF ( para_Op%spectral .AND. para_Op%n_Op /= para_Op%spectral_Op ) THEN
        write(out_unitp,*) 'para_Op%ComOp%nb_vp_spec',para_Op%ComOp%nb_vp_spec
        IF (para_Op%cplx) THEN

          CALL alloc_NParray(Cmat1,(/para_Op%nb_tot,para_Op%nb_tot/),     &
                          'Cmat1','sub_Spectral_Op')
          CALL alloc_NParray(Cmat2,(/para_Op%nb_tot,para_Op%nb_tot/),     &
                          'Cmat2','sub_Spectral_Op')

          Cmat1(:,:) = matmul(para_Op%Cmat,para_Op%ComOp%Cvp_spec)
          Cmat2(:,:) = transpose(para_Op%ComOp%Cvp_spec)
          para_Op%Cmat(:,:) = matmul(Cmat2,Cmat1)
          Cmat1(:,:) = para_Op%Cmat(:,:)
          para_Op%nb_tot = para_Op%ComOp%nb_vp_spec
          CALL dealloc_array(para_Op%Cmat,'para_Op%Cmat',name_sub)
          CALL alloc_array(para_Op%Cmat,(/para_Op%nb_tot,para_Op%nb_tot/),&
                          'para_Op%Cmat',name_sub)

          para_Op%Cmat(:,:) =                                           &
              Cmat1(para_Op%ComOp%liste_spec,para_Op%ComOp%liste_spec)

          CALL dealloc_NParray(Cmat1,'Cmat1',name_sub)
          CALL dealloc_NParray(Cmat2,'Cmat2',name_sub)
        ELSE
          CALL alloc_NParray(Rmat1,(/para_Op%nb_tot,para_Op%nb_tot/),     &
                          'Rmat1',name_sub)
          CALL alloc_NParray(Rmat2,(/para_Op%nb_tot,para_Op%nb_tot/),     &
                          'Rmat2',name_sub)

          Rmat1(:,:) = matmul(para_Op%Rmat,para_Op%ComOp%Rvp_spec)
          Rmat2(:,:) = transpose(para_Op%ComOp%Rvp_spec)
          para_Op%Rmat(:,:) = matmul(Rmat2,Rmat1)
          Rmat1(:,:) = para_Op%Rmat(:,:)
          para_Op%nb_tot = para_Op%ComOp%nb_vp_spec
          CALL dealloc_array(para_Op%Rmat,'para_Op%Rmat',name_sub)
          CALL alloc_array(para_Op%Rmat,(/para_Op%nb_tot,para_Op%nb_tot/),&
                          'para_Op%Rmat',name_sub)

          para_Op%Rmat(:,:) =                                           &
              Rmat1(para_Op%ComOp%liste_spec,para_Op%ComOp%liste_spec)

          CALL dealloc_NParray(Rmat1,'Rmat1',name_sub)
          CALL dealloc_NParray(Rmat2,'Rmat2',name_sub)

        END IF
        IF (para_Op%pack_Op)  THEN
          CALL dealloc_array(para_Op%ind_Op,"para_Op%ind_Op",name_sub)
          CALL dealloc_array(para_Op%dim_Op,"para_Op%dim_Op",name_sub)
        END IF
      END IF
!     --------------------------------------------------------


!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*)
         write(out_unitp,*) 'END ',name_sub
       END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_Spectral_Op
!====================================================================
!
!     para_Op%para_ReadOp%para_FileGrid%Type_FileGrid == 0
!      and
!     para_Op%para_ReadOp%para_FileGrid%Save_MemGrid = .F.
!
!=====================================================================
      SUBROUTINE sub_MatOp_WITH_FileGrid_type0(para_Op)
      USE mod_system
      USE mod_nDindex
      USE mod_SetOp
      IMPLICIT NONE

!----- Operator variables --------------------------------------------
      TYPE (param_Op)  :: para_Op


!------ quadrature points and weight -----------------------------
      real (kind=Rkind) :: WnD

      real (kind=Rkind) :: Qdyn(para_Op%mole%nb_var)
      real (kind=Rkind) :: Qact(para_Op%mole%nb_act1)
      real (kind=Rkind), allocatable :: mat1Im(:,:)
      real (kind=Rkind), allocatable :: mat1(:,:)
      real (kind=Rkind), allocatable :: mat2(:,:)
      real (kind=Rkind), allocatable :: mat3(:,:)


!------ for td0b ...         -------------------------------------
      real (kind=Rkind), allocatable    :: VecQ(:)
      integer                           :: kmem
      real (kind=Rkind), allocatable    :: td0b(:,:)
      TYPE (param_d0MatOp), allocatable :: d0MatOpd0bWrho(:,:)
      TYPE (param_d0MatOp)              :: d0MatOp
      integer                           :: type_Op
      integer                           :: iterm_Op,iterm
!----- for the  memory ---------------------------------
      integer                       :: nreste,nplus,nb_blocks

!----- divers ----------------------------------------------------
      integer  :: nb_ba,nb_bie
      integer  :: n

      integer  :: i,k
      integer  :: i1,i2,f1,f2
      integer  :: i1_h,i2_h,i_h,ib1

      real (kind=Rkind) :: Hinter


      integer :: nDGridI

!----- for debuging --------------------------------------------------
      integer :: err,err_mem,memory
      character (len=*), parameter :: name_sub="sub_MatOp_WITH_FileGrid_type0"
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (.NOT. para_Op%alloc_mat)                                      &
             CALL alloc_para_Op(para_Op,Mat=.TRUE.,Grid=.FALSE.)
      IF (para_Op%para_ReadOp%para_FileGrid%Type_FileGrid /= 0 .OR.     &
          para_Op%para_ReadOp%para_FileGrid%Save_MemGrid_done .OR. para_Op%mat_done) RETURN

      para_Op%mat_done = .TRUE.
      n = para_Op%nb_tot
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub,para_Op%n_Op
        write(out_unitp,*) 'nb_act1,nb_var',para_Op%mole%nb_act1,               &
                                   para_Op%mole%nb_var
        write(out_unitp,*) 'nb_tot',para_Op%nb_tot
        write(out_unitp,*)
        CALL write_param_Op(para_Op)
      END IF
!-----------------------------------------------------------

!     ----------------------------------------------------------------
!     - test on nb_bi : nb_bi>0
      IF (para_Op%nb_bi <= 0) THEN
        write(out_unitp,*) ' ERROR : nb_bi MUST be > 0'
        write(out_unitp,*) 'nb_bi',para_Op%nb_bi
        write(out_unitp,*) 'STOP in ',name_sub
        write(out_unitp,*)
        CALL write_param_Op(para_Op)
        STOP
      END IF
!     ----------------------------------------------------------------

      !-------------------------------------------------------------
      !-     memories allocation: td0b, Opd0bWrho mat1, mat2, mat3
      !-------------------------------------------------------------
      nb_ba = para_Op%nb_ba
      nb_bie = para_Op%nb_bie


      IF (para_Op%name_Op == 'H') THEN
        type_Op = para_Op%para_PES%Type_HamilOp ! H
        IF (type_Op /= 1) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '    Type_HamilOp MUST be equal to 1 for HADA or cHAC'
          write(out_unitp,*) '    CHECK your data!!'
          STOP
        END IF
      ELSE
        type_Op = 0
      END IF

      IF (para_Op%cplx) THEN
        CALL alloc_NParray(mat1Im,(/nb_ba,nb_ba/),'mat1Im',name_sub)
      END IF

      CALL alloc_NParray(mat1,(/nb_ba,nb_ba/),'mat1',name_sub)
      IF (para_Op%ComOp%contrac_ba_ON_HAC) THEN
        CALL alloc_NParray(mat2,(/nb_ba,nb_ba/),'mat2',name_sub)
        CALL alloc_NParray(mat3,(/nb_ba,nb_ba/),'mat3',name_sub)
      END IF

      ! selected the optimal value of kmem, as function of max_mem and mem_tot
      !kmem = min(10,para_Op%nb_qa)
      kmem = para_Op%nb_qa/10
      IF (kmem == 0) kmem = para_Op%nb_qa
      allocate(d0MatOpd0bWrho(kmem,nb_ba),stat=err_mem)
      memory = kmem*nb_ba
      CALL error_memo_allo(err_mem,memory,'d0MatOpd0bWrho',name_sub,'d0MatOp')
      DO i=1,nb_ba
      DO k=1,kmem
        CALL Init_d0MatOp(d0MatOpd0bWrho(k,i),type_Op,0,nb_bie,         &
                                            JRot=0,cplx=para_Op%cplx)
      END DO
      END DO

      nb_blocks = para_Op%nb_qa/kmem
      IF (mod(para_Op%nb_qa,kmem) /= 0) nb_blocks = nb_blocks + 1
      IF (print_level>0 .OR. debug) write(out_unitp,*) 'number of blocks',nb_blocks

      CALL alloc_NParray(td0b,(/nb_ba,kmem/),'td0b',name_sub)
      CALL alloc_NParray(VecQ,(/kmem/),'VecQ',name_sub)

      CALL Init_d0MatOp(d0MatOp,para_Op%param_TypeOp,nb_bie)

      !--------------------------------------------------------

      !-- built the Operator matrix ---------------------------
      !--------------------------------------------------------
      IF (para_Op%cplx) THEN
        para_Op%Cmat(:,:) = CZERO
      ELSE
        para_Op%Rmat(:,:) = ZERO
      END IF


!     - Hmin: Vmin and Hmax: the largest diagonal element of H -----
      para_Op%Hmin = huge(ONE)
      para_Op%Hmax = -huge(ONE)


!-------------------------------------------------------------
!     - Real part of the Operator ---------------------------
!-------------------------------------------------------------

!     --------------------------------------------------------
!     - loop of kmem block -----------------------------------
      nreste = para_Op%nb_qa


      nDGridI = 0

 98   CONTINUE

      IF (print_level>-1) write(out_unitp,'(a)',ADVANCE='no') 'Block of ReOpd0b(:,:) (%): '
      CALL flush_perso(out_unitp)

      nplus = min(kmem,nreste)
      DO k=1,nplus

        nDGridI = nDGridI + 1
        IF (mod(k,max(1,int(nplus/10))) == 0 .AND. print_level>-1) THEN
          write(out_unitp,'(a,i3)',ADVANCE='no') ' -',                        &
                      int(real(k,kind=Rkind)*HUNDRED/real(nplus,kind=Rkind))
          CALL flush_perso(out_unitp)
        END IF

        CALL sub_reading_Op(nDGridI,para_Op%nb_qa,                  &
                                d0MatOp,para_Op%n_Op,                   &
                                Qdyn,para_Op%mole%nb_var,Qact,          &
                                WnD,para_Op%ComOp)

        ! -- WARNING: it has to be done before calc_td0b_OpRVd0bW !!
        iterm = d0MatOp%derive_term_TO_iterm(0,0)
        DO i1_h=1,para_Op%nb_bie
          para_Op%Hmin = min(para_Op%Hmin,d0MatOp%ReVal(i1_h,i1_h,iterm))
          para_Op%Hmax = max(para_Op%Hmax,d0MatOp%ReVal(i1_h,i1_h,iterm))
        END DO

        CALL calc_td0b_OpRVd0bW(nDGridI,k,td0b,d0MatOpd0bWrho,          &
                                WnD,kmem,d0MatOp,para_Op,               &
                                para_Op%BasisnD)

      END DO
      IF (print_level>-1) write(out_unitp,'(a)',ADVANCE='yes') ' - End'
      CALL flush_perso(out_unitp)

      !DO iterm_Op=1,d0MatOpd0bWrho(1,1)%nb_term ! one term
      iterm_Op=1

        !================================================================
        ! loop on i1_h, i2_h
        !================================================================
        DO i2_h=1,para_Op%nb_bie
        DO i1_h=1,para_Op%nb_bie

          i1 = sum(para_Op%ComOp%nb_ba_ON_HAC(1:i1_h-1))
          i2 = sum(para_Op%ComOp%nb_ba_ON_HAC(1:i2_h-1))
          f1 = para_Op%ComOp%nb_ba_ON_HAC(i1_h)
          f2 = para_Op%ComOp%nb_ba_ON_HAC(i2_h)

          DO ib1=1,para_Op%nb_ba
            DO k=1,nplus
              VecQ(k) = d0MatOpd0bWrho(k,ib1)%ReVal(i1_h,i2_h,iterm_Op)
            END DO
            mat1(:,ib1) = matmul(td0b(:,1:nplus),VecQ(1:nplus))
          END DO
          !CALL Write_Mat(mat1,out_unitp,5)

          IF (para_Op%ComOp%contrac_ba_ON_HAC) THEN
            mat2 = matmul(mat1,para_Op%ComOp%d0Cba_ON_HAC(:,:,i2_h) )
            mat3 = transpose(para_Op%ComOp%d0Cba_ON_HAC(:,:,i1_h))
            mat1 = matmul(mat3,mat2)
          END IF
          !CALL Write_Mat(mat1,out_unitp,5)

          IF (para_Op%cplx) THEN
            DO ib1=1,para_Op%nb_ba
              DO k=1,nplus
                VecQ(k) = d0MatOpd0bWrho(k,ib1)%ImVal(i1_h,i2_h)
              END DO
              mat1Im(:,ib1) = matmul(td0b(:,1:nplus),VecQ(1:nplus))
            END DO
            !CALL Write_Mat(mat1Im,out_unitp,5)

            IF (para_Op%ComOp%contrac_ba_ON_HAC) THEN
              mat2 = matmul(mat1Im,para_Op%ComOp%d0Cba_ON_HAC(:,:,i2_h) )
              mat3 = transpose(para_Op%ComOp%d0Cba_ON_HAC(:,:,i1_h))
              mat1Im = matmul(mat3,mat2)
            END IF
            !CALL Write_Mat(mat1Im,out_unitp,5)
          END IF


          IF (para_Op%cplx) THEN
            para_Op%Cmat(i1+1:i1+f1 , i2+1:i2+f2) =                     &
                    para_Op%Cmat(i1+1:i1+f1 , i2+1:i2+f2) +             &
                 cmplx(mat1(1:f1 , 1:f2),mat1Im(1:f1 , 1:f2),kind=Rkind)
          ELSE
            para_Op%Rmat(i1+1:i1+f1 , i2+1:i2+f2) =                     &
                  para_Op%Rmat(i1+1:i1+f1 , i2+1:i2+f2) +               &
                              mat1(1:f1 , 1:f2)
          END IF


        END DO
        END DO
        !================================================================
        ! END loop on i1_h, i2_h
        !================================================================
      !END DO

      nreste = nreste - nplus
      IF (nreste .GT. 0) GOTO 98
!     - END of loop of kmem block ----------------------------
!     --------------------------------------------------------
      close(para_Op%ComOp%file_HADA%unit)
      CALL flush_perso(out_unitp)


      !- determination of Hmax --------------------------------
      IF (para_Op%cplx) THEN
        DO i=1,para_Op%nb_tot
          para_Op%Hmax = max(para_Op%Hmax,real(para_Op%Cmat(i,i)))
        END DO
       ELSE
        DO i=1,para_Op%nb_tot
          para_Op%Hmax = max(para_Op%Hmax,para_Op%Rmat(i,i))
        END DO
      END IF
      !--------------------------------------------------------

      !--------------------------------------------------------
      CALL dealloc_NParray(VecQ,     'VecQ',     name_sub)
      CALL dealloc_NParray(td0b,     'td0b',     name_sub)

      DO i=1,nb_ba
      DO k=1,kmem
        CALL dealloc_d0MatOp(d0MatOpd0bWrho(k,i))
      END DO
      END DO
      deallocate(d0MatOpd0bWrho,stat=err_mem)
      memory = kmem*nb_ba
      CALL error_memo_allo(err_mem,memory,'d0MatOpd0bWrho',name_sub,'d0MatOp')

      CALL dealloc_d0MatOp(d0MatOp)

      IF (allocated(mat1Im)) THEN
        CALL dealloc_NParray(mat1Im,'mat1Im',name_sub)
      END IF
      IF (allocated(mat1)) THEN
        CALL dealloc_NParray(mat1,'mat1',name_sub)
      END IF
      IF (allocated(mat2)) THEN
        CALL dealloc_NParray(mat2,'mat2',name_sub)
      END IF
      IF (allocated(mat3)) THEN
        CALL dealloc_NParray(mat3,'mat3',name_sub)
      END IF
      !---------------------------------------------------------------

      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) para_Op%name_Op,' non symetrized'
        IF (para_Op%cplx) THEN
          CALL Write_Mat(para_Op%Cmat,out_unitp,3)
        ELSE
          CALL Write_Mat(para_Op%Rmat,out_unitp,5)
        END IF
        write(out_unitp,*)
      END IF


!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*)
         write(out_unitp,*) 'END ',name_sub
       END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_MatOp_WITH_FileGrid_type0

      SUBROUTINE sub_MatOpVibRot_WITH_FileGrid_type0(para_Op)
      USE mod_system
      USE mod_nDindex
      USE mod_SetOp
      IMPLICIT NONE

!----- Operator variables --------------------------------------------
      TYPE (param_Op)  :: para_Op

!------ quadrature points and weight -----------------------------
      real (kind=Rkind) :: WnD

      real (kind=Rkind) :: Qdyn(para_Op%mole%nb_var)
      real (kind=Rkind) :: Qact(para_Op%mole%nb_act1)

      real (kind=Rkind), allocatable :: mat2(:,:)
      real (kind=Rkind), allocatable :: mat3(:,:)
      real (kind=Rkind), allocatable :: VecQ(:)
      real (kind=Rkind), allocatable :: VecVib(:,:),EneVib(:)


!------ for td0b ...         -------------------------------------
      integer                           :: kmem
      real (kind=Rkind), allocatable    :: td0b(:,:)
      TYPE (param_d0MatOp), allocatable :: d0MatOpd0bWrho(:,:)
      TYPE (param_d0MatOp)              :: MatRV
      integer                           :: type_Op


!----- for the  memory ---------------------------------
      integer                       :: nreste,nplus,nb_blocks

!----- divers ----------------------------------------------------
      integer  :: nb_ba,nb_bie,nb_Op,i_Op

      integer  :: n

      integer  :: i,k,KRot,JRot,ibRot,jbRot
      integer  :: i1,i2,f1,f2
      integer  :: i1_h,i2_h,i_h,ib1
      integer  :: J1,J2,iterm,iterm_BasisRot,iterm_Op,nb_term_BasisRot
      real (kind=Rkind) :: Val_BasisRot


      real (kind=Rkind) :: Hinter

      TYPE (param_d0MatOp) :: d0MatOp

      integer :: nDGridI
      logical :: spectral = .TRUE.
      !logical :: spectral = .FALSE.

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub="sub_MatOpVibRot_WITH_FileGrid_type0"
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (.NOT. para_Op%alloc_mat)                                      &
             CALL alloc_para_Op(para_Op,Mat=.TRUE.,Grid=.FALSE.)
      IF (para_Op%para_ReadOp%para_FileGrid%Type_FileGrid /= 0 .OR.                 &
          para_Op%para_ReadOp%para_FileGrid%Save_MemGrid_done .OR. para_Op%mat_done) RETURN

      para_Op%mat_done = .TRUE.
      n = para_Op%nb_tot
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub,para_Op%n_Op
        write(out_unitp,*) 'nb_act1,nb_var',para_Op%mole%nb_act1,               &
                                   para_Op%mole%nb_var
        write(out_unitp,*) 'nb_tot',para_Op%nb_tot
        write(out_unitp,*)
        !CALL write_param_Op(para_Op)
      END IF
!-----------------------------------------------------------

      JRot = para_Op%Para_Tnum%JJ
      spectral = spectral .AND. (JRot > 0)
      CALL init_RotBasis_Param(para_Op%BasisnD%RotBasis,Jrot)

!     ----------------------------------------------------------------
!     - test on nb_bi : nb_bi>0
      IF (para_Op%nb_bi <= 0) THEN
        write(out_unitp,*) ' ERROR : nb_bi MUST be > 0'
        write(out_unitp,*) 'nb_bi',para_Op%nb_bi
        write(out_unitp,*) 'STOP in ',name_sub
        write(out_unitp,*)
        CALL write_param_Op(para_Op)
        STOP
      END IF
!     ----------------------------------------------------------------


      nb_ba  = para_Op%nb_ba
      nb_bie = para_Op%nb_bie

      IF (para_Op%name_Op == 'H') THEN
        type_Op = para_Op%para_PES%Type_HamilOp ! H
        IF (type_Op /= 1) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '    Type_HamilOp MUST be equal to 1 for HADA or cHAC'
          write(out_unitp,*) '    CHECK your data!!'
          STOP
        END IF
      ELSE
        type_Op = 0
      END IF

!-------------------------------------------------------------
!-     memories allocation: td0b, Opd0bWrho matRV, mat2, mat3
!-------------------------------------------------------------
      CALL Init_d0MatOp(MatRV,type_Op,0,nb_ba,JRot=JRot,cplx=para_Op%cplx)

      IF (para_Op%ComOp%contrac_ba_ON_HAC) THEN
        CALL alloc_NParray(mat2,(/nb_ba,nb_ba/),'mat2',name_sub)
        CALL alloc_NParray(mat3,(/nb_ba,nb_ba/),'mat3',name_sub)
      END IF
      CALL alloc_NParray(VecVib,(/nb_ba,nb_ba/),'VecVib',name_sub)
      CALL alloc_NParray(EneVib,(/nb_ba/),'EneVib',name_sub)


      ! selected the optimal value of kmem, as function of max_mem and mem_tot
      kmem = min(10,para_Op%nb_qa)
      !kmem = para_Op%nb_qa
      allocate(d0MatOpd0bWrho(kmem,nb_ba),stat=err_mem)
      memory = kmem*nb_ba
      CALL error_memo_allo(err_mem,memory,'d0MatOpd0bWrho',name_sub,'d0MatOp')
      DO i=1,nb_ba
      DO k=1,kmem
        CALL Init_d0MatOp(d0MatOpd0bWrho(k,i),type_Op,0,nb_bie,         &
                                            JRot=JRot,cplx=para_Op%cplx)
      END DO
      END DO

      nb_blocks = para_Op%nb_qa/kmem
      IF (mod(para_Op%nb_qa,kmem) /= 0) nb_blocks = nb_blocks + 1
      IF (print_level>0) write(out_unitp,*) 'number of blocks',nb_blocks

      CALL alloc_NParray(td0b,(/nb_ba,kmem/),'td0b',name_sub)
      CALL alloc_NParray(VecQ,(/kmem/),'VecQ',name_sub)

      CALL Init_d0MatOp(d0MatOp,type_Op,para_Op%mole%nb_act1,nb_bie,    &
                                             JRot=JRot,cplx=para_Op%cplx)
      !--------------------------------------------------------


      !-- built the Operator matrix ---------------------------
      !--------------------------------------------------------
      IF (para_Op%cplx) THEN
        para_Op%Cmat(:,:) = CZERO
      ELSE
        para_Op%Rmat(:,:) = ZERO
      END IF


      !- Hmin: Vmin and Hmax: the largest diagonal element of H -----
      para_Op%Hmin =  huge(ONE)
      para_Op%Hmax = -huge(ONE)


      !--------------------------------------------------------
      !- loop of kmem block -----------------------------------
      nreste = para_Op%nb_qa


      nDGridI = 0

 98   CONTINUE
      IF (print_level>-1) write(out_unitp,'(a)',ADVANCE='no') 'Block of ReOpd0b(:,:) (%): '
      CALL flush_perso(out_unitp)

      nplus = min(kmem,nreste)
      DO k=1,nplus

        nDGridI = nDGridI + 1
        IF (mod(k,max(1,int(nplus/10))) == 0 .AND. print_level>-1) THEN
          write(out_unitp,'(a,i3)',ADVANCE='no') ' -',                  &
                           int(real(k,kind=Rkind)*HUNDRED/real(nplus,kind=Rkind))
          CALL flush_perso(out_unitp)
        END IF

        CALL sub_reading_Op(nDGridI,para_Op%nb_qa,                  &
                                d0MatOp,para_Op%n_Op,                   &
                                Qdyn,para_Op%mole%nb_var,Qact,          &
                                WnD,para_Op%ComOp)

        iterm = d0MatOp%derive_term_TO_iterm(0,0)
        DO i1_h=1,para_Op%nb_bie
          para_Op%Hmin = min(para_Op%Hmin,d0MatOp%ReVal(i1_h,i1_h,iterm))
          para_Op%Hmax = max(para_Op%Hmax,d0MatOp%ReVal(i1_h,i1_h,iterm))
        END DO

        CALL calc_td0b_OpRVd0bW(nDGridI,k,td0b,d0MatOpd0bWrho,          &
                                WnD,kmem,d0MatOp,para_Op,               &
                                para_Op%BasisnD)

      END DO


      IF (print_level>-1) write(out_unitp,'(a)',ADVANCE='yes') ' - End'
      CALL flush_perso(out_unitp)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      DO iterm_Op=1,MatRV%nb_term

        !================================================================
        ! loop on i1_h, i2_h
        !================================================================
        DO i2_h=1,para_Op%nb_bie
        DO i1_h=1,para_Op%nb_bie

          DO ib1=1,para_Op%nb_ba
            DO k=1,nplus
              VecQ(k) = d0MatOpd0bWrho(k,ib1)%ReVal(i1_h,i2_h,iterm_Op)
            END DO
            MatRV%ReVal(:,ib1,iterm_Op) = matmul(td0b(:,1:nplus),VecQ(1:nplus))
          END DO
          !write(6,*) '===================================================='
          !write(6,*) 'MatRV%ReVal',iterm_Op,MatRV%derive_termQact(:,iterm_Op)
          !CALL Write_Mat(MatRV%ReVal(:,:,iterm_Op),out_unitp,5)


          ! Rotational contribution
          J1       = MatRV%derive_termQact(1,iterm_Op)
          J2       = MatRV%derive_termQact(2,iterm_Op)
          iterm_BasisRot = para_Op%BasisnD%RotBasis%tab_der_TO_iterm(J1,J2)
          !write(6,*) 'J1,J2',J1,J2,'iterm_Op,iterm_BasisRot',iterm_Op,iterm_BasisRot

          DO ibRot=1,para_Op%nb_bRot
          DO jbRot=1,para_Op%nb_bRot
            Val_BasisRot = para_Op%BasisnD%RotBasis%tab_RotOp(ibRot,jbRot,iterm_BasisRot)
            IF (abs(Val_BasisRot) < ONETENTH**9) CYCLE

            i1 = (ibRot-1)*para_Op%nb_baie + sum(para_Op%ComOp%nb_ba_ON_HAC(1:i1_h-1))
            i2 = (jbRot-1)*para_Op%nb_baie + sum(para_Op%ComOp%nb_ba_ON_HAC(1:i2_h-1))
            f1 = para_Op%ComOp%nb_ba_ON_HAC(i1_h)
            f2 = para_Op%ComOp%nb_ba_ON_HAC(i2_h)

            IF (debug) THEN
              write(6,*) 'J1,J2',J1,J2
              write(6,*) 'i1_h,ibRot,i1+1:i1+f1',i1_h,ibRot,i1+1,i1+f1
              write(6,*) 'i2_h,jbRot,i2+1:i2+f2',i2_h,jbRot,i2+1,i2+f2
            END IF

            IF (para_Op%cplx) THEN
              para_Op%Cmat(i1+1:i1+f1 , i2+1:i2+f2) =                   &
                                para_Op%Cmat(i1+1:i1+f1 , i2+1:i2+f2) + &
                          MatRV%ReVal(1:f1 , 1:f2,iterm_Op)*Val_BasisRot
            ELSE
              para_Op%Rmat(i1+1:i1+f1 , i2+1:i2+f2) =                   &
                                para_Op%Rmat(i1+1:i1+f1 , i2+1:i2+f2) + &
                          MatRV%ReVal(1:f1 , 1:f2,iterm_Op)*Val_BasisRot
            END IF
          END DO
          END DO
          !---- END LOOP on the rotational basis function

        END DO
        END DO
        !================================================================
        ! END loop on i1_h, i2_h
        !================================================================
      END DO

      nreste = nreste - nplus
      IF (nreste .GT. 0) GOTO 98
      !- END of loop of kmem block ----------------------------
      !--------------------------------------------------------
      CALL flush_perso(out_unitp)

      !- determination of Hmax --------------------------------
      DO i=1,para_Op%nb_tot
        IF (para_Op%cplx) THEN
          Hinter = real(para_Op%Cmat(i,i))
        ELSE
          Hinter = para_Op%Rmat(i,i)
        END IF
        IF (Hinter > para_Op%Hmax) para_Op%Hmax = Hinter
      END DO
      !--------------------------------------------------------

      !--------------------------------------------------------
      CALL dealloc_d0MatOp(MatRV)

      IF (allocated(mat2)) THEN
        CALL dealloc_NParray(mat2,'mat2',name_sub)
      END IF
      IF (allocated(mat3)) THEN
        CALL dealloc_NParray(mat3,'mat3',name_sub)
      END IF

      IF (allocated(VecVib)) THEN
        CALL dealloc_NParray(VecVib,'VecVib',name_sub)
      END IF

      IF (allocated(EneVib)) THEN
        CALL dealloc_NParray(EneVib,'EneVib',name_sub)
      END IF

      CALL dealloc_NParray(VecQ,     'VecQ',     name_sub)
      CALL dealloc_NParray(td0b,     'td0b',     name_sub)

      DO i=1,nb_ba
      DO k=1,kmem
        CALL dealloc_d0MatOp(d0MatOpd0bWrho(k,i))
      END DO
      END DO
      deallocate(d0MatOpd0bWrho,stat=err_mem)
      memory = kmem*nb_ba
      CALL error_memo_allo(err_mem,memory,'d0MatOpd0bWrho',name_sub,'d0MatOp')

      CALL dealloc_d0MatOp(d0MatOp)
      !--------------------------------------------------------

      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) para_Op%name_Op,' non symetrized'
        IF (para_Op%cplx) THEN
          CALL Write_Mat(para_Op%Cmat,out_unitp,3)
        ELSE
          CALL Write_Mat(para_Op%Rmat,out_unitp,5)
        END IF
        write(out_unitp,*)
      END IF


!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*)
         write(out_unitp,*) 'END ',name_sub
       END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_MatOpVibRot_WITH_FileGrid_type0


      SUBROUTINE sub_MatOpVibRot_WITH_FileGrid_type0_old(para_Op)
      USE mod_system
      USE mod_nDindex
      USE mod_Constant
      USE mod_SetOp
      IMPLICIT NONE

!----- Operator variables --------------------------------------------
      TYPE (param_Op)  :: para_Op

!------ quadrature points and weight -----------------------------
      real (kind=Rkind) :: WnD

      real (kind=Rkind) :: Qdyn(para_Op%mole%nb_var)
      real (kind=Rkind) :: Qact(para_Op%mole%nb_act1)

      real (kind=Rkind), allocatable :: mat2(:,:)
      real (kind=Rkind), allocatable :: mat3(:,:)
      real (kind=Rkind), allocatable :: VecQ(:)
      real (kind=Rkind), allocatable :: VecVib(:,:),EneVib(:)


!------ for td0b ...         -------------------------------------
      integer                           :: kmem
      real (kind=Rkind), allocatable    :: td0b(:,:)
      TYPE (param_d0MatOp), allocatable :: d0MatOpd0bWrho(:,:)
      TYPE (param_d0MatOp)              :: MatRV

      integer                           :: type_Op


!----- for the  memory ---------------------------------
      integer                       :: nreste,nplus,nb_blocks

!----- divers ----------------------------------------------------
      integer  :: nb_ba,nb_bie,nb_Op,i_Op

      integer  :: n

      integer  :: i,k,KRot,JRot,ibRot,jbRot
      integer  :: i1,i2,f1,f2
      integer  :: i1_h,i2_h,i_h,ib1
      integer  :: J1,J2,iterm,iterm_BasisRot,iterm_Op
      real (kind=Rkind) :: Val_BasisRot


      real (kind=Rkind) :: Hinter,auTOcm_inv

      TYPE (param_d0MatOp) :: d0MatOp

      integer :: nDGridI
      logical :: spectral = .TRUE.
      !logical :: spectral = .FALSE.

!----- for debuging --------------------------------------------------
      integer :: err,err_mem,memory
      character (len=*), parameter :: name_sub="sub_MatOpVibRot_WITH_FileGrid_type0_old"
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (.NOT. para_Op%alloc_mat)                                      &
             CALL alloc_para_Op(para_Op,Mat=.TRUE.,Grid=.FALSE.)
      IF (para_Op%para_ReadOp%para_FileGrid%Type_FileGrid /= 0 .OR.                 &
          para_Op%para_ReadOp%para_FileGrid%Save_MemGrid_done .OR. para_Op%mat_done) RETURN

      para_Op%mat_done = .TRUE.
      n = para_Op%nb_tot
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub,para_Op%n_Op
        write(out_unitp,*) 'nb_act1,nb_var',para_Op%mole%nb_act1,               &
                                   para_Op%mole%nb_var
        write(out_unitp,*) 'nb_tot',para_Op%nb_tot
        write(out_unitp,*)
        CALL write_param_Op(para_Op)
      END IF
!-----------------------------------------------------------

      JRot = para_Op%Para_Tnum%JJ
      spectral = spectral .AND. (JRot > 0)
      IF (JRot > 0) THEN
        CALL init_RotBasis_Param(para_Op%BasisnD%RotBasis,Jrot)
        CALL Write_RotBasis_Param(para_Op%BasisnD%RotBasis)
      END IF

!     ----------------------------------------------------------------
!     - test on nb_bi : nb_bi>0
      IF (para_Op%nb_bi <= 0) THEN
        write(out_unitp,*) ' ERROR : nb_bi MUST be > 0'
        write(out_unitp,*) 'nb_bi',para_Op%nb_bi
        write(out_unitp,*) 'STOP in ',name_sub
        write(out_unitp,*)
        CALL write_param_Op(para_Op)
        STOP
      END IF
!     ----------------------------------------------------------------

      nb_ba  = para_Op%nb_ba
      nb_bie = para_Op%nb_bie

      IF (para_Op%name_Op == 'H') THEN
        type_Op = para_Op%para_PES%Type_HamilOp ! H
        IF (type_Op /= 1) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '    Type_HamilOp MUST be equal to 1 for HADA or cHAC'
          write(out_unitp,*) '    CHECK your data!!'
          STOP
        END IF
      ELSE
        type_Op = 0
      END IF

!-------------------------------------------------------------
!-     memories allocation: td0b, Opd0bWrho matRV, mat2, mat3
!-------------------------------------------------------------
      CALL Init_d0MatOp(MatRV,type_Op,0,nb_ba,JRot=JRot,cplx=para_Op%cplx)

      IF (para_Op%ComOp%contrac_ba_ON_HAC) THEN
        CALL alloc_NParray(mat2,(/nb_ba,nb_ba/),'mat2',name_sub)
        CALL alloc_NParray(mat3,(/nb_ba,nb_ba/),'mat3',name_sub)
      END IF
      CALL alloc_NParray(VecVib,(/nb_ba,nb_ba/),'VecVib',name_sub)
      CALL alloc_NParray(EneVib,(/nb_ba/),'EneVib',name_sub)


      ! selected the optimal value of kmem, as function of max_mem and mem_tot
      kmem = para_Op%nb_qa
      allocate(d0MatOpd0bWrho(para_Op%nb_qa,nb_ba),stat=err)
      memory = kmem*nb_ba
      CALL error_memo_allo(err,memory,'d0MatOpd0bWrho',name_sub,'d0MatOp')
      DO i=1,nb_ba
      DO k=1,kmem
        CALL Init_d0MatOp(d0MatOpd0bWrho(k,i),type_Op,0,nb_bie,         &
                                            JRot=JRot,cplx=para_Op%cplx)
      END DO
      END DO

      nb_blocks = para_Op%nb_qa/kmem
      IF (mod(para_Op%nb_qa,kmem) /= 0) nb_blocks = nb_blocks + 1
      IF (print_level>0) write(out_unitp,*) 'number of blocks',nb_blocks

      CALL alloc_NParray(td0b,(/nb_ba,kmem/),'td0b',name_sub)

      CALL Init_d0MatOp(d0MatOp,type_Op,para_Op%mole%nb_act1,nb_bie,    &
                                             JRot=JRot,cplx=para_Op%cplx)
      !--------------------------------------------------------


      !-- built the Operator matrix ---------------------------
      !--------------------------------------------------------
      IF (para_Op%cplx) THEN
        para_Op%Cmat(:,:) = CZERO
      ELSE
        para_Op%Rmat(:,:) = ZERO
      END IF


      !- Hmin: Vmin and Hmax: the largest diagonal element of H -----
      para_Op%Hmin =  huge(ONE)
      para_Op%Hmax = -huge(ONE)


      !--------------------------------------------------------
      !- loop of kmem block -----------------------------------
      nreste = para_Op%nb_qa


      nDGridI = 0

 98   CONTINUE
      IF (print_level>-1) write(out_unitp,'(a)',ADVANCE='no') 'Block of ReOpd0b(:,:) (%): '
      CALL flush_perso(out_unitp)

      nplus = min(kmem,nreste)
      DO k=1,nplus

        nDGridI = nDGridI + 1
        IF (mod(k,max(1,int(nplus/10))) == 0 .AND. print_level>-1) THEN
          write(out_unitp,'(a,i3)',ADVANCE='no') ' -',                  &
                   int(real(k,kind=Rkind)*HUNDRED/real(nplus,kind=Rkind))
          CALL flush_perso(out_unitp)
        END IF

        CALL sub_reading_Op(nDGridI,para_Op%nb_qa,                  &
                                d0MatOp,para_Op%n_Op,                   &
                                Qdyn,para_Op%mole%nb_var,Qact,          &
                                WnD,para_Op%ComOp)

        iterm = d0MatOp%derive_term_TO_iterm(0,0)
        DO i1_h=1,para_Op%nb_bie
          para_Op%Hmin = min(para_Op%Hmin,d0MatOp%ReVal(i1_h,i1_h,iterm))
          para_Op%Hmax = max(para_Op%Hmax,d0MatOp%ReVal(i1_h,i1_h,iterm))
        END DO

        CALL calc_td0b_OpRVd0bW(nDGridI,k,td0b,d0MatOpd0bWrho,          &
                                WnD,kmem,d0MatOp,para_Op,               &
                                para_Op%BasisnD)

      END DO


      IF (print_level>-1) write(out_unitp,'(a)',ADVANCE='yes') ' - End'
      IF (print_level>-1) write(out_unitp,'(a)',ADVANCE='no') 'Block of ReMatOp(:,:) (%): '
      CALL flush_perso(out_unitp)


      DO i2_h=1,para_Op%nb_bie
        IF (mod(i2_h,max(1,int(para_Op%nb_bie/10))) == 0 .AND.          &
                          print_level>-1 .AND. para_Op%nb_bie /= 0) THEN
            write(out_unitp,'(a,i3)',ADVANCE='no') ' -',                &
              int(real(i2_h,kind=Rkind)*HUNDRED/real(para_Op%nb_bie,kind=Rkind))
            CALL flush_perso(out_unitp)
        END IF
      DO i1_h=1,para_Op%nb_bie
        DO i_Op=1,d0MatOpd0bWrho(1,1)%nb_term


          !$OMP parallel default(none)                               &
          !$OMP shared(para_Op,nplus,i1_h,i2_h,i_Op,print_level)     &
          !$OMP shared(MatRV,td0b,d0MatOpd0bWrho,out_unitp,kmem)     &
          !$OMP private(ib1,k,VecQ)                                  &
          !$OMP num_threads(MatOp_maxth)
          CALL alloc_NParray(VecQ,(/kmem/),'VecQ',name_sub)
          !$OMP do
          DO ib1=1,para_Op%nb_ba
            IF (mod(ib1,max(1,int(para_Op%nb_ba/10))) == 0 .AND.        &
                          print_level>-1 .AND. para_Op%nb_bie == 0) THEN
              write(out_unitp,'(a,i3)',ADVANCE='no') ' -',              &
                        int(real(ib1,kind=Rkind)*HUNDRED/real(para_Op%nb_ba,kind=Rkind))
              CALL flush_perso(out_unitp)
            END IF

            DO k=1,nplus
              VecQ(k) = d0MatOpd0bWrho(k,ib1)%ReVal(i1_h,i2_h,i_Op)
            END DO
            MatRV%ReVal(:,ib1,i_Op) = matmul(td0b(:,:),VecQ(1:nplus))
          END DO
          !$OMP end do
          CALL dealloc_NParray(VecQ,'VecQ',name_sub)
          !$OMP end parallel

          IF (para_Op%ComOp%contrac_ba_ON_HAC) THEN
            mat2 = matmul(MatRV%ReVal(:,:,i_Op),para_Op%ComOp%d0Cba_ON_HAC(:,:,i2_h) )
            mat3 = transpose(para_Op%ComOp%d0Cba_ON_HAC(:,:,i1_h))
            MatRV%ReVal(:,:,i_Op) = matmul(mat3,mat2)
          END IF
        END DO

        IF (para_Op%cplx) THEN
          !$OMP parallel default(none)                             &
          !$OMP shared(para_Op,nplus,i1_h,i2_h,i_Op,print_level)   &
          !$OMP shared(MatRV,td0b,d0MatOpd0bWrho,out_unitp,kmem)   &
          !$OMP private(ib1,k,VecQ)                                &
          !$OMP num_threads(MatOp_maxth)
          CALL alloc_NParray(VecQ,(/kmem/),'VecQ',name_sub)
          !$OMP do
          DO ib1=1,para_Op%nb_ba
            IF (mod(ib1,max(1,int(para_Op%nb_ba/10))) == 0 .AND.        &
                          print_level>-1 .AND. para_Op%nb_bie == 0) THEN
              write(out_unitp,'(a,i3)',ADVANCE='no') ' -',              &
                        int(real(ib1,kind=Rkind)*HUNDRED/real(para_Op%nb_ba,kind=Rkind))
              CALL flush_perso(out_unitp)
            END IF
            DO k=1,nplus
              VecQ(k) = d0MatOpd0bWrho(k,ib1)%ImVal(i1_h,i2_h)
            END DO
            MatRV%ImVal(:,ib1) = matmul(td0b(:,1:nplus),VecQ(1:nplus))
          END DO
          !$OMP end do
          CALL dealloc_NParray(VecQ,'VecQ',name_sub)
          !$OMP end parallel

          IF (para_Op%ComOp%contrac_ba_ON_HAC) THEN
            mat2 = matmul(MatRV%ImVal(:,:),para_Op%ComOp%d0Cba_ON_HAC(:,:,i2_h) )
            mat3 = transpose(para_Op%ComOp%d0Cba_ON_HAC(:,:,i1_h))
            MatRV%ImVal(:,:) = matmul(mat3,mat2)
          END IF
        END IF

        IF (JRot > 0 .AND. para_Op%name_Op == 'H') THEN
          i_Op = MatRV%derive_term_TO_iterm(0,0)
          CALL diagonalization(MatRV%ReVal(:,:,i_Op),EneVib,VecVib,nb_ba,2,1,.TRUE.)
          para_Op%ComOp%ZPE     = Get_ZPE(EneVib)
          para_Op%ComOp%Set_ZPE = .TRUE.
          auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')
          write(out_unitp,*) 'ZPE (cm-1)',para_Op%ComOp%ZPE*auTOcm_inv
          write(out_unitp,*) 'Ene (J=0)',(EneVib(1:min(10,nb_ba))-para_Op%ComOp%ZPE)*auTOcm_inv
        END IF

        IF (spectral) THEN
          DO i_Op=1,MatRV%nb_term
            MatRV%ReVal(:,:,i_Op) = matmul(MatRV%ReVal(:,:,i_Op),VecVib)
            MatRV%ReVal(:,:,i_Op) = matmul(transpose(VecVib),MatRV%ReVal(:,:,i_Op))
          END DO
          IF (para_Op%cplx) THEN
            MatRV%ImVal(:,:) = matmul(MatRV%ImVal(:,:),VecVib)
            MatRV%ImVal(:,:) = matmul(transpose(VecVib),MatRV%ImVal(:,:))
          END IF
        END IF


        IF (debug) THEN
          DO i_Op=1,MatRV%nb_term
            write(out_unitp,*) 'matRV,i_Op',i_Op,                       &
                                'der_term',MatRV%derive_termQact(:,i_Op)
            CALL Write_Mat(MatRV%ReVal(:,:,i_Op),out_unitp,5)
          END DO
        END IF


        !---- LOOP on the rotational basis function
        ! Vibrational contribution
        iterm_Op       = MatRV%derive_term_TO_iterm(0,0)
        DO ibRot=1,para_Op%nb_bRot

          i1 = (ibRot-1)*para_Op%nb_baie + sum(para_Op%ComOp%nb_ba_ON_HAC(1:i1_h-1))
          i2 = (ibRot-1)*para_Op%nb_baie + sum(para_Op%ComOp%nb_ba_ON_HAC(1:i2_h-1))
          f1 = para_Op%ComOp%nb_ba_ON_HAC(i1_h)
          f2 = para_Op%ComOp%nb_ba_ON_HAC(i2_h)

          IF (para_Op%cplx) THEN
            para_Op%Cmat(i1+1:i1+f1 , i2+1:i2+f2) =                     &
                                para_Op%Cmat(i1+1:i1+f1 , i2+1:i2+f2) + &
                              cmplx(MatRV%ReVal(1:f1 , 1:f2,iterm_Op),  &
                                    MatRV%ImVal(1:f1 , 1:f2),kind=Rkind)
          ELSE
            para_Op%Rmat(i1+1:i1+f1 , i2+1:i2+f2) =                     &
                                para_Op%Rmat(i1+1:i1+f1 , i2+1:i2+f2) + &
                                MatRV%ReVal(1:f1 , 1:f2,iterm_Op)
          END IF
        END DO

        !---- END LOOP on the rotational basis function
        IF (MatRV%JRot > 0 .AND. para_Op%name_Op == 'H') THEN
        !---- LOOP on the rotational basis function
        ! Rotaional contribution
        DO J1=-3,-1
        DO J2=-3,-1
          iterm_Op       = MatRV%derive_term_TO_iterm(J1,J2)
          iterm_BasisRot = para_Op%BasisnD%RotBasis%tab_der_TO_iterm(J1,J2)
          !write(6,*) 'J1,J2',J1,J2,'iterm_Op,iterm_BasisRot',iterm_Op,iterm_BasisRot

          DO ibRot=1,para_Op%nb_bRot
          DO jbRot=1,para_Op%nb_bRot
            Val_BasisRot = para_Op%BasisnD%RotBasis%tab_RotOp(ibRot,jbRot,iterm_BasisRot)
            IF (abs(Val_BasisRot) < ONETENTH**9) CYCLE

            i1 = (ibRot-1)*para_Op%nb_baie + sum(para_Op%ComOp%nb_ba_ON_HAC(1:i1_h-1))
            i2 = (jbRot-1)*para_Op%nb_baie + sum(para_Op%ComOp%nb_ba_ON_HAC(1:i2_h-1))
            f1 = para_Op%ComOp%nb_ba_ON_HAC(i1_h)
            f2 = para_Op%ComOp%nb_ba_ON_HAC(i2_h)

            IF (debug) THEN
              write(6,*) 'J1,J2',J1,J2
              write(6,*) 'i1_h,ibRot,i1+1:i1+f1',i1_h,ibRot,i1+1,i1+f1
              write(6,*) 'i2_h,jbRot,i2+1:i2+f2',i2_h,jbRot,i2+1,i2+f2
            END IF

            IF (para_Op%cplx) THEN
              para_Op%Cmat(i1+1:i1+f1 , i2+1:i2+f2) =                   &
                                para_Op%Cmat(i1+1:i1+f1 , i2+1:i2+f2) + &
                          MatRV%ReVal(1:f1 , 1:f2,iterm_Op)*Val_BasisRot
            ELSE
              para_Op%Rmat(i1+1:i1+f1 , i2+1:i2+f2) =                   &
                                para_Op%Rmat(i1+1:i1+f1 , i2+1:i2+f2) + &
                          MatRV%ReVal(1:f1 , 1:f2,iterm_Op)*Val_BasisRot
            END IF
          END DO
          END DO
        END DO
        END DO
        !---- END LOOP on the rotational basis function

        !---- LOOP on the rotational basis function
        ! Coriolis contribution
        DO J1=-3,-1
          iterm_Op       = MatRV%derive_term_TO_iterm(J1,0)
          iterm_BasisRot = para_Op%BasisnD%RotBasis%tab_der_TO_iterm(J1,0)
          DO ibRot=1,para_Op%nb_bRot
          DO jbRot=1,para_Op%nb_bRot
            Val_BasisRot = -para_Op%BasisnD%RotBasis%tab_RotOp(ibRot,jbRot,iterm_BasisRot)
            ! the minus sign is comming from i * i (i=sqrt(-1))

            IF (abs(Val_BasisRot) < ONETENTH**9) CYCLE

            i1 = (ibRot-1)*para_Op%nb_baie + sum(para_Op%ComOp%nb_ba_ON_HAC(1:i1_h-1))
            i2 = (jbRot-1)*para_Op%nb_baie + sum(para_Op%ComOp%nb_ba_ON_HAC(1:i2_h-1))
            f1 = para_Op%ComOp%nb_ba_ON_HAC(i1_h)
            f2 = para_Op%ComOp%nb_ba_ON_HAC(i2_h)

            IF (debug) THEN
              write(6,*) 'J1',J1
              write(6,*) 'i1_h,ibRot,i1+1:i1+f1',i1_h,ibRot,i1+1,i1+f1
              write(6,*) 'i2_h,jbRot,i2+1:i2+f2',i2_h,jbRot,i2+1,i2+f2
            END IF

            IF (para_Op%cplx) THEN
              para_Op%Cmat(i1+1:i1+f1 , i2+1:i2+f2) =                   &
                  para_Op%Cmat(i1+1:i1+f1 , i2+1:i2+f2) +               &
                              MatRV%ReVal(1:f1 , 1:f2,iterm_Op)*Val_BasisRot
            ELSE
              para_Op%Rmat(i1+1:i1+f1 , i2+1:i2+f2) =                   &
                  para_Op%Rmat(i1+1:i1+f1 , i2+1:i2+f2) +               &
                              MatRV%ReVal(1:f1 , 1:f2,iterm_Op)*Val_BasisRot
            END IF
          END DO
          END DO
        END DO
        !---- END LOOP on the rotational basis function
        END IF

      END DO
      END DO



      IF (print_level>-1) write(out_unitp,'(a)',ADVANCE='yes') ' - End'
      CALL flush_perso(out_unitp)

      nreste = nreste - nplus
      IF (nreste .GT. 0) GOTO 98
      !- END of loop of kmem block ----------------------------
      !--------------------------------------------------------
      CALL flush_perso(out_unitp)

      !- determination of Hmax --------------------------------
      DO i=1,para_Op%nb_tot
        IF (para_Op%cplx) THEN
          Hinter = real(para_Op%Cmat(i,i))
        ELSE
          Hinter = para_Op%Rmat(i,i)
        END IF
        IF (Hinter > para_Op%Hmax) para_Op%Hmax = Hinter
      END DO
      !--------------------------------------------------------

      !--------------------------------------------------------
      CALL dealloc_d0MatOp(MatRV)

      IF (allocated(mat2)) THEN
        CALL dealloc_NParray(mat2,'mat2',name_sub)
      END IF
      IF (allocated(mat3)) THEN
        CALL dealloc_NParray(mat3,'mat3',name_sub)
      END IF

      IF (allocated(VecVib)) THEN
        CALL dealloc_NParray(VecVib,'VecVib',name_sub)
      END IF

      IF (allocated(EneVib)) THEN
        CALL dealloc_NParray(EneVib,'EneVib',name_sub)
      END IF

      CALL dealloc_NParray(td0b,     'td0b',     name_sub)

      DO i=1,nb_ba
      DO k=1,kmem
        CALL dealloc_d0MatOp(d0MatOpd0bWrho(k,i))
      END DO
      END DO
      deallocate(d0MatOpd0bWrho,stat=err)
      memory = kmem*nb_ba
      CALL error_memo_allo(err,memory,'d0MatOpd0bWrho',name_sub,'d0MatOp')

      CALL dealloc_d0MatOp(d0MatOp)
      !--------------------------------------------------------

      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) para_Op%name_Op,' non symetrized'
        IF (para_Op%cplx) THEN
          CALL Write_Mat(para_Op%Cmat,out_unitp,3)
        ELSE
          CALL Write_Mat(para_Op%Rmat,out_unitp,5)
        END IF
        write(out_unitp,*)
      END IF


!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*)
         write(out_unitp,*) 'END ',name_sub
       END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_MatOpVibRot_WITH_FileGrid_type0_old

      SUBROUTINE sub_MatOpVibRot_WITH_FileGrid_type0_v1(para_Op)
      USE mod_system
      USE mod_nDindex
      USE mod_SetOp
      IMPLICIT NONE

!----- Operator variables --------------------------------------------
      TYPE (param_Op)  :: para_Op

!------ quadrature points and weight -----------------------------
      real (kind=Rkind) :: WnD

      real (kind=Rkind) :: Qdyn(para_Op%mole%nb_var)
      real (kind=Rkind) :: Qact(para_Op%mole%nb_act1)

      real (kind=Rkind), allocatable :: mat2(:,:)
      real (kind=Rkind), allocatable :: mat3(:,:)
      real (kind=Rkind), allocatable :: VecQ(:)
      real (kind=Rkind), allocatable :: VecVib(:,:),EneVib(:)


!------ for td0b ...         -------------------------------------
      integer                           :: kmem
      real (kind=Rkind), allocatable    :: td0b(:,:)
      TYPE (param_d0MatOp), allocatable :: d0MatOpd0bWrho(:,:)
      TYPE (param_d0MatOp)              :: MatRV

      integer                           :: type_Op


!----- for the  memory ---------------------------------
      integer                       :: nreste,nplus,nb_blocks

!----- divers ----------------------------------------------------
      integer  :: nb_ba,nb_bie,nb_Op,i_Op

      integer  :: n

      integer  :: i,k,KRot,JRot,ibRot,jbRot
      integer  :: i1,i2,f1,f2
      integer  :: i1_h,i2_h,i_h,ib1
      integer  :: J1,J2,iterm,iterm_BasisRot,iterm_Op,nb_term_BasisRot
      real (kind=Rkind) :: Val_BasisRot


      real (kind=Rkind) :: Hinter,auTOcm_inv

      TYPE (param_d0MatOp) :: d0MatOp

      integer :: nDGridI
      logical :: spectral = .TRUE.
      !logical :: spectral = .FALSE.

!----- for debuging --------------------------------------------------
      integer :: err,err_mem,memory
      character (len=*), parameter :: name_sub="sub_MatOpVibRot_WITH_FileGrid_type0_v1"
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (.NOT. para_Op%alloc_mat)                                      &
             CALL alloc_para_Op(para_Op,Mat=.TRUE.,Grid=.FALSE.)
      IF (para_Op%para_ReadOp%para_FileGrid%Type_FileGrid /= 0 .OR.                 &
          para_Op%para_ReadOp%para_FileGrid%Save_MemGrid_done .OR. para_Op%mat_done) RETURN

      para_Op%mat_done = .TRUE.
      n = para_Op%nb_tot
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub,para_Op%n_Op
        write(out_unitp,*) 'nb_act1,nb_var',para_Op%mole%nb_act1,               &
                                   para_Op%mole%nb_var
        write(out_unitp,*) 'nb_tot',para_Op%nb_tot
        write(out_unitp,*)
        !CALL write_param_Op(para_Op)
      END IF
!-----------------------------------------------------------

      JRot = para_Op%Para_Tnum%JJ
      spectral = spectral .AND. (JRot > 0)
      CALL init_RotBasis_Param(para_Op%BasisnD%RotBasis,Jrot)

!     ----------------------------------------------------------------
!     - test on nb_bi : nb_bi>0
      IF (para_Op%nb_bi <= 0) THEN
        write(out_unitp,*) ' ERROR : nb_bi MUST be > 0'
        write(out_unitp,*) 'nb_bi',para_Op%nb_bi
        write(out_unitp,*) 'STOP in ',name_sub
        write(out_unitp,*)
        CALL write_param_Op(para_Op)
        STOP
      END IF
!     ----------------------------------------------------------------

      nb_ba  = para_Op%nb_ba
      nb_bie = para_Op%nb_bie

      IF (para_Op%name_Op == 'H') THEN
        type_Op = 1
      ELSE
        type_Op = 0
      END IF

!-------------------------------------------------------------
!-     memories allocation: td0b, Opd0bWrho matRV, mat2, mat3
!-------------------------------------------------------------
      CALL Init_d0MatOp(MatRV,type_Op,0,nb_ba,JRot=JRot,cplx=para_Op%cplx)

      IF (para_Op%ComOp%contrac_ba_ON_HAC) THEN
        CALL alloc_NParray(mat2,(/nb_ba,nb_ba/),'mat2',name_sub)
        CALL alloc_NParray(mat3,(/nb_ba,nb_ba/),'mat3',name_sub)
      END IF
      CALL alloc_NParray(VecVib,(/nb_ba,nb_ba/),'VecVib',name_sub)
      CALL alloc_NParray(EneVib,(/nb_ba/),'EneVib',name_sub)


      ! selected the optimal value of kmem, as function of max_mem and mem_tot
      kmem = para_Op%nb_qa
      allocate(d0MatOpd0bWrho(para_Op%nb_qa,nb_ba),stat=err)
      memory = kmem*nb_ba
      CALL error_memo_allo(err,memory,'d0MatOpd0bWrho',name_sub,'d0MatOp')
      DO i=1,nb_ba
      DO k=1,kmem
        CALL Init_d0MatOp(d0MatOpd0bWrho(k,i),type_Op,0,nb_bie,         &
                                            JRot=JRot,cplx=para_Op%cplx)
      END DO
      END DO

      nb_blocks = para_Op%nb_qa/kmem
      IF (mod(para_Op%nb_qa,kmem) /= 0) nb_blocks = nb_blocks + 1
      IF (print_level>0) write(out_unitp,*) 'number of blocks',nb_blocks

      CALL alloc_NParray(td0b,(/nb_ba,kmem/),'td0b',name_sub)

      CALL Init_d0MatOp(d0MatOp,type_Op,para_Op%mole%nb_act1,nb_bie,    &
                                             JRot=JRot,cplx=para_Op%cplx)
      !--------------------------------------------------------


      !-- built the Operator matrix ---------------------------
      !--------------------------------------------------------
      IF (para_Op%cplx) THEN
        para_Op%Cmat(:,:) = CZERO
      ELSE
        para_Op%Rmat(:,:) = ZERO
      END IF


      !- Hmin: Vmin and Hmax: the largest diagonal element of H -----
      para_Op%Hmin =  huge(ONE)
      para_Op%Hmax = -huge(ONE)


      !--------------------------------------------------------
      !- loop of kmem block -----------------------------------
      nreste = para_Op%nb_qa


      nDGridI = 0

 98   CONTINUE
      IF (print_level>-1) write(out_unitp,'(a)',ADVANCE='no') 'Block of ReOpd0b(:,:) (%): '
      CALL flush_perso(out_unitp)

      nplus = min(kmem,nreste)
      DO k=1,nplus

        nDGridI = nDGridI + 1
        IF (mod(k,max(1,int(nplus/10))) == 0 .AND. print_level>-1) THEN
          write(out_unitp,'(a,i3)',ADVANCE='no') ' -',                  &
                   int(real(k,kind=Rkind)*HUNDRED/real(nplus,kind=Rkind))
          CALL flush_perso(out_unitp)
        END IF

        CALL sub_reading_Op(nDGridI,para_Op%nb_qa,                  &
                                d0MatOp,para_Op%n_Op,                   &
                                Qdyn,para_Op%mole%nb_var,Qact,          &
                                WnD,para_Op%ComOp)

        iterm = d0MatOp%derive_term_TO_iterm(0,0)
        DO i1_h=1,para_Op%nb_bie
          para_Op%Hmin = min(para_Op%Hmin,d0MatOp%ReVal(i1_h,i1_h,iterm))
          para_Op%Hmax = max(para_Op%Hmax,d0MatOp%ReVal(i1_h,i1_h,iterm))
        END DO

        CALL calc_td0b_OpRVd0bW(nDGridI,k,td0b,d0MatOpd0bWrho,          &
                                WnD,kmem,d0MatOp,para_Op,               &
                                para_Op%BasisnD)

      END DO


      IF (print_level>-1) write(out_unitp,'(a)',ADVANCE='yes') ' - End'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      IF (print_level>-1) write(out_unitp,'(a)',ADVANCE='no') 'Block of ReMatOp(:,:) (%): '
      CALL flush_perso(out_unitp)


      DO i2_h=1,para_Op%nb_bie
        IF (mod(i2_h,max(1,int(para_Op%nb_bie/10))) == 0 .AND.          &
                          print_level>-1 .AND. para_Op%nb_bie /= 0) THEN
            write(out_unitp,'(a,i3)',ADVANCE='no') ' -',                &
              int(real(i2_h,kind=Rkind)*HUNDRED/real(para_Op%nb_bie,kind=Rkind))
            CALL flush_perso(out_unitp)
        END IF
      DO i1_h=1,para_Op%nb_bie
        DO i_Op=1,d0MatOpd0bWrho(1,1)%nb_term


          !$OMP parallel default(none)                                  &
          !$OMP shared(para_Op,nplus,i1_h,i2_h,i_Op,print_level)        &
          !$OMP shared(MatRV,td0b,d0MatOpd0bWrho,out_unitp,kmem)        &
          !$OMP private(ib1,k,VecQ)                                     &
          !$OMP num_threads(MatOp_maxth)
          CALL alloc_NParray(VecQ,(/kmem/),'VecQ',name_sub)
          !$OMP do
          DO ib1=1,para_Op%nb_ba
            IF (mod(ib1,max(1,int(para_Op%nb_ba/10))) == 0 .AND.        &
                          print_level>-1 .AND. para_Op%nb_bie == 0) THEN
              write(out_unitp,'(a,i3)',ADVANCE='no') ' -',              &
                   int(real(ib1,kind=Rkind)*HUNDRED/real(para_Op%nb_ba,kind=Rkind))
              CALL flush_perso(out_unitp)
            END IF

            DO k=1,nplus
              VecQ(k) = d0MatOpd0bWrho(k,ib1)%ReVal(i1_h,i2_h,i_Op)
            END DO
            MatRV%ReVal(:,ib1,i_Op) = matmul(td0b(:,:),VecQ(1:nplus))
          END DO
          !$OMP end do
          CALL dealloc_NParray(VecQ,'VecQ',name_sub)
          !$OMP end parallel

          IF (para_Op%ComOp%contrac_ba_ON_HAC) THEN
            mat2 = matmul(MatRV%ReVal(:,:,i_Op),para_Op%ComOp%d0Cba_ON_HAC(:,:,i2_h) )
            mat3 = transpose(para_Op%ComOp%d0Cba_ON_HAC(:,:,i1_h))
            MatRV%ReVal(:,:,i_Op) = matmul(mat3,mat2)
          END IF
        END DO

        IF (para_Op%cplx) THEN
          !$OMP parallel default(none)                             &
          !$OMP shared(para_Op,nplus,i1_h,i2_h,i_Op,print_level)   &
          !$OMP shared(MatRV,td0b,d0MatOpd0bWrho,out_unitp,kmem)   &
          !$OMP private(ib1,k,VecQ)                                &
          !$OMP num_threads(MatOp_maxth)
          CALL alloc_NParray(VecQ,(/kmem/),'VecQ',name_sub)
          !$OMP do
          DO ib1=1,para_Op%nb_ba
            IF (mod(ib1,max(1,int(para_Op%nb_ba/10))) == 0 .AND.        &
                          print_level>-1 .AND. para_Op%nb_bie == 0) THEN
              write(out_unitp,'(a,i3)',ADVANCE='no') ' -',              &
                        int(real(ib1,kind=Rkind)*HUNDRED/real(para_Op%nb_ba,kind=Rkind))
              CALL flush_perso(out_unitp)
            END IF
            DO k=1,nplus
              VecQ(k) = d0MatOpd0bWrho(k,ib1)%ImVal(i1_h,i2_h)
            END DO
            MatRV%ImVal(:,ib1) = matmul(td0b(:,1:nplus),VecQ(1:nplus))
          END DO
          !$OMP end do
          CALL dealloc_NParray(VecQ,'VecQ',name_sub)
          !$OMP end parallel

          IF (para_Op%ComOp%contrac_ba_ON_HAC) THEN
            mat2 = matmul(MatRV%ImVal(:,:),para_Op%ComOp%d0Cba_ON_HAC(:,:,i2_h) )
            mat3 = transpose(para_Op%ComOp%d0Cba_ON_HAC(:,:,i1_h))
            MatRV%ImVal(:,:) = matmul(mat3,mat2)
          END IF
        END IF

!        IF (JRot > 0 .AND. para_Op%name_Op == 'H') THEN
!          i_Op = MatRV%derive_term_TO_iterm(0,0)
!          CALL diagonalization(MatRV%ReVal(:,:,i_Op),EneVib,VecVib,nb_ba,2,1,.TRUE.)
!          para_Op%ComOp%ZPE     = Get_ZPE(EneVib)
!          para_Op%ComOp%Set_ZPE = .TRUE.
!          auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')
!          write(out_unitp,*) 'ZPE',para_Op%ComOp%ZPE*auTOcm_inv
!          write(out_unitp,*) 'Ene (J=0)',(EneVib(1:min(10,nb_ba))-para_Op%ComOp%ZPE)*auTOcm_inv
!        END IF
!
!        IF (spectral) THEN
!          DO i_Op=1,MatRV%nb_term
!            MatRV%ReVal(:,:,i_Op) = matmul(MatRV%ReVal(:,:,i_Op),VecVib)
!            MatRV%ReVal(:,:,i_Op) = matmul(transpose(VecVib),MatRV%ReVal(:,:,i_Op))
!          END DO
!          IF (para_Op%cplx) THEN
!            MatRV%ImVal(:,:) = matmul(MatRV%ImVal(:,:),VecVib)
!            MatRV%ImVal(:,:) = matmul(transpose(VecVib),MatRV%ImVal(:,:))
!          END IF
!        END IF


        IF (debug) THEN
          DO i_Op=1,MatRV%nb_term
            write(out_unitp,*) 'matRV,i_Op',i_Op,                       &
                                'der_term',MatRV%derive_termQact(:,i_Op)
            CALL Write_Mat(MatRV%ReVal(:,:,i_Op),out_unitp,5)
          END DO
        END IF

        !---- LOOP on the rotational basis function
        ! Rotational contribution
        nb_term_BasisRot = 0
        IF (para_Op%name_Op == 'H') nb_term_BasisRot = para_Op%BasisnD%RotBasis%nb_term

        DO iterm_BasisRot=0,nb_term_BasisRot
          J1       = para_Op%BasisnD%RotBasis%tab_iterm_TO_der(1,iterm_BasisRot)
          J2       = para_Op%BasisnD%RotBasis%tab_iterm_TO_der(2,iterm_BasisRot)
          iterm_Op = MatRV%derive_term_TO_iterm(J1,J2)
          !write(6,*) 'J1,J2',J1,J2,'iterm_Op,iterm_BasisRot',iterm_Op,iterm_BasisRot

          DO ibRot=1,para_Op%nb_bRot
          DO jbRot=1,para_Op%nb_bRot
            Val_BasisRot = para_Op%BasisnD%RotBasis%tab_RotOp(ibRot,jbRot,iterm_BasisRot)
            IF (abs(Val_BasisRot) < ONETENTH**9) CYCLE

            i1 = (ibRot-1)*para_Op%nb_baie + sum(para_Op%ComOp%nb_ba_ON_HAC(1:i1_h-1))
            i2 = (jbRot-1)*para_Op%nb_baie + sum(para_Op%ComOp%nb_ba_ON_HAC(1:i2_h-1))
            f1 = para_Op%ComOp%nb_ba_ON_HAC(i1_h)
            f2 = para_Op%ComOp%nb_ba_ON_HAC(i2_h)

            IF (debug) THEN
              write(6,*) 'J1,J2',J1,J2
              write(6,*) 'i1_h,ibRot,i1+1:i1+f1',i1_h,ibRot,i1+1,i1+f1
              write(6,*) 'i2_h,jbRot,i2+1:i2+f2',i2_h,jbRot,i2+1,i2+f2
            END IF

            IF (para_Op%cplx) THEN
              para_Op%Cmat(i1+1:i1+f1 , i2+1:i2+f2) =                   &
                                para_Op%Cmat(i1+1:i1+f1 , i2+1:i2+f2) + &
                          MatRV%ReVal(1:f1 , 1:f2,iterm_Op)*Val_BasisRot
            ELSE
              para_Op%Rmat(i1+1:i1+f1 , i2+1:i2+f2) =                   &
                                para_Op%Rmat(i1+1:i1+f1 , i2+1:i2+f2) + &
                          MatRV%ReVal(1:f1 , 1:f2,iterm_Op)*Val_BasisRot
            END IF
          END DO
          END DO
        END DO
        END DO
        !---- END LOOP on the rotational basis function


      END DO



      IF (print_level>-1) write(out_unitp,'(a)',ADVANCE='yes') ' - End'
      CALL flush_perso(out_unitp)

      nreste = nreste - nplus
      IF (nreste .GT. 0) GOTO 98
      !- END of loop of kmem block ----------------------------
      !--------------------------------------------------------
      CALL flush_perso(out_unitp)

      !- determination of Hmax --------------------------------
      DO i=1,para_Op%nb_tot
        IF (para_Op%cplx) THEN
          Hinter = real(para_Op%Cmat(i,i))
        ELSE
          Hinter = para_Op%Rmat(i,i)
        END IF
        IF (Hinter > para_Op%Hmax) para_Op%Hmax = Hinter
      END DO
      !--------------------------------------------------------

      !--------------------------------------------------------
      CALL dealloc_d0MatOp(MatRV)

      IF (allocated(mat2)) THEN
        CALL dealloc_NParray(mat2,'mat2',name_sub)
      END IF
      IF (allocated(mat3)) THEN
        CALL dealloc_NParray(mat3,'mat3',name_sub)
      END IF

      IF (allocated(VecVib)) THEN
        CALL dealloc_NParray(VecVib,'VecVib',name_sub)
      END IF

      IF (allocated(EneVib)) THEN
        CALL dealloc_NParray(EneVib,'EneVib',name_sub)
      END IF

      CALL dealloc_NParray(td0b,     'td0b',     name_sub)

      DO i=1,nb_ba
      DO k=1,kmem
        CALL dealloc_d0MatOp(d0MatOpd0bWrho(k,i))
      END DO
      END DO
      deallocate(d0MatOpd0bWrho,stat=err)
      memory = kmem*nb_ba
      CALL error_memo_allo(err,memory,'d0MatOpd0bWrho',name_sub,'d0MatOp')

      CALL dealloc_d0MatOp(d0MatOp)
      !--------------------------------------------------------

      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) para_Op%name_Op,' non symetrized'
        IF (para_Op%cplx) THEN
          CALL Write_Mat(para_Op%Cmat,out_unitp,3)
        ELSE
          CALL Write_Mat(para_Op%Rmat,out_unitp,5)
        END IF
        write(out_unitp,*)
      END IF


!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*)
         write(out_unitp,*) 'END ',name_sub
       END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_MatOpVibRot_WITH_FileGrid_type0_v1


!================================================================
!
!     Read Operator (matrix form)
!
!================================================================
      SUBROUTINE sub_Read_MatOp(para_Op)
      USE mod_system
      USE mod_SetOp
      IMPLICIT NONE

!----- Operator variables --------------------------------------------
      TYPE (param_Op)  :: para_Op


!----- divers ----------------------------------------------------

      integer  :: i,j,k,n


      integer           :: nbcol,partial
      logical           :: diag_only,print_Op
      real (kind=Rkind) :: conv,Opij
      integer           :: idum,err
      namelist /read_Op/ nbcol,diag_only,conv,print_Op,partial

!----- for debuging --------------------------------------------------
      character (len=*), parameter ::name_sub='sub_Read_MatOp'
      logical,parameter :: debug=.FALSE.
!     logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------

      IF (.NOT. para_Op%alloc_mat)                                      &
             CALL alloc_para_Op(para_Op,Mat=.TRUE.,Grid=.FALSE.)
      IF (para_Op%mat_done) RETURN

      para_Op%mat_done = .TRUE.
      print_Op = .FALSE.
      n = para_Op%nb_tot
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING sub_Read_MatOp',para_Op%n_Op
        write(out_unitp,*) 'nb_act1,nb_var',para_Op%mole%nb_act1,               &
                                   para_Op%mole%nb_var
        write(out_unitp,*) 'read MatOp',para_Op%n_Op,para_Op%read_Op
        write(out_unitp,*) 'nb_tot',para_Op%nb_tot
        write(out_unitp,*)
        CALL write_param_Op(para_Op)
      END IF
!-----------------------------------------------------------

        nbcol = 5
        diag_only = .FALSE.
        partial   = 0
        conv  = ONE
        read(in_unitp,read_Op)
        IF (print_level>-1) write(out_unitp,read_Op)
        IF (diag_only .AND. partial > 0) THEN
          write(out_unitp,*) ' ERROR in the read_Op namelist of ',name_sub
          write(out_unitp,*) ' you cannot have diag_only=t and partial > 0'
          STOP
        END IF
        IF (para_Op%cplx) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' CANNOT read complex operator yet!!!'
          STOP
        ELSE
          IF (diag_only) THEN
            para_Op%Rmat(:,:) = ZERO
            DO i=1,para_Op%nb_tot
              read(in_unitp,*) idum,para_Op%Rmat(i,i)
            END DO
          ELSE IF (partial > 0) THEN
            para_Op%Rmat(:,:) = ZERO
            DO k=1,partial
              read(in_unitp,*) i,j,Opij
              para_Op%Rmat(i,j) = Opij
              para_Op%Rmat(j,i) = Opij
            END DO

          ELSE
            CALL Read_Mat(para_Op%Rmat,in_unitp,para_Op%nb_tot,err)
            IF (err /= 0) THEN
              write(out_unitp,*) 'ERROR in ',name_sub
              write(out_unitp,*) ' reading the matrix "para_Op%Rmat"'
              STOP
            END IF
          END IF
          IF (conv /= ONE) para_Op%Rmat = para_Op%Rmat * conv
        END IF

        IF (debug .OR. print_Op) THEN
          write(out_unitp,*) 'Read: ',para_Op%name_Op
          IF (para_Op%cplx) THEN
            CALL Write_Mat(para_Op%Cmat,out_unitp,3)
          ELSE
            CALL Write_Mat(para_Op%Rmat,out_unitp,5)
          END IF
        END IF
        IF (debug) THEN
          write(out_unitp,*)
          write(out_unitp,*) 'END ',name_sub
        END IF
        CALL flush_perso(out_unitp)
      END SUBROUTINE sub_Read_MatOp
!===============================================================================
!
!     not(para_Op%para_FileGrid%Type_FileGrid = 0
!      and
!     para_Op%para_FileGrid%Save_MemGrid  = .F.)
!
!     second version for MatOp_omp=2,1
!===============================================================================
      SUBROUTINE sub_MatOp_direct2(para_Op)
      USE mod_system
!$    USE omp_lib, only : OMP_GET_THREAD_NUM

      USE mod_SetOp
      USE mod_ana_psi
      USE mod_MPI
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op), intent(inout) :: para_Op

!------ working parameters --------------------------------
      integer       :: i,nb_thread


!----- for debuging --------------------------------------------------
      character (len=*), parameter ::name_sub='sub_MatOp_direct2'
      logical, parameter :: debug=.FALSE.
!      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
        write(out_unitp,*) 'Build matrix of ',para_Op%nb_tot
        CALL flush_perso(out_unitp)
      END IF

      IF (.NOT. para_Op%alloc_mat)                                      &
                     CALL alloc_para_Op(para_Op,Mat=.TRUE.,Grid=.FALSE.)
      IF (para_Op%mat_done) RETURN

      para_Op%Make_mat = .FALSE.

!     - scaling of Op ---------------------------------------
      para_Op%E0     = ZERO
      para_Op%Esc    = ONE
      para_Op%scaled = .FALSE.
!-----------------------------------------------------------
      IF (MatOp_omp /= 2) THEN
        nb_thread = 1
      ELSE
        nb_thread = MatOp_maxth
      END IF
      IF (print_level>-1) write(out_unitp,*) 'nb_thread in ',name_sub,' : ',nb_thread


!       ----------------------------------------------------------
!       - build H and H0
        IF (print_level > -1) THEN
          write(out_unitp,'(a)')              'MatOp(:,i) (%): [--0-10-20-30-40-50-60-70-80-90-100]'
          write(out_unitp,'(a)',ADVANCE='no') 'MatOp(:,i) (%): ['
          CALL flush_perso(out_unitp)
        END IF
        
        !$OMP parallel do default(none)                    &
        !$OMP shared(para_Op,print_level,out_unitp,MPI_id) &
        !$OMP private(i)                                   &
        !$OMP num_threads(nb_thread)
        DO i=1,para_Op%nb_tot
          !$ !write(out_unitp,*) "thread",omp_get_thread_num(),"doing",i ; CALL flush_perso(out_unitp)
          CALL sub_OpBasisFi(para_Op,i)

          IF (mod(i,max(1,int(para_Op%nb_tot/10))) == 0 .AND. print_level > -1         &
              .AND. MPI_id==0) THEN
            write(out_unitp,'(a)',ADVANCE='no') '---'
            CALL flush_perso(out_unitp)
          END IF

        END DO
        !$OMP end parallel do

        IF (print_level > -1 .AND. MPI_id==0) THEN
          write(out_unitp,'(a)',ADVANCE='yes') '----]'
          CALL flush_perso(out_unitp)
        END IF

      IF (debug) THEN
        write(out_unitp,*) para_Op%name_Op,' non-symmetrized'
        IF (para_Op%cplx) THEN
          CALL Write_Mat(para_Op%Cmat,out_unitp,3)
        ELSE
          CALL Write_Mat(para_Op%Rmat,out_unitp,5)
        END IF
      END IF

!     ----------------------------------------------------------
      para_Op%Make_mat = .TRUE.

!----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'END ',name_sub
       END IF
!----------------------------------------------------------

      END SUBROUTINE sub_MatOp_direct2
      
      SUBROUTINE sub_MatOp_direct1(para_Op)
      USE mod_system
      USE mod_SetOp
      USE mod_OpPsi
      USE mod_psi_set_alloc
      USE mod_psi_Op
      USE mod_MPI
      IMPLICIT NONE


!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op), intent(inout) :: para_Op

!------ for OpenMP --------------------------------
      TYPE (param_psi), allocatable :: psi(:),Hpsi(:)
      integer       :: i,ib,ib0,n

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter ::name_sub='sub_MatOp_direct1'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
        write(out_unitp,*) 'Build matrix of ',para_Op%nb_tot
        CALL flush_perso(out_unitp)
      END IF

      IF (.NOT. para_Op%alloc_mat)                                      &
             CALL alloc_para_Op(para_Op,Mat=.TRUE.,Grid=.FALSE.)
      IF (para_Op%mat_done) RETURN

      para_Op%Make_mat = .FALSE.
      !- scaling of Op ---------------------------------------
      para_Op%E0     = ZERO
      para_Op%Esc    = ONE
      para_Op%scaled = .FALSE.
!-----------------------------------------------------------

!---- initialization -------------------------------------
      n = min(100,para_Op%nb_tot)
      IF(MPI_id==0) THEN

        CALL alloc_NParray(psi, (/n/),"psi", name_sub)
        CALL alloc_NParray(Hpsi,(/n/),"Hpsi",name_sub)

        DO i=1,n
          CALL init_psi(psi(i),para_Op,para_Op%cplx)
          CALL init_psi(Hpsi(i),para_Op,para_Op%cplx)
          psi(i) = ZERO
        END DO

!     ----------------------------------------------------------
!       - build H and H0
        IF (print_level > -1) THEN
          write(out_unitp,'(a)')              'MatOp(:,i) (%): [--0-10-20-30-40-50-60-70-80-90-100]'
          write(out_unitp,'(a)',ADVANCE='no') 'MatOp(:,i) (%): ['
          CALL flush_perso(out_unitp)
        END IF
      ENDIF ! for MPI=0

      ib = 0
      DO
        IF (ib >= para_Op%nb_tot) EXIT

        IF (para_Op%cplx) THEN
          ib0 = ib
          DO i=1,n
            ib = ib + 1
            IF (ib > para_Op%nb_tot) EXIT
            IF(MPI_id==0) THEN
              psi(i)%CvecB(:) = CZERO
              psi(i)%CvecB(ib) = CONE
              CALL Set_symab_OF_psiBasisRep(psi(i))
            ENDIF
          END DO

          CALL sub_TabOpPsi(Psi,HPsi,para_Op)

          ib = ib0
          DO i=1,n
            ib = ib + 1
            IF (ib > para_Op%nb_tot) EXIT

            !> Cmat assigned here
            IF(MPI_id==0) para_Op%Cmat(:,ib)  = Hpsi(i)%CvecB(:)

            IF (mod(ib,max(1,int(para_Op%nb_tot/10))) == 0 .AND. print_level > -1      &
                .AND. MPI_id==0) THEN
              write(out_unitp,'(a)',ADVANCE='no') '---'
              CALL flush_perso(out_unitp)
            END IF
          END DO
        ELSE
          ib0 = ib
          DO i=1,n
            ib = ib + 1
            IF (ib > para_Op%nb_tot) EXIT

            IF(MPI_id==0) THEN
              psi(i)%RvecB(:) = ZERO
              psi(i)%RvecB(ib) = ONE
              CALL Set_symab_OF_psiBasisRep(psi(i))
            ENDIF
          END DO

          CALL sub_TabOpPsi(Psi,HPsi,para_Op)

          ib = ib0
          DO i=1,n
            ib = ib + 1
            IF (ib > para_Op%nb_tot) EXIT

            !> Rmat assigned here
            IF(MPI_id==0) para_Op%Rmat(:,ib)  = Hpsi(i)%RvecB(:)

            IF (mod(ib,max(1,int(para_Op%nb_tot/10))) == 0 .AND. print_level > -1      &
                .AND. MPI_id==0) THEN
              write(out_unitp,'(a)',ADVANCE='no') '---'
              CALL flush_perso(out_unitp)
            END IF
          END DO
        END IF
      END DO

      IF (print_level > -1 .AND. MPI_id==0) THEN
        write(out_unitp,'(a)',ADVANCE='yes') '----]'
        CALL flush_perso(out_unitp)
      END IF

      IF (debug) THEN
        write(out_unitp,*) para_Op%name_Op,' non-symmetrized'
        IF (para_Op%cplx) THEN
          CALL Write_Mat(para_Op%Cmat,out_unitp,3)
        ELSE
          CALL Write_Mat(para_Op%Rmat,out_unitp,5)
        END IF
      END IF

!---- deallocation -------------------------------------
      IF(MPI_id==0) THEN
        DO i=1,n
          CALL dealloc_psi(Hpsi(i))
          CALL dealloc_psi(psi(i))
        END DO
        CALL dealloc_NParray(psi ,"psi", name_sub)
        CALL dealloc_NParray(Hpsi,"Hpsi",name_sub)
      ENDIF
!     ----------------------------------------------------------
      para_Op%Make_mat = .TRUE.

!----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'END ',name_sub
       END IF
!----------------------------------------------------------
      END SUBROUTINE sub_MatOp_direct1
!=======================================================================================

      SUBROUTINE sub_MatOp_direct1_old(para_Op)
      USE mod_system
!$    USE omp_lib, only : OMP_GET_THREAD_NUM

      USE mod_SetOp
      USE mod_OpPsi
      USE mod_psi_set_alloc

      USE mod_psi_Op
      USE mod_MPI
      IMPLICIT NONE


!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op), intent(inout) :: para_Op

!------ for OpenMP --------------------------------
      TYPE (param_psi), pointer :: psi(:),Hpsi(:)
      integer       :: i,nb_thread,ith,n

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter ::name_sub='sub_MatOp_direct1'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
        write(out_unitp,*) 'Build matrix of ',para_Op%nb_tot
        CALL flush_perso(out_unitp)
      END IF

      IF (.NOT. para_Op%alloc_mat)                                      &
             CALL alloc_para_Op(para_Op,Mat=.TRUE.,Grid=.FALSE.)
      IF (para_Op%mat_done) RETURN

      para_Op%Make_mat = .FALSE.
!     - scaling of Op ---------------------------------------
        para_Op%E0     = ZERO
        para_Op%Esc    = ONE
        para_Op%scaled = .FALSE.
!-----------------------------------------------------------
      IF (MatOp_omp /= 1) THEN
        nb_thread = 1
      ELSE
        nb_thread = MatOp_maxth
      END IF
      IF (print_level>-1) write(out_unitp,*) 'nb_thread in ',name_sub,' : ',nb_thread

!---- initialization -------------------------------------
      nullify(psi)
      CALL alloc_array(psi, (/nb_thread/),"psi", name_sub)
      nullify(Hpsi)
      CALL alloc_array(Hpsi,(/nb_thread/),"Hpsi",name_sub)

      DO ith=1,nb_thread
        CALL init_psi(psi(ith),para_Op,para_Op%cplx)
        CALL init_psi(Hpsi(ith),para_Op,para_Op%cplx)
        psi(ith) = ZERO
      END DO

!     ----------------------------------------------------------
!       - build H and H0

        IF (print_level > -1) THEN
          write(out_unitp,'(a)')              'MatOp(:,i) (%): [--0-10-20-30-40-50-60-70-80-90-100]'
          write(out_unitp,'(a)',ADVANCE='no') 'MatOp(:,i) (%): ['
          CALL flush_perso(out_unitp)
        END IF

        n=para_Op%nb_tot/nb_thread
        !$OMP parallel do &
        !$OMP default(none) &
        !$OMP shared(para_Op,psi,Hpsi,print_level,out_unitp,MPI_id) &
        !$OMP private(i,ith) &
        !$OMP num_threads(nb_thread)
        DO i=1,para_Op%nb_tot
          ith = 1
          !$ ith = omp_get_thread_num()+1
          !write(out_unitp,*) 'i,ith',i,ith ; CALL flush_perso(out_unitp)

          IF (psi(ith)%cplx) THEN
            psi(ith)%CvecB(:) = CZERO
            psi(ith)%CvecB(i) = CONE
            CALL Set_symab_OF_psiBasisRep(psi(ith))
            CALL sub_OpPsi(psi(ith),Hpsi(ith),para_Op)
            para_Op%Cmat(:,i)  = Hpsi(ith)%CvecB(:)
          ELSE
            psi(ith)%RvecB(:) = ZERO
            psi(ith)%RvecB(i) = ONE
            CALL Set_symab_OF_psiBasisRep(psi(ith))

            CALL sub_OpPsi(psi(ith),Hpsi(ith),para_Op)

            para_Op%Rmat(:,i)  = Hpsi(ith)%RvecB(:)
          END IF
          IF (mod(i,max(1,int(para_Op%nb_tot/10))) == 0 .AND. print_level > -1         &
              .AND. MPI_id==0) THEN
            write(out_unitp,'(a)',ADVANCE='no') '---'
            CALL flush_perso(out_unitp)
          END IF

        END DO
        !$OMP end parallel do

        IF (print_level > -1 .AND. MPI_id==0) THEN
          write(out_unitp,'(a)',ADVANCE='yes') '----]'
          CALL flush_perso(out_unitp)
        END IF

      IF (debug) THEN
        write(out_unitp,*) para_Op%name_Op,' non-symmetrized'
        IF (para_Op%cplx) THEN
          CALL Write_Mat(para_Op%Cmat,out_unitp,3)
        ELSE
          CALL Write_Mat(para_Op%Rmat,out_unitp,5)
        END IF
      END IF

!---- deallocation -------------------------------------
      DO ith=1,nb_thread
        CALL dealloc_psi(Hpsi(ith))
        CALL dealloc_psi(psi(ith))
      END DO
      CALL dealloc_array(psi ,"psi", name_sub)
      CALL dealloc_array(Hpsi,"Hpsi",name_sub)
!     ----------------------------------------------------------
      para_Op%Make_mat = .TRUE.

!----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'END ',name_sub
       END IF
!----------------------------------------------------------



      END SUBROUTINE sub_MatOp_direct1_old

      SUBROUTINE sub_MatOp_direct1_Overlap(para_Op)
      USE mod_system
      USE mod_SetOp
      USE mod_OpPsi
      USE mod_psi_set_alloc
      USE mod_psi_B_TO_G

      USE mod_psi_Op
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op), intent(inout) :: para_Op

!------ for OpenMP --------------------------------
      TYPE (param_psi) :: psi,Hpsi
      integer       :: i,j,n
      real (kind=Rkind) :: max_Sii,max_Sij

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter ::name_sub='sub_MatOp_direct1_Overlap'
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
        write(out_unitp,*) 'Build matrix of ',para_Op%nb_tot
        CALL flush_perso(out_unitp)
      END IF

      IF (.NOT. para_Op%alloc_mat)                                      &
             CALL alloc_para_Op(para_Op,Mat=.TRUE.,Grid=.FALSE.)
      IF (para_Op%mat_done) RETURN

      para_Op%Make_mat = .FALSE.
!     - scaling of Op ---------------------------------------
        para_Op%E0     = ZERO
        para_Op%Esc    = ONE
        para_Op%scaled = .FALSE.
!-----------------------------------------------------------

!---- initialization -------------------------------------
        CALL init_psi(psi,para_Op,para_Op%cplx)
        CALL init_psi(Hpsi,para_Op,para_Op%cplx)
        psi = ZERO

!     ----------------------------------------------------------
!       - build H and H0
        IF (print_level>-1) write(out_unitp,'(a)',ADVANCE='no') 'MatOp(:,i) (%): 0'
        CALL flush_perso(out_unitp)
        para_Op%Rmat(:,:) = ZERO
        DO i=1,para_Op%nb_tot

          IF (psi%cplx) THEN
            psi%CvecB(:) = CZERO
            psi%CvecB(i) = CONE
            CALL Set_symab_OF_psiBasisRep(psi)
            CALL sub_OpPsi(psi,Hpsi,para_Op)
            para_Op%Cmat(:,i)  = Hpsi%CvecB(:)
          ELSE
            psi%RvecB(:) = ZERO
            psi%RvecB(i) = ONE

            CALL Set_symab_OF_psiBasisRep(psi)
            !write(6,*) 'i',i,psi%RvecB ; flush(6)

            Hpsi = psi
            CALL sub_PsiBasisRep_TO_GridRep(Hpsi)
            !write(6,*) 'i',i,Hpsi%RvecG ; flush(6)
            CALL sub_PsiGridRep_TO_BasisRep(Hpsi)
            !write(6,*) 'i',i,Hpsi%RvecB ; flush(6)

            para_Op%Rmat(:,i)  = Hpsi%RvecB(:)
          END IF

          IF (mod(i,max(1,int(para_Op%nb_tot/10))) == 0 .AND. print_level>-1) THEN
            write(out_unitp,'(a,i3)',ADVANCE='no') ' -',                        &
              int(real(i,kind=Rkind)*HUNDRED/real(para_Op%nb_tot,kind=Rkind))
            CALL flush_perso(out_unitp)
          END IF

        END DO
        IF (print_level>-1) write(out_unitp,'(a)',ADVANCE='yes') ' - 100'
        CALL flush_perso(out_unitp)

!       - analysis of the overlap matrix
        IF (para_Op%cplx) THEN
          CALL sub_ana_cplxS(para_Op%Cmat,para_Op%nb_tot,max_Sii,max_Sij,.TRUE.)
        ELSE
          CALL sub_ana_S(para_Op%Rmat,para_Op%nb_tot,max_Sii,max_Sij,.TRUE.)
        END IF


      IF (debug) THEN
        IF (para_Op%cplx) THEN
          CALL Write_Mat(para_Op%Cmat,out_unitp,3)
        ELSE
          CALL Write_Mat(para_Op%Rmat,out_unitp,5)
        END IF
      END IF

!---- deallocation -------------------------------------
        CALL dealloc_psi(Hpsi)
        CALL dealloc_psi(psi)

!     ----------------------------------------------------------
      para_Op%Make_mat = .TRUE.

!----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'END ',name_sub
       END IF
!----------------------------------------------------------



      END SUBROUTINE sub_MatOp_direct1_Overlap
      SUBROUTINE sub_MatOp_Grid(para_Op)
      USE mod_system
      USE mod_SetOp
      USE mod_OpPsi
      USE mod_psi_set_alloc

      USE mod_psi_Op
      IMPLICIT NONE


!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op), intent(inout) :: para_Op

!------ for OpenMP --------------------------------
      TYPE (param_psi) :: psi,Hpsi
      integer          :: iq,nq

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter ::name_sub='sub_MatOp_Grid'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
        write(out_unitp,*) 'Build matrix of ',para_Op%nb_tot
        CALL flush_perso(out_unitp)
      END IF

!     - scaling of Op ---------------------------------------
        para_Op%E0     = ZERO
        para_Op%Esc    = ONE
        para_Op%scaled = .FALSE.
!-----------------------------------------------------------

!---- initialization -------------------------------------
        CALL init_psi(psi,para_Op,para_Op%cplx)
        CALL init_psi(Hpsi,para_Op,para_Op%cplx)
        psi = ZERO

      nq = psi%nb_qaie


!     ----------------------------------------------------------
!       - build H and H0
        IF (print_level > -1) THEN
          write(out_unitp,'(a)')              'MatOp(:,i) (%): [--0-10-20-30-40-50-60-70-80-90-100]'
          write(out_unitp,'(a)',ADVANCE='no') 'MatOp(:,i) (%): ['
          CALL flush_perso(out_unitp)
        END IF

        DO iq=1,nq

          IF (para_Op%cplx) THEN

            psi%CvecG(:)  = CZERO
            psi%CvecG(iq) = CONE

            CALL sub_OpPsi(Psi,HPsi,para_Op)

            para_Op%Cmat(:,iq)  = Hpsi%CvecG(:)

          ELSE

            psi%RvecG(:)  = ZERO
            psi%RvecG(iq) = ONE

            CALL sub_OpPsi(Psi,HPsi,para_Op)

            para_Op%Rmat(:,iq)  = Hpsi%RvecG(:)

          END IF

          IF (mod(iq,max(1,int(nq/10))) == 0 .AND. print_level > -1) THEN
            write(out_unitp,'(a)',ADVANCE='no') '---'
            CALL flush_perso(out_unitp)
          END IF


        END DO

        IF (print_level > -1) THEN
          write(out_unitp,'(a)',ADVANCE='yes') '----]'
          CALL flush_perso(out_unitp)
        END IF

      IF (debug) THEN
        write(out_unitp,*) para_Op%name_Op,' non-symmetrized'
        IF (para_Op%cplx) THEN
          CALL Write_Mat(para_Op%Cmat,out_unitp,3)
        ELSE
          CALL Write_Mat(para_Op%Rmat,out_unitp,5)
        END IF
      END IF

!---- deallocation -------------------------------------
      CALL dealloc_psi(Hpsi)
      CALL dealloc_psi(psi)
!     ----------------------------------------------------------
      para_Op%Make_mat = .TRUE.

!----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'END ',name_sub
       END IF
!----------------------------------------------------------



      END SUBROUTINE sub_MatOp_Grid
      SUBROUTINE sub_MatOp_Overlap_SG4(para_Op)
      USE mod_system
      USE mod_SetOp
      USE mod_OpPsi
      USE mod_psi_set_alloc
      USE mod_psi_B_TO_G

      USE mod_psi_Op
      USE mod_basis_BtoG_GtoB_SGType4
      IMPLICIT NONE


!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op), intent(inout) :: para_Op

!------ for OpenMP --------------------------------
      TYPE (param_psi) :: psi

      TYPE(Type_SmolyakRep)     :: SRep_i,SRep_j,Srep_w

      integer       :: i,j,n
      real (kind=Rkind) :: max_Sii,max_Sij
      logical, parameter :: Grid=.TRUE.
      !logical, parameter :: Grid=.FALSE.

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter ::name_sub='sub_MatOp_Overlap_SG4'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
        write(out_unitp,*) 'Build matrix of ',para_Op%nb_tot
        CALL flush_perso(out_unitp)
      END IF

      IF (.NOT. para_Op%alloc_mat)                                      &
             CALL alloc_para_Op(para_Op,Mat=.TRUE.,Grid=.FALSE.)
      IF (para_Op%mat_done) RETURN

      para_Op%Make_mat = .FALSE.
!     - scaling of Op ---------------------------------------
        para_Op%E0     = ZERO
        para_Op%Esc    = ONE
        para_Op%scaled = .FALSE.
!-----------------------------------------------------------

!---- initialization -------------------------------------
        CALL init_psi(psi,para_Op,para_Op%cplx)
        CALL alloc_psi(psi,BasisRep=.TRUE.)

!     ----------------------------------------------------------
!       - build H and H0
        IF (print_level>-1) write(out_unitp,'(a)',ADVANCE='no') 'MatOp(:,i) (%): 0'
        CALL flush_perso(out_unitp)
        para_Op%Rmat(:,:) = ZERO

          IF (para_Op%cplx) THEN
            STOP 'not yet cplx in sub_MatOp_Overlap_SG4'
          END IF


        IF (Grid) THEN
          Srep_w = Set_weight_TO_SmolyakRep(                            &
             para_Op%BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval, &
             para_Op%BasisnD%tab_basisPrimSG)
        END IF

        DO i=1,para_Op%nb_tot
          psi%RvecB(:)   = ZERO
          psi%RvecB(i)   = ONE

          CALL tabPackedBasis_TO_SmolyakRepBasis(SRep_i,psi%RvecB,      &
                                       para_Op%BasisnD%tab_basisPrimSG, &
                     para_Op%BasisnD%nDindB,para_Op%BasisnD%para_SGType2)

          IF (Grid) THEN
            CALL BSmolyakRep_TO3_GSmolyakRep(SRep_i,                    &
                                          para_Op%BasisnD%para_SGType2, &
                                        para_Op%BasisnD%tab_basisPrimSG)
          END IF

          DO j=1,para_Op%nb_tot
            psi%RvecB(:)   = ZERO
            psi%RvecB(j)   = ONE

            CALL tabPackedBasis_TO_SmolyakRepBasis(SRep_j,psi%RvecB,    &
                                       para_Op%BasisnD%tab_basisPrimSG, &
                     para_Op%BasisnD%nDindB,para_Op%BasisnD%para_SGType2)

            IF (Grid) THEN
               CALL BSmolyakRep_TO3_GSmolyakRep(SRep_j,                 &
                                          para_Op%BasisnD%para_SGType2, &
                                        para_Op%BasisnD%tab_basisPrimSG)

               para_Op%Rmat(j,i) = dot_product_SmolyakRep_Grid(SRep_j,  &
                                 SRep_i,Srep_w,para_Op%BasisnD%WeightSG)
            ELSE
               para_Op%Rmat(j,i) = dot_product_SmolyakRep_Basis(SRep_j, &
                                        SRep_i,para_Op%BasisnD%WeightSG)
            END IF
          END DO

          IF (mod(i,max(1,int(para_Op%nb_tot/10))) == 0 .AND. print_level>-1) THEN
            write(out_unitp,'(a,i3)',ADVANCE='no') ' -',                        &
              int(real(i,kind=Rkind)*HUNDRED/real(para_Op%nb_tot,kind=Rkind))
            CALL flush_perso(out_unitp)
          END IF

        END DO
        IF (print_level>-1) write(out_unitp,'(a)',ADVANCE='yes') ' - 100'
        CALL flush_perso(out_unitp)

!       - analysis of the overlap matrix
        IF (para_Op%cplx) THEN
          CALL sub_ana_cplxS(para_Op%Cmat,para_Op%nb_tot,max_Sii,max_Sij,.TRUE.)
        ELSE
          CALL sub_ana_S(para_Op%Rmat,para_Op%nb_tot,max_Sii,max_Sij,.TRUE.)
        END IF


      IF (debug) THEN
        IF (para_Op%cplx) THEN
          CALL Write_Mat(para_Op%Cmat,out_unitp,3)
        ELSE
          CALL Write_Mat(para_Op%Rmat,out_unitp,5)
        END IF
      END IF

!---- deallocation -------------------------------------
        CALL dealloc_psi(psi)

!     ----------------------------------------------------------
      para_Op%Make_mat = .TRUE.

!----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'END ',name_sub
       END IF
!----------------------------------------------------------



      END SUBROUTINE sub_MatOp_Overlap_SG4
      SUBROUTINE sub_MatOp_V_SG4(para_Op)
      USE mod_system
      USE mod_SetOp
      USE mod_OpPsi
      USE mod_psi_set_alloc
      USE mod_psi_B_TO_G

      USE mod_psi_Op
      USE mod_basis_BtoG_GtoB_SGType4
      IMPLICIT NONE


!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op), intent(inout) :: para_Op

!------ for OpenMP --------------------------------
      TYPE (param_psi) :: psi

      TYPE(Type_SmolyakRep)     :: SRep_i,SRep_j,Srep_w,Srep_V

      integer       :: i,j,n
      real (kind=Rkind) :: max_Sii,max_Sij

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter ::name_sub='sub_MatOp_V_SG4'
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
        write(out_unitp,*) 'Build matrix of ',para_Op%nb_tot
        CALL flush_perso(out_unitp)
      END IF

      IF (.NOT. para_Op%alloc_mat)                                      &
             CALL alloc_para_Op(para_Op,Mat=.TRUE.,Grid=.FALSE.)
      IF (para_Op%mat_done) RETURN

      para_Op%Make_mat = .FALSE.
!     - scaling of Op ---------------------------------------
        para_Op%E0     = ZERO
        para_Op%Esc    = ONE
        para_Op%scaled = .FALSE.
!-----------------------------------------------------------

!---- initialization -------------------------------------
        CALL init_psi(psi,para_Op,para_Op%cplx)
        CALL alloc_psi(psi,BasisRep=.TRUE.)


      IF (para_Op%BasisnD%para_SGType2%nb0 /= 1) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*)
        write(out_unitp,*) 'nb0 > 1',para_Op%BasisnD%para_SGType2%nb0
        STOP 'nb0 /= 1'
      END IF
!     ----------------------------------------------------------
!       - build H and H0
        IF (print_level>-1) write(out_unitp,'(a)',ADVANCE='no') 'MatOp(:,i) (%): 0'
        CALL flush_perso(out_unitp)
        para_Op%Rmat(:,:) = ZERO

          IF (para_Op%cplx) THEN
            STOP 'not yet cplx in sub_MatOp_Overlap_SG4'
          END IF

        CALL alloc2_SmolyakRep(SRep_V, &
                       para_Op%BasisnD%para_SGType2%nDind_SmolyakRep, &
                       para_Op%BasisnD%tab_basisPrimSG,grid=.TRUE.)

        CALL tabR2bis_TO_SmolyakRep1(SRep_V,para_Op%OpGrid(1)%Grid(:,1,1))

        Srep_w = Set_weight_TO_SmolyakRep(                            &
             para_Op%BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval, &
             para_Op%BasisnD%tab_basisPrimSG) * SRep_V




        DO i=1,para_Op%nb_tot
          psi%RvecB(:)   = ZERO
          psi%RvecB(i)   = ONE

          CALL tabPackedBasis_TO_SmolyakRepBasis(SRep_i,psi%RvecB,      &
                                       para_Op%BasisnD%tab_basisPrimSG, &
                     para_Op%BasisnD%nDindB,para_Op%BasisnD%para_SGType2)

           CALL BSmolyakRep_TO3_GSmolyakRep(SRep_i,                      &
                                           para_Op%BasisnD%para_SGType2, &
                                          para_Op%BasisnD%tab_basisPrimSG)

          DO j=1,para_Op%nb_tot
            psi%RvecB(:)   = ZERO
            psi%RvecB(j)   = ONE

            CALL tabPackedBasis_TO_SmolyakRepBasis(SRep_j,psi%RvecB,    &
                                       para_Op%BasisnD%tab_basisPrimSG, &
                     para_Op%BasisnD%nDindB,para_Op%BasisnD%para_SGType2)

             CALL BSmolyakRep_TO3_GSmolyakRep(SRep_j,                    &
                                           para_Op%BasisnD%para_SGType2, &
                                         para_Op%BasisnD%tab_basisPrimSG)

          para_Op%Rmat(j,i) = dot_product_SmolyakRep_Grid(SRep_j,SRep_i,&
                                        Srep_w,para_Op%BasisnD%WeightSG)

          END DO

          IF (mod(i,max(1,int(para_Op%nb_tot/10))) == 0 .AND. print_level>-1) THEN
            write(out_unitp,'(a,i3)',ADVANCE='no') ' -',                        &
              int(real(i,kind=Rkind)*HUNDRED/real(para_Op%nb_tot,kind=Rkind))
            CALL flush_perso(out_unitp)
          END IF

        END DO
        IF (print_level>-1) write(out_unitp,'(a)',ADVANCE='yes') ' - 100'
        CALL flush_perso(out_unitp)

!       - analysis of the overlap matrix
        IF (para_Op%cplx) THEN
          CALL sub_ana_cplxS(para_Op%Cmat,para_Op%nb_tot,max_Sii,max_Sij,.TRUE.)
        ELSE
          CALL sub_ana_S(para_Op%Rmat,para_Op%nb_tot,max_Sii,max_Sij,.TRUE.)
        END IF


      IF (debug) THEN
        IF (para_Op%cplx) THEN
          CALL Write_Mat(para_Op%Cmat,out_unitp,3)
        ELSE
          CALL Write_Mat(para_Op%Rmat,out_unitp,5)
        END IF
      END IF

!---- deallocation -------------------------------------
        CALL dealloc_psi(psi)

!     ----------------------------------------------------------
      para_Op%Make_mat = .TRUE.

!----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'END ',name_sub
       END IF
!----------------------------------------------------------



      END SUBROUTINE sub_MatOp_V_SG4
      SUBROUTINE sub_MatOp_OpExact_SG4(para_Op)
      USE mod_system
      USE mod_SetOp
      USE mod_OpPsi
      USE mod_psi_set_alloc
      USE mod_psi_B_TO_G

      USE mod_psi_Op
      USE mod_basis_BtoG_GtoB_SGType4
      IMPLICIT NONE


!----- variables pour la namelist minimum ----------------------------
 TYPE (param_Op), intent(inout) :: para_Op


 integer       :: n

 ! local variables
 TYPE (CoordType), pointer :: mole
 TYPE (basis),   pointer :: BasisnD


 integer                :: ib,jb,ipb,jpb,i,iG,ith,nb,iBSRep
 integer, allocatable   :: tab_l(:)

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter ::name_sub='sub_MatOp_OpExact_SG4'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
  mole    => para_Op%mole
  BasisnD => para_Op%BasisnD

      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
        write(out_unitp,*) 'Build matrix of ',para_Op%nb_tot
        CALL flush_perso(out_unitp)
      END IF

      IF (.NOT. para_Op%alloc_mat)                                      &
             CALL alloc_para_Op(para_Op,Mat=.TRUE.,Grid=.FALSE.)
      IF (para_Op%mat_done) RETURN

      para_Op%Make_mat = .FALSE.


!---- initialization -------------------------------------
      para_Op%Rmat(:,:) = ZERO
!     ----------------------------------------------------------

   !--------------------------------------------------------------
   !-- For the initialization of tab_l(:) and the use of ADD_ONE_TO_nDindex in the parallel loop
   CALL alloc_NParray(tab_l,(/ BasisnD%para_SGType2%nDind_SmolyakRep%ndim /),'tabl_l',name_sub)
   ith = 0
   tab_l(:) = BasisnD%para_SGType2%nDval_init(:,ith+1)
   !--------------------------------------------------------------

   !write(6,*) 'ith,tab_l(:)',ith,':',tab_l

   ! we are not using the parallel do, to be able to use the correct initialized tab_l with nDval_init
   DO iG=BasisnD%para_SGType2%iG_th(ith+1),BasisnD%para_SGType2%fG_th(ith+1)
     !write(6,*) 'iG',iG ; flush(6)

     CALL ADD_ONE_TO_nDindex(BasisnD%para_SGType2%nDind_SmolyakRep,tab_l,iG=iG)

     nb = BasisnD%para_SGType2%tab_nb_OF_SRep(iG)

     IF (iG == 1) THEN
       iBSRep = 0
     ELSE
       iBSRep = BasisnD%para_SGType2%tab_Sum_nb_OF_SRep(iG-1)
     END IF

     !write(6,*) 'iG,iBSRep,nb,WeightSG',iG,iBSRep,nb,BasisnD%WeightSG(iG)

     DO ib=1,nb
       ipb = BasisnD%para_SGType2%tab_iB_OF_SRep_TO_iB(iBSRep+ib)
       IF (ipb == 0) CYCLE
       !write(6,*) 'ib,ipb',ib,ipb

       DO jb=1,nb
         jpb = BasisnD%para_SGType2%tab_iB_OF_SRep_TO_iB(iBSRep+jb)
         IF (jpb == 0) CYCLE

         !write(6,*) 'jb,jpb',jb,jpb

         para_Op%Rmat(ipb,jpb) = para_Op%Rmat(ipb,jpb) + BasisnD%WeightSG(iG)

       END DO
     END DO

     !write(6,*) 'iG done:',iG ; flush(6)
   END DO
   CALL dealloc_NParray(tab_l,'tabl_l',name_sub)

   write(out_unitp,*) '# 1',count(para_Op%Rmat == 1)
   write(out_unitp,*) '# 0',count(para_Op%Rmat == 0)
   write(out_unitp,*) '# (diff 0 and diff 1)',size(para_Op%Rmat)-       &
                     (count(para_Op%Rmat == 0)+count(para_Op%Rmat == 1))

   CALL Write_Mat(para_Op%Rmat,out_unitp,5)

   write(6,*) 'for MatOp**2'
   para_Op%Rmat = matmul(para_Op%Rmat,para_Op%Rmat)

   write(out_unitp,*) '# 1',count(para_Op%Rmat == 1)
   write(out_unitp,*) '# 0',count(para_Op%Rmat == 0)
   write(out_unitp,*) '# (diff 0 and diff 1)',size(para_Op%Rmat)-       &
                   (count(para_Op%Rmat == 0)+count(para_Op%Rmat == 1))

   CALL Write_Mat(para_Op%Rmat,out_unitp,5)
!     ----------------------------------------------------------
      para_Op%Make_mat = .TRUE.

!----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'END ',name_sub
       END IF
!----------------------------------------------------------

      END SUBROUTINE sub_MatOp_OpExact_SG4
      
!===============================================================================     
      SUBROUTINE sub_OpBasisFi(para_Op,i)
      USE mod_system
      USE mod_SetOp
      USE mod_OpPsi
      USE mod_psi_set_alloc

      USE mod_psi_Op
      USE mod_MPI
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op), intent(inout) :: para_Op
      integer, intent(in)            :: i

!------ working parameters --------------------------------
      TYPE (param_psi)   :: psi,Hpsi


!----- for debuging --------------------------------------------------
      character (len=*), parameter ::name_sub='sub_OpBasisFi'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
        write(out_unitp,*) 'Build Op(:,i) ',para_Op%nb_tot
        CALL flush_perso(out_unitp)
      END IF

      IF (.NOT. para_Op%alloc_mat) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' The matrix of Op HAS to be allocated'
          write(out_unitp,*) ' CHECK the fortran source!'
          STOP
      END IF

      !write(out_unitp,*) 'in total memory: ',para_mem%mem_tot
!------ initialization -------------------------------------
      CALL init_psi(psi,para_Op,para_Op%cplx)
      CALL init_psi(Hpsi,para_Op,para_Op%cplx)
      psi = ZERO
      !write(out_unitp,*) 'psi%nb_qa',psi%nb_qa
!     ----------------------------------------------------------
      IF (psi%cplx) THEN
        psi%CvecB(i) = CONE

        CALL Set_symab_OF_psiBasisRep(psi)

        CALL sub_OpPsi(psi,Hpsi,para_Op)
        para_Op%Cmat(:,i)  = Hpsi%CvecB(:) !< Cmat calculated
      ELSE
        psi%RvecB(i) = ONE
        CALL Set_symab_OF_psiBasisRep(psi)
        CALL sub_OpPsi(psi,Hpsi,para_Op)
        para_Op%Rmat(:,i)  = Hpsi%RvecB(:) !< Rmat calculated
      END IF
      CALL dealloc_psi(Hpsi)
      CALL dealloc_psi(psi)
!     ----------------------------------------------------------
       !write(out_unitp,*) 'out total memory: ',para_mem%mem_tot
!----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'END ',name_sub
       END IF
!----------------------------------------------------------

      END SUBROUTINE sub_OpBasisFi
!===============================================================================     

END MODULE Mod_MatOp
