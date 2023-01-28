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
MODULE mod_OpPsi
USE mod_OpPsi_SG4
USE mod_OpPsi_SG4_MPI

PRIVATE
PUBLIC :: sub_PsiOpPsi,sub_OpPsi,sub_scaledOpPsi,sub_OpiPsi,sub_TabOpPsi
PUBLIC :: sub_PsiDia_TO_PsiAdia_WITH_MemGrid
!PUBLIC :: sub_TabOpPsi_OF_ONEDP_FOR_SGtype4
!PUBLIC :: sub_TabOpPsi_OF_ONEGDP_WithOp_FOR_SGtype4
PUBLIC :: sub_psiHitermPsi
PUBLIC :: sub_OpBasis_OneBF,sub_OpBasis_OneCBF

CONTAINS
!======================================================
!
!     E = <Psi | Op | Psi>
!
!======================================================
      FUNCTION skip_term(derOp,derOp_FROM_Qdyn)
      USE mod_system
      IMPLICIT NONE

      logical :: skip_term
!----- variables pour la namelist minimum ----------------------------
      integer, intent(in) :: derOp(2),derOp_FROM_Qdyn(2)

      integer :: derOp_FROM_Qdyn_loc(2)
      logical :: not_skip,li11_22,li12_21

      !----- for debuging ----------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='skip_term'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'derOp          : ',derOp
        write(out_unitp,*) 'derOp_FROM_Qdyn: ',derOp_FROM_Qdyn
      END IF


      li11_22 = ( derOp(1) == min(0,derOp_FROM_Qdyn(1)) ) .AND.         &
                ( derOp(2) == min(0,derOp_FROM_Qdyn(2)) )

      li12_21 = ( derOp(1) == min(0,derOp_FROM_Qdyn(2)) ) .AND.         &
                ( derOp(2) == min(0,derOp_FROM_Qdyn(1)) )
      not_skip = li11_22 .OR. li12_21

      !derOp_FROM_Qdyn_loc(:) = derOp_FROM_Qdyn(:)
      !WHERE (derOp_FROM_Qdyn_loc(:) > 0) derOp_FROM_Qdyn_loc(:) = 0

      !not_skip = ( sum(abs(derOp_FROM_Qdyn_loc-derOp)) == 0 ) .OR.      &
      !           ( sum(abs(cshift(derOp_FROM_Qdyn_loc,1)-derOp)) == 0)

      skip_term = .NOT. not_skip

      IF (debug) THEN
        write(out_unitp,*) 'skip_term: ',(.NOT. not_skip)
        write(out_unitp,*) 'END ',name_sub
      END IF
      !write(out_unitp,*) 'skip_term: ',derOp,derOp_FROM_Qdyn,(.NOT. not_skip)

      END FUNCTION skip_term


!======================================================
!
!     E = <Psi | Op | Psi>
!
!======================================================
      SUBROUTINE sub_PsiOpPsi(E,Psi,OpPsi,para_Op,iOp)
      USE mod_system
      USE mod_psi,             ONLY : param_psi,ecri_psi,Overlap_psi1_psi2
      USE mod_SetOp,           ONLY : param_Op,write_param_Op
      USE mod_MPI
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)                            :: para_Op
      integer,              intent(in), optional :: iOp
      complex (kind=Rkind), intent(inout)        :: E

!----- variables for the WP ----------------------------------------
      TYPE (param_psi)    :: Psi
      TYPE (param_psi)    :: OpPsi

  !----- for debuging --------------------------------------------------
  logical, parameter :: debug = .FALSE.
  !logical, parameter :: debug = .TRUE.
  !---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING sub_PsiOpPsi',para_Op%nb_tot,     &
                               para_Op%nb_ba,para_Op%nb_bi,para_Op%nb_be
        write(out_unitp,*) 'nb_act1',para_Op%mole%nb_act1
        write(out_unitp,*) 'nb_var',para_Op%mole%nb_var
        write(out_unitp,*) 'symab of para_Op,psi ',para_Op%symab,psi%symab
        IF (present(iOp)) write(out_unitp,*) 'iOp',iOp
        write(out_unitp,*)
        CALL write_param_Op(para_Op)
        IF (allocated(para_Op%Cmat)) CALL Write_Mat(para_Op%Cmat,out_unitp,3)
        IF (allocated(para_Op%Rmat)) CALL Write_Mat(para_Op%Rmat,out_unitp,5)
        write(out_unitp,*)
        write(out_unitp,*) 'PsiBasisRep'
        CALL ecri_psi(Psi=Psi)
        flush(out_unitp)
      END IF
!-----------------------------------------------------------

      IF (present(iOp)) THEN
        CALL sub_OpiPsi(Psi,OpPsi,para_Op,iOp)
      ELSE
        CALL sub_OpPsi(Psi,OpPsi,para_Op)
      END IF

      IF (debug) THEN
        write(out_unitp,*) 'Op.PsiBasisRep'
        CALL ecri_psi(Psi=OpPsi)
      END IF

      IF(keep_MPI) CALL Overlap_psi1_psi2(E,Psi,OpPsi)

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'E (au)',E
         write(out_unitp,*) 'END sub_PsiOpPsi'
       END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_PsiOpPsi

!===============================================================================
      SUBROUTINE sub_OpPsi(Psi,OpPsi,para_Op,derOp,With_Grid,TransfoOp)
      USE mod_system
      USE mod_psi,             ONLY : param_psi,ecri_psi,dealloc_psi
      USE mod_SetOp,           ONLY : param_Op,write_param_Op
      USE mod_MPI_aux
      IMPLICIT NONE

      !----- variables pour la namelist minimum ------------------------
      TYPE (param_Op)  :: para_Op
      integer, intent(in), optional :: derOp(2)
      logical, intent(in), optional :: With_Grid,TransfoOp
      integer          :: n


      !----- variables for the WP --------------------------------------
      TYPE (param_psi)   :: Psi,OpPsi
      integer            :: OpPsi_symab

      integer :: derOp_loc(2)
      logical :: With_Grid_loc,TransfoOp_loc

      TYPE (param_psi)   :: Psi_loc


      !----- for debuging ----------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_OpPsi'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      !-----------------------------------------------------------------
      n = para_Op%nb_tot
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub,' ',n
        write(out_unitp,*) 'nb_bie,nb_baie',para_Op%nb_bie,para_Op%nb_baie
        write(out_unitp,*) 'symab of para_Op,psi ',para_Op%symab,psi%symab
        write(out_unitp,*) 'para_Op%mat_done',para_Op%mat_done
        flush(out_unitp)
        CALL write_param_Op(para_Op)
        IF (allocated(para_Op%Cmat) .AND. para_Op%mat_done) CALL Write_Mat(para_Op%Cmat,out_unitp,3)
        IF (allocated(para_Op%Rmat) .AND. para_Op%mat_done) CALL Write_Mat(para_Op%Rmat,out_unitp,5)
        write(out_unitp,*)
        write(out_unitp,*) 'Psi'
        CALL ecri_psi(Psi=Psi)
        flush(out_unitp)
      END IF
      !-----------------------------------------------------------------
      IF (present(derOp)) THEN
        derOp_loc(:) = derOp(:)
      ELSE
        derOp_loc(:) = 0
      END IF

      IF (present(With_Grid)) THEN
        With_Grid_loc = With_Grid
      ELSE
        With_Grid_loc = .FALSE.
      END IF

      IF(keep_MPI) THEN
        IF (.NOT. psi%BasisnD%dnGGRep) With_Grid_loc = .FALSE.
      ENDIF

      IF(openmpi .AND. MPI_scheme/=1) CALL MPI_Bcast_(With_Grid_loc,size1_MPI,root_MPI)

      IF (present(TransfoOp)) THEN
        TransfoOp_loc = TransfoOp
      ELSE
        TransfoOp_loc = .FALSE.
      END IF
      !-------------------------------------------------------------------------

      IF (para_Op%para_ReadOp%Op_Transfo .AND. TransfoOp_loc) THEN
        CALL sub_PrimOpPsi(Psi,Psi_loc,para_Op,derOp_loc,With_Grid_loc)
        CALL sub_scaledOpPsi(Psi,Psi_loc,para_Op%para_ReadOp%E0_Transfo,ONE)
        CALL sub_PrimOpPsi(Psi_loc,OpPsi,para_Op,derOp_loc,With_Grid_loc)
        CALL sub_scaledOpPsi(Psi_loc,OpPsi,para_Op%para_ReadOp%E0_Transfo,ONE)
        CALL dealloc_psi(Psi_loc,delete_all=.TRUE.)
      ELSE
        ! MPI action here
        CALL sub_PrimOpPsi(Psi,OpPsi,para_Op,derOp_loc,With_Grid_loc)
      END IF


!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'OpPsiBasisRep'
        CALL ecri_psi(Psi=OpPsi)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_OpPsi
!===============================================================================

!===============================================================================
      SUBROUTINE sub_PrimOpPsi(Psi,OpPsi,para_Op,derOp,With_Grid,pot_only)
      USE mod_system
      USE mod_SymAbelian,     ONLY : Calc_symab1_EOR_symab2
      USE mod_psi,            ONLY : param_psi,ecri_psi,copy_psi2TOpsi1,&
                                     alloc_psi,dealloc_psi,             &
                                     Set_symab_OF_psiBasisRep,          &
                                     sub_PsiGridRep_TO_BasisRep
      USE mod_SetOp,          ONLY : param_Op,write_param_Op,read_OpGrid_OF_Op
      USE mod_OpPsi_SG4_MPI,  ONLY : sub_TabOpPsi_FOR_SGtype4_SRB_MPI
      USE mod_MPI_aux
      IMPLICIT NONE

      !----- variables pour la namelist minimum ------------------------
      TYPE (param_Op)  :: para_Op
      integer, intent(in), optional :: derOp(2)
      logical, intent(in), optional :: With_Grid,pot_only
      integer                       :: n


      !----- variables for the WP --------------------------------------
      TYPE (param_psi)   :: Psi,OpPsi
      integer :: OpPsi_symab

      integer :: derOp_loc(2)
      logical :: With_Grid_loc,pot_only_loc
      TYPE (param_psi)   :: RPsi,ROpPsi
      TYPE (param_psi)   :: RCPsi(2),RCOpPsi(2)
      TYPE (param_psi)   :: RROpPsi(1)

      logical :: direct_KEO,SGtype4


      !----- for debuging ------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_PrimOpPsi'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      !-------------------------------------------------------------------------
      n = para_Op%nb_tot
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub,' ',n
        write(out_unitp,*) 'nb_bie,nb_baie',para_Op%nb_bie,para_Op%nb_baie
        write(out_unitp,*) 'symab of para_Op,psi ',para_Op%symab,psi%symab
        write(out_unitp,*) 'para_Op%mat_done',para_Op%mat_done
        write(out_unitp,*) 'para_Op%... %Save_MemGrid_done',            &
                    para_Op%para_ReadOp%para_FileGrid%Save_MemGrid_done
        flush(out_unitp)
        !CALL write_param_Op(para_Op)
        !IF (allocated(para_Op%Cmat) .AND. para_Op%mat_done) CALL Write_Mat(para_Op%Cmat,out_unitp,3)
        !IF (allocated(para_Op%Rmat) .AND. para_Op%mat_done) CALL Write_Mat(para_Op%Rmat,out_unitp,5)
        write(out_unitp,*)
        write(out_unitp,*) 'Psi'
        CALL ecri_psi(Psi=Psi)
        flush(out_unitp)
      END IF
      !-------------------------------------------------------------------------

      IF (present(derOp)) THEN
        derOp_loc(:) = derOp(:)
      ELSE
        derOp_loc(:) = 0
      END IF

      IF (present(With_Grid)) THEN
        With_Grid_loc = With_Grid
      ELSE
        With_Grid_loc = .FALSE.
      END IF

      IF(keep_MPI) THEN
        IF (.NOT. psi%BasisnD%dnGGRep) With_Grid_loc = .FALSE.
      ENDIF

      IF(openmpi .AND. MPI_scheme/=1) CALL MPI_Bcast_(With_Grid_loc,size1_MPI,root_MPI)

      IF (present(pot_only)) THEN
        pot_only_loc = pot_only .AND. (para_Op%n_Op ==0) ! para_Op has to be H
      ELSE
        pot_only_loc = .FALSE.
      END IF

      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      !     - test on nb_bi : nb_bi>0
      IF (para_Op%nb_bi <= 0) THEN
        write(out_unitp,*) ' ERROR : nb_bi MUST be > 0'
        write(out_unitp,*) 'nb_bi',para_Op%nb_bi
        write(out_unitp,*) 'STOP in ',name_sub
        STOP
      END IF
      !-------------------------------------------------------------------------

      !-------------------------------------------------------------------------
      ! special case: to deal with TabPsi
      SGtype4    = (para_Op%BasisnD%SparseGrid_type == 4)

      direct_KEO = para_Op%para_ReadOp%para_FileGrid%Save_MemGrid
      IF (SGtype4) THEN
        direct_KEO = para_Op%BasisnD%dnGGRep
      ELSE
        direct_KEO = direct_KEO .AND. para_Op%BasisnD%dnGGRep
      END IF
      direct_KEO = direct_KEO .AND. (para_Op%type_Op == 10)
      direct_KEO = direct_KEO .AND. para_Op%direct_KEO

!SGtype4=.FALSE.

      IF (debug) write(out_unitp,*) 'SGtype4,direct_KEO',SGtype4,direct_KEO
      flush(out_unitp)

      ! Number of Hamiltonian operation counter
      para_Op%nb_OpPsi = para_Op%nb_OpPsi + 1

      IF (SGtype4) THEN
        IF(Psi%SRB_MPI) THEN

        IF(openmpi) THEN
          CALL sub_TabOpPsi_FOR_SGtype4_SRB_MPI([Psi],RROpPsi,para_Op)
          OpPsi=RROpPsi(1)  ! find a way to avoid it
        ENDIF

        ELSE IF (Psi%cplx) THEN

          IF(keep_MPI) RCPsi = Psi

          IF(openmpi) THEN
            CALL sub_TabOpPsi_FOR_SGtype4_MPI(RCPsi,RCOpPsi,para_Op)
          ELSE
            CALL sub_TabOpPsi_FOR_SGtype4(RCPsi,RCOpPsi,para_Op)

          ENDIF
          IF(keep_MPI) OpPsi = RCOpPsi

          CALL dealloc_psi(RCPsi(1),  delete_all=.TRUE.)
          CALL dealloc_psi(RCPsi(2),  delete_all=.TRUE.)
          CALL dealloc_psi(RCOpPsi(1),delete_all=.TRUE.)
          CALL dealloc_psi(RCOpPsi(2),delete_all=.TRUE.)

        ELSE
          IF(openmpi) THEN
            CALL sub_TabOpPsi_FOR_SGtype4_MPI( [Psi] ,RROpPsi,para_Op)
          ELSE
            CALL sub_TabOpPsi_FOR_SGtype4( [Psi] ,RROpPsi,para_Op)
          ENDIF
          IF(keep_MPI) OpPsi = RROpPsi(1)
          CALL dealloc_psi(RROpPsi(1),delete_all=.TRUE.)
        END IF
      ELSE

        IF(keep_MPI) THEN
          OpPsi = Psi    ! For the allocation of OpPsi

          IF (para_Op%mat_done) THEN
            CALL sub_OpPsi_WITH_MatOp(Psi,OpPsi,para_Op)
          ELSE
            ! With_Grid_loc F
            IF (With_Grid_loc) THEN
              IF (.NOT. Psi%GridRep) THEN
                write(out_unitp,*) ' ERROR in ',name_sub
                write(out_unitp,*) '      Psi%GridRep=F and With_Grid=T'
                write(out_unitp,*) 'The wavepacket MUST be on the grid'
                write(out_unitp,*) 'CHECK the fortran!'
                STOP
              END IF

              IF (.NOT. para_Op%para_ReadOp%para_FileGrid%Save_MemGrid_done) THEN
                CALL read_OpGrid_OF_Op(para_Op)
              END IF
              OpPsi = ZERO

              IF (para_Op%type_Op == 10 .AND. .NOT. pot_only_loc) THEN
                IF (psi%cplx) THEN

                  CALL copy_psi2TOpsi1(RPsi,Psi,alloc=.FALSE.)
                  RPsi%cplx = .FALSE.
                  CALL alloc_psi(RPsi,GridRep=.TRUE.)
                  ROpPsi = RPsi

                  ! Real part
                  RPsi%RvecG(:) = Real(Psi%CvecG(:),kind=Rkind)
                  CALL sub_OpPsi_WITH_MemGrid_BGG_Hamil10(RPsi,ROpPsi,para_Op,derOp_loc,.TRUE.,pot_only=pot_only_loc)
                  OpPsi%CvecG(:) = cmplx(ROpPsi%RvecG,kind=Rkind)

                  ! Imaginary part
                  RPsi%RvecG(:) = aimag(Psi%CvecG(:))
                  CALL sub_OpPsi_WITH_MemGrid_BGG_Hamil10(RPsi,ROpPsi,para_Op,derOp_loc,.TRUE.,pot_only=pot_only_loc)
                  OpPsi%CvecG(:) = OpPsi%CvecG + EYE * ROpPsi%RvecG

                  CALL dealloc_psi(RPsi)
                  CALL dealloc_psi(ROpPsi)

                ELSE
                  CALL sub_OpPsi_WITH_MemGrid_BGG_Hamil10(Psi,OpPsi,para_Op,derOp_loc,.TRUE.,pot_only=pot_only_loc )
                END IF ! for psi%cplx
              ELSE ! for para_Op%type_Op /= 10
                CALL sub_OpPsi_WITH_MemGrid_BGG(Psi,OpPsi,para_Op,derOp_loc,.TRUE.,pot_only=pot_only_loc)
              END IF

              !STOP 'OpPsi with Grid'
            ELSE ! for .NOT. With_Grid_loc
              !--- allocate OpPsi and psi%GridRep  -----------------------------
              Psi%GridRep=.TRUE.
              CALL alloc_psi(Psi)
              OpPsi%GridRep=.TRUE.
              CALL alloc_psi(OpPsi)
              OpPsi = ZERO

              !-----------------------------------------------------------------
              IF (para_Op%para_ReadOp%para_FileGrid%Save_MemGrid) THEN

                IF (.NOT. para_Op%para_ReadOp%para_FileGrid%Save_MemGrid_done) THEN
                   CALL read_OpGrid_OF_Op(para_Op)
                END IF

                ! psi%BasisnD%dnGGRep T
                ! para_Op%type_Op 1
                IF (psi%BasisnD%dnGGRep) THEN
                  IF (para_Op%type_Op == 10) THEN
                    IF (psi%cplx) THEN

                      CALL copy_psi2TOpsi1(RPsi,Psi,alloc=.FALSE.)
                      RPsi%cplx = .FALSE.
                      CALL alloc_psi(RPsi,BasisRep=.TRUE.)
                      ROpPsi = RPsi

                      ! Real part
                      RPsi%RvecB(:) = Real(Psi%CvecB(:),kind=Rkind)
                      !write(out_unitp,*) 'Real part of Psi'
                      !CALL ecri_psi(Psi=RPsi)
                      CALL sub_OpPsi_WITH_MemGrid_BGG_Hamil10(RPsi,ROpPsi,para_Op,derOp_loc,.FALSE.,pot_only=pot_only_loc)
                      !write(out_unitp,*) 'Real part of OpPsi'
                      !CALL ecri_psi(Psi=ROpPsi)
                      OpPsi%CvecG(:) = cmplx(ROpPsi%RvecG,kind=Rkind)

                      ! Imaginary part
                      RPsi%RvecB(:) = aimag(Psi%CvecB(:))
                      !write(out_unitp,*) 'Imag part of Psi'
                      !CALL ecri_psi(Psi=RPsi)
                      CALL sub_OpPsi_WITH_MemGrid_BGG_Hamil10(RPsi,ROpPsi,para_Op,derOp_loc,.FALSE.,pot_only=pot_only_loc)
                      !write(out_unitp,*) 'Imag part of OpPsi'
                      !CALL ecri_psi(Psi=ROpPsi)
                      OpPsi%CvecG(:) = OpPsi%CvecG + EYE * ROpPsi%RvecG

                      CALL dealloc_psi(RPsi)
                      CALL dealloc_psi(ROpPsi)


                    ELSE
                      CALL sub_OpPsi_WITH_MemGrid_BGG_Hamil10(Psi,OpPsi,para_Op,derOp_loc,.FALSE.,pot_only=pot_only_loc)
                    END IF
                  ELSE
                    CALL sub_OpPsi_WITH_MemGrid_BGG(Psi,OpPsi,para_Op,derOp_loc,.FALSE.,pot_only=pot_only_loc)
                  END IF
                ELSE ! for .NOT. psi%BasisnD%dnGGRep
                  CALL sub_OpPsi_WITH_MemGrid(Psi,OpPsi,para_Op,derOp_loc,pot_only=pot_only_loc)
                END IF ! for psi%BasisnD%dnGGRep

              ELSE
                IF (para_Op%para_ReadOp%para_FileGrid%Type_FileGrid == 0) THEN
                  CALL sub_OpPsi_WITH_FileGrid_type0(Psi,OpPsi,para_Op,derOp_loc,pot_only=pot_only_loc)
                ELSE ! Type_FileGrid 1 or 2
                  IF (psi%BasisnD%dnGGRep) THEN
                    CALL sub_OpPsi_WITH_FileGrid_type12_BGG(Psi,OpPsi,para_Op,derOp_loc,.FALSE.,pot_only=pot_only_loc)
                  ELSE
                    CALL sub_OpPsi_WITH_FileGrid_type12(Psi,OpPsi,para_Op,derOp_loc,pot_only=pot_only_loc)
                  END IF
                END IF

              END IF

              !- the projection of PsiGridRep on PsiBasisRep -------------------
              CALL sub_PsiGridRep_TO_BasisRep(OpPsi)
              nb_mult_OpPsi = nb_mult_OpPsi + nb_mult_GTOB

              !--------------------------------------------------------
              !--------------------------------------------------------
              Psi%GridRep=.FALSE.
              CALL alloc_psi(Psi)
              OpPsi%GridRep=.FALSE.
              CALL alloc_psi(OpPsi)
              !--------------------------------------------------------

            END IF ! for With_Grid_loc
          END IF ! for para_Op%mat_done
        ENDIF ! for keep_MPI
      END IF

!=====================================================================
!
!      Symmetrization of  Op.Psi (if psi is on the basis)
!
!=====================================================================
      IF (debug) write(out_unitp,*) 'symab of para_Op,psi ',para_Op%symab,psi%symab

      IF(keep_MPI) THEN
        OpPsi_symab = Calc_symab1_EOR_symab2(para_Op%symab,psi%symab)
        CALL Set_symab_OF_psiBasisRep(OpPsi,OpPsi_symab)
      ENDIF
      IF (debug) write(out_unitp,*) 'OpPsi_symab',OpPsi%symab

!write(out_unitp,*) 'para_Op,psi symab ',para_Op%symab,psi%symab
!write(out_unitp,*) 'OpPsi_symab',OpPsi%symab
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'OpPsiBasisRep'
        CALL ecri_psi(Psi=OpPsi)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_PrimOpPsi

      SUBROUTINE sub_OpBasis_OneBF(Psi,OpPsi,para_Op,i)
      USE mod_system
      USE mod_psi,     ONLY : param_psi,Set_symab_OF_psiBasisRep
      USE mod_SetOp
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_psi), intent(inout) :: Psi,OpPsi
      TYPE (param_Op),  intent(in)    :: para_Op
      integer,          intent(in)    :: i        ! index of the basis function

!----- for debuging --------------------------------------------------
      character (len=*), parameter ::name_sub='sub_OpBasis_OneBF'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
        write(out_unitp,*) 'Build Op(:,i) ',para_Op%nb_tot
        flush(out_unitp)
      END IF

      psi = ZERO
      IF (psi%cplx) THEN
        psi%CvecB(i) = CONE
      ELSE
        psi%RvecB(i) = ONE
      END IF

      CALL Set_symab_OF_psiBasisRep(psi)
      CALL sub_OpPsi(psi,OpPsi,para_Op)

!----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'END ',name_sub
       END IF
!----------------------------------------------------------

      END SUBROUTINE sub_OpBasis_OneBF
      SUBROUTINE sub_OpBasis_OneCBF(Psi,OpPsi,para_Op,i)
      USE mod_system
      USE mod_psi,     ONLY : param_psi,Set_symab_OF_psiBasisRep,dealloc_psi
      USE mod_SetOp
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_psi), intent(inout) :: Psi,OpPsi
      TYPE (param_Op),  intent(in)    :: para_Op
      integer,          intent(in)    :: i        ! index of the basis function

!------ working parameters --------------------------------
      integer :: nb,nbc
      real (kind=Rkind),    allocatable :: RBC(:)
      complex (kind=Rkind), allocatable :: CBC(:)
      TYPE (param_psi)                  :: Temp_OpPsi

!----- for debuging --------------------------------------------------
      character (len=*), parameter ::name_sub='sub_OpBasis_OneCBF'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
        write(out_unitp,*) 'Build Op(:,i) ',para_Op%nb_tot
        flush(out_unitp)
      END IF

      CALL init_psi(Temp_OpPsi,para_Op,para_Op%cplx)
      Temp_OpPsi%nb_tot = para_Op%nb_tot_ini

      nbc = para_Op%nb_tot
      nb  = psi%nb_tot

      IF (para_Op%cplx) THEN
        CALL alloc_NParray(CBC,[NBC], 'CBC', name_sub)
        CBC(i) = CONE
        CALL CVecBC_TO_CvecB(CBC,psi%CvecB,nbc,nb,psi%BasisnD)

        CALL Set_symab_OF_psiBasisRep(psi)
        CALL sub_OpPsi(psi,Temp_OpPsi,para_Op)

        CALL CVecB_TO_CvecBC(Temp_OpPsi%CvecB,OpPsi%CvecB,nb,nbc,psi%BasisnD)

        CALL dealloc_NParray(CBC, 'CBC', name_sub)

      ELSE
        CALL alloc_NParray(RBC,[NBC], 'RBC', name_sub)
        RBC(i) = ONE
        CALL RVecBC_TO_RvecB(RBC,psi%RvecB,nbc,nb,psi%BasisnD)

        CALL Set_symab_OF_psiBasisRep(psi)
        CALL sub_OpPsi(psi,Temp_OpPsi,para_Op)

        CALL RVecB_TO_RvecBC(Temp_OpPsi%RvecB,OpPsi%RvecB,nb,nbc,psi%BasisnD)

        CALL dealloc_NParray(RBC, 'RBC', name_sub)

      END IF

      CALL dealloc_psi(Temp_OpPsi)
!     ----------------------------------------------------------
       !write(out_unitp,*) 'out total memory: ',para_mem%mem_tot
!----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'END ',name_sub
       END IF
!----------------------------------------------------------

END SUBROUTINE sub_OpBasis_OneCBF

!=======================================================================================

!=======================================================================================
!     | OpPsi> = Op | Psi>
!=======================================================================================
      SUBROUTINE sub_TabOpPsi(TabPsi,TabOpPsi,para_Op,derOp,With_Grid,TransfoOp)
      USE mod_system
      USE mod_psi,        ONLY : param_psi,ecri_psi,dealloc_psi
      USE mod_SetOp,      ONLY : param_Op,write_param_Op
      IMPLICIT NONE

      !----- variables pour la namelist minimum ------------------------
      TYPE (param_Op)  :: para_Op
      integer, intent(in), optional :: derOp(2)
      logical, intent(in), optional :: With_Grid,TransfoOp
      integer          :: n


      !----- variables for the WP --------------------------------------
      TYPE (param_psi)   :: TabPsi(:),TabOpPsi(:)
      integer :: OpPsi_symab

      integer :: i,derOp_loc(2)
      logical :: With_Grid_loc,TransfoOp_loc

      TYPE (param_psi)   :: Psi_loc


      !----- for debuging ----------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_TabOpPsi'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      !-----------------------------------------------------------------
      n = para_Op%nb_tot
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub,' ',n
        write(out_unitp,*) 'nb_bie,nb_baie',para_Op%nb_bie,para_Op%nb_baie
        write(out_unitp,*) 'para_Op%mat_done',para_Op%mat_done
        flush(out_unitp)
        CALL write_param_Op(para_Op)
        IF (allocated(para_Op%Cmat) .AND. para_Op%mat_done) CALL Write_Mat(para_Op%Cmat,out_unitp,3)
        IF (allocated(para_Op%Rmat) .AND. para_Op%mat_done) CALL Write_Mat(para_Op%Rmat,out_unitp,5)
        write(out_unitp,*)
        write(out_unitp,*) 'TabPsi'
        DO i=1,size(TabPsi)
          CALL ecri_psi(Psi=TabPsi(i))
        END DO
        flush(out_unitp)
      END IF
      !-----------------------------------------------------------------
!CALL Check_mem()

      !-----------------------------------------------------------------
      IF (present(derOp)) THEN
        derOp_loc(:) = derOp(:)
      ELSE
        derOp_loc(:) = 0
      END IF

      IF (present(With_Grid)) THEN
        With_Grid_loc = With_Grid
      ELSE
        With_Grid_loc = .FALSE.
      END IF

      IF (present(TransfoOp)) THEN
        TransfoOp_loc = TransfoOp
      ELSE
        TransfoOp_loc = .FALSE.
      END IF

      IF (para_Op%para_ReadOp%Op_Transfo .AND. TransfoOp_loc) THEN
        DO i=1,size(TabPsi)
          CALL sub_PrimOpPsi(TabPsi(i),Psi_loc,para_Op,derOp_loc,With_Grid_loc)
          CALL sub_scaledOpPsi(TabPsi(i),Psi_loc,para_Op%para_ReadOp%E0_Transfo,ONE)
          CALL sub_PrimOpPsi(Psi_loc,TabOpPsi(i),para_Op,derOp_loc,With_Grid_loc)
          CALL sub_scaledOpPsi(Psi_loc,TabOpPsi(i),para_Op%para_ReadOp%E0_Transfo,ONE)
          CALL dealloc_psi(Psi_loc,delete_all=.TRUE.)
        END DO
      ELSE
        ! MPI action here
        CALL sub_PrimTabOpPsi(TabPsi,TabOpPsi,para_Op,derOp_loc,With_Grid_loc)
      END IF

!CALL UnCheck_mem()
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'TabOpPsiBasisRep'
        DO i=1,size(TabPsi)
          CALL ecri_psi(Psi=TabOpPsi(i))
        END DO
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_TabOpPsi
!=======================================================================================

!=======================================================================================
      SUBROUTINE sub_PrimTabOpPsi(TabPsi,TabOpPsi,para_Op,derOp,With_Grid)
      USE mod_system
      USE mod_SymAbelian,    ONLY : Calc_symab1_EOR_symab2
      USE mod_psi,           ONLY : param_psi,alloc_psi,dealloc_psi,    &
                                    ecri_psi,Set_symab_OF_psiBasisRep,  &
                                    sub_PsiGridRep_TO_BasisRep
      USE mod_SetOp,         ONLY : param_Op,read_OpGrid_OF_Op,Write_FileGrid
      USE mod_MPI_aux
      IMPLICIT NONE

      !----- variables pour la namelist minimum ------------------------
      TYPE (param_Op)     :: para_Op
      integer, intent(in) :: derOp(2)
      logical, intent(in) :: With_Grid
      integer             :: n
      integer             :: size_TabPsi

      !----- variables for the WP --------------------------------------
      TYPE (param_psi)   :: TabPsi(:),TabOpPsi(:)
      integer :: OpPsi_symab

      integer :: i
      logical :: direct_KEO,SGtype4


      !----- for debuging ----------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_PrimTabOpPsi'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      !-----------------------------------------------------------------
      If(keep_MPI) size_TabPsi=size(TabPsi)

      IF(openmpi .AND. MPI_scheme/=1) CALL MPI_Bcast_(size_TabPsi,size1_MPI,root_MPI)

      n = para_Op%nb_tot
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub,' ',n
        write(out_unitp,*) 'nb_bie,nb_baie',para_Op%nb_bie,para_Op%nb_baie
        write(out_unitp,*) 'para_Op%mat_done',para_Op%mat_done
        flush(out_unitp)
        write(out_unitp,*)
        write(out_unitp,*) 'TabPsi'
        !DO i=1,size(TabPsi)
        DO i=1,size_TabPsi
          CALL ecri_psi(Psi=TabPsi(i))
        END DO
        flush(out_unitp)
      END IF
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      !- test on nb_bi : nb_bi>0
      IF (para_Op%nb_bi <= 0) THEN
        write(out_unitp,*) ' ERROR : nb_bi MUST be > 0'
        write(out_unitp,*) 'nb_bi',para_Op%nb_bi
        write(out_unitp,*) 'STOP in ',name_sub
        STOP
      END IF
      !-----------------------------------------------------------------

      ! special case: to deal with TabPsi
      SGtype4    = (para_Op%BasisnD%SparseGrid_type == 4)

      direct_KEO = para_Op%para_ReadOp%para_FileGrid%Save_MemGrid
      IF (SGtype4) THEN
        direct_KEO = para_Op%BasisnD%dnGGRep
      ELSE
        direct_KEO = direct_KEO .AND. para_Op%BasisnD%dnGGRep
      END IF
      direct_KEO = direct_KEO .AND. (para_Op%type_Op == 10)
      direct_KEO = direct_KEO .AND. para_Op%direct_KEO

      !IF (SGtype4 .AND. direct_KEO) THEN

      IF (SGtype4) THEN
        ! para_Op%BasisnD%SparseGrid_type=4
        IF(openmpi) THEN
          CALL sub_TabOpPsi_FOR_SGtype4_MPI(TabPsi,TabOpPsi,para_Op)
        ELSE
          CALL sub_TabOpPsi_FOR_SGtype4(TabPsi,TabOpPsi,para_Op)
        ENDIF
        !para_Op%nb_OpPsi = para_Op%nb_OpPsi + size(TabPsi)
        para_Op%nb_OpPsi=para_Op%nb_OpPsi+size_TabPsi
        RETURN
      END IF

      IF (debug) write(out_unitp,*) 'SGtype4,direct_KEO',SGtype4,direct_KEO
      flush(out_unitp)

      ! direct_KEO false by defult
      ! note this part is not modified for MPI yet
      IF (direct_KEO) THEN
        IF (.NOT. para_Op%para_ReadOp%para_FileGrid%Save_MemGrid_done) THEN
          CALL read_OpGrid_OF_Op(para_Op)
        END IF

        IF (SGtype4) THEN
          IF(openmpi) THEN
            CALL sub_TabOpPsi_FOR_SGtype4_MPI(TabPsi,TabOpPsi,para_Op)
          ELSE
            CALL sub_TabOpPsi_FOR_SGtype4(TabPsi,TabOpPsi,para_Op)
          ENDIF
          ! para_Op%nb_OpPsi = para_Op%nb_OpPsi + size(TabPsi)
          para_Op%nb_OpPsi = para_Op%nb_OpPsi+size_TabPsi
        ELSE
          CALL sub_TabOpPsi_WITH_MemGrid_BGG_Hamil10(TabPsi,TabOpPsi,para_Op,derOp,.FALSE.)
          ! para_Op%nb_OpPsi = para_Op%nb_OpPsi + size(TabPsi)
          para_Op%nb_OpPsi = para_Op%nb_OpPsi+size_TabPsi

          !DO i=1,size(TabPsi)
          DO i=1,size_TabPsi

            !- the projection of PsiGridRep on PsiBasisRep -------------------
            CALL sub_PsiGridRep_TO_BasisRep(TabOpPsi(i))
            nb_mult_OpPsi = nb_mult_OpPsi + nb_mult_GTOB

            !--------------------------------------------------------
            TabPsi(i)%GridRep=.FALSE.
            CALL alloc_psi(TabPsi(i))
            TabOpPsi(i)%GridRep=.FALSE.
            CALL alloc_psi(TabOpPsi(i))
            !--------------------------------------------------------

            OpPsi_symab = Calc_symab1_EOR_symab2(para_Op%symab,TabPsi(i)%symab)
            CALL Set_symab_OF_psiBasisRep(TabOpPsi(i),OpPsi_symab)

            !write(out_unitp,*) 'para_Op,psi symab ',para_Op%symab,TabPsi(i)%symab
            !write(out_unitp,*) 'OpPsi_symab',TabOpPsi(i)%symab

          END DO
        END IF

      ELSE
        !DO i=1,size(TabPsi)
        DO i=1,size_TabPsi
          CALL sub_PrimOpPsi(TabPsi(i),TabOpPsi(i),para_Op,derOp,With_Grid)
        END DO
      END IF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'OpPsiBasisRep'
        write(out_unitp,*) 'TabOpPsiBasisRep'
        DO i=1,size(TabOpPsi)
          CALL ecri_psi(Psi=TabOpPsi(i))
        END DO
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_PrimTabOpPsi
!=======================================================================================

!=======================================================================================
      SUBROUTINE sub_OpPsi_WITH_MatOp(Psi,OpPsi,para_Op)
      USE mod_system
      USE mod_psi,        ONLY : param_psi,ecri_psi
      USE mod_SetOp,      ONLY : param_Op,write_param_Op
      IMPLICIT NONE

      !----- variables pour la namelist minimum ------------------------
      TYPE (param_Op)  :: para_Op

      !----- variables for the WP --------------------------------------
      TYPE (param_psi)   :: Psi,OpPsi

      integer :: i,ki,k,n

      !----- for debuging ----------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_OpPsi_WITH_MatOp'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      !-----------------------------------------------------------------
      n = para_Op%nb_tot
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub,' ',n
        write(out_unitp,*) 'nb_bie,nb_baie',para_Op%nb_bie,para_Op%nb_baie
        write(out_unitp,*) 'nb_act1',para_Op%mole%nb_act1
        write(out_unitp,*) 'nb_var',para_Op%mole%nb_var
        write(out_unitp,*) 'para_Op%mat_done',para_Op%mat_done
        flush(out_unitp)
        CALL write_param_Op(para_Op)
        IF (allocated(para_Op%Cmat) .AND. para_Op%mat_done) CALL Write_Mat(para_Op%Cmat,out_unitp,3)
        IF (allocated(para_Op%Rmat) .AND. para_Op%mat_done) CALL Write_Mat(para_Op%Rmat,out_unitp,5)
        write(out_unitp,*)
        write(out_unitp,*) 'PsiBasisRep'
        CALL ecri_psi(Psi=Psi)
        !write(out_unitp,*) 'ini OpPsiBasisRep'
        !CALL ecri_psi(Psi=OpPsi)
        flush(out_unitp)
      END IF
      !-----------------------------------------------------------------

      IF (para_Op%mat_done) THEN

      !=================================================================
      !     Op.Psi with the Op as a matrix
      !=================================================================
        !- Calculation of the matrix (only once) -----------------------
        !CALL sub_MatOp(para_Op,.FALSE.)

        !- Op psi (with a matrix) --------------------------------------
        !---------------------------------------------------------------
        IF (para_Op%pack_Op .AND. para_Op%ratio_pack < para_Op%tol_nopack) THEN
          !- Op is packed ----------------------------------------------
          IF (Psi%cplx .AND. para_Op%cplx) THEN
            DO i=1,Psi%nb_tot
              OpPsi%CvecB(i) = CZERO
              DO k=1,para_Op%dim_Op(i)
                ki = para_Op%ind_Op(k,i)
                OpPsi%CvecB(i) = OpPsi%CvecB(i) +                       &
                                        para_Op%Cmat(i,ki)*Psi%CvecB(ki)
              END DO
            END DO
          ELSE IF (Psi%cplx .AND. .NOT. para_Op%cplx) THEN
            DO i=1,Psi%nb_tot
              OpPsi%CvecB(i) = CZERO
              DO k=1,para_Op%dim_Op(i)
                ki = para_Op%ind_Op(k,i)
                OpPsi%CvecB(i) = OpPsi%CvecB(i) +                       &
                                        para_Op%Rmat(i,ki)*Psi%CvecB(ki)
              END DO
            END DO
          ELSE
            DO i=1,Psi%nb_tot
              OpPsi%RvecB(i) = ZERO
              DO k=1,para_Op%dim_Op(i)
                ki = para_Op%ind_Op(k,i)
                OpPsi%RvecB(i) = OpPsi%RvecB(i) +                       &
                                        para_Op%Rmat(i,ki)*Psi%RvecB(ki)
              END DO
            END DO
          END IF
        ELSE
          !- Op is NOT packed ------------------------------------------
          IF (Psi%cplx) THEN
            IF (para_Op%cplx) THEN
              OpPsi%CvecB(:) = matmul(para_Op%Cmat,Psi%CvecB)
            ELSE
              OpPsi%CvecB(:) = matmul(para_Op%Rmat,Psi%CvecB)
            END IF
          ELSE
            OpPsi%RvecB(:) = matmul(para_Op%Rmat,Psi%RvecB)
          END IF
        END IF
      ELSE
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'Make_Mat=.FALSE. is not possible in this subroutine'
        write(out_unitp,*) 'Check the fortran'
        STOP
      END IF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'OpPsiBasisRep'
        CALL ecri_psi(Psi=OpPsi)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_OpPsi_WITH_MatOp
!=======================================================================================

      SUBROUTINE sub_OpPsi_WITH_MemGrid(Psi,OpPsi,para_Op,derOp,pot_only)
      USE mod_system
!$    USE omp_lib,      only : OMP_GET_THREAD_NUM
      USE mod_psi,      ONLY : param_psi,ecri_psi,                      &
                               alloc_psi,alloc_array,dealloc_array,     &
                               sub_d0d1d2PsiBasisRep_TO_GridRep,        &
                               sub_PsiBasisRep_TO_GridRep
      USE mod_SetOp,           ONLY : param_Op,write_param_Op
      IMPLICIT NONE

      !----- variables pour la namelist minimum ------------------------
      TYPE (param_Op)     :: para_Op
      integer, intent(in) :: derOp(2)
      logical, intent(in) :: pot_only

      integer          :: n

      integer          :: iOp

      !----- variables for the WP --------------------------------------
      TYPE (param_psi)   :: Psi,OpPsi
      integer :: OpPsi_symab
      TYPE (param_psi), pointer :: tab_Psi(:)


      !----- working variables -----------------------------------------
      integer :: i1_bi,i2_bi
      integer :: i_qa,iqi1,fqi1,iqi2,fqi2
      integer :: i,ki,k

      integer                  :: nio,error,iterm


      !------ for OpenMP -----------------------------------------------
      TYPE (param_psi), pointer :: thread_Psi(:),thread_OpPsi(:)
      integer       :: nb_thread,ith,ithread


      !----- for debuging ----------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_OpPsi_WITH_MemGrid'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      !-----------------------------------------------------------------
      n = para_Op%nb_tot
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub,' ',n
        write(out_unitp,*) 'nb_bie,nb_baie',para_Op%nb_bie,para_Op%nb_baie
        write(out_unitp,*) 'nb_act1',para_Op%mole%nb_act1
        write(out_unitp,*) 'nb_var',para_Op%mole%nb_var
        write(out_unitp,*) 'para_Op%mat_done',para_Op%mat_done
        flush(out_unitp)
        CALL write_param_Op(para_Op)
        write(out_unitp,*)
        write(out_unitp,*) 'PsiBasisRep'
        CALL ecri_psi(Psi=Psi)
        write(out_unitp,*) 'ini OpPsiBasisRep'
        CALL ecri_psi(Psi=OpPsi)
        flush(out_unitp)
      END IF
      !-----------------------------------------------------------------


      IF (para_Op%para_ReadOp%para_FileGrid%Save_MemGrid_done) THEN

        IF (OpPsi_omp /= 1) THEN
          nb_thread = 1
        ELSE
          nb_thread = OpPsi_maxth
        END IF
        IF (debug) write(out_unitp,*) 'nb_thread in ',name_sub,' : ',nb_thread

        !-----------------------------------------------------------------
        IF (nb_thread == 1) THEN

          DO iterm=1,para_Op%nb_term

            IF (para_Op%OpGrid(iterm)%grid_zero) CYCLE
            IF ( skip_term(derOp,para_Op%derive_termQdyn(:,iterm)) ) CYCLE
            IF (pot_only .AND. iterm /= para_Op%derive_term_TO_iterm(0,0)) CYCLE


            IF (para_Op%OpGrid(iterm)%grid_cte) THEN

              CALL sub_d0d1d2PsiBasisRep_TO_GridRep(Psi,para_Op%derive_termQdyn(:,iterm))

              DO i1_bi=1,para_Op%nb_bie
              DO i2_bi=1,para_Op%nb_bie
                iqi1 = 1         + (i1_bi-1) * Psi%nb_qa
                fqi1 = Psi%nb_qa + (i1_bi-1) * Psi%nb_qa
                iqi2 = 1         + (i2_bi-1) * Psi%nb_qa
                fqi2 = Psi%nb_qa + (i2_bi-1) * Psi%nb_qa


                IF (Psi%cplx) THEN
                  OpPsi%CvecG(iqi1:fqi1) = OpPsi%CvecG(iqi1:fqi1) +     &
                     Psi%CvecG(iqi2:fqi2) *                             &
                     para_Op%OpGrid(iterm)%Mat_cte(i1_bi,i2_bi)
                ELSE
                  OpPsi%RvecG(iqi1:fqi1) = OpPsi%RvecG(iqi1:fqi1) +     &
                     Psi%RvecG(iqi2:fqi2) *                             &
                     para_Op%OpGrid(iterm)%Mat_cte(i1_bi,i2_bi)
                END IF
              END DO
              END DO
            ELSE
              CALL sub_d0d1d2PsiBasisRep_TO_GridRep(Psi,para_Op%derive_termQdyn(:,iterm))

              !write(out_unitp,*) 'iterm,Grid',iterm,Grid
              DO i1_bi=1,para_Op%nb_bie
              DO i2_bi=1,para_Op%nb_bie
                iqi1 = 1         + (i1_bi-1) * Psi%nb_qa
                fqi1 = Psi%nb_qa + (i1_bi-1) * Psi%nb_qa
                iqi2 = 1         + (i2_bi-1) * Psi%nb_qa
                fqi2 = Psi%nb_qa + (i2_bi-1) * Psi%nb_qa

                IF (Psi%cplx) THEN
                  OpPsi%CvecG(iqi1:fqi1) = OpPsi%CvecG(iqi1:fqi1) +     &
                                                Psi%CvecG(iqi2:fqi2) *  &
                                para_Op%OpGrid(iterm)%Grid(:,i1_bi,i2_bi)
                ELSE
                  OpPsi%RvecG(iqi1:fqi1) = OpPsi%RvecG(iqi1:fqi1) +     &
                                               Psi%RvecG(iqi2:fqi2)  *  &
                                para_Op%OpGrid(iterm)%Grid(:,i1_bi,i2_bi)
                END IF
              END DO
              END DO

            END IF

          END DO

        ELSE
          !---- initialization -------------------------------------
          nullify(thread_Psi)
          CALL alloc_array(thread_Psi,[nb_thread],                    &
                          "thread_Psi",name_sub)
          nullify(thread_OpPsi)
          CALL alloc_array(thread_OpPsi,[nb_thread],                  &
                          "thread_OpPsi",name_sub)

          DO ith=1,nb_thread
            thread_Psi(ith) = Psi
            thread_Psi(ith)%GridRep=.TRUE.
            CALL alloc_psi(thread_Psi(ith))

            thread_OpPsi(ith) = Psi
            thread_OpPsi(ith)%GridRep=.TRUE.
            CALL alloc_psi(thread_OpPsi(ith))
          END DO
          !-------------------------------------------------------------
          !-------------------------------------------------------------

          !$OMP parallel do default(none) &
          !$OMP shared(para_Op,thread_Psi,thread_OpPsi,OpPsi,derOp,pot_only) &
          !$OMP private(iterm,ith) &
          !$OMP num_threads(nb_thread)
          DO iterm=1,para_Op%nb_term
            IF ( skip_term(derOp,para_Op%derive_termQdyn(:,iterm)) ) CYCLE
            IF (pot_only .AND. iterm /= para_Op%derive_term_TO_iterm(0,0)) CYCLE


            ith = 1
            !$ ith = omp_get_thread_num()+1

            ! in this subroutine Psi,OpPsi must be allocated with GridRep=.TRUE.
            CALL sub_itermOpPsi_GridRep(thread_Psi(ith),thread_OpPsi(ith),iterm,para_Op)

            !$OMP critical(CRIT_sub_OpPsi)
            IF (OpPsi%cplx) THEN
              OpPsi%CvecG(:) = OpPsi%CvecG(:) + thread_OpPsi(ith)%CvecG(:)
            ELSE
              OpPsi%RvecG(:) = OpPsi%RvecG(:) + thread_OpPsi(ith)%RvecG(:)
            END IF
            !$OMP end critical(CRIT_sub_OpPsi)
          END DO
          !$OMP end parallel do

          !---- deallocation -------------------------------------
          CALL dealloc_array(thread_Psi,  "thread_Psi",  name_sub)
          CALL dealloc_array(thread_OpPsi,"thread_OpPsi",name_sub)

        END IF

        IF (para_Op%cplx) THEN
        IF ( .NOT. skip_term(derOp, [0,0] ) ) THEN
        IF (.NOT. para_Op%imOpGrid(1)%grid_zero) THEN
          CALL sub_PsiBasisRep_TO_GridRep(Psi)

          IF (para_Op%imOpGrid(1)%grid_cte) THEN
            DO i1_bi=1,para_Op%nb_bie
            DO i2_bi=1,para_Op%nb_bie
              iqi1 = 1         + (i1_bi-1) * Psi%nb_qa
              fqi1 = Psi%nb_qa + (i1_bi-1) * Psi%nb_qa
              iqi2 = 1         + (i2_bi-1) * Psi%nb_qa
              fqi2 = Psi%nb_qa + (i2_bi-1) * Psi%nb_qa

              OpPsi%CvecG(iqi1:fqi1) = OpPsi%CvecG(iqi1:fqi1) +       &
                                     Psi%CvecG(iqi2:fqi2) *           &
                          EYE * para_Op%imOpGrid(1)%Mat_cte(i1_bi,i2_bi)

            END DO
            END DO
          ELSE

            DO i1_bi=1,para_Op%nb_bie
            DO i2_bi=1,para_Op%nb_bie
              iqi1 = 1         + (i1_bi-1) * Psi%nb_qa
              fqi1 = Psi%nb_qa + (i1_bi-1) * Psi%nb_qa
              iqi2 = 1         + (i2_bi-1) * Psi%nb_qa
              fqi2 = Psi%nb_qa + (i2_bi-1) * Psi%nb_qa

              OpPsi%CvecG(iqi1:fqi1) = OpPsi%CvecG(iqi1:fqi1) +         &
                                           Psi%CvecG(iqi2:fqi2) * EYE * &
                                  para_Op%imOpGrid(1)%Grid(:,i1_bi,i2_bi)
            END DO
            END DO

          END IF

        END IF
        END IF
        END IF

      ELSE
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'Save_MemGrid_done=.FALSE. is not possible in this subroutine'
        write(out_unitp,*) 'Check the fortran!'
        STOP
      END IF


!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'OpPsiGridRep'
        CALL ecri_psi(Psi=OpPsi)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_OpPsi_WITH_MemGrid

      SUBROUTINE sub_OpPsi_WITH_MemGrid_BGG(Psi,OpPsi,para_Op,derOp,With_Grid,pot_only)
      USE mod_system
      USE mod_basis_BtoG_GtoB, ONLY : DerivOp_TO_CVecG,DerivOp_TO_RVecG
      USE mod_psi,             ONLY : param_psi,ecri_psi,ecri_init_psi,sub_PsiBasisRep_TO_GridRep
      USE mod_SetOp,           ONLY : param_Op,write_param_Op
      IMPLICIT NONE

      !----- variables pour la namelist minimum ------------------------
      TYPE (param_Op)     :: para_Op
      integer, intent(in) :: derOp(2)
      logical, intent(in) :: With_Grid
      logical, intent(in) :: pot_only

      integer          :: n

      integer          :: iOp

      !----- variables for the WP --------------------------------------
      TYPE (param_psi)   :: Psi,OpPsi
      integer :: OpPsi_symab
      real (kind=Rkind), allocatable       :: RG1(:)
      complex (kind=Rkind), allocatable    :: CG1(:)


      !----- working variables -----------------------------------------
      integer :: i1_bi,i2_bi
      integer :: i_qa,iqi1,fqi1,iqi2,fqi2
      integer :: i,ki,k

      integer                  :: nio,error,iterm


      !----- for debuging ----------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_OpPsi_WITH_MemGrid_BGG'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      !-----------------------------------------------------------------
      n = para_Op%nb_tot
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub,' ',n
        write(out_unitp,*) 'nb_bie,nb_baie',para_Op%nb_bie,para_Op%nb_baie
        write(out_unitp,*) 'nb_act1',para_Op%mole%nb_act1
        write(out_unitp,*) 'nb_var',para_Op%mole%nb_var
        write(out_unitp,*) 'para_Op...%Save_MemGrid_done',para_Op%para_ReadOp%para_FileGrid%Save_MemGrid_done
        !flush(out_unitp)
        !CALL write_param_Op(para_Op)
        write(out_unitp,*)
        write(out_unitp,*) 'PsiBasisRep'
        CALL ecri_init_psi(Psi=Psi)
        !CALL ecri_psi(Psi=Psi)
        flush(out_unitp)
      END IF
      !-----------------------------------------------------------------

      IF (para_Op%para_ReadOp%para_FileGrid%Save_MemGrid_done) THEN

        IF (.NOT. With_Grid) THEN
          nb_mult_OpPsi = 0
          CALL sub_PsiBasisRep_TO_GridRep(Psi)
          IF (debug) write(out_unitp,*) 'PsiGridRep'
          IF (debug) CALL ecri_psi(Psi=Psi)

          nb_mult_OpPsi = nb_mult_OpPsi + nb_mult_BTOG
          IF (debug) THEN
            write(out_unitp,*) 'PsiGridRep done'
            write(out_unitp,*) 'PsiBasisRep'
            !CALL ecri_init_psi(Psi=Psi)
            CALL ecri_psi(Psi=Psi)
            flush(out_unitp)
          END IF
        END IF

        CALL sub_sqRhoOVERJac_Psi(Psi,para_Op,inv=.FALSE.)

        IF (Psi%cplx) THEN
          CALL alloc_NParray(CG1,[Psi%nb_qa],"CG1",name_sub)
        ELSE
          CALL alloc_NParray(RG1,[Psi%nb_qa],"RG1",name_sub)
        END IF

        DO iterm=1,para_Op%nb_term

          IF (para_Op%OpGrid(iterm)%grid_zero) CYCLE
          IF (pot_only .AND. iterm /= para_Op%derive_term_TO_iterm(0,0)) CYCLE
          IF ( skip_term(derOp,para_Op%derive_termQdyn(:,iterm)) ) CYCLE


          IF (para_Op%OpGrid(iterm)%grid_cte) THEN

            DO i1_bi=1,para_Op%nb_bie
            DO i2_bi=1,para_Op%nb_bie
              iqi1 = 1         + (i1_bi-1) * Psi%nb_qa
              fqi1 = Psi%nb_qa + (i1_bi-1) * Psi%nb_qa
              iqi2 = 1         + (i2_bi-1) * Psi%nb_qa
              fqi2 = Psi%nb_qa + (i2_bi-1) * Psi%nb_qa

              IF (Psi%cplx) THEN
                CG1 = Psi%CvecG(iqi2:fqi2)
                CALL DerivOp_TO_CVecG(CG1,Psi%nb_qa,para_Op%BasisnD,  &
                                      para_Op%derive_termQdyn(:,iterm))

                OpPsi%CvecG(iqi1:fqi1) = OpPsi%CvecG(iqi1:fqi1) +     &
                      CG1 * para_Op%OpGrid(iterm)%Mat_cte(i1_bi,i2_bi)
              ELSE
                RG1 = Psi%RvecG(iqi2:fqi2)

                !write(6,*) iterm,'RG1>',RG1

                CALL DerivOp_TO_RVecG(RG1,Psi%nb_qa,para_Op%BasisnD,  &
                                      para_Op%derive_termQdyn(:,iterm))

                !write(6,*) iterm,'Op(iterm)RG1>',RG1*para_Op%OpGrid(iterm)%Mat_cte(i1_bi,i2_bi)
                OpPsi%RvecG(iqi1:fqi1) = OpPsi%RvecG(iqi1:fqi1) +     &
                      RG1 * para_Op%OpGrid(iterm)%Mat_cte(i1_bi,i2_bi)

              END IF
            END DO
            END DO
          ELSE

            !write(out_unitp,*) 'iterm,Grid',iterm,Grid
            DO i1_bi=1,para_Op%nb_bie
            DO i2_bi=1,para_Op%nb_bie
              iqi1 = 1         + (i1_bi-1) * Psi%nb_qa
              fqi1 = Psi%nb_qa + (i1_bi-1) * Psi%nb_qa
              iqi2 = 1         + (i2_bi-1) * Psi%nb_qa
              fqi2 = Psi%nb_qa + (i2_bi-1) * Psi%nb_qa

              IF (Psi%cplx) THEN
                CG1 = Psi%CvecG(iqi2:fqi2)
                CALL DerivOp_TO_CVecG(CG1,Psi%nb_qa,para_Op%BasisnD,  &
                                      para_Op%derive_termQdyn(:,iterm))

                OpPsi%CvecG(iqi1:fqi1) = OpPsi%CvecG(iqi1:fqi1) +     &
                       CG1 * para_Op%OpGrid(iterm)%Grid(:,i1_bi,i2_bi)
              ELSE
                RG1 = Psi%RvecG(iqi2:fqi2)

                CALL DerivOp_TO_RVecG(RG1,Psi%nb_qa,para_Op%BasisnD,  &
                                      para_Op%derive_termQdyn(:,iterm))

                OpPsi%RvecG(iqi1:fqi1) = OpPsi%RvecG(iqi1:fqi1) +     &
                         RG1* para_Op%OpGrid(iterm)%Grid(:,i1_bi,i2_bi)

              END IF
            END DO
            END DO

          END IF

        END DO

        IF (para_Op%cplx) THEN
        IF (.NOT. para_Op%imOpGrid(1)%grid_zero) THEN

          IF (para_Op%imOpGrid(1)%grid_cte) THEN
            DO i1_bi=1,para_Op%nb_bie
            DO i2_bi=1,para_Op%nb_bie
              iqi1 = 1         + (i1_bi-1) * Psi%nb_qa
              fqi1 = Psi%nb_qa + (i1_bi-1) * Psi%nb_qa
              iqi2 = 1         + (i2_bi-1) * Psi%nb_qa
              fqi2 = Psi%nb_qa + (i2_bi-1) * Psi%nb_qa

              OpPsi%CvecG(iqi1:fqi1) = OpPsi%CvecG(iqi1:fqi1) +       &
                                     Psi%CvecG(iqi2:fqi2) *           &
                          EYE * para_Op%imOpGrid(1)%Mat_cte(i1_bi,i2_bi)

            END DO
            END DO
          ELSE

            DO i1_bi=1,para_Op%nb_bie
            DO i2_bi=1,para_Op%nb_bie
              iqi1 = 1         + (i1_bi-1) * Psi%nb_qa
              fqi1 = Psi%nb_qa + (i1_bi-1) * Psi%nb_qa
              iqi2 = 1         + (i2_bi-1) * Psi%nb_qa
              fqi2 = Psi%nb_qa + (i2_bi-1) * Psi%nb_qa

              OpPsi%CvecG(iqi1:fqi1) = OpPsi%CvecG(iqi1:fqi1) +         &
                                           Psi%CvecG(iqi2:fqi2) * EYE * &
                                  para_Op%imOpGrid(1)%Grid(:,i1_bi,i2_bi)
            END DO
            END DO

          END IF

        END IF
        END IF

        IF (Psi%cplx) THEN
          CALL dealloc_NParray(CG1,"CG1",name_sub)
        ELSE
          CALL dealloc_NParray(RG1,"RG1",name_sub)
        END IF


      ELSE
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'Save_MemGrid_done=.FALSE. is not possible in this subroutine'
        write(out_unitp,*) 'Check the fortran!'
        STOP
      END IF

      CALL sub_sqRhoOVERJac_Psi(Psi,para_Op,inv=.TRUE.)
      CALL sub_sqRhoOVERJac_Psi(OpPsi,para_Op,inv=.TRUE.)

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'OpPsiGridRep'
        CALL ecri_init_psi(Psi=OpPsi)
        !CALL ecri_psi(Psi=OpPsi)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------
      END SUBROUTINE sub_OpPsi_WITH_MemGrid_BGG



      SUBROUTINE sub_OpPsi_WITH_MemGrid_BGG_Hamil10(Psi,OpPsi,para_Op,derOp,With_Grid,pot_only)
      USE mod_system
      USE mod_Coord_KEO,       ONLY : get_d0GG
      USE mod_basis,           ONLY : rec_Qact
      USE mod_basis_BtoG_GtoB, ONLY : DerivOp_TO_RVecG

      USE mod_psi,             ONLY : param_psi,ecri_psi,sub_PsiBasisRep_TO_GridRep
      USE mod_SetOp,           ONLY : param_Op,write_param_Op
      IMPLICIT NONE


      !----- variables pour la namelist minimum ------------------------
      TYPE (param_Op)     :: para_Op
      integer, intent(in) :: derOp(2)
      logical, intent(in) :: With_Grid
      logical, intent(in) :: pot_only

      integer          :: n

      integer          :: iOp

      !----- variables for the WP --------------------------------------
      TYPE (param_psi)   :: Psi,OpPsi
      integer :: OpPsi_symab
      real (kind=Rkind), allocatable       :: RG1(:)
      real (kind=Rkind), allocatable       :: derRGi(:,:),derRGj(:,:)

      complex (kind=Rkind), allocatable    :: CG1(:)


      !----- working variables -----------------------------------------
      integer :: i1_bi,i2_bi
      integer :: i_qa,iqi1,fqi1,iqi2,fqi2
      integer :: i,j,ki,k,iq,iblock,nb_block,block_size,iq1,iq2,nb_thread

      integer                  :: nio,error,iterm,derive_termQdyn(2)

      real (kind=Rkind), allocatable       :: Qact(:)
      real (kind=Rkind), allocatable       :: GGiq(:,:,:)
      !logical :: direct_KEO = .FALSE.

      !----- for debuging ----------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_OpPsi_WITH_MemGrid_BGG_Hamil10'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      !-----------------------------------------------------------------
      n = para_Op%nb_tot
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub,' ',n
        write(out_unitp,*) 'nb_bie,nb_baie',para_Op%nb_bie,para_Op%nb_baie
        write(out_unitp,*) 'nb_act1',para_Op%mole%nb_act1
        write(out_unitp,*) 'nb_var',para_Op%mole%nb_var
        write(out_unitp,*) 'para_Op...%Save_MemGrid_done',para_Op%para_ReadOp%para_FileGrid%Save_MemGrid_done
        !flush(out_unitp)
        !CALL write_param_Op(para_Op)
        write(out_unitp,*)
        write(out_unitp,*) 'PsiBasisRep'
        !CALL ecri_psi(Psi=Psi)
        flush(out_unitp)
      END IF
      !-----------------------------------------------------------------

     IF (BasisTOGrid_omp == 0) THEN
       nb_thread = 1
     ELSE
       nb_thread = Grid_maxth
     END IF

     IF (pot_only) STOP 'pot_only is not possible in sub_OpPsi_WITH_MemGrid_BGG_Hamil10'

      IF (para_Op%direct_KEO) THEN
        block_size = 1500
        nb_block   = Psi%nb_qa/block_size
        IF (mod(Psi%nb_qa,block_size) > 0) nb_block = nb_block + 1

        CALL alloc_NParray(GGiq,[block_size,para_Op%mole%nb_act,para_Op%mole%nb_act],&
                          'GGiq',name_sub)
      ELSE
        block_size = Psi%nb_qa
        nb_block   = 1
      END IF

      IF (para_Op%para_ReadOp%para_FileGrid%Save_MemGrid_done) THEN

        IF (.NOT. With_Grid) THEN
          nb_mult_OpPsi = 0
          CALL sub_PsiBasisRep_TO_GridRep(Psi)
          nb_mult_OpPsi = nb_mult_OpPsi + nb_mult_BTOG
          IF (debug) THEN
            write(out_unitp,*) 'PsiGridRep done'
            write(out_unitp,*) 'PsiBasisRep'
            !CALL ecri_psi(Psi=Psi)
            flush(out_unitp)
          END IF
        END IF

        CALL sub_sqRhoOVERJac_Psi(Psi,para_Op,inv=.FALSE.)

        IF (Psi%cplx) THEN
          STOP 'cplx in sub_OpPsi_WITH_MemGrid_BGG_Hamil10'
          CALL alloc_NParray(CG1,[Psi%nb_qa],"CG1",name_sub)
          OpPsi%CvecG = CZERO
        ELSE
          CALL alloc_NParray(derRGi,[Psi%nb_qa,para_Op%nb_Qact],"derRGi",name_sub)
          CALL alloc_NParray(derRGj,[Psi%nb_qa,para_Op%nb_Qact],"derRGj",name_sub)
          OpPsi%RvecG = ZERO
        END IF

        IF (para_Op%nb_bie /= 1) STOP 'nb_bie /= 1 in sub_OpPsi_WITH_MemGrid_BGG_Hamil10'

        DO i1_bi=1,para_Op%nb_bie
        DO i2_bi=1,para_Op%nb_bie
          iqi1 = 1         + (i1_bi-1) * Psi%nb_qa
          fqi1 = Psi%nb_qa + (i1_bi-1) * Psi%nb_qa
          iqi2 = 1         + (i2_bi-1) * Psi%nb_qa
          fqi2 = Psi%nb_qa + (i2_bi-1) * Psi%nb_qa

  IF (Psi%cplx) THEN
STOP 'cplx in sub_OpPsi_WITH_MemGrid_BGG_Hamil10'
  ELSE

     ! first derivatives of sqrt(J/rho)*psi,  in derRGi(:,i)
     DO i=1,para_Op%nb_Qact
       derive_termQdyn(:) = [para_Op%mole%liste_QactTOQdyn(i),0]
       derRGi(:,i) = Psi%RvecG(iqi2:fqi2)

       CALL DerivOp_TO_RVecG(derRGi(:,i),Psi%nb_qa,para_Op%BasisnD,  &
                             derive_termQdyn)
     END DO
     derRGj(:,:) = ZERO

     DO iblock=1,nb_block
        iq1 = (iblock-1)*block_size + 1
        iq2 = min(iblock*block_size,Psi%nb_qa)
        IF (para_Op%direct_KEO) THEN

          !$OMP parallel default(none)                 &
          !$OMP shared(para_Op,GGiq,iq1,iq2)           &
          !$OMP private(iq,Qact)                       &
          !$OMP num_threads(nb_thread)
          CALL alloc_NParray(Qact,[para_Op%mole%nb_var],'Qact',name_sub)
          !$OMP  do
          DO iq=iq1,iq2
            CALL Rec_Qact(Qact,para_Op%para_AllBasis%BasisnD,iq,para_Op%mole)
            CALL get_d0GG(Qact,para_Op%para_Tnum,para_Op%mole,d0GG=GGiq(iq-iq1+1,:,:),def=.TRUE.)
          END DO
          !$OMP end do
          CALL dealloc_NParray(Qact,'Qact',name_sub)
          !$OMP end parallel

        END IF

        DO j=1,para_Op%nb_Qact

          ! sum over i, GG(j,i) derRGi(:,i) => derRGj(:,j)
          DO i=1,para_Op%nb_Qact

            IF (para_Op%direct_KEO) THEN
              derRGj(iq1:iq2,j) = derRGj(iq1:iq2,j) + derRGi(iq1:iq2,i) * GGiq(1:iq2-iq1+1,i,j)
            ELSE
              iterm = para_Op%derive_term_TO_iterm(j,i)
              IF (para_Op%OpGrid(iterm)%grid_zero) CYCLE
              IF ( skip_term(derOp,para_Op%derive_termQdyn(:,iterm)) ) CYCLE

              IF (para_Op%OpGrid(iterm)%grid_cte) THEN
                derRGj(iq1:iq2,j) = derRGj(iq1:iq2,j) + derRGi(iq1:iq2,i) * para_Op%OpGrid(iterm)%Mat_cte(i1_bi,i2_bi)
              ELSE
                derRGj(iq1:iq2,j) = derRGj(iq1:iq2,j) + derRGi(iq1:iq2,i) * para_Op%OpGrid(iterm)%Grid(iq1:iq2,i1_bi,i2_bi)
              END IF
            END IF

          END DO
        END DO
     END DO
     CALL dealloc_NParray(derRGi,"derRGi",name_sub)

     DO j=1,para_Op%nb_Qact

       ! multiply by Jac
       IF (allocated(para_Op%para_AllBasis%basis_ext%Jac)) THEN
        derRGj(:,j) = derRGj(:,j) * para_Op%para_AllBasis%basis_ext%Jac
       END IF

       ! derivative with respect to Qact_j
       derive_termQdyn(:) = [para_Op%mole%liste_QactTOQdyn(j),0]
       CALL DerivOp_TO_RVecG(derRGj(:,j),Psi%nb_qa,para_Op%BasisnD,derive_termQdyn)

       ! add each term to OpPsi%RvecG
       OpPsi%RvecG(iqi1:fqi1) = OpPsi%RvecG(iqi1:fqi1) + derRGj(:,j)

     END DO
     CALL dealloc_NParray(derRGj,"derRGj",name_sub)

     IF (allocated(para_Op%para_AllBasis%basis_ext%Jac)) THEN
      OpPsi%RvecG(iqi1:fqi1) = -HALF * OpPsi%RvecG(iqi1:fqi1) / para_Op%para_AllBasis%basis_ext%Jac
     ELSE
      OpPsi%RvecG(iqi1:fqi1) = -HALF * OpPsi%RvecG(iqi1:fqi1)
    END IF

     ! add the potential
     iterm = para_Op%derive_term_TO_iterm(0,0)
     IF (.NOT. para_Op%OpGrid(iterm)%grid_zero .AND.            &
         .NOT. skip_term(derOp,para_Op%derive_termQdyn(:,iterm)) ) THEN

       IF (para_Op%OpGrid(iterm)%grid_cte) THEN
         OpPsi%RvecG(iqi1:fqi1) = OpPsi%RvecG(iqi1:fqi1) +      &
             Psi%RvecG(iqi2:fqi2) * para_Op%OpGrid(iterm)%Mat_cte(i1_bi,i2_bi)
       ELSE
         OpPsi%RvecG(iqi1:fqi1) = OpPsi%RvecG(iqi1:fqi1) +      &
                Psi%RvecG(iqi2:fqi2) * para_Op%OpGrid(iterm)%Grid(:,i1_bi,i2_bi)
       END IF
     END IF

  END IF

           IF (para_Op%cplx) THEN
           stop 'Op  cplx'
!            IF (.NOT. para_Op%imOpGrid(1)%grid_zero) THEN
!
!            IF (para_Op%imOpGrid(1)%grid_cte) THEN
!
!                OpPsi%CvecG(iqi1:fqi1) = OpPsi%CvecG(iqi1:fqi1) +       &
!                                       Psi%CvecG(iqi2:fqi2) *           &
!                            EYE * para_Op%imOpGrid(1)%Mat_cte(i1_bi,i2_bi)
!
!            ELSE
!
!                OpPsi%CvecG(iqi1:fqi1) = OpPsi%CvecG(iqi1:fqi1) +         &
!                                             Psi%CvecG(iqi2:fqi2) * EYE * &
!                                    para_Op%imOpGrid(1)%Grid(:,i1_bi,i2_bi)
!
!           END IF
!           END IF
          END IF

        END DO
        END DO

        IF (Psi%cplx) THEN
          CALL dealloc_NParray(CG1,"CG1",name_sub)
        END IF


      ELSE
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'Save_MemGrid_done=.FALSE. is not possible in this subroutine'
        write(out_unitp,*) 'Check the fortran!'
        STOP
      END IF

      CALL sub_sqRhoOVERJac_Psi(Psi,para_Op,inv=.TRUE.)
      CALL sub_sqRhoOVERJac_Psi(OpPsi,para_Op,inv=.TRUE.)

      IF (allocated(GGiq)) CALL dealloc_NParray(GGiq,'GGiq',name_sub)
!CALL UnCheck_mem()

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'OpPsiGridRep'
        CALL ecri_psi(Psi=OpPsi)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_OpPsi_WITH_MemGrid_BGG_Hamil10

 !SUBROUTINE sub_TabOpPsi_WITH_MemGrid_BGG_Hamil10(Psi,OpPsi,para_Op,derOp,With_Grid,pot_only)
 SUBROUTINE sub_TabOpPsi_WITH_MemGrid_BGG_Hamil10(Psi,OpPsi,para_Op,derOp,With_Grid)
 USE mod_system
 USE mod_Coord_KEO,               ONLY : get_d0GG

 USE mod_basis,                   ONLY : rec_Qact
 USE mod_basis_BtoG_GtoB,         ONLY : DerivOp_TO_RVecG
 USE mod_basis_BtoG_GtoB_SGType4, ONLY : Type_SmolyakRep,Write_SmolyakRep,   &
                                         alloc_SmolyakRep,dealloc_SmolyakRep,&
                                         tabR2bis_TO_SmolyakRep1

 USE mod_psi,                     ONLY : param_psi,ecri_psi,alloc_psi,dealloc_psi,&
                                         sub_PsiBasisRep_TO_GridRep

 USE mod_SetOp,                   ONLY : param_Op,write_param_Op

 IMPLICIT NONE

 !----- variables pour la namelist minimum ------------------------
 TYPE (param_Op)     :: para_Op
 integer, intent(in) :: derOp(2)
 logical, intent(in) :: With_Grid
 !logical, intent(in) :: pot_only

 integer          :: n

 integer          :: iOp

 !----- variables for the WP --------------------------------------
 TYPE (param_psi)   :: Psi(:),OpPsi(:)
 integer :: OpPsi_symab
 real (kind=Rkind), allocatable       :: RG1(:)
 real (kind=Rkind), allocatable       :: derRGi(:,:,:),derRGj(:,:,:)
 TYPE(Type_SmolyakRep)                :: SRep ! temp variable

 complex (kind=Rkind), allocatable    :: CG1(:)


 !----- working variables -----------------------------------------
 integer :: i1_bi,i2_bi
 integer :: i_qa,iqi1,fqi1,iqi2,fqi2
 integer :: itab,i,j,ki,k,iq,iblock,nb_block,block_size,iq1,iq2,nb_thread

 integer                  :: nio,error,iterm,derive_termQdyn(2)

 real (kind=Rkind), allocatable       :: Qact(:)
 real (kind=Rkind), allocatable       :: GGiq(:,:,:)
 !logical :: direct_KEO = .FALSE.

 !----- for debuging ----------------------------------------------
 integer :: err_mem,memory
 character (len=*), parameter :: name_sub='sub_TabOpPsi_WITH_MemGrid_BGG_Hamil10'
 logical, parameter :: debug = .FALSE.
 !logical, parameter :: debug = .TRUE.
 !-----------------------------------------------------------------
 n = para_Op%nb_tot
 IF (debug) THEN
   write(out_unitp,*) 'BEGINNING ',name_sub,' ',n
   write(out_unitp,*) 'nb_bie,nb_baie',para_Op%nb_bie,para_Op%nb_baie
   write(out_unitp,*) 'nb_act1',para_Op%mole%nb_act1
   write(out_unitp,*) 'nb_var',para_Op%mole%nb_var
   write(out_unitp,*) 'para_Op...%Save_MemGrid_done',para_Op%para_ReadOp%para_FileGrid%Save_MemGrid_done
   !flush(out_unitp)
   !CALL write_param_Op(para_Op)
   write(out_unitp,*)
   write(out_unitp,*) 'PsiBasisRep'
   DO itab=1,size(Psi)
     CALL ecri_psi(Psi=Psi(itab))
   END DO
   flush(out_unitp)
 END IF
 !-----------------------------------------------------------------
 IF (Psi(1)%cplx) STOP 'cplx in sub_TabOpPsi_WITH_MemGrid_BGG_Hamil10'
 IF (para_Op%nb_bie /= 1) STOP 'nb_bie /= 1 in sub_TabOpPsi_WITH_MemGrid_BGG_Hamil10'

 IF (BasisTOGrid_omp == 0) THEN
   nb_thread = 1
 ELSE
   nb_thread = Grid_maxth
 END IF

 IF (para_Op%direct_KEO) THEN
   block_size = 10000
   nb_block   = Psi(1)%nb_qa/block_size
   IF (mod(Psi(1)%nb_qa,block_size) > 0) nb_block = nb_block + 1
 ELSE
   block_size = Psi(1)%nb_qa
   nb_block   = 1
 END IF

 CALL alloc_NParray(derRGi,[Psi(1)%nb_qa,para_Op%nb_Qact,size(Psi)],"derRGi",name_sub)

 nb_mult_OpPsi = 0
 DO itab=1,size(Psi)

   IF (.NOT. With_Grid) THEN
     Psi(itab)%GridRep=.TRUE.
     CALL alloc_psi(Psi(itab))
     OpPsi(itab)%GridRep=.TRUE.
     CALL alloc_psi(OpPsi(itab))

     CALL sub_PsiBasisRep_TO_GridRep(Psi(itab))
     nb_mult_OpPsi = nb_mult_OpPsi + nb_mult_BTOG
   END IF

   IF (debug) THEN
     write(out_unitp,*) 'PsiGridRep done'
     write(out_unitp,*) 'PsiBasisRep'
     CALL ecri_psi(Psi=Psi(itab))
     flush(out_unitp)
   END IF

   !CALL ecri_psi(Psi=Psi(itab),ecri_GridRep=.TRUE.)

   CALL sub_sqRhoOVERJac_Psi(Psi(itab),para_Op,inv=.FALSE.)

 END DO

 IF (.NOT. allocated(para_Op%para_AllBasis%basis_ext%Jac)) THEN
   write(out_unitp,*) ' ERROR in ',name_sub
   write(out_unitp,*) ' ....%basis_ext%Jac(:) is not allocated '
   write(out_unitp,*) ' Set JacSave = .TRUE., around line 186 of sub_HSH_harm.f90.'
   STOP
 END IF

 ! first the potential (-2.Jac.V.Psi^G)
 iterm = para_Op%derive_term_TO_iterm(0,0)
 IF (.NOT. para_Op%OpGrid(iterm)%grid_zero .AND.                    &
     .NOT. skip_term(derOp,para_Op%derive_termQdyn(:,iterm)) ) THEN

   IF (para_Op%OpGrid(iterm)%grid_cte) THEN
     DO itab=1,size(Psi)
       OpPsi(itab)%RvecG(:) = -TWO*para_Op%para_AllBasis%basis_ext%Jac * &
             para_Op%OpGrid(iterm)%Mat_cte(1,1) * Psi(itab)%RvecG(:)
     END DO
   ELSE
     DO itab=1,size(Psi)
       OpPsi(itab)%RvecG(:) = -TWO*para_Op%para_AllBasis%basis_ext%Jac * &
              para_Op%OpGrid(iterm)%Grid(:,1,1) * Psi(itab)%RvecG(:)
     END DO
   END IF
 ELSE
   OpPsi(itab)%RvecG(:) = ZERO
 END IF

 !write(out_unitp,*) 'sqRhoOVERJac',para_Op%para_AllBasis%basis_ext%sqRhoOVERJac(:)
 !write(out_unitp,*) 'Jac',para_Op%para_AllBasis%basis_ext%Jac(:)
 !write(out_unitp,*) 'V',para_Op%OpGrid(iterm)%Grid(:,1,1)

 !Transfert sqRhoOVERJac, Jac and the potential in Smolyak rep (Grid)
 IF (debug) THEN
   CALL alloc_SmolyakRep(SRep,                               &
    para_Op%para_AllBasis%BasisnD%para_SGType2%nDind_SmolyakRep%Tab_nDval, &
                para_Op%para_AllBasis%BasisnD%tab_basisPrimSG,grid=.TRUE.)
 END IF

 ! then the derivatives of sqrt(J/rho)*psi,  in derRGi(:,i)
 ! after we don't need psi on the grid
 DO itab=1,size(Psi)

   IF (debug) THEN
     write(out_unitp,*) 'Psi * sqRhoOVERJac'
     SRep = Psi(itab)%RvecG
     CALL Write_SmolyakRep(SRep)
   END IF

   DO i=1,para_Op%nb_Qact
     derive_termQdyn(:) = [para_Op%mole%liste_QactTOQdyn(i),0]
     derRGi(:,i,itab) = Psi(itab)%RvecG(:)
     CALL DerivOp_TO_RVecG(derRGi(:,i,itab),Psi(itab)%nb_qa,para_Op%BasisnD,  &
                           derive_termQdyn)
     IF (debug) THEN
       write(out_unitp,*) 'dQi ',i
       CALL tabR2bis_TO_SmolyakRep1(SRep,derRGi(:,i,itab))
       CALL Write_SmolyakRep(SRep)
     END IF

   END DO
   Psi(itab)%GridRep=.FALSE.
   CALL alloc_psi(Psi(itab))
 END DO

 CALL alloc_NParray(derRGj,[Psi(1)%nb_qa,para_Op%nb_Qact,size(Psi)],"derRGj",name_sub)
 derRGj(:,:,:) = ZERO


 CALL alloc_NParray(GGiq,[block_size,para_Op%mole%nb_act,para_Op%mole%nb_act],&
                   'GGiq',name_sub)

 DO iblock=1,nb_block
    iq1 = (iblock-1)*block_size + 1
    iq2 = min(iblock*block_size,Psi(1)%nb_qa)
    IF (para_Op%direct_KEO) THEN

      !$OMP parallel default(none)                 &
      !$OMP shared(para_Op,GGiq,iq1,iq2)  &
      !$OMP private(iq,Qact)                       &
      !$OMP num_threads(nb_thread)
      CALL alloc_NParray(Qact,[para_Op%mole%nb_var],'Qact',name_sub)
      !$OMP  do
      DO iq=iq1,iq2
        CALL Rec_Qact(Qact,para_Op%para_AllBasis%BasisnD,iq,para_Op%mole)
        CALL get_d0GG(Qact,para_Op%para_Tnum,para_Op%mole,d0GG=GGiq(iq-iq1+1,:,:),def=.TRUE.)
        !write(out_unitp,*) 'iq,Gij',iq,GGiq(iq-iq1+1,:,:)
      END DO
      !$OMP end do
      CALL dealloc_NParray(Qact,'Qact',name_sub)
      !$OMP end parallel
    END IF

    DO itab=1,size(Psi)
    DO j=1,para_Op%nb_Qact

      ! sum over i, GG(j,i) derRGi(:,i) => derRGj(:,j)
      DO i=1,para_Op%nb_Qact

        derRGj(iq1:iq2,j,itab) = derRGj(iq1:iq2,j,itab) +           &
                      derRGi(iq1:iq2,i,itab) * GGiq(1:iq2-iq1+1,i,j)

      END DO
    END DO
    END DO
 END DO
 CALL dealloc_NParray(derRGi,"derRGi",name_sub)
 CALL dealloc_NParray(GGiq,'GGiq',name_sub)

 DO itab=1,size(Psi)
   DO j=1,para_Op%nb_Qact

     ! multiply by Jac
     derRGj(:,j,itab) = derRGj(:,j,itab) * para_Op%para_AllBasis%basis_ext%Jac

     ! derivative with respect to Qact_j
     derive_termQdyn(:) = [para_Op%mole%liste_QactTOQdyn(j),0]
     CALL DerivOp_TO_RVecG(derRGj(:,j,itab),Psi(itab)%nb_qa,para_Op%BasisnD,derive_termQdyn)

     IF (debug) THEN
       write(out_unitp,*) 'dQj ',j
       CALL tabR2bis_TO_SmolyakRep1(SRep,derRGj(:,j,itab))
       CALL Write_SmolyakRep(SRep)
     END IF

     ! add each term to OpPsi%RvecG
     OpPsi(itab)%RvecG(:) = OpPsi(itab)%RvecG(:) + derRGj(:,j,itab)

   END DO
 END DO
 CALL dealloc_NParray(derRGj,"derRGj",name_sub)


 DO itab=1,size(Psi)
   OpPsi(itab)%RvecG(:) = -HALF * OpPsi(itab)%RvecG(:) / para_Op%para_AllBasis%basis_ext%Jac
   CALL sub_sqRhoOVERJac_Psi(OpPsi(itab),para_Op,inv=.TRUE.)

   IF (debug) THEN
     write(out_unitp,*) 'OpPsi Grid '
     CALL tabR2bis_TO_SmolyakRep1(SRep,OpPsi(itab)%RvecG)
     CALL Write_SmolyakRep(SRep)
   END IF

 END DO

 CALL dealloc_SmolyakRep(SRep)
!CALL UnCheck_mem()

 !-----------------------------------------------------------
 IF (debug) THEN
   write(out_unitp,*) 'OpPsiGridRep'
   DO itab=1,size(Psi)
     CALL ecri_psi(Psi=OpPsi(itab))
   END DO
   write(out_unitp,*)
   write(out_unitp,*) 'END ',name_sub
 END IF
 END SUBROUTINE sub_TabOpPsi_WITH_MemGrid_BGG_Hamil10

      SUBROUTINE sub_OpPsi_WITH_FileGrid_type12_BGG(Psi,OpPsi,para_Op,derOp,With_Grid,pot_only)
      USE mod_system
      USE mod_basis_BtoG_GtoB, ONLY : DerivOp_TO_CVecG,DerivOp_TO_RVecG

      USE mod_psi,             ONLY : param_psi,ecri_psi,alloc_psi,dealloc_psi,&
                                      sub_PsiBasisRep_TO_GridRep

      USE mod_SetOp,           ONLY : param_Op
      USE mod_OpGrid,          ONLY : sub_ReadSeq_Grid_iterm,sub_ReadDir_Grid_iterm
      IMPLICIT NONE

      !----- variables pour la namelist minimum ------------------------
      TYPE (param_Op)     :: para_Op
      integer, intent(in) :: derOp(2)
      logical, intent(in) :: With_Grid
      logical, intent(in) :: pot_only

      integer          :: n

      integer          :: iOp

      !----- variables for the WP --------------------------------------
      TYPE (param_psi)   :: Psi,OpPsi
      integer :: OpPsi_symab
      real (kind=Rkind), allocatable       :: RG1(:)
      complex (kind=Rkind), allocatable    :: CG1(:)

      real (kind=Rkind), allocatable :: Grid(:,:,:) ! grid when Save_Grid_iterm=t


      !----- working variables -----------------------------------------
      integer :: i1_bi,i2_bi
      integer :: i_qa,iqi1,fqi1,iqi2,fqi2
      integer :: i,ki,k

      integer                  :: nio,error,iterm


      !----- for debuging ----------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_OpPsi_WITH_FileGrid_type12_BGG'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      !-----------------------------------------------------------------
      n = para_Op%nb_tot
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub,' ',n
        write(out_unitp,*) 'nb_bie,nb_baie',para_Op%nb_bie,para_Op%nb_baie
        write(out_unitp,*) 'nb_act1',para_Op%mole%nb_act1
        write(out_unitp,*) 'nb_var',para_Op%mole%nb_var
        write(out_unitp,*) 'para_Op...%Save_MemGrid_done',              &
                     para_Op%para_ReadOp%para_FileGrid%Save_MemGrid_done
        !flush(out_unitp)
        !CALL write_param_Op(para_Op)
        write(out_unitp,*)
        write(out_unitp,*) 'PsiBasisRep'
        CALL ecri_psi(Psi=Psi)
        flush(out_unitp)
      END IF
      !-----------------------------------------------------------------

      IF (.NOT. para_Op%para_ReadOp%para_FileGrid%Save_MemGrid_done .AND. &
             (para_Op%para_ReadOp%para_FileGrid%Type_FileGrid == 1 .OR.   &
              para_Op%para_ReadOp%para_FileGrid%Type_FileGrid == 2)) THEN

        CALL alloc_NParray(Grid,                                          &
                       [para_Op%nb_qa,para_Op%nb_bie,para_Op%nb_bie], &
                        'Grid',name_sub)

        IF (.NOT. With_Grid) THEN
          nb_mult_OpPsi = 0
          CALL sub_PsiBasisRep_TO_GridRep(Psi)
          nb_mult_OpPsi = nb_mult_OpPsi + nb_mult_BTOG
          IF (debug) THEN
            write(out_unitp,*) 'PsiGridRep done'
            write(out_unitp,*) 'PsiBasisRep'
            CALL ecri_psi(Psi=Psi)
            flush(out_unitp)
          END IF
        END IF

        CALL sub_sqRhoOVERJac_Psi(Psi,para_Op,inv=.FALSE.)


        IF (Psi%cplx) THEN
          CALL alloc_NParray(CG1,[Psi%nb_qa],"CG1",name_sub)
        ELSE
          CALL alloc_NParray(RG1,[Psi%nb_qa],"RG1",name_sub)
        END IF

        DO iterm=1,para_Op%nb_term

          IF (para_Op%OpGrid(iterm)%grid_zero) CYCLE
          IF (pot_only .AND. iterm /= para_Op%derive_term_TO_iterm(0,0)) CYCLE
          IF ( skip_term(derOp,para_Op%derive_termQdyn(:,iterm)) ) CYCLE


          IF (para_Op%OpGrid(iterm)%grid_cte) THEN

            DO i1_bi=1,para_Op%nb_bie
            DO i2_bi=1,para_Op%nb_bie
              iqi1 = 1         + (i1_bi-1) * Psi%nb_qa
              fqi1 = Psi%nb_qa + (i1_bi-1) * Psi%nb_qa
              iqi2 = 1         + (i2_bi-1) * Psi%nb_qa
              fqi2 = Psi%nb_qa + (i2_bi-1) * Psi%nb_qa

              IF (Psi%cplx) THEN
                CG1 = Psi%CvecG(iqi2:fqi2)
                CALL DerivOp_TO_CVecG(CG1,Psi%nb_qa,para_Op%BasisnD,  &
                                      para_Op%derive_termQdyn(:,iterm))

                OpPsi%CvecG(iqi1:fqi1) = OpPsi%CvecG(iqi1:fqi1) +     &
                      CG1 * para_Op%OpGrid(iterm)%Mat_cte(i1_bi,i2_bi)
              ELSE
                RG1 = Psi%RvecG(iqi2:fqi2)
                CALL DerivOp_TO_RVecG(RG1,Psi%nb_qa,para_Op%BasisnD,  &
                                      para_Op%derive_termQdyn(:,iterm))

                OpPsi%RvecG(iqi1:fqi1) = OpPsi%RvecG(iqi1:fqi1) +     &
                      RG1 * para_Op%OpGrid(iterm)%Mat_cte(i1_bi,i2_bi)


              END IF
            END DO
            END DO
          ELSE

            Grid(:,:,:) = ZERO

            !$OMP critical(CRIT2_sub_OpPsi_WITH_FileGrid_type12_BGG)
            IF (para_Op%OpGrid(iterm)%file_Grid%seq) THEN   ! sequential acces file
              CALL sub_ReadSeq_Grid_iterm(Grid,para_Op%OpGrid(iterm))
            ELSE  ! direct acces file
              CALL sub_ReadDir_Grid_iterm(Grid,para_Op%OpGrid(iterm))
            END IF
            !$OMP end critical(CRIT2_sub_OpPsi_WITH_FileGrid_type12_BGG)

            !write(out_unitp,*) 'iterm,Grid',iterm,Grid
            DO i1_bi=1,para_Op%nb_bie
            DO i2_bi=1,para_Op%nb_bie
              iqi1 = 1         + (i1_bi-1) * Psi%nb_qa
              fqi1 = Psi%nb_qa + (i1_bi-1) * Psi%nb_qa
              iqi2 = 1         + (i2_bi-1) * Psi%nb_qa
              fqi2 = Psi%nb_qa + (i2_bi-1) * Psi%nb_qa

              IF (Psi%cplx) THEN
                CG1 = Psi%CvecG(iqi2:fqi2)
                CALL DerivOp_TO_CVecG(CG1,Psi%nb_qa,para_Op%BasisnD,  &
                                      para_Op%derive_termQdyn(:,iterm))

                OpPsi%CvecG(iqi1:fqi1) = OpPsi%CvecG(iqi1:fqi1) +     &
                       CG1 * Grid(:,i1_bi,i2_bi)
              ELSE
                RG1 = Psi%RvecG(iqi2:fqi2)

                CALL DerivOp_TO_RVecG(RG1,Psi%nb_qa,para_Op%BasisnD,  &
                                      para_Op%derive_termQdyn(:,iterm))

                OpPsi%RvecG(iqi1:fqi1) = OpPsi%RvecG(iqi1:fqi1) +     &
                         RG1 * Grid(:,i1_bi,i2_bi)

              END IF
            END DO
            END DO

          END IF

        END DO

        IF (para_Op%cplx) THEN
        IF (.NOT. para_Op%imOpGrid(1)%grid_zero) THEN

          IF (para_Op%imOpGrid(1)%grid_cte) THEN
            DO i1_bi=1,para_Op%nb_bie
            DO i2_bi=1,para_Op%nb_bie
              iqi1 = 1         + (i1_bi-1) * Psi%nb_qa
              fqi1 = Psi%nb_qa + (i1_bi-1) * Psi%nb_qa
              iqi2 = 1         + (i2_bi-1) * Psi%nb_qa
              fqi2 = Psi%nb_qa + (i2_bi-1) * Psi%nb_qa

              OpPsi%CvecG(iqi1:fqi1) = OpPsi%CvecG(iqi1:fqi1) +       &
                                     Psi%CvecG(iqi2:fqi2) *           &
                          EYE * para_Op%imOpGrid(1)%Mat_cte(i1_bi,i2_bi)

            END DO
            END DO
          ELSE

            Grid(:,:,:) = ZERO
            !$OMP critical(CRIT3_sub_OpPsi_WITH_FileGrid_type12_BGG)
            IF (para_Op%imOpGrid(1)%file_Grid%seq) THEN   ! sequential acces file
              CALL sub_ReadSeq_Grid_iterm(Grid,para_Op%imOpGrid(1))
            ELSE  ! direct acces file
              CALL sub_ReadDir_Grid_iterm(Grid,para_Op%imOpGrid(1))
            END IF
            !$OMP end critical(CRIT3_sub_OpPsi_WITH_FileGrid_type12_BGG)

            DO i1_bi=1,para_Op%nb_bie
            DO i2_bi=1,para_Op%nb_bie
              iqi1 = 1         + (i1_bi-1) * Psi%nb_qa
              fqi1 = Psi%nb_qa + (i1_bi-1) * Psi%nb_qa
              iqi2 = 1         + (i2_bi-1) * Psi%nb_qa
              fqi2 = Psi%nb_qa + (i2_bi-1) * Psi%nb_qa

              OpPsi%CvecG(iqi1:fqi1) = OpPsi%CvecG(iqi1:fqi1) +         &
                        Psi%CvecG(iqi2:fqi2) * EYE * Grid(:,i1_bi,i2_bi)
            END DO
            END DO

          END IF

        END IF
        END IF

        IF (Psi%cplx) THEN
          CALL dealloc_NParray(CG1,"CG1",name_sub)
        ELSE
          CALL dealloc_NParray(RG1,"RG1",name_sub)
        END IF

        CALL dealloc_NParray(Grid,'Grid',name_sub)


      ELSE
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'Save_MemGrid_done=.TRUE. is not possible in this subroutine'
        write(out_unitp,*) 'Check the fortran!'
        STOP
      END IF

      CALL sub_sqRhoOVERJac_Psi(  Psi,para_Op,inv=.TRUE.)
      CALL sub_sqRhoOVERJac_Psi(OpPsi,para_Op,inv=.TRUE.)

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'OpPsiGridRep'
        CALL ecri_psi(Psi=OpPsi)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------
      END SUBROUTINE sub_OpPsi_WITH_FileGrid_type12_BGG



      SUBROUTINE sub_OpPsi_WITH_FileGrid_type12(Psi,OpPsi,para_Op,derOp,pot_only)
      USE mod_system
      USE mod_psi,          ONLY : param_psi,ecri_psi,                  &
                                   sub_d0d1d2PsiBasisRep_TO_GridRep,    &
                                   sub_PsiBasisRep_TO_GridRep

      USE mod_SetOp,        ONLY : param_Op,write_param_Op
      USE mod_OpGrid,       ONLY : sub_ReadSeq_Grid_iterm,sub_ReadDir_Grid_iterm
      IMPLICIT NONE

      !----- variables pour la namelist minimum ------------------------
      TYPE (param_Op)     :: para_Op
      integer, intent(in) :: derOp(2)
      logical, intent(in) :: pot_only

      integer          :: n

      !----- variables for the WP --------------------------------------
      TYPE (param_psi)   :: Psi,OpPsi


      !----- working variables -----------------------------------------
      integer :: i1_bi,i2_bi
      integer :: i_qa,iqi1,fqi1,iqi2,fqi2

      real (kind=Rkind), allocatable :: Grid(:,:,:) ! grid when Save_Grid_iterm=t



      integer                  :: nio,error,iterm


      !----- for debuging ----------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_OpPsi_WITH_FileGrid_type12'
      logical, parameter :: debug = .FALSE.
!      logical, parameter :: debug = .TRUE.
      !-----------------------------------------------------------------
      n = para_Op%nb_tot
      !write(out_unitp,*) 'BEGINNING ',name_sub,' ',n
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub,' ',n
        write(out_unitp,*) 'nb_bie,nb_baie',para_Op%nb_bie,para_Op%nb_baie
        write(out_unitp,*) 'nb_act1',para_Op%mole%nb_act1
        write(out_unitp,*) 'nb_var',para_Op%mole%nb_var
        write(out_unitp,*) 'para_Op%mat_done',para_Op%mat_done
        flush(out_unitp)
        CALL write_param_Op(para_Op)
        write(out_unitp,*)
        write(out_unitp,*) 'PsiBasisRep'
        CALL ecri_psi(Psi=Psi)
        write(out_unitp,*) 'ini OpPsiBasisRep'
        CALL ecri_psi(Psi=OpPsi)
        flush(out_unitp)
      END IF
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
!     - test on nb_bi : nb_bi>0
      IF (para_Op%nb_bi <= 0) THEN
        write(out_unitp,*) ' ERROR : nb_bi MUST be > 0'
        write(out_unitp,*) 'nb_bi',para_Op%nb_bi
        write(out_unitp,*) 'STOP in ',name_sub
        STOP
      END IF
      !-----------------------------------------------------------------

      IF (para_Op%para_ReadOp%para_FileGrid%Type_FileGrid /= 0 .AND.    &
              .NOT. para_Op%para_ReadOp%para_FileGrid%Save_MemGrid) THEN
      !=================================================================
      !       with grid on a sequential/direct access file
      !=================================================================

        CALL alloc_NParray(Grid,                                          &
                       [para_Op%nb_qa,para_Op%nb_bie,para_Op%nb_bie], &
                        'Grid',name_sub)

        !-----------------------------------------------------------------
        DO iterm=1,para_Op%nb_term

          IF (para_Op%OpGrid(iterm)%grid_zero) CYCLE
          IF ( skip_term(derOp,para_Op%derive_termQdyn(:,iterm)) ) CYCLE
          IF (pot_only .AND. iterm /= para_Op%derive_term_TO_iterm(0,0)) CYCLE

          IF (para_Op%OpGrid(iterm)%grid_cte) THEN

            CALL sub_d0d1d2PsiBasisRep_TO_GridRep(Psi,para_Op%derive_termQdyn(:,iterm))

            DO i1_bi=1,para_Op%nb_bie
            DO i2_bi=1,para_Op%nb_bie
              iqi1 = 1         + (i1_bi-1) * Psi%nb_qa
              fqi1 = Psi%nb_qa + (i1_bi-1) * Psi%nb_qa
              iqi2 = 1         + (i2_bi-1) * Psi%nb_qa
              fqi2 = Psi%nb_qa + (i2_bi-1) * Psi%nb_qa


              IF (Psi%cplx) THEN
                OpPsi%CvecG(iqi1:fqi1) = OpPsi%CvecG(iqi1:fqi1) +     &
                   Psi%CvecG(iqi2:fqi2) *                             &
                   para_Op%OpGrid(iterm)%Mat_cte(i1_bi,i2_bi)
              ELSE
                OpPsi%RvecG(iqi1:fqi1) = OpPsi%RvecG(iqi1:fqi1) +     &
                   Psi%RvecG(iqi2:fqi2) *                             &
                   para_Op%OpGrid(iterm)%Mat_cte(i1_bi,i2_bi)
              END IF
            END DO
            END DO
          ELSE
            CALL sub_d0d1d2PsiBasisRep_TO_GridRep(Psi,para_Op%derive_termQdyn(:,iterm))

            Grid(:,:,:) = ZERO

            !$OMP critical(CRIT2_sub_OpPsi)
            IF (para_Op%OpGrid(iterm)%file_Grid%seq) THEN   ! sequential acces file
              CALL sub_ReadSeq_Grid_iterm(Grid,para_Op%OpGrid(iterm))
            ELSE  ! direct acces file
              CALL sub_ReadDir_Grid_iterm(Grid,para_Op%OpGrid(iterm))
            END IF
            !$OMP end critical(CRIT2_sub_OpPsi)


            !write(out_unitp,*) 'iterm,Grid',iterm,Grid
            DO i1_bi=1,para_Op%nb_bie
            DO i2_bi=1,para_Op%nb_bie
              iqi1 = 1         + (i1_bi-1) * Psi%nb_qa
              fqi1 = Psi%nb_qa + (i1_bi-1) * Psi%nb_qa
              iqi2 = 1         + (i2_bi-1) * Psi%nb_qa
              fqi2 = Psi%nb_qa + (i2_bi-1) * Psi%nb_qa

              IF (Psi%cplx) THEN
                OpPsi%CvecG(iqi1:fqi1) = OpPsi%CvecG(iqi1:fqi1) +     &
                   Psi%CvecG(iqi2:fqi2) * Grid(:,i1_bi,i2_bi)
              ELSE
                OpPsi%RvecG(iqi1:fqi1) = OpPsi%RvecG(iqi1:fqi1) +     &
                   Psi%RvecG(iqi2:fqi2) * Grid(:,i1_bi,i2_bi)
              END IF
            END DO
            END DO

          END IF

        END DO


        IF (para_Op%cplx ) THEN
        IF ( .NOT. skip_term(derOp, [0,0] ) ) THEN
        IF (.NOT. para_Op%imOpGrid(1)%grid_zero) THEN
          CALL sub_PsiBasisRep_TO_GridRep(Psi)

          IF (para_Op%imOpGrid(1)%grid_cte) THEN
            DO i1_bi=1,para_Op%nb_bie
            DO i2_bi=1,para_Op%nb_bie
              iqi1 = 1         + (i1_bi-1) * Psi%nb_qa
              fqi1 = Psi%nb_qa + (i1_bi-1) * Psi%nb_qa
              iqi2 = 1         + (i2_bi-1) * Psi%nb_qa
              fqi2 = Psi%nb_qa + (i2_bi-1) * Psi%nb_qa

              OpPsi%CvecG(iqi1:fqi1) = OpPsi%CvecG(iqi1:fqi1) +         &
                                     Psi%CvecG(iqi2:fqi2) *             &
                          EYE * para_Op%imOpGrid(1)%Mat_cte(i1_bi,i2_bi)

            END DO
            END DO
          ELSE
            Grid(:,:,:) = ZERO
            !$OMP critical(CRIT3_sub_OpPsi)
            IF (para_Op%imOpGrid(1)%file_Grid%seq) THEN   ! sequential acces file
              CALL sub_ReadSeq_Grid_iterm(Grid,para_Op%imOpGrid(1))
            ELSE  ! direct acces file
              CALL sub_ReadDir_Grid_iterm(Grid,para_Op%imOpGrid(1))
            END IF
            !$OMP end critical(CRIT3_sub_OpPsi)

            DO i1_bi=1,para_Op%nb_bie
            DO i2_bi=1,para_Op%nb_bie
              iqi1 = 1         + (i1_bi-1) * Psi%nb_qa
              fqi1 = Psi%nb_qa + (i1_bi-1) * Psi%nb_qa
              iqi2 = 1         + (i2_bi-1) * Psi%nb_qa
              fqi2 = Psi%nb_qa + (i2_bi-1) * Psi%nb_qa

              OpPsi%CvecG(iqi1:fqi1) = OpPsi%CvecG(iqi1:fqi1) +         &
                        Psi%CvecG(iqi2:fqi2) * EYE * Grid(:,i1_bi,i2_bi)
            END DO
            END DO

          END IF

        END IF
        END IF
        END IF

        CALL dealloc_NParray(Grid,'Grid',name_sub)

      ELSE
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'Save_MemGrid=.TRUE. or Type_FileGrid=0 is not possible in this subroutine'
        write(out_unitp,*) 'Check the fortran!'
        STOP
      END IF
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'OpPsiGridRep'
        CALL ecri_psi(Psi=OpPsi)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------
      !write(out_unitp,*) 'END ',name_sub

      END SUBROUTINE sub_OpPsi_WITH_FileGrid_type12


      SUBROUTINE  sub_OpPsi_WITH_FileGrid_type0(Psi,OpPsi,para_Op,derOp,pot_only)
      USE mod_system
      USE mod_PrimOp,  ONLY : param_d0MatOp,Init_d0MatOp,dealloc_d0MatOp
      USE mod_psi,     ONLY : param_psi,ecri_psi,alloc_psi,dealloc_psi, &
                              alloc_array,dealloc_array,                &
                              sub_d0d1d2PsiBasisRep_TO_GridRep,         &
                              sub_PsiBasisRep_TO_GridRep

      USE mod_SetOp,   ONLY : param_Op,write_param_Op

      IMPLICIT NONE

      !----- variables pour la namelist minimum ------------------------
      TYPE (param_Op)     :: para_Op
      integer, intent(in) :: derOp(2)
      logical, intent(in) :: pot_only

      integer          :: n

      integer          :: iOp

      !----- variables for the WP --------------------------------------
      TYPE (param_psi)   :: Psi,OpPsi
      integer :: OpPsi_symab
      TYPE (param_psi), pointer :: tab_Psi(:)


      !----- working variables -----------------------------------------

      integer :: i1_bi,i2_bi
      integer :: i_qa,iqi1,fqi1,iqi2,fqi2
      integer :: i,ki,k

      real (kind=Rkind) :: WnD
      real (kind=Rkind) :: Qdyn(para_Op%mole%nb_var)
      real (kind=Rkind) :: Qact(para_Op%mole%nb_act1)

      TYPE (param_d0MatOp) :: d0MatOp
      integer              :: type_Op,iterm_Op,id1,id2

      real (kind=Rkind), allocatable :: Grid(:,:,:) ! grid when Save_Grid_iterm=t


      character (len=Line_len) :: name_file,name_file_th
      character (len=Name_len) :: name_term,name_th_num,name_n_Op
      integer                  :: nio,error,iterm
      integer                  :: lrecl_Grid_iterm


      !------ for OpenMP -----------------------------------------------
      TYPE (param_psi), pointer :: thread_Psi(:),thread_OpPsi(:)
      integer       :: nb_thread,ith,ithread
      logical       :: file_is_para


      !----- for debuging ----------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_OpPsi_WITH_FileGrid_type0'
      logical, parameter :: debug = .FALSE.
!      logical, parameter :: debug = .TRUE.
      !-----------------------------------------------------------------
      n = para_Op%nb_tot
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub,' ',n
        write(out_unitp,*) 'nb_bie,nb_baie',para_Op%nb_bie,para_Op%nb_baie
        write(out_unitp,*) 'nb_act1',para_Op%mole%nb_act1
        write(out_unitp,*) 'nb_var',para_Op%mole%nb_var
        write(out_unitp,*) 'para_Op%mat_done',para_Op%mat_done
        flush(out_unitp)
        CALL write_param_Op(para_Op)
        write(out_unitp,*)
        write(out_unitp,*) 'PsiBasisRep'
        CALL ecri_psi(Psi=Psi)
        write(out_unitp,*) 'ini OpPsiBasisRep'
        CALL ecri_psi(Psi=OpPsi)
        flush(out_unitp)
      END IF
      !-----------------------------------------------------------------



       IF (para_Op%para_ReadOp%para_FileGrid%Type_FileGrid == 0 .AND.   &
              .NOT. para_Op%para_ReadOp%para_FileGrid%Save_MemGrid) THEN
      !=====================================================================
      !
      !       with Type_FileGrid=0 grid (normal SH_HADA file)
      !       Read SH_HADA file each time (for huge grid)
      !
      !=====================================================================

        nullify(tab_Psi)
        CALL alloc_array(tab_Psi,[para_Op%nb_term],"tab_Psi",name_sub)
        DO iterm=1,para_Op%nb_term
          tab_Psi(iterm) = Psi
          tab_Psi(iterm)%GridRep=.TRUE.
          CALL alloc_psi(tab_Psi(iterm))
          !- calculation of d0d1d2psi as a function derive_termQact
          CALL sub_d0d1d2PsiBasisRep_TO_GridRep(tab_Psi(iterm),           &
                                     para_Op%derive_termQdyn(:,iterm))

        END DO
        IF (para_Op%cplx) CALL sub_PsiBasisRep_TO_GridRep(Psi)


        IF (para_Op%name_Op == 'H') THEN
          type_Op = para_Op%para_ReadOp%Type_HamilOp ! H
          IF (type_Op /= 1) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) '    Type_HamilOp MUST be equal to 1 for HADA or cHAC'
            write(out_unitp,*) '    CHECK your data!!'
            STOP
          END IF
        ELSE
          type_Op = 0
        END IF

        CALL Init_d0MatOp(d0MatOp,para_Op%param_TypeOp,para_Op%nb_bie)

        DO i_qa=1,para_Op%nb_qa

          CALL sub_reading_Op(i_qa,para_Op%nb_qa,d0MatOp,para_Op%n_Op,  &
                          Qdyn,para_Op%mole%nb_var,para_Op%mole%nb_act1, &
                          Qact,WnD,para_Op%file_grid)

          DO i1_bi=1,para_Op%nb_bie
          DO i2_bi=1,para_Op%nb_bie
            iqi1 = i_qa      + (i1_bi-1) * Psi%nb_qa
            fqi1 = i_qa      + (i1_bi-1) * Psi%nb_qa
            iqi2 = i_qa      + (i2_bi-1) * Psi%nb_qa
            fqi2 = i_qa      + (i2_bi-1) * Psi%nb_qa

            IF (Psi%cplx) THEN
              DO iterm=1,para_Op%nb_term
                IF ( skip_term(derOp,para_Op%derive_termQdyn(:,iterm)) ) CYCLE
                IF (pot_only .AND. iterm /= para_Op%derive_term_TO_iterm(0,0)) CYCLE

                id1 = para_Op%derive_termQact(1,iterm)
                id2 = para_Op%derive_termQact(2,iterm)
                iterm_Op = d0MatOp%derive_term_TO_iterm(id1,id2)

                OpPsi%CvecG(iqi1:fqi1) = OpPsi%CvecG(iqi1:fqi1) +       &
                                      tab_Psi(iterm)%CvecG(iqi2:fqi2) * &
                                      d0MatOp%ReVal(i1_bi,i2_bi,iterm_Op)
              END DO
            ELSE
              DO iterm=1,para_Op%nb_term
                IF ( skip_term(derOp,para_Op%derive_termQdyn(:,iterm)) ) CYCLE
                IF (pot_only .AND. iterm /= para_Op%derive_term_TO_iterm(0,0)) CYCLE

                id1 = para_Op%derive_termQact(1,iterm)
                id2 = para_Op%derive_termQact(2,iterm)
                iterm_Op = d0MatOp%derive_term_TO_iterm(id1,id2)

                OpPsi%RvecG(iqi1:fqi1) = OpPsi%RvecG(iqi1:fqi1) +       &
                                      tab_Psi(iterm)%RvecG(iqi2:fqi2) * &
                                      d0MatOp%ReVal(i1_bi,i2_bi,iterm_Op)
              END DO
            END IF

            IF (para_Op%cplx) THEN
            IF ( skip_term(derOp,[0,0]) ) THEN
               OpPsi%CvecG(iqi1:fqi1) = OpPsi%CvecG(iqi1:fqi1) +      &
                                        Psi%CvecG(iqi2:fqi2) *        &
                                       EYE * d0MatOp%ImVal(i1_bi,i2_bi)
            END IF
            END IF
          END DO
          END DO

        END DO

        CALL dealloc_d0MatOp(d0MatOp)

        DO iterm=1,para_Op%nb_term
          CALL dealloc_psi(tab_Psi(iterm))
        END DO
        CALL dealloc_array(tab_Psi,"tab_Psi","sub_OpPsi")

      ELSE
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'Save_MemGrid=.TRUE. or Type_FileGrid=1,2 is not possible in this subroutine'
        write(out_unitp,*) 'Check the fortran!'
      END IF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'OpPsiGridRep'
        CALL ecri_psi(Psi=OpPsi)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_OpPsi_WITH_FileGrid_type0

      SUBROUTINE sub_itermOpPsi_GridRep(Psi,OpPsi,iterm,para_Op)
      USE mod_system
      USE mod_SetOp,    ONLY : param_Op,write_param_Op
      USE mod_psi,      ONLY : param_psi,ecri_psi,sub_d0d1d2PsiBasisRep_TO_GridRep
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)              :: para_Op
      integer,         intent(in)  :: iterm

!----- variables for the WP ----------------------------------------
      TYPE (param_psi)      :: Psi
      TYPE (param_psi)      :: OpPsi


!----- working variables -------------------------------------------
      integer :: i1_bi,i2_bi
      integer :: i_qa,iqi1,fqi1,iqi2,fqi2


!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING sub_itermOpPsi_GridRep'
        write(out_unitp,*) 'iterm',iterm
        write(out_unitp,*) 'PsiBasisRep'
        CALL ecri_psi(Psi=Psi)
        write(out_unitp,*) 'ini OpPsiBasisRep'
        CALL ecri_psi(Psi=OpPsi)
      END IF
!-----------------------------------------------------------

      IF (para_Op%para_ReadOp%para_FileGrid%Save_MemGrid_done) THEN
!=====================================================================
!
!       with in memory Grid
!
!=====================================================================

!       Psi%GridRep=.TRUE.
!       CALL alloc_psi(Psi)
!       OpPsi%GridRep=.TRUE.
!       CALL alloc_psi(OpPsi)

        OpPsi = ZERO
!-------------------------------------------------------------
!-------------------------------------------------------------
!       - calculation of d0d1d2psi as a function derive_termQact
        IF (para_Op%OpGrid(iterm)%grid_zero) RETURN

        CALL sub_d0d1d2PsiBasisRep_TO_GridRep(Psi,para_Op%derive_termQdyn(:,iterm))

        DO i1_bi=1,para_Op%nb_bie
        DO i2_bi=1,para_Op%nb_bie
          iqi1 = 1         + (i1_bi-1) * Psi%nb_qa
          fqi1 = Psi%nb_qa + (i1_bi-1) * Psi%nb_qa
          iqi2 = 1         + (i2_bi-1) * Psi%nb_qa
          fqi2 = Psi%nb_qa + (i2_bi-1) * Psi%nb_qa


          IF (para_Op%OpGrid(iterm)%grid_cte) THEN
            IF (Psi%cplx) THEN
              OpPsi%CvecG(iqi1:fqi1) = OpPsi%CvecG(iqi1:fqi1) +         &
                 Psi%CvecG(iqi2:fqi2) *                                 &
                   para_Op%OpGrid(iterm)%Mat_cte(i1_bi,i2_bi)
            ELSE
              OpPsi%RvecG(iqi1:fqi1) = OpPsi%RvecG(iqi1:fqi1) +         &
                 Psi%RvecG(iqi2:fqi2) *                                 &
                   para_Op%OpGrid(iterm)%Mat_cte(i1_bi,i2_bi)
            END IF
          ELSE
            IF (Psi%cplx) THEN
              OpPsi%CvecG(iqi1:fqi1) = OpPsi%CvecG(iqi1:fqi1) +         &
                 Psi%CvecG(iqi2:fqi2) *                                 &
                   para_Op%OpGrid(iterm)%Grid(:,i1_bi,i2_bi)
            ELSE
              OpPsi%RvecG(iqi1:fqi1) = OpPsi%RvecG(iqi1:fqi1) +         &
                 Psi%RvecG(iqi2:fqi2) *                                 &
                   para_Op%OpGrid(iterm)%Grid(:,i1_bi,i2_bi)
            END IF
          END IF
        END DO
        END DO

      ELSE ! para_Op%para_ReadOp%para_FileGrid%Save_MemGrid_done = .FLASE.
         write(out_unitp,*) 'ERROR in  sub_itermOpPsi_GridRep'
         write(out_unitp,*) 'Save_MemGrid_done = .FALSE. is not possible'
         STOP
      END IF

      END SUBROUTINE sub_itermOpPsi_GridRep

!=======================================================================================
!>     OpPsi = (OpPsi - E0.Psi) / Esc
!>     Note the symab
!=======================================================================================
      SUBROUTINE sub_scaledOpPsi(Psi,OpPsi,E0,Esc)
      USE mod_system
      USE mod_psi,   ONLY : param_psi,ecri_psi
      IMPLICIT NONE

!----- for the scaling -------------------------------------------
      real (kind=Rkind)  :: E0,Esc

!----- variables for the WP ----------------------------------------
      TYPE (param_psi)   :: OpPsi,Psi


!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'sub_scaledOpPsi'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'E0,Esc',E0,Esc
        write(out_unitp,*) 'Psi'
        CALL ecri_psi(Psi=Psi)
        write(out_unitp,*) 'OpPsi'
        CALL ecri_psi(Psi=OpPsi)
      END IF

      IF(keep_MPI) THEN
        IF (Psi%cplx) THEN
          OpPsi%CvecB(:) = (OpPsi%CvecB(:) - E0*Psi%CvecB(:))/Esc
        ELSE
          OpPsi%RvecB(:) = (OpPsi%RvecB(:) - E0*Psi%RvecB(:))/Esc
        END IF
      ENDIF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'sub_scaledOpPsi'
        CALL ecri_psi(Psi=OpPsi)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_scaledOpPsi
!=======================================================================================

!=======================================================================================
!     | Op(i).Psi> = Op | Psi>
!     just for ONE term (iop) of the op
!
!=======================================================================================
      SUBROUTINE sub_OpiPsi(Psi,OpPsi,para_Op,iOp)
      USE mod_system
      USE mod_SetOp,   ONLY : param_Op,write_param_Op,alloc_para_Op,read_OpGrid_OF_Op
      USE mod_psi,     ONLY : param_psi,ecri_psi,alloc_psi,dealloc_psi, &
                              sub_d0d1d2PsiBasisRep_TO_GridRep,         &
                              sub_PsiGridRep_TO_BasisRep
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)  :: para_Op
      integer          :: n

      integer          :: iOp

!----- variables for the WP ----------------------------------------
      TYPE (param_psi)   :: Psi,OpPsi



!----- working variables -------------------------------------------
      integer :: i1_bi,i2_bi
      integer :: iqi1,fqi1,iqi2,fqi2



!     logical       :: old = .TRUE.
      logical       :: old = .FALSE.

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_OpiPsi'
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
      n = para_Op%nb_tot
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub,' ',n
        write(out_unitp,*) 'nb_bie,nb_baie',para_Op%nb_bie,para_Op%nb_baie
        write(out_unitp,*) 'nb_act1',para_Op%mole%nb_act1
        write(out_unitp,*) 'nb_var',para_Op%mole%nb_var
        write(out_unitp,*)
        CALL write_param_Op(para_Op)
        IF (allocated(para_Op%Cmat)) CALL Write_Mat(para_Op%Cmat,out_unitp,3)
        IF (allocated(para_Op%Rmat)) CALL Write_Mat(para_Op%Rmat,out_unitp,5)
        write(out_unitp,*)
        write(out_unitp,*) 'PsiBasisRep'
        CALL ecri_psi(Psi=Psi)
        write(out_unitp,*) 'ini OpPsiBasisRep'
        CALL ecri_psi(Psi=OpPsi)
      END IF
!-----------------------------------------------------------

!-----------------------------------------------------------
!     ----------------------------------------------------------------
!     - test on nb_bi : nb_bi>0
      IF (para_Op%nb_bi <= 0) THEN
        write(out_unitp,*) ' ERROR : nb_bi MUST be > 0'
        write(out_unitp,*) 'nb_bi',para_Op%nb_bi
        write(out_unitp,*) 'STOP in ',name_sub
        STOP
      END IF
!     ----------------------------------------------------------------


      para_Op%nb_OpPsi = para_Op%nb_OpPsi + 1

      !- For the allocation of OpPsi ----------------------------------
      OpPsi = Psi


      IF (para_Op%para_ReadOp%para_FileGrid%Save_MemGrid) THEN
!=====================================================================
!
!       with in memory (the grid can be readed)
!
!=====================================================================

      Psi%GridRep=.TRUE.
      CALL alloc_psi(Psi)
      OpPsi%GridRep=.TRUE.
      CALL alloc_psi(OpPsi)
      OpPsi = ZERO


!-------------------------------------------------------------
!     --- allocate OpGrid ------------------------------------
      IF (.NOT. para_Op%alloc_Grid)                                     &
             CALL alloc_para_Op(para_Op,Mat=.FALSE.,Grid=.TRUE.)
!-------------------------------------------------------------

!-------------------------------------------------------------
      IF (old .OR. psi%nb_baie>psi%nb_tot) THEN
         write(out_unitp,*) ' Impossible with old OpPsi : ',name_sub
         STOP
      END IF
!-------------------------------------------------------------

!-------------------------------------------------------------
!     --- read  OpGrid ---------------------------------------
      IF (para_Op%para_ReadOp%para_FileGrid%Read_FileGrid .AND.         &
         .NOT. para_Op%para_ReadOp%para_FileGrid%Save_MemGrid_done) THEN
        CALL read_OpGrid_OF_Op(para_Op)
      END IF
!-------------------------------------------------------------


!       write(out_unitp,*) 'para_Op%derive_termQact',iOp,para_Op%derive_termQact(:,iOp)
!       - calculation of d0d1d2psi as a function derive_termQact

        IF (.NOT. para_Op%OpGrid(iOp)%grid_zero) THEN
          CALL sub_d0d1d2PsiBasisRep_TO_GridRep(Psi,para_Op%derive_termQdyn(:,iOp))

          IF (para_Op%OpGrid(iOp)%grid_cte) THEN

            DO i1_bi=1,para_Op%nb_bie
            DO i2_bi=1,para_Op%nb_bie
              iqi1 = 1         + (i1_bi-1) * Psi%nb_qa
              fqi1 = Psi%nb_qa + (i1_bi-1) * Psi%nb_qa
              iqi2 = 1         + (i2_bi-1) * Psi%nb_qa
              fqi2 = Psi%nb_qa + (i2_bi-1) * Psi%nb_qa

              IF (Psi%cplx) THEN
                OpPsi%CvecG(iqi1:fqi1) = OpPsi%CvecG(iqi1:fqi1) +       &
                      Psi%CvecG(iqi2:fqi2) * para_Op%OpGrid(iOp)%Mat_cte(i1_bi,i2_bi)
              ELSE
                OpPsi%RvecG(iqi1:fqi1) = OpPsi%RvecG(iqi1:fqi1) +       &
                      Psi%RvecG(iqi2:fqi2) * para_Op%OpGrid(iOp)%Mat_cte(i1_bi,i2_bi)
              END IF
            END DO
            END DO
          ELSE
            DO i1_bi=1,para_Op%nb_bie
            DO i2_bi=1,para_Op%nb_bie
              iqi1 = 1         + (i1_bi-1) * Psi%nb_qa
              fqi1 = Psi%nb_qa + (i1_bi-1) * Psi%nb_qa
              iqi2 = 1         + (i2_bi-1) * Psi%nb_qa
              fqi2 = Psi%nb_qa + (i2_bi-1) * Psi%nb_qa

              IF (Psi%cplx) THEN
                OpPsi%CvecG(iqi1:fqi1) = OpPsi%CvecG(iqi1:fqi1) +         &
                   Psi%CvecG(iqi2:fqi2) *                                 &
                     para_Op%OpGrid(iOp)%Grid(:,i1_bi,i2_bi)
              ELSE
                OpPsi%RvecG(iqi1:fqi1) = OpPsi%RvecG(iqi1:fqi1) +         &
                   Psi%RvecG(iqi2:fqi2) *                                 &
                     para_Op%OpGrid(iOp)%Grid(:,i1_bi,i2_bi)
              END IF
            END DO
            END DO
            END IF
        ELSE
          IF (Psi%cplx) THEN
            OpPsi%CvecG(:) = CZERO
          ELSE
            OpPsi%RvecG(:) = ZERO
          END IF
        END IF

!     - the projection of PsiGridRep on PsiBasisRep -------------------
!     write(out_unitp,*) 'OpPsiGridRep'
!     CALL ecri_psi(Psi=OpPsi)
      CALL sub_PsiGridRep_TO_BasisRep(OpPsi)
!     --------------------------------------------------------

!     --------------------------------------------------------
      Psi%GridRep=.FALSE.
      CALL alloc_psi(Psi)
      OpPsi%GridRep=.FALSE.
      CALL alloc_psi(OpPsi)
!     --------------------------------------------------------

!=====================================================================
!
!       END H.Psi
!
!=====================================================================
      ELSE ! para_Op%para_ReadOp%para_FileGrid%Save_MemGrid = .FLASE.
         write(out_unitp,*) 'ERROR in ',name_sub
         write(out_unitp,*) 'Save_MemGrid = .FALSE. is not possible'
         STOP
      END IF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'OpPsiBasisRep'
        CALL ecri_psi(Psi=OpPsi)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------


      END SUBROUTINE sub_OpiPsi

      SUBROUTINE sub_sqRhoOVERJac_Psi(Psi,para_Op,inv)
      USE mod_system
      USE mod_psi,       ONLY : param_psi,ecri_psi
      USE mod_SetOp,     ONLY : param_Op,write_param_Op

      IMPLICIT NONE

      !----- variables for the WP --------------------------------------
      TYPE (param_psi), intent(inout)   :: Psi

      !----- variables pour la namelist minimum ------------------------
      TYPE (param_Op), intent(in)     :: para_Op
      logical,         intent(in)     :: inv


      !----- working variables -----------------------------------------
      integer          :: n
      integer :: i2_bi
      integer :: i_qa,iqi2,fqi2

      !----- for debuging ----------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_sqRhoOVERJac_Psi'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      !-----------------------------------------------------------------

      IF (para_Op%para_Tnum%nrho /= 0) RETURN

      n = para_Op%nb_tot
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub,' ',n
        write(out_unitp,*) 'nb_bie,nb_baie',para_Op%nb_bie,para_Op%nb_baie
        write(out_unitp,*) 'nb_act1',para_Op%mole%nb_act1
        write(out_unitp,*) 'nb_var',para_Op%mole%nb_var
        write(out_unitp,*)
        write(out_unitp,*) 'PsiBasisRep'
        CALL ecri_psi(Psi=Psi)
        flush(out_unitp)
      END IF
      !-----------------------------------------------------------------

      IF (.NOT. allocated(para_Op%para_AllBasis%basis_ext%sqRhoOVERJac)) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' sqRhoOVERJac MUST be on allocated!!!'
         write(out_unitp,*) ' ... You have to force it in "sub_HSOp_inact"!!'
         write(out_unitp,*) ' Check the fortran !!'
         STOP
      END IF

      IF (Psi%cplx) THEN
       IF (.NOT. allocated(Psi%CvecG)) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' psi MUST be on the grid'
         write(out_unitp,*) ' Check the fortran !!'
         STOP
       END IF

        IF (inv) THEN
          DO i2_bi=1,para_Op%nb_bie
            iqi2 = 1         + (i2_bi-1) * Psi%nb_qa
            fqi2 = Psi%nb_qa + (i2_bi-1) * Psi%nb_qa

            Psi%CvecG(iqi2:fqi2) = Psi%CvecG(iqi2:fqi2) / para_Op%para_AllBasis%basis_ext%sqRhoOVERJac(:)

          END DO
        ELSE
          DO i2_bi=1,para_Op%nb_bie
            iqi2 = 1         + (i2_bi-1) * Psi%nb_qa
            fqi2 = Psi%nb_qa + (i2_bi-1) * Psi%nb_qa

            Psi%CvecG(iqi2:fqi2) = Psi%CvecG(iqi2:fqi2) * para_Op%para_AllBasis%basis_ext%sqRhoOVERJac(:)

          END DO
        END IF

      ELSE

        IF (.NOT. allocated(Psi%RvecG)) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' psi MUST be on the grid'
          write(out_unitp,*) ' Check the fortran !!'
          STOP
        END IF

        IF (inv) THEN

          DO i2_bi=1,para_Op%nb_bie
            iqi2 = 1         + (i2_bi-1) * Psi%nb_qa
            fqi2 = Psi%nb_qa + (i2_bi-1) * Psi%nb_qa

            Psi%RvecG(iqi2:fqi2) = Psi%RvecG(iqi2:fqi2) / para_Op%para_AllBasis%basis_ext%sqRhoOVERJac(:)

          END DO
        ELSE

          DO i2_bi=1,para_Op%nb_bie
            iqi2 = 1         + (i2_bi-1) * Psi%nb_qa
            fqi2 = Psi%nb_qa + (i2_bi-1) * Psi%nb_qa

            Psi%RvecG(iqi2:fqi2) = Psi%RvecG(iqi2:fqi2) * para_Op%para_AllBasis%basis_ext%sqRhoOVERJac(:)

          END DO

        END IF


      END IF


!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'PsiGridRep * sqrt(Rho/Jac)'
        CALL ecri_psi(Psi=Psi)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------
      END SUBROUTINE sub_sqRhoOVERJac_Psi

SUBROUTINE sub_PsiDia_TO_PsiAdia_WITH_MemGrid(Psi,para_H)
  USE mod_system
  USE mod_psi,      ONLY : param_psi,ecri_psi,                                  &
                           get_CVec_OF_psi_AT_ind_a,get_RVec_OF_psi_AT_ind_a,   &
                           set_CVec_OF_psi_AT_ind_a,set_RVec_OF_psi_AT_ind_a

  USE mod_SetOp,    ONLY : param_Op,write_param_Op
  USE mod_param_SGType2
  IMPLICIT NONE

  !----- variables pour la namelist minimum ------------------------
  TYPE (param_Op)  :: para_H ! the operator is the Hamiltonian

  integer          :: iOp

  !----- variables for the WP --------------------------------------
  TYPE (param_psi)   :: Psi


  !----- working variables -----------------------------------------
  integer :: i1_bi
  integer :: i,ki,k,iq,iqbi

  integer              :: iterm

  real (kind=Rkind)    :: EigenVec(para_H%nb_bie,para_H%nb_bie)
  real (kind=Rkind)    :: V(para_H%nb_bie,para_H%nb_bie)
  real (kind=Rkind)    :: EigenVal(para_H%nb_bie)
  real (kind=Rkind)    :: Rpsi_iq(para_H%nb_bie)
  complex (kind=Rkind) :: Cpsi_iq(para_H%nb_bie)
  logical              :: GridDone
  TYPE (OldParam)      :: OldPara

  !----- for debuging ----------------------------------------------
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub='sub_PsiDia_TO_PsiAdia_WITH_MemGrid'
  logical, parameter :: debug = .FALSE.
  !logical, parameter :: debug = .TRUE.
  !-----------------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*) 'nb_bie,nb_baie',para_H%nb_bie,para_H%nb_baie
    flush(out_unitp)
    !CALL write_param_Op(para_H)
    write(out_unitp,*)
    write(out_unitp,*) 'PsiDia'
    CALL ecri_psi(Psi=Psi)
    flush(out_unitp)
  END IF
  !-----------------------------------------------------------------

  IF (Psi%cplx) THEN
    GridDone = allocated(Psi%CvecG)
  ELSE
    GridDone = allocated(Psi%RvecG)
  END IF
  IF (debug) THEN
    write(out_unitp,*) 'GridDone',GridDone
    write(out_unitp,*) 'allo RvecG',allocated(Psi%RvecG)
    write(out_unitp,*) 'allo CvecG',allocated(Psi%CvecG)
  END IF

  IF (para_H%para_ReadOp%para_FileGrid%Save_MemGrid_done .AND. GridDone) THEN

    iterm = para_H%derive_term_TO_iterm(0,0) ! The potential term index

    IF (debug) THEN
      write(out_unitp,*) 'iterm    ',iterm
      write(out_unitp,*) 'grid_zero',para_H%OpGrid(iterm)%grid_zero
      write(out_unitp,*) 'grid_cte ',para_H%OpGrid(iterm)%grid_cte
    END IF

    IF (para_H%OpGrid(iterm)%grid_zero) THEN
      CONTINUE ! wp_adia = wp_dia (nothing to be done)
    ELSE IF (para_H%OpGrid(iterm)%grid_cte) THEN

      CALL diagonalization(para_H%OpGrid(iterm)%Mat_cte,        &
                     EigenVal,EigenVec,para_H%nb_bie,2,1,.TRUE.)

      IF (Psi%cplx) THEN

        DO iq=1,Psi%nb_qa
          CALL get_CVec_OF_psi_AT_ind_a(Cpsi_iq,Psi,iq,OldPara=OldPara)
          Cpsi_iq(:) = matmul(transpose(EigenVec),Cpsi_iq)
          CALL set_CVec_OF_psi_AT_ind_a(Cpsi_iq,Psi,iq,OldPara=OldPara)
        END DO
      ELSE

        DO iq=1,Psi%nb_qa
          CALL get_RVec_OF_psi_AT_ind_a(Rpsi_iq,Psi,iq,OldPara=OldPara)
          Rpsi_iq(:) = matmul(transpose(EigenVec),Rpsi_iq)
          CALL set_RVec_OF_psi_AT_ind_a(Rpsi_iq,Psi,iq,OldPara=OldPara)
        END DO

      END IF

    ELSE

      IF (Psi%cplx) THEN

        DO iq=1,Psi%nb_qa

          V(:,:) = para_H%OpGrid(iterm)%Grid(iq,:,:)

          CALL diagonalization(V,EigenVal,EigenVec,para_H%nb_bie,2,1,.TRUE.)

          CALL get_CVec_OF_psi_AT_ind_a(Cpsi_iq,Psi,iq,OldPara=OldPara)
          Cpsi_iq(:) = matmul(transpose(EigenVec),Cpsi_iq)
          CALL set_CVec_OF_psi_AT_ind_a(Cpsi_iq,Psi,iq,OldPara=OldPara)
        END DO
      ELSE

        DO iq=1,Psi%nb_qa

          V(:,:) = para_H%OpGrid(iterm)%Grid(iq,:,:)
          CALL diagonalization(V,EigenVal,EigenVec,para_H%nb_bie,2,1,.TRUE.)

          CALL get_RVec_OF_psi_AT_ind_a(Rpsi_iq,Psi,iq,OldPara=OldPara)
          Rpsi_iq(:) = matmul(transpose(EigenVec),Rpsi_iq)
          CALL set_RVec_OF_psi_AT_ind_a(Rpsi_iq,Psi,iq,OldPara=OldPara)
        END DO

      END IF

    END IF

  ELSE
    write(out_unitp,*) 'ERROR in ',name_sub
    write(out_unitp,*) 'Save_MemGrid_done=.FALSE. is not possible in this subroutine'
    write(out_unitp,*) '  Save_MemGrid_done',para_H%para_ReadOp%para_FileGrid%Save_MemGrid_done
    write(out_unitp,*) 'OR'
    write(out_unitp,*) 'The wave function is not on the grid:'
    write(out_unitp,*) '  alloc Psi%CvecG',allocated(Psi%CvecG)
    write(out_unitp,*) '  alloc Psi%RvecG',allocated(Psi%RvecG)
    write(out_unitp,*) 'Check the fortran!'
    STOP
  END IF
  !-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'PsiAdia'
    CALL ecri_psi(Psi=Psi)
    write(out_unitp,*)
    write(out_unitp,*) 'END ',name_sub
  END IF
  !-----------------------------------------------------------

END SUBROUTINE sub_PsiDia_TO_PsiAdia_WITH_MemGrid
      SUBROUTINE sub_psiHitermPsi(Psi,iPsi,info,para_H)
      USE mod_system
      USE mod_SetOp,    ONLY : param_Op
      USE mod_psi,      ONLY : param_psi,dealloc_psi
      IMPLICIT NONE

      TYPE (param_psi)               :: Psi
      integer                        :: iPsi
      character (len=*)              :: info
      TYPE (param_Op)                :: para_H



      !----- local variables --------------------------------------
      TYPE (param_psi)           :: OpPsi


      complex(kind=Rkind) :: avOp
      integer             :: iOp
!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING sub_psiHitermPsi'
         write(out_unitp,*) 'ipsi,info',iPsi,info
       END IF
!-----------------------------------------------------------

       OpPsi = psi  ! for the initialization

       !for H
       IF (para_H%name_Op /= 'H') STOP 'wrong Operator !!'
       write(out_unitp,*) 'nb_Term',para_H%nb_Term ; flush(out_unitp)
       DO iOp=1,para_H%nb_Term

         CALL sub_PsiOpPsi(avOp,Psi,OpPsi,para_H,iOp)
         write(out_unitp,"(i0,a,i0,a,2(i0,1x),2a,f15.9,1x,f15.9)") iPsi,  &
          ' H(',iOp,') der[',para_H%derive_termQact(:,iOp),']: ',info,avOp

         flush(out_unitp)

       END DO

       CALL dealloc_psi(OpPsi)

!----------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'END sub_psiHitermPsi'
          flush(out_unitp)
        END IF
!----------------------------------------------------------


        end subroutine sub_psiHitermPsi
END MODULE mod_OpPsi
