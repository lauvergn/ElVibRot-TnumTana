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
      SUBROUTINE sub_qa_bhe(para_AllOp)
      USE mod_system
      USE mod_Op
      USE mod_MPI
      IMPLICIT NONE

!=====================================================================
!     variables
!=====================================================================

!----- variables for the construction of H ---------------------------
      TYPE (param_AllOp), intent(inout) :: para_AllOp

!----- active parameters -------------------------------------------------
      integer                          :: nb_act,nb_act1,nb_inact2n
      logical, allocatable             :: Grid_cte(:)

      real (kind=Rkind) :: max_Sii,max_Sij
      real (kind=Rkind) :: max_Hii,max_Hij


!------ working variables ---------------------------------
      integer :: i,j,iq,iqf,nb_thread,nb_Qtransfo
      integer :: iOp,iterm,id1,id2,i_e
      logical :: lformatted,freq_only

      TYPE (OldParam) :: OldPara

!----- for debuging --------------------------------------------------
      !integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_qa_bhe'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

      nb_act     = para_AllOp%tab_Op(1)%mole%nb_act
      nb_act1    = para_AllOp%tab_Op(1)%mole%nb_act1
      nb_inact2n = para_AllOp%tab_Op(1)%mole%nb_inact2n
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub

        write(out_unitp,*) 'nb_qa',para_AllOp%tab_Op(1)%nb_qa

        CALL RecWrite_basis(para_AllOp%tab_Op(1)%para_AllBasis%BasisnD)

        write(out_unitp,*) 'pot_act,HarD,pot_cplx',                     &
                 para_AllOp%tab_Op(1)%para_PES%pot_act,                 &
                 para_AllOp%tab_Op(1)%para_PES%HarD,                    &
                 para_AllOp%tab_Op(1)%para_PES%pot_cplx

        write(out_unitp,*) 'nb_act1',nb_act1
        write(out_unitp,*) 'nb_inact2n',nb_inact2n
        write(out_unitp,*) 'nb_inact',para_AllOp%tab_Op(1)%mole%nb_inact
      END IF
!-----------------------------------------------------------

      !----- built tables ----------------------------------------------
      !Define the volume element (nrho of Basis => nrho of Tnum
      CALL nrho_Basis_TO_nhro_Tnum(para_AllOp%tab_Op(1)%para_AllBasis,  &
                                   para_AllOp%tab_Op(1)%mole)

      !----- zero of max... ------------------------------------------------
      para_AllOp%tab_Op(1)%para_PES%min_pot =  huge(ONE)
      para_AllOp%tab_Op(1)%para_PES%max_pot = -huge(ONE)

      DO iOp=1,para_AllOp%nb_Op
        IF (para_AllOp%tab_Op(iOp)%n_op == -1 .AND.                     &
                para_AllOp%tab_Op(1)%BasisnD%SparseGrid_type == 4) CYCLE ! for S

        IF (.NOT. para_AllOp%tab_Op(iOp)%alloc_Grid) THEN

          CALL alloc_NParray(Grid_cte,(/ para_AllOp%tab_Op(iOp)%nb_term /),&
                            "Grid_cte",name_sub)
          Grid_cte(:) = .FALSE.

          IF (para_AllOp%tab_Op(iOp)%para_Tnum%Gcte .AND. para_AllOp%tab_Op(iOp)%n_op == 0) THEN
            iterm = para_AllOp%tab_Op(iOp)%derive_term_TO_iterm(0,0)
            Grid_cte(:)     = .TRUE.  ! KEO
            Grid_cte(iterm) = .FALSE. ! potential
          ELSE IF (para_AllOp%tab_Op(iOp)%para_Tnum%NonGcteRange(1) > 0 .AND. &
                   para_AllOp%tab_Op(iOp)%para_Tnum%NonGcteRange(2) > 0 .AND. &
                   para_AllOp%tab_Op(iOp)%n_op == 0) THEN
            Grid_cte(:)     = .TRUE.  ! KEO

            iterm = para_AllOp%tab_Op(iOp)%derive_term_TO_iterm(0,0)
            Grid_cte(iterm)     = .FALSE. ! pot

            DO i=para_AllOp%tab_Op(iOp)%para_Tnum%NonGcteRange(1),      &
                 para_AllOp%tab_Op(iOp)%para_Tnum%NonGcteRange(2)
              iterm = para_AllOp%tab_Op(iOp)%derive_term_TO_iterm(i,0)
              Grid_cte(iterm)     = .FALSE. ! KEO
            DO j=para_AllOp%tab_Op(iOp)%para_Tnum%NonGcteRange(1),      &
                 para_AllOp%tab_Op(iOp)%para_Tnum%NonGcteRange(2)
                iterm = para_AllOp%tab_Op(iOp)%derive_term_TO_iterm(i,j)
                Grid_cte(iterm)     = .FALSE. ! KEO
              END DO
            END DO

          END IF
          IF (para_AllOp%tab_Op(iOp)%para_Tnum%Gdiago .AND. para_AllOp%tab_Op(iOp)%n_op == 0) THEN
            DO i=1,para_AllOp%tab_Op(1)%mole%nb_act
            DO j=1,para_AllOp%tab_Op(1)%mole%nb_act
              IF (i /= j) THEN
                iterm = para_AllOp%tab_Op(iOp)%derive_term_TO_iterm(i,j)
                Grid_cte(iterm)     = .TRUE. ! off diag terms
              END IF
            END DO
            END DO
          END IF

          IF (para_AllOp%tab_Op(iOp)%para_Tnum%Gcte .AND. para_AllOp%tab_Op(iOp)%name_Op(1:3) == "Mu_") THEN
            Grid_cte(:) = .TRUE.  ! KEO, rotation, Mu_xx,Mu_xy...
          END IF
          IF (para_AllOp%tab_Op(iOp)%para_Tnum%Gcte .AND. para_AllOp%tab_Op(iOp)%name_Op(1:5) == "Corr_") THEN
            Grid_cte(:) = .TRUE.  ! KEO, Coriolis
          END IF

          !grid will be allocated in the action part
          CALL alloc_para_Op(para_AllOp%tab_Op(iOp),Grid=.TRUE.,Mat=.FALSE.,Grid_cte=Grid_cte)

          CALL dealloc_NParray(Grid_cte,"Grid_cte",name_sub)
        END IF
      END DO

      !----- Transfert the constant KEO to Mate_cte -----------------
      IF (para_AllOp%tab_Op(1)%para_Tnum%Gcte) THEN
        DO iOp=1,para_AllOp%nb_Op

          IF (para_AllOp%tab_Op(iOp)%n_op == 0) THEN

            DO iterm=1,para_AllOp%tab_Op(iOp)%nb_term

              IF (.NOT. para_AllOp%tab_Op(iOp)%OpGrid(iterm)%Grid_cte) CYCLE

              id1 = para_AllOp%tab_Op(iOp)%derive_termQact(1,iterm)
              id2 = para_AllOp%tab_Op(iOp)%derive_termQact(2,iterm)
              IF (id1 /= 0 .AND. id2 /= 0) THEN ! f2 for G
                IF (id1 == id2) THEN
                  DO i_e=1,para_AllOp%tab_Op(iOp)%para_PES%nb_elec
                    para_AllOp%tab_Op(iOp)%OpGrid(iterm)%Mat_cte(i_e,i_e) = &
                      -HALF* para_AllOp%tab_Op(iOp)%para_Tnum%Gref(id1,id2)
                  END DO
                ELSE
                 DO i_e=1,para_AllOp%tab_Op(iOp)%para_PES%nb_elec
                    para_AllOp%tab_Op(iOp)%OpGrid(iterm)%Mat_cte(i_e,i_e) = &
                      - para_AllOp%tab_Op(iOp)%para_Tnum%Gref(id1,id2)
                  END DO
                END IF
              END IF

            END DO

          END IF
          IF (para_AllOp%tab_Op(iOp)%name_Op(1:3) == "Mu_") THEN
            STOP 'Rot cte not yet!'
          END IF
          IF (para_AllOp%tab_Op(iOp)%name_Op(1:5) == "Corr_") THEN
            STOP 'Corriolis cte not yet!'
          END IF

        END DO
      END IF

      IF (para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Test_Grid) THEN
        IF(MPI_id==0) write(out_unitp,*) ' TEST:  Operators at Qdyn0'
      ELSE
        IF (print_level > 0) write(out_unitp,*) 'Grid qact Veff T1 T2'
      END IF

      IF (para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Type_FileGrid /= 0) THEN
        CALL Set_File_OF_tab_Op(para_AllOp%tab_Op)
      END IF

      IF (para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Type_FileGrid == 4 .OR. &
          para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Read_FileGrid) GOTO 999

      !- test ---------------------------------------------------------
      IF (para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Test_Grid) THEN
        nb_Qtransfo = para_AllOp%tab_Op(1)%mole%nb_Qtransfo
        iq=0
        freq_only = .FALSE.

        CALL sub_HSOp_inact(iq,freq_only,para_AllOp,max_Sii,max_Sij,    &
               para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Test_Grid,OldPara)
        
        write(out_unitp,*)
        write(out_unitp,*)
        CALL time_perso('sub_qa_bhe')
        write(out_unitp,*) ' VIB: END ',name_sub
        write(out_unitp,*) '================================================='
        write(out_unitp,*)
        STOP 'test ini_act'
      END IF

      !----------------------------------------------------------------
      !-- Check for a restart on HADA file ----------------------------
      IF (para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Last_GridPoint < 1 .OR.    &
          para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Last_GridPoint >           &
          para_AllOp%tab_Op(1)%nb_qa) THEN
              para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Last_GridPoint =       &
                          para_AllOp%tab_Op(1)%nb_qa
      END IF

      IF (para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%First_GridPoint < 1 .OR.  &
          para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%First_GridPoint >         &
          para_AllOp%tab_Op(1)%nb_qa) THEN
              para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%First_GridPoint = 1
      END IF

      iqf = 0
      IF (para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Save_FileGrid) THEN
        IF (para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Type_FileGrid == 0) THEN
          IF (para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Restart_Grid) THEN
            write(out_unitp,*) '----------------------------------------'
            write(out_unitp,*) 'Restart_Grid=t'
            CALL check_HADA(iqf,para_AllOp%tab_Op(1)%ComOp)
            IF (iqf > para_AllOp%tab_Op(1)%nb_qa) iqf = 0
            iqf = max(para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%First_GridPoint,iqf+1)
            write(out_unitp,*) 'First new grid point:',iqf
            write(out_unitp,*) '----------------------------------------'
          ELSE
            iqf = para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%First_GridPoint
            lformatted = para_AllOp%tab_Op(1)%ComOp%file_HADA%formatted
            IF (print_level > 1) write(out_unitp,*) 'file_HADA%formatted',lformatted
            para_AllOp%tab_Op(1)%ComOp%file_HADA%nb_thread = Grid_maxth
            CALL file_delete(para_AllOp%tab_Op(1)%ComOp%file_HADA)
            para_AllOp%tab_Op(1)%ComOp%file_HADA%formatted = lformatted
          END IF
        ELSE
          CALL Open_File_OF_tab_Op(para_AllOp%tab_Op)
          iqf = 0
        END IF
      END IF

      IF (print_level > 1) THEN
         write(out_unitp,*) 'num_grid',                                  &
         para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%First_GridPoint,&
         para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Last_GridPoint
         write(out_unitp,*) 'num_grid iqf',iqf
      END IF

      IF (para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Save_MemGrid_done .AND.    &
          para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Type_FileGrid /= 0) GOTO 999

      !-----------------------------------------------------
      !-- Multidimensional loop ----------------------------
      IF (print_level > 1) write(out_unitp,*) 'nb_thread in ',name_sub,' : ',Grid_maxth

      IF (.NOT. para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Test_Grid .AND.    &
         print_level > 0 .AND. para_AllOp%tab_Op(1)%nb_qa > max_nb_G_FOR_print) THEN
        write(out_unitp,'(a)') 'Grid (%): [--0-10-20-30-40-50-60-70-80-90-100]'
        write(out_unitp,'(a)',ADVANCE='no') 'Grid (%): ['
        CALL flush_perso(out_unitp)
      END IF


      max_Sii = ZERO
      max_Sij = ZERO

!$OMP   PARALLEL &
!$OMP   DEFAULT(NONE) &
!$OMP   SHARED(para_AllOp,max_Sii,max_Sij,iqf) &
!$OMP   PRIVATE(iq,freq_only,OldPara) &
!$OMP   NUM_THREADS(Grid_maxth)

!$OMP   DO SCHEDULE(STATIC)
        DO iq=1,iqf-1
write(out_unitp,*) 'freq_only'
          freq_only = .TRUE.
          CALL sub_HSOp_inact(iq,freq_only,para_AllOp,max_Sii,max_Sij,  &
               para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Test_Grid,OldPara)

        END DO
!$OMP   END DO


!$OMP   DO SCHEDULE(STATIC)
        DO iq=max(1,iqf),para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Last_GridPoint

          freq_only = .FALSE.
          CALL sub_HSOp_inact(iq,freq_only,para_AllOp,max_Sii,max_Sij,  &
               para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Test_Grid,OldPara)


        END DO
!$OMP   END DO

!$OMP   END PARALLEL

      IF (.NOT. para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Test_Grid .AND.  &
         print_level > 0 .AND. para_AllOp%tab_Op(1)%nb_qa > max_nb_G_FOR_print) THEN
         IF(MPI_id==0) write(out_unitp,'(a)',ADVANCE='yes') '----]'
      END IF

      DO iOp=1,para_AllOp%nb_Op
        IF (associated(para_AllOp%tab_Op(iOp)%OpGrid)) THEN
          DO iterm=1,size(para_AllOp%tab_Op(iOp)%OpGrid)
            para_AllOp%tab_Op(iOp)%OpGrid(iterm)%Grid_done = .TRUE.
          END DO
        END IF
        IF (associated(para_AllOp%tab_Op(iOp)%imOpGrid)) THEN
          DO iterm=1,size(para_AllOp%tab_Op(iOp)%imOpGrid)
            para_AllOp%tab_Op(iOp)%imOpGrid(iterm)%Grid_done = .TRUE.
          END DO
        END IF
      END DO

      CALL flush_perso(out_unitp)
      !- END multidimentional loop ---------------------------------------
      !-------------------------------------------------------------------

 999  CONTINUE

      !-------------------------------------------------------------------
      !- Analysis of the grid (zero or constant terms)
      DO iOp=1,para_AllOp%nb_Op
        CALL Analysis_OpGrid_OF_Op(para_AllOp%tab_Op(iOp))
      END DO
      !-------------------------------------------------------------------

      !-------------------------------------------------------------------
      ! close files
      CALL Close_File_OF_tab_Op(para_AllOp%tab_Op)
      !-------------------------------------------------------------------

      !-------------------------------------------------------------------

      ! write dnTError
      IF (associated(para_AllOp%tab_Op(1)%mole%tab_Cart_transfo)) THEN
      IF (para_AllOp%tab_Op(1)%mole%tab_Cart_transfo(1)%CartesianTransfo%check_dnT) THEN
        IF(MPI_id==0) write(out_unitp,*) ' Error det(dnT-1) ?',                        &
                para_AllOp%tab_Op(1)%mole%tab_Cart_transfo(1)%CartesianTransfo%dnTErr(:)
      END IF
      END IF
      !-------------------------------------------------------------------

      !-------------------------------------------------------------------
      ! test the number of elements for the RPH transfo
      IF (associated(para_AllOp%tab_Op(1)%mole%RPHTransfo) .AND. MPI_id==0) THEN
        write(out_unitp,*) '------------------------------'
        write(out_unitp,*) 'Number of RPH points (active coordiantes)', &
          size(para_AllOp%tab_Op(1)%mole%RPHTransfo%tab_RPHpara_AT_Qact1)
        write(out_unitp,*) '------------------------------'
      END IF
      !-------------------------------------------------------------------

      !-------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
      !-------------------------------------------------------------------

      IF (.NOT. para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Read_FileGrid) THEN
      IF (para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Last_GridPoint <       &
          para_AllOp%tab_Op(1)%nb_qa .OR.                                       &
          para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%First_GridPoint > 1) THEN
        write(out_unitp,*) 'WARNING : the grid is not completed'
        write(out_unitp,*) para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%First_GridPoint, &
                           para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Last_GridPoint
        write(out_unitp,*)
        write(out_unitp,*)
        CALL time_perso('sub_qa_bhe')
        write(out_unitp,*) ' VIB: END ',name_sub
        write(out_unitp,*) '================================================='
        write(out_unitp,*)
        write(out_unitp,*) 'STOP in ',name_sub
        STOP
      END IF
      END IF

      END SUBROUTINE sub_qa_bhe

