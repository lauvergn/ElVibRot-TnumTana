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
      SUBROUTINE sub_qa_bhe(para_AllOp)
      USE mod_system
      USE mod_Op
      IMPLICIT NONE

!=====================================================================
!
!     variables
!
!=====================================================================

!----- variables for the construction of H ---------------------------
      TYPE (param_AllOp) :: para_AllOp

!----- active parameters -------------------------------------------------
      integer   :: nb_act,nb_act1,nb_inact2n
      logical, allocatable           :: Grid_cte(:)

      real (kind=Rkind) :: max_Sii,max_Sij
      real (kind=Rkind) :: max_Hii,max_Hij


!------ working variables ---------------------------------
      integer :: iq,iqf,nb_thread,nb_Qtransfo
      integer :: iOp,iterm,id1,id2,i_e
      logical :: lformatted,freq_only

      TYPE (OldParam) :: OldPara

      real (kind=Rkind), allocatable  :: Qact_min(:)
      real (kind=Rkind)               :: Pot_min

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
            Grid_cte(:) = .TRUE.  ! KEO
            Grid_cte(1) = .FALSE. ! potential
          END IF
          IF (para_AllOp%tab_Op(iOp)%para_Tnum%Gcte .AND. para_AllOp%tab_Op(iOp)%name_Op(1:3) == "Mu_") THEN
            Grid_cte(:) = .TRUE.  ! KEO, rotation, Mu_xx,Mu_xy...
          END IF
          IF (para_AllOp%tab_Op(iOp)%para_Tnum%Gcte .AND. para_AllOp%tab_Op(iOp)%name_Op(1:5) == "Corr_") THEN
            Grid_cte(:) = .TRUE.  ! KEO, Coriolis
          END IF

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
        write(out_unitp,*) ' TEST:  Operators at Qdyn0'
      ELSE
        IF (print_level > 0) write(out_unitp,*) 'Grid qact Veff T1 T2'
      END IF

      IF (para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Type_FileGrid /= 0) THEN
        CALL Set_File_OF_tab_Op(para_AllOp%tab_Op)
      END IF

      IF (para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Type_FileGrid == 4) RETURN
      IF (para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Read_FileGrid) RETURN



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

      IF (para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Type_FileGrid == 0) THEN

        iqf = 0
        IF (para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Restart_Grid) THEN
          CALL check_HADA(iqf,para_AllOp%tab_Op(1)%ComOp)
          IF (iqf > para_AllOp%tab_Op(1)%nb_qa) iqf = 0
          iqf = max(para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%First_GridPoint,iqf+1)
        ELSE
          iqf = para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%First_GridPoint
          lformatted = para_AllOp%tab_Op(1)%ComOp%file_HADA%formatted
          IF (print_level > 1) write(out_unitp,*) 'file_HADA%formatted',lformatted
          IF (Grid_omp == 0) THEN
            nb_thread = 1
          ELSE
            nb_thread = Grid_maxth
          END IF
          para_AllOp%tab_Op(1)%ComOp%file_HADA%nb_thread = nb_thread
          CALL file_delete(para_AllOp%tab_Op(1)%ComOp%file_HADA)
          para_AllOp%tab_Op(1)%ComOp%file_HADA%formatted = lformatted
        END IF
      ELSE

        CALL Open_File_OF_tab_Op(para_AllOp%tab_Op)

        iqf = 0
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
      IF (Grid_omp == 0) THEN
        nb_thread = 1
      ELSE
        nb_thread = Grid_maxth
      END IF
      IF (print_level > 1) write(out_unitp,*) 'nb_thread in ',name_sub,' : ',nb_thread

      IF (.NOT. para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Test_Grid .AND.    &
         print_level > 0 .AND. para_AllOp%tab_Op(1)%nb_qa > max_nb_G_FOR_print) THEN
        write(out_unitp,'(a)') 'Grid_HADA (%): [--0-10-20-30-40-50-60-70-80-90-100]'
        write(out_unitp,'(a)',ADVANCE='no') 'Grid_HADA (%): ['
        CALL flush_perso(out_unitp)
      END IF

      Pot_min = huge(ONE)
      CALL alloc_NParray(Qact_min,(/para_AllOp%tab_Op(1)%mole%nb_var/),'Qact_min',name_sub)

!$OMP   PARALLEL &
!$OMP   DEFAULT(NONE) &
!$OMP   SHARED(para_AllOp,max_Sii,max_Sij,iqf,Pot_min,Qact_min) &
!$OMP   PRIVATE(iq,freq_only,OldPara) &
!$OMP   NUM_THREADS(nb_thread)

!$OMP   DO SCHEDULE(STATIC)
        DO iq=1,iqf-1

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

!        !$OMP  CRITICAL (sub_qa_bhe_CRIT)
!        IF (minval(para_AllOp%tab_Op(1)%OpGrid(1)%Grid(iq,:,:)) < pot_min) THEN
!          Pot_min = minval(para_AllOp%tab_Op(1)%OpGrid(1)%Grid(iq,:,:))
!          CALL Rec_Qact(Qact_min,                                       &
!                        para_AllOp%tab_Op(1)%para_AllBasis%BasisnD,iq,  &
!                        para_AllOp%tab_Op(1)%mole)
!        END IF
!        !$OMP END CRITICAL (sub_qa_bhe_CRIT)

        END DO
!$OMP   END DO

!$OMP   END PARALLEL

      IF (.NOT. para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Test_Grid .AND.  &
         print_level > 0 .AND. para_AllOp%tab_Op(1)%nb_qa > max_nb_G_FOR_print) THEN
        write(out_unitp,'(a)',ADVANCE='yes') '----]'
      END IF
      CALL flush_perso(out_unitp)
      !- END multidimentional loop ---------------------------------------
      !-------------------------------------------------------------------

 999  CONTINUE

      !-------------------------------------------------------------------
      !- Analysis of the grid (zero or constant terms)
      DO iOp=1,para_AllOp%nb_Op
        CALL Analysis_OpGrid_OF_Op(para_AllOp%tab_Op(iOp))
        !CALL write_param_Op(para_AllOp%tab_Op(iOp))
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
        write(out_unitp,*) ' Error det(dnT-1) ?',                       &
           para_AllOp%tab_Op(1)%mole%tab_Cart_transfo(1)%CartesianTransfo%dnTErr(:)
      END IF
      END IF
      !-------------------------------------------------------------------

      !-------------------------------------------------------------------
      ! test the number of elements for the RPH transfo
      IF (associated(para_AllOp%tab_Op(1)%mole%RPHTransfo)) THEN
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

      END SUBROUTINE sub_qa_bhe

      FUNCTION compar_GridPoint(Q1,Q2,n)
      USE mod_system
      IMPLICIT NONE
        integer :: compar_GridPoint
        integer, intent(in) :: n
        real(kind=Rkind), intent(in) :: Q1(n),Q2(n)
        integer :: i

        DO i=1,n
          IF (abs(Q1(i)-Q2(i)) < ONETENTH**10) THEN
            compar_GridPoint = 0
          ELSE
            IF (Q1(i) < Q2(i)) THEN
              compar_GridPoint = -1
            ELSE
              compar_GridPoint =  1
            END IF
            EXIT
          END IF
        END DO
      END FUNCTION compar_GridPoint

      SUBROUTINE Set_paraPRH(mole,para_Tnum,BasisnD)
      USE mod_system
      USE mod_Coord_KEO
      USE mod_PrimOp
      USE mod_basis
      IMPLICIT NONE

!
!=====================================================================
!
!     variables
!
!=====================================================================

!----- variables for the construction of H ---------------------------
      TYPE (basis)   :: BasisnD
      TYPE (zmatrix) :: mole
      TYPE (Tnum)    :: para_Tnum


!------ working variables ---------------------------------
      integer :: nq_part,iq,iq_list,nb_act1_RPH,nb_inact21_RPH,it,comp,iq_list_small

      real (kind=Rkind), allocatable :: Qact1_fromBasisnD(:)

      real (kind=Rkind) :: Qdyn(mole%nb_var)
      real (kind=Rkind) :: Qact(mole%nb_var),auTOcm_inv
      real (kind=Rkind), allocatable :: List_Qact1(:,:),List_tmp_Qact1(:,:)

      logical :: Find_in_List,tab_skip_transfo(mole%nb_Qtransfo)
      integer :: compar_GridPoint ! function

      TYPE (OldParam) :: OldPara ! to iq, ...



!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Set_paraPRH'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
      IF (mole%tab_Qtransfo(mole%itRPH)%skip_transfo) RETURN

      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')


      nb_act1_RPH    = mole%RPHTransfo%nb_act1
      nb_inact21_RPH = mole%RPHTransfo%nb_inact21

      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub

        write(out_unitp,*) 'nb_qa',get_nq_FROM_basis(BasisnD)

        !CALL RecWrite_basis(BasisnD)

        write(out_unitp,*) 'nb_act1_RPH',nb_act1_RPH
        write(out_unitp,*) 'nb_inact21_RPH',nb_inact21_RPH

        CALL Write_RPHTransfo(mole%RPHTransfo)

      END IF

      DO it=mole%nb_Qtransfo-1,mole%itRPH+1,-1
        tab_skip_transfo(it) = mole%tab_Qtransfo(it)%skip_transfo
        mole%tab_Qtransfo(it)%skip_transfo = .TRUE.
      END DO

      ! for RPHpara_AT_Qref
      IF (.NOT. associated(mole%RPHTransfo%RPHpara_AT_Qref)) THEN
        IF (debug) write(out_unitp,*) ' RPHpara_AT_Qref'
        CALL alloc_array(mole%RPHTransfo%RPHpara_AT_Qref,(/ 1 /),       &
                        'mole%RPHTransfo%RPHpara_AT_Qref',name_sub)

        CALL get_Qact(Qact,mole%ActiveTransfo) ! rigid, flexible coordinates

        CALL Set_RPHpara_AT_Qact1(mole%RPHTransfo%RPHpara_AT_Qref(1),   &
                                  Qact,para_Tnum,mole,mole%RPHTransfo)
        mole%RPHTransfo%init_Qref = .TRUE.

        write(out_unitp,*) ' Frequencies, normal modes at the reference geometry'

        write(out_unitp,11) Qact(1:nb_act1_RPH), &
               mole%RPHTransfo%RPHpara_AT_Qref(1)%dnEHess%d0(:)*auTOcm_inv
 11     format(' frequencies : ',30f10.4)
        write(out_unitp,*) 'dnQopt'
        CALL Write_dnVec(mole%RPHTransfo%RPHpara_AT_Qref(1)%dnQopt)
        write(out_unitp,*) 'dnC_inv'
        CALL Write_dnMat(mole%RPHTransfo%RPHpara_AT_Qref(1)%dnC_inv)
        CALL flush_perso(out_unitp)

      END IF

      ! for tab_RPHpara_AT_Qact1
      !----------------------------------------------------------------
      !--- First the number of grid points ----------------------------
      CALL time_perso('Grid RPH')
      write(out_unitp,*) 'Grid RPH'

      CALL alloc_NParray(Qact1_fromBasisnD,(/ nb_act1_RPH /),'Qact1_fromBasisnD',name_sub)

      nq_part = get_nq_FROM_basis(BasisnD)/100
      DO iq=1,get_nq_FROM_basis(BasisnD)

        IF (debug) THEN
        IF (mod(iq,nq_part)==0) THEN
          write(out_unitp,*) 'iq,nq',iq,get_nq_FROM_basis(BasisnD)
          CALL flush_perso(out_unitp)
        END IF
        END IF

        CALL Rec_Qact(Qact,BasisnD,iq,mole,OldPara)
        !write(6,*) 'iq,size(List_Qact1,dim=2),Qact',iq,size(List_Qact1,dim=2),Qact ; flush(6)
        CALL Qact_TO_Qdyn_FROM_ActiveTransfo(Qact,Qdyn,mole%ActiveTransfo)
        !write(6,*) 'Qdyn',Qdyn

        Qact1_fromBasisnD(:) = Qdyn(mole%RPHTransfo%list_QactTOQdyn(1:nb_act1_RPH))
        !write(6,*) 'Qact1_fromBasisnD',Qact1_fromBasisnD


        IF (.NOT. allocated(List_Qact1)) THEN ! first point
          CALL alloc_NParray(List_Qact1,(/ nb_act1_RPH, 1 /),'List_tmp_Qact1',name_sub)
          List_Qact1(:,1) = Qact1_fromBasisnD(:)
        ELSE

          Find_in_List = .FALSE.
          iq_list_small = 0
          DO iq_list=1,size(List_Qact1,dim=2)
            !write(6,*) 'coucou iq,iq_list',iq,iq_list ; flush(6)
            comp = compar_GridPoint(List_Qact1(:,iq_list),Qact1_fromBasisnD,nb_act1_RPH)
            IF (comp == -1) iq_list_small = iq_list
            Find_in_List = (comp == 0)
            !write(6,*) 'coucou iq,iq_list,comp',iq,iq_list,comp ; flush(6)

            !Find_in_List = (sum(abs(List_Qact1(:,iq_list)-Qact1_fromBasisnD)) <= ONETENTH**5)
            IF (Find_in_List) EXIT
          END DO

          IF (.NOT. Find_in_List) THEN ! add the new point in the list

            CALL alloc_NParray(List_tmp_Qact1,(/ nb_act1_RPH, size(List_Qact1,dim=2)+1 /),'List_tmp_Qact1',name_sub)

            ! find the position to add the point
            IF (iq_list_small == 0) THEN ! add in the first point
              List_tmp_Qact1(:,1) = Qact1_fromBasisnD(:)
              List_tmp_Qact1(:,2:size(List_Qact1,dim=2)+1) = List_Qact1(:,:)
            ELSE IF (iq_list_small == size(List_Qact1,dim=2)) THEN ! add in the last point
              List_tmp_Qact1(:,1:iq_list_small) = List_Qact1(:,1:iq_list_small)
              List_tmp_Qact1(:,  iq_list_small+1) = Qact1_fromBasisnD(:)
            ELSE
              List_tmp_Qact1(:,1:iq_list_small) = List_Qact1(:,1:iq_list_small)
              List_tmp_Qact1(:,  iq_list_small+1) = Qact1_fromBasisnD(:)
              List_tmp_Qact1(:,  iq_list_small+2:size(List_Qact1,dim=2)+1) = List_Qact1(:,iq_list_small+1:size(List_Qact1,dim=2))
            END IF
            CALL dealloc_NParray(List_Qact1,'List_Qact1',name_sub)

            CALL alloc_NParray(List_Qact1,(/ nb_act1_RPH, size(List_Qact1,dim=2)+1 /),'List_Qact1',name_sub)
            List_Qact1 = List_tmp_Qact1
            CALL dealloc_NParray(List_tmp_Qact1,'List_tmp_Qact1',name_sub)

          END IF


        END IF

      END DO
      CALL time_perso('Grid RPH')
      write(out_unitp,*) 'nb of Qact1 grid points',size(List_Qact1,dim=2)
      CALL flush_perso(out_unitp)
      !----------------------------------------------------------------
      !----------------------------------------------------------------

      !----------------------------------------------------------------
      !---- allocation of tab_RPHpara_AT_Qact1 ------------------------
      mole%RPHTransfo%nb_Qa = size(List_Qact1,dim=2)
      CALL alloc_array(mole%RPHTransfo%tab_RPHpara_AT_Qact1,(/ mole%RPHTransfo%nb_Qa /),&
                      'mole%RPHTransfo%tab_RPHpara_AT_Qact1',name_sub)
      !----------------------------------------------------------------


      !----------------------------------------------------------------
      !---- Multidimensional loop -------------------------------------
        DO iq_list=1,size(List_Qact1,dim=2)

          CALL get_Qact(Qact,mole%ActiveTransfo) ! rigid, flexible???? coordinates
          Qact(1:nb_act1_RPH) = List_Qact1(:,iq_list)

          write(out_unitp,*) 'new RPH point',iq_list
          CALL flush_perso(out_unitp)

          CALL Set_RPHpara_AT_Qact1(                                    &
                          mole%RPHTransfo%tab_RPHpara_AT_Qact1(iq_list),&
                                    Qact,para_Tnum,mole,mole%RPHTransfo)

        END DO

      write(out_unitp,*) 'END Grid RPH'
      CALL flush_perso(out_unitp)
      mole%RPHTransfo%init = .TRUE.

      DO it=mole%nb_Qtransfo-1,mole%itRPH+1,-1
        mole%tab_Qtransfo(it)%skip_transfo = tab_skip_transfo(it)
      END DO

      CALL dealloc_NParray(List_Qact1,'List_Qact1',name_sub)
      CALL dealloc_NParray(Qact1_fromBasisnD,'Qact1_fromBasisnD',name_sub)

!     -------------------------------------------------------
      IF (debug) THEN
        CALL Write_RPHTransfo(mole%RPHTransfo)
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -------------------------------------------------------
      END SUBROUTINE Set_paraPRH
