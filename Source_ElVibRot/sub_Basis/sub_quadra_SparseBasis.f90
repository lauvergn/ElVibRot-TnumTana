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
      RECURSIVE SUBROUTINE RecSparseGrid_ForDP_type1(basis_SG,          &
                          para_Tnum,mole,para_ReadOp)
      USE mod_system
      USE mod_dnSVM
      USE mod_nDindex
      USE mod_Coord_KEO
      USE mod_PrimOp
      USE mod_basis
      USE mod_Op
      USE mod_Auto_Basis
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole
      TYPE (Tnum)    :: para_Tnum

!----- for the basis set ----------------------------------------------
      TYPE (basis), intent(inout) :: basis_SG

!----- variables for the construction of H ---------------------------
      TYPE (param_ReadOp) :: para_ReadOp


!---------------------------------------------------------------------
      TYPE (basis)        :: basis_DPsave,basis_cuba
      logical             :: cuba,SG,DP
      !logical             :: tab_SG_for_primitive(basis_SG%nb_basis)


      integer             :: i,ib,nbb,ibasis,ibi,nq_DP,nq,nb,nq_cuba
      TYPE (Type_nDindex) :: nDindL
      TYPE (Type_SymAbelian), pointer    :: P_SymAbelian_save => null()

      integer             :: L,Lmin,Lmax,i_SG,DeltaL,nq_iSG,nq_SG
      real(kind=Rkind)    :: Effi
      integer :: tab_time(8) = 0
      integer :: nb_var
      character (len=Name_len)   :: basis_name

      real(kind=Rkind)    :: OldNormB
      real (kind=Rkind), allocatable :: wrho(:)


!----- function ------------------------------------------------------
      real (kind=Rkind) :: binomial
!---------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='RecSparseGrid_ForDP_type1'
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
!        write(out_unitp,*) '--------------------------'
!        CALL RecWrite_basis(basis_SG,.TRUE.)
!        write(out_unitp,*) '--------------------------'
      END IF
!-----------------------------------------------------------
      IF (debug) basis_SG%print_info_OF_basisDP    = .TRUE.
      Effi = TEN

      Lmax = basis_SG%L_SparseGrid
      Lmin = max(0,Lmax-basis_SG%nb_basis+1)


      IF (basis_SG%print_info_OF_basisDP .AND. print_level > -1) THEN
        write(out_unitp,*) '================================================='
        write(out_unitp,*) '======== SPARSE GRID type1 (coucou) ============='
        write(out_unitp,*) '======== OR cubature (hermite only) ============='
        write(out_unitp,*) '================================================='
        write(out_unitp,*) '- Sparse Grid, Lmin,Lmax:',Lmin,Lmax
      END IF

      ! first cubature with HO basis
      cuba                                = basis_SG%SparseGrid_With_Cuba
      basis_cuba%opt_param                = basis_SG%opt_param
      basis_cuba%ndim                     = basis_SG%nb_basis
      basis_cuba%active                   = .TRUE.
      basis_cuba%name                     = "cuba_HO"
      basis_cuba%type                     = 2000
      basis_cuba%check_basis              = .FALSE.
      basis_cuba%check_nq_OF_basis        = .FALSE.
      basis_cuba%print_info_OF_basisDP    = .FALSE.
      basis_cuba%dnGGRep                  = .FALSE.
      CALL alloc_init_basis(basis_cuba)

      nb_var = size(basis_SG%Tabder_Qdyn_TO_Qbasis)
      CALL alloc_NParray(basis_cuba%Tabder_Qdyn_TO_Qbasis,(/ nb_var-1 /), &
                      'basis_cuba%Tabder_Qdyn_TO_Qbasis',name_sub,(/ 0 /))
      basis_cuba%Tabder_Qdyn_TO_Qbasis(:) = basis_SG%Tabder_Qdyn_TO_Qbasis(:)

      nq = 0
      DO i=1,basis_SG%nb_basis
        basis_name = basis_SG%tab_Pbasis(i)%Pbasis%name
        CALL string_uppercase_TO_lowercase(basis_name)
        !write(out_unitp,*) 'tab_basis',i,'basis_name',basis_name
        IF (basis_name /= "ho" .AND. basis_name /= "hm" .AND.           &
            basis_name /= "hermite") THEN
          cuba = .FALSE.
          EXIT
        END IF
        basis_cuba%iQdyn(i)      = basis_SG%tab_Pbasis(i)%Pbasis%iQdyn(1)
        basis_cuba%Q0(i)         = basis_SG%tab_Pbasis(i)%Pbasis%Q0(1)
        basis_cuba%scaleQ(i)     = basis_SG%tab_Pbasis(i)%Pbasis%scaleQ(1)
        basis_cuba%opt_Q0(i)     = basis_SG%tab_Pbasis(i)%Pbasis%opt_Q0(1)
        basis_cuba%opt_scaleQ(i) = basis_SG%tab_Pbasis(i)%Pbasis%opt_scaleQ(1)

        nq = max(nq,Get_nq_FROM_l_OF_PrimBasis(basis_SG%L_SparseGrid,   &
                                         basis_SG%tab_Pbasis(i)%Pbasis))
      END DO
      basis_cuba%nb = int(basis_SG%Norm_OF_nDindB) + 1   ! to be change

      CALL Set_nq_OF_basis(basis_cuba,nq)
      IF (print_level > -1) write(out_unitp,*) 'cuba,nq(1D),nb(1D): ',cuba,nq,basis_cuba%nb
      IF (cuba) THEN
        CALL construct_primitive_basis(basis_cuba)

        nq_cuba = get_nq_FROM_basis(basis_cuba)
        IF (nq_cuba < 1) CALL Set_nq_OF_basis(basis_cuba,huge(1))
      ELSE
        nq_cuba = huge(1)
        CALL Set_nq_OF_basis(basis_cuba,huge(1))
      END IF
      nq_cuba = get_nq_FROM_basis(basis_cuba)

      ! save the informations of basis_SG to basis_DPsave
      CALL basis2TObasis1(basis_DPsave,basis_SG)

      basis_DPsave%L_SparseGrid             = -1
      basis_DPsave%SparseGrid_type          = 0

      CALL Set_nq_OF_basis(basis_DPsave,nq=0)
      CALL Basis_Grid_ParamTOBasis_Grid_Param_init(basis_DPsave%Basis_Grid_Para)

      basis_DPsave%check_basis              = .FALSE.
      basis_DPsave%check_nq_OF_basis        = .FALSE.
      basis_DPsave%print_info_OF_basisDP    = .FALSE.

      DO ib=1,basis_SG%nb_basis
        basis_DPsave%tab_Pbasis(ib)%Pbasis%check_basis              = .FALSE.
        basis_DPsave%tab_Pbasis(ib)%Pbasis%check_nq_OF_basis        = .FALSE.
        basis_DPsave%tab_Pbasis(ib)%Pbasis%print_info_OF_basisDP    = .FALSE.
      END DO



      IF (basis_SG%L_SparseGrid < 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        CALL RecWrite_basis(basis_SG,.TRUE.)
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' L_SparseGrid < 0 !!',basis_SG%L_SparseGrid
        STOP
      END IF

      ! Calculation of the number and the index of the sum of SG basis-sets
!                            MaxCoupling=basis_SG%MaxCoupling_OF_nDindB, &
      nDindL%NormWithInit = .FALSE.
      CALL init_nDindexPrim(nDindL,ndim=basis_SG%nb_basis,              &
                            nDsize=(/ (Lmax+1,i=1,basis_SG%nb_basis) /),&
                            nDinit=(/ (0,i=1,basis_SG%nb_basis) /),     &
                            type_OF_nDindex=0,                          &
                            MinNorm=real(Lmin,kind=Rkind),              &
                            MaxNorm=real(Lmax,kind=Rkind))
      !CALL Write_nDindex(nDindL)
      !IF (debug) CALL Write_nDindex(nDindL)
      basis_SG%nb_SG = nDindL%Max_nDI
      IF (basis_SG%print_info_OF_basisDP .AND. print_level > -1 .AND. MPI_id==0) THEN
        write(out_unitp,*) '- Number of grids: ',nDindL%Max_nDI
        write(out_unitp,*) '================================================='
      END IF

        IF (basis_SG%print_info_OF_basisDP .AND. print_level > -1 .AND. MPI_id==0) THEN
          write(out_unitp,*) '================================================='
          write(out_unitp,*) '=====Set-up SG primtive basis sets==============='
        END IF
        CALL alloc_array(basis_SG%tab_basisPrimSG,                      &
                                          (/ basis_SG%nb_basis,Lmax /), &
                        'basis_SG%tab_basisPrimSG',name_sub,(/ 1,0 /))

        DO L=0,Lmax
          IF (basis_SG%print_info_OF_basisDP .AND. print_level > -1 .AND. MPI_id==0)    &
            write(out_unitp,*) '================================================='

          DO ib=1,basis_SG%nb_basis

            IF (basis_DPsave%tab_Pbasis(ib)%Pbasis%Nested > 0 .AND.     &
                basis_DPsave%tab_Pbasis(ib)%Pbasis%nq_max_Nested == -1) THEN
              write(out_unitp,*) ' ERROR you must set the "nq_max_Nested" parameter'
              STOP
            END IF

            CALL basis2TObasis1(basis_SG%tab_basisPrimSG(ib,L),         &
                                basis_DPsave%tab_Pbasis(ib)%Pbasis)

            basis_SG%tab_basisPrimSG(ib,L)%L_SparseGrid    = L
            basis_SG%tab_basisPrimSG(ib,L)%packed          = .TRUE.
            basis_SG%tab_basisPrimSG(ib,L)%L_SparseBasis   = int(basis_SG%Norm_OF_nDindB)
            basis_SG%tab_basisPrimSG(ib,L)%Norm_OF_nDindB  = basis_SG%Norm_OF_nDindB

            IF (L > 1) THEN
              IF (basis_SG%tab_basisPrimSG(ib,L-1)%nqSG_SMALLER_nqDP) THEN
                basis_SG%tab_basisPrimSG(ib,L)%SparseGrid_type = basis_SG%SparseGrid_type
              ELSE
                basis_SG%tab_basisPrimSG(ib,L)%SparseGrid_type = 0
              END IF
            END IF

            CALL RecAuto_basis(para_Tnum,mole,                          &
                               basis_SG%tab_basisPrimSG(ib,L),          &
                               para_ReadOp)

            !CALL RecWriteMini_basis(basis_SG%tab_basisPrimSG(ib,L) )
          END DO

          IF (basis_SG%print_info_OF_basisDP .AND. print_level > -1) THEN
            write(out_unitp,*) 'L,SG:  ',L,                                     &
              ((basis_SG%tab_basisPrimSG(ib,L)%SparseGrid_type > 0),ib=1,basis_SG%nb_basis)
            write(out_unitp,*) 'L,nq(:)',L,                                     &
      (get_nq_FROM_basis(basis_SG%tab_basisPrimSG(ib,L)),ib=1,basis_SG%nb_basis)
            write(out_unitp,*) 'L,nb(:)',L,                                     &
                      (basis_SG%tab_basisPrimSG(ib,L)%nb,ib=1,basis_SG%nb_basis)
            write(out_unitp,*) '================================================='
          END IF
          CALL flush_perso(out_unitp)
        END DO

        ! find the optimal grid point number : cubature, SG or DP
        ! grid point number for DP basis (with Lmax)
        ! first check if nq_DP will be too huge !!
        Effi = real(huge(1),kind=Rkind)
        DO ib=1,basis_SG%nb_basis
          Effi = Effi / real(get_nq_FROM_basis(basis_SG%tab_basisPrimSG(ib,Lmax)),kind=Rkind)
          IF (Effi < ONE) THEN ! nq_DP will be too large
            Effi = -ONE
            EXIT
          END IF
        END DO

        IF (Effi > ZERO) THEN
          nq_DP = 1
          DO ib=1,basis_SG%nb_basis
            nq_DP = nq_DP * get_nq_FROM_basis(basis_SG%tab_basisPrimSG(ib,Lmax))
          END DO
        ELSE ! nq_DP will be too large
          nq_DP = huge(1)
        END IF

        IF (basis_SG%print_info_OF_basisDP .AND. print_level > -1) THEN
          IF(MPI_id==0) write(out_unitp,*) '=END Set-up SG primtive basis sets==============='
          IF(MPI_id==0) write(out_unitp,*) '================================================='
          CALL flush_perso(out_unitp)
        END IF


!      ! Just the grid point number for the SG:
!      nq_SG = 0
!      DO i_SG=1,basis_SG%nb_SG
!
!        nq_DP = 1
!        DO ib=1,basis_SG%nb_basis
!          L = nDindL%Tab_nDval(ib,i_SG)
!          nq_DP = nq_DP * get_nq_FROM_basis(basis_SG%tab_basisPrimSG(ib,L))
!        END DO
!
!        nq_SG = nq_SG + nq_DP
!
!      END DO
!      write(out_unitp,*) 'SparseGrid number         : ',basis_SG%nb_SG
!      write(out_unitp,*) 'Grid point number (for SG): ',nq_SG
!      STOP

      ! BEGINNING Sparse Grid
      ! Set up the weights and the SG basis-sets
      CALL alloc_array(basis_SG%tab_PbasisSG,(/ basis_SG%nb_SG /),      &
                      'basis_SG%tab_PbasisSG',name_sub)
      CALL alloc_NParray(basis_SG%WeightSG,(/ basis_SG%nb_SG /),          &
                      'basis_SG%WeightSG',name_sub)

      !para_mem%mem_debug = .TRUE.

      DO i_SG=1,basis_SG%nb_SG

        DeltaL = Lmax - sum(nDindL%Tab_nDval(:,i_SG))
        IF (DeltaL < 0) STOP 'DeltaL < 0'
        IF (DeltaL > basis_SG%nb_basis -1) STOP 'DeltaL > nb_basis-1'
        IF (mod(DeltaL,2) == 0) THEN
          basis_SG%WeightSG(i_SG) =  binomial(basis_SG%nb_basis-1,deltaL)
        ELSE
          basis_SG%WeightSG(i_SG) = -binomial(basis_SG%nb_basis-1,deltaL)
        END IF

        IF (basis_SG%print_info_OF_basisDP .AND. print_level > -1) THEN
          write(out_unitp,*) '================================================='
          write(out_unitp,*) '==== i_SG',i_SG,':',nDindL%Tab_nDval(:,i_SG)
          write(out_unitp,*) '================================================='
        END IF

        CALL alloc_array(basis_SG%tab_PbasisSG(i_SG)%Pbasis,            &
                        'basis_SG%tab_PbasisSG(i_SG)%Pbasis',name_sub)

        DO ib=1,basis_SG%nb_basis
          basis_DPsave%tab_Pbasis(ib)%Pbasis%L_SparseGrid = nDindL%Tab_nDval(ib,i_SG)

          CALL Set_nq_OF_basis(basis_DPsave%tab_Pbasis(ib)%Pbasis,nq=0)
          CALL Basis_Grid_ParamTOBasis_Grid_Param_init(basis_DPsave%tab_Pbasis(ib)%Pbasis%Basis_Grid_Para)
        END DO

        CALL basis2TObasis1(basis_SG%tab_PbasisSG(i_SG)%Pbasis,basis_DPsave)
        ! here tab_Pbasis(ib)%Pbasis is also copied. It has to be removed !!
        IF (associated(basis_SG%tab_PbasisSG(i_SG)%Pbasis%tab_Pbasis)) THEN
          DO ib=1,basis_SG%tab_PbasisSG(i_SG)%Pbasis%nb_basis
            CALL dealloc_basis(basis_SG%tab_PbasisSG(i_SG)%Pbasis%tab_Pbasis(ib)%Pbasis)
            CALL dealloc_array(basis_SG%tab_PbasisSG(i_SG)%Pbasis%tab_Pbasis(ib)%Pbasis,&
                              'basis_SG%tab_PbasisSG(i_SG)%Pbasis%tab_Pbasis(ib)%Pbasis',name_sub)
          END DO
          CALL dealloc_array(basis_SG%tab_PbasisSG(i_SG)%Pbasis%tab_Pbasis,&
                            'basis_SG%tab_PbasisSG(i_SG)%Pbasis%tab_Pbasis',name_sub)
        END IF

        CALL alloc_array(basis_SG%tab_PbasisSG(i_SG)%Pbasis%tab_Pbasis,&
                                              (/ basis_SG%nb_basis /), &
                        'basis_SG%tab_PbasisSG(i_SG)%Pbasis%tab_Pbasis',name_sub)

        DO ib=1,basis_SG%nb_basis
          L = nDindL%Tab_nDval(ib,i_SG)

          basis_SG%tab_PbasisSG(i_SG)%Pbasis%tab_Pbasis(ib)%Pbasis => &
                   basis_SG%tab_basisPrimSG(ib,L)
          !CALL RecWriteMini_basis(basis_SG%tab_PbasisSG(i_SG)%Pbasis%tab_Pbasis(ib)%Pbasis)

        END DO
        basis_SG%tab_PbasisSG(i_SG)%Pbasis%tab_basis_done   = .TRUE.
        basis_SG%tab_PbasisSG(i_SG)%Pbasis%tab_basis_linked = .TRUE.

        !CALL RecWriteMini_basis(basis_SG%tab_PbasisSG(i_SG)%Pbasis)

        !basis_SG%tab_PbasisSG(i_SG)%Pbasis%print_info_OF_basisDP = .TRUE.

        CALL RecAuto_basis(para_Tnum,mole,                              &
                           basis_SG%tab_PbasisSG(i_SG)%Pbasis,          &
                           para_ReadOp)

        IF (i_SG == 1) THEN
          !-- copy tab_PbasisSG(i_SG)%Pbasis%nDindB in basis_SG%nDindB --------- !! utile ???
          basis_SG%nDindB = basis_SG%tab_PbasisSG(1)%Pbasis%nDindB
          !write(out_unitp,*) 'nDindB'
          !CALL Write_nDindex(basis_SG%nDindB)
          basis_SG%nb = basis_SG%nDindB%Max_nDI

          ! transfert of tab_symab ---------------------------------------
          CALL SymAbelian1_TO_SymAbelian2(basis_SG%tab_PbasisSG(1)%Pbasis%  &
                                         P_SymAbelian,basis_SG%P_SymAbelian)

          CALL SymAbelian1_TO_SymAbelian2(basis_SG%tab_PbasisSG(1)%Pbasis%  &
                                         P_SymAbelian,P_SymAbelian_save)

        END IF
        CALL dealloc_nDindex(basis_SG%tab_PbasisSG(i_SG)%Pbasis%nDindB)
        CALL dealloc_array(basis_SG%tab_PbasisSG(i_SG)%Pbasis%nDindB,     &
                          'basis_SG%tab_PbasisSG(i_SG)%Pbasis%nDindB',name_sub)

        basis_SG%tab_PbasisSG(i_SG)%Pbasis%nDindB => basis_SG%nDindB

        CALL dealloc_SymAbelian(basis_SG%tab_PbasisSG(i_SG)%Pbasis%P_SymAbelian)
        basis_SG%tab_PbasisSG(i_SG)%Pbasis%P_SymAbelian => basis_SG%P_SymAbelian


        IF (basis_SG%print_info_OF_basisDP .AND. print_level > -1) THEN
          write(out_unitp,*) '------------------------------------------------'
          write(out_unitp,*) '- grid number, weight',i_SG,basis_SG%WeightSG(i_SG)
          write(out_unitp,*) '- L(:) ',nDindL%Tab_nDval(:,i_SG)
          write(out_unitp,*) '- nq(:)',                                 &
                 (get_nq_FROM_basis(basis_SG%tab_PbasisSG(i_SG)%Pbasis% &
                           tab_Pbasis(i)%Pbasis),i=1,basis_SG%nb_basis)
          write(out_unitp,*) '- nq',get_nq_FROM_basis(basis_SG%tab_PbasisSG(i_SG)%Pbasis)
          write(out_unitp,*) '------------------------------------------------'
        END IF
        !IF (debug) CALL RecWrite_basis(basis_SG%tab_PbasisSG(i_SG)%Pbasis,write_all=.TRUE.)

        IF (basis_SG%print_info_OF_basisDP .AND. print_level > -1) THEN
          write(out_unitp,*) '================================================='
          write(out_unitp,*) '================================================='
          write(out_unitp,*) '================================================='
        END IF

      END DO


      !-- Check is primitive basis sets are set up ---------------------
      DO i_SG=1,basis_SG%nb_SG
        IF (.NOT. basis_SG%tab_PbasisSG(i_SG)%Pbasis%primitive_done) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  The primitive basis sets are not set up !!'
          write(out_unitp,*) '  i_SG: ',i_SG
          STOP
        END IF
      END DO
      basis_SG%primitive_done = .TRUE.

      ! number of total SG grid points
      nq_SG = 0
      DO i_SG=1,basis_SG%nb_SG
        nq_SG = nq_SG + get_nq_FROM_basis(basis_SG%tab_PbasisSG(i_SG)%Pbasis)
      END DO
      CALL Set_nq_OF_basis(basis_SG,nq_SG)

      !-- Packed the basis if needed -------------------------------
      !write(out_unitp,*) 'CALL1 pack_basis from ',name_sub
      !write(out_unitp,*) 'Pack basis_SG? ',basis_SG%packed
      CALL pack_basis(basis_SG,sortX=.TRUE.)
      !basis_DPsave%SparseGrid               = .TRUE. !????
      basis_DPsave%SparseGrid_type = 1
      nq_SG = get_nq_FROM_basis(basis_SG)

      !para_mem%mem_debug = .FALSE.


      ! TEST for the use of SG or DP basis or cuba
      !nq_cuba = huge(1)
      IF (print_level > -1) THEN
        IF (nq_DP == huge(1)) THEN
          write(out_unitp,*) 'Grid point number (for DP)     : huge(1)'
        ELSE
          write(out_unitp,*) 'Grid point number (for DP)     : ',nq_DP
        END IF
        IF (nq_cuba == huge(1)) THEN
          write(out_unitp,*) 'Grid point number (for cuba)   : huge(1)'
        ELSE
          write(out_unitp,*) 'Grid point number (for cuba)   : ',nq_cuba
        END IF
          write(out_unitp,*) 'Smolyak Grid number            : ',basis_SG%nb_SG
          write(out_unitp,*) 'Grid point number (for smolyak): ',get_nq_FROM_basis(basis_SG)
      END IF

      cuba = .FALSE.
      DP   = .FALSE.
      SG   = .FALSE.
      IF (nq_SG   < nq_DP .AND. nq_SG <  nq_cuba .AND. basis_SG%SparseGrid_With_Smolyak)  SG = .TRUE.
      IF (nq_cuba < nq_DP .AND. nq_cuba <= nq_SG .AND. basis_SG%SparseGrid_With_Cuba)     cuba = .TRUE.
      DP = (.NOT. SG) .AND. (.NOT. cuba) .AND. basis_SG%SparseGrid_With_DP

      IF (print_level > -1) write(out_unitp,*) 'cuba,SG,DP: ',cuba,SG,DP

      Effi = real(nq_DP,kind=Rkind)/real(nq_SG,kind=Rkind)
      IF (print_level > -1) write(out_unitp,*) 'Efficiency (SG)                : ',Effi
      Effi = real(nq_DP,kind=Rkind)/real(nq_cuba,kind=Rkind)
      IF (print_level > -1) write(out_unitp,*) 'Efficiency (cuba)              : ',Effi

      IF (cuba) THEN
        CALL basis2TObasis1(basis_SG,basis_cuba)

        !-- Packed the basis if needed -------------------------------
        !write(out_unitp,*) 'CALL2 pack_basis from ',name_sub
        CALL pack_basis(basis_SG,sortX=.TRUE.)
        !STOP 'cuba'
      ELSE IF (DP) THEN
        IF (print_level > -1) write(out_unitp,*) 'WARNNING: Efficiencies are <= 1.'
        IF (print_level > -1) write(out_unitp,*) '   The Smolyak grid or cubature are not used!'

        DO ib=1,basis_SG%nb_basis
          CALL basis2TObasis1(basis_DPsave%tab_Pbasis(ib)%Pbasis,     &
                              basis_SG%tab_basisPrimSG(ib,Lmax))

        END DO
        basis_DPsave%tab_basis_done           = .TRUE.

        DO L=0,Lmax
        DO ib=1,basis_SG%nb_basis
           CALL dealloc_basis(basis_SG%tab_basisPrimSG(ib,L))
         END DO
         END DO
         CALL dealloc_array(basis_SG%tab_basisPrimSG,                 &
                           'basis_SG%tab_basisPrimSG',name_sub)

        basis_DPsave%L_SparseGrid             = Lmax
        basis_DPsave%print_info_OF_basisDP    = .FALSE.
        basis_DPsave%SparseGrid_type          = 0

        CALL RecAuto_basis(para_Tnum,mole,basis_DPsave,para_ReadOp)

        CALL basis2TObasis1(basis_SG,basis_DPsave)

        CALL SymAbelian1_TO_SymAbelian2(basis_DPsave%P_SymAbelian,        &
                                        P_SymAbelian_save)
        !-- Packed the basis if needed -------------------------------
        !write(out_unitp,*) 'CALL2 pack_basis from ',name_sub
        CALL pack_basis(basis_SG,sortX=.TRUE.)
      ELSE ! SG=t
        ! nothing
      END IF


      CALL SymAbelian1_TO_SymAbelian2(P_SymAbelian_save,basis_SG%P_SymAbelian)
      CALL Set_nbPERsym_FROM_SymAbelian(basis_SG%P_SymAbelian)


      basis_SG%nqSG_SMALLER_nqDP = .TRUE.
      !basis_SG%nqSG_SMALLER_nqDP = (nq_SG < nq_DP)

      CALL dealloc_basis(basis_DPsave)
      CALL dealloc_basis(basis_cuba)

      IF (basis_SG%print_info_OF_basisDP .AND. print_level > -1) THEN
        CALL Write_SymAbelian(basis_SG%P_SymAbelian)
        CALL Write_nDindex(basis_SG%nDindB)
        write(out_unitp,*) '================================================='
        write(out_unitp,*) '======== END SPARSE GRID ========================'
        write(out_unitp,*) '================================================='
      END IF
      CALL flush_perso(out_unitp)

      !CALL RecWriteMini_basis(basis_SG)

!-----------------------------------------------------------
      IF (debug) THEN
        !CALL RecWrite_basis(basis_SG)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------
      END SUBROUTINE RecSparseGrid_ForDP_type1

      RECURSIVE SUBROUTINE RecSparseGrid_ForDP_type2(basis_SG,          &
                          para_Tnum,mole,para_ReadOp)
      USE mod_system
      USE mod_dnSVM
      USE mod_nDindex
      USE mod_Coord_KEO
      USE mod_PrimOp
      USE mod_basis
      USE mod_Op
      USE mod_Auto_Basis
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole
      TYPE (Tnum)    :: para_Tnum

!----- for the basis set ----------------------------------------------
      TYPE (basis), intent(inout) :: basis_SG

!----- variables for the construction of H ---------------------------
      TYPE (param_ReadOp) :: para_ReadOp




!---------------------------------------------------------------------
      real (kind=Rkind) :: binomial
!---------------------------------------------------------------------

      integer             :: LB,L,Lmin,Lmax,i_SG,DeltaL,nq_iSG,nq_SG,ib,nb,i
      integer             :: iq,nq,nqq,ndim
      integer             :: A,B,LG_L,LB_L,n
      integer             :: nDsize(basis_SG%nb_basis)

      TYPE (Type_IntVec), allocatable :: tab_i_TO_l(:)
      real (kind=Rkind), allocatable :: wrho(:)

      TYPE (Type_nDindex)        :: nDind_SmolyakRep_temp
      integer                    :: i_SGm1
      integer, allocatable       :: nDval(:)
      real (kind=Rkind), allocatable :: WeightSG1(:),WeightSG2(:)

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='RecSparseGrid_ForDP_type2'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'basis_SG%nb_basis',basis_SG%nb_basis
!        write(out_unitp,*) '--------------------------'
!        CALL RecWrite_basis(basis_SG,.TRUE.)
!        write(out_unitp,*) '--------------------------'
      END IF
!-----------------------------------------------------------

      IF (debug) basis_SG%print_info_OF_basisDP    = .TRUE.
      Lmax = basis_SG%L_SparseGrid
      Lmin = max(0,Lmax-basis_SG%nb_basis+1)
      LB   = basis_SG%L_SparseBasis

      IF (basis_SG%print_info_OF_basisDP .AND. print_level > -1) THEN
        write(out_unitp,*) '================================================='
        write(out_unitp,*) '======== SPARSE GRID type2 (coucou) ============='
        write(out_unitp,*) '================================================='
        write(out_unitp,*) '- Sparse Grid, packed   :',basis_SG%packed
        write(out_unitp,*) '- Sparse Grid, Lmin,Lmax:',Lmin,Lmax
      END IF

      DO ib=1,basis_SG%nb_basis
        basis_SG%tab_Pbasis(ib)%Pbasis%check_basis              = .TRUE.
        basis_SG%tab_Pbasis(ib)%Pbasis%check_nq_OF_basis        = .FALSE.
        basis_SG%tab_Pbasis(ib)%Pbasis%print_info_OF_basisDP    = .FALSE.
        basis_SG%tab_Pbasis(ib)%Pbasis%With_L                   = .TRUE.
      END DO


      IF (basis_SG%L_SparseGrid < 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        CALL RecWrite_basis(basis_SG,.TRUE.)
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' L_SparseGrid < 0 !!',basis_SG%L_SparseGrid
        STOP
      END IF

      allocate(tab_i_TO_l(basis_SG%nb_basis))

      IF (basis_SG%print_info_OF_basisDP .AND. print_level > -1) THEN
        write(out_unitp,*) '================================================='
        write(out_unitp,*) '=====Set-up SG primtive basis sets==============='
      END IF

      CALL alloc_array(basis_SG%tab_basisPrimSG,                      &
                                         (/Lmax,basis_SG%nb_basis/),  &
                      'basis_SG%tab_basisPrimSG',name_sub, (/0,1/) )


      DO ib=1,basis_SG%nb_basis
        DO L=0,Lmax

          IF (debug) THEN
            write(out_unitp,*) '================================================='
            write(out_unitp,*) '===L,ib: ',L,ib,'==============='
          END IF

          CALL basis2TObasis1(basis_SG%tab_basisPrimSG(L,ib),           &
                                         basis_SG%tab_Pbasis(ib)%Pbasis)
          ! change nb, nq (function of L)
          IF (basis_SG%tab_basisPrimSG(L,ib)%nb_basis < 1) THEN
            LG_L = L ! old
            LB_L = LB !old
            !write(out_unitp,*) 'primitive basis'
          ELSE
            basis_SG%tab_basisPrimSG(L,ib)%L_TO_nq%A = 0
            CALL init_Basis_L_TO_n(basis_SG%tab_basisPrimSG(L,ib)%L_TO_nq,Lmax=L)
            LG_L = get_n_FROM_Basis_L_TO_n(basis_SG%tab_basisPrimSG(L,ib)%L_TO_nq,L)

            !CALL Write_Basis_L_TO_n(basis_SG%tab_basisPrimSG(L,ib)%L_TO_nq)


            basis_SG%tab_basisPrimSG(L,ib)%L_TO_nb%A = 0
            CALL init_Basis_L_TO_n(basis_SG%tab_basisPrimSG(L,ib)%L_TO_nb,Lmax=LB)
            LB_L = get_n_FROM_Basis_L_TO_n(basis_SG%tab_basisPrimSG(L,ib)%L_TO_nb,LB)
            !CALL Write_Basis_L_TO_n(basis_SG%tab_basisPrimSG(L,ib)%L_TO_nb)
            !write(out_unitp,*) 'L ,LG_L',L,LG_L
            !write(out_unitp,*) 'LB,LB_L',LB,min(LG_L,LB_L)
            !write(out_unitp,*) 'not primitive basis'
          END IF

          basis_SG%tab_basisPrimSG(L,ib)%L_SparseGrid    = LG_L
          basis_SG%tab_basisPrimSG(L,ib)%L_SparseBasis   = min(LG_L,LB_L)
          basis_SG%tab_basisPrimSG(L,ib)%Norm_OF_nDindB  = min(LG_L,LB_L)
          basis_SG%tab_basisPrimSG(L,ib)%Type_OF_nDindB  = 0

          basis_SG%tab_basisPrimSG(L,ib)%packed = .TRUE.

          basis_SG%tab_basisPrimSG(L,ib)%print_info_OF_basisDP =        &
                                          basis_SG%print_info_OF_basisDP
          basis_SG%tab_basisPrimSG(L,ib)%print_info_OF_basisDP = .FALSE.

          IF (.NOT. basis_SG%check_nq_OF_basis)                         &
               basis_SG%tab_basisPrimSG(L,ib)%check_nq_OF_basis = .FALSE.
          IF (.NOT. basis_SG%check_basis)                               &
                    basis_SG%tab_basisPrimSG(L,ib)%check_basis = .FALSE.

          CALL RecAuto_basis(para_Tnum,mole,                            &
                             basis_SG%tab_basisPrimSG(L,ib),para_ReadOp)

          CALL sort_basis(basis_SG%tab_basisPrimSG(L,ib))

          nDsize(ib) = basis_SG%tab_basisPrimSG(L,ib)%nb
          IF (debug) write(out_unitp,*) 'L,ib',L,ib,'nb:',                      &
                                       basis_SG%tab_basisPrimSG(L,ib)%nb

          IF (debug) write(out_unitp,*) 'primtive basis sets of SG,L,ib',L,ib,' done'
          CALL flush_perso(out_unitp)

        END DO
        IF(MPI_id==0) THEN
          write(out_unitp,*) ib,'nb(L)',(basis_SG%tab_basisPrimSG(L,ib)%nb,L=0,Lmax)
          write(out_unitp,*) ib,'nq(L)',(get_nb_FROM_basis(basis_SG%tab_basisPrimSG(L,ib)),L=0,Lmax)
        ENDIF
      END DO
      basis_SG%primitive_done = .TRUE.
      IF (basis_SG%print_info_OF_basisDP .AND. print_level > -1 .AND. MPI_id==0) THEN
        write(out_unitp,*) '=END Set-up SG primtive basis sets==============='
        write(out_unitp,*) '================================================='
        CALL flush_perso(out_unitp)
      END IF


      DO ib=1,basis_SG%nb_basis
        nb = basis_SG%tab_basisPrimSG(Lmax,ib)%nb
        CALL alloc_dnSVM(tab_i_TO_l(ib),nb)
        IF (basis_SG%tab_basisPrimSG(Lmax,ib)%nb_basis < 1) THEN
          tab_i_TO_l(ib)%vec(:) = basis_SG%tab_basisPrimSG(Lmax,ib)%nDindB%Tab_L(:)
        ELSE
          DO i=1,nb
            n = basis_SG%tab_basisPrimSG(Lmax,ib)%nDindB%Tab_L(i)
            tab_i_TO_l(ib)%vec(i) = get_L_FROM_Basis_L_TO_n(basis_SG%tab_basisPrimSG(Lmax,ib)%L_TO_nb,n)
          END DO
        END IF
      END DO

      ! for the Basis functions -----------------------------------------
      IF (basis_SG%print_info_OF_basisDP .AND. print_level > -1 .AND. MPI_id==0) THEN
        write(out_unitp,*) '============ Set nDindB'
        CALL flush_perso(out_unitp)
      END IF
      CALL dealloc_nDindex(basis_SG%nDindB)
      CALL init_nDindexPrim(basis_SG%nDindB,basis_SG%nb_basis,nDsize,   &
                        type_OF_nDindex=3,Lmax=LB,tab_i_TO_l=tab_i_TO_l)
      basis_SG%nb = basis_SG%nDindB%Max_nDI
      IF (debug) THEN
        CALL Write_nDindex(basis_SG%nDindB)
      ELSE IF (basis_SG%print_info_OF_basisDP .AND. print_level > 1 .AND. MPI_id==0) THEN
        DO i=1,basis_SG%nDindB%Max_nDI
          write(out_unitp,*) 'ib,tab_L',i,basis_SG%nDindB%Tab_nDval(:,i)
        END DO
        CALL flush_perso(out_unitp)
      END IF

      IF (basis_SG%print_info_OF_basisDP .AND. print_level > -1 .AND. MPI_id==0) THEN
        write(out_unitp,*) '============ Set nDindB: done' ; CALL flush_perso(out_unitp)
        CALL flush_perso(out_unitp)
      END IF

      IF (basis_SG%print_info_OF_basisDP .AND. print_level > -1) THEN
        write(out_unitp,*) '============ Set nDind_SmolyakRep' ; CALL flush_perso(out_unitp)
        CALL flush_perso(out_unitp)
      END IF
      ! for the Smolyak grids -----------------------------------------
      CALL dealloc_SGType2(basis_SG%para_SGType2)

      CALL init_nDindexPrim(basis_SG%para_SGType2%nDind_SmolyakRep, &
                         basis_SG%nb_basis,nDsize,type_OF_nDindex=-4, &
                         Lmin=Lmin,Lmax=Lmax)
      IF (debug) CALL Write_nDindex(basis_SG%para_SGType2%nDind_SmolyakRep)

      basis_SG%nb_SG              = basis_SG%para_SGType2%nDind_SmolyakRep%Max_nDI
      basis_SG%para_SGType2%nb_SG = basis_SG%para_SGType2%nDind_SmolyakRep%Max_nDI

      ! for the Smolyak Weights ----------------------------------------
      CALL alloc_NParray(basis_SG%WeightSG,(/basis_SG%nb_SG/),          &
                        'basis_SG%WeightSG',name_sub)
      DO i_SG=1,basis_SG%nb_SG
        DeltaL = Lmax - sum(basis_SG%para_SGType2%nDind_SmolyakRep%Tab_nDval(:,i_SG))
        IF (DeltaL < 0) STOP 'DeltaL < 0'
        IF (DeltaL > basis_SG%nb_basis -1) STOP 'DeltaL > nb_basis-1'
        IF (mod(DeltaL,2) == 0) THEN
          basis_SG%WeightSG(i_SG) =  binomial(basis_SG%nb_basis-1,deltaL)
        ELSE
          basis_SG%WeightSG(i_SG) = -binomial(basis_SG%nb_basis-1,deltaL)
        END IF
      END DO

      IF (debug) THEN
        CALL alloc_NParray(nDval,(/basis_SG%para_SGType2%nDind_SmolyakRep%ndim/),&
                          'nDval',name_sub)
        DO i_SG=1,basis_SG%nb_SG
          CALL calc_nDindex(basis_SG%para_SGType2%nDind_SmolyakRep,i_SG,nDval)
          write(out_unitp,*) 'i_SG,nDval,coef',i_SG,nDval(:),basis_SG%WeightSG(i_SG)
        END DO
        CALL dealloc_NParray(nDval,'nDval',name_sub)
      END IF
      IF (basis_SG%print_info_OF_basisDP .AND. print_level > -1 .AND. MPI_id==0) THEN
        write(out_unitp,*) '============ Set nDind_SmolyakRep: done'
        CALL flush_perso(out_unitp)
      END IF

      IF (basis_SG%print_info_OF_basisDP .AND. print_level > -1 .AND. MPI_id==0) THEN
        write(out_unitp,*) '============ Set para_SGType2%nDind_DPG'
        CALL flush_perso(out_unitp)
      END IF

      ! for the number grid points --------------------------------------
      CALL alloc_NParray(basis_SG%para_SGType2%nDind_DPG,(/ basis_SG%nb_SG /),&
                        'basis_SG%para_SGType2%nDind_DPG',name_sub)

      CALL alloc_NParray(basis_SG%para_SGType2%tab_Sum_nq_OF_SRep,      &
                         (/ basis_SG%nb_SG /),                          &
                        'basis_SG%para_SGType2%tab_Sum_nq_OF_SRep',name_sub)
      basis_SG%para_SGType2%tab_Sum_nq_OF_SRep(:) = 0

      CALL alloc_NParray(basis_SG%para_SGType2%tab_nq_OF_SRep,(/ basis_SG%nb_SG /),&
                        'basis_SG%para_SGType2%tab_nq_OF_SRep',name_sub)
      basis_SG%para_SGType2%tab_nq_OF_SRep(:) = 0

      nqq = 0
      DO i_SG=1,basis_SG%nb_SG
        DO ib=1,basis_SG%nb_basis
          L = basis_SG%para_SGType2%nDind_SmolyakRep%Tab_nDval(ib,i_SG)
          nDsize(ib) = get_nq_FROM_basis(basis_SG%tab_basisPrimSG(L,ib))
        END DO
        CALL init_nDindexPrim(basis_SG%para_SGType2%nDind_DPG(i_SG),     &
                              basis_SG%nb_basis,nDsize,type_OF_nDindex=-1)
        nq = basis_SG%para_SGType2%nDind_DPG(i_SG)%Max_nDI
        nqq = nqq + nq

        basis_SG%para_SGType2%tab_Sum_nq_OF_SRep(i_SG) = nqq
        basis_SG%para_SGType2%tab_nq_OF_SRep(i_SG)     = nq


        IF (debug) THEN
          write(out_unitp,*) 'nq at i_SG',i_SG,':',                     &
            basis_SG%para_SGType2%nDind_SmolyakRep%Tab_nDval(:,i_SG),':',nq
        END IF
      END DO

      CALL Set_nq_OF_basis(basis_SG,nqq)
      IF (debug) write(out_unitp,*) 'nqq',nqq

      IF (basis_SG%print_info_OF_basisDP .AND. print_level > -1 .AND. MPI_id==0) THEN
        write(out_unitp,*) '============ Set para_SGType2%nDind_DPG: done'
        CALL flush_perso(out_unitp)
      END IF


      ! set of nrho ---------------------------------------
      iq = 1
      DO ib=1,basis_SG%nb_basis
        ndim = basis_SG%tab_basisPrimSG(Lmax,ib)%ndim
        basis_SG%nrho(iq:iq+ndim-1) = basis_SG%tab_basisPrimSG(Lmax,ib)%nrho(:)
        iq = iq + ndim
      END DO

      !-- Packed the basis if needed -------------------------------
      IF (basis_SG%print_info_OF_basisDP .AND. print_level > -1) THEN
        write(out_unitp,*) '============ Set pack_basis',basis_SG%packed
        CALL flush_perso(out_unitp)
      END IF
      CALL pack_basis(basis_SG,sortX=.TRUE.)
      IF (basis_SG%print_info_OF_basisDP .AND. print_level > -1) THEN
        write(out_unitp,*) '============ Set pack_basis: done'
        CALL flush_perso(out_unitp)
      END IF

      IF (basis_SG%print_info_OF_basisDP .AND. print_level > -1) THEN
        write(out_unitp,*) '============ Set_SymAbelian_OF_BasisDP'
        CALL flush_perso(out_unitp)
      END IF
      CALL Set_SymAbelian_OF_BasisDP(basis_SG)
      IF (basis_SG%print_info_OF_basisDP .AND. print_level > -1) THEN
        write(out_unitp,*) '============ Set_SymAbelian_OF_BasisDP: done'
        CALL flush_perso(out_unitp)
      END IF

      IF (debug) THEN
        CALL Write_SymAbelian(basis_SG%P_SymAbelian)
        write(out_unitp,*) '==== nDindB ====================================='
        CALL Write_nDindex(basis_SG%nDindB)
        write(out_unitp,*) '==== nDind_SmolyakRep ========================='
        CALL Write_nDindex(basis_SG%para_SGType2%nDind_SmolyakRep)
        CALL flush_perso(out_unitp)
      END IF
      IF (basis_SG%print_info_OF_basisDP .AND. print_level > -1) THEN
        write(out_unitp,*) '================================================='
        write(out_unitp,*) '== number of DP grids (nb_SG):',basis_SG%nb_SG
        write(out_unitp,*) '================================================='
        write(out_unitp,*) '======== END SPARSE GRID ========================'
        write(out_unitp,*) '================================================='
      END IF
      CALL flush_perso(out_unitp)

      DO ib=1,size(tab_i_TO_l)
        CALL dealloc_dnSVM(tab_i_TO_l(ib))
      END DO
      deallocate(tab_i_TO_l)

!-----------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(basis_SG)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------
      END SUBROUTINE RecSparseGrid_ForDP_type2

!=======================================================================================
      RECURSIVE SUBROUTINE RecSparseGrid_ForDP_type4(basis_SG,          &
                                             para_Tnum,mole,para_ReadOp)
      USE mod_system
      USE mod_dnSVM
      USE mod_nDindex
      USE mod_Coord_KEO
      USE mod_PrimOp
      USE mod_basis
      USE mod_Op
      USE mod_Auto_Basis
      USE mod_basis_BtoG_GtoB_SGType4
      USE mod_module_DInd
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole
      TYPE (Tnum)    :: para_Tnum

!----- for the basis set ----------------------------------------------
      TYPE (basis), intent(inout) :: basis_SG

!----- variables for the construction of H ---------------------------
      TYPE (param_ReadOp) :: para_ReadOp

      integer             :: LB,L,Lmin,Lmax,i_SG,DeltaL,nq_iSG,nq_SG,ib,nb,i
      integer             :: iq,nq,nqq,ndim,nbb
      integer             :: A,B,LG_L,LB_L,n
      integer             :: nDsize(basis_SG%nb_basis)
      logical             :: with_L2max
      real(kind=Rkind)    :: SG4_Mat_size,Mat_size


      TYPE (Type_IntVec), allocatable :: tab_i_TO_l(:)
      real (kind=Rkind),  allocatable :: wrho(:)

      integer       :: tab_l(basis_SG%nb_basis)
      integer       :: tab_nb(basis_SG%nb_basis)
      integer       :: tab_ib(basis_SG%nb_basis)
      integer       :: tab_nq(basis_SG%nb_basis)

      integer       :: nDNum_OF_Lmax(basis_SG%nb_basis)
      integer       :: L1maxB,L2maxB,L1maxG,L2maxG


      logical       :: Print_basis
      TYPE (basis)  :: basis_temp

      character (len=:), allocatable :: fformat
!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='RecSparseGrid_ForDP_type4'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'basis_SG%nb_basis',basis_SG%nb_basis
!        write(out_unitp,*) '--------------------------'
!        CALL RecWrite_basis(basis_SG,.TRUE.)
!        write(out_unitp,*) '--------------------------'
      END IF
!-----------------------------------------------------------

      Print_basis = basis_SG%print_info_OF_basisDP .AND. print_level > -1 .OR. debug
      Print_basis = Print_basis .AND. MPI_id==0

      Lmax = basis_SG%L_SparseGrid
      Lmin = max(0,Lmax-basis_SG%nb_basis+1)
      LB   = basis_SG%L_SparseBasis

      IF (Print_basis) THEN
        write(out_unitp,*) '================================================='
        write(out_unitp,*) '======== SPARSE GRID type4 (coucou) ============='
        write(out_unitp,*) '================================================='
        write(out_unitp,*) '- Sparse Grid, packed   :',basis_SG%packed
        write(out_unitp,*) '- Sparse Grid, Lmin,Lmax:',Lmin,Lmax
        write(out_unitp,*) '- L1_SparseBasis        :',basis_SG%para_SGType2%L1_SparseBasis
        write(out_unitp,*) '- L1_SparseGrid         :',basis_SG%para_SGType2%L1_SparseGrid
        write(out_unitp,*) '- L2_SparseBasis        :',basis_SG%para_SGType2%L2_SparseBasis
        write(out_unitp,*) '- L2_SparseGrid         :',basis_SG%para_SGType2%L2_SparseGrid
      END IF

      DO ib=1,basis_SG%nb_basis
        basis_SG%tab_Pbasis(ib)%Pbasis%check_basis              = .TRUE.
        basis_SG%tab_Pbasis(ib)%Pbasis%check_nq_OF_basis        = .FALSE.
        basis_SG%tab_Pbasis(ib)%Pbasis%print_info_OF_basisDP    = .FALSE.
        basis_SG%tab_Pbasis(ib)%Pbasis%With_L                   = .TRUE.
        !CALL RecWrite_basis(basis_SG%tab_Pbasis(ib)%Pbasis,write_all=.TRUE.)

        IF (basis_SG%tab_Pbasis(ib)%Pbasis%auto_basis) THEN
          CALL AutoParam_basis(basis_SG%tab_Pbasis(ib)%Pbasis,para_Tnum, &
                               mole,para_ReadOp)
        END IF

      END DO


      IF (basis_SG%L_SparseGrid < 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        CALL RecWrite_basis(basis_SG,.TRUE.)
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' L_SparseGrid < 0 !!',basis_SG%L_SparseGrid
        STOP
      END IF

      IF (Print_basis) THEN
        write(out_unitp,*) '================================================='
        write(out_unitp,*) '=====Set-up SG primtive basis sets==============='
      END IF

      CALL alloc_array(basis_SG%tab_basisPrimSG,                      &
                                          (/Lmax,basis_SG%nb_basis/), &
                      'basis_SG%tab_basisPrimSG',name_sub, (/0,1/) )

      DO ib=1,basis_SG%nb_basis
        DO L=0,Lmax

          IF (debug) THEN
            write(out_unitp,*) '================================================='
            write(out_unitp,*) '===L,ib: ',L,ib,'==============='
          END IF

          CALL basis2TObasis1(basis_SG%tab_basisPrimSG(L,ib),           &
                                       basis_SG%tab_Pbasis(ib)%Pbasis)

          IF (L < Lmax) THEN
            basis_SG%tab_basisPrimSG(L,ib)%contrac          = .FALSE.
            basis_SG%tab_basisPrimSG(L,ib)%contrac_analysis = .FALSE.
            basis_SG%tab_basisPrimSG(L,ib)%auto_contrac     = .FALSE.
          END IF

          ! change nb, nq (function of L)
          IF (basis_SG%tab_basisPrimSG(L,ib)%nb_basis < 1) THEN
            LG_L = L ! old
            LB_L = min(L,LB) !old
            IF (debug) write(out_unitp,*) 'primitive basis'
          ELSE
            basis_SG%tab_basisPrimSG(L,ib)%L_TO_nq%A = 0
            CALL init_Basis_L_TO_n(basis_SG%tab_basisPrimSG(L,ib)%L_TO_nq,Lmax=L)
            LG_L = get_n_FROM_Basis_L_TO_n(basis_SG%tab_basisPrimSG(L,ib)%L_TO_nq,L)
            !CALL Write_Basis_L_TO_n(basis_SG%tab_basisPrimSG(L,ib)%L_TO_nq)

            basis_SG%tab_basisPrimSG(L,ib)%L_TO_nb%A = 0
            CALL init_Basis_L_TO_n(basis_SG%tab_basisPrimSG(L,ib)%L_TO_nb,Lmax=L)
            LB_L = get_n_FROM_Basis_L_TO_n(basis_SG%tab_basisPrimSG(L,ib)%L_TO_nb,L)
            !CALL Write_Basis_L_TO_n(basis_SG%tab_basisPrimSG(L,ib)%L_TO_nb)
            IF (debug) write(out_unitp,*) 'L ,LG_L',L,LG_L
            IF (debug) write(out_unitp,*) 'LB,LB_L',LB,LB_L
            IF (Print_basis) write(out_unitp,*) 'not primitive basis'
          END IF

          basis_SG%tab_basisPrimSG(L,ib)%L_SparseGrid    = LG_L
          basis_SG%tab_basisPrimSG(L,ib)%L_SparseBasis   = LG_L !LB_L
          basis_SG%tab_basisPrimSG(L,ib)%Norm_OF_nDindB  = basis_SG%tab_basisPrimSG(L,ib)%L_SparseBasis
          basis_SG%tab_basisPrimSG(L,ib)%Type_OF_nDindB  = 0

          basis_SG%tab_basisPrimSG(L,ib)%packed = .TRUE.

          basis_SG%tab_basisPrimSG(L,ib)%print_info_OF_basisDP =        &
                                            basis_SG%print_info_OF_basisDP
          basis_SG%tab_basisPrimSG(L,ib)%print_info_OF_basisDP = .FALSE.

          IF (.NOT. basis_SG%check_nq_OF_basis)                         &
                 basis_SG%tab_basisPrimSG(L,ib)%check_nq_OF_basis = .FALSE.
          IF (.NOT. basis_SG%check_basis)                               &
                      basis_SG%tab_basisPrimSG(L,ib)%check_basis = .FALSE.

          CALL RecAuto_basis(para_Tnum,mole,                            &
                             basis_SG%tab_basisPrimSG(L,ib),para_ReadOp)

          CALL sort_basis(basis_SG%tab_basisPrimSG(L,ib))

          IF (L == Lmax .AND. allocated(basis_SG%tab_basisPrimSG(L,ib)%Rvec)) THEN
            basis_SG%tab_Pbasis(ib)%Pbasis%Rvec = basis_SG%tab_basisPrimSG(L,ib)%Rvec
          END IF

          nDsize(ib) = basis_SG%tab_basisPrimSG(L,ib)%nb


          IF (debug) THEN
            CALL RecWrite_basis(basis_SG%tab_basisPrimSG(L,ib),.TRUE.)
            write(out_unitp,*) 'primtive basis sets of SG,L,ib',L,ib,' done'
            CALL flush_perso(out_unitp)
          END IF

        END DO
        IF(MPI_id==0) THEN
          write(out_unitp,*) ib,'nb(L)',(basis_SG%tab_basisPrimSG(L,ib)%nb,L=0,Lmax)
          write(out_unitp,*) ib,'nq(L)',(get_nq_FROM_basis(basis_SG%tab_basisPrimSG(L,ib)),L=0,Lmax)
          ! was get_nq_FROM_basis
        ENDIF
      END DO
      basis_SG%primitive_done = .TRUE.
      IF (Print_basis) THEN
        write(out_unitp,*) '=END Set-up SG primtive basis sets==============='
        write(out_unitp,*) '================================================='
        CALL flush_perso(out_unitp)
      END IF

      ! for the Basis functions -----------------------------------------
      IF (Print_basis) THEN
        write(out_unitp,*) '============ Set nDindB'
        CALL flush_perso(out_unitp)
      END IF
      L1maxB = basis_SG%para_SGType2%L1_SparseBasis
      L2maxB = basis_SG%para_SGType2%L2_SparseBasis
      L1maxG = basis_SG%para_SGType2%L1_SparseGrid
      L2maxG = basis_SG%para_SGType2%L2_SparseGrid
      CALL dealloc_SGType2(basis_SG%para_SGType2) !! why ???
      basis_SG%para_SGType2%nb0 = 1


      ! for nDind_SmolyakRep%nDNum_OF_Lmax and nDindB%nDNum_OF_Lmax
      DO ib=1,basis_SG%nb_basis
        nDNum_OF_Lmax(ib) = basis_SG%tab_Pbasis(ib)%Pbasis%para_SGType2%Num_OF_Lmax
      END DO
      IF(MPI_id==0) THEN
        write(out_unitp,*) 'L1maxB, L2maxB (basis)',L1maxB,L2maxB
        write(out_unitp,*) 'nDNum_OF_Lmax',nDNum_OF_Lmax
      END IF

      allocate(tab_i_TO_l(basis_SG%nb_basis))
      DO ib=1,basis_SG%nb_basis
        nb = basis_SG%tab_basisPrimSG(Lmax,ib)%nb
        CALL alloc_dnSVM(tab_i_TO_l(ib),nb)
        IF (basis_SG%tab_basisPrimSG(Lmax,ib)%nb_basis < 1) THEN
          tab_i_TO_l(ib)%vec(:) = basis_SG%tab_basisPrimSG(Lmax,ib)%nDindB%Tab_L(:)
        ELSE
          DO i=1,nb
            n = basis_SG%tab_basisPrimSG(Lmax,ib)%nDindB%Tab_L(i)
            tab_i_TO_l(ib)%vec(i) = get_L_FROM_Basis_L_TO_n(basis_SG%tab_basisPrimSG(Lmax,ib)%L_TO_nb,n)
          END DO
        END IF
      END DO

      CALL dealloc_nDindex(basis_SG%nDindB)
      IF (count(nDNum_OF_Lmax == 0) == basis_SG%nb_basis .AND.          &
          basis_SG%MaxCoupling_OF_nDindB >= basis_SG%nb_basis) THEN
        !basis_SG%nDindB%packed = .FALSE.
        basis_SG%nDindB%packed = .TRUE. ! with false the mapping is too long !!
        CALL init_nDindexPrim(basis_SG%nDindB,basis_SG%nb_basis,nDsize, &
                              type_OF_nDindex=5,Lmax=LB,                &
                             MaxCoupling=basis_SG%MaxCoupling_OF_nDindB,&
                              tab_i_TO_l=tab_i_TO_l)
      ELSE
        basis_SG%nDindB%packed = .TRUE.
        CALL init_nDindexPrim(basis_SG%nDindB,                          &
                             basis_SG%nb_basis,nDsize,type_OF_nDindex=5,&
                             Lmax=LB,nDNum_OF_Lmax=nDNum_OF_Lmax,       &
                             L1max=L1maxB,L2max=L2maxB,                 &
                             MaxCoupling=basis_SG%MaxCoupling_OF_nDindB,&
                             tab_i_TO_l=tab_i_TO_l )
      END IF

      basis_SG%nb = basis_SG%nDindB%Max_nDI

      DO ib=1,size(tab_i_TO_l)
        CALL dealloc_dnSVM(tab_i_TO_l(ib))
      END DO
      deallocate(tab_i_TO_l)

      CALL dealloc_nDInd(basis_SG%nDindB%Tab_DInd)

      IF (debug) THEN
        CALL Write_nDindex(basis_SG%nDindB)
        CALL flush_perso(out_unitp)

      ELSE IF (Print_basis) THEN
        IF (allocated(basis_SG%nDindB%Tab_nDval)) THEN
          DO i=1,min(100,basis_SG%nDindB%Max_nDI)
            IF(MPI_id==0) write(out_unitp,*) 'ib,tab_L',i,basis_SG%nDindB%Tab_nDval(:,i)
          END DO
        ELSE
          CALL init_nDval_OF_nDindex(basis_SG%nDindB,tab_ib)
          DO i=1,min(100,basis_SG%nDindB%Max_nDI)
            CALL ADD_ONE_TO_nDindex(basis_SG%nDindB,tab_ib,iG=i)
            IF(MPI_id==0) write(out_unitp,*) 'ib,tab_L',i,tab_ib
          END DO
        END IF

        IF (basis_SG%nDindB%Max_nDI > 100) THEN
          IF(MPI_id==0) write(out_unitp,*) 'ib,tab_L .....'
        END IF
        CALL flush_perso(out_unitp)
      END IF

      IF (Print_basis) THEN
        write(out_unitp,*) '============ Set nDindB: done'
        CALL flush_perso(out_unitp)
      END IF

      IF (Print_basis) THEN
        write(out_unitp,*) '============ Set nDind_SmolyakRep'
        CALL flush_perso(out_unitp)
      END IF

      IF (count(nDNum_OF_Lmax == 0) == basis_SG%nb_basis .AND.          &
          basis_SG%MaxCoupling_OF_nDindB >= basis_SG%nb_basis) THEN

        !basis_SG%para_SGType2%nDind_SmolyakRep%packed = .FALSE.
        basis_SG%para_SGType2%nDind_SmolyakRep%packed = .TRUE.
        CALL init_nDindexPrim(basis_SG%para_SGType2%nDind_SmolyakRep,   &
                            basis_SG%nb_basis,nDsize,type_OF_nDindex=-5,&
                            Lmin=Lmin,Lmax=Lmax,                        &
                            MaxCoupling=basis_SG%MaxCoupling_OF_nDindB, &
                            nDinit=(/ (0,i=1,basis_SG%nb_basis) /) )
      ELSE
        IF (MPI_id==0) write(out_unitp,*) 'L1maxG, L2maxG (grid)',L1maxG,L2maxG

        basis_SG%para_SGType2%nDind_SmolyakRep%packed = .TRUE.
        CALL init_nDindexPrim(basis_SG%para_SGType2%nDind_SmolyakRep,   &
                           basis_SG%nb_basis,nDsize,type_OF_nDindex=-5, &
                           Lmin=0,Lmax=Lmax,nDNum_OF_Lmax=nDNum_OF_Lmax,&
                           L1max=L1maxG,L2max=L2maxG,                   &
                           MaxCoupling=basis_SG%MaxCoupling_OF_nDindB,  &
                           nDinit=(/ (0,i=1,basis_SG%nb_basis) /) )
      END IF

      IF (debug) CALL Write_nDindex(basis_SG%para_SGType2%nDind_SmolyakRep)

      basis_SG%nb_SG              = basis_SG%para_SGType2%nDind_SmolyakRep%Max_nDI
      basis_SG%para_SGType2%nb_SG = basis_SG%para_SGType2%nDind_SmolyakRep%Max_nDI

      ! for the Smolyak Weights ----------------------------------------
      CALL alloc_NParray(basis_SG%WeightSG,(/basis_SG%nb_SG/),          &
                        'basis_SG%WeightSG',name_sub)

      CALL calc_Weight_OF_SRep(basis_SG%WeightSG,basis_SG%para_SGType2%nDind_SmolyakRep)
      !CALL unpack_nDindex(basis_SG%para_SGType2%nDind_SmolyakRep)

      IF (Print_basis) THEN
        write(out_unitp,*) '============ Set nDind_SmolyakRep: done'
        CALL flush_perso(out_unitp)
      END IF

      IF (Print_basis) THEN
        write(out_unitp,*) '============ Set para_SGType2%nDind_DPG and para_SGType2%nDind_DPB'
        write(out_unitp,*) 'nb_SG:',basis_SG%nb_SG
        CALL flush_perso(out_unitp)
      END IF
      ! for the number grid points --------------------------------------

      CALL alloc_NParray(basis_SG%para_SGType2%tab_Sum_nq_OF_SRep,      &
                         (/ basis_SG%nb_SG /),                          &
                        'basis_SG%para_SGType2%tab_Sum_nq_OF_SRep',name_sub)
      basis_SG%para_SGType2%tab_Sum_nq_OF_SRep(:) = 0

      CALL alloc_NParray(basis_SG%para_SGType2%tab_nq_OF_SRep,(/ basis_SG%nb_SG /),&
                        'basis_SG%para_SGType2%tab_nq_OF_SRep',name_sub)
      basis_SG%para_SGType2%tab_nq_OF_SRep(:) = 0

      CALL alloc_NParray(basis_SG%para_SGType2%tab_Sum_nb_OF_SRep,      &
                         (/ basis_SG%nb_SG /),                          &
                        'basis_SG%para_SGType2%tab_Sum_nb_OF_SRep',name_sub)
      basis_SG%para_SGType2%tab_Sum_nb_OF_SRep(:) = 0

      CALL alloc_NParray(basis_SG%para_SGType2%tab_nb_OF_SRep,(/ basis_SG%nb_SG /),&
                        'basis_SG%para_SGType2%tab_nb_OF_SRep',name_sub)
      basis_SG%para_SGType2%tab_nb_OF_SRep(:) = 0

      nbb         = 0
      nqq         = 0
      SG4_Mat_size = ZERO
      Mat_size     = real(basis_SG%nb,kind=Rkind)**2

      CALL init_nDval_OF_nDindex(basis_SG%para_SGType2%nDind_SmolyakRep,tab_l)
      DO i_SG=1,basis_SG%nb_SG
        CALL ADD_ONE_TO_nDindex(basis_SG%para_SGType2%nDind_SmolyakRep,tab_l,iG=i_SG)

        tab_nq(:) = getbis_tab_nq(tab_l,basis_SG%tab_basisPrimSG)
        tab_nb(:) = getbis_tab_nb(tab_l,basis_SG%tab_basisPrimSG)

        nq          = product(tab_nq)
        nqq         = nqq         + nq

        basis_SG%para_SGType2%tab_Sum_nq_OF_SRep(i_SG) = nqq
        basis_SG%para_SGType2%tab_nq_OF_SRep(i_SG)     = nq

        nb  = product(tab_nb)
        nbb         = nbb         + nb
        basis_SG%para_SGType2%tab_Sum_nb_OF_SRep(i_SG) = nbb
        basis_SG%para_SGType2%tab_nb_OF_SRep(i_SG)     = nb

        SG4_Mat_size = SG4_Mat_size + basis_SG%WeightSG(i_SG)*real(nb,kind=Rkind)**2


        IF (debug) THEN
          fformat = '(a,i0,a,' // int_TO_char(basis_SG%nb_basis) // '(1x,i0),a,i0)'

          write(out_unitp,fformat) 'nb at i_SG ',i_SG,' : ',tab_nb(:),' : ',nb
          write(out_unitp,fformat) 'nq at i_SG ',i_SG,' : ',tab_nq(:),' : ',nq

          deallocate(fformat)
        END IF
      END DO

      IF(MPI_id==0) THEN
        write(out_unitp,*) ' max nq nb:',maxval(basis_SG%para_SGType2%tab_nq_OF_SRep), &
                                         maxval(basis_SG%para_SGType2%tab_nb_OF_SRep)
        CALL flush_perso(out_unitp)
      ENDIF

      CALL Set_nq_OF_basis(basis_SG,nqq)

      CALL Set_nDval_init_FOR_SG4(basis_SG%para_SGType2,version=1)

      ! save mapping table on certain threads only
      CALL Set_tables_FOR_SmolyakRepBasis_TO_tabPackedBasis(basis_SG)
      !CALL unpack_nDindex(basis_SG%nDindB)

      IF (Print_basis) THEN
        write(out_unitp,*) '============ Set para_SGType2%nDind_DPG and para_SGType2%nDind_DPB: done'
        CALL flush_perso(out_unitp)
      END IF

      IF(MPI_id==0) THEN
        write(out_unitp,*) 'nbb          (Smolyak Rep)',nbb
        write(out_unitp,*) 'nqq          (Smolyak Rep)',nqq
        write(out_unitp,*) 'SG4_Mat_size (Smolyak Rep)',SG4_Mat_size
        write(out_unitp,*) '    Mat_size (Smolyak Rep)',Mat_size
        write(out_unitp,*) 'Mat_size/SG4_Mat_size     ',Mat_size/SG4_Mat_size
      END IF
      CALL flush_perso(out_unitp)

      ! set of nrho ---------------------------------------
      iq = 1
      DO ib=1,basis_SG%nb_basis
        ndim = basis_SG%tab_basisPrimSG(Lmax,ib)%ndim
        basis_SG%nrho(iq:iq+ndim-1) = basis_SG%tab_basisPrimSG(Lmax,ib)%nrho(:)
        iq = iq + ndim
      END DO

      !-- Packed the basis if needed -------------------------------
      IF (Print_basis) THEN
        write(out_unitp,*) '============ pack_basis',basis_SG%packed
        CALL flush_perso(out_unitp)
      END IF
      CALL pack_basis(basis_SG,sortX=.TRUE.)
      IF (Print_basis) THEN
        write(out_unitp,*) '============ pack_basis: done'
        CALL flush_perso(out_unitp)
      END IF

      IF (Print_basis) THEN
        write(out_unitp,*) '============ Set_SymAbelian_OF_BasisDP'
        CALL flush_perso(out_unitp)
      END IF
      CALL Set_SymAbelian_OF_BasisDP(basis_SG)
      IF (Print_basis) THEN
        write(out_unitp,*) '============ Set_SymAbelian_OF_BasisDP: done'
        CALL flush_perso(out_unitp)
      END IF


      IF (debug) THEN
        CALL Write_SymAbelian(basis_SG%P_SymAbelian)
        write(out_unitp,*) '==== nDindB ====================================='
        CALL Write_nDindex(basis_SG%nDindB)
        write(out_unitp,*) '==== nDind_SmolyakRep ========================='
        CALL Write_nDindex(basis_SG%para_SGType2%nDind_SmolyakRep)
        CALL flush_perso(out_unitp)
      END IF

      IF (Print_basis) THEN
        write(out_unitp,*) '================================================='
        write(out_unitp,*) '== number of DP grids (nb_SG):',basis_SG%nb_SG
        write(out_unitp,*) '================================================='
        write(out_unitp,*) '======== END SPARSE GRID ========================'
        write(out_unitp,*) '================================================='
        CALL flush_perso(out_unitp)
      END IF

!-----------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(basis_SG)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------
      END SUBROUTINE RecSparseGrid_ForDP_type4
!=======================================================================================
