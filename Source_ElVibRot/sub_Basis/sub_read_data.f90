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

!================================================================
! ++    basis set reading (1D+nD)
!
!       type_basis_set,n_contrac,A,B ...
!================================================================
      SUBROUTINE read_basis5(BasisnD,mole)
      USE mod_system
      USE mod_basis
      use mod_Coord_KEO, only: zmatrix, alloc_array, alloc_nparray, dealloc_nparray
      IMPLICIT NONE

      !----- for the active basis set ------------------------------------
      TYPE (basis)               :: BasisnD

      !-- for Tnum and mole ----------------------------------------------
      TYPE (zmatrix), intent(in) :: mole


      integer       :: i,i0,i1,nb_act1_test,i_Qdyn,i_Qbasis
      integer       :: liste_var_basis(mole%nb_var)
      logical       :: check_ba
      TYPE (basis)  :: BasisnD_loc

      character (len=Name_len)   :: name

!---------------------------------------------------------------------
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'read_basis5'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
      END IF
!---------------------------------------------------------------------
      CALL alloc_array(BasisnD_loc%tab_Pbasis,(/ mole%nb_act1+1 /),     &
                      'BasisnD_loc%tab_Pbasis',name_sub)
      DO i=1,size(BasisnD_loc%tab_Pbasis)
        CALL alloc_array(BasisnD_loc%tab_Pbasis(i)%Pbasis,              &
                        'BasisnD_loc%tab_Pbasis(i)%Pbasis',name_sub)
      END DO


      write(out_unitp,*) '---------------------------------------'
      write(out_unitp,*) '----------- BasisnD -------------------'
      write(out_unitp,*) '---------------------------------------'
      write(out_unitp,*)
      write(out_unitp,*)


      ! parameter for BasisnD
      BasisnD_loc%type           = 1
      BasisnD_loc%name           = 'direct_prod'
      BasisnD_loc%contrac        = .FALSE.

!---------------------------------------------------------------
!    - read the parameters of all basis set --------------------
      !read(in_unitp,*) name
      name='basisnd'
      CALL string_uppercase_TO_lowercase(name)
      !write(out_unitp,*) 'name ',name

      IF (name .EQ. "nd" .OR. name .EQ. "basisnd") THEN
!        - read the parameters of all nD-basis set ---------------
         nb_act1_test = 0
         i = 1
         BasisnD_loc%nb_basis = 0
         BasisnD_loc%opt_param = 0
         DO WHILE (nb_act1_test < mole%nb_act1)

           CALL RecRead5_Basis(BasisnD_loc%tab_Pbasis(i)%Pbasis,mole)

           IF (BasisnD_loc%tab_Pbasis(i)%Pbasis%active) THEN
             BasisnD_loc%opt_param = BasisnD_loc%opt_param +            &
                              BasisnD_loc%tab_Pbasis(i)%Pbasis%opt_param

             nb_act1_test = nb_act1_test + BasisnD_loc%tab_Pbasis(i)%Pbasis%ndim
             i = i + 1
           END IF
           !write(out_unitp,*) ' nb_act1_test: ',nb_act1_test,mole%nb_act1
         END DO
         BasisnD_loc%nb_basis = i-1

         IF (nb_act1_test /= mole%nb_act1) THEN
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) 'nb_act1_test /= nb_act1 !!',nb_act1_test,mole%nb_act1
           STOP
         END IF
      ELSE
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) '   The old way to define the basis-set is not possible'
        STOP
      END IF

!     - check if the active coordinates are associated with a basis
      liste_var_basis(:) = 0
      DO i=1,BasisnD_loc%nb_basis
        DO i_Qbasis=1,BasisnD_loc%tab_Pbasis(i)%Pbasis%ndim
          i_Qdyn = BasisnD_loc%tab_Pbasis(i)%Pbasis%iQdyn(i_Qbasis)
          IF (liste_var_basis(i_Qdyn) == 1) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' The coordinate (iQdyn',i_Qdyn,                 &
                        ') is associated with TWO basis sets !!!'
            write(out_unitp,*) ' CHECK your DATA'
            STOP
          END IF
          liste_var_basis(i_Qdyn) = 1
        END DO
      END DO

      check_ba = .TRUE.
      DO i=1,mole%nb_var
        IF (mole%ActiveTransfo%list_act_OF_Qdyn(i) == 1 .AND.           &
                                           liste_var_basis(i) /= 1) THEN
          check_ba = .FALSE.
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' NO basis set for the ACTIVE coordinate',i
          write(out_unitp,*) '   list_act_OF_Qdyn',mole%ActiveTransfo%list_act_OF_Qdyn(:)
          write(out_unitp,*) '   liste_var_basis',liste_var_basis
          write(out_unitp,*) ' CHECK your DATA'
        END IF
      END DO
      IF (.NOT. check_ba) STOP

      ! set ndim, iQdyn, nrho, ....
      BasisnD_loc%ndim = 0
      DO i=1,BasisnD_loc%nb_basis
        BasisnD_loc%ndim = BasisnD_loc%ndim +                           &
                            BasisnD_loc%tab_Pbasis(i)%Pbasis%ndim
      END DO

      CALL alloc_init_basis(BasisnD_loc)

      ! transfert of iQdyn, nrho --------------------------------
      i0 = 0
      DO i=1,BasisnD_loc%nb_basis
        IF (BasisnD_loc%tab_Pbasis(i)%Pbasis%ndim == 0) CYCLE
        i1 = i0 + BasisnD_loc%tab_Pbasis(i)%Pbasis%ndim
        BasisnD_loc%iQdyn(i0+1:i1) = BasisnD_loc%tab_Pbasis(i)%Pbasis%iQdyn(:)
        BasisnD_loc%nrho(i0+1:i1)  = BasisnD_loc%tab_Pbasis(i)%Pbasis%nrho(:)
        i0 = i1
      END DO

      ! Set Tabder_Qdyn_TO_Qbasis(:)
      CALL alloc_NParray(BasisnD_loc%Tabder_Qdyn_TO_Qbasis,               &
                                                   (/ mole%nb_var /),   &
                      'BasisnD_loc%Tabder_Qdyn_TO_Qbasis',name_sub,(/ 0 /))
      BasisnD_loc%Tabder_Qdyn_TO_Qbasis(:) = 0
      DO i=1,BasisnD_loc%ndim
        BasisnD_loc%Tabder_Qdyn_TO_Qbasis(BasisnD_loc%iQdyn(i)) = i
      END DO

      BasisnD_loc%active            = .TRUE.
      write(out_unitp,*)
      IF (BasisnD_loc%cplx)       write(out_unitp,*) 'BasisnD is COMPLEX'
      IF (.NOT. BasisnD_loc%cplx) write(out_unitp,*) 'BasisnD is REAL'
      write(out_unitp,*)
      write(out_unitp,*) 'Number of active basis sets:',BasisnD_loc%nb_basis

      IF (BasisnD_loc%nb_basis == 1) THEN
        write(out_unitp,*) 'WARNING: ONE layer of basis has been removed!!'
        CALL basis2TObasis1(BasisnD,BasisnD_loc%tab_Pbasis(1)%Pbasis)
      ELSE
        CALL basis2TObasis1(BasisnD,BasisnD_loc)
      END IF

      IF (BasisnD%MaxCoupling_OF_nDindB < 1) BasisnD%MaxCoupling_OF_nDindB = BasisnD%nb_basis

      write(out_unitp,*) 'Tabder_Qdyn_TO_Qbasis',BasisnD%Tabder_Qdyn_TO_Qbasis(1:mole%nb_var)
      write(out_unitp,*)
      write(out_unitp,*) 'BasisnD%opt_param',BasisnD%opt_param
      write(out_unitp,*)

      CALL dealloc_basis(BasisnD_loc)
      !CALL RecWriteMini_basis(BasisnD)

      write(out_unitp,*) '---------------------------------------'
      write(out_unitp,*) '---------------------------------------'
      write(out_unitp,*) '---------------------------------------'

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BasisnD'
        CALL RecWrite_basis(BasisnD,.TRUE.)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------
      end subroutine read_basis5
      RECURSIVE SUBROUTINE RecRead5_Basis(basis_temp,mole)

      USE mod_system
      USE mod_basis
      use mod_Coord_KEO, only: zmatrix
      IMPLICIT NONE

!----- for the active basis set ---------------------------------------
      TYPE (basis)  :: basis_temp

!----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole


      integer       :: i,i0,i1

!---------------------------------------------------------------------
      character (len=*), parameter :: name_sub='RecRead5_Basis'
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unitp,*)
      write(out_unitp,*) 'BEGINNING ',name_sub
      END IF
!---------------------------------------------------------------------

!     - read the namelist -------------------------------
      CALL read5_basis_nD(basis_temp,mole)

      IF (basis_temp%nb_basis > 0) THEN ! Direct_product, SparseBasis
        basis_temp%active = .TRUE.
        CALL alloc_array(basis_temp%tab_Pbasis,(/ basis_temp%nb_basis /),&
                        'basis_temp%tab_Pbasis',name_sub)
        basis_temp%opt_param = 0
        DO i=1,basis_temp%nb_basis
          CALL alloc_array(basis_temp%tab_Pbasis(i)%Pbasis,             &
                          'basis_temp%tab_Pbasis(i)%Pbasis',name_sub)
          IF (basis_temp%packed) basis_temp%tab_Pbasis(i)%Pbasis%packed = .TRUE.
          IF (basis_temp%With_L) basis_temp%tab_Pbasis(i)%Pbasis%With_L = .TRUE.
          CALL RecRead5_Basis(basis_temp%tab_Pbasis(i)%Pbasis,mole)

          basis_temp%opt_param = basis_temp%opt_param +                 &
                              basis_temp%tab_Pbasis(i)%Pbasis%opt_param

          basis_temp%active = basis_temp%active .AND. basis_temp%tab_Pbasis(i)%Pbasis%active
          IF (.NOT. basis_temp%active) EXIT
        END DO

        IF (basis_temp%active) THEN
          basis_temp%ndim = 0
          DO i=1,basis_temp%nb_basis
            basis_temp%ndim = basis_temp%ndim +                         &
                            basis_temp%tab_Pbasis(i)%Pbasis%ndim
          END DO


          CALL alloc_init_basis(basis_temp)

          ! transfert of iQdyn, nrho, --------------------------------
          i0 = 0
          DO i=1,basis_temp%nb_basis
            IF (basis_temp%tab_Pbasis(i)%Pbasis%ndim == 0) CYCLE
            i1 = i0 + basis_temp%tab_Pbasis(i)%Pbasis%ndim
            basis_temp%iQdyn(i0+1:i1) = basis_temp%tab_Pbasis(i)%Pbasis%iQdyn(:)
            basis_temp%nrho(i0+1:i1)  = basis_temp%tab_Pbasis(i)%Pbasis%nrho(:)
            i0 = i1
          END DO
        END IF

      END IF

      IF (basis_temp%active) THEN

        CALL alloc_NParray(basis_temp%Tabder_Qdyn_TO_Qbasis,(/ mole%nb_var /), &
                        'basis_temp%Tabder_Qdyn_TO_Qbasis',name_sub,(/ 0 /))
        basis_temp%Tabder_Qdyn_TO_Qbasis(:) = 0
        ! now basis_temp%iQdyn(:) and Tabder_Qdyn_TO_Qbasis(:) are set-up
        DO i=1,basis_temp%ndim
          basis_temp%Tabder_Qdyn_TO_Qbasis(basis_temp%iQdyn(i)) = i
        END DO
        IF (print_level > 1) write(out_unitp,*) 'Tabder_Qdyn_TO_Qbasis',&
                        basis_temp%Tabder_Qdyn_TO_Qbasis(1:mole%nb_var)

      ELSE
        CALL dealloc_basis(basis_temp)
      END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(basis_temp)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------

      END SUBROUTINE RecRead5_Basis
!================================================================
! ++    read basis_nD
!
!        Now nD<50D (change max_dim = 50)
!
!================================================================
      SUBROUTINE read5_basis_nD(basis_temp,mole)
      USE mod_system
      use mod_Constant
      use mod_Coord_KEO, only: zmatrix
      USE mod_basis
      IMPLICIT NONE

!----- for the active basis set ---------------------------------------
      integer              :: nb,nq,nbc,nqc
      integer              :: Nested,nq_max_Nested
      integer              :: SparseGrid_type
      logical              :: SparseGrid,With_L
      logical              :: SparseGrid_With_Cuba    ! When 2 or more are true, the program choses the optimal one
      logical              :: SparseGrid_With_Smolyak ! When only one is true, the program tries to use only one
      logical              :: SparseGrid_With_DP      ! Remark: when only SparseGrid_With_Cuba=T, and the grid does not exit the program stops
      integer              :: L_TO_n_type,max_nb,max_nq
      integer              :: L_SparseGrid,L_TO_nq_A,L_TO_nq_B,L_TO_nq_C,Lexpo_TO_nq
      integer              :: L1_SparseGrid,L2_SparseGrid,Num_OF_Lmax
      integer              :: Type_OF_nDindB,MaxCoupling_OF_nDindB,nDinit_OF_nDindB
      integer              :: nb_OF_MinNorm_OF_nDindB,Div_nb_TO_Norm_OF_nDindB
      logical              :: contrac_WITH_nDindB
      real (kind=Rkind)    :: Norm_OF_nDindB,weight_OF_nDindB
      integer              :: L_SparseBasis,L_TO_nb_A,L_TO_nb_B,Lexpo_TO_nb,L1_SparseBasis,L2_SparseBasis
      integer, allocatable :: Tab_L_TO_n(:)
      logical              :: read_L_TO_n

      integer            :: ndim,nb_basis

      logical            :: packed,dnBBRep,contrac,contrac_analysis,read_contrac_file
      logical            :: auto_basis,auto_contrac,POGridRep,POGridRep_polyortho,make_cubature
      logical            :: restart_make_cubature
      TYPE (REAL_WU)     :: max_ene_contrac
      integer            :: max_nbc,min_nbc,nqPLUSnbc_TO_nqc
      integer            :: auto_contrac_type1_TO,auto_contrac_type21_TO
      character (len=Line_len) :: name_contrac_file

      logical            :: cplx,BasisEl

      integer, parameter :: max_dim = 50
      integer            :: iQact(max_dim)
      integer            :: iQsym(max_dim)
      integer            :: iQdyn(max_dim)
      integer            :: symab(3),index_symab(3)

      character (len=Name_len) :: name
      character (len=Name_longlen) :: dummy_name

      real (kind=Rkind)  :: cte(20,max_dim),A(max_dim),B(max_dim),Q0(max_dim),scaleQ(max_dim)
      integer :: opt_A(max_dim),opt_B(max_dim),opt_Q0(max_dim),opt_scaleQ(max_dim)


      TYPE (basis)       :: basis_temp

!----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

      integer       :: i

      NAMELIST /basis_nD/iQact,iQsym,iQdyn,name,                        &
                         nb,nq,nbc,nqc,contrac,contrac_analysis,cte,cplx,   &
                 auto_basis,A,B,Q0,scaleQ,opt_A,opt_B,opt_Q0,opt_scaleQ,&
                         symab,index_symab,                             &
                         L_TO_n_type,                                   &
                         L_SparseGrid,L_TO_nq_A,L_TO_nq_B,L_TO_nq_C,    &
                         L1_SparseGrid,L2_SparseGrid,Num_OF_Lmax,       &
                         Lexpo_TO_nq,Lexpo_TO_nb,max_nb,max_nq,         &
                         L_SparseBasis,L_TO_nb_A,L_TO_nb_B,read_L_TO_n, &
                         L1_SparseBasis,L2_SparseBasis,                 &
                         SparseGrid,SparseGrid_type,With_L,             &
                         SparseGrid_With_Cuba,SparseGrid_With_Smolyak,  &
                         SparseGrid_With_DP, &
                         Nested,nq_max_Nested,                          &
                         Type_OF_nDindB,Norm_OF_nDindB,weight_OF_nDindB,&
                         nb_OF_MinNorm_OF_nDindB,Div_nb_TO_Norm_OF_nDindB,&
                         MaxCoupling_OF_nDindB,nDinit_OF_nDindB,contrac_WITH_nDindB,   &
                         packed,dnBBRep,                             &
                         name_contrac_file,auto_contrac,                &
                         make_cubature,restart_make_cubature,           &
                         POGridRep_polyortho,POGridRep,nb_basis,        &
                      max_ene_contrac,max_nbc,min_nbc,nqPLUSnbc_TO_nqc, &
                         auto_contrac_type1_TO,auto_contrac_type21_TO

!----- for debuging --------------------------------------------------
      integer :: err_io
      character (len=*), parameter :: name_sub='read5_basis_nD'
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
      END IF
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!------- read the namelist --------------------------------------
      name                     = "0"
      packed                   = .FALSE.
      dnBBRep                  = .FALSE.
      contrac                  = .FALSE.
      contrac_analysis         = .FALSE.
      make_cubature            = .FALSE.
      restart_make_cubature    = .FALSE.
      auto_contrac             = .FALSE.
      POGridRep                = .FALSE.
      POGridRep_polyortho      = .FALSE.
      max_nbc                  = 0
      min_nbc                  = 0
      auto_contrac_type1_TO    = 100
      auto_contrac_type21_TO   = 200
      name_contrac_file        = " "
      max_ene_contrac          = REAL_WU(TEN**4,'cm-1','E') ! 10000 cm-1
      cplx                     = .FALSE.

      cte(:,:)                 = ZERO
      A(:)                     = ZERO
      B(:)                     = ZERO
      Q0(:)                    = ZERO
      scaleQ(:)                = ZERO
      opt_A(:)                 = 0
      opt_B(:)                 = 0
      opt_Q0(:)                = 0
      opt_scaleQ(:)            = 0
      auto_basis               = .FALSE.

      iQact(:)                 = 0
      iQsym(:)                 = 0
      iQdyn(:)                 = 0
      nb                       = 0
      nq                       = 0
      Nested                   = 0
      nq_max_Nested            = -1
      nbc                      = -1
      nqc                      = -1
      symab(:)                 = -1
      index_symab(:)           = -1
      nqPLUSnbc_TO_nqc         = 2
      nb_basis                 = 0

      L_TO_n_type              = 0   ! for both nq and nb
      max_nq                   = huge(1)
      L_TO_nq_A                = 1   ! nq(L) = L_TO_nq_A + L_TO_nq_B * L**expo
      L_TO_nq_B                = 1   ! nq(L) = L_TO_nq_A + L_TO_nq_B * L**expo
      Lexpo_TO_nq              = 1
      L_TO_nq_C                = 0
      L_SparseGrid             = -1
      L1_SparseGrid            = huge(1)
      L2_SparseGrid            = huge(1)
      Num_OF_Lmax              = 0

      max_nb                   = huge(1)
      L_TO_nb_A                = -1   ! nb(L) = L_TO_nb_A + L_TO_nb_B * L**expo
      L_TO_nb_B                = -1   ! nb(L) = L_TO_nb_A + L_TO_nb_B * L**expo
      Lexpo_TO_nb              = -1
      read_L_TO_n              = .FALSE.
      L_SparseBasis            = -1
      L1_SparseBasis           = huge(1)
      L2_SparseBasis           = huge(1)

      SparseGrid               = .FALSE.
      SparseGrid_type          = -1
      With_L                   = basis_temp%With_L  ! because, it can be set up in RecRead5_Basis

      SparseGrid_With_Cuba     = .TRUE.
      SparseGrid_With_Smolyak  = .TRUE.
      SparseGrid_With_DP       = .TRUE.
      Type_OF_nDindB           = 1
      Norm_OF_nDindB           = huge(ONE)
      weight_OF_nDindB         = ONE
      nDinit_OF_nDindB         = 1
      MaxCoupling_OF_nDindB    = -1
      nb_OF_MinNorm_OF_nDindB  = 1
      Div_nb_TO_Norm_OF_nDindB = 1
      contrac_WITH_nDindB      = .FALSE.

      read(in_unitp,basis_nD,IOSTAT=err_io)
      IF (err_io < 0) THEN
        write(out_unitp,basis_nD)
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) '  while reading the namelist "basis_nD"'
        write(out_unitp,*) ' end of file or end of record'
        write(out_unitp,*) ' Probably, you forget a basis set ...'
        write(out_unitp,*) ' Check your data !!'
        STOP
      END IF
      IF (err_io > 0) THEN
        write(out_unitp,basis_nD)
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) '  while reading the namelist "basis_nD"'
        write(out_unitp,*) ' Probably, some arguments of namelist are wrong.'
        write(out_unitp,*) ' Check your data !!'
        STOP
      END IF

      IF (print_level > 1 .OR. debug) THEN
         write(out_unitp,basis_nD)
      ELSE
         write(out_unitp,*) 'Basis name : ',name
      END IF

      IF (count(cte /= ZERO) > 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' Do not use cte(:,:) to define the scaling or the range parameters'
        write(out_unitp,*) ' Instead, use A and B or Q0 and scaleQ'
        write(out_unitp,*) ' CHECK your data'
        STOP
      END IF
      packed  = basis_temp%packed .OR. packed .OR. contrac .OR. auto_contrac .OR. (nb_basis < 1)

      ! Here only iQact(:) will be set-up, although iQdyn are read from the data
      IF ( count(iQsym(:) > 0) > 0 ) THEN
        write(out_unitp,*) ' WARNNING in ',name_sub
        write(out_unitp,*) '  You sould use iQdyn instead of iQsym in your basis'
      END IF
      ndim = count(iQact(:) > 0)
      IF (ndim == 0) THEN
        ndim = count(iQdyn(:) > 0)
        IF (ndim == 0) THEN
          ndim = count(iQsym(:) > 0)
          iQdyn(:) = iQsym(:)
        END IF
        DO i=1,ndim
          iQact(i) = mole%ActiveTransfo%list_QdynTOQact(iQdyn(i))
        END DO
      END IF
      ! Here only iQact(:) is set-up
      basis_temp%active = .TRUE.
      DO i=1,ndim
        basis_temp%active = basis_temp%active .AND. (iQact(i) <= mole%nb_act1)
      END DO

      IF (.NOT. basis_temp%active) RETURN

      basis_temp%name                   = name

      CALL string_uppercase_TO_lowercase(basis_temp%name)

      IF (trim(adjustl(basis_temp%name)) == 'el') THEN
        BasisEl = .TRUE.
      ELSE
        BasisEl = .FALSE.
      END IF


      IF (ndim <= 0 .AND. nb_basis == 0 .AND. .NOT. BasisEl) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' The primitive basis has no coordinates !'
        write(out_unitp,*) ' Specified the active coordinates with iQact(:) or iQdyn(:)'
        write(out_unitp,*) ' or "nb_basis" for direct-product basis'
        write(out_unitp,*) ' STOP in ',name_sub
        write(out_unitp,basis_nD)
        STOP
      END IF

      IF (auto_contrac_type1_TO /= 0 .AND. auto_contrac_type1_TO /= 100) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' The possible values for "auto_contrac_type1_TO" are: '
         write(out_unitp,*) ' "0" or "100"'
         write(out_unitp,*) ' Your value is:',auto_contrac_type1_TO
         STOP
      END IF
      IF (auto_contrac_type21_TO /= 0 .AND. auto_contrac_type21_TO /= 100 .AND. &
          auto_contrac_type21_TO /= 20 .AND. auto_contrac_type21_TO /= 200) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' the possible values for "auto_contrac_type21_TO" are:'
         write(out_unitp,*) ' "0", "100", "20" or "200"'
         write(out_unitp,*) ' Your value is:',auto_contrac_type21_TO
         STOP
      END IF

      read_contrac_file = (name_contrac_file /= " ")

      IF (POGridRep_polyortho) POGridRep = .TRUE.
      IF (.NOT. POGridRep) nqPLUSnbc_TO_nqc = 0
      IF (nqPLUSnbc_TO_nqc < 0) nqPLUSnbc_TO_nqc = 0

      basis_temp%ndim                   = ndim
      basis_temp%nb                     = nb
      basis_temp%nbc                    = nbc
      basis_temp%nqc                    = nqc
      basis_temp%nqPLUSnbc_TO_nqc       = nqPLUSnbc_TO_nqc
      basis_temp%nq_max_Nested          = nq_max_Nested
      basis_temp%Nested                 = Nested
      CALL Set_nq_OF_basis(basis_temp,nq)

      IF (L_SparseBasis /= -1 .AND. Norm_OF_nDindB == huge(ONE)) THEN
        Type_OF_nDindB = 0
        Norm_OF_nDindB = real(L_SparseBasis,kind=Rkind)
      END IF
      IF (Norm_OF_nDindB /= huge(ONE) .AND. L_SparseBasis == -1) THEN
        L_SparseBasis = int(Norm_OF_nDindB)
      END IF

      IF (SparseGrid .AND. SparseGrid_type == -1) THEN
        SparseGrid_type = 1 !with SG (old way)
      ELSE IF (SparseGrid .AND. SparseGrid_type == 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) '  SparseGrid=T (with SG) and SparseGrid_type=0 (without SG)'
        write(out_unitp,*) '  Do not use the "SparseGrid" variable, ...'
        write(out_unitp,*) '  ... use only the "SparseGrid_type"'
        write(out_unitp,*) '  Check your data!!'
        STOP
      ELSE IF (.NOT. SparseGrid .AND. SparseGrid_type == -1) THEN
        SparseGrid_type = 0 ! without SG
      END IF
      SparseGrid = (SparseGrid_type > 0)

      IF (SparseGrid_type > 0) THEN

        IF (L_SparseGrid < 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  You cannot set-up a sparse grid with L_SparseGrid < 0'
          write(out_unitp,*) '  SparseGrid_type,L_SparseGrid',SparseGrid_type,L_SparseGrid
          STOP
        END IF

        IF (L_SparseBasis > L_SparseGrid) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  L_SparseBasis > L_SparseGrid : it is impossible'
          write(out_unitp,*) '  L_SparseBasis,L_SparseGrid',L_SparseBasis,L_SparseGrid
          STOP
        END IF

      END IF

      basis_temp%SparseGrid_type             = SparseGrid_type
      basis_temp%L_SparseGrid                = L_SparseGrid
      basis_temp%L_SparseBasis               = L_SparseBasis


      IF (Num_OF_Lmax < 0 .OR. Num_OF_Lmax > 2) Num_OF_Lmax = 0
      basis_temp%para_SGType2%Num_OF_Lmax    = Num_OF_Lmax

      IF (L1_SparseGrid  == huge(1) .AND. L1_SparseBasis < huge(1)) L1_SparseGrid  = L1_SparseBasis
      IF (L1_SparseBasis == huge(1) .AND. L1_SparseGrid  < huge(1)) L1_SparseBasis = L1_SparseGrid
      IF (L1_SparseBasis > L1_SparseGrid) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  L1_SparseBasis > L1_SparseGrid : it is impossible'
          write(out_unitp,*) '  L1_SparseBasis,L1_SparseGrid',L1_SparseBasis,L1_SparseGrid
          STOP
      END IF
      basis_temp%para_SGType2%L1_SparseGrid  = L1_SparseGrid
      basis_temp%para_SGType2%L1_SparseBasis = L1_SparseBasis

      IF (L2_SparseGrid  == huge(1) .AND. L2_SparseBasis < huge(1)) L2_SparseGrid  = L2_SparseBasis
      IF (L2_SparseBasis == huge(1) .AND. L2_SparseGrid  < huge(1)) L2_SparseBasis = L2_SparseGrid
      IF (L2_SparseBasis > L2_SparseGrid) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  L2_SparseBasis > L2_SparseGrid : it is impossible'
          write(out_unitp,*) '  L2_SparseBasis,L2_SparseGrid',L2_SparseBasis,L2_SparseGrid
          STOP
      END IF
      basis_temp%para_SGType2%L2_SparseGrid  = L2_SparseGrid
      basis_temp%para_SGType2%L2_SparseBasis = L2_SparseBasis

      IF (max_nb > max_nq ) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  max_nb MUST be <= max_nq'
          write(out_unitp,*) '  max_nb, max_nq',max_nb, max_nq
          STOP
      END IF


        IF (read_L_TO_n) THEN
          CALL alloc_NParray(Tab_L_TO_n,(/ 10 /),"Tab_L_TO_n",name_sub,(/ 0 /))

          read(in_unitp,*,IOSTAT=err_io) dummy_name,Tab_L_TO_n
          IF (err_io /= 0) THEN
            write(out_unitp,basis_nD)
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) '  while reading  "Tab_L_TO_n" for L_TO_nb'
            write(out_unitp,*) ' Probably, you some intergers are missing ...'
            write(out_unitp,*) ' => The line has to be like that (with 11 integers):'
            write(out_unitp,*) ' l_to_nb 1 2 3 4 5 6   6 6 6 6 6'
            write(out_unitp,*) ' Check your data !!'
            STOP
          END IF
          write(out_unitp,'(2a,11(x,i0))') 'dummy_name ',dummy_name,Tab_L_TO_n
          CALL Set_Basis_L_TO_n(basis_temp%L_TO_nb,max_n=max_nb,Tab_L_TO_n=Tab_L_TO_n)



          read(in_unitp,*,IOSTAT=err_io) dummy_name,Tab_L_TO_n
          IF (err_io /= 0) THEN
            write(out_unitp,basis_nD)
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) '  while reading "Tab_L_TO_n" for L_TO_nq'
            write(out_unitp,*) ' Probably, you some intergers are missing ...'
            write(out_unitp,*) ' => The line has to like that (with 11 integers):'
            write(out_unitp,*) ' l_to_nq 1 2 3 4 5 6   6 6 6 6 6'
            write(out_unitp,*) ' Check your data !!'
            STOP
          END IF
          write(out_unitp,'(2a,11(x,i0))') 'dummy_name ',dummy_name,Tab_L_TO_n
          CALL Set_Basis_L_TO_n(basis_temp%L_TO_nq,max_n=max_nq,Tab_L_TO_n=Tab_L_TO_n)

          CALL dealloc_NParray(Tab_L_TO_n,"Tab_L_TO_n",name_sub)

        ELSE
          IF (L_TO_nb_A   == -1) L_TO_nb_A   = L_TO_nq_A
          IF (L_TO_nb_B   == -1) L_TO_nb_B   = L_TO_nq_B
          IF (Lexpo_TO_nb == -1) Lexpo_TO_nb = Lexpo_TO_nq

          IF (L_TO_nb_B > L_TO_nq_B .OR. L_TO_nb_A > L_TO_nq_A .OR. Lexpo_TO_nb > Lexpo_TO_nq) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) '  L_TO_nb_B MUST be <= L_TO_nq_B'
            write(out_unitp,*) '  L_TO_nb_B,L_TO_nq_B',L_TO_nb_B,L_TO_nq_B

            write(out_unitp,*) '       and'

            write(out_unitp,*) '  L_TO_nb_A MUST be <= L_TO_nq_A'
            write(out_unitp,*) '  L_TO_nb_A,L_TO_nq_A',L_TO_nb_A,L_TO_nq_A

            write(out_unitp,*) '       and'

            write(out_unitp,*) '  Lexpo_TO_nb MUST be <= Lexpo_TO_nq'
            write(out_unitp,*) '  Lexpo_TO_nb,Lexpo_TO_nq',Lexpo_TO_nb,Lexpo_TO_nq

            STOP
          END IF
          IF (L_TO_nq_C == 0) L_TO_nq_C = L_TO_nq_B

          CALL Set_Basis_L_TO_n(basis_temp%L_TO_nb,L_TO_nb_A,L_TO_nb_B, &
                  expo=Lexpo_TO_nb,max_n=max_nb,L_TO_n_type=L_TO_n_type)
          CALL Set_Basis_L_TO_n(basis_temp%L_TO_nq,L_TO_nq_A,L_TO_nq_B, &
             L_TO_nq_C,Lexpo_TO_nq,max_n=max_nq,L_TO_n_type=L_TO_n_type)


        END IF

        basis_temp%SparseGrid_With_Cuba    = SparseGrid_With_Cuba
        basis_temp%SparseGrid_With_Smolyak = SparseGrid_With_Smolyak
        basis_temp%SparseGrid_With_DP      = SparseGrid_With_DP

        basis_temp%With_L                  = With_L
        IF (SparseGrid_type == 2) basis_temp%With_L = .TRUE.
        IF (SparseGrid_type == 4) basis_temp%With_L = .TRUE.

        ! Be carefull we can deal with only one type of Sparse Grid in the same calculation
!        IF ( SGtype > 0 .AND. SparseGrid_type > 0) THEN
!        IF ( (SparseGrid_type-SGtype) /= 0) THEN
!          write(out_unitp,*) ' ERROR in ',name_sub
!          write(out_unitp,*) ' There are two Smolyak types',SGtype,SparseGrid_type
!          write(out_unitp,*) ' ... in the same calculation'
!          write(out_unitp,*) ' You can use only one type!'
!          write(out_unitp,*) ' Check your data !'
!          STOP
!        END IF
!        END IF
        IF (SparseGrid_type == 1 .OR. SparseGrid_type == 2 .OR. SparseGrid_type == 4) SGtype = SparseGrid_type


        IF (print_level > 1) THEN
          write(out_unitp,*) ' Parameters: nq(L) = A + B * L**expo'
          CALL Write_Basis_L_TO_n(basis_temp%L_TO_nq)
          write(out_unitp,*) ' Parameters: nb(L) = A + B * L**expo'
          CALL Write_Basis_L_TO_n(basis_temp%L_TO_nb)
        END IF


      basis_temp%Type_OF_nDindB           = Type_OF_nDindB
      basis_temp%Norm_OF_nDindB           = Norm_OF_nDindB
      basis_temp%weight_OF_nDindB         = weight_OF_nDindB
      basis_temp%nDinit_OF_nDindB         = nDinit_OF_nDindB
      basis_temp%nb_OF_MinNorm_OF_nDindB  = nb_OF_MinNorm_OF_nDindB
      basis_temp%Div_nb_TO_Norm_OF_nDindB = Div_nb_TO_Norm_OF_nDindB
      basis_temp%contrac_WITH_nDindB      = contrac_WITH_nDindB
      IF (nb_basis == 0) MaxCoupling_OF_nDindB = 1
      IF (nb_basis > 0 .AND. MaxCoupling_OF_nDindB < 1) MaxCoupling_OF_nDindB = nb_basis
      basis_temp%MaxCoupling_OF_nDindB    = MaxCoupling_OF_nDindB

      basis_temp%dnBBRep                  = dnBBRep
      basis_temp%packed                   = packed
      basis_temp%contrac                  = contrac .OR. auto_contrac .OR. contrac_analysis
      basis_temp%contrac_analysis         = contrac_analysis
      basis_temp%auto_contrac             = auto_contrac
      basis_temp%make_cubature            = make_cubature
      basis_temp%restart_make_cubature    = restart_make_cubature
      basis_temp%max_ene_contrac          = convRWU_TO_R(max_ene_contrac)
      basis_temp%max_nbc                  = max_nbc
      basis_temp%min_nbc                  = min_nbc
      basis_temp%auto_contrac_type1_TO    = auto_contrac_type1_TO
      basis_temp%auto_contrac_type21_TO   = auto_contrac_type21_TO
      basis_temp%POGridRep                = POGridRep
      basis_temp%POGridRep_polyortho      = POGridRep_polyortho
      basis_temp%read_contrac_file        = read_contrac_file
      basis_temp%file_contrac%name        = name_contrac_file
      basis_temp%cplx                     = cplx
      basis_temp%nb_basis                 = nb_basis


      CALL Set_ReadsymabOFSymAbelian(basis_temp%P_SymAbelian,symab(1))

      IF (ndim == 0) THEN
         IF (.NOT. associated(basis_temp%nDindB)) THEN
           CALL alloc_array(basis_temp%nDindB,"basis_temp%nDindB",name_sub)
         END IF
         basis_temp%nDindB%packed          = basis_temp%packed .OR. (Type_OF_nDindB == 0)
      ELSE
        CALL alloc_init_basis(basis_temp)
        basis_temp%nDindB%packed          = basis_temp%packed .OR. (Type_OF_nDindB == 0)
        ! now basis_temp%iQdyn(:)
        DO i=1,ndim
          basis_temp%iQdyn(i) = mole%ActiveTransfo%list_QactTOQdyn(iQact(i))
          IF (basis_temp%iQdyn(i) > mole%nb_var) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' iQdyn(i) is larger than nb_var',               &
                                    i,basis_temp%iQdyn(i),mole%nb_var
            STOP
          END IF
        END DO

        IF (print_level > 1) THEN
          write(out_unitp,*) 'auto_basis:    ',auto_basis
          write(out_unitp,*) 'A,B,Q0,scaleQ: ',A(1:ndim),B(1:ndim),Q0(1:ndim),scaleQ(1:ndim)
          write(out_unitp,*) 'opt of A,B,Q0,scaleQ: ',opt_A(1:ndim),opt_B(1:ndim),opt_Q0(1:ndim),opt_scaleQ(1:ndim)
        END IF
        basis_temp%A(:)      = ZERO
        basis_temp%B(:)      = ZERO
        basis_temp%Q0(:)     = ZERO
        basis_temp%scaleQ(:) = ONE

        DO i=1,ndim
          IF (A(i) /= B(i) .AND. scaleQ(i) > ZERO) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' You give the range ("A" and "B"):',A(i),B(i)
            write(out_unitp,*) '   and also the scaling factors ("scaleQ" and "Q0"):',scaleQ(i),Q0(i)
            write(out_unitp,*) ' You have to chose the range or the scaling factors'
            write(out_unitp,*) ' CHECK your data'
            write(out_unitp,basis_nD)
            STOP
          ELSE IF (A(i) > B(i) .AND. scaleQ(i) == ZERO) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' The range ("A" and "B") is :',A(i),B(i)
            write(out_unitp,*) '   "A" MUST be < "B" '
            write(out_unitp,*) ' CHECK your data'
            write(out_unitp,basis_nD)
            STOP
          ELSE IF (A(i) /= B(i) .AND. scaleQ(i) == ZERO) THEN
            basis_temp%A(i)      = A(i)
            basis_temp%B(i)      = B(i)
            basis_temp%Q0(i)     = ZERO
            basis_temp%scaleQ(i) = ONE
            basis_temp%opt_A(i)  = opt_A(i)
            basis_temp%opt_B(i)  = opt_B(i)
          ELSE IF (A(i) == B(i) .AND. scaleQ(i) > ZERO) THEN
            basis_temp%A(i)          = ZERO
            basis_temp%B(i)          = ZERO
            basis_temp%Q0(i)         = Q0(i)
            basis_temp%scaleQ(i)     = scaleQ(i)
            basis_temp%opt_Q0(i)     = opt_Q0(i)
            basis_temp%opt_scaleQ(i) = opt_scaleQ(i)
          ELSE IF (A(i) == B(i) .AND. scaleQ(i) == ZERO .AND. Q0(i) /= ZERO) THEN
            basis_temp%A(i)          = ZERO
            basis_temp%B(i)          = ZERO
            basis_temp%Q0(i)         = Q0(i)
            basis_temp%scaleQ(i)     = ONE
            basis_temp%opt_Q0(i)     = opt_Q0(i)
            basis_temp%opt_scaleQ(i) = opt_scaleQ(i)
          ELSE
            basis_temp%A(i)      = ZERO
            basis_temp%B(i)      = ZERO
            basis_temp%Q0(i)     = ZERO
            basis_temp%scaleQ(i) = ONE
          END IF
        END DO

        basis_temp%auto_basis = (auto_basis .AND. ndim == 1)
        IF (auto_basis .AND. ndim == 1 .AND. associated(mole%NMTransfo)) THEN
          IF (.NOT. mole%tab_Qtransfo(mole%itNM)%skip_transfo) THEN
          IF (associated(mole%NMTransfo%Q0_HObasis) .AND. associated(mole%NMTransfo%scaleQ_HObasis)) THEN

            basis_temp%Q0(1)     = mole%NMTransfo%Q0_HObasis(basis_temp%iQdyn(1))
            basis_temp%scaleQ(1) = mole%NMTransfo%scaleQ_HObasis(basis_temp%iQdyn(1))
            write(out_unitp,*) 'Q0,scaleQ (from auto_basis): ',basis_temp%Q0,basis_temp%scaleQ
            basis_temp%auto_basis = .FALSE.
          ELSE
            write(out_unitp,*) ' WARNING in ',name_sub
            write(out_unitp,*) '  auto_basis=t and ...'
            write(out_unitp,*) '   NMTransfo%Q0_HObasis or NMTransfo%scaleQ_HObasis are not associated'
            write(out_unitp,*) '      The Q0 and scaleQ are not modified !!'
            write(out_unitp,*) ' CHECK your data'
            !write(out_unitp,basis_nD)
            !STOP
          END IF
          END IF
        ELSE IF (auto_basis .AND. ndim == 1 .AND. .NOT. associated(mole%NMTransfo)) THEN
          write(out_unitp,*) ' WARNING in ',name_sub
          write(out_unitp,*) '  auto_basis=t and mole%NMTransfo is not associated'
          write(out_unitp,*) ' CHECK your data'
          !write(out_unitp,basis_nD)
          !STOP
        END IF

        basis_temp%opt_param = count(basis_temp%opt_A /= 0) +           &
          count(basis_temp%opt_B /= 0) + count(basis_temp%opt_Q0 /= 0) +&
                                       count(basis_temp%opt_scaleQ /= 0)
        IF (print_level > 1)                                            &
                  write(out_unitp,*) 'Parameter(s) to be optimized?: ', &
                                                   basis_temp%opt_param

      END IF


!---------------------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(basis_temp)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
      CALL flush_perso(out_unitp)
!---------------------------------------------------------------------
      end subroutine read5_basis_nD
!================================================================
!       analysis of a string
!       output : nb_word,word (the i th word of name)
!================================================================
      SUBROUTINE analysis_name(name,word,i,nb_word)
      USE mod_system
      IMPLICIT NONE


!     - analysis of the basis name --------------------------
      character (len=*)  :: word,name
      character          :: ch
      integer            :: i,nb_word
      integer            :: iw,ic,icw
      logical            :: blank


      iw = 0
      icw = 0
      blank = .TRUE.
!     write(out_unitp,*) 'analysis_name: ',name,len(name)
      DO ic=1,len(name)
        ch = name(ic:ic)
        IF (ch .EQ. " ") THEN
          IF (.NOT. blank) THEN
            iw = iw + 1
            blank = .TRUE.
          END IF
        ELSE
          IF (iw .EQ. i-1) THEN
            icw = icw + 1
            word(icw:icw) = ch
          END IF
          blank = .FALSE.
        END IF
!       write(out_unitp,*) 'analysis_name: ',ic,ch,blank,iw
      END DO

      nb_word = iw
!     write(out_unitp,*) 'analysis_name: ',name,':',nb_word
!     write(out_unitp,*) 'analysis_name: ',i,word


      end subroutine analysis_name
