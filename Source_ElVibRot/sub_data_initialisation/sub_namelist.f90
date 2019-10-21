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
! ++    read inactive parameters (for cHAC or HADA methods)
!
!       nb_quadrature : grid quadrature points
!       max_excit     : (or max_herm) max degree for the hermite polynomias
!       max_ene_h     : energy limit (in cm-1 for the imput then in au)
!       num           : .TRUE. if numerical derivatives
!       step          : step for numerical derivatives
!       auTOcm_inv    : conversion factor hartree (au) to cm-1
!
!================================================================
!
      SUBROUTINE read_inactive(Basis2n,ComOp,mole)
      USE mod_system
      USE mod_nDindex
      USE mod_Constant, only : REAL_WU, convRWU_TO_R
      use mod_Coord_KEO, only: zmatrix, alloc_array, dealloc_array, &
                               set_rphtransfo, tnum, alloc_nparray
      USE mod_basis
      USE mod_Op
      IMPLICIT NONE

!----- for the basis set ----------------------------------------------
      TYPE (basis)          :: Basis2n
      TYPE (param_ComOp)    :: ComOp
      TYPE (zmatrix), intent(inout) :: mole

!-----------------------------------------------------------

!----- working variables -----------------------------------

      integer           :: nb_quadrature,max_excit,max_coupling,n_h
      real (kind=Rkind) :: step
      logical           :: num,ADA
      TYPE (REAL_WU)    :: max_ene_h



      integer, parameter :: max_inact2n = 500
      integer       :: i,k,nb_inact21
      integer       :: tab_nq(max_inact2n)
      integer       :: tab_nb(max_inact2n)
      integer       :: nDinit(max_inact2n)
      logical       :: SparseGrid
      integer       :: isort,L_SparseGrid

      logical       :: H0_sym,gradTOpot0,diabatic_freq
      integer       :: Qinact2n_sym(max_inact2n)
      integer       :: Qinact2n_eq(max_inact2n,max_inact2n)

      logical       :: non_adia,contrac_ba_ON_HAC
      integer       :: max_nb_ba_ON_HAC

      NAMELIST /inactives/nb_quadrature,max_excit,max_coupling,         &
                          isort,                                        &
                          max_ene_h,num,step,n_h,                       &
                          ADA,tab_nq,tab_nb,                            &
                          SparseGrid,L_SparseGrid,                      &
                          H0_sym,diabatic_freq,Qinact2n_sym,Qinact2n_eq,&
                          gradTOpot0,                                   &
                          non_adia,contrac_ba_ON_HAC,max_nb_ba_ON_HAC


!------- read the inactive namelist ----------------------------
!      logical, parameter :: debug=.TRUE.
      logical, parameter :: debug=.FALSE.
      character (len=*), parameter :: name_sub='read_inactive'
!      -----------------------------------------------------------------
      IF(MPI_id==0) write(out_unitp,*) 'INACTIVES PARAMETERS'
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
      END IF


      ! get nb_inact21 from mole%nb_inact2n or mole%RPHTransfo%nb_inact21
      IF (associated(mole%RPHTransfo)) THEN
        nb_inact21 = mole%RPHTransfo%nb_inact21
      ELSE
        nb_inact21 = mole%nb_inact2n
      END IF

      IF (nb_inact21 <= 0) RETURN
      IF (nb_inact21 > max_inact2n) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'max_inact2n is too small',max_inact2n
        write(out_unitp,*) 'It should be larger than "nb_inact21"',nb_inact21
        STOP
      END IF

      tab_nq(:) = 0
      tab_nb(:) = 0

      Qinact2n_sym(:)  = 0
      Qinact2n_eq(:,:) = 0
      H0_sym           = .FALSE.
      diabatic_freq    = .FALSE.
      gradTOpot0       = .FALSE.

      step             = ONETENTH**5
      num              = .TRUE. ! not used anymore
      ADA              = .FALSE.

      max_ene_h        = REAL_WU(huge(ONE),   'au','E') ! huge (ua)
      nb_quadrature    = 10
      max_excit        = -1
      max_coupling     = mole%nb_inact2n
      n_h              = -1 ! all channels
      SparseGrid       = .FALSE.
      L_SparseGrid     = -1
      isort            = 1 ! 1: sort the HADA basis (energy)
                           ! 2: sort the HADA basis in term of excitation


      contrac_ba_ON_HAC = .FALSE.
      max_nb_ba_ON_HAC  = huge(1)

      read(in_unitp,inactives)

      IF (step < epsilon(ONE)*TEN**3) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'step is too small',step
        write(out_unitp,*) 'It should be larger than',epsilon(ONE)*TEN**3
        STOP
      END IF

      IF (.NOT. associated(mole%RPHTransfo)) THEN

        IF (product(tab_nq(1:nb_inact21)) == 0)  tab_nq(:) = nb_quadrature

        IF (product(tab_nb(1:nb_inact21)) == 0 ) tab_nb(:) = max_excit + 1
        ! here tab_nb can be 0 if max_excit=-1, therefore, we test again
        IF (product(tab_nb(1:nb_inact21)) == 0 ) tab_nb(:) = 1


        IF (max_excit == -1) THEN
          max_excit = sum(tab_nb(1:nb_inact21))-nb_inact21
        END IF

        write(6,*) 'max_excit',max_excit
        write(6,*) 'tab_nb(:)',tab_nb(1:nb_inact21)

       nDinit(:) = 0
       CALL alloc_array(Basis2n%nDindB,'Basis2n%nDindB',name_sub)
       IF (isort == 1) THEN ! sort with energy
         CALL init_nDindexPrim(Basis2n%nDindB,nb_inact21,               &
             tab_nb(1:nb_inact21),nDinit(1:nb_inact21),nb_OF_MinNorm=0, &
                     MaxNorm=convRWU_TO_R(max_ene_h),type_OF_nDindex=0, &
                                   With_nDindex=.FALSE.)
       ELSE IF (isort == 2) THEN ! sort with excitation
         CALL init_nDindexPrim(Basis2n%nDindB,nb_inact21,               &
             tab_nb(1:nb_inact21),nDinit(1:nb_inact21),nb_OF_MinNorm=0, &
                                      Lmax=max_excit,type_OF_nDindex=0, &
                                      With_nDindex=.FALSE.)
       ELSE
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) '   the isort value ',isort,' cannot be used.'
         write(out_unitp,*) '   check your data!!'
         STOP
       END IF
       !CALL Write_nDindex(Basis2n%nDindB)

       CALL init_nDindexPrim(Basis2n%nDindG,nb_inact21,                 &
             tab_nq(1:nb_inact21),type_OF_nDindex=0,With_nDindex=.FALSE.)
       CALL Write_nDindex(Basis2n%nDindG)

        ComOp%ADA               = ADA
        ComOp%contrac_ba_ON_HAC = contrac_ba_ON_HAC
        ComOp%max_nb_ba_ON_HAC  = max_nb_ba_ON_HAC
        ComOp%max_ene_ON_HAC    = convRWU_TO_R(max_ene_h)


        Basis2n%nb_basis            = nb_inact21
        IF (SparseGrid) THEN
          Basis2n%SparseGrid_type   = 3
        ELSE
          Basis2n%SparseGrid_type   = 0
        END IF
        Basis2n%L_SparseGrid        = L_SparseGrid

        CALL alloc_tab_Pbasis_OF_basis(Basis2n)

        IF (SparseGrid .AND. isort /= 2) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' With SparseGrid, you MUST use isort=2'
          write(out_unitp,*) '   SparseGrid,isort: ',SparseGrid,isort
          write(out_unitp,*) 'Check your data!'
          STOP
        END IF

        IF (Basis2n%SparseGrid_type == 3 .AND. Basis2n%L_SparseGrid < 1) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' SparseGrid_type /= 3 and L_SparseGrid < 1'
          write(out_unitp,*) ' You should increase L_SparseGrid!'
          STOP
        END IF
      END IF

      IF (nb_inact21 > 0 .AND. .NOT. associated(mole%RPHTransfo)) THEN
        ! here it is for the HADA/cHAC models (not RPH)
        IF (associated(mole%RPHTransfo_inact2n)) THEN
          CALL dealloc_array(mole%RPHTransfo_inact2n,                   &
                            'mole%RPHTransfo_inact2n',name_sub)
        END IF

        CALL alloc_array(mole%RPHTransfo_inact2n,                       &
                        'mole%RPHTransfo_inact2n',name_sub)

        CALL Set_RPHTransfo(mole%RPHTransfo_inact2n,                    &
                                mole%ActiveTransfo%list_act_OF_Qdyn,    &
                                       gradTOpot0,diabatic_freq,step,   &
                             H0_sym,H0_sym,Qinact2n_sym(1:nb_inact21),  &
                             Qinact2n_eq(1:nb_inact21,1:nb_inact21))

      ELSE IF (nb_inact21 > 0 .AND. associated(mole%RPHTransfo)) THEN
        IF (mole%RPHTransfo%option == 0) THEN
          ! here we set up not list_act_OF_Qdyn because it is already
          !     done in type_var_analysis

          CALL Set_RPHTransfo(mole%RPHTransfo,                          &
            gradTOpot0=gradTOpot0,diabatic_freq=diabatic_freq,step=step,&
                                    purify_hess=H0_sym,eq_hess=H0_sym,  &
                              Qinact2n_sym=Qinact2n_sym(1:nb_inact21),  &
                     Qinact2n_eq=Qinact2n_eq(1:nb_inact21,1:nb_inact21))
        END IF
      END IF


      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF

      END SUBROUTINE read_inactive
!
!================================================================
! ++    read the active parameters
!
!       max_ene1       : energy limit (in cm-1 for the imput then in au)
!       test           : .TRUE.    if a test on ONE active grid point
!       auTOcm_inv     : conversion factor hartree (au) to cm-1
!
!================================================================
!
      SUBROUTINE read_active(para_Tnum,mole,para_AllBasis,              &
                             ComOp,para_ReadOp,para_PES)

      USE mod_system
      USE mod_nDindex
      USE mod_Constant, only : REAL_WU,convRWU_TO_R
      USE mod_PrimOp
      USE mod_Op
      USE mod_basis
      USE mod_Auto_Basis
      IMPLICIT NONE

!----- for the zmatrix and Tnum --------------------------------------
      TYPE (Tnum)    :: para_Tnum
      TYPE (zmatrix) :: mole

!----- for the basis set ----------------------------------------------
      TYPE (param_AllBasis) :: para_AllBasis

!----- variables for the construction of H ---------------------------
      TYPE (param_ComOp)  :: ComOp
      TYPE (param_ReadOp) :: para_ReadOp
      TYPE (param_PES)    :: para_PES


      logical       :: test
      logical       :: comput_S

      logical       :: lect,restart
      logical       :: read_Grid,restart_Grid
      logical       :: Save_FileGrid,Keep_FileGrid,Save_MemGrid,Seq_FileGrid

      integer       :: num_grid_i,num_grid_f
      integer       :: JJ,Type_HamilOp

      logical       :: pot_only,T_only,pack_Op,read_Op,make_MatOp,direct_KEO,direct_ScalOp

      character (len=Name_len) :: name0

      integer           :: direct,Type_FileGrid
      real (kind=Rkind) :: tol_pack,tol_nopack

      logical           :: Op_Transfo
      TYPE (REAL_WU)    :: E0_Transfo


!     - working variables -----
      integer       :: ib,ind_b(mole%nb_act1+1)
      integer       :: i,j,k,i_term,nb_bi


!-------variables for the file names -------------------------------------
       character (len=Line_len) :: name_HADA
       logical                  :: formatted_HADA
       character (len=Line_len) :: name_Grid
       logical                  :: formatted_Grid
!-------------------------------------------------------------------------

      NAMELIST /actives/comput_S,test,                                  &
                        read_Grid,lect,restart_Grid,restart,            &
                        Save_FileGrid,Keep_FileGrid,Save_MemGrid,       &
                        Seq_FileGrid,Type_FileGrid,                     &
                        num_grid_i,num_grid_f,                          &
                        pot_only,T_only,read_Op,                        &
                        name_HADA,formatted_HADA,                       &
                        name_Grid,formatted_Grid,                       &
                        JJ,Type_HamilOp,direct_KEO,direct_ScalOp,       &
                        direct,make_MatOp,pack_Op,tol_pack,tol_nopack,  &
                        Op_Transfo,E0_Transfo


!------- test on max_HADA and n_h ---------------------------------
      write(out_unitp,*) ' ACTIVES PARAMETERS'


      IF (print_level > 0) write(out_unitp,*) 'BEGINNING read_active'
      IF (print_level > 0) write(out_unitp,*) 'nb_act1',mole%nb_act1


!------- read the active namelist ----------------------------

      make_MatOp      = .FALSE.
      tol_pack        = ONETENTH**7
      tol_nopack      = NINE/TEN
      pack_Op         = .FALSE.
      read_Op         = .FALSE.
      test            = .TRUE.
      comput_S        = .FALSE.
      Type_HamilOp    = 1
      direct_KEO      = .FALSE.
      direct_ScalOp   = .FALSE.

      direct          = 0
      Type_FileGrid   = 0
      Read_Grid       = .FALSE.
      Save_FileGrid   = .FALSE.
      Save_MemGrid    = .FALSE.
      Keep_FileGrid   = .FALSE.
      Seq_FileGrid    = .FALSE.
      lect            = .FALSE.
      Restart_Grid    = .FALSE.
      restart         = .FALSE.
      num_grid_i      = 0
      num_grid_f      = 0

      JJ              = -1

      name_HADA     = 'SH_HADA'
      formatted_HADA= .TRUE.
      name_Grid     = 'SH_HADA'
      formatted_Grid= .TRUE.


      pot_only      = .FALSE.
      T_only        = .FALSE.

      Op_Transfo    = .FALSE.
      E0_Transfo    = REAL_WU(ZERO,   'cm-1','E') !ZERO with cm-1 as default unit

      read(in_unitp,actives)
      IF (direct == 0 .OR. read_Op) make_MatOp = .TRUE.
      IF (print_level > 1) write(out_unitp,actives)


      IF (name_HADA /= 'SH_HADA' .OR. .NOT. formatted_HADA) THEN
        write(out_unitp,*) 'ERROR in read_active'
        write(out_unitp,*) '  You should not use the parameters (name_HADA, formatted_HADA) '
        write(out_unitp,*) '    in read_active namelist, '
        write(out_unitp,*) '  name_Grid and formatted_Grid'
        STOP
        name_Grid      = name_HADA
        formatted_Grid = formatted_HADA
      END IF

      IF (lect .OR. restart) THEN
        write(out_unitp,*) 'ERROR in read_active'
        write(out_unitp,*) '  You should not use the parameters (lect or restart) '
        write(out_unitp,*) '    in read_active namelist, '
        write(out_unitp,*) '  instead, you must use: Read_Grid and Restart_Grid'
        STOP
      END IF

      IF (JJ /= -1) THEN
        para_Tnum%JJ                = JJ
        para_Tnum%With_Cart_Transfo = (JJ>0) .AND. mole%Cart_transfo
      END IF

!     - for ComOp ----------------------------
      nb_bi = get_nb_bi_FROM_AllBasis(para_AllBasis)
      CALL All2_param_TO_ComOp(ComOp,para_AllBasis,mole,nb_bi,          &
                               para_PES%nb_elec,name_Grid,formatted_Grid)

      para_PES%Type_HamilOp  = Type_HamilOp
      para_PES%direct_KEO    = direct_KEO
      para_PES%direct_ScalOp = direct_ScalOp



      IF (print_level > 0) write(out_unitp,*) ' END read_active'

!     - Copy parameters in para_ReadOp --------

      SELECT CASE (direct)
      !         direct=0    => Make_Mat=T, SaveFile_Grid=T, SaveMem_Grid=F
      !         direct=1    => Make_Mat=F, SaveFile_Grid=T, SaveMem_Grid=T
      !         direct=2    => Make_Mat=F, SaveFile_Grid=F, SaveMem_Grid=T
      !         direct=3    => Make_Mat=F, SaveFile_Grid=T, SaveMem_Grid=F (for huge grid, like cHAC)
      CASE Default ! id case(0)
        para_ReadOp%Make_Mat = .TRUE.
        CALL init_FileGrid(para_ReadOp%para_FileGrid,                   &
                           Type_FileGrid=0,Save_FileGrid =.TRUE.,       &
                           Keep_FileGrid=.TRUE.,Save_MemGrid=.FALSE.,   &
                           Read_FileGrid=Read_Grid,                     &
                           Restart_Grid=Restart_Grid,test_Grid=test,    &
                           First_GridPoint=num_grid_i,                  &
                           Last_GridPoint=num_grid_f,                   &
                           Base_FileName_Grid=name_Grid)
      CASE (0)
        para_ReadOp%Make_Mat = .TRUE.
        CALL init_FileGrid(para_ReadOp%para_FileGrid,                   &
                           Type_FileGrid=0,Save_FileGrid =.TRUE.,       &
                           Formatted_FileGrid=formatted_Grid,           &
                           Keep_FileGrid=.TRUE.,Save_MemGrid=.FALSE.,   &
                           Read_FileGrid=Read_Grid,                     &
                           Restart_Grid=Restart_Grid,test_Grid=test,    &
                           First_GridPoint=num_grid_i,                  &
                           Last_GridPoint=num_grid_f,                   &
                           Base_FileName_Grid=name_Grid)
      CASE (1)
        para_ReadOp%Make_Mat = .FALSE.
        CALL init_FileGrid(para_ReadOp%para_FileGrid,                   &
                           Type_FileGrid=0,Save_FileGrid =.TRUE.,       &
                           Formatted_FileGrid=formatted_Grid,           &
                           Keep_FileGrid=.TRUE.,Save_MemGrid=.TRUE.,    &
                           Read_FileGrid=Read_Grid,                     &
                           Restart_Grid=Restart_Grid,test_Grid=test,    &
                           First_GridPoint=num_grid_i,                  &
                           Last_GridPoint=num_grid_f,                   &
                           Base_FileName_Grid=name_Grid)
      CASE (2)
        para_ReadOp%Make_Mat = .FALSE.
        CALL init_FileGrid(para_ReadOp%para_FileGrid,                   &
                           Type_FileGrid=0,Save_FileGrid =.FALSE.,      &
                           Formatted_FileGrid=formatted_Grid,           &
                           Keep_FileGrid=.FALSE.,Save_MemGrid=.TRUE.,   &
                           Read_FileGrid=Read_Grid,                     &
                           Restart_Grid=Restart_Grid,test_Grid=test,    &
                           First_GridPoint=num_grid_i,                  &
                           Last_GridPoint=num_grid_f,                   &
                           Base_FileName_Grid=name_Grid)

        IF (Type_FileGrid /= 0) THEN
          para_ReadOp%para_FileGrid%Save_MemGrid       = Save_MemGrid
          para_ReadOp%para_FileGrid%Type_FileGrid      = Type_FileGrid
          para_ReadOp%para_FileGrid%Keep_FileGrid      = Keep_FileGrid
          para_ReadOp%para_FileGrid%Formatted_FileGrid = .FALSE.
          para_ReadOp%para_FileGrid%Save_FileGrid      = .TRUE.
        END IF

      CASE (3)
        para_ReadOp%Make_Mat = .FALSE.
        CALL init_FileGrid(para_ReadOp%para_FileGrid,                   &
                           Type_FileGrid=0,Save_FileGrid =.TRUE.,       &
                           Formatted_FileGrid=formatted_Grid,           &
                           Keep_FileGrid=.TRUE.,Save_MemGrid=.FALSE.,   &
                           Read_FileGrid=Read_Grid,                     &
                           Restart_Grid=Restart_Grid,test_Grid=test,    &
                           First_GridPoint=num_grid_i,                  &
                           Last_GridPoint=num_grid_f,                   &
                           Base_FileName_Grid=name_Grid)

      CASE (4) ! for SG4, it enable to use directKEO and also direct pot
        para_ReadOp%Make_Mat = .FALSE.
        CALL init_FileGrid(para_ReadOp%para_FileGrid,                   &
                           Type_FileGrid=4,Save_FileGrid=Save_FileGrid, &
                           Formatted_FileGrid=.FALSE.,                  &
                           Keep_FileGrid=Keep_FileGrid,                 &
                           Save_MemGrid=Save_MemGrid,                   &
                           Read_FileGrid=Read_Grid,                     &
                           Restart_Grid=Restart_Grid,test_Grid=test,    &
                           First_GridPoint=num_grid_i,                  &
                           Last_GridPoint=num_grid_f,                   &
                           Base_FileName_Grid=name_Grid)
      END SELECT

      !CALL Write_FileGrid(para_ReadOp%para_FileGrid)

      IF (make_MatOp) para_ReadOp%Make_Mat = .TRUE.

      IF (direct < 2) THEN
        para_ReadOp%comput_S = comput_S
      ELSE
        para_ReadOp%comput_S = .FALSE.
      END IF

      para_ReadOp%pack_Op         = pack_Op .AND. para_ReadOp%Make_Mat
      para_ReadOp%read_Op         = read_Op
      para_ReadOp%tol_pack        = tol_pack
      para_ReadOp%tol_nopack      = tol_nopack

      para_ReadOp%pot_only        = pot_only
      para_ReadOp%T_only          = T_only

      para_ReadOp%Op_Transfo      = Op_Transfo
      IF(MPI_id==0) write(out_unitp,*) 'para_ReadOp%Op_Transfo',para_ReadOp%Op_Transfo
      IF (Op_Transfo) THEN
        !write(out_unitp,*) 'E0_Transfo',E0_Transfo

        para_ReadOp%E0_Transfo      = convRWU_TO_R(E0_Transfo)
        para_ReadOp%degree_Transfo  = 2
        CALL alloc_NParray(para_ReadOp%Poly_Transfo,(/ 2 /),            &
                          'para_ReadOp%Poly_Transfo','read_active',(/ 0 /))
        para_ReadOp%Poly_Transfo = (/ ZERO,ZERO,ONE /) ! x^2
        write(out_unitp,*) 'para_ReadOp%E0_Transfo',para_ReadOp%E0_Transfo
        !write(out_unitp,*) 'l u bounds',ubound(para_ReadOp%Poly_Transfo),lbound(para_ReadOp%Poly_Transfo)

      END IF

      IF (JJ > -1) THEN
        para_ReadOp%nb_bRot         = 2*JJ+1
      ELSE
        para_ReadOp%nb_bRot         = 1
      END IF
      write(out_unitp,*) 'The number of rotational basis is:',para_ReadOp%nb_bRot

      END SUBROUTINE read_active
