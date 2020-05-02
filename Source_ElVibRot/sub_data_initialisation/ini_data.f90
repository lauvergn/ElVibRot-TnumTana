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


!=======================================================================================
!     INITIALIZATION of most of parameters :
!        mole, Operators....
!
! >const_phys: type constant contains all the physical constants
!              *Note: may not safe with this name
!              setup with calling 'sub_constantes', 
!              *Note: a namelist 'constants' is required in input file for this 
!                     subroutine if call with sub_constantes(const_phys,.TRUE.)
! >para_OTF: type param_OTF
! >para_Tnum: type Tnum, which contains type param_PES_FromTnum and sum_opnd
! >mole: type CoordType
! >para_AllBasis: type param_AllBasis, contians four Basis pointers 
!                 BasisnD,Basis2n,BasisElec,BasisRot
!                 setup with calling alloc_AllBasis(para_AllBasis)
! >para_AllOp: type param_AllOp, for the construction of H
!=======================================================================================
      SUBROUTINE ini_data(const_phys,para_Tnum,mole,                    &
                          para_AllBasis,para_AllOp,                     &
                          para_ana,para_intensity,intensity_only,       &
                          para_propa)

      use mod_system
      USE mod_dnSVM,     only : Type_dnMat
      USE mod_Constant,  only : constant, sub_constantes, REAL_WU
      USE mod_Coord_KEO, only : CoordType, Tnum, get_Qact0, read_RefGeom
      USE mod_PrimOp,    only : write_typeop, param_typeop,init_typeop, &
                                Finalize_tnumtana_coord_primop,         &
                                derive_termqact_to_derive_termqdyn,     &
                                param_otf, PrimOp_t, write_PrimOp
      USE mod_basis
      USE mod_Set_paraRPH
      USE mod_Op
      USE mod_analysis
      USE mod_propa
      USE mod_Auto_Basis
      USE mod_MPI
      IMPLICIT NONE

!----- physical and mathematical constants ----------------------------
      TYPE (constant)        :: const_phys

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole
      TYPE (Tnum)      :: para_Tnum

!----- variables for the active and inactive namelists ----------------

!----- for the basis set ----------------------------------------------
      TYPE (param_AllBasis) :: para_AllBasis

!----- variables pour la namelist analyse ----------------------------
      TYPE (param_ana)       :: para_ana
      TYPE (param_intensity) :: para_intensity
      logical                :: intensity_only

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa)     :: para_propa
      integer                :: nb_vp_specWP

!----- Operators -------------------------------------------
      TYPE (param_AllOp)     :: para_AllOp

!----- working variables ---------------------------------------------


!----- Parameters to set the operators --------------------------------
      TYPE (param_ReadOp)            :: para_ReadOp

      integer                        :: i,j,rk,rl,i_term,iOp,it
      integer                        :: nq,nb_ba,nb_bi,nb_be,nb_bie
      logical                        :: spectral_H
      real (kind=Rkind), allocatable :: Qana(:),Qact(:)

!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'ini_data'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
      END IF
!---------------------------------------------------------------------

!=====================================================================
!=====================================================================
!=====================================================================

!---------------------------------------------------------------------
!------ read or set up the physical constants ------------------------
      CALL sub_constantes(const_phys,.TRUE.)

      CALL flush_perso(out_unitp)
      IF(MPI_id==0) THEN
        write(out_unitp,*) "============================================================"
        write(out_unitp,*) "============================================================"
        write(out_unitp,*) "=== COORDINATES (TNUM) ====================================="
        write(out_unitp,*) "============================================================"
        write(out_unitp,*) "============================================================"
      ENDIF 
      CALL flush_perso(out_unitp)

!---------------------------------------------------------------------
!------- read the coordinates ....     -------------------------------
      CALL Read_CoordType(mole,para_Tnum,const_phys)

!---------------------------------------------------------------------
!----- Read the namelist:  minimum -----------------------------------
      CALL read_RefGeom(mole,para_Tnum)
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!--------------------- TO finalize the coordinates (NM) and the KEO ----
!     If needed, Tana must be done after auto_basis, otherwise nrho(:) could be wrong
      CALL Finalize_TnumTana_Coord_PrimOp(para_Tnum,mole,para_ReadOp%PrimOp_t,Tana=.FALSE.)
!-----------------------------------------------------------------------

      IF(MPI_id==0) THEN
        write(out_unitp,*) "============================================================"
        write(out_unitp,*) "============================================================"
        write(out_unitp,*) "=== END COORDINATES (TNUM) ================================="
        write(out_unitp,*) "============================================================"
        write(out_unitp,*) "============================================================"
        CALL flush_perso(out_unitp)
      ENDIF

!---------------------------------------------------------------------
!------- read basis -- -----------------------------------------------
!     ----------------------------------------------------------------
!     Calculation of the quadrature and weight points
!     and calculations of the basis (d0b) and their dervivatives (d1b d2b)
!     on the grid points.
!     ----------------------------------------------------------------
      IF(MPI_id==0) THEN 
        write(out_unitp,*) "============================================================"
        write(out_unitp,*) "============================================================"
        write(out_unitp,*) "=== BASIS =================================================="
        write(out_unitp,*) "============================================================"
        write(out_unitp,*) "============================================================"
        CALL flush_perso(out_unitp)
      ENDIF

      ! allocate para_AllBasis, but no big mem at this point
      CALL alloc_AllBasis(para_AllBasis)

      CALL flush_perso(out_unitp)
      IF(MPI_id==0) THEN
        write(out_unitp,*) '================================================='
        write(out_unitp,*) '== READ ACTIVE BASIS ============================'
        write(out_unitp,*) '================================================='
        CALL flush_perso(out_unitp)
      ENDIF

      ! basis information in BasisnD
      CALL read_basis5(para_AllBasis%BasisnD,mole)

      IF(MPI_id==0) THEN
        write(out_unitp,*) '================================================='
        write(out_unitp,*) '== END READ ACTIVE BASIS ========================'
        write(out_unitp,*) '================================================='

        write(out_unitp,*) '================================================='
        write(out_unitp,*) '== INACTIVE BASIS ==============================='
        write(out_unitp,*) '================================================='
        CALL flush_perso(out_unitp)
      ENDIF

      CALL read_inactive(para_AllBasis,mole)

      IF (mole%nb_inact2n == 0) THEN
        CALL dealloc_basis(para_AllBasis%Basis2n)
        para_AllBasis%Basis2n%nb = 1
      ELSE

        !------------------------------------------------------------
        ! List of basis functions
        CALL sub2_ind_harm(para_AllBasis%Basis2n,para_ReadOp%PrimOp_t,para_Tnum,mole)
        !------------------------------------------------------------

        !------------------------------------------------------------
        ! Grid points and basis functions (hermite) for inactive coordinates (2n)
        !------------------------------------------------------------
        IF (para_AllBasis%Basis2n%SparseGrid_type == 3) THEN
          ! For Sparse Basis and SparseGrid
          CALL sub_quadra_SparseBasis2n(para_AllBasis%Basis2n,mole)
        ELSE
          ! For direct product grid
          CALL sub_quadra_inact(para_AllBasis%Basis2n,mole)
        END IF
      END IF

      IF(MPI_id==0) THEN
        write(out_unitp,*) '================================================='
        write(out_unitp,*) '== END INACTIVE BASIS ==========================='
        write(out_unitp,*) '================================================='
        CALL flush_perso(out_unitp)
      ENDIF
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      IF(MPI_id==0) THEN 
        write(out_unitp,*) "============================================================"
        write(out_unitp,*) "============================================================"
        write(out_unitp,*) "=== END BASIS =============================================="
        write(out_unitp,*) "============================================================"
        write(out_unitp,*) "============================================================"
        CALL flush_perso(out_unitp)
      ENDIF

      CALL read_active(para_Tnum,mole,para_ReadOp)

!---------------------------------------------------------------------
!------- read the parameter to analyze wave functions ----------------
      CALL alloc_NParray(Qana,(/ mole%nb_var /),"Qana",name_sub)
      CALL get_Qact0(Qana,mole%ActiveTransfo)

      CALL read_analyse(para_ana,Qana)

      CALL dealloc_NParray(Qana,"Qana",name_sub)

      IF (para_ana%VibRot .AND. para_ana%JJmax <= 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) '  VibRot=t and JJmax<1'
        write(out_unitp,*) ' It is impossible, you have to:'
        write(out_unitp,*) '(i)     Set VibRot=f'
        write(out_unitp,*) '(ii) or Set VibRot=t and JJmax > 0'
        STOP
      END IF
      IF (para_ana%VibRot .AND. para_ana%JJmax > 0 .AND. Para_Tnum%JJ > 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) '  VibRot=t and ...'
        write(out_unitp,*) '  para_ana%JJmax > 0 and Para_Tnum%JJ > 0'
        STOP
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!     -- reading parameters for intensity ----------------------------
      para_intensity%Temp     = -ONE
      IF (para_ana%intensity) THEN
        CALL read_intensity(para_intensity)
        IF (para_intensity%l_IntVR .AND.para_Tnum%JJ == 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '   You have incompatible parameters in: '
          write(out_unitp,*) ' l_IntVR in the namelist "intensity"',para_intensity%l_IntVR
          write(out_unitp,*) ' JJ in the namelist "variable"',para_Tnum%JJ
          write(out_unitp,*) ' You should have: l_IntVR=f and JJ=0 '
          write(out_unitp,*) ' or : l_IntVR=t and JJ>0'
          STOP
        END IF
      END IF
      para_Tnum%Inertia = para_intensity%Inertia

      IF (para_ana%Temp < ZERO) para_ana%Temp = para_intensity%Temp
      IF (para_ana%Temp < ZERO) para_ana%Temp = 298_Rkind
!---------------------------------------------------------------------
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!     -- reading the WP propagation ----------------------------------
      IF (para_ana%propa) THEN
        para_propa%control = para_ana%control !< for controling the calcutation of Hmin
        para_propa%max_ana = para_ana%max_ana
        CALL read_propagation(para_propa,mole%nb_act1,nb_vp_specWP)
        spectral_H = (para_propa%type_WPpropa == 10)

        IF (para_propa%spectral) THEN
           para_ReadOp%spectral     = .TRUE.
           para_ana%Spectral_ScalOp = .TRUE.
           para_AllBasis%basis_ext%nb_vp_spec = nb_vp_specWP
           IF (nb_vp_specWP > 0) THEN
             CALL alloc_NParray(para_AllBasis%basis_ext%liste_spec,     &
                                                        [nb_vp_specWP], &
                               'para_AllBasis%basis_ext%liste_spec',name_sub)
             para_AllBasis%basis_ext%liste_spec(:) = [ (i,i=1,nb_vp_specWP) ]
           END IF
        END IF

        para_propa%para_WP0%WP0n_h =                                    &
                            min(get_nb_bi_FROM_AllBasis(para_AllBasis), &
                                para_propa%para_WP0%WP0n_h)

        IF (para_propa%Write_WPAdia) para_ana%ana_psi%adia = .TRUE.
        para_ana%ana_psi%file_Psi  = para_propa%file_WP

        para_propa%ana_psi = para_ana%ana_psi

      ELSE
        spectral_H = .FALSE.
      END IF

      IF (para_ana%davidson .OR. para_ana%arpack .OR. para_ana%filter) THEN

        IF (para_ana%davidson) THEN
           para_propa%name_WPpropa      = 'Davidson'
           para_propa%type_WPpropa      = 33
        ELSE IF (para_ana%arpack) THEN
           para_propa%name_WPpropa      = 'Arpack'
           para_propa%type_WPpropa      = 33
        ELSE IF (para_ana%filter) THEN
           para_propa%name_WPpropa      = 'Filter'
           para_propa%type_WPpropa      = 33
        END IF

        para_propa%file_WP%name                        = 'file_WP'
        para_propa%file_WP%formatted                   = para_ana%formatted_file_WP
        para_propa%para_Davidson%formatted_file_readWP = para_ana%formatted_file_WP ! to set the same default as in para_ana%formatted_file_WP
        para_propa%para_Davidson%formatted_file_WP     = para_ana%formatted_file_WP ! idem

        CALL read_Davidson(para_propa%para_Davidson,para_propa)

       para_propa%file_WP%name    = para_propa%para_Davidson%name_file_saveWP

       para_ana%ana_psi%file_Psi  = para_propa%file_WP

      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------

!---------------------------------------------------------------------
!---------------------------------------------------------------------
!     -- Some parameters for the operators ----------------------------
      para_ReadOp%calc_scalar_Op = para_propa%with_field .OR.           &
                                   para_ana%intensity .OR. para_ana%NLO

      IF (para_ReadOp%nb_scalar_Op > 0) para_ReadOp%calc_scalar_Op = .TRUE.

      IF ((para_propa%with_field .OR. para_ana%intensity) .AND.         &
             para_ReadOp%nb_scalar_Op == 0) para_ReadOp%nb_scalar_Op = 3

      IF (.NOT. para_ana%davidson .AND. .NOT. para_ana%arpack .AND.     &
          .NOT. para_ana%filter   .AND. para_ana%CRP <= 0 .AND.         &
          .NOT. para_ana%propa)             para_ReadOp%make_Mat = .TRUE.
!---------------------------------------------------------------------
!---------------------------------------------------------------------

!=====================================================================
!
!      End of data reading
!
!=====================================================================
      CALL flush_perso(out_unitp)
      IF(MPI_id==0) THEN
        write(out_unitp,*) "============================================================"
        write(out_unitp,*) "============================================================"
        write(out_unitp,*) "=== AUTO BASIS ============================================="
        write(out_unitp,*) "============================================================"
        write(out_unitp,*) "============================================================"
        CALL flush_perso(out_unitp)
      ENDIF

      CALL Auto_basis(para_Tnum,mole,para_AllBasis,para_ReadOp)

      IF (para_AllBasis%BasisnD%SparseGrid_type == 4) THEN
        para_AllBasis%BasisnD%para_SGType2%nb0 = para_ReadOp%nb_elec ! to be changed
      END IF

      nb_bie = get_nb_bi_FROM_AllBasis(para_AllBasis) * para_ReadOp%nb_elec
      nb_ba  = get_nb_FROM_Basis(para_AllBasis%BasisnD)
      CALL alloc_basis_ext2n(para_AllBasis%basis_ext2n,nb_ba,nb_bie)


      IF (para_Tnum%Tana) THEN
        CALL alloc_NParray(Qact,(/ mole%nb_var /),"Qact",name_sub)
        CALL get_Qact0(Qact,mole%ActiveTransfo)
        CALL nrho_Basis_TO_nhro_Tnum(para_AllBasis,mole)
        CALL compute_analytical_KEO(para_Tnum%TWOxKEO,mole,para_Tnum,Qact)
        CALL comparison_G_FROM_Tnum_Tana(para_Tnum%ExpandTWOxKEO,mole,para_Tnum,Qact)
        CALL dealloc_NParray(Qact,"Qact",name_sub)
      END IF

      IF(MPI_id==0) THEN
        !CALL RecWrite_basis(para_AllBasis%BasisnD,write_all=.TRUE.) ; stop
        write(out_unitp,*) "============================================================"
        write(out_unitp,*) "============================================================"
        IF (debug) THEN
          CALL RecWrite_basis(para_AllBasis%BasisnD)
          CALL write_basis_ext2n(para_AllBasis%basis_ext2n)
        END IF
        write(out_unitp,*) "============================================================"
        write(out_unitp,*) "============================================================"
        write(out_unitp,*) "=== END AUTO BASIS ========================================="
        write(out_unitp,*) "============================================================"
        write(out_unitp,*) "============================================================"
        CALL flush_perso(out_unitp)
      ENDIF

       !write(out_unitp,*) 'pack ?',para_AllBasis%BasisnD%packed_done
       !IF (para_AllBasis%BasisnD%packed_done) THEN
       !  CALL sub_MatOFdnSX_basis(para_AllBasis%BasisnD)
       !END IF
       !STOP

!=====================================================================
!=====================================================================
      nq    = get_nq_FROM_basis(para_AllBasis%BasisnD)
      nb_be = para_ReadOp%nb_elec
      nb_bi = get_nb_FROM_basis(para_AllBasis%Basis2n)
      nb_ba = get_nb_FROM_basis(para_AllBasis%BasisnD)
      CALL MemoryEstimation(nb_ba,nq,mole%nb_act,nb_bi*nb_be,           &
                            para_propa%para_Davidson%nb_WP)
!=====================================================================
!=====================================================================


!=====================================================================
!     initialization of AllOp,
!     i=1 => for H       : n_op = 0
!     i=2 => for S       : n_op = -1
!     i=3 => for Dipx    : n_op = 1 (or ScalOp1)
!     i=4 => for Dipy    : n_op = 2 (or ScalOp2)
!     i=5 => for Dipz    : n_op = 3 (or ScalOp3)
!     i=5 => for ScalOp4 : n_op = 4
!     ....
      IF(MPI_id==0) THEN
        write(out_unitp,*) "============================================================"
        write(out_unitp,*) "============================================================"
        write(out_unitp,*) "====== List of Operators ==================================="
        write(out_unitp,*)
        write(out_unitp,*) 'para_ReadOp%nb_scalar_Op : ',para_ReadOp%nb_scalar_Op
      ENDIF

      IF (para_ReadOp%nb_scalar_Op > 27) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) 'You have defined too many Operator'
        write(out_unitp,*) ' nb_scalar_Op must be < 28',para_ReadOp%nb_scalar_Op
        STOP
      END IF
      IF (para_ReadOp%calc_scalar_Op .AND. para_ReadOp%nb_scalar_Op < 1) THEN
        para_ReadOp%nb_scalar_Op = 3
        write(out_unitp,*) ' WARNING in ',name_sub
        write(out_unitp,*) 'calc_scalar_Op=t and nb_scalar_Op < 1'
        write(out_unitp,*) ' You MUST set nb_scalar_Op in the namelis "minimun"'
      END IF
      IF (para_ReadOp%calc_scalar_Op .AND. para_ReadOp%nb_scalar_Op < 3) THEN
        write(out_unitp,*) ' WARNING in ',name_sub
        write(out_unitp,*) 'calc_scalar_Op=t and nb_scalar_Op < 3'
      END IF
      para_ReadOp%calc_scalar_Op = (para_ReadOp%nb_scalar_Op > 0)

      para_AllOp%nb_Op = para_ReadOp%nb_scalar_Op + 2  ! for H and S

      IF (debug) write(out_unitp,*) 'para_AllOp%nb_Op        : ',para_AllOp%nb_Op

      CALL alloc_array(para_AllOp%tab_Op,(/ para_AllOp%nb_Op /),        &
                      'para_AllOp%tab_Op',name_sub)

      iOp = 1 ! => for H
      IF (para_ana%Read_zpe) THEN
        para_AllOp%tab_Op(iOp)%Set_ZPE = .TRUE.
        para_AllOp%tab_Op(iOp)%ZPE     = para_ana%Ezpe
      END IF
      IF (para_ana%VibRot) Para_Tnum%JJ = para_ana%JJmax

      CALL All_param_TO_para_H(para_Tnum,mole,para_AllBasis,            &
                               para_ReadOp,para_AllOp%tab_Op(iOp))

      IF (spectral_H) THEN
        para_AllOp%tab_Op(iOp)%spectral    = .TRUE.
        para_AllOp%tab_Op(iOp)%spectral_Op = para_AllOp%tab_Op(iOp)%n_Op
         para_AllOp%tab_Op(iOp)%para_ReadOp%make_Mat = .TRUE.
      END IF
      para_AllOp%tab_Op(iOp)%symab      = 0 ! totally symmetric
      IF (para_ana%VibRot) Para_Tnum%JJ = 0
      IF(MPI_id==0) THEN
        write(out_unitp,*) 'para_ReadOp%Make_Mat',para_ReadOp%Make_Mat
        write(out_unitp,*) 'para_H%Make_Mat',para_AllOp%tab_Op(iOp)%Make_Mat
      ENDIF

      IF (debug) CALL Write_TypeOp(para_AllOp%tab_Op(iOp)%param_TypeOp)

      iOp = 2 ! for S
      CALL param_Op1TOparam_Op2(para_AllOp%tab_Op(1),para_AllOp%tab_Op(iOp))
      para_AllOp%tab_Op(iOp)%name_Op     = 'S'
      para_AllOp%tab_Op(iOp)%n_Op        = -1

      CALL Init_TypeOp(para_AllOp%tab_Op(iOp)%param_TypeOp,             &
                       type_Op=0,nb_Qact=mole%nb_act1,cplx=.FALSE.,     &
                       JRot=Para_Tnum%JJ,direct_KEO=.FALSE.,direct_ScalOp=.FALSE.)
      CALL derive_termQact_TO_derive_termQdyn(                          &
                              para_AllOp%tab_Op(iOp)%derive_termQdyn,   &
                              para_AllOp%tab_Op(iOp)%derive_termQact,   &
                              mole%ActiveTransfo%list_QactTOQdyn)

      para_AllOp%tab_Op(iOp)%symab    = 0 ! totally symmetric
      para_AllOp%tab_Op(iOp)%spectral = para_ana%Spectral_ScalOp


      DO i=1,para_ReadOp%nb_scalar_Op  ! for scalar operators (Dip)
        iOp = iOp + 1
        CALL param_Op1TOparam_Op2(para_AllOp%tab_Op(2),                 &
                                  para_AllOp%tab_Op(iOp))
        para_AllOp%tab_Op(iOp)%n_Op    = i
        para_AllOp%tab_Op(iOp)%name_Op = 'OpScal' // int_TO_char(i)

        CALL Init_TypeOp(para_AllOp%tab_Op(iOp)%param_TypeOp,           &
                         type_Op=0,nb_Qact=mole%nb_act1,cplx=.FALSE.,   &
                         JRot=Para_Tnum%JJ,                             &
                         direct_KEO=.FALSE.,direct_ScalOp=para_ReadOp%direct_ScalOp)
        CALL derive_termQact_TO_derive_termQdyn(                        &
                              para_AllOp%tab_Op(iOp)%derive_termQdyn,   &
                              para_AllOp%tab_Op(iOp)%derive_termQact,   &
                              mole%ActiveTransfo%list_QactTOQdyn)

        para_AllOp%tab_Op(iOp)%symab    = -1  ! the symmetry is not used
        para_AllOp%tab_Op(iOp)%spectral = para_ana%Spectral_ScalOp
      END DO

      IF (para_ReadOp%nb_scalar_Op == 3) THEN ! dipole moment operators
        para_AllOp%tab_Op(3)%name_Op = 'Dipx'
        para_AllOp%tab_Op(4)%name_Op = 'Dipy'
        para_AllOp%tab_Op(5)%name_Op = 'Dipz'
      END IF

      IF(MPI_id==0) THEN
        DO i=1,para_AllOp%nb_Op
          write(out_unitp,*) i,'Operator name: ',trim(para_AllOp%tab_Op(i)%name_Op)
        END DO

        write(out_unitp,*)
        write(out_unitp,*) "============================================================"
        write(out_unitp,*) "============================================================"

        write(out_unitp,*) "============================================================"
        write(out_unitp,*) "============================================================"
        write(out_unitp,*) "====== Finalize RPH transfo (Tnum) and ... ================="
        write(out_unitp,*) "====== ... EneH0 of the basis sets ========================="
        write(out_unitp,*)
      ENDIF

      IF (associated(mole%RPHTransfo)) THEN
        CALL Set_paraPRH(mole,para_Tnum,para_AllBasis%BasisnD)
      END IF

      IF (para_propa%para_Davidson%NewVec_type == 4 .OR. para_ana%CRP > 0) THEN
        CALL RecSet_EneH0(para_Tnum,mole,para_AllBasis%BasisnD,para_ReadOp)
      END IF

      IF(MPI_id==0) THEN
        write(out_unitp,*)
        write(out_unitp,*) "============================================================"
        write(out_unitp,*) "============================================================"
      ENDIF

!=====================================================================
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'nb_Op',para_AllOp%nb_Op
        DO i=1,para_AllOp%nb_Op
          CALL write_param_Op(para_AllOp%tab_Op(i))

        END DO
        DO i=1,para_AllBasis%BasisnD%nb_basis
          write(out_unitp,*) 'basis',i
          CALL RecWrite_basis(para_AllBasis%BasisnD)
        END DO
        write(out_unitp,*) 'END ',name_sub
      END IF
      !---------------------------------------------------------------------
      END SUBROUTINE ini_data
!=======================================================================================

!===============================================================================
SUBROUTINE MemoryEstimation(nb,nq,nb_Q,nb_channels,nb_psi)
USE mod_system
USE mod_MPI_Aid
IMPLICIT NONE

integer, intent(in) :: nb,nq,nb_Q,nb_channels,nb_psi

integer :: GridMem,BasisMem
integer :: MappingSG4Meme,PotMem,KEO_type1_Mem,KEO_type10_Mem,psi_Mem
integer :: Mem,nb_psi_loc
character (len=2) :: MemUnit


nb_psi_loc = max(nb_psi,1)


GridMem  = nq*Rkind
BasisMem = nb*Rkind

MappingSG4Meme = nq*Ikind

psi_Mem        = BasisMem * nb_channels
PotMem         = GridMem  * nb_channels**2
KEO_type1_Mem  = PotMem   * (nb_Q+1)*(nb_Q+2)/2 ! F2+F1+vep
KEO_type10_Mem = GridMem  * (nb_Q**2 + 2) ! size of G + jac+rho
              ! We suppose, the KEO are the same on each channel !! Pb...

IF(MPI_id==0) THEN
  write(out_unitp,*) "============================================================"
  write(out_unitp,*) "============================================================"
  write(out_unitp,*) "====== Memory psi and H ===================================="
ENDIF

IF(MPI_id==0) write(out_unitp,*) "------------------------------------------------------------"
Mem = psi_Mem
CALL convertMem(Mem,MemUnit)
IF(MPI_id==0) write(out_unitp,'(a,i0,1x,a)') "One psi: ",Mem,MemUnit
Mem = psi_Mem * 4
CALL convertMem(Mem,MemUnit)
IF(MPI_id==0) write(out_unitp,'(a,i0,1x,a)') "One psi's with Davidson (num_resetH=1): ",Mem,MemUnit
IF(MPI_id==0) write(out_unitp,*) "------------------------------------------------------------"

IF (nb_psi_loc > 0) THEN
  Mem = psi_Mem * nb_psi_loc
  CALL convertMem(Mem,MemUnit)
  IF(MPI_id==0) write(out_unitp,'(i0,a,i0,1x,a)') nb_psi_loc," psi's: ",Mem,MemUnit
  Mem = psi_Mem * nb_psi_loc * 4
  CALL convertMem(Mem,MemUnit)
  IF(MPI_id==0) write(out_unitp,'(i0,a,i0,1x,a)') nb_psi_loc," psi's with Davidson (num_resetH=1): ",Mem,MemUnit
  IF(MPI_id==0) write(out_unitp,*) "------------------------------------------------------------"
END IF

IF(MPI_id==0) write(out_unitp,*) "====== Memory for Type 1 (F2+F1+Vep+V) ====================="
IF(MPI_id==0) write(out_unitp,*) "-SG4 Full direct --"
Mem = MappingSG4Meme
CALL convertMem(Mem,MemUnit)
IF(MPI_id==0) write(out_unitp,'(a,i0,1x,a)') "H memory (mapping):       ",Mem,MemUnit
IF(MPI_id==0) write(out_unitp,*) "-SG4 KEO direct --"
Mem = MappingSG4Meme+PotMem
CALL convertMem(Mem,MemUnit)
IF(MPI_id==0) write(out_unitp,'(a,i0,1x,a)') "H memory (mapping+V):     ",Mem,MemUnit
IF(MPI_id==0) write(out_unitp,*) "-SG4 --"
Mem = MappingSG4Meme+KEO_type1_Mem ! PotMem is not here because is already counted in KEO_type1_Mem
CALL convertMem(Mem,MemUnit)
IF(MPI_id==0) write(out_unitp,'(a,i0,1x,a)') "H memory (mapping+V+KEO): ",Mem,MemUnit

IF(MPI_id==0) write(out_unitp,*) "====== Memory for Type 10 (G+V) ============================"
IF(MPI_id==0) write(out_unitp,*) "-SG4 Full direct --"
Mem = MappingSG4Meme
CALL convertMem(Mem,MemUnit)
IF(MPI_id==0) write(out_unitp,'(a,i0,1x,a)') "H memory (mapping):       ",Mem,MemUnit
IF(MPI_id==0) write(out_unitp,*) "-SG4 KEO direct --"
Mem = MappingSG4Meme+PotMem
CALL convertMem(Mem,MemUnit)
IF(MPI_id==0) write(out_unitp,'(a,i0,1x,a)') "H memory (mapping+V):     ",Mem,MemUnit
IF(MPI_id==0) write(out_unitp,*) "-SG4 --"
Mem = MappingSG4Meme+PotMem+KEO_type10_Mem
CALL convertMem(Mem,MemUnit)
IF(MPI_id==0) write(out_unitp,'(a,i0,1x,a)') "H memory (mapping+V+KEO): ",Mem,MemUnit

IF(MPI_id==0) THEN
  write(out_unitp,*) "============================================================"
  write(out_unitp,*) "============================================================"
ENDIF

END SUBROUTINE MemoryEstimation
