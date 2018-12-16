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


!=====================================================================
!
!     INITIALIZATION of most of parameters :
!        mole, WP0, Operators....
!
!=====================================================================
SUBROUTINE ini_data(const_phys,                                         &
                          para_OTF,                                     &
                          para_Tnum,mole,                               &
                          para_AllBasis,BasisnD_Save,                   &
                          para_PES,ComOp,para_AllOp,                    &
                          para_ana,para_intensity,intensity_only,       &
                          para_propa,WP0)

      use mod_system,    only : rkind, out_unitp, flush_perso, alloc_nparray,&
                                dealloc_nparray, one, zero, alloc_array,     &
                                int_to_char
      USE mod_dnSVM,     only : Type_dnMat
      USE mod_Constant,  only : constant, sub_constantes, REAL_WU
      USE mod_Coord_KEO, only : zmatrix, Tnum, get_Qact0, read_RefGeom
      use mod_PrimOp,    only : param_otf, param_pes, write_typeop, param_typeop,&
                                finalyze_tnumtana_coord_primop, init_typeop,     &
                                derive_termqact_to_derive_termqdyn
      USE mod_basis
      USE mod_Op
      USE mod_analysis
      USE mod_propa
      USE mod_psi_set_alloc
      USE mod_Auto_Basis
      IMPLICIT NONE


!----- On the fly parameters (at this time for gaussian) -------------
      TYPE (param_OTF) :: para_OTF

!----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole
      TYPE (Tnum)    :: para_Tnum

!----- variables for the active and inactive namelists ----------------

!----- for the basis set ----------------------------------------------
      TYPE (param_AllBasis) :: para_AllBasis
      TYPE (basis) :: BasisnD_Save

!----- variables pour la namelist analyse ----------------------------
      TYPE (param_ana)       :: para_ana
      TYPE (param_intensity) :: para_intensity
      logical                :: intensity_only

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_psi)   :: WP0
      integer            :: nb_vp_specWP

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_PES)   :: para_PES

!----- variables for the construction of H ---------------------------
      TYPE (param_AllOp)  :: para_AllOp
      TYPE (param_ComOp)  :: ComOp
      TYPE (param_ReadOp) :: para_ReadOp

!----- physical and mathematical constants ----------------------------
      TYPE (constant)             :: const_phys


!----- working variables ---------------------------------------------
      integer       :: i,j,rk,rl,i_term,iOp,it
      real (kind=Rkind), allocatable :: Qana(:)


      TYPE(Type_dnMat) :: dnGG
      integer       :: nq


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
      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "=== COORDINATES (TNUM) ====================================="
      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "============================================================"
      CALL flush_perso(out_unitp)

!---------------------------------------------------------------------

!---------------------------------------------------------------------
!------- read the coordinates ....     -------------------------------
      CALL Read_mole(mole,para_Tnum,const_phys)

!---------------------------------------------------------------------

!---------------------------------------------------------------------
!----- Read the namelist:  minimum -----------------------------------
      CALL read_RefGeom(mole,para_Tnum)
!---------------------------------------------------------------------

!-----------------------------------------------------------------------
!--------------------- TO finalize the coordinates (NM) and the KEO ----
      CALL Finalyze_TnumTana_Coord_PrimOp(para_Tnum,mole,para_PES)
!-----------------------------------------------------------------------

      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "=== END COORDINATES (TNUM) ================================="
      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "============================================================"
      CALL flush_perso(out_unitp)



!---------------------------------------------------------------------
!------- read basis -- -----------------------------------------------
!     ----------------------------------------------------------------
!     Calculation of the quadrature and weight points
!     and calculations of the basis (d0b) and their dervivatives (d1b d2b)
!     on the grid points.
!     ----------------------------------------------------------------
      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "=== BASIS =================================================="
      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "============================================================"
      CALL flush_perso(out_unitp)

      CALL alloc_AllBasis(para_AllBasis)


      CALL flush_perso(out_unitp)
      write(out_unitp,*) '================================================='
      write(out_unitp,*) '== READ ACTIVE BASIS ============================'
      write(out_unitp,*) '================================================='
      CALL flush_perso(out_unitp)

      CALL read_basis5(para_AllBasis%BasisnD,mole)
      CALL basis2TObasis1(BasisnD_Save,para_AllBasis%BasisnD)


      CALL flush_perso(out_unitp)
      write(out_unitp,*) '================================================='
      write(out_unitp,*) '== END READ ACTIVE BASIS ========================'
      write(out_unitp,*) '================================================='

      write(out_unitp,*) '================================================='
      write(out_unitp,*) '== INACTIVE BASIS ==============================='
      write(out_unitp,*) '================================================='
      CALL flush_perso(out_unitp)

      CALL read_inactive(para_AllBasis%Basis2n,ComOp,mole)

      IF (mole%nb_inact2n == 0) THEN
        CALL dealloc_basis(para_AllBasis%Basis2n)
        para_AllBasis%Basis2n%nb = 1
      ELSE

        CALL sub2_ind_harm(para_AllBasis%Basis2n,para_PES,para_Tnum,mole)
!       ------------------------------------------------------------

!       ------------------------------------------------------------
!       Grid points and basis functions (hermite) for inactive coordinates (2n)
!       ------------------------------------------------------------
        IF (para_AllBasis%Basis2n%SparseGrid_type == 3) THEN
!         For Sparse Basis and SparseGrid
          CALL sub_quadra_SparseBasis2n(para_AllBasis%Basis2n,mole)
        ELSE
!         For direct product grid
          CALL sub_quadra_inact(para_AllBasis%Basis2n,mole)
        END IF
      END IF
      write(out_unitp,*) '================================================='
      write(out_unitp,*) '== END INACTIVE BASIS ==========================='
      write(out_unitp,*) '================================================='
      CALL flush_perso(out_unitp)

!---------------------------------------------------------------------

!---------------------------------------------------------------------
      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "=== END BASIS =============================================="
      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "============================================================"
      CALL flush_perso(out_unitp)

      CALL read_active(para_Tnum,mole, para_AllBasis,ComOp,para_ReadOp,para_PES)

!---------------------------------------------------------------------
!------- read the parameter to analyze wave functions ----------------
      CALL alloc_NParray(Qana,(/ mole%nb_act /),"Qana",name_sub)
      CALL get_Qact0(Qana,mole%ActiveTransfo)

      CALL read_analyse(para_ana,Qana)

      CALL dealloc_NParray(Qana,"Qana",name_sub)

      IF (para_ana%VibRot .AND. para_ana%JJmax <= 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) '  VibRot=t and JJmax<1'
        write(out_unitp,*) ' It impossible, you have to:'
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
      IF (para_ana%Temp < ZERO) para_ana%Temp = real(298,kind=Rkind)

!     -- reading the WP propagation ---------------------------------
      IF (para_ana%propa) THEN
        para_propa%control = para_ana%control
        para_propa%max_ana = para_ana%max_ana
        CALL read_propagation(para_propa,mole%nb_act1,nb_vp_specWP)

        IF (para_propa%spectral) THEN
           para_ReadOp%spectral     = .TRUE.
           para_ana%Spectral_ScalOp = .TRUE.
           ComOp%nb_vp_spec         = nb_vp_specWP
           IF (nb_vp_specWP > 0) THEN
             CALL alloc_array(ComOp%liste_spec,(/ nb_vp_specWP /),      &
                             'ComOp%liste_spec',name_sub)
             ComOp%liste_spec(:) = (/ (i,i=1,nb_vp_specWP) /)
           END IF
        END IF

        para_propa%para_WP0%WP0n_h =                                    &
                            min(get_nb_bi_FROM_AllBasis(para_AllBasis), &
                                para_propa%para_WP0%WP0n_h)

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

      END IF

      para_PES%calc_scalar_Op = para_propa%with_field .OR.              &
                                    para_ana%intensity .OR. para_ana%NLO

      IF (para_PES%nb_scalar_Op > 0) para_PES%calc_scalar_Op = .TRUE.

      IF ((para_propa%with_field .OR. para_ana%intensity) .AND.         &
                   para_PES%nb_scalar_Op == 0) para_PES%nb_scalar_Op = 3

      IF (.NOT. para_ana%davidson .AND. .NOT. para_ana%arpack .AND.     &
          .NOT. para_ana%filter .AND. .NOT. para_ana%propa) para_ReadOp%make_Mat = .TRUE.

!---------------------------------------------------------------------

!=====================================================================
!
!      End of data reading
!
!=====================================================================
      CALL flush_perso(out_unitp)
      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "=== AUTO BASIS ============================================="
      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "============================================================"
      CALL flush_perso(out_unitp)

      CALL Auto_basis(para_Tnum,mole,para_AllBasis,ComOp,para_PES,para_ReadOp)


      !CALL RecWrite_basis(para_AllBasis%BasisnD,write_all=.TRUE.) ; stop
      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "============================================================"
      IF (debug) CALL RecWrite_basis(para_AllBasis%BasisnD)
      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "=== END AUTO BASIS ========================================="
      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "============================================================"
      CALL flush_perso(out_unitp)

       !write(out_unitp,*) 'pack ?',para_AllBasis%BasisnD%packed_done
       !IF (para_AllBasis%BasisnD%packed_done) THEN
       !  CALL sub_MatOFdnSX_basis(para_AllBasis%BasisnD)
       !END IF
       !STOP


!=====================================================================
!     initialization of AllOp,
!     i=1 => for H       : n_op = 0
!     i=2 => for S       : n_op = -1
!     i=3 => for Dipx    : n_op = 1 (or ScalOp1)
!     i=4 => for Dipy    : n_op = 2 (or ScalOp2)
!     i=5 => for Dipz    : n_op = 3 (or ScalOp3)
!     i=5 => for ScalOp4 : n_op = 4
!     ....
      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "====== List of Operators ==================================="
      write(out_unitp,*)
      write(out_unitp,*) 'para_PES%nb_scalar_Op : ',para_PES%nb_scalar_Op

      IF (debug) write(out_unitp,*) 'para_PES%nb_scalar_Op : ',para_PES%nb_scalar_Op
      IF (para_PES%nb_scalar_Op > 27) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) 'You have defined too many Operator'
        write(out_unitp,*) ' nb_scalar_Op must be < 28',para_PES%nb_scalar_Op
        STOP
      END IF
      IF (para_PES%calc_scalar_Op .AND. para_PES%nb_scalar_Op < 1) THEN
        para_PES%nb_scalar_Op = 3
        write(out_unitp,*) ' WARNING in ',name_sub
        write(out_unitp,*) 'calc_scalar_Op=t and nb_scalar_Op < 1'
        write(out_unitp,*) ' You MUST set nb_scalar_Op in the namelis "minimun"'
      END IF
      IF (para_PES%calc_scalar_Op .AND. para_PES%nb_scalar_Op < 3) THEN
        write(out_unitp,*) ' WARNING in ',name_sub
        write(out_unitp,*) 'calc_scalar_Op=t and nb_scalar_Op < 3'
      END IF

      para_AllOp%nb_Op = para_PES%nb_scalar_Op + 2  ! for H and S

      IF (debug) write(out_unitp,*) 'para_AllOp%nb_Op        : ',para_AllOp%nb_Op

      CALL alloc_array(para_AllOp%tab_Op,(/ para_AllOp%nb_Op /),        &
                      'para_AllOp%tab_Op',name_sub)


      iOp = 1 ! => for H
      IF (para_ana%Read_zpe) THEN
        ComOp%Set_ZPE = .TRUE.
        ComOp%ZPE     = para_ana%Ezpe
      END IF
      IF (para_ana%VibRot) Para_Tnum%JJ = para_ana%JJmax

      CALL All_param_TO_para_H(para_Tnum,mole,                          &
                               para_AllBasis,                           &
                               ComOp,para_PES,para_ReadOp,              &
                               para_AllOp%tab_Op(iOp))
      para_AllOp%tab_Op(iOp)%symab      = 0 ! totally symmetric
      IF (para_ana%VibRot) Para_Tnum%JJ = 0

      IF (debug) CALL Write_TypeOp(para_AllOp%tab_Op(iOp)%param_TypeOp)

      iOp = 2 ! for S
      CALL param_Op1TOparam_Op2(para_AllOp%tab_Op(1),para_AllOp%tab_Op(iOp))
      para_AllOp%tab_Op(iOp)%name_Op     = 'S'
      para_AllOp%tab_Op(iOp)%n_Op        = -1

      CALL Init_TypeOp(para_AllOp%tab_Op(iOp)%param_TypeOp,             &
                       type_Op=0,nb_Qact=mole%nb_act1,cplx=.FALSE.,     &
                       JRot=Para_Tnum%JJ,direct_KEO=.FALSE.)
      CALL derive_termQact_TO_derive_termQdyn(                          &
                              para_AllOp%tab_Op(iOp)%derive_termQdyn,   &
                              para_AllOp%tab_Op(iOp)%derive_termQact,   &
                              mole%ActiveTransfo%list_QactTOQdyn)

      para_AllOp%tab_Op(iOp)%symab    = 0 ! totally symmetric
      para_AllOp%tab_Op(iOp)%spectral = para_ana%Spectral_ScalOp


      DO i=1,para_PES%nb_scalar_Op  ! for scalar operators (Dip)
        iOp = iOp + 1
        CALL param_Op1TOparam_Op2(para_AllOp%tab_Op(2),                 &
                                  para_AllOp%tab_Op(iOp))
        para_AllOp%tab_Op(iOp)%n_Op    = i
        para_AllOp%tab_Op(iOp)%name_Op = 'OpScal' // int_TO_char(i)

        CALL Init_TypeOp(para_AllOp%tab_Op(iOp)%param_TypeOp,           &
                         type_Op=0,nb_Qact=mole%nb_act1,cplx=.FALSE.,   &
                         JRot=Para_Tnum%JJ,direct_KEO=.FALSE.)
        CALL derive_termQact_TO_derive_termQdyn(                        &
                              para_AllOp%tab_Op(iOp)%derive_termQdyn,   &
                              para_AllOp%tab_Op(iOp)%derive_termQact,   &
                              mole%ActiveTransfo%list_QactTOQdyn)

        para_AllOp%tab_Op(iOp)%symab    = -1  ! the symmetry is not used
        para_AllOp%tab_Op(iOp)%spectral = para_ana%Spectral_ScalOp
      END DO

      IF (para_PES%nb_scalar_Op == 3) THEN ! dipole moment operators
        para_AllOp%tab_Op(3)%name_Op = 'Dipx'
        para_AllOp%tab_Op(4)%name_Op = 'Dipy'
        para_AllOp%tab_Op(5)%name_Op = 'Dipz'
      END IF



      DO i=1,para_AllOp%nb_Op
        write(out_unitp,*) i,'Operator name: ',trim(para_AllOp%tab_Op(i)%name_Op)
      END DO

      write(out_unitp,*)
      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "============================================================"


      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "====== Finalyze RPH transfo (Tnum) and ... ================="
      write(out_unitp,*) "====== ... EneH0 of the basis sets ========================="

      write(out_unitp,*)
      IF (associated(mole%RPHTransfo)) THEN
        CALL Set_paraPRH(mole,para_Tnum,para_AllBasis%BasisnD)
      END IF

      IF (para_propa%para_Davidson%NewVec_type == 4) THEN
        CALL RecSet_EneH0(para_Tnum,mole,para_AllBasis%BasisnD,         &
                          para_PES,para_ReadOp,ComOp)
      END IF

      write(out_unitp,*)
      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "============================================================"

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
