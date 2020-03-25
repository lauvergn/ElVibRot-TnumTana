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
MODULE mod_nDGridFit
      USE mod_system
      USE mod_nDindex, only : Type_nDindex
      IMPLICIT NONE

      TYPE param_nDGrid ! it mays change in the futur (more like "basis" type)

      TYPE (Type_nDindex) :: nDindG  ! enable to use multidimensional index for the basis
      real (kind=Rkind), pointer :: Q0(:) => null()
      real (kind=Rkind), pointer :: stepQ(:) => null()

      TYPE (param_file)        :: Grid_FOR_Fit_file1
      TYPE (param_file)        :: Grid_FOR_Fit_file2
      character (len=Line_len) :: name_Grid = ''

      integer :: MinCoupling = 0
      integer :: MaxCoupling = 4
      integer :: nb_DeltaQ   = 20
      integer :: nb_G        = 0


      real (kind=Rkind) :: MinNorm = ZERO
      real (kind=Rkind) :: MaxNorm = TWO

      END TYPE param_nDGrid

      PRIVATE
      PUBLIC  :: sub_nDGrid_nDfit


      CONTAINS

      SUBROUTINE dealloc_nDGrid(para_nDGrid)
      USE mod_nDindex, only : dealloc_nDindex
      IMPLICIT NONE

      TYPE (param_nDGrid), intent(inout) :: para_nDGrid

!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'dealloc_nDGrid'
!      logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------

      CALL dealloc_nDindex(para_nDGrid%nDindG)
      IF (associated(para_nDGrid%Q0))                                   &
                CALL dealloc_array(para_nDGrid%Q0,'Q0',name_sub)
      IF (associated(para_nDGrid%stepQ))                                &
                CALL dealloc_array(para_nDGrid%stepQ,'stepQ',name_sub)

      para_nDGrid%MinNorm      = ZERO
      para_nDGrid%MaxNorm      = TWO
      para_nDGrid%MaxCoupling  = 4
      para_nDGrid%MinCoupling  = 0
      para_nDGrid%nb_DeltaQ    = 20
      para_nDGrid%nb_G         = 0

      CALL file_close(para_nDGrid%Grid_FOR_Fit_file1)
      CALL file_close(para_nDGrid%Grid_FOR_Fit_file2)

      END SUBROUTINE dealloc_nDGrid


      SUBROUTINE sub_nDGrid_nDfit()
      USE mod_system
      USE mod_dnSVM
      USE mod_Constant
      USE mod_Coord_KEO
      USE mod_PrimOp
      IMPLICIT NONE

!
!=====================================================================
!
!     variables
!
!=====================================================================
!

!----- physical and mathematical constants ---------------------------
      TYPE (constant) :: const_phys

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole
      TYPE (Tnum)    :: para_Tnum
      real (kind=Rkind), pointer :: hCC(:,:)
      TYPE(Type_dnMat) :: dnGG

!----- variables for the construction of H ----------------------------
      TYPE (param_PES)    :: para_PES


!----- variables divers ----------------------------------------------
      logical                    :: Grid,Fit,Grid1_TO_Grid2,Transfo_fit,Analysis
      TYPE (param_nDGrid)        :: para_nDGrid
      TYPE (param_nDFit)         :: para_nDFit
      real (kind=Rkind)          :: val_nDfit
      real (kind=Rkind), pointer :: Q0(:),Q(:)
      character (len=Line_len)   :: name_Grid,name_Fit
      real (kind=Rkind), allocatable :: Qact(:)
      real (kind=Rkind)          :: auTOenergy


      namelist /Grid_Fit/ Grid,Fit,Grid1_TO_Grid2,Transfo_fit,name_Grid,name_Fit,Analysis

!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'sub_nDGrid_nDfit'
!      logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
!para_mem%mem_debug = .TRUE.
!para_mem%mem_print = .FALSE.

!=====================================================================
!=====================================================================
!=====================================================================
!=====================================================================


!---------------------------------------------------------------------
!------ read or set up the physical constants ------------------------
      CALL sub_constantes(const_phys,.TRUE.)

      auTOenergy = get_Conv_au_TO_unit('E',' ',WorkingUnit=.FALSE.)
      CALL flush_perso(out_unitp)
      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "=== COORDINATES (TNUM) ====================================="
      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "============================================================"
      CALL flush_perso(out_unitp)
!=======================================================================
!=======================================================================
      CALL versionEVRT(.TRUE.)

      !-----------------------------------------------------------------
      !     - read the coordinate tansformations :
      !     -   CoordType, polysperical, bunch...
      !     ------------------------------------------------------------
      CALL Read_CoordType(mole,para_Tnum,const_phys)
      !     ------------------------------------------------------------
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !     - read coordinate values -----------------------------------
      !     ------------------------------------------------------------
      CALL read_RefGeom(mole,para_Tnum)
      !     ------------------------------------------------------------
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !     ---- TO finalize the coordinates (NM) and the KEO ----------
      !     ------------------------------------------------------------
      CALL Finalyze_TnumTana_Coord_PrimOp(para_Tnum,mole,para_PES)
      !-----------------------------------------------------------------

      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "=== END COORDINATES (TNUM) ================================="
      write(out_unitp,*) "============================================================"
      write(out_unitp,*) "============================================================"
      CALL flush_perso(out_unitp)


      CALL alloc_NParray(Qact,(/mole%nb_var/),'Qact',name_sub)
      CALL get_Qact0(Qact,mole%ActiveTransfo)

      nullify(Q0)
      CALL alloc_array(Q0,(/ mole%nb_act /),'Q0',name_sub)
      Q0(:) = Qact(1:mole%nb_act)

!=====================================================================
        name_Grid = "Grid_FOR_Fit"
        name_Fit  = "Param_FOR_Fit-col"

        Grid           = .FALSE.
        Fit            = .FALSE.
        Grid1_TO_Grid2 = .FALSE.
        Transfo_fit    = .FALSE.
        Analysis       = .FALSE.
        read(in_unitp,Grid_Fit,IOSTAT=err_read)
        IF (err_read < 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' End-of-file or End-of-record'
          write(out_unitp,*) ' The namelist "Grid_Fit" is probably absent'
          write(out_unitp,*) ' check your data!'
          write(out_unitp,*) ' ERROR in ',name_sub
          STOP
        ELSE IF (err_read > 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' Some parameter name of the namelist "Grid_Fit" are probaly wrong'
          write(out_unitp,*) ' check your data!'
          write(out_unitp,Grid_Fit)
          write(out_unitp,*) ' ERROR in ',name_sub
          STOP
        END IF
        IF (debug) write(out_unitp,Grid_Fit)

        para_nDGrid%name_Grid = name_Grid
        para_nDFit%Analysis   = Analysis
        para_nDFit%name_Fit   = name_Fit

        IF (Grid1_TO_Grid2) THEN
          CALL sub_nGrid1_TO_nGrid2(Q0)
          STOP
        END IF

        IF (Grid) THEN
          CALL sub_nDGrid(para_nDGrid,Qact,mole,para_PES)
        ELSE
          CALL sub_nDGrid_WiTHOUT_calc(para_nDGrid,Qact,mole,para_PES)
        END IF

        IF (Fit) THEN
          CALL sub_nDFit(para_nDFit,para_nDGrid,Qact,mole,para_PES)
        END IF

        IF (Transfo_fit) THEN
          CALL nDFit1_TO_TnDFit2()
        END IF

        IF (Analysis) THEN
          CALL Read_Analysis(para_nDFit%para_Analysis,Q0)
          CALL Analysis_nDFitW(para_nDFit,auTOenergy)
        END IF

        IF (para_PES%nDfit_Op) THEN
          para_PES%para_nDFit_V%Param_Fit_file%name = para_nDFit%Param_Fit_file%name

          nullify(Q)
          CALL alloc_array(Q,(/ mole%nb_act /),'Q',name_sub)

          Q(:) = Q0(:)
          Q(1) = Q0(1) + ZERO

          CALL sub_nDFunc_FROM_nDFit(val_nDfit,Q,para_PES%para_nDFit_V)
          write(out_unitp,*) 'val_nDfit',val_nDfit

          CALL dealloc_array(Q,'Q',name_sub)
        END IF

!=====================================================================
      CALL dealloc_table_at(const_phys%mendeleev)

      CALL dealloc_CoordType(mole)
      IF (associated(para_Tnum%Gref)) THEN
        CALL dealloc_array(para_Tnum%Gref,"para_Tnum%Gref",name_sub)
      END IF

      CALL dealloc_array(Q0,'Q0',name_sub)


      write(out_unitp,*) 'mem_tot',para_mem%mem_tot

      END SUBROUTINE sub_nDGrid_nDfit
      SUBROUTINE sub_nDGrid(para_nDGrid,Qact,mole,para_PES)
      USE mod_system
      USE mod_dnSVM
      USE mod_nDindex
      USE mod_Constant,  only : get_Conv_au_TO_unit
      USE mod_Coord_KEO, only : CoordType, Tnum, get_Qact0
      USE mod_PrimOp
      IMPLICIT NONE

      TYPE (param_nDGrid), intent(inout) :: para_nDGrid

!=====================================================================
!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType)      :: mole
      TYPE (Tnum)         :: para_Tnum
      real (kind=Rkind)   :: Qact(:)

!----- variables for the construction of H ---------------------------
      TYPE (param_PES)    :: para_PES
!=====================================================================


!----- variables divers ----------------------------------------------
      integer :: nderivE,MinCoupling,MaxCoupling,nb_DeltaQ
      real (kind=Rkind) :: MinNorm,MaxNorm,conv_ene
      real (kind=Rkind), pointer :: nDweight(:)

      integer, pointer :: nDsize(:)
      integer, pointer :: nDinit(:)
      integer :: i,ic,iGP,isign,nb_couplings,iGPtot

      integer           :: iOp,nioGrid1,nioGrid2

      integer                  :: idum
      character (len=Name_len) :: name_dum
      logical                  :: nDsize_read,nb_Gonly,add_Grid

      TYPE (param_dnMatOp) :: dnMatOp(1)


      TYPE (Type_nDindex),pointer :: Tab_PointGridSign(:) ! enable to change the sing of the DelatQ

      namelist /nDGrid/ MinNorm,MaxNorm,MinCoupling,MaxCoupling,        &
                        nDsize_read,nb_DeltaQ,nb_Gonly,add_Grid

!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'sub_nDGrid'
      logical, parameter :: debug=.FALSE.
!      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      CALL alloc_array(para_nDGrid%Q0,(/mole%nb_act/),                  &
                     'para_nDGrid%Q0',name_sub)
      para_nDGrid%Q0(:) = Qact(1:mole%nb_act)
      CALL alloc_array(para_nDGrid%stepQ,(/mole%nb_act/),               &
                      'para_nDGrid%stepQ',name_sub)
      para_nDGrid%stepQ(:) = ZERO

      nullify(nDsize)
      CALL alloc_array(nDsize,(/mole%nb_act/),'nDsize',name_sub)

      nullify(nDinit)
      CALL alloc_array(nDinit,(/mole%nb_act/),'nDinit',name_sub)
      nDinit(:) = 1
      nullify(nDweight)
      CALL alloc_array(nDweight,(/mole%nb_act/),'nDweight',name_sub)

      CALL Init_Tab_OF_dnMatOp(dnMatOp,mole%nb_act,para_PES%nb_elec,    &
                               nderiv=0,cplx=.FALSE.,JRot=0)

!=====================================================================
      conv_ene = get_Conv_au_TO_unit('E',' ',WorkingUnit=.FALSE.)


      para_nDGrid%Grid_FOR_Fit_file1%name = adjustl(trim(para_nDGrid%name_Grid)) // '1'
      para_nDGrid%Grid_FOR_Fit_file2%name = adjustl(trim(para_nDGrid%name_Grid)) // '2'

      CALL file_open(para_nDGrid%Grid_FOR_FIT_file1,nioGrid1)
      CALL file_open(para_nDGrid%Grid_FOR_FIT_file2,nioGrid2)

      add_Grid = .TRUE.
      DO ; IF (.NOT. add_Grid) EXIT  ! for the end of the loop

        MinNorm     = ZERO
        MaxNorm     = TWO
        MaxCoupling = 4
        MinCoupling = 0
        nb_DeltaQ   = 20
        nDsize_read = .FALSE.
        nb_Gonly    = .FALSE.
        add_Grid    = .FALSE.
        read(in_unitp,nDGrid,IOSTAT=err_read)
        IF (err_read < 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' End-of-file or End-of-record'
          write(out_unitp,*) ' The namelist "nDGrid" is probably absent'
          write(out_unitp,*) ' check your data!'
          write(out_unitp,*) ' ERROR in ',name_sub
          STOP
        ELSE IF (err_read > 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' Some parameter name of the namelist "nDGrid" are probaly wrong'
          write(out_unitp,*) ' check your data!'
          write(out_unitp,nDGrid)
          write(out_unitp,*) ' ERROR in ',name_sub
          STOP
        END IF
        IF (debug) write(out_unitp,nDGrid)
        read(in_unitp,*) nDweight(:)
        read(in_unitp,*) para_nDGrid%stepQ(:)
        IF (nDsize_read) read(in_unitp,*) nDsize(:)

        IF (debug) write(out_unitp,*) 'nDweight',nDweight(:)
        IF (debug) write(out_unitp,*) 'stepQ',para_nDGrid%stepQ(:)
        CALL flush_perso(out_unitp)


        para_nDGrid%MinNorm      = MinNorm
        para_nDGrid%MaxNorm      = MaxNorm
        para_nDGrid%MaxCoupling  = MaxCoupling
        para_nDGrid%MinCoupling  = MinCoupling
        para_nDGrid%nb_DeltaQ    = nb_DeltaQ
        IF (.NOT. nDsize_read) nDsize(:) = nb_DeltaQ
        IF (debug) write(out_unitp,*) 'nDsize',nDsize(:)
        CALL flush_perso(out_unitp)

        CALL init_nDindexPrim(para_nDGrid%nDindG,mole%nb_act,nDsize,    &
                              nDweight=nDweight,type_OF_nDindex=0,      &
                              nDinit=nDinit,                            &
                              MaxNorm=MaxNorm,MinNorm=MinNorm,          &
                              MaxCoupling=MaxCoupling,MinCoupling=MinCoupling)
        CALL sort_nDindex(para_nDGrid%nDindG)
        para_nDGrid%nDindG%Tab_nDval(:,:) = para_nDGrid%nDindG%Tab_nDval(:,:) - 1
        CALL Write_nDindex(para_nDGrid%nDindG)

        nullify(Tab_PointGridSign)
        CALL alloc_array(Tab_PointGridSign,(/MaxCoupling/),             &
                        'Tab_PointGridSign',name_sub)
        DO ic=1,MaxCoupling
          CALL init_nDindexPrim(Tab_PointGridSign(ic),ic,(/ (2,i=1,ic) /), &
                                type_OF_nDindex=1)
          Tab_PointGridSign(ic)%packed = .TRUE.
          CALL pack_nDindex(Tab_PointGridSign(ic))
          WHERE (Tab_PointGridSign(ic)%Tab_nDval(:,:) == 1)
            Tab_PointGridSign(ic)%Tab_nDval(:,:) = -1
          ELSEWHERE
            Tab_PointGridSign(ic)%Tab_nDval(:,:) = 1
          END WHERE
          !CALL Write_nDindex(Tab_PointGridSign(ic))
        END DO

        !number of points
        iGPtot = 0
        DO iGP=1,para_nDGrid%nDindG%Max_nDI
          nb_couplings=count(para_nDGrid%nDindG%Tab_nDval(:,iGP) > 0)
          IF (nb_couplings < MinCoupling) CYCLE

          IF (nb_couplings == 0) THEN
            iGPtot = iGPtot + 1
          ELSE
            iGPtot = iGPtot + Tab_PointGridSign(nb_couplings)%Max_nDI
          END IF
        END DO
        write(out_unitp,*) "nb_G",iGPtot
        write(out_unitp,*) "======================================"

        IF (nb_Gonly) CYCLE

        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
        iGPtot = 0
        CALL get_Qact0(Qact,mole%ActiveTransfo)
        DO iGP=1,para_nDGrid%nDindG%Max_nDI
          nb_couplings=count(para_nDGrid%nDindG%Tab_nDval(:,iGP) > 0)
          IF (nb_couplings < MinCoupling) CYCLE

          IF (nb_couplings == 0) THEN
            Qact(1:mole%nb_act) = para_nDGrid%Q0(:)

            CALL Set_ZERO_TO_Tab_OF_dnMatOp(dnMatOp)
            CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole,para_Tnum,para_PES,nderiv=0)

            iGPtot = iGPtot +1
            isign = 0
            nDinit(:) = para_nDGrid%nDindG%Tab_nDval(:,iGP)

            IF (debug) THEN
              write(out_unitp,"(a,3i7,a,10i3)") 'iGPtot,iGp,isign,Qact',&
                                       iGPtot,iGp,isign,' : ',nDinit(:)
              write(out_unitp,*) 'Grid',iGPtot,Qact(1:mole%nb_act),     &
                             dnMatOp(1)%tab_dnMatOp(:,:,1)%d0*conv_ene, &
                (dnMatOp(iOp)%tab_dnMatOp(:,:,1)%d0,iOp=3,size(dnMatOp))
            END IF

            write(nioGrid1,*) 'Grid1',iGPtot,nDinit(:),                 &
                                      dnMatOp(1)%tab_dnMatOp(:,:,1)%d0, &
                (dnMatOp(iOp)%tab_dnMatOp(:,:,1)%d0,iOp=3,size(dnMatOp))
            write(nioGrid2,*) 'Grid2',iGPtot,Qact(1:mole%nb_act),       &
                                      dnMatOp(1)%tab_dnMatOp(:,:,1)%d0, &
                (dnMatOp(iOp)%tab_dnMatOp(:,:,1)%d0,iOp=3,size(dnMatOp))
          ELSE

            DO isign=1,Tab_PointGridSign(nb_couplings)%Max_nDI
              iGPtot = iGPtot +1
              ic = 0
              nDinit(:) = para_nDGrid%nDindG%Tab_nDval(:,iGP)
              DO i=1,mole%nb_act
                IF (para_nDGrid%nDindG%Tab_nDval(i,iGP) > 0) THEN
                  ic = ic +1
                  nDinit(i) = para_nDGrid%nDindG%Tab_nDval(i,iGP) *     &
                     Tab_PointGridSign(nb_couplings)%Tab_nDval(ic,isign)
                END IF
              END DO
              Qact(1:mole%nb_act) = para_nDGrid%Q0(:) +                 &
                       real(nDinit(:),kind=Rkind) * para_nDGrid%stepQ(:)

              CALL Set_ZERO_TO_Tab_OF_dnMatOp(dnMatOp)
              CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole,para_Tnum,para_PES,nderiv=0)

              IF (debug) THEN
                write(out_unitp,"(a,3i7,a,10i3)") 'iGPtot,iGp,isign,Qact',&
                                       iGPtot,iGp,isign,' : ',nDinit(:)
                write(out_unitp,*) 'Grid',iGPtot,Qact(1:mole%nb_act),   &
                             dnMatOp(1)%tab_dnMatOp(:,:,1)%d0*conv_ene, &
                (dnMatOp(iOp)%tab_dnMatOp(:,:,1)%d0,iOp=3,size(dnMatOp))
              END IF

              write(nioGrid1,*) 'Grid1',iGPtot,nDinit(:),               &
                                      dnMatOp(1)%tab_dnMatOp(:,:,1)%d0, &
                (dnMatOp(iOp)%tab_dnMatOp(:,:,1)%d0,iOp=3,size(dnMatOp))
              write(nioGrid2,*) 'Grid2',iGPtot,Qact(1:mole%nb_act),     &
                                      dnMatOp(1)%tab_dnMatOp(:,:,1)%d0, &
                (dnMatOp(iOp)%tab_dnMatOp(:,:,1)%d0,iOp=3,size(dnMatOp))
            END DO
          END IF
        END DO
        para_nDGrid%nb_G = iGPtot
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"

      END DO

      CALL file_close(para_nDGrid%Grid_FOR_Fit_file1)
      CALL file_close(para_nDGrid%Grid_FOR_Fit_file2)

      IF (nb_Gonly) STOP 'nb_Gonly=t'


      CALL dealloc_array(nDsize,'nDsize',name_sub)
      CALL dealloc_array(nDinit,'nDinit',name_sub)
      DO ic=1,MaxCoupling
        CALL dealloc_nDindex(Tab_PointGridSign(ic))
      END DO
      CALL dealloc_array(Tab_PointGridSign,'Tab_PointGridSign',name_sub)

      CALL dealloc_Tab_OF_dnMatOp(dnMatOp)

      END SUBROUTINE sub_nDGrid
      SUBROUTINE sub_nDGrid_WiTHOUT_calc(para_nDGrid,Qact,mole,para_PES)
      USE mod_system
      USE mod_Constant,  only : get_Conv_au_TO_unit
      USE mod_Coord_KEO, only : CoordType, Tnum, get_Qact0
      USE mod_PrimOp
      IMPLICIT NONE

      TYPE (param_nDGrid), intent(inout) :: para_nDGrid

!=====================================================================
!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType)      :: mole
      TYPE (Tnum)         :: para_Tnum
      real (kind=Rkind)   :: Qact(:)

!----- variables for the construction of H ---------------------------
      TYPE (param_PES)    :: para_PES
!=====================================================================


!----- variables divers ----------------------------------------------
      real(kind=Rkind), pointer :: Q(:),val(:)
      integer :: iGPtot

      integer           :: nioGrid1,nioGrid2

      integer                  :: idum
      character (len=Name_len) :: name_dum
      real(kind=Rkind)         :: auTOenergy


!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'sub_nDGrid_WiTHOUT_calc'
!      logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
       auTOenergy = get_Conv_au_TO_unit('E',' ',WorkingUnit=.FALSE.)

       CALL get_Qact0(Qact,mole%ActiveTransfo)


       CALL alloc_array(para_nDGrid%Q0,(/mole%nb_act/),                 &
                       'para_nDGrid%Q0',name_sub)
        para_nDGrid%Q0(:) = Qact(1:mole%nb_act)

        nullify(Q)
        CALL alloc_array(Q,(/mole%nb_act/),'Q',name_sub)
        nullify(val)
        CALL alloc_array(val,(/4/),'val',name_sub)

        para_nDGrid%Grid_FOR_Fit_file1%name = adjustl(trim(para_nDGrid%name_Grid)) // '1'
        para_nDGrid%Grid_FOR_Fit_file2%name = adjustl(trim(para_nDGrid%name_Grid)) // '2'

        CALL file_open(para_nDGrid%Grid_FOR_FIT_file2,nioGrid2)
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
        iGPtot = 0
        DO
          !read(nioGrid2,*,IOSTAT=err_read) name_dum,idum,Q(:),val(:)
          read(nioGrid2,*,IOSTAT=err_read) name_dum,idum,Q(:)

          IF (err_read < 0) EXIT
          !write(6,*) name_dum,idum,Q(:),val(1)/auTOenergy,val(2:4)
          iGPtot = iGPtot + 1
        END DO
        para_nDGrid%nb_G = iGPtot

        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
        CALL file_close(para_nDGrid%Grid_FOR_Fit_file2)

        CALL dealloc_array(Q,'Q',name_sub)
        CALL dealloc_array(val,'val',name_sub)


      END SUBROUTINE sub_nDGrid_WiTHOUT_calc
      SUBROUTINE sub_nDFit(para_nDFit,para_nDGrid,Qact,mole,para_PES)
      USE mod_system
      USE mod_dnSVM
      USE mod_nDindex
      USE mod_Constant
      USE mod_Coord_KEO
      USE mod_PrimOp
      IMPLICIT NONE

      TYPE (param_nDGrid), intent(in)   :: para_nDGrid
      TYPE (param_nDFit), intent(inout) :: para_nDFit

!=====================================================================
!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType)      :: mole
      TYPE (Tnum)         :: para_Tnum
      real (kind=Rkind)   :: Qact(:)

!----- variables for the construction of H ---------------------------
      TYPE (param_PES)    :: para_PES
!=====================================================================


!----- variables divers ----------------------------------------------

      integer, pointer         :: nDinit(:)
      integer, pointer         :: nDvalG(:),nDvalB(:)
      real (kind=Rkind), pointer :: Q(:)
      integer                  :: i,ic,iGP,iB,jB,idum,nb_B,nb_WB
      character (len=Name_len) :: name_dum
      character (len=Name_len) :: name_file
      character (len=Name_len) :: name_int

      integer           :: nioGrid1,nioGrid2
      integer           :: nioFit

      real (kind=Rkind), pointer :: F(:)   ! size nb_B ! intermediate
      real (kind=Rkind), pointer :: A(:,:) ! size nb_B,nb_B
      real (kind=Rkind), pointer :: B(:)   ! size nb_B ! result
      real (kind=Rkind), pointer :: val(:) ! number of data which can be fit
      real (kind=Rkind), pointer :: w(:)
      real (kind=Rkind), pointer :: vv(:,:)
      real (kind=Rkind)          :: wmin,wmax,conv_ene,Weight_iGP


!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'sub_nDFit'
!      logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
!
!=====================================================================
        conv_ene = get_Conv_au_TO_unit('E',' ',WorkingUnit=.FALSE.)

        CALL Read_nDFit(para_nDFit,Qact(1:mole%nb_act))

        nullify(val)
        CALL alloc_array(val,(/para_nDFit%nb_val/),'val',name_sub)
        val(:) = ZERO

        nullify(Q)
        CALL alloc_array(Q,(/mole%nb_act/),'Q',name_sub)

        nullify(nDinit)
        CALL alloc_array(nDinit,(/mole%nb_act/),'nDinit',name_sub)
        nDinit(:) = 1

        nullify(nDvalB)
        CALL alloc_array(nDvalB,(/mole%nb_act/),'nDvalB',name_sub)
        nullify(nDvalG)
        CALL alloc_array(nDvalG,(/mole%nb_act/),'nDvalG',name_sub)


        CALL init_nDindexPrim(para_nDFit%nDindB,mole%nb_act,para_nDFit%nDsize,     &
                              nDweight=para_nDFit%nDweight,type_OF_nDindex=0,      &
                              nDinit=nDinit,                            &
                              MinNorm=para_nDFit%MinNorm,MaxNorm=para_nDFit%MaxNorm,&
                              MinCoupling=para_nDFit%MinCoupling,MaxCoupling=para_nDFit%MaxCoupling)
        CALL sort_nDindex(para_nDFit%nDindB)
        para_nDFit%nDindB%Tab_nDval(:,:) = para_nDFit%nDindB%Tab_nDval(:,:) - 1
        CALL Write_nDindex(para_nDFit%nDindB)
        nb_B = para_nDFit%nDindB%Max_nDI

        nullify(F)
        CALL alloc_array(F,(/nb_B/),'F',name_sub)
        nullify(B)
        CALL alloc_array(B,(/nb_B/),'B',name_sub)
        B(:) = ZERO
        nullify(A)
        CALL alloc_array(A,(/nb_B,nb_B/),'A',name_sub)
        A(:,:) = ZERO
        nullify(w)
        CALL alloc_array(w,(/nb_B/),'w',name_sub)
        nullify(vv)
        CALL alloc_array(vv,(/nb_B,nb_B/),'vv',name_sub)

        write(out_unitp,*) "======================================"
        write(out_unitp,*) "=== FIT =============================="
        write(out_unitp,*) "======================================"
        CALL file_open(para_nDGrid%Grid_FOR_FIT_file2,nioGrid2)

        DO iGP=1,para_nDGrid%nb_G

          read(nioGrid2,*) name_dum,idum,Q(:),val(:)

          IF (para_nDFit%Col_FOR_WeightOFFit > 0) THEN
            Weight_iGP = exp(-para_nDFit%Scal_FOR_WeightOFFit *         &
                                    val(para_nDFit%Col_FOR_WeightOFFit))
          ELSE
            Weight_iGP = ONE
          END IF

          DO iB=1,nb_B
            F(iB) = nDFunct_WITH_Q(Q,iB,para_nDFit)
          END DO

           DO iB=1,nb_B
             DO jB=1,nb_B
               A(jB,iB) = A(jB,iB) + Weight_iGP*F(iB)*F(jB)
             END DO
             B(iB)      = B(iB)    + Weight_iGP*F(iB)*val(para_nDFit%ind_val)
           END DO

        END DO

        IF (para_nDFit%svd) THEN
          !une facon.... SVD
          !write(6,*) 'a',a
          !write(6,*) 'b',b
          CALL SVDCMP(A,nb_B,nb_B,w,vv,nb_B)
          !Find maximum singular value
          wmax = maxval(w)
          !Define "small"
          wmin = wmax * para_nDFit%epsi
          !Zero the "small" singular values
          WHERE (w<wmin) w = ZERO

          F(:) = B(:)
          CALL SVBKSB(a,w,vv,nb_B,nb_B,F,b,nb_B)
        ELSE
          STOP 'not SVD, not yet'
        END IF
        CALL file_close(para_nDGrid%Grid_FOR_Fit_file2)
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"


        para_nDFit%ndim = mole%nb_act
        CALL ReadWrite_nDFitW(para_nDFit,.FALSE.,B=B,conv_ene=conv_ene)

        ! check the fit
        CALL sub_ChecknDFit2(para_nDFit,para_nDGrid)


        ! check the fit
        CALL Analysis_nDFit(para_nDFit,conv_ene=conv_ene)
        ! check the fit
        !CALL Analysis_nDFitW(para_nDFit,conv_ene=conv_ene)


        CALL dealloc_array(Q,'Q',name_sub)
        CALL dealloc_array(nDinit,'nDinit',name_sub)
        CALL dealloc_array(nDvalB,'nDvalB',name_sub)
        CALL dealloc_array(nDvalG,'nDvalG',name_sub)

        CALL dealloc_array(F,'F',name_sub)
        CALL dealloc_array(B,'B',name_sub)
        CALL dealloc_array(A,'A',name_sub)
        CALL dealloc_array(w,'w',name_sub)
        CALL dealloc_array(vv,'vv',name_sub)


      END SUBROUTINE sub_nDFit

      SUBROUTINE sub_ChecknDFit2(para_nDFit,para_nDGrid)
      USE mod_system
      USE mod_Constant
      USE mod_PrimOp
      IMPLICIT NONE

      TYPE (param_nDGrid), intent(in)   :: para_nDGrid
      TYPE (param_nDFit), intent(inout) :: para_nDFit


      integer                    :: iGP,iB
      real (kind=Rkind), pointer :: val(:),Q(:)

      integer                    :: idum
      character (len=Name_len)   :: name_dum

      integer           :: nioGrid1,nioGrid2
      real (kind=Rkind) :: val_fit,max_err,RMS_err


!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'sub_ChecknDFit2'
!      logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------

!
!=====================================================================
        nullify(val)
        CALL alloc_array(val,(/para_nDFit%ind_val/),'val',name_sub)
        val(:) = ZERO

        nullify(Q)
        CALL alloc_array(Q,(/para_nDFit%nDindB%ndim/),'Q',name_sub)


        write(out_unitp,*) "======================================"
        write(out_unitp,*) "=== CHECKING FIT ====================="
        write(out_unitp,*) "======================================"

        CALL file_open(para_nDGrid%Grid_FOR_FIT_file2,nioGrid2)

        max_err = ZERO
        RMS_err = ZERO

        DO iGP=1,para_nDGrid%nb_G

          read(nioGrid2,*) name_dum,idum,Q(:),val(:)

          CALL sub_nDFunc_FROM_nDFit(val_fit,Q,para_nDFit)


          IF (abs(val(para_nDFit%ind_val)-val_fit) > FIVE*ONETENTH**4) THEN
            IF (para_nDFit%ind_val == 1) THEN
              write(out_unitp,*) 'iG,val,val_fit,diff (cm-1)',iGP,        &
                                      val(para_nDFit%ind_val),val_fit,    &
                                    abs(val(para_nDFit%ind_val)-val_fit)* &
                                      get_Conv_au_TO_unit('E','cm-1')
            ELSE
              write(out_unitp,*) 'iG,val,val_fit,diff',iGP,               &
                                      val(para_nDFit%ind_val),val_fit,    &
                                      abs(val(para_nDFit%ind_val)-val_fit)
            END IF
          END IF


          val_fit = abs(val(para_nDFit%ind_val)-val_fit)
          IF (val_fit > max_err) max_err = val_fit
          RMS_err = RMS_err + val_fit**2
        END DO
        CALL file_close(para_nDGrid%Grid_FOR_Fit_file2)


        RMS_err = sqrt(RMS_err/para_nDGrid%nb_G)

        IF (para_nDFit%ind_val == 1) THEN
          write(out_unitp,*) 'Max_error (cm-1): ',max_err*              &
                                         get_Conv_au_TO_unit('E','cm-1')

          write(out_unitp,*) 'RMS_error (cm-1): ',RMS_err*              &
                                         get_Conv_au_TO_unit('E','cm-1')

        ELSE
          write(out_unitp,*) 'Max_error: ',max_err
          write(out_unitp,*) 'RMS_error: ',RMS_err
        END IF

        write(out_unitp,*) "======================================"
        write(out_unitp,*) "=== END CHECKING FIT ================="
        write(out_unitp,*) "======================================"


        CALL dealloc_array(val,'val',name_sub)
        CALL dealloc_array(Q,'Q',name_sub)


      END SUBROUTINE sub_ChecknDFit2

      SUBROUTINE sub_nGrid1_TO_nGrid2(Q0)
      USE mod_system
      IMPLICIT NONE

      real (kind=Rkind) :: Q0(:)

      integer, pointer           :: nDvalG1(:),nDvalG2(:)
      integer, pointer           :: list1(:)
      real (kind=Rkind), pointer :: val(:),Q1(:),Q2(:)



      integer           :: nioGrid1,nioGrid2,i
      integer                    :: idum
      character (len=Name_len)   :: name_dum

      ! for the namelist
      integer                    :: nb_act1,nb_act2
      character (len=Line_len)   :: name_Grid1,name_Grid2
      character (len=Line_len)   :: name_file

      namelist /nGrid1_TO_nGrid2/ nb_act1,name_Grid1,name_Grid2

      ! for the namelist


!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'sub_nGrid1_TO_nGrid2'
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------

!
!=====================================================================
        nb_act2 = size(Q0)

        write(out_unitp,*) 'Q0',Q0(:)

        nb_act1    = 0
        name_Grid1 = ''
        name_Grid2 = ''
        read(in_unitp,nGrid1_TO_nGrid2,IOSTAT=err_read)
        IF (err_read < 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' End-of-file or End-of-record'
          write(out_unitp,*) ' The namelist "nGrid1_TO_nGrid2" is probably absent'
          write(out_unitp,*) ' ERROR in ',name_sub
          STOP
        ELSE IF (err_read > 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' Some parameter name of the namelist "nGrid1_TO_nGrid2" are probaly wrong'
          write(out_unitp,*) ' It should never append !!'
          write(out_unitp,*) ' Check the fortran'
          write(out_unitp,*) ' ERROR in ',name_sub
          STOP
        END IF
        IF (debug) write(out_unitp,nGrid1_TO_nGrid2)

        nullify(val)
        CALL alloc_array(val,(/4/),'val',name_sub)
        nullify(Q1)
        CALL alloc_array(Q1,(/nb_act1/),'Q1',name_sub)
        nullify(Q2)
        CALL alloc_array(Q2,(/nb_act2/),'Q2',name_sub)

        nullify(nDvalG1)
        CALL alloc_array(nDvalG1,(/nb_act1/),'nDvalG1',name_sub)
        nullify(list1)
        CALL alloc_array(list1,(/nb_act1/),'list1',name_sub)
        nullify(nDvalG2)
        CALL alloc_array(nDvalG2,(/nb_act2/),'nDvalG2',name_sub)


        read(in_unitp,*) list1(:) ! for list of active variables of grid 1
        IF (debug) write(out_unitp,*) list1(:)
        CALL flush_perso(out_unitp)


        write(out_unitp,*) "======================================"
        write(out_unitp,*) "=== READING THE GRIDS ================"
        write(out_unitp,*) "======================================"
        CALL flush_perso(out_unitp)

        ! for the file with the grid points
        name_file = adjustl(trim(name_Grid1)) // '2'
        CALL file_open2(name_file,nioGrid1)
        name_file = adjustl(trim(name_Grid2)) // '2'
        CALL file_open2(name_file,nioGrid2)

        DO

          read(nioGrid1,*,IOSTAT=err_read) name_dum,idum,Q1(:),val(:)
          IF (err_read /= 0) EXIT

          Q2(:) = Q0(:)
          DO i=1,nb_act1
            Q2(list1(i)) = Q1(i)
          END DO
          write(nioGrid2,*) name_dum,idum,Q2(:),val(:)


        END DO
        close(nioGrid1)
        close(nioGrid2)

        ! for the file with the basis indices
        name_file = adjustl(trim(name_Grid1)) // '1'
        CALL file_open2(name_file,nioGrid1)
        name_file = adjustl(trim(name_Grid2)) // '1'
        CALL file_open2(name_file,nioGrid2)

        DO

          read(nioGrid1,*,IOSTAT=err_read) name_dum,idum,nDvalG1(:),val(:)
          IF (err_read /= 0) EXIT

          nDvalG2(:) = 1
          DO i=1,nb_act1
            nDvalG2(list1(i)) = nDvalG1(i)
          END DO
          write(nioGrid2,*) name_dum,idum,nDvalG2(:),val(:)


        END DO
        close(nioGrid1)
        close(nioGrid2)


        CALL dealloc_array(val,    'val',name_sub)
        CALL dealloc_array(Q1,     'Q1',name_sub)
        CALL dealloc_array(Q2,     'Q2',name_sub)
        CALL dealloc_array(nDvalG1,'nDvalG1',name_sub)
        CALL dealloc_array(list1,  'list1',name_sub)
        CALL dealloc_array(nDvalG2,'nDvalG2',name_sub)

      END SUBROUTINE sub_nGrid1_TO_nGrid2

END MODULE mod_nDGridFit
