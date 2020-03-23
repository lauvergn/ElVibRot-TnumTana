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
!=============================================================
!
!      Build the direct product basis set
!
!=============================================================
      SUBROUTINE sub_quadra_SparseBasis2n(SparseBasis,mole)
      USE mod_system
      USE mod_nDindex
      USE mod_Coord_KEO
      USE mod_basis
      IMPLICIT NONE

!---------------------------------------------------------------------
      TYPE (basis)  :: SparseBasis

!----- variables for the inactive namelist ----------------
      TYPE (CoordType)     :: mole

!---------------------------------------------------------------------
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      integer       :: iq_modif


!----  For Sparse Grid -----------------------------------------------
      integer :: ndim,L_SparseGrid,nb_grid,nb_basis,nb_points,nb_remove
      TYPE (basis), pointer  :: tab_basis_loc(:,:)

      integer :: ib,L,iq,iq_loc,tab_iq(1),iqf,iqi,nq,nqL
      integer, allocatable :: tab_iq_loc_TO_iq(:,:,:)
      real (kind=Rkind), allocatable :: x(:)
      real (kind=Rkind), allocatable :: w(:)
      real (kind=Rkind), allocatable :: w2(:)

      real (kind=Rkind) :: x_loc
      integer, allocatable :: list_nDgrid_points(:,:)
      integer, allocatable :: list2_nDgrid_points(:,:)
      integer :: nb2_points,iqnD,iq2nD,ident,iqnD_id

      integer :: ig,minL,maxL
      integer, allocatable :: max_nbL_basis(:)
      integer, allocatable :: tab_nb_Grid_per_L(:)

      integer, allocatable :: indL_per_nDgrid(:,:)     ! indgrid_per_nDgrid(nb_basis,nb_grid)
      integer, allocatable :: indgrid_per_nDgrid(:,:)  ! indgrid_per_nDgrid(nb_basis,nb_grid)

      integer           :: normL,DeltaL
      real (kind=Rkind) :: C
      integer, allocatable  :: min_nq(:)
      integer, allocatable  :: max_nq(:)
      integer, allocatable  :: ind_nq(:)

      real (kind=Rkind) :: rnb_points
      integer :: err_sub


!----- function ------------------------------------------------------
      integer           :: calc_ind_n
      integer           :: L_TO_nq
      real (kind=Rkind) :: binomial
      logical           :: inferior
!---------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_quadra_SparseBasis2n'
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------
      nullify(tab_basis_loc)


      nb_basis     = SparseBasis%nb_basis
      ndim         = SparseBasis%nb_basis
      L_SparseGrid = SparseBasis%L_SparseGrid

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_basis,L_SparseGrid',nb_basis,L_SparseGrid
      END IF

!-----------------------------------------------------------
!     allocation of SparseBasis

      write(out_unitp,*) '================================================='
      write(out_unitp,*) '== SPARSE BASIS ================================='
      write(out_unitp,*) '================================================='
      SparseBasis%packed            = .FALSE.
      SparseBasis%packed_done       = .FALSE.


      SparseBasis%nb                    = SparseBasis%nDindB%max_nDI

      SparseBasis%MaxCoupling_OF_nDindB = nb_basis
      SparseBasis%ndim                  = nb_basis

      CALL alloc_tab_Pbasis_OF_basis(SparseBasis)
      CALL alloc_init_basis(SparseBasis)

!        DO ib=1,SparseBasis%nb_basis
!          write(out_unitp,*) '--------------------------'
!          write(out_unitp,*) 'tab_basis',ib
!          CALL RecWrite_basis(SparseBasis%tab_Pbasis(ib)%Pbasis)
!        END DO
!        write(out_unitp,*) '--------------------------'
!-----------------------------------------------------------

!----------------------------------------------------------------------------
      write(out_unitp,*) '-------------------------------------------------'
      write(out_unitp,*) '-------------------------------------------------'
      write(out_unitp,*) '-------------------------------------------------'
      write(out_unitp,*) '    Basis: Sparse Basis'
      write(out_unitp,*) 'nb_basis,L_SparseGrid',nb_basis,L_SparseGrid
      IF (SparseBasis%cplx) THEN
        write(out_unitp,*) ' STOP the basis is complex!!'
        STOP
      END IF

!----------------------------------------------------------------------------
!     1st: calculation of grid points and weight for the 1D-grid
!         => in tab_basis_loc

      CALL alloc_array(tab_basis_loc,(/ L_SparseGrid,nb_basis /),       &
                      'tab_basis_loc',name_sub,(/0,1/))
      DO ib=1,nb_basis
      DO L=0,L_SparseGrid

        tab_basis_loc(L,ib)%nb            = 0
        tab_basis_loc(L,ib)%ndim          = 1
        tab_basis_loc(L,ib)%packed        = .TRUE.
        tab_basis_loc(L,ib)%packed_done   = .TRUE.

        CALL Set_Basis_L_TO_n(tab_basis_loc(L,ib)%L_TO_nq,A=1,B=1,C=1,expo=1,L_TO_n_type=0)

        tab_basis_loc(L,ib)%L_SparseBasis = SparseBasis%nDindB%Lmax

        nq = Get_nq_FROM_l_OF_PrimBasis(L,tab_basis_loc(L,ib))

        CALL Set_nq_OF_basis(tab_basis_loc(L,ib),nq)
        CALL alloc_init_basis(tab_basis_loc(L,ib))
!        write(out_unitp,*) 'packed',tab_basis_loc(L,ib)%packed
!        write(out_unitp,*) 'asso x,w',associated(tab_basis_loc(L,ib)%x),associated(tab_basis_loc(L,ib)%w)

        CALL alloc_xw_OF_basis(tab_basis_loc(L,ib))
        CALL hercom(nq,tab_basis_loc(L,ib)%x,tab_basis_loc(L,ib)%w)
        tab_basis_loc(L,ib)%wrho(:) = tab_basis_loc(L,ib)%w(:)  ! useless !
        tab_basis_loc(L,ib)%rho(:)  = ONE                       ! useless !

!       write(out_unitp,*) '--------------------------'
!       write(out_unitp,*) 'basis',ib,L
!       CALL RecWrite_basis(tab_basis_loc(L,iib))
!       CALL flush_perso(out_unitp)
      END DO
      END DO
      write(out_unitp,*) 'L, 1D max_nq',L_SparseGrid,                   &
       (get_nq_FROM_basis(tab_basis_loc(L_SparseGrid,ib)),ib=1,nb_basis)
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
!     2d: 1D-grid in SparseBasis%tab_basis(:) (from tab_basis_loc)
!       2a: first local grid x = union_of(tab_basis_loc(:,ib)%x
!       2b: number of the duplicate values
!       2c: x in SparseBasis%tab_basis
!       2d: d0b, d1b, d2b in SparseBasis%tab_Pbasis%Pbasis
!
      DO ib=1,nb_basis
!       2a first local grid x = union_of(tab_basis_loc(:,ib)%x
        nq = 0
        DO L=0,L_SparseGrid
          nq = nq + get_nq_FROM_basis(tab_basis_loc(L,ib))
        END DO
        CALL alloc_NParray(x,(/ nq /),"x",name_sub)
        iqi = 0
        DO L=0,L_SparseGrid
          nqL  = get_nq_FROM_basis(tab_basis_loc(L,ib))
          iqf = iqi + nqL
          x(iqi+1:iqf) = tab_basis_loc(L,ib)%x(1,:)
          iqi = iqf
        END DO

!       2b: number of the duplicate values
        CALL trie_tab(nq,x,nq) ! sort the 1D-grid:x
        nb_remove = 0
        DO iq=1,nq-1
          IF (x(iq) == x(iq+1)) nb_remove = nb_remove + 1
        END DO

!       2c: x in SparseBasis%tab_Pbasis
        SparseBasis%tab_Pbasis(ib)%Pbasis%packed            = .TRUE.
        SparseBasis%tab_Pbasis(ib)%Pbasis%packed_done       = .TRUE.


        SparseBasis%tab_Pbasis(ib)%Pbasis%ndim     = 1
        SparseBasis%tab_Pbasis(ib)%Pbasis%nb       = SparseBasis%nDindB%Lmax + 1
        SparseBasis%tab_Pbasis(ib)%Pbasis%name     = "Hm"
        SparseBasis%tab_Pbasis(ib)%Pbasis%type     = 20
        CALL Set_nq_OF_basis(SparseBasis%tab_Pbasis(ib)%Pbasis,nq-nb_remove)

        SparseBasis%tab_Pbasis(ib)%Pbasis%nbc      =                    &
                                    SparseBasis%tab_Pbasis(ib)%Pbasis%nb
        SparseBasis%tab_Pbasis(ib)%Pbasis%nqc      =                    &
                    get_nq_FROM_basis(SparseBasis%tab_Pbasis(ib)%Pbasis)
        SparseBasis%tab_Pbasis(ib)%Pbasis%contrac  = .FALSE.

        CALL alloc_init_basis(SparseBasis%tab_Pbasis(ib)%Pbasis)
        CALL alloc_xw_OF_basis(SparseBasis%tab_Pbasis(ib)%Pbasis)
        CALL alloc_dnb_OF_basis(SparseBasis%tab_Pbasis(ib)%Pbasis)

        iQ = mole%nb_act1 + ib
        SparseBasis%tab_Pbasis(ib)%Pbasis%iQdyn(1) = mole%ActiveTransfo%list_QactTOQdyn(iQ)
        CALL alloc_NParray(                                              &
              SparseBasis%tab_Pbasis(ib)%Pbasis%Tabder_Qdyn_TO_Qbasis, &
                                                    (/ mole%nb_var /), &
             'SparseBasis%tab_Pbasis(ib)%Pbasis%Tabder_Qdyn_TO_Qbasis',&
                                                      name_sub,(/ 0 /))
        SparseBasis%tab_Pbasis(ib)%Pbasis%Tabder_Qdyn_TO_Qbasis(:)  = 0
        SparseBasis%tab_Pbasis(ib)%Pbasis%Tabder_Qdyn_TO_Qbasis(iQ) = 1

        SparseBasis%iQdyn(ib) = iQ

        SparseBasis%tab_Pbasis(ib)%Pbasis%x(1,1) = x(1)
        iqi = 1
        DO iq=2,nq
          IF (x(iq) == x(iq-1)) CYCLE
          iqi = iqi + 1
          SparseBasis%tab_Pbasis(ib)%Pbasis%x(1,iqi) = x(iq)
        END DO
        SparseBasis%tab_Pbasis(ib)%Pbasis%w(:)    = ZERO
        SparseBasis%tab_Pbasis(ib)%Pbasis%rho(:)  = ONE
        SparseBasis%tab_Pbasis(ib)%Pbasis%wrho(:) = ZERO

!       2d: d0b, d1b, d2b in SparseBasis%tab_Pbasis
        nq = get_nq_FROM_basis(SparseBasis%tab_Pbasis(ib)%Pbasis)
        CALL d0d1d2poly_Hermite_exp_grille(                             &
                        SparseBasis%tab_Pbasis(ib)%Pbasis%x(1,:),       &
                        SparseBasis%tab_Pbasis(ib)%Pbasis%dnRGB%d0(:,:),     &
                        SparseBasis%tab_Pbasis(ib)%Pbasis%dnRGB%d1(:,:,1),   &
                        SparseBasis%tab_Pbasis(ib)%Pbasis%dnRGB%d2(:,:,1,1), &
                        SparseBasis%tab_Pbasis(ib)%Pbasis%nb,nq,        &
                        .TRUE.,.FALSE.,ONE)
!       write(out_unitp,*) '--------------------------'
!       write(out_unitp,*) 'tab_Pbasis',ib
!       CALL RecWrite_basis(SparseBasis%tab_Pbasis(ib)%Pbasis)
!       CALL flush_perso(out_unitp)


        CALL dealloc_NParray(x,"x",name_sub)

      END DO
      write(out_unitp,*) '1D-basis: done'
      CALL flush_perso(out_unitp)
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
!     3d: relation between iq_loc and iq: tab_basis_loc(L,ib)%x(1,iq_loc) and SparseBasis%tab_Pbasis(ib)%x(1,iq)
!         ==> Tab_iq_loc_TO_iq(0:Lmax,1:nb_basis,1:max(nq))

      nq = 0
      DO ib=1,nb_basis
      DO L=0,L_SparseGrid
        nq = max(nq,get_nq_FROM_basis(tab_basis_loc(L,ib)))
      END DO
      END DO

      CALL alloc_NParray(tab_iq_loc_TO_iq,(/ L_SparseGrid,nb_basis,nq /), &
                        "tab_iq_loc_TO_iq",name_sub,(/ 0,1,1 /))
      tab_iq_loc_TO_iq(:,:,:) = 0
      DO ib=1,nb_basis
!     write(out_unitp,*) 'ib,x ',ib,L,SparseBasis%tab_Pbasis(ib)%Pbasis%x(1,:)
      DO L=0,L_SparseGrid
!       write(out_unitp,*) 'ib,L,x_loc',ib,L,tab_basis_loc(L,ib)%x(1,:)
        DO iq_loc=1,get_nq_FROM_basis(tab_basis_loc(L,ib))
          x_loc = tab_basis_loc(L,ib)%x(1,iq_loc)
          tab_iq = maxloc(SparseBasis%tab_Pbasis(ib)%Pbasis%x(1,:),             &
                          SparseBasis%tab_Pbasis(ib)%Pbasis%x(1,:) == x_loc)

          tab_iq_loc_TO_iq(L,ib,iq_loc) = tab_iq(1)
        END DO
!       write(out_unitp,*) 'ib,L,tab_loc',ib,L,
!    *            tab_iq_loc_TO_iq(L,ib,1:get_nq_FROM_basis(tab_basis_loc(L,ib)))
      END DO
      END DO
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
!     4th: the number of direct_product grids: nb_grid
!     then the index in L for each grid: indL_per_nDgrid
      CALL alloc_NParray(max_nbL_basis,(/nb_basis/),                      &
                        "max_nbL_basis",name_sub)
      CALL alloc_NParray(tab_nb_Grid_per_L,(/L_SparseGrid/),              &
                        "tab_nb_Grid_per_L",name_sub,(/0/))
!     max_nbL_basis(:) = SparseBasis%tab_Pbasis(:)%Pbasis%nb -1
!     where(max_nbL_basis(:) > L_SparseGrid)
!    *         max_nbL_basis(:) = L_SparseGrid
      max_nbL_basis(:) = L_SparseGrid+1


!     write(out_unitp,*) 'max_nbL_basis',max_nbL_basis

      minL = max(0,L_SparseGrid - nb_basis + 1)
      maxL = L_SparseGrid

      CALL alloc_NParray(max_nq,(/nb_basis/),"max_nq",name_sub)
      CALL alloc_NParray(min_nq,(/nb_basis/),"min_nq",name_sub)
      CALL alloc_NParray(ind_nq,(/nb_basis/),"ind_nq",name_sub)
      max_nq(:) = maxL
      min_nq(:) = 0
      ind_nq(:) = -1
      tab_nb_Grid_per_L(:) = 0
      DO L=minL,maxL
        nb_grid = 0
        CALL Rec_calc_nb_HierarchicalRep(nb_grid,L,                     &
                           nb_basis,ind_nq,min_nq,max_nq,nb_basis)
        tab_nb_Grid_per_L(L) = nb_grid
      END DO
      nb_grid = sum(tab_nb_Grid_per_L)
      write(out_unitp,*) 'tab_nb_Grid_per_L',tab_nb_Grid_per_L
      write(out_unitp,*) 'nb_grid',nb_grid
      max_nq(:) = maxL
      min_nq(:) = 0
      ind_nq(:) = -1
      CALL alloc_NParray(indL_per_nDgrid,(/nb_basis,nb_grid/),            &
                      "indL_per_nDgrid",name_sub)
      ig = 0
      DO L=minL,maxL
        CALL Rec_calc_ind_HierarchicalRep(ig,L,                         &
                           nb_basis,ind_nq,                             &
                           indL_per_nDgrid,nb_grid,                     &
                           min_nq,max_nq,nb_basis)
      END DO
      CALL flush_perso(out_unitp)

      CALL dealloc_NParray(max_nq,"max_nq",name_sub)
      CALL dealloc_NParray(min_nq,"min_nq",name_sub)
      CALL dealloc_NParray(ind_nq,"ind_nq",name_sub)

!----------------------------------------------------------------------------



!----------------------------------------------------------------------------
!     5th: number of grid points: nb_points (with duplicate values)
!          and the index in q for each grid: indgrid_per_nDgrid
      CALL alloc_NParray(indgrid_per_nDgrid,(/nb_basis,nb_grid/),         &
                      "indgrid_per_nDgrid",name_sub)
      indgrid_per_nDgrid(:,:)= 0

      nb_points = 0
      DO ig=1,nb_grid
        !write(out_unitp,*) 'ig,indL_per_nDgrid',ig,indL_per_nDgrid(:,ig)
        normL = sum(indL_per_nDgrid(:,ig))
        DeltaL = SparseBasis%L_SparseGrid - normL
        IF (DeltaL < 0) STOP 'DeltaL < 0'
        IF (DeltaL > nb_basis -1) STOP 'DeltaL > nb_basis-1'

        DO ib=1,nb_basis
          L = indL_per_nDgrid(ib,ig)
          indgrid_per_nDgrid(ib,ig) = get_nq_FROM_basis(tab_basis_loc(L,ib))
        END DO


        nb_points = nb_points + product(indgrid_per_nDgrid(:,ig))
!       write(out_unitp,*) 'ig,indL   ',ig,indL_per_nDgrid(:,ig)
!       write(out_unitp,*) 'ig,indgrid',ig,indgrid_per_nDgrid(:,ig)
      END DO
      write(out_unitp,*) 'nb_points (with duplicates): ',nb_points
      CALL flush_perso(out_unitp)
!----------------------------------------------------------------------------

!----------------------------------------------------------------------------
!     6th: multidimentional index points and weight (with duplicate values)
!         => temporary tables
      CALL alloc_NParray(list_nDgrid_points,(/ nb_basis,nb_points /),     &
                        "list_nDgrid_points",name_sub)
      CALL alloc_NParray(w,(/ nb_points /),"w",name_sub)
      nb2_points = nb_points



!     write(out_unitp,*) '------------------------------------------------'
!     write(out_unitp,*) '--   normL E [L_SparseGrid-D+1 ... L_SparseGrid]'
!     write(out_unitp,*) ' L_SparseGrid-D+1',L_SparseGrid-nb_basis+1
!     write(out_unitp,*) ' L_SparseGrid',L_SparseGrid
!     write(out_unitp,*) '------------------------------------------------'


      CALL alloc_NParray(max_nq,(/ nb_basis /),"max_nq",name_sub)
      CALL alloc_NParray(ind_nq,(/ nb_basis /),"ind_nq",name_sub)
      nb_points = 0
      DO ig=1,nb_grid
        max_nq(:) = indgrid_per_nDgrid(:,ig)
        ind_nq(:) = 0

        normL = sum(indL_per_nDgrid(:,ig))
        DeltaL = L_SparseGrid - normL
!       write(out_unitp,*) 'For grid ',ig,' normL',normL,'L',
!    *                                      indL_per_nDgrid(:,ig)
!       write(out_unitp,*) 'For grid ',ig,' nb_points',product(max_nq(:))

        IF (DeltaL < 0) STOP 'DeltaL < 0'
        IF (DeltaL > nb_basis -1) STOP 'DeltaL > nb_basis-1'
        C = (-ONE)**DeltaL * binomial(nb_basis-1,deltaL)
!       write(out_unitp,*) 'For grid ',ig,' C:',C
!       CALL flush_perso(out_unitp)


        DO iq=1,product(max_nq(:))

          nb_points = nb_points + 1
!         -----------------------------------------------------
!         determine l'indice modifie (i_modif) de la variable
!         et les indices de la frequence ind_freq en fonction de i_point
          iq_modif = calc_ind_n(ind_nq,iq,max_nq,nb_basis,1)

          w(nb_points) = C
          DO ib=1,nb_basis
            L = indL_per_nDgrid(ib,ig)
!           write(out_unitp,*) 'L',L,tab_basis_loc(L,ib)%x(1,ind_nq(ib))
            w(nb_points) = w(nb_points) *                               &
                                      tab_basis_loc(L,ib)%w(ind_nq(ib))

            list_nDgrid_points(ib,nb_points) =                          &
                      tab_iq_loc_TO_iq(L,ib,ind_nq(ib))
          END DO
!         write(out_unitp,*) 'ig,iq,x:',ig,iq,
!    *        (tab_basis_loc(indL_per_nDgrid(ib,ig),ib)%x(1,ind_nq(ib)),
!    *         ib=1,nb_basis)
!         write(out_unitp,*) 'ig,iq,ind_nq,w:',ig,iq,
!    *                  list_nDgrid_points(:,nb_points),w(nb_points)
!         CALL flush_perso(out_unitp)


        END DO
      END DO
      CALL dealloc_NParray(max_nq,"max_nq",name_sub)
      CALL dealloc_NParray(ind_nq,"ind_nq",name_sub)

      IF (nb_points /= nb2_points) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' You MUST have nb_points = nb2_points',     &
                    nb_points,nb2_points
        write(out_unitp,*) ' Check the source !!!!'
        STOP
      END IF
      write(out_unitp,*) 'list_nDgrid_points: done'
      CALL flush_perso(out_unitp)
!     -----------------------------------------------------------------------
!     6b: sort the list_nDgrid_points
      CALL alloc_NParray(list2_nDgrid_points,(/ nb_basis,nb_points /),    &
                        "list2_nDgrid_points",name_sub)
      CALL alloc_NParray(w2,(/ nb_points /),"w2",name_sub)
      list2_nDgrid_points(:,:) = 0
      w2(:) = ZERO

      DO ib=1,nb_basis

        nq = get_nq_FROM_basis(SparseBasis%tab_Pbasis(ib)%Pbasis)
        CALL alloc_NParray(max_nq,(/ nq /),"max_nq",name_sub)
        CALL alloc_NParray(ind_nq,(/ nq /),"ind_nq",name_sub)
        ind_nq(:) = 0
        iq = 1
        ind_nq(iq) = 1
        max_nq(iq) = count(list_nDgrid_points(ib,:) == iq)
        DO iq=2,nq
          ind_nq(iq) = sum(max_nq(1:iq-1)) + 1
          max_nq(iq) = count(list_nDgrid_points(ib,:) == iq)
        END DO
!       write(out_unitp,*) 'ib,max_nq',ib,max_nq
!       write(out_unitp,*) 'ib,ind_nq',ib,ind_nq

        DO iqnD=1,nb_points
          iq = list_nDgrid_points(ib,iqnD)
          iq2nD = ind_nq(iq)
          ind_nq(iq) = ind_nq(iq) + 1
          list2_nDgrid_points(:,iq2nD) = list_nDgrid_points(:,iqnD)
          w2(iq2nD) = w(iqnD)
        END DO
        w(:) = w2(:)
        list_nDgrid_points(:,:) = list2_nDgrid_points(:,:)

        CALL dealloc_NParray(max_nq,"max_nq",name_sub)
        CALL dealloc_NParray(ind_nq,"ind_nq",name_sub)
      END DO
      CALL dealloc_NParray(list2_nDgrid_points,"list2_nDgrid_points",name_sub)
      CALL dealloc_NParray(w2,"w2", name_sub)

      write(out_unitp,*) 'sort list_nDgrid_points: done'
      CALL flush_perso(out_unitp)

      nb_remove = 0
      iqnD_id = 0
      iqnD = 1
!     write(out_unitp,*) 'iqnD',iqnD,':',list_nDgrid_points(:,iqnD)
      DO iqnD=2,nb_points
!       write(out_unitp,*) 'iqnD',iqnD,':',list_nDgrid_points(:,iqnD)
        ident = sum(abs(list_nDgrid_points(:,iqnD)-                     &
                        list_nDgrid_points(:,iqnD-1)))

        IF (ident == 0) THEN
          IF (iqnD_id == 0) iqnD_id = iqnD-1
          nb_remove = nb_remove + 1
          w(iqnD_id) = w(iqnD_id) + w(iqnD)
          w(iqnD) = ZERO
        ELSE
          iqnD_id = 0
        END IF
      END DO
      write(out_unitp,*) 'nb_remove',nb_remove
      CALL flush_perso(out_unitp)
!     -----------------------------------------------------------------------
!     6c: list_nDgrid_points, w in SparseBasis
      nq   = nb_points-nb_remove
      CALL Set_nq_OF_basis(SparseBasis,nq)
      SparseBasis%ndim = nb_basis
      CALL alloc_xw_OF_basis(SparseBasis)

      CALL alloc_NParray(list2_nDgrid_points,(/SparseBasis%nb_basis,nq/), &
                      "list2_nDgrid_points",name_sub)

      write(out_unitp,*) 'SparseBasis: nq',nq
      write(out_unitp,*) 'Efficiency',                                  &
           real(L_SparseGrid+1,kind=Rkind)**nb_basis/real(nq,kind=Rkind)
      CALL flush_perso(out_unitp)

      iq2nD = 0
      DO iqnD=1,nb_points
        IF (w(iqnD) /= ZERO) THEN
          iq2nD = iq2nD + 1
          SparseBasis%w(iq2nD) = w(iqnD)
          list2_nDgrid_points(:,iq2nD) = list_nDgrid_points(:,iqnD)
          DO ib=1,nb_basis
            iq = list_nDgrid_points(ib,iqnD)
            SparseBasis%x(ib,iq2nD) = SparseBasis%tab_Pbasis(ib)%Pbasis%x(1,iq)
          END DO
!         write(out_unitp,*) 'iqnD',iq2nD,'w,x:',SparseBasis%x(:,iq2nD),
!    *                                   SparseBasis%w(iq2nD)
        END IF
      END DO

      SparseBasis%wrho(:) = SparseBasis%w(:)
      SparseBasis%rho(:)  = ONE

      CALL dealloc_NParray(list_nDgrid_points,"list_nDgrid_points",name_sub)
      CALL dealloc_NParray(w,"w",name_sub)
!----------------------------------------------------------------------------
      CALL dealloc_nDindex(SparseBasis%nDindG)
      CALL init_nDindex_typeTAB(SparseBasis%nDindG,                     &
                                 nb_basis,list2_nDgrid_points,nq,err_sub)
      IF (err_sub /= 0) THEN
        write(*,*) ' ERROR in ',name_sub
        STOP ' from init_nDindex_typeTAB'
      END IF


      SparseBasis%nDindB%Tab_nDval(:,:) = SparseBasis%nDindB%Tab_nDval(:,:) + 1

      CALL dealloc_NParray(list2_nDgrid_points,"list2_nDgrid_points",name_sub)
!----------------------------------------------------------------------------
!     free memories
!----------------------------------------------------------------------------
      DO ib=1,nb_basis
      DO L=0,L_SparseGrid
        CALL dealloc_basis(tab_basis_loc(L,ib))
      END DO
      END DO

      CALL dealloc_array(tab_basis_loc,     'tab_basis_loc',     name_sub)

      CALL dealloc_NParray(tab_iq_loc_TO_iq,  'tab_iq_loc_TO_iq',  name_sub)
      CALL dealloc_NParray(indgrid_per_nDgrid,'indgrid_per_nDgrid',name_sub)
      CALL dealloc_NParray(indL_per_nDgrid,   'indL_per_nDgrid',   name_sub)
      CALL dealloc_NParray(max_nbL_basis,     'max_nbL_basis',     name_sub)
      CALL dealloc_NParray(tab_nb_Grid_per_L, 'tab_nb_Grid_per_L', name_sub)
      write(out_unitp,*) 'free memory: done'
      CALL flush_perso(out_unitp)
!----------------------------------------------------------------------

!-----------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(SparseBasis)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------
      write(out_unitp,*) ' END Basis: Sparse Basis'
      write(out_unitp,*) '-------------------------------------------------'
      write(out_unitp,*) '-------------------------------------------------'
      write(out_unitp,*) '-------------------------------------------------'

      write(out_unitp,*) '================================================='
      write(out_unitp,*) '== END SPARSE BASIS ============================='
      write(out_unitp,*) '================================================='


      END SUBROUTINE sub_quadra_SparseBasis2n

      RECURSIVE SUBROUTINE Rec_calc_nb_HierarchicalRep(i_list,sum_val,  &
                                i_val,ind_val,min_val,max_val,dim_val)
      USE mod_system
      IMPLICIT NONE


      integer, intent(inout)     :: i_list
      integer, intent(in)      :: sum_val,i_val,dim_val
      integer, intent(in)      :: min_val(dim_val)
      integer, intent(in)      :: max_val(dim_val)
      integer, intent(inout)   :: ind_val(dim_val)


      integer :: min_i_val,max_i_val,i

      !--------------------------------------------------------------
      character (len=*), parameter ::                                   &
                                name_sub='Rec_calc_nb_HierarchicalRep'
      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'i_list',i_list
        write(out_unitp,*) 'sum_val',sum_val
        write(out_unitp,*) 'i_val',i_val
        write(out_unitp,*) 'min_val',min_val(:)
        write(out_unitp,*) 'max_val',max_val(:)
        write(out_unitp,*) 'ind_val',ind_val(:)
      END IF
      !--------------------------------------------------------------

      IF (sum(min_val(:)) > sum_val) THEN
        i_list = i_list + 0
        RETURN
      END IF
      IF (i_val == 1) THEN
        i_list = i_list + 1
        ind_val(1) = sum_val
        !write(out_unitp,*) name_sub,',',i_list,':',ind_val(:)
      ELSE
        min_i_val = max(0,min_val(i_val))
        max_i_val = min(sum_val,max_val(i_val))
        !write(out_unitp,*) 'i_val..',i_val,min_i_val,max_i_val
        DO i=min_i_val,max_i_val
           ind_val(i_val) = i
           CALL Rec_calc_nb_HierarchicalRep(i_list,sum_val-i,           &
                                i_val-1,ind_val,min_val,max_val,dim_val)
        END DO
      END IF

      IF (debug) THEN
        write(out_unitp,*) name_sub,',',i_list,':',ind_val(:)
        write(out_unitp,*) 'END ',name_sub
      END IF

      END SUBROUTINE Rec_calc_nb_HierarchicalRep
      RECURSIVE SUBROUTINE Rec_calc_ind_HierarchicalRep(i_list,sum_val, &
                                i_val,ind_val,                          &
                                list_ind_val,nb_list,                   &
                                min_val,max_val,dim_val)
      USE mod_system
      IMPLICIT NONE


      integer, intent(inout)   :: i_list
      integer, intent(in)      :: sum_val,i_val,dim_val,nb_list
      integer, intent(in)      :: min_val(dim_val)
      integer, intent(in)      :: max_val(dim_val)
      integer, intent(inout)   :: ind_val(dim_val)
      integer, intent(inout)   :: list_ind_val(dim_val,nb_list)


      integer :: min_i_val,max_i_val,i

      !--------------------------------------------------------------
      character (len=*), parameter ::                                   &
                                name_sub='Rec_calc_ind_HierarchicalRep'
      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'i_list',i_list
        write(out_unitp,*) 'sum_val',sum_val
        write(out_unitp,*) 'i_val',i_val
        write(out_unitp,*) 'min_val',min_val(:)
        write(out_unitp,*) 'max_val',max_val(:)
        write(out_unitp,*) 'ind_val',ind_val(:)
      END IF
      !--------------------------------------------------------------

      IF (sum(min_val(:)) > sum_val) THEN
        i_list = i_list + 0
        RETURN
      END IF
      IF (i_val == 1) THEN
        i_list = i_list + 1
        ind_val(1) = sum_val
        list_ind_val(:,i_list) = ind_val(:)
        !write(out_unitp,*) name_sub,',',i_list,':',ind_val(:)
      ELSE
        min_i_val = max(0,min_val(i_val))
        max_i_val = min(sum_val,max_val(i_val))
        !write(out_unitp,*) 'i_val..',i_val,min_i_val,max_i_val
        DO i=min_i_val,max_i_val
           ind_val(i_val) = i
           CALL Rec_calc_ind_HierarchicalRep(i_list,sum_val-i,          &
                                i_val-1,ind_val,list_ind_val,nb_list,   &
                                min_val,max_val,dim_val)
        END DO
      END IF

      IF (debug) THEN
        write(out_unitp,*) name_sub,',',i_list,':',ind_val(:)
        write(out_unitp,*) 'END ',name_sub
      END IF

      END SUBROUTINE Rec_calc_ind_HierarchicalRep
      integer FUNCTION calc_ind_n(ind_n,i_point,n,ndim,ind0)
      USE mod_system
      IMPLICIT NONE

       integer ndim
       integer ind_n(ndim)
       integer n(ndim)
       integer i_modif,ind0,i_point,i
!      ind0 : parametre de debut de tableau ( 0 ou 1 )
!       =>  ind0 = 0  => ind_n(0,0,0,...0)
!       =>  ind0 = 1  => ind_n(1,1,1,...1)


!----- for debuging --------------------------------------------------
      logical   debug
      parameter (debug=.FALSE.)
!     parameter (debug=.TRUE.)
!---------------------------------------------------------------------
      IF (debug) THEN
      write(out_unitp,*)
      write(out_unitp,*) 'BEGINNING calc_ind_n'
      write(out_unitp,*) 'i_point',i_point
      write(out_unitp,*) 'ndim,ind0',ndim,ind0
      write(out_unitp,*) 'n',n
      END IF
!---------------------------------------------------------------------

        IF (i_point .EQ. 1) THEN
!         initialisation des ind_n(i)
          DO i=1,ndim
            ind_n(i)=ind0
          END DO
          ind_n(ndim)=ind0-1
        END IF


!       determine les indices (ind_n) en fonction de i_point
        ind_n(ndim) = ind_n(ndim)+1
        i_modif = ndim
        DO i=ndim,1,-1
          IF (ind_n(i) .EQ. (n(i)+ind0) ) THEN
            ind_n(i)=ind0
            IF (i .NE. 0) ind_n(i-1)=ind_n(i-1)+1
            i_modif = i-1
          END IF
        END DO

!----- for debuging --------------------------------------------------
      IF (debug) THEN
      write(out_unitp,*)
      write(out_unitp,*) 'i_modif',i_modif
      write(out_unitp,*) 'ind_n',ind_n
      write(out_unitp,*) 'END calc_ind_n'
      END IF
!---------------------------------------------------------------------

        calc_ind_n = i_modif

        end function calc_ind_n
