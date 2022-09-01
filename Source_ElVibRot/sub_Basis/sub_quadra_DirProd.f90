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
      SUBROUTINE sub_DirProd_basis(basis_DP)
      USE mod_system
      USE mod_nDindex
      USE mod_basis
      IMPLICIT NONE

!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis), intent(inout) :: basis_DP

!---------------------------------------------------------------------
!---------------------------------------------------------------------

!---------------------------------------------------------------------

      integer           :: ib,jb,iq,ibasis,iqi,ibi,nb_basis,ndim,nq
      integer           :: i0,i1
      integer           :: i,iv,nDI_DP

      integer           :: nb_var,symab,symab_ibasis
      logical           :: packed_loc,primitive_done
      real (kind=Rkind), allocatable :: wrho(:)

      integer, allocatable :: tab_nq(:)
      TYPE (basis)         :: basis_temp


!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_DirProd_basis'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) '--------------------------'
        CALL RecWrite_basis(basis_DP)
        write(out_unitp,*) '--------------------------'
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------------------

     IF (basis_DP%SparseGrid_type /= 0) THEN
       write(out_unitp,*) ' ERROR in ',name_sub
       write(out_unitp,*) ' You cannot use ',name_sub,' subroutine with SPARSE GRID'
       write(out_unitp,*) ' Check the Fortran !!'
       STOP
     END IF
!-----------------------------------------------------------------------
     nb_basis                       = basis_DP%nb_basis
     IF (debug) basis_DP%print_info_OF_basisDP = .TRUE.
     !--- check if the basis is complex => STOP -----------------------
     DO i=1,basis_DP%nb_basis
        basis_DP%cplx = basis_DP%tab_Pbasis(i)%Pbasis%cplx
        IF (basis_DP%cplx) EXIT
     END DO



      IF (basis_DP%cplx) THEN
        write(out_unitp,*) ' STOP the basis is complex!!'
        STOP
      END IF

      !--- check if primitive_done=.TRUE. -----------------------
      DO i=1,basis_DP%nb_basis
        IF (.NOT. basis_DP%tab_Pbasis(i)%Pbasis%primitive_done) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' primitive_done is .FALSE. for the ibasis:',i
          STOP
        END IF
      END DO
      basis_DP%primitive_done = .TRUE.

      !write(out_unitp,*) 'basis_DP%tab_Pbasis(:) packed?',(basis_DP%tab_Pbasis(i)%Pbasis%packed,i=1,nb_basis)


      basis_DP%opt_A(:)      = 0
      basis_DP%opt_B(:)      = 0
      basis_DP%opt_Q0(:)     = 0
      basis_DP%opt_scaleQ(:) = 0

      !--- number of basis functions -----------------------------------
      CALL calc_nDindB_ForDP(basis_DP)

      basis_DP%nb_init = basis_DP%nb
      IF (debug) CALL Write_nDindex(basis_DP%nDindB)
      IF (debug) write(out_unitp,*) '      nb',basis_DP%nb
      IF (debug)                                                        &
       write(out_unitp,*) 'Tabder_Qdyn_TO_Qbasis',basis_DP%Tabder_Qdyn_TO_Qbasis
      CALL flush_perso(out_unitp)


      ! set of tab_symab ---------------------------------------
      CALL Set_SymAbelian_OF_BasisDP(basis_DP)
      IF (basis_DP%print_info_OF_basisDP .AND. print_level > -1)        &
                                         write(out_unitp,*) 'symab done'
      CALL flush_perso(out_unitp)

      ! set of nrho ---------------------------------------
      iq = 1
      DO i=1,basis_DP%nb_basis
        ndim = basis_DP%tab_Pbasis(i)%Pbasis%ndim
        IF (ndim == 0) CYCLE
        basis_DP%nrho(iq:iq+ndim-1) = basis_DP%tab_Pbasis(i)%Pbasis%nrho(:)
        iq = iq + ndim
      END DO

      CALL alloc_NParray(tab_nq,[nb_basis],"tab_nq",name_sub)
      DO i=1,basis_DP%nb_basis
        tab_nq(i) = get_nq_FROM_basis(basis_DP%tab_Pbasis(i)%Pbasis)
      END DO
      IF (debug) write(out_unitp,*) 'tab_nq',tab_nq
      CALL dealloc_nDindex(basis_DP%nDindG)
      CALL init_nDindexPrim(nDindex=basis_DP%nDindG,ndim=nb_basis,    &
                            nDsize=tab_nq)

      CALL dealloc_NParray(tab_nq,"tab_nq",name_sub)

      CALL Set_nq_OF_basis(basis_DP,basis_DP%nDindG%max_nDI)
      CALL Basis_Grid_ParamTOBasis_Grid_Param_init(basis_DP%Basis_Grid_Para)

      IF (debug) CALL Write_nDindex(basis_DP%nDindG)

      nq = get_nq_FROM_basis(basis_DP)
      IF (debug) write(out_unitp,*) '      nq or nqaie',nq
      IF (debug) write(out_unitp,*) ' DP basis: grid done'

      !--- for the grid points -----------------------------------------

      !-- Packed the basis if necessary --------------------------------
      !write(out_unitp,*) 'CALL pack_basis from ',name_sub
      !write(out_unitp,*) 'packed?',basis_DP%packed
      CALL pack_basis(basis_DP)
      IF (debug) write(out_unitp,*) ' DP basis: pack done'

      CALL sub_dnGB_TO_dnBB(basis_DP)
      IF (debug) write(out_unitp,*) ' DP basis: dnBB done'


      !-------------------------------------------------
      !- d1b => dnBGG%d1 and  d2b => dnBGG%d2 ------------
      CALL sub_dnGB_TO_dnGG(basis_DP)
      IF (debug) write(out_unitp,*) ' DP basis: dnGG done'
      !-------------------------------------------------


!     - check the overlap matrix -----------------------------
      !write(out_unitp,*) 'CALL check_ortho_basis from ',name_sub
      !write(out_unitp,*) 'check_basis?',basis_DP%check_basis
      CALL check_ortho_basis(basis_DP)
      IF (debug) write(out_unitp,*) ' DP basis: check done'


      !-------------------------------------------------
      !- wrho ------------
      !CALL alloc_NParray(wrho,[nq],"wrho",name_sub)
      !DO iq=1,nq
      !  wrho(iq) = Rec_WrhonD(basis_DP,iq)
      !END DO
      !CALL alloc_NParray(basis_DP%wrho,[nq],"basis_DP%wrho",name_sub)
      !basis_DP%wrho(:) = wrho(:)
      !CALL dealloc_NParray(wrho,"wrho",name_sub)

      !IF (debug) write(out_unitp,*) ' DP basis: wrho done'
      !-------------------------------------------------

!-----------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(basis_DP)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_DirProd_basis

      RECURSIVE SUBROUTINE calc_nDindB_ForDP(basis_DP)
      USE mod_system
      USE mod_nDindex
      USE mod_dnSVM
      USE mod_basis
      IMPLICIT NONE

!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis), intent(inout) :: basis_DP

!---------------------------------------------------------------------
!---------------------------------------------------------------------

!---------------------------------------------------------------------

      integer             :: ib,jb,iq,ibasis,iqi,ibi,nb_basis
      integer             :: i0,i1
      real (kind=Rkind)   :: weight(basis_DP%nb_basis)
      integer             :: tab_nb(basis_DP%nb_basis)

      real (kind=Rkind)   :: Norm
      TYPE (Type_nDindex) :: nDindB_loc
      integer             :: i,ndim,iv,nDI_DP

      integer             :: ibb,nnb,nb2,newnb2,nb_var,ndim_OF_ib,ivec,nb_vec,val_vec
      integer, allocatable    :: nDval_ref(:)
      integer, allocatable    :: nDval(:)


!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='calc_nDindB_ForDP'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

      nb_basis = basis_DP%nb_basis

      !--- for the basis functions -------------------------------------
      CALL dealloc_nDindex(basis_DP%nDindB)
      CALL alloc_nDindex(basis_DP%nDindB,ndim=nb_basis)

      IF (basis_DP%Type_OF_nDindB == 0) THEN
        CALL alloc_array(basis_DP%nDindB%tab_nDNorm,[nb_basis],     &
                        'basis_DP%nDindB%tab_nDNorm',name_sub)

        DO ib=1,nb_basis

          CALL sort_basis(basis_DP%tab_Pbasis(ib)%Pbasis)
          !CALL Write_nDindex(basis_DP%tab_Pbasis(ib)%Pbasis%nDindB)

          ndim_OF_ib = basis_DP%tab_Pbasis(ib)%Pbasis%nDindB%max_nDI
          IF (debug) write(out_unitp,*) 'ib,ndim_OF_ib',ib,ndim_OF_ib
          CALL alloc_dnSVM(basis_DP%nDindB%tab_nDNorm(ib),              &
                                                  nb_var_vec=ndim_OF_ib)
          basis_DP%nDindB%tab_nDNorm(ib)%d0(:) =                        &
                       basis_DP%tab_Pbasis(ib)%Pbasis%nDindB%Tab_Norm(:)
          IF (debug) write(out_unitp,*) 'ib,tab_nDNorm(ib)%d0',ib,basis_DP%nDindB%tab_nDNorm(ib)%d0(:)
          CALL flush_perso(out_unitp)

        END DO

      END IF

      DO ib=1,nb_basis
        IF (allocated(basis_DP%tab_Pbasis(ib)%Pbasis%nDindB%nDweight)) THEN
          weight(ib) = basis_DP%tab_Pbasis(ib)%Pbasis%nDindB%nDweight(1)
        ELSE
          weight(ib) = ONE
        END IF
        tab_nb(ib) = basis_DP%tab_Pbasis(ib)%Pbasis%nb
      END DO

      CALL init_nDindexPrim(nDindex=basis_DP%nDindB,ndim=nb_basis,      &
                            Type_OF_nDindex=basis_DP%Type_OF_nDindB,    &
                            MaxNorm=basis_DP%Norm_OF_nDindB,            &
                            nDsize=tab_nb,nDweight=weight,              &
                            MaxCoupling=basis_DP%MaxCoupling_OF_nDindB)

      IF (debug) CALL Write_nDindex(basis_DP%nDindB,'nDindB: ')

      basis_DP%nb      = basis_DP%nDindB%max_nDI
      IF (debug) write(out_unitp,*) 'basis_DP%nb',basis_DP%nb
      CALL flush_perso(out_unitp)


      ! set-up Tab_OF_Tabnb2(:) for SparseBasis and type_OF_nDindex=0
      IF (.NOT. associated(basis_DP%Tab_OF_Tabnb2) ) THEN
        CALL alloc_array(basis_DP%Tab_OF_Tabnb2,[nb_basis],                     &
                        "basis_DP%Tab_OF_Tabnb2",name_sub)
      END IF


      IF (basis_DP%nDindB%Type_OF_nDindex == 0) THEN

        DO ibasis=1,nb_basis
          IF (ibasis == 1) THEN
            nnb = 1
            CALL alloc_IntVec(basis_DP%Tab_OF_Tabnb2(ibasis),nnb)
            basis_DP%Tab_OF_Tabnb2(ibasis)%vec(1) = maxval(basis_DP%nDindB%Tab_nDval(ibasis,:))

            nb2 = basis_DP%tab_Pbasis(ibasis)%Pbasis%nb
            newnb2 = count(basis_DP%nDindB%Tab_nDNorm(ibasis)%d0(:) <= basis_DP%nDindB%MaxNorm)
            !IF (nb2 /= newnb2) STOP 'ERROR nb2 /= newnb2'

          ELSE
            nnb = sum(basis_DP%Tab_OF_Tabnb2(ibasis-1)%vec(:))
            CALL alloc_IntVec(basis_DP%Tab_OF_Tabnb2(ibasis),nnb)

            CALL alloc_NParray(nDval,    [ibasis-1],"nDval",    name_sub)
            CALL alloc_NParray(nDval_ref,[ibasis-1],"nDval_ref",name_sub)

            ivec    = 1
            !val_vec = basis_DP%nDindB%Tab_nDval(ibasis-1,1)
            nDval_ref(:) = basis_DP%nDindB%Tab_nDval(1:ibasis-1,1)
            DO nDI_DP=1,basis_DP%nDindB%max_nDI
              nDval(:) = basis_DP%nDindB%Tab_nDval(1:ibasis-1,nDI_DP)
              IF (sum(abs(nDval-nDval_ref)) /= 0) THEN
                nDval_ref(:) = nDval(:)
                ivec = ivec + 1
              END IF
              basis_DP%Tab_OF_Tabnb2(ibasis)%vec(ivec) =                 &
                                basis_DP%nDindB%Tab_nDval(ibasis,nDI_DP)
             !write(out_unitp,*) ' ibasis,ivec,shape vec',ibasis,ivec,shape(basis_DP%Tab_OF_Tabnb2(ibasis)%vec)
            END DO

            CALL dealloc_NParray(nDval,    "nDval",    name_sub)
            CALL dealloc_NParray(nDval_ref,"nDval_ref",name_sub)

          END IF
          !write(out_unitp,*) 'ibasis,size Tab_OF_Tabnb2(ibasis)%vec',ibasis,size(basis_DP%Tab_OF_Tabnb2(ibasis)%vec)
          !write(out_unitp,*) 'ibasis,Tab_OF_Tabnb2(ibasis)%vec',ibasis,basis_DP%Tab_OF_Tabnb2(ibasis)%vec
          IF ((basis_DP%print_info_OF_basisDP .AND. print_level > -1) .OR. debug) &
                             write(out_unitp,'(a,i4,":",2i6)')          &
                             'ibasis,size(vec),sum(vec)',ibasis,        &
                             basis_DP%Tab_OF_Tabnb2(ibasis)%nb_var_vec, &
                             sum(basis_DP%Tab_OF_Tabnb2(ibasis)%vec)
          CALL flush_perso(out_unitp)
        END DO

      ELSE   ! for true direct-product (not SparseBasis)

        nnb = basis_DP%nb
        DO ibasis=nb_basis,1,-1

          nb2 = basis_DP%tab_Pbasis(ibasis)%Pbasis%nb
          nnb = nnb / nb2
          CALL alloc_IntVec(basis_DP%Tab_OF_Tabnb2(ibasis),nnb)
          basis_DP%Tab_OF_Tabnb2(ibasis)%vec(:) = nb2

          IF ((basis_DP%print_info_OF_basisDP .AND. print_level > -1).OR. debug) &
                             write(out_unitp,'(a,i4,":",2i6)')                   &
                             'ibasis,size(vec),sum(vec)',ibasis,        &
                             basis_DP%Tab_OF_Tabnb2(ibasis)%nb_var_vec, &
                             sum(basis_DP%Tab_OF_Tabnb2(ibasis)%vec)
          CALL flush_perso(out_unitp)
        END DO
      END IF

      ! verification ....
      DO ibasis=1,nb_basis
        IF (ibasis == nb_basis) THEN
          nnb = basis_DP%nb
        ELSE
          nnb = basis_DP%Tab_OF_Tabnb2(ibasis+1)%nb_var_vec
        END IF
        IF (sum(basis_DP%Tab_OF_Tabnb2(ibasis)%vec) /= nnb) THEN
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) 'Wrong Tab_OF_Tabnb2 !!'
          write(out_unitp,*) 'ibasis,Tab_OF_Tabnb2',ibasis,                     &
                                 basis_DP%Tab_OF_Tabnb2(ibasis)%vec
          write(out_unitp,*) 'sum(Tab_OF_Tabnb2(ibasis)%vec)',                  &
                              sum(basis_DP%Tab_OF_Tabnb2(ibasis)%vec)
          write(out_unitp,*) 'Tab_OF_Tabnb2(ibasis+1)%nb_var_vec or basis_DP%nb', &
                         nnb
          write(out_unitp,*) 'weight(:)',basis_DP%nDindB%nDweight(:)
          write(out_unitp,*) 'weight(:)',(basis_DP%tab_Pbasis(i)%Pbasis%nDindB%nDweight(1),i=1,nb_basis)
          STOP
        END IF
      END DO

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

      END SUBROUTINE calc_nDindB_ForDP
