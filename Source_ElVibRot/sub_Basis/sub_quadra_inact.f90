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
!      determination des tous les Hn(xi)=poly_h(n,i)
!      + les derivees 1er et 2d.
!
!=============================================================
      SUBROUTINE sub_quadra_inact(Basis2n,mole)
      USE mod_system
      USE mod_nDindex
      use mod_Coord_KEO, only: CoordType
      USE mod_basis
      IMPLICIT NONE

!----- variables for the inactive namelist ----------------
      TYPE (CoordType)     :: mole

!----- for the inactive basis sets ----------------------------------
      TYPE (Basis) :: Basis2n

!---------------------------------------------------------------------
      logical               :: num
      real (kind=Rkind)     :: step
      integer               :: i,i0,i1,iQ,nb,nq
      TYPE (Type_nDindex)   :: nDindB

      integer               :: tab_nq(Basis2n%nb_basis)


!---------------------------------------------------------------------
      character (len=*), parameter :: name_sub='sub_quadra_inact'
      integer :: err_mem,memory
!      logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) ' SparseGrid_type :',Basis2n%SparseGrid_type
        write(out_unitp,*) ' nb_inact2n or nb_basis',Basis2n%nb_basis
        write(out_unitp,*) ' tab_nq',Basis2n%nDindG%nDsize(:)
      END IF
!---------------------------------------------------------------------
      num              = .FALSE.
      step             = ZERO
      tab_nq(:)        = Basis2n%nDindG%nDsize(:)


      ! parameter for Basis2n%tab_Pbasis
      CALL alloc_tab_Pbasis_OF_basis(Basis2n)

      DO i=1,Basis2n%nb_basis
        nb = maxval(Basis2n%nDindB%Tab_nDval(i,:))+1
        nq = tab_nq(i)

        Basis2n%tab_Pbasis(i)%Pbasis%primitive_done    = .TRUE.
        Basis2n%tab_Pbasis(i)%Pbasis%ndim              = 1
        Basis2n%tab_Pbasis(i)%Pbasis%name              = "Hm"
        Basis2n%tab_Pbasis(i)%Pbasis%type              = 20
        Basis2n%tab_Pbasis(i)%Pbasis%nb                = nb

        CALL Set_nq_OF_basis(Basis2n%tab_Pbasis(i)%Pbasis,nq)
        Basis2n%tab_Pbasis(i)%Pbasis%nbc               = nb
        Basis2n%tab_Pbasis(i)%Pbasis%nqc               = nq
        Basis2n%tab_Pbasis(i)%Pbasis%contrac           = .FALSE.

        CALL alloc_init_basis(Basis2n%tab_Pbasis(i)%Pbasis)

        CALL sub_quadra_hermite(Basis2n%tab_Pbasis(i)%Pbasis,-1,num,step)
        tab_nq(i) = get_nq_FROM_basis(Basis2n%tab_Pbasis(i)%Pbasis)
        Basis2n%tab_Pbasis(i)%Pbasis%A(:)              = ZERO
        Basis2n%tab_Pbasis(i)%Pbasis%B(:)              = ZERO
        Basis2n%tab_Pbasis(i)%Pbasis%Q0(:)             = ZERO
        Basis2n%tab_Pbasis(i)%Pbasis%scaleQ(:)         = ONE


        CALL init_nDindexPrim(Basis2n%tab_Pbasis(i)%Pbasis%nDindB,1,(/ nb /))
        Basis2n%tab_Pbasis(i)%Pbasis%nDindB%nDweight(1)     = ONE
        Basis2n%tab_Pbasis(i)%Pbasis%nDindB%Type_OF_nDindex = 1
        Basis2n%tab_Pbasis(i)%Pbasis%nDindB%MaxNorm         = ZERO

        iQ = mole%nb_act1 + i
        Basis2n%tab_Pbasis(i)%Pbasis%iQdyn(1) =                         &
                                  mole%ActiveTransfo%list_QactTOQdyn(iQ)
        CALL alloc_NParray(Basis2n%tab_Pbasis(i)%Pbasis%Tabder_Qdyn_TO_Qbasis, &
             (/mole%nb_var/),"...%Tabder_Qdyn_TO_Qbasis",name_sub,(/0/))
        Basis2n%tab_Pbasis(i)%Pbasis%Tabder_Qdyn_TO_Qbasis(:) = 0
        Basis2n%tab_Pbasis(i)%Pbasis%Tabder_Qdyn_TO_Qbasis(iQ) = 1

      END DO

      ! parameter for Basis2n
      Basis2n%type                  = 1
      Basis2n%name                  = 'direct_prod'
      Basis2n%contrac               = .FALSE.
      Basis2n%ndim                  = Basis2n%nb_basis
      Basis2n%packed                = .FALSE.
      Basis2n%packed_done           = .FALSE.
      Basis2n%primitive_done        = .TRUE.
      Basis2n%check_basis           = .FALSE.
      Basis2n%print_info_OF_basisDP = .FALSE.
      Basis2n%check_nq_OF_basis     = .FALSE.
      Basis2n%MaxCoupling_OF_nDindB = Basis2n%nb_basis


      CALL alloc_init_basis(Basis2n)

      ! transfert of iQdyn, nrho --------------------------------
      i0 = 0
      DO i=1,Basis2n%nb_basis
        i1 = i0 + Basis2n%tab_Pbasis(i)%Pbasis%ndim
        Basis2n%iQdyn(i0+1:i1) = Basis2n%tab_Pbasis(i)%Pbasis%iQdyn(:)
        Basis2n%nrho(i0+1:i1)  = Basis2n%tab_Pbasis(i)%Pbasis%nrho(:)
        i0 = i1
      END DO
      ! Set Tabder_Qdyn_TO_Qbasis(:)
      CALL alloc_NParray(Basis2n%Tabder_Qdyn_TO_Qbasis,(/ mole%nb_var /), &
                        "Basis2n%Tabder_Qdyn_TO_Qbasis",name_sub,(/0/))
      Basis2n%Tabder_Qdyn_TO_Qbasis(:) = 0
      DO i=1,Basis2n%ndim
        Basis2n%Tabder_Qdyn_TO_Qbasis(Basis2n%iQdyn(i)) = i
      END DO


      ! table of basis functions
      nDindB = Basis2n%nDindB
      CALL sub_DirProd_basis(Basis2n)
      Basis2n%nDindB = nDindB
      CALL dealloc_nDindex(nDindB)
      Basis2n%nDindB%Tab_nDval(:,:) = Basis2n%nDindB%Tab_nDval(:,:) + 1

      Basis2n%nb     = Basis2n%nDindB%max_nDI


      IF (debug) THEN
        write(out_unitp,*) 'Basis2n functions:',Basis2n%nb
        CALL Write_nDindex(Basis2n%nDindB)
        CALL flush_perso(out_unitp)
      END IF


      ! table of grid points
      CALL dealloc_nDindex(Basis2n%nDindG)
      CALL init_nDindexPrim(nDindex=Basis2n%nDindG,ndim=Basis2n%nb_basis,&
                            type_OF_nDindex=0,nDsize=tab_nq)

      IF (debug) THEN
        write(out_unitp,*) 'Basis2n grids:',get_nq_FROM_basis(Basis2n)
        CALL Write_nDindex(Basis2n%nDindG)
        CALL flush_perso(out_unitp)
      END IF

      IF (Basis2n%nDindG%max_nDI /= get_nq_FROM_basis(Basis2n)) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' nq /= max_nDI'
        write(out_unitp,*) ' ... nq:     ',get_nq_FROM_basis(Basis2n)
        write(out_unitp,*) ' ... max_nDI:',Basis2n%nDindG%max_nDI
        STOP
      END IF


      write(out_unitp,*) 'Basis2n%iQdyn',Basis2n%iQdyn(:)
      write(out_unitp,*) 'Basis2n%Tabder_Qdyn_TO_Qbasis',Basis2n%Tabder_Qdyn_TO_Qbasis(:)
!----------------------------------------------------------------------------


!---------------------------------------------------------------------
      IF (debug) THEN
        DO i=1,Basis2n%nb_basis
          write(out_unitp,*) 'basis (active order)',i
          write(out_unitp,*)
          CALL RecWrite_basis(Basis2n%tab_Pbasis(i)%Pbasis)
        END DO
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------

      end subroutine sub_quadra_inact

