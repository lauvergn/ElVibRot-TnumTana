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
!      determination des tous les Ln(xi)=serie_fourier(n,i)
!      + les derivees 1er et 2d.
!
!=============================================================
      SUBROUTINE sub_quadra_Wigner(base)
      USE mod_system
      USE mod_nDindex
      USE mod_basis
      IMPLICIT NONE

!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis), intent(inout)      :: base


!---------------------------------------------------------------------
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      real (kind=Rkind), allocatable :: xm(:)
      real (kind=Rkind), allocatable :: wm(:)
      real (kind=Rkind), allocatable :: xk(:)
      real (kind=Rkind), allocatable :: wk(:)
      real (kind=Rkind), allocatable :: xj(:)
      real (kind=Rkind), allocatable :: wj(:)

      integer           :: max_j,max_k,max_m,max_qj,max_qm,max_qk
      integer           :: j,k,m
      integer           :: ib,iq,nb,nq
      real (kind=Rkind) :: xq(3)
      real (kind=Rkind) :: d0,d1(3),d2(3,3)

      integer           :: nio


!---------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      integer :: err_io,err_mem,memory
      character (len=*), parameter :: name_sub='sub_quadra_Wigner'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------
       nq = get_nq_FROM_basis(base)

       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'nb,nq',base%nb,nq
         write(out_unitp,*) 'L_SparseBasis',base%L_SparseBasis
         write(out_unitp,*) 'L_SparseGrid',base%L_SparseGrid
       END IF
!-----------------------------------------------------------
      base%packed            = .TRUE.
      base%packed_done       = .TRUE.

     IF (base%nb <= 0 .AND. base%Norm_OF_nDindB < 0) STOP 'ERROR nb<=0'

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!      nb and nq
!----------------------------------------------------------------------------
      nb = base%nb
      write(out_unitp,*) '    Basis: Wigner (3D-rotation)'
      write(out_unitp,*) '      old nb_Wigner',nb

      IF (base%xPOGridRep_done) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'xPOGridRep_done=t and a 3D-basis is impossible'
        write(out_unitp,*) 'CHECK the source'
        STOP
      END IF

      IF (nb > 0) THEN
        max_j = int(sqrt(real(nb,kind=Rkind)))-1
        IF ( (max_j+1)**2 .NE. nb ) max_j = max_j+1
        nb = (max_j+1)**2
      ELSE
        max_j = Get_nb_FROM_l_OF_PrimBasis(base%L_SparseBasis,base)-1
        nb = (max_j+1)**2
      END IF
      max_m = 2*max_j+1
      max_k = 2*max_j+1

      write(out_unitp,*) '      new nb_Wigner (without symmetry): ',nb
      write(out_unitp,*) '      max_j,max_k,max_m',max_j,max_k,max_m

      IF (base%L_SparseGrid > -1) THEN
         max_qj = Get_nq_FROM_l_OF_PrimBasis(base%L_SparseGrid,base)-1
         max_qj = max(max_qj,max_j)
      ELSE
         max_qj = -1
         max_qj = max_j
      END IF
      write(out_unitp,*) 'max_qj',max_qj


      write(out_unitp,*) '      old nb_quadra',nq

      max_qj = int(sqrt(real(nq,kind=Rkind)))
      IF (max_qj <= max_j) max_qj = max_j+1
      max_qm = 2*max_qj+1
      max_qk = max_qm

      nq = max_qj*max_qm*max_qk

      write(out_unitp,*) '      new nb_quadra',nq
      write(out_unitp,*) '      max_qj,max_qk,max_qm',max_qj,max_qk,max_qm

      CALL Set_nq_OF_basis(base,nq)

!----------------------------------------------------------------------------

      CALL alloc_xw_OF_basis(base)

      CALL alloc_NParray(xm,[max_qm],'xm',name_sub)
      CALL alloc_NParray(wm,[max_qm],'wm',name_sub)
      CALL alloc_NParray(xk,[max_qk],'xk',name_sub)
      CALL alloc_NParray(wk,[max_qk],'wk',name_sub)
      CALL alloc_NParray(xj,[max_qj],'xj',name_sub)
      CALL alloc_NParray(wj,[max_qj],'wj',name_sub)
      !--------------------------------------------------------------
      !----- grid and weight calculation ----------------------------
      CALL gauss_fourier(xk,wk,max_qk)
      CALL gauss_fourier(xm,wm,max_qm)
      CALL gauleg(-ONE,ONE,xj,wj,max_qj)

      !- tranformation of cos(th) in th --------------------------------
      wj(:) = wj(:)/sqrt(ONE-xj**2)
      xj(:) = acos(xj)

      iq = 0
      DO j=1,max_qj
      DO k=1,max_qk
      DO m=1,max_qm
          iq = iq + 1
          base%x(1,iq)  = xj(j)
          base%x(2,iq)  = xk(k)
          base%x(3,iq)  = xm(m)
          base%w(iq)    = wj(j)*wm(k)*wm(m)
          base%rho(iq)  = sin(xj(j))
          base%wrho(iq) = base%w(iq) * base%rho(iq)
      END DO
      END DO
      END DO

      IF (debug) THEN
        write(out_unitp,*) 'grid for the Ylm. nq:',nq
        DO iq=1,nq
          write(out_unitp,*) base%x(:,iq)
        END DO
      END IF

      !----------------------------------------------------------------
      ! nrho(1) has to be changed
      base%nrho(1) = 2  ! we are working with theta, rho(theta)=sin(theta) instead of one

      !----------------------------------------------------------------
      CALL alloc_dnb_OF_basis(base)

      CALL dealloc_nDindex(base%nDindB)
      base%nDindB%packed = .TRUE.

      IF (base%With_L) THEN
        STOP 'With_L not yet!'
      ELSE
        CALL init_nDindexPrim(base%nDindB,ndim=1,nDsize=[base%nb],  &
                              nDweight=[base%weight_OF_nDindB]     )
      END IF

      ib  = 0
      DO j = 0,max_j
      DO k = 1,2*j+1
      DO m = 1,2*j+1
        ! here a test on j,k,m !!!!
        ib = ib + 1
        base%tab_ndim_index(:,ib) = [j,k,m]
        CALL Set_symabOFSymAbelian_AT_ib(base%P_SymAbelian,ib,-1)

        ! for nDindB
        base%nDindB%Tab_Norm(ib) = j*base%weight_OF_nDindB

        DO iq=1,nq
          xq(:) = base%x(:,iq)
          d0 = ZERO ; d1 = ZERO ; d2 = ZERO
          !CALL d0d1d2Ylm(d0,d1,d2,xq,ibb,num,step)
          base%dnRGB%d0(iq,ib)     = d0
          base%dnRGB%d1(iq,ib,:)   = d1(:)
          base%dnRGB%d2(iq,ib,:,:) = d2(:,:)
        END DO
      END DO
      END DO
      END DO

      IF (ib /= base%nb) STOP 'ib /= nb : ERROR in sub_quadra_Wigner'

      CALL dealloc_NParray(xm,'xm',name_sub)
      CALL dealloc_NParray(wm,'wm',name_sub)
      CALL dealloc_NParray(xk,'xk',name_sub)
      CALL dealloc_NParray(wk,'wk',name_sub)
      CALL dealloc_NParray(xj,'xj',name_sub)
      CALL dealloc_NParray(wj,'wj',name_sub)
!-----------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(base)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

END SUBROUTINE sub_quadra_Wigner
