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
      TYPE (basis), intent(in)      :: base


!---------------------------------------------------------------------
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      integer       :: nb,ib,ibb,iq,iqm,iql,nbl,ibl,ibm
      integer       :: symab,symab_l,symab_mfourier,m,Read_symab

      real (kind=Rkind), allocatable :: xm(:)
      real (kind=Rkind), allocatable :: wm(:)
      real (kind=Rkind), allocatable :: xl(:)
      real (kind=Rkind), allocatable :: wl(:)

      integer           :: max_l,max_m,max_ql,max_qm,nq,ibbl
      real (kind=Rkind) :: xq(2)
      real (kind=Rkind) :: d0,d1(2),d2(2,2)
      real (kind=Rkind) :: s,c

      character (len=3) :: name_i
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
!      test sur nb_fourier et nb_quadra
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
        max_l = int(sqrt(real(nb,kind=Rkind)))-1
        IF ( (max_l+1)**2 .NE. nb ) max_l = max_l+1
        nb = (max_l+1)**2
      ELSE
        max_l = Get_nb_FROM_l_OF_PrimBasis(base%L_SparseBasis,base)-1
        nb = (max_l+1)**2
      END IF
      max_m = 2*max_l+1

      write(out_unitp,*) '      new nb_Wigner (without symmetry): ',nb
      write(out_unitp,*) '      max_l, max_m',max_l,max_m

      IF (base%L_SparseGrid > -1) THEN
         max_ql = Get_nq_FROM_l_OF_PrimBasis(base%L_SparseGrid,base)-1
         max_ql = max(max_ql,max_l)
      ELSE
         max_ql = -1
         max_ql = max_l
      END IF
      write(out_unitp,*) 'max_ql',max_ql


      IF (lebedev < 0) THEN
        write(out_unitp,*) '      old nb_quadra',nq
        write(out_unitp,*) '  not lebedev',lebedev

        max_ql = int(sqrt(real(nq,kind=Rkind)))
        IF (max_ql <= max_l) max_ql = max_l+1
        max_qm = 2*max_ql+1
        IF (mod(max_qm,2) .EQ. 1) max_qm = max_qm + 1
        nq = max_ql*max_qm

        write(out_unitp,*) '      new nb_quadra',nq
        write(out_unitp,*) '      max_ql, max_qm',max_ql,max_qm
      END IF

      CALL Set_nq_OF_basis(base,nq)

      ! calculation of base%nb  without symmetry
        base%nb = 0
        DO ibl = 0,max_l
        DO ibm = 1,2*ibl+1
          base%nb = base%nb + 1
        END DO
        END DO
      write(out_unitp,*) '      new nb_Ylm (with symmetry): ',base%nb

!----------------------------------------------------------------------------

      CALL alloc_xw_OF_basis(base)

        CALL alloc_NParray(xm,(/ max_qm /),'xm',name_sub)
        CALL alloc_NParray(wm,(/ max_qm /),'wm',name_sub)
        CALL alloc_NParray(xl,(/ max_qm /),'xl',name_sub)
        CALL alloc_NParray(wl,(/ max_qm /),'wl',name_sub)

        !--------------------------------------------------------------
        !----- grid and weight calculation ----------------------------
        CALL gauss_fourier(xm,wm,max_qm)
        CALL gauleg(-ONE,ONE,xl,wl,max_ql)

        !- tranformation cos(th) => th --------------------------------
        DO iql=1,max_ql
          c  = xl(iql)
          s = sqrt(ONE-c*c)
          xl(iql) = acos(c)
          wl(iql) = wl(iql)/s
        END DO

        iq = 0
        DO iql=1,max_ql
        DO iqm=1,max_qm
          iq = iq + 1
          base%x(1,iq)  = xl(iql)
          base%x(2,iq)  = xm(iqm)
          base%w(iq)    = wl(iql)*wm(iqm)
          base%rho(iq)  = sin(xl(iql))
          base%wrho(iq) = base%w(iq) * base%rho(iq)
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
      CALL alloc_SymAbelian(base%P_SymAbelian,base%nb)
      Read_symab = Get_Read_symabOFSymAbelian(base%P_SymAbelian)
      CALL Set_ReadsymabOFSymAbelian(base%P_SymAbelian,Read_symab)

      CALL alloc_dnb_OF_basis(base)

      CALL dealloc_nDindex(base%nDindB)
      base%nDindB%packed = .TRUE.

      IF (base%With_L) THEN
          CALL init_nDindexPrim(base%nDindB,1,(/ base%nb /))
          base%nDindB%With_L      = .TRUE.
          base%nDindB%Tab_L(:)    = -1
          base%nDindB%Tab_Norm(:) = -ONE

          ibb  = 0
          ibbl = -1
          DO ibl = 0,max_l
            IF (isyml >= 0 .AND. mod(ibl,2) /= isyml) CYCLE
            IF (ibl+1 <= base%L_TO_nb%A) THEN
              ibbl = 0
            ELSE
              IF (mod(ibl- base%L_TO_nb%A,base%L_TO_nb%B) == 0) THEN
                ibbl = ibbl + 1
              END IF
          END IF
          DO ibm = 1,2*ibl+1
            ibb = ibb + 1
            base%nDindB%Tab_L(ibb)    = ibbl
            base%nDindB%Tab_Norm(ibb) = real(ibbl,kind=Rkind)
          END DO
          END DO
          !write(out_unitp,*) 'base%nDindB%Tab_L',base%nDindB%Tab_L
          !STOP 'With_L'
      ELSE
        CALL init_nDindexPrim(base%nDindB,ndim=1,nDsize=(/ base%nb /),  &
                              nDweight=(/ base%weight_OF_nDindB /)     )
      END IF

      ib  = 0
      ibb = 0
      DO ibl = 0,max_l
      DO ibm = 1,2*ibl+1
        ibb = ibb + 1
        IF (isyml >= 0 .AND. mod(ibl,2) /= isyml) CYCLE
        ib = ib + 1  ! for the symmetry
        base%tab_ndim_index(:,ib) = (/ ibl,ibm /)

        ! for symab (test)
        !m              = ibm/2
        !symab_l        = 0
        !symab_mfourier = 0
        !IF (mod(ibl-m,2) == 1) symab_l        = 2 ! odd, Plm
        !IF (mod(ibm,2) == 0)   symab_mfourier = 1 ! odd, m (the true m)
        !symab = Calc_symab1_EOR_symab2(symab_l,symab_mfourier)
        !CALL Set_symabOFSymAbelian_AT_ib(base%P_SymAbelian,ib,symab)
        !write(out_unitp,*) 'ib,symab_l,symab_mfourier',ib,symab_l,symab_mfourier
        !write(out_unitp,*) 'ib,symab',ib,symab

        SELECT CASE (Read_symab)
        CASE (-1)
          CALL Set_symabOFSymAbelian_AT_ib(base%P_SymAbelian,ib,-1)
        CASE (0,1,2,3,4,5,6,7)
          IF (mod(ibl,2) == 0) THEN
            CALL Set_symabOFSymAbelian_AT_ib(base%P_SymAbelian,ib,0)
          ELSE
            CALL Set_symabOFSymAbelian_AT_ib(base%P_SymAbelian,ib,Read_symab)
          END IF
        CASE DEFAULT
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  it should never append. The error should come from'
          write(out_unitp,*) ' "Set_ReadsymabOFSymAbelian" subroutine'
          write(out_unitp,*) ' CHECK the fortran!!'
          STOP
        END SELECT

        ! for nDindB
        base%nDindB%Tab_Norm(ib) = real(ibl,kind=Rkind)*base%weight_OF_nDindB

        DO iq=1,nq
          xq(:) = base%x(:,iq)
          CALL d0d1d2Ylm(d0,d1,d2,xq,ibb,num,step)
          base%dnRGB%d0(iq,ib)     = d0
          base%dnRGB%d1(iq,ib,:)   = d1(:)
          base%dnRGB%d2(iq,ib,:,:) = d2(:,:)
        END DO
      END DO
      END DO

      CALL Set_nbPERsym_FROM_SymAbelian(base%P_SymAbelian)

      IF (debug) CALL Write_SymAbelian(base%P_SymAbelian)

      IF (ib /= base%nb) STOP 'ib /= nb : ERROR in sub_quadra_Wigner'

      IF (lebedev < 0) THEN
        CALL dealloc_NParray(xm,'xm',name_sub)
        CALL dealloc_NParray(wm,'wm',name_sub)
        CALL dealloc_NParray(xl,'xl',name_sub)
        CALL dealloc_NParray(wl,'wl',name_sub)
      END IF
!-----------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(base)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

END SUBROUTINE sub_quadra_Wigner
