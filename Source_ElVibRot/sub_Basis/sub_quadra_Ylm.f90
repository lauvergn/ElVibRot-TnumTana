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

!=============================================================
!
!      determination des tous les Ln(xi)=serie_fourier(n,i)
!      + les derivees 1er et 2d.
!
!=============================================================
      SUBROUTINE sub_quadra_Ylm(base,isyml,isymm)
      USE mod_system
      USE mod_nDindex
      USE mod_basis
      IMPLICIT NONE

!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis)      :: base
      logical           :: num
      real (kind=Rkind) :: step
      integer           :: isyml,isymm


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

      integer           :: Lebedev = 1
      character (len=3) :: name_i
      integer           :: nio
      TYPE (param_file) :: Lebedev_file
      real (kind=Rkind) :: theta,phi,thetaR


!---------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      integer :: err_io,err_mem,memory
      character (len=*), parameter :: name_sub='sub_quadra_Ylm'
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
       num   = .FALSE.

      base%packed            = .TRUE.
      base%packed_done       = .TRUE.


     IF (base%nb <= 0 .AND. base%Norm_OF_nDindB < 0) STOP 'ERROR nb<=0'

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!      test sur nb_fourier et nb_quadra
!----------------------------------------------------------------------------
      nb = base%nb
      write(out_unitp,*) '    Basis: Ylm'
      write(out_unitp,*) '      old nb_Ylm',nb

      IF (base%xPOGridRep_done) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'xPOGridRep_done=t and a 2D-basis is impossible'
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

      write(out_unitp,*) '      new nb_Ylm (without symmetry): ',nb
      write(out_unitp,*) '      max_l, max_m',max_l,max_m

      IF (base%L_SparseGrid > -1) THEN
         max_ql = Get_nq_FROM_l_OF_PrimBasis(base%L_SparseGrid,base)-1
         max_ql = max(max_ql,max_l)
      ELSE
         max_ql = -1
         max_ql = max_l
      END IF
      write(out_unitp,*) 'max_ql',max_ql

      lebedev = 2*max_ql+1
      !lebedev = -1

      IF (lebedev >= 0) THEN

        write(out_unitp,*) '      old nb_quadra',nq
        write(out_unitp,*) '      max_ql',max_ql
        write(out_unitp,*) '      lebedev',lebedev

        DO
          write(name_i,'(i3)') lebedev
          IF (name_i(1:1) == ' ') name_i(1:1) = '0'
          IF (name_i(2:2) == ' ') name_i(2:2) = '0'

          Lebedev_file%name = trim(EVRT_path) //                        &
           '/Internal_data/Lebedev-grid/lebedev_' // trim(name_i) // '.txt'
          !write(out_unitp,*) 'Lebedev_file%name: ',Lebedev_file%name
          CALL file_open(Lebedev_file,nio,old=.TRUE.,err_file=err_io)
          IF (err_io == 0) THEN
            write(out_unitp,*) ' Lebedev parameter: ',lebedev
            ! first the grid file and then the number of grid points
            iq = 0
            DO
              read(nio,*,iostat=err_io)
              IF (err_io /= 0) EXIT
              iq=iq+1
            END DO
            nq = iq
            EXIT
          ELSE
            lebedev = lebedev + 2
            IF (lebedev > 131) THEN
              write(out_unitp,*) 'WARNING in ',name_sub
              write(out_unitp,*) 'the lebedev parameter is too large',lebedev
              write(out_unitp,*) ' => Do not use sparse grid (L_SparseGrid=-1)'
              lebedev = -1
              EXIT
            END IF
          END IF
        END DO
        close(nio)
        write(out_unitp,*) '      new nb_quadra',nq
      END IF

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

      ! calculation of base%nb with symmetry
      write(out_unitp,*) 'isyml,isymm',isyml,isymm
      IF (isyml >= 0) THEN
        base%nb = 0
        DO ibl = 0,max_l
        IF (mod(ibl,2) /= isyml) CYCLE
        write(out_unitp,*) 'ibl',ibl
        DO ibm = 1,2*ibl+1
          base%nb = base%nb + 1
        END DO
        END DO
      ELSE ! without symmetry
        base%nb = 0
        DO ibl = 0,max_l
        DO ibm = 1,2*ibl+1
          base%nb = base%nb + 1
        END DO
        END DO
      END IF
      write(out_unitp,*) '      new nb_Ylm (with symmetry): ',base%nb

!----------------------------------------------------------------------------

      CALL alloc_xw_OF_basis(base)

      IF (lebedev < 0) THEN
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
      ELSE
        CALL file_open(Lebedev_file,nio)

        ! rotation of all points of 0.02 rad in theta to avoid both poles (0 and pi).
        ! and also of 0.02 rad in phi to avoid the (0 and pi)
        !CALL random_number(thetaR)
        !thetaR = thetaR/FIVE ! now thetaR E [-0.2:0.2]
        thetaR = TWO/TEN**2
        write(out_unitp,*) 'thetaR',thetaR
        DO iq=1,nq
          read(nio,*,iostat=err_io) phi,theta,base%wrho(iq)
          IF (err_io /= 0) STOP 'iq > nq STOP in Ylm'
          !write(out_unitp,*) 'iq,theta,phi,w,(ori)',iq,theta,phi, base%wrho(iq)

          theta = theta / 180._Rkind*pi
          phi   = phi   / 180._Rkind*pi

          base%x(1,iq) = acos(cos(theta)*cos(thetaR)-cos(phi)*sin(theta)*sin(thetaR))
          base%x(2,iq) = thetaR + atan2(cos(phi)*cos(thetaR)*sin(theta)+cos(theta)*sin(thetaR), sin(phi)*sin(theta))

          base%wrho(iq) = base%wrho(iq) * FOUR*pi
          base%rho(iq)  = sin(base%x(1,iq))
          base%w(iq)    = base%wrho(iq) / base%rho(iq)

          !write(out_unitp,*) 'iq,theta,phi,w',iq,base%x(:,iq) * 180._Rkind/pi, base%wrho(iq)

        END DO

        close(nio)
      END IF
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
          !write(6,*) 'base%nDindB%Tab_L',base%nDindB%Tab_L
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
        !write(6,*) 'ib,symab_l,symab_mfourier',ib,symab_l,symab_mfourier
        !write(6,*) 'ib,symab',ib,symab

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

      CALL Write_SymAbelian(base%P_SymAbelian)

      IF (ib /= base%nb) STOP 'ib /= nb : ERROR in sub_quadra_Ylm'

      IF (lebedev < 0) THEN
        CALL dealloc_NParray(xm,'xm',name_sub)
        CALL dealloc_NParray(wm,'wm',name_sub)
        CALL dealloc_NParray(xl,'xl',name_sub)
        CALL dealloc_NParray(wl,'wl',name_sub)
      END IF
!-----------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(base)
        write(out_unitp,*) 'END sub_quadra_Ylm'
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_quadra_Ylm

!=============================================================
!
!      Yj1k1(th1,0)Yj2k2(th2,phi2)
!      k1=-k2
!      + les derivees 1er et 2d.
!
!=============================================================
      SUBROUTINE sub_quadra_ABplusCD(base)
      USE mod_system
      USE mod_basis
      IMPLICIT NONE

!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis)  :: base


!---------------------------------------------------------------------
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      integer       :: nb,ib,iq
      integer       :: il1,il2,im,mmm,iql1,iql2,iqm

      real (kind=Rkind), allocatable :: xm(:)
      real (kind=Rkind), allocatable :: wm(:)
      real (kind=Rkind), allocatable :: xl1(:)
      real (kind=Rkind), allocatable :: wl1(:)
      real (kind=Rkind), allocatable :: xl2(:)
      real (kind=Rkind), allocatable :: wl2(:)

      integer                    :: max_l,max_m,max_ql,max_qm,nq
      real (kind=Rkind) :: xq(3)
      real (kind=Rkind) :: s,c
      real (kind=Rkind) :: d0pl1m,d1pl1m,d2pl1m
      real (kind=Rkind) :: d0pl2m,d1pl2m,d2pl2m
      real (kind=Rkind) :: d0fm,d1fm,d2fm,d3fm

!---------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_quadra_ABplusCD'
      logical,parameter :: debug=.FALSE.
!     logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------
       nq = get_nq_FROM_basis(base)

       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'nb,nq',base%nb,nq
       END IF
!-----------------------------------------------------------

      base%packed            = .TRUE.
      base%packed_done       = .TRUE.

      IF (base%nb <= 0) STOP 'ERROR in sub_quadra_ABplusCD: nb<=0'

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
!      test sur nb_fourier et nb_quadra
!----------------------------------------------------------------------------
      max_l = base%nb
      write(out_unitp,*) '    Basis: Coll AB+CD'
      write(out_unitp,*) '      old max_l',max_l

      IF (base%xPOGridRep_done) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'xPOGridRep_done=t and a 3D-basis is impossible'
        write(out_unitp,*) 'CHECK the source'
        STOP
      END IF

      nb = 0
      DO il1=0,max_l
      DO il2=0,max_l
      DO im=-max_l,max_l
        IF (abs(im) <= il1 .AND. abs(im) <= il2) THEN
          nb = nb + 1
        END IF
      END DO
      END DO
      END DO

      base%nb = nb

      nq = (max_l+1)**2 * (2*max_l+1)
      max_ql  = max_l+1
      max_qm  = 2*max_l+1

      write(out_unitp,*) '      new nb       ',nb
      write(out_unitp,*) '      new nb_quadra',nq
      CALL Set_nq_OF_basis(base,nq)


!----------------------------------------------------------------------------

      CALL alloc_xw_OF_basis(base)

      CALL alloc_NParray(xm, (/ max_qm /),'xm', name_sub)
      CALL alloc_NParray(wm, (/ max_qm /),'wm', name_sub)
      CALL alloc_NParray(xl1,(/ max_qm /),'xl1',name_sub)
      CALL alloc_NParray(wl1,(/ max_qm /),'wl1',name_sub)
      CALL alloc_NParray(xl2,(/ max_qm /),'xl2',name_sub)
      CALL alloc_NParray(wl2,(/ max_qm /),'wl2',name_sub)

!----------------------------------------------------------------------------
!----- grid and weight calculation ------------------------------------------
      CALL gauss_fourier(xm,wm,max_qm)
      xm(:) = xm(:) + pi
      CALL gauleg(-ONE,ONE,xl1,wl1,max_ql)
      CALL gauleg(-ONE,ONE,xl2,wl2,max_ql)


!     - tranformation cos(th) => th -----------------------------------------
      DO iql1=1,max_ql
        c  = xl1(iql1)
        s = sqrt(ONE-c*c)

        xl1(iql1) = acos(c)
        wl1(iql1) = wl1(iql1)/s
        xl2(iql1) = acos(c)
        wl2(iql1) = wl2(iql1)/s
      END DO

      ! nrho(1) has to be changed
      base%nrho(1) = 2  ! it means, the variable, x, is substituted by cos(theta)
      base%nrho(2) = 2  ! it means, the variable, x, is substituted by cos(theta)

      iq = 0
      DO iql1=1,max_ql
      DO iql2=1,max_ql
      DO iqm=1,max_qm
        iq = iq + 1
        base%x(1,iq)  = xl1(iql1)
        base%x(2,iq)  = xl2(iql2)
        base%x(3,iq)  = xm(iqm)

        base%w(iq)    = wl1(iql1) * wl2(iql2) * wm(iqm)
        base%rho(iq)  = sin(xl1(iql1)) * sin(xl2(iql2))

        base%wrho(iq) = base%w(iq) * base%rho(iq)
      END DO
      END DO
      END DO

!----------------------------------------------------------------------------
      CALL alloc_dnb_OF_basis(base)

      ib  = 0
      DO il1 = 0,max_l
      DO il2 = 0,max_l
      DO mmm=1,2*max_l+1
        im = mmm/2
        IF (abs(im) <= il1 .AND. abs(im) <=il2) THEN
          ib = ib + 1
          !write(out_unitp,*) 'il1,il2,im,mmm',il1,il2,im,mmm
          !CALL flush_perso(out_unitp)
          DO iq=1,nq
            xq(:) = base%x(:,iq)

            CALL d0d1d2Plm(d0pl1m,d1pl1m,d2pl1m,xq(1),il1,im)
            CALL d0d1d2Plm(d0pl2m,d1pl2m,d2pl2m,xq(2),il2,im)
            CALL d0d1d2d3fourier(xq(3),d0fm,d1fm,d2fm,d3fm,mmm)


            base%dnRGB%d0(iq,ib)     = d0pl1m * d0pl2m * d0fm

            base%dnRGB%d1(iq,ib,1)   = d1pl1m * d0pl2m * d0fm
            base%dnRGB%d1(iq,ib,2)   = d0pl1m * d1pl2m * d0fm
            base%dnRGB%d1(iq,ib,3)   = d0pl1m * d0pl2m * d1fm

            base%dnRGB%d2(iq,ib,1,1)   = d2pl1m * d0pl2m * d0fm
            base%dnRGB%d2(iq,ib,1,2)   = d1pl1m * d1pl2m * d0fm
            base%dnRGB%d2(iq,ib,1,3)   = d1pl1m * d0pl2m * d1fm

            base%dnRGB%d2(iq,ib,2,1)   = d1pl1m * d1pl2m * d0fm
            base%dnRGB%d2(iq,ib,2,2)   = d0pl1m * d2pl2m * d0fm
            base%dnRGB%d2(iq,ib,2,3)   = d0pl1m * d1pl2m * d1fm

            base%dnRGB%d2(iq,ib,3,1)   = d1pl1m * d0pl2m * d1fm
            base%dnRGB%d2(iq,ib,3,2)   = d0pl1m * d1pl2m * d1fm
            base%dnRGB%d2(iq,ib,3,3)   = d0pl1m * d0pl2m * d2fm

        END DO

        END IF
      END DO
      END DO
      END DO

      IF (ib /= base%nb) STOP 'ib /= nb : ERROR in sub_quadra_ABplusCD'

      CALL dealloc_NParray(xm, 'xm', name_sub)
      CALL dealloc_NParray(wm, 'wm', name_sub)
      CALL dealloc_NParray(xl1,'xl1',name_sub)
      CALL dealloc_NParray(wl1,'wl1',name_sub)
      CALL dealloc_NParray(xl2,'xl2',name_sub)
      CALL dealloc_NParray(wl2,'wl2',name_sub)
!-----------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(base)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_quadra_ABplusCD
