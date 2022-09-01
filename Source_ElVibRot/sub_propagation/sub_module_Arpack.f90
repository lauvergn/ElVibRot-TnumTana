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
MODULE mod_Arpack
USE mod_Constant
IMPLICIT NONE

PRIVATE
PUBLIC :: sub_propagation_Arpack,sub_propagation_Arpack_Sym

CONTAINS

!=======================================================================================
      SUBROUTINE sub_propagation_Arpack(psi,Ene,nb_diago,max_diago,     &
                                        para_H,para_propa)
      USE mod_system
      USE mod_psi,    ONLY : param_psi,alloc_psi,trie_psi,dealloc_psi,  &
                             param_WP0
      USE mod_Op
      USE mod_propa
      USE mod_MPI_aux
      IMPLICIT NONE

      !----- Operator: Hamiltonian ----------------------------
      TYPE (param_Op)   :: para_H

      !----- WP, energy ... -----------------------------------
      TYPE (param_propa)        :: para_propa
      TYPE (param_WP0)          :: para_WP0

      integer                   :: nb_diago,max_diago
      TYPE (param_psi)          :: psi(max_diago)

      real (kind=Rkind)         :: Ene(max_diago)
      logical                   :: With_Grid,With_Basis

      !------ working parameters --------------------------------
      TYPE (param_psi)               :: psi_loc,Hpsi_loc
      real (kind=Rkind)              :: ZPE
      real (kind=Rkind), allocatable :: Evec(:)

!     %-----------------------------%
!     | Define leading dimensions   |
!     | for all arrays.             |
!     | MAXN:   Maximum dimension   |
!     |         of the A allowed.   |
!     | MAXNEV: Maximum NEV allowed |
!     | MAXNCV: Maximum NCV allowed |
!     %-----------------------------%

      integer  (kind=4)      ::   maxn, maxnev, maxncv, ldv

!     %--------------%
!     | Local Arrays |
!     %--------------%

      real(kind=Rkind), allocatable :: v(:,:),workl(:),workd(:),workev(:),d(:,:),resid(:),ax(:)
      logical, allocatable          :: select(:)
      integer  (kind=4)             :: iparam(11), ipntr(14)

!     %---------------%
!     | Local Scalars |
!     %---------------%

      character (len=1)   :: bmat
      character (len=2)   :: which

      integer (kind=4)    :: ido, n, nev, ncv, lworkl, info, ierr, j
      integer (kind=4)    :: nconv, maxitr, mode, ishfts
      logical             :: rvec,first
      logical             :: if_deq0=.FALSE. ! for MPI

      real(kind=Rkind)    :: tol, sigmar, sigmai

      integer (kind=4)    :: it

!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%

      real(kind=Rkind) ::  dnrm2,dlapy2
      external         dnrm2, dlapy2, daxpy

      TYPE(param_file)     :: Log_file
      integer              :: iunit
      logical              :: cplx
      complex (kind=Rkind) :: cplxE

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_propagation_Arpack'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) ' nb_diago',nb_diago
        write(out_unitp,*) ' max_diago',max_diago
        write(out_unitp,*) ' para_Davidson',para_propa%para_Davidson
        write(out_unitp,*)
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

      With_Grid  = para_propa%para_Davidson%With_Grid
      With_Basis = .NOT. para_propa%para_Davidson%With_Grid
      cplx = para_H%cplx
      IF (cplx) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' cplx is not yet possible in this subroutine!'
        STOP
      END IF

      CALL init_psi(psi_loc,para_H,cplx)
      CALL alloc_psi(psi_loc,BasisRep=With_Basis,GridRep=With_Grid)
      CALL init_psi(Hpsi_loc,para_H,cplx)
      CALL alloc_psi(Hpsi_loc,BasisRep=With_Basis,GridRep=With_Grid)

!------ initialization -------------------------------------
      Log_file%name='Arpack.log'
      IF(MPI_id==0) CALL file_open(Log_file,iunit)

      IF (With_Grid) THEN
        n = int(para_H%nb_qaie,kind=4)
      ELSE
        n = int(para_H%nb_tot,kind=4)
      END IF
!#if(run_MPI)
!      CALL MPI_Bcast(n,size1_MPI,MPI_Integer4,root_MPI,MPI_COMM_WORLD,MPI_err)
!#endif
      IF(openmpi .AND. MPI_scheme/=1) CALL MPI_Bcast_(n,size1_MPI,root_MPI)

      IF (nb_diago == 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' nb_diago=0 is not possible with ARPACK'
        STOP
      END IF

      nev =  int(nb_diago,kind=4)
      ncv =  2*nev+10
      ncv =  3*nev+10 
      !was ncv =  2*nev+10
      
      ncv=MIN(N,ncv)  ! prevent infor=-3 error:
                      ! NCV must be greater than NEV and less than or equal to N.

      maxn   = n
      ldv    = maxn
      maxnev = nev
      maxncv = ncv

      allocate(select(maxncv))
      allocate(ax(maxn))
      allocate(d(maxncv,3))
      allocate(resid(maxn)) !< initial guess vector when infor=1

      allocate(v(ldv,maxncv))
      allocate(workd(3*maxn))
      allocate(workev(3*maxncv))
      allocate(workl(3*maxncv*maxncv+6*maxncv))

      bmat  = 'I'
      which = 'SR'
      !which = 'SM'

!     %--------------------------------------------------%
!     | The work array WORKL is used in DSAUPD as        |
!     | workspace.  Its dimension LWORKL is set as       |
!     | illustrated below.  The parameter TOL determines |
!     | the stopping criterion.  If TOL<=0, machine      |
!     | precision is used.  The variable IDO is used for |
!     | reverse communication and is initially set to 0. |
!     | Setting INFO=0 indicates that a random vector is |
!     | generated in DSAUPD to start the Arnoldi         |
!     | iteration.                                       |
!     %--------------------------------------------------%

      lworkl = 3*ncv**2+6*ncv
      tol    = ZERO
      !tol    = ONETENTH**6

      ido    = 0

      info   = 1
      IF (info /= 0) THEN
        IF(keep_MPI) CALL ReadWP0_Arpack(psi_loc,nb_diago,max_diago,                   &
                                            para_propa%para_Davidson,para_H%cplx)
        IF (para_propa%para_Davidson%With_Grid) THEN
          IF(keep_MPI) resid(:) = psi_loc%RvecG(:)
        ELSE
          IF(keep_MPI) resid(:) = psi_loc%RvecB(:)
        END IF
      END IF


      DO j=1,max_diago
        CALL init_psi(psi(j),para_H,cplx)
      END DO
      DO j=1,nb_diago
        CALL alloc_psi(psi(j),BasisRep=With_Basis,GridRep=With_Grid)
      END DO

!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of DSAUPD is used     |
!     | (IPARAM(7) = 1).  All these options may be        |
!     | changed by the user. For details, see the         |
!     | documentation in DSAUPD.                          |
!     %---------------------------------------------------%

      ishfts = 1
      maxitr = int(min(max_diago,para_propa%para_Davidson%max_it),kind=4)
      mode   = 1

      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode

!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%
      it = 0
      DO
        it = it +1
!       %---------------------------------------------%
!       | Repeatedly call the routine DNAUPD and take |
!       | actions indicated by parameter IDO until    |
!       | either convergence is indicated or maxitr   |
!       | has been exceeded.                          |
!       %---------------------------------------------%
#if __ARPACK == 1
        IF(keep_MPI) call dnaupd(ido, bmat, n, which, nev, tol, resid,                 &
                                  ncv, v, ldv, iparam, ipntr, workd, workl, lworkl,    &
                                  info)
!#if(run_MPI)
!        CALL MPI_Bcast(info,size1_MPI,MPI_Integer4,root_MPI,MPI_COMM_WORLD,MPI_err)
!        CALL MPI_Bcast(ido, size1_MPI,MPI_Integer4,root_MPI,MPI_COMM_WORLD,MPI_err)
!#endif
       IF(openmpi .AND. MPI_scheme/=1) THEN
         CALL MPI_Bcast_(info,size1_MPI,root_MPI)
         CALL MPI_Bcast_(ido ,size1_MPI,root_MPI)
       ENDIF

#else
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) ' The ARPACK library is not present!'
        write(out_unitp,*) 'Use Arpack=f and Davidson=t'
        STOP 'ARPACK has been removed'
#endif
        !write(out_unitp,*) 'it,ido',it,ido

        IF (abs(ido) /= 1) EXIT

!       %--------------------------------------%
!       | Perform matrix vector multiplication |
!       |              y <--- OP*x             |
!       | The user should supply his/her own   |
!       | matrix vector multiplication routine |
!       | here that takes workd(ipntr(1)) as   |
!       | the input, and return the result to  |
!       | workd(ipntr(2)).                     |
!       %--------------------------------------%

        CALL sub_OpV1_TO_V2_Arpack(workd(ipntr(1):ipntr(1)-1+n),                       &
                                   workd(ipntr(2):ipntr(2)-1+n),                       &
                                   psi_loc,Hpsi_loc,para_H,cplxE,para_propa,int(n,4))

        IF(MPI_id==0) THEN
          write(iunit,*) 'Arpack <psi H psi>:',                                        &
                  dot_product(workd(ipntr(1):ipntr(1)-1+n),workd(ipntr(2):ipntr(2)-1+n))
          CALL flush_perso(iunit)
        ENDIF
      END DO

      IF(MPI_id==0) THEN
        write(iunit,*) 'End Arpack ' 
        CALL flush_perso(iunit)
        CALL file_close(Log_file)
      ENDIF
!---------------------------------------------------------------------------------------

!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%

      if ( info .lt. 0 ) then
!
!       %--------------------------%
!       | Error message. Check the |
!       | documentation in DSAUPD. |
!       %--------------------------%

        write(out_unitp,*)
        write(out_unitp,*) ' Error with _saupd, info = ', info,' from ', MPI_id
        write(out_unitp,*) ' Check documentation in _naupd '
        write(out_unitp,*) ' '

      else

!       %-------------------------------------------%
!       | No fatal errors occurred.                 |
!       | Post-Process using DSEUPD.                |
!       |                                           |
!       | Computed eigenvalues may be extracted.    |
!       |                                           |
!       | Eigenvectors may also be computed now if  |
!       | desired.  (indicated by rvec = .true.)    |
!       %-------------------------------------------%

        rvec = .true.

#if __ARPACK == 1
        ! consider p-arpark 
        IF(keep_MPI) call dneupd(rvec, 'A', select, d, d(1,2), v, ldv,                 &
                                  sigmar, sigmai, workev, bmat, n, which, nev, tol,    &
                                  resid, ncv, v, ldv, iparam, ipntr, workd, workl,     &
                                  lworkl, ierr )
!#if(run_MPI)
!        CALL MPI_Bcast(ierr,size1_MPI,MPI_Integer4,root_MPI,MPI_COMM_WORLD,MPI_err)
!#endif
       IF(openmpi .AND. MPI_scheme/=1) CALL MPI_Bcast_(ierr,size1_MPI,root_MPI)
#else
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) ' The ARPACK library is not present!'
        write(out_unitp,*) 'Use Arpack=f and Davidson=t'
        STOP 'ARPACK has been removed'
#endif
!       %-----------------------------------------------%
!       | The real part of the eigenvalue is returned   |
!       | in the first column of the two dimensional    |
!       | array D, and the imaginary part is returned   |
!       | in the second column of D.  The corresponding |
!       | eigenvectors are returned in the first NEV    |
!       | columns of the two dimensional array V if     |
!       | requested.  Otherwise, an orthogonal basis    |
!       | for the invariant subspace corresponding to   |
!       | the eigenvalues in D is returned in V.        |
!       %-----------------------------------------------%

        if ( ierr .ne. 0) then

!         %------------------------------------%
!         | Error condition:                   |
!         | Check the documentation of _neupd  |
!         %------------------------------------%

          write(out_unitp,*)
          write(out_unitp,*) ' Error with _seupd, info = ', ierr
          write(out_unitp,*) ' Check the documentation of _neupd. '
          write(out_unitp,*)

        else

          first  = .true.
          nconv  = iparam(5)
!#if(run_MPI)
!          CALL MPI_Bcast(nconv,size1_MPI,MPI_Integer4,root_MPI,MPI_COMM_WORLD,MPI_err)
!#endif
          IF(openmpi .AND. MPI_scheme/=1) CALL MPI_Bcast_(nconv,size1_MPI,root_MPI)
          do j=1, nconv

!           %---------------------------%
!           | Compute the residual norm |
!           |                           |
!           |   ||  A*x - lambda*x ||   |
!           |                           |
!           | for the NCONV accurately  |
!           | computed eigenvalues and  |
!           | eigenvectors.  (iparam(5) |
!           | indicates how many are    |
!           | accurate to the requested |
!           | tolerance)                |
!           %---------------------------%

            IF(openmpi .AND. MPI_scheme/=1) THEN
              IF(MPI_id==0) THEN
                IF(d(j,2)==zero) THEN
                  if_deq0=.TRUE.
                ELSE
                  if_deq0=.FALSE.
                ENDIF
              ENDIF
              !CALL MPI_Bcast(if_deq0,size1_MPI,MPI_logical,root_MPI,MPI_COMM_WORLD,MPI_err)
              CALL MPI_Bcast_(if_deq0,size1_MPI,root_MPI)
            ENDIF
            
            temp_logi=.FALSE.
            IF(openmpi .AND. MPI_scheme/=1) THEN
              if(if_deq0) temp_logi=.TRUE.
            ELSE
              IF(d(j,2) == zero) temp_logi=.TRUE.
            ENDIF
            IF(temp_logi) THEN

!             %--------------------%
!             | Ritz value is real |
!             %--------------------%

              CALL sub_OpV1_TO_V2_Arpack(v(:,j),ax,psi(j),Hpsi_loc,                    &
                                         para_H,cplxE,para_propa,int(n,4))
#if __ARPACK == 1
              call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
              d(j,3) = dnrm2(n, ax, 1)
              d(j,3) = d(j,3) / abs(d(j,1))
#else
              write(out_unitp,*) 'ERROR in ',name_sub
              write(out_unitp,*) ' The ARPACK library is not present!'
              write(out_unitp,*) 'Use Arpack=f and Davidson=t'
              STOP 'ARPACK has been removed'
#endif


              Ene(j)          = d(j,1)
              psi(j)%CAvOp    = cmplx(d(j,1),ZERO,kind=Rkind)
              psi(j)%IndAvOp  = para_H%n_Op  ! it should be 0
              psi(j)%convAvOp = .TRUE.
            else if (first) then

!             %------------------------%
!             | Ritz value is complex. |
!             | Residual of one Ritz   |
!             | value of the conjugate |
!             | pair is computed.      |
!             %------------------------%

              !call av(nx, v(1,j), ax)
              CALL sub_OpV1_TO_V2_Arpack(v(:,j),ax,psi(j),Hpsi_loc,                    &
                                         para_H,cplxE,para_propa,int(n,4))

#if __ARPACK == 1
              call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
              call daxpy(n, d(j,2), v(1,j+1), 1, ax, 1)
              d(j,3) = dnrm2(n, ax, 1)
#else
              write(out_unitp,*) 'ERROR in ',name_sub
              write(out_unitp,*) ' The ARPACK library is not present!'
              write(out_unitp,*) 'Use Arpack=f and Davidson=t'
              STOP 'ARPACK has been removed'
#endif


              Ene(j)          = d(j,1)
              psi(j)%CAvOp    = cmplx(d(j,1),d(j,2),kind=Rkind)
              psi(j)%IndAvOp  = para_H%n_Op  ! it should be 0
              psi(j)%convAvOp = .TRUE.


              !call av(nx, v(1,j+1), ax)
              CALL sub_OpV1_TO_V2_Arpack(v(:,j+1),ax,psi(j+1),Hpsi_loc,&
                                        para_H,cplxE,para_propa,int(n,4))

              Ene(j+1)          = d(j,1)
              psi(j+1)%CAvOp    = cmplx(d(j,1),-d(j,2),kind=Rkind)
              psi(j+1)%IndAvOp  = para_H%n_Op  ! it should be 0
              psi(j+1)%convAvOp = .TRUE.

#if __ARPACK == 1
              call daxpy(n, -d(j,2), v(1,j), 1, ax, 1)
              call daxpy(n, -d(j,1), v(1,j+1), 1, ax, 1)
              d(j,3) = dlapy2( d(j,3), dnrm2(n, ax, 1) )
              d(j,3) = d(j,3) / dlapy2(d(j,1),d(j,2))
              d(j+1,3) = d(j,3)
              first = .false.
#else
              write(out_unitp,*) 'ERROR in ',name_sub
              write(out_unitp,*) ' The ARPACK library is not present!'
              write(out_unitp,*) 'Use Arpack=f and Davidson=t'
              STOP 'ARPACK has been removed'
#endif
            else
              first = .true.
            end if

            IF (debug) write(out_unitp,*) 'j,cplx ene ?',j,         &
                             psi(j)%CAvOp * get_Conv_au_TO_unit('E','cm-1')

          END DO ! for j=1, nconv

!         %-----------------------------%
!         | Display computed residuals. |
!         %-----------------------------%
#if __ARPACK == 1
          ! bug here for some compiler, diable the output currently
!          IF(MPI_id==0) call dmout(6, nconv, 3, d, maxncv, -6,                     &
!                     'Ritz values (Real,Imag) and relative residuals')
#else
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) ' The ARPACK library is not present!'
          write(out_unitp,*) 'Use Arpack=f and Davidson=t'
          STOP 'ARPACK has been removed'
#endif
        end if ! for ierr .ne. 0

!       %-------------------------------------------%
!       | Print additional convergence information. |
!       %-------------------------------------------%

        if ( info .eq. 1) then
          print *, ' '
          print *, ' Maximum number of iterations reached.'
          print *, ' '
        else if ( info .eq. 3) then
          print *, ' '
          print *, ' No shifts could be applied during implicit',    &
                   ' Arnoldi update, try increasing NCV.'
          print *, ' '
        end if

        print *, ' '
        print *, ' _NDRV1 '
        print *, ' ====== '
        print *, ' '
        print *, ' Size of the matrix is ', n
        print *, ' The number of Ritz values requested is ', nev
        print *, ' The number of Arnoldi vectors generated (NCV) is ', ncv
        print *, ' What portion of the spectrum: ', which
        print *, ' The number of converged Ritz values is ',nconv
        print *, ' The number of Implicit Arnoldi update',             &
                 ' iterations taken is ', iparam(3)
        print *, ' The number of OP*x is ', iparam(9)
        print *, ' The convergence criterion is ', tol
        print *, ' '

      end if ! for info .lt. 0

      nb_diago = nconv

      IF(keep_MPI) CALL trie_psi(psi,Ene,nb_diago)

      if(allocated(v)) deallocate(v)
      deallocate(workl)
      deallocate(workd)
      deallocate(d)
      deallocate(resid)
      deallocate(ax)
      deallocate(select)

!     %---------------------------%
!     | Done with program dndrv1. |
!     %---------------------------%

      CALL dealloc_psi(Hpsi_loc,delete_all=.TRUE.)
      CALL dealloc_psi(psi_loc,delete_all=.TRUE.)


!----------------------------------------------------------
       IF (debug) THEN
        CALL alloc_NParray(Evec,[nb_diago],'Evec',name_sub)
        Evec(:) = real( psi(1:nb_diago)%CAvOp,kind=Rkind)
        ZPE = minval(Evec)
        DO j=1,nb_diago
           write(out_unitp,*) j,Evec(j)*get_Conv_au_TO_unit('E','cm-1'),   &
                             (Evec(j)-ZPE)*get_Conv_au_TO_unit('E','cm-1')
        END DO
        CALL dealloc_NParray(Evec,'Evec',name_sub)
         write(out_unitp,*) 'END ',name_sub
       END IF
!----------------------------------------------------------
      END SUBROUTINE sub_propagation_Arpack
!=======================================================================================
!
!=======================================================================================
      SUBROUTINE sub_OpV1_TO_V2_Arpack(V1,V2,psi1,psi2,                 &
                                       para_H,cplxE,para_propa,n)
      USE mod_system
      USE mod_Op
      USE mod_psi,    ONLY : param_psi,Overlap_psi1_psi2,Set_symab_OF_psiBasisRep
      USE mod_propa
      IMPLICIT NONE

      !----- Operator: Hamiltonian ----------------------------
      TYPE (param_Op)   :: para_H
      integer (kind=4)  :: n
      ! was integer     :: n
      real(kind=Rkind)  :: V1(n),V2(n)

      !----- WP, energy ... -----------------------------------
      TYPE (param_propa)        :: para_propa
      logical                   :: With_Grid,With_Basis,cplx


      TYPE (param_psi)          :: psi1,psi2 ! psi2 = Hpsi
      complex(kind=Rkind)       :: cplxE



!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_OpV1_TO_V2_Arpack'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) ' para_Davidson',para_propa%para_Davidson
        write(out_unitp,*)
        write(out_unitp,*) 'n',n
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

      With_Grid  = para_propa%para_Davidson%With_Grid
      With_Basis = .NOT. para_propa%para_Davidson%With_Grid
      cplx       = para_H%cplx


      IF (With_Grid) THEN
        IF(keep_MPI) psi1%RvecG(:) = V1(:)/sqrt(para_H%BasisnD%wrho(:))
      ELSE
        IF(keep_MPI) psi1%RvecB(:) = V1(:)
      END IF

      CALL sub_OpPsi(psi1,psi2,para_H,With_Grid=With_Grid)

      CALL Set_symab_OF_psiBasisRep(psi2,para_propa%para_Davidson%symab)


      IF (debug) THEN
        CALL Overlap_psi1_psi2(cplxE,psi1,psi2,With_Grid=With_Grid)
        write(out_unitp,*) 'Arpack <psi H psi>:',cplxE
      END IF

      IF (With_Grid) THEN
        IF(keep_MPI) V2(:) = psi2%RvecG(:)*sqrt(para_H%BasisnD%wrho(:))
      ELSE
        IF(keep_MPI) V2(:) = psi2%RvecB(:)
      END IF

!----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'END ',name_sub
       END IF
!----------------------------------------------------------
      END SUBROUTINE sub_OpV1_TO_V2_Arpack
!=======================================================================================

!=======================================================================================
      SUBROUTINE sub_propagation_Arpack_Sym(psi,Ene,nb_diago,max_diago, &
                                          para_H,para_propa)
      USE mod_system
      USE mod_psi,    ONLY : param_psi,alloc_psi,dealloc_psi,param_WP0

      USE mod_Op
      USE mod_propa
      USE mod_MPI_aux
      IMPLICIT NONE

      !----- Operator: Hamiltonian ----------------------------
      TYPE (param_Op)   :: para_H

      !----- WP, energy ... -----------------------------------
      TYPE (param_propa)        :: para_propa
      TYPE (param_WP0)          :: para_WP0

      integer                   :: nb_diago,max_diago
      TYPE (param_psi)          :: psi(max_diago)

      real (kind=Rkind)         :: Ene(max_diago)


      !------ working parameters --------------------------------
      TYPE (param_psi)               :: psi_loc,Hpsi_loc
      real (kind=Rkind)              :: ZPE
      real (kind=Rkind), allocatable :: Vec(:,:),Evec(:)
      complex (kind=Rkind)           :: cplxE

!     %-----------------------------%
!     | Define leading dimensions   |
!     | for all arrays.             |
!     | MAXN:   Maximum dimension   |
!     |         of the A allowed.   |
!     | MAXNEV: Maximum NEV allowed |
!     | MAXNCV: Maximum NCV allowed |
!     %-----------------------------%

      integer (kind=4)       ::   maxn, maxnev, maxncv, ldv

!     %--------------%
!     | Local Arrays |
!     %--------------%

      real(kind=Rkind), allocatable :: v(:,:),workl(:),workd(:),d(:,:),resid(:),ax(:)
      logical, allocatable          :: select(:)
      integer (kind=4)              :: iparam(11), ipntr(11)

!     %---------------%
!     | Local Scalars |
!     %---------------%

      character (len=1)   :: bmat
      character (len=2)   :: which

      integer (kind=4)    :: ido, n, nev, ncv, lworkl, info, ierr, j
      integer (kind=4)    :: nconv, maxitr, mode, ishfts
      logical             :: rvec
      real(kind=Rkind)    :: tol, sigma
      logical             :: if_deq0=.FALSE. !< for MPI 

!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%

      real(kind=Rkind), external ::  dnrm2
      external                       daxpy


      TYPE(param_file)  :: Log_file
      integer           :: iunit
      logical           :: cplx


!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_propagation_Arpack_Sym'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) ' nb_diago',nb_diago
        write(out_unitp,*) ' max_diago',max_diago
        write(out_unitp,*) ' para_Davidson',para_propa%para_Davidson
        write(out_unitp,*)
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

      cplx = para_H%cplx
      IF (cplx) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' cplx is not yet possible in this subroutine!'
        STOP
      END IF

      CALL init_psi(psi_loc,para_H,cplx)
      CALL alloc_psi(psi_loc)
      CALL init_psi(Hpsi_loc,para_H,cplx)
      CALL alloc_psi(Hpsi_loc)

!------ initialization -------------------------------------
      Log_file%name='Arpack.log'
      IF(MPI_id==0) CALL file_open(Log_file,iunit)

      n = para_H%nb_tot
!#if(run_MPI)
!      CALL MPI_Bcast(n,size1_MPI,MPI_Integer4,root_MPI,MPI_COMM_WORLD,MPI_err)
!#endif
      IF(openmpi .AND. MPI_scheme/=1) CALL MPI_Bcast_(n,size1_MPI,root_MPI)

      IF (nb_diago == 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' nb_diago=0 is not possible with ARPACK'
        STOP
      END IF
      nev =  nb_diago
      ncv =  2*nev+10
      
      ncv=MIN(N,ncv)  ! prevent infor=-3 error:
                      ! NCV must be greater than NEV and less than or equal to N.

      maxn   = n
      ldv    = maxn
      maxnev = nev
      maxncv = ncv

      ! allocate: OK
      allocate(v(ldv,maxncv))
      allocate( workl(maxncv*(maxncv+8)) )
      allocate(workd(3*maxn))
      allocate(d(maxncv,2))
      allocate(resid(maxn))
      allocate(ax(maxn))
      allocate(select(maxncv))

      bmat  = 'I'
      which = 'SM' ! smalest magnitude
      !which = 'SA' ! algebraically smalest

!     %--------------------------------------------------%
!     | The work array WORKL is used in DSAUPD as        |
!     | workspace.  Its dimension LWORKL is set as       |
!     | illustrated below.  The parameter TOL determines |
!     | the stopping criterion.  If TOL<=0, machine      |
!     | precision is used.  The variable IDO is used for |
!     | reverse communication and is initially set to 0. |
!     | Setting INFO=0 indicates that a random vector is |
!     | generated in DSAUPD to start the Arnoldi         |
!     | iteration.                                       |
!     %--------------------------------------------------%

      lworkl = ncv*(ncv+8) ! ok
      tol    = ZERO
      ido    = 0

      info   = 0
      info   = 1
      IF (info /= 0) THEN
        IF(keep_MPI) CALL ReadWP0_Arpack(psi_loc,nb_diago,max_diago,                   &
                                          para_propa%para_Davidson,para_H%cplx)
        IF (para_propa%para_Davidson%With_Grid) THEN
          IF(keep_MPI) resid(:) = psi_loc%RvecG(:)
        ELSE
          IF(keep_MPI) resid(:) = psi_loc%RvecB(:)
        END IF
      END IF


      DO j=1,max_diago
        CALL init_psi(psi(j),para_H,cplx)
      END DO
      DO j=1,nb_diago
        CALL alloc_psi(psi(j))
      END DO


!     %---------------------------------------------------%
!     | This program uses exact shifts with respect to    |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of DSAUPD is used     |
!     | (IPARAM(7) = 1).  All these options may be        |
!     | changed by the user. For details, see the         |
!     | documentation in DSAUPD.                          |
!     %---------------------------------------------------%

      ishfts = 1
      maxitr = min(max_diago,para_propa%para_Davidson%max_it)
      mode   = 1

      iparam(1) = ishfts
      iparam(3) = maxitr
      iparam(7) = mode

!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%

      DO
!       %---------------------------------------------%
!       | Repeatedly call the routine DSAUPD and take |
!       | actions indicated by parameter IDO until    |
!       | either convergence is indicated or maxitr   |
!       | has been exceeded.                          |
!       %---------------------------------------------%
#if __ARPACK == 1
        IF(keep_MPI) call dsaupd(ido, bmat, n, which, nev, tol, resid,             &
                                  ncv, v, ldv, iparam, ipntr, workd, workl,        &
                                  lworkl, info )
!#if(run_MPI)
!        CALL MPI_Bcast(info,size1_MPI, MPI_Integer4,root_MPI,MPI_COMM_WORLD,MPI_err)
!        CALL MPI_Bcast(ido, size1_MPI, MPI_Integer4,root_MPI,MPI_COMM_WORLD,MPI_err)
!        CALL MPI_Bcast(ipntr,INT(11,Ikind),MPI_Integer4,root_MPI,MPI_COMM_WORLD,MPI_err)
!#endif
        IF(openmpi .AND. MPI_scheme/=1) THEN
          CALL MPI_Bcast_(info,size1_MPI,root_MPI)
          CALL MPI_Bcast_(ido ,size1_MPI,root_MPI)
          CALL MPI_Bcast_(ipntr,INT(11,kind=MPI_INTEGER_KIND),root_MPI)
        ENDIF

#else
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) ' The ARPACK library is not present!'
          write(out_unitp,*) 'Use Arpack=f and Davidson=t'
          STOP 'ARPACK has been removed'
#endif

         IF (abs(ido) /= 1) EXIT

!        %--------------------------------------%
!        | Perform matrix vector multiplication |
!        |              y <--- OP*x             |
!        | The user should supply his/her own   |
!        | matrix vector multiplication routine |
!        | here that takes workd(ipntr(1)) as   |
!        | the input, and return the result to  |
!        | workd(ipntr(2)).                     |
!        %--------------------------------------%

         CALL sub_OpV1_TO_V2_Arpack(workd(ipntr(1):ipntr(1)-1+n),       &
                                    workd(ipntr(2):ipntr(2)-1+n),       &
                                    psi_loc,Hpsi_loc,para_H,cplxE,para_propa,int(n,kind=4))

         !psi_loc%RvecB(:) = workd(ipntr(1):ipntr(1)-1+n)
         !CALL sub_OpPsi(psi_loc,Hpsi_loc,para_H)
         !workd(ipntr(2):ipntr(2)-1+n) = Hpsi_loc%RvecB(:)
         
        !IF(keep_MPI) psi_loc%RvecB(:) = workd(ipntr(1):ipntr(1)-1+n)
        !CALL sub_OpPsi(psi_loc,Hpsi_loc,para_H)
        !IF(keep_MPI) THEN
        !  workd(ipntr(2):ipntr(2)-1+n) = Hpsi_loc%RvecB(:)
        
        IF(MPI_id==0) THEN  
          write(iunit,*) 'Arpack <psi H psi>:',                          &
           dot_product(workd(ipntr(1):ipntr(1)-1+n),workd(ipntr(2):ipntr(2)-1+n))
          CALL flush_perso(iunit)
        ENDIF
      END DO

      IF(MPI_id==0) THEN
        write(iunit,*) 'End Arpack ' ; CALL flush_perso(iunit)
        CALL file_close(Log_file)
      ENDIF

!----------------------------------------------------------


!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%

      if ( info .lt. 0 ) then
!
!       %--------------------------%
!       | Error message. Check the |
!       | documentation in DSAUPD. |
!       %--------------------------%

         write(out_unitp,*)
         write(out_unitp,*) ' Error with _saupd, info = ', info,' from ', MPI_id
         write(out_unitp,*) ' Check documentation in _saupd '
         write(out_unitp,*) ' '

      else

!       %-------------------------------------------%
!       | No fatal errors occurred.                 |
!       | Post-Process using DSEUPD.                |
!       |                                           |
!       | Computed eigenvalues may be extracted.    |
!       |                                           |
!       | Eigenvectors may also be computed now if  |
!       | desired.  (indicated by rvec = .true.)    |
!       %-------------------------------------------%

        rvec = .true.

#if __ARPACK == 1
        IF(keep_MPI) call dseupd(rvec, 'All', select, d, v, ldv, sigma,                &
                                  bmat, n, which, nev, tol, resid, ncv, v, ldv,        &
                                  iparam, ipntr, workd, workl, lworkl, ierr )
!#if(run_MPI)
!        CALL MPI_Bcast(ierr,size1_MPI,MPI_Integer4,root_MPI,MPI_COMM_WORLD,MPI_err)
!#endif
        IF(openmpi .AND. MPI_scheme/=1) CALL MPI_Bcast_(ierr,size1_MPI,root_MPI)

#endif
!       %----------------------------------------------%
!       | Eigenvalues are returned in the first column |
!       | of the two dimensional array D and the       |
!       | corresponding eigenvectors are returned in   |
!       | the first NEV columns of the two dimensional |
!       | array V if requested.  Otherwise, an         |
!       | orthogonal basis for the invariant subspace  |
!       | corresponding to the eigenvalues in D is     |
!       | returned in V.                               |
!       %----------------------------------------------%

        if ( ierr /= 0) then

!         %------------------------------------%
!         | Error condition:                   |
!         | Check the documentation of DSEUPD. |
!         %------------------------------------%

          write(out_unitp,*)
          write(out_unitp,*) ' Error with _seupd, info = ', ierr
          write(out_unitp,*) ' Check the documentation of _seupd. '
          write(out_unitp,*)
          STOP

        else

          nconv =  iparam(5)
!#if(run_MPI)
!          CALL MPI_Bcast(nconv,size1_MPI,MPI_Integer4,root_MPI,MPI_COMM_WORLD,MPI_err)
!#endif
          IF(openmpi .AND. MPI_scheme/=1) CALL MPI_Bcast_(nconv,size1_MPI,root_MPI)

          DO j=1, nconv

!           %---------------------------%
!           | Compute the residual norm |
!           |                           |
!           |   ||  A*x - lambda*x ||   |
!           |                           |
!           | for the NCONV accurately  |
!           | computed eigenvalues and  |
!           | eigenvectors.  (iparam(5) |
!           | indicates how many are    |
!           | accurate to the requested |
!           | tolerance)                |
!           %---------------------------%

            !call av(nx, v(1,j), ax)
            CALL sub_OpV1_TO_V2_Arpack(v(:,j),ax,psi(j),Hpsi_loc,&
                                       para_H,cplxE,para_propa,int(n,kind=4))

            !psi(j)%RvecB(:) = v(:,j)
            !CALL sub_OpPsi(psi(j),Hpsi_loc,para_H)
            !ax(:) = Hpsi_loc%RvecB(:)
            
            !IF(keep_MPI) psi(j)%RvecB(:) = v(:,j)
            !CALL sub_OpPsi(psi(j),Hpsi_loc,para_H)
            !IF(keep_MPI) ax(:) = Hpsi_loc%RvecB(:)

#if __ARPACK == 1
            IF(keep_MPI) THEN
              call daxpy(n, -d(j,1), v(:,j), 1, ax, 1)
              d(j,2) = dnrm2(n, ax, 1)
              d(j,2) = d(j,2) / abs(d(j,1))
            END IF
#endif
            IF(keep_MPI) THEN
              IF (debug) write(out_unitp,*) 'j,ene ?',j,              &
                             d(j,1) * get_Conv_au_TO_unit('E','cm-1')

              Ene(j)          = d(j,1)
              psi(j)%CAvOp    = Ene(j)
              psi(j)%IndAvOp  = para_H%n_Op  ! it should be 0
              psi(j)%convAvOp = .TRUE.
            ENDIF ! for keep_MPI
          END DO

!         %-------------------------------%
!         | Display computed residuals    |
!         %-------------------------------%
#if __ARPACK == 1
          ! bug here for some compiler, diable the output currently
!          IF(keep_MPI) call dmout(6, nconv, 2, d, maxncv, -6,                         &
!                        'Ritz values and relative residuals')
#endif
        end if


!       %------------------------------------------%
!       | Print additional convergence information |
!       %------------------------------------------%

        if ( info .eq. 1) then
          print *, ' '
          print *, ' Maximum number of iterations reached.'
          print *, ' '
        else if ( info .eq. 3) then
          print *, ' '
          print *, ' No shifts could be applied during implicit',     &
                   ' Arnoldi update, try increasing NCV.'
          print *, ' '
        end if

        print *, ' '
        print *, ' _SDRV1 '
        print *, ' ====== '
        print *, ' '
        print *, ' Size of the matrix is ', n
        print *, ' The number of Ritz values requested is ', nev
        print *, ' The number of Arnoldi vectors generated',           &
                ' (NCV) is ', ncv
        print *, ' What portion of the spectrum: ', which
        print *, ' The number of converged Ritz values is ',           &
                   nconv
        print *, ' The number of Implicit Arnoldi update',             &
                 ' iterations taken is ', iparam(3)
        print *, ' The number of OP*x is ', iparam(9)
        print *, ' The convergence criterion is ', tol
        print *, ' '

      end if

      nb_diago = nconv


      CALL dealloc_psi(Hpsi_loc,delete_all=.TRUE.)
      CALL dealloc_psi(psi_loc,delete_all=.TRUE.)
      deallocate(v)
      deallocate(workl)
      deallocate(workd)
      deallocate(d)
      deallocate(resid)
      deallocate(ax)
      deallocate(select)

!----------------------------------------------------------
       IF (debug) THEN
        CALL alloc_NParray(Evec,[nb_diago],'Evec',name_sub)
        Evec(:) = real( psi(1:nb_diago)%CAvOp,kind=Rkind)
        ZPE = minval(Evec)
        DO j=1,nb_diago
           write(out_unitp,*) j,Evec(j)*get_Conv_au_TO_unit('E','cm-1'),   &
                             (Evec(j)-ZPE)*get_Conv_au_TO_unit('E','cm-1')
        END DO
        CALL dealloc_NParray(Evec,'Evec',name_sub)
         write(out_unitp,*) 'END ',name_sub
       END IF
!----------------------------------------------------------

      END SUBROUTINE sub_propagation_Arpack_Sym
!=======================================================================================

 SUBROUTINE ReadWP0_Arpack(psi,nb_diago,max_diago,para_Davidson,cplx)
 USE mod_system
 USE mod_psi,    ONLY : param_psi,norm2_psi,renorm_psi,dealloc_psi,     &
                        Set_symab_OF_psiBasisRep,copy_psi2TOpsi1,       &
                        sub_PsiBasisRep_TO_GridRep,alloc_psi,           &
                        param_WP0,sub_read_psi0,set_random_psi, ecri_init_psi

 USE mod_propa,  ONLY : param_Davidson
 IMPLICIT NONE


 TYPE (param_Davidson) :: para_Davidson
 logical               :: cplx

 !----- WP, energy ... -----------------------------------
 TYPE (param_WP0)          :: para_WP0

 integer                   :: nb_diago,max_diago
 TYPE (param_psi)          :: psi
 TYPE (param_psi)          :: psi0(max_diago)

 integer :: i

 !----- for debuging --------------------------------------------------
 integer :: err_mem,memory
 character (len=*), parameter :: name_sub='ReadWP0_Arpack'
 logical, parameter :: debug=.FALSE.
 !logical, parameter :: debug=.TRUE.
 !-----------------------------------------------------------

 IF (debug) THEN
   write(out_unitp,*) 'BEGINNING ',name_sub
   write(out_unitp,*) ' nb_diago',nb_diago
   write(out_unitp,*) ' max_diago',max_diago

   write(out_unitp,*) ' para_Davidson',para_Davidson
   write(out_unitp,*)
   CALL flush_perso(out_unitp)
 END IF
 !-----------------------------------------------------------

 !------ read guess vectors ---------------------------------
 IF (nb_diago > 0 .OR. para_Davidson%read_WP) THEN

   DO i=1,max_diago
     CALL copy_psi2TOpsi1(psi0(i),psi,alloc=.FALSE.)
   END DO

   para_WP0%nb_WP0              = nb_diago
   para_WP0%read_file           = para_Davidson%read_WP
   IF (para_Davidson%nb_readWP_OF_List > 0) THEN
     para_WP0%read_listWP0        = .FALSE.
   ELSE
     para_WP0%read_listWP0        = para_Davidson%read_listWP
   END IF
   para_WP0%WP0cplx             = cplx
   para_WP0%file_WP0%name       = para_Davidson%name_file_readWP
   para_WP0%file_WP0%formatted  = para_Davidson%formatted_file_readWP

   IF(keep_MPI) CALL sub_read_psi0(psi0,para_WP0,max_diago,                            &
                                    symab=para_Davidson%symab,ortho=.TRUE.)

   nb_diago = para_WP0%nb_WP0
   para_Davidson%nb_WP0 = para_WP0%nb_WP0

   IF (nb_diago < 1) THEN
     write(out_unitp,*) ' ERROR while reading the vector(s)'
     write(out_unitp,*) ' ... the number is 0'
     write(out_unitp,*) ' Probably, in the namelist (davidson)...'
     write(out_unitp,*) '   ... you select a wrong symmetry',para_Davidson%symab
     STOP
   END IF

   psi = ZERO
   DO i=1,nb_diago
     psi = psi + psi0(i)
     CALL norm2_psi(psi0(i))
     IF (debug) write(out_unitp,*) '   norm^2 of psi0(i)',i,psi0(i)%norm2
     IF ( abs(psi0(i)%norm2-ONE) > ONETENTH**8) THEN
       write(out_unitp,*) ' ERROR while reading the vector(s)'
       write(out_unitp,*) ' ... the norm^2 of psi0(i) is /= 1',i,psi0(i)%norm2
       STOP
     END IF
   END DO
   CALL Set_symab_OF_psiBasisRep(psi,para_Davidson%symab)
   CALL renorm_psi(psi,BasisRep=.TRUE.)


 ELSE
   para_Davidson%nb_WP0 = 0
   nb_diago = 1
   CALL alloc_psi(psi)
   CALL Set_Random_psi(psi)
   CALL Set_symab_OF_psiBasisRep(psi,para_Davidson%symab)
   CALL renorm_psi(psi,BasisRep=.TRUE.)
 END IF

 DO i=1,max_diago
   CALL dealloc_psi(psi0(i),delete_all=.TRUE.)
 END DO

 !----------------------------------------------------------
 IF (debug) THEN
   write(out_unitp,*) 'psi, nb_diago',nb_diago
   CALL ecri_init_psi(psi)
   write(out_unitp,*) 'END ',name_sub
 END IF
 !----------------------------------------------------------

 END SUBROUTINE ReadWP0_Arpack

END MODULE mod_Arpack
