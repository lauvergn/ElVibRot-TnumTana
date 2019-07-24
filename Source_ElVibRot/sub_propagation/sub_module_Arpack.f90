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
MODULE mod_Arpack
USE mod_Constant
IMPLICIT NONE

PRIVATE
PUBLIC :: sub_propagation_Arpack,sub_propagation_Arpack_Sym

CONTAINS

      SUBROUTINE sub_propagation_Arpack(psi,Ene,nb_diago,max_diago,     &
                                        para_H,para_propa)
      USE mod_system
      USE mod_Op

      USE mod_psi_set_alloc
      USE mod_psi_SimpleOp
      USE mod_ana_psi
      USE mod_psi_B_TO_G
      USE mod_psi_Op
      USE mod_param_WP0
      USE mod_propa
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
      real (kind=Rkind), allocatable :: Vec(:,:),Evec(:)

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

      real(kind=Rkind)    ::  tol, sigmar, sigmai

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
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
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
      DO j=1,max_diago
        CALL init_psi(psi(j),para_H,cplx)
      END DO
      DO j=1,nb_diago
        CALL alloc_psi(psi(j),BasisRep=With_Basis,GridRep=With_Grid)
      END DO





!------ initialization -------------------------------------
      Log_file%name='Arpack.log'
      CALL file_open(Log_file,iunit)



      IF (With_Grid) THEN
        n = int(para_H%nb_qaie,kind=4)
      ELSE
        n = int(para_H%nb_tot,kind=4)
      END IF
      IF (nb_diago == 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' nb_diago=0 is not possible with ARPACK'
        STOP
      END IF

      nev =  int(nb_diago,kind=4)
      ncv =  2*nev+10

      maxn   = n
      ldv    = maxn
      maxnev = nev
      maxncv = ncv

      allocate(select(maxncv))


      allocate(ax(maxn))
      allocate(d(maxncv,3))
      allocate(resid(maxn))

      allocate(v(ldv,maxncv))
      allocate(workd(3*maxn))
      allocate(workev(3*maxncv))
      allocate( workl(3*maxncv*maxncv+6*maxncv) )

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
      info   = 0

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

      DO
!        %---------------------------------------------%
!        | Repeatedly call the routine DNAUPD and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
#if __ARPACK == 1
         call dnaupd ( ido, bmat, n, which, nev, tol, resid,            &
              ncv, v, ldv, iparam, ipntr, workd, workl, lworkl,         &
              info )
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
                              psi_loc,Hpsi_loc,para_H,cplxE,para_propa,int(n))

         write(iunit,*) 'Arpack <psi H psi>:',                          &
          dot_product(workd(ipntr(1):ipntr(1)-1+n),workd(ipntr(2):ipntr(2)-1+n))
         CALL flush_perso(iunit)

      END DO

      write(iunit,*) 'End Arpack ' ; CALL flush_perso(iunit)
      CALL file_close(Log_file)


!----------------------------------------------------------


!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%

      if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | Error message. Check the |
!        | documentation in DSAUPD. |
!        %--------------------------%

         write(out_unitp,*)
         write(out_unitp,*) ' Error with _saupd, info = ', info
         write(out_unitp,*) ' Check documentation in _naupd '
         write(out_unitp,*) ' '

      else

!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DSEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        %-------------------------------------------%

         rvec = .true.

#if __ARPACK == 1
         call dneupd ( rvec, 'A', select, d, d(1,2), v, ldv,            &
              sigmar, sigmai, workev, bmat, n, which, nev, tol,         &
              resid, ncv, v, ldv, iparam, ipntr, workd, workl,          &
              lworkl, ierr )
#else
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) ' The ARPACK library is not present!'
          write(out_unitp,*) 'Use Arpack=f and Davidson=t'
          STOP 'ARPACK has been removed'
#endif
!        %-----------------------------------------------%
!        | The real part of the eigenvalue is returned   |
!        | in the first column of the two dimensional    |
!        | array D, and the imaginary part is returned   |
!        | in the second column of D.  The corresponding |
!        | eigenvectors are returned in the first NEV    |
!        | columns of the two dimensional array V if     |
!        | requested.  Otherwise, an orthogonal basis    |
!        | for the invariant subspace corresponding to   |
!        | the eigenvalues in D is returned in V.        |
!        %-----------------------------------------------%

         if ( ierr .ne. 0) then

!            %------------------------------------%
!            | Error condition:                   |
!            | Check the documentation of _neupd  |
!            %------------------------------------%

             write(out_unitp,*)
             write(out_unitp,*) ' Error with _seupd, info = ', ierr
             write(out_unitp,*) ' Check the documentation of _neupd. '
             write(out_unitp,*)

         else

             first  = .true.
             nconv  = iparam(5)
             do j=1, nconv

!               %---------------------------%
!               | Compute the residual norm |
!               |                           |
!               |   ||  A*x - lambda*x ||   |
!               |                           |
!               | for the NCONV accurately  |
!               | computed eigenvalues and  |
!               | eigenvectors.  (iparam(5) |
!               | indicates how many are    |
!               | accurate to the requested |
!               | tolerance)                |
!               %---------------------------%

                if (d(j,2) == zero)  then

!                  %--------------------%
!                  | Ritz value is real |
!                  %--------------------%

                   CALL sub_OpV1_TO_V2_Arpack(v(:,j),ax,psi(j),Hpsi_loc,&
                                              para_H,cplxE,para_propa,int(n))

                   call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                   d(j,3) = dnrm2(n, ax, 1)
                   d(j,3) = d(j,3) / abs(d(j,1))

                   Ene(j)          = d(j,1)
                   psi(j)%CAvOp    = cmplx(d(j,1),ZERO,kind=Rkind)
                   psi(j)%IndAvOp  = para_H%n_Op  ! it should be 0
                   psi(j)%convAvOp = .TRUE.
                else if (first) then

!                  %------------------------%
!                  | Ritz value is complex. |
!                  | Residual of one Ritz   |
!                  | value of the conjugate |
!                  | pair is computed.      |
!                  %------------------------%

                   !call av(nx, v(1,j), ax)
                   CALL sub_OpV1_TO_V2_Arpack(v(:,j),ax,psi(j),Hpsi_loc,&
                                             para_H,cplxE,para_propa,int(n))

#if __ARPACK == 1
                   call daxpy(n, -d(j,1), v(1,j), 1, ax, 1)
                   call daxpy(n, d(j,2), v(1,j+1), 1, ax, 1)
#else
             write(out_unitp,*) 'ERROR in ',name_sub
             write(out_unitp,*) ' The ARPACK library is not present!'
             write(out_unitp,*) 'Use Arpack=f and Davidson=t'
             STOP 'ARPACK has been removed'
#endif

                   d(j,3) = dnrm2(n, ax, 1)

                   Ene(j)          = d(j,1)
                   psi(j)%CAvOp    = cmplx(d(j,1),d(j,2),kind=Rkind)
                   psi(j)%IndAvOp  = para_H%n_Op  ! it should be 0
                   psi(j)%convAvOp = .TRUE.


                   !call av(nx, v(1,j+1), ax)
                   CALL sub_OpV1_TO_V2_Arpack(v(:,j+1),ax,psi(j+1),Hpsi_loc,&
                                             para_H,cplxE,para_propa,int(n))

                   Ene(j+1)          = d(j,1)
                   psi(j+1)%CAvOp    = cmplx(d(j,1),-d(j,2),kind=Rkind)
                   psi(j+1)%IndAvOp  = para_H%n_Op  ! it should be 0
                   psi(j+1)%convAvOp = .TRUE.

#if __ARPACK == 1
                   call daxpy(n, -d(j,2), v(1,j), 1, ax, 1)
                   call daxpy(n, -d(j,1), v(1,j+1), 1, ax, 1)
#else
             write(out_unitp,*) 'ERROR in ',name_sub
             write(out_unitp,*) ' The ARPACK library is not present!'
             write(out_unitp,*) 'Use Arpack=f and Davidson=t'
             STOP 'ARPACK has been removed'
#endif
                   d(j,3) = dlapy2( d(j,3), dnrm2(n, ax, 1) )
                   d(j,3) = d(j,3) / dlapy2(d(j,1),d(j,2))
                   d(j+1,3) = d(j,3)
                   first = .false.
                else
                   first = .true.
                end if

                IF (debug) write(out_unitp,*) 'j,cplx ene ?',j,         &
                             psi(j)%CAvOp * get_Conv_au_TO_unit('E','cm-1')

            END DO

!            %-----------------------------%
!            | Display computed residuals. |
!            %-----------------------------%
#if __ARPACK == 1
             call dmout(6, nconv, 3, d, maxncv, -6,                     &
                       'Ritz values (Real,Imag) and relative residuals')
#else
             write(out_unitp,*) 'ERROR in ',name_sub
             write(out_unitp,*) ' The ARPACK library is not present!'
             write(out_unitp,*) 'Use Arpack=f and Davidson=t'
             STOP 'ARPACK has been removed'
#endif
          end if

!        %-------------------------------------------%
!        | Print additional convergence information. |
!        %-------------------------------------------%

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

      end if

      nb_diago = nconv

      CALL trie_psi(psi,Ene,nb_diago)

      deallocate(v)
      deallocate(workl)
      deallocate(workd)
      deallocate(d)
      deallocate(resid)
      deallocate(ax)
      deallocate(select)

!     %---------------------------%
!     | Done with program dndrv1. |
!     %---------------------------%



      IF (debug) THEN
        CALL sub_build_MatOp(psi,nb_diago,para_H,.TRUE.,.FALSE.)
        CALL alloc_NParray(vec,(/ nb_diago,nb_diago /),'vec',name_sub)
        CALL alloc_NParray(Evec,(/ nb_diago /),'Evec',name_sub)
        CALL diagonalization(para_H%Rmat,Evec,Vec,nb_diago,2,1,.TRUE.)
        ZPE = minval(Evec)
        DO j=1,nb_diago
           write(out_unitp,*) j,Evec(j)*get_Conv_au_TO_unit('E','cm-1'),   &
                             (Evec(j)-ZPE)*get_Conv_au_TO_unit('E','cm-1')
        END DO
        CALL dealloc_NParray(vec,'vec',name_sub)
        CALL dealloc_NParray(Evec,'Evec',name_sub)
      END IF

      CALL dealloc_psi(Hpsi_loc,delete_all=.TRUE.)
      CALL dealloc_psi(psi_loc,delete_all=.TRUE.)


!----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'END ',name_sub
       END IF
!----------------------------------------------------------
      END SUBROUTINE sub_propagation_Arpack

      SUBROUTINE sub_OpV1_TO_V2_Arpack(V1,V2,psi1,psi2,                 &
                                       para_H,cplxE,para_propa,n)
      USE mod_system
      USE mod_Op

      USE mod_psi_set_alloc
      USE mod_psi_SimpleOp
      USE mod_ana_psi
      USE mod_psi_B_TO_G
      USE mod_psi_Op,         ONLY : Overlap_psi1_psi2
      USE mod_param_WP0
      USE mod_propa
      IMPLICIT NONE

      !----- Operator: Hamiltonian ----------------------------
      TYPE (param_Op)   :: para_H
      integer           :: n
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
        psi1%RvecG(:) = V1(:)/sqrt(para_H%BasisnD%wrho(:))
      ELSE
        psi1%RvecB(:) = V1(:)
      END IF

      CALL sub_OpPsi(psi1,psi2,para_H,With_Grid=With_Grid)

      IF (debug) THEN
        CALL Overlap_psi1_psi2(cplxE,psi1,psi2,With_Grid=With_Grid)
        write(out_unitp,*) 'Arpack <psi H psi>:',cplxE
      END IF

      IF (With_Grid) THEN
        V2(:) = psi2%RvecG(:)*sqrt(para_H%BasisnD%wrho(:))
      ELSE
        V2(:) = psi2%RvecB(:)
      END IF

!----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'END ',name_sub
       END IF
!----------------------------------------------------------
      END SUBROUTINE sub_OpV1_TO_V2_Arpack

      SUBROUTINE sub_propagation_Arpack_Sym(psi,Ene,nb_diago,max_diago, &
                                          para_H,para_propa)
      USE mod_system
      USE mod_Op

      USE mod_psi_set_alloc
      USE mod_psi_SimpleOp
      USE mod_ana_psi
      USE mod_psi_B_TO_G
      USE mod_param_WP0
      USE mod_propa
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
      TYPE (param_psi)          :: psi_loc,Hpsi_loc
      real (kind=Rkind) :: ZPE
      real (kind=Rkind), allocatable :: Vec(:,:),Evec(:)

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

!     %-----------------------------%
!     | BLAS & LAPACK routines used |
!     %-----------------------------%

      real(kind=Rkind) ::  dnrm2
      external         dnrm2, daxpy


      TYPE(param_file)  :: Log_file
      integer           :: iunit
      logical           :: cplx


!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_propagation_Arpack_Sym'
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
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
      DO j=1,max_diago
        CALL init_psi(psi(j),para_H,cplx)
      END DO
      DO j=1,nb_diago
        CALL alloc_psi(psi(j))
      END DO





!------ initialization -------------------------------------
      Log_file%name='Arpack.log'
      CALL file_open(Log_file,iunit)



      n = para_H%nb_tot
      IF (nb_diago == 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' nb_diago=0 is not possible with ARPACK'
        STOP
      END IF
      nev =  nb_diago
      ncv =  nev+10

      maxn   = n
      ldv    = maxn
      maxnev = nev
      maxncv = ncv


      allocate(v(ldv,maxncv))
      allocate( workl(maxncv*(maxncv+8)) )
      allocate(workd(3*maxn))
      allocate(d(maxncv,2))
      allocate(resid(maxn))
      allocate(ax(maxn))
      allocate(select(maxncv))


      bmat  = 'I'
      which = 'SM'

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

      lworkl = ncv*(ncv+8)
      tol    = ZERO
      ido    = 0
      info   = 0

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
!        %---------------------------------------------%
!        | Repeatedly call the routine DSAUPD and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
#if __ARPACK == 1
         call dsaupd ( ido, bmat, n, which, nev, tol, resid,            &
                       ncv, v, ldv, iparam, ipntr, workd, workl,        &
                       lworkl, info )
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

         psi_loc%RvecB(:) = workd(ipntr(1):ipntr(1)-1+n)
         CALL sub_OpPsi(psi_loc,Hpsi_loc,para_H)
         workd(ipntr(2):ipntr(2)-1+n) = Hpsi_loc%RvecB(:)

         write(iunit,*) 'Arpack <psi H psi>:',                          &
          dot_product(workd(ipntr(1):ipntr(1)-1+n),workd(ipntr(2):ipntr(2)-1+n))
         CALL flush_perso(iunit)

      END DO


      write(iunit,*) 'End Arpack ' ; CALL flush_perso(iunit)
      CALL file_close(Log_file)


!----------------------------------------------------------


!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%

      if ( info .lt. 0 ) then
!
!        %--------------------------%
!        | Error message. Check the |
!        | documentation in DSAUPD. |
!        %--------------------------%

         write(out_unitp,*)
         write(out_unitp,*) ' Error with _saupd, info = ', info
         write(out_unitp,*) ' Check documentation in _saupd '
         write(out_unitp,*) ' '

      else

!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using DSEUPD.                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        %-------------------------------------------%

         rvec = .true.

#if __ARPACK == 1
         call dseupd ( rvec, 'All', select, d, v, ldv, sigma,           &
              bmat, n, which, nev, tol, resid, ncv, v, ldv,             &
              iparam, ipntr, workd, workl, lworkl, ierr )
#endif
!        %----------------------------------------------%
!        | Eigenvalues are returned in the first column |
!        | of the two dimensional array D and the       |
!        | corresponding eigenvectors are returned in   |
!        | the first NEV columns of the two dimensional |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%

         if ( ierr .ne. 0) then

!            %------------------------------------%
!            | Error condition:                   |
!            | Check the documentation of DSEUPD. |
!            %------------------------------------%

             write(out_unitp,*)
             write(out_unitp,*) ' Error with _seupd, info = ', ierr
             write(out_unitp,*) ' Check the documentation of _seupd. '
             write(out_unitp,*)
             STOP

         else

             nconv =  iparam(5)
             DO j=1, nconv

!               %---------------------------%
!               | Compute the residual norm |
!               |                           |
!               |   ||  A*x - lambda*x ||   |
!               |                           |
!               | for the NCONV accurately  |
!               | computed eigenvalues and  |
!               | eigenvectors.  (iparam(5) |
!               | indicates how many are    |
!               | accurate to the requested |
!               | tolerance)                |
!               %---------------------------%

                psi(j)%RvecB(:) = v(:,j)
                CALL sub_OpPsi(psi(j),Hpsi_loc,para_H)
                ax(:) = Hpsi_loc%RvecB(:)

#if __ARPACK == 1
                call daxpy(n, -d(j,1), v(:,j), 1, ax, 1)
#endif
                d(j,2) = dnrm2(n, ax, 1)
                d(j,2) = d(j,2) / abs(d(j,1))

                IF (debug) write(out_unitp,*) 'j,ene ?',j,              &
                             d(j,1) * get_Conv_au_TO_unit('E','cm-1')

                Ene(j)          = d(j,1)
                psi(j)%CAvOp    = Ene(j)
                psi(j)%IndAvOp  = para_H%n_Op  ! it should be 0
                psi(j)%convAvOp = .TRUE.

             END DO

!            %-------------------------------%
!            | Display computed residuals    |
!            %-------------------------------%
#if __ARPACK == 1
             call dmout(6, nconv, 2, d, maxncv, -6,                     &
                        'Ritz values and relative residuals')
#endif
         end if


!        %------------------------------------------%
!        | Print additional convergence information |
!        %------------------------------------------%

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


      IF (debug) THEN
        CALL sub_build_MatOp(psi,nb_diago,para_H,.TRUE.,.FALSE.)
        CALL alloc_NParray(vec,(/ nb_diago,nb_diago /),'vec',name_sub)
        CALL alloc_NParray(Evec,(/ nb_diago /),'Evec',name_sub)
        CALL diagonalization(para_H%Rmat,Evec,Vec,nb_diago,2,1,.TRUE.)
        ZPE = minval(Evec)
        DO j=1,nb_diago
           write(out_unitp,*) j,Evec(j)*get_Conv_au_TO_unit('E','cm-1'),&
                          (Evec(j)-ZPE)*get_Conv_au_TO_unit('E','cm-1')
        END DO
        CALL dealloc_NParray(vec,'vec',name_sub)
        CALL dealloc_NParray(Evec,'Evec',name_sub)
      END IF


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
         write(out_unitp,*) 'END ',name_sub
       END IF
!----------------------------------------------------------
      END SUBROUTINE sub_propagation_Arpack_Sym

END MODULE mod_Arpack
