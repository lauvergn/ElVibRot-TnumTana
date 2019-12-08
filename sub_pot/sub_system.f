c
C================================================================
C    calc_Op : calculation of the potential and dipolar matrices
c    mat_V(nb_be,nb_be) and mat_dip(nb_be,nb_be,3)
c    nb_be : nb of elctronic surface
c    Q are the coordinates in active order or syl order
c    dipolar calculation if calc_dip = T
C================================================================
      SUBROUTINE calcN_op(mat_V,mat_imV,mat_ScalOp,nb_be,nb_ScalOp,
     *                   Qdyn,nb_var,mole,
     *                   calc_ScalOp,pot_cplx)

      USE mod_system
      USE mod_Tnum
      IMPLICIT NONE

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

      integer           :: nb_be,nb_ScalOp,nb_var
      logical           :: calc_ScalOp,pot_cplx
      real (kind=Rkind) :: mat_V(nb_be,nb_be),mat_imV(nb_be,nb_be)
      real (kind=Rkind) :: mat_ScalOp(nb_be,nb_be,nb_ScalOp)
      real (kind=Rkind) :: Qdyn(nb_var)
      real (kind=Rkind) :: Qsym(nb_var)

      real (kind=Rkind) :: pot0,im_pot0
      real (kind=Rkind) :: dip(3)
      integer  :: it


      !write(6,*) 'Qdyn',Qdyn
      IF (nb_be == 1 ) THEN
        mat_V(1,1) = pot0(Qdyn)
        IF (pot_cplx) mat_imV(1,1) = im_pot0(Qsym)
        IF (calc_ScalOp) THEN
          CALL sub_dipole(dip,Qsym)
          mat_ScalOp(1,1,1) = dip(1)
          mat_ScalOp(1,1,2) = dip(2)
          mat_ScalOp(1,1,3) = dip(3)
        END IF
      END IF

      RETURN
      END
C================================================================
C    fonction im_pot0(x) imaginary part of pot0
C================================================================
       real*8 FUNCTION im_pot0(Qsym0)


       real*8 Qsym0(1)

       im_pot0 = 0.0

       END
C================================================================
C    fonction pot_rest(x) rest of the DL : pot0 + v2 + pot_rest
C================================================================
       real*8 FUNCTION pot_rest(Qact,Delta_Qact,nb_inact2n)

       real*8 Qact(1)
       integer nb_inact2n
       real*8 Delta_Qact(nb_inact2n)

       pot_rest = 0.0

      END
C================================================================
C    subroutine calculant le gradient
C================================================================
       SUBROUTINE d0d1d2_g(d0g,d1g,d2g,Qsym0,mole,deriv,num,step)

      USE mod_Tnum
      implicit none

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

      real (kind=8) :: d0g(mole%nb_inact2n)
      real (kind=8) :: d1g(mole%nb_inact2n,mole%nb_act1)
      real (kind=8) :: d2g(mole%nb_inact2n,mole%nb_act1,mole%nb_act1)

      real (kind=8) :: Qsym0(mole%nb_var)
      real (kind=8) :: step
      logical       :: deriv,num




      real (kind=8) :: Qact(mole%nb_act1)
      real (kind=8) :: d0f
      real (kind=8) :: d1f(mole%nb_act1)
      real (kind=8) :: d2f(mole%nb_act1,mole%nb_act1)
      integer       :: i,k,vi,ibo

      real (kind=8) :: d0gzmt(mole%nb_var)

c----- for debuging ----------------------------------
c     logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
c---------------------------------------------------------------------
      IF (debug) THEN
      write(6,*)
      write(6,*) 'BEGINNING d0d1d2_g'
      write(6,*) 'nb_var',mole%nb_var
      write(6,*) 'nb_act1',mole%nb_act1
      write(6,*) 'nb_inact22,nb_inact21',mole%nb_inact22,mole%nb_inact21
      write(6,*) 'nb_inact2n',mole%nb_inact2n
      write(6,*) 'deriv',deriv
      END IF

c---------------------------------------------------------------------
       Qact(1) = Qsym0(mole%liste_QactTOQsym(1))
c---------------------------------------------------------------------

      d0g(:) = 0.d0

c---------------------------------------------------------------------
       IF (debug) THEN
         write(6,*) 'd0g at Qact:',Qact
         write(6,*) d0g(:)
         write(6,*) 'END d0d1d2_g'
       END IF
c---------------------------------------------------------------------

       RETURN
       END
C================================================================
C    subroutine calculant la matrice hessienne
C    en fonction de x=cos(theta)
C================================================================
       SUBROUTINE sub_hessian (h)

       real*8 h

       h = 0.d0


       RETURN
       END
       SUBROUTINE H0_sym(h)

       real*8 h

       RETURN
       END
C================================================================
C    subroutine for the hessian matrix
C================================================================
      SUBROUTINE d0d1d2_h(d0h,d1h,d2h,Qsym,mole,deriv,num,step)

      USE mod_system
      USE mod_Tnum
      implicit none

c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

       real(kind=Rkind) :: Qsym(mole%nb_var)

       real(kind=Rkind)  :: step
       logical :: deriv,num

       real(kind=Rkind) :: d0h(mole%nb_inact2n,mole%nb_inact2n)
       real(kind=Rkind) :: d1h(mole%nb_inact2n,mole%nb_inact2n)
       real(kind=Rkind) :: d2h(mole%nb_inact2n,mole%nb_inact2n)


c----- for debuging ----------------------------------
      logical :: debug
      parameter (debug=.FALSE.)
c      parameter (debug=.TRUE.)
c---------------------------------------------------------------------
      STOP 'd0d1d2_h'

       END
C================================================================
C    analytical derivative (Qeq Qeq' Qeq" Qeq'") calculation
c    for the variable i_qsym
C================================================================
      SUBROUTINE calc_dnQflex(iq,dnQflex,Qact,nb_act,nderiv,it)
      USE mod_system
      USE mod_dnSVM
      IMPLICIT NONE

       integer :: iq,nb_act
       real (kind=Rkind) ::  Qact(nb_act)

       integer :: nderiv,it ! it number of the transformation

       TYPE (Type_dnS)   :: dnQflex

       real (kind=Rkind) :: d0req,d1req,d2req,d3req
       integer :: i
       integer, parameter :: nb_var=12
       integer, parameter :: liste_QactTOQsym(nb_var)=(/(i,i=1,nb_var)/)
       integer, parameter :: liste_QsymTOQact(nb_var)=(/(i,i=1,nb_var)/)

c----- for debuging ----------------------------------
      character (len=*), parameter :: name_sub='dnQflex'
      logical, parameter :: debug=.FALSE.
c     logical, parameter :: debug=.TRUE.
c----- for debuging ----------------------------------


c---------------------------------------------------------------------
      IF (debug) THEN
        write(6,*) 'BEGINNING ',name_sub
        write(6,*) 'nb_act',nb_act
        write(6,*) 'iq',iq
        CALL flush_perso(6)
      END IF
c---------------------------------------------------------------------

c---------------------------------------------------------------------
c      Qact value. Rq: only ONE active variable is possible
c---------------------------------------------------------------------
       IF (nb_act /= 1) THEN
         write(6,*) ' ERROR in ',name_sub
         write(6,*) ' the number of Active variable'
         write(6,*) ' should be 1. But nb_act =',nb_act
         STOP
       END IF

c---------------------------------------------------------------------
c---------------------------------------------------------------------

c---------------------------------------------------------------------
       dnQflex%d0                          = ZERO
       IF (nderiv >= 1)  dnQflex%d1(:)     = ZERO
       IF (nderiv >= 2)  dnQflex%d2(:,:)   = ZERO
       IF (nderiv >= 3)  dnQflex%d3(:,:,:) = ZERO
c---------------------------------------------------------------------

       CALL  subd0d1d2d3_Qeq(iq,
     *                         d0req,d1req,d2req,d3req,
     *                         liste_QactTOQsym,liste_QsymTOQact,
     *                         Qact,nderiv)

       IF (nderiv >= 0) dnQflex%d0 = d0req
       IF (nderiv >= 1) dnQflex%d1 = d1req
       IF (nderiv >= 2) dnQflex%d2 = d2req
       IF (nderiv == 3) dnQflex%d3 = d3req


       !write(6,*) 'Qflex',iq,Qact(1),dnQflex%d0

c---------------------------------------------------------------------
      IF (debug) THEN
        write(6,*) 'dnQflex : ',Qact
        CALL write_dnS(dnQflex,nderiv)
        write(6,*) 'END ',name_sub
        CALL flush_perso(6)
      END IF
c---------------------------------------------------------------------

      END
c================================================================
c    analytical derivative (Qeq Qeq' Qeq" Qeq'") calculation
c    for the variable i_qsym
c================================================================
      SUBROUTINE d0d1d2d3_Qeq(i_qsym,
     *                        d0req,d1req,d2req,d3req,
     *                        Qsym,mole,nderiv)

       
      USE mod_system
      USE mod_Tnum
      implicit none
         
c----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole

       integer :: i_qsym
       real(kind=Rkind) :: Qsym(mole%nb_var)

       integer :: nderiv
       integer :: i

       real(kind=Rkind) :: d0req
       real(kind=Rkind) :: d1req(mole%nb_act)
       real(kind=Rkind) :: d2req(mole%nb_act,mole%nb_act)
       real(kind=Rkind) :: d3req(mole%nb_act,mole%nb_act,mole%nb_act)

       real(kind=Rkind) :: d0zeq
       real(kind=Rkind) :: d1zeq(mole%nb_act)
       real(kind=Rkind) :: d2zeq(mole%nb_act,mole%nb_act)
       real(kind=Rkind) :: d3zeq(mole%nb_act,mole%nb_act,mole%nb_act)

      real (kind=Rkind) :: pi2,pi23
      parameter (pi2 =pi+pi)
      parameter (pi23 =pi2/3.d0)


       d0req = ZERO
       IF (nderiv .GE. 1)  d1req = ZERO
       IF (nderiv .GE. 2)  d2req = ZERO
       IF (nderiv .GE. 3)  d3req = ZERO
       STOP 'd0d1d2d3_Qeq'

      END
C================================================================
C    Fonction, pot0, Bowman potential
C================================================================
      real*8 FUNCTION pot0(Qxyz)
      USE mod_system
      implicit none

       integer, parameter :: ncc=18
       real (kind=8) :: Qxyz(ncc),nrg
       real (kind=8) :: xyz(6,3)

       logical, save :: begin = .TRUE.


c---------------------------------------------------------------
c      initialization (only once)
c$OMP    CRITICAL (pot0_CRIT)
       IF (begin) THEN
         CALL prepot()
         begin=.FALSE.
       END IF
c$OMP    END CRITICAL (pot0_CRIT)
c      END initialization
c---------------------------------------------------------------

       xyz(:,:) = transpose(reshape(Qxyz, (/3,6/) ))

       CALL calcpot(nrg,xyz)

       pot0 = nrg

       !write(6,*) 'Qxyz',Qxyz,nrg
       !STOP

       END
!********************************************************************************
!------------Note by Xinchuan Huang on 2005-04-08, updated in April, 2009--------
! Potential Energy Surface for CH3OH
! Intrinsic version 4, CANNOT describe HO...CH3 dissociation
!
! Citation/Reference for this PES
!  Joel M. Bowman, Xinchuan Huang, Nicholas C. Handy, Stuart Carter, J. Phys. Chem. A, 111, 7317 (2007).  
!
!  19,315 ab initio single points, computed at ccsd(t)/aug-cc-pvtz level with MOLPRO 2002.6
!
!  Original fitting codes were written by Bastiaan J. Bramms on 2004-01-29, using energies only
!  Codes were modified and applied to fit by Xinchuan Huang    ---- Finished on 2005-04-07
!
!     No.of.Points     Energy Cut    RMS fitting error:
!         6,371           500 cm-1     0.46 cm-1
!         8,742         2,500 cm-1     1.42 cm-1
!        10,564         5,000 cm-1     2.56 cm-1
!        13,269        10,000 cm-1     5.26 cm-1
!        15,038        15,000 cm-1     6.67 cm-1
!        16,538        20,000 cm-1     8.35 cm-1
!        18,308        30,000 cm-1     9.91 cm-1
!
!  pgf90 & ifc/ifort compatible
!  one data file is required :  ch3oh.pes4.coeff.dat
!
!  The Cs-symmetry minimum geometry on PES-4 in the order of C O H H H H is: (in bohr)
!          X                  Y                 Z
!C   0.000000000000000    0.00000000000000   0.00000000000000D+000
!O   2.69384866444459     0.00000000000000   0.00000000000000D+000
!H   3.25464844199020     0.00000000000000   1.72748736457417
!H  -0.590685613953303    0.00000000000000  -1.97213587090349
!H  -0.770030152497506    1.68581929552297   0.919398475020390
!H  -0.770030152497506   -1.68581929552297   0.919398475020390
!
!  The saddle point on PES-4 is found at 342.113517704766 cm-1 above the minimum
!  its cartesian coordinates in the order of C O H H H H is: (in bohr)
!          X                  Y                 Z
!C   0.00000000000000   0.00000000000000   0.00000000000000D+000
!O   2.70168258977627   0.00000000000000   0.00000000000000D+000
!H   3.27720953588083   0.00000000000000  -1.71828642403232
!H -0.773235753505252   0.00000000000000  -1.91454906757083
!H -0.686338573590350   1.67973201897574   0.984326851977137
!H -0.686338573590350  -1.67973201897574   0.984326851977136
!
!  prepot() should be called before first calling calcpot(V,xx)
!    xx(6,3) is cartesian array containing C O H H H H coor, in bohr
!    returned V is potential in hartree, taking Cs-sym minimum as zero potential reference
!
!  Permanent contact :  Joel M. Bowman    E-mail: jmbowma@emory.edu
!  Techinical questions may send to :  xinchuan@hotmail.com 
!-----------------------End of Note----------------------------------------------
!********************************************************************************
        subroutine prepot()
        use mod_file
        implicit none

        integer i,j,k,i1,j1,k1,m,mr,nio

        double precision dc0(0:5,0:5),dw0(0:5,0:5)
        double precision coef(0:3337)

        common/NCOE/m,mr
        common/ch3ohcoe/dc0,dw0,coef

! for complete 6th-order fit
        m=3250 ; mr=22
        CALL file_open2('ch3oh.pes4/ch3oh.pes4.coeff.dat',nio)
!       open(20,file='ch3oh.pes4/ch3oh.pes4.coeff.dat',status='old')

!        read(20,*)
!        read(20,*)dc0
        read(nio,*)
!        read(20,*)dw0
        read(nio,*)
        read(nio,*)(coef(i1),i1=0,m+4*mr-1)
!        write(*,*)(coef(i1),i1=0,m+4*mr-1)
        close(nio)

        return
        end 

!***********************************************************
      subroutine calcpot(V,xyz)

      implicit none

      double precision, dimension(6,3) :: xyz
      double precision, dimension(0:2,0:5) :: xn
      double precision, dimension(0:2,0:5) :: gf0
      double precision V

      integer i,j

      do i=1,3
        xn(i-1,4:5)=xyz(1:2,i)
        xn(i-1,0:3)=xyz(3:6,i)
      end do

      call potshell(xn,V)

! set MIN to zero-point of potential  ---- PES-4 min
      V=V+115.562385144095d0

!      V=V*219474.63067d0
      return
      end subroutine calcpot

!***********************************************************
      subroutine potshell(xn,V)

      implicit none
      
      integer m,mr

      double precision dc0(0:5,0:5),dw0(0:5,0:5)
      double precision coef(0:3337), f0, gf0(0:2,0:5)

      common/NCOE/m,mr
      common/ch3ohcoe/dc0,dw0,coef

      double precision V, xn(0:2,0:5)
      double precision d0(0:5,0:5), dw(0:5,0:5)
      double precision vec(0:m+4*mr-1)

!      call getfit(m,mr,dc0,dw0,coef,xn,f0,gf0,vec)

      call getvec(m,mr,xn,vec)
      f0=dot_product(coef,vec)

      V=f0

      return
      end 

!****************************************************************
       
      subroutine getd0 (nk, r0, d0)
      implicit none
      integer, parameter :: wp=selected_real_kind(12,300)
      integer nk
      double precision r0(0:nk-1,0:nk-1), d0(0:nk-1,0:nk-1)
      integer i, j
      do i = 0, nk-1
       d0(i,i) = 0
       do j = i+1, nk-1
        d0(i,j) = dexp(-r0(i,j)/3)
        d0(j,i) = d0(i,j)
       enddo
      enddo
      return
      end
!      subroutine getfit (ms, mr, coef, xn, f0, gf0, vec)
!      implicit none
!      integer, parameter :: wp=selected_real_kind(12,300)
!      integer nk, ms, mr
!      parameter (nk=6)
!      double precision coef(0:ms+4*mr-1), xn(0:2,0:nk-1), f0, 
!     $ gf0(0:2,0:nk-1), vec(0:ms+4*mr-1)
!      double precision dd
!      parameter (dd=1.0e-6)
!      integer i, j
!      double precision xn1(0:2,0:nk-1), t0, t1
!!----------------------------------------------------------------------
!      do i = 0, 2
!       do j = 0, nk-1
!        xn1 = xn ; xn1(i,j) = xn1(i,j)-dd
!        call getvec (ms, mr, xn, vec)
!        t0 = dot_product(coef,vec)
!        xn1 = xn ; xn1(i,j) = xn1(i,j)+dd
!        call getvec (ms, mr, xn, vec)
!        t1 = dot_product(coef,vec)
!! Note that gf0 will be the negative gradient
!        gf0(i,j) = -(t1-t0)/(2*dd)
!       enddo
!      enddo
!! f0 is last, so vec will contain a sensible return value
!      call getvec (ms, mr, xn, vec)
!      f0 = dot_product(coef,vec)
!      return
!      end
      subroutine getr0 (nk, xn, r0)
      implicit none
      integer, parameter :: wp=selected_real_kind(12,300)
      integer nk
      double precision xn(0:2,0:nk-1), r0(0:nk-1,0:nk-1)
      integer i, j
      do i = 0, nk-1
       r0(i,i) = 0
       do j = i+1, nk-1
        r0(i,j) = sqrt((xn(0,j)-xn(0,i))**2+(xn(1,j)-xn(1,i))**2+ 
     $      (xn(2,j)-xn(2,i))**2)
        r0(j,i) = r0(i,j)
       enddo
      enddo
      return
      end
      subroutine getrvec (ms, r, vec)
      implicit none
      integer, parameter :: wp=selected_real_kind(12,300)
! version for X4YZ
      integer nk, ms
      parameter (nk=6)
      double precision r(0:nk-1,0:nk-1), vec(0:ms-1)
      integer i, j
      double precision x(0:3), r1(0:nk-1,0:nk-1), t0, t1
!-----------------------------------------------------------------------
! Test for compatibility
      if (.not.(ms.eq.1.or.ms.eq.5)) then
       stop 'getrvec - wrong dimension'
      endif
! Computation
      x = 0
      do i = 0, nk-1
       do j = 0, nk-1
        if (i.eq.j) then
         r1(i,j) = 0
        else
         r1(i,j) = dexp(-r(i,j))/r(i,j)
        endif
       enddo
      enddo
! XX distance
      t0 = 0
      do i = 0, nk-3
       do j = i+1, nk-3
        t0 = t0+r1(i,j)
       enddo
      enddo
      x(0) = t0/((nk-2)*(nk-3)/2)
! XY and XZ distances
      t0 = 0 ; t1 = 0
      do i = 0, nk-3
       t0 = t0+r1(i,nk-2)
       t1 = t1+r1(i,nk-1)
      enddo
      x(1) = t0/(nk-2)
      x(2) = t1/(nk-2)
! YZ distance
      x(3) = r1(nk-2,nk-1)
! set vec
      vec(0) = 1
      if (5.le.ms) then
       vec(1:4) = x
      endif
      return
      end
      subroutine getvec (ms, mr, xn, vec)
      implicit none
      integer, parameter :: wp=selected_real_kind(12,300)
! version for H4CO (reordering of CH3OH).
      integer nk, ms, mr
      parameter (nk=6)
      double precision xn(0:2,0:nk-1), vec(0:ms+4*mr-1)
      integer k, l
      double precision rvec(0:4), d0(0:nk-1,0:nk-1), r0(0:nk-1,0:nk-1)
!-----------------------------------------------------------------------
      vec = 0
      call getr0 (nk, xn, r0)
      call getd0 (nk, r0, d0)
      call getv_x4yz (ms, d0, vec(0:ms-1))
      call getrvec (5, r0, rvec)
      do l = 0, mr-1
       do k = 0, 3
        vec(ms+4*l+k) = rvec(k+1)*vec(l)
       enddo
      enddo
      return
      end
      subroutine getv_x4yz (m, d, vec)
      implicit none
      integer, parameter :: wp=selected_real_kind(12,300)
      integer, parameter :: nk=6
      integer m
      double precision d(0:nk-1,0:nk-1), vec(0:m-1)
! version for molecule X4YZ
! MolienSeries(0:8): 1 4 17 65 230 736 2197 6093 15864
! #Primaries(1:8):   4 4 4 3 0 0 0 0
! #Secondaries(1:8): 0 3 13 32 62 129 221 335
      integer, parameter :: l0=1, l1=l0+4, l2=l1+17, l3=l2+65,l4=l3+230,
     $ l5=l4+736, l6=l5+2197, l7=l6+6093-221
!! We haven't incorporated the secondaries at degree 7
      integer, parameter :: np1=4, np2=4, np3=4, np4=3, np5=0, np6=0, 
     $   np7=0
      integer, parameter :: ns1=0, ns2=3, ns3=13, ns4=32, ns5=62, 
     $   ns6=129, ns7=221
      double precision x(0:np1-1), y(0:np2-1), z(0:np3-1), u(0:np4-1)
      double precision x2(0:np1-1), x3(0:np1-1),x4(0:np1-1),x5(0:np1-1),
     $   x6(0:np1-1), x7(0:np1-1), y2(0:np2-1), y3(0:np2-1),z2(0:np3-1)
      double precision ys(0:ns2-1),zs(0:ns3-1),us(0:ns4-1),
     $   vs(0:ns5-1), ws(0:ns6-1), w7s(0:ns7-1)
      integer mdeg, i, j, k, l, i0, i1, i2, i3, j0, k0
      double precision d2(0:nk-1,0:nk-1), d3(0:nk-1,0:nk-1), d4(0:nk-1,
     $ 0:nk-1),d5(0:nk-1,0:nk-1), d6(0:nk-1,0:nk-1), d7(0:nk-1,0:nk-1)
      double precision t0
      double precision pol2, pol3, pol4, pol5, pol6, pol7
      pol2(t0) = t0**2
      pol3(t0) = t0**3
      pol4(t0) = t0**4
      pol5(t0) = t0**5
      pol6(t0) = t0**6
      pol7(t0) = t0**7
!-----------------------------------------------------------------------
! Test for compatibility, set mdeg
      select case (m)
      case (l0)
       mdeg = 0
      case (l1)
       mdeg = 1
      case (l2)
       mdeg = 2
      case (l3)
       mdeg = 3
      case (l4)
       mdeg = 4
      case (l5)
       mdeg = 5
      case (l6)
       mdeg = 6
      case (l7)
       mdeg = 7
      case default
       stop 'getv - wrong dimension'
      endselect
! auxiliary distances
      do i = 0, nk-1
       do j = i+1, nk-1
        d2(i,j) = pol2(d(i,j))
        d2(j,i) = d2(i,j)
        d3(i,j) = pol3(d(i,j))
        d3(j,i) = d3(i,j)
        d4(i,j) = pol4(d(i,j))
        d4(j,i) = d4(i,j)
        d5(i,j) = pol5(d(i,j))
        d5(j,i) = d5(i,j)
        d6(i,j) = pol6(d(i,j))
        d6(j,i) = d6(i,j)
        d7(i,j) = pol7(d(i,j))
        d7(j,i) = d7(i,j)
       enddo
      enddo
! Primary Invariants
!      x = 0 ; y = 0 ; z = 0 ; u = 0 ; v = 0 ; w = 0 ; w7 = 0
       x = 0 ; y = 0 ; z = 0 ; u = 0 

      do i0 = 0, 3
       t0 = 0
       do i1 = 0, 3
       if (i1.ne.i0) then
        t0 = t0+d(i0,i1)/3
        x(0) = x(0)+d(i0,i1)/12
        y(0) = y(0)+d2(i0,i1)/12
        z(0) = z(0)+d3(i0,i1)/12
       endif
       enddo
       y(1) = y(1)+pol2(t0)/4
       z(1) = z(1)+pol3(t0)/4
       u(0) = u(0)+pol4(t0)/4
      enddo
      x(1) = sum(d(0:3,4))/4
      y(2) = sum(d2(0:3,4))/4
      z(2) = sum(d3(0:3,4))/4
      u(1) = sum(d4(0:3,4))/4
      x(2) = sum(d(0:3,5))/4
      y(3) = sum(d2(0:3,5))/4
      z(3) = sum(d3(0:3,5))/4
      u(2) = sum(d4(0:3,5))/4
      x(3) = d(4,5)
! Required powers
      do i = 0, np1-1
       x2(i) = pol2(x(i))
       x3(i) = pol3(x(i))
       x4(i) = pol4(x(i))
       x5(i) = pol5(x(i))
       x6(i) = pol6(x(i))
       x7(i) = pol7(x(i))
      enddo
      do i = 0, np2-1
       y2(i) = pol2(y(i))
       y3(i) = pol3(y(i))
      enddo
      do i = 0, np3-1
       z2(i) = pol2(z(i))
      enddo
! Secondary Invariants
!! ys(0:2), zs(0:12), us(0:31), vs(0:61), ws(0:128), w7s(0:220)
!! reducible: us(0:5), vs(0:38), ws(0:127), w7s(0:..)
      ys = 0 ; zs = 0 ; us = 0 ; vs = 0 ; ws = 0 ; w7s = 0
! Irreducible secondaries
      do i0 = 0, 3
       do i1 = 0, 3
       if (i1.ne.i0) then
        us(6) = us(6)+d4(i0,i1)/12
        vs(39) = vs(39)+d5(i0,i1)/12
        do i2 = 0, 3
        if (i2.ne.i1.and.i2.ne.i0) then
         zs(0) = zs(0)+d2(i0,i1)*d(i0,i2)/24
        endif
        enddo
       endif
       enddo
      enddo
      j0 = 4 ; k0 = 5
      do i0 = 0, 3
       ys(2) = ys(2)+d(i0,j0)*d(i0,k0)/4
       zs(8) = zs(8)+d2(i0,j0)*d(i0,k0)/4
       zs(11) = zs(11)+d(i0,j0)*d2(i0,k0)/4
       us(20) = us(20)+d3(i0,j0)*d(i0,k0)/4
       us(27) = us(27)+d2(i0,j0)*d2(i0,k0)/4
       us(30) = us(30)+d(i0,j0)*d3(i0,k0)/4
       do i1 = 0, 3
       if (i1.ne.i0) then
        ys(0) = ys(0)+d(i0,i1)*d(i0,j0)/12
        ys(1) = ys(1)+d(i0,i1)*d(i0,k0)/12
        zs(1) = zs(1)+d2(i0,i1)*d(i0,j0)/12
        zs(3) = zs(3)+d(i0,i1)*d2(i0,j0)/12
        zs(4) = zs(4)+d(i0,i1)*d(i0,j0)*d(i1,j0)/12
        zs(5) = zs(5)+d2(i0,i1)*d(i0,k0)/12
        zs(7) = zs(7)+d(i0,i1)*d(i0,j0)*d(i0,k0)/12
        zs(9) = zs(9)+d(i0,i1)*d(i1,j0)*d(i0,k0)/12
        zs(10) = zs(10)+d(i0,i1)*d2(i0,k0)/12
        zs(12) = zs(12)+d(i0,i1)*d(i0,k0)*d(i1,k0)/12
        us(7) = us(7)+d3(i0,i1)*d(i0,j0)/12
        us(10) = us(10)+d2(i0,i1)*d2(i0,j0)/12
        us(12) = us(12)+d(i0,i1)*d3(i0,j0)/12
        us(13) = us(13)+d2(i0,i1)*d(i0,j0)*d(i1,j0)/12
        us(14) = us(14)+d3(i0,i1)*d(i0,k0)/12
        us(17) = us(17)+d2(i0,i1)*d(i0,j0)*d(i0,k0)/12
        us(19) = us(19)+d(i0,i1)*d2(i0,j0)*d(i0,k0)/12
        us(21) = us(21)+d2(i0,i1)*d(i1,j0)*d(i0,k0)/12
        us(23) = us(23)+d(i0,i1)*d(i0,j0)*d(i1,j0)*d(i0,k0)/12
        us(24) = us(24)+d2(i0,i1)*d2(i0,k0)/12
        us(26) = us(26)+d(i0,i1)*d(i0,j0)*d2(i0,k0)/12
        us(28) = us(28)+d(i0,i1)*d(i1,j0)*d2(i0,k0)/12
        us(29) = us(29)+d(i0,i1)*d3(i0,k0)/12
        us(31) = us(31)+d2(i0,i1)*d(i0,k0)*d(i1,k0)/12
        vs(40) = vs(40)+d4(i0,i1)*d(i0,j0)/12
        vs(42) = vs(42)+d3(i0,i1)*d2(i0,j0)/12
        vs(44) = vs(44)+d2(i0,i1)*d3(i0,j0)/12
        vs(45) = vs(45)+d(i0,i1)*d3(i0,j0)*d(i1,j0)/12
        vs(46) = vs(46)+d4(i0,i1)*d(i0,k0)/12
        vs(48) = vs(48)+d3(i0,i1)*d(i0,j0)*d(i0,k0)/12
        vs(50) = vs(50)+d2(i0,i1)*d2(i0,j0)*d(i0,k0)/12
        vs(52) = vs(52)+d2(i0,i1)*d(i0,j0)*d(i1,j0)*d(i0,k0)/12
        vs(53) = vs(53)+d(i0,i1)*d2(i0,j0)*d(i1,j0)*d(i0,k0)/12
        vs(54) = vs(54)+d3(i0,i1)*d2(i0,k0)/12
        vs(56) = vs(56)+d2(i0,i1)*d(i0,j0)*d2(i0,k0)/12
        vs(57) = vs(57)+d2(i0,i1)*d(i1,j0)*d2(i0,k0)/12
        vs(58) = vs(58)+d(i0,i1)*d(i0,j0)*d(i1,j0)*d2(i0,k0)/12
        vs(59) = vs(59)+d2(i0,i1)*d3(i0,k0)/12
        vs(60) = vs(60)+d(i0,i1)*d(i1,j0)*d3(i0,k0)/12
        vs(61) = vs(61)+d(i0,i1)*d3(i0,k0)*d(i1,k0)/12
        do i2 = 0, 3
        if (i2.ne.i1.and.i2.ne.i0) then
         zs(2) = zs(2)+d(i0,i1)*d(i0,i2)*d(i0,j0)/24
         zs(6) = zs(6)+d(i0,i1)*d(i0,i2)*d(i0,k0)/24
         us(8) = us(8)+d2(i0,i1)*d(i0,i2)*d(i0,j0)/24
         us(9) = us(9)+d2(i0,i1)*d(i1,i2)*d(i0,j0)/24
         us(11) = us(11)+d(i0,i1)*d(i0,i2)*d2(i0,j0)/24
         us(15) = us(15)+d2(i0,i1)*d(i0,i2)*d(i0,k0)/24
         us(16) = us(16)+d2(i0,i1)*d(i1,i2)*d(i0,k0)/24
         us(18) = us(18)+d(i0,i1)*d(i0,i2)*d(i0,j0)*d(i0,k0)/24
         us(22) = us(22)+d(i0,i1)*d(i0,i2)*d(i1,j0)*d(i0,k0)/24
         us(25) = us(25)+d(i0,i1)*d(i0,i2)*d2(i0,k0)/24
         vs(41) = vs(41)+d3(i0,i1)*d(i0,i2)*d(i0,j0)/24
         vs(43) = vs(43)+d2(i0,i1)*d(i0,i2)*d2(i0,j0)/24
         vs(47) = vs(47)+d3(i0,i1)*d(i0,i2)*d(i0,k0)/24
         vs(49) = vs(49)+d2(i0,i1)*d(i0,i2)*d(i0,j0)*d(i0,k0)/24
         vs(51) = vs(51)+d2(i0,i1)*d(i0,i2)*d(i1,j0)*d(i0,k0)/24
         vs(55) = vs(55)+d2(i0,i1)*d(i0,i2)*d2(i0,k0)/24
         ws(128) = ws(128)+d2(i0,i1)*d2(i0,i2)*d(i1,j0)*d(i0,k0)/24
        endif
        enddo
       endif
       enddo
      enddo
! Reducible secondaries
! us(0:5)
      k = 0
      do j = 0, 2
       do i = 0, j-1
        us(k) = ys(i)*ys(j) ; k = k+1
       enddo
       us(k) = pol2(ys(j)) ; k = k+1
      enddo
! vs(0:38)
      do j = 0, 12
       do i = 0, 2
        vs(3*j+i) = ys(i)*zs(j)
       enddo
      enddo
! ws(0:127)
      ws(0) = pol2(zs(0))
      ws(1) = zs(0)*zs(1)
      ws(2) = zs(0)*zs(2)
      ws(3) = zs(0)*zs(3)
      ws(4) = zs(1)*zs(3)
      ws(5) = zs(2)*zs(3)
      ws(6) = zs(0)*zs(4)
      ws(7) = zs(1)*zs(4)
      ws(8) = zs(3)*zs(4)
      ws(9) = pol2(zs(4))
      ws(10) = zs(0)*zs(5)
      ws(11) = zs(0)*zs(6)
      ws(12) = zs(1)*zs(6)
      ws(13) = zs(0)*zs(7)
      ws(14) = zs(1)*zs(7)
      ws(15) = zs(0)*zs(8)
      ws(16) = zs(1)*zs(8)
      ws(17) = zs(2)*zs(8)
      ws(18) = zs(4)*zs(8)
      ws(19) = zs(0)*zs(9)
      ws(20) = zs(1)*zs(9)
      ws(21) = zs(2)*zs(9)
      ws(22) = zs(3)*zs(9)
      ws(23) = zs(4)*zs(9)
      ws(24) = zs(0)*zs(10)
      ws(25) = zs(1)*zs(10)
      ws(26) = zs(4)*zs(10)
      ws(27) = zs(5)*zs(10)
      ws(28) = zs(6)*zs(10)
      ws(29) = zs(0)*zs(11)
      ws(30) = zs(1)*zs(11)
      ws(31) = zs(2)*zs(11)
      ws(32) = zs(3)*zs(11)
      ws(33) = zs(4)*zs(11)
      ws(34) = zs(5)*zs(11)
      ws(35) = zs(6)*zs(11)
      ws(36) = zs(7)*zs(11)
      ws(37) = zs(8)*zs(11)
      ws(38) = zs(0)*zs(12)
      ws(39) = zs(1)*zs(12)
      ws(40) = zs(2)*zs(12)
      ws(41) = zs(3)*zs(12)
      ws(42) = zs(4)*zs(12)
      ws(43) = zs(5)*zs(12)
      ws(44) = zs(7)*zs(12)
      ws(45) = zs(8)*zs(12)
      ws(46) = zs(9)*zs(12)
      ws(47) = zs(10)*zs(12)
      ws(48) = zs(11)*zs(12)
      ws(49) = pol2(zs(12))
      do j = 0, 25
       do i = 0, 2
        ws(50+3*j+i) = ys(i)*us(6+j)
       enddo
      enddo
! w7s(0:..) !! to follow
! Compute vec(0:*).  The code was created using these parameters:
! MolienSeries(0:*): 1 4 17 65 230 736 2197 6093
! #Primaries(1:*):   4 4 4 3 0 0 0
! #Secondaries(1:*): 0 3 13 32 62 129 221
! The MolienSeries partial sums (allowed size of vec) are:
! Molien Sums(0:*):  1 5 22 87 317 1053 3250 9343 9347
! constant term
      vec(0) = 1
! first degree terms
      if (1.le.mdeg) then
       vec(1) = x(0)
       vec(2) = x(1)
       vec(3) = x(2)
       vec(4) = x(3)
      endif
! second degree terms
      if (2.le.mdeg) then
       vec(5) = x2(0)
       vec(6:8) = x(0)*vec(2:4)
       vec(9) = x2(1)
       vec(10:11) = x(1)*vec(3:4)
       vec(12) = x2(2)
       vec(13:13) = x(2)*vec(4:4)
       vec(14) = x2(3)
       vec(15:14) = x(3)*vec(5:4)
       vec(15) = y(0)
       vec(16) = y(1)
       vec(17) = y(2)
       vec(18) = y(3)
       vec(19:21) = ys(0:2)
      endif
! third degree terms
      if (3.le.mdeg) then
       vec(22) = x3(0)
       vec(23:25) = x2(0)*vec(2:4)
       vec(26:38) = x(0)*vec(9:21)
       vec(39) = x3(1)
       vec(40:41) = x2(1)*vec(3:4)
       vec(42:51) = x(1)*vec(12:21)
       vec(52) = x3(2)
       vec(53:53) = x2(2)*vec(4:4)
       vec(54:61) = x(2)*vec(14:21)
       vec(62) = x3(3)
       vec(63:62) = x2(3)*vec(5:4)
       vec(63:69) = x(3)*vec(15:21)
       vec(70) = z(0)
       vec(71) = z(1)
       vec(72) = z(2)
       vec(73) = z(3)
       vec(74:86) = zs(0:12)
      endif
! fourth degree terms
      if (4.le.mdeg) then
       vec(87) = x4(0)
       vec(88:90) = x3(0)*vec(2:4)
       vec(91:103) = x2(0)*vec(9:21)
       vec(104:151) = x(0)*vec(39:86)
       vec(152) = x4(1)
       vec(153:154) = x3(1)*vec(3:4)
       vec(155:164) = x2(1)*vec(12:21)
       vec(165:199) = x(1)*vec(52:86)
       vec(200) = x4(2)
       vec(201:201) = x3(2)*vec(4:4)
       vec(202:209) = x2(2)*vec(14:21)
       vec(210:234) = x(2)*vec(62:86)
       vec(235) = x4(3)
       vec(236:235) = x3(3)*vec(5:4)
       vec(236:242) = x2(3)*vec(15:21)
       vec(243:259) = x(3)*vec(70:86)
       vec(260) = y2(0)
       vec(261:266) = y(0)*vec(16:21)
       vec(267) = y2(1)
       vec(268:272) = y(1)*vec(17:21)
       vec(273) = y2(2)
       vec(274:277) = y(2)*vec(18:21)
       vec(278) = y2(3)
       vec(279:281) = y(3)*vec(19:21)
       vec(282) = u(0)
       vec(283) = u(1)
       vec(284) = u(2)
       vec(285:316) = us(0:31)
      endif
! fifth degree terms
      if (5.le.mdeg) then
       vec(317) = x5(0)
       vec(318:320) = x4(0)*vec(2:4)
       vec(321:333) = x3(0)*vec(9:21)
       vec(334:381) = x2(0)*vec(39:86)
       vec(382:546) = x(0)*vec(152:316)
       vec(547) = x5(1)
       vec(548:549) = x4(1)*vec(3:4)
       vec(550:559) = x3(1)*vec(12:21)
       vec(560:594) = x2(1)*vec(52:86)
       vec(595:711) = x(1)*vec(200:316)
       vec(712) = x5(2)
       vec(713:713) = x4(2)*vec(4:4)
       vec(714:721) = x3(2)*vec(14:21)
       vec(722:746) = x2(2)*vec(62:86)
       vec(747:828) = x(2)*vec(235:316)
       vec(829) = x5(3)
       vec(830:829) = x4(3)*vec(5:4)
       vec(830:836) = x3(3)*vec(15:21)
       vec(837:853) = x2(3)*vec(70:86)
       vec(854:910) = x(3)*vec(260:316)
       vec(911:927) = y(0)*vec(70:86)
       vec(928:944) = y(1)*vec(70:86)
       vec(945:961) = y(2)*vec(70:86)
       vec(962:978) = y(3)*vec(70:86)
       vec(979:981) = z(0)*ys(0:2)
       vec(982:984) = z(1)*ys(0:2)
       vec(985:987) = z(2)*ys(0:2)
       vec(988:990) = z(3)*ys(0:2)
       vec(991:1052) = vs(0:61)
      endif
! sixth degree terms
      if (6.le.mdeg) then
       vec(1053) = x6(0)
       vec(1054:1056) = x5(0)*vec(2:4)
       vec(1057:1069) = x4(0)*vec(9:21)
       vec(1070:1117) = x3(0)*vec(39:86)
       vec(1118:1282) = x2(0)*vec(152:316)
       vec(1283:1788) = x(0)*vec(547:1052)
       vec(1789) = x6(1)
       vec(1790:1791) = x5(1)*vec(3:4)
       vec(1792:1801) = x4(1)*vec(12:21)
       vec(1802:1836) = x3(1)*vec(52:86)
       vec(1837:1953) = x2(1)*vec(200:316)
       vec(1954:2294) = x(1)*vec(712:1052)
       vec(2295) = x6(2)
       vec(2296:2296) = x5(2)*vec(4:4)
       vec(2297:2304) = x4(2)*vec(14:21)
       vec(2305:2329) = x3(2)*vec(62:86)
       vec(2330:2411) = x2(2)*vec(235:316)
       vec(2412:2635) = x(2)*vec(829:1052)
       vec(2636) = x6(3)
       vec(2637:2636) = x5(3)*vec(5:4)
       vec(2637:2643) = x4(3)*vec(15:21)
       vec(2644:2660) = x3(3)*vec(70:86)
       vec(2661:2717) = x2(3)*vec(260:316)
       vec(2718:2859) = x(3)*vec(911:1052)
       vec(2860) = y3(0)
       vec(2861:2866) = y2(0)*vec(16:21)
       vec(2867:2916) = y(0)*vec(267:316)
       vec(2917) = y3(1)
       vec(2918:2922) = y2(1)*vec(17:21)
       vec(2923:2966) = y(1)*vec(273:316)
       vec(2967) = y3(2)
       vec(2968:2971) = y2(2)*vec(18:21)
       vec(2972:3010) = y(2)*vec(278:316)
       vec(3011) = y3(3)
       vec(3012:3014) = y2(3)*vec(19:21)
       vec(3015:3049) = y(3)*vec(282:316)
       vec(3050) = z2(0)
       vec(3051:3066) = z(0)*vec(71:86)
       vec(3067) = z2(1)
       vec(3068:3082) = z(1)*vec(72:86)
       vec(3083) = z2(2)
       vec(3084:3097) = z(2)*vec(73:86)
       vec(3098) = z2(3)
       vec(3099:3111) = z(3)*vec(74:86)
       vec(3112:3114) = u(0)*ys(0:2)
       vec(3115:3117) = u(1)*ys(0:2)
       vec(3118:3120) = u(2)*ys(0:2)
       vec(3121:3249) = ws(0:128)
      endif
! seventh degree terms
      if (7.le.mdeg) then
       vec(3250) = x7(0)
       vec(3251:3253) = x6(0)*vec(2:4)
       vec(3254:3266) = x5(0)*vec(9:21)
       vec(3267:3314) = x4(0)*vec(39:86)
       vec(3315:3479) = x3(0)*vec(152:316)
       vec(3480:3985) = x2(0)*vec(547:1052)
       vec(3986:5446) = x(0)*vec(1789:3249)
       vec(5447) = x7(1)
       vec(5448:5449) = x6(1)*vec(3:4)
       vec(5450:5459) = x5(1)*vec(12:21)
       vec(5460:5494) = x4(1)*vec(52:86)
       vec(5495:5611) = x3(1)*vec(200:316)
       vec(5612:5952) = x2(1)*vec(712:1052)
       vec(5953:6907) = x(1)*vec(2295:3249)
       vec(6908) = x7(2)
       vec(6909:6909) = x6(2)*vec(4:4)
       vec(6910:6917) = x5(2)*vec(14:21)
       vec(6918:6942) = x4(2)*vec(62:86)
       vec(6943:7024) = x3(2)*vec(235:316)
       vec(7025:7248) = x2(2)*vec(829:1052)
       vec(7249:7862) = x(2)*vec(2636:3249)
       vec(7863) = x7(3)
       vec(7864:7863) = x6(3)*vec(5:4)
       vec(7864:7870) = x5(3)*vec(15:21)
       vec(7871:7887) = x4(3)*vec(70:86)
       vec(7888:7944) = x3(3)*vec(260:316)
       vec(7945:8086) = x2(3)*vec(911:1052)
       vec(8087:8476) = x(3)*vec(2860:3249)
       vec(8477:8493) = y2(0)*vec(70:86)
       vec(8494:8618) = y(0)*vec(928:1052)
       vec(8619:8635) = y2(1)*vec(70:86)
       vec(8636:8743) = y(1)*vec(945:1052)
       vec(8744:8760) = y2(2)*vec(70:86)
       vec(8761:8851) = y(2)*vec(962:1052)
       vec(8852:8868) = y2(3)*vec(70:86)
       vec(8869:8942) = y(3)*vec(979:1052)
       vec(8943:8977) = z(0)*vec(282:316)
       vec(8978:9012) = z(1)*vec(282:316)
       vec(9013:9047) = z(2)*vec(282:316)
       vec(9048:9082) = z(3)*vec(282:316)
       vec(9083:9095) = u(0)*zs(0:12)
       vec(9096:9108) = u(1)*zs(0:12)
       vec(9109:9121) = u(2)*zs(0:12)
!! vec(9122:9342) = w7s(0:220)
!! Secondaries at degree 7 yet to follow
      endif
      return
      end
      
       
      
      

C================================================================
C    analytical derivative (Qeq Qeq' Qeq" Qeq'") calculation
c    for the variable i_qsym
c    Ordre des coordonnees (i_qym):
c
c    1 => CO
c    2 => OH
c    3 => HOC
c    4,5 => ROH1, OCH1
c    6   => torsion: 1/3(phi1+phi2+phi3)  (active)
c    7,8 => ROH2, OCH2
c    9   => a*(phi2-phi3)
c    10,11 => ROH3, OCH3
c    12    => b*(2phi1-phi2-phi3)
c
c MODIFIER POUR AVOIR ANGLE DIEDRE EN CA
C================================================================
      SUBROUTINE subd0d1d2d3_Qeq(i_qsym,
     *                           d0Qop,d1Qop,d2Qop,d3Qop,
     *                           liste_QactTOQsym,liste_QsymTOQact,
     *                           Qsym0,nderiv)
      IMPLICIT NONE

       integer nb_var,nb_act,nb_act1,nb_inact
       parameter (nb_var=12)
       parameter (nb_act=12)
       parameter (nb_inact=11)
       parameter (nb_act1=1)
       integer liste_QactTOQsym(nb_var)
       integer liste_QsymTOQact(nb_var)

       integer i_qact,i_qsym
       real*8  Qact1,Qsym0(nb_var)

       integer nderiv

       real*8  d0Qop,d1Qop,d2Qop,d3Qop


      real*8 pi,pi2,pi23
      parameter (pi = 3.141592653589793238462643383279
     *                          50288419716939937511d0)
      parameter (pi2 =pi+pi)
      parameter (pi23 =pi2/3.d0)

       integer i,k,kl,ih,ifo
       integer i_sym_act,i_qsym_inact


       real*8  d0,d1,d2,d3

       integer   max_points,nb_inactb
       parameter (max_points=3)
       parameter (nb_inactb=12)
       real*8 F(max_points,nb_inactb)
       integer nn(nb_inactb)
       data nn /2,2,2,2,2,0,2,2,2,3,3,3/
       integer ntyp(nb_inactb)
       data ntyp/43,43,43,43,43,0, 28,28,28,27,27,27 /


       !parameters obtained from the minimun and the TS of Bowman pot
       data (F(k, 1 ),k=1, 2 ) /6.7622955939d0,6.9426396d-3/
       data (F(k, 2 ),k=1, 2 ) /4.5474547926d0,-3.6563152d-3/
       data (F(k, 3 ),k=1, 2 ) /4.7358865624d0,8.2400303d-3/
       data (F(k, 4 ),k=1, 2 ) /5.1761441937d0,-8.80d-4/
       data (F(k, 5 ),k=1, 2 ) /4.8211810965d0,2.3425337d-3/

        data (F(k, 7 ),k=1, 2 ) /-7.7010927d-3,-7.0239467d-3/
        data (F(k, 8 ),k=1, 2 ) /-9.78875318d-2,-3.28873691d-2/
        data (F(k, 9 ),k=1, 2 ) /3.8577602d-2,2.23729421d-2/

        data (F(k,10 ),k=1, 3 ) /0.d0, 7.7010927d-3,-7.0239467d-3/
        data (F(k,11 ),k=1, 3 ) /0.d0, 9.78875318d-2,-3.28873691d-2/
        data (F(k,12 ),k=1, 3 ) /0.d0,3.8577602d-2,-2.23729421d-2/

c----- for debuging ----------------------------------
      logical debug
      parameter (debug=.FALSE.)
c     parameter (debug=.TRUE.)
c----- for debuging ----------------------------------

       SAVE F,nn

c---------------------------------------------------------------------
      IF (debug) THEN
        write(6,*) 'BEGINNING d0d1d2d3_Qeq'
        write(6,*) 'nb_act,nb_var',nb_act,nb_var
        write(6,*) 'i_qsym',i_qsym
      END IF
c---------------------------------------------------------------------


c---------------------------------------------------------------------
       i_qact = nb_act1
       i_sym_act = liste_QactTOQsym(1)
       Qact1  = Qsym0(i_sym_act)
c      write(6,*) 'i_qact i_sym_act Qact1',i_qact,i_sym_act,Qact1

c---------------------------------------------------------------------


c      --------------------------------------------------------------
       d0Qop = 0.d0
       d1Qop = 0.d0
       d2Qop = 0.d0
       d3Qop = 0.d0
       DO i=1,nn(i_Qsym)
            IF (ntyp(i_Qsym) == 43) THEN
              ! cos(3 k phi), i=1 => if=1, i=2 => if=7
              ! cos(3 k phi), i=1 => if= 6*(i-1) + 1
              ifo = 6*i-5
            ELSE IF (ntyp(i_Qsym) == 27) THEN
              ! cos(k phi), i=1 => if=1, i=2 => if=3
              ! cos(k phi), i=1 => if= 2*(i-1) + 1
              ifo = 2*i-1
            ELSE IF (ntyp(i_Qsym) == 28) THEN
              ! sin(k phi), i=1 => if=2, i=2 => if=4
              ! sin(k phi), i=1 => if= 2*i
              ifo = 2*i
            ELSE
              write(6,*) ' ERROR in d0d1d2d3_Qeq'
              write(6,*) ' ntyp unknown for i_qsym',ntyp(i_Qsym),i_Qsym
              STOP
            END IF
            CALL d0d1d2d3fourier(Qact1,d0,d1,d2,d3,ifo)
            !write(6,*) 'i_Qsym,i,if',i_Qsym,i,ifo,d0

            d0Qop = d0Qop + d0 * F(i,i_Qsym)
            IF (nderiv > 0) d1Qop = d1Qop + d1 * F(i,i_Qsym)
            IF (nderiv > 1) d2Qop = d2Qop + d2 * F(i,i_Qsym)
            IF (nderiv > 2) d3Qop = d3Qop + d3 * F(i,i_Qsym)
       END DO
       IF (i_Qsym == 12) d0Qop = d0Qop + sqrt(0.5d0)*4.d0/3.d0*pi
       !write(6,*) 'd0Qop : ',i_Qsym,Qact1,d0Qop

c---------------------------------------------------------------------
      IF (debug) THEN
        write(6,*) 'd0Qop : ',Qact1,d0Qop
        write(6,*) 'd1Qop : ',Qact1,d1Qop
        write(6,*) 'd2Qop : ',Qact1,d2Qop
        write(6,*) 'd3Qop : ',Qact1,d3Qop
        write(6,*) 'END d0d1d2d3_Qeq'
      END IF
c---------------------------------------------------------------------

       RETURN
       END
C================================================================
C    The tri component of the dipole moment.
C================================================================
       SUBROUTINE sub_dipole(dip,Q)

       real*8 Q(21)
       real*8 dip(3)

c     ----------------------------------------------------------------
      real*8 pi,pi2
      parameter (pi = 3.141592653589793238462643383279
     *                          50288419716939937511d0)
      parameter (pi2 =pi+pi)

      real*8 sq2pi,sqpi
c     ----------------------------------------------------------------

       real*8 x,z
       integer kl


       character*14 nom
       logical exist

       integer max_points,nn
       parameter (max_points=30)
       real*8 Fx(max_points)
       real*8 Fy(max_points)
       real*8 Fz(max_points)

       logical begin
       data begin/.true./
       SAVE begin,Fx,Fy,Fz,nn

       dip(:) = 0.d0

       END
