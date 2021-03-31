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
      MODULE mod_BFGS
      USE mod_system

      IMPLICIT NONE

      TYPE param_BFGS

      integer           :: max_iteration       = 10

      real (kind=Rkind) :: max_grad            = 0.000015_Rkind
      real (kind=Rkind) :: RMS_grad            = 0.000010_Rkind

      real (kind=Rkind) :: max_step            = 0.000060_Rkind
      real (kind=Rkind) :: RMS_step            = 0.000040_Rkind

      logical           :: calc_hessian        = .FALSE.

      real (kind=Rkind), pointer :: hessian_inv_init(:,:) => null()


      END TYPE param_BFGS

      CONTAINS

      SUBROUTINE Read_param_BFGS(para_BFGS)
      TYPE (param_BFGS), intent(inout) :: para_BFGS

      integer           :: max_iteration

      real (kind=Rkind) :: max_grad
      real (kind=Rkind) :: RMS_grad

      real (kind=Rkind) :: max_step
      real (kind=Rkind) :: RMS_step
      logical           :: calc_hessian


      integer :: err_io
      NAMELIST /BFGS/ max_iteration,max_grad,RMS_grad,max_step,RMS_step,calc_hessian

        max_iteration       = 10

        max_grad            = 0.000015_Rkind
        RMS_grad            = 0.000010_Rkind

        max_step            = 0.000060_Rkind
        RMS_step            = 0.000040_Rkind
        calc_hessian        = .FALSE.

        read(in_unitp,BFGS,IOSTAT=err_io)
        IF (err_io /= 0) THEN
           write(out_unitp,*) ' WARNING in Read_param_BFGS'
           write(out_unitp,*) '  while reading the "BFGS" namelist'
           write(out_unitp,*) ' end of file or end of record'
           write(out_unitp,*) ' Check your data !!'
           STOP
        END IF
        IF (print_level > 1) write(out_unitp,BFGS)

        para_BFGS%max_iteration       = max_iteration

        para_BFGS%max_grad            = max_grad
        para_BFGS%RMS_grad            = RMS_grad

        para_BFGS%max_step            = max_step
        para_BFGS%RMS_step            = RMS_step

        para_BFGS%calc_hessian        = calc_hessian

      END SUBROUTINE Read_param_BFGS

      SUBROUTINE Write_param_BFGS(para_BFGS)
      TYPE (param_BFGS), intent(in)   :: para_BFGS

      write(out_unitp,*) '  WRITE param_BFGS'
      write(out_unitp,*)
      write(out_unitp,*) '  max_iteration  ',para_BFGS%max_iteration
      write(out_unitp,*)
      write(out_unitp,*) '  max_grad       ',para_BFGS%max_grad
      write(out_unitp,*) '  RMS_grad       ',para_BFGS%RMS_grad
      write(out_unitp,*) '  max_step       ',para_BFGS%max_step
      write(out_unitp,*) '  RMS_step       ',para_BFGS%RMS_step
      write(out_unitp,*)
      write(out_unitp,*) '  END WRITE param_BFGS'

      END SUBROUTINE Write_param_BFGS

      SUBROUTINE Sub_BFGS(BasisnD,xOpt_min,SQ,nb_Opt,                   &
                          para_Tnum,mole,PrimOp,Qact,para_BFGS)

      USE mod_system
      USE mod_Coord_KEO
      USE mod_PrimOp
      USE mod_basis
      USE mod_Op
      USE mod_Auto_Basis
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole
      TYPE (Tnum)    :: para_Tnum
      logical        :: Cart_Transfo_save
      real (kind=Rkind), intent(inout) :: Qact(:)


!----- for the basis set ----------------------------------------------
      TYPE (basis)          :: BasisnD

!----- variables pour la namelist minimum ----------------------------
      TYPE (PrimOp_t)  :: PrimOp
      integer          :: nb_scalar_Op
      logical          :: calc_scalar_Op

!----- variables for the construction of H ---------------------------
      TYPE (param_ReadOp)         :: para_ReadOp
      logical                     :: Save_FileGrid,Save_MemGrid


!----- local variables -----------------------------------------------
!----- variables for the construction of H ---------------------------
      TYPE (param_AllOp)          :: para_AllOp_loc

!----- for the basis set ----------------------------------------------
      TYPE (param_AllBasis)       :: para_AllBasis_loc
      TYPE (basis)                :: basis_temp

!----- for the optimization -------------------------------------------
      TYPE (param_BFGS) :: para_BFGS
      integer, intent(in) :: nb_Opt
      real (kind=Rkind), intent(inout) :: xOpt_min(nb_Opt),SQ(nb_Opt)

!---------- working variables ----------------------------------------
  TYPE (param_dnMatOp) :: dnMatOp(1)
  integer        :: nderiv_alloc

!---------------------------------------------------------------------
!      logical,parameter :: debug= .FALSE.
      logical,parameter :: debug= .TRUE.
      character (len=*), parameter :: name_sub='Sub_BFGS'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'xopt_min',xopt_min(:)
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------

        Qact(1:nb_Opt) = xopt_min(:)


        write(out_unitp,*) 'Qact',Qact
        !---------JML-------------------------------------
        write(out_unitp,*) 'mole%nb_act', mole%nb_act
        write(out_unitp,*) 'RMS_grad', para_BFGS%RMS_grad
        write(out_unitp,*) 'RMS_step', para_BFGS%RMS_step
        write(out_unitp,*) 'max_iteration',  para_BFGS%max_iteration

        IF (para_BFGS%calc_hessian) THEN
          write(out_unitp,*) ' The initial hessian is calculated'
          !-------- allocation -----------------------------------------------
          CALL Init_Tab_OF_dnMatOp(dnMatOp,nb_Opt,PrimOp%nb_elec,nderiv=2)
          CALL alloc_array(para_BFGS%hessian_inv_init,(/ nb_Opt,nb_Opt /),  &
                          'para_BFGS%hessian_inv_init',name_sub)
          !-------- end allocation --------------------------------------------

          !----- Hessian ------------------------------------

          CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole,para_Tnum,PrimOp)

          CALL Get_Hess_FROM_Tab_OF_dnMatOp(para_BFGS%hessian_inv_init,dnMatOp) ! for the ground state

          CALL Write_Mat(para_BFGS%hessian_inv_init,out_unitp,5)

          !-------- deallocation ---------------------------------------------
          CALL dealloc_Tab_OF_dnMatOp(dnMatOp)
          !-------- end deallocation -----------------------------------------
        END IF

        !-------- allocation -----------------------------------------------
        CALL Init_Tab_OF_dnMatOp(dnMatOp,nb_Opt,PrimOp%nb_elec,nderiv=1)
        !-------- end allocation --------------------------------------------

        CALL dfpmin_new(Qact,dnMatOp(1)%tab_dnMatOp(:,:,1),             &
                        mole,PrimOp,para_Tnum,para_BFGS,              &
          para_BFGS%RMS_grad,para_BFGS%RMS_step,para_BFGS%max_iteration)

        xopt_min(:) = Qact(1:nb_Opt)

        !---------JMLend-------------------------------------
        !--------------------------------------------------

        !--------------------------------------------------
        ! this subroutine print the  matrix of derived type.
        CALL Write_MatOFdnS(dnMatOp(1)%tab_dnMatOp(:,:,1))
        !--------------------------------------------------


        !-------- deallocation ---------------------------------------------
        CALL dealloc_Tab_OF_dnMatOp(dnMatOp)
        IF (associated(para_BFGS%hessian_inv_init)) THEN
          CALL dealloc_array(para_BFGS%hessian_inv_init,                &
                            'para_BFGS%hessian_inv_init',name_sub)
        END IF
        !-------- end deallocation -----------------------------------------


!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


      END SUBROUTINE Sub_BFGS

!!!!!!!!!!!!!!!!!!!JML!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!---------------------------------------------------------------------------
SUBROUTINE dfpmin_new(Qact,MatdnE,mole,PrimOp,para_Tnum,para_BFGS,    &
                      gtol,tolx,itmax)
!---------------------------------------------------------------------------
!
 USE mod_system
 USE mod_Coord_KEO
 USE mod_PrimOp
 USE mod_basis
 USE mod_Op
 USE mod_Auto_Basis
!
!
!  Given a starting point p(1:n) that is a vector of length n, Broyden-Fletxer-
! Goldfrab-Shano varian of Davidon-Fletcher-Powell minimization is performed
! on a subroutine dfunc(...0), using this gradient as calculated by a routine dfunc(...1).
! The convergence requirement on zeroing the gradient by is input as gtol.
! Returned quantities are p(1:n) (the location of the minimum, iter (the number
! of iterations that were performed), and fret (the minimum value of the
! function). The routine lnsrch is called to perform approximate line
! minimizations, Parameters ITMAX
! is the maximum allowed number of iterations; STPMX is the scaled maximun step
! lengh allowed in the line search; TOLX is the convergence criterium in the
! x values
!
 implicit none

 real (kind=Rkind), pointer :: p(:) => null()
 real (kind=Rkind), intent(in) ::  gtol, tolx
 integer, intent(in) :: itmax
 integer :: n
 logical :: check
 real (kind=Rkind), parameter :: EPS=epsilon(p), STPMX=100
 real (kind=Rkind) :: xxxg, den, fp, temp, sum, stpmax, test, fret, fac, fae, fad, sumdg, sumxi
 integer :: its, iun,i,j
 real (kind=Rkind), allocatable :: g(:),dg(:),pnew(:),hdg(:),hessin(:,:),xi(:)

 real (kind=Rkind), target :: Qact(:)
 TYPE(Type_dnS)    :: MatdnE(:,:)
 TYPE (CoordType)  :: mole
 TYPE (PrimOp)     :: PrimOp
 TYPE (Tnum)       :: para_Tnum
 TYPE (param_BFGS) :: para_BFGS


 p => Qact(1:mole%nb_act)
 n =  mole%nb_act
 !write(out_unitp,*) 'n',n
 !write(out_unitp,*) 'gtol,tolx,itmax',gtol,tolx,itmax
 allocate (g(n),hdg(n),pnew(n),dg(n),xi(n),hessin(n,n))

 call dfunc(p,g,fp,MatdnE,mole,PrimOp,para_Tnum,1) ! Calculate starting function value and gradient.
 write(out_unitp,*)
 write(out_unitp,*) 'Iteration=    0 '
 write(out_unitp,*) ' Geometry:'
 write(out_unitp,*) ' p = ', p
 write(out_unitp,*) ' Energy = ',fp

 call proescvec(g,g,xxxg,n)
 xxxg=sqrt(xxxg/n)
!!!!!!!!!!!!!!
  test=ZERO            ! Test of convergence on zero gradient
  den=max(fp,ONE)
  do i=1,n
   temp=abs(g(i))*max(abs(p(i)),ONE)/den
   if (temp > test) test=temp
  end do
!!!!!!!!!!!!!!
 write(out_unitp,*) ' RMS Gradient = ',xxxg
 write(out_unitp,*) ' Test on gradient convergence = ', test
 call flush(out_unitp)

 IF (associated(para_BFGS%hessian_inv_init)) THEN
   write(out_unitp,*) ' The initial hessian is transfered'
   hessin(:,:) = para_BFGS%hessian_inv_init(:,:)
 ELSE
   hessin(:,:) = ZERO
   do i=1,n
     hessin(i,i)=ONE    !!!!!!!!!! Here Eckart has 1.0d0/a20(i)
   end do
 END IF

 xi=-g
 call proescvec(p,p,sum,n)
 stpmax=STPMX*max(sqrt(sum),real(n,kind=Rkind))
 do its=1,ITMAX                 ! Main loop over the iterations.
  call lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,check,MatdnE,mole,PrimOp,para_Tnum)
  ! The new function evaluation occurs in lnsrch; save the value in fp for the
  ! next line search. It is usually safe to ignore the value of check.
  fp=fret
  xi=pnew-p    ! update the line direction
  p=pnew       ! and the current point
  test=ZERO
  do i=1,n     ! Test the convergence in Delta(x)
   temp=abs(xi(i))/max(abs(p(i)),ONE)
   if (temp > test) test=temp
  end do
  dg=g                  ! Save the old gradient
  call dfunc(p,g,fret,MatdnE,mole,PrimOp,para_Tnum,1)  ! and get a new gradient

  test=ZERO            ! Test of convergence on zero gradient
  den=max(fret,ONE)
  do i=1,n
   temp=abs(g(i))*max(abs(p(i)),ONE)/den
   if (temp > test) test=temp
  end do
  call proescvec(g,g,xxxg,n)
  xxxg=sqrt(xxxg/n)
  IF (print_level > 0 .OR. test < tolx) THEN
    write(out_unitp,*)
    write(out_unitp,'(a,i4)') ' Iteration= ',its
    write(out_unitp,*) ' Geometry: '
    write(out_unitp,*) ' p = ', p
    write(out_unitp,*) ' Energy',fret
    write(out_unitp,*) ' RMS Gradient',xxxg
    write(out_unitp,*) ' Test on gradient convergence = ', test
    call flush(out_unitp)
  END IF
   if (test < tolx) then  !!! Testing what happen if this part is removed
   write(out_unitp,*) ' Geometry coordinates converged !! RMS step criteria'
   call flush(out_unitp)
   deallocate (g,hdg,pnew,dg,xi,hessin)
   return
  end if
!
  if ( test < gtol ) then
   write(out_unitp,*) ' Gradient converged !!'
   write(out_unitp,*) ' Optimized energy: ',fret
   write(out_unitp,*) ' Optimized RMS gradient: ', xxxg, test
   write(out_unitp,*) ' Optimized geometry: '
   write(out_unitp,*) ' p = ', p
   deallocate (g,hdg,pnew,dg,xi,hessin)
   return
  end if
!
  dg=g-dg          ! Compute diference of gradient
  do i=1,n         ! and difference times current matrix.
   call proescvec(hessin(:,i),dg,hdg(i),n)
  end do
  call proescvec(dg,xi,fac,n)  ! Calculate dot products for the denominators.
  call proescvec(dg,hdg,fae,n)
  call proescvec(dg,dg,sumdg,n)
  call proescvec(xi,xi,sumxi,n)
  if (fac**2 > EPS*sumdg*sumxi) then ! Skip update if fac is not
   fac=ONE/fac                        ! sufficiently positive.
   fad=ONE/fae
   do i=1,n         ! The vector that makes BFGS different from DFP.
    dg(i)=fac*xi(i)-fad*hdg(i)
   end do
   do i=1,n         ! The BFGS updating formula:
    do j=1,n
     hessin(i,j)=hessin(i,j)+fac*xi(i)*xi(j)-fad*hdg(i)*hdg(j)+fae*dg(i)*dg(j)
    end do
   end do
  end if
  do i=1,n          ! Now calculate the next direction to go,
   call proescvec(hessin(:,i),g,xi(i),n)
  end do
  xi=-xi
 end do             ! and go back for another iteration.
 deallocate (g,hdg,pnew,dg,xi,hessin)
 write(out_unitp,*) 'Too many iterations in dfpmin'
 write(out_unitp,*) ' energy: ',fret
 write(out_unitp,*) ' Geometry: '
 write(out_unitp,*) ' p = ', p
 stop
 return
 end subroutine dfpmin_new

!----------------------------------------------------------------------
  subroutine lnsrch(n,xold,fold,g,p,x,f,stpmax,check,MatdnE,mole,PrimOp,para_Tnum)
!---------------------------------------------------------------------
!
 USE mod_system
 USE mod_Coord_KEO
 USE mod_PrimOp
 USE mod_basis
 USE mod_Op
 USE mod_Auto_Basis
!
 implicit none
 LOGICAL :: check
 integer :: n,i
 real (kind=Rkind) :: g(n),p(n),x(n),xold(n)
 real (kind=Rkind), PARAMETER :: ALF=ONETENTH**4,TOLX=ONETENTH**9
 real (kind=Rkind) :: fold,f,stpmax,sum,slope,test,temp,alamin,alam,tmplam
 real (kind=Rkind) :: rhs1, rhs2, a, b, alam2, disc, f2
!
 TYPE(Type_dnS)   :: MatdnE(:,:)
 TYPE (CoordType) :: mole
 TYPE (PrimOp)    :: PrimOp
 TYPE (Tnum)      :: para_Tnum
!
 check=.false.
 call proescvec(p,p,sum,n)
 sum=sqrt(sum)
!
! write(out_unitp,*) 'sum=', sum, 'stpmax=', stpmax
!
 call flush(out_unitp)
 if(sum.gt.stpmax)then
  do i=1,n
   p(i)=p(i)*(stpmax/sum)
  end do
 endif
 call proescvec(g,p,slope,n)
 test=ZERO
 do i=1,n
  temp=abs(p(i))/max(abs(xold(i)),ONE)
  if(temp.gt.test)test=temp
 end do
 alamin=TOLX/test
 alam=ONE
1 continue
 do  i=1,n
  x(i)=xold(i)+alam*p(i)
 end do
 call dfunc(x,g,f,MatdnE,mole,PrimOp,para_Tnum,0)
 if(alam < alamin)then
  x=xold
  check=.true.
  return
 else if(f <= fold+ALF*alam*slope) then
  return
 else
  if(alam == ONE)then
    tmplam=-slope/(TWO*(f-fold-slope))
   else
    rhs1=f-fold-alam*slope
    rhs2=f2-fold-alam2*slope
    a=(rhs1/alam**2-rhs2/alam2**2)/(alam-alam2)
    b=(-alam2*rhs1/alam**2+alam*rhs2/alam2**2)/(alam-alam2)
    if(a == ZERO)then
     tmplam=-slope/(TWO*b)
    else
     disc=b*b-THREE*a*slope
     if (disc < ZERO) then
      tmplam=HALF*alam
     else if (b < ZERO) then
      tmplam=(-b+sqrt(disc))/(THREE*a)
     else
      tmplam=-slope/(b+sqrt(disc))
     end if
    endif
    if(tmplam > HALF*alam) tmplam=HALF*alam
   endif
  endif
  alam2=alam
  f2=f
  alam=max(tmplam,ONETENTH*alam)
  goto 1
  end subroutine
!  (C) Copr. 1986-92 Numerical Recipes Software Bc21D#,#5,15!".
!
!---------------------------------------------------------------------------
  SUBROUTINE dfunc(xt,df,f,MatdnE,mole,PrimOp,para_Tnum,nderiv_dnE)
!---------------------------------------------------------------------------
 USE mod_system
 USE mod_Coord_KEO
 USE mod_PrimOp
 USE mod_basis
 USE mod_Op
 USE mod_Auto_Basis
!
! calculation of the gradient at xt
!
 IMPLICIT none
 integer, intent(in) :: nderiv_dnE
 TYPE(Type_dnS)      :: MatdnE(:,:)
 TYPE (CoordType)    :: mole
 TYPE (PrimOp)       :: PrimOp
 TYPE (Tnum)         :: para_Tnum

 integer :: i
 real (kind=Rkind)   :: Qact(mole%nb_var)

 real(kind=Rkind), intent(in) :: xt(mole%nb_act)

 real(kind=Rkind),intent(inout) :: df(mole%nb_act), f
!
! write(out_unitp,*) 'dfunc subroutine',mole%nb_act,nderiv_dnE
! write(out_unitp,*) 'xt = ', xt

 Qact(:) = ZERO
 Qact(1:mole%nb_act)=xt

! write(out_unitp,*) 'Qact = ',Qact
! flush(out_unitp)
!
! The subroutine below enables to calculate the energy, the gradient and/or the hessian
! nderiv_dnE = 0 : => energy only
! nderiv_dnE = 1 : => energy and gradient (here nderiv_dnE = 1)
! nderiv_dnE = 2 : => energy, gradient and hessian
! the derivatives are done with respect to the number of active coordinates

 MatdnE(1,1)%d0    = ZERO
 IF (nderiv_dnE > 0) MatdnE(1,1)%d1(:) = ZERO


 CALL dnOp_grid(Qact,MatdnE,nderiv_dnE,mole,para_Tnum,PrimOp)

! the results are in the matrix of derived type MatdnE(:,:)
!   Remark: we have a matrix because ElVibRot can deal with several diabatic electronic states
!
!  For one electronic state (which is the only possibility with on-the-fly calculation)
!    the energy   is in : MatdnE(1,1)%d0
!    the gradient is in : MatdnE(1,1)%d1(:)     (if nderiv_dnE >= 1)
!    the hessian  is in : MatdnE(1,1)%d2(:,:)   (if nderiv_dnE >= 2)

!  Remark, the allocation of MatdnE have to be done with "nderiv_alloc" larger than "nderiv_dnE"
!     With nderiv_alloc=2 and nderiv_dnE=1, you can calculate the energy and the gradient
!     and update the hessian in MatdnE(1,1)%d2(:,:), if you want.

!  write(out_unitp,*) ' MatdnE(1,1)%d0 = ', MatdnE(1,1)%d0
!  write(out_unitp,*) ' MatdnE(1,1)%d1 = ', MatdnE(1,1)%d1
!
  f=MatdnE(1,1)%d0
  if (nderiv_dnE >= 1) df=MatdnE(1,1)%d1(1:mole%nb_act)

!  write(out_unitp,*) ' Energy = ',f
!  write(out_unitp,*) ' Active modes gradient = ', df
!  write(out_unitp,*) 'end dfunc subroutine'
!  flush(out_unitp)

 END SUBROUTINE dfunc
! ------------------------------------------------------------------------
      subroutine proescvec(vec1,vec2,ps,ndim)
! ------------------------------------------------------------------------
!     producte escalar vec1*vec2
!
      implicit none
      integer :: ndim, i
      real(kind=Rkind) :: vec1(ndim),vec2(ndim)
      real(kind=Rkind) :: ps
!
      intent (in) ndim,vec1,vec2
      intent (inout) ps
!
      ps=ZERO
      do i=1,ndim
       ps=ps+vec1(i)*vec2(i)
      end do
!
      return
      end subroutine
!

END MODULE mod_BFGS
