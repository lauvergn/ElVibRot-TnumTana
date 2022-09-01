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
MODULE mod_Filter
USE mod_system
USE mod_Constant
USE mod_psi, ONLY : param_psi,alloc_psi,dealloc_psi
IMPLICIT NONE
TYPE param_filter

   integer :: filter_type = 1    ! rectangular in [A:B]

   real (kind=Rkind) :: A=-ONE   ! for type 1,3
   real (kind=Rkind) :: B= ONE   ! for type 1,3

   real (kind=Rkind) :: beta = TEN  ! for type 3, smooth rectangular in [A:B]

   real (kind=Rkind) :: l = 1    ! for type 4, sine in [A:B]


   real (kind=Rkind) :: sigma = ONETENTH ! for type 2 (gaussian)
   real (kind=Rkind) :: E0    = ZERO     ! for type 2 (gaussian)

END TYPE param_filter

PRIVATE
PUBLIC :: sub_FilterDiagonalization,sub_GaussianFilterDiagonalization,sub_GaussianFilterDiagonalization_v0

CONTAINS

      SUBROUTINE sub_GaussianFilterDiagonalization(psi,Ene,nb_diago,max_diago, &
                                                   para_H,para_propa)
      USE mod_system
      USE mod_psi,    ONLY : param_psi,alloc_psi,Set_Random_psi,        &
                             Set_symab_OF_psiBasisRep,renorm_psi,       &
                             Overlap_psi1_psi2,norm2_psi,dealloc_psi

      USE mod_Op
      USE mod_propa
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ----------------------------
      integer                      :: nb_diago,max_diago
      TYPE (param_propa)           :: para_propa
      TYPE (param_psi)             :: psi(max_diago)
      real (kind=Rkind)            :: Ene(max_diago)
!------ working parameters --------------------------------

      real (kind=Rkind), allocatable :: x_cheby(:),w_cheby(:),d0P_cheby(:,:)
      real (kind=Rkind), allocatable :: sqw_cheby(:),P0_cheby(:),P1_cheby(:),P2_cheby(:)

      integer :: nb,nq,nb0,nstep

      real (kind=Rkind), allocatable :: wfl(:),wfl_cheby(:)


      real (kind=Rkind) :: phi_j,Delta_Lambda,Lambda

      integer       :: i,n,jorth,jsave

      TYPE (param_psi)             :: q1 ! intial vector (random ?)
      TYPE (param_psi)             :: z(para_propa%para_Davidson%Lmax_filter) ! intial vector (random ?)
      TYPE (param_psi)             :: Hz(para_propa%para_Davidson%Lmax_filter) ! intial vector (random ?)

      TYPE (param_psi)             :: Tnq1(para_propa%para_Davidson%Mmax_filter+1) ! intial vector (random ?)
      TYPE (param_psi)             :: g,w1 ! working vector

      real (kind=Rkind), allocatable :: f(:,:)
      real (kind=Rkind) :: acEj,acDeltaj,sigma,ff_filter,El,A_filter,B_filter,filter_err
      TYPE (param_filter) :: filter
      integer :: j,k,l,JJ


      real (kind=Rkind), allocatable :: H(:,:),Vec(:,:)

      logical                   :: Conv,convergeEne(max_diago),convergeResi(max_diago)
      real (kind=Rkind)         :: non_hermitic,epsi,auTOcm_inv
      real (kind=Rkind)         :: RS,a,max_Sii,max_Sij
      real (kind=Rkind)         :: tab_norm2g(max_diago)
      complex (kind=Rkind)      :: Overlap
      integer :: m0,mf,nb_Vec_IN_Window,nb_ConvVec_IN_Window,DeltaL
      real (kind=Rkind)         :: filter_err_thresh = ONETENTH**8

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_FilterDiagonalization'
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      para_propa%para_Davidson%E0_filter =                              &
                   para_propa%para_Davidson%E0_filter + para_propa%Hmin
      nb_diago = para_propa%para_Davidson%L_filter
      DeltaL   = para_propa%para_Davidson%L_filter
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Hmin,Hmax     : ',para_propa%Hmin,para_propa%Hmax
        write(out_unitp,*) 'E0_filter (ua): ',para_propa%para_Davidson%E0_filter
        write(out_unitp,*) 'L_filter      : ',para_propa%para_Davidson%L_filter
        write(out_unitp,*) 'Lmax_filter   : ',para_propa%para_Davidson%Lmax_filter
        write(out_unitp,*) 'M_filter      : ',para_propa%para_Davidson%M_filter
        write(out_unitp,*) 'Mmax_filter   : ',para_propa%para_Davidson%Mmax_filter
        !write(out_unitp,*) 'nb_diago      : ',nb_diago
        write(out_unitp,*) 'max_diago     : ',max_diago

      END IF
!-----------------------------------------------------------
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')
!-----------------------------------------------------------

      write(out_unitp,*) ' Propagation: ',para_propa%name_WPpropa
      CALL Set_ZPE_OF_Op(para_H,ZPE=para_propa%Hmin,forced=.TRUE.)
      write(out_unitp,*) 'ZPE (cm-1)',para_H%ZPE * auTOcm_inv

      ! change Hmin and Hmax to be sure that the spectral range is between Hmin and Hmax.
      para_propa%Hmin = para_propa%Hmin - ONETENTH**2 * (para_propa%Hmax - para_propa%Hmin)
      para_propa%Hmax = para_propa%Hmax + ONETENTH**2 * (para_propa%Hmax - para_propa%Hmin)

!     - scaling of H ---------------------------------------
      para_propa%para_poly%deltaE = para_propa%Hmax - para_propa%Hmin
      para_propa%para_poly%E0     = para_propa%Hmin + HALF * para_propa%para_poly%deltaE
      para_propa%para_poly%Esc    = HALF * para_propa%para_poly%deltaE

      para_H%scaled = .TRUE.
      para_H%E0     = para_propa%para_poly%E0
      para_H%Esc    = para_propa%para_poly%Esc

      write(out_unitp,*) ' deltaE,E0,Esc: ',para_propa%para_poly%deltaE,&
                                            para_H%E0,para_H%Esc
!-----------------------------------------------------------

     write(out_unitp,*) 'W_filter (ua)  : ',para_propa%para_Davidson%W_filter
     write(out_unitp,*) 'W_filter (cm-1): ',para_propa%para_Davidson%W_filter*&
                                                                   auTOcm_inv


     Delta_Lambda = para_propa%para_Davidson%W_filter /                 &
                    real(para_propa%para_Davidson%L_filter-1,kind=Rkind)

!     Delta_Lambda = para_propa%para_Davidson%W_filter *TWO /            &
!                    real(para_propa%para_Davidson%L_filter-1,kind=Rkind)

     write(out_unitp,*) ' Delta_Lambda (ua)  : ',Delta_Lambda
     write(out_unitp,*) ' Delta_Lambda (cm-1): ',Delta_Lambda*auTOcm_inv

     para_propa%para_Davidson%LambdaMin = max(para_H%ZPE,         &
         para_propa%para_Davidson%E0_filter - HALF*para_propa%para_Davidson%W_filter)
     para_propa%para_Davidson%LambdaMax =                               &
         para_propa%para_Davidson%LambdaMin + para_propa%para_Davidson%W_filter
     write(out_unitp,*) ' LambdaMin,LambdaMax (cm-1): ',                &
                        para_propa%para_Davidson%LambdaMin * auTOcm_inv,&
                        para_propa%para_Davidson%LambdaMax * auTOcm_inv
     CALL flush_perso(out_unitp)

!-----------------------------------------------------------
! define the number of chebychev polynomials automatically
!-----------------------------------------------------------
sigma    = para_propa%para_Davidson%W_filter /para_H%Esc / para_propa%para_Davidson%L_filter
El       = (para_propa%para_Davidson%E0_filter-para_H%E0)/para_H%Esc ! to scale the energy between [-1,1]
A_filter = El-Delta_Lambda/para_H%Esc/TWO
B_filter = El+Delta_Lambda/para_H%Esc/TWO
nq       = pi / abs(acos(B_filter)-acos(A_filter))/3
nb0      = nq
nb       = nb0
write(out_unitp,*) ' nb,nq init: ',nb,nq
CALL flush_perso(out_unitp)

CALL alloc_NParray(x_cheby,[nq],'x_cheby',name_sub)
CALL alloc_NParray(w_cheby,[nq],'w_cheby',name_sub)
CALL alloc_NParray(sqw_cheby,[nq],'sqw_cheby',name_sub)
CALL alloc_NParray(P0_cheby,[nq],'P0_cheby',name_sub)
CALL alloc_NParray(P1_cheby,[nq],'P1_cheby',name_sub)
CALL alloc_NParray(P2_cheby,[nq],'P2_cheby',name_sub)

DO
  nb = nb*11/10
  CALL dealloc_NParray(x_cheby,'x_cheby',name_sub)
  CALL dealloc_NParray(w_cheby,'w_cheby',name_sub)
  CALL dealloc_NParray(sqw_cheby,'sqw_cheby',name_sub)

  CALL dealloc_NParray(P0_cheby,'P0_cheby',name_sub)
  CALL dealloc_NParray(P1_cheby,'P1_cheby',name_sub)
  CALL dealloc_NParray(P2_cheby,'P2_cheby',name_sub)

  ! first the chebychev polynomials + grid/weight
  !nb = para_propa%para_Davidson%M_filter+1
  nq = nb
  CALL alloc_NParray(x_cheby,[nq],'x_cheby',name_sub)
  CALL alloc_NParray(w_cheby,[nq],'w_cheby',name_sub)
  CALL alloc_NParray(sqw_cheby,[nq],'sqw_cheby',name_sub)
  !CALL alloc_NParray(d0P_cheby,[nq,nb],'d0P_cheby',name_sub)
  CALL alloc_NParray(P0_cheby,[nq],'P0_cheby',name_sub)
  CALL alloc_NParray(P1_cheby,[nq],'P1_cheby',name_sub)
  CALL alloc_NParray(P2_cheby,[nq],'P2_cheby',name_sub)

  CALL gauss_chebyWeight(x_cheby,w_cheby,nq)
  sqw_cheby = ONE/sqrt(sqrt(1-x_cheby*x_cheby)) * sqrt(TWO/pi)
  !write(out_unitp,*) 'x_cheby',x_cheby
  !write(out_unitp,*) 'w_cheby',w_cheby
  !CALL d0poly_chebyWeight_grid(x_cheby,d0P_cheby,nb,nq)
  !DO i=1,nb
  !DO j=i,nb
  !  write(out_unitp,*) 'Sij',i,j,dot_product(d0P_cheby(:,j),w_cheby*d0P_cheby(:,i))
  !  CALL flush_perso(out_unitp)
  !END DO
  !END DO
  !-----------------------------------------------------------
  ! then one filter (at window center)
  CALL Set_filter(filter,2,E0=El,sigma=sigma)
  !CALL Set_filter(filter,1,A=A_filter,B=B_filter)
  !CALL Set_filter(filter,3,A=A_filter,B=B_filter,beta=FOUR)
  !CALL Set_filter(filter,4,A=A_filter,B=B_filter)

  CALL alloc_NParray(f,[nb,1],'f',name_sub)
  CALL alloc_NParray(wfl,[nq],'wfl',name_sub)
  CALL alloc_NParray(wfl_cheby,[nq],'wfl',name_sub)
  wfl_cheby = ZERO

  ! l-filter values x weight on the grid
  DO k=1,nq
    wfl(k) = f_filter(x_cheby(k),filter) * w_cheby(k)
  END DO

  P0_cheby(:) = ONE
  j=1
  f(j,1) = dot_product(wfl,P0_cheby*sqw_cheby/sqrt(TWO))

  IF (debug) wfl_cheby(:) = f(j,1)*P0_cheby*sqw_cheby/sqrt(TWO)


  P1_cheby(:) = x_cheby
  j=2
  f(j,1) = dot_product(wfl,P1_cheby*sqw_cheby)*KDF(j,nb)
  IF (debug) wfl_cheby(:) = wfl_cheby(:) + f(j,1)*P1_cheby*sqw_cheby



  DO j=3,nb
    P2_cheby = TWO*x_cheby* P1_cheby - P0_cheby
    f(j,1) = dot_product(wfl,P2_cheby*sqw_cheby)*KDF(j,nb)
    IF (debug) wfl_cheby(:) = wfl_cheby(:) + f(j,1)*P2_cheby*sqw_cheby

    P0_cheby = P1_cheby
    P1_cheby = P2_cheby
  END DO

  DO k=1,nq
    write(66,*) x_cheby(k),wfl(k),wfl_cheby(k),log10(abs(wfl(k)-wfl_cheby(k)))
  END DO

  filter_err = sum(abs(f(nb-3:nb,1)))/THREE

  CALL dealloc_NParray(f,'f',name_sub)
  CALL dealloc_NParray(wfl,'wfl',name_sub)
  CALL dealloc_NParray(wfl_cheby,'wfl',name_sub)



  write(out_unitp,*) 'nb,filter_err',nb,filter_err
  CALL flush_perso(out_unitp)
  IF (filter_err < filter_err_thresh .OR. nb > para_propa%para_Davidson%M_filter) EXIT



END DO
nb = size(P0_cheby)

  !-----------------------------------------------------------
! then the filter(s) are projected on the chebyshev polynomials
  CALL alloc_NParray(f,[nb,para_propa%para_Davidson%L_filter],'f',name_sub)
  CALL alloc_NParray(wfl,[nq],'wfl',name_sub)
  sigma       = para_propa%para_Davidson%W_filter /para_H%Esc / para_propa%para_Davidson%L_filter


  DO l=1,para_propa%para_Davidson%L_filter
    El = para_propa%para_Davidson%LambdaMin + real(l-1,kind=Rkind)*Delta_Lambda

    El = (El-para_H%E0)/para_H%Esc ! to scale the energy between [-1,1]
    A_filter  = El-Delta_Lambda/para_H%Esc/TWO
    B_filter  = El+Delta_Lambda/para_H%Esc/TWO

    CALL Set_filter(filter,2,E0=El,sigma=sigma)
    !CALL Set_filter(filter,1,A=A_filter,B=B_filter)
    !CALL Set_filter(filter,3,A=A_filter,B=B_filter,beta=FOUR)
    !CALL Set_filter(filter,4,A=A_filter,B=B_filter)

    IF (debug) write(out_unitp,*) 'l,El',l,El
    CALL flush_perso(out_unitp)


    ! l-filter values x weight on the grid
    DO k=1,nq
      wfl(k) = f_filter(x_cheby(k),filter) * w_cheby(k)
    END DO


    P0_cheby(:) = ONE
    j=1
    f(j,l) = dot_product(wfl,P0_cheby*sqw_cheby/sqrt(TWO))

    P1_cheby(:) = x_cheby
    j=2
    f(j,l) = dot_product(wfl,P1_cheby*sqw_cheby)*KDF(j,nb)

    DO j=3,nb
      P2_cheby = TWO*x_cheby* P1_cheby - P0_cheby
      f(j,l) = dot_product(wfl,P2_cheby*sqw_cheby)*KDF(j,nb)
      P0_cheby = P1_cheby
      P1_cheby = P2_cheby
    END DO

    IF (debug) write(out_unitp,*) 'l,end err of f',l,sum(abs(f(nb-10:nb,l)))/TEN
    CALL flush_perso(out_unitp)
    !write(out_unitp,*) '========================================='
    !write(out_unitp,*) 'f(:,l)',l
    !write(out_unitp,'(10f10.6)') ,f(:,l)

  END DO
  CALL dealloc_NParray(x_cheby,'x_cheby',name_sub)
  CALL dealloc_NParray(w_cheby,'w_cheby',name_sub)
  CALL dealloc_NParray(sqw_cheby,'sqw_cheby',name_sub)

  CALL dealloc_NParray(P0_cheby,'P0_cheby',name_sub)
  CALL dealloc_NParray(P1_cheby,'P1_cheby',name_sub)
  CALL dealloc_NParray(P2_cheby,'P2_cheby',name_sub)

  IF (debug) write(out_unitp,*) 'chebychev coef: done',nb
  CALL flush_perso(out_unitp)
!-----------------------------------------------------------

!- vector initialization:  q1 -----------
  CALL init_psi(q1,para_H,para_H%cplx)
  CALL init_psi(g,para_H,para_H%cplx)

  CALL alloc_psi(q1)
  CALL Set_Random_psi(q1)
  CALL Set_symab_OF_psiBasisRep(q1,para_propa%para_Davidson%symab)
  CALL renorm_psi(q1,BasisRep=.TRUE.)
!- vector initialization:  q1 -----------

!----------------------------------------------------------

     !Conv = .FALSE.
     !DO

       !- chebychev recursion -------------------------------------------
       CALL sub_chebychev_recursion(Tnq1,q1,0,nb-1,para_H)
       IF (debug) write(out_unitp,*) 'chebychev recursion: done',nb
       CALL flush_perso(out_unitp)

       !- Coefficient of the chebychev expansion (with the filter)  -----

       !- z vectors --------------------------------------------
       CALL sub_newZ_vectors_withf(z,Tnq1,nb,para_H,para_propa,f)
       nb_diago = count(abs(z(:)%CAvOp)>ONETENTH**9)
       IF (debug) write(out_unitp,*) 'ortho z vectors: done',nb_diago
       CALL flush_perso(out_unitp)
       !- z vectors -------------------------------------------

        !- diagonalization -----------------------------------
        CALL alloc_NParray(H,  [nb_diago,nb_diago],'H',  name_sub)
        CALL alloc_NParray(Vec,[nb_diago,nb_diago],'Vec',name_sub)

        DO j=1,nb_diago
        DO i=1,nb_diago
          CALL Overlap_psi1_psi2(Overlap,z(i),z(j))
          H(i,j) = real(Overlap,kind=Rkind)
        END DO
        END DO
        !IF (debug) write(out_unitp,*) 'S matrix: done'
        CALL flush_perso(out_unitp)
        CALL sub_ana_S(H,nb_diago,max_Sii,max_Sij,.TRUE.)
        !IF (debug) CALL Write_Mat(H,out_unitp,5)



        DO j=1,nb_diago
          !H.z(j) calc
          CALL sub_OpPsi(z(j),Hz(j),para_H)                            ! => H.z(j)

          DO i=1,nb_diago
            CALL Overlap_psi1_psi2(Overlap,z(i),Hz(j))
            H(i,j) = real(Overlap,kind=Rkind)
          END DO
        END DO
        !IF (debug) write(out_unitp,*) 'H matrix: done'
        CALL flush_perso(out_unitp)

        CALL sub_hermitic_H(H,nb_diago,non_hermitic,para_H%sym_Hamil)
        !IF (debug) CALL Write_Mat(H,out_unitp,5)

        IF (non_hermitic > FOUR*ONETENTH**4) THEN
          If(MPI_id==0) write(out_unitp,*) 'WARNING: non_hermitic is BIG'
          If(MPI_id==0) write(out_unitp,31) non_hermitic
 31       format(' Hamiltonien: ',f16.12,' au')
        ELSE
          If(MPI_id==0) write(out_unitp,51) non_hermitic*auTOcm_inv
 51       format(' Hamiltonien: ',f16.12,' cm-1')
        END IF
        epsi = max(para_propa%para_Davidson%conv_resi,                    &
                   TEN**para_propa%para_Davidson%conv_hermitian *       &
                                               non_hermitic)

        IF (para_H%sym_Hamil) THEN
          CALL diagonalization(H,Ene(1:nb_diago),Vec,nb_diago,3,1,.FALSE.)
        ELSE
          CALL diagonalization(H,Ene(1:nb_diago),Vec,nb_diago,4,1,.FALSE.)
        END IF

         write(out_unitp,21) Ene(1:nb_diago)*auTOcm_inv
         write(out_unitp,21) (Ene(1:nb_diago)-para_H%ZPE)*auTOcm_inv
 21      format(' Filter: ',50(1x,f18.4))


        !----------------------------------------------------------
        !- residual vector ---------------------------
        IF (debug) write(out_unitp,*) 'residual'
        CALL flush_perso(out_unitp)
        DO j=1,nb_diago
          g = ZERO
          DO i=1,nb_diago
            w1 = Hz(i) - z(i) * Ene(j)
            g = g + w1 * Vec(i,j)
          END DO
          CALL norm2_psi(g)
          tab_norm2g(j) = sqrt(g%norm2)
          convergeResi(j) = tab_norm2g(j) < epsi

        END DO
        write(out_unitp,41) 'tab_norm2g(:):      ',tab_norm2g(1:nb_diago)
        write(out_unitp,42) 'convergenceResi(:): ',convergeResi(1:nb_diago)
 41     format(a,100(1x,e9.2))
 42     format(a,100(1x,l9))
        IF (debug) write(out_unitp,*) 'residual: done'

        nb_Vec_IN_Window     = 0
        nb_ConvVec_IN_Window = 0
        DO j=1,nb_diago
          IF (Ene(j) >= para_propa%para_Davidson%LambdaMin .AND.         &
              Ene(j) <= para_propa%para_Davidson%LambdaMax) THEN

            nb_Vec_IN_Window = nb_Vec_IN_Window + 1
            IF (convergeResi(j)) nb_ConvVec_IN_Window = nb_ConvVec_IN_Window + 1
          END IF
        END DO
        Conv = (nb_ConvVec_IN_Window == nb_Vec_IN_Window) .AND. nb_ConvVec_IN_Window > 0

        write(out_unitp,*)  ' Converged levels:               ',count(convergeResi(1:nb_diago))
        write(out_unitp,*)  ' Converged levels in the windox: ',nb_ConvVec_IN_Window
        write(out_unitp,*)  ' Levels in the window:           ',nb_Vec_IN_Window
        write(out_unitp,*)  ' Convergence ?:                  ',Conv
        CALL flush_perso(out_unitp)
        !- residual vector and convergence ------------------------
        !----------------------------------------------------------
        !IF (Conv .OR. mf >= para_propa%para_Davidson%Mmax_filter) EXIT

        CALL dealloc_NParray(H,  'H',  name_sub)
        CALL dealloc_NParray(Vec,'Vec',name_sub)

STOP
        m0 = mf + 1
        mf = min(mf + para_propa%para_Davidson%DeltaM_filter,           &
                                   para_propa%para_Davidson%Mmax_filter)

        IF (nb_diago == para_propa%para_Davidson%L_filter) THEN
          para_propa%para_Davidson%L_filter = para_propa%para_Davidson%L_filter + DeltaL
        END IF
        IF (para_propa%para_Davidson%L_filter > para_propa%para_Davidson%Lmax_filter) THEN
          para_propa%para_Davidson%L_filter = para_propa%para_Davidson%Lmax_filter
        END IF
        IF (para_propa%para_Davidson%L_filter > max_diago) THEN
          para_propa%para_Davidson%L_filter = max_diago
        END IF
        STOP
      !END DO

      !----------------------------------------------------------
      !! save converged vectors and vectors in the window
      jsave = 0
      DO j=1,nb_diago
        IF (convergeResi(j) .OR.                                        &
                    (Ene(j) >= para_propa%para_Davidson%LambdaMin .AND. &
                     Ene(j) <= para_propa%para_Davidson%LambdaMax)) THEN
          jsave = jsave + 1

          CALL init_psi(psi(jsave),para_H,para_H%cplx)
          psi(jsave) = ZERO
          DO i=1,nb_diago
            psi(jsave) = psi(jsave) + Vec(i,j) * z(i)
          END DO
          psi(jsave)%CAvOp    = Ene(j)
          psi(jsave)%IndAvOp  = para_H%n_Op  ! it should be 0
          psi(jsave)%convAvOp = convergeResi(j)
        END IF
      END DO
      nb_diago = jsave
      IF (nb_diago < 1) THEN
        write(out_unitp,*) 'WARNING in ',name_sub
        write(out_unitp,*) ' There is no vector in the filter range!!'
        !STOP 'filter diago'
      END IF

     CALL dealloc_psi(w1)
     CALL dealloc_psi(g)
     CALL dealloc_psi(q1)
     DO i=1,size(z)
       CALL dealloc_psi(z(i))
     END DO
     DO i=1,size(Hz)
       CALL dealloc_psi(Hz(i))
     END DO
     DO i=lbound(Tnq1,dim=1),ubound(Tnq1,dim=1)
       CALL dealloc_psi(Tnq1(i))
     END DO
     CALL dealloc_NParray(H,  'H',  name_sub)
     CALL dealloc_NParray(Vec,'Vec',name_sub)

!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------


      END SUBROUTINE sub_GaussianFilterDiagonalization

      ! ONE filter, but several initial vector
      SUBROUTINE sub_BlockFilterDiagonalization(psi,Ene,nb_diago,max_diago, &
                                                   para_H,para_propa)
      USE mod_system
      USE mod_psi,    ONLY : param_psi,alloc_psi,Set_Random_psi,        &
                             Set_symab_OF_psiBasisRep,renorm_psi,       &
                             Overlap_psi1_psi2,norm2_psi,dealloc_psi
      USE mod_Op
      USE mod_propa
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ----------------------------
      integer                      :: nb_diago,max_diago
      TYPE (param_propa)           :: para_propa
      TYPE (param_psi)             :: psi(max_diago)
      real (kind=Rkind)            :: Ene(max_diago)
!------ working parameters --------------------------------

      real (kind=Rkind), allocatable :: x_cheby(:),w_cheby(:),d0P_cheby(:,:)
      integer :: nb,nq

      real (kind=Rkind), allocatable :: wfl(:)


      real (kind=Rkind) :: phi_j,Delta_Lambda,Lambda

      integer       :: i,n,jorth,jsave

      TYPE (param_psi)             :: q1 ! intial vector (random ?)
      TYPE (param_psi)             :: z(para_propa%para_Davidson%Lmax_filter) ! intial vector (random ?)
      TYPE (param_psi)             :: Hz(para_propa%para_Davidson%Lmax_filter) ! intial vector (random ?)

      TYPE (param_psi)             :: Tnq1(para_propa%para_Davidson%Mmax_filter+1) ! intial vector (random ?)
      TYPE (param_psi)             :: g,w1 ! working vector

      real (kind=Rkind), allocatable :: f(:,:)
      real (kind=Rkind) :: acEj,acDeltaj,sigma,ff_filter,El
      integer :: j,k,l,JJ


      real (kind=Rkind), allocatable :: H(:,:),Vec(:,:)

      logical                   :: Conv,convergeEne(max_diago),convergeResi(max_diago)
      real (kind=Rkind)         :: non_hermitic,epsi,auTOcm_inv
      real (kind=Rkind)         :: RS,a,max_Sii,max_Sij
      real (kind=Rkind)         :: tab_norm2g(max_diago)
      complex (kind=Rkind)      :: Overlap
      integer :: m0,mf,nb_Vec_IN_Window,nb_ConvVec_IN_Window,DeltaL

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_BlockFilterDiagonalization'
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      para_propa%para_Davidson%E0_filter =                              &
                   para_propa%para_Davidson%E0_filter + para_propa%Hmin
      nb_diago = para_propa%para_Davidson%L_filter
      DeltaL   = para_propa%para_Davidson%L_filter
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Hmin,Hmax     : ',para_propa%Hmin,para_propa%Hmax
        write(out_unitp,*) 'E0_filter (ua): ',para_propa%para_Davidson%E0_filter
        write(out_unitp,*) 'L_filter      : ',para_propa%para_Davidson%L_filter
        write(out_unitp,*) 'Lmax_filter   : ',para_propa%para_Davidson%Lmax_filter
        write(out_unitp,*) 'M_filter      : ',para_propa%para_Davidson%M_filter
        write(out_unitp,*) 'Mmax_filter   : ',para_propa%para_Davidson%Mmax_filter
        !write(out_unitp,*) 'nb_diago      : ',nb_diago
        write(out_unitp,*) 'max_diago     : ',max_diago

      END IF
!-----------------------------------------------------------
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')
!-----------------------------------------------------------

      write(out_unitp,*) ' Propagation: ',para_propa%name_WPpropa
      CALL Set_ZPE_OF_Op(para_H,ZPE=para_propa%Hmin,forced=.TRUE.)
      write(out_unitp,*) 'ZPE (cm-1)',para_H%ZPE * auTOcm_inv

      ! change Hmin and Hmax to be sure that the spectral range is between Hmin and Hmax.
      para_propa%Hmin = para_propa%Hmin - ONETENTH**2 * (para_propa%Hmax - para_propa%Hmin)
      para_propa%Hmax = para_propa%Hmax + ONETENTH**2 * (para_propa%Hmax - para_propa%Hmin)

!     - scaling of H ---------------------------------------
      para_propa%para_poly%deltaE = para_propa%Hmax - para_propa%Hmin
      para_propa%para_poly%E0     = para_propa%Hmin + HALF * para_propa%para_poly%deltaE
      para_propa%para_poly%Esc    = HALF * para_propa%para_poly%deltaE

      para_H%scaled = .TRUE.
      para_H%E0     = para_propa%para_poly%E0
      para_H%Esc    = para_propa%para_poly%Esc

      write(out_unitp,*) ' deltaE,E0,Esc: ',para_propa%para_poly%deltaE,&
                                            para_H%E0,para_H%Esc
!-----------------------------------------------------------

     write(out_unitp,*) 'W_filter (ua)  : ',para_propa%para_Davidson%W_filter
     write(out_unitp,*) 'W_filter (cm-1): ',para_propa%para_Davidson%W_filter*&
                                                                   auTOcm_inv


     Delta_Lambda = para_propa%para_Davidson%W_filter /                 &
                    real(para_propa%para_Davidson%L_filter-1,kind=Rkind)
     write(out_unitp,*) ' Delta_Lambda (ua)  : ',Delta_Lambda
     write(out_unitp,*) ' Delta_Lambda (cm-1): ',Delta_Lambda*auTOcm_inv

     para_propa%para_Davidson%LambdaMin = max(para_H%ZPE,         &
         para_propa%para_Davidson%E0_filter - HALF*para_propa%para_Davidson%W_filter)
     para_propa%para_Davidson%LambdaMax =                               &
         para_propa%para_Davidson%LambdaMin + para_propa%para_Davidson%W_filter
     write(out_unitp,*) ' LambdaMin,LambdaMax (cm-1): ',                &
                        para_propa%para_Davidson%LambdaMin * auTOcm_inv,&
                        para_propa%para_Davidson%LambdaMax * auTOcm_inv

!-----------------------------------------------------------
! first the chebychev polynomials + grid/weight
  nb = para_propa%para_Davidson%M_filter+1
  nq = para_propa%para_Davidson%M_filter+1
  CALL alloc_NParray(x_cheby,[nq],'x_cheby',name_sub)
  CALL alloc_NParray(w_cheby,[nq],'w_cheby',name_sub)
  CALL alloc_NParray(d0P_cheby,[nq,nb],'d0P_cheby',name_sub)
  CALL gauss_chebyWeight(x_cheby,w_cheby,nq)
  !write(out_unitp,*) 'x_cheby',x_cheby
  !write(out_unitp,*) 'w_cheby',w_cheby
  CALL d0poly_chebyWeight_grid(x_cheby,d0P_cheby,nb,nq)
  !DO i=1,nb
  !DO j=i,nb
  !  write(out_unitp,*) 'Sij',i,j,dot_product(d0P_cheby(:,j),w_cheby*d0P_cheby(:,i))
  !  CALL flush_perso(out_unitp)
  !END DO
  !END DO
!-----------------------------------------------------------
!STOP
!-----------------------------------------------------------
! then the filter(s) are projected on the chebyshev polynomials
  sigma = para_propa%para_Davidson%W_filter /para_H%Esc
  CALL alloc_NParray(f,[nb,para_propa%para_Davidson%L_filter],'f',name_sub)
  CALL alloc_NParray(wfl,[nq],'wfl',name_sub)

  DO l=1,para_propa%para_Davidson%L_filter
    El = para_propa%para_Davidson%LambdaMin + real(l-1,kind=Rkind)*Delta_Lambda

    El = (El-para_H%E0)/para_H%Esc ! to scale the energy between [-1,1]
    IF (debug) write(out_unitp,*) 'l,El',l,El
    CALL flush_perso(out_unitp)


    ! l-filter values x weight on the grid
    DO k=1,nq
      wfl(k) = f_filter_gauss(x_cheby(k),El,sigma) * w_cheby(k)
    END DO

    DO j=1,nb
      f(j,l) = dot_product(wfl,d0P_cheby(:,j))
    END DO

    IF (debug) write(out_unitp,*) 'l,end err of f',l,sum(abs(f(nb-10:nb,l)))/TEN
    CALL flush_perso(out_unitp)
    !write(out_unitp,*) '========================================='
    !write(out_unitp,*) 'f(:,l)',l
    !write(out_unitp,'(10f10.6)') ,f(:,l)

    !DO k=1,nq
    !  ff_filter = dot_product(f(:,l),d0P_cheby(k,:))
    !  wfl(k) = f_filter_gauss(x_cheby(k),El,sigma)
    !  write(out_unitp6,*) x_cheby(k),wfl(k),ff_filter,log10(abs(ff_filter-wfl(k)))
    !END DO

  END DO
  IF (debug) write(out_unitp,*) 'chebychev coef: done',nb
  CALL flush_perso(out_unitp)
!-----------------------------------------------------------

!- vector initialization:  q1 -----------
  CALL init_psi(q1,para_H,para_H%cplx)
  CALL init_psi(g,para_H,para_H%cplx)

  CALL alloc_psi(q1)
  CALL Set_Random_psi(q1)
  CALL Set_symab_OF_psiBasisRep(q1,para_propa%para_Davidson%symab)
  CALL renorm_psi(q1,BasisRep=.TRUE.)
!- vector initialization:  q1 -----------



!----------------------------------------------------------

     !Conv = .FALSE.
     !DO

       !- chebychev recursion -------------------------------------------
       CALL sub_chebychev_recursion(Tnq1,q1,0,nb-1,para_H)
       IF (debug) write(out_unitp,*) 'chebychev recursion: done',nb
       CALL flush_perso(out_unitp)

       !- Coefficient of the chebychev expansion (with the filter)  -----

       !- z vectors --------------------------------------------
       CALL sub_newZ_vectors_withf(z,Tnq1,nb,para_H,para_propa,f)
       nb_diago = count(abs(z(:)%CAvOp)>ONETENTH**9)
       IF (debug) write(out_unitp,*) 'ortho z vectors: done',nb_diago
       CALL flush_perso(out_unitp)
       !- z vectors -------------------------------------------

        !- diagonalization -----------------------------------
        CALL alloc_NParray(H,  [nb_diago,nb_diago],'H',  name_sub)
        CALL alloc_NParray(Vec,[nb_diago,nb_diago],'Vec',name_sub)

        DO j=1,nb_diago
        DO i=1,nb_diago
          CALL Overlap_psi1_psi2(Overlap,z(i),z(j))
          H(i,j) = real(Overlap,kind=Rkind)
        END DO
        END DO
        !IF (debug) write(out_unitp,*) 'S matrix: done'
        CALL flush_perso(out_unitp)
        CALL sub_ana_S(H,nb_diago,max_Sii,max_Sij,.TRUE.)
        !IF (debug) CALL Write_Mat(H,out_unitp,5)



        DO j=1,nb_diago
          !H.z(j) calc
          CALL sub_OpPsi(z(j),Hz(j),para_H)                            ! => H.z(j)

          DO i=1,nb_diago
            CALL Overlap_psi1_psi2(Overlap,z(i),Hz(j))
            H(i,j) = real(Overlap,kind=Rkind)
          END DO
        END DO
        !IF (debug) write(out_unitp,*) 'H matrix: done'
        CALL flush_perso(out_unitp)

        CALL sub_hermitic_H(H,nb_diago,non_hermitic,para_H%sym_Hamil)
        !IF (debug) CALL Write_Mat(H,out_unitp,5)

        IF (non_hermitic > FOUR*ONETENTH**4) THEN
          If(MPI_id==0) write(out_unitp,*) 'WARNING: non_hermitic is BIG'
          If(MPI_id==0) write(out_unitp,31) non_hermitic
 31       format(' Hamiltonien: ',f16.12,' au')
        ELSE
          If(MPI_id==0) write(out_unitp,51) non_hermitic*auTOcm_inv
 51       format(' Hamiltonien: ',f16.12,' cm-1')
        END IF
        epsi = max(para_propa%para_Davidson%conv_resi,                    &
                   TEN**para_propa%para_Davidson%conv_hermitian *       &
                                               non_hermitic)

        IF (para_H%sym_Hamil) THEN
          CALL diagonalization(H,Ene(1:nb_diago),Vec,nb_diago,3,1,.FALSE.)
        ELSE
          CALL diagonalization(H,Ene(1:nb_diago),Vec,nb_diago,4,1,.FALSE.)
        END IF

         write(out_unitp,21) Ene(1:nb_diago)*auTOcm_inv
         write(out_unitp,21) (Ene(1:nb_diago)-para_H%ZPE)*auTOcm_inv
 21      format(' Filter: ',50(1x,f18.4))


        !----------------------------------------------------------
        !- residual vector ---------------------------
        IF (debug) write(out_unitp,*) 'residual'
        CALL flush_perso(out_unitp)
        DO j=1,nb_diago
          g = ZERO
          DO i=1,nb_diago
            w1 = Hz(i) - z(i) * Ene(j)
            g = g + w1 * Vec(i,j)
          END DO
          CALL norm2_psi(g)
          tab_norm2g(j) = sqrt(g%norm2)
          convergeResi(j) = tab_norm2g(j) < epsi

        END DO
        write(out_unitp,41) 'tab_norm2g(:):      ',tab_norm2g(1:nb_diago)
        write(out_unitp,42) 'convergenceResi(:): ',convergeResi(1:nb_diago)
 41     format(a,100(1x,e9.2))
 42     format(a,100(1x,l9))
        IF (debug) write(out_unitp,*) 'residual: done'

        nb_Vec_IN_Window     = 0
        nb_ConvVec_IN_Window = 0
        DO j=1,nb_diago
          IF (Ene(j) >= para_propa%para_Davidson%LambdaMin .AND.         &
              Ene(j) <= para_propa%para_Davidson%LambdaMax) THEN

            nb_Vec_IN_Window = nb_Vec_IN_Window + 1
            IF (convergeResi(j)) nb_ConvVec_IN_Window = nb_ConvVec_IN_Window + 1
          END IF
        END DO
        Conv = (nb_ConvVec_IN_Window == nb_Vec_IN_Window) .AND. nb_ConvVec_IN_Window > 0

        write(out_unitp,*)  ' Converged levels:               ',count(convergeResi(1:nb_diago))
        write(out_unitp,*)  ' Converged levels in the windox: ',nb_ConvVec_IN_Window
        write(out_unitp,*)  ' Levels in the window:           ',nb_Vec_IN_Window
        write(out_unitp,*)  ' Convergence ?:                  ',Conv
        CALL flush_perso(out_unitp)
        !- residual vector and convergence ------------------------
        !----------------------------------------------------------
        !IF (Conv .OR. mf >= para_propa%para_Davidson%Mmax_filter) EXIT

        CALL dealloc_NParray(H,  'H',  name_sub)
        CALL dealloc_NParray(Vec,'Vec',name_sub)

STOP
        m0 = mf + 1
        mf = min(mf + para_propa%para_Davidson%DeltaM_filter,           &
                                   para_propa%para_Davidson%Mmax_filter)

        IF (nb_diago == para_propa%para_Davidson%L_filter) THEN
          para_propa%para_Davidson%L_filter = para_propa%para_Davidson%L_filter + DeltaL
        END IF
        IF (para_propa%para_Davidson%L_filter > para_propa%para_Davidson%Lmax_filter) THEN
          para_propa%para_Davidson%L_filter = para_propa%para_Davidson%Lmax_filter
        END IF
        IF (para_propa%para_Davidson%L_filter > max_diago) THEN
          para_propa%para_Davidson%L_filter = max_diago
        END IF
        STOP
      !END DO

      !----------------------------------------------------------
      !! save converged vectors and vectors in the window
      jsave = 0
      DO j=1,nb_diago
        IF (convergeResi(j) .OR.                                        &
                    (Ene(j) >= para_propa%para_Davidson%LambdaMin .AND. &
                     Ene(j) <= para_propa%para_Davidson%LambdaMax)) THEN
          jsave = jsave + 1

          CALL init_psi(psi(jsave),para_H,para_H%cplx)
          psi(jsave) = ZERO
          DO i=1,nb_diago
            psi(jsave) = psi(jsave) + Vec(i,j) * z(i)
          END DO
          psi(jsave)%CAvOp    = Ene(j)
          psi(jsave)%IndAvOp  = para_H%n_Op  ! it should be 0
          psi(jsave)%convAvOp = convergeResi(j)
        END IF
      END DO
      nb_diago = jsave
      IF (nb_diago < 1) THEN
        write(out_unitp,*) 'WARNING in ',name_sub
        write(out_unitp,*) ' There is no vector in the filter range!!'
        !STOP 'filter diago'
      END IF

     CALL dealloc_psi(w1)
     CALL dealloc_psi(g)
     CALL dealloc_psi(q1)
     DO i=1,size(z)
       CALL dealloc_psi(z(i))
     END DO
     DO i=1,size(Hz)
       CALL dealloc_psi(Hz(i))
     END DO
     DO i=lbound(Tnq1,dim=1),ubound(Tnq1,dim=1)
       CALL dealloc_psi(Tnq1(i))
     END DO
     CALL dealloc_NParray(H,  'H',  name_sub)
     CALL dealloc_NParray(Vec,'Vec',name_sub)

!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------


      END SUBROUTINE sub_BlockFilterDiagonalization


      SUBROUTINE sub_GaussianFilterDiagonalization_v0(psi,Ene,nb_diago,max_diago, &
                                                   para_H,para_propa)
      USE mod_system
      USE mod_psi,    ONLY : param_psi,alloc_psi,Set_Random_psi,        &
                             Set_symab_OF_psiBasisRep,renorm_psi,       &
                             Overlap_psi1_psi2,norm2_psi,dealloc_psi
      USE mod_Op
      USE mod_propa
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ----------------------------
      integer                      :: nb_diago,max_diago
      TYPE (param_propa)           :: para_propa
      TYPE (param_psi)             :: psi(max_diago)
      real (kind=Rkind)            :: Ene(max_diago)
!------ working parameters --------------------------------


      real (kind=Rkind) :: phi_j,Delta_Lambda,Lambda

      integer       :: i,n,jorth,jsave

      TYPE (param_psi)             :: q1 ! intial vector (random ?)
      TYPE (param_psi)             :: z(para_propa%para_Davidson%Lmax_filter) ! intial vector (random ?)
      TYPE (param_psi)             :: Hz(para_propa%para_Davidson%Lmax_filter) ! intial vector (random ?)

      TYPE (param_psi)             :: Tnq1(0:para_propa%para_Davidson%Mmax_filter) ! intial vector (random ?)
      TYPE (param_psi)             :: g,w1 ! working vector

      real (kind=Rkind), allocatable :: f(:,:)
      real (kind=Rkind) :: acEj,acDeltaj,sigma,ff_filter,El
      integer :: j,k,l,JJ


      real (kind=Rkind), allocatable :: H(:,:),Vec(:,:)

      logical                   :: Conv,convergeEne(max_diago),convergeResi(max_diago)
      real (kind=Rkind)         :: non_hermitic,epsi,auTOcm_inv
      real (kind=Rkind)         :: RS,a,max_Sii,max_Sij
      real (kind=Rkind)         :: tab_norm2g(max_diago)
      complex (kind=Rkind)      :: Overlap
      integer :: m0,mf,nb_Vec_IN_Window,nb_ConvVec_IN_Window,DeltaL

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_GaussianFilterDiagonalization_v0'
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      para_propa%para_Davidson%E0_filter =                              &
                   para_propa%para_Davidson%E0_filter + para_propa%Hmin
      nb_diago = para_propa%para_Davidson%L_filter
      DeltaL   = para_propa%para_Davidson%L_filter
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Hmin,Hmax     : ',para_propa%Hmin,para_propa%Hmax
        write(out_unitp,*) 'E0_filter (ua): ',para_propa%para_Davidson%E0_filter
        write(out_unitp,*) 'L_filter      : ',para_propa%para_Davidson%L_filter
        write(out_unitp,*) 'Lmax_filter   : ',para_propa%para_Davidson%Lmax_filter
        write(out_unitp,*) 'M_filter      : ',para_propa%para_Davidson%M_filter
        write(out_unitp,*) 'Mmax_filter   : ',para_propa%para_Davidson%Mmax_filter
        !write(out_unitp,*) 'nb_diago      : ',nb_diago
        write(out_unitp,*) 'max_diago     : ',max_diago

      END IF
!-----------------------------------------------------------
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')
!-----------------------------------------------------------

      write(out_unitp,*) ' Propagation: ',para_propa%name_WPpropa
      CALL Set_ZPE_OF_Op(para_H,ZPE=para_propa%Hmin,forced=.TRUE.)
      write(out_unitp,*) 'ZPE (cm-1)',para_H%ZPE * auTOcm_inv

      ! change Hmin and Hmax to be sure that the spectral range is between Hmin and Hmax.
      para_propa%Hmin = para_propa%Hmin - ONETENTH**2 * (para_propa%Hmax - para_propa%Hmin)
      para_propa%Hmax = para_propa%Hmax + ONETENTH**2 * (para_propa%Hmax - para_propa%Hmin)

!     - scaling of H ---------------------------------------
      para_propa%para_poly%deltaE = para_propa%Hmax - para_propa%Hmin
      para_propa%para_poly%E0     = para_propa%Hmin + HALF * para_propa%para_poly%deltaE
      para_propa%para_poly%Esc    = HALF * para_propa%para_poly%deltaE

      para_H%scaled = .TRUE.
      para_H%E0     = para_propa%para_poly%E0
      para_H%Esc    = para_propa%para_poly%Esc

      write(out_unitp,*) ' deltaE,E0,Esc: ',para_propa%para_poly%deltaE,&
                                            para_H%E0,para_H%Esc
!-----------------------------------------------------------

      !- vector initialization:  q1 + others -----------
      CALL init_psi(q1,para_H,para_H%cplx)
      CALL init_psi(g,para_H,para_H%cplx)

      CALL alloc_psi(q1)
      CALL Set_Random_psi(q1)
      CALL Set_symab_OF_psiBasisRep(q1,para_propa%para_Davidson%symab)
      CALL renorm_psi(q1,BasisRep=.TRUE.)
      !- vector initialization:  q1 + others -----------

     write(out_unitp,*) 'W_filter (ua)  : ',para_propa%para_Davidson%W_filter
     write(out_unitp,*) 'W_filter (cm-1): ',para_propa%para_Davidson%W_filter*&
                                                                   auTOcm_inv


     Delta_Lambda = para_propa%para_Davidson%W_filter /                 &
                    real(para_propa%para_Davidson%L_filter-1,kind=Rkind)
     write(out_unitp,*) ' Delta_Lambda (ua)  : ',Delta_Lambda
     write(out_unitp,*) ' Delta_Lambda (cm-1): ',Delta_Lambda*auTOcm_inv

     para_propa%para_Davidson%LambdaMin = max(para_H%ZPE,         &
         para_propa%para_Davidson%E0_filter - HALF*para_propa%para_Davidson%W_filter)
     para_propa%para_Davidson%LambdaMax =                               &
         para_propa%para_Davidson%LambdaMin + para_propa%para_Davidson%W_filter
     write(out_unitp,*) ' LambdaMin,LambdaMax (cm-1): ',                &
                        para_propa%para_Davidson%LambdaMin * auTOcm_inv,&
                        para_propa%para_Davidson%LambdaMax * auTOcm_inv
!----------------------------------------------------------

     m0 = 0
     mf = para_propa%para_Davidson%M_filter
     Conv = .FALSE.
     JJ = mf+1  ! J=K+1   in Guo paper

     sigma = para_propa%para_Davidson%W_filter /para_H%Esc / para_propa%para_Davidson%L_filter
     CALL alloc_NParray(f,[mf,para_propa%para_Davidson%L_filter],'f',name_sub,[0,1])

     DO

       !- chebychev recursion -------------------------------------------
       CALL sub_chebychev_recursion(Tnq1,q1,m0,mf,para_H)
       IF (debug) write(out_unitp,*) 'chebychev recursion: done',mf
       CALL flush_perso(out_unitp)

       !- Coefficient of the chebychev expansion (with the filter)  -----
       DO l=1,para_propa%para_Davidson%L_filter
         El = para_propa%para_Davidson%LambdaMin + real(l-1,kind=Rkind)*Delta_Lambda

         El = (El-para_H%E0)/para_H%Esc ! to scale the energy between [-1,1]

         IF (debug) write(out_unitp,*) 'l,El',l,El
         CALL flush_perso(out_unitp)
         DO k=0,mf
           f(k,l) = ZERO

           acDeltaj = pi/real(JJ,kind=Rkind)
           acEj = -HALF*acDeltaj
           !write(out_unitp,*) 'acDeltaj,acEj',acDeltaj,acEj
           DO j=1,JJ
             acEj = acEj + acDeltaj

             f(k,l) = f(k,l) + f_filter_gauss(cos(acEj),El,sigma)*cos(real(k,kind=Rkind)*acEj)

           END DO
           f(k,l) = f(k,l) / real(JJ,kind=Rkind)
           IF (k > 0) THEN
             f(k,l) = f(k,l) * TWO
           END IF
         END DO
         IF (debug) write(out_unitp,*) 'l,end err of f',l,sum(abs(f(mf-10:mf,l)))/TEN
         CALL flush_perso(out_unitp)
         !write(out_unitp,*) '========================================='
         !write(out_unitp,*) 'f(:,l)',l
         !write(out_unitp,'(10f10.6)') ,f(:,l)

         acDeltaj = pi/real(JJ,kind=Rkind)
         acEj = -HALF*acDeltaj
         DO j=1,JJ
           acEj = acEj + acDeltaj

           ff_filter = f(0,l)

           DO k=1,para_propa%para_Davidson%M_filter
             ff_filter = ff_filter + f(k,l) * cos(real(k,kind=Rkind)*acEj)
           END DO
           write(66,*) cos(acEj),f_filter_gauss(cos(acEj),El,sigma),    &
             ff_filter,log10(abs(ff_filter-f_filter_gauss(cos(acEj),El,sigma)))
           !write(67,*) cos(acEj),log10(f_filter_gauss(cos(acEj),El,sigma)),&
           !                                            log10(ff_filter)
         END DO
         STOP
       END DO
       IF (debug) write(out_unitp,*) 'chebychev coef: done',mf
       CALL flush_perso(out_unitp)

       !- z vectors --------------------------------------------
       CALL sub_Z_vectors_withf(z,Tnq1,mf,para_H,para_propa,f)
       nb_diago = count(abs(z(:)%CAvOp)>ONETENTH**9)
       IF (debug) write(out_unitp,*) 'ortho z vectors: done',nb_diago
       CALL flush_perso(out_unitp)
       !- z vectors -------------------------------------------

        !- diagonalization -----------------------------------
        CALL alloc_NParray(H,  [nb_diago,nb_diago],'H',  name_sub)
        CALL alloc_NParray(Vec,[nb_diago,nb_diago],'Vec',name_sub)

        DO j=1,nb_diago
        DO i=1,nb_diago
          CALL Overlap_psi1_psi2(Overlap,z(i),z(j))
          H(i,j) = real(Overlap,kind=Rkind)
        END DO
        END DO
        !IF (debug) write(out_unitp,*) 'S matrix: done'
        CALL flush_perso(out_unitp)
        CALL sub_ana_S(H,nb_diago,max_Sii,max_Sij,.TRUE.)
        !IF (debug) CALL Write_Mat(H,out_unitp,5)



        DO j=1,nb_diago
          !H.z(j) calc
          CALL sub_OpPsi(z(j),Hz(j),para_H)                            ! => H.z(j)

          DO i=1,nb_diago
            CALL Overlap_psi1_psi2(Overlap,z(i),Hz(j))
            H(i,j) = real(Overlap,kind=Rkind)
          END DO
        END DO
        !IF (debug) write(out_unitp,*) 'H matrix: done'
        CALL flush_perso(out_unitp)

        CALL sub_hermitic_H(H,nb_diago,non_hermitic,para_H%sym_Hamil)
        !IF (debug) CALL Write_Mat(H,out_unitp,5)

        IF (non_hermitic > FOUR*ONETENTH**4) THEN
          If(MPI_id==0) write(out_unitp,*) 'WARNING: non_hermitic is BIG'
          If(MPI_id==0) write(out_unitp,31) non_hermitic
 31       format(' Hamiltonien: ',f16.12,' au')
        ELSE
          If(MPI_id==0) write(out_unitp,51) non_hermitic*auTOcm_inv
 51       format(' Hamiltonien: ',f16.12,' cm-1')
        END IF
        epsi = max(para_propa%para_Davidson%conv_resi,                    &
                   TEN**para_propa%para_Davidson%conv_hermitian *       &
                                               non_hermitic)

        IF (para_H%sym_Hamil) THEN
          CALL diagonalization(H,Ene(1:nb_diago),Vec,nb_diago,3,1,.FALSE.)
        ELSE
          CALL diagonalization(H,Ene(1:nb_diago),Vec,nb_diago,4,1,.FALSE.)
        END IF

         write(out_unitp,21) Ene(1:nb_diago)*auTOcm_inv
         write(out_unitp,21) (Ene(1:nb_diago)-para_H%ZPE)*auTOcm_inv
 21      format(' Filter: ',50(1x,f18.4))


        !----------------------------------------------------------
        !- residual vector ---------------------------
        IF (debug) write(out_unitp,*) 'residual'
        CALL flush_perso(out_unitp)
        DO j=1,nb_diago
          g = ZERO
          DO i=1,nb_diago
            w1 = Hz(i) - z(i) * Ene(j)
            g = g + w1 * Vec(i,j)
          END DO
          CALL norm2_psi(g)
          tab_norm2g(j) = sqrt(g%norm2)
          convergeResi(j) = tab_norm2g(j) < epsi

        END DO
        write(out_unitp,41) 'tab_norm2g(:):      ',tab_norm2g(1:nb_diago)
        write(out_unitp,42) 'convergenceResi(:): ',convergeResi(1:nb_diago)
 41     format(a,100(1x,e9.2))
 42     format(a,100(1x,l9))
        IF (debug) write(out_unitp,*) 'residual: done'

        nb_Vec_IN_Window     = 0
        nb_ConvVec_IN_Window = 0
        DO j=1,nb_diago
          IF (Ene(j) >= para_propa%para_Davidson%LambdaMin .AND.         &
              Ene(j) <= para_propa%para_Davidson%LambdaMax) THEN

            nb_Vec_IN_Window = nb_Vec_IN_Window + 1
            IF (convergeResi(j)) nb_ConvVec_IN_Window = nb_ConvVec_IN_Window + 1
          END IF
        END DO
        Conv = (nb_ConvVec_IN_Window == nb_Vec_IN_Window) .AND. nb_ConvVec_IN_Window > 0

        write(out_unitp,*)  ' Converged levels:               ',count(convergeResi(1:nb_diago))
        write(out_unitp,*)  ' Converged levels in the windox: ',nb_ConvVec_IN_Window
        write(out_unitp,*)  ' Levels in the window:           ',nb_Vec_IN_Window
        write(out_unitp,*)  ' Convergence ?:                  ',Conv
        CALL flush_perso(out_unitp)
        !- residual vector and convergence ------------------------
        !----------------------------------------------------------
        IF (Conv .OR. mf >= para_propa%para_Davidson%Mmax_filter) EXIT

        CALL dealloc_NParray(H,  'H',  name_sub)
        CALL dealloc_NParray(Vec,'Vec',name_sub)
        m0 = mf + 1
        mf = min(mf + para_propa%para_Davidson%DeltaM_filter,           &
                                   para_propa%para_Davidson%Mmax_filter)

        IF (nb_diago == para_propa%para_Davidson%L_filter) THEN
          para_propa%para_Davidson%L_filter = para_propa%para_Davidson%L_filter + DeltaL
        END IF
        IF (para_propa%para_Davidson%L_filter > para_propa%para_Davidson%Lmax_filter) THEN
          para_propa%para_Davidson%L_filter = para_propa%para_Davidson%Lmax_filter
        END IF
        IF (para_propa%para_Davidson%L_filter > max_diago) THEN
          para_propa%para_Davidson%L_filter = max_diago
        END IF
        STOP
      END DO

      !----------------------------------------------------------
      !! save converged vectors and vectors in the window
      jsave = 0
      DO j=1,nb_diago
        IF (convergeResi(j) .OR.                                        &
                    (Ene(j) >= para_propa%para_Davidson%LambdaMin .AND. &
                     Ene(j) <= para_propa%para_Davidson%LambdaMax)) THEN
          jsave = jsave + 1

          CALL init_psi(psi(jsave),para_H,para_H%cplx)
          psi(jsave) = ZERO
          DO i=1,nb_diago
            psi(jsave) = psi(jsave) + Vec(i,j) * z(i)
          END DO
          psi(jsave)%CAvOp    = Ene(j)
          psi(jsave)%IndAvOp  = para_H%n_Op  ! it should be 0
          psi(jsave)%convAvOp = convergeResi(j)
        END IF
      END DO
      nb_diago = jsave
      IF (nb_diago < 1) THEN
        write(out_unitp,*) 'WARNING in ',name_sub
        write(out_unitp,*) ' There is no vector in the filter range!!'
        !STOP 'filter diago'
      END IF

     CALL dealloc_psi(w1)
     CALL dealloc_psi(g)
     CALL dealloc_psi(q1)
     DO i=1,size(z)
       CALL dealloc_psi(z(i))
     END DO
     DO i=1,size(Hz)
       CALL dealloc_psi(Hz(i))
     END DO
     DO i=lbound(Tnq1,dim=1),ubound(Tnq1,dim=1)
       CALL dealloc_psi(Tnq1(i))
     END DO
     CALL dealloc_NParray(H,  'H',  name_sub)
     CALL dealloc_NParray(Vec,'Vec',name_sub)

!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------


      END SUBROUTINE sub_GaussianFilterDiagonalization_v0

!================================================================
!
!          Filter diagonalization (Cheby)
!          Hmin and Hmax are the parameter to scale H
!
!================================================================
      SUBROUTINE sub_FilterDiagonalization(psi,Ene,nb_diago,max_diago, &
                                           para_H,para_propa)
      USE mod_system
      USE mod_psi,    ONLY : param_psi,alloc_psi,Set_Random_psi,        &
                             Set_symab_OF_psiBasisRep,renorm_psi,       &
                             Overlap_psi1_psi2,norm2_psi,dealloc_psi
      USE mod_Op
      USE mod_propa
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ----------------------------
      integer                      :: nb_diago,max_diago
      TYPE (param_propa)           :: para_propa
      TYPE (param_psi)             :: psi(max_diago)
      real (kind=Rkind)            :: Ene(max_diago)
!------ working parameters --------------------------------


      real (kind=Rkind) :: phi_j,Delta_Lambda,Lambda
      !real (kind=Rkind) :: phi_j,Delta_Lambda,Lambda,LambdaMin,LambdaMax

      integer       :: i,j,n,jorth,jsave

      TYPE (param_psi)             :: q1 ! intial vector (random ?)
      TYPE (param_psi)             :: z(para_propa%para_Davidson%Lmax_filter) ! intial vector (random ?)
      TYPE (param_psi)             :: Hz(para_propa%para_Davidson%Lmax_filter) ! intial vector (random ?)

      TYPE (param_psi)             :: Tnq1(0:para_propa%para_Davidson%Mmax_filter) ! intial vector (random ?)
      TYPE (param_psi)             :: g,w1 ! working vector

      real (kind=Rkind), allocatable :: H(:,:),Vec(:,:)

      logical                   :: Conv,convergeEne(max_diago),convergeResi(max_diago)
      real (kind=Rkind)         :: non_hermitic,epsi,auTOcm_inv
      real (kind=Rkind)         :: RS,a,max_Sii,max_Sij
      real (kind=Rkind)         :: tab_norm2g(max_diago)
      complex (kind=Rkind)      :: Overlap
      integer :: m0,mf,nb_Vec_IN_Window,nb_ConvVec_IN_Window,DeltaL

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_FilterDiagonalization'
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      para_propa%para_Davidson%E0_filter =                              &
                   para_propa%para_Davidson%E0_filter + para_propa%Hmin
      nb_diago = para_propa%para_Davidson%L_filter
      DeltaL   = para_propa%para_Davidson%L_filter
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Hmin,Hmax     : ',para_propa%Hmin,para_propa%Hmax
        write(out_unitp,*) 'E0_filter (ua): ',para_propa%para_Davidson%E0_filter
        write(out_unitp,*) 'L_filter      : ',para_propa%para_Davidson%L_filter
        write(out_unitp,*) 'Lmax_filter   : ',para_propa%para_Davidson%Lmax_filter
        write(out_unitp,*) 'M_filter      : ',para_propa%para_Davidson%M_filter
        write(out_unitp,*) 'Mmax_filter   : ',para_propa%para_Davidson%Mmax_filter
        !write(out_unitp,*) 'nb_diago      : ',nb_diago
        write(out_unitp,*) 'max_diago     : ',max_diago

      END IF
!-----------------------------------------------------------
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')
!-----------------------------------------------------------

      write(out_unitp,*) ' Propagation: ',para_propa%name_WPpropa
      CALL Set_ZPE_OF_Op(para_H,ZPE=para_propa%Hmin,forced=.TRUE.)
      write(out_unitp,*) 'ZPE (cm-1)',para_H%ZPE * auTOcm_inv

      ! change Hmin and Hmax to be sure that the spectral range is between Hmin and Hmax.
      para_propa%Hmin = para_propa%Hmin - ONETENTH**2 * (para_propa%Hmax - para_propa%Hmin)
      para_propa%Hmax = para_propa%Hmax + ONETENTH**2 * (para_propa%Hmax - para_propa%Hmin)

!     - scaling of H ---------------------------------------
      para_propa%para_poly%deltaE = para_propa%Hmax - para_propa%Hmin
      para_propa%para_poly%E0     = para_propa%Hmin + HALF * para_propa%para_poly%deltaE
      para_propa%para_poly%Esc    = HALF * para_propa%para_poly%deltaE

      para_H%scaled = .TRUE.
      para_H%E0     = para_propa%para_poly%E0
      para_H%Esc    = para_propa%para_poly%Esc

      write(out_unitp,*) ' deltaE,E0,Esc: ',para_propa%para_poly%deltaE,&
                                            para_H%E0,para_H%Esc
!-----------------------------------------------------------

      !- vector initialization:  q1 + others -----------
      CALL init_psi(q1,para_H,para_H%cplx)
      CALL init_psi(g,para_H,para_H%cplx)

      CALL alloc_psi(q1)
      IF (q1%cplx) THEN
        DO i=1,q1%nb_tot
          CALL random_number(a)
          q1%CvecB(i) = cmplx(a,ZERO,kind=Rkind)
        END DO
      ELSE
        DO i=1,q1%nb_tot
          CALL random_number(a)
          q1%RvecB(i) = a
        END DO
      END IF
      CALL Set_symab_OF_psiBasisRep(q1,para_propa%para_Davidson%symab)
      CALL renorm_psi(q1,BasisRep=.TRUE.)
      !- vector initialization:  q1 + others -----------

     write(out_unitp,*) 'W_filter (ua)  : ',para_propa%para_Davidson%W_filter
     write(out_unitp,*) 'W_filter (cm-1): ',para_propa%para_Davidson%W_filter*&
                                                                   auTOcm_inv


     Delta_Lambda = para_propa%para_Davidson%W_filter /                 &
                    real(para_propa%para_Davidson%L_filter-1,kind=Rkind)
     write(out_unitp,*) ' Delta_Lambda (ua)  : ',Delta_Lambda
     write(out_unitp,*) ' Delta_Lambda (cm-1): ',Delta_Lambda*auTOcm_inv

     para_propa%para_Davidson%LambdaMin = max(para_H%ZPE,         &
         para_propa%para_Davidson%E0_filter - HALF*para_propa%para_Davidson%W_filter)
     para_propa%para_Davidson%LambdaMax =                               &
         para_propa%para_Davidson%LambdaMin + para_propa%para_Davidson%W_filter
     write(out_unitp,*) ' LambdaMin,LambdaMax (cm-1): ',                &
                        para_propa%para_Davidson%LambdaMin * auTOcm_inv,&
                        para_propa%para_Davidson%LambdaMax * auTOcm_inv
!----------------------------------------------------------

     m0 = 0
     mf = para_propa%para_Davidson%M_filter
     Conv = .FALSE.

     DO

       !- chebychev recursion -------------------------------------------
       CALL sub_chebychev_recursion(Tnq1,q1,m0,mf,para_H)
       IF (debug) write(out_unitp,*) 'chebychev recursion: done',mf
       CALL flush_perso(out_unitp)

       !- z vectors --------------------------------------------
       CALL sub_Z_vectors(z,Tnq1,m0,mf,para_H,para_propa)
       nb_diago = count(abs(z(:)%CAvOp)>ONETENTH**9)
       IF (debug) write(out_unitp,*) 'ortho z vectors: done',nb_diago
       CALL flush_perso(out_unitp)
       !- z vectors -------------------------------------------

        !- diagonalization -----------------------------------
        CALL alloc_NParray(H,  [nb_diago,nb_diago],'H',  name_sub)
        CALL alloc_NParray(Vec,[nb_diago,nb_diago],'Vec',name_sub)

        DO j=1,nb_diago
        DO i=1,nb_diago
          CALL Overlap_psi1_psi2(Overlap,z(i),z(j))
          H(i,j) = real(Overlap,kind=Rkind)
        END DO
        END DO
        !IF (debug) write(out_unitp,*) 'S matrix: done'
        CALL flush_perso(out_unitp)
        CALL sub_ana_S(H,nb_diago,max_Sii,max_Sij,.TRUE.)
        !IF (debug) CALL Write_Mat(H,out_unitp,5)



        DO j=1,nb_diago
          !H.z(j) calc
          CALL sub_OpPsi(z(j),Hz(j),para_H)                            ! => H.z(j)

          DO i=1,nb_diago
            CALL Overlap_psi1_psi2(Overlap,z(i),Hz(j))
            H(i,j) = real(Overlap,kind=Rkind)
          END DO
        END DO
        !IF (debug) write(out_unitp,*) 'H matrix: done'
        CALL flush_perso(out_unitp)

        CALL sub_hermitic_H(H,nb_diago,non_hermitic,para_H%sym_Hamil)
        !IF (debug) CALL Write_Mat(H,out_unitp,5)

        IF (non_hermitic > FOUR*ONETENTH**4) THEN
          If(MPI_id==0) write(out_unitp,*) 'WARNING: non_hermitic is BIG'
          If(MPI_id==0) write(out_unitp,31) non_hermitic
 31       format(' Hamiltonien: ',f16.12,' au')
        ELSE
          If(MPI_id==0) write(out_unitp,51) non_hermitic*auTOcm_inv
 51       format(' Hamiltonien: ',f16.12,' cm-1')
        END IF
        epsi = max(para_propa%para_Davidson%conv_resi,                    &
                   TEN**para_propa%para_Davidson%conv_hermitian *       &
                                               non_hermitic)

        IF (para_H%sym_Hamil) THEN
          CALL diagonalization(H,Ene(1:nb_diago),Vec,nb_diago,3,1,.FALSE.)
        ELSE
          CALL diagonalization(H,Ene(1:nb_diago),Vec,nb_diago,4,1,.FALSE.)
        END IF

         write(out_unitp,21) Ene(1:nb_diago)*auTOcm_inv
         write(out_unitp,21) (Ene(1:nb_diago)-para_H%ZPE)*auTOcm_inv
 21      format(' Filter: ',50(1x,f18.4))


        !----------------------------------------------------------
        !- residual vector ---------------------------
        IF (debug) write(out_unitp,*) 'residual'
        CALL flush_perso(out_unitp)
        DO j=1,nb_diago
          g = ZERO
          DO i=1,nb_diago
            w1 = Hz(i) - z(i) * Ene(j)
            g = g + w1 * Vec(i,j)
          END DO
          CALL norm2_psi(g)
          tab_norm2g(j) = sqrt(g%norm2)
          convergeResi(j) = tab_norm2g(j) < epsi

        END DO
        write(out_unitp,41) 'tab_norm2g(:):      ',tab_norm2g(1:nb_diago)
        write(out_unitp,42) 'convergenceResi(:): ',convergeResi(1:nb_diago)
 41     format(a,100(1x,e9.2))
 42     format(a,100(1x,l9))
        IF (debug) write(out_unitp,*) 'residual: done'

        nb_Vec_IN_Window     = 0
        nb_ConvVec_IN_Window = 0
        DO j=1,nb_diago
          IF (Ene(j) >= para_propa%para_Davidson%LambdaMin .AND.         &
              Ene(j) <= para_propa%para_Davidson%LambdaMax) THEN

            nb_Vec_IN_Window = nb_Vec_IN_Window + 1
            IF (convergeResi(j)) nb_ConvVec_IN_Window = nb_ConvVec_IN_Window + 1
          END IF
        END DO
        Conv = (nb_ConvVec_IN_Window == nb_Vec_IN_Window) .AND. nb_ConvVec_IN_Window > 0

        write(out_unitp,*)  ' Converged levels:               ',count(convergeResi(1:nb_diago))
        write(out_unitp,*)  ' Converged levels in the windox: ',nb_ConvVec_IN_Window
        write(out_unitp,*)  ' Levels in the window:           ',nb_Vec_IN_Window
        write(out_unitp,*)  ' Convergence ?:                  ',Conv
        CALL flush_perso(out_unitp)
        !- residual vector and convergence ------------------------
        !----------------------------------------------------------
        IF (Conv .OR. mf >= para_propa%para_Davidson%Mmax_filter) EXIT

        CALL dealloc_NParray(H,  'H',  name_sub)
        CALL dealloc_NParray(Vec,'Vec',name_sub)
        m0 = mf + 1
        mf = min(mf + para_propa%para_Davidson%DeltaM_filter,           &
                                   para_propa%para_Davidson%Mmax_filter)

        IF (nb_diago == para_propa%para_Davidson%L_filter) THEN
          para_propa%para_Davidson%L_filter = para_propa%para_Davidson%L_filter + DeltaL
        END IF
        IF (para_propa%para_Davidson%L_filter > para_propa%para_Davidson%Lmax_filter) THEN
          para_propa%para_Davidson%L_filter = para_propa%para_Davidson%Lmax_filter
        END IF
        IF (para_propa%para_Davidson%L_filter > max_diago) THEN
          para_propa%para_Davidson%L_filter = max_diago
        END IF
      END DO

      !----------------------------------------------------------
      !! save converged vectors and vectors in the window
      jsave = 0
      DO j=1,nb_diago
        IF (convergeResi(j) .OR.                                        &
                    (Ene(j) >= para_propa%para_Davidson%LambdaMin .AND. &
                     Ene(j) <= para_propa%para_Davidson%LambdaMax)) THEN
          jsave = jsave + 1

          CALL init_psi(psi(jsave),para_H,para_H%cplx)
          psi(jsave) = ZERO
          DO i=1,nb_diago
            psi(jsave) = psi(jsave) + Vec(i,j) * z(i)
          END DO
          psi(jsave)%CAvOp    = Ene(j)
          psi(jsave)%IndAvOp  = para_H%n_Op  ! it should be 0
          psi(jsave)%convAvOp = convergeResi(j)
        END IF
      END DO
      nb_diago = jsave
      IF (nb_diago < 1) THEN
        write(out_unitp,*) 'WARNING in ',name_sub
        write(out_unitp,*) ' There is no vector in the filter range!!'
        !STOP 'filter diago'
      END IF

     CALL dealloc_psi(w1)
     CALL dealloc_psi(g)
     CALL dealloc_psi(q1)
     DO i=1,size(z)
       CALL dealloc_psi(z(i))
     END DO
     DO i=1,size(Hz)
       CALL dealloc_psi(Hz(i))
     END DO
     DO i=lbound(Tnq1,dim=1),ubound(Tnq1,dim=1)
       CALL dealloc_psi(Tnq1(i))
     END DO
     CALL dealloc_NParray(H,  'H',  name_sub)
     CALL dealloc_NParray(Vec,'Vec',name_sub)

!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------


      END SUBROUTINE sub_FilterDiagonalization

      SUBROUTINE sub_FilterDiagonalization_v1(psi,Ene,nb_diago,max_diago, &
                                           para_H,para_propa)
      USE mod_system
      USE mod_psi,    ONLY : param_psi,alloc_psi,Set_Random_psi,        &
                             Set_symab_OF_psiBasisRep,renorm_psi,       &
                             Overlap_psi1_psi2,norm2_psi,dealloc_psi,  &
                             sub_Lowdin
      USE mod_Op
      USE mod_propa
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ----------------------------
      integer                      :: nb_diago,max_diago
      TYPE (param_propa)           :: para_propa
      TYPE (param_psi)             :: psi(max_diago)
      real (kind=Rkind)            :: Ene(max_diago)
!------ working parameters --------------------------------


      real (kind=Rkind) :: phi_j(para_propa%para_Davidson%L_filter)
      real (kind=Rkind) :: Delta_Lambda,Lambda,LambdaMin,LambdaMax
      integer       :: i,j,n,jorth,jsave

      TYPE (param_psi)             :: q1 ! intial vector (random ?)
      TYPE (param_psi)             :: z(para_propa%para_Davidson%L_filter) ! intial vector (random ?)
      TYPE (param_psi)             :: Hz(para_propa%para_Davidson%L_filter) ! intial vector (random ?)

      TYPE (param_psi)             :: Tnq1(0:para_propa%para_Davidson%M_filter) ! intial vector (random ?)
      TYPE (param_psi)             :: g,w1 ! working vector

      real (kind=Rkind), allocatable :: H(:,:),Vec(:,:)

      logical                   :: convergeEne(max_diago),convergeResi(max_diago)
      real (kind=Rkind)         :: non_hermitic,epsi,auTOcm_inv
      real (kind=Rkind)         :: RS,a,max_Sii,max_Sij
      real (kind=Rkind)         :: tab_norm2g(max_diago)
      complex (kind=Rkind)      :: Overlap

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_FilterDiagonalization_v1'
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      para_propa%para_Davidson%E0_filter =                               &
                    para_propa%para_Davidson%E0_filter + para_propa%Hmin
      nb_diago = para_propa%para_Davidson%L_filter
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Hmin,Hmax     : ',para_propa%Hmin,para_propa%Hmax
        write(out_unitp,*) 'E0_filter (ua): ',para_propa%para_Davidson%E0_filter
        write(out_unitp,*) 'L_filter      : ',para_propa%para_Davidson%L_filter
        write(out_unitp,*) 'M_filter      : ',para_propa%para_Davidson%M_filter
        write(out_unitp,*) 'nb_diago      : ',nb_diago
        write(out_unitp,*) 'max_diago     : ',max_diago

      END IF
!-----------------------------------------------------------
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')
!-----------------------------------------------------------

      write(out_unitp,*) ' Propagation: ',para_propa%name_WPpropa

!     - scaling of H ---------------------------------------
      para_propa%para_poly%deltaE = para_propa%Hmax - para_propa%Hmin
      para_propa%para_poly%E0     = para_propa%Hmin + HALF * para_propa%para_poly%deltaE
      para_propa%para_poly%Esc = HALF * para_propa%para_poly%deltaE

      para_H%scaled = .TRUE.
      para_H%E0     = para_propa%para_poly%E0
      para_H%Esc    = para_propa%para_poly%Esc

      write(out_unitp,*) ' deltaE,E0,Esc: ',para_propa%para_poly%deltaE,&
                                            para_H%E0,para_H%Esc
!-----------------------------------------------------------

      !- vector initialization:  q1, z(:), Tnq1(:), psi(:) -----------
      CALL init_psi(q1,para_H,para_H%cplx)
      CALL init_psi(g,para_H,para_H%cplx)

      CALL alloc_psi(q1)
      IF (q1%cplx) THEN
        DO i=1,q1%nb_tot
          CALL random_number(a)
          q1%CvecB(i) = cmplx(a,ZERO,kind=Rkind)
        END DO
      ELSE
        DO i=1,q1%nb_tot
          CALL random_number(a)
          q1%RvecB(i) = a
        END DO
      END IF
      CALL Set_symab_OF_psiBasisRep(q1,para_propa%para_Davidson%symab)
      CALL renorm_psi(q1,BasisRep=.TRUE.)

      DO j=1,nb_diago
        CALL init_psi(psi(j),para_H,para_H%cplx)
        CALL init_psi(z(j),para_H,para_H%cplx)
        CALL init_psi(Hz(j),para_H,para_H%cplx)
      END DO
      DO n=0,para_propa%para_Davidson%M_filter
        CALL init_psi(Tnq1(n),para_H,para_H%cplx)
      END DO
      !- vector initialization:  q1, z(:), Tnq1(:), psi(:) -----------



     para_propa%para_Davidson%W_filter =                                &
                    real(para_propa%para_Davidson%L_filter,kind=Rkind)* &
                                        para_propa%para_poly%deltaE /   &
          (1.2_Rkind*real(para_propa%para_Davidson%M_filter,kind=Rkind))
     para_propa%para_Davidson%W_filter = 200._Rkind / auTOcm_inv ! 200 cm-1
     write(out_unitp,*) 'W_filter (ua)',para_propa%para_Davidson%W_filter
     write(out_unitp,*) 'W_filter (cm-1)',para_propa%para_Davidson%W_filter*&
                                                            auTOcm_inv


     Delta_Lambda = para_propa%para_Davidson%W_filter /                 &
                    real(para_propa%para_Davidson%L_filter-1,kind=Rkind)
     write(out_unitp,*) ' Delta_Lambda (ua)  : ',Delta_Lambda
     write(out_unitp,*) ' Delta_Lambda (cm-1): ',Delta_Lambda* auTOcm_inv

     LambdaMin = para_propa%para_Davidson%E0_filter -                   &
                                 HALF*para_propa%para_Davidson%W_filter
     IF (LambdaMin < para_propa%Hmin) LambdaMin = para_propa%Hmin
     LambdaMax = LambdaMin + para_propa%para_Davidson%W_filter
     write(out_unitp,*) ' LambdaMin,LambdaMax (cm-1): ',                &
                           LambdaMin * auTOcm_inv,LambdaMax * auTOcm_inv

     DO j=1,para_propa%para_Davidson%L_filter
       Lambda = LambdaMin + real(j-1,kind=Rkind)*Delta_Lambda
       !write(out_unitp,*) 'j,Lambda (cm-1)',j,Lambda * auTOcm_inv
       phi_j(j) = acos((Lambda-para_H%E0)/para_H%Esc)
     END DO
     IF (debug) write(out_unitp,*) 'phi_j(:): done'
     CALL flush_perso(out_unitp)

!----------------------------------------------------------

      !- chebychev recursion -------------------------------------------
      Tnq1(0)  = q1                                                ! => q1
      CALL sub_OpPsi(q1,Tnq1(1),para_H)                            ! => H.q1
      CALL sub_scaledOpPsi(q1,Tnq1(1),para_H%E0,para_H%Esc)        ! scaling

      DO n=2,para_propa%para_Davidson%M_filter
        CALL sub_OpPsi(Tnq1(n-1),w1,para_H)                        ! => Ti-1.q1
        CALL sub_scaledOpPsi(Tnq1(n-1),w1,para_H%E0,para_H%Esc)    ! scaling
        Tnq1(n) = TWO * w1
        Tnq1(n) = Tnq1(n) - Tnq1(n-2)
      END DO
     IF (debug) write(out_unitp,*) 'chebychev recursion: done'
     CALL flush_perso(out_unitp)

      !- chebychev recursion -------------------------------------------

      !- z vectors --------------------------------------------
      DO j=1,para_propa%para_Davidson%L_filter
        z(j) = Tnq1(0)
        DO n=1,para_propa%para_Davidson%M_filter
          z(j) = z(j) + TWO*cos(real(n,kind=rkind)*phi_j(j)) * Tnq1(n)
        END DO
      END DO
     IF (debug) write(out_unitp,*) 'z vectors: done'
     CALL flush_perso(out_unitp)
      !- z vectors -------------------------------------------

      !- z vectors (orthonormalized) -----------------------------------
      CALL sub_Lowdin(z,para_propa%para_Davidson%L_filter)
      nb_diago = count(abs(z(:)%CAvOp)>ONETENTH**8)


     IF (debug) write(out_unitp,*) 'ortho z vectors: done',nb_diago
     CALL flush_perso(out_unitp)
      !- z vectors -------------------------------------------

      !- diagonalization -----------------------------------
      CALL alloc_NParray(H,  [nb_diago,nb_diago],'H',  name_sub)
      CALL alloc_NParray(Vec,[nb_diago,nb_diago],'Vec',name_sub)

      DO j=1,nb_diago
      DO i=1,nb_diago
        CALL Overlap_psi1_psi2(Overlap,z(i),z(j))
        H(i,j) = real(Overlap,kind=Rkind)
      END DO
      END DO
      !IF (debug) write(out_unitp,*) 'H matrix: done'
      CALL flush_perso(out_unitp)
      CALL sub_ana_S(H,nb_diago,max_Sii,max_Sij,.TRUE.)

      !IF (debug) CALL Write_Mat(H,out_unitp,5)



      DO j=1,nb_diago
        !H.z(j) calc
        CALL sub_OpPsi(z(j),Hz(j),para_H)                            ! => H.z(j)
        !CALL sub_scaledOpPsi(z(j),Hz(j),para_H%E0,para_H%Esc)    ! scaling

        DO i=1,nb_diago
          CALL Overlap_psi1_psi2(Overlap,z(i),Hz(j))
          H(i,j) = real(Overlap,kind=Rkind)
        END DO
      END DO
      !IF (debug) write(out_unitp,*) 'H matrix: done'
      CALL flush_perso(out_unitp)

      CALL sub_hermitic_H(H,nb_diago,non_hermitic,para_H%sym_Hamil)
      !IF (debug) CALL Write_Mat(H,out_unitp,5)

      IF (non_hermitic > FOUR*ONETENTH**4) THEN
        If(MPI_id==0) write(out_unitp,*) 'WARNING: non_hermitic is BIG'
        If(MPI_id==0) write(out_unitp,31) non_hermitic
 31     format(' Hamiltonien: ',f16.12,' au')
      ELSE
        If(MPI_id==0) write(out_unitp,51) non_hermitic*auTOcm_inv
 51     format(' Hamiltonien: ',f16.12,' cm-1')
      END IF
      epsi = max(para_propa%para_Davidson%conv_resi,                    &
                   TEN**para_propa%para_Davidson%conv_hermitian *       &
                                               non_hermitic)

      CALL Set_ZPE_OF_Op(para_H,ZPE=para_propa%Hmin,forced=.TRUE.)
      write(out_unitp,*) 'ZPE',para_H%ZPE

        IF (para_H%sym_Hamil) THEN
          CALL diagonalization(H,Ene(1:nb_diago),Vec,nb_diago,3,1,.FALSE.)
        ELSE
          CALL diagonalization(H,Ene(1:nb_diago),Vec,nb_diago,4,1,.FALSE.)
        END IF
      !write(out_unitp,*) Ene(1:nb_diago)*auTOcm_inv

       write(out_unitp,21) Ene(1:nb_diago)*auTOcm_inv
       write(out_unitp,21) (Ene(1:nb_diago)-para_H%ZPE)*auTOcm_inv

 21   format(' Filter: ',50(1x,f18.4))


      !----------------------------------------------------------
      !- residual vector ---------------------------
      IF (debug) write(out_unitp,*) 'residual'
      IF (debug) CALL flush_perso(out_unitp)
      DO j=1,nb_diago
          g = ZERO
          DO i=1,nb_diago
            w1 = Hz(i) - z(i) * Ene(j)
            g = g + w1 * Vec(i,j)
          END DO
          CALL norm2_psi(g)
          tab_norm2g(j) = sqrt(g%norm2)
          convergeResi(j) = tab_norm2g(j) < epsi

      END DO
      write(out_unitp,41) 'tab_norm2g          ',tab_norm2g(1:nb_diago)
      write(out_unitp,42)  'convergenceResi(:): ',convergeResi(1:nb_diago)
 41   format(a,100(1x,e9.2))
 42   format(a,100(1x,l9))
      IF (debug) write(out_unitp,*) 'residual: done'
      write(out_unitp,*)  ' number of converged levels: ',count(convergeResi(1:nb_diago))

      IF (debug) CALL flush_perso(out_unitp)
      !- residual vector and convergence ------------------------
      !----------------------------------------------------------

      jsave = 0
      DO j=1,nb_diago
        IF (convergeResi(j) .OR. (Ene(j) >= LambdaMin .AND. Ene(j) <= LambdaMax)) THEN
          jsave = jsave + 1
          psi(jsave) = ZERO
          DO i=1,nb_diago
            psi(jsave) = psi(jsave) + Vec(i,j) * z(i)
          END DO
          psi(jsave)%CAvOp    = Ene(j)
          psi(jsave)%IndAvOp  = para_H%n_Op  ! it should be 0
          psi(jsave)%convAvOp = convergeResi(j)
        END IF
      END DO
      nb_diago = jsave
      IF (nb_diago < 1) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) ' There is no vector in the filter range!!'
        STOP 'filter diago'
      END IF


!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------


      END SUBROUTINE sub_FilterDiagonalization_v1
      SUBROUTINE sub_FilterDiagonalization_v0(psi,Ene,nb_diago,max_diago, &
                                           para_H,para_propa)
      USE mod_system
      USE mod_psi,    ONLY : param_psi,alloc_psi,Set_Random_psi,        &
                             Set_symab_OF_psiBasisRep,renorm_psi,       &
                             Overlap_psi1_psi2,norm2_psi,dealloc_psi
      USE mod_Op
      USE mod_propa
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_H

!----- variables for the WP propagation ----------------------------
      integer                      :: nb_diago,max_diago
      TYPE (param_propa)           :: para_propa
      TYPE (param_psi)             :: psi(max_diago)
      real (kind=Rkind)            :: Ene(max_diago)
!------ working parameters --------------------------------


      real (kind=Rkind) :: phi_j(para_propa%para_Davidson%L_filter)
      real (kind=Rkind) :: Delta_Lambda,Lambda
      integer       :: i,j,n,jorth

      TYPE (param_psi)             :: q1 ! intial vector (random ?)
      TYPE (param_psi)             :: z(para_propa%para_Davidson%L_filter) ! intial vector (random ?)
      TYPE (param_psi)             :: Hz(para_propa%para_Davidson%L_filter) ! intial vector (random ?)

      TYPE (param_psi)             :: Tnq1(0:para_propa%para_Davidson%M_filter) ! intial vector (random ?)
      TYPE (param_psi)             :: g,w1 ! working vector

      real (kind=Rkind), allocatable :: H(:,:),Vec(:,:)

      logical                   :: convergeEne(max_diago),convergeResi(max_diago)
      real (kind=Rkind)         :: non_hermitic,epsi,auTOcm_inv
      real (kind=Rkind)         :: RS,a,max_Sii,max_Sij
      real (kind=Rkind)         :: tab_norm2g(max_diago)
      complex (kind=Rkind)      :: Overlap

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_FilterDiagonalization_v0'
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      para_propa%para_Davidson%E0_filter =                               &
                    para_propa%para_Davidson%E0_filter + para_propa%Hmin
      nb_diago = para_propa%para_Davidson%L_filter
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Hmin,Hmax     : ',para_propa%Hmin,para_propa%Hmax
        write(out_unitp,*) 'E0_filter (ua): ',para_propa%para_Davidson%E0_filter
        write(out_unitp,*) 'L_filter      : ',para_propa%para_Davidson%L_filter
        write(out_unitp,*) 'M_filter      : ',para_propa%para_Davidson%M_filter
        write(out_unitp,*) 'nb_diago      : ',nb_diago
        write(out_unitp,*) 'max_diago     : ',max_diago

      END IF
!-----------------------------------------------------------
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')
!-----------------------------------------------------------

      write(out_unitp,*) ' Propagation: ',para_propa%name_WPpropa

!     - scaling of H ---------------------------------------
      para_propa%para_poly%deltaE = para_propa%Hmax - para_propa%Hmin
      para_propa%para_poly%E0     = para_propa%Hmin + HALF * para_propa%para_poly%deltaE
      para_propa%para_poly%Esc = HALF * para_propa%para_poly%deltaE

      para_H%scaled = .TRUE.
      para_H%E0     = para_propa%para_poly%E0
      para_H%Esc    = para_propa%para_poly%Esc

      write(out_unitp,*) ' deltaE,E0,Esc: ',para_propa%para_poly%deltaE,&
                                            para_H%E0,para_H%Esc
!-----------------------------------------------------------

      !- vector initialization:  q1, z(:), Tnq1(:), psi(:) -----------
      CALL init_psi(q1,para_H,para_H%cplx)
      CALL init_psi(g,para_H,para_H%cplx)

      CALL alloc_psi(q1)
      IF (q1%cplx) THEN
        DO i=1,q1%nb_tot
          CALL random_number(a)
          q1%CvecB(i) = cmplx(a,ZERO,kind=Rkind)
        END DO
      ELSE
        DO i=1,q1%nb_tot
          CALL random_number(a)
          q1%RvecB(i) = a
        END DO
      END IF
      CALL Set_symab_OF_psiBasisRep(q1,para_propa%para_Davidson%symab)
      CALL renorm_psi(q1,BasisRep=.TRUE.)

      DO j=1,nb_diago
        CALL init_psi(psi(j),para_H,para_H%cplx)
        CALL init_psi(z(j),para_H,para_H%cplx)
        CALL init_psi(Hz(j),para_H,para_H%cplx)
      END DO
      DO n=0,para_propa%para_Davidson%M_filter
        CALL init_psi(Tnq1(n),para_H,para_H%cplx)
      END DO
      !- vector initialization:  q1, z(:), Tnq1(:), psi(:) -----------



     para_propa%para_Davidson%W_filter =                                &
                    real(para_propa%para_Davidson%L_filter,kind=Rkind)* &
                                        para_propa%para_poly%deltaE /   &
          (1.2_Rkind*real(para_propa%para_Davidson%M_filter,kind=Rkind))
     para_propa%para_Davidson%W_filter = 200._Rkind / auTOcm_inv ! 200 cm-1
     write(out_unitp,*) 'W_filter (ua)',para_propa%para_Davidson%W_filter
     write(out_unitp,*) 'W_filter (cm-1)',para_propa%para_Davidson%W_filter*&
                                                               auTOcm_inv


     Delta_Lambda = para_propa%para_Davidson%W_filter /                 &
                    real(para_propa%para_Davidson%L_filter-1,kind=Rkind)
     write(out_unitp,*) ' Delta_Lambda (ua)  : ',Delta_Lambda
     write(out_unitp,*) ' Delta_Lambda (cm-1): ',Delta_Lambda*auTOcm_inv

     Lambda = para_propa%para_Davidson%E0_filter -                      &
                                 HALF*para_propa%para_Davidson%W_filter
     IF (Lambda < para_propa%Hmin) Lambda = para_propa%Hmin
     DO j=1,para_propa%para_Davidson%L_filter
       write(out_unitp,*) 'j,Lambda (cm-1)',j,Lambda * auTOcm_inv
       phi_j(j) = acos((Lambda-para_H%E0)/para_H%Esc)
       Lambda = Lambda + Delta_Lambda
     END DO
     IF (debug) write(out_unitp,*) 'phi_j(:): done'
     CALL flush_perso(out_unitp)

!----------------------------------------------------------



      !- chebychev recursion -------------------------------------------
      Tnq1(0)  = q1                                                ! => q1
      CALL sub_OpPsi(q1,Tnq1(1),para_H)                            ! => H.q1
      CALL sub_scaledOpPsi(q1,Tnq1(1),para_H%E0,para_H%Esc)        ! scaling

      DO n=2,para_propa%para_Davidson%M_filter
        CALL sub_OpPsi(Tnq1(n-1),w1,para_H)                        ! => Ti-1.q1
        CALL sub_scaledOpPsi(Tnq1(n-1),w1,para_H%E0,para_H%Esc)    ! scaling
        Tnq1(n) = TWO * w1
        Tnq1(n) = Tnq1(n) - Tnq1(n-2)
      END DO
     IF (debug) write(out_unitp,*) 'chebychev recursion: done'
     CALL flush_perso(out_unitp)

      !- chebychev recursion -------------------------------------------

      !- z vectors --------------------------------------------
      DO j=1,para_propa%para_Davidson%L_filter
        z(j) = Tnq1(0)
        DO n=1,para_propa%para_Davidson%M_filter
          z(j) = z(j) + TWO*cos(real(n,kind=rkind)*phi_j(j)) * Tnq1(n)
        END DO
      END DO
     IF (debug) write(out_unitp,*) 'z vectors: done'
     CALL flush_perso(out_unitp)
      !- z vectors -------------------------------------------


      !- z vectors (orthonormalized) -----------------------------------
      jorth = 0
      DO j=1,para_propa%para_Davidson%L_filter
        CALL renorm_psi(z(j))
        write(out_unitp,*) 'j,norm',j,z(j)%norm2

        DO i=1,jorth
          CALL norm2_psi(z(i))
          !write(out_unitp,*) '    i,norm',i,z(i)%norm2
          CALL Overlap_psi1_psi2(Overlap,z(j),z(i))
          !write(out_unitp,*) '    j,i,S(j,i)',j,i,real(Overlap,kind=Rkind)

          z(j) = z(j) - z(i) * real(Overlap,kind=Rkind)
        END DO
        DO i=1,jorth
          CALL Overlap_psi1_psi2(Overlap,z(j),z(i))
          RS = real(Overlap,kind=Rkind)
          z(j) = z(j) - z(i) * RS
        END DO
        CALL norm2_psi(z(j))
        !write(out_unitp,*) 'j,norm',j,z(j)%norm2
        IF (z(j)%norm2 > ONETENTH**6) THEN
          jorth = jorth + 1
          z(j) = z(j) * (ONE/sqrt(z(j)%norm2))
          IF (jorth < j) z(jorth) = z(j)
        ELSE
          z(j) = ZERO
        END IF
        write(out_unitp,*) 'jorth',jorth
        write(out_unitp,*)

      END DO
      nb_diago = jorth
     IF (debug) write(out_unitp,*) 'ortho z vectors: done',jorth
     CALL flush_perso(out_unitp)
      !- z vectors -------------------------------------------

      !- diagonalization -----------------------------------
      CALL alloc_NParray(H,  [nb_diago,nb_diago],'H',  name_sub)
      CALL alloc_NParray(Vec,[nb_diago,nb_diago],'Vec',name_sub)

      DO j=1,nb_diago
      DO i=1,nb_diago
        CALL Overlap_psi1_psi2(Overlap,z(i),z(j))
        H(i,j) = real(Overlap,kind=Rkind)
      END DO
      END DO
      !IF (debug) write(out_unitp,*) 'H matrix: done'
      CALL flush_perso(out_unitp)
      CALL sub_ana_S(H,nb_diago,max_Sii,max_Sij,.TRUE.)

      !IF (debug) CALL Write_Mat(H,out_unitp,5)



      DO j=1,nb_diago
        !H.z(j) calc
        CALL sub_OpPsi(z(j),Hz(j),para_H)                            ! => H.z(j)
        !CALL sub_scaledOpPsi(z(j),Hz(j),para_H%E0,para_H%Esc)    ! scaling

        DO i=1,nb_diago
          CALL Overlap_psi1_psi2(Overlap,z(i),Hz(j))
          H(i,j) = real(Overlap,kind=Rkind)
        END DO
      END DO
      !IF (debug) write(out_unitp,*) 'H matrix: done'
      CALL flush_perso(out_unitp)

      CALL sub_hermitic_H(H,nb_diago,non_hermitic,para_H%sym_Hamil)
      !IF (debug) CALL Write_Mat(H,out_unitp,5)

      IF (non_hermitic > FOUR*ONETENTH**4) THEN
        If(MPI_id==0) write(out_unitp,*) 'WARNING: non_hermitic is BIG'
        If(MPI_id==0) write(out_unitp,31) non_hermitic
 31     format(' Hamiltonien: ',f16.12,' au')
      ELSE
        If(MPI_id==0) write(out_unitp,51) non_hermitic*auTOcm_inv
 51     format(' Hamiltonien: ',f16.12,' cm-1')
      END IF
      epsi = max(para_propa%para_Davidson%conv_resi,                    &
                   TEN**para_propa%para_Davidson%conv_hermitian *       &
                                               non_hermitic)

        IF (para_H%sym_Hamil) THEN
          CALL diagonalization(H,Ene(1:nb_diago),Vec,nb_diago,3,1,.FALSE.)
        ELSE
          CALL diagonalization(H,Ene(1:nb_diago),Vec,nb_diago,4,1,.FALSE.)
        END IF
      !write(out_unitp,*) Ene(1:nb_diago)*auTOcm_inv

       write(out_unitp,21) Ene(1:nb_diago)*auTOcm_inv
       write(out_unitp,21) (Ene(1:nb_diago)-para_propa%Hmin)*auTOcm_inv

 21   format(' Filter: ',50(1x,f18.4))


      !----------------------------------------------------------
      !- residual vector ---------------------------
      IF (debug) write(out_unitp,*) 'residual'
      IF (debug) CALL flush_perso(out_unitp)
      DO j=1,nb_diago
          g = ZERO
          psi(j) = ZERO
          DO i=1,nb_diago
            psi(j) = psi(j) + Vec(i,j) * z(i)
            w1 = Hz(i) - z(i) * Ene(j)
            g = g + w1 * Vec(i,j)
          END DO
          CALL norm2_psi(g)
          tab_norm2g(j) = sqrt(g%norm2)
          convergeResi(j) = tab_norm2g(j) < epsi

      END DO
      write(out_unitp,41) 'tab_norm2g          ',tab_norm2g(1:nb_diago)
      write(out_unitp,42)  'convergenceResi(:): ',convergeResi(1:nb_diago)
 41   format(a,100(1x,e9.2))
 42   format(a,100(1x,l9))
      IF (debug) write(out_unitp,*) 'residual: done'
      write(out_unitp,*)  ' number of converged levels: ',count(convergeResi(1:nb_diago))

      IF (debug) CALL flush_perso(out_unitp)
      !- residual vector and convergence ------------------------
      !----------------------------------------------------------

      DO j=1,nb_diago
         psi(j)%CAvOp    = Ene(j)
         psi(j)%IndAvOp  = para_H%n_Op  ! it should be 0
         psi(j)%convAvOp = convergeResi(j)
      END DO


!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------


      END SUBROUTINE sub_FilterDiagonalization_v0

     SUBROUTINE sub_chebychev_recursion(Tnq1,q1,m0,mf,para_Op)
      USE mod_system
      USE mod_psi,     ONLY : param_psi,dealloc_psi
      USE mod_Op
      IMPLICIT NONE

      ! Operator (Hamiltonian)
      TYPE (param_Op)   :: para_Op

      !-----vector ----------------------------
      TYPE (param_psi), intent(inout)    :: Tnq1(0:)   ! vector for the chebychev recursion
      TYPE (param_psi), intent(in)       :: q1        ! initial vector
      integer, intent(in)                :: m0,mf     ! recursion form m0 to mf

!------ working parameters --------------------------------

      integer                      :: n,m0_loc
      TYPE (param_psi)             :: w1 ! working vector

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_chebychev_recursion'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'm0,mf              : ',m0,mf
        write(out_unitp,*) 'ubound,lbound Tnq1 : ',lbound(Tnq1,dim=1),ubound(Tnq1,dim=1)
      END IF
!-----------------------------------------------------------

     IF (lbound(Tnq1,dim=1) /= 0 .OR. lbound(Tnq1,dim=1) > m0 .OR.      &
                                         ubound(Tnq1,dim=1) < mf) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' incompatible parameters:'
        write(out_unitp,*) 'm0,mf              : ',m0,mf
        write(out_unitp,*) 'lbound,ubound Tnq1 : ',lbound(Tnq1,dim=1),ubound(Tnq1,dim=1)
        write(out_unitp,*) 'lbound(Tnq1) must be = 0 and <= m0'
        write(out_unitp,*) 'ubound(Tnq1) must be >= mf'
        STOP 'chebychev recursion'
     END IF

     m0_loc = m0
     !m0_loc = 0 ! for debugging

     !- chebychev recursion -------------------------------------------
     IF (print_level>-1) write(out_unitp,'(a)',ADVANCE='no') 'cheby rec (%): '
     CALL flush_perso(out_unitp)

     IF (m0_loc == 0) THEN
        Tnq1(0)  = q1                                              ! => T0=q1
     END IF
     IF (m0_loc <= 1) THEN
        CALL sub_OpPsi(q1,Tnq1(1),para_Op)                         ! => T1=H.q1
        CALL sub_scaledOpPsi(q1,Tnq1(1),para_Op%E0,para_Op%Esc)    ! scaling
     END IF

     DO n=max(2,m0_loc),mf
       CALL sub_OpPsi(Tnq1(n-1),w1,para_Op)                         ! => H.Tn-1
       CALL sub_scaledOpPsi(Tnq1(n-1),w1,para_Op%E0,para_Op%Esc)    ! scaling
       Tnq1(n) = TWO * w1
       Tnq1(n) = Tnq1(n) - Tnq1(n-2)

       IF (mod(n-m0,max(1,int((mf-m0)/10))) == 0 .AND. print_level>-1) THEN
         write(out_unitp,'(a,i3)',ADVANCE='no') ' -',(n-m0)*100/(mf-m0)
         CALL flush_perso(out_unitp)
       END IF
     END DO
     IF (print_level>-1) write(out_unitp,'(a)',ADVANCE='yes') ' end'
     CALL flush_perso(out_unitp)
     CALL dealloc_psi(w1)

      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF

     END SUBROUTINE sub_chebychev_recursion

     SUBROUTINE sub_Z_vectors(z,Tnq1,m0,mf,para_H,para_propa)
      USE mod_system
      USE mod_psi,    ONLY : param_psi,sub_Lowdin,sub_Schmidt
      USE mod_Op
      USE mod_propa
      IMPLICIT NONE


      ! Operator (Hamiltonian)
      TYPE (param_Op)              :: para_H
      TYPE (param_propa)           :: para_propa

      !-----vector ----------------------------
      TYPE (param_psi), intent(inout)    :: z(:)      ! z vectors
      TYPE (param_psi), intent(inout)    :: Tnq1(0:)  ! vectors for the chebychev recursion
      integer, intent(in)                :: m0,mf     ! recursion form m0 to mf

!------ working parameters --------------------------------

      integer                      :: j,n,m0_loc
      real (kind=Rkind) :: phi_j,Delta_Lambda,Lambda

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_Z_vectors'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'm0,mf              : ',m0,mf
        write(out_unitp,*) 'ubound,lbound Tnq1 : ',lbound(Tnq1,dim=1),ubound(Tnq1,dim=1)
      END IF
!-----------------------------------------------------------

     IF (lbound(Tnq1,dim=1) /= 0 .OR. lbound(Tnq1,dim=1) > m0 .OR.      &
                                         ubound(Tnq1,dim=1) < mf) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' incompatible parameters:'
        write(out_unitp,*) 'm0,mf              : ',m0,mf
        write(out_unitp,*) 'lbound,ubound Tnq1 : ',lbound(Tnq1,dim=1),ubound(Tnq1,dim=1)
        write(out_unitp,*) 'lbound(Tnq1) must be = 0 and <= m0'
        write(out_unitp,*) 'ubound(Tnq1) must be >= mf'
        STOP 'chebychev recursion'
     END IF

     Delta_Lambda = para_propa%para_Davidson%W_filter /                 &
                    real(para_propa%para_Davidson%L_filter-1,kind=Rkind)

     ! we cannot use m0, because the z-vectors are orthonormalized
     m0_loc = 0

     !- z vectors --------------------------------------------
     IF (m0_loc == 0) THEN
       DO j=1,para_propa%para_Davidson%L_filter
         z(j) = Tnq1(0)
       END DO
     END IF
     DO j=1,para_propa%para_Davidson%L_filter
       Lambda = para_propa%para_Davidson%LambdaMin +                    &
                                      real(j-1,kind=Rkind)*Delta_Lambda
       phi_j = acos((Lambda-para_H%E0)/para_H%Esc)
       DO n=max(1,m0_loc),mf
         z(j) = z(j) + TWO*cos(real(n,kind=rkind)*phi_j) * Tnq1(n)
       END DO
     END DO

     !- z vectors (orthonormalized) -----------------------------------
     CALL sub_Lowdin(z,para_propa%para_Davidson%L_filter)
     CALL sub_Schmidt(z,para_propa%para_Davidson%L_filter) ! to have a better identity

     IF (debug) THEN
       write(out_unitp,*) 'END ',name_sub
     END IF

     END SUBROUTINE sub_Z_vectors

     SUBROUTINE sub_Z_vectors_withf(z,Tnq1,mf,para_H,para_propa,f)
      USE mod_system
      USE mod_psi,    ONLY : param_psi,sub_Lowdin,sub_Schmidt
      USE mod_Op
      USE mod_propa
      IMPLICIT NONE


      ! Operator (Hamiltonian)
      TYPE (param_Op)              :: para_H
      TYPE (param_propa)           :: para_propa

      !-----vector ----------------------------
      TYPE (param_psi), intent(inout)    :: z(:)      ! z vectors
      TYPE (param_psi), intent(inout)    :: Tnq1(0:)  ! vectors for the chebychev recursion
      integer, intent(in)                :: mf     ! recursion form m0 to mf
      real (kind=Rkind), allocatable     :: f(:,:)

!------ working parameters --------------------------------

      integer                      :: l,k

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_Z_vectors_withf'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'mf              : ',mf
        write(out_unitp,*) 'ubound,lbound Tnq1 : ',lbound(Tnq1,dim=1),ubound(Tnq1,dim=1)
      END IF
!-----------------------------------------------------------

     IF (lbound(Tnq1,dim=1) /= 0 .OR. ubound(Tnq1,dim=1) < mf) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' incompatible parameters:'
        write(out_unitp,*) 'mf              : ',mf
        write(out_unitp,*) 'lbound,ubound Tnq1 : ',lbound(Tnq1,dim=1),ubound(Tnq1,dim=1)
        write(out_unitp,*) 'lbound(Tnq1) must be = 0'
        write(out_unitp,*) 'ubound(Tnq1) must be >= mf'
        STOP 'chebychev recursion'
     END IF

     !- z vectors --------------------------------------------
     DO l=1,para_propa%para_Davidson%L_filter
       z(l) = Tnq1(0)
       z(l) = ZERO
       DO k=0,mf
         z(l) = z(l) + f(k,l) * Tnq1(k)
       END DO
     END DO

     !- z vectors (orthonormalized) -----------------------------------
     CALL sub_Lowdin(z,para_propa%para_Davidson%L_filter)
     CALL sub_Schmidt(z,para_propa%para_Davidson%L_filter) ! to have a better identity

     IF (debug) THEN
       write(out_unitp,*) 'END ',name_sub
     END IF

     END SUBROUTINE sub_Z_vectors_withf

     SUBROUTINE sub_newZ_vectors_withf(z,Tnq1,nb,para_H,para_propa,f)
      USE mod_system
      USE mod_psi,    ONLY : param_psi,sub_Lowdin,sub_Schmidt
      USE mod_Op
      USE mod_propa
      IMPLICIT NONE

      ! Operator (Hamiltonian)
      TYPE (param_Op)              :: para_H
      TYPE (param_propa)           :: para_propa

      !-----vector ----------------------------
      TYPE (param_psi), intent(inout)    :: z(:)     ! z vectors
      TYPE (param_psi), intent(inout)    :: Tnq1(:)  ! vectors for the chebychev recursion
      integer, intent(in)                :: nb       ! recursion form 1 to nb
      real (kind=Rkind), allocatable     :: f(:,:)

!------ working parameters --------------------------------

      integer                      :: l,k

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_newZ_vectors_withf'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb                 : ',nb
        write(out_unitp,*) 'ubound,lbound Tnq1 : ',lbound(Tnq1,dim=1),ubound(Tnq1,dim=1)
        DO l=1,para_propa%para_Davidson%L_filter
          write(out_unitp,*) 'f(:,l)',l,f(:,l)
        END DO
      END IF
!-----------------------------------------------------------

     IF (lbound(Tnq1,dim=1) /= 1 .OR. ubound(Tnq1,dim=1) < nb) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' incompatible parameters:'
        write(out_unitp,*) 'nb                 : ',nb
        write(out_unitp,*) 'lbound,ubound Tnq1 : ',lbound(Tnq1,dim=1),ubound(Tnq1,dim=1)
        write(out_unitp,*) 'lbound(Tnq1) must be = 1'
        write(out_unitp,*) 'ubound(Tnq1) must be >= nb'
        STOP 'chebychev recursion'
     END IF

     !- z vectors --------------------------------------------
     DO l=1,para_propa%para_Davidson%L_filter
       z(l) = Tnq1(1) ! for the allocation of z(l)
       z(l) = ZERO
       DO k=1,nb
         z(l) = z(l) + f(k,l) * Tnq1(k)
       END DO
     END DO

     !- z vectors (orthonormalized) -----------------------------------
     CALL sub_Lowdin(z,para_propa%para_Davidson%L_filter)
     CALL sub_Schmidt(z,para_propa%para_Davidson%L_filter) ! to have a better identity

     IF (debug) THEN
       write(out_unitp,*) 'END ',name_sub
     END IF

     END SUBROUTINE sub_newZ_vectors_withf

     FUNCTION f_filter_gauss(E,El,sigma)
      USE mod_system
      IMPLICIT NONE

      real (kind=Rkind), intent(in) :: E,El,sigma
      real (kind=Rkind)             :: f_filter_gauss

      !----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='f_filter_gauss'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      !-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'E,El,sigma              : ',E,El,sigma
      END IF
      !-----------------------------------------------------------

      f_filter_gauss = exp(-((E-El)/sigma)**2)

      !-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'DE/sigma                : ',(E-El)/sigma
        write(out_unitp,*) 'f_filter                : ',f_filter_gauss
        write(out_unitp,*) 'END ',name_sub
      END IF
      !-----------------------------------------------------------

     END FUNCTION f_filter_gauss

     FUNCTION KDF(j,nb) ! Kenerl damping factor
      USE mod_system
      IMPLICIT NONE

      integer, intent(in) :: j,nb
      real (kind=Rkind)             :: KDF

      real (kind=Rkind)             :: r,r0

      !----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='KDF'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      !-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'j              : ',j
      END IF
      !-----------------------------------------------------------

      !KDF = ONE
      !RETURN
      IF (j < 2) THEN
        KDF = ONE
      ELSE
        r0 = pi/real(nb-1,kind=Rkind)
        r = pi*real(j-1,kind=Rkind)/real(nb-1,kind=Rkind)
        KDF = (cos(r)*real(nb-j,kind=Rkind)+sin(r)*sin(r0)/cos(r0))/real(nb-1,kind=Rkind)
        !r  = pi*real(j-1,kind=Rkind)/real(nb,kind=Rkind)
        !KDF = (sin(r)/r)**2

      END IF

      !-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'KDF                : ',KDF
        write(out_unitp,*) 'END ',name_sub
      END IF
      !-----------------------------------------------------------

     END FUNCTION KDF

  FUNCTION f_filter(E,filter)
   USE mod_system
   IMPLICIT NONE

   real (kind=Rkind), intent(in)   :: E
   TYPE (param_filter), intent(in) :: filter
   real (kind=Rkind)               :: f_filter

   !----- for debuging --------------------------------------------------
   integer :: err_mem,memory
   character (len=*), parameter :: name_sub='f_filter'
   logical, parameter :: debug=.FALSE.
   !logical, parameter :: debug=.TRUE.
   !-----------------------------------------------------------
   IF (debug) THEN
     write(out_unitp,*) 'BEGINNING ',name_sub
     write(out_unitp,*) 'filter              : ',filter
   END IF
   !-----------------------------------------------------------

   SELECT CASE (filter%filter_type)
   CASE (1) ! rectangular

     IF (E >= filter%A .AND. E <= filter%B) THEN
       f_filter = ONE
     ELSE
       f_filter = ZERO
     END IF

   CASE (2) ! gaussian

     f_filter = exp(-((E-filter%E0)/filter%sigma)**2)

   CASE (3) ! smooth rectangular

     f_filter = HALF*( tanh(filter%beta* (E-filter%A)/(filter%B-filter%A) )-&
                       tanh(filter%beta* (E-filter%B)/(filter%B-filter%A) ))

   CASE (4) ! sine in [A:B]

     IF (E >= filter%A .AND. E <= filter%B) THEN
       !f_filter = sin( (E-filter%A)*pi/(filter%B-filter%A) )
       f_filter = HALF*(ONE+cos((E-filter%E0)*TWO*pi/(filter%B-filter%A)))
     ELSE
       f_filter = ZERO
     END IF

   CASE DEFAULT ! rectangular

     IF (E >= filter%A .AND. E <= filter%B) THEN
       f_filter = ONE
     ELSE
       f_filter = ZERO
     END IF

   END SELECT

   !-----------------------------------------------------------
   IF (debug) THEN
     write(out_unitp,*) 'f_filter                : ',f_filter
     write(out_unitp,*) 'END ',name_sub
   END IF
   !-----------------------------------------------------------

  END FUNCTION f_filter
  SUBROUTINE Set_filter(filter,filter_type,A,B,beta,sigma,E0)
   USE mod_system
   IMPLICIT NONE

   integer, intent(in) :: filter_type

   real (kind=Rkind), intent(in),optional :: A,B
   real (kind=Rkind), intent(in),optional :: beta

   real (kind=Rkind), intent(in),optional :: sigma,E0


   TYPE (param_filter), intent(inout) :: filter

   !----- for debuging --------------------------------------------------
   integer :: err_mem,memory
   character (len=*), parameter :: name_sub='Set_filter'
   logical, parameter :: debug=.FALSE.
   !logical, parameter :: debug=.TRUE.
   !-----------------------------------------------------------
   IF (debug) THEN
     write(out_unitp,*) 'BEGINNING ',name_sub
     write(out_unitp,*) 'filter_type         : ',filter_type
     write(out_unitp,*) 'filter              : ',filter
   END IF
   !-----------------------------------------------------------

   filter%filter_type = filter_type

   SELECT CASE (filter%filter_type)
   CASE (1) ! rectangular

     IF (present(A) .AND. present(B)) THEN
       filter%A = A
       filter%B = B
     ELSE
       write(out_unitp,*) 'ERROR in ',name_sub
       write(out_unitp,*) ' A or B are not present',present(A),present(B)
       STOP 'in set_filter'
     END IF

   CASE (2) ! gaussian

     IF (present(E0) .AND. present(sigma)) THEN
       filter%E0    = E0
       filter%sigma = sigma
     ELSE
       write(out_unitp,*) 'ERROR in ',name_sub
       write(out_unitp,*) ' E0 or sigma are not present',present(E0),present(sigma)
       STOP 'in set_filter'
     END IF

   CASE (3) ! smooth rectangular

     IF (present(A) .AND. present(B)) THEN
       filter%A = A
       filter%B = B
     ELSE
       write(out_unitp,*) 'ERROR in ',name_sub
       write(out_unitp,*) ' A or B are not present',present(A),present(B)
       STOP 'in set_filter'
     END IF
     IF (present(beta)) filter%beta = beta

   CASE (4) ! sine in [A:B]

     IF (present(A) .AND. present(B)) THEN
       filter%A  = A
       filter%B  = B
       filter%E0 = (A+B)/TWO
     ELSE
       write(out_unitp,*) 'ERROR in ',name_sub
       write(out_unitp,*) ' A or B are not present',present(A),present(B)
       STOP 'in set_filter'
     END IF

   CASE DEFAULT ! rectangular
     write(out_unitp,*) 'ERROR in ',name_sub
     write(out_unitp,*) ' Unknown filter_type: ',filter_type
     STOP 'in set_filter'

   END SELECT

   !-----------------------------------------------------------
   IF (debug) THEN
     write(out_unitp,*) 'filter                : ',filter
     write(out_unitp,*) 'END ',name_sub
   END IF
   !-----------------------------------------------------------

  END SUBROUTINE Set_filter


END MODULE mod_Filter
