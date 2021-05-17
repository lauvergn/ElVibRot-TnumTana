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
!=====================================================================
! POGridRep_basis
!=====================================================================
MODULE BasisMakeGrid
USE mod_system
IMPLICIT NONE

  TYPE param_SimulatedAnnealing

  integer           :: nb_mc_tot           =  100000
  integer           :: nb_mc_partial       =  100

  integer           :: TempInit_type       =  1
  real (kind=Rkind) :: Tmax                = -ONE
  real (kind=Rkind) :: Tmin                =  ONETENTH**7
  real (kind=Rkind) :: DeltaT              =  ZERO

  real (kind=Rkind) :: RangeScal           =  0.8_Rkind
  real (kind=Rkind) :: RangeScalInit       =  1._Rkind

  logical           :: With_RangeInit      = .FALSE.
  real (kind=Rkind) :: RangeInit           = 1._Rkind

  integer           :: TempScheduling_type =  2 ! 1: linear, 2: geometrical ...
  real (kind=Rkind) :: ExpCoolParam        =  0.95_Rkind

  logical           :: ResetTemp           = .TRUE.
  real (kind=Rkind) :: ResetTempScal       =  ONE/THREE

  integer           :: Restart_Opt         =  0


  END TYPE param_SimulatedAnnealing
  TYPE param_Grid_FOR_SA

  integer           :: type_weight         =  0

  integer           :: type_grid           =  0
  logical           :: ReOriented_grid     = .TRUE.
  logical           :: ReCentered_grid     = .TRUE.


  integer           :: nb_mc_partial       =  100


  integer           :: TempInit_type       =  1
  real (kind=Rkind) :: Tmax                = -ONE
  real (kind=Rkind) :: Tmin                =  ONETENTH**7
  real (kind=Rkind) :: DeltaT              =  ZERO

  real (kind=Rkind) :: RangeScal           =  0.8_Rkind
  real (kind=Rkind) :: RangeScalInit       =  1._Rkind

  logical           :: With_RangeInit      = .FALSE.
  real (kind=Rkind) :: RangeInit           = 1._Rkind

  integer           :: TempScheduling_type =  2 ! 1: linear, 2: geometrical ...
  real (kind=Rkind) :: ExpCoolParam        =  0.95_Rkind


  real (kind=Rkind) :: ResetTempScal       =  ONE/THREE

  integer           :: Restart_Opt         =  0


  END TYPE param_Grid_FOR_SA

CONTAINS

      SUBROUTINE Read_param_SimulatedAnnealing(para_SimulatedAnnealing)
      TYPE (param_SimulatedAnnealing), intent(inout) :: para_SimulatedAnnealing

        integer :: nb_mc_tot     = 1000
        integer :: nb_mc_partial = 100

        real (kind=Rkind) :: Tmax          = -ONE
        real (kind=Rkind) :: Tmin          =  ONETENTH**7
        real (kind=Rkind) :: DeltaT        =  ZERO
        real (kind=Rkind) :: ResetTempScal = ONE/THREE


        real (kind=Rkind) :: ExpCoolParam = 0.95_Rkind

        real (kind=Rkind) :: RangeScal     = 0.8_Rkind
        real (kind=Rkind) :: RangeScalInit = 1._Rkind
        logical           :: With_RangeInit = .FALSE.
        real (kind=Rkind) :: RangeInit     = 1._Rkind

        logical :: ResetTemp               = .TRUE.
        integer :: TempScheduling_type     = 2 ! 1: linear, 2: geometrical ...
        integer :: TempInit_type           = 1

        integer :: Restart_Opt     = 0

        integer :: err_io
        NAMELIST /SimulatedAnnealing/nb_mc_tot,nb_mc_partial,           &
                                        Tmax,Tmin,DeltaT,TempInit_type, &
                                               RangeScal,RangeScalInit, &
                                              With_RangeInit,RangeInit, &
                                      TempScheduling_type,ExpCoolParam, &
                                    ResetTemp,ResetTempScal,Restart_Opt

        nb_mc_tot           =  1000
        nb_mc_partial       =  100

        TempInit_type       =  1
        Tmax                = -ONE
        Tmin                =  ONETENTH**7
        DeltaT              =  ZERO

        RangeScal           =  0.8_Rkind
        RangeScalInit       =  1._Rkind

        TempScheduling_type =  2 ! 1: linear, 2: geometrical ...
        ExpCoolParam        =  0.95_Rkind

        ResetTemp           = .TRUE.
        ResetTempScal       =  ONE/THREE

        Restart_Opt         =  0

        read(in_unitp,SimulatedAnnealing,IOSTAT=err_io)
        IF (err_io /= 0) THEN
           write(out_unitp,*) ' WARNING in Read_param_SimulatedAnnealing'
           write(out_unitp,*) '  while reading the "SimulatedAnnealing" namelist'
           write(out_unitp,*) ' end of file or end of record'
           write(out_unitp,*) ' Check your data !!'
           STOP
        END IF
        IF (print_level > 1) write(out_unitp,SimulatedAnnealing)

        para_SimulatedAnnealing%nb_mc_tot           =  nb_mc_tot
        para_SimulatedAnnealing%nb_mc_partial       =  nb_mc_partial

        para_SimulatedAnnealing%TempInit_type       =  TempInit_type
        para_SimulatedAnnealing%Tmax                =  Tmax
        para_SimulatedAnnealing%Tmin                =  Tmin
        para_SimulatedAnnealing%DeltaT              =  DeltaT

        para_SimulatedAnnealing%With_RangeInit      =  With_RangeInit
        para_SimulatedAnnealing%RangeInit           =  RangeInit
        IF (With_RangeInit) RangeScalInit           =  ONE

        para_SimulatedAnnealing%RangeScal           =  RangeScal
        para_SimulatedAnnealing%RangeScalInit       =  RangeScalInit



        para_SimulatedAnnealing%TempScheduling_type =  TempScheduling_type
        para_SimulatedAnnealing%ExpCoolParam        =  ExpCoolParam

        para_SimulatedAnnealing%ResetTemp           =  ResetTemp
        para_SimulatedAnnealing%ResetTempScal       =  ResetTempScal

        para_SimulatedAnnealing%Restart_Opt         =  Restart_Opt

      END SUBROUTINE Read_param_SimulatedAnnealing

      SUBROUTINE Write_param_SimulatedAnnealing(para_SimulatedAnnealing)
      TYPE (param_SimulatedAnnealing), intent(in)   :: para_SimulatedAnnealing

      write(out_unitp,*) '  WRITE param_SimulatedAnnealing'
      write(out_unitp,*)
      write(out_unitp,*) '  nb_mc_tot          ',para_SimulatedAnnealing%nb_mc_tot
      write(out_unitp,*) '  nb_mc_partial      ',para_SimulatedAnnealing%nb_mc_partial
      write(out_unitp,*)
      write(out_unitp,*) '  TempInit_type      ',para_SimulatedAnnealing%TempInit_type
      write(out_unitp,*) '  Tmax               ',para_SimulatedAnnealing%Tmax
      write(out_unitp,*) '  Tmin               ',para_SimulatedAnnealing%Tmin
      write(out_unitp,*) '  DeltaT             ',para_SimulatedAnnealing%DeltaT
      write(out_unitp,*)
      write(out_unitp,*) '  RangeScal          ',para_SimulatedAnnealing%RangeScal
      write(out_unitp,*) '  RangeScalInit      ',para_SimulatedAnnealing%RangeScalInit

      write(out_unitp,*) '  With_RangeInit     ',para_SimulatedAnnealing%With_RangeInit
      write(out_unitp,*) '  RangeInit          ',para_SimulatedAnnealing%RangeInit

      write(out_unitp,*)
      write(out_unitp,*) '  TempScheduling_type',para_SimulatedAnnealing%TempScheduling_type
      write(out_unitp,*) '  ExpCoolParam       ',para_SimulatedAnnealing%ExpCoolParam
      write(out_unitp,*)
      write(out_unitp,*) '  ResetTemp          ',para_SimulatedAnnealing%ResetTemp
      write(out_unitp,*) '  ResetTempScal      ',para_SimulatedAnnealing%ResetTempScal
      write(out_unitp,*)
      write(out_unitp,*) '  Restart_Opt        ',para_SimulatedAnnealing%Restart_Opt
      write(out_unitp,*)
      write(out_unitp,*) '  END WRITE param_SimulatedAnnealing'

      END SUBROUTINE write_param_SimulatedAnnealing


      SUBROUTINE POGridRep_basis(basis_POGridRep,nb0)
      USE mod_system
      USE mod_Coord_KEO
      USE mod_basis
      IMPLICIT NONE

!----- for the basis set ----------------------------------------------
      TYPE (basis), intent(inout)   :: basis_POGridRep
      integer,      intent(in)      :: nb0  ! basis function number before contraction

!----- local variables
      TYPE (basis)                   :: basis_temp
      real (kind=Rkind), allocatable :: matX_basis(:,:)
      real (kind=Rkind), allocatable :: RdiagX(:)
      real (kind=Rkind), allocatable :: RvpX(:,:)

      integer                        :: i,j,k,nb_var,nq
      integer                        :: nqc0
      real (kind=Rkind)              :: a
      logical                        :: periodic = .FALSE.
      integer                        :: iweight = 1


!-------------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'POGridRep_basis'
!---------------------------------------------------------------------
      IF (basis_POGridRep%ndim > 1 .OR. .NOT. basis_POGridRep%POGridRep) RETURN

      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) ' The Contracted basis'
        !CALL RecWrite_basis(basis_POGridRep)
      END IF
!---------------------------------------------------------------------
      ! The POGridRep basis set contains the nb contracted basis functions on
      ! the primitive grid (nq)

      nqc0 = basis_POGridRep%nqc

      CALL basis2TObasis1(basis_temp,basis_POGridRep)

      CALL alloc_NParray(matX_basis,(/ basis_temp%nb,basis_temp%nb /), &
                        "matX_basis",name_sub)
      CALL alloc_NParray(RvpX,(/ basis_temp%nb,basis_temp%nb /),       &
                        "RvpX",name_sub)
      CALL alloc_NParray(RdiagX,(/ basis_temp%nb /),"RdiagX",name_sub)

      !---------------------------------------------------------------
      ! make the matrix of Q
      IF (periodic) THEN
        CALL sub_matperioX_basis(basis_temp,matX_basis)
      ELSE
        CALL sub_matX_basis(basis_temp,matX_basis)
      END IF

      !---------------------------------------------------------------
      ! digonalization of Q
      CALL sub_diago_H(matX_basis,RdiagX,RvpX,basis_temp%nb,.TRUE.)

      IF (periodic) THEN
        RdiagX(:) = TWO * acos(RdiagX(:)) - PI
      END IF

      !---------------------------------------------------------------
      ! POGridRep basis on the primitive grid points
!     IF (debug) THEN
!       IF (allocated(basis_temp%Rvec))                             &
!               CALL dealloc_array(basis_temp%Rvec,'Rvec',name_sub)
!       CALL alloc_array(basis_temp%Rvec,                            &
!                    (/ basis_temp%nb,basis_temp%nb /),'Rvec',name_sub)
!       basis_temp%Rvec = RvpX
!       CALL sub_contraction_basis(basis_temp,.TRUE.)
!       write(out_unitp,*) 'POGridRep basis on the grid'
!       DO i=1,get_nq_FROM_basis(basis_temp)
!         write(out_unitp,*) i,basis_temp%x(:,i),basis_temp%dnRGB%d0(i,:)
!       END DO
!       CALL flush_perso(out_unitp)
!       STOP
!     END IF

      !---------------------------------------------------------------
      ! Eigenvalues of Q =>  the POGridRep grid points
      IF (debug) THEN
         write(out_unitp,*) ' Eigenvalues of Q: POGridRep grid points'
         DO i=1,basis_temp%nb
           write(out_unitp,*) i,RdiagX(i)
         END DO
         write(out_unitp,*) ' RvpX:'
         CALL Write_Mat(RvpX,out_unitp,5)
         !write(out_unitp,*) 'Id?',matmul(transpose(RvpX),RvpX) ; stop
         CALL flush_perso(out_unitp)
      END IF

      !---------------------------------------------------------------
      ! Primitive basis on the POGridRep grid points

      CALL basis2TObasis1(basis_POGridRep,basis_temp,init_only=.TRUE.)

      basis_POGridRep%nb = nb0
      CALL Set_nq_OF_basis(basis_POGridRep,nqc0)
      CALL alloc_xw_OF_basis(basis_POGridRep)

      ! need to unscaled x, because the scaling has to be done in construct_primitive_basis
      basis_POGridRep%x(1,:) = ( RdiagX(:) - basis_POGridRep%Q0(1) ) *   &
                                               basis_POGridRep%scaleQ(1)
      basis_POGridRep%xPOGridRep_done = .TRUE.

      CALL construct_primitive_basis(basis_POGridRep)
      !CALL RecWrite_basis(basis_POGridRep,write_all=.TRUE.)

      !---------------------------------------------------------------
      ! first transformation nb0 => nb with nq
      ! => contracted basis functions on the POGridRep grid points
      ! Remark: If Rvec is not allocated, we are using a primitive basis
      IF (allocated(basis_temp%Rvec)) THEN

        CALL alloc_NParray(basis_POGridRep%Rvec,                       &
                        (/ basis_POGridRep%nb,basis_POGridRep%nb /), &
                        "basis_POGridRep%Rvec",name_sub)
        basis_POGridRep%Rvec = basis_temp%Rvec
        !CALL Write_Mat(basis_temp%Rvec,out_unitp,5) ; stop
        basis_POGridRep%nbc = nqc0
        CALL sub_contraction_basis(basis_POGridRep,.TRUE.)
        IF (debug) THEN
          write(out_unitp,*) 'Eigenvectors on the POGridRep grid:'
          DO i=1,get_nq_FROM_basis(basis_POGridRep)
            write(out_unitp,*) i,basis_POGridRep%x(:,i),basis_POGridRep%dnRGB%d0(i,:)
          END DO
          CALL flush_perso(out_unitp)
        END IF
      END IF
      nq = get_nq_FROM_basis(basis_POGridRep)
      IF (iweight == 1) THEN
        ! the POGridRep weights
        DO i=1,nq
          basis_POGridRep%w(i) = ONE / sum(basis_POGridRep%dnRGB%d0(i,:)**2)
        END DO
        basis_POGridRep%wrho(:) = basis_POGridRep%w(:) * basis_POGridRep%rho(:)
      ELSE
        !---------------------------------------------------------------
        ! second transformation nbc => nbc with nqc
        ! => POGridRep basis on the POGridRep grid points
        ! enable to calculate  the POGridRep weights
        !write(out_unitp,*) 'shape Rvec',shape(basis_POGridRep%Rvec)
        !write(out_unitp,*) 'shape RvpX',shape(RvpX)
        !CALL flush_perso(out_unitp)
        IF (allocated(basis_POGridRep%Rvec))  THEN
          CALL dealloc_NParray(basis_POGridRep%Rvec,"basis_POGridRep%Rvec",name_sub)
        END IF
        CALL alloc_NParray(basis_POGridRep%Rvec,(/ nq,nq /),              &
                        "basis_POGridRep%Rvec",name_sub)
        basis_POGridRep%Rvec = RvpX
        CALL sub_contraction_basis(basis_POGridRep,.TRUE.)
        IF (debug) THEN
          write(out_unitp,*) 'POGridRep basis on the POGridRep grid'
          DO i=1,nq
            write(out_unitp,*) i,basis_POGridRep%x(:,i),basis_POGridRep%dnRGB%d0(i,:)
          END DO
          CALL flush_perso(out_unitp)
        END IF

        ! the POGridRep weights
        DO i=1,nq
          basis_POGridRep%w(i) = ONE/basis_POGridRep%dnRGB%d0(i,i)**2
        END DO
        basis_POGridRep%wrho(:) = basis_POGridRep%w(:) * basis_POGridRep%rho(:)

      END IF
      IF (debug) THEN
         write(out_unitp,*) ' Weight on POGridRep grid'
         DO i=1,nq
           write(out_unitp,*) i,basis_POGridRep%x(:,i),basis_POGridRep%w(i)
        END DO
        CALL flush_perso(out_unitp)
      END IF
      !- d1b => d1BasisRep and  d2b => d2BasisRep ------------
      CALL sub_dnGB_TO_dnBB(basis_POGridRep)

      !---------------------------------------------------------------
      !- check the overlap matrix -----------------------------
      CALL check_ortho_basis(basis_POGridRep,test_stop=.FALSE.)

      !---------------------------------------------------------------
      ! The basis is orthogonalized (it should not append !!!!)
      CALL ortho_basis_schmidt(basis_POGridRep)

      IF (debug) THEN
         write(out_unitp,*) ' basis on POGridRep grid'
         DO i=1,nq
           write(out_unitp,*) i,basis_POGridRep%x(:,i),basis_POGridRep%dnRGB%d0(i,:)
        END DO
        CALL flush_perso(out_unitp)
      END IF

      CALL dealloc_NParray(matX_basis,"matX_basis",name_sub)
      CALL dealloc_NParray(RvpX,"RvpX",name_sub)
      CALL dealloc_NParray(RdiagX,"RdiagX",name_sub)
      CALL dealloc_basis(basis_temp)

!=====================================================================
!---------------------------------------------------------------------
      IF (debug) THEN
        !write(out_unitp,*) ' POGridRep basis'
        !CALL RecWrite_basis(basis_POGridRep)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------

      END SUBROUTINE POGridRep_basis
      SUBROUTINE POGridRep2_basis(basis_POGridRep,nb0)
      USE mod_system
      USE mod_nDindex
      USE mod_Coord_KEO
      USE mod_basis
      IMPLICIT NONE

!----- for the basis set ----------------------------------------------
      TYPE (basis), intent(inout)   :: basis_POGridRep
      integer,      intent(in)      :: nb0  ! basis function number before contraction

!----- local variables
      TYPE (basis)               :: basis_temp
      real (kind=Rkind), allocatable :: matX_basis(:,:)
      real (kind=Rkind), allocatable :: RdiagX(:)
      real (kind=Rkind), allocatable :: RvpX(:,:)
      integer                    :: i,j,k
      integer                    :: nqc0,nq,nq1
      real (kind=Rkind)          :: a
      logical                    :: periodic = .FALSE.
      integer                    :: iweight = 1


!-------------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
!      logical, parameter :: debug = .FALSE.
      logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'POGridRep2_basis'
!---------------------------------------------------------------------
      IF (basis_POGridRep%ndim > 1 .OR. .NOT. basis_POGridRep%POGridRep) RETURN

      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) ' The Contracted basis'
        CALL RecWrite_basis(basis_POGridRep)
      END IF
!---------------------------------------------------------------------

         ! The POGridRep basis set contains the nb contracted basis functions on
         ! the primitive grid (nq)

         nqc0 = basis_POGridRep%nqc

         CALL basis2TObasis1(basis_temp,basis_POGridRep)

         CALL alloc_NParray(matX_basis,(/ basis_temp%nb,basis_temp%nb /), &
                                                 "matX_basis",name_sub)
         CALL alloc_NParray(RvpX,(/ basis_temp%nb,basis_temp%nb /),       &
                                                        "RvpX",name_sub)
         CALL alloc_NParray(RdiagX,(/ basis_temp%nb /),"RdiagX",name_sub)

         !---------------------------------------------------------------
         ! make the matrix of Q
         IF (periodic) THEN
           CALL sub_matperioX_basis(basis_temp,matX_basis)
         ELSE
           CALL sub_matX_basis(basis_temp,matX_basis)
         END IF
         !---------------------------------------------------------------
         ! digonalization of Q
         CALL sub_diago_H(matX_basis,RdiagX,RvpX,basis_temp%nb,.TRUE.)

         IF (periodic) THEN
           RdiagX(:) = TWO * acos(RdiagX(:)) - PI
         END IF

         !---------------------------------------------------------------
         ! POGridRep basis on the primitive grid points
!         IF (debug) THEN
!           IF (allocated(basis_temp%Rvec))                             &
!                   CALL dealloc_NParray(basis_temp%Rvec,'Rvec',name_sub)
!           CALL alloc_NParray(basis_temp%Rvec,                            &
!                      (/ basis_temp%nb,basis_temp%nb /),'Rvec',name_sub)
!           basis_temp%Rvec = RvpX
!           CALL sub_contraction_basis(basis_temp,.TRUE.)
!           write(out_unitp,*) 'POGridRep basis on the grid'
!           DO i=1,get_nq_FROM_basis(basis_temp)
!             write(out_unitp,*) i,basis_temp%x(:,i),basis_temp%dnRGB%d0(i,:)
!           END DO
!           CALL flush_perso(out_unitp)
!           STOP
!         END IF

         !---------------------------------------------------------------
         ! Eigenvalues of Q =>  the POGridRep grid points
         IF (debug) THEN
            write(out_unitp,*) ' Eigenvalues of Q: POGridRep grid points'
            DO i=1,basis_temp%nb
              write(out_unitp,*) i,RdiagX(i)
            END DO
            write(out_unitp,*) ' RvpX:'
            CALL Write_Mat(RvpX,out_unitp,5)
            !write(out_unitp,*) 'Id?',matmul(transpose(RvpX),RvpX) ; stop
            CALL flush_perso(out_unitp)
         END IF

         !---------------------------------------------------------------
         ! Primitive basis augmented with the POGridRep grid points
         CALL basis2TObasis1(basis_POGridRep,basis_temp)
         CALL dealloc_basis(basis_POGridRep)
         basis_POGridRep%nDindB = basis_temp%nDindB
         basis_POGridRep%nb     = nb0 ! nb before contraction
         nq = get_nq_FROM_basis(basis_temp)
         CALL Set_nq_OF_basis(basis_POGridRep,nq)
         CALL alloc_init_basis(basis_POGridRep)
         CALL alloc_xw_OF_basis(basis_POGridRep)
         basis_POGridRep%iQdyn  = basis_temp%iQdyn
         basis_POGridRep%A      = basis_temp%A
         basis_POGridRep%B      = basis_temp%B
         basis_POGridRep%Q0     = basis_temp%Q0
         basis_POGridRep%scaleQ = basis_temp%scaleQ

         ! grid points + POGridRep points
         nq1 = get_nq_FROM_basis(basis_POGridRep)
         basis_POGridRep%x(1,1:nq)     = basis_temp%x(1,:)
         basis_POGridRep%x(1,nq+1:nq1) = RdiagX(:)
         ! need to unscaled x, because the scaling is done in construct_basis
         basis_POGridRep%x = ( basis_POGridRep%x - basis_POGridRep%Q0(1) ) *     &
                                                  basis_POGridRep%scaleQ(1)
         !write(out_unitp,*) 'unscaled basis_POGridRep%x',basis_POGridRep%x
         basis_POGridRep%xPOGridRep_done = .TRUE.
         CALL construct_primitive_basis(basis_POGridRep)

         nq = get_nq_FROM_basis(basis_temp)
         basis_POGridRep%w(:)    = ZERO
         basis_POGridRep%w(1:nq) = basis_temp%w(:)
         basis_POGridRep%wrho(:) = basis_POGridRep%w(:) * basis_POGridRep%rho(:)
         write(out_unitp,*) 'nb,nq',basis_POGridRep%nb,get_nq_FROM_basis(basis_POGridRep)
         write(out_unitp,*) 'nbc,nqc',basis_POGridRep%nbc,basis_POGridRep%nqc

         !---------------------------------------------------------------
         ! first transformation nb0 => nb with nq
         ! => contracted basis functions on the POGridRep grid points
         ! Remark: If Rvec is not allocated, we are using a primitive basis
         IF (allocated(basis_temp%Rvec)) THEN
           CALL alloc_NParray(basis_POGridRep%Rvec,                       &
                           (/ basis_POGridRep%nb,basis_POGridRep%nb /), &
                                        "basis_POGridRep%Rvec",name_sub)
           basis_POGridRep%Rvec = basis_temp%Rvec
           !CALL Write_Mat(basis_temp%Rvec,out_unitp,5) ; stop
           CALL sub_contraction_basis(basis_POGridRep,.TRUE.)
           IF (debug) THEN
             write(out_unitp,*) 'Eigenvectors on the POGridRep grid:'
             DO i=1,get_nq_FROM_basis(basis_POGridRep)
               write(out_unitp,*) i,basis_POGridRep%x(:,i),basis_POGridRep%dnRGB%d0(i,:)
             END DO
             CALL flush_perso(out_unitp)
           END IF
         END IF
GOTO 99
         IF (iweight == 1) THEN
           ! the POGridRep weights
           DO i=1,get_nq_FROM_basis(basis_POGridRep)
             basis_POGridRep%w(i) = ZERO
             DO k=1,get_nq_FROM_basis(basis_POGridRep)
               basis_POGridRep%w(i) = basis_POGridRep%w(i) + basis_POGridRep%dnRGB%d0(i,k)**2
             END DO
             basis_POGridRep%w(i) = ONE/basis_POGridRep%w(i)
           END DO
           basis_POGridRep%wrho(:) = basis_POGridRep%w(:) * basis_POGridRep%rho(:)

         ! Ortho base POGridRep on POGridRep grid
!         DO j=1,get_nq_FROM_basis(basis_POGridRep) ! base
!         DO i=1,get_nq_FROM_basis(basis_POGridRep) ! grid
!           a = ZERO
!           DO k=1,get_nq_FROM_basis(basis_POGridRep)
!           a = a + basis_POGridRep%dnRGB%d0(i,k)*basis_POGridRep%dnRGB%d0(j,k)
!           END DO
!           a = a * basis_POGridRep%w(j)
!           write(out_unitp,*) 'u_j(xi)',j,i,a
!         END DO
!         END DO
!         STOP

         ELSE
           !---------------------------------------------------------------
           ! second transformation nbc => nbc with nqc
           ! => POGridRep basis on the POGridRep grid points
           ! enable to calculate  the POGridRep weights
           write(out_unitp,*) 'shape Rvec',shape(basis_POGridRep%Rvec)
           write(out_unitp,*) 'shape RvpX',shape(RvpX)
           CALL flush_perso(out_unitp)
           IF (allocated(basis_POGridRep%Rvec))  THEN
             CALL dealloc_NParray(basis_POGridRep%Rvec,"basis_POGridRep%Rvec",name_sub)
           END IF
           nq = get_nq_FROM_basis(basis_POGridRep)
           CALL alloc_NParray(basis_POGridRep%Rvec,(/ nq,nq /),           &
                           "basis_POGridRep%Rvec",name_sub)
           basis_POGridRep%Rvec = RvpX
           CALL sub_contraction_basis(basis_POGridRep,.TRUE.)
           IF (debug) THEN
             write(out_unitp,*) 'POGridRep basis on the POGridRep grid'
             DO i=1,get_nq_FROM_basis(basis_POGridRep)
               write(out_unitp,*) i,basis_POGridRep%x(:,i),basis_POGridRep%dnRGB%d0(i,:)
             END DO
             CALL flush_perso(out_unitp)
           END IF

           ! the POGridRep weights
           DO i=1,nq
             basis_POGridRep%w(i) = ONE/basis_POGridRep%dnRGB%d0(i,i)**2
           END DO
           basis_POGridRep%wrho(:) = basis_POGridRep%w(:) * basis_POGridRep%rho(:)

         END IF
         IF (debug) THEN
            write(out_unitp,*) ' Weight on POGridRep grid'
            DO i=1,get_nq_FROM_basis(basis_POGridRep)
              write(out_unitp,*) i,basis_POGridRep%x(:,i),basis_POGridRep%w(i)
           END DO
           CALL flush_perso(out_unitp)
         END IF
         99 CONTINUE
         !---------------------------------------------------------------
         ! The basis is orthogonalized (it should not append !!!!)
         !CALL ortho_basis(basis_POGridRep)
         CALL ortho_basis_schmidt(basis_POGridRep)

         !---------------------------------------------------------------
!        - check the overlap matrix -----------------------------
         CALL check_ortho_basis(basis_POGridRep)

         IF (debug) THEN
            write(out_unitp,*) ' basis on POGridRep grid'
            DO i=1,get_nq_FROM_basis(basis_POGridRep)
              write(out_unitp,*) i,basis_POGridRep%x(:,i),basis_POGridRep%dnRGB%d0(i,:)
           END DO
           CALL flush_perso(out_unitp)
         END IF

         CALL dealloc_NParray(matX_basis,"matX_basis",name_sub)
         CALL dealloc_NParray(RvpX,"RvpX",name_sub)
         CALL dealloc_NParray(RdiagX,"RdiagX",name_sub)
         CALL dealloc_basis(basis_temp)

!=====================================================================
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) ' POGridRep basis'
        CALL RecWrite_basis(basis_POGridRep)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------

      END SUBROUTINE POGridRep2_basis
!=====================================================================
! POGridRep_basis
!=====================================================================
      SUBROUTINE sub_make_polyorthobasis(base,Rvp)
      USE mod_system
      USE mod_dnSVM
      USE mod_basis
      IMPLICIT NONE

!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis), intent(inout)   :: base
      real (kind=Rkind), intent(in) :: Rvp(base%nb,base%nb)

!     - working variables --------------------------------------------
      integer       :: i,j,nb0,nq
      real (kind=Rkind), allocatable :: d0bc(:,:)
      real (kind=Rkind), allocatable :: d1bc(:,:,:)
      real (kind=Rkind), allocatable :: d2bc(:,:,:,:)
      real (kind=Rkind), allocatable :: d0b(:,:)

      real (kind=Rkind) :: x0,Sjj,Sii,Sij,phi0max
      real (kind=Rkind), allocatable :: Wsuzy(:)
      real (kind=Rkind), allocatable :: V(:)
      real (kind=Rkind), allocatable :: rho(:)
      real (kind=Rkind), allocatable :: y(:)

!---------------------------------------------------------------------
      integer :: Get_symabOFSymAbelianOFBasis_AT_ib ! function

!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_make_polyorthobasis'
      logical,parameter :: debug=.FALSE.
!     logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL RecWrite_basis(base)
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------


      base%packed            = .TRUE.
      base%packed_done       = .TRUE.

       IF (base%ndim /= 1) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) ' ndim /= 1!',base%ndim
        STOP
       END IF

      IF (base%cplx) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' I cannot contract a COMPLEX basis set !'
        write(out_unitp,*) ' Not yet implemented'
        STOP
      END IF

      IF (Get_symabOFSymAbelianOFBasis_AT_ib(base,1) /= -1) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' I cannot use the symmetry'
        write(out_unitp,*) ' Not yet implemented'
        STOP
      END IF

!     - contraction -----------------------------------
      nq = get_nq_FROM_basis(base)
      CALL alloc_NParray(d0bc, (/ nq,base%nbc /),"d0bc",name_sub)
      CALL alloc_NParray(Wsuzy,(/ nq /),"Wsuzy",name_sub)
      CALL alloc_NParray(d1bc, (/ nq,base%nbc,base%ndim /),"d1bc",name_sub)
      CALL alloc_NParray(d2bc, (/ nq,base%nbc,base%ndim,base%ndim /),"d2bc",name_sub)

      IF (allocated(base%Rvec))  THEN
        CALL dealloc_NParray(base%Rvec,"base%Rvec",name_sub)
      END IF
      CALL alloc_NParray(base%Rvec,(/ base%nb,base%nb /),"base%Rvec",name_sub)
      base%Rvec(:,:) = ZERO
      nb0 = base%nb


      d0bc(:,:)     = ZERO
      d1bc(:,:,:)   = ZERO
      d2bc(:,:,:,:) = ZERO
      ! the first vector
      j=1
      !DO i=1,max(base%nb/10,base%nbc)
      DO i=1,base%nb
         d0bc(:,j)     = d0bc(:,j)     + Rvp(i,j) * base%dnRGB%d0(:,i)
         d1bc(:,j,:)   = d1bc(:,j,:)   + Rvp(i,j) * base%dnRGB%d1(:,i,:)
         d2bc(:,j,:,:) = d2bc(:,j,:,:) + Rvp(i,j) * base%dnRGB%d2(:,i,:,:)
      END DO
      phi0max = maxval(abs(d0bc(:,j)))
      Wsuzy(:) = -d1bc(:,j,1)/d0bc(:,j)
      DO i=1,nq
        IF (abs(d0bc(i,j))/phi0max > 0.0001_Rkind)                      &
        write(out_unitp,*) base%x(1,i),Wsuzy(i),d2bc(i,j,1,1)/d0bc(i,j)
      END DO
      Sjj = dot_product(d0bc(:,j),base%wrho(:)*d0bc(:,j))
      d0bc(:,j)     = d0bc(:,j)     / sqrt(Sjj)
      d1bc(:,j,1)   = d1bc(:,j,1)   / sqrt(Sjj)
      d2bc(:,j,1,1) = d2bc(:,j,1,1) / sqrt(Sjj)
      write(out_unitp,*) 'Sjj',j,Sjj


      ! recursion for the other vectors + orthonormalization
      x0 = dot_product(d0bc(:,1),base%x(1,:)*base%wrho(:)*d0bc(:,1))
      DO j=2,base%nbc
         d0bc(:,j)     =                           d0bc(:,j-1) * (base%x(1,:)-x0)
         d1bc(:,j,1)   = d0bc(:,j-1) +           d1bc(:,j-1,1) * (base%x(1,:)-x0)
         d2bc(:,j,1,1) = d1bc(:,j-1,1) * TWO + d2bc(:,j-1,1,1) * (base%x(1,:)-x0)

         DO i=1,j-1
           Sij = dot_product(d0bc(:,i),base%wrho(:)*d0bc(:,j))
           Sii = dot_product(d0bc(:,i),base%wrho(:)*d0bc(:,i))
           d0bc(:,j)     = d0bc(:,j)    *Sii - d0bc(:,i)    *Sij
           d1bc(:,j,1)   = d1bc(:,j,1)  *Sii - d1bc(:,i,1)  *Sij
           d2bc(:,j,1,1) = d2bc(:,j,1,1)*Sii - d2bc(:,i,1,1)*Sij

           Sij = dot_product(d0bc(:,i),base%wrho(:)*d0bc(:,j))
           d0bc(:,j)     = d0bc(:,j)    *Sii - d0bc(:,i)    *Sij
           d1bc(:,j,1)   = d1bc(:,j,1)  *Sii - d1bc(:,i,1)  *Sij
           d2bc(:,j,1,1) = d2bc(:,j,1,1)*Sii - d2bc(:,i,1,1)*Sij
         END DO
         Sjj = dot_product(d0bc(:,j),base%wrho(:)*d0bc(:,j))
         d0bc(:,j)     = d0bc(:,j)     / sqrt(Sjj)
         d1bc(:,j,1)   = d1bc(:,j,1)   / sqrt(Sjj)
         d2bc(:,j,1,1) = d2bc(:,j,1,1) / sqrt(Sjj)

         write(out_unitp,*) 'Sjj',j,Sjj
      END DO
!      -------------------------------------------------

      CALL alloc_NParray(d0b,(/ nq,base%nb /),"d0b",name_sub)
      d0b(:,:) = base%dnRGB%d0

!     - transfert to base% ----------------------------
      CALL dealloc_dnMat(base%dnRGB)
      CALL alloc_dnMat(base%dnRGB,nq,base%nbc,base%ndim,nderiv=2)

!     write(out_unitp,*) ' mv dnbc TO dnb'
      base%dnRGB%d0 = d0bc
      base%dnRGB%d1 = d1bc
      base%dnRGB%d2 = d2bc



      CALL Set_tab_SymAbelian(base%P_SymAbelian,base%nbc,symab=-1)
      base%nb  = base%nbc

      CALL ortho_basis_schmidt(base)
      CALL ortho_basis_schmidt(base)

      DO j=1,base%nbc
      DO i=1,nb0
        base%Rvec(i,j) = dot_product(d0b(:,i)*base%wrho(:),d0bc(:,j))
        !write(out_unitp,*) 'Rvec',j,i,base%Rvec(i,j)
      END DO
      END DO

      CALL dealloc_NParray(d0b,"d0b",name_sub)
      CALL dealloc_NParray(Wsuzy,"Wsuzy",name_sub)

      CALL dealloc_NParray(d0bc,"d0bc",name_sub)
      CALL dealloc_NParray(d1bc,"d1bc",name_sub)
      CALL dealloc_NParray(d2bc,"d2bc",name_sub)
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        CALL RecWrite_basis(base)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE sub_make_polyorthobasis

!=============================================================
!
!      check the overlap matrix of a basis
!
!=============================================================
      SUBROUTINE ortho_basis(basis_temp)
      USE mod_system
      USE mod_basis
      IMPLICIT NONE
!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis)  :: basis_temp

!---------- working variables ----------------------------------------
      real (kind=Rkind), allocatable :: tbasiswrho(:,:)
      real (kind=Rkind) :: matS(basis_temp%nb,basis_temp%nb)
      real (kind=Rkind) :: Rvp(basis_temp%nb,basis_temp%nb)
      real (kind=Rkind) :: Rdiag(basis_temp%nb)

      integer           ::  i,j,ierr,ind_i,nq
      real (kind=Rkind) ::  Sij,Sii,perm
      real (kind=Rkind) :: max_Sii,max_Sij

!---------------------------------------------------------------------
      integer :: err_mem,memory
      logical,parameter :: debug= .FALSE.
!     logical,parameter :: debug= .TRUE.
      character (len=*), parameter :: name_sub='ortho_basis'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL RecWrite_basis(basis_temp)
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------

      IF (basis_temp%cplx) THEN
        STOP 'ortho_basis in complex: not yet!'
      ELSE
        nq = get_nq_FROM_basis(basis_temp)
        CALL alloc_NParray(tbasiswrho,(/basis_temp%nb,nq/),"tbasiswrho",name_sub)

        DO i=1,basis_temp%nb
          tbasiswrho(i,:) = basis_temp%dnRGB%d0(:,i) * basis_temp%wrho(:)
        END DO
        matS(:,:) = matmul(tbasiswrho,basis_temp%dnRGB%d0)
        CALL dealloc_NParray(tbasiswrho,"tbasiswrho",name_sub)

        CALL Write_Mat(matS,out_unitp,5)

        DO i=1,basis_temp%nb
          Sii = matS(i,i) -ONE
          IF (abs(Sii) > max_Sii) max_Sii = abs(Sii)
          DO j=i+1,basis_temp%nb
            Sij = matS(i,j)
            IF (abs(Sij) > max_Sij) max_Sij = abs(Sij)
          END DO
        END DO

        write(out_unitp,11) max_Sii,max_Sij
 11     format('     Max Overlap:',2f15.12)
        IF (max_Sii > ONETENTH**5 .OR. max_Sij > ONETENTH**5) THEN
          write(out_unitp,*) ' WARNING in ',name_sub
          write(out_unitp,*) ' the basis is not orthonormal !!'
        END IF

        CALL diagonalization(matS,Rdiag,Rvp,basis_temp%nb,2,0,.FALSE.)

        write(out_unitp,*) 'Rvp:'
        CALL Write_Mat(Rvp,out_unitp,5)


        IF (max_Sii > ONETENTH**5 .OR. max_Sij > ONETENTH**5) THEN
          IF (allocated(basis_temp%Rvec))  THEN
            CALL dealloc_NParray(basis_temp%Rvec,'basis_temp%Rvec',name_sub)
          END IF
          CALL alloc_NParray(basis_temp%Rvec,                             &
                                      (/ basis_temp%nb,basis_temp%nb /),&
                          'basis_temp%Rvec',name_sub)
          basis_temp%nbc = basis_temp%nb
          basis_temp%Rvec = Rvp
          CALL sub_contraction_basis(basis_temp,.TRUE.)

          DO i=1,basis_temp%nb
            basis_temp%dnRGB%d0(:,i) = basis_temp%dnRGB%d0(:,i) / sqrt(Rdiag(i))
            basis_temp%dnRGB%d1(:,i,:) = basis_temp%dnRGB%d1(:,i,:) / sqrt(Rdiag(i))
            basis_temp%dnRGB%d2(:,i,:,:) = basis_temp%dnRGB%d2(:,i,:,:) / sqrt(Rdiag(i))
          END DO

          CALL dealloc_NParray(basis_temp%Rvec,'basis_temp%Rvec',name_sub)
        END IF

      END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(basis_temp)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


      end subroutine ortho_basis
      SUBROUTINE ortho_basis_schmidt(basis_temp)
      USE mod_system
      USE mod_basis
      IMPLICIT NONE
!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis)  :: basis_temp

!---------- working variables ----------------------------------------
      integer           ::  i,j
      real (kind=Rkind) ::  Sij,Sii,Sjj

!---------------------------------------------------------------------
      logical,parameter :: debug= .FALSE.
!     logical,parameter :: debug= .TRUE.
      character (len=*), parameter :: name_sub='ortho_basis_schmidt'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL RecWrite_basis(basis_temp)
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------

       CALL check_ortho_basis(basis_temp,.FALSE.)

      IF (basis_temp%cplx) THEN
        STOP 'ortho_basis in complex: not yet!'
      ELSE

        DO i=1,basis_temp%nb
          DO j=1,i-1
            Sij = dot_product(basis_temp%dnRGB%d0(:,i)*basis_temp%wrho(:),   &
                              basis_temp%dnRGB%d0(:,j))
            Sjj = dot_product(basis_temp%dnRGB%d0(:,j)*basis_temp%wrho(:),   &
                              basis_temp%dnRGB%d0(:,j))

            basis_temp%dnRGB%d0(:,i)     = Sjj * basis_temp%dnRGB%d0(:,i) -       &
                                      Sij * basis_temp%dnRGB%d0(:,j)
            basis_temp%dnRGB%d1(:,i,:)   = Sjj * basis_temp%dnRGB%d1(:,i,:) -     &
                                      Sij * basis_temp%dnRGB%d1(:,j,:)
            basis_temp%dnRGB%d2(:,i,:,:) = Sjj * basis_temp%dnRGB%d2(:,i,:,:) -   &
                                      Sij * basis_temp%dnRGB%d2(:,j,:,:)
          END DO
          Sii = dot_product(basis_temp%dnRGB%d0(:,i)*basis_temp%wrho(:),     &
                            basis_temp%dnRGB%d0(:,i))
          Sii = sqrt(Sii)
          basis_temp%dnRGB%d0(:,i)     = basis_temp%dnRGB%d0(:,i) / Sii
          basis_temp%dnRGB%d1(:,i,:)   = basis_temp%dnRGB%d1(:,i,:) / Sii
          basis_temp%dnRGB%d2(:,i,:,:) = basis_temp%dnRGB%d2(:,i,:,:) / Sii
        END DO

      END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(basis_temp)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


      end subroutine ortho_basis_schmidt
!=============================================================
!
!      Calculation of <d0b(:,i) | x | d0b(:,j) >
!
!=============================================================
      SUBROUTINE sub_matX_basis(basis_temp,matX_basis)
      USE mod_system
      USE mod_basis
      IMPLICIT NONE
!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis)  :: basis_temp
      real (kind=Rkind) :: matX_basis(basis_temp%nb,basis_temp%nb)

!---------- working variables ----------------------------------------
      integer       ::  i,nq
      real (kind=Rkind), allocatable :: tbasiswrhox(:,:)

!---------------------------------------------------------------------
      logical,parameter :: debug= .FALSE.
!     logical,parameter :: debug= .TRUE.
      character (len=*), parameter :: name_sub='sub_matX_basis'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL RecWrite_basis(basis_temp)
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
      IF (basis_temp%ndim > 1) THEN
        write(out_unitp,*) ' ERROR in',name_sub
        write(out_unitp,*) ' This transformation is not possible for basis with ndim>1'
        STOP
      END IF

      IF (basis_temp%cplx) THEN
        write(out_unitp,*) ' ERROR in',name_sub
        write(out_unitp,*) ' This transformation is not possible for complexe basis'
        STOP
      ELSE

        nq = get_nq_FROM_basis(basis_temp)
        CALL alloc_NParray(tbasiswrhox,(/basis_temp%nb,nq/),"tbasiswrhox",name_sub)

        DO i=1,basis_temp%nb
          tbasiswrhox(i,:) = basis_temp%dnRGB%d0(:,i) *                 &
                                    basis_temp%wrho(:)*basis_temp%x(1,:)
        END DO
        matX_basis(:,:) = matmul(tbasiswrhox,basis_temp%dnRGB%d0)
        CALL dealloc_NParray(tbasiswrhox,"tbasiswrho",name_sub)
        !CALL Write_Mat(matX_basis,out_unitp,5)
      END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'matX_basis:'
        CALL Write_Mat(matX_basis,out_unitp,5)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


      END SUBROUTINE sub_matX_basis
      SUBROUTINE sub_matperioX_basis(basis_temp,matX_basis)
      USE mod_system
      USE mod_basis
      IMPLICIT NONE
!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis)  :: basis_temp
      real (kind=Rkind) :: matX_basis(basis_temp%nb,basis_temp%nb)

!---------- working variables ----------------------------------------
      integer       ::  i,nq
      real (kind=Rkind), allocatable :: tbasiswrhox(:,:)

!---------------------------------------------------------------------
      logical,parameter :: debug= .FALSE.
!     logical,parameter :: debug= .TRUE.
      character (len=*), parameter :: name_sub='sub_matperioX_basis'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL RecWrite_basis(basis_temp)
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
      write(out_unitp,*) 'BEGINNING ',name_sub
      IF (basis_temp%ndim > 1) THEN
        write(out_unitp,*) ' ERROR in',name_sub
        write(out_unitp,*) ' This transformation is not possible for basis with ndim>1'
        STOP
      END IF

      IF (basis_temp%cplx) THEN
        write(out_unitp,*) ' ERROR in',name_sub
        write(out_unitp,*) ' This transformation is not possible for complexe basis'
        STOP
      ELSE
        nq = get_nq_FROM_basis(basis_temp)
        CALL alloc_NParray(tbasiswrhox,(/basis_temp%nb,nq/),"tbasiswrhox",name_sub)

        DO i=1,basis_temp%nb
          tbasiswrhox(i,:) = basis_temp%dnRGB%d0(:,i) *                 &
                      basis_temp%wrho(:)*cos(HALF*(basis_temp%x(1,:)+PI))
        END DO
        matX_basis(:,:) = matmul(tbasiswrhox,basis_temp%dnRGB%d0)

        CALL dealloc_NParray(tbasiswrhox,"tbasiswrho",name_sub)
        !CALL Write_Mat(matX_basis,out_unitp,5)
      END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'matX_basis:'
        CALL Write_Mat(matX_basis,out_unitp,5)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------
      write(out_unitp,*) 'END ',name_sub

      END SUBROUTINE sub_matperioX_basis
      SUBROUTINE sub_MatOFdnSX_basis(basis_temp)
      USE mod_system
      USE mod_dnSVM
      USE mod_basis
      IMPLICIT NONE
!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis)      :: basis_temp
      TYPE (Type_dnS)   :: MatOFdnSX(basis_temp%nb,basis_temp%nb)
      TYPE (Type_dnS)   :: EigVecdnS(basis_temp%nb,basis_temp%nb)


!---------- working variables ----------------------------------------
      integer           :: i,ib,i1,i2,nq
      real (kind=Rkind),allocatable :: tbasiswrhox(:,:)
      real (kind=Rkind) :: MatX(basis_temp%nb,basis_temp%nb)
      real (kind=Rkind) :: R0

!---------------------------------------------------------------------
!      logical,parameter :: debug= .FALSE.
      logical,parameter :: debug= .TRUE.
      character (len=*), parameter :: name_sub='sub_MatOFdnSX_basis'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL RecWrite_basis(basis_temp)
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
      IF (basis_temp%cplx) THEN
        write(out_unitp,*) ' ERROR in',name_sub
        write(out_unitp,*) ' This transformation is not possible for complexe basis'
        STOP
      END IF
      CALL alloc_MatOFdnS(MatOFdnSX,basis_temp%ndim,1)
      CALL sub_ZERO_TO_MatOFdnS(MatOFdnSX)
      R0 = ONE/real(basis_temp%ndim,kind=Rkind)

      nq = get_nq_FROM_basis(basis_temp)
      CALL alloc_NParray(tbasiswrhox,(/basis_temp%nb,nq/),"tbasiswrhox",name_sub)


      DO i=1,basis_temp%ndim
        DO ib=1,basis_temp%nb
          tbasiswrhox(ib,:) = basis_temp%dnRGB%d0(:,ib) *        &
                        basis_temp%wrho(:)*basis_temp%x(i,:)
        END DO
        MatX(:,:) = matmul(tbasiswrhox,basis_temp%dnRGB%d0)
        CALL dealloc_NParray(tbasiswrhox,"tbasiswrho",name_sub)
        !write(out_unitp,*) 'matX:',i
        !CALL Write_Mat(matX,out_unitp,5)


        DO i1=1,basis_temp%nb
        DO i2=1,basis_temp%nb
          MatOFdnSX(i1,i2)%d0    = MatOFdnSX(i1,i2)%d0 + R0 * MatX(i1,i2)
          MatOFdnSX(i1,i2)%d1(i) = MatX(i1,i2)
        END DO
        END DO

      END DO

      DO i=1,basis_temp%ndim
        DO i1=1,basis_temp%nb
        DO i2=1,basis_temp%nb
          MatOFdnSX(i1,i2)%d1(i) = MatOFdnSX(i1,i2)%d1(i) - MatOFdnSX(i1,i2)%d0
        END DO
        END DO
      END DO
      write(out_unitp,*) 'MatOFdnSX:'
      CALL Write_MatOFdnS(MatOFdnSX)


      CALL DIAG_MatOFdnS(MatOFdnSX,EigVecdnS,type_diago=3)


      DO i=1,basis_temp%ndim
        DO i1=1,basis_temp%nb
        DO i2=1,basis_temp%nb
          MatX(i1,i2) = MatOFdnSX(i1,i2)%d1(i)
        END DO
        END DO
      write(out_unitp,*) 'mat for i:',i
      CALL Write_Mat(MatX,6,5)

      END DO

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Diag MatOFdnSX:'
        CALL Write_MatOFdnS(MatOFdnSX)
        write(out_unitp,*) 'EigVecdnS:'
        CALL Write_MatOFdnS(EigVecdnS)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------
STOP

      END SUBROUTINE sub_MatOFdnSX_basis
!=============================================================
!
!      Transformation of basis:
!        sqrt(w(:)*rho(:))*d0b(:,i) TO d0b(:,i)
!        idem for d1b and d2b
!
!=============================================================
      SUBROUTINE sub_sqrtwrhoBasis_TO_basis(basis_temp)
      USE mod_system
      USE mod_basis
      IMPLICIT NONE
!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis)  :: basis_temp

!---------- working variables ----------------------------------------
      integer       ::  i,j,k

!---------------------------------------------------------------------
      logical,parameter :: debug= .FALSE.
!     logical,parameter :: debug= .TRUE.
      character (len=*), parameter :: name_sub='sub_sqrtwrhoBasis_TO_basis'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL RecWrite_basis(basis_temp)
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------

      basis_temp%wrho(:) = sqrt(basis_temp%wrho(:))

      IF (basis_temp%cplx) THEN
        DO i=1,basis_temp%nb
          basis_temp%dnCGB%d0(:,i)         = basis_temp%dnCGB%d0(:,i)    *&
                                    cmplx(basis_temp%wrho(:),kind=Rkind)
          DO j=1,basis_temp%ndim
            basis_temp%dnCGB%d1(:,i,j)     = basis_temp%dnCGB%d1(:,i,j)   *&
                                     cmplx(basis_temp%wrho(:),kind=Rkind)
            DO k=1,basis_temp%ndim
              basis_temp%dnCGB%d2(:,i,j,k) = basis_temp%dnCGB%d2(:,i,j,k) *&
                                       cmplx(basis_temp%wrho(:),kind=Rkind)
            END DO
          END DO
        END DO
      ELSE
        DO i=1,basis_temp%nb
          basis_temp%dnRGB%d0(:,i) = basis_temp%dnRGB%d0(:,i) * basis_temp%wrho(:)
          DO j=1,basis_temp%ndim
            basis_temp%dnRGB%d1(:,i,j) = basis_temp%dnRGB%d1(:,i,j) * basis_temp%wrho(:)
            DO k=1,basis_temp%ndim
              basis_temp%dnRGB%d2(:,i,j,k) = basis_temp%dnRGB%d2(:,i,j,k) * basis_temp%wrho(:)
            END DO
          END DO
        END DO
      END IF
      basis_temp%w(:) = ONE
      basis_temp%wrho(:) = ONE

!---------------------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(basis_temp)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------

      END SUBROUTINE sub_sqrtwrhoBasis_TO_basis

      SUBROUTINE Weight_OF_grid(w,d0RGB,nb,nq)
      USE mod_system
      USE mod_basis
      IMPLICIT NONE
!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      integer,           intent(in)    :: nb,nq
      real (kind=Rkind), intent(inout) :: w(nq)
      real (kind=Rkind), intent(in)    :: d0RGB(nq,nb)



!---------- working variables ----------------------------------------
      real (kind=Rkind) :: A(nb**2,nq)
      real (kind=Rkind) :: B(nb**2)

      real (kind=Rkind) :: AtA(nq,nq)
      real (kind=Rkind) :: AtB(nq)

      integer           ::  i,j,ij

!---------------------------------------------------------------------
      logical,parameter :: debug= .FALSE.
!     logical,parameter :: debug= .TRUE.
      character (len=*), parameter :: name_sub='Weight_OF_grid'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------

      ! matrices of the linear (rectangular) system A.W=B
      ij = 0
      B(:) = ZERO
      DO i=1,nb
      DO j=1,nb
        ij = ij + 1
        IF (i == j) B(ij) = ONE
        A(ij,:) = d0RGB(:,i)*d0RGB(:,j)
      END DO
      END DO

      ! matrices of the linear (square nq*nq) system AtA.W=AtB
      AtA(:,:) = matmul(transpose(A),A)
      AtB(:)   = matmul(transpose(A),B)

      ! Solve the linear system AtA.W=AtB
      CALL Linear_Sys(AtA,AtB,w,nq)

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'w',w
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


      END SUBROUTINE Weight_OF_grid

      SUBROUTINE Weight2_OF_grid_basis(basis_temp)
      USE mod_system
      USE mod_basis
      IMPLICIT NONE
!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis)  :: basis_temp

!---------- working variables ----------------------------------------
      real (kind=Rkind), allocatable :: A(:,:)
      real (kind=Rkind) :: B(basis_temp%nb**2)

      real (kind=Rkind),allocatable :: AtA(:,:)
      real (kind=Rkind),allocatable :: AtB(:)

      integer           ::  i,j,k,ij,nq
      logical           :: check_basis_save

!---------------------------------------------------------------------
      logical,parameter :: debug= .FALSE.
!     logical,parameter :: debug= .TRUE.
      character (len=*), parameter :: name_sub='Weight2_OF_grid_basis'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL RecWrite_basis(basis_temp)
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
      nq = get_nq_FROM_basis(basis_temp)
      ! write(out_unitp,*) 'nb,nq',basis_temp%nb,nq
      IF (.NOT. allocated(basis_temp%x) ) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' You cannot use this subroutine the grid is not done!'
        STOP
      END IF

      CALL Weight_OF_grid(basis_temp%w,basis_temp%dnRGB%d0,basis_temp%nb,nq)

!      CALL alloc_NParray(A,  (/basis_temp%nb**2,nq/),"A",  name_sub)
!      CALL alloc_NParray(AtA,(/nq,nq/),              "AtA",name_sub)
!      CALL alloc_NParray(AtB,(/nq/),                 "AtB",name_sub)
!
!
!
!      ! matrices of the linear (rectangular) system A.W=B
!      ij = 0
!      B(:) = ZERO
!      DO i=1,basis_temp%nb
!      DO j=1,basis_temp%nb
!        ij = ij + 1
!        IF (i == j) B(ij) = ONE
!        DO k=1,nq
!          A(ij,k) = basis_temp%dnRGB%d0(k,i)*basis_temp%dnRGB%d0(k,j)
!        END DO
!      END DO
!      END DO
!
!      ! matrices of the linear (square nq*nq) system AtA.W=AtB
!      AtA(:,:) = matmul(transpose(A),A)
!      AtB(:)   = matmul(transpose(A),B)
!
!
!      ! Solve the linear system AtA.W=AtB
!      CALL Linear_Sys(AtA,AtB,basis_temp%w(iqi:iqe),nq)
      basis_temp%wrho(:) = basis_temp%w(:) * basis_temp%rho(:)

!
!      CALL dealloc_NParray(A,  "A",  name_sub)
!      CALL dealloc_NParray(AtA,"AtA",name_sub)
!      CALL dealloc_NParray(AtB,"AtB",name_sub)

!---------------------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(basis_temp)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


      END SUBROUTINE Weight2_OF_grid_basis
      SUBROUTINE Weight_OF_grid_basis(basis_temp)
      USE mod_system
      USE mod_basis
      IMPLICIT NONE
!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis)  :: basis_temp

!---------- working variables ----------------------------------------
      real (kind=Rkind), allocatable :: A(:,:)
      real (kind=Rkind) :: B(basis_temp%nb*(basis_temp%nb+1)/2)

      real (kind=Rkind), allocatable :: AtA(:,:)
      real (kind=Rkind), allocatable :: AtB(:)

      integer           ::  i,j,k,ij,nq
      logical           :: check_basis_save

!---------------------------------------------------------------------
      logical,parameter :: debug= .FALSE.
      !logical,parameter :: debug= .TRUE.
      character (len=*), parameter :: name_sub='Weight_OF_grid_basis'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        !CALL RecWrite_basis(basis_temp)
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
      nq = get_nq_FROM_basis(basis_temp)
      ! write(out_unitp,*) 'nb,nq',basis_temp%nb,nq
      IF (.NOT. allocated(basis_temp%x) ) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' You cannot use this subroutine the grid is not done!'
        STOP
      END IF

      !write(out_unitp,*) 'nq',nq ; flush(out_unitp)
      CALL alloc_NParray(A, (/ basis_temp%nb*(basis_temp%nb+1)/2,nq /), &
                        "A",  name_sub)
      CALL alloc_NParray(AtA,(/nq,nq/),"AtA",name_sub)
      CALL alloc_NParray(AtB,(/nq/),   "AtB",name_sub)

      ! matrices of the linear (rectangular) system A.W=B
      ij = 0
      B(:) = ZERO
      DO i=1,basis_temp%nb
      DO j=i,basis_temp%nb
        ij = ij + 1
        IF (i == j) B(ij) = ONE
        DO k=1,nq
          A(ij,k) = basis_temp%dnRGB%d0(k,i)*basis_temp%dnRGB%d0(k,j)
        END DO
      END DO
      END DO

      ! matrices of the linear (square nq*nq) system AtA.W=AtB
      AtA(:,:) = matmul(transpose(A),A)
      AtB(:)   = matmul(transpose(A),B)


      ! Solve the linear system AtA.W=AtB
      CALL Linear_Sys(AtA,AtB,basis_temp%w,nq)
      IF (debug) write(out_unitp,*) 'w',basis_temp%w
      WHERE (basis_temp%w < ZERO)                                       &
                          basis_temp%w(:) = basis_temp%w(:) * 0.9_Rkind
      basis_temp%wrho(:) = basis_temp%w(:) * basis_temp%rho(:)



      CALL dealloc_NParray(A,  "A",  name_sub)
      CALL dealloc_NParray(AtA,"AtA",name_sub)
      CALL dealloc_NParray(AtB,"AtB",name_sub)

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'w',basis_temp%w
        !CALL RecWrite_basis(basis_temp)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


      END SUBROUTINE Weight_OF_grid_basis
      SUBROUTINE ConstantWeight_OF_grid_basis(basis_temp)
      USE mod_system
      USE mod_basis
      IMPLICIT NONE
!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis)  :: basis_temp

!---------- working variables ----------------------------------------
      real (kind=Rkind) :: SumIntExa,SumIntGrid

      integer           ::  i,j,k,ij,nq
      logical           :: check_basis_save

!---------------------------------------------------------------------
      logical,parameter :: debug= .FALSE.
!     logical,parameter :: debug= .TRUE.
      character (len=*), parameter :: name_sub='ConstantWeight_OF_grid_basis'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL RecWrite_basis(basis_temp)
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
      nq = get_nq_FROM_basis(basis_temp)
      ! write(out_unitp,*) 'nb,nq',basis_temp%nb,nq

      IF (.NOT. allocated(basis_temp%x) ) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' You cannot use this subroutine the grid is not done!'
        STOP
      END IF

      ! matrices of the linear (rectangular) system A.W=B
      SumIntExa  = real(basis_temp%nb,kind=Rkind)
      SumIntGrid = ZERO
      DO i=1,basis_temp%nb
      DO j=i,basis_temp%nb
        SumIntGrid = SumIntGrid +                                       &
         dot_product(basis_temp%dnRGB%d0(:,i),basis_temp%dnRGB%d0(:,j))
      END DO
      END DO

      basis_temp%w(:) = SumIntExa/SumIntGrid

      basis_temp%wrho(:) = basis_temp%w(:) * basis_temp%rho(:)

!---------------------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(basis_temp)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------

      END SUBROUTINE ConstantWeight_OF_grid_basis
      SUBROUTINE Make_grid_basis(basis_cuba)
      USE mod_system
      USE mod_basis
      IMPLICIT NONE
!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis), intent(inout) :: basis_cuba


!---------- working variables ----------------------------------------
      TYPE (basis)      :: basis_temp


      real (kind=Rkind) :: NormA


      !integer           :: i,j,k,ij,imc,iq,mc_max,Lq,idum
      integer           :: i,iq,nq,Lq,idum

      logical           :: check_basis_save,err_cuba

      integer                  :: nio,err_io
      TYPE (param_file)        :: cubature_file
      character (len=Name_len) :: name_i,name_j


      TYPE (param_SimulatedAnnealing) :: SA_para



!---------------------------------------------------------------------
      logical,parameter :: debug= .FALSE.
!     logical,parameter :: debug= .TRUE.
      character (len=*), parameter :: name_sub='Make_grid_basis'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL RecWrite_basis(basis_cuba)
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
      if ( basis_cuba%read_para_cubature ) then
        write(out_unitp,*) ' READ SA parameters'
        CALL Read_param_SimulatedAnnealing(SA_para)
        CALL Write_param_SimulatedAnnealing(SA_para)
        STOP
      end if

      nq = basis_cuba%nqc
      ! Save the old grid in basis_temp
      !CALL basis2TObasis1(basis_temp,basis_cuba)
      write(out_unitp,*) 'basis_cuba%nb',basis_cuba%nb
      write(out_unitp,*) 'basis_cuba%nq',nq

!      NormA = Norm_OF_grid_basis(basis_cuba,1,err_cuba)
!      write(out_unitp,*) ' Optimal norm',NormA
!      DO iq=1,nq
!        write(out_unitp,*) iq,basis_cuba%x(:,iq),basis_cuba%w(iq)
!      END DO
!
!      write(out_unitp,*) 'NormA',NormA

      Lq = int(basis_cuba%nDindB%MaxNorm)
      CALL Write_int_IN_char(basis_cuba%ndim,name_i)
      CALL Write_int_IN_char(Lq,             name_j)
      cubature_file%name = trim(name_i) // 'D_deg' // trim(name_j)
      write(out_unitp,*) 'cubature_file%name: ',cubature_file%name

!      IF (basis_cuba%Read_make_cubature) THEN
!        STOP 'Read_make_cubature: not yet'
!      ELSE
        SA_para%nb_mc_tot           = 1000000
        SA_para%nb_mc_partial       = 0             ! ?????

        SA_para%TempInit_type       =  1
        SA_para%Tmax                = -ONE
        SA_para%Tmin                =  ONETENTH**40
        SA_para%DeltaT              =  ZERO

        SA_para%RangeScal           =  0.8_Rkind
        SA_para%RangeScalInit       =  1._Rkind

        SA_para%With_RangeInit      = .FALSE.
        SA_para%RangeInit           = 1._Rkind

        SA_para%TempScheduling_type =  2 ! 1: linear, 2: geometrical (exp cooling) ...
        SA_para%ExpCoolParam        =  0.995_Rkind
        SA_para%ExpCoolParam        =  0.999_Rkind
        !SA_para%ExpCoolParam        =  0.99_Rkind

        SA_para%ResetTemp           = .TRUE.
        SA_para%ResetTempScal       =  ONE/THREE

        SA_para%Restart_Opt         =  3

!      END IF

      IF (basis_cuba%Restart_make_cubature) THEN

        CALL Set_nq_OF_basis(basis_cuba,basis_cuba%nqc)
        CALL alloc_dnb_OF_basis(basis_cuba)
        CALL alloc_xw_OF_basis(basis_cuba)
        basis_cuba%rho(:) = ONE

        nq = get_nq_FROM_basis(basis_cuba)

        CALL file_open(cubature_file,nio,err_file=err_io)
        read(nio,*,iostat=err_io) idum
        IF (idum /= nq) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' Wrong cubature file'
          write(out_unitp,*) ' nq from the file and the basis are different'
          write(out_unitp,*) ' nq:',idum,nq
          STOP
        END IF
        DO iq=1,nq
          read(nio,*) idum,basis_cuba%x(:,iq),basis_cuba%w(iq)
        END DO
        CALL file_close(cubature_file)

      ELSE
        CALL Make_grid_basis_SimulatedAnnealing(basis_cuba,1,SA_para,restart=.FALSE.)
      END IF
      SA_para%RangeScalInit = TEN**6
      DO i=1,SA_para%Restart_Opt
        CALL Make_grid_basis_SimulatedAnnealing(basis_cuba,1,SA_para,restart=.TRUE.)
      END DO


      basis_cuba%rho(:) = ONE
      basis_cuba%wrho(:) = basis_cuba%w(:)


      NormA = Norm_OF_grid_basis(basis_cuba,1,err_cuba)
      write(out_unitp,*) ' Optimal norm',NormA
      write(out_unitp,*) 'av',((sum(basis_cuba%x(i,:)) /                &
                               real(nq,kind=Rkind)),i=1,basis_cuba%ndim)
      CALL flush_perso(out_unitp)

      CALL file_open(cubature_file,nio,err_file=err_io)
      write(nio,*,iostat=err_io) nq
      DO iq=1,nq
        write(nio,*) iq,basis_cuba%x(:,iq),basis_cuba%w(iq)
      END DO
      write(nio,*) ' Rkind',Rkind
      write(nio,*) ' Optimal norm',NormA
      CALL file_close(cubature_file)


      STOP
      ! check the calulation
      check_basis_save       = basis_cuba%check_basis
      basis_cuba%check_basis = .TRUE.
      CALL check_ortho_basis(basis_cuba,test_stop=.TRUE.)
      basis_cuba%check_basis = check_basis_save

      !write(out_unitp,*) 'new w',basis_temp%w

!---------------------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(basis_cuba)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


      END SUBROUTINE Make_grid_basis
      SUBROUTINE Make_grid_basis_SimulatedAnnealing(basis_cuba,type_weight,SA_para,restart)
      USE mod_system
      USE mod_basis
      IMPLICIT NONE
!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis),                    intent(inout) :: basis_cuba
      integer,                         intent(in)    :: type_weight
      TYPE (param_SimulatedAnnealing), intent(in)    :: SA_para
      logical,                         intent(in)    :: restart


!---------- working variables ----------------------------------------
      real (kind=Rkind) :: x0(basis_cuba%ndim,basis_cuba%nqc)
      real (kind=Rkind) :: xmin(basis_cuba%ndim,basis_cuba%nqc),Norm_min,Norm_max

      real (kind=Rkind) :: DNorm,NormA,NormB,Temp,Temp_max,PTemp,DTemp
      real (kind=Rkind) :: x(basis_cuba%ndim),xi,xav,S

      real (kind=Rkind) :: QA(basis_cuba%ndim)
      real (kind=Rkind) :: QB(basis_cuba%ndim)
      real (kind=Rkind) :: SQ(basis_cuba%ndim)
      real (kind=Rkind) :: QAv(basis_cuba%ndim)

      integer           :: i,j,k,ij,imc,iq,nb_Norm_min,nb_block_WithoutMin,option_rand_grid

      logical           :: check_basis_save,err_cuba,Norm_down


!---------------------------------------------------------------------
      logical,parameter :: debug= .FALSE.
!     logical,parameter :: debug= .TRUE.
      character (len=*), parameter :: name_sub='Make_grid_basis_SimulatedAnnealing'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL RecWrite_basis(basis_cuba)
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
      option_rand_grid = 0
      Norm_down        = .FALSE.


      DO i=1,basis_cuba%ndim
        QA(i) = minval(basis_cuba%x(i,:))
        QB(i) = maxval(basis_cuba%x(i,:))
        QAv(i)= (QA(i)+QB(i))*HALF
        SQ(i) = (QB(i)-QA(i))/TWO
      END DO

      IF (restart) THEN

        Norm_min  = Norm_OF_grid_basis(basis_cuba,type_weight,err_cuba) ! initial norm
        write(out_unitp,*) 'Initial Norm',Norm_min

        SQ(:)     = SQ(:)*min(ONE,Norm_min*SA_para%RangeScalInit)

        x0(:,:)   = basis_cuba%x(:,:)


      ELSE
        SQ(:) = SQ(:)*SA_para%RangeScalInit

        ! because, nq has changed (nq->nqc)
        CALL Set_nq_OF_basis(basis_cuba,basis_cuba%nqc)
        CALL alloc_dnb_OF_basis(basis_cuba)
        CALL alloc_xw_OF_basis(basis_cuba)

        DO i=1,basis_cuba%ndim
          x0(i,:) = QAv(i)
        END DO
        CALL ReOriented_grid(x0,basis_cuba%ndim,basis_cuba%nqc)
        CALL ReCentered_grid(x0,basis_cuba%ndim,basis_cuba%nqc)

      END IF
      write(out_unitp,*) 'QA',QA
      write(out_unitp,*) 'QB',QB
      write(out_unitp,*) 'SQ',SQ

      Norm_min  = Norm_OF_grid_basis(basis_cuba,type_weight,err_cuba)
      xmin(:,:) = x0(:,:)
      Norm_max  = ZERO


      write(out_unitp,'(a,i0,a)') 'Initial Temperature, with ',SA_para%nb_mc_tot/10,' evaluations'
      flush(out_unitp)

      write(out_unitp,'(a)') 'Eval (%): [--0-10-20-30-40-50-60-70-80-90-100]'
      write(out_unitp,'(a)',ADVANCE='no') 'Eval (%): ['
      CALL flush_perso(out_unitp)
      NormA  = ZERO
      ! first find the average Energy (Norm), then Temp
      DO imc=1,SA_para%nb_mc_tot/10
        err_cuba = .TRUE.
        DO WHILE (err_cuba)
          CALL Random_grid(basis_cuba%x,x0,SQ,QA,QB,basis_cuba%ndim,basis_cuba%nqc,option_rand_grid)
          NormB = Norm_OF_grid_basis(basis_cuba,type_weight,err_cuba)
        END DO
        !write(out_unitp,*) 'NormB',imc,NormB
        NormA = NormA + NormB
        IF (NormB > Norm_max) Norm_max = NormB

        IF (NormB < Norm_min) THEN
          xmin(:,:) = basis_cuba%x(:,:)
          Norm_down = .TRUE.
          Norm_min  = NormB
        END IF

        IF (mod(imc,max(1,int(SA_para%nb_mc_tot/100))) == 0 .AND. MPI_id==0) THEN
          write(out_unitp,'(a)',ADVANCE='no') '---'
          flush(out_unitp)
        END IF

        IF (NormB < SA_para%Tmin) EXIT
      END DO
      IF(MPI_id==0) write(out_unitp,'(a)',ADVANCE='yes') '----]'

      NormA = NormA / real(SA_para%nb_mc_tot/10,kind=Rkind)
      write(out_unitp,*) 'Min, Average, Max Norm',Norm_min,NormA,Norm_max
      flush(out_unitp)

      IF (NormB < SA_para%Tmin) RETURN

      basis_cuba%x(:,:) = xmin(:,:)
      Norm_min          = Norm_OF_grid_basis(basis_cuba,type_weight,err_cuba)
      x0(:,:)           = basis_cuba%x(:,:)
      NormB             = Norm_min
      write(out_unitp,*) ' Norm (cubature)',xmin,Norm_min

      Temp_max = Norm_max-Norm_min
      Temp_max = NormA

      Temp     = Temp_max
      write(out_unitp,*) 'Average Norm, Temp',NormA,Temp
      CALL flush_perso(out_unitp)

      DTemp               = Temp_max/real(SA_para%nb_mc_tot,kind=Rkind)
      imc                 = 1
      nb_Norm_min         = 0
      nb_block_WithoutMin = 0

      DO

        err_cuba = .TRUE.
        DO WHILE (err_cuba)
          CALL Random_grid(basis_cuba%x,x0,SQ,QA,QB,basis_cuba%ndim,basis_cuba%nqc,option_rand_grid)
          NormA = Norm_OF_grid_basis(basis_cuba,type_weight,err_cuba)
        END DO
        DNorm = NormA - NormB
        !write(out_unitp,*) 'Norm',imc,NormA
        CALL flush_perso(out_unitp)

        IF ( NormA < Norm_min) THEN
          nb_Norm_min = nb_Norm_min + 1
          xmin(:,:)   = basis_cuba%x(:,:)
          Norm_down   = .TRUE.
          Norm_min    = NormA
          !write(out_unitp,*) ' imc, Temp, Norm_min',imc,Temp,Norm_min
        END IF


        IF ( DNorm < ZERO) THEN
          x0(:,:)     = basis_cuba%x(:,:)
          NormB       = NormA
          !write(out_unitp,*) ' imc, Temp, Norm (cubature)',imc,Temp,NormA
          !write(out_unitp,*) ' accepted, DNorm',DNorm
        ELSE
          CALL random_number(PTemp)
          !Ptemp = genrand64_real1()
          IF ( PTemp < exp(-DNorm/Temp) ) THEN
            x0(:,:)   = basis_cuba%x(:,:)
            NormB     = NormA
            !write(out_unitp,*) ' imc, Temp, Norm (propa)',imc,Temp,NormA
            !write(out_unitp,*) ' accepted, Proba',PTemp,exp(-DNorm/Temp)
          END IF
        END IF
        imc = imc + 1

        SELECT CASE(SA_para%TempScheduling_type)
        CASE (1)
          Temp = Temp - DTemp        ! linear cooling
        CASE (2)
          Temp = SA_para%ExpCoolParam * Temp ! Exponential cooling
        CASE DEFAULT
          Temp = SA_para%ExpCoolParam * Temp ! Exponential cooling
        END SELECT
        !write(out_unitp,*) 'imc,Temp,Temp_max,ResetTempScal,Norm_min',imc,Temp,Temp_max,ResetTempScal

        IF (Temp < SA_para%ResetTempScal*Temp_max) THEN

          write(out_unitp,*) 'imc,Temp_max,nb_Norm_min,Norm_min',imc,   &
                                           Temp_max,nb_Norm_min,Norm_min
          CALL flush_perso(out_unitp)

          SQ(:)               = SQ(:)*SA_para%RangeScal
          Temp_max            = (ONE-SA_para%ResetTempScal)*Temp_max
          Temp                = Temp_max
          DTemp               = Temp_max/real(max(10,SA_para%nb_mc_tot-imc),kind=Rkind)

          x0(:,:)             = xmin(:,:)
          IF (nb_Norm_min == 0) THEN
            nb_block_WithoutMin = nb_block_WithoutMin + 1
          ELSE
            nb_block_WithoutMin = 0
            nb_Norm_min         = 0
          END IF
        END IF

        IF (nb_block_WithoutMin > 20) THEN
          write(out_unitp,*) 'restart because, nb_block_WithoutMin > 20'
          flush(out_unitp)
          EXIT
        END IF
        IF (Temp < SA_para%Tmin .OR. imc > SA_para%nb_mc_tot) EXIT

      END DO

      ! optimal values.
      basis_cuba%x(:,:) = xmin(:,:)


!---------------------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(basis_cuba)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


      END SUBROUTINE Make_grid_basis_SimulatedAnnealing

      FUNCTION Norm_OF_grid_basis(basis_cuba,type_weight,err_cuba)
      USE mod_system
      USE mod_basis
      IMPLICIT NONE
!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      real (kind=Rkind)           :: Norm_OF_grid_basis
      TYPE (basis), intent(inout) :: basis_cuba
      logical, intent(inout)      :: err_cuba
      integer, intent(in)         :: type_weight



!---------- working variables ----------------------------------------
      real (kind=Rkind), allocatable :: basiswrho(:,:)
      real (kind=Rkind) :: S(basis_cuba%nb,basis_cuba%nb)

      real (kind=Rkind) :: Norm,Norm_max,Norm_A,Norm_ini

      integer           :: i,j,k,ij
      integer           :: nq

!---------------------------------------------------------------------
      logical,parameter :: debug= .FALSE.
!     logical,parameter :: debug= .TRUE.
      character (len=*), parameter :: name_sub='Norm_OF_grid_basis'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL RecWrite_basis(basis_cuba)
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
      nq = get_nq_FROM_basis(basis_cuba)

      DO i=1,basis_cuba%nb
      DO j=1,nq
        basis_cuba%dnRGB%d0(j,i) = Rec_d0bnD_AT_Q(basis_cuba,i,basis_cuba%x(:,j))
      END DO
      END DO

      ! then we calculate the weight
      !write(out_unitp,*) 'old',basis_cuba%w
      IF (type_weight == 0) THEN
        CALL ConstantWeight_OF_grid_basis(basis_cuba)
      ELSE
        CALL Weight_OF_grid_basis(basis_cuba)
        !CALL Weight2_OF_grid_basis(basis_cuba)
      END IF
      err_cuba = .FALSE.

      CALL alloc_NParray(basiswrho,(/ nq,basis_cuba%nb /),                &
                        "basiswrho",name_sub)

      ! basiswrho
      DO i=1,basis_cuba%nb
        basiswrho(:,i) = basis_cuba%dnRGB%d0(:,i)*basis_cuba%w(:)
      END DO
      S(:,:) = matmul(transpose(basiswrho),basis_cuba%dnRGB%d0)
      DO i=1,basis_cuba%nb
        S(i,i) = S(i,i)-ONE
      END DO
      Norm_OF_grid_basis = sqrt(sum(S*S)/real(basis_cuba%nb**2,kind=Rkind))

      !Norm_OF_grid_basis = sum(S*S)*HALF
      !Norm_OF_grid_basis = sqrt(sum(S*S)*HALF)

      RETURN


      Norm     = ZERO
      Norm_max = ZERO
      DO i=1,basis_cuba%nb
        Norm_A = (dot_product(basis_cuba%dnRGB%d0(:,i),basiswrho(:,i))-ONE)**2
        IF (Norm_A > Norm_max) Norm_max = Norm_A
        Norm = Norm + Norm_A
        DO j=i+1,basis_cuba%nb
        !DO j=1,basis_cuba%nb
          Norm_A = dot_product(basis_cuba%dnRGB%d0(:,i),basiswrho(:,j))**2
          IF (Norm_A > Norm_max) Norm_max = Norm_A
          Norm = Norm + Norm_A
        END DO
      END DO
      Norm     = Norm + Norm_ini
      Norm_max = Norm_max + Norm_ini

      !Norm_OF_grid_basis = sqrt(Norm)
      !Norm_OF_grid_basis = sqrt(Norm)/real(basis_cuba%nb*(basis_cuba%nb+1)/2,kind=Rkind)
      Norm_OF_grid_basis = Norm
      !Norm_OF_grid_basis = Norm_max


!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) ' Norm (cubature)',Norm
        CALL RecWrite_basis(basis_cuba)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------

      END FUNCTION Norm_OF_grid_basis

      SUBROUTINE ReCentered_grid(x,ndim,nqc)
      USE mod_system
      USE mod_basis
      IMPLICIT NONE


      integer :: ndim,nqc
      real (kind=Rkind) :: x(ndim,nqc)


      integer           :: i

!---------------------------------------------------------------------
      logical,parameter :: debug= .FALSE.
!     logical,parameter :: debug= .TRUE.
      character (len=*), parameter :: name_sub='ReCentered_grid'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        DO i=1,nqc
          write(out_unitp,*) 'i,x(:,i)',i,x(:,i)
        END DO
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------

      DO i=1,ndim
        x(i,:) = x(i,:) - sum(x(i,:))/real(nqc,kind=Rkind)
      END DO

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Recentered grid'
        DO i=1,nqc
          write(out_unitp,*) 'i,x(:,i)',i,x(:,i)
        END DO
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


      END SUBROUTINE ReCentered_grid
      SUBROUTINE ReOriented_grid(x,ndim,nqc)
      USE mod_system
      USE mod_basis
      IMPLICIT NONE


      integer :: ndim,nqc
      real (kind=Rkind) :: x(ndim,nqc)


      integer           :: i
      real (kind=Rkind) :: NR,C,S,MatRot(ndim,ndim)

      integer :: option = 1

!---------------------------------------------------------------------
      logical,parameter :: debug= .FALSE.
!      logical,parameter :: debug= .TRUE.
      character (len=*), parameter :: name_sub='ReOriented_grid'
!---------------------------------------------------------------------
      IF (ndim == 1) RETURN

      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        DO i=1,nqc
          write(out_unitp,*) 'i,x(:,i)',i,x(:,i)
        END DO
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------

      IF (option == 0) THEN

        DO i=1,ndim-1
          IF (i > nqc) CYCLE
          x(1:ndim-i,i) = ZERO ! to avoid rotation
        END DO

      ELSE
        MatRot(:,:) = ZERO
        DO i=1,ndim
          MatRot(i,i) = ONE
        END DO
        NR = sqrt(dot_product(x(:,1),x(:,1)))
        IF (NR < ONETENTH**5) RETURN
        C  = x(1,1) / NR
        S  = sqrt(ONE-C*C)
        MatRot(1,1:2) = (/  C, S /)
        MatRot(2,1:2) = (/ -S, C /)

        DO i=1,nqc
          x(:,i) = matmul(MatRot,x(:,i))
        END DO
      END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'ReOriented grid'
        DO i=1,nqc
          write(out_unitp,*) 'i,x(:,i)',i,x(:,i)
        END DO
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


      END SUBROUTINE ReOriented_grid
      SUBROUTINE Random_grid(x,x0,SQ,QA,QB,ndim,nqc,option)
      USE mod_system
      USE mod_basis
      IMPLICIT NONE


      integer :: ndim,nqc,option
      real (kind=Rkind) :: x(ndim,nqc)

      real (kind=Rkind) :: x0(ndim,nqc),SQ(ndim),QA(ndim),QB(ndim)



      integer           :: iq,i,k

      real (kind=Rkind) :: xi,S


!---------------------------------------------------------------------
      logical,parameter :: debug= .FALSE.
!     logical,parameter :: debug= .TRUE.
      character (len=*), parameter :: name_sub='Random_grid'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        DO i=1,nqc
          write(out_unitp,*) 'i,x(:,i)',i,x(:,i)
        END DO
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
  IF (option == 0) THEN

    DO iq=1,nqc
    DO i=1,ndim
      DO k=1,1000
        CALL random_number(xi)
        xi = TWO*xi-ONE
        !S=minval( (/ QB(i)-x0(i,iq),x0(i,iq)-QA(i),SQ(i) /) )
        S=SQ(i)
        x(i,iq) = x0(i,iq) + xi*S
        IF (x(i,iq) <= QB(i) .AND. x(i,iq) >= QA(i)) EXIT
      END DO
    END DO
    END DO

  ELSE

    DO iq=1,nqc
    DO i=1,ndim
      CALL random_number(xi)
      !xi = genrand64_real1()
      xi = TWO*xi-ONE
      S=SQ(i)
      x(i,iq) = x0(i,iq) + xi*S
    END DO
    END DO

  END IF

  CALL ReOriented_grid(x,ndim,nqc)
  CALL ReCentered_grid(x,ndim,nqc)


!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Random grid'
        DO i=1,nqc
          write(out_unitp,*) 'i,x(:,i)',i,x(:,i)
        END DO
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


      END SUBROUTINE Random_grid
END MODULE BasisMakeGrid
