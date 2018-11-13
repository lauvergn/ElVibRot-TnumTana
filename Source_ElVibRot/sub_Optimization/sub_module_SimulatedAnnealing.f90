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

      MODULE mod_SimulatedAnnealing
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

      SUBROUTINE Sub_SimulatedAnnealing(BasisnD,xOpt_min,Norm_min,      &
                                        SQ,nb_Opt,                      &
                                     para_Tnum,mole,ComOp,para_PES,Qact,&
                                        para_SimulatedAnnealing)

      USE mod_system
      USE mod_Coord_KEO
      USE mod_PrimOp
      USE mod_basis
      USE mod_Op
      USE mod_Auto_Basis
      IMPLICIT NONE

!----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole
      TYPE (Tnum)    :: para_Tnum
      logical        :: Cart_Transfo_save
      real (kind=Rkind), intent(inout) :: Qact(:)

!----- for the basis set ----------------------------------------------
      TYPE (basis)          :: BasisnD

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_PES) :: para_PES
      integer          :: nb_scalar_Op
      logical          :: calc_scalar_Op

!----- variables for the construction of H ---------------------------
      TYPE (param_ComOp)          :: ComOp
      TYPE (param_ReadOp)         :: para_ReadOp
      logical                     :: Save_FileGrid,Save_MemGrid


!----- local variables -----------------------------------------------
!----- variables for the construction of H ---------------------------
      TYPE (param_AllOp)          :: para_AllOp_loc

!----- for the basis set ----------------------------------------------
      TYPE (param_AllBasis)       :: para_AllBasis_loc
      TYPE (basis)                :: basis_temp

!----- for the optimization -------------------------------------------
      TYPE (param_SimulatedAnnealing) :: para_SimulatedAnnealing
      integer, intent(in) :: nb_Opt
      real (kind=Rkind), intent(inout) :: xOpt_min(nb_Opt),SQ(nb_Opt)

!---------- working variables ----------------------------------------
      real (kind=Rkind) :: Norm_min,Norm_max

      real (kind=Rkind) :: DNorm,NormA,NormB
      real (kind=Rkind) :: Temp,Temp_max,PTemp,DTemp


      real (kind=Rkind) :: xi,xav
      real (kind=Rkind), pointer :: x(:)
      real (kind=Rkind), pointer :: x0(:)
      real (kind=Rkind), pointer :: QA(:)
      real (kind=Rkind), pointer :: QB(:)
      real (kind=Rkind), pointer :: QAv(:)

      real (kind=Rkind), pointer :: xOpt(:)

      integer           :: i,imc,nb_Norm_min
      integer           :: nb_Accepted_DNorm,nb_Accepted_proba,nb_NotAccepted_proba


!---------------------------------------------------------------------
      logical,parameter :: debug= .FALSE.
!     logical,parameter :: debug= .TRUE.
      character (len=*), parameter :: name_sub='Sub_SimulatedAnnealing'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------

      ! allocation x ...
      nullify(x)
      nullify(x0)
      nullify(QA)
      nullify(QB)
      nullify(QAv)
      nullify(xOpt)

      CALL alloc_array(x,   (/ nb_Opt /),'x',   name_sub)
      CALL alloc_array(x0,  (/ nb_Opt /),'x0',  name_sub)
      CALL alloc_array(xOpt,(/ nb_Opt /),'xOpt',name_sub)
      CALL alloc_array(QA,  (/ nb_Opt /),'QA',  name_sub)
      CALL alloc_array(QB,  (/ nb_Opt /),'QB',  name_sub)
      CALL alloc_array(QAv, (/ nb_Opt /),'Qav', name_sub)
      ! end allocation

      CALL Write_param_SimulatedAnnealing(para_SimulatedAnnealing)

      print_level = 0
      CALL Sub_Energ_OF_ParamBasis(Norm_min,xOpt_min,nb_Opt,BasisnD,    &
                                   para_Tnum,mole,ComOp,para_PES,Qact)
      write(out_unitp,*) 'Initial param',xOpt_min,' Energy',Norm_min
          !STOP 'coucou'


      QA(:) = xOpt_min(:) - SQ(:)*HALF
      QB(:) = xOpt_min(:) + SQ(:)*HALF
      QAv(:)= xOpt_min(:)
      write(out_unitp,*) 'QA,QB',QA,QB


      Norm_max = Norm_min
      NormA    = ZERO
      print_level = -1
      ! first find the average Energy (Norm), then Temp
      DO imc=1,para_SimulatedAnnealing%nb_mc_tot/10
          CALL random_number(x)
          xOpt(:) = QA(:) + x(:)*SQ(:)

          CALL Sub_Energ_OF_ParamBasis(NormB,xOpt,nb_Opt,BasisnD,       &
                                       para_Tnum,mole,ComOp,para_PES,Qact)

          !write(6,*) 'NormB',imc,NormB
          NormA = NormA + NormB
          IF (NormB > Norm_max) Norm_max = NormB

          IF (NormB < Norm_min) THEN
            xOpt_min(:) = xOpt(:)
            Norm_min = NormB
            write(out_unitp,*) 'NormB',imc,NormB
          END IF
      END DO
      NormA = NormA / real(para_SimulatedAnnealing%nb_mc_tot/10,kind=Rkind)
      write(out_unitp,*) 'Min, Average, Max Energy',Norm_min,NormA,Norm_max

      x0(:)  = xOpt_min(:)
      write(out_unitp,*) ' Energy_min',xOpt_min,Norm_min

      Temp_max = Norm_max-Norm_min
      Temp     = Temp_max
      write(out_unitp,*) 'Average Energy, Temp',NormA,Temp
      CALL flush_perso(out_unitp)

      DTemp       = Temp_max/real(max(10,para_SimulatedAnnealing%nb_mc_tot),kind=Rkind)

      imc = 1
      nb_Norm_min = 0
      DO

        DO i=1,nb_Opt
          DO
            CALL random_number(xi)
            xi = TWO*xi-ONE
            xOpt(i) = x0(i) + xi*SQ(i)
            IF (xOpt(i) <= QB(i) .AND. xOpt(i) >= QA(i)) EXIT
          END DO
        END DO

        CALL Sub_Energ_OF_ParamBasis(NormA,xOpt,nb_Opt,BasisnD,         &
                                     para_Tnum,mole,ComOp,para_PES,Qact)

        DNorm = NormA - NormB
        !write(out_unitp,*) 'Norm',imc,xOpt,NormA
        CALL flush_perso(out_unitp)

        IF ( NormA < Norm_min) THEN
          nb_Norm_min   = nb_Norm_min + 1
          xOpt_min(:)   = xOpt(:)
          Norm_min      = NormA
          write(out_unitp,*) ' imc, Temp, nb_Energy_min, Energy_min',imc,Temp,nb_Norm_min,Norm_min
        END IF


        IF (DNorm < ZERO) THEN
          x0(:)       = xOpt(:)
          NormB       = NormA
        ELSE
          CALL random_number(PTemp)
          IF ( PTemp < exp(-DNorm/Temp) ) THEN
            x0(:)       = xOpt(:)
            NormB       = NormA
            !write(out_unitp,*) ' imc, Temp, Norm (propa)',imc,Temp,NormA
            !write(out_unitp,*) ' accepted, Proba',PTemp,exp(-DNorm/Temp)
          END IF
        END IF
        imc = imc + 1
        SELECT CASE (para_SimulatedAnnealing%TempScheduling_type)
        CASE (1)
          !write(6,*) 'LinCoolParam'
          Temp = Temp - DTemp
        CASE (2)
          !write(6,*) 'ExpCoolParam'
          Temp = para_SimulatedAnnealing%ExpCoolParam * Temp
        CASE DEFAULT
          !write(6,*) 'ExpCoolParam'
          Temp = para_SimulatedAnnealing%ExpCoolParam * Temp
        END SELECT

        IF (para_SimulatedAnnealing%ResetTemp .AND.                     &
             Temp < para_SimulatedAnnealing%ResetTempScal*Temp_max) THEN

          write(out_unitp,*) 'Temp_max,nb_Energy_min,Energy_min',           &
                                           Temp_max,nb_Norm_min,Norm_min
          CALL flush_perso(out_unitp)
          SQ(:)       = SQ(:)*para_SimulatedAnnealing%RangeScal
          IF (nb_Norm_min == 0) THEN ! another scaling+no reset Tmax and DTemp
            SQ(:) = SQ(:)*para_SimulatedAnnealing%RangeScal ! another scaling
          ELSE
            Temp_max    = (ONE-para_SimulatedAnnealing%ResetTempScal)*Temp_max
            DTemp       = Temp_max/real(max(10,para_SimulatedAnnealing%nb_mc_tot-imc),kind=Rkind)
          END IF


          Temp        = Temp_max
          x0(:)       = xOpt_min(:)
          nb_Norm_min = 0
        END IF

        IF (Temp < para_SimulatedAnnealing%Tmin .OR. imc > para_SimulatedAnnealing%nb_mc_tot) EXIT
      END DO
      write(out_unitp,*) 'Optimal param',xOpt_min,' Energy',Norm_min

      !deallocation
      CALL dealloc_array(x,   'x',   name_sub)
      CALL dealloc_array(x0,  'x0',  name_sub)
      CALL dealloc_array(xOpt,'xOpt',name_sub)
      CALL dealloc_array(QA,  'QA',  name_sub)
      CALL dealloc_array(QB,  'QB',  name_sub)
      CALL dealloc_array(QAv, 'Qav', name_sub)
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


      END SUBROUTINE Sub_SimulatedAnnealing

      SUBROUTINE Sub_SimulatedAnnealing_cuba(BasisnD,xOpt_min,Norm_min, &
                                             SQ,QA,QB,nb_Opt,           &
                                    para_Tnum,mole,ComOp,para_PES,Qact, &
                                             para_SimulatedAnnealing)

      USE mod_system
      USE mod_Coord_KEO
      USE mod_PrimOp
      USE mod_basis
      USE mod_Op
      USE mod_Auto_Basis
      IMPLICIT NONE

!----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole
      TYPE (Tnum)    :: para_Tnum
      logical        :: Cart_Transfo_save
      real (kind=Rkind), intent(inout) :: Qact(:)

!----- for the basis set ----------------------------------------------
      TYPE (basis)          :: BasisnD

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_PES) :: para_PES
      integer          :: nb_scalar_Op
      logical          :: calc_scalar_Op

!----- variables for the construction of H ---------------------------
      TYPE (param_ComOp)          :: ComOp
      TYPE (param_ReadOp)         :: para_ReadOp
      logical                     :: Save_FileGrid,Save_MemGrid

!----- local variables -----------------------------------------------
!----- variables for the construction of H ---------------------------
      TYPE (param_AllOp)          :: para_AllOp_loc

!----- for the basis set ----------------------------------------------
      TYPE (param_AllBasis)       :: para_AllBasis_loc
      TYPE (basis)                :: basis_temp

!----- for the optimization -------------------------------------------
      TYPE (param_SimulatedAnnealing)  :: para_SimulatedAnnealing
      integer, intent(in)              :: nb_Opt
      real (kind=Rkind), intent(inout) :: xOpt_min(nb_Opt),SQ(nb_Opt),  &
                                          QA(nb_Opt),QB(nb_Opt)
!---------- working variables ----------------------------------------
      real (kind=Rkind) :: Norm_min,Norm_max

      real (kind=Rkind) :: DNorm,NormA,NormB
      real (kind=Rkind) :: Temp,Temp_max,PTemp,DTemp


      real (kind=Rkind) :: xi,xav
      real (kind=Rkind), pointer :: x(:)
      real (kind=Rkind), pointer :: x0(:)
      real (kind=Rkind), pointer :: QAv(:)

      real (kind=Rkind), pointer :: xOpt(:)

      integer           :: i,k,imc,imc_block,nb_Norm_min,nb_block_WithoutMin
      integer           :: nb_Accepted_DNorm,nb_Accepted_proba,nb_NotAccepted_proba

!---------------------------------------------------------------------
      logical,parameter :: debug= .FALSE.
!     logical,parameter :: debug= .TRUE.
      character (len=*), parameter :: name_sub='Sub_SimulatedAnnealing_cuba'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------

      ! allocation x ...
      nullify(x)
      nullify(x0)
      nullify(QAv)
      nullify(xOpt)

      CALL alloc_array(x,   (/ nb_Opt /),'x',   name_sub)
      CALL alloc_array(x0,  (/ nb_Opt /),'x0',  name_sub)
      CALL alloc_array(xOpt,(/ nb_Opt /),'xOpt',name_sub)
      CALL alloc_array(QAv, (/ nb_Opt /),'Qav', name_sub)
      ! end allocation

      CALL Write_param_SimulatedAnnealing(para_SimulatedAnnealing)


      print_level = 0
      QAv(:)= xOpt_min(:)
      x0(:) = xOpt_min(:)
      DO i=1,nb_Opt
        DO k=1,1000
          CALL random_number(xi)
          xi = TWO*xi-ONE
          xOpt(i) = x0(i) + xi*SQ(i)/TWO
          IF (xOpt(i) <= QB(i) .AND. xOpt(i) >= QA(i)) EXIT
        END DO
      END DO
      CALL Sub_Energ_OF_ParamBasis(Norm_min,xOpt,nb_Opt,BasisnD,        &
                                   para_Tnum,mole,ComOp,para_PES,Qact)

      write(out_unitp,*) 'Initial param',xOpt,' Energy',Norm_min


      write(out_unitp,*) 'SQ',SQ
      write(out_unitp,*) 'QA',QA
      write(out_unitp,*) 'QB',QB


      Norm_max = Norm_min
      NormA    = ZERO
      print_level = -1
      ! first find the average Energy (Norm), then Temp
      DO imc=1,para_SimulatedAnnealing%nb_mc_tot/10
        DO i=1,nb_Opt
          DO k=1,1000
            CALL random_number(xi)
            xi = TWO*xi-ONE
            xOpt(i) = x0(i) + xi*SQ(i)
            IF (xOpt(i) <= QB(i) .AND. xOpt(i) >= QA(i)) EXIT
          END DO
        END DO

        CALL Sub_Energ_OF_ParamBasis(NormB,xOpt,nb_Opt,BasisnD,         &
                                     para_Tnum,mole,ComOp,para_PES,Qact)

        !write(6,*) 'NormB',imc,NormB
        NormA = NormA + NormB
        IF (NormB > Norm_max) Norm_max = NormB

        IF (NormB < Norm_min) THEN
          xOpt_min(:) = xOpt(:)
          Norm_min = NormB
        END IF
      END DO
      NormA = NormA / real(para_SimulatedAnnealing%nb_mc_tot/10,kind=Rkind)
      write(out_unitp,*) 'Min, Average, Max Energy',Norm_min,NormA,Norm_max

      x0(:)   = xOpt_min(:)
      Xopt(:) = xOpt_min(:)
      CALL Sub_Energ_OF_ParamBasis(NormB,xOpt,nb_Opt,BasisnD,         &
                                   para_Tnum,mole,ComOp,para_PES,Qact)
      Norm_min = NormB
      write(out_unitp,*) ' Energy_min',xOpt_min,Norm_min

      Temp_max = Norm_max-Norm_min
      Temp_max = NormA
      Temp     = Temp_max
      write(out_unitp,*) 'Average Energy, Temp',NormA,Temp
      CALL flush_perso(out_unitp)

      DTemp       = Temp_max/real(max(10,para_SimulatedAnnealing%nb_mc_tot),kind=Rkind)

      nb_Accepted_DNorm    = 0
      nb_Accepted_proba    = 0
      nb_NotAccepted_proba = 0
      imc                  = 0
      nb_Norm_min          = 0
      nb_block_WithoutMin  = 0
      imc_block            = 0
      !SQ(:) = SQ(:)/TEN
      write(6,*) 'NormB,SQ(1)',NormB,SQ(1)
      DO
        imc       = imc       + 1
        imc_block = imc_block + 1

        DO i=1,nb_Opt
          DO k=1,1000
            CALL random_number(xi)
            xi = TWO*xi-ONE
            xOpt(i) = x0(i) + xi*SQ(i)
            IF (xOpt(i) <= QB(i) .AND. xOpt(i) >= QA(i)) EXIT
          END DO
        END DO


        CALL Sub_Energ_OF_ParamBasis(NormA,xOpt,nb_Opt,BasisnD,         &
                                     para_Tnum,mole,ComOp,para_PES,Qact)

        !write(6,*) 'imc,Xopt,NormA,SQ',imc,xOpt,NormA,SQ(1)

        DNorm = NormA - NormB
        !write(out_unitp,*) 'Norm',imc,xOpt,NormA
        CALL flush_perso(out_unitp)

        IF ( NormA < Norm_min) THEN
          nb_Norm_min   = nb_Norm_min + 1
          xOpt_min(:)   = xOpt(:)
          Norm_min      = NormA
          !write(out_unitp,*) ' imc, Temp, nb_Energy_min, Energy_min',imc,Temp,nb_Norm_min,Norm_min
        END IF

        IF (DNorm < ZERO) THEN
          x0(:)       = xOpt(:)
          NormB       = NormA
          !write(out_unitp,*) ' accepted, DNorm',imc,DNorm
          nb_Accepted_DNorm = nb_Accepted_DNorm + 1
        ELSE
          CALL random_number(PTemp)
          IF ( PTemp < exp(-DNorm/Temp) ) THEN
            x0(:)       = xOpt(:)
            NormB       = NormA
            !write(out_unitp,*) ' imc, Temp, Norm (propa)',imc,Temp,NormA
            !write(out_unitp,*) ' accepted, Proba',imc,PTemp,exp(-DNorm/Temp)
            nb_Accepted_proba = nb_Accepted_proba + 1
          ELSE
            nb_NotAccepted_proba = nb_NotAccepted_proba + 1
            !write(out_unitp,*) ' Not accepted, Proba',imc,PTemp,exp(-DNorm/Temp)
          END IF
        END IF


        IF (imc_block > para_SimulatedAnnealing%nb_mc_partial) THEN
        SELECT CASE (para_SimulatedAnnealing%TempScheduling_type)
        CASE (1)
          !write(6,*) 'LinCoolParam'
          Temp = Temp - DTemp
        CASE (2)
          !write(6,*) 'ExpCoolParam'
          Temp = para_SimulatedAnnealing%ExpCoolParam * Temp
        CASE DEFAULT
          !write(6,*) 'ExpCoolParam'
          Temp = para_SimulatedAnnealing%ExpCoolParam * Temp
        END SELECT
        END IF

        !write(out_unitp,*) 'imc,Temp',imc,Temp,Temp_max,para_SimulatedAnnealing%ResetTempScal
        !write(out_unitp,*) 'imc,Temp',imc,Temp,Temp_max*para_SimulatedAnnealing%ResetTempScal

        !write(out_unitp,*) 'test ?',(Temp < para_SimulatedAnnealing%ResetTempScal*Temp_max)


        IF (para_SimulatedAnnealing%ResetTemp .AND.                     &
             Temp < para_SimulatedAnnealing%ResetTempScal*Temp_max) THEN

          write(out_unitp,*) 'imc,Temp_max,nb_Energy_min,Energy_min',   &
                                       imc,Temp_max,nb_Norm_min,Norm_min
          write(out_unitp,*) 'imc,Temp_max,nb_Accepted_...',imc,Temp_max,&
                nb_Accepted_DNorm,nb_Accepted_proba,nb_NotAccepted_proba
          CALL flush_perso(out_unitp)
          SQ(:)       = SQ(:)*para_SimulatedAnnealing%RangeScal
          Temp_max    = (ONE-para_SimulatedAnnealing%ResetTempScal)*Temp_max
          DTemp       = Temp_max/real(max(10,para_SimulatedAnnealing%nb_mc_tot-imc),kind=Rkind)
          Temp        = Temp_max
          x0(:)       = xOpt_min(:)
          IF (nb_Norm_min == 0) THEN
            nb_block_WithoutMin = nb_block_WithoutMin + 1
          ELSE
            nb_block_WithoutMin = 0
            nb_Norm_min         = 0
          END IF

          nb_Accepted_DNorm    = 0
          nb_Accepted_proba    = 0
          nb_NotAccepted_proba = 0
          imc_block            = 0

        END IF

        IF (Temp < para_SimulatedAnnealing%Tmin .OR. imc > para_SimulatedAnnealing%nb_mc_tot) EXIT
      END DO
      write(out_unitp,*) 'Optimal param',xOpt_min,' Energy',Norm_min

      !deallocation
      CALL dealloc_array(x,   'x',   name_sub)
      CALL dealloc_array(x0,  'x0',  name_sub)
      CALL dealloc_array(xOpt,'xOpt',name_sub)
      CALL dealloc_array(QAv, 'Qav', name_sub)
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


      END SUBROUTINE Sub_SimulatedAnnealing_cuba


      END MODULE mod_SimulatedAnnealing

