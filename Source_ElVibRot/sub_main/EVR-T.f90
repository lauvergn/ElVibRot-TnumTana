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
!%LATEX-USER-DOC-Driver
!
!\begin{itemize}
!  \item
!    Time independent calculations: energy levels, spectrum with intensities...
!  \item
!    Time dependent calculations (wavepacket propagations):
!    propagations with time dependent field, relaxation, optimal control ...
!\end{itemize}
!
!The main originality concerns the use of numerical kinetic energy operator (Tnum),
! which enables a large flexibility in the choice of the curvilinear coordinates.
!
!\section(Input file)
!
!The input file has four mains sections:
!\begin{itemize}
!  \item
!   SYSTEM and CONSTANTS, which define general parameters for parallelization, printing levels, energy unit, physical constants...
!  \item
!   COORDINATES, which defines the curvilinear coordinates, the coordinates transformations and some aspects of the physical models (constraints....). This section is part of Tnum.
!  \item
!   OPERATORS and BASIS SETS, which define parameters of scalar operators (potential, dipole moments...) and the active and inactive basis set (contracted).
!  \item
!   ANALYSIS, which defines parameters for time dependent (including optimal control) or independent calculations, intensities.
!\end{itemize}
!%END-LATEX-USER-DOC-Driver
!===========================================================================
      PROGRAM ElVibRot
      USE mod_system
!$    USE omp_lib, only : omp_get_max_threads
      USE mod_nDGridFit
      IMPLICIT NONE

      logical  :: intensity_only,analysis_only,Popenmp
      integer  :: PMatOp_omp,POpPsi_omp,PBasisTOGrid_omp,PGrid_omp,optimization
      integer  :: maxth,PMatOp_maxth,POpPsi_maxth,PBasisTOGrid_maxth,PGrid_maxth
      integer  :: PSG4_omp,PSG4_maxth
      integer (kind=ILkind)  :: max_mem
      integer  :: printlevel,err
      logical  :: test,EVR,cart,nDfit,nDGrid,mem_debug
      logical  :: GridTOBasis_test,OpPsi_test,main_test
      character (len=Name_longlen) :: EneFormat
      character (len=Name_longlen) :: RMatFormat
      character (len=Name_longlen) :: CMatFormat

      namelist /system/ max_mem,mem_debug,test,printlevel,Popenmp,      &
                          PSG4_omp,PSG4_maxth,                          &
                          PMatOp_omp,PMatOp_maxth,                      &
                          POpPsi_omp,POpPsi_maxth,                      &
                          PBasisTOGrid_omp,PBasisTOGrid_maxth,          &
                          PGrid_omp,PGrid_maxth,                        &
                          RMatFormat,CMatFormat,EneFormat,              &
                          intensity_only,analysis_only,EVR,cart,        &
                          GridTOBasis_test,OpPsi_test,                  &
                          optimization,nDfit,nDGrid,                    &
                          main_test,                                    &
                          EVRT_path,base_FileName


        intensity_only   = .FALSE.
        analysis_only    = .FALSE.
        test             = .FALSE.
        cart             = .FALSE.
        GridTOBasis_test = .FALSE.
        OpPsi_test       = .FALSE.
        EVR              = .FALSE.   ! ElVibRot (default)
        nDfit            = .FALSE.
        nDGrid           = .FALSE.
        main_test        = .FALSE.
        optimization     = 0

        maxth              = 1
        !$ maxth           = omp_get_max_threads()
        Popenmp            = .TRUE.
        PMatOp_omp         = 0
        PMatOp_maxth       = maxth
        POpPsi_omp         = 0
        POpPsi_maxth       = maxth
        PBasisTOGrid_omp   = 0
        PBasisTOGrid_maxth = maxth
        PGrid_omp          = 1
        PGrid_maxth        = maxth

        PSG4_omp           = 1
        PSG4_maxth         = maxth

        max_mem          = 4000000000_ILkind/Rkind ! 4GO
        mem_debug        = .FALSE.
        printlevel       = 0

        EneFormat        = "f18.10"
        RMatFormat       = "f18.10"
        CMatFormat       = "f15.7"

        CALL versionEVRT(.TRUE.)

        read(in_unitp,system,IOSTAT=err)
        IF (err < 0) THEN
          write(out_unitp,*) ' ERROR in ElVibRot (main program)'
          write(out_unitp,*) ' End-of-file or End-of-record'
          write(out_unitp,*) ' The namelist "system" is probably absent'
          write(out_unitp,*) ' check your data!'
          write(out_unitp,*) ' ERROR in ElVibRot (main program)'
          STOP
        ELSE IF (err > 0) THEN
          write(out_unitp,*) ' ERROR in ElVibRot (main program)'
          write(out_unitp,*) ' Some parameter name of the namelist "system" are probaly wrong'
          write(out_unitp,*) ' check your data!'
          write(out_unitp,system)
          write(out_unitp,*) ' ERROR in ElVibRot (main program)'
          STOP
        END IF

        para_mem%mem_debug = mem_debug

        !IF (analysis_only .OR. GridTOBasis_test .OR. OpPsi_test .OR. cart .OR. main_test) EVR=.FALSE.
        EVR = .NOT. (analysis_only .OR. GridTOBasis_test .OR.            &
                     OpPsi_test .OR. cart .OR. main_test .OR. nDfit .OR. &
                     nDGrid .OR. optimization /= 0 .OR. analysis_only)

        IF (printlevel > 1) write(out_unitp,system)

        para_EVRT_calc%optimization     = optimization
        para_EVRT_calc%EVR              = EVR
        para_EVRT_calc%analysis_only    = analysis_only
        para_EVRT_calc%intensity_only   = intensity_only
        para_EVRT_calc%cart             = cart
        para_EVRT_calc%GridTOBasis_test = GridTOBasis_test
        para_EVRT_calc%OpPsi_test       = OpPsi_test

        para_EVRT_calc%nDfit            = nDfit
        para_EVRT_calc%nDGrid           = nDGrid
        para_EVRT_calc%main_test        = main_test

        print_level = printlevel ! print_level is in mod_system.mod

        EneIO_format  = EneFormat
        RMatIO_format = RMatFormat
        CMatIO_format = "'('," // trim(adjustl(CMatFormat)) //      &
                    ",' +i'," // trim(adjustl(CMatFormat)) // ",')'"


        openmp              = Popenmp ! openmp is in mod_system.mod
        IF (.NOT. openmp) THEN
           MatOp_omp          = 0
           OpPsi_omp          = 0
           BasisTOGrid_omp    = 0
           Grid_omp           = 0
           SG4_omp            = 0

           MatOp_maxth        = 1
           OpPsi_maxth        = 1
           BasisTOGrid_maxth  = 1
           Grid_maxth         = 1
           SG4_maxth          = 1
        ELSE
           MatOp_omp          = PMatOp_omp
           OpPsi_omp          = POpPsi_omp
           BasisTOGrid_omp    = PBasisTOGrid_omp
           Grid_omp           = PGrid_omp
           SG4_omp            = PSG4_omp

           IF (MatOp_omp > 0) THEN
             MatOp_maxth        = min(PMatOp_maxth,maxth)
           ELSE
             MatOp_maxth        = 1
           END IF

           IF (OpPsi_omp > 0) THEN
             OpPsi_maxth        = min(POpPsi_maxth,maxth)
           ELSE
             OpPsi_maxth        = 1
           END IF

           IF (BasisTOGrid_omp > 0) THEN
             BasisTOGrid_maxth  = min(PBasisTOGrid_maxth,maxth)
           ELSE
             BasisTOGrid_maxth  = 1
           END IF

           IF (Grid_omp > 0) THEN
             Grid_maxth         = min(PGrid_maxth,maxth)
           ELSE
             Grid_maxth         = 1
           END IF

           IF (SG4_omp > 0) THEN
             SG4_maxth         = PSG4_maxth
           ELSE
             SG4_maxth         = 1
           END IF

        END IF

        write(out_unitp,*) '========================================='
        write(out_unitp,*) 'OpenMP parameters:'
        write(out_unitp,*) 'Max number of threads:           ',maxth
        write(out_unitp,*) 'MatOp_omp,      MatOp_maxth      ',MatOp_omp,MatOp_maxth
        write(out_unitp,*) 'OpPsi_omp,      OpPsi_maxth      ',OpPsi_omp,OpPsi_maxth
        write(out_unitp,*) 'BasisTOGrid_omp,BasisTOGrid_maxth',BasisTOGrid_omp,BasisTOGrid_maxth
        write(out_unitp,*) 'Grid_omp,       Grid_maxth       ',Grid_omp,Grid_maxth
        write(out_unitp,*) 'SG4_omp,        SG4_maxth        ',SG4_omp,SG4_maxth
        write(out_unitp,*) '========================================='

        para_mem%max_mem    = max_mem/Rkind
        write(out_unitp,*) '========================================='
        write(out_unitp,*) '========================================='
        IF (para_EVRT_calc%optimization /= 0) THEN
          write(out_unitp,*) ' Optimization calculation'
          write(out_unitp,*) '========================================='
          CALL sub_Optimization_OF_VibParam(max_mem)

        ELSE IF (para_EVRT_calc%nDfit .OR. para_EVRT_calc%nDGrid) THEN
          write(out_unitp,*) ' nDfit or nDGrid calculation'
          write(out_unitp,*) '========================================='
          CALL sub_nDGrid_nDfit()

        ELSE IF (para_EVRT_calc%EVR) THEN
          write(out_unitp,*) ' ElVibRot calculation'
          write(out_unitp,*) '========================================='
          CALL vib(max_mem,test,intensity_only)

        ELSE IF (para_EVRT_calc%cart) THEN
          write(out_unitp,*) ' cart calculation'
          write(out_unitp,*) '========================================='
          CALL sub_cart(max_mem)

        ELSE IF (para_EVRT_calc%GridTOBasis_test) THEN
          write(out_unitp,*) ' sub_GridTOBasis calculation'
          write(out_unitp,*) '========================================='
          CALL sub_GridTOBasis_test(max_mem)

        ELSE IF (para_EVRT_calc%OpPsi_test) THEN
          write(out_unitp,*) ' OpPsi calculation'
          write(out_unitp,*) '========================================='
          CALL Sub_OpPsi_test(max_mem)

        ELSE IF (para_EVRT_calc%analysis_only) THEN
          write(out_unitp,*) ' WP analysis calculation'
          write(out_unitp,*) '========================================='
          CALL sub_analysis_only(max_mem)

        ELSE IF (para_EVRT_calc%main_test) THEN
          write(out_unitp,*) ' Smolyat test calculation'
          write(out_unitp,*) '========================================='
          CALL sub_main_Smolyak_test()

        ELSE
          write(out_unitp,*) ' ElVibRot calculation (default)'
          write(out_unitp,*) '========================================='
          CALL vib(max_mem,test,intensity_only)
        END IF



        write(out_unitp,*) '========================================='
        write(out_unitp,*) '========================================='
        END PROGRAM ElVibRot

