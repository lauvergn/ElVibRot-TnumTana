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
SUBROUTINE init_EVR_new()
   USE mod_EVR
      IMPLICIT NONE
      logical  :: intensity_only,analysis_only,Popenmp,Popenmpi
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
      character (len=Line_len)     :: base_FileName = ''

      integer        :: ith


      ! parameters for system setup
      ! make sure to be prepared in file
      namelist /system/ max_mem,mem_debug,test,printlevel,              &
                          Popenmp,Popenmpi,                             &
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
                          EVRT_path,File_path,base_FileName

        !> initialize MPI
        !> id=0 to be the master
        !---------------------------------------------------------------------------------
#if(run_MPI)
        CALL ini_MPI()
        Popenmpi           = .TRUE.  !< True to run with MPI
        Popenmp            = .FALSE.
#endif

#if(run_openMP)
        Popenmp            = .TRUE.  !< True to run openMP
        Popenmpi           = .FALSE.
#endif

        intensity_only     = .FALSE.
        analysis_only      = .FALSE.
        test               = .FALSE.
        cart               = .FALSE.
        GridTOBasis_test   = .FALSE.
        OpPsi_test         = .FALSE.   !< True for test of action
        EVR                = .FALSE.   ! ElVibRot (default)
        nDfit              = .FALSE.
        nDGrid             = .FALSE.
        main_test          = .FALSE.
        optimization       = 0

        maxth              = 1
        !$ maxth           = omp_get_max_threads()

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
        printlevel       = -1

        EneFormat        = "f18.10"
        RMatFormat       = "f18.10"
        CMatFormat       = "f15.7"

        ! version and copyright statement
        CALL versionEVRT(.TRUE.)
        write(out_unitp,*)
        IF(Popenmpi) THEN
          CALL time_perso('MPI start, initial time')
          CALL MPI_ini_messaage()
        ENDIF

        !> automatically decide the reading of namelist, from file or shell
        !> NOTE: remember to use vib to ensure "rm namelist" to prevent the
        !> reading of old namelist
        CALL file_open2('namelist',in_unitp,old=.TRUE.,err_file=err)
        !in_unitp=10
        !open(in_unitp,file='namelist',STATUS='OLD',IOSTAT=err)
        IF(err/=0) THEN
          write(out_unitp,*) ' ERROR in init_EVR'
          write(out_unitp,*) '   "namelist" file does not exist.'
          STOP 'ERROR: namelist file does not exist'
        END IF
        read(in_unitp,system,IOSTAT=err)

        IF (err < 0) THEN
          write(out_unitp,*) ' ERROR in init_EVR'
          write(out_unitp,*) ' End-of-file or End-of-record'
          write(out_unitp,*) ' The namelist "system" is probably absent'
          write(out_unitp,*) ' check your data!'
          write(out_unitp,*) ' ERROR in init_EVR'
          STOP 'ERROR in init_EVR'
        ELSE IF (err > 0) THEN
          write(out_unitp,*) ' ERROR in init_EVR'
          write(out_unitp,*) ' Some parameter names of the namelist "system" are probaly wrong'
          write(out_unitp,*) ' check your data!'
          write(out_unitp,system)
          write(out_unitp,*) ' ERROR in init_EVR'
          STOP 'ERROR in init_EVR'
        END IF

        IF (base_FileName /= "" .AND. File_path /= "") THEN
          write(out_unitp,*) ' ERROR in init_EVR'
          write(out_unitp,*) ' base_FileName and File_path are both set!!'
          write(out_unitp,*) ' You MUST define only File_path.'
          write(out_unitp,*) ' check your data!'
          write(out_unitp,system)
          write(out_unitp,*) ' ERROR in init_EVR'
          STOP 'ERROR in init_EVR'
        ELSE IF (base_FileName /= "") THEN
          File_path = base_FileName
        END IF

        para_mem%mem_debug = mem_debug

        EVR = .NOT. (analysis_only .OR. GridTOBasis_test .OR.            &
                     OpPsi_test .OR. cart .OR. main_test .OR. nDfit .OR. &
                     nDGrid .OR. optimization /= 0 .OR. analysis_only)

       IF (.NOT. EVR) THEN
          write(out_unitp,*) ' ERROR in init_EVR'
          write(out_unitp,*) ' This subroutine is used only to initialized ElVibRot dynamics'
          write(out_unitp,*) ' => EVR=T'
          write(out_unitp,system)
          write(out_unitp,*) ' ERROR in init_EVR'
          STOP 'ERROR in init_EVR'
       END IF
       IF (intensity_only) THEN
          write(out_unitp,*) ' ERROR in init_EVR'
          write(out_unitp,*) ' This subroutine is used only to initialized ElVibRot dynamics'
          write(out_unitp,*) ' intensity_only CANNOT be true'
          write(out_unitp,system)
          write(out_unitp,*) ' ERROR in init_EVR'
          STOP 'ERROR in init_EVR'
       END IF
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
        openmpi             = Popenmpi

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

        IF(MPI_id==0 .AND. .NOT. openmpi) THEN
          write(out_unitp,*) '========================================='
          write(out_unitp,*) 'OpenMP parameters:'
          write(out_unitp,*) 'Max number of threads:           ',maxth
          write(out_unitp,*) 'MatOp_omp,      MatOp_maxth      ',MatOp_omp,MatOp_maxth
          write(out_unitp,*) 'OpPsi_omp,      OpPsi_maxth      ',OpPsi_omp,OpPsi_maxth
          write(out_unitp,*) 'BasisTOGrid_omp,BasisTOGrid_maxth',BasisTOGrid_omp,BasisTOGrid_maxth
          write(out_unitp,*) 'Grid_omp,       Grid_maxth       ',Grid_omp,Grid_maxth
          write(out_unitp,*) 'SG4_omp,        SG4_maxth        ',SG4_omp,SG4_maxth
          write(out_unitp,*) '========================================='

          write(out_unitp,*) '========================================='
          write(out_unitp,*) 'File_path: ',trim(adjustl(File_path))
          write(out_unitp,*) '========================================='
        ENDIF ! for MPI_id=0

        para_mem%max_mem    = max_mem/Rkind

        close(in_unitp)

  allocate(tab_EVRT(maxth))
  DO ith=1,maxth
    CALL file_open2('namelist',in_unitp,old=.TRUE.,err_file=err)


        IF(MPI_id==0) THEN
          write(out_unitp,*) '========================================='
          write(out_unitp,*) '========================================='
          write(out_unitp,*) ' ElVibRot calculation',ith
          write(out_unitp,*) '========================================='
        ENDIF ! for MPI_id=0

        ! from vib.f90 file

        IF(MPI_id==0) THEN
          write(out_unitp,*)
          write(out_unitp,*) '================================================='
          write(out_unitp,*) ' VIB: BEGINNING ini_data'
          CALL time_perso('ini_data ini')
          write(out_unitp,*)
        ENDIF
        !CALL system_mem_usage(memory_RSS,'before ini_data')

        CALL   ini_data(tab_EVRT(ith)%const_phys,                       &
                        tab_EVRT(ith)%para_Tnum,tab_EVRT(ith)%mole,     &
                        tab_EVRT(ith)%para_AllBasis,                    &
                        tab_EVRT(ith)%para_AllOp,tab_EVRT(ith)%para_ana,&
                        tab_EVRT(ith)%para_intensity,intensity_only,    &
                        tab_EVRT(ith)%para_propa)

        !CALL system_mem_usage(memory_RSS,'after ini_data')

        IF(MPI_id==0) THEN
          write(out_unitp,*)
          CALL time_perso('ini_data end')
          write(out_unitp,*)
          write(out_unitp,*) ' VIB: END ini_data'
          write(out_unitp,*) '================================================='
          write(out_unitp,*)
        ENDIF


        IF(MPI_id==0) THEN
          write(out_unitp,*) '========================================='
          write(out_unitp,*) ith
          write(out_unitp,*) '========================================='
        ENDIF ! for MPI_id=0
     close(in_unitp)
     CALL flush_perso(out_unitp)
  END DO

END SUBROUTINE init_EVR_new


SUBROUTINE get_nb_var_new(nb_var)
  USE mod_EVR
!$ USE omp_lib, only : OMP_GET_THREAD_NUM
  IMPLICIT NONE

  integer,          intent(out) :: nb_var


  integer :: ith

  ith=1
  !$ ith=OMP_GET_THREAD_NUM()+1

  nb_var = tab_EVRT(ith)%mole%nb_var

END SUBROUTINE get_nb_var_new
SUBROUTINE get_TnumRefGeom_Q0_new(Q0,nb_Q0,Qdyn0)
   USE mod_EVR
!$ USE omp_lib, only : OMP_GET_THREAD_NUM

  IMPLICIT NONE

  integer,          intent(in)    :: nb_Q0
  real(kind=Rkind), intent(inout) :: Q0(nb_Q0)
  logical,          intent(in)    :: Qdyn0


  integer :: ith

  ith=1
  !$ ith=OMP_GET_THREAD_NUM()+1

  IF (nb_Q0 /= tab_EVRT(ith)%mole%nb_var) THEN
    write(out_unitp,*) ' ERROR in get_TnumRefGeom_Q0'
    write(out_unitp,*) ' The size of Q0 is different from mole%nb_var'
    write(out_unitp,*) '   nb_Q0 ',nb_Q0
    write(out_unitp,*) '   nb_var',tab_EVRT(ith)%mole%nb_var
    STOP 'ERROR in get_TnumRefGeom_Q0'
  END IF

  IF (Qdyn0) THEN
    Q0(:) = tab_EVRT(ith)%mole%ActiveTransfo%Qdyn0(:)
  ELSE ! Qact order
    Q0(:) = tab_EVRT(ith)%mole%ActiveTransfo%Qact0(:)
  END IF

END SUBROUTINE get_TnumRefGeom_Q0_new
SUBROUTINE Modify_TnumRefGeom_Q0_new(Q0,nb_Q0,Qdyn0)
   USE mod_EVR
!$ USE omp_lib, only : OMP_GET_THREAD_NUM
  IMPLICIT NONE

  integer,          intent(in) :: nb_Q0
  real(kind=Rkind), intent(in) :: Q0(nb_Q0)
  logical,          intent(in) :: Qdyn0


  ! local variables
  real(kind=Rkind) :: Q1(nb_Q0)
  integer          :: ith


  ith=1
  !$ ith=OMP_GET_THREAD_NUM()+1

  IF (nb_Q0 /= tab_EVRT(ith)%mole%nb_var) THEN
    write(out_unitp,*) ' ERROR in Modify_TnumRefGeom_Q0_new'
    write(out_unitp,*) ' The size of Q0 is different from mole%nb_var'
    write(out_unitp,*) '   nb_Q0 ',nb_Q0
    write(out_unitp,*) '   nb_var',tab_EVRT(ith)%mole%nb_var
    STOP 'ERROR in Modify_TnumRefGeom_Q0'
  END IF

  IF (Qdyn0) THEN
    tab_EVRT(ith)%mole%ActiveTransfo%Qdyn0(:) =  Q0(:)
    CALL Qdyn_TO_Qact_FROM_ActiveTransfo(Q0,Q1,tab_EVRT(ith)%mole%ActiveTransfo)
    tab_EVRT(ith)%mole%ActiveTransfo%Qact0(:) =  Q1(:)
  ELSE
    tab_EVRT(ith)%mole%ActiveTransfo%Qact0(:) =  Q0(:)
    CALL Qact_TO_Qdyn_FROM_ActiveTransfo(Q0,Q1,tab_EVRT(ith)%mole%ActiveTransfo)
    tab_EVRT(ith)%mole%ActiveTransfo%Qdyn0(:) =  Q1(:)
  END IF

END SUBROUTINE Modify_TnumRefGeom_Q0_new

SUBROUTINE get_nb_nq_new(nb,nq)
   USE mod_EVR
!$ USE omp_lib, only : OMP_GET_THREAD_NUM
  IMPLICIT NONE

  integer,          intent(out) :: nb,nq

  integer :: ith

  ith=1
  !$ ith=OMP_GET_THREAD_NUM()+1

  nb = get_nb_FROM_basis(tab_EVRT(ith)%para_AllBasis%BasisnD)
  nq = get_nq_FROM_basis(tab_EVRT(ith)%para_AllBasis%BasisnD)

END SUBROUTINE get_nb_nq_new

SUBROUTINE levels_EVR_new(EigenVal,EigenVecB,EigenVecG,RhoWeight,nb,nq,nb_vec)
   USE mod_EVR
!$ USE omp_lib, only : OMP_GET_THREAD_NUM

      USE mod_FullPropa
      USE mod_FullControl
      USE mod_Davidson
      USE mod_Filter
      USE mod_Arpack
      USE mod_Op
      USE mod_analysis
      USE mod_fullanalysis
      USE mod_Auto_Basis
      USE mod_psi
      USE mod_MPI_aux
      IMPLICIT NONE

      integer,           intent(in)    :: nb,nq             ! numbers of basis functions and grid points
      integer,           intent(inout) :: nb_vec            ! number of eigenvectors/eigenvalues
      real (kind=Rkind), intent(inout) :: EigenVal(nb)      ! eigenvalues.               EigenVal(i)
      real (kind=Rkind), intent(inout) :: EigenVecB(nb,nb)  ! eigenvectors on the basis. EigenVecB(:,i)
      real (kind=Rkind), intent(inout) :: EigenVecG(nq,nb)  ! eigenvectors on the grid.  EigenVecG(:,i)
      real (kind=Rkind), intent(inout) :: RhoWeight(nq)     ! rho(Q).Weight(Q), on the grid points.

!----- variables for the construction of H ---------------------------------------------
      TYPE (param_Op), pointer    :: para_H      => null()
      TYPE (param_Op), pointer    :: para_S      => null()
      TYPE (param_Op), pointer    :: para_Dip(:) => null()
      integer                     :: iOp
      real (kind=Rkind)           :: max_Sii,max_Sij

!----- variables for the WP propagation ------------------------------------------------
      TYPE (param_psi)            :: WP0tmp,MuWP0

!----- for Davidson diagonalization ----------------------------------------------------
      integer                     :: nb_diago
      integer                     :: max_diago
      TYPE (param_psi),  allocatable  :: Tab_Psi(:)
      real (kind=Rkind), allocatable  :: Ene0(:)

      integer :: Get_nbPERsym_FROM_SymAbelianOFAllBasis ! function

!----- variables divers ----------------------------------------------------------------
      integer           :: i,ip,i_baie,f_baie,id,nb_ScalOp
      real (kind=Rkind) :: T,DE,Ep,Em,Q,fac,zpe,pop
      logical           :: print_mat
      integer           :: err
      integer           :: err_mem,memory
      real (kind=Rkind) :: part_func ! function
      integer           :: ith


!para_mem%mem_debug=.TRUE.
!---------------------------------------------------------------------------------------

  ith=1
  !$ ith=OMP_GET_THREAD_NUM()+1
  write(out_unitp,*) 'ith',ith ; flush(out_unitp)

      para_H => tab_EVRT(ith)%para_AllOp%tab_Op(1)

!---------------------------------------------------------------------------------------
!      Grids (V, T, Dip) calculations
!---------------------------------------------------------------------------------------
      !> turn off the allocation of grid, to be done in action
      !> sub_qa_bhe should be refined later
!#if(run_MPI)
!      Grid_allco=.FALSE.
!#endif
      IF(MPI_id==0 .AND. print_level > -1) THEN
        write(out_unitp,*)
        write(out_unitp,*) '================================================='
        write(out_unitp,*) ' VIB: BEGINNING sub_qa_bhe'
        CALL time_perso('sub_qa_bhe ini')
        write(out_unitp,*)
      ENDIF

      IF (tab_EVRT(ith)%para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Read_FileGrid) THEN ! test only for H
        IF(MPI_id==0 .AND. print_level > -1) THEN
          write(out_unitp,*)
          write(out_unitp,*) ' The grid of operators (S V Veff T1 and T2) will be read'
          write(out_unitp,*)
        ENDIF
        IF (tab_EVRT(ith)%para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Test_Grid) THEN
          IF(MPI_id==0) write(out_unitp,*) 'Test_Grid=.TRUE. => STOP in vib'
          STOP
        END IF
      END IF

      IF (tab_EVRT(ith)%para_ana%VibRot) THEN
         tab_EVRT(ith)%para_Tnum%JJ = tab_EVRT(ith)%para_ana%JJmax
         tab_EVRT(ith)%para_Tnum%With_Cart_Transfo = (tab_EVRT(ith)%para_Tnum%JJ>0) .AND. tab_EVRT(ith)%mole%Cart_transfo
      END IF

      ! calculate potential in action
      CALL sub_qa_bhe(tab_EVRT(ith)%para_AllOp)

      IF (tab_EVRT(ith)%para_ana%VibRot) THEN
         tab_EVRT(ith)%para_Tnum%JJ = 0
         tab_EVRT(ith)%para_Tnum%With_Cart_Transfo = (tab_EVRT(ith)%para_Tnum%JJ>0) .AND. tab_EVRT(ith)%mole%Cart_transfo
      END IF

      IF(MPI_id==0 .AND. print_level > -1) THEN
        write(out_unitp,*)
        CALL time_perso('sub_qa_bhe end')
        write(out_unitp,*) ' VIB: END sub_qa_bhe'
        write(out_unitp,*) '================================================='
        write(out_unitp,*)
      ENDIF

!#if(run_MPI)
!      Grid_allco=.True.
!#endif
!---------------------------------------------------------------------------------------
!      contraction of the active basis set with HADA basis
!---------------------------------------------------------------------------------------
      IF (para_H%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC .AND. para_H%nb_bi>1) THEN
        IF(MPI_id==0 .AND. print_level > -1) THEN
          write(out_unitp,*)
          write(out_unitp,*) '================================================='
          write(out_unitp,*) ' VIB: BEGINNING HADA contraction',          &
                                   tab_EVRT(ith)%para_AllBasis%BasisnD%nb
          CALL time_perso('HADA contraction')
          write(out_unitp,*)
          write(out_unitp,*)
        ENDIF

        CALL sub_MatOp_HADA(para_H,tab_EVRT(ith)%para_ana,                  &
                        tab_EVRT(ith)%para_intensity,tab_EVRT(ith)%para_AllOp,  &
                                                   tab_EVRT(ith)%const_phys)

        IF(MPI_id==0 .AND. print_level > -1) THEN
          write(out_unitp,*)
          write(out_unitp,*)
          CALL time_perso('HADA contraction')
          write(out_unitp,*) ' VIB: END HADA contraction'
          write(out_unitp,*) '================================================='
          write(out_unitp,*)
        ENDIF

        tab_EVRT(ith)%para_AllOp%tab_Op(:)%nb_tot     =                 &
               sum(para_H%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(:)) *  &
               para_H%para_ReadOp%nb_elec
        tab_EVRT(ith)%para_AllOp%tab_Op(:)%nb_tot_ini = tab_EVRT(ith)%para_AllOp%tab_Op(:)%nb_tot
      END IF

!=====================================================================
!       => Time-independent calculation
!=====================================================================

        !================================================================
        !================================================================
        !================================================================
        !===== Tune the number of threads (for SG4) =====================
        !================================================================
        max_diago = max(10,tab_EVRT(ith)%para_propa%para_Davidson%nb_WP,para_H%nb_tot/10)
        max_diago = min(max_diago,10,para_H%nb_tot)
        !CALL Tune_SG4threads_HPsi(para_H%cplx,max_diago,para_H)

        !================================================================
        !===== build S and/or H if necessary ============================
        !================================================================
        IF (para_H%Make_Mat) THEN
          IF(MPI_id==0 .AND. print_level > -1) THEN
            write(out_unitp,*)
            write(out_unitp,*) '================================================='
            write(out_unitp,*) ' VIB: BEGINNING sub_matOp',para_H%nb_tot
            CALL time_perso('sub_matOp: H and S')
            write(out_unitp,*)
            write(out_unitp,*) 'para_S...%comput_S',tab_EVRT(ith)%para_AllOp%tab_Op(2)%para_ReadOp%comput_S
            write(out_unitp,*)
          ENDIF

          IF (tab_EVRT(ith)%para_AllOp%tab_Op(2)%para_ReadOp%comput_S) THEN

            para_S => tab_EVRT(ith)%para_AllOp%tab_Op(2)
            CALL sub_MatOp(para_S,tab_EVRT(ith)%para_ana%print)

            !- analysis of the overlap matrix
            CALL sub_ana_S(para_S%Rmat,para_S%nb_tot,max_Sii,max_Sij,.TRUE.)

            CALL dealloc_para_Op(para_S)
            nullify(para_S)
          END IF ! for para_AllOp%tab_Op(2)%para_ReadOp%comput_S

          CALL sub_MatOp(para_H,tab_EVRT(ith)%para_ana%print)

          IF(MPI_id==0 .AND. print_level > -1) THEN
            write(out_unitp,*)
            write(out_unitp,*)
            CALL time_perso('sub_matOp: H and S')
            write(out_unitp,*) ' VIB: END sub_matOp'
            write(out_unitp,*) '================================================='
            write(out_unitp,*)
          ENDIF
        END IF ! for para_H%Make_Mat
        !================================================================
        !================================================================

        !================================================================
        !===== Hmax calculation (for filter diagonalization)
        !================================================================
        IF (tab_EVRT(ith)%para_ana%filter .AND. tab_EVRT(ith)%para_propa%auto_Hmax) THEN
          IF(MPI_id==0 .AND. print_level > -1) THEN
            write(out_unitp,*)
            write(out_unitp,*) '================================================'
            write(out_unitp,*) ' VIB: Hmax and Hmin calculation'
            CALL time_perso('sub_Hmax ini2')
            write(out_unitp,*)
            write(out_unitp,*)
          ENDIF

        IF(.NOT. openmpi) CALL sub_Hmax(tab_EVRT(ith)%para_propa,para_H)

          IF(MPI_id==0 .AND. print_level > -1) THEN
            write(out_unitp,*)
            write(out_unitp,*)
            CALL time_perso('sub_Hmax end2')
            write(out_unitp,*) ' VIB: END Hmax and Hmin calculation'
            write(out_unitp,*) '================================================'
            write(out_unitp,*)
          ENDIF
        END IF ! for para_ana%filter .AND. para_propa%auto_Hmax

        !================================================================
        !===== Diagonalisation ==========================================
        !================================================================
        IF (tab_EVRT(ith)%para_ana%davidson .OR. tab_EVRT(ith)%para_ana%arpack .OR. tab_EVRT(ith)%para_ana%filter) THEN
          IF(MPI_id==0 .AND. print_level > -1) THEN
            write(out_unitp,*)
            write(out_unitp,*) '================================================='
            write(out_unitp,*) ' VIB: BEGINNING ITERATIVE DIAGONALIZATION'
            CALL time_perso('sub_Iterative_Diago')
            write(out_unitp,*)
          ENDIF

          IF (tab_EVRT(ith)%para_propa%para_Davidson%max_WP == 0) THEN
            max_diago = max(1000,tab_EVRT(ith)%para_propa%para_Davidson%nb_WP,          &
                            para_H%nb_tot/10)
          ELSE
            max_diago = tab_EVRT(ith)%para_propa%para_Davidson%max_WP
          END IF
          IF (Get_nbPERsym_FROM_SymAbelianOFAllBasis(tab_EVRT(ith)%para_AllBasis,       &
                               tab_EVRT(ith)%para_propa%para_Davidson%symab) == 0) THEN
            max_diago = min(max_diago,para_H%nb_tot)
          ELSE
            max_diago = min(max_diago,para_H%nb_tot,                      &
                  Get_nbPERsym_FROM_SymAbelianOFAllBasis(tab_EVRT(ith)%para_AllBasis, &
                                        tab_EVRT(ith)%para_propa%para_Davidson%symab))
          END IF
          tab_EVRT(ith)%para_propa%para_Davidson%max_WP = max_diago

          nb_diago = min(tab_EVRT(ith)%para_propa%para_Davidson%nb_WP,para_H%nb_tot,max_diago)
!#if(run_MPI)
!          CALL MPI_Bcast(nb_diago,size1_MPI,Int_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
!#endif
          IF(openmpi .AND. MPI_scheme/=1) CALL MPI_Bcast_(nb_diago,size1_MPI,root_MPI)

          CALL alloc_NParray(Tab_Psi,(/ max_diago /),"Tab_Psi","vib")
          CALL alloc_NParray(Ene0,(/max_diago/),"Ene0","vib")

          IF (tab_EVRT(ith)%para_ana%davidson) THEN

            CALL sub_propagation_Davidson(Tab_Psi,Ene0,nb_diago,max_diago,             &
                                         para_H,tab_EVRT(ith)%para_propa%para_Davidson,&
                                         tab_EVRT(ith)%para_propa)

          ELSE IF (tab_EVRT(ith)%para_ana%arpack) THEN ! arpack=t
            !CALL sub_propagation_Arpack(Tab_Psi,Ene0,nb_diago,max_diago,  &
            !                            para_H,tab_EVRT(ith)%para_propa)
            CALL sub_propagation_Arpack_Sym(Tab_Psi,Ene0,nb_diago,max_diago,&
                                           para_H,tab_EVRT(ith)%para_propa)

          ELSE ! filter diagonalization
            CALL sub_GaussianFilterDiagonalization(Tab_Psi,Ene0,nb_diago,max_diago,&
                                           para_H,tab_EVRT(ith)%para_propa)

            !CALL sub_GaussianFilterDiagonalization_v0(Tab_Psi,Ene0,nb_diago,max_diago,&
            !                               para_H,tab_EVRT(ith)%para_propa)

            !CALL sub_FilterDiagonalization(Tab_Psi,Ene0,nb_diago,max_diago,&
            !                               para_H,tab_EVRT(ith)%para_propa)

          END IF

          CALL dealloc_NParray(Ene0,"Ene0","vib")
          para_H%spectral = .TRUE.

          IF(MPI_id==0 .AND. print_level > -1) THEN
            write(out_unitp,*)
            write(out_unitp,*)
            CALL time_perso('sub_Iterative_Diago')
            write(out_unitp,*) ' VIB: END ITERATIVE DIAGONALIZATION'
            write(out_unitp,*) '================================================='
            write(out_unitp,*)
          ENDIF

        ELSE ! for para_ana%davidson .OR. para_ana%arpack .OR. para_ana%filter
          IF(MPI_id==0 .AND. print_level > -1) THEN
            write(out_unitp,*)
            write(out_unitp,*) '================================================='
            write(out_unitp,*) ' VIB: BEGINNING DIAGONALIZATION'
            CALL time_perso('sub_diago_H')
            write(out_unitp,*)
            write(out_unitp,*)
          ENDIF

          nb_diago  = para_H%nb_tot
          max_diago = para_H%nb_tot

          IF(keep_MPI) CALL alloc_NParray(Tab_Psi,(/ max_diago /),"Tab_Psi","vib")

          IF(keep_MPI) THEN
            DO i=1,max_diago
              CALL init_psi(Tab_psi(i),para_H,para_H%cplx)
            END DO
          ENDIF

          para_H%diago = .TRUE.
          IF(keep_MPI) CALL alloc_para_Op(para_H,Grid=.FALSE.,Mat=.TRUE.)
          IF (para_H%cplx) THEN
            IF(keep_MPI) THEN
              CALL sub_diago_CH(para_H%Cmat,para_H%Cdiag,para_H%Cvp,      &
                                para_H%nb_tot)

              nb_diago = count(abs(para_H%Cdiag(:)-para_H%Cdiag(1))< tab_EVRT(ith)%para_ana%max_ene)
              nb_diago = min(nb_diago,tab_EVRT(ith)%para_ana%max_ana)
              IF (nb_diago < 3) nb_diago = 10
              IF (nb_diago > para_H%nb_tot) nb_diago = para_H%nb_tot
              DO i=1,nb_diago
                CALL alloc_psi(Tab_Psi(i))
                Tab_Psi(i)%CvecB(:)    = para_H%Cvp(:,i)
                Tab_psi(i)%CAvOp       = para_H%Cdiag(i)
                Tab_psi(i)%IndAvOp     = para_H%n_Op  ! it should be 0
                Tab_psi(i)%convAvOp    = .TRUE.
              END DO
              para_H%para_AllBasis%basis_ext%Cvp_spec    => para_H%Cvp
            ENDIF
          ELSE
            IF(keep_MPI) THEN
              CALL sub_diago_H(para_H%Rmat,para_H%Rdiag,para_H%Rvp,       &
                               para_H%nb_tot,para_H%sym_Hamil)

              IF(print_level > -1) THEN
                write(out_unitp,*) 'HMin,HMax (ua)  : ',(/ minval(para_H%Rdiag),maxval(para_H%Rdiag) /)
                write(out_unitp,*) 'HMin,HMax (cm-1): ',   &
                  (/ minval(para_H%Rdiag),maxval(para_H%Rdiag) /)*get_Conv_au_TO_unit('E','cm-1')
              END IF

              nb_diago = count((para_H%Rdiag(:)-para_H%Rdiag(1))< tab_EVRT(ith)%para_ana%max_ene)
              IF (tab_EVRT(ith)%para_ana%max_ana > 0) nb_diago = min(nb_diago,tab_EVRT(ith)%para_ana%max_ana)
              IF (nb_diago < 3) nb_diago = 10
              IF (nb_diago > para_H%nb_tot) nb_diago = para_H%nb_tot
              DO i=1,nb_diago
                CALL alloc_psi(Tab_Psi(i))
                Tab_Psi(i)%RvecB(:)    = para_H%Rvp(:,i)
                Tab_psi(i)%CAvOp       = para_H%Rdiag(i)
                Tab_psi(i)%IndAvOp     = para_H%n_Op  ! it should be 0
                Tab_psi(i)%convAvOp    = .TRUE.
                CALL Set_symab_OF_psiBasisRep(Tab_psi(i))
              END DO
              para_H%para_AllBasis%basis_ext%Rvp_spec    => para_H%Rvp
            ENDIF ! for MPI_id=0
          END IF ! for para_H%cplx
          para_H%para_AllBasis%basis_ext%nb_vp_spec  = nb_diago

          IF (associated(tab_EVRT(ith)%para_AllBasis%basis_ext%liste_spec)) THEN
            CALL dealloc_NRarray(tab_EVRT(ith)%para_AllBasis%basis_ext%liste_spec,&
                                                            "basis_ext%liste_spec","vib")
          END IF
          CALL alloc_NRarray(tab_EVRT(ith)%para_AllBasis%basis_ext%liste_spec,&
                                             [nb_diago],"basis_ext%liste_spec","vib")
          tab_EVRT(ith)%para_AllBasis%basis_ext%liste_spec(:) = (/ (i,i=1,nb_diago) /)

          IF (associated(para_H%Rmat)) THEN
            CALL dealloc_array(para_H%Rmat,"para_H%Rmat","vib")
          END IF
          IF (associated(para_H%Cmat)) THEN
            CALL dealloc_array(para_H%Cmat,"para_H%Cmat","vib")
            nullify(para_H%Cmat)
          END IF
          tab_EVRT(ith)%para_ana%max_ana = nb_diago

          IF(print_level > -1) THEN
            write(out_unitp,*)
            write(out_unitp,*)
            CALL time_perso('sub_diago_H')
            write(out_unitp,*) ' VIB: END DIAGONALIZATION'
            write(out_unitp,*) '================================================='
            write(out_unitp,*)
          ENDIF
        END IF ! for para_ana%davidson .OR. para_ana%arpack .OR. para_ana%filter
        CALL flush_perso(out_unitp)
        !===============================================================
        !===============================================================

         IF (nb /= get_nb_FROM_basis(tab_EVRT(ith)%para_AllBasis%BasisnD) .OR.        &
             nq /= get_nq_FROM_basis(tab_EVRT(ith)%para_AllBasis%BasisnD)) THEN

           write(out_unitp,*) ' ERROR in levels_EVR_new'
           write(out_unitp,*) ' inconsistant nb values',nb,get_nb_FROM_basis(tab_EVRT(ith)%para_AllBasis%BasisnD)
           write(out_unitp,*) '   or'
           write(out_unitp,*) ' inconsistant nq values',nq,get_nq_FROM_basis(tab_EVRT(ith)%para_AllBasis%BasisnD)
           STOP 'ERROR in levels_EVR'

         END IF

         nb_vec = nb_diago

         DO i=1,nb_diago
           EigenVal(i)    = para_H%Rdiag(i)
           EigenVecB(:,i) = Tab_Psi(i)%RvecB(:)

           CALL sub_PsiBasisRep_TO_GridRep(Tab_Psi(i))
           EigenVecG(:,i) = Tab_Psi(i)%RvecG(:)
           CALL alloc_psi(Tab_Psi(i),BasisRep=.TRUE.,GridRep=.FALSE.)

         END DO

         !- calculation of WrhonD ------------------------------
         DO i=1,nq
           RhoWeight(i) = Rec_WrhonD(tab_EVRT(ith)%para_AllBasis%BasisnD,i)
         END DO

        !===============================================================
        !===============================================================
        IF(MPI_id==0 .AND. print_level > -1) THEN
          write(out_unitp,*)
          write(out_unitp,*) '================================================='
          write(out_unitp,*) ' VIB: BEGINNING WAVE FUNCTION ANALYSIS'
          CALL time_perso('sub_analyse ini')
          write(out_unitp,*)
        ENDIF

        IF(keep_MPI) CALL sub_analyse(Tab_Psi,nb_diago,para_H,         &
                           tab_EVRT(ith)%para_ana,tab_EVRT(ith)%para_intensity, &
                               tab_EVRT(ith)%para_AllOp,tab_EVRT(ith)%const_phys)
        CALL flush_perso(out_unitp)

        IF (.NOT. para_H%cplx .AND. tab_EVRT(ith)%para_ana%VibRot) THEN
          CALL sub_VibRot(Tab_Psi,tab_EVRT(ith)%para_ana%max_ana,para_H,tab_EVRT(ith)%para_ana)
        END IF
        CALL flush_perso(out_unitp)

        !===============================================================
        ! Spectral representation of operator
        !===============================================================
        print_mat = (tab_EVRT(ith)%para_ana%MaxWP_TO_Write_MatOp >= nb_diago)
        IF (tab_EVRT(ith)%para_ana%Spectral_ScalOp) THEN
          DO iOp=1,size(tab_EVRT(ith)%para_AllOp%tab_Op)

            IF (tab_EVRT(ith)%para_AllOp%tab_Op(iOp)%n_Op == -1) CYCLE ! S
            IF (tab_EVRT(ith)%para_AllOp%tab_Op(iOp)%spectral) THEN
              write(out_unitp,*) '==========================================='
              write(out_unitp,*) ' Spectral representation of: ',        &
                                 trim(tab_EVRT(ith)%para_AllOp%tab_Op(iOp)%name_Op)

              CALL sub_build_MatOp(Tab_Psi,nb_diago,                     &
                                   tab_EVRT(ith)%para_AllOp%tab_Op(iOp),.TRUE.,print_mat)

              write(out_unitp,*) '==========================================='
              CALL flush_perso(out_unitp)
            END IF
          END DO
        END IF

        IF(MPI_id==0 .AND. print_level > -1) THEN
          write(out_unitp,*)
          write(out_unitp,*)
          CALL time_perso('sub_analyse end')
          write(out_unitp,*) ' VIB: END WAVE FUNCTION ANALYSIS'
          write(out_unitp,*) '================================================='
          write(out_unitp,*)
        ENDIF

        !===============================================================
        !===============================================================

        !===============================================================
        ! Deallocation of Tab_Psi
        !===============================================================
        IF(keep_MPI) THEN
          DO i=1,size(Tab_Psi)
            CALL dealloc_psi(Tab_Psi(i))
          END DO
          CALL dealloc_NParray(Tab_Psi,"Tab_Psi","vib")
        ENDIF
        !===============================================================
        !===============================================================

!=====================================================================
!       => end Time-independent calculation
!=====================================================================
!=====================================================================
!=====================================================================

      IF (.NOT. para_H%cplx .AND. tab_EVRT(ith)%para_ana%intensity) THEN
        IF(MPI_id==0 .AND. print_level > -1) THEN
          write(out_unitp,*)
          write(out_unitp,*) '================================================'
          write(out_unitp,*) ' VIB: BEGINNING sub_intensity',para_H%nb_tot,tab_EVRT(ith)%para_ana%max_ana
          CALL time_perso('sub_intensity')
          write(out_unitp,*)
          write(out_unitp,*)
        ENDIF
!         ----- file for restart or for changed the temprature --------------
        CALL file_open(tab_EVRT(ith)%para_intensity%file_resart_int,tab_EVRT(ith)%nio_res_int)

        IF (tab_EVRT(ith)%intensity_only) THEN
          read(tab_EVRT(ith)%nio_res_int,*,IOSTAT=err) para_H%nb_tot,tab_EVRT(ith)%para_ana%max_ana
          read(tab_EVRT(ith)%nio_res_int,*,IOSTAT=err)

          para_H%diago = .TRUE.
          CALL alloc_para_Op(para_H)

          CALL Read_Vec(para_H%Rdiag,tab_EVRT(ith)%nio_res_int,5,err)
          IF (err /= 0) THEN
            write(out_unitp,*) 'ERROR in levels_EVR_new'
            write(out_unitp,*) ' reading the vector "para_H%Rdiag"'
            STOP
          END IF

          IF (tab_EVRT(ith)%intensity_only) THEN
            write(out_unitp,*)
            Q =  part_func(para_H%Rdiag,size(para_H%Rdiag),tab_EVRT(ith)%para_ana%Temp,tab_EVRT(ith)%const_phys)
            fac = tab_EVRT(ith)%const_phys%Eh / (tab_EVRT(ith)%const_phys%k * tab_EVRT(ith)%para_ana%Temp)
            zpe = minval(para_H%Rdiag(:))
            write(out_unitp,*) 'population at T, Q',tab_EVRT(ith)%para_ana%Temp,Q
            write(out_unitp,*) 'Energy level (',tab_EVRT(ith)%const_phys%ene_unit,') pop:'
            DO i=1,size(para_H%Rdiag)
              pop = exp(-(para_H%Rdiag(i)-zpe)*fac)
              write(out_unitp,10) i,para_H%Rdiag(i) * tab_EVRT(ith)%const_phys%auTOenergy,&
                    (para_H%Rdiag(i) - zpe) * tab_EVRT(ith)%const_phys%auTOenergy,pop/Q
10             format(i4,x,3f20.5)
            END DO
          END IF

          read(tab_EVRT(ith)%nio_res_int,*,IOSTAT=err)

          CALL Read_Mat(para_H%Rvp,tab_EVRT(ith)%nio_res_int,5,err)
          IF (err /= 0) THEN
            write(out_unitp,*) 'ERROR in levels_EVR_new'
            write(out_unitp,*) ' reading the matrix "para_H%Rvp"'
            STOP
          END IF
        ELSE ! for intensity_only
          write(out_unitp,*) 'write restart file for intensity: ',      &
                          tab_EVRT(ith)%para_intensity%file_resart_int%name
          CALL flush_perso(out_unitp)
          write(tab_EVRT(ith)%nio_res_int,*) para_H%nb_tot,tab_EVRT(ith)%para_ana%max_ana
          write(tab_EVRT(ith)%nio_res_int,*) 'ene'
          CALL Write_Vec(para_H%Rdiag,tab_EVRT(ith)%nio_res_int,5,Rformat='e30.23')
          write(tab_EVRT(ith)%nio_res_int,*) 'psi'
          CALL flush_perso(out_unitp)
          CALL Write_Mat(para_H%Rvp,tab_EVRT(ith)%nio_res_int,5,Rformat='e30.23')
          CALL flush_perso(tab_EVRT(ith)%nio_res_int)
        END IF ! for intensity_only
        CALL flush_perso(out_unitp)
!         -------------------------------------------------------------------

        iOp = 2
        para_Dip => tab_EVRT(ith)%para_AllOp%tab_Op(iOp+1:iOp+3)
        CALL sub_intensity(para_Dip,                                    &
            tab_EVRT(ith)%para_ana%print,para_H,tab_EVRT(ith)%para_ana%max_ana, &
                          tab_EVRT(ith)%para_intensity,tab_EVRT(ith)%const_phys,&
                          tab_EVRT(ith)%intensity_only,tab_EVRT(ith)%nio_res_int)

        IF(MPI_id==0 .AND. print_level > -1) THEN
          write(out_unitp,*)
          write(out_unitp,*)
          CALL time_perso('sub_intensity')
          write(out_unitp,*) ' VIB: END sub_intensity'
          write(out_unitp,*) '================================================'
          write(out_unitp,*)
        ENDIF
        close(tab_EVRT(ith)%nio_res_int)
        nullify(para_Dip)
      END IF !for .NOT. para_H%cplx .AND. para_ana%intensity
      CALL flush_perso(out_unitp)

      IF (.NOT. para_H%cplx .AND. tab_EVRT(ith)%para_ana%Psi_ScalOp) THEN
        IF(MPI_id==0 .AND. print_level > -1) THEN
          write(out_unitp,*)
          write(out_unitp,*) '================================================'
          write(out_unitp,*) ' VIB: BEGINNING sub_AnalysePsy_ScalOp',para_H%nb_tot,tab_EVRT(ith)%para_ana%max_ana
          CALL time_perso('sub_AnalysePsy_ScalOp')
          write(out_unitp,*)
          write(out_unitp,*) 'nb_scalar_Op',para_H%para_ReadOp%nb_scalar_Op
        ENDIF

        iOp = 2
        nb_ScalOp = para_H%para_ReadOp%nb_scalar_Op
        para_Dip => tab_EVRT(ith)%para_AllOp%tab_Op(iOp+1:iOp+nb_ScalOp)
        CALL sub_AnalysePsy_ScalOp(para_Dip,nb_ScalOp,para_H,tab_EVRT(ith)%para_ana%max_ana)

        IF(MPI_id==0 .AND. print_level > -1) THEN
          write(out_unitp,*)
          write(out_unitp,*)
          CALL time_perso('sub_AnalysePsy_ScalOp')
          write(out_unitp,*) ' VIB: END sub_AnalysePsy_ScalOp'
          write(out_unitp,*) '================================================'
          write(out_unitp,*)
        ENDIF
        nullify(para_Dip)
      END IF
      CALL flush_perso(out_unitp)

      IF (.NOT. para_H%cplx .AND. tab_EVRT(ith)%para_ana%NLO) THEN
        IF(MPI_id==0 .AND. print_level > -1) THEN
          write(out_unitp,*)
          write(out_unitp,*) '================================================'
          write(out_unitp,*) ' VIB: BEGINNING sub_NLO',para_H%nb_tot,tab_EVRT(ith)%para_ana%max_ana
          CALL time_perso('sub_NLO')
          write(out_unitp,*)
          write(out_unitp,*)
        ENDIF

        iOp = 2
        para_Dip => tab_EVRT(ith)%para_AllOp%tab_Op(iOp+1:iOp+3)
        CALL sub_NLO(para_Dip,tab_EVRT(ith)%para_ana%print,para_H,tab_EVRT(ith)%para_ana%max_ana, &
                     tab_EVRT(ith)%para_intensity)

        IF(MPI_id==0 .AND. print_level > -1) THEN
          write(out_unitp,*)
          write(out_unitp,*)
          CALL time_perso('sub_NLO')
          write(out_unitp,*) ' VIB: END sub_NLO'
          write(out_unitp,*) '================================================'
          write(out_unitp,*)
        ENDIF
        nullify(para_Dip)
      END IF
      CALL flush_perso(out_unitp)

!=====================================================================
!=====================================================================

!=====================================================================
!       deallocated memories
!=====================================================================
      IF ( allocated(Tab_Psi) ) THEN
        DO i=1,size(Tab_Psi)
          CALL dealloc_psi(Tab_Psi(i))
        END DO
        CALL dealloc_NParray(Tab_Psi,"Tab_Psi","vib")
      END IF

      IF (associated(tab_EVRT(ith)%para_AllOp%tab_Op)) THEN
        DO i=1,size(tab_EVRT(ith)%para_AllOp%tab_Op)
          CALL dealloc_para_Op(tab_EVRT(ith)%para_AllOp%tab_Op(i),keep_init=.TRUE.)
        END DO
      END IF

      IF(print_level > -1) THEN
        write(out_unitp,*) 'mem_tot,max_mem_used',para_mem%mem_tot,para_mem%max_mem_used
        write(out_unitp,*) 'nb_alloc,nb_dealloc',para_mem%nb_alloc,para_mem%nb_dealloc
        write(out_unitp,*) '================================================'
        IF(openmpi) THEN
          write(out_unitp,*) ' ElVibRot-Tnum AU REVOIR!!!', ' from ', MPI_id
        ELSE
          write(out_unitp,*) ' ElVibRot-Tnum AU REVOIR!!!'
        ENDIF
        write(out_unitp,*) '================================================'
      END IF

   END SUBROUTINE levels_EVR_new


SUBROUTINE init_EVR()
   USE mod_EVR
      IMPLICIT NONE
      logical  :: intensity_only,analysis_only,Popenmp,Popenmpi
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
      character (len=Line_len)     :: base_FileName = ''

      integer        :: ith


      ! parameters for system setup
      ! make sure to be prepared in file
      namelist /system/ max_mem,mem_debug,test,printlevel,              &
                          Popenmp,Popenmpi,                             &
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
                          EVRT_path,File_path,base_FileName

        !> initialize MPI
        !> id=0 to be the master
        !---------------------------------------------------------------------------------
#if(run_MPI)
        CALL ini_MPI()
        Popenmpi           = .TRUE.  !< True to run with MPI
        Popenmp            = .FALSE.
#endif

#if(run_openMP)
        Popenmp            = .TRUE.  !< True to run openMP
        Popenmpi           = .FALSE.
#endif

        intensity_only     = .FALSE.
        analysis_only      = .FALSE.
        test               = .FALSE.
        cart               = .FALSE.
        GridTOBasis_test   = .FALSE.
        OpPsi_test         = .FALSE.   !< True for test of action
        EVR                = .FALSE.   ! ElVibRot (default)
        nDfit              = .FALSE.
        nDGrid             = .FALSE.
        main_test          = .FALSE.
        optimization       = 0

        maxth              = 1
        !$ maxth           = omp_get_max_threads()

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
        printlevel       = -1

        EneFormat        = "f18.10"
        RMatFormat       = "f18.10"
        CMatFormat       = "f15.7"

        ! version and copyright statement
        CALL versionEVRT(.TRUE.)
        write(out_unitp,*)
        IF(Popenmpi) THEN
          CALL time_perso('MPI start, initial time')
          CALL MPI_ini_messaage()
        ENDIF

        !> automatically decide the reading of namelist, from file or shell
        !> NOTE: remember to use vib to ensure "rm namelist" to prevent the
        !> reading of old namelist
        CALL file_open2('namelist',in_unitp,old=.TRUE.,err_file=err)
        !in_unitp=10
        !open(in_unitp,file='namelist',STATUS='OLD',IOSTAT=err)
        IF(err/=0) THEN
          write(out_unitp,*) ' ERROR in init_EVR'
          write(out_unitp,*) '   "namelist" file does not exist.'
          STOP 'ERROR: namelist file does not exist'
        END IF
        read(in_unitp,system,IOSTAT=err)

        IF (err < 0) THEN
          write(out_unitp,*) ' ERROR in init_EVR'
          write(out_unitp,*) ' End-of-file or End-of-record'
          write(out_unitp,*) ' The namelist "system" is probably absent'
          write(out_unitp,*) ' check your data!'
          write(out_unitp,*) ' ERROR in init_EVR'
          STOP 'ERROR in init_EVR'
        ELSE IF (err > 0) THEN
          write(out_unitp,*) ' ERROR in init_EVR'
          write(out_unitp,*) ' Some parameter names of the namelist "system" are probaly wrong'
          write(out_unitp,*) ' check your data!'
          write(out_unitp,system)
          write(out_unitp,*) ' ERROR in init_EVR'
          STOP 'ERROR in init_EVR'
        END IF

        IF (base_FileName /= "" .AND. File_path /= "") THEN
          write(out_unitp,*) ' ERROR in init_EVR'
          write(out_unitp,*) ' base_FileName and File_path are both set!!'
          write(out_unitp,*) ' You MUST define only File_path.'
          write(out_unitp,*) ' check your data!'
          write(out_unitp,system)
          write(out_unitp,*) ' ERROR in init_EVR'
          STOP 'ERROR in init_EVR'
        ELSE IF (base_FileName /= "") THEN
          File_path = base_FileName
        END IF

        para_mem%mem_debug = mem_debug

        EVR = .NOT. (analysis_only .OR. GridTOBasis_test .OR.            &
                     OpPsi_test .OR. cart .OR. main_test .OR. nDfit .OR. &
                     nDGrid .OR. optimization /= 0 .OR. analysis_only)

       IF (.NOT. EVR) THEN
          write(out_unitp,*) ' ERROR in init_EVR'
          write(out_unitp,*) ' This subroutine is used only to initialized ElVibRot dynamics'
          write(out_unitp,*) ' => EVR=T'
          write(out_unitp,system)
          write(out_unitp,*) ' ERROR in init_EVR'
          STOP 'ERROR in init_EVR'
       END IF
       IF (intensity_only) THEN
          write(out_unitp,*) ' ERROR in init_EVR'
          write(out_unitp,*) ' This subroutine is used only to initialized ElVibRot dynamics'
          write(out_unitp,*) ' intensity_only CANNOT be true'
          write(out_unitp,system)
          write(out_unitp,*) ' ERROR in init_EVR'
          STOP 'ERROR in init_EVR'
       END IF
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


        openmp                = Popenmp 
        openmpi               = Popenmpi

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

        IF(MPI_id==0 .AND. .NOT. openmpi) THEN
          write(out_unitp,*) '========================================='
          write(out_unitp,*) 'OpenMP parameters:'
          write(out_unitp,*) 'Max number of threads:           ',maxth
          write(out_unitp,*) 'MatOp_omp,      MatOp_maxth      ',MatOp_omp,MatOp_maxth
          write(out_unitp,*) 'OpPsi_omp,      OpPsi_maxth      ',OpPsi_omp,OpPsi_maxth
          write(out_unitp,*) 'BasisTOGrid_omp,BasisTOGrid_maxth',BasisTOGrid_omp,BasisTOGrid_maxth
          write(out_unitp,*) 'Grid_omp,       Grid_maxth       ',Grid_omp,Grid_maxth
          write(out_unitp,*) 'SG4_omp,        SG4_maxth        ',SG4_omp,SG4_maxth
          write(out_unitp,*) '========================================='

          write(out_unitp,*) '========================================='
          write(out_unitp,*) 'File_path: ',trim(adjustl(File_path))
          write(out_unitp,*) '========================================='
        ENDIF ! for MPI_id=0

        para_mem%max_mem    = max_mem/Rkind

        IF(MPI_id==0) THEN
          write(out_unitp,*) '========================================='
          write(out_unitp,*) '========================================='
          write(out_unitp,*) ' ElVibRot calculation'
          write(out_unitp,*) '========================================='
        ENDIF ! for MPI_id=0

        ! from vib.f90 file

        IF(MPI_id==0) THEN
          write(out_unitp,*)
          write(out_unitp,*) '================================================='
          write(out_unitp,*) ' VIB: BEGINNING ini_data'
          CALL time_perso('ini_data ini')
          write(out_unitp,*)
        ENDIF
        !CALL system_mem_usage(memory_RSS,'before ini_data')

        CALL   ini_data(para_EVRT%const_phys,                              &
                        para_EVRT%para_Tnum,para_EVRT%mole,                &
                        para_EVRT%para_AllBasis,                           &
                        para_EVRT%para_AllOp,para_EVRT%para_ana,           &
                        para_EVRT%para_intensity,para_EVRT%intensity_only, &
                        para_EVRT%para_propa)

        !CALL system_mem_usage(memory_RSS,'after ini_data')

        IF(MPI_id==0) THEN
          write(out_unitp,*)
          CALL time_perso('ini_data end')
          write(out_unitp,*)
          write(out_unitp,*) ' VIB: END ini_data'
          write(out_unitp,*) '================================================='
          write(out_unitp,*)
        ENDIF


        IF(MPI_id==0) THEN
          write(out_unitp,*) '========================================='
          write(out_unitp,*) '========================================='
        ENDIF ! for MPI_id=0

END SUBROUTINE init_EVR
SUBROUTINE get_nb_var(nb_var)
  USE mod_EVR
  IMPLICIT NONE

  integer,          intent(out) :: nb_var

  nb_var = para_EVRT%mole%nb_var

END SUBROUTINE get_nb_var
SUBROUTINE get_TnumRefGeom_Q0(Q0,nb_Q0,Qdyn0)
  USE mod_EVR
  IMPLICIT NONE

  integer,          intent(in)    :: nb_Q0
  real(kind=Rkind), intent(inout) :: Q0(nb_Q0)
  logical,          intent(in)    :: Qdyn0


  IF (nb_Q0 /= para_EVRT%mole%nb_var) THEN
    write(out_unitp,*) ' ERROR in get_TnumRefGeom_Q0'
    write(out_unitp,*) ' The size of Q0 is different from mole%nb_var'
    write(out_unitp,*) '   nb_Q0 ',nb_Q0
    write(out_unitp,*) '   nb_var',para_EVRT%mole%nb_var
    STOP 'ERROR in get_TnumRefGeom_Q0'
  END IF

  IF (Qdyn0) THEN
    Q0(:) = para_EVRT%mole%ActiveTransfo%Qdyn0(:)
  ELSE ! Qact order
    Q0(:) = para_EVRT%mole%ActiveTransfo%Qact0(:)
  END IF

END SUBROUTINE get_TnumRefGeom_Q0
SUBROUTINE Modify_TnumRefGeom_Q0(Q0,nb_Q0,Qdyn0)
  USE mod_EVR
  IMPLICIT NONE

  integer,          intent(in) :: nb_Q0
  real(kind=Rkind), intent(in) :: Q0(nb_Q0)
  logical,          intent(in) :: Qdyn0


  ! local variables
  real(kind=Rkind) :: Q1(nb_Q0)


  IF (nb_Q0 /= para_EVRT%mole%nb_var) THEN
    write(out_unitp,*) ' ERROR in Modify_TnumRefGeom_Q0'
    write(out_unitp,*) ' The size of Q0 is different from mole%nb_var'
    write(out_unitp,*) '   nb_Q0 ',nb_Q0
    write(out_unitp,*) '   nb_var',para_EVRT%mole%nb_var
    STOP 'ERROR in Modify_TnumRefGeom_Q0'
  END IF

  IF (Qdyn0) THEN
    para_EVRT%mole%ActiveTransfo%Qdyn0(:) =  Q0(:)
    CALL Qdyn_TO_Qact_FROM_ActiveTransfo(Q0,Q1,para_EVRT%mole%ActiveTransfo)
    para_EVRT%mole%ActiveTransfo%Qact0(:) =  Q1(:)
  ELSE
    para_EVRT%mole%ActiveTransfo%Qact0(:) =  Q0(:)
    CALL Qact_TO_Qdyn_FROM_ActiveTransfo(Q0,Q1,para_EVRT%mole%ActiveTransfo)
    para_EVRT%mole%ActiveTransfo%Qdyn0(:) =  Q1(:)
  END IF

END SUBROUTINE Modify_TnumRefGeom_Q0

SUBROUTINE get_nb_nq(nb,nq)
  USE mod_EVR
  IMPLICIT NONE

  integer,          intent(out) :: nb,nq

  nb = get_nb_FROM_basis(para_EVRT%para_AllBasis%BasisnD)
  nq = get_nq_FROM_basis(para_EVRT%para_AllBasis%BasisnD)

END SUBROUTINE get_nb_nq

SUBROUTINE levels_EVR(EigenVal,EigenVecB,EigenVecG,RhoWeight,nb,nq,nb_vec)
   USE mod_EVR

      USE mod_FullPropa
      USE mod_FullControl
      USE mod_Davidson
      USE mod_Filter
      USE mod_Arpack
      USE mod_Op
      USE mod_analysis
      USE mod_fullanalysis
      USE mod_Auto_Basis
      USE mod_psi
      USE mod_MPI_aux
      IMPLICIT NONE

      integer,           intent(in)    :: nb,nq             ! numbers of basis functions and grid points
      integer,           intent(inout) :: nb_vec            ! number of eigenvectors/eigenvalues
      real (kind=Rkind), intent(inout) :: EigenVal(nb)      ! eigenvalues.               EigenVal(i)
      real (kind=Rkind), intent(inout) :: EigenVecB(nb,nb)  ! eigenvectors on the basis. EigenVecB(:,i)
      real (kind=Rkind), intent(inout) :: EigenVecG(nq,nb)  ! eigenvectors on the grid.  EigenVecG(:,i)
      real (kind=Rkind), intent(inout) :: RhoWeight(nq)     ! rho(Q).Weight(Q), on the grid points.

!----- variables for the construction of H ---------------------------------------------
      TYPE (param_Op),       pointer :: para_H       => null()
      TYPE (param_Op),       pointer :: para_S       => null()
      TYPE (param_Op),       pointer :: para_Dip(:)  => null()
      integer                        :: iOp
      real (kind=Rkind)              :: max_Sii,max_Sij

      TYPE (CoordType),      pointer :: mole          => null() ! true POINTER
      TYPE (Tnum),           pointer :: para_Tnum     => null() ! true POINTER
      TYPE (param_AllBasis), pointer :: para_AllBasis => null() ! true POINTER

!----- for Davidson diagonalization ----------------------------------------------------
      integer                     :: nb_diago
      integer                     :: max_diago
      TYPE (param_psi),pointer    :: Tab_Psi(:) => null()
      real (kind=Rkind),allocatable  :: Ene0(:)

      integer :: Get_nbPERsym_FROM_SymAbelianOFAllBasis ! function

!----- variables divers ----------------------------------------------------------------
      integer           :: i,ip,i_baie,f_baie,id,nb_ScalOp
      real (kind=Rkind) :: T,DE,Ep,Em,Q,fac,zpe,pop
      logical           :: print_mat
      integer           :: err
      integer           :: err_mem,memory
      real (kind=Rkind) :: part_func ! function
      integer           :: nio_res_int


!para_mem%mem_debug=.TRUE.
!---------------------------------------------------------------------------------------

      para_H        => para_EVRT%para_AllOp%tab_Op(1)
      mole          => para_H%mole
      para_Tnum     => para_H%para_Tnum
      para_AllBasis => para_H%para_AllBasis

!---------------------------------------------------------------------------------------
!      Grids (V, T, Dip) calculations
!---------------------------------------------------------------------------------------
      !> turn off the allocation of grid, to be done in action
      !> sub_qa_bhe should be refined later
!#if(run_MPI)
!      Grid_allco=.FALSE.
!#endif
      IF(MPI_id==0 .AND. print_level > -1) THEN
        write(out_unitp,*)
        write(out_unitp,*) '================================================='
        write(out_unitp,*) ' VIB: BEGINNING sub_qa_bhe'
        CALL time_perso('sub_qa_bhe ini')
        write(out_unitp,*)
      ENDIF

      IF (para_H%para_ReadOp%para_FileGrid%Read_FileGrid) THEN ! test only for H
        IF(MPI_id==0 .AND. print_level > -1) THEN
          write(out_unitp,*)
          write(out_unitp,*) ' The grid of operators (S V Veff T1 and T2) will be read'
          write(out_unitp,*)
        ENDIF
        IF (para_H%para_ReadOp%para_FileGrid%Test_Grid) THEN
          IF(MPI_id==0) write(out_unitp,*) 'Test_Grid=.TRUE. => STOP in vib'
          STOP
        END IF
      END IF

      IF (para_EVRT%para_ana%VibRot) THEN
         para_Tnum%JJ = para_EVRT%para_ana%JJmax
         para_Tnum%With_Cart_Transfo = (para_Tnum%JJ>0) .AND. mole%Cart_transfo
      END IF

      ! calculate potential in action
      CALL sub_qa_bhe(para_EVRT%para_AllOp)

      IF (para_EVRT%para_ana%VibRot) THEN
         para_Tnum%JJ = 0
         para_Tnum%With_Cart_Transfo = (para_Tnum%JJ>0) .AND. mole%Cart_transfo
      END IF

      IF(MPI_id==0 .AND. print_level > -1) THEN
        write(out_unitp,*)
        CALL time_perso('sub_qa_bhe end')
        write(out_unitp,*) ' VIB: END sub_qa_bhe'
        write(out_unitp,*) '================================================='
        write(out_unitp,*)
      ENDIF

!#if(run_MPI)
!      Grid_allco=.True.
!#endif
!---------------------------------------------------------------------------------------
!      contraction of the active basis set with HADA basis
!---------------------------------------------------------------------------------------
      IF (para_H%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC .AND. para_H%nb_bi>1) THEN
        IF(MPI_id==0 .AND. print_level > -1) THEN
          write(out_unitp,*)
          write(out_unitp,*) '================================================='
          write(out_unitp,*) ' VIB: BEGINNING HADA contraction',          &
                                   para_AllBasis%BasisnD%nb
          CALL time_perso('HADA contraction')
          write(out_unitp,*)
          write(out_unitp,*)
        ENDIF

        CALL sub_MatOp_HADA(para_H,para_EVRT%para_ana,para_EVRT%para_intensity,para_EVRT%para_AllOp,  &
                            para_EVRT%const_phys)

        IF(MPI_id==0 .AND. print_level > -1) THEN
          write(out_unitp,*)
          write(out_unitp,*)
          CALL time_perso('HADA contraction')
          write(out_unitp,*) ' VIB: END HADA contraction'
          write(out_unitp,*) '================================================='
          write(out_unitp,*)
        ENDIF

        para_EVRT%para_AllOp%tab_Op(:)%nb_tot     =                     &
                sum(para_H%para_AllBasis%basis_ext2n%nb_ba_ON_HAC(:)) * &
                    para_H%para_ReadOp%nb_elec
        para_EVRT%para_AllOp%tab_Op(:)%nb_tot_ini = para_EVRT%para_AllOp%tab_Op(:)%nb_tot
      END IF

!=====================================================================
!       => Time-independent calculation
!=====================================================================

        !================================================================
        !================================================================
        !================================================================
        !===== Tune the number of threads (for SG4) =====================
        !================================================================
        max_diago = max(10,para_EVRT%para_propa%para_Davidson%nb_WP,para_H%nb_tot/10)
        max_diago = min(max_diago,10,para_H%nb_tot)
        !CALL Tune_SG4threads_HPsi(para_H%cplx,max_diago,para_H)

        !================================================================
        !===== build S and/or H if necessary ============================
        !================================================================
        IF (para_H%Make_Mat) THEN
          IF(MPI_id==0 .AND. print_level > -1) THEN
            write(out_unitp,*)
            write(out_unitp,*) '================================================='
            write(out_unitp,*) ' VIB: BEGINNING sub_matOp',para_H%nb_tot
            CALL time_perso('sub_matOp: H and S')
            write(out_unitp,*)
            write(out_unitp,*) 'para_S...%comput_S',para_EVRT%para_AllOp%tab_Op(2)%para_ReadOp%comput_S
            write(out_unitp,*)
          ENDIF

          IF (para_EVRT%para_AllOp%tab_Op(2)%para_ReadOp%comput_S) THEN

            para_S => para_EVRT%para_AllOp%tab_Op(2)
            CALL sub_MatOp(para_S,para_EVRT%para_ana%print)

            !- analysis of the overlap matrix
            CALL sub_ana_S(para_S%Rmat,para_S%nb_tot,max_Sii,max_Sij,.TRUE.)

            CALL dealloc_para_Op(para_S)
            nullify(para_S)
          END IF ! for para_EVRT%para_AllOp%tab_Op(2)%para_ReadOp%comput_S

          CALL sub_MatOp(para_H,para_EVRT%para_ana%print)


          IF(MPI_id==0 .AND. print_level > -1) THEN
            write(out_unitp,*)
            write(out_unitp,*)
            CALL time_perso('sub_matOp: H and S')
            write(out_unitp,*) ' VIB: END sub_matOp'
            write(out_unitp,*) '================================================='
            write(out_unitp,*)
          ENDIF
        END IF ! for para_H%Make_Mat
        !================================================================
        !================================================================

        !================================================================
        !===== Hmax calculation (for filter diagonalization)
        !================================================================
        IF (para_EVRT%para_ana%filter .AND. para_EVRT%para_propa%auto_Hmax) THEN
          IF(MPI_id==0 .AND. print_level > -1) THEN
            write(out_unitp,*)
            write(out_unitp,*) '================================================'
            write(out_unitp,*) ' VIB: Hmax and Hmin calculation'
            CALL time_perso('sub_Hmax ini2')
            write(out_unitp,*)
            write(out_unitp,*)
          ENDIF

          IF(.NOT. openmpi) CALL sub_Hmax(para_EVRT%para_propa,para_H)

          IF(MPI_id==0 .AND. print_level > -1) THEN
            write(out_unitp,*)
            write(out_unitp,*)
            CALL time_perso('sub_Hmax end2')
            write(out_unitp,*) ' VIB: END Hmax and Hmin calculation'
            write(out_unitp,*) '================================================'
            write(out_unitp,*)
          ENDIF
        END IF ! for para_EVRT%para_ana%filter .AND. para_EVRT%para_propa%auto_Hmax

        !================================================================
        !===== Diagonalisation ==========================================
        !================================================================
        IF (para_EVRT%para_ana%davidson .OR. para_EVRT%para_ana%arpack .OR. para_EVRT%para_ana%filter) THEN
          IF(MPI_id==0 .AND. print_level > -1) THEN
            write(out_unitp,*)
            write(out_unitp,*) '================================================='
            write(out_unitp,*) ' VIB: BEGINNING ITERATIVE DIAGONALIZATION'
            CALL time_perso('sub_Iterative_Diago')
            write(out_unitp,*)
          ENDIF

          IF (para_EVRT%para_propa%para_Davidson%max_WP == 0) THEN
            max_diago = max(1000,para_EVRT%para_propa%para_Davidson%nb_WP,          &
                            para_H%nb_tot/10)
          ELSE
            max_diago = para_EVRT%para_propa%para_Davidson%max_WP
          END IF
          IF (Get_nbPERsym_FROM_SymAbelianOFAllBasis(para_AllBasis,       &
                               para_EVRT%para_propa%para_Davidson%symab) == 0) THEN
            max_diago = min(max_diago,para_H%nb_tot)
          ELSE
            max_diago = min(max_diago,para_H%nb_tot,                      &
                  Get_nbPERsym_FROM_SymAbelianOFAllBasis(para_AllBasis, &
                                        para_EVRT%para_propa%para_Davidson%symab))
          END IF
          para_EVRT%para_propa%para_Davidson%max_WP = max_diago

          nb_diago = min(para_EVRT%para_propa%para_Davidson%nb_WP,para_H%nb_tot,max_diago)
!#if(run_MPI)
!          CALL MPI_Bcast(nb_diago,size1_MPI,Int_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
!#endif
          IF(openmpi .AND. MPI_scheme/=1) CALL MPI_Bcast_(nb_diago,size1_MPI,root_MPI)

          nullify(Tab_Psi)
          CALL alloc_array(Tab_Psi,(/ max_diago /),"Tab_Psi","vib")
          CALL alloc_NParray(Ene0,(/max_diago/),"Ene0","vib")

          IF (para_EVRT%para_ana%davidson) THEN

            CALL sub_propagation_Davidson(Tab_Psi,Ene0,nb_diago,max_diago,             &
                                          para_H,para_EVRT%para_propa%para_Davidson,   &
                                          para_EVRT%para_propa)

          ELSE IF (para_EVRT%para_ana%arpack) THEN ! arpack=t
            !CALL sub_propagation_Arpack(Tab_Psi,Ene0,nb_diago,max_diago,  &
            !                            para_H,para_EVRT%para_propa)
            CALL sub_propagation_Arpack_Sym(Tab_Psi,Ene0,nb_diago,max_diago,&
                                           para_H,para_EVRT%para_propa)

          ELSE ! filter diagonalization
            CALL sub_GaussianFilterDiagonalization(Tab_Psi,Ene0,nb_diago,max_diago,&
                                           para_H,para_EVRT%para_propa)

            !CALL sub_GaussianFilterDiagonalization_v0(Tab_Psi,Ene0,nb_diago,max_diago,&
            !                               para_H,para_EVRT%para_propa)

            !CALL sub_FilterDiagonalization(Tab_Psi,Ene0,nb_diago,max_diago,&
            !                               para_H,para_EVRT%para_propa)

          END IF

          CALL dealloc_NParray(Ene0,"Ene0","vib")
          para_H%spectral = .TRUE.

          IF(MPI_id==0 .AND. print_level > -1) THEN
            write(out_unitp,*)
            write(out_unitp,*)
            CALL time_perso('sub_Iterative_Diago')
            write(out_unitp,*) ' VIB: END ITERATIVE DIAGONALIZATION'
            write(out_unitp,*) '================================================='
            write(out_unitp,*)
          ENDIF

        ELSE ! for para_EVRT%para_ana%davidson .OR. para_EVRT%para_ana%arpack .OR. para_EVRT%para_ana%filter
          IF(MPI_id==0 .AND. print_level > -1) THEN
            write(out_unitp,*)
            write(out_unitp,*) '================================================='
            write(out_unitp,*) ' VIB: BEGINNING DIAGONALIZATION'
            CALL time_perso('sub_diago_H')
            write(out_unitp,*)
            write(out_unitp,*)
          ENDIF

          nb_diago  = para_H%nb_tot
          max_diago = para_H%nb_tot

          nullify(Tab_Psi)
          IF(keep_MPI) CALL alloc_array(Tab_Psi,(/ max_diago /),"Tab_Psi","vib")

          IF(keep_MPI) THEN
            DO i=1,max_diago
              CALL init_psi(Tab_psi(i),para_H,para_H%cplx)
            END DO
          ENDIF

          para_H%diago = .TRUE.
          IF(keep_MPI) CALL alloc_para_Op(para_H,Grid=.FALSE.,Mat=.TRUE.)
          IF (para_H%cplx) THEN
            IF(keep_MPI) THEN
              CALL sub_diago_CH(para_H%Cmat,para_H%Cdiag,para_H%Cvp,      &
                                para_H%nb_tot)

              nb_diago = count(abs(para_H%Cdiag(:)-para_H%Cdiag(1))< para_EVRT%para_ana%max_ene)
              nb_diago = min(nb_diago,para_EVRT%para_ana%max_ana)
              IF (nb_diago < 3) nb_diago = 10
              IF (nb_diago > para_H%nb_tot) nb_diago = para_H%nb_tot
              DO i=1,nb_diago
                CALL alloc_psi(Tab_Psi(i))
                Tab_Psi(i)%CvecB(:)    = para_H%Cvp(:,i)
                Tab_psi(i)%CAvOp       = para_H%Cdiag(i)
                Tab_psi(i)%IndAvOp     = para_H%n_Op  ! it should be 0
                Tab_psi(i)%convAvOp    = .TRUE.
              END DO
              para_H%para_AllBasis%basis_ext%Cvp_spec    => para_H%Cvp
            ENDIF
          ELSE
            IF(keep_MPI) THEN
              CALL sub_diago_H(para_H%Rmat,para_H%Rdiag,para_H%Rvp,       &
                               para_H%nb_tot,para_H%sym_Hamil)

              IF(print_level > -1) THEN
                write(out_unitp,*) 'HMin,HMax (ua)  : ',(/ minval(para_H%Rdiag),maxval(para_H%Rdiag) /)
                write(out_unitp,*) 'HMin,HMax (cm-1): ',   &
                  (/ minval(para_H%Rdiag),maxval(para_H%Rdiag) /)*get_Conv_au_TO_unit('E','cm-1')
              END IF

              nb_diago = count((para_H%Rdiag(:)-para_H%Rdiag(1))< para_EVRT%para_ana%max_ene)
              IF (para_EVRT%para_ana%max_ana > 0) nb_diago = min(nb_diago,para_EVRT%para_ana%max_ana)
              IF (nb_diago < 3) nb_diago = 10
              IF (nb_diago > para_H%nb_tot) nb_diago = para_H%nb_tot
              DO i=1,nb_diago
                CALL alloc_psi(Tab_Psi(i))
                Tab_Psi(i)%RvecB(:)    = para_H%Rvp(:,i)
                Tab_psi(i)%CAvOp       = para_H%Rdiag(i)
                Tab_psi(i)%IndAvOp     = para_H%n_Op  ! it should be 0
                Tab_psi(i)%convAvOp    = .TRUE.
                CALL Set_symab_OF_psiBasisRep(Tab_psi(i))
              END DO
              para_H%para_AllBasis%basis_ext%Rvp_spec    => para_H%Rvp
            ENDIF ! for keep_MPI
          END IF ! for para_H%cplx
          para_H%para_AllBasis%basis_ext%nb_vp_spec  = nb_diago

          IF (associated(para_AllBasis%basis_ext%liste_spec))           &
               CALL dealloc_NRarray(para_AllBasis%basis_ext%liste_spec, &
                                                 "basis_ext%liste_spec","vib")
          CALL alloc_NRarray(para_AllBasis%basis_ext%liste_spec,        &
                               [nb_diago],"basis_ext%liste_spec","vib")
          para_AllBasis%basis_ext%liste_spec(:) = (/ (i,i=1,nb_diago) /)


          IF (associated(para_H%Rmat)) THEN
            CALL dealloc_array(para_H%Rmat,"para_H%Rmat","vib")
          END IF
          IF (associated(para_H%Cmat)) THEN
            CALL dealloc_array(para_H%Cmat,"para_H%Cmat","vib")
            nullify(para_H%Cmat)
          END IF
          para_EVRT%para_ana%max_ana = nb_diago

          IF(MPI_id==0 .AND. print_level > -1) THEN
            write(out_unitp,*)
            write(out_unitp,*)
            CALL time_perso('sub_diago_H')
            write(out_unitp,*) ' VIB: END DIAGONALIZATION'
            write(out_unitp,*) '================================================='
            write(out_unitp,*)
          ENDIF
        END IF ! for para_EVRT%para_ana%davidson .OR. para_EVRT%para_ana%arpack .OR. para_EVRT%para_ana%filter
        CALL flush_perso(out_unitp)
        !===============================================================
        !===============================================================

         IF (nb /= get_nb_FROM_basis(para_AllBasis%BasisnD) .OR.        &
             nq /= get_nq_FROM_basis(para_AllBasis%BasisnD)) THEN

           write(out_unitp,*) ' ERROR in levels_EVR'
           write(out_unitp,*) ' inconsistant nb values',nb,get_nb_FROM_basis(para_AllBasis%BasisnD)
           write(out_unitp,*) '   or'
           write(out_unitp,*) ' inconsistant nq values',nq,get_nq_FROM_basis(para_AllBasis%BasisnD)
           STOP 'ERROR in levels_EVR'

         END IF

         nb_vec = nb_diago

         DO i=1,nb_diago
           EigenVal(i)    = para_H%Rdiag(i)
           EigenVecB(:,i) = Tab_Psi(i)%RvecB(:)

           CALL sub_PsiBasisRep_TO_GridRep(Tab_Psi(i))
           EigenVecG(:,i) = Tab_Psi(i)%RvecG(:)
           CALL alloc_psi(Tab_Psi(i),BasisRep=.TRUE.,GridRep=.FALSE.)

         END DO

         !- calculation of WrhonD ------------------------------
         DO i=1,nq
           RhoWeight(i) = Rec_WrhonD(para_AllBasis%BasisnD,i)
         END DO

        !===============================================================
        !===============================================================
        IF(MPI_id==0 .AND. print_level > -1) THEN
          write(out_unitp,*)
          write(out_unitp,*) '================================================='
          write(out_unitp,*) ' VIB: BEGINNING WAVE FUNCTION ANALYSIS'
          CALL time_perso('sub_analyse ini')
          write(out_unitp,*)
        ENDIF

        IF(keep_MPI) CALL sub_analyse(Tab_Psi,nb_diago,para_H,         &
                           para_EVRT%para_ana,para_EVRT%para_intensity, &
                               para_EVRT%para_AllOp,para_EVRT%const_phys)

        CALL flush_perso(out_unitp)

        IF (.NOT. para_H%cplx .AND. para_EVRT%para_ana%VibRot) THEN
          CALL sub_VibRot(Tab_Psi,para_EVRT%para_ana%max_ana,para_H,para_EVRT%para_ana)
        END IF
        CALL flush_perso(out_unitp)

        !===============================================================
        ! Spectral representation of operator
        !===============================================================
        print_mat = (para_EVRT%para_ana%MaxWP_TO_Write_MatOp >= nb_diago)
        IF (para_EVRT%para_ana%Spectral_ScalOp) THEN
          DO iOp=1,size(para_EVRT%para_AllOp%tab_Op)

            IF (para_EVRT%para_AllOp%tab_Op(iOp)%n_Op == -1) CYCLE ! S
            IF (para_EVRT%para_AllOp%tab_Op(iOp)%spectral) THEN
              write(out_unitp,*) '==========================================='
              write(out_unitp,*) ' Spectral representation of: ',        &
                                 trim(para_EVRT%para_AllOp%tab_Op(iOp)%name_Op)

              CALL sub_build_MatOp(Tab_Psi,nb_diago,                     &
                                   para_EVRT%para_AllOp%tab_Op(iOp),.TRUE.,print_mat)

              write(out_unitp,*) '==========================================='
              CALL flush_perso(out_unitp)
            END IF
          END DO
        END IF

        IF(MPI_id==0 .AND. print_level > -1) THEN
          write(out_unitp,*)
          write(out_unitp,*)
          CALL time_perso('sub_analyse end')
          write(out_unitp,*) ' VIB: END WAVE FUNCTION ANALYSIS'
          write(out_unitp,*) '================================================='
          write(out_unitp,*)
        ENDIF

        !===============================================================
        !===============================================================

        !===============================================================
        ! Deallocation of Tab_Psi
        !===============================================================
        IF(keep_MPI) THEN
          DO i=1,size(Tab_Psi)
            CALL dealloc_psi(Tab_Psi(i))
          END DO
          CALL dealloc_array(Tab_Psi,"Tab_Psi","vib")
        ENDIF
        !===============================================================
        !===============================================================

!=====================================================================
!       => end Time-independent calculation
!=====================================================================
!=====================================================================
!=====================================================================
!RETURN
write(out_unitp,*) 'intensity',para_EVRT%para_ana%intensity ; flush(out_unitp)

      IF (.NOT. para_H%cplx .AND. para_EVRT%para_ana%intensity) THEN
        IF(MPI_id==0 .AND. print_level > -1) THEN
          write(out_unitp,*)
          write(out_unitp,*) '================================================'
          write(out_unitp,*) ' VIB: BEGINNING sub_intensity',para_H%nb_tot,para_EVRT%para_ana%max_ana
          CALL time_perso('sub_intensity')
          write(out_unitp,*)
          write(out_unitp,*)
        ENDIF
        !----- file for restart or for changed the temprature --------------
        CALL file_open(para_EVRT%para_intensity%file_resart_int,nio_res_int)

        IF (para_EVRT_calc%intensity_only) THEN
          read(nio_res_int,*,IOSTAT=err) para_H%nb_tot,para_EVRT%para_ana%max_ana
          read(nio_res_int,*,IOSTAT=err)

          para_H%diago = .TRUE.
          CALL alloc_para_Op(para_H)

          CALL Read_Vec(para_H%Rdiag,nio_res_int,5,err)
          IF (err /= 0) THEN
            write(out_unitp,*) 'ERROR in vib'
            write(out_unitp,*) ' reading the vector "para_H%Rdiag"'
            STOP
          END IF

          IF (para_EVRT_calc%intensity_only) THEN
            write(out_unitp,*)
            Q =  part_func(para_H%Rdiag,size(para_H%Rdiag),para_EVRT%para_ana%Temp,para_EVRT%const_phys)
            fac = para_EVRT%const_phys%Eh / (para_EVRT%const_phys%k * para_EVRT%para_ana%Temp)
            zpe = minval(para_H%Rdiag(:))
            write(out_unitp,*) 'population at T, Q',para_EVRT%para_ana%Temp,Q
            write(out_unitp,*) 'Energy level (',para_EVRT%const_phys%ene_unit,') pop:'
            DO i=1,size(para_H%Rdiag)
              pop = exp(-(para_H%Rdiag(i)-zpe)*fac)
              write(out_unitp,10) i,para_H%Rdiag(i) * para_EVRT%const_phys%auTOenergy,&
                    (para_H%Rdiag(i) - zpe) * para_EVRT%const_phys%auTOenergy,pop/Q
10             format(i4,x,3f20.5)
            END DO
          END IF

          read(nio_res_int,*,IOSTAT=err)

          CALL Read_Mat(para_H%Rvp,nio_res_int,5,err)
          IF (err /= 0) THEN
            write(out_unitp,*) 'ERROR in vib'
            write(out_unitp,*) ' reading the matrix "para_H%Rvp"'
            STOP
          END IF
        ELSE ! for intensity_only
          write(out_unitp,*) 'write restart file for intensity: ',      &
                           para_EVRT%para_intensity%file_resart_int%name
          CALL flush_perso(out_unitp)
          write(nio_res_int,*) para_H%nb_tot,para_EVRT%para_ana%max_ana
          write(nio_res_int,*) 'ene'
          CALL Write_Vec(para_H%Rdiag,nio_res_int,5,Rformat='e30.23')
          write(nio_res_int,*) 'psi'
          CALL flush_perso(out_unitp)
          CALL Write_Mat(para_H%Rvp,nio_res_int,5,Rformat='e30.23')
          CALL flush_perso(nio_res_int)
        END IF ! for intensity_only
        CALL flush_perso(out_unitp)
        !-------------------------------------------------------------------

        iOp = 2
        para_Dip => para_EVRT%para_AllOp%tab_Op(iOp+1:iOp+3)
        CALL sub_intensity(para_Dip,                                  &
                           para_EVRT%para_ana%print,para_H,para_EVRT%para_ana%max_ana,    &
                           para_EVRT%para_intensity,para_EVRT%const_phys,para_EVRT_calc%intensity_only,nio_res_int)

        IF(MPI_id==0 .AND. print_level > -1) THEN
          write(out_unitp,*)
          write(out_unitp,*)
          CALL time_perso('sub_intensity')
          write(out_unitp,*) ' VIB: END sub_intensity'
          write(out_unitp,*) '================================================'
          write(out_unitp,*)
        ENDIF
        close(nio_res_int)
        nullify(para_Dip)
      END IF !for .NOT. para_H%cplx .AND. para_EVRT%para_ana%intensity
      CALL flush_perso(out_unitp)
RETURN
      IF (.NOT. para_H%cplx .AND. para_EVRT%para_ana%Psi_ScalOp) THEN
        IF(MPI_id==0 .AND. print_level > -1) THEN
          write(out_unitp,*)
          write(out_unitp,*) '================================================'
          write(out_unitp,*) ' VIB: BEGINNING sub_AnalysePsy_ScalOp',para_H%nb_tot,para_EVRT%para_ana%max_ana
          CALL time_perso('sub_AnalysePsy_ScalOp')
          write(out_unitp,*)
          write(out_unitp,*) 'nb_scalar_Op',para_H%para_ReadOp%nb_scalar_Op
        ENDIF

        iOp = 2
        nb_ScalOp = para_H%para_ReadOp%nb_scalar_Op
        para_Dip => para_EVRT%para_AllOp%tab_Op(iOp+1:iOp+nb_ScalOp)
        CALL sub_AnalysePsy_ScalOp(para_Dip,nb_ScalOp,para_H,para_EVRT%para_ana%max_ana)

        IF(MPI_id==0 .AND. print_level > -1) THEN
          write(out_unitp,*)
          write(out_unitp,*)
          CALL time_perso('sub_AnalysePsy_ScalOp')
          write(out_unitp,*) ' VIB: END sub_AnalysePsy_ScalOp'
          write(out_unitp,*) '================================================'
          write(out_unitp,*)
        ENDIF
        nullify(para_Dip)
      END IF
      CALL flush_perso(out_unitp)

      IF (.NOT. para_H%cplx .AND. para_EVRT%para_ana%NLO) THEN
        IF(MPI_id==0 .AND. print_level > -1) THEN
          write(out_unitp,*)
          write(out_unitp,*) '================================================'
          write(out_unitp,*) ' VIB: BEGINNING sub_NLO',para_H%nb_tot,para_EVRT%para_ana%max_ana
          CALL time_perso('sub_NLO')
          write(out_unitp,*)
          write(out_unitp,*)
        ENDIF

        iOp = 2
        para_Dip => para_EVRT%para_AllOp%tab_Op(iOp+1:iOp+3)
        CALL sub_NLO(para_Dip,para_EVRT%para_ana%print,para_H,para_EVRT%para_ana%max_ana, &
                     para_EVRT%para_intensity)

        IF(MPI_id==0 .AND. print_level > -1) THEN
          write(out_unitp,*)
          write(out_unitp,*)
          CALL time_perso('sub_NLO')
          write(out_unitp,*) ' VIB: END sub_NLO'
          write(out_unitp,*) '================================================'
          write(out_unitp,*)
        ENDIF
        nullify(para_Dip)
      END IF
      CALL flush_perso(out_unitp)

!=====================================================================
!=====================================================================

!=====================================================================
!       deallocated memories
!=====================================================================
      IF ( associated(Tab_Psi) ) THEN
        DO i=1,size(Tab_Psi)
          CALL dealloc_psi(Tab_Psi(i))
        END DO
        CALL dealloc_array(Tab_Psi,"Tab_Psi","vib")
      END IF

      IF (associated(para_EVRT%para_AllOp%tab_Op)) THEN
        DO i=1,size(para_EVRT%para_AllOp%tab_Op)
          CALL dealloc_para_Op(para_EVRT%para_AllOp%tab_Op(i),keep_init=.TRUE.)
        END DO
      END IF

      IF(print_level > -1) THEN
        write(out_unitp,*) 'mem_tot,max_mem_used',para_mem%mem_tot,para_mem%max_mem_used
        write(out_unitp,*) 'nb_alloc,nb_dealloc',para_mem%nb_alloc,para_mem%nb_dealloc
        write(out_unitp,*) '================================================'
        IF(openmpi) THEN
          write(out_unitp,*) ' ElVibRot-Tnum AU REVOIR!!!', ' from ', MPI_id
        ELSE
          write(out_unitp,*) ' ElVibRot-Tnum AU REVOIR!!!'
        ENDIF
        write(out_unitp,*) '================================================'
      END IF

   END SUBROUTINE levels_EVR


SUBROUTINE Finalize_EVR()
  USE mod_EVR
  IMPLICIT NONE

      CALL dealloc_table_at(para_EVRT%const_phys%mendeleev)

      CALL dealloc_CoordType(para_EVRT%mole)
      IF (associated(para_EVRT%para_Tnum%Gref)) THEN
        CALL dealloc_array(para_EVRT%para_Tnum%Gref,"para_Tnum%Gref","vib")
      END IF
      !CALL dealloc_Tnum(para_EVRT%para_Tnum)

      CALL dealloc_para_AllOp(para_EVRT%para_AllOp)
      CALL dealloc_para_ana(para_EVRT%para_ana)
      CALL dealloc_param_propa(para_EVRT%para_propa)
      CALL dealloc_psi(para_EVRT%WP0)
      CALL dealloc_AllBasis(para_EVRT%para_AllBasis)

      IF(openmpi) THEN
        CALL end_MPI()
        CALL time_perso('MPI closed')
      ENDIF
      close(in_unitp)

END SUBROUTINE Finalize_EVR

