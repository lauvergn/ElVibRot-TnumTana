DIR_EVRT:=$(shell pwd)
#===============================================================================
#===============================================================================
## Compiler? Possible values: ifort; gfortran; pgf90 (v17),mpifort
# F90 = mpifort
  F90 = gfortran
# F90 = nagfor
# F90 = ifort
# F90 = pgf90

## parallel_make=1 to enable parallel make
## parallel_make=0 for fast debug make, no parallel
parallel_make=0

## Optimize? Empty: default No optimization; 0: No Optimization; 1 Optimzation
OPT = 0
#
## OpenMP? Empty: default with OpenMP; 0: No OpenMP; 1 with OpenMP
OMP = 1
#
## force the default integer (without kind) during the compillation.
## default 4: , INT=8 (for kind=8)
INT = 4
#
## Arpack? Empty: default No Arpack; 0: without Arpack; 1 with Arpack
ARPACK = 0
## CERFACS? Empty: default No CERFACS; 0: without CERFACS; 1 with CERFACS
CERFACS = 0
## Lapack/blas/mkl? Empty: default with Lapack; 0: without Lapack; 1 with Lapack
LAPACK = 1
## Quantum Model Lib (QMLib) Empty: default with QMLib; 0: without QMLib; 1 with QMLib
## Always QML=1
QML = 1
#
## extension for the "sub_system." file. Possible values: f; f90 or $(EXTFextern)
## if $(EXTFextern) is empty, the default is f
extf = $(EXTFextern)
## Some compilers (like PGF90) do not have inverse hyperbolic functions: atanh, asinh, acosh
# NVHYP  = 1 : with intrinsic inverse hyperbolic functions
# NVHYP  = 0 : with external inverse hyperbolic functions (without intrinsic ones)
INVHYP  = 1
#
## Operating system, OS? automatic using uname:
OS=$(shell uname)
#
#===============================================================================
# External Libraries directory (QML ...)
ExternalLibDIR=Ext_Lib
#===============================================================================
#
#===============================================================================
# turn off ARPACK when using pgf90
ifeq ($(F90),pgf90)
  ARPACK = 0
endif

# setup for mpifort
ifeq ($(F90),mpifort)
  ## MPI compiled with: gfortran or ifort
  MPICORE := $(shell ompi_info | grep 'Fort compiler:' | sed 's/Fort compiler://g' | sed 's/ //g')
  OMP = 0
  ifeq ($(INT),8)
    ARPACK = 0 ## temp here, disable ARPACK for 64-bit case
  endif
endif

# obj directory
ifeq ($(F90),mpifort)
  obj_dir = obj/obj_$(F90)_$(MPICORE)_$(INT)
else
  obj_dir = obj/obj_$(F90)_$(INT)
endif
## turn off ARPACK
#=================================================================================
# External pot for the library: libpot.a,
# with epxort variable (POTDIRextern) or with explicit name
ExternalDIR := $(POTDIRextern)
# Example pot Bowman (new lib)
#ExternalDIR := /u/lauvergn/trav/ElVibRot-Tnum/exa_work/exa_clathrate-Bowman/H2-clathrate-PES-ompNewLib/ver1
# Valiron pot
#ExternalDIR := /Users/lauvergn/Documents/Papiers/Clathrate/Papier1-2016/Papier1-2016/Pot_Valiron/code_PES_V08
#ExternalDIR := /Users/lauvergn/git/PotV08_6D/EXEC_6D-v3
DIRLIB := -L$(ExternalDIR)
PESLIB := -lpot
# If ExternalDIR is empty, PESLIB must be empty
ifeq  ($(strip $(ExternalDIR)),)
  PESLIB =
  DIRLIB =
endif
#===============================================================================
#
#===============================================================================
# Quantum Model Lib (ECAM)
#QMLibDIR := /Users/lauvergn/git/QuantumModelLib
QMLibDIR := $(DIR_EVRT)/$(ExternalLibDIR)/QuantumModelLib
QMLModDIR := $(QMLibDIR)/OBJ/obj_$(F90)_omp$(OMP)
QMLIB := -L$(QMLibDIR) -lQMLib_$(F90)_omp$(OMP)
QMLibDIR_full := $(QMLibDIR)/libQMLib_$(F90)_omp$(OMP).a

# dnSVM Lib
#dnSVMLibDIR := /Users/lauvergn/git/AD_dnSVM
dnSVMLibDIR := $(DIR_EVRT)/$(ExternalLibDIR)/dnSVMLib
dnSVMLIB := -L$(dnSVMLibDIR) -lAD_dnSVM
dnSVMLibDIR_full := $(dnSVMLibDIR)/libAD_dnSVM.a
AD_dnSVMModDIR   := $(dnSVMLibDIR)/OBJ/obj_$(F90)_omp$(OMP)
#===============================================================================
#
#===============================================================================
# If EXTFextern is empty, extf must be empty
ifeq  ($(strip $(EXTFextern)),)
  extf = f
endif
#===============================================================================

#===============================================================================
# We cannot use ARPACK without lapack
ifeq ($(LAPACK),0)
  ARPACK = 0
endif
#===============================================================================
#
CompC=gcc

#===============================================================================
# nag compillation (nagfor)
#===============================================================================
ifeq ($(F90),nagfor)
   # for c++ preprocessing
   ifeq ($(OMP),1)
     CPPpre  = -fpp -Drun_openMP=1
   else
     CPPpre  = -fpp
   endif

   # omp management
   ifeq ($(OMP),0)
      OMPFLAG =
   else
      OMPFLAG = -openmp
   endif
   # opt management
   ifeq ($(OPT),1)
      F90FLAGS = -O4  $(OMPFLAG) -o -compatible -kind=byte -Ounroll=4 -s
   else
      #F90FLAGS = -O0 $(OMPFLAG) -g -C=all -mtrace=all
      #  -C=undefined is not compatible with: (i) -framework Accelerate or lapack lib (ii) openmp
      #with -mtrace=all add information on the memmory allocation/deallocation.
      # The option -C=dangling causes troubles since it is used => -C=all cannot be used.
      # -gline is not compatible with openmp
      ifeq ($(OMP),0)
        ifeq ($(LAPACK),0)
          F90FLAGS = -O0            -g -gline -kind=byte -C -C=alias -C=intovf -C=undefined
        else
          F90FLAGS = -O0 $(OMPFLAG) -g -gline -kind=byte -C -C=alias -C=intovf
        endif
      else
          F90FLAGS = -O0 $(OMPFLAG) -g        -kind=byte -C -C=alias -C=intovf
      endif
   endif

   ifeq ($(LAPACK),1)
     F90LIB = -framework Accelerate
   else
     F90LIB =
   endif

   ifeq ($(INT),8)
     F90FLAGS := $(F90FLAGS) -i8 -Dint8=1
   endif

   F90_VER = $(shell $(F90) -V 3>&1 1>&2 2>&3 | head -1 )

endif

#=================================================================================
# ifort compillation v17 with mkl
#=================================================================================
ifeq ($(F90),ifort)
   # for c++ preprocessing
   ifeq ($(OMP),1)
     CPPpre  = -cpp -Drun_openMP=1
   else
     CPPpre  = -cpp
   endif

   # omp management
   ifeq ($(OMP),0)
      OMPFLAG =
   else
      OMPFLAG = -qopenmp
   endif
   # opt management
   ifeq ($(OPT),1)
      F90FLAGS = -O  $(OMPFLAG) -parallel -g -traceback -assume realloc_lhs
   else
      F90FLAGS = -O0 $(OMPFLAG) -debug -g -traceback -check all -ftrapuv -assume realloc_lhs
   endif

   ifeq ($(LAPACK),1)
     F90LIB = -mkl -lpthread
     #F90LIB = $(MKLROOT)/lib/libmkl_lapack95_ilp64.a $(MKLROOT)/lib/libmkl_core.a $(MKLROOT)/lib/libmkl_blas95_ilp64.a -lpthread
   else
     F90LIB = -lpthread
   endif

   ifeq ($(INT),8)
     F90FLAGS := $(F90FLAGS) -i8 -Dint8=1
   endif
   MOD_FLAGS := -M$(QMLModDIR) -M$(AD_dnSVMModDIR)

   F90_VER = $(shell $(F90) --version | head -1 )

endif
#=================================================================================

#=================================================================================
# pgf90 compillation
#=================================================================================
ifeq ($(F90),pgf90)
   # for c++ preprocessing
   ifeq ($(OMP),1)
     CPPpre  = -Mpreprocess -Drun_openMP=1
   else
     CPPpre  = -Mpreprocess
   endif

   # omp management
   ifeq ($(OMP),0)
      OMPFLAG =
   else
      OMPFLAG = -mp=allcores
      F90LIB = -lpthread
   endif
   # opt management
   ifeq ($(OPT),1)
      F90FLAGS = -O $(OMPFLAG) -fast -Mallocatable=03
   else
      F90FLAGS = -O0 $(OMPFLAG)      -Mallocatable=03 -Mbounds -Mchkstk -g
   endif

   ifeq ($(LAPACK),1)
     F90LIB += -lblas -llapack
   else
     F90LIB +=
   endif

   F90_VER = $(shell $(F90) --version | head -2 | tail -1 )


endif
#===============================================================================

#===============================================================================
# gfortran (osx and linux)
#ifeq ($(F90),gfortran)
#===============================================================================
ifeq ($(F90),$(filter $(F90),gfortran gfortran-8))
   # for c++ preprocessing
   ifeq ($(OMP),1)
     CPPpre  = -cpp -Drun_openMP=1
   else
     CPPpre  = -cpp
   endif

   # omp management
   ifeq ($(OMP),0)
      OMPFLAG =
   else
      OMPFLAG = -fopenmp
   endif
   # OS management
   ifeq ($(LAPACK),1)
     ifeq ($(OS),Darwin)    # OSX
        # OSX libs (included lapack+blas)
        F90LIB = -framework Accelerate
        CompC  = gcc
     else                   # Linux
        # linux libs
        F90LIB = -llapack -lblas
        #
        # linux libs with mkl and with openmp
        #F90LIB = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread
        # linux libs with mkl and without openmp
        #F90LIB = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential
     endif
   else
    F90LIB =
   endif
   #
   # opt management
   # -finit-local-zero
   ifeq ($(OPT),1)
      F90FLAGS = -O5 -g -fbacktrace $(OMPFLAG) -funroll-loops -ftree-vectorize -falign-loops=16
      CFLAGS   = -O5 -g             $(OMPFLAG) -funroll-loops -ftree-vectorize -falign-loops=16
   else
      #F90FLAGS = -O0 -g -fbacktrace $(OMPFLAG) -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized -Wconversion -Wconversion-extra
      #F90FLAGS = -O0 -g -fbacktrace $(OMPFLAG) -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized -Wunused
       F90FLAGS = -O0 -g -fbacktrace $(OMPFLAG) -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized
       CFLAGS   = -O0 -g             $(OMPFLAG) -fwhole-file -Wuninitialized
      #F90FLAGS = -O0 -fbounds-check -Wuninitialized
   endif
   # integer kind management
   ifeq ($(INT),8)
      F90FLAGS := $(F90FLAGS) -fdefault-integer-8 -Dint8=1
   endif

   MOD_FLAGS := -I $(QMLModDIR) -I $(AD_dnSVMModDIR)


   F90_VER = $(shell $(F90) --version | head -1 )

endif
#=================================================================================

#=================================================================================
# mpifort (osx and linux)
# ifeq ($(F90),mpifort)
# for gfortran core
#=================================================================================
ifeq ($(F90),mpifort)
   # for c++ preprocessing
   CPPpre = -cpp -Drun_MPI=1

   # omp management
   ifeq ($(OMP),0)
      OMPFLAG =
   else
      $(error disable openMP when using MPI!)
   endif
   # OS management
   ifeq ($(LAPACK),1)
     ifeq ($(OS),Darwin)    # OSX
        # OSX libs (included lapack+blas)
        F90LIB = -framework Accelerate
     else                   # Linux
        # linux libs
        F90LIB = -llapack -lblas
        # linux libs with mkl and without openmp
        #F90LIB = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential
     endif
   else
    F90LIB =
   endif
   #
   # opt management
   ifeq ($(MPICORE), gfortran)
     CPPpre += -Drun_MPI_gfortran=1
     ifeq ($(OPT),1)
        F90FLAGS = -O5 -g -fbacktrace $(OMPFLAG) -funroll-loops -ftree-vectorize -falign-loops=16
     else
        #F90FLAGS = -O0 -g -fbacktrace $(OMPFLAG) -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized -Wconversion -Wconversion-extra
        #F90FLAGS = -O0 -g -fbacktrace $(OMPFLAG) -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized -Wunused
         F90FLAGS = -O0 -g -fbacktrace $(OMPFLAG) -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized
        #F90FLAGS = -O0 -fbounds-check -Wuninitialized
     endif
   else
     CPPpre += -Drun_MPI_ifort=1
     ifeq ($(OPT),1)
       F90FLAGS =  -O5 -g #-check all -fpe0 -warn -traceback -debug extended
     else
       F90FLAGS =  -O0 -g #-check all -fpe0 -warn -traceback -debug extended
     endif
   endif
   # integer kind management
   ifeq ($(INT),8)
      ifeq ($(MPICORE),gfortran)
         F90FLAGS += -fdefault-integer-8
      else
         F90FLAGS += -i8
      endif
   endif
endif
F90FLAGS := $(F90FLAGS)   $(EXTMOD)

ifeq ($(F90),mpifort)
  F90_VER = $(shell ompi_info | grep 'Open MPI:' | sed 's/ //g' )
else
  F90_VER = $(shell $(F90) --version | head -1 )
endif

GIT_Branch := $(shell git status | grep "On branch")

#===============================================================================
#===============================================================================
$(info ************************************************************************)
$(info ***********OS:               $(OS))
$(info ***********git:              $(GIT_Branch))
$(info ***********COMPILER:         $(F90))
$(info ***********OPTIMIZATION:     $(OPT))
$(info ***********COMPILER VERSION: $(F90_VER))
ifeq ($(F90),mpifort)
$(info ***********COMPILED with:    $(MPICORE))
endif
$(info ***********OpenMP:           $(OMPFLAG))
$(info ***********Arpack:           $(ARPACK))
$(info ***********CERFACS:          $(CERFACS))
$(info ***********Lapack:           $(LAPACK))
$(info ***********QMLib:            $(QMLIB))
$(info ***********DIR of QMLib:     $(QMLibDIR))
$(info ***********F90FLAGS:         $(F90FLAGS))
$(info ***********F90LIB:           $(F90LIB))
$(info ***********subsystem file:   sub_system.$(extf))
$(info ***********DIR of potlib.a:  $(ExternalDIR))
$(info ***********potLib:           $(PESLIB))
$(info ***********INVHYP:           $(INVHYP))
$(info ***********MOD_FLAGS:        $(MOD_FLAGS))
$(info ************************************************************************)

F90_FLAGS = $(F90) $(F90FLAGS) $(MOD_FLAGS)
LYNK90 = $(F90_FLAGS)

#===============================================================================
# Arpack library
#===============================================================================
ifeq ($(ARPACK),1)
  # Arpack management with the OS
  ifeq ($(OS),Darwin)    # OSX
    #ARPACKLIB=/Users/chen/Linux/Software/ARPACK/libarpack_MAC.a
    ARPACKLIB=/Users/lauvergn/trav/ARPACK/libarpack_OSX.a
  else                   # Linux
    ifeq ($(F90), mpifort)
      ifeq ($(MPICORE), gfortran)
        ARPACKLIB=/u/achen/Software/ARPACK/libarpack_Linux_gfortran.a
      else ifeq ($(MPICORE), ifort)
        ARPACKLIB=/u/achen/Software/ARPACK/libarpack_Linux_ifort.a
      endif
    else ifeq ($(F90), gfortran)
      ARPACKLIB=/u/achen/Software/ARPACK/libarpack_Linux_gfortran.a
    else ifeq ($(F90), ifort)
      ARPACKLIB=/u/achen/Software/ARPACK/libarpack_Linux_ifort.a
    endif
    #ARPACKLIB=/usr/lib64/libarpack.a
    ARPACKLIB=/userTMP/lauvergn/EVR/ARPACK_DML/libarpack_Linux.a
  endif
else
  ARPACKLIB =
endif
#=================================================================================
#=================================================================================
#LIBS := $(DIRLIB) $(QMLIB) $(PESLIB) $(ARPACKLIB) $(F90LIB)
LIBS := $(DIRLIB) $(QMLIB) $(dnSVMLIB) $(PESLIB) $(ARPACKLIB) $(F90LIB)
 LYNKFLAGS = $(LIBS)

#=================================================================================
#=================================================================================
DIR_EVRT:=$(shell pwd)
OBJ = $(DIR_EVRT)/$(obj_dir)

TNUM_ver:=$(shell awk '/Tnum/ {print $$3}' $(DIR_EVRT)/version-EVR-T)
TANA_ver:=$(shell awk '/Tana/ {print $$3}' $(DIR_EVRT)/version-EVR-T)
EVR_ver:=$(shell awk '/EVR/ {print $$3}' $(DIR_EVRT)/version-EVR-T)
$(info ***********************************************************************)
$(info ***************** TNUM_ver: $(TNUM_ver))
$(info ***************** TANA_ver: $(TANA_ver))
$(info ****************** EVR_ver: $(EVR_ver))
$(info ***********************************************************************)
$(info ***********************************************************************)
$(info ************ run UnitTests: make UT)
$(info ********** clean UnitTests: make clean_UT)
$(info ***********************************************************************)
#
#
CPPSHELL = -D__COMPILE_DATE="\"$(shell date +"%a %e %b %Y - %H:%M:%S")\"" \
           -D__COMPILE_HOST="\"$(shell hostname -s)\"" \
           -D__COMPILER="'$(F90)'" \
           -D__COMPILER_VER="'$(F90_VER)'" \
           -D__COMPILER_OPT="'$(F90FLAGS)'" \
           -D__COMPILER_LIBS="'$(F90LIB)'" \
           -D__EVRTPATH="'$(DIR_EVRT)'" \
           -D__EVR_VER="'$(EVR_ver)'" \
           -D__TNUM_VER="'$(TNUM_ver)'" \
           -D__TANA_VER="'$(TANA_ver)'"
CPPSHELL_ARPACK  = -D__ARPACK="$(ARPACK)"
CPPSHELL_CERFACS = -D__CERFACS="$(CERFACS)"
CPPSHELL_INVHYP  = -D__INVHYP="$(INVHYP)"
CPPSHELL_LAPACK  = -D__LAPACK="$(LAPACK)"
CPPSHELL_QML     = -D__QML="$(QML)"
#CPPSHELL_GIT     = -D__GIT="'master'"
CPPSHELL_GIT     = -D__GIT="'$(GIT_Branch)'"
#==========================================
# the different programs
#  vib:  make or make EVR
#  Tnum:  make tnum or make Tnum or make tnum-dist
#  work (program):  make work
#
VIBEXE  = vib.exe
VIBMAIN = EVR-T
VIBDIRm = sub_main
#Variables for documentation generation
#DOCPREFIX = doc
DOCGEN = ./scripts/docgen.pl
DOCINDEXGEN = ./scripts/doc_indexgen.pl
DOCINDEXFILE = doc/reference/index.store
DOCINDEX = --index=$(DOCINDEXFILE)
REFPATH = doc/reference
HTML = $(patsubst sub_module/%.f90, $(REFPATH)/%.html, $(wildcard sub_module/*.f90))
#
#
PhysConstEXE  = PhysConst.exe
PhysConstMAIN = PhysicalConstants_Main
#
KEOTESTEXE  = TEST_TnumTana.exe
KEOTEST     = TEST_TnumTana
#
Main_TnumTana_FDriverEXE=Main_TnumTana_FDriver.exe
Main_TnumTana_cDriverEXE=Main_TnumTana_cDriver.exe
#
TNUMEXE  = Tnum90.exe
TNUMMAIN = Tnum90
#
TNUMMCTDHEXE = Tnum90_MCTDH.exe
TNUMMCTDHMAIN = Tnum90_MCTDH
#
TNUM_MiddasCppEXE  = Tnum90_MidasCpp.exe
TNUM_MiddasCppMAIN = Tnum90_MidasCpp
#
GWPEXE = gauss.exe
GWPMAIN = Gauss_numlH
#
WORKEXE  = work.exe
#WORKMAIN = Tnum90_AverageHessian
WORKMAIN = CurviRPH

#==========================================
EXE = $(VIBEXE) $(TNUMEXE) $(TNUMDISTEXE) $(GWPEXE) $(WORKEXE)
#==========================================
#
#Lib and its directories
DirLib     = $(DIR_EVRT)/Source_Lib

DirSys     = $(DirLib)/sub_system
DirMath    = $(DirLib)/sub_communf90/sub_math
DirIO      = $(DirLib)/sub_communf90/sub_io
DirSHTOOLS = $(DirLib)/sub_communf90/SHTOOLS
DirdnSVM   = $(DirLib)/sub_dnSVM
DirnDind   = $(DirLib)/sub_nDindex

DirMod     = $(DirLib)/sub_module

#Physical constants, masses
DIRPhyCte = $(DIR_EVRT)/Source_PhysicalConstants

#Tnum + Tana + Coordinates
DirTNUM   = $(DIR_EVRT)/Source_TnumTana_Coord

#Primitive Operators
DIRPrimOp = $(DIR_EVRT)/Source_PrimOperator
DirPot    = $(DIR_EVRT)/sub_pot

#ElVibRot + Optimization + fit
DirEVR    = $(DIR_EVRT)/Source_ElVibRot

DIRvib     = $(DirEVR)/sub_main
DIR1       = $(DirEVR)/sub_data_initialisation
DIRba      = $(DirEVR)/sub_Basis
DIRbaSG4   = $(DirEVR)/sub_Basis/sub_Basis_SG4
DIRWP      = $(DirEVR)/sub_WP
DIR2       = $(DirEVR)/sub_inactive
DIR5       = $(DirEVR)/sub_active
DIROp      = $(DirEVR)/sub_Operator
DIRana     = $(DirEVR)/sub_analysis
DIRpropa   = $(DirEVR)/sub_propagation
DIRSmolyak = $(DirEVR)/sub_Smolyak_test
DIROpt     = $(DirEVR)/sub_Optimization
DIRCRP     = $(DirEVR)/sub_CRP

DIRUT      = $(DIR_EVRT)/UnitTests

#============================================================================
#Libs, Minimize Only list: OK
# USE mod_system
Obj_Primlib  = \
  $(OBJ)/sub_module_NumParameters.o $(OBJ)/sub_module_MPI.o \
  $(OBJ)/sub_module_memory.o $(OBJ)/sub_module_string.o \
  $(OBJ)/sub_module_memory_Pointer.o $(OBJ)/sub_module_memory_NotPointer.o \
  $(OBJ)/sub_module_file.o $(OBJ)/sub_module_RW_MatVec.o $(OBJ)/mod_Frac.o \
  $(OBJ)/sub_module_system.o \
  $(OBJ)/sub_module_MPI_aux.o

Obj_math =\
   $(OBJ)/sub_diago.o $(OBJ)/sub_trans_mat.o $(OBJ)/sub_math_util.o $(OBJ)/sub_integration.o \
   $(OBJ)/sub_polyortho.o $(OBJ)/sub_function.o $(OBJ)/sub_derive.o $(OBJ)/sub_pert.o \
   $(OBJ)/sub_fft.o \
   $(OBJ)/Wigner3j.o

Obj_io = $(OBJ)/sub_io.o

# dnSVM, Minimize Only list: OK
# USE mod_dnSVM
Obj_dnSVM = \
  $(OBJ)/sub_module_dnS.o $(OBJ)/sub_module_VecOFdnS.o $(OBJ)/sub_module_MatOFdnS.o \
  $(OBJ)/sub_module_dnV.o $(OBJ)/sub_module_dnM.o $(OBJ)/sub_module_IntVM.o \
  $(OBJ)/sub_module_dnSVM.o

Obj_FiniteDiff = $(OBJ)/mod_FiniteDiff.o

# nDindex, Minimize Only list: OK
# USE mod_mod_nDindex and mod_module_DInd
Obj_nDindex  = $(OBJ)/sub_module_DInd.o $(OBJ)/sub_module_nDindex.o

# nDfit, Minimize Only list: OK
Obj_nDfit    = $(OBJ)/sub_module_nDfit.o

Obj_lib  = $(Obj_Primlib) $(Obj_math) $(Obj_io) $(Obj_dnSVM) \
           $(Obj_FiniteDiff) $(Obj_nDindex) $(Obj_nDfit)


#============================================================================

#============================================================================
#Smolyak test
Obj_Smolyak_test = \
  $(OBJ)/sub_Smolyak_DInd.o $(OBJ)/sub_Smolyak_ba.o $(OBJ)/sub_Smolyak_RDP.o \
  $(OBJ)/sub_Smolyak_module.o $(OBJ)/sub_Smolyak_test.o
#============================================================================

#============================================================================
#Physical constant, Minimize Only list: OK
#USE mod_constant
Obj_PhyCte = $(OBJ)/sub_module_RealWithUnit.o $(OBJ)/sub_module_Atom.o $(OBJ)/sub_module_constant.o
#============================================================================


#============================================================================
#KEO objects
#
#TanaPrim objects, Minimize Only list: OK
Obj_TanaPrim = $(OBJ)/sub_module_Tana_OpEl.o \
  $(OBJ)/sub_module_Tana_Op1D.o $(OBJ)/sub_module_Tana_OpnD.o \
  $(OBJ)/sub_module_Tana_SumOpnD.o $(OBJ)/sub_module_Tana_VecSumOpnD.o \
  $(OBJ)/sub_module_Tana_PiEulerRot.o


Obj_Coord = \
  $(OBJ)/Lib_QTransfo.o \
  $(OBJ)/BunchPolyTransfo.o $(OBJ)/ZmatTransfo.o $(OBJ)/QTOXanaTransfo.o $(OBJ)/CartesianTransfo.o \
  $(OBJ)/OneDTransfo.o $(OBJ)/ThreeDTransfo.o $(OBJ)/TwoDTransfo.o $(OBJ)/Rot2CoordTransfo.o \
  $(OBJ)/FlexibleTransfo.o $(OBJ)/GeneTransfo.o \
  $(OBJ)/HyperSpheTransfo.o $(OBJ)/LinearNMTransfo.o $(OBJ)/RectilinearNM_Transfo.o \
  $(OBJ)/sub_freq.o $(OBJ)/RPHTransfo.o $(OBJ)/RPHQMLTransfo.o $(OBJ)/ProjectTransfo.o \
  $(OBJ)/ActiveTransfo.o $(OBJ)/Qtransfo.o \
  $(OBJ)/Calc_Tab_dnQflex.o

#Minimize Only list: OK
Obj_Tnum = \
  $(OBJ)/sub_module_Tnum.o $(OBJ)/sub_module_paramQ.o \
  $(OBJ)/calc_f2_f1Q.o $(OBJ)/Sub_X_TO_Q_ana.o $(OBJ)/sub_dnDetGG_dnDetg.o $(OBJ)/sub_dnRho.o \
  $(OBJ)/calc_dng_dnGG.o $(OBJ)/sub_export_KEO.o

#Tana objects
Obj_Tana = \
  $(OBJ)/sub_module_Tana_vec_operations.o $(OBJ)/sub_module_Tana_op.o \
  $(OBJ)/sub_module_Tana_Export_KEO.o \
  $(OBJ)/sub_module_Tana_NumKEO.o $(OBJ)/sub_module_Tana_keo.o

#Minimize Only list: OK
Obj_TnumTana = $(OBJ)/calc_f2_f1Q_num.o $(OBJ)/sub_module_Tana_Tnum.o

#Minimize Only list: OK
Obj_Coord_KEO = $(Obj_TanaPrim) $(Obj_Coord) $(Obj_Tnum) $(Obj_Tana) $(Obj_TnumTana) $(OBJ)/sub_module_Coord_KEO.o
#============================================================================

#============================================================================
#Primitive Operators, Minimize Only list: OK
Obj_PrimOperator = \
   $(OBJ)/sub_module_SimpleOp.o $(OBJ)/sub_module_OnTheFly_def.o \
	 $(OBJ)/mod_CAP.o $(OBJ)/mod_HStep.o\
	 $(OBJ)/sub_PrimOp_def.o \
   $(OBJ)/sub_onthefly.o $(OBJ)/sub_PrimOp_RPH.o $(OBJ)/sub_PrimOp.o \
   $(OBJ)/sub_system.o $(OBJ)/read_para.o
#============================================================================

#===========================================================
#  For Tnum/Tana + Primitive Operators only
Obj_KEO_PrimOp= \
  $(Obj_lib) $(Obj_PhyCte) $(OBJ)/versionEVR-T.o \
  $(Obj_Coord_KEO) $(Obj_PrimOperator) $(OBJ)/Module_ForTnumTana_Driver.o $(OBJ)/TnumTana_Lib.o

#============================================================


#============================================================================
#Objects for ElVibRot
Obj_main   =  $(OBJ)/vib.o $(OBJ)/versionEVR-T.o $(OBJ)/cart.o \
  $(OBJ)/sub_main_Optimization.o $(OBJ)/sub_main_nDfit.o

Obj_EVR-Mod =  $(OBJ)/EVR-Module.o $(OBJ)/versionEVR-T.o

#Minimize Only list: sub_module_RotBasis, sub_module_basis_Grid_Param, sub_SymAbelian
#... sub_module_Basis_LTO_n,
Obj_module =  \
 $(OBJ)/sub_module_RotBasis.o $(OBJ)/sub_module_basis_Grid_Param.o $(OBJ)/sub_module_Basis_LTO_n.o \
 $(OBJ)/sub_SymAbelian.o \
 $(OBJ)/sub_module_param_SGType2.o $(OBJ)/sub_module_param_RD.o \
 $(OBJ)/sub_module_basis_set_alloc.o \
 $(OBJ)/sub_module_basis_RCVec_SG4.o \
 $(OBJ)/sub_module_basis_BtoG_GtoB_SG4_MPI.o $(OBJ)/sub_module_basis_BtoG_GtoB_SG4.o\
 $(OBJ)/sub_module_basis_BtoG_GtoB_SGType2.o \
 $(OBJ)/sub_module_basis_BtoG_GtoB_MPI.o $(OBJ)/sub_module_basis_BtoG_GtoB.o \
 $(OBJ)/sub_module_basis.o \
 $(OBJ)/sub_module_BasisMakeGrid.o \
 $(OBJ)/sub_module_poly.o $(OBJ)/sub_module_GWP.o \
 $(OBJ)/sub_module_cart.o


Obj_Basis = \
 $(OBJ)/sub_read_data.o \
 $(OBJ)/sub_quadra_inact.o \
 $(OBJ)/sub_basis_El.o \
 $(OBJ)/sub_quadra_herm.o $(OBJ)/sub_quadra_laguerre.o $(OBJ)/sub_quadra_legendre.o \
 $(OBJ)/sub_quadra_fourier.o $(OBJ)/sub_quadra_box.o $(OBJ)/sub_quadra_SincDVR.o $(OBJ)/sub_quadra_ft.o \
 $(OBJ)/sub_quadra_Ylm.o $(OBJ)/sub_quadra_Wigner.o $(OBJ)/sub_quadra_DirProd.o \
 $(OBJ)/sub_quadra_SparseBasis2n.o \
 $(OBJ)/sub_SymAbelian_OF_Basis.o

Obj_WP = \
 $(OBJ)/sub_module_type_ana_psi.o $(OBJ)/sub_module_param_WP0.o \
 $(OBJ)/sub_module_psi_set_alloc.o \
 $(OBJ)/sub_module_psi_B_TO_G.o \
 $(OBJ)/sub_module_ana_psi.o  $(OBJ)/sub_module_ana_psi_MPI.o \
 $(OBJ)/sub_module_psi_Op.o $(OBJ)/sub_module_psi_Op_MPI.o \
 $(OBJ)/sub_module_psi_io.o \
 $(OBJ)/mod_psi.o

Obj_propagation = \
 $(OBJ)/sub_module_field.o $(OBJ)/sub_module_ExactFact.o \
 $(OBJ)/sub_module_propagation.o $(OBJ)/sub_module_propagation_MPI.o \
 $(OBJ)/sub_module_propa_march_SG4.o \
 $(OBJ)/sub_module_propa_march_MPI.o $(OBJ)/sub_module_propa_march.o \
 $(OBJ)/sub_module_Filter.o \
 $(OBJ)/sub_module_Davidson_MPI.o $(OBJ)/sub_module_Davidson.o \
 $(OBJ)/sub_module_Arpack.o \
 $(OBJ)/sub_propagation.o $(OBJ)/sub_Hmax_MPI.o $(OBJ)/sub_Hmax.o $(OBJ)/sub_control.o \
 $(OBJ)/sub_TF_autocorr.o

ifeq ($(CERFACS),1)
  # CERFACS management
  Obj_CRP = $(OBJ)/sub_CRP.o $(OBJ)/CERFACS_lib.o $(OBJ)/QMRPACK_lib.o
else
  Obj_CRP = $(OBJ)/sub_CRP.o $(OBJ)/QMRPACK_lib.o
endif


Obj_inactive = \
 $(OBJ)/sub_HST_harm.o $(OBJ)/sub_inactive_harmo.o \
 $(OBJ)/sub_changement_de_var.o $(OBJ)/sub_ana_HS.o

Obj_active = \
 $(OBJ)/sub_Grid_SG4.o $(OBJ)/sub_ini_act_harm.o $(OBJ)/sub_lib_act.o \
 $(OBJ)/sub_diago_H.o $(OBJ)/sub_paraRPH.o

Obj_Operator = \
 $(OBJ)/sub_module_OpGrid.o $(OBJ)/sub_module_ReadOp.o $(OBJ)/sub_module_SetOp.o \
 $(OBJ)/sub_OpPsi_SG4.o $(OBJ)/sub_OpPsi_SG4_MPI.o \
 $(OBJ)/sub_OpPsi_MPI.o $(OBJ)/sub_OpPsi.o \
 $(OBJ)/sub_MatOp.o $(OBJ)/sub_lib_Op.o \
 $(OBJ)/sub_module_Op.o

Obj_analysis = \
 $(OBJ)/sub_module_analysis.o $(OBJ)/sub_analyse.o \
 $(OBJ)/sub_NLO.o $(OBJ)/sub_VibRot.o $(OBJ)/sub_intensity.o


Obj_Optimization = \
  $(OBJ)/sub_module_SimulatedAnnealing.o $(OBJ)/sub_module_BFGS.o $(OBJ)/sub_module_Optimization.o

Obj_ini = \
 $(OBJ)/ini_data.o $(OBJ)/sub_namelist.o $(OBJ)/nb_harm.o


Obj_Basis_WP_Op_propa = \
 $(OBJ)/sub_Auto_Basis.o $(OBJ)/sub_quadra_SparseBasis.o
#
Obj_EVRT =\
  $(Obj_lib) $(Obj_PhyCte) $(Obj_Coord_KEO) $(Obj_PrimOperator) \
  $(Obj_module) $(Obj_Basis) \
  $(Obj_WP) \
  $(Obj_Operator) \
  $(Obj_inactive) $(Obj_active) \
  $(Obj_propagation) $(Obj_CRP) \
  $(Obj_analysis) \
  $(Obj_Basis_WP_Op_propa) \
  $(Obj_Optimization) \
  $(Obj_ini) $(Obj_main) \
  $(Obj_Smolyak_test) \
  $(LIBARPACK)
#
#===============================================
#==============================================
#ElVibRot:

#make all : EVR
.PHONY: all evr EVR libEVR libevr
evr EVR all :obj vib $(VIBEXE)
	@echo "EVR OK"
libEVR libevr: obj $(OBJ)/libEVR.a
	@echo "libEVR.a OK"

#============================================================================
# All tnum/Tana ...
.PHONY: Tnum_FDriver Tnum_cDriver libTnum libTnum.a keotest
.PHONY: tnum Tnum tnum-dist Tnum-dist Tnum_MCTDH Tnum_MidasCpp Midas midas

Tnum_FDriver: obj qml $(Main_TnumTana_FDriverEXE)
	@echo "Main_TnumTana_FDriver OK"
Tnum_cDriver: obj qml $(Main_TnumTana_cDriverEXE)
	@echo "Main_TnumTana_cDriver OK"
#
libTnum libTnum.a: obj qml $(OBJ)/libTnum.a
	@echo "libTnum.a OK"
#
keotest: obj qml $(KEOTESTEXE)
	@echo "TEST_TnumTana OK"

tnum Tnum tnum-dist Tnum-dist: obj qml $(TNUMEXE)
	@echo "Tnum OK"
#
Tnum_MCTDH: obj qml $(TNUMMCTDHEXE)
	@echo "Tnum_MCTDH OK"
#
#TNUM_MiddasCppEXE
Tnum_MidasCpp Midas midas: obj qml $(TNUM_MiddasCppEXE)
	@echo "Tnum_MidasCpp OK"
#
.PHONY: Tana_test
Tana_test: Tana_test.exe
Tana_test.exe: obj qml $(Obj_lib) $(OBJ)/libTnum.a $(OBJ)/Tana_test.o
	$(LYNK90)   -o Tana_test.exe $(OBJ)/Tana_test.o $(OBJ)/libTnum.a $(LYNKFLAGS)
$(OBJ)/Tana_test.o: $(DirTNUM)/sub_main/Tana_test.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/sub_main/Tana_test.f90
#============================================================================
# Some all programs
.PHONY: gauss GWP work
gauss GWP: obj $(GWPEXE)
	@echo "GWP OK"
#
work:obj $(WORKEXE)
	@echo "work OK"
#============================================================================
# Physical Constants
.PHONY: PhysConst
PhysConst: obj $(PhysConstEXE)
	@echo "Physical Constants OK"
#============================================================================
# Unitary tests
.PHONY: ut UT UnitTests
ut UT UnitTests: UT_Frac UT_PhysConst UT_Tnum UT_HNO3 UT_HCN UT_HCN-WP
#
.PHONY: UT_Frac ut_frac
UT_Frac ut_frac : UnitTests_Frac.exe
	@echo "---------------------------------------"
	@echo "Unitary tests for the Frac module"
	@./UnitTests_Frac.exe > $(DIRUT)/res_UT_Frac ; awk -f $(DIRUT)/frac.awk $(DIRUT)/res_UT_Frac
	@echo "---------------------------------------"
UnitTests_Frac.exe: obj $(Obj_lib) $(OBJ)/UnitTests_Frac.o
	$(LYNK90)   -o UnitTests_Frac.exe $(OBJ)/UnitTests_Frac.o $(Obj_lib) $(LYNKFLAGS)
$(OBJ)/UnitTests_Frac.o: $(DIRUT)/UnitTests_Frac.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRUT)/UnitTests_Frac.f90
#
.PHONY: UT_PhysConst ut_physconst
UT_PhysConst ut_physconst: PhysConst
	@echo "---------------------------------------"
	@echo "Unitary tests for the PhysConst module"
	@cd Examples/exa_PhysicalConstants ; ./run_tests > $(DIRUT)/res_UT_PhysConst ; $(DIRUT)/PhysConst.sh $(DIRUT)/res_UT_PhysConst
	@echo "---------------------------------------"
#
.PHONY: UT_Tnum ut_Tnum UT_tnum ut_tnum
UT_Tnum ut_Tnum UT_tnum ut_tnum: Tnum
	@echo "---------------------------------------"
	@echo "Unitary tests for the Tnum"
	@cd UnitTests/Tnum_UT ; ./run_tests
	@echo "---------------------------------------"
#
.PHONY: UT_HNO3 ut_hno3
UT_HNO3 ut_hno3: EVR
	@echo "---------------------------------------"
	@echo "Unitary tests for the HNO3 ElVibRot calculations"
ifeq ($(F90),mpifort)
else
	@cd UnitTests/HNO3_UT ; ./run_tests small
endif
	@echo "---------------------------------------"
#
.PHONY: UT_HCN ut_hcn
UT_HCN ut_hcn: EVR
	@echo "---------------------------------------"
	@echo "Unitary tests for the HCN (diago) ElVibRot calculations"
ifeq ($(F90),mpifort)
	@cd UnitTests/HCN_UT_MPI ; ./run_tests small
else
	@cd UnitTests/HCN_UT ; ./run_tests small
endif
	@echo "---------------------------------------"
#
.PHONY: UT_HCN-WP ut_hcn-wp
UT_HCN-WP ut_hcn-wp: EVR
	@echo "---------------------------------------"
ifeq ($(F90),mpifort)
	@echo "Unitary tests for the pyrazine (propagation) ElVibRot calculations"
	@cd UnitTests/PYR-WP_UT_MPI ; ./run_tests small
else
	@echo "Unitary tests for the HCN (propagation) ElVibRot calculations"
	@cd UnitTests/HCN-WP_UT ; ./run_tests small
endif
	@echo "---------------------------------------"
#===============================================
#===============================================
#
# QML
QMLObjDIR   :=OBJ/obj_$(F90)_omp$(OMP)
QMLMODFILE= $(QMLObjDIR)/adiachannels_basis_m.mod $(QMLObjDIR)/irc_m.mod $(QMLObjDIR)/opt_m.mod $(QMLObjDIR)/model_m.mod

.PHONY: qml QML
qml QML: $(QMLibDIR) $(QMLibDIR_full)
	@echo "make qml library"
$(QMLibDIR_full): $(QMLibDIR)
	cd $(QMLibDIR) ; make lib
#cd $(QMLibDIR) ; make ; cp $(QMLMODFILE) $(OBJ)

$(QMLibDIR):
	cd $(ExternalLibDIR) ; ./get_QML.sh
	test -d $(QMLibDIR) || exit 1
.PHONY: clean_qml clean_QML
clean_qml clean_QML:
	cd $(ExternalLibDIR) ; rm -rf QuantumModelLib*
#
# dnS libraries
#
dnSVMObjDIR   :=OBJ/obj_$(F90)_omp$(OMP)

dnSMODFILE= $(dnSVMObjDIR)/addnsvm_m.mod $(dnSVMObjDIR)/addnsvm_dns_m.mod \
            $(dnSVMObjDIR)/addnsvm_dnmat_m.mod $(dnSVMObjDIR)/addnsvm_dnpoly_m.mod
.PHONY: dns dnS
dns dnS: $(dnSVMLibDIR) $(dnSVMLibDIR_full)
	@echo "make dnS library"
$(dnSVMLibDIR_full): $(dnSVMLibDIR)
	@echo "make dnS library"
	cd $(dnSVMLibDIR) ; make lib
#cd $(dnSVMLibDIR) ; make lib ; cp $(dnSMODFILE) $(OBJ)
$(dnSVMLibDIR):
	cd $(ExternalLibDIR) ; ./get_dnSVM.sh
	test -d $(dnSVMLibDIR) || exit 1

.PHONY: clean_dns clean_dnS
clean_dns clean_dnS:
	cd $(ExternalLibDIR) ; rm -rf AD_dnSVM* dnSVMLib
#
##################################################################################


# obj directory
.PHONY: obj
obj:
	@echo "=> obj directory: $(obj_dir)"
	@mkdir -p $(obj_dir)

# vib script
.PHONY: vib
vib:
	@echo "make vib script"
	./scripts/make_vib.sh $(DIR_EVRT) $(F90)
	chmod a+x vib

.PHONY: clean_UT
clean_UT:
	@cd UnitTests ; ./clean
	@echo "UnitTests cleaned"

# clean
.PHONY: clean
clean: clean_example clean_dnS clean_qml
	rm -f *.lst $(OBJ)/*.o *.mod *.MOD $(OBJ)/*.mod $(OBJ)/*.MOD $(EXE) *.exe $(OBJ)/*.a vib
	rm -f *.lst $(DIR_EVRT)/obj/*/*.o $(DIR_EVRT)/obj/*/*.mod $(DIR_EVRT)/obj/*/*.MOD $(EXE) *.exe $(DIR_EVRT)/obj/*/*.a vib
	rm -rf *.dSYM
	rm -f .DS_Store */.DS_Store */*/.DS_Store */*/*/.DS_Store
	@cd sub_pot                              ; rm -f sub_system.f sub_system.f90
	@cd Source_TnumTana_Coord/sub_operator_T ; rm -f calc_f2_f1Q.f90 Calc_Tab_dnQflex.f90 Sub_X_TO_Q_ana.f90
	@echo "  done remove the system dependent files (sub_system.f, calc_f2_f1Q.f90 ...) "
	@cd Examples/exa_hcn-dist ; ./clean
	@cd Examples/exa_TnumDriver ; ./clean
	@cd Examples/exa_direct-dist ; ./clean
	@cd Examples/exa_TnumTana_Coord-dist ; ./clean
	@cd Examples/exa_PhysicalConstants ; ./clean
	@cd UnitTests ; ./clean
	@cd Examples/exa_TanaCheck ; ./clean
	@cd Working_tests/exa_TnumTana-dist ; ./clean
	@cd Working_tests/exa_TnumPoly-dist ; ./clean
	@cd Working_tests/examples_Tana ; ./clean
	@cd Working_tests/exa_hcn-test ; ./clean
	@echo "  done cleaning up the example directories"
	@cd doc_ElVibRot-TnumTana && ./clean
	@rm -f doc/reference/index.store
	@echo "  done cleaning up the documentation"
#===============================================
#===============================================
#
$(VIBEXE): obj $(Obj_EVRT) $(OBJ)/$(VIBMAIN).o $(QMLibDIR_full) $(dnSVMLibDIR_full)
	@echo EVR-T
	$(LYNK90)   -o $(VIBEXE) $(Obj_EVRT) $(OBJ)/$(VIBMAIN).o $(LYNKFLAGS)
#	if test $(F90) = "pgf90" ; then mv $(VIBEXE) $(VIBEXE)2 ; echo "export OMP_STACKSIZE=50M" > $(VIBEXE) ; echo $(DIR_EVRT)/$(VIBEXE)2 >> $(VIBEXE) ; chmod a+x $(VIBEXE) ; fi
#===============================================
#
$(OBJ)/libTnum.a: obj $(Obj_KEO_PrimOp) $(QMLibDIR_full) $(dnSVMLibDIR_full)
	ar cr $(OBJ)/libTnum.a   $(Obj_KEO_PrimOp)
$(OBJ)/libEVR.a:obj $(Obj_EVRT) $(OBJ)/EVR_Module.o $(OBJ)/EVR_driver.o $(QMLibDIR_full) $(dnSVMLibDIR_full)
	ar cr $(OBJ)/libEVR.a $(Obj_EVRT)  $(OBJ)/EVR_Module.o $(OBJ)/EVR_driver.o
$(KEOTESTEXE): obj $(OBJ)/libTnum.a $(OBJ)/$(KEOTEST).o
	$(LYNK90)   -o $(KEOTESTEXE) $(OBJ)/$(KEOTEST).o $(OBJ)/libTnum.a $(LYNKFLAGS)

#Main_TnumTana_FDriver
$(Main_TnumTana_FDriverEXE): obj $(OBJ)/libTnum.a $(QMLibDIR_full) $(dnSVMLibDIR_full) $(OBJ)/Main_TnumTana_FDriver.o
	$(LYNK90)   -o $(Main_TnumTana_FDriverEXE) $(OBJ)/Main_TnumTana_FDriver.o $(OBJ)/libTnum.a $(LYNKFLAGS) $(dnSVMLIB)
$(Main_TnumTana_cDriverEXE): obj $(OBJ)/libTnum.a $(QMLibDIR_full) $(dnSVMLibDIR_full) $(OBJ)/Main_TnumTana_cDriver.o
	cp $(OBJ)/libTnum.a $(OBJ)/libTnumForcDriver.a
	ar d $(OBJ)/libTnumForcDriver.a sub_integration.o
	$(CompC) -o $(Main_TnumTana_cDriverEXE) $(CFLAGS) $(OBJ)/Main_TnumTana_cDriver.o $(OBJ)/libTnumForcDriver.a $(LYNKFLAGS) $(dnSVMLIB) -lgfortran -lm
#
$(TNUMEXE): obj $(OBJ)/libTnum.a $(QMLibDIR_full) $(dnSVMLibDIR_full) $(OBJ)/$(TNUMMAIN).o
	$(LYNK90)   -o $(TNUMEXE) $(OBJ)/$(TNUMMAIN).o $(OBJ)/libTnum.a $(LYNKFLAGS) $(dnSVMLIB)
#
$(TNUMMCTDHEXE): obj $(OBJ)/libTnum.a $(QMLibDIR_full) $(dnSVMLibDIR_full) $(OBJ)/$(TNUMMCTDHMAIN).o
	$(LYNK90)   -o $(TNUMMCTDHEXE) $(OBJ)/$(TNUMMCTDHMAIN).o $(OBJ)/libTnum.a $(LYNKFLAGS) $(dnSVMLIB)
# TNUM_MiddasCppEXE  = Tnum90_MidasCpp.exe
# TNUM_MiddasCppMAIN = Tnum90_MidasCpp
$(TNUM_MiddasCppEXE): obj $(OBJ)/libTnum.a $(QMLibDIR_full) $(dnSVMLibDIR_full) $(OBJ)/$(TNUM_MiddasCppMAIN).o
	$(LYNK90)   -o $(TNUM_MiddasCppEXE) $(OBJ)/$(TNUM_MiddasCppMAIN).o  $(OBJ)/libTnum.a $(LYNKFLAGS) $(dnSVMLIB)
#
$(GWPEXE): obj $(Obj_All) $(OBJ)/$(GWPMAIN).o
	$(LYNK90)   -o $(GWPEXE) $(Obj_All) $(OBJ)/$(GWPMAIN).o  $(LYNKFLAGS)
#
$(WORKEXE): obj $(Obj_KEO_PrimOp) $(QMLibDIR_full) $(dnSVMLibDIR_full) $(OBJ)/$(WORKMAIN).o
	$(LYNK90)   -o $(WORKEXE) $(Obj_KEO_PrimOp) $(OBJ)/$(WORKMAIN).o  $(LYNKFLAGS) $(dnSVMLIB)
#$(Obj_Primlib)
#===============================================PhysConst:obj $(PhysConstEXE)
#$(PhysConstEXE): obj $(Obj_lib) $(Obj_PhyCte) $(OBJ)/$(PhysConstMAIN).o
#	$(LYNK90)   -o $(PhysConstEXE) $(Obj_lib) $(Obj_PhyCte) $(OBJ)/$(PhysConstMAIN).o  $(LYNKFLAGS)
$(PhysConstEXE): obj $(Obj_Primlib) $(Obj_math) $(Obj_io) $(Obj_PhyCte) $(OBJ)/$(PhysConstMAIN).o
	ar cr $(OBJ)/libPhysConst.a $(Obj_Primlib) $(Obj_math) $(Obj_io) $(Obj_PhyCte)
	$(LYNK90)   -o $(PhysConstEXE) $(Obj_Primlib) $(Obj_math) $(Obj_io) $(Obj_PhyCte) $(OBJ)/$(PhysConstMAIN).o  $(LYNKFLAGS)
#===================================================================================
# lib
$(OBJ)/sub_module_MPI.o:$(DirSys)/sub_module_MPI.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DirSys)/sub_module_MPI.f90
$(OBJ)/sub_module_NumParameters.o:$(DirSys)/sub_module_NumParameters.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirSys)/sub_module_NumParameters.f90
$(OBJ)/mod_Frac.o:$(DirSys)/mod_Frac.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirSys)/mod_Frac.f90
$(OBJ)/sub_module_memory.o:$(DirSys)/sub_module_memory.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirSys)/sub_module_memory.f90
$(OBJ)/sub_module_memory_Pointer.o:$(DirSys)/sub_module_memory_Pointer.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirSys)/sub_module_memory_Pointer.f90
$(OBJ)/sub_module_memory_NotPointer.o:$(DirSys)/sub_module_memory_NotPointer.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirSys)/sub_module_memory_NotPointer.f90
$(OBJ)/sub_module_file.o:$(DirSys)/sub_module_file.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirSys)/sub_module_file.f90
$(OBJ)/sub_module_string.o:$(DirSys)/sub_module_string.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DirSys)/sub_module_string.f90
$(OBJ)/sub_module_RW_MatVec.o:$(DirSys)/sub_module_RW_MatVec.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPSHELL)  -c $(DirSys)/sub_module_RW_MatVec.f90
$(OBJ)/sub_module_system.o:$(DirSys)/sub_module_system.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) $(CPPSHELL) $(CPPSHELL_GIT) -c $(DirSys)/sub_module_system.f90
$(OBJ)/sub_module_MPI_aux.o:$(DirSys)/sub_module_MPI_aux.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DirSys)/sub_module_MPI_aux.f90
###
$(OBJ)/sub_module_DInd.o:$(DirnDind)/sub_module_DInd.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirnDind)/sub_module_DInd.f90
$(OBJ)/sub_module_nDindex.o:$(DirnDind)/sub_module_nDindex.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirnDind)/sub_module_nDindex.f90
###
$(OBJ)/sub_module_dnS.o:$(DirdnSVM)/sub_module_dnS.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirdnSVM)/sub_module_dnS.f90
$(OBJ)/sub_module_VecOFdnS.o:$(DirdnSVM)/sub_module_VecOFdnS.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirdnSVM)/sub_module_VecOFdnS.f90
$(OBJ)/sub_module_MatOFdnS.o:$(DirdnSVM)/sub_module_MatOFdnS.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirdnSVM)/sub_module_MatOFdnS.f90
$(OBJ)/sub_module_dnV.o:$(DirdnSVM)/sub_module_dnV.f90 $(dnSVMLibDIR_full)
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirdnSVM)/sub_module_dnV.f90
$(OBJ)/sub_module_dnM.o:$(DirdnSVM)/sub_module_dnM.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirdnSVM)/sub_module_dnM.f90
$(OBJ)/sub_module_IntVM.o:$(DirdnSVM)/sub_module_IntVM.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirdnSVM)/sub_module_IntVM.f90
$(OBJ)/sub_module_dnSVM.o:$(DirdnSVM)/sub_module_dnSVM.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirdnSVM)/sub_module_dnSVM.f90
#
$(OBJ)/mod_FiniteDiff.o:$(DirSys)/mod_FiniteDiff.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirSys)/mod_FiniteDiff.f90
#
$(OBJ)/sub_module_nDfit.o:$(DirMod)/sub_module_nDfit.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirMod)/sub_module_nDfit.f90
#
#===================================================================================
# Physical constants
$(OBJ)/sub_module_RealWithUnit.o:$(DIRPhyCte)/sub_module_RealWithUnit.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIRPhyCte)/sub_module_RealWithUnit.f90
$(OBJ)/sub_module_Atom.o:$(DIRPhyCte)/sub_module_Atom.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIRPhyCte)/sub_module_Atom.f90
$(OBJ)/sub_module_constant.o:$(DIRPhyCte)/sub_module_constant.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRPhyCte)/sub_module_constant.f90
$(OBJ)/$(PhysConstMAIN).o:$(DIRPhyCte)/$(PhysConstMAIN).f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRPhyCte)/$(PhysConstMAIN).f90
#
#===================================================================================
# Tnum Tana Coord
##
# TanaPrim files
$(OBJ)/sub_module_Tana_OpEl.o:$(DirTNUM)/Tana/sub_module_Tana_OpEl.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Tana/sub_module_Tana_OpEl.f90
$(OBJ)/sub_module_Tana_Op1D.o:$(DirTNUM)/Tana/sub_module_Tana_Op1D.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Tana/sub_module_Tana_Op1D.f90
$(OBJ)/sub_module_Tana_OpnD.o:$(DirTNUM)/Tana/sub_module_Tana_OpnD.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Tana/sub_module_Tana_OpnD.f90
$(OBJ)/sub_module_Tana_SumOpnD.o:$(DirTNUM)/Tana/sub_module_Tana_SumOpnD.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Tana/sub_module_Tana_SumOpnD.f90
$(OBJ)/sub_module_Tana_VecSumOpnD.o:$(DirTNUM)/Tana/sub_module_Tana_VecSumOpnD.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Tana/sub_module_Tana_VecSumOpnD.f90
$(OBJ)/sub_module_Tana_PiEulerRot.o:$(DirTNUM)/Tana/sub_module_Tana_PiEulerRot.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Tana/sub_module_Tana_PiEulerRot.f90
#
# Coordinates , Qtransfo, zmat...
#
$(OBJ)/Lib_QTransfo.o:$(DirTNUM)/Qtransfo/Lib_QTransfo.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Qtransfo/Lib_QTransfo.f90
$(OBJ)/CartesianTransfo.o:$(DirTNUM)/Qtransfo/CartesianTransfo.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Qtransfo/CartesianTransfo.f90
$(OBJ)/BunchPolyTransfo.o:$(DirTNUM)/Qtransfo/BunchPolyTransfo.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Qtransfo/BunchPolyTransfo.f90
$(OBJ)/ZmatTransfo.o:$(DirTNUM)/Qtransfo/ZmatTransfo.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Qtransfo/ZmatTransfo.f90
$(OBJ)/QTOXanaTransfo.o:$(DirTNUM)/Qtransfo/QTOXanaTransfo.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Qtransfo/QTOXanaTransfo.f90
$(OBJ)/OneDTransfo.o:$(DirTNUM)/Qtransfo/OneDTransfo.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Qtransfo/OneDTransfo.f90
$(OBJ)/ThreeDTransfo.o:$(DirTNUM)/Qtransfo/ThreeDTransfo.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Qtransfo/ThreeDTransfo.f90
$(OBJ)/TwoDTransfo.o:$(DirTNUM)/Qtransfo/TwoDTransfo.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Qtransfo/TwoDTransfo.f90
$(OBJ)/Rot2CoordTransfo.o:$(DirTNUM)/Qtransfo/Rot2CoordTransfo.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Qtransfo/Rot2CoordTransfo.f90
$(OBJ)/FlexibleTransfo.o:$(DirTNUM)/Qtransfo/FlexibleTransfo.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Qtransfo/FlexibleTransfo.f90
$(OBJ)/GeneTransfo.o:$(DirTNUM)/Qtransfo/GeneTransfo.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Qtransfo/GeneTransfo.f90
$(OBJ)/HyperSpheTransfo.o:$(DirTNUM)/Qtransfo/HyperSpheTransfo.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Qtransfo/HyperSpheTransfo.f90
$(OBJ)/LinearNMTransfo.o:$(DirTNUM)/Qtransfo/LinearNMTransfo.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Qtransfo/LinearNMTransfo.f90
$(OBJ)/ProjectTransfo.o:$(DirTNUM)/Qtransfo/ProjectTransfo.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Qtransfo/ProjectTransfo.f90
$(OBJ)/RectilinearNM_Transfo.o:$(DirTNUM)/Qtransfo/RectilinearNM_Transfo.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Qtransfo/RectilinearNM_Transfo.f90
$(OBJ)/RPHTransfo.o:$(DirTNUM)/Qtransfo/RPHTransfo.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Qtransfo/RPHTransfo.f90
$(OBJ)/RPHQMLTransfo.o:$(DirTNUM)/Qtransfo/RPHQMLTransfo.f90 $(QMLibDIR_full) $(dnSVMLibDIR_full)
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Qtransfo/RPHQMLTransfo.f90
$(OBJ)/ActiveTransfo.o:$(DirTNUM)/Qtransfo/ActiveTransfo.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Qtransfo/ActiveTransfo.f90
$(OBJ)/Qtransfo.o:$(DirTNUM)/Qtransfo/Qtransfo.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DirTNUM)/Qtransfo/Qtransfo.f90
$(OBJ)/sub_freq.o:$(DirTNUM)/Qtransfo/sub_freq.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Qtransfo/sub_freq.f90
$(OBJ)/Calc_Tab_dnQflex.o:$(DirTNUM)/sub_operator_T/Calc_Tab_dnQflex.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/sub_operator_T/Calc_Tab_dnQflex.f90
#
$(OBJ)/sub_module_Tnum.o:$(DirTNUM)/sub_module_Tnum.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DirTNUM)/sub_module_Tnum.f90
$(OBJ)/sub_module_paramQ.o:$(DirTNUM)/sub_module_paramQ.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DirTNUM)/sub_module_paramQ.f90
$(OBJ)/sub_module_Coord_KEO.o:$(DirTNUM)/sub_module_Coord_KEO.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/sub_module_Coord_KEO.f90
#
#
# Tnum files
$(OBJ)/calc_f2_f1Q.o:$(DirTNUM)/sub_operator_T/calc_f2_f1Q.f90
	sed "s/zmatrix/CoordType/" $(DirTNUM)/sub_operator_T/calc_f2_f1Q.f90 > $(DirTNUM)/sub_operator_T/calc_f2_f1Q.i
	sed "s/Write_mole/Write_CoordType/" $(DirTNUM)/sub_operator_T/calc_f2_f1Q.i > $(DirTNUM)/sub_operator_T/calc_f2_f1Q.f90
	@echo Warning the calc_f2_f1Q.f90 file has been modified.
	rm $(DirTNUM)/sub_operator_T/calc_f2_f1Q.i
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/sub_operator_T/calc_f2_f1Q.f90
$(OBJ)/Sub_X_TO_Q_ana.o:$(DirTNUM)/sub_operator_T/Sub_X_TO_Q_ana.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/sub_operator_T/Sub_X_TO_Q_ana.f90
#
# copy Sub_X_TO_Q_ana_save.f90, Calc_Tab_dnQflex_save.f90 and calc_f2_f1Q_save.f90
# to Sub_X_TO_Q_ana.f90, Calc_Tab_dnQflex.f90 and calc_f2_f1Q.f90 when they are not present
#
$(DirTNUM)/sub_operator_T/Sub_X_TO_Q_ana.f90:
	cp $(DirTNUM)/sub_operator_T/Sub_X_TO_Q_ana_save.f90 $(DirTNUM)/sub_operator_T/Sub_X_TO_Q_ana.f90
$(DirTNUM)/sub_operator_T/Calc_Tab_dnQflex.f90:
	cp $(DirTNUM)/sub_operator_T/Calc_Tab_dnQflex_save.f90 $(DirTNUM)/sub_operator_T/Calc_Tab_dnQflex.f90
$(DirTNUM)/sub_operator_T/calc_f2_f1Q.f90:
	cp $(DirTNUM)/sub_operator_T/calc_f2_f1Q_save.f90 $(DirTNUM)/sub_operator_T/calc_f2_f1Q.f90
#
$(OBJ)/calc_f2_f1Q_num.o:$(DirTNUM)/Tnum/calc_f2_f1Q_num.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Tnum/calc_f2_f1Q_num.f90
$(OBJ)/calc_dng_dnGG.o:$(DirTNUM)/Tnum/calc_dng_dnGG.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Tnum/calc_dng_dnGG.f90
$(OBJ)/sub_export_KEO.o:$(DirTNUM)/Tnum/sub_export_KEO.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Tnum/sub_export_KEO.f90
$(OBJ)/sub_dnDetGG_dnDetg.o:$(DirTNUM)/Tnum/sub_dnDetGG_dnDetg.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Tnum/sub_dnDetGG_dnDetg.f90
$(OBJ)/sub_dnRho.o:$(DirTNUM)/Tnum/sub_dnRho.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DirTNUM)/Tnum/sub_dnRho.f90
#
# Tana files
$(OBJ)/sub_module_Tana_Export_KEO.o:$(DirTNUM)/Tana/sub_module_Tana_Export_KEO.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Tana/sub_module_Tana_Export_KEO.f90
$(OBJ)/sub_module_Tana_vec_operations.o:$(DirTNUM)/Tana/sub_module_Tana_vec_operations.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Tana/sub_module_Tana_vec_operations.f90
$(OBJ)/sub_module_Tana_op.o:$(DirTNUM)/Tana/sub_module_Tana_op.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Tana/sub_module_Tana_op.f90
$(OBJ)/sub_module_Tana_keo.o:$(DirTNUM)/Tana/sub_module_Tana_keo.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Tana/sub_module_Tana_keo.f90
$(OBJ)/sub_module_Tana_NumKEO.o:$(DirTNUM)/Tana/sub_module_Tana_NumKEO.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Tana/sub_module_Tana_NumKEO.f90
$(OBJ)/sub_module_Tana_Tnum.o:$(DirTNUM)/Tana/sub_module_Tana_Tnum.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Tana/sub_module_Tana_Tnum.f90
#
#===================================================================================
# Primitive Operators
$(OBJ)/sub_module_SimpleOp.o:$(DIRPrimOp)/sub_module_SimpleOp.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRPrimOp)/sub_module_SimpleOp.f90
$(OBJ)/sub_PrimOp_def.o:$(DIRPrimOp)/sub_PrimOp_def.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRPrimOp)/sub_PrimOp_def.f90
$(OBJ)/sub_module_OnTheFly_def.o:$(DIRPrimOp)/sub_module_OnTheFly_def.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRPrimOp)/sub_module_OnTheFly_def.f90
$(OBJ)/mod_CAP.o:$(DIRPrimOp)/mod_CAP.f90
	cd $(OBJ) ; $(F90_FLAGS)  -c $(DIRPrimOp)/mod_CAP.f90
$(OBJ)/mod_HStep.o:$(DIRPrimOp)/mod_HStep.f90
	cd $(OBJ) ; $(F90_FLAGS)  -c $(DIRPrimOp)/mod_HStep.f90
$(OBJ)/sub_PrimOp_RPH.o:$(DIRPrimOp)/sub_PrimOp_RPH.f90
	cd $(OBJ) ; $(F90_FLAGS)  -c $(DIRPrimOp)/sub_PrimOp_RPH.f90
$(OBJ)/sub_PrimOp.o:$(DIRPrimOp)/sub_PrimOp.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) $(CPPSHELL_QML)  -c $(DIRPrimOp)/sub_PrimOp.f90
$(OBJ)/sub_onthefly.o:$(DIRPrimOp)/sub_onthefly.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRPrimOp)/sub_onthefly.f90
$(OBJ)/Module_ForTnumTana_Driver.o:$(DirTNUM)/Module_ForTnumTana_Driver.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Module_ForTnumTana_Driver.f90
$(OBJ)/TnumTana_Lib.o:$(DirTNUM)/TnumTana_Lib.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/TnumTana_Lib.f90
#
#===================================================================================
# mains TnumTana_PrimOp
$(OBJ)/$(TNUMMAIN).o:$(DirTNUM)/$(TNUMMAIN).f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/$(TNUMMAIN).f90
$(OBJ)/$(TNUMDISTMAIN).o:$(DirTNUM)/$(TNUMDISTMAIN).f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/$(TNUMDISTMAIN).f90
$(OBJ)/$(TNUMMCTDHMAIN).o:$(DirTNUM)/$(TNUMMCTDHMAIN).f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/$(TNUMMCTDHMAIN).f90
#Emil change: Could not compile without inserting the following two lines of code
$(OBJ)/$(TNUM_MiddasCppMAIN).o:$(DirTNUM)/$(TNUM_MiddasCppMAIN).f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/$(TNUM_MiddasCppMAIN).f90
$(OBJ)/$(KEOTEST).o:$(DirTNUM)/$(KEOTEST).f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/$(KEOTEST).f90
$(OBJ)/Main_TnumTana_FDriver.o:$(DirTNUM)/Main_TnumTana_FDriver.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/Main_TnumTana_FDriver.f90
$(OBJ)/Main_TnumTana_cDriver.o:$(DirTNUM)/Main_TnumTana_cDriver.c
	cd $(OBJ) ; $(CompC) $(CFLAGS)  -c $(DirTNUM)/Main_TnumTana_cDriver.c
#
#===============================================================================
# module
$(OBJ)/sub_module_poly.o:$(DirMod)/sub_module_poly.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirMod)/sub_module_poly.f90
$(OBJ)/sub_module_GWP.o:$(DirMod)/sub_module_GWP.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirMod)/sub_module_GWP.f90
$(OBJ)/sub_module_cart.o:$(DirMod)/sub_module_cart.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirMod)/sub_module_cart.f90
#
#===============================================================================
# sub_Basis :
$(OBJ)/sub_module_RotBasis.o:$(DIRba)/sub_module_RotBasis.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRba)/sub_module_RotBasis.f90
$(OBJ)/sub_module_basis_Grid_Param.o:$(DIRba)/sub_module_basis_Grid_Param.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRba)/sub_module_basis_Grid_Param.f90
$(OBJ)/sub_module_Basis_LTO_n.o:$(DIRba)/sub_module_Basis_LTO_n.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRba)/sub_module_Basis_LTO_n.f90
$(OBJ)/sub_SymAbelian.o:$(DIRba)/sub_SymAbelian/sub_SymAbelian.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRba)/sub_SymAbelian/sub_SymAbelian.f90
$(OBJ)/sub_module_param_SGType2.o:$(DIRba)/sub_module_param_SGType2.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIRba)/sub_module_param_SGType2.f90
$(OBJ)/sub_module_param_RD.o:$(DIRba)/sub_ReducedDensity/sub_module_param_RD.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRba)/sub_ReducedDensity/sub_module_param_RD.f90
$(OBJ)/sub_module_basis_set_alloc.o:$(DIRba)/sub_module_basis_set_alloc.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRba)/sub_module_basis_set_alloc.f90

$(OBJ)/sub_module_basis_RCVec_SG4.o:$(DIRbaSG4)/sub_module_basis_RCVec_SG4.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRbaSG4)/sub_module_basis_RCVec_SG4.f90
$(OBJ)/sub_module_basis_BtoG_GtoB_SG4.o:$(DIRbaSG4)/sub_module_basis_BtoG_GtoB_SG4.f90
	cd $(OBJ) ; $(F90_FLAGS)  $(CPPpre) -c $(DIRbaSG4)/sub_module_basis_BtoG_GtoB_SG4.f90
$(OBJ)/sub_module_basis_BtoG_GtoB_SG4_MPI.o:$(DIRbaSG4)/sub_module_basis_BtoG_GtoB_SG4_MPI.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIRbaSG4)/sub_module_basis_BtoG_GtoB_SG4_MPI.f90

$(OBJ)/sub_module_basis_BtoG_GtoB_SGType2.o:$(DIRba)/sub_module_basis_BtoG_GtoB_SGType2.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRba)/sub_module_basis_BtoG_GtoB_SGType2.f90
$(OBJ)/sub_module_basis_BtoG_GtoB_MPI.o:$(DIRba)/sub_module_basis_BtoG_GtoB_MPI.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIRba)/sub_module_basis_BtoG_GtoB_MPI.f90
$(OBJ)/sub_module_basis_BtoG_GtoB.o:$(DIRba)/sub_module_basis_BtoG_GtoB.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRba)/sub_module_basis_BtoG_GtoB.f90
$(OBJ)/sub_module_basis.o:$(DIRba)/sub_module_basis.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRba)/sub_module_basis.f90
$(OBJ)/sub_Auto_Basis.o:$(DIRba)/sub_Auto_Basis.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIRba)/sub_Auto_Basis.f90
$(OBJ)/sub_read_data.o:$(DIRba)/sub_read_data.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRba)/sub_read_data.f90
$(OBJ)/sub_basis_El.o:$(DIRba)/sub_basis_El.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRba)/sub_basis_El.f90
$(OBJ)/sub_quadra_inact.o:$(DIRba)/sub_quadra_inact.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRba)/sub_quadra_inact.f90
$(OBJ)/sub_quadra_herm.o:$(DIRba)/sub_quadra_herm.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRba)/sub_quadra_herm.f90
$(OBJ)/sub_quadra_laguerre.o:$(DIRba)/sub_quadra_laguerre.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRba)/sub_quadra_laguerre.f90
$(OBJ)/sub_quadra_legendre.o:$(DIRba)/sub_quadra_legendre.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRba)/sub_quadra_legendre.f90
$(OBJ)/sub_quadra_fourier.o:$(DIRba)/sub_quadra_fourier.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRba)/sub_quadra_fourier.f90
$(OBJ)/sub_quadra_box.o:$(DIRba)/sub_quadra_box.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRba)/sub_quadra_box.f90
$(OBJ)/sub_quadra_SincDVR.o:$(DIRba)/sub_quadra_SincDVR.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRba)/sub_quadra_SincDVR.f90
$(OBJ)/sub_quadra_ft.o:$(DIRba)/sub_quadra_ft.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRba)/sub_quadra_ft.f90
$(OBJ)/sub_quadra_Ylm.o:$(DIRba)/sub_quadra_Ylm.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRba)/sub_quadra_Ylm.f90
$(OBJ)/sub_quadra_Wigner.o:$(DIRba)/sub_quadra_Wigner.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRba)/sub_quadra_Wigner.f90
$(OBJ)/sub_quadra_DirProd.o:$(DIRba)/sub_quadra_DirProd.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRba)/sub_quadra_DirProd.f90
$(OBJ)/sub_quadra_SparseBasis.o:$(DIRba)/sub_quadra_SparseBasis.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIRba)/sub_quadra_SparseBasis.f90
$(OBJ)/sub_module_BasisMakeGrid.o:$(DIRba)/sub_module_BasisMakeGrid.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRba)/sub_module_BasisMakeGrid.f90
$(OBJ)/sub_quadra_SparseBasis2n.o:$(DIRba)/sub_quadra_SparseBasis2n.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRba)/sub_quadra_SparseBasis2n.f90
$(OBJ)/sub_SymAbelian_OF_Basis.o:$(DIRba)/sub_SymAbelian_OF_Basis.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRba)/sub_SymAbelian_OF_Basis.f90
#
#===============================================================================
# sub_WP
$(OBJ)/sub_module_type_ana_psi.o:$(DIRWP)/sub_module_type_ana_psi.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRWP)/sub_module_type_ana_psi.f90
$(OBJ)/sub_module_param_WP0.o:$(DIRWP)/sub_module_param_WP0.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRWP)/sub_module_param_WP0.f90
$(OBJ)/sub_module_psi_set_alloc.o:$(DIRWP)/sub_module_psi_set_alloc.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIRWP)/sub_module_psi_set_alloc.f90
$(OBJ)/sub_module_ana_psi.o:$(DIRWP)/sub_module_ana_psi.f90
	cd $(OBJ) ; $(F90_FLAGS) -c $(DIRWP)/sub_module_ana_psi.f90
$(OBJ)/sub_module_ana_psi_MPI.o:$(DIRWP)/sub_module_ana_psi_MPI.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIRWP)/sub_module_ana_psi_MPI.f90
$(OBJ)/sub_module_psi_B_TO_G.o:$(DIRWP)/sub_module_psi_B_TO_G.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRWP)/sub_module_psi_B_TO_G.f90
$(OBJ)/sub_module_psi_Op.o:$(DIRWP)/sub_module_psi_Op.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRWP)/sub_module_psi_Op.f90
$(OBJ)/sub_module_psi_Op_MPI.o:$(DIRWP)/sub_module_psi_Op_MPI.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIRWP)/sub_module_psi_Op_MPI.f90
$(OBJ)/sub_module_psi_io.o:$(DIRWP)/sub_module_psi_io.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIRWP)/sub_module_psi_io.f90
$(OBJ)/mod_psi.o:$(DIRWP)/mod_psi.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIRWP)/mod_psi.f90
$(OBJ)/sub_ana_psi.o:$(DIRWP)/sub_ana_psi.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRWP)/sub_ana_psi.f90
#
#===============================================================================
# sub_propagation sub_module_ExactFact
$(OBJ)/sub_module_field.o:$(DIRpropa)/sub_module_field.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRpropa)/sub_module_field.f90
$(OBJ)/sub_module_ExactFact.o:$(DIRpropa)/sub_module_ExactFact.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRpropa)/sub_module_ExactFact.f90
$(OBJ)/sub_module_propagation.o:$(DIRpropa)/sub_module_propagation.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRpropa)/sub_module_propagation.f90
$(OBJ)/sub_module_propagation_MPI.o:$(DIRpropa)/sub_module_propagation_MPI.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIRpropa)/sub_module_propagation_MPI.f90
$(OBJ)/sub_module_propa_march.o:$(DIRpropa)/sub_module_propa_march.f90
	cd $(OBJ) ; $(F90_FLAGS) -c $(DIRpropa)/sub_module_propa_march.f90
$(OBJ)/sub_module_propa_march_MPI.o:$(DIRpropa)/sub_module_propa_march_MPI.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIRpropa)/sub_module_propa_march_MPI.f90
$(OBJ)/sub_module_propa_march_SG4.o:$(DIRpropa)/sub_module_propa_march_SG4.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRpropa)/sub_module_propa_march_SG4.f90
$(OBJ)/sub_module_Filter.o:$(DIRpropa)/sub_module_Filter.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRpropa)/sub_module_Filter.f90
$(OBJ)/sub_module_Davidson.o:$(DIRpropa)/sub_module_Davidson.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRpropa)/sub_module_Davidson.f90
$(OBJ)/sub_module_Davidson_MPI.o:$(DIRpropa)/sub_module_Davidson_MPI.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIRpropa)/sub_module_Davidson_MPI.f90
$(OBJ)/sub_module_Arpack.o:$(DIRpropa)/sub_module_Arpack.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) $(CPPSHELL_ARPACK)  -c $(DIRpropa)/sub_module_Arpack.f90
$(OBJ)/sub_propagation.o:$(DIRpropa)/sub_propagation.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIRpropa)/sub_propagation.f90
$(OBJ)/sub_Hmax.o:$(DIRpropa)/sub_Hmax.f90
	cd $(OBJ) ; $(F90_FLAGS) -c $(DIRpropa)/sub_Hmax.f90
$(OBJ)/sub_Hmax_MPI.o:$(DIRpropa)/sub_Hmax_MPI.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIRpropa)/sub_Hmax_MPI.f90
$(OBJ)/sub_control.o:$(DIRpropa)/sub_control.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRpropa)/sub_control.f90
$(OBJ)/sub_TF_autocorr.o:$(DIRpropa)/sub_TF_autocorr.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRpropa)/sub_TF_autocorr.f90
#
#===============================================================================
# sub_CRP:
$(OBJ)/sub_CRP.o:$(DIRCRP)/sub_CRP.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) $(CPPSHELL_CERFACS) $(CPPSHELL_ARPACK) -c $(DIRCRP)/sub_CRP.f90
$(OBJ)/CERFACS_lib.o:$(DIRCRP)/CERFACS_lib.f
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRCRP)/CERFACS_lib.f
$(OBJ)/QMRPACK_lib.o:$(DIRCRP)/QMRPACK_lib.f
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRCRP)/QMRPACK_lib.f
#===============================================================================
#Operator ....
$(OBJ)/sub_module_OpGrid.o:$(DIROp)/sub_module_OpGrid.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIROp)/sub_module_OpGrid.f90
$(OBJ)/sub_module_ReadOp.o:$(DIROp)/sub_module_ReadOp.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIROp)/sub_module_ReadOp.f90
$(OBJ)/sub_module_SetOp.o:$(DIROp)/sub_module_SetOp.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIROp)/sub_module_SetOp.f90
#
$(OBJ)/sub_HST_harm.o:$(DIR2)/sub_HST_harm.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIR2)/sub_HST_harm.f90
$(OBJ)/sub_inactive_harmo.o:$(DIR2)/sub_inactive_harmo.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIR2)/sub_inactive_harmo.f90
$(OBJ)/sub_changement_de_var.o:$(DIR2)/sub_changement_de_var.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIR2)/sub_changement_de_var.f90
$(OBJ)/sub_ana_HS.o:$(DIR2)/sub_ana_HS.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIR2)/sub_ana_HS.f90
#
# sub_active
$(OBJ)/sub_ini_act_harm.o:$(DIR5)/sub_ini_act_harm.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIR5)/sub_ini_act_harm.f90
$(OBJ)/sub_Grid_SG4.o:$(DIR5)/sub_Grid_SG4.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIR5)/sub_Grid_SG4.f90
$(OBJ)/sub_lib_act.o:$(DIR5)/sub_lib_act.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIR5)/sub_lib_act.f90
$(OBJ)/sub_diago_H.o:$(DIR5)/sub_diago_H.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIR5)/sub_diago_H.f90
$(OBJ)/sub_paraRPH.o:$(DIR5)/sub_paraRPH.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIR5)/sub_paraRPH.f90
#
# Operator
$(OBJ)/sub_MatOp.o:$(DIROp)/sub_MatOp.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIROp)/sub_MatOp.f90
$(OBJ)/sub_OpPsi_SG4.o:$(DIROp)/sub_OpPsi_SG4.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIROp)/sub_OpPsi_SG4.f90
$(OBJ)/sub_OpPsi_SG4_MPI.o:$(DIROp)/sub_OpPsi_SG4_MPI.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIROp)/sub_OpPsi_SG4_MPI.f90
$(OBJ)/sub_OpPsi.o:$(DIROp)/sub_OpPsi.f90
	cd $(OBJ) ; $(F90_FLAGS) -c $(DIROp)/sub_OpPsi.f90
$(OBJ)/sub_OpPsi_MPI.o:$(DIROp)/sub_OpPsi_MPI.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIROp)/sub_OpPsi_MPI.f90
$(OBJ)/sub_lib_Op.o:$(DIROp)/sub_lib_Op.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIROp)/sub_lib_Op.f90
$(OBJ)/sub_module_Op.o:$(DIROp)/sub_module_Op.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIROp)/sub_module_Op.f90
#
#===============================================================================
# sub_analysis
$(OBJ)/sub_module_analysis.o:$(DIRana)/sub_module_analysis.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIRana)/sub_module_analysis.f90
$(OBJ)/sub_analyse.o:$(DIRana)/sub_analyse.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRana)/sub_analyse.f90
$(OBJ)/sub_NLO.o:$(DIRana)/sub_NLO.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRana)/sub_NLO.f90
$(OBJ)/sub_VibRot.o:$(DIRana)/sub_VibRot.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRana)/sub_VibRot.f90
$(OBJ)/sub_intensity.o:$(DIRana)/sub_intensity.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRana)/sub_intensity.f90
#===============================================================================
# sub_Optimization
$(OBJ)/sub_main_Optimization.o:$(DIROpt)/sub_main_Optimization.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIROpt)/sub_main_Optimization.f90
$(OBJ)/sub_module_SimulatedAnnealing.o:$(DIROpt)/sub_module_SimulatedAnnealing.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIROpt)/sub_module_SimulatedAnnealing.f90
$(OBJ)/sub_module_BFGS.o:$(DIROpt)/sub_module_BFGS.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIROpt)/sub_module_BFGS.f90
$(OBJ)/sub_module_Optimization.o:$(DIROpt)/sub_module_Optimization.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIROpt)/sub_module_Optimization.f90
#
#===============================================================================
# sub_data_initialisation
$(OBJ)/ini_data.o:$(DIR1)/ini_data.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIR1)/ini_data.f90
$(OBJ)/sub_namelist.o:$(DIR1)/sub_namelist.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIR1)/sub_namelist.f90
$(OBJ)/nb_harm.o:$(DIR1)/nb_harm.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIR1)/nb_harm.f90
#
#===============================================================================
# mains
$(OBJ)/$(VIBMAIN).o:$(DIRvib)/$(VIBMAIN).f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIRvib)/$(VIBMAIN).f90
#
$(OBJ)/$(GWPMAIN).o:$(DIRvib)/$(GWPMAIN).f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRvib)/$(GWPMAIN).f90
$(OBJ)/$(WORKMAIN).o:$(DirTNUM)/sub_main/$(WORKMAIN).f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirTNUM)/sub_main/$(WORKMAIN).f90
#
$(OBJ)/EVR_Module.o:$(DIRvib)/EVR_Module.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIRvib)/EVR_Module.f90
$(OBJ)/EVR_driver.o:$(DIRvib)/EVR_driver.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIRvib)/EVR_driver.f90

$(OBJ)/vib.o:$(DIRvib)/vib.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) $(CPPSHELL_ARPACK)  -c $(DIRvib)/vib.f90
$(OBJ)/versionEVR-T.o:$(DIRvib)/versionEVR-T.f90
	cd $(OBJ) ; $(F90_FLAGS) $(CPPpre) -c $(DIRvib)/versionEVR-T.f90
$(OBJ)/cart.o:$(DIRvib)/cart.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRvib)/cart.f90
$(OBJ)/sub_main_nDfit.o:$(DIRvib)/sub_main_nDfit.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRvib)/sub_main_nDfit.f90
#
#===============================================================================
# sub_Smolyak_test
$(OBJ)/sub_Smolyak_DInd.o:$(DIRSmolyak)/sub_Smolyak_DInd.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRSmolyak)/sub_Smolyak_DInd.f90
$(OBJ)/sub_Smolyak_RDP.o:$(DIRSmolyak)/sub_Smolyak_RDP.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRSmolyak)/sub_Smolyak_RDP.f90
$(OBJ)/sub_Smolyak_ba.o:$(DIRSmolyak)/sub_Smolyak_ba.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRSmolyak)/sub_Smolyak_ba.f90
$(OBJ)/sub_Smolyak_module.o:$(DIRSmolyak)/sub_Smolyak_module.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRSmolyak)/sub_Smolyak_module.f90
$(OBJ)/sub_Smolyak_test.o:$(DIRSmolyak)/sub_Smolyak_test.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DIRSmolyak)/sub_Smolyak_test.f90
#
#===============================================================================
#
#===============================================================================
# system (molecule) dependent subroutines and fonctions (potential ....)
# sub_system.o
$(OBJ)/sub_system.o:$(DirPot)/sub_system.$(extf)
	sed "s/zmatrix/CoordType/" $(DirPot)/sub_system.$(extf) > $(DirPot)/sub_system.i
	@echo Warning the sub_system.$(extf) file has been modified.
	mv $(DirPot)/sub_system.i $(DirPot)/sub_system.$(extf)
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirPot)/sub_system.$(extf)
$(OBJ)/read_para.o:$(DirPot)/read_para.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirPot)/read_para.f90
#
# when the file sub_system.f or sub_system.f90 are not present, ...
#  ... we need the copies from sub_system_save.f or sub_system_save.f90.
$(DirPot)/sub_system.f:
	cp $(DirPot)/sub_system_save.f $(DirPot)/sub_system.f
$(DirPot)/sub_system.f90:
	cp $(DirPot)/sub_system_save.f90 $(DirPot)/sub_system.f90
#
#===============================================================================
$(OBJ)/sub_diago.o:$(DirMath)/sub_diago.f90
	cd $(OBJ) ; $(F90_FLAGS)  $(CPPpre) $(CPPSHELL_LAPACK) -c $(DirMath)/sub_diago.f90
$(OBJ)/sub_trans_mat.o:$(DirMath)/sub_trans_mat.f90
	cd $(OBJ) ; $(F90_FLAGS)  $(CPPpre) $(CPPSHELL_LAPACK) -c $(DirMath)/sub_trans_mat.f90
$(OBJ)/sub_math_util.o:$(DirMath)/sub_math_util.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirMath)/sub_math_util.f90
$(OBJ)/sub_integration.o:$(DirMath)/sub_integration.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirMath)/sub_integration.f90
$(OBJ)/sub_SparseBG.o:$(DirMath)/sub_SparseBG.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirMath)/sub_SparseBG.f90
$(OBJ)/sub_polyortho.o:$(DirMath)/sub_polyortho.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirMath)/sub_polyortho.f90
$(OBJ)/sub_function.o:$(DirMath)/sub_function.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirMath)/sub_function.f90
$(OBJ)/sub_derive.o:$(DirMath)/sub_derive.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirMath)/sub_derive.f90
$(OBJ)/sub_pert.o:$(DirMath)/sub_pert.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirMath)/sub_pert.f90
$(OBJ)/sub_fft.o:$(DirMath)/sub_fft.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirMath)/sub_fft.f90
$(OBJ)/Wigner3j.o:$(DirSHTOOLS)/Wigner3j.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirSHTOOLS)/Wigner3j.f90
#
$(OBJ)/sub_io.o:$(DirIO)/sub_io.f90
	cd $(OBJ) ; $(F90_FLAGS)   -c $(DirIO)/sub_io.f90

install: all doc
	@echo "Installing documentation in doc"
	@mkdir -p doc/img
	@mkdir -p doc/reference
#doc: $(HTML)
#	@perl $(DOCGEN) $(DOCINDEX) sub_module/sub_module_deftype_TnumTana.f90 $(REFPATH)/def.html # second pass
#	@perl $(DOCINDEXGEN) $(DOCINDEXFILE) $(REFPATH)/lib_index.html
#	@rm -f doc/reference/index.store
doxy:
	@echo "Installing documentation with doxygen"
	@cd doc_ElVibRot-TnumTana ; doxygen PhysConst_doxygen_settings

$(HTML) : $(REFPATH)/%.html : sub_module/%.f90
	@perl $(DOCGEN) $(DOCINDEX) $< $@

#=======================================================================================
#=======================================================================================
#add dependence for parallelization

ifeq ($(parallel_make),1)
  include ./dependency.mk
endif

#=======================================================================================
#=======================================================================================
ifeq ($(F90),mpifort)
$(info ***********************************************************************)
$(info ********** run MPI example: make example)
$(info ******** clean MPI example: make clean_example)
$(info ***********************************************************************)
endif

# test

ifeq ($(OMP),1)
  parall=openMP
  parall_name=_openMP
else
  parall=NaN
  parall_name=_noparall
endif

.PHONY: example
example:
ifeq ($(F90),mpifort)
	@cd ./Working_tests/MPI_tests ; ./MPI_test.sh
else
	@cd ./Working_tests/MPI_tests ; ./openMP_test.sh
endif

# clean test results
.PHONY: clean_example
clean_example:
	@echo "clean MPI examples"
	@rm -rf ./Working_tests/MPI_tests/*/result
	@echo "removed ./Working_tests/MPI_tests/*/result"
ifeq ($(F90),mpifort)
	@rm -rf ./Working_tests/MPI_tests/*/MPI_test.log
	@echo "removed ./Working_tests/MPI_tests/*/MPI_test.log"
endif
ifeq ($(F90),$(filter $(F90), gfortran ifort pgf90))
	@rm -rf ./Working_tests/MPI_tests/*/openMP_test.log
	@echo "removed ./Working_tests/MPI_tests/*/openMP_test.log"
endif
	@echo "MPI examples cleaned"
