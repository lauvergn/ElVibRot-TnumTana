#=================================================================================
#=================================================================================
# Compiler?
#Possible values: (Empty: gfortran)
#                gfortran (version: 9.0 linux and osx)
# F90 = mpifort
 FC = gfortran
#
## parallel_make=1 to enable parallel make
## parallel_make=0 for fast debug make, no parallel
parallel_make=0
# Optimize? Empty: default No optimization; 0: No Optimization; 1 Optimzation
OPT = 0
## OpenMP? Empty: default with OpenMP; 0: No OpenMP; 1 with OpenMP
OMP = 1
## Lapack/blas/mkl? Empty: default with Lapack; 0: without Lapack; 1 with Lapack
LAPACK = 1
## force the default integer (without kind) during the compillation.
## default 4: , INT=8 (for kind=8)
INT = 4
#
## Arpack? Empty: default No Arpack; 0: without Arpack; 1 with Arpack
ARPACK = 0
## CERFACS? Empty: default No CERFACS; 0: without CERFACS; 1 with CERFACS
CERFACS = 0
## extension for the "sub_system." file. Possible values: f; f90 or $(EXTFextern)
## if $(EXTFextern) is empty, the default is f
extf = $(EXTFextern)
## how to get external libraries;  "loc" (default): from local zip file, Empty or something else (v0.5): from github
EXTLIB_TYPE = loc
#=================================================================================
#=================================================================================
ifeq ($(FC),)
  FFC      := gfortran
else
  FFC      := $(FC)
endif
ifeq ($(OPT),)
  OOPT      := 1
else
  OOPT      := $(OPT)
endif
ifeq ($(OMP),)
  OOMP      := 1
else
  OOMP      := $(OMP)
endif
ifeq ($(LAPACK),)
  LLAPACK      := 1
else
  LLAPACK      := $(LAPACK)
endif
#===============================================================================
# turn off ARPACK when using pgf90
ifeq ($(FFC),pgf90)
  ARPACK = 0
endif
#===============================================================================
# setup for mpifort
ifeq ($(FFC),mpifort)
  ## MPI compiled with: gfortran or ifort
  MPICORE := $(shell ompi_info | grep 'Fort compiler:' | awk '{print $3}')
  OOMP = 0
  ifeq ($(INT),8)
    ARPACK = 0 ## temp here, disable ARPACK for 64-bit case
  endif
endif
#===============================================================================
#
# Operating system, OS? automatic using uname:
OS :=$(shell uname)

# about EVRT, path, versions ...:
EVRT_path:= $(shell pwd)
TNUM_ver:=$(shell awk '/Tnum/ {print $$3}' $(EVRT_path)/version-EVR-T)
TANA_ver:=$(shell awk '/Tana/ {print $$3}' $(EVRT_path)/version-EVR-T)
EVR_ver:=$(shell awk '/EVR/ {print $$3}' $(EVRT_path)/version-EVR-T)

# Extension for the object directory and the library
ifeq ($(FFC),mpifort)
  ext_obj:=_$(FFC)_$(MPICORE)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)_int$(INT)
else
  ext_obj:=_$(FFC)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)_int$(INT)
endif
extlib_obj:=_$(FFC)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)



OBJ_DIR = obj/obj$(ext_obj)
$(info ***********OBJ_DIR:            $(OBJ_DIR))
$(shell [ -d $(OBJ_DIR) ] || mkdir -p $(OBJ_DIR))
MOD_DIR=$(OBJ_DIR)
#SRC_DIR=SRC
#MAIN_DIR=APP
#TESTS_DIR=Tests
UT_DIR      = $(EVRT_path)/UnitTests

# library name
EVRTLIBA=libEVRT$(ext_obj).a
#=================================================================================
# cpp preprocessing
CPPSHELL = -D__COMPILE_DATE="\"$(shell date +"%a %e %b %Y - %H:%M:%S")\"" \
           -D__COMPILE_HOST="\"$(shell hostname -s)\"" \
           -D__COMPILER="'$(FFC)'" \
           -D__COMPILER_VER="'$(FC_VER)'" \
           -D__EVRTPATH="'$(EVRT_path)'" \
           -D__EVR_VER="'$(EVR_ver)'" \
           -D__TNUM_VER="'$(TNUM_ver)'" \
           -D__TANA_VER="'$(TANA_ver)'"
CPPSHELL_ARPACK  = -D__ARPACK="$(ARPACK)"
CPPSHELL_CERFACS = -D__CERFACS="$(CERFACS)"
CPPSHELL_LAPACK  = -D__LAPACK="$(LLAPACK)"

#===============================================================================
#
#===============================================================================
# external lib (QML, AD_dnSVM ...)
ifeq ($(ExtLibDIR),)
  ExtLibDIR := $(EVRT_path)/Ext_Lib
endif

QML_DIR    = $(ExtLibDIR)/QuantumModelLib
QMLMOD_DIR = $(QML_DIR)/OBJ/obj$(extlib_obj)
QMLLIBA    = $(QML_DIR)/libQMLib$(extlib_obj).a

AD_DIR    = $(ExtLibDIR)/AD_dnSVM
ADMOD_DIR = $(AD_DIR)/OBJ/obj$(extlib_obj)
ADLIBA    = $(AD_DIR)/libAD_dnSVM$(extlib_obj).a

QD_DIR    = $(ExtLibDIR)/QDUtilLib
QDMOD_DIR = $(QD_DIR)/OBJ/obj$(extlib_obj)
QDLIBA    = $(QD_DIR)/libQD$(extlib_obj).a

EXTLib     = $(QMLLIBA) $(ADLIBA) $(QDLIBA)
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
# Arpack library
#===============================================================================
ifeq ($(ARPACK),1)
  # Arpack management with the OS
  ifeq ($(OS),Darwin)    # OSX
    #EXTLib += /Users/chen/Linux/Software/ARPACK/libarpack_MAC.a
    EXTLib += /Users/lauvergn/trav/ARPACK/libarpack_OSX.a
  else                   # Linux
    ifeq ($(F90), mpifort)
      ifeq ($(MPICORE), gfortran)
        EXTLib += /u/achen/Software/ARPACK/libarpack_Linux_gfortran.a
      else ifeq ($(MPICORE), ifort)
        EXTLib += /u/achen/Software/ARPACK/libarpack_Linux_ifort.a
      endif
    else ifeq ($(F90), gfortran)
      EXTLib += /u/achen/Software/ARPACK/libarpack_Linux_gfortran.a
    else ifeq ($(F90), ifort)
      EXTLib += /u/achen/Software/ARPACK/libarpack_Linux_ifort.a
    endif
    #ARPACKLIB=/usr/lib64/libarpack.a
    EXTLib += /userTMP/lauvergn/EVR/ARPACK_DML/libarpack_Linux.a
  endif
endif
#===============================================================================
#
CompC=gcc
#===============================================================================

#===============================================================================
# gfortran (osx and linux)
#ifeq ($(F90),gfortran)
#===============================================================================
ifeq ($(F90),$(filter $(F90),gfortran gfortran-8))

  # opt management
  ifeq ($(OOPT),1)
    FFLAGS = -O5 -g -fbacktrace -funroll-loops -ftree-vectorize -falign-loops=16
  else
    FFLAGS = -Og -g -fbacktrace -fcheck=all -fwhole-file -fcheck=pointer -Wuninitialized -finit-real=nan -finit-integer=nan
  endif

  # integer kind management
  ifeq ($(INT),8)
    FFLAGS += -fdefault-integer-8 -Dint8=1
  endif

  # omp management
  ifeq ($(OOMP),1)
    FFLAGS += -fopenmp
  endif
  FFLAGS0 := $(FFLAGS)


  # where to store the .mod files
  FFLAGS +=-J$(MOD_DIR)

  # where to look the .mod files
  FFLAGS += -I$(QMLMOD_DIR) -I$(ADMOD_DIR) -I$(QDMOD_DIR)

  # some cpreprocessing
  FFLAGS += -cpp $(CPPSHELL)
  ifeq ($(OMP),1)
    FFLAGS += -Drun_openMP=1
  endif

  FLIB   = $(EXTLib)
  # OS management
  ifeq ($(LLAPACK),1)
    ifeq ($(OS),Darwin)    # OSX
      # OSX libs (included lapack+blas)
      FLIB += -framework Accelerate
    else                   # Linux
      # linux libs
      FLIB += -llapack -lblas
      #
      # linux libs with mkl and with openmp
      #FLIB = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread
      # linux libs with mkl and without openmp
      #FLIB = -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential
    endif
  endif

  FC_VER = $(shell $(FFC) --version | head -1 )

endif
#=================================================================================

#===============================================================================
#===============================================================================
$(info ************************************************************************)
$(info ***********OS:               $(OS))
$(info ***********COMPILER:         $(FFC))
$(info ***********OPTIMIZATION:     $(OOPT))
$(info ***********COMPILER VERSION: $(FC_VER))
ifeq ($(FFC),mpifort)
$(info ***********COMPILED with:    $(MPICORE))
endif
$(info ***********OpenMP:           $(OOMP))
$(info ***********Arpack:           $(ARPACK))
$(info ***********CERFACS:          $(CERFACS))
$(info ***********Lapack:           $(LLAPACK))
$(info ***********FFLAGS0:          $(FFLAGS0))
$(info ***********FLIB:             $(FLIB))
$(info ***********subsystem file:   sub_system.$(extf))
$(info ************************************************************************)
$(info ************************************************************************)
$(info ***************** TNUM_ver: $(TNUM_ver))
$(info ***************** TANA_ver: $(TANA_ver))
$(info ****************** EVR_ver: $(EVR_ver))
$(info ************************************************************************)
$(info ************************************************************************)
$(info ************ run UnitTests: make UT)
$(info ********** clean UnitTests: make clean_UT)
$(info ************************************************************************)

#==========================================
#VPATH = Source_Lib/sub_system Source_Lib/sub_nDindex Source_Lib/sub_dnSVM Source_Lib/sub_module \
        Source_Lib/sub_communf90/sub_math Source_Lib/sub_communf90/sub_io \
        Source_PhysicalConstants
VPATH = Source_Lib/sub_system Source_Lib/sub_nDindex Source_Lib/sub_dnSVM Source_Lib/sub_module \
        Source_Lib/sub_communf90/sub_math \
        Source_PhysicalConstants

#Primlib_SRCFILES  = \
  sub_module_NumParameters.f90 sub_module_MPI.f90 \
  sub_module_memory.f90 sub_module_string.f90 \
  sub_module_memory_Pointer.f90 sub_module_memory_NotPointer.f90 \
  sub_module_file.f90 sub_module_RW_MatVec.f90 mod_Frac.f90 \
  sub_module_system.f90 \
  sub_module_MPI_aux.f90
Primlib_SRCFILES  = \
  sub_module_MPI.f90 \
  sub_module_memory.f90 sub_module_memory_Pointer.f90 sub_module_memory_NotPointer.f90 \
  sub_module_system.f90 \
  sub_module_MPI_aux.f90

#math_SRCFILES =\
   sub_diago.f90 sub_trans_mat.f90 sub_math_util.f90 sub_integration.f90 \
   sub_polyortho.f90 sub_function.f90 sub_derive.f90 sub_pert.f90 \
   sub_fft.f90
math_SRCFILES =\
   sub_math_util.f90 sub_integration.f90 \
   sub_polyortho.f90 sub_function.f90 sub_derive.f90 sub_pert.f90 \
   sub_fft.f90

#io_SRCFILES = sub_io.f90
io_SRCFILES =

dnSVM_SRCFILES = \
  sub_module_dnS.f90 sub_module_VecOFdnS.f90 sub_module_MatOFdnS.f90 \
  sub_module_dnV.f90 sub_module_dnM.f90 sub_module_IntVM.f90 \
  sub_module_dnSVM.f90

FiniteDiff_SRCFILES = mod_FiniteDiff.f90

# nDindex, Minimize Only list: OK
# USE mod_mod_nDindex and mod_module_DInd
nDindex_SRCFILES  = sub_module_DInd.f90 sub_module_nDindex.f90

# nDfit, Minimize Only list: OK
nDfit_SRCFILES    = sub_module_nDfit.f90


#============================================================================
#Physical constants
#USE mod_constant
PhysConst_SRCFILES = sub_module_RealWithUnit.f90 sub_module_Atom.f90 sub_module_constant.f90
PhysConstEXE       = PhysConst.exe
PhysConstMAIN      = PhysicalConstants_Main
#============================================================================

SRCFILES= $(Primlib_SRCFILES) $(math_SRCFILES) $(io_SRCFILES) $(dnSVM_SRCFILES) $(FiniteDiff_SRCFILES)  \
          $(nDindex_SRCFILES) $(nDfit_SRCFILES) \
          $(PhysConst_SRCFILES)

OBJ0=${SRCFILES:.f90=.o}
OBJ=$(addprefix $(OBJ_DIR)/, $(OBJ0))
$(info ************ OBJ: $(OBJ))
#
#============================================================================
# Physical Constants
.PHONY: PhysConst
PhysConst: $(PhysConstEXE)
	@echo "Physical Constants OK"
#
$(PhysConstEXE): $(EVRTLIBA) $(OBJ_DIR)/$(PhysConstMAIN).o
	$(FFC) $(FFLAGS) -o $(PhysConstEXE) $(OBJ_DIR)/$(PhysConstMAIN).o $(EVRTLIBA) $(FLIB)
#
.PHONY: UT_PhysConst ut_physconst
UT_PhysConst ut_physconst: $(PhysConstEXE)
	@echo "---------------------------------------"
	@echo "Unitary tests for the PhysConst module"
	cd Examples/exa_PhysicalConstants ; ./run_tests > $(UT_DIR)/res_UT_PhysConst ; $(UT_DIR)/PhysConst.sh $(UT_DIR)/res_UT_PhysConst
	@echo "---------------------------------------"
#===============================================
#===============================================
.PHONY: clean_UT
clean_UT:
	@cd UnitTests ; ./clean
	@echo "UnitTests cleaned"
#===============================================
#============= Library: libEVRT.a  =============
#===============================================
.PHONY: lib
lib: $(EVRTLIBA)

$(EVRTLIBA): $(OBJ)
	ar -cr $(EVRTLIBA) $(OBJ)
	@echo "  done Library: "$(EVRTLIBA)
#
#===============================================
#============= compilation =====================
#===============================================
$(OBJ_DIR)/%.o: %.f90
	@echo "  compile: " $<
	$(FFC) $(FFLAGS) -o $@ -c $<
#===============================================
#================ cleaning =====================
.PHONY: clean cleanall
clean:
	rm -f $(OBJ_DIR)/*/*.o $(OBJ_DIR)/*.o
	rm -f lib*.a
	rm -f *.log 
	rm -f TEST*.x
	@echo "  done cleaning"

cleanall : clean clean_extlib
	rm -fr OBJ/obj* OBJ/*mod build
	rm -f lib*.a
	rm -f *.exe
	rm -f TESTS/res* TESTS/*log
	@echo "  done all cleaning"
#===============================================
#=== external libraries ========================
# AD_dnSVM + QML Lib
#===============================================
#
$(QMLLIBA):
	@test -d $(ExtLibDIR) || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(QML_DIR) || (cd $(ExtLibDIR) ; ./get_QML.sh $(EXTLIB_TYPE))
	@test -d $(QML_DIR) || (echo $(QML_DIR) "does not exist" ; exit 1)
	cd $(QML_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) ExtLibDIR=$(ExtLibDIR)
	@echo "  done " $(QDLIBA) " in QML"
#
$(ADLIBA):
	@test -d $(ExtLibDIR) || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(AD_DIR) || (cd $(ExtLibDIR) ; ./get_dnSVM.sh  $(EXTLIB_TYPE))
	@test -d $(AD_DIR) || (echo $(AD_DIR) "does not exist" ; exit 1)
	cd $(AD_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) ExtLibDIR=$(ExtLibDIR)
	@echo "  done " $(AD_DIR) " in QML"
#
$(QDLIBA):
	@test -d $(ExtLibDIR) || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(QD_DIR) || (cd $(ExtLibDIR) ; ./get_QDUtilLib.sh $(EXTLIB_TYPE))
	@test -d $(QD_DIR) || (echo $(QD_DIR) "does not exist" ; exit 1)
	cd $(QD_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) ExtLibDIR=$(ExtLibDIR)
	@echo "  done " $(QDLIBA) " in QML"
##
.PHONY: clean_extlib
clean_extlib:
	cd $(ExtLibDIR) ; ./cleanlib
#=======================================================================================
#=======================================================================================
#add dependence for parallelization
$(OBJ): $(QMLLIBA) $(ADLIBA) $(QDLIBA)
ifeq ($(parallel_make),1)
  include ./dependency.mk
endif