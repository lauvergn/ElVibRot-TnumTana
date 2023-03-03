#=================================================================================
#=================================================================================
# Compiler?
#Possible values: (Empty: gfortran)
#                gfortran (version: 9.0 linux and osx)
# F90 = mpifort
 FC = gfortran
#
# Optimize? Empty: default No optimization; 0: No Optimization; 1 Optimzation
OPT = 1
## OpenMP? Empty: default with OpenMP; 0: No OpenMP; 1 with OpenMP
OMP = 1
## Lapack/blas/mkl? Empty: default with Lapack; 0: without Lapack; 1 with Lapack
LAPACK = 1
## force the default integer (without kind) during the compillation.
## default 4: , INT=8 (for kind=8)
INT = 4
## extension for the "sub_system." file. Possible values: f; f90
extf = f
#
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
# setup for mpifort
ifeq ($(FFC),mpifort)
  ## MPI compiled with: gfortran or ifort
  MPICORE := $(shell ompi_info | grep 'Fort compiler:' | awk '{print $3}')
  OOMP = 0
endif
#===============================================================================
#
# Operating system, OS? automatic using uname:
OS :=$(shell uname)

# about EVRT, path, versions ...:
LOC_path:= $(shell pwd)

# Extension for the object directory and the library
ifeq ($(FFC),mpifort)
  extlibwi_obj:=_$(FFC)_$(MPICORE)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)_int$(INT)
else
  extlibwi_obj:=_$(FFC)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)_int$(INT)
endif
extlib_obj:=_$(FFC)_opt$(OOPT)_omp$(OOMP)_lapack$(LLAPACK)

OBJ_DIR = obj/obj$(extlibwi_obj)
$(info ***********OBJ_DIR:            $(OBJ_DIR))
$(shell [ -d $(OBJ_DIR) ] || mkdir -p $(OBJ_DIR))
MOD_DIR=$(OBJ_DIR)
#
# library name
LIBA=libEVR$(extlibwi_obj).a
#=================================================================================
# cpp preprocessing
CPPSHELL_LAPACK  = -D__LAPACK="$(LLAPACK)"

#===============================================================================
#
#===============================================================================
# external lib (QML, AD_dnSVM ...)
ifeq ($(ExtLibDIR),)
  ExtLibDIR := $(LOC_path)/Ext_Lib
endif
QD_DIR            = $(ExtLibDIR)/QDUtilLib
QDMOD_DIR         = $(QD_DIR)/OBJ/obj$(extlib_obj)
QDLIBA            = $(QD_DIR)/libQD$(extlib_obj).a

AD_DIR            = $(ExtLibDIR)/AD_dnSVM
ADMOD_DIR         = $(AD_DIR)/OBJ/obj$(extlib_obj)
ADLIBA            = $(AD_DIR)/libAD_dnSVM$(extlib_obj).a

QML_DIR           = $(ExtLibDIR)/QuantumModelLib
QMLMOD_DIR        = $(QML_DIR)/OBJ/obj$(extlib_obj)
QMLLIBA           = $(QML_DIR)/libQMLib$(extlib_obj).a

FOREVRT_DIR       = $(ExtLibDIR)/FOR_EVRT
FOREVRTMOD_DIR    = $(FOREVRT_DIR)/obj/obj$(extlibwi_obj)
FOREVRTLIBA       = $(FOREVRT_DIR)/libFOR_EVRT$(extlibwi_obj).a

CONSTPHYS_DIR     = $(ExtLibDIR)/ConstPhys
CONSTPHYSMOD_DIR  = $(CONSTPHYS_DIR)/obj/obj$(extlibwi_obj)
CONSTPHYSLIBA     = $(CONSTPHYS_DIR)/libPhysConst$(extlibwi_obj).a

KEO_DIR           = $(ExtLibDIR)/Coord_KEO_PrimOp
KEOMOD_DIR        = $(KEO_DIR)/obj/obj$(extlibwi_obj)
KEOLIBA           = $(KEO_DIR)/libCoord_KEO_PrimOp$(extlibwi_obj).a

EXTLib     = $(KEOLIBA) $(CONSTPHYSLIBA) $(FOREVRTLIBA) $(QMLLIBA) $(ADLIBA) $(QDLIBA)
#===============================================================================
#
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
  FFLAGS += -I$(KEOMOD_DIR) -I$(CONSTPHYSMOD_DIR) -I$(FOREVRTMOD_DIR) -I$(QMLMOD_DIR) -I$(ADMOD_DIR) -I$(QDMOD_DIR)

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
#=================================================================================
#=================================================================================
# ifort compillation v17 v18 with mkl
#=================================================================================
ifeq ($(FFC),ifort)

  # opt management
  ifeq ($(OOPT),1)
      #F90FLAGS = -O -parallel -g -traceback
      FFLAGS = -O  -g -traceback
  else
      FFLAGS = -O0 -check all -g -traceback
  endif

  # where to store the modules
  FFLAGS +=-module $(MOD_DIR)

  # omp management
  ifeq ($(OOMP),1)
    FFLAGS += -qopenmp
  endif

  # where to look the .mod files
  FFLAGS += -I$(KEOMOD_DIR) -I$(CONSTPHYSMOD_DIR) -I$(FOREVRTMOD_DIR) -I$(QMLMOD_DIR) -I$(ADMOD_DIR) -I$(QDMOD_DIR)

  # some cpreprocessing
  FFLAGS += -cpp $(CPPSHELL)
  ifeq ($(OMP),1)
    FFLAGS += -Drun_openMP=1
  endif

  FLIB    = $(EXTLib)
  ifeq ($(LLAPACK),1)
    #FLIB += -mkl -lpthread
    #FLIB += -qmkl -lpthread
    FLIB +=  ${MKLROOT}/lib/libmkl_blas95_ilp64.a ${MKLROOT}/lib/libmkl_lapack95_ilp64.a ${MKLROOT}/lib/libmkl_intel_ilp64.a \
             ${MKLROOT}/lib/libmkl_intel_thread.a ${MKLROOT}/lib/libmkl_core.a -liomp5 -lpthread -lm -ldl
  else
    FLIB += -lpthread
  endif

  FC_VER = $(shell $(F90) --version | head -1 )

endif
#=================================================================================
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
$(info ***********Lapack:           $(LLAPACK))
$(info ***********FFLAGS0:          $(FFLAGS0))
$(info ***********FLIB:             $(FLIB))
$(info ***********ExtLibDIR:        $(ExtLibDIR))
$(info ************************************************************************)
$(info ************************************************************************)

#==========================================
VPATH = Source_ElVibRot/sub_Basis Source_ElVibRot/sub_Basis/sub_Basis_SG4 \
  Source_ElVibRot/sub_Basis/sub_ReducedDensity Source_ElVibRot/sub_Basis/sub_SymAbelian \
  Source_ElVibRot/sub_CRP Source_ElVibRot/sub_GWP Source_ElVibRot/sub_Operator \
  Source_ElVibRot/sub_Optimization Source_ElVibRot/sub_Smolyak_test Source_ElVibRot/sub_WP \
  Source_ElVibRot/sub_active Source_ElVibRot/sub_analysis Source_ElVibRot/sub_data_initialisation Source_ElVibRot/sub_inactive \
  Source_ElVibRot/sub_main Source_ElVibRot/sub_propagation Source_ElVibRot/sub_rotation sub_pot



#SRCFILES=  $(basis_SRCFILES) $(main_SRCFILES) $(EVR-Mod_SRCFILES) 
include ./f90list.mk

OBJ0=${SRCFILES:.f90=.o}
OBJ0 += QMRPACK_lib.o

OBJ=$(addprefix $(OBJ_DIR)/, $(OBJ0))
$(info ************ OBJ: $(OBJ))
#
#===============================================
#============= Several mains ===================
#===============================================
#===============================================
#==============================================
#ElVibRot:
VIBEXE  = vib.exe
VIBMAIN = EVR-T


#make all : EVR
.PHONY: all evr EVR libEVR libevr lib
evr EVR all :obj vib $(VIBEXE)
	@echo "EVR OK"
lib libEVR libevr: $(LIBA)
	@echo $(LIBA) " OK"
#
# vib script
.PHONY: vib
vib:
	@echo "make vib script"
	./scripts/make_vib.sh $(LOC_path) $(FFC)
	chmod a+x vib
#
$(VIBEXE): $(OBJ_DIR)/$(VIBMAIN).o $(OBJ_DIR)/sub_system.o $(LIBA) $(EXTLib)
	$(FFC) $(FFLAGS) -o $(VIBEXE) $(OBJ_DIR)/$(VIBMAIN).o $(OBJ_DIR)/sub_system.o $(LIBA) $(FLIB)
	@echo EVR-T
#===============================================
#============= TESTS ===========================
#===============================================
#===============================================
#============= Library: lib_FOR_EVRT.a  ========
#===============================================
$(LIBA): $(OBJ) $(EXTLib)
	@echo "  LIBA from OBJ files"
	ar -cr $(LIBA) $(OBJ)
	@echo "  done Library: "$(LIBA)
#
#===============================================
#============= make sub_system =================
#=============  with the .f or .f90 extention ==
#===============================================
sub_pot/sub_system.$(extf): sub_pot/sub_system_save.$(extf)
	cp sub_pot/sub_system_save.$(extf) sub_pot/sub_system.$(extf)
#
#===============================================
#============= compilation =====================
#===============================================
$(OBJ_DIR)/%.o: %.f90
	@echo "  compile: " $<
	$(FFC) $(FFLAGS) -o $@ -c $<
$(OBJ_DIR)/%.o: %.f
	@echo "  compile: " $<
	$(FFC) $(FFLAGS) -o $@ -c $<
#===============================================
#================ cleaning =====================
.PHONY: clean cleanall
clean:
	rm -f  $(OBJ_DIR)/*.o
	rm -f *.log 
	rm -f vib.exe
	@echo "  done cleaning"

cleanall : clean clean_extlib
	rm -fr obj/* build
	rm -f lib*.a
	rm -f *.exe
	rm -f TESTS/res* TESTS/*log
	@echo "  done all cleaning"
#===============================================
#================ zip and copy the directory ===
ExtLibSAVEDIR := /Users/lauvergn/git/Ext_Lib
BaseName := EVR
.PHONY: zip
zip: cleanall
	test -d $(ExtLibSAVEDIR) || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	cd $(ExtLibSAVEDIR) ; rm -rf $(BaseName)_devloc
	mkdir $(ExtLibSAVEDIR)/$(BaseName)_devloc
	cp -r * $(ExtLibSAVEDIR)/$(BaseName)_devloc
	cd $(ExtLibSAVEDIR) ; zip -r Save_$(BaseName)_devloc.zip $(BaseName)_devloc
	cd $(ExtLibSAVEDIR) ; rm -rf $(BaseName)_devloc
	@echo "  done zip"
#===============================================
#=== external libraries ========================
# AD_dnSVM + QML Libs ...
#===============================================
#
$(KEOLIBA):
	@test -d $(ExtLibDIR)     || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(KEO_DIR) || (cd $(ExtLibDIR) ; ./get_Coord_KEO_PrimOp.sh $(EXTLIB_TYPE))
	@test -d $(KEO_DIR) || (echo $(KEO_DIR) "does not exist" ; exit 1)
	cd $(KEO_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) ExtLibDIR=$(ExtLibDIR) INT=$(INT)
	@echo "  done " $(KEO_DIR) " in "$(BaseName)
#	ln -s $(KEO_DIR)/sub_pot sub_pot
#
$(CONSTPHYSLIBA):
	@test -d $(ExtLibDIR)     || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(CONSTPHYS_DIR) || (cd $(ExtLibDIR) ; ./get_ConstPhys.sh $(EXTLIB_TYPE))
	@test -d $(CONSTPHYS_DIR) || (echo $(CONSTPHYS_DIR) "does not exist" ; exit 1)
	cd $(CONSTPHYS_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) ExtLibDIR=$(ExtLibDIR) INT=$(INT)
	@echo "  done " $(CONSTPHYS_DIR) " in "$(BaseName)
#
$(FOREVRTLIBA):
	@test -d $(ExtLibDIR)   || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(FOREVRT_DIR) || (cd $(ExtLibDIR) ; ./get_FOR_EVRT.sh $(EXTLIB_TYPE))
	@test -d $(FOREVRT_DIR) || (echo $(FOREVRT_DIR) "does not exist" ; exit 1)
	cd $(FOREVRT_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) ExtLibDIR=$(ExtLibDIR) INT=$(INT)
	@echo "  done " $(FOREVRTLIBA) " in "$(BaseName)
#
$(QMLLIBA):
	@test -d $(ExtLibDIR) || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(QML_DIR)   || (cd $(ExtLibDIR) ; ./get_QML.sh $(EXTLIB_TYPE))
	@test -d $(QML_DIR)   || (echo $(QML_DIR) "does not exist" ; exit 1)
	cd $(QML_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) ExtLibDIR=$(ExtLibDIR)
	@echo "  done " $(QDLIBA) " in "$(BaseName)
#
$(ADLIBA):
	@test -d $(ExtLibDIR) || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(AD_DIR)    || (cd $(ExtLibDIR) ; ./get_dnSVM.sh  $(EXTLIB_TYPE))
	@test -d $(AD_DIR)    || (echo $(AD_DIR) "does not exist" ; exit 1)
	cd $(AD_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) ExtLibDIR=$(ExtLibDIR)
	@echo "  done " $(AD_DIR) " in "$(BaseName)
#
$(QDLIBA):
	@test -d $(ExtLibDIR) || (echo $(ExtLibDIR) "does not exist" ; exit 1)
	@test -d $(QD_DIR)    || (cd $(ExtLibDIR) ; ./get_QDUtilLib.sh $(EXTLIB_TYPE))
	@test -d $(QD_DIR)    || (echo $(QD_DIR) "does not exist" ; exit 1)
	cd $(QD_DIR) ; make lib FC=$(FFC) OPT=$(OOPT) OMP=$(OOMP) LAPACK=$(LLAPACK) ExtLibDIR=$(ExtLibDIR)
	@echo "  done " $(QDLIBA) " in "$(BaseName)
##
.PHONY: clean_extlib
clean_extlib:
	cd $(ExtLibDIR) ; ./cleanlib
#=======================================================================================
#=======================================================================================
#add dependence for parallelization
#$(OBJ):                     $(EXTLib)
#	@echo "OBJ with EXTLib"
$(OBJ) : $(EXTLib)

#$(OBJ_DIR)/$(VIBMAIN).o:    $(LIBA)
#	@echo "done VIBMAIN.o"

#$(OBJ_DIR)/sub_Auto_Basis.o : $(OBJ_DIR)/sub_module_basis.o
#$(OBJ_DIR)/sub_module_basis.o: $(OBJ_DIR)/sub_module_RotBasis.o $(OBJ_DIR)/sub_module_basis_Grid_Param.o \
#    $(OBJ_DIR)/sub_module_Basis_LTO_n.o $(OBJ_DIR)/sub_SymAbelian.o $(OBJ_DIR)/sub_module_param_SGType2.o

#===============================================
mod_auto_basis = $(OBJ_DIR)/sub_Auto_Basis.o
mod_basis_btog_gtob_sgtype4 = $(OBJ_DIR)/sub_module_basis_BtoG_GtoB_SG4.o
mod_basis_btog_gtob_sgtype4_mpi = $(OBJ_DIR)/sub_module_basis_BtoG_GtoB_SG4_MPI.o
mod_basis_rcvec_sgtype4 = $(OBJ_DIR)/sub_module_basis_RCVec_SG4.o
mod_param_rd = $(OBJ_DIR)/sub_module_param_RD.o
mod_symabelian = $(OBJ_DIR)/sub_SymAbelian.o
basismakegrid = $(OBJ_DIR)/sub_module_BasisMakeGrid.o
mod_basis_l_to_n = $(OBJ_DIR)/sub_module_Basis_LTO_n.o
mod_rotbasis_param = $(OBJ_DIR)/sub_module_RotBasis.o
mod_basis = $(OBJ_DIR)/sub_module_basis.o
mod_basis_btog_gtob = $(OBJ_DIR)/sub_module_basis_BtoG_GtoB.o
mod_basis_btog_gtob_mpi = $(OBJ_DIR)/sub_module_basis_BtoG_GtoB_MPI.o
mod_basis_btog_gtob_sgtype2 = $(OBJ_DIR)/sub_module_basis_BtoG_GtoB_SGType2.o
mod_basis_grid_param = $(OBJ_DIR)/sub_module_basis_Grid_Param.o
mod_basis_set_alloc = $(OBJ_DIR)/sub_module_basis_set_alloc.o
mod_param_sgtype2 = $(OBJ_DIR)/sub_module_param_SGType2.o
mod_crp = $(OBJ_DIR)/sub_CRP.o
mod_matop = $(OBJ_DIR)/sub_MatOp.o
mod_oppsi = $(OBJ_DIR)/sub_OpPsi.o
mod_oppsi_mpi = $(OBJ_DIR)/sub_OpPsi_MPI.o
mod_oppsi_sg4 = $(OBJ_DIR)/sub_OpPsi_SG4.o
mod_oppsi_sg4_mpi = $(OBJ_DIR)/sub_OpPsi_SG4_MPI.o
mod_op = $(OBJ_DIR)/sub_module_Op.o
mod_opgrid = $(OBJ_DIR)/sub_module_OpGrid.o
mod_readop = $(OBJ_DIR)/sub_module_ReadOp.o
mod_setop = $(OBJ_DIR)/sub_module_SetOp.o
mod_bfgs = $(OBJ_DIR)/sub_module_BFGS.o
mod_optimization = $(OBJ_DIR)/sub_module_Optimization.o
mod_simulatedannealing = $(OBJ_DIR)/sub_module_SimulatedAnnealing.o
mod_smolyak_dind = $(OBJ_DIR)/sub_Smolyak_DInd.o
mod_smolyak_rdp = $(OBJ_DIR)/sub_Smolyak_RDP.o
mod_smolyak_ba = $(OBJ_DIR)/sub_Smolyak_ba.o
mod_smolyak_test = $(OBJ_DIR)/sub_Smolyak_module.o
mod_psi = $(OBJ_DIR)/mod_psi.o
mod_ana_psi = $(OBJ_DIR)/sub_module_ana_psi.o
mod_ana_psi_mpi = $(OBJ_DIR)/sub_module_ana_psi_MPI.o
mod_param_wp0 = $(OBJ_DIR)/sub_module_param_WP0.o
mod_psi_b_to_g = $(OBJ_DIR)/sub_module_psi_B_TO_G.o
mod_psi_op = $(OBJ_DIR)/sub_module_psi_Op.o
mod_psi_op_mpi = $(OBJ_DIR)/sub_module_psi_Op_MPI.o
mod_psi_io = $(OBJ_DIR)/sub_module_psi_io.o
mod_psi_set_alloc = $(OBJ_DIR)/sub_module_psi_set_alloc.o
mod_type_ana_psi = $(OBJ_DIR)/sub_module_type_ana_psi.o
mod_set_pararph = $(OBJ_DIR)/sub_paraRPH.o
mod_fullanalysis = $(OBJ_DIR)/sub_analyse.o
mod_analysis = $(OBJ_DIR)/sub_module_analysis.o
mod_evr = $(OBJ_DIR)/EVR_Module.o
mod_ndgridfit = $(OBJ_DIR)/sub_main_nDfit.o
constants = $(OBJ_DIR)/constants.o
system = $(OBJ_DIR)/system.o
tdpes = $(OBJ_DIR)/tdPES.o
mod_hmax_mpi = $(OBJ_DIR)/sub_Hmax_MPI.o
mod_fullcontrol = $(OBJ_DIR)/sub_control.o
mod_arpack = $(OBJ_DIR)/sub_module_Arpack.o
mod_davidson = $(OBJ_DIR)/sub_module_Davidson.o
mod_davidson_mpi = $(OBJ_DIR)/sub_module_Davidson_MPI.o
mod_exactfact = $(OBJ_DIR)/sub_module_ExactFact.o
mod_filter = $(OBJ_DIR)/sub_module_Filter.o
mod_field = $(OBJ_DIR)/sub_module_field.o
mod_march = $(OBJ_DIR)/sub_module_propa_march.o
mod_march_mpi = $(OBJ_DIR)/sub_module_propa_march_MPI.o
mod_march_sg4 = $(OBJ_DIR)/sub_module_propa_march_SG4.o
mod_propa = $(OBJ_DIR)/sub_module_propagation.o
mod_propa_mpi = $(OBJ_DIR)/sub_module_propagation_MPI.o
mod_fullpropa = $(OBJ_DIR)/sub_propagation.o
#===============================================
$(OBJ_DIR)/sub_Auto_Basis.o : \
          $(mod_coord_keo) \
          $(mod_primop) \
          $(mod_ndindex) \
          $(mod_basis) \
          $(mod_op) \
          $(mod_mpi) \
          $(mod_system) \
          $(basismakegrid) \
          $(mod_constant) \
          $(mod_analysis) \
          $(mod_propa) \
          $(mod_psi) \
          $(mod_fullanalysis)
$(OBJ_DIR)/sub_module_basis_BtoG_GtoB_SG4.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_basis_set_alloc) \
          $(mod_param_sgtype2) \
          $(mod_basis_rcvec_sgtype4) \
          $(mod_mpi_aux) \
          $(mod_basis_btog_gtob_sgtype4_mpi)
$(OBJ_DIR)/sub_module_basis_BtoG_GtoB_SG4_MPI.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_basis_set_alloc) \
          $(mod_param_sgtype2) \
          $(mod_basis_rcvec_sgtype4) \
          $(mod_mpi_aux) \
          $(mod_mpi)
$(OBJ_DIR)/sub_module_basis_RCVec_SG4.o : \
          $(mod_system)
$(OBJ_DIR)/sub_module_param_RD.o : \
          $(mod_system) \
          $(mod_ndindex)
$(OBJ_DIR)/sub_SymAbelian.o : \
          $(mod_system)
$(OBJ_DIR)/sub_SymAbelian_OF_Basis.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_basis)
$(OBJ_DIR)/sub_basis_El.o : \
          $(mod_system) \
          $(mod_basis)
$(OBJ_DIR)/sub_module_BasisMakeGrid.o : \
          $(mod_system) \
          $(mod_coord_keo) \
          $(mod_basis) \
          $(mod_ndindex) \
          $(mod_dnsvm)
$(OBJ_DIR)/sub_module_Basis_LTO_n.o : \
          $(mod_system)
$(OBJ_DIR)/sub_module_RotBasis.o : \
          $(mod_system) \
          $(mod_ndindex)
$(OBJ_DIR)/sub_module_basis.o : \
          $(mod_ndindex) \
          $(mod_coord_keo) \
          $(mod_rotbasis_param) \
          $(mod_basis_grid_param) \
          $(mod_basis_l_to_n) \
          $(mod_symabelian) \
          $(mod_param_sgtype2) \
          $(mod_basis_set_alloc) \
          $(mod_basis_btog_gtob) \
          $(mod_system) \
          $(mod_dnsvm) \
          $(mod_basis_btog_gtob_sgtype4)
$(OBJ_DIR)/sub_module_basis_BtoG_GtoB.o : \
          $(mod_system) \
          $(mod_basis_set_alloc) \
          $(mod_basis_btog_gtob_sgtype2) \
          $(mod_param_sgtype2) \
          $(mod_basis_btog_gtob_mpi) \
          $(mod_basis_btog_gtob_sgtype4) \
          $(mod_basis_rcvec_sgtype4) \
          $(mod_mpi_aux)
$(OBJ_DIR)/sub_module_basis_BtoG_GtoB_MPI.o : \
          $(mod_system) \
          $(mod_basis_btog_gtob_sgtype4) \
          $(mod_mpi_aux)
$(OBJ_DIR)/sub_module_basis_BtoG_GtoB_SGType2.o : \
          $(mod_system) \
          $(mod_basis_set_alloc) \
          $(mod_module_dind)
$(OBJ_DIR)/sub_module_basis_Grid_Param.o : \
          $(mod_system)
$(OBJ_DIR)/sub_module_basis_set_alloc.o : \
          $(mod_system) \
          $(mod_dnsvm) \
          $(mod_ndindex) \
          $(mod_rotbasis_param) \
          $(mod_basis_grid_param) \
          $(mod_symabelian) \
          $(mod_param_sgtype2) \
          $(mod_basis_l_to_n) \
          $(mod_param_rd) \
          $(mod_mpi)
$(OBJ_DIR)/sub_module_param_SGType2.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_mpi_aux)
$(OBJ_DIR)/sub_quadra_DirProd.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_basis) \
          $(mod_dnsvm)
$(OBJ_DIR)/sub_quadra_SincDVR.o : \
          $(mod_system) \
          $(mod_basis)
$(OBJ_DIR)/sub_quadra_SparseBasis.o : \
          $(mod_system) \
          $(mod_dnsvm) \
          $(mod_ndindex) \
          $(mod_coord_keo) \
          $(mod_primop) \
          $(mod_basis) \
          $(mod_op) \
          $(mod_auto_basis) \
          $(mod_basis_btog_gtob_sgtype4) \
          $(mod_module_dind)
$(OBJ_DIR)/sub_quadra_SparseBasis2n.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_coord_keo) \
          $(mod_basis)
$(OBJ_DIR)/sub_quadra_Wigner.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_basis)
$(OBJ_DIR)/sub_quadra_Ylm.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_basis)
$(OBJ_DIR)/sub_quadra_box.o : \
          $(mod_system) \
          $(mod_basis)
$(OBJ_DIR)/sub_quadra_fourier.o : \
          $(mod_system) \
          $(mod_basis)
$(OBJ_DIR)/sub_quadra_ft.o : \
          $(mod_system) \
          $(mod_basis)
$(OBJ_DIR)/sub_quadra_herm.o : \
          $(mod_system) \
          $(mod_basis) \
          $(basismakegrid) \
          $(mod_ndindex) \
          $(mod_dnsvm)
$(OBJ_DIR)/sub_quadra_inact.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_coord_keo) \
          $(mod_basis)
$(OBJ_DIR)/sub_quadra_laguerre.o : \
          $(mod_system) \
          $(mod_basis)
$(OBJ_DIR)/sub_quadra_legendre.o : \
          $(mod_system) \
          $(mod_basis)
$(OBJ_DIR)/sub_read_data.o : \
          $(mod_system) \
          $(mod_basis) \
          $(mod_coord_keo) \
          $(mod_constant)
$(OBJ_DIR)/sub_CRP.o : \
          $(mod_system) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_basis) \
          $(mod_op) \
          $(mod_psi) \
          $(mod_realwithunit) \
          $(mod_dnsvm) \
          $(mod_ndindex)
$(OBJ_DIR)/sub_calc_crp_P_lanczos-withNAG.o : \
          $(mod_op) \
          $(mod_system) \
          $(mod_psi)
$(OBJ_DIR)/sub_MatOp.o : \
          $(mod_system) \
          $(mod_constant) \
          $(mod_setop) \
          $(mod_psi) \
          $(mod_oppsi) \
          $(mod_ndindex) \
          $(mod_ana_psi) \
          $(mod_basis_btog_gtob_sgtype4)
$(OBJ_DIR)/sub_OpPsi.o : \
          $(mod_oppsi_sg4) \
          $(mod_oppsi_sg4_mpi) \
          $(mod_system) \
          $(mod_psi) \
          $(mod_setop) \
          $(mod_mpi) \
          $(mod_mpi_aux) \
          $(mod_symabelian) \
          $(mod_basis_btog_gtob) \
          $(mod_coord_keo) \
          $(mod_basis) \
          $(mod_basis_btog_gtob_sgtype4) \
          $(mod_opgrid) \
          $(mod_primop) \
          $(mod_param_sgtype2)
$(OBJ_DIR)/sub_OpPsi_MPI.o : \
          $(mod_system) \
          $(mod_psi)
$(OBJ_DIR)/sub_OpPsi_SG4.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_coord_keo) \
          $(mod_basis_set_alloc) \
          $(mod_basis_rcvec_sgtype4) \
          $(mod_basis_btog_gtob_sgtype4) \
          $(mod_psi) \
          $(mod_setop) \
          $(mod_basis) \
          $(mod_symabelian) \
          $(mod_psi_set_alloc) \
          $(mod_mpi_aux) \
          $(mod_primop)
$(OBJ_DIR)/sub_OpPsi_SG4_MPI.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_coord_keo) \
          $(mod_symabelian) \
          $(mod_basis_set_alloc) \
          $(mod_basis_btog_gtob_sgtype4) \
          $(mod_psi) \
          $(mod_setop) \
          $(mod_mpi_aux) \
          $(mod_oppsi_sg4) \
          $(mod_basis_btog_gtob_sgtype4_mpi) \
          $(mod_psi_set_alloc) \
          $(mod_mpi)
$(OBJ_DIR)/sub_lib_Op.o : \
          $(mod_system) \
          $(mod_primop) \
          $(mod_basis) \
          $(mod_setop)
$(OBJ_DIR)/sub_module_Op.o : \
          $(mod_setop) \
          $(mod_readop) \
          $(mod_matop) \
          $(mod_oppsi)
$(OBJ_DIR)/sub_module_OpGrid.o : \
          $(mod_system) \
          $(mod_basis_btog_gtob_sgtype4) \
          $(mod_mpi)
$(OBJ_DIR)/sub_module_ReadOp.o : \
          $(mod_system) \
          $(mod_opgrid) \
          $(mod_primop)
$(OBJ_DIR)/sub_module_SetOp.o : \
          $(mod_system) \
          $(mod_primop) \
          $(mod_basis) \
          $(mod_opgrid) \
          $(mod_readop) \
          $(mod_mpi) \
          $(mod_param_sgtype2) \
          $(mod_psi)
$(OBJ_DIR)/sub_main_Optimization.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_dnsvm) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_primop) \
          $(mod_basis) \
          $(basismakegrid) \
          $(mod_psi) \
          $(mod_propa) \
          $(mod_op) \
          $(mod_analysis) \
          $(mod_auto_basis) \
          $(mod_optimization)
$(OBJ_DIR)/sub_module_BFGS.o : \
          $(mod_system) \
          $(mod_dnsvm) \
          $(mod_coord_keo) \
          $(mod_primop) \
          $(mod_basis) \
          $(mod_op) \
          $(mod_auto_basis)
$(OBJ_DIR)/sub_module_Optimization.o : \
          $(mod_system) \
          $(mod_simulatedannealing) \
          $(mod_bfgs)
$(OBJ_DIR)/sub_module_SimulatedAnnealing.o : \
          $(mod_system) \
          $(mod_coord_keo) \
          $(mod_primop) \
          $(mod_basis) \
          $(mod_op) \
          $(mod_auto_basis)
$(OBJ_DIR)/sub_Smolyak_DInd.o : \
          $(mod_system)
$(OBJ_DIR)/sub_Smolyak_RDP.o : \
          $(mod_system) \
          $(mod_smolyak_ba) \
          $(mod_smolyak_dind)
$(OBJ_DIR)/sub_Smolyak_ba.o : \
          $(mod_smolyak_dind) \
          $(mod_system)
$(OBJ_DIR)/sub_Smolyak_module.o : \
          $(mod_smolyak_dind) \
          $(mod_smolyak_rdp) \
          $(mod_smolyak_ba) \
          $(mod_system)
$(OBJ_DIR)/sub_Smolyak_test.o : \
          $(mod_system) \
          $(mod_smolyak_dind) \
          $(mod_smolyak_rdp) \
          $(mod_smolyak_ba) \
          $(mod_smolyak_test)
$(OBJ_DIR)/mod_psi.o : \
          $(mod_param_wp0) \
          $(mod_type_ana_psi) \
          $(mod_psi_set_alloc) \
          $(mod_ana_psi) \
          $(mod_psi_io) \
          $(mod_psi_b_to_g) \
          $(mod_psi_op) \
          $(mod_ana_psi_mpi) \
          $(mod_system) \
          $(mod_coord_keo) \
          $(mod_basis) \
          $(mod_mpi_aux)
$(OBJ_DIR)/sub_module_ana_psi.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_constant) \
          $(mod_type_ana_psi) \
          $(mod_basis) \
          $(mod_psi_set_alloc) \
          $(mod_psi_b_to_g) \
          $(mod_dnsvm) \
          $(mod_param_sgtype2) \
          $(mod_param_rd) \
          $(mod_basis_btog_gtob_sgtype4) \
          $(mod_mpi_aux)
$(OBJ_DIR)/sub_module_ana_psi_MPI.o : \
          $(mod_system) \
          $(mod_psi_set_alloc) \
          $(mod_basis) \
          $(mod_param_sgtype2) \
          $(mod_ana_psi) \
          $(mod_mpi_aux)
$(OBJ_DIR)/sub_module_param_WP0.o : \
          $(mod_system) \
          $(mod_coord_keo)
$(OBJ_DIR)/sub_module_psi_B_TO_G.o : \
          $(mod_basis) \
          $(mod_system) \
          $(mod_psi_set_alloc)
$(OBJ_DIR)/sub_module_psi_Op.o : \
          $(mod_basis) \
          $(mod_system) \
          $(mod_psi_set_alloc) \
          $(mod_mpi_aux) \
          $(mod_ana_psi)
$(OBJ_DIR)/sub_module_psi_Op_MPI.o : \
          $(mod_basis) \
          $(mod_system) \
          $(mod_psi_set_alloc) \
          $(mod_mpi_aux) \
          $(mod_param_sgtype2) \
          $(mod_basis_btog_gtob_sgtype4) \
          $(mod_psi_op)
$(OBJ_DIR)/sub_module_psi_io.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_basis) \
          $(mod_psi_set_alloc) \
          $(mod_ana_psi) \
          $(mod_psi_op) \
          $(mod_param_wp0)
$(OBJ_DIR)/sub_module_psi_set_alloc.o : \
          $(mod_system) \
          $(mod_basis) \
          $(mod_type_ana_psi) \
          $(mod_mpi_aux) \
          $(mod_mpi)
$(OBJ_DIR)/sub_module_type_ana_psi.o : \
          $(mod_system)
$(OBJ_DIR)/sub_Grid_SG4.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_op)
$(OBJ_DIR)/sub_diago_H.o : \
          $(mod_system) \
          $(mod_constant)
$(OBJ_DIR)/sub_ini_act_harm.o : \
          $(mod_system) \
          $(mod_op) \
          $(mod_primop)
$(OBJ_DIR)/sub_lib_act.o : \
          $(mod_system) \
          $(mod_primop) \
          $(mod_op) \
          $(mod_coord_keo) \
          $(mod_basis)
$(OBJ_DIR)/sub_paraRPH.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_dnsvm) \
          $(mod_constant) \
          $(mod_primop) \
          $(mod_basis)
$(OBJ_DIR)/sub_NLO.o : \
          $(mod_system) \
          $(mod_coord_keo) \
          $(mod_basis) \
          $(mod_op) \
          $(mod_analysis)
$(OBJ_DIR)/sub_VibRot.o : \
          $(mod_system) \
          $(mod_coord_keo) \
          $(mod_primop) \
          $(mod_basis) \
          $(mod_psi) \
          $(mod_op) \
          $(mod_analysis)
$(OBJ_DIR)/sub_analyse.o : \
          $(mod_constant) \
          $(mod_system) \
          $(mod_coord_keo) \
          $(mod_basis) \
          $(mod_param_rd) \
          $(mod_psi) \
          $(mod_op) \
          $(mod_analysis) \
          $(mod_ndindex)
$(OBJ_DIR)/sub_intensity.o : \
          $(mod_system) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_basis) \
          $(mod_op) \
          $(mod_analysis)
$(OBJ_DIR)/sub_module_analysis.o : \
          $(mod_system) \
          $(mod_constant) \
          $(mod_psi) \
          $(mod_crp) \
          $(mod_coord_keo)
$(OBJ_DIR)/ini_data.o : \
          $(mod_system) \
          $(mod_dnsvm) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_primop) \
          $(mod_cap) \
          $(mod_basis) \
          $(mod_set_pararph) \
          $(mod_readop) \
          $(mod_op) \
          $(mod_analysis) \
          $(mod_propa) \
          $(mod_auto_basis) \
          $(mod_mpi_aux)
$(OBJ_DIR)/nb_harm.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_primop) \
          $(mod_basis)
$(OBJ_DIR)/sub_namelist.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_basis) \
          $(mod_op) \
          $(mod_primop) \
          $(mod_cap) \
          $(mod_hstep)
$(OBJ_DIR)/sub_HST_harm.o : \
          $(mod_system) \
          $(mod_dnsvm) \
          $(mod_coord_keo) \
          $(mod_basis) \
          $(mod_op) \
          $(mod_primop) \
          $(mod_constant) \
          $(mod_ndindex)
$(OBJ_DIR)/sub_ana_HS.o : \
          $(mod_system)
$(OBJ_DIR)/sub_changement_de_var.o : \
          $(mod_system) \
          $(mod_basis)
$(OBJ_DIR)/sub_inactive_harmo.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_primop) \
          $(mod_basis)
$(OBJ_DIR)/EVR-T.o : \
          $(mod_system) \
          $(mod_ndgridfit)
$(OBJ_DIR)/EVR_Module.o : \
          $(mod_system) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_primop) \
          $(mod_basis) \
          $(mod_psi) \
          $(mod_propa) \
          $(mod_op) \
          $(mod_analysis)
$(OBJ_DIR)/EVR_driver.o : \
          $(mod_evr) \
          $(mod_fullpropa) \
          $(mod_fullcontrol) \
          $(mod_davidson) \
          $(mod_filter) \
          $(mod_arpack) \
          $(mod_op) \
          $(mod_analysis) \
          $(mod_fullanalysis) \
          $(mod_auto_basis) \
          $(mod_psi) \
          $(mod_mpi_aux)
$(OBJ_DIR)/Gauss_numlH.o : \
          $(mod_system) \
          $(mod_coord_keo) \
          $(mod_poly) \
          $(mod_gwp) \
          $(mod_propa)
$(OBJ_DIR)/cart.o : \
          $(mod_system) \
          $(mod_dnsvm) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_cart)
$(OBJ_DIR)/sub_main_nDfit.o : \
          $(mod_system) \
          $(mod_ndindex) \
          $(mod_dnsvm) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_primop)
$(OBJ_DIR)/vib.o : \
          $(mod_system) \
          $(mod_constant) \
          $(mod_coord_keo) \
          $(mod_primop) \
          $(mod_basis) \
          $(mod_psi) \
          $(mod_propa) \
          $(mod_fullpropa) \
          $(mod_fullcontrol) \
          $(mod_davidson) \
          $(mod_filter) \
          $(mod_arpack) \
          $(mod_crp) \
          $(mod_op) \
          $(mod_analysis) \
          $(mod_fullanalysis) \
          $(mod_auto_basis) \
          $(mod_mpi_aux)
$(OBJ_DIR)/main.o : \
          $(system) \
          $(tdpes)
$(OBJ_DIR)/tdPES.o : \
          $(system) \
          $(constants)
$(OBJ_DIR)/sub_Hmax.o : \
          $(mod_system) \
          $(mod_op) \
          $(mod_psi) \
          $(mod_ana_psi_mpi) \
          $(mod_propa) \
          $(mod_fullpropa) \
          $(mod_davidson) \
          $(mod_hmax_mpi) \
          $(mod_mpi_aux) \
          $(mod_march)
$(OBJ_DIR)/sub_Hmax_MPI.o : \
          $(mod_system) \
          $(mod_setop) \
          $(mod_mpi_aux)
$(OBJ_DIR)/sub_TF_autocorr.o : \
          $(mod_system) \
          $(mod_constant) \
          $(mod_propa)
$(OBJ_DIR)/sub_control.o : \
          $(mod_system) \
          $(mod_constant) \
          $(mod_psi) \
          $(mod_op) \
          $(mod_field) \
          $(mod_propa) \
          $(mod_fullpropa) \
          $(mod_march)
$(OBJ_DIR)/sub_module_Arpack.o : \
          $(mod_constant) \
          $(mod_system) \
          $(mod_psi) \
          $(mod_op) \
          $(mod_propa) \
          $(mod_mpi_aux)
$(OBJ_DIR)/sub_module_Davidson.o : \
          $(mod_constant) \
          $(mod_mpi) \
          $(mod_system) \
          $(mod_psi) \
          $(mod_op) \
          $(mod_propa) \
          $(mod_propa_mpi) \
          $(mod_davidson_mpi) \
          $(mod_ana_psi_mpi) \
          $(mod_psi_op_mpi) \
          $(mod_mpi_aux) \
          $(mod_basis)
$(OBJ_DIR)/sub_module_Davidson_MPI.o : \
          $(mod_constant) \
          $(mod_mpi) \
          $(mod_system) \
          $(mod_ana_psi_mpi) \
          $(mod_psi) \
          $(mod_psi_op_mpi) \
          $(mod_propa) \
          $(mod_mpi_aux)
$(OBJ_DIR)/sub_module_ExactFact.o : \
          $(mod_system) \
          $(mod_basis) \
          $(mod_psi) \
          $(mod_op) \
          $(mod_field)
$(OBJ_DIR)/sub_module_Filter.o : \
          $(mod_system) \
          $(mod_constant) \
          $(mod_psi) \
          $(mod_op) \
          $(mod_propa)
$(OBJ_DIR)/sub_module_LinearSystem.o : \
          $(mod_constant) \
          $(mod_mpi) \
          $(mod_system) \
          $(mod_psi) \
          $(mod_propa) \
          $(mod_op)
$(OBJ_DIR)/sub_module_field.o : \
          $(mod_system)
$(OBJ_DIR)/sub_module_propa_march.o : \
          $(mod_system) \
          $(mod_constant) \
          $(mod_psi) \
          $(mod_propa) \
          $(mod_field) \
          $(mod_march_mpi) \
          $(mod_march_sg4) \
          $(mod_op) \
          $(mod_mpi_aux) \
          $(mod_basis) \
          $(mod_oppsi_sg4)
$(OBJ_DIR)/sub_module_propa_march_MPI.o : \
          $(mod_system) \
          $(mod_psi) \
          $(mod_propa) \
          $(mod_march_sg4) \
          $(mod_op) \
          $(mod_psi_op_mpi) \
          $(mod_psi_set_alloc) \
          $(mod_oppsi_mpi) \
          $(mod_propa_mpi) \
          $(mod_basis_btog_gtob_sgtype4) \
          $(mod_basis_btog_gtob_sgtype4_mpi) \
          $(mod_mpi_aux) \
          $(mod_basis_set_alloc) \
          $(mod_param_sgtype2) \
          $(mod_ndindex) \
          $(mod_symabelian) \
          $(mod_psi_op)
$(OBJ_DIR)/sub_module_propa_march_SG4.o : \
          $(mod_system) \
          $(mod_psi) \
          $(mod_propa) \
          $(mod_ndindex) \
          $(mod_coord_keo) \
          $(mod_basis_set_alloc) \
          $(mod_basis_btog_gtob_sgtype4) \
          $(mod_op) \
          $(mod_oppsi_sg4) \
          $(mod_mpi_aux)
$(OBJ_DIR)/sub_module_propagation.o : \
          $(mod_system) \
          $(mod_constant) \
          $(mod_psi) \
          $(mod_field) \
          $(mod_op) \
          $(mod_exactfact) \
          $(mod_mpi_aux) \
          $(mod_realwithunit) \
          $(mod_type_ana_psi) \
          $(mod_coord_keo)
$(OBJ_DIR)/sub_module_propagation_MPI.o : \
          $(mod_propa) \
          $(mod_mpi_aux) \
          $(mod_system) \
          $(mod_op) \
          $(mod_psi_set_alloc) \
          $(mod_psi_op_mpi)
$(OBJ_DIR)/sub_propagation.o : \
          $(mod_constant) \
          $(mod_mpi) \
          $(mod_system) \
          $(mod_op) \
          $(mod_propa) \
          $(mod_psi) \
          $(mod_field) \
          $(mod_march)

