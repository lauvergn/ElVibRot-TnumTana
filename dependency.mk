#mod_MPI
lib_dep_mod_MPI=$(OBJ)/sub_module_string.o $(OBJ)/sub_module_memory_NotPointer.o       \
                $(OBJ)/sub_module_MPI_aux.o
$(lib_dep_mod_MPI):$(OBJ)/sub_module_MPI.o

#mod_memory
lib_dep_mod_memory=$(OBJ)/sub_module_memory_Pointer.o $(OBJ)/sub_module_string.o       \
                   $(OBJ)/sub_module_memory_NotPointer.o
$(lib_dep_mod_memory):$(OBJ)/sub_module_memory.o

# mod_system
lib_dep_mod_system=$(OBJ)/Wigner3j.o $(OBJ)/sub_fft.o $(OBJ)/sub_pert.o                \
                   $(OBJ)/sub_io.o $(OBJ)/sub_derive.o $(OBJ)/sub_function.o           \
                   $(OBJ)/sub_polyortho.o $(OBJ)/sub_integration.o                     \
                   $(OBJ)/sub_trans_mat.o $(OBJ)/sub_diago.o                           \
                   $(OBJ)/ThreeDTransfo.o $(OBJ)/TwoDTransfo.o                         \
                   $(OBJ)/OneDTransfo.o $(OBJ)/QTOXanaTransfo.o $(OBJ)/ZmatTransfo.o   \
                   $(OBJ)/CartesianTransfo.o $(OBJ)/Lib_QTransfo.o                     \
                   $(OBJ)/sub_module_Tana_OpnD.o $(OBJ)/sub_module_Tana_OpEl.o         \
                   $(OBJ)/sub_module_Atom.o $(OBJ)/sub_module_IntVM.o                  \
                   $(OBJ)/sub_module_dnV.o $(OBJ)/sub_module_MatOFdnS.o                \
                   $(OBJ)/sub_module_VecOFdnS.o $(OBJ)/sub_module_dnS.o                \
                   $(OBJ)/sub_module_MPI_aux.o $(OBJ)/sub_module_dnM.o                 \
                   $(OBJ)/sub_module_constant.o $(OBJ)/sub_module_Tana_Op1D.o          \
                   $(OBJ)/sub_module_Tana_SumOpnD.o $(OBJ)/sub_module_Tana_VecSumOpnD.o\
                   $(OBJ)/sub_module_Tana_PiEulerRot.o $(OBJ)/FlexibleTransfo.o        \
                   $(OBJ)/BunchPolyTransfo.o $(OBJ)/HyperSpheTransfo.o                 \
                   $(OBJ)/Rot2CoordTransfo.o $(OBJ)/GeneTransfo.o                      \
                   $(OBJ)/LinearNMTransfo.o $(OBJ)/RectilinearNM_Transfo.o             \
                   $(OBJ)/sub_freq.o $(OBJ)/RPHTransfo.o $(OBJ)/ActiveTransfo.o        \
                   $(OBJ)/Qtransfo.o $(OBJ)/Sub_X_TO_Q_ana.o                           \
                   $(OBJ)/sub_dnDetGG_dnDetg.o $(OBJ)/sub_module_Tana_vec_operations.o \
                   $(OBJ)/sub_module_Tana_op.o $(OBJ)/sub_module_Tana_NumKEO.o         \
                   $(OBJ)/sub_module_SimpleOp.o $(OBJ)/sub_module_OnTheFly_def.o       \
                   $(OBJ)/read_para.o $(OBJ)/sub_module_RotBasis.o                     \
                   $(OBJ)/sub_module_basis_Grid_Param.o                                \
                   $(OBJ)/sub_module_Basis_LTO_n.o $(OBJ)/sub_SymAbelian.o             \
                   $(OBJ)/sub_module_param_SGType2.o $(OBJ)/sub_module_param_RD.o      \
                   $(OBJ)/sub_module_basis_RCVec_SG4.o $(OBJ)/sub_module_poly.o        \
                   $(OBJ)/sub_quadra_box.o $(OBJ)/sub_quadra_ft.o                      \
                   $(OBJ)/sub_quadra_Ylm.o $(OBJ)/sub_quadra_DirProd.o                 \
                   $(OBJ)/sub_quadra_SparseBasis2n.o                                   \
                   $(OBJ)/sub_SymAbelian_OF_Basis.o $(OBJ)/sub_module_type_ana_psi.o   \
                   $(OBJ)/sub_module_param_WP0.o $(OBJ)/sub_module_psi_set_alloc.o     \
                   $(OBJ)/sub_module_OpGrid.o                                          \
                   $(OBJ)/sub_inactive_harmo.o $(OBJ)/sub_changement_de_var.o          \
                   $(OBJ)/sub_ana_HS.o $(OBJ)/sub_diago_H.o $(OBJ)/sub_paraRPH.o       \
                   $(OBJ)/sub_module_field.o $(OBJ)/sub_module_propa_march_SG4.o       \
                   $(OBJ)/sub_TF_autocorr.o $(OBJ)/sub_module_analysis.o               \
                   $(OBJ)/sub_NLO.o $(OBJ)/sub_CRP.o $(OBJ)/sub_VibRot.o               \
                   $(OBJ)/sub_intensity.o $(OBJ)/sub_quadra_SparseBasis.o              \
                   $(OBJ)/sub_module_Optimization.o $(OBJ)/sub_module_BFGS.o           \
                   $(OBJ)/sub_module_SimulatedAnnealing.o $(OBJ)/ini_data.o            \
                   $(OBJ)/sub_namelist.o $(OBJ)/nb_harm.o $(OBJ)/versionEVR-T.o        \
                   $(OBJ)/vib.o $(OBJ)/cart.o $(OBJ)/sub_main_nDfit.o                  \
                   $(OBJ)/sub_main_Optimization.o  $(OBJ)/sub_Smolyak_DInd.o           \
                   $(OBJ)/sub_Smolyak_ba.o $(OBJ)/sub_module_cart.o                    \
                   $(OBJ)/sub_Smolyak_RDP.o $(OBJ)/sub_Smolyak_test.o                  \
                   $(OBJ)/$(VIBMAIN).o $(OBJ)/QMRPACK_lib.o $(OBJ)/EVR_Module.o        \
                   $(OBJ)/sub_math_util.o $(OBJ)/Calc_Tab_dnQflex.o                    \
                   $(OBJ)/sub_module_basis_BtoG_GtoB_SG4_MPI.o $(OBJ)/sub_OpPsi_MPI.o  \
                   $(OBJ)/sub_Hmax_MPI.o $(OBJ)/sub_module_propa_march_MPI.o           \
                   $(OBJ)/sub_module_ana_psi_MPI.o $(OBJ)/mod_CAP.o $(OBJ)/mod_HStep.o \
                   $(OBJ)/sub_quadra_SincDVR.o
$(lib_dep_mod_system):$(OBJ)/sub_module_system.o

#mod_EVR
lib_dep_mod_EVR=$(OBJ)/EVR_driver.o
$(lib_dep_mod_EVR):$(OBJ)/EVR_Module.o

#mod_CRP
lib_dep_mod_CRP=$(OBJ)/versionEVR-T.o $(OBJ)/sub_module_analysis.o $(OBJ)/EVR_driver.o \
                $(OBJ)/QMRPACK_lib.o
$(lib_dep_mod_CRP):$(OBJ)/sub_CRP.o

#mod_Coord_KEO
lib_dep_mod_Coord_KEO=$(OBJ)/sub_Auto_Basis.o $(OBJ)/sub_PrimOp_def.o                  \
                      $(OBJ)/sub_module_basis.o $(OBJ)/sub_quadra_SparseBasis2n.o      \
                      $(OBJ)/cart.o $(OBJ)/nb_harm.o $(OBJ)/sub_main_Optimization.o    \
                      $(OBJ)/sub_main_nDfit.o $(OBJ)/EVR_Module.o
$(lib_dep_mod_Coord_KEO):$(OBJ)/sub_module_Coord_KEO.o

#mod_Constant
lib_dep_mod_Constant=$(OBJ)/sub_analyse.o $(OBJ)/sub_freq.o $(OBJ)/sub_diago_H.o       \
                     $(OBJ)/sub_module_analysis.o $(OBJ)/EVR_Module.o                  \
                     $(OBJ)/sub_CRP.o $(OBJ)/sub_module_Davidson_MPI.o
$(lib_dep_mod_Constant):$(OBJ)/sub_module_constant.o

#mod_NumParameters
lib_dep_mod_NumParameters=$(OBJ)/sub_module_memory.o $(OBJ)/sub_module_MPI.o           \
                          $(OBJ)/sub_module_RealWithUnit.o
$(lib_dep_mod_NumParameters):$(OBJ)/sub_module_NumParameters.o

#mod_Tana_Op1D
lib_dep_mod_Tana_Op1D=$(OBJ)/sub_module_Tana_OpnD.o $(OBJ)/sub_module_Tana_OpnD.o      \
                      $(OBJ)/sub_module_Tana_PiEulerRot.o $(OBJ)/sub_module_Tana_op.o  \
                      $(OBJ)/sub_module_Tana_SumOpnD.o
$(lib_dep_mod_Tana_Op1D):$(OBJ)/sub_module_Tana_Op1D.o

#mod_dnSVM
lib_dep_mod_dnSVM=$(OBJ)/Lib_QTransfo.o $(OBJ)/sub_module_DInd.o                       \
                  $(OBJ)/BunchPolyTransfo.o $(OBJ)/QTOXanaTransfo.o                    \
                  $(OBJ)/OneDTransfo.o $(OBJ)/ThreeDTransfo.o $(OBJ)/TwoDTransfo.o     \
                  $(OBJ)/Rot2CoordTransfo.o $(OBJ)/FlexibleTransfo.o                   \
                  $(OBJ)/GeneTransfo.o $(OBJ)/HyperSpheTransfo.o                       \
                  $(OBJ)/RectilinearNM_Transfo.o $(OBJ)/LinearNMTransfo.o              \
                  $(OBJ)/RPHTransfo.o $(OBJ)/ActiveTransfo.o $(OBJ)/Qtransfo.o         \
                  $(OBJ)/sub_dnDetGG_dnDetg.o $(OBJ)/sub_module_SimpleOp.o             \
                  $(OBJ)/sub_module_cart.o $(OBJ)/sub_math_util.o                      \
                  $(OBJ)/mod_FiniteDiff.o $(OBJ)/Calc_Tab_dnQflex.o
$(lib_dep_mod_dnSVM):$(OBJ)/sub_module_dnSVM.o

#mod_dnM
lib_dep_mod_dnM=$(OBJ)/sub_module_dnSVM.o
$(lib_dep_mod_dnM):$(OBJ)/sub_module_dnM.o

#mod_module_DInd
lib_dep_mod_module_DInd=$(OBJ)/sub_module_nDindex.o
$(lib_dep_mod_module_DInd):$(OBJ)/sub_module_DInd.o

#mod_nDindex
lib_dep_mod_nDindex=$(OBJ)/sub_module_nDfit.o $(OBJ)/sub_module_RotBasis.o             \
                    $(OBJ)/sub_module_param_SGType2.o $(OBJ)/sub_module_param_RD.o     \
                    $(OBJ)/nb_harm.o $(OBJ)/sub_main_Optimization.o                    \
                    $(OBJ)/sub_main_nDfit.o $(OBJ)/sub_module_basis_BtoG_GtoB_SG4_MPI.o
$(lib_dep_mod_nDindex):$(OBJ)/sub_module_nDindex.o

#mod_file
lib_dep_mod_file=$(OBJ)/sub_module_RW_MatVec.o
$(lib_dep_mod_file):$(OBJ)/sub_module_file.o

#mod_string
lib_dep_mod_string=$(OBJ)/sub_module_system.o $(OBJ)/sub_module_RealWithUnit.o         \
                   $(OBJ)/mod_Frac.o $(OBJ)/sub_module_file.o
$(lib_dep_mod_string):$(OBJ)/sub_module_string.o

#mod_RW_MatVec
lib_dep_mod_RW_MatVec=$(OBJ)/sub_module_system.o
$(lib_dep_mod_RW_MatVec):$(OBJ)/sub_module_RW_MatVec.o

#mod_Lib_QTransfo
lib_dep_mod_Lib_QTransfo=$(OBJ)/CartesianTransfo.o $(OBJ)/ZmatTransfo.o                \
                         $(OBJ)/BunchPolyTransfo.o
$(lib_dep_mod_Lib_QTransfo):$(OBJ)/Lib_QTransfo.o

#mod_nDFit
lib_dep_mod_nDFit=$(OBJ)/sub_module_Tnum.o $(OBJ)/sub_PrimOp_def.o                     \
                  $(OBJ)/sub_PrimOp_RPH.o
$(lib_dep_mod_nDFit):$(OBJ)/sub_module_nDfit.o

#mod_Tnum
lib_dep_mod_Tnum=$(OBJ)/calc_f2_f1Q.o $(OBJ)/sub_module_paramQ.o $(OBJ)/sub_dnRho.o    \
                 $(OBJ)/sub_export_KEO.o $(OBJ)/sub_system.o                           \
                 $(OBJ)/sub_module_Tana_NumKEO.o $(OBJ)/sub_module_paramQ.o
$(lib_dep_mod_Tnum):$(OBJ)/sub_module_Tnum.o

#mod_paramQ
lib_dep_mod_paramQ=$(OBJ)/calc_dng_dnGG.o $(OBJ)/sub_module_Tana_Export_KEO.o          \
                   $(OBJ)/sub_system.o
$(lib_dep_mod_paramQ):$(OBJ)/sub_module_paramQ.o

#mod_dnGG_dng
lib_dep_mod_dnGG_dng=$(OBJ)/sub_export_KEO.o
$(lib_dep_mod_dnGG_dng):$(OBJ)/calc_dng_dnGG.o

#mod_Tana_NumKEO
lib_dep_mod_Tana_NumKEO=$(OBJ)/sub_module_Tana_Tnum.o $(OBJ)/sub_module_Tana_keo.o
$(lib_dep_mod_Tana_NumKEO):$(OBJ)/sub_module_Tana_NumKEO.o

#mod_dnGG_dng
lib_dep_mod_dnGG_dng=$(OBJ)/calc_f2_f1Q_num.o
$(lib_dep_mod_dnGG_dng):$(OBJ)/calc_dng_dnGG.o

#mod_export_KEO
lib_dep_mod_export_KEO=$(OBJ)/sub_module_Coord_KEO.o
$(lib_dep_mod_export_KEO):$(OBJ)/sub_export_KEO.o

#mod_f2f2Vep
lib_dep_mod_f2f2Vep=$(OBJ)/sub_module_Tana_Tnum.o
$(lib_dep_mod_f2f2Vep):$(OBJ)/calc_f2_f1Q_num.o

#mod_Tana_Tnum
lib_dep_mod_Tana_Tnum=$(OBJ)/sub_module_Coord_KEO.o
$(lib_dep_mod_Tana_Tnum):$(OBJ)/sub_module_Tana_Tnum.o

#mod_PrimOp_def
lib_dep_mod_PrimOp_def=$(OBJ)/sub_PrimOp.o $(OBJ)/sub_onthefly.o $(OBJ)/sub_PrimOp_RPH.o
$(lib_dep_mod_PrimOp_def):$(OBJ)/sub_PrimOp_def.o

#mod_OTF_def
lib_dep_mod_OTF_def=$(OBJ)/sub_PrimOp_def.o
$(lib_dep_mod_OTF_def):$(OBJ)/sub_module_OnTheFly_def.o

#mod_param_SGType2
lib_dep_mod_param_SGType2=$(OBJ)/sub_module_basis_set_alloc.o
$(lib_dep_mod_param_SGType2):$(OBJ)/sub_module_param_SGType2.o

#mod_OTF
lib_dep_mod_OTF=$(OBJ)/sub_PrimOp.o $(OBJ)/sub_PrimOp_RPH.o
$(lib_dep_mod_OTF):$(OBJ)/sub_onthefly.o

#mod_basis_set_alloc
lib_dep_mod_basis_set_alloc=$(OBJ)/sub_module_basis_BtoG_GtoB_SGType2.o                \
                            $(OBJ)/sub_module_basis_BtoG_GtoB.o                        \
                            $(OBJ)/sub_module_basis_BtoG_GtoB_SG4.o                    \
                            $(OBJ)/sub_module_basis_BtoG_GtoB_SG4_MPI.o
$(lib_dep_mod_basis_set_alloc):$(OBJ)/sub_module_basis_set_alloc.o

#mod_basis_BtoG_GtoB
lib_dep_mod_basis_BtoG_GtoB=$(OBJ)/sub_module_basis.o
$(lib_dep_mod_basis_BtoG_GtoB):$(OBJ)/sub_module_basis_BtoG_GtoB.o

#mod_basis_BtoG_GtoB_SGType2
lib_dep_mod_basis_BtoG_GtoB_SGType2=$(OBJ)/sub_module_basis_BtoG_GtoB.o
$(lib_dep_mod_basis_BtoG_GtoB_SGType2):$(OBJ)/sub_module_basis_BtoG_GtoB_SGType2.o

#mod_basis
lib_dep_mod_basis=$(OBJ)/sub_module_BasisMakeGrid.o $(OBJ)/sub_read_data.o             \
                  $(OBJ)/sub_quadra_inact.o $(OBJ)/sub_basis_El.o                      \
                  $(OBJ)/sub_quadra_legendre.o $(OBJ)/sub_quadra_laguerre.o            \
                  $(OBJ)/sub_quadra_herm.o $(OBJ)/sub_quadra_fourier.o                 \
                  $(OBJ)/sub_quadra_box.o $(OBJ)/sub_quadra_ft.o                       \
                  $(OBJ)/sub_quadra_Ylm.o $(OBJ)/sub_quadra_DirProd.o                  \
                  $(OBJ)/sub_SymAbelian_OF_Basis.o $(OBJ)/sub_module_psi_set_alloc.o   \
                  $(OBJ)/sub_changement_de_var.o $(OBJ)/sub_quadra_SparseBasis2n.o     \
                  $(OBJ)/sub_inactive_harmo.o $(OBJ)/sub_paraRPH.o $(OBJ)/nb_harm.o    \
                  $(OBJ)/sub_module_psi_Op_MPI.o $(OBJ)/sub_quadra_SincDVR.o
$(lib_dep_mod_basis):$(OBJ)/sub_module_basis.o

#mod_poly
lib_dep_mod_poly=$(OBJ)/sub_module_GWP.o
$(lib_dep_mod_poly):$(OBJ)/sub_module_poly.o

#mod_basis_BtoG_GtoB_SGType4
lib_dep_mod_basis_BtoG_GtoB_SGType4=$(OBJ)/sub_module_ComOp.o                          \
                                    $(OBJ)/sub_module_OpGrid.o                         \
                                    $(OBJ)/sub_module_basis_BtoG_GtoB.o
$(lib_dep_mod_basis_BtoG_GtoB_SGType4):$(OBJ)/sub_module_basis_BtoG_GtoB_SG4.o

#mod_basis_BtoG_GtoB_SGType4_MPI
lib_dep_mod_basis_BtoG_GtoB_SGType4_MPI=$(OBJ)/sub_OpPsi_SG4_MPI.o                     \
                                        $(OBJ)/sub_module_basis_BtoG_GtoB_SG4.o        \
                                        $(OBJ)/sub_module_ComOp.o                      \
                                        $(OBJ)/sub_module_OpGrid.o                     \
                                        $(OBJ)/sub_module_basis_BtoG_GtoB.o
$(lib_dep_mod_basis_BtoG_GtoB_SGType4_MPI):$(OBJ)/sub_module_basis_BtoG_GtoB_SG4_MPI.o

#BasisMakeGrid
lib_dep_BasisMakeGrid=$(OBJ)/sub_quadra_herm.o
$(lib_dep_BasisMakeGrid):$(OBJ)/sub_module_BasisMakeGrid.o

#mod_psi_set_alloc
lib_dep_mod_psi_set_alloc=$(OBJ)/sub_module_psi_B_TO_G.o $(OBJ)/sub_module_ana_psi.o   \
                          $(OBJ)/sub_module_SetOp.o $(OBJ)/mod_psi.o                   \
                          $(OBJ)/sub_module_ana_psi_MPI.o                              \
                          $(OBJ)/sub_module_psi_Op_MPI.o
$(lib_dep_mod_psi_set_alloc):$(OBJ)/sub_module_psi_set_alloc.o

#mod_ana_psi
lib_dep_mod_ana_psi=$(OBJ)/sub_module_psi_io.o $(OBJ)/sub_module_psi_Op.o              \
                    $(OBJ)/mod_psi.o $(OBJ)/sub_module_ana_psi_MPI.o
$(lib_dep_mod_ana_psi):$(OBJ)/sub_module_ana_psi.o

#mod_ana_psi_MPI
lib_dep_mod_ana_psi_MPI=$(OBJ)/sub_module_Davidson.o $(OBJ)/sub_module_Davidson_MPI.o
$(lib_dep_mod_ana_psi_MPI):$(OBJ)/sub_module_ana_psi_MPI.o

#mod_psi_Op
lib_dep_mod_psi_Op=$(OBJ)/sub_module_psi_io.o $(OBJ)/sub_OpPsi_SG4.o                   \
                   $(OBJ)/sub_OpPsi_SG4_MPI.o $(OBJ)/EVR_Module.o                      \
                   $(OBJ)/sub_module_psi_Op_MPI.o
$(lib_dep_mod_psi_Op):$(OBJ)/sub_module_psi_Op.o

#mod_psi_Op_MPI
lib_dep_mod_psi_Op_MPI=$(OBJ)/sub_module_propagation_MPI.o
$(lib_dep_mod_psi_Op_MPI):$(OBJ)/sub_module_psi_Op_MPI.o

#mod_OpGrid
lib_dep_mod_OpGrid=$(OBJ)/sub_module_ReadOp.o $(OBJ)/sub_module_SetOp.o
$(lib_dep_mod_OpGrid):$(OBJ)/sub_module_OpGrid.o

#mod_SetOp
lib_dep_mod_SetOp=$(OBJ)/sub_OpPsi_SG4.o $(OBJ)/sub_OpPsi_SG4_MPI.o $(OBJ)/sub_MatOp.o \
                  $(OBJ)/sub_module_Op.o $(OBJ)/sub_lib_Op.o $(OBJ)/sub_Hmax_MPI.o
$(lib_dep_mod_SetOp):$(OBJ)/sub_module_SetOp.o

#mod_OpPsi
lib_dep_mod_OpPsi=$(OBJ)/sub_MatOp.o
$(lib_dep_mod_OpPsi):$(OBJ)/sub_OpPsi.o

#mod_OpPsi_MPI
lib_dep_mod_OpPsi_MPI=$(OBJ)/sub_module_propa_march.o
$(lib_dep_mod_OpPsi_MPI):$(OBJ)/sub_OpPsi_MPI.o

#mod_OpPsi_SG4
lib_dep_mod_OpPsi_SG4=$(OBJ)/sub_OpPsi.o sub_module_basis_BtoG_GtoB_SG4.o              \
                      $(OBJ)/sub_OpPsi_SG4_MPI.o
$(lib_dep_mod_OpPsi_SG4):$(OBJ)/sub_OpPsi_SG4.o

#mod_OpPsi_SG4_MPI
lib_dep_mod_OpPsi_SG4_MPI=$(OBJ)/sub_OpPsi.o
$(lib_dep_mod_OpPsi_SG4_MPI):$(OBJ)/sub_OpPsi_SG4_MPI.o

#mod_Op
lib_dep_mod_Op=$(OBJ)/sub_HST_harm.o $(OBJ)/sub_Grid_SG4.o $(OBJ)/sub_ini_act_harm.o   \
               $(OBJ)/sub_lib_act.o $(OBJ)/sub_module_ExactFact.o                      \
               $(OBJ)/sub_module_Davidson.o $(OBJ)/sub_module_Filter.o                 \
               $(OBJ)/sub_module_propagation.o $(OBJ)/sub_module_Arpack.o              \
               $(OBJ)/sub_Hmax.o $(OBJ)/sub_propagation.o $(OBJ)/sub_NLO.o             \
               $(OBJ)/sub_CRP.o $(OBJ)/sub_VibRot.o $(OBJ)/sub_analyse.o               \
               $(OBJ)/sub_Auto_Basis.o $(OBJ)/sub_intensity.o                          \
               $(OBJ)/sub_quadra_SparseBasis.o $(OBJ)/sub_module_SimulatedAnnealing.o  \
               $(OBJ)/sub_module_BFGS.o $(OBJ)/ini_data.o $(OBJ)/sub_namelist.o        \
               $(OBJ)/QMRPACK_lib.o
$(lib_dep_mod_Op):$(OBJ)/sub_module_Op.o

#mod_psi_io
lib_dep_mod_psi_io=$(OBJ)/sub_analyse.o $(OBJ)/mod_psi.o
$(lib_dep_mod_psi_io):$(OBJ)/sub_module_psi_io.o

#mod_MatOp
lib_dep_mod_MatOp=$(OBJ)/sub_module_Op.o
$(lib_dep_mod_MatOp):$(OBJ)/sub_MatOp.o

#mod_propa
lib_dep_mod_propa=$(OBJ)/sub_module_propa_march.o $(OBJ)/sub_module_propa_march_SG4.o  \
                  $(OBJ)/sub_module_propa_march_SG4.o $(OBJ)/sub_module_propa_march.o  \
                  $(OBJ)/sub_TF_autocorr.o $(OBJ)/vib.o $(OBJ)/sub_module_Filter.o     \
                  $(OBJ)/sub_main_Optimization.o $(OBJ)/sub_module_Davidson.o          \
                  $(OBJ)/sub_module_Arpack.o $(OBJ)/ini_data.o $(OBJ)/EVR_Module.o     \
                  $(OBJ)/sub_module_propagation_MPI.o $(OBJ)/sub_module_Davidson_MPI.o \
                  $(OBJ)/sub_module_propa_march_MPI.o
$(lib_dep_mod_propa):$(OBJ)/sub_module_propagation.o

#mod_propa_MPI
lib_dep_mod_propa_MPI=$(OBJ)/sub_module_Davidson.o $(OBJ)/sub_module_propa_march.o
$(lib_dep_mod_propa_MPI):$(OBJ)/sub_module_propagation_MPI.o

#mod_march_SG4
lib_dep_mod_march_SG4=$(OBJ)/sub_module_propa_march.o                                  \
                      $(OBJ)/sub_module_propa_march_MPI.o
$(lib_dep_mod_march_SG4):$(OBJ)/sub_module_propa_march_SG4.o

#mod_march
lib_dep_mod_march=$(OBJ)/sub_propagation.o
$(lib_dep_mod_march):$(OBJ)/sub_module_propa_march.o

#mod_march_MPI
lib_dep_mod_march_MPI=$(OBJ)/sub_module_propa_march.o
$(lib_dep_mod_march_MPI):$(OBJ)/sub_module_propa_march_MPI.o

#mod_FullPropa
lib_dep_mod_FullPropa=$(OBJ)/sub_control.o $(OBJ)/sub_Hmax.o $(OBJ)/EVR_driver.o
$(lib_dep_mod_FullPropa):$(OBJ)/sub_propagation.o

#mod_Smolyak_DInd
lib_dep_mod_Smolyak_DInd=$(OBJ)/sub_Smolyak_module.o $(OBJ)/sub_Smolyak_ba.o
$(lib_dep_mod_Smolyak_DInd):$(OBJ)/sub_Smolyak_DInd.o

#mod_Atom
lib_dep_mod_Atom=$(OBJ)/sub_module_constant.o
$(lib_dep_mod_Atom):$(OBJ)/sub_module_Atom.o

#mod_Tana_OpEl
lib_dep_mod_Tana_OpEl=$(OBJ)/sub_module_Tana_OpnD.o $(OBJ)/sub_module_Tana_Op1D.o      \
                      $(OBJ)/sub_module_Tana_SumOpnD.o
$(lib_dep_mod_Tana_OpEl):$(OBJ)/sub_module_Tana_OpEl.o

#mod_dnV
lib_dep_mod_dnV=$(OBJ)/sub_module_dnM.o
$(lib_dep_mod_dnV):$(OBJ)/sub_module_dnV.o

#mod_Tana_sum_opnd
lib_dep_mod_Tana_sum_opnd=$(OBJ)/sub_module_Tana_VecSumOpnD.o                          \
                          $(OBJ)/sub_module_Tana_PiEulerRot.o
$(lib_dep_mod_Tana_sum_opnd):$(OBJ)/sub_module_Tana_SumOpnD.o

#mod_Tana_OpnD
lib_dep_mod_Tana_OpnD=$(OBJ)/sub_module_Tana_vec_operations.o                          \
                      $(OBJ)/sub_module_Tana_SumOpnD.o                                 \
                      $(OBJ)/sub_module_Tana_PiEulerRot.o $(OBJ)/sub_module_Tana_op.o
$(lib_dep_mod_Tana_OpnD):$(OBJ)/sub_module_Tana_OpnD.o

#mod_Tana_VecSumOpnD
lib_dep_mod_Tana_VecSumOpnD=$(OBJ)/BunchPolyTransfo.o                                  \
                            $(OBJ)/sub_module_Tana_PiEulerRot.o
$(lib_dep_mod_Tana_VecSumOpnD):$(OBJ)/sub_module_Tana_VecSumOpnD.o

#mod_BunchPolyTransfo
lib_dep_mod_BunchPolyTransfo=$(OBJ)/Qtransfo.o $(OBJ)/sub_module_Tana_vec_operations.o
$(lib_dep_mod_BunchPolyTransfo):$(OBJ)/BunchPolyTransfo.o

#mod_QTransfo
lib_dep_mod_QTransfo=$(OBJ)/sub_module_Tnum.o
$(lib_dep_mod_QTransfo):$(OBJ)/Qtransfo.o

#mod_Tana_vec_operations
lib_dep_mod_Tana_vec_operations=(OBJ)/sub_module_Tana_op.o $(OBJ)/sub_module_Tana_op.o
$(lib_dep_mod_Tana_vec_operations):$(OBJ)/sub_module_Tana_vec_operations.o

#mod_dnRho
lib_dep_mod_dnRho=$(OBJ)/sub_module_Tana_NumKEO.o $(OBJ)/calc_dng_dnGG.o
$(lib_dep_mod_dnRho):$(OBJ)/sub_dnRho.o

#mod_Tana_write_mctdh
lib_dep_mod_Tana_write_mctdh=$(OBJ)/sub_module_Tana_keo.o
$(lib_dep_mod_Tana_write_mctdh):$(OBJ)/sub_module_Tana_Export_KEO.o

#mod_PrimOp
lib_dep_mod_PrimOp=$(OBJ)/sub_inactive_harmo.o $(OBJ)/sub_module_SetOp.o               \
                   $(OBJ)/sub_paraRPH.o $(OBJ)/nb_harm.o $(OBJ)/sub_main_Optimization.o\
                   $(OBJ)/sub_main_nDfit.o $(OBJ)/EVR_Module.o                         \
                   $(OBJ)/sub_module_ReadOp.o
$(lib_dep_mod_PrimOp):$(OBJ)/sub_PrimOp.o

#mod_SimulatedAnnealing
lib_dep_mod_SimulatedAnnealing=$(OBJ)/sub_module_Optimization.o
$(lib_dep_mod_SimulatedAnnealing):$(OBJ)/sub_module_SimulatedAnnealing.o

#mod_FullPropa
lib_dep_mod_FullPropa=$(OBJ)/vib.o
$(lib_dep_mod_FullPropa):$(OBJ)/sub_propagation.o

#mod_FullControl
lib_dep_mod_FullControl=$(OBJ)/vib.o $(OBJ)/EVR_driver.o
$(lib_dep_mod_FullControl):$(OBJ)/sub_control.o

#mod_Smolyak_ba
lib_dep_mod_Smolyak_ba=$(OBJ)/sub_Smolyak_RDP.o
$(lib_dep_mod_Smolyak_ba):$(OBJ)/sub_Smolyak_ba.o

#mod_Smolyak_RDP
lib_dep_mod_Smolyak_RDP=$(OBJ)/sub_Smolyak_module.o
$(lib_dep_mod_Smolyak_RDP):$(OBJ)/sub_Smolyak_RDP.o

#mod_Smolyak_test
lib_dep_mod_Smolyak_test=$(OBJ)/sub_Smolyak_test.o
$(lib_dep_mod_Smolyak_test):$(OBJ)/sub_Smolyak_module.o

#mod_nDGridFit
lib_dep_mod_nDGridFit=$(OBJ)/$(VIBMAIN).o
$(lib_dep_mod_nDGridFit):$(OBJ)/sub_main_nDfit.o

#mod_Auto_Basis
lib_dep_mod_Auto_Basis=$(OBJ)/sub_quadra_SparseBasis.o                                 \
                       $(OBJ)/sub_module_SimulatedAnnealing.o                          \
                       $(OBJ)/sub_module_BFGS.o $(OBJ)/sub_namelist.o                  \
                       $(OBJ)/ini_data.o $(OBJ)/vib.o $(OBJ)/EVR_driver.o
$(lib_dep_mod_Auto_Basis):$(OBJ)/sub_Auto_Basis.o

#mod_Optimization
lib_dep_mod_Optimization=$(OBJ)/sub_main_Optimization.o
$(lib_dep_mod_Optimization):$(OBJ)/sub_module_Optimization.o

#mod_BFGS
lib_dep_mod_BFGS=$(OBJ)/sub_module_Optimization.o
$(lib_dep_mod_BFGS):$(OBJ)/sub_module_BFGS.o

#mod_ExactFact
lib_dep_mod_ExactFact=$(OBJ)/sub_module_propagation.o
$(lib_dep_mod_ExactFact):$(OBJ)/sub_module_ExactFact.o

#mod_psi_B_TO_G
lib_dep_mod_psi_B_TO_G=$(OBJ)/sub_module_ana_psi.o
$(lib_dep_mod_psi_B_TO_G):$(OBJ)/sub_module_psi_B_TO_G.o

#mod_CartesianTransfo
lib_dep_mod_CartesianTransfo=$(OBJ)/Qtransfo.o
$(lib_dep_mod_CartesianTransfo):$(OBJ)/CartesianTransfo.o

#mod_Tana_Sum_OpnD
lib_dep_mod_Tana_Sum_OpnD=$(OBJ)/sub_module_Tana_VecSumOpnD.o
$(lib_dep_mod_Tana_Sum_OpnD):$(OBJ)/sub_module_Tana_SumOpnD.o

#mod_Davidson
lib_dep_mod_Davidson=$(OBJ)/sub_Hmax.o
$(lib_dep_mod_Davidson):$(OBJ)/sub_module_Davidson.o

#mod_Davidson_MPI
lib_dep_mod_Davidson_MPI=$(OBJ)/sub_module_Davidson.o
$(lib_dep_mod_Davidson_MPI):$(OBJ)/sub_module_Davidson_MPI.o

#mod_Hmax_MPI
lib_dep_mod_Hmax_MPI=$(OBJ)/sub_Hmax.o
$(lib_dep_mod_Hmax_MPI):$(OBJ)/sub_Hmax_MPI.o

#mod_param_RD
lib_dep_mod_param_RD=$(OBJ)/sub_module_basis_set_alloc.o
$(lib_dep_mod_param_RD):$(OBJ)/sub_module_param_RD.o

#mod_VecOFdnS
lib_dep_mod_VecOFdnS=$(OBJ)/sub_module_MatOFdnS.o $(OBJ)/sub_module_dnSVM.o
$(lib_dep_mod_VecOFdnS):$(OBJ)/sub_module_VecOFdnS.o

#mod_analysis
lib_dep_mod_analysis=$(OBJ)/sub_analyse.o $(OBJ)/sub_NLO.o $(OBJ)/sub_VibRot.o         \
                     $(OBJ)/sub_intensity.o
$(lib_dep_mod_analysis):$(OBJ)/sub_module_analysis.o

#mod_fullanalysis
lib_dep_mod_fullanalysis=$(OBJ)/sub_Auto_Basis.o
$(lib_dep_mod_fullanalysis):$(OBJ)/sub_analyse.o

#mod_dnS
lib_dep_mod_dnS=$(OBJ)/sub_module_VecOFdnS.o $(OBJ)/sub_module_dnV.o
$(lib_dep_mod_dnS):$(OBJ)/sub_module_dnS.o

#mod_param_WP0
lib_dep_mod_param_WP0=$(OBJ)/mod_psi.o
$(lib_dep_mod_param_WP0):$(OBJ)/sub_module_param_WP0.o

#mod_MatOFdnS
lib_dep_mod_MatOFdnS=$(OBJ)/sub_module_dnSVM.o
$(lib_dep_mod_MatOFdnS):$(OBJ)/sub_module_MatOFdnS.o

#mod_psi
lib_dep_mod_psi=$(OBJ)/sub_module_SetOp.o $(OBJ)/sub_OpPsi_MPI.o                       \
                $(OBJ)/sub_module_propa_march_MPI.o $(OBJ)/sub_module_Davidson_MPI.o
$(lib_dep_mod_psi):$(OBJ)/mod_psi.o

#mod_PrimOp_RPH
lib_dep_mod_PrimOp_RPH=$(OBJ)/sub_PrimOp.o
$(lib_dep_mod_PrimOp_RPH):$(OBJ)/sub_PrimOp_RPH.o

#mod_IntVM
lib_dep_mod_IntVM=$(OBJ)/sub_module_dnSVM.o
$(lib_dep_mod_IntVM):$(OBJ)/sub_module_IntVM.o
