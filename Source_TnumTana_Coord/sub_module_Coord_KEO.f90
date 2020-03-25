!===========================================================================
!===========================================================================
!This file is part of Tnum-Tana.
!
!    Tnum-Tana is a free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    Tnum-Tana is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ElVibRot.  If not, see <http://www.gnu.org/licenses/>.
!
!    Copyright 2015  David Lauvergnat
!      Tnum is written David Lauvergnat [1]
!      Tana is written by Mamadou Ndong [1] and David Lauvergnat [1]
!         with contributions
!          Emil Lund klinting (coupling with MidasCpp) [3]'
!
![1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France
![3]: Department of Chemistry, Aarhus University, DK-8000 Aarhus C, Denmark
!
!===========================================================================
!===========================================================================
MODULE mod_Coord_KEO


  USE mod_Constant

  USE mod_Lib_QTransfo,    ONLY : Write_dnx
  USE mod_freq,            ONLY : gaussian_width
  USE mod_ActiveTransfo,   ONLY : get_Qact,get_Qact0,Set_AllActive,     &
                                  Qact_TO_Qdyn_FROM_ActiveTransfo,      &
                                  Qdyn_TO_Qact_FROM_ActiveTransfo,      &
                                  Qinact2n_TO_Qact_FROM_ActiveTransfo
  USE mod_RPHTransfo,      ONLY : Type_RPHpara_AT_Qact1,Type_RPHTransfo,&
                                  alloc_array,dealloc_array,            &
                                  alloc_rphpara_at_qact1,switch_rph,    &
                                  write_rphtransfo,set_rphtransfo,      &
                                  write_rphpara_at_qact1
  USE mod_CartesianTransfo, ONLY: calc_dnteckart,calc_dntxdnxin_to_dnxout,&
                                  calc_eckartrot,dnmwx_multiref
  USE mod_export_KEO,      ONLY : export3_MCTDH_T
  USE mod_Tnum,            ONLY : Tnum,CoordType,zmatrix,param_PES_FromTnum,&
                                  Read_CoordType,write_coordtype,           &
                                  coordtype1tocoordtype2,dealloc_coordtype, &
                                  sub_coordtype_to_pararph,                 &
                                  sub_pararph_to_coordtype,                 &
                                  type_var_analysis_of_coordtype,           &
                                  CoordTypeRPH_TO_CoordTypeFlex,            &
                                  Set_OptimizationPara_FROM_CoordType

  USE mod_paramQ,          ONLY : sub_dnFCC_TO_dnFcurvi,sub_QactTOdnx,  &
                                  sub_qacttoqit,sub_qplusdq_to_cart,    &
                               read_RefGeom,sub_QactTOd0x,sub_d0xTOQact,&
                                  Set_paramQ_FOR_optimization,          &
                                  Write_Cartg98, Write_XYZ

  USE mod_dnRho,           ONLY : sub3_dnrho_ana,Write_Rho
  USE mod_dnGG_dng,        ONLY : get_d0GG,get_dng_dnGG
  USE mod_dnDetGG_dnDetg,  ONLY : sub3_dndetgg
  USE mod_f2f2Vep,         ONLY : calc3_f2_f1q_num,calc3_f2_f1q_numtay0qinact2n

  USE mod_Tana_keo,        ONLY : compute_analytical_keo
  USE mod_Tana_Tnum,       ONLY : comparison_g_from_tnum_tana
  USE mod_Tana_Sum_OpnD,   ONLY : sum_opnd,write_op,delete_op,Expand_Sum_OpnD_TO_Sum_OpnD
  IMPLICIT NONE

END MODULE mod_Coord_KEO
