MODULE Module_ForTnumTana_Driver
  USE mod_system,                ONLY : Rkind,out_unitp,print_level
  USE mod_Constant,              ONLY : constant
  USE mod_Coord_KEO,             ONLY : zmatrix,Tnum,Read_mole,read_RefGeom,sub_QactTOd0x,sub_d0xTOQact
  USE mod_PrimOp,                ONLY : param_PES,Finalyze_TnumTana_Coord_PrimOp
  IMPLICIT NONE

  TYPE (constant)  :: const_phys
  TYPE (zmatrix)   :: mole
  TYPE (Tnum)      :: para_Tnum
  TYPE (param_PES) :: para_PES

  integer          :: Init = 0 ! Initialization is not done

END MODULE Module_ForTnumTana_Driver
