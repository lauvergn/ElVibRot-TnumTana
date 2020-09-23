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

MODULE mod_propa_MPI

PRIVATE
PUBLIC :: MPI_Bcast_param_Davidson,Calc_AutoCorr_SR_MPI

CONTAINS
!=======================================================================================
!> for boardcast derived types: param_Davidson   
!=======================================================================================
  SUBROUTINE MPI_Bcast_param_Davidson(param_DS)
    USE mod_propa,ONLY:param_Davidson
    USE mod_MPI_aux
    IMPLICIT NONE
    
    TYPE(param_Davidson),intent(in)       :: param_DS
    
!        TYPE(MPI_Datatype)                    :: type_MPI(47)
!        TYPE(MPI_Datatype)                    :: param_Davidson_MPI
    Integer(Kind=MPI_INTEGER_KIND)        :: type_MPI(47)
    Integer(Kind=MPI_INTEGER_KIND)        :: param_Davidson_MPI
    
    Integer(Kind=MPI_ADDRESS_KIND)        :: disp(47)
    Integer(Kind=MPI_ADDRESS_KIND)        :: base
    Integer(Kind=MPI_INTEGER_KIND)        :: block_length(47)
    Integer(Kind=MPI_INTEGER_KIND)        :: n_count(3)
    Integer                               :: ii

#if(run_MPI)

    CAll MPI_GET_ADDRESS(param_DS%num_resetH,            disp(1),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%num_checkS,            disp(2),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%residual_max_nb,       disp(3),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%max_it,                disp(4),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%nb_WP,                 disp(5),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%max_WP,                disp(6),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%num_LowestWP,          disp(7),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%nb_WP0,                disp(8),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%nb_readWP,             disp(9),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%nb_readWP_OF_List,     disp(10),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%save_interal,          disp(11),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%save_max_nb,           disp(12),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%symab,                 disp(13),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%conv_hermitian,        disp(14),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%NewVec_type,           disp(15),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%L_filter,              disp(16),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%Lmax_filter,           disp(17),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%M_filter,              disp(18),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%DeltaM_filter,         disp(19),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%Mmax_filter,           disp(20),MPI_err) 
    
    CAll MPI_GET_ADDRESS(param_DS%one_by_one,            disp(21),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%read_WP,               disp(22),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%read_listWP,           disp(23),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%precond,               disp(24),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%formatted_file_readWP, disp(25),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%save_all,              disp(26),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%formatted_file_WP,     disp(27),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%all_lower_states,      disp(28),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%lower_states,          disp(29),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%project_WP0,           disp(30),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%Hmin_propa,            disp(31),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%Hmax_propa,            disp(32),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%With_Grid,             disp(33),MPI_err) 

    CAll MPI_GET_ADDRESS(param_DS%precond_tol,           disp(34),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%save_max_ene,          disp(35),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%scaled_max_ene,        disp(36),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%thresh_project,        disp(37),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%Max_ene,               disp(38),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%RMS_ene,               disp(39),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%conv_ene,              disp(40),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%RMS_resi,              disp(41),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%conv_resi,             disp(42),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%E0_filter,             disp(43),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%W_filter,              disp(44),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%LambdaMin,             disp(45),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%conv_resi,             disp(46),MPI_err) 
    CAll MPI_GET_ADDRESS(param_DS%LambdaMax,             disp(47),MPI_err) 
    
    n_count(1)=20
    n_count(2)=33
    n_count(3)=47
    
    base=disp(1)
    DO ii=1,n_count(3)
      disp(ii)=disp(ii)-base
      block_length(ii)=1
    ENDDO 
    
    DO ii=1,n_count(1)
      type_MPI(ii)=Int_MPI
    ENDDO

    DO ii=n_count(1)+1,n_count(2)
      type_MPI(ii)=MPI_Logical
    ENDDO
    
    DO ii=n_count(2)+1,n_count(3)
      type_MPI(ii)=Real_MPI
    ENDDO

    CALL MPI_TYPE_CREATE_STRUCT(n_count(3),block_length,disp,type_MPI,                 &
                                param_Davidson_MPI,MPI_err) 
    CALL MPI_TYPE_COMMIT(param_Davidson_MPI,MPI_err) 
    CALL MPI_Bcast(param_DS,size1_MPI,param_Davidson_MPI,root_MPI,                     &
                   MPI_COMM_WORLD,MPI_err) 

#endif
  ENDSUBROUTINE MPI_Bcast_param_Davidson
!=======================================================================================


!=======================================================================================
!< calculate auto-correcetion function on Smolyak rep. 
!=======================================================================================
  FUNCTION Calc_AutoCorr_SR_MPI(psi0,psi,para_propa,TT,Write_AC)
    USE mod_system
    USE mod_Op,           ONLY:param_Op
    USE mod_propa,        ONLY:param_propa,Write_AutoCorr
    USE mod_psi_set_alloc,ONLY:param_psi
    USE mod_psi_Op_MPI,   ONLY:Overlap_psi1_psi2_SRB_MPI,Overlap_psi1_psi2_SRG_MPI
    IMPLICIT NONE

    Complex(kind=Rkind)                           :: Calc_AutoCorr_SR_MPI
    TYPE(param_psi),                intent(in)    :: psi0
    TYPE(param_psi),                intent(in)    :: psi
    TYPE(param_propa),              intent(in)    :: para_propa
    Real(kind=Rkind),               intent(in)    :: TT
    Logical,optional,               intent(in)    :: Write_AC

    Complex(kind=Rkind)                           :: cdot
    Logical                                       :: Write_AC_loc

#if(run_MPI)

    IF(present(Write_AC)) THEN
      Write_AC_loc=Write_AC
    ELSE
      Write_AC_loc=.FALSE.
    ENDIF

    IF(psi%SRG_MPI) THEN
      CALL Overlap_psi1_psi2_SRG_MPI(cdot,psi0,psi)
    ELSE IF(psi%SRB_MPI) THEN
      CALL Overlap_psi1_psi2_SRB_MPI(cdot,psi0,psi)
    ENDIF

    IF(Write_AC_loc) THEN
      CALL Write_AutoCorr(para_propa%file_autocorr%unit,TT,cdot)
    ENDIF

    Calc_AutoCorr_SR_MPI=cdot

#endif
  ENDFUNCTION Calc_AutoCorr_SR_MPI
!=======================================================================================

ENDMODULE mod_propa_MPI
