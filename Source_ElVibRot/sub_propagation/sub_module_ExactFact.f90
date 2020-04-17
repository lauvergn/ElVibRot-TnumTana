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
MODULE mod_ExactFact
  !USE mod_system
  !USE mod_Constant
  !USE mod_basis
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: sub_ExactFact_analysis


CONTAINS

!================================================================
!
!     Exact factorization analysis
!
!================================================================
SUBROUTINE sub_ExactFact_analysis(T,psi,ana_psi,para_H,Tmax,deltaT,para_field)
  USE mod_system
  USE mod_basis
  USE mod_psi,    ONLY : param_psi,param_ana_psi
  USE mod_Op
  USE mod_field


  IMPLICIT NONE


  real (kind=Rkind),    intent(in)           :: T       ! time
  real (kind=Rkind),    intent(in)           :: Tmax,deltaT ! Tmax, deltaT: Time step
  TYPE (param_psi),     intent(inout)        :: psi

  TYPE (param_ana_psi), intent(inout)        :: ana_psi


  integer :: iterm_pot
!----- for the operator ----------------------------
  TYPE (param_Op),    intent(in)           :: para_H
  TYPE (param_field), intent(in), optional :: para_field

!- for debuging --------------------------------------------------
  character (len=*), parameter :: name_sub='sub_ExactFact_analysis'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
!-------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*) 'Time',T
    write(out_unitp,*)
   END IF
!-------------------------------------------------------

  ! normally, the psi is known on the basis ( psi%CvecB(:) ) and on the grid ( psi%CvecG(:) )
  IF (.NOT. allocated(psi%CvecG) ) THEN
    write(out_unitp,*) 'ERROR in ',name_sub
    write(out_unitp,*) '  psi%CvecG is not allocated'
    write(out_unitp,*) '  => check the source!!'
    STOP
  END IF

  ! the diabatic potential is in para_H%OpGrid(iterm_pot)%Grid(1:nb_qa,1:nb_be,1:nb_be).
  ! Normally, iterm_pot=1 for the potential, but it is better to use para_H%derive_term_TO_iterm(0,0)
  IF (.NOT. associated(para_H%OpGrid)) THEN
    write(out_unitp,*) 'ERROR in ',name_sub
    write(out_unitp,*) '  para_H%OpGrid is not associated.'
    write(out_unitp,*) '  the operators have to be stored in memory'
    write(out_unitp,*) '  Use direct=2 in &active namelist.'
    write(out_unitp,*) '  => check your data!!'
    STOP
  END IF
  iterm_pot = para_H%derive_term_TO_iterm(0,0)
  IF (.NOT. associated(para_H%OpGrid(iterm_pot)%Grid)) THEN
    write(out_unitp,*) 'ERROR in ',name_sub
    write(out_unitp,*) '  para_H%OpGrid(iterm_pot)%Grid is not associated.'
    write(out_unitp,*) '  iterm_pot',iterm_pot
    write(out_unitp,*) '  => check the source!!'
    STOP
  END IF

  SELECT CASE (ana_psi%ExactFact)
  CASE (1)
    IF (T == ZERO) THEN
      CALL sub_ExactFact_analysis_gV(psi,para_H,Tmax,deltaT)
    END IF
    CALL sub_ExactFact_analysis_option1(T,psi,ana_psi,para_H)
  CASE (2)
    CALL sub_ExactFact_analysis_option2(T,psi,ana_psi,para_H)
    STOP 'EF option 2: not yet'
  CASE Default
    write(out_unitp,*) 'ERROR in ',name_sub
    write(out_unitp,*) '  The option for the Exact Factorization are: 1 or 2'
    write(out_unitp,*) '  The option value (ExactFact) is',ana_psi%ExactFact
    write(out_unitp,*) '  => Check the ExactFact value in the &analyse'
    STOP 'EF : no default'
  END SELECT

!----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'END ',name_sub
    CALL flush_perso(out_unitp)
  END IF
!----------------------------------------------------------

END SUBROUTINE sub_ExactFact_analysis
SUBROUTINE sub_ExactFact_analysis_option2(T,psi,ana_psi,para_H)
  USE mod_system
  USE mod_basis
  USE mod_psi,    ONLY : param_psi,param_ana_psi,dealloc_psi,           &
                         sub_PsiBasisRep_TO_GridRep,                    &
                         sub_PsiGridRep_TO_BasisRep,                    &
                         sub_d0d1d2PsiBasisRep_TO_GridRep
  USE mod_Op
  USE mod_field
  IMPLICIT NONE


  real (kind=Rkind),    intent(in)           :: T       ! time
  TYPE (param_psi),     intent(inout)        :: psi

  TYPE (param_ana_psi), intent(inout)        :: ana_psi

!----- for the operator ----------------------------
  TYPE (param_Op),    intent(in)           :: para_H


!-- working parameters --------------------------------
  integer              :: ie,je,iqe,jqe,iq,iterm_pot

  !real (kind=Rkind)    :: Wrho(psi%nb_qa)
  real (kind=Rkind)    :: grid(psi%nb_act1)
  integer              :: iact1,idyn,nio
  complex (kind=Rkind) :: d0psi(psi%nb_qa,psi%nb_be)
  complex (kind=Rkind) :: d1psi(psi%nb_qa,psi%nb_be,psi%nb_act1)
  complex (kind=Rkind) :: dtpsi(psi%nb_qa,psi%nb_be)
  complex (kind=Rkind) :: d1dtpsi(psi%nb_qa,psi%nb_be,psi%nb_act1)
  complex (kind=Rkind) :: d1dtpsi_2(psi%nb_qa,psi%nb_be,psi%nb_act1)

  complex (kind=Rkind) :: d0psi2(psi%nb_qa,psi%nb_be)


  TYPE (param_psi)     :: dpsi


!- for debuging --------------------------------------------------
  character (len=*), parameter :: name_sub='sub_ExactFact_analysis_option2'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
!-------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*) 'Time',T
    write(out_unitp,*)
   END IF
!-------------------------------------------------------

  IF (T == ZERO) THEN
    CALL file_open2(name_file='EF_parameter_dpsi',iunit=nio,lformatted=.TRUE.,append=.FALSE.)
    write(nio,*) '# T, iq, grid(ia), Psi(ie), dPsi(ie)/dT, dPsi(ie)/dQia, d^2Psi(ie)/dTdQia'
  ELSE
    CALL file_open2(name_file='EF_parameter_dpsi',iunit=nio,lformatted=.TRUE.,append=.TRUE.)
  END IF

  ! no derivative
  d0psi(:,:) = reshape(psi%CvecG,shape=(/ psi%nb_qa,psi%nb_be /))

  ! time derivative
  CALL sub_OpPsi(psi,dpsi,para_H) ! H.psi
  CALL sub_PsiBasisRep_TO_GridRep(dpsi) ! put H.psi on the grid
  dtpsi(:,:) = reshape(dpsi%CvecG,shape=(/ psi%nb_qa,psi%nb_be /))
  dtpsi(:,:) = -EYE*dtpsi(:,:) ! -i H.psi

  ! Q derivatives
  DO iact1=1,psi%nb_act1
    dpsi = psi
    idyn = para_H%mole%liste_QactTOQdyn(iact1)
    CALL sub_d0d1d2PsiBasisRep_TO_GridRep(dpsi,tab_derQdyn=(/ idyn,0 /) ) ! put d./dQ psi on the grid
    d1psi(:,:,iact1) = reshape(dpsi%CvecG,shape=(/ psi%nb_qa,psi%nb_be /))
  END DO

  ! cross Q-time derivatives (we use the time derivative, dtpsi(:,:))
  DO iact1=1,psi%nb_act1
    dpsi%CvecG = reshape(dtpsi,shape=(/ psi%nb_qa*psi%nb_be /))
    CALL sub_PsiGridRep_TO_BasisRep(dpsi) ! put dtpsi on the basis

    idyn = para_H%mole%liste_QactTOQdyn(iact1)
    CALL sub_d0d1d2PsiBasisRep_TO_GridRep(dpsi,tab_derQdyn=(/ idyn,0 /) ) ! put d./dQ psi on the grid
    d1dtpsi(:,:,iact1) = reshape(dpsi%CvecG,shape=(/ psi%nb_qa,psi%nb_be /))
  END DO


  ! print the informations
  DO iq=1,psi%nb_qa
    CALL Rec_Qact(Grid(:),psi%BasisnD,iq,para_H%mole)
    write(nio,*) T,iq,Grid(:),d0psi(iq,:),dtpsi(iq,:),d1psi(iq,:,:),d1dtpsi(iq,:,:)
  END DO
  write(nio,*)


  !---------------------------------------------------------------------
  CALL dealloc_psi(dpsi,delete_all=.TRUE.)
  close(nio)

!----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'END ',name_sub
    CALL flush_perso(out_unitp)
  END IF
!----------------------------------------------------------

END SUBROUTINE sub_ExactFact_analysis_option2
SUBROUTINE sub_ExactFact_analysis_gV(psi,para_H,Tmax,deltaT)
  USE mod_system
  USE mod_basis
  USE mod_psi,    ONLY : param_psi
  USE mod_Op
  USE mod_field
  IMPLICIT NONE

  TYPE (param_psi),     intent(inout)        :: psi

  real (kind=Rkind),    intent(in)           :: Tmax,deltaT ! Tmax, deltaT: Time step

!----- for the operator ----------------------------
  TYPE (param_Op),    intent(in)           :: para_H


!-- working parameters --------------------------------
  integer              :: ie,je,iqe,jqe,iq,iterm_pot
  integer, allocatable :: tab_nq(:)
  real (kind=Rkind)    :: Qact0(psi%nb_act1),d0GG(psi%nb_act1,psi%nb_act1)

  real (kind=Rkind)    :: Wrho
  real (kind=Rkind)    :: grid(psi%nb_act1)
  integer              :: iact1,idyn,nio

!- for debuging --------------------------------------------------
  character (len=*), parameter :: name_sub='sub_ExactFact_analysis_gV'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
!-------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*)
   END IF
!-------------------------------------------------------

    CALL file_open2(name_file='EF_parameter_gV',iunit=nio,lformatted=.TRUE.,append=.FALSE.)

    write(nio,*) 'nb_be',psi%nb_be,' number of electronic surfaces'
    write(nio,*) 'nb_act1',psi%nb_act1,' number of active coordinates'
    write(nio,*) 'nb_qa',psi%nb_qa,' full number of grid points'
    IF (psi%BasisnD%nb_basis == 0) THEN ! not a direct product, just one basis set
      write(nio,*) 'tab_nq(:)',psi%nb_qa,' number of grid points per basis set'
    ELSE
      CALL alloc_NParray(tab_nq,(/ psi%BasisnD%nb_basis /),'tab_nq',name_sub)
      CALL get_tab_nq_OF_Qact(tab_nq,psi%BasisnD)
      write(nio,*) 'tab_nq(:)',tab_nq,' number of grid points per basis set'
      CALL dealloc_NParray(tab_nq,'tab_nq',name_sub)
    END IF
    write(nio,*) 'Tmax,deltaT',Tmax,deltaT

    !write the metric tensor
    IF (associated(para_H%para_Tnum%Gref)) THEN
      d0GG = para_H%para_Tnum%Gref(1:psi%nb_act1,1:psi%nb_act1)
    ELSE
      CALL get_Qact0(Qact0,para_H%mole%ActiveTransfo)
      CALL get_d0GG(Qact0,para_H%para_Tnum,para_H%mole,d0GG,def=.TRUE.)
    END IF
    write(nio,*) ' metric Tensor (nb_act1 x nb_act1)',psi%nb_act1
    CALL Write_Mat(d0GG,nio,psi%nb_act1)
    write(nio,*)

    ! set the grid and the diabatic potential
    ! this is written only for T=0, hence we write ZERO
    iterm_pot = para_H%derive_term_TO_iterm(0,0)
    write(nio,*) '# T, iq, Wrho, grid(ia), DiabPot(je,ie)'
    DO iq=1,psi%nb_qa
      CALL Rec_Qact(Grid(:),psi%BasisnD,iq,para_H%mole)
      Wrho = Rec_WrhonD(psi%BasisnD,iq)
      write(nio,*) ZERO,iq,Wrho,Grid(:),para_H%OpGrid(iterm_pot)%Grid(iq,:,:)
    END DO
    write(nio,*)
    close(nio)

!----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'END ',name_sub
    CALL flush_perso(out_unitp)
  END IF
!----------------------------------------------------------

END SUBROUTINE sub_ExactFact_analysis_gV
SUBROUTINE sub_ExactFact_analysis_option1(T,psi,ana_psi,para_H)
  USE mod_system
  USE mod_basis
  USE mod_psi,    ONLY : param_psi,param_ana_psi,dealloc_psi,           &
                         sub_PsiBasisRep_TO_GridRep,                    &
                         sub_PsiGridRep_TO_BasisRep,                    &
                         sub_d0d1d2PsiBasisRep_TO_GridRep
  USE mod_Op
  USE mod_field
  IMPLICIT NONE


  real (kind=Rkind),    intent(in)           :: T       ! time
  TYPE (param_psi),     intent(inout)        :: psi

  TYPE (param_ana_psi), intent(inout)        :: ana_psi

!----- for the operator ----------------------------
  TYPE (param_Op),    intent(in)           :: para_H


!-- working parameters --------------------------------
  integer              :: ie,je,iqe,jqe,iq,iterm_pot

  !real (kind=Rkind)    :: Wrho(psi%nb_qa)
  real (kind=Rkind)    :: grid(psi%nb_act1)
  integer              :: iact1,idyn,nio
  complex (kind=Rkind) :: d0psi(psi%nb_qa,psi%nb_be)
  complex (kind=Rkind) :: d1psi(psi%nb_qa,psi%nb_be,psi%nb_act1)
  complex (kind=Rkind) :: dtpsi(psi%nb_qa,psi%nb_be)
  complex (kind=Rkind) :: d1dtpsi(psi%nb_qa,psi%nb_be,psi%nb_act1)
  complex (kind=Rkind) :: d1dtpsi_2(psi%nb_qa,psi%nb_be,psi%nb_act1)

  TYPE (param_psi)     :: dpsi


!- for debuging --------------------------------------------------
  character (len=*), parameter :: name_sub='sub_ExactFact_analysis_option1'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
!-------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*) 'Time',T
    write(out_unitp,*)
   END IF
!-------------------------------------------------------

  IF (T == ZERO) THEN
    CALL file_open2(name_file='EF_parameter_dpsi',iunit=nio,lformatted=.TRUE.,append=.FALSE.)
    write(nio,*) '# T, iq, grid(ia), Psi(ie), dPsi(ie)/dT, dPsi(ie)/dQia, d^2Psi(ie)/dTdQia'
  ELSE
    CALL file_open2(name_file='EF_parameter_dpsi',iunit=nio,lformatted=.TRUE.,append=.TRUE.)
  END IF

  ! no derivative
  d0psi(:,:) = reshape(psi%CvecG,shape=(/ psi%nb_qa,psi%nb_be /))

  ! time derivative
  CALL sub_OpPsi(psi,dpsi,para_H) ! H.psi
  CALL sub_PsiBasisRep_TO_GridRep(dpsi) ! put H.psi on the grid
  dtpsi(:,:) = reshape(dpsi%CvecG,shape=(/ psi%nb_qa,psi%nb_be /))
  dtpsi(:,:) = -EYE*dtpsi(:,:) ! -i H.psi

  ! Q derivatives
  DO iact1=1,psi%nb_act1
    dpsi = psi
    idyn = para_H%mole%liste_QactTOQdyn(iact1)
    CALL sub_d0d1d2PsiBasisRep_TO_GridRep(dpsi,tab_derQdyn=(/ idyn,0 /) ) ! put d./dQ psi on the grid
    d1psi(:,:,iact1) = reshape(dpsi%CvecG,shape=(/ psi%nb_qa,psi%nb_be /))
  END DO

  ! cross Q-time derivatives (we use the time derivative, dtpsi(:,:))
  DO iact1=1,psi%nb_act1
    dpsi%CvecG = reshape(dtpsi,shape=(/ psi%nb_qa*psi%nb_be /))
    CALL sub_PsiGridRep_TO_BasisRep(dpsi) ! put dtpsi on the basis

    idyn = para_H%mole%liste_QactTOQdyn(iact1)
    CALL sub_d0d1d2PsiBasisRep_TO_GridRep(dpsi,tab_derQdyn=(/ idyn,0 /) ) ! put d./dQ psi on the grid
    d1dtpsi(:,:,iact1) = reshape(dpsi%CvecG,shape=(/ psi%nb_qa,psi%nb_be /))
  END DO


  ! print the informations
  DO iq=1,psi%nb_qa
    CALL Rec_Qact(Grid(:),psi%BasisnD,iq,para_H%mole)
    write(nio,*) T,iq,Grid(:),d0psi(iq,:),dtpsi(iq,:),d1psi(iq,:,:),d1dtpsi(iq,:,:)
  END DO
  write(nio,*)


  !---------------------------------------------------------------------
  CALL dealloc_psi(dpsi,delete_all=.TRUE.)
  close(nio)

!----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'END ',name_sub
    CALL flush_perso(out_unitp)
  END IF
!----------------------------------------------------------

END SUBROUTINE sub_ExactFact_analysis_option1
SUBROUTINE sub_ExactFact_analysis_v1(T,psi,ana_psi,para_H,Tmax,deltaT,para_field)
  USE mod_system
  USE mod_basis
  USE mod_psi,    ONLY : param_psi,param_ana_psi,dealloc_psi,           &
                         sub_PsiBasisRep_TO_GridRep,                    &
                         sub_PsiGridRep_TO_BasisRep,                    &
                         sub_d0d1d2PsiBasisRep_TO_GridRep
  USE mod_Op
  USE mod_field
  IMPLICIT NONE


  real (kind=Rkind),    intent(in)           :: T       ! time
  real (kind=Rkind),    intent(in)           :: Tmax,deltaT ! Tmax, deltaT: Time step
  TYPE (param_psi),     intent(inout)        :: psi

  TYPE (param_ana_psi), intent(inout)        :: ana_psi

!----- for the operator ----------------------------
  TYPE (param_Op),    intent(in)           :: para_H
  TYPE (param_field), intent(in), optional :: para_field


!-- working parameters --------------------------------
  integer              :: ie,je,iqe,jqe,iq,iterm_pot
  integer, allocatable :: tab_nq(:)
  real (kind=Rkind)    :: Qact0(psi%nb_act1),d0GG(psi%nb_act1,psi%nb_act1)

  real (kind=Rkind)    :: Wrho(psi%nb_qa)
  real (kind=Rkind)    :: grid(psi%nb_act1,psi%nb_qa)
  integer              :: iact1,idyn,nio
  complex (kind=Rkind) :: d0psi(psi%nb_qa,psi%nb_be)
  complex (kind=Rkind) :: d1psi(psi%nb_qa,psi%nb_be,psi%nb_act1)
  complex (kind=Rkind) :: dtpsi(psi%nb_qa,psi%nb_be)
  complex (kind=Rkind) :: d1dtpsi(psi%nb_qa,psi%nb_be,psi%nb_act1)
  complex (kind=Rkind) :: d1dtpsi_2(psi%nb_qa,psi%nb_be,psi%nb_act1)

  TYPE (param_psi)     :: dpsi,ddpsi


!- for debuging --------------------------------------------------
  character (len=*), parameter :: name_sub='sub_ExactFact_analysis_v1'
  logical, parameter :: debug=.FALSE.
  !logical, parameter :: debug=.TRUE.
!-------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*) 'Time',T
    write(out_unitp,*)
   END IF
!-------------------------------------------------------

  ! normally, the psi is known on the basis ( psi%CvecB(:) ) and on the grid ( psi%CvecG(:) )
  IF (.NOT. allocated(psi%CvecG) ) THEN
    write(out_unitp,*) 'ERROR in ',name_sub
    write(out_unitp,*) '  psi%CvecG is not allocated'
    write(out_unitp,*) '  => check the source!!'
    STOP
  END IF

  ! the diabatic potential is in para_H%OpGrid(iterm_pot)%Grid(1:nb_qa,1:nb_be,1:nb_be).
  ! Normally, iterm_pot=1 for the potential, but it is better to use para_H%derive_term_TO_iterm(0,0)
  IF (.NOT. associated(para_H%OpGrid)) THEN
    write(out_unitp,*) 'ERROR in ',name_sub
    write(out_unitp,*) '  para_H%OpGrid is not associated.'
    write(out_unitp,*) '  the operators have to be stored in memory'
    write(out_unitp,*) '  Use direct=2 in &active namelist.'
    write(out_unitp,*) '  => check your data!!'
    STOP
  END IF
  iterm_pot = para_H%derive_term_TO_iterm(0,0)
  IF (.NOT. associated(para_H%OpGrid(iterm_pot)%Grid)) THEN
    write(out_unitp,*) 'ERROR in ',name_sub
    write(out_unitp,*) '  para_H%OpGrid(iterm_pot)%Grid is not associated.'
    write(out_unitp,*) '  iterm_pot',iterm_pot
    write(out_unitp,*) '  => check the source!!'
    STOP
  END IF

  IF (T == ZERO) THEN
    CALL file_open2(name_file='EF_parameter_gV',iunit=nio,lformatted=.TRUE.,append=.FALSE.)

    write(nio,*) 'nb_be',psi%nb_be,' number of electronic surfaces'
    write(nio,*) 'nb_act1',psi%nb_act1,' number of active coordinates'
    write(nio,*) 'nb_qa',psi%nb_qa,' full number of grid points'
    IF (psi%BasisnD%nb_basis == 0) THEN ! not a direct product, just one basis set
      write(nio,*) 'tab_nq(:)',psi%nb_qa,' number of grid points per basis set'
    ELSE
      CALL alloc_NParray(tab_nq,(/ psi%BasisnD%nb_basis /),'tab_nq',name_sub)
      CALL get_tab_nq_OF_Qact(tab_nq,psi%BasisnD)
      write(nio,*) 'tab_nq(:)',tab_nq,' number of grid points per basis set'
      CALL dealloc_NParray(tab_nq,'tab_nq',name_sub)
    END IF
    write(nio,*) 'Tmax,deltaT',Tmax,deltaT

    !write the metric tensor
    IF (associated(para_H%para_Tnum%Gref)) THEN
      d0GG = para_H%para_Tnum%Gref(1:psi%nb_act1,1:psi%nb_act1)
    ELSE
      CALL get_Qact0(Qact0,para_H%mole%ActiveTransfo)
      CALL get_d0GG(Qact0,para_H%para_Tnum,para_H%mole,d0GG,def=.TRUE.)
    END IF
    write(nio,*) ' metric Tensor (nb_act1 x nb_act1)',psi%nb_act1
    CALL Write_Mat(d0GG,nio,psi%nb_act1)
    write(nio,*)

    ! set the grid and the diabatic potential
    write(nio,*) '# T, iq, Wrho, grid(ia), DiabPot(je,ie)'
    DO iq=1,psi%nb_qa
      CALL Rec_Qact(Grid(:,iq),psi%BasisnD,iq,para_H%mole)
      Wrho(iq) = Rec_WrhonD(psi%BasisnD,iq)
      write(nio,*) T,iq,Wrho(iq),Grid(:,iq),para_H%OpGrid(iterm_pot)%Grid(iq,:,:)
    END DO
    write(nio,*)
    close(nio)

    CALL file_open2(name_file='EF_parameter_dpsi',iunit=nio,lformatted=.TRUE.,append=.FALSE.)
    write(nio,*) '# T, iq, grid(ia), Psi(ie), dPsi(ie)/dT, dPsi(ie)/dQia, d^2Psi(ie)/dTdQia'
  ELSE
    CALL file_open2(name_file='EF_parameter_dpsi',iunit=nio,lformatted=.TRUE.,append=.TRUE.)
  END IF


  ! no derivative
  d0psi(:,:) = reshape(psi%CvecG,shape=(/ psi%nb_qa,psi%nb_be /))

  ! time derivative
  CALL sub_OpPsi(psi,dpsi,para_H) ! H.psi
  CALL sub_PsiBasisRep_TO_GridRep(dpsi) ! put H.psi on the grid
  dtpsi(:,:) = reshape(dpsi%CvecG,shape=(/ psi%nb_qa,psi%nb_be /))
  dtpsi(:,:) = -EYE*dtpsi(:,:) ! -i H.psi

  ! Q derivatives
  DO iact1=1,psi%nb_act1
    dpsi = psi
    idyn = para_H%mole%liste_QactTOQdyn(iact1)
    CALL sub_d0d1d2PsiBasisRep_TO_GridRep(dpsi,tab_derQdyn=(/ idyn,0 /) ) ! put d./dQ psi on the grid
    d1psi(:,:,iact1) = reshape(dpsi%CvecG,shape=(/ psi%nb_qa,psi%nb_be /))
  END DO

  ! cross Q-time derivatives (we use the time derivative, dtpsi(:,:))
  DO iact1=1,psi%nb_act1
    dpsi%CvecG = reshape(dtpsi,shape=(/ psi%nb_qa*psi%nb_be /))
    CALL sub_PsiGridRep_TO_BasisRep(dpsi) ! put dtpsi on the basis

    idyn = para_H%mole%liste_QactTOQdyn(iact1)
    CALL sub_d0d1d2PsiBasisRep_TO_GridRep(dpsi,tab_derQdyn=(/ idyn,0 /) ) ! put d./dQ psi on the grid
    d1dtpsi(:,:,iact1) = reshape(dpsi%CvecG,shape=(/ psi%nb_qa,psi%nb_be /))
  END DO


  ! print the informations
  DO iq=1,psi%nb_qa
    CALL Rec_Qact(Grid(:,iq),psi%BasisnD,iq,para_H%mole)
    write(nio,*) T,iq,Grid(:,iq),d0psi(iq,:),dtpsi(iq,:),d1psi(iq,:,:),d1dtpsi(iq,:,:)
  END DO
  write(nio,*)


  !---------------------------------------------------------------------
  CALL dealloc_psi(dpsi,delete_all=.TRUE.)
  close(nio)

!----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'END ',name_sub
    CALL flush_perso(out_unitp)
  END IF
!----------------------------------------------------------

END SUBROUTINE sub_ExactFact_analysis_v1
SUBROUTINE sub_ExactFact_analysis_v0(T,psi,ana_psi,para_H,para_field)
  USE mod_system
  USE mod_basis
  USE mod_psi,    ONLY : param_psi,param_ana_psi

  USE mod_Op,     ONLY : param_Op,sub_PsiOpPsi
  USE mod_field,  ONLY : param_field,sub_dnE
  IMPLICIT NONE


  real (kind=Rkind),    intent(in)           :: T      ! time
  TYPE (param_psi),     intent(inout)        :: psi

  TYPE (param_ana_psi), intent(inout)        :: ana_psi

!----- for the operator ----------------------------
  TYPE (param_Op),    intent(in)           :: para_H
  TYPE (param_field), intent(in), optional :: para_field


!-- working parameters --------------------------------
  integer           :: ie,je,iqe,jqe,iq,iterm_pot
  real (kind=Rkind), allocatable :: chi2(:),chi(:)
  real (kind=Rkind), allocatable :: EF_Scal_Pot(:) ! exact factorisation scalar potential
  real (kind=Rkind), allocatable :: Qact1(:)

!- for debuging --------------------------------------------------
  character (len=*), parameter :: name_sub='sub_ExactFact_analysis_v0'
  logical, parameter :: debug=.FALSE.
! logical, parameter :: debug=.TRUE.
!-------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*) 'Time',T
    write(out_unitp,*)
   END IF
!-------------------------------------------------------

! normally, the psi is known on the basis ( psi%CvecB(:) ) and on the grid ( psi%CvecG(:) )
  IF (.NOT. allocated(psi%CvecG) ) THEN
    write(out_unitp,*) 'ERROR in ',name_sub
    write(out_unitp,*) '  psi%CvecG is not allocated'
    write(out_unitp,*) '  => check the source!!'
    STOP
  END IF

  ! the diabatic potential is in para_H%OpGrid(iterm_pot)%Grid(1:nb_qa,1:nb_be,1:nb_be).
  ! Normally, iterm_pot=1 for the potential, but it is better to use para_H%derive_term_TO_iterm(0,0)
  IF (.NOT. associated(para_H%OpGrid)) THEN
    write(out_unitp,*) 'ERROR in ',name_sub
    write(out_unitp,*) '  para_H%OpGrid is not associated.'
    write(out_unitp,*) '  the operators have to be stored in memory'
    write(out_unitp,*) '  Use direct=2 in &active namelist.'
    write(out_unitp,*) '  => check your data!!'
    STOP
  END IF
  iterm_pot = para_H%derive_term_TO_iterm(0,0)
  IF (.NOT. associated(para_H%OpGrid(iterm_pot)%Grid)) THEN
    write(out_unitp,*) 'ERROR in ',name_sub
    write(out_unitp,*) '  para_H%OpGrid(iterm_pot)%Grid is not associated.'
    write(out_unitp,*) '  iterm_pot',iterm_pot
    write(out_unitp,*) '  => check the source!!'
    STOP
  END IF


  !---------------------------------------------------------------------
  ! allocation of chi2(:) and chi(:)
  CALL alloc_NParray(chi2, (/ psi%nb_qa /),'chi2',name_sub)
  CALL alloc_NParray(chi,  (/ psi%nb_qa /),'chi' ,name_sub)

  ! calculation of chi2(:) and chi(:)
  iqe = 0
  chi2(:) = ZERO
  DO ie=1,psi%nb_be
    chi2(:) =  chi2(:) + abs(psi%CvecG(iqe+1:iqe+psi%nb_qa))**2
    iqe = iqe+psi%nb_qa
  END DO
  chi(:) = sqrt(chi2(:))
  !---------------------------------------------------------------------

  !---------------------------------------------------------------------
  ! calculation of EF_Scal_Pot(:)
  !---------------------------------------------------------------------

  ! allocation of EF_Scal_Pot(:)
  CALL alloc_NParray(EF_Scal_Pot,(/ psi%nb_qa /),'EF_Scal_Pot' ,name_sub)

  iqe = 0
  EF_Scal_Pot(:) = ZERO
  DO ie=1,psi%nb_be
    jqe = 0
    DO je=1,psi%nb_be

      EF_Scal_Pot(:) =  EF_Scal_Pot(:) +                                                  &
         real(psi%CvecG(iqe+1:iqe+psi%nb_qa)*psi%CvecG(jqe+1:jqe+psi%nb_qa),kind=Rkind) * &
         para_H%OpGrid(iterm_pot)%Grid(:,ie,je)

      jqe = jqe+psi%nb_qa
    END DO
    iqe = iqe+psi%nb_qa
  END DO
  !EF_Scal_Pot(:) =  EF_Scal_Pot(:)/chi2(:)

  ! Write the EF_Scal_Pot(:)
  CALL alloc_NParray(Qact1,(/ para_H%mole%nb_act1 /),'Qact1',name_sub)
  DO iq=1,psi%nb_qa
    CALL Rec_Qact(Qact1,psi%BasisnD,iq,para_H%mole)
    write(out_unitp,*) 'EF_Scal_Pot',T,Qact1,chi2(iq),                          &
           (para_H%OpGrid(iterm_pot)%Grid(iq,ie,ie),ie=1,psi%nb_be),    &
           EF_Scal_Pot(iq),EF_Scal_Pot(iq)/chi2(iq)
  END DO
  write(out_unitp,*)
  CALL dealloc_NParray(Qact1,'Qact1',name_sub)
  !---------------------------------------------------------------------


  !---------------------------------------------------------------------
  ! deallocation of ....
  CALL dealloc_NParray(chi2,'chi2',name_sub)
  CALL dealloc_NParray(chi, 'chi' ,name_sub)
  CALL dealloc_NParray(EF_Scal_Pot,'EF_Scal_Pot' ,name_sub)
  !---------------------------------------------------------------------

!----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'END ',name_sub
  END IF
!----------------------------------------------------------
END SUBROUTINE sub_ExactFact_analysis_v0


END MODULE mod_ExactFact

