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
MODULE mod_Set_paraRPH
IMPLICIT NONE

PRIVATE
PUBLIC Set_paraPRH
CONTAINS
      SUBROUTINE Set_paraPRH(mole,para_Tnum,BasisnD)
      USE mod_system
      USE mod_nDindex
      USE mod_dnSVM
      USE mod_Constant
      USE mod_PrimOp
      USE mod_basis
      IMPLICIT NONE

!
!=====================================================================
!
!     variables
!
!=====================================================================

!----- variables for the construction of H ---------------------------
      TYPE (basis)   :: BasisnD
      TYPE (CoordType) :: mole
      TYPE (Tnum)    :: para_Tnum


!------ working variables ---------------------------------
      integer :: ib,nq_part,iq,iq_list,nb_act1_RPH,nb_inact21_RPH,it,comp,iq_list_small


      real (kind=Rkind) :: Qdyn(mole%nb_var)
      real (kind=Rkind) :: Qact(mole%nb_var)

      logical :: Find_in_List,tab_skip_transfo(mole%nb_Qtransfo),RPHCoord_IN_OneBasis



!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Set_paraPRH'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
      IF (mole%tab_Qtransfo(mole%itRPH)%skip_transfo) RETURN


      nb_act1_RPH    = mole%RPHTransfo%nb_act1
      nb_inact21_RPH = mole%RPHTransfo%nb_inact21

      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub

        write(out_unitp,*) 'nb_qa',get_nq_FROM_basis(BasisnD)

        !CALL RecWrite_basis(BasisnD)

        write(out_unitp,*) 'nb_act1_RPH',nb_act1_RPH
        write(out_unitp,*) 'nb_inact21_RPH',nb_inact21_RPH

        CALL Write_RPHTransfo(mole%RPHTransfo)

      END IF

      DO it=mole%nb_Qtransfo-1,mole%itRPH+1,-1
        tab_skip_transfo(it) = mole%tab_Qtransfo(it)%skip_transfo
        mole%tab_Qtransfo(it)%skip_transfo = .TRUE.
      END DO

      ! for RPHpara_AT_Qref
      IF (.NOT. associated(mole%RPHTransfo%RPHpara_AT_Qref)) THEN
        IF (debug) write(out_unitp,*) ' RPHpara_AT_Qref'
        CALL alloc_array(mole%RPHTransfo%RPHpara_AT_Qref,(/ 1 /),       &
                        'mole%RPHTransfo%RPHpara_AT_Qref',name_sub)

        CALL get_Qact(Qact,mole%ActiveTransfo) ! rigid, flexible coordinates

        CALL Set_RPHpara_AT_Qact1(mole%RPHTransfo%RPHpara_AT_Qref(1),   &
                                  Qact,para_Tnum,mole,mole%RPHTransfo)
        mole%RPHTransfo%init_Qref = .TRUE.

        write(out_unitp,*) ' Frequencies, normal modes at the reference geometry'

        write(out_unitp,11) Qact(1:nb_act1_RPH), &
               mole%RPHTransfo%RPHpara_AT_Qref(1)%dnEHess%d0(:)*get_Conv_au_TO_unit('E','cm-1')
 11     format(' frequencies : ',30f10.4)
        write(out_unitp,*) 'dnQopt'
        CALL Write_dnVec(mole%RPHTransfo%RPHpara_AT_Qref(1)%dnQopt)
        write(out_unitp,*) 'dnC_inv'
        CALL Write_dnMat(mole%RPHTransfo%RPHpara_AT_Qref(1)%dnC_inv)
        CALL flush_perso(out_unitp)

      END IF

      ! Check if the nb_act1_RPH coordinates belong to one basis set (primitive ?)
      ! 1) RPHTransfo MUST be the 2d transformation after the active one.
      !write(out_unitp,*) 'asso RPH, itRPH,nb_Qtransfo',associated(mole%RPHTransfo),mole%itRPH,mole%nb_Qtransfo
      RPHCoord_IN_OneBasis = associated(mole%RPHTransfo) .AND. (mole%itRPH == mole%nb_Qtransfo-1)

      RPHCoord_IN_OneBasis = RPHCoord_IN_OneBasis .AND.                 &
        (count(mole%RPHTransfo%list_act_OF_Qdyn(1:nb_act1_RPH) == 1) == nb_act1_RPH)
      !write(out_unitp,*) 'list_act_OF_Qdyn',mole%RPHTransfo%list_act_OF_Qdyn


      ! 2) basis functions of BasisnD are defined as a product (BasisnD%nb_basis > 0)
      !  => if true, RPHCoord_IN_OneBasis CAN be true
      RPHCoord_IN_OneBasis = RPHCoord_IN_OneBasis .AND. (BasisnD%nb_basis > 0)

      ! 3) Check nb_act1_RPH coordinates belong to one primitive basis set
      IF (RPHCoord_IN_OneBasis) THEN
        DO ib=1,BasisnD%nb_basis
          !write(out_unitp,*) 'ib,iQdyn',ib,':',BasisnD%tab_Pbasis(ib)%Pbasis%iQdyn(:)
          IF (BasisnD%tab_Pbasis(ib)%Pbasis%ndim == nb_act1_RPH) THEN
            IF (all(BasisnD%tab_Pbasis(ib)%Pbasis%iQdyn ==              &
                  mole%RPHTransfo%list_QactTOQdyn(1:nb_act1_RPH)) ) EXIT
          END IF
        END DO
        RPHCoord_IN_OneBasis = RPHCoord_IN_OneBasis .AND. (ib <= BasisnD%nb_basis)
      END IF

      IF (RPHCoord_IN_OneBasis) THEN
        CALL Set_paraPRH_ONEBasis(mole,para_Tnum,BasisnD,ib)
      ELSE
        CALL Set_paraPRH_gene(mole,para_Tnum,BasisnD)
      END IF

      CALL flush_perso(out_unitp)
      mole%RPHTransfo%init = .TRUE.

      DO it=mole%nb_Qtransfo-1,mole%itRPH+1,-1
        mole%tab_Qtransfo(it)%skip_transfo = tab_skip_transfo(it)
      END DO


!     -------------------------------------------------------
      IF (debug) THEN
        CALL Write_RPHTransfo(mole%RPHTransfo)
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -------------------------------------------------------
      END SUBROUTINE Set_paraPRH
      SUBROUTINE Set_paraPRH_OneBasis(mole,para_Tnum,BasisnD,ib)
      USE mod_system
      USE mod_nDindex
      USE mod_dnSVM
      USE mod_Constant
      USE mod_PrimOp
      USE mod_basis
      IMPLICIT NONE

!
!=====================================================================
!
!     variables
!
!=====================================================================

!----- variables for the construction of H ---------------------------
      TYPE (basis)   :: BasisnD
      TYPE (CoordType) :: mole
      TYPE (Tnum)    :: para_Tnum
      integer, intent(in) :: ib ! index of the basis in tab_Pbasis(:) or tab_basisPrimSG(:,:)


!------ working variables ---------------------------------
      integer :: L,iq,iq_list,nb_act1_RPH,nb_inact21_RPH,it,comp,iq_list_small,nq

      real (kind=Rkind), allocatable :: Qact1_fromBasisnD(:)

      real (kind=Rkind) :: Qdyn(mole%nb_var)
      real (kind=Rkind) :: Qact(mole%nb_var)
      real (kind=Rkind), allocatable :: List_Qact1(:,:),List_tmp_Qact1(:,:)

      logical :: Find_in_List,iqLoop_end


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Set_paraPRH_OneBasis'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

      nb_act1_RPH    = mole%RPHTransfo%nb_act1
      nb_inact21_RPH = mole%RPHTransfo%nb_inact21

      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'ib',ib

        !CALL RecWrite_basis(BasisnD)

        write(out_unitp,*) 'nb_act1_RPH',nb_act1_RPH
        write(out_unitp,*) 'nb_inact21_RPH',nb_inact21_RPH

        CALL Write_RPHTransfo(mole%RPHTransfo)

      END IF

      ! for tab_RPHpara_AT_Qact1
      !----------------------------------------------------------------
      !--- First the number of grid points ----------------------------
      CALL time_perso('Grid RPH')
      write(out_unitp,*) 'Grid RPH'

      CALL alloc_NParray(Qact1_fromBasisnD,(/ nb_act1_RPH /),'Qact1_fromBasisnD',name_sub)

      iqLoop_end = .FALSE.
      iq = 1
      L  = 0
      SELECT CASE (BasisnD%SparseGrid_type)
      CASE (0) ! normal direct product basis/grid
        nq = get_nq_FROM_basis(BasisnD%tab_Pbasis(ib)%Pbasis) ! it's used one for SparseGrid_type=0
      CASE (1) ! First Smolyak
           STOP ' SG1 not yet'
      CASE (2,4) ! First Smolyak
        nq = get_nq_FROM_basis(BasisnD%tab_basisPrimSG(L,ib))
      CASE Default
         STOP ' no default'
      END SELECT
      write(out_unitp,*) 'L,ib,nq',L,ib,nq ; flush(out_unitp)



      DO

        SELECT CASE (BasisnD%SparseGrid_type)
        CASE (0) ! normal direct product basis/grid
          Qact1_fromBasisnD = BasisnD%tab_Pbasis(ib)%Pbasis%x(:,iq)
          iq = iq + 1
          iqLoop_end = (iq > nq)
        CASE (1) ! First Smolyak
           STOP ' SG1 not yet'
        CASE (2,4) ! First Smolyak

          Qact1_fromBasisnD = BasisnD%tab_basisPrimSG(L,ib)%x(:,iq)
          iq = iq + 1

          IF (iq > nq) THEN
            iq = 1
            L  = L + 1
            IF (L <= BasisnD%L_SparseGrid)                              &
                   nq = get_nq_FROM_basis(BasisnD%tab_basisPrimSG(L,ib))
          END IF

          iqLoop_end = (L > BasisnD%L_SparseGrid)

        CASE Default
           STOP ' no default'
        END SELECT

        !write(out_unitp,*) 'L,iq,Qact1_fromBasisnD',L,iq-1,':',Qact1_fromBasisnD


        IF (.NOT. allocated(List_Qact1)) THEN ! first point
          CALL alloc_NParray(List_Qact1,(/ nb_act1_RPH, 1 /),'List_tmp_Qact1',name_sub)
          List_Qact1(:,1) = Qact1_fromBasisnD(:)
        ELSE

          Find_in_List = .FALSE.
          iq_list_small = 0
          DO iq_list=1,size(List_Qact1,dim=2)
            comp = compar_GridPoint(List_Qact1(:,iq_list),Qact1_fromBasisnD,nb_act1_RPH)
            IF (comp == -1) iq_list_small = iq_list
            Find_in_List = (comp == 0)

            IF (Find_in_List) EXIT
          END DO

          IF (.NOT. Find_in_List) THEN ! add the new point in the list

            CALL alloc_NParray(List_tmp_Qact1,(/ nb_act1_RPH, size(List_Qact1,dim=2)+1 /),'List_tmp_Qact1',name_sub)

            ! find the position to add the point
            IF (iq_list_small == 0) THEN ! add in the first point
              List_tmp_Qact1(:,1) = Qact1_fromBasisnD(:)
              List_tmp_Qact1(:,2:size(List_Qact1,dim=2)+1) = List_Qact1(:,:)
            ELSE IF (iq_list_small == size(List_Qact1,dim=2)) THEN ! add in the last point
              List_tmp_Qact1(:,1:iq_list_small)   = List_Qact1(:,1:iq_list_small)
              List_tmp_Qact1(:,  iq_list_small+1) = Qact1_fromBasisnD(:)
            ELSE
              List_tmp_Qact1(:,1:iq_list_small)   = List_Qact1(:,1:iq_list_small)
              List_tmp_Qact1(:,  iq_list_small+1) = Qact1_fromBasisnD(:)
              List_tmp_Qact1(:,  iq_list_small+2:size(List_Qact1,dim=2)+1) = List_Qact1(:,iq_list_small+1:size(List_Qact1,dim=2))
            END IF
            CALL dealloc_NParray(List_Qact1,'List_Qact1',name_sub)

            CALL alloc_NParray(List_Qact1,shape(List_tmp_Qact1),'List_Qact1',name_sub)
            List_Qact1(:,:) = List_tmp_Qact1
            CALL dealloc_NParray(List_tmp_Qact1,'List_tmp_Qact1',name_sub)

          END IF


        END IF

        IF (iqLoop_end) EXIT


      END DO
      CALL time_perso('Grid RPH')
      write(out_unitp,*) 'nb of Qact1 grid points',size(List_Qact1,dim=2)
      CALL flush_perso(out_unitp)

      !----------------------------------------------------------------
      !----------------------------------------------------------------

      !----------------------------------------------------------------
      !---- allocation of tab_RPHpara_AT_Qact1 ------------------------
      mole%RPHTransfo%nb_Qa = size(List_Qact1,dim=2)
      CALL alloc_array(mole%RPHTransfo%tab_RPHpara_AT_Qact1,(/ mole%RPHTransfo%nb_Qa /),&
                      'mole%RPHTransfo%tab_RPHpara_AT_Qact1',name_sub)
      !----------------------------------------------------------------


      !----------------------------------------------------------------
      !---- Multidimensional loop -------------------------------------
        DO iq_list=1,size(List_Qact1,dim=2)

          CALL get_Qact(Qact,mole%ActiveTransfo) ! rigid, flexible???? coordinates
          Qact(1:nb_act1_RPH) = List_Qact1(:,iq_list)

          write(out_unitp,*) 'new RPH point',iq_list
          CALL flush_perso(out_unitp)

          CALL Set_RPHpara_AT_Qact1(                                    &
                          mole%RPHTransfo%tab_RPHpara_AT_Qact1(iq_list),&
                           Qact,para_Tnum,mole,mole%RPHTransfo)

        END DO

      write(out_unitp,*) 'END Grid RPH'
      CALL flush_perso(out_unitp)

      CALL dealloc_NParray(List_Qact1,'List_Qact1',name_sub)
      CALL dealloc_NParray(Qact1_fromBasisnD,'Qact1_fromBasisnD',name_sub)

!     -------------------------------------------------------
      IF (debug) THEN
        CALL Write_RPHTransfo(mole%RPHTransfo)
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -------------------------------------------------------
      END SUBROUTINE Set_paraPRH_OneBasis
      SUBROUTINE Set_paraPRH_gene(mole,para_Tnum,BasisnD)
      USE mod_system
      USE mod_nDindex
      USE mod_dnSVM
      USE mod_Constant
      USE mod_PrimOp
      USE mod_basis
      IMPLICIT NONE

!
!=====================================================================
!
!     variables
!
!=====================================================================

!----- variables for the construction of H ---------------------------
      TYPE (basis)   :: BasisnD
      TYPE (CoordType) :: mole
      TYPE (Tnum)    :: para_Tnum


!------ working variables ---------------------------------
      integer :: nq_part,iq,iq_list,nb_act1_RPH,nb_inact21_RPH,it,comp,iq_list_small

      real (kind=Rkind), allocatable :: Qact1_fromBasisnD(:)

      real (kind=Rkind) :: Qdyn(mole%nb_var)
      real (kind=Rkind) :: Qact(mole%nb_var),auTOcm_inv
      real (kind=Rkind), allocatable :: List_Qact1(:,:),List_tmp_Qact1(:,:)

      logical :: Find_in_List

      TYPE (OldParam) :: OldPara ! to iq, ...



!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Set_paraPRH_gene'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------

      nb_act1_RPH    = mole%RPHTransfo%nb_act1
      nb_inact21_RPH = mole%RPHTransfo%nb_inact21

      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub

        write(out_unitp,*) 'nb_qa',get_nq_FROM_basis(BasisnD)

        !CALL RecWrite_basis(BasisnD)

        write(out_unitp,*) 'nb_act1_RPH',nb_act1_RPH
        write(out_unitp,*) 'nb_inact21_RPH',nb_inact21_RPH

        CALL Write_RPHTransfo(mole%RPHTransfo)

      END IF

      ! for tab_RPHpara_AT_Qact1
      !----------------------------------------------------------------
      !--- First the number of grid points ----------------------------
      CALL time_perso('Grid RPH')
      write(out_unitp,*) 'Grid RPH'

      CALL alloc_NParray(Qact1_fromBasisnD,(/ nb_act1_RPH /),'Qact1_fromBasisnD',name_sub)

      nq_part = get_nq_FROM_basis(BasisnD)/100
      DO iq=1,get_nq_FROM_basis(BasisnD)

        IF (debug) THEN
        IF (mod(iq,nq_part)==0) THEN
          write(out_unitp,*) 'iq,nq',iq,get_nq_FROM_basis(BasisnD)
          CALL flush_perso(out_unitp)
        END IF
        END IF

        Qact(:) = ZERO
        CALL Rec_Qact(Qact,BasisnD,iq,mole,OldPara)
        !write(out_unitp,*) 'iq,size(List_Qact1,dim=2),Qact',iq,size(List_Qact1,dim=2),Qact
        !flush(out_unitp)
        CALL Qact_TO_Qdyn_FROM_ActiveTransfo(Qact,Qdyn,mole%ActiveTransfo)
        !write(out_unitp,*) 'Qdyn',Qdyn
        !flush(out_unitp)

        Qact1_fromBasisnD(:) = Qdyn(mole%RPHTransfo%list_QactTOQdyn(1:nb_act1_RPH))
        !write(out_unitp,*) 'Qact1_fromBasisnD',Qact1_fromBasisnD


        IF (.NOT. allocated(List_Qact1)) THEN ! first point
          CALL alloc_NParray(List_Qact1,(/ nb_act1_RPH, 1 /),'List_tmp_Qact1',name_sub)
          List_Qact1(:,1) = Qact1_fromBasisnD(:)
        ELSE

          Find_in_List = .FALSE.
          iq_list_small = 0
          DO iq_list=1,size(List_Qact1,dim=2)
            comp = compar_GridPoint(List_Qact1(:,iq_list),Qact1_fromBasisnD,nb_act1_RPH)
            IF (comp == -1) iq_list_small = iq_list
            Find_in_List = (comp == 0)
            IF (Find_in_List) EXIT
          END DO

          IF (.NOT. Find_in_List) THEN ! add the new point in the list

            CALL alloc_NParray(List_tmp_Qact1,(/ nb_act1_RPH, size(List_Qact1,dim=2)+1 /),'List_tmp_Qact1',name_sub)

            ! find the position to add the point
            IF (iq_list_small == 0) THEN ! add in the first point
              List_tmp_Qact1(:,1) = Qact1_fromBasisnD(:)
              List_tmp_Qact1(:,2:size(List_Qact1,dim=2)+1) = List_Qact1(:,:)
            ELSE IF (iq_list_small == size(List_Qact1,dim=2)) THEN ! add in the last point
              List_tmp_Qact1(:,1:iq_list_small) = List_Qact1(:,1:iq_list_small)
              List_tmp_Qact1(:,  iq_list_small+1) = Qact1_fromBasisnD(:)
            ELSE
              List_tmp_Qact1(:,1:iq_list_small) = List_Qact1(:,1:iq_list_small)
              List_tmp_Qact1(:,  iq_list_small+1) = Qact1_fromBasisnD(:)
              List_tmp_Qact1(:,  iq_list_small+2:size(List_Qact1,dim=2)+1) = List_Qact1(:,iq_list_small+1:size(List_Qact1,dim=2))
            END IF
            CALL dealloc_NParray(List_Qact1,'List_Qact1',name_sub)

            CALL alloc_NParray(List_Qact1,shape(List_tmp_Qact1),'List_Qact1',name_sub)
            List_Qact1 = List_tmp_Qact1
            CALL dealloc_NParray(List_tmp_Qact1,'List_tmp_Qact1',name_sub)

          END IF


        END IF

      END DO
      CALL time_perso('Grid RPH')
      write(out_unitp,*) 'nb of Qact1 grid points',size(List_Qact1,dim=2)
      CALL flush_perso(out_unitp)
      !----------------------------------------------------------------
      !----------------------------------------------------------------

      !----------------------------------------------------------------
      !---- allocation of tab_RPHpara_AT_Qact1 ------------------------
      mole%RPHTransfo%nb_Qa = size(List_Qact1,dim=2)
      CALL alloc_array(mole%RPHTransfo%tab_RPHpara_AT_Qact1,(/ mole%RPHTransfo%nb_Qa /),&
                      'mole%RPHTransfo%tab_RPHpara_AT_Qact1',name_sub)
      !----------------------------------------------------------------


      !----------------------------------------------------------------
      !---- Multidimensional loop -------------------------------------
        DO iq_list=1,size(List_Qact1,dim=2)

          CALL get_Qact(Qact,mole%ActiveTransfo) ! rigid, flexible???? coordinates
          Qact(1:nb_act1_RPH) = List_Qact1(:,iq_list)

          write(out_unitp,*) 'new RPH point',iq_list
          CALL flush_perso(out_unitp)

          CALL Set_RPHpara_AT_Qact1(                                    &
                          mole%RPHTransfo%tab_RPHpara_AT_Qact1(iq_list),&
                          Qact,para_Tnum,mole,mole%RPHTransfo)

        END DO

      write(out_unitp,*) 'END Grid RPH'
      CALL flush_perso(out_unitp)

      CALL dealloc_NParray(List_Qact1,'List_Qact1',name_sub)
      CALL dealloc_NParray(Qact1_fromBasisnD,'Qact1_fromBasisnD',name_sub)

!     -------------------------------------------------------
      IF (debug) THEN
        CALL Write_RPHTransfo(mole%RPHTransfo)
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -------------------------------------------------------
      END SUBROUTINE Set_paraPRH_gene

      FUNCTION compar_GridPoint(Q1,Q2,n)
      USE mod_system
      IMPLICIT NONE
        integer :: compar_GridPoint
        integer, intent(in) :: n
        real(kind=Rkind), intent(in) :: Q1(n),Q2(n)
        integer :: i

        DO i=1,n
          IF (abs(Q1(i)-Q2(i)) < ONETENTH**10) THEN
            compar_GridPoint = 0
          ELSE
            IF (Q1(i) < Q2(i)) THEN
              compar_GridPoint = -1
            ELSE
              compar_GridPoint =  1
            END IF
            EXIT
          END IF
        END DO
      END FUNCTION compar_GridPoint

END MODULE mod_Set_paraRPH
