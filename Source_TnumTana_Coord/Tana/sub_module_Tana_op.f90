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
   MODULE mod_Tana_op
   !Description:
   USE mod_system
   USE mod_Tana_OpEl
   USE mod_Tana_Op1D
   USE mod_Tana_OpnD
   USE mod_Tana_sum_opnd
   USE mod_Tana_VecSumOpnD
   USE mod_Tana_PiEulerRot
   USE mod_Tana_vec_operations
   IMPLICIT NONE
   PRIVATE
   PUBLIC ::  get_opLi, get_opL1, get_opL2, get_keo_for_Qactiv,         &
              get_opKEO, get_opKEO_subsyst, add_Vextr, add_Vextr_new,   &
              Get_F2_F1_FROM_TWOxKEO

   CONTAINS 

   !> @description: Defines the total angular momentum in terms of
   !!               derivative operator of the Euler's angles in the SF frame
   !!               for a subsystem
   !!               which is defined with only two Euler's angles (beta, gamma).
   !!               Eq A10 of Ndong et al J. Chem. Phys. 136,034107 (2012).
   !> @param:       opJ    The  vector (type: vec_sum_opnd)
   !> @param:       F1_sum Sum of elementary op 
   !> @param:       fbeta   an elementary op  which contains the needed
   !!                      information on \beta or ub coordinate
   !> @param:       fgamma   an elementary op  which contains the needed
   !!                      information on \gamma 
   !> @param:       dag   Logical, if present and = true, the adjoint of Li will be obtained
   !!                        by  vector transpose time a matrix
   SUBROUTINE get_opJ_projected_into_ref_frameEq170(opJ, F1_sum, fbeta, fgamma, dag)
     type(vec_sum_opnd),      intent(inout)      :: opJ
     type(sum_opnd),          intent(in)         :: F1_sum
     type(opel),              intent(in)         :: fbeta
     type(opel),              intent(in)         :: fgamma
     logical, optional,       intent(in)         :: dag

     type(sum_opnd), allocatable         :: M_opnd(:,:)
     type(sum_opnd), allocatable         :: M_opnd_tr(:,:)
     type(vec_sum_opnd)                  :: V

     logical                         :: dag_loc
     integer                         :: error
     character (len = *), parameter  :: routine_name='get_opJ_projected_into_ref_frameEq170'

     IF (present(dag)) THEN
       dag_loc = dag
     ELSE
       dag_loc = .FALSE.
     END IF

     if(abs(fbeta%idq) /= 7 .or. fgamma%idq /= 8) then
       write(out_unitp,*) ' ERROR in',routine_name
       write(out_unitp,*) "Data structure of idq of fbeta or fgamma are not correct"
       STOP
     end if
     if(fgamma%idf /= 1 .or. fbeta%idf /= 1) then
       write(out_unitp,*) ' ERROR in',routine_name
       write(out_unitp,*) "The elementary operators should be the Id"
       STOP
     end if 

     call allocate_op(opJ, 3)
     call allocate_op(V, 3)

     V%vec_sum(1) = F1_sum
     IF (dag_loc) THEN
       V%vec_sum(2) = get_Pq_dag(fbeta)
       V%vec_sum(3) = get_Pq_dag(fgamma)
     ELSE
       V%vec_sum(2) = get_Pq(fbeta)
       V%vec_sum(3) = get_Pq(fgamma)
     END IF


     CALL alloc_NParray(M_opnd,(/3,3/),'M_opnd',routine_name)

     M_opnd(1,1) = get_cot(fbeta)
     M_opnd(1,1)%Cn(1) = -CONE
     M_opnd(1,2) = CZERO
     M_opnd(1,3) = get_sin(fbeta,-1)

     M_opnd(2,1) = CZERO
     if(fbeta%idq == 7) then
       M_opnd(2,2) = CONE
     else
       M_opnd(2,2) = get_sin(fbeta)
       M_opnd(2,2)%Cn(1) = -CONE
     end if
     M_opnd(2,3) = CZERO

     M_opnd(3,1) = CONE
     M_opnd(3,2) = CZERO
     M_opnd(3,3) = CZERO


     if(dag_loc) then
       M_opnd_tr = Transpose_Mat_OF_sum_opnd(M_opnd)

       call V_times_M_opnd_in_Vres(V, M_opnd_tr, opJ)

     else
       call M_opnd_times_V_in_Vres(M_opnd, V, opJ)
     end if

     call delete_op(V)
     CALL dealloc_NParray(M_opnd,'M_opnd',routine_name)
     CALL dealloc_NParray(M_opnd_tr,'M_opnd_tr',routine_name)

   END SUBROUTINE get_opJ_projected_into_ref_frameEq170



   !! @description: Defines the total angular momentum in terms of
   !!               derivative operator of the Euler's angles.
   !!               Eq A7 of Ndong et al. J. Chem. Phys. 136,034107 (2012).
   !! @param:       OpJ      The  vector (type: vec_sum_opnd)
   !! @param:       falpha   an elementary op  which contains the needed
   !!                         information on \alpha or u coordinate
   !! @param:       fbeta    an elementary op  which contains the needed
   !!                         information on \beta or ub coordinate
   !! @param:       fgamma   an elementary op  which contains the needed
   !!                         information on \gamma
   !! @param:       dag   Logical, if present and = true, the adjoint of Li will be obtained
   !!                        by  vector transpose time a matrix
   SUBROUTINE get_opJ_projected_into_ref_frame(opJ, falpha, fbeta, fgamma, dag)
     type(vec_sum_opnd),      intent(inout)      :: opJ
     type(opel),              intent(in)         :: falpha
     type(opel),              intent(in)         :: fbeta
     type(opel),              intent(in)         :: fgamma
     logical, optional,       intent(in)         :: dag

     type(sum_opnd), allocatable         :: M_opnd(:,:)
     type(sum_opnd), allocatable         :: M_opnd_tr(:,:)
     type(vec_sum_opnd)                  :: V

     logical                         :: dag_loc
     integer                         :: error
     character (len = Name_longlen)  :: routine_name='get_opJ_projected_into_ref_frame'

     IF (present(dag)) THEN
       dag_loc = dag
     ELSE
       dag_loc = .FALSE.
     END IF

     if(falpha%idq /= 6 .or. abs(fbeta%idq) /= 7 .or. fgamma%idq /= 8 ) then
       write(out_unitp,*) ' ERROR in',routine_name
       write(out_unitp,*) "Data structure of idq in falpha or fbeta or fgamma are not correct"
       STOP
     end if
     if(falpha%idf /= 1 .or. fbeta%idf /= 1 .OR. fgamma%idf /= 1 ) then
       write(out_unitp,*) ' ERROR in',routine_name
       write(out_unitp,*) "The elementary operators should be the Id"
       STOP
     end if


     call allocate_op(opJ, 3)
     call allocate_op(V,   3)

     IF (dag_loc) THEN
       V%vec_sum(1) = get_Pq_dag(falpha)
       V%vec_sum(2) = get_Pq_dag(fbeta)
       V%vec_sum(3) = get_Pq_dag(fgamma)
     ELSE
       V%vec_sum(1) = get_Pq(falpha)
       V%vec_sum(2) = get_Pq(fbeta)
       V%vec_sum(3) = get_Pq(fgamma)
     END IF


     CALL alloc_NParray(M_opnd,(/3,3/),'M_opnd',routine_name)

     if(fbeta%idq == 7) then
       M_opnd(1,2) = get_sin(falpha)
       M_opnd(1,2)%Cn(1) = -CONE
       M_opnd(2,2) = get_cos(falpha)
     else
       M_opnd(1,2) = get_sin(falpha) * get_sin(fbeta)
       M_opnd(2,2) = get_cos(falpha) * get_sin(fbeta)
       M_opnd(2,2)%Cn(1) = -CONE
     end if

     M_opnd(1,1) = get_cos(falpha) * get_cot(fbeta)
     M_opnd(1,1)%Cn(1) = -CONE
     M_opnd(1,3) = get_cos(falpha) * get_sin(fbeta,-1)


     M_opnd(2,1) = get_sin(falpha) * get_cot(fbeta)
     M_opnd(2,1)%Cn(1) = -CONE
     M_opnd(2,3) = get_sin(falpha) * get_sin(fbeta,-1)

     M_opnd(3,1) = CONE
     M_opnd(3,2) = CZERO
     M_opnd(3,3) = CZERO

     if(dag_loc) then
       M_opnd_tr = Transpose_Mat_OF_sum_opnd(M_opnd)

       call V_times_M_opnd_in_Vres(V, M_opnd_tr, opJ)

     else
       call M_opnd_times_V_in_Vres(M_opnd, V, opJ)
     end if


     call delete_op(V)
     CALL dealloc_NParray(M_opnd,'M_opnd',routine_name)
     CALL dealloc_NParray(M_opnd_tr,'M_opnd_tr',routine_name)

   END SUBROUTINE get_opJ_projected_into_ref_frame


   !! @description: Defines the angular momentum Li in the reference frame (not used)
   !!               Initially, they are defined in the local BF frame
   !! @param:       L      The  vector (type: vec_sum_opnd)
   !! @param:       falpha   an elementary op  which contains the needed
   !!                         information on \alpha or u coordinate
   !! @param:       fbeta    an elementary op  which contains the needed
   !!                         information on \beta or ub coordinate
   !! @param:       fgamma   an elementary op  which contains the needed
   !!                         information on \gamma
   !! @param:       dag   Logical, if present and = true, the adjoint of Li will be obtained
   !!                        by  vector transpose time a matrix
   SUBROUTINE get_opLi_projected_into_ref_frame(opLi, falpha, fbeta, fgamma, dag)
     type(vec_sum_opnd),      intent(inout)      :: opLi
     type(opel),              intent(in)         :: falpha
     type(opel),              intent(in)         :: fbeta
     type(opel),              intent(in)         :: fgamma
     logical, optional,       intent(in)         :: dag

     type(sum_opnd), allocatable         :: M_opnd(:,:)
     type(sum_opnd), allocatable         :: M_opnd_tr(:,:)
     type(vec_sum_opnd)                  :: V

     logical                         :: dag_loc
     integer                         :: error
     character (len=*), parameter    :: routine_name='get_opLi_projected_into_ref_frame'

     IF (present(dag)) THEN
       dag_loc = dag
     ELSE
       dag_loc = .FALSE.
     END IF

     if(falpha%idq /= 6 .or. abs(fbeta%idq) /= 7 .or. fgamma%idq /= 8 ) then
       write(out_unitp,*) ' ERROR in',routine_name
       write(out_unitp,*) "Data structure of idq in falpha or fbeta or fgamma are not correct"
       STOP
     end if
     if(falpha%idf /= 1 .or. fbeta%idf /= 1 .OR. fgamma%idf /= 1 ) then
       write(out_unitp,*) ' ERROR in',routine_name
       write(out_unitp,*) "The elementary operators should be the Id"
       STOP
     end if


     call allocate_op(opLi, 3)
     call allocate_op(V, 3)

     IF (dag_loc) THEN
       V%vec_sum(1) = get_Pq_dag(falpha)
       V%vec_sum(2) = get_Pq_dag(fbeta)
       V%vec_sum(3) = get_Pq_dag(fgamma)
     ELSE
       V%vec_sum(1) = get_Pq(falpha)
       V%vec_sum(2) = get_Pq(fbeta)
       V%vec_sum(3) = get_Pq(fgamma)
     END IF


     CALL alloc_NParray(M_opnd,(/3,3/),'M_opnd',routine_name)

     if(fbeta%idq == 7) then
       M_opnd(1,2) = get_sin(falpha)
       M_opnd(1,2)%Cn(1) = -CONE
       M_opnd(2,2) = get_cos(falpha)
     else
       M_opnd(1,2) = get_sin(falpha) * get_sin(fbeta)
       M_opnd(2,2) = get_cos(falpha) * get_sin(fbeta)
       M_opnd(2,2)%Cn(1) = -CONE
     end if

     M_opnd(1,1) = get_cos(falpha) * get_cot(fbeta)
     M_opnd(1,1)%Cn(1) = -CONE
     M_opnd(1,3) = get_cos(falpha) * get_sin(fbeta,-1)


     M_opnd(2,1) = get_sin(falpha) * get_cot(fbeta)
     M_opnd(2,1)%Cn(1) = -CONE
     M_opnd(2,3) = get_sin(falpha) * get_sin(fbeta,-1)

     M_opnd(3,1) = CONE
     M_opnd(3,2) = CZERO
     M_opnd(3,3) = CZERO

     if(dag_loc) then
       M_opnd_tr = Transpose_Mat_OF_sum_opnd(M_opnd)

       call V_times_M_opnd_in_Vres(V, M_opnd_tr, opLi)

     else
       call M_opnd_times_V_in_Vres(M_opnd, V, opLi)
     end if

     call delete_op(V)
     CALL dealloc_NParray(M_opnd_tr,'M_opnd_tr',routine_name)
     CALL dealloc_NParray(M_opnd,'M_opnd',routine_name)

   END SUBROUTINE get_opLi_projected_into_ref_frame

   !! @description: Defines the total angular momentum in terms of
   !!               derivative operator of the Euler's angles.
   !! @param:       J      The  vector (type: vec_sum_opnd)
   !! @param:       F1el   an elementary op  which contains the needed
   !!                      information on \alpha or u coordinate
   !! @param:       F2el   an elementary op  which contains the needed
   !!                      information on \beta or ub coordinate
   !! @param:       F3el   an elementary op  which contains the needed
   !!                      information on \gamma 
   !! @param:       dag   Logical, if present and = true, the adjoint of Li will be obtained
   !!                        by  vector transpose time a matrix
   SUBROUTINE get_opJ_projected_into_BF(opJ, F1el, F2el, F3el, dag)
     type(vec_sum_opnd),      intent(inout)      :: opJ
     type(opel),              intent(in)         :: F1el
     type(opel),              intent(in)         :: F2el
     type(opel), optional,    intent(in)         :: F3el
     logical, optional,       intent(in)         :: dag

     type(sum_opnd), allocatable         :: M_opnd(:,:)
     type(sum_opnd), allocatable         :: M_opnd_tr(:,:)
     type(vec_sum_opnd)              :: V


     logical                         :: dag_loc
     integer                         :: i, j
     integer                         :: error
     character (len = *), parameter  :: routine_name='get_opJ_projected_into_BF'

     IF (present(dag)) THEN
       dag_loc = dag
     ELSE
       dag_loc = .FALSE.
     END IF

     if(F1el%idq /=6 .or. (F2el%idq /= 7 .and. F2el%idq /=-7)) then
       write(out_unitp,*) ' ERROR in',routine_name
       write(out_unitp,*) "Data structure of idq in F1el or in F2el is not correct"
       STOP
     end if        
     if(present(F3el)) then
       if(F3el%idq /= 8 ) then
         write(out_unitp,*) ' ERROR in',routine_name
         write(out_unitp,*) "Data structure of idq in F3el is not correct"
         STOP
       end if
     end if        
     if(F1el%idf /= 1 .or. F2el%idf /= 1) then
       write(out_unitp,*) ' ERROR in',routine_name
       write(out_unitp,*) "The elementary operators should be the Id"
       STOP
     end if 
     if(present(F3el)) then
       if(F3el%idf /= 1 ) then
         write(out_unitp,*) ' ERROR in',routine_name
         write(out_unitp,*) "The elementary operators should be the Id"
         STOP
       end if
     end if 

     if(present(F3el)) then
       call allocate_op(opJ, 3)
       call allocate_op(V, 3)

       IF (dag_loc) THEN
         V%vec_sum(1) = get_Pq_dag(F1el) ! Pqb
         V%vec_sum(2) = get_Pq_dag(F2el) ! Pqb
         V%vec_sum(3) = get_Pq_dag(F3el) ! Pqg
       ELSE
         V%vec_sum(1) = get_Pq(F1el) ! Pqb
         V%vec_sum(2) = get_Pq(F2el) ! Pqb
         V%vec_sum(3) = get_Pq(F3el) ! Pqg
       END IF


       CALL alloc_NParray(M_opnd,(/3,3/),'M_opnd',routine_name)

       if(F2el%idq == 7) then

         M_opnd(1,2) = get_sin(F3el) ! sin(g)

         M_opnd(2,2) = get_cos(F3el) ! cos(g)

       else

         M_opnd(1,2) = get_sin(F2el) * get_sin(F3el) ! sin(b) * sin(g)
         M_opnd(1,2)%Cn(1) = -CONE

         M_opnd(2,2) = get_sin(F2el) * get_cos(F3el) ! sin(b) * cos(g)
         M_opnd(2,2)%Cn(1) = -CONE

       end if

       M_opnd(3,1) = czero
       M_opnd(3,2) = czero
       M_opnd(3,3) = cone

       M_opnd(1,1) = get_sin(F2el,-1) * get_cos(F3el) ! sin(b)^-1 * cos(g)
       M_opnd(1,1)%Cn(1) = -CONE

       !CALL get_F1_1d_times_F2_nd_to_Fres_nd(get_cos(F3el), &
       !          (get_cos(F2el) * get_sin(F2el,-1) ), M_opnd(1,3))
       M_opnd(1,3) = (get_cos(F2el) * get_sin(F2el,-1) ) * get_cos(F3el) ! cot(b) * cos(g)

       M_opnd(2,1) = get_sin(F2el,-1) * get_sin(F3el) ! sin(b)^-1 * sin(g)

       M_opnd(2,3) = (get_cos(F2el) * get_sin(F2el,-1) ) * get_sin(F3el) ! cot(b) * sin(g)
       M_opnd(2,3)%Cn(1) = -CONE

       if(dag_loc) then
         M_opnd_tr = Transpose_Mat_OF_sum_opnd(M_opnd)

         call V_times_M_opnd_in_Vres(V, M_opnd_tr, opJ)

       else
         call M_opnd_times_V_in_Vres(M_opnd, V, opJ)
       end if
     else
       write(out_unitp,*) ' ERROR in ',routine_name
       write(out_unitp,*) "Less than three Euler's angles is not yet taken into account"
       STOP
     end if

     CALL dealloc_NParray(M_opnd,'M_opnd',routine_name)
     CALL dealloc_NParray(M_opnd_tr,'M_opnd_tr',routine_name)

   END SUBROUTINE get_opJ_projected_into_BF




   !! @description: Defines the total angular momentum in terms of
   !!               derivative operator of the Euler's angles for a subsystem
   !!               which is defined with only two Euler's angles.
   !!               Eq A9 (A6) from Ndong et al. J. Chem. Phys. 136,034107 (2012)
   !!               Remark: wrong minus sign if thev(1,3) matrix element
   !! @param:       J      The  vector (type: vec_sum_opnd)
   !! @param:       F1_sum Sum of elementary op 
   !! @param:       fbeta   an elementary op  which contains the needed
   !!                      information on \beta or ub coordinate
   !! @param:       fgamma   an elementary op  which contains the needed
   !!                      information on \gamma 
   !! @param:       dag   Logical, if present and = true, the adjoint of Li will be obtained
   !!                        by  vector transpose time a matrix
   SUBROUTINE get_opJ_projected_into_BFEq171(opJ, F1_sum, fbeta, fgamma, dag)
     type(vec_sum_opnd),      intent(inout)      :: opJ
     type(sum_opnd),          intent(in)         :: F1_sum
     type(opel),              intent(in)         :: fbeta
     type(opel),              intent(in)         :: fgamma
     logical, optional,       intent(in)         :: dag

     type(sum_opnd), allocatable         :: M_opnd(:,:)
     type(sum_opnd), allocatable         :: M_opnd_tr(:,:)
     type(vec_sum_opnd)                  :: V

     logical                         :: dag_loc
     integer                         :: error
     character (len=*), parameter :: routine_name='get_opJ_projected_into_BFEq171'

     IF (present(dag)) THEN
       dag_loc = dag
     ELSE
       dag_loc = .FALSE.
     END IF

     if((fbeta%idq /= 7 .and. fbeta%idq /=-7).or. fgamma%idq /= 8) then
       write(out_unitp,*) ' ERROR in',routine_name
       write(out_unitp,*) "TData structure of idq of fbeta or fgamma is not correct"
       STOP
     end if        
     if(fgamma%idf /= 1 .or. fbeta%idf /= 1) then
       write(out_unitp,*) ' ERROR in',routine_name
       write(out_unitp,*) "The elementary operators should be the Id"
       STOP
     end if 

     call allocate_op(opJ, 3)
     call allocate_op(V, 3)


     V%vec_sum(1) = F1_sum
     IF (dag_loc .AND. fbeta%idq == 7) THEN
       V%vec_sum(2) = get_Pq_dag(fbeta)
       V%vec_sum(3) = get_Pq_dag(fgamma)
     ELSE
       V%vec_sum(2) = get_Pq(fbeta)
       V%vec_sum(3) = get_Pq(fgamma)
     END IF

     CALL alloc_NParray(M_opnd,(/3,3/),'M_opnd',routine_name)

     M_opnd(1,1) = get_sin(fbeta,-1) * get_cos(fgamma)
     M_opnd(1,1)%Cn(1) = -CONE
     M_opnd(1,3) = get_cot(fbeta) * get_cos(fgamma)

     M_opnd(2,1) = get_sin(fbeta,-1) * get_sin(fgamma)
     M_opnd(2,3) = get_cot(fbeta) * get_sin(fgamma)
     M_opnd(2,3)%Cn(1) = -CONE

     if(fbeta%idq == 7) then
       M_opnd(1,2) = get_sin(fgamma)
       M_opnd(2,2) = get_cos(fgamma)
     else
       M_opnd(1,2) = get_sin(fbeta) * get_sin(fgamma)
       M_opnd(1,2)%Cn(1) = -CONE
       M_opnd(2,2) = get_sin(fbeta) * get_cos(fgamma)
       M_opnd(2,2)%Cn(1) = -CONE
     end if

     M_opnd(3,1) = CZERO
     M_opnd(3,2) = CZERO
     M_opnd(3,3) = CONE

     if(dag_loc) then

       M_opnd_tr = Transpose_Mat_OF_sum_opnd(M_opnd)

       call V_times_M_opnd_in_Vres(V, M_opnd_tr, opJ)

     else
       call M_opnd_times_V_in_Vres(M_opnd, V, opJ)
     end if


     call delete_op(V)
     CALL dealloc_NParray(M_opnd,'M_opnd',routine_name)
     CALL dealloc_NParray(M_opnd_tr,'M_opnd_tr',routine_name)

   END SUBROUTINE get_opJ_projected_into_BFEq171


   !! @description: Defines the vector Li, the partial angular momentum
   !!               associated to Ri
   !! @param:       L      The  vector (type: vec_sum_opnd)
   !! @param:       theta   an elementary op  which contains the needed
   !!                      information on \theta or u or \beta coordinate
   !! @param:       phi   an elementary op  which contains the needed
   !!                      information on \phi or \alpha
   !! @param:       index_L  Integer corresponding to the index of L
   !! @param:       dag   Logical, if present and = true, the adjoint of Li will be obtained
   !!                        by  vector transpose time a matrix
   SUBROUTINE get_opLi(L, theta, phi, index_L, dag, Li)
     type(vec_sum_opnd),      intent(inout)      :: L
     type(opel),              intent(in)         :: theta
     type(opel),              intent(in)         :: phi
     integer,                 intent(in)         :: index_L
     logical, optional,       intent(in)         :: dag,Li

     type(sum_opnd), allocatable     :: M_opnd(:,:)
     type(sum_opnd), allocatable     :: M_opnd_tr(:,:)
     type(vec_sum_opnd)              :: V

     logical                         :: dag_loc,Li_loc
     integer                         :: error
     character (len=*), parameter    :: routine_name='get_opLi'

     IF (present(dag)) THEN
       dag_loc = dag
     ELSE
       dag_loc = .FALSE.
     END IF

     IF (present(Li)) THEN
       Li_loc = Li
     ELSE
       Li_loc = .FALSE.
     END IF

     IF (Li_loc) THEN
       call allocate_op(L,3)
       L%vec_sum(1) = get_Lx(theta)
       L%vec_sum(2) = get_Ly(theta)
       L%vec_sum(3) = get_Lz(phi)
       RETURN
     END IF

     if((phi%idq /=4 .and. phi%idq /=6) .or. &
       & (theta%idq /= 3 .and. theta%idq /=-3  .and.&
       &  theta%idq /= 7 .and. theta%idq /=-7)) then
       write(out_unitp,*) ' ERROR in',routine_name
       write(out_unitp,*) "Data structure of idq in theta or in phi is not correct"
       STOP
     end if
     if(theta%idf /= 1 .or. phi%idf /= 1) then
       write(out_unitp,*) ' ERROR in',routine_name
       write(out_unitp,*) "The elementary operators should be the Id"
       STOP
     end if 
     if(index_L < 3) then
       write(out_unitp,*) ' ERROR in',routine_name
       write(out_unitp,*) "The routine evaluates only L_i, i = 3, Ndim_L-1"
       STOP
     end if

     call allocate_op(L, 3)

     call allocate_op(V, 2)
     if(dag_loc) then
       V%vec_sum(1) = get_Pq_dag(theta)
       V%vec_sum(2) = get_Pq_dag(phi)
     else
       V%vec_sum(1) = get_Pq(theta)
       V%vec_sum(2) = get_Pq(phi)
     end if

     CALL alloc_NParray(M_opnd,(/3,2/),'M_opnd',routine_name)

     if(theta%idq == 3 .or. theta%idq == 7) then

       M_opnd(1,1) = get_sin(phi)
       M_opnd(1,1)%Cn(1) = -CONE
       M_opnd(1,2) = get_cos(phi) * get_cot(theta)
       M_opnd(1,2)%Cn(1) = -CONE

       M_opnd(2,1) = get_cos(phi)
       M_opnd(2,2) = get_sin(phi) * get_cot(theta)
       M_opnd(2,2)%Cn(1) = -CONE

       M_opnd(3,1) = czero
       M_opnd(3,2) = cone

     else

       M_opnd(1,1) = get_sin(phi) * get_sin(theta)
       M_opnd(1,1)%Cn(1) =  CONE
       M_opnd(1,2) = get_cos(phi) * get_cot(theta)
       M_opnd(1,2)%Cn(1) = -CONE


       M_opnd(2,1) = get_cos(phi) * get_sin(theta)
       M_opnd(2,1)%Cn(1) = -CONE
       M_opnd(2,2) = get_sin(phi) * get_cot(theta)
       M_opnd(2,2)%Cn(1) = -CONE

       M_opnd(3,1) = czero
       M_opnd(3,2) = cone

     end if

     if(dag_loc) then
       M_opnd_tr = Transpose_Mat_OF_sum_opnd(M_opnd)

       call V_times_M_opnd_in_Vres(V, M_opnd_tr, L)

     else
       call M_opnd_times_V_in_Vres(M_opnd, V, L)
     end if

     call delete_op(V)
     CALL dealloc_NParray(M_opnd_tr,'M_opnd_tr',routine_name)
     CALL dealloc_NParray(M_opnd,'M_opnd',routine_name)

   END SUBROUTINE get_opLi

   !! @description: Defines the vector L1, the partial angular momentum
   !!               associated to R1
   !! @param:       L1     The  vector (type: vec_sum_opnd)
   !! @param:       J      The  total angular momentum (type: vec_sum_opnd)
   !! @param:       L_all  All L from 2 to ndim_sytem-1
   !! @param:       index_L  Integer corresponding to the index of L
   SUBROUTINE get_opL1(L1, J, L_all, index_L)
     type(vec_sum_opnd),      intent(inout)      :: L1
     type(vec_sum_opnd),      intent(in)         :: J
     type(vec_sum_opnd),      intent(in)         :: L_all(:)
     integer,                 intent(in)         :: index_L

     type(vec_sum_opnd)              :: Vtmp

     integer                         :: i
     logical                         :: minus
     character (len=*), parameter :: routine_name='get_opL1'

     if(index_L /= 1) then
       write(out_unitp,*) ' ERROR in',routine_name
       write(out_unitp,*) "The routine can be call only for index_L = 1"
       STOP
     end if

     !write(out_unitp,*) 'get_opL1'

     call copy_F1_into_F2(J, Vtmp)
     !write(out_unitp,*) 'J'
     !CALL write_op(J)
     !write(out_unitp,*) 'END J'

     do i = 2, size(L_all)
       !write(out_unitp,*) 'L(i)',i
       !CALL write_op(L_all(i))
       !write(out_unitp,*) 'END L(i)',i
       call V1_plus_V2_in_Vres(V1 = Vtmp, V2 = L_all(i), Vres = L1, minus = .true.)
       call copy_F1_into_F2(L1, Vtmp)
     end do
     call delete_op(Vtmp)

     !write(out_unitp,*) 'END get_opL1'


   END SUBROUTINE get_opL1

   !! @description: Evaluates the vector L2, the partial angular momentum
   !!               associated to R2
   !! @param:       L      The  vector (type: vec_sum_opnd)
   !! @param:       Fel   an elementary op  which contains the needed
   !!                      information on \theta or u coordinate
   !! @param:       Jz     The z component of the  total angular 
   !!               momentum (type: sum_opnd)
   !! @param:       L_all  All L from 2 to ndim_sytem-1
   !! @param:       index_L  Integer corresponding to the index of L
   SUBROUTINE get_opL2(L, Fel, Jz, Lz_all, index_L, dag)
     type(vec_sum_opnd),      intent(inout)      :: L
     type(opel),              intent(in)         :: Fel
     type(sum_opnd),          intent(in)         :: Jz
     type(sum_opnd),          intent(in)         :: Lz_all(:)
     integer,                 intent(in)         :: index_L
     logical, optional,       intent(in)         :: dag

     type(sum_opnd)                  :: F_sum_tmp

     logical                         :: dag_loc
     integer                         :: i
     logical                         :: minus
     character (len=*), parameter :: routine_name='get_opL2'

     IF (present(dag)) THEN
       dag_loc = dag
     ELSE
       dag_loc = .FALSE.
     END IF


     if(index_L /= 2) then
       write(out_unitp,*) ' ERROR in ',routine_name
       write(out_unitp,*) "The routine can be call only for index_L = 2"
       STOP
     end if

     if(Fel%idf /= 1) then
       write(out_unitp,*) ' ERROR in ',routine_name
       write(out_unitp,*) "The elementary operators should be the Id"
       STOP
     end if 

     call allocate_op(L, 3)

     F_sum_tmp = Jz

     do i = 3, size(Lz_all)
       call get_F1_plus_F2_to_F_sum_nd(F1_sum_nd = F_sum_tmp, F2_sum_nd = Lz_all(i),&
                          & Fres_sum_nd = L%vec_sum(1), minus = .true.)

       F_sum_tmp = L%vec_sum(1)
     end do
     L%vec_sum(3) = F_sum_tmp


     if(Fel%idq == 3 .or. Fel%idq == 7) then

       L%vec_sum(1) = get_cot(Fel) * F_sum_tmp
       L%vec_sum(1)%Cn(:) = -L%vec_sum(1)%Cn(:)

       if (dag_loc) then
         L%vec_sum(2) = get_Pq_dag(Fel)
       else
         L%vec_sum(2) = get_Pq(Fel)
       end if

     else if(Fel%idq == -3 .or. Fel%idq == -7) then
       L%vec_sum(1) = get_cot(Fel) * F_sum_tmp
       L%vec_sum(1)%Cn(:) = -L%vec_sum(1)%Cn(:)

       if (dag_loc) then
         L%vec_sum(2) =  get_Pq(Fel) * get_sin(Fel)
       else
         L%vec_sum(2) =  get_sin(Fel) * get_Pq(Fel)
       end if
       L%vec_sum(2)%Cn(:) = -L%vec_sum(2)%Cn(:)

     else
       write(out_unitp,*) ' ERROR in',routine_name
       write(out_unitp,*) "In Fel, the idq should be 3, -3, 7 or -7"
       STOP
     end if

     call delete_op(F_sum_tmp)


   END SUBROUTINE get_opL2

   !! @description: Evaluates the vector L1, the partial angular momentum
   !!               associated to R1 in the case of subsystem with
   !!               only one euler angle
   !! @param:       L1    The  vector (type: vec_sum_opnd)
   !! @param:       beta  an elementary op  which contains the needed
   !!                      information on \theta (u) or \beta (ub) coordinate
   !! @param:       J_a   The z component of the  total angular 
   !!               momentum (type: sum_opnd)
   SUBROUTINE get_opL1_beta(L1, beta, J_a, dag)
     type(vec_sum_opnd),      intent(inout)      :: L1
     type(opel),              intent(in)         :: beta
     type(sum_opnd),          intent(in)         :: J_a
     logical, optional,       intent(in)         :: dag

     logical :: dag_loc

     character (len=*), parameter :: routine_name='get_opL1_beta'


     IF (present(dag)) THEN
       dag_loc = dag
     ELSE
       dag_loc = .FALSE.
     END IF

     if(beta%idf /= 1 .AND. abs(beta%idq) == 3 ) then
       write(out_unitp,*) ' ERROR in',routine_name
       write(out_unitp,*) "The elementary operators should be the Id"
       STOP
     end if 

     call allocate_op(L1, 3)


     L1%vec_sum(3) = J_a

     L1%vec_sum(1)    = get_cot(beta) * J_a
     L1%vec_sum(1)%Cn = -L1%vec_sum(1)%Cn

     if(beta%idq == 3 .or. beta%idq == 7) then


       if (dag_loc) then
         L1%vec_sum(2) = get_Pq_dag(beta)
       else
         L1%vec_sum(2) = get_Pq(beta)
       end if

     else if(beta%idq == -3 .or. beta%idq == -7) then

       if (dag_loc) then
         L1%vec_sum(2) = get_Pq(beta) * get_sin(beta)
       else
         L1%vec_sum(2) = get_sin(beta) * get_Pq(beta)
       end if
       L1%vec_sum(2)%Cn = -L1%vec_sum(2)%Cn
     else
       write(out_unitp,*) ' ERROR in',routine_name
       write(out_unitp,*) "In beta, the idq should be 3, -3, 7 or -7"
       STOP
     end if

     !write(out_unitp,*) ' L1_beta',dag_loc
     !CALL write_op(L1)
   END SUBROUTINE get_opL1_beta

   !! @description: Defines the conjugate momentum operator
   !! @param:       P          The output  The input function  
   !! @param:       FRel       an elementary op  which contains the needed
   !!                          information on R coordinate
   !! @param:       L          Partial angular momentum
   !! @param:       E          unit vector in spherical coordinates
   SUBROUTINE get_opPi(P, FRel, L, E)

     type(vec_sum_opnd),      intent(inout)      :: P
     type(vec_sum_opnd),      intent(in)         :: L
     type(vec_sum_opnd),      intent(in)         :: E
     type(opel),              intent(in)         :: FRel

     type(vec_sum_opnd)                  :: PRi_Ei 
     type(vec_sum_opnd)                  :: Ei_cross_Li 
     type(vec_sum_opnd)                  :: V_tmp 
     type(sum_opnd)                      :: F_sum_nd

     !logical, parameter           :: debug=.TRUE.
     logical, parameter           :: debug=.FALSE.
     character (len=*), parameter :: routine_name='get_opPi'

     IF (debug) THEN
       write(out_unitp,*) ' BEGINNING ',routine_name
       write(out_unitp,*) ' L'
       CALL write_op(L,header=.TRUE.)
       write(out_unitp,*) ' E'
       CALL write_op(E,header=.TRUE.)
       write(out_unitp,*) ' FRel'
       CALL write_op(FRel,header=.TRUE.)
       CALL flush_perso(out_unitp)
     END IF

     if(FRel%idq /= 2) then
       write(out_unitp,*) ' ERROR in ',routine_name
       write(out_unitp,*) "Data structure of idq of FRel should be 2"
       STOP
     end if


     if(FRel%idf /= 1) then
       write(out_unitp,*) ' ERROR in ',routine_name
       write(out_unitp,*) "The elementary operators should be the Id"
       STOP
     end if 

     F_sum_nd = get_Pq(FRel)   ! PqR

     !PRi_Ei = F_sum_nd * E
     call F_sum_nd_times_V_in_Vres(F_sum_nd = F_sum_nd, V = E, Vres = PRi_Ei)

     call V1_cross_V2_in_Vres(V1 = E, V2 = L, Vres = Ei_cross_Li)


     F_sum_nd = get_Q(FRel,-1)  ! R^-1

     !V_tmp = F_sum_nd * Ei_cross_Li
     call F_sum_nd_times_V_in_Vres(F_sum_nd = F_sum_nd, V = Ei_cross_Li, Vres = V_tmp)

     call V1_plus_V2_in_Vres(V1 = PRi_Ei, V2 = V_tmp, Vres = P, minus = .true.)

    call delete_op(V_tmp)
    call delete_op(Ei_cross_Li)
    call delete_op(PRi_Ei)
    call delete_op(F_sum_nd)

     IF (debug) THEN
       write(out_unitp,*) ' P'
       CALL write_op(P,header=.TRUE.)
       write(out_unitp,*) ' END ',routine_name
       CALL flush_perso(out_unitp)
     END IF

   END SUBROUTINE get_opPi

   !! @description: Defines the adjoint conjugate momentum operator
   !! @param:       P_dag          The output  The input function  
   !! @param:       FRel       an elementary op  which contains the needed
   !!                          information on R coordinate
   !! @param:       L          Partial angular momentum
   !! @param:       E          unit vector in spherical coordinates
   SUBROUTINE get_opPi_dagger(P_dag, FRel, L_dag, E)

     type(vec_sum_opnd),      intent(inout)      :: P_dag
     type(vec_sum_opnd),      intent(in)         :: L_dag
     type(vec_sum_opnd),      intent(in)         :: E
     type(opel),              intent(in)         :: FRel

     type(vec_sum_opnd)                  :: PRi_Ei 
     type(vec_sum_opnd)                  :: Li_cross_Ei 
     type(vec_sum_opnd)                  :: V_tmp 
     type(sum_opnd)                      :: F_sum_nd 

     character (len=*), parameter :: routine_name='get_opPi_dagger'

     IF (FRel%idq /= 2) THEN
       write(out_unitp,*) ' ERROR in',routine_name
       write(out_unitp,*) "Data structure of idq of FRel should be 2"
       STOP
     END IF

     IF (FRel%idf /= 1) THEN
       write(out_unitp,*) ' ERROR in',routine_name
       write(out_unitp,*) "The elementary operators should be the Id"
       STOP
     END IF

     F_sum_nd = get_Pq_dag(FRel)   ! PqR_dag

     call F_sum_nd_times_V_in_Vres(F_sum_nd = F_sum_nd, V = E, Vres = PRi_Ei)
     call V1_cross_V2_in_Vres(V1 = L_dag, V2 = E, Vres = Li_cross_Ei)

     F_sum_nd = get_Q(FRel,-1)     ! R^-1

     call F_sum_nd_times_V_in_Vres(F_sum_nd = F_sum_nd, V = Li_cross_Ei, Vres = V_tmp)
     call V1_plus_V2_in_Vres(V1 = PRi_Ei, V2 = V_tmp, Vres = P_dag)

     call delete_op(V_tmp)
     call delete_op(Li_cross_Ei)
     call delete_op(PRi_Ei)
     call delete_op(F_sum_nd)

   END SUBROUTINE get_opPi_dagger

   !! @description: Defines the spherical unit vectors ei.
   !!               Eq 10, from Ndong et al. J. Chem. Phys. V136, pp034107, 2012
   !! @param:       V      The unit vector (type: vec_sum_opnd)
   !! @param:       theta       an elementary op  which contains the needed
   !!                          information on \theta or u coordinate
   !! @param:       phi       an elementary op  which contains the needed
   !!                          information on \phi  coordinate
   SUBROUTINE get_unit_vector_Ei(V, theta, phi, index_v)
     type(vec_sum_opnd),      intent(inout)      :: V
     type(opel),              intent(in)         :: theta
     type(opel),              intent(in)         :: phi
     integer,                 intent(in)         :: index_v

     character (len = *), parameter :: routine_name='get_unit_vector_Ei'

     call allocate_op(V, 3)

     if (index_v == 1) then

       V%vec_sum(1) = czero
       V%vec_sum(2) = czero
       V%vec_sum(3) = cone

     else if (index_v == 2) then ! it deals automatically idq=3 or -3

       V%vec_sum(1) = get_sin(theta)
       V%vec_sum(2) = get_zero(theta)
       V%vec_sum(3) = get_cos(theta)

     else

       V%vec_sum(1) =  get_sin(theta) * get_cos(phi)
       V%vec_sum(2) =  get_sin(theta) * get_sin(phi)
       V%vec_sum(3) =  get_cos(theta)

     end if

   END SUBROUTINE get_unit_vector_Ei

   !! @description: Determines the Euler matrix rotation in a given frame (not used!!)
   !! @param:       F_system   The  data structure of the subsystem (type: TYpe_BFtransfo)
   !! @param:       mat_R   The  matrix (type: sum_opnd)
   !! @param:       mat_RTranspo   The  transpose of Mat_R (type: sum_opnd)
   SUBROUTINE project_Pi(F_system, Pi, Pi_dag)
     USE mod_BunchPolyTransfo, only : Type_BFTransfo

     type(Type_BFtransfo),     intent(inout)      :: F_system
     type(vec_sum_opnd),       intent(inout)      :: Pi
     type(vec_sum_opnd),       intent(inout)      :: Pi_dag

     type(sum_opnd), allocatable         :: D_a(:,:)
     type(sum_opnd), allocatable         :: D_b(:,:)
     type(sum_opnd), allocatable         :: D_g(:,:)
     type(sum_opnd), allocatable         :: D_trans(:,:)
     type(vec_sum_opnd)                  :: Pi_tmp
     type(opel)                          :: falpha,fbeta,fgamma

     integer                         :: error
     character (len=*), parameter    :: routine_name='project_Pi'


     if(F_system%euler(3)) then

       fgamma = F_system%QEuler(3)

       D_g = get_MatRotz(fgamma) ! gamma

       call M_opnd_times_V_in_Vres(D_g, Pi, Pi_tmp)
       call copy_F1_into_F2(Pi_tmp, Pi)

       D_trans = Transpose_Mat_OF_sum_opnd(D_g)

       call V_times_M_opnd_in_Vres(Pi_dag, D_trans, Pi_tmp)
       call copy_F1_into_F2(Pi_tmp, Pi_dag)
     end if

     if(F_system%euler(2)) then

       fbeta = F_system%QEuler(2)

       D_b = get_MatRoty(fbeta) ! beta

       call M_opnd_times_V_in_Vres(D_b, Pi, Pi_tmp)
       call copy_F1_into_F2(Pi_tmp, Pi)

       D_trans = Transpose_Mat_OF_sum_opnd(D_b)

       call V_times_M_opnd_in_Vres(Pi_dag, D_trans, Pi_tmp)
       call copy_F1_into_F2(Pi_tmp, Pi_dag)
     end if
     if(F_system%euler(1)) then

       falpha = F_system%QEuler(1)

       D_a = get_MatRotz(falpha) ! alpha

       call M_opnd_times_V_in_Vres(D_a, Pi, Pi_tmp)
       call copy_F1_into_F2(Pi_tmp, Pi)

       D_trans = Transpose_Mat_OF_sum_opnd(D_b)

       call V_times_M_opnd_in_Vres(Pi_dag, D_trans, Pi_tmp)
       call copy_F1_into_F2(Pi_tmp, Pi_dag)
     end if

     CALL dealloc_NParray(D_a,'D_a',routine_name)
     CALL dealloc_NParray(D_b,'D_b',routine_name)
     CALL dealloc_NParray(D_g,'D_g',routine_name)
     CALL dealloc_NParray(D_trans,'D_trans',routine_name)
     call delete_op(Pi_tmp)

   END SUBROUTINE project_Pi

   

   !! @description: Determines the Euler matrix rotation in a given frame
   !! @param:       F_system   The  data structure of the subsystem (type: TYpe_BFtransfo)
   !! @param:       mat_R   The  matrix (type: sum_opnd)
   !! @param:       mat_RTranspo   The  transpose of Mat_R (type: sum_opnd)
   SUBROUTINE get_Euler_MatRot(F_system, Mat_R, Mat_RTranspo)
     USE mod_BunchPolyTransfo, only : Type_BFTransfo

     type(Type_BFtransfo),         intent(inout)      :: F_system
     type(sum_opnd), allocatable,  intent(inout)      :: mat_R(:,:)
     type(sum_opnd), allocatable,  intent(inout)      :: mat_RTranspo(:,:)

     type(sum_opnd), allocatable         :: M_opnd(:,:)
     type(sum_opnd), allocatable         :: D_a(:,:)
     type(sum_opnd), allocatable         :: D_b(:,:)
     type(sum_opnd), allocatable         :: D_g(:,:)

     type(opel)                      :: falpha
     type(opel)                      :: fbeta
     type(opel)                      :: fgamma

     integer                         :: error
     character (len=*), parameter :: routine_name = 'get_Euler_MatRot'


     ! TTF
     if (compare_tab(F_system%euler, (/.true., .true., .false./))) then

       falpha = F_system%QEuler(1)
       fbeta  = F_system%QEuler(2)

       D_a = get_MatRotz(falpha) ! alpha
       D_b = get_MatRoty(fbeta) ! beta

       CALL alloc_NParray(Mat_R,(/3,3/),'Mat_R',routine_name)
       CALL M1_times_M2_in_Mres(D_a, D_b, Mat_R)

     ! FTT
     else if (compare_tab(F_system%euler, (/.false., .true., .true./))) then

       fbeta  = F_system%QEuler(2)
       fgamma = F_system%QEuler(3)

       D_b = get_MatRoty(fbeta) ! beta
       D_g = get_MatRotz(fgamma) ! gamma

       CALL alloc_NParray(Mat_R,(/3,3/),'Mat_R',routine_name)

       call M1_times_M2_in_Mres(D_b, D_g, Mat_R)

     !TTT
     else if (compare_tab(F_system%euler, (/.true., .true., .true./))) then

       falpha = F_system%QEuler(1)
       fbeta  = F_system%QEuler(2)
       fgamma = F_system%QEuler(3)

       D_a = get_MatRotz(falpha) ! alpha
       D_b = get_MatRoty(fbeta) ! beta
       D_g = get_MatRotz(fgamma) ! gamma


       CALL alloc_NParray(M_opnd,(/3,3/),'M_opnd',routine_name)
       call M1_times_M2_in_Mres(D_b, D_g, M_opnd)

       CALL alloc_NParray(Mat_R,(/3,3/),'Mat_R',routine_name)
       call M1_times_M2_in_Mres(D_a, M_opnd, Mat_R)

     end if

     CALL dealloc_NParray(M_opnd,'M_opnd',routine_name)
     CALL dealloc_NParray(D_a,'D_a',routine_name)
     CALL dealloc_NParray(D_b,'D_b',routine_name)
     CALL dealloc_NParray(D_g,'D_g',routine_name)



     Mat_RTranspo = Transpose_Mat_OF_sum_opnd(Mat_R)

     
   END SUBROUTINE get_Euler_MatRot
   FUNCTION get_MatRotz(angle) result(mat_R)
     type(opel), intent(in)                       :: angle

     type(sum_opnd), allocatable                  :: mat_R(:,:)

     character (len=*), parameter :: routine_name = 'get_MatRotz'


       CALL alloc_NParray(mat_R,(/3,3/),'mat_R',routine_name)

       mat_R(1,1) = get_cos(angle)
       mat_R(1,2) = get_sin(angle)
       mat_R(1,2)%Cn(1) = -CONE
       mat_R(1,3) = CZERO

       mat_R(2,1) = get_sin(angle)
       mat_R(2,2) = get_cos(angle)
       mat_R(2,3) = CZERO

       mat_R(3,1) = CZERO
       mat_R(3,2) = CZERO
       mat_R(3,3) = CONE

   END FUNCTION get_MatRotz
   FUNCTION get_MatRoty(angle) result(mat_R)
     type(opel), intent(in)                       :: angle

     type(sum_opnd), allocatable                  :: mat_R(:,:)

     character (len=*), parameter :: routine_name = 'get_MatRoty'


       CALL alloc_NParray(mat_R,(/3,3/),'mat_R',routine_name)

       mat_R(1,1) = get_cos(angle)
       mat_R(1,2) = CZERO
       mat_R(1,3) = get_sin(angle)

       mat_R(2,1) = CZERO
       mat_R(2,2) = CONE
       mat_R(2,3) = CZERO

       mat_R(3,1) = get_sin(angle)
       mat_R(3,1)%Cn(1) = -CONE
       mat_R(3,2) = CZERO
       mat_R(3,3) = get_cos(angle)

   END FUNCTION get_MatRoty


   !! @description: get the conjugate mommentum associated with the vectors
   !!               which belong to the system F_system  including those its subsystems
   !! @param:     F_system       The  data structure of the system (type: TYpe_BFtransfo)
   !! @param:     P_Euler        The  data structure of the conjugate momenta
   !! @param:     Pi_BF          Array of conjugate momenta of the subsystem (type: vec_sum_opnd)
   !! @param:     Pi_dag_BF      Array of the  adjoint of Pi_BF      (type: vec_sum_opnd)  
   !! @param:     zero_Pi_BF     Logical array to check if Pi_BF is set to zero or not
   RECURSIVE SUBROUTINE get_Pi_subsyst(F_system, P_Euler, Pi_BF, Pi_dag_BF,&
                                       zero_Pi_BF)
     USE mod_BunchPolyTransfo, only : Type_BFTransfo


     type(Type_BFtransfo),         intent(in)      :: F_system
     type(Type_PiEulerRot),        intent(in)      :: P_Euler(:)
     type(vec_sum_opnd),           intent(inout)   :: Pi_BF(:)
     type(vec_sum_opnd),           intent(inout)   :: Pi_dag_BF(:)
     logical,                      intent(inout)   :: zero_Pi_BF(:)


     integer                                     :: nvec, nsub_syst
     integer                                     :: i, j, n, error, n_size
     logical                                     :: true_BF
     character (len=*), parameter  :: routine_name='get_Pi_subsyst'

     nsub_syst = 0
     do i=1, F_system%nb_vect
       if(F_system%tab_BFTransfo(i)%frame) nsub_syst = nsub_syst+1
     end do
     nvec = F_system%nb_vect-nsub_syst+1

     if (compare_tab(F_system%euler, (/.false., .false., .false./))) then
       true_BF = .true.
     else
       true_BF = .false. 
     end if

     if(.not.true_BF) then !project onto the frame of the container 
       do i = 1, nvec
         call copy_F1_into_F2(P_Euler(F_system%listVFr(i))%Pi, &
         &                           Pi_BF(F_system%listVFr(i)))
         call copy_F1_into_F2(P_Euler(F_system%listVFr(i))%Pidag,&
         &                    Pi_dag_BF(F_system%listVFr(i)))
         zero_Pi_BF(F_system%listVFr(i)) = .false.
       end do
     end if

     do i=1, F_system%nb_vect
       if(F_system%tab_BFTransfo(i)%frame) then
         call get_Pi_subsyst(F_system%tab_BFTransfo(i), P_Euler, Pi_BF, &
         &                                             Pi_dag_BF, zero_Pi_BF)
       end if
     end do 

   END SUBROUTINE get_Pi_subsyst

   !! @description: Defines the total angular momentum in terms of
   !! @param:     F_system       The  data structure of the subsystem (type: TYpe_BFtransfo)
   !! @param:     mat_R          The  matrix (type: sum_opnd)
   !! @param:     mat_RTranspo   The  transpose of Mat_R (type: sum_opnd)
   !! @param:     Pi             The  conjugate momenta of the subsystem (type: vec_sum_opnd)
   !! @param:     Pi_dag         The  adjoint of Pi      (type: vec_sum_opnd)  
   RECURSIVE SUBROUTINE Mat_Rot_times_Pi_subsyst(F_system, Mat_R, Mat_RTranspo, &
   &                                             P_Euler)
     USE mod_BunchPolyTransfo, only : Type_BFTransfo

     type(Type_BFtransfo),         intent(inout)      :: F_system
     type(sum_opnd),               intent(inout)      :: mat_R(:,:)
     type(sum_opnd),               intent(inout)      :: mat_RTranspo(:,:)
     type(Type_PiEulerRot),        intent(inout)      :: P_Euler(:)

     type(vec_sum_opnd)                          :: Pi_BF
     type(vec_sum_opnd)                          :: Pi_dag_BF
     integer                                     :: nvec, nsub_syst
     integer                                     :: i, j, n, error, n_size
     logical                                     :: true_BF
     character (len =*), parameter  :: routine_name='Mat_Rot_times_Pi_subsyst'


     nsub_syst = 0
     do i=1, F_system%nb_vect
       if(F_system%tab_BFTransfo(i)%frame) nsub_syst = nsub_syst+1
     end do
     nvec = F_system%nb_vect-nsub_syst+1
     n = size(P_Euler)

     if (compare_tab(F_system%euler, (/.false., .false., .false./))) then
       true_BF = .true.
     else
       true_BF = .false. 
     end if

     if(.not.true_BF) then !project onto the frame of the container 
       do i = 1, nvec
         call M_opnd_times_V_in_Vres(Mat_R, P_Euler(F_system%listVFr(i))%Pi, Pi_BF)
         call copy_F1_into_F2(Pi_BF, P_Euler(F_system%listVFr(i))%Pi)
         call V_times_M_opnd_in_Vres(P_Euler(F_system%listVFr(i))%Pidag, Mat_RTranspo, Pi_dag_BF)
         call copy_F1_into_F2(Pi_dag_BF, P_Euler(F_system%listVFr(i))%Pidag)
       end do
       call delete_op(Pi_BF)
       call delete_op(Pi_dag_BF)
     else
       return
     end if
     do i=1, F_system%nb_vect
       if(F_system%tab_BFTransfo(i)%frame) then
         call Mat_Rot_times_Pi_subsyst(F_system%tab_BFTransfo(i), Mat_R, Mat_RTranspo, &
         &                                             P_Euler)
       end if
     end do 

   END SUBROUTINE Mat_Rot_times_Pi_subsyst

   SUBROUTINE get_Lz_F_system_parent(F_system_parent, Liz)
     USE mod_BunchPolyTransfo, only : Type_BFTransfo

     type(Type_BFtransfo),        intent(in)          :: F_system_parent
     type(sum_opnd),              intent(inout)       :: Liz(:)

     type(opel)                  :: phi,x,y
     integer                     :: i, j
     character (len=*), parameter :: routine_name='get_Lz_F_system_parent'


      if(.not.F_system_parent%tab_BFTransfo(1)%frame) then
        write(out_unitp,*) ' ERROR in',routine_name
        write(out_unitp,*) "This routine can be called only "
        write(out_unitp,*) "  if frame=T for the first vector"
        STOP
      end if

       j = 0
       do i = 1, F_system_parent%nb_vect
         if(.not.F_system_parent%tab_BFTransfo(i)%frame) then
           j=j+1
           if(F_system_parent%tab_BFTransfo(i)%cart) then

             x = F_system_parent%tab_BFTransfo(i)%Qvec(1)
             y = F_system_parent%tab_BFTransfo(i)%Qvec(2)

             Liz(j) = get_Q(x) * get_Pq(y)

             CALL F1_nd_MINUS_TO_Fres_sum_nd( get_Q(y) * get_Pq(x) ,Liz(j))

           else

             phi = F_system_parent%tab_BFTransfo(i)%Qvec(3)

             Liz(j) = get_Pq(phi)

           end if
         end if
       end do

   END  SUBROUTINE get_Lz_F_system_parent


   !! @description: Defines the KEO
   !! @param:       TWOxKEO The output
   !! @param:       nvec    Integer corresponding to the size of the system
   RECURSIVE SUBROUTINE get_opKEO(F_system, TWOxKEO, P_Euler, M_mass_out, &
                                  scalar_PiPj, F_system_parent)
     USE mod_BunchPolyTransfo, only : Type_BFTransfo

     type(Type_BFtransfo),           intent(inout)       :: F_system
     type(sum_opnd),                 intent(inout)       :: TWOxKEO
     type(Type_PiEulerRot),          intent(inout)       :: P_Euler(:)
     type(sum_opnd),                 intent(in)          :: M_mass_out(:,:)
     logical,                        intent(inout)       :: scalar_PiPj(:,:)
     type(Type_BFtransfo), optional, intent(in)          :: F_system_parent

     integer                         :: i, j
     integer                         :: i_syst, j_syst
     integer                         :: nsub_syst
     integer                         :: index_tmp, i_v1_ref


     !logical, parameter           :: debug=.TRUE.
     logical, parameter           :: debug=.FALSE.
     character (len = *), parameter :: routine_name= 'get_opKEO'


     IF (debug) THEN
       write(out_unitp,*) 'BEGINNING ',routine_name
     END IF

     write(out_unitp,*) 'entree S_(', F_system%tab_num_frame,')'
     CALL flush_perso(out_unitp)
     nsub_syst = 0
     do i=1, F_system%nb_vect
       if(F_system%tab_BFTransfo(i)%frame) nsub_syst = nsub_syst+1
     end do

     do i_syst = F_system%nb_vect, 1, -1
       if(F_system%tab_BFTransfo(i_syst)%frame) then
         call get_opKEO(F_system%tab_BFTransfo(i_syst), TWOxKEO, P_Euler,   &
                        M_mass_out, scalar_PiPj, F_system)
       end if
     end do

     write(out_unitp,*) 'sub_system S_(', F_system%tab_num_frame,')'
     CALL flush_perso(out_unitp)

     if(F_system%frame) then

       if (compare_tab(F_system%euler, (/.false., .true., .true./))) then

         IF (present(F_system_parent)) THEN
           call get_opKEO_subsyst_2euler(F_system, P_Euler, M_mass_out,   &
                                       scalar_PiPj, F_system_parent)
         ELSE
           STOP 'ERROR in get_opKEO: F_system_parent is absent'
         END IF

       ELSE

         IF (F_system%nb_vect_tot == 1) THEN

           IF (present(F_system_parent)) THEN
             CALL get_opKEO_subsyst_nvectot1(F_system, P_Euler,M_mass_out,&
                                             scalar_PiPj,F_system_parent)
           ELSE
             STOP 'ERROR in get_opKEO: F_system_parent is absent'
           END IF

         ELSE
           CALL get_opKEO_subsyst(F_system, P_Euler,M_mass_out,scalar_PiPj)
!           IF (present(F_system_parent)) THEN
!             CALL get_opKEO_subsyst(F_system, P_Euler,M_mass_out,       &
!                                    scalar_PiPj,F_system_parent)
!           ELSE
!             STOP 'ERROR in get_opKEO: F_system_parent is absent'
!           END IF

         END IF
       END IF

       write(out_unitp,*) 'system S_r, r=', F_system%tab_num_frame
       write(out_unitp,*) 'size before simplify', size(F_system%KEO%sum_prod_op1d)
       CALL flush_perso(out_unitp)

       call Simplify_Sum_OpnD(F_system%KEO,Expand_Sin2=.TRUE.)

       write(out_unitp,*) 'size after simplify', size(F_system%KEO%sum_prod_op1d)
       CALL flush_perso(out_unitp)
       CALL F1_sum_nd_PLUS_TO_Fres_sum_nd(F_system%KEO,TWOxKEO)

     end if
     !call Simplify_Sum_OpnD(TWOxKEO,Expand_Sin2=.TRUE.)

     IF (debug) THEN
       CALL write_op(F_system%keo)
       write(out_unitp,*) 'END ',routine_name
     END IF

   END SUBROUTINE get_opKEO


   !! @param:       KOE          The output 
   !! @param:       nvec    Integer corresponding to the size of the system
   !! @param:       M_mass  The mass matrix
   !RECURSIVE SUBROUTINE get_opKEO_subsyst(F_system, P_Euler, M_mass_out, scalar_PiPj, &
   !                                       F_system_parent)
   RECURSIVE SUBROUTINE get_opKEO_subsyst(F_system, P_Euler, M_mass_out, scalar_PiPj)
     USE mod_BunchPolyTransfo, only : Type_BFTransfo

     type(Type_BFtransfo),                 intent(inout)      :: F_system
     type(Type_PiEulerRot),                intent(inout)      :: P_Euler(:)
     type(sum_opnd),                       intent(in)         :: M_mass_out(:,:)
     logical,                              intent(inout)      :: scalar_PiPj(:,:)
     !type(Type_BFtransfo),                 intent(in)         :: F_system_parent

     type(sum_opnd)                  :: RR_inv,MPRPR,PRPR,PiPi,MPiPi
     type(sum_opnd)                  :: LiLi,JLi,LiJ,JJ,L1L1


     type(vec_sum_opnd), pointer         :: Pi_BF(:)
     type(vec_sum_opnd), pointer         :: Pi_dag_BF(:)
     type(vec_sum_opnd), pointer         :: L(:)
     type(vec_sum_opnd), pointer         :: L_dag(:)
     type(sum_opnd), pointer             :: Lz(:)
     type(sum_opnd), pointer             :: Liz_parent(:)
     type(vec_sum_opnd)                  :: L1_bis
     type(vec_sum_opnd)                  :: L1dag_bis
     type(sum_opnd),     allocatable     :: Mat_R(:,:)
     type(sum_opnd),     allocatable     :: Mat_RTranspo(:,:)
     type(sum_opnd),     allocatable     :: Mres(:,:)

     type(vec_sum_opnd)                  :: V1_tmp
     type(vec_sum_opnd)                  :: V2_tmp

     type(vec_sum_opnd)                  :: V_sum_J
     type(vec_sum_opnd)                  :: V_sum_Jdag
     type(sum_opnd)                      :: keo_Pi

     type(sum_opnd)                      :: Ja_sum_subsyst

     type(opel)                          :: Ja_el
     type(opel)                          :: Jb_el
     type(opel)                          :: Jg_el

     integer                         :: index_L
     integer                         :: i, j, k1, k2
     integer                         :: k3, k4, k, n, iv
     integer                         :: nvec
     integer                         :: nsub_syst, nvec_tot
     integer                         :: nvec_parent, nsub_syst_parent
     integer                         :: error
     logical                         :: true_BF
     logical                         :: alloc_subsyst
     logical                         :: euler
     logical, pointer                :: zero_Pi_BF(:)
     real(kind=Rkind)                :: sum_op, prod_op

     !logical, parameter           :: debug=.TRUE.
     logical, parameter           :: debug=.FALSE.
     character (len = *), parameter :: routine_name= 'get_opKEO_subsyst'

     IF (debug) THEN
       write(out_unitp,*) 'BEGINNING ',routine_name
     END IF

     nullify(Pi_BF)
     nullify(Pi_dag_BF)
     nullify(L)
     nullify(L_dag)
     nullify(Lz)
     nullify(Liz_parent)

     nullify(zero_Pi_BF)

     nsub_syst = 0
     do i=1, F_system%nb_vect
       if(F_system%tab_BFTransfo(i)%frame) nsub_syst = nsub_syst+1
     end do
     nvec = F_system%nb_vect-nsub_syst+1
     nvec_tot = F_system%nb_vect_tot
     if (compare_tab(F_system%euler, (/.false., .false., .false./))) then
       true_BF = .true.
     else
       true_BF = .false.
     end if

     do i = 1, nvec
       if(.not.associated(P_Euler(F_system%listVFr(i))%Tab_num_Frame)) then

         CALL alloc_array(P_Euler(F_system%listVFr(i))%Tab_num_Frame,shape(F_system%tab_num_Frame), &
                         'P_Euler(F_system%listVFr(i))%Tab_num_Frame',routine_name)
         P_Euler(F_system%listVFr(i))%Tab_num_Frame(:) = F_system%Tab_num_Frame(:)
       end if
     end do

     !!allocation
     CALL alloc_array(L,    (/nvec/),'L',routine_name)
     CALL alloc_array(L_dag,(/nvec/),'L_dag',routine_name)
     if(nsub_syst > 0) then
       if(.not.F_system%tab_BFTransfo(1)%frame) then
         CALL alloc_array(Lz,(/nvec+nsub_syst/),'Lz',routine_name)
       else
         CALL alloc_array(Lz,(/nvec/),'Lz',routine_name)
       end if
     else
       CALL alloc_array(Lz,(/nvec/),'Lz',routine_name)
     end if
     CALL alloc_array(Pi_BF,    shape(P_Euler),'Pi_BF',routine_name)
     CALL alloc_array(Pi_dag_BF,shape(P_Euler),'Pi_dag_BF',routine_name)
     CALL alloc_array(zero_Pi_BF,shape(P_Euler),'zero_Pi_BF',routine_name)

     ! Initiliation of the elementary operators, J and Jdag
     write(out_unitp,*) 'init elementaries op. for S_(',F_system%tab_num_frame,')'

     Ja_el = F_system%QEuler(1) ! alpha
     Jb_el = F_system%QEuler(2) ! beta or cos(beta)
     Jg_el = F_system%QEuler(3) ! gamma

     if(true_BF) then

       call allocate_op(F_system%J,    3)
       call allocate_op(F_system%Jdag, 3)

       F_system%J%vec_sum(1) = get_Jx(Ja_el)
       F_system%J%vec_sum(2) = get_Jy(Jb_el)
       F_system%J%vec_sum(3) = get_Jz(Jg_el)

       F_system%Jdag%vec_sum(1) = get_Jx(Ja_el)
       F_system%Jdag%vec_sum(2) = get_Jy(Jb_el)
       F_system%Jdag%vec_sum(3) = get_Jz(Jg_el)

     else
       call get_opJ_projected_into_BF(F_system%J,    Ja_el, Jb_el, F3el = Jg_el)
       call get_opJ_projected_into_BF(F_system%Jdag, Ja_el, Jb_el, F3el = Jg_el, dag = .true.)
     end if


     do i = 1, size(Lz)
       Lz(i) = czero
     end do

     !get the corresponding Euler matrix rotation
     do i = 1, size(Pi_BF)
       call allocate_op(Pi_BF(i),3)
       CALL zero_TO_vec_sum_opnd(Pi_BF(i))
       call allocate_op(Pi_dag_BF(i),3)
       CALL zero_TO_vec_sum_opnd(Pi_dag_BF(i))

       zero_Pi_BF(i) = .true.

     end do
     if(.not.true_BF) then
       call get_Euler_MatRot(F_system, Mat_R, Mat_RTranspo)
     end if

     ! Initiliation of the Li (i>2) operators
     !computation of the Pi, i>2 in their local frame
     call allocate_op(V1_tmp, 3)
     CALL zero_TO_vec_sum_opnd(V1_tmp)

     call allocate_op(V2_tmp, 3)
     CALL zero_TO_vec_sum_opnd(V2_tmp)


     iv = 1
     do i = 1, F_system%nb_vect
       if(.not.F_system%tab_BFTransfo(i)%frame) then
         iv=iv+1
         if(i>1) then
           if(F_system%tab_BFTransfo(i)%cart) then

             V1_tmp%vec_sum(1) = get_Q(F_system%tab_BFTransfo(i)%Qvec(1)) !x
             V1_tmp%vec_sum(2) = get_Q(F_system%tab_BFTransfo(i)%Qvec(2)) !y
             V1_tmp%vec_sum(3) = get_Q(F_system%tab_BFTransfo(i)%Qvec(3)) !z

             Pi_BF(F_system%listVFr(iv))%vec_sum(1) = get_Pq(F_system%tab_BFTransfo(i)%Qvec(1)) !Px
             Pi_BF(F_system%listVFr(iv))%vec_sum(2) = get_Pq(F_system%tab_BFTransfo(i)%Qvec(2)) !Py
             Pi_BF(F_system%listVFr(iv))%vec_sum(3) = get_Pq(F_system%tab_BFTransfo(i)%Qvec(3)) !Pz

             call V1_cross_V2_in_Vres(V1_tmp, Pi_BF(F_system%listVFr(iv)), L(iv))
             call V1_cross_V2_in_Vres(Pi_BF(F_system%listVFr(iv)), V1_tmp, L_dag(iv))

             ! L_dag(iv) = -L_dag(iv)
             call V1_plus_V2_in_Vres(V2_tmp, L_dag(iv), V1_tmp, minus=.true.)
             call copy_F1_into_F2(V1_tmp,  L_dag(iv))

             call copy_F1_into_F2(Pi_BF(F_system%listVFr(iv)), Pi_dag_BF(F_system%listVFr(iv)))

           else
             IF (iv == 2) THEN
               index_L = iv+1
             ELSE
               index_L = iv
             END IF

             call get_opLi(L = L(iv),     theta = F_system%tab_BFTransfo(i)%Qvec(2), &
                                          phi   = F_system%tab_BFTransfo(i)%Qvec(3), &
                                  index_L = index_L, Li=F_system%tab_BFTransfo(i)%Li)

             call get_opLi(L = L_dag(iv), theta = F_system%tab_BFTransfo(i)%Qvec(2), &
                                          phi   = F_system%tab_BFTransfo(i)%Qvec(3), &
                     index_L = index_L, dag = .true., Li=F_system%tab_BFTransfo(i)%Li)

             call get_opPi(       Pi_BF(F_system%listVFr(iv)),                       &
                           FRel=F_system%tab_BFTransfo(i)%Qvec(1), L=L(iv),          &
                           E=F_system%tab_BFTransfo(i)%Unit_Vector)
             call get_opPi_dagger(Pi_dag_BF(F_system%listVFr(iv)),                   &
                           FRel=F_system%tab_BFTransfo(i)%Qvec(1), L_dag=L_dag(iv),  &
                           E=F_system%tab_BFTransfo(i)%Unit_Vector)

           end if
           Lz(iv) = L(iv)%vec_sum(3)
           zero_Pi_BF(F_system%listVFr(iv)) = .false.
         end if
       end if
     end do

     if(nsub_syst>0) then
       if(.not.F_system%tab_BFTransfo(1)%frame) then
         j = 0
         do i = 1, F_system%nb_vect
           if(F_system%tab_BFTransfo(i)%frame) then
             j = j+1
             Lz(nvec+j) = get_Pq(F_system%tab_BFTransfo(i)%QEuler(1)) ! alpha Eq A7 (Tana2012)
           end if
         end do
       end if
     end if

     if(F_system%nb_vect >0) then
       if(.not.F_system%tab_BFTransfo(1)%frame) then

         call get_opL2(L = L(2),     Fel = F_system%tab_BFTransfo(1)%Qvec(2), &
                      Jz = F_system%J%vec_sum(3), Lz_all = Lz, index_L = 2)

         call get_opL2(L = L_dag(2), Fel = F_system%tab_BFTransfo(1)%Qvec(2), &
                      Jz = F_system%Jdag%vec_sum(3), Lz_all = Lz, index_L = 2, dag = .true.)

         call get_opPi(       Pi_BF(F_system%listVFr(2)),     F_system%tab_BFTransfo(1)%Qvec(1), L(2),    &
                              E=F_system%tab_BFTransfo(1)%Unit_Vector)
         call get_opPi_dagger(Pi_dag_BF(F_system%listVFr(2)), F_system%tab_BFTransfo(1)%Qvec(1), L_dag(2),&
                              E=F_system%tab_BFTransfo(1)%Unit_Vector)

         zero_Pi_BF(F_system%listVFr(2)) = .false.
       end if
     end if

     if(nvec == 1) then
       CALL copy_F1_into_F2(F_system%J,    L(1))
       CALL copy_F1_into_F2(F_system%Jdag, L_dag(1))
     else
       call get_opL1(L1 = L(1),     J = F_system%J,    L_all = L,     index_L = 1)
       call get_opL1(L1 = L_dag(1), J = F_system%Jdag, L_all = L_dag, index_L = 1)
     end if

     CALL copy_F1_into_F2(L(1),     L1_bis)
     CALL copy_F1_into_F2(L_dag(1), L1dag_bis)

     CALL allocate_op(V_sum_J, 3)
     CALL zero_TO_vec_sum_opnd(V_sum_J)
     CALL allocate_op(V_sum_Jdag, 3)
     CALL zero_TO_vec_sum_opnd(V_sum_Jdag)

     if(nsub_syst > 0) then
       do i = 1, F_system%nb_vect
         if(F_system%tab_BFTransfo(i)%frame) then
           CALL V1_PLUS_TO_Vres(F_system%tab_BFTransfo(i)%J, V_sum_J)
           CALL V1_PLUS_TO_Vres(F_system%tab_BFTransfo(i)%Jdag, V_sum_Jdag)
         end if
       end do
       !L1_bis = L1_bis - V_sum_J
       CALL V1_MINUS_TO_Vres(V_sum_J,    L1_bis)
       CALL V1_MINUS_TO_Vres(V_sum_Jdag, L1dag_bis)

     end if

     call get_opPi(       Pi_BF(F_system%listVFr(1)),     F_system%Qvec(1), L1_bis,   &
                          E=F_system%Unit_Vector)
     call get_opPi_dagger(Pi_dag_BF(F_system%listVFr(1)), F_system%Qvec(1), L1dag_bis,&
                          E=F_system%Unit_Vector)

     zero_Pi_BF(F_system%listVFr(1)) = .false.
     do i = 1, F_system%nb_vect
       if(F_system%tab_BFTransfo(i)%frame) then
         call get_Pi_subsyst(F_system%tab_BFTransfo(i), P_Euler, Pi_BF, &
         &                   Pi_dag_BF, zero_Pi_BF)
       end if
     end do

     !Coupling terms between subsystems
     keo_Pi = czero
     do i = 1, size(Pi_BF)
     do j = 1, size(Pi_BF)
         if(i/=j  .and. abs(M_mass_out(i,j)%Cn(1)) >1.0e-13_Rkind &
           &      .and. .not.zero_Pi_BF(i) &
           &      .and. .not.zero_Pi_BF(j)) then
           if(.not.compare_tab(P_euler(i)%tab_num_frame, &
             &       P_euler(j)%tab_num_frame) .and. &
             &       .not.scalar_PiPj(i,j)) then

             call  V1_scalar_V2_in_F_sum_nd(Pi_dag_BF(i),Pi_BF(j), PiPi)
             MPiPi = M_mass_out(i,j) * PiPi

             CALL F1_sum_nd_PLUS_TO_Fres_sum_nd(MPiPi,keo_Pi)

             scalar_PiPj(i,j) = .true.
           end if

         end if
     end do
     end do
     if(.not.true_BF) then !project onto the frame of the container
       do i = 1, nvec
         call M_opnd_times_V_in_Vres(Mat_R, Pi_BF(F_system%listVFr(i)), &
         & P_Euler(F_system%listVFr(i))%Pi)
         call V_times_M_opnd_in_Vres(Pi_dag_BF(F_system%listVFr(i)), &
         &                Mat_RTranspo, P_Euler(F_system%listVFr(i))%Pidag)
       end do
       do i = 1, F_system%nb_vect
         if(F_system%tab_BFTransfo(i)%frame) then
           call Mat_Rot_times_Pi_subsyst(F_system%tab_BFTransfo(i), &
           &                   Mat_R, Mat_RTranspo, P_Euler)
         end if
       end do
     else
       do i = 1, nvec
         call copy_F1_into_F2(Pi_BF(F_system%listVFr(i)), &
         &  P_Euler(F_system%listVFr(i))%Pi)
         call copy_F1_into_F2(Pi_dag_BF(F_system%listVFr(i)), &
         &  P_Euler(F_system%listVFr(i))%Pidag)
       end do
     end if

     do i = 1, nvec
     do j = 1, nvec
       scalar_PiPj(F_system%listVFr(i), F_system%listVFr(j)) = .true.
     end do
     end do
     CALL dealloc_array(zero_Pi_BF,'zero_Pi_BF',routine_name)

     CALL dealloc_NParray(Mat_R,'Mat_R',routine_name)
     CALL dealloc_NParray(Mat_RTranspo,'Mat_RTranspo',routine_name)

     !!! Starting the computation: Diagonal terms
     L1L1         = czero ! L1L1
     F_system%KEO = keo_Pi ! keo

     iv = 1
     do i = 1, F_system%nb_vect
       if(.not.F_system%tab_BFTransfo(i)%frame) then
         iv=iv+1
         if(i>1) then

           if(F_system%tab_BFTransfo(i)%cart) then

             call  V1_scalar_V2_in_F_sum_nd(Pi_dag_BF(F_system%listVFr(iv)), &
                                            Pi_BF(F_system%listVFr(iv)), PiPi)

             call  V1_scalar_V2_in_F_sum_nd(L_dag(iv),L(iv), LiLi)
             CALL F1_sum_nd_PLUS_TO_Fres_sum_nd(LiLi,L1L1)

           else

             call Li_scalar_Li_from_Eq75(theta = F_system%tab_BFTransfo(i)%Qvec(2), &
             &                           phi   = F_system%tab_BFTransfo(i)%Qvec(3), &
             &                           LiLi  =  LiLi)
             CALL F1_sum_nd_PLUS_TO_Fres_sum_nd(LiLi,L1L1)


             PiPi = get_Pq_dag (F_system%tab_BFTransfo(i)%Qvec(1)) *       &
                    get_Pq(F_system%tab_BFTransfo(i)%Qvec(1)) ! PR^dag * PR

             RR_inv = get_Q(F_system%tab_BFTransfo(i)%Qvec(1),-2) !1/R^2

             CALL F1_sum_nd_PLUS_TO_Fres_sum_nd( RR_inv*LiLi ,PiPi)

           end if
           MPiPi = F_system%M_mass(iv,iv) * PiPi
           CALL F1_sum_nd_PLUS_TO_Fres_sum_nd(MPiPi,F_system%KEO)

         end if
       end if
     end do

     !!L2dag_L2
     if(F_system%nb_vect >0) then
       if(.not. F_system%tab_BFTransfo(1)%frame) then

         CALL V1_scalar_V2_in_F_sum_nd(L_dag(2),L(2),LiLi)
         CALL F1_sum_nd_PLUS_TO_Fres_sum_nd(LiLi,L1L1)

         PiPi = get_Pq_dag (F_system%tab_BFTransfo(1)%Qvec(1)) *       &
                get_Pq(F_system%tab_BFTransfo(1)%Qvec(1)) ! PR^dag * PR

         RR_inv = get_Q(F_system%tab_BFTransfo(1)%Qvec(1),-2) !1/R^2

         CALL F1_sum_nd_PLUS_TO_Fres_sum_nd( RR_inv*LiLi ,PiPi)

         MPiPi = F_system%M_mass(2,2) * PiPi
         CALL F1_sum_nd_PLUS_TO_Fres_sum_nd(MPiPi,F_system%KEO)

       end if
     end if

     !!L1dag_L1
     do i = 2, nvec
     do j = 2, nvec
       if(i/=j) then
           CALL V1_scalar_V2_in_F_sum_nd(L_dag(i),L(j),LiLi)
           CALL F1_sum_nd_PLUS_TO_Fres_sum_nd(LiLi,L1L1)
       end if
     end do
     end do

     !!Lidag_J and Jdag_Li
     do i = 2, nvec

       call V1_scalar_V2_in_F_sum_nd(F_system%Jdag,L(i),JLi)
       CALL F1_sum_nd_MINUS_TO_Fres_sum_nd(JLi,L1L1)

       call V1_scalar_V2_in_F_sum_nd(L_dag(i),F_system%J,LiJ)
       CALL F1_sum_nd_MINUS_TO_Fres_sum_nd(LiJ,L1L1)

     end do

     !!L1dag \times Jsubsyst and J_dag_subsyst \times J_subsyst
     if(nsub_syst > 0) then
       call V1_scalar_V2_in_F_sum_nd(L_dag(1),V_sum_J, LiJ)
       CALL F1_sum_nd_MINUS_TO_Fres_sum_nd(LiJ,L1L1)

       call V1_scalar_V2_in_F_sum_nd(V_sum_Jdag,L(1),JLi)
       CALL F1_sum_nd_MINUS_TO_Fres_sum_nd(JLi,L1L1)

       do i = 1, F_system%nb_vect
       do j = 1, F_system%nb_vect
           if (F_system%tab_BFTransfo(i)%frame .and. F_system%tab_BFTransfo(j)%frame) then
             if(i==j) then
               if(compare_tab(F_system%tab_BFTransfo(i)%euler, (/.true., .true., .true./))) then

                 call Jdag_scalar_J_from_Eq122(falpha = F_system%tab_BFTransfo(i)%QEuler(1), &
                 &                             fbeta  = F_system%tab_BFTransfo(i)%QEuler(2), &
                 &                             fgamma = F_system%tab_BFTransfo(i)%QEuler(3), &
                 &                             JJ = LiLi)

               else if(compare_tab(F_system%tab_BFTransfo(i)%euler, (/.true., .true., .false./))) then
!
                 call Li_scalar_Li_from_Eq75(theta = F_system%tab_BFTransfo(i)%QEuler(2), &
                 &                           phi   = F_system%tab_BFTransfo(i)%QEuler(1), &
                 &                           LiLi  = LiLi)

               else if(compare_tab(F_system%tab_BFTransfo(i)%euler, (/.false., .true., .true./))) then
!
                 Ja_sum_subsyst = F_system%tab_BFTransfo(i)%J%vec_sum(3)

                 call Jdag_scalar_J_from_Eq171(F1_sum = Ja_sum_subsyst, &
                 &                             fbeta  = F_system%tab_BFTransfo(i)%QEuler(2), &
                 &                             fgamma = F_system%tab_BFTransfo(i)%QEuler(3), &
                 &                             JJ = LiLi)

                 call delete_op(Ja_sum_subsyst)

               else if(compare_tab(F_system%tab_BFTransfo(i)%euler, (/.false., .true., .false./))) then

                 call  V1_scalar_V2_in_F_sum_nd(F_system%tab_BFTransfo(i)%Jdag,&
                                                F_system%tab_BFTransfo(i)%J,LiLi)
               end if
             else
               call V1_scalar_V2_in_F_sum_nd(F_system%tab_BFTransfo(i)%Jdag, &
               &                             F_system%tab_BFTransfo(j)%J, LiLi)
             end if
             CALL F1_sum_nd_PLUS_TO_Fres_sum_nd(LiLi,L1L1)

           end if
       end do
       end do

       CALL V1_MINUS_TO_Vres(V_sum_J,   L(1))
       CALL V1_MINUS_TO_Vres(V_sum_Jdag,L_dag(1))

     end if

     if(true_BF) then
       call V1_scalar_V2_in_F_sum_nd(F_system%Jdag, F_system%J,JJ)
     else
       call Jdag_scalar_J_from_Eq122(falpha = Ja_el, &
       &                             fbeta  = jb_el, &
       &                             fgamma = jg_el, &
       &                             JJ = JJ)
     end if
     CALL F1_sum_nd_PLUS_TO_Fres_sum_nd(JJ,L1L1)

      iv = 1
      PiPi = get_Pq_dag (F_system%Qvec(1)) * get_Pq(F_system%Qvec(1)) ! PR^dag * PR

      RR_inv = get_Q(F_system%Qvec(1),-2) !1/R^2

     CALL F1_sum_nd_PLUS_TO_Fres_sum_nd( RR_inv*L1L1 ,PiPi)

     MPiPi = F_system%M_mass(iv,iv) * PiPi

     CALL F1_sum_nd_PLUS_TO_Fres_sum_nd(MPiPi,F_system%KEO)

     !get total J in the ref. frame
     if(.not.true_BF) then
       call get_opJ_projected_into_ref_frame(F_system%J,    Ja_el, Jb_el, Jg_el)
       call get_opJ_projected_into_ref_frame(F_system%Jdag, Ja_el, Jb_el, Jg_el, dag = .true.)
      end if

     !!! off diagonal terms
     do i = 1, nvec
     do j = 1, nvec
         if(i/=j  .and. abs(F_system%M_mass(i,j)%Cn(1)) >1.0e-13_Rkind) then
           call V1_scalar_V2_in_F_sum_nd(Pi_dag_BF(F_system%listVFr(i)),&
           &                             Pi_BF(F_system%listVFr(j)),PiPi)

           MPiPi = F_system%M_mass(i,j) * PiPi
           CALL F1_sum_nd_PLUS_TO_Fres_sum_nd(MPiPi,F_system%KEO)

         end if
     end do
     end do

     call delete_op(V1_tmp)
     call delete_op(V2_tmp)
     call delete_op(V_sum_J)
     call delete_op(V_sum_Jdag)
     call delete_op(keo_Pi)
     call delete_op(Ja_sum_subsyst)

     call delete_op(RR_inv)
     call delete_op(MPRPR)
     call delete_op(PRPR)
     call delete_op(PiPi)
     call delete_op(MPiPi)
     call delete_op(LiLi)
     call delete_op(JLi)
     call delete_op(LiJ)
     call delete_op(JJ)
     call delete_op(L1L1)

     IF (debug) THEN
       CALL write_op(F_system%keo)
       write(out_unitp,*) 'END ',routine_name
     END IF

   END SUBROUTINE get_opKEO_subsyst

   !! @param:       KOE          The output
   !! @param:       nvec    Integer corresponding to the size of the system
   !! @param:       M_mass  The mass matrix
   RECURSIVE SUBROUTINE get_opKEO_subsyst_nvectot1(F_system, P_Euler, M_mass_out, scalar_PiPj, &
                                                   F_system_parent)
     USE mod_BunchPolyTransfo, only : Type_BFTransfo

     type(Type_BFtransfo),                 intent(inout)      :: F_system
     type(Type_PiEulerRot),                intent(inout)      :: P_Euler(:)
     type(sum_opnd),                       intent(in)         :: M_mass_out(:,:)
     logical,                              intent(inout)      :: scalar_PiPj(:,:)
     type(Type_BFtransfo),                 intent(in)         :: F_system_parent

     type(sum_opnd)                          :: RR_inv,PiPi,LiLi
     type(sum_opnd)                          :: Ja_sum, Ja_sum_dag
     type(vec_sum_opnd), allocatable         :: E(:)
     type(vec_sum_opnd), allocatable         :: L(:)
     type(vec_sum_opnd), allocatable         :: L_dag(:)
     type(sum_opnd), allocatable             :: Liz_parent(:)
     integer                                 :: i,j,n_size,nvec,nvec_parent

     logical, parameter           :: debug=.TRUE.
     !logical, parameter           :: debug=.FALSE.
     character (len = *), parameter :: routine_name= 'get_opKEO_subsyst_nvectot1'

     IF (debug) THEN
       write(out_unitp,*) 'BEGINNING ',routine_name
     END IF

     IF (F_system%nb_vect_tot /= 1) THEN
       write(out_unitp,*) ' ERROR in ',routine_name
       write(out_unitp,*) ' F_system%nb_vect_tot > 1 is not possible in this subroutine'
       write(out_unitp,*) ' F_system%nb_vect_tot: ',F_system%nb_vect_tot
       STOP
     END IF


     nvec = 1

     n_size = size(F_system%tab_num_Frame)
     do i = 1, nvec
       if(.not.associated(P_Euler(F_system%listVFr(i))%Tab_num_Frame)) then

         CALL alloc_array(P_Euler(F_system%listVFr(i))%Tab_num_Frame,(/n_size/), &
                         'P_Euler(F_system%listVFr(i))%Tab_num_Frame',routine_name)
         P_Euler(F_system%listVFr(i))%Tab_num_Frame(:) = F_system%Tab_num_Frame(:)

       end if
     end do

       write(out_unitp,*) 'processing S_(',F_system%tab_num_frame,')'
       CALL flush_perso(out_unitp)
       CALL alloc_NParray(E,    (/nvec/),'E',routine_name)
       CALL alloc_NParray(L,    (/nvec/),'L',routine_name)
       CALL alloc_NParray(L_dag,(/nvec/),'L_dag',routine_name)

       if( compare_tab(F_system%euler, (/.true., .true., .false./))) then

         CALL Li_scalar_Li_from_Eq75(theta = F_system%QEuler(2), phi = F_system%QEuler(1), LiLi = LiLi)

         PiPi = get_Pq_dag(F_system%Qvec(1)) * get_Pq(F_system%Qvec(1)) ! PR^dag * PR

         RR_inv = get_Q(F_system%Qvec(1),-2) !1/R^2

         CALL F1_sum_nd_PLUS_TO_Fres_sum_nd( RR_inv*LiLi ,PiPi)

         F_system%KEO = F_system%M_mass(1,1) * PiPi

         CALL delete_op(LiLi)
         CALL delete_op(PiPi)
         CALL delete_op(RR_inv)

         !projection onto the ref. frame
         call get_opLi(L = F_system%J,    theta = F_system%QEuler(2), phi = F_system%QEuler(1), index_L = 3)
         call get_opLi(L = F_system%Jdag, theta = F_system%QEuler(2), phi = F_system%QEuler(1), index_L = 3, dag = .true.)
         call copy_F1_into_F2(F_system%J,    L(1))
         call copy_F1_into_F2(F_system%Jdag, L_dag(1))

         call get_unit_vector_Ei(E(1), F_system%QEuler(2), F_system%QEuler(1), 3) ! the F_system%Unit_vector cannot be used

         call get_opPi(       P_Euler(F_system%listVFr(1))%Pi,    F_system%Qvec(1),&
                              L(1),E(1))
         call get_opPi_dagger(P_Euler(F_system%listVFr(1))%Pidag, F_system%Qvec(1),&
                              L_dag(1), E(1))

       else if(compare_tab(F_system%euler, (/.false., .false., .false./))) then

         STOP 'FFF'

       else if(compare_tab(F_system%euler, (/.false., .true., .false./))) then

         nvec_parent = count( .NOT. F_system_parent%tab_BFTransfo(:)%frame )
         IF (nvec_parent >=1 ) THEN
           CALL alloc_NParray(Liz_parent,(/nvec_parent/),'Liz_parent',routine_name)
           CALL get_Lz_F_system_parent(F_system_parent, Liz_parent)
         END IF

         if(compare_tab(F_system_parent%euler, (/.false., .false., .false./))) then ! true BF

           Ja_sum     = get_Jz(F_system_parent%QEuler(3))
           Ja_sum_dag = get_Jz(F_system_parent%QEuler(3))

         else if(compare_tab(F_system_parent%euler, (/.true., .true., .true./))) then

           Ja_sum     = get_Pq(F_system_parent%QEuler(3))
           Ja_sum_dag = get_Pq_dag(F_system_parent%QEuler(3))

         else if(compare_tab(F_system_parent%euler, (/.false., .true., .true./))) then

           Ja_sum     = get_Pq(F_system_parent%QEuler(3))
           Ja_sum_dag = get_Pq_dag(F_system_parent%QEuler(3))

         end if

         DO i = F_system_parent%nb_vect, 2, -1
           IF (F_system_parent%tab_BFTransfo(i)%frame) THEN
             CALL F1_sum_nd_MINUS_TO_Fres_sum_nd(F_system_parent%tab_BFTransfo(i)%J%vec_sum(3) ,Ja_sum)
             CALL F1_sum_nd_MINUS_TO_Fres_sum_nd(F_system_parent%tab_BFTransfo(i)%Jdag%vec_sum(3) ,Ja_sum_dag)
           END IF
         END DO

         j = 0
         DO i = 1, F_system_parent%nb_vect
           IF (.not.F_system_parent%tab_BFTransfo(i)%frame) THEN
             j =j+1
             CALL F1_sum_nd_MINUS_TO_Fres_sum_nd(Liz_parent(j) ,Ja_sum)
             CALL F1_sum_nd_MINUS_TO_Fres_sum_nd(Liz_parent(j) ,Ja_sum_dag)
           END IF
         END DO

         call get_opL1_beta(L1 = L(1),     beta = F_system%QEuler(2), J_a = Ja_sum)
         call get_opL1_beta(L1 = L_dag(1), beta = F_system%QEuler(2), J_a = Ja_sum_dag, dag=.true.)

         !projection onto the ref. frame
         call copy_F1_into_F2(L(1),     F_system%J)
         call copy_F1_into_F2(L_dag(1), F_system%Jdag)

         call get_unit_vector_Ei(E(1), F_system%QEuler(2), F_system%QEuler(1), 2) ! the second angle is not used


         call get_opPi(       P_Euler(F_system%listVFr(1))%Pi,    F_system%Qvec(1), L(1),     E(1))
         call get_opPi_dagger(P_Euler(F_system%listVFr(1))%Pidag, F_system%Qvec(1), L_dag(1), E(1))


         call  V1_scalar_V2_in_F_sum_nd(P_Euler(F_system%listVFr(1))%Pidag,   &
                                        P_Euler(F_system%listVFr(1))%Pi, PiPi)

         F_system%keo = F_system%M_mass(1,1) * PiPi

         CALL delete_op(PiPi)
         CALL delete_op(Ja_sum)
         CALL delete_op(Ja_sum_dag)

       end if

       CALL dealloc_NParray(L,'L',routine_name)
       CALL dealloc_NParray(L_dag,'L_dag',routine_name)
       CALL dealloc_NParray(E,'E',routine_name)
       CALL dealloc_NParray(Liz_parent,'Liz_parent',routine_name)

     IF (debug) THEN
       CALL write_op(F_system%keo)
       write(out_unitp,*) 'END ',routine_name
     END IF

   END SUBROUTINE get_opKEO_subsyst_nvectot1

   !! @param:       KOE          The output 
   !! @param:       nvec    Integer corresponding to the size of the system
   !! @param:       M_mass  The mass matrix
   RECURSIVE SUBROUTINE get_opKEO_subsyst_2euler(F_system, &
                                         P_Euler, M_mass_out, scalar_PiPj, &
                                         F_system_parent)
     USE mod_BunchPolyTransfo, only : Type_BFTransfo

     type(Type_BFtransfo),         intent(inout)      :: F_system
     type(Type_PiEulerRot),        intent(inout)      :: P_Euler(:) 
     type(sum_opnd),               intent(in)         :: M_mass_out(:,:) 
     logical,                      intent(inout)      :: scalar_PiPj(:,:) 
     type(Type_BFtransfo),         intent(in)         :: F_system_parent

     type(vec_sum_opnd), allocatable :: L(:)
     type(vec_sum_opnd), allocatable :: L_dag(:)
     type(sum_opnd),     allocatable :: Lz(:)
     type(sum_opnd),     allocatable :: Liz_parent(:)
     type(vec_sum_opnd), allocatable :: Pi_BF(:)
     type(vec_sum_opnd), allocatable :: Pi_dag_BF(:)
     type(sum_opnd),     allocatable :: Mat_R(:,:)
     type(sum_opnd),     allocatable :: Mat_RTranspo(:,:)
     type(sum_opnd)                  :: RR_inv,MPRPR,PRPR,PiPi,MPiPi
     type(sum_opnd)                  :: LiLi,JLi,LiJ,JJ,L1L1

     type(vec_sum_opnd)                  :: V1_tmp 
     type(vec_sum_opnd)                  :: V2_tmp 
     type(vec_sum_opnd)                  :: V_sum_Li 
     type(vec_sum_opnd)                  :: V_sum_Lidag 
     type(vec_sum_opnd)                  :: V_sum_J 
     type(vec_sum_opnd)                  :: V_sum_Jdag 

     type(sum_opnd)                      :: Ja_sum ,Ja_sum_dag
     type(sum_opnd)                      :: Ja_sum_subsyst 
     type(sum_opnd)                      :: Jdag_J_sub 


     integer                         :: i, j, n, k, iv, index_L
     integer                         :: k1, k2, K3
     integer                         :: nvec
     integer                         :: nsub_syst, nvec_tot
     integer                         :: nvec_parent
     integer                         :: error
     logical                         :: parent_true_BF
     logical                         :: alloc_subsyst
     logical, allocatable            :: zero_Pi_BF(:)

     !logical, parameter           :: debug=.TRUE.
     logical, parameter           :: debug=.FALSE.
     character (len = *), parameter  :: routine_name='get_opKEO_subsyst_2euler'

     IF (debug) THEN
       write(out_unitp,*) 'BEGINNING ',routine_name
     END IF

     nsub_syst = 0
     do i=1, F_system%nb_vect
       if(F_system%tab_BFTransfo(i)%frame) nsub_syst = nsub_syst+1
     end do

     nvec = F_system%nb_vect-nsub_syst+1
     nvec_tot = F_system%nb_vect_tot

     do i = 1, nvec
       if(.not. associated(P_Euler(F_system%listVFr(i))%Tab_num_Frame)) then
          CALL alloc_array(P_Euler(F_system%listVFr(i))%Tab_num_Frame,shape(F_system%tab_num_Frame), &
                          'P_Euler(F_system%listVFr(i))%Tab_num_Frame',routine_name)
          P_Euler(F_system%listVFr(i))%Tab_num_Frame(:) = F_system%Tab_num_Frame(:)
       end if
     end do

     if(.not.compare_tab(F_system%euler, (/.false., .true., .true./))) then
       write(out_unitp,*) ' ERROR in',routine_name
       write(out_unitp,*) "This routine can be call only for a subsystem"
       write(out_unitp,*) "which has two Euler's angles"
       STOP
     end if

     parent_true_BF = .true.
     do i = 1, size(F_system_parent%euler)
       if(F_system_parent%euler(i)) then
         parent_true_BF = .false. 
         exit
       end if
     end do


     CALL alloc_NParray(L,(/nvec/),'L',routine_name)
     CALL alloc_NParray(L_dag,(/nvec/),'L_dag',routine_name)
     if(nsub_syst > 0) then
       CALL alloc_NParray(Lz,(/nvec+nsub_syst/),'Lz',routine_name)
     else
       CALL alloc_NParray(Lz,(/nvec/),'Lz',routine_name)
     end if
     CALL alloc_NParray(Pi_BF,     shape(P_Euler),'Pi_BF',routine_name)
     CALL alloc_NParray(Pi_dag_BF, shape(P_Euler),'Pi_dag_BF',routine_name)
     CALL alloc_NParray(zero_Pi_BF,shape(P_Euler),'zero_Pi_BF',routine_name)

     !Computation the elementary operators
     write(out_unitp,*) 'get elementaries op. S_(',F_system%tab_num_frame,')' 
     write(out_unitp,*) 'with euler(1) = false' 
     CALL flush_perso(out_unitp)
     iv = 1

     nvec_parent = count( .NOT. F_system_parent%tab_BFTransfo(:)%frame )
     IF (nvec_parent >= 1 ) THEN
       CALL alloc_NParray(Liz_parent,(/nvec_parent/),'Liz_parent',routine_name)
       CALL get_Lz_F_system_parent(F_system_parent, Liz_parent)
     END IF


     if (parent_true_BF) then

       Ja_sum =     get_Jz(F_system_parent%QEuler(3))
       Ja_sum_dag = get_Jz(F_system_parent%QEuler(3))

     else if(compare_tab(F_system_parent%euler, (/.true., .true., .true./)) .OR. &
            compare_tab(F_system_parent%euler, (/.false., .true., .true./))) then

        Ja_sum     = get_Pq(    F_system_parent%QEuler(3))
        Ja_sum_dag = get_Pq_dag(F_system_parent%QEuler(3))


     else if(compare_tab(F_system_parent%euler, (/.true., .true., .false./))) then
       write(out_unitp,*) 'Case not possible'
       stop
     end if
    ! write(out_unitp,*) 'euler parent', F_system_parent%euler 
    ! write(out_unitp,*) 'euler parent', F_system_parent%tab_num_frame 

     do i =  F_system_parent%nb_vect, 2, -1
       if(F_system_parent%tab_BFTransfo(i)%frame) then
         CALL F1_sum_nd_MINUS_TO_Fres_sum_nd( F_system_parent%tab_BFTransfo(i)%J%vec_sum(3) ,   Ja_sum)
         CALL F1_sum_nd_MINUS_TO_Fres_sum_nd( F_system_parent%tab_BFTransfo(i)%Jdag%vec_sum(3) ,Ja_sum_dag)

       end if
     end do

     j = 0
     do i = F_system_parent%nb_vect, 2, -1
       if(.not.F_system_parent%tab_BFTransfo(i)%frame) then
         j = j+1
         CALL F1_sum_nd_MINUS_TO_Fres_sum_nd( Liz_parent(j) ,Ja_sum)
         CALL F1_sum_nd_MINUS_TO_Fres_sum_nd( Liz_parent(j) ,Ja_sum_dag)
       end if
     end do
     CALL dealloc_NParray(Liz_parent,'Liz_parent',routine_name)

     !write(out_unitp,*) 'Ja_sum'
     !call write_op(Ja_sum, out_unitp, header=.true.)
     !write(out_unitp,*) 'END Ja_sum'

     call get_opJ_projected_into_BFEq171(F_system%J,    Ja_sum,     F_system%QEuler(2),F_system%QEuler(3))
     call get_opJ_projected_into_BFEq171(F_system%Jdag, Ja_sum_dag, F_system%QEuler(2),F_system%QEuler(3), dag = .true.)

     do i = 1, size(Pi_BF)
       call allocate_op(Pi_BF(i),3)
       call allocate_op(Pi_dag_BF(i),3)
       zero_Pi_BF(i) = .true.
       do j = 1, 3
         Pi_BF(i)%vec_sum(j)     = czero
         Pi_dag_BF(i)%vec_sum(j) = czero
       end do
     end do
     call get_Euler_MatRot(F_system, Mat_R, Mat_RTranspo)

     do i = 1, size(Lz)
       Lz(i) = czero
     end do
     if(nvec>=3) then
       write(out_unitp,*) 'computation of Pi, i>=3 for S_(',F_system%tab_num_frame,')' 
       write(out_unitp,*) 'with euler(1) = false' 
     end if

     ! Initiliation of the Li (i>2) operators
     !computation of the Pi, i>2 in their local frame
     call allocate_op(V1_tmp, 3)
     CALL zero_TO_vec_sum_opnd(V1_tmp)
     call allocate_op(V2_tmp, 3)
     CALL zero_TO_vec_sum_opnd(V2_tmp)

     iv = 1
     do i = 1, F_system%nb_vect
       if(.not.F_system%tab_BFTransfo(i)%frame) then
         iv=iv+1
         if(i>1) then
           if(F_system%tab_BFTransfo(i)%cart) then

             V1_tmp%vec_sum(1) = get_Q(F_system%tab_BFTransfo(i)%Qvec(1)) !x
             V1_tmp%vec_sum(2) = get_Q(F_system%tab_BFTransfo(i)%Qvec(2)) !y
             V1_tmp%vec_sum(3) = get_Q(F_system%tab_BFTransfo(i)%Qvec(3)) !z

             Pi_BF(F_system%listVFr(iv))%vec_sum(1) = get_Pq(F_system%tab_BFTransfo(i)%Qvec(1)) !Px
             Pi_BF(F_system%listVFr(iv))%vec_sum(2) = get_Pq(F_system%tab_BFTransfo(i)%Qvec(2)) !Py
             Pi_BF(F_system%listVFr(iv))%vec_sum(3) = get_Pq(F_system%tab_BFTransfo(i)%Qvec(3)) !Pz

             call V1_cross_V2_in_Vres(V1_tmp, Pi_BF(F_system%listVFr(iv)), L(iv))
             call V1_cross_V2_in_Vres(Pi_BF(F_system%listVFr(iv)), V1_tmp, L_dag(iv))

             ! L_dag(iv) = -L_dag(iv)
             call V1_plus_V2_in_Vres(V2_tmp, L_dag(iv), V1_tmp, minus=.true.)
             call copy_F1_into_F2(V1_tmp,  L_dag(iv))

             call copy_F1_into_F2(Pi_BF(F_system%listVFr(iv)), Pi_dag_BF(F_system%listVFr(iv)))

           else
             IF (iv == 2) THEN
               index_L = iv+1
             ELSE
               index_L = iv
             END IF

             call get_opLi(L = L(iv),     theta = F_system%tab_BFTransfo(i)%Qvec(2), &
                                          phi   = F_system%tab_BFTransfo(i)%Qvec(3), &
                                  index_L = index_L, Li=F_system%tab_BFTransfo(i)%Li)

             call get_opLi(L = L_dag(iv), theta = F_system%tab_BFTransfo(i)%Qvec(2), &
                                          phi   = F_system%tab_BFTransfo(i)%Qvec(3), &
                    index_L = index_L, dag = .true., Li=F_system%tab_BFTransfo(i)%Li)

             call get_opPi(       Pi_BF(F_system%listVFr(iv)),                       &
                           FRel=F_system%tab_BFTransfo(i)%Qvec(1), L=L(iv),          &
                           E=F_system%tab_BFTransfo(i)%Unit_Vector)
             call get_opPi_dagger(Pi_dag_BF(F_system%listVFr(iv)),                   &
                           FRel=F_system%tab_BFTransfo(i)%Qvec(1), L_dag=L_dag(iv),  &
                           E=F_system%tab_BFTransfo(i)%Unit_Vector)

           end if
           Lz(iv) = L(iv)%vec_sum(3)
           zero_Pi_BF(F_system%listVFr(iv)) = .false.
         end if
       end if
     end do

     if(nsub_syst > 0) then
       j = 0
       do i = 1, F_system%nb_vect
         if(F_system%tab_BFTransfo(i)%frame) then
           j = j+1
           Lz(j+nvec) = F_system%tab_BFTransfo(i)%J%vec_sum(3)
         end if
       end do
     end if

     if(F_system%nb_vect >0) then
       if(.not.F_system%tab_BFTransfo(1)%frame) then
         write(out_unitp,*) 'computation of P2, for S_(',F_system%tab_num_frame,')' 
         write(out_unitp,*) 'with euler(1) = false' 

         call get_opL2(L = L(2),     Fel = F_system%tab_BFTransfo(1)%Qvec(2), &
                       Jz = F_system%J%vec_sum(3), Lz_all = Lz, index_L = 2)
         call get_opL2(L = L_dag(2), Fel = F_system%tab_BFTransfo(1)%Qvec(2), &
                       Jz = F_system%Jdag%vec_sum(3), Lz_all = Lz, index_L = 2, dag = .true.)

         call get_opPi(       Pi_BF(F_system%listVFr(2)),     F_system%tab_BFTransfo(1)%Qvec(1), L(2),    &
                       E=F_system%tab_BFTransfo(1)%Unit_Vector)
         call get_opPi_dagger(Pi_dag_BF(F_system%listVFr(2)), F_system%tab_BFTransfo(1)%Qvec(1), L_dag(2),&
                       E=F_system%tab_BFTransfo(1)%Unit_Vector)

         zero_Pi_Bf(F_system%listVFr(2)) = .false.
       end if
     end if

     if(nvec == 1) then
       call copy_F1_into_F2(F_system%J,    L(1))
       call copy_F1_into_F2(F_system%Jdag, L_dag(1))
     else
       call get_opL1(L1 = L(1),     J = F_system%J,    L_all = L,     index_L = 1)
       call get_opL1(L1 = L_dag(1), J = F_system%Jdag, L_all = L_dag, index_L = 1)
     end if


     call allocate_op(V_sum_J, 3)
     CALL zero_TO_vec_sum_opnd(V_sum_J)
     call allocate_op(V_sum_Jdag, 3)
     CALL zero_TO_vec_sum_opnd(V_sum_Jdag)

     DO i = 1, F_system%nb_vect
       IF (F_system%tab_BFTransfo(i)%frame) THEN

         CALL V1_PLUS_TO_Vres( F_system%tab_BFTransfo(i)%J,    V_sum_J )
         CALL V1_PLUS_TO_Vres( F_system%tab_BFTransfo(i)%Jdag, V_sum_Jdag)

       END IF
     END DO


     !!! Diagonal terms
     !Compute sum1
     L1L1  = CZERO
     F_system%KEO = CZERO

     iv = 1

     do i = 1, F_system%nb_vect
       if(.not.F_system%tab_BFTransfo(i)%frame) then
         iv=iv+1
         if(i>1) then
           if(F_system%tab_BFTransfo(i)%cart) then

             CALL V1_scalar_V2_in_F_sum_nd(L_dag(iv), L(iv), LiLi)

             CALL V1_scalar_V2_in_F_sum_nd(Pi_dag_BF(F_system%listVFr(iv)), &
                                               Pi_BF(F_system%listVFr(iv)), PiPi)

           else

             CALL Li_scalar_Li_from_Eq75(theta = F_system%tab_BFTransfo(i)%Qvec(2), &
             &                           phi   = F_system%tab_BFTransfo(i)%Qvec(3), &
             &                           LiLi  = LiLi)

             PiPi = get_Pq_dag (F_system%tab_BFTransfo(i)%Qvec(1)) *       &
                    get_Pq(F_system%tab_BFTransfo(i)%Qvec(1)) ! PR^dag * PR

             RR_inv = get_Q(F_system%tab_BFTransfo(i)%Qvec(1),-2) !1/R^2

             CALL F1_sum_nd_PLUS_TO_Fres_sum_nd( RR_inv*LiLi ,PiPi)

           end if

           CALL F1_sum_nd_PLUS_TO_Fres_sum_nd(LiLi,L1L1)

           MPiPi = F_system%M_mass(iv,iv) * PiPi

           CALL F1_sum_nd_PLUS_TO_Fres_sum_nd(MPiPi,F_system%KEO)

         end if
       end if
     end do

     !! L2dag_L2 and M.Pr^dag* Pr
     iv=2
     if(F_system%nb_vect>0) then
       if( .not.F_system%tab_BFTransfo(1)%frame) then

         CALL V1_scalar_V2_in_F_sum_nd(L_dag(iv), L(iv), LiLi)
         CALL F1_sum_nd_PLUS_TO_Fres_sum_nd(LiLi, L1L1)

         PiPi = get_Pq_dag (F_system%tab_BFTransfo(1)%Qvec(1)) * &
                get_Pq(F_system%tab_BFTransfo(1)%Qvec(1))  ! PR^dag * PR

         RR_inv = get_Q(F_system%tab_BFTransfo(1)%Qvec(1),-2) !1/R^2

         CALL F1_sum_nd_PLUS_TO_Fres_sum_nd( RR_inv*LiLi ,PiPi)

         MPiPi = F_system%M_mass(iv,iv) * PiPi

         CALL F1_sum_nd_PLUS_TO_Fres_sum_nd(MPiPi,F_system%KEO)

       end if
     end if

     !!L1dag_L1
     do i = 2, nvec
     do j = 2, nvec
         if(i/=j) then
           CALL V1_scalar_V2_in_F_sum_nd(L_dag(i), L(j), LiLi)
           CALL F1_sum_nd_PLUS_TO_Fres_sum_nd(LiLi, L1L1)
         end if
     end do
     end do

     do i = 2, nvec

       CALL V1_scalar_V2_in_F_sum_nd(F_system%Jdag, L(i), JLi)
       CALL F1_sum_nd_MINUS_TO_Fres_sum_nd(JLi, L1L1)

       CALL V1_scalar_V2_in_F_sum_nd(L_dag(i), F_system%J, LiJ)
       CALL F1_sum_nd_MINUS_TO_Fres_sum_nd(LiJ, L1L1)

     end do

     if(nsub_syst > 0) then
       CALL V1_scalar_V2_in_F_sum_nd(L_dag(1), V_sum_J, LiJ)
       CALL F1_sum_nd_MINUS_TO_Fres_sum_nd(LiJ, L1L1)

       call V1_scalar_V2_in_F_sum_nd(V_sum_Jdag, L(1), JLi)
       CALL F1_sum_nd_MINUS_TO_Fres_sum_nd(JLi, L1L1)

       call Jdag_scalarJ_subsystem(F_system, Jdag_J_sub)
       CALL F1_sum_nd_PLUS_TO_Fres_sum_nd(Jdag_J_sub, L1L1)

       CALL V1_MINUS_TO_Vres( V_sum_J,   L(1) )
       CALL V1_MINUS_TO_Vres( V_sum_Jdag,L_dag(1) )

     end if

     !write(out_unitp,*) 'L1L1'
     !CALL write_op(L1L1)
     !write(out_unitp,*) 'END L1L1'
     !STOP

     write(out_unitp,*) 'computation of P1, for S_(',F_system%tab_num_frame,')'
     write(out_unitp,*) 'with euler(1) = false'

     call get_opPi(       Pi_BF(F_system%listVFr(1)),     F_system%Qvec(1), L(1),    &
                          E=F_system%Unit_Vector)
     call get_opPi_dagger(Pi_dag_BF(F_system%listVFr(1)), F_system%Qvec(1), L_dag(1),&
                          E=F_system%Unit_Vector)

     zero_Pi_Bf(F_system%listVFr(1)) = .false.

     do i = 1, F_system%nb_vect  
       if(F_system%tab_BFTransfo(i)%frame) then
         call get_Pi_subsyst(F_system%tab_BFTransfo(i), P_Euler, Pi_BF, &
         &                   Pi_dag_BF, zero_Pi_BF)
       end if
     end do

! for coupled vectors, if they does not been taking into account (scalar_PiPj(i,j))

       do i = 1, size(Pi_BF)
       do j = 1, size(Pi_BF)

         if (i/=j .and. .not.zero_Pi_BF(i) .and. .not.zero_Pi_BF(j)) then
         if (.not.compare_tab(P_euler(i)%tab_num_frame,P_euler(j)%tab_num_frame) .and. &
                   .not.scalar_PiPj(i,j)) then

           call  V1_scalar_V2_in_F_sum_nd(Pi_dag_BF(i), Pi_BF(j), PiPi)

     !MPiPi = F_system%M_mass(F_system%listVFr(i), F_system%listVFr(j)) * PiPi ! ne marche pas
     !MPiPi = F_system%M_mass(i,j) * PiPi ! ne marche pas

           MPiPi = M_mass_out(i,j) * PiPi

           CALL F1_sum_nd_PLUS_TO_Fres_sum_nd(MPiPi,F_system%KEO)

           scalar_PiPj(i,j) = .true.
         end if
         end if
       end do
       end do

       do i = 1, nvec
         call M_opnd_times_V_in_Vres(Mat_R, Pi_BF(F_system%listVFr(i)), &
         &                           P_Euler(F_system%listVFr(i))%Pi)

         call V_times_M_opnd_in_Vres(Pi_dag_BF(F_system%listVFr(i)), &
       &                             Mat_RTranspo, P_Euler(F_system%listVFr(i))%Pidag)
       end do
    
       do i = 1, F_system%nb_vect  
         if(F_system%tab_BFTransfo(i)%frame) then
           call Mat_Rot_times_Pi_subsyst(F_system%tab_BFTransfo(i), &
           &                   Mat_R, Mat_RTranspo, P_Euler)
         end if
       end do

      do i = 1, nvec
      do j = 1, nvec
          scalar_PiPj(F_system%listVFr(i), F_system%listVFr(j)) = .true.
      end do
      end do

      CALL dealloc_NParray(zero_Pi_BF,'zero_Pi_BF',routine_name)
      CALL dealloc_NParray(Mat_R,'Mat_R',routine_name)
      CALL dealloc_NParray(Mat_RTranspo,'Mat_RTranspo',routine_name)
! END for coupled vectors, if they does not been taking into account (scalar_PiPj(i,j))


     !! end J1dagJ1, P1dag.P1
     call Jdag_scalar_J_from_Eq171(F1_sum = Ja_sum,             &
     &                             fbeta  = F_system%QEuler(2), &
     &                             fgamma = F_system%QEuler(3), &
     &                             JJ = JJ)

     CALL F1_sum_nd_PLUS_TO_Fres_sum_nd(JJ,L1L1)

     iv = 1
     PiPi = get_Pq_dag(F_system%Qvec(1)) * get_Pq(F_system%Qvec(1)) ! PR^dag * PR


     RR_inv = get_Q(F_system%Qvec(1),-2) !1/R^2

     CALL F1_sum_nd_PLUS_TO_Fres_sum_nd( RR_inv*L1L1 ,PiPi)

     MPiPi = F_system%M_mass(iv,iv) * PiPi

     CALL F1_sum_nd_PLUS_TO_Fres_sum_nd(MPiPi,F_system%KEO)


     !get total J in the ref. frame
     call get_opJ_projected_into_ref_frameEq170(F_system%J,    Ja_sum,     &
                                                F_system%QEuler(2),F_system%QEuler(3))
     call get_opJ_projected_into_ref_frameEq170(F_system%Jdag, Ja_sum_dag, &
                                                F_system%QEuler(2),F_system%QEuler(3), dag = .true.)
     call delete_op(Ja_sum)
     call delete_op(Ja_sum_dag)


     !!! off diagonal terms
     !Compute sum1 

     do i = 1, nvec
     do j = 1, nvec
       if(i/=j  .and. abs(F_system%M_mass(i,j)%Cn(1)) >1.0e-13_Rkind) then

         CALL V1_scalar_V2_in_F_sum_nd(Pi_dag_BF(F_system%listVFr(i)),  &
                                       Pi_BF(F_system%listVFr(j)), PiPi)
         MPiPi = F_system%M_mass(i,j) * PiPi
         CALL F1_sum_nd_PLUS_TO_Fres_sum_nd(MPiPi,F_system%KEO)

       end if
     end do
     end do

     CALL dealloc_NParray(Pi_BF,       'Pi_BF',routine_name)
     CALL dealloc_NParray(Pi_dag_BF,   'Pi_dag_BF',routine_name)
     CALL dealloc_NParray(L,           'L',routine_name)
     CALL dealloc_NParray(Lz,          'Lz',routine_name)
     CALL dealloc_NParray(L_dag,       'L_dag',routine_name)

     call delete_op(RR_inv)
     call delete_op(MPRPR)
     call delete_op(PRPR)
     call delete_op(PiPi)
     call delete_op(MPiPi)
     call delete_op(LiLi)
     call delete_op(JLi)
     call delete_op(LiJ)
     call delete_op(JJ)
     call delete_op(L1L1)

     call delete_op(V1_tmp)
     call delete_op(V2_tmp)
     call delete_op(V_sum_J)
     call delete_op(V_sum_Jdag)
     call delete_op(V_sum_Li)
     call delete_op(V_sum_Lidag)

     call delete_op(Ja_sum)
     call delete_op(Ja_sum_subsyst)
     call delete_op(Jdag_J_sub)

     IF (debug) THEN
       CALL write_op(F_system%keo)
       write(out_unitp,*) 'END ',routine_name
     END IF

   END SUBROUTINE get_opKEO_subsyst_2euler

 !! @description:  Extract the terms of the KEO for active coordinates
 !! @param.in:     F_sum_nd   Kinetic operator energy (type: sum_opnd). 
 !! @param.in:     list_QdyntoQact New list of the coordinate. The active
 !!                coordinates are classified before the non active coordinates
  subroutine get_KEO_for_Qactiv(TWOxKEO, constraint, Qval,tabQpoly_Qel,tabQact_Qel, &
                                list_Qactiv,list_QpolytoQact)

   type(sum_opnd),      intent(inout)             :: TWOxKEO
   logical,             intent(in)                :: constraint
   integer,             intent(in)                :: list_QpolytoQact(:)
   integer,             intent(in)                :: list_Qactiv(:)
   real(kind=Rkind),    intent(in)                :: Qval(:) ! Qact order ???
   TYPE(opel),          intent(in)                :: tabQpoly_Qel(:)
   TYPE(opel),          intent(inout)             :: tabQact_Qel(:)


   TYPE(opel), allocatable    :: tab_Qel_loc(:)
   integer                    :: i, k
   integer                    :: indexq
   complex(kind=Rkind)        :: opval

   integer :: Pq(2),JJ(2),LL(2),nb_var,nb_act

   !logical, parameter :: debug = .TRUE.
   logical, parameter :: debug = .FALSE.
   character (len=*), parameter :: routine_name='get_keo_for_Qactiv'

   if(.not.allocated(TWOxKEO%sum_prod_op1d)) then
     write(out_unitp,*) ' ERROR in',routine_name
     write(out_unitp,*) "TWOxKEO should be allocated"
     STOP
   end if

   nb_var = size(list_QpolytoQact)
   nb_act = count(list_Qactiv == 1 .OR. list_Qactiv == 21)

   IF (debug) THEN
     write(out_unitp,*) ' BEGINNING in ',routine_name
     write(out_unitp,*) 'nb_var',nb_var
     write(out_unitp,*) 'nb_act',nb_act
     write(out_unitp,*) 'Qval',Qval
     write(out_unitp,*) 'list_Qactiv',list_Qactiv
     write(out_unitp,*) 'list_QpolytoQact',list_QpolytoQact
     write(out_unitp,*) 'before removing inactive terms'
     CALL write_op(TWOxKEO,header=.TRUE.)
   END IF

   DO i=1,size(tabQpoly_Qel)
     IF (i <= nb_var) THEN
       tabQact_Qel(list_QpolytoQact(i)) = tabQpoly_Qel(i)
       tabQact_Qel(list_QpolytoQact(i))%indexQ = list_QpolytoQact(i) ! it has to be changed
     ELSE
       tabQact_Qel(i) = tabQpoly_Qel(i)
     END IF
   END DO

   DO i = 1, size(TWOxKEO%sum_prod_op1d)
   DO k = 1, size(TWOxKEO%sum_prod_op1d(i)%prod_op1d)
     indexq = get_indexQ_OF_Op1D(TWOxKEO%sum_prod_op1d(i)%prod_op1d(k))
     IF (indexq <= nb_var) then ! we don't change the index for Jx,Jy,Jz
       TWOxKEO%sum_prod_op1d(i)%prod_op1d(k)%prod_opel(:)%indexq =      &
                                               list_QpolytoQact(indexq)
     END IF
   END DO
   END DO


   IF (debug) THEN
     write(out_unitp,*) ' Change indexq (active order)'
     CALL write_op(TWOxKEO,header=.TRUE.)
   END IF

   IF (nb_act == nb_var) THEN
     IF (debug) THEN
       write(out_unitp,*) ' no constraint'
       write(out_unitp,*) ' END in',routine_name
     END IF
     RETURN
   END IF


   ! remove all terms with PQi PQj where Qi or Qj are inactive coordinates
   DO i = 1, size(TWOxKEO%sum_prod_op1d)
     CALL get_pqJL_OF_OpnD(Pq,JJ,LL,TWOxKEO%sum_prod_op1d(i))
     !write(out_unitp,*)
     !write(out_unitp,*) 'i (sum), Pq',i,Pq
     !write(out_unitp,*) 'i (sum), #Pq inact',i,count(Pq>nb_act)
     !CALL write_op(TWOxKEO%sum_prod_op1d(i))

     IF ( count(Pq>nb_act) > 0 ) THEN
       TWOxKEO%Cn(i)            = czero
       TWOxKEO%sum_prod_op1d(i) = czero
     END IF
     !CALL write_op(TWOxKEO%sum_prod_op1d(i))

   END DO
   CALL remove_opzero_in_F_sum_nd(TWOxKEO, 'TWOxKEO '// routine_name)

   IF (debug) THEN
     write(out_unitp,*) 'After removing inactive terms'
     CALL write_op(TWOxKEO,header=.TRUE.)
   END IF

   !write(out_unitp,*) 'Transfert inactive coef'

   ! calculation of Op1D for inactive coordinates and transfert it to Cn(i)
   ! Then, these Op1D are removed
   DO i = 1, size(TWOxKEO%sum_prod_op1d)
     DO k = 1, size(TWOxKEO%sum_prod_op1d(i)%prod_op1d)

       indexq = get_indexQ_OF_Op1D(TWOxKEO%sum_prod_op1d(i)%prod_op1d(k))

       IF (indexq > nb_act .AND. indexq <= nb_var) THEN
         !write(out_unitp,*) 'inactiv coord, i (sum),k (op1d)',i,k
         !write(out_unitp,*) 'inactiv coord, indexq',indexq
         CALL get_NumVal_Op1D(opval,Qval(indexq),TWOxKEO%sum_prod_op1d(i)%prod_op1d(k))
         !write(out_unitp,*) 'old Cn',i,k,TWOxKEO%Cn(i)

         !write(out_unitp,*) 'inactiv coord, opval',opval

         TWOxKEO%Cn(i) = TWOxKEO%Cn(i) * opval
         TWOxKEO%sum_prod_op1d(i)%prod_op1d(k) = cone ! IdOp
         !write(out_unitp,*) 'new Cn',i,k,TWOxKEO%Cn(i)

       END IF

     END DO
     !write(out_unitp,*)
     !write(out_unitp,*) 'i (sum)',i
     !CALL write_op(TWOxKEO%sum_prod_op1d(i))

     call Simplify_OpnD(TWOxKEO%sum_prod_op1d(i)) ! remove Id op
   END DO

   call Simplify_Sum_OpnD(TWOxKEO,Expand_Sin2=.TRUE.)


   IF (debug) THEN
     write(out_unitp,*) 'After new coef'
     CALL write_op(TWOxKEO,header=.TRUE.)
     write(out_unitp,*) ' END ',routine_name
   END IF

 END subroutine get_keo_for_Qactiv

   !! @description: Defines the rho(Q) associated to the volume element
   !! @param:       F_system     The vector/frame (recursive)
   !! @param:       rho          The output, rho
   RECURSIVE FUNCTION get_rho(F_system,alfa) RESULT (rho)
     USE mod_BunchPolyTransfo, only : Type_BFTransfo

     type(opnd)                                       :: rho
     type(Type_BFtransfo),        intent(in)          :: F_system
     TYPE(Frac_t),      optional, intent(in)          :: alfa

     TYPE(Frac_t)               :: alfa_loc
     integer                         :: i
     character (len =*), parameter   :: routine_name='get_rho'

     IF (present(alfa)) THEN
       alfa_loc = alfa
     ELSE
       alfa_loc = 1
     END IF

     IF ( compare_tab(F_system%euler, (/.FALSE., .FALSE., .FALSE./)) ) THEN ! true BF
       rho = get_rho_OF_Q(F_system%QVec(1),alfa_loc) *                 &
             get_rho_OF_Q(F_system%QVec(2),alfa_loc) *                 &
             get_rho_OF_Q(F_system%QVec(3),alfa_loc)
     ELSE
       rho = get_rho_OF_Q(F_system%QVec(1),alfa_loc) *                 &
             get_rho_OF_Q(F_system%QVec(2),alfa_loc) *                 &
             get_rho_OF_Q(F_system%QVec(3),alfa_loc) *                 &
             get_rho_OF_Q(F_system%QEuler(1),alfa_loc) *               &
             get_rho_OF_Q(F_system%QEuler(2),alfa_loc) *               &
             get_rho_OF_Q(F_system%QEuler(3),alfa_loc)
     END IF

     do i=1,F_system%nb_vect
       rho = rho * get_rho(F_system%tab_BFTransfo(i),alfa_loc)
     end do

   END FUNCTION get_rho
   !! @description: Defines the Jacobian (deformation part)
   !! @param:       F_system     The vector/frame (recursive)
   !! @param:       Jac          The output, Jacobian
   RECURSIVE FUNCTION get_Jac(F_system,alfa) RESULT (Jac)
     USE mod_BunchPolyTransfo, only : Type_BFTransfo

     type(Type_BFtransfo),        intent(in)          :: F_system
     type(opnd)                                       :: Jac
     TYPE(Frac_t),      optional, intent(in)          :: alfa

     TYPE(Frac_t)               :: alfa_loc
     integer                         :: i
     character (len =*), parameter   :: routine_name='get_Jac'

     IF (present(alfa)) THEN
       alfa_loc = alfa
     ELSE
       alfa_loc = 1
     END IF

     IF ( compare_tab(F_system%euler, (/.FALSE., .FALSE., .FALSE./)) ) THEN
       Jac = get_Jac_OF_Q(F_system%QVec(1),alfa_loc) *                 &
             get_Jac_OF_Q(F_system%QVec(2),alfa_loc) *                 &
             get_Jac_OF_Q(F_system%QVec(3),alfa_loc)
     ELSE
       Jac = get_Jac_OF_Q(F_system%QVec(1),alfa_loc) *                 &
             get_Jac_OF_Q(F_system%QVec(2),alfa_loc) *                 &
             get_Jac_OF_Q(F_system%QVec(3),alfa_loc) *                 &
             get_Jac_OF_Q(F_system%QEuler(1),alfa_loc) *               &
             get_Jac_OF_Q(F_system%QEuler(2),alfa_loc) *               &
             get_Jac_OF_Q(F_system%QEuler(3),alfa_loc)
     END IF

     do i=1,F_system%nb_vect
       Jac = Jac * get_Jac(F_system%tab_BFTransfo(i),alfa_loc)
     end do

   END FUNCTION get_Jac
   !! @description: Calculates the extra potential term and add to the keo
   !! @param:       keo       The kinetic energy operator (inout)
   SUBROUTINE add_Vextr_new(F_system,keo,tab_Qel,nrho,nb_Qel)
     USE mod_BunchPolyTransfo, only : Type_BFTransfo

     type(Type_BFtransfo), intent(in)                 :: F_system
     integer,              intent(in)                 :: nrho,nb_Qel
     type(sum_opnd),       intent(inout)              :: keo
     TYPE(opel),           intent(in)                 :: tab_Qel(:)

     type(opnd)                 :: Jac,rho
     type(opnd)                 :: x,x_inv
     TYPE(Frac_t)               :: Fp12,Fm12
     integer                    :: i

     SELECT CASE(nrho)
     CASE (0) ! rho=Jac => the KEO is Euclidean
       CONTINUE
     CASE (1) ! rho=1
       Fp12  = Frac_t( 1,2)
       Fm12  = Frac_t(-1,2)
       x     = CONE
       x_inv = CONE

       DO i=1,nb_Qel
         x     = x     * get_Jac_OF_Q(tab_Qel(i),Fm12)
         x_inv = x_inv * get_Jac_OF_Q(tab_Qel(i),Fp12)
       END DO

       !x     = get_Jac(F_system, Fm12) ! sqrt(1/Jac)
       !x_inv = get_Jac(F_system, Fp12) ! sqrt(Jac/1)

       keo = x_inv * keo * x
     CASE (2,3)
       Fp12 = Frac_t( 1,2)
       Fm12 = Frac_t(-1,2)
       x     = CONE
       x_inv = CONE

       DO i=1,nb_Qel
         x     = x     * get_rho_OF_Q(tab_Qel(i),Fp12) * get_Jac_OF_Q(tab_Qel(i),Fm12)
         x_inv = x_inv * get_rho_OF_Q(tab_Qel(i),Fm12) * get_Jac_OF_Q(tab_Qel(i),Fp12)
       END DO

       !x     = get_rho(F_system, Fp12) * get_Jac(F_system, Fm12) ! sqrt(rho/Jac)
       !x_inv = get_rho(F_system, Fm12) * get_Jac(F_system, Fp12) ! sqrt(Jac/rho)

       keo = x_inv * keo * x

     CASE DEFAULT
       STOP 'nrho /= 0,1,2,3'
     END SELECT

     write(out_unitp,*) 'x (sqrt(rho/jac)'
     CALL write_op(x)
     write(out_unitp,*) 'x_inv (sqrt(jac/rho)'
     CALL write_op(x_inv)

     call delete_op(x)
     call delete_op(x_inv)


   END SUBROUTINE add_Vextr_new
   SUBROUTINE Get_F2_F1_FROM_TWOxKEO(F_system,TWOxKEO,ExpandTWOxKEO,tabQact_Qel,nb_act,nb_var,nrho)
     USE mod_BunchPolyTransfo, only : Type_BFTransfo

     type(Type_BFtransfo), intent(in)                 :: F_system
     integer,              intent(in)                 :: nb_act,nb_var,nrho
     type(sum_opnd),       intent(in)                 :: TWOxKEO
     TYPE(opel),           intent(in)                 :: tabQact_Qel(nb_var)
     type(sum_opnd),       intent(inout)              :: ExpandTWOxKEO

     type(sum_opnd),       allocatable  :: Gana(:,:)
     type(sum_opnd),       allocatable  :: d1Gana(:,:) ! d./dQj Gana(i,j) !! same j for dQj and the index

     type(sum_opnd)                     :: Gij
     type(opnd)                         :: dQi

     integer :: i,indexQ
     integer                    :: pq(2),JJ(2),LL(2),iG,jG


     type(opnd)                 :: rho,rho_inv
     type(opnd), allocatable    :: d1lnrho(:)
     type(Sum_OF_op1d)          :: d1Op1D ! it will contain a sum of Op1D

     TYPE(Frac_t)          :: Fm1


     logical, parameter :: debug = .TRUE.
     !logical, parameter :: debug = .FALSE.
     character (len=*), parameter  :: routine_name='Get_F2_F1_FROM_TWOxKEO'

     IF (debug) THEN
       write(out_unitp,*) ' BEGINNING ',routine_name
       CALL write_op(TWOxKEO,header=.TRUE.)

       DO i=1,size(tabQact_Qel)
         CALL write_op(tabQact_Qel(i))
       END DO
       CALL flush_perso(out_unitp)
     END IF

     ! first extract G(ij)
     CALL alloc_NParray(Gana,(/nb_act+3,nb_act+3/),'Gana',routine_name)
     CALL C_TO_Mat_OF_sum_opnd(Gana,CZERO)

     DO i = 1, size(TWOxKEO%sum_prod_op1d)

       CALL get_pqJL_OF_OpnD(pq,JJ,LL,TWOxKEO%sum_prod_op1d(i))

       iG = 0
       jG = 0
       !write(out_unitp,*) 'term:',i
       IF (pq(1) > 0 .AND. pq(2) > 0) THEN ! def
         !write(out_unitp,*) 'def'
         iG = pq(1)
         jG = pq(2)
       ELSE IF (pq(1) > 0 .AND. JJ(1) > 0) THEN ! cor
         !write(out_unitp,*) 'cor'
         iG = pq(1)
         jG = JJ(1) -(nb_var-nb_act)
       ELSE IF (JJ(1) > 0 .AND. JJ(2) > 0) THEN ! rot
         !write(out_unitp,*) 'rot'
         iG = JJ(1) -(nb_var-nb_act)
         jG = JJ(2) -(nb_var-nb_act)
       ELSE IF (JJ(1) == 0 .AND. pq(1) > 0) THEN ! rot
         !write(out_unitp,*) 'pq^1'
         !CALL write_op(TWOxKEO%sum_prod_op1d(i),header=.TRUE.)
         iG = pq(1)
         jG = 0
       ELSE IF(JJ(1) > 0 .AND. pq(1) == 0) THEN ! rot
         !write(out_unitp,*) 'J^1'
         !CALL write_op(TWOxKEO%sum_prod_op1d(i),header=.TRUE.)
         iG = JJ(1) -(nb_var-nb_act)
         jG = 0
       END IF

       !write(out_unitp,*) i,'pq,JJ',pq,JJ
       !write(out_unitp,*) i,'iG,jG',iG,jG

       IF (iG > nb_act+3 .OR. jG > nb_act+3 .OR. iG < 0 .OR. jG < 0) THEN
         write(out_unitp,*) ' ERROR in ',routine_name
         write(out_unitp,*) ' iG or jG have a wrong range'
         write(out_unitp,*) ' iG, jG',iG,jG
         write(out_unitp,*) 'range: [1:',nb_act+3,'] or '
         write(out_unitp,*) 'iG = 0 and iG = 0 for the vep'
         write(out_unitp,*) 'CHECK the FORTRAN'
         STOP
       END IF

       IF (iG > 0 .AND. jG > 0) THEN

         Gij       = TWOxKEO%sum_prod_op1d(i)
         Gij%Cn(1) = TWOxKEO%Cn(i)
         !CALL write_op(Gij)

         CALL Change_PQ_OF_OpnD_TO_Id_OF_OnD(Gij%sum_prod_op1d(1))
         !CALL write_op(Gij)

         IF (iG == jG) THEN
           CALL F1_sum_nd_PLUS_TO_Fres_sum_nd(Gij,Gana(iG,jG)) !  adds Gij to Gana(iG,jG)
         ELSE  ! because we get both (iG,jG) and (jG,iG) elements in TWOxKEO
           Gij%Cn(1) = Gij%Cn(1) * CHALF

           CALL F1_sum_nd_PLUS_TO_Fres_sum_nd(Gij,Gana(iG,jG)) !  adds Gij to Gana(iG,jG)

           CALL F1_sum_nd_PLUS_TO_Fres_sum_nd(Gij,Gana(jG,iG)) !  adds Gij to Gana(jG,iG)

         END IF

       END IF
       !write(out_unitp,*) 'i,iG, jG',i,iG,jG
     END DO
     CALL delete_op(Gij)

     !DO iG=1,nb_act+3
     !DO jG=1,nb_act+3
     !  write(out_unitp,*) 'Gana,iG, jG',iG,jG
     !  CALL write_op(Gana(iG,jG))
     !END DO
     !END DO

     write(out_unitp,*) '================================================'
     write(out_unitp,*) '=========== der of Gana ========================'
     write(out_unitp,*) '================================================'


     ! 2d the 1st derivative of G(ij), just for the deformation
     CALL alloc_NParray(d1Gana,(/nb_act+3,nb_act+3/),'Gana',routine_name)
     CALL C_TO_Mat_OF_sum_opnd(d1Gana,CZERO)

     DO iG=1,nb_act
     DO jG=1,nb_act
       CALL Der1_OF_Sum_OpnD_TO_Sum_OpnD(Gana(iG,jG),jG,d1Gana(iG,jG))
       !write(out_unitp,*) 'Gana,iG, jG',iG,jG
       !CALL write_op(Gana(iG,jG))
       !write(out_unitp,*) 'd1Gana,iG, jG (der)',iG,jG
       !CALL write_op(d1Gana(iG,jG))
     END DO
     END DO

     write(out_unitp,*) '================================================'
     write(out_unitp,*) '=========== d1lnrho ============================'
     write(out_unitp,*) '================================================'

     CALL alloc_NParray(d1lnrho,(/nb_act/),'d1lnrho',routine_name)
     DO i=1,size(d1lnrho)
       d1lnrho(i) = CZERO
     END DO
     rho     = CONE
     rho_inv = CONE
     Fm1     = Frac_t(-1,1)
     SELECT CASE(nrho)
     CASE (0) ! rho=Jac => the KEO is Euclidean
        !rho         = get_Jac(F_system)
        !rho_inv     = get_Jac(F_system,Fm1)
        DO i=1,size(tabQact_Qel)
          rho     = rho     * get_Jac_OF_Q(tabQact_Qel(i))
          rho_inv = rho_inv * get_Jac_OF_Q(tabQact_Qel(i),Fm1)
        END DO

     CASE (1) ! rho=1
       rho     = cone
       rho_inv = cone

     CASE (2,3)
       !rho         = get_rho(F_system)
       !rho_inv     = get_rho(F_system,Fm1)
        DO i=1,size(tabQact_Qel)
          rho     = rho     * get_rho_OF_Q(tabQact_Qel(i))
          rho_inv = rho_inv * get_rho_OF_Q(tabQact_Qel(i),Fm1)
        END DO
     CASE DEFAULT
       STOP 'nrho /= 0,1,2,3'
     END SELECT


     DO i=1,size(rho%prod_op1d)
       indexQ = get_indexQ_OF_Op1D(rho%prod_op1d(i))
       d1Op1D = Der1_OF_d0Op1D(rho%prod_op1d(i))

       IF (size(d1Op1D%Sum_op1D) > 1) STOP 'size > 1 !!!'
       IF (indexQ > 0 .AND. indexQ <= nb_act) THEN
         d1lnrho(indexQ) = d1Op1D%Sum_op1D(1) * rho_inv%prod_op1d(i)
       END IF

     END DO

     write(out_unitp,*) 'rho'
     CALL write_op(rho)
     write(out_unitp,*) 'rho_inv'
     CALL write_op(rho_inv)
     write(out_unitp,*) 'd1lnrho'
     DO i=1,size(d1lnrho)
       write(out_unitp,*) 'd1lnrho',i
       CALL flush_perso(out_unitp)
       CALL write_op(d1lnrho(i))
     END DO


     write(out_unitp,*) '================================================'
     write(out_unitp,*) '=========== F2+F1 =============================='
     write(out_unitp,*) '================================================'
     ExpandTWOxKEO = CZERO
     DO iG=1,nb_act
     DO jG=1,nb_act
       Gij = Gana(iG,jG) * get_Pq(tabQact_Qel(iG)) * get_Pq(tabQact_Qel(jG))

       CALL F1_sum_nd_PLUS_TO_Fres_sum_nd(Gij,ExpandTWOxKEO) ! for iG /= jG, with the double sum, it is added twice (normal)

     END DO
     END DO

     DO iG=1,nb_act
     DO jG=1,nb_act
       dQi = get_Pq(tabQact_Qel(iG)) * get_Id(tabQact_Qel(iG),EYE)  ! d./dq = i Pq

       Gij = Gana(iG,jG) * d1lnrho(jG) * dQi
       CALL F1_sum_nd_MINUS_TO_Fres_sum_nd(Gij,ExpandTWOxKEO)

       Gij = d1Gana(iG,jG) * dQi
       CALL F1_sum_nd_MINUS_TO_Fres_sum_nd(Gij,ExpandTWOxKEO)

     END DO
     END DO

     call delete_op(dQi)
     call delete_op(rho_inv)
     call delete_op(rho)
     call delete_op(d1Op1D)
     CALL dealloc_NParray(d1lnrho,'d1lnrho',routine_name)
     call delete_op(Gij)
     CALL dealloc_NParray(Gana,'Gana',routine_name)
     CALL dealloc_NParray(d1Gana,'d1Gana',routine_name)

     IF (debug) THEN
       CALL write_op(ExpandTWOxKEO,header=.TRUE.)
       write(out_unitp,*) ' END ',routine_name
       CALL flush_perso(out_unitp)
     END IF

   END SUBROUTINE Get_F2_F1_FROM_TWOxKEO
   SUBROUTINE add_Vextr(keo)
     type(sum_opnd),       intent(inout)                :: keo

     type(sum_opnd)                 :: V_extr,RiRj_term

     integer                        :: i, j, k
     integer                        :: k1, k2, j1,j2
     logical                        :: l1_qact, l2_qact

     integer                        :: indexq1, indexq2
     complex (kind=Rkind)           :: coeff1, coeff2
     TYPE(Frac_t)                   :: alfa

     call init_to_opzero(F_sum_nd = V_extr)

     do i = 1, size(keo%sum_prod_op1d)
       indexq1 = 0
       indexq2 = 0
       do j = 1, size(keo%sum_prod_op1d(i)%prod_op1d)
         IF (get_idq_OF_Op1D(keo%sum_prod_op1d(i)%prod_op1d(j)) /= 2) CYCLE

         do k = 1, size(keo%sum_prod_op1d(i)%prod_op1d(j)%prod_opel)
           alfa = keo%sum_prod_op1d(i)%prod_op1d(j)%prod_opel(k)%alfa
           if(keo%sum_prod_op1d(i)%prod_op1d(j)%prod_opel(k)%idf == 4 .and. &
             & indexq1 == 0 .and. alfa == 1) then
             indexq1 = keo%sum_prod_op1d(i)%prod_op1d(j)%prod_opel(k)%indexq
             k1      = k
             j1      = j
             coeff1  = keo%sum_prod_op1d(i)%prod_op1d(j)%prod_opel(k)%coeff
           else if(keo%sum_prod_op1d(i)%prod_op1d(j)%prod_opel(k)%idf == 4 .and. &
             & indexq2 == 0 .and. alfa ==1) then
             indexq2 = keo%sum_prod_op1d(i)%prod_op1d(j)%prod_opel(k)%indexq
             k2      = k
             j2      = j
             coeff2  = keo%sum_prod_op1d(i)%prod_op1d(j)%prod_opel(k)%coeff

           end if
         end do
       end do

       if(indexq1 /= 0 .and. indexq2 /= 0 .and. indexq1 /= indexq2) then
         !write(out_unitp,*) 'indexq1=', indexq1, 'indexq2=', indexq2

         RiRj_term = keo%sum_prod_op1d(i)
         RiRj_term%Cn(1) = -keo%Cn(i)

         !!! substitue P_Q (idf=4) => 1/Q (idf=2) (just for R coordinates idq=2)
         RiRj_term%sum_prod_op1d(1)%prod_op1d(j1)%prod_opel(k1) = set_opel( &
                       idf = 2, idq = 2, alfa = -1, indexq = indexq1, coeff = coeff1)

         RiRj_term%sum_prod_op1d(1)%prod_op1d(j2)%prod_opel(k2) = set_opel( &
                       idf = 2, idq = 2, alfa = -1, indexq = indexq2, coeff = coeff2)

         CALL F1_sum_nd_PLUS_TO_Fres_sum_nd(RiRj_term,V_extr) !  adds RiRj_term to V_extr

       end if
     end do

     call F1_sum_nd_PLUS_TO_Fres_sum_nd(V_extr,KEO) !  adds V_extr to the KEO

     call delete_op(V_extr)
     call delete_op(RiRj_term)


   END SUBROUTINE add_Vextr
 END MODULE mod_Tana_op
