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
 module mod_Tana_NumKEO
 use mod_system
 USE mod_Tnum,     only : CoordType,Tnum
 USE mod_dnRho ! all

 USE mod_Tana_OpEl
 USE mod_Tana_Op1D
 USE mod_Tana_OpnD
 USE mod_Tana_sum_opnd
 USE mod_Tana_VecSumOpnD
 IMPLICIT NONE

 PRIVATE
 PUBLIC  :: get_NumG_WITH_AnaKEO, get_Numf2f1vep_WITH_AnaKEO

 CONTAINS

 !! @description: Write the KEO in a certain order
 !! @param:     F_sum_nd   Kinetic operator energy (type: sum_opnd).
 !! @param:     nvec       Number of vector in.
 !! @param:     index_q1   Index of the vecor which is
 !!                parallel to the BF-z frame.
  subroutine get_NumG_WITH_AnaKEO(TWOxKEO,Qval,mole,Gana,vep)

   type(sum_opnd),             intent(inout)             :: TWOxKEO
   real(kind=Rkind),           intent(inout)             :: Qval(:)
   real(kind=Rkind),           intent(inout)             :: Gana(:,:)
   real(kind=Rkind),           intent(inout)             :: vep
   TYPE (CoordType),           intent(in)                :: mole

   complex(kind=Rkind)        :: opval
   integer                    :: i,nb_act,nb_var
   integer                    :: error
   integer                    :: pq(2),JJ(2),LL(2),iG,jG

   logical, parameter :: debug = .FALSE.
   !logical, parameter :: debug = .TRUE.
   character (len=*), parameter :: routine_name='get_NumG_WITH_AnaKEO'

   IF (debug) THEN
     write(out_unitp,*) 'BEGINNING ',routine_name
     CALL write_op(TWOxKEO)
     CALL flush_perso(out_unitp)
   END IF


   IF( .NOT. allocated(TWOxKEO%sum_prod_op1d) ) THEN
     write(out_unitp,*) ' ERROR in ',routine_name
     write(out_unitp,*) ' TWOxKEO%sum_prod_op1d is not allocated'
     write(out_unitp,*) ' CHECK the fortran!!'
     STOP
   END IF

   nb_act = mole%nb_act
   nb_var = size(Qval)

   Gana(:,:) = zero
   vep       = ZERO

   DO i = 1, size(TWOxKEO%sum_prod_op1d)
     IF (debug) write(out_unitp,*)
     IF (debug) write(out_unitp,*) "==========================================="
     IF (debug) write(out_unitp,*) "============== term",i,"========================"

     CALL get_NumVal_OpnD(opval,Qval,TWOxKEO%sum_prod_op1d(i))
     opval = opval * TWOxKEO%Cn(i)

     IF (debug) write(out_unitp,*) 'term:',i,' TWOxKEO%Cn(i):',TWOxKEO%Cn(i)
     IF (debug) write(out_unitp,*) 'term:',i,' opval:',opval


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

     !write(out_unitp,*) 'pq,JJ',pq,JJ
     !write(out_unitp,*) 'iG,jG',iG,jG

     IF (iG > nb_act+3 .OR. jG > nb_act+3 .OR. iG < 0 .OR. jG < 0) THEN
       write(out_unitp,*) ' ERROR in ',routine_name
       write(out_unitp,*) ' iG or jG have a wrong range'
       write(out_unitp,*) ' iG, jG',iG,jG
       write(out_unitp,*) 'range: [1:',nb_act+3,'] or '
       write(out_unitp,*) 'iG = 0 and iG = 0 for the vep'
       write(out_unitp,*) 'CHECK the FORTRAN'
       STOP
     END IF

     IF (iG == 0 .AND. jG == 0) THEN ! vep
       !write(out_unitp,*) 'i,iG,jG (vep)',i,iG,jG
       !CALL write_op(TWOxKEO%sum_prod_op1d(i))
       vep = vep + HALF * real(opval,kind=Rkind)
     ELSE IF (iG > 0 .AND. jG > 0) THEN ! Gdef
       IF (iG /= jG) opval = opval * CHALF ! because we get both (iG,jG) and (jG,iG) elements
       Gana(iG,jG) = Gana(iG,jG) + real(opval,kind=Rkind)
       Gana(jG,iG) = Gana(iG,jG)

       IF (debug) write(out_unitp,*) 'i,iG, jG',i,iG,jG,Gana(iG,jG)

       !IF (iG /= jG .AND. iG <=nb_act .AND. jG <= nb_act) THEN
       !  write(out_unitp,*) ' G(iG,jG)',iG,jG,TWOxKEO%Cn(i)
       !  CALL write_op(TWOxKEO%sum_prod_op1d(i))
       !END IF

     END IF
     !write(out_unitp,*) 'i,iG, jG',i,iG,jG
   END DO

   IF (debug) THEN
     write(out_unitp,*) 'vep',vep
     write(out_unitp,*) 'G of Tana  '
     CALL write_mat(Gana,out_unitp,4)
     write(out_unitp,*) 'END ',routine_name
   END IF

 end subroutine get_NumG_WITH_AnaKEO

  subroutine get_Numf2f1vep_WITH_AnaKEO(TWOxKEO,Qval,mole,para_Tnum,f2,f1,vep,rho)
   USE mod_dnSVM

   type(sum_opnd),             intent(inout)             :: TWOxKEO
   real(kind=Rkind),           intent(inout)             :: Qval(:)
   real(kind=Rkind),           intent(inout)             :: f2(:,:),f1(:)
   real(kind=Rkind),           intent(inout)             :: vep,rho
   TYPE (CoordType),           intent(in)                :: mole
   TYPE (Tnum),                intent(in)                :: para_Tnum


   complex(kind=Rkind)        :: opval
   integer                    :: i,nb_var,nb_act
   integer                    :: error
   integer                    :: pq(2),JJ(2),LL(2),iG,jG,nb_J
   TYPE(Type_dnS)             :: dnrho,dnJac


   logical, parameter :: debug = .FALSE.
   !logical, parameter :: debug = .TRUE.
   character (len=*), parameter :: routine_name='get_Numf2f1vep_WITH_AnaKEO'

   IF (debug) THEN
     write(out_unitp,*) 'BEGINNING ',routine_name
     CALL write_op(TWOxKEO)
   END IF

   IF( .NOT. allocated(TWOxKEO%sum_prod_op1d) ) THEN
     write(out_unitp,*) ' ERROR in ',routine_name
     write(out_unitp,*) ' TWOxKEO%sum_prod_op1d is not allocated'
     write(out_unitp,*) ' CHECK the fortran!!'
     STOP
   END IF

   nb_act = mole%nb_act
   nb_var = size(Qval)

   CALL alloc_dnS(dnrho, nb_var_deriv=nb_act, nderiv=0)
   CALL alloc_dnS(dnJac, nb_var_deriv=nb_act, nderiv=0)

   CALL sub3_dnrho(dnrho,dnjac,Qval,mole,0,para_Tnum%num_x,para_Tnum%stepT,para_Tnum%nrho)
   rho = dnrho%d0

   CALL dealloc_dnS(dnrho)

   f2(:,:) = zero
   f1(:)   = zero
   vep     = ZERO

   DO i = 1, size(TWOxKEO%sum_prod_op1d)

     CALL get_NumVal_OpnD(opval,Qval,TWOxKEO%sum_prod_op1d(i))
     opval = opval * TWOxKEO%Cn(i)

     CALL get_pqJL_OF_OpnD(pq,JJ,LL,TWOxKEO%sum_prod_op1d(i))

     nb_J = count(JJ>0)
     iG = 0
     jG = 0
     IF (debug) write(out_unitp,*) 'term:',i,'pq',pq,'JJ',JJ,'LL',LL
     IF (pq(1) > 0 .AND. pq(2) > 0) THEN ! def
       IF (debug) write(out_unitp,*) 'def'
       iG = pq(1)
       jG = pq(2)
     ELSE IF (pq(1) > 0 .AND. JJ(1) > 0) THEN ! cor
       IF (debug) write(out_unitp,*) 'cor'
       iG = pq(1)
       jG = JJ(1) -(nb_var-nb_act)
     ELSE IF (JJ(1) > 0 .AND. JJ(2) > 0) THEN ! rot
       IF (debug) write(out_unitp,*) 'rot'
       iG = JJ(1) -(nb_var-nb_act)
       jG = JJ(2) -(nb_var-nb_act)
     ELSE IF (JJ(1) == 0 .AND. pq(1) > 0) THEN ! f1
       IF (debug) write(out_unitp,*) 'pq^1'
       iG = pq(1)
       jG = 0
     ELSE IF (JJ(2) == 0 .AND. pq(2) > 0) THEN ! f1
       IF (debug) write(out_unitp,*) 'pq^1'
       iG = pq(2)
       jG = 0
     ELSE IF(JJ(1) > 0 .AND. pq(1) == 0) THEN ! rot/cor
       IF (debug) write(out_unitp,*) 'J^1'
       iG = JJ(1) -(nb_var-nb_act)
       jG = 0
     END IF

     !IF (nb_J == 0) write(out_unitp,*) 'i (sum)',i,' Pq',pq
     !write(out_unitp,*) 'pq,JJ',pq,JJ
     !write(out_unitp,*) 'iG,jG',iG,jG

     IF (iG > nb_act+3 .OR. jG > nb_act+3 .OR. iG < 0 .OR. jG < 0) THEN
       write(out_unitp,*) ' ERROR in ',routine_name
       write(out_unitp,*) ' iG or jG have a wrong range'
       write(out_unitp,*) ' iG, jG',iG,jG
       write(out_unitp,*) 'range: [1:',nb_act+3,'] or '
       write(out_unitp,*) 'iG = 0 and iG = 0 for the vep'
       write(out_unitp,*) 'CHECK the FORTRAN'
       STOP
     END IF

     IF (iG == 0 .AND. jG == 0) THEN ! vep
       IF (debug)  write(out_unitp,*) 'add vep',iG
       vep = vep + HALF * real(opval,kind=Rkind)
     ELSE IF (iG > 0 .AND. jG > 0 .AND. nb_J == 0) THEN ! f2
       IF (debug)  write(out_unitp,*) 'add f2',iG
       ! the (-) is comming from Pq^2 = (-EYE d/dq)^2 = - d/dq ^2
       f2(iG,jG) = f2(iG,jG) -HALF * real(opval,kind=Rkind)
       f2(jG,iG) = f2(iG,jG)
     ELSE IF (iG > 0 .AND. jG == 0 .AND. nb_J == 0) THEN ! f1
       IF (debug) write(out_unitp,*) 'add f1',iG
       ! the (-EYE) is comming from Pq = -EYE d/dq
       f1(iG) = f1(iG) + HALF * real(-EYE*opval,kind=Rkind)
     END IF
   END DO

   IF (debug) THEN
     write(out_unitp,*) 'vep',vep
     write(out_unitp,*) 'f1 of Tana  '
     CALL write_vec(f1,out_unitp,4)
     write(out_unitp,*) 'f2 of Tana  '
     CALL write_mat(f2,out_unitp,4)
     write(out_unitp,*) 'END ',routine_name
   END IF

 end subroutine get_Numf2f1vep_WITH_AnaKEO

 end module mod_Tana_NumKEO
