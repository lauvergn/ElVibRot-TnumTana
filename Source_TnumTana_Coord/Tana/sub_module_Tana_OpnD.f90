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
 module mod_Tana_OpnD
 use mod_system
 use mod_Tana_OpEl ! all
 use mod_Tana_Op1D ! all
 IMPLICIT NONE
 PRIVATE

      !-----------------------------------------------------------!
      !                        OPND                               !
      !-----------------------------------------------------------!
      ! Product of 1d operators Type
      !! @description: Definition of a type of a product 1d operators (nd_operator)
      !!   like $ \hat{P}_{q_{k_1}} q_{k_2}^\alpha \cos^\alpha q_{k_2}, \,\,\,\, $
      !!    $ \sin^\alpha q_{k_1}\hat{P}_{q_{k_3}}, \,\,\,\, $ $ \cdots $.
      !!               This type is used for the analytical
      !!               computation of the KEO
      !! @param: prod_op1d        Array of 1d operators
      TYPE opnd
          type(op1d), allocatable     :: prod_op1d(:)
      END TYPE opnd

      INTERFACE alloc_NParray
        MODULE PROCEDURE alloc_NParray_OF_OpnDdim1
      END INTERFACE
      INTERFACE dealloc_NParray
        MODULE PROCEDURE dealloc_NParray_OF_OpnDdim1
      END INTERFACE
      INTERFACE check_NParray
        MODULE PROCEDURE check_NParray_OF_OpnDdim1
      END INTERFACE

  !!@description: Generic routine that check the allocated array of operator types
  interface check_allocate_op
    module procedure check_allocate_opnd
  end interface

  !!@description: Generic routine that deletes an array of operator types
  interface delete_op
    module procedure delete_opnd
  end interface

  !!@description: Generic routine that initializes a variable of operator type to zero 
  interface init_to_opzero
    module procedure init_opzero_opnd
  end interface
 
  !!@description: Generic routine that allocate variables of operator type  
  interface allocate_op
    module procedure allocate_opnd
  end interface

   INTERFACE write_op
     module procedure  write_opnd
   END INTERFACE

  !!@description: Generic routine that compares two nd-operators
  interface compare_op
    module procedure compare_F1_nd_and_F2_nd
  end interface


  !!@description: Generic routine that copy a operator F1 to another operator F2
  interface copy_F1_into_F2
    module procedure  copy_F1_nd_into_F2_nd, &
                      copy_F1_el_into_F2_nd, copy_F1_1d_into_F2_nd
  end interface

  INTERFACE assignment (=)
     MODULE PROCEDURE OpnD2_TO_OpnD1,Op1D2_TO_OpnD1,OpEl2_TO_OpnD1,R_TO_OpnD1,C_TO_OpnD1
  END INTERFACE

   INTERFACE get_F1_times_F2_to_F_nd
     module procedure get_F1_1d_times_F2_1d_to_Fres_nd, get_F1_1d_times_F2_nd_to_Fres_nd, &
                      get_F1_nd_times_F2_1d_to_Fres_nd, get_F1_nd_times_F2_nd_to_Fres_nd
   END INTERFACE

   INTERFACE operator (*)
      MODULE PROCEDURE F1_1d_times_F2_1d_TO_OpnD,F1_1d_times_F2_nd_TO_OpnD,F1_nd_times_F2_1d_TO_OpnD,F1_nd_times_F2_nd_TO_OpnD
   END INTERFACE

   PUBLIC  :: opnd, allocate_op, delete_op, check_allocate_op, write_op, compare_op
   PUBLIC  :: init_to_opzero, present_op_zero_in_F_nd, Set_coeff_OF_OpnD_TO_ONE, Simplify_OpnD
   PUBLIC  :: alloc_NParray, dealloc_NParray, check_NParray
   PUBLIC  :: get_F1_times_F2_to_F_nd, copy_F1_into_F2, assignment (=), operator (*)

   PUBLIC  :: get_coeff_OF_OpnD, get_sin, get_cos, get_Pq, get_Pq_dag, get_Id
   PUBLIC  :: get_Jac_OF_Q, get_rho_OF_Q, get_Q, get_Jx, get_Jy, get_Jz
   PUBLIC  :: get_Lx, get_Ly, get_Lz, get_zero, get_cot

   PUBLIC  :: Change_PQ_OF_OpnD_TO_Id_OF_OnD, Expand_OpnD_TO_SumOpnD
   PUBLIC  :: Export_Latex_Opnd, Export_Midas_Opnd
   PUBLIC  :: Export_MCTDH_Opnd, Export_VSCF_Opnd
   PUBLIC  :: get_NumVal_OpnD, get_pq_OF_OpnD, get_pqJL_OF_OpnD, set_indexQ_OF_OpnD



  contains

      SUBROUTINE check_NParray_OF_OpnDdim1(tab,name_var,name_sub)
      IMPLICIT NONE

      type(opnd), allocatable, intent(in) :: tab(:)

      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'check_NParray_OF_OpnDdim1'
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       IF (.NOT. allocated(tab))                                        &
                 CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       IF (size(tab) < 1) THEN
         write(out_unitp,*) ' ERROR in ',name_sub_alloc
         write(out_unitp,*) '   Called from ',name_sub
         write(out_unitp,*) ' Size of ',trim(name_var),' is wrong!!'
         STOP
       END IF

      END SUBROUTINE check_NParray_OF_OpnDdim1

      SUBROUTINE alloc_NParray_OF_OpnDdim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      type(opnd), allocatable, intent(inout) :: tab(:)
      integer,                 intent(in)    :: tab_ub(:)
      integer, optional,       intent(in)    :: tab_lb(:)
      character (len=*),       intent(in)    :: name_var,name_sub

      integer, parameter :: ndim=1
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_NParray_OF_OpnDdim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       CALL dealloc_NParray_OF_OpnDdim1(tab,name_var,name_sub_alloc // '  from ' // name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)

         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1)),stat=err_mem)
       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'OpnD')

      END SUBROUTINE alloc_NParray_OF_OpnDdim1
      SUBROUTINE dealloc_NParray_OF_OpnDdim1(tab,name_var,name_sub)
      IMPLICIT NONE

      type(opnd), allocatable, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

      integer :: i
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_NParray_OF_OpnDdim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       IF (.NOT. allocated(tab)) RETURN

       DO i=lbound(tab,dim=1),ubound(tab,dim=1)
         CALL delete_op(tab(i))
       END DO

       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'OpnD')

      END SUBROUTINE dealloc_NParray_OF_OpnDdim1


subroutine check_allocate_opnd(F_nd)

   type(opnd),           intent(in)    :: F_nd

   character (len=*), parameter :: routine_name="check_allocate_opnd"

   CALL check_NParray(F_nd%prod_op1d,'F_nd%prod_op1d',routine_name)

 end subroutine check_allocate_opnd

 !! @description: Deallocated a nd operator 
 !! @param:    F_nd    The 1d operator (type: opnd). 
 subroutine delete_opnd(F_nd)

   type(opnd),           intent(inout)    :: F_nd

   character (len=*), parameter :: routine_name="delete_opnd"

   if(allocated(F_nd%prod_op1d)) then
     CALL dealloc_NParray(F_nd%prod_op1d,'F_nd%prod_op1d',routine_name)
   end if

 end subroutine delete_opnd

 !! @description: Allocated a nd operator 
 !! @param:    F_nd    The nd operator (type: opnd). 
 !! @param:       ndim     Size of F_nd%prod_op1d. 
 subroutine allocate_opnd(F_nd, ndim)

   type(opnd),           intent(inout)    :: F_nd
   integer,              intent(in)       :: ndim

   character (len=*), parameter :: routine_name="allocate_opnd"

   CALL delete_opnd(F_nd)

   CALL alloc_NParray(F_nd%prod_op1d,(/ndim/),'F_nd%prod_op1d',routine_name)

 end subroutine allocate_opnd

 !! @description: Initialized a nd operator to zero
 !! @param:       F_nd    The nd operator (type: opnd). 
 subroutine init_opzero_opnd(F_nd)

   type(opnd),           intent(inout)    :: F_nd

   character (len=*), parameter :: routine_name="init_opzero_nd"

   call allocate_opnd(F_nd,1)
   call init_to_opzero(F_nd%prod_op1d(1))
 end subroutine init_opzero_opnd

 FUNCTION get_Id(Q_El,coeff)
   type(OpEl),           intent(in)               :: Q_El
   complex (kind=Rkind), intent(in), optional     :: coeff
   type(OpnD)                                     :: get_Id

   character (len=*), parameter   :: routine_name='get_Id'

   IF (present(coeff)) THEN
     get_Id = set_opel(1,Q_El%idq,1,Q_El%indexq,coeff)
   ELSE
     get_Id = set_opel(1,Q_El%idq,1,Q_El%indexq,cone)
   END IF

 END FUNCTION get_Id

 FUNCTION get_zero(Q_El)
   type(OpEl),      intent(in)       :: Q_El
   type(OpnD)                        :: get_zero

   integer             :: indexQ
   character (len=*), parameter   :: routine_name='get_Id'

   get_zero = set_opel(0,Q_El%idq,1,Q_El%indexq,cone)

 END FUNCTION get_zero

 FUNCTION get_Pq_dag(Q_El)
   type(OpEl),      intent(in)            :: Q_El

   type(OpnD)                        :: get_Pq_dag

   character (len=*), parameter   :: routine_name='get_Pq_dag'


     SELECT CASE (Q_El%idq) ! Pqdag = Pq
     CASE(1,-3,4,6,-7,8) ! cart, u=cos(th),phi,alpha, ub=cos(beta), gamma
       get_Pq_dag = set_opel(4,Q_El%idq,1,Q_El%indexq,cone)

     CASE(2) ! R (special case for  (R^-2 d/dR R^2)

       get_Pq_dag =  set_opel(2,Q_El%idq,-2,Q_El%indexq,cone) * &
                    (set_opel(4,Q_El%idq, 1,Q_El%indexq,cone) * &
                     set_opel(2,Q_El%idq, 2,Q_El%indexq,cone) )

       !get_Pq_dag = set_opel(4,Q_El%idq,1,Q_El%indexq,cone)

     CASE(3,7) ! th,beta   =>  Pqdag = sin^-1 Pq  Sin
       get_Pq_dag =  set_opel(6,Q_El%idq,-1,Q_El%indexq,cone) * &
                    (set_opel(4,Q_El%idq, 1,Q_El%indexq,cone) * &
                     set_opel(6,Q_El%idq, 1,Q_El%indexq,cone) )

     CASE DEFAULT
       write(out_unitp,*) 'ERROR in ',routine_name
       write(out_unitp,*) 'idq=',Q_El%idq
       write(out_unitp,*) "This idq does not exist"
       STOP
     END SELECT

 END FUNCTION get_Pq_dag
 FUNCTION get_Pq(Q_El)
   type(OpEl),      intent(in)       :: Q_El
   type(OpnD)                        :: get_Pq

   character (len=*), parameter   :: routine_name='get_Pq'

   SELECT CASE (Q_El%idq) ! Pq
   CASE(1,2,3,-3,4,6,7,-7,8) ! cart, R, th, u=cos(th),phi,alpha, beta, ub=cos(beta), gamma

     get_Pq = set_opel(4,Q_El%idq,1,Q_El%indexq,cone)

   CASE DEFAULT
     write(out_unitp,*) 'ERROR in ',routine_name
     write(out_unitp,*) 'idq=',Q_El%idq
     write(out_unitp,*) "This idq does not exist"
     STOP
   END SELECT

 END FUNCTION get_Pq
 FUNCTION get_Jx(Q_El)
   type(OpEl),      intent(in)       :: Q_El
   type(OpnD)                        :: get_Jx

   character (len=*), parameter   :: routine_name='get_Jx'

   get_Jx = set_opel(9,5,1,Q_El%indexq,cone)

 END FUNCTION get_Jx
 FUNCTION get_Jy(Q_El)
   type(OpEl),      intent(in)       :: Q_El
   type(OpnD)                        :: get_Jy

   character (len=*), parameter   :: routine_name='get_Jy'

   get_Jy = set_opel(10,5,1,Q_El%indexq,cone)

 END FUNCTION get_Jy
 FUNCTION get_Jz(Q_El)
   type(OpEl),      intent(in)       :: Q_El
   type(OpnD)                        :: get_Jz

   character (len=*), parameter   :: routine_name='get_Jz'

   get_Jz = set_opel(11,5,1,Q_El%indexq,cone)

 END FUNCTION get_Jz
 FUNCTION get_Lx(Q_El)
   type(OpEl),      intent(in)       :: Q_El
   type(OpnD)                        :: get_Lx

   character (len=*), parameter   :: routine_name='get_Lx'


   get_Lx = set_opel(idf=24,idq=5,alfa=1,indexq=Q_El%indexq,coeff=cone)

 END FUNCTION get_Lx
 FUNCTION get_Ly(Q_El)
   type(OpEl),      intent(in)       :: Q_El
   type(OpnD)                        :: get_Ly

   character (len=*), parameter   :: routine_name='get_Ly'

   get_Ly = set_opel(idf=25,idq=5,alfa=1,indexq=Q_El%indexq,coeff=cone)

 END FUNCTION get_Ly
 FUNCTION get_Lz(Q_El)
   type(OpEl),      intent(in)       :: Q_El
   type(OpnD)                        :: get_Lz

   character (len=*), parameter   :: routine_name='get_Lz'

   get_Lz = set_opel(idf=26,idq=5,alfa=1,indexq=Q_El%indexq,coeff=cone)

 END FUNCTION get_Lz
 FUNCTION get_Q(Q_El,alfa)
   type(OpEl),      intent(in)       :: Q_El
   integer, intent(in), optional     :: alfa
   type(OpnD)                        :: get_Q

   character (len=*), parameter   :: routine_name='get_Q'

   SELECT CASE (Q_El%idq)
   CASE(1,2,-3,-7) ! x,y,z,R,u,ub (u and ub should be removed)
     IF (present(alfa)) THEN
       get_Q = set_opel(2,Q_El%idq,alfa,Q_El%indexq,cone)
     ELSE
       get_Q = set_opel(2,Q_El%idq,1,Q_El%indexq,cone)
     END IF

   CASE DEFAULT
     write(out_unitp,*) 'ERROR in ',routine_name
     write(out_unitp,*) 'idq=',Q_El%idq
     write(out_unitp,*) "This idq does not exist"
     STOP
   END SELECT

 END FUNCTION get_Q
 FUNCTION get_cos(Q_El,alfa)
   type(OpEl),      intent(in)       :: Q_El
   integer, intent(in), optional     :: alfa
   type(OpnD)                        :: get_cos

   character (len=*), parameter   :: routine_name='get_cos'

   SELECT CASE (Q_El%idq)
   CASE(3,4,6,7,8) ! th ,phi,alpha, beta, gamma => cos(q)
     IF (present(alfa)) THEN
       get_cos = set_opel(5,Q_El%idq,alfa,Q_El%indexq,cone)
     ELSE
       get_cos = set_opel(5,Q_El%idq,1,Q_El%indexq,cone)
     END IF

   CASE(-3,-7) ! u=cos(th), ub=cos(beta) => q
     IF (present(alfa)) THEN
       get_cos = set_opel(2,Q_El%idq,alfa,Q_El%indexq,cone)
     ELSE
       get_cos = set_opel(2,Q_El%idq,1,Q_El%indexq,cone)
     END IF
   CASE DEFAULT
     write(out_unitp,*) 'ERROR in ',routine_name
     write(out_unitp,*) 'idq=',Q_El%idq
     write(out_unitp,*) "This idq does not exist"
     STOP
   END SELECT

 END FUNCTION get_cos
 FUNCTION get_cot(Q_El,alfa)
   type(OpEl),      intent(in)       :: Q_El
   integer, intent(in), optional     :: alfa
   type(OpnD)                        :: get_cot

   character (len=*), parameter   :: routine_name='get_cot'

   SELECT CASE (Q_El%idq)
   CASE(3,4,6,7,8) ! th ,phi,alpha, beta, gamma => cos(q)
     IF (present(alfa)) THEN
       get_cot = set_opel(5,Q_El%idq, alfa,Q_El%indexq,cone) *   &
                 set_opel(6,Q_El%idq,-alfa,Q_El%indexq,cone)
     ELSE
       get_cot = set_opel(5,Q_El%idq, 1,Q_El%indexq,cone) *      &
                 set_opel(6,Q_El%idq,-1,Q_El%indexq,cone)

     END IF

   CASE(-3,-7) ! u=cos(th), ub=cos(beta) => q
     IF (present(alfa)) THEN
       get_cot = set_opel(2,Q_El%idq, alfa,Q_El%indexq,cone) *   &
                 set_opel(3,Q_El%idq,-alfa,Q_El%indexq,cone)
     ELSE
       get_cot = set_opel(2,Q_El%idq, 1,Q_El%indexq,cone) *      &
                 set_opel(3,Q_El%idq,-1,Q_El%indexq,cone)
     END IF
   CASE DEFAULT
     write(out_unitp,*) 'ERROR in ',routine_name
     write(out_unitp,*) 'idq=',Q_El%idq
     write(out_unitp,*) "This idq does not exist"
     STOP
   END SELECT

 END FUNCTION get_cot
 FUNCTION get_sin(Q_El,alfa)
   type(OpEl),      intent(in)       :: Q_El
   integer, intent(in), optional     :: alfa
   type(OpnD)                        :: get_sin

   character (len=*), parameter   :: routine_name='get_sin'

   SELECT CASE (Q_El%idq) !
   CASE(3,4,6,7,8) ! th ,phi,alpha, beta, gamma => sin(q)
     IF (present(alfa)) THEN
       get_sin = set_opel(6,Q_El%idq,alfa,Q_El%indexq,cone)
     ELSE
       get_sin = set_opel(6,Q_El%idq,1,Q_El%indexq,cone)
     END IF

   CASE(-3,-7) ! u=cos(th), ub=cos(beta) =>  sqrt(1-q^2)
     IF (present(alfa)) THEN
       get_sin = set_opel(3,Q_El%idq,alfa,Q_El%indexq,cone)
     ELSE
       get_sin = set_opel(3,Q_El%idq,1,Q_El%indexq,cone)
     END IF
   CASE DEFAULT
     write(out_unitp,*) 'ERROR in ',routine_name
     write(out_unitp,*) 'idq=',Q_El%idq
     write(out_unitp,*) "This idq does not exist"
     STOP
   END SELECT

 END FUNCTION get_sin

 function get_rho_OF_Q(Q,alfa) RESULT(rho)
 type(opel), intent(in)                   :: Q
 type(OpnD)                               :: rho
 TYPE(FracInteger), optional, intent(in)  :: alfa

 TYPE(FracInteger)               :: alfa_loc
 integer :: idf
 character (len = *), parameter :: routine_name = 'get_rho_OF_Q'

  alfa_loc = 1
  SELECT CASE(Q%idq)
   CASE(1) ! x,y,z
     idf = 1

   CASE(2) ! R
     idf = 1

   CASE(-3) ! u=cos(th)
     idf = 1
   CASE(3) ! th
     idf = 6

   CASE(4) ! phi
     idf = 1

   CASE(6) ! alpha
     idf = 1

   CASE(7) ! beta
     idf = 6
   CASE(-7) ! ub=cos(beta)
     idf = 1

   CASE(8) ! gamma
     idf = 1

   CASE DEFAULT ! if a Q is not initialized => rho = 1
      idf = 1
   END SELECT

   IF (present(alfa)) THEN
     alfa_loc = alfa_loc * alfa
   END IF

   rho = set_opel(idf, Q%idq, alfa_loc%num, Q%indexq, cone, alfa_den=alfa_loc%den)

 end function get_rho_OF_Q
 function get_Jac_OF_Q(Q,alfa) RESULT(jac)
 type(opel), intent(in)                   :: Q
 type(OpnD)                               :: Jac
 TYPE(FracInteger), optional, intent(in)  :: alfa

 TYPE(FracInteger)               :: alfa_loc

 integer :: idf
 character (len = *), parameter :: routine_name = 'get_Jac_OF_Q'

  alfa_loc = 1
  SELECT CASE(Q%idq)
   CASE(1) ! x,y,z
     idf = 1

   CASE(2) ! R
     idf  = 2
     alfa_loc = 2

   CASE(-3) ! u=cos(th)
     idf = 1
   CASE(3) ! th
     idf = 6

   CASE(4) ! phi
     idf = 1

   CASE(6) ! alpha
     idf = 1

   CASE(7) ! beta
     idf = 6
   CASE(-7) ! ub=cos(beta)
     idf = 1

   CASE(8) ! gamma
     idf = 1

   CASE DEFAULT ! if a Q is not initialized => Jac = 1
     idf = 1
   END SELECT

   IF (present(alfa)) THEN
     alfa_loc = alfa_loc * alfa
   END IF

   Jac = set_opel(idf, Q%idq, alfa_loc%num, Q%indexq, cone, alfa_den=alfa_loc%den)

 end function get_Jac_OF_Q

 subroutine Export_Latex_Opnd(Fnd,tab_Qname,FndName)
   type(opnd),                       intent(in)      :: Fnd
   character (len = :), allocatable, intent(inout)     :: FndName
   character(len=*),                 intent(in)       :: tab_Qname(:)

   ! local variables
   character (len = :), allocatable     :: qname,FndName_loc
   character (len = :), allocatable     :: F1dName
   integer :: j,m
   character (len = *), parameter       :: mult = ' '

   character (len = *), parameter :: routine_name = 'Export_Latex_Opnd'


   !CALL write_op(Fnd)
   FndName_loc = String_TO_String('')

   IF (size(Fnd%prod_op1d) > 0) THEN
     DO j=1,size(Fnd%prod_op1d)
       m = get_indexQ_OF_Op1D( Fnd%prod_op1d(j) )
       if (m <= size(tab_Qname) .AND. m > 0) then
         qname = String_TO_String(trim(tab_Qname(m)))
       else
         qname = String_TO_String(trim(Fnd%prod_op1d(j)%prod_opel(1)%opname))
       end if

       CALL Export_Latex_Op1D(Fnd%prod_op1d(j),qname,F1dName)
       FndName_loc = String_TO_String( FndName_loc // mult // F1dName)
     END DO
     IF (allocated(F1dName)) deallocate(F1dName)

   ELSE
     FndName_loc = String_TO_String('')
   END IF

   FndName = FndName_loc

   IF (allocated(FndName_loc))   deallocate(FndName_loc)
   IF (allocated(qname))         deallocate(qname)

 end subroutine Export_Latex_Opnd

 subroutine Export_Midas_Opnd(Fnd, tab_Qname, FndName)
   type(opnd),                       intent(in)      :: Fnd
   character (len = :), allocatable, intent(inout)   :: FndName
   character(len=*),                 intent(in)      :: tab_Qname(:)

   ! local variables
   character (len = :), allocatable     :: qname, Qdispname !Emil new
   character (len = :), allocatable     :: F1dName,FndName_loc
   integer :: j,m
   character (len = *), parameter       :: mult = ' '
   character (len = Name_len)           :: cindexq

   character (len = *), parameter       :: routine_name = 'Export_Midas_Opnd'

   !CALL write_op(Fnd)
   FndName_loc = String_TO_String('')

   IF (size(Fnd%prod_op1d) > 0) THEN
     DO j = 1, size(Fnd%prod_op1d)
       m = get_indexQ_OF_Op1D( Fnd%prod_op1d(j) ) -1

       Qname = String_TO_String('(Q' // int_TO_char(m) // ')')
       Qdispname = String_TO_String('Q+Q' // int_TO_char(m) // '_ref') !Emil new

       CALL Export_Midas_Op1D(Fnd%prod_op1d(j), Qdispname, F1dName) !Emil change
       !CALL Export_Midas_Op1D(Fnd%prod_op1d(j), 'Q', F1dName) ! old one DML
       FndName_loc = String_TO_String( FndName_loc // mult // F1dName // Qname)
     END DO
     IF (allocated(F1dName)) deallocate(F1dName)
     IF (allocated(Qdispname)) deallocate(Qdispname) !Emil new

   ELSE
     FndName_loc = String_TO_String('')
   END IF

   FndName = FndName_loc

   IF (allocated(FndName_loc))   deallocate(FndName_loc)
   IF (allocated(qname))         deallocate(qname)

 end subroutine Export_Midas_Opnd

 subroutine Export_MCTDH_Opnd(Fnd,FndName,nb_act)
   type(opnd),                       intent(inout)   :: Fnd
   integer,                          intent(in)      :: nb_act
   character (len = :), allocatable, intent(inout)   :: FndName


   ! local variables
   character (len = :), allocatable     :: F1dName,FndName_loc
   character (len = :), allocatable     :: SepCoord

   integer :: j,m,idq
   logical :: First_idq5

   character (len = *), parameter :: routine_name = 'Export_MCTDH_Opnd'


   !CALL write_op(Fnd)
   FndName_loc = String_TO_String('')

   First_idq5 = .TRUE.
   IF (size(Fnd%prod_op1d) > 0) THEN
     DO j=1,size(Fnd%prod_op1d)
       idq = get_idq_OF_Op1D( Fnd%prod_op1d(j) )
       m = get_indexQ_OF_Op1D( Fnd%prod_op1d(j) )
       IF (idq == 5) THEN ! Jx,Jy,Jz
         m = nb_act+1
         IF (First_idq5) THEN
           SepCoord = String_TO_String(' |' // int_TO_char(m) // '   ',ltrim=.FALSE.)
           First_idq5 = .FALSE.
         ELSE
           SepCoord = String_TO_String('*')
         END IF
       ELSE
         SepCoord = String_TO_String(' |' // int_TO_char(m) // '   ',ltrim=.FALSE.)
       END IF

       CALL Export_MCTDH_Op1D(Fnd%prod_op1d(j),F1dName)
       FndName_loc = String_TO_String( FndName_loc // SepCoord // F1dName)
     END DO
     IF (allocated(F1dName))    deallocate(F1dName)
     IF (allocated(SepCoord))   deallocate(SepCoord)

   END IF

   FndName = FndName_loc

   IF (allocated(FndName_loc))   deallocate(FndName_loc)

 end subroutine Export_MCTDH_Opnd

 subroutine Export_VSCF_Opnd(Fnd,tab_Qname,FndName)
   type(opnd),                       intent(in)      :: Fnd
   character (len = :), allocatable, intent(inout)   :: FndName
   character(len=*),                 intent(in)      :: tab_Qname(:)


   !local variables
   character (len = :), allocatable     :: qname,FndName_loc
   character (len = :), allocatable     :: F1dName
   integer :: j,m
   character (len = *), parameter       :: mult = ' '

   character (len = *), parameter :: routine_name = 'Export_VSCF_Opnd'


   !CALL write_op(Fnd)
   FndName_loc = String_TO_String('')

   IF (size(Fnd%prod_op1d) > 0) THEN
     DO j=1,size(Fnd%prod_op1d)
       m = get_indexQ_OF_Op1D( Fnd%prod_op1d(j) )
       qname = String_TO_String('Q' // int_TO_char(m) )

       CALL Export_VSCF_Op1D(Fnd%prod_op1d(j),qname,F1dName)
       FndName_loc = String_TO_String( FndName_loc // mult // F1dName)
     END DO
     IF (allocated(F1dName)) deallocate(F1dName)

   ELSE
     FndName_loc = String_TO_String('')
   END IF

   FndName = FndName_loc

   IF (allocated(FndName_loc))   deallocate(FndName_loc)
   IF (allocated(qname))         deallocate(qname)


 end subroutine Export_VSCF_Opnd

   !! @description: Write an array of nd operators,
   !! @param:       F_nd      The operator (type: opnd).
   !! @param:       header    If present write comment in the biginning of the
   !!                         of the file
   SUBROUTINE write_opnd(F_nd,  i_file, header, append, close_file)
     type(opnd),                intent(in)       :: F_nd
     integer, optional,         intent(in)       :: i_file
     logical, optional,         intent(in)       :: header
     logical, optional,         intent(in)       :: append
     logical, optional,         intent(in)       :: close_file

     integer                        :: error
     integer                        :: i
     integer                        :: i_open
     logical                        :: header_loc
     integer             :: pq(2),J(2),L(2)

     character (len=*), parameter :: routine_name='write_opnd'

     header_loc = .FALSE.
     if (present(header)) header_loc = header

     if(present(i_file)) then
       i_open = i_file
     else
       i_open = out_unitp
     end if

     if (header_loc) then
       write(i_open, *)
       write(i_open, *)
       write(i_open, '(40x, A)') '========== writing a nd operator ======== '
       write(i_open, *)
       write(i_open, '(A)')  "   idf          iqd        alfa       indexq          op_name                 coeff"
     end if

     CALL get_pqJL_OF_OpnD(pq,J,L,F_nd)
     IF (pq(1) == 7 .AND. pq(2) == 8 .OR. pq(1) == 8 .AND. pq(2) == 7 .OR. &
         pq(1) == 0 .AND. pq(2) == 8 .OR. pq(1) == 8 .AND. pq(2) == 0) THEN
     do i = 1, size(F_nd%prod_op1d)
       call write_op(F_nd%prod_op1d(i), i_open)
     end do
     END IF
   END SUBROUTINE write_opnd

 !! @description: Copy a opnd operator F1_nd to 
 !!               another sum_opnd operator F2_nd
 !! @param:     F1_nd    The operator which will be copied 
 !! @param:    F2_nd    The operator in which F1_sum_nd will be copied 
 subroutine copy_F1_nd_into_F2_nd(F1_nd, F2_nd)

   type(opnd),           intent(in)    :: F1_nd
   type(opnd),           intent(inout) :: F2_nd

   integer                    :: i
   integer                    :: ndim_op1d, ndim_opel
   character (len=*), parameter :: routine_name="copy_F1_nd_into_F2_nd"

   call delete_op(F2_nd)
   ndim_op1d = size(F1_nd%prod_op1d)
   call allocate_op(F2_nd,ndim_op1d)
   do i = 1, ndim_op1d
     call copy_F1_into_F2(F1_nd%prod_op1d(i), F2_nd%prod_op1d(i))
   end do
 end subroutine copy_F1_nd_into_F2_nd

 subroutine OpnD2_TO_OpnD1(OpnD1,OpnD2)

   type(opnd),           intent(in)    :: OpnD2
   type(opnd),           intent(inout) :: OpnD1

   integer                    :: i
   character (len=*), parameter :: routine_name="OpnD2_TO_OpnD1"

   call allocate_op(OpnD1, size(OpnD2%prod_op1d) )

   do i = 1, size(OpnD2%prod_op1d)
     OpnD1%prod_op1d(i) = OpnD2%prod_op1d(i)
   end do

 end subroutine OpnD2_TO_OpnD1

 !! @description: Copy a 1d operator F1_1d to 
 !!               another a 1d operator F2_nd
 !! @param:     F1_1d    The operator which will be copied 
 !! @param:    F2_nd    The operator in which F1_1d will be copied 
 subroutine copy_F1_1d_into_F2_nd(F1_1d, F2_nd)

   type(op1d),           intent(in)    :: F1_1d
   type(opnd),           intent(inout) :: F2_nd

   character (len=*), parameter :: routine_name="copy_F1_1d_into_F2_nd"

   call delete_op(F2_nd)
   call allocate_opnd(F2_nd,1)
   call copy_F1_into_F2(F1_1d, F2_nd%prod_op1d(1))
 end subroutine copy_F1_1d_into_F2_nd
 subroutine Op1D2_TO_OpnD1(OpnD1,Op1D2)

   type(op1d),           intent(in)    :: Op1D2
   type(opnd),           intent(inout) :: OpnD1

   character (len=*), parameter :: routine_name="Op1D2_TO_OpnD1"

   call allocate_op(OpnD1, 1 )

   OpnD1%prod_op1d(1) = Op1D2

 end subroutine Op1D2_TO_OpnD1

 !! @description: Copy an elementary operator F1_el to 
 !!               another a nd operator F2_nd
 !! @param:     F1_el    The operator which will be copied 
 !! @param:    F2_nd    The operator in which F1_el will be copied 
 subroutine copy_F1_el_into_F2_nd(F1_el, F2_nd)

   type(opel),           intent(in)    :: F1_el
   type(opnd),           intent(inout) :: F2_nd

   character (len=*), parameter :: routine_name="copy_F1_el_into_F2_nd"

   call delete_op(F2_nd)
   call allocate_opnd(F2_nd,1)
   call copy_F1_into_F2(F1_el, F2_nd%prod_op1d(1))
 end subroutine copy_F1_el_into_F2_nd

 subroutine OpEl2_TO_OpnD1(OpnD1,OpEl2)

   type(opEl),           intent(in)    :: OpEl2
   type(opnd),           intent(inout) :: OpnD1

   character (len=*), parameter :: routine_name="OpEl2_TO_OpnD1"

   call allocate_op(OpnD1, 1 )

   OpnD1%prod_op1d(1) = OpEl2

 end subroutine OpEl2_TO_OpnD1

 subroutine R_TO_OpnD1(OpnD1,R)

   real (kind=Rkind),    intent(in)    :: R
   type(opnd),           intent(inout) :: OpnD1

   character (len=*), parameter :: routine_name="R_TO_OpnD1"

   call allocate_op(OpnD1, 1 )

   OpnD1%prod_op1d(1) = R

 end subroutine R_TO_OpnD1
 subroutine C_TO_OpnD1(OpnD1,C)

   complex (kind=Rkind),    intent(in) :: C
   type(opnd),           intent(inout) :: OpnD1

   character (len=*), parameter :: routine_name="C_TO_OpnD1"

   call allocate_op(OpnD1, 1 )

   OpnD1%prod_op1d(1) = C

 end subroutine C_TO_OpnD1
   !! @description: Does the product of two 1d operators.,
   !!               the indexq of the two operator can be
   !!               different. The result is saved the result in Fres_nd.
   !! @param.in:    F1_1d First  1d operator (type: op1d).
   !! @param.in:    F2_1d second 1d operator (type: op1d).
   !! @param.out:   Fres_nd  array of nd operators
   !!               the result is saved (type: opnd).

   SUBROUTINE get_F1_1d_times_F2_1d_to_Fres_nd(F1_1d, F2_1d, Fres_nd)
     type(op1d),       intent(in)       :: F1_1d
     type(op1d),       intent(in)       :: F2_1d
     type(opnd),       intent(inout)    :: Fres_nd

     integer                    :: iF1_opzero, iF2_opzero

     !logical, parameter           :: debug=.TRUE.
     logical, parameter           :: debug=.FALSE.
     character (len=*), parameter :: routine_name='get_F1_1d_times_F2_1d_to_Fres_nd'

   IF (debug) THEN
     write(out_unitp,*) 'BEGINNING ',routine_name
     CALL write_op(F1_1d,header=.TRUE.)
     CALL write_op(F2_1d,header=.TRUE.)
   END IF

   CALL check_allocate_op(F1_1d)
   CALL check_allocate_op(F2_1d)

   call present_op_zero_in_F_1d(F1_1d, iF1_opzero,'F1_1d from '//routine_name)
   if(iF1_opzero /= -1 ) then
     Fres_nd = F1_1d%prod_opel(iF1_opzero)
     IF (debug) THEN
       CALL write_op(Fres_nd,header=.TRUE.)
       write(out_unitp,*) 'END ',routine_name
     END IF
     return
   end  if

   call present_op_zero_in_F_1d(F2_1d, iF2_opzero,'F2_1d from '//routine_name)
   if(iF2_opzero /= -1 ) then
     Fres_nd = F2_1d%prod_opel(iF2_opzero)
     IF (debug) THEN
       CALL write_op(Fres_nd,header=.TRUE.)
       write(out_unitp,*) 'END ',routine_name
     END IF
     return
   end  if

   call check_idq_in_F_1d(F1_1d, 'F1_1d from '//routine_name)
   call check_idq_in_F_1d(F2_1d, 'F2_1d from '//routine_name)

   if(compare_indexq(F1_1d, F2_1d)) then
     call allocate_op(Fres_nd, 1)
     Fres_nd%prod_op1d(1) = F1_1d_times_F2_1d(F1_1d , F2_1d)
   else
     call allocate_op(Fres_nd, 2)
     Fres_nd%prod_op1d(1) = F1_1d
     Fres_nd%prod_op1d(2) = F2_1d
   end if

   call Simplify_OpnD(Fres_nd)

   IF (debug) THEN
     CALL write_op(Fres_nd,header=.TRUE.)
     write(out_unitp,*) 'END ',routine_name
   END IF

   END SUBROUTINE get_F1_1d_times_F2_1d_to_Fres_nd

   function F1_1d_times_F2_1d_TO_OpnD(F1_1d, F2_1d) result(Fres_nD)
     type(opnd)    :: Fres_nD
     type(op1d),       intent(in)       :: F1_1d
     type(op1d),       intent(in)       :: F2_1d

     integer                    :: i_opzero
     integer                    :: ndim1, ndim2

   character (len=*), parameter :: routine_name='F1_1d_times_F2_1d_TO_OpnD'

   CALL get_F1_1d_times_F2_1d_to_Fres_nd(F1_1d, F2_1d, Fres_nd)

   end function F1_1d_times_F2_1d_TO_OpnD

   function F1_1d_times_F2_nd_TO_OpnD(F1_1d, F2_nd) result(Fres_nD)
     type(opnd)    :: Fres_nD
     type(op1d),       intent(in)       :: F1_1d
     type(opnd),       intent(in)       :: F2_nd

   character (len=*), parameter :: routine_name='F1_1d_times_F2_nd_TO_OpnD'


   CALL get_F1_1d_times_F2_nd_to_Fres_nd(F1_1d, F2_nd, Fres_nd)

   end function F1_1d_times_F2_nd_TO_OpnD
   function F1_nd_times_F2_1d_TO_OpnD(F1_nd, F2_1d) result(Fres_nD)
     type(opnd)    :: Fres_nD
     type(op1d),       intent(in)       :: F2_1d
     type(opnd),       intent(in)       :: F1_nd

   character (len=*), parameter :: routine_name='F1_nd_times_F2_1d_TO_OpnD'

   CALL get_F1_nd_times_F2_1d_to_Fres_nd(F1_nd, F2_1d, Fres_nd)

   end function F1_nd_times_F2_1d_TO_OpnD
   function F1_nd_times_F2_nd_TO_OpnD(F1_nd, F2_nd) result(Fres_nD)
     type(opnd)    :: Fres_nD
     type(opnd),       intent(in)       :: F2_nd
     type(opnd),       intent(in)       :: F1_nd

   character (len=*), parameter :: routine_name='F1_nd_times_F2_nd_TO_OpnD'

   CALL get_F1_nd_times_F2_nd_to_Fres_nd(F1_nd, F2_nd, Fres_nd)

   end function F1_nd_times_F2_nd_TO_OpnD

   !! @description: Does the product of an 1d operator F1_1d with
   !!               an nd operator F2_nd
   !!               and save the result in Fres_nd.
   !! @param.in:    F1_1d the 1d operator (type: op1d).
   !! @param.in:    F2_nd_the nd operator (type: opnd).
   !! @param.out:   Fres_nd  Array of nd operators in which
   !!               the result is saved (type: opnd).
   SUBROUTINE get_F1_1d_times_F2_nd_to_Fres_nd(F1_1d, F2_nd, Fres_nd)
     type(op1d),       intent(in)       :: F1_1d
     type(opnd),       intent(in)       :: F2_nd
     type(opnd),       intent(inout)    :: Fres_nd

     integer                    :: error
     integer                    :: iF1_opzero
     integer                    :: iF2_opzero, jF2_opzero
     integer                    :: i, p
     !logical, parameter           :: debug=.TRUE.
     logical, parameter           :: debug=.FALSE.
     character (len=*), parameter :: routine_name='get_F1_1d_times_F2_nd_to_Fres_nd'

   IF (debug) THEN
     write(out_unitp,*) 'BEGINNING ',routine_name
     CALL write_op(F1_1d,header=.TRUE.)
     CALL write_op(F2_nd,header=.TRUE.)
   END IF

   CALL check_allocate_op(F1_1d)
   CALL check_allocate_op(F2_nd)


   call present_op_zero_in_F_1d(F1_1d, iF1_opzero,'F1_1d from '//routine_name)
   if(iF1_opzero /= -1 ) then
     Fres_nd = F1_1d%prod_opel(iF1_opzero)
     IF (debug) THEN
       CALL write_op(Fres_nd,header=.TRUE.)
       write(out_unitp,*) 'END ',routine_name
     END IF
     return
   end  if

   call present_op_zero_in_F_nd(F2_nd, iF2_opzero,jF2_opzero, 'F2_nd from '//routine_name)
   if(jF2_opzero /= -1 ) then
     Fres_nd = F2_nd%prod_op1d(iF2_opzero)%prod_opel(jF2_opzero)

     IF (debug) THEN
       CALL write_op(Fres_nd,header=.TRUE.)
       write(out_unitp,*) 'END ',routine_name
     END IF
     return
   end if

   call check_idq_in_F_1d(F1_1d, 'F1_1d from '//routine_name)
   do i = 1, size(F2_nd%prod_op1d)
     call check_idq_in_F_1d(F2_nd%prod_op1d(i), 'F2_nd from '//routine_name)
   end do
   p = 0
   do i = 1, size(F2_nd%prod_op1d)
     if(compare_indexq(F1_1d, F2_nd%prod_op1d(i))) then
       p = i
       exit
     end if
   end do


   if(p /= 0) then
     call allocate_op(Fres_nd, size(F2_nd%prod_op1d))
     do i = 1, size(F2_nd%prod_op1d)
       if (i == p) then
         Fres_nd%prod_op1d(i) = F1_1d_times_F2_1d( F1_1d , F2_nd%prod_op1d(i) )
       else
         Fres_nd%prod_op1d(i) = F2_nd%prod_op1d(i)
       end if
     end do
   else
     call allocate_op(Fres_nd, size(F2_nd%prod_op1d)+1)

     Fres_nd%prod_op1d(1) = F1_1d

     do i = 1, size(F2_nd%prod_op1d)
       Fres_nd%prod_op1d(i+1) = F2_nd%prod_op1d(i)
     end do

   end if

   call Simplify_OpnD(Fres_nd)

   IF (debug) THEN
     CALL write_op(Fres_nd,header=.TRUE.)
     write(out_unitp,*) 'END ',routine_name
   END IF

   END SUBROUTINE get_F1_1d_times_F2_nd_to_Fres_nd

   !! @description: Does the product of an nd operator F1_nd with
   !!               an 1d operator F2_1d
   !!               and save the result in Fres_nd.
   !! @param.in:    F1_nd the 1d operator (type: opnd).
   !! @param.in:    F2_1d_the nd operator (type: op1d).
   !! @param.out:   Fres_nd  Array of nd operators in which
   !!               the result is saved (type: opnd).
   SUBROUTINE get_F1_nd_times_F2_1d_to_Fres_nd(F1_nd, F2_1d, Fres_nd)
     type(opnd),       intent(in)       :: F1_nd
     type(op1d),       intent(in)       :: F2_1d
     type(opnd),       intent(inout)    :: Fres_nd

     integer                    :: iF1_opzero, jF1_opzero
     integer                    :: iF2_opzero
     integer                    :: i
     integer                    :: p
     !logical, parameter           :: debug=.TRUE.
     logical, parameter           :: debug=.FALSE.
     character (len=*), parameter :: routine_name='get_F1_nd_times_F2_1d_to_Fres_nd'

     IF (debug) THEN
       write(out_unitp,*) ' BEGINNING ',routine_name
     CALL write_op(F1_nd,header=.TRUE.)
     CALL write_op(F2_1d,header=.TRUE.)
       CALL flush_perso(out_unitp)
     END IF

     CALL check_allocate_op(F1_nd)
     CALL check_allocate_op(F2_1d)

     call present_op_zero_in_F_1d(F2_1d, iF2_opzero,'F2_1d from '//routine_name)
     if(iF2_opzero /= -1 ) then
       Fres_nd = F2_1d%prod_opel(iF2_opzero)
       IF (debug) THEN
         CALL write_op(Fres_nd,header=.TRUE.)
         write(out_unitp,*) 'END ',routine_name
         CALL flush_perso(out_unitp)
       END IF
       return
     end  if
     call present_op_zero_in_F_nd(F1_nd, iF1_opzero,jF1_opzero, 'F1_nd from '//routine_name)
     if(jF1_opzero /= -1 ) then
       Fres_nd = F1_nd%prod_op1d(iF1_opzero)%prod_opel(jF1_opzero)
       IF (debug) THEN
         CALL write_op(Fres_nd,header=.TRUE.)
         write(out_unitp,*) 'END ',routine_name
         CALL flush_perso(out_unitp)
       END IF
       return
     end if

     call check_idq_in_F_1d(F2_1d, 'F2_1d from '//routine_name)
     do i = 1, size(F1_nd%prod_op1d)
       call check_idq_in_F_1d(F1_nd%prod_op1d(i), 'F1_nd from '//routine_name)
     end do

     p = 0
     do i = 1, size(F1_nd%prod_op1d)
       if(compare_indexq(F1_nd%prod_op1d(i), F2_1d)) then
         p = i
         exit
       end if
     end do

     if(p /= 0) then
       call allocate_op(Fres_nd, size(F1_nd%prod_op1d))
       do i = 1, size(F1_nd%prod_op1d)
         if (i == p) then
           Fres_nd%prod_op1d(i) = F1_1d_times_F2_1d( F1_nd%prod_op1d(i) , F2_1d )
         else
           Fres_nd%prod_op1d(i) = F1_nd%prod_op1d(i)
         end if
       end do
     else
       call allocate_op(Fres_nd, size(F1_nd%prod_op1d)+1)

       do i = 1, size(F1_nd%prod_op1d)
           Fres_nd%prod_op1d(i) = F1_nd%prod_op1d(i)
       end do

       Fres_nd%prod_op1d(size(F1_nd%prod_op1d)+1) = F2_1d

     end if
     call Simplify_OpnD(Fres_nd)

     IF (debug) THEN
       CALL write_op(Fres_nd,header=.TRUE.)
       write(out_unitp,*) 'END ',routine_name
       CALL flush_perso(out_unitp)
     END IF
   END SUBROUTINE get_F1_nd_times_F2_1d_to_Fres_nd

   !! @description: Does the product of an nd operator F1_nd with
   !!               an nd operator F2_nd
   !!               and save the result in Fres_nd.
   !! @param.in:    F1_nd the 1d operator (type: opnd).
   !! @param.in:    F2_nd_the nd operator (type: opnd).
   !! @param.out:   Fres_nd  Array of nd operators in which
   !!               the result is saved (type: opnd).
   SUBROUTINE get_F1_nd_times_F2_nd_to_Fres_nd(F1_nd, F2_nd, Fres_nd)
     type(opnd),       intent(in)       :: F1_nd
     type(opnd),       intent(in)       :: F2_nd
     type(opnd),       intent(inout)    :: Fres_nd

     type(opnd)                 :: Ftmp_nd
     integer                    :: iF1_opzero, iF2_opzero
     integer                    :: jF1_opzero, jF2_opzero
     integer                    :: i, j

     !logical, parameter           :: debug=.TRUE.
     logical, parameter           :: debug=.FALSE.
     character (len=*), parameter :: routine_name='get_F1_nd_times_F2_nd_to_Fres_nd'

     IF (debug) THEN
       write(out_unitp,*) ' BEGINNING ',routine_name
       CALL flush_perso(out_unitp)
     END IF

   CALL check_allocate_op(F1_nd)
   CALL check_allocate_op(F2_nd)


     call present_op_zero_in_F_nd(F1_nd, iF1_opzero,jF1_opzero, 'F1_nd from '//routine_name)
     if(jF1_opzero /= -1 ) then
       Fres_nd = F1_nd%prod_op1d(iF1_opzero)%prod_opel(jF1_opzero)
       IF (debug) THEN
         CALL write_op(Fres_nd,header=.true.)
         write(out_unitp,*) ' END ',routine_name
         CALL flush_perso(out_unitp)
       END IF
       return
     end if

     call present_op_zero_in_F_nd(F2_nd, iF2_opzero,jF2_opzero, 'F2_nd from '//routine_name)
     if(jF2_opzero /= -1 ) then
       Fres_nd = F2_nd%prod_op1d(iF2_opzero)%prod_opel(jF2_opzero)
       IF (debug) THEN
         CALL write_op(Fres_nd,header=.true.)
         write(out_unitp,*) ' END ',routine_name
         CALL flush_perso(out_unitp)
       END IF
       return
     end if

     do i = 1, size(F1_nd%prod_op1d)
       call check_idq_in_F_1d(F1_nd%prod_op1d(i), 'F1_nd from '//routine_name)
     end do
     do i = 1, size(F2_nd%prod_op1d)
       call check_idq_in_F_1d(F2_nd%prod_op1d(i), 'F2_nd from '//routine_name)
     end do

     Ftmp_nd = F2_nd
     do i = 1, size(F1_nd%prod_op1d)
       call get_F1_1d_times_F2_nd_to_Fres_nd(F1_nd%prod_op1d(i),Ftmp_nd, Fres_nd)
       Ftmp_nd = Fres_nd
     end do

     call delete_Op(Ftmp_nd)

     CALL Simplify_OpnD(Fres_nd)

     IF (debug) THEN
       CALL write_op(Fres_nd,header=.true.)
       write(out_unitp,*) ' END ',routine_name
       CALL flush_perso(out_unitp)
     END IF

   END SUBROUTINE get_F1_nd_times_F2_nd_to_Fres_nd



 !! @description:       Check if op_zero is a componant of F_nd
 !! @param:  F_nd       The nd operator (type: op1d).
 !! @param:  i_opzero   The index of the componant of F_nd%prod_op1d
 !!                        which is equal to op_zero.
 !!                        If F_nd doesn't contain a zero, i_opzero
 !!                        will be initialized to -1
 !! @param:  j_opzero   The index of the componant of
 !!                        F_nd%prod_op1d(i_opzero)%prod_opel
 !!                        which is equal to op_zero.
 !!                        If F_nd doesn't contain a zero, j_opzero
 !!                        will be initialized to -1.
 !! @param:    string   Message text. It just help to localize the problem.
 !! @param:   l_zero   Optional out parameter. If l_zero = .true., that
 !!                        means that F_nd cointains a op_zero component.
 subroutine present_op_zero_in_F_nd(F_nd, i_opzero, j_opzero, string, l_zero)
   type(opnd),          intent(in)       :: F_nd
   integer,             intent(inout)    :: i_opzero
   integer,             intent(inout)    :: j_opzero
   character (len = *), intent(in)       :: string
   logical,  optional,   intent(inout)   :: l_zero

   integer                    :: i
   character (len=*), parameter :: routine_name='present_op_zero_in_F_nd'

   CALL check_allocate_op(F_nd)

   i_opzero = -1
   j_opzero = -1
   if(present(l_zero)) l_zero = .false.
   do i = 1, size(F_nd%prod_op1d)
     call present_op_zero_in_F_1d(F_nd%prod_op1d(i), &
                                 & j_opzero, string)
     if(j_opzero /= -1 ) then
       i_opzero = i
       if(present(l_zero)) l_zero = .true.
       exit
     end if
   end do
 end subroutine present_op_zero_in_F_nd

 !! @description: Simplify a nd operator by removind all Id operators
 !! @param:    F_nd    The nd operator (type: opnd). It will be
 !!                            overwritten it necessary.
 !! @param:   l_Id     Optional out parameter. If l_Id = .true., that
 !!                        means that Fin_nd = Id operator. In this case
 !!                        the size of the arrays data structures of Fout_nd
 !!                        is = 1
 subroutine remove_Idop_in_F_nd(F_nd,  l_Id)

   type(opnd),           intent(inout)    :: F_nd
   logical,  optional,   intent(inout)    :: l_Id

   type(opnd)                 :: Ftmp_nd

   integer                    :: i, ii,j, p
   integer                    :: n, ndim
   integer                    :: error, n_opId
   complex(kind=Rkind)        :: coeff
   logical, allocatable       :: l_Id1(:)

   !logical, parameter           :: debug=.TRUE.
   logical, parameter           :: debug=.FALSE.
   character (len=*), parameter :: routine_name='remove_Idop_in_F_nd'


   IF (debug) THEN
     write(out_unitp,*) 'BEGINNING ',routine_name
     CALL write_op(F_nd,header=.TRUE.)
   END IF

   CALL check_allocate_op(F_nd)

   if(present(l_Id)) l_Id = .false.

   CALL alloc_NParray(l_Id1,shape(F_nd%prod_op1d),'l_Id1',routine_name)
   l_Id1(:) = .FALSE.

   n = 0
   coeff = cone
   do i = 1, size(F_nd%prod_op1d)
     call remove_Idop_in_F_1d(F_nd%prod_op1d(i), l_Id1(i))
     if(l_Id1(i)) then
       n = n+1
       coeff = coeff * F_nd%prod_op1d(i)%prod_opel(1)%coeff
     end if
   end do

   if(n == size(F_nd%prod_op1d)) then ! because all Op1D are Id op (1).
     F_nd = set_opel(idf=1, idq=1, alfa=1, indexq=1, coeff=coeff)
   else
     Ftmp_nd = F_nd
     CALL allocate_op(F_nd, size(Ftmp_nd%prod_op1d)-n)
     ii = 0
     DO i=1,size(Ftmp_nd%prod_op1d)
       IF ( .NOT. l_Id1(i) ) THEN
         ii = ii + 1
         F_nd%prod_op1d(ii) = Ftmp_nd%prod_op1d(i)
       END IF
     END DO
     F_nd%prod_op1d(1)%prod_opel(1)%coeff = F_nd%prod_op1d(1)%prod_opel(1)%coeff * coeff

     call delete_opnd(Ftmp_nd)

   end if

   CALL dealloc_NParray(l_Id1,'l_Id1',routine_name)

   IF (debug) THEN
     CALL write_op(F_nd,header=.TRUE.)
     write(out_unitp,*) 'END ',routine_name
   END IF

 end subroutine remove_Idop_in_F_nd

 !! @description:  Compare if two operator of type opnd are identical
 !!                Return true if they are identical
 !! @param:     F_1nd     The first opnd operator (type: opnd).
 !! @param:     F_2nd     The second opnd operator (type: opnd).
 LOGICAL FUNCTION compare_F1_nd_and_F2_nd(F1_nd, F2_nd)

   type(opnd),           intent(in)       :: F1_nd
   type(opnd),           intent(in)       :: F2_nd

   integer                    :: i, j
   character (len=*), parameter :: routine_name='compare_F1_nd_and_F2_nd'

   CALL check_allocate_op(F1_nd)
   CALL check_allocate_op(F2_nd)


   compare_F1_nd_and_F2_nd = (size(F1_nd%prod_op1d) == size(F2_nd%prod_op1d))
   IF (.NOT. compare_F1_nd_and_F2_nd) RETURN

   DO i = 1,size(F1_nd%prod_op1d)
     compare_F1_nd_and_F2_nd = compare_op(F1_nd%prod_op1d(i),F2_nd%prod_op1d(i))
     IF (.NOT. compare_F1_nd_and_F2_nd) EXIT
   END DO
 END FUNCTION compare_F1_nd_and_F2_nd

 !! @description: Permuts two operators of the type of op_nd
 !!               Depending on their index of Pq
 !!               Return true if the indexes of the Pq in F2_nd
 !!               are smaller than there in F1_nd
 !! @param:     F_1nd     The first opnd operator (type: opnd).
 !! @param:     F_2nd     The second opnd operator (type: opnd).
 LOGICAL FUNCTION permut_F1_nd_and_F2_nd(F1_nd, F2_nd)

   type(opnd),           intent(in)       :: F1_nd
   type(opnd),           intent(in)       :: F2_nd

   integer                    :: i, j, p
   integer                    :: k1, k2
   integer                    :: k3, k4
   integer                    :: idq1_F1, idq2_F1
   integer                    :: idq1_F2, idq2_F2
   character (len=*), parameter :: routine_name='permut_F1_nd_and_F2_nd'

   CALL check_allocate_op(F1_nd)
   CALL check_allocate_op(F2_nd)


   permut_F1_nd_and_F2_nd = .false.
   k1 = 0; k2 = 0
   k3 = 0; k4 = 0
   do i = 1,size(F1_nd%prod_op1d)
     do j = 1,size(F1_nd%prod_op1d(i)%prod_opel)
       if(F1_nd%prod_op1d(i)%prod_opel(j)%idf == 4 .and. k1 == 0) then
         k1 = F1_nd%prod_op1d(i)%prod_opel(j)%indexq
         idq1_F1 = F1_nd%prod_op1d(i)%prod_opel(j)%idq
       else if (F1_nd%prod_op1d(i)%prod_opel(j)%idf == 4 .and. k2 == 0) then
         k2 = F1_nd%prod_op1d(i)%prod_opel(j)%indexq
         idq2_F1 = F1_nd%prod_op1d(i)%prod_opel(j)%idq
       end if
     end do
   end do
   if(k1 == 0 .or. k2 == 0) then
     return
   end if
   do i = 1,size(F2_nd%prod_op1d)
     do j = 1,size(F2_nd%prod_op1d(i)%prod_opel)
       if(F2_nd%prod_op1d(i)%prod_opel(j)%idf == 4 .and. k3 == 0) then
         k3 = F2_nd%prod_op1d(i)%prod_opel(j)%indexq
         idq1_F2 = F2_nd%prod_op1d(i)%prod_opel(j)%idq
       else if (F2_nd%prod_op1d(i)%prod_opel(j)%idf == 4 .and. k4 == 0) then
         k4 = F2_nd%prod_op1d(i)%prod_opel(j)%indexq
         idq2_F2 = F2_nd%prod_op1d(i)%prod_opel(j)%idq
       end if
     end do
   end do
   if(k3 == 0 .or. k4 == 0) then
     return
   end if
   if(k1>k3  .and. k2 > k4) then
      permut_F1_nd_and_F2_nd = .true.
   else if(k1<=k3  .and. k2 > k4) then
      permut_F1_nd_and_F2_nd = .true.
   else if(k1>=k3 .and. k2<k4) then
      permut_F1_nd_and_F2_nd = .true.
   else if(k1>k3) then
      permut_F1_nd_and_F2_nd = .true.
   end if
 END FUNCTION permut_F1_nd_and_F2_nd

 SUBROUTINE Change_PQ_OF_OpnD_TO_Id_OF_OnD(F_OpnD)
   type(opnd),           intent(inout)       :: F_OpnD ! Product of Op1D

   integer                     :: i

   character (len=*), parameter   :: routine_name='Change_PQ_OF_OpnD_TO_Id_OF_OnD'

   DO i=1,size(F_OpnD%prod_op1d)

     CALL Change_PQ_OF_Op1D_TO_Id_OF_Op1D(F_OpnD%prod_op1D(i))

   END DO

   CALL Simplify_OpnD(F_OpnD)
 END SUBROUTINE Change_PQ_OF_OpnD_TO_Id_OF_OnD

 subroutine Sort_OpnD(FOpnD)
   type(opnd),       intent(inout)       :: FOpnD

   type(opnd)       :: temp_OpnD
   integer          :: i,j,n,idqi,idqj,indexQi,indexQj
   integer          :: pq(2),JJ(2),LL(2),nb_J,i_JJ,ii
   logical          :: Li
   type(op1D)       :: tempOp

   !logical, parameter           :: debug=.TRUE.
   logical, parameter           :: debug=.FALSE.
   character (len =*), parameter :: routine_name='Sort_OpnD'

   IF (debug) THEN
     write(out_unitp,*) ' BEGINNING ',routine_name
    write(out_unitp,*) ' Unsorted OpnD',size(FOpnD%prod_op1d)
     CALL write_op(FOpnD,header=.TRUE.)
     CALL flush_perso(out_unitp)
   END IF
   n = size(FOpnD%prod_op1D)

   ! first put the Jx, Jy or Jz operators at the end, because we cannot sort them (JxJy /= JyJx)
   CALL get_pqJL_OF_OpnD(pq,JJ,LL,FOpnD)
   IF (JJ(1) > 0 .AND. JJ(2) > 0)THEN
     IF (JJ(1) == JJ(2) ) THEN
       nb_J = 1
     ELSE
       nb_J = 2
     END IF
   ELSE IF (JJ(1) > 0 .AND. JJ(2) == 0)THEN
     nb_J = 1
   ELSE
     nb_J = 0
   END IF
   !write(6,*) 'n,nb_J',n,nb_J
   IF (nb_J > 0) THEN
     temp_OpnD = FOpnD

     i_JJ = n
     ii   = n-nb_J
     DO i=n,1,-1
       idqi = get_idq_OF_Op1D(temp_OpnD%prod_op1D(i))
        Li = (idqi == 5) ! Jx,Jy,Jz
        !write(6,*) 'idqi,Li',idqi,Li ; flush(6)
        !write(6,*) 'i,ii,i_JJ',i,ii,i_JJ ; flush(6)
        IF (Li) THEN
          FOpnD%prod_op1D(i_JJ) = temp_OpnD%prod_op1D(i)
          i_JJ = i_JJ-1
        ELSE
          FOpnD%prod_op1D(ii) = temp_OpnD%prod_op1D(i)
          ii = ii-1
        END IF
     END DO

     CALL delete_op(temp_OpnD)
   END IF


   ! sort the op1D
   DO i=1,n-nb_J
   DO j=i+1,n-nb_J
     indexQi = get_indexQ_OF_Op1D(FOpnD%prod_op1D(i))
     indexQj = get_indexQ_OF_Op1D(FOpnD%prod_op1D(j))

     IF ( indexQi > indexQj ) THEN
       tempOp             = FOpnD%prod_op1D(j)
       FOpnD%prod_op1D(j) = FOpnD%prod_op1D(i)
       FOpnD%prod_op1D(i) = tempOp
     END IF
   END DO
   END DO

   IF (debug) THEN
     write(out_unitp,*) ' Sorted OpnD'
     CALL write_op(FOpnD,header=.TRUE.)
     write(out_unitp,*) ' END ',routine_name
     CALL flush_perso(out_unitp)
   END IF

 END SUBROUTINE Sort_OpnD

 SUBROUTINE Simplify_OpnD(FOpnD)

   type(opnd),               intent(inout)       :: FOpnD


   integer                    :: i_opzero, j_opzero



   !logical, parameter           :: debug=.TRUE.
   logical, parameter           :: debug=.FALSE.
   character (len =*), parameter :: routine_name='Simplify_OpnD'

   IF (debug) THEN
     write(out_unitp,*) ' BEGINNING ',routine_name
    write(out_unitp,*) ' Unsimplified OpnD',size(FOpnD%prod_op1d)
     CALL write_op(FOpnD,header=.TRUE.)
     CALL flush_perso(out_unitp)
   END IF

   ! 1) check if it is zero
   CALL present_op_zero_in_F_nd(FOpnD, i_opzero, j_opzero, 'from' // routine_name)
   IF (i_opzero > 0) THEN
     FOpnD = czero
   ELSE

     CALL Sort_OpnD(FOpnD)

     CALL remove_Idop_in_F_nd(FOpnD)
   END IF


   IF (debug) THEN
     write(out_unitp,*) ' Simplified OpnD',size(FOpnD%prod_op1d)
     CALL write_op(FOpnD,header=.TRUE.)
     write(out_unitp,*) ' END ',routine_name
     CALL flush_perso(out_unitp)
   END IF

 END SUBROUTINE Simplify_OpnD



subroutine Expand_OpnD_TO_SumOpnD(FOpnD,SumOpnD)

   type(opnd), allocatable,  intent(inout)    :: SumOpnD(:) ! it will contain a sum of OpnD
   type(opnd),               intent(in)       :: FOpnD

   integer                    :: pq(2),JJ(2),LL(2),k,ip,jp
   type(opnd)                 :: temp_OpnD,temp_OpnD_ij

   integer                    :: i,j,ij,ndim_i,ndim_j
   type(Sum_OF_op1d)          :: temp_SumOp1D_i,temp_SumOp1D_j



   !logical, parameter           :: debug=.TRUE.
   logical, parameter           :: debug=.FALSE.
   character (len =*), parameter :: routine_name='Expand_OpnD_TO_SumOpnD'

   IF (debug) THEN
     write(out_unitp,*) ' BEGINNING ',routine_name
     CALL write_op(FOpnD,header=.TRUE.)
     CALL flush_perso(out_unitp)
   END IF

   ! 1) Determined in which prod_op1d(:) the P operator are located
   ip = 0
   jp = 0
   DO k=1,size(FOpnD%prod_op1d)
     CALL get_pqJL_OF_Op1D(pq,JJ,LL,FOpnD%prod_op1d(k))
     IF (pq(1) /= 0) THEN
       IF (ip == 0) THEN
         ip = k
       ELSE IF (jp == 0) THEN
         jp = k
       ELSE
         write(out_unitp,*) ' ERROR in ',routine_name
         write(out_unitp,*) ' More than TWO Pq operators'
         write(out_unitp,*) ' ip,jp',ip,jp
         STOP
       END IF
     END IF
   END DO

  !2) Differente cases in function of ip and jp values
  IF (ip == 0 .AND. jp == 0) THEN
    CALL alloc_NParray(SumOpnD,(/1/),'SumOpnD',routine_name)
    SumOpnD(1) = FOpnD
  ELSE IF (ip > 0 .AND. jp == 0) THEN
    temp_OpnD = FOpnD
    temp_OpnD%prod_op1d(ip) = cone

    IF (debug) THEN
      write(out_unitp,*) ' temp_OpnD ',ip
      CALL write_op(temp_OpnD)
      write(out_unitp,*) ' FOpnD%prod_op1d(ip) ',ip
      CALL write_op(FOpnD%prod_op1d(ip))
    END IF


    CALL Expand_Op1D_TO_SumOp1D( FOpnD%prod_op1d(ip), temp_SumOp1D_i)

    IF (debug) THEN
      write(out_unitp,*) ' temp_SumOp1D_i',ip
      CALL write_op(temp_SumOp1D_i)
    END IF

    ndim_i = size(temp_SumOp1D_i%Sum_op1D)
    CALL alloc_NParray(SumOpnD,(/ndim_i/),'SumOpnD',routine_name)
    DO i=1,ndim_i
      CALL get_F1_times_F2_to_F_nd(temp_OpnD,temp_SumOp1D_i%Sum_op1D(i),SumOpnD(i))
    END DO
  ELSE
    temp_OpnD = FOpnD
    temp_OpnD%prod_op1d(ip) = cone
    temp_OpnD%prod_op1d(jp) = cone

    CALL Expand_Op1D_TO_SumOp1D(FOpnD%prod_op1d(ip),temp_SumOp1D_i)
    CALL Expand_Op1D_TO_SumOp1D(FOpnD%prod_op1d(jp),temp_SumOp1D_j)

    ndim_i = size(temp_SumOp1D_i%Sum_op1D)
    ndim_j = size(temp_SumOp1D_j%Sum_op1D)

    CALL alloc_NParray(SumOpnD,(/ndim_i*ndim_j/),'SumOpnD',routine_name)
    ij = 0
    DO i=1,ndim_i
    DO j=1,ndim_j
      CALL get_F1_times_F2_to_F_nd(temp_SumOp1D_i%Sum_op1D(i),          &
                                   temp_SumOp1D_j%Sum_op1D(j),temp_OpnD_ij)
      ij = ij+1
      CALL get_F1_times_F2_to_F_nd(temp_OpnD,temp_OpnD_ij,SumOpnD(ij))
    END DO
    END DO

  END IF

  CALL delete_op(temp_OpnD)
  CALL delete_op(temp_SumOp1D_i)
  CALL delete_op(temp_SumOp1D_j)

  IF (debug) THEN
    write(out_unitp,*) ' Expand of OpnD',size(SumOpnD)
    DO ij=1,size(SumOpnD)
      write(out_unitp,*) ' Term in the sum: ',ij
      CALL write_op(SumOpnD(ij),header=.TRUE.)
    END DO
    write(out_unitp,*) ' END ',routine_name
    CALL flush_perso(out_unitp)
  END IF

 end subroutine Expand_OpnD_TO_SumOpnD

   FUNCTION get_coeff_OF_OpnD(FnD)
     type(OpnD),                intent(in)       :: FnD
     complex (kind=Rkind)                        :: get_coeff_OF_OpnD

     complex (kind=Rkind)                        :: temp
     integer :: i
     character (len=*), parameter   :: routine_name='get_coeff_OF_OpnD'

     temp = cone
     DO i=1,size(FnD%prod_op1d)
       temp = temp * get_coeff_OF_Op1D(FnD%prod_op1d(i))
     END DO
     get_coeff_OF_OpnD = temp

   END FUNCTION get_coeff_OF_OpnD

   SUBROUTINE Set_coeff_OF_OpnD_TO_ONE(FnD)
     type(OpnD),                intent(inout)       :: FnD

     integer :: i
     character (len=*), parameter   :: routine_name='Set_coeff_OF_OpnD_TO_ONE'

     DO i=1,size(FnD%prod_op1d)
       CALL Set_coeff_OF_Op1D_TO_ONE(FnD%prod_op1d(i))
     END DO

   END SUBROUTINE Set_coeff_OF_OpnD_TO_ONE

   SUBROUTINE get_pq_OF_OpnD(pq1,pq2,nb_pq,nb_J,nb_L,nb_var,F_nd)
     type(opnd),                intent(in)       :: F_nd
     integer,                   intent(in)       :: nb_var

     integer,                   intent(inout)    :: pq1,pq2,nb_pq,nb_J,nb_L

     integer    :: pq(2),J(2),L(2)

     integer             :: i
     character (len=*), parameter   :: routine_name='get_pq_OF_OpnD'

     CALL get_pqJL_OF_OpnD(pq,J,L,F_nd)
     nb_pq = count(pq /= 0)
     nb_J  = count(J /= 0)
     nb_L  = count(L /= 0)
     pq1   = pq(1)
     pq2   = pq(2)
     IF (nb_J > 0 .AND. nb_pq == 0) THEN
       pq1 = J(1)
       pq2 = J(2)
     ELSE IF (nb_J > 0 .AND. nb_pq == 1) THEN
       pq2 = J(1)
     END IF

   END SUBROUTINE get_pq_OF_OpnD

   SUBROUTINE get_pqJL_OF_OpnD(pq,J,L,F_nd)
     type(opnd),                intent(in)       :: F_nd
     integer,                   intent(inout)    :: pq(2),J(2),L(2)

     integer             :: pq_1D(2),J_1D(2),L_1D(2)
     integer             :: i
     character (len=*), parameter   :: routine_name='get_pqJL_OF_OpnD'

     pq   = 0
     J    = 0
     L    = 0
     DO i=1,size(F_nd%prod_op1d)

       CALL get_pqJL_OF_Op1D(pq_1D,J_1D,L_1D,F_nd%prod_op1d(i))

       CALL set_pqORJORL(pq,pq_1D)
       CALL set_pqORJORL(J,J_1D)
       CALL set_pqORJORL(L,L_1D)

     END DO

   END SUBROUTINE get_pqJL_OF_OpnD


   SUBROUTINE set_indexQ_OF_OpnD(F_nd)
     type(opnd),                intent(inout)       :: F_nd

     integer             :: i
     character (len=*), parameter   :: routine_name='set_indexQ_OF_OpnD'

     DO i=1,size(F_nd%prod_op1d)
       !CALL set_indexQ_OF_Op1D(F_nd%prod_op1d(i))
     END DO
   END SUBROUTINE set_indexQ_OF_OpnD

   !! @description: numerical calculation of an 1D Op,
   !! @param:       F_1d      The operator (type: op1d).
   !! @param:       Qval      value of the coordinate associated with F_1d
   !! @param:       ValOpEl   value of the operator, F_1d
   SUBROUTINE get_NumVal_OpnD(ValOp,Qval,F_nd)
     type(opnd),                intent(in)       :: F_nd
     real(kind=Rkind),          intent(in)       :: Qval(:)
     complex(kind=Rkind),       intent(inout)    :: ValOp


     complex(kind=Rkind)    :: ValOp1D
     integer                :: i,iQval
     character (len=*), parameter   :: routine_name='get_NumVal_OpnD'

     ValOp = CONE
     DO i =1,size(F_nd%prod_op1d)
       iQval = get_indexQ_OF_Op1D(F_nd%prod_op1d(i))
       IF (iQval == 0) THEN
         CALL get_NumVal_Op1D(ValOp1D,zero,F_nd%prod_op1d(i))
         ValOp = ValOp * ValOp1D
       ELSE IF (iQval <= size(Qval)) THEN
         CALL get_NumVal_Op1D(ValOp1D,Qval(iQval),F_nd%prod_op1d(i))
         ValOp = ValOp * ValOp1D
       END IF
     END DO
   END SUBROUTINE get_NumVal_OpnD

 end module mod_Tana_OpnD
