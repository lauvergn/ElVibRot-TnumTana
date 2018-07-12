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
!      with contributions of Mamadou Ndong
!
!===========================================================================
!===========================================================================

 module mod_Tana_OpEl
   !! @description: This module defines the data structures. It contains also
   !!               the some standard routine that initilize, 
   !!               allocate and delete the data structure
 USE mod_system
 USE mod_string
 IMPLICIT NONE
! PRIVATE

        !-----------------------------------------------------------!
        !                OpEl                                       !
        !-----------------------------------------------------------!

        !! @description:   Information about an elementary
        !!                 operator.
        !!                 This type is used for the analytical
        !!                 computation of the KEO
        !!                 like <span class="eq"> \hat{P}_q, \,\, q^\alpha,
        !!                 \,\, \cos^\alpha q, \,\, \sin^\alpha q </span>
        !! @param: idf     Integer which identifies the an elementary op.
        !! @param: idq     Integer which identifies the coordinate-dependent
        !!                 of the op. in the BF frame, (R,\theta\phi).
        !!                 Ex: idq = 1 means that the elementary op
        !!                 depends on R. If idq = 2, it depends on theta.
        !!                 if idq = 4, that would mean that Cartesian
        !!                 coordinates are used.
        !! @param: alfa      numerator of the power of the elementary op.
        !! @param: alfa_den  denominator of the power of the elementary op.
        !! @param: indexq  the index of the coordinate in the BF.
        !! @param: opname  name of to the  elementary op.
        !! @param: coeff   Coefficient in front of the operator.
        TYPE OpEl
          integer                   :: idf        = 0
          integer                   :: idq        = 0
          integer                   :: iv         = 0
          TYPE(FracInteger)         :: alfa       = FracInteger(0,1)

          integer                   :: indexq     = 0
          character(len = Name_len) :: opname     = 'init0'
          complex(kind = Rkind)     :: coeff      = CZERO
        END TYPE OpEl

      INTERFACE alloc_NParray
        MODULE PROCEDURE alloc_NParray_OF_OpEldim1
      END INTERFACE
      INTERFACE dealloc_NParray
        MODULE PROCEDURE dealloc_NParray_OF_OpEldim1
      END INTERFACE
      INTERFACE check_NParray
        MODULE PROCEDURE check_NParray_OF_OpEldim1
      END INTERFACE

  !!@description: Generic routine that compares the index of the coordinate on
  !!              which depends two 1d-operators 
  interface compare_indexq
    module procedure compare_indexq_F1el_F2el
  end interface
 
  !!@description: Generic routine that compares two 1d-operators
  interface compare_op
    module procedure compare_F1el_F2el
  end interface

  !!@description: Generic routine that copy a operator F1 to another operator F2
  interface copy_F1_into_F2
    module procedure copy_F1_el_into_F2_el
  end interface

   INTERFACE write_op
     module procedure write_opel
   END INTERFACE

   INTERFACE assignment (=)
     MODULE PROCEDURE R_TO_OpEl,C_TO_OpEl
   END INTERFACE

   INTERFACE operator (*)
     MODULE PROCEDURE R_times_OpEl,OpEl_times_R,C_times_OpEl,OpEl_times_C
   END INTERFACE

  contains

      SUBROUTINE alloc_NParray_OF_OpEldim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      type(opel), allocatable, intent(out) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=1
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_NParray_OF_OpEldim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       CALL dealloc_NParray_OF_OpEldim1(tab,name_var,name_sub_alloc // ' from ' // name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)

         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1)),stat=err_mem)
       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'OpEl')

      END SUBROUTINE alloc_NParray_OF_OpEldim1
      SUBROUTINE dealloc_NParray_OF_OpEldim1(tab,name_var,name_sub)
      IMPLICIT NONE

      type(opel), allocatable, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_NParray_OF_OpEldim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       IF (.NOT. allocated(tab)) RETURN

       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'OpEl')

      END SUBROUTINE dealloc_NParray_OF_OpEldim1

      SUBROUTINE check_NParray_OF_OpEldim1(tab,name_var,name_sub)
      IMPLICIT NONE

      type(opel), allocatable, intent(in) :: tab(:)

      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'check_NParray_OF_OpEldim1'
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

      END SUBROUTINE check_NParray_OF_OpEldim1

 function idq_TO_CoordName(idq,err_el)
 character (len=:), allocatable :: idq_TO_CoordName
 integer, intent(in) :: idq
 integer, intent(inout), optional :: err_el

 integer:: err_el_loc
 character (len = *), parameter :: routine_name = 'idq_TO_CoordName'

 err_el_loc = 0
 SELECT CASE(idq)
   CASE(1)
     idq_TO_CoordName = String_TO_String('cart')

   CASE(2)
     idq_TO_CoordName = String_TO_String('r')

   CASE(3)
     idq_TO_CoordName = String_TO_String('t')
   CASE(-3)
     idq_TO_CoordName = String_TO_String('u')

   CASE(4)
     idq_TO_CoordName = String_TO_String('f')

   CASE(6)
     idq_TO_CoordName = String_TO_String('a')

   CASE(7)
     idq_TO_CoordName = String_TO_String('b')
   CASE(-7)
     idq_TO_CoordName = String_TO_String('ub')

   CASE(8)
     idq_TO_CoordName = String_TO_String('g')

   CASE DEFAULT
     write(out_unitp,*) 'ERROR in ',routine_name
     write(out_unitp,*) 'idq=', idq
     write(out_unitp,*) "This idq does not exist"
     err_el_loc = 1
   END SELECT

   IF (present(err_el)) THEN
     err_el = err_el_loc
   ELSE
     STOP 'wrong idq'
   END IF

 end function idq_TO_CoordName

 !! @description: Initialized a list of elementary operators.
 subroutine init_list_opel()
   type(opel)               :: Fel
   integer                  :: error
   character (len = *), parameter :: routine_name = 'init_list_opel'

   Fel%idf =  0; Fel%opname = 'op_zero'
   Fel%idf =  1; Fel%opname = 'Id'
   Fel%idf =  2; Fel%opname = 'q^alfa'
   Fel%idf =  3; Fel%opname = 'qus^alfa' ! sqrt(1-q^2)^alfa

   Fel%idf =  4; Fel%opname = 'Pq'   ! -I*hb * dq

   Fel%idf =  5; Fel%opname = 'cosq^alfa'
   Fel%idf =  6; Fel%opname = 'sinq^alfa'
   Fel%idf =  7; Fel%opname = 'tanq^alfa'
   Fel%idf =  8; Fel%opname = 'cotq^alfa'
   Fel%idf =  9; Fel%opname = 'J_x'
   Fel%idf = 10; Fel%opname = 'J_y'
   Fel%idf = 11; Fel%opname = 'J_z'

   Fel%idf = 12; Fel%opname = 'Pq_q^alfa'

   Fel%idf = 13; Fel%opname = 'q^alfa_Pq'

   Fel%idf = 14; Fel%opname = 'Pq_cosq^alfa'

   Fel%idf = 15; Fel%opname = 'cosq^alfa_Pq'

   Fel%idf = 16; Fel%opname = 'Pq_sinq^alfa'

   Fel%idf = 17; Fel%opname = 'sinq^alfa_Pq'

   Fel%idf = 18; Fel%opname = 'Pq_tanq^alfa'

   Fel%idf = 19; Fel%opname = 'tanq^alfa_Pq'

   Fel%idf = 20; Fel%opname = 'Pq_cotq^alfa'

   Fel%idf = 21; Fel%opname = 'cotq^alfa_Pq'

   Fel%idf = 22; Fel%opname = 'Pu_qus^alfa'

   Fel%idf = 23; Fel%opname = 'qus^alfa_Pu'

   Fel%idf = 24; Fel%opname = 'L_x'
   Fel%idf = 25; Fel%opname = 'L_y'
   Fel%idf = 26; Fel%opname = 'L_z'

 end subroutine init_list_opel

 !! @description: Initialize the elementary operator
 !! @param: Fel     elementary operator which will be initialized
 !! @param: idf     id number of Fel
 !! @param: idq     id coordinate on which depends Fel
 !! @param: alfa    Power of the standard operator
 !! @param: indexq  index of the coordinate
 !! @param: coeff   coefficient value
 !! @param: qval    value given to the coordinate
 subroutine get_opel(Fel, idf, idq, alfa, indexq, coeff, err_el)
   type(opel),                   intent(inout)   :: Fel
   integer,                      intent(in)      :: idf
   integer,                      intent(in)      :: idq
   TYPE(FracInteger),            intent(in)      :: alfa
   integer,                      intent(in)      :: indexq
   complex(kind = Rkind),        intent(in)      :: coeff
   integer, optional,            intent(inout)   :: err_el

   TYPE(FracInteger) :: alfa_loc
   integer :: err_el_loc
   character (len = Name_len)     ::       calfa
   character (len = Name_len)     ::       cindex
   character(len = Name_len)      ::       opname

   character (len = *), parameter :: routine_name = 'get_opel'

   alfa_loc = alfa

   CALL frac_simplification(alfa_loc%num,alfa_loc%den)

   IF ( (idf > 11 .AND. idf < 24).OR. idf==7 .OR. idf==8) STOP 'get_opel idf > 11 =7,8'

   Fel%idf      = idf
   Fel%idq      = idq
   Fel%alfa     = alfa_loc
   Fel%indexq   = indexq
   Fel%coeff    = coeff

   IF ((idf == 4 .OR. idf == 9 .OR. idf == 10 .OR. idf == 11 .OR. &
                      idf ==24 .OR. idf == 25 .OR. idf == 26 ) .AND. &
        alfa_loc%den /= 1) THEN
     write(out_unitp,*) 'ERROR in ',routine_name
     write(out_unitp,*) 'idf=', idf
     write(out_unitp,*) 'alfa_loc%den=', alfa_loc%den
     write(out_unitp,*) ' The operator contains a P or a J or L and alfa_loc%den /= 1'
     err_el_loc = 1
   END IF


   !special cases:
   IF (abs(coeff) == zero) Fel%idf = 0 ! zero

   IF (alfa_loc == 1) THEN
     calfa = ''
   ELSE
     calfa = '^(' // frac_TO_string(alfa_loc) // ')'
   END IF

   IF (alfa_loc == 0 .AND. idf >= 12 .AND. idf <= 23) THEN
     Fel%idf  = 4 ! Pq
     Fel%alfa = 1
     calfa    = ''
   END IF


   err_el_loc = 0
   select case (Fel%idf)
     case(0) 
       Fel%opname   = 'op_zero'
       Fel%coeff    = zero
       Fel%alfa     = 1

     case(1) 
       Fel%opname   = 'Id'
       Fel%alfa     = 1

     case(2) 
       IF (idq /= 1 .AND. idq /= 2 .AND. idq /=  -3 .AND. idq /= -7) THEN
         write(out_unitp,*) 'ERROR in ',routine_name
         write(out_unitp,*) 'idq=', idq
         write(out_unitp,*) 'idf=', idf
         write(out_unitp,*) "for idf = 2, idq should be:"
         write(out_unitp,*) "   1 (cart.) or 2 (=r) or -3 (=ut) or -7 (=ub)"
         err_el_loc = 1
       END IF

       if(idq == 1) then
         Fel%opname = 'qcart'//calfa
        else if(idq == 2) then
          Fel%opname = 'r'//calfa
       else if(idq == -3) then
         Fel%opname = 'u'//calfa
       else if(idq == -7) then
         Fel%opname = 'ub'//calfa
       end if

     case(3) 
       if (idq == -3) then
         Fel%opname = 'qus'//calfa
       else if(idq == -7) then
         Fel%opname = 'qubs'//calfa
       else
         write(out_unitp,*) 'ERROR in ',routine_name
         write(out_unitp,*) 'idq=', idq
         write(out_unitp,*) 'idf=', idf
         write(out_unitp,*) "for idf = 3, idq should be:"
         write(out_unitp,*) "   -3 (=ut) or -7 (=ub)"
         err_el_loc = 1
       end if

     case(4) 
       if(idq == 1) then
         Fel%opname = 'Pqcart' // calfa
       else if(idq == 2) then
         Fel%opname = 'Pqr' // calfa
       else if(idq == 3) then
         Fel%opname = 'Pqt' // calfa
       else if(idq == 4) then
         Fel%opname = 'Pqf' // calfa
       else if(idq == -3) then
         Fel%opname = 'Pqu' // calfa
       else if(idq == 6) then
         Fel%opname = 'Pqa' // calfa
       else if(idq == 7) then
         Fel%opname = 'Pqb' // calfa
       else if(idq == -7) then
         Fel%opname = 'Pqub' // calfa
       else if(idq == 8) then
         Fel%opname = 'Pqg' // calfa
       else
         write(out_unitp,*) 'ERROR in ',routine_name
         write(out_unitp,*) 'idq=', idq
         write(out_unitp,*) 'idf=', idf
         write(out_unitp,*) "This idq does not exist"
         err_el_loc = 1
       end if

      case(5)
         if(idq == 3) then 
           Fel%opname = 'cosqt'//calfa
         else if(idq == 4) then
           Fel%opname = 'cosqf'//calfa
         else if(idq == 6) then
           Fel%opname = 'cosqa'//calfa
         else if(idq == 7) then
           Fel%opname = 'cosqb'//calfa
         else if(idq == 8) then
           Fel%opname = 'cosqg'//calfa
         else
           write(out_unitp,*) 'ERROR in ',routine_name
           write(out_unitp,*) 'idq=', idq
           write(out_unitp,*) 'idf=', idf
           write(out_unitp,*) "for idf = 5, idq should be:"
           write(out_unitp,*) "   3 (=\theta) or 4 (=\phi) or"
           write(out_unitp,*) "   6 (=\alpha or 7 (=\beta) or 7 (=\gamma)"
           err_el_loc = 1
         end if

     case(6) 
         if(idq == 3) then 
           Fel%opname = 'sinqt'//calfa
         else if(idq == 4) then
           Fel%opname = 'sinqf'//calfa
         else if(idq == 6) then
           Fel%opname = 'sinqa'//calfa
         else if(idq == 7) then
           Fel%opname = 'sinqb'//calfa
         else if(idq == 8) then
           Fel%opname = 'sinqg'//calfa
         else
           write(out_unitp,*) 'ERROR in ',routine_name
           write(out_unitp,*) 'idq=', idq
           write(out_unitp,*) 'idf=', idf
           write(out_unitp,*) "for idf = 6, idq should be:"
           write(out_unitp,*) "   3 (=\theta) or 4 (=\phi) or"
           write(out_unitp,*) "   6 (=\alpha or 7 (=\beta) or"
           write(out_unitp,*) "   7 (=\gamma) for idf = 6"
           err_el_loc = 1
         end if

     case(7) 
         if(idq == 3) then
           Fel%opname = 'tanqt'//calfa
         else if(idq == 4) then
           Fel%opname = 'tanqf'//calfa
         else if(idq == 6) then
           Fel%opname = 'tanqa'//calfa
         else if(idq == 7) then
           Fel%opname = 'tanqb'//calfa
         else if(idq == 8) then
           Fel%opname = 'tanqg'//calfa
         else
           write(out_unitp,*) 'ERROR in ',routine_name
           write(out_unitp,*) 'idq=', idq
           write(out_unitp,*) 'idf=', idf
           write(out_unitp,*) "for idf = 7, idq should be:"
           write(out_unitp,*) "   3 (=\theta) or 4 (=\phi) or"
           write(out_unitp,*) "   6 (=\alpha or 7 (=\beta) or"
           write(out_unitp,*) "   7 (=\gamma) for idf = 6"
           err_el_loc = 1
         end if

     case(8) 
         if(idq == 3) then 
           Fel%opname = 'cotqt'//calfa
         else if(idq == 4) then
           Fel%opname = 'cotqf'//calfa
         else if(idq == 6) then
           Fel%opname = 'cotqa'//calfa
         else if(idq == 7) then
           Fel%opname = 'cotqb'//calfa
         else if(idq == 8) then
           Fel%opname = 'cotqg'//calfa
         else
           write(out_unitp,*) 'ERROR in ',routine_name
           write(out_unitp,*) 'idq=', idq
           write(out_unitp,*) 'idf=', idf
           write(out_unitp,*) "for idf = 8, idq should be:"
           write(out_unitp,*) "   3 (=\theta) or 4 (=\phi) or"
           write(out_unitp,*) "   6 (=\alpha or 7 (=\beta) or"
           write(out_unitp,*) "   7 (=\gamma) for idf = 6"
           err_el_loc = 1
         end if

     case(9) 
       Fel%idq    = 5
       Fel%opname = 'J_x'//calfa
     case(10) 
       Fel%idq    = 5
       Fel%opname = 'J_y'//calfa
     case(11) 
       Fel%idq    = 5
       Fel%opname = 'J_z'//calfa

     case(12) 
       if(idq == 1) then 
         Fel%opname = 'Pqcart_qcart'//calfa
        else if(idq == 2) then
         Fel%opname = 'Pqr_r'//calfa
       else if(idq == -3) then
         Fel%opname = 'Pqu_u'//calfa
       else
         write(out_unitp,*) 'ERROR in ',routine_name
         write(out_unitp,*) 'idq=', idq
         write(out_unitp,*) 'idf=', idf
         write(out_unitp,*) "for idf = 12, idq should be:"
         write(out_unitp,*) "   1 (cart.) or 2 (=r) or -3 (=u)"
         err_el_loc = 1
         !STOP
       end if

     case(13) 
       if(idq == 1) then 
         Fel%opname = 'qcart'//calfa//'_Pqcart'
        else if(idq == 2) then
         Fel%opname = 'r'//calfa//'_Pqr'
       else if(idq == -3) then
         Fel%opname = 'u'//calfa//'_Pqu'
       else
         write(out_unitp,*) 'ERROR in ',routine_name
         write(out_unitp,*) 'idq=', idq
         write(out_unitp,*) 'idf=', idf
         write(out_unitp,*) "for idf = 13, idq should be:"
         write(out_unitp,*) "   1 (cart.) or 2 (=r) or -3 (=u)"
         err_el_loc = 1
         !STOP
       end if

     case(14) 
       if(idq == 3) then 
         Fel%opname = 'Pqt_cosqt'//calfa
       else if(idq == 4) then
         Fel%opname = 'Pqf_cosqf'//calfa
       else
         write(out_unitp,*) 'ERROR in ',routine_name
         write(out_unitp,*) 'idq=', idq
         write(out_unitp,*) 'idf=', idf
         write(out_unitp,*) "for idf = 14, idq should be:"
         write(out_unitp,*) "   3 (=\theta) or 4 (=\phi)"
         err_el_loc = 1
         !STOP
       end if

     case(15) 
       if(idq == 3) then 
         Fel%opname = 'cosqt'//calfa//'_Pqt'
       else if(idq == 4) then
         Fel%opname = 'cosqf'//calfa//'_Pqf'
       else
         write(out_unitp,*) 'ERROR in ',routine_name
         write(out_unitp,*) 'idq=', idq
         write(out_unitp,*) 'idf=', idf
         write(out_unitp,*) "for idf = 15, idq should be:"
         write(out_unitp,*) "   3 (=\theta) or 4 (=\phi)"
         err_el_loc = 1
         !STOP
       end if

     case(16) 
       if(idq == 3) then 
         Fel%opname = 'Pqt_sinqt'//calfa
       else if(idq == 4) then
         Fel%opname = 'Pqf_sinqf'//calfa
       else
         write(out_unitp,*) 'ERROR in ',routine_name
         write(out_unitp,*) 'idq=', idq
         write(out_unitp,*) 'idf=', idf
         write(out_unitp,*) "for idf = 16, idq should be:"
         write(out_unitp,*) "   3 (=\theta) or 4 (=\phi)"
         err_el_loc = 1
         !STOP
       end if

     case(17) 
       if(idq == 3) then 
         Fel%opname = 'sinqt'//calfa//'_Pqt'
       else if(idq == 4) then
         Fel%opname = 'sinqf'//calfa//'_Pqf'
       else
         write(out_unitp,*) 'ERROR in ',routine_name
         write(out_unitp,*) 'idq=', idq
         write(out_unitp,*) 'idf=', idf
         write(out_unitp,*) "for idf = 17, idq should be:"
         write(out_unitp,*) "   3 (=\theta) or 4 (=\phi)"
         err_el_loc = 1
         !STOP
       end if

     case(18) 
       if(idq == 3) then 
         Fel%opname = 'Pqt_tanqt'//calfa
       else if(idq == 4) then
         Fel%opname = 'Pqf_tanqf'//calfa
       else
         write(out_unitp,*) 'ERROR in ',routine_name
         write(out_unitp,*) 'idq=', idq
         write(out_unitp,*) 'idf=', idf
         write(out_unitp,*) "for idf = 18, idq should be:"
         write(out_unitp,*) "   3 (=\theta) or 4 (=\phi)"
         err_el_loc = 1
         !STOP
       end if

     case(19) 
       if(idq == 3) then 
         Fel%opname = 'tanqt'//calfa//'_Pqt'
       else if(idq == 4) then
         Fel%opname = 'tanqf'//calfa//'_Pqf'
       else
         write(out_unitp,*) 'ERROR in ',routine_name
         write(out_unitp,*) 'idq=', idq
         write(out_unitp,*) 'idf=', idf
         write(out_unitp,*) "for idf = 19, idq should be:"
         write(out_unitp,*) "   3 (=\theta) or 4 (=\phi)"
         err_el_loc = 1
         !STOP
       end if

     case(20) 
       if(idq == 3) then 
         Fel%opname = 'Pqt_cotqt'//calfa
       else if(idq == 4) then
         Fel%opname = 'Pqf_cotqf'//calfa
       else
         write(out_unitp,*) 'ERROR in ',routine_name
         write(out_unitp,*) 'idq=', idq
         write(out_unitp,*) 'idf=', idf
         write(out_unitp,*) "for idf = 20, idq should be:"
         write(out_unitp,*) "   3 (=\theta) or 4 (=\phi)"
         err_el_loc = 1
         !STOP
       end if

     case(21) 
       if(idq == 3) then 
         Fel%opname = 'cotqt'//calfa//'_Pqt'
       else if(idq == 4) then
         Fel%opname = 'cotqf'//calfa//'_Pqf'
       else
         write(out_unitp,*) 'ERROR in ',routine_name
         write(out_unitp,*) 'idq=', idq
         write(out_unitp,*) 'idf=', idf
         write(out_unitp,*) "for idf = 21, idq should be:"
         write(out_unitp,*) "   3 (=\theta) or 4 (=\phi)"
         err_el_loc = 1
         !STOP
       end if

     case(22) 
       if(idq /= -3) then
         write(out_unitp,*) 'ERROR in ',routine_name
         write(out_unitp,*) 'idq=', idq
         write(out_unitp,*) 'idf=', idf
         write(out_unitp,*) "for idf = 22, idq should be: -3 (=u)"
         err_el_loc = 1
         !STOP
       else
         Fel%opname = 'Pqu_qus'//calfa
       end if 

     case(23) 
       if(idq /= -3) then
         write(out_unitp,*) 'ERROR in ',routine_name
         write(out_unitp,*) 'idq=', idq
         write(out_unitp,*) 'idf=', idf
         write(out_unitp,*) "for idf = 23, idq should be: -3 (=u)"
         err_el_loc = 1

       else
         Fel%opname = 'qus'//calfa//'_Pqu'
       end if 

     case(24) 
       Fel%idq    = -5
       Fel%opname = 'L_x'//calfa
     case(25) 
       Fel%idq    = -5
       Fel%opname = 'L_y'//calfa
     case(26) 
       Fel%idq    = -5
       Fel%opname = 'L_z'//calfa

     case default
       write(out_unitp,*) 'ERROR in ',routine_name
       write(out_unitp,*) 'idf=', idf
       write(out_unitp,*) "This idf is not registered for an elementary operator"
       write(out_unitp,*) "   illegal value of idf . (Internal Bug)"
       err_el_loc = 1
     end select

     IF (present(err_el)) THEN
       err_el = err_el_loc
     ELSE
       IF (err_el_loc /= 0) STOP
     END IF

 end subroutine get_opel

 function set_opel(idf, idq, alfa, indexq, coeff, alfa_den, err_el)
   type(opel)                                    :: set_opel
   integer,                      intent(in)      :: idf
   integer,                      intent(in)      :: idq
   integer,                      intent(in)      :: alfa
   integer,                      intent(in)      :: indexq
   complex(kind = Rkind),        intent(in)      :: coeff
   integer, optional,            intent(in)      :: alfa_den
   integer, optional,            intent(inout)   :: err_el

   TYPE(FracInteger) :: alfa_loc


   IF (present(alfa_den)) THEN
     alfa_loc = FracInteger(alfa,alfa_den)
   ELSE
     alfa_loc = FracInteger(alfa,1)
   END IF

   IF (present(err_el)) THEN
     CALL get_opel(set_opel, idf, idq, alfa_loc, indexq, coeff, err_el)
   ELSE
     CALL get_opel(set_opel, idf, idq, alfa_loc, indexq, coeff)
   END IF

 end function set_opel

 subroutine Export_MCTDH_OpEl(Fel,FelName)
   type(opel),          intent(inout)                   :: Fel

   character (len = :), allocatable, intent(inout) :: FelName

   character (len = :), allocatable     :: PName,SQName,QName

   character (len = *), parameter :: routine_name = 'Export_MCTDH_OpEl'

   IF (allocated(FelName)) deallocate(FelName)


   !CALL write_op(Fel)

   QName     = String_TO_String('q')
   PName     = String_TO_String('dq')
   SQName    = String_TO_String('qs')

   select case (Fel%idf)
     case(0) ! 0
       FelName = String_TO_String('0')
     case(1) ! Id
       FelName = String_TO_String('1')

     case(2) ! q^alfa
       FelName  = Qnamealfa_MCTDH(Qname,Fel%alfa)

     case(3) ! sqrt(1-Q^2)^alfa
       !FelName = String_TO_String('\sqrt{1-' // FuncQName // '^2}' )
       !IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // calfa)

       FelName  = Qnamealfa_MCTDH(SQname,Fel%alfa)

     case(4) ! PQ^alfa
       FelName  = Qnamealfa_MCTDH(PName,Fel%alfa)
       Fel%coeff = (-EYE)**Fel%alfa%num * Fel%coeff ! Fel%alfa%den MUST be equal to 1

     case(5) ! cos(Q)^alfa
       FelName  = Qnamealfa_MCTDH('cos',Fel%alfa)

     case(6) ! sin(Q)^alfa
       FelName = Qnamealfa_MCTDH('sin',Fel%alfa)

     case(7) ! tan(Q)^alfa
       FelName = Qnamealfa_MCTDH('tan',Fel%alfa)

     case(8) ! cot(Q)^alfa
       FelName = Qnamealfa_MCTDH('cot',Fel%alfa)

     case(9) ! Jx
       IF (Fel%alfa == 1) THEN
         FelName = String_TO_String('Jx')
       ELSE ! no more than alfa=2
         FelName = String_TO_String('Jx*Jx')
       END IF

     case(10) ! Jy
       IF (Fel%alfa == 1) THEN
         FelName = String_TO_String('Jy')
         Fel%coeff = -EYE*Fel%coeff

       ELSE ! no more than alfa=2
         FelName = String_TO_String('Jy*Jy')
         Fel%coeff = -Fel%coeff
       END IF
     case(11) ! Jz  (because Jz .psi = k Psi
       FelName  = Qnamealfa_MCTDH(Qname,Fel%alfa)

     case default
           write(out_unitp,*) 'ERROR in ',routine_name
           write(out_unitp,*) 'idf=', Fel%idf
           write(out_unitp,*) "This idf is not registered for an elementary operator"
           write(out_unitp,*) "   illegal value of idf . (Internal Bug)"
           STOP
     end select

     IF (allocated(QName))     deallocate(QName)
     IF (allocated(PName))     deallocate(PName)
     IF (allocated(SQName))    deallocate(SQName)

 end subroutine Export_MCTDH_OpEl
FUNCTION Qnamealfa_MCTDH(Qname,alfa)
  character(len=:), allocatable     :: Qnamealfa_MCTDH
  character(len=*), intent(in)      :: Qname
  TYPE(FracInteger), intent(in)     :: alfa
  real(kind=Rkind)                  :: ralfa


  IF (allocated(Qnamealfa_MCTDH)) deallocate(Qnamealfa_MCTDH)


  IF (alfa == 1) THEN
    Qnamealfa_MCTDH = String_TO_String(Qname)
  ELSE IF (alfa%den == 1) THEN
    Qnamealfa_MCTDH = String_TO_String(Qname // "^" // int_TO_char(alfa%num) )
  ELSE
    ralfa = Frac_TO_real(alfa)
    Qnamealfa_MCTDH = String_TO_String(Qname // "^" // real_TO_char(ralfa) )
  END IF

END FUNCTION Qnamealfa_MCTDH

 subroutine Export_Latex_OpEl(Fel,qname,FelName)
   type(opel),          intent(in)                   :: Fel
   character (len =*), intent(in)                    :: qname

   character (len = :), allocatable, intent(inout) :: FelName

   character (len = :), allocatable     :: PName,SQName

   character (len = *), parameter :: routine_name = 'Export_Latex_OpEl'

   IF (allocated(FelName)) deallocate(FelName)


   FelName = String_TO_String('')
   !CALL write_op(Fel)

   PName     = String_TO_String('\hat{P}_{' // qname // '}')
   SQName    = String_TO_String('v' // qname(2:len_trim(qname)) )

   select case (Fel%idf)
     case(0) ! 0
       FelName = String_TO_String('0')
     case(1) ! Id
       FelName = String_TO_String('1')

     case(2) ! q^alfa
       FelName  = Qnamealfa_Latex(Qname,Fel%alfa)

     case(3) ! sqrt(1-Q^2)^alfa
       !FelName = String_TO_String('\sqrt{1-' // FuncQName // '^2}' )
       !IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // calfa)

       FelName  = Qnamealfa_Latex(SQname,Fel%alfa)

     case(4) ! PQ^alfa
       FelName  = Qnamealfa_Latex(PName,Fel%alfa)

     case(5) ! cos(Q)^alfa
       FelName  = fnamealfa_Latex('\cos',Qname,Fel%alfa)

     case(6) ! sin(Q)^alfa
       FelName = fnamealfa_Latex('\sin',Qname,Fel%alfa)

     case(7) ! tan(Q)^alfa
       FelName = fnamealfa_Latex('\tan',Qname,Fel%alfa)

     case(8) ! cot(Q)^alfa
       FelName = fnamealfa_Latex('\cot',Qname,Fel%alfa)

     case(9) ! Jx
       FelName = Qnamealfa_Latex('\hat{J}_x',Fel%alfa)

       !FelName = String_TO_String('\hat{J}_x')
       !IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // calfa)
     case(10) ! Jy
       FelName = Qnamealfa_Latex('\hat{J}_y',Fel%alfa)

     case(11) ! Jz
       FelName = Qnamealfa_Latex('\hat{J}_z',Fel%alfa)

     case(12) ! PQ Q^alfa
       FelName = Qnamealfa_Latex(Qname,Fel%alfa)
       FelName = String_TO_String(PName // FelName)

     case(13) ! Q^alfa PQ
       FelName = Qnamealfa_Latex(Qname,Fel%alfa)
       FelName = String_TO_String(FelName // PName)

     case(14) ! PQ cos(Q)^alfa
       FelName = fnamealfa_Latex('\cos',Qname,Fel%alfa)
       FelName = String_TO_String(PName // FelName)

     case(15) ! cos(Q)^alfa PQ
       FelName = fnamealfa_Latex('\cos',Qname,Fel%alfa)
       FelName = String_TO_String(FelName // PName)

     case(16)  ! PQ sin(Q)^alfa
       FelName = fnamealfa_Latex('\sin',Qname,Fel%alfa)
       FelName = String_TO_String(PName // FelName)

     case(17) ! sin(Q)^alfa PQ
       FelName = fnamealfa_Latex('\sin',Qname,Fel%alfa)
       FelName = String_TO_String(FelName // PName)

     case(18)   ! PQ tan(Q)^alfa
       FelName = fnamealfa_Latex('\tan',Qname,Fel%alfa)
       FelName = String_TO_String(PName // FelName)


     case(19) ! tan(Q)^alfa PQ
       FelName  = fnamealfa_Latex('\tan',Qname,Fel%alfa)
       FelName = String_TO_String(FelName // PName)

     case(20)   ! PQ cot(Q)^alfa
      FelName  = fnamealfa_Latex('\cot',Qname,Fel%alfa)
       FelName = String_TO_String(PName // FelName)

     case(21)  ! cot(Q)^alfa PQ
       FelName  = fnamealfa_Latex('\cot',Qname,Fel%alfa)
       FelName = String_TO_String(FelName // PName)

     case(22)   ! PQ sqrt(1-Q^2)^alfa = PQ sQ^alfa
       !FelName = String_TO_String(PName // '\sqrt{1-' // FuncQName // '^2}' )
       !IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // calfa)

       FelName = Qnamealfa_Latex(SQname,Fel%alfa)
       FelName = String_TO_String(PName // FelName)

     case(23) ! sqrt(1-Q^2)^alfa PQ
       !FelName = String_TO_String('\sqrt{1-' // FuncQName // '^2}' )
       !IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // calfa)
       !FelName = String_TO_String(FelName // PName)

       FelName = Qnamealfa_Latex(SQname,Fel%alfa)
       FelName = String_TO_String(FelName // PName)

     case(24) ! Lx
       FelName = String_TO_String('\hat{L}_x')
     case(25) ! Ly
       FelName = String_TO_String('\hat{L}_y')
     case(26) ! Lz
       FelName = String_TO_String('\hat{L}_z')

     case default
           write(out_unitp,*) 'ERROR in ',routine_name
           write(out_unitp,*) 'idf=', Fel%idf
           write(out_unitp,*) "This idf is not registered for an elementary operator"
           write(out_unitp,*) "   illegal value of idf . (Internal Bug)"
           STOP
     end select

     IF (allocated(PName))     deallocate(PName)
     IF (allocated(SQName))    deallocate(SQName)

 end subroutine Export_Latex_OpEl
FUNCTION fnamealfa_Latex(fname,Qname,alfa)
  character(len=:), allocatable     :: fnamealfa_Latex
  character(len=*), intent(in)      :: fname,Qname
  TYPE(FracInteger), intent(in)     :: alfa

  character(len=:), allocatable     :: fnamealfa_loc
  TYPE(FracInteger)                 :: alfa_loc


  IF (allocated(fnamealfa_loc)) deallocate(fnamealfa_loc)

  IF (alfa > 0) THEN
    fnamealfa_loc = String_TO_String(fname)
    alfa_loc = alfa
  ELSE
    fnamealfa_loc = String_TO_String( "\frac{1}{" // fname)
    alfa_loc = FracInteger(-alfa%num,alfa%den)
  END IF

  IF (alfa /= 0) THEN
    fnamealfa_loc = String_TO_String(fnamealfa_loc // "^{"// frac_TO_string(alfa_loc) // "}")
  END IF

  fnamealfa_loc = String_TO_String(fnamealfa_loc // '\left(' // QName // '\right)')

  IF (alfa < 0) THEN
    fnamealfa_loc = String_TO_String(fnamealfa_loc // "}")
  END IF

  fnamealfa_Latex = fnamealfa_loc

  deallocate(fnamealfa_loc)

END FUNCTION fnamealfa_Latex

FUNCTION Qnamealfa_Latex(Qname,alfa)
  character(len=:), allocatable     :: Qnamealfa_Latex
  character(len=*), intent(in)      :: Qname
  TYPE(FracInteger),intent(in)      :: alfa
  character(len=:), allocatable     :: Qnamealfa_loc

  TYPE(FracInteger)                 :: alfa_loc

  IF (allocated(Qnamealfa_Latex)) deallocate(Qnamealfa_Latex)

  IF (alfa > 0) THEN
    Qnamealfa_loc = String_TO_String("{" // Qname // "}")
    alfa_loc = alfa
  ELSE
    Qnamealfa_loc = String_TO_String( "\frac{1}{" // "{" // Qname // "}")
    alfa_loc = FracInteger(-alfa%num,alfa%den)
  END IF

  IF (alfa /= 0) THEN
    Qnamealfa_loc = String_TO_String(Qnamealfa_loc // "^{"// frac_TO_string(alfa_loc) // "}")
  END IF

  IF (alfa < 0) THEN
    Qnamealfa_loc = String_TO_String(Qnamealfa_loc // "}")
  END IF

  Qnamealfa_Latex = Qnamealfa_loc

  deallocate(Qnamealfa_loc)

END FUNCTION Qnamealfa_Latex
 subroutine Export_Midas_OpEl(Fel, Qname, FelName)
   type(opel),          intent(in)                    :: Fel
   character (len =*),  intent(in)                    :: Qname

   character (len = :), allocatable, intent(inout)    :: FelName

   character (len = :), allocatable                   :: PName, SQName, FuncQName
   character (len = Name_len)                         :: calfa

   character (len = *), parameter                     :: routine_name = 'Export_Midas_OpEl'

   IF (allocated(FelName)) deallocate(FelName)

   FelName = String_TO_String('')
   !CALL write_op(Fel)

   calfa = '^(' // frac_TO_string(Fel%alfa) // ')'

   PName        = String_TO_String('(DDQ)')
   FuncQName    = String_TO_String(Qname)
   SQName       = String_TO_String('v' // Qname(2:len_trim(Qname)) )


!    allocate(character(len=len_trim('(DDQ)')) :: PName)
!    PName     = '(DDQ)'
!    allocate(character(len=len_trim(Qname)) :: FuncQName)
!    FuncQName = trim(Qname)
!    allocate(character(len=('v' // Qname(2:len_trim(Qname)))) :: SQName)
!    SQName    = 'v' // Qname(2:len_trim(Qname))


   select case (Fel%idf)
     case(0) ! 0
       FelName = String_TO_String('0')
     case(1) ! Id
       FelName = String_TO_String('1')

     case(2) ! q^alfa
       IF (Fel%alfa /= 1) THEN
         FelName = String_TO_String(FuncQName // trim(calfa))
       ELSE
         FelName = String_TO_String(QName)
       END IF

     case(3) ! sqrt(1-Q^2)^alfa
       !FelName = String_TO_String('sqrt(1-' // FuncQName // '^2)' )
       !IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // trim(calfa))

       IF (Fel%alfa /= 1) THEN
         FelName = String_TO_String('(' // SQName // ')')
         FelName = String_TO_String(FelName // trim(calfa))
       ELSE
         FelName = String_TO_String(SQName)
       END IF


     case(4) ! PQ^alfa
       FelName = String_TO_String(PName)
       FelName = String_TO_String(FelName // trim(calfa))

     case(5) ! cos(Q)^alfa
       FelName = String_TO_String('cos')
       IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // trim(calfa))
       FelName = String_TO_String(FelName // '(' // FuncQName // ')')

     case(6) ! sin(Q)^alfa
       FelName = String_TO_String('sin')
       IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // trim(calfa))
       FelName = String_TO_String(FelName // '(' // FuncQName // ')')

     case(7) ! tan(Q)^alfa
       FelName = String_TO_String('tan')
       IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // trim(calfa))
       FelName = String_TO_String(FelName // '(' // FuncQName // ')')

     case(8) ! cot(Q)^alfa
       FelName = String_TO_String('cot')
       IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // trim(calfa))
       FelName = String_TO_String(FelName // '(' // FuncQName // ')')

     case(9) ! Jx
       FelName = String_TO_String('Jx')
       IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // trim(calfa))
     case(10) ! Jy
       FelName = String_TO_String('Jy')
       IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // trim(calfa))
     case(11) ! Jz
       FelName = String_TO_String('Jz')
        IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // trim(calfa))

     case(12) ! PQ Q^alfa
       IF (Fel%alfa /= 1) THEN
         FelName = String_TO_String(PName // FuncQName // trim(calfa))
       ELSE
         FelName = String_TO_String(PName // QName )
       END IF

     case(13) ! Q^alfa PQ
       IF (Fel%alfa /= 1) THEN
         FelName = String_TO_String(FuncQName // trim(calfa) // PName)
       ELSE
         FelName = String_TO_String(QName // PName)
       END IF

     case(14) ! PQ cos(Q)^alfa
       FelName = String_TO_String(PName // 'cos')
       IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // trim(calfa))
       FelName = String_TO_String(FelName // FuncQName)

     case(15) ! cos(Q)^alfa PQ
       FelName = String_TO_String('cos')
       IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // trim(calfa))
       FelName = String_TO_String(FelName // FuncQName // PName)

     case(16)  ! PQ sin(Q)^alfa
       FelName = String_TO_String(PName // 'sin')
       IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // trim(calfa))
       FelName = String_TO_String(FelName // FuncQName)

     case(17) ! sin(Q)^alfa PQ
       FelName = String_TO_String('sin')
       IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // trim(calfa))
       FelName = String_TO_String(FelName // FuncQName // PName)

     case(18)   ! PQ tan(Q)^alfa
       FelName = String_TO_String(PName // 'tan')
       IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // trim(calfa))
       FelName = String_TO_String(FelName // FuncQName)

     case(19) ! tan(Q)^alfa PQ
       FelName = String_TO_String('tan')
       IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // trim(calfa))
       FelName = String_TO_String(FelName // FuncQName // PName)

     case(20)   ! PQ cot(Q)^alfa
       FelName = String_TO_String(PName // 'cot')
       IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // trim(calfa))
       FelName = String_TO_String(FelName // FuncQName)

     case(21)  ! cot(Q)^alfa PQ
       FelName = String_TO_String('cot')
       IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // trim(calfa))
       FelName = String_TO_String(FelName // FuncQName // PName)

     case(22)   ! PQ sqrt(1-Q^2)^alfa
       !FelName = String_TO_String(PName // 'sqrt(1-' // FuncQName // '^2)' )
       !IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // trim(calfa))

       FelName = String_TO_String(PName)
       IF (Fel%alfa /= 1) THEN
         FelName = String_TO_String(FelName // '(' // SQName // ')')
         FelName = String_TO_String(FelName // trim(calfa))
       ELSE
         FelName = String_TO_String(FelName // SQName)
       END IF

     case(23) ! sqrt(1-Q^2)^alfa PQ
       !FelName = String_TO_String('sqrt(1-' // FuncQName // '^2)' )
       !IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // trim(calfa))
       !FelName = String_TO_String(FelName // PName)

       IF (Fel%alfa /= 1) THEN
         FelName = String_TO_String('(' // SQName // ')')
         FelName = String_TO_String(FelName // trim(calfa))
       ELSE
         FelName = String_TO_String(SQName)
       END IF
       FelName = String_TO_String(FelName // PName)


     case default
           write(out_unitp,*) 'ERROR in ',routine_name
           write(out_unitp,*) 'idf=', Fel%idf
           write(out_unitp,*) "This idf is not registered for an elementary operator"
           write(out_unitp,*) "   illegal value of idf . (Internal Bug)"
           STOP
     end select


     IF (allocated(PName))     deallocate(PName)
     IF (allocated(SQName))    deallocate(SQName)
     IF (allocated(FuncQName)) deallocate(FuncQName)


 end subroutine Export_Midas_OpEl
 subroutine Export_VSCF_OpEl(Fel,qname,FelName)
   type(opel),          intent(in)                   :: Fel
   character (len =*), intent(in)                    :: qname

   character (len = :), allocatable, intent(inout) :: FelName

   character (len = :), allocatable     :: PName,SQName,FuncQName
   character (len = :), allocatable     :: calfa

   character (len = *), parameter :: routine_name = 'Export_VSCF_OpEl'

   IF (allocated(FelName)) deallocate(FelName)


   FelName = String_TO_String('')
   !CALL write_op(Fel)

   calfa        = String_TO_String( '^(' // frac_TO_string(Fel%alfa) // ')' )
   PName        = String_TO_String('P_' // qname)
   FuncQName    = String_TO_String('(' // qname // ')' )
   SQName       = String_TO_String('v' // qname(2:len_trim(qname)) )

   select case (Fel%idf)
     case(0) ! 0
       FelName = String_TO_String('0')
     case(1) ! Id
       FelName = String_TO_String('1')

     case(2) ! q^alfa
       IF (Fel%alfa /= 1) THEN
         FelName = String_TO_String(FuncQName // calfa)
       ELSE
         FelName = String_TO_String(QName)
       END IF

     case(3) ! sqrt(1-Q^2)^alfa
       !FelName = String_TO_String('sqrt(1-' // FuncQName // '^2)' )
       !IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // calfa)

       IF (Fel%alfa /= 1) THEN
         FelName = String_TO_String('(' // SQName // ')')
         FelName = String_TO_String(FelName // calfa)
       ELSE
         FelName = String_TO_String(SQName)
       END IF


     case(4) ! PQ^alfa
       FelName = String_TO_String(PName)
       IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // calfa)

     case(5) ! cos(Q)^alfa
       FelName = String_TO_String('cos')
       IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // calfa)
       FelName = String_TO_String(FelName // FuncQName)

     case(6) ! sin(Q)^alfa
       FelName = String_TO_String('sin')
       IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // calfa)
       FelName = String_TO_String(FelName // FuncQName)

     case(7) ! tan(Q)^alfa
       FelName = String_TO_String('tan')
       IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // calfa)
       FelName = String_TO_String(FelName // FuncQName)

     case(8) ! cot(Q)^alfa
       FelName = String_TO_String('cot')
       IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // calfa)
       FelName = String_TO_String(FelName // FuncQName)

     case(9) ! Jx
       FelName = String_TO_String('Jx')
       IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // calfa)
     case(10) ! Jy
       FelName = String_TO_String('Jy')
       IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // calfa)
     case(11) ! Jz
       FelName = String_TO_String('Jz')
        IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // calfa)

     case(12) ! PQ Q^alfa
       IF (Fel%alfa /= 1) THEN
         FelName = String_TO_String(PName // FuncQName // calfa)
       ELSE
         FelName = String_TO_String(PName // QName )
       END IF

     case(13) ! Q^alfa PQ
       IF (Fel%alfa /= 1) THEN
         FelName = String_TO_String(FuncQName // calfa // PName)
       ELSE
         FelName = String_TO_String(QName // PName)
       END IF

     case(14) ! PQ cos(Q)^alfa
       FelName = String_TO_String(PName // 'cos')
       IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // calfa)
       FelName = String_TO_String(FelName // FuncQName)

     case(15) ! cos(Q)^alfa PQ
       FelName = String_TO_String('cos')
       IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // calfa)
       FelName = String_TO_String(FelName // FuncQName // PName)

     case(16)  ! PQ sin(Q)^alfa
       FelName = String_TO_String(PName // 'sin')
       IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // calfa)
       FelName = String_TO_String(FelName // FuncQName)

     case(17) ! sin(Q)^alfa PQ
       FelName = String_TO_String('sin')
       IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // calfa)
       FelName = String_TO_String(FelName // FuncQName // PName)

     case(18)   ! PQ tan(Q)^alfa
       FelName = String_TO_String(PName // 'tan')
       IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // calfa)
       FelName = String_TO_String(FelName // FuncQName)

     case(19) ! tan(Q)^alfa PQ
       FelName = String_TO_String('tan')
       IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // calfa)
       FelName = String_TO_String(FelName // FuncQName // PName)

     case(20)   ! PQ cot(Q)^alfa
       FelName = String_TO_String(PName // 'cot')
       IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // calfa)
       FelName = String_TO_String(FelName // FuncQName)

     case(21)  ! cot(Q)^alfa PQ
       FelName = String_TO_String('cot')
       IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // calfa)
       FelName = String_TO_String(FelName // FuncQName // PName)

     case(22)   ! PQ sqrt(1-Q^2)^alfa
       !FelName = String_TO_String(PName // 'sqrt(1-' // FuncQName // '^2)' )
       !IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // calfa)

       FelName = String_TO_String(PName)
       IF (Fel%alfa /= 1) THEN
         FelName = String_TO_String(FelName // '(' // SQName // ')')
         FelName = String_TO_String(FelName // calfa)
       ELSE
         FelName = String_TO_String(FelName // SQName)
       END IF

     case(23) ! sqrt(1-Q^2)^alfa PQ
       !FelName = String_TO_String('sqrt(1-' // FuncQName // '^2)' )
       !IF (Fel%alfa /= 1) FelName = String_TO_String(FelName // calfa)
       !FelName = String_TO_String(FelName // PName)

       IF (Fel%alfa /= 1) THEN
         FelName = String_TO_String('(' // SQName // ')')
         FelName = String_TO_String(FelName // calfa)
       ELSE
         FelName = String_TO_String(SQName)
       END IF
       FelName = String_TO_String(FelName // PName)

     case(24) ! Lx
       FelName = String_TO_String('Lx')
     case(25) ! Ly
       FelName = String_TO_String('Ly')
     case(26) ! Lz
       FelName = String_TO_String('Lz')

     case default
           write(out_unitp,*) 'ERROR in ',routine_name
           write(out_unitp,*) 'idf=', Fel%idf
           write(out_unitp,*) "This idf is not registered for an elementary operator"
           write(out_unitp,*) "   illegal value of idf . (Internal Bug)"
           STOP
     end select

     IF (allocated(calfa))     deallocate(calfa)
     IF (allocated(PName))     deallocate(PName)
     IF (allocated(SQName))    deallocate(SQName)
     IF (allocated(FuncQName)) deallocate(FuncQName)


 end subroutine Export_VSCF_OpEl
 !! @description: Compares the idq of two elementary operators 
 !! @param:    F1_el     The first elementary operator. 
 !! @param:    F2_el     The second elementary operator. 
 logical FUNCTION compare_indexq_F1el_F2el(F1_el, F2_el)

   type(opel),           intent(in)       :: F1_el
   type(opel),           intent(in)       :: F2_el

   character (len=*), parameter :: routine_name='compare_indexq_F1el_F2el'

   compare_indexq_F1el_F2el = (F1_el%indexq == F2_el%indexq) .OR.       &
                            (F2_el%idf == 1) .OR. (F1_el%idf == 1) .OR. &
                            (F2_el%idf == 0) .OR. (F1_el%idf == 0)

 END FUNCTION compare_indexq_F1el_F2el

 logical FUNCTION compare_F1el_F2el(F1_el, F2_el)

   type(opel),           intent(in)       :: F1_el
   type(opel),           intent(in)       :: F2_el

   character (len=*), parameter :: routine_name='compare_F1el_F2el'

   compare_F1el_F2el = compare_indexq_F1el_F2el(F1_el,F2_el)  .AND.     &
                                    (F2_el%idf  == F1_el%idf) .AND.     &
                                    (F2_el%alfa == F1_el%alfa)

 END FUNCTION compare_F1el_F2el

 !! @description: Copy a elementary operator F1_el to 
 !!               another elementary operator F2_el
 !! @param:     F1_el    The operator which will be copied 
 !! @param:     F2_el    The operator in which F1_el will be copied 
 !! @param:     idf      The new id of the operation, optional
 !! @param:     idq      The new id of the coordinate, optional
 !! @param:     alfa     The new power of the operator, optional
 !! @param:     coeff    The new coefficient in front of the operator, optional
 subroutine copy_F1_el_into_F2_el(F1_el, F2_el, idf, idq, alfa, frac_alfa, indexq, coeff)

   type(opel),                    intent(in)    :: F1_el
   type(opel),                    intent(inout) :: F2_el
   integer, optional,             intent(in)    :: idf
   integer, optional,             intent(in)    :: idq
   integer, optional,             intent(in)    :: alfa
   TYPE(FracInteger), optional,   intent(in)    :: frac_alfa
   integer, optional,             intent(in)    :: indexq
   complex(kind=Rkind), optional, intent(in)    :: coeff

   TYPE(FracInteger) :: alfa_loc
   character (len=*), parameter :: routine_name="copy_F1_el_into_F2_el"

   call get_opel(Fel = F2_el, idf = F1_el%idf, idq = F1_el%idq, &
   & alfa = F1_el%alfa, indexq = F1_el%indexq, coeff = F1_el%coeff)

   if(present(idf)) then
     call get_opel(Fel = F2_el, idf = idf, idq = F2_el%idq, &
     & alfa = F2_el%alfa, indexq =  F2_el%indexq, coeff = F2_el%coeff)
   end if
   if(present(idq)) then
     call get_opel(Fel = F2_el, idf = F2_el%idf, idq = idq, &
     & alfa = F2_el%alfa, indexq =  F2_el%indexq, coeff = F2_el%coeff)
   end if

   if(present(alfa) .AND. present(frac_alfa)) then
     STOP 'both alfa and frac_alfa are present'
   ELSE if(present(alfa)) THEN
     alfa_loc = alfa
     call get_opel(Fel = F2_el, idf = F2_el%idf, idq = F2_el%idq, &
     & alfa = alfa_loc, indexq =  F2_el%indexq, coeff = F2_el%coeff)
   ELSE if(present(frac_alfa)) THEN
     call get_opel(Fel = F2_el, idf = F2_el%idf, idq = F2_el%idq, &
     & alfa = frac_alfa, indexq =  F2_el%indexq, coeff = F2_el%coeff)
   end if

   if(present(indexq)) then
     call get_opel(Fel = F2_el, idf = F2_el%idf, idq = F2_el%idq, &
     & alfa = F2_el%alfa, indexq =  indexq, coeff = F2_el%coeff)
   end if

   if(present(coeff)) then
     call get_opel(Fel = F2_el, idf = F2_el%idf, idq = F2_el%idq, &
     & alfa = F2_el%alfa, indexq =  F2_el%indexq, coeff = coeff)
   end if

 end subroutine copy_F1_el_into_F2_el

 subroutine R_TO_OpEl(OpEl1,R)

   type(opel),                 intent(inout)  :: OpEl1
   real(kind=Rkind),           intent(in)     :: R

   TYPE(FracInteger) :: alfa_loc
   character (len=*), parameter :: routine_name="R_TO_OpEl"

   IF (R == ZERO) THEN ! zero
     alfa_loc = 0
     call get_opel(Fel = OpEl1, idf = 0, idq = 0, alfa = alfa_loc, indexq = 0, &
                   coeff = czero)
   ELSE ! identity Op x R
     alfa_loc = 1
     call get_opel(Fel = OpEl1, idf = 1, idq = 0, alfa = alfa_loc, indexq = 0, &
                   coeff = cmplx(R,zero,kind=Rkind) )
   END IF

 end subroutine R_TO_OpEl

 subroutine C_TO_OpEl(OpEl1,C)

   type(opel),                 intent(inout)  :: OpEl1
   complex(kind=Rkind),        intent(in)     :: C

   TYPE(FracInteger) :: alfa_loc
   character (len=*), parameter :: routine_name="C_TO_OpEl"

   IF (C == CZERO) THEN
     alfa_loc = 0
     call get_opel(Fel = OpEl1, idf = 0, idq = 0, alfa = alfa_loc, indexq = 0, &
                   coeff = czero)
   ELSE ! identity Op x C
     alfa_loc = 1
     call get_opel(Fel = OpEl1, idf = 1, idq = 0, alfa = alfa_loc, indexq = 0, &
                   coeff = C)
   END IF

 end subroutine C_TO_OpEl

   function R_times_OpEl(R, FOpEl) result(Fres)
     type(opel)    :: Fres
     real (kind = Rkind),    intent(in) :: R
     type(opel),       intent(in)       :: FOpEl

   character (len=*), parameter :: routine_name='R_times_OpEl'

   Fres = FOpEl
   Fres%coeff = FOpEl%coeff * cmplx(R,zero,kind=Rkind)

   end function R_times_OpEl
   function OpEl_times_R(FOpEl,R) result(Fres)
     type(opel)    :: Fres
     real (kind = Rkind),    intent(in) :: R
     type(opel),       intent(in)       :: FOpEl

   character (len=*), parameter :: routine_name='OpEl_times_R'

   Fres = FOpEl
   Fres%coeff = FOpEl%coeff * cmplx(R,zero,kind=Rkind)

   end function OpEl_times_R
   function C_times_OpEl(C, FOpEl) result(Fres)
     type(opel)    :: Fres
     complex (kind = Rkind),    intent(in) :: C
     type(opel),       intent(in)       :: FOpEl

   character (len=*), parameter :: routine_name='C_times_OpEl'

   Fres = FOpEl
   Fres%coeff = FOpEl%coeff * C

   end function C_times_OpEl
   function OpEl_times_C(FOpEl,C) result(Fres)
     type(opel)    :: Fres
     complex (kind = Rkind),    intent(in) :: C
     type(opel),       intent(in)       :: FOpEl

   character (len=*), parameter :: routine_name='OpEl_times_C'

   Fres       = FOpEl
   Fres%coeff = FOpEl%coeff * C

   end function OpEl_times_C

 subroutine Sort_TabOpEl(TabOpEl)

   type(opel),                 intent(inout)  :: TabOpEl(:)


   integer :: i,j,n
   type(opel) :: tempOpEl
   character (len=*), parameter :: routine_name="Sort_TabOpEl"

   n = size(TabOpEl)



   DO i=1,n
   DO j=i+1,n
     IF ( TabOpEl(i)%idf > TabOpEl(j)%idf ) THEN
       tempOpEl   = TabOpEl(j)
       TabOpEl(j) = TabOpEl(i)
       TabOpEl(i) = tempOpEl
     END IF
   END DO
   END DO

 end subroutine Sort_TabOpEl

 subroutine Merge_TabOpEl(TabOpEl)

   type(opel),                 intent(inout)  :: TabOpEl(:)


   integer :: i,j,n
   type(opel) :: tempOpEl
   character (len=*), parameter :: routine_name="Merge_TabOpEl"

   n = size(TabOpEl)

   DO i=1,n-1
     IF (TabOpEl(i)%idf >= 11 .AND. TabOpEl(i)%idf <= 23) CYCLE
     IF (TabOpEl(i)%idf == TabOpEl(i+1)%idf) THEN
       ! we use get_opel, to have the correct opname
       CALL get_opel(TabOpEl(i+1),                                      &
                     idf    = TabOpEl(i)%idf,                           &
                     idq    = TabOpEl(i)%idq,                           &
                     alfa   = TabOpEl(i+1)%alfa  + TabOpEl(i)%alfa,     &
                     indexq = TabOpEl(i)%indexq,                        &
                     coeff  = TabOpEl(i+1)%coeff * TabOpEl(i)%coeff )

       TabOpEl(i) = cone
     END IF
   END DO

 end subroutine Merge_TabOpEl

   !! @description: Write an elementary ,
   !! @param:       Fel       The operator (type: opel).
   !! @param:       filename  Name of the output file.
   !! @param:       header    If present write comment in the biginning of the
   !!                         of the file
   SUBROUTINE write_opel(Fel,  i_file, header, append, close_file)
     type(opel),                intent(in)       :: Fel
     integer, optional,         intent(in)       :: i_file
     logical, optional,         intent(in)       :: header
     logical, optional,         intent(in)       :: append
     logical, optional,         intent(in)       :: close_file

     integer                        :: error
     integer                        :: i_open
     logical                        :: header_loc

     character (len=*), parameter   :: routine_name='write_opel'

     header_loc = .FALSE.
     if (present(header)) header_loc = header


     if(present(i_file)) then
       i_open = i_file
     else
       i_open = out_unitp
     end if

     if (header_loc) then
       write(i_open, '(A)')  "   idf          iqd        alfa       alfa_den       &
              &indexq           op_name                  coeff"
     end if
     write(i_open,"(3x, 5(I2,10x), 3x, A, (E13.4,' Ix ',E13.4))" ) Fel%idf, Fel%idq, Fel%alfa%num, &
                 Fel%alfa%den, Fel%indexq, Fel%opname,Fel%coeff
   END SUBROUTINE write_opel

   FUNCTION get_coeff_OF_OpEl(Fel)
     type(opel),                intent(in)       :: Fel
     complex (kind=Rkind)                        :: get_coeff_OF_OpEl

     character (len=*), parameter   :: routine_name='get_coeff_OF_OpEl'

     get_coeff_OF_OpEl = Fel%coeff

   END FUNCTION get_coeff_OF_OpEl

   SUBROUTINE Set_coeff_OF_OpEl_TO_ONE(Fel)
     type(opel),                intent(inout)       :: Fel

     character (len=*), parameter   :: routine_name='Set_coeff_OF_OpEl_TO_ONE'

     Fel%coeff = cone

   END SUBROUTINE Set_coeff_OF_OpEl_TO_ONE

   SUBROUTINE get_pqJL_OF_OpEl(pq,J,L,Fel)
     type(opel),                intent(in)       :: Fel
     integer,                   intent(inout)    :: pq(2),J(2),L(2)

     character (len=*), parameter   :: routine_name='get_pqJL_OF_OpEl'

   pq = 0
   L  = 0
   J  = 0
   select case (Fel%idf)
     case(4)
       pq(1) = Fel%indexq
       IF (Fel%alfa == 2) pq(2) = Fel%indexq
     case(9,10,11) ! Jx, Jy, Jz
       J(1) = Fel%indexq
       IF (Fel%alfa == 2) J(2) = Fel%indexq
     case(12,13,14,15,16,17,18,19,20,21,22,23)
       pq(1) = Fel%indexq
       ! no alpha on pq
     case(24) ! Lx
       L(1) = -4
       IF (Fel%alfa == 2) L(2) = -4
     case(25) ! Ly
       L(1) = -5
       IF (Fel%alfa == 2) L(2) = -5
     case(26) ! Lz
       L(1) = -6
       IF (Fel%alfa == 2) L(2) = -6
     end select

   END SUBROUTINE get_pqJL_OF_OpEl

   !! @description: numerical calculation of an elementary Op,
   !! @param:       Fel       The operator (type: opel).
   !! @param:       Qval      value of the coordinate associated with Fel
   !! @param:       ValOpEl   value of the operator, Fel
   SUBROUTINE get_NumVal_OpEl(ValOpEl,Qval,Fel)
     type(opel),                intent(in)       :: Fel
     real(kind=Rkind),          intent(in)       :: Qval
     complex(kind=Rkind),       intent(inout)    :: ValOpEl

     real(kind=Rkind)       :: ralfa,RvalOp

     character (len=*), parameter   :: routine_name='get_NumVal_OpEl'

   ralfa = Frac_TO_Real(Fel%alfa)

   select case (Fel%idf)
     case(0) ! zero
       RValOp = zero
     case(1) ! id
       RValOp = one
     case(2) ! qval**alfa
       RValOp = qval**ralfa
     case(3)
       RValOp = sqrt(one-qval**2)**ralfa
     case(4)
       ValOpEl = Fel%coeff**ralfa
     case(5)
       RValOp = cos(qval)**ralfa
     case(6)
       RValOp = sin(qval)**ralfa
     case(7)
       RValOp = tan(qval)**ralfa
     case(8)
       RValOp = ONE / tan(qval)**ralfa
     case(9,10,11)
       RValOp = ONE
     case(12)
       RValOp = qval**ralfa
     case(13)
       RValOp = qval**ralfa
     case(14)
       RValOp = cos(qval)**ralfa
     case(15)
       RValOp = cos(qval)**ralfa
     case(16)
       RValOp = sin(qval)**ralfa
     case(17)
       RValOp = sin(qval)**ralfa
     case(18,19)
       RValOp = tan(qval)**ralfa
     case(20,21)
       RValOp = ONE/tan(qval)**ralfa
     case(22,23)
       RValOp = sqrt(one-qval**2)**ralfa
     case(24,25,26)
       RValOp = ONE
     case default
       write(out_unitp,*) 'ERROR in ',routine_name
       write(out_unitp,*) 'idf=',Fel%idf
       write(out_unitp,*) "This idf is not registered for an elementary operator"
       write(out_unitp,*) "   illegal value of idf . (Internal Bug)"
       STOP
     end select

     IF (Fel%idf /= 4) THEN ! already done for idf=4
        ValOpEl =  Fel%coeff * real(RvalOp,kind=Rkind)
     END IF

   END SUBROUTINE get_NumVal_OpEl
   SUBROUTINE Change_PQ_OF_OpEl_TO_Id_OF_OpEl(F_OpEl)
     type(opel),                intent(inout)       :: F_OpEl    ! PQ or Jx ...=> Id or PQ Q^alfa => Q^alfa

     type(opel)       :: tmp_OpEl

   character (len=*), parameter   :: routine_name='Change_PQ_OF_OpEl_TO_Id_OF_OpEl'


   select case (F_OpEl%idf)
   case(0,1,2,3,5,6,7,8)
     CONTINUE ! nothing
  case(4,9,10,11) ! PQ, Jx,Jy,Jz  (we keep the coef)
     F_OpEl = F_OpEl%coeff

   case(12,13) ! PQ Q^alfa   => Q^alfa
     CALL copy_F1_el_into_F2_el(F_OpEl, tmp_OpEl, idf=2)
     F_OpEl = tmp_OpEl
   case(14,15) ! PQ cos(Q)^alfa => cos(Q)^alfa
     CALL copy_F1_el_into_F2_el(F_OpEl, tmp_OpEl, idf=5)
     F_OpEl = tmp_OpEl

   case(16,17) ! PQ sin(Q)^alfa
     CALL copy_F1_el_into_F2_el(F_OpEl, tmp_OpEl, idf=6)
     F_OpEl = tmp_OpEl

   case(18,19) ! PQ tan(Q)^alfa
     CALL copy_F1_el_into_F2_el(F_OpEl, tmp_OpEl, idf=7)
     F_OpEl = tmp_OpEl

   case(20,21) ! PQ cot(Q)^alfa
     CALL copy_F1_el_into_F2_el(F_OpEl, tmp_OpEl, idf=8)
     F_OpEl = tmp_OpEl

   case(22,23) ! PQ qs(Q)^alfa Rq: qs(Q)=sqrt(1-q^2)
     CALL copy_F1_el_into_F2_el(F_OpEl, tmp_OpEl, idf=3)
     F_OpEl = tmp_OpEl

   case(24,25,26) ! Lx,Ly,Lz
     write(out_unitp,*) 'ERROR in ',routine_name
     write(out_unitp,*) 'idf=', F_OpEl%idf
     write(out_unitp,*) "Impossible to deal with Lx, Ly, Lz !!"
     STOP

   case default
     write(out_unitp,*) 'ERROR in ',routine_name
     write(out_unitp,*) 'idf=', F_OpEl%idf
     write(out_unitp,*) "This idf is not registered for an elementary operator"
     write(out_unitp,*) "   illegal value of idf . (Internal Bug)"
     STOP
   end select

   END SUBROUTINE Change_PQ_OF_OpEl_TO_Id_OF_OpEl
   SUBROUTINE Split_OpEl_TO_SplitOpEl(F_OpEl,SplitOpEl)
     type(opel),                intent(in)       :: F_OpEl
     type(opel),  allocatable, intent(inout)     :: SplitOpEl(:) ! product of new OpEl
                                                                 ! PQ Q^alfa => PQ   and   Q^alfa

   integer :: idf1,idf2
   character (len=*), parameter   :: routine_name='Split_OpEl_TO_SplitOpEl'


   CALL dealloc_NParray(SplitOpEl,'SplitOpEl in ',routine_name)

   select case (F_OpEl%idf)
   case(0,1,2,3,4,5,6,7,8,9,10,11,  24,25,26)
     idf1 = F_OpEl%idf ; idf2 = -1

   case(12) ! PQ Q^alfa
     idf1 = 4 ; idf2 = 2
   case(13) ! Q^alfa PQ
     idf1 = 2 ; idf2 = 4

   case(14) ! PQ cos(Q)^alfa
     idf1 = 4 ; idf2 = 5
   case(15) ! cos(Q)^alfa PQ
     idf1 = 5 ; idf2 = 4

   case(16) ! PQ sin(Q)^alfa
     idf1 = 4 ; idf2 = 6
   case(17) ! sin(Q)^alfa PQ
     idf1 = 6 ; idf2 = 4

   case(18) ! PQ tan(Q)^alfa
     idf1 = 4 ; idf2 = 7
   case(19) ! tan(Q)^alfa PQ
     idf1 = 7 ; idf2 = 4

   case(20) ! PQ cot(Q)^alfa
     idf1 = 4 ; idf2 = 8
   case(21) ! cot(Q)^alfa PQ
     idf1 = 8 ; idf2 = 4

   case(22) ! PQ qs(Q)^alfa Rq: qs(Q)=sqrt(1-q^2)
     idf1 = 4 ; idf2 = 3
   case(23) ! qs(Q)^alfa PQ
     idf1 = 3 ; idf2 = 4

   case default
     write(out_unitp,*) 'ERROR in ',routine_name
     write(out_unitp,*) 'idf=', F_OpEl%idf
     write(out_unitp,*) "This idf is not registered for an elementary operator"
     write(out_unitp,*) "   illegal value of idf . (Internal Bug)"
     STOP
   end select

     IF (idf2 == -1) THEN
       CALL alloc_NParray(SplitOpEl,(/1/),'SplitOpEl in ',routine_name)
       SplitOpEl(1) = F_OpEl

     ELSE IF (idf1 == 4) THEN
       CALL alloc_NParray(SplitOpEl,(/2/),'SplitOpEl in',routine_name)

       CALL get_opel(SplitOpEl(1),                                      &
                     idf    = 4,                                        &
                     idq    = F_OpEl%idq,                               &
                     alfa   = FracInteger(1,1),                         &
                     indexq = F_OpEl%indexq,                            &
                     coeff  = F_OpEl%coeff)  ! PQ

       CALL get_opel(SplitOpEl(2),                                      &
                     idf    = idf2,                                     &
                     idq    = F_OpEl%idq,                               &
                     alfa   = F_OpEl%alfa,                              &
                     indexq = F_OpEl%indexq,                            &
                     coeff  = CONE)  ! f(Q)^alfa with idf2
     ELSE ! idf2 == 4
       CALL alloc_NParray(SplitOpEl,(/2/),'SplitOpEl in',routine_name)

       CALL get_opel(SplitOpEl(1),                                      &
                     idf    = idf1,                                     &
                     idq    = F_OpEl%idq,                               &
                     alfa   = F_OpEl%alfa,                              &
                     indexq = F_OpEl%indexq,                            &
                     coeff  = CONE)  ! f(Q)^alfa with idf1

       CALL get_opel(SplitOpEl(2),                                      &
                     idf    = 4,                                        &
                     idq    = F_OpEl%idq,                               &
                     alfa   = FracInteger(1,1),                         &
                     indexq = F_OpEl%indexq,                            &
                     coeff  = F_OpEl%coeff)  ! PQ

     END IF

   END SUBROUTINE Split_OpEl_TO_SplitOpEl

   SUBROUTINE Der1_OF_d0OpEl_TO_d1OpEl(d0OpEl,d1OpEl)
     type(opel),                intent(in)       :: d0OpEl
     type(opel), allocatable, intent(inout)      :: d1OpEl(:) ! product of OpEl
                                                             ! cos(Q)^alfa => -Sin(Q) * cos(Q)^alfa-1
!  It doesn't work when OpEl contains a P or J or L

   TYPE(FracInteger)    :: alfa_prim
   real(kind=Rkind)     :: ralfa
   complex(kind=Rkind)  :: coeff
   character (len=*), parameter   :: routine_name='Der1_OF_d0OpEl_TO_d1OpEl'

   CALL alloc_NParray(d1OpEl,(/2/),'d1OpEl in ',routine_name)

   ralfa = Frac_TO_Real(d0OpEl%alfa)


! For most of the function
!                        d1OpEl(1)           d1OpEl(2)
! Op(Q) = f(Q)^alfa =>   f'(Q)     *    ( alfa. f(Q)^(alfa-1)  )

   coeff = d0OpEl%coeff * cmplx(ralfa,zero,kind=Rkind)

   alfa_prim = d0OpEl%alfa-1
   CALL get_opel(d1OpEl(2),                                             &
                     idf    = d0OpEl%idf,                               &
                     idq    = d0OpEl%idq,                               &
                     alfa   = alfa_prim,                                &
                     indexq = d0OpEl%indexq,                            &
                     coeff  = coeff)

   select case (d0OpEl%idf)
   case(0,1) ! constant
     d1OpEl(1) = CZERO

   case(2) ! q^alfa => alfa * q^(alfa-1)
     d1OpEl(1) = CONE

   case(3) ! sqrt(1-q^2)^alfa => -alfa x * sqrt(1-q^2)^alfa^(alfa-2)
     CALL get_opel(d1OpEl(1),                                           &
                     idf    = 2,                                        &
                     idq    = d0OpEl%idq,                               &
                     alfa   = FracInteger(1,1),                         &
                     indexq = d0OpEl%indexq,                            &
                     coeff  = cone)

     coeff = -d0OpEl%coeff * cmplx(ralfa,zero,kind=Rkind)

     alfa_prim = d0OpEl%alfa-2
     CALL get_opel(d1OpEl(2),                                           &
                     idf    = d0OpEl%idf,                               &
                     idq    = d0OpEl%idq,                               &
                     alfa   = alfa_prim,                                &
                     indexq = d0OpEl%indexq,                            &
                     coeff  = coeff)

   case(5) ! cos^alfa => -sin   *    alfa cos^(alfa-1)
     CALL get_opel(d1OpEl(1),                                           &
                     idf    = 6,                                        &
                     idq    = d0OpEl%idq,                               &
                     alfa   = FracInteger(1,1),                         &
                     indexq = d0OpEl%indexq,                            &
                     coeff  = -cone)

   case(6) ! sin^alfa => cos   *    alfa sin^(alfa-1)
     CALL get_opel(d1OpEl(1),                                           &
                     idf    = 5,                                        &
                     idq    = d0OpEl%idq,                               &
                     alfa   = FracInteger(1,1),                         &
                     indexq = d0OpEl%indexq,                            &
                     coeff  = cone)

   case(7) ! tan^alfa => cos^-2   *    alfa tan^(alfa-1)
     CALL get_opel(d1OpEl(1),                                           &
                     idf    = 5,                                        &
                     idq    = d0OpEl%idq,                               &
                     alfa   = FracInteger(-2,1),                        &
                     indexq = d0OpEl%indexq,                            &
                     coeff  = cone)

   case(8) ! cot^alfa => -sin^-2   *    alfa cot^(alfa-1)
     CALL get_opel(d1OpEl(1),                                           &
                     idf    = 6,                                        &
                     idq    = d0OpEl%idq,                               &
                     alfa   = FracInteger(-2,1),                        &
                     indexq = d0OpEl%indexq,                            &
                     coeff  = -cone)

   case(9,10,11) ! Jx,Jy,Jz
     d1OpEl(1) = d0OpEl
     d1OpEl(2) = CONE


   case default
     write(out_unitp,*) 'ERROR in ',routine_name
     write(out_unitp,*) 'idf=', d0OpEl%idf
     write(out_unitp,*) "This idf is not registered for an elementary operator"
     write(out_unitp,*) "   illegal value of idf . (Internal Bug)"
     STOP
   end select



   END SUBROUTINE Der1_OF_d0OpEl_TO_d1OpEl

 end module mod_Tana_OpEl
