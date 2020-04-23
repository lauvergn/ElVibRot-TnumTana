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
      PROGRAM Tana_test
      USE mod_system
      USE mod_Constant
      USE mod_Coord_KEO
      USE mod_Tnum
      USE mod_PrimOp
      IMPLICIT NONE


!     - parameters for para_Tnum -----------------------
      TYPE (constant)  :: const_phys
      TYPE (zmatrix)   :: mole
      TYPE (Tnum)      :: para_Tnum
      TYPE (param_PES) :: para_PES

      TYPE(sum_opnd)             :: TWOxKEO


      integer :: nderiv,nb_act

      real (kind=Rkind), allocatable :: Qact(:)

!     ------------------------------------------------------

!     - working parameters ------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Tana_test'

CALL test_String_TO_Sum_OpnD()
STOP
!===========================================================
!===========================================================
!      !para_mem%mem_debug = .TRUE.
!      CALL versionEVRT(.TRUE.)
!
!      !-----------------------------------------------------------------
!      !     - read the coordinate tansformations :
!      !     -   zmatrix, polysperical, bunch...
!      !     ------------------------------------------------------------
!      CALL Read_mole(mole,para_Tnum,const_phys)
!      !     ------------------------------------------------------------
!      !-----------------------------------------------------------------
!
!      !-----------------------------------------------------------------
!      !     - read coordinate values -----------------------------------
!      !     ------------------------------------------------------------
!      CALL read_RefGeom(mole,para_Tnum)
!      !     ------------------------------------------------------------
!      !-----------------------------------------------------------------
!
!      !-----------------------------------------------------------------
!      !     ---- TO finalize the coordinates (NM) and the KEO ----------
!      !     ------------------------------------------------------------
!      CALL Finalyze_TnumTana_Coord_PrimOp(para_Tnum,mole,para_PES)
!      !-----------------------------------------------------------------
!===========================================================
!===========================================================

!-------------------------------------------------
! test Tana library:
! 1) Elementary operator (OpEl)
!CALL test_Opel()

! 2) 1D operator (Op1D)
!CALL test_Op1D()

! 3) sum of 1D operator (Sum_OF_Op1D):
!  derivative (result in Sum_OF_Op1D)
!CALL test_Der_OF_Op1D()
!   Sum_OF_Op1D subroutines
CALL test_SumOp1D()
!    Expand Op1D
!CALL test_Expand_Op1D()

! 4) OpnD:
!CALL test_OpnD()

! 5) Sum_OpnD:
!CALL test_Sum_OpnD()
!-------------------------------------------------

  CALL dealloc_zmat(mole)

  write(out_unitp,*) 'END Tana_test'

END PROGRAM Tana_test

SUBROUTINE test_Opel()
  USE mod_system
  USE mod_Constant
  USE mod_Coord_KEO
  USE mod_Tana_OpEl
  IMPLICIT NONE


  TYPE(OpEl)              :: F1el,F2el
  TYPE(OpEl), allocatable :: tab_Fel(:)


  integer :: pq(2),JJ(2),LL(2),idf,idq,alfa,err_Op
  integer :: tab_idq(0:27) = (/ 1,1,2,-3,2,  3,3,3,3,               &
    5,5,5,  2,2,  3,3,3,3,3,3,3,3,  -3,-3,  -5,-5,-5,  2 /)

  logical :: first
  character (len = :), allocatable :: FelName
  integer :: i,j,n,ndim

!     ------------------------------------------------------

!     - working parameters ------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='test_Opel'


!-------------------------------------------------
! test Tana library:
! 1) Elementary operator (OpEl)
!-------------------------------------------------
write(out_unitp,*) '===================================='
write(out_unitp,*) '======== Elementry Op =============='
alfa=0
DO idq=1,8
  write(out_unitp,*) '===================================='
  write(out_unitp,*) 'idq',idq
  F1el = set_opel(4, idq, alfa, indexq=1, coeff=cone,err_el=err_Op)
  IF (err_Op == 0) CALL write_op(F1el,header=.TRUE.)
END DO

idq=-3
write(out_unitp,*) '===================================='
write(out_unitp,*) 'idq',idq
F1el = set_opel(4, idq, alfa, indexq=1, coeff=cone,err_el=err_Op)
IF (err_Op == 0) CALL write_op(F1el,header=.TRUE.)

idq=-7
write(out_unitp,*) '===================================='
write(out_unitp,*) 'idq',idq
F1el = set_opel(4, idq, alfa, indexq=1, coeff=cone,err_el=err_Op)
IF (err_Op == 0) CALL write_op(F1el,header=.TRUE.)
stop

write(out_unitp,*) '===================================='
write(out_unitp,*) '======== Elementry Op =============='
DO idf=0,ubound(tab_idq,dim=1)
  first = .TRUE.
  write(out_unitp,*) '===================================='
  write(out_unitp,*) 'idf',idf

  DO alfa=-2,2
    F1el = set_opel(idf, tab_idq(idf), alfa, indexq=1, coeff=cone,err_el=err_Op)
    CALL write_op(F1el,header=.TRUE.)
    CALL Export_VSCF_OpEl(F1el,'Q0',FelName)
    write(out_unitp,*) idf,tab_idq(idf),alfa,'FelName: ',FelName

    CALL Split_OpEl_TO_SplitOpEl(F1el,tab_Fel)
    write(out_unitp,'(a,"(",f12.6," +Ix ",f12.6,")*")',advance='no') 'Split_Fel: ',product(tab_Fel(:)%coeff)
    DO i=1,size(tab_Fel)-1
      CALL Export_VSCF_OpEl(tab_Fel(i),'Q0',FelName)
      write(out_unitp,'(a)',advance='no') FelName // ' * '
    END DO
    CALL Export_VSCF_OpEl(tab_Fel(size(tab_Fel)),'Q0',FelName)
    write(out_unitp,'(a)',advance='yes') FelName

    CALL get_pqJL_OF_OpEl(pq,JJ,LL,F1el)
    IF (pq(1) == 0) THEN
      CALL Der1_OF_d0OpEl_TO_d1OpEl(F1el,tab_Fel)
      write(out_unitp,'(a,"(",f12.6," +Ix ",f12.6,")*")',advance='no') 'Der1_Of_Fel: ',product(tab_Fel(:)%coeff)
      DO i=1,size(tab_Fel)-1
        CALL Export_VSCF_OpEl(tab_Fel(i),'Q0',FelName)
        write(out_unitp,'(a)',advance='no') FelName // ' * '
      END DO
      CALL Export_VSCF_OpEl(tab_Fel(size(tab_Fel)),'Q0',FelName)
      write(out_unitp,'(a)',advance='yes') FelName
    END IF

    !CALL write_opel(F1el, header=first)
    first = .FALSE.
    write(out_unitp,*)
  END DO
END DO


END SUBROUTINE test_Opel
SUBROUTINE test_Op1D()
  USE mod_system
  USE mod_Constant
  USE mod_Coord_KEO
  USE mod_Tana_OpEl
  USE mod_Tana_Op1D
  IMPLICIT NONE


      TYPE(OpEl)             :: F1el
      TYPE(Op1D)             :: F11D,F21D

      character (len = :), allocatable :: FelName
      integer :: i,idq,idf

!     - working parameters ------------------------------------------
      integer :: err_mem,memory,err_Op
      character (len=*), parameter :: name_sub='test_Op1D'


!-------------------------------------------------
! 1) 1D operator (Op1D)
!-------------------------------------------------
write(out_unitp,*) '===================================='
write(out_unitp,*) '======== 1D Op ====================='
idq=2
idf=2

write(out_unitp,*)
write(out_unitp,*) '------------------------------'
F11D = cone
F1el = set_opel(idf, idq, alfa=2, indexq=1, coeff=cone,err_el=err_Op)
F11D = F11D * F1el
CALL Export_VSCF_Op1D(F11D,'Q0',FelName)
write(out_unitp,*) 'Op1D = 1 x Q^2 ?'
write(out_unitp,'(a,a)') '1DOp: ',FelName
CALL flush_perso(out_unitp)

write(out_unitp,*)
write(out_unitp,*) '------------------------------'
write(out_unitp,*) 'Op1D = Op1D x Q^-1 = Q ?'
F1el = set_opel(idf, idq, alfa=-1, indexq=1, coeff=cone,err_el=err_Op)
F11D = F11D * F1el
CALL Export_VSCF_Op1D(F11D,'Q0',FelName)
write(out_unitp,'(a,a)') '1DOp: ',FelName
CALL flush_perso(out_unitp)

write(out_unitp,*)
write(out_unitp,*) '------------------------------'
write(out_unitp,*) 'Op1D = P x Op1D = P Q ?'
F1el = set_opel(4, idq, alfa=1, indexq=1, coeff=cone,err_el=err_Op)
F11D = F1el * F11D
CALL Export_VSCF_Op1D(F11D,'Q0',FelName)
write(out_unitp,'(a,a)') '1DOp: ',FelName
CALL flush_perso(out_unitp)

write(out_unitp,*)
write(out_unitp,*) '------------------------------'
write(out_unitp,*) 'Op1D = (Q^-2xP) x Op1D = Q^-2P P Q ?'
F1el = set_opel(13, idq, alfa=-2, indexq=1, coeff=cone,err_el=err_Op)
F11D = F1el * F11D
CALL Export_VSCF_Op1D(F11D,'Q0',FelName)
write(out_unitp,'(a,a)') '1DOp: ',FelName
CALL flush_perso(out_unitp)

write(out_unitp,*)
write(out_unitp,*) '------------------------------'
write(out_unitp,*) 'Op1D = 1 x Op1D = Q^-2P P Q ?'
F1el = set_opel(1, idq, alfa=2, indexq=1, coeff=cone,err_el=err_Op)
F11D = F11D * F1el
CALL Export_VSCF_Op1D(F11D,'Q0',FelName)
write(out_unitp,'(a,a)') '1DOp: ',FelName
CALL flush_perso(out_unitp)

write(out_unitp,*)
write(out_unitp,*) '------------------------------'
write(out_unitp,*) 'Simplify(Op1D) = Op1D = Q^-2P P Q ?'
CALL Simplify_Op1D(F11D)
CALL Export_VSCF_Op1D(F11D,'Q0',FelName)
write(out_unitp,'(a,a)') 'Simplify 1DOp: ',FelName
CALL flush_perso(out_unitp)

write(out_unitp,*)
write(out_unitp,*) '------------------------------'
write(out_unitp,*) 'Split(Op1D) = Q^-2 x P^2 x Q ?'
CALL Split_Op1D_TO_SplitOp1D(F11D,F21D) ; F11D=F21D
CALL Export_VSCF_Op1D(F11D,'Q0',FelName)
write(out_unitp,'(a,a)') 'Split+simplify 1DOp: ',FelName
CALL flush_perso(out_unitp)

CALL delete_op(F11D)
CALL delete_op(F21D)


END SUBROUTINE test_Op1D
SUBROUTINE test_Der_OF_Op1D()
  USE mod_system
  USE mod_Constant
  USE mod_Coord_KEO
  USE mod_Tana_OpEl
  USE mod_Tana_Op1D
  IMPLICIT NONE


  TYPE(OpEl)             :: F1el
  TYPE(Op1D)             :: F11D
  TYPE(Sum_OF_Op1D)      :: tab_F1D

  character (len = :), allocatable :: FelName
  integer :: i,idq


  ! - working parameters ------------------------------------------
  integer :: err_mem,memory,err_Op
  character (len=*), parameter :: name_sub='test_Der_OF_Op1D'

write(out_unitp,*) '===================================='
write(out_unitp,*) '========der1 of 1D Op =============='
idq=3


F11D = cone
F1el = set_opel(5, idq, alfa=2, indexq=1, coeff=cone,err_el=err_Op)
F11D = F11D * F1el
CALL Export_VSCF_Op1D(F11D,'Q0',FelName)
write(out_unitp,*) 'Op1D = 1 x cos(Q)^2 ?'
write(out_unitp,'(a,a)') '1DOp: ',FelName
CALL flush_perso(out_unitp)

write(out_unitp,*) 'Op1D = Op1D x sin(Q)^-1 = cos^2 sin^-1 ?'
F1el = set_opel(6, idq, alfa=-1, indexq=1, coeff=cone,err_el=err_Op)
F11D = F11D * F1el
CALL Export_VSCF_Op1D(F11D,'Q0',FelName)
write(out_unitp,'(a,a)') '1DOp: ',FelName
CALL flush_perso(out_unitp)

write(out_unitp,*) ' First derivative:'
tab_F1D = Der1_OF_d0Op1D(F11D)
DO i=1,size(tab_F1D%Sum_Op1D)
  CALL Export_VSCF_Op1D(tab_F1D%Sum_Op1D(i),'Q0',FelName)
  write(out_unitp,*) i,product(tab_F1D%Sum_Op1D(i)%prod_opel(:)%coeff),' * ',FelName
END DO
CALL delete_op(tab_F1D)


write(out_unitp,*) '===================================='
write(out_unitp,*) '========der2 of 1D Op =============='

write(out_unitp,*) ' Second derivative:'
tab_F1D = Der2_OF_d0Op1D(F11D)

DO i=1,size(tab_F1D%Sum_Op1D)
  CALL Export_VSCF_Op1D(tab_F1D%Sum_Op1D(i),'Q0',FelName)
  write(out_unitp,*) i,product(tab_F1D%Sum_Op1D(i)%prod_opel(:)%coeff),' * ',FelName
END DO

CALL delete_op(F11D)
CALL delete_op(tab_F1D)

END SUBROUTINE test_Der_OF_Op1D

SUBROUTINE test_SumOp1D()
  USE mod_system
  USE mod_Constant

  USE mod_Coord_KEO
  USE mod_Tana_OpEl
  USE mod_Tana_Op1D
  IMPLICIT NONE

  TYPE(Op1D)             :: F11D,F21D
  TYPE(Sum_OF_Op1D)      :: Sum1_1D,Sum2_1D,Sum3_1D


  character (len = :), allocatable :: FelName
  integer :: i,idq

!     ------------------------------------------------------

!     - working parameters ------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='test_SumOp1D'

write(out_unitp,*) '===================================='
write(out_unitp,*) '======== Simplify sin^2 ============'
idq=-3
write(out_unitp,*) '------------------------------------'
F11D = set_opel(3, idq, alfa=9, indexq=1, coeff=cone) * & ! sin^9
       set_opel(2, idq, alfa=3, indexq=1, coeff=cone) * & ! cos^3
       set_opel(4, idq, alfa=1, indexq=1, coeff=cone) * & ! P
       set_opel(3, idq, alfa=2, indexq=1, coeff=cone) * & ! sin^2
       set_opel(2, idq, alfa=1, indexq=1, coeff=cone)     ! cos
CALL Export_VSCF_Op1D(F11D,'Q0',FelName)
write(out_unitp,*) 'Op1D = su^9 * u^3 * P * su^2 * u ?'
write(out_unitp,'(a,a)') '1DOp: ',FelName

CALL Expand_Sin2_IN_Op1D_TO_SumOp1D(F11D,Sum1_1D)

write(out_unitp,*) 'Expand_Sin2_IN_Op1D_TO_SumOp1D'
CALL write_op(Sum1_1D)

write(out_unitp,*) '===================================='
write(out_unitp,*) '======== Sum OF 1D Op =============='

idq=2
write(out_unitp,*) '------------------------------------'
F11D = set_opel(2, idq, alfa=-2, indexq=1, coeff=cone) * &
       set_opel(2, idq, alfa= 3, indexq=1, coeff=cone)
CALL Export_VSCF_Op1D(F11D,'Q0',FelName)
write(out_unitp,*) 'Op1D = Q^-2 * Q^3 ?'
write(out_unitp,'(a,a)') '1DOp: ',FelName

write(out_unitp,*) '------------------------------------'
F21D = set_opel(2, idq, alfa=-2, indexq=1, coeff=cone) * &
       set_opel(4, idq, alfa= 1, indexq=1, coeff=cone)
CALL Export_VSCF_Op1D(F21D,'Q0',FelName)
write(out_unitp,*) 'Op1D = Q^-2 * P  ?'
write(out_unitp,'(a,a)') '1DOp: ',FelName

CALL flush_perso(out_unitp)

write(out_unitp,*) '------------------------------------'
Sum1_1D = F11D + F11D
write(out_unitp,*) 'F11D+F11D'
CALL write_op(Sum1_1D)

write(out_unitp,*) '------------------------------------'
Sum1_1D = F11D + F21D
write(out_unitp,*) 'F11D+F21D'
CALL write_op(Sum1_1D)

write(out_unitp,*) '------------------------------------'
F11D = set_opel(2, idq, alfa=-2, indexq=1, coeff=cone)
Sum1_1D = Sum1_1D + F11D
write(out_unitp,*) 'Sum1_1D+F11D(Q^-2)'
CALL write_op(Sum1_1D)

write(out_unitp,*) '------------------------------------'
F11D = set_opel(2, idq, alfa=2, indexq=1, coeff=cone)
Sum1_1D = F11D*Sum1_1D
write(out_unitp,*) 'F1D(Q^2) * Sum1_1D'
CALL write_op(Sum1_1D)

write(out_unitp,*) '------------------------------------'
F11D = set_opel(2, idq, alfa=-2, indexq=1, coeff=cone)
Sum1_1D = Sum1_1D*F11D
write(out_unitp,*) 'Sum1_1D * F11D(Q^-2)'
CALL write_op(Sum1_1D)

write(out_unitp,*) '------------------------------------'
F11D = set_opel(2, idq, alfa=-2, indexq=1, coeff=cone)
F21D = set_opel(4, idq, alfa=2, indexq=1, coeff=cone)
Sum2_1D = F11D+F21D
write(out_unitp,*) 'Sum2_1D'
CALL write_op(Sum2_1D)
Sum3_1D = Sum1_1D+Sum2_1D
write(out_unitp,*) 'Sum3_1D'
CALL write_op(Sum3_1D)

write(out_unitp,*) '------------------------------------'
Sum1_1D = Sum1_1D+Sum2_1D
write(out_unitp,*) 'Sum1_1D = Sum1_1D+Sum2_1D'
CALL write_op(Sum1_1D)

write(out_unitp,*) '------------------------------------'
F11D = set_opel(2, idq, alfa=-2, indexq=1, coeff=cone)
F21D = set_opel(2, idq, alfa=2, indexq=1, coeff=cone)
Sum1_1D = F11D
write(out_unitp,*) 'Sum1_1D'
CALL write_op(Sum1_1D)
Sum2_1D = F11D+F21D
write(out_unitp,*) 'Sum2_1D'
CALL write_op(Sum2_1D)
Sum1_1D = Sum1_1D + F21D * Der1_OF_d0SumOp1D(Sum2_1D)
write(out_unitp,*) 'Sum1_1D = Sum1_1D + F21D * Der1_OF_d0SumOp1D(Sum2_1D)'
CALL write_op(Sum1_1D)

END SUBROUTINE test_SumOp1D
SUBROUTINE test_Expand_Op1D()
  USE mod_system
  USE mod_Constant

  USE mod_Coord_KEO
  USE mod_Tana_OpEl
  USE mod_Tana_Op1D
  IMPLICIT NONE

  TYPE(Op1D)             :: F11D
  TYPE(Sum_OF_Op1D)      :: tab_F1D


  character (len = :), allocatable :: FelName
  integer :: i,n1,n2,n3,idq
  integer, parameter :: ex = 2
  integer :: list_exP2(3,8) = reshape((/ 0,0,0,  0,0,ex,  0,ex,0,  ex,0,0,  0,ex,ex,  ex,0,ex,  ex,ex,0,  ex,ex,ex /),(/3,8/))

!     - working parameters ------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='test_Expand_Op1D'

write(out_unitp,*) '===================================='
write(out_unitp,*) '======== Expand 1D Op =============='
idq=2

write(out_unitp,*) '------------------------------------'
F11D = set_opel(2, idq, alfa=-2, indexq=1, coeff=cone) * &
       set_opel(2, idq, alfa= 3, indexq=1, coeff=cone)
CALL Export_VSCF_Op1D(F11D,'Q0',FelName)
write(out_unitp,*) 'Op1D = Q^-2 * P * Q^3 ?'
write(out_unitp,'(a,a)') '1DOp: ',FelName
CALL flush_perso(out_unitp)

CALL Expand_Op1D_TO_SumOp1D(F11D,tab_F1D)
CALL write_op(tab_F1D)

write(out_unitp,*) '------------------------------------'
F11D = set_opel(2, idq, alfa=-2, indexq=1, coeff=cone) * &
       set_opel(4, idq, alfa= 1, indexq=1, coeff=cone) * &
       set_opel(2, idq, alfa= 3, indexq=1, coeff=cone)
CALL Export_VSCF_Op1D(F11D,'Q0',FelName)
write(out_unitp,*) 'Op1D = Q^-2 * P * Q^3 ?'
write(out_unitp,'(a,a)') '1DOp: ',FelName
CALL flush_perso(out_unitp)

CALL Expand_Op1D_TO_SumOp1D(F11D,tab_F1D)
CALL write_op(tab_F1D)

write(out_unitp,*) '------------------------------------'
F11D = set_opel(4, idq, alfa= 1, indexq=1, coeff=cone) * &
       set_opel(2, idq, alfa= 3, indexq=1, coeff=cone)
CALL Export_VSCF_Op1D(F11D,'Q0',FelName)
write(out_unitp,*) 'Op1D = P * Q^3 ?'
write(out_unitp,'(a,a)') '1DOp: ',FelName
CALL flush_perso(out_unitp)

CALL Expand_Op1D_TO_SumOp1D(F11D,tab_F1D)
CALL write_op(tab_F1D)

write(out_unitp,*) '------------------------------------'
F11D = set_opel(2, idq, alfa=-2, indexq=1, coeff=cone) * &
       set_opel(4, idq, alfa= 1, indexq=1, coeff=cone)
CALL Export_VSCF_Op1D(F11D,'Q0',FelName)
write(out_unitp,*) 'Op1D = Q^-2 * P ?'
write(out_unitp,'(a,a)') '1DOp: ',FelName
CALL flush_perso(out_unitp)

CALL Expand_Op1D_TO_SumOp1D(F11D,tab_F1D)
CALL write_op(tab_F1D)

write(out_unitp,*) '------------------------------------'
F11D = set_opel(2, idq, alfa= 2, indexq=1, coeff=cone) * &
       set_opel(4, idq, alfa= 1, indexq=1, coeff=cone) * &
       set_opel(2, idq, alfa=-3, indexq=1, coeff=cone) * &
       set_opel(4, idq, alfa= 1, indexq=1, coeff=cone) * &
       set_opel(2, idq, alfa= 6, indexq=1, coeff=cone)

CALL Export_VSCF_Op1D(F11D,'Q0',FelName)
write(out_unitp,*) 'Op1D = Q^2 * P * Q^-3 * P * Q^6 ?'
write(out_unitp,'(a,a)') '1DOp: ',FelName
CALL flush_perso(out_unitp)

CALL Expand_Op1D_TO_SumOp1D(F11D,tab_F1D)
CALL write_op(tab_F1D)

DO i=1,size(list_exP2(1,:))
  write(out_unitp,*) '------------------------------------'
  write(out_unitp,*) '-F1D = x^n1 * P * x^n2 * P * x^n3 --'
  write(out_unitp,*) 'n1,n2,n3',list_exP2(:,i)
  n1 = list_exP2(1,i)
  n2 = list_exP2(2,i)
  n3 = list_exP2(3,i)

  F11D = set_opel(2, idq, alfa=n1, indexq=1, coeff=cone) * &
         set_opel(4, idq, alfa= 1, indexq=1, coeff=cone) * &
         set_opel(2, idq, alfa=n2, indexq=1, coeff=cone) * &
         set_opel(4, idq, alfa= 1, indexq=1, coeff=cone) * &
         set_opel(2, idq, alfa=n3, indexq=1, coeff=cone)

  CALL Export_VSCF_Op1D(F11D,'Q0',FelName)
  write(out_unitp,'(a,a)') '1DOp: ',FelName
  CALL flush_perso(out_unitp)
  write(out_unitp,*) ' Expand of F11D:'
  write(out_unitp,'(a,i0,a,i0,a)') '       ',1,'* x^',n1+n2+n3,' * P^2'
  write(out_unitp,'(a,i0,a,i0,a)') '       ',(n2 + 2*n3),'* x^',n1+n2+n3-1,' * P'
  write(out_unitp,'(a,i0,a,i0)') '       ',n3*(-1 + n2 + n3),'* x^',n1+n2+n3-2

  CALL Expand_Op1D_TO_SumOp1D(F11D,tab_F1D)
  CALL write_op(tab_F1D)

END DO


END SUBROUTINE test_Expand_Op1D
SUBROUTINE test_OpnD()
  USE mod_system
  USE mod_Constant

  USE mod_Coord_KEO
  USE mod_Tana_OpEl
  USE mod_Tana_Op1D
  USE mod_Tana_OpnD
  IMPLICIT NONE


      TYPE(Op1D)              :: F11D,F21D,F31D
      TYPE(Sum_OF_Op1D)       :: tab_F1D
      TYPE(OpnD)              :: FOpnD
      TYPE(OpnD), allocatable :: SumOpnD(:)

      integer :: i,idq,idf,indexq

!     - working parameters ------------------------------------------
      integer :: err_mem,memory,err_Op
      character (len=*), parameter :: name_sub='test_OpnD'


!-------------------------------------------------
! 1) 1D operator (Op1D)
!-------------------------------------------------
write(out_unitp,*) '===================================='
write(out_unitp,*) '======== nD Op ====================='
write(out_unitp,*) '------------------------------------'
idq=2
indexq=1
F11D = set_opel(2, idq, alfa=-2, indexq=indexq, coeff=cone) * &
       set_opel(4, idq, alfa= 1, indexq=indexq, coeff=cone) * &
       set_opel(2, idq, alfa= 3, indexq=indexq, coeff=cone)

indexq=2
F21D = set_opel(2, idq, alfa=-2, indexq=indexq, coeff=cone) * &
       set_opel(4, idq, alfa= 1, indexq=indexq, coeff=cone) * &
       set_opel(2, idq, alfa= 3, indexq=indexq, coeff=cone)

CALL get_F1_times_F2_to_F_nd(F11D,F21D,FOpnD)

write(out_unitp,*) '------------------------------------'
write(out_unitp,*) '------------------------------------'
write(out_unitp,*) 'FOpnD'
CALL write_op(FOpnD)
write(out_unitp,*) '------------------------------------'

CALL Expand_OpnD_TO_SumOpnD(FOpnD,SumOpnD)

write(out_unitp,*) '------------------------------------'
write(out_unitp,*) 'SumOpnD'
DO i=1,size(SumOpnD)
   CALL write_op(SumOpnD(i))
END DO
write(out_unitp,*) '------------------------------------'
write(out_unitp,*) '------------------------------------'
write(out_unitp,*) '------------------------------------'
CALL delete_op(F11D)
CALL delete_op(F21D)
CALL delete_op(F31D)
CALL delete_op(FOpnD)
CALL dealloc_NParray(SumOpnD,'SumOpnD',name_sub)

END SUBROUTINE test_OpnD
SUBROUTINE test_Sum_OpnD()
  USE mod_system
  USE mod_Constant

  USE mod_Coord_KEO
  USE mod_Tana_OpEl
  USE mod_Tana_Op1D
  USE mod_Tana_OpnD
  USE mod_Tana_Sum_OpnD
  IMPLICIT NONE

      TYPE(Sum_OpnD)              :: SumOpnD
      TYPE(Sum_OpnD)              :: ExpSumOpnD

      TYPE(Op1D)                  :: F11D,F21D,F31D

      integer :: i,idq,idf,indexq

!     - working parameters ------------------------------------------
      integer :: err_mem,memory,err_Op
      character (len=*), parameter :: name_sub='test_Sum_OpnD'


!-------------------------------------------------
! 1) 1D operator (Op1D)
!-------------------------------------------------
write(out_unitp,*) '===================================='
write(out_unitp,*) '======== Sum nD Op ================='
write(out_unitp,*) '------------------------------------'
CALL allocate_op(SumOpnD,4)

idq=2
indexq=1
F11D = set_opel(2, idq, alfa=-2, indexq=indexq, coeff=cone) * &
       set_opel(4, idq, alfa= 1, indexq=indexq, coeff=cone) * &
       set_opel(2, idq, alfa= 3, indexq=indexq, coeff=cone)

indexq=2
F21D = set_opel(2, idq, alfa=-2, indexq=indexq, coeff=cone) * &
       set_opel(4, idq, alfa= 1, indexq=indexq, coeff=cone) * &
       set_opel(2, idq, alfa= 3, indexq=indexq, coeff=cone)

CALL get_F1_times_F2_to_F_nd(F11D,F21D,SumOpnD%sum_prod_op1d(1))
SumOpnD%Cn(1) = ONE

idq=2
indexq=1
F11D = set_opel(2, idq, alfa=-2, indexq=indexq, coeff=cone) * &
       set_opel(4, idq, alfa= 1, indexq=indexq, coeff=cone) * &
       set_opel(2, idq, alfa= 3, indexq=indexq, coeff=cone)

indexq=3
F21D = set_opel(2, idq, alfa=-2, indexq=indexq, coeff=cone) * &
       set_opel(2, idq, alfa= 1, indexq=indexq, coeff=cone) * &
       set_opel(2, idq, alfa= 3, indexq=indexq, coeff=cone)

CALL get_F1_times_F2_to_F_nd(F11D,F21D,SumOpnD%sum_prod_op1d(2))
SumOpnD%Cn(2) = TWO

idq=2
indexq=1
F11D = set_opel(2, idq, alfa=-2, indexq=indexq, coeff=cone) * &
       set_opel(2, idq, alfa= 1, indexq=indexq, coeff=cone) * &
       set_opel(2, idq, alfa= 3, indexq=indexq, coeff=cone)

indexq=2
F21D = set_opel(2, idq, alfa=-2, indexq=indexq, coeff=cone) * &
       set_opel(2, idq, alfa= 1, indexq=indexq, coeff=cone) * &
       set_opel(2, idq, alfa= 3, indexq=indexq, coeff=cone)

CALL get_F1_times_F2_to_F_nd(F11D,F21D,SumOpnD%sum_prod_op1d(3))
SumOpnD%Cn(3) = THREE


idq=2
indexq=2
F11D = set_opel(2, idq, alfa=-2, indexq=indexq, coeff=cone) * &
       set_opel(4, idq, alfa= 1, indexq=indexq, coeff=cone) * &
       set_opel(2, idq, alfa= 3, indexq=indexq, coeff=cone) * &
       set_opel(4, idq, alfa= 1, indexq=indexq, coeff=cone) * &
       set_opel(2, idq, alfa= 3, indexq=indexq, coeff=cone)

indexq=3
F21D = set_opel(2, idq, alfa=-2, indexq=indexq, coeff=cone) * &
       set_opel(2, idq, alfa= 3, indexq=indexq, coeff=cone)

CALL get_F1_times_F2_to_F_nd(F11D,F21D,SumOpnD%sum_prod_op1d(4))
SumOpnD%Cn(4) = FOUR

write(out_unitp,*) '------------------------------------'
write(out_unitp,*) '------------------------------------'
write(out_unitp,*) 'SumOpnD'
CALL write_op(SumOpnD)
write(out_unitp,*) '------------------------------------'

CALL Expand_Sum_OpnD_TO_Sum_OpnD(SumOpnD,ExpSumOpnD)

write(out_unitp,*) '------------------------------------'
write(out_unitp,*) 'ExpSumOpnD'
write(out_unitp,*) '------------------------------------'
CALL write_op(ExpSumOpnD)
write(out_unitp,*) '------------------------------------'
write(out_unitp,*) '------------------------------------'
write(out_unitp,*) '------------------------------------'

CALL delete_op(F11D)
CALL delete_op(F21D)
CALL delete_op(F31D)

END SUBROUTINE test_Sum_OpnD
SUBROUTINE test_String_TO_Sum_OpnD()
  USE mod_system
  USE mod_Constant

  USE mod_Coord_KEO
  USE mod_Tana_OpEl
  USE mod_Tana_Op1D
  USE mod_Tana_OpnD
  USE mod_Tana_Sum_OpnD
  USE mod_Tana_write_mctdh
  IMPLICIT NONE

      TYPE(Sum_OpnD)              :: SumOpnD
      TYPE(Sum_OpnD)              :: ExpSumOpnD

      TYPE(Op1D)                  :: F1D

      TYPE(OpEl)                  :: F1el
      TYPE(OpnD)                  :: F1nD
      TYPE(sum_opnd)              :: KEO_MCTDH

      character (len = :), allocatable     :: String_OpnD

      integer :: i,idq,idf,indexq

      !- working parameters ------------------------------------------
      integer :: err_mem,memory,err_Op
      character (len=*), parameter :: name_sub='test_String_TO_Sum_OpnD'

!  CALL StringMCTDH_TO_OpEl(F1el,'qs^-1',indexq=3)
!  CALL write_op(F1el,header=.TRUE.)
!  write(out_unitp,*) '===================================='
!
!  CALL StringMCTDH_TO_Op1d(F1D,'sin^-1', indexq=1)
!  CALL write_op(F1D,header=.TRUE.)
!
!  CALL StringMCTDH_TO_Op1d(F1D,'q^-1*dq*q', indexq=2)
!  CALL write_op(F1D,header=.TRUE.)
!
!  write(out_unitp,*) '===================================='
!  write(out_unitp,*) '       q^-1*dq*q^2*dq*q^-1'
!  CALL StringMCTDH_TO_Op1d(F1D,'q^-1*dq*q^2*dq*q^-1', indexq=2)
!  CALL write_op(F1D,header=.TRUE.)
!  write(out_unitp,*) '===================================='
!
!  write(out_unitp,*) '===================================='
!  String_OpnD = '-1.4437526746680724d-005 |1   q^-1 |3   q*qs^-1 |4   q^-1 |5   dq*q*qs |6   sin*dq'
!  write(out_unitp,*) String_OpnD
!
!  CALL StringMCTDH_TO_Opnd(F1nD,String_OpnD,nb_act=6)
!  CALL write_op(F1nD,header=.TRUE.)
!  write(out_unitp,*) '===================================='
!
!  CALL delete_op(F1D)
!  CALL delete_op(F1nD)

  CALL read_keo_mctdh_form(nb_act=6,keo=KEO_MCTDH,io=in_unitp)
  CALL write_op(KEO_MCTDH,header=.TRUE.)


END SUBROUTINE test_String_TO_Sum_OpnD
