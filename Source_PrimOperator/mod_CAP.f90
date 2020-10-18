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
   MODULE mod_CAP
   USE mod_system
   IMPLICIT NONE
   PRIVATE

     TYPE CAP_t

        integer                       :: Type_CAP             = 1   ! 1:  A*B * x**n_exp
                                                                    ! 2: Woods-Saxon: 2*A/(1+Exp(-x))
        character (len=Name_longlen)  :: Name_Cap             = ""  ! 1 => "x^n" or "xn"
                                                                    ! 2 => "Woods-Saxon" or "WS"
        integer                       :: n_exp                = 2
        real (kind=Rkind)             :: A                    = ONE
        real (kind=Rkind)             :: B                    = ONE

                                                              ! x = (Q-Q0)/LQ
        real (kind=Rkind)             :: Q0                   = ZERO
        real (kind=Rkind)             :: LQ                   = ONE
        integer                       :: ind_Q                = 1       ! index of the coordinate (active order)

        integer                       :: iOp                  = 0       ! index of the Operator

      CONTAINS
        PROCEDURE, PRIVATE, PASS(CAP1) :: CAP2_TO_CAP1
        GENERIC,   PUBLIC  :: assignment(=) => CAP2_TO_CAP1
      END TYPE CAP_t

    PUBLIC :: CAP_t, Read_CAP, write_CAP, dealloc_CAP, calc_CAP

  CONTAINS

  SUBROUTINE write_CAP(CAP)
  IMPLICIT NONE

      TYPE (CAP_t) :: CAP


      write(out_unitp,*) ' BEGINNING write_CAP'

      IF (CAP%Type_CAP > 0) THEN
        write(out_unitp,*) 'Type_CAP > 0'
        write(out_unitp,*) 'CAP(Q)'
        write(out_unitp,*) ' ^'
        write(out_unitp,*) ' |                   /'
        write(out_unitp,*) ' |                  /'
        write(out_unitp,*) ' |                 /'
        write(out_unitp,*) ' |----------------/..........> Q'
        write(out_unitp,*) '                 Q0'
        write(out_unitp,*) ' CAP(Q)=0 when Q<= Q0'
        write(out_unitp,*) ' x=+(Q-0)'
      ELSE
        write(out_unitp,*) 'Type_CAP < 0'
        write(out_unitp,*) 'CAP(Q)'
        write(out_unitp,*) ' ^'
        write(out_unitp,*) ' |             \'
        write(out_unitp,*) ' |              \'
        write(out_unitp,*) ' |               \'
        write(out_unitp,*) ' |................\----------> Q'
        write(out_unitp,*) '                  Q0'
        write(out_unitp,*) ' CAP(Q)=0 when Q>= Q0'
        write(out_unitp,*) ' x=-(Q-0)'

      END IF

      write(out_unitp,*) 'Name_CAP',CAP%Name_CAP
      write(out_unitp,*) 'Type_CAP',CAP%Type_CAP
      SELECT CASE (CAP%Type_CAP)
      CASE(1,-1)
        write(out_unitp,*) 'CAP(Q)=A * B*x^',int_TO_char(CAP%n_exp)
        write(out_unitp,*) '   with B=(',int_TO_char(CAP%n_exp),'+1)/2'
        write(out_unitp,*) ' ref: A. Vibok and G. G. Balint-Kurti, J. Phys. Chem. (1992) 96, 8712-8719'
      CASE(2,-2)
        write(out_unitp,*) 'CAP(Q)=A * 2/(1-exp(-x))'
        write(out_unitp,*) ' ref: T. Seideman and W. H. Miller, J. Chem. Phys. (1992) 96, 4412'
        write(out_unitp,*) ' ref: R. D. Woods and D. S. Saxon, Phys. Rev. (1954) 95, 577'
      CASE(3,-3)
        write(out_unitp,*) 'CAP(Q)=A * B*exp(-2/x)'
        write(out_unitp,*) '   with B=13.22'
        write(out_unitp,*) ' ref: A. Vibok and G. G. Balint-Kurti, J. Phys. Chem. (1992), 96, 8712-8719'
      CASE default
        STOP 'ERROR in write_CAP: no default'
      END SELECT

      write(out_unitp,*) 'n_exp   ',CAP%n_exp
      write(out_unitp,*) 'A       ',CAP%A
      write(out_unitp,*) 'B       ',CAP%B

      write(out_unitp,*) 'ind_Q (active order)   ',CAP%ind_Q
      write(out_unitp,*) 'Q0      ',CAP%Q0
      write(out_unitp,*) 'LQ      ',CAP%LQ

      write(out_unitp,*) 'iOp     ',CAP%iOp

    write(out_unitp,*) ' END write_CAP'
    CALL flush_perso(out_unitp)
  END SUBROUTINE write_CAP
  SUBROUTINE CAP2_TO_CAP1(CAP1,CAP2)
  IMPLICIT NONE
      CLASS (CAP_t), intent(inout) :: CAP1
      TYPE (CAP_t),  intent(in)    :: CAP2

      !write(out_unitp,*) ' BEGINNING CAP2_TO_CAP1'

      CAP1%Type_CAP            = CAP2%Type_CAP
      CAP1%Name_CAP            = CAP2%Name_CAP
      CAP1%n_exp               = CAP2%n_exp
      CAP1%A                   = CAP2%A
      CAP1%B                   = CAP2%B

      CAP1%Q0                  = CAP2%Q0
      CAP1%LQ                  = CAP2%LQ
      CAP1%ind_Q               = CAP2%ind_Q

      CAP1%iOp                 = CAP2%iOp

     !write(out_unitp,*) ' END CAP2_TO_CAP1'
     !CALL flush_perso(out_unitp)
  END SUBROUTINE CAP2_TO_CAP1
  SUBROUTINE Read_CAP(CAP_in)
  IMPLICIT NONE
      CLASS (CAP_t),    intent(inout) :: CAP_in

      character (len=Name_longlen)  :: Name_Cap
      integer             :: Type_CAP,n_exp,ind_Q
      real(kind=Rkind)    :: A,Q0,LQ
      integer             :: err_read

      namelist / CAP / Type_CAP,Name_CAP,n_exp,A,Q0,LQ,ind_Q

      Type_CAP             = 0
      Name_Cap             = ""
      n_exp                = 2
      A                    = ONE
      Q0                   = ZERO
      LQ                   = ONE
      ind_Q                = -1
      read(in_unitp,CAP,IOSTAT=err_read)
      IF (err_read < 0) THEN
        write(out_unitp,*) ' ERROR in Read_CAP'
        write(out_unitp,*) ' End-of-file or End-of-record'
        write(out_unitp,*) ' The namelist "CAP" is probably absent'
        write(out_unitp,*) ' check your data!'
        write(out_unitp,*) ' ERROR in Read_CAP'
        STOP ' ERROR in Read_CAP'
      ELSE IF (err_read > 0) THEN
        write(out_unitp,*) ' ERROR in Read_CAP'
        write(out_unitp,*) ' Some parameter name of the namelist "CAP" are probaly wrong'
        write(out_unitp,*) ' check your data!'
        write(out_unitp,CAP)
        write(out_unitp,*) ' ERROR in Read_CAP'
        STOP ' ERROR in Read_CAP'
      END IF
      IF (print_level > 1) write(out_unitp,CAP)

      CALL Init_CAP(CAP_in,Type_CAP,Name_Cap,n_exp,A,Q0,LQ,ind_Q)

      CALL Write_CAP(CAP_in)

  END SUBROUTINE Read_CAP
  SUBROUTINE Init_CAP(CAP,Type_CAP,Name_Cap,n_exp,A,Q0,LQ,ind_Q)
  IMPLICIT NONE
      CLASS (CAP_t),                intent(inout) :: CAP
      character (len=Name_longlen), intent(in)    :: Name_Cap
      integer,                      intent(in)    :: Type_CAP,n_exp,ind_Q
      real(kind=Rkind),             intent(in)    :: A,Q0,LQ

      if ( len(adjustl(trim(Name_Cap))) == 0 .AND. Type_CAP == 0 ) then
        CAP%Type_CAP = 1
      else if ( Type_CAP == 0 ) then
        CAP%Name_Cap = Name_Cap
        CALL string_uppercase_TO_lowercase(CAP%Name_Cap)
        SELECT CASE (CAP%Name_Cap)
        CASE("xn","x^n","-xn","-x^n","+xn","+x^n")
          CAP%Type_CAP = 1
        CASE("ws","woods-saxon","-ws","-woods-saxon","+ws","+woods-saxon")
          CAP%Type_CAP = 2
        CASE("exp","-exp","+exp")
          CAP%Type_CAP = 3
        CASE default
          STOP 'ERROR in Init_CAP: no Name_Cap default'
        END SELECT
        if (CAP%Name_Cap(1:1) == "-") CAP%Type_CAP = -CAP%Type_CAP
      end if

      CAP = CAP_t(Type_CAP=Type_CAP,n_exp=n_exp,A=A,Q0=Q0,LQ=LQ,ind_Q=ind_Q)

      SELECT CASE (CAP%Type_CAP)
      CASE(-1,1) ! as function of n_exp, B= 1 3/2 2 5/2 ...
        CAP%B = ONE + (n_exp-1)*HALF
        CAP%Name_Cap = "x^n"
      CASE(-2,2)
        CAP%B = ZERO
        CAP%Name_Cap = "Woods-Saxon"
      CASE(-3,3)
        CAP%B = 13.22_Rkind
        CAP%Name_Cap = "exp"
      CASE default
        STOP 'ERROR in Init_CAP: no Type_CAP default'
      END SELECT

      if ( CAP%Type_CAP > 0 ) then
        CAP%Name_Cap = '+' // trim(CAP%Name_Cap)
      else
        CAP%Name_Cap = '-' // trim(CAP%Name_Cap)
      end if

      !CALL Write_CAP(CAP)

  END SUBROUTINE Init_CAP
  SUBROUTINE dealloc_CAP(CAP)
  IMPLICIT NONE
      CLASS (CAP_t), intent(inout) :: CAP

      !write(out_unitp,*) ' BEGINNING dealloc_CAP'

        CAP%Type_CAP             = 1       ! 1:  A*(B*x)**n_exp
        CAP%Name_Cap             = "x^n"   ! 1:  A*(B*x)**n_exp
        CAP%n_exp                = 2
        CAP%A                    = ONE
        CAP%B                    = THREE/TWO

                                           ! x = (Q-Q0)/LQ
        CAP%Q0                   = ZERO
        CAP%LQ                   = ONE
        CAP%ind_Q                = 1

        CAP%iOp                  = 0
     !write(out_unitp,*) ' END dealloc_CAP'
     !CALL flush_perso(out_unitp)
  END SUBROUTINE dealloc_CAP

  FUNCTION calc_CAP(CAP,Q)
  IMPLICIT NONE
      real (kind=Rkind)               :: calc_CAP
      CLASS (CAP_t),    intent(in)    :: CAP
      real(kind=Rkind), intent(in)    :: Q(:)

      real(kind=Rkind)    :: x

      calc_CAP = ZERO

      SELECT CASE (CAP%Type_CAP)
      CASE(1) ! x^n
        IF ( Q(CAP%ind_Q) > CAP%Q0 ) THEN
          x = (Q(CAP%ind_Q)-CAP%Q0)/CAP%LQ
          calc_CAP = CAP%A * CAP%B * x**CAP%n_exp
        END IF
      CASE(-1) ! x^n
        IF ( Q(CAP%ind_Q) < CAP%Q0 ) THEN
          x = -(Q(CAP%ind_Q)-CAP%Q0)/CAP%LQ
          calc_CAP = CAP%A * CAP%B * x**CAP%n_exp
        END IF

      CASE(2) ! WS
        x = (Q(CAP%ind_Q)-CAP%Q0)/CAP%LQ
        calc_CAP = TWO*CAP%A/(ONE+Exp(-x))
      CASE(-2)  ! WS
        x = -(Q(CAP%ind_Q)-CAP%Q0)/CAP%LQ
        calc_CAP = TWO*CAP%A/(ONE+Exp(-x))

      CASE(3)  ! exp
        IF ( Q(CAP%ind_Q) > CAP%Q0 ) THEN
          x = (Q(CAP%ind_Q)-CAP%Q0)/CAP%LQ
          calc_CAP = CAP%A * CAP%B*exp(-TWO/x)
        END IF
      CASE(-3) ! exp
        IF ( Q(CAP%ind_Q) < CAP%Q0 ) THEN
          x = -(Q(CAP%ind_Q)-CAP%Q0)/CAP%LQ
          calc_CAP = CAP%A * CAP%B*exp(-TWO/x)
        END IF

      CASE default
        STOP 'ERROR in calc_CAP: no default'
      END SELECT

     !write(out_unitp,*) ' END calc_CAP'
     !CALL flush_perso(out_unitp)
  END FUNCTION calc_CAP

  END MODULE mod_CAP
