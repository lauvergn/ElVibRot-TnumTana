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
   MODULE mod_HStep
   USE mod_system
   IMPLICIT NONE
   PRIVATE

     TYPE HStep_t

        integer           :: Type_HStep          = 1  ! 1:  A*B * x**n_exp

        real (kind=Rkind) :: Q0                   = ZERO
        integer           :: ind_Q                = 1       ! index of the coordinate

        integer           :: iOp                  = 0       ! index of the Operator

      CONTAINS
        PROCEDURE, PRIVATE, PASS(HStep1) :: HStep2_TO_HStep1
        GENERIC,   PUBLIC  :: assignment(=) => HStep2_TO_HStep1
      END TYPE HStep_t

    PUBLIC :: HStep_t, Read_HStep, write_HStep, dealloc_HStep, calc_HStep

  CONTAINS

  SUBROUTINE write_HStep(HStep)
  IMPLICIT NONE

      TYPE (HStep_t) :: HStep

      write(out_unitp,*) ' BEGINNING write_HStep'

      IF (HStep%Type_HStep > 0) THEN
        write(out_unitp,*) 'Type_HStep > 0'
        write(out_unitp,*) 'HStep(Q)'
        write(out_unitp,*) ' ^'
        write(out_unitp,*) ' 1                ------------------'
        write(out_unitp,*) ' |                |'
        write(out_unitp,*) ' |                |'
        write(out_unitp,*) ' |                |'
        write(out_unitp,*) ' |----------------|..................> Q'
        write(out_unitp,*) '                 Q0'
        write(out_unitp,*) ' HStep(Q)=0 when Q<= Q0'
        write(out_unitp,*) ' HStep(Q)=1 when Q> Q0'

      ELSE
        write(out_unitp,*) 'Type_HStep < 0'
        write(out_unitp,*) 'HStep(Q)'
        write(out_unitp,*) ' ^'
        write(out_unitp,*) ' 1-----------------'
        write(out_unitp,*) ' |                |'
        write(out_unitp,*) ' |                |'
        write(out_unitp,*) ' |                |'
        write(out_unitp,*) ' |................|-----------------> Q'
        write(out_unitp,*) '                  Q0'
        write(out_unitp,*) ' HStep(Q)=0 when Q>= Q0'
        write(out_unitp,*) ' HStep(Q)=1 when Q< Q0'
      END IF

      write(out_unitp,*) 'Type_HStep',HStep%Type_HStep
      write(out_unitp,*) 'ind_Q   ',HStep%ind_Q
      write(out_unitp,*) 'Q0      ',HStep%Q0

      write(out_unitp,*) 'iOp     ',HStep%iOp
      write(out_unitp,*) ' END write_HStep'
      CALL flush_perso(out_unitp)
  END SUBROUTINE write_HStep
  SUBROUTINE HStep2_TO_HStep1(HStep1,HStep2)
  IMPLICIT NONE
      CLASS (HStep_t), intent(inout) :: HStep1
      TYPE (HStep_t),  intent(in)    :: HStep2

      !write(out_unitp,*) ' BEGINNING HStep2_TO_HStep1'

      HStep1%Type_HStep         = HStep2%Type_HStep

      HStep1%Q0                  = HStep2%Q0
      HStep1%ind_Q               = HStep2%ind_Q

      HStep1%iOp                 = HStep2%iOp

     !write(out_unitp,*) ' END HStep2_TO_HStep1'
     !CALL flush_perso(out_unitp)
  END SUBROUTINE HStep2_TO_HStep1
  SUBROUTINE Read_HStep(HStep_in)
  IMPLICIT NONE
      CLASS (HStep_t),    intent(inout) :: HStep_in


      integer             :: Type_HStep,n_exp,ind_Q
      real(kind=Rkind)    :: A,Q0,LQ
      integer             :: err_read

      namelist / HStep / Type_HStep,Q0,ind_Q

      Type_HStep          = 1       ! 1:  A*B*x**n_exp
      Q0                   = ZERO
      ind_Q                = -1
      read(in_unitp,HStep,IOSTAT=err_read)
      IF (err_read < 0) THEN
        write(out_unitp,*) ' ERROR in Read_HStep'
        write(out_unitp,*) ' End-of-file or End-of-record'
        write(out_unitp,*) ' The namelist "HStep" is probably absent'
        write(out_unitp,*) ' check your data!'
        write(out_unitp,*) ' ERROR in Read_HStep'
        STOP ' ERROR in Read_HStep'
      ELSE IF (err_read > 0) THEN
        write(out_unitp,*) ' ERROR in Read_HStep'
        write(out_unitp,*) ' Some parameter name of the namelist "HStep" are probaly wrong'
        write(out_unitp,*) ' check your data!'
        write(out_unitp,HStep)
        write(out_unitp,*) ' ERROR in Read_HStep'
        STOP ' ERROR in Read_HStep'
      END IF
      IF (print_level > 1) write(out_unitp,HStep)

      CALL Init_HStep(HStep_in,Type_HStep,Q0,ind_Q)

      CALL Write_HStep(HStep_in)

  END SUBROUTINE Read_HStep
  SUBROUTINE Init_HStep(HStep,Type_HStep,Q0,ind_Q)
  IMPLICIT NONE
      CLASS (HStep_t), intent(inout) :: HStep
      integer,          intent(in)    :: Type_HStep,ind_Q
      real(kind=Rkind), intent(in)    :: Q0


      HStep = HStep_t(Type_HStep=Type_HStep,Q0=Q0,ind_Q=ind_Q)

      !CALL Write_HStep(HStep)

  END SUBROUTINE Init_HStep
  SUBROUTINE dealloc_HStep(HStep)
  IMPLICIT NONE
      CLASS (HStep_t), intent(inout) :: HStep

      !write(out_unitp,*) ' BEGINNING dealloc_HStep'

        HStep%Type_HStep          = 1       ! 1:  A*(B*x)**n_exp

        HStep%Q0                   = ZERO
        HStep%ind_Q                = 1

        HStep%iOp                  = 0
     !write(out_unitp,*) ' END dealloc_HStep'
     !CALL flush_perso(out_unitp)
  END SUBROUTINE dealloc_HStep

  FUNCTION calc_HStep(HStep,Q)
  IMPLICIT NONE
      real (kind=Rkind)               :: calc_HStep
      CLASS (HStep_t),  intent(in)    :: HStep
      real(kind=Rkind), intent(in)    :: Q(:)

      integer          :: option = 1
      real(kind=Rkind) :: Scale  = FIVE

      real(kind=Rkind) :: x

      SELECT CASE (option)
      CASE (0)
        IF ((Q(HStep%ind_Q) <= HStep%Q0 .AND. HStep%Type_HStep > 0) .OR.                             &
            (Q(HStep%ind_Q) >= HStep%Q0 .AND. HStep%Type_HStep < 0) ) THEN
          calc_HStep = ZERO
        ELSE
          calc_HStep = ONE
        END IF
      CASE (1)
        x = Scale*(Q(HStep%ind_Q)-HStep%Q0)
        IF (HStep%Type_HStep > 0) THEN
          calc_HStep = HALF*(ONE+tanh(x))
        else
          calc_HStep = HALF*(ONE+tanh(-x))
        END IF
      END SELECT



     !write(out_unitp,*) ' END calc_HStep'
     !CALL flush_perso(out_unitp)
  END FUNCTION calc_HStep

  END MODULE mod_HStep
