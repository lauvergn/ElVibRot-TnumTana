!===========================================================================
!===========================================================================
!This file is part of ElVibRot.
!
!    ElVibRot is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ElVibRot is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ElVibRot.  If not, see <http://www.gnu.org/licenses/>.
!
!    Copyright 2015  David Lauvergnat
!      with contributions of Mamadou Ndong, Josep Maria Luis
!
!    ElVibRot includes:
!        - Tnum-Tana under the GNU LGPL3 license
!        - Somme subroutines of John Burkardt under GNU LGPL license
!             http://people.sc.fsu.edu/~jburkardt/
!        - Somme subroutines of SHTOOLS written by Mark A. Wieczorek under BSD license
!             http://shtools.ipgp.fr
!===========================================================================
!===========================================================================

      MODULE mod_Optimization
      USE mod_system
      USE mod_SimulatedAnnealing
      USE mod_BFGS

      IMPLICIT NONE


        TYPE param_Optimization

        character (len=Name_len) :: Optimization_method = 'SimulatedAnnealing'
        character (len=Name_len) :: Optimization_param  = 'coordbasis'
        logical                  :: FinalEnergy         = .TRUE.
        logical                  :: Grad                = .FALSE.
        logical                  :: Freq                = .FALSE.
        logical                  :: ReadRange           = .FALSE.

        TYPE (param_SimulatedAnnealing) :: para_SimulatedAnnealing
        TYPE (param_BFGS)               :: para_BFGS


        END TYPE param_Optimization

      CONTAINS

      SUBROUTINE Read_param_Optimization(para_Optimization)
      TYPE (param_Optimization), intent(inout) :: para_Optimization

        character (len=Name_len) :: Optimization_method = 'SimulatedAnnealing'
        character (len=Name_len) :: Optimization_param  = 'coordbasis'
        logical                  :: FinalEnergy         = .TRUE.
        logical                  :: Freq                = .FALSE.
        logical                  :: Grad                = .FALSE.
        logical                  :: ReadRange           = .FALSE.

        integer :: err_io
        NAMELIST /Optimization/ Optimization_method,Optimization_param, &
                                FinalEnergy,Freq,Grad,ReadRange

        Optimization_method = 'SimulatedAnnealing'
        Optimization_param  = 'coordbasis'
        FinalEnergy         = .TRUE. ! redo the optimal energy
        Freq                = .FALSE. ! Calculation of the frequencies (geometry optimization)
        Grad                = .FALSE. ! Calculation of the gradiant
        ReadRange           = .FALSE.

        read(in_unitp,Optimization,IOSTAT=err_io)
        IF (err_io /= 0) THEN
          write(out_unitp,*) ' WARNING in Read_param_Optimization'
          write(out_unitp,*) '  while reading the "Optimization" namelist'
          write(out_unitp,*) ' end of file or end of record'
          write(out_unitp,*) ' Check your data !!'
          STOP
        END IF
        IF (print_level > 1) write(out_unitp,Optimization)

        CALL string_uppercase_TO_lowercase(Optimization_method)
        para_Optimization%Optimization_method    = Optimization_method
        CALL string_uppercase_TO_lowercase(Optimization_param)
        para_Optimization%Optimization_param     = Optimization_param

        para_FOR_optimization%Optimization_param = Optimization_param
        !IF (Optimization_param == 'geometry') FinalEnergy = .FALSE.
        IF (Optimization_param /= 'geometry') Freq = .FALSE.

        para_Optimization%FinalEnergy            = FinalEnergy
        para_Optimization%Freq                   = Freq
        para_Optimization%Grad                   = Grad
        para_Optimization%ReadRange              = ReadRange

        SELECT CASE (para_Optimization%Optimization_method)
        CASE ('simulatedannealing','sa')
          CALL Read_param_SimulatedAnnealing(                           &
                              para_Optimization%para_SimulatedAnnealing)
          !CALL Read_param_BFGS(para_Optimization%para_BFGS)

        CASE ('bfgs')
          CALL Read_param_BFGS(para_Optimization%para_BFGS)

        CASE DEFAULT
          write(out_unitp,*) 'ERROR in Read_param_Optimization'
          write(out_unitp,*) 'Wrong Optimization_method: ',             &
                                    para_Optimization%Optimization_method
          write(out_unitp,*) ' Check your data!!'
          STOP
        END SELECT

      END SUBROUTINE Read_param_Optimization
      SUBROUTINE Write_param_Optimization(para_Optimization)
      TYPE (param_Optimization), intent(in) :: para_Optimization

      write(out_unitp,*) 'WRITE param_Optimization'
      write(out_unitp,*)
      write(out_unitp,*) 'Optimization_method        ',                 &
                                    para_Optimization%Optimization_method
      write(out_unitp,*) 'Optimization_param         ',                 &
                                    para_Optimization%Optimization_param
      write(out_unitp,*) 'ReadRange                  ',                 &
                                    para_Optimization%ReadRange
      write(out_unitp,*) 'FinalEnergy                ',                 &
                                    para_Optimization%FinalEnergy
      write(out_unitp,*) 'Grad, Freq                 ',                 &
                           para_Optimization%Grad,para_Optimization%Freq

      SELECT CASE (para_Optimization%Optimization_method)
      CASE ('simulatedannealing','sa')
        CALL Write_param_SimulatedAnnealing(                            &
                              para_Optimization%para_SimulatedAnnealing)
      CASE ('bfgs')
        CALL Write_param_BFGS(para_Optimization%para_BFGS)
      CASE DEFAULT
        write(out_unitp,*) 'WARNING in Write_param_param_Optimization'
        write(out_unitp,*) 'No parameter for this Optimization_method', &
                                     para_Optimization%Optimization_method
      END SELECT
      write(out_unitp,*) 'END WRITE param_Optimization'

      END SUBROUTINE Write_param_Optimization

      END MODULE mod_Optimization

