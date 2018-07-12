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

!================================================================
!     NLO
!================================================================
      SUBROUTINE sub_NLO(para_Dip,print_Op,para_H,nb_ana,para_intensity)

      USE mod_system
      USE mod_Tnum
      USE mod_basis
      USE mod_Op
      USE mod_constant
      USE mod_analysis
      IMPLICIT NONE


!----- Operator variables ----------------------------------------------
      TYPE (param_Op)  :: para_Dip(3),para_H
      logical          :: print_Op
      integer          :: nb_ana

!----- variables pour la namelist analyse ------------------------------
      TYPE (param_intensity) :: para_intensity



!---- variable for the Z-matrix ----------------------------------------
      TYPE (zmatrix), pointer    :: mole
      TYPE (Tnum), pointer       :: para_Tnum


!----- working variables -----------------------------
      real (kind=Rkind), allocatable ::    Mat_Aif(:,:)
      integer       ::    i,k



!----- for debuging --------------------------------------------------
      integer   :: err
      logical, parameter :: debug=.FALSE.
!      logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'sub_NLO'
!-----------------------------------------------------------
      mole       => para_H%mole
      para_Tnum  => para_H%para_Tnum

      write(out_unitp,*) 'BEGINNING ',name_sub
      IF (debug) THEN
        write(out_unitp,*) 'nb_ana',nb_ana
        write(out_unitp,*) 'Rvp',shape(para_H%Rvp)
!       CALL Write_Mat(para_H%Rvp,out_unitp,5)
        write(out_unitp,*)
      END IF
!-----------------------------------------------------------


      IF (nb_ana > 0) THEN
        CALL alloc_NParray(Mat_Aif,(/ nb_ana,nb_ana /),'Mat_Aif',name_sub)
        Mat_Aif(:,:) = ZERO
      ELSE
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' nb_ana <=0',nb_ana
        STOP
      END IF
      CALL flush_perso(out_unitp)



!     -------------------------------------------
      write(out_unitp,*) ' ene (ua): ',nb_ana
      DO i=1,nb_ana
        write(out_unitp,*) i,para_H%Rdiag(i)
      END DO
      write(out_unitp,*) ' END ene',nb_ana
      write(out_unitp,*)

      write(out_unitp,*) '==================================================='
      write(out_unitp,*) '==================================================='
      write(out_unitp,*) ' Calculation of "Mat_Aif(:,:)": '

      Mat_Aif(:,:) = ZERO
      DO k=1,3
        Mat_Aif(:,:) = Mat_Aif(:,:) + para_Dip(k)%Rmat(1:nb_ana,1:nb_ana)**2
      END DO
      CALL Write_Mat(Mat_Aif,out_unitp,5,Rformat='e30.23')
      write(out_unitp,*) '==================================================='
      write(out_unitp,*) '==================================================='
      CALL flush_perso(out_unitp)


      CALL dealloc_NParray(Mat_Aif,'Mat_Aif',name_sub)


!----------------------------------------------------------
      IF (debug) THEN
      END IF
      write(out_unitp,*) 'END ',name_sub
!----------------------------------------------------------

      end subroutine sub_NLO
