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
      MODULE mod_HyperSpheTransfo
      use mod_system
      USE mod_dnSVM
      IMPLICIT NONE

      PRIVATE
      !!@description: TODO
      !!@param: TODO
      TYPE Type_HyperSpheTransfo
        integer          :: nb_HyperSphe      = 0
        integer, pointer :: list_HyperSphe(:) => null()
      END TYPE Type_HyperSpheTransfo

      PUBLIC :: Type_HyperSpheTransfo, Read_HyperSpheTransfo, Calc_HyperSpheTransfo

      CONTAINS

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE Read_HyperSpheTransfo(HyperSpheTransfo,nb_Qin)

      TYPE (Type_HyperSpheTransfo), intent(inout) :: HyperSpheTransfo
      integer, intent(in) :: nb_Qin

      integer :: i,ii,it,err
      integer :: list_HyperSphe(nb_Qin)

      character (len=*), parameter :: name_sub='Read_HyperSpheTransfo'


      read(in_unitp,*,IOSTAT=err) list_HyperSphe(:)
      IF (err /= 0) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) '  while reading "list_HyperSphe"'
         write(out_unitp,*) ' end of file or end of record'
         write(out_unitp,*) ' Check your data !!'
         STOP
      END IF

      HyperSpheTransfo%nb_HyperSphe = count(list_HyperSphe /= 0)

      IF (HyperSpheTransfo%nb_HyperSphe < 2) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' nb_HyperSphe is smaller than 2 !!'
        write(out_unitp,*) ' Check your data !!'
        write(out_unitp,*) 'list_HyperSphe: ',HyperSpheTransfo%list_HyperSphe(:)
        STOP
      END IF

      CALL alloc_array(HyperSpheTransfo%list_HyperSphe,                 &
                                   (/HyperSpheTransfo%nb_HyperSphe/),   &
                      "HyperSpheTransfo%list_HyperSphe",name_sub)
      HyperSpheTransfo%list_HyperSphe(:) = 0

      IF (count(list_HyperSphe == 1) == HyperSpheTransfo%nb_HyperSphe) THEN

        ii = 0
        DO i=1,nb_Qin
          IF (list_HyperSphe(i) /= 0) THEN
            ii = ii + 1
            HyperSpheTransfo%list_HyperSphe(ii) = i
          END IF
        END DO

      ELSE ! enables to chose the order R,th or th,R

      DO i=1,HyperSpheTransfo%nb_HyperSphe ! check is the list includes 1 2 3 ... (only once)
        IF (count(list_HyperSphe == i) /= 1) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' Wrong list: '
          write(out_unitp,*) 'list_HyperSphe: ',list_HyperSphe
          write(out_unitp,*) ' Check your data !!'
          STOP
        END IF
      END DO

        DO i=1,nb_Qin
          IF (list_HyperSphe(i) /= 0) THEN
            HyperSpheTransfo%list_HyperSphe(list_HyperSphe(i)) = i
          END IF
        END DO


      END IF


      write(out_unitp,*) 'nb_HyperSphe: ',                              &
                 HyperSpheTransfo%nb_HyperSphe
      write(out_unitp,*) 'list_HyperSphe: ',HyperSpheTransfo%list_HyperSphe(:)

      END SUBROUTINE Read_HyperSpheTransfo

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE calc_HyperSpheTransfo(dnQin,dnQout,HyperSpheTransfo,   &
                                       nderiv,inTOout)

        TYPE (Type_dnVec), intent(inout)         :: dnQin,dnQout
        TYPE (Type_HyperSpheTransfo), intent(in) :: HyperSpheTransfo
        integer, intent(in)                      :: nderiv
        logical, intent(in)                      :: inTOout


        TYPE (Type_dnS) :: dnRho,dnTheta,dnCosTheta,dnSinTheta
        TYPE (Type_dnS) :: dnX,dnY
        integer :: i1,i2
!----- for debuging ----------------------------------
       character (len=*),parameter :: name_sub='calc_HyperSpheTransfo'
       logical, parameter :: debug=.FALSE.
!       logical, parameter :: debug=.TRUE.
!----- for debuging ----------------------------------


!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'dnQin'
        CALL Write_dnSVM(dnQin,nderiv)
      END IF
!---------------------------------------------------------------------

      IF (HyperSpheTransfo%nb_HyperSphe > 2) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) '  This subroutine works nb_HyperSphe = 2'
        write(out_unitp,*) '  nb_HyperSphe: ',HyperSpheTransfo%nb_HyperSphe
        write(out_unitp,*) ' Check your data !!'
        STOP
      END IF

      CALL check_alloc_dnVec(dnQin,'dnQin',name_sub)
      CALL check_alloc_dnVec(dnQout,'dnQout',name_sub)

      CALL alloc_dnSVM(dnX,dnQin%nb_var_deriv,nderiv)
      CALL alloc_dnSVM(dnY,dnQin%nb_var_deriv,nderiv)
      CALL alloc_dnSVM(dnRho,dnQin%nb_var_deriv,nderiv)
      CALL alloc_dnSVM(dnTheta,dnQin%nb_var_deriv,nderiv)
      CALL alloc_dnSVM(dnCosTheta,dnQin%nb_var_deriv,nderiv)
      CALL alloc_dnSVM(dnSinTheta,dnQin%nb_var_deriv,nderiv)

      IF (inTOout) THEN
        CALL sub_dnVec1_TO_dnVec2(dnQin,dnQout,nderiv)

        i1 = HyperSpheTransfo%list_HyperSphe(1)
        i2 = HyperSpheTransfo%list_HyperSphe(2)
        !write(out_unitp,*) 'i1,i2',i1,i2
        CALL sub_dnVec_TO_dnS(dnQin,dnRho,i1,nderiv)
        CALL sub_dnVec_TO_dnS(dnQin,dnTheta,i2,nderiv)

        CALL sub_dnS1_TO_dntR2(dnTheta,dnCosTheta,2,nderiv)
        CALL sub_dnS1_TO_dntR2(dnTheta,dnSinTheta,3,nderiv)

        CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnRho,dnCosTheta,dnX,nderiv)
        CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnRho,dnSinTheta,dnY,nderiv)

        CALL sub_dnS_TO_dnVec(dnX,dnQout,i1,nderiv)
        CALL sub_dnS_TO_dnVec(dnY,dnQout,i2,nderiv)

       !write(out_unitp,*) 'hyper rho,theta',dnQin%d0(i1)/pi*180._Rkind,dnQin%d0(i2)/pi*180._Rkind
       !write(out_unitp,*) 'hyper x,y',dnQout%d0(i1)/pi*180._Rkind,dnQout%d0(i2)/pi*180._Rkind
      ELSE
        CALL sub_dnVec1_TO_dnVec2(dnQout,dnQin,nderiv)

        i1 = HyperSpheTransfo%list_HyperSphe(1)
        i2 = HyperSpheTransfo%list_HyperSphe(2)
        write(out_unitp,*) 'hyper i1,i2',i1,i2
        write(out_unitp,*) 'hyper x,y',dnQout%d0(i1),dnQout%d0(i2)

        CALL sub_dnVec_TO_dnS(dnQout,dnX,i1,nderiv)
        CALL sub_dnVec_TO_dnS(dnQout,dnY,i2,nderiv)

        dnRho%d0 = sqrt(dnX%d0**2 + dnY%d0**2)
        dnTheta%d0 = atan2(dnY%d0,dnX%d0)
        CALL dihedral_range(dnTheta%d0,2)

        CALL sub_dnS_TO_dnVec(dnRho,dnQin,i1,nderiv)
        CALL sub_dnS_TO_dnVec(dnTheta,dnQin,i2,nderiv)

        write(out_unitp,*) 'hyper rho,theta',dnQin%d0(i1),dnQin%d0(i2)
        write(out_unitp,*) 'hyper rho,theta',dnQin%d0(i1),dnQin%d0(i2)/pi*180._Rkind

      END IF

      CALL dealloc_dnSVM(dnRho)
      CALL dealloc_dnSVM(dnTheta)
      CALL dealloc_dnSVM(dnCosTheta)
      CALL dealloc_dnSVM(dnSinTheta)
      CALL dealloc_dnSVM(dnX)
      CALL dealloc_dnSVM(dnY)

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------
      END SUBROUTINE calc_HyperSpheTransfo

      END MODULE mod_HyperSpheTransfo

