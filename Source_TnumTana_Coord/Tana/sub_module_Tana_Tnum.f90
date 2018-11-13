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

MODULE mod_Tana_Tnum
   !!@description: 
   IMPLICIT NONE

   PRIVATE
   PUBLIC :: comparison_G_FROM_Tnum_Tana

   CONTAINS

   SUBROUTINE comparison_G_FROM_Tnum_Tana(TWOxKEO,mole,para_Tnum,Qact)
      USE mod_system
      USE mod_Tnum
      USE mod_paramQ
      USE mod_Tana_PiEulerRot
      USE mod_Tana_sum_opnd
      USE mod_Tana_op
      USE mod_Tana_NumKEO
      USE mod_Tana_write_mctdh
      USE mod_dnSVM
      USE mod_dnGG_dng
      USE mod_f2f2Vep
      IMPLICIT NONE

      TYPE(sum_opnd),        intent(inout)        :: TWOxKEO
      TYPE (zmatrix),        intent(inout)        :: mole
      TYPE (Tnum),           intent(inout)        :: para_Tnum
      real (kind=Rkind),     intent(inout)        :: Qact(:)


      real (kind=Rkind),pointer  :: Gana(:,:)
      real (kind=Rkind),pointer  :: f2_ana(:,:),f1_ana(:)
      real (kind=Rkind), pointer :: Tdef2(:,:)
      real (kind=Rkind), pointer :: Tdef1(:)

      real (kind=Rkind), pointer :: Tcor2(:,:)
      real (kind=Rkind), pointer :: Tcor1(:)
      real (kind=Rkind), pointer :: Trot(:,:)

      TYPE(Type_dnMat)               :: dng,dnGG
      real(kind=Rkind)               :: vep,vep_ana,rho,rho_ana,max_error,vep_error,maxval_f1

     character (len=*), parameter    :: routine_name='comparison_G_FROM_Tnum_Tana'

!===========================================================
!===========================================================

      write(out_unitp,*) '================================================='
      write(out_unitp,*) ' BEGINNING ',routine_name

      nullify(Gana)
      CALL alloc_array(Gana,(/mole%ndimG,mole%ndimG/),'Gana',routine_name)
      CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,2)
      CALL alloc_dnSVM(dng,mole%ndimG,mole%ndimG,mole%nb_act,2)

      CALL get_NumG_WITH_AnaKEO(TWOxKEO,Qact,mole,Gana,vep_ana)
      CALL get_dnGG_vep(Qact,para_Tnum,mole,dnGG,vep,rho,2)

      IF (maxval(abs(Gana-dnGG%d0))/maxval(abs(dnGG%d0)) > ONETENTH**10) THEN
        write(out_unitp,*) 'G of Tana  '
        CALL write_mat(Gana,out_unitp,4)
        write(out_unitp,*) 'G of Tnum  '
        CALL write_mat(dnGG%d0,out_unitp,4)
        write(out_unitp,*) 'Difference of G'
        CALL write_mat(Gana-dnGG%d0,out_unitp,4)
      END IF
      write(out_unitp,*) '         max diff G: ',maxval(abs(Gana-dnGG%d0))
      write(out_unitp,*) 'Relative max diff G: ',maxval(abs(Gana-dnGG%d0))/maxval(abs(dnGG%d0))
      write(out_unitp,*)

      max_error = maxval(abs(Gana-dnGG%d0))/maxval(abs(dnGG%d0))

      CALL dealloc_array(Gana,'Gana',routine_name)
      CALL dealloc_dnSVM(dnGG)
      CALL dealloc_dnSVM(dng)

      nullify(f2_ana)
      nullify(f1_ana)
      CALL alloc_array(f2_ana,(/mole%nb_act,mole%nb_act/),'f2_ana',routine_name)
      CALL alloc_array(f1_ana,(/mole%nb_act/),'f1_ana',routine_name)

      nullify(Tdef2)
      nullify(Tdef1)
      nullify(Tcor2)
      nullify(Tcor1)
      nullify(Trot)
      CALL alloc_array(Tdef2,(/ mole%nb_act,mole%nb_act /),'Tdef2',routine_name)
      CALL alloc_array(Tdef1,(/ mole%nb_act /),            'Tdef1',routine_name)
      CALL alloc_array(Tcor2,(/ mole%nb_act,3 /),          'Tcor2',routine_name)
      CALL alloc_array(Tcor1,(/ 3 /),                      'Tcor1',routine_name)
      CALL alloc_array(Trot, (/ 3,3 /),                    'Trot', routine_name)

      para_Tnum%Tana = .FALSE.
      CALL   calc3_f2_f1Q_num(Qact,                                  &
                              Tdef2,Tdef1,vep,rho,                   &
                              Tcor2,Tcor1,Trot,                      &
                              para_Tnum,mole)
      para_Tnum%Tana = .TRUE.
      CALL get_Numf2f1vep_WITH_AnaKEO(TWOxKEO,Qact,mole,       &
                                f2_ana,f1_ana,vep_ana,rho_ana)

      IF (vep < ONETENTH**6) THEN
        vep_error = abs(vep-vep_ana)
      ELSE
        vep_error = abs(vep-vep_ana)/vep
      END IF
      vep_error = vep_error / TEN

      maxval_f1 = maxval(abs(Tdef1))
      IF (maxval_f1 < ONETENTH**6) maxval_f1 = ONE

      write(out_unitp,*) '         max diff f2: ',maxval(abs(f2_ana-Tdef2))
      write(out_unitp,*) 'Relative max diff f2: ',maxval(abs(f2_ana-Tdef2))/maxval(abs(Tdef2))
      write(out_unitp,*) '         max diff f1: ',maxval(abs(f1_ana-Tdef1))
      write(out_unitp,*) 'Relative max diff f1: ',maxval(abs(f1_ana-Tdef1))/maxval_f1
      write(out_unitp,*) '        max diff vep: ',abs(vep-vep_ana)
      write(out_unitp,*) '       vep from Tana: ',vep_ana
      write(out_unitp,*) '       vep from Tnum: ',vep

      write(out_unitp,*)
      IF (mole%nb_act == mole%nb_var) max_error = max( max_error, vep_error )
      max_error = max( max_error, maxval(abs(f2_ana-Tdef2))/maxval(abs(Tdef2)) )
      max_error = max( max_error, maxval(abs(f1_ana-Tdef1))/maxval_f1 )

      write(out_unitp,'(a,e8.2)') '         max error: ',max_error

      CALL dealloc_array(Tdef2, 'Tdef2', routine_name)
      CALL dealloc_array(Tdef1, 'Tdef1', routine_name)
      CALL dealloc_array(Tcor2, 'Tcor2', routine_name)
      CALL dealloc_array(Tcor1, 'Tcor1', routine_name)
      CALL dealloc_array(Trot,  'Trot',  routine_name)
      CALL dealloc_array(f2_ana,'f2_ana',routine_name)
      CALL dealloc_array(f1_ana,'f1_ana',routine_name)


      write(out_unitp,*) ' END ',routine_name
      write(out_unitp,*) '================================================='


   END SUBROUTINE comparison_G_FROM_Tnum_Tana


END MODULE mod_Tana_Tnum
