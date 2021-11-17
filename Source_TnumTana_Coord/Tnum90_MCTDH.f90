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
      PROGRAM Tnum_f90
      USE mod_system
      USE mod_dnSVM
      USE mod_Constant
      USE mod_Coord_KEO
      USE mod_PrimOp
      IMPLICIT NONE


!     - parameters for para_Tnum -----------------------
      TYPE (constant)  :: const_phys
      TYPE (CoordType) :: mole
      TYPE (Tnum)      :: para_Tnum
      TYPE (PrimOp_t)  :: PrimOp

      real (kind=Rkind) :: vep,rho
      real (kind=Rkind), pointer :: Tdef2(:,:) => null()
      real (kind=Rkind), pointer :: Tdef1(:) => null()

      real (kind=Rkind), pointer :: Tcor2(:,:) => null()
      real (kind=Rkind), pointer :: Tcor1(:) => null()
      real (kind=Rkind), pointer :: Trot(:,:) => null()

      integer :: nderiv
      real (kind=Rkind), allocatable :: Qact(:),Qxyz(:)

!     - working parameters ------------------------------------------
      integer :: i,n

      character (len=*), parameter :: name_sub='Tnum90_MCTDH'

!===========================================================
!===========================================================
      !para_mem%mem_debug = .TRUE.
      CALL versionEVRT(.TRUE.)
      !print_level=-1
      print_level=2
      !-----------------------------------------------------------------
      !     - read the coordinate transformations :
      !     -   zmatrix, polysperical, bunch...
      !     ------------------------------------------------------------
      para_Tnum%MCTDHForm = .TRUE.
      CALL Read_CoordType(mole,para_Tnum,const_phys)
      !     ------------------------------------------------------------
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !     - read coordinate values -----------------------------------
      !     ------------------------------------------------------------
      CALL read_RefGeom(mole,para_Tnum)
      !     ------------------------------------------------------------
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !     ---- TO finalize the coordinates (NM) and the KEO ----------
      !     ------------------------------------------------------------
      CALL time_perso('Tnum90_MCTDH')
      CALL Finalize_TnumTana_Coord_PrimOp(para_Tnum,mole,PrimOp)
      CALL time_perso('Tnum90_MCTDH')

      !-----------------------------------------------------------------
      n=10! (several) evaluations (for cpu time)
!===========================================================
!===========================================================

      CALL alloc_NParray(Qact,(/ mole%nb_var /),'Qact',name_sub)
      CALL get_Qact0(Qact,mole%ActiveTransfo)

!-------------------------------------------------
!     - Cartesian coordinates --------------------
!     --------------------------------------------

       write(out_unitp,*) "======================================"
       write(out_unitp,*) "======================================"
       write(out_unitp,*) "======================================"
       write(out_unitp,*) "======================================"
       mole%WriteCC = .FALSE.
       CALL alloc_NParray(Qxyz,(/ mole%ncart_act /),'Qxyz',name_sub)

       write(out_unitp,*) "======================================"
       CALL time_perso('Qact->Qxyz')

       DO i=1,n
         CALL sub_QactTOd0x(Qxyz,Qact,mole,Gcenter=.FALSE.)
       END DO

       CALL time_perso('Qact->Qxyz')
       write(out_unitp,*) "======================================"

       mole%WriteCC = .TRUE.
       CALL sub_QactTOd0x(Qxyz,Qact,mole,Gcenter=.FALSE.)

       write(out_unitp,*) ' Cartesian coordinates:'
       CALL Write_Cartg98(Qxyz,mole)

       mole%WriteCC = .FALSE.

       CALL dealloc_NParray(Qxyz,'Qxyz',name_sub)
       write(out_unitp,*) "======================================"
       write(out_unitp,*) "======================================"
       write(out_unitp,*) "======================================"
       write(out_unitp,*) "======================================"
       write(out_unitp,*) n,' evaluation of sub_QactTOd0x'
       write(out_unitp,*) "======================================"

!-------------------------------------------------
!-------------------------------------------------

!-------------------------------------------------
!-------------------------------------------------
!     - calculation of f2, f1, vep, rho ----------
!     --------------------------------------------

       CALL alloc_array(Tdef2,(/ mole%nb_act,mole%nb_act /),'Tdef2',name_sub)
       CALL alloc_array(Tdef1,(/ mole%nb_act /),            'Tdef1',name_sub)
       CALL alloc_array(Tcor2,(/ mole%nb_act,3 /),          'Tcor2',name_sub)
       CALL alloc_array(Tcor1,(/ 3 /),                      'Tcor1',name_sub)
       CALL alloc_array(Trot, (/ 3,3 /),                    'Trot', name_sub)

       write(out_unitp,*) "======================================"
       write(out_unitp,*) "======================================"
       write(out_unitp,*) "======================================"
       write(out_unitp,*) n,' evaluation of calc3_f2_f1Q_num'
       write(out_unitp,*) "======================================"
       para_Tnum%WriteT = .FALSE.
       CALL time_perso('f2')
       DO i=1,n-1

        CALL    calc3_f2_f1Q_num(Qact,                                  &
                                 Tdef2,Tdef1,vep,rho,                   &
                                 Tcor2,Tcor1,Trot,                      &
                                 para_Tnum,mole)
       END DO
       para_Tnum%WriteT = .TRUE.
       CALL    calc3_f2_f1Q_num(Qact,                                   &
                                 Tdef2,Tdef1,vep,rho,                   &
                                 Tcor2,Tcor1,Trot,                      &
                                 para_Tnum,mole)
       para_Tnum%WriteT = .FALSE.
       CALL time_perso('f2')
       write(out_unitp,*) "======================================"
       write(out_unitp,*) "======================================"
       write(out_unitp,*) "======================================"
       write(out_unitp,*) n,' evaluation of calc3_f2_f1Q_num'
       write(out_unitp,*) "======================================"

!-------------------------------------------------
!      FOR MCTDH (the calc3_f2_f1Q_num subroutine can be supressed)
       IF (para_Tnum%MCTDHform) CALL export3_MCTDH_T(Qact,para_Tnum,mole)
!-------------------------------------------------

       CALL dealloc_array(Tdef2,'Tdef2',name_sub)
       CALL dealloc_array(Tdef1,'Tdef1',name_sub)
       CALL dealloc_array(Tcor2,'Tcor2',name_sub)
       CALL dealloc_array(Tcor1,'Tcor1',name_sub)
       CALL dealloc_array(Trot, 'Trot', name_sub)

       CALL dealloc_CoordType(mole)
       CALL dealloc_NParray(Qact,'Qact',name_sub)

       write(out_unitp,*) 'mem_tot,max_mem_used',para_mem%mem_tot,para_mem%max_mem_used
       write(out_unitp,*) 'nb_alloc,nb_dealloc',para_mem%nb_alloc,para_mem%nb_dealloc
       write(out_unitp,*) 'END Tnum'

      END PROGRAM Tnum_f90
