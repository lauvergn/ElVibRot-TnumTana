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
      use mod_system
      use mod_dnSVM
      use mod_Constant
      use mod_Coord_KEO

      IMPLICIT NONE

!     - parameters for para_Tnum -----------------------
      TYPE (constant)  :: const_phys
      TYPE (CoordType) :: mole
      TYPE (Tnum)      :: para_Tnum

      real (kind=Rkind) :: vep,rho
      real (kind=Rkind), pointer :: Tdef2(:,:)=>null()
      real (kind=Rkind), pointer :: Tdef1(:)=>null()
      real (kind=Rkind), pointer :: Tcor2(:,:)=>null()
      real (kind=Rkind), pointer :: Tcor1(:)=>null()
      real (kind=Rkind), pointer :: Trot(:,:)=>null()

      TYPE(Type_dnMat) :: dng,dnGG
      TYPE(Type_dnVec) :: dnx

!     ------------------------------------------------------

!     - for the coordinate values ----------------------------------
      real (kind=Rkind), allocatable :: Qact(:)

!     - working parameters ------------------------------------------
      integer :: nderiv,err_mem,memory,err_read

      character (len=*), parameter :: name_sub='TEST_TnumTana'


!=======================================================================
!=======================================================================
      CALL versionEVRT(.TRUE.)
      print_level=2

      !-----------------------------------------------------------------
      !     - read the coordinate transformations :
      !     -   zmatrix, polysperical, bunch...
      !     ------------------------------------------------------------
      CALL Read_CoordType(mole,para_Tnum,const_phys)
      !     ------------------------------------------------------------
      !-----------------------------------------------------------------

      IF (associated(mole%NMTransfo) .OR. associated(mole%RPHTransfo)) THEN
        write(out_unitp,*) "ERROR: This test program cannot be used with"
        write(out_unitp,*) "Normal modes (NM) or RPH"
        STOP
      END IF

      !-----------------------------------------------------------------
      !     - read coordinate values -----------------------------------
      !     ------------------------------------------------------------
      CALL read_RefGeom(mole,para_Tnum)
      !     ------------------------------------------------------------
      !-----------------------------------------------------------------
!=======================================================================
!=======================================================================

!===========================================================
!===========================================================

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
        CALL time_perso('sub_QactTOdnx')

        nderiv = 0
        CALL alloc_dnSVM(dnx,mole%ncart,mole%nb_act,nderiv)
        write(out_unitp,*) "======================================"
        CALL sub_QactTOdnx(Qact,dnx,mole,nderiv,.FALSE.)
        write(out_unitp,*) 'dnx: ',mole%ncart
        CALL write_dnx(1,mole%ncart,dnx,nderiv)

        CALL Write_Cartg98(dnx%d0,mole)

        CALL dealloc_dnSVM(dnx)
        CALL time_perso('sub_QactTOdnx')
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
!-------------------------------------------------
!-------------------------------------------------


      CALL dealloc_CoordType(mole)
      CALL dealloc_NParray(Qact,'Qact',name_sub)


      write(out_unitp,*) 'END ',name_sub

      END PROGRAM Tnum_f90
