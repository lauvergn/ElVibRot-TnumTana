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
      USE mod_Coord_KEO
      IMPLICIT NONE


!     - parameters for para_Tnum -----------------------
      TYPE (CoordType)   :: mole
      TYPE (Tnum)      :: para_Tnum
      TYPE (PrimOp_t)  :: PrimOp

      TYPE(Type_dnMat) :: dng,dnGG
      TYPE(Type_dnS), pointer :: GOFdnS(:,:),EigenVecOFdnS(:,:)
      TYPE(Type_dnS)   :: dnS

      integer :: nderiv,nb_act
!     ------------------------------------------------------


!     - for the coordinate values ----------------------------------
      TYPE (param_Q)    :: para_Q

!     - working parameters ------------------------------------------
      integer :: i,j,n,ndim
      integer :: err_mem,memory

!=======================================================================
!=======================================================================
      CALL versionEVRT(.TRUE.)

      !-----------------------------------------------------------------
      !     - read the coordinate tansformations :
      !     -   zmatrix, polysperical, bunch...
      !     ------------------------------------------------------------
      CALL Read_mole(mole,para_Tnum)
      !     ------------------------------------------------------------
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !     - read coordinate values -----------------------------------
      !     ------------------------------------------------------------
      CALL read_RefGeom(para_Q,mole,para_Tnum)
      !     ------------------------------------------------------------
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !     ---- TO finalize the coordinates (NM) and the KEO ----------
      !     ------------------------------------------------------------
      CALL Finalize_TnumTana_Coord_PrimOp(para_Q%Qact,para_Tnum,mole,PrimOp)
      !-----------------------------------------------------------------
!=======================================================================
!=======================================================================


!-------------------------------------------------
!-------------------------------------------------
!     FOR G and g metric tensors
      CALL alloc_dnSVM(dng ,mole%ndimG,mole%ndimG,mole%nb_act,2)
      CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,2)

      write(out_unitp,*) "======================================"
      write(out_unitp,*) "======================================"
      write(out_unitp,*) "======================================"
      write(out_unitp,*) "======================================"
      CALL time_perso('G and g')

      para_Tnum%WriteT    = .TRUE.
      nderiv              = 2
      CALL get_dng_dnGG(para_Q%Qact,para_Tnum,mole,dng,dnGG,nderiv)

      write(out_unitp,*) 'matrix G:',mole%nb_act
      CALL Write_dnSVM(dnGG)


      nb_act = mole%nb_act

      nullify(GOFdnS)
      nullify(EigenVecOFdnS)
      CALL alloc_array(GOFdnS,(/nb_act,nb_act/),'GOFdnS','main')
      CALL alloc_array(EigenVecOFdnS,(/nb_act,nb_act/),'EigenVecOFdnS','main')
      CALL alloc_MatOFdnS(GOFdnS,nb_act,nderiv)
      CALL alloc_MatOFdnS(EigenVecOFdnS,nb_act,nderiv)
      CALL alloc_dnS(dnS,nb_act,nderiv)

      DO i=1,nb_act
      DO j=1,nb_act
         CALL sub_dnMat_TO_dnS(dnGG,i,j,GOFdnS(i,j))
      END DO
      END DO
      write(out_unitp,*) 'matrix G:',mole%nb_act
      CALL Write_MatOFdnS(GOFdnS)


      CALL DIAG_MatOFdnS(GOFdnS,EigenVecOFdnS,type_diago=4)

      !diagonal G matrix
      write(out_unitp,*) 'Diagonal G:'
      CALL Write_MatOFdnS(GOFdnS)

      CALL sub_ZERO_TO_dnS(dnS)
      DO i=1,nb_act
      DO j=1,nb_act
         IF (j == i) CYCLE
         CALL sub_ABSdnS1_PLUS_dnS2_TO_dnS2(GOFdnS(i,j),dnS)
      END DO
      END DO
      write(out_unitp,*) 'non-Diagonal G:?'
      CALL Write_dnS(dnS)



      CALL time_perso('G and g')
      write(out_unitp,*) "======================================"
      write(out_unitp,*) "======================================"
      write(out_unitp,*) "======================================"
      write(out_unitp,*) "======================================"

      CALL dealloc_MatOFdnS(EigenVecOFdnS)
      CALL dealloc_MatOFdnS(GOFdnS)
      CALL dealloc_array(GOFdnS,'GOFdnS','main')
      CALL dealloc_array(EigenVecOFdnS,'EigenVecOFdnS','main')
      CALL dealloc_dnSVM(dnS)
      CALL dealloc_dnSVM(dng)
      CALL dealloc_dnSVM(dnGG)
!-------------------------------------------------
!-------------------------------------------------

      CALL dealloc_CoordType(mole)
      CALL dealloc_param_Q(para_Q)


      write(out_unitp,*) 'END Tnum'

      end program Tnum_f90
