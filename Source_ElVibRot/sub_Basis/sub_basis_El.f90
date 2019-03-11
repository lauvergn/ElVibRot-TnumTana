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

!=============================================================
!
!      Basis set for several diabatic PES
!        For this basis nq must be equal to nb
!
!=============================================================
      SUBROUTINE sub_Basis_El(base)
      USE mod_system
      USE mod_basis
      IMPLICIT NONE

!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis)      :: base



!---------------------------------------------------------------------
!---------------------------------------------------------------------

!---------------------------------------------------------------------


      integer           :: i,Read_symab

!---------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      integer :: err_io,err_mem,memory
      character (len=*), parameter :: name_sub='sub_Basis_El'
      !logical,parameter :: debug=.FALSE.
      logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'nb',base%nb
       END IF
!-----------------------------------------------------------
       IF (NewBasisEl) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) '  A electronic basis is already defined'
         write(out_unitp,*) '   CHECK your data!!'
         STOP
       END IF

     IF (base%nb <= 0) STOP 'ERROR nb<=0'
     NewBasisEl = .TRUE.

     ! here nq=nb or nq=0
     CALL Set_nq_OF_basis(base,0)

!----------------------------------------------------------------------------

      !CALL alloc_xw_OF_basis(base)
      !base%x(:,1) = (/ (i,i=1,base%nb) /)
      !base%w(:)   = ONE/real(base%nb,kind=Rkind)

      !CALL alloc_dnb_OF_basis(base)
      !CALL mat_id(base%dnRGB%d0,base%nb,base%nb)
      base%primitive_done = .TRUE.
      base%packed_done    = .TRUE.

      CALL dealloc_nDindex(base%nDindB)
      base%nDindB%packed = .TRUE.
      CALL init_nDindexPrim(base%nDindB,1,(/ base%nb /))
      base%nDindB%With_L      = .TRUE.
      base%nDindB%Tab_L(:)    = 0
      base%nDindB%Tab_Norm(:) = ZERO





      CALL alloc_SymAbelian(base%P_SymAbelian,base%nb)
      Read_symab = Get_Read_symabOFSymAbelian(base%P_SymAbelian)
      CALL Set_ReadsymabOFSymAbelian(base%P_SymAbelian,Read_symab)
      DO i=1,base%nb
        CALL Set_symabOFSymAbelian_AT_ib(base%P_SymAbelian,i,-1)
      END DO
      CALL Set_nbPERsym_FROM_SymAbelian(base%P_SymAbelian)
      !CALL Write_SymAbelian(base%P_SymAbelian)


!-----------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(base)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_Basis_El
