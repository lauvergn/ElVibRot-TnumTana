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
     RECURSIVE SUBROUTINE Set_SymAbelian_OF_BasisDP(basis_DP)
      USE mod_system
      USE mod_nDindex
      USE mod_basis
      IMPLICIT NONE

      TYPE (basis), intent(inout) :: basis_DP

!---------------------------------------------------------------------
      integer           :: ib,ibi,ibasis,LG
      integer           :: symab,symab_ibasis
      integer           :: nDval(basis_DP%nDindB%ndim)

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Set_SymAbelian_OF_BasisDP'
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL Write_SymAbelian(basis_DP%P_SymAbelian)
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------------------

      IF (basis_DP%nb_basis < 1) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' This subroutine works only when nb_basis > 0'
        write(out_unitp,*) 'nb_basis: ',basis_DP%nb_basis
        STOP
      END IF
      CALL alloc_SymAbelian(basis_DP%P_SymAbelian,basis_DP%nb)
      CALL Set_ReadsymabOFSymAbelian(basis_DP%P_SymAbelian,Read_symab=0)

      SELECT CASE (basis_DP%SparseGrid_type)
      CASE (0) ! Direct product

        DO ib=1,basis_DP%nb
          CALL calc_nDindex(basis_DP%nDindB,ib,nDval)
          symab = 0
          DO ibasis=1,basis_DP%nb_basis
            ibi = nDval(ibasis)
            symab_ibasis = Get_symabOFSymAbelian_AT_ib(                   &
                      basis_DP%tab_Pbasis(ibasis)%Pbasis%P_SymAbelian,ibi)
            symab = Calc_symab1_EOR_symab2(symab,symab_ibasis)
            IF (symab == -1) EXIT
          END DO
          CALL Set_symabOFSymAbelian_AT_ib(basis_DP%P_SymAbelian,ib,symab)

          IF (debug) write(out_unitp,*) 'ib,nDval',ib,nDval, &
                                        'symab',WriteTOstring_symab(symab)
        END DO
        CALL Set_nbPERsym_FROM_SymAbelian(basis_DP%P_SymAbelian)

      CASE (1) ! Sparse basis (Smolyak 1st implementation)
        ! recursive call for the first Smolyak grid
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' It should not be called with SparseGrid_type=1'
        STOP

!        CALL Set_SymAbelian_OF_BasisDP(basis_DP%tab_PbasisSG(1)%Pbasis)
!
!        CALL SymAbelian1_TO_SymAbelian2(basis_DP%P_SymAbelian,          &
!                          basis_DP%tab_PbasisSG(1)%Pbasis%P_SymAbelian)
!        CALL Set_nbPERsym_FROM_SymAbelian(basis_DP%P_SymAbelian)

      CASE (2,4) ! Sparse basis (Smolyak 2d or 4th implementation)

        LG = basis_DP%L_SparseGrid

        CALL init_nDval_OF_nDindex(basis_DP%nDindB,nDval)
        DO ib=1,basis_DP%nb
          CALL ADD_ONE_TO_nDindex(basis_DP%nDindB,nDval)
          !CALL calc_nDindex(basis_DP%nDindB,ib,nDval)

          symab = 0
          DO ibasis=1,basis_DP%nb_basis
            ibi = nDval(ibasis)
            symab_ibasis = Get_symabOFSymAbelian_AT_ib(                 &
                   basis_DP%tab_basisPrimSG(LG,ibasis)%P_SymAbelian,ibi)
            symab = Calc_symab1_EOR_symab2(symab,symab_ibasis)
            IF (symab == -1) EXIT
          END DO
          CALL Set_symabOFSymAbelian_AT_ib(basis_DP%P_SymAbelian,ib,symab)

          IF (debug) write(out_unitp,*) 'ib,nDval',ib,nDval, &
                                        'symab',WriteTOstring_symab(symab)
        END DO
        CALL Set_nbPERsym_FROM_SymAbelian(basis_DP%P_SymAbelian)


      CASE DEFAULT
        write(out_unitp,*) ' ERROR in',name_sub
        write(out_unitp,*) ' WRONG SparseGrid_type',basis_DP%SparseGrid_type
        write(out_unitp,*) ' The possibilities are: 0, 1, 2, 4'
        STOP
      END SELECT

      IF (basis_DP%nb_basis > 0) THEN ! FOR DP

      END IF

!-----------------------------------------------------------
      IF (debug) THEN
        CALL Write_SymAbelian(basis_DP%P_SymAbelian)
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------------------


      END SUBROUTINE Set_SymAbelian_OF_BasisDP

      integer FUNCTION Get_symabOFSymAbelianOFBasis_AT_ib(basis_set,ib)
         USE mod_system
         USE mod_basis
         IMPLICIT NONE

         TYPE (basis), intent(in) :: basis_set
         integer, intent(in) :: ib


         character (len=*), parameter :: name_sub='Get_symabOFSymAbelianOFBasis_AT_ib'

         Get_symabOFSymAbelianOFBasis_AT_ib =                           &
                    Get_symabOFSymAbelian_AT_ib(basis_set%P_SymAbelian,ib)

      END FUNCTION Get_symabOFSymAbelianOFBasis_AT_ib

      integer FUNCTION Get_nbPERsym_FROM_SymAbelianOFBasis(basis_set,symab)
      USE mod_system
      USE mod_basis
      IMPLICIT NONE

      TYPE (basis), intent(inout) :: basis_set
      integer, intent(in) :: symab

      character (len=*), parameter :: name_sub='Get_nbPERsym_FROM_SymAbelianOFBasis'

      Get_nbPERsym_FROM_SymAbelianOFBasis =                             &
               Get_nbPERsym_FROM_SymAbelian(basis_set%P_SymAbelian,symab)

      END FUNCTION Get_nbPERsym_FROM_SymAbelianOFBasis
      integer FUNCTION Get_nbPERsym_FROM_SymAbelianOFAllBasis(AllBasis,symab)
      USE mod_system
      USE mod_basis
      IMPLICIT NONE

      TYPE (param_AllBasis), intent(inout) :: AllBasis
      integer, intent(in) :: symab

      integer :: Get_nbPERsym_FROM_SymAbelianOFBasis ! function


      character (len=*), parameter :: name_sub='Get_nbPERsym_FROM_SymAbelianOFAllBasis'

      IF (AllBasis%nb_be > 1 .OR. AllBasis%nb_bi > 1) THEN
        Get_nbPERsym_FROM_SymAbelianOFAllBasis =                        &
              Get_nbPERsym_FROM_SymAbelianOFBasis(AllBasis%basisnD,-1)* &
                                           AllBasis%nb_be*AllBasis%nb_bi
      ELSE
        Get_nbPERsym_FROM_SymAbelianOFAllBasis =                        &
            Get_nbPERsym_FROM_SymAbelianOFBasis(AllBasis%basisnD,symab)
      END IF

      END FUNCTION Get_nbPERsym_FROM_SymAbelianOFAllBasis
