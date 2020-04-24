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


!=====================================================================
!
!      Determination of the nD-harmonic functions: nb_bi
!      and the corresponding table: ind_harm(nb_var+1,nb_bi)
!      as a function of max_excit (the degree maximal of the hermite polynomia)
!      and the harmonic frequencies (at the Qdyn0).
!
!      The nD-functions are sorted (or not) relatively to the harmonic energy.
!
!=====================================================================
      SUBROUTINE sub2_ind_harm(Basis2n,para_PES,para_Tnum,mole)
      use mod_system
      USE mod_nDindex
      use mod_Constant,  only: get_conv_au_to_unit
      USE mod_Coord_KEO, only: CoordType, Tnum, gaussian_width, get_Qact0
      use mod_PrimOp,    only: param_pes, sub_freq2_rph
      USE mod_basis
      IMPLICIT NONE

!----- for the Basis2n -----------------------------------------------
      TYPE (basis)   :: Basis2n

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole
      TYPE (Tnum)    :: para_Tnum

!----- for the PES ---------------------------------------------------
      TYPE (param_PES)   :: para_PES

!------ working variables -------------------------------------------

      real (kind=Rkind) :: norm2

      real (kind=Rkind) :: d0ehess(mole%nb_inact2n)
      real (kind=Rkind) :: d0Qeq(mole%nb_inact2n)
      real (kind=Rkind) :: fe(mole%nb_inact2n+1)
      real (kind=Rkind) :: d0hess(mole%nb_inact2n,mole%nb_inact2n)
      real (kind=Rkind) :: d0g(mole%nb_inact2n),pot0_corgrad

      real (kind=Rkind) :: d0c(mole%nb_inact2n,mole%nb_inact2n),        &
                           d0c_inv(mole%nb_inact2n,mole%nb_inact2n),    &
                           A(mole%nb_inact2n,mole%nb_inact2n)

      real (kind=Rkind) :: Qact(mole%nb_var)


      logical       :: deriv,num
      integer       :: i,k,nderiv
      integer       :: min_i(mole%nb_inact2n)

      integer :: nb_somme_herm,i_point,nb_coupling,n_h
      real (kind=Rkind) :: ene_freq,ZPE,ene,auTOcm_inv

!---------------------------------------------------------------------
!---------------------------------------------------------------------


!---------------------------------------------------------------------
      character (len=*), parameter :: name_sub='sub2_ind_harm'
      !logical, parameter :: debug=.TRUE.
      logical, parameter :: debug=.FALSE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) ' BEGINNING ',name_sub
        write(out_unitp,*) '  Qdyn0',mole%ActiveTransfo%Qdyn0
        write(out_unitp,*) '  nb_inact2n',mole%nb_inact2n
        write(out_unitp,*) '  nb_var',mole%nb_var
      END IF
!---------------------------------------------------------------------
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')

      IF (mole%nb_inact2n == 0 .AND. para_PES%nb_elec == 1) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' nb_inact2n (nb_inact21+nb_inact22)= 0'
        write(out_unitp,*) ' no harmonic inactive variables !!'
        STOP
      END IF
!------------------------------------------------------------

      CALL get_Qact0(Qact,mole%ActiveTransfo)

      CALL sub_freq2_RPH(d0ehess,d0c,d0c_inv,norm2,d0hess,d0Qeq,d0g,    &
                         pot0_corgrad,                                  &
                         Qact,para_Tnum,mole,mole%RPHTransfo_inact2n)

      write(out_unitp,*) 'freq',d0ehess(:)*auTOcm_inv
      write(out_unitp,*) 'd0Qeq',d0Qeq
      write(out_unitp,*) 'd0c'
      CALL Write_Mat(d0c,out_unitp,5)
      write(out_unitp,*)
!-----------------------------------------------------------------
      n_h = Basis2n%nDindB%Max_nDI

      CALL gaussian_width(mole%nb_inact2n,A,d0c)

      IF (.NOT. Basis2n%nDindB%With_L) THEN
        Basis2n%nDindB%nDweight = d0ehess ! change the weight with the frequnecies
      END IF
      n_h = Basis2n%nDindB%Max_nDI
      CALL init_nDindexPrim(Basis2n%nDindB,mole%nb_inact2n,min_i,With_init=.FALSE.)
      CALL sort_nDindex(Basis2n%nDindB)
      IF (n_h > 0) Basis2n%nDindB%Max_nDI = min(n_h,Basis2n%nDindB%Max_nDI)
      !CALL write_nDindex(Basis2n%nDindB)
!-----------------------------------------------------------------


!     -----------------------------------------------------
!     write label of the anharmonic basis functions
      ZPE = HALF*sum(d0ehess)
      DO i=1,Basis2n%nDindB%max_nDI
        Ene = ZPE + sum(real(Basis2n%nDindB%Tab_nDval(:,i),kind=Rkind)*d0ehess(:))
        write(out_unitp,'(i4,2(1x,f16.4))',advance='no') i,             &
                                     Ene*auTOcm_inv,(Ene-ZPE)*auTOcm_inv
        DO k=1,mole%nb_inact2n
          write(out_unitp,'(i4)',advance='no') Basis2n%nDindB%Tab_nDval(k,i)
        END DO
        write(out_unitp,*)
      END DO
!     -----------------------------------------------------


      IF (Basis2n%nDindB%max_nDI < 1) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' nb_bi < 1 !!',Basis2n%nDindB%max_nDI
        STOP
      END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        CALL write_nDindex(Basis2n%nDindB,'Basis2n%nDindB')
        write(out_unitp,*) ' number of nD inactive harmonic functions, nb_bi: ', &
                    Basis2n%nDindB%max_nDI
        write(out_unitp,*) ' END ',name_sub
      END IF
!---------------------------------------------------------------------

      END SUBROUTINE sub2_ind_harm
