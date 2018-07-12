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
!     RoVibrational levels calculation
!================================================================
      SUBROUTINE sub_VibRot(Tab_Psi,nb_psi,para_H,para_ana)

      USE mod_system
      USE mod_Tnum
      USE mod_basis
      !USE mod_psi
      USE mod_psi_set_alloc
      USE mod_psi_Op
      USE mod_Op
      USE mod_OpPsi
      USE mod_constant
      USE mod_analysis
      IMPLICIT NONE


!----- variables for the WP propagation ----------------------------
      integer            :: nb_psi
      TYPE (param_psi)   :: Tab_Psi(nb_psi)

!----- Operator variables ----------------------------------------------
      TYPE (param_Op)  :: para_H
      logical          :: print_Op

!----- variables pour la namelist analyse ------------------------------
      TYPE (param_ana)           :: para_ana



!---- variable for the Z-matrix ----------------------------------------
      TYPE (zmatrix), pointer    :: mole
      TYPE (Tnum), pointer       :: para_Tnum

!----- working variables -----------------------------
      TYPE (param_d0MatOp)              :: MatRV
      TYPE (param_psi)                  :: OpPsi
      integer                           :: iv,jv

      integer           :: nb_bVR
      real (kind=Rkind), allocatable ::    H_VR(:,:)
      real (kind=Rkind), allocatable ::    Ene_VR(:)
      real (kind=Rkind), allocatable ::    Vec_VR(:,:)
      real (kind=Rkind) ::    rho_V(nb_psi)

      integer              :: nb_bRot,JRot,iterm_Op,iterm_BasisRot
      integer              :: J1,J2,i1,i2,f1,f2,ibRot,jbRot,nb_ana
      integer              :: i,jR,nb_shift,type_Op

      real (kind=Rkind)    :: Val_BasisRot,non_hermitic,auTOcm_inv
      complex (kind=Rkind) :: C_over



!----- for debuging --------------------------------------------------
      integer   :: err
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'sub_VibRot'
!-----------------------------------------------------------
      mole       => para_H%mole
      para_Tnum  => para_H%para_Tnum

      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) 'nb_psi',nb_psi

      IF (debug) THEN
        write(out_unitp,*)
      END IF
      CALL flush_perso(out_unitp)

      !-----------------------------------------------------------
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')
      JRot = para_ana%JJmax
      para_H%Mat_done = .FALSE.

      type_Op = para_H%para_PES%Type_HamilOp ! H
      IF (type_Op /= 1) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) '    Type_HamilOp MUST be equal to 1 here!!'
        write(out_unitp,*) '    CHECK your data!!'
        STOP
      END IF

      CALL Init_d0MatOp(MatRV,type_Op,0,nb_psi,JRot=JRot,cplx=para_H%cplx)
      para_H%Make_Mat = .FALSE.
      DO iterm_Op=1,MatRV%nb_term
        DO iv=1,nb_psi
          IF (debug) THEN
            write(out_unitp,*) '======================================='
            write(out_unitp,*) '======================================='
            write(out_unitp,*) '======================================='
            write(out_unitp,*) "iterm_Op,iv",iterm_Op,iv
            CALL flush_perso(out_unitp)
          END IF

          CALL sub_OpPsi(Tab_Psi(iv),OpPsi,para_H,MatRV%derive_termQact(:,iterm_Op))
write(6,*) 'coucou oppsi' ; flush(6)
          DO jv=1,nb_psi
            CALL Overlap_psi1_psi2(C_over,Tab_Psi(jv),OpPsi)
            !write(6,*) 'jv,iv,C_over',jv,iv,C_over
            MatRV%ReVal(jv,iv,iterm_Op) = real(C_over,kind=Rkind)
            IF (MatRV%cplx) MatRV%ImVal(jv,iv) = aimag(C_over)
          END DO
        END DO
      END DO
      IF (debug) CALL Write_d0MatOp(MatRV)
      CALL flush_perso(out_unitp)

      !SET Rotational basis
      CALL init_RotBasis_Param(para_H%BasisnD%RotBasis,Jrot)
      IF (debug) CALL Write_RotBasis_Param(para_H%BasisnD%RotBasis)
      CALL flush_perso(out_unitp)

      nb_bRot = para_H%BasisnD%RotBasis%nb_Rot
      nb_bVR  = nb_psi*nb_bRot


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      CALL alloc_NParray(H_VR,  (/ nb_bVR,nb_bVR /),'H_VR',  name_sub)
      CALL alloc_NParray(Vec_VR,(/ nb_bVR,nb_bVR /),'Vec_VR',name_sub)
      CALL alloc_NParray(Ene_VR,(/ nb_bVR /),       'Ene_VR',name_sub)

      H_VR(:,:) = ZERO

      DO iterm_Op=1,MatRV%nb_term

          ! Rotational contribution
          J1       = MatRV%derive_termQact(1,iterm_Op)
          J2       = MatRV%derive_termQact(2,iterm_Op)
          iterm_BasisRot = para_H%BasisnD%RotBasis%tab_der_TO_iterm(J1,J2)
          !write(6,*) 'J1,J2',J1,J2,'iterm_Op,iterm_BasisRot',iterm_Op,iterm_BasisRot

          DO ibRot=1,nb_bRot
          DO jbRot=1,nb_bRot
            Val_BasisRot = para_H%BasisnD%RotBasis%tab_RotOp(ibRot,jbRot,iterm_BasisRot)
            IF (abs(Val_BasisRot) < ONETENTH**9) CYCLE

            i1 = (ibRot-1)*nb_psi
            i2 = (jbRot-1)*nb_psi
            f1 = nb_psi
            f2 = nb_psi

            !IF (debug) THEN
            !  write(out_unitp,*) 'J1,J2',J1,J2
            !  write(out_unitp,*) 'ibRot,i1+1:i1+f1',ibRot,i1+1,i1+f1
            !  write(out_unitp,*) 'jbRot,i2+1:i2+f2',jbRot,i2+1,i2+f2
            !END IF

            H_VR(i1+1:i1+f1 , i2+1:i2+f2) =                             &
                                        H_VR(i1+1:i1+f1 , i2+1:i2+f2) + &
                          MatRV%ReVal(1:f1 , 1:f2,iterm_Op)*Val_BasisRot
          END DO
          END DO
          !---- END LOOP on the rotational basis function


      END DO

      IF (debug) THEN
        write(out_unitp,*) 'H_VR'
        CALL Write_Mat(H_VR,out_unitp,5)
      END IF

      CALL sub_hermitic_H(H_VR,nb_bVR,non_hermitic,para_H%sym_Hamil)

      IF (non_hermitic >= FOUR/TEN**4) THEN
        write(out_unitp,*) 'WARNING: non_hermitic is BIG'
        write(out_unitp,31) non_hermitic
 31     format(' Hamiltonien: ',f16.12,' au')
      ELSE
        IF (print_level>-1) write(out_unitp,21) non_hermitic*auTOcm_inv
 21     format(' Hamiltonien: ',f16.12,' cm-1')
      END IF
      CALL flush_perso(out_unitp)


      IF (para_H%sym_Hamil) THEN
        CALL diagonalization(H_VR,Ene_VR,Vec_VR,nb_bVR,3,1,.TRUE.)
      ELSE
        CALL diagonalization(H_VR,Ene_VR,Vec_VR,nb_bVR,4,1,.TRUE.)
      END IF
      nb_shift = count(Ene_VR(:) <= para_H%ComOp%ZPE)
      IF (nb_shift > 0) THEN
        write(out_unitp,*) 'WARNING the vectors 1 to ',nb_shift,'have negative energies',Ene_VR(1:nb_shift)
        write(out_unitp,*) '=> They will be shifted'
        Vec_VR = cshift(Vec_VR,shift=nb_shift,dim=2)
        Ene_VR = cshift(Ene_VR,shift=nb_shift)
      END IF

      IF (debug) THEN
        write(out_unitp,*) 'Vec_VR (in column)'
        CALL Write_Mat(Vec_VR,out_unitp,5)
      END IF


      write(out_unitp,*) 'ZPE',para_H%ComOp%ZPE*auTOcm_inv
      nb_ana = min(nb_bVR,nb_bRot*2)
      write(out_unitp,'(A,i4,30f15.6)') 'Ene RV',JRot,                  &
                         (Ene_VR(1:nb_bRot)-para_H%ComOp%ZPE)*auTOcm_inv
      write(out_unitp,'(A,i4,30f15.6)') 'Ene RV',JRot,                  &
          (Ene_VR(nb_bRot+1:nb_ana)-real(Tab_psi(2)%CAvOp,kind=Rkind))* &
                                                             auTOcm_inv

      write(out_unitp,*) 'Ene RV (all), J:',JRot
      DO i=1,min(nb_bVR,10*nb_bRot)
        write(out_unitp,'(A,i5,f15.6)') ' levR:',i,                     &
                                 (Ene_VR(i)-para_H%ComOp%ZPE)*auTOcm_inv
        rho_V(:) = ZERO
        DO jR=1,nb_bRot
          i1 = (jR-1)*nb_psi + 1
          f1 = (jR-1)*nb_psi + nb_psi
          rho_V(:) = rho_V(:) + abs(Vec_VR( i1:f1 ,i))**2
        END DO
        write(out_unitp,'(A,10f5.2)') ' densVib:',rho_V(1:min(10,nb_psi))

      END DO
      CALL flush_perso(out_unitp)

      CALL dealloc_NParray(H_VR,  'H_VR',  name_sub)
      CALL dealloc_NParray(Vec_VR,'Vec_VR',name_sub)
      CALL dealloc_NParray(Ene_VR,'Ene_VR',name_sub)

      CALL dealloc_d0MatOp(MatRV)

!----------------------------------------------------------
      IF (debug) THEN
      END IF
      write(out_unitp,*) 'END ',name_sub
      CALL flush_perso(out_unitp)

!----------------------------------------------------------

      end subroutine sub_VibRot
