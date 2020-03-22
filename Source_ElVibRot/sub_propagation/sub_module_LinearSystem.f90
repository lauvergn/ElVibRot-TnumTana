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

MODULE mod_LinearSystem
USE mod_Constant
USE mod_MPI 
IMPLICIT NONE

PRIVATE
PUBLIC :: sub_LinearSystem

CONTAINS
!===============================================================================
!
!  Solve with an iterative method: Op|TabPsi> = |TabOpPsi>
!     => find TabPsi from TabOpPsi
!
!===============================================================================
SUBROUTINE sub_LinearSystem(TabOpPsi,TabPsi,Op,para_propa)
      USE mod_system
      USE mod_Op
      USE mod_psi_set_alloc
      USE mod_psi_SimpleOp
      USE mod_ana_psi,        ONLY : norm2_psi
      USE mod_psi_Op,         ONLY : sub_LCpsi_TO_psi
      USE mod_psi_io,         ONLY : sub_save_psi
      USE mod_propa,          ONLY : param_propa,param_Davidson
      USE mod_MPI
      IMPLICIT NONE

      !----- psi ---------------------------------------------
      TYPE (param_psi),        intent(in)    :: OpTabPsi(:)
      TYPE (param_psi),        intent(inout) :: TabPsi(:)

      !----- Operator: Hamiltonian ----------------------------
      TYPE (param_Op),         intent(in)    :: Op

      !----- WP, energy ... -----------------------------------
      TYPE (param_propa),      intent(in)    :: para_propa


      !----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_LinearSystem'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      !-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
        CALL flush_perso(out_unitp)
      END IF
      !-----------------------------------------------------------

      CALL sub_LinearSystem_Jacobi(TabOpPsi,TabPsi,Op,para_propa)

      !----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
      !----------------------------------------------------------

 END SUBROUTINE sub_LinearSystem

 ! solve iteratively Op.psi = OpPsi with:
 !   psi(k) = psi(k-1) + D0^-1(OpPsi - Op.
 SUBROUTINE sub_LinearSystem_Jacobi(TabOpPsi,TabPsi,para_H,Ene,para_propa)
      USE mod_system
      USE mod_Op
      USE mod_psi_set_alloc
      USE mod_psi_SimpleOp
      USE mod_ana_psi,        ONLY : norm2_psi
      USE mod_psi_Op,         ONLY : sub_LCpsi_TO_psi
      USE mod_psi_io,         ONLY : sub_save_psi
      USE mod_propa,          ONLY : param_propa,param_Davidson
      USE mod_MPI
      IMPLICIT NONE

      !----- psi ---------------------------------------------
      TYPE (param_psi), intent(in)    :: OpTabPsi(:)
      TYPE (param_psi), intent(inout) :: TabPsi(:)

      !----- Operator: Hamiltonian ----------------------------
      TYPE (param_Op),  intent(in)    :: para_H

      !----- WP, energy ... -----------------------------------
      TYPE (param_propa), intent(in)  :: para_propa




      !------ working parameters --------------------------------
      TYPE(param_file)  :: Log_file
      integer           :: iunit

      TYPE (param_psi) :: HDiagInv


      !----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_LinearSystem'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      !-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
        CALL flush_perso(out_unitp)
      END IF


      !-----------------------------------------------------------
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')
      auTOene    = get_Conv_au_TO_WriteUnit('E',WriteUnit)

      !------ initialization -------------------------------------
      Log_file%name='LinearSystem.log'
      CALL file_open(Log_file,iunit)


      ! first the diagonal part of Ene-H
      CALL init_psi(HDiagInv,para_H,Op%cplx)
      IF (allocated(para_H%BasisnD%EneH0)) THEN
        write(out_unitp,*) 'precon /= 1. DML'
        IF (OpInvDiag_part%cplx) THEN
          HDiagInv%CvecB(:) = ONE/(Ene-para_H%BasisnD%EneH0(:)) ! approximation of 1/(Ene-H(i,i))
        ELSE
          HDiagInv%RvecB(:) = ONE/(Ene-para_H%BasisnD%EneH0(:)) ! approximation of 1/(Ene-H(i,i))
        END IF
      ELSE
        write(out_unitp,*) 'precon = 1. DML'
        STOP 'ERROR in sub_LinearSystem: Op%BasisnD%EneH0 is not allocated!!'
      END IF

      DO it=1,max_it
      END DO

      !para_propa%file_WP%formatted = .TRUE.
      para_propa%file_WP%name                    = para_propa%para_Davidson%name_file_saveWP
      para_propa%para_Davidson%formatted_file_WP = para_propa%file_WP%formatted

      !CALL time_perso('Davidson psi0')
      IF(MPI_id==0) THEN
        CALL alloc_NParray(psi0,shape(psi),'psi0',name_sub)
        CALL alloc_NParray(Hpsi,shape(psi),'Hpsi',name_sub)
      ENDIF

      DO i=1,max_diago
        IF(MPI_id==0) CALL init_psi( psi(i),para_H,para_H%cplx)
        IF(MPI_id==0) CALL init_psi(Hpsi(i),para_H,para_H%cplx)
      END DO
      IF (With_Grid) THEN
        DO i=1,max_diago
          Hpsi(i)%BasisRep = .FALSE.
          Hpsi(i)%GridRep  = .TRUE.
        END DO
      END IF

      IF(MPI_id==0) THEN
        CALL init_psi(g,para_H,para_H%cplx)
        CALL alloc_psi(g,      BasisRep=With_Basis,GridRep=With_Grid)
        ! read guss on master
        CALL ReadWP0_Davidson(psi,psi0,Vec0,nb_diago,max_diago,   &
                              para_propa%para_Davidson,para_H%cplx)

        ! save the nb_diago wp
        CALL sub_save_psi(psi,nb_diago,para_propa%file_WP)
        write(out_unitp,*) '  sub_save_psi: psi done'
        RealTime = Delta_RealTime(DavidsonTime)
        CALL flush_perso(out_unitp)
      ENDIF ! for MPI_id==0

      !CALL time_perso('Davidson psi0')
!para_mem%mem_debug = .TRUE.
      !===================================================================
      !===================================================================
      ! LOOP
      !===================================================================
      !===================================================================
      nb_diago        = max(1,nb_diago) ! number of Eign value
      epsi            = para_propa%para_Davidson%conv_resi
      norm2g          = HUNDRED * epsi
      conv_Ene        = HUNDRED * epsi
      tab_norm2g(:)   = norm2g
      convergeEne(:)  = .FALSE.
      convergeResi(:) = .FALSE.
      converge(:)     = .FALSE.

      it      = 0
      ndim    = nb_diago
      ndim0   = 0
      iresidu = 0
      conv    = .FALSE.

      IF(MPI_id==0) THEN
        CALL alloc_NParray(vec,(/ndim,ndim/),"vec",name_sub)

        IF (MatOp_omp /= 2) THEN
          nb_thread = 1
        ELSE
          nb_thread = MatOp_maxth
        END IF
        write(out_unitp,*) 'nb_thread in ',name_sub,' : ',nb_thread
        write(out_unitp,*) 'Beginning Davidson iteration'
        CALL flush_perso(out_unitp)

        write(iunit,*) 'Beginning Davidson iteration' ; CALL flush_perso(iunit)
      ENDIF

      !--------------------------------------------------------------------------------
      ! loop for davidson with maximum iter number para_propa%para_Davidson%max_it
      ! careful about the exit when appling MPI
      !--------------------------------------------------------------------------------
      DO it=0,para_propa%para_Davidson%max_it

        IF(MPI_id==0) THEN
          write(out_unitp,*) '--------------------------------------------------'
          write(iunit,*) 'Davidson iteration',it ; CALL flush_perso(iunit)
        ENDIF
        save_WP = .FALSE.
        !CALL time_perso('Beginining it')

        !---------------------------------------------------------------
        !- Hpsi(:) -----------------------------------------------------
        IF (debug) write(out_unitp,*) 'Hpsi(:)',it,ndim,ndim0
        IF (debug) CALL flush_perso(out_unitp)


        CALL sub_MakeHPsi_Davidson(it,psi(1:ndim),Hpsi(1:ndim),Ene,ndim0, &
                                   para_H,para_propa%para_Davidson,iunit)
        !CALL time_perso('MakeHPsi done')

        IF (debug) write(out_unitp,*) 'Hpsi(:) done',it,ndim,ndim0
        IF (debug) CALL flush_perso(out_unitp)
        !- Hpsi(:) -----------------------------------------------------
        !---------------------------------------------------------------

        !---------------------------------------------------------------
        !- build H
        IF (debug) write(out_unitp,*) 'H mat',it,ndim,ndim0
        IF (debug) CALL flush_perso(out_unitp)

        IF(MPI_id==0) THEN
          ! H built from psi and Hpsi
          ! add MPI for this subroutine later
          CALL sub_MakeH_Davidson(it,psi(1:ndim),Hpsi(1:ndim),H,para_propa%para_Davidson)
        ENDIF
        ndim0 = ndim
        !CALL time_perso('MakeH done')

        IF(MPI_id==0) THEN
          ! if symmetric
          CALL sub_hermitic_H(H,ndim,non_hermitic,para_H%sym_Hamil)
          IF (debug) CALL Write_Mat(H,out_unitp,5)

          IF (non_hermitic > FOUR*ONETENTH**4) THEN
            write(out_unitp,*) 'WARNING: non_hermitic is BIG'
            write(out_unitp,31) non_hermitic
31          format(' Hamiltonien: ',f16.12,' au')
          ELSE
            write(out_unitp,51) non_hermitic*auTOcm_inv
51          format(' Hamiltonien: ',f16.12,' cm-1')
          END IF
        ENDIF ! for MPI_id==0

        IF (para_H%sym_Hamil) THEN
          epsi = max(para_propa%para_Davidson%conv_resi,                &
                   TEN**para_propa%para_Davidson%conv_hermitian *       &
                                               non_hermitic)
        ELSE
          epsi = para_propa%para_Davidson%conv_resi
        END IF
        !CALL time_perso('MakeH done')

        IF (debug) write(out_unitp,*) 'H mat',it,ndim,ndim0
        CALL flush_perso(out_unitp)
        !- build H
        !----------------------------------------------------------

        !----------------------------------------------------------
        !- diagonalization
        IF (debug) write(out_unitp,*) 'diago',it,ndim,ndim0
        IF (debug) CALL flush_perso(out_unitp)

        IF(MPI_id==0) THEN
          CALL dealloc_NParray(vec,"vec",name_sub)
          CALL alloc_NParray(vec,(/ndim,ndim/),"vec",name_sub)
        ENDIF
        Ene(:) = ZERO

        ! write(out_unitp,*) 'ndim',ndim
        ! write(out_unitp,*) 'shape ..',shape(H),shape(Vec),shape(Ene),shape(trav)
        IF (para_H%sym_Hamil) THEN
          !IF(MPI_id==0) CALL diagonalization(H,Ene(1:ndim),Vec,ndim,3,1,.FALSE.)
          ! consider the MPI of diagonalization
          IF(MPI_id==0) CALL diagonalization(H,Ene(1:ndim),Vec,ndim,3,1,.True.)
          !CALL diagonalization(H,Ene(1:ndim),Vec,ndim,2,1,.FALSE.)
        ELSE
          !IF(MPI_id==0) CALL diagonalization(H,Ene(1:ndim),Vec,ndim,4,1,.FALSE.)
          IF(MPI_id==0) CALL diagonalization(H,Ene(1:ndim),Vec,ndim,4,1,.True.)
        END IF

        IF (it == 0 .OR. (it > 1 .AND.                                    &
              mod(it-1,para_propa%para_Davidson%num_resetH) == 0) ) THEN
          EneRef(:) = Ene(:)
        END IF

        IF (it == 0) Ene0(:) = Ene(:) + conv_Ene
        !CALL time_perso('diago done')
        IF (debug) write(out_unitp,*) 'diago',it,ndim,ndim0
        IF (debug) CALL flush_perso(out_unitp)
        !- diagonalization
        !----------------------------------------------------------


        !----------------------------------------------------------
        ! Save vec(:) on vec0(:)
        IF (debug) write(out_unitp,*) 'selec',it,ndim,ndim0
        IF (debug) CALL flush_perso(out_unitp)

        IF(MPI_id==0) THEN
          CALL sub_projec_Davidson(Ene,VecToBeIncluded,nb_diago,min_Ene,para_H%para_PES%min_pot,  &
                                   psi,psi0,Vec,Vec0,para_propa%para_Davidson,it,.TRUE.)
          !CALL time_perso('projec done')

          !> MPI note:  para_H%ComOp%ZPE is updated just on maaster now
          IF (para_H%para_ReadOp%Op_Transfo) THEN
            CALL Set_ZPE_OF_ComOp(para_H%ComOp,Ene(1:count(VecToBeIncluded)),forced=.TRUE.)
          ELSE
            CALL Set_ZPE_OF_ComOp(para_H%ComOp,Ene(1:count(VecToBeIncluded)),Ene_min=min_Ene,forced=.TRUE.)
          END IF
          ZPE = para_H%ComOp%ZPE

          IF (debug) write(out_unitp,*) 'selec',it,ndim,ndim0
          IF (debug) CALL flush_perso(out_unitp)
          ! CALL Write_Mat(Vec,out_unitp,5)
          ! Save vec(:) on vec0(:)
          !----------------------------------------------------------

          IF (para_propa%para_Davidson%all_lower_states) THEN
            nb_diago = count((Ene(1:count(VecToBeIncluded))-ZPE) <= para_propa%para_Davidson%max_Ene)
            VecToBeIncluded = .FALSE.
            VecToBeIncluded(1:nb_diago) = .TRUE.
          END IF

          IF (debug) write(out_unitp,*) 'nb_diago,ZPE,min_Ene',nb_diago,ZPE,min_Ene
          IF (debug) write(out_unitp,*) 'it,Ene(:)',it,Ene(1:ndim)*auTOene

          DEne(1:nb_diago) = Ene(1:nb_diago)-Ene0(1:nb_diago)
          conv_Ene = maxval(abs(DEne(1:nb_diago)))
          write(out_unitp,41) 'convergence (it, norm2g/epsi, conv_Ene): ',&
                                                   it,norm2g/epsi,conv_Ene
          IF (para_propa%para_Davidson%Hmax_propa) THEN
            write(out_unitp,21) it,ndim,norm2g/epsi,iresidu,              &
                             -Ene(1:nb_diago)*auTOene
          ELSE
            write(out_unitp,21) it,ndim,norm2g/epsi,iresidu,              &
                              Ene(1:nb_diago)*auTOene
          END IF

          DO j=1,nb_diago
            convergeEne(j) = abs(DEne(j)) < para_propa%para_Davidson%conv_Ene
          END DO
          write(out_unitp,41) 'it Diff Ene (' // trim(WriteUnit) // '):    ',it, &
                           DEne(1:nb_diago) * auTOene
          write(out_unitp,42) 'it convergenceEne(:):  ',it,               &
                                                    convergeEne(1:nb_diago)
          CALL flush_perso(out_unitp)

          write(iunit,21) it,ndim,norm2g/epsi,iresidu,Ene(1:ndim)*auTOene
          CALL flush_perso(iunit)

          !----------------------------------------------------------
          !- residual vector ---------------------------
          !-  and convergence --------------------------
          IF (debug) write(out_unitp,*) 'residual',it,ndim,ndim0
          IF (debug) CALL flush_perso(out_unitp)
          norm2g    = -ONE
          fresidu   = 0
          ! time consuming in MakeResidual_Davidson
          DO j=1,ndim
            IF (.NOT. converge(j) .AND. VecToBeIncluded(j)) THEN
              CALL MakeResidual_Davidson(j,g,psi,Hpsi,Ene,Vec)

              CALL norm2_psi(g)
              tab_norm2g(j) = sqrt(g%norm2)
              IF (fresidu == 0) fresidu = j
              IF (tab_norm2g(j) > norm2g) THEN
                iresidu = j
                norm2g = tab_norm2g(iresidu)
              END IF

              convergeResi(j) = tab_norm2g(j) < epsi
            END IF
            converge(j) = (convergeEne(j) .AND. convergeResi(j))

          END DO
          conv = all(converge(1:nb_diago))
        ENDIF ! for MPI_id==0

#if(run_MPI)
        CALL MPI_BCAST(conv,size1_MPI,MPI_LOGICAL,root_MPI,MPI_COMM_WORLD,MPI_err)
#endif

        IF(MPI_id==0) THEN
          Ene0(1:nb_diago) = Ene(1:nb_diago)
          write(out_unitp,41) 'it tab_norm2g          ',it,tab_norm2g(1:nb_diago)
          write(out_unitp,42) 'it convergenceResi(:): ',it,convergeResi(1:nb_diago)
41        format(a,i3,100(1x,e9.2))
42        format(a,i3,100(1x,l9))

          !CALL time_perso('residual done')

          IF (debug) write(out_unitp,*) 'residual',it,ndim,ndim0
          IF (debug) CALL flush_perso(out_unitp)
          !- residual vector and convergence ------------------------
          !----------------------------------------------------------

          !----------------------------------------------------------
          !- convergence ? ------------------------------------------
          nb_conv_states   = count( converge(1:nb_diago) )
          nb_unconv_states = nb_diago-nb_conv_states
          write(out_unitp,*) 'it, conv',it,conv
          write(out_unitp,*) 'it, nb   converged state(s)',it,nb_conv_states
          write(out_unitp,*) 'it, nb unconverged state(s)',it,nb_unconv_states
          CALL flush_perso(out_unitp)
          !- convergence ? ------------------------------------------
          !----------------------------------------------------------

          CALL sub_NewVec_Davidson(it,psi,Hpsi,Ene,Ene0,EneRef,Vec,       &
                             converge,VecToBeIncluded,nb_diago,max_diago, &
                                   para_propa%para_Davidson,fresidu,ndim, &
                                   para_H%para_ReadOp%Op_Transfo,para_H%para_ReadOp%E0_Transfo)
          !CALL time_perso('NewVec done')

          nb_added_states = ndim-ndim0
          save_WP = (ndim == max_diago) .OR. conv .OR.                    &
                    it == para_propa%para_Davidson%max_it .OR.            &
           (it > 0 .AND. mod(it,para_propa%para_Davidson%num_resetH) == 0)
           Save_WP = Save_WP .AND. .NOT. Hmin_OR_Hmax
           !- new vectors --------------------------------------------
           !----------------------------------------------------------

          write(out_unitp,*) 'ndim,max_diago',ndim,max_diago
          write(out_unitp,*) 'save_WP,conv',save_WP,conv
          IF (debug) THEN
            write(out_unitp,*) 'it,max_it',it,para_propa%para_Davidson%max_it, &
                                     (it == para_propa%para_Davidson%max_it)
            write(out_unitp,*) 'mod(it,para_propa%para_Davidson%num_resetH)',  &
                                mod(it,para_propa%para_Davidson%num_resetH)
          END IF
          CALL flush_perso(out_unitp)

          !----------------------------------------------------------
          !- save psi(:) on file
          IF (save_WP) THEN

            CALL sub_projec_Davidson(Ene,VecToBeIncluded,nb_diago,        &
                                     min_Ene,para_H%para_PES%min_pot,     &
                                     psi,psi0,Vec,Vec0,para_propa%para_Davidson,it,.TRUE.)

            IF (para_H%para_ReadOp%Op_Transfo) THEN
              CALL Set_ZPE_OF_ComOp(para_H%ComOp,Ene(1:count(VecToBeIncluded)),forced=.TRUE.)
            ELSE
              CALL Set_ZPE_OF_ComOp(para_H%ComOp,Ene(1:count(VecToBeIncluded)),Ene_min=min_Ene,forced=.TRUE.)
            END IF
            ZPE = para_H%ComOp%ZPE
            !CALL time_perso('projec done')

            write(out_unitp,*) 'save psi(:)',it,ndim,ndim0
            CALL flush_perso(out_unitp)

            !- check the orthogonality ------------------------
            CALL sub_MakeS_Davidson(it,psi(1:ndim),With_Grid,debug)
            !CALL time_perso('MakeS done')


            CALL sub_LCpsi_TO_psi(psi,Vec,ndim0,nb_diago)
            write(out_unitp,*) '  sub_LCpsi_TO_psi: psi done',ndim0,nb_diago
            CALL flush_perso(out_unitp)

            ! move the new vectors (nb_added_states), after the nb_diago ones
            DO i=1,nb_added_states
              psi(nb_diago+i) = psi(ndim0+i)
            END DO
            write(out_unitp,*) '  move psi done'
            CALL flush_perso(out_unitp)

            ! deallocation
            DO i=nb_diago+nb_added_states+1,ndim
              !write(out_unitp,*) 'dealloc psi(i), Hpsi(i)',i
              CALL dealloc_psi(psi(i))
            END DO
            write(out_unitp,*) '  deallocation psi done'
            CALL flush_perso(out_unitp)

            ! save the nb_diago wp
            If(MPI_id==0) THEN
              CALL sub_save_psi(psi,nb_diago,para_propa%file_WP)
              write(out_unitp,*) '  sub_save_psi: psi done'
              CALL flush_perso(out_unitp)
            ENDIF

            CALL sub_LCpsi_TO_psi(Hpsi,Vec,ndim0,nb_diago)
            write(out_unitp,*) '  sub_LCpsi_TO_psi: Hpsi done',ndim0,nb_diago
            CALL flush_perso(out_unitp)

            DO i=nb_diago+1,ndim
              !write(out_unitp,*) 'dealloc psi(i), Hpsi(i)',i
              CALL dealloc_psi(Hpsi(i))
            END DO
            write(out_unitp,*) '  deallocation Hpsi done'
            CALL flush_perso(out_unitp)



            CALL mat_id(Vec,ndim0,ndim0)

            ndim0 = nb_diago
            ndim  = nb_diago+nb_added_states


            IF (allocated(Vec0))  THEN
              CALL dealloc_NParray(Vec0,"Vec0",name_sub)
            END IF
            CALL alloc_NParray(Vec0,(/ndim0,ndim0/),"Vec0",name_sub)
            CALL mat_id(Vec0,ndim0,ndim0)

            IF (allocated(H))  THEN
              CALL dealloc_NParray(H,"H",name_sub)
            END IF
          ELSE
  !          IF (para_H%para_ReadOp%Op_Transfo) THEN
  !            CALL Set_ZPE_OF_ComOp(para_H%ComOp,Ene(1:count(VecToBeIncluded)),forced=.TRUE.)
  !          ELSE
  !            CALL Set_ZPE_OF_ComOp(para_H%ComOp,Ene(1:count(VecToBeIncluded)),Ene_min=min_Ene,forced=.TRUE.)
  !          END IF
  !          ZPE = para_H%ComOp%ZPE
  !
  !          CALL sub_save_LCpsi(psi,Vec,ndim0,nb_diago,para_propa%file_WP)
  !          CALL time_perso('save_LCpsi done')
          END IF
          !- save psi(:) on file
          !----------------------------------------------------------

          RealTime = Delta_RealTime(DavidsonTime)
          IF (RealTime < TEN) THEN
            write(out_unitp,'(a,i0,a,i0)') 'At Davidson iteration: ',it,', Delta Real time (ms): ',int(10**3*RealTime)
          ELSE
            write(out_unitp,'(a,i0,a,i0)') 'At Davidson iteration: ',it,', Delta Real time (s): ',int(RealTime)
          END IF
          CALL flush_perso(out_unitp)
        ENDIF ! for MPI_id==0

        IF (conv) EXIT

      END DO ! for it=0,para_propa%para_Davidson%max_it
      write(out_unitp,*) '--------------------------------------------------'

      !===================================================================
      !===================================================================
      ! END LOOP
      !===================================================================
      !===================================================================
      write(iunit,*) 'End Davidson ' ; CALL flush_perso(iunit)

 21   format(' Davidson: ',2(i5,1x),e10.3,i5,1x,50(1x,f18.4))

      write(out_unitp,*)
      write(out_unitp,*) '==========================================='
      write(out_unitp,*) '==========================================='
      IF (conv) THEN
        IF(MPI_id==0) write(out_unitp,*) ' Davidson has converged after ',it,' iterations'
      ELSE
        write(out_unitp,*) ' WARNNING: Davidson has NOT converged after',it,' iterations'
      END IF
      CALL flush_perso(out_unitp)

      IF (para_H%para_ReadOp%Op_Transfo) THEN
        ! The energies have to be recalculate without T(Op)
        para_H%para_ReadOp%Op_Transfo = .FALSE.
        CALL sub_MakeHPsi_Davidson(it,psi(1:nb_diago),Hpsi(1:nb_diago),Ene,0, &
                                   para_H,para_propa%para_Davidson,iunit)
        DO j=1,nb_diago
          g = Hpsi(j) - psi(j) * Ene(j)
          CALL norm2_psi(g)
          tab_norm2g(j) = sqrt(g%norm2)
          write(out_unitp,*) 'lev:',j,Ene(j)*auTOene,tab_norm2g(j)
        END DO
      END IF
      CALL file_close(Log_file)

      IF (.NOT. Hmin_OR_Hmax) THEN
        CALL Set_ZPE_OF_ComOp(para_H%ComOp,Ene(1:nb_diago),             &
                              Ene_min=min_Ene,forced=.TRUE.)
        ZPE = para_H%ComOp%ZPE
      END IF

      DO j=1,nb_diago
        psi(j)%CAvOp    = Ene(j)
        psi(j)%IndAvOp  = para_H%n_Op  ! it should be 0
        psi(j)%convAvOp = convergeEne(j) .AND. convergeResi(j)
      END DO

      IF(MPI_id==0) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'Number of Hamiltonian operations (H I psi >)',para_H%nb_OpPsi
        write(out_unitp,*)
      ENDIF

      IF (para_propa%para_Davidson%Hmax_propa) THEN
        para_propa%Hmax = -Ene(1)
        para_H%Hmax     = -Ene(1)
        write(out_unitp,*) 'Hmax (ua)  : ',para_propa%Hmax
        write(out_unitp,*) 'Hmax (cm-1): ',para_propa%Hmax*auTOene
      END IF
      IF (para_propa%para_Davidson%Hmin_propa) THEN
        para_propa%Hmin = Ene(1)
        para_H%Hmin     = Ene(1)
        write(out_unitp,*) 'Hmin (ua)  : ',para_propa%Hmin
        write(out_unitp,*) 'Hmin (cm-1): ',para_propa%Hmin*auTOene
      END IF

      IF(MPI_id==0) THEN
        write(out_unitp,*)
        write(out_unitp,*) '==========================================='
        write(out_unitp,*) '==========================================='
        CALL flush_perso(out_unitp)
      ENDIF

      !----------------------------------------------------------
      IF (allocated(Vec))  THEN
        CALL dealloc_NParray(Vec,"Vec",name_sub)
      END IF
      IF (allocated(H))  THEN
        CALL dealloc_NParray(H,"H",name_sub)
      END IF
      IF (allocated(Vec0))  THEN
        CALL dealloc_NParray(Vec0,"Vec0",name_sub)
      END IF

      IF(MPI_id==0) THEN
        DO i=1,size(Hpsi,dim=1)
          CALL dealloc_psi(Hpsi(i),delete_all=.TRUE.)
        END DO
        CALL dealloc_psi(g,delete_all=.TRUE.)

        DO i=nb_diago+1,max_diago
          CALL dealloc_psi(psi(i),delete_all=.TRUE.)
        END DO
        DO i=1,max_diago
          CALL dealloc_psi(psi0(i),delete_all=.TRUE.)
        END DO

        CALL dealloc_NParray(psi0,'psi0',name_sub)
        CALL dealloc_NParray(Hpsi,'Hpsi',name_sub)

        write(out_unitp,*) 'total memory',para_mem%mem_tot
      ENDIF

      !----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
      !----------------------------------------------------------

 END SUBROUTINE sub_LinearSystem_Jacobi

END MODULE mod_LinearSystem
