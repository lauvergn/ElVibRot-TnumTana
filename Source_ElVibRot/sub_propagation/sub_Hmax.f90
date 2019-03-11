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

      SUBROUTINE sub_Hmax(para_propa,para_H)
      USE mod_system
      USE mod_Op
      USE mod_psi_set_alloc
      USE mod_psi_SimpleOp
      USE mod_ana_psi
      USE mod_propa
      USE mod_FullPropa
      USE mod_Davidson
      IMPLICIT NONE

!----- Operator variables --------------------------------------------
      TYPE (param_Op)    :: para_H

!----- variables for the WP ----------------------------------------
      TYPE (param_psi)   :: psi,Hpsi
      TYPE (param_psi)   :: WP0,WP
      TYPE (param_propa) :: para_propa,para_propa_loc

!------ active Matrix H ------------------------------------------
      real (kind=Rkind)    :: Hinter
      complex (kind=Rkind) :: Emax

!----- quadrature points and weight -----------------------------
      real (kind=Rkind)    :: WnD

!----- working parameters --------------------------------------------
      integer           :: i_qa
      integer           :: i_h,i1_h,i,k_term
      integer           :: nio
      real (kind=Rkind) :: a
      !logical           :: relax = .TRUE.
      logical           :: relax = .FALSE.

      real (kind=Rkind) :: Qdyn(para_H%mole%nb_var)
      real (kind=Rkind) :: Qact(para_H%mole%nb_act1)
      TYPE (param_d0MatOp) :: d0MatOp
      integer              :: type_Op


      real (kind=Rkind), pointer :: Grid(:,:,:)

      integer                    :: nb_diago,max_diago
      TYPE (param_psi),pointer   :: Tab_Psi(:) => null()
      real (kind=Rkind),allocatable  :: Ene0(:)


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub = 'sub_Hmax'
      integer :: err_mem,memory
!      logical,parameter :: debug=.FALSE.
     logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------
      write(out_unitp,*) 'BEGINNING ',name_sub,' ',para_H%nb_tot
      write(out_unitp,*) 'Hmin,Hmax',para_H%Hmin,para_H%Hmax
      IF (debug) THEN

      END IF
!-----------------------------------------------------------

!     --------------------------------------------------------
!     -  Hmin ------------------------------------------------
!     --------------------------------------------------------
      IF (para_H%para_ReadOp%para_FileGrid%Save_MemGrid_done) THEN
        ! Minimal value of Veff
        para_H%Hmin = minval(para_H%OpGrid(1)%Grid(:,:,:))
      ELSE ! the grids are on a file

!       -  Hmin ------------------------------------------------
        para_H%Hmin = huge(ONE)

        IF (para_H%para_ReadOp%para_FileGrid%Type_FileGrid == 0) THEN
          type_Op = para_H%para_PES%Type_HamilOp ! H
          IF (type_Op /= 1) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) '    Type_HamilOp MUST be equal to 1 for the usual SHADA file'
            write(out_unitp,*) '    CHECK your data!!'
            STOP
          END IF

          CALL Init_d0MatOp(d0MatOp,para_H%param_TypeOp,para_H%nb_bie)

          !- loop on i_qa -----------------------------------------
          DO i_qa=1,para_H%nb_qa

            CALL sub_reading_Op(i_qa,para_H%nb_qa,d0MatOp,para_H%n_Op,  &
                                Qdyn,para_H%mole%nb_var,Qact,WnD,para_H%ComOp)

            k_term = d0MatOp%derive_term_TO_iterm(0,0)
            DO i1_h=1,para_H%nb_bie
              para_H%Hmin = min(para_H%Hmin,d0MatOp%ReVal(i1_h,i1_h,k_term))
            END DO

          END DO
          CALL dealloc_d0MatOp(d0MatOp)

        ELSE
          nullify(Grid)
          CALL alloc_array(Grid,                                        &
                          (/para_H%nb_qa,para_H%nb_bie,para_H%nb_bie/), &
                          'Grid',name_sub)

          DO k_term=1,para_H%nb_term
            IF (para_H%OpGrid(k_term)%derive_termQact(1) /= 0 .OR.      &
                para_H%OpGrid(k_term)%derive_termQact(2) /= 0) CYCLE

            IF (para_H%OpGrid(k_term)%grid_cte) THEN
               para_H%Hmin = minval(para_H%OpGrid(k_term)%Mat_cte(:,:))
            ELSE

              !$OMP critical(CRIT_sub_Hmax)
              IF (para_H%OpGrid(k_term)%file_Grid%seq) THEN   ! sequential acces file
                CALL sub_ReadSeq_Grid_iterm(Grid,para_H%OpGrid(k_term))
              ELSE  ! direct acces file
                CALL sub_ReadDir_Grid_iterm(Grid,para_H%OpGrid(k_term))
              END IF
              !$OMP end critical(CRIT_sub_Hmax)
            END IF

            EXIT ! we exit the loop because we find the potential (k_term)
          END DO
          para_H%Hmin = minval(Grid)

          CALL dealloc_array(Grid,'Grid',name_sub)

        END IF

!       --------------------------------------------------------
      END IF
      IF (debug) write(out_unitp,*) 'Hmin: ',para_H%Hmin
      CALL flush_perso(out_unitp)

!     --------------------------------------------------------
!     -  Hmax ------------------------------------------------
!     --------------------------------------------------------
      para_H%Hmax = -huge(ONE)

      CALL init_psi(psi,para_H,para_H%cplx)
      IF (debug) write(out_unitp,*) 'nb_tot',psi%nb_tot
      CALL alloc_psi(psi,BasisRep =.TRUE.)

      psi = ZERO
      IF (para_H%cplx) THEN
        psi%CvecB(psi%nb_tot) = ONE
      ELSE
        psi%RvecB(psi%nb_tot) = ONE
      END IF
      psi%symab = -1

      CALL renorm_psi(psi,BasisRep=.TRUE.)
      IF (debug) write(out_unitp,*) 'psi%symab',psi%symab

      Hpsi = psi

      CALL sub_PsiOpPsi(Emax,psi,Hpsi,para_H)
      para_H%Hmax = Real(Emax,kind=Rkind)
      IF (debug) write(out_unitp,*) 'Hmax: ',para_H%Hmax
      CALL flush_perso(out_unitp)

      CALL dealloc_psi(psi)
      CALL dealloc_psi(Hpsi)
!       --------------------------------------------------------

      write(out_unitp,*) 'non-auto : Hmin,Hmax',para_H%Hmin,para_H%Hmax
      write(out_unitp,*) 'nb_tot',para_H%nb_tot
      CALL flush_perso(out_unitp)

      IF (para_H%read_Op) THEN
        IF (para_H%cplx) THEN
          para_H%Hmax = para_H%Cmat(para_H%nb_tot,para_H%nb_tot)
          para_H%Hmin = para_H%Cmat(1,1)
        ELSE
          para_H%Hmax = para_H%Rmat(para_H%nb_tot,para_H%nb_tot)
          para_H%Hmin = para_H%Rmat(1,1)
        END IF
        para_propa%Hmax = para_H%Hmax
        para_propa%Hmin = para_H%Hmin
        write(out_unitp,*) 'read_Op : Hmin,Hmax',para_H%Hmin,para_H%Hmax
        RETURN
      END IF

      IF (para_H%spectral) THEN
        IF (para_H%cplx) THEN
          para_H%Hmax = para_H%Cdiag(para_H%ComOp%nb_vp_spec)
          para_H%Hmin = para_H%Cdiag(1)
        ELSE
          para_H%Hmax = para_H%Rdiag(para_H%ComOp%nb_vp_spec)
          para_H%Hmin = para_H%Rdiag(1)
        END IF
        para_propa%Hmax = para_H%Hmax
        para_propa%Hmin = para_H%Hmin
        write(out_unitp,*) 'spectral : Hmin,Hmax',para_H%Hmin,para_H%Hmax
        RETURN
      END IF

      IF (para_propa%auto_Hmax) THEN

        IF (para_H%cplx) relax = .TRUE. ! Because Davidson doesn't work (yet) with complex H
relax = .TRUE.
        IF (relax) THEN
          CALL sub_Auto_HmaxHmin_relax(para_propa,para_H)
        ELSE
          write(6,*) 'Davidson Hmin Hmax' ; flush(6)
          para_propa%para_Davidson%nb_WP            = 0
          para_propa%para_Davidson%lower_states     = .TRUE.
          para_propa%para_Davidson%project_WP0      = .FALSE.
          para_propa%para_Davidson%all_lower_states = .FALSE.
          para_propa%para_Davidson%max_it           = 200
          para_propa%para_Davidson%num_resetH       = para_propa%para_Davidson%max_it
          para_propa%para_Davidson%num_checkS       = para_propa%para_Davidson%max_it
          para_propa%para_Davidson%NewVec_type      = 4
          para_propa%para_Davidson%read_WP          = .FALSE.
          para_propa%para_Davidson%conv_ene         = ONETENTH**6
          para_propa%para_Davidson%conv_resi        = ONETENTH**4
          para_propa%para_Davidson%symab            = -1


          nb_diago = 0
          max_diago = para_propa%para_Davidson%max_it + 1
          nullify(Tab_Psi)
          CALL alloc_array(Tab_Psi,(/ max_diago /),"Tab_Psi",name_sub)
          CALL alloc_NParray(Ene0,(/max_diago/),"Ene0",name_sub)


          para_propa%para_Davidson%Hmin_propa       = .TRUE.
          para_propa%para_Davidson%Hmax_propa       = .FALSE.
          para_propa%para_Davidson%name_file_saveWP = 'file_WP_Hmin'

          CALL sub_propagation_Davidson(Tab_Psi,Ene0,nb_diago,max_diago,&
                                        para_H,para_propa)

          nb_diago = 0
          para_propa%para_Davidson%Hmin_propa       = .FALSE.
          para_propa%para_Davidson%Hmax_propa       = .TRUE.
          para_propa%para_Davidson%name_file_saveWP = 'file_WP_Hmax'

          CALL sub_propagation_Davidson(Tab_Psi,Ene0,nb_diago,max_diago,&
                                        para_H,para_propa)

          CALL dealloc_array(Tab_Psi,"Tab_Psi",name_sub)
          CALL dealloc_NParray(Ene0,"Ene0",name_sub)
          nullify(Tab_Psi)

          write(6,*) 'END Davidson Hmin Hmax' ; flush(6)


        END IF
      END IF

      para_propa%Hmax = para_H%Hmax
      para_propa%Hmin = para_H%Hmin
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
      END IF
      write(out_unitp,*) 'END ',name_sub
!-----------------------------------------------------------

      END SUBROUTINE sub_Hmax
      SUBROUTINE sub_Auto_HmaxHmin_relax(para_propa,para_H)
      USE mod_system
      USE mod_Op
      USE mod_psi_set_alloc
      USE mod_psi_SimpleOp
      USE mod_ana_psi
      USE mod_propa
      USE mod_FullPropa
      IMPLICIT NONE

!----- Operator variables --------------------------------------------
      TYPE (param_Op)  :: para_H

!----- variables for the WP ----------------------------------------
      TYPE (param_psi)   :: WP0,WP
      TYPE (param_propa) :: para_propa,para_propa_loc

!----- working parameters --------------------------------------------
      complex (kind=Rkind) :: Emax
      integer           :: i
      real (kind=Rkind) :: a

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub = 'sub_Auto_HmaxHmin_relax'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------
      write(out_unitp,*) 'BEGINNING ',name_sub,' ',para_H%nb_tot
      write(out_unitp,*) 'Hmin,Hmax',para_H%Hmin,para_H%Hmax
      IF (debug) THEN

      END IF
!-----------------------------------------------------------

        para_propa_loc%Hmin           = para_H%Hmin
        para_propa_loc%Hmax           = para_H%Hmax
        para_propa_loc%para_poly%Hmin = para_H%Hmin
        para_propa_loc%para_poly%Hmax = para_H%Hmax

        para_propa_loc%WPTmax         = TEN**6
        para_propa_loc%WPdeltaT       = ONE
        para_propa_loc%nb_micro       = 1

        para_propa_loc%para_poly%max_poly     = 20
        para_propa_loc%para_poly%npoly        = 20

        CALL alloc_array(para_propa_loc%work_WP,                        &
                                  (/para_propa_loc%para_poly%npoly+6/), &
                        "para_propa_loc%work_WP",name_sub,(/0/))

        CALL alloc_array(para_propa_loc%para_poly%coef_poly,(/2/),      &
                        "para_propa_loc%para_poly%coef_poly",name_sub)

        para_propa_loc%para_poly%poly_tol     = ONETENTH**6
        para_propa_loc%para_poly%DHmax        = ZERO

        para_propa_loc%write_iter        = (print_level > 1)
        para_propa_loc%n_WPecri          = 1
        para_propa_loc%WPpsi2            = .FALSE.
        para_propa_loc%WPpsi             = .FALSE.
        para_propa_loc%file_autocorr%name= 'WP_auto'


        CALL init_psi(WP,para_H,para_H%cplx)
        WP%GridRep = .TRUE.
        CALL alloc_psi(WP)

        CALL init_psi(WP0,para_H,para_H%cplx)
        WP0%BasisRep = .TRUE.
        CALL alloc_psi(WP0)

        !---- for Hmax -----------------------------------------
        para_propa_loc%name_WPpropa      = 'Emax'
        para_propa_loc%file_WP%name     = make_FileName('file_WP_Hmax')
        para_propa_loc%ana_psi%file_Psi = para_propa_loc%file_WP
        para_propa_loc%ana_psi%propa    = .TRUE.

        WP0 = ZERO
        DO i=1,WP0%nb_tot

          CALL random_number(a)
          IF (WP0%cplx) THEN
           WP0%CvecB(WP0%nb_tot+1-i) = cmplx(a-HALF,ZERO,kind=Rkind)
          ELSE
            WP0%RvecB(WP0%nb_tot+1-i) = a-HALF
          END IF
        END DO
        WP0%symab = -1
        CALL renorm_psi(WP0,BasisRep=.TRUE.)

        para_propa_loc%type_WPpropa = -3
        CALL sub_propagation3(Emax,WP0,WP,para_H,para_propa_loc)
        para_H%Hmax = Real(Emax,kind=Rkind)

        !---- END for Hmax -------------------------------------


        !---- for Hmin -----------------------------------------
        para_propa_loc%name_WPpropa       = 'Emin'
        para_propa_loc%file_WP%name       = make_FileName('file_WP_Hmin')
        para_propa_loc%ana_psi%file_Psi   = para_propa_loc%file_WP
        para_propa_loc%para_poly%poly_tol = ONETENTH**8
        para_propa_loc%WPdeltaT           = ONE

        WP0 = ZERO
        DO i=1,max(WP0%nb_tot/100,1)
          CALL random_number(a)
          IF (WP0%cplx) THEN
           WP0%CvecB(i) = cmplx(a-HALF,ZERO,kind=Rkind)
          ELSE
            WP0%RvecB(i) = a-HALF
          END IF
        END DO
        WP0%symab = -1
        CALL renorm_psi(WP0,BasisRep=.TRUE.)

        para_propa_loc%type_WPpropa = 3
        CALL sub_propagation3(Emax,WP0,WP,para_H,para_propa_loc)
        para_H%Hmin = Real(Emax,kind=Rkind)
        !---- END for Hmin -------------------------------------


        CALL dealloc_psi(WP0)
        CALL dealloc_psi(WP)

        CALL dealloc_param_propa(para_propa_loc)

        para_propa%para_poly%DHmax        = ZERO

        write(out_unitp,*) 'auto : Hmin,Hmax',para_H%Hmin,para_H%Hmax

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
      END IF
      write(out_unitp,*) 'END ',name_sub
!-----------------------------------------------------------

      END SUBROUTINE sub_Auto_HmaxHmin_relax

