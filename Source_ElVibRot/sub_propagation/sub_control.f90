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
MODULE mod_FullControl
IMPLICIT NONE

 PRIVATE
 PUBLIC :: sub_nonOpt_control,sub_Opt_control
 CONTAINS
 !================================================================
 !     control optimal : WP propagation with field (type 25)
 !
 !     Just one propagation ...
 !
 !     INTERFACE IN vib.f
 !================================================================
 SUBROUTINE sub_nonOpt_control(para_AllOp,para_propa)
      USE mod_system
      USE mod_Constant
      USE mod_psi,    ONLY : param_psi,alloc_psi,alloc_array,dealloc_array, &
                             renorm_psi,ecri_psi,sub_PsiBasisRep_TO_GridRep
      USE mod_Op
      USE mod_field
      USE mod_propa
      USE mod_FullPropa
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa

!----- for the control --------------------------------------------
      integer :: nb_WP,nb_WPba
      TYPE (param_psi),  pointer :: tab_WP0(:),tab_WPt(:),tab_WP(:)
      real (kind=Rkind), pointer :: Obj(:)
      real (kind=Rkind)              :: SObj
      integer :: it

!----- variables for H ---------------------------------------------
      TYPE (param_AllOp), target :: para_AllOp
      TYPE (param_Op), pointer   :: para_H

!----- for printing --------------------------------------------------
      logical :: print_cont


!------ working variables ---------------------------------
      integer           :: i,j,jt,iDip,iOp
      real (kind=Rkind) :: c,s,T

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
      para_H => para_AllOp%tab_Op(1)
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING sub_nonOpt_control'
        write(out_unitp,*)
        CALL write_param_Op(para_H)
      END IF
!-----------------------------------------------------------


      para_propa%Hmax = para_propa%Hmax + para_propa%para_poly%DHmax

      para_propa%para_poly%Hmin = para_propa%Hmin
      para_propa%para_poly%Hmax = para_propa%Hmax


      write(out_unitp,*) 'Tmax,DeltaT (ua)=> ',                         &
             para_propa%WPTmax,para_propa%WPdeltaT
      write(out_unitp,*) 'Tmax,DeltaT (fs)=> ',                         &
                       para_propa%WPTmax*get_Conv_au_TO_unit('t','fs'), &
                     para_propa%WPdeltaT*get_Conv_au_TO_unit('t','fs')
      write(out_unitp,*) 'Tmax,DeltaT (ps)=> ',                         &
                       para_propa%WPTmax*get_Conv_au_TO_unit('t','ps'), &
                     para_propa%WPdeltaT*get_Conv_au_TO_unit('t','ps')
      write(out_unitp,*) '... DeltaE,Emax (cm-1)',                      &
            TWO*pi/para_propa%WPTmax * get_Conv_au_TO_unit('E','cm-1'), &
          TWO*pi/para_propa%WPdeltaT * get_Conv_au_TO_unit('E','cm-1')

!     - for initialization of field variables -----------
      CALL init0_field(para_propa%para_field,para_propa%WPTmax)
      CALL read_field(para_propa%para_field)
      IF (para_propa%para_field%type == 'grid') THEN
        CALL init_field_grid(para_propa%para_field,                     &
                             para_propa%WPTmax,para_propa%WPdeltaT)
        para_propa%para_field%type_der = 2
      END IF
      CALL save_field_grid(para_propa%para_field)

!     - dipole moment on the BasisRep basis ---------------------
      iOp = 3
      DO i=iOp,iOp+2
        iDip = para_AllOp%tab_Op(i)%n_Op
        IF (para_propa%para_field%pola_xyz(iDip)) THEN
          write(out_unitp,*) 'Control with ',trim(para_AllOp%tab_Op(i)%name_Op),&
                            para_AllOp%tab_Op(i)%n_Op
        END IF
      END DO


!     - for the optimal control -------------------------

      nb_WP   = para_propa%para_control%nb_WP
      nb_WPba = para_propa%para_control%nb_WPba
      it      = 0
      nullify(tab_WP0)
      CALL alloc_array(tab_WP0,(/nb_WP/),"tab_WP0","sub_nonOpt_control")
      nullify(tab_WPt)
      CALL alloc_array(tab_WPt,(/nb_WP/),"tab_WPt","sub_nonOpt_control")
      nullify(tab_WP)
      CALL alloc_array(tab_WP,(/2*nb_WP/),"tab_WP","sub_nonOpt_control")

      DO i=1,nb_WP
        CALL init_psi(tab_WP0(i),para_H,.TRUE.)
        CALL alloc_psi(tab_WP0(i))
        CALL init_psi(tab_WPt(i),para_H,.TRUE.)
        CALL alloc_psi(tab_WPt(i))
        tab_WP0(i) = ZERO
        tab_WPt(i) = ZERO
      END DO

      DO i=1,nb_WP
        IF (para_propa%para_control%tab_WP0(i) > tab_WP0(i)%nb_tot)     &
              THEN
          write(out_unitp,*) ' ERROR in sub_nonOpt_control'
          write(out_unitp,*) ' tab_WP0(i) > nb_tot',i,                          &
                para_propa%para_control%tab_WP0(i),tab_WP0(i)%nb_tot
          STOP
        END IF
      END DO


!     --------------------------------------------------
      IF (para_propa%para_control%gate) THEN
        IF (debug) THEN
          write(out_unitp,*) 'nb_WP,nb_WPba',nb_WP,nb_WPba
          DO i=1,nb_WP
            write(out_unitp,*) ' #WP0',i,para_propa%para_control%Mgate0(i,:)
            write(out_unitp,*) ' #WPt',i,para_propa%para_control%Mgatet(i,:)
          END DO
        END IF

        DO i=1,nb_WP
        DO j=1,nb_WPba
          tab_WP0(i)%CvecB(para_propa%para_control%tab_WP0(j)) =      &
             para_propa%para_control%Mgate0(i,j)
          tab_WPt(i)%CvecB(para_propa%para_control%tab_WP0(j)) =      &
             para_propa%para_control%Mgatet(i,j)
        END DO
        END DO
      ELSE
        DO i=1,nb_WP
          IF (para_propa%para_control%tab_WPt(i) > tab_WPt(i)%nb_tot)   &
             THEN
            write(out_unitp,*) ' ERROR in sub_nonOpt_control'
            write(out_unitp,*) ' tab_WPt(i) > nb_tot',i,                        &
                para_propa%para_control%tab_WPt(i),tab_WPt(i)%nb_tot
            STOP
          END IF
          tab_WPt(i)%CvecB(para_propa%para_control%tab_WPt(i))=ONE
          tab_WP0(i)%CvecB(para_propa%para_control%tab_WP0(i))=ONE
        END DO
      END IF

      DO i=1,nb_WP

        T = ZERO
        write(out_unitp,*) 'WP0 (BasisRep)',i
        CALL renorm_psi(tab_WP0(i))
        CALL ecri_psi(T=T,psi=tab_WP0(i))

        write(out_unitp,*) 'WPt (BasisRep)',i
        CALL renorm_psi(tab_WPt(i))
        CALL ecri_psi(T=T,psi=tab_WPt(i))

      END DO
!     --------------------------------------------------

      DO i=1,nb_WP
        tab_WP(nb_WP+i) = tab_WPt(i)
        tab_WP(i)       = tab_WP0(i)
      END DO

      nullify(Obj)
      CALL alloc_array(Obj,(/ nb_WP /),"Obj","sub_nonOpt_control")
      Obj(:) = ONETENTH**2

      print_cont = para_propa%para_control%post_control

      para_propa%WPdeltaT = abs(para_propa%WPdeltaT)
      para_propa%para_field%deltaT = para_propa%WPdeltaT

      DO i=1,nb_WP
        tab_WP(i) = tab_WP0(i)
        tab_WP(nb_WP+i) = tab_WPt(i)
      END DO

      CALL sub_propagation25(tab_WP(1:nb_WP),nb_WP,                     &
                               tab_WP0,nb_WP,                           &
                               tab_WPt,nb_WP,                           &
                               print_cont,                              &
                               para_propa%para_field,.FALSE.,Obj,       &
                               para_H,para_AllOp%tab_Op(iOp:iOp+2),     &
                               para_propa)



      T = para_propa%WPTmax

      !- calculation of the objectif ----------------------
      CALL calc_fidelity(nb_WP,SObj,Obj,                                &
                           tab_WP(1:nb_WP),tab_WP(nb_WP+1:2*nb_WP),     &
                           T,it,para_propa%para_control%alpha,          &
                           .TRUE.)

      !- write WP -----------------------------------------
      IF (print_cont) THEN
        DO j=1,nb_WP

          write(out_unitp,*) 'WP (BasisRep) at T=',j,T
          CALL ecri_psi(T=T,psi=tab_WP(j))

          !for the writting of the wp (GridRep)
          write(out_unitp,*) 'WP (GridRep) at T=',j,T
          CALL sub_PsiBasisRep_TO_GridRep(tab_WP(j))
          CALL ecri_psi(T=T,psi=tab_WP(j))
        END DO ! j loop (nb_WP)
      END IF

      nullify(para_H)
      CALL dealloc_array(tab_WP0,"tab_WP0","sub_nonOpt_control")
      CALL dealloc_array(tab_WPt,"tab_WPt","sub_nonOpt_control")
      CALL dealloc_array(tab_WP, "tab_WP", "sub_nonOpt_control")
      CALL dealloc_array(Obj,    "Obj",    "sub_nonOpt_control")

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END sub_nonOpt_control'
      END IF
!-----------------------------------------------------------
      END SUBROUTINE sub_nonOpt_control
!================================================================
!     control optimal : WP propagation with field (type 25)
!
!     Le nouveau champ EST utilise au cours de la propagation
!
!     INTERFACE IN vib.f
!================================================================
SUBROUTINE sub_Opt_control(para_AllOp,para_propa)
      USE mod_system
      USE mod_Constant
      USE mod_psi,    ONLY : param_psi,alloc_psi,alloc_array,dealloc_array, &
                            renorm_psi,norm2_psi,sub_PsiBasisRep_TO_GridRep,&
                             sub_analyze_tab_psi,ecri_psi,sub_read_psi0
      USE mod_Op
      USE mod_field
      USE mod_propa
      USE mod_FullPropa
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_field) :: para_field_new

!----- for the control --------------------------------------------
      integer :: nb_WP,nb_WPba
      TYPE (param_psi), pointer  :: tab_WP0(:),tab_WPt(:),tab_WP(:)
      TYPE (param_psi), pointer  :: tab_WP_save(:)
      real (kind=Rkind), pointer :: Obj(:)
      real (kind=Rkind)              :: SObj,SObj_new,SObj1
      real (kind=Rkind)              :: epsi_obj,scal_obj
      integer                        :: it,nb_iter
      real (kind=Rkind)              :: alpha0

!----- variables for H ---------------------------------------------
      TYPE (param_AllOp), target :: para_AllOp
      TYPE (param_Op), pointer   :: para_H

!----- for printing --------------------------------------------------
      logical :: print_cont
      logical :: make_field

!------ working variables ---------------------------------
      integer           :: i,j,jt,iDip,iOp
      real (kind=Rkind) :: c,s,T

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_Opt_control'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
      para_H => para_AllOp%tab_Op(1)
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
        CALL write_param_Op(para_H)
      END IF
!-----------------------------------------------------------


      para_propa%Hmax = para_propa%Hmax + para_propa%para_poly%DHmax

      para_propa%para_poly%Hmin = para_propa%Hmin
      para_propa%para_poly%Hmax = para_propa%Hmax

      write(out_unitp,*) 'Tmax,DeltaT (ua)=> ',                         &
             para_propa%WPTmax,para_propa%WPdeltaT
      write(out_unitp,*) 'Tmax,DeltaT (fs)=> ',                         &
                       para_propa%WPTmax*get_Conv_au_TO_unit('t','fs'), &
                     para_propa%WPdeltaT*get_Conv_au_TO_unit('t','fs')
      write(out_unitp,*) 'Tmax,DeltaT (ps)=> ',                         &
                       para_propa%WPTmax*get_Conv_au_TO_unit('t','ps'), &
                     para_propa%WPdeltaT*get_Conv_au_TO_unit('t','ps')
      write(out_unitp,*) '... DeltaE,Emax (cm-1)',                      &
            TWO*pi/para_propa%WPTmax * get_Conv_au_TO_unit('E','cm-1'), &
          TWO*pi/para_propa%WPdeltaT * get_Conv_au_TO_unit('E','cm-1')


!     - for initialization of field variables -----------
      CALL init0_field(para_propa%para_field,para_propa%WPTmax)
      CALL read_field(para_propa%para_field)
      CALL init_field_grid(para_propa%para_field,                       &
                           para_propa%WPTmax,para_propa%WPdeltaT)
      CALL save_field_grid(para_propa%para_field)
      para_propa%para_field%type_der = 2
      IF (debug) CALL print_field_grid(para_propa%para_field)

      CALL init0_field(para_field_new,para_propa%WPTmax)
      para_field_new%type = 'grid'
      CALL init_field_grid(para_field_new,                              &
                           para_propa%WPTmax,para_propa%WPdeltaT)
      para_field_new%type_der = 2
      para_field_new%grid_E = para_propa%para_field%grid_E
      para_field_new%pola_xyz(:) = para_propa%para_field%pola_xyz(:)
      IF (debug) CALL print_field_grid(para_field_new)


!     - dipole moment on the BasisRep basis ---------------------
      iOp = 3
      DO i=iOp,iOp+2
        iDip = para_AllOp%tab_Op(i)%n_Op
        IF (para_propa%para_field%pola_xyz(iDip)) THEN
          write(out_unitp,*) 'Control with ',trim(para_AllOp%tab_Op(i)%name_Op),&
                            para_AllOp%tab_Op(i)%n_Op
        END IF
      END DO


!     - for the optimal control -------------------------

      nb_WP   = para_propa%para_control%nb_WP
      nb_WPba = para_propa%para_control%nb_WPba

      nullify(tab_WP0)
      CALL alloc_array(tab_WP0,    (/  nb_WP/),"tab_WP0",    name_sub)
      nullify(tab_WPt)
      CALL alloc_array(tab_WPt,    (/  nb_WP/),"tab_WPt",    name_sub)
      nullify(tab_WP)
      CALL alloc_array(tab_WP,     (/2*nb_WP/),"tab_WP",     name_sub)
      nullify(tab_WP_save)
      CALL alloc_array(tab_WP_save,(/  nb_WP/),"tab_WP_save",name_sub)


      DO i=1,nb_WP
        CALL init_psi(tab_WP0(i),para_H,.TRUE.)
        CALL alloc_psi(tab_WP0(i))
        tab_WP0(i) = ZERO

        CALL init_psi(tab_WPt(i),para_H,.TRUE.)
        CALL alloc_psi(tab_WPt(i))
        tab_WPt(i) = ZERO
      END DO

      DO i=1,nb_WP
        tab_WP(i)       = tab_WP0(i)
        tab_WP(nb_WP+i) = tab_WP0(i)
      END DO

      DO i=1,nb_WP
        IF (para_propa%para_control%tab_WP0(i) > tab_WP0(i)%nb_tot)     &
              THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' tab_WP0(i) > nb_tot',i,                          &
                para_propa%para_control%tab_WP0(i),tab_WP0(i)%nb_tot
          STOP
        END IF
      END DO


!     --------------------------------------------------
      IF (para_propa%para_WP0%New_Read_WP0) THEN
        para_propa%para_WP0%nb_WP0  = nb_WPba
        CALL sub_read_psi0(tab_WP,para_propa%para_WP0,nb_WPba)
        DO j=1,nb_WPba
          CALL norm2_psi(tab_WP(j))
          write(out_unitp,*) 'norm tab_WP(j): ',j,tab_WP(j)%norm2
        END DO

        DO i=1,nb_WP
        DO j=1,nb_WPba
          tab_WP0(i) = tab_WP0(i) +                                     &
            para_propa%para_control%Mgate0(i,j) * tab_WP(j)
          tab_WPt(i) = tab_WPt(i) +                                     &
            para_propa%para_control%Mgatet(i,j) * tab_WP(j)
        END DO
        END DO

      ELSE
        IF (para_propa%para_control%gate) THEN
          DO i=1,nb_WP
          DO j=1,nb_WPba
            tab_WP0(i)%CvecB(para_propa%para_control%tab_WP0(j)) =      &
               para_propa%para_control%Mgate0(i,j)
            tab_WPt(i)%CvecB(para_propa%para_control%tab_WP0(j)) =      &
               para_propa%para_control%Mgatet(i,j)
          END DO
          END DO
        ELSE
          DO i=1,nb_WP
            IF (para_propa%para_control%tab_WPt(i) > tab_WPt(i)%nb_tot) THEN
              write(out_unitp,*) ' ERROR in ',name_sub
              write(out_unitp,*) ' tab_WPt(i) > nb_tot',i,              &
                  para_propa%para_control%tab_WPt(i),tab_WPt(i)%nb_tot
              STOP
            END IF
            tab_WPt(i)%CvecB(para_propa%para_control%tab_WPt(i))=ONE
            tab_WP0(i)%CvecB(para_propa%para_control%tab_WP0(i))=ONE
          END DO
        END IF
      END IF

      DO i=1,nb_WP
        CALL renorm_psi(tab_WP0(i))
        CALL renorm_psi(tab_WPt(i))
      END DO

      para_propa%ana_psi%T = T
      CALL sub_analyze_tab_psi(tab_WP0,para_propa%ana_psi,adia=.FALSE.)
      CALL sub_analyze_tab_psi(tab_WPt,para_propa%ana_psi,adia=.FALSE.)

!     --------------------------------------------------

      DO i=1,nb_WP
        tab_WP(nb_WP+i) = tab_WPt(i)
        tab_WP(i)       = tab_WP0(i)
      END DO

      nullify(Obj)
      CALL alloc_array(Obj,(/ nb_WP /),"Obj",name_sub)
      Obj(:) = ONETENTH**2

!      print_cont = .TRUE.
      print_cont =.FALSE.


      IF (para_propa%para_control%max_iter > 0) THEN
!       - backward propagation of the target -------------
        para_propa%WPdeltaT = -abs(para_propa%WPdeltaT)
        para_propa%para_field%deltaT = para_propa%WPdeltaT
        para_field_new%deltaT        = para_propa%WPdeltaT

        CALL sub_propagation25(tab_WP(nb_WP+1:2*nb_WP),nb_WP,           &
                               tab_WP0,nb_WP,                           &
                               tab_WPt,nb_WP,                           &
                               print_cont,                              &
                               para_field_new,.FALSE.,Obj,              &
                               para_H,para_AllOp%tab_Op(iOp:iOp+2),     &
                               para_propa)


!       - the target or the WP are saved -------------------
        IF (para_propa%WPdeltaT < ZERO) THEN
          DO i=1,nb_WP ! the target
            tab_WP_save(i) = tab_WP(nb_WP+i)
          END DO
        ELSE ! the WP
          DO i=1,nb_WP
            tab_WP_save(i) = tab_WP(i)
          END DO
        END IF
      END IF


!------- iterations --------------------------------------------
      alpha0 = para_propa%para_control%alpha
      SObj = ZERO
      it = 0
!     print_cont = .TRUE.
      print_cont =.FALSE.
      DO ! control loop
        IF (para_propa%para_control%max_iter <1) EXIT

!       - calculation of the objectif ----------------------
        T = para_propa%WPTmax
        IF (para_propa%WPdeltaT < ZERO) T = ZERO
        CALL calc_fidelity(nb_WP,SObj_new,Obj,                          &
                           tab_WP(1:nb_WP),tab_WP(nb_WP+1:2*nb_WP),     &
                           T,it,para_propa%para_control%alpha,          &
                           .TRUE.)

        CALL flush_perso(out_unitp)

        epsi_obj = ONETENTH**2
        IF (.NOT. para_propa%para_control%krotov) epsi_obj = FIVE
        scal_obj = 1.2_Rkind
        IF ((SObj_new - SObj) < -epsi_obj) THEN
          para_propa%para_control%alpha =                               &
                             para_propa%para_control%alpha*scal_obj
          para_field_new%grid_E = para_propa%para_field%grid_E

!         - the target or the WP are restored ----------------
          IF (para_propa%WPdeltaT < ZERO) THEN
            DO i=1,nb_WP
              tab_WP(i) = tab_WP_save(i)
              tab_WP(nb_WP+i) = tab_WPt(i)
            END DO
          ELSE
            DO i=1,nb_WP
              tab_WP(nb_WP+i) = tab_WP_save(i)
              tab_WP(i) = tab_WP0(i)
            END DO
          END IF


        ELSE
!         - the target or the WP are saved -------------------
          IF (para_propa%WPdeltaT < ZERO) THEN
            DO i=1,nb_WP ! the target
              tab_WP_save(i) = tab_WP(nb_WP+i)
              tab_WP(i)      = tab_WP0(i)
            END DO
          ELSE ! the WP
            DO i=1,nb_WP
              tab_WP_save(i)  = tab_WP(i)
              tab_WP(nb_WP+i) = tab_WPt(i)
            END DO
          END IF

          para_propa%para_field%grid_E = para_field_new%grid_E
          SObj = SObj_new
          CALL save_field_grid(para_propa%para_field)
          para_propa%WPdeltaT = -para_propa%WPdeltaT
          para_propa%para_field%deltaT = para_propa%WPdeltaT
          para_field_new%deltaT        = para_propa%WPdeltaT
          it = it+1
        END IF

        IF (para_propa%para_control%alpha >                             &
                  para_propa%para_control%Max_alpha) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' alpha is too large',                             &
                 para_propa%para_control%alpha
          STOP
        END IF

        IF (SObj_new > para_propa%para_control%conv                     &
               .OR. it > para_propa%para_control%max_iter) EXIT

!       - WP0 or WPt in WP -------------------------------
        IF (print_cont) THEN
          write(out_unitp,*) 'it,DeltaT',it,para_propa%WPdeltaT
          T = para_propa%WPTmax
          IF (para_propa%WPdeltaT < ZERO) T = ZERO
          DO j=1,2*nb_WP
            write(out_unitp,*) 'WP (BasisRep) ',j
            CALL ecri_psi(T=T,psi=tab_WP(j))
          END DO
        END IF
        CALL sub_propagation25(tab_WP(1:2*nb_WP),2*nb_WP,               &
                               tab_WP0,nb_WP,                           &
                               tab_WPt,nb_WP,                           &
                               print_cont,                              &
                               para_field_new,.TRUE.,Obj,               &
                               para_H,para_AllOp%tab_Op(iOp:iOp+2),     &
                               para_propa)

        IF (print_cont) THEN
          DO j=1,2*nb_WP
            write(out_unitp,*) 'WP (BasisRep) ',j
            CALL ecri_psi(T=T,psi=tab_WP(j))
          END DO
        END IF


      END DO ! control iteration
!------- iterations --------------------------------------------

!     - forward propagation of the WP (after CV) ------
      IF (para_propa%para_control%max_iter <1 .OR.                      &
            para_propa%para_control%post_control) THEN

        print_cont = para_propa%para_control%post_control

        para_propa%WPdeltaT = abs(para_propa%WPdeltaT)
        para_propa%para_field%deltaT = para_propa%WPdeltaT
        para_field_new%deltaT        = para_propa%WPdeltaT

        DO i=1,nb_WP
          tab_WP(i) = tab_WP0(i)
          tab_WP(nb_WP+i) = tab_WPt(i)
        END DO

        CALL sub_propagation25(tab_WP(1:nb_WP),nb_WP,                   &
                               tab_WP0,nb_WP,                           &
                               tab_WPt,nb_WP,                           &
                               print_cont,                              &
                               para_field_new,.FALSE.,Obj,              &
                               para_H,para_AllOp%tab_Op(iOp:iOp+2),     &
                               para_propa)



        T = para_propa%WPTmax
        IF (para_propa%WPdeltaT < ZERO) T = ZERO


!       - calculation of the objectif ----------------------
        CALL calc_fidelity(nb_WP,SObj_new,Obj,                          &
                           tab_WP(1:nb_WP),tab_WP(nb_WP+1:2*nb_WP),     &
                           T,it,para_propa%para_control%alpha,          &
                           .TRUE.)

!       - write WP -----------------------------------------
        IF (print_cont) THEN
          DO j=1,nb_WP

            write(out_unitp,*) 'WP (BasisRep) at T=',j,T
            CALL ecri_psi(T=T,psi=tab_WP(j))

!           for the writting of the wp (GridRep)
            write(out_unitp,*) 'WP (GridRep) at T=',j,T
            CALL sub_PsiBasisRep_TO_GridRep(tab_WP(j))
            CALL ecri_psi(T=T,psi=tab_WP(j))
          END DO ! j loop (nb_WP)
        END IF
      END IF


      ! - deallocate WP ------------------------------------
      CALL dealloc_array(tab_WP0,    "tab_WP0",    name_sub)
      CALL dealloc_array(tab_WPt,    "tab_WPt",    name_sub)
      CALL dealloc_array(tab_WP,     "tab_WP",     name_sub)
      CALL dealloc_array(tab_WP_save,"tab_WP_save",name_sub)

      CALL dealloc_array(Obj,        "Obj",        name_sub)

      nullify(para_H)
      CALL dealloc_param_field(para_field_new)
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------
      END SUBROUTINE sub_Opt_control
!================================================================
!
!    25 : nOD propagation with
!         a time dependant pulse in Hamiltonian (W(t))
!         H is the square matrix (dimension n)
!         Hmin and Hmax are the parameter to scale H
!         WP and WPt are associated wth the inital WP and the target WP
!
!    for the control
!================================================================
      SUBROUTINE sub_propagation25(WP,nb_WP,WP0,nb_WP0,WPt,nb_WPt,      &
                                   print_Op,                            &
                                   para_field_new,make_field,Obj0,      &
                                   para_H,para_Dip,para_propa)
      USE mod_system
      USE mod_Constant
      USE mod_psi,    ONLY : param_psi,renorm_psi,ecri_psi
      USE mod_Op
      USE mod_propa
      USE mod_march
      USE mod_field
      IMPLICIT NONE

!----- Operator : H and Dip(:) ---------------------------------------
      TYPE (param_Op)   :: para_H
      TYPE (param_Op)   :: para_Dip(3)

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa

      integer            :: nb_WP,nb_WP0,nb_WPt
      TYPE (param_psi)   :: WP(nb_WP),WP0(nb_WP0),WPt(nb_WPt)

!----- for the fidelity or objectif --------------------------------
      real (kind=Rkind)  :: Obj0(nb_WPt)
      real (kind=Rkind) :: SObj,Obj(nb_WPt),alpha

!----- for the field --------------------------------------------------
      TYPE (param_field) :: para_field_new
      logical :: make_field
!----- for the field --------------------------------------------------

!----- for printing --------------------------------------------------
      logical :: print_Op, print_WP

!------ working parameters --------------------------------
      integer       :: i,j,it,it_max
      real (kind=Rkind) :: T      ! time

      integer  ::   nioWP
      logical :: BasisRep,GridRep


!----- for debuging --------------------------------------------------
       logical, parameter :: debug=.FALSE.
       !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING sub_propagation25'
        write(out_unitp,*) 'Tmax,deltaT',para_propa%WPTmax,para_propa%WPdeltaT
        write(out_unitp,*) 'Hmin,Hmax',para_propa%para_poly%Hmin,               &
                                para_propa%para_poly%Hmax
        write(out_unitp,*)
        write(out_unitp,*) 'nb_ba,nb_qa',WP(1)%nb_ba,WP(1)%nb_qa
        write(out_unitp,*) 'nb_bi',WP(1)%nb_bi
        write(out_unitp,*)

        DO j=1,nb_WP
          CALL renorm_psi(WP(j),GridRep=.FALSE.,BasisRep=.TRUE.)

          write(out_unitp,*) 'WP0 (BasisRep)',j
          CALL ecri_psi(T=ZERO,psi=WP(j),                               &
                        ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.,      &
                        ecri_psi2=.FALSE.)

        END DO ! j loop (nb_WP)

      END IF
!-----------------------------------------------------------

      BasisRep = WP0(1)%BasisRep
      GridRep  = WP0(1)%GridRep

!-----------------------------------------------------------
      IF (print_Op .OR. debug) THEN
        write(out_unitp,*) ' Propagation ',para_propa%name_WPpropa
        IF (para_propa%para_field%pola_xyz(1)) write(out_unitp,*) 'with Dipx'
        IF (para_propa%para_field%pola_xyz(2)) write(out_unitp,*) 'with Dipy'
        IF (para_propa%para_field%pola_xyz(3)) write(out_unitp,*) 'with Dipz'
      END IF

!     - parameters for poly (cheby and nOD) ... ------------
      CALL initialisation1_poly(para_propa%para_poly,                   &
                                para_propa%WPdeltaT,                    &
                                para_propa%type_WPpropa)

!     - scaling of H ---------------------------------------
      para_H%scaled = .TRUE.
      para_H%E0     = para_propa%para_poly%E0
      para_H%Esc    = para_propa%para_poly%Esc

      para_propa%march_error = .FALSE.
!-----------------------------------------------------------
!------- propagation loop ---------------------------------
      T = ZERO
      IF (para_propa%WPdeltaT < 0) T = para_propa%WPTmax
      it     = 0
      it_max = para_propa%WPTmax/abs(para_propa%WPdeltaT)
      CALL file_open(para_propa%file_WP,nioWP)
      para_propa%test_max_norm = .FALSE.
!     ------------------------------------------------------
      DO

!       - Build of the field ------------------------------
        IF (make_field)                                                 &
           CALL build_field(T,WP,nb_WP,para_field_new,Obj0,             &
                            para_Dip,para_propa)

        print_WP = (print_Op .AND. mod(it,para_propa%n_WPecri) == 0)

        IF (print_Op .OR. debug) THEN
          CALL sub_analyze_WP_OpWP(T,WP,nb_WP,para_H,para_propa,             &
                                   para_field=para_propa%para_field)
          IF (nb_WP == nb_WPt) THEN
            alpha=ZERO
            CALL calc_fidelity(nb_WPt,SObj,Obj,WP,WPt,T,0,alpha,.TRUE.)
          END IF

        END IF
        IF (it == it_max .OR. para_propa%march_error) EXIT  ! exit the propagation loop

        IF (make_field .AND. para_propa%WPdeltaT > 0) THEN

          CALL march_gene(T,WP(1:nb_WP/2),WP(1:nb_WP/2),nb_WP/2,print_WP,&
                          para_H,para_propa,                            &
                          para_Dip,para_field_new)
          CALL march_gene(T,WP(nb_WP/2+1:nb_WP),WP(nb_WP/2+1:nb_WP),nb_WP/2,print_WP,&
                          para_H,para_propa,                            &
                          para_Dip,para_propa%para_field)

        ELSE IF (make_field .AND. para_propa%WPdeltaT < 0) THEN

          CALL march_gene(T,WP(1:nb_WP/2),WP(1:nb_WP/2),nb_WP/2,print_WP,&
                          para_H,para_propa,                            &
                          para_Dip,para_propa%para_field)
          CALL march_gene(T,WP(nb_WP/2+1:nb_WP),WP(nb_WP/2+1:nb_WP),nb_WP/2,print_WP,&
                          para_H,para_propa,                            &
                          para_Dip,para_field_new)
        ELSE

          CALL march_gene(T,WP(:),WP(:),nb_WP,print_WP,                 &
                          para_H,para_propa,                            &
                          para_Dip,para_propa%para_field)
        END IF

         it = it + 1
         T = real(it,kind=Rkind) * para_propa%WPdeltaT
         IF (para_propa%WPdeltaT < ZERO) T = T + para_propa%WPTmax

      END DO ! loop on the Time iteration
      close(nioWP)
!-----------------------------------------------------------
!-----------------------------------------------------------



!----------------------------------------------------------
!     - write the final WP---------------------------------
      IF (print_Op .OR. debug) THEN
        DO j=1,nb_WP
          CALL renorm_psi(WP(j),GridRep=.FALSE.,BasisRep=.TRUE.)

          write(out_unitp,*) 'WP (BasisRep)',j,' at T=',T
          CALL ecri_psi(T=T,psi=WP(j),                                  &
                        ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.,      &
                        ecri_psi2=.FALSE.)

        END DO ! j loop (nb_WP)
      END IF
      CALL flush_perso(out_unitp)

!----------------------------------------------------------

      IF (para_propa%march_error) THEN
        write(out_unitp,*) ' ERROR in sub_propagation25'
        write(out_unitp,*) ' March: norm too large, no convergence...'
        IF (para_propa%test_max_norm)                                   &
         write(out_unitp,*) ' the norm2 is too large! ',WP(:)%norm2
        STOP
      END IF



!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END sub_propagation25'
      END IF
!----------------------------------------------------------


      end subroutine sub_propagation25

!================================================================
!    construct the new field
!    for a given value of T
!================================================================
      SUBROUTINE build_field(T,WP,nb_WP,para_field_new,Obj0,            &
                             para_Dip,para_propa)
      USE mod_system
      USE mod_psi,    ONLY : param_psi,dealloc_psi,Overlap_psi1_psi2
      USE mod_Op
      USE mod_field
      USE mod_propa
      IMPLICIT NONE

!----- variables pour la namelist minimum ----------------------------
      TYPE (param_Op)   :: para_Dip(3)

!----- variables for the WP propagation ----------------------------
      TYPE (param_propa) :: para_propa
      TYPE (param_field) :: para_field_new
      real (kind=Rkind)      :: E_new(3)

      integer            :: nb_WP
      TYPE (param_psi)   :: WP(nb_WP)
      real (kind=Rkind)      :: Obj0(nb_WP/2)
      real (kind=Rkind)      :: Obj0_w(nb_WP/2)


!------ working parameters --------------------------------
      TYPE (param_psi)     :: w2
      integer              :: i,j,k,jt,ip
      real (kind=Rkind)    :: T      ! time
      complex (kind=Rkind) :: S,avMu
      real (kind=Rkind)    :: alpha_j,Obj_min
      logical          :: add

      integer  ::   nioWP

!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING build_field'
        write(out_unitp,*)
        write(out_unitp,*) para_field_new%type
        write(out_unitp,*) para_field_new%allo_grid
        CALL flush_perso(out_unitp)
       END IF
!-----------------------------------------------------------

!      - Build of the field ------------------------------
       IF (nb_WP > 1 .AND. mod(nb_WP,2) == 0) THEN
         IF (para_propa%para_control%Obj_TO_alpha) THEN
           Obj0_w(:) = max(ONETENTH**4,Obj0(:))
           Obj_min = minval(Obj0_w(:))
         ELSE
           Obj0_w(:) = ONE
           Obj_min = minval(Obj0_w(:))
         END IF

         E_new(:) = ZERO
         DO j=1,nb_WP/2
           jt = j+nb_WP/2
           IF (para_propa%para_control%Turinici) THEN
             S = ONE
           ELSE
             CALL Overlap_psi1_psi2(S,WP(j),WP(jt))
           END IF
           alpha_j = Obj0_w(j)/Obj_min * para_propa%para_control%alpha
           DO ip=1,3
             IF (.NOT. para_propa%para_field%pola_xyz(ip)) CYCLE
             CALL sub_OpPsi(WP(j),w2,para_Dip(ip))
             CALL Overlap_psi1_psi2(avMu,WP(jt),w2)
             E_new(ip) = E_new(ip) - aimag(S*avMu)/alpha_j
           END DO
!          write(out_unitp,*) 'T,S,avMu',j,T,S,avMu
         END DO
         DO ip=1,3
           IF (.NOT. para_propa%para_field%pola_xyz(ip)) CYCLE
           IF (para_propa%para_control%envelopp) THEN
             E_new(ip) = E_new(ip) *                                    &
                   envelopp(T,para_propa%para_control%Tenvelopp)
           END IF
         END DO
         add = para_propa%para_control%Krotov
!        write(out_unitp,*) 'T,S,E_new,add,alpha_j',T,S,E_new(:),add,alpha_j
         CALL EatT_TO_para_field(E_new,T,para_field_new,add)
       END IF


       IF (T == para_propa%para_field%grid_T(1) .AND.                   &
             para_propa%para_control%Obj_TO_alpha) THEN
         write(out_unitp,*) 'T Log(Obj/min_Obj)',T,int(log10(Obj0_w(:)/Obj_min))
       END IF

       CALL dealloc_psi(w2)

!----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'END build_field'
       END IF
!----------------------------------------------------------


      end subroutine build_field
!================================================================
!     Fidelity calculation
!
!================================================================
      SUBROUTINE calc_fidelity(nb_WP,SObj,Obj,tab_WP,tab_WPt,           &
                               T,it,alpha,print_fid)
      USE mod_system
      USE mod_psi,    ONLY : param_psi,Overlap_psi1_psi2
      IMPLICIT NONE

!-------------------------------------------------------------------------

!----- for the control --------------------------------------------
      integer              :: nb_WP
      TYPE (param_psi)     :: tab_WP(nb_WP),tab_WPt(nb_WP)
      real (kind=Rkind)    :: Obj(nb_WP)
      real (kind=Rkind)    :: SObj

!----- for printing --------------------------------------------------
      logical :: print_fid


!------ working variables ---------------------------------
      integer  :: j,it
      real (kind=Rkind) :: alpha,T
      complex (kind=Rkind) :: Over

!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING calc_fidelity'
        write(out_unitp,*) 'nb_WP',nb_WP
      END IF
!-----------------------------------------------------------


!     - calculation of the objectif ----------------------
      DO j=1,nb_WP
        CALL Overlap_psi1_psi2(Over,tab_WP(j),tab_WPt(j))
        Obj(j) = abs(Over)**2
      END DO

      SObj = sum(Obj(1:nb_WP))/nb_WP
      IF (print_fid .OR. debug)                                         &
           write(out_unitp,11) ' it,T,alpha,J:',it,T,alpha,SObj,Obj(1:nb_WP)
 11       format(a,i4,1x,f12.2,1x,f8.2,20(1x,f6.4))

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING calc_fidelity'
      END IF
!-----------------------------------------------------------
      end subroutine calc_fidelity
END MODULE mod_FullControl

