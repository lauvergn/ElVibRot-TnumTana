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
      MODULE mod_field
      USE mod_system
      IMPLICIT NONE
        integer, parameter :: max_pulse = 110

        !!@description: TODO
        !!@param: TODO
        TYPE param_field

          logical :: init0,notinit0

          integer :: nb_pola = 0
          real (kind=Rkind), pointer :: dnE(:,:) => null() ! derivative of the field
          logical :: pola_xyz(3) = (/ .FALSE., .FALSE.,.FALSE. /)

          integer :: nb_pulse = 0
          real (kind=Rkind) :: w(3,max_pulse),E0(3,max_pulse)
          real (kind=Rkind) :: t_cent(3,max_pulse),sigma(3,max_pulse)
          real (kind=Rkind) :: phase(3,max_pulse),t1(3,max_pulse)
          real (kind=Rkind) :: wmin,wmax,stepw
          character (len=Name_len) :: type,type_init_grid

          integer :: n_der    ! derivative order
          integer :: max_der  ! if n_der > max_der => dnE=0
          integer :: type_der ! derivative kind



          logical :: init_grid,init_spline,allo_grid
          TYPE (param_file) :: file
          integer :: nb_T
          real (kind=Rkind) :: Tmin,Tmax,DeltaT
          real (kind=Rkind), pointer :: grid_T(:) => null()    ! grid_T(nb_T)
          real (kind=Rkind), pointer :: grid_E(:,:) => null()  ! grid_E(nb_T,nb_pola)
          real (kind=Rkind), pointer :: grid_E2(:,:) => null() ! for the splin

        END TYPE param_field
      CONTAINS
      !====================================
      ! initialization of param_field
      !====================================
      !!@description: TODO
      !!@param: TODO
       SUBROUTINE init0_field(field,Tmax)

       TYPE (param_field) :: field
       real (kind=Rkind) :: Tmax

         field%init0          = .TRUE.
         field%notinit0       = .FALSE.

         field%nb_pola        = 0
         nullify(field%dnE)
         field%pola_xyz(:)    = .FALSE.

         field%n_der          = 0
         field%max_der        = -1
         field%type_der       = 0

         field%nb_pulse       = 1
         field%w              = ZERO
         field%E0             = ZERO
         field%phase          = ZERO
         field%t1             = Tmax
         field%t_cent         = Tmax/2.d0
         field%sigma          = ZERO
         field%wmin           = ZERO
         field%wmax           = ZERO
         field%stepw          = ZERO
         field%type           = 'cos'
         field%type_init_grid = 'cos'

         field%init_grid      = .FALSE.
         field%init_spline    = .FALSE.
         field%allo_grid      = .FALSE.
         field%file%name      = 'file_field'
         field%nb_T           = 0
         field%Tmin           = ZERO
         field%Tmax           = ZERO
         field%DeltaT         = ZERO
         nullify(field%grid_T)
         nullify(field%grid_E)
         nullify(field%grid_E2)

       END SUBROUTINE init0_field

!      ==========================================================
!
!     check if init0 has been done
!
!      ==========================================================
      !!@description: TODO
      !!@param: TODO
       SUBROUTINE check_init0_field(A,name_A,name_sub)
        TYPE (param_field), intent(in) :: A
        character (len=*), intent(in) :: name_A
        character (len=*), intent(in) :: name_sub

        IF ( (A%init0 .EQV. A%notinit0) .OR.                            &
             (A%notinit0 .AND. .NOT. A%init0) ) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) name_A,' has NOT been initiated with "init0_field"'
          write(out_unitp,*) ' CHECK the source!!!!!'
          STOP
        END IF
      END SUBROUTINE check_init0_field
      !====================================
      ! read of the param_field
      !====================================
      !!@description: TODO
      !!@param: TODO
       SUBROUTINE read_field(para_field)

       TYPE (param_field) :: para_field

       integer                :: nb_pulse
       real (kind=Rkind)      :: w(3,max_pulse),E0(3,max_pulse)
       real (kind=Rkind)      :: t_cent(3,max_pulse)
       real (kind=Rkind)      :: sigma(3,max_pulse),tt(3,max_pulse) ! sigma = tt (tt should not be used)
       real (kind=Rkind)      :: phase(3,max_pulse),t1(3,max_pulse)
       real (kind=Rkind)      :: wmax,wmin,stepw
       character (len=Name_len) :: type,type_init_grid
       integer            :: max_der
       logical            :: pola_xyz(3)
       integer            :: i_pola

       real (kind=Rkind) :: cte(0:10)

       NAMELIST /field/ cte,max_der,t1,w,phase,E0,t_cent,tt,sigma,      &
                        wmin,wmax,                                      &
                        stepw,nb_pulse,type,type_init_grid,pola_xyz


       CALL check_init0_field(para_field,'para_field','read_field')

       pola_xyz(:) = .FALSE.
       pola_xyz(3) = .TRUE.
       cte(:)    = ZERO
       nb_pulse  = 1
       w         = ZERO
       wmin      = ZERO
       wmax      = ZERO
       stepw     = ZERO
       E0        = ZERO
       phase     = ZERO
       t1        = para_field%t1
       t_cent    = para_field%t_cent
       sigma     = ZERO
       tt        = ZERO
       max_der   = -1
       type      = 'cos'
       type_init_grid = 'cos'
       read(in_unitp,field)
       write(out_unitp,field)


       IF (nb_pulse > max_pulse) THEN
         write(out_unitp,*) ' ERROR in read_field'
         write(out_unitp,*) ' the number of pulse (cos) is too large'
         write(out_unitp,*) ' nb_pulse > max_pulse',nb_pulse,max_pulse
         STOP
       END IF

       para_field%nb_pola = count(pola_xyz(:))

       IF (para_field%nb_pola == 0) THEN
         write(out_unitp,*) ' WARNING in read_field'
         write(out_unitp,*) ' no polarisation is set up'
       END IF

       IF (sum(abs(cte(:))) /= ZERO) THEN
         E0 = cte(0)
         w  = cte(1)
         t1 = cte(2)

         wmax  = cte(5)
         stepw = cte(6)
         IF (stepw > ZERO) wmin = cte(1)
       END IF

       IF (stepw > ZERO .AND. para_field%nb_pola > 1) THEN
         write(out_unitp,*) ' ERROR in read_field'
         write(out_unitp,*) ' the scan on the frequency CAN NOT be done with '
         write(out_unitp,*) ' several polarisations!'
         write(out_unitp,*) ' stepw,nb_pola',stepw,para_field%nb_pola
         STOP
       END IF

       write(out_unitp,*) 'E0',E0(:,1:nb_pulse)
       write(out_unitp,*) 'w',w(:,1:nb_pulse)
       write(out_unitp,*) 't1',t1(:,1:nb_pulse)
       write(out_unitp,*) 'IF t<t1 => E0.cos(w.t)'
       IF (stepw > ZERO) THEN
         write(out_unitp,*) 'wmin,wmax,stepw',wmin,wmax,stepw
         write(out_unitp,*) 'scan [wmin, wmin+stepw, wmin+2stepw... wmax]'
       END IF



        para_field%max_der     = max_der
        para_field%pola_xyz(:) = pola_xyz(:)

        para_field%nb_pulse    = nb_pulse
        para_field%w           = w
        para_field%E0          = E0
        para_field%phase       = phase
        para_field%t1          = t1
        IF (sum(abs(sigma)) > ZERO) THEN
          para_field%sigma       = sigma
        ELSE
          IF (sum(abs(tt)) > ZERO) THEN
            para_field%sigma       = tt
            write(out_unitp,*) ' WARNNIG in read_field'
            write(out_unitp,*) ' tt SHOULD NOT be use. Use sigma instead'
          ELSE
            para_field%sigma = ZERO
          END IF
        END IF
        para_field%t_cent      = t_cent
        para_field%wmin        = wmin
        para_field%wmax        = wmax
        para_field%stepw       = stepw
        para_field%type        = type
        para_field%type_init_grid   = type_init_grid

        para_field%init_grid   = .FALSE.
        para_field%init_spline = .FALSE.
        para_field%allo_grid   = .FALSE.
        para_field%file%name   = 'file_field'
        para_field%nb_T        = 0
        para_field%Tmin        = ZERO
        para_field%Tmax        = ZERO

       END SUBROUTINE read_field
!      ==========================================================
!
!     deallocate param_field
!
!      ==========================================================
      !!@description:  deallocate param_field
      !!@param: TODO
       SUBROUTINE dealloc_param_field(para_field)
        TYPE (param_field), intent(inout) :: para_field

        IF ( associated(para_field%dnE) ) THEN
          CALL dealloc_array(para_field%dnE,                            &
                            "para_field%dnE","dealloc_param_field")
        END IF

        IF ( associated(para_field%grid_T) ) THEN
          CALL dealloc_array(para_field%grid_T,                         &
                            "para_field%grid_T","dealloc_param_field")
        END IF

        IF ( associated(para_field%grid_E) ) THEN
          CALL dealloc_array(para_field%grid_E,                         &
                            "para_field%grid_E","dealloc_param_field")
        END IF

        IF ( associated(para_field%grid_E2) ) THEN
          CALL dealloc_array(para_field%grid_E2,                        &
                            "para_field%grid_E2","dealloc_param_field")
        END IF

      END SUBROUTINE dealloc_param_field


!===============================================
!     initialization of the field on the grid
!===============================================
      SUBROUTINE  init_field_grid(para_field,WPTmax,WPdeltaT)
      USE mod_system
      IMPLICIT NONE

      type (param_field) :: para_field

      integer :: i,j,k,npt,nio
      real (kind=Rkind) :: sigma,T,WPTmax,WPdeltaT
      real (kind=Rkind) :: ch
      real (kind=Rkind) :: T0,Tf

      !real (kind=Rkind) :: dnE,dnEcos,envelopp,dnEgauss

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING init_field_grid'
        write(out_unitp,*) 'npt,WPTmax,WPdeltaT',npt,WPTmax,WPdeltaT
        write(out_unitp,*) 'type ',para_field%type
        write(out_unitp,*) 'allo_grid ',para_field%allo_grid
        write(out_unitp,*) 'init_grid ',para_field%init_grid
        write(out_unitp,*) 'type_init_grid ',para_field%type_init_grid
        write(out_unitp,*) 'nb_pulse',para_field%nb_pulse
        write(out_unitp,*) 'E0',para_field%E0(:,1:para_field%nb_pulse)
        write(out_unitp,*) 'w',para_field%w(:,1:para_field%nb_pulse)
        write(out_unitp,*) 't1',para_field%t1(:,1:para_field%nb_pulse)
        CALL flush_perso(out_unitp)
      END IF


      IF (para_field%type .NE. 'grid') RETURN

      npt = int(WPTmax/abs(WPdeltaT))
      para_field%nb_T = npt
      write(out_unitp,*) 'Grid init_field_grid, npt: ',npt

      IF (.NOT. para_field%allo_grid) THEN
        para_field%allo_grid = .TRUE.
        CALL alloc_array(para_field%grid_T,(/npt/),                     &
                        "para_field%grid_T","init_field_grid",(/0/))
        CALL alloc_array(para_field%grid_E,(/npt,3/),                   &
                        "para_field%grid_E","init_field_grid",(/0,1/))
        para_field%grid_E(:,:) = ZERO
        para_field%grid_T(:)   = ZERO
      END IF

      IF (.NOT. para_field%init_grid) THEN
        para_field%init_grid = .TRUE.
        T = ZERO
        CALL file_open(para_field%file,nio)
        DO i=0,npt
          T = real(i,kind=Rkind)*WPdeltaT
          para_field%grid_T(i) = T
          IF (para_field%type_init_grid == 'grid' .OR.                  &
              para_field%type_init_grid == 'read' ) THEN
            read(nio,*) sigma,para_field%grid_E(i,1:3)
            IF (abs(sigma-T) > ONETENTH**4) THEN
              write(out_unitp,*) ' ERROR in init_field_grid'
              write(out_unitp,*) ' The grid point associated with T',sigma
              write(out_unitp,*) ' of the file:',para_field%file%name
              write(out_unitp,*) ' does not match with T',T
              write(out_unitp,*) ' nb_T',npt
              write(out_unitp,*) ' DeltaT',WPdeltaT
              write(out_unitp,*) 'grid_T',para_field%grid_T(1:10)
              STOP
            END IF
          ELSE IF (para_field%type_init_grid == 'cos_env' .OR.          &
              para_field%type_init_grid == 'cos' ) THEN
            DO k=1,3
              ch = ZERO
              DO j=1,para_field%nb_pulse
                ch = ch +                                               &
                   para_field%E0(k,j)*dnEcos(0,T,para_field%w(k,j)) *   &
                   envelopp(T,para_field%t1(k,j))
              END DO
              para_field%grid_E(i,k) = ch
            END DO
          ELSE IF (para_field%type_init_grid == 'cos_exp') THEN
            DO k=1,3
              ch = ZERO
              DO j=1,para_field%nb_pulse
                ch = ch +                                               &
                      dnEgauss(0,T,para_field%w(k,j),                   &
                               para_field%t_cent(k,j),                  &
                               para_field%sigma(k,j),                   &
                               para_field%phase(k,j) ) *                &
                        para_field%E0(k,j)
              END DO
              para_field%grid_E(i,k) = ch
            END DO
          ELSE IF (para_field%type_init_grid == 'cos_whitoutenv' ) THEN
            DO k=1,3
              ch = ZERO
              DO j=1,para_field%nb_pulse
                T0 = para_field%t_cent(k,j) - para_field%t1(k,j)/HALF
                Tf = para_field%t_cent(k,j) + para_field%t1(k,j)/HALF
                IF (T > t0 .AND. T<Tf) ch = ch +                        &
                   para_field%E0(k,j)*dnEcos(0,T,para_field%w(k,j))
              END DO
              para_field%grid_E(i,k) = ch
            END DO
          ELSE ! error
            write(out_unitp,*) ' ERROR in init_field_grid'
            write(out_unitp,*) ' WRONG type_init_grid: ',                       &
                             para_field%type_init_grid
            write(out_unitp,*) ' Possibilities are: '
            write(out_unitp,*) '    read, grid, cos_env, cos, gauss'
            STOP
          END IF
          IF (debug) write(out_unitp,*) ' field : ',i,T,para_field%grid_E(i,:)
          CALL flush_perso(out_unitp)
!         T = T + abs(WPdeltaT)
        END DO
        close(nio)
      END IF

      IF (para_field%max_der > 2) para_field%max_der = 2
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END init_field_grid'
      END IF
!-----------------------------------------------------------

      end subroutine  init_field_grid
!===============================================
!     save the filed on a grid
!===============================================
      SUBROUTINE  save_field_grid(para_field)
      implicit none

      type (param_field) :: para_field

      integer :: i,nio


      IF (para_field%type .NE. 'grid') RETURN


      CALL file_open(para_field%file,nio)
      DO i=0,para_field%nb_T
        write(nio,11) para_field%grid_T(i),para_field%grid_E(i,:)
 11     format(e30.20,3(1X,f30.20))
      END DO
      close(nio)


      end subroutine  save_field_grid
      SUBROUTINE  print_field_grid(para_field)
      implicit none

      type (param_field) :: para_field

      integer :: i,nio


      IF (para_field%type .NE. 'grid') THEN
        write(out_unitp,*) ' WARNING in print_field_grid'
        write(out_unitp,*) ' CANNOT print the field on a grid'
        RETURN
      ELSE
        write(out_unitp,*) ' ====================================='
        write(out_unitp,*) ' Field on a grid'
        write(out_unitp,*) ' ====================================='
        DO i=0,para_field%nb_T
        write(out_unitp,11) para_field%grid_T(i),para_field%grid_E(i,:)
 11     format(e30.20,3(1X,f30.20))
        END DO
        write(out_unitp,*) ' ====================================='
      END IF


      end subroutine  print_field_grid
!===============================================
!     find the the grid number, it, and store E
!     in grid_E(it)
!===============================================
      SUBROUTINE  EatT_TO_para_field(E,t,para_field,add)
      implicit none

      real (kind=Rkind)      :: E(3)
      real (kind=Rkind)      :: t
      type (param_field) :: para_field
      logical            :: add


      integer       :: it
      real (kind=Rkind) :: DeltaT

      IF (.NOT. para_field%allo_grid) THEN
        write(out_unitp,*) ' ERROR in EatT_TO_para_field'
        write(out_unitp,*) ' the temporal grids are not allocated!'
        STOP
      END IF


!     find the index assiociated with t
      DeltaT = para_field%grid_T(1)-para_field%grid_T(0)
      it = int(t/DeltaT)
      IF (abs(t-para_field%grid_T(it)) > 1.d-4) THEN
        write(out_unitp,*) ' ERROR in EatT_TO_para_field'
        write(out_unitp,*) ' I cannot find the grid point associated with T',t
        write(out_unitp,*) ' it,grid_T(it)',it,para_field%grid_T(it)
        write(out_unitp,*) ' nb_T',para_field%nb_T
        write(out_unitp,*) ' DeltaT',DeltaT
        write(out_unitp,*) 'grid_T',para_field%grid_T(1:10)
        write(out_unitp,*) 'grid_Ex',para_field%grid_E(1:10,1)
        write(out_unitp,*) 'grid_Ey',para_field%grid_E(1:10,2)
        write(out_unitp,*) 'grid_Ez',para_field%grid_E(1:10,3)
        STOP
      END IF

      IF (add) THEN
        para_field%grid_E(it,:) = para_field%grid_E(it,:) + E(:)
      ELSE
        para_field%grid_E(it,:) = E(:)
      END IF

      end subroutine  EatT_TO_para_field
!===============================================
!
!     calculation of dnE (general subroutine)
!     chose several kind of field :
!     cos, cos*exp, cos*env(sin), grid
!
!===============================================
      SUBROUTINE sub_dnE(dnE,n,t,para_field)
      implicit none

      type (param_field) :: para_field

      integer :: n
      integer            :: k ! polarisation
      integer            :: j
      real (kind=Rkind) :: t
      real (kind=Rkind) :: dnE(3)


!     write(out_unitp,*) 'para_field%w',para_field%w
!     write(out_unitp,*) 'para_field%E0',para_field%E0

      dnE(:) = ZERO
      IF (para_field%type .EQ. 'cos') THEN
        DO k=1,3
          DO j=1,para_field%nb_pulse
            IF (t < para_field%t1(k,j)) THEN
              dnE(k) = dnE(k) + dnEcos(n,t,para_field%w(k,j)) *         &
                           para_field%E0(k,j)
            END IF
          END DO
        END DO
      ELSE IF (para_field%type .EQ. 'cos_env') THEN
        DO k=1,3
          DO j=1,para_field%nb_pulse
            IF (t < para_field%t1(k,j)) THEN
              dnE(k) = dnE(k) +                                         &
                dnEcos_env(n,t,para_field%t1(k,j),para_field%w(k,j)) *  &
                        para_field%E0(k,j)
            END IF
          END DO
        END DO
      ELSE IF (para_field%type .EQ. 'cos_exp') THEN
        DO k=1,3
          DO j=1,para_field%nb_pulse
            IF (t < para_field%t1(k,j)) THEN
               dnE(k)=dnE(k) +                                          &
                               dnEgauss(n,t,para_field%w(k,j),          &
                                        para_field%t_cent(k,j),         &
                                        para_field%sigma(k,j),          &
                                        para_field%phase(k,j)) *        &
                        para_field%E0(k,j)
            END IF
          END DO
        END DO
      ELSE IF (para_field%type .EQ. 'grid') THEN
!       - para_field%E0 is already used in dnEgrid -----
        DO k=1,3
          dnE(k) = dnEgrid(n,t,k,para_field)
        END DO
      END IF

      end subroutine  sub_dnE
!===============================================
!     the field is on a grid
!===============================================
      FUNCTION  dnEgrid(n,t,k,para_field)
      USE mod_system
      IMPLICIT NONE

      real (kind=Rkind)      :: dnEgrid
      real (kind=Rkind)      :: E0,Ep,Epp,Em,Emm,d0E,d1E,d2E
      integer            :: n
      integer            :: k ! polarisation
      real (kind=Rkind)      :: t
      type (param_field) :: para_field


      integer, save :: it
      real (kind=Rkind), save :: DeltaT

      IF (.NOT. para_field%allo_grid) THEN
        write(out_unitp,*) ' ERROR in dnEgrid'
        write(out_unitp,*) ' the temporal grids are not allocated!'
        STOP
      END IF


!     find the index assiociated with t
      DeltaT = para_field%grid_T(1)-para_field%grid_T(0)
      it = int(t/DeltaT)
      IF (it < 0 .OR. it > para_field%nb_T) THEN
        write(out_unitp,*) ' ERROR in dnEgrid'
        write(out_unitp,*) ' it is <0 or > nb_T !!!'
        write(out_unitp,*) ' it,nb_T',it,para_field%nb_T
        write(out_unitp,*) ' T,DeltaT',t,DeltaT
        write(out_unitp,*) 'grid_T',para_field%grid_T(1:10)
        write(out_unitp,*) 'grid_E',para_field%grid_E(1:10,1)
        STOP
      END IF
      IF (abs(t-para_field%grid_T(it)) > 1.d-4) THEN
        write(out_unitp,*) ' ERROR in dnEgrid'
        write(out_unitp,*) ' I cannot find the grid point associated with T',t
        write(out_unitp,*) ' it,grid_T(it)',it,para_field%grid_T(it)
        write(out_unitp,*) ' nb_T',para_field%nb_T
        write(out_unitp,*) ' DeltaT',DeltaT
        write(out_unitp,*) 'grid_T',para_field%grid_T(1:10)
        write(out_unitp,*) 'grid_Ex',para_field%grid_E(1:10,1)
        write(out_unitp,*) 'grid_Ey',para_field%grid_E(1:10,2)
        write(out_unitp,*) 'grid_Ez',para_field%grid_E(1:10,3)
        STOP
      END IF

      d0E = para_field%grid_E(it,k)
      d1E = ZERO
      d2E = ZERO
      E0 = d0E

      IF (para_field%type_der == 1) THEN
        IF (it >0 .AND. it < para_field%nb_T) THEN
          Ep  = para_field%grid_E(it+1,k)
          Em  = para_field%grid_E(it-1,k)
          d1E = (Ep-Em) / (TWO*DeltaT)
          d2E = ( Ep+Em-E0-E0) / DeltaT**2
        ELSE IF (it == 0) THEN
          Ep  = para_field%grid_E(it+1,k)
          Epp = para_field%grid_E(it+2,k)
          d1E = (FOUR*Ep - Epp - THREE*E0) / (TWO*DeltaT)
          d2E = (Epp+E0-Ep-Ep) / DeltaT**2
        ELSE IF (para_field%type_der == 2) THEN
          Em  = para_field%grid_E(it-1,k)
          Emm = para_field%grid_E(it-2,k)
          d1E = (FOUR*Em - Emm - THREE*E0) / (-TWO*DeltaT)
          d2E = (Emm+E0-Em-Em) / DeltaT**2
        END IF
      ELSE ! type=2
        IF (para_field%DeltaT>0) THEN
          IF (it == 0) THEN
            d1E = ZERO
            d2E = ZERO
          ELSE IF (it == 1) THEN
            Em  = para_field%grid_E(it-1,k)
            d1E = (E0-Em)/DeltaT
            d2E = ZERO
          ELSE
            Em  = para_field%grid_E(it-1,k)
            Emm = para_field%grid_E(it-2,k)
            d1E = (FOUR*Em - Emm - THREE*E0) / (-TWO*DeltaT)
            d2E = (Emm+E0-Em-Em) / DeltaT**2
          END IF
        ELSE ! DeltaT<0
          IF (it == para_field%nb_T) THEN
            d1E = ZERO
            d2E = ZERO
          ELSE IF (it == para_field%nb_T-1) THEN
            Ep  = para_field%grid_E(it+1,k)
            d1E = (Ep-E0)/DeltaT
            d2E = ZERO
          ELSE
            Ep  = para_field%grid_E(it+1,k)
            Epp = para_field%grid_E(it+2,k)
            d1E = (FOUR*Ep - Epp - THREE*E0) / (TWO*DeltaT)
            d2E = (Epp+E0-Ep-Ep) / DeltaT**2
          END IF
        END IF
      END IF

      dnEgrid = ZERO
      IF (n == 0) dnEgrid = d0E
      IF (n == 1) dnEgrid = d1E
      IF (n == 2) dnEgrid = d2E

      end function  dnEgrid
!===============================================
!     sin(t * Pi/tmax)**2      * cos(w t) =
!     1/2cos(w t) -1/4 cos( 2w_env-w t) -1/4 cos( 2w_env+w t) with w_env = Pi/tmax
!===============================================
      FUNCTION  dnEcos_env(n,t,tmax,w)
      USE mod_system
      IMPLICIT NONE

      integer :: n
      real (kind=Rkind) :: t,w,ph,tmax
      real (kind=Rkind) :: dnE,dnEcos_env

      real (kind=Rkind) :: w_env,w1,w2

      integer :: k


      w_env = Pi / tmax

      w1 = w_env+w_env - w
      w2 = w_env+w_env + w

!     pour k=0
      dnE = HALF  * dnEcos(n,t,w) -                                     &
            HALF*HALF * ( dnEcos(n,t,w1) + dnEcos(n,t,w2) )

      dnEcos_env = dnE

      end function  dnEcos_env
!===============================================
!     sin(t * Pi/tmax)**2      * cos(w t)
!     1/2(1-cos(t *2Pi/tmax) ) * cos(w t)
!===============================================
      FUNCTION  dnEcos_env_old(n,t,tmax,w)
      USE mod_system
      IMPLICIT NONE

      integer :: n
      real (kind=Rkind) :: t,w,ph,tmax
      real (kind=Rkind) :: dnE,dnEcos_env_old

      real (kind=Rkind) :: w_env,a

      ! function
      real (kind=Rkind) :: combi

      integer :: k


      w_env = (pi+pi) / tmax

!     pour k=0
      dnE = HALF*(ONE-dnEcos(0,t,w_env)) *                              &
                        dnEcos(n,t,w)

      DO k=1,n
       dnE = dnE -HALF * combi(n,k) *                                   &
                        dnEcos(k,t,w_env) *                             &
                        dnEcos(n-k,t,w)
      END DO


      dnEcos_env_old = dnE

      end function  dnEcos_env_old
!===============================================
!     cos(w t)
!===============================================
      FUNCTION  dnEcos(n,t,w)
      USE mod_system
      IMPLICIT NONE

      integer :: n
      real (kind=Rkind) :: t,w,ph
      real (kind=Rkind) :: dnEcos
      real (kind=Rkind) :: dnE

      integer :: k


      k = (n+1)/2

      ph = pi/TWO
      ph = ZERO

!     write(out_unitp,*) 'n,mod(n,2),k',n,mod(n,2),k
      IF ( n .EQ. 0) THEN
        dnE = cos(w*t+ph)
!       write(out_unitp,*) 'n,cos,w**0',n

      ELSE IF ( mod(n,2) .EQ. 0) THEN
        dnE = ((-ONE)**k)*cos(w*t+ph)*w**n
!       write(out_unitp,*) 'n,k,cos,w**n*f',n,k,(-ONE)**k

      ELSE
        dnE = ((-ONE)**k)*sin(w*t+ph)*w**n
!       write(out_unitp,*) 'n,k,sin,w**n*f',n,k,(-ONE)**k

      END IF

      dnEcos = dnE

      end function  dnEcos
!==================================================================
      FUNCTION  dnEgauss(n,t,w,t_cent,sigma,ph)
      USE mod_system
      IMPLICIT NONE


      integer :: n
      real (kind=Rkind) :: t,w,ph,t_cent,sigma
      real    (kind=Rkind) :: dnEgauss
      real    (kind=Rkind) :: dnE,fac
      real (kind=Rkind) :: dt,t0,wt,sigma2,w2


      dt = t-t_cent
      wt = w*t
      sigma2 = sigma*sigma
      w2 = w*w


      IF (sigma == ZERO) THEN
        fac = ZERO
      ELSE
        dnE = exp(-(dt/sigma)**2)*cos(wt+ph)

        IF (n .EQ. 0) fac = ONE
        IF (n .EQ. 1) fac = -TWO*dt/sigma2 -w*tan(wt+ph)
        IF (n .EQ. 2) fac = (FOUR*dt**2-sigma2*(TWO+sigma2*w2)+         &
                          FOUR*dt*sigma2*w*tan(wt+ph) )/sigma**4
        IF (n .EQ. 3) fac =(-TWO*dt*(FOUR*dt**2-THREE*sigma2*           &
                            (TWO+sigma2*w2))+                           &
                     sigma2*w*(-TWELVE*dt**2+sigma2*(SIX+sigma2*w2))*   &
                     tan(wt+ph))/                                       &
                      sigma**6
        IF (n .EQ. 4) fac = w**4+TWELVE/sigma**4+TWELVE*w2/sigma2-      &
                      48._Rkind*dt**2/sigma**6-24._Rkind*w2*dt**2/sigma**4+     &
                      16._Rkind*dt**4/sigma**8+(-48._Rkind*w*dt/sigma**4-       &
                      EIGHT*w**3*dt/sigma**2+32._Rkind*dt**3*w/sigma**6)*   &
                      tan(wt+ph)
        IF (n .EQ. 5) fac=-120._Rkind*dt/sigma**6-120._Rkind*w2*dt/sigma**4-    &
                      TEN*w**4*dt/sigma2+160._Rkind*dt**3/sigma**8+         &
                      80._Rkind*w2*dt**3/sigma**6-32._Rkind*dt**5/sigma**10+    &
                      (-w**5-60._Rkind*w/sigma**4-20._Rkind*w**3/sigma2+        &
                      240._Rkind*w*dt**2/sigma**6+40._Rkind*w**3*dt**2/sigma**4-&
                      80._Rkind*w*dt**4/sigma**8)*tan(w*t+ph)
        IF (n .EQ. 6) fac=-w**6-120._Rkind/sigma**6-180._Rkind*w2/sigma**4-        &
                      30._Rkind*w**4/sigma2+720._Rkind*dt**2/sigma**8+          &
                     720._Rkind*w2*dt**2/sigma**6+60._Rkind*w**4*dt**2/sigma**4-&
                      480._Rkind*dt**4/sigma**10-240._Rkind*w2*dt**4/sigma**8+  &
                      64._Rkind*dt**6/sigma**12+(720._Rkind*w*dt/sigma**6+      &
                      240._Rkind*w**3*dt/sigma**4+TWELVE*w**5*dt/sigma2-    &
                      960._Rkind*w*dt**3/sigma**8+192._Rkind*w*dt**5/sigma**10)*&
                      tan(w*t+ph)
        IF (n .ge. 7) fac=ZERO
      END IF

      dnEgauss = dnE*fac

      end function  dnEgauss
!===============================================
!     sin(t * Pi/t1)**2
!===============================================
      FUNCTION  envelopp(t,t1)
      USE mod_system
      IMPLICIT NONE

      real (kind=Rkind) :: envelopp
      real (kind=Rkind) :: t,t1

      envelopp = sin(t/t1 * Pi) **2

      end function  envelopp

      END MODULE mod_field

