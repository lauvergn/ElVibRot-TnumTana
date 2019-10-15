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
      MODULE mod_OpGrid

      USE mod_system
      USE mod_basis_BtoG_GtoB_SGType4
      IMPLICIT NONE

       TYPE param_FileGrid

          logical                    :: Read_FileGrid      = .FALSE.   ! Read the grid from a file

          logical                    :: Save_FileGrid      = .TRUE.   ! Save the grid in a file
          logical                    :: Save_FileGrid_done = .FALSE.  ! T, if the grid is save in a file
          logical                    :: Formatted_FileGrid = .TRUE.   ! format of the file

          integer                    :: Type_FileGrid      = 0        ! 0 in normal SH_HADA file
                                                                      ! 1 unformatted sequential acces
                                                                      ! 2 unformatted direct acces
                                                                      ! 4 unformatted sequential acces (for SG4)

          logical                    :: Keep_FileGrid      = .TRUE.   ! Keep the file

          logical                    :: Save_MemGrid       = .FALSE.  ! Save the grid in memory
          logical                    :: Save_MemGrid_done  = .FALSE.  ! T, if the grid is save in memory

          character (len=line_len)   :: Base_FileName_Grid = "SH_HADA" ! base name of grid file

          !Remarks:
          !         direct=0    => Make_Mat=T, SaveFile_Grid=T, SaveMem_Grid=F
          !         direct=1    => Make_Mat=F, SaveFile_Grid=T, SaveMem_Grid=T
          !         direct=2    => Make_Mat=F, SaveFile_Grid=F, SaveMem_Grid=T
          !         direct=3    => Make_Mat=F, SaveFile_Grid=T, SaveMem_Grid=F (for huge grid, like cHAC)

           logical                  :: Test_Grid       = .TRUE.    ! test calculation on one active grid point (Qdyn0)
           logical                  :: Restart_Grid    = .FALSE.   ! if t => restarting the coupled adiabatic calculations
           integer                  :: First_GridPoint = 0         ! if =0 calculation for all active grid points
           integer                  :: Last_GridPoint  = 0         ! if =0 calculation for all active grid points

        END TYPE param_FileGrid

        TYPE param_OpGrid
          integer                    :: nb_qa              =  0
          integer                    :: nb_bie             =  0

          integer, pointer           :: derive_termQact(:) => null() ! (2)
          integer, pointer           :: derive_termQdyn(:) => null() ! (2)
          real (kind=Rkind), pointer :: Grid(:,:,:)        => null() ! (nb_qa,nb_bie,nb_bie)
          real (kind=Rkind), pointer :: Mat_cte(:,:)       => null() ! (nb_bie,nb_bie). Used if grid_cte=.true.
          logical                    :: grid_zero          =  .FALSE.! TRUE if the grid is ZERO
          logical                    :: grid_cte           =  .FALSE.! TRUE if the grid is constante
          logical                    :: ana_grid           =  .FALSE.! ???
          logical                    :: cplx               =  .FALSE.! TRUE if it the imaginary part of the complex grid
          logical                    :: alloc_Grid         =  .FALSE.
          logical                    :: Grid_done          =  .FALSE.

          integer                    :: iq_min             =  0
          integer                    :: iq_max             =  0
          real (kind=Rkind)          :: Op_min             =  huge(ONE)
          real (kind=Rkind)          :: Op_max             = -huge(ONE)

          TYPE(Type_SmolyakRep)      :: SRep                         ! Smolyak Rep (SG4)

          TYPE (param_file)          :: file_Grid                    ! file of the grid
          TYPE (param_FileGrid)      :: para_FileGrid

        END TYPE param_OpGrid

      INTERFACE assignment (=)
        MODULE PROCEDURE para_FileGrid2TOpara_FileGrid1
      END INTERFACE

      INTERFACE alloc_array
        MODULE PROCEDURE alloc_array_OF_OpGriddim1
      END INTERFACE

      INTERFACE dealloc_array
        MODULE PROCEDURE dealloc_array_OF_OpGriddim1
      END INTERFACE

      CONTAINS

      SUBROUTINE para_FileGrid2TOpara_FileGrid1(para_FileGrid1,para_FileGrid2)

      TYPE (param_FileGrid), intent(inout) :: para_FileGrid1
      TYPE (param_FileGrid), intent(in)    :: para_FileGrid2

      para_FileGrid1%Read_FileGrid      = para_FileGrid2%Read_FileGrid

      para_FileGrid1%Save_FileGrid      = para_FileGrid2%Save_FileGrid
      para_FileGrid1%Save_FileGrid_done = para_FileGrid2%Save_FileGrid_done
      para_FileGrid1%Formatted_FileGrid = para_FileGrid2%Formatted_FileGrid
      para_FileGrid1%Type_FileGrid      = para_FileGrid2%Type_FileGrid
      para_FileGrid1%Keep_FileGrid      = para_FileGrid2%Keep_FileGrid
      para_FileGrid1%Base_FileName_Grid = para_FileGrid2%Base_FileName_Grid

      para_FileGrid1%Save_MemGrid       = para_FileGrid2%Save_MemGrid
      para_FileGrid1%Save_MemGrid_done  = para_FileGrid2%Save_MemGrid_done

      para_FileGrid1%Test_Grid          = para_FileGrid2%Test_Grid
      para_FileGrid1%Restart_Grid       = para_FileGrid2%Restart_Grid
      para_FileGrid1%First_GridPoint    = para_FileGrid2%First_GridPoint
      para_FileGrid1%Last_GridPoint     = para_FileGrid2%Last_GridPoint

      END SUBROUTINE para_FileGrid2TOpara_FileGrid1

      SUBROUTINE init_FileGrid(para_FileGrid,Type_FileGrid,             &
                               Read_FileGrid,Restart_Grid,Test_Grid,    &
                               First_GridPoint,Last_GridPoint,          &
                               Save_FileGrid,Formatted_FileGrid,        &
                               Keep_FileGrid,Save_MemGrid,              &
                               Base_FileName_Grid)

      TYPE (param_FileGrid), intent(inout) :: para_FileGrid
      logical, optional :: Save_FileGrid,Keep_FileGrid,Save_MemGrid,Formatted_FileGrid
      logical, optional :: Read_FileGrid,Restart_Grid,Test_Grid

      integer, optional :: Type_FileGrid
      integer, optional :: First_GridPoint,Last_GridPoint

      character (len=*), optional  :: Base_FileName_Grid


      IF (present(Read_FileGrid)) THEN
        para_FileGrid%Read_FileGrid    = Read_FileGrid
      ELSE
        para_FileGrid%Read_FileGrid    = .FALSE.   ! in SH_HADA file
      END IF

      IF (present(Restart_Grid)) THEN
        para_FileGrid%Restart_Grid    = Restart_Grid
      ELSE
        para_FileGrid%Restart_Grid    = .FALSE.   ! in SH_HADA file
      END IF
      IF (present(Test_Grid)) THEN
        para_FileGrid%Test_Grid    = Test_Grid
      ELSE
        para_FileGrid%Test_Grid    = .TRUE.   ! in SH_HADA file
      END IF
      IF (present(First_GridPoint)) THEN
        para_FileGrid%First_GridPoint    = First_GridPoint
      ELSE
        para_FileGrid%First_GridPoint    = 0
      END IF
      IF (present(Last_GridPoint)) THEN
        para_FileGrid%Last_GridPoint    = Last_GridPoint
      ELSE
        para_FileGrid%Last_GridPoint    = 0
      END IF

      IF (present(Save_FileGrid)) THEN
        para_FileGrid%Save_FileGrid    = Save_FileGrid
      ELSE
        para_FileGrid%Save_FileGrid    = .TRUE.   ! in SH_HADA file
      END IF
      para_FileGrid%Save_FileGrid_done = .FALSE.  ! in SH_HADA file

      IF (present(Formatted_FileGrid)) THEN
        para_FileGrid%Formatted_FileGrid    = Formatted_FileGrid
      ELSE
        para_FileGrid%Formatted_FileGrid    = .TRUE.
      END IF

      IF (present(Type_FileGrid)) THEN
        para_FileGrid%Type_FileGrid    = Type_FileGrid
      ELSE
        para_FileGrid%Type_FileGrid    = 0        ! 0 in normal SH_HADA file
                                                  ! 1 unformatted sequential acces
                                                  ! 2 unformatted direct acces
                                                  ! 4 unformatted sequential acces (for SG4)
      END IF
      IF (para_FileGrid%Type_FileGrid /= 0) THEN
        para_FileGrid%Formatted_FileGrid    = .FALSE.
      END IF

      IF (present(Keep_FileGrid)) THEN
        para_FileGrid%Keep_FileGrid    = Keep_FileGrid
      ELSE
        para_FileGrid%Keep_FileGrid    = .TRUE.   ! Keep the SH_HADA file
      END IF

      IF (present(Save_MemGrid)) THEN
        para_FileGrid%Save_MemGrid     = Save_MemGrid
      ELSE
        para_FileGrid%Save_MemGrid     = .FALSE.  ! in param_OpGrid stucture
      END IF
      para_FileGrid%Save_MemGrid_done  = .FALSE.  ! in param_OpGrid stucture


      IF (present(Base_FileName_Grid)) THEN
        para_FileGrid%Base_FileName_Grid    = make_FileName(Base_FileName_Grid)
      ELSE
        para_FileGrid%Base_FileName_Grid    = make_FileName("SH_HADA")
      END IF

      END SUBROUTINE init_FileGrid
      SUBROUTINE Write_FileGrid(para_FileGrid)

      TYPE (param_FileGrid), intent(in) :: para_FileGrid

      write(out_unitp,*) 'BEGINNING Write_FileGrid'

      write(out_unitp,*) 'Save_FileGrid      ',para_FileGrid%Save_FileGrid
      write(out_unitp,*) 'Save_FileGrid_done ',para_FileGrid%Save_FileGrid_done
      write(out_unitp,*) 'Read_FileGrid      ',para_FileGrid%Read_FileGrid
      write(out_unitp,*) 'Formatted_FileGrid ',para_FileGrid%Formatted_FileGrid

      write(out_unitp,*) 'Type_FileGrid      ',para_FileGrid%Type_FileGrid
      write(out_unitp,*) 'Keep_FileGrid      ',para_FileGrid%Keep_FileGrid

      write(out_unitp,*) 'Save_MemGrid       ',para_FileGrid%Save_MemGrid
      write(out_unitp,*) 'Save_MemGrid_done  ',para_FileGrid%Save_MemGrid_done

      write(out_unitp,*) 'Base_FileName_Grid: ',trim(adjustl(para_FileGrid%Base_FileName_Grid))

      write(out_unitp,*) 'Test_Grid           ',para_FileGrid%Test_Grid
      write(out_unitp,*) 'Restart_Grid        ',para_FileGrid%Restart_Grid
      write(out_unitp,*) 'First_GridPoint     ',para_FileGrid%First_GridPoint
      write(out_unitp,*) 'Last_GridPoint      ',para_FileGrid%Last_GridPoint

      write(out_unitp,*) 'END Write_FileGrid'

      END SUBROUTINE Write_FileGrid


      SUBROUTINE alloc_OpGrid(OpGrid,nb_qa,nb_bie,                      &
                              derive_termQact,derive_termQdyn,SmolyakRep,nb_SG,info)

          TYPE (param_OpGrid), intent(inout) :: OpGrid
          integer,             intent(in)    :: nb_qa,nb_bie
          integer,             intent(in)    :: derive_termQact(2)
          integer,             intent(in)    :: derive_termQdyn(2)
          logical,             intent(in)    :: SmolyakRep
          integer,             intent(in)    :: nb_SG
          character (len=*),   intent(in)    :: info

          character (len=Name_longlen) :: info2
          integer :: err
          integer :: err_mem,memory
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='alloc_OpGrid'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'nb_qa,nb_bie',nb_qa,nb_bie,nb_SG
         write(out_unitp,*) 'derive_termQact(:)',derive_termQact(:)
         write(out_unitp,*) 'derive_termQdyn(:)',derive_termQdyn(:)
         write(out_unitp,*) 'grid_cte',OpGrid%grid_cte
         write(out_unitp,*) 'Save_MemGrid',OpGrid%para_FileGrid%Save_MemGrid
         CALL flush_perso(out_unitp)
       END IF

       OpGrid%nb_qa      = nb_qa
       OpGrid%nb_bie     = nb_bie

       CALL alloc_array(OpGrid%derive_termQact,(/2/),                &
                       "OpGrid%derive_termQact",name_sub)
       OpGrid%derive_termQact(:) = derive_termQact(:)

       CALL alloc_array(OpGrid%derive_termQdyn,(/2/),                &
                       "OpGrid%derive_termQdyn",name_sub)
       OpGrid%derive_termQdyn(:) = derive_termQdyn(:)

       IF (nb_bie <1) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nb_bie <1',nb_bie
          STOP
       END IF

       info2 = name_sub // ' of ' // trim(info)
       CALL alloc_array(OpGrid%Mat_cte,(/nb_bie,nb_bie/),            &
                       "OpGrid%Mat_cte",info2)
       OpGrid%Mat_cte(:,:) = ZERO
       IF (debug) write(out_unitp,*) info2,'Mat_cte(:,:)',size(OpGrid%Mat_cte)

       IF (.NOT. OpGrid%grid_cte .AND. OpGrid%para_FileGrid%Save_MemGrid) THEN
         CALL alloc_array(OpGrid%Grid,(/nb_qa,nb_bie,nb_bie/),       &
                         "OpGrid%Grid",info2)
         OpGrid%Grid(:,:,:) = ZERO

         IF (print_level > -1) write(out_unitp,*) info2,size(OpGrid%Grid)

         IF (SmolyakRep) THEN
           IF (print_level > -1) write(out_unitp,*) info2 // ': OpGrid%SRep allocated'
           CALL alloc_SmolyakRep_only(OpGrid%SRep,nb_SG,delta=.FALSE.,grid=.TRUE.,nb0=nb_bie)
         END IF
       END IF
       CALL flush_perso(out_unitp)

       IF (debug) THEN
         write(out_unitp,*) 'END ',name_sub
         CALL flush_perso(out_unitp)
       END IF

      END SUBROUTINE alloc_OpGrid

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      SUBROUTINE dealloc_OpGrid(OpGrid,keep_FileGrid)
          TYPE (param_OpGrid), intent(inout) :: OpGrid
          logical, intent(in) :: keep_FileGrid

          character (len=line_len)   :: FileName_Grid


          character (len=*), parameter :: name_sub='dealloc_OpGrid'

          IF (associated(OpGrid%derive_termQact)) THEN
            CALL dealloc_array(OpGrid%derive_termQact,                  &
                              "OpGrid%derive_termQact",name_sub)
          END IF

          IF (associated(OpGrid%derive_termQdyn)) THEN
            CALL dealloc_array(OpGrid%derive_termQdyn,                  &
                              "OpGrid%derive_termQdyn",name_sub)
          END IF

          IF (associated(OpGrid%Grid)) THEN
            CALL dealloc_array(OpGrid%Grid,"OpGrid%Grid",name_sub)
          END IF
          IF (associated(OpGrid%Mat_cte)) THEN
            CALL dealloc_array(OpGrid%Mat_cte,"OpGrid%Mat_cte",name_sub)
          END IF

          OpGrid%nb_qa      = 0
          OpGrid%nb_bie     = 0
          OpGrid%grid_zero  = .FALSE.
          OpGrid%grid_cte   = .FALSE.
          OpGrid%ana_grid   = .FALSE.
          OpGrid%cplx       = .FALSE.

          IF (.NOT. keep_FileGrid) THEN
            CALL file_delete(OpGrid%file_Grid)
          END IF

          CALL dealloc_SmolyakRep(OpGrid%SRep) ! for SG4


      END SUBROUTINE dealloc_OpGrid

      SUBROUTINE Write_OpGrid(OpGrid)

      TYPE (param_OpGrid), intent(inout) :: OpGrid


      write(out_unitp,*) 'BEGINNING Write_OpGrid'

      write(out_unitp,*) 'nb_qa,nb_bie      ',OpGrid%nb_qa,OpGrid%nb_bie
      IF (associated(OpGrid%derive_termQact)) THEN
        write(out_unitp,*) 'derive_termQact   ',OpGrid%derive_termQact
      ELSE
        write(out_unitp,*) 'derive_termQact: not associated'

      END IF
      IF (associated(OpGrid%derive_termQdyn)) THEN
        write(out_unitp,*) 'derive_termQdyn   ',OpGrid%derive_termQdyn
      ELSE
        write(out_unitp,*) 'derive_termQdyn: not associated'
      END IF
      write(out_unitp,*) 'grid_zero         ',OpGrid%grid_zero
      write(out_unitp,*) 'grid_cte          ',OpGrid%grid_cte
      write(out_unitp,*) 'ana_grid          ',OpGrid%ana_grid
      write(out_unitp,*) 'cplx              ',OpGrid%cplx
      write(out_unitp,*) 'alloc_Grid        ',OpGrid%alloc_Grid
      write(out_unitp,*) 'Grid_done         ',OpGrid%Grid_done
      write(out_unitp,*) 'asso Grid         ',associated(OpGrid%Grid)
      write(out_unitp,*) 'asso Mat_cte      ',associated(OpGrid%Mat_cte)
      write(out_unitp,*) 'alloc Smolyak Rep ',allocated(OpGrid%SRep%SmolyakRep)
      write(out_unitp,*) 'iq_min,Op_min     ',OpGrid%iq_min,OpGrid%Op_min
      write(out_unitp,*) 'iq_max,Op_max     ',OpGrid%iq_max,OpGrid%Op_max
      CALL flush_perso(out_unitp)

      IF (associated(OpGrid%Mat_cte)) &
         write(out_unitp,*) 'Mat_cte           ',OpGrid%Mat_cte
      CALL flush_perso(out_unitp)

      write(out_unitp,*) 'END Write_OpGrid'
      CALL flush_perso(out_unitp)

      END SUBROUTINE Write_OpGrid

      SUBROUTINE alloc_array_OF_OpGriddim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      TYPE (param_OpGrid), pointer, intent(inout) :: tab(:)
      integer,                      intent(in)    :: tab_ub(:)
      integer, optional,            intent(in)    :: tab_lb(:)
      character (len=*),            intent(in)    :: name_var,name_sub

      integer, parameter :: ndim=1
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_OpGriddim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (associated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)

         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1)),stat=err_mem)
       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'param_OpGrid')

      END SUBROUTINE alloc_array_OF_OpGriddim1
      SUBROUTINE dealloc_array_OF_OpGriddim1(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (param_OpGrid), pointer, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_OpGriddim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN
       IF (.NOT. associated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'param_OpGrid')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_OpGriddim1

      SUBROUTINE Set_file_OF_OpGrid(OpGrid,Type_FileGrid,iOp,           &
                                      Base_FileName_Grid,name_Op,nb_bie)

      TYPE (param_OpGrid), pointer, intent(inout)   :: OpGrid(:)
      integer, intent(in)                           :: iOp,nb_bie,Type_FileGrid
      character (len=*), intent(in)                 :: name_Op
      character (len=Line_len), intent(in)          :: Base_FileName_Grid


      integer :: iterm,frecl
      real (kind=Rkind) :: Mat(nb_bie,nb_bie) ! for the direct acces file

      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Set_file_OF_OpGrid'

      !write(6,*) ' in ',name_sub,' asso OpGrid ',associated(OpGrid), 'name_Op: ',name_Op
      IF (.NOT. associated(OpGrid)) RETURN
      IF (size(OpGrid) < 1) RETURN

      DO iterm=1,size(OpGrid)
        IF (Grid_omp == 0) THEN
          OpGrid(iterm)%file_Grid%nb_thread = 1
        ELSE
          OpGrid(iterm)%file_Grid%nb_thread = Grid_maxth
        END IF
        OpGrid(iterm)%para_fileGrid%Type_FileGrid = Type_FileGrid

        SELECT CASE (OpGrid(iterm)%para_fileGrid%Type_FileGrid)
        CASE (0) ! normal SH_HADA file
          OpGrid(iterm)%file_Grid%seq = .TRUE.
          ! formatted is already defined
          STOP

        CASE (1) ! sequential
          OpGrid(iterm)%file_Grid%seq       = .TRUE.
          OpGrid(iterm)%file_Grid%formatted = .FALSE.
          OpGrid(iterm)%file_Grid%init      = .TRUE.

          IF (OpGrid(iterm)%cplx) THEN
            OpGrid(iterm)%file_Grid%name = trim(adjustl(Base_FileName_Grid)) // &
                                               "_im" // trim(adjustl(name_Op))
          ELSE
            OpGrid(iterm)%file_Grid%name = trim(adjustl(Base_FileName_Grid)) // &
                     "_" // trim(adjustl(name_Op)) // int_TO_char(iterm)
          END IF

        CASE (2) ! direct
          OpGrid(iterm)%file_Grid%seq       = .FALSE.
          OpGrid(iterm)%file_Grid%formatted = .FALSE.
          OpGrid(iterm)%file_Grid%init      = .TRUE.
          OpGrid(iterm)%file_Grid%nb_thread = 0 ! its means only one file, but it can be used with sevral threads

          INQUIRE(iolength=frecl) Mat(:,:)
          OpGrid(iterm)%file_Grid%frecl = frecl

          IF (OpGrid(iterm)%cplx) THEN
            OpGrid(iterm)%file_Grid%name = trim(adjustl(Base_FileName_Grid)) // &
                                               "_im" // trim(adjustl(name_Op))
          ELSE
            OpGrid(iterm)%file_Grid%name = trim(adjustl(Base_FileName_Grid)) // &
                         "_" // trim(adjustl(name_Op)) // int_TO_char(iterm)
          END IF

        CASE (4) ! sequential for SG4 (here it just the base name for each Smolyak term)
          OpGrid(iterm)%file_Grid%seq       = .TRUE.
          OpGrid(iterm)%file_Grid%formatted = .FALSE.
          OpGrid(iterm)%file_Grid%init      = .TRUE.

          IF (OpGrid(iterm)%cplx) THEN
            OpGrid(iterm)%file_Grid%name = trim(adjustl(Base_FileName_Grid)) // &
                                               "_im" // trim(adjustl(name_Op))
          ELSE
            OpGrid(iterm)%file_Grid%name = trim(adjustl(Base_FileName_Grid)) // &
                     "_" // trim(adjustl(name_Op)) // int_TO_char(iterm)
          END IF

        CASE default ! normal SH_HADA file
          STOP
        END SELECT

      END DO

      END SUBROUTINE Set_file_OF_OpGrid

      SUBROUTINE Open_file_OF_OpGrid(OpGrid,nio)

      TYPE (param_OpGrid), pointer, intent(inout) :: OpGrid(:)
      integer, intent(inout)                      :: nio

      integer :: iterm
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Open_file_OF_OpGrid'

      IF (.NOT. associated(OpGrid)) RETURN
      IF (size(OpGrid) < 1) RETURN

      DO iterm=1,size(OpGrid)

        SELECT CASE (OpGrid(iterm)%para_FileGrid%Type_FileGrid)
        CASE (0) ! normal SH_HADA file
          CONTINUE ! nothing the file is opened/closed at each grid points
        CASE (1,2) ! sequential

          CALL file_open(OpGrid(iterm)%file_Grid,nio)

        CASE (4) ! for SG4: one file per Smolyak term
          CONTINUE ! nothing the files will be open when need
        CASE default ! normal SH_HADA file
          CONTINUE ! nothing the file is opened/closed at each grid points
        END SELECT


      END DO

      END SUBROUTINE Open_file_OF_OpGrid
      SUBROUTINE Close_file_OF_OpGrid(OpGrid)

      TYPE (param_OpGrid), pointer, intent(inout) :: OpGrid(:)

      integer :: iterm
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Close_file_OF_OpGrid'

      IF (.NOT. associated(OpGrid)) RETURN
      IF (size(OpGrid) < 1) RETURN

      DO iterm=1,size(OpGrid)

        SELECT CASE (OpGrid(iterm)%para_FileGrid%Type_FileGrid)
        CASE (0) ! normal SH_HADA file
          CONTINUE ! nothing the file is opened/closed at each grid points
        CASE (1,2) ! sequential

          CALL file_close(OpGrid(iterm)%file_Grid)

        CASE (4) ! for SG4: one file per Smolyak term
          CONTINUE ! nothing the files will be opened/closed when need
        CASE default ! normal SH_HADA file
          CONTINUE ! nothing the file is opened/closed at each grid points
        END SELECT


      END DO

      END SUBROUTINE Close_file_OF_OpGrid


      SUBROUTINE sub_ReadDir_Grid_iterm(Grid,OpGrid)
      USE mod_system
      IMPLICIT NONE

      real (kind=Rkind), intent(inout) :: Grid(:,:,:) ! grid when Save_Grid_iterm=t
      TYPE (param_OpGrid) :: OpGrid


      integer :: i_qa,lrecl_Grid_iterm,nio,error,nb_qa

      !-- Read the direct acces file, then delete it
      lrecl_Grid_iterm = OpGrid%file_Grid%frecl
      nb_qa = size(Grid,dim=1)

      CALL file_open(OpGrid%file_Grid,                            &
                      nio,lformatted=.FALSE.,seq=.FALSE.,               &
                      lrecl=lrecl_Grid_iterm)
       DO i_qa=1,nb_qa
         read(nio,REC=i_qa,iostat=error) Grid(i_qa,:,:)

         IF (error /= 0) THEN
           write(out_unitp,*) ' ERROR in sub_ReadDir_Grid_iterm'
           write(out_unitp,*) ' Impossible to read the file: ',         &
                  OpGrid%file_Grid%name
           write(out_unitp,*) ' i_qa,nb_qa',i_qa,nb_qa

           STOP
         END IF
       END DO
       CALL file_close(OpGrid%file_Grid)

      END SUBROUTINE sub_ReadDir_Grid_iterm
      SUBROUTINE sub_ReadSeq_Grid_iterm(Grid,OpGrid)
      USE mod_system
      IMPLICIT NONE

      real (kind=Rkind), intent(inout) :: Grid(:,:,:) ! grid when Save_Grid_iterm=t
      TYPE (param_OpGrid) :: OpGrid


      integer :: i_qa,nio,error,ithread,nb_qa
      logical :: file_is_para

      CALL file_open(OpGrid%file_Grid,nio,lformatted=.FALSE.)
      nb_qa = size(Grid,dim=1)
      ! for parallel calculation of the Grid file
      IF (OpGrid%file_Grid%nb_thread > 1) THEN
        !write(out_unitp,*) ' OMP: read one term of HADA file'
        ithread      = 0
        file_is_para = .TRUE.
        nio = OpGrid%file_Grid%tab_unit(ithread)
      ELSE
        nio = OpGrid%file_Grid%unit
        !write(out_unitp,*) ' non OMP: read one term of HADA file'
      END IF

      DO i_qa=1,nb_qa
         read(nio,iostat=error) Grid(i_qa,:,:)

         IF (error > 0 .OR. error < 0 .AND. .NOT. file_is_para) THEN
           write(out_unitp,*) ' ERROR in sub_ReadSeq_Grid_iterm'
           write(out_unitp,*) ' ERROR on file: ',OpGrid%file_Grid%name
           write(out_unitp,*) ' i_qa,nb_qa',i_qa,nb_qa
           STOP
         ELSE IF (error < 0 .AND. file_is_para) THEN ! end of file
           ithread = ithread + 1
           nio = OpGrid%file_Grid%tab_unit(ithread)

           read(nio) Grid(i_qa,:,:)

         END IF
         !write(out_unitp,*) 'i_qa,Grid',i_qa,Grid(i_qa,:,:)

       END DO

       CALL file_close(OpGrid%file_Grid)

      END SUBROUTINE sub_ReadSeq_Grid_iterm

!================================================================
!     Analysis of the grid (zero or constant terms)
!================================================================
      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      SUBROUTINE Analysis_OpGrid(OpGrid,n_Op)

      TYPE (param_OpGrid), pointer, intent(inout) :: OpGrid(:)
      integer,                      intent(in)    :: n_Op

      integer           :: k,k_term,iq
      real (kind=Rkind) :: Op_temp


      character (len=*), parameter :: name_sub='Analysis_OpGrid'

      IF (.NOT. associated(OpGrid) ) RETURN
      IF (size(OpGrid) < 1 ) RETURN

      IF (print_level>-1) THEN
        write(out_unitp,*)'--------------------------------------------------------------'
        write(out_unitp,*)'n_Op,k_term,derive_term,cte,zero,    minval,    maxval,dealloc'
        CALL flush_perso(out_unitp)
      END IF

      DO k_term=1,size(OpGrid)

        IF (associated(OpGrid(k_term)%Grid) .AND. OpGrid(k_term)%Grid_done) THEN
          IF (.NOT. OpGrid(k_term)%grid_cte) THEN
            OpGrid(k_term)%Mat_cte(:,:) = OpGrid(k_term)%Grid(1,:,:)
            OpGrid(k_term)%grid_cte = .TRUE.
            DO iq=1,OpGrid(k_term)%nb_qa
              OpGrid(k_term)%grid_cte =                                 &
                 OpGrid(k_term)%grid_cte .AND.                          &
                 (sum(abs(OpGrid(k_term)%Grid(iq,:,:) -                 &
                    OpGrid(k_term)%Grid(1,:,:))) < ONETENTH**12)
              IF (.NOT. OpGrid(k_term)%grid_cte) EXIT
            END DO

            OpGrid(k_term)%Op_min = huge(ONE)
            OpGrid(k_term)%Op_max = -huge(ONE)
            OpGrid(k_term)%iq_min = 0
            OpGrid(k_term)%iq_max = 0
            DO iq=1,OpGrid(k_term)%nb_qa

              DO k=1,OpGrid(k_term)%nb_bie
                Op_temp = OpGrid(k_term)%Grid(iq,k,k)

                IF (Op_temp < OpGrid(k_term)%Op_min) THEN
                  OpGrid(k_term)%Op_min = Op_temp
                  OpGrid(k_term)%iq_min = iq
                END IF
                IF (Op_temp > OpGrid(k_term)%Op_max) THEN
                  OpGrid(k_term)%Op_max = Op_temp
                  OpGrid(k_term)%iq_max = iq
                END IF
              END DO
            END DO

          END IF

          OpGrid(k_term)%grid_zero = (sum(abs(                          &
            OpGrid(k_term)%Mat_cte(:,:))) < ONETENTH**12) .AND.         &
            OpGrid(k_term)%grid_cte

          IF (OpGrid(k_term)%grid_cte .AND. associated(OpGrid(k_term)%Grid)) THEN
            CALL dealloc_array(OpGrid(k_term)%Grid,"OpGrid%Grid",name_sub)
          END IF

          IF (print_level>-1) THEN
            write(out_unitp,'(i5,x,i6,2x,2i5,l3,x,l4,x,2e11.2,x,l3)')   &
                   n_Op,k_term,OpGrid(k_term)%derive_termQact(:),       &
                   OpGrid(k_term)%grid_cte,OpGrid(k_term)%grid_zero,    &
                   OpGrid(k_term)%Op_min,OpGrid(k_term)%Op_max,         &
                   (.NOT. associated(OpGrid(k_term)%Grid))
          END IF

        ELSE

          OpGrid(k_term)%grid_zero = OpGrid(k_term)%grid_cte .AND.      &
                 (sum(abs(OpGrid(k_term)%Mat_cte(:,:))) < ONETENTH**12)

          IF (print_level>-1) THEN
            write(out_unitp,'(i5,x,i6,2x,2i5,l3,x,l4,24x,l3)') n_Op,k_term,    &
                   OpGrid(k_term)%derive_termQact(:),                   &
                   OpGrid(k_term)%grid_cte,OpGrid(k_term)%grid_zero,    &
                   (.NOT. associated(OpGrid(k_term)%Grid))

          END IF
        END IF

      END DO
      IF (print_level>-1) THEN
        write(out_unitp,*)'--------------------------------------------------------------'
        CALL flush_perso(out_unitp)
      END IF

      END SUBROUTINE Analysis_OpGrid


      END MODULE mod_OpGrid

