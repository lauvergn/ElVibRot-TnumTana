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
!
!=====================================================================
!  save on ONE file the matrices :
!       S_bhe, H_bhe, Veff_bhe, im_Veff_bhe
!       T1_bhe, T2_bhe ....
!=====================================================================
      SUBROUTINE sub_Save_GridFile_AllOp(iq,d0MatOp,nb_Op,           &
                                         Qdyn,nb_var,Qact,nb_act1,   &
                                         para_AllOp,w)


      USE mod_system
      USE mod_PrimOp
      USE mod_Op
      IMPLICIT NONE

      TYPE (param_AllOp) :: para_AllOp
      integer       :: nb_var,nb_act1
      real (kind=Rkind), intent(inout)        :: Qdyn(nb_var),Qact(nb_act1)


      integer            :: iq

!----- variables for calc_Op (mat_V, mat_imV and mat_ScalOp) ----------------
      integer                  :: nb_Op
      TYPE (param_d0MatOp)     :: d0MatOp(nb_Op)


      real (kind=Rkind) :: w



      character (len=name_len) :: name_Op

      integer       :: nio
      integer       :: iOp,n_Op
      integer       :: i1,i2,nb_bie
      integer       :: iterm,nb_allterms,k_term
      logical       :: Save_FileGrid_Op

      integer       :: nb_thread,ithread

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_Save_GridFile_AllOp'
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'iq',iq
         write(out_unitp,*) 'w',w
         write(out_unitp,*) 'nb_Op',para_AllOp%nb_Op
         DO iOp=1,para_AllOp%nb_Op
           write(out_unitp,*) 'iOp,n_Op,nb_term',                              &
                iOp,para_AllOp%tab_Op(iOp)%n_Op,para_AllOp%tab_Op(iOp)%nb_term
         END DO
       END IF
!-----------------------------------------------------------
       nb_bie          = para_AllOp%tab_Op(1)%nb_bie

      !---------------------------------------------------------
      !allocation of OpGrid
      !---------------------------------------------------------
      DO iOp=1,para_AllOp%nb_Op
        IF (para_AllOp%tab_Op(iOp)%n_Op == -1) CYCLE  ! for S

        !--- alloc OpGrid ------------------------------------
        IF (.NOT. associated(para_AllOp%tab_Op(iOp)%OpGrid))            &
           CALL alloc_para_Op(para_AllOp%tab_Op(iOp),Mat=.FALSE.,Grid=.TRUE.)
      END DO

      !---------------------------------------------------------
      !whitout saving on SH_HADA file
      !---------------------------------------------------------
      IF (para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Type_FileGrid /= 0) THEN !test on H=tab_Op(1)

        DO iOp=1,para_AllOp%nb_Op

          IF (para_AllOp%tab_Op(iOp)%n_Op == -1) THEN ! S=tab_Op(2), useless
            para_AllOp%tab_Op(iOp)%para_ReadOp%para_FileGrid%Save_FileGrid      = .FALSE.
            para_AllOp%tab_Op(iOp)%para_ReadOp%para_FileGrid%Save_FileGrid_done = .FALSE.
            para_AllOp%tab_Op(iOp)%OpGrid(1)%para_FileGrid%Save_FileGrid        = .FALSE.
            para_AllOp%tab_Op(iOp)%OpGrid(1)%para_FileGrid%Save_FileGrid_done   = .FALSE.

            CYCLE
          END IF

          Save_FileGrid_Op = .TRUE.
          DO k_term=1,para_AllOp%tab_Op(iOp)%nb_term

            Save_FileGrid_Op = Save_FileGrid_Op .AND.                   &
               para_AllOp%tab_Op(iOp)%OpGrid(k_term)%para_FileGrid%Save_FileGrid

            IF (.NOT. para_AllOp%tab_Op(iOp)%OpGrid(k_term)%para_FileGrid%Save_FileGrid) CYCLE
            IF (para_AllOp%tab_Op(iOp)%OpGrid(k_term)%Grid_cte) CYCLE

            para_AllOp%tab_Op(iOp)%OpGrid(k_term)%para_FileGrid%Save_FileGrid_done = .TRUE.

            i1 = para_AllOp%tab_Op(iOp)%derive_termQact(1,k_term)
            i2 = para_AllOp%tab_Op(iOp)%derive_termQact(2,k_term)
            IF ((i1 < 0 .OR. i2 < 0) .AND. para_AllOp%tab_Op(1)%para_Tnum%JJ == 0) CYCLE
            iterm = d0MatOp(iOp)%derive_term_TO_iterm(i1,i2)

!$OMP       CRITICAL (sub_Save_GridFile_AllOp_CRIT1)

            nio = file_GetUnit(para_AllOp%tab_Op(iOp)%OpGrid(k_term)%file_Grid)

            IF (para_AllOp%tab_Op(iOp)%OpGrid(k_term)%para_FileGrid%Type_FileGrid == 1) THEN
              ! sequential access
              write(nio) d0MatOp(iOp)%ReVal(:,:,iterm)
              flush(nio)
            ELSE ! direct access (Type_FileGrid=2)
              write(nio,rec=iq) d0MatOp(iOp)%ReVal(:,:,iterm)
            END IF

!$OMP       END CRITICAL (sub_Save_GridFile_AllOp_CRIT1)

          END DO

          IF (para_AllOp%tab_Op(iOp)%cplx) THEN

            Save_FileGrid_Op = Save_FileGrid_Op .AND.                   &
                para_AllOp%tab_Op(iOp)%ImOpGrid(1)%para_FileGrid%Save_FileGrid

            IF (.NOT. para_AllOp%tab_Op(iOp)%ImOpGrid(1)%para_FileGrid%Save_FileGrid) CYCLE
            IF (para_AllOp%tab_Op(iOp)%ImOpGrid(1)%Grid_cte) CYCLE

            para_AllOp%tab_Op(iOp)%ImOpGrid(1)%para_FileGrid%Save_FileGrid_done = .TRUE.
!$OMP       CRITICAL (sub_Save_GridFile_AllOp_CRIT2)
            ! sequential access
            IF (para_AllOp%tab_Op(iOp)%ImOpGrid(1)%para_FileGrid%Type_FileGrid == 1) THEN
              ithread      = 0
!$            ithread      = OMP_GET_THREAD_NUM()

              IF (Grid_omp == 0) THEN
                nb_thread = 1
              ELSE
                nb_thread = Grid_maxth
              END IF
              !write(out_unitp,*) 'nb_thread in sub_saving5_Op_all: ',nb_thread

              IF (nb_thread > 1) THEN
                nio = para_AllOp%tab_Op(iOp)%imOpGrid(1)%file_Grid%tab_unit(ithread)
              ELSE
                nio = para_AllOp%tab_Op(iOp)%imOpGrid(1)%file_Grid%unit
              END IF
              write(nio) d0MatOp(iOp)%ImVal(:,:)
              flush(nio)
            ELSE! direct access (Type_FileGrid=2)
              nio = para_AllOp%tab_Op(iOp)%imOpGrid(1)%file_Grid%unit
              write(nio,rec=iq) d0MatOp(iOp)%ImVal(:,:)
            END IF
!$OMP       END CRITICAL (sub_Save_GridFile_AllOp_CRIT2)

          END IF
          para_AllOp%tab_Op(iOp)%para_ReadOp%para_FileGrid%Save_FileGrid_done = Save_FileGrid_Op

        END DO


      ELSE IF (para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Save_FileGrid) THEN ! para_AllOp%tab_Op(1)%para_ReadOp%para_FileGrid%Type_FileGrid=0

        !---------------------------------------------------------
        !Saving on SH_HADA file
        !---------------------------------------------------------
        !---------------------------------------------------------
        nb_allterms = 0
        DO iOp=1,para_AllOp%nb_Op
          nb_allterms = nb_allterms + para_AllOp%tab_Op(iOp)%nb_term
          IF (para_AllOp%tab_Op(iOp)%cplx) nb_allterms = nb_allterms + 1
        END DO

!$OMP   CRITICAL (sub_Save_GridFile_AllOp_CRIT3)
        ithread      = 0
!$      ithread      = OMP_GET_THREAD_NUM()
        IF (para_AllOp%tab_Op(1)%ComOp%file_HADA%nb_thread > 1 .AND. iq == 1) THEN
          CALL file_open(para_AllOp%tab_Op(1)%ComOp%file_HADA,nio,      &
             lformatted=para_AllOp%tab_Op(1)%ComOp%file_HADA%formatted)
         nio = para_AllOp%tab_Op(1)%ComOp%file_HADA%unit
         IF (para_AllOp%tab_Op(1)%ComOp%file_HADA%formatted) THEN

           write(nio,*) '- Beginning_th ',                              &
                        para_AllOp%tab_Op(1)%ComOp%file_HADA%nb_thread, &
                       '-------------------'
         ELSE
            !write(nio) 'Beginning_th'
            write(nio) para_AllOp%tab_Op(1)%ComOp%file_HADA%nb_thread
         END IF
         CALL file_close(para_AllOp%tab_Op(1)%ComOp%file_HADA)
       END IF

       IF (para_AllOp%tab_Op(1)%ComOp%file_HADA%nb_thread > 1) THEN
         CALL file_open2(                                               &
          para_AllOp%tab_Op(1)%ComOp%file_HADA%tab_name_th(ithread),nio,&
              lformatted=para_AllOp%tab_Op(1)%ComOp%file_HADA%formatted,&
                     append=.TRUE.)
       ELSE
         CALL file_open2(                                               &
            para_AllOp%tab_Op(1)%ComOp%file_HADA%name,nio,              &
              lformatted=para_AllOp%tab_Op(1)%ComOp%file_HADA%formatted,&
                     append=.TRUE.)
       END IF
       !---------------------------------------------------------
       nb_var  = para_AllOp%tab_Op(1)%mole%nb_var
       nb_act1 = para_AllOp%tab_Op(1)%mole%nb_act1
       nb_bie  = para_AllOp%tab_Op(1)%nb_bie

       IF (para_AllOp%tab_Op(1)%ComOp%file_HADA%formatted) THEN

          !- write parameters at iq --------------------------------
          write(nio,*) '- Beginning ---------------------------------'
          write(nio,*) 'iq: ',iq
          write(nio,*) 'nb_bie:',nb_bie
          write(nio,*) 'w: ',w
          write(nio,*) 'Qdyn: ',nb_var
          CALL Write_Vec(Qdyn,nio,5,Rformat='e30.23')
          write(nio,*) 'Qact: ',nb_act1
          CALL Write_Vec(Qact(1:nb_act1),nio,5,Rformat='e30.23')
          write(nio,*) 'JJ: ',para_AllOp%tab_Op(1)%para_Tnum%JJ
          write(nio,*) 'pot_cplx: ',para_AllOp%tab_Op(1)%cplx
          write(nio,*) 'calc_scalar_Op: ',para_AllOp%tab_Op(1)%para_PES%calc_scalar_Op
          write(nio,*) 'nb_allterms: ',nb_allterms

          DO iOp=1,para_AllOp%nb_Op
            para_AllOp%tab_Op(iOp)%para_ReadOp%para_FileGrid%Save_FileGrid_done = .TRUE.

            n_Op = para_AllOp%tab_Op(iOp)%n_Op

            name_Op = '-' // trim(para_AllOp%tab_Op(iOp)%name_Op) // '-----'

            !write(6,*) 'para_AllOp%tab_Op(iOp)%nb_term',iOp,para_AllOp%tab_Op(iOp)%nb_term
            DO k_term=1,para_AllOp%tab_Op(iOp)%nb_term
              para_AllOp%tab_Op(iOp)%OpGrid(k_term)%para_FileGrid%Save_FileGrid_done = .TRUE.

              i1 = para_AllOp%tab_Op(iOp)%derive_termQact(1,k_term)
              i2 = para_AllOp%tab_Op(iOp)%derive_termQact(2,k_term)
              IF ((i1 < 0 .OR. i2 < 0) .AND. para_AllOp%tab_Op(1)%para_Tnum%JJ == 0) CYCLE
              iterm = d0MatOp(iOp)%derive_term_TO_iterm(i1,i2)

              write(nio,*) name_Op,n_Op,i1,i2,.FALSE.
              CALL Write_Mat(d0MatOp(iOp)%ReVal(:,:,iterm),nio,5,Rformat='e30.23')
            END DO

            IF (para_AllOp%tab_Op(iOp)%cplx) THEN
              para_AllOp%tab_Op(iOp)%ImOpGrid(1)%para_FileGrid%Save_FileGrid_done = .TRUE.

              name_Op = '-im_' // trim(para_AllOp%tab_Op(iOp)%name_Op)  &
                         // '-----'
              write(nio,*) name_Op,n_Op,0,0,.TRUE.
              CALL Write_Mat(d0MatOp(iOp)%ImVal(:,:),nio,5,Rformat='e30.23')
            END IF
          END DO

          write(nio,*) '- End ',ithread,                                &
                      para_AllOp%tab_Op(1)%ComOp%file_HADA%nb_thread,   &
                      '----------------------'

       ELSE


          !- write parameters at iq --------------------------------
          write(nio) para_AllOp%tab_Op(1)%ComOp%file_HADA%nb_thread
          write(nio) iq
          write(nio) nb_bie
          write(nio) w
          write(nio) nb_var
          write(nio) Qdyn(:)
          write(nio) nb_act1
          write(nio) Qact(1:nb_act1)
          write(nio) para_AllOp%tab_Op(1)%para_Tnum%JJ
          write(nio) para_AllOp%tab_Op(1)%cplx
          write(nio) para_AllOp%tab_Op(1)%para_PES%calc_scalar_Op
          write(nio) nb_allterms

          DO iOp=1,para_AllOp%nb_Op
            para_AllOp%tab_Op(iOp)%para_ReadOp%para_FileGrid%Save_FileGrid_done = .TRUE.

            n_Op = para_AllOp%tab_Op(iOp)%n_Op

            DO k_term=1,para_AllOp%tab_Op(iOp)%nb_term
              para_AllOp%tab_Op(iOp)%OpGrid(k_term)%para_FileGrid%Save_FileGrid_done = .TRUE.

              i1 = para_AllOp%tab_Op(iOp)%derive_termQact(1,k_term)
              i2 = para_AllOp%tab_Op(iOp)%derive_termQact(2,k_term)
              IF ((i1 < 0 .OR. i2 < 0) .AND. para_AllOp%tab_Op(1)%para_Tnum%JJ == 0) CYCLE
              iterm = d0MatOp(iOp)%derive_term_TO_iterm(i1,i2)

              write(nio) n_Op,i1,i2,.FALSE.
              write(nio) d0MatOp(iOp)%ReVal(:,:,iterm)
            END DO

            IF (para_AllOp%tab_Op(iOp)%cplx) THEN
              para_AllOp%tab_Op(iOp)%ImOpGrid(1)%para_FileGrid%Save_FileGrid_done = .TRUE.
              write(nio) n_Op,0,0,.TRUE.
              write(nio) d0MatOp(iOp)%ImVal(:,:)
            END IF
          END DO

          write(nio) para_AllOp%tab_Op(1)%ComOp%file_HADA%nb_thread

       END IF

       !---------------------------------------------------------
       close(nio)
       !---------------------------------------------------------
!$OMP  END CRITICAL (sub_Save_GridFile_AllOp_CRIT3)
     END IF


       IF (debug) THEN
         write(out_unitp,*) 'END ',name_sub
       END IF


      END SUBROUTINE sub_Save_GridFile_AllOp
      SUBROUTINE sub_Save_GridMem_AllOp(iq,d0MatOp,nb_Op,para_AllOp)

      USE mod_system
      USE mod_PrimOp
      USE mod_Op
      IMPLICIT NONE

      TYPE (param_AllOp) :: para_AllOp

!----- variables for calc_Op (mat_V, mat_imV and mat_ScalOp) ----------------
      integer                  :: nb_Op
      TYPE (param_d0MatOp)     :: d0MatOp(nb_Op)

      integer            :: iq


      integer       :: iOp,i1,i2
      integer       :: k_term,iterm
      logical       :: SaveGrid_Op


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_Save_GridMem_AllOp'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'iq',iq
         write(out_unitp,*) 'nb_Op',para_AllOp%nb_Op
         DO iOp=1,para_AllOp%nb_Op
           write(out_unitp,*) 'iOp,n_Op,nb_term',                       &
                iOp,para_AllOp%tab_Op(iOp)%n_Op,para_AllOp%tab_Op(iOp)%nb_term
           write(out_unitp,*) 'asso OpGrid',iOp, &
                               associated(para_AllOp%tab_Op(iOp)%OpGrid)
           write(out_unitp,*) 'para_ReadOp%para_FileGrid%Save_MemGrid', &
                 para_AllOp%tab_Op(iOp)%para_ReadOp%para_FileGrid%Save_MemGrid
         END DO
       END IF
!-----------------------------------------------------------

      !---------------------------------------------------------
      !whitout saving on memory
      !---------------------------------------------------------
      DO iOp=1,para_AllOp%nb_Op
        IF (para_AllOp%tab_Op(iOp)%n_Op == -1) CYCLE  ! for S


        !--- alloc OpGrid ------------------------------------
        IF (.NOT. associated(para_AllOp%tab_Op(iOp)%OpGrid)) THEN
           CALL alloc_para_Op(para_AllOp%tab_Op(iOp),Mat=.FALSE.,Grid=.TRUE.)

!           DO k_term=1,para_AllOp%tab_Op(iOp)%nb_term
!             write(out_unitp,*) 'para_AllOp%tab_Op(iOp)%OpGrid(k_term)%para_FileGrid',iOp,k_term
!             CALL Write_FileGrid(para_AllOp%tab_Op(iOp)%OpGrid(k_term)%para_FileGrid)
!             write(out_unitp,*) 'asso para_AllOp%tab_Op(iOp)%OpGrid(k_term)%Grid',iOp,k_term, &
!                      associated(para_AllOp%tab_Op(iOp)%OpGrid(k_term)%Grid)
!           END DO
         END IF
      END DO



        DO iOp=1,para_AllOp%nb_Op

          IF (para_AllOp%tab_Op(iOp)%n_Op == -1) THEN !S=tab_Op(2), NO STORAGE: useless
            para_AllOp%tab_Op(iOp)%para_ReadOp%para_FileGrid%Save_MemGrid      = .FALSE.
            para_AllOp%tab_Op(iOp)%para_ReadOp%para_FileGrid%Save_MemGrid_done = .FALSE.
            para_AllOp%tab_Op(iOp)%OpGrid(1)%para_FileGrid%Save_MemGrid        = .FALSE.
            para_AllOp%tab_Op(iOp)%OpGrid(1)%para_FileGrid%Save_MemGrid_done   = .FALSE.
            CYCLE
          END IF

          SaveGrid_Op = .TRUE.
          DO k_term=1,para_AllOp%tab_Op(iOp)%nb_term

            SaveGrid_Op = SaveGrid_Op .AND.                             &
               para_AllOp%tab_Op(iOp)%OpGrid(k_term)%para_FileGrid%Save_MemGrid

            IF (para_AllOp%tab_Op(iOp)%OpGrid(k_term)%Grid_cte) CYCLE

            IF (para_AllOp%tab_Op(iOp)%OpGrid(k_term)%para_FileGrid%Save_MemGrid) THEN

              para_AllOp%tab_Op(iOp)%OpGrid(k_term)%para_FileGrid%Save_MemGrid_done = .TRUE.
              i1 = para_AllOp%tab_Op(iOp)%derive_termQact(1,k_term)
              i2 = para_AllOp%tab_Op(iOp)%derive_termQact(2,k_term)
              IF ((i1 < 0 .OR. i2 < 0) .AND. para_AllOp%tab_Op(1)%para_Tnum%JJ == 0) CYCLE
              iterm = d0MatOp(iOp)%derive_term_TO_iterm(i1,i2)

              para_AllOp%tab_Op(iOp)%OpGrid(k_term)%Grid(iq,:,:) =      &
                                           d0MatOp(iOp)%ReVal(:,:,iterm)

            END IF

          END DO

          IF (para_AllOp%tab_Op(iOp)%cplx) THEN

            SaveGrid_Op = SaveGrid_Op .AND.                           &
                 para_AllOp%tab_Op(iOp)%ImOpGrid(1)%para_FileGrid%Save_MemGrid

            IF (para_AllOp%tab_Op(iOp)%ImOpGrid(1)%Grid_cte) CYCLE


            IF (para_AllOp%tab_Op(iOp)%imOpGrid(1)%para_FileGrid%Save_MemGrid) THEN

              para_AllOp%tab_Op(iOp)%imOpGrid(1)%para_FileGrid%Save_MemGrid_done = .TRUE.
              para_AllOp%tab_Op(iOp)%imOpGrid(1)%Grid(iq,:,:) =         &
                                                 d0MatOp(iOp)%ImVal(:,:)

            END IF

          END IF
          para_AllOp%tab_Op(iOp)%para_ReadOp%para_FileGrid%Save_MemGrid_done = SaveGrid_Op
        END DO

       IF (debug) THEN
         write(out_unitp,*) 'END ',name_sub
       END IF


      END SUBROUTINE sub_Save_GridMem_AllOp

!=====================================================================
!
!  read on files the matrices nb_term * Op_bhe and imOp_bhe
!
!=====================================================================

      SUBROUTINE sub_reading_Op(iq,nb_qa,d0MatOp,n_Op,Qdyn,nb_var,Qact,w,ComOp)


      USE mod_system
      USE mod_Op
      IMPLICIT NONE

      TYPE (param_d0MatOp)     :: d0MatOp

      integer, intent(in)                     :: iq,nb_qa,nb_var,n_Op
      TYPE (param_ComOp), intent(inout)       :: ComOp

      real (kind=Rkind), intent(inout)        :: Qdyn(nb_var),Qact(ComOp%nb_act1)
      real (kind=Rkind), intent(inout)        :: w



      integer            :: n_Op_lect,nb_Op
      real (kind=Rkind), allocatable  :: work_bhe(:,:)

      integer            :: id1,id2

      character (len=Name_len) :: name1,name2

      integer       :: i,i_term
      integer       :: n1,n2,n3,iqr

      logical       :: pot_cplx,calc_scalar_Op,cplx
      integer       :: JJ

      logical, save                  :: file_is_para = .FALSE.
      integer, save                  :: nb_thread_file = 0
      integer                        :: iocond,err
      integer, save                  :: nio,ithread

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_reading_Op'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'iq',iq
         write(out_unitp,*) 'nb_bi',ComOp%nb_bi
         write(out_unitp,*) 'nb_var,nb_act1,nb_bie',nb_var,ComOp%nb_act1,d0MatOp%nb_bie
         write(out_unitp,*) 'formatted_HADA',ComOp%file_HADA%formatted
         write(out_unitp,*)
       END IF
!-----------------------------------------------------------

      CALL alloc_NParray(work_bhe,(/d0MatOp%nb_bie,d0MatOp%nb_bie/),    &
                        'work_bhe',name_sub)

      work_bhe(:,:) = ZERO

      IF (iq == 1) THEN
        !write(out_unitp,*) 'read SH_HADA file: ',ComOp%file_HADA%name
        ComOp%file_HADA%nb_thread = 0
        CALL file_open(ComOp%file_HADA,nio,                             &
                       lformatted=ComOp%file_HADA%formatted)

        IF (ComOp%file_HADA%formatted) THEN
          read(nio,*) name1,name2

          IF (debug) write(out_unitp,*) 'init,name2: ',name2

          CALL file_close(ComOp%file_HADA)

          ! for parallel calculation of SH_HADA file (the file is split in SH_HADA.0, SH_HADA.1 ...)
          IF (trim(name2) == 'Beginning_th') THEN
            !write(out_unitp,*) ' OMP calc of HADA file'
            ithread      = 0
            file_is_para = .TRUE.
            CALL file_open(ComOp%file_HADA,nio,                         &
                                   lformatted=ComOp%file_HADA%formatted)
            read(nio,*) name1,name2,nb_thread_file
            CALL file_close(ComOp%file_HADA)

            IF (debug) write(out_unitp,*) 'file_is_para,nb_thread_file',file_is_para,nb_thread_file

            ComOp%file_HADA%nb_thread = nb_thread_file
          ELSE
            file_is_para = .FALSE.
          END IF

          CALL file_open(ComOp%file_HADA,nio,                           &
                                   lformatted=ComOp%file_HADA%formatted)

          IF (file_is_para) THEN
            ithread = 0
            nio = ComOp%file_HADA%tab_unit(ithread)
          ELSE
            nio = ComOp%file_HADA%unit
          END IF

        ELSE  ! unformatted file
          read(nio) nb_thread_file
          IF (debug) write(out_unitp,*) 'init,nb_thread_file: ',nb_thread_file

          CALL file_close(ComOp%file_HADA)

          ! for parallel calculation of SH_HADA file (the file is split in SH_HADA.0, SH_HADA.1 ...)
          IF (nb_thread_file > 1) THEN
            !write(out_unitp,*) ' OMP calc of HADA file'
            ithread      = 0
            file_is_para = .TRUE.

            IF (debug) write(out_unitp,*) 'file_is_para,nb_thread_file',file_is_para,nb_thread_file
            ComOp%file_HADA%nb_thread = nb_thread_file
          ELSE
            ComOp%file_HADA%nb_thread = 0
            file_is_para = .FALSE.
          END IF

          CALL file_open(ComOp%file_HADA,nio,                           &
                         lformatted=ComOp%file_HADA%formatted)

          IF (file_is_para) THEN
            ithread = 0
            nio = ComOp%file_HADA%tab_unit(ithread)
          ELSE
            nio = ComOp%file_HADA%unit
          END IF

        END IF
      END IF



      IF (ComOp%file_HADA%formatted) THEN

!       - read parameters at iq --------------------------------
        read(nio,*,iostat=iocond) name1,name2
        IF (debug) write(out_unitp,*) 'iq,iocond,name2',iq,iocond,name2
        IF (debug) write(out_unitp,*) 'iq,ithread,nio',iq,ithread,nio

        IF (iocond > 0 .OR. iocond < 0 .AND. .NOT. file_is_para) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' Problem with the SH_HADA file'
          STOP
        ELSE IF (iocond < 0 .AND. file_is_para) THEN
          ithread = ithread + 1
          nio = ComOp%file_HADA%tab_unit(ithread)
          read(nio,*) name1,name2
        END IF


        IF (trim(name2) /= 'Beginning') THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' The beginning of the record is not correct'
          write(out_unitp,*)  trim(name2),' instead of "Beginning"'
          write(out_unitp,*) ' Probably, you should restart with Read_Grid=f'
          STOP
        END IF


        read(nio,*) name1,iqr
        !write(out_unitp,*) 'name1,iqr',name1,iqr
        IF (iqr /= iq) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' iq of HADA file .NE. iq prog',iqr,iq
          write(out_unitp,*) ' Check restart parameter: restart and num_grid'
          write(out_unitp,*) '  OR'
          write(out_unitp,*) ' Restart the calculation with Read_Grid=f'
          STOP
        END IF
        read(nio,*) name1,n3
!       write(out_unitp,*) name1,n3
        IF (n3 /= d0MatOp%nb_bie) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' iq=',iq
          write(out_unitp,*) ' nb_bie(file) /= nb_bie(data)',n3,d0MatOp%nb_bie
          write(out_unitp,*) ' Restart the calculation with Read_Grid=f'
          STOP
        END IF


        read(nio,*) name1,w
!       write(out_unitp,*) name1,w

        read(nio,*) name1,n1
!       write(out_unitp,*) name1,n1
        IF (n1 /= nb_var) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nb_var of HADA file .NE. nb_var of data',        &
                      n1,nb_var
          write(out_unitp,*) ' Check the HADA file'
          STOP
        END IF
        CALL Read_RVec(Qdyn,nio,5,err)
        IF (err /= 0) THEN
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) ' reading the vector "Qdyn"'
          write(out_unitp,*) ' grid point, iq:',iq
          STOP
        END IF

        read(nio,*) name1,n2
        IF (n2 /= ComOp%nb_act1) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nb_act1 of HADA file .NE. nb_act1 of data',      &
                       n2,ComOp%nb_act1
          write(out_unitp,*) ' Check the HADA file'
          STOP
        END IF
        CALL Read_RVec(Qact,nio,5,err)
        IF (err /= 0) THEN
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) ' reading the vector "Qact"'
          write(out_unitp,*) ' grid point, iq:',iq
          STOP
        END IF
        !write(out_unitp,*) 'Qact',Qact

        read(nio,*) name1,JJ
        read(nio,*) name1,pot_cplx
        read(nio,*) name1,calc_scalar_Op
        read(nio,*) name1,nb_Op
!       write(out_unitp,*) name1,nb_Op
!       ---------------------------------------------------------


        DO i=1,nb_Op

          !write(out_unitp,*) name_sub,',iqr,i,nb_Op,i_term',iqr,i,nb_Op,i_term
          !call flush_perso(out_unitp)
          read(nio,*) name1,n_OP_lect,id1,id2,cplx
          !write(out_unitp,*) name_sub,name1,n_OP_lect,derive_lect(:),cplx
          !call flush_perso(out_unitp)

          CALL Read_RMat(work_bhe,nio,5,err)
          IF (err /= 0) THEN
            write(out_unitp,*) 'ERROR in ',name_sub
            write(out_unitp,*) ' reading the matrix "work_bhe"'
            write(out_unitp,*) ' grid point, iq:',iq
            STOP
          END IF


          IF (n_OP_lect == n_OP) THEN
            IF (cplx) THEN
              d0MatOp%Imval(:,:) = work_bhe(:,:)
            ELSE
              i_term = d0MatOp%derive_term_TO_iterm(id1,id2)
              d0MatOp%ReVal(:,:,i_term) = work_bhe(:,:)
            END IF
          END IF
!         write(out_unitp,*) name_sub,', i_term',i_term


!IF (i_term == 1) write(777,*) Qact,min(work_bhe(1,1),0.05_Rkind)


        END DO

        read(nio,*) name1,name2
        IF (trim(name2) /= 'End') THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' The end of the record is not correct'
          write(out_unitp,*)  trim(name2),' instead of "End"'
          write(out_unitp,*) ' Probably, you should restart with Read_Grid=f'
          STOP
        END IF

      ELSE

!       - read parameters at iq --------------------------------
        read(nio,iostat=iocond) n3
!        write(out_unitp,*) 'iq,iocond,nb_thread',iq,iocond,n3
!        IF (n3 > 1) THEN
!          write(out_unitp,*) 'name_SHADA',ComOp%file_HADA%tab_name_th(ithread)
!        ELSE
!          write(out_unitp,*) 'name_SHADA',ComOp%file_HADA%name
!        END IF
!        CALL flush_perso(out_unitp)
        IF (iocond > 0 .OR. iocond < 0 .AND. .NOT. file_is_para) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' Problem with the SH_HADA file'
          STOP
        ELSE IF (iocond < 0 .AND. file_is_para) THEN
          ithread = ithread + 1
          nio = ComOp%file_HADA%tab_unit(ithread)
          read(nio,iostat=iocond) n3
        END IF

        read(nio) iqr
!       write(out_unitp,*) iqr,iq
        IF (iqr /= iq) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' iq of HADA file .NE. iq prog',iqr,iq
          write(out_unitp,*) ' Check restart parameter: restart and num_grid'
          write(out_unitp,*) '  OR'
          write(out_unitp,*) ' Restart the calculation with Read_Grid=f'
          STOP
        END IF
        read(nio) n3
!       write(out_unitp,*) n3,nb_bie
        IF (n3 /= d0MatOp%nb_bie) THEN
          write(out_unitp,*) ' ERROR in sub_reading_Mat'
          write(out_unitp,*) ' iq=',iq
          write(out_unitp,*) ' nb_bie(file) .NE. nb_bie(data)',n3,d0MatOp%nb_bie
          write(out_unitp,*) ' Restart the calculation with Read_Grid=f'
          STOP
        END IF

        read(nio) w
!       write(out_unitp,*) w

        read(nio) n1
!       write(out_unitp,*) n1,nb_var
        IF (n1 /= nb_var) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nb_var of HADA file .NE. nb_var of data',        &
                      n1,nb_var
          write(out_unitp,*) ' Check the HADA file'
          STOP
        END IF
        read(nio) Qdyn
!       write(out_unitp,*) Qdyn
        read(nio) n2
!       write(out_unitp,*) n2,nb_act1
        IF (n2 /= ComOp%nb_act1) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nb_act1 of HADA file .NE. nb_act1 of data',      &
                       n2,ComOp%nb_act1
          write(out_unitp,*) ' Check the HADA file'
          STOP
        END IF
        read(nio) Qact
!       write(out_unitp,*) Qact


        read(nio) JJ
!       write(out_unitp,*) JJ
        read(nio) pot_cplx
!       write(out_unitp,*) pot_cplx
        read(nio) calc_scalar_Op
!       write(out_unitp,*) calc_scalar_Op
        read(nio) nb_Op
!       write(out_unitp,*) nb_Op
!       ---------------------------------------------------------

        DO i=1,nb_Op

          read(nio) n_OP_lect,id1,id2,cplx
!         write(out_unitp,*) name_sub,n_OP_lect,id1,id2,cplx

          read(nio) work_bhe

          IF (n_OP_lect == n_OP) THEN
            IF (cplx) THEN
              d0MatOp%ImVal(:,:) = work_bhe(:,:)
            ELSE
              i_term = d0MatOp%derive_term_TO_iterm(id1,id2)
              d0MatOp%ReVal(:,:,i_term) = work_bhe(:,:)
            END IF
          END IF

        END DO

        read(nio) n3
      END IF


      IF (iq == nb_qa) CALL file_close(ComOp%file_HADA)

      CALL dealloc_NParray(work_bhe,'work_bhe',name_sub)

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Qact,w',Qact,w
        write(out_unitp,*) 'Op',n_Op
        CALL Write_d0MatOp(d0MatOp)
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_reading_Op

      SUBROUTINE sub_ReadDir_TO_SaveSeq_Grid(para_AllOp)
      USE mod_system
      USE mod_Op
      IMPLICIT NONE

      TYPE (param_AllOp) :: para_AllOp


      integer :: i_qa,nb_qa,nb_bie,iOp,iterm,lrecl_Grid_iterm,nio,error

      real (kind=Rkind), allocatable :: Grid(:,:,:) ! grid when Save_Grid_iterm=t


      nb_qa  = para_AllOp%tab_Op(1)%nb_qa
      nb_bie = para_AllOp%tab_Op(1)%nb_bie

      CALL alloc_NParray(Grid,(/ nb_qa,nb_bie,nb_bie /),                  &
                                   'Grid','sub_ReadDir_TO_SaveSeq_Grid')

      !-- Read / write the grids ---------------------------------------
      DO iOp=1,para_AllOp%nb_Op

        IF (para_AllOp%tab_Op(iOp)%n_Op == -1) EXIT

        DO iterm=1,para_AllOp%tab_Op(iOp)%nb_term

          IF (para_AllOp%tab_Op(iOp)%OpGrid(iterm)%Grid_cte) CYCLE


          !-- Read the direct acces file, then delete it
          lrecl_Grid_iterm = para_AllOp%tab_Op(iOp)%OpGrid(iterm)%file_Grid%frecl
          CALL file_open(para_AllOp%tab_Op(iOp)%OpGrid(iterm)%file_Grid,&
                         nio,lformatted=.FALSE.,seq=.FALSE.,            &
                         lrecl=lrecl_Grid_iterm)
          DO i_qa=1,nb_qa
            read(nio,REC=i_qa,iostat=error) Grid(i_qa,:,:)

            IF (error /= 0) THEN
              write(out_unitp,*) ' ERROR in sub_ReadDir_TO_SaveSeq_Grid'
              write(out_unitp,*) ' Impossible to read the file: ',      &
                     para_AllOp%tab_Op(iOp)%OpGrid(iterm)%file_Grid%name
              STOP
            END IF
          END DO
          CALL file_delete(para_AllOp%tab_Op(iOp)%OpGrid(iterm)%file_Grid)


          !-- Open and write the sequential acces file, then close it
          para_AllOp%tab_Op(iOp)%OpGrid(iterm)%file_Grid%frecl = 0
          para_AllOp%tab_Op(iOp)%OpGrid(iterm)%file_Grid%seq   = .TRUE.
          CALL file_open(para_AllOp%tab_Op(iOp)%OpGrid(iterm)%file_Grid,&
                         nio,lformatted=.FALSE.)

          DO i_qa=1,nb_qa
            write(nio,iostat=error) Grid(i_qa,:,:)
            IF (error /= 0) THEN
              write(out_unitp,*) ' ERROR in sub_ReadDir_TO_SaveSeq_Grid'
              write(out_unitp,*) ' Impossible to write the file: ',     &
                    para_AllOp%tab_Op(iOp)%OpGrid(iterm)%file_Grid%name
              STOP
            END IF
          END DO

          CALL file_close(para_AllOp%tab_Op(iOp)%OpGrid(iterm)%file_Grid)

        END DO

      END DO

      CALL dealloc_NParray(Grid,'Grid','sub_ReadDir_TO_SaveSeq_Grid')

      END SUBROUTINE sub_ReadDir_TO_SaveSeq_Grid

!=====================================================================
!
!  calculation of d0bnD d1bnD(:) and d2bnD(:,:)
!  with nb_act1 derivatives
!
!=====================================================================
      SUBROUTINE d0d1d2bnDQact(d0b,d1b,d2b,BasisnD,iq,ib,mole)
      USE mod_system
      USE mod_basis
      USE mod_Tnum
      implicit none

      !----- variables for the Basis and quadrature points -----------------
      TYPE (Basis) :: BasisnD

      !----- for the zmatrix  --------------------------------------
      TYPE (zmatrix) :: mole

      integer           :: iq,ib
      real (kind=Rkind) :: d0b
      real (kind=Rkind) :: d1b(mole%nb_act1)
      real (kind=Rkind) :: d2b(mole%nb_act1,mole%nb_act1)

      !------ working variables ---------------------------------
      integer           :: i,j
      real (kind=Rkind) :: d0b_loc
      real (kind=Rkind), allocatable  :: d1b_loc(:)
      real (kind=Rkind), allocatable  :: d2b_loc(:,:)
      integer, allocatable            :: iQact(:)

       !----- for debuging --------------------------------------------------
       logical, parameter :: debug = .FALSE.
       !logical, parameter :: debug = .TRUE.
       !-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING d0d1d2bnDQact'
         write(out_unitp,*) 'iq,ib',iq,ib
       END IF
       !-----------------------------------------------------------

       CALL alloc_NParray(iQact,  (/mole%nb_act1/),'iQact',  'd0d1d2bnDQact')
       CALL alloc_NParray(d1b_loc,shape(d1b),      'd1b_loc','d0d1d2bnDQact')
       CALL alloc_NParray(d2b_loc,shape(d2b),      'd2b_loc','d0d1d2bnDQact')

       ! in Rec_d0d1d2bnD, the derivatives are in basis order
       CALL Rec_d0d1d2bnD(d0b_loc,d1b_loc,d2b_loc,BasisnD,iq,ib)

       iQact(:) = mole%ActiveTransfo%list_QdynTOQact(BasisnD%iQdyn(:))

       ! Put the derivatives in Qact order
       d0b      = d0b_loc
       DO i=1,mole%nb_act1
         d1b(iQact(i)) = d1b_loc(i)
         DO j=1,mole%nb_act1
           d2b(iQact(i),iQact(j)) = d2b_loc(i,j)
         END DO
       END DO

       CALL dealloc_NParray(iQact,  'iQact',  'd0d1d2bnDQact')
       CALL dealloc_NParray(d1b_loc,'d1b_loc','d0d1d2bnDQact')
       CALL dealloc_NParray(d2b_loc,'d2b_loc','d0d1d2bnDQact')

      !-------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) ' d0b',d0b
        write(out_unitp,*) ' d1b',d1b
        write(out_unitp,*) ' d2b',d2b
        write(out_unitp,*)
        write(out_unitp,*) 'END d0d1d2bnDQact'
      END IF
      !-------------------------------------------------------


      END SUBROUTINE d0d1d2bnDQact

!=====================================================================
!
!  check : read the HADA file and determine the last grid point (iqf)
!
!=====================================================================
      SUBROUTINE check_HADA(iqf,ComOp)

      USE mod_system
      USE mod_Op
      IMPLICIT NONE

      TYPE (param_ComOp)   :: ComOp


      integer :: i,iq,iqf
      integer :: ios

      integer             :: nio
      character (len=132) :: name
      character (len=3)   :: name3
      character (len=Name_len)  :: nameA20


      CALL file_open(ComOp%file_HADA,nio,                               &
                     lformatted=ComOp%file_HADA%formatted)

      iqf = 0

      IF (ComOp%file_HADA%formatted) THEN

 10     CONTINUE

        read(nio,11,END=20) name
 11     format(A132)
!       write(out_unitp,*) name
        read(name,*,IOSTAT=ios) name3,iq
        IF (ios /= 0) GOTO 10

        IF (name3 == "iq:") iqf=iq
!       write(out_unitp,*) name3,iqf

        GOTO 10

 20     CONTINUE

      ELSE
        nameA20 = 'iq:'

        i = 0
 50     CONTINUE
        i = i + 1

        read(nio,END=60,ERR=50) name3,iq

        IF (name3 == nameA20) iqf=iq
!       write(out_unitp,*) name3,iq,iqf

        GOTO 50

 60     CONTINUE

      END IF



      close(ComOp%file_HADA%unit)


      write(out_unitp,*) 'check_HADA: last grid point, iqf=',iqf

      end subroutine check_HADA

