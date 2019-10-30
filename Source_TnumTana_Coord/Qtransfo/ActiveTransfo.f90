!===========================================================================
!===========================================================================
!This file is part of Tnum-Tana.
!
!    Tnum-Tana is a free software: you can redistribute it and/or modify
!    it under the terms of the GNU Lesser General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    Tnum-Tana is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Lesser General Public License for more details.
!
!    You should have received a copy of the GNU Lesser General Public License
!    along with ElVibRot.  If not, see <http://www.gnu.org/licenses/>.
!
!    Copyright 2015  David Lauvergnat
!      with contributions of Mamadou Ndong
!
!===========================================================================
!===========================================================================
MODULE mod_ActiveTransfo
      use mod_system
      use mod_dnSVM, only: alloc_array, dealloc_array, type_dnvec,   &
                           type_dns, write_dnsvm, alloc_dnsvm,       &
                           set_zero_to_dnsvm, sub_dns_to_dnvec,      &
                           dealloc_dnsvm
      IMPLICIT NONE

      PRIVATE

      !! @description: TODO
      !! @param: nb_var TODO
      !! @param: TODO
      !! @param: TODO
      TYPE Type_ActiveTransfo
          integer :: nb_var      = 0
          integer :: nb_act      = 0
          integer :: nb_act1     = 0
          integer :: nb_inact2n  = 0
          integer :: nb_inact21  = 0
          integer :: nb_inact22  = 0
          integer :: nb_inact20  = 0
          integer :: nb_inact    = 0
          integer :: nb_inact31  = 0
          integer :: nb_rigid0   = 0
          integer :: nb_rigid100 = 0
          integer :: nb_rigid    = 0
          real (kind=Rkind), pointer  :: Qdyn0(:)            => null() ! value of rigid coordinates (Qdyn order)
          real (kind=Rkind), pointer  :: Qact0(:)            => null() ! value of rigid coordinates (Qact order)
          integer, pointer            :: list_act_OF_Qdyn(:) => null()  ! "active" transfo
          integer, pointer            :: list_QactTOQdyn(:)  => null() ! "active" transfo
          integer, pointer            :: list_QdynTOQact(:)  => null() ! "active" transfo
      END TYPE Type_ActiveTransfo

      INTERFACE alloc_array
        ! for RPHTransfo
        MODULE PROCEDURE alloc_array_OF_ActiveTransfodim0
      END INTERFACE
      INTERFACE dealloc_array
        ! for RPHTransfo
        MODULE PROCEDURE dealloc_array_OF_ActiveTransfodim0
      END INTERFACE

      PUBLIC :: Type_ActiveTransfo, alloc_ActiveTransfo, dealloc_ActiveTransfo, ActiveTransfo1TOActiveTransfo2
      PUBLIC :: alloc_array, dealloc_array
      PUBLIC :: Read_ActiveTransfo, Read2_ActiveTransfo, Write_ActiveTransfo
      PUBLIC :: calc_ActiveTransfo
      PUBLIC :: get_Qact, get_Qact0, Set_AllActive
      PUBLIC :: Qact_TO_Qdyn_FROM_ActiveTransfo, Qdyn_TO_Qact_FROM_ActiveTransfo, Qinact2n_TO_Qact_FROM_ActiveTransfo

      CONTAINS

!================================================================
!      Subroutines for the Active Transfo:
!       alloc_ActiveTransfo
!       dealloc_ActiveTransfo
!       Read_ActiveTransfo
!       Check_ActiveTransfo
!       calc_Activetransfo
!================================================================
      !!@description:  Subroutines for the Active Transfo:
      !!       alloc_ActiveTransfo
      !!       dealloc_ActiveTransfo
      !!       Read_ActiveTransfo
      !!       Check_ActiveTransfo
      !!       calc_Activetransfo
      !!@param: TODO
      SUBROUTINE alloc_ActiveTransfo(ActiveTransfo,nb_var)

      TYPE (Type_ActiveTransfo), intent(inout) :: ActiveTransfo
      integer, intent(in) :: nb_var

      character (len=*), parameter :: name_sub='alloc_ActiveTransfo'

      ActiveTransfo%nb_var = nb_var

      CALL alloc_array(ActiveTransfo%list_act_OF_Qdyn,(/nb_var/),       &
                      "ActiveTransfo%list_act_OF_Qdyn",name_sub)
      ActiveTransfo%list_act_OF_Qdyn(:) = 0
      CALL alloc_array(ActiveTransfo%list_QactTOQdyn,(/nb_var/),        &
                      "ActiveTransfo%list_QactTOQdyn",name_sub)
      ActiveTransfo%list_QactTOQdyn(:) = 0
      CALL alloc_array(ActiveTransfo%list_QdynTOQact,(/nb_var/),        &
                      "ActiveTransfo%list_QdynTOQact",name_sub)
      ActiveTransfo%list_QdynTOQact(:) = 0

      CALL alloc_array(ActiveTransfo%Qdyn0,(/nb_var/),                  &
                      "ActiveTransfo%Qdyn0",name_sub)
      ActiveTransfo%Qdyn0(:) = ZERO

      CALL alloc_array(ActiveTransfo%Qact0,(/nb_var/),                  &
                      "ActiveTransfo%Qact0",name_sub)
      ActiveTransfo%Qact0(:) = ZERO

      END SUBROUTINE alloc_ActiveTransfo
      !-----------------------------------------------------------------------

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE dealloc_ActiveTransfo(ActiveTransfo)

      TYPE (Type_ActiveTransfo), intent(inout) :: ActiveTransfo

      character (len=*), parameter :: name_sub='dealloc_ActiveTransfo'


      IF (associated(ActiveTransfo%list_act_OF_Qdyn)) THEN
        CALL dealloc_array(ActiveTransfo%list_act_OF_Qdyn,              &
                          "ActiveTransfo%list_act_OF_Qdyn",name_sub)
      END IF
      IF (associated(ActiveTransfo%list_QactTOQdyn)) THEN
        CALL dealloc_array(ActiveTransfo%list_QactTOQdyn,               &
                          "ActiveTransfo%list_QactTOQdyn",name_sub)
      END IF
      IF (associated(ActiveTransfo%list_QdynTOQact)) THEN
        CALL dealloc_array(ActiveTransfo%list_QdynTOQact,               &
                          "ActiveTransfo%list_QdynTOQact",name_sub)
      END IF

      IF (associated(ActiveTransfo%Qdyn0))  THEN
        CALL dealloc_array(ActiveTransfo%Qdyn0,                         &
                          "ActiveTransfo%Qdyn0",name_sub)
      END IF

      IF (associated(ActiveTransfo%Qact0))  THEN
        CALL dealloc_array(ActiveTransfo%Qact0,                         &
                          "ActiveTransfo%Qact0",name_sub)
      END IF

      ActiveTransfo%nb_var      = 0
      ActiveTransfo%nb_act      = 0
      ActiveTransfo%nb_act1     = 0
      ActiveTransfo%nb_inact2n  = 0
      ActiveTransfo%nb_inact2n  = 0
      ActiveTransfo%nb_inact2n  = 0
      ActiveTransfo%nb_inact20  = 0
      ActiveTransfo%nb_inact    = 0
      ActiveTransfo%nb_inact31  = 0
      ActiveTransfo%nb_rigid0   = 0
      ActiveTransfo%nb_rigid100 = 0
      ActiveTransfo%nb_rigid    = 0

      END SUBROUTINE dealloc_ActiveTransfo

    SUBROUTINE alloc_array_OF_ActiveTransfodim0(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (Type_ActiveTransfo), pointer, intent(inout) :: tab

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=0
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_ActiveTransfodim0'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (associated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       memory = 1
       allocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Type_ActiveTransfo')

      END SUBROUTINE alloc_array_OF_ActiveTransfodim0
      SUBROUTINE dealloc_array_OF_ActiveTransfodim0(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (Type_ActiveTransfo), pointer, intent(inout) :: tab
      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_ActiveTransfodim0'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN
       IF (.NOT. associated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = 1
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Type_ActiveTransfo')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_ActiveTransfodim0


      !!@description: TODO
      !!@param: TODO
      SUBROUTINE Read_ActiveTransfo(ActiveTransfo,nb_Qin)
      USE mod_MPI

      TYPE (Type_ActiveTransfo), intent(inout) :: ActiveTransfo
      integer, intent(in) :: nb_Qin

      integer :: err
      character (len=*), parameter :: name_sub='Read_ActiveTransfo'

      CALL alloc_ActiveTransfo(ActiveTransfo,nb_Qin)

      read(in_unitp,*,IOSTAT=err) ActiveTransfo%list_act_OF_Qdyn(:)
      IF(MPI_id==0) write(out_unitp,*) 'list_act_OF_Qdyn or type_var',                 &
                                        ActiveTransfo%list_act_OF_Qdyn(:)
      IF (err /= 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) '  while reading "list_act_OF_Qdyn"'
        write(out_unitp,*) ' end of file or end of record'
        write(out_unitp,*) ' Check your data !!'
        STOP
      END IF


      CALL flush_perso(out_unitp)

      END SUBROUTINE Read_ActiveTransfo

      SUBROUTINE Read2_ActiveTransfo(ActiveTransfo,nb_Qin)

      TYPE (Type_ActiveTransfo), intent(inout) :: ActiveTransfo
      integer, intent(in) :: nb_Qin

      character (len=Name_len) :: name_int
      integer :: i,nb_Qact
      integer, pointer :: list_Qact(:)

      integer :: err_io
      character (len=*), parameter :: name_sub='Read2_ActiveTransfo'

      CALL alloc_ActiveTransfo(ActiveTransfo,nb_Qin)

      read(in_unitp,*,IOSTAT=err_io) ActiveTransfo%list_act_OF_Qdyn(:)
      IF(MPI_id==0) write(out_unitp,*) 'list_act_OF_Qdyn or type_var',                 &
                                        ActiveTransfo%list_act_OF_Qdyn(:)
      IF (err_io /= 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) '  while reading "list_act_OF_Qdyn"'
        write(out_unitp,*) ' end of file or end of record'
        write(out_unitp,*) ' Check your data !!'
        STOP
      END IF

      IF (count(ActiveTransfo%list_act_OF_Qdyn(:) == 1) /=0) THEN
        write(out_unitp,*) ' WARNNING in ',name_sub
        write(out_unitp,*) ' You have already defined active coordinates in'
        write(out_unitp,*) 'list_act_OF_Qdyn(:)',ActiveTransfo%list_act_OF_Qdyn
      END IF

      nullify(list_Qact)
      CALL alloc_array(list_Qact,(/nb_Qin/),'list_Qact',name_sub)
      list_Qact(:) = 0

      DO i=1,nb_Qin
        CALL read_name_advNo(in_unitp,name_int,err_io)

        IF (len_trim(name_int) == 0) EXIT
        !write(out_unitp,*) 'i,err_io',i,err_io
        !write(out_unitp,*) 'i,name_int',i,name_int
        read(name_int,*) list_Qact(i)
        IF (err_io /= 0) EXIT ! end of the liste

      END DO
      write(out_unitp,*) 'list_Qact_order',list_Qact(:)
      CALL flush_perso(out_unitp)

      ! modify ActiveTransfo%list_act_OF_Qdyn with list_Qact
      DO i=1,count(list_Qact(:) > 0)
        ActiveTransfo%list_act_OF_Qdyn(list_Qact(i)) = 1
      END DO
      write(out_unitp,*) 'New list_act_OF_Qdyn(:)',ActiveTransfo%list_act_OF_Qdyn


      IF (count(ActiveTransfo%list_act_OF_Qdyn(:) == 1) == 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' There is no active coordiantes!'
        write(out_unitp,*) 'list_act_OF_Qdyn(:)',ActiveTransfo%list_act_OF_Qdyn
        write(out_unitp,*) 'Check your data!'
        STOP
      END IF


      CALL dealloc_array(list_Qact,'list_Qact',name_sub)


      CALL flush_perso(out_unitp)

      END SUBROUTINE Read2_ActiveTransfo

      SUBROUTINE Write_ActiveTransfo(ActiveTransfo)

      TYPE (Type_ActiveTransfo), pointer, intent(in) :: ActiveTransfo

      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Write_ActiveTransfo'

      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*) 'asso ActiveTransfo:',associated(ActiveTransfo)
      IF (associated(ActiveTransfo) .AND. MPI_id==0) THEN
        write(out_unitp,*) 'nb_var:        ',ActiveTransfo%nb_var
        write(out_unitp,*) 'nb_act:        ',ActiveTransfo%nb_act
        write(out_unitp,*) 'nb_act1:       ',ActiveTransfo%nb_act1
        write(out_unitp,*) 'nb_inact2n:    ',ActiveTransfo%nb_inact2n
        write(out_unitp,*) 'nb_inact21:    ',ActiveTransfo%nb_inact21
        write(out_unitp,*) 'nb_inact22:    ',ActiveTransfo%nb_inact22
        write(out_unitp,*) 'nb_inact20:    ',ActiveTransfo%nb_inact20
        write(out_unitp,*) 'nb_inact:      ',ActiveTransfo%nb_inact
        write(out_unitp,*) 'nb_inact31:    ',ActiveTransfo%nb_inact31
        write(out_unitp,*) 'nb_rigid0:     ',ActiveTransfo%nb_rigid0
        write(out_unitp,*) 'nb_rigid100:   ',ActiveTransfo%nb_rigid100
        write(out_unitp,*) 'nb_rigid:      ',ActiveTransfo%nb_rigid

        IF (associated(ActiveTransfo%list_act_OF_Qdyn)) THEN
          write(out_unitp,*) 'list_act_OF_Qdyn or type_var: ',ActiveTransfo%list_act_OF_Qdyn(:)
        ELSE
          write(out_unitp,*) 'asso list_act_OF_Qdyn?   F'
        END IF

        IF (associated(ActiveTransfo%list_QactTOQdyn)) THEN
          write(out_unitp,*) 'list_QactTOQdyn: ',ActiveTransfo%list_QactTOQdyn(:)
        ELSE
          write(out_unitp,*) 'asso list_QactTOQdyn?   F'
        END IF

        IF (associated(ActiveTransfo%list_QdynTOQact)) THEN
          write(out_unitp,*) 'list_QdynTOQact: ',ActiveTransfo%list_QdynTOQact(:)
        ELSE
          write(out_unitp,*) 'asso list_QdynTOQact?   F'
        END IF

        IF (associated(ActiveTransfo%Qdyn0)) THEN
          write(out_unitp,*) ' Rigid coordinate values (Qdyn order):',  &
                                                 ActiveTransfo%Qdyn0(:)
        ELSE
          write(out_unitp,*) 'asso Qdyn0?   F'
        END IF

        IF (associated(ActiveTransfo%Qact0)) THEN
          write(out_unitp,*) 'Qact0: ',ActiveTransfo%Qact0(:)
        ELSE
          write(out_unitp,*) 'asso Qact0?   F'
        END IF

      END IF
      write(out_unitp,*) 'END ',name_sub

      END SUBROUTINE Write_ActiveTransfo

      SUBROUTINE calc_ActiveTransfo(dnQact,dnQdyn,ActiveTransfo,nderiv,inTOout)
      IMPLICIT NONE

        TYPE (Type_dnVec), intent(inout)        :: dnQact,dnQdyn
        TYPE (Type_ActiveTransfo), intent(in)   :: ActiveTransfo
        integer, intent(in)                     :: nderiv
        logical                                 :: inTOout


        TYPE (Type_dnS)    :: dnQ
        integer :: typ_var_act,i_Qdyn,i_Qact,nb_act1


!      -----------------------------------------------------------------
!       logical, parameter :: debug=.TRUE.
       logical, parameter :: debug=.FALSE.
       character (len=*), parameter :: name_sub='calc_ActiveTransfo'
!      -----------------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'nderiv',nderiv
         write(out_unitp,*) 'nb_var_deriv',dnQact%nb_var_deriv
         write(out_unitp,*) 'nb_act',ActiveTransfo%nb_act
         write(out_unitp,*) 'asso Qact0 ?',associated(ActiveTransfo%Qact0)
         IF (inTOOut) THEN
           write(out_unitp,*) 'dnQact'
           CALL Write_dnSVM(dnQact)
         ELSE
           write(out_unitp,*) 'dnQdyn'
           CALL Write_dnSVM(dnQdyn)
         END IF
         write(out_unitp,*)
       END IF
!      -----------------------------------------------------------------

       dnQ%nb_var_deriv = dnQact%nb_var_deriv
       dnQ%nderiv       = nderiv
       nb_act1          = ActiveTransfo%nb_act1

       IF (inTOout) THEN ! Qact => Qdyn (with the derivatives)
         CALL alloc_dnSVM(dnQ)

         DO i_Qact=1,ActiveTransfo%nb_var

           CALL Set_ZERO_TO_dnSVM(dnQ,nderiv)

           i_Qdyn      = ActiveTransfo%list_QactTOQdyn(i_Qact)
           typ_var_act = ActiveTransfo%list_act_OF_Qdyn(i_Qdyn)

           SELECT CASE (typ_var_act)
           CASE (1,-1,21,22,31)
             ! active coordinate
             dnQ%d0 = dnQact%d0(i_Qact)
             IF ( nderiv >= 1 ) dnQ%d1(i_Qact) = ONE
           CASE (20)
             ! inactive coordinate : flexible constraints
             CALL calc_dnQflex(i_Qdyn,dnQ,dnQact%d0(1:nb_act1),nb_act1,nderiv,-1)
           CASE (200)
             ! inactive coordinate : flexible constraints
             ! nderiv MUST be 0
             CALL calc_dnQflex(i_Qdyn,dnQ,dnQact%d0(1:nb_act1),nb_act1,0,-1)
           CASE (0,100)
             ! inactive coordinate : rigid0 and rigid100
             dnQ%d0 = ActiveTransfo%Qact0(i_Qact)
           CASE default
             write(out_unitp,*) ' ERROR in ',name_sub
             write(out_unitp,*) ' I do not know this variable type:',typ_var_act
             write(out_unitp,*) ' Check your data!!'
             STOP
           END SELECT

           CALL sub_dnS_TO_dnVec(dnQ,dnQact,i_Qact)
           CALL sub_dnS_TO_dnVec(dnQ,dnQdyn,i_Qdyn)

         END DO
         CALL dealloc_dnSVM(dnQ)

       ELSE ! Qdyn => Qact (without the derivatives)
         IF (nderiv > 0) THEN
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) ' you cannot use this subroutine with inTOout=f and nderiv>0'
           write(out_unitp,*) '    Qdyn => Qact  '
           write(out_unitp,*) ' Check your the fortran!!'
           STOP
         END IF

         CALL Qdyn_TO_Qact_FROM_ActiveTransfo(dnQdyn%d0,dnQact%d0,ActiveTransfo)

       END IF


!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'dnQdyn'
        CALL Write_dnSVM(dnQdyn)
        write(out_unitp,*) 'dnQact'
        CALL Write_dnSVM(dnQact)
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -----------------------------------------------------------------

      END SUBROUTINE calc_ActiveTransfo

      SUBROUTINE get_Qact(Qact,ActiveTransfo,With_All)
      IMPLICIT NONE

        real (kind=Rkind), intent(inout) :: Qact(:)
        TYPE (Type_ActiveTransfo), intent(in)   :: ActiveTransfo
        logical, intent(in), optional :: With_All


        TYPE (Type_dnS)    :: dnQ
        integer :: typ_var_act,i_Qdyn,i_Qact,nb_act1
        logical :: With_All_loc


!      -----------------------------------------------------------------
!       logical, parameter :: debug=.TRUE.
       logical, parameter :: debug=.FALSE.
       character (len=*), parameter :: name_sub='get_Qact'
!      -----------------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'nb_act',ActiveTransfo%nb_act
         write(out_unitp,*) 'asso Qact0 ?',associated(ActiveTransfo%Qact0)
         write(out_unitp,*) 'Qact',Qact(:)
         write(out_unitp,*)
         CALL flush_perso(out_unitp)
       END IF
!      -----------------------------------------------------------------

       dnQ%nb_var_deriv = ActiveTransfo%nb_act
       dnQ%nderiv       = 0
       nb_act1          = ActiveTransfo%nb_act1

       IF (present(With_All)) THEN
         With_All_loc = With_All
       ELSE
         With_All_loc = .TRUE.
       END IF

       CALL alloc_dnSVM(dnQ)

       DO i_Qact=1,ActiveTransfo%nb_var

         CALL Set_ZERO_TO_dnSVM(dnQ,nderiv=0)

         i_Qdyn      = ActiveTransfo%list_QactTOQdyn(i_Qact)
         typ_var_act = ActiveTransfo%list_act_OF_Qdyn(i_Qdyn)

         SELECT CASE (typ_var_act)
         CASE (1,-1,21,22,31)
           ! active coordinate, nothing here, because it the Qact coord
           ! except if With_All_loc=.TRUE.
           IF (With_All_loc) Qact(i_Qact) = ActiveTransfo%Qact0(i_Qact)

         CASE (20)
           ! inactive coordinate : flexible constraints
           CALL calc_dnQflex(i_Qdyn,dnQ,Qact(1:nb_act1),nb_act1,0,-1)
           Qact(i_Qact) = dnQ%d0
         CASE (200)
           ! inactive coordinate : flexible constraints
           ! nderiv MUST be 0
           CALL calc_dnQflex(i_Qdyn,dnQ,Qact(1:nb_act1),nb_act1,0,-1)
           Qact(i_Qact) = dnQ%d0

         CASE (0,100)
           ! inactive coordinate : rigid0 and rigid100
           Qact(i_Qact) = ActiveTransfo%Qact0(i_Qact)
         CASE default
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) ' I do not know this variable type:',typ_var_act
           write(out_unitp,*) ' Check your data!!'
           STOP
         END SELECT

       END DO
       CALL dealloc_dnSVM(dnQ)


!     -----------------------------------------------------------------
      IF (debug) THEN
         write(out_unitp,*) 'Qact',Qact(:)
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -----------------------------------------------------------------

      END SUBROUTINE get_Qact

      SUBROUTINE get_Qact0(Qact0,ActiveTransfo)
      IMPLICIT NONE

        real (kind=Rkind), intent(inout)      :: Qact0(:)
        TYPE (Type_ActiveTransfo), intent(in) :: ActiveTransfo


        TYPE (Type_dnS)    :: dnQ
        integer :: typ_var_act,i_Qdyn,i_Qact,nb_act1


!      -----------------------------------------------------------------
!       logical, parameter :: debug=.TRUE.
       logical, parameter :: debug=.FALSE.
       character (len=*), parameter :: name_sub='get_Qact0'
!      -----------------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'nb_act',ActiveTransfo%nb_act
         write(out_unitp,*) 'asso Qact0 ?',associated(ActiveTransfo%Qact0)
         write(out_unitp,*) 'Qact0',Qact0(:)
         write(out_unitp,*)
         CALL flush_perso(out_unitp)
       END IF
!      -----------------------------------------------------------------

       dnQ%nb_var_deriv = ActiveTransfo%nb_act
       dnQ%nderiv       = 0
       nb_act1          = ActiveTransfo%nb_act1

       CALL alloc_dnSVM(dnQ)

       DO i_Qact=1,size(Qact0)

         CALL Set_ZERO_TO_dnSVM(dnQ,nderiv=0)

         i_Qdyn      = ActiveTransfo%list_QactTOQdyn(i_Qact)
         typ_var_act = ActiveTransfo%list_act_OF_Qdyn(i_Qdyn)

         SELECT CASE (typ_var_act)
         CASE (1,-1,21,22,31)
           ! active coordinate, nothing here, because it the Qact coord
           ! except if With_All_loc=.TRUE.
           Qact0(i_Qact) = ActiveTransfo%Qact0(i_Qact)

         CASE (20)
           ! inactive coordinate : flexible constraints
           CALL calc_dnQflex(i_Qdyn,dnQ,Qact0(1:nb_act1),nb_act1,0,-1)
           Qact0(i_Qact) = dnQ%d0
         CASE (200)
           ! inactive coordinate : flexible constraints
           ! nderiv MUST be 0
           CALL calc_dnQflex(i_Qdyn,dnQ,Qact0(1:nb_act1),nb_act1,0,-1)
           Qact0(i_Qact) = dnQ%d0

         CASE (0,100)
           ! inactive coordinate : rigid0 and rigid100
           Qact0(i_Qact) = ActiveTransfo%Qact0(i_Qact)
         CASE default
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) ' I do not know this variable type:',typ_var_act
           write(out_unitp,*) ' Check your data!!'
           STOP
         END SELECT

       END DO
       CALL dealloc_dnSVM(dnQ)


!     -----------------------------------------------------------------
      IF (debug) THEN
         write(out_unitp,*) 'Qact0',Qact0(:)
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -----------------------------------------------------------------

      END SUBROUTINE get_Qact0

      SUBROUTINE Set_AllActive(dnQact)
      IMPLICIT NONE

        TYPE (Type_dnVec), intent(inout)        :: dnQact


        integer :: i
        real (kind=Rkind) :: Qact(dnQact%nb_var_vec)


!      -----------------------------------------------------------------
!       logical, parameter :: debug=.TRUE.
       logical, parameter :: debug=.FALSE.
       character (len=*), parameter :: name_sub='Set_AllActive'
!      -----------------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'dnQact'
         CALL Write_dnSVM(dnQact)
         write(out_unitp,*)
       END IF
!      -----------------------------------------------------------------

       Qact(:) = dnQact%d0(:)

       CALL Set_ZERO_TO_dnSVM(dnQact)

       IF (dnQact%nderiv > 0) THEN
         DO i=1,dnQact%nb_var_deriv
           dnQact%d1(i,i) = ONE
         END DO
       END IF
       dnQact%d0(:) =        Qact(:)

!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'dnQact'
        CALL Write_dnSVM(dnQact)
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -----------------------------------------------------------------

      END SUBROUTINE Set_AllActive

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE ActiveTransfo1TOActiveTransfo2(ActiveTransfo1,ActiveTransfo2)

!      for the Activerix and Tnum --------------------------------------
      TYPE (Type_ActiveTransfo), intent(in)    :: ActiveTransfo1
      TYPE (Type_ActiveTransfo), intent(inout) :: ActiveTransfo2

      character (len=*), parameter ::                                   &
                             name_sub = 'ActiveTransfo1TOActiveTransfo2'

      CALL dealloc_ActiveTransfo(ActiveTransfo2)

      ActiveTransfo2%nb_var      = ActiveTransfo1%nb_var
      ActiveTransfo2%nb_act      = ActiveTransfo1%nb_act
      ActiveTransfo2%nb_act1     = ActiveTransfo1%nb_act1
      ActiveTransfo2%nb_inact2n  = ActiveTransfo1%nb_inact2n
      ActiveTransfo2%nb_inact21  = ActiveTransfo1%nb_inact21
      ActiveTransfo2%nb_inact22  = ActiveTransfo1%nb_inact22
      ActiveTransfo2%nb_inact20  = ActiveTransfo1%nb_inact20
      ActiveTransfo2%nb_inact    = ActiveTransfo1%nb_inact
      ActiveTransfo2%nb_inact31  = ActiveTransfo1%nb_inact31
      ActiveTransfo2%nb_rigid0   = ActiveTransfo1%nb_rigid0
      ActiveTransfo2%nb_rigid100 = ActiveTransfo1%nb_rigid100
      ActiveTransfo2%nb_rigid    = ActiveTransfo1%nb_rigid

      CALL alloc_ActiveTransfo(ActiveTransfo2,ActiveTransfo1%nb_var)

      IF (associated(ActiveTransfo1%list_act_OF_Qdyn))                  &
        ActiveTransfo2%list_act_OF_Qdyn(:)  = ActiveTransfo1%list_act_OF_Qdyn(:)

      IF (associated(ActiveTransfo1%list_QactTOQdyn))                   &
        ActiveTransfo2%list_QactTOQdyn(:)   = ActiveTransfo1%list_QactTOQdyn(:)

      IF (associated(ActiveTransfo1%list_QdynTOQact))                   &
        ActiveTransfo2%list_QdynTOQact(:)   = ActiveTransfo1%list_QdynTOQact(:)

      IF (associated(ActiveTransfo1%Qdyn0))                             &
        ActiveTransfo2%Qdyn0(:)   = ActiveTransfo1%Qdyn0(:)

      IF (associated(ActiveTransfo1%Qact0))                              &
        ActiveTransfo2%Qact0(:)   = ActiveTransfo1%Qact0(:)

!     write(out_unitp,*) 'END ActiveTransfo1TOActiveTransfo2'

      END SUBROUTINE ActiveTransfo1TOActiveTransfo2

!
!=====================================================================
!
! ++   transfert Qact to Qdyn with the list ActiveTransfo
!       and      Qdyn to Qact with the list ActiveTransfo
!
!=====================================================================
!
      SUBROUTINE Qact_TO_Qdyn_FROM_ActiveTransfo(Qact,Qdyn,ActiveTransfo)
      IMPLICIT NONE


      real (kind=Rkind), intent(in)         :: Qact(:)
      real (kind=Rkind), intent(inout)      :: Qdyn(:)

      TYPE (Type_ActiveTransfo), intent(in) :: ActiveTransfo

!---------------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING Qact_TO_Qdyn_FROM_ActiveTransfo'
        write(out_unitp,*) 'Qact',Qact
      END IF
!---------------------------------------------------------------------

      Qdyn(:) = Qact(ActiveTransfo%list_QdynTOQact)

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Qdyn',Qdyn
        write(out_unitp,*) 'END Qact_TO_Qdyn_FROM_ActiveTransfo'
      END IF
!---------------------------------------------------------------------

      END SUBROUTINE Qact_TO_Qdyn_FROM_ActiveTransfo
      SUBROUTINE Qdyn_TO_Qact_FROM_ActiveTransfo(Qdyn,Qact,ActiveTransfo)
      IMPLICIT NONE

      real (kind=Rkind), intent(in)         :: Qdyn(:)
      real (kind=Rkind), intent(inout)      :: Qact(:)

      TYPE (Type_ActiveTransfo), intent(in) :: ActiveTransfo

!---------------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING Qdyn_TO_Qact_FROM_ActiveTransfo'
        write(out_unitp,*) 'Qdyn',Qdyn
      END IF
!---------------------------------------------------------------------

      Qact(:) = Qdyn(ActiveTransfo%list_QactTOQdyn)

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Qact',Qact
        write(out_unitp,*) 'END Qact_TO_Qdyn_FROM_ActiveTransfo'
      END IF
!---------------------------------------------------------------------

      END SUBROUTINE Qdyn_TO_Qact_FROM_ActiveTransfo
      SUBROUTINE Qinact2n_TO_Qact_FROM_ActiveTransfo(Qinact2n,Qact,ActiveTransfo)
      IMPLICIT NONE


      real (kind=Rkind), intent(inout) :: Qinact2n(:)
      real (kind=Rkind), intent(inout) :: Qact(:)

      TYPE (Type_ActiveTransfo), intent(in)   :: ActiveTransfo

      integer i_Qact,i_Qdyn,i_Q2n
!---------------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING Qinact2n_TO_Qact_FROM_ActiveTransfo'
        write(out_unitp,*) 'Qinact2n',Qinact2n
        write(out_unitp,*) 'Qact',Qact
      END IF
!---------------------------------------------------------------------

      ! we use Qdyn order because
      i_Q2n = 0
      DO i_Qdyn=1,ActiveTransfo%nb_var
        IF (ActiveTransfo%list_act_OF_Qdyn(i_Qdyn) == 21 .OR.           &
            ActiveTransfo%list_act_OF_Qdyn(i_Qdyn) == 22 .OR.           &
            ActiveTransfo%list_act_OF_Qdyn(i_Qdyn) == 31 ) THEN

          i_Q2n  = i_Q2n + 1
          i_Qact = ActiveTransfo%list_QdynTOQact(i_Qdyn)

          IF (i_Q2n > size(Qinact2n)) THEN
            write(out_unitp,*) ' ERROR Qinact2n_TO_Qact_FROM_ActiveTransfo'
            write(out_unitp,*) ' i_Q2n > size(Qinact2n)',i_Q2n,size(Qinact2n)
            STOP
          END IF
          Qact(i_Qact) = Qinact2n(i_Q2n)

        END IF
      END DO

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Qact',Qact
        write(out_unitp,*) 'END Qinact2n_TO_Qact_FROM_ActiveTransfo'
      END IF
!---------------------------------------------------------------------

      END SUBROUTINE Qinact2n_TO_Qact_FROM_ActiveTransfo

END MODULE mod_ActiveTransfo
