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

MODULE mod_Tana_keo
   use mod_system
   use mod_Tnum,   only: zmatrix, tnum, write_mole
   use mod_ActiveTransfo,  only : qact_to_qdyn_from_activetransfo
   USE mod_paramQ
   USE mod_Tana_PiEulerRot
   USE mod_Tana_sum_opnd
   USE mod_Tana_op
   USE mod_Tana_NumKEO
   USE mod_Tana_write_mctdh

   IMPLICIT NONE

   PRIVATE
   PUBLIC :: compute_analytical_KEO

   CONTAINS

   SUBROUTINE compute_analytical_KEO(TWOxKEO,mole, para_Tnum, Qact)
      USE mod_Tana_OpEl , ONLY : opel, assignment(=)
      USE mod_Tana_op,    ONLY : add_Vextr_new, Get_F2_F1_FROM_TWOxKEO
      IMPLICIT NONE

      TYPE(sum_opnd),        intent(inout)        :: TWOxKEO
      TYPE (zmatrix),        intent(inout)        :: mole
      TYPE (Tnum),           intent(inout)        :: para_Tnum
      real (kind=Rkind),     intent(inout)        :: Qact(:)

      type(Type_PiEulerRot), pointer    :: P_euler(:)
      TYPE(sum_opnd), pointer           :: M_mass_out(:,:)

      integer, pointer                  :: list_Qactiv(:)
      integer, pointer                  :: list_QpolytoQact(:)
      integer, pointer                  :: list_QactTOQpoly(:)
      real(kind=Rkind), pointer         :: tab_Q(:)
      character (len=Name_len), pointer :: tab_Qname(:)
      TYPE(opel)                        :: tabQpoly_Qel(mole%nb_var),tabQact_Qel(mole%nb_var)

      integer                           :: ndim
      logical                           :: constraint
      logical, pointer                  :: scalar_PiPj(:,:)
      real (kind=Rkind)                 :: Qdyn(mole%nb_var)
      logical :: new = .FALSE.
      logical :: With_Li

!     - working parameters ------------------------------------------
      integer :: iQpoly,iQprim,iQact,Qpoly_type
      integer :: i,n, j,k, i_transfo,nio
      integer :: nb_act, i_var
      logical :: frame,poly
      integer :: nb_terms_KEO_withoutVep,nb_terms_KEO_withVep

     !logical, parameter :: debug = .TRUE.
     logical, parameter :: debug = .FALSE.
     character (len=*), parameter  :: routine_name='compute_analytical_KEO'


!===========================================================
!===========================================================

      IF (debug) THEN
        write(out_unitp,*) '================================================='
        write(out_unitp,*) ' BEGINNING Tana'
        CALL flush_perso(out_unitp)
      END IF

      nullify(M_mass_out)
      poly = .false.
      i_transfo = -1
      do i = 1, size(mole%tab_Qtransfo)
        if(mole%tab_Qtransfo(i)%name_transfo .eq. 'poly') then
          poly = .true.
          i_transfo = i
          exit
        end if
      end do
      if (.not. poly .and. para_Tnum%Tana) then
        CALL Write_mole(mole,.TRUE.)
        write(out_unitp,*) ' ERROR in ',routine_name
        write(out_unitp,*) "Tana works only with the polyspherical coordinates"
        write(out_unitp,*) " Check your data input"
        STOP
      end if
      frame =  mole%tab_Qtransfo(i_transfo)%BFTransfo%frame
      if (.not.frame) then
        CALL Write_mole(mole,.TRUE.)
        write(out_unitp,*) ' ERROR in ',routine_name
        write(out_unitp,*) "The first vector should define a frame"
        write(out_unitp,*) " its corresponding data structure frame should be true"
        STOP
      end if


      CALL Qact_TO_Qdyn_FROM_ActiveTransfo(Qact,Qdyn,mole%ActiveTransfo)

      nullify(tab_Q)
      CALL alloc_array(tab_Q,shape(Qdyn),'tab_Q',routine_name)
      nullify(list_Qactiv)
      CALL alloc_array(list_Qactiv,shape(Qdyn),'list_Qactiv',routine_name)
      nullify(list_QpolytoQact)
      CALL alloc_array(list_QpolytoQact,shape(Qdyn),'list_QpolytoQact',routine_name)
      nullify(list_QactTOQpoly)
      CALL alloc_array(list_QactTOQpoly,shape(Qdyn),'list_QactTOQpoly',routine_name)
      nullify(tab_Qname)
      CALL alloc_array(tab_Qname,shape(Qdyn),Name_len,'tab_Qname',routine_name)


      ndim = size(mole%tab_Qtransfo(1)%BunchTransfo%M_Tana(:,1))
      nullify(P_euler)
      CALL alloc_array(P_euler,(/ndim/),"P_euler",routine_name)

      nullify(scalar_PiPj)
      CALL alloc_array(scalar_PiPj,(/ndim,ndim/),'scalar_PiPj',routine_name)
      scalar_PiPj(:,:) = .false.


      !-----------------------------------------------------------------
      !      DML 9/12/2012
      !      modification to use a different ordering than the polysperical one
      !      Lists have been used in "BFTransfo":list_Qpoly_TO_Qprim and list_Qprim_TO_Qpoly
      !      Now:
      !      - the order of "Qdyn" is "Qprim"
      !      - the order of "tab_Q" and "list_Qactiv" is "Qpoly"
      !      - "list_QpolytoQact" is in fact "list_QpolytoQact"
      !      - "list_QactTOQpoly" is in fact "list_QactTOQpoly"
      !-----------------------------------------------------------------
      DO iQpoly = 1, size(tab_Q)
        iQprim = mole%tab_Qtransfo(i_transfo)%BFTransfo%list_Qpoly_TO_Qprim(iQpoly)

        tab_Q(iQpoly) = Qdyn(iQprim)

        list_Qactiv(iQpoly)     = mole%ActiveTransfo%list_act_of_Qdyn(iQprim)

        list_QpolytoQact(iQpoly) = mole%ActiveTransfo%list_QdynTOQact(iQprim)

        list_QactTOQpoly(list_QpolytoQact(iQpoly) ) = iQpoly
      END DO

      nb_act     = count(mole%ActiveTransfo%list_act_of_Qdyn == 1 .OR.  &
                         mole%ActiveTransfo%list_act_of_Qdyn == 21)
      constraint = (count(mole%ActiveTransfo%list_act_of_Qdyn /= 1 .AND.&
                          mole%ActiveTransfo%list_act_of_Qdyn /= 21) > 0 )
      IF (nb_act /= mole%nb_act) THEN
        write(out_unitp,*) ' ERROR in ',routine_name
        write(out_unitp,*) "  mole%nb_act from mole is not equal to nb_act"
        write(out_unitp,*) '  nb_act,mole%nb_act',nb_act,mole%nb_act
        STOP
      END IF
      IF (mole%nb_act == 0) THEN
        write(out_unitp,*) ' ERROR in ',routine_name
        write(out_unitp,*) "  there is no active coordinates"
        write(out_unitp,*) '  mole%nb_act',mole%nb_act
        STOP
      END IF
      !-----------------------------------------------------------------

      write(out_unitp,*) '================================================='
      write(out_unitp,*) 'vectors indices in their subsystem to the indice in th BF'
      CALL flush_perso(out_unitp)
      i_var = 1
      call  iv_system_to_iv_BF(mole%tab_Qtransfo(i_transfo)%BFTransfo,i_var)
      call init_tab_num_frame_Peuler(mole%tab_Qtransfo(i_transfo)%BFTransfo,P_euler)

      write(out_unitp,*) '================================================='
      write(out_unitp,*) 'Initialization of the coordinates associated with each subsystem'
      CALL flush_perso(out_unitp)
      i_var = 1
      With_Li = .FALSE.
      call extract_qval_F_system(mole%tab_Qtransfo(i_transfo)%BFTransfo, &
      &                     tab_Q, list_Qactiv, tab_Qname, tabQpoly_Qel, i_var , With_Li)

      DO i_var=1,mole%nb_var
        write(out_unitp,*) tab_Qname(i_var)
        CALL write_op(tabQpoly_Qel(i_var))
      END DO

      write(out_unitp,*) '================================================='
      write(out_unitp,*) 'Initialization of the M_mass matrix of each subsystem'
      CALL flush_perso(out_unitp)

      do k = 1, size(mole%tab_Qtransfo(1)%BunchTransfo%M_Tana(:,1))
      do j = 1, size(mole%tab_Qtransfo(1)%BunchTransfo%M_Tana(:,1))
        if (abs(mole%tab_Qtransfo(1)%BunchTransfo%M_Tana(k,j)) < 1.0e-13_Rkind) &
              mole%tab_Qtransfo(1)%BunchTransfo%M_Tana(k,j) = zero
      end do
      end do


      call extract_bloc_matrix(mole%tab_Qtransfo(i_transfo)%BFTransfo, &
      &               mole%tab_Qtransfo(1)%BunchTransfo%M_Tana)
      call transform_M_mass_to_M_mass_opnd(mole%tab_Qtransfo(1)%BunchTransfo%M_Tana, &
      & M_mass_out)

      TWOxKEO = CZERO
      write(out_unitp,*) '================================================='
      write(out_unitp,*) ' Computation of the 2xKEO (in full dimension)'
      CALL flush_perso(out_unitp)

      call  get_opKEO(mole%tab_Qtransfo(i_transfo)%BFTransfo,  TWOxKEO, &
                      P_Euler, M_mass_out, scalar_PiPj)

      nb_terms_KEO_withoutVep = size(TWOxKEO%sum_prod_op1d)

      write(out_unitp,*) '================================================='

      IF (debug) THEN
        write(out_unitp,*) ' Write 2xKEO (in full dimension)'
        call write_op(TWOxKEO,header=.TRUE.)
      END IF

      IF ( para_Tnum%nrho == 1 .OR. para_Tnum%nrho == 2 ) THEN

        ! this version with Qact, when get_KEO_for_Qactiv is BEFORE add_Vextr_new
        !CALL add_Vextr_new(mole%tab_Qtransfo(i_transfo)%BFTransfo,TWOxKEO,&
        !                   tabQact_Qel,para_Tnum%nrho,mole%nb_act)

        ! this version with Qpoly, when get_KEO_for_Qactiv is AFTER add_Vextr_new
        CALL add_Vextr_new(mole%tab_Qtransfo(i_transfo)%BFTransfo,TWOxKEO,&
                           tabQpoly_Qel,para_Tnum%nrho,mole%nb_var)
      END IF
      nb_terms_KEO_withVep = size(TWOxKEO%sum_prod_op1d)

      write(out_unitp,*) '================================================='
      write(out_unitp,*) ' Add the extra potential term  KEO'
      write(out_unitp,*) '================================================='
      write(out_unitp,*) 'number of terms before adding Vextr=',nb_terms_KEO_withoutVep
      write(out_unitp,*) 'number of terms after  adding Vextr=',nb_terms_KEO_withVep
      CALL flush_perso(out_unitp)


      write(out_unitp,*) '================================================='
      write(out_unitp,*) ' Computation of the 2xKEO (in reduced dimension)'
      CALL flush_perso(out_unitp)
      CALL get_KEO_for_Qactiv(TWOxKEO, constraint,Qact,tabQpoly_Qel,tabQact_Qel, &
                              list_Qactiv,list_QpolytoQact)
      write(out_unitp,*) '================================================='



      IF (With_Li) THEN
        CALL write_keo_VSCFform(mole, para_Tnum%TWOxKEO, out_unitp, tab_Qname, para_Tnum%JJ)
        STOP
      ELSE
        !new = .TRUE.
        new = .FALSE.
        IF (new) THEN
          write(out_unitp,*) '================================================='
          write(out_unitp,*) ' GET F2 F1 (in reduced dimension)'
          CALL flush_perso(out_unitp)
          CALL  Get_F2_F1_FROM_TWOxKEO(mole%tab_Qtransfo(i_transfo)%BFTransfo,&
                                     TWOxKEO,para_Tnum%ExpandTWOxKEO,       &
                                     tabQact_Qel,mole%nb_act,mole%nb_var,para_Tnum%nrho)
          IF (debug) CALL write_op(para_Tnum%ExpandTWOxKEO,header=.TRUE.)
          write(out_unitp,*) '================================================='
        ELSE
          !!!! The expansion MUST be done after removing inactive terms.
          !!!! Otherwise, the KEO can be not hermitian !
          write(out_unitp,*) '================================================='
          write(out_unitp,*) ' Expand 2xKEO (in reduced dimension)'
          CALL flush_perso(out_unitp)
          CALL Expand_Sum_OpnD_TO_Sum_OpnD(TWOxKEO,para_Tnum%ExpandTWOxKEO)
          write(out_unitp,*) 'number of terms after the expansion=',size(para_Tnum%ExpandTWOxKEO%sum_prod_op1d)
          IF (debug) CALL write_op(para_Tnum%ExpandTWOxKEO,header=.TRUE.)
          write(out_unitp,*) '================================================='
        END IF
      END IF


      write(out_unitp,*) '================================================='
      write(out_unitp,*) ' output of the analytical  KEO'
      write(out_unitp,*) ' MCTDH    Form: ',para_Tnum%MCTDHForm
      write(out_unitp,*) ' VSCF     Form: ',para_Tnum%VSCFForm
      write(out_unitp,*) ' MidasCpp Form: ',para_Tnum%MidasCppForm
      write(out_unitp,*) ' LaTex    Form: ',para_Tnum%LaTexForm
      write(out_unitp,*) '================================================='
      CALL flush_perso(out_unitp)

      tab_Qname(:) = tab_Qname(list_QactTOQpoly(:)) ! to change the order du to the "constraints"

      IF (para_Tnum%MCTDHForm) THEN
        write(out_unitp,*) '================================================='
        write(out_unitp,*) "output MCTDH format"
        write(out_unitp,*) '-------------------------------------------------'
        !call write_keo_mctdh_form(mole, para_Tnum%ExpandTWOxKEO, out_unitp, tab_Qname, para_Tnum%JJ)
        call write_keo_mctdh_form(mole, para_Tnum%TWOxKEO, out_unitp, tab_Qname, para_Tnum%JJ)

        write(out_unitp,*) '================================================='
      END IF

      IF (para_Tnum%VSCFForm) THEN
        write(out_unitp,*) '================================================='
        write(out_unitp,*) 'VSCF form'
        write(out_unitp,*) '-------------------------------------------------'
        CALL write_keo_VSCFform(mole, para_Tnum%TWOxKEO, out_unitp, tab_Qname, para_Tnum%JJ)
        !CALL write_keo_VSCFform(mole, para_Tnum%ExpandTWOxKEO, out_unitp, tab_Qname, para_Tnum%JJ)
        write(out_unitp,*) '================================================='
      END IF

      IF (para_Tnum%MidasCppForm) THEN
        write(out_unitp,*) '================================================='
        write(out_unitp,*) 'MidasCpp formatted file'
        write(out_unitp,*) '-------------------------------------------------'
        CALL file_open2(name_file = 'MidasCpp_KEO.mop', iunit = nio)
        !CALL write_keo_MidasCppForm(mole, TWOxKEO, nio, tab_Qname, 0)
        CALL write_keo_MidasCppForm(mole, para_Tnum%ExpandTWOxKEO, nio, tab_Qname, 0)
        close(nio)
        write(out_unitp,*) '================================================='
        write(out_unitp,*) 'MidasCpp formatted KEO: T = '
        write(out_unitp,*) '-------------------------------------------------'
        !CALL write_keo_MidasCppForm(mole,TWOxKEO, out_unitp, tab_Qname, 0)
        CALL write_keo_MidasCppForm(mole, para_Tnum%ExpandTWOxKEO, out_unitp, tab_Qname, 0)
        write(out_unitp,*) '================================================='
      END IF

      IF (para_Tnum%LaTexForm) THEN
        write(out_unitp,*) '================================================='
        write(out_unitp,*) 'LaTex file'
        write(out_unitp,*) '-------------------------------------------------'
        CALL file_open2(name_file = 'Eq_KEO.tex', iunit = nio)
        CALL write_keo_Latexform(mole, para_Tnum%ExpandTWOxKEO, nio, tab_Qname, para_Tnum%JJ)
        close(nio)
      END IF

      ! deallocation ...

      CALL dealloc_array(P_euler,"P_euler",routine_name)

      CALL dealloc_array(M_mass_out,'M_mass_out',routine_name)

      CALL dealloc_array(tab_Q,'tab_Q',routine_name)
      CALL dealloc_array(list_Qactiv,'list_Qactiv',routine_name)

      CALL dealloc_array(list_QpolytoQact,'list_QpolytoQact',routine_name)
      CALL dealloc_array(list_QactTOQpoly,'list_QactTOQpoly',routine_name)
      CALL dealloc_array(tab_Qname,'tab_Qname',routine_name)

      CALL dealloc_array(scalar_PiPj,'scalar_PiPj',routine_name)

      IF (debug) THEN
        write(out_unitp,*) ' END Tana'
        write(out_unitp,*) '================================================='
        CALL flush_perso(out_unitp)
      END IF


   END SUBROUTINE compute_analytical_KEO
   SUBROUTINE compute_analytical_KEO_old(TWOxKEO,mole, para_Tnum, Qact)
      USE mod_Tana_OpEl , ONLY : opel, assignment(=)
      USE mod_Tana_op,    ONLY : add_Vextr_new, Get_F2_F1_FROM_TWOxKEO
      IMPLICIT NONE

      TYPE(sum_opnd),        intent(inout)        :: TWOxKEO
      TYPE (zmatrix),        intent(inout)        :: mole
      TYPE (Tnum),           intent(inout)        :: para_Tnum
      real (kind=Rkind),     intent(inout)        :: Qact(:)

      type(Type_PiEulerRot), pointer    :: P_euler(:)
      TYPE(sum_opnd), pointer           :: M_mass_out(:,:)

      integer, pointer                  :: list_Qactiv(:)
      integer, pointer                  :: list_QpolytoQact(:)
      integer, pointer                  :: list_QactTOQpoly(:)
      real(kind=Rkind), pointer         :: tab_Q(:)
      character (len=Name_len), pointer :: tab_Qname(:)
      TYPE(opel)                        :: tabQpoly_Qel(mole%nb_var),tabQact_Qel(mole%nb_var)

      integer                           :: ndim
      logical                           :: constraint
      logical, pointer                  :: scalar_PiPj(:,:)
      real (kind=Rkind)                 :: Qdyn(mole%nb_var)
      logical :: new = .FALSE.

!     - working parameters ------------------------------------------
      integer :: iQpoly,iQprim,iQact,Qpoly_type
      integer :: i,n, j,k, i_transfo,nio
      integer :: nb_act, i_var
      logical :: frame,poly,With_Li
      integer :: nb_terms_KEO_withoutVep,nb_terms_KEO_withVep

     !logical, parameter :: debug = .TRUE.
     logical, parameter :: debug = .FALSE.
     character (len=*), parameter  :: routine_name='compute_analytical_KEO_old'


!===========================================================
!===========================================================
      nullify(M_mass_out)

      write(out_unitp,*) '================================================='
      write(out_unitp,*) ' BEGINNING Tana'
      CALL flush_perso(out_unitp)
      poly = .false.
      i_transfo = -1
      do i = 1, size(mole%tab_Qtransfo)
        if(mole%tab_Qtransfo(i)%name_transfo .eq. 'poly') then
          poly = .true.
          i_transfo = i
          exit
        end if
      end do
      if (.not. poly .and. para_Tnum%Tana) then
        CALL Write_mole(mole,.TRUE.)
        write(out_unitp,*) ' ERROR in ',routine_name
        write(out_unitp,*) "Tana works only with the polyspherical coordinates"
        write(out_unitp,*) " Check your data input"
        STOP
      end if
      frame =  mole%tab_Qtransfo(i_transfo)%BFTransfo%frame
      if (.not.frame) then
        CALL Write_mole(mole,.TRUE.)
        write(out_unitp,*) ' ERROR in ',routine_name
        write(out_unitp,*) "The first vector should define a frame"
        write(out_unitp,*) " its corresponding data structure frame should be true"
        STOP
      end if


      CALL Qact_TO_Qdyn_FROM_ActiveTransfo(Qact,Qdyn,mole%ActiveTransfo)

      nullify(tab_Q)
      CALL alloc_array(tab_Q,shape(Qdyn),'tab_Q',routine_name)
      nullify(list_Qactiv)
      CALL alloc_array(list_Qactiv,shape(Qdyn),'list_Qactiv',routine_name)
      nullify(list_QpolytoQact)
      CALL alloc_array(list_QpolytoQact,shape(Qdyn),'list_QpolytoQact',routine_name)
      nullify(list_QactTOQpoly)
      CALL alloc_array(list_QactTOQpoly,shape(Qdyn),'list_QactTOQpoly',routine_name)
      nullify(tab_Qname)
      CALL alloc_array(tab_Qname,shape(Qdyn),Name_len,'tab_Qname',routine_name)


      ndim = size(mole%tab_Qtransfo(1)%BunchTransfo%M_Tana(:,1))
      nullify(P_euler)
      CALL alloc_array(P_euler,(/ndim/),"P_euler",routine_name)

      nullify(scalar_PiPj)
      CALL alloc_array(scalar_PiPj,(/ndim,ndim/),'scalar_PiPj',routine_name)
      scalar_PiPj(:,:) = .false.


      !-----------------------------------------------------------------
      !      DML 9/12/2012
      !      modification to use a different ordering than the polysperical one
      !      Lists have been used in "BFTransfo":list_Qpoly_TO_Qprim and list_Qprim_TO_Qpoly
      !      Now:
      !      - the order of "Qdyn" is "Qprim"
      !      - the order of "tab_Q" and "list_Qactiv" is "Qpoly"
      !      - "list_QpolytoQact" is in fact "list_QpolytoQact"
      !      - "list_QactTOQpoly" is in fact "list_QactTOQpoly"
      !-----------------------------------------------------------------
      DO iQpoly = 1, size(tab_Q)
        iQprim = mole%tab_Qtransfo(i_transfo)%BFTransfo%list_Qpoly_TO_Qprim(iQpoly)

        tab_Q(iQpoly) = Qdyn(iQprim)

        list_Qactiv(iQpoly)     = mole%ActiveTransfo%list_act_of_Qdyn(iQprim)

        list_QpolytoQact(iQpoly) = mole%ActiveTransfo%list_QdynTOQact(iQprim)

        list_QactTOQpoly(list_QpolytoQact(iQpoly) ) = iQpoly
      END DO

      nb_act     = count(mole%ActiveTransfo%list_act_of_Qdyn == 1 .OR.  &
                         mole%ActiveTransfo%list_act_of_Qdyn == 21)
      constraint = (count(mole%ActiveTransfo%list_act_of_Qdyn /= 1 .AND.&
                          mole%ActiveTransfo%list_act_of_Qdyn /= 21) > 0 )
      IF (nb_act /= mole%nb_act) THEN
        write(out_unitp,*) ' ERROR in ',routine_name
        write(out_unitp,*) "  mole%nb_act from mole is not equal to nb_act"
        write(out_unitp,*) '  nb_act,mole%nb_act',nb_act,mole%nb_act
        STOP
      END IF
      IF (mole%nb_act == 0) THEN
        write(out_unitp,*) ' ERROR in ',routine_name
        write(out_unitp,*) "  there is no active coordinates"
        write(out_unitp,*) '  mole%nb_act',mole%nb_act
        STOP
      END IF
      !-----------------------------------------------------------------

      write(out_unitp,*) '================================================='
      write(out_unitp,*) 'vectors indices in their subsystem to the indice in th BF'
      CALL flush_perso(out_unitp)
      i_var = 1
      call  iv_system_to_iv_BF(mole%tab_Qtransfo(i_transfo)%BFTransfo,i_var)
      call init_tab_num_frame_Peuler(mole%tab_Qtransfo(i_transfo)%BFTransfo,P_euler)

      write(out_unitp,*) '================================================='
      write(out_unitp,*) 'Initialization of the coordinates associated with each subsystem'
      CALL flush_perso(out_unitp)
      i_var = 1
      With_Li = .FALSE.
      call extract_qval_F_system(mole%tab_Qtransfo(i_transfo)%BFTransfo, &
      &                     tab_Q, list_Qactiv, tab_Qname, tabQpoly_Qel, i_var, With_Li)

      DO i_var=1,mole%nb_var
        write(out_unitp,*) tab_Qname(i_var)
        CALL write_op(tabQpoly_Qel(i_var))
      END DO

      write(out_unitp,*) '================================================='
      write(out_unitp,*) 'Initialization of the M_mass matrix of each subsystem'
      CALL flush_perso(out_unitp)

      do k = 1, size(mole%tab_Qtransfo(1)%BunchTransfo%M_Tana(:,1))
      do j = 1, size(mole%tab_Qtransfo(1)%BunchTransfo%M_Tana(:,1))
        if (abs(mole%tab_Qtransfo(1)%BunchTransfo%M_Tana(k,j)) < 1.0e-13_Rkind) &
              mole%tab_Qtransfo(1)%BunchTransfo%M_Tana(k,j) = zero
      end do
      end do


      call extract_bloc_matrix(mole%tab_Qtransfo(i_transfo)%BFTransfo, &
      &               mole%tab_Qtransfo(1)%BunchTransfo%M_Tana)
      call transform_M_mass_to_M_mass_opnd(mole%tab_Qtransfo(1)%BunchTransfo%M_Tana, &
      & M_mass_out)

      TWOxKEO = CZERO
      write(out_unitp,*) '================================================='
      write(out_unitp,*) ' Computation of the 2xKEO (in full dimension)'
      CALL flush_perso(out_unitp)

      call  get_opKEO(mole%tab_Qtransfo(i_transfo)%BFTransfo,  TWOxKEO, &
                      P_Euler, M_mass_out, scalar_PiPj)

      nb_terms_KEO_withoutVep = size(TWOxKEO%sum_prod_op1d)
      IF ( para_Tnum%nrho == 1 .OR. para_Tnum%nrho == 2 ) THEN
        CALL add_Vextr_new(mole%tab_Qtransfo(i_transfo)%BFTransfo,TWOxKEO,&
                           tabQact_Qel,para_Tnum%nrho,mole%nb_act)
      END IF
      nb_terms_KEO_withVep = size(TWOxKEO%sum_prod_op1d)

      write(out_unitp,*) '================================================='
      write(out_unitp,*) ' Add the extra potential term  KEO'
      write(out_unitp,*) '================================================='
      write(out_unitp,*) 'number of terms before adding Vextr=',nb_terms_KEO_withoutVep
      write(out_unitp,*) 'number of terms after  adding Vextr=',nb_terms_KEO_withVep
      CALL flush_perso(out_unitp)

      IF (debug) THEN
        write(out_unitp,*) ' Write 2xKEO (in full dimension)'
        call write_op(TWOxKEO,header=.TRUE.)
      END IF

      write(out_unitp,*) '================================================='
      write(out_unitp,*) ' Computation of the 2xKEO (in reduced dimension)'
      CALL flush_perso(out_unitp)
      CALL get_KEO_for_Qactiv(TWOxKEO, constraint,Qact,tabQpoly_Qel,tabQact_Qel, &
                              list_Qactiv,list_QpolytoQact)
      write(out_unitp,*) '================================================='

      !new = .TRUE.
      new = .FALSE.
      IF (new) THEN
      write(out_unitp,*) '================================================='
      write(out_unitp,*) ' GET F2 F1 (in reduced dimension)'
      CALL flush_perso(out_unitp)
      CALL  Get_F2_F1_FROM_TWOxKEO(mole%tab_Qtransfo(i_transfo)%BFTransfo,&
                                   TWOxKEO,para_Tnum%ExpandTWOxKEO,       &
                                   tabQact_Qel,mole%nb_act,mole%nb_var,para_Tnum%nrho)
      IF (debug) CALL write_op(para_Tnum%ExpandTWOxKEO,header=.TRUE.)
      write(out_unitp,*) '================================================='
      ELSE
      !!!! The expansion MUST be done after removing inactive terms.
      !!!! Otherwise, the KEO can be not hermitian !
      write(out_unitp,*) '================================================='
      write(out_unitp,*) ' Expand 2xKEO (in reduced dimension)'
      CALL flush_perso(out_unitp)
      CALL Expand_Sum_OpnD_TO_Sum_OpnD(TWOxKEO,para_Tnum%ExpandTWOxKEO)
      write(out_unitp,*) 'number of terms after the expansion=',size(para_Tnum%ExpandTWOxKEO%sum_prod_op1d)
      IF (debug) CALL write_op(para_Tnum%ExpandTWOxKEO,header=.TRUE.)
      write(out_unitp,*) '================================================='
      END IF




      write(out_unitp,*) '================================================='
      write(out_unitp,*) ' output of the analytical  KEO'
      write(out_unitp,*) ' MCTDH    Form: ',para_Tnum%MCTDHForm
      write(out_unitp,*) ' VSCF     Form: ',para_Tnum%VSCFForm
      write(out_unitp,*) ' MidasCpp Form: ',para_Tnum%MidasCppForm
      write(out_unitp,*) ' LaTex    Form: ',para_Tnum%LaTexForm
      write(out_unitp,*) '================================================='
      CALL flush_perso(out_unitp)

      tab_Qname(:) = tab_Qname(list_QactTOQpoly(:)) ! to change the order du to the "constraints"

      IF (para_Tnum%MCTDHForm) THEN
        write(out_unitp,*) '================================================='
        write(out_unitp,*) "output MCTDH format"
        write(out_unitp,*) '-------------------------------------------------'
        !call write_keo_mctdh_form(mole, para_Tnum%ExpandTWOxKEO, out_unitp, tab_Qname, para_Tnum%JJ)
        call write_keo_mctdh_form(mole, para_Tnum%TWOxKEO, out_unitp, tab_Qname, para_Tnum%JJ)

        write(out_unitp,*) '================================================='
      END IF

      IF (para_Tnum%VSCFForm) THEN
        write(out_unitp,*) '================================================='
        write(out_unitp,*) 'VSCF form'
        write(out_unitp,*) '-------------------------------------------------'
        CALL write_keo_VSCFform(mole, para_Tnum%ExpandTWOxKEO, out_unitp, tab_Qname, para_Tnum%JJ)
        write(out_unitp,*) '================================================='
      END IF

      IF (para_Tnum%MidasCppForm) THEN
        write(out_unitp,*) '================================================='
        write(out_unitp,*) 'MidasCpp formatted file'
        write(out_unitp,*) '-------------------------------------------------'
        CALL file_open2(name_file = 'MidasCpp_KEO.mop', iunit = nio)
        !CALL write_keo_MidasCppForm(mole, TWOxKEO, nio, tab_Qname, 0)
        CALL write_keo_MidasCppForm(mole, para_Tnum%ExpandTWOxKEO, nio, tab_Qname, 0)
        close(nio)
        write(out_unitp,*) '================================================='
        write(out_unitp,*) 'MidasCpp formatted KEO: T = '
        write(out_unitp,*) '-------------------------------------------------'
        !CALL write_keo_MidasCppForm(mole,TWOxKEO, out_unitp, tab_Qname, 0)
        CALL write_keo_MidasCppForm(mole, para_Tnum%ExpandTWOxKEO, out_unitp, tab_Qname, 0)
        write(out_unitp,*) '================================================='
      END IF

      IF (para_Tnum%LaTexForm) THEN
        write(out_unitp,*) '================================================='
        write(out_unitp,*) 'LaTex file'
        write(out_unitp,*) '-------------------------------------------------'
        CALL file_open2(name_file = 'Eq_KEO.tex', iunit = nio)
        CALL write_keo_Latexform(mole, para_Tnum%ExpandTWOxKEO, nio, tab_Qname, para_Tnum%JJ)
        close(nio)
      END IF

      ! deallocation ...

      CALL dealloc_array(P_euler,"P_euler",routine_name)

      CALL dealloc_array(M_mass_out,'M_mass_out',routine_name)

      CALL dealloc_array(tab_Q,'tab_Q',routine_name)
      CALL dealloc_array(list_Qactiv,'list_Qactiv',routine_name)

      CALL dealloc_array(list_QpolytoQact,'list_QpolytoQact',routine_name)
      CALL dealloc_array(list_QactTOQpoly,'list_QactTOQpoly',routine_name)
      CALL dealloc_array(tab_Qname,'tab_Qname',routine_name)

      CALL dealloc_array(scalar_PiPj,'scalar_PiPj',routine_name)

   END SUBROUTINE compute_analytical_KEO_old

   !! @description: Relation between the indices of the variable in
   !!               the BF frame and their indice in the intermediate frames
   !! @param:       F_system    The  variable in which all needed
   !!                          information will be saved (type: system)
   !!@param:        i_var input index of the first variable of the subsystem
   recursive subroutine iv_system_to_iv_BF(F_system, i_var)
     USE mod_BunchPolyTransfo, only : Type_BFTransfo

     type(Type_BFTransfo),         intent(inout)      :: F_system
     integer,                      intent(inout)      :: i_var
  
     integer                            :: i, j, iv, ivF
     integer                            :: nvec, nsub_syst 
     character (len=*), parameter       :: routine_name='iv_system_to_iv_BF'
  
       nsub_syst = 0
       do i=1, F_system%nb_vect
         if(F_system%tab_BFTransfo(i)%frame) nsub_syst = nsub_syst+1
       end do
       nvec = F_system%nb_vect-nsub_syst+1     
       CALL alloc_array(F_system%listVFr,(/nvec/),'F_system%listVFr',routine_name)
       ivF = 1         
       F_system%listVFr(ivF) = i_var
       i_var = i_var+1
       ivF = ivF+1
       do iv = 1, F_system%nb_vect
         if(F_system%tab_BFTransfo(iv)%frame) then
           call iv_system_to_iv_BF(F_system%tab_BFTransfo(iv), i_var)
         else  
           F_system%listVFr(ivF) = i_var
           i_var = i_var+1
           ivF = ivF+1
         end if    
       end do
   end subroutine iv_system_to_iv_BF



   !! @description: Initialize the tab_num_frame of the P_euler
   !!               data structure
   !! @param:       F_system    The  variable in which all needed
   !!                          information will be saved (type: BFtranfo)
   !!@param:        P_euler
   recursive subroutine init_tab_num_frame_Peuler(F_system, P_euler)
     USE mod_BunchPolyTransfo, only : Type_BFTransfo

     type(Type_BFTransfo),         intent(inout)      :: F_system
     type(Type_PiEulerRot), intent(inout)             :: P_euler(:)
  
     integer                            :: i, j, iv 
     integer                            :: nvec, nsub_syst 
     integer                            :: n_size
     character (len=*), parameter       :: routine_name='init_tab_num_frame_Peuler'

       nsub_syst = 0
       do i=1, F_system%nb_vect
         if(F_system%tab_BFTransfo(i)%frame) nsub_syst = nsub_syst+1
       end do
       nvec = F_system%nb_vect-nsub_syst+1     
       n_size = size(F_system%tab_num_Frame)
       do i = 1, nvec
         CALL alloc_array(P_Euler(F_system%listVFr(i))%Tab_num_Frame,(/n_size/),&
                                 'F_system%listVFr(i))%Tab_num_Frame',routine_name)
         do j = 1, n_size
            P_Euler(F_system%listVFr(i))%Tab_num_Frame(j) = F_system%Tab_num_Frame(j)
          end do
       end do
       do iv = 1, F_system%nb_vect
         if(F_system%tab_BFTransfo(iv)%frame) then
           call init_tab_num_frame_Peuler(F_system%tab_BFTransfo(iv), &
           &     P_Euler)
         end if    
       end do
   end subroutine init_tab_num_frame_Peuler
   


   !! @description: Relation between the indices of the variable in
   !!               the BF frame and their indice in the intermediate frames
   !! @param:       F_system    The  variable in which all needed
   !!                          information will be saved (type: system)
   !!@param:        i_var input index of the first variable of the subsystem
   recursive subroutine iv_BF_to_iv_subsystem(F_system, tab_iv, i_var)
     USE mod_BunchPolyTransfo, only : Type_BFTransfo

     type(Type_BFTransfo),         intent(inout)      :: F_system
     integer,                      intent(inout)      :: tab_iv(:)
     integer,                      intent(inout)      :: i_var
  
     integer                            :: i,  iv, ivF
     integer                            :: nvec, nsub_syst 
     character (len=*), parameter       :: routine_name='iv_system_to_iv_BF'
  
      ! if(.not.allocated(tab_v)) then
      !   write(out_unitp,*) ' ERROR in',routine_name
      !   write(out_unitp,*) "you should allocate tab_v before calling this subroutine"
      ! end if   
       ivF = 1
       tab_iv(i_var) = ivF
       i_var = i_var+1
       ivF = ivF+1
       do iv = 1, F_system%nb_vect
         if(F_system%tab_BFTransfo(iv)%frame) then
           call iv_BF_to_iv_subsystem(F_system%tab_BFTransfo(iv), tab_iv, i_var)
         else  
           tab_iv(i_var) = ivF
           i_var = i_var+1
           ivF = ivF+1
         end if    
       end do
   end subroutine iv_BF_to_iv_subsystem

   !! @description: Fill the model F_system para with the data
   !!               in tab_Q
   !! @param:       F_system    The  variable in which all needed
   !!                          information will be saved (type: system)
   recursive subroutine extract_qval_F_system(F_system, tab_Q, tab_Qactiv, &
                                              tab_Qname, tab_Qel, i_var, with_Li)
     USE mod_Tana_OpEl , ONLY : opel, assignment(=)
     USE mod_BunchPolyTransfo, only : Type_BFTransfo

     type(Type_BFTransfo),            intent(inout)      :: F_system
     real(kind = Rkind),              intent(in)         :: tab_Q(:)
     integer,                         intent(in)         :: tab_Qactiv(:)
     character(len=*),                intent(inout)      :: tab_Qname(:)
     integer,                         intent(inout)      :: i_var
     TYPE(opel),                      intent(inout)      :: tab_Qel(:)
     logical,                         intent(inout)      :: with_Li

     integer                           :: i, j, k, i_syst, iv
     integer                           :: ivF, nb_var
     integer                           :: nvec, nsub_syst
     integer                           :: i_len, iq
     character (len = Name_longlen)    :: ci
     character (len=*), parameter      :: routine_name='extract_qval_F_system'

      nsub_syst = 0
      do i=1, F_system%nb_vect
        if(F_system%tab_BFTransfo(i)%frame) nsub_syst = nsub_syst+1
      end do
      nvec = F_system%nb_vect-nsub_syst+1     
      nb_var = max(1,3*nvec-3)
      write(out_unitp,*) '************************************'
      write(out_unitp,*) 'begin sub_system', F_system%tab_num_frame
      write(out_unitp,*) 'nvec=', nvec
      write(out_unitp,*) 'nsub_syt=', nsub_syst
      do i = 1, size(F_system%euler)
        if(F_system%euler(i)) nb_var = nb_var + 1
      end do        
     ! if(F_system%nb_vect > 0) then
     !   if(F_system%tab_BFTransfo(1)%frame)  nb_var = nb_var + 1
     ! end if
     ! if(compare_tab(F_system%euler, (/.false., .false., .false./))) then
     if( F_system%nb_vect >1) then
       if(F_system%tab_BFTransfo(1)%frame) then
         do i=2, F_system%nb_vect
           if(.not.F_system%tab_BFTransfo(i)%frame) then
             nb_var = nb_var + 1
             exit
           end if
         end do
       end if
     end if

     ivF = 1         
     i_len = len(trim(F_system%name_frame))
     tab_Qel(i_var) = F_system%Qvec(1)

     tab_Qname(i_var) = "R_1^{"//trim(F_system%name_frame(2:i_len))//'}'
     write(out_unitp,*) 'i_var_F= R', ivF, 'i_var_BF=', i_var,'qval=',tab_Q(i_var)
     i_var = i_var+1
      ivF = ivF+1
      iq = 1
      do iv = 1, F_system%nb_vect
        if(F_system%tab_BFTransfo(iv)%frame) then
          call extract_qval_F_system(F_system%tab_BFTransfo(iv), tab_Q, tab_Qactiv, &
          & tab_Qname, tab_Qel, i_var ,with_Li)
        else 
         with_Li = with_Li .OR. F_system%tab_BFTransfo(iv)%Li
         iq = iq+1 
         call write_int_in_char(iq, ci)
          if(iv == 1) then
            tab_Qel(i_var) = F_system%tab_BFTransfo(iv)%Qvec(1)

            tab_Qname(i_var) = "R_"//trim(ci)//"^{"//trim(F_system%name_frame(2:i_len))//'}'
            write(out_unitp,*) 'i_var_F= R', ivF, 'i_var_BF=', i_var,'qval=',tab_Q(i_var)

            i_var = i_var+1
            ivF = ivF+1

            tab_Qel(i_var) = F_system%tab_BFTransfo(iv)%Qvec(2)

            write(out_unitp,*) 'i_var_F th=', ivF, 'i_var_BF=', i_var,'qval=',tab_Q(i_var)
            if(F_system%cos_th) then
              tab_Qname(i_var) = &
              &  "u_"//trim(ci)//"^{"//trim(F_system%name_frame(2:i_len))//'}'
            else
              tab_Qname(i_var) = &
              &"\theta_"//trim(ci)//"^{"//trim(F_system%name_frame(2:i_len))//'}'
            end if
            i_var = i_var+1
            ivF = ivF+1
          else
            tab_Qel(i_var+0) = F_system%tab_BFTransfo(iv)%Qvec(1)
            tab_Qel(i_var+1) = F_system%tab_BFTransfo(iv)%Qvec(2)
            tab_Qel(i_var+2) = F_system%tab_BFTransfo(iv)%Qvec(3)

            write(out_unitp,*) 'i_var_F R=', ivF, 'i_var_BF=', i_var, 'qval=',tab_Q(i_var)
            if(F_system%tab_BFTransfo(iv)%cart) then
              tab_Qname(i_var) = "x_"//trim(ci)//"^{"//trim(F_system%name_frame(2:i_len))//'}'
              i_var = i_var+1
              ivF = ivF+1
              write(out_unitp,*) 'i_var_F y=', ivF, 'i_var_BF=', i_var,'qval=',tab_Q(i_var)
              tab_Qname(i_var) = "y_"//trim(ci)//"^{"//trim(F_system%name_frame(2:i_len))//'}'
              i_var = i_var+1
              ivF = ivF+1
              write(out_unitp,*) 'i_var_F z=', ivF, 'i_var_BF=', i_var,'qval=',tab_Q(i_var)
              tab_Qname(i_var) = "z_"//trim(ci)//"^{"//trim(F_system%name_frame(2:i_len))//'}'
              i_var = i_var+1
              ivF = ivF+1
            else
              tab_Qname(i_var) = "R_"//trim(ci)//"^{"//trim(F_system%name_frame(2:i_len))//'}'
              i_var = i_var+1
              ivF = ivF+1
              write(out_unitp,*) 'i_var_F th=', ivF, 'i_var_BF=', i_var,'qval=',tab_Q(i_var)
              if(F_system%cos_th) then
                tab_Qname(i_var) = "u_"//trim(ci)//"^{"//trim(F_system%name_frame(2:i_len))//'}'
              else
                tab_Qname(i_var) = "\theta_"//trim(ci)//"^{"//trim(F_system%name_frame(2:i_len))//'}'
              end if
              i_var = i_var+1
              ivF = ivF+1
              write(out_unitp,*) 'i_var_F phi=', ivF, 'i_var_BF=', i_var,'qval=',tab_Q(i_var)
              tab_Qname(i_var) = "\varphi_"//trim(ci)//"^{"//trim(F_system%name_frame(2:i_len))//'}'
              i_var = i_var+1
              ivF = ivF+1
            end if
          end if
        end if
      end do
      do i = 1, size(F_system%euler)
        if(F_system%euler(i)) then
          tab_Qel(i_var) = F_system%QEuler(i)
          write(out_unitp,*) 'i_var_F euler =', ivF, 'i_var_BF=', i_var, 'qval=',tab_Q(i_var)
          if(i==1) then
            tab_Qname(i_var) = "\alpha^{"//trim(F_system%name_frame(2:i_len))//'}'
          else if(i==2) then
            if(F_system%cos_th) then
              tab_Qname(i_var) = "u_\beta^{"//trim(F_system%name_frame(2:i_len))//'}'
            else
              tab_Qname(i_var) = "\beta^{"//trim(F_system%name_frame(2:i_len))//'}'
            end if
          else
            tab_Qname(i_var) = "\gamma^{"//trim(F_system%name_frame(2:i_len))//'}'
          end if
          i_var = i_var+1
          ivF = ivF+1
        end if
      end do
      write(out_unitp,*) 'end sub_system', F_system%tab_num_frame
      write(out_unitp,*) '************************************'
      write(out_unitp,*)
      write(out_unitp,*)
   end subroutine extract_qval_F_system


   !! @description: Extract the bloc of each subsystem from variable M_mass
   !! @param:       F_system    The  variable in which all needed
   !!                          information will be saved (type: system)
   !! @param:       M_mass      matrix of mass of the global system
   Recursive subroutine  extract_bloc_matrix(F_system, M_mass)
     USE mod_BunchPolyTransfo, only : Type_BFTransfo

     type(Type_BFTransfo),                intent(inout)      :: F_system
     real(kind = Rkind),pointer                              :: M_mass(:,:)

     integer                         :: nvec, nvec_subsyst
     integer                         :: i, j, ivF, i_var
     integer                         :: nb_var, nsub_syst
     integer                         :: iv
     integer                         :: i_BF, j_BF
     character (len=*), parameter    :: routine_name='extract_bloc_matrix'
     
     if(associated(F_system%M_mass)) then
       do iv=1, F_system%nb_vect
         IF (associated(F_system%tab_BFTransfo(iv)%M_mass)) THEN
           CALL dealloc_array(F_system%tab_BFTransfo(iv)%M_mass,        &
                             'F_system%tab_BFTransfo(iv)%M_mass',routine_name)
         END IF
       end do
       CALL dealloc_array(F_system%M_mass,'F_system%M_mass',routine_name)
     end if

      if(F_system%frame) then
        nsub_syst = 0
        do i=1, F_system%nb_vect
          if(F_system%tab_BFTransfo(i)%frame) nsub_syst = nsub_syst+1
        end do
        nvec = F_system%nb_vect-nsub_syst+1
        nb_var = max(1,3*nvec-3)
        do i = 1, size(F_system%euler)
          if(F_system%euler(i)) nb_var = nb_var + 1
        end do        

        CALL alloc_array(F_system%M_mass,(/nvec, nvec/),'F_system%M_mass',routine_name)

       do  i = 1, nvec
         i_BF = F_system%listVFr(i)
         do j = 1, nvec
           j_BF = F_system%listVFr(j)
           if(M_mass(i_BF,j_BF) /= zero) then
             F_system%M_mass(i,j) = cmplx(M_mass(i_BF,j_BF),zero,kind=Rkind)
           else
             F_system%M_mass(i,j) = CZERO
           end if
         end do  
       end do  
     end if
     do iv = 1, F_system%nb_vect
       if(F_system%tab_BFTransfo(iv)%frame) then
         call extract_bloc_matrix(F_system%tab_BFTransfo(iv), M_mass)
       end if
     end do
   end subroutine extract_bloc_matrix

   !! @description: Transforms the mass matrix in the new data structure
   !! @param:       M_mass_in       matrix of mass of the global system, input
   !!                               paramater
   !! @param:       M_mass_out      matrix of mass of the global system, input
   !!                               paramater
   subroutine  transform_M_mass_to_M_mass_opnd(M_mass_in, M_mass_out)
     real(kind = Rkind), intent(in)  :: M_mass_in(:,:)
     type(sum_opnd), pointer         :: M_mass_out(:,:)

     integer                         :: i, j
     character (len=*), parameter    :: routine_name='transform_M_mass_to_M_mass_opnd'

     if(associated(M_mass_out)) then
       do i=1,size(M_mass_out(:,1))
         do j=1,size(M_mass_out(1,:))
           call delete_op(M_mass_out(i,j))
         end do
       end do
       CALL dealloc_array(M_mass_out,'M_mass_out',routine_name)
     end if
     CALL alloc_array(M_mass_out,shape(M_mass_in),'M_mass_out',routine_name)

     do i = 1, size(M_mass_in(:,1))
     do j = 1, size(M_mass_in(1,:))
       if(abs(M_mass_in(i,j)) > ONETENTH**13) then
         M_mass_out(i,j) = M_mass_in(i,j)
       else
         M_mass_out(i,j) = CZERO
       end if
     end do
     end do  
     
   end subroutine transform_M_mass_to_M_mass_opnd

END MODULE mod_Tana_keo
