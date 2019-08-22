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

   MODULE mod_PrimOp

    !USE mod_Coord_KEO

    USE mod_nDFit
    USE mod_PrimOp_def
    USE mod_OTF_def
    USE mod_OTF
    USE mod_SimpleOp

   IMPLICIT NONE

   PUBLIC

   PRIVATE   dnOp_num_grid_v2, calc4_NM_TO_sym, calc5_NM_TO_sym

   PRIVATE   get_hess_k, Set_RPHpara_AT_Qact1_opt2, Set_RPHpara_AT_Qact1_opt01
   !PRIVATE   calc_freq, calc_freq_block, calc_freq_WITH_d0c, calc_freqNM, calc_freq_width
   !PRIVATE   H0_symmetrization, sort_with_Tab

   CONTAINS

      SUBROUTINE Sub_init_dnOp(mole,para_Tnum,para_PES)
      USE mod_system
      USE mod_SimpleOp,   only : param_d0MatOp,Init_d0MatOp,dealloc_d0MatOp
      USE mod_PrimOp_def, only : param_PES
      USE mod_Coord_KEO,  only : zmatrix,Tnum
      IMPLICIT NONE

!----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix)   :: mole
      TYPE (Tnum)      :: para_Tnum


      TYPE (param_PES) :: para_PES

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------

!----- working variables ----------------------------------------
      TYPE (param_d0MatOp), allocatable :: d0MatOp(:)
      integer                           :: k,nb_Op
      real (kind=Rkind)                 :: Qact(mole%nb_var)
      integer                           :: err_io
      character (len=Name_longlen)      :: name_dum

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Sub_init_dnOp'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'para_PES%nb_scalar_Op ',para_PES%nb_scalar_Op
         CALL flush_perso(out_unitp)
       END IF
!-----------------------------------------------------------
     allocate(d0MatOp(para_PES%nb_scalar_Op+2))

      nb_Op = size(d0MatOp)

      CALL Init_d0MatOp(d0MatOp(1),para_PES%Type_HamilOp,mole%nb_act,   &
                        para_PES%nb_elec,JRot=para_Tnum%JJ,             &
                        cplx=para_PES%pot_cplx,direct_KEO=para_PES%direct_KEO) ! H

      DO k=2,nb_Op
        CALL Init_d0MatOp(d0MatOp(k),0,mole%nb_act,para_PES%nb_elec,    &
                          JRot=para_Tnum%JJ,cplx=.FALSE.,direct_KEO=.FALSE.) ! Scalar Operator
      END DO


      Qact(:) = mole%ActiveTransfo%Qact0(:)

      IF (para_PES%QMLib) THEN
        IF (debug) write(out_unitp,*) 'Initialization with Quantum Model Lib'

#if __QML == 1
        CALL sub_Init_Qmodel(mole%nb_act,para_PES%nb_elec,'read_model',.FALSE.,0)
#else
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) ' The "Quantum Model Lib" (QML) library is not present!'
        write(out_unitp,*) '  Qmodel cannot be intialized!'
        write(out_unitp,*) 'Use another potential/model'
        STOP 'QML is not present'
#endif

        IF (allocated(para_PES%Qit_TO_QQMLib)) THEN
          CALL dealloc_NParray(para_PES%Qit_TO_QQMLib,'Qit_TO_QQMLib',name_sub)
        END IF
        CALL alloc_NParray(para_PES%Qit_TO_QQMLib,(/ mole%nb_act /),'Qit_TO_QQMLib',name_sub)
        para_PES%Qit_TO_QQMLib(:) = (/ (k,k=1,mole%nb_act) /)

        IF (para_PES%pot_itQtransfo == mole%nb_Qtransfo-1) THEN ! Qdyn Coord
          read(in_unitp,*,IOSTAT=err_io) name_dum,para_PES%Qit_TO_QQMLib
          IF (err_io /= 0) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) '  while reading "Qit_TO_QQMLib"'
            write(out_unitp,*) ' end of file or end of record'
            write(out_unitp,*) ' Probably, you have forgotten the list of integers ...'
            write(out_unitp,*) ' Check your data !!'
            STOP
          END IF
        END IF
      ELSE
        CALL get_d0MatOp_AT_Qact(Qact,d0MatOp,mole,para_Tnum,para_PES)
      END IF

      DO k=1,nb_Op
        CALL dealloc_d0MatOp(d0MatOp(k))
      END DO
      deallocate(d0MatOp)

      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF


      END SUBROUTINE Sub_init_dnOp

!================================================================
!    subroutine enables to calculate the energy, gradient and hessian
!     On-the-fly (with gaussian or gamess) or not
!
!    input : coordinates Qact (used in the dynamics) unit (bohr, radian)
!            nderiv = 0 (pot only : d0pot)
!            nderiv = 1 (pot + gradient only : d0pot, d1pot)
!            nderiv = 2 (pot + gradient + hessian only : d0pot, d1pot, d2pot)
!
!    output : dnE%d..  dnMu(:)%d...
!================================================================
      SUBROUTINE get_d0MatOp_AT_Qact(Qact,d0MatOp,mole,para_Tnum,para_PES)
      USE mod_system
      USE mod_dnSVM
      use mod_nDFit, only: sub_ndfunc_from_ndfit
      USE mod_Coord_KEO
      USE mod_SimpleOp
      USE mod_PrimOp_def
      USE mod_OTF
      IMPLICIT NONE

!----- for the zmatrix and Tnum --------------------------------------
      TYPE (Tnum)    :: para_Tnum
      TYPE (zmatrix) :: mole

!----- for Qact ... ---------------------------------------------
      real (kind=Rkind), intent(inout) :: Qact(:)
      TYPE (param_PES) :: para_PES

!----- input output variables ----------------------------------------
      integer           :: nb_Op
      TYPE (param_d0MatOp), intent(inout) :: d0MatOp(:)

      integer           :: nderivE,nderivS
      integer           :: nderivImE,nderivScal


!----- working variables -------------------------------------------------
      integer             :: nderivScal_loc


      real (kind=Rkind)   :: Qxyz(mole%ncart_act)
      TYPE(Type_dnVec)    :: dnXin
      real (kind=Rkind)   :: d0Scal_loc1(para_PES%nb_elec,para_PES%nb_elec,para_PES%nb_scalar_Op)
      real (kind=Rkind)   :: d0Scal_loc2(para_PES%nb_elec,para_PES%nb_elec,para_PES%nb_scalar_Op)
      real (kind=Rkind)   :: d0T(3,3) ! for the Eckart rotation matrix

      logical             :: Gcenter,Cart_transfo

      integer             :: i,i1,i2,ie,je,io,iOpE,itermE,iOpS,iOpScal,itermS,iOp,iterm

!     - for the conversion gCC -> gzmt=d1pot -----------
      TYPE(Type_dnS), allocatable :: MatdnECC(:,:)
      TYPE(Type_dnS), allocatable :: MatdnScalCC(:,:,:)

      real (kind=Rkind) :: mat_V(para_PES%nb_elec,para_PES%nb_elec)
      real (kind=Rkind) :: mat_imV(para_PES%nb_elec,para_PES%nb_elec)
      real (kind=Rkind) :: mat_ScalOp(para_PES%nb_elec,para_PES%nb_elec,para_PES%nb_scalar_Op)

      ! for HarD
      real (kind=Rkind) :: Vinact
      real (kind=Rkind), allocatable :: Qact1(:)
      real (kind=Rkind), allocatable :: Qinact21(:)
      real (kind=Rkind), allocatable :: Qit(:)
      real (kind=Rkind) :: Qdyn(mole%nb_var)

      integer :: iQa,iQ,iQact1,iQinact21

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='get_d0MatOp_AT_Qact'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      nb_Op = size(d0MatOp)
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Qact',Qact
        write(out_unitp,*) 'nb_Op',nb_Op
        write(out_unitp,*) 'nb_scalar_Op',para_PES%nb_scalar_Op
        write(out_unitp,*) 'calc_scalar_Op',para_PES%calc_scalar_Op
        write(out_unitp,*) 'pot_cplx',para_PES%pot_cplx
        write(out_unitp,*) 'pot_itQtransfo',para_PES%pot_itQtransfo
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

      !----------------------------------------------------------------
      IF (.NOT. para_PES%Read_OnTheFly_only) THEN
        CALL sub_QactTOQit(Qact,Qit,para_PES%pot_itQtransfo,mole,.FALSE.)
        !write(6,*) 'Qact',Qact
        !write(6,*) 'Qit',Qit
      ELSE
        ! why this allocation ????
        IF (allocated(Qit)) CALL dealloc_NParray(Qit,'Qit',name_sub)
        CALL alloc_NParray(Qit,(/mole%ncart_act/),'Qit',name_sub)
        Qit(:) = ZERO
      END IF


      !----------------------------------------------------------------

!     ----------------------------------------------------------------
      iOpE    = 1
      itermE  = d0MatOp(iOpE)%derive_term_TO_iterm(0,0)
      nderivE = 0
      iOpS    = iOpE
      iOpScal = iOpE

      !------------ The Overlap -------------------------------------
      IF (nb_Op >= 2) THEN
        iOpS    = iOpS + 1
        iOpScal = iOpS
        itermS  = d0MatOp(iOpS)%derive_term_TO_iterm(0,0)
        nderivS = 0
        DO ie=1,para_PES%nb_elec
          d0MatOp(iOpS)%ReVal(ie,ie,itermS) = ONE
        END DO
      END IF
      !----------------------------------------------------------------


      IF (nb_Op >= 3) THEN
        nderivScal = 0
        iOpScal    = iOpScal + 1
      ELSE
        nderivScal = -1
      END IF

      IF (para_PES%OnTheFly) THEN

        IF (nderivScal > -1 .AND. para_PES%nb_scalar_Op < 3) THEN
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) 'nderivScal > -1 and nb_scalar_Op < 3'
          write(out_unitp,*) 'nderivScal (mu)',nderivScal
          write(out_unitp,*) 'nb_scalar_Op',para_PES%nb_scalar_Op
          write(out_unitp,*) 'With on-the-fly calculation,'
          write(out_unitp,*) ' nb_scalar_Op MUST be >= 2 !'
          STOP
        END IF

        allocate(MatdnECC(para_PES%nb_elec,para_PES%nb_elec))
        allocate(MatdnScalCC(para_PES%nb_elec,para_PES%nb_elec,para_PES%nb_scalar_Op))

        CALL dnOp_grid_OnTheFly(Qit,MatdnECC,nderivE,                   &
                                MatdnScalCC,nderivScal,                 &
                                mole,para_PES)


        !write(77,*) Qact(1:mole%nb_act),MatdnECC%d0,MatdnScalCC%d0

        !----------------------------------------------------------------
        !- then conversion: CC=>Q
        DO ie=1,para_PES%nb_elec
        DO je=1,para_PES%nb_elec
          d0MatOp(iOpE)%ReVal(:,:,itermE) = MatdnECC(:,:)%d0
        END DO
        END DO
        CALL dealloc_MatOFdnS(MatdnECC)
        deallocate(MatdnECC)


        !- then conversion: CC=>Q
        IF (nderivScal > -1) THEN
          DO i=1,para_PES%nb_scalar_Op
            iterm = d0MatOp(iOpScal-1+i)%derive_term_TO_iterm(0,0)
            d0MatOp(iOpScal-1+i)%ReVal(:,:,iterm) = MatdnScalCC(:,:,i)%d0
          END DO
          DO i=1,para_PES%nb_scalar_Op
            CALL dealloc_MatOFdnS(MatdnScalCC(:,:,i))
          END DO
        END IF
        deallocate(MatdnScalCC)
        !----------------------------------------------------------------

        !----------------------------------------------------------------
        DO ie=1,para_PES%nb_elec
          d0MatOp(iOpE)%ReVal(ie,ie,itermE) =                           &
                       d0MatOp(iOpE)%ReVal(ie,ie,itermE) - para_PES%pot0
        END DO
        !----------------------------------------------------------------

      ELSE

          IF (para_PES%QMLib) THEN
            IF (debug) THEN
               write(out_unitp,*) 'With Quantum Model Lib'
               write(out_unitp,*) 'QQMLib',Qit(para_PES%Qit_TO_QQMLib)
            END IF
#if __QML == 1
            CALL sub_Qmodel_V(d0MatOp(iOpE)%ReVal(1,1,itermE),Qit(para_PES%Qit_TO_QQMLib))
#else
            write(out_unitp,*) 'ERROR in ',name_sub
            write(out_unitp,*) ' The "Quantum Model Lib" (QML) library is not present!'
            write(out_unitp,*) 'Use another potential/model'
            STOP 'QML is not present'
#endif

          ELSE IF (para_PES%nDfit_Op) THEN
            IF (debug) write(out_unitp,*) 'With nDFit'
            IF (para_PES%nb_elec > 1) STOP 'ERROR nb_elec > 1 with nDFit'

            ! potential
            CALL sub_nDFunc_FROM_nDFit(d0MatOp(iOpE)%ReVal(1,1,itermE), &
                                       Qit,para_PES%para_nDFit_V)

            ! Scalar Op
            IF (para_PES%calc_scalar_Op) THEN
              DO i=1,para_PES%nb_scalar_Op
                iterm = d0MatOp(iOpScal-1+i)%derive_term_TO_iterm(0,0)
                CALL sub_nDFunc_FROM_nDFit(                             &
                                 d0MatOp(iOpScal-1+i)%ReVal(1,1,iterm), &
                             Qit,para_PES%para_nDFit_Scalar_Op(i))
              END DO
            END IF
          ELSE
            IF (debug) THEN
              write(out_unitp,*) 'With calcN_op'
              CALL flush_perso(out_unitp)
            END IF

            CALL calcN_op(d0MatOp(iOpE)%ReVal(:,:,itermE),              &
                          mat_imV,mat_ScalOp,                           &
                          para_PES%nb_elec,para_PES%nb_scalar_Op,       &
                          Qit,size(Qit),                                &
                          mole,para_PES%calc_scalar_Op,para_PES%pot_cplx)

            IF (d0MatOp(iOpE)%cplx) THEN
              d0MatOp(iOpE)%ImVal(:,:) = mat_imV(:,:)
            END IF

            DO i=1,para_PES%nb_scalar_Op
              iterm = d0MatOp(iOpScal-1+i)%derive_term_TO_iterm(0,0)
              d0MatOp(iOpScal-1+i)%ReVal(:,:,iterm) = mat_ScalOp(:,:,i)
            END DO

            !----------------------------------------------------------------
            DO ie=1,para_PES%nb_elec
             d0MatOp(iOpE)%ReVal(ie,ie,itermE) =                        &
                       d0MatOp(iOpE)%ReVal(ie,ie,itermE) - para_PES%pot0
            END DO
            !----------------------------------------------------------------
          END IF
          IF (para_PES%HarD .AND. associated(mole%RPHTransfo) .AND. para_PES%nb_elec == 1) THEN
            !here it should be Qin of RPH (therefore Qdyn ?????)
            CALL Qact_TO_Qdyn_FROM_ActiveTransfo(Qact,Qdyn,mole%ActiveTransfo)
            !write(6,*) 'test HARD without HAC'

            ! transfert the dnQin coordinates: type21 in dnVecQin and ....
            !   the other (active, rigid ..) in dnQout
            iQinact21 = 0
            iQact1    = 0

            CALL alloc_NParray(Qact1,(/mole%RPHTransfo%nb_act1/),       &
                              "Qact1",name_sub)

            CALL alloc_NParray(Qinact21,(/mole%RPHTransfo%nb_inact21/), &
                              "Qinact21",name_sub)

            DO iQ=1,mole%nb_var
              IF (mole%RPHTransfo%list_act_OF_Qdyn(iQ) == 21) THEN
                iQinact21 = iQinact21 + 1
                Qinact21(iQinact21) = Qdyn(iQ)
              ELSE IF (mole%RPHTransfo%list_act_OF_Qdyn(iQ) == 1) THEN
                iQact1 = iQact1 + 1
                Qact1(iQact1) = Qdyn(iQ)
              END IF
            END DO
            !write(6,*) 'Qact1',Qact1
            !write(6,*) 'Qinact21',Qinact21

            ! find the iQa from tab_RPHpara_AT_Qact1
            DO iQa=1,mole%RPHTransfo%nb_Qa
              IF (sum(abs(Qact1-mole%RPHTransfo%tab_RPHpara_AT_Qact1(iQa)%Qact1)) < ONETENTH**5) EXIT
            END DO
            IF (iQa > mole%RPHTransfo%nb_Qa) THEN
              IF (sum(abs(Qact1-mole%RPHTransfo%RPHpara_AT_Qref(1)%Qact1)) < ONETENTH**5) THEN
                Vinact = HALF*sum(mole%RPHTransfo%RPHpara_AT_Qref(1)%dnehess%d0(:)*Qinact21(:)**2)
              ELSE
                CALL Write_RPHTransfo(mole%RPHTransfo)
                write(out_unitp,*) ' Qact1',Qact1(:)
                write(out_unitp,*) 'ERROR in ',name_sub
                write(out_unitp,*) ' I cannot find Qact1(:) in tab_RPHpara_AT_Qact1'
                write(out_unitp,*) '  or  in tab_RPHpara_AT_Qref(1)'
                STOP
              END IF
           ELSE
             Vinact = HALF*sum(mole%RPHTransfo%tab_RPHpara_AT_Qact1(iQa)%dnehess%d0(:)*Qinact21(:)**2)
           END IF

            !write(6,*) 'iQa',iQa
            !write(6,*) 'Qinact21',Qinact21(:)
            !write(6,*) 'dnehess',mole%RPHTransfo%tab_RPHpara_AT_Qact1(iQa)%dnehess%d0(:)
            !write(6,*) 'Vinact',Vinact

            DO ie=1,para_PES%nb_elec
              d0MatOp(iOpE)%ReVal(ie,ie,itermE) =                       &
                              d0MatOp(iOpE)%ReVal(ie,ie,itermE) + Vinact
            END DO

            CALL dealloc_NParray(Qact1,"Qact1",name_sub)
            CALL dealloc_NParray(Qinact21,"Qinact21",name_sub)

          END IF

      END IF

      DO ie=1,para_PES%nb_elec
        para_PES%min_pot = min(para_PES%min_pot,d0MatOp(iOpE)%ReVal(ie,ie,itermE))
        para_PES%max_pot = max(para_PES%max_pot,d0MatOp(iOpE)%ReVal(ie,ie,itermE))
      END DO

      !----------------------------------------------------------------
      IF (mole%Rot_Dip_with_EC .AND. nderivScal > -1 .AND.              &
                                         para_PES%nb_scalar_Op > 2) THEN
        CALL alloc_dnSVM(dnXin,mole%ncart,mole%nb_act,nderiv=nderivScal)

        CALL sub_QactTOdnx(Qact,dnXin,mole,nderiv=nderivScal,           &
                                    Gcenter=.TRUE.,Cart_transfo=.FALSE.)

        IF (debug) write(out_unitp,*) ' WARNING dip(:) are rotated'

        ! initial rotation of the dipole moment
        d0T(:,:) = mole%tab_Cart_transfo(1)%CartesianTransfo%Rot_initial

        DO i=1,para_PES%nb_scalar_Op
          iterm = d0MatOp(iOpScal-1+i)%derive_term_TO_iterm(0,0)
          d0Scal_loc1(:,:,i) = d0MatOp(iOpScal-1+i)%ReVal(:,:,iterm)
        END DO
        DO ie=1,para_PES%nb_elec
        DO je=1,para_PES%nb_elec
          d0Scal_loc2(ie,je,:) = matmul(d0T,d0Scal_loc1(ie,je,:))
        END DO
        END DO
        ! End initial rotation

        ! Eckart rotation of the dipole moment
        CALL calc_EckartRot(dnXin,d0T,                                  &
                 mole%tab_Cart_transfo(1)%CartesianTransfo,Qact)

        ! rotation of the dipole moment
        DO ie=1,para_PES%nb_elec
        DO je=1,para_PES%nb_elec
          d0Scal_loc1(ie,je,:) = matmul(d0T,d0Scal_loc2(ie,je,:))
        END DO
        END DO

        DO i=1,para_PES%nb_scalar_Op
          iterm = d0MatOp(iOpScal-1+i)%derive_term_TO_iterm(0,0)
          d0MatOp(iOpScal-1+i)%ReVal(:,:,iterm) = d0Scal_loc1(:,:,i)
        END DO
        ! End Eckart rotation

        CALL dealloc_dnSVM(dnXin)

      END IF
      !----------------------------------------------------------------

      IF (allocated(Qit)) CALL dealloc_NParray(Qit,'Qit',name_sub)


!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'd0MatOp(:)'
        CALL Write_Tab_OF_d0MatOp(d0MatOp)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

     END SUBROUTINE get_d0MatOp_AT_Qact
     SUBROUTINE get_dnMatOp_AT_Qact(Qact,Tab_dnMatOp,mole,para_Tnum,para_PES,nderiv)
      USE mod_system
      USE mod_dnSVM
      use mod_nDFit, only: sub_ndfunc_from_ndfit
      USE mod_Coord_KEO
      USE mod_SimpleOp
      USE mod_PrimOp_def
      USE mod_OTF
      IMPLICIT NONE

!----- for the zmatrix and Tnum --------------------------------------
      TYPE (Tnum)    :: para_Tnum
      TYPE (zmatrix) :: mole

!----- for Qact ... ---------------------------------------------
      real (kind=Rkind), intent(inout) :: Qact(:)
      TYPE (param_PES) :: para_PES

!----- input output variables ----------------------------------------
      TYPE (param_dnMatOp), intent(inout) :: Tab_dnMatOp(:)
      integer, intent(in), optional       :: nderiv


      integer           :: nderivE,nderivS,nderiv_loc
      integer           :: nderivImE,nderivScal


!----- working variables -------------------------------------------------
      integer             :: nderivScal_loc
      integer             :: nb_Op


      real (kind=Rkind)   :: Qxyz(mole%ncart_act)
      TYPE(Type_dnVec)    :: dnXin,dnXout
      TYPE(Type_dnS)      :: dnScal_loc1(para_PES%nb_scalar_Op)
      TYPE(Type_dnS)      :: dnScal_loc2(para_PES%nb_scalar_Op)
      TYPE(Type_dnS)      :: dnT(3,3) ! for the Eckart rotation matrix
      TYPE(Type_dnS)      :: dnXref(3,mole%nat_act)

      logical             :: Gcenter,Cart_transfo

      integer             :: i,i1,i2,ie,je,io,iOpE,itermE,iOpS,iOpScal,itermS,iOp,iterm

!     - for the conversion gCC -> gzmt=d1pot -----------
      TYPE(Type_dnS) :: MatdnECC(para_PES%nb_elec,para_PES%nb_elec)
      TYPE(Type_dnS) :: MatdnScalCC(para_PES%nb_elec,para_PES%nb_elec,para_PES%nb_scalar_Op)

      real (kind=Rkind) :: mat_V(para_PES%nb_elec,para_PES%nb_elec)
      real (kind=Rkind) :: mat_imV(para_PES%nb_elec,para_PES%nb_elec)
      real (kind=Rkind) :: mat_ScalOp(para_PES%nb_elec,para_PES%nb_elec,para_PES%nb_scalar_Op)

      ! for HarD
      real (kind=Rkind) :: Vinact
      real (kind=Rkind), allocatable :: Qact1(:)
      real (kind=Rkind), allocatable :: Qinact21(:)
      real (kind=Rkind), allocatable :: Qit(:)
      real (kind=Rkind) :: Qdyn(mole%nb_var)

      integer :: iQa,iQ,iQact1,iQinact21

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='get_dnMatOp_AT_Qact'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      nb_Op = size(Tab_dnMatOp)

      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Qact',Qact
        write(out_unitp,*) 'nb_Op',nb_Op
        write(out_unitp,*) 'nb_scalar_Op',para_PES%nb_scalar_Op
        write(out_unitp,*) 'calc_scalar_Op',para_PES%calc_scalar_Op
        write(out_unitp,*) 'pot_cplx',para_PES%pot_cplx
        write(out_unitp,*) 'pot_itQtransfo',para_PES%pot_itQtransfo
      END IF
!-----------------------------------------------------------

      IF (present(nderiv)) THEN
        nderiv_loc = nderiv
      ELSE
        nderiv_loc = huge(1)
      END IF

      !----------------------------------------------------------------
      IF (.NOT. para_PES%Read_OnTheFly_only) THEN
        CALL sub_QactTOQit(Qact,Qit,para_PES%pot_itQtransfo,mole,.FALSE.)
        !write(6,*) 'Qact',Qact
        !write(6,*) 'Qit',Qit
      ELSE
        ! why this allocation ????
        IF (allocated(Qit)) CALL dealloc_NParray(Qit,'Qit',name_sub)
        CALL alloc_NParray(Qit,(/mole%ncart_act/),'Qit',name_sub)
        Qit(:) = ZERO
      END IF


      !----------------------------------------------------------------

!     ----------------------------------------------------------------
      iOpE    = 1
      itermE  = Tab_dnMatOp(iOpE)%derive_term_TO_iterm(0,0)
      nderivE = min(nderiv_loc,Tab_dnMatOp(iOpE)%nderiv)
      iOpS    = iOpE
      iOpScal = iOpE
      !------------ The Overlap -------------------------------------
      IF (nb_Op >= 2) THEN
        iOpS    = iOpS + 1
        iOpScal = iOpS
        itermS  = Tab_dnMatOp(iOpS)%derive_term_TO_iterm(0,0)
        nderivS = min(nderiv_loc,Tab_dnMatOp(iOpS)%nderiv)
        DO ie=1,para_PES%nb_elec
          Tab_dnMatOp(iOpS)%tab_dnMatOp(ie,ie,itermE)%d0 = ONE
        END DO
      END IF
      !----------------------------------------------------------------


      IF (nb_Op >= 3) THEN
        iOpScal    = iOpScal + 1
        nderivScal = min(nderiv_loc,Tab_dnMatOp(iOpScal)%nderiv)
      ELSE
        nderivScal =  -1
      END IF

      IF (para_PES%Read_OnTheFly_only .OR. (para_PES%OnTheFly .AND.     &
          (nderivE == 0 .OR. .NOT. para_PES%deriv_WITH_FiniteDiff))) THEN

        IF (nderivScal > -1 .AND. para_PES%nb_scalar_Op < 3) THEN
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) 'nderivScal > -1 and nb_scalar_Op < 3'
          write(out_unitp,*) 'nderivScal (mu)',nderivScal
          write(out_unitp,*) 'nb_scalar_Op',para_PES%nb_scalar_Op
          write(out_unitp,*) 'With on-the-fly calculation,'
          write(out_unitp,*) ' nb_scalar_Op MUST be >= 2 !'
          STOP
        END IF
        CALL dnOp_grid_OnTheFly(Qit,MatdnECC,nderivE,                   &
                                MatdnScalCC,nderivScal,                 &
                                mole,para_PES)


        !write(77,*) Qact(1:mole%nb_act),MatdnECC%d0,MatdnScalCC%d0

        !----------------------------------------------------------------
        !- then conversion: CC=>Q
        DO ie=1,para_PES%nb_elec
        DO je=1,para_PES%nb_elec
         CALL sub_dnFCC_TO_dnFcurvi(Qact,MatdnECC(ie,je),             &
                       Tab_dnMatOp(iOpE)%tab_dnMatOp(ie,je,itermE),mole)
        END DO
        END DO
        CALL dealloc_MatOFdnS(MatdnECC)

        !- then conversion: CC=>Q
        IF (nderivScal > -1) THEN
          DO i=1,para_PES%nb_scalar_Op
            iterm = Tab_dnMatOp(iOpScal-1+i)%derive_term_TO_iterm(0,0)
            DO ie=1,para_PES%nb_elec
            DO je=1,para_PES%nb_elec
              CALL sub_dnFCC_TO_dnFcurvi(Qact,MatdnScalCC(ie,je,i),   &
                      Tab_dnMatOp(iOpScal-1+i)%tab_dnMatOp(ie,je,iterm),mole)

            END DO
            END DO
          END DO
          DO i=1,3
            CALL dealloc_MatOFdnS(MatdnScalCC(:,:,i))
          END DO
        END IF
        !----------------------------------------------------------------

        !----------------------------------------------------------------
        DO ie=1,para_PES%nb_elec
          Tab_dnMatOp(iOpE)%tab_dnMatOp(ie,ie,itermE)%d0 =              &
                       Tab_dnMatOp(iOpE)%tab_dnMatOp(ie,ie,itermE)%d0 - &
                                                          para_PES%pot0
        END DO
        !----------------------------------------------------------------

      ELSE

        IF (nderivE == 0 ) THEN
          IF (para_PES%nDfit_Op) THEN
            IF (debug) write(out_unitp,*) 'With nDFit'
            IF (para_PES%nb_elec > 1) STOP 'ERROR nb_elec > 1 with nDFit'

            ! potential
            CALL sub_nDFunc_FROM_nDFit(                                 &
                          Tab_dnMatOp(iOpE)%tab_dnMatOp(1,1,itermE)%d0, &
                                       Qit,para_PES%para_nDFit_V)

            ! Scalar Op
            IF (para_PES%calc_scalar_Op .AND. nb_Op >=3) THEN
              DO i=1,para_PES%nb_scalar_Op
                iterm = Tab_dnMatOp(iOpScal-1+i)%derive_term_TO_iterm(0,0)
                CALL sub_nDFunc_FROM_nDFit(                             &
                        Tab_dnMatOp(iOpScal-1+i)%tab_dnMatOp(1,1,iterm)%d0, &
                             Qit,para_PES%para_nDFit_Scalar_Op(i))
              END DO
            END IF
          ELSE
            IF (debug) write(out_unitp,*) 'With calcN_op'

            CALL calcN_op(Tab_dnMatOp(iOpE)%tab_dnMatOp(:,:,itermE)%d0, &
                          mat_imV,mat_ScalOp,                           &
                          para_PES%nb_elec,para_PES%nb_scalar_Op,       &
                          Qit,size(Qit),                     &
                          mole,para_PES%calc_scalar_Op,para_PES%pot_cplx)

            IF (Tab_dnMatOp(iOpE)%cplx) THEN
              Tab_dnMatOp(iOpE)%Im_dnMatOp(:,:)%d0 = mat_imV(:,:)
            END IF

            IF (nb_Op >=3) THEN
              DO i=1,para_PES%nb_scalar_Op
                iterm = Tab_dnMatOp(iOpScal-1+i)%derive_term_TO_iterm(0,0)
                Tab_dnMatOp(iOpScal-1+i)%tab_dnMatOp(:,:,iterm)%d0 = mat_ScalOp(:,:,i)
              END DO
            END IF

            !----------------------------------------------------------------
            DO ie=1,para_PES%nb_elec
              Tab_dnMatOp(iOpE)%tab_dnMatOp(ie,ie,itermE)%d0 =          &
                       Tab_dnMatOp(iOpE)%tab_dnMatOp(ie,ie,itermE)%d0 - &
                                           para_PES%pot0
            END DO
            !----------------------------------------------------------------
          END IF
          IF (para_PES%HarD .AND. associated(mole%RPHTransfo) .AND. para_PES%nb_elec == 1) THEN
            !here it should be Qin of RPH (therefore Qdyn ?????)
            CALL Qact_TO_Qdyn_FROM_ActiveTransfo(Qact,Qdyn,mole%ActiveTransfo)
            !write(6,*) 'test HARD without HAC'

            ! transfert the dnQin coordinates: type21 in dnVecQin and ....
            !   the other (active, rigid ..) in dnQout
            iQinact21 = 0
            iQact1    = 0

            CALL alloc_NParray(Qact1,(/mole%RPHTransfo%nb_act1/),       &
                              "Qact1",name_sub)

            CALL alloc_NParray(Qinact21,(/mole%RPHTransfo%nb_inact21/), &
                              "Qinact21",name_sub)

            DO iQ=1,mole%nb_var
              IF (mole%RPHTransfo%list_act_OF_Qdyn(iQ) == 21) THEN
                iQinact21 = iQinact21 + 1
                Qinact21(iQinact21) = Qdyn(iQ)
              ELSE IF (mole%RPHTransfo%list_act_OF_Qdyn(iQ) == 1) THEN
                iQact1 = iQact1 + 1
                Qact1(iQact1) = Qdyn(iQ)
              END IF
            END DO
            !write(6,*) 'Qact1',Qact1
            !write(6,*) 'Qinact21',Qinact21

            ! find the iQa from tab_RPHpara_AT_Qact1
            DO iQa=1,mole%RPHTransfo%nb_Qa
              IF (sum(abs(Qact1-mole%RPHTransfo%tab_RPHpara_AT_Qact1(iQa)%Qact1)) < ONETENTH**5) EXIT
            END DO
            IF (iQa > mole%RPHTransfo%nb_Qa) THEN
              IF (sum(abs(Qact1-mole%RPHTransfo%RPHpara_AT_Qref(1)%Qact1)) < ONETENTH**5) THEN
                Vinact = HALF*sum(mole%RPHTransfo%RPHpara_AT_Qref(1)%dnehess%d0(:)*Qinact21(:)**2)
              ELSE
                CALL Write_RPHTransfo(mole%RPHTransfo)
                write(out_unitp,*) ' Qact1',Qact1(:)
                write(out_unitp,*) 'ERROR in ',name_sub
                write(out_unitp,*) ' I cannot find Qact1(:) in tab_RPHpara_AT_Qact1'
                write(out_unitp,*) '  or  in tab_RPHpara_AT_Qref(1)'
                STOP
              END IF
           ELSE
             Vinact = HALF*sum(mole%RPHTransfo%tab_RPHpara_AT_Qact1(iQa)%dnehess%d0(:)*Qinact21(:)**2)
           END IF



            !write(6,*) 'iQa',iQa
            !write(6,*) 'Qinact21',Qinact21(:)
            !write(6,*) 'dnehess',mole%RPHTransfo%tab_RPHpara_AT_Qact1(iQa)%dnehess%d0(:)
            !write(6,*) 'Vinact',Vinact


            DO ie=1,para_PES%nb_elec
              Tab_dnMatOp(iOpE)%tab_dnMatOp(ie,ie,itermE)%d0 =          &
                       Tab_dnMatOp(iOpE)%tab_dnMatOp(ie,ie,itermE)%d0 + Vinact
            END DO

            CALL dealloc_NParray(Qact1,"Qact1",name_sub)
            CALL dealloc_NParray(Qinact21,"Qinact21",name_sub)

          END IF


        ELSE ! with finite difference (even for on-the-fly calculation)
          IF (debug) write(out_unitp,*) 'With numerical derivatives'

            CALL dnOp_num_grid_v2(Qact,Tab_dnMatOp,mole,para_Tnum,para_PES,nderiv_loc)

        END IF
      END IF

      DO ie=1,para_PES%nb_elec
        para_PES%min_pot = min(para_PES%min_pot,Tab_dnMatOp(iOpE)%tab_dnMatOp(ie,ie,itermE)%d0)
        para_PES%max_pot = max(para_PES%max_pot,Tab_dnMatOp(iOpE)%tab_dnMatOp(ie,ie,itermE)%d0)
      END DO

      !----------------------------------------------------------------
      IF (mole%Rot_Dip_with_EC .AND. nderivScal > -1 .AND.              &
                                         para_PES%nb_scalar_Op > 2) THEN
        CALL alloc_dnSVM(dnXin,mole%ncart,mole%nb_act,nderiv=nderivScal)

        CALL sub_QactTOdnx(Qact,dnXin,mole,nderiv=nderivScal,           &
                                    Gcenter=.TRUE.,Cart_transfo=.FALSE.)

        IF (debug) write(out_unitp,*) ' WARNING dip(:) are rotated'

        ! initial rotation of the dipole moment
        CALL alloc_MatOFdnS(dnT,dnXin%nb_var_deriv,nderivScal)
        CALL sub_ZERO_TO_MatOFdnS(dnT)
        dnT(:,:)%d0 = mole%tab_Cart_transfo(1)%CartesianTransfo%Rot_initial

        CALL alloc_VecOFdnS(dnScal_loc1,dnXin%nb_var_deriv,nderivScal)
        CALL alloc_VecOFdnS(dnScal_loc2,dnXin%nb_var_deriv,nderivScal)

        DO ie=1,para_PES%nb_elec
        DO je=1,para_PES%nb_elec
          DO i=1,para_PES%nb_scalar_Op
            iterm = Tab_dnMatOp(iOpScal-1+i)%derive_term_TO_iterm(0,0)
            CALL sub_dnS1_TO_dnS2(                                      &
                          Tab_dnMatOp(iOpScal-1+i)%tab_dnMatOp(ie,je,iterm), &
                                              dnScal_loc1(i),nderivScal)
          END DO

          CALL Mat1OFdnS_MUL_Vec2OFdnS_TO_Vec3OFdnS(dnT,dnScal_loc1,    &
                                                dnScal_loc2,nderivScal)

          DO i=1,para_PES%nb_scalar_Op
            iterm = Tab_dnMatOp(iOpScal-1+i)%derive_term_TO_iterm(0,0)
             CALL sub_dnS1_TO_dnS2(dnScal_loc2(i),                      &
                 Tab_dnMatOp(iOpScal-1+i)%tab_dnMatOp(ie,je,iterm),nderivScal)
          END DO

        END DO
        END DO
        ! End initial rotation

        ! Eckart rotation of the dipole moment
        !  1st: rotate dnXin with the initial rotation
        CALL alloc_dnSVM(dnXout,dnXin%nb_var_vec,dnXin%nb_var_deriv,nderivScal)
        CALL calc_dnTxdnXin_TO_dnXout(dnXin,dnT,dnXout,                 &
                   mole%tab_Cart_transfo(1)%CartesianTransfo,nderivScal)
        CALL sub_dnVec1_TO_dnVec2(dnXout,dnXin)

        !  2d: calculation of the dnXref
        CALL alloc_MatOFdnS(dnXref,dnXin%nb_var_deriv,nderivScal)

        CALL dnMWX_MultiRef(dnXref,                                     &
                             mole%tab_Cart_transfo(1)%CartesianTransfo, &
              Qact(1:mole%tab_Cart_transfo(1)%CartesianTransfo%nb_Qact),&
                                                                  dnXin)

        !  3d: calculation of the dnT (Eckart rotational matrix)
        CALL calc_dnTEckart(dnXin,dnT,dnXref,                           &
                   mole%tab_Cart_transfo(1)%CartesianTransfo,nderivScal)
        IF (debug) write(out_unitp,*) 'dnT%d0',Qact(1:mole%nb_act1),dnT(:,:)%d0
        CALL dealloc_MatOFdnS(dnXref)

        ! 4th: rotation of the dipole moment
        DO ie=1,para_PES%nb_elec
        DO je=1,para_PES%nb_elec
          DO i=1,para_PES%nb_scalar_Op
            iterm = Tab_dnMatOp(iOpScal-1+i)%derive_term_TO_iterm(0,0)
            CALL sub_dnS1_TO_dnS2(                                      &
                          Tab_dnMatOp(iOpScal-1+i)%tab_dnMatOp(ie,je,iterm), &
                                              dnScal_loc1(i),nderivScal)
          END DO


          CALL Mat1OFdnS_MUL_Vec2OFdnS_TO_Vec3OFdnS(dnT,dnScal_loc1,    &
                                                 dnScal_loc2,nderivScal)

          DO i=1,para_PES%nb_scalar_Op
            iterm = Tab_dnMatOp(iOpScal-1+i)%derive_term_TO_iterm(0,0)
             CALL sub_dnS1_TO_dnS2(dnScal_loc2(i),                      &
                 Tab_dnMatOp(iOpScal-1+i)%tab_dnMatOp(ie,je,iterm),nderivScal)
          END DO

        END DO
        END DO
        ! End initial rotation


        CALL dealloc_VecOFdnS(dnScal_loc1)
        CALL dealloc_VecOFdnS(dnScal_loc2)
        CALL dealloc_MatOFdnS(dnT)
        CALL dealloc_dnSVM(dnXin)
        CALL dealloc_dnSVM(dnXout)

      END IF
      !----------------------------------------------------------------

      IF (allocated(Qit)) CALL dealloc_NParray(Qit,'Qit',name_sub)

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'tab_dnMatOp'
        CALL Write_Tab_OF_dnMatOp(Tab_dnMatOp)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE get_dnMatOp_AT_Qact


      SUBROUTINE dnOp_num_grid_v2(Qact,Tab_dnMatOp,                     &
                                  mole,para_Tnum,para_PES,nderiv)
      USE mod_system
      !$ USE omp_lib, only : OMP_GET_THREAD_NUM
      USE mod_dnSVM
      USE mod_Coord_KEO
      USE mod_SimpleOp
      USE mod_PrimOp_def
      USE mod_OTF
      IMPLICIT NONE

!----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix)   :: mole
      TYPE (Tnum)      :: para_Tnum

      real (kind=Rkind), intent(inout) :: Qact(:)
      TYPE (param_PES) :: para_PES

      TYPE (param_dnMatOp), intent(inout) :: Tab_dnMatOp(:)
      integer, intent(in)                 :: nderiv



      real (kind=Rkind), allocatable    :: Qact_th(:,:)
      TYPE (param_d0MatOp), allocatable :: d0MatOp_th(:,:)

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------


!----- pour les derivees ---------------------------------------------
      real (kind=Rkind) ::    step2,stepp,step24


!----- working variables ---------------------------------------------
      integer           :: i,j,k,ie,je,nderiv_loc,ith,iOp,nb_Op
      real (kind=Rkind) :: vi,vj,poti,potij
      integer           :: nb_thread


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='dnOp_num_grid_v2'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'Qact',Qact
         write(out_unitp,*) 'stepOp',para_PES%stepOp
         write(out_unitp,*) 'nb_elec',para_PES%nb_elec

         CALL flush_perso(out_unitp)
       END IF
!-----------------------------------------------------------
      nderiv_loc = min(nderiv,Tab_dnMatOp(1)%tab_dnMatOp(1,1,1)%nderiv) ! We assume that all nderiv are identical!!

      IF (para_PES%stepOp <= ZERO) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' stepOp is <= 0',para_PES%stepOp
        STOP
      END IF
      step2 = ONE/(para_PES%stepOp*para_PES%stepOp)
      step24 = step2*HALF*HALF
      stepp = ONE/(para_PES%stepOp+para_PES%stepOp)

      IF (Grid_omp == 0) THEN
        nb_thread = 1
      ELSE
        nb_thread = Grid_maxth
      END IF
      IF (print_level > 1) write(out_unitp,*) 'nb_thread in ',name_sub,' : ',nb_thread
      CALL flush_perso(out_unitp)



      nb_Op = size(Tab_dnMatOp)

      allocate(d0MatOp_th(nb_Op,nb_thread))
      DO ith=1,nb_thread
        CALL Init_Tab_OF_d0MatOp(d0MatOp_th(:,ith),mole%nb_act,para_PES%nb_elec, &
                                 para_PES%Type_HamilOp,JRot=para_Tnum%JJ, &
                                 cplx=para_PES%pot_cplx,direct_KEO=para_PES%direct_KEO) ! H
      END DO

      !-- pot0 Qact(i) ------------------
      CALL get_d0MatOp_AT_Qact(Qact,d0MatOp_th(:,1),mole,para_Tnum,para_PES)

      DO k=1,nb_Op
        CALL d0MatOp_TO_dnMatOp(d0MatOp_th(k,1),Tab_dnMatOp(k),(/0,0/))
        IF (nderiv_loc > 1) THEN ! diagonal hessian (the -2f(0) contribution)
          DO i=1,mole%nb_act
            CALL d0MatOp_TO_dnMatOp(d0MatOp_th(k,1),Tab_dnMatOp(k),(/i,i/))
            CALL WeightDer_dnMatOp(Tab_dnMatOp(k),-TWO,(/i,i/))
          END DO
        END IF
      END DO


      allocate(Qact_th(size(Qact),nb_thread))
      DO ith=1,nb_thread
        Qact_th(:,ith) = Qact(:)
      END DO

!$OMP   PARALLEL &
!$OMP   DEFAULT(NONE) &
!$OMP   SHARED(mole,para_Tnum,para_PES,stepp,step2) &
!$OMP   SHARED(tab_dnMatOp,nderiv_loc,Qact,Qact_th,d0MatOp_th,nb_Op) &
!$OMP   PRIVATE(i,ie,je,ith) &
!$OMP   NUM_THREADS(nb_thread)

!$OMP   DO SCHEDULE(STATIC)
!-----------------------------------------------------------------
!----- d/Qqi et d2/dQi2 of pot0 ----------------------------------
!-----------------------------------------------------------------
      DO i=1,mole%nb_act

        ith = 1
        !$ ith = omp_get_thread_num()+1

        !-- pot0 Qact(i)+step -------------
        Qact_th(i,ith) = Qact(i) + para_PES%stepOp

        CALL get_d0MatOp_AT_Qact(Qact_th(:,ith),d0MatOp_th(:,ith),mole,para_Tnum,para_PES)


        DO k=1,size(Tab_dnMatOp)
          IF (nderiv_loc > 0) THEN ! gradient
            CALL d0MatOp_TO_dnMatOp(d0MatOp_th(k,ith),Tab_dnMatOp(k),(/i,0/))
          END IF
          IF (nderiv_loc > 1) THEN ! hessian (diagonal)
            CALL d0MatOp_wADDTO_dnMatOp(d0MatOp_th(k,ith),Tab_dnMatOp(k),(/i,i/),ONE)
          END IF
        END DO

        !-- pot0 Qact(i)-step -------------
        Qact_th(i,ith) = Qact(i) - para_PES%stepOp

        CALL get_d0MatOp_AT_Qact(Qact_th(:,ith),d0MatOp_th(:,ith),mole,para_Tnum,para_PES)

        DO k=1,size(Tab_dnMatOp)
          IF (nderiv_loc > 0) THEN ! gradient
            CALL d0MatOp_wADDTO_dnMatOp(d0MatOp_th(k,ith),Tab_dnMatOp(k),(/i,0/),-ONE)
            CALL WeightDer_dnMatOp(Tab_dnMatOp(k),stepp,(/i,0/))
          END IF
          IF (nderiv_loc > 1) THEN ! diagonal hessian
            CALL d0MatOp_wADDTO_dnMatOp(d0MatOp_th(k,ith),Tab_dnMatOp(k),(/i,i/),ONE)
            CALL WeightDer_dnMatOp(Tab_dnMatOp(k),step2,(/i,i/))
          END IF
        END DO

        Qact_th(i,ith) = Qact(i)

      END DO
!-----------------------------------------------------------------
!----- end d/Qqi and d2/dQi2 of pot0 -----------------------------
!-----------------------------------------------------------------
!$OMP   END DO
!$OMP   END PARALLEL

      IF (nderiv_loc > 1) THEN ! hessian (off diagonal)
!$OMP   PARALLEL &
!$OMP   DEFAULT(NONE) &
!$OMP   SHARED(mole,para_Tnum,para_PES,stepp,step2,step24) &
!$OMP   SHARED(tab_dnMatOp,nderiv_loc,Qact,Qact_th,d0MatOp_th,nb_Op) &
!$OMP   PRIVATE(i,j,ie,je,ith) &
!$OMP   NUM_THREADS(nb_thread)

!$OMP   DO SCHEDULE(STATIC)
!-----------------------------------------------------------------
!----- d2/dQidQj of pot0 (4 points) -----------------------
!      d2/dQidQj = ( v(Qi+,Qj+)+v(Qi-,Qj-)-v(Qi-,Qj+)-v(Qi+,Qj-) )/(4*s*s)
!-----------------------------------------------------------------
      DO i=1,mole%nb_act
      DO j=i+1,mole%nb_act

        ith = 1
        !$ ith = omp_get_thread_num()+1

        !-- pot0 at Qact(i)+step Qact(j)+step
        Qact_th(i,ith) = Qact(i) + para_PES%stepOp
        Qact_th(j,ith) = Qact(j) + para_PES%stepOp

        CALL get_d0MatOp_AT_Qact(Qact_th(:,ith),d0MatOp_th(:,ith),mole,para_Tnum,para_PES)

        DO k=1,size(Tab_dnMatOp)
          CALL d0MatOp_TO_dnMatOp(d0MatOp_th(k,ith),Tab_dnMatOp(k),(/i,j/))
        END DO

        !-- pot0 at Qact(i)-step Qact(j)-step
        Qact_th(i,ith) = Qact(i) - para_PES%stepOp
        Qact_th(j,ith) = Qact(j) - para_PES%stepOp

        CALL get_d0MatOp_AT_Qact(Qact_th(:,ith),d0MatOp_th(:,ith),mole,para_Tnum,para_PES)

        DO k=1,size(Tab_dnMatOp)
          CALL d0MatOp_wADDTO_dnMatOp(d0MatOp_th(k,ith),Tab_dnMatOp(k),(/i,j/),ONE)
        END DO


        !-- pot0 at Qact(i)-step Qact(j)+step
        Qact_th(i,ith) = Qact(i) - para_PES%stepOp
        Qact_th(j,ith) = Qact(j) + para_PES%stepOp

        CALL get_d0MatOp_AT_Qact(Qact_th(:,ith),d0MatOp_th(:,ith),mole,para_Tnum,para_PES)

        DO k=1,size(Tab_dnMatOp)
          CALL d0MatOp_wADDTO_dnMatOp(d0MatOp_th(k,ith),Tab_dnMatOp(k),(/i,j/),-ONE)
        END DO


        !-- pot0 at Qact(i)+step Qact(j)-step
        Qact_th(i,ith) = Qact(i) + para_PES%stepOp
        Qact_th(j,ith) = Qact(j) - para_PES%stepOp

        CALL get_d0MatOp_AT_Qact(Qact_th(:,ith),d0MatOp_th(:,ith),mole,para_Tnum,para_PES)

        DO k=1,size(Tab_dnMatOp)
          CALL d0MatOp_wADDTO_dnMatOp(d0MatOp_th(k,ith),Tab_dnMatOp(k),(/i,j/),-ONE)
          CALL WeightDer_dnMatOp(Tab_dnMatOp(k),step24,(/i,j/))
          CALL dnMatOp2Der_TO_dnMatOp1Der(Tab_dnMatOp(k),(/j,i/),Tab_dnMatOp(k),(/i,j/))
        END DO

        Qact_th(i,ith) = Qact(i)
        Qact_th(j,ith) = Qact(j)

      END DO
      END DO
!-----------------------------------------------------------------
!----- end d2/dQidQj of pot0 -------------------------------------
!-----------------------------------------------------------------
!$OMP   END DO
!$OMP   END PARALLEL
      END IF

      DO ith=1,nb_thread
      DO k=1,nb_Op
        CALL dealloc_d0MatOp(d0MatOp_th(k,ith))
      END DO
      END DO
      deallocate(d0MatOp_th)

      deallocate(Qact_th)

      IF (debug) THEN
        write(out_unitp,*) 'Tab_dnMatOp'
        CALL Write_Tab_OF_dnMatOp(Tab_dnMatOp)
        write(out_unitp,*) 'END ',name_sub
      END IF

      END SUBROUTINE dnOp_num_grid_v2


      SUBROUTINE TnumKEO_TO_tab_dnH(Qact,Tab_dnH,mole,para_Tnum)
      USE mod_system
      USE mod_dnSVM
      USE mod_Coord_KEO, only : zmatrix,Tnum,get_dng_dnGG, sub3_dnrho_ana, &
                                calc3_f2_f1Q_num, sub3_dndetGG
      USE mod_SimpleOp
      IMPLICIT NONE

!----- for the zmatrix and Tnum --------------------------------------
      TYPE (Tnum)    :: para_Tnum
      TYPE (zmatrix) :: mole

!----- for Qact ... ---------------------------------------------
      real (kind=Rkind), intent(inout) :: Qact(:)

!----- input output variables ----------------------------------------
      TYPE (param_dnMatOp), intent(inout) :: Tab_dnH


!----- working variables -------------------------------------------------
      real (kind=Rkind) ::                                              &
                T2(mole%nb_act,mole%nb_act),T1(mole%nb_act),vep,rho,    &
                Tcor2(mole%nb_act,3),Tcor1(3),Trot(3,3)

      integer             :: i,j,ie,iterm
      TYPE(Type_dnS)      :: dnrho


!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='TnumKEO_TO_Tab_dnH'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Qact',Qact
      END IF
!-----------------------------------------------------------

  IF (Tab_dnH%type_Op == 1) THEN

      CALL calc3_f2_f1Q_num(Qact,T2,T1,vep,rho,Tcor2,Tcor1,Trot, &
                            para_Tnum,mole)

      IF (para_Tnum%nrho == 0) THEN
        Tab_dnH%Jac = rho ! because in calc3_f2_f1Q_num and when nrho=0, rho=jac

        CALL alloc_dnSVM(dnrho,mole%nb_act,nderiv=0)
        CALL sub3_dnrho_ana(dnrho,Qact,mole,0)
        Tab_dnH%rho = dnrho%d0 ! here we should have the good rho for the basis set
      END IF


      DO ie=1,Tab_dnH%nb_bie
        ! T2
        DO i=1,Tab_dnH%nb_Qact
        DO j=i,Tab_dnH%nb_Qact
          iterm = Tab_dnH%derive_term_TO_iterm(j,i)
          Tab_dnH%tab_dnMatOp(ie,ie,iterm)%d0 = T2(j,i)
        END DO
        END DO

        ! T1
        DO i=1,Tab_dnH%nb_Qact
          iterm = Tab_dnH%derive_term_TO_iterm(i,0)
          Tab_dnH%tab_dnMatOp(ie,ie,iterm)%d0 = T1(i)
        END DO

        ! vep is added
        iterm = Tab_dnH%derive_term_TO_iterm(0,0)
        Tab_dnH%tab_dnMatOp(ie,ie,iterm)%d0 =                           &
                               Tab_dnH%tab_dnMatOp(ie,ie,iterm)%d0 + vep

        ! rot
        DO i=-3,-1
        DO j=i,-1
          iterm = Tab_dnH%derive_term_TO_iterm(j,i)
          Tab_dnH%tab_dnMatOp(ie,ie,iterm)%d0 = Trot(-j,-i)
        END DO
        END DO

        ! cor
        DO i=-3,-1
          iterm = Tab_dnH%derive_term_TO_iterm(i,0)
          Tab_dnH%tab_dnMatOp(ie,ie,iterm)%d0 = Tcor1(-i)
          DO j=1,Tab_dnH%nb_Qact
            iterm = Tab_dnH%derive_term_TO_iterm(j,i)
            Tab_dnH%tab_dnMatOp(ie,ie,iterm)%d0 = Tcor2(j,-i)
          END DO
        END DO

      END DO
  ELSE
    write(out_unitp,*) 'ERROR in ',name_sub
    write(out_unitp,*) ' This type of Op is not possible for the KEO',Tab_dnH%type_Op
    write(out_unitp,*) ' The possibilities are: 1 or 10'
    write(out_unitp,*) '    .... 10 not yet !!!!!!'

    write(out_unitp,*) '   Check the fortran!!'
    STOP
  END IF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Tab_dnH'
        CALL Write_dnMatOp(Tab_dnH)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE TnumKEO_TO_tab_dnH

      SUBROUTINE TnumKEO_TO_tab_d0H(Qact,d0MatH,mole,para_Tnum)
      USE mod_system
      USE mod_dnSVM
      USE mod_Coord_KEO, only : zmatrix,Tnum,get_dng_dnGG, sub3_dnrho_ana, &
                                calc3_f2_f1Q_num, sub3_dndetGG
      USE mod_SimpleOp
      IMPLICIT NONE

!----- for the zmatrix and Tnum --------------------------------------
      TYPE (Tnum)    :: para_Tnum
      TYPE (zmatrix) :: mole

!----- for Qact ... ---------------------------------------------
      real (kind=Rkind), intent(inout) :: Qact(:)

!----- input output variables ----------------------------------------
      TYPE (param_d0MatOp), intent(inout) :: d0MatH


!----- working variables -------------------------------------------------
      real (kind=Rkind) ::                                              &
                T2(mole%nb_act,mole%nb_act),T1(mole%nb_act),vep,rho,    &
                Tcor2(mole%nb_act,3),Tcor1(3),Trot(3,3)

      TYPE(Type_dnMat) :: dnGG


      integer          :: i,j,ie,iterm
      TYPE(Type_dnS)   :: dnrho


!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='TnumKEO_TO_Tab_d0H'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Qact',Qact
      END IF
!-----------------------------------------------------------
  IF (d0MatH%type_Op == 1) THEN
      !----------------------------------------------------------------
      CALL calc3_f2_f1Q_num(Qact,T2,T1,vep,rho,Tcor2,Tcor1,Trot, &
                            para_Tnum,mole)

      IF (para_Tnum%nrho == 0) THEN
        d0MatH%Jac = rho ! because in calc3_f2_f1Q_num and when nrho=0, rho=jac

        CALL alloc_dnSVM(dnrho,mole%nb_act,nderiv=0)
        CALL sub3_dnrho_ana(dnrho,Qact,mole,0)
        d0MatH%rho = dnrho%d0 ! here we should have the good rho for the basis set
        CALL dealloc_dnSVM(dnrho)

      END IF

      !write(6,*) 'nrho,Qact,vep',para_Tnum%nrho,Qact(1:mole%nb_act),vep
      !IF (d0MatH%ReVal(1,1,1) < -0.01_Rkind) write(6,*) 'V',d0MatH%ReVal(1,1,1)

      DO ie=1,d0MatH%nb_bie
        ! T2
        DO i=1,d0MatH%nb_Qact
        DO j=i,d0MatH%nb_Qact
          iterm = d0MatH%derive_term_TO_iterm(j,i)
          d0MatH%ReVal(ie,ie,iterm) = T2(j,i)
        END DO
        END DO

        ! T1
        DO i=1,d0MatH%nb_Qact
          iterm = d0MatH%derive_term_TO_iterm(i,0)
          d0MatH%ReVal(ie,ie,iterm) = T1(i)
        END DO

        ! vep is added
        iterm = d0MatH%derive_term_TO_iterm(0,0)
        d0MatH%ReVal(ie,ie,iterm) = d0MatH%ReVal(ie,ie,iterm) + vep

        IF (d0MatH%JRot > 0) THEN
          ! rot
          DO i=-3,-1
          DO j=i,-1
            iterm = d0MatH%derive_term_TO_iterm(j,i)
            d0MatH%ReVal(ie,ie,iterm) = Trot(-j,-i)
          END DO
          END DO

          ! cor
          DO i=-3,-1
            iterm = d0MatH%derive_term_TO_iterm(i,0)
            d0MatH%ReVal(ie,ie,iterm) = Tcor1(-i)
            DO j=1,d0MatH%nb_Qact
              iterm = d0MatH%derive_term_TO_iterm(j,i)
              d0MatH%ReVal(ie,ie,iterm) = Tcor2(j,-i)
            END DO
          END DO
        END IF

      END DO
  ELSE IF (d0MatH%type_Op == 10) THEN

      CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,0)
      CALL get_dng_dnGG(Qact,para_Tnum,mole,dnGG=dnGG,nderiv=0)

      IF (.NOT. d0MatH%direct_KEO) THEN

      DO ie=1,d0MatH%nb_bie

        DO i=1,d0MatH%nb_Qact
        DO j=1,d0MatH%nb_Qact
          iterm = d0MatH%derive_term_TO_iterm(j,i)
          d0MatH%ReVal(ie,ie,iterm) = dnGG%d0(j,i)
        END DO
        END DO

        IF (d0MatH%JRot > 0) THEN
          DO iterm=1,3*d0MatH%nb_Qact*2
            i = d0MatH%derive_termQact(1,iterm)
            IF (i < 0) i = d0MatH%nb_Qact-i
            j = d0MatH%derive_termQact(2,iterm)
            IF (j < 0) j = d0MatH%nb_Qact-j

            d0MatH%ReVal(ie,ie,iterm) = dnGG%d0(j,i)
          END DO
        END IF

      END DO
      END IF

      CALL alloc_dnSVM(dnrho,mole%nb_act,nderiv=0)

      ! here dnrho contains dnJac
      CALL sub3_dndetGG(dnrho,dnGG,0,mole%masses,mole%Mtot_inv,mole%ncart)
      d0MatH%Jac = dnrho%d0

      ! now the true dnrho
      CALL sub3_dnrho_ana(dnrho,Qact,mole,0)
      d0MatH%rho = dnrho%d0 ! here we should have the good rho for the basis set

      CALL dealloc_dnSVM(dnrho)
      CALL dealloc_dnSVM(dnGG)


  ELSE
    write(out_unitp,*) 'ERROR in ',name_sub
    write(out_unitp,*) ' This type of Op is not possible for the KEO',d0MatH%type_Op
    write(out_unitp,*) ' The possibilities are: 1 or 10'
    write(out_unitp,*) '   Check the fortran!!'
    STOP
  END IF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'd0MatH'
          CALL Write_d0MatOp(d0MatH)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------
!STOP 'TnumKEO_TO_tab_d0H'
      END SUBROUTINE TnumKEO_TO_tab_d0H

!=============================================================
!
!     frequency calculations at Qact
!
!=============================================================
      SUBROUTINE sub_freq_AT_Qact(freq,Qact,para_Tnum,mole,para_PES,print_freq,d0h_opt)
      USE mod_system
      USE mod_dnSVM
      USE mod_Constant
      USE mod_Coord_KEO, only : zmatrix,Tnum,get_dng_dnGG, calc_freq
      USE mod_SimpleOp
      USE mod_PrimOp_def
      IMPLICIT NONE

!----- for the zmatrix and Tnum --------------------------------------
      TYPE (Tnum)        :: para_Tnum
      TYPE (zmatrix)     :: mole
      TYPE (param_PES)   :: para_PES
      real (kind=Rkind), intent(inout) :: Qact(:)
      logical, intent(in), optional :: print_freq
      real (kind=Rkind), optional :: d0h_opt(:,:)


      real (kind=Rkind), intent(inout) :: freq(mole%nb_act)



      TYPE (param_dnMatOp) :: dnMatOp(1)

      TYPE(Type_dnMat) :: dnGG

      real (kind=Rkind), allocatable :: d0c_inv(:,:)
      real (kind=Rkind), allocatable :: d0c_ini(:,:)
      real (kind=Rkind), allocatable :: d0k(:,:)
      real (kind=Rkind), allocatable :: d0h(:,:),d0grad(:)
      real (kind=Rkind), allocatable :: d0c(:,:)
      real (kind=Rkind) :: norme
      integer           :: i,i2
      logical           :: print_freq_loc

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_freq_AT_Qact'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
     IF (present(print_freq)) THEN
       print_freq_loc = print_freq
     ELSE
       print_freq_loc = .FALSE.
     END IF

      IF (debug .OR. print_freq_loc) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Qact',Qact
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------
      CALL alloc_NParray(d0c,    (/ mole%nb_act,mole%nb_act /),"d0c",    name_sub)
      CALL alloc_NParray(d0k,    (/ mole%nb_act,mole%nb_act /),"d0k",    name_sub)
      CALL alloc_NParray(d0h,    (/ mole%nb_act,mole%nb_act /),"d0h",    name_sub)
      CALL alloc_NParray(d0grad, (/ mole%nb_act /),            "d0grad", name_sub)
      CALL alloc_NParray(d0c_inv,(/ mole%nb_act,mole%nb_act /),"d0c_inv",name_sub)
      CALL alloc_NParray(d0c_ini,(/ mole%nb_act,mole%nb_act /),"d0c_ini",name_sub)


      !----- Hessian ------------------------------------
      IF (present(d0h_opt)) THEN
        d0h = d0h_opt
      ELSE
        CALL Init_Tab_OF_dnMatOp(dnMatOp,mole%nb_act,1,nderiv=2)

        CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole,para_Tnum,para_PES)

        CALL Get_Hess_FROM_Tab_OF_dnMatOp(d0h,dnMatOp)
        CALL Get_Grad_FROM_Tab_OF_dnMatOp(d0grad,dnMatOp)

        IF (debug .OR. print_freq_loc) THEN
          write(out_unitp,*) 'Energy:',Get_Scal_FROM_Tab_OF_dnMatOp(dnMatOp)
          write(out_unitp,*) 'Gradient:'
          CALL Write_VecMat(d0grad,out_unitp,5)
          write(out_unitp,*) 'Hessian:'
          CALL Write_VecMat(d0h,out_unitp,5)
        END IF

        CALL dealloc_Tab_OF_dnMatOp(dnMatOp)
      END IF

      !----- kinetic part ---------------------------------
      CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,0)

      para_Tnum%WriteT    = .FALSE.
      CALL get_dng_dnGG(Qact,para_Tnum,mole,dnGG=dnGG,nderiv=0)

      d0c_ini(:,:) = ZERO
      d0k = dnGG%d0(1:mole%nb_act,1:mole%nb_act)

      CALL dealloc_dnSVM(dnGG)

      IF (debug .OR. print_freq_loc) THEN
        write(out_unitp,*) 'Curvilinear kinetic part:'
        CALL Write_VecMat(d0k,out_unitp,5)
      END IF

      !----- frequencies ---------------------------------
      CALL calc_freq(mole%nb_act,d0h,d0k,freq,d0c,d0c_inv,norme,d0c_ini,.FALSE.)

      IF (debug .OR. print_freq_loc) THEN
        write(out_unitp,*) 'd0c (NM?):'
        CALL Write_VecMat(d0c,out_unitp,5)
      END IF

      CALL dealloc_NParray(d0c,    "d0c",    name_sub)
      CALL dealloc_NParray(d0k,    "d0k",    name_sub)
      CALL dealloc_NParray(d0h,    "d0h",    name_sub)
      CALL dealloc_NParray(d0grad, "d0grad", name_sub)
      CALL dealloc_NParray(d0c_inv,"d0c_inv",name_sub)
      CALL dealloc_NParray(d0c_ini,"d0c_ini",name_sub)

      IF (debug .OR. print_freq_loc) THEN
        write(out_unitp,*) 'ZPE (cm-1): ',HALF*sum(freq(:))*get_Conv_au_TO_unit('E','cm-1')
        write(out_unitp,*) 'ZPE   (eV): ',HALF*sum(freq(:))*get_Conv_au_TO_unit('E','eV')
        write(out_unitp,*) 'ZPE   (au): ',HALF*sum(freq(:))

        DO i=1,size(freq),3
          i2 = min(i+2,mole%nb_act)
          write(out_unitp,'(a,i0,"-",i0,3(x,f0.4))') 'frequencies (cm-1): ',  &
                                 i,i2,freq(i:i2) * get_Conv_au_TO_unit('E','cm-1')
        END DO

        write(out_unitp,*) 'END ',name_sub
      END IF
      CALL flush_perso(out_unitp)
!-----------------------------------------------------------

      END SUBROUTINE sub_freq_AT_Qact
!=====================================================================
!
! ++   calculation gaussian_width and freq with curvilinear coordinates
!
!=====================================================================
      SUBROUTINE calc3_NM_TO_sym(Qact,mole,para_Tnum,para_PES,hCC,l_hCC)
      USE mod_system
      USE mod_dnSVM
      USE mod_Constant, only : get_Conv_au_TO_unit
      USE mod_Coord_KEO
      USE mod_SimpleOp
      USE mod_PrimOp_def
      IMPLICIT NONE

      TYPE (zmatrix) :: mole,mole_1
      TYPE (Tnum)    :: para_Tnum

      real (kind=Rkind), intent(inout) :: Qact(:)
      TYPE (param_PES)   :: para_PES
      real (kind=Rkind)  :: hCC(mole%ncart_act,mole%ncart_act)
      logical            :: l_hCC  ! if .TRUE. hCC is already calculated (for PVSCF)


      TYPE(Type_dnMat) :: dnGG

      TYPE(Type_dnS)       :: dnECC(1,1),dnE(1,1)
      TYPE (param_dnMatOp) :: dnMatOp(1)

      real (kind=Rkind), allocatable :: d0k(:,:),d0k_save(:,:),d0h(:,:),d0grad(:)


      real (kind=Rkind), allocatable :: d0c_inv(:,:),d0c_ini(:,:)
      real (kind=Rkind), allocatable :: d0c(:,:),d0eh(:),ScalePara(:),tab_sort(:)
      real (kind=Rkind), allocatable :: mat_inv(:,:),mat(:,:)

      real (kind=Rkind) :: max_freq,norme,step = ONETENTH**3

      real (kind=Rkind) :: auTOcm_inv

      integer :: i,j,k,i_act,i_sym,k_act,k_sym,ierr,nb_NM

      logical                  :: Read_OnTheFly_only,OnTheFly
      character (len=Line_len) :: name_FChk


!      -----------------------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug=.FALSE.
!      logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'calc3_NM_TO_sym'
!      -----------------------------------------------------------------
       IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
        write(out_unitp,*) 'Qdyn0 =',mole%ActiveTransfo%Qdyn0(:)
        write(out_unitp,*)
!       CALL Write_mole(mole)
!       write(out_unitp,*)
       END IF
!      -----------------------------------------------------------------
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')

      write(out_unitp,*) '========================================='
      write(out_unitp,*) '========== calc3_NM_TO_sym =============='
      write(out_unitp,*) '========================================='
      write(out_unitp,*)

      Qact = mole%ActiveTransfo%Qact0(:)

      write(out_unitp,*) '========================================='
      write(out_unitp,*) '==== hessian and kinetic matrices ======='
      write(out_unitp,*) '========================================='

      IF (mole%NMTransfo%hessian_read .AND. mole%NMTransfo%k_read) THEN
        nb_NM = ubound(mole%NMTransfo%d0k,dim=1)
        CALL alloc_NParray(d0k,(/nb_NM,nb_NM/),"d0k",name_sub)
        d0k(:,:) = mole%NMTransfo%d0k(:,:)

        IF (nb_NM /= ubound(mole%NMTransfo%d0h,dim=1)) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' The dimension of d0k and d0h are different!'
          write(out_unitp,*) ' shape(d0k):',shape(mole%NMTransfo%d0k)
          write(out_unitp,*) ' shape(d0h):',shape(mole%NMTransfo%d0h)
          write(out_unitp,*) ' Check your data !!'
          STOP
        END IF
        CALL alloc_NParray(d0h,(/nb_NM,nb_NM/),"d0h",name_sub)
        d0h(:,:) = mole%NMTransfo%d0h(:,:)

        CALL alloc_NParray(d0grad,(/nb_NM/),"d0grad",name_sub)
        d0grad(:) = ZERO


      ELSE ! both are false
        !- create mole_1 (type=-1 => type=1)
        CALL mole1TOmole2(mole,mole_1)
        DO i=1,mole_1%nb_var
          IF (mole_1%ActiveTransfo%list_act_OF_Qdyn(i) == -1)           &
                            mole_1%ActiveTransfo%list_act_OF_Qdyn(i) = 1
        END DO
        IF (debug) THEN
          write(out_unitp,*) 'mole_1:'
          CALL Write_mole(mole_1)
        END IF

        !- calc G_1
        CALL alloc_dnSVM(dnGG,mole_1%ndimG,mole_1%ndimG,mole_1%nb_act,0)

        CALL get_dng_dnGG(Qact,para_Tnum,mole_1,dnGG=dnGG,nderiv=0)

        nb_NM = mole_1%nb_act

        CALL alloc_NParray(d0k,(/nb_NM,nb_NM/),"d0k",name_sub)
        d0k(:,:) = dnGG%d0(1:nb_NM,1:nb_NM)

        CALL dealloc_dnSVM(dnGG)

        !- calculation of the hessian (mole_1)
        CALL alloc_NParray(d0h,(/nb_NM,nb_NM/),"d0h",name_sub)
        d0h(:,:) = ZERO
        CALL alloc_NParray(d0grad,(/nb_NM/),"d0grad",name_sub)
        d0grad(:) = ZERO

        CALL alloc_MatOFdnS(dnECC,mole_1%ncart_act,2)
        CALL alloc_MatOFdnS(dnE,nb_NM,2)

        IF (mole%NMTransfo%hessian_old) THEN
          IF (mole%NMTransfo%hessian_onthefly) THEN

            mole%NMTransfo%hessian_cart = .TRUE.
            ! save on-the-fly parameters
            name_FChk          = para_PES%para_OTF%file_FChk%name
            Read_OnTheFly_only = para_PES%Read_OnTheFly_only
            OnTheFly           = para_PES%OnTheFly

            ! set-up on-the-fly parameters to read the hessian
            para_PES%OnTheFly                = .TRUE.
            para_PES%Read_OnTheFly_only      = .TRUE.
            para_PES%para_OTF%file_FChk%name = mole%NMTransfo%file_hessian%name

            write(out_unitp,*) 'Read ab initio hessian from file: ',    &
                                  trim(para_PES%para_OTF%file_FChk%name)

            !CALL  dnOp_grid(Qact,dnE,2,mole_1,para_Tnum,para_PES)
            CALL Init_Tab_OF_dnMatOp(dnMatOp,nb_NM,1,nderiv=2)
            CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole_1,para_Tnum,para_PES)
            CALL Get_Hess_FROM_Tab_OF_dnMatOp(d0h,dnMatOp)
            CALL Get_Grad_FROM_Tab_OF_dnMatOp(d0grad,dnMatOp)
            CALL dealloc_Tab_OF_dnMatOp(dnMatOp)


            ! restore the on-the-fly parameters
            para_PES%para_OTF%file_FChk%name = name_FChk
            para_PES%Read_OnTheFly_only      = Read_OnTheFly_only
            para_PES%OnTheFly                = OnTheFly
            !STOP 'coucou'
          ELSE
            IF (mole%NMTransfo%hessian_cart) THEN
              write(out_unitp,*) 'Old hessian : mole_1%ncart_act',mole_1%ncart_act
              dnECC(1,1)%d1(:)   = ZERO
              IF (.NOT. l_hCC) THEN
                dnECC(1,1)%d2(:,:)   = ZERO
                CALL sub_hessian(dnECC(1,1)%d2)
              ELSE
                dnECC(1,1)%d2(:,:) = hCC(:,:)
              END IF
              CALL sub_dnFCC_TO_dnFcurvi(Qact,dnECC(1,1),dnE(1,1),mole_1)
              d0h(:,:) = dnE(1,1)%d2(:,:)
            ELSE
              write(out_unitp,*) 'hessian : mole_1%nb_act',mole_1%nb_act
              CALL sub_hessian(d0h)
            END IF
          END IF
        ELSE

          !CALL  dnOp_grid(Qact,dnE,2,mole_1,para_Tnum,para_PES)
          CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole_1,para_Tnum,para_PES)
          CALL Get_Hess_FROM_Tab_OF_dnMatOp(d0h,dnMatOp)
          CALL Get_Grad_FROM_Tab_OF_dnMatOp(d0grad,dnMatOp)
          CALL dealloc_Tab_OF_dnMatOp(dnMatOp)

        END IF

        IF (debug) THEN
          write(out_unitp,*) 'Qref (Qact)',Qact
          write(out_unitp,*) 'pot_Qref',dnE%d0
        END IF

        write(out_unitp,*) 'gradient:'
        DO i=1,nb_NM
          write(out_unitp,*) 'Q',i,d0grad(i)
        END DO

        CALL dealloc_MatOFdnS(dnE)
        CALL dealloc_MatOFdnS(dnECC)
        CALL dealloc_zmat(mole_1)
        CALL dealloc_NParray(d0grad,"d0grad",name_sub)

      END IF

      !write with high precision to be able to read them
      write(out_unitp,*) 'hessian matrix'
      write(out_unitp,*) nb_NM,5
      CALL Write_Mat(d0h,out_unitp,5,Rformat='e20.13')
      write(out_unitp,*) 'kinetic matrix'
      write(out_unitp,*) nb_NM,5
      CALL Write_Mat(d0k,out_unitp,5,Rformat='e20.13')


      ! parameters for uncoupled HO
      CALL alloc_NParray(ScalePara,(/ nb_NM /),"ScalePara",name_sub)
      DO i=1,nb_NM
        ScalePara(i) = sqrt(sqrt(abs(d0h(i,i)/d0k(i,i))))
        !write(out_unitp,*) 'i,d0h,d0k',i,d0h(i,i),d0k(i,i)
        !write(out_unitp,*) 'i,ScalePara(i)',i,ScalePara(i)
      END DO
      ! parameters for uncoupled HO

!     - frequencies
      CALL alloc_NParray(d0c_inv,(/ nb_NM,nb_NM /),"d0c_inv",name_sub)
      CALL alloc_NParray(d0c_ini,(/ nb_NM,nb_NM /),"d0c_ini",name_sub)
      CALL alloc_NParray(d0c,(/ nb_NM,nb_NM /),"d0c",name_sub)
      CALL alloc_NParray(d0eh,(/ nb_NM /),"d0eh",name_sub)

      IF (mole%NMTransfo%purify_hess) THEN

        CALL H0_symmetrization(d0h,nb_NM,mole%NMTransfo%Qact1_sym,      &
                         mole%NMTransfo%dim_equi,mole%NMTransfo%tab_equi)
        write(out_unitp,*) 'purified hessian matrix'
        write(out_unitp,*) nb_NM,5
        CALL Write_Mat(d0h,out_unitp,5,Rformat='e20.13')

        CALL H0_symmetrization(d0k,nb_NM,mole%NMTransfo%Qact1_sym,      &
                        mole%NMTransfo%dim_equi,mole%NMTransfo%tab_equi)
        write(out_unitp,*) 'purified K (kinetic) matrix'
        write(out_unitp,*) nb_NM,5
        CALL Write_Mat(d0k,out_unitp,5,Rformat='e20.13')

      END IF
      write(out_unitp,*) '========================================='
      write(out_unitp,*)


      d0c_ini(:,:) = ZERO
      write(out_unitp,*) '========================================='
      write(out_unitp,*) '======== frequencies ===================='
      write(out_unitp,*) '========================================='

      CALL alloc_NParray(d0k_save,(/ nb_NM,nb_NM /),"d0k_save",name_sub)
      d0k_save(:,:) = d0k(:,:)
      IF (mole%NMTransfo%d0c_read) THEN
        d0c(:,:) = mole%NMTransfo%d0c(:,:)
        CALL dealloc_array(mole%NMTransfo%d0c,"mole%NMTransfo%d0c",name_sub)

        CALL calc_freq_WITH_d0c(nb_NM,d0h,d0k_save,d0eh,                &
                                d0c,d0c_inv,norme)
      ELSE
        CALL calc_freq(nb_NM,d0h,d0k_save,d0eh,                         &
                       d0c,d0c_inv,norme,d0c_ini,.FALSE.)

        !write with high precision to be able to read it
        write(out_unitp,*) 'd0c'
        write(out_unitp,*) nb_NM,5
        CALL Write_Mat(d0c,out_unitp,5,Rformat='e20.13')
      END IF

      write(out_unitp,*) '========================================='
      IF (associated(mole%NMTransfo%Qact1_sym)) THEN
        IF (debug) write(out_unitp,*) '   d0eh,d0c,d0c_inv after "sort_with_Tab"'

        CALL alloc_NParray(tab_sort,(/ nb_NM /),"tab_sort",name_sub)
        tab_sort(:) = ZERO
        max_freq = maxval(d0eh(:))
        DO i=1,nb_NM
          k = maxloc(abs(d0c(:,i)),dim=1)
          tab_sort(i) = real(mole%NMTransfo%Qact1_sym(k),kind=Rkind)
          IF (tab_sort(i) > ZERO) tab_sort(i) = d0eh(i) + tab_sort(i) * max_freq
        END DO
        write(out_unitp,*) 'tab_sort: ',tab_sort(:)

        CALL sort_with_Tab(d0c,d0c_inv,d0eh,tab_sort,nb_NM)

        CALL dealloc_NParray(tab_sort,"tab_sort",name_sub)

      ELSE
        IF (debug) write(out_unitp,*) '   d0eh,d0c,d0c_inv'
      END IF

      IF (debug) THEN
        write(out_unitp,*) 'frequencies (cm-1): ',d0eh(:)*auTOcm_inv
        write(out_unitp,*) 'd0c'
        CALL Write_Mat(d0c,out_unitp,5)
        write(out_unitp,*) 'd0c_inv'
        CALL Write_Mat(d0c_inv,out_unitp,5)
      END IF

      IF (mole%NMTransfo%k_Half) THEN
        write(out_unitp,*) '==========================='
        d0k_save = matmul(transpose(d0c),matmul(d0k,d0c))
        write(out_unitp,*) 'new d0k'
        CALL Write_Mat(d0k_save,out_unitp,5)
        DO i=1,nb_NM
          d0c(:,i)     = d0c(:,i)     / sqrt(d0k_save(i,i))
          d0c_inv(i,:) = d0c_inv(i,:) * sqrt(d0k_save(i,i))
        END DO

        write(out_unitp,*) '==========================='
      END IF

      mole%NMTransfo%nb_NM = nb_NM

      CALL alloc_array(mole%NMTransfo%d0c,(/ nb_NM,nb_NM /),            &
                      "mole%NMTransfo%d0c",name_sub)
      mole%NMTransfo%d0c(:,:)     = d0c(:,:)

      CALL alloc_array(mole%NMTransfo%d0c_inv,(/ nb_NM,nb_NM /),        &
                      "mole%NMTransfo%d0c_inv",name_sub)
      mole%NMTransfo%d0c_inv(:,:) = d0c_inv(:,:)

      CALL alloc_array(mole%NMTransfo%d0eh,(/ nb_NM /),                 &
                      "mole%NMTransfo%d0eh",name_sub)
      mole%NMTransfo%d0eh(:)      = d0eh(:)

      write(out_unitp,*) 'frequencies (cm-1):',d0eh(:)*auTOcm_inv
      !write(out_unitp,*) 'all scaling Gaussian  :',sqrt(abs(d0eh))
      !write(out_unitp,*) 'd0c'
      !CALL Write_Mat(d0c,out_unitp,5)
      !write(out_unitp,*) 'd0c_inv'
      !CALL Write_Mat(d0c_inv,out_unitp,5)

      write(out_unitp,*) '========================================='
      write(out_unitp,*)

      CALL alloc_NParray(mat,    (/ mole%nb_var,mole%nb_var /),"mat",name_sub)
      CALL alloc_NParray(mat_inv,(/ mole%nb_var,mole%nb_var /),"mat_inv",name_sub)


      CALL mat_id(mat,mole%nb_var,mole%nb_var)
      CALL mat_id(mat_inv,mole%nb_var,mole%nb_var)

      DO i_act=1,nb_NM
      DO k_act=1,nb_NM
        i_sym = mole%ActiveTransfo%list_QactTOQdyn(i_act)
        k_sym = mole%ActiveTransfo%list_QactTOQdyn(k_act)
        mat(i_sym,k_sym)     = d0c_inv(k_act,i_act)
        mat_inv(i_sym,k_sym) = d0c(k_act,i_act)
      END DO
      END DO

      IF (print_level > 0) THEN
        write(out_unitp,*) '========================================='
        write(out_unitp,*) '======== Basis parameters ==============='
        write(out_unitp,*) '========================================='
        write(out_unitp,*) '==========================='
        write(out_unitp,*) 'Parameters for the uncoupled HO basis'
        write(out_unitp,*) ' Qdyn0 =',mole%ActiveTransfo%Qdyn0
        write(out_unitp,*) '==========================='
        write(out_unitp,*) ' Hm basis set with Qdyn'
        DO i=1,nb_NM
          i_sym = mole%ActiveTransfo%list_QactTOQdyn(i)
          write(out_unitp,*) '&basis_nD iQdyn(1)=',i_sym,' name="Hm" nq=5 nb=5 Q0=', &
            mole%ActiveTransfo%Qdyn0(i_sym),' scaleQ=',ScalePara(i),' /'
        END DO
        write(out_unitp,*) '==========================='
      END IF


      mole%ActiveTransfo%Qdyn0(:) = matmul(mat_inv,mole%ActiveTransfo%Qdyn0(:))

      CALL Qdyn_TO_Qact_FROM_ActiveTransfo(mole%ActiveTransfo%Qdyn0,    &
                                           mole%ActiveTransfo%Qact0,    &
                                           mole%ActiveTransfo)

      Qact  = mole%ActiveTransfo%Qact0(:)

      IF (print_level > 0) THEN
        write(out_unitp,*) '==========================='
        write(out_unitp,*) 'Parameters for the HO basis (Normal Modes)'
        write(out_unitp,*) ' New Qdyn0, mole with the normal modes'
        write(out_unitp,*) ' Qdyn0 =',mole%ActiveTransfo%Qdyn0
        write(out_unitp,*) '==========================='
        write(out_unitp,*) ' Hm basis set with New Qdyn'
        DO i=1,nb_NM
          i_sym = mole%ActiveTransfo%list_QactTOQdyn(i)
          IF (mole%NMTransfo%k_Half) THEN
            write(out_unitp,*) '&basis_nD iQdyn(1)=',i_sym,' name="Hm" nq=5 nb=5 Q0=', &
              mole%ActiveTransfo%Qdyn0(i_sym),' scaleQ=',sqrt(d0k_save(i,i)),'/'
          ELSE
            write(out_unitp,*) '&basis_nD iQdyn(1)=',i_sym,' name="Hm" nq=5 nb=5 Q0=', &
              mole%ActiveTransfo%Qdyn0(i_sym),' scaleQ=1. /'
          END IF
        END DO
        write(out_unitp,*) '==========================='
      END IF

      IF (allocated(d0k_save))  THEN
        CALL dealloc_NParray(d0k_save,"d0k_save",name_sub)
      END IF

      mole%tab_Qtransfo(mole%itNM)%LinearTransfo%mat     = mat
      mole%tab_Qtransfo(mole%itNM)%LinearTransfo%mat_inv = mat_inv

      write(out_unitp,*) '========================================='
      write(out_unitp,*)

      CALL dealloc_NParray(mat,    "mat",    name_sub)
      CALL dealloc_NParray(mat_inv,"mat_inv",name_sub)

      write(out_unitp,*) '========================================='
      write(out_unitp,*) '======== New G matrix ==================='
      write(out_unitp,*) '========================================='

      DO i=1,mole%nb_var
        IF (mole%ActiveTransfo%list_act_OF_Qdyn(i) == -1)               &
                            mole%ActiveTransfo%list_act_OF_Qdyn(i) = 100
      END DO
      CALL type_var_analysis(mole) ! Qdyn0 => Qact0 is done in this subroutine
      ! but Qact has to be changed
      Qact(:) = mole%ActiveTransfo%Qact0

!     -----------------------------------------------------------------
!     - calc G_1
      CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,0)

      CALL get_dng_dnGG(Qact,para_Tnum,mole,dnGG=dnGG,nderiv=0)
      DO i=1,mole%nb_act
        d0eh(i) = dnGG%d0(i,i)
      END DO
      write(out_unitp,*) 'new d0GG',mole%nb_act
      CALL Write_Mat(dnGG%d0,out_unitp,5)

      CALL dealloc_dnSVM(dnGG)
!     -----------------------------------------------------------------
      write(out_unitp,*) '========================================='
      write(out_unitp,*)

      IF (print_level > 0) THEN
        write(out_unitp,*) '========================================='
        write(out_unitp,*) '======== Uncoupled Harmonic Hamiltonian ='
        write(out_unitp,*) '========================================='
        write(out_unitp,*) 'H(Qact(:)) = Sum_i H1D(Qact(i))'
        write(out_unitp,*) 'Active Freq (cm-1)',d0eh(1:mole%nb_act)*auTOcm_inv

        DO i=1,mole%nb_act
        write(out_unitp,*) 'H1D(Qact_i)',i,':',d0eh(i),'*1/2(-d./dQact_i^2+dQact_i^2)'
        END DO

        write(out_unitp,*) '========================================='
      END IF

      CALL dealloc_NParray(d0c_inv,  "d0c_inv",  name_sub)
      CALL dealloc_NParray(d0c_ini,  "d0c_ini",  name_sub)
      CALL dealloc_NParray(d0c,      "d0c",      name_sub)
      CALL dealloc_NParray(d0eh,     "d0eh",     name_sub)
      CALL dealloc_NParray(d0k,      "d0k",      name_sub)
      CALL dealloc_NParray(d0h,      "d0h",      name_sub)
      CALL dealloc_NParray(ScalePara,"ScalePara",name_sub)


!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*)
        !CALL Write_mole(mole)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF

!     -----------------------------------------------------------------

      END SUBROUTINE calc3_NM_TO_sym
      SUBROUTINE calc4_NM_TO_sym(Qact,mole,para_Tnum,para_PES,hCC,l_hCC)
      USE mod_system
      USE mod_dnSVM
      USE mod_Constant, only : get_Conv_au_TO_unit
      USE mod_Coord_KEO
      USE mod_SimpleOp
      USE mod_PrimOp_def
      IMPLICIT NONE

      TYPE (zmatrix) :: mole,mole_1
      TYPE (Tnum)    :: para_Tnum

      real (kind=Rkind), intent(inout) :: Qact(:)
      TYPE (param_PES) :: para_PES
      real (kind=Rkind), optional :: hCC(mole%ncart_act,mole%ncart_act)
      logical, optional           :: l_hCC  ! if .TRUE. hCC is already calculated (for PVSCF)


      TYPE(Type_dnMat) :: dnGG

      real (kind=Rkind), allocatable :: d0k(:,:),d0k_save(:,:),d0h(:,:)
      real (kind=Rkind), allocatable :: d0k_PerBlock(:,:),d0h_PerBlock(:,:)


      real (kind=Rkind), allocatable :: d0c_inv(:,:),d0c_ini(:,:)
      real (kind=Rkind), allocatable :: d0c(:,:),d0eh(:),d0eh_all(:)
      real (kind=Rkind), allocatable :: mat_inv(:,:),mat(:,:)
      real (kind=Rkind), allocatable :: ScalePara(:),ScalePara_NM(:)

      real (kind=Rkind) :: hCC_loc(mole%ncart_act,mole%ncart_act)
      logical           :: l_hCC_loc  ! if .TRUE. hCC is already calculated (for PVSCF)

      real (kind=Rkind) :: norm,max_freq


      integer :: k,i_act,i_sym,k_act,k_sym,ierr

      integer :: nb_NM
      integer, allocatable :: Ind_Coord_PerBlock(:)
      integer, allocatable :: nb_PerBlock(:),Ind_Coord_AtBlock(:)
      integer :: i_Block,nb_Block
      integer ::i,j,iQ,jQ,i2
      real (kind=Rkind) ::  auTOcm_inv

!      -----------------------------------------------------------------
      integer :: err_mem,memory
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'calc4_NM_TO_sym'
!      -----------------------------------------------------------------
       IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
        write(out_unitp,*) 'Qdyn0 =',mole%ActiveTransfo%Qdyn0(:)
        write(out_unitp,*)
!       CALL Write_mole(mole)
!       write(out_unitp,*)
       END IF

      !-----------------------------------------------------------------
      write(out_unitp,*) '========================================='
      write(out_unitp,*) '========== calc4_NM_TO_sym =============='
      write(out_unitp,*) '========================================='
      write(out_unitp,*)
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')

      IF (present(l_hCC) .AND. present(hCC)) THEN
        l_hCC_loc    = l_hCC
        IF (l_hCC_loc) hCC_loc(:,:) = hCC(:,:)
      END IF

      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      ! analysis the number of block ... => nb_NM
      CALL alloc_NParray(Ind_Coord_PerBlock,(/ mole%nb_var /),          &
                        "Ind_Coord_PerBlock",name_sub)

      ! parameters for uncoupled HO
      CALL alloc_NParray(ScalePara   ,(/mole%nb_var/),"ScalePara",name_sub)
      CALL alloc_NParray(ScalePara_NM,(/mole%nb_var/),"ScalePara_NM",name_sub)
      ScalePara(:)    = ONE
      ScalePara_NM(:) = ONE

      ! to store temporaly mat and mat_inv
      CALL alloc_NParray(mat,    (/ mole%nb_var,mole%nb_var /),"mat",name_sub)
      CALL alloc_NParray(mat_inv,(/ mole%nb_var,mole%nb_var /),"mat_inv",name_sub)

      CALL mat_id(mat,mole%nb_var,mole%nb_var)
      CALL mat_id(mat_inv,mole%nb_var,mole%nb_var)

      CALL alloc_NParray(d0eh_all,(/mole%nb_var/),"d0eh_all",name_sub)
      d0eh_all(:) = ZERO

      IF (associated(mole%NMTransfo%Qact1_sym)) THEN
        Ind_Coord_PerBlock(:) = mole%NMTransfo%Qact1_sym(:)

        ! first count the blocks
        nb_Block = 0
        DO
          i = maxval(Ind_Coord_PerBlock)
          IF ( i /= -Huge(1)) THEN
             nb_Block = nb_Block + 1
             WHERE (Ind_Coord_PerBlock == i) Ind_Coord_PerBlock = -Huge(1)
          ELSE
             EXIT
          END IF
        END DO
        Ind_Coord_PerBlock(:) = mole%NMTransfo%Qact1_sym(:)

        ! then, count the number of coordinates per block
        CALL alloc_NParray(nb_PerBlock,      (/ nb_Block /),"nb_PerBlock",      name_sub)
        CALL alloc_NParray(Ind_Coord_AtBlock,(/ nb_Block /),"Ind_Coord_AtBlock",name_sub)
        Ind_Coord_AtBlock(:) = 0
        nb_PerBlock(:)       = 0

        DO i_Block=1,nb_Block

          i = maxval(Ind_Coord_PerBlock)
          Ind_Coord_AtBlock(i_Block) = i

          IF (i /= 0) nb_PerBlock(i_Block) = count(Ind_Coord_PerBlock == i)

          WHERE (Ind_Coord_PerBlock == i) Ind_Coord_PerBlock = -Huge(1)

        END DO
        Ind_Coord_PerBlock(:) = mole%NMTransfo%Qact1_sym(:)

      ELSE
        Ind_Coord_PerBlock(:) = 1
        nb_Block    = 1
        CALL alloc_NParray(nb_PerBlock,      (/ nb_Block /),"nb_PerBlock",      name_sub)
        CALL alloc_NParray(Ind_Coord_AtBlock,(/ nb_Block /),"Ind_Coord_AtBlock",name_sub)

        Ind_Coord_AtBlock(1) = 1
        nb_PerBlock(1) = count(Ind_Coord_PerBlock == 1)
      END IF
      nb_NM = sum(nb_PerBlock)


      write(out_unitp,*) '  nb_Block            ',nb_Block
      write(out_unitp,*) '  nb_PerBlock(:)      ',nb_PerBlock(:)
      write(out_unitp,*) '  Ind_Coord_AtBlock(:)',Ind_Coord_AtBlock(:)
      write(out_unitp,*) '  nb_NM',nb_NM
      CALL flush_perso(out_unitp)
      !-----------------------------------------------------------------
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      DO i_Block=1,nb_Block

        IF (Ind_Coord_AtBlock(i_Block) == 0) CYCLE

        nb_NM = nb_PerBlock(i_Block)
        write(out_unitp,*) '========================================='
        write(out_unitp,*) '=========       Block: ',i_Block
        write(out_unitp,*) '========= nb_PerBlock: ',nb_PerBlock(i_Block)
        write(out_unitp,*) '========================================='
        write(out_unitp,*)
        CALL flush_perso(out_unitp)

        !- create mole_1 (type=-1 => type=1)
        CALL mole1TOmole2(mole,mole_1)
        ! a changer (utilisation de Qread_TO_Qact !!!
        DO i=1,mole_1%nb_var
          IF (Ind_Coord_PerBlock(i) == Ind_Coord_AtBlock(i_Block)) THEN
            mole_1%ActiveTransfo%list_act_OF_Qdyn(i) = 1
          ELSE IF (Ind_Coord_PerBlock(i) /= 0) THEN
            mole_1%ActiveTransfo%list_act_OF_Qdyn(i) = 100
          END IF
        END DO
        CALL type_var_analysis(mole_1)
        write(out_unitp,*) 'mole_1%...list_act_OF_Qdyn',                &
                                mole_1%ActiveTransfo%list_act_OF_Qdyn(:)

        CALL Qdyn_TO_Qact_FROM_ActiveTransfo(mole%ActiveTransfo%Qdyn0,  &
                                             mole%ActiveTransfo%Qact0,  &
                                             mole%ActiveTransfo)

        Qact = mole_1%ActiveTransfo%Qact0(:)
        IF (print_level > 1) write(out_unitp,*) 'Qact',Qact
        CALL flush_perso(out_unitp)


        IF (debug) THEN
          write(out_unitp,*) 'mole_1:'
          CALL Write_mole(mole_1)
          CALL flush_perso(out_unitp)
        END IF

        CALL alloc_NParray(d0c,     (/nb_NM,nb_NM/),"d0c",     name_sub)
        CALL alloc_NParray(d0c_inv, (/nb_NM,nb_NM/),"d0c_inv", name_sub)

        CALL alloc_NParray(d0k,     (/nb_NM,nb_NM/),"d0k",     name_sub)
        CALL alloc_NParray(d0h,     (/nb_NM,nb_NM/),"d0h",     name_sub)
        CALL alloc_NParray(d0k_save,(/nb_NM,nb_NM/),"d0k_save",name_sub)
        CALL alloc_NParray(d0c_ini, (/nb_NM,nb_NM/),"d0c_ini", name_sub)
        CALL alloc_NParray(d0eh,    (/ nb_NM /),    "d0eh",    name_sub)

        CALL get_hess_k(d0k,d0h,nb_NM,Qact,mole_1,para_Tnum,          &
                        Ind_Coord_AtBlock(i_Block),Ind_Coord_PerBlock,  &
                        para_PES,hCC_loc,l_hCC_loc)

        ! scaleQ for the uncoupled HO
        iQ = 0
        DO i=1,mole%nb_var
          IF (Ind_Coord_PerBlock(i) /= Ind_Coord_AtBlock(i_Block) .OR.  &
              abs(mole_1%ActiveTransfo%list_act_OF_Qdyn(i)) /= 1) CYCLE
          iQ = iQ + 1
          ScalePara(i) = sqrt(sqrt(abs(d0h(iQ,iQ)/d0k(iQ,iQ))))
        END DO


        d0k_save(:,:) = d0k(:,:)
        d0c_ini(:,:)  = ZERO
        CALL calc_freq(nb_NM,d0h,d0k_save,d0eh,d0c,d0c_inv,norm,d0c_ini,.FALSE.)

        DO i=1,mole_1%nb_act,3
          i2 = min(i+2,mole_1%nb_act)
          write(out_unitp,'("frequencies (cm-1): ",i0,"-",i0,3(x,f0.4))') &
                          i,i2,d0eh(i:i2)*auTOcm_inv
        END DO
        CALL flush_perso(out_unitp)

        ! frequencies
        iQ = 0
        DO i=1,mole%nb_var
          IF (Ind_Coord_PerBlock(i) == Ind_Coord_AtBlock(i_Block)) THEN
            iQ = iQ + 1
            d0eh_all(i) = d0eh(iQ)
          END IF
        END DO

        d0k_save = matmul(transpose(d0c),matmul(d0k,d0c))
        IF (print_level > 1 .OR. debug) THEN
          write(out_unitp,*) '==========================='
          write(out_unitp,*) 'new d0k'
          CALL Write_Mat(d0k_save,out_unitp,5)
          write(out_unitp,*) '==========================='
          CALL flush_perso(out_unitp)
        END IF

        ! change d0c, d0c_inv
        IF (mole%NMTransfo%k_Half) THEN
          DO i=1,nb_NM
            d0c(:,i)     = d0c(:,i)     / sqrt(d0k_save(i,i))
            d0c_inv(i,:) = d0c_inv(i,:) * sqrt(d0k_save(i,i))
          END DO

          ! scaleQ for the coupled HO (NM)
          iQ = 0
          DO i=1,mole%nb_var
            IF (Ind_Coord_PerBlock(i) /= Ind_Coord_AtBlock(i_Block) .OR.  &
              abs(mole_1%ActiveTransfo%list_act_OF_Qdyn(i)) /= 1) CYCLE
            iQ = iQ + 1
            ScalePara_NM(i) = sqrt(d0k_save(iQ,iQ))
            !write(6,*) 'i,iQ,d0h,d0k',i,iQ,d0h(iQ,iQ),d0k(iQ,iQ)
            !write(6,*) 'i,iQ,ScalePara_NM(i)',i,iQ,ScalePara_NM(i)
          END DO
        END IF


        IF (debug) THEN
          !write with high precision to be able to read it
          write(out_unitp,*) 'd0c'
          write(out_unitp,*) nb_NM,5
          CALL Write_Mat(d0c,out_unitp,5,Rformat='e20.13')
        END IF

        IF (Ind_Coord_AtBlock(i_block) > 0) THEN
        iQ = 0
        DO i=1,mole%nb_var
          IF (Ind_Coord_PerBlock(i) /= Ind_Coord_AtBlock(i_Block) .OR.  &
              abs(mole_1%ActiveTransfo%list_act_OF_Qdyn(i)) /= 1) CYCLE
          iQ = iQ + 1

          jQ = 0
          DO j=1,mole%nb_var
            IF (Ind_Coord_PerBlock(j) /= Ind_Coord_AtBlock(i_Block) .OR.&
              abs(mole_1%ActiveTransfo%list_act_OF_Qdyn(j)) /= 1) CYCLE
            jQ = jQ + 1
            mat(i,j)     = d0c_inv(jQ,iQ)
            mat_inv(i,j) = d0c(jQ,iQ)
          END DO

        END DO
        END IF
        mole%ActiveTransfo%Qdyn0 = mole_1%ActiveTransfo%Qdyn0

        CALL dealloc_NParray(d0k_save,"d0k_save",name_sub)
        CALL dealloc_NParray(d0c_ini, "d0c_ini", name_sub)
        CALL dealloc_NParray(d0eh,    "d0eh",    name_sub)
        CALL dealloc_NParray(d0k,     "d0k",     name_sub)
        CALL dealloc_NParray(d0h,     "d0h",     name_sub)
        CALL dealloc_NParray(d0c,     "d0c",     name_sub)
        CALL dealloc_NParray(d0c_inv, "d0c_inv", name_sub)

        CALL dealloc_zmat(mole_1)


        write(out_unitp,*) '========================================='
        write(out_unitp,*) '===== END Block: ',i_Block
        write(out_unitp,*) '========================================='
        CALL flush_perso(out_unitp)

      END DO

      DO i=1,mole%nb_act,3
        i2 = min(i+2,mole%nb_act)
        write(out_unitp,'("frequencies (cm-1): ",i0,"-",i0,3(x,f0.4))') &
                          i,i2,d0eh_all(i:i2)*auTOcm_inv
      END DO
      CALL alloc_array(mole%NMTransfo%d0eh,(/ nb_NM /),                 &
                      "mole%NMTransfo%d0eh",name_sub)

      mole%NMTransfo%d0eh(1:nb_NM) = d0eh_all(1:nb_NM)
      mole%NMTransfo%nb_NM         = nb_NM

      !write(out_unitp,*) ' Mat'
      !CALL Write_Mat(mat,out_unitp,5)
      !write(out_unitp,*) ' Mat_inv'
      !CALL Write_Mat(mat_inv,out_unitp,5)

      IF (print_level > 1 .OR. debug) THEN
        write(out_unitp,*) '==========================='
        write(out_unitp,*) 'Parameters for the uncoupled HO basis'
        write(out_unitp,*) '==========================='
        write(out_unitp,*) ' HO basis set with Qdyn'
        DO i=1,mole%nb_var
          IF (Ind_Coord_PerBlock(i) /= 0) THEN
            write(out_unitp,'(a,i0,a,f0.3,a,f0.3,a)')                   &
                 '  &basis_nD iQdyn(1)= ',i,' name="Hm" nq=5 nb=5 Q0=', &
                mole%ActiveTransfo%Qdyn0(i),' scaleQ=',ScalePara(i),' /'
          END IF
        END DO
        write(out_unitp,*) '==========================='
      END IF


      mole%ActiveTransfo%Qdyn0(:) = matmul(mat_inv,mole%ActiveTransfo%Qdyn0(:))

      CALL Qdyn_TO_Qact_FROM_ActiveTransfo(mole%ActiveTransfo%Qdyn0,    &
                                           mole%ActiveTransfo%Qact0,    &
                                           mole%ActiveTransfo)



      CALL alloc_array(mole%NMTransfo%Q0_HObasis,(/mole%nb_var/),       &
                      "mole%NMTransfo%Q0_HObasis",name_sub)
      mole%NMTransfo%Q0_HObasis(:) = mole%ActiveTransfo%Qdyn0(:)

      CALL alloc_array(mole%NMTransfo%scaleQ_HObasis,(/mole%nb_var/),   &
                      "mole%NMTransfo%scaleQ_HObasis",name_sub)
      mole%NMTransfo%scaleQ_HObasis(:) = ScalePara_NM(:)

      write(out_unitp,*) '==========================='
      write(out_unitp,*) 'New Qact0',mole%ActiveTransfo%Qact0(:)
      write(out_unitp,*) '==========================='

      IF (print_level > 1 .OR. debug) THEN

        write(out_unitp,*) '==========================='
        write(out_unitp,*) 'Parameters for the HO basis (Normal Modes)'

        write(out_unitp,*) '==========================='
        write(out_unitp,*) ' Hm basis set with New Qdyn'

        DO i=1,mole%nb_var
          IF (Ind_Coord_PerBlock(i) /= 0) THEN
            write(out_unitp,'(a,i0,a,f0.3,a,e11.4,a)')                &
               '  &basis_nD iQdyn(1)= ',i,' name="Hm" nq=5 nb=5 Q0=', &
             mole%ActiveTransfo%Qdyn0(i),' scaleQ=',ScalePara_NM(i),' /'
          END IF
        END DO
        write(out_unitp,*) '==========================='
      END IF

      IF (print_level > 1 .OR. debug) THEN
        write(out_unitp,*) '=================================='
        write(out_unitp,*) '=== Mat of Linear transformation ='
        write(out_unitp,*) '=================================='
        write(out_unitp,*) " &Coord_transfo name_transfo='linear' check_LinearTransfo=f /"
        write(out_unitp,*) 'Mat of linear Transfo for the Normal modes'
        write(out_unitp,*) 5
        CALL Write_Mat(mat,out_unitp,5,Rformat='e20.13')
        write(out_unitp,*) '=================================='
        write(out_unitp,*) '=================================='
      END IF


      IF (print_level > 2 .OR. debug) THEN
        write(out_unitp,*) '=================================='
        write(out_unitp,*) '=================================='
        write(out_unitp,*) 'transpose(mat_inv) or d0c'
        write(out_unitp,*) 'Each column corresponds to one normal mode'
        CALL Write_Mat(transpose(mat_inv),out_unitp,5,Rformat='f10.3')
        write(out_unitp,*) '=================================='
        write(out_unitp,*) '=================================='
        write(out_unitp,*) '=================================='
        write(out_unitp,*) '=================================='
        write(out_unitp,*) 'transpose(mat) or d0c_inv'
        write(out_unitp,*) 'Each line corresponds to one normal mode'
        CALL Write_Mat(transpose(mat),out_unitp,5,Rformat='f14.7')
        write(out_unitp,*) '=================================='
        write(out_unitp,*) '=================================='
      END IF

      mole%tab_Qtransfo(mole%itNM)%LinearTransfo%mat(:,:)     = mat(:,:)
      mole%tab_Qtransfo(mole%itNM)%LinearTransfo%mat_inv(:,:) = mat_inv(:,:)


      CALL dealloc_NParray(mat,    "mat",name_sub)
      CALL dealloc_NParray(mat_inv,"mat_inv",name_sub)

      CALL dealloc_NParray(nb_PerBlock,"nb_PerBlock",name_sub)
      CALL dealloc_NParray(Ind_Coord_AtBlock,"Ind_Coord_AtBlock",name_sub)
      CALL dealloc_NParray(ScalePara,"ScalePara",name_sub)
      CALL dealloc_NParray(ScalePara_NM,"ScalePara_NM",name_sub)
      CALL dealloc_NParray(d0eh_all,"d0eh_all",name_sub)
      CALL dealloc_NParray(Ind_Coord_PerBlock,"Ind_Coord_PerBlock",name_sub)

      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      DO i=1,mole%nb_var
        IF (mole%ActiveTransfo%list_act_OF_Qdyn(i) == -1)               &
                            mole%ActiveTransfo%list_act_OF_Qdyn(i) = 100
      END DO
      CALL type_var_analysis(mole) ! Qdyn0 => Qact0 is done in this subroutine
      ! but Qact has to be changed

      CALL get_Qact0(Qact,mole%ActiveTransfo)

      !-----------------------------------------------------------------
      !- calc G_1
      IF (print_level > 1 .OR. debug) THEN
        write(out_unitp,*) '========================================='
        write(out_unitp,*) '======== New G matrix ==================='
        write(out_unitp,*) '========================================='

        CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,0)

        CALL get_dng_dnGG(Qact,para_Tnum,mole,dnGG=dnGG,nderiv=0)

        write(out_unitp,*) 'New d0GG',mole%nb_act
        CALL Write_Mat(dnGG%d0,out_unitp,5)
        write(out_unitp,*) '========================================='
        write(out_unitp,*) '========================================='

        CALL dealloc_dnSVM(dnGG)
      END IF
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*)
        !CALL Write_mole(mole)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF

!     -----------------------------------------------------------------

      END SUBROUTINE calc4_NM_TO_sym
      SUBROUTINE calc5_NM_TO_sym(Qact,mole,para_Tnum,para_PES,hCC,l_hCC)
      USE mod_system
      USE mod_dnSVM
      USE mod_Constant, only : get_Conv_au_TO_unit
      USE mod_Coord_KEO
      USE mod_SimpleOp
      USE mod_PrimOp_def
      IMPLICIT NONE

      TYPE (zmatrix) :: mole,mole_1
      TYPE (Tnum)    :: para_Tnum

      real (kind=Rkind), intent(inout) :: Qact(:)
      TYPE (param_PES) :: para_PES
      real (kind=Rkind), optional :: hCC(mole%ncart_act,mole%ncart_act)
      logical, optional           :: l_hCC  ! if .TRUE. hCC is already calculated (for PVSCF)


      TYPE(Type_dnMat) :: dnGG

      real (kind=Rkind), allocatable :: d0k(:,:),d0k_save(:,:),d0h(:,:)
      real (kind=Rkind), allocatable :: d0k_PerBlock(:,:),d0h_PerBlock(:,:)


      real (kind=Rkind), allocatable :: d0c_inv(:,:),d0c_ini(:,:)
      real (kind=Rkind), allocatable :: d0c(:,:),d0eh(:),d0eh_all(:)
      real (kind=Rkind), allocatable :: mat_inv(:,:),mat(:,:)
      real (kind=Rkind), allocatable :: ScalePara(:),ScalePara_NM(:)

      real (kind=Rkind) :: hCC_loc(mole%ncart_act,mole%ncart_act)
      logical           :: l_hCC_loc  ! if .TRUE. hCC is already calculated (for PVSCF)

      real (kind=Rkind) :: norm,max_freq


      integer :: k,i_act,i_sym,k_act,k_sym,ierr

      integer :: nb_NM
      integer, allocatable :: Ind_Coord_PerBlock(:)
      integer, allocatable :: nb_PerBlock(:),Ind_Coord_AtBlock(:)
      integer :: i_Block,nb_Block
      integer ::i,j,iQ,jQ,i2
      real (kind=Rkind) ::  auTOcm_inv

!      -----------------------------------------------------------------
      integer :: err_mem,memory
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'calc5_NM_TO_sym'
!      -----------------------------------------------------------------
       IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
        write(out_unitp,*) 'Qdyn0 =',mole%ActiveTransfo%Qdyn0(:)
        write(out_unitp,*)
!       CALL Write_mole(mole)
!       write(out_unitp,*)
       END IF

      !-----------------------------------------------------------------
      write(out_unitp,*) '========================================='
      write(out_unitp,*) '========== calc5_NM_TO_sym =============='
      write(out_unitp,*) '========================================='
      write(out_unitp,*)
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')

      IF (present(l_hCC) .AND. present(hCC)) THEN
        l_hCC_loc    = l_hCC
        IF (l_hCC_loc) hCC_loc(:,:) = hCC(:,:)
      END IF

      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      ! analysis the number of block ... => nb_NM
      CALL alloc_NParray(Ind_Coord_PerBlock,(/ mole%nb_var /),          &
                        "Ind_Coord_PerBlock",name_sub)

      ! parameters for uncoupled HO
      CALL alloc_NParray(ScalePara   ,(/mole%nb_var/),"ScalePara",name_sub)
      CALL alloc_NParray(ScalePara_NM,(/mole%nb_var/),"ScalePara_NM",name_sub)
      ScalePara(:)    = ONE
      ScalePara_NM(:) = ONE

      ! to store temporaly mat and mat_inv
      CALL alloc_NParray(mat,    (/ mole%nb_var,mole%nb_var /),"mat",name_sub)
      CALL alloc_NParray(mat_inv,(/ mole%nb_var,mole%nb_var /),"mat_inv",name_sub)

      CALL mat_id(mat,mole%nb_var,mole%nb_var)
      CALL mat_id(mat_inv,mole%nb_var,mole%nb_var)

      CALL alloc_NParray(d0eh_all,(/mole%nb_var/),"d0eh_all",name_sub)
      d0eh_all(:) = ZERO

      IF (.NOT. associated(mole%NMTransfo%Qact1_sym)) THEN
        CALL alloc_array(mole%NMTransfo%Qact1_sym,(/mole%nb_var/),      &
                        "mole%NMTransfo%Qact1_sym",name_sub)
        mole%NMTransfo%Qact1_sym(:) = mole%ActiveTransfo%list_act_OF_Qdyn(:)
      END IF

      IF (associated(mole%NMTransfo%Qact1_sym)) THEN
        Ind_Coord_PerBlock(:) = mole%NMTransfo%Qact1_sym(:)

        ! first count the blocks
        nb_Block = 0
        DO
          i = maxval(Ind_Coord_PerBlock)
          IF ( i /= -Huge(1)) THEN
             nb_Block = nb_Block + 1
             WHERE (Ind_Coord_PerBlock == i) Ind_Coord_PerBlock = -Huge(1)
          ELSE
             EXIT
          END IF
        END DO
        Ind_Coord_PerBlock(:) = mole%NMTransfo%Qact1_sym(:)

        ! then, count the number of coordinates per block
        CALL alloc_NParray(nb_PerBlock,      (/ nb_Block /),"nb_PerBlock",      name_sub)
        CALL alloc_NParray(Ind_Coord_AtBlock,(/ nb_Block /),"Ind_Coord_AtBlock",name_sub)
        Ind_Coord_AtBlock(:) = 0
        nb_PerBlock(:)       = 0

        DO i_Block=1,nb_Block

          i = maxval(Ind_Coord_PerBlock)
          Ind_Coord_AtBlock(i_Block) = i

          nb_PerBlock(i_Block) = count(Ind_Coord_PerBlock == i)

          WHERE (Ind_Coord_PerBlock == i) Ind_Coord_PerBlock = -Huge(1)

        END DO
        Ind_Coord_PerBlock(:) = mole%NMTransfo%Qact1_sym(:)

      ELSE
        Ind_Coord_PerBlock(:) = 1
        nb_Block    = 1
        CALL alloc_NParray(nb_PerBlock,      (/ nb_Block /),"nb_PerBlock",      name_sub)
        CALL alloc_NParray(Ind_Coord_AtBlock,(/ nb_Block /),"Ind_Coord_AtBlock",name_sub)

        Ind_Coord_AtBlock(1) = 1
        nb_PerBlock(1) = count(Ind_Coord_PerBlock == 1)
      END IF
      nb_NM = sum(nb_PerBlock)


      write(out_unitp,*) '  nb_Block            ',nb_Block
      write(out_unitp,*) '  nb_PerBlock(:)      ',nb_PerBlock(:)
      write(out_unitp,*) '  Ind_Coord_AtBlock(:)',Ind_Coord_AtBlock(:)
      write(out_unitp,*) '  nb_NM',nb_NM

      !-----------------------------------------------------------------
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      DO i_Block=1,nb_Block

        nb_NM = nb_PerBlock(i_Block)
        write(out_unitp,*) '========================================='
        write(out_unitp,*) '=========       Block: ',i_Block
        write(out_unitp,*) '========= nb_PerBlock: ',nb_PerBlock(i_Block)
        write(out_unitp,*) '========= Ind_AtBlock: ',Ind_Coord_AtBlock(i_Block)
        write(out_unitp,*) '========================================='
        write(out_unitp,*)

        IF (Ind_Coord_AtBlock(i_Block) == 0 .OR. Ind_Coord_AtBlock(i_Block) == 100) THEN

        ELSE
          !- create mole_1 (type=-1 => type=1)
          CALL mole1TOmole2(mole,mole_1)
          ! a changer (utilisation de Qread_TO_Qact !!!
          DO i=1,mole_1%nb_var
            IF (Ind_Coord_PerBlock(i) == Ind_Coord_AtBlock(i_Block)) THEN
              mole_1%ActiveTransfo%list_act_OF_Qdyn(i) = 1
            ELSE IF (Ind_Coord_PerBlock(i) /= 0) THEN
              mole_1%ActiveTransfo%list_act_OF_Qdyn(i) = 100
            END IF
          END DO
          CALL type_var_analysis(mole_1)
          write(out_unitp,*) 'mole_1%...list_act_OF_Qdyn',                &
                                  mole_1%ActiveTransfo%list_act_OF_Qdyn(:)

          CALL Qdyn_TO_Qact_FROM_ActiveTransfo(mole%ActiveTransfo%Qdyn0,  &
                                               mole%ActiveTransfo%Qact0,  &
                                               mole%ActiveTransfo)

          Qact = mole_1%ActiveTransfo%Qact0(:)
          IF (print_level > 1) write(out_unitp,*) 'Qact',Qact


          IF (debug) THEN
            write(out_unitp,*) 'mole_1:'
            CALL Write_mole(mole_1)
          END IF

          CALL alloc_NParray(d0c,     (/nb_NM,nb_NM/),"d0c",     name_sub)
          CALL alloc_NParray(d0c_inv, (/nb_NM,nb_NM/),"d0c_inv", name_sub)

          CALL alloc_NParray(d0k,     (/nb_NM,nb_NM/),"d0k",     name_sub)
          CALL alloc_NParray(d0h,     (/nb_NM,nb_NM/),"d0h",     name_sub)
          CALL alloc_NParray(d0k_save,(/nb_NM,nb_NM/),"d0k_save",name_sub)
          CALL alloc_NParray(d0c_ini, (/nb_NM,nb_NM/),"d0c_ini", name_sub)
          CALL alloc_NParray(d0eh,    (/ nb_NM /),    "d0eh",    name_sub)


          CALL get_hess_k(d0k,d0h,nb_NM,Qact,mole_1,para_Tnum,          &
                          Ind_Coord_AtBlock(i_Block),Ind_Coord_PerBlock,  &
                          para_PES,hCC_loc,l_hCC_loc)

          ! scaleQ for the uncoupled HO
          iQ = 0
          DO i=1,mole%nb_var
            IF (Ind_Coord_PerBlock(i) /= Ind_Coord_AtBlock(i_Block) .OR.  &
                abs(mole_1%ActiveTransfo%list_act_OF_Qdyn(i)) /= 1) CYCLE
            iQ = iQ + 1
            ScalePara(i) = sqrt(sqrt(abs(d0h(iQ,iQ)/d0k(iQ,iQ))))
          END DO


          d0k_save(:,:) = d0k(:,:)
          d0c_ini(:,:)  = ZERO
          CALL calc_freq(nb_NM,d0h,d0k_save,d0eh,d0c,d0c_inv,norm,d0c_ini,.FALSE.)

          DO i=1,mole_1%nb_act,3
            i2 = min(i+2,mole_1%nb_act)
            write(out_unitp,'("frequencies (cm-1): ",i0,"-",i0,3(x,f0.4))') &
                            i,i2,d0eh(i:i2)*auTOcm_inv
          END DO

          ! frequencies
          iQ = 0
          DO i=1,mole%nb_var
            IF (Ind_Coord_PerBlock(i) == Ind_Coord_AtBlock(i_Block)) THEN
              iQ = iQ + 1
              d0eh_all(i) = d0eh(iQ)
            END IF
          END DO

          d0k_save = matmul(transpose(d0c),matmul(d0k,d0c))
          IF (print_level > 1 .OR. debug) THEN
            write(out_unitp,*) '==========================='
            write(out_unitp,*) 'new d0k'
            CALL Write_Mat(d0k_save,out_unitp,5)
            write(out_unitp,*) '==========================='
          END IF

          ! change d0c, d0c_inv
          IF (mole%NMTransfo%k_Half) THEN
            DO i=1,nb_NM
              d0c(:,i)     = d0c(:,i)     / sqrt(d0k_save(i,i))
              d0c_inv(i,:) = d0c_inv(i,:) * sqrt(d0k_save(i,i))
            END DO

            ! scaleQ for the coupled HO (NM)
            iQ = 0
            DO i=1,mole%nb_var
              IF (Ind_Coord_PerBlock(i) /= Ind_Coord_AtBlock(i_Block) .OR.  &
                abs(mole_1%ActiveTransfo%list_act_OF_Qdyn(i)) /= 1) CYCLE
              iQ = iQ + 1
              ScalePara_NM(i) = sqrt(d0k_save(iQ,iQ))
              !write(6,*) 'i,iQ,d0h,d0k',i,iQ,d0h(iQ,iQ),d0k(iQ,iQ)
              !write(6,*) 'i,iQ,ScalePara_NM(i)',i,iQ,ScalePara_NM(i)
            END DO
          END IF


          IF (debug) THEN
            !write with high precision to be able to read it
            write(out_unitp,*) 'd0c'
            write(out_unitp,*) nb_NM,5
            CALL Write_Mat(d0c,out_unitp,5,Rformat='e20.13')
          END IF

          IF (Ind_Coord_AtBlock(i_block) > 0) THEN
            iQ = 0
            DO i=1,mole%nb_var
              IF (Ind_Coord_PerBlock(i) /= Ind_Coord_AtBlock(i_Block) .OR.  &
                  abs(mole_1%ActiveTransfo%list_act_OF_Qdyn(i)) /= 1) CYCLE
              iQ = iQ + 1

              jQ = 0
              DO j=1,mole%nb_var
                IF (Ind_Coord_PerBlock(j) /= Ind_Coord_AtBlock(i_Block) .OR.&
                  abs(mole_1%ActiveTransfo%list_act_OF_Qdyn(j)) /= 1) CYCLE
                jQ = jQ + 1
                mat(i,j)     = d0c_inv(jQ,iQ)
                mat_inv(i,j) = d0c(jQ,iQ)
              END DO

            END DO
          END IF
          mole%ActiveTransfo%Qdyn0 = mole_1%ActiveTransfo%Qdyn0

          CALL dealloc_NParray(d0k_save,"d0k_save",name_sub)
          CALL dealloc_NParray(d0c_ini, "d0c_ini", name_sub)
          CALL dealloc_NParray(d0eh,    "d0eh",    name_sub)
          CALL dealloc_NParray(d0k,     "d0k",     name_sub)
          CALL dealloc_NParray(d0h,     "d0h",     name_sub)
          CALL dealloc_NParray(d0c,     "d0c",     name_sub)
          CALL dealloc_NParray(d0c_inv, "d0c_inv", name_sub)

          CALL dealloc_zmat(mole_1)
        END IF

        write(out_unitp,*) '========================================='
        write(out_unitp,*) '===== END Block: ',i_Block
        write(out_unitp,*) '========================================='
        CALL flush_perso(out_unitp)

      END DO

      DO i=1,mole%nb_act,3
        i2 = min(i+2,mole%nb_act)
        write(out_unitp,'("frequencies (cm-1): ",i0,"-",i0,3(x,f0.4))') &
                          i,i2,d0eh_all(i:i2)*auTOcm_inv
      END DO
      CALL alloc_array(mole%NMTransfo%d0eh,(/ nb_NM /),                 &
                      "mole%NMTransfo%d0eh",name_sub)

      mole%NMTransfo%d0eh(:)      = d0eh_all(:)
      mole%NMTransfo%nb_NM        = nb_NM

      !write(out_unitp,*) ' Mat'
      !CALL Write_Mat(mat,out_unitp,5)
      !write(out_unitp,*) ' Mat_inv'
      !CALL Write_Mat(mat_inv,out_unitp,5)

      IF (print_level > 1 .OR. debug) THEN
        write(out_unitp,*) '==========================='
        write(out_unitp,*) 'Parameters for the uncoupled HO basis'
        write(out_unitp,*) '==========================='
        write(out_unitp,*) ' HO basis set with Qdyn'
        DO i=1,mole%nb_var
          IF (Ind_Coord_PerBlock(i) /= 0) THEN
            write(out_unitp,'(a,i0,a,f0.3,a,f0.3,a)')                   &
                 '  &basis_nD iQdyn(1)= ',i,' name="Hm" nq=5 nb=5 Q0=', &
                mole%ActiveTransfo%Qdyn0(i),' scaleQ=',ScalePara(i),' /'
          END IF
        END DO
        write(out_unitp,*) '==========================='
      END IF


      mole%ActiveTransfo%Qdyn0(:) = matmul(mat_inv,mole%ActiveTransfo%Qdyn0(:))

      CALL Qdyn_TO_Qact_FROM_ActiveTransfo(mole%ActiveTransfo%Qdyn0,    &
                                           mole%ActiveTransfo%Qact0,    &
                                           mole%ActiveTransfo)



      CALL alloc_array(mole%NMTransfo%Q0_HObasis,(/mole%nb_var/),       &
                      "mole%NMTransfo%Q0_HObasis",name_sub)
      mole%NMTransfo%Q0_HObasis(:) = mole%ActiveTransfo%Qdyn0(:)

      CALL alloc_array(mole%NMTransfo%scaleQ_HObasis,(/mole%nb_var/),   &
                      "mole%NMTransfo%scaleQ_HObasis",name_sub)
      mole%NMTransfo%scaleQ_HObasis(:) = ScalePara_NM(:)

      write(out_unitp,*) '==========================='
      write(out_unitp,*) 'New Qact0',mole%ActiveTransfo%Qact0(:)
      write(out_unitp,*) '==========================='

      IF (print_level > 1 .OR. debug) THEN

        write(out_unitp,*) '==========================='
        write(out_unitp,*) 'Parameters for the HO basis (Normal Modes)'

        write(out_unitp,*) '==========================='
        write(out_unitp,*) ' Hm basis set with New Qdyn'

        DO i=1,mole%nb_var
          IF (Ind_Coord_PerBlock(i) /= 0) THEN
            write(out_unitp,'(a,i0,a,f0.3,a,e11.4,a)')                &
               '  &basis_nD iQdyn(1)= ',i,' name="Hm" nq=5 nb=5 Q0=', &
             mole%ActiveTransfo%Qdyn0(i),' scaleQ=',ScalePara_NM(i),' /'
          END IF
        END DO
        write(out_unitp,*) '==========================='
      END IF

      IF (print_level > 1 .OR. debug) THEN
        write(out_unitp,*) '=================================='
        write(out_unitp,*) '=== Mat of Linear transformation ='
        write(out_unitp,*) '=================================='
        write(out_unitp,*) " &Coord_transfo name_transfo='linear' check_LinearTransfo=f /"
        write(out_unitp,*) 'Mat of linear Transfo for the Normal modes'
        write(out_unitp,*) 5
        CALL Write_Mat(mat,out_unitp,5,Rformat='e20.13')
        write(out_unitp,*) '=================================='
        write(out_unitp,*) '=================================='
      END IF


      IF (print_level > 2 .OR. debug) THEN
        write(out_unitp,*) '=================================='
        write(out_unitp,*) '=================================='
        write(out_unitp,*) 'transpose(mat_inv) or d0c'
        write(out_unitp,*) 'Each column corresponds to one normal mode'
        CALL Write_Mat(transpose(mat_inv),out_unitp,5,Rformat='f10.3')
        write(out_unitp,*) '=================================='
        write(out_unitp,*) '=================================='
      END IF

      mole%tab_Qtransfo(mole%itNM)%LinearTransfo%mat(:,:)     = mat(:,:)
      mole%tab_Qtransfo(mole%itNM)%LinearTransfo%mat_inv(:,:) = mat_inv(:,:)


      CALL dealloc_NParray(mat,    "mat",name_sub)
      CALL dealloc_NParray(mat_inv,"mat_inv",name_sub)

      CALL dealloc_NParray(nb_PerBlock,"nb_PerBlock",name_sub)
      CALL dealloc_NParray(Ind_Coord_AtBlock,"Ind_Coord_AtBlock",name_sub)
      CALL dealloc_NParray(ScalePara,"ScalePara",name_sub)
      CALL dealloc_NParray(ScalePara_NM,"ScalePara_NM",name_sub)
      CALL dealloc_NParray(d0eh_all,"d0eh_all",name_sub)
      CALL dealloc_NParray(Ind_Coord_PerBlock,"Ind_Coord_PerBlock",name_sub)

      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      DO i=1,mole%nb_var
        IF (mole%ActiveTransfo%list_act_OF_Qdyn(i) == -1)               &
                            mole%ActiveTransfo%list_act_OF_Qdyn(i) = 100
      END DO
      CALL type_var_analysis(mole) ! Qdyn0 => Qact0 is done in this subroutine
      ! but Qact has to be changed

      CALL get_Qact0(Qact,mole%ActiveTransfo)

      !-----------------------------------------------------------------
      !- calc G_1
      IF (print_level > 1 .OR. debug) THEN
        write(out_unitp,*) '========================================='
        write(out_unitp,*) '======== New G matrix ==================='
        write(out_unitp,*) '========================================='

        CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,0)

        CALL get_dng_dnGG(Qact,para_Tnum,mole,dnGG=dnGG,nderiv=0)

        write(out_unitp,*) 'New d0GG',mole%nb_act
        CALL Write_Mat(dnGG%d0,out_unitp,5)
        write(out_unitp,*) '========================================='
        write(out_unitp,*) '========================================='

        CALL dealloc_dnSVM(dnGG)
      END IF
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*)
        !CALL Write_mole(mole)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF

!     -----------------------------------------------------------------

      END SUBROUTINE calc5_NM_TO_sym
      SUBROUTINE get_hess_k(d0k,d0h,nb_NM,Qact,mole,para_Tnum,          &
                            Ind_Coord_AtBlock,Ind_Coord_PerBlock,       &
                            para_PES,hCC,l_hCC)
      USE mod_system
      USE mod_dnSVM
      USE mod_Coord_KEO, only : zmatrix,Tnum,get_dng_dnGG, sub_dnFCC_TO_dnFcurvi

      USE mod_SimpleOp
      USE mod_PrimOp_def
      IMPLICIT NONE

      integer           :: nb_NM
      real (kind=Rkind) :: d0k(nb_NM,nb_NM),d0h(nb_NM,nb_NM)


      TYPE (zmatrix)     :: mole
      TYPE (Tnum)        :: para_Tnum

      real (kind=Rkind), intent(inout) :: Qact(:)
      TYPE (param_PES)   :: para_PES
      real (kind=Rkind)  :: hCC(mole%ncart_act,mole%ncart_act)
      logical            :: l_hCC  ! if .TRUE. hCC is already calculated (for PVSCF)
      integer :: Ind_Coord_AtBlock,Ind_Coord_PerBlock(mole%nb_var)
      real (kind=Rkind) :: d0grad(nb_NM)


      TYPE(Type_dnMat) :: dnGG

      TYPE(Type_dnS)       :: dnECC(1,1),dnE(1,1)
      TYPE (param_dnMatOp) :: dnMatOp(1)

      real (kind=Rkind) :: d0eh(nb_NM),d0ch(nb_NM,nb_NM),d0hh(nb_NM,nb_NM)


      integer :: i,j,ib,jb,k,i_act,i_sym,k_act,k_sym,ierr

      logical                  :: Read_OnTheFly_only,OnTheFly,calc_scalar_Op
      integer                  :: nb_scalar_Op

      character (len=Line_len) :: name_FChk


      !-----------------------------------------------------------------
      integer :: err_mem,memory
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'get_hess_k'
      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
        write(out_unitp,*) 'hessian_read,k_read',mole%NMTransfo%hessian_read,mole%NMTransfo%k_read
        write(out_unitp,*) 'hessian_old',mole%NMTransfo%hessian_old
        write(out_unitp,*) 'hessian_onthefly',mole%NMTransfo%hessian_onthefly
        write(out_unitp,*) 'hessian_cart',mole%NMTransfo%hessian_cart
        write(out_unitp,*)
        !CALL Write_mole(mole)
        !write(out_unitp,*)
        CALL flush_perso(out_unitp)
      END IF
      !-----------------------------------------------------------------

      Qact = mole%ActiveTransfo%Qact0(:)
      IF (debug) write(out_unitp,*) 'Qact =',Qact

      IF (print_level > 1 .OR. debug) THEN
        write(out_unitp,*) '========================================='
        write(out_unitp,*) '========================================='
        write(out_unitp,*) '========================================='
        write(out_unitp,*) '==== hessian and kinetic matrices ======='
        write(out_unitp,*) '========================================='
        CALL flush_perso(out_unitp)
      END IF

      IF (mole%NMTransfo%hessian_read .AND. mole%NMTransfo%k_read) THEN
        ib = 0
        DO i=1,mole%nb_var
          IF (Ind_Coord_PerBlock(i) == Ind_Coord_AtBlock) THEN
            ib = ib + 1
            jb = 0
            DO j=1,mole%nb_var
              IF (Ind_Coord_PerBlock(j) == Ind_Coord_AtBlock) THEN
                jb = jb + 1
                d0k(ib,jb) = mole%NMTransfo%d0k(i,j)
                d0h(ib,jb) = mole%NMTransfo%d0h(i,j)
              END IF
            END DO
          END IF
        END DO

      ELSE ! both are false
        !- calc G_1
        CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,0)

        CALL get_dng_dnGG(Qact,para_Tnum,mole,dnGG=dnGG,nderiv=0)

        IF (nb_NM /= mole%nb_act) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nb_NM /= mole%nb_act',nb_NM,mole%nb_act
          STOP
        END IF

        d0k(:,:) = dnGG%d0(1:nb_NM,1:nb_NM)

        CALL dealloc_dnSVM(dnGG)


        !- calculation of the hessian (mole)
        IF (mole%NMTransfo%hessian_old) THEN
          IF (mole%NMTransfo%hessian_onthefly) THEN
            mole%NMTransfo%hessian_cart = .TRUE.
            ! save on-the-fly parameters
            name_FChk          = para_PES%para_OTF%file_FChk%name
            Read_OnTheFly_only = para_PES%Read_OnTheFly_only
            OnTheFly           = para_PES%OnTheFly

            ! set-up on-the-fly parameters to read the hessian
            para_PES%OnTheFly                = .TRUE.
            para_PES%Read_OnTheFly_only      = .TRUE.
            para_PES%para_OTF%file_FChk%name = mole%NMTransfo%file_hessian%name

            write(out_unitp,*) 'Read ab initio hessian from file: ',    &
                                  trim(para_PES%para_OTF%file_FChk%name)

            CALL Init_Tab_OF_dnMatOp(dnMatOp,nb_NM,1,nderiv=2)
            CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole,para_Tnum,para_PES)
            CALL Get_Hess_FROM_Tab_OF_dnMatOp(d0h,dnMatOp)
            CALL Get_Grad_FROM_Tab_OF_dnMatOp(d0grad,dnMatOp)
            CALL dealloc_Tab_OF_dnMatOp(dnMatOp)

            ! restore the on-the-fly parameters
            para_PES%para_OTF%file_FChk%name = name_FChk
            para_PES%Read_OnTheFly_only      = Read_OnTheFly_only
            para_PES%OnTheFly                = OnTheFly
          ELSE
            IF (mole%NMTransfo%hessian_cart) THEN
              CALL alloc_MatOFdnS(dnE,nb_NM,2)
              CALL alloc_MatOFdnS(dnECC,mole%ncart_act,2)

              write(out_unitp,*) 'Old hessian : mole%ncart_act',mole%ncart_act
              dnECC(1,1)%d1(:)   = ZERO
              IF (.NOT. l_hCC) THEN
                dnECC(1,1)%d2(:,:)   = ZERO
                CALL sub_hessian(dnECC(1,1)%d2)
              ELSE
                dnECC(1,1)%d2(:,:) = hCC(:,:)
              END IF
              CALL sub_dnFCC_TO_dnFcurvi(Qact,dnECC(1,1),dnE(1,1),mole)
              d0h(:,:) = dnE(1,1)%d2(:,:)

              CALL dealloc_MatOFdnS(dnE)
              CALL dealloc_MatOFdnS(dnECC)
            ELSE
              write(out_unitp,*) 'hessian : mole%nb_act',mole%nb_act
              CALL sub_hessian(d0h)
            END IF
            d0grad(:) = ZERO
          END IF
        ELSE

          nb_scalar_Op            = para_PES%nb_scalar_Op
          para_PES%nb_scalar_Op   = 0
          calc_scalar_Op          = para_PES%calc_scalar_Op
          para_PES%calc_scalar_Op = .FALSE.

          CALL Init_Tab_OF_dnMatOp(dnMatOp,nb_NM,1,nderiv=2)
          CALL get_dnMatOp_AT_Qact(Qact,dnMatOp,mole,para_Tnum,para_PES)

          CALL Get_Hess_FROM_Tab_OF_dnMatOp(d0h,dnMatOp)
          CALL Get_Grad_FROM_Tab_OF_dnMatOp(d0grad,dnMatOp)
          CALL dealloc_Tab_OF_dnMatOp(dnMatOp)

          para_PES%nb_scalar_Op   = nb_scalar_Op
          para_PES%calc_scalar_Op = calc_scalar_Op


        END IF

        IF (debug) THEN
          write(out_unitp,*) 'Qref (Qact)',Qact
          write(out_unitp,*) 'pot_Qref',dnE%d0
          CALL flush_perso(out_unitp)
        END IF

        IF (para_Tnum%WriteT .OR. debug .OR. print_level > 0) THEN
          write(out_unitp,*) 'gradient:'
          DO i=1,nb_NM
            write(out_unitp,'(a,x,i0,x,f12.6)') 'Q',i,d0grad(i)
          END DO
          CALL flush_perso(out_unitp)
        END IF

      END IF

      !write with high precision to be able to read them
      IF (para_Tnum%WriteT .OR. debug .OR. print_level > 1) THEN
        write(out_unitp,*) 'hessian matrix (not purified)'
        write(out_unitp,*) nb_NM,5
        !CALL Write_Mat(d0h,out_unitp,5,Rformat='e20.13')
        CALL Write_Mat(d0h,out_unitp,5,Rformat='f12.6')
        write(out_unitp,*) 'kinetic matrix  (not purified)'
        write(out_unitp,*) nb_NM,5
        !CALL Write_Mat(d0k,out_unitp,5,Rformat='e20.13')
        CALL Write_Mat(d0k,out_unitp,5,Rformat='f12.6')
        CALL flush_perso(out_unitp)
      END IF

      IF (debug) THEN
        d0hh = d0h
        CALL diagonalization(d0hh,d0eh,d0ch,nb_NM,2,1,.TRUE.)
        DO i=1,nb_NM,3
          write(out_unitp,*) i,'d0eh:',d0eh(i:min(i+2,nb_NM))
        END DO
        CALL flush_perso(out_unitp)
      END IF


      IF (print_level > 1) THEN
        write(out_unitp,*) '========================================='
        write(out_unitp,*) '==== END hessian and kinetic matrices ==='
        write(out_unitp,*) '========================================='
        write(out_unitp,*) '========================================='
        write(out_unitp,*)
        CALL flush_perso(out_unitp)
      END IF

!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF

!     -----------------------------------------------------------------

      END SUBROUTINE get_hess_k

      !=============================================================
      !     Harmonic potential part : pot2
      !     h(nb_inact2,nb_inact2)  : hessian matrix
      !     deltaQ(nb_inact2)       : displacement
      !=============================================================
      FUNCTION pot2(h,deltaQ,nb_inact2)
      USE mod_system
      IMPLICIT NONE

      real (kind=Rkind) :: pot2

      integer           :: nb_inact2
      real (kind=Rkind) :: deltaQ(nb_inact2)
      real (kind=Rkind) :: h(nb_inact2,nb_inact2)
      real (kind=Rkind) :: v2

      integer       :: i,j

      !----- for debuging ----------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      !---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING pot2'
        write(out_unitp,*) '  h ',h
        write(out_unitp,*) '  deltaQ',deltaQ
      END IF
      !---------------------------------------------------------------------

      v2 = ZERO
      DO i=1,nb_inact2
        v2 = v2 + h(i,i)*deltaQ(i)*deltaQ(i)
      END DO
      v2 = v2*HALF

      DO i=1,nb_inact2
        DO j=i+1,nb_inact2
          v2 = v2 + h(i,j)*deltaQ(i)*deltaQ(j)
        END DO
      END DO

      pot2 = v2

      !---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) '  v2 ',v2
        write(out_unitp,*) 'END pot2'
      END IF
      !---------------------------------------------------------------------

      END FUNCTION pot2

      SUBROUTINE Set_RPHpara_AT_Qact1(RPHpara_AT_Qact1,             &
                                      Qact,para_Tnum,mole,RPHTransfo)
      USE mod_system
      USE mod_Coord_KEO
      IMPLICIT NONE


      !----- for the zmatrix and Tnum --------------------------------------
      TYPE (Type_RPHpara_AT_Qact1), intent(inout) :: RPHpara_AT_Qact1

      TYPE (Tnum)             :: para_Tnum
      TYPE (zmatrix)          :: mole
      TYPE (Type_RPHTransfo)  :: RPHTransfo

      real (kind=Rkind), intent(inout) :: Qact(:)

      real (kind=Rkind)                :: det,over
      integer :: iQ,ndim
      real (kind=Rkind), allocatable   :: vecNM1(:),vecNM2(:)

      TYPE (Type_RPHpara_AT_Qact1) :: RPHpara_AT_0PlusStep,RPHpara_AT_0MinusStep


      !-----------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Set_RPHpara_AT_Qact1'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      !-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL flush_perso(out_unitp)
      END IF
      !-----------------------------------------------------------

      IF (RPHTransfo%option == 2) THEN

        CALL Set_RPHpara_AT_Qact1_opt2(RPHpara_AT_Qact1,                &
                                        Qact,para_Tnum,mole,RPHTransfo)

      ELSE ! option 0 ou 1

!        IF (RPHTransfo%nb_act1 == 1 .AND. abs(Qact(1)) < ONETENTH**5) THEN
!
!          ndim = RPHTransfo%nb_inact21
!          !CALL RPHpara1_AT_Qact1_TO_RPHpara2_AT_Qact1(RPHpara_AT_Qact1,RPHpara_AT_0PlusStep)
!          !CALL RPHpara1_AT_Qact1_TO_RPHpara2_AT_Qact1(RPHpara_AT_Qact1,RPHpara_AT_0MinusStep)
!
!          Qact(1) = ONETENTH**2
!          CALL Set_RPHpara_AT_Qact1_opt01(RPHpara_AT_0PlusStep,         &
!                                          Qact,para_Tnum,mole,RPHTransfo)
!          Qact(1) = -ONETENTH**2
!          CALL Set_RPHpara_AT_Qact1_opt01(RPHpara_AT_0MinusStep,         &
!                                          Qact,para_Tnum,mole,RPHTransfo)
!
!          ! copy for the allocation ...
!          CALL RPHpara1_AT_Qact1_TO_RPHpara2_AT_Qact1(RPHpara_AT_0PlusStep,RPHpara_AT_Qact1)
!          RPHpara_AT_Qact1%Qact1(:) = ZERO
!
!          ! average of Qopt
!          CALL dnVec2_wPLUS_dnVec3_TO_dnVec1(RPHpara_AT_Qact1%dnQopt,1,  &
!                                         RPHpara_AT_0PlusStep%dnQopt,1,HALF,&
!                                        RPHpara_AT_0MinusStep%dnQopt,1,HALF,&
!                                               ndim)
!          ! average of dnC_inv
!          CALL dnMat1_PLUS_dnMat2_TO_dnMat3(RPHpara_AT_0PlusStep%dnC_inv,&
!                                            RPHpara_AT_0MinusStep%dnC_inv,&
!                                            RPHpara_AT_Qact1%dnC_inv,     &
!                                            w1=HALF,w2=HALF)
!
!          CALL dealloc_RPHpara_AT_Qact1(RPHpara_AT_0PlusStep)
!          CALL dealloc_RPHpara_AT_Qact1(RPHpara_AT_0MinusStep)
!
!          !STOP 'Qact(1) = ZERO'
!        ELSE
          CALL Set_RPHpara_AT_Qact1_opt01(RPHpara_AT_Qact1,               &
                                          Qact,para_Tnum,mole,RPHTransfo)
!        END IF
     END IF


     CALL Det_OF_m1(RPHpara_AT_Qact1%dnC_inv%d0,det,RPHTransfo%nb_inact21,0)
     IF (debug) write(out_unitp,*) 'det of dnC_inv',det


    !check the sign with respect to iref-1 and changes it when negative
    CALL alloc_NParray(vecNM1,(/ RPHTransfo%nb_inact21 /),'vecNM1',name_sub)
    CALL alloc_NParray(vecNM2,(/ RPHTransfo%nb_inact21 /),'vecNM2',name_sub)

    DO iq=1,RPHTransfo%nb_inact21

      vecNM1 = RPHTransfo%RPHpara_AT_Qref(1)%dnC_inv%d0(iq,:)
      vecNM1 = vecNM1 / sqrt(dot_product(vecNM1,vecNM1))

      vecNM2 = RPHpara_AT_Qact1%dnC_inv%d0(iq,:)
      vecNM2 = vecNM2 / sqrt(dot_product(vecNM2,vecNM2))

      over = dot_product(vecNM1,vecNM2)

      IF (debug) write(out_unitp,*) 'over',iq,over

    END DO
    CALL dealloc_NParray(vecNM1,'vecNM1',name_sub)
    CALL dealloc_NParray(vecNM2,'vecNM2',name_sub)

     IF (debug) THEN
        CALL Write_RPHpara_AT_Qact1(RPHpara_AT_Qact1)
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
     END IF

     END SUBROUTINE Set_RPHpara_AT_Qact1

      SUBROUTINE Set_RPHpara_AT_Qact1_opt2(RPHpara_AT_Qact1,            &
                                           Qact,para_Tnum,mole,RPHTransfo)
      USE mod_system
      USE mod_dnSVM
      USE mod_Constant, only : get_Conv_au_TO_unit
      USE mod_Coord_KEO
      IMPLICIT NONE


      !----- for the zmatrix and Tnum --------------------------------------
      TYPE (Type_RPHpara_AT_Qact1), intent(inout) :: RPHpara_AT_Qact1
      integer :: nb_act1,nb_inact21

      TYPE (Tnum)             :: para_Tnum
      TYPE (zmatrix)          :: mole
      TYPE (Type_RPHTransfo)  :: RPHTransfo

      real (kind=Rkind), intent(inout) :: Qact(:)



      !------ for the frequencies -------------------------------
      integer               :: nderiv
      real (kind=Rkind)     :: auTOcm_inv
      integer               :: i,iact,idyn,RPHoption,iref,nb_ref
      real (kind=Rkind)     :: Qdyn(mole%nb_var)

      integer               :: iact1,iq,jq,iQinact21,jQinact21
      integer               :: listNM_selected(mole%nb_var)


      TYPE (Type_dnS), pointer   :: dnSwitch(:)
      TYPE (Type_dnS)            :: dnW1
      real (kind=Rkind)          :: sc,det

      TYPE (Type_dnVec)                  :: dnQact
      real (kind=Rkind), allocatable     :: QrefQact(:,:)     ! QrefQact(nb_Qact1,nb_ref)

      !-----------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Set_RPHpara_AT_Qact1_opt2'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      !-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL flush_perso(out_unitp)
      END IF
      !-----------------------------------------------------------
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')

      nderiv     = 3
      IF (para_Tnum%nrho == 0 .OR. para_Tnum%nrho == 10 .OR. para_Tnum%nrho == 20) nderiv = 2

      nb_act1    = RPHTransfo%nb_act1
      nb_inact21 = RPHTransfo%nb_inact21
      nb_ref     = RPHTransfo%RPHpara2%nb_ref

      CALL alloc_RPHpara_AT_Qact1(RPHpara_AT_Qact1,nb_act1,nb_inact21,nderiv)


      !here it should be Qin of RPH (therefore Qdyn ?????)
      CALL Qact_TO_Qdyn_FROM_ActiveTransfo(Qact,Qdyn,mole%ActiveTransfo)

      RPHpara_AT_Qact1%Qact1(:) = Qdyn(RPHTransfo%list_QactTOQdyn(1:nb_act1))

      ! 1st: dnQact (derivatives, just for the active coordinates)
      CALL alloc_dnSVM(dnQact,  nb_act1,nb_act1,           nderiv)
      dnQact%d0(:) = RPHpara_AT_Qact1%Qact1(:)
      CALL Set_AllActive(dnQact)

      ! 2d: the reference Qact
      CALL alloc_NParray(QrefQact,(/ nb_act1,nb_ref /),'QrefQact',name_sub)
      QrefQact(:,:) = RPHTransfo%RPHpara2%QoutRef(1:nb_act1,:)

      ! 3d: dnSwitch
      sc = TWO ! to be changed, from Read_RPHpara2
      nullify(dnSwitch)
      CALL alloc_array(dnSwitch,(/nb_ref/),"dnSwitch",name_sub)
      CALL alloc_VecOFdnS(dnSwitch,nb_act1,nderiv)
      CALL Switch_RPH(dnSwitch,dnQact,QrefQact,sc,nderiv)
      !write(6,*) 'dnSwitch(:)',dnSwitch(:)%d0

      CALL alloc_dnSVM(dnW1,  nb_act1,           nderiv)

      ! 4th: dnQopt
      !old (one ref)
      CALL sub_ZERO_TO_dnVec(RPHpara_AT_Qact1%dnQopt)
      DO iQinact21=1,nb_inact21
        CALL sub_ZERO_TO_dnS(dnW1)
        DO iref=1,nb_ref
          !dnW1 = dnW1 + dnSwitch(iref)*RPHTransfo%RPHpara2%QoutRef(nb_act1+iQinact21,iref)
          CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnW1,ONE,                    &
             dnSwitch(iref),RPHTransfo%RPHpara2%QoutRef(nb_act1+iQinact21,iref),&
                                           dnW1)
        END DO
        CALL sub_dnS_TO_dnVec(dnW1,RPHpara_AT_Qact1%dnQopt,iQinact21)
      END DO
      !write(99,*) 'Qact,Qopt',dnQact%d0(:),RPHpara_AT_Qact1%dnQopt%d0

      !5th: dnC_inv
      CALL sub_ZERO_TO_dnMat(RPHpara_AT_Qact1%dnC_inv)
      listNM_selected(:) = 0
      DO iact1=1,nb_act1
        listNM_selected(RPHTransfo%RPHpara2%listNM_act1(iact1)) = 1
      END DO

      iQinact21 = 0
      DO iq=1,nb_act1+nb_inact21
        IF (listNM_selected(iq) /= 0) CYCLE
        iQinact21 = iQinact21 + 1

        DO jQinact21=1,nb_inact21

          CALL sub_ZERO_TO_dnS(dnW1)
          DO iref=1,nb_ref
            CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnW1,ONE,                  &
                dnSwitch(iref),RPHTransfo%RPHpara2%CinvRef(iq,nb_act1+jQinact21,iref), &
                                             dnW1)
          END DO
          CALL sub_dnS_TO_dnMat(dnW1,RPHpara_AT_Qact1%dnC_inv,iQinact21,jQinact21)


        END DO
      END DO

      CALL dealloc_dnSVM(dnQact)
      CALL dealloc_NParray(QrefQact,'QrefQact',name_sub)

      CALL dealloc_VecOFdnS(dnSwitch)
      CALL dealloc_array(dnSwitch,"dnSwitch",name_sub)
      nullify(dnSwitch)

      CALL dealloc_dnSVM(dnW1)


      ! just dnC%d0
      CALL inv_m1_TO_m2(RPHpara_AT_Qact1%dnC_inv%d0,RPHpara_AT_Qact1%dnC%d0, &
                        nb_inact21,0,ZERO)

     IF (debug) THEN
        CALL Write_RPHpara_AT_Qact1(RPHpara_AT_Qact1)
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
     END IF

     END SUBROUTINE Set_RPHpara_AT_Qact1_opt2

      SUBROUTINE Set_RPHpara_AT_Qact1_opt01(RPHpara_AT_Qact1,           &
                                      Qact,para_Tnum,mole,RPHTransfo)
      USE mod_system
      USE mod_dnSVM
      USE mod_Constant, only : get_Conv_au_TO_unit
      USE mod_Coord_KEO
      IMPLICIT NONE


      !----- for the zmatrix and Tnum --------------------------------------
      TYPE (Type_RPHpara_AT_Qact1), intent(inout) :: RPHpara_AT_Qact1
      integer :: nb_act1,nb_inact21

      TYPE (Tnum)             :: para_Tnum
      TYPE (zmatrix)          :: mole
      TYPE (Type_RPHTransfo)  :: RPHTransfo

      real (kind=Rkind), intent(inout) :: Qact(:)



      !------ for the frequencies -------------------------------
      TYPE (Type_dnMat)     :: dnC,dnC_inv      ! derivative with respect to Qact1
      TYPE (Type_dnVec)     :: dnQeq            ! derivative with respect to Qact1
      TYPE (Type_dnVec)     :: dnEHess          ! derivative with respect to Qact1
      TYPE (Type_dnVec)     :: dnGrad           ! derivative with respect to Qact1
      TYPE (Type_dnMat)     :: dnHess           ! derivative with respect to Qact1
      TYPE (Type_dnS)       :: dnLnN            ! derivative with respect to Qact1
      integer               :: nderiv
      real (kind=Rkind)     :: pot0_corgrad,stepp,step_loc,vi,auTOcm_inv
      integer               :: i,iact,idyn,RPHoption
      real (kind=Rkind)     :: Qdyn(mole%nb_var)

      !-----------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Set_RPHpara_AT_Qact1_opt01'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      !-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL flush_perso(out_unitp)
      END IF
      !-----------------------------------------------------------
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')

      nderiv     = 3
      IF (para_Tnum%nrho == 0 .OR. para_Tnum%nrho == 10 .OR. para_Tnum%nrho == 20) nderiv = 2

      step_loc = RPHTransfo%step
      stepp    = ONE/(step_loc+step_loc)

      nb_act1    = RPHTransfo%nb_act1
      nb_inact21 = RPHTransfo%nb_inact21

      CALL alloc_RPHpara_AT_Qact1(RPHpara_AT_Qact1,nb_act1,nb_inact21,nderiv)



      RPHoption = RPHTransfo%option
      CALL Sub_paraRPH_TO_mole(mole) ! switch back mole
      mole%tab_Qtransfo(mole%itRPH)%skip_transfo = .TRUE. ! we have to skip RPH transfo because, ...
                                                          ! this subroutine calculates the RPH parameters

      CALL alloc_dnSVM(dnC,    nb_inact21,nb_inact21,nb_act1,nderiv)
      CALL alloc_dnSVM(dnC_inv,nb_inact21,nb_inact21,nb_act1,nderiv)
      CALL alloc_dnSVM(dnQeq,  nb_inact21,nb_act1,           nderiv)
      CALL alloc_dnSVM(dnEHess,nb_inact21,nb_act1,           nderiv)
      CALL alloc_dnSVM(dnHess, nb_inact21,nb_inact21,nb_act1,nderiv)
      CALL alloc_dnSVM(dnGrad, nb_inact21,nb_act1,           nderiv)
      CALL alloc_dnSVM(dnLnN,  nb_act1,                      nderiv)

      !here it should be Qin of RPH (therefore Qdyn ?????)
      CALL Qact_TO_Qdyn_FROM_ActiveTransfo(Qact,Qdyn,mole%ActiveTransfo)

      RPHpara_AT_Qact1%Qact1(:) = Qdyn(RPHTransfo%list_QactTOQdyn(1:nb_act1))

      CALL sub_dnfreq_8p(RPHpara_AT_Qact1%dnQopt,RPHpara_AT_Qact1%dnC,  &
                        RPHpara_AT_Qact1%dnLnN,RPHpara_AT_Qact1%dnEHess,&
                        RPHpara_AT_Qact1%dnhess,dnGrad,                 &
                        RPHpara_AT_Qact1%dnC_inv,pot0_corgrad,          &
                        Qact,para_Tnum,mole,RPHTransfo,nderiv,.FALSE.)

     IF (debug) THEN
       write(out_unitp,*) 'dnC_inv'
       CALL Write_dnMat(RPHpara_AT_Qact1%dnC_inv,nderiv=0)
     END IF

     IF (nderiv == 3) THEN
       DO i=1,nb_act1

         idyn = RPHTransfo%list_QactTOQdyn(i)
         iact = mole%liste_QsymTOQact(idyn)
         vi = Qact(iact)

         !-- frequencies calculation at Qact(i)+step -------------
         Qact(iact)      = vi + step_loc

         CALL sub_dnfreq_8p(dnQeq,dnC,dnLnN,dnEHess,dnhess,dnGrad,dnC_inv, &
                         pot0_corgrad,Qact,    &
                         para_Tnum,mole,RPHTransfo,nderiv,.FALSE.)

         RPHpara_AT_Qact1%dnLnN%d3(:,:,i)         = dnLnN%d2(:,:)

         RPHpara_AT_Qact1%dnC%d3(:,:,:,:,i)       = dnC%d2(:,:,:,:)
         RPHpara_AT_Qact1%dnC_inv%d3(:,:,:,:,i)   = dnC_inv%d2(:,:,:,:)
         RPHpara_AT_Qact1%dnhess%d3(:,:,:,:,i)    = dnHess%d2(:,:,:,:)

         RPHpara_AT_Qact1%dnEHess%d3(:,:,:,i)     = dnEHess%d2(:,:,:)
         RPHpara_AT_Qact1%dnQopt%d3(:,:,:,i)      = dnQeq%d2(:,:,:)

         !-- frequencies calculation at Qact(i)-step -------------
         Qact(iact)      = vi - step_loc

         CALL sub_dnfreq_8p(dnQeq,dnC,dnLnN,dnEHess,dnhess,dnGrad,dnC_inv, &
                            pot0_corgrad,Qact,    &
                            para_Tnum,mole,RPHTransfo,nderiv,.FALSE.)

         RPHpara_AT_Qact1%dnLnN%d3(:,:,i)         =                     &
                                (RPHpara_AT_Qact1%dnLnN%d3(:,:,i)      -&
                                                    dnLnN%d2(:,:))*stepp

         RPHpara_AT_Qact1%dnC%d3(:,:,:,:,i)       =                     &
                                (RPHpara_AT_Qact1%dnC%d3(:,:,:,:,i)    -&
                                                  dnC%d2(:,:,:,:))*stepp
         RPHpara_AT_Qact1%dnC_inv%d3(:,:,:,:,i)   =                     &
                                (RPHpara_AT_Qact1%dnC_inv%d3(:,:,:,:,i)-&
                                              dnC_inv%d2(:,:,:,:))*stepp
         RPHpara_AT_Qact1%dnhess%d3(:,:,:,:,i)    =                     &
                                (RPHpara_AT_Qact1%dnhess%d3(:,:,:,:,i) -&
                                               dnHess%d2(:,:,:,:))*stepp

         RPHpara_AT_Qact1%dnEHess%d3(:,:,:,i)     =                     &
                                (RPHpara_AT_Qact1%dnEHess%d3(:,:,:,i)  -&
                                                dnEHess%d2(:,:,:))*stepp
         RPHpara_AT_Qact1%dnQopt%d3(:,:,:,i)      =                     &
                               (RPHpara_AT_Qact1%dnQopt%d3(:,:,:,i)    -&
                                                  dnQeq%d2(:,:,:))*stepp

         Qact(iact)      = vi
       END DO
     END IF

     write(out_unitp,11) Qact(1:RPHTransfo%nb_act1),                    &
                               RPHpara_AT_Qact1%dnEHess%d0(:)*auTOcm_inv
 11  format(' frequencies : ',30f10.4)

     CALL dealloc_dnSVM(dnC)
     CALL dealloc_dnSVM(dnC_inv)
     CALL dealloc_dnSVM(dnQeq)
     CALL dealloc_dnSVM(dnEHess)
     CALL dealloc_dnSVM(dnHess)
     CALL dealloc_dnSVM(dnGrad)
     CALL dealloc_dnSVM(dnLnN)

     mole%tab_Qtransfo(mole%itRPH)%skip_transfo = .FALSE.
     CALL Sub_mole_TO_paraRPH(mole)
     RPHTransfo%option = RPHoption

     IF (debug) THEN
        CALL Write_RPHpara_AT_Qact1(RPHpara_AT_Qact1)
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
     END IF

     END SUBROUTINE Set_RPHpara_AT_Qact1_opt01
!
!=============================================================
!
!     derivative frequency calculations
!
!=============================================================

      SUBROUTINE sub_dnfreq_8p(dnQeq,dnC,dnLnN,dnEHess,dnHess,dnGrad,dnC_inv,&
                               pot0_corgrad,Qact,                      &
                               para_Tnum,mole,RPHTransfo,nderiv,test)
      USE mod_system
      USE mod_dnSVM
      USE mod_Constant, only : get_Conv_au_TO_unit
      USE mod_Coord_KEO
      IMPLICIT NONE

!----- for the zmatrix and Tnum --------------------------------------
      TYPE (Tnum)    :: para_Tnum
      TYPE (zmatrix) :: mole

      real (kind=Rkind), intent(inout) :: Qact(:)

!----- variables for the active and inactive namelists ----------------
      TYPE (Type_RPHTransfo)  :: RPHTransfo
!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
      integer :: nderiv

!------ for the frequencies -------------------------------
        TYPE (Type_dnMat)     :: dnC,dnC_inv      ! derivative with respect to Qact1
        TYPE (Type_dnVec)     :: dnQeq            ! derivative with respect to Qact1
        TYPE (Type_dnVec)     :: dnEHess          ! derivative with respect to Qact1
        TYPE (Type_dnVec)     :: dnGrad           ! derivative with respect to Qact1
        TYPE (Type_dnMat)     :: dnHess           ! derivative with respect to Qact1
        TYPE (Type_dnS)       :: dnLnN            ! derivative with respect to Qact1

      real (kind=Rkind) :: pot0_corgrad,pot0_corgrad2


!----- pour les derivees ---------------------------------------------
      real (kind=Rkind) ::    step,step2,stepp,step24
      real (kind=Rkind) ::    d1


!----- for testing ---------------------------------------------------
      logical :: test


!----- working variables ---------------------------------------------
      integer           :: i,j,k,nb_inact21,nb_act1
      real (kind=Rkind) :: vi,vj

      real (kind=Rkind) ::  mat0(RPHTransfo%nb_inact21,RPHTransfo%nb_inact21)
      real (kind=Rkind) ::  mat1(RPHTransfo%nb_inact21,RPHTransfo%nb_inact21)
      real (kind=Rkind) ::  mat2(RPHTransfo%nb_inact21,RPHTransfo%nb_inact21)
      real (kind=Rkind) ::  vec1(RPHTransfo%nb_inact21),vec0(RPHTransfo%nb_inact21)
      real (kind=Rkind) ::  vec2(RPHTransfo%nb_inact21)

      real (kind=Rkind) ::  mat2p(RPHTransfo%nb_inact21,RPHTransfo%nb_inact21)
      real (kind=Rkind) ::  mat2m(RPHTransfo%nb_inact21,RPHTransfo%nb_inact21)
      real (kind=Rkind) ::  mat22p(RPHTransfo%nb_inact21,RPHTransfo%nb_inact21)
      real (kind=Rkind) ::  mat22m(RPHTransfo%nb_inact21,RPHTransfo%nb_inact21)

      real (kind=Rkind) ::  mat2_s(RPHTransfo%nb_inact21,RPHTransfo%nb_inact21)
      real (kind=Rkind) ::  mat2_s2(RPHTransfo%nb_inact21,RPHTransfo%nb_inact21)



      real (kind=Rkind) ::  vec0p(RPHTransfo%nb_inact21)
      real (kind=Rkind) ::  vec0m(RPHTransfo%nb_inact21)
      real (kind=Rkind) ::  vec02p(RPHTransfo%nb_inact21)
      real (kind=Rkind) ::  vec02m(RPHTransfo%nb_inact21)

      real (kind=Rkind) ::  vec0_s(RPHTransfo%nb_inact21)
      real (kind=Rkind) ::  vec0_s2(RPHTransfo%nb_inact21)

      real (kind=Rkind) ::  auTOcm_inv


!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_dnfreq_8p'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Qact',Qact
        write(out_unitp,*) 'purify_hess,eq_hess',                       &
                              RPHTransfo%purify_hess,RPHTransfo%eq_hess
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')

      step       = RPHTransfo%step
      step2      = step * HALF

      nb_inact21 = RPHTransfo%nb_inact21
      nb_act1    = RPHTransfo%nb_act1

      IF (RPHTransfo%step <= ZERO) THEN
        write(out_unitp,*) ' ERROR : RPHTransfo%step is < zero'
        STOP
      END IF


      CALL check_alloc_dnMat(dnC,'dnC',name_sub)
      CALL check_alloc_dnMat(dnC_inv,'dnC_inv',name_sub)
      CALL check_alloc_dnVec(dnQeq,'dnQeq',name_sub)
      CALL check_alloc_dnVec(dnEHess,'dnEHess',name_sub)
      CALL check_alloc_dnMat(dnHess,'dnHess',name_sub)
      CALL check_alloc_dnVec(dnGrad,'dnGrad',name_sub)
      CALL check_alloc_dnS(dnLnN,'dnLnN',name_sub)

      IF (nderiv == 0) THEN
        CALL sub_freq2_RPH(dnEHess%d0,dnC%d0,dnC_inv%d0,        &
                       dnLnN%d0,dnHess%d0,dnQeq%d0,                     &
                       dnGrad%d0,pot0_corgrad,                          &
                       Qact,para_Tnum,mole,RPHTransfo)

        IF (debug) THEN
          write(out_unitp,*) 'dnQeq%d0',dnQeq%d0(:)
          write(out_unitp,*) 'freq',dnEHess%d0(:)*auTOcm_inv
          write(out_unitp,*) 'dnC'
          CALL Write_dnMat(dnC)
          write(out_unitp,*) 'dnHess'
          CALL Write_dnMat(dnHess)
          write(out_unitp,*) 'END ',name_sub
        END IF
        RETURN
      END IF


!-----------------------------------------------------------------
!----- frequencies calculation at Qact --------------------------
!-----------------------------------------------------------------

      CALL sub_freq2_RPH(dnEHess%d0,dnC%d0,dnC_inv%d0,                  &
                         dnLnN%d0,dnHess%d0,dnQeq%d0,                   &
                         dnGrad%d0,pot0_corgrad,                        &
                         Qact,para_Tnum,mole,RPHTransfo)

!-----------------------------------------------------------------
!----- end frequencies calculation at Qact ----------------------
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!----- d/Qqi et d2/dQi2 of frequencies ---------------------------
!-----------------------------------------------------------------
      DO i=1,RPHTransfo%nb_act1

        vi = Qact(i)

!       -- frequencies calculation at Qact(i)+step -------------
        Qact(i) = vi + step

        CALL sub_freq2_RPH(vec1,mat1,mat2,                      &
                       d1,mat0,vec0,                                    &
                       vec2,pot0_corgrad2,                              &
                       Qact,para_Tnum,mole,RPHTransfo)

        mat2p   = mat2-dnC_inv%d0
        vec0p   = vec0-dnQeq%d0


!       -- frequencies calculation at Qact(i)-step -------------
        Qact(i) = vi - step

        CALL sub_freq2_RPH(vec1,mat1,mat2,                      &
                       d1,mat0,vec0,                                    &
                       vec2,pot0_corgrad2,                              &
                       Qact,para_Tnum,mole,RPHTransfo)

        mat2m   = mat2-dnC_inv%d0
        vec0m   = vec0-dnQeq%d0


!       -- frequencies calculation at Qact(i)+step -------------
        Qact(i) = vi + step2

        CALL sub_freq2_RPH(vec1,mat1,mat2,                      &
                       d1,mat0,vec0,                                    &
                       vec2,pot0_corgrad2,                              &
                       Qact,para_Tnum,mole,RPHTransfo)

        mat22p   = mat2-dnC_inv%d0
        vec02p   = vec0-dnQeq%d0


!       -- frequencies calculation at Qact(i)-step -------------
        Qact(i) = vi - step2

        CALL sub_freq2_RPH(vec1,mat1,mat2,                      &
                       d1,mat0,vec0,                                    &
                       vec2,pot0_corgrad2,                              &
                       Qact,para_Tnum,mole,RPHTransfo)

        mat22m   = mat2-dnC_inv%d0
        vec02m   = vec0-dnQeq%d0



        dnC_inv%d1(:,:,i)   = (EIGHT*(mat22p-mat22m)-(mat2p-mat2m))/(SIX*step)
        dnC_inv%d2(:,:,i,i) = (16._Rkind*(mat22p+mat22m)-(mat2p+mat2m)) / (THREE*step*step)

        dnQeq%d1(:,i)     = (EIGHT*(vec02p-vec02m)-(vec0p-vec0m))/(SIX*step)
        dnQeq%d2(:,i,i)   = (16._Rkind*(vec02p+vec02m)-(vec0p+vec0m)) / (THREE*step*step)

        Qact(i) = vi
      END DO


!-----------------------------------------------------------------
!----- end d/Qqi and d2/dQi2 of frequencies ----------------------
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!----- d2/dQidQj of frequencies (4 points) -----------------------
!      d2/dQidQj = ( v(Qi+,Qj+)+v(Qi-,Qj-)-v(Qi-,Qj+)-v(Qi+,Qj-) )/(4*s*s)
!-----------------------------------------------------------------
      DO i=1,RPHTransfo%nb_act1
      DO j=i+1,RPHTransfo%nb_act1

        vi = Qact(i)
        vj = Qact(j)


!       -- frequencies calculation at Qact(i)+step Qact(j)+step
        Qact(i) = vi + step
        Qact(j) = vj + step
        CALL sub_freq2_RPH(vec1,mat1,mat2,                      &
                       d1,mat0,vec0,                                    &
                       vec2,pot0_corgrad2,                              &
                       Qact,para_Tnum,mole,RPHTransfo)

        mat2_s   = mat2
        vec0_s   = vec0



!       -- frequencies calculation at Qact(i)-step Qact(j)-step
        Qact(i) = vi - step
        Qact(j) = vj - step

        CALL sub_freq2_RPH(vec1,mat1,mat2,                      &
                       d1,mat0,vec0,                                    &
                       vec2,pot0_corgrad2,                              &
                       Qact,para_Tnum,mole,RPHTransfo)

        mat2_s   = mat2_s + mat2
        vec0_s   = vec0_s + vec0


!       -- frequencies calculation at Qact(i)-step Qact(j)+step
        Qact(i) = vi - step
        Qact(j) = vj + step

        CALL sub_freq2_RPH(vec1,mat1,mat2,                      &
                       d1,mat0,vec0,                                    &
                       vec2,pot0_corgrad2,                              &
                       Qact,para_Tnum,mole,RPHTransfo)

        mat2_s   = mat2_s - mat2
        vec0_s   = vec0_s - vec0

!       -- frequencies calculation at Qact(i)+step Qact(j)-step
        Qact(i) = vi + step
        Qact(j) = vj - step

        CALL sub_freq2_RPH(vec1,mat1,mat2,                      &
                       d1,mat0,vec0,                                    &
                       vec2,pot0_corgrad2,                              &
                       Qact,para_Tnum,mole,RPHTransfo)

        mat2_s   = mat2_s - mat2
        vec0_s   = vec0_s - vec0


!       -- frequencies calculation at Qact(i)+step Qact(j)+step
        Qact(i) = vi + step2
        Qact(j) = vj + step2
        CALL sub_freq2_RPH(vec1,mat1,mat2,                      &
                       d1,mat0,vec0,                                    &
                       vec2,pot0_corgrad2,                              &
                       Qact,para_Tnum,mole,RPHTransfo)

        mat2_s2   = mat2
        vec0_s2   = vec0



!       -- frequencies calculation at Qact(i)-step Qact(j)-step
        Qact(i) = vi - step2
        Qact(j) = vj - step2

        CALL sub_freq2_RPH(vec1,mat1,mat2,                      &
                       d1,mat0,vec0,                                    &
                       vec2,pot0_corgrad2,                              &
                       Qact,para_Tnum,mole,RPHTransfo)

        mat2_s2   = mat2_s2 + mat2
        vec0_s2   = vec0_s2 + vec0


!       -- frequencies calculation at Qact(i)-step Qact(j)+step
        Qact(i) = vi - step2
        Qact(j) = vj + step2

        CALL sub_freq2_RPH(vec1,mat1,mat2,                      &
                       d1,mat0,vec0,                                    &
                       vec2,pot0_corgrad2,                              &
                       Qact,para_Tnum,mole,RPHTransfo)

        mat2_s2   = mat2_s2 - mat2
        vec0_s2   = vec0_s2 - vec0

!       -- frequencies calculation at Qact(i)+step Qact(j)-step
        Qact(i) = vi + step2
        Qact(j) = vj - step2

        CALL sub_freq2_RPH(vec1,mat1,mat2, d1,mat0,vec0,vec2,           &
                           pot0_corgrad2,                               &
                           Qact,para_Tnum,mole,RPHTransfo)

        mat2_s2   = mat2_s2 - mat2
        vec0_s2   = vec0_s2 - vec0


!       -- d2/dQi/dQj -----------------------------------------

        dnC_inv%d2(:,:,i,j) = (16._Rkind*mat2_s2 - mat2_s)/(step*step*TWELVE)
        dnQeq%d2(:,i,j)     = (16._Rkind*vec0_s2 - vec0_s)/(step*step*TWELVE)

        dnC_inv%d2(:,:,j,i) = dnC_inv%d2(:,:,i,j)
        dnQeq%d2(:,j,i)     = dnQeq%d2(:,i,j)


        Qact(i) = vi
        Qact(j) = vj
      END DO
      END DO
!-----------------------------------------------------------------
!----- end d2/dQidQj of frequencies ------------------------------
!-----------------------------------------------------------------

!-----------------------------------------------------------

       IF (debug .OR. test) THEN
         write(out_unitp,11)                         &
                  Qact(1:RPHTransfo%nb_act1),dnEHess%d0(:)*auTOcm_inv
 11      format(' frequencies : ',30f10.4)
         write(out_unitp,*) 'dnQeq'
         CALL Write_dnVec(dnQeq)
         write(out_unitp,*) 'dnC_inv'
         CALL Write_dnMat(dnC_inv)
       END IF

       IF (debug) THEN
         write(out_unitp,*) 'END ',name_sub
       END IF
      CALL flush_perso(out_unitp)
!-----------------------------------------------------------

      END SUBROUTINE sub_dnfreq_8p
      SUBROUTINE sub_dnfreq_4p(dnQeq,dnC,dnLnN,dnEHess,dnHess,dnGrad,dnC_inv,&
                               pot0_corgrad,Qact,                       &
                               para_Tnum,mole,RPHTransfo,nderiv,test)
      USE mod_system
      USE mod_dnSVM
      USE mod_Constant, only : get_Conv_au_TO_unit
      USE mod_Coord_KEO
      IMPLICIT NONE

!----- for the zmatrix and Tnum --------------------------------------
      TYPE (Tnum)             :: para_Tnum
      TYPE (zmatrix)          :: mole
      TYPE (Type_RPHTransfo)  :: RPHTransfo

      real (kind=Rkind), intent(inout) :: Qact(:)

!----- variables for the active and inactive namelists ----------------

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
      integer :: nderiv

!------ for the frequencies -------------------------------
      TYPE (Type_dnMat)     :: dnC,dnC_inv      ! derivative with respect to Qact1
      TYPE (Type_dnVec)     :: dnQeq            ! derivative with respect to Qact1
      TYPE (Type_dnVec)     :: dnEHess          ! derivative with respect to Qact1
      TYPE (Type_dnVec)     :: dnGrad           ! derivative with respect to Qact1
      TYPE (Type_dnMat)     :: dnHess           ! derivative with respect to Qact1
      TYPE (Type_dnS)       :: dnLnN            ! derivative with respect to Qact1

      real (kind=Rkind) :: pot0_corgrad,pot0_corgrad2

      real (kind=Rkind)  :: Qact1(RPHTransfo%nb_act1)


!----- pour les derivees ---------------------------------------------
      real (kind=Rkind) ::    step,step2,stepp,step24
      real (kind=Rkind) ::    d1


!----- for testing ---------------------------------------------------
      logical :: test


!----- working variables ---------------------------------------------
      integer           :: i,j,k,nb_inact21,nb_act1
      real (kind=Rkind) :: vi,vj

      real (kind=Rkind) ::  mat0(RPHTransfo%nb_inact21,RPHTransfo%nb_inact21)
      real (kind=Rkind) ::  mat1(RPHTransfo%nb_inact21,RPHTransfo%nb_inact21)
      real (kind=Rkind) ::  mat2(RPHTransfo%nb_inact21,RPHTransfo%nb_inact21)
      real (kind=Rkind) ::  vec1(RPHTransfo%nb_inact21),vec0(RPHTransfo%nb_inact21)
      real (kind=Rkind) ::  vec2(RPHTransfo%nb_inact21)

      real (kind=Rkind) ::  auTOcm_inv

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_dnfreq_4p'
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Qact',Qact
        write(out_unitp,*) 'purify_hess,eq_hess',                       &
                              RPHTransfo%purify_hess,RPHTransfo%eq_hess
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')

      step     = RPHTransfo%step
      step2    = ONE/(step*step)
      step24   = step2*HALF*HALF
      stepp    = ONE/(step+step)

      nb_inact21 = RPHTransfo%nb_inact21
      nb_act1    = RPHTransfo%nb_act1

      IF (RPHTransfo%step <= ZERO) THEN
        write(out_unitp,*) ' ERROR : RPHTransfo%step is <= to zero'
        STOP
      END IF


      CALL check_alloc_dnMat(dnC,'dnC',name_sub)
      CALL check_alloc_dnMat(dnC_inv,'dnC_inv',name_sub)
      CALL check_alloc_dnVec(dnQeq,'dnQeq',name_sub)
      CALL check_alloc_dnVec(dnEHess,'dnEHess',name_sub)
      CALL check_alloc_dnMat(dnHess,'dnHess',name_sub)
      CALL check_alloc_dnVec(dnGrad,'dnGrad',name_sub)
      CALL check_alloc_dnS(dnLnN,'dnLnN',name_sub)


      IF (nderiv == 0) THEN
        CALL sub_freq2_RPH(dnEHess%d0,dnC%d0,dnC_inv%d0,                &
                           dnLnN%d0,dnHess%d0,dnQeq%d0,                 &
                           dnGrad%d0,pot0_corgrad,                      &
                           Qact,para_Tnum,mole,RPHTransfo)

        IF (debug) THEN
          write(out_unitp,*) 'dnQeq%d0',dnQeq%d0(:)
          write(out_unitp,*) 'freq',dnEHess%d0(:)*auTOcm_inv
          write(out_unitp,*) 'dnC'
          CALL Write_dnMat(dnC)
          write(out_unitp,*) 'dnHess'
          CALL Write_dnMat(dnHess)
          write(out_unitp,*) 'END ',name_sub
        END IF
        RETURN
      END IF


!-----------------------------------------------------------------
!----- frequencies calculation at Qact --------------------------
!-----------------------------------------------------------------

      CALL sub_freq2_RPH(dnEHess%d0,dnC%d0,dnC_inv%d0,                  &
                         dnLnN%d0,dnHess%d0,dnQeq%d0,dnGrad%d0,         &
                         pot0_corgrad,Qact,para_Tnum,mole,RPHTransfo)

!-----------------------------------------------------------------
!----- end frequencies calculation at Qact ----------------------
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!----- d/Qqi et d2/dQi2 of frequencies ---------------------------
!-----------------------------------------------------------------
      DO i=1,RPHTransfo%nb_act1

        vi = Qact(i)

!       -- frequencies calculation at Qact(i)+step -------------
        Qact(i) = vi + step

        CALL sub_freq2_RPH(vec1,mat1,mat2,d1,mat0,vec0,vec2,    &
                       pot0_corgrad2,Qact,para_Tnum,mole,RPHTransfo)

        dnC_inv%d1(:,:,i)   = mat2
        dnC_inv%d2(:,:,i,i) = mat2

        dnLnN%d1(i)         = d1
        dnLnN%d2(i,i)       = d1

        dnC%d1(:,:,i)       = mat1(:,:)
        dnC%d2(:,:,i,i)     = mat1(:,:)

        dnHess%d1(:,:,i)    = mat0(:,:)
        dnHess%d2(:,:,i,i)  = mat0(:,:)

        dnQeq%d1(:,i)       = vec0(:)
        dnQeq%d2(:,i,i)     = vec0(:)

!       -- frequencies calculation at Qact(i)-step -------------
        Qact(i) = vi - step

        CALL sub_freq2_RPH(vec1,mat1,mat2,d1,mat0,vec0,vec2,    &
                       pot0_corgrad2,Qact,para_Tnum,mole,RPHTransfo)

        dnC_inv%d1(:,:,i)   = (dnC_inv%d1(:,:,i)   - mat2)*stepp
        dnC_inv%d2(:,:,i,i) = (dnC_inv%d2(:,:,i,i) + mat2 -TWO*dnC_inv%d0)*step2

        dnLnN%d1(i)        = (dnLnN%d1(i)   - d1)*stepp
        dnLnN%d2(i,i)      = (dnLnN%d2(i,i) + d1 -dnLnN%d0-dnLnN%d0)*step2

        dnC%d1(:,:,i)      = ( dnC%d1(:,:,i) -  mat1(:,:) ) * stepp
        dnC%d2(:,:,i,i)    = (dnC%d2(:,:,i,i)+mat1(:,:)-dnC%d0(:,:)-dnC%d0(:,:))*step2

        dnHess%d1(:,:,i)   = ( dnHess%d1(:,:,i) -  mat0(:,:) ) * stepp
        dnHess%d2(:,:,i,i) = (dnHess%d2(:,:,i,i)+mat0(:,:)-               &
                                   dnHess%d0(:,:)-dnHess%d0(:,:))*step2

        dnQeq%d1(:,i)      = ( dnQeq%d1(:,i) - vec0(:) ) * stepp
        dnQeq%d2(:,i,i)    = (dnQeq%d2(:,i,i)+vec0(:)-dnQeq%d0(:)-dnQeq%d0(:))*step2


        Qact(i) = vi
      END DO


!-----------------------------------------------------------------
!----- end d/Qqi and d2/dQi2 of frequencies ----------------------
!-----------------------------------------------------------------

!-----------------------------------------------------------------
!----- d2/dQidQj of frequencies (4 points) -----------------------
!      d2/dQidQj = ( v(Qi+,Qj+)+v(Qi-,Qj-)-v(Qi-,Qj+)-v(Qi+,Qj-) )/(4*s*s)
!-----------------------------------------------------------------
      DO i=1,RPHTransfo%nb_act1
      DO j=i+1,RPHTransfo%nb_act1

        vi = Qact(i)
        vj = Qact(j)


!       -- frequencies calculation at Qact(i)+step Qact(j)+step
        Qact(i) = vi + step
        Qact(j) = vj + step
        CALL sub_freq2_RPH(vec1,mat1,mat2,d1,mat0,vec0,vec2,    &
                       pot0_corgrad2,Qact,para_Tnum,mole,RPHTransfo)

        dnC_inv%d2(:,:,i,j) = mat2
        dnLnN%d2(i,j)       = d1
        dnC%d2(:,:,i,j)     = mat1(:,:)
        dnHess%d2(:,:,i,j)  = mat0(:,:)
        dnQeq%d2(:,i,j)     = vec0(:)

!       -- frequencies calculation at Qact(i)-step Qact(j)-step
        Qact(i) = vi - step
        Qact(j) = vj - step

        CALL sub_freq2_RPH(vec1,mat1,mat2,d1,mat0,vec0,vec2,    &
                       pot0_corgrad2,Qact,para_Tnum,mole,RPHTransfo)

        dnC_inv%d2(:,:,i,j) = dnC_inv%d2(:,:,i,j) + mat2
        dnLnN%d2(i,j)       = dnLnN%d2(i,j)       + d1
        dnC%d2(:,:,i,j)     = dnC%d2(:,:,i,j)     + mat1
        dnHess%d2(:,:,i,j)  = dnHess%d2(:,:,i,j)  + mat0
        dnQeq%d2(:,i,j)     = dnQeq%d2(:,i,j)     + vec0

!       -- frequencies calculation at Qact(i)-step Qact(j)+step
        Qact(i) = vi - step
        Qact(j) = vj + step

        CALL sub_freq2_RPH(vec1,mat1,mat2,d1,mat0,vec0,vec2,    &
                       pot0_corgrad2,Qact,para_Tnum,mole,RPHTransfo)

        dnC_inv%d2(:,:,i,j) = dnC_inv%d2(:,:,i,j) - mat2
        dnLnN%d2(i,j)       = dnLnN%d2(i,j)       - d1
        dnC%d2(:,:,i,j)     = dnC%d2(:,:,i,j)     - mat1
        dnHess%d2(:,:,i,j)  = dnHess%d2(:,:,i,j)  - mat0
        dnQeq%d2(:,i,j)     = dnQeq%d2(:,i,j)     - vec0

!       -- frequencies calculation at Qact(i)+step Qact(j)-step
        Qact(i) = vi + step
        Qact(j) = vj - step

        CALL sub_freq2_RPH(vec1,mat1,mat2,d1,mat0,vec0,vec2,    &
                       pot0_corgrad2,Qact,para_Tnum,mole,RPHTransfo)

        dnC_inv%d2(:,:,i,j) = dnC_inv%d2(:,:,i,j) - mat2
        dnLnN%d2(i,j)       = dnLnN%d2(i,j)       - d1
        dnC%d2(:,:,i,j)     = dnC%d2(:,:,i,j)     - mat1
        dnHess%d2(:,:,i,j)  = dnHess%d2(:,:,i,j)  - mat0
        dnQeq%d2(:,i,j)     = dnQeq%d2(:,i,j)     - vec0


!       -- d2/dQi/dQj -----------------------------------------

        dnC_inv%d2(:,:,i,j) = dnC_inv%d2(:,:,i,j) * step24
        dnC_inv%d2(:,:,j,i) = dnC_inv%d2(:,:,i,j)

        dnLnN%d2(i,j)       = dnLnN%d2(i,j)       * step24
        dnLnN%d2(j,i)       = dnLnN%d2(i,j)

        dnC%d2(:,:,i,j)     = dnC%d2(:,:,i,j)     * step24
        dnC%d2(:,:,j,i)     = dnC%d2(:,:,i,j)

        dnHess%d2(:,:,i,j)  = dnHess%d2(:,:,i,j)  * step24
        dnHess%d2(:,:,j,i)  = dnHess%d2(:,:,i,j)

        dnQeq%d2(:,i,j)     = dnQeq%d2(:,i,j)     * step24
        dnQeq%d2(:,j,i)     = dnQeq%d2(:,i,j)


        Qact(i) = vi
        Qact(j) = vj
      END DO
      END DO
!-----------------------------------------------------------------
!----- end d2/dQidQj of frequencies ------------------------------
!-----------------------------------------------------------------


      DO i=1,RPHTransfo%nb_act1
        dnLnN%d1(i) = dnLnN%d1(i)/dnLnN%d0
      END DO
      DO i=1,RPHTransfo%nb_act1
      DO j=1,RPHTransfo%nb_act1
        dnLnN%d2(i,j) = dnLnN%d2(i,j)/dnLnN%d0 - dnLnN%d1(i)*dnLnN%d1(j)
      END DO
      END DO

!-----------------------------------------------------------
       IF (print_level == 1) write(out_unitp,11)                        &
                                            Qact(1:RPHTransfo%nb_act1), &
                                            dnEHess%d0(:)*auTOcm_inv
 11    format(' frequencies : ',30f10.4)

       IF (print_level > 1) write(out_unitp,*) ' frequencies : ',       &
                                            Qact(1:RPHTransfo%nb_act1), &
                                            dnEHess%d0(:)*auTOcm_inv

       IF (debug .OR. test) THEN
         write(out_unitp,*) 'dnQeq'
         CALL Write_dnVec(dnQeq)
         write(out_unitp,*) 'dnHess'
         CALL Write_dnMat(dnHess)
         write(out_unitp,*) 'dnC_inv'
         CALL Write_dnMat(dnC_inv)
       END IF

       IF (debug) THEN
         write(out_unitp,*) 'END ',name_sub
       END IF
      CALL flush_perso(out_unitp)
!-----------------------------------------------------------

      END SUBROUTINE sub_dnfreq_4p
!=============================================================
!
!     frequency calculations along Qact
!
!=============================================================
      SUBROUTINE sub_freq2_RPH(d0ehess,d0c,d0c_inv,                     &
                               norme,d0hess,d0Qeq,d0g,pot0_corgrad,     &
                               Qact,para_Tnum,mole,RPHTransfo)
      USE mod_system
      USE mod_dnSVM
      USE mod_Constant, only : get_Conv_au_TO_unit
      USE mod_Coord_KEO
      IMPLICIT NONE

      !----- for the zmatrix and Tnum --------------------------------------
      TYPE (Tnum)        :: para_Tnum
      TYPE (zmatrix)     :: mole

      real (kind=Rkind), intent(inout) :: Qact(:)

      real (kind=Rkind)  :: rho,vep

!----- variables for the active and inactive namelists ----------------
      TYPE (Type_RPHTransfo)  :: RPHTransfo

!-------------------------------------------------------------------------
!-------------------------------------------------------------------------



!------ pour les frequences -------------------------------

       real (kind=Rkind) :: d0ehess(RPHTransfo%nb_inact21)
       real (kind=Rkind) :: d0ek(RPHTransfo%nb_inact21)
       real (kind=Rkind) :: d0g(RPHTransfo%nb_inact21)
       real (kind=Rkind) :: d0Qeq(RPHTransfo%nb_inact21)
       real (kind=Rkind) :: d0hess(RPHTransfo%nb_inact21,RPHTransfo%nb_inact21)
       real (kind=Rkind) :: d0c(RPHTransfo%nb_inact21,RPHTransfo%nb_inact21)
       real (kind=Rkind) :: d0c_inv(RPHTransfo%nb_inact21,RPHTransfo%nb_inact21)
       real (kind=Rkind) :: norme

       real (kind=Rkind) :: pot0_corgrad


!----- working variables ---------------------------------------------
!----- variables pour les derivees -----------------------------------
      logical       :: deriv,num
      integer       :: i,j,i_Qdyn,nderiv
      logical       :: special

      integer       :: nb_act1,nb_inact21
      real (kind=Rkind) :: Qdyn(mole%nb_var)


      real (kind=Rkind) :: a,d0req
      real (kind=Rkind) :: auTOcm_inv

      real (kind=Rkind), allocatable ::                                 &
        d1req(:),d2req(:,:),d3req(:,:,:),                               &
        d1g(:,:),d2g(:,:,:),                                            &
        d0h(:,:),d1hess(:,:,:),d2hess(:,:,:,:),                         &
        d0k(:,:),                                                       &
        d0hess_inv(:,:),trav1(:),NonDiag_Scaling(:)

      TYPE(Type_dnMat) :: dnGG

!----- for debuging --------------------------------------------------
       integer :: err_mem,memory
       character (len=*), parameter :: name_sub = 'sub_freq2_RPH'
       logical, parameter :: debug = .FALSE.
       !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'Qact',Qact
         write(out_unitp,*) 'RPHTransfo%step',RPHTransfo%step
         CALL flush_perso(out_unitp)
       END IF
!-----------------------------------------------------------
      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')

      nb_act1    = RPHTransfo%nb_act1
      nb_inact21 = RPHTransfo%nb_inact21
      IF (.NOT. associated(RPHTransfo%C_ini)) THEN
        CALL alloc_array(RPHTransfo%C_ini,(/nb_inact21,nb_inact21/),    &
                          "RPHTransfo%C_ini",name_sub)
        RPHTransfo%C_ini(:,:) = ZERO
      END IF
      IF (debug) THEN
        write(out_unitp,*) 'RPHTransfo%C_ini'
        CALL Write_Mat(RPHTransfo%C_ini,out_unitp,4)
        CALL flush_perso(out_unitp)
      END IF

!-----------------------------------------------------------------
!--------- Qact => Qdyn ------------------------------------------
! we need Qdyn because, we calculate, the hessian, grandient with Qdyn coord
!-----------------------------------------------------------------
       CALL Qact_TO_Qdyn_FROM_ActiveTransfo(Qact,Qdyn,mole%ActiveTransfo)
!-----------------------------------------------------------------


!-----------------------------------------------------------------
!-----------------------------------------------------------------
!----- d0Qeq, d0g and d0h at Qdyn --------------------------------
!      Qdyn and Qact are also modified
!-----------------------------------------------------------------
      nderiv = 0
      deriv  = .FALSE.

      CALL alloc_NParray(d1req,(/ nb_act1 /),"d1req",name_sub)
      CALL alloc_NParray(d2req,(/ nb_act1,nb_act1 /),"d2req",name_sub)
      CALL alloc_NParray(d3req,(/ nb_act1,nb_act1,nb_act1 /),"d3req",name_sub)

      DO i=1,nb_inact21

        i_Qdyn = mole%ActiveTransfo%list_QactTOQdyn(nb_act1+i)

        CALL d0d1d2d3_Qeq(i_Qdyn,d0req,d1req,d2req,d3req,Qdyn,mole,nderiv)

        IF (debug) write(out_unitp,*) 'i_Qdyn,i,d0req',i_Qdyn,i,d0req
        CALL flush_perso(out_unitp)
        d0Qeq(i)             = d0req
        Qdyn(i_Qdyn)         = d0req
        Qact(nb_act1+i)      = d0req

      END DO

      CALL dealloc_NParray(d1req,"d1req",name_sub)
      CALL dealloc_NParray(d2req,"d2req",name_sub)
      CALL dealloc_NParray(d3req,"d3req",name_sub)

      !------ The gradient ----------------------------------
      CALL alloc_NParray(d1g,(/ nb_inact21,nb_act1 /),"d1g",name_sub)
      CALL alloc_NParray(d2g,(/ nb_inact21,nb_act1,nb_act1 /),"d2g",name_sub)

      CALL d0d1d2_g(d0g,d1g,d2g,Qdyn,mole,.FALSE.,.FALSE.,RPHTransfo%step)

      CALL dealloc_NParray(d1g,"d1g",name_sub)
      CALL dealloc_NParray(d2g,"d2g",name_sub)

      !------ The hessian ----------------------------------
      CALL alloc_NParray(d1hess,(/ nb_inact21,nb_inact21,nb_act1 /),    &
                        "d1hess",name_sub)
      CALL alloc_NParray(d2hess,(/nb_inact21,nb_inact21,nb_act1,nb_act1/),&
                        "d2hess",name_sub)

      CALL d0d1d2_h(d0hess,d1hess,d2hess,Qdyn,mole,.FALSE.,.FALSE.,RPHTransfo%step)

      CALL dealloc_NParray(d1hess,"d1hess",name_sub)
      CALL dealloc_NParray(d2hess,"d2hess",name_sub)


      !-----------------------------------------------------------------
      !- the gardient is taken into account for d0Qeq -------------
      CALL alloc_NParray(d0h,(/nb_inact21,nb_inact21/),"d0h",name_sub)
      d0h(:,:) = d0hess(:,:)

      IF (RPHTransfo%gradTOpot0) THEN
        CALL alloc_NParray(d0hess_inv,(/nb_inact21,nb_inact21/),"d0hess_inv",name_sub)
        CALL alloc_NParray(trav1,(/nb_inact21/),"trav1",name_sub)

        CALL inv_m1_TO_m2(d0h,d0hess_inv,nb_inact21,0,ZERO) ! not SVD
        trav1(:)     = matmul(d0hess_inv,d0g)
        pot0_corgrad = -HALF*dot_product(d0g,trav1)
        d0Qeq(:)     = d0Qeq(:) - trav1(:)
        d0g(:)       = ZERO

        !-- merge d0Qeq(:) with Qact(:)
        CALL Qinact2n_TO_Qact_FROM_ActiveTransfo(d0Qeq,Qact,mole%ActiveTransfo)

        CALL dealloc_NParray(d0hess_inv,"d0hess_inv",name_sub)
        CALL dealloc_NParray(trav1,"trav1",name_sub)
      ELSE
        pot0_corgrad = ZERO
      END IF

      !-----------------------------------------------------------------
      !------ The kinetic part -------------------------------
      CALL alloc_NParray(d0k,(/nb_inact21,nb_inact21/),"d0k",name_sub)

      CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,0)

      CALL get_dng_dnGG(Qact,para_Tnum,mole,dnGG=dnGG,nderiv=0)

      d0k(:,:) = dnGG%d0(nb_act1+1:nb_act1+nb_inact21,                  &
                         nb_act1+1:nb_act1+nb_inact21)

      CALL dealloc_dnSVM(dnGG)


!-----------------------------------------------------------------
!     --- frequencies and normal modes calculation ....
      IF (debug) THEN
        write(out_unitp,*) 'd0hess,d0k'
        CALL Write_Mat(d0hess,out_unitp,4)
        CALL Write_Mat(d0k,out_unitp,4)
        CALL flush_perso(out_unitp)
      END IF

      d0h(:,:) = d0hess(:,:)

      IF (RPHTransfo%purify_hess .OR. RPHTransfo%eq_hess) THEN
        CALL H0_symmetrization(d0h,nb_inact21,                        &
                               RPHTransfo%Qinact2n_sym,               &
                               RPHTransfo%dim_equi,RPHTransfo%tab_equi)
        CALL H0_symmetrization(d0k,nb_inact21,                        &
                               RPHTransfo%Qinact2n_sym,               &
                               RPHTransfo%dim_equi,RPHTransfo%tab_equi)
        IF (debug) THEN
          write(out_unitp,*) 'sym : d0hess,d0k'
          CALL Write_Mat(d0h,out_unitp,4)
          CALL Write_Mat(d0k,out_unitp,4)
        END IF

        CALL calc_freq_block(nb_inact21,d0h,d0k,d0ehess,                &
                             d0c,d0c_inv,norme,RPHTransfo%C_ini,        &
                             RPHTransfo%diabatic_freq,RPHTransfo%Qinact2n_sym)

        !CALL calc_freq(nb_inact21,d0h,d0k,d0ehess,d0c,d0c_inv,norme,    &
        !               RPHTransfo%C_ini,RPHTransfo%diabatic_freq)

      ELSE
        CALL calc_freq(nb_inact21,d0h,d0k,d0ehess,d0c,d0c_inv,norme,    &
                       RPHTransfo%C_ini,RPHTransfo%diabatic_freq)

      END IF

      CALL dealloc_NParray(d0h,"d0h",name_sub)
      CALL dealloc_NParray(d0k,"d0k",name_sub)

!-----------------------------------------------------------------

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'norm: ',norme
        write(out_unitp,*) 'freq : ',d0ehess(:)*auTOcm_inv
        write(out_unitp,*) 'd0c : '
        CALL Write_Mat(d0c,out_unitp,4)
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------
      END SUBROUTINE sub_freq2_RPH

      SUBROUTINE Finalyze_TnumTana_Coord_PrimOp(para_Tnum,mole,para_PES,Tana)
      USE mod_system
      USE mod_dnSVM
      USE mod_Constant, only : get_Conv_au_TO_unit
      USE mod_Coord_KEO
      USE mod_SimpleOp
      USE mod_PrimOp_def

      IMPLICIT NONE

!----- for the zmatrix and Tnum --------------------------------------
      TYPE (Tnum)        :: para_Tnum
      TYPE (zmatrix)     :: mole
      TYPE (param_PES)   :: para_PES
      logical, optional  :: Tana


!-------------------------------------------------------------------------
!-------------------------------------------------------------------------
      real (kind=Rkind) :: Qact(mole%nb_var)
      TYPE(Type_dnMat)  :: dnGG
      TYPE (Type_dnVec) :: dnx

      integer           :: NM_TO_sym_ver=4
      real (kind=Rkind), allocatable :: hCC(:,:),GGdef(:,:)
      logical :: Gref,tab_skip_transfo(mole%nb_Qtransfo),Tana_loc
      integer :: iQa,nb_act1_RPH,nb_inact21_RPH,it,nderiv
      real (kind=Rkind) :: auTOcm_inv

!----- for debuging --------------------------------------------------
       integer :: err_mem,memory
       character (len=*), parameter :: name_sub = 'Finalyze_TnumTana_Coord_PrimOp'
       logical, parameter :: debug = .FALSE.
       !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       !IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         CALL flush_perso(out_unitp)
       !END IF
!-----------------------------------------------------------

      auTOcm_inv = get_Conv_au_TO_unit('E','cm-1')

      IF (present(Tana)) THEN
        Tana_loc = Tana
      ELSE
        Tana_loc = .TRUE.
      END IF

!-----------------------------------------------------------------------
!--------------------- TO finalize the coordinates (NM) and the KEO ----

  CALL Sub_PES_FromTnum_TO_PES(para_PES,para_Tnum%para_PES_FromTnum)
  IF (para_PES%nb_scalar_Op > 0) para_PES%calc_scalar_Op = .TRUE.

  !-----------------------------------------------------------------
  ! initialization of the scalar operators
  CALL Sub_init_dnOp(mole,para_Tnum,para_PES)
  !-----------------------------------------------------------------


  !----- calc and transfert NM to LinearTransfo%mat if needed ---------------
  IF (associated(mole%NMTransfo)) THEN
      IF (.NOT. mole%tab_Qtransfo(mole%itNM)%skip_transfo) THEN

        CALL get_Qact0(Qact,mole%ActiveTransfo)

        IF (NM_TO_sym_ver == 3) THEN

          CALL alloc_NParray(hCC,(/ mole%ncart_act,mole%ncart_act /),"hCC",name_sub)
          CALL calc3_NM_TO_sym(Qact,mole,para_Tnum,para_PES,hCC,.FALSE.)
          CALL dealloc_NParray(hCC,"hCC",name_sub)

        ELSE
          CALL calc4_NM_TO_sym(Qact,mole,para_Tnum,para_PES)
        END IF
        IF (print_level > 1) CALL sub_QplusDQ_TO_Cart(Qact,mole)

      END IF
  END IF

  !----- set RPH transfo of Qref -----------------------------------
  IF (associated(mole%RPHTransfo)) THEN
    IF (.NOT. mole%tab_Qtransfo(mole%itRPH)%skip_transfo .AND.          &
                                       mole%RPHTransfo%option /= 0) THEN

      CALL get_Qact0(Qact,mole%ActiveTransfo)

      !we cannot set RPHpara_AT_Qref for option=0,
      ! because the RPHTransfo parameters will be readed in "inactive" namelist

      DO it=mole%nb_Qtransfo-1,mole%itRPH+1,-1
          tab_skip_transfo(it) = mole%tab_Qtransfo(it)%skip_transfo
          mole%tab_Qtransfo(it)%skip_transfo = .TRUE.
      END DO

        ! for RPHpara_AT_Qref
        IF (.NOT. associated(mole%RPHTransfo%RPHpara_AT_Qref)) THEN
          nb_act1_RPH    = mole%RPHTransfo%nb_act1
          nb_inact21_RPH = mole%RPHTransfo%nb_inact21
          CALL alloc_array(mole%RPHTransfo%RPHpara_AT_Qref,(/ 1 /),       &
                          'mole%RPHTransfo%RPHpara_AT_Qref',name_sub)
        END IF

        CALL Set_RPHpara_AT_Qact1(mole%RPHTransfo%RPHpara_AT_Qref(1),     &
                                  Qact,para_Tnum,mole,mole%RPHTransfo)

        mole%RPHTransfo%init_Qref = .TRUE.

        DO it=mole%nb_Qtransfo-1,mole%itRPH+1,-1
          mole%tab_Qtransfo(it)%skip_transfo = tab_skip_transfo(it)
        END DO


!        ! new Qact0/Qdyn0
!        ! first Delta_Qdyn0
!        mole%ActiveTransfo%Qdyn0(mole%RPHTransfo%nb_act1+1:                        &
!                             mole%RPHTransfo%nb_act1+mole%RPHTransfo%nb_inact21) = &
!                               mole%ActiveTransfo%Qdyn0(mole%RPHTransfo%nb_act1+1: &
!                             mole%RPHTransfo%nb_act1+mole%RPHTransfo%nb_inact21) - &
!                                     mole%RPHTransfo%RPHpara_AT_Qref(1)%dnQopt%d0
!     write(6,*) 'delta, mole%ActiveTransfo%Qdyn0',mole%ActiveTransfo%Qdyn0 ; flush(6)
!        !Qdyn0 = matmul(Delta_Qdyn0, ...dnC%d0)
!        mole%ActiveTransfo%Qdyn0(mole%RPHTransfo%nb_act1+1:                        &
!                             mole%RPHTransfo%nb_act1+mole%RPHTransfo%nb_inact21) = &
!                    matmul(mole%ActiveTransfo%Qdyn0(mole%RPHTransfo%nb_act1+1:     &
!                             mole%RPHTransfo%nb_act1+mole%RPHTransfo%nb_inact21),  &
!                             mole%RPHTransfo%RPHpara_AT_Qref(1)%dnC%d0)
!
!     write(6,*) 'mole%ActiveTransfo%Qdyn0',mole%ActiveTransfo%Qdyn0 ; flush(6)


        CALL Qdyn_TO_Qact_FROM_ActiveTransfo(mole%ActiveTransfo%Qdyn0,  &
                                             mole%ActiveTransfo%Qact0,  &
                                             mole%ActiveTransfo)

        write(out_unitp,*) 'New Qact0',mole%ActiveTransfo%Qact0

        write(out_unitp,*) ' Frequencies, normal modes at the reference geometry'

        write(out_unitp,11) Qact(1:mole%RPHTransfo%RPHpara_AT_Qref(1)%nb_act1), &
               mole%RPHTransfo%RPHpara_AT_Qref(1)%dnEHess%d0(:)*auTOcm_inv
 11     format(' frequencies : ',30f10.4)

        write(out_unitp,*) 'dnQopt'
        CALL Write_dnVec(mole%RPHTransfo%RPHpara_AT_Qref(1)%dnQopt,nderiv=0)
        write(out_unitp,*) 'dnC_inv'
        CALL Write_dnMat(mole%RPHTransfo%RPHpara_AT_Qref(1)%dnC_inv,nderiv=0)
        CALL flush_perso(out_unitp)


       IF (debug) CALL Write_RPHTransfo(mole%RPHTransfo)

    END IF
  END IF

  !----- Gcte if needed --------------------------------------------
  Gref = .TRUE.
  IF (associated(mole%RPHTransfo)) THEN
    Gref = Gref .AND. associated(mole%RPHTransfo%RPHpara_AT_Qref)
  END IF

  IF (Gref) THEN
    CALL get_Qact0(Qact,mole%ActiveTransfo)
    IF (print_level > 1) write(out_unitp,*) ' para_Tnum%Gcte'

    CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,0)

    IF (para_PES%QMLib) THEN
      CALL alloc_NPArray(GGdef,(/ mole%nb_act,mole%nb_act /),'GGdef',name_sub)

#if __QML == 1
      CALL get_Qmodel_GGdef(GGdef)
#else
      write(out_unitp,*) 'ERROR in ',name_sub
      write(out_unitp,*) ' The "Quantum Model Lib" (QML) library is not present!'
      write(out_unitp,*) 'Use another potential/model'
      STOP 'QML is not present'
#endif

      dnGG%d0(:,:) = ZERO
      dnGG%d0(1:mole%nb_act,1:mole%nb_act) = GGdef

      CALL dealloc_NPArray(GGdef,'GGdef',name_sub)
    ELSE
      CALL get_dng_dnGG(Qact,para_Tnum,mole,dnGG=dnGG,nderiv=0)
    END IF
    IF (para_Tnum%Gcte) THEN

      CALL alloc_array(para_Tnum%Gref,(/ mole%ndimG,mole%ndimG /),    &
                      'para_Tnum%Gref',name_sub)

      para_Tnum%Gref(:,:) = dnGG%d0(:,:)
    END IF

    IF (print_level > 1) THEN
      write(out_unitp,*) ' dnGG%d0'
      CALL Write_Mat(dnGG%d0,out_unitp,5)
    END IF

    CALL dealloc_dnSVM(dnGG)
  END IF

      !----- Tana if needed --------------------------------------------
      IF (para_Tnum%Tana .AND. Tana_loc) THEN
        write(out_unitp,*)
        write(out_unitp,*) '================================================='
        write(out_unitp,*) ' BEGINNING Tana'
        CALL time_perso('Tana')
        write(out_unitp,*) 'nrho_OF_Qdyn(:)',mole%nrho_OF_Qdyn(:)
        write(out_unitp,*) 'nrho_OF_Qact(:)',mole%nrho_OF_Qact(:)

        !IF (print_level > 1) write(out_unitp,*) ' para_Tnum%Tana'
        CALL compute_analytical_KEO(para_Tnum%TWOxKEO,mole,para_Tnum,Qact)
        IF (debug) CALL write_op(para_Tnum%TWOxKEO,header=.TRUE.)

        CALL comparison_G_FROM_Tnum_Tana(para_Tnum%ExpandTWOxKEO,mole,para_Tnum,Qact)

        write(out_unitp,*)
        CALL time_perso('Tana')
        write(out_unitp,*) ' END Tana'
        write(out_unitp,*) '================================================='

      END IF

      write(out_unitp,*) '================================================='
      write(out_unitp,*) '=== Reference geometry (not recenter) ==========='
      CALL alloc_dnSVM(dnx,mole%ncart,mole%nb_act,nderiv=0)

      CALL get_Qact0(Qact,mole%ActiveTransfo)
      CALL sub_QactTOdnx(Qact,dnx,mole,nderiv=0,Gcenter=.FALSE.,WriteCC=.TRUE.)

      CALL dealloc_dnSVM(dnx)
      write(out_unitp,*) '================================================='

!-----------------------------------------------------------
      !IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      !END IF
!-----------------------------------------------------------

      END SUBROUTINE Finalyze_TnumTana_Coord_PrimOp


   END MODULE mod_PrimOp

