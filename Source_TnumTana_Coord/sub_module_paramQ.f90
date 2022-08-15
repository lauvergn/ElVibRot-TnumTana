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
!      Tnum is written David Lauvergnat [1]
!      Tana is written by Mamadou Ndong [1] and David Lauvergnat [1]
!         with contributions
!          Emil Lund klinting (coupling with MidasCpp) [3]'
!
![1]: Institut de Chimie Physique, UMR 8000, CNRS-Université Paris-Saclay, France
![3]: Department of Chemistry, Aarhus University, DK-8000 Aarhus C, Denmark
!
!===========================================================================
!===========================================================================
MODULE mod_paramQ
  use mod_system
  USE mod_dnSVM
  use mod_Lib_QTransfo,     only: write_cart, calc_cross_product,       &
                                  write_dnx, sub3_dnx_at1
  use mod_ActiveTransfo,    only: qact_to_qdyn_from_activetransfo
  use mod_CartesianTransfo, only: alloc_cartesiantransfo,               &
           Set_P_Axis_CartesianTransfo, Set_Eckart_CartesianTransfo,    &
                               centre_masse, write_cartesiantransfo,    &
                             sub_dnxmassweight, sub3_dncentre_masse,    &
                     calc_cartesiantransfo_new, sub_dnxnomassweight,    &
                                              sub3_nodncentre_masse
  use mod_Qtransfo,         only: write_Qtransfo, calc_Qtransfo
  use mod_Tnum,             only: tnum, zmatrix, write_mole,            &
                                  CoordType, Write_CoordType

  USE mod_Constant,         ONLY: real_wu,get_conv_au_to_unit,          &
                                  rwu_write,convRWU_TO_R_WITH_WorkingUnit
  IMPLICIT NONE

  INTERFACE sub_QactTOdnx
    MODULE PROCEDURE sub_QactTOdnx_zmatrix,sub_QactTOdnx_CoordType
  END INTERFACE
  INTERFACE sub_QactTOd0x
    MODULE PROCEDURE sub_QactTOd0x_zmatrix,sub_QactTOd0x_CoordType
  END INTERFACE

  INTERFACE Write_Cartg98
    MODULE PROCEDURE Write_Cartg98_zmatrix,Write_Cartg98_CoordType
  END INTERFACE
  INTERFACE analyze_dnx
    MODULE PROCEDURE analyze_dnx_zmatrix,analyze_dnx_CoordType
  END INTERFACE
  INTERFACE read_RefGeom
    MODULE PROCEDURE read_RefGeom_CoordType,read_RefGeom_zmatrix
  END INTERFACE


  PRIVATE
  PUBLIC :: read_RefGeom, Get_Qread
  PUBLIC :: sub_QactTOQit, sub_QinRead_TO_Qact
  PUBLIC :: sub_QplusDQ_TO_Cart, sub_QactTOdnMWx, sub_QactTOdnx, sub_QactTOd0x, sub_d0xTOQact
  PUBLIC :: Write_d0Q, Write_Q_WU, Write_Cartg98, Write_XYZ
  PUBLIC :: analyze_dnx, sub_dnFCC_TO_dnFcurvi, write_dnx
  PUBLIC :: Set_paramQ_FOR_optimization

CONTAINS

!=======================================================================================
!      Read reference geometry
!=======================================================================================
      SUBROUTINE read_RefGeom_zmatrix(mole,para_Tnum)
      USE mod_MPI
      IMPLICIT NONE

!----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix)     :: mole
      TYPE (Tnum)        :: para_Tnum

      CALL read_RefGeom_CoordType(mole%CoordType,para_Tnum)

      END SUBROUTINE read_RefGeom_zmatrix
      SUBROUTINE read_RefGeom_CoordType(mole,para_Tnum)
      USE mod_MPI
      IMPLICIT NONE

!----- for the zmatrix and Tnum --------------------------------------
      TYPE (CoordType)   :: mole
      TYPE (Tnum)        :: para_Tnum

      logical            :: read_Qact0,read_Qdyn0,read_Qsym0
      logical            :: read_xyz0,read_xyz0_with_dummy
      logical            :: read_nameQ
      integer            :: read_itQ0transfo,read_itQtransfo_OF_Qin0
      character (len=Name_len) :: name,unit



!----- The coordinates which are read --------------------------------
      character (len=Name_len), pointer :: QnameRead(:)
      real(kind=Rkind), allocatable     :: Qread(:)
      real(kind=Rkind)                  :: Qact(mole%nb_var)
      real(kind=Rkind)                  :: Qdyn(mole%nb_var)


      logical           :: opt,pot_act,pot_cart,HarD,pot_cplx,OnTheFly
      TYPE(Type_dnVec)  :: dnx

      integer           :: pot_itQtransfo
      !real (kind=Rkind) :: pot0
      TYPE (REAL_WU)    :: pot0

      integer           :: nb_elec,nb_scalar_Op,nb_CAP

      logical                  :: deriv_WITH_FiniteDiff  = .FALSE.
      logical                  :: nDfit_Op               = .FALSE.
      logical                  :: QMLib                  = .FALSE.
      character (len=Line_len) :: BaseName_nDfit_file

      real (kind=Rkind) :: ex(3),nx,ey(3),ny,ez(3),nz
      real (kind=Rkind) :: VG(3)


      integer           :: iref,i,nb_t,type_Qin,type_Qread,nc1,nc2,nc3
      character (len=Line_len) :: info_Qread

      !-----------------------------------------------------------------

      NAMELIST /minimum/ read_itQ0transfo,read_Qsym0,read_Qdyn0,read_xyz0, &
                         read_nameQ,unit,read_xyz0_with_dummy,read_Qact0,  &
                         pot0,pot_act,pot_cart,pot_itQtransfo,pot_cplx,    &
                         HarD,nb_elec,nb_scalar_Op,nb_CAP,                 &
                         OnTheFly,nDfit_Op,QMLib,BaseName_nDfit_file,      &
                         deriv_WITH_FiniteDiff,                            &
                         opt

      !-----------------------------------------------------------------
      integer :: err_mem,memory,err_io
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='read_RefGeom'
      !-----------------------------------------------------------------


      write(out_unitp,*) 'BEGINNING ',name_sub

!------- allocation of Q... -----------------------------------------


      IF (.NOT. associated(mole%name_Qdyn))  THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' mole%name_Qdyn is not associated in mole!!'
        write(out_unitp,*) ' Check the source !'
        STOP
      END IF

      IF (print_level > 1) THEN
        write(out_unitp,*) '===================================='
        write(out_unitp,*) 'nb_var',mole%nb_var
      END IF

!------- read the namelist minimum -----------------------------
      read_Qsym0           = .FALSE.
      read_Qdyn0           = .FALSE.
      read_Qact0           = .FALSE.
      read_xyz0            = .FALSE.
      read_xyz0_with_dummy = .TRUE.
      read_nameQ           = .FALSE.
      read_itQ0transfo     = -1
      unit                 = 'au'


      nb_scalar_Op = 0
      nb_CAP       = 0
      pot_cplx     = .FALSE.
      nb_elec      = 0

      pot_act              = .FALSE.
      pot_cart             = .FALSE.
      pot_itQtransfo       = -1
      IF (associated(mole%RPHTransfo)) THEN
        HarD               = .FALSE.
      ELSE
        HarD               = .TRUE.
      END IF
      pot0                 = REAL_WU(ZERO,'au','E')

      nDfit_Op             = .FALSE.
      BaseName_nDfit_file  = ""
      OnTheFly             = .FALSE.
      QMLib                = .FALSE.

      deriv_WITH_FiniteDiff = .FALSE.
      opt                   = .FALSE.

      read(in_unitp,minimum,IOSTAT=err_io)
      IF (err_io < 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) '  while reading the namelist "minimum"'
        write(out_unitp,*) ' end of file or end of record'
        write(out_unitp,*) ' Probably, you have forgotten the namelist ...'
        write(out_unitp,*) ' Check your data !!'
        STOP
      END IF
      IF (err_io > 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) '  while reading the namelist "minimum"'
        write(out_unitp,*) ' Probably, some arguments of namelist are wrong.'
        write(out_unitp,*) ' Check your data !!'
        STOP
      END IF
      IF (print_level > 1) write(out_unitp,minimum)

      CALL string_uppercase_TO_lowercase(unit)

      IF (unit /= 'au' .AND. unit /= 'bohr' .AND. unit /= 'angs') THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) '  while reading the namelist "minimum"'
        write(out_unitp,*) '  The unit is wrong, unit: "',trim(unit),'"'
        write(out_unitp,*) '  The possible values are: "au" or "bohr" or "angs"'
        write(out_unitp,*) ' Check your data !!'
        STOP
      END IF


      IF (.NOT. allocated(mole%opt_Qdyn)) THEN
        CALL alloc_NParray(mole%opt_Qdyn,(/ mole%nb_var /),'mole%opt_Qdyn',name_sub)
      END IF
      mole%opt_Qdyn(:) = 0
      IF (opt) THEN
         DO i=1,mole%nb_var
          IF (mole%ActiveTransfo%list_act_OF_Qdyn(i) == 1) THEN
            mole%opt_Qdyn(i) = 1
          END IF
        END DO
      END IF

      IF(MPI_id==0) THEN
        write(out_unitp,*)  '------------------------------------------------------'
        write(out_unitp,*)  '--- Coordinates used for the operators ---------------'
        write(out_unitp,*)  '------------------------------------------------------'
      ENDIF

      IF (mole%nb_Qtransfo == -1) THEN
        pot_itQtransfo = 0
        pot_cart       = .FALSE.
        pot_act        = .FALSE.
      END IF

      IF ( (pot_act  .AND. pot_cart) .OR.                               &
           (pot_act  .AND. pot_itQtransfo /= -1) .OR.                   &
           (pot_cart .AND. pot_itQtransfo /= -1) ) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) '(pot_act=t and pot_cart=t) .OR. ...'
        write(out_unitp,*) 'pot_act',pot_act
        write(out_unitp,*) 'pot_cart',pot_cart
        write(out_unitp,*) 'pot_itQtransfo ',pot_itQtransfo
        write(out_unitp,*) ' You have to chose between these options'
        write(out_unitp,*) ' Check your data !'
        STOP
      END IF
      IF (OnTheFly) pot_itQtransfo = 0 ! cart

      IF(pot_itQtransfo == -1) THEN
        IF (pot_cart) THEN
          pot_itQtransfo = 0                      ! Qcart
        ELSE IF (pot_act) THEN
          pot_itQtransfo = mole%nb_Qtransfo       ! Qact
        ELSE
          pot_itQtransfo = mole%nb_Qtransfo-1     ! Qdyn
        END IF

      END IF
      IF (pot_itQtransfo < 0 .OR. pot_itQtransfo > mole%nb_Qtransfo) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,'(a,i0,a)') ' pot_itQtransfo is out of the range [0:',mole%nb_Qtransfo,']'
        write(out_unitp,*) ' Check your data !'
        STOP
      END IF

      IF (nb_elec < 1 .AND. .NOT. QMLib) nb_elec = 1
      para_Tnum%para_PES_FromTnum%opt            = opt
      para_Tnum%para_PES_FromTnum%pot0           = convRWU_TO_R_WITH_WorkingUnit(pot0)
      para_Tnum%para_PES_FromTnum%HarD           = HarD
      para_Tnum%para_PES_FromTnum%nb_elec        = nb_elec
      para_Tnum%para_PES_FromTnum%pot_cplx       = pot_cplx
      para_Tnum%para_PES_FromTnum%OnTheFly       = OnTheFly
      para_Tnum%para_PES_FromTnum%pot_itQtransfo = pot_itQtransfo
      para_Tnum%para_PES_FromTnum%nb_scalar_Op   = nb_scalar_Op
      para_Tnum%para_PES_FromTnum%nb_CAP         = nb_CAP

      para_Tnum%para_PES_FromTnum%deriv_WITH_FiniteDiff  = deriv_WITH_FiniteDiff

      para_Tnum%para_PES_FromTnum%nDfit_Op = nDfit_Op
      IF (nDfit_Op) THEN
       !STOP 'nDfit_Op is not possible anymore!!'
       write(out_unitp,*)  'BaseName_nDfit_file: ',trim(adjustl(BaseName_nDfit_file))

       para_Tnum%para_PES_FromTnum%BaseName_nDfit_file = BaseName_nDfit_file
       IF (len_trim(BaseName_nDfit_file) == 0) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nDfit_Op=.TRUE. and BaseName_nDfit_file is empty'
         write(out_unitp,*) ' Check your data !'
         STOP
       END IF
       para_Tnum%para_PES_FromTnum%nDFit_V_name_Fit =                   &
                           trim(adjustl(BaseName_nDfit_file)) // "-col1"

       allocate(para_Tnum%para_PES_FromTnum%nDFit_Scalar_Op_name_Fit(nb_scalar_Op))
       DO i=1,nb_scalar_Op
         para_Tnum%para_PES_FromTnum%nDFit_Scalar_Op_name_Fit(i) =      &
            trim(adjustl(BaseName_nDfit_file)) // "-col" // int_TO_char(i)
       END DO
      END IF

      para_Tnum%para_PES_FromTnum%QMLib = QMLib

      IF(MPI_id==0) THEN
        write(out_unitp,*) 'nb_scalar_Op',nb_scalar_Op
        write(out_unitp,*) 'nb_CAP',nb_CAP
        write(out_unitp,*) 'pot_itQtransfo',pot_itQtransfo

        IF (para_Tnum%para_PES_FromTnum%pot_itQtransfo == 0) THEN
          write(out_unitp,*) 'Operators (pot...) with Cartesian coordinates'
        ELSE IF (para_Tnum%para_PES_FromTnum%pot_itQtransfo == 1) THEN
          write(out_unitp,*) 'Operators (pot...) with primitive (zmat...) coordinates'
        ELSE IF (para_Tnum%para_PES_FromTnum%pot_itQtransfo == mole%nb_Qtransfo-1) THEN
          write(out_unitp,*) 'Operators (pot...) with Qdyn coordinates'
        ELSE IF (para_Tnum%para_PES_FromTnum%pot_itQtransfo == mole%nb_Qtransfo) THEN
          write(out_unitp,*) 'Operators (pot...) with Qact coordinates'
        END IF

        IF (QMLib .AND. para_Tnum%para_PES_FromTnum%pot_itQtransfo /= mole%nb_Qtransfo-1) THEN
          write(out_unitp,*)  ' WARNING QMLib=.TRUE. and its coordinates are not Qdyn!!'
        END IF

        write(out_unitp,*)  '------------------------------------------------------'
      ENDIF

      !=================================================================
      !=================================================================
      !=================================================================
      IF(MPI_id==0) THEN
        write(out_unitp,*)  '------------------------------------------------------'
        write(out_unitp,*)  '--- Coordinates used for the reference geometry ------'
        write(out_unitp,*)  '------------------------------------------------------'
      ENDIF

      ! first defined how to read the reference geometry:
      ! - with read_itQ0transfo    or
      ! - with read_Qsym0, read_xyz0, ...
      IF (read_Qsym0) read_Qdyn0 = .TRUE.

      IF ( (read_Qdyn0 .AND. read_Qact0) .OR.                           &
           (read_Qdyn0 .AND. read_xyz0) .OR.                            &
           (read_Qdyn0 .AND. read_itQ0transfo /= -1) .OR.               &
           (read_Qact0 .AND. read_xyz0) .OR.                            &
           (read_Qact0 .AND. read_itQ0transfo /= -1) .OR.               &
           (read_xyz0  .AND. read_itQ0transfo /= -1)) THEN

        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) '(read_Qdyn0=t and read_xyz0=t) .OR. ...'
        write(out_unitp,*) 'read_Qdyn0 .OR. read_Qsym0',read_Qdyn0
        write(out_unitp,*) 'read_Qact0',read_Qact0
        write(out_unitp,*) 'read_xyz0 ',read_xyz0
        write(out_unitp,*) 'read_itQ0transfo ',read_itQ0transfo
        write(out_unitp,*) ' You have to chose between these options'
        write(out_unitp,*) ' Check your data !'
        STOP
      END IF


      read_itQtransfo_OF_Qin0 = read_itQ0transfo

      IF (read_itQtransfo_OF_Qin0 == -1) THEN ! old way with read_Qsym0 or read_xyz0 ....
        IF (read_Qdyn0) THEN
          IF (print_level > 1 .OR. debug) write(out_unitp,*) ' Read Qdyn0 coordinates:'
          read_itQtransfo_OF_Qin0 = mole%nb_Qtransfo-1
        ELSE IF (read_xyz0) THEN
          IF (print_level > 1 .OR. debug) write(out_unitp,*) ' Read xyz0 coordinates:'
          read_itQtransfo_OF_Qin0 = 0
        ELSE IF (read_Qact0) THEN
          IF (print_level > 1 .OR. debug) write(out_unitp,*) ' Read Qact0 coordinates:'
          read_itQtransfo_OF_Qin0 = mole%nb_Qtransfo
        ELSE
          IF (print_level > 1 .OR. debug) write(out_unitp,*) ' Read Qdyn0 coordinates:'
          read_itQtransfo_OF_Qin0 = mole%nb_Qtransfo-1
        END IF

          !IF (print_level > 1 .OR. debug) write(out_unitp,*) ' Read Qprim0 (zmat, ...) coordinates:'
          !read_itQtransfo_OF_Qin0 = mole%itPrim
      END IF


      ! check if 0<= read_itQtransfo_OF_Qin0 <= mole%nb_Qtransfo
      IF (read_itQtransfo_OF_Qin0 < 0 .OR.                              &
          read_itQtransfo_OF_Qin0 > mole%nb_Qtransfo) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' read_itQ0transfo (or read_itQtransfo_OF_Qin0) ...'
        write(out_unitp,'(a,i0,a)') '... is out of the range [0:',mole%nb_Qtransfo,']'
        write(out_unitp,*) ' Check your data !'
        STOP
      END IF
      read_xyz0 = (read_itQtransfo_OF_Qin0 == 0)
      IF (print_level > 1 .OR. debug) write(out_unitp,*) 'read_itQtransfo_OF_Qin0',read_itQtransfo_OF_Qin0
      CALL flush_perso(out_unitp)
      ! defined the "info" from read_itQtransfo_OF_Qin0
      IF (read_itQtransfo_OF_Qin0 == mole%nb_Qtransfo) THEN ! Qact
        info_Qread = ' Read Qact0 coordinates:'
      ELSE IF (read_itQtransfo_OF_Qin0 == mole%nb_Qtransfo-1) THEN ! Qdyn
        info_Qread = ' Read Qdyn0 coordinates:'
      ELSE IF (read_itQtransfo_OF_Qin0 == mole%itPrim) THEN ! Qprim (zmat, poly ...)
        info_Qread = ' Read Qprim0 coordinates:'
      ELSE IF (read_itQtransfo_OF_Qin0 == 0) THEN ! Qprim (zmat, poly ...)
        info_Qread = ' Read xyz0 coordinates:'
      ELSE
        info_Qread = ' Read Qin0 coordinates, from transfo_it' //       &
                            int_TO_char(read_itQtransfo_OF_Qin0) // ':'
      END IF
      IF(MPI_id==0) write(out_unitp,*) info_Qread
      ! ----------------------------------------------

      ! ----------------------------------------------
      ! read the coordinates + conversion (angs,deg => bohr, radian)
      IF (read_itQtransfo_OF_Qin0 == 0) THEN ! special case for Cartesian coordinates

        IF (read_xyz0_with_dummy) THEN
          CALL alloc_NParray(Qread,(/ mole%tab_Qtransfo(1)%nb_Qout /),'Qread',name_sub)
          Qread(:) = ZERO

          CALL Get_Qread(Qread(1:3*mole%nat0),                        &
                         mole%tab_Qtransfo(1)%name_Qout,              &
                         mole%tab_Qtransfo(1)%type_Qout,              &
                         read_nameQ,unit,read_xyz0,info=info_Qread)


        ELSE
          CALL alloc_NParray(Qread,(/ mole%tab_Qtransfo(1)%nb_Qout /),'Qread',name_sub)
          Qread(:) = ZERO

          CALL Get_Qread(Qread(1:3*mole%nat_act),                     &
                         mole%tab_Qtransfo(1)%name_Qout,              &
                         mole%tab_Qtransfo(1)%type_Qout,              &
                         read_nameQ,unit,read_xyz0,info=info_Qread)
        END IF

      ELSE

        CALL alloc_NParray(Qread,(/mole%tab_Qtransfo(read_itQtransfo_OF_Qin0)%nb_Qin /),'Qread',name_sub)

        CALL Get_Qread(Qread,                                         &
                mole%tab_Qtransfo(read_itQtransfo_OF_Qin0)%name_Qin,  &
                mole%tab_Qtransfo(read_itQtransfo_OF_Qin0)%type_Qin,  &
                       read_nameQ,unit,read_xyz0,info=info_Qread)

      END IF
      ! ----------------------------------------------

      CALL sub_QinRead_TO_Qact(Qread,Qact,mole,read_itQtransfo_OF_Qin0)
      CALL Qact_TO_Qdyn_FROM_ActiveTransfo(Qact,Qdyn,mole%ActiveTransfo)

 111  format(a,1x,f15.6)
      IF(MPI_id==0) write(out_unitp,*) 'Qdyn0 coordinates (not transformed): [bohr]/[rad or cos(angle)]'
      DO i=1,mole%nb_var
        name = mole%tab_Qtransfo(mole%nb_Qtransfo)%name_Qout(i)
        IF(MPI_id==0) write(out_unitp,111) name,Qdyn(i)
      END DO

      IF(MPI_id==0) write(out_unitp,*) 'Qact0 coordinates (not transformed): [bohr]/[rad or cos(angle)]'
      DO i=1,mole%nb_var
        name = mole%tab_Qtransfo(mole%nb_Qtransfo)%name_Qin(i)
        IF(MPI_id==0) write(out_unitp,111) name,Qact(i)
      END DO

      IF(MPI_id==0) CALL Write_Q_WU(Qdyn,                                              &
                                    mole%tab_Qtransfo(mole%nb_Qtransfo)%name_Qout,     &
                                    mole%tab_Qtransfo(mole%nb_Qtransfo)%type_Qout,     &
                                    info='Qdyn0 coordinates (transformed):')

      IF(MPI_id==0) CALL Write_Q_WU(Qact,                                              &
                                    mole%tab_Qtransfo(mole%nb_Qtransfo)%name_Qin,      &
                                    mole%tab_Qtransfo(mole%nb_Qtransfo)%type_Qin,      &
                                    info='Qact0 coordinates (transformed):')

      ! Transfert the rigid values in ActiveTransfo%Qdyn0 and ActiveTransfo%Qact0
      mole%tab_Qtransfo(mole%nb_Qtransfo)%ActiveTransfo%Qdyn0(:) =  Qdyn(:)
      mole%tab_Qtransfo(mole%nb_Qtransfo)%ActiveTransfo%Qact0(:) =  Qact(:)
      mole%tab_Qtransfo(mole%nb_Qtransfo)%print_done = .FALSE.
      IF (debug) CALL Write_Qtransfo(mole%tab_Qtransfo(mole%nb_Qtransfo))
      ! END Transfert the rigid values in ActiveTransfo%Qdyn0

      write(out_unitp,*)  '------------------------------------------------------'
      !=================================================================
      !=================================================================
      !=================================================================


!======================================================================
!     IF Cart_transfo=t
!======================================================================
      IF (mole%tab_Qtransfo(1)%name_transfo == 'zmat'  .AND.            &
              mole%tab_Qtransfo(1)%ZmatTransfo%New_Orient .AND.         &
         sum(abs(mole%tab_Qtransfo(1)%ZmatTransfo%vAt1)) == ZERO .AND.  &
         sum(abs(mole%tab_Qtransfo(1)%ZmatTransfo%vAt2)) == ZERO .AND.  &
         sum(abs(mole%tab_Qtransfo(1)%ZmatTransfo%vAt3)) == ZERO .AND.  &
         read_xyz0) THEN

        nc1 = mole%tab_Qtransfo(1)%ZmatTransfo%ind_zmat(1,1)
        nc2 = mole%tab_Qtransfo(1)%ZmatTransfo%ind_zmat(1,2)
        nc3 = mole%tab_Qtransfo(1)%ZmatTransfo%ind_zmat(1,3)

        mole%tab_Qtransfo(1)%ZmatTransfo%vAt1(:) = Qread(nc1:nc1+2)
        mole%tab_Qtransfo(1)%ZmatTransfo%vAt2(:) = Qread(nc2:nc2+2)
        mole%tab_Qtransfo(1)%ZmatTransfo%vAt3(:) = Qread(nc3:nc3+2)
        write(out_unitp,*) 'vAt1', mole%tab_Qtransfo(1)%ZmatTransfo%vAt1
        write(out_unitp,*) 'vAt2', mole%tab_Qtransfo(1)%ZmatTransfo%vAt2
        write(out_unitp,*) 'vAt3', mole%tab_Qtransfo(1)%ZmatTransfo%vAt3
      END IF


      ! Transfert/calculate the reference geometries to the CartesianTransfo
      ! => initial rotation matrix
      IF (mole%Cart_transfo) THEN
        write(out_unitp,*) '===================================='
        write(out_unitp,*) '==== CartesianTransfo =============='
        CALL flush_perso(out_unitp)

        CALL alloc_dnSVM(dnx,mole%ncart,mole%nb_act,nderiv=0)


        ! The reference geometries are not read in "Read_CartesianTransfo"
        IF (.NOT. mole%tab_Cart_transfo(1)%CartesianTransfo%ReadRefGeometry) THEN

          IF (mole%tab_Cart_transfo(1)%CartesianTransfo%nb_RefGeometry > 1) THEN
            ! Here several reference geometries => They MUST be calculated with sub_QactTOdnx

            DO iref=1,mole%tab_Cart_transfo(1)%CartesianTransfo%nb_RefGeometry

              CALL sub_QactTOdnx(Qact,dnx,mole,0,Gcenter=.FALSE.,Cart_Transfo=.FALSE.)

              mole%tab_Cart_transfo(1)%CartesianTransfo%Qxyz(:,:,iref) =     &
                 reshape(dnx%d0(1:mole%ncart_act),(/ 3,mole%ncart_act/3 /) )

            END DO

          ELSE
            ! Here ONE reference geometry
            iref = 1
            IF (read_xyz0) THEN

              mole%tab_Cart_transfo(1)%CartesianTransfo%Qxyz(:,:,iref) =       &
                  reshape(Qread(1:mole%ncart_act),(/ 3,mole%ncart_act/3 /) )

            ELSE

              IF (.NOT. associated(mole%tab_Cart_transfo(1)%CartesianTransfo%Qxyz)) THEN
                STOP 'CartesianTransfo%Qxyz NOT allocated'
              END IF
              CALL sub_QactTOdnx(Qact,dnx,mole,0,Gcenter=.FALSE.,Cart_Transfo=.FALSE.)
              !CALL Write_XYZ(dnx%d0,mole)

              mole%tab_Cart_transfo(1)%CartesianTransfo%Qxyz(:,:,iref) =     &
                 reshape(dnx%d0(1:mole%ncart_act),(/ 3,mole%nat_act /) )

            END IF

          END IF

        END IF
        CALL dealloc_dnSVM(dnx)

        IF (mole%tab_Cart_transfo(1)%CartesianTransfo%P_Axis_Ref) THEN
          ! this is valid also for Eckart
          CALL Set_P_Axis_CartesianTransfo(mole%tab_Cart_transfo(1)%CartesianTransfo)
        ELSE
          IF (mole%tab_Cart_transfo(1)%CartesianTransfo%Eckart .OR.        &
             mole%tab_Cart_transfo(1)%CartesianTransfo%MultiRefEckart) THEN
            CALL Set_Eckart_CartesianTransfo(mole%tab_Cart_transfo(1)%CartesianTransfo)
          END IF
        END IF

        IF (mole%tab_Cart_transfo(1)%CartesianTransfo%New_Orient) THEN
          CALL Set_NewOrient_CartesianTransfo(mole)
        END IF

        IF (debug) CALL Write_CartesianTransfo(mole%tab_Cart_transfo(1)%CartesianTransfo)

      END IF
!======================================================================
!     END Cart_transfo=t
!======================================================================

      IF(MPI_id==0) THEN
        IF (print_level > 1)                                              &
            write(out_unitp,*) '===================================='
        write(out_unitp,*) 'END ',name_sub
      ENDIF

      CALL dealloc_NParray(QRead,"QRead",name_sub)

      END SUBROUTINE read_RefGeom_CoordType
!=======================================================================================

      SUBROUTINE sub_QactTOQit(Qact,Qit,it_QTransfo,mole,print_Qtransfo)
      USE mod_system
      USE mod_dnSVM
      USE mod_Tnum
      IMPLICIT NONE


      TYPE (CoordType), intent(in)                  :: mole
      integer, intent(in)                           :: it_QTransfo
      real (kind=Rkind), intent(in)                 :: Qact(:)
      real (kind=Rkind), intent(inout), allocatable :: Qit(:)

      logical, intent(in), optional                 :: print_Qtransfo


!     - working variables -------------------------
      integer           :: it,nderiv
      TYPE (Type_dnVec) :: dnQin,dnQout,dnx
      logical           :: print_Qtransfo_loc,Cart_transfo

!     -----------------------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='sub_QactTOQit'
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'it_QTransfo',it_QTransfo
        !write(out_unitp,*) 'Qact =',Qact
        write(out_unitp,*)
        !CALL Write_mole(mole)
        write(out_unitp,*)
        CALL flush_perso(out_unitp)
      END IF
!     -----------------------------------------------------------------
!      IF (size(Qact) == 0) THEN
!        write(out_unitp,*) 'ERROR in ',name_sub
!        write(out_unitp,*) ' the size of Qact(:) is zero!'
!        write(out_unitp,*) ' Check the Frantran source!!'
!        STOP
!      END IF

      IF (allocated(Qit)) CALL dealloc_NParray(Qit,'Qit',name_sub)

      IF (it_QTransfo < 0 .OR. it_QTransfo > mole%nb_Qtransfo) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) ' it_QTransfo has a wrong value: ',it_QTransfo
        write(out_unitp,'(a,i0,a)') ' it must be between [0:',mole%nb_Qtransfo,']'
        write(out_unitp,*) ' Check the Frantran source!!'
        STOP
      END IF

      nderiv = 0
      print_Qtransfo_loc = .FALSE.
      IF (present(print_Qtransfo)) print_Qtransfo_loc = print_Qtransfo
      print_Qtransfo_loc = print_Qtransfo_loc .OR. debug


      IF (it_QTransfo == 0) THEN ! Qcart ! (special case for Cart_transfo)

        CALL alloc_dnSVM(dnx,mole%ncart,mole%nb_act,nderiv=0)

        Cart_transfo = mole%Cart_transfo
        IF (mole%Rot_Dip_with_EC) Cart_transfo = .FALSE.

        CALL sub_QactTOdnx(Qact,dnx,mole,nderiv=0,Gcenter=.FALSE.,      &
                                              Cart_Transfo=Cart_transfo)

        CALL alloc_NParray(Qit,(/mole%ncart/),'Qit',name_sub)
        Qit(:) = dnx%d0(1:mole%ncart)

        IF (print_Qtransfo_loc .OR. debug) THEN
          CALL Write_d0Q(it,'Qin (Qact)',Qact,6)
          CALL Write_d0Q(0, 'Qit (Cart)',Qit,3)
          CALL flush_perso(out_unitp)
        END IF

        CALL dealloc_dnSVM(dnx)


      ELSE IF (it_QTransfo == mole%nb_Qtransfo) THEN    ! Qact only
        ! it enables to add the constraints (rigid, flexible)

        it=mole%nb_Qtransfo
        CALL alloc_dnSVM(dnQin,mole%tab_Qtransfo(it)%nb_Qout,mole%nb_act,nderiv)
        dnQin%d0(:) = Qact(:)


        CALL alloc_dnSVM(dnQout,mole%tab_Qtransfo(it)%nb_Qout,mole%nb_act,nderiv)

        CALL calc_Qtransfo(dnQin,dnQout,mole%tab_Qtransfo(it),nderiv,.TRUE.)

        IF (print_Qtransfo_loc .OR. debug) THEN

          write(out_unitp,*) '-----------------------------------------'
          write(out_unitp,*) 'name_transfo',it,mole%tab_Qtransfo(it)%name_transfo

          CALL Write_d0Q(it,'Qin  (Qact)',dnQin%d0 ,6)
          CALL Write_d0Q(it,'Qout (Qdyn)',dnQout%d0,6)
          CALL flush_perso(out_unitp)

        END IF

        ! Here we have to use dnQin, because we want to set-up Qact (with the constraints)
        CALL alloc_NParray(Qit,(/dnQin%nb_var_vec/),'Qit',name_sub)
        Qit(:) = dnQin%d0(:)

        CALL dealloc_dnSVM(dnQin)
        CALL dealloc_dnSVM(dnQout)


      ELSE

        CALL alloc_dnSVM(dnQin,mole%tab_Qtransfo(mole%nb_Qtransfo)%nb_Qout,mole%nb_act,nderiv)
        dnQin%d0(:) = Qact(:)

        DO it=mole%nb_Qtransfo,it_QTransfo+1,-1 ! we add 1 to it_QTransfo,
                ! because it_QTransfo is set for Qin and at a given iteration we get Qout

          CALL alloc_dnSVM(dnQout,mole%tab_Qtransfo(it)%nb_Qout,mole%nb_act,nderiv)

          CALL calc_Qtransfo(dnQin,dnQout,                              &
                                    mole%tab_Qtransfo(it),nderiv,.TRUE.)

          IF (print_Qtransfo_loc .OR. debug) THEN
            write(out_unitp,*) '-----------------------------------------'
            write(out_unitp,*) 'name_transfo',it,mole%tab_Qtransfo(it)%name_transfo
            CALL Write_d0Q(it,'Qin ',dnQin%d0 ,6)
            CALL Write_d0Q(it,'Qout',dnQout%d0,6)
            CALL flush_perso(out_unitp)

          END IF

          CALL dealloc_dnSVM(dnQin)
          CALL alloc_dnSVM(dnQin,mole%tab_Qtransfo(it)%nb_Qout,mole%nb_act,nderiv)
          CALL sub_dnVec1_TO_dnVec2(dnQout,dnQin)
          CALL dealloc_dnSVM(dnQout)

        END DO

        ! Here we have to use dnQin, because dnQout is transferred in dnQin and then it is deallocated.
        CALL alloc_NParray(Qit,(/dnQin%nb_var_vec/),'Qit',name_sub)
        Qit(:) = dnQin%d0(:)


        CALL dealloc_dnSVM(dnQin)


      END IF

!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
        write(out_unitp,*)
        CALL flush_perso(out_unitp)
      END IF
!     -----------------------------------------------------------------

      END SUBROUTINE sub_QactTOQit

      SUBROUTINE sub_QinRead_TO_Qact(Qread,Qact,mole,it_QinRead)
      USE mod_system
      USE mod_dnSVM
      USE mod_Tnum
      IMPLICIT NONE


      integer           :: it_QinRead
      TYPE (CoordType)  :: mole
      real (kind=Rkind) :: Qread(:)
      real (kind=Rkind) :: Qact(:)

      real (kind=Rkind) :: Rot_initial(3,3),Qat1(3)


      !- working variables -------------------------
      integer :: it,nb_act,i,ic1,ic
      integer :: it_QoutRead

      TYPE (Type_dnVec) :: dnQin,dnQout

      !-----------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='sub_QinRead_TO_Qact'
      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'it_QinRead,mole%nb_Qtransfo',it_QinRead,mole%nb_Qtransfo
        write(out_unitp,*) 'Qread =',Qread
        write(out_unitp,*)
        !CALL Write_mole(mole)
        write(out_unitp,*)
      END IF
      !-----------------------------------------------------------------

      ! since it is going from out to in, it is better to use it_QoutRead (= it_QinRead+1)
      it_QoutRead = it_QinRead + 1

      IF (it_QoutRead == mole%nb_Qtransfo+1) THEN ! read_Qact0
        Qact(:) = Qread(:)
      ELSE
        it = it_QoutRead
        nb_act = mole%tab_Qtransfo(it_QoutRead)%nb_act

        CALL alloc_dnSVM(dnQout,mole%tab_Qtransfo(it)%nb_Qout,nb_act,0)

        dnQout%d0(1:size(Qread)) = Qread(:)

        DO it=it_QoutRead,mole%nb_Qtransfo

          CALL alloc_dnSVM(dnQin,mole%tab_Qtransfo(it)%nb_Qin,nb_act,0)

          IF (debug) THEN
            CALL Write_d0Q(it,'Qout ' // trim(adjustl(mole%tab_Qtransfo(it)%name_transfo)),dnQout%d0,6)
            write(out_unitp,*) 'Qout ',it,mole%tab_Qtransfo(it)%name_transfo,dnQout%d0
            CALL flush_perso(out_unitp)
          END IF

          CALL calc_Qtransfo(dnQin,dnQout,mole%tab_Qtransfo(it),0,inTOout=.FALSE.)

          IF (debug) THEN
            CALL Write_d0Q(it,'Qin  ' // trim(adjustl(mole%tab_Qtransfo(it)%name_transfo)),dnQin%d0,6)
            CALL flush_perso(out_unitp)
          END IF

          CALL dealloc_dnSVM(dnQout)
          CALL alloc_dnSVM(dnQout,mole%tab_Qtransfo(it)%nb_Qin,nb_act,0)

          CALL sub_dnVec1_TO_dnVec2(dnQin,dnQout,nderiv=0)
          CALL dealloc_dnSVM(dnQin)

        END DO

        Qact(:) = dnQout%d0(1:size(Qact))
        CALL dealloc_dnSVM(dnQout)
      END IF



      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Qact',Qact(:)
        write(out_unitp,*) 'END ',name_sub
        write(out_unitp,*)
      END IF
      CALL flush_perso(out_unitp)
      !-----------------------------------------------------------------

      END SUBROUTINE sub_QinRead_TO_Qact

      SUBROUTINE Set_NewOrient_CartesianTransfo(mole)
      USE mod_system
      USE mod_dnSVM
      USE mod_Tnum
      IMPLICIT NONE

      TYPE (CoordType),   intent(inout)      :: mole



!     - working variables -------------------------
      integer           :: iref,iat
      real (kind=Rkind) :: Qxyz(mole%ncart)
      logical           :: vAti_EQ_ZERO


!     -----------------------------------------------------------------
      !logical, parameter :: debug = .FALSE.
      logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='Set_NewOrient_CartesianTransfo'
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
      END IF
!     -----------------------------------------------------------------

      vAti_EQ_ZERO = all(mole%tab_Cart_transfo(1)%CartesianTransfo%vAt1 == ZERO) .AND. &
                     all(mole%tab_Cart_transfo(1)%CartesianTransfo%vAt2 == ZERO) .AND. &
                     all(mole%tab_Cart_transfo(1)%CartesianTransfo%vAt3 == ZERO)

      IF (.NOT. vAti_EQ_ZERO) THEN

        CALL sub_vAtiTOexeyez(mole%tab_Cart_transfo(1)%CartesianTransfo%TransVect(:,1),mole)

        mole%tab_Cart_transfo(1)%CartesianTransfo%Qxyz(:,:,:)   = ZERO
        mole%tab_Cart_transfo(1)%CartesianTransfo%MWQxyz(:,:,:) = ZERO

      ELSE

        IF (mole%ncart /= mole%ncart_act+3) THEN
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) '  ncart /= ncart_act+3',mole%ncart,mole%ncart_act
          write(out_unitp,*) ' => dummy atoms are present. It is not possible with this option.'
          write(out_unitp,*) ' => Choose New_Orient with vAt1,vAt3,vAt3. '
          write(out_unitp,*) '  CHECK your data!!'
          STOP
        END IF

        DO iref=1,mole%tab_Cart_transfo(1)%CartesianTransfo%nb_RefGeometry

          Qxyz(1:mole%ncart_act) = reshape(mole%tab_Cart_transfo(1)%CartesianTransfo%Qxyz(:,:,iref), shape=(/ mole%ncart_act /))

          CALL sub_QxyzTOexeyez(Qxyz,mole%tab_Cart_transfo(1)%CartesianTransfo%TransVect(:,iref),mole)

          mole%tab_Cart_transfo(1)%CartesianTransfo%Qxyz(:,:,iref) = reshape(Qxyz(1:mole%ncart_act), shape=(/ 3,mole%nat_act /))

          DO iat=1,mole%tab_Cart_transfo(1)%CartesianTransfo%nat_act
            mole%tab_Cart_transfo(1)%CartesianTransfo%MWQxyz(:,iat,iref) = &
             mole%tab_Cart_transfo(1)%CartesianTransfo%Qxyz(:,iat,iref) *  &
               sqrt(mole%tab_Cart_transfo(1)%CartesianTransfo%masses_at(iat))
          END DO

        END DO
      END IF

      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
        write(out_unitp,*)
      END IF
      CALL flush_perso(out_unitp)

      END SUBROUTINE Set_NewOrient_CartesianTransfo

      SUBROUTINE sub_QxyzTOexeyez(Qxyz,VT,mole)
      USE mod_system
      USE mod_dnSVM
      USE mod_Tnum
      IMPLICIT NONE


      TYPE (CoordType)  :: mole
      real (kind=Rkind) :: Qxyz(:),VT(:)



!     - working variables -------------------------
      logical           :: case1
      integer :: i,it,nb_act,ncart,nc1,nc2,nc3
      TYPE (Type_dnVec) :: dnQin,dnQout
      real (kind=Rkind) :: ex(3),nx,ey(3),ny,ez(3),nz

!     -----------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='sub_QxyzTOexeyez'
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Qxyz (not mass-weigthed, Angs)'
        DO i=1,mole%nat_act
          write(out_unitp,*) Qxyz(3*i-2:3*i) * get_Conv_au_TO_unit("L","Angs")
        END DO
        write(out_unitp,*) 'Qxyz (not mass-weigthed, Bohr)'
        DO i=1,mole%nat_act
          write(out_unitp,*) Qxyz(3*i-2:3*i)
        END DO
        write(out_unitp,*) 'num_transfo',mole%tab_Qtransfo(1)%num_transfo
        write(out_unitp,*) 'name_transfo ',mole%tab_Qtransfo(1)%name_transfo
      END IF
!     -----------------------------------------------------------------

      IF (.NOT. mole%Cart_transfo) RETURN
      IF (.NOT. associated(mole%tab_Cart_transfo)) RETURN


      VT(:) = ZERO

      SELECT CASE (mole%tab_Qtransfo(1)%name_transfo)
      CASE ('zmat')
        IF (debug) write(out_unitp,*) 'zmat'
        nc1 = mole%tab_Qtransfo(1)%ZmatTransfo%ind_zmat(1,1)
        nc2 = mole%tab_Qtransfo(1)%ZmatTransfo%ind_zmat(1,2)
        nc3 = mole%tab_Qtransfo(1)%ZmatTransfo%ind_zmat(1,3)

        VT(:) = Qxyz(nc1:nc1+2)
        DO i=1,mole%nat_act
          Qxyz(3*i-2:3*i) = Qxyz(3*i-2:3*i)-VT(:)
        END DO

        IF (debug) write(out_unitp,*) 'zmat, nc1,nc2,nc3',nc1,nc2,nc3

        IF (nc1 <= mole%ncart .AND. nc2 <= mole%ncart .AND. nc3 <= mole%ncart) THEN

          ez(:) = Qxyz(nc2:nc2+2)
          ez(:) = ez(:)/sqrt(dot_product(ez,ez))

          case1 = (mole%tab_Qtransfo(1)%ZmatTransfo%ind_zmat(2,3) ==    &
                 mole%tab_Qtransfo(1)%ZmatTransfo%ind_zmat(1,1) )

          IF (case1) THEN
            ex(:) = Qxyz(nc3:nc3+2)
          ELSE
            ex(:) = Qxyz(nc3:nc3+2)-Qxyz(nc2:nc2+2)
          END IF
          ex(:) = ex(:) - ez(:)*dot_product(ez,ex)
          ex(:) = ex(:)/sqrt(dot_product(ex,ex))

          CALL calc_cross_product(ez,nz,ex,nx,ey,ny)
        ELSE
          ex(:) = (/ONE,ZERO,ZERO/)
          ey(:) = (/ZERO,ONE,ZERO/)
          ez(:) = (/ZERO,ZERO,ONE/)
        END IF

      CASE ('bunch','bunch_poly')

        it = 1
        nb_act = mole%tab_Qtransfo(it)%nb_act
        CALL alloc_dnSVM(dnQout,mole%tab_Qtransfo(it)%nb_Qout,nb_act,0)

        dnQout%d0(1:size(Qxyz)) = Qxyz(:)

       CALL alloc_dnSVM(dnQin,mole%tab_Qtransfo(it)%nb_Qin,nb_act,0)

        CALL Write_d0Q(it,'Qxyz ' // trim(adjustl(mole%tab_Qtransfo(it)%name_transfo)),dnQout%d0,3)
        CALL calc_Qtransfo(dnQin,dnQout,mole%tab_Qtransfo(it),0,inTOout=.FALSE.)
        DO i=1,3*mole%tab_Qtransfo(it)%BunchTransfo%nb_vect,3
          write(out_unitp,*) 'QVect',int(i/3)+1,                        &
                    sqrt(dot_product(dnQin%d0(i:i+2),dnQin%d0(i:i+2))), &
                                                         dnQin%d0(i:i+2)
        END DO

        ez(:) = dnQin%d0(1:3)
        ez(:) = ez(:)/sqrt(dot_product(ez,ez))

        ex(:) = dnQin%d0(4:6)
        ex(:) = ex(:) - ez(:)*dot_product(ez,ex)
        ex(:) = ex(:)/sqrt(dot_product(ex,ex))

        CALL calc_cross_product(ez,nz,ex,nx,ey,ny)

        CALL dealloc_dnSVM(dnQout)
        CALL dealloc_dnSVM(dnQin)

      CASE ('QTOX_ana')
        write(out_unitp,*) 'QTOX_ana: correct orientation ???'
        ez(:) = Qxyz(1:3)
        ez(:) = ez(:)/sqrt(dot_product(ez,ez))

        ex(:) = Qxyz(4:6)
        ex(:) = ex(:) - ez(:)*dot_product(ez,ex)
        ex(:) = ex(:)/sqrt(dot_product(ex,ex))

        CALL calc_cross_product(ez,nz,ex,nx,ey,ny)
      CASE default ! ERROR: wrong transformation !
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) '  Wrong transformation !!'
        write(out_unitp,*) 'name_transfo',mole%tab_Qtransfo(1)%name_transfo
        write(out_unitp,*) '  CHECK the fortran!!'
        STOP
      END SELECT

      mole%tab_Cart_transfo(1)%CartesianTransfo%Rot_initial(:,1) = ex(:)
      mole%tab_Cart_transfo(1)%CartesianTransfo%Rot_initial(:,2) = ey(:)
      mole%tab_Cart_transfo(1)%CartesianTransfo%Rot_initial(:,3) = ez(:)

!=================================================


!=================================================
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'ex',ex(:)
        write(out_unitp,*) 'ey',ey(:)
        write(out_unitp,*) 'ez',ez(:)
        write(out_unitp,*) 'END ',name_sub
        write(out_unitp,*)
      END IF
      CALL flush_perso(out_unitp)

!     -----------------------------------------------------------------
!=================================================

      END SUBROUTINE sub_QxyzTOexeyez
      !Initial rotation with vAt1, vAt2, vAt3
      SUBROUTINE sub_vAtiTOexeyez(VT,mole)
      USE mod_system
      USE mod_dnSVM
      USE mod_Tnum
      IMPLICIT NONE


      TYPE (CoordType)  :: mole
      real (kind=Rkind) :: VT(:)



!     - working variables -------------------------
      logical           :: case1
      integer :: i,it,nb_act,ncart,nc1,nc2,nc3
      TYPE (Type_dnVec) :: dnQin,dnQout
      real (kind=Rkind) :: ex(3),nx,ey(3),ny,ez(3),nz
      logical :: vAti_EQ_ZERO

!     -----------------------------------------------------------------
      !logical, parameter :: debug = .FALSE.
      logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='sub_vAtiTOexeyez'
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'vAt1',mole%tab_Cart_transfo(1)%CartesianTransfo%vAt1
        write(out_unitp,*) 'vAt2',mole%tab_Cart_transfo(1)%CartesianTransfo%vAt2
        write(out_unitp,*) 'vAt3',mole%tab_Cart_transfo(1)%CartesianTransfo%vAt3
      END IF
!     -----------------------------------------------------------------


      vAti_EQ_ZERO = all(mole%tab_Cart_transfo(1)%CartesianTransfo%vAt1 == ZERO) .AND. &
                     all(mole%tab_Cart_transfo(1)%CartesianTransfo%vAt2 == ZERO) .AND. &
                     all(mole%tab_Cart_transfo(1)%CartesianTransfo%vAt3 == ZERO)

     IF (vAti_EQ_ZERO) THEN
       write(out_unitp,*) 'ERROR in ',name_sub
       write(out_unitp,*) '  vAt1=0 and vAt2=0 and vAt3=0 '
       write(out_unitp,*) '  CHECK you data!!'
       STOP
     END IF

      VT(:) = mole%tab_Cart_transfo(1)%CartesianTransfo%vAt1(:)

      SELECT CASE (mole%tab_Qtransfo(1)%name_transfo)
      CASE ('zmat')
        IF (debug) write(out_unitp,*) 'zmat'
        nc1 = mole%tab_Qtransfo(1)%ZmatTransfo%ind_zmat(1,1)
        nc2 = mole%tab_Qtransfo(1)%ZmatTransfo%ind_zmat(1,2)
        nc3 = mole%tab_Qtransfo(1)%ZmatTransfo%ind_zmat(1,3)


        IF (debug) write(out_unitp,*) 'zmat, nc1,nc2,nc3',nc1,nc2,nc3

        IF (nc1 <= mole%ncart .AND. nc2 <= mole%ncart .AND. nc3 <= mole%ncart) THEN

          ez(:) = mole%tab_Cart_transfo(1)%CartesianTransfo%vAt2(:)-VT(:)
          ez(:) = ez(:)/sqrt(dot_product(ez,ez))

          case1 = (mole%tab_Qtransfo(1)%ZmatTransfo%ind_zmat(2,3) ==    &
                 mole%tab_Qtransfo(1)%ZmatTransfo%ind_zmat(1,1) )

          IF (case1) THEN
            ex(:) = mole%tab_Cart_transfo(1)%CartesianTransfo%vAt3(:)-VT(:)
          ELSE
            ex(:) = mole%tab_Cart_transfo(1)%CartesianTransfo%vAt3(:)-mole%tab_Cart_transfo(1)%CartesianTransfo%vAt2(:)
          END IF
          ex(:) = ex(:) - ez(:)*dot_product(ez,ex)
          ex(:) = ex(:)/sqrt(dot_product(ex,ex))

          CALL calc_cross_product(ez,nz,ex,nx,ey,ny)
        ELSE
          ex(:) = (/ONE,ZERO,ZERO/)
          ey(:) = (/ZERO,ONE,ZERO/)
          ez(:) = (/ZERO,ZERO,ONE/)
        END IF

      CASE ('bunch','bunch_poly')

        STOP ' not yet with bunch or bunch_poly'

      CASE ('QTOX_ana')

        STOP ' not yet with QTOX_ana'


      CASE default ! ERROR: wrong transformation !
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) '  Wrong transformation !!'
        write(out_unitp,*) 'name_transfo',mole%tab_Qtransfo(1)%name_transfo
        write(out_unitp,*) '  CHECK the fortran!!'
        STOP
      END SELECT

      mole%tab_Cart_transfo(1)%CartesianTransfo%Rot_initial(:,1) = ex(:)
      mole%tab_Cart_transfo(1)%CartesianTransfo%Rot_initial(:,2) = ey(:)
      mole%tab_Cart_transfo(1)%CartesianTransfo%Rot_initial(:,3) = ez(:)

!=================================================


!=================================================
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'ex',ex(:)
        write(out_unitp,*) 'ey',ey(:)
        write(out_unitp,*) 'ez',ez(:)
        write(out_unitp,*) 'END ',name_sub
        write(out_unitp,*)
      END IF
      CALL flush_perso(out_unitp)

!     -----------------------------------------------------------------
!=================================================

      END SUBROUTINE sub_vAtiTOexeyez
      SUBROUTINE sub_Qxyz0TORot(Qxyz,Rot_initial,mole)
      USE mod_system
      USE mod_dnSVM
      USE mod_Tnum
      IMPLICIT NONE


      TYPE (CoordType)  :: mole
      real (kind=Rkind) :: Qxyz(:)
      real (kind=Rkind) :: Rot_initial(3,3)

!     - working variables -------------------------
      logical           :: case1
      integer :: i,it,nb_act,ncart,nc1,nc2,nc3
      TYPE (Type_dnVec) :: dnQin,dnQout
      real (kind=Rkind) :: ex(3),nx,ey(3),ny,ez(3),nz

!     -----------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='sub_Qxyz0TORot'
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Qxyz =',Qxyz
        write(out_unitp,*) 'num_transfo',mole%tab_Qtransfo(1)%num_transfo
        write(out_unitp,*) 'name_transfo ',mole%tab_Qtransfo(1)%name_transfo
      END IF
!     -----------------------------------------------------------------

      ncart = min(size(Qxyz),size(mole%d0sm))

      SELECT CASE (mole%tab_Qtransfo(1)%name_transfo)
      CASE ('zmat')
        nc1 = mole%tab_Qtransfo(1)%ZmatTransfo%ind_zmat(1,1)
        nc2 = mole%tab_Qtransfo(1)%ZmatTransfo%ind_zmat(1,2)
        nc3 = mole%tab_Qtransfo(1)%ZmatTransfo%ind_zmat(1,3)

        ez(:) = Qxyz(nc2:nc2+2)-Qxyz(nc1:nc1+2)
        ez(:) = ez(:)/sqrt(dot_product(ez,ez))

        case1 = (mole%tab_Qtransfo(1)%ZmatTransfo%ind_zmat(2,3) ==      &
                 mole%tab_Qtransfo(1)%ZmatTransfo%ind_zmat(1,1) )

        IF (case1) THEN
          ex(:) = Qxyz(nc3:nc3+2)-Qxyz(nc1:nc1+2)
        ELSE
          ex(:) = Qxyz(nc3:nc3+2)-Qxyz(nc2:nc2+2)
        END IF
        ex(:) = ex(:) - ez(:)*dot_product(ez,ex)
        ex(:) = ex(:)/sqrt(dot_product(ex,ex))

        CALL calc_cross_product(ez,nz,ex,nx,ey,ny)

      CASE ('bunch','bunch_poly')

        it = 1
        nb_act = mole%tab_Qtransfo(it)%nb_act
        CALL alloc_dnSVM(dnQout,mole%tab_Qtransfo(it)%nb_Qout,nb_act,0)

        dnQout%d0(1:size(Qxyz)) = Qxyz(:)

       CALL alloc_dnSVM(dnQin,mole%tab_Qtransfo(it)%nb_Qin,nb_act,0)

        IF (debug) THEN
          CALL Write_d0Q(it,'Qxyz ' // trim(adjustl(mole%tab_Qtransfo(it)%name_transfo)),dnQout%d0,3)
        END IF
        CALL calc_Qtransfo(dnQin,dnQout,mole%tab_Qtransfo(it),0,inTOout=.FALSE.)

        IF (debug) THEN
          DO i=1,3*mole%tab_Qtransfo(it)%BunchTransfo%nb_vect,3
            write(out_unitp,*) 'QVect',int(i/3)+1,                      &
                    sqrt(dot_product(dnQin%d0(i:i+2),dnQin%d0(i:i+2))), &
                                                         dnQin%d0(i:i+2)
          END DO
        END IF

        ez(:) = dnQin%d0(1:3)
        ez(:) = ez(:)/sqrt(dot_product(ez,ez))

        ex(:) = dnQin%d0(4:6)
        ex(:) = ex(:) - ez(:)*dot_product(ez,ex)
        ex(:) = ex(:)/sqrt(dot_product(ex,ex))

        CALL calc_cross_product(ez,nz,ex,nx,ey,ny)

        CALL dealloc_dnSVM(dnQout)
        CALL dealloc_dnSVM(dnQin)

      CASE ('QTOX_ana')
        write(out_unitp,*) 'QTOX_ana: correct orientation ???'
        ez(:) = Qxyz(1:3)
        ez(:) = ez(:)/sqrt(dot_product(ez,ez))

        ex(:) = Qxyz(4:6)
        ex(:) = ex(:) - ez(:)*dot_product(ez,ex)
        ex(:) = ex(:)/sqrt(dot_product(ex,ex))

        CALL calc_cross_product(ez,nz,ex,nx,ey,ny)

      CASE default ! ERROR: wrong transformation !
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) '  Wrong transformation !!'
        write(out_unitp,*) 'name_transfo',mole%tab_Qtransfo(1)%name_transfo
        write(out_unitp,*) '  CHECK the fortran!!'
        STOP
      END SELECT

      Rot_initial(:,1) = ex(:)
      Rot_initial(:,2) = ey(:)
      Rot_initial(:,3) = ez(:)

!=================================================


!=================================================
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Rotational matrix'
        CALL Write_Mat(Rot_initial,out_unitp,5)
        write(out_unitp,*) 'END ',name_sub
        write(out_unitp,*)
      END IF
      CALL flush_perso(out_unitp)

!     -----------------------------------------------------------------
!=================================================

      END SUBROUTINE sub_Qxyz0TORot


      SUBROUTINE sub_QplusDQ_TO_Cart(Qact,mole)
      USE mod_system
      USE mod_dnSVM
      USE mod_Tnum
      IMPLICIT NONE


      TYPE (CoordType) :: mole
      real (kind=Rkind), intent(inout) :: Qact(:)
      TYPE (Type_dnVec) :: dnx0
      TYPE (Type_dnVec) :: dnx

      real (kind=Rkind) :: a0,Norm
      integer           :: Z_act(mole%nat)

      integer           :: i,iZ,iQ,niofreq
      TYPE (param_file) :: file_freq


!     -----------------------------------------------------------------
      integer :: nderiv_debug = 0
      logical, parameter :: debug = .FALSE.
!      logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='sub_QplusDQ_TO_Cart'
!     -----------------------------------------------------------------

!     -----------------------------------------------------------------
      ! Some initializations
!     -----------------------------------------------------------------
      Z_act(:) = -1
      iZ = 0
      DO i=1,mole%nat
        IF (mole%Z(i) > 0) THEN
          iZ = iZ + 1
          Z_act(iZ) = mole%Z(i)
        END IF
      END DO
      a0 = get_Conv_au_TO_unit("L","Angs")
!     -----------------------------------------------------------------
!     -----------------------------------------------------------------


      CALL alloc_dnSVM(dnx0,mole%ncart,mole%nb_act,0)
      CALL alloc_dnSVM(dnx,mole%ncart,mole%nb_act,0)

      Qact = mole%ActiveTransfo%Qact0
      CALL sub_QactTOdnx(Qact,dnx0,mole,0,Gcenter=.TRUE.)

      write(out_unitp,*) '=============================================='
      write(out_unitp,*) '= XYZ format (reference geometry) ============'
      write(out_unitp,*) mole%nat_act
      write(out_unitp,*)

      iZ = 0
      DO i=1,mole%ncart_act,3
        iZ = iZ + 1
        write(out_unitp,112) Z_act(iZ),dnx0%d0(i:i+2)*a0
 112    format(2x,i5,3(2x,f12.5))

      END DO

      write(out_unitp,*) '= END XYZ format ============================='
      write(out_unitp,*) '=============================================='

      file_freq%name='freq.xyz'
      CALL file_open(file_freq,niofreq)
      ! loop on all the coordinates (active order)
      DO iQ=1,mole%nb_var

        Qact = mole%ActiveTransfo%Qact0
        Qact(iQ) = Qact(iQ) + ONETENTH
        CALL sub_QactTOdnx(Qact,dnx,mole,0,Gcenter=.TRUE.)
        dnx%d0 = dnx%d0 - dnx0%d0 ! dxyz
        Norm = sqrt(dot_product(dnx%d0,dnx%d0))
        dnx%d0 = dnx%d0/Norm
        write(niofreq,*) mole%nat_act
        write(niofreq,*) '  Coord: ',iQ

        iZ = 0
        DO i=1,mole%ncart_act,3
          iZ = iZ + 1
          write(niofreq,113) Z_act(iZ),dnx0%d0(i:i+2)*a0,0,dnx%d0(i:i+2)*a0
 113      format(2x,i5,3(2x,f12.5),i5,3(2x,f12.5))

        END DO
      END DO
      CALL file_close(file_freq)

      Qact = mole%ActiveTransfo%Qact0

      END SUBROUTINE sub_QplusDQ_TO_Cart

!================================================================
!       conversion d0Q (zmat,poly, bunch ...) => d0x (mass weighted)
!================================================================
      SUBROUTINE sub_QactTOdnMWx(Qact,dnMWx,mole,nderiv,Gcenter,Cart_Transfo,WriteCC)
      USE mod_system
      USE mod_dnSVM
      USE mod_Tnum
      IMPLICIT NONE


      real (kind=Rkind), intent(in)             :: Qact(:)
      TYPE (CoordType),  intent(in)             :: mole
      TYPE (Type_dnVec), intent(inout)          :: dnMWx
      integer,           intent(in)             :: nderiv
      logical,           intent(in)             :: Gcenter
      logical,           intent(in),   optional :: Cart_Transfo,WriteCC


!     - working variables -------------------------
      logical           :: Cart_Transfo_loc,WriteCC_loc


!     -----------------------------------------------------------------
      integer :: nderiv_debug = 1
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='sub_QactTOdnMWx'
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nderiv',nderiv
        write(out_unitp,*) 'ncart',mole%ncart
        write(out_unitp,*) 'Qact =',Qact
        write(out_unitp,*)
        !CALL Write_CoordType(mole)
        write(out_unitp,*)
        CALL write_dnx(1,mole%ncart,dnMWx,nderiv_debug)
      END IF
!     -----------------------------------------------------------------

      IF (present(WriteCC)) THEN
        WriteCC_loc = WriteCC
      ELSE
        WriteCC_loc = mole%WriteCC
      END IF

      IF (present(Cart_Transfo)) THEN
        Cart_Transfo_loc = Cart_Transfo
      ELSE
        Cart_Transfo_loc = mole%Cart_transfo
      END IF

      CALL sub_QactTOdnx(Qact,dnMWx,mole,nderiv,Gcenter,Cart_Transfo_loc,WriteCC_loc)

      CALL sub_dnxMassWeight(dnMWx,mole%d0sm,mole%ncart,mole%ncart_act,nderiv)


      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Mass Weighted Cartessian coordinates'
        CALL write_dnx(1,mole%ncart,dnMWx,nderiv_debug)
        write(out_unitp,*) 'END ',name_sub
        write(out_unitp,*)
      END IF
      !-----------------------------------------------------------------

      END SUBROUTINE sub_QactTOdnMWx
!================================================================
!       conversion d0Q (zmat,poly, bunch ...) => d0x (not mass weighted)
!================================================================
      SUBROUTINE sub_QactTOdnx_zmatrix(Qact,dnx,mole,                   &
                                     nderiv,Gcenter,Cart_Transfo,WriteCC)
      USE mod_system
      USE mod_dnSVM
      USE mod_Tnum
      IMPLICIT NONE

      real (kind=Rkind), intent(in)           :: Qact(:)
      TYPE (zmatrix),    intent(in)           :: mole
      TYPE (Type_dnVec), intent(inout)        :: dnx
      integer,           intent(in)           :: nderiv
      logical,           intent(in)           :: Gcenter
      logical,           intent(in), optional :: Cart_Transfo,WriteCC

      IF (present(Cart_Transfo)) THEN
        IF (present(WriteCC)) THEN
          CALL sub_QactTOdnx_CoordType(Qact,dnx,mole%CoordType,nderiv,Gcenter,&
                               Cart_Transfo=Cart_Transfo,WriteCC=WriteCC)
        ELSE
          CALL sub_QactTOdnx_CoordType(Qact,dnx,mole%CoordType,nderiv,Gcenter,&
                                              Cart_Transfo=Cart_Transfo)
        END IF
      ELSE
        IF (present(WriteCC)) THEN
          CALL sub_QactTOdnx_CoordType(Qact,dnx,mole%CoordType,nderiv,Gcenter,&
                                                        WriteCC=WriteCC)
        ELSE
          CALL sub_QactTOdnx_CoordType(Qact,dnx,mole%CoordType,nderiv,Gcenter)
        END IF
      END IF

      END SUBROUTINE sub_QactTOdnx_zmatrix
      RECURSIVE SUBROUTINE sub_QactTOdnx_CoordType(Qact,dnx,mole,       &
                                     nderiv,Gcenter,Cart_Transfo,WriteCC)
      USE mod_system
      USE mod_dnSVM
      USE mod_Tnum
      IMPLICIT NONE

      real (kind=Rkind), intent(in)           :: Qact(:)
      TYPE (CoordType),  intent(in)           :: mole
      TYPE (Type_dnVec), intent(inout)        :: dnx
      integer,           intent(in)           :: nderiv
      logical,           intent(in)           :: Gcenter
      logical,           intent(in), optional :: Cart_Transfo,WriteCC

!     - working variables -------------------------
      TYPE (Type_dnVec) :: dnQin,dnQout
      real (kind=Rkind) :: Qacti,Qactj
      real (kind=Rkind) :: step2,step24,stepp
      integer           :: i,j
      integer           :: it,ic,icG,nb_act,iii
      logical           :: Gcenter_loc,Cart_Transfo_loc,GCenter_done,WriteCC_loc
      real (kind=Rkind) :: Qact_loc(size(Qact))


!     -----------------------------------------------------------------
      integer :: nderiv_debug = 1
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='sub_QactTOdnx'
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nderiv',nderiv
        write(out_unitp,*) 'ncart',mole%ncart
        write(out_unitp,*) 'Qact =',Qact
        write(out_unitp,*)
        !CALL Write_CoordType(mole)
        write(out_unitp,*)
        CALL write_dnx(1,mole%ncart,dnx,nderiv_debug)
      END IF
!     -----------------------------------------------------------------

      IF (size(Qact) /= mole%nb_var .AND. size(Qact) /= mole%nb_act) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) ' the size of Qact(:) is not mole%nb_var or '
        write(out_unitp,*) ' the size of Qact(:) is not mole%nb_act!'
        write(out_unitp,*) ' the size of Qact(:): ',size(Qact)
        write(out_unitp,*) ' mole%nb_var:         ',mole%nb_var
        write(out_unitp,*) ' mole%nb_act:         ',mole%nb_act
        write(out_unitp,*) ' Check the Fortran source!!'
        STOP
      END IF

      IF (present(WriteCC)) THEN
        WriteCC_loc = WriteCC
      ELSE
        WriteCC_loc = mole%WriteCC
      END IF

      IF (present(Cart_Transfo)) THEN
        Cart_Transfo_loc = Cart_Transfo
      ELSE
        Cart_Transfo_loc = mole%Cart_transfo
      END IF

      IF (debug) write(out_unitp,*) 'Cart_Transfo_loc, mole%Cart_transfo',Cart_Transfo_loc,mole%Cart_transfo

      IF (mole%num_x .AND. nderiv > 0) THEN
        step2 = ONE/(mole%stepQ*mole%stepQ)
        step24 = step2/FOUR
        stepp = ONE/(mole%stepQ+mole%stepQ)
        IF (nderiv >= 3) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nderiv > 2 is impossible with numerical derivatives'
          write(out_unitp,*) ' nderiv: ',nderiv
          STOP
        END IF
        IF (nderiv >= 1) THEN ! first and second (diagonal) derivatives

         Qact_loc(:) = Qact(:)

         DO i=1,mole%nb_act

            Qacti = Qact_loc(i)

            Qact_loc(i) = Qacti + mole%stepQ
            CALL sub_QactTOdnx(Qact_loc,dnx,mole,0,Gcenter,Cart_Transfo_loc,WriteCC=.FALSE.)
            dnx%d1(:,i) = dnx%d0(:)
            IF (nderiv == 2) dnx%d2(:,i,i) = dnx%d0(:)

            Qact_loc(i) = Qacti - mole%stepQ
            CALL sub_QactTOdnx(Qact_loc,dnx,mole,0,Gcenter,Cart_Transfo_loc,WriteCC=.FALSE.)
            dnx%d1(:,i) = (dnx%d1(:,i) - dnx%d0(:)) * stepp
            IF (nderiv == 2) dnx%d2(:,i,i) = dnx%d2(:,i,i) + dnx%d0

            Qact_loc(i) = Qacti

          END DO ! end first and second (diagonal) derivatives
        END IF

        IF (nderiv == 2) THEN ! second derivatives (crossing term)
          DO i=1,mole%nb_act
          DO j=i+1,mole%nb_act

            Qacti = Qact_loc(i)
            Qactj = Qact_loc(j)

            Qact_loc(i) = Qacti + mole%stepQ
            Qact_loc(j) = Qactj + mole%stepQ
            CALL sub_QactTOdnx(Qact_loc,dnx,mole,0,Gcenter,Cart_Transfo_loc,WriteCC=.FALSE.)
            dnx%d2(:,i,j) = dnx%d0(:)

            Qact_loc(i) = Qacti - mole%stepQ
            Qact_loc(j) = Qactj - mole%stepQ
            CALL sub_QactTOdnx(Qact_loc,dnx,mole,0,Gcenter,Cart_Transfo_loc,WriteCC=.FALSE.)
            dnx%d2(:,i,j) = dnx%d2(:,i,j) + dnx%d0(:)

            Qact_loc(i) = Qacti - mole%stepQ
            Qact_loc(j) = Qactj + mole%stepQ
            CALL sub_QactTOdnx(Qact_loc,dnx,mole,0,Gcenter,Cart_Transfo_loc,WriteCC=.FALSE.)
            dnx%d2(:,i,j) = dnx%d2(:,i,j) - dnx%d0(:)

            Qact_loc(i) = Qacti - mole%stepQ
            Qact_loc(j) = Qactj - mole%stepQ
            CALL sub_QactTOdnx(Qact_loc,dnx,mole,0,Gcenter,Cart_Transfo_loc,WriteCC=.FALSE.)
            dnx%d2(:,i,j) = dnx%d2(:,i,j) - dnx%d0(:)

            dnx%d2(:,i,j) = dnx%d2(:,i,j) * step24
            dnx%d2(:,j,i) = dnx%d2(:,i,j)

            Qact_loc(i) = Qacti
            Qact_loc(j) = Qactj

          END DO
          END DO
        END IF ! end second derivatives (crossing term)

        ! no derivative values
        CALL sub_QactTOdnx(Qact_loc,dnx,mole,0,Gcenter,Cart_Transfo_loc,WriteCC=WriteCC_loc)

        IF (nderiv == 2) THEN
          DO i=1,mole%nb_act
            dnx%d2(:,i,i) = ( dnx%d2(:,i,i) - TWO*dnx%d0(:) ) * step2
          END DO
        END IF

      ELSE

        it = mole%nb_Qtransfo
        nb_act = mole%tab_Qtransfo(mole%nb_Qtransfo)%nb_act

        IF (WriteCC_loc .OR. debug) THEN
          CALL Write_d0Q(it,'Qact',Qact(1:mole%nb_act),6)
        END IF

        CALL alloc_dnSVM(dnQin,mole%tab_Qtransfo(it)%nb_Qout,nb_act,nderiv)
        dnQin%d0(1:size(Qact)) = Qact(:)
        IF (WriteCC_loc .OR. debug) CALL Write_d0Q(it,'Qin (Qact)',dnQin%d0,6)

        DO it=mole%nb_Qtransfo,1,-1
          IF (mole%tab_Qtransfo(it)%skip_transfo) CYCLE
          IF (WriteCC_loc .OR. debug) write(out_unitp,*) 'name_transfo',it,&
                                  mole%tab_Qtransfo(it)%name_transfo
          CALL flush_perso(out_unitp)

          IF (WriteCC_loc .OR. debug) CALL Write_d0Q(it,'Qin ',dnQin%d0,6)

          CALL calc_Qtransfo(dnQin,dnQout,mole%tab_Qtransfo(it),nderiv,.TRUE.)

          IF (WriteCC_loc .OR. debug) THEN
            IF (it == 1) THEN
              CALL Write_d0Q(it,'Qxyz',dnQout%d0,3)
            ELSE
              CALL Write_d0Q(it,'Qout',dnQout%d0,6)
            END IF
          END IF

          CALL dealloc_dnSVM(dnQin)
          CALL alloc_dnSVM(dnQin,mole%tab_Qtransfo(it)%nb_Qout,nb_act,nderiv)
          CALL sub_dnVec1_TO_dnVec2(dnQout,dnQin)
          CALL dealloc_dnSVM(dnQout)

        END DO

        it = 0
        CALL sub_dnVec1_TO_dnVec2(dnQin,dnx)
        CALL dealloc_dnSVM(dnQin)

        !=================================================
        IF((WriteCC_loc .OR. debug) .AND. MPI_id==0) THEN
          write(out_unitp,*) ' Cartesian coordinates without the Cartesian Transformation (au):'
          CALL write_dnx(1,mole%ncart,dnx,nderiv_debug)
          write(out_unitp,*) ' Cartesian coordinates without the Cartesian Transformation (ang):'
          CALL Write_Cartg98(dnx%d0,mole)
          CALL flush_perso(out_unitp)
        END IF
        !=================================================


        !=================================================
        IF (Cart_Transfo_loc) THEN

          IF (debug) write(out_unitp,*) ' calc_CartesianTransfo_new?',Cart_Transfo_loc

          CALL calc_CartesianTransfo_new(dnx,dnx,                      &
                            mole%tab_Cart_transfo(1)%CartesianTransfo, &
                            Qact,nderiv,.TRUE.)

          !=================================================
          IF((WriteCC_loc .OR. debug) .AND. MPI_id==0)THEN
            write(out_unitp,*) ' Cartesian coordinates after the Cartesian Transformation (au):'
            CALL write_dnx(1,mole%ncart,dnx,nderiv_debug)
            write(out_unitp,*) ' Cartesian coordinates after the Cartesian Transformation (ang):'
            CALL Write_Cartg98(dnx%d0,mole)
            CALL flush_perso(out_unitp)
          END IF
          !=================================================

        END IF

        !=================================================
        IF (Gcenter .AND. mole%Centered_ON_CoM) THEN
          icG = mole%ncart-2
          CALL sub3_dncentre_masse(mole%ncart_act,mole%nb_act,mole%ncart, &
                                     dnx,mole%masses,mole%Mtot_inv,icG,nderiv)
        END IF

      END IF

      !=================================================
      ! for partial hessian (pvscf)
      DO ic=1,mole%ncart_act,3
        IF (mole%active_masses(ic) == 0) CALL sub3_dnx_AT1(dnx,ic,nderiv)
      END DO
      !=================================================


!=================================================
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Cartessian coordinates (au)'
        CALL write_dnx(1,mole%ncart,dnx,nderiv_debug)
        write(out_unitp,*) 'END ',name_sub
        write(out_unitp,*)
      END IF
!     -----------------------------------------------------------------
!=================================================

      END SUBROUTINE sub_QactTOdnx_CoordType
      SUBROUTINE sub_QactTOd0x_zmatrix(Qxyz,Qact,mole,Gcenter)
      USE mod_system
      USE mod_dnSVM
      IMPLICIT NONE

      TYPE (zmatrix),     intent(in)     :: mole
      real (kind=Rkind),  intent(in)     :: Qact(:)
      real (kind=Rkind),  intent(inout)  :: Qxyz(mole%ncart_act)
      logical,            intent(in)     :: Gcenter

      CALL sub_QactTOd0x_CoordType(Qxyz,Qact,mole%CoordType,Gcenter)

      END SUBROUTINE sub_QactTOd0x_zmatrix
      SUBROUTINE sub_QactTOd0x_CoordType(Qxyz,Qact,mole,Gcenter)
      USE mod_system
      USE mod_dnSVM
      IMPLICIT NONE

      TYPE (CoordType),   intent(in)     :: mole
      real (kind=Rkind),  intent(in)     :: Qact(:)
      real (kind=Rkind),  intent(inout)  :: Qxyz(mole%ncart_act)
      logical,            intent(in)     :: Gcenter


      TYPE (Type_dnVec) :: dnx
      integer :: nderiv

!     -----------------------------------------------------------------
      integer :: nderiv_debug = 0
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!     -----------------------------------------------------------------
      nderiv = 0
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING sub_QactTOd0x'
        write(out_unitp,*) 'nderiv',nderiv
        write(out_unitp,*) 'ncart',mole%ncart
        write(out_unitp,*) 'Qact =',Qact
        write(out_unitp,*)
        CALL Write_CoordType(mole)
        write(out_unitp,*)
      END IF
!     -----------------------------------------------------------------


      CALL alloc_dnSVM(dnx,mole%ncart,mole%nb_act,nderiv)

      CALL sub_QactTOdnx(Qact,dnx,mole,nderiv,Gcenter)


!=================================================
!=================================================
!     d0x => Qxyz

      Qxyz(:) = dnx%d0(1:mole%ncart_act)
!=================================================
!=================================================



!=================================================
!=================================================
!     -----------------------------------------------------------------
      IF (debug) THEN
       CALL write_dnx(1,mole%ncart,dnx,nderiv_debug)
       write(out_unitp,*) 'END QTOd0x'
       write(out_unitp,*)
      END IF
!     -----------------------------------------------------------------
!=================================================
      CALL dealloc_dnSVM(dnx)

      END SUBROUTINE sub_QactTOd0x_CoordType

      SUBROUTINE sub_d0xTOQact(Qxyz,Qact,mole)
      USE mod_system
      USE mod_dnSVM
      USE mod_Tnum
      IMPLICIT NONE


      TYPE (CoordType)  :: mole
      real (kind=Rkind) :: Qxyz(:)
      real (kind=Rkind) :: Qact(:)

      real (kind=Rkind) :: Rot_initial(3,3),Qat1(3)


      !- working variables -------------------------
      integer :: it,nb_act,i,ic1,ic
      integer :: it_QoutRead,it_QinRead

      TYPE (Type_dnVec) :: dnQin,dnQout

      !-----------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='sub_d0xTOQact'
      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'mole%nb_Qtransfo',mole%nb_Qtransfo
        write(out_unitp,*) 'Qxyz =',Qxyz(:)
        write(out_unitp,*)
        !CALL Write_CoordType(mole)
        write(out_unitp,*)
      END IF
      !-----------------------------------------------------------------

      ! since it is going from out to in, it is better to use it_QoutRead (= it_QinRead+1)
      it_QinRead  = 0 ! Qread=Qxyz => it_QinRead=0
      it_QoutRead = it_QinRead + 1

        it = it_QoutRead
        nb_act = mole%tab_Qtransfo(it_QoutRead)%nb_act


        CALL alloc_dnSVM(dnQout,mole%tab_Qtransfo(it)%nb_Qout,nb_act,0)

        dnQout%d0(1:size(Qxyz)) = Qxyz(:)

        DO it=it_QoutRead,mole%nb_Qtransfo

          CALL alloc_dnSVM(dnQin,mole%tab_Qtransfo(it)%nb_Qin,nb_act,0)

          IF (debug) THEN
            CALL Write_d0Q(it,'Qout ' // trim(adjustl(mole%tab_Qtransfo(it)%name_transfo)),dnQout%d0,6)
            write(out_unitp,*) 'Qout ',it,mole%tab_Qtransfo(it)%name_transfo,dnQout%d0
            CALL flush_perso(out_unitp)
          END IF

          CALL calc_Qtransfo(dnQin,dnQout,mole%tab_Qtransfo(it),0,inTOout=.FALSE.)

          IF (debug) THEN
            CALL Write_d0Q(it,'Qin  ' // trim(adjustl(mole%tab_Qtransfo(it)%name_transfo)),dnQin%d0,6)
            CALL flush_perso(out_unitp)
          END IF

          CALL dealloc_dnSVM(dnQout)
          CALL alloc_dnSVM(dnQout,mole%tab_Qtransfo(it)%nb_Qin,nb_act,0)

          CALL sub_dnVec1_TO_dnVec2(dnQin,dnQout,nderiv=0)
          CALL dealloc_dnSVM(dnQin)

        END DO

        Qact(:) = dnQout%d0(1:size(Qact))
        CALL dealloc_dnSVM(dnQout)



      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Qact',Qact(:)
        write(out_unitp,*) 'END ',name_sub
        write(out_unitp,*)
      END IF
      CALL flush_perso(out_unitp)
      !-----------------------------------------------------------------

      END SUBROUTINE sub_d0xTOQact

      SUBROUTINE Write_d0Q(it,name_info,d0Q,iblock)
      USE mod_system
      IMPLICIT NONE


      integer, intent(in)           :: it,iblock
      character (len=*), intent(in) :: name_info
      real (kind=Rkind), intent(in) :: d0Q(:)

!     - working variables -------------------------
      integer           :: i,iend

      !-----------------------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='Write_d0Q'
      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
      END IF
      !-----------------------------------------------------------------

      IF(keep_MPI) THEN
        write(out_unitp,*) '-----------------------------------------'
        DO i=1,size(d0Q),iblock
          iend = min(size(d0Q),i+iblock-1)
          write(out_unitp,'(a,a,i0,1x,6(1x,f0.4))') name_info,',it_Qtransfo: ',it,d0Q(i:iend)
        END DO
        write(out_unitp,*) '-----------------------------------------'
        CALL flush_perso(out_unitp)
      ENDIF

      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
        write(out_unitp,*)
      END IF
      !-----------------------------------------------------------------

      END SUBROUTINE Write_d0Q

!=======================================================================================
      SUBROUTINE Write_Q_WU(Q,name_Q,type_Q,info)
      USE mod_system
      IMPLICIT NONE

      real (kind=Rkind), intent(in)           :: Q(:)
      integer, intent(in)                     :: type_Q(:)
      character (len=Name_len), intent(in)    :: name_Q(:)
      character (len=*), intent(in), optional :: info

!     - working variables -------------------------
      integer           :: i
      TYPE (REAL_WU)    :: QWU

      !-----------------------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='Write_Q_WU'
      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'i,name_Q(i),type_Q(i),Q(i)'
        DO i=1,size(Q)
          write(out_unitp,*) i,name_Q(i),type_Q(i),Q(i)
        END DO
      END IF
      !-----------------------------------------------------------------

      write(out_unitp,*) '-----------------------------------------'
      IF (present(info) .AND. MPI_id==0) write(out_unitp,*) info

      DO i=1,size(Q)

        SELECT CASE (type_Q(i))
        CASE (-3)
          QWU = REAL_WU(acos(Q(i)),'rad','angle')
        CASE (3,4)
          QWU = REAL_WU(Q(i),'rad','angle')
        CASE (1,2)
          QWU = REAL_WU(Q(i),'bohr','L')
        CASE default
          QWU = REAL_WU(Q(i),'','no_dim')
        END SELECT

        write(out_unitp,'(a,i0,5x,a,5x,a)') name_Q(i),i,             &
                  RWU_Write(QWU,WithUnit=.TRUE.,WorkingUnit=.TRUE.),&
                  RWU_Write(QWU,WithUnit=.TRUE.,WorkingUnit=.FALSE.)
      END DO
      write(out_unitp,*) '-----------------------------------------'
      CALL flush_perso(out_unitp)

      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
        write(out_unitp,*)
      END IF
      !-----------------------------------------------------------------

      END SUBROUTINE Write_Q_WU
!=======================================================================================

      SUBROUTINE Get_Qread(Q,name_Q,type_Q,read_nameQ,unit,read_xyz0,info)
      USE mod_system
      IMPLICIT NONE


      real (kind=Rkind), intent(inout)        :: Q(:)
      character (len=Name_len), intent(inout) :: name_Q(:)
      integer, intent(in)                     :: type_Q(:)
      logical, intent(in)                     :: read_nameQ,read_xyz0
      character (len=Name_len), intent(in)    :: unit

      character (len=*), intent(in), optional :: info

      !- working variables -------------------------
      integer           :: i,k,err_ioQ
      TYPE (REAL_WU)    :: QWU
      character (len=Line_len) :: Read_name
      character (len=Name_len) :: unit_Q

      !-----------------------------------------------------------------
      integer :: err_mem,memory,err_io
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='Get_Qread'
      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'type_Q',type_Q
        write(out_unitp,*) 'unit: ',unit
        write(out_unitp,*) 'read_xyz0',read_xyz0
        write(out_unitp,*) 'info',info

      END IF
      !-----------------------------------------------------------------

      IF (read_xyz0) THEN

        DO i=1,size(Q)/3

          read(in_unitp,*,IOSTAT=err_io) name_Q(3*i-2),Q(3*i-2:3*i)
          !write(out_unitp,*) name_Q(3*i-2),Q(3*i-2:3*i)*.52d0

          IF (err_io /= 0) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) '  while reading the Cartessian reference geometry ...'
            write(out_unitp,*) '   ... just after the namelist "minimum".'
            write(out_unitp,'(a,i0,a,i0,a)') '  Trying to read the atom:',i,' among ',size(Q)/3,'.'
            write(out_unitp,*) ' Check your data !!'
            STOP
          END IF

          name_Q(3*i-0) = "Z_" // trim(adjustl(name_Q(3*i-2)))
          name_Q(3*i-1) = "Y_" // trim(adjustl(name_Q(3*i-2)))
          name_Q(3*i-2) = "X_" // trim(adjustl(name_Q(3*i-2)))
        END DO

        ! conversion of unit if needed
        IF (unit == 'angs' ) THEN
          DO i=1,size(Q)
            SELECT CASE (type_Q(i))
            CASE (3,4)
              QWU = REAL_WU(Q(i),'°',    'angle')
            CASE (1,2)
              QWU = REAL_WU(Q(i),'Angs', 'L')
            CASE default
              QWU = REAL_WU(Q(i),'',     'no_dim')
            END SELECT
            Q(i) = convRWU_TO_R_WITH_WorkingUnit(QWU)

            !write(out_unitp,*) 'i,QWU, conv',i,QWU,Q(i)


          END DO
        END IF

      ELSE
        DO i=1,size(Q)
           ! read the first word: it can be the variable name or its value
           CALL read_name_advNo(in_unitp,Read_name,err_io)
           !write(6,*) i,'Read_name: ',Read_name
           ! try to read its value
           read(Read_name,*,IOSTAT=err_ioQ) Q(i)
           !write(6,*) i,'Read_name: ',Read_name,'err_ioQ',err_ioQ

           IF (err_ioQ /= 0) THEN ! an error, it should be the variable name or a true error
             name_Q(i) = trim(adjustl(Read_name))

             !write(6,*) i,'name_Q(i): ',name_Q(i)

             !now we read the value
             CALL read_name_advNo(in_unitp,Read_name,err_io)
             !write(6,*) i,'Read_name: ',Read_name
             read(Read_name,*,IOSTAT=err_ioQ) Q(i)
             !write(6,*) i,'Read_name: ',Read_name,'err_ioQ',err_ioQ ; flush(6)
             IF (err_ioQ /= 0) THEN
               write(out_unitp,*) ' ERROR in ',name_sub
               write(out_unitp,*) '  while reading the curvilinear reference geometry '
               write(out_unitp,*) '   ... just after the namelist "minimum"'
               write(out_unitp,*) ' error with the value name: ',trim(adjustl(Read_name))
               write(out_unitp,*) ' the variable name:         ',trim(adjustl(name_Q(i)))
               write(out_unitp,*) ' Check your data !!'
               STOP
             END IF
           END IF
           ! Normally, the value is readed. Try to read the unit
           IF (err_io < 0) THEN ! end-of-line ?
             Read_name = ''
           ELSE
             CALL read_name_advNo(in_unitp,Read_name,err_io)
           END IF

           IF(MPI_id==0) THEN
             write(out_unitp,*) i,name_Q(i),':',Q(i),':',trim(adjustl(Read_name))
             write(out_unitp,*) i,'type_Q(i) :',type_Q(i)
           ENDIF

           IF (len_trim(Read_name) > 0) THEN
             IF (trim(adjustl(Read_name)) == '°' .OR.                           &
                 trim(adjustl(Read_name)) == 'rad') THEN
               QWU = REAL_WU(Q(i),trim(adjustl(Read_name)),'angle')
               IF (type_Q(i) /= 3 .AND. type_Q(i) /= 4 .AND. type_Q(i) /= 0) THEN
                 write(out_unitp,*) ' ERROR in ',name_sub
                 write(out_unitp,*) '  The unit and type_Q(i) are incompatible.'
                 write(out_unitp,*) '    unit     ',trim(adjustl(Read_name))
                 write(out_unitp,*) '    type_Q(i)',type_Q(i)
                 write(out_unitp,*) ' The compatible values are:'
                 write(out_unitp,*) '  angs or bohr => type_Q(i)=1,2 or 0'
                 write(out_unitp,*) '  ° or rad     => type_Q(i)=3,4 or 0'
                 write(out_unitp,*) ' Check your data !!'
                 STOP 'ERROR in Get_Qread: Wrong unit'
               END IF
             ELSE IF (trim(adjustl(Read_name)) == 'Angs' .OR.                   &
                      trim(adjustl(Read_name)) == 'angs' .OR.                   &
                      trim(adjustl(Read_name)) == 'bohr') THEN
               QWU = REAL_WU(Q(i),trim(adjustl(Read_name)),'L')
               IF (type_Q(i) /= 1 .AND. type_Q(i) /= 2 .AND. type_Q(i) /= 0) THEN
                 write(out_unitp,*) ' ERROR in ',name_sub
                 write(out_unitp,*) '  The unit and type_Q(i) are incompatible.'
                 write(out_unitp,*) '    unit     ',trim(adjustl(Read_name))
                 write(out_unitp,*) '    type_Q(i)',type_Q(i)
                 write(out_unitp,*) ' The compatible values are:'
                 write(out_unitp,*) '  angs or bohr => type_Q(i)=1,2 or 0'
                 write(out_unitp,*) '  ° or rad     => type_Q(i)=3,4 or 0'
                 write(out_unitp,*) ' Check your data !!'
                 STOP 'ERROR in Get_Qread: Wrong unit'
               END IF
             ELSE
               write(out_unitp,*) ' ERROR in ',name_sub
               write(out_unitp,*) '  The unit is wrong: ',trim(adjustl(Read_name))
               write(out_unitp,*) ' The possible values are:'
               write(out_unitp,*) '  angs or bohr'
               write(out_unitp,*) '  ° or rad'
               write(out_unitp,*) ' Check your data !!'
               STOP 'ERROR in Get_Qread: Wrong unit'
             END IF
             ! write(6,*) 'read with unit:'
             ! SELECT CASE (type_Q(i))
             ! CASE (3,4)
             !   QWU = REAL_WU(Q(i),trim(adjustl(Read_name)),'angle')
             ! CASE (1,2)
             !   QWU = REAL_WU(Q(i),trim(adjustl(Read_name)),'L')
             ! CASE default
             !   QWU = REAL_WU(Q(i),trim(adjustl(Read_name)),'no_dim')
             ! END SELECT

           ELSE IF (unit == 'angs' ) THEN ! angs + degree
             SELECT CASE (type_Q(i))
             CASE (3,4)
               QWU = REAL_WU(Q(i),'°',    'angle')
             CASE (1,2)
               QWU = REAL_WU(Q(i),'Angs', 'L')
             CASE default
               QWU = REAL_WU(Q(i),'',     'no_dim')
             END SELECT

          ELSE ! bohr + radian
             SELECT CASE (type_Q(i))
             CASE (3,4)
               QWU = REAL_WU(Q(i),'rad',  'angle')
             CASE (1,2)
               QWU = REAL_WU(Q(i),'bohr', 'L')
             CASE default
               QWU = REAL_WU(Q(i),'',     'no_dim')
             END SELECT

          END IF

          Q(i) = convRWU_TO_R_WITH_WorkingUnit(QWU)
          IF(MPI_id==0) write(out_unitp,*) i,QWU,'working value Q(i)',Q(i)

        END DO
      END IF

      !-----------------------------------------------------------------
      IF (debug) THEN
        IF (present(info)) THEN
          CALL Write_Q_WU(Q,name_Q,type_Q,info)
        ELSE
          CALL Write_Q_WU(Q,name_Q,type_Q)
        END IF
        write(out_unitp,*) 'END ',name_sub
        write(out_unitp,*)
      END IF
      !-----------------------------------------------------------------

      END SUBROUTINE Get_Qread


!================================================================
!       Write Cartesian coordinates (for gaussian)
!================================================================
      SUBROUTINE Write_Cartg98_zmatrix(d0x,mole)
      USE mod_system
      IMPLICIT NONE

      TYPE (zmatrix),    intent(in) :: mole
      real (kind=Rkind), intent(in) :: d0x(:)

        CALL Write_Cartg98_CoordType(d0x,mole%CoordType)

      END SUBROUTINE Write_Cartg98_zmatrix

      SUBROUTINE Write_Cartg98_CoordType(d0x,mole)
      USE mod_system
      IMPLICIT NONE

      TYPE (CoordType),  intent(in) :: mole
      real (kind=Rkind), intent(in) :: d0x(:)

      real (kind=Rkind) :: a0
      integer           :: Z_act(mole%nat)

      integer       :: i,iZ,ncart

      ncart = size(d0x)


      Z_act(:) = -1
      iZ = 0
      DO i=1,mole%nat
        IF (mole%Z(i) > 0) THEN
          iZ = iZ + 1
          Z_act(iZ) = mole%Z(i)
        END IF
      END DO

      a0 = get_Conv_au_TO_unit("L","Angs")

      iZ = 0
      write(out_unitp,*) '=============================================='
      write(out_unitp,*) '= Gaussian CC ================================'
      DO i=1,ncart,3
        iZ = iZ + 1
        write(out_unitp,111) Z_act(iZ),0,d0x(i+0:i+2)*a0
 111    format(1x,2(1x,i5),3(2x,f20.9))

      END DO
      write(out_unitp,*) '= END Gaussian CC ============================'
      write(out_unitp,*) '=============================================='

      write(out_unitp,*) '=============================================='
      write(out_unitp,*) '= XYZ format ================================='
      write(out_unitp,*) ncart/3
      write(out_unitp,*)

      iZ = 0
      DO i=1,ncart,3
        iZ = iZ + 1
        write(out_unitp,112) Z_act(iZ),d0x(i+0)*a0,d0x(i+1)*a0,d0x(i+2)*a0
 112    format(2x,i5,3(2x,f20.9))

      END DO
      write(out_unitp,*) '= END XYZ format ============================='
      write(out_unitp,*) '=============================================='

      END SUBROUTINE Write_Cartg98_CoordType

!================================================================
!       Write Cartesian coordinates (xyz format)
!================================================================
      SUBROUTINE Write_XYZ(d0x,mole,unit,io_unit)
      USE mod_system
      IMPLICIT NONE

      TYPE (CoordType),            intent(in) :: mole
      real (kind=Rkind),           intent(in) :: d0x(mole%ncart)
      character (len=*), optional, intent(in) :: unit
      integer,           optional, intent(in) :: io_unit

      real (kind=Rkind) :: a0
      integer           :: Z_act(mole%nat)

      integer       :: i,iZ,io_unit_loc

      IF (present(io_unit)) THEN
        io_unit_loc = io_unit
      ELSE
        io_unit_loc = out_unitp
      END IF


      Z_act(:) = -1
      iZ = 0
      DO i=1,mole%nat
        IF (mole%Z(i) > 0) THEN
          iZ = iZ + 1
          Z_act(iZ) = mole%Z(i)
        END IF
      END DO

      write(io_unit_loc,*) '=============================================='

      IF (present(unit)) THEN
        a0 = get_Conv_au_TO_unit("L",unit)
        write(io_unit_loc,*) '= XYZ format (',unit,') =========================='
      ELSE
        a0 = get_Conv_au_TO_unit("L","Angs")
        write(io_unit_loc,*) '= XYZ format (Angs) =========================='
      END IF

      write(io_unit_loc,*) mole%nat_act
      write(io_unit_loc,*)

      iZ = 0
      DO i=1,mole%ncart_act,3
        iZ = iZ + 1
        write(io_unit_loc,112) Z_act(iZ),d0x(i+0)*a0,d0x(i+1)*a0,d0x(i+2)*a0
 112    format(2x,i5,3(2x,f20.16))
!112    format(2x,i5,3(2x,f20.9)) ! Emil change: Needed more digits for coordinate conversion between Tana and MidasCpp

      END DO
      write(io_unit_loc,*) '= END XYZ format ============================='
      write(io_unit_loc,*) '=============================================='

      END SUBROUTINE Write_XYZ

      SUBROUTINE analyze_dnx_zmatrix(dnx,Qact,mole)
      USE mod_system
      USE mod_dnSVM
      IMPLICIT NONE

      TYPE (Type_dnVec) :: dnx
      TYPE (zmatrix)    :: mole
      real (kind=Rkind) :: Qact(:)

      CALL analyze_dnx_CoordType(dnx,Qact,mole%CoordType)

      END SUBROUTINE analyze_dnx_zmatrix

      SUBROUTINE analyze_dnx_CoordType(dnx,Qact,mole)
      USE mod_system
      USE mod_dnSVM
      IMPLICIT NONE

      TYPE (Type_dnVec) :: dnx
      TYPE (CoordType)    :: mole
      real (kind=Rkind) :: Qact(:)

      integer :: i,j,k
      real (kind=Rkind) :: d
      character (len=*), parameter :: name_sub = 'analyze_dnx'
      logical, parameter :: debug = .FALSE.
!      logical, parameter :: debug = .TRUE.


        CALL check_alloc_dnVec(dnx,'dnx',name_sub)

        write(out_unitp,*) 'BEGINNING in ',name_sub

        IF (3*mole%nat_act > dnx%nb_var_vec) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' 3*nat_act > nb_var_vec',3*mole%nat_act,dnx%nb_var_vec
          write(out_unitp,*) ' Check the fortran !!!!'
          STOP
        END IF

        IF (debug) THEN
          DO i=1,mole%nat_act
          DO j=i+1,mole%nat_act
             d = (dnx%d0(3*i-2)-dnx%d0(3*j-2))**2 +                     &
                 (dnx%d0(3*i-1)-dnx%d0(3*j-1))**2 +                     &
                 (dnx%d0(3*i-0)-dnx%d0(3*j-0))**2
             write(out_unitp,*) 'mass weighted distances: ',i,j,sqrt(d)
          END DO
          END DO
        END IF

        write(out_unitp,*) ' d0x Mass weighted'
        CALL write_dnx(1,dnx%nb_var_vec,dnx,0)

        DO i=1,mole%ncart_act
          dnx%d0(i) = dnx%d0(i)/mole%d0sm(i)
        END DO

        IF (debug) THEN
          DO i=1,mole%nat_act
          DO j=i+1,mole%nat_act
             d = (dnx%d0(3*i-2)-dnx%d0(3*j-2))**2 +                     &
                 (dnx%d0(3*i-1)-dnx%d0(3*j-1))**2 +                     &
                 (dnx%d0(3*i-0)-dnx%d0(3*j-0))**2
             write(out_unitp,*) 'distances: ',i,j,sqrt(d)
          END DO
          END DO
        END IF
        write(out_unitp,*) ' d0x NOT Mass weighted'
        CALL write_dnx(1,dnx%nb_var_vec,dnx,0)
        DO i=1,mole%ncart_act
          dnx%d0(i) = dnx%d0(i)*mole%d0sm(i)
        END DO


        write(out_unitp,*) ' d0x NOT Mass weighted and NOT recentered/CM'
        CALL sub_QactTOdnx(Qact,dnx,mole,0,.FALSE.)
        CALL write_dnx(1,dnx%nb_var_vec,dnx,0)

        write(out_unitp,*) 'END ',name_sub

      END SUBROUTINE analyze_dnx_CoordType


      SUBROUTINE sub_dnFCC_TO_dnFcurvi(Qact,dnFCC,dnFcurvi,mole)

      USE mod_system
      USE mod_dnSVM
      USE mod_Tnum
      IMPLICIT NONE

      TYPE (CoordType),  intent(in)    :: mole
      real (kind=Rkind), intent(in)    :: Qact(:)


      TYPE(Type_dnS), intent(in)    :: dnFCC
      TYPE(Type_dnS), intent(inout) :: dnFcurvi


      real(kind=Rkind)  :: work(mole%nb_act,mole%ncart_act)
      TYPE(Type_dnVec)  :: dnx
      logical           :: Gcenter
      integer           :: i,j,k,l,ncc
      integer           :: nderiv

!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'sub_dnFCC_TO_dnFcurvi'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
        write(out_unitp,*) 'Val, grad and hessian in CC'
        CALL Write_dnSVM(dnFCC)
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

      nderiv = dnFCC%nderiv
      IF (.NOT. dnFcurvi%alloc) CALL alloc_dnS(dnFcurvi,mole%nb_act,nderiv)
      nderiv = min(dnFCC%nderiv,dnFcurvi%nderiv)

      Gcenter = .FALSE.

      CALL alloc_dnSVM(dnx,mole%ncart,mole%nb_act,nderiv)

      CALL sub_QactTOdnx(Qact,dnx,mole,nderiv,Gcenter)

      dnFcurvi%d0 = dnFCC%d0
      ncc = mole%ncart_act

      IF (nderiv >= 1) THEN
        DO i=1,dnFcurvi%nb_var_deriv
          dnFcurvi%d1(i) = dot_product(dnx%d1(1:ncc,i),dnFCC%d1(:))
        END DO
      END IF

      IF (nderiv == 2) THEN

        DO i=1,dnFcurvi%nb_var_deriv
        DO k=1,ncc
          work(i,k) = dot_product(dnx%d1(1:ncc,i),dnFCC%d2(k,:))
        END DO
        END DO

        DO i=1,dnFcurvi%nb_var_deriv
        DO j=1,i
          dnFcurvi%d2(i,j) = dot_product(dnx%d2(1:ncc,i,j),dnFCC%d1(:)) + &
                             dot_product(dnx%d1(1:ncc,j),work(i,:))

          dnFcurvi%d2(j,i) = dnFcurvi%d2(i,j)
        END DO
        END DO


      END IF

      CALL dealloc_dnSVM(dnx)
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'E, grad and hessian in zmt'
        CALL Write_dnSVM(dnFcurvi)
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

      END SUBROUTINE sub_dnFCC_TO_dnFcurvi


      SUBROUTINE Set_paramQ_FOR_optimization(Qact,mole,Set_Val)
      USE mod_system
      IMPLICIT NONE


!----- for the CoordType and Tnum --------------------------------------
      real (kind=Rkind), intent(inout) :: Qact(:)
      TYPE (CoordType),  intent(inout) :: mole
      integer, intent(in)              :: Set_Val

      integer :: nopt,ib,i_RVec,i,i1,i2,iQdyn
!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'Set_paramQ_FOR_optimization'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'i_OptParam ',para_FOR_optimization%i_OptParam
        write(out_unitp,*) 'Qact',Qact
      END IF
!---------------------------------------------------------------------
      IF (count(mole%opt_Qdyn > 0) < 1 ) RETURN

      nopt = mole%nb_act1
      i1 = para_FOR_optimization%i_OptParam+1
      i2 = para_FOR_optimization%i_OptParam+nopt
      IF (debug) write(out_unitp,*) 'nopt geometry',nopt
      IF (nopt > 0) THEN
        para_FOR_optimization%nb_OptParam =                             &
                                para_FOR_optimization%nb_OptParam + nopt
        !write(out_unitp,*) ' size opt_Qdyn',size(mole%opt_Qdyn)
        DO i=1,size(mole%opt_Qdyn)
          IF (mole%opt_Qdyn(i) == 1) mole%opt_Qdyn(i) = 5
        END DO

        IF (Set_Val == -1) THEN

          i = 0
          DO iQdyn=1,size(mole%opt_Qdyn)
            IF (mole%opt_Qdyn(iQdyn) /= 0) THEN
              i = i + 1
              IF (i > nopt) THEN
                write(out_unitp,*) ' ERROR in ',name_sub
                write(out_unitp,*) '  the value of i is larger than  nopt',i,nopt
                write(out_unitp,*) ' Check the source !!'
                STOP
              END IF
              para_FOR_optimization%Val_RVec(i1+i-1) = Qact(i)
              para_FOR_optimization%Opt_RVec(i1+i-1) = mole%opt_Qdyn(iQdyn)
            END IF
          END DO

        ELSE IF (Set_Val == 1) THEN
          i = 0
          DO iQdyn=1,size(mole%opt_Qdyn)
            IF (mole%opt_Qdyn(iQdyn) /= 0) THEN
              i = i + 1
              IF (i > nopt) THEN
                write(out_unitp,*) ' ERROR in ',name_sub
                write(out_unitp,*) '  the value of i is larger than  nopt',i,nopt
                write(out_unitp,*) ' Check the source !!'
                STOP
              END IF
              Qact(i) = para_FOR_optimization%Val_RVec(i1+i-1)
            END IF
          END DO
        END IF
        para_FOR_optimization%i_OptParam = i2

      END IF

      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'nb_OptParam ',para_FOR_optimization%nb_OptParam
        write(out_unitp,*) 'Val_RVec ',para_FOR_optimization%Val_RVec

        write(out_unitp,*) 'Qact ',Qact
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF

      END SUBROUTINE Set_paramQ_FOR_optimization

END MODULE mod_paramQ
