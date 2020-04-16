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
MODULE mod_OTF
   USE mod_system
   USE mod_dnSVM
   use mod_OTF_def,    only: assignment(=),param_otf
   use mod_PrimOp_def, only: param_pes, write_param_pes
   USE mod_Constant
   IMPLICIT NONE

   PRIVATE
   PUBLIC   dnOp_grid_OnTheFly, Read_GradHess_Molpro, Read_dnDipCC_Gauss
   PUBLIC   Read_dnPolarizabilityCC_Gauss
   PUBLIC   Read_hess_Fchk

   CONTAINS

!================================================================
!    subroutine enables to calculate the energy, gradient and hessian
!     On-the-fly (with gaussian or gamess)
!
!    You should not use this subroutine directly.
!    Instead, one must use the "sub_Operator/dnOp_grid" subroutine
!
!================================================================
      SUBROUTINE dnOp_grid_OnTheFly(Qxyz,MatdnECC,nderivE,              &
                                    MatdnMuCC,nderivMu,                 &
                                    mole,para_PES)
      use mod_Coord_KEO,  only: assignment(=),CoordType

      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole

!----- for Qdyn Qact ... ---------------------------------------------
      real (kind=Rkind)   :: Qxyz(mole%ncart_act)
      TYPE (param_PES)    :: para_PES

!----- input output variables ----------------------------------------
      integer        :: nderivE,nderivMu
      TYPE(Type_dnS) :: MatdnECC(para_PES%nb_elec,para_PES%nb_elec)
      TYPE(Type_dnS) :: MatdnMuCC(para_PES%nb_elec,para_PES%nb_elec,3)

      integer        :: err_calc


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='dnOp_grid_OnTheFly'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Qxyz'
        write(out_unitp,'(3(1x,f16.6))') Qxyz
        write(out_unitp,*) 'nderivE,nderivMu',nderivE,nderivMu
        write(out_unitp,*)
      END IF
!-----------------------------------------------------------

      !----------------------------------------------------------------
      IF (.NOT. para_PES%OnTheFly) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) ' This subroutine works only with on-the-fly calculation'
        write(out_unitp,*) ' It should never append!!'
        write(out_unitp,*) '  CHECK the fortran'
        STOP
      END IF
      !----------------------------------------------------------------

      !----------------------------------------------------------------
      SELECT CASE (para_PES%para_OTF%ab_initio_prog)
      CASE ('g03','g09')
          IF (para_PES%nb_elec /= 1) STOP 'Yet we connot use gaussian whit nb_elec>1'
          IF (debug) write(out_unitp,*) 'With pot_mu_onthefly_gauss'
          CALL pot_mu_onthefly_gauss(Qxyz,MatdnECC(1,1),nderivE,        &
                                     MatdnMuCC(1,1,:),nderivMu,         &
                                     mole,para_PES,err_calc)
          IF (err_calc /= 0) MatdnECC(1,1)%d0 = para_PES%pot0 + ONE
      CASE ('gamess','gamess2014')
          IF (para_PES%nb_elec /= 1) STOP 'Yet we connot use gamess whit nb_elec>1'
          IF (debug) write(out_unitp,*) 'With pot_mu_onthefly_gamess'
          CALL pot_mu_onthefly_gamess(Qxyz,MatdnECC(1,1),nderivE,       &
                                      MatdnMuCC(1,1,:),nderivMu,        &
                                      mole,para_PES)
      CASE ('generic')
          IF (debug) write(out_unitp,*) 'With pot_mu_onthefly_generic'
          CALL onthefly_generic(Qxyz,MatdnECC,nderivE,                  &
                                MatdnMuCC,nderivMu,                     &
                                mole,para_PES)
      CASE default ! ERROR: wrong program !
          CALL write_param_PES(para_PES)
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' The ab initio program is UNKNOWN ',      &
                      trim(para_PES%para_OTF%ab_initio_prog)
        STOP
      END SELECT
      !----------------------------------------------------------------

IF (.NOT. para_PES%Read_OnTheFly_only) THEN
  ! remove the ab-initio files
  CALL file_delete(para_PES%para_OTF%file_data)
  CALL file_delete(para_PES%para_OTF%file_log)
  CALL file_delete(para_PES%para_OTF%file_FChk)
  CALL file_delete(para_PES%para_OTF%file_pun)
END IF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'MatdnECC'
        CALL Write_MatOFdnS(MatdnECC)
        IF (nderivMu > -1) THEN
          write(out_unitp,*) 'MatdnMuCC(:,:,1): MatDipX'
          CALL Write_MatOFdnS(MatdnMuCC(:,:,1))
          write(out_unitp,*) 'MatdnMuCC(:,:,2): MatDipY'
          CALL Write_MatOFdnS(MatdnMuCC(:,:,2))
          write(out_unitp,*) 'MatdnMuCC(:,:,3): MatDipZ'
          CALL Write_MatOFdnS(MatdnMuCC(:,:,3))
        END IF
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE dnOp_grid_OnTheFly
      SUBROUTINE pot_mu_onthefly_gauss(Qxyz,dnECC,nderivE,dnMuCC,nderivMu,&
                                        mole,para_PES,err_calc)
      use mod_Coord_KEO,  only: assignment(=),CoordType

      USE mod_system
      USE mod_dnSVM
      USE mod_PrimOp_def
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType)    :: mole

!----- for Qdyn Qact ... ---------------------------------------------
      real (kind=Rkind) :: Qxyz(mole%ncart_act)
      TYPE (param_PES)  :: para_PES
      integer, optional :: err_calc


!----- input output variables ----------------------------------------
      integer           :: nderivE,nderivMu
      TYPE(Type_dnS)    :: dnECC,dnMuCC(3)

!----- working variables ----------------------------------------
      integer       :: err
      integer       :: i


!----- function -------------------------------------------------
!----------------------------------------------------------------

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='pot_mu_onthefly_gauss'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING',name_sub
        write(out_unitp,*) 'Qxyz',Qxyz(:)
        write(out_unitp,*) 'nderivE,nderivMu',nderivE,nderivMu
        write(out_unitp,*) 'Read_OnTheFly_only',para_PES%Read_OnTheFly_only
        write(out_unitp,*)
      END IF
!-----------------------------------------------------------

      !-----------------------------------------------------------------
      IF (.NOT. para_PES%Read_OnTheFly_only) THEN
        IF (present(err_calc)) THEN
         CALL Calc_EneDip_WITH_gauss(Qxyz,nderivE,nderivMu,            &
                                        mole,para_PES,para_PES%para_OTF,err_calc)
        ELSE
          CALL Calc_EneDip_WITH_gauss(Qxyz,nderivE,nderivMu,            &
                                        mole,para_PES,para_PES%para_OTF)
        END IF
      END IF
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !- read the energy from the file energy, gradient, hessian
      CALL Read_dnECC_Gauss(dnECC,para_PES%para_OTF%file_FChk%name,     &
                                                 nderivE,mole%ncart_act)

      !-----------------------------------------------------------------
      !- read the Dipole Moment from the file Test.FChk
      IF (nderivMu > -1) THEN
        CALL Read_dnDipCC_Gauss(dnMuCC,para_PES%para_OTF%file_FChk%name,&
                                nderivMu,mole%ncart_act)
      END IF
      !-----------------------------------------------------------------



!-----------------------------------------------------------
      IF (debug) THEN
         write(out_unitp,*) 'dnECC'
         CALL Write_dnSVM(dnECC)
         IF (nderivMu > -1) THEN
           write(out_unitp,*) 'dnMuCC(1): DipX'
           CALL Write_dnSVM(dnMuCC(1))
           write(out_unitp,*) 'dnMuCC(2): DipY'
           CALL Write_dnSVM(dnMuCC(2))
           write(out_unitp,*) 'dnMuCC(:): DipZ'
           CALL Write_dnSVM(dnMuCC(3))
         END IF
         write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE pot_mu_onthefly_gauss

      SUBROUTINE Calc_EneDip_WITH_gauss(Qxyz,nderivE,nderivDip,&
                                        mole,para_PES,para_OTF,err_calc)

      USE mod_system
      use mod_Coord_KEO,  only: assignment(=),CoordType
      USE mod_PrimOp_def
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType)    :: mole

!----- for Qdyn Qact ... ---------------------------------------------
      real (kind=Rkind) :: Qxyz(mole%ncart_act)
      TYPE (param_PES)  :: para_PES
      TYPE (param_OTF)  :: para_OTF
      integer, optional :: err_calc

!----- for the files -------------------------------------------------
      integer :: nio,nio_header,nio_footer

!----- input output variables ----------------------------------------
      integer        :: nderivE,nderivDip

!----- working variables ----------------------------------------
      integer       :: err,err_calc_loc
      integer                       :: i,iq,Z

      character (len=Name_longlen)  :: ab_level
      character (len=Name_longlen)  :: ab_level_temp

      character (len=Name_longlen)  :: labelR
      logical                       :: located
      character (len=Line_len)      :: line

!----- function -------------------------------------------------
!----------------------------------------------------------------

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Calc_EneDip_WITH_gauss'
      logical, parameter :: debug=.FALSE.
!      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING',name_sub
        write(out_unitp,*) 'Qxyz',Qxyz(:)
        write(out_unitp,*) 'nderivE,nderivDip',nderivE,nderivDip
        write(out_unitp,*)
      END IF
      !----------------------------------------------------------------


      !----------------------------------------------------------------
      !- gaussian input file ------------------------------------------
      CALL file_open(para_OTF%file_data,nio)

      !----------------------------------------------------------------
      !- add a header ------------------------------------------
      IF (para_OTF%header) THEN
        CALL file_open(para_OTF%file_header,nio_header)

        DO
          read(nio_header,21,end=998) line
 21       format(132A)
          write(nio,*) trim(adjustl(line))
        END DO
 998    close(nio_header)
      END IF
      !----------------------------------------------------------------


      ab_level = trim(para_OTF%ab_initio_meth) // ' ' //                &
                                          trim(para_OTF%ab_initio_basis)
      ab_level_temp = ab_level
      CALL string_uppercase_TO_lowercase(ab_level_temp)

      write(nio,*) '%chk=',trim(para_OTF%file_name)

      IF (nderivE == 0) write(nio,*) '#n ',ab_level
      IF (nderivE == 1) write(nio,*) '#n force ',ab_level
      IF (nderivE == 2) write(nio,*) '#n freq=noraman ',ab_level

      IF (index(ab_level_temp,'conver') == 0) THEN ! check if scf=conver=...) is present
        write(nio,*) '# unit=(au,rad) nosymm scf=tight FormCheck=All'
      ELSE
        write(nio,*) '# unit=(au,rad) nosymm FormCheck=All'
      END IF
      write(nio,*) '# integral=(grid=ultrafine) '
      IF (nderivDip > -1) write(nio,*) '# density=current'
      write(nio,*) '# iop(1/18=40,2/9=111,2/11=1,7/33=1,1/33=1)'
      write(nio,*) '# iop(3/27=30,11/27=30)'

      write(nio,*) ' '
      write(nio,*) ' xx '
      write(nio,*) ' '
      write(nio,*) mole%charge," ",mole%multiplicity

      iq=0
      DO i=1,mole%nat_act
        Z = mole%Z(i)
        IF (mole%Z(i) == 0) Z=-1
        write(nio,13) Z,0,Qxyz(iq+1:iq+3)
 13     format(i5,1x,i5,3(1x,f20.15))
        iq = iq+3
      END DO
      write(nio,*) ' '

      IF (para_OTF%footer) THEN
        CALL file_open(para_OTF%file_footer,nio_footer)

        DO
          read(nio_footer,22,end=997) line
 22       format(132A)
          write(nio,*) trim(adjustl(line))
        END DO
 997    close(nio_footer)
      END IF

      close(nio)
      !----------------------------------------------------------------



      !----------------------------------------------------------------
      !- gaussian execution -------------------------------------------
      CALL file_delete(para_OTF%file_log)

      !CALL system(para_OTF%commande_unix)
      CALL EXECUTE_COMMAND_LINE(para_OTF%commande_unix)


      located = .FALSE.
      CALL file_open(para_OTF%file_log,nio,append=.TRUE.)
      backspace(nio,err=999)
      read(nio,'(a32)',err=999) labelR
      !write(out_unitp,*) 'last line: ',labelR
      located = verify(labelR,' Normal termination of Gaussian') == 0
      close(nio)

 999  CONTINUE
      IF (present(err_calc)) THEN
        IF (.NOT. located) THEN
          err_calc = 1
        ELSE
          err_calc = 0
        END IF
      ELSE
        IF (.NOT. located) THEN
          write(out_unitp,*) 'ERROR in the gaussian execution'
          write(out_unitp,*) 'no line: "Normal termination of Gaussian"'
          write(out_unitp,*) 'last line: ',labelR
          write(out_unitp,*) 'unix command:',para_OTF%commande_unix
          STOP
        END IF
      END IF

      IF (debug) THEN
        write(out_unitp,*) 'END',name_sub
      END IF

      END SUBROUTINE Calc_EneDip_WITH_gauss

      SUBROUTINE Read_hess_Fchk(dnFCC,fchk_name,nderiv,ncart_act)

      USE mod_system
      USE mod_dnSVM
      implicit none

!----- input output variables ----------------------------------------
      integer                  :: nderiv,ncart_act
      character (len=Line_len) :: fchk_name
      TYPE (param_file)        :: file_FChk
      TYPE(Type_dnS)           :: dnFCC

!----- for the files -------------------------------------------------
      integer :: nio,nio_header


      logical                    :: zmt

      integer       :: i,j,k,l,vi,vj,iq

      character (len=Name_len)  :: ab_level
      character (len=Name_longlen)  :: labelR
      integer             :: Z
      character (len=Name_len)  :: name1_i,name2_i,name3_i
      integer             :: i1,i2,i3
      logical             :: located
      integer             :: err

      character (len=Line_len) :: line


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Read_hess_Fchk'
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------

      IF (debug)  write(out_unitp,*) 'BEGINNING',name_sub

      CALL Read_dnECC_Gauss(dnFCC,fchk_name,nderiv,ncart_act)

      IF (debug) write(out_unitp,*) 'END',name_sub

      END SUBROUTINE Read_hess_Fchk
      SUBROUTINE Read_dnECC_Gauss(dnECC,fchk_name,nderiv,ncart_act)

      USE mod_system
      USE mod_dnSVM
      implicit none

!----- input output variables ----------------------------------------
      integer                  :: nderiv,ncart_act
      character (len=Line_len) :: fchk_name
      TYPE (param_file)        :: file_FChk
      TYPE(Type_dnS)           :: dnECC

!----- for the files -------------------------------------------------
      integer :: nio

      integer                   :: i,j
      character (len=Name_len)  :: name1_i
      logical                   :: located
      integer                   :: err


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Read_dnECC_Gauss'
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING',name_sub
        write(out_unitp,*) 'nderiv',nderiv
        write(out_unitp,*) 'fchk_name: ',fchk_name
        write(out_unitp,*)
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

      CALL alloc_dnSVM(dnECC,ncart_act,nderiv)
      CALL Set_ZERO_TO_dnSVM(dnECC,nderiv)

!----------------------------------------------------------------------
!     ----------------------------------------------------------------
!     - read the energy from the file energy
      file_FChk%name = fchk_name

      CALL file_open(file_FChk,nio)

      CALL Find_Label(nio,'Total Energy',located)
      IF (debug) write(out_unitp,*) 'located: Total Energy',located
      IF (located) THEN
        read(nio,*,iostat=err) name1_i,dnECC%d0
      ELSE
        err = -1
      END IF

      IF (.NOT. located .OR. err /=0) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'I cannot find the energy in :',file_FChk%name
        write(out_unitp,*) 'located,err',located,err
        STOP
      END IF

      close(nio)
!     ----------------------------------------------------------------
!----------------------------------------------------------------------

      CALL file_open(file_FChk,nio)
      IF (nderiv >= 1) THEN

        CALL Find_Label(nio,'Cartesian Gradient',located)
        IF (debug) write(out_unitp,*) 'located: Cartesian Gradient',located
        IF (located) THEN
          read(nio,*,iostat=err)
          read(nio,*,iostat=err) dnECC%d1(:)
        END IF
        IF (.NOT. located .OR. err /=0) THEN
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) 'I cannot find the Gradient in :',file_FChk%name
          write(out_unitp,*) 'located,err',located,err
          STOP
        END IF
      END IF

      IF (nderiv == 2) THEN

        CALL Find_Label(nio,'Cartesian Force Constants',located)
        IF (debug) write(out_unitp,*) 'located: Cartesian Force Constants (hessian)',located
        IF (located) THEN
          read(nio,*,iostat=err)
          read(nio,*,iostat=err) ((dnECC%d2(i,j),i=1,j),j=1,ncart_act)
        END IF
        IF (.NOT. located .OR. err /=0) THEN
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) 'I cannot find the hessian in :',file_FChk%name
          write(out_unitp,*) 'located,err',located,err
          STOP
        END IF
!       CALL Write_Mat(dnECC%d2,out_unitp,5)

        DO j=1,ncart_act
        DO i=1,j-1
          dnECC%d2(j,i) = dnECC%d2(i,j)
        END DO
        END DO

      END IF
      close(nio)

!     ----------------------------------------------------------------
!----------------------------------------------------------------------

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'dnECC'
        CALL Write_dnS(dnECC,nderiv)
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------


      END SUBROUTINE Read_dnECC_Gauss
      SUBROUTINE Read_dnDipCC_Gauss(dnDipCC,fchk_name,nderiv,ncart_act)

      USE mod_system
      USE mod_dnSVM
      implicit none

!----- input output variables ----------------------------------------
      integer, intent(in)      :: nderiv,ncart_act
      character (len=Line_len) :: fchk_name
      TYPE (param_file)        :: file_FChk
      TYPE(Type_dnS)           :: dnDipCC(3)

!----- for the files -------------------------------------------------
      integer :: nio

      integer                   :: i,j,nderiv_loc
      character (len=Name_len)  :: name1_i
      logical                   :: located
      integer                   :: err

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Read_dnDipCC_Gauss'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (nderiv > 1) THEN
         nderiv_loc = 1
      ELSE
         nderiv_loc = nderiv
      END IF
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING',name_sub
        write(out_unitp,*) 'nderiv',nderiv_loc
        write(out_unitp,*) 'fchk_name: ',fchk_name
        write(out_unitp,*)
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------
      IF (nderiv_loc < 0) RETURN

      CALL alloc_VecOFdnS(dnDipCC,ncart_act,nderiv_loc)
      CALL sub_ZERO_TO_VecOFdnS(dnDipCC,nderiv_loc)

!----------------------------------------------------------------------
!     ----------------------------------------------------------------
!     - read the Dipole Moment from the file Test.FChk
      file_FChk%name = fchk_name

      CALL file_open(file_FChk,nio)
      CALL Find_Label(nio,'Dipole Moment',located)
      IF (debug) write(out_unitp,*) 'Located: Dipole Moment ',located
      IF (located) THEN
        read(nio,*,iostat=err)
        read(nio,*,iostat=err) dnDipCC(:)%d0
      ELSE
        err = -1
      END IF

      IF (.NOT. located .OR. err /=0) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'I cannot find the Dipole Moment in :',fchk_name
        write(out_unitp,*) 'located,err',located,err
        STOP
      END IF
      close(nio)

!     - read the gradient of the Dipole Moment
      IF (nderiv_loc >= 1) THEN
        CALL file_open(file_FChk,nio)

        CALL Find_Label(nio,'Dipole Derivatives',located)
        IF (debug) write(out_unitp,*) 'Located: Dipole Derivatives ',located
        IF (located) THEN

          read(nio,*,iostat=err)
          read(nio,*,iostat=err) ((dnDipCC(i)%d1(j),i=1,3),j=1,ncart_act)

        END IF
        IF (.NOT. located .OR. err /=0) THEN
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) 'I cannot find the Dipole Derivatives in :',fchk_name
          write(out_unitp,*) 'located,err',located,err
          STOP
        END IF
        close(nio)
      END IF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'dnDipCC(:)'
        CALL Write_VecOFdnS(dnDipCC)
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------


      END SUBROUTINE Read_dnDipCC_Gauss

     SUBROUTINE Read_dnPolarizabilityCC_Gauss(dnPolarCC,fchk_name,nderiv,ncart_act)

      USE mod_system
      USE mod_dnSVM
      implicit none

!----- input output variables ----------------------------------------
      integer, intent(in)      :: nderiv,ncart_act
      character (len=Line_len) :: fchk_name
      TYPE (param_file)        :: file_FChk
      TYPE(Type_dnS)           :: dnPolarCC(6) ! αxx, αxy, αyy, αxz, αyz, αzz

!----- for the files -------------------------------------------------
      integer :: nio

      integer                   :: i,j,nderiv_loc
      character (len=Name_len)  :: name1_i
      logical                   :: located
      integer                   :: err

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Read_dnPolarizabilityCC_Gauss'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (nderiv > 1) THEN
         nderiv_loc = 1
      ELSE
         nderiv_loc = nderiv
      END IF
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING',name_sub
        write(out_unitp,*) 'nderiv',nderiv_loc
        write(out_unitp,*) 'fchk_name: ',fchk_name
        write(out_unitp,*)
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------
      IF (nderiv_loc < 0) RETURN

      CALL alloc_VecOFdnS(dnPolarCC,ncart_act,nderiv_loc)
      CALL sub_ZERO_TO_VecOFdnS(dnPolarCC,nderiv_loc)

!----------------------------------------------------------------------
!     ----------------------------------------------------------------
!     - read the Dipole Moment from the file Test.FChk
      file_FChk%name = fchk_name

      CALL file_open(file_FChk,nio)
      CALL Find_Label(nio,'Polarizability',located)
      IF (debug) write(out_unitp,*) 'Located: Polarizability ',located
      IF (located) THEN
        read(nio,*,iostat=err)
        read(nio,*,iostat=err) dnPolarCC(:)%d0
      ELSE
        err = -1
      END IF

      IF (.NOT. located .OR. err /=0) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'I cannot find the Polarizability in :',fchk_name
        write(out_unitp,*) 'located,err',located,err
        !STOP
      END IF
      close(nio)

!     - read the gradient of the Dipole Moment
      IF (nderiv_loc >= 1) THEN
        CALL file_open(file_FChk,nio)

        CALL Find_Label(nio,'Polarizability Derivatives',located)
        IF (debug) write(out_unitp,*) 'Located: Polarizability Derivatives ',located
        IF (located) THEN

          read(nio,*,iostat=err)
          read(nio,*,iostat=err) ((dnPolarCC(i)%d1(j),i=1,6),j=1,ncart_act)

        END IF
        IF (.NOT. located .OR. err /=0) THEN
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) 'I cannot find the Polarizability Derivatives in :',fchk_name
          write(out_unitp,*) 'located,err',located,err
          !STOP
        END IF
        close(nio)
      END IF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'dnPolarCC(:)'
        CALL Write_VecOFdnS(dnPolarCC)
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------


      END SUBROUTINE Read_dnPolarizabilityCC_Gauss

      SUBROUTINE pot_mu_onthefly_gamess(Qxyz,dnECC,nderivE,dnMuCC,nderivMu, &
                                         mole,para_PES)

      USE mod_system
      USE mod_dnSVM
      USE mod_Coord_KEO,  only: assignment(=),CoordType
      USE mod_PrimOp_def
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole

!----- for Qdyn Qact ... ---------------------------------------------
      real (kind=Rkind) :: Qxyz(mole%ncart_act)
      TYPE (param_PES)  :: para_PES

!----- input output variables ----------------------------------------
      integer        :: nderivE,nderivMu
      TYPE(Type_dnS) :: dnECC,dnMuCC(3)


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='pot_mu_onthefly_gamess'
      logical, parameter :: debug=.FALSE.
!      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING',name_sub
        write(out_unitp,*) 'Qxyz',Qxyz
        write(out_unitp,*) 'nderivE,nderivMu',nderivE,nderivMu
        write(out_unitp,*)
      END IF
!-----------------------------------------------------------

      IF (.NOT. para_PES%Read_OnTheFly_only) THEN

        CALL Calc_EneDip_WITH_gamess(Qxyz,nderivE,nderivMu,             &
                                        mole,para_PES,para_PES%para_OTF)

      END IF



      !-----------------------------------------------------------------
      !- read the energy from the file energy, gradient, hessian
      CALL Read_dnECC_Gamess(dnECC,para_PES%para_OTF%file_log%name,     &
                 para_PES%para_OTF%file_pun%name,nderivE,mole%ncart_act)

      !-----------------------------------------------------------------
      !- read the Dipole Moment from the file xx.pun
      IF (nderivMu > -1) THEN
        CALL Read_dnDipCC_Gamess(dnMuCC,para_PES%para_OTF%file_pun%name,&
                                 nderivMu,mole%ncart_act)

      END IF
      !-----------------------------------------------------------------


!-----------------------------------------------------------
      IF (debug) THEN
         write(out_unitp,*) 'dnECC'
         CALL Write_dnSVM(dnECC)
         IF (nderivMu > -1) THEN
           write(out_unitp,*) 'dnMuCC(1): DipX'
           CALL Write_dnSVM(dnMuCC(1))
           write(out_unitp,*) 'dnMuCC(2): DipY'
           CALL Write_dnSVM(dnMuCC(2))
           write(out_unitp,*) 'dnMuCC(:): DipZ'
           CALL Write_dnSVM(dnMuCC(3))
         END IF
         write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------


      END SUBROUTINE pot_mu_onthefly_gamess
      SUBROUTINE Calc_EneDip_WITH_gamess(Qxyz,nderivE,nderivDip,&
                                        mole,para_PES,para_OTF)

      USE mod_system
      USE mod_Coord_KEO,  only: assignment(=),CoordType
      USE mod_PrimOp_def
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType)    :: mole

!----- for Qdyn Qact ... ---------------------------------------------
      real (kind=Rkind) :: Qxyz(mole%ncart_act)
      TYPE (param_PES)  :: para_PES
      TYPE (param_OTF)  :: para_OTF

!----- for the files -------------------------------------------------
      integer :: nio,nio_header

!----- input output variables ----------------------------------------
      integer        :: nderivE,nderivDip

!----- working variables ----------------------------------------
      integer       :: err
      integer                       :: i,iq,Z

      character (len=Name_longlen)  :: ab_level
      character (len=Name_longlen)  :: labelR
      logical                       :: located
      character (len=Line_len)      :: line

!----- function -------------------------------------------------
!----------------------------------------------------------------

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Calc_EneDip_WITH_gamess'
      logical, parameter :: debug=.FALSE.
!      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING',name_sub
        write(out_unitp,*) 'Qxyz',Qxyz(:)
        write(out_unitp,*) 'nderivE,nderivDip',nderivE,nderivDip
        write(out_unitp,*)
      END IF
      !----------------------------------------------------------------

      !- gamess input file ------------------------------------------
      ! write(out_unitp,*) 'file_data gamess ',para_OTF%file_data

      CALL file_open(para_OTF%file_data,nio)

      IF (para_OTF%header) THEN
        CALL file_open(para_OTF%file_header,nio_header)

        DO
          read(nio_header,21,end=998) line
 21       format(132A)
          write(nio,*) trim(adjustl(line))
        END DO
 998    close(nio_header)
      END IF

      write(nio,*) '$CONTRL'
      write(nio,*) '  EXETYP=RUN UNITS=BOHR NOSYM=1'
      write(nio,*) '  COORD=UNIQUE '
      write(nio,*) '  ',para_OTF%ab_initio_meth
      write(nio,*) '  ICHARG=',mole%charge
      write(nio,*) '  MULT=',mole%multiplicity
      IF (nderivE == 0) write(nio,*) '  RUNTYP=energy'
      IF (nderivE == 1) write(nio,*) '  RUNTYP=GRADIENT'
      IF (nderivE == 2) write(nio,*) '  RUNTYP=HESSIAN'
      write(nio,*) '$END'
      write(nio,*) '$BASIS'
      write(nio,*) '  ',para_OTF%ab_initio_basis
      write(nio,*) '$END'
      write(nio,*) '$DATA'
      write(nio,*) '  xx'
      write(nio,*) 'C1'

      !write(out_unitp,*) 'Z',mole%Z(:)
      iq=0
      DO i=1,mole%nat_act
        Z = mole%Z(i)
        IF (mole%Z(i) == 0) Z=0
        write(nio,13) 'at ',real(Z,kind=Rkind),Qxyz(iq+1:iq+3)
 13     format(a,1x,f12.6,3(1x,f20.15))
        iq = iq+3
      END DO
      write(nio,*) '$END'

      close(nio)
      !---------------------------------------------------------------

!       - gamess execution -------------------------------------------
        CALL file_delete(para_OTF%file_log)

        !CALL system(para_OTF%commande_unix)
        CALL EXECUTE_COMMAND_LINE(para_OTF%commande_unix)


        located = .FALSE.
        CALL file_open(para_OTF%file_log,nio)
        CALL Find_Label(nio,' ddikick.x: exited gracefully.',located)
        close(nio)

        !located = .FALSE.
        !CALL file_open(para_OTF%file_log,nio,append=.TRUE.)
        !backspace(nio,err=999)
        !read(nio,'(a32)',err=999) labelR
        !located = verify(labelR,'ddikick.x: exited gracefully.') == 0
        !close(nio)

 999    IF (.NOT. located) THEN
          write(out_unitp,*) 'ERROR in the gamess execution'
          write(out_unitp,*) 'no line: "ddikick.x: exited gracefully."'
          write(out_unitp,*) 'last line: ',labelR
          write(out_unitp,*) 'unix command:',para_OTF%commande_unix
          STOP
        END IF


      IF (debug) THEN
        write(out_unitp,*) 'END',name_sub
      END IF

      END SUBROUTINE Calc_EneDip_WITH_gamess
      SUBROUTINE Read_dnECC_Gamess(dnECC,log_name,pun_name,nderiv,ncart_act)

      USE mod_system
      USE mod_dnSVM
      implicit none

!----- input output variables ----------------------------------------
      integer                  :: nderiv,ncart_act
      character (len=Line_len) :: log_name,pun_name
      TYPE (param_file)        :: file_log,file_pun
      TYPE(Type_dnS)           :: dnECC

!----- for the files -------------------------------------------------
      integer :: nio

      integer                   :: i,j,k,iq,idum,jdum,nbligne,nbreste
      real (kind=Rkind)         :: RZ
      character (len=Name_len)  :: name1_i
      logical                   :: located
      integer                   :: err


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Read_dnECC_Gamess'
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING',name_sub
        write(out_unitp,*) 'nderiv',nderiv
        write(out_unitp,*) 'log_name: ',log_name
        write(out_unitp,*)
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

      CALL alloc_dnSVM(dnECC,ncart_act,nderiv)
      CALL Set_ZERO_TO_dnSVM(dnECC,nderiv)


      !----------------------------------------------------------------
      !- read the energy from the file energy
      file_log%name = log_name
      CALL file_open(file_log,nio)

      CALL NFind_Label(nio,'TOTAL ENERGY =',located,37)
      !CALL NFind_Label(nio,'E(MP2)=',located,20)
      !CALL NFind_Label(nio,'FINAL R-PM3 ENERGY IS',located,22)
      IF (located) THEN
        read(nio,*,iostat=err) dnECC%d0
      ELSE
        err = -1
      END IF
!     write(out_unitp,*) 'located,err',located,err

      IF (.NOT. located .OR. err /=0) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'I cannot find the energy in :',log_name
        write(out_unitp,*) 'located,err',located,err
        STOP
      END IF
      close(nio)
      !----------------------------------------------------------------

      !----------------------------------------------------------------
      !- read the gradient and the hessian
      file_pun%name = pun_name
      CALL file_open(file_pun,nio)

      IF (nderiv >= 1) THEN

        CALL NFind_Label(nio,'$GRAD',located,6)
        IF (debug) write(out_unitp,*) 'located gradient',located
        IF (located) THEN
          read(nio,*,iostat=err)
          read(nio,*,iostat=err)
          iq = 0
          DO i=1,ncart_act/3
            read(nio,*,iostat=err) name1_i,RZ,dnECC%d1(iq+1:iq+3)
            iq = iq + 3
          END DO
        END IF
        IF (.NOT. located .OR. err /= 0) THEN
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) 'I cannot find the Gradient in :',pun_name
          write(out_unitp,*) 'located,err',located,err
          STOP
        END IF

      END IF

      IF (nderiv == 2) THEN

        CALL NFind_Label(nio,' $HESS',located,6)
        IF (debug) write(out_unitp,*) 'located hessian',located
        IF (located) THEN
          read(nio,*,iostat=err)
          read(nio,*,iostat=err)
          nbligne = int(ncart_act/5)
          nbreste = ncart_act-5*nbligne
          DO  i=1,ncart_act
            DO  j=0,nbligne-1
              read(nio,41,iostat=err) idum,jdum,(dnECC%d2(i,j*5+k),k=1,5)
 41           format (i2,i3,5E16.9)
              !Recommended by ifort, W>=D+7, was 5E15.9
            END DO
            IF (nbreste > 0 ) THEN
              j = nbligne
              read(nio,41,iostat=err) idum,jdum,(dnECC%d2(i,j*5+k),k=1,nbreste)
            END IF
          END DO
        END IF
        IF (.NOT. located .OR. err /=0) THEN
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) 'I cannot find the hessian in :',pun_name
          write(out_unitp,*) 'located,err',located,err
          STOP
        END IF

      END IF
      close(nio)



      IF (debug) THEN
        write(out_unitp,*) 'dnECC'
        CALL Write_dnS(dnECC,nderiv)
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------


      END SUBROUTINE Read_dnECC_Gamess
      SUBROUTINE Read_dnDipCC_Gamess(dnDipCC,pun_name,nderiv,ncart_act)

      USE mod_system
      USE mod_dnSVM
      implicit none

!----- input output variables ----------------------------------------
      integer                  :: nderiv,ncart_act
      character (len=Line_len) :: pun_name
      TYPE (param_file)        :: file_pun
      TYPE(Type_dnS)           :: dnDipCC(3)

!----- for the files -------------------------------------------------
      integer :: nio
      real(kind=Rkind)         :: conv

      integer                   :: i,j
      character (len=Name_len)  :: name1_i
      logical                   :: located
      integer                   :: err
      real(kind=Rkind)          :: convDebyeTOau,a0

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Read_dnDipCC_Gamess'
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING',name_sub
        write(out_unitp,*) 'nderiv',nderiv
        write(out_unitp,*) 'pun_name: ',pun_name
        write(out_unitp,*)
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------
      IF (nderiv < 0) RETURN

      a0            =     get_Conv_au_TO_unit("L" ,"Angs")
      convDebyeTOau = ONE/get_Conv_au_TO_unit("QL","D") ! Debye

      CALL alloc_VecOFdnS(dnDipCC,ncart_act,nderiv)
      CALL sub_ZERO_TO_VecOFdnS(dnDipCC,nderiv)

      !----------------------------------------------------------------
      !- read the Dipole Moment from the file xx.pun
      file_pun%name = pun_name

      CALL file_open(file_pun,nio)
      CALL NFind_Label(nio,'DIPOLE',located,7)
      IF (located) THEN
        read(nio,*,iostat=err) dnDipCC(:)%d0
      ELSE
        err = -1
      END IF

      IF (.NOT. located .OR. err /=0) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'I cannot find the Dipole Moment in :',pun_name
        write(out_unitp,*) 'located,err',located,err
        STOP
      END IF
      close(nio)

      dnDipCC(:)%d0 = dnDipCC(:)%d0 * convDebyeTOau ! because mu is in Debye


!     - read the gradient of the Dipole Moment
      IF (nderiv >= 1) THEN
        conv = convDebyeTOau * a0 ! gradient in debye / Angstrom
        CALL file_open(file_pun,nio)

        CALL NFind_Label(nio,'$DIPDR',located,7)
        IF (debug) write(out_unitp,*) 'Located Dipole Derivatives ?',located
        IF (located) THEN
          read(nio,*,iostat=err)
          DO i=1,ncart_act
            read(nio,'(1x,3e15.8)',iostat=err) (dnDipCC(j)%d1(i),j=1,3)
          END DO

          DO j=1,3
            dnDipCC(j)%d1(:) = dnDipCC(j)%d1(:) * conv
          END DO

          close(nio)
        END IF
      END IF


!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'dnDipCC(:)'
        CALL Write_VecOFdnS(dnDipCC,nderiv)
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------


      END SUBROUTINE Read_dnDipCC_Gamess

!================================================================
!    subroutine enables to calculate the energy, gradient and hessian
!    directly with gaussian
!
!    input : coordinates Qdyn (used in the dynamics) unit (bohr, radian)
!            nderiv = 0 (pot only : d0pot)
!            nderiv = 1 (pot + gradient only : d0pot, d1pot)
!            nderiv = 2 (pot + gradient + hessian only : d0pot, d1pot, d2pot)
!
!    ouput : d0pot, d1pot, d2pot
!================================================================
      SUBROUTINE Read_GradHess_Molpro(dnFCC,outm_name,nderiv,ncart_act)

      USE mod_system
      USE mod_dnSVM
      implicit none

!----- input output variables ----------------------------------------
      integer                  :: nderiv,ncart_act
      character (len=Line_len) :: outm_name
      TYPE (param_file)        :: file_outm
      TYPE(Type_dnS)           :: dnFCC

!----- for the files -------------------------------------------------
      integer :: nio


      integer                    :: i,j,k,l,vi,vj,iq
      integer                    :: icart,icart_ini,iblock,nb_blocks,idum

      character (len=Name_len)      :: ab_level
      character (len=Name_longlen)  :: labelR
      integer             :: Z
      character (len=Name_len)  :: name1_i,name2_i,name3_i
      integer             :: i1,i2,i3
      logical             :: located
      integer             :: err

      character (len=Line_len) :: line


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Read_hess_Molpro'
      logical, parameter :: debug=.FALSE.
!      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING',name_sub
        write(out_unitp,*) 'nderiv',nderiv
        write(out_unitp,*) 'outm_name: ',outm_name
        write(out_unitp,*)
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

!----------------------------------------------------------------------
!     ----------------------------------------------------------------
!     - read the energy from the file energy
!     Not easy with molpro
!     ----------------------------------------------------------------
!----------------------------------------------------------------------

 31   format(1x,5D15.8)
      CALL alloc_dnSVM(dnFCC,ncart_act,nderiv)
      file_outm%name = outm_name

      IF (nderiv >=1) THEN ! the gradient

        dnFCC%d1(:)   = ZERO

        CALL file_open(file_outm,nio)
        CALL FFind_Label(nio,                                           &
         ' Atom          dE/dx               dE/dy               dE/dz',&
                         located,'(a60)')
        write(out_unitp,*) 'located Gradient',located
        CALL flush_perso(out_unitp)
        IF (located) THEN
          read(nio,*,iostat=err)
          DO icart=1,ncart_act,3
            read(nio,*,iostat=err) idum,dnFCC%d1(icart:icart+2)
          END DO
        ELSE
          err = 0
          write(out_unitp,*) 'We cannot find the gradient!'
          write(out_unitp,*) '=> The gradient is assumed to be zero'
          dnFCC%d1(:) = ZERO
        END IF
        IF (err /= 0) THEN
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) 'Problem while reading the gradient :',file_outm%name
          write(out_unitp,*) 'located,err',located,err
          STOP
        END IF
        IF (debug) THEN
          write(out_unitp,31) dnFCC%d1(:)
        END IF
        close(nio)
        CALL flush_perso(out_unitp)
      END IF

      IF (nderiv == 2) THEN

        dnFCC%d2(:,:) = ZERO

        CALL file_open(file_outm,nio)
        CALL FFind_Label(nio,' Force Constants',located,'(a16)')
        write(out_unitp,*) 'located hessian',located
        CALL flush_perso(out_unitp)
        IF (located) THEN
          nb_blocks = ncart_act/5
          IF (nb_blocks*5 < ncart_act) nb_blocks = nb_blocks + 1
          icart_ini = 1
          DO iblock=1,nb_blocks
            read(nio,*,iostat=err)
            DO icart=icart_ini,ncart_act
              IF (icart < icart_ini+5) THEN
                read(nio,*,iostat=err) labelR,dnFCC%d2(icart,icart_ini:icart)
              ELSE
                read(nio,*,iostat=err) labelR,dnFCC%d2(icart,icart_ini:icart_ini+4)
              END IF
              !write(out_unitp,*) 'iblock,icart,labelR',iblock,icart,labelR
            END DO
            icart_ini = icart_ini + 5
          END DO
        ELSE
          err = -1
        END IF
        IF (.NOT. located .OR. err /=0) THEN
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) 'I cannot find the hessian in :',file_outm%name
          write(out_unitp,*) 'located,err',located,err
          STOP
        END IF
        !CALL Write_Mat(dnFCC%d2,out_unitp,5)
        DO j=1,ncart_act
        DO i=1,j-1
          dnFCC%d2(i,j) = dnFCC%d2(j,i)
        END DO
        END DO
        IF (debug) THEN
          CALL Write_Mat(dnFCC%d2,out_unitp,5)
          CALL flush_perso(out_unitp)
        END IF

        close(nio)
      END IF
!     ----------------------------------------------------------------
!----------------------------------------------------------------------
        CALL flush_perso(out_unitp)
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------


      END SUBROUTINE Read_GradHess_Molpro
      SUBROUTINE onthefly_generic(Qxyz,MatdnECC,nderivE,MatdnMuCC,nderivMu, &
                                  mole,para_PES)

      USE mod_system
      USE mod_dnSVM
      USE mod_Coord_KEO,  only: assignment(=),CoordType
      USE mod_PrimOp_def
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType)    :: mole

!----- for Qdyn Qact ... ---------------------------------------------
      real (kind=Rkind) :: Qxyz(mole%ncart_act)
      TYPE (param_PES)  :: para_PES


!----- input output variables ----------------------------------------
      integer        :: nderivE,nderivMu
      TYPE(Type_dnS) :: MatdnECC(para_PES%nb_elec,para_PES%nb_elec)
      TYPE(Type_dnS) :: MatdnMuCC(para_PES%nb_elec,para_PES%nb_elec,3)

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='onthefly_generic'
      logical, parameter :: debug=.FALSE.
!      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING',name_sub
        write(out_unitp,*) 'Qxyz',Qxyz
        write(out_unitp,*) 'nderivE,nderivMu',nderivE,nderivMu
        write(out_unitp,*) 'Read_OnTheFly_only',para_PES%Read_OnTheFly_only
        write(out_unitp,*)
      END IF
!-----------------------------------------------------------
      IF (nderivE >= 1 .OR. nderivMu >= 1) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'The reading of the gradient and the hessian'
        write(out_unitp,*) 'are not ready with "generic" on-the-fly calculation'
        STOP
      END IF

      !-----------------------------------------------------------------
      IF (.NOT. para_PES%Read_OnTheFly_only) THEN

        CALL Calc_EneDip_WITH_generic(Qxyz,nderivE,nderivMu,            &
                                        mole,para_PES,para_PES%para_OTF)
      END IF
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !- read the energy from the file energy, gradient, hessian
      CALL Read_dnECC_generic(MatdnECC,para_PES%para_OTF%file_log%name, &
                                nderivE,para_PES%nb_elec,mole%ncart_act)
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !- read the Dipole Moment from the file Test.FChk
      IF (nderivMu > -1) THEN

        CALL Read_dnDipCC_generic(MatdnMuCC,para_PES%para_OTF%file_log%name,&
                               nderivMu,para_PES%nb_elec,mole%ncart_act)

      END IF
      !-----------------------------------------------------------------



!-----------------------------------------------------------
      IF (debug) THEN
         write(out_unitp,*) 'MatdnECC'
         CALL Write_MatOFdnS(MatdnECC)
         IF (nderivMu > -1) THEN
           write(out_unitp,*) 'MatdnMuCC(:,:,1): MatDipX'
           CALL Write_MatOFdnS(MatdnMuCC(:,:,1))
           write(out_unitp,*) 'MatdnMuCC(:,:,2): MatDipY'
           CALL Write_MatOFdnS(MatdnMuCC(:,:,2))
           write(out_unitp,*) 'MatdnMuCC(:,:,3): MatDipZ'
           CALL Write_MatOFdnS(MatdnMuCC(:,:,3))
         END IF
         write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE onthefly_generic

      SUBROUTINE Calc_EneDip_WITH_generic(Qxyz,nderivE,nderivDip,       &
                                          mole,para_PES,para_OTF)

      USE mod_system
      USE mod_Coord_KEO,  only: assignment(=),CoordType
      USE mod_PrimOp_def
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType)    :: mole

!----- for Qdyn Qact ... ---------------------------------------------
      real (kind=Rkind) :: Qxyz(mole%ncart_act)
      TYPE (param_PES)  :: para_PES
      TYPE (param_OTF)  :: para_OTF

!----- for the files -------------------------------------------------
      integer :: nio,nio_header

!----- input output variables ----------------------------------------
      integer        :: nderivE,nderivDip

!----- working variables ----------------------------------------
      integer       :: err
      integer                       :: i,iq,Z

      character (len=Name_longlen)  :: ab_level
      character (len=Name_longlen)  :: labelR
      logical                       :: located
      character (len=Line_len)      :: line

!----- function -------------------------------------------------
!----------------------------------------------------------------

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Calc_EneDip_WITH_generic'
      logical, parameter :: debug=.FALSE.
!      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING',name_sub
        write(out_unitp,*) 'Qxyz',Qxyz(:)
        write(out_unitp,*) 'nderivE,nderivDip',nderivE,nderivDip
        write(out_unitp,*)
      END IF
      !----------------------------------------------------------------

      !----------------------------------------------------------------
      CALL file_open(para_PES%para_OTF%file_data,nio)
      write(nio,*) 'geom'

      iq=0
      DO i=1,mole%nat_act
        Z = mole%Z(i)
        IF (mole%Z(i) == 0) Z=0
          !mass = get_mass_Tnum(mole%mendeleev,Z=Z,name=name_at)
          !write(nio,13) trim(name_at),Z,Qxyz(iq+1),Qxyz(iq+2),Qxyz(iq+3)
        write(nio,13) "At ",Z,Qxyz(iq+1),Qxyz(iq+2),Qxyz(iq+3)

 13     format(a,1x,i4,3(1x,f20.15))
        iq = iq+3
      END DO

      write(nio,*) 'end geom'
      close(nio)
      !----------------------------------------------------------------



      !----------------------------------------------------------------
      !- system call => ab initio calulation --------------------------
      !CALL system(para_PES%para_OTF%commande_unix // " " // int_TO_char(para_PES%nb_elec) )
      CALL EXECUTE_COMMAND_LINE(para_PES%para_OTF%commande_unix // " " // int_TO_char(para_PES%nb_elec) )

      located = .FALSE.
      CALL file_open(para_PES%para_OTF%file_log,nio,append=.TRUE.)
      backspace(nio,err=999)
      read(nio,'(a32)',err=999) labelR
      located = verify(labelR,'calculation done') == 0
      close(nio)

 999  IF (.NOT. located) THEN
        write(out_unitp,*) 'ERROR in the generic execution'
        write(out_unitp,*) 'no line: "calculation done"'
        write(out_unitp,*) 'last line: ',labelR
        write(out_unitp,*) 'unix command:',para_PES%para_OTF%commande_unix
        STOP
      END IF

      IF (debug) THEN
        write(out_unitp,*) 'END',name_sub
      END IF

      END SUBROUTINE Calc_EneDip_WITH_generic
      SUBROUTINE Read_dnECC_generic(MatdnECC,log_name,nderiv,nb_elec,ncart_act)

      USE mod_system
      USE mod_dnSVM
      implicit none

!----- input output variables ----------------------------------------
      integer                  :: nderiv,ncart_act,nb_elec
      character (len=Line_len) :: log_name
      TYPE (param_file)        :: file_log
      TYPE(Type_dnS)           :: MatdnECC(nb_elec,nb_elec)

!----- for the files -------------------------------------------------
      integer :: nio

      integer                   :: i1,i2,icount
      logical                   :: located
      integer                   :: err
      real(kind=Rkind)          :: val


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Read_dnECC_generic'
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING',name_sub
        write(out_unitp,*) 'nderiv',nderiv
        write(out_unitp,*) 'log_name: ',log_name
        write(out_unitp,*)
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

      CALL alloc_MatOFdnS(MatdnECC,ncart_act,nderiv)
      CALL sub_ZERO_TO_MatOFdnS(MatdnECC,nderiv)


      !----------------------------------------------------------------
      !- read the energy
      file_log%name = log_name
      CALL file_open(file_log,nio)

      CALL NFind_Label(nio,'energy',located,len('energy'))
      IF (located) THEN
        icount = 0
        DO
          read(nio,*,iostat=err) i1,i2,val
          IF (err /= 0) THEN
            ! end energy or error
            IF (icount /= 0) err = 0
            EXIT
          END IF
          icount = icount + 1
          MatdnECC(i1,i2)%d0 = val
        END DO
      ELSE
        err = -1
      END IF
      !write(out_unitp,*) 'located,err',located,err
      !write(out_unitp,*) 'MatdnECC(:,:)%d0',MatdnECC(:,:)%d0

      IF (.NOT. located .OR. err /=0) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'I cannot find the energy in :',log_name
        write(out_unitp,*) 'located,err',located,err
        STOP
      END IF
      close(nio)
      !----------------------------------------------------------------

      IF (nderiv >= 1) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'The reading of the gradient and the hessian'
        write(out_unitp,*) 'are not ready with "generic" on-the-fly calculation'
        STOP
      END IF


!     ----------------------------------------------------------------
!----------------------------------------------------------------------

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'MatdnECC'
        CALL Write_MatOFdnS(MatdnECC,nderiv)
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------


      END SUBROUTINE Read_dnECC_generic
      SUBROUTINE Read_dnDipCC_generic(MatdnDipCC,log_name,nderiv,nb_elec,ncart_act)

      USE mod_system
      USE mod_dnSVM
      implicit none

!----- input output variables ----------------------------------------
      integer                  :: nderiv,nb_elec,ncart_act
      character (len=Line_len) :: log_name
      TYPE (param_file)        :: file_log
      TYPE(Type_dnS) :: MatdnDipCC(nb_elec,nb_elec,3)

!----- for the files -------------------------------------------------
      integer :: nio

      integer                   :: i,i1,i2,icount
      real(kind=Rkind)          :: val
      logical                   :: located
      integer                   :: err

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Read_dnDipCC_generic'
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING',name_sub
        write(out_unitp,*) 'nderiv',nderiv
        write(out_unitp,*) 'log_name: ',log_name
        write(out_unitp,*)
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------
      IF (nderiv < 0) RETURN

      DO i=1,3
        CALL alloc_MatOFdnS(MatdnDipCC(:,:,i),ncart_act,nderiv)
        CALL Sub_ZERO_TO_MatOFdnS(MatdnDipCC(:,:,i))
      END DO



      !----------------------------------------------------------------
      !- read the Dipole Moment
      file_log%name = log_name
      CALL file_open(file_log,nio)

      CALL NFind_Label(nio,'mux',located,len('mux'))
      IF (located) THEN
        icount = 0
        DO
          read(nio,*,iostat=err) i1,i2,val
          IF (err /= 0) THEN
            ! end mux or error
            IF (icount /= 0) err = 0
            EXIT
          END IF
          icount = icount + 1
          MatdnDipCC(i1,i2,1)%d0 = val
        END DO
      ELSE
        err = -1
      END IF

      IF (.NOT. located .OR. err /=0) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'I cannot find the "mux" in :',log_name
        write(out_unitp,*) 'located,err',located,err
        STOP
      END IF

      CALL NFind_Label(nio,'muy',located,len('muy'))
      IF (located) THEN
        DO
          read(nio,*,iostat=err) i1,i2,val
          IF (err /= 0) THEN
            ! end muy or error
            IF (icount /= 0) err = 0
            EXIT
          END IF
          icount = icount + 1
          MatdnDipCC(i1,i2,2)%d0 = val
        END DO
      ELSE
        err = -1
      END IF

      IF (.NOT. located .OR. err /=0) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'I cannot find the "muy" in :',log_name
        write(out_unitp,*) 'located,err',located,err
        STOP
      END IF

      CALL NFind_Label(nio,'muz',located,len('muz'))
      IF (located) THEN
        DO
          read(nio,*,iostat=err) i1,i2,val
          IF (err /= 0) THEN
            ! end muz or error
            IF (icount /= 0) err = 0
            EXIT
          END IF
          icount = icount + 1
          MatdnDipCC(i1,i2,3)%d0 = val
        END DO
      ELSE
        err = -1
      END IF

      IF (.NOT. located .OR. err /=0) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'I cannot find the "muz" in :',log_name
        write(out_unitp,*) 'located,err',located,err
        STOP
      END IF

      close(nio)
      !----------------------------------------------------------------

      IF (nderiv >= 1) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'The reading of the gradient and the hessian'
        write(out_unitp,*) 'are not ready with "generic" on-the-fly calculation'
        STOP
      END IF

!-----------------------------------------------------------
      IF (debug) THEN
         write(out_unitp,*) 'MatdnDipCC(:,:,1)'
         CALL Write_MatOFdnS(MatdnDipCC(:,:,i))
         write(out_unitp,*) 'MatdnDipCC(:,:,2)'
         CALL Write_MatOFdnS(MatdnDipCC(:,:,i))
         write(out_unitp,*) 'MatdnDipCC(:,:,3)'
         CALL Write_MatOFdnS(MatdnDipCC(:,:,i))
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------


      END SUBROUTINE Read_dnDipCC_generic

     SUBROUTINE Read_OnTheFly_OF_PES(para_PES)
      USE mod_system
      USE mod_PrimOp_def
      IMPLICIT NONE

      TYPE (param_PES)    :: para_PES

!     - for the molecule -------------------------------
      integer :: charge,multiplicity
      logical :: header,footer
      character (len=Name_len)      :: file_name_OTF
      character (len=Line_len)      :: commande_unix
      character (len=Name_longlen)  :: ab_initio_meth,ab_initio_basis
      character (len=Name_longlen)  :: ab_initio_methEne,ab_initio_basisEne
      character (len=Name_longlen)  :: ab_initio_methDip,ab_initio_basisdip

      character (len=Name_len)      :: ab_initio_prog

      NAMELIST /OnTheFly/ charge,multiplicity,                          &
                          header,footer,file_name_OTF,                  &
                          commande_unix,ab_initio_prog,                 &
                          ab_initio_meth,ab_initio_basis,               &
                          ab_initio_methEne,ab_initio_basisEne,         &
                          ab_initio_methDip,ab_initio_basisDip

!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Read_OnTheFly_OF_PES'
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
       END IF

      charge = 0
      multiplicity = -1
      header = .FALSE.
      footer   = .FALSE.
      file_name_OTF      = 'xx'
      ab_initio_meth     = ' hf '
      ab_initio_basis    = ' sto-3g '
      ab_initio_methEne  = ''
      ab_initio_basisEne = ''
      ab_initio_methDip  = ''
      ab_initio_basisDip = ''
      ab_initio_prog     = 'g03'
      commande_unix      = 'g03.run xx >err'


      read(in_unitp,OnTheFly,IOSTAT=err_read)
      IF (err_read < 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' End-of-file or End-of-record'
        write(out_unitp,*) ' The namelist "OnTheFly" is probably absent'
        write(out_unitp,*) ' check your data!'
        write(out_unitp,*) ' ERROR in ',name_sub
        STOP
      ELSE IF (err_read > 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' Some parameter names of the namelist "OnTheFly" are probaly wrong'
        write(out_unitp,*) ' check your data!'
        write(out_unitp,OnTheFly)
        write(out_unitp,*) ' ERROR in ',name_sub
        STOP
      END IF
      IF (debug) write(out_unitp,OnTheFly)

      para_PES%para_OTF%charge       = charge
      para_PES%para_OTF%multiplicity = multiplicity

        SELECT CASE (ab_initio_prog)
        CASE ('g03','g98','G98','G03','gaussian')
          para_PES%para_OTF%file_data%name=trim(file_name_OTF)   // '.com'
          para_PES%para_OTF%file_log%name=trim(file_name_OTF)    // '.log'
          para_PES%para_OTF%file_header%name=trim(file_name_OTF) //'.header'
          para_PES%para_OTF%file_footer%name=trim(file_name_OTF)   //'.footer'
          para_PES%para_OTF%file_FChk%name='Test.FChk'
        CASE ('gamess')
          para_PES%para_OTF%file_data%name=trim(file_name_OTF) // '.inpg'
          para_PES%para_OTF%file_log%name=trim(file_name_OTF) // '.outg'
          para_PES%para_OTF%file_header%name=trim(file_name_OTF)//'.header'
          para_PES%para_OTF%file_footer%name=trim(file_name_OTF)//'.footer'
          para_PES%para_OTF%file_pun%name=trim(file_name_OTF) // '.dat'
        CASE ('gamess2014')
          para_PES%para_OTF%file_data%name=trim(file_name_OTF) // '.inp'
          para_PES%para_OTF%file_log%name=trim(file_name_OTF) // '.log'
          para_PES%para_OTF%file_header%name=trim(file_name_OTF)//'.header'
          para_PES%para_OTF%file_footer%name=trim(file_name_OTF)//'.footer'
          para_PES%para_OTF%file_pun%name=trim(file_name_OTF) // '.dat'
        CASE ('generic')
          para_PES%para_OTF%file_data%name=trim(file_name_OTF) // '.evrti'
          para_PES%para_OTF%file_log%name=trim(file_name_OTF) // '.evrto'
          para_PES%para_OTF%file_header%name=trim(file_name_OTF)//'.header'
          para_PES%para_OTF%file_footer%name=trim(file_name_OTF)//'.footer'
          para_PES%para_OTF%file_pun%name=trim(file_name_OTF) // '.dat'
        CASE ('g09','G09') ! particular case because the keyword formchk is obsolet
          para_PES%para_OTF%file_data%name=trim(file_name_OTF)   // '.com'
          para_PES%para_OTF%file_log%name=trim(file_name_OTF)    // '.log'
          para_PES%para_OTF%file_header%name=trim(file_name_OTF) //'.header'
          para_PES%para_OTF%file_footer%name=trim(file_name_OTF)   //'.footer'
          para_PES%para_OTF%file_FChk%name=trim(file_name_OTF)   //'.fchk'
        CASE default ! ERROR: wrong program !
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' The ab initio program is UNKNOWN ',ab_initio_prog
          STOP
        END SELECT
        write(out_unitp,*) ' Files for the OTF'
        write(out_unitp,*) trim(para_PES%para_OTF%file_data%name)
        write(out_unitp,*) trim(para_PES%para_OTF%file_log%name)
        write(out_unitp,*) trim(para_PES%para_OTF%file_Fchk%name)
        write(out_unitp,*) trim(para_PES%para_OTF%file_pun%name)
        write(out_unitp,*) ' Program for the OTF ',ab_initio_prog
        write(out_unitp,*) ' Unix script for the OTF ',commande_unix
        para_PES%para_OTF%header          = header
        para_PES%para_OTF%footer            = footer
        para_PES%para_OTF%file_name       = file_name_OTF
        para_PES%para_OTF%ab_initio_prog  = ab_initio_prog
        para_PES%para_OTF%commande_unix   = commande_unix

        para_PES%para_OTF_DIP = para_PES%para_OTF


        IF (len_trim(ab_initio_methEne) == 0) THEN
          ab_initio_methEne = ab_initio_meth
        END IF
        IF (len_trim(ab_initio_basisEne) == 0) THEN
          ab_initio_basisEne = ab_initio_basis
        END IF

        IF (len_trim(ab_initio_methDip) == 0) THEN
          ab_initio_methDip = ab_initio_meth
        END IF
        IF (len_trim(ab_initio_basisDip) == 0) THEN
          ab_initio_basisDip = ab_initio_basis
        END IF
        para_PES%para_OTF%ab_initio_meth      = ab_initio_methEne
        para_PES%para_OTF%ab_initio_basis     = ab_initio_basisEne
        para_PES%para_OTF_Dip%ab_initio_meth  = ab_initio_methDip
        para_PES%para_OTF_Dip%ab_initio_basis = ab_initio_basisDip

        para_PES%levelEne_EQ_levelDip    =                              &
                        (ab_initio_methEne .EQ. ab_initio_methDip) .AND.&
                        (ab_initio_basisEne .EQ. ab_initio_basisDip)


       IF (debug) THEN
         write(out_unitp,*) 'END ',name_sub
       END IF

      END SUBROUTINE Read_OnTheFly_OF_PES


      SUBROUTINE Find_Label(nio,label,located)
      USE mod_system
      IMPLICIT NONE
      character (len=*) :: label
      logical :: located
      integer :: nio

      character (len=len(label)) :: labelR
      character (len=name_len) :: format_label

      integer, save :: i_line = 0

      located = .FALSE.
      format_label='(A' // int_TO_char(len(label)) // ')'

      !write(out_unitp,*) 'format_label',format_label
!     write(out_unitp,*) 'label to find:',label,located
      DO
        i_line = i_line + 1

        read(nio,format_label,end=999,eor=888,err=888,advance='no') labelR
 11     format(A40)


        located = (labelR .EQ. label)
        !located = verify(label,labelR) == 0
        !write(out_unitp,*) i_line,located,labelR
        IF (located) EXIT
        read(nio,*,end=999,err=888)

 888    CONTINUE
      END DO


 999  CONTINUE
      i_line = 0

      !write(out_unitp,*) 'Find_Label: ',label,located
      END SUBROUTINE Find_Label

      SUBROUTINE NFind_Label(nio,label,located,iformat)
      USE mod_system
      IMPLICIT NONE
      character (len=*) :: label
      integer :: iformat
      character (len=Name_len) :: fformat
      logical :: located
      integer :: nio

      character (len=iformat) :: labelR
      integer, save :: i_line = 0


      located = .FALSE.

      IF (iformat < 1) THEN
        write(out_unitp,*) ' ERROR in NFind_Label'
        write(out_unitp,*) ' iformat < 1',iformat
        STOP
      END IF
      fformat='(A' // int_TO_char(iformat) // ')'
      !write(out_unitp,*) iformat,fformat
      !write(out_unitp,*) 'label to find:',label,located
      DO
        i_line = i_line + 1
        read(nio,fformat,end=999,eor=888,err=888,advance='no') labelR

        located = verify(label,labelR) == 0
        !write(out_unitp,*) i_line,located,labelR
        IF (located) EXIT
        read(nio,*,end=999,err=888)

 888    CONTINUE
      END DO


 999  CONTINUE
      i_line = 0
      end subroutine NFind_Label

      SUBROUTINE FFind_Label(nio,label,located,fformat)
      USE mod_system
      IMPLICIT NONE
      character (len=*)   :: label
      character (len=*)   :: fformat
      logical             :: located
      integer             :: nio

      character (len=256) :: labelR
      integer, save       :: i_line = 0


      located = .FALSE.

      !write(out_unitp,*) 'label to find:',label,located
      DO
        i_line = i_line + 1
        read(nio,fformat,end=999,eor=888,err=888,advance='no') labelR


        located = verify(label,labelR) == 0
        !write(out_unitp,*) i_line,located,labelR
        IF (located) EXIT
        read(nio,*,end=999,err=888)

 888    CONTINUE
      END DO


 999  CONTINUE
      i_line = 0
      IF (located) read(nio,*)

      END SUBROUTINE FFind_Label

END MODULE mod_OTF
