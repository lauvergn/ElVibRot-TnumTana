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
      PROGRAM Tnum_f90
      use mod_system
      use mod_dnSVM
      use mod_Constant
      ! in the use mod_Coord_KEO, we have to use "only", because "calc_freq" is
      !   a subroutine in mod_Coord_KEO and also a variable in the namelist.
      use mod_Coord_KEO,  ONLY: CoordType,Tnum,Read_CoordType,              &
                                read_RefGeom,get_Qact0,sub_QactTOdnx,       &
                                Write_Cartg98,Write_dnx,calc3_f2_f1Q_num,   &
                                get_dng_dnGG,sub_QplusDQ_TO_Cart,           &
                                sub_dnFCC_TO_dnFcurvi,dealloc_CoordType
      use mod_PrimOp

      IMPLICIT NONE

!     - parameters for para_Tnum -----------------------
      TYPE (constant)  :: const_phys
      TYPE (CoordType) :: mole
      TYPE (Tnum)      :: para_Tnum
      TYPE (param_PES) :: para_PES

      real (kind=Rkind) :: vep,rho
      real (kind=Rkind), pointer :: Tdef2(:,:)=>null()
      real (kind=Rkind), pointer :: Tdef1(:)=>null()
      real (kind=Rkind), pointer :: Tcor2(:,:)=>null()
      real (kind=Rkind), pointer :: Tcor1(:)=>null()
      real (kind=Rkind), pointer :: Trot(:,:)=>null()

      TYPE(Type_dnMat) :: dng,dnGG
      TYPE(Type_dnVec) :: dnx

      TYPE(Type_dnS), pointer :: MatdnE(:,:)=>null()
      TYPE(Type_dnS), pointer :: MatdnImE(:,:)=>null()
      TYPE(Type_dnS), pointer :: MatdnScalOp(:,:,:)=>null()
      TYPE (param_dnMatOp), allocatable :: Tab_dnMatOp(:)


      TYPE(Type_dnS)   :: dnFCC,dnFcurvi
      TYPE(Type_dnS)   :: dnMuCC(3),dnMucurvi(3)
      TYPE(Type_dnS)   :: dnPolarCC(6),dnPolarcurvi(6)

      character (len=Line_len) :: outm_name,fchk_name
      real (kind=Rkind), pointer :: d0c_inv(:,:)=>null()
      real (kind=Rkind), pointer :: d0c_ini(:,:)=>null()
      real (kind=Rkind), pointer :: d0k(:,:)=>null()
      real (kind=Rkind), pointer :: d0c(:,:)=>null()
      real (kind=Rkind), pointer :: d0eh(:)=>null()
      real (kind=Rkind), pointer :: freq(:)=>null()
      real (kind=Rkind), allocatable :: grad(:),hess(:,:)

      real (kind=Rkind)  :: norme
      character (len=50) :: wqi,mqi


      integer :: nderiv,i,j,i1,i2,icart,idum
      character (len=Name_longlen) :: name_i,name_j
!     ------------------------------------------------------

!     - for the coordinate values ----------------------------------
      real (kind=Rkind), allocatable :: Qact(:)

!     - working parameters ------------------------------------------
      integer :: err_mem,memory,err_read

      logical :: calc_QTOx,calc_Tnum,calc_gG
      logical :: calc_grad,calc_hessian
      logical :: OnTheFly,calc_freq
      integer :: nderivGg,n_eval
      character (len=*), parameter :: name_sub='Tnum_f90'


      NAMELIST /calculation/ calc_QTOx,calc_Tnum,calc_gG,nderivGg,      &
                             calc_freq,OnTheFly,n_eval,                 &
                             calc_grad,calc_hessian,outm_name,fchk_name

!=======================================================================
!=======================================================================
      CALL versionEVRT(.TRUE.)
      print_level=2

      !CALL sub_constantes(const_phys,Read_Namelist=.FALSE.)

      !-----------------------------------------------------------------
      !     - read the coordinate transformations :
      !     -   zmatrix, polysperical, bunch...
      !     ------------------------------------------------------------
      CALL Read_CoordType(mole,para_Tnum,const_phys)
      !     ------------------------------------------------------------
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !     - read coordinate values -----------------------------------
      !     ------------------------------------------------------------
      CALL read_RefGeom(mole,para_Tnum)
      !     ------------------------------------------------------------
      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      !     ---- TO finalize the coordinates (NM) and the KEO ----------
      !     ------------------------------------------------------------
      CALL Finalize_TnumTana_Coord_PrimOp(para_Tnum,mole,para_PES)
      !-----------------------------------------------------------------
!=======================================================================
!=======================================================================

!===========================================================
!===========================================================

      write(out_unitp,*) "======================================"
      calc_QTOx    = .TRUE.
      calc_Tnum    = .TRUE.
      calc_gG      = .FALSE.
      nderivGg     = 2
      calc_grad    = .FALSE.
      calc_hessian = .FALSE.
      calc_freq    = .FALSE.
      OnTheFly     = .FALSE.
      outm_name    = ''
      fchk_name    = ''
      n_eval       = 1
      read(in_unitp,calculation,IOSTAT=err_read)
      write(out_unitp,calculation)
      IF (err_read /= 0) THEN
        write(out_unitp,*) ' NO namelist "calculation"!'
        write(out_unitp,*) ' => calc_QTOx=t and calc_Tnum=t'

      END IF
      write(out_unitp,*) "======================================"

!===========================================================
!===========================================================

      CALL alloc_NParray(Qact,(/ mole%nb_var /),'Qact',name_sub)
      CALL get_Qact0(Qact,mole%ActiveTransfo)



!-------------------------------------------------
!     - On The Fly calculation -------------------
!     --------------------------------------------
      IF (OnTheFly) THEN
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======== OnTheFly ===================="
        write(out_unitp,*) "======================================"

        nderiv = 2

        allocate(Tab_dnMatOp(para_PES%nb_scalar_Op+2))
        CALL Init_Tab_OF_dnMatOp(Tab_dnMatOp,mole%nb_act,para_PES%nb_elec, &
                                 nderiv,cplx=para_PES%pot_cplx,JRot=para_Tnum%JJ) ! H



        CALL get_dnMatOp_AT_Qact(Qact,Tab_dnMatOp,mole,para_Tnum,para_PES)

        write(out_unitp,*) "Energy: ",Get_Scal_FROM_Tab_OF_dnMatOp(Tab_dnMatOp,1)
        write(out_unitp,*) "Dipole Moments: ",Get_Scal_FROM_Tab_OF_dnMatOp(Tab_dnMatOp,3),&
         Get_Scal_FROM_Tab_OF_dnMatOp(Tab_dnMatOp,4),Get_Scal_FROM_Tab_OF_dnMatOp(Tab_dnMatOp,5)

        IF (nderiv > 0) THEN
          allocate(Grad(mole%nb_act))
          CALL Get_Grad_FROM_Tab_OF_dnMatOp(Grad,Tab_dnMatOp,1)
          write(out_unitp,*) "Grad of E: ",Grad
          CALL Get_Grad_FROM_Tab_OF_dnMatOp(Grad,Tab_dnMatOp,3)
          write(out_unitp,*) "Grad of Dipx: ",Grad
          CALL Get_Grad_FROM_Tab_OF_dnMatOp(Grad,Tab_dnMatOp,4)
          write(out_unitp,*) "Grad of Dipy: ",Grad
          CALL Get_Grad_FROM_Tab_OF_dnMatOp(Grad,Tab_dnMatOp,5)
          write(out_unitp,*) "Grad of Dipz: ",Grad
          deallocate(Grad)
        END IF
        IF (nderiv > 1) THEN
          allocate(hess(mole%nb_act,mole%nb_act))
          CALL Get_Hess_FROM_Tab_OF_dnMatOp(hess,Tab_dnMatOp,1)
          write(out_unitp,*) "Hessian of E: "
          CALL Write_Mat(hess,out_unitp,5)
          deallocate(hess)
        END IF

        !write(out_unitp,*) "======================================"
        !CALL Write_Tab_OF_dnMatOp(Tab_dnMatOp)
        !write(out_unitp,*) "======================================"


        CALL dealloc_Tab_OF_dnMatOp(Tab_dnMatOp)
        deallocate(Tab_dnMatOp)


        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
      END IF
!-------------------------------------------------
!-------------------------------------------------


!-------------------------------------------------
!     - Cartesian coordinates --------------------
!     --------------------------------------------
      IF (calc_QTOx) THEN
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
        CALL time_perso('sub_QactTOdnx')

        nderiv = 0
        CALL alloc_dnSVM(dnx,mole%ncart,mole%nb_act,nderiv)
        write(out_unitp,*) "======================================"
        DO i=1,n_eval-1
          CALL sub_QactTOdnx(Qact,dnx,mole,nderiv,.FALSE.)
        END DO
        write(out_unitp,*) 'dnx: ',mole%ncart
        mole%WriteCC = .TRUE.
        CALL sub_QactTOdnx(Qact,dnx,mole,nderiv,.FALSE.)
        mole%WriteCC = .FALSE.

        CALL write_dnx(1,mole%ncart,dnx,nderiv)

        CALL Write_Cartg98(dnx%d0,mole)

        CALL dealloc_dnSVM(dnx)
        CALL time_perso('sub_QactTOdnx')
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
      END IF
!-------------------------------------------------
!-------------------------------------------------

!-------------------------------------------------
!-------------------------------------------------
!     - calculation of f2, f1, vep, rho ----------
!     --------------------------------------------
      IF (calc_Tnum) THEN
        CALL alloc_array(Tdef2,(/ mole%nb_act,mole%nb_act /),'Tdef2',name_sub)
        CALL alloc_array(Tdef1,(/ mole%nb_act /),            'Tdef1',name_sub)
        CALL alloc_array(Tcor2,(/ mole%nb_act,3 /),          'Tcor2',name_sub)
        CALL alloc_array(Tcor1,(/ 3 /),                      'Tcor1',name_sub)
        CALL alloc_array(Trot, (/ 3,3 /),                    'Trot', name_sub)

        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "====== calc3_f2_f1Q_num =============="
        write(out_unitp,*) "======================================"
        CALL time_perso('calc3_f2_f1Q_num')

        DO i=1,n_eval
          para_Tnum%WriteT = (i == 1)  ! write only when i=1
          CALL  calc3_f2_f1Q_num(Qact,                                  &
                                 Tdef2,Tdef1,vep,rho,                   &
                                 Tcor2,Tcor1,Trot,                      &
                                 para_Tnum,mole)
        END DO

        CALL time_perso('calc3_f2_f1Q_num')
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"

        CALL dealloc_array(Tdef2,'Tdef2',name_sub)
        CALL dealloc_array(Tdef1,'Tdef1',name_sub)
        CALL dealloc_array(Tcor2,'Tcor2',name_sub)
        CALL dealloc_array(Tcor1,'Tcor1',name_sub)
        CALL dealloc_array(Trot, 'Trot', name_sub)
      END IF

!-------------------------------------------------
!-------------------------------------------------
!     FOR G and g metric tensors
      IF (calc_gG) THEN
        CALL alloc_dnSVM(dng ,mole%ndimG,mole%ndimG,mole%nb_act,nderivGg)
        CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,nderivGg)

        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "====== get_dng_dnGG =================="
        write(out_unitp,*) "======================================"
        CALL time_perso('get_dng_dnGG')

        DO i=1,n_eval
          para_Tnum%WriteT    = (i == 1) ! write only when i=1
          CALL get_dng_dnGG(Qact,para_Tnum,mole,dng,dnGG,nderiv=nderivGg)
        END DO
        write(out_unitp,*) ' dng'
        CALL Write_dnSVM(dng,0)
        write(out_unitp,*) ' dnG'
        CALL Write_dnSVM(dnGG,0)

        CALL time_perso('get_dng_dnGG')
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"

        CALL dealloc_dnSVM(dng)
        CALL dealloc_dnSVM(dnGG)
      END IF
!-------------------------------------------------
!-------------------------------------------------

!-------------------------------------------------
!-------------------------------------------------
      !!! frequencies
      IF (calc_freq .AND. .NOT. calc_hessian) THEN
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "====== sub_freq_AT_Qact =============="
        write(out_unitp,*) "======================================"
        CALL alloc_array(freq,(/ mole%nb_act /),"freq",name_sub)


        CALL sub_freq_AT_Qact(freq,Qact,para_Tnum,mole,para_PES,print_freq=.TRUE.)

        write(out_unitp,*) 'ZPE (cm-1): ',HALF*sum(freq(:))*get_Conv_au_TO_unit('E','cm-1')
        write(out_unitp,*) 'ZPE   (eV): ',HALF*sum(freq(:))*get_Conv_au_TO_unit('E','eV')
        write(out_unitp,*) 'ZPE   (au): ',HALF*sum(freq(:))

        DO i=1,mole%nb_act,3
          i2 = min(i+2,mole%nb_act)
          write(out_unitp,'("frequencies (cm-1): ",i0,"-",i0,3(1x,f0.4))') &
                         i,i2,freq(i:i2)* get_Conv_au_TO_unit('E','cm-1')
        END DO

        CALL dealloc_array(freq,"freq",name_sub)

        CALL sub_QplusDQ_TO_Cart(Qact,mole)

        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
      END IF
!-------------------------------------------------
!-------------------------------------------------
      !!! hessian gradient tranformation from cartessian to curvilinear
      IF (calc_hessian .OR. calc_grad) THEN
        nderiv = 1
        IF (calc_hessian) nderiv = 2

        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======= grad/hessian ================="

        para_Tnum%WriteT    = .TRUE.

        IF (nderiv == 2) THEN
          IF (len_trim(outm_name) > 0) THEN
            write(out_unitp,*) 'read FCC from molpro file:',outm_name
            CALL Read_GradHess_Molpro(dnFCC,outm_name,nderiv,mole%ncart_act)
          ELSE IF (len_trim(fchk_name) > 0) THEN
            write(out_unitp,*) 'read FCC from gaussian file:',fchk_name
            CALL Read_hess_Fchk(dnFCC,fchk_name,nderiv,mole%ncart_act)
            CALL Read_dnDipCC_Gauss(dnMuCC,fchk_name,nderiv,mole%ncart_act)
            CALL Read_dnPolarizabilityCC_Gauss(dnPolarCC,fchk_name,nderiv,mole%ncart_act)
          ELSE
            write(out_unitp,*) ' ERROR it is not possible to read ...'
            write(out_unitp,*) ' ... the hessian and gradient from the input file'
            STOP
          END IF
        ELSE ! nderiv=1
          write(out_unitp,*) 'read FCC from the input file:'
          CALL flush_perso(out_unitp)
          CALL alloc_dnSVM(dnFCC,mole%ncart_act,nderiv)

          ! read the gradient
          read(in_unitp,*,IOSTAT=err_read)
          DO icart=1,mole%ncart_act,3
            read(in_unitp,*,iostat=err_read) name_i,dnFCC%d1(icart:icart+2)
          END DO
          !DO icart=1,mole%ncart_act
          !  read(in_unitp,*,iostat=err_read) name_i,dnFCC%d1(icart)
          !END DO

          IF (err_read /= 0) THEN
            write(out_unitp,*) ' ERROR while reading the gradient'
            write(out_unitp,*) ' => check your data!!'
            STOP
          END IF
        END IF

        CALL sub_dnFCC_TO_dnFcurvi(Qact,dnFCC,dnFcurvi,mole)
        write(out_unitp,*) 'Energy=',dnFcurvi%d0
        write(out_unitp,*) 'Gradient in cuvilinear coordinates'
        DO i=1,mole%nb_act
          write(out_unitp,"(i4,1x,a,2f10.6)") i,                        &
                       mole%tab_Qtransfo(mole%nb_Qtransfo)%name_Qin(i), &
                       Qact(i),dnFcurvi%d1(i)
        END DO
        IF (nderiv == 2) THEN
          write(out_unitp,*) 'Curvilinear hessian:'
          CALL Write_VecMat(dnFcurvi%d2,out_unitp,5)
        END IF

        IF (nderiv == 2) THEN
          DO i=1,size(dnMuCC)
            CALL sub_dnFCC_TO_dnFcurvi(Qact,dnMuCC(i),dnMucurvi(i),mole)
          END DO

          write(out_unitp,*) 'Dipole moment:',dnMuCC(:)%d0
          write(out_unitp,*) 'Gradient of the Dipole moment (curvi):'
          DO i=1,mole%nb_act
            write(out_unitp,*) i,Qact(i),(dnMucurvi(j)%d1(i),j=1,size(dnMuCC))
          END DO
        END IF

        IF (nderiv == 2) THEN
          DO i=1,size(dnPolarCC)
            CALL sub_dnFCC_TO_dnFcurvi(Qact,dnPolarCC(i),dnPolarcurvi(i),mole)
          END DO

          write(out_unitp,*) 'Polarizability:',dnPolarCC(:)%d0
          write(out_unitp,*) 'Gradient of the Polarizability (curvi):'
          DO i=1,mole%nb_act
            write(out_unitp,*) i,Qact(i),(dnPolarCC(j)%d1(i),j=1,size(dnPolarCC))
          END DO
        END IF


        IF (calc_freq) THEN
          CALL alloc_array(freq,(/ mole%nb_act /),"freq",name_sub)

          CALL sub_freq_AT_Qact(freq,Qact,para_Tnum,mole,para_PES,d0h_opt=dnFcurvi%d2)

          write(out_unitp,*) 'ZPE (cm-1): ',HALF*sum(freq(:))*get_Conv_au_TO_unit('E','cm-1')
          write(out_unitp,*) 'ZPE   (eV): ',HALF*sum(freq(:))*get_Conv_au_TO_unit('E','eV')
          write(out_unitp,*) 'ZPE   (au): ',HALF*sum(freq(:))

          DO i=1,mole%nb_act,3
            i2 = min(i+2,mole%nb_act)
            write(out_unitp,'("frequencies (cm-1): ",i0,"-",i0,3(1x,f0.4))') &
                         i,i2,freq(i:i2)* get_Conv_au_TO_unit('E','cm-1')
          END DO

          CALL dealloc_array(freq,"freq",name_sub)
        END IF


        write(out_unitp,*) '======================================================'
        write(out_unitp,*) '======================================================'
        write(out_unitp,*) '======================================================'



        CALL dealloc_dnSVM(dnFcurvi)
        CALL dealloc_dnSVM(dnFCC)
      END IF
!-------------------------------------------------
!-------------------------------------------------

      CALL dealloc_CoordType(mole)
      CALL dealloc_NParray(Qact,'Qact',name_sub)


      write(out_unitp,*) 'END Tnum'

      END PROGRAM Tnum_f90
