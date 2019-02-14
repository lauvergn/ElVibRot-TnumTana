      PROGRAM Tnum_f90
      use mod_system
      use mod_dnSVM
      use mod_Constant
      use mod_Coord_KEO

      IMPLICIT NONE

!     - parameters for para_Tnum -----------------------
      TYPE (constant)  :: const_phys
      TYPE (zmatrix)   :: mole
      TYPE (Tnum)      :: para_Tnum

      real (kind=Rkind) :: vep,rho
      real (kind=Rkind), pointer :: Tdef2(:,:)=>null()
      real (kind=Rkind), pointer :: Tdef1(:)=>null()
      real (kind=Rkind), pointer :: Tcor2(:,:)=>null()
      real (kind=Rkind), pointer :: Tcor1(:)=>null()
      real (kind=Rkind), pointer :: Trot(:,:)=>null()

      TYPE(Type_dnMat) :: dng,dnGG
      TYPE(Type_dnVec) :: dnx

!     ------------------------------------------------------

!     - for the coordinate values ----------------------------------
      real (kind=Rkind), allocatable :: Qact(:)

!     - working parameters ------------------------------------------
      integer :: nderiv,err_mem,memory,err_read

      character (len=*), parameter :: name_sub='TEST_TnumTana'


!=======================================================================
!=======================================================================
      CALL versionEVRT(.TRUE.)
      print_level=2

      !-----------------------------------------------------------------
      !     - read the coordinate transformations :
      !     -   zmatrix, polysperical, bunch...
      !     ------------------------------------------------------------
      CALL Read_mole(mole,para_Tnum,const_phys)
      !     ------------------------------------------------------------
      !-----------------------------------------------------------------

      IF (associated(mole%NMTransfo) .OR. associated(mole%RPHTransfo)) THEN
        write(out_unitp,*) "ERROR: This test program cannot be used with"
        write(out_unitp,*) "Normal modes (NM) or RPH"
        STOP
      END IF

      !-----------------------------------------------------------------
      !     - read coordinate values -----------------------------------
      !     ------------------------------------------------------------
      CALL read_RefGeom(mole,para_Tnum)
      !     ------------------------------------------------------------
      !-----------------------------------------------------------------
!=======================================================================
!=======================================================================

!===========================================================
!===========================================================

!===========================================================
!===========================================================

      CALL alloc_NParray(Qact,(/ mole%nb_var /),'Qact',name_sub)
      CALL get_Qact0(Qact,mole%ActiveTransfo)


!-------------------------------------------------
!     - Cartesian coordinates --------------------
!     --------------------------------------------
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
        CALL time_perso('sub_QactTOdnx')

        nderiv = 0
        CALL alloc_dnSVM(dnx,mole%ncart,mole%nb_act,nderiv)
        write(out_unitp,*) "======================================"
        CALL sub_QactTOdnx(Qact,dnx,mole,nderiv,.FALSE.)
        write(out_unitp,*) 'dnx: ',mole%ncart
        CALL write_dnx(1,mole%ncart,dnx,nderiv)

        CALL Write_Cartg98(dnx%d0,mole)

        CALL dealloc_dnSVM(dnx)
        CALL time_perso('sub_QactTOdnx')
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
!-------------------------------------------------
!-------------------------------------------------


      CALL dealloc_zmat(mole)
      CALL dealloc_NParray(Qact,'Qact',name_sub)


      write(out_unitp,*) 'END ',name_sub

      END PROGRAM Tnum_f90
