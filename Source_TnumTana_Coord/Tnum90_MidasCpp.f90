      PROGRAM Tnum90_MidasCpp
      USE mod_system
      USE mod_Tnum
      USE mod_Tana_keo
      USE mod_Tana_Tnum
      USE mod_dnGG_dng
      USE mod_PrimOp_def
      USE mod_OTF
      USE mod_PrimOp
      USE mod_Lib_QTransfo, only : write_dnx
      IMPLICIT NONE


!     - parameters for para_Tnum -----------------------
      TYPE (constant)  :: const_phys
      TYPE (zmatrix)   :: mole
      TYPE (Tnum)      :: para_Tnum
      TYPE (param_PES) :: para_PES

      TYPE(Type_dnMat) :: dnGG

      TYPE(sum_opnd)   :: TWOxKEO,ExpandTWOxKEO


      real (kind=Rkind), allocatable :: Qact(:)
      real (kind=Rkind), allocatable :: Qxyz(:)
!     - working parameters ------------------------------------------
      integer :: nada,i,j,n,ndim
!     ------------------------------------------------------


      NAMELIST /NewQ/ nada


!     - working parameters ------------------------------------------
      integer :: err_mem,memory,err_io
      character (len=*), parameter :: name_sub='Tnum90_MidasCpp'

!===========================================================
!===========================================================
      !para_mem%mem_debug = .TRUE.
      CALL versionEVRT(.TRUE.)

      !-----------------------------------------------------------------
      !     - read the coordinate tansformations :
      !     -   zmatrix, polysperical, bunch...
      !     ------------------------------------------------------------
      CALL Read_mole(mole,para_Tnum,const_phys)
      para_Tnum%MidasCppForm = .TRUE.
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
      para_Tnum%Tana =.FALSE.
      CALL Finalyze_TnumTana_Coord_PrimOp(para_Tnum,mole,para_PES)
      !-----------------------------------------------------------------
!===========================================================
!===========================================================

      CALL alloc_NParray(Qact,(/ mole%nb_var /),'Qact',name_sub)
      CALL get_Qact0(Qact,mole%ActiveTransfo) ! important when constraints (rigid, flexible are added)

!-------------------------------------------------
!     - Cartesian coordinates of the reference geometry
!     --------------------------------------------

       CALL alloc_NParray(Qxyz,(/ mole%ncart /),'Qxyz',name_sub)


       CALL sub_QactTOd0x(Qxyz,Qact,mole,Gcenter=.FALSE.)

       !write(out_unitp,*) 'Qxyz: ',Qxyz
       CALL Write_XYZ(Qxyz,mole)
!-------------------------------------------------
!-------------------------------------------------


!-------------------------------------------------
!  Evaluation of Qact TO xyz (Once)
!-------------------------------------------------
         read(in_unitp,NewQ,IOSTAT=err_io)
         IF (err_io == 0) THEN
           read(in_unitp,*,IOSTAT=err_io) Qact
           IF (err_io == 0) THEN
             CALL sub_QactTOd0x(Qxyz,Qact,mole,Gcenter=.FALSE.)
             CALL Write_XYZ(Qxyz,mole,unit='bohr',io_unit=out_unitp)
           END IF
         END IF
!-------------------------------------------------
!-------------------------------------------------

!-------------------------------------------------
       para_Tnum%Tana =.TRUE.
       IF (para_Tnum%Tana .AND. err_io /= 0) THEN
         write(out_unitp,*) "======================================"
         write(out_unitp,*) "======================================"
         write(out_unitp,*) "======================================"
         write(out_unitp,*) "======================================"
         CALL time_perso('Tana')

         CALL compute_analytical_KEO(TWOxKEO,mole,para_Tnum,Qact)

         IF (print_level > 2) CALL write_sum_opnd(TWOxKEO,header=.TRUE.)

         write(out_unitp,*) '================================================='
         write(out_unitp,*) ' Expand 2xKEO (in reduced dimension)'
         CALL Expand_Sum_OpnD_TO_Sum_OpnD(TWOxKEO, ExpandTWOxKEO)
         IF (print_level > 2)  CALL write_sum_opnd(ExpandTWOxKEO,header=.TRUE.)
         write(out_unitp,*) '================================================='

         CALL comparison_G_FROM_Tnum_Tana(ExpandTWOxKEO,mole,para_Tnum,Qact)

         CALL delete_op(TWOxKEO)
         CALL delete_op(ExpandTWOxKEO)


         ! calculation of the G matrix. Then print the diagonal elements
         CALL alloc_dnSVM(dnGG,mole%ndimG,mole%ndimG,mole%nb_act,nderiv=0)

         para_Tnum%WriteT    = .FALSE.
         CALL get_dng_dnGG(Qact,para_Tnum,mole,dnGG=dnGG,nderiv=0)

         write(out_unitp,*) 'Coordinate, value, GQQ'
         DO i=1,mole%nb_act
           write(out_unitp,*) String_TO_String('Q' // int_TO_char(i-1) ),Qact(i),dnGG%d0(i,i)
         END DO

         CALL dealloc_dnSVM(dnGG)

         CALL time_perso('Tana')
         write(out_unitp,*) "======================================"
         write(out_unitp,*) "======================================"
         write(out_unitp,*) "======================================"
         write(out_unitp,*) "======================================"
       END IF
!-------------------------------------------------


       CALL dealloc_zmat(mole)
       CALL dealloc_NParray(Qact,'Qact',name_sub)
       CALL dealloc_NParray(Qxyz,'Qxyz',name_sub)

       write(out_unitp,*) 'END ',name_sub

      END PROGRAM Tnum90_MidasCpp
