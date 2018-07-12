      PROGRAM Tnum_f90
      USE mod_system
      USE mod_Tnum
      USE mod_Tana_keo
      USE mod_Tana_Tnum
      USE mod_export_KEO
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

      real (kind=Rkind) :: vep,rho
      real (kind=Rkind), pointer :: Tdef2(:,:) => null()
      real (kind=Rkind), pointer :: Tdef1(:) => null()

      real (kind=Rkind), pointer :: Tcor2(:,:) => null()
      real (kind=Rkind), pointer :: Tcor1(:) => null()
      real (kind=Rkind), pointer :: Trot(:,:) => null()


      TYPE(Type_dnVec) :: dnx

      integer :: nderiv
      real (kind=Rkind), allocatable :: Qact(:)

!     - working parameters ------------------------------------------
      integer :: i,n

      character (len=*), parameter :: name_sub='Tnum90_MCTDH'

     !CALL test_FracInteger()
     !STOP
!===========================================================
!===========================================================
      !para_mem%mem_debug = .TRUE.
      CALL versionEVRT(.TRUE.)
      print_level=0
      !-----------------------------------------------------------------
      !     - read the coordinate tansformations :
      !     -   zmatrix, polysperical, bunch...
      !     ------------------------------------------------------------
      para_Tnum%LaTeXForm = .TRUE.
      CALL Read_mole(mole,para_Tnum,const_phys)
      para_Tnum%MCTDHForm = .TRUE.
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
      CALL time_perso('Tnum90_MCTDH')
      CALL Finalyze_TnumTana_Coord_PrimOp(para_Tnum,mole,para_PES)
      CALL time_perso('Tnum90_MCTDH')

      !-----------------------------------------------------------------
      n=1! (several) evaluations (for cpu time)
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
       nderiv = 0
       mole%WriteCC = .TRUE.
       CALL alloc_dnSVM(dnx,mole%ncart,mole%nb_act,nderiv)
       write(out_unitp,*) "======================================"
       CALL time_perso('dnx')

       DO i=1,n
         CALL sub_QactTOdnx(Qact,dnx,mole,nderiv,.FALSE.)
       END DO

       CALL time_perso('dnx')
       write(out_unitp,*) "======================================"

       write(out_unitp,*) 'dnx: ',mole%ncart
       CALL write_dnx(1,mole%ncart,dnx,nderiv)

       CALL sub_QactTOdnx(Qact,dnx,mole,nderiv,.TRUE.)

       write(out_unitp,*) 'dnx (mass weighted): ',mole%ncart
       CALL write_dnx(1,mole%ncart,dnx,nderiv)

       CALL sub_dnxNOMassWeight(dnx,mole%d0sm,mole%ncart,mole%ncart_act,nderiv)
       write(out_unitp,*) ' Cartesian coordinates with Eckart (au):'
       CALL write_dnx(1,mole%ncart,dnx,nderiv)
       write(out_unitp,*) ' Cartesian coordinates  with Eckart (ang):'
       CALL Write_Cartg98(dnx%d0,mole)

       mole%WriteCC = .FALSE.
       CALL dealloc_dnSVM(dnx)

       write(out_unitp,*) "======================================"
       write(out_unitp,*) "======================================"
       write(out_unitp,*) "======================================"
       write(out_unitp,*) "======================================"
       write(out_unitp,*) n,' evaluation of sub_QactTOdnx'
       write(out_unitp,*) "======================================"
!-------------------------------------------------
!-------------------------------------------------

!-------------------------------------------------
!-------------------------------------------------
!     - calculation of f2, f1, vep, rho ----------
!     --------------------------------------------

       CALL alloc_array(Tdef2,(/ mole%nb_act,mole%nb_act /),'Tdef2',name_sub)
       CALL alloc_array(Tdef1,(/ mole%nb_act /),            'Tdef1',name_sub)
       CALL alloc_array(Tcor2,(/ mole%nb_act,3 /),          'Tcor2',name_sub)
       CALL alloc_array(Tcor1,(/ 3 /),                      'Tcor1',name_sub)
       CALL alloc_array(Trot, (/ 3,3 /),                    'Trot', name_sub)

       write(out_unitp,*) "======================================"
       write(out_unitp,*) "======================================"
       write(out_unitp,*) "======================================"
       write(out_unitp,*) n,' evaluation of calc3_f2_f1Q_num'
       write(out_unitp,*) "======================================"
       para_Tnum%WriteT = .FALSE.
       CALL time_perso('f2')
       DO i=1,n-1

        CALL    calc3_f2_f1Q_num(Qact,                                  &
                                 Tdef2,Tdef1,vep,rho,                   &
                                 Tcor2,Tcor1,Trot,                      &
                                 para_Tnum,mole)
       END DO
       para_Tnum%WriteT = .TRUE.
       CALL    calc3_f2_f1Q_num(Qact,                                   &
                                 Tdef2,Tdef1,vep,rho,                   &
                                 Tcor2,Tcor1,Trot,                      &
                                 para_Tnum,mole)
       para_Tnum%WriteT = .FALSE.
       CALL time_perso('f2')
       write(out_unitp,*) "======================================"
       write(out_unitp,*) "======================================"
       write(out_unitp,*) "======================================"
       write(out_unitp,*) n,' evaluation of calc3_f2_f1Q_num'
       write(out_unitp,*) "======================================"

!-------------------------------------------------
!      FOR MCTDH (the calc3_f2_f1Q_num subroutine can be supressed)
       IF (para_Tnum%MCTDHform) CALL export3_MCTDH_T(Qact,para_Tnum,mole)
!-------------------------------------------------

       CALL dealloc_array(Tdef2,'Tdef2',name_sub)
       CALL dealloc_array(Tdef1,'Tdef1',name_sub)
       CALL dealloc_array(Tcor2,'Tcor2',name_sub)
       CALL dealloc_array(Tcor1,'Tcor1',name_sub)
       CALL dealloc_array(Trot, 'Trot', name_sub)

       CALL dealloc_zmat(mole)
       CALL dealloc_NParray(Qact,'Qact',name_sub)

       write(out_unitp,*) 'mem_tot,max_mem_used',para_mem%mem_tot,para_mem%max_mem_used
       write(out_unitp,*) 'nb_alloc,nb_dealloc',para_mem%nb_alloc,para_mem%nb_dealloc
       write(out_unitp,*) 'END Tnum'

      END PROGRAM Tnum_f90
