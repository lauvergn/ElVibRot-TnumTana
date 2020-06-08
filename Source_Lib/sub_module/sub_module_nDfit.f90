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
      MODULE mod_nDFit
      USE mod_system
      USE mod_nDindex
      IMPLICIT NONE

      PRIVATE

      TYPE nDFitAna_t

      ! range for the grid (1D, 2D)
      real (kind=Rkind), allocatable :: A(:)
      real (kind=Rkind), allocatable :: B(:)
      real (kind=Rkind), allocatable :: Step(:)
      integer,           allocatable :: nq(:)
      integer,           allocatable :: coord_list(:)

      logical :: Grid1D = .TRUE.
      logical :: Grid2D = .TRUE.

      logical :: Minimum   = .FALSE.
      logical :: all_coord = .TRUE.

      CONTAINS
        PROCEDURE, PRIVATE, PASS(nDFitAna1) :: nDFitAna2_TO_nDFitAna1
        GENERIC,   PUBLIC  :: assignment(=) => nDFitAna2_TO_nDFitAna1
      END TYPE nDFitAna_t


      TYPE param_nDFit ! it mays change in the futur (more like "basis" type)

      TYPE (Type_nDindex)            :: nDindB  ! enable to use multidimensional index for the basis
      integer                        :: ndim  = 0             ! size of Q0, nDsize, nDweight
      real (kind=Rkind), allocatable :: Q0(:)
      integer,           allocatable :: nDsize(:)
      real (kind=Rkind), allocatable :: nDweight(:)
      integer,           allocatable :: ntyp(:)

      integer                        :: nb_WB = 0
      real (kind=Rkind), allocatable :: B(:)
      integer,           allocatable :: nDvalB(:,:)
      integer                        :: MinCoupling    = 0
      integer                        :: MaxCoupling    = 4
      integer                        :: max_nb         = 10
      integer                        :: ind_val        = 1
      integer                        :: nb_val         = 1
      integer                        :: MR_order       = -1

      real (kind=Rkind)              :: MinNorm        = ZERO
      real (kind=Rkind)              :: MaxNorm        = FOUR
      real (kind=Rkind)              :: epsi           = ONETENTH**10
      real (kind=Rkind)              :: epsi_inter     = ONETENTH**10
      logical                        :: svd            = .TRUE.

      logical                        :: Analysis       = .FALSE.
      TYPE (nDFitAna_t)              :: nDFitAna

      integer                        :: Col_FOR_WeightOFFit  = 0 ! it is not used
      real (kind=Rkind)              :: Scal_FOR_WeightOFFit = 200._Rkind

      TYPE (param_file)              :: Param_Fit_file
      character (len=Line_len)       :: name_Fit = ''

      integer                        :: nb_fit = 0
      TYPE (param_nDFit), pointer    :: Tab_para_nDFit(:) => null()

      CONTAINS
        PROCEDURE, PRIVATE, PASS(nDFit1) :: nDFit2_TO_nDFit1
        GENERIC,   PUBLIC  :: assignment(=) => nDFit2_TO_nDFit1
      END TYPE param_nDFit

      PUBLIC :: nDFitAna_t, Read_Analysis
      PUBLIC :: param_nDFit, dealloc_nDFit,Read_nDFit,ReadWrite_nDFitW, &
                Write_nDFit, Analysis_nDFit, Analysis_nDFitW,           &
                Read_FOR_nDFit1_TO_TnDFit2, sub_nDFunc_FROM_nDFit,      &
                sub_ONLYnDFunc_FROM_nDFit, nDFunct_WITH_Q, nDFunct_WITH_tQ

      CONTAINS

      RECURSIVE SUBROUTINE nDFit2_TO_nDFit1(nDFit1,nDFit2)
      USE mod_system
      IMPLICIT NONE

      CLASS (param_nDFit), intent(inout) :: nDFit1
      TYPE (param_nDFit),  intent(in)    :: nDFit2

      integer :: i
!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'nDFit2_TO_nDFit1'
!      logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------

      !write(out_unitp,*) 'SUBROUTINE ',trim(name_sub)

      IF (allocated(nDFit2%Q0))       nDFit1%Q0       = nDFit2%Q0
      IF (allocated(nDFit2%nDsize))   nDFit1%nDsize   = nDFit2%nDsize
      IF (allocated(nDFit2%nDweight)) nDFit1%nDweight = nDFit2%nDweight
      IF (allocated(nDFit2%ntyp))     nDFit1%ntyp     = nDFit2%ntyp
      IF (allocated(nDFit2%B))        nDFit1%B        = nDFit2%B
      IF (allocated(nDFit2%nDvalB))   nDFit1%nDvalB   = nDFit2%nDvalB

      nDFit1%nDindB               = nDFit2%nDindB
      nDFit1%ndim                 = nDFit2%ndim
      nDFit1%nb_WB                = nDFit2%nb_WB

      nDFit1%MinCoupling          = nDFit2%MinCoupling
      nDFit1%MaxCoupling          = nDFit2%MaxCoupling
      nDFit1%max_nb               = nDFit2%max_nb
      nDFit1%ind_val              = nDFit2%ind_val
      nDFit1%nb_val               = nDFit2%nb_val
      nDFit1%MR_order             = nDFit2%MR_order

      nDFit1%MinNorm              = nDFit2%MinNorm
      nDFit1%MaxNorm              = nDFit2%MaxNorm
      nDFit1%epsi                 = nDFit2%epsi
      nDFit1%epsi_inter           = nDFit2%epsi_inter
      nDFit1%svd                  = nDFit2%svd

      nDFit1%Analysis             = nDFit2%Analysis
      nDFit1%nDFitAna             = nDFit2%nDFitAna

      nDFit1%Col_FOR_WeightOFFit  = nDFit2%Col_FOR_WeightOFFit
      nDFit1%Scal_FOR_WeightOFFit = nDFit2%Scal_FOR_WeightOFFit

      nDFit1%Param_Fit_file       = nDFit2%Param_Fit_file
      nDFit1%name_Fit             = nDFit2%name_Fit

      nDFit1%nb_fit               = nDFit2%nb_fit

      IF (associated(nDFit2%Tab_para_nDFit)) THEN
        allocate(nDFit1%Tab_para_nDFit(size(nDFit2%Tab_para_nDFit)))
        DO i=1,size(nDFit2%Tab_para_nDFit)
          nDFit1%Tab_para_nDFit(i) = nDFit2%Tab_para_nDFit(i)
        END DO
      END IF

      END SUBROUTINE nDFit2_TO_nDFit1
      RECURSIVE SUBROUTINE dealloc_nDFit(nDFit)
      USE mod_system
      IMPLICIT NONE

      CLASS (param_nDFit), intent(inout) :: nDFit

      integer :: i
!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'dealloc_nDFit'
!      logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------

      !write(out_unitp,*) 'SUBROUTINE ',trim(name_sub)

      IF (allocated(nDFit%Q0))       deallocate(nDFit%Q0)
      IF (allocated(nDFit%nDsize))   deallocate(nDFit%nDsize)
      IF (allocated(nDFit%nDweight)) deallocate(nDFit%nDweight)
      IF (allocated(nDFit%ntyp))     deallocate(nDFit%ntyp)
      IF (allocated(nDFit%B))        deallocate(nDFit%B)
      IF (allocated(nDFit%nDvalB))   deallocate(nDFit%nDvalB)

      CALL dealloc_nDindex(nDFit%nDindB)

      nDFit%ndim                 = 0
      nDFit%nb_WB                = 0

      nDFit%MinCoupling          = 0
      nDFit%MaxCoupling          = 4
      nDFit%max_nb               = 10
      nDFit%ind_val              = 1
      nDFit%nb_val               = 1
      nDFit%MR_order             = -1

      nDFit%MinNorm              = ZERO
      nDFit%MaxNorm              = FOUR
      nDFit%epsi                 = ONETENTH**10
      nDFit%epsi_inter           = ONETENTH**10
      nDFit%svd                  = .TRUE.

      nDFit%Analysis             = .FALSE.
      CALL dealloc_nDFitAna(nDFit%nDFitAna)

      nDFit%Col_FOR_WeightOFFit  = 0
      nDFit%Scal_FOR_WeightOFFit = 200._Rkind

      CALL file_dealloc(nDFit%Param_Fit_file)
      nDFit%name_Fit             = ''

      nDFit%nb_fit               = 0

      IF (associated(nDFit%Tab_para_nDFit)) THEN
        DO i=1,size(nDFit%Tab_para_nDFit)
          CALL dealloc_nDFit(nDFit%Tab_para_nDFit(i))
        END DO
        deallocate(nDFit%Tab_para_nDFit)
        nullify(nDFit%Tab_para_nDFit)
      END IF

      END SUBROUTINE dealloc_nDFit
      SUBROUTINE nDFitAna2_TO_nDFitAna1(nDFitAna1,nDFitAna2)
      USE mod_system
      IMPLICIT NONE

      CLASS (nDFitAna_t), intent(inout) :: nDFitAna1
      TYPE (nDFitAna_t),  intent(in)    :: nDFitAna2

!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'nDFitAna2_TO_nDFitAna1'
!      logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------

      !write(out_unitp,*) 'SUBROUTINE ',trim(name_sub)

      ! range for the grid (1D, 2D)
      IF (allocated(nDFitAna2%A))          nDFitAna1%A          = nDFitAna2%A
      IF (allocated(nDFitAna2%B))          nDFitAna1%B          = nDFitAna2%B
      IF (allocated(nDFitAna2%Step))       nDFitAna1%Step       = nDFitAna2%Step
      IF (allocated(nDFitAna2%nq))         nDFitAna1%nq         = nDFitAna2%nq
      IF (allocated(nDFitAna2%coord_list)) nDFitAna1%coord_list = nDFitAna2%coord_list

      nDFitAna1%Grid1D    = nDFitAna2%Grid1D
      nDFitAna1%Grid2D    = nDFitAna2%Grid2D
      nDFitAna1%Minimum   = nDFitAna2%Minimum
      nDFitAna1%all_coord = nDFitAna2%all_coord

      END SUBROUTINE nDFitAna2_TO_nDFitAna1
      SUBROUTINE dealloc_nDFitAna(nDFitAna)
      USE mod_system
      IMPLICIT NONE

      CLASS (nDFitAna_t), intent(inout) :: nDFitAna

!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'dealloc_nDFitAna'
!      logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------

      !write(out_unitp,*) 'SUBROUTINE ',trim(name_sub)

      ! range for the grid (1D, 2D)
      IF (allocated(nDFitAna%A))          deallocate(nDFitAna%A)
      IF (allocated(nDFitAna%B))          deallocate(nDFitAna%B)
      IF (allocated(nDFitAna%Step))       deallocate(nDFitAna%Step)
      IF (allocated(nDFitAna%nq))         deallocate(nDFitAna%nq)
      IF (allocated(nDFitAna%coord_list)) deallocate(nDFitAna%coord_list)

      nDFitAna%Grid1D    = .TRUE.
      nDFitAna%Grid2D    = .TRUE.
      nDFitAna%Minimum   = .FALSE.
      nDFitAna%all_coord = .TRUE.

      END SUBROUTINE dealloc_nDFitAna
      SUBROUTINE Read_Analysis(nDFitAna,Q0)
      USE mod_system
      IMPLICIT NONE

      TYPE (nDFitAna_t), intent(inout) :: nDFitAna
      real (kind=Rkind), intent(in)        :: Q0(:)

      integer                    :: i,ndim
      character (len=Name_len)   :: name_dum,name_int


!----- For the namelist ----------------------------------------------
      logical :: Grid1D
      logical :: Grid2D
      logical :: Minimum,DeltaRange_Read,all_coord

      namelist /Analysis/ Grid1D,Grid2D,Minimum,DeltaRange_Read,all_coord


!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'Read_Analysis'
!      logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
!
!=====================================================================
        DeltaRange_Read         = .TRUE.

        Grid1D                  = .TRUE.
        Grid2D                  = .TRUE.
        Minimum                 = .TRUE.
        all_coord               = .TRUE.

        read(in_unitp,Analysis,IOSTAT=err_read)
        IF (err_read < 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' End-of-file or End-of-record'
          write(out_unitp,*) ' The namelist "Analysis" is probably absent'
          write(out_unitp,*) ' check your data!'
          write(out_unitp,*) ' ERROR in ',name_sub
          STOP
        ELSE IF (err_read > 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' Some parameter name of the namelist "Analysis" are probaly wrong'
          write(out_unitp,*) ' check your data!'
          write(out_unitp,Analysis)
          write(out_unitp,*) ' ERROR in ',name_sub
          STOP
        END IF
        IF (debug) write(out_unitp,Analysis)


       ndim = size(Q0)
       IF (size(Q0) < 1) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) '  The size of Q0(:) is < 1'
            write(out_unitp,*) '  shap(Q0)',shape(Q0)
            write(out_unitp,*) '  It should never append'
            write(out_unitp,*) '  Check the fortran !'
            write(out_unitp,*) ' ERROR in ',name_sub
            STOP
        END IF


        nDFitAna%Grid1D              = Grid1D
        nDFitAna%Grid2D              = Grid2D
        nDFitAna%Minimum             = Minimum
        nDFitAna%all_coord           = all_coord


        CALL alloc_NParray(nDFitAna%nq,shape(Q0),'nDFitAna%nq',name_sub)
        CALL alloc_NParray(nDFitAna%A,shape(Q0),'nDFitAna%A',name_sub)
        CALL alloc_NParray(nDFitAna%B,shape(Q0),'nDFitAna%B',name_sub)
        CALL alloc_NParray(nDFitAna%Step,shape(Q0),'nDFitAna%Step',name_sub)
        CALL alloc_NParray(nDFitAna%coord_list,shape(Q0),'nDFitAna%coord_list',name_sub)

        nDFitAna%A(:)  = ZERO
        nDFitAna%B(:)  = ZERO
        nDFitAna%nq(:) = 0

        IF (.NOT. all_coord) THEN
          nDFitAna%coord_list(:) = 0

          DO i=1,ndim
           CALL read_name_advNo(in_unitp,name_int,err_read)

           IF (len_trim(name_int) == 0) EXIT
             !write(out_unitp,*) 'i,err_io',i,err_io
             !write(out_unitp,*) 'i,name_int',i,name_int
             read(name_int,*) nDFitAna%coord_list(i)
             IF (err_read /= 0) EXIT ! end of the liste

          END DO
          write(out_unitp,*) 'coord_list',nDFitAna%coord_list(:)
        ELSE
          nDFitAna%coord_list(:) = (/ (i,i=1,ndim) /)
        END IF

        read(in_unitp,*) name_dum,nDFitAna%nq(:)
        read(in_unitp,*) name_dum,nDFitAna%A(:)
        read(in_unitp,*) name_dum,nDFitAna%B(:)

        nDFitAna%Step(:) = (nDFitAna%B(:)-nDFitAna%A(:)) /              &
                                 real(nDFitAna%nq(:)-1,kind=Rkind)


        IF (DeltaRange_Read) THEN
           nDFitAna%A(:) = Q0 + nDFitAna%A(:)
           nDFitAna%B(:) = Q0 + nDFitAna%B(:)
        END IF

        write(out_unitp,*) 'nq:',nDFitAna%nq(:)
        CALL Write_VecMat(nDFitAna%A,out_unitp,5,name_info='range A:')
        CALL Write_VecMat(nDFitAna%B,out_unitp,5,name_info='range B:')
        CALL Write_VecMat(nDFitAna%Step,out_unitp,5,name_info='Step:')

        CALL flush_perso(out_unitp)


      END SUBROUTINE Read_Analysis

      SUBROUTINE Read_nDFit(para_nDFit,Q0)
      USE mod_system
      IMPLICIT NONE

      real (kind=Rkind), intent(in)     :: Q0(:)
      TYPE (param_nDFit), intent(inout) :: para_nDFit

      integer           :: nb_act


!----- For the namelist ----------------------------------------------
      integer           :: MinCoupling,MaxCoupling
      real (kind=Rkind) :: MinNorm,MaxNorm,conv_col,max_b,Weight_iGP
      logical           :: svd,ntyp_read,Analysis
      real (kind=Rkind) :: epsi,epsi_inter
      integer           :: ind_val,nb_val,nb_G,max_nb,MR_order
      integer           :: Col_FOR_WeightOFFit
      real (kind=Rkind) :: Scal_FOR_WeightOFFit

      namelist /nDFit/ MinNorm,MaxNorm,MinCoupling,MaxCoupling,         &
                       ind_val,nb_val,svd,ntyp_read,                    &
                       epsi,epsi_inter,max_nb,                          &
                       Col_FOR_WeightOFFit,Scal_FOR_WeightOFFit,MR_order,&
                       Analysis


!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'Read_nDFit'
!      logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
!
!=====================================================================
        MR_order             = -1  ! order of the multimode representation (-1: not use)
        MaxNorm              = FOUR
        MinNorm              = 0
        MaxCoupling          = 4
        MinCoupling          = 0
        svd                  = .TRUE.
        epsi                 = ONETENTH**10
        epsi_inter           = ONETENTH**5
        nb_val               = 1
        ind_val              = 1
        max_nb               = 10
        Col_FOR_WeightOFFit  = 0
        Scal_FOR_WeightOFFit = 200._Rkind
        ntyp_read            = .FALSE.
        Analysis             = .FALSE.

        read(in_unitp,nDFit,IOSTAT=err_read)
        IF (err_read < 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' End-of-file or End-of-record'
          write(out_unitp,*) ' The namelist "nDFit" is probably absent'
          write(out_unitp,*) ' check your data!'
          write(out_unitp,*) ' ERROR in ',name_sub
          STOP
        ELSE IF (err_read > 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' Some parameter name of the namelist "nDFit" are probaly wrong'
          write(out_unitp,*) ' check your data!'
          write(out_unitp,nDFit)
          write(out_unitp,*) ' ERROR in ',name_sub
          STOP
        END IF
        IF (debug) write(out_unitp,nDFit)

        IF (nb_val < 1 .OR. ind_val < 1 .OR. ind_val > nb_val) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nb_val or ind_val are probaly wrong'
          write(out_unitp,*) ' nb_val, ind_val: ',nb_val,ind_val
          write(out_unitp,*) ' check your data!'
          write(out_unitp,nDFit)
          write(out_unitp,*) ' ERROR in ',name_sub
          STOP
        END IF

        para_nDFit%MR_order              = MR_order
        para_nDFit%nb_val                = nb_val
        para_nDFit%ind_val               = ind_val
        para_nDFit%MinNorm               = MinNorm
        para_nDFit%MaxNorm               = MaxNorm
        para_nDFit%MinCoupling           = MinCoupling
        para_nDFit%MaxCoupling           = MaxCoupling
        para_nDFit%svd                   = svd
        para_nDFit%max_nb                = max_nb
        para_nDFit%epsi                  = epsi
        para_nDFit%epsi_inter            = epsi_inter
        para_nDFit%Col_FOR_WeightOFFit   = Col_FOR_WeightOFFit
        para_nDFit%Scal_FOR_WeightOFFit  = Scal_FOR_WeightOFFit
        para_nDFit%Analysis              = Analysis

       para_nDFit%ndim                   = size(Q0)

       IF (para_nDFit%ndim < 1) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) '  The size of Q0(:) is < 1'
            write(out_unitp,*) '  It should never append'
            write(out_unitp,*) '  Check the fortran !'
            write(out_unitp,*) ' ERROR in ',name_sub
            STOP
        END IF


        CALL alloc_NParray(para_nDFit%Q0,shape(Q0),'para_nDFit%Q0',name_sub)
        para_nDFit%Q0(:) = Q0(:)

        CALL alloc_NParray(para_nDFit%nDweight,shape(Q0),                 &
                        'para_nDFit%nDweight',name_sub)
        CALL alloc_NParray(para_nDFit%nDsize,shape(Q0),                   &
                        'para_nDFit%nDsize',name_sub)
        CALL alloc_NParray(para_nDFit%ntyp,shape(Q0),                     &
                        'para_nDFit%ntyp',name_sub)
        para_nDFit%ntyp(:) = 15 ! polynomial

        read(in_unitp,*) para_nDFit%nDweight(:)
        read(in_unitp,*) para_nDFit%nDsize(:)
        IF (ntyp_read) read(in_unitp,*) para_nDFit%ntyp

        write(out_unitp,*) para_nDFit%nDweight(:)
        write(out_unitp,*) para_nDFit%nDsize(:)
        write(out_unitp,*) para_nDFit%ntyp(:)

        IF (Analysis) THEN
          CALL Read_Analysis(para_nDFit%nDFitAna,Q0)
        END IF

        CALL flush_perso(out_unitp)


      END SUBROUTINE Read_nDFit

      SUBROUTINE Write_nDFit(para_nDFit)
      USE mod_system
      IMPLICIT NONE

      TYPE (param_nDFit), intent(in) :: para_nDFit

!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'Write_nDFit'
!      logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------

      write(out_unitp,*)  '============================================================'
      write(out_unitp,*)  'BEGINING Write_nDFit'
      write(out_unitp,*)  '  ---------------------------------------------------------'

      write(out_unitp,*)  'ndim',para_nDFit%ndim
      write(out_unitp,*)  'allo Q0',allocated(para_nDFit%Q0)
      IF (allocated(para_nDFit%Q0)) write(out_unitp,*)  'Q0',para_nDFit%Q0
      write(out_unitp,*)  'allo nDsize',allocated(para_nDFit%nDsize)
      IF (allocated(para_nDFit%nDsize)) write(out_unitp,*)  'nDsize',para_nDFit%nDsize
      write(out_unitp,*)  'allo nDweight',allocated(para_nDFit%nDweight)
      IF (allocated(para_nDFit%nDweight)) write(out_unitp,*)  'nDweight',para_nDFit%nDweight
      write(out_unitp,*)  'allo ntyp',allocated(para_nDFit%ntyp)
      IF (allocated(para_nDFit%ntyp)) write(out_unitp,*)  'ntyp',para_nDFit%ntyp

      write(out_unitp,*)  'nDindB'
      CALL Write_nDindex(para_nDFit%nDindB,'para_nDFit')
      write(out_unitp,*)  '  ---------------------------------------------------------'


      write(out_unitp,*)  'nb_WB',para_nDFit%nb_WB
      write(out_unitp,*)  'allo B',allocated(para_nDFit%B)
      IF (allocated(para_nDFit%nDvalB)) write(out_unitp,*)  'nDvalB',para_nDFit%nDvalB
      write(out_unitp,*)  'allo nDvalB',allocated(para_nDFit%nDvalB)
      IF (allocated(para_nDFit%nDvalB)) write(out_unitp,*)  'nDvalB',para_nDFit%nDvalB
      write(out_unitp,*)  '  ---------------------------------------------------------'

      write(out_unitp,*)  'MinCoupling,MaxCoupling',para_nDFit%MinCoupling,para_nDFit%MaxCoupling
      write(out_unitp,*)  'max_nb',para_nDFit%max_nb
      write(out_unitp,*)  'ind_val',para_nDFit%ind_val
      write(out_unitp,*)  'nb_val',para_nDFit%nb_val
      write(out_unitp,*)  'MR_order',para_nDFit%MR_order
      write(out_unitp,*)  'MinNorm,MaxNorm',para_nDFit%MinNorm,para_nDFit%MaxNorm
      write(out_unitp,*)  'epsi,epsi_inter',para_nDFit%epsi,para_nDFit%epsi_inter
      write(out_unitp,*)  'svd',para_nDFit%svd
      write(out_unitp,*)  '  ---------------------------------------------------------'



      write(out_unitp,*)  'Analysis',para_nDFit%Analysis
      !CALL Write_nDindex(para_nDFit%nDFitAna)
      write(out_unitp,*)  '  ---------------------------------------------------------'

      write(out_unitp,*)  'Col_FOR_WeightOFFit',para_nDFit%Col_FOR_WeightOFFit
      write(out_unitp,*)  'Scal_FOR_WeightOFFit',para_nDFit%Scal_FOR_WeightOFFit
      write(out_unitp,*)  '  ---------------------------------------------------------'


      write(out_unitp,*)  'name_Fit: ',trim(adjustl(para_nDFit%name_Fit))
      CALL file_Write(para_nDFit%Param_Fit_file)
      write(out_unitp,*)  '  ---------------------------------------------------------'


      write(out_unitp,*)  'nb_fit: ',para_nDFit%nb_fit
      !TYPE (param_nDFit), pointer :: Tab_para_nDFit(:) => null()
      write(out_unitp,*)  '  ---------------------------------------------------------'

      write(out_unitp,*)  'END Write_nDFit'
      write(out_unitp,*)  '============================================================'

      CALL flush_perso(out_unitp)


      END SUBROUTINE Write_nDFit

      RECURSIVE SUBROUTINE ReadWrite_nDFitW(para_nDFit,ReadData,B,Tab_nDval,conv_ene)
      USE mod_system
      IMPLICIT NONE

      logical, intent(in)               :: ReadData  ! if true read, else write
      TYPE (param_nDFit), intent(inout) :: para_nDFit
      real (kind=Rkind), optional       :: B(:)
      integer, optional                 :: Tab_nDval(:,:)
      real (kind=Rkind), optional       :: conv_ene


      integer                    :: i,iB,idum,nioFit,max_b
      character (len=Name_len)   :: name_dum,name_int
      real (kind=Rkind)          :: conv_col,conv_ene_loc


      !-----------------------------------------------------------------
      ! for the namelist
      integer                  :: MinCoupling,MaxCoupling
      real (kind=Rkind)        :: MinNorm,MaxNorm,max_valB
      integer                  :: ndim,nb_B,nb_WB,MR_order
      logical                  :: svd
      real (kind=Rkind)        :: epsi,epsi_inter
      integer                  :: ind_val,nb_val,max_nb
      integer                  :: Col_FOR_WeightOFFit
      real (kind=Rkind)        :: Scal_FOR_WeightOFFit
      integer                  :: nb_Fit

      namelist /nDFitW/ MinNorm,MaxNorm,MinCoupling,MaxCoupling,        &
                        ind_val,nb_val,svd,                             &
                        epsi,epsi_inter,max_nb,ndim,MR_order,           &
                        Col_FOR_WeightOFFit,Scal_FOR_WeightOFFit,       &
                        nb_B,nb_WB,                                     &
                        nb_Fit
      ! for the namelist
      !-----------------------------------------------------------------


!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'ReadWrite_nDFitW'
!      logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
!
!=====================================================================
      IF (present(conv_ene)) THEN
        conv_ene_loc = conv_ene
      ELSE
        conv_ene_loc = ONE
      END IF


      IF (ReadData) THEN
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "=== READ PARAM FIT ==================="
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "  File: ",para_nDFit%Param_Fit_file%name
        IF (len_trim(para_nDFit%Param_Fit_file%name) == 0) STOP 'file name of para_nDFit is empty!!'

        CALL file_open(para_nDFit%Param_Fit_file,nioFit,old=.TRUE.)

        MR_order     = -1
        ndim         =  0
        MinCoupling  =  0
        MinNorm      =  ZERO
        nb_Fit       =  0
        read(nioFit,nDFitW,IOSTAT=err_read)
        IF (err_read < 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' End-of-file or End-of-record'
          write(out_unitp,*) ' The namelist "nDFitW" is probably absent'
          write(out_unitp,*) ' from the file: ',trim(adjustl(para_nDFit%Param_Fit_file%name))
          write(out_unitp,*) ' ERROR in ',name_sub
          STOP
        ELSE IF (err_read > 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' Some parameter name of the namelist "nDFitW" are probaly wrong'
          write(out_unitp,*) ' in the file: ',trim(adjustl(para_nDFit%Param_Fit_file%name))
          write(out_unitp,*) ' It should never append !!'
          write(out_unitp,*) ' Check the fortran'
          write(out_unitp,*) ' ERROR in ',name_sub
          STOP
        END IF
        IF (debug) write(out_unitp,nDFitW)
        para_nDFit%nb_Fit = nb_Fit

        IF (nb_Fit > 0) THEN
          allocate(para_nDFit%Tab_para_nDFit(para_nDFit%nb_Fit))
          DO i=1,para_nDFit%nb_Fit
            read(nioFit,*) para_nDFit%Tab_para_nDFit(i)%name_Fit
            para_nDFit%Tab_para_nDFit(i)%Param_Fit_file%name =          &
                                    para_nDFit%Tab_para_nDFit(i)%name_Fit
          END DO
          CALL file_close(para_nDFit%Param_Fit_file)
          DO i=1,para_nDFit%nb_Fit
            CALL ReadWrite_nDFitW(para_nDFit%Tab_para_nDFit(i),.TRUE.)
          END DO
          para_nDFit%ndim = para_nDFit%Tab_para_nDFit(1)%ndim
          para_nDFit%Q0   = para_nDFit%Tab_para_nDFit(1)%Q0
        ELSE
          write(out_unitp,*) "  ndim (nb_act): ",ndim
          write(out_unitp,*) "  nb read functions: ",nb_WB

          IF (ndim < 1) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' Some parameter name of the namelist "nDFitW" are probaly wrong'
            write(out_unitp,*) ' in the file: ',trim(adjustl(para_nDFit%Param_Fit_file%name))
            write(out_unitp,*) '  Check the file and add the (correct) value of ndim'
            write(out_unitp,*) ' ERROR in ',name_sub
            STOP
          END IF

          IF (.NOT. allocated(para_nDFit%Q0)) THEN
            CALL alloc_NParray(para_nDFit%Q0,(/ndim/),                  &
                              'para_nDFit%Q0',name_sub)
          END IF
          IF (.NOT. allocated(para_nDFit%nDweight)) THEN
            CALL alloc_NParray(para_nDFit%nDweight,(/ndim/),            &
                              'para_nDFit%nDweight',name_sub)
          END IF
          IF (.NOT. allocated(para_nDFit%nDsize)) THEN
            CALL alloc_NParray(para_nDFit%nDsize,(/ndim/),              &
                              'para_nDFit%nDsize',name_sub)
          END IF
          IF (.NOT. allocated(para_nDFit%ntyp)) THEN
            CALL alloc_NParray(para_nDFit%ntyp,(/ndim/),                &
                              'para_nDFit%ntyp',name_sub)
          END IF

          read(nioFit,*) name_dum,para_nDFit%Q0 ! for Q0
          read(nioFit,*) name_dum,para_nDFit%nDweight ! for nDweight
          read(nioFit,*) name_dum,para_nDFit%nDsize ! for nDsize
          read(nioFit,*) name_dum,para_nDFit%ntyp ! for ntyp
          write(out_unitp,*) 'Q0       ',para_nDFit%Q0
          write(out_unitp,*) 'nDweight ',para_nDFit%nDweight
          write(out_unitp,*) 'nDsize   ',para_nDFit%nDsize
          write(out_unitp,*) 'ntyp     ',para_nDFit%ntyp
          CALL flush_perso(out_unitp)

          para_nDFit%nb_WB = nb_WB
          para_nDFit%ndim  = ndim

          para_nDFit%MR_order              = MR_order
          para_nDFit%nb_val                = nb_val
          para_nDFit%ind_val               = ind_val
          para_nDFit%MinNorm               = MinNorm
          para_nDFit%MaxNorm               = MaxNorm
          para_nDFit%MinCoupling           = MinCoupling
          para_nDFit%MaxCoupling           = MaxCoupling
          para_nDFit%svd                   = svd
          para_nDFit%max_nb                = max_nb
          para_nDFit%epsi                  = epsi
          para_nDFit%epsi_inter            = epsi_inter
          para_nDFit%Col_FOR_WeightOFFit   = Col_FOR_WeightOFFit
          para_nDFit%Scal_FOR_WeightOFFit  = Scal_FOR_WeightOFFit

          CALL flush_perso(out_unitp)

          CALL alloc_NParray(para_nDFit%B,(/nb_WB/),                    &
                            'para_nDFit%B',name_sub)
          CALL alloc_NParray(para_nDFit%nDvalB,(/ndim,nb_WB/),          &
                            'para_nDFit%nDvalB',name_sub)

          DO iB=1,nb_WB
            read(nioFit,*) idum,para_nDFit%nDvalB(:,iB),para_nDFit%B(iB)
          END DO

          CALL file_close(para_nDFit%Param_Fit_file)
        END IF
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
      ELSE
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "=== EXPORT PARAM FIT ================="
        write(out_unitp,*) "======================================"
        IF (.NOT. present(b)) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' To export B.'
          write(out_unitp,*) ' B MUST be present in input argument list of the called subroutine'
          write(out_unitp,*) '  Check the fortran'
          write(out_unitp,*) ' ERROR in ',name_sub
          STOP
        END IF

        IF (para_nDFit%ind_val == 1) THEN
          conv_col = conv_ene_loc
        ELSE
          conv_col = ONE
        END IF

        nb_B  = size(B)
        max_valB = maxval(abs(b))
        write(out_unitp,*) 'Largest b value',max_valB

        nb_WB = count(abs(b(:)/max_valB) > para_nDFit%epsi_inter)

        write(out_unitp,*) 'number of large b value',nb_WB

        IF (.NOT. present(Tab_nDval)) THEN
          DO iB=1,nb_B
            IF (abs(b(ib)/max_valB) < para_nDFit%epsi_inter) THEN
                b(iB) = ZERO
            ELSE
              IF (count(para_nDFit%nDindB%Tab_nDval(:,iB)/=0)==1) THEN
                write(out_unitp,*) iB,para_nDFit%nDindB%Tab_nDval(:,iB),b(ib)*conv_col
              END IF
            END IF
          END DO
        END IF

        ndim                 = para_nDFit%ndim
        MR_order             = para_nDFit%MR_order
        nb_val               = para_nDFit%nb_val
        ind_val              = para_nDFit%ind_val
        MinNorm              = para_nDFit%MinNorm
        MaxNorm              = para_nDFit%MaxNorm
        MinCoupling          = para_nDFit%MinCoupling
        MaxCoupling          = para_nDFit%MaxCoupling
        svd                  = para_nDFit%svd
        max_nb               = para_nDFit%max_nb
        epsi                 = para_nDFit%epsi
        epsi_inter           = para_nDFit%epsi_inter
        Col_FOR_WeightOFFit  = para_nDFit%Col_FOR_WeightOFFit
        Scal_FOR_WeightOFFit = para_nDFit%Scal_FOR_WeightOFFit
        nb_Fit               = 0

        CALL Write_int_IN_char(ind_val,name_int)
        para_nDFit%Param_Fit_file%name =                                &
           trim(adjustl(para_nDFit%name_fit)) // trim(adjustl(name_int))
        write(out_unitp,*) 'name_fit_file: ',trim(para_nDFit%Param_Fit_file%name)

        CALL file_open(para_nDFit%Param_Fit_file,nioFit)


        write(nioFit,nDFitW)
        write(nioFit,*) 'Q0 ',para_nDFit%Q0(:)
        write(nioFit,*) 'nDweight ',para_nDFit%nDweight(:)
        write(nioFit,*) 'nDsize ',para_nDFit%nDsize(:)
        write(nioFit,*) 'ntyp ',para_nDFit%ntyp(:)

        IF (present(Tab_nDval)) THEN
          DO iB=1,nb_B
            IF (b(iB) /= ZERO) THEN
              write(nioFit,*) iB,Tab_nDval(:,iB),b(iB)
            END IF
          END DO
        ELSE
          DO iB=1,nb_B
            IF (b(iB) /= ZERO) THEN
              write(nioFit,*) iB,para_nDFit%nDindB%Tab_nDval(:,iB),b(iB)
            END IF
          END DO
        END IF


        CALL file_close(para_nDFit%Param_Fit_file)
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"
      END IF


      END SUBROUTINE ReadWrite_nDFitW

      RECURSIVE SUBROUTINE Analysis_nDFitW(para_nDFit,conv_ene)
      USE mod_system
      IMPLICIT NONE

      TYPE (param_nDFit), intent(inout) :: para_nDFit
      real (kind=Rkind), intent(in)     :: conv_ene

!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'Analysis_nDFitW'
!      logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
!
!=====================================================================


        write(out_unitp,*) "======================================"
        write(out_unitp,*) "=== ANALYSIS PARAM FIT ==============="
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "  File: ",para_nDFit%Param_Fit_file%name
        IF (para_nDFit%Param_Fit_file%name == '') THEN
          para_nDFit%Param_Fit_file%name = para_nDFit%name_Fit
        END IF


        IF (allocated(para_nDFit%B))                                   &
               CALL dealloc_NParray(para_nDFit%B,'para_nDFit%B',name_sub)
        IF (allocated(para_nDFit%nDvalB))                              &
           CALL dealloc_NParray(para_nDFit%nDvalB,'para_nDFit%nDvalB',name_sub)

        CALL ReadWrite_nDFitW(para_nDFit,.TRUE.)

        IF (para_nDFit%Analysis) THEN
          CALL Analysis_nDFit(para_nDFit,conv_ene)
        END IF

        write(out_unitp,*) "======================================"
        write(out_unitp,*) "=== END ANALYSIS PARAM FIT ==========="
        write(out_unitp,*) "======================================"


      END SUBROUTINE Analysis_nDFitW

      SUBROUTINE Analysis_nDFit(para_nDFit,conv_ene)
      USE mod_system
      IMPLICIT NONE

      TYPE (param_nDFit), intent(inout) :: para_nDFit
      real (kind=Rkind), intent(in)     :: conv_ene


      integer                    :: ii,i,iB,ndim_coord
      integer                    :: ic1,ic2,ic3, i1,i2,i3 ,iq1,iq2,iq3

      real (kind=Rkind)          :: conv_col,val_nDfit
      real (kind=Rkind)          :: Q(para_nDFit%ndim)
      TYPE (param_file)          :: Grid1D_file
      TYPE (param_file)          :: Grid2D_file

      real (kind=Rkind)          :: val0,Q10,Q20,Q30
      real (kind=Rkind)          :: val_min,Q1_min,Q2_min,Q3_min,val_min_1D
      real (kind=Rkind)          :: val_max,Q1_max,Q2_max,Q3_max

!----- for debuging --------------------------------------------------
      integer :: err_nio,nio
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'Analysis_nDFit'
!      logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
!
!=====================================================================


        write(out_unitp,*) "======================================"
        write(out_unitp,*) "=== ANALYSIS PARAM FIT ==============="
        write(out_unitp,*) "======================================"

        IF (para_nDFit%ind_val == 1) THEN
          conv_col = conv_ene
        ELSE
          conv_col = ONE
        END IF



        DO iB=1,para_nDFit%nb_WB
          i1 = 0
          i2 = 0
          DO i=1,para_nDFit%ndim
            IF (para_nDFit%nDvalB(i,iB) /=0) THEN
              IF (i1 /= 0) THEN
                i2 = i
              ELSE
                i1 = i
              END IF
            END IF
            IF (i2 /= 0) EXIT
          END DO
          !nb_couplings = count(para_nDFit%nDvalB(:,iB) /=0)

          IF (i1 /= 0 .AND. i2 == 0) THEN
            write(out_unitp,*) 'i1,ndDval(i1,.),B',i1,                  &
             '(',para_nDFit%nDvalB(i1,iB),')',para_nDFit%B(iB)*conv_col
          ELSE IF (i1 /= 0 .AND. i2 /= 0) THEN
            write(out_unitp,*) 'i1,ndDval(i1,.)...,B',i1,i2,            &
              '(',para_nDFit%nDvalB(i1,iB),para_nDFit%nDvalB(i2,iB),')',&
                para_nDFit%B(iB)*conv_col
          END IF

        END DO
        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"


        IF (para_nDFit%Analysis .AND.                                   &
            (para_nDFit%nDFitAna%Grid1D .OR.                       &
             para_nDFit%nDFitAna%Grid2D) ) THEN
          write(out_unitp,*) "======================================"
          write(out_unitp,*) "=== GRID PARAM FIT ==================="
          write(out_unitp,*) "======================================"
          Q(:) = para_nDFit%Q0(:)
          CALL sub_nDFunc_FROM_nDFit(val0,Q,para_nDFit)
          val_min_1D = val0
          ndim_coord = count(para_nDFit%nDFitAna%coord_list /= 0)

          IF (para_nDFit%nDFitAna%Grid1D) THEN
            Grid1D_file%name="Grid1D"
            CALL file_open(Grid1D_file,nio)
            ii = 0

            DO ic1=1,ndim_coord
              i1 = para_nDFit%nDFitAna%coord_list(ic1)
              val_min = huge(ONE)
              val_max = -huge(ONE)
              write(out_unitp,*) "====1D Grid",i1,"====================="
              Q(:) = para_nDFit%Q0(:)
              Q10  = para_nDFit%Q0(i1)
              Q(i1) = para_nDFit%nDFitAna%A(i1)
              DO iq1=1,para_nDFit%nDFitAna%nq(i1)
                CALL sub_nDFunc_FROM_nDFit(val_nDfit,Q,para_nDFit)
                IF (val_nDfit < val_min) THEN
                   val_min = val_nDfit
                   Q1_min  = Q(i1)
                END IF
                IF (val_nDfit > val_max) THEN
                   val_max = val_nDfit
                   Q1_max  = Q(i1)
                END IF
                write(nio,*) ii,i1,Q(i1),val_nDfit*conv_col
                Q(i1) = Q(i1) + para_nDFit%nDFitAna%Step(i1)
              END DO
              write(nio,*)
              write(nio,*)
              ii = ii + 1
              IF (val_min < val_min_1D) val_min_1D = val_min
              write(out_unitp,*) 'Q0   val0   ',i1,Q10,val0*conv_col
              write(out_unitp,*) 'Qmin val_min',i1,Q1_min,val_min*conv_col
              write(out_unitp,*) 'Qmax val_max',i1,Q1_max,val_max*conv_col
              IF (val_min < val0) write(out_unitp,*) 'WARNNING: val_min < val0'
              CALL flush_perso(out_unitp)
            END DO
            write(out_unitp,*) 'val_min_1D',val_min_1D*conv_col

            CALL file_close(Grid1D_file)
          END IF
          IF (para_nDFit%nDFitAna%Grid2D) THEN
            Grid2D_file%name="Grid2D"
            CALL file_open(Grid2D_file,nio)
            ii = 0

            DO ic1=1,ndim_coord
            DO ic2=ic1+1,ndim_coord
              i1 = para_nDFit%nDFitAna%coord_list(ic1)
              i2 = para_nDFit%nDFitAna%coord_list(ic2)

              write(out_unitp,*) "====2D Grid",i1,i2,"================"
              val_min = huge(ONE)
              val_max = -huge(ONE)
              Q(:) = para_nDFit%Q0(:)
              Q10  = para_nDFit%Q0(i1)
              Q20  = para_nDFit%Q0(i2)

              Q(i1) = para_nDFit%nDFitAna%A(i1)

              DO iq1=1,para_nDFit%nDFitAna%nq(i1)

                Q(i2) = para_nDFit%nDFitAna%A(i2)
                DO iq2=1,para_nDFit%nDFitAna%nq(i2)
                  CALL sub_nDFunc_FROM_nDFit(val_nDfit,Q,para_nDFit)
                  write(nio,*) ii,i1,i2,Q(i1),Q(i2),val_nDfit*conv_col

                  IF (val_nDfit < val_min) THEN
                     val_min = val_nDfit
                     Q1_min  = Q(i1)
                     Q2_min  = Q(i2)
                  END IF
                  IF (val_nDfit > val_max) THEN
                     val_max = val_nDfit
                     Q1_max  = Q(i1)
                     Q2_max  = Q(i2)
                  END IF
                  Q(i2) = Q(i2) + para_nDFit%nDFitAna%Step(i2)
                END DO
                Q(i1) = Q(i1) + para_nDFit%nDFitAna%Step(i1)
                write(nio,*)
              END DO
              write(out_unitp,*) 'Q0   val0   ',i1,i2,Q10,Q20,val0*conv_col
              write(out_unitp,*) 'Qmin val_min',i1,i2,Q1_min,Q2_min,val_min*conv_col
              write(out_unitp,*) 'Qmax val_max',i1,i2,Q1_max,Q2_max,val_max*conv_col
              IF (val_min < val_min_1D) write(out_unitp,*) 'WARNNING: val_min < val_min_1D'
              CALL flush_perso(out_unitp)

              write(nio,*)
              ii = ii + 1
            END DO
            END DO

            CALL file_close(Grid2D_file)
          END IF

          IF (para_nDFit%nDFitAna%Grid2D) THEN
            ii = 0

            DO ic1=1,ndim_coord
            DO ic2=ic1+1,ndim_coord
            DO ic3=ic2+1,ndim_coord
              i1 = para_nDFit%nDFitAna%coord_list(ic1)
              i2 = para_nDFit%nDFitAna%coord_list(ic2)
              i3 = para_nDFit%nDFitAna%coord_list(ic3)

              val_min = huge(ONE)
              val_max = -huge(ONE)
              Q(:) = para_nDFit%Q0(:)
              Q10  = para_nDFit%Q0(i1)
              Q20  = para_nDFit%Q0(i2)
              Q30  = para_nDFit%Q0(i2)

              Q(i1) = para_nDFit%nDFitAna%A(i1)

              DO iq1=1,para_nDFit%nDFitAna%nq(i1)

                Q(i2) = para_nDFit%nDFitAna%A(i2)
                DO iq2=1,para_nDFit%nDFitAna%nq(i2)

                  Q(i3) = para_nDFit%nDFitAna%A(i3)
                  DO iq3=1,para_nDFit%nDFitAna%nq(i3)

                    CALL sub_nDFunc_FROM_nDFit(val_nDfit,Q,para_nDFit)

                    IF (val_nDfit < val_min) THEN
                       val_min = val_nDfit
                       Q1_min  = Q(i1)
                       Q2_min  = Q(i2)
                       Q3_min  = Q(i3)
                    END IF
                    IF (val_nDfit > val_max) THEN
                       val_max = val_nDfit
                       Q1_max  = Q(i1)
                       Q2_max  = Q(i2)
                       Q3_max  = Q(i3)
                    END IF

                    Q(i3) = Q(i3) + para_nDFit%nDFitAna%Step(i3)
                  END DO
                  Q(i2) = Q(i2) + para_nDFit%nDFitAna%Step(i2)
                END DO
                Q(i1) = Q(i1) + para_nDFit%nDFitAna%Step(i1)
              END DO
              IF (val_min < val_min_1D) THEN
                write(out_unitp,*) "====3D Grid",i1,i2,i3,"============="
                write(out_unitp,*) 'Q0   val0   ',i1,i2,i3,Q10,Q20,Q30,val0*conv_col
                write(out_unitp,*) 'Qmin val_min',i1,i2,i3,Q1_min,Q2_min,Q3_min,val_min*conv_col
                write(out_unitp,*) 'Qmax val_max',i1,i2,i3,Q1_max,Q2_max,Q3_max,val_max*conv_col
                write(out_unitp,*) 'WARNNING: val_min < val_min_1D'
                CALL flush_perso(out_unitp)
              END IF

              ii = ii + 1
            END DO
            END DO
            END DO

          END IF


          write(out_unitp,*) "======================================"
          write(out_unitp,*) "======================================"
        END IF

      END SUBROUTINE Analysis_nDFit


      SUBROUTINE Read_FOR_nDFit1_TO_TnDFit2()
      USE mod_system
      IMPLICIT NONE

      TYPE (param_nDFit) :: para_nDFit1
      TYPE (param_nDFit) :: para_nDFit2


      integer                    :: i,iB,idum,nioFit,nb_Qact
      integer                    :: nb_coupling1,nb_coupling_act,nb_coupling_inact
      character (len=Name_len)   :: name_dum,name_int
      character (len=Line_len)   :: name_Fit1,name_Fit2
      integer, allocatable       :: nDinit(:)
      integer, allocatable       :: list_Qact(:)

      real (kind=Rkind)          :: Norm1
      logical                    :: keep

!----- For the namelist ----------------------------------------------
      integer           :: MinCoupling,MaxCoupling
      real (kind=Rkind) :: MinNorm,MaxNorm,conv_col,max_b,Weight_iGP
      logical           :: svd,ntyp_read,Keep_act_inact_couplings
      real (kind=Rkind) :: epsi,epsi_inter
      integer           :: ind_val,nb_val,nb_G,max_nb,MR_order
      integer           :: Col_FOR_WeightOFFit
      real (kind=Rkind) :: Scal_FOR_WeightOFFit

      namelist /nDFit/ MinNorm,MaxNorm,MinCoupling,MaxCoupling,         &
                       ind_val,nb_val,svd,ntyp_read,                    &
                       epsi,epsi_inter,max_nb,                          &
                       Col_FOR_WeightOFFit,Scal_FOR_WeightOFFit,MR_order,&
                       name_Fit1,name_Fit2,Keep_act_inact_couplings
      ! for the namelist
      !-----------------------------------------------------------------


!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'nDFit1_TO_TnDFit2'
!      logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------

      !write(out_unitp,*) 'SUBROUTINE ',trim(name_sub)
      !-- Read the new parameters
        MR_order             = -1  ! order of the multimode representation (-1: not use)
        MaxNorm              = FOUR
        MinNorm              = 0
        MaxCoupling          = 4
        MinCoupling          = 0
        svd                  = .TRUE.
        epsi                 = ONETENTH**10
        epsi_inter           = ONETENTH**5
        nb_val               = 1
        ind_val              = 1
        max_nb               = 10
        Col_FOR_WeightOFFit  = 0
        Scal_FOR_WeightOFFit = 200._Rkind
        ntyp_read            = .FALSE.
        name_Fit2            = ''
        name_Fit1            = "Param_FOR_Fit-col"
        Keep_act_inact_couplings = .TRUE.

        read(in_unitp,nDFit,IOSTAT=err_read)
        IF (err_read < 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' End-of-file or End-of-record'
          write(out_unitp,*) ' The namelist "nDFit" is probably absent'
          write(out_unitp,*) ' check your data!'
          write(out_unitp,*) ' ERROR in ',name_sub
          STOP
        ELSE IF (err_read > 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' Some parameter name of the namelist "nDFit" are probaly wrong'
          write(out_unitp,*) ' check your data!'
          write(out_unitp,nDFit)
          write(out_unitp,*) ' ERROR in ',name_sub
          STOP
        END IF
        IF (debug) write(out_unitp,nDFit)
!=====================================================================

      !-- First read the parameters of the grid
      CALL Write_int_IN_char(ind_val,name_int)
      para_nDFit1%name_Fit =                               &
           trim(adjustl(name_Fit1)) // trim(adjustl(name_int))
      para_nDFit1%Param_Fit_file%name = para_nDFit1%name_Fit
      write(out_unitp,*) 'name_fit_file: ',trim(para_nDFit1%Param_Fit_file%name)
      CALL ReadWrite_nDFitW(para_nDFit1,.TRUE.)

        write(out_unitp,*) "======================================"
        write(out_unitp,*) "=== EXPORT TRANSFORMED PARAM FIT ====="
        write(out_unitp,*) "======================================"

        para_nDFit2%ndim                 = para_nDFit1%ndim
        para_nDFit2%MR_order             = MR_order
        para_nDFit2%nb_val               = para_nDFit1%nb_val
        para_nDFit2%ind_val              = para_nDFit1%ind_val
        para_nDFit2%MinNorm              = MinNorm
        para_nDFit2%MaxNorm              = MaxNorm
        para_nDFit2%MinCoupling          = MinCoupling
        para_nDFit2%MaxCoupling          = MaxCoupling
        para_nDFit2%svd                  = para_nDFit1%svd
        para_nDFit2%max_nb               = max_nb
        para_nDFit2%epsi                 = epsi
        para_nDFit2%epsi_inter           = epsi_inter
        para_nDFit2%Col_FOR_WeightOFFit  = para_nDFit1%Col_FOR_WeightOFFit
        para_nDFit2%Scal_FOR_WeightOFFit = para_nDFit1%Scal_FOR_WeightOFFit


        CALL alloc_NParray(para_nDFit2%Q0,(/para_nDFit2%ndim/),         &
                          'para_nDFit2%Q0',name_sub)
        para_nDFit2%Q0(:) = para_nDFit1%Q0(:)
        CALL alloc_NParray(nDinit,(/para_nDFit2%ndim/),'nDinit',name_sub)
        nDinit(:) = 1
        CALL alloc_NParray(list_Qact,(/para_nDFit2%ndim/),'list_Qact',name_sub)
        list_Qact(:) = 0

        CALL alloc_NParray(para_nDFit2%nDweight,(/para_nDFit2%ndim/),   &
                          'para_nDFit2%nDweight',name_sub)
        CALL alloc_NParray(para_nDFit2%nDsize,(/para_nDFit2%ndim/),     &
                          'para_nDFit2%nDsize',name_sub)
        CALL alloc_NParray(para_nDFit2%ntyp,(/para_nDFit2%ndim/),       &
                          'para_nDFit2%ntyp',name_sub)
        para_nDFit2%ntyp(:) = para_nDFit1%ntyp(:)

        read(in_unitp,*) para_nDFit2%nDweight(:)
        read(in_unitp,*) para_nDFit2%nDsize(:)

        write(out_unitp,*) para_nDFit2%nDweight(:)
        write(out_unitp,*) para_nDFit2%nDsize(:)

        DO i=1,para_nDFit2%ndim
          CALL read_name_advNo(in_unitp,name_int,err_read)

          IF (len_trim(name_int) == 0) EXIT
          !write(out_unitp,*) 'i,err_io',i,err_io
          !write(out_unitp,*) 'i,name_int',i,name_int
          read(name_int,*) list_Qact(i)
          IF (err_read /= 0) EXIT ! end of the list

        END DO
        nb_Qact = count(list_Qact(:) > 0)
        write(out_unitp,*) 'list_Qact',list_Qact(:)


        CALL flush_perso(out_unitp)


        CALL init_nDindexPrim(para_nDFit2%nDindB,                       &
                              para_nDFit2%ndim,para_nDFit2%nDsize,      &
                              nDweight=para_nDFit2%nDweight,            &
                              type_OF_nDindex=0,                        &
                              nDinit=nDinit,                            &
                              MinNorm=para_nDFit2%MinNorm,              &
                              MaxNorm=para_nDFit2%MaxNorm,              &
                              MinCoupling=para_nDFit2%MinCoupling,      &
                              MaxCoupling=para_nDFit2%MaxCoupling)
        CALL sort_nDindex(para_nDFit2%nDindB)
        para_nDFit2%nDindB%Tab_nDval(:,:) = para_nDFit2%nDindB%Tab_nDval(:,:) - 1
        CALL Write_nDindex(para_nDFit2%nDindB)

        para_nDFit2%Param_Fit_file%name = name_Fit2
        para_nDFit2%name_Fit            = name_Fit2


        !First the new size
        para_nDFit2%nb_WB = 0
        DO iB=1,para_nDFit1%nb_WB

          keep = .TRUE.
          DO i=1,para_nDFit1%ndim
            keep = para_nDFit1%nDvalB(i,iB)+1 <= para_nDFit2%nDsize(i)
            IF (.NOT. keep) EXIT
          END DO
          nb_coupling1      = count(para_nDFit1%nDvalB(:,iB) > 0)
          nb_coupling_act   = count(para_nDFit1%nDvalB(list_Qact(1:nb_Qact),iB) > 0)
          nb_coupling_inact = nb_coupling1 - nb_coupling_act
          !write(out_unitp,*) 'iB,tab',iB,':',para_nDFit1%nDvalB(:,iB)
          !write(out_unitp,*) 'nb_coupling_act,nb_coupling_inact',nb_coupling_act,nb_coupling_inact

          !    act         inact    keep
          !1    0           0        +T
          !2    1           0        +T
          !3    0           1         F
          !4    1           1        +T or F  (T, if Keep_act_inact_couplings = T)

          IF (nb_coupling_inact == 0) THEN
            keep = keep .AND. .TRUE.             ! 1 and 2
          ELSE
            IF (nb_coupling_act == 0) THEN
              keep = .FALSE.                     ! 3
            ELSE
              IF (Keep_act_inact_couplings) THEN
                keep = keep .AND. .TRUE.         ! 4
              ELSE
                keep = .FALSE.                   ! 4
              END IF
            END IF
          END IF

          Norm1 = sum(real(para_nDFit1%nDvalB(:,iB),kind=Rkind)*para_nDFit1%nDweight)

          IF (Norm1 <= para_nDFit2%MaxNorm .AND.                        &
              Norm1 >= para_nDFit2%MinNorm .AND.                        &
              nb_coupling1 <= para_nDFit2%MaxCoupling .AND.             &
              nb_coupling1 >= para_nDFit2%MinCoupling .AND. keep) THEN

              para_nDFit2%nb_WB = para_nDFit2%nb_WB + 1

          ELSE
            para_nDFit1%B(iB) = ZERO
          END IF
        END DO

        CALL ReadWrite_nDFitW(para_nDFit2,.FALSE.,para_nDFit1%B,para_nDFit1%nDvalB)


        write(out_unitp,*) "======================================"
        write(out_unitp,*) "======================================"

        ! check the fit
        CALL Analysis_nDFitW(para_nDFit2,conv_ene=219475._Rkind)


        CALL dealloc_NParray(nDinit,   'nDinit',name_sub)
        CALL dealloc_NParray(list_Qact,'list_Qact',name_sub)

      END SUBROUTINE Read_FOR_nDFit1_TO_TnDFit2

      SUBROUTINE sub_nDFunc_FROM_nDFit(val_nDfit,Q,para_nDFit)
      USE mod_system
      IMPLICIT NONE

      real (kind=Rkind), intent(inout)       :: val_nDfit ! value of the function
      real (kind=Rkind), intent(in)          :: Q(:) ! value of the variables
      TYPE (param_nDFit), intent(inout)      :: para_nDFit



      real (kind=Rkind) :: val_nDfiti
      integer           :: i

!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'sub_nDFunc_FROM_nDFit'
!      logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------

!=====================================================================

!=====================================================================
!      initialization (only once)
!$OMP CRITICAL (nDFunc_FROM_nDFit_CRIT)
      IF (para_nDFit%nb_WB == 0 .AND. para_nDFit%nb_Fit == 0) THEN
        CALL ReadWrite_nDFitW(para_nDFit,ReadData=.TRUE.)
      END IF
!$OMP END CRITICAL (nDFunc_FROM_nDFit_CRIT)
!      END initialization
!=====================================================================

      IF (para_nDFit%nb_Fit == 0) THEN
        CALL sub_ONLYnDFunc_FROM_nDFit(val_nDfit,Q,para_nDFit)
      ELSE
        val_nDfit = ZERO
        DO i=1,para_nDFit%nb_Fit
          CALL sub_ONLYnDFunc_FROM_nDFit(val_nDfiti,Q,para_nDFit%Tab_para_nDFit(i))
          val_nDfit = val_nDfit +  val_nDfiti
        END DO
      END IF

      END SUBROUTINE sub_nDFunc_FROM_nDFit
      SUBROUTINE sub_ONLYnDFunc_FROM_nDFit(val_nDfit,Q,para_nDFit)
      USE mod_system
      IMPLICIT NONE

      real (kind=Rkind), intent(inout)       :: val_nDfit ! value of the function
      real (kind=Rkind), intent(in)          :: Q(:) ! value of the variables
      TYPE (param_nDFit), intent(inout)      :: para_nDFit

      integer           :: iB
      real (kind=Rkind) :: tQ(size(Q)) ! value of the variables
      real (kind=Rkind) :: t2Q(size(Q)) ! value of the variables
      real (kind=Rkind) :: a


!----- for debuging --------------------------------------------------
      integer :: err_read
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = 'sub_ONLYnDFunc_FROM_nDFit'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) " BEGINNING ",name_sub
        write(out_unitp,*) "Fit file name: ",trim(para_nDFit%Param_Fit_file%name)
        CALL Write_VecMat(Q,out_unitp,5,name_info='Q:')
      END IF


      IF (para_nDFit%nb_WB == 0 .OR. para_nDFit%nb_Fit /= 0) THEN
        write(out_unitp,*) " ERROR in ",name_sub
        write(out_unitp,*) "nb_WB=0 or nb_Fit /= 0 !",para_nDFit%nb_WB,para_nDFit%nb_Fit
        write(out_unitp,*) "It should not append"
        write(out_unitp,*) "CHECK the fortran!!"
        STOP
      END IF

      a=TWO
      tQ(1:para_nDFit%ndim)  = Q(1:para_nDFit%ndim)-para_nDFit%Q0(:)
      !t2Q(1:para_nDFit%ndim) = tanh(a*tQ(1:para_nDFit%ndim))/a


      val_nDfit = ZERO
      DO iB=1,para_nDFit%nb_WB

        !IF (count(para_nDFit%nDvalB(:,iB)  >  1) /= 0 .AND.             &
        !    count(para_nDFit%ntyp(:) /= 15) == 0) THEN
        !  val_nDfit = val_nDfit + para_nDFit%B(iB) *                    &
        !   nDFunct_WITH_tQ(t2Q(1:para_nDFit%ndim),iB,para_nDFit,para_nDFit%nDvalB(:,iB))
        !ELSE
          val_nDfit = val_nDfit + para_nDFit%B(iB) *                    &
           nDFunct_WITH_tQ(tQ(1:para_nDFit%ndim),iB,para_nDFit,para_nDFit%nDvalB(:,iB))
        !END IF
      END DO


      IF (debug) THEN
        CALL Write_VecMat(tQ,out_unitp,5,name_info='tQ:')
        write(out_unitp,*) 'Norm tQ',sqrt(dot_product(tQ,tQ))
        write(out_unitp,*) 'val_nDfit',val_nDfit
        write(out_unitp,*) " END ",name_sub
        CALL flush_perso(out_unitp)
      END IF

      END SUBROUTINE sub_ONLYnDFunc_FROM_nDFit
      FUNCTION nDFunct_WITH_Q(Q,inD,para_nDFit,nDvalB)
      USE mod_system
      IMPLICIT NONE


      ! type of the returned value
      real (kind=Rkind) :: nDFunct_WITH_Q

      real (kind=Rkind), intent(in)    :: Q(:) ! input variables
      integer, intent(in)              :: inD  ! index of the function
      TYPE (param_nDFit), intent(in)   :: para_nDFit ! definition of the function
      integer, intent(in), optional    :: nDvalB(:) ! input variables
      real (kind=Rkind)                :: V
      real (kind=Rkind)                :: tQ(size(Q)) ! input variables


      tQ(:) = Q(:)-para_nDFit%Q0(:)

      IF (present(nDvalB)) THEN
        nDFunct_WITH_Q = nDFunct_WITH_tQ(tQ,inD,para_nDFit,nDvalB)
      ELSE
        nDFunct_WITH_Q = nDFunct_WITH_tQ(tQ,inD,para_nDFit)
      END IF

      END FUNCTION nDFunct_WITH_Q
      FUNCTION nDFunct_WITH_tQ(tQ,inD,para_nDFit,nDvalB)
      USE mod_system
      IMPLICIT NONE


      ! type of the returned value
      real (kind=Rkind) :: nDFunct_WITH_tQ

      real (kind=Rkind), intent(in)    :: tQ(:) ! transformed Q, such tQ=Q-Q0
      integer, intent(in)              :: inD  ! index of the function
      TYPE (param_nDFit), intent(in)   :: para_nDFit ! definition of the function
      integer, intent(in), optional    :: nDvalB(:) ! input variables
      real (kind=Rkind) :: V
      real (kind=Rkind) :: poly_Hermite ! function
      real (kind=Rkind) :: Funct_1D     ! function

      integer :: i

      V = ONE
      IF (present(nDvalB)) THEN
        !write(out_unitp,*) 'nDvalB',nDvalB
        DO i=1,size(tQ)
          V = V * Funct_1D(tQ(i),nDvalB(i),para_nDFit%ntyp(i),0)
        END DO

      ELSE
        !write(out_unitp,*) 'nDvalB',para_nDFit%nDindB%Tab_nDval(:,inD)
        DO i=1,size(tQ)
          V = V * Funct_1D(tQ(i),para_nDFit%nDindB%Tab_nDval(i,inD),    &
                                                   para_nDFit%ntyp(i),0)
        END DO
      END IF
      nDFunct_WITH_tQ = V


      END FUNCTION nDFunct_WITH_tQ

      END MODULE mod_nDFit
