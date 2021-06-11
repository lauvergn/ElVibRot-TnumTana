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
MODULE mod_CRP
  USE mod_system
  IMPLICIT NONE

  TYPE CRP_Eckart_t
    real (kind=Rkind) :: V0 = 0.0156_Rkind   ! Baloitcha values
    real (kind=Rkind) :: L  = ONE            !  //
    real (kind=Rkind) :: m  = 1060._Rkind    !  //
  END TYPE CRP_Eckart_t
  TYPE CRP_Channel_AT_TS_t
    real (kind=Rkind)              :: EneTS = 0.0105_Rkind
    real (kind=Rkind)              :: w1    = 0.015625_Rkind
    real (kind=Rkind), allocatable :: w(:)
    integer                        :: option = 1
    integer                        :: nb_channels_added = 1
  END TYPE CRP_Channel_AT_TS_t
  TYPE param_CRP
    real (kind=Rkind) :: Ene    = ZERO            ! Total energy for CRP
    real (kind=Rkind) :: DEne   = ZERO            ! Energy increment for the CRP
    integer           :: nb_Ene = 1               ! Number of CRP calculation

    integer           :: iOp_CAP_Reactif = 3      ! Operator index of reactif CAP
    integer           :: iOp_CAP_Product = 4      ! Operator index of product CAP

    integer           :: iOp_Flux_Reactif = 5      ! Operator index of reactif flux
    integer           :: iOp_Flux_Product = 6      ! Operator index of product flux

    character (len=Name_len) :: CRP_Type            = 'lanczos'
    integer                  :: KS_max_it           = 100
    real (kind=Rkind)        :: KS_accuracy         = ONETENTH**5
    character (len=Name_len) :: LinSolv_Type        = 'MatInv'
    integer                  :: LinSolv_max_it      = 100
    real (kind=Rkind)        :: LinSolv_accuracy    = ONETENTH**7
    character (len=Name_len) :: Preconditioner_Type = 'Identity'
    logical                  :: FluxOp_test         = .FALSE.

    logical                  :: With_Eckart         = .FALSE.
    logical                  :: Read_Channel_AT_TS  = .FALSE.

    logical                  :: Build_MatOp         = .FALSE.


    TYPE (CRP_Eckart_t)        :: Eckart
    TYPE (CRP_Channel_AT_TS_t) :: Channel_AT_TS

  END TYPE param_CRP

CONTAINS

SUBROUTINE read_CRP(para_CRP,ny)
USE mod_system
USE mod_Constant
IMPLICIT NONE

  !----- variables pour la namelist analyse ----------------------------
  TYPE (param_CRP),     intent(inout)  :: para_CRP
  integer,              intent(in)     :: ny


  TYPE (REAL_WU) :: Ene,DEne
  integer        :: nb_Ene

  character (len=Name_len) :: CRP_Type            = 'Lanczos'
  integer                  :: KS_max_it           = 100
  real (kind=Rkind)        :: KS_accuracy         = ONETENTH**5
  character (len=Name_len) :: LinSolv_Type        = 'MatInv'
  integer                  :: LinSolv_max_it      = 100
  real (kind=Rkind)        :: LinSolv_accuracy    = ONETENTH**7
  character (len=Name_len) :: Preconditioner_Type = 'Identity'
  logical                  :: FluxOp_test         = .FALSE.

    TYPE (CRP_Eckart_t)    :: Eckart ! to be able to compar with Eckart CRP
    logical                :: With_Eckart

    logical                :: Read_Channel

  !----- for debuging --------------------------------------------------
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub = "read_CRP"
  !logical, parameter :: debug=.FALSE.
  logical, parameter :: debug=.TRUE.
  !-----------------------------------------------------------

  NAMELIST /CRP/Ene,DEne,nb_Ene,CRP_Type,                               &
                KS_max_it,KS_accuracy,                                  &
                LinSolv_Type,LinSolv_max_it,LinSolv_accuracy,           &
                Preconditioner_Type,FluxOp_test,                        &
                Eckart,With_Eckart,Read_Channel

  Ene                 = REAL_WU(ZERO,'cm-1','E')
  DEne                = REAL_WU(ZERO,'cm-1','E')
  nb_Ene              = 1

  CRP_Type            = 'lanczos'
  KS_max_it           = 100
  KS_accuracy         = ONETENTH**5
  LinSolv_Type        = 'QMR'
  LinSolv_max_it      = 100
  LinSolv_accuracy    = ONETENTH**7
  Preconditioner_Type = 'Diag'
  FluxOp_test         = .FALSE.
  With_Eckart         = .FALSE.
  Eckart              = CRP_Eckart_t(V0=0.0156_Rkind,L=ONE,m=1060._Rkind)
  Read_Channel        = .FALSE.

  read(in_unitp,CRP)
  write(out_unitp,CRP)

  para_CRP%With_Eckart = With_Eckart
  IF (With_Eckart) para_CRP%Eckart = Eckart

  CALL string_uppercase_TO_lowercase(CRP_Type)
  CALL string_uppercase_TO_lowercase(LinSolv_Type)
  CALL string_uppercase_TO_lowercase(Preconditioner_Type)

  IF (print_level > 0) write(out_unitp,CRP)
  write(out_unitp,*)

  para_CRP%Ene                  = convRWU_TO_R_WITH_WorkingUnit(Ene)
  para_CRP%DEne                 = convRWU_TO_R_WITH_WorkingUnit(DEne)
  para_CRP%nb_Ene               = nb_Ene
  para_CRP%CRP_Type             = CRP_Type
  para_CRP%KS_max_it            = KS_max_it
  para_CRP%KS_accuracy          = KS_accuracy
  para_CRP%LinSolv_Type         = LinSolv_Type
  para_CRP%LinSolv_max_it       = LinSolv_max_it
  para_CRP%LinSolv_accuracy     = LinSolv_accuracy
  para_CRP%Preconditioner_Type  = Preconditioner_Type
  para_CRP%FluxOp_test          = FluxOp_test
  para_CRP%Read_Channel_AT_TS   = Read_Channel

  para_CRP%Build_MatOp          = (CRP_type == 'withmat')                  .OR. &
                                  (CRP_type == 'withmat_flux')             .OR. &
                                  (CRP_type == 'withmatspectral')          .OR. &
         (CRP_type == 'lanczos'        .AND. LinSolv_Type == 'matinv')     .OR. &
         (CRP_type == 'lanczos'        .AND. LinSolv_Type == 'matlinsolv') .OR. &
         (CRP_type == 'lanczos_arpack' .AND. LinSolv_Type == 'matinv')     .OR. &
         (CRP_type == 'lanczos_arpack' .AND. LinSolv_Type == 'matlinsolv')

  IF (debug) write(out_unitp,*) 'E,DE,nb_E   : ',para_CRP%Ene,para_CRP%DEne,para_CRP%nb_Ene

  IF (Read_Channel) CALL Read_Channel_AT_TS(para_CRP%Channel_AT_TS,ny)

  write(out_unitp,*)
  CALL flush_perso(out_unitp)

END SUBROUTINE read_CRP
!================================================================
!     CRP
!================================================================
      SUBROUTINE sub_CRP(tab_Op,nb_Op,print_Op,para_CRP)
      USE mod_system
      USE mod_Coord_KEO
      USE mod_basis
      USE mod_Op
      IMPLICIT NONE


!----- Operator variables ----------------------------------------------
      integer,           intent(in)    :: nb_Op
      TYPE (param_Op),   intent(inout) :: tab_Op(nb_Op)
      logical,           intent(in)    :: print_Op

      TYPE (param_CRP),  intent(in)    :: para_CRP


      integer                           :: i
      real (kind=Rkind)                 :: Ene
      complex(kind=Rkind), allocatable  :: GuessVec(:)

!----- for debuging --------------------------------------------------
      integer   :: err
      logical, parameter :: debug=.FALSE.
!      logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'sub_CRP'
!-----------------------------------------------------------
      write(out_unitp,*) 'BEGINNING ',name_sub
      IF (debug) THEN
        write(out_unitp,*) 'shape tab_op',shape(tab_Op)
        CALL flush_perso(out_unitp)
        write(out_unitp,*)
      END IF
!-----------------------------------------------------------

      IF (para_CRP%FluxOp_test .AND. nb_Op < 6) Then
        write(out_unitp,*) ' The number of operator is wrong'
        write(out_unitp,*) ' nb_Op=',nb_Op
        write(out_unitp,*) ' For testing the flux, you MUST have 6 or more operators.'
        write(out_unitp,*)
        STOP ' ERROR in sub_CRP: wrong operator number'
      END IF
      IF (nb_Op < 4) THEN
        write(out_unitp,*) ' The number of operator is wrong'
        write(out_unitp,*) ' nb_Op=',nb_Op
        write(out_unitp,*) ' You MUST have 4 or more operators.'
        write(out_unitp,*) ' You HAVE to set-up: '
        write(out_unitp,*) '   - nb_scalar_Op=2 in the &minimum namelist'
        write(out_unitp,*) '  or '
        write(out_unitp,*) '   - nb_CAP=2 in the &active namelist'
        write(out_unitp,*)
        STOP ' ERROR in sub_CRP: wrong operator number'
      END IF

      IF (para_CRP%Build_MatOp) THEN
        CALL sub_MatOp(tab_Op(1),print_Op) ! H
        DO i=3,nb_Op ! for the CAP
          IF (i == para_CRP%iOp_CAP_Reactif .OR. i == para_CRP%iOp_CAP_Product)   &
              CALL sub_MatOp(tab_Op(i),print_Op)
        END DO
      END IF
      IF (tab_Op(1)%Partial_MatOp) STOP 'STOP the Matrices are incomplete'


      SELECT CASE (para_CRP%CRP_type)
      CASE ('withmat') ! old one
        CALL sub_CRP_BasisRep_WithMat(tab_Op,nb_Op,print_Op,para_CRP)

      CASE ('withmat_flux')
        CALL sub_CRP_BasisRep_WithMat_flux(tab_Op,nb_Op,print_Op,para_CRP)

      CASE ('withmatspectral') ! old one
        CALL sub_CRP_BasisRep_WithMatSpectral(tab_Op,nb_Op,print_Op,para_CRP)
        !CALL sub_CRP_BasisRep_WithMatSpectral_old(tab_Op,nb_Op,print_Op,para_CRP)

      CASE ('lanczos') ! lanczos (Lucien Dupuy)

!$OMP   PARALLEL &
!$OMP   DEFAULT(NONE) &
!$OMP   SHARED(para_CRP,tab_Op,nb_Op) &
!$OMP   PRIVATE(i,Ene,GuessVec) &
!$OMP   NUM_THREADS(CRP_maxth)

        CALL alloc_NParray(GuessVec, [tab_Op(1)%nb_tot], 'GuessVec', name_sub)
        GuessVec(:) = ZERO

!$OMP   DO SCHEDULE(STATIC)

        DO i = 0, para_CRP%nb_Ene-1
          Ene = para_CRP%Ene+real(i,kind=Rkind)*para_CRP%DEne
          CALL calc_crp_P_lanczos(tab_Op, nb_Op,para_CRP,Ene,GuessVec)
        END DO

!$OMP   END DO

        CALL dealloc_NParray(GuessVec, 'GuessVec', name_sub)

!$OMP   END PARALLEL

      CASE ('lanczos_arpack') ! lanczos (Lucien Dupuy)

        ! OpenMP ne fonctionne pas avec Arpack
        DO i = 0, para_CRP%nb_Ene-1
          Ene = para_CRP%Ene+real(i,kind=Rkind)*para_CRP%DEne
          CALL calc_crp_IRL(tab_Op, nb_Op,para_CRP,Ene)
        END DO

      END SELECT

!----------------------------------------------------------
      IF (debug) THEN
      END IF
      write(out_unitp,*) 'END ',name_sub
!----------------------------------------------------------

      end subroutine sub_CRP
      SUBROUTINE sub_CRP_BasisRep_WithMat(tab_Op,nb_Op,print_Op,para_CRP)

      USE mod_system
      USE mod_Constant
      USE mod_Coord_KEO
      USE mod_basis
      USE mod_Op
      IMPLICIT NONE


!----- Operator variables ----------------------------------------------
      integer,            intent(in)      :: nb_Op
      TYPE (param_Op)                     :: tab_Op(nb_Op)
      logical,            intent(in)      :: print_Op
      TYPE (param_CRP),   intent(in)      :: para_CRP

      !real (kind=Rkind) :: CRP_Ene,CRP_DEne
      !integer           :: nb_CRP_Ene

!---- variable for the Z-matrix ----------------------------------------
      TYPE (CoordType), pointer  :: mole
      TYPE (Tnum), pointer       :: para_Tnum


!----- working variables -----------------------------
      integer       ::    i,j,k,ie
      real (kind=Rkind), allocatable :: EneH(:),Vec(:,:) ! for debuging


      complex (kind=Rkind), allocatable :: G(:,:)
      complex (kind=Rkind), allocatable :: Ginv(:,:)
      complex (kind=Rkind), allocatable :: gGgG(:,:)
      complex (kind=Rkind) :: CRP
      real (kind=Rkind) :: Ene
      TYPE(REAL_WU)     :: RWU_E



!----- for debuging --------------------------------------------------
      integer   :: err
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'sub_CRP_BasisRep_WithMat'
!-----------------------------------------------------------
      mole       => tab_Op(1)%mole
      para_Tnum  => tab_Op(1)%para_Tnum

      write(out_unitp,*) 'BEGINNING ',name_sub
      IF (debug) THEN
        write(out_unitp,*) 'shape tab_op',shape(tab_Op)
        CALL flush_perso(out_unitp)
        write(out_unitp,*)
      END IF
!-----------------------------------------------------------

      write(out_unitp,*) 'nb_tot of H',tab_Op(1)%nb_tot
      CALL flush_perso(out_unitp)

      IF (debug) THEN
        write(out_unitp,*) 'shape H',shape(tab_Op(1)%Rmat)
        CALL alloc_NParray(Vec,shape(tab_Op(1)%Rmat),'Vec',name_sub)
        CALL alloc_NParray(EneH,shape(tab_Op(1)%Rmat(:,1)),'EneH',name_sub)

        CALL sub_diago_H(tab_Op(1)%Rmat,EneH,Vec,tab_Op(1)%nb_tot,.TRUE.)
        write(out_unitp,*) 'Ene (ua)',EneH(1:min(10,tab_Op(1)%nb_tot))

        CALL dealloc_NParray(Vec,'Vec',name_sub)
        CALL dealloc_NParray(EneH,'EneH',name_sub)
      END IF


      CALL alloc_NParray(Ginv,shape(tab_Op(1)%Rmat),'Ginv',name_sub)
      CALL alloc_NParray(G,shape(tab_Op(1)%Rmat),'G',name_sub)
      CALL alloc_NParray(gGgG,shape(tab_Op(1)%Rmat),'gGgG',name_sub)


      write(out_unitp,*) 'Ginv calc'
      Ginv(:,:) = -tab_Op(1)%Rmat + EYE*HALF * (tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat+ &
                                                tab_Op(para_CRP%iOp_CAP_Product)%Rmat)

      DO i=1,tab_Op(1)%nb_tot
        Ginv(i,i) = Ginv(i,i) + para_CRP%Ene-para_CRP%DEne
      END DO
      Ene = para_CRP%Ene-para_CRP%DEne

      DO ie=1,para_CRP%nb_Ene

        DO i=1,tab_Op(1)%nb_tot
          Ginv(i,i) = Ginv(i,i) + para_CRP%DEne
        END DO
        Ene = Ene + para_CRP%DEne

! write(6,*) 'Ginv'
!        DO i=1,tab_Op(1)%nb_tot
!        DO j=1,tab_Op(1)%nb_tot
!           write(6,*) i,j,Ginv(j,i)
!        END DO
!        END DO

        CALL inv_m1_TO_m2_cplx(Ginv,G,tab_Op(1)%nb_tot,0,ZERO)
        !Ginv = matmul(Ginv,G)
        !DO i=1,tab_Op(1)%nb_tot
        !  Ginv(i,i) = Ginv(i,i) - CONE
        !END DO
        !write(out_unitp,*) 'id diff ?',maxval(abs(Ginv))

! write(6,*) 'G'
!        DO i=1,tab_Op(1)%nb_tot
!        DO j=1,tab_Op(1)%nb_tot
!           write(6,*) i,j,G(j,i)
!        END DO
!        END DO

        gGgG(:,:) = matmul(tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat,               &
           matmul(G,matmul(tab_Op(para_CRP%iOp_CAP_Product)%Rmat,conjg(G))))

! write(6,*) 'gammaR.G.gammaP.G*'
!        DO i=1,tab_Op(1)%nb_tot
!        DO j=1,tab_Op(1)%nb_tot
!           write(6,*) i,j,gGgG(j,i)
!        END DO
!        END DO
!STOP

        CRP = ZERO
        DO i=1,tab_Op(1)%nb_tot
          CRP = CRP + gGgG(i,i)
        END DO

        RWU_E  = REAL_WU(Ene,'au','E')

        if (para_CRP%With_Eckart) then
          write(out_unitp,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                            real(CRP,kind=Rkind),aimag(CRP),CRP_Eckart(Ene,para_CRP%Eckart),&
                            real(CRP,kind=Rkind)-CRP_Eckart(Ene,para_CRP%Eckart)
        else
          write(out_unitp,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                            real(CRP,kind=Rkind),aimag(CRP)
        end if
        CALL flush_perso(out_unitp)

      END DO


      CALL dealloc_NParray(Ginv,'Ginv',name_sub)
      CALL dealloc_NParray(gGgG,'gGgG',name_sub)
      CALL dealloc_NParray(G,'G',name_sub)
!----------------------------------------------------------
      IF (debug) THEN
      END IF
      write(out_unitp,*) 'END ',name_sub
      CALL flush_perso(out_unitp)
!----------------------------------------------------------

END SUBROUTINE sub_CRP_BasisRep_WithMat
SUBROUTINE sub_CRP_BasisRep_WithMatSpectral(tab_Op,nb_Op,print_Op,para_CRP)

      USE mod_system
      USE mod_Constant
      USE mod_Coord_KEO
      USE mod_basis
      USE mod_Op
      IMPLICIT NONE


!----- Operator variables ----------------------------------------------
      integer,            intent(in)      :: nb_Op
      TYPE (param_Op)                     :: tab_Op(nb_Op)
      logical,            intent(in)      :: print_Op
      TYPE (param_CRP),   intent(in)      :: para_CRP

      !real (kind=Rkind) :: CRP_Ene,CRP_DEne
      !integer           :: nb_CRP_Ene

!---- variable for the Z-matrix ----------------------------------------
      TYPE (CoordType), pointer  :: mole
      TYPE (Tnum), pointer       :: para_Tnum


!----- working variables -----------------------------
      integer       ::    i,j,k,ie

      complex (kind=Rkind), allocatable :: G(:,:)
      complex (kind=Rkind), allocatable :: Ginv(:,:)
      complex (kind=Rkind), allocatable :: VecPGinv(:,:)
      complex (kind=Rkind), allocatable :: ValPGinv(:)
      complex (kind=Rkind), allocatable :: ValPG(:)
      complex (kind=Rkind), allocatable :: Vec1(:)
      complex (kind=Rkind) :: CRP
      real (kind=Rkind) :: Ene
      TYPE(REAL_WU)     :: RWU_E



!----- for debuging --------------------------------------------------
      integer   :: err
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'sub_CRP_BasisRep_WithMatSpectral'
!-----------------------------------------------------------
      mole       => tab_Op(1)%mole
      para_Tnum  => tab_Op(1)%para_Tnum

      write(out_unitp,*) 'BEGINNING ',name_sub
      IF (debug) THEN
        write(out_unitp,*) 'shape tab_op',shape(tab_Op)
        CALL flush_perso(out_unitp)
        write(out_unitp,*)
      END IF
!-----------------------------------------------------------

      write(out_unitp,*) 'nb_tot of H',tab_Op(1)%nb_tot
      CALL flush_perso(out_unitp)

      CALL alloc_NParray(Ginv,    shape(tab_Op(1)%Rmat),'Ginv',    name_sub)
      CALL alloc_NParray(VecPGinv,shape(tab_Op(1)%Rmat),'VecPGinv',name_sub)
      CALL alloc_NParray(ValPGinv,[tab_Op(1)%nb_tot],   'ValPGinv',name_sub)
      CALL alloc_NParray(ValPG,   [tab_Op(1)%nb_tot],   'ValPG',   name_sub)

      write(out_unitp,*) 'Ginv calc' ; CALL flush_perso(out_unitp)
      Ginv(:,:) = -tab_Op(1)%Rmat + EYE*HALF * (tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat+ &
                                                tab_Op(para_CRP%iOp_CAP_Product)%Rmat)

      write(out_unitp,*) 'Ginv diago' ; CALL flush_perso(out_unitp)
      CALL sub_diago_CH(Ginv,ValPGinv,VecPGinv,tab_Op(1)%nb_tot)
      write(out_unitp,*) 'Ginv diago: done' ; CALL flush_perso(out_unitp)

      CALL dealloc_NParray(Ginv,'Ginv',name_sub)


      CALL alloc_NParray(Vec1,[tab_Op(1)%nb_tot],'Vec1',name_sub)


      Ene = para_CRP%Ene-para_CRP%DEne

      DO ie=1,para_CRP%nb_Ene

        Ene = Ene + para_CRP%DEne

        ValPG(:) = ONE/(ValPGinv+Ene)

        CRP = ZERO
        DO i=1,tab_Op(1)%nb_tot

          Vec1 = conjg(matmul(VecPGinv,ValPG * VecPGinv(i,:)))

          Vec1 = matmul(tab_Op(para_CRP%iOp_CAP_Product)%Rmat,Vec1)
          Vec1 = matmul(Vec1,VecPGinv) ! equvivalent to matmul(transpose(VecPGinv),Vec1)

          Vec1 = matmul(VecPGinv,ValPG * Vec1)
          CRP = CRP + sum(tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat(i,:)*Vec1) ! Cannot use dot_product because of conjg of CAP_Reactif

        END DO

        RWU_E  = REAL_WU(Ene,'au','E')

        if (para_CRP%With_Eckart) then
          write(out_unitp,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                            real(CRP,kind=Rkind),aimag(CRP),CRP_Eckart(Ene,para_CRP%Eckart),&
                            real(CRP,kind=Rkind)-CRP_Eckart(Ene,para_CRP%Eckart)
        else
          write(out_unitp,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                            real(CRP,kind=Rkind),aimag(CRP)
        end if
        CALL flush_perso(out_unitp)

      END DO

      CALL dealloc_NParray(VecPGinv,'VecPGinv',name_sub)
      CALL dealloc_NParray(ValPGinv,'ValPGinv',name_sub)
      CALL dealloc_NParray(ValPG,   'ValPG',   name_sub)
      CALL dealloc_NParray(Vec1,   'Vec1',   name_sub)

!----------------------------------------------------------
      IF (debug) THEN
      END IF
      write(out_unitp,*) 'END ',name_sub
      CALL flush_perso(out_unitp)
!----------------------------------------------------------

END SUBROUTINE sub_CRP_BasisRep_WithMatSpectral
SUBROUTINE sub_CRP_BasisRep_WithMatSpectral_old(tab_Op,nb_Op,print_Op,para_CRP)

      USE mod_system
      USE mod_Constant
      USE mod_Coord_KEO
      USE mod_basis
      USE mod_Op
      IMPLICIT NONE


!----- Operator variables ----------------------------------------------
      integer,            intent(in)      :: nb_Op
      TYPE (param_Op)                     :: tab_Op(nb_Op)
      logical,            intent(in)      :: print_Op
      TYPE (param_CRP),   intent(in)      :: para_CRP

      !real (kind=Rkind) :: CRP_Ene,CRP_DEne
      !integer           :: nb_CRP_Ene

!---- variable for the Z-matrix ----------------------------------------
      TYPE (CoordType), pointer  :: mole
      TYPE (Tnum), pointer       :: para_Tnum


!----- working variables -----------------------------
      integer       ::    i,j,k,ie

      complex (kind=Rkind), allocatable :: G(:,:)
      complex (kind=Rkind), allocatable :: Ginv(:,:)
      complex (kind=Rkind), allocatable :: VecPGinv(:,:)
      complex (kind=Rkind), allocatable :: ValPGinv(:)
      complex (kind=Rkind), allocatable :: ValPG(:)
      complex (kind=Rkind), allocatable :: Vec1(:),Vec2(:)
      complex (kind=Rkind), allocatable :: Mat1(:,:),Mat2(:,:)

      complex (kind=Rkind), allocatable :: CAP_reactif(:,:),CAP_product(:,:)


      complex (kind=Rkind), allocatable :: gGgG(:,:)
      complex (kind=Rkind) :: CRP
      real (kind=Rkind) :: Ene
      TYPE(REAL_WU)     :: RWU_E



!----- for debuging --------------------------------------------------
      integer   :: err
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'sub_CRP_BasisRep_WithMatSpectral_old'
!-----------------------------------------------------------
      mole       => tab_Op(1)%mole
      para_Tnum  => tab_Op(1)%para_Tnum

      write(out_unitp,*) 'BEGINNING ',name_sub
      IF (debug) THEN
        write(out_unitp,*) 'shape tab_op',shape(tab_Op)
        CALL flush_perso(out_unitp)
        write(out_unitp,*)
      END IF
!-----------------------------------------------------------

      write(out_unitp,*) 'nb_tot of H',tab_Op(1)%nb_tot
      CALL flush_perso(out_unitp)

      CALL alloc_NParray(Ginv,    shape(tab_Op(1)%Rmat),'Ginv',    name_sub)
      CALL alloc_NParray(VecPGinv,shape(tab_Op(1)%Rmat),'VecPGinv',name_sub)
      CALL alloc_NParray(ValPGinv,[tab_Op(1)%nb_tot],   'ValPGinv',name_sub)
      CALL alloc_NParray(ValPG,   [tab_Op(1)%nb_tot],   'ValPG',   name_sub)


      CALL alloc_NParray(gGgG,shape(tab_Op(1)%Rmat),'gGgG',name_sub)


      write(out_unitp,*) 'Ginv calc' ; CALL flush_perso(out_unitp)
      Ginv(:,:) = -tab_Op(1)%Rmat + EYE*HALF * (tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat+ &
                                                tab_Op(para_CRP%iOp_CAP_Product)%Rmat)

      write(out_unitp,*) 'Ginv diago' ; CALL flush_perso(out_unitp)
      CALL sub_diago_CH(Ginv,ValPGinv,VecPGinv,tab_Op(1)%nb_tot)
      write(out_unitp,*) 'Ginv diago: done' ; CALL flush_perso(out_unitp)

!tesy diago
      Ginv(:,:) = -tab_Op(1)%Rmat + EYE*HALF * (tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat+ &
                                                 tab_Op(para_CRP%iOp_CAP_Product)%Rmat)
      Ginv = matmul(transpose(VecPGinv),matmul(Ginv,VecPGinv))
      !CALL Write_Mat(Ginv,6,7)
      DO i=1,tab_Op(1)%nb_tot
           Ginv(i,i) = Ginv(i,i)-ValPGinv(i)
      END DO
      write(6,*) 'diago?',maxval(abs(Ginv))

      Ginv = ZERO
      DO i=1,tab_Op(1)%nb_tot
           Ginv(i,i) = ValPGinv(i)
      END DO
      Ginv = matmul(VecPGinv,matmul(Ginv,transpose(VecPGinv)))
      Ginv(:,:) = Ginv - (-tab_Op(1)%Rmat + EYE*HALF * (tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat+ &
                                          tab_Op(para_CRP%iOp_CAP_Product)%Rmat))
      write(6,*) 'diago?',maxval(abs(Ginv))
      stop

      CALL dealloc_NParray(Ginv,'Ginv',name_sub)

      CALL alloc_NParray(CAP_Reactif,shape(tab_Op(1)%Rmat),'CAP_Reactif',name_sub)
      CAP_Reactif = matmul(transpose(VecPGinv),matmul(tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat,VecPGinv))
      !write(out_unitp,*) 'CAP_Reactif'
      !CALL Write_Mat(CAP_Reactif,6,5)

      CALL alloc_NParray(CAP_Product,shape(tab_Op(1)%Rmat),'CAP_Product',name_sub)
      CAP_Product = matmul(transpose(VecPGinv),matmul(tab_Op(para_CRP%iOp_CAP_Product)%Rmat,VecPGinv))
      !write(out_unitp,*) 'CAP_Product'
      !CALL Write_Mat(CAP_Product,6,5)

      CALL alloc_NParray(Mat1,shape(tab_Op(1)%Rmat),'Mat1',name_sub)
      CALL alloc_NParray(Mat2,shape(tab_Op(1)%Rmat),'Mat2',name_sub)

      !CAP_Reactif = matmul(VecPGinv,matmul(CAP_Reactif,transpose(VecPGinv))) ! v1
      CAP_Reactif = matmul(VecPGinv,CAP_Reactif) ! v2

      CAP_Product = matmul(VecPGinv,matmul(CAP_Product,transpose(VecPGinv))) !v1,v2


      ValPGinv(:) = ValPGinv(:) + para_CRP%Ene-para_CRP%DEne
      Ene = para_CRP%Ene-para_CRP%DEne

      DO ie=1,para_CRP%nb_Ene


        ValPGinv(:) = ValPGinv(:) + para_CRP%DEne
        Ene = Ene + para_CRP%DEne

        ValPG(:) = ONE/ValPGinv

        Mat1 = CZERO
        Mat2 = CZERO
        DO i=1,tab_Op(1)%nb_tot
          Mat1(i,i) = conjg(ValPG(i))
          Mat2(i,i) = ValPG(i)
        END DO

        Mat1 = matmul(conjg(VecPGinv),matmul(Mat1,transpose(conjg(VecPGinv)))) ! v1,v2
        !Mat2 = matmul(VecPGinv,matmul(Mat2,transpose(VecPGinv)))                !v1
        Mat2 = matmul(Mat2,transpose(VecPGinv))                !v2



        gGgG(:,:) = matmul(matmul(CAP_Reactif,Mat2),matmul(CAP_Product,Mat1))

        !gGgG(:,:) = matmul(tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat,               &
        !   matmul(G,matmul(tab_Op(para_CRP%iOp_CAP_Product)%Rmat,conjg(G))))


        CRP = ZERO
        DO i=1,tab_Op(1)%nb_tot
!           Vec1 = ValPG(:) * CAP_Product(:,i) * conjg(ValPG(i))
! CRP = CRP + dot_product(CAP_Reactif(:,i),Vec1)
          CRP = CRP + gGgG(i,i)
        END DO

        RWU_E  = REAL_WU(Ene,'au','E')

        if (para_CRP%With_Eckart) then
          write(out_unitp,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                            real(CRP,kind=Rkind),aimag(CRP),CRP_Eckart(Ene,para_CRP%Eckart),&
                            real(CRP,kind=Rkind)-CRP_Eckart(Ene,para_CRP%Eckart)
        else
          write(out_unitp,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                            real(CRP,kind=Rkind),aimag(CRP)
        end if
        CALL flush_perso(out_unitp)

      END DO

      CALL dealloc_NParray(CAP_Reactif,'CAP_Reactif',name_sub)
      CALL dealloc_NParray(CAP_Product,'CAP_Product',name_sub)
      CALL dealloc_NParray(VecPGinv,'VecPGinv',name_sub)
      CALL dealloc_NParray(ValPGinv,'ValPGinv',name_sub)
      CALL dealloc_NParray(ValPG,   'ValPG',   name_sub)

      CALL dealloc_NParray(gGgG,'gGgG',name_sub)
!----------------------------------------------------------
      IF (debug) THEN
      END IF
      write(out_unitp,*) 'END ',name_sub
      CALL flush_perso(out_unitp)
!----------------------------------------------------------

END SUBROUTINE sub_CRP_BasisRep_WithMatSpectral_old
SUBROUTINE sub_CRP_BasisRep_WithMat_flux(tab_Op,nb_Op,print_Op,para_CRP)

      USE mod_system
      USE mod_Constant
      USE mod_Coord_KEO
      USE mod_basis
      USE mod_Op
      IMPLICIT NONE


!----- Operator variables ----------------------------------------------
      integer,            intent(in)      :: nb_Op
      TYPE (param_Op)                     :: tab_Op(nb_Op)
      logical,            intent(in)      :: print_Op
      TYPE (param_CRP),   intent(in)      :: para_CRP

      !real (kind=Rkind) :: CRP_Ene,CRP_DEne
      !integer           :: nb_CRP_Ene

!---- variable for the Z-matrix ----------------------------------------
      TYPE (CoordType), pointer  :: mole
      TYPE (Tnum), pointer       :: para_Tnum


!----- working variables -----------------------------
      integer       ::    i,k,ie


      complex (kind=Rkind), allocatable :: G(:,:)
      complex (kind=Rkind), allocatable :: gGgG(:,:)
      complex (kind=Rkind) :: CRP
      real (kind=Rkind) :: Ene
      TYPE(REAL_WU)     :: RWU_E

      real (kind=Rkind), allocatable :: mEYE_FluxOpReactif_mat(:,:)
      real (kind=Rkind), allocatable :: mEYE_FluxOpProduct_mat(:,:)



!----- for debuging --------------------------------------------------
      integer   :: err
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'sub_CRP_BasisRep_WithMat_flux'
!-----------------------------------------------------------
      mole       => tab_Op(1)%mole
      para_Tnum  => tab_Op(1)%para_Tnum

      write(out_unitp,*) 'BEGINNING ',name_sub
      IF (debug) THEN
        write(out_unitp,*) 'shape tab_op',shape(tab_Op)
        CALL flush_perso(out_unitp)
        write(out_unitp,*)
      END IF
!-----------------------------------------------------------

      write(out_unitp,*) 'nb_tot of H',tab_Op(1)%nb_tot


      CALL alloc_NParray(G,shape(tab_Op(1)%Rmat),'G',name_sub)
      CALL alloc_NParray(gGgG,shape(tab_Op(1)%Rmat),'gGgG',name_sub)

      CALL alloc_NParray(mEYE_FluxOpProduct_mat,shape(tab_Op(1)%Rmat), &
                        'mEYE_FluxOpProduct_mat',name_sub)
      CALL alloc_NParray(mEYE_FluxOpReactif_mat,shape(tab_Op(1)%Rmat), &
                        'mEYE_FluxOpReactif_mat',name_sub)

      DO i=3,nb_Op ! for the flux and the cap
        ! Here,  we don't calculate the flux operator, but -i.FluxOp = [H,HStep]
        ! Because the corresponding matrix is real.
        IF (i == para_CRP%iOp_Flux_Reactif) Then
          write(out_unitp,*) 'Op name: ',tab_Op(i)%name_Op
          CALL sub_MatOp(tab_Op(i),print_Op)
          CALL FluxOp_Mat_v0(tab_Op(1),tab_Op(i),mEYE_FluxOpReactif_mat)
          STOP
        END IF
        IF (i == para_CRP%iOp_Flux_Product) Then
          write(out_unitp,*) 'Op name: ',tab_Op(i)%name_Op
          CALL sub_MatOp(tab_Op(i),print_Op)
          CALL FluxOp_Mat(tab_Op(1),tab_Op(i),mEYE_FluxOpProduct_mat)
        END IF

      END DO

      Ene = para_CRP%Ene
      DO ie=1,para_CRP%nb_Ene

        CALL G_Mat(tab_Op(1),tab_Op(para_CRP%iOp_CAP_Reactif),                  &
                             tab_Op(para_CRP%iOp_CAP_Product),Ene,G)

        gGgG(:,:) = matmul(mEYE_FluxOpReactif_mat,                              &
                    matmul(G,matmul(mEYE_FluxOpProduct_mat,conjg(G))))

        CRP = ZERO
        DO i=1,tab_Op(1)%nb_tot
          CRP = CRP + gGgG(i,i)
        END DO

        RWU_E  = REAL_WU(Ene,'au','E')

        if (para_CRP%With_Eckart) then
          write(out_unitp,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                            real(CRP,kind=Rkind),aimag(CRP),CRP_Eckart(Ene,para_CRP%Eckart),&
                            real(CRP,kind=Rkind)-CRP_Eckart(Ene,para_CRP%Eckart)

        else
          write(out_unitp,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                            real(CRP,kind=Rkind),aimag(CRP)
        end if
        CALL flush_perso(out_unitp)

        Ene = Ene + para_CRP%DEne


      END DO


      CALL dealloc_NParray(gGgG,'gGgG',name_sub)
      CALL dealloc_NParray(G,'G',name_sub)

      CALL dealloc_NParray(mEYE_FluxOpProduct_mat,'mEYE_FluxOpProduct_mat',name_sub)
      CALL dealloc_NParray(mEYE_FluxOpReactif_mat,'mEYE_FluxOpReactif_mat',name_sub)


!----------------------------------------------------------
      IF (debug) THEN
      END IF
      write(out_unitp,*) 'END ',name_sub
      CALL flush_perso(out_unitp)
!----------------------------------------------------------

END SUBROUTINE sub_CRP_BasisRep_WithMat_flux
SUBROUTINE calc_crp_p_lanczos(tab_Op,nb_Op,para_CRP,Ene,GuessVec)

      USE mod_Constant
      USE mod_Op
      implicit none

!----- Operator variables ----------------------------------------------
      integer,             intent(in)      :: nb_Op
      TYPE (param_Op),     intent(inout)   :: tab_Op(nb_Op)

      TYPE (param_CRP),    intent(in)      :: para_CRP
      real(kind=Rkind),    intent(in)      :: Ene
      complex(kind=Rkind), intent(inout)   :: GuessVec(:)

!----- working variables -----------------------------
      TYPE(REAL_WU)     :: RWU_E

      ! Calculate the inverse matrix explicitly? For debuging. It is equivalent to sub_CRP_BasisRep_WithMat
      logical, parameter :: Inv = .FALSE.

      ! Size of Hamiltonian matrix
      integer :: ncooked
      ! Loop integers
      integer :: nks, mks, i, j
      ! Vectors for Pmult
      complex(kind=Rkind),  allocatable :: Krylov_vectors(:,:),h(:,:)
      real (kind=Rkind),    allocatable :: Eigvals(:)
      complex(kind=Rkind),  allocatable :: EVec(:,:)

      complex(kind=Rkind) y, len
      real(kind=Rkind) ranr, rani

      ! Operator variables ----------------------------------------------
      logical           :: print_Op


      integer           ::    k,ie
      real (kind=Rkind) :: CRP, oldcrp,crp2,DeltaCRP
      complex (kind=Rkind), allocatable :: G(:,:)
      complex (kind=Rkind), allocatable :: Ginv(:,:)
      complex (kind=Rkind), allocatable :: gGgG(:,:)
      complex (kind=Rkind), allocatable :: M1(:)

       ! temporary variables for the LU decomposition
       integer,             allocatable :: indx(:)
       complex(kind=Rkind), allocatable :: trav(:)
       complex(kind=Rkind)              :: d

       integer :: lwork,ierr
       real (kind=Rkind), ALLOCATABLE :: work(:)

       TYPE (param_time) :: CRP_Time
       real(kind=Rkind)  :: RealTime

!----- for debuging --------------------------------------------------
      integer   :: err
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'calc_crp_p_lanczos'
!-----------------------------------------------------------


      ncooked   = tab_Op(1)%nb_tot


    ! If need be, generate explicit matrix representation of operators
    IF ( Inv ) THEN

!$OMP SINGLE
      CALL sub_MatOp(tab_Op(1),print_Op=.FALSE.) ! H
      DO i=3,nb_Op ! for the CAP
        IF (i == para_CRP%iOp_CAP_Reactif .OR. i == para_CRP%iOp_CAP_Product)   &
           CALL sub_MatOp(tab_Op(i),print_Op=.FALSE.)
      END DO
!$OMP END SINGLE

      CALL alloc_NParray(Ginv,shape(tab_Op(1)%Rmat),'Ginv',name_sub)
      CALL alloc_NParray(G,   shape(tab_Op(1)%Rmat),'G',   name_sub)
      CALL alloc_NParray(gGgG,shape(tab_Op(1)%Rmat),'gGgG',name_sub)

      Ginv(:,:) = -tab_Op(1)%Rmat + EYE*HALF * (tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat+ &
                                                tab_Op(para_CRP%iOp_CAP_Product)%Rmat)

      DO i=1,tab_Op(1)%nb_tot
        Ginv(i,i) = Ginv(i,i) + Ene
      END DO

      CALL inv_m1_TO_m2_cplx(Ginv,G,tab_Op(1)%nb_tot,0,ZERO)

      gGgG(:,:) = matmul(tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat,                 &
         matmul(G,matmul(tab_Op(para_CRP%iOp_CAP_Product)%Rmat,conjg(G))))

      crp2 = ZERO
      do i=1,tab_Op(1)%nb_tot
        crp2 = crp2 + real( gGgG(i,i), kind=Rkind)
      end do

      CALL dealloc_NParray(Ginv,'Ginv',name_sub)
      CALL dealloc_NParray(gGgG,'gGgG',name_sub)
      CALL dealloc_NParray(G,'G',name_sub)

      write(out_unitp,*) 'CRP at E (ua)', Ene, crp,'CRP with explicit inversion =', crp2
    ELSE
      RealTime = Delta_RealTime(CRP_Time)

      CALL alloc_NParray(Krylov_vectors,[tab_Op(1)%nb_tot,para_CRP%KS_max_it], &
                        'Krylov_vectors',name_sub,tab_lb=[1,0])

      CALL alloc_NParray(h,      [para_CRP%KS_max_it,para_CRP%KS_max_it],'h',name_sub)
      CALL alloc_NParray(Eigvals,[para_CRP%KS_max_it],                   'Eigvals',name_sub)

      ! Generate first Krylov vector randomly or from a guess (previous energy iteration)
      IF (size(GuessVec) /= tab_Op(1)%nb_tot) THEN
       write(out_unitp,*) ' ERROR in',name_sub
       write(out_unitp,*) '  The GuessVec size is wrong: ',size(GuessVec)
       write(out_unitp,*) '  H%nb_tot:                   ',tab_Op(1)%nb_tot
       write(out_unitp,*) ' CHECK the fortran source !'
       STOP ' ERROR in calc_crp_p_lanczos: The GuessVec size is wrong.'
      END IF
      IF (sqrt(dot_product(GuessVec,GuessVec)) == 0) THEN
        write(out_unitp,*) '  Random vector'
        CALL Random_CplxVec(GuessVec)
      END IF
      Krylov_vectors(:,0) = GuessVec


      IF (para_CRP%LinSolv_type == 'matinv') THEN

        CALL alloc_NParray(Ginv,shape(tab_Op(1)%Rmat),'Ginv',name_sub)
        CALL alloc_NParray(G,   shape(tab_Op(1)%Rmat),'G',   name_sub)
        CALL alloc_NParray(gGgG,shape(tab_Op(1)%Rmat),'gGgG',name_sub)

        Ginv(:,:) = -tab_Op(1)%Rmat + EYE*HALF * (tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat+ &
                                                  tab_Op(para_CRP%iOp_CAP_Product)%Rmat)

       DO i=1,tab_Op(1)%nb_tot
          Ginv(i,i) = Ginv(i,i) + Ene
       END DO

        CALL inv_m1_TO_m2_cplx(Ginv,G,tab_Op(1)%nb_tot,0,ZERO)

        gGgG(:,:) = matmul(tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat,               &
           matmul(G,matmul(tab_Op(para_CRP%iOp_CAP_Product)%Rmat,conjg(G))))

       CALL dealloc_NParray(Ginv,'Ginv',name_sub)
       CALL dealloc_NParray(G,   'G',   name_sub)
     ELSE IF (para_CRP%LinSolv_type == 'matlinsolv') THEN

       CALL alloc_NParray(Ginv,shape(tab_Op(1)%Rmat),'Ginv',name_sub)
       CALL alloc_NParray(trav,[tab_Op(1)%nb_tot],'trav',name_sub)
       CALL alloc_NParray(indx,[tab_Op(1)%nb_tot],'indx',name_sub)

       Ginv(:,:) = -tab_Op(1)%Rmat + EYE*HALF * (tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat+ &
                                                 tab_Op(para_CRP%iOp_CAP_Product)%Rmat)

       DO i=1,tab_Op(1)%nb_tot
          Ginv(i,i) = Ginv(i,i) + Ene
       END DO

       CALL ludcmp_cplx(Ginv,tab_Op(1)%nb_tot,trav,indx,d)
       !now in Ginv we have its LU decomposition

       CALL dealloc_NParray(trav,'trav',name_sub)
     ELSE IF (para_CRP%LinSolv_type == 'qmr' .OR. para_CRP%LinSolv_type == 'gmres') THEN
       CALL alloc_NParray(M1,[tab_Op(1)%nb_tot],'M1',name_sub)

       IF (allocated(tab_Op(1)%BasisnD%EneH0)) THEN
         M1(:) = ONE/(Ene-tab_Op(1)%BasisnD%EneH0(:)) ! approximation of 1/(Ene-H(i,i))
         write(out_unitp,*) 'precon /= 1. DML'
       ELSE
         M1(:)        = CONE
         write(out_unitp,*) 'precon = 1. DML'
       END IF
       !M1(:)        = CONE
       !write(out_unitp,*) 'precon = 1. DML'
     END IF

      ! Begin Lanczos scheme
      oldcrp = ZERO
      do nks=1,para_CRP%KS_max_it

         IF (print_level > 1) then
           write(out_unitp,*) '######################'
           write(out_unitp,*) '# in KS iterations, n=',nks
           write(out_unitp,*) '# before p_multiply'
           call flush_perso(out_unitp)
         end if

         SELECT CASE ( para_CRP%LinSolv_type )

         CASE ( 'matinv' )

            Krylov_vectors(:,nks) = matmul(gGgG,Krylov_vectors(:,nks-1))

         CASE ( 'matlinsolv')

            call p_multiplyLU(Krylov_vectors(:,nks-1),Krylov_vectors(:,nks),    &
                              tab_Op,nb_Op,Ene,tab_Op(1)%nb_tot,Ginv,indx,      &
                              para_CRP%iOp_CAP_Reactif,para_CRP%iOp_CAP_Product)

         CASE ( 'qmr' )

            call p_multiplyQMR(Krylov_vectors(:,nks-1),Krylov_vectors(:,nks),   &
                         tab_Op,nb_Op,Ene,ncooked,M1,para_CRP%LinSolv_accuracy, &
                         para_CRP%iOp_CAP_Reactif,para_CRP%iOp_CAP_Product)

         CASE ( 'gmres' )
#if __CERFACS == 1
            call p_multiplyGMRES(Krylov_vectors(:,nks-1),Krylov_vectors(:,nks), &
                         tab_Op,nb_Op,Ene,ncooked,M1,para_CRP%LinSolv_accuracy, &
                         para_CRP%iOp_CAP_Reactif,para_CRP%iOp_CAP_Product)

#else
           write(out_unitp,*) ' ERROR in',name_sub
           write(out_unitp,*) '  CERFACS GMRES is not implemented.'
           write(out_unitp,*) '  You have to choose between: "MatInv" or "QMR".'
           STOP ' ERROR CERFACS GMRES is not implemented'
#endif
         CASE Default
           write(out_unitp,*) ' ERROR in',name_sub
           write(out_unitp,*) '  No Default for LinSolv_type:',para_CRP%LinSolv_type
           write(out_unitp,*) '  You have to choose between: "MatInv" or "QMR".'
           STOP ' ERROR No Default for LinSolv_type'
         END SELECT
         call flush_perso(out_unitp)

         ! Calculate matrix
         IF (debug) write(out_unitp,*) '# in KS iterations, buiding h'
         do mks = 0, nks-1
            h(mks+1, nks) = dot_product(Krylov_vectors(:,mks),Krylov_vectors(:,nks))
            h(nks, mks+1) = conjg(h(mks+1,nks))
         end do
         IF (debug) write(out_unitp,*) '# in KS iterations, h'
         IF (debug) CALL Write_Mat(h(1:nks,1:nks),out_unitp,5)

         ! Orthogonalize vectors
         IF (debug) write(out_unitp,*) '# in KS iterations: Orthogonalize the vectors'
         do mks = 0, nks-1
            y = dot_product(Krylov_vectors(:,mks),Krylov_vectors(:,nks))
            Krylov_vectors(:,nks) = Krylov_vectors(:,nks) - y * Krylov_vectors(:,mks)
         end do
         ! Normalize vector
         CALL ReNorm_CplxVec(Krylov_vectors(:,nks))

         if (nks > 1) then

            Eigvals(:) = ZERO
            IF (allocated(EVec)) CALL dealloc_NParray(EVec,'EVec',name_sub)
            CALL alloc_NParray(EVec,[nks,nks],'EVec',name_sub)
            CALL diagonalization_HerCplx(h(1:nks,1:nks),Eigvals,EVec,nks,3,.FALSE.,.FALSE.)
            IF (debug) write(out_unitp,*) '# in KS iterations, Eigvals',Eigvals(1:nks)
            !write(out_unitp,*) '# in KS iterations, Eigvals',Eigvals(1:nks)

            crp = ZERO
            do i=1,nks
              crp = crp + h(i,i)
            end do

            DeltaCRP = oldcrp-crp
            if (abs(DeltaCRP) < para_CRP%KS_accuracy) EXIT
            oldcrp = crp
         endif

      !write(out_unitp,*) 'Krylov_vectors(:,0)',Krylov_vectors(:,0)
      !write(out_unitp,*) 'Krylov_vectors(:,1)',Krylov_vectors(:,1) ; stop

      end do

      !actual_iterations = nks
      write(out_unitp,*) '# in KS iterations, n=',nks
      write(out_unitp,*) 'accuracy: ',DeltaCRP
      IF (nks > para_CRP%KS_max_it .OR. DeltaCRP >= para_CRP%KS_accuracy) THEN
        write(out_unitp,*) 'CRP diago, minval: ',sum(Eigvals(1:para_CRP%KS_max_it)),&
                                         minval(abs(Eigvals(1:para_CRP%KS_max_it)))
        write(out_unitp,*) 'WARNING: Lanczos did not converged'
        nks = para_CRP%KS_max_it
      ELSE
        write(out_unitp,*) 'CRP diago, minval: ',sum(Eigvals(1:nks)),minval(abs(Eigvals(1:nks)))
      END IF

      RWU_E  = REAL_WU(Ene,'au','E')
      if (para_CRP%With_Eckart) then
        write(out_unitp,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                                     CRP,CRP_Eckart(Ene,para_CRP%Eckart),       &
                                     CRP-CRP_Eckart(Ene,para_CRP%Eckart)

      else
        write(out_unitp,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                          CRP
      end if
      RealTime = Delta_RealTime(CRP_Time)
      IF (debug .OR. print_Op .OR. print_level > 0) Then
        write(out_unitp,*) 'CRP Energy iteration: Delta Real Time',RealTime
      END IF
      CALL flush_perso(out_unitp)

      !CALL Random_CplxVec(GuessVec)
      GuessVec(:) = ZERO
      do mks = 0, nks-1
        GuessVec(:) = GuessVec(:) + Krylov_vectors(:,mks)* sum(EVec(mks+1,:))
      end do
      CALL ReNorm_CplxVec(GuessVec)
      IF (allocated(EVec)) CALL dealloc_NParray(EVec,'EVec',name_sub)


    end if

    IF (allocated(gGgG))           CALL dealloc_NParray(gGgG,'gGgG',name_sub)

    IF (allocated(M1))             CALL dealloc_NParray(M1,'M1',name_sub)
    IF (allocated(Krylov_vectors)) CALL dealloc_NParray(Krylov_vectors,'Krylov_vectors',name_sub)
    IF (allocated(h))              CALL dealloc_NParray(h,'h',name_sub)
    IF (allocated(Eigvals))        CALL dealloc_NParray(Eigvals,'Eigvals',name_sub)

    IF (allocated(indx))           CALL dealloc_NParray(indx,'indx',name_sub)
    IF (allocated(Ginv))           CALL dealloc_NParray(Ginv,'Ginv',name_sub)


END SUBROUTINE calc_crp_p_lanczos

SUBROUTINE calc_crp_IRL(tab_Op,nb_Op,para_CRP,Ene)

!---------------- Attempt at using ARPACK routine ---------------!
!                  Implicitly Restarted Lanczos                  !
!                  Scheme without shift-invert                   !
!----------------------------------------------------------------!

  USE mod_Constant
  USE mod_Op
  IMPLICIT NONE

!----- Operator variables ----------------------------------------------
  INTEGER,            INTENT(in)      :: nb_Op
  TYPE (param_Op),    INTENT(inout)   :: tab_Op(nb_Op)

  TYPE (param_CRP),   INTENT(in)      :: para_CRP
  REAL(kind=Rkind),   INTENT(in)      :: Ene

!----- working variables -----------------------------
  TYPE(REAL_WU)     :: RWU_E

      ! Calculate the inverse matrix explicitly? For debuging. It is equivalent to sub_CRP_BasisRep_WithMat
  LOGICAL, PARAMETER :: Inv = .FALSE.

      ! Size of Hamiltonian matrix
  INTEGER :: ncooked
  ! number of modes
  INTEGER :: ny
  ! number of possible excitations
  INTEGER :: nv

      ! Loop integers
  INTEGER :: i, j


  ! y harmonic frequency
  REAL (kind=Rkind) ::  wy

  ! zero point energy
  REAL (kind=Rkind) :: zpe


      !    IRL local arrays     !
  INTEGER           iparam(11), ipntr(14)
  LOGICAL, ALLOCATABLE ::  SELECT(:)
  COMPLEX (kind=Rkind), ALLOCATABLE :: &
       &                  ax(:), d(:), &
       &                  v(:,:), workd(:), &
       &                  workev(:), resid(:), &
       &                  workl(:)
  REAL (kind=Rkind), ALLOCATABLE ::   rwork(:), rd(:,:)

     !    IRL  local scalars    !
  CHARACTER         bmat*1, which*2
  INTEGER           ido, n, nx, nev, ncv, lworkl, info, &
       &                  ierr, nconv, maxitr, ishfts, mode
  INTEGER ldv
  COMPLEX (kind=Rkind)   sigma
  REAL (kind=Rkind)   tol
  LOGICAL           rvec

     ! BLAS & LAPACK routines used by IRL !
  REAL (kind=Rkind)  dznrm2 , dlapy2
  EXTERNAL          dznrm2 , zaxpy , dlapy2

      ! Operator variables ----------------------------------------------
  LOGICAL           :: print_Op

  INTEGER           ::    k,ie
  REAL (kind=Rkind) :: CRP, oldcrp,crp2,DeltaCRP
  COMPLEX (kind=Rkind), ALLOCATABLE :: G(:,:)
  COMPLEX (kind=Rkind), ALLOCATABLE :: Ginv(:,:)
  COMPLEX (kind=Rkind), ALLOCATABLE :: gGgG(:,:)
  COMPLEX (kind=Rkind), ALLOCATABLE :: M1(:)

       ! temporary variables for the LU decomposition
  INTEGER,             ALLOCATABLE :: indx(:)
  COMPLEX(kind=Rkind), ALLOCATABLE :: trav(:)
  COMPLEX(kind=Rkind)              :: dLU

  INTEGER :: lwork
  REAL (kind=Rkind), ALLOCATABLE :: work(:)

  TYPE (param_time) :: CRP_Time
  real(kind=Rkind)  :: RealTime

!----- for debuging --------------------------------------------------
  INTEGER   :: err
  LOGICAL, PARAMETER :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
  CHARACTER (len=*), PARAMETER :: name_sub = 'calc_crp_IRL'
!-----------------------------------------------------------

  RealTime = Delta_RealTime(CRP_Time)

!======================= END OF HEADER ==========================!

  ncooked   = tab_Op(1)%nb_tot

  IF (para_CRP%LinSolv_type == 'matinv') THEN

     CALL alloc_NParray(Ginv,SHAPE(tab_Op(1)%Rmat),'Ginv',name_sub)
     CALL alloc_NParray(G,   SHAPE(tab_Op(1)%Rmat),'G',   name_sub)
     CALL alloc_NParray(gGgG,SHAPE(tab_Op(1)%Rmat),'gGgG',name_sub)

     Ginv(:,:) = -tab_Op(1)%Rmat + EYE*HALF * (tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat+ &
          tab_Op(para_CRP%iOp_CAP_Product)%Rmat)

     DO i=1,tab_Op(1)%nb_tot
        Ginv(i,i) = Ginv(i,i) + Ene
     END DO

     CALL inv_m1_TO_m2_cplx(Ginv,G,tab_Op(1)%nb_tot,0,ZERO)

     gGgG(:,:) = MATMUL(tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat,               &
          MATMUL(G,MATMUL(tab_Op(para_CRP%iOp_CAP_Product)%Rmat,CONJG(G))))

     CALL dealloc_NParray(Ginv,'Ginv',name_sub)
     CALL dealloc_NParray(G,   'G',   name_sub)

      open(unit=17, file='P.dat')
      DO i=1,ncooked
         DO j=1, ncooked
            WRITE(17,*) gGgG(i,j)
         END DO
      END DO
      CLOSE (17)


     CRP = ZERO
     DO i=1,tab_Op(1)%nb_tot
        CRP = CRP + gGgG(i,i)
     END DO

     RWU_E  = REAL_WU(Ene,'au','E')

     WRITE(out_unitp,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
          CRP

     CALL dealloc_NParray(gGgG,   'gGgG',   name_sub)


  ELSE IF (para_CRP%LinSolv_type == 'matlinsolv') THEN

     CALL alloc_NParray(Ginv,SHAPE(tab_Op(1)%Rmat),'Ginv',name_sub)
     CALL alloc_NParray(trav,[tab_Op(1)%nb_tot],'trav',name_sub)
     CALL alloc_NParray(indx,[tab_Op(1)%nb_tot],'indx',name_sub)

     Ginv(:,:) = -tab_Op(1)%Rmat + EYE*HALF * (tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat+ &
          tab_Op(para_CRP%iOp_CAP_Product)%Rmat)

     DO i=1,tab_Op(1)%nb_tot
        Ginv(i,i) = Ginv(i,i) + Ene
     END DO

     CALL ludcmp_cplx(Ginv,tab_Op(1)%nb_tot,trav,indx,dLU)
         !now in Ginv we have its LU decomposition


     CALL dealloc_NParray(trav,'trav',name_sub)

  ELSE IF (para_CRP%LinSolv_type == 'qmr' .OR. para_CRP%LinSolv_type == 'gmres') THEN
     CALL alloc_NParray(M1,[tab_Op(1)%nb_tot],'M1',name_sub)

     IF (ALLOCATED(tab_Op(1)%BasisnD%EneH0)) THEN
        M1(:) = ONE/(Ene-tab_Op(1)%BasisnD%EneH0(:)) ! approximation of 1/(Ene-H(i,i))
        WRITE(out_unitp,*) 'precon /= 1. DML'
     ELSE
        M1(:)        = CONE
        WRITE(out_unitp,*) 'precon = 1. DML'
     END IF
  END IF
  CALL flush_perso(out_unitp)

!     %--------------------------------------------------%
!     | The number N(=NX*NX) is the dimension of the     |
!     | matrix.  A standard eigenvalue problem is        |
!     | solved (BMAT = 'I').  NEV is the number of       |
!     | eigenvalues to be approximated.  The user can    |
!     | modify NX, NEV, NCV, WHICH to solve problems of  |
!     | different sizes, and to get different parts of   |
!     | the spectrum.  However, The following            |
!     | conditions must be satisfied:                    |
!     |                   N <= MAXN                      |
!     |                 NEV <= MAXNEV                    |
!     |           NEV + 2 <= NCV <= MAXNCV               |
!     %--------------------------------------------------%

  n = ncooked
  ldv = n

      ! As we know the number of opened channels at a given energy,
      ! we know how many eigenvalues we want ???

  nev = ChannelNumber_AT_TS(Ene,para_CRP,tab_Op(1))

  ncv   = 2*nev ! recommended in manual

         ! array allocation for IRL
  ALLOCATE ( ax(n), d(ncv), &
       &     v(ldv,ncv), workd(3*n), &
       &     workev(3*ncv), resid(n), &
       &     workl(3*ncv*ncv+5*ncv), &
       &     SELECT(ncv), &
       &     rwork(ncv), rd(ncv,3) &
       & )

      ! Problem type: AX = l x
  bmat  = 'I'
      ! We want highest eigenvalues
  which = 'LM'
!
!     %---------------------------------------------------%
!     | The work array WORKL is used in ZNAUPD  as         |
!     | workspace.  Its dimension LWORKL is set as        |
!     | illustrated below.  The parameter TOL determines  |
!     | the stopping criterion. If TOL<=0, machine        |
!     | precision is used.  The variable IDO is used for  |
!     | reverse communication, and is initially set to 0. |
!     | Setting INFO=0 indicates that a random vector is  |
!     | generated to start the ARNOLDI iteration.         |
!     %---------------------------------------------------%
!
  lworkl  = 3*ncv**2+5*ncv
  tol    = para_CRP%KS_accuracy
  ido    = 0
  info   = 0
!
!     %---------------------------------------------------%
!     | This program uses exact shift with respect to     |
!     | the current Hessenberg matrix (IPARAM(1) = 1).    |
!     | IPARAM(3) specifies the maximum number of Arnoldi |
!     | iterations allowed.  Mode 1 of ZNAUPD  is used     |
!     | (IPARAM(7) = 1). All these options can be changed |
!     | by the user. For details see the documentation in |
!     | ZNAUPD .                                           |
!     %---------------------------------------------------%
!
  ishfts = 1
!  max iterations:
  maxitr = para_CRP%KS_max_it
  mode   = 1
!
  iparam(1) = ishfts
  iparam(3) = maxitr
  iparam(7) = mode

!----------------------------------------------------------------!
!-----------------------Begin Lanczos scheme---------------------!
!----------------------------------------------------------------!
!----------------------------------------------------------------!
  DO        ! Stop is handled by reverse communication  !
!----------------------------------------------------------------!
!
!     %-------------------------------------------%
!     | M A I N   L O O P (Reverse communication) |
!     %-------------------------------------------%
!
!
!        %---------------------------------------------%
!        | Repeatedly call the routine ZNAUPD  and take |
!        | actions indicated by parameter IDO until    |
!        | either convergence is indicated or maxitr   |
!        | has been exceeded.                          |
!        %---------------------------------------------%
!
#if __ARPACK == 1
     CALL znaupd  ( ido, bmat, n, which, nev, tol, resid, ncv,&
          &        v, ldv, iparam, ipntr, workd, workl, lworkl,&
          &        rwork,info )
#else
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) ' The ARPACK library is not present!'
        write(out_unitp,*) "Use CRP_Type='lanczos' instead of CRP_Type='lanczos_Arpack'"
        write(out_unitp,*) '  or recompile ElVibRot with ARPACK = 1 (makefile)'
        STOP 'ARPACK has been removed'
#endif
!
     IF (ido .EQ. -1 .OR. ido .EQ. 1) THEN
!
!           %-------------------------------------------%
!           | Perform matrix vector multiplication      |
!           |                y <--- OP*x                |
!           | takes workd(ipntr(1)) as the input vector,|
!           |  and return the matrix vector product     |
!           |         to workd(ipntr(2)).               |
!           %-------------------------------------------%
!
        SELECT CASE ( para_CRP%LinSolv_type )

        CASE ( 'matinv' )

           workd(ipntr(2):ipntr(2)+n-1) = MATMUL(gGgG,workd(ipntr(1):ipntr(1)+n-1))

        CASE ( 'matlinsolv')

           CALL p_multiplyLU(workd(ipntr(1):ipntr(1)+n-1),workd(ipntr(2):ipntr(2)+n-1),    &
                &  tab_Op,nb_Op,Ene,tab_Op(1)%nb_tot,Ginv,indx,      &
                &  para_CRP%iOp_CAP_Reactif,para_CRP%iOp_CAP_Product)

        CASE ( 'qmr' )

           CALL p_multiplyQMR(workd(ipntr(1):ipntr(1)+n-1),workd(ipntr(2):ipntr(2)+n-1),   &
                &  tab_Op,nb_Op,Ene,ncooked,M1,para_CRP%LinSolv_accuracy, &
                &  para_CRP%iOp_CAP_Reactif,para_CRP%iOp_CAP_Product)

        CASE ( 'gmres' )
#if __CERFACS == 1
           CALL p_multiplyGMRES(workd(ipntr(1):ipntr(1)+n-1),workd(ipntr(2):ipntr(2)+n-1), &
                &  tab_Op,nb_Op,Ene,ncooked,M1,para_CRP%LinSolv_accuracy, &
                &  para_CRP%iOp_CAP_Reactif,para_CRP%iOp_CAP_Product)

#else
           WRITE(out_unitp,*) ' ERROR in',name_sub
           WRITE(out_unitp,*) '  CERFACS GMRES is not implemented.'
           WRITE(out_unitp,*) '  You have to choose between: "MatInv" or "QMR".'
           STOP ' ERROR CERFACS GMRES is not implemented'
#endif
        CASE Default
           WRITE(out_unitp,*) ' ERROR in',name_sub
           WRITE(out_unitp,*) '  No Default for LinSolv_type:',para_CRP%LinSolv_type
           WRITE(out_unitp,*) '  You have to choose between: "MatInv" or "QMR".'
           STOP ' ERROR No Default for LinSolv_type'
        END SELECT
!
!           %-----------------------------------------%
!           | L O O P   B A C K to call ZNAUPD  again. |
!           %-----------------------------------------%
!
     ELSE
!
!     %----------------------------------------%
!     | Either we have convergence or there is |
!     | an error.                              |
!     %----------------------------------------%
!
        IF ( info .LT. 0 ) THEN
!
!        %--------------------------%
!        | Error message, check the |
!        | documentation in ZNAUPD   |
!        %--------------------------%
!
           PRINT *, ' '
           PRINT *, ' Error with _naupd, info = ', info
           PRINT *, ' Check the documentation of _naupd'
           PRINT *, ' '
!
        ELSE
!
!        %-------------------------------------------%
!        | No fatal errors occurred.                 |
!        | Post-Process using ZNEUPD .                |
!        |                                           |
!        | Computed eigenvalues may be extracted.    |
!        |                                           |
!        | Eigenvectors may also be computed now if  |
!        | desired.  (indicated by rvec = .true.)    |
!        %-------------------------------------------%
!
           rvec = .FALSE.
!
#if __ARPACK == 1
           CALL zneupd  (rvec, 'A', SELECT, d, v, ldv, sigma, &
                &        workev, bmat, n, which, nev, tol, resid, ncv, &
                &        v, ldv, iparam, ipntr, workd, workl, lworkl, &
                &        rwork, ierr)
#else
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) ' The ARPACK library is not present!'
        write(out_unitp,*) "Use CRP_Type='lanczos' instead of CRP_Type='lanczos_Arpack'"
        write(out_unitp,*) '  or recompile ElVibRot with ARPACK = 1 (makefile)'
        STOP 'ARPACK has been removed'
#endif

!
!        %----------------------------------------------%
!        | Eigenvalues are returned in the one          |
!        | dimensional array D.  The corresponding      |
!        | eigenvectors are returned in the first NCONV |
!        | (=IPARAM(5)) columns of the two dimensional  |
!        | array V if requested.  Otherwise, an         |
!        | orthogonal basis for the invariant subspace  |
!        | corresponding to the eigenvalues in D is     |
!        | returned in V.                               |
!        %----------------------------------------------%
!
           IF ( ierr .NE. 0) THEN
!
!           %------------------------------------%
!           | Error condition:                   |
!           | Check the documentation of ZNEUPD . |
!           %------------------------------------%
!
              PRINT *, ' '
              PRINT *, ' Error with _neupd, info = ', ierr
              PRINT *, ' Check the documentation of _neupd. '
              PRINT *, ' '
!
           ELSE
!
              nconv = iparam(5)
              DO j=1, nconv
!
!               %---------------------------%
!               | Compute the residual norm |
!               |                           |
!               |   ||  A*x - lambda*x ||   |
!               |                           |
!               | for the NCONV accurately  |
!               | computed eigenvalues and  |
!               | eigenvectors.  (iparam(5) |
!               | indicates how many are    |
!               | accurate to the requested |
!               | tolerance)                |
!               %---------------------------%
!
                 SELECT CASE ( para_CRP%LinSolv_type )

                 CASE ( 'matinv' )

                    ax = MATMUL(gGgG,v(:,j))

                 CASE ( 'matlinsolv')

                    CALL p_multiplyLU(v(1,j),ax,    &
                         &  tab_Op,nb_Op,Ene,tab_Op(1)%nb_tot,Ginv,indx,      &
                         &  para_CRP%iOp_CAP_Reactif,para_CRP%iOp_CAP_Product)

                 CASE ( 'qmr' )

                    CALL p_multiplyQMR(v(1,j),ax,   &
                         &  tab_Op,nb_Op,Ene,ncooked,M1,para_CRP%LinSolv_accuracy, &
                         &  para_CRP%iOp_CAP_Reactif,para_CRP%iOp_CAP_Product)

                 CASE ( 'gmres' )
#if __CERFACS == 1
                    CALL p_multiplyGMRES(v(1,j),ax, &
                         &  tab_Op,nb_Op,Ene,ncooked,M1,para_CRP%LinSolv_accuracy, &
                         &  para_CRP%iOp_CAP_Reactif,para_CRP%iOp_CAP_Product)

#else
                    WRITE(out_unitp,*) ' ERROR in',name_sub
                    WRITE(out_unitp,*) '  CERFACS GMRES is not implemented.'
                    WRITE(out_unitp,*) '  You have to choose between: "MatInv" or "QMR".'
                    STOP ' ERROR CERFACS GMRES is not implemented'
#endif
                 CASE Default
                    WRITE(out_unitp,*) ' ERROR in',name_sub
                    WRITE(out_unitp,*) '  No Default for LinSolv_type:',para_CRP%LinSolv_type
                    WRITE(out_unitp,*) '  You have to choose between: "MatInv" or "QMR".'
                    STOP ' ERROR No Default for LinSolv_type'
                 END SELECT

                 CALL zaxpy (n, -d(j), v(1,j), 1, ax, 1)
                 rd(j,1) = real(d(j))
                 rd(j,2) = aimag (d(j))
                 rd(j,3) = dznrm2 (n, ax, 1)
                 rd(j,3) = rd(j,3) / dlapy2 (rd(j,1),rd(j,2))
              END DO
!
!            %-----------------------------%
!            | Display computed residuals. |
!            %-----------------------------%
!
#if __ARPACK == 1
              CALL dmout (6, nconv, 3, rd, ncv, -6, &
                   &            'Ritz values (Real, Imag) and relative residuals')
#else
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) ' The ARPACK library is not present!'
        write(out_unitp,*) "Use CRP_Type='lanczos' instead of CRP_Type='lanczos_Arpack'"
        write(out_unitp,*) '  or recompile ElVibRot with ARPACK = 1 (makefile)'
        STOP 'ARPACK has been removed'
#endif
           END IF
!
!        %-------------------------------------------%
!        | Print additional convergence information. |
!        %-------------------------------------------%
!
           IF ( info .EQ. 1) THEN
              PRINT *, ' '
              PRINT *, ' Maximum number of iterations reached.'
              PRINT *, ' Eigvals '
              DO i=1, nev
                 PRINT *, REAL(D(i))
              END DO
              PRINT *, ' '
              nconv = 0 ! DML 18/01/2021
           ELSE IF ( info .EQ. 3) THEN
              PRINT *, ' '
              PRINT *, ' No shifts could be applied during implicit', &
                   & ' Arnoldi update, try increasing NCV.'
              PRINT *, ' '
           END IF

           PRINT *, ' '
           PRINT *, '_NDRV1'
           PRINT *, '====== '
           PRINT *, ' '
           PRINT *, ' Size of the matrix is ', n
           PRINT *, ' The number of Ritz values requested is ', nev
           PRINT *, ' The number of Arnoldi vectors generated',&
                &   ' (NCV) is ', ncv
           PRINT *, ' What portion of the spectrum: ', which
           PRINT *, ' The number of converged Ritz values is ',&
                &   nconv
           PRINT *, ' The number of Implicit Arnoldi update',&
                &   ' iterations taken is ', iparam(3)
           PRINT *, ' The number of OP*x is ', iparam(9)
           PRINT *, ' The convergence criterion is ', tol
           PRINT *, ' '

        END IF
!
!     %---------------------------%
!     |           Done.           |
!     %---------------------------%
!
        EXIT
     END IF
!----------------------------------------------------------------!
  END DO         !        End of main loop         !
!----------------------------------------------------------------!

  CRP = ZERO
  DO i=1,nconv
     CRP = CRP + real(D(i))
  END DO

  RWU_E  = REAL_WU(Ene,'au','E')

  IF (para_CRP%With_Eckart) THEN
     WRITE(out_unitp,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                                     CRP,CRP_Eckart(Ene,para_CRP%Eckart),       &
                                     CRP-CRP_Eckart(Ene,para_CRP%Eckart)

  ELSE
     WRITE(out_unitp,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
          CRP
  END IF
  RealTime = Delta_RealTime(CRP_Time)
  IF (debug .OR. print_Op .OR. print_level > 0) Then
    write(out_unitp,*) 'CRP Energy iteration: Delta Real Time',RealTime
  END IF
  CALL flush_perso(out_unitp)

  DEALLOCATE ( ax, d, &
       &     v, workd, &
       &     workev, resid, &
       &     workl, &
       &     SELECT, &
       &     rwork, rd &
       & )

  IF (ALLOCATED(gGgG))           CALL dealloc_NParray(gGgG,'gGgG',name_sub)
  IF (ALLOCATED(M1))             CALL dealloc_NParray(M1,'M1',name_sub)
  IF (ALLOCATED(indx))           CALL dealloc_NParray(indx,'indx',name_sub)
  IF (ALLOCATED(Ginv))           CALL dealloc_NParray(Ginv,'Ginv',name_sub)


END SUBROUTINE calc_crp_IRL
SUBROUTINE calc_crp_p_lanczos_old(tab_Op,nb_Op,para_CRP,Ene)

      USE mod_Constant
      USE mod_Op
      implicit none

!----- Operator variables ----------------------------------------------
      integer,            intent(in)      :: nb_Op
      TYPE (param_Op),    intent(inout)   :: tab_Op(nb_Op)

      TYPE (param_CRP),   intent(in)      :: para_CRP
      real(kind=Rkind),   intent(in)      :: Ene

!----- working variables -----------------------------
      TYPE(REAL_WU)     :: RWU_E

      ! Calculate the inverse matrix explicitly? For debuging. It is equivalent to sub_CRP_BasisRep_WithMat
      logical, parameter :: Inv = .FALSE.

      ! Size of Hamiltonian matrix
      integer :: ncooked
      ! Loop integers
      integer :: nks, mks, i, j
      ! Vectors for Pmult
      complex(kind=Rkind), dimension(:,:), allocatable:: Krylov_vectors,h
      real (kind=Rkind),    allocatable :: Eigvals(:)
      complex(kind=Rkind),  allocatable :: EVec(:,:)

      complex(kind=Rkind) y, len
      real(kind=Rkind) ranr, rani

      ! Operator variables ----------------------------------------------
      logical           :: print_Op


      integer           ::    k,ie
      real (kind=Rkind) :: CRP, oldcrp,crp2,DeltaCRP
      complex (kind=Rkind), allocatable :: G(:,:)
      complex (kind=Rkind), allocatable :: Ginv(:,:)
      complex (kind=Rkind), allocatable :: gGgG(:,:)
      complex (kind=Rkind), allocatable :: M1(:)

       ! temporary variables for the LU decomposition
       integer,             allocatable :: indx(:)
       complex(kind=Rkind), allocatable :: trav(:)
       complex(kind=Rkind)              :: d

       integer :: lwork,ierr
       real (kind=Rkind), ALLOCATABLE :: work(:)

!----- for debuging --------------------------------------------------
      integer   :: err
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'calc_crp_p_lanczos_old'
!-----------------------------------------------------------


      ncooked   = tab_Op(1)%nb_tot


    ! If need be, generate explicit matrix representation of operators
    IF ( Inv ) THEN

      CALL sub_MatOp(tab_Op(1),print_Op=.FALSE.) ! H
      DO i=3,nb_Op ! for the CAP
        IF (i == para_CRP%iOp_CAP_Reactif .OR. i == para_CRP%iOp_CAP_Product)   &
           CALL sub_MatOp(tab_Op(i),print_Op=.FALSE.)
      END DO

      CALL alloc_NParray(Ginv,shape(tab_Op(1)%Rmat),'Ginv',name_sub)
      CALL alloc_NParray(G,   shape(tab_Op(1)%Rmat),'G',   name_sub)
      CALL alloc_NParray(gGgG,shape(tab_Op(1)%Rmat),'gGgG',name_sub)

      Ginv(:,:) = -tab_Op(1)%Rmat + EYE*HALF * (tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat+ &
                                                tab_Op(para_CRP%iOp_CAP_Product)%Rmat)

      DO i=1,tab_Op(1)%nb_tot
        Ginv(i,i) = Ginv(i,i) + Ene
      END DO

      CALL inv_m1_TO_m2_cplx(Ginv,G,tab_Op(1)%nb_tot,0,ZERO)

      gGgG(:,:) = matmul(tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat,                 &
         matmul(G,matmul(tab_Op(para_CRP%iOp_CAP_Product)%Rmat,conjg(G))))

      crp2 = ZERO
      do i=1,tab_Op(1)%nb_tot
        crp2 = crp2 + real( gGgG(i,i), kind=Rkind)
      end do

      CALL dealloc_NParray(Ginv,'Ginv',name_sub)
      CALL dealloc_NParray(gGgG,'gGgG',name_sub)
      CALL dealloc_NParray(G,'G',name_sub)

      write(out_unitp,*) 'CRP at E (ua)', Ene, crp,'CRP with explicit inversion =', crp2
    ELSE

      CALL alloc_NParray(Krylov_vectors,[tab_Op(1)%nb_tot,para_CRP%KS_max_it], &
                        'Krylov_vectors',name_sub,tab_lb=[1,0])

      CALL alloc_NParray(h,      [para_CRP%KS_max_it,para_CRP%KS_max_it],'h',name_sub)
      CALL alloc_NParray(Eigvals,[para_CRP%KS_max_it],                   'Eigvals',name_sub)

      ! Generate first Krylov vector randomly
      do i =1,tab_Op(1)%nb_tot
         CALL random_number(ranr)
         CALL random_number(rani)
         Krylov_vectors(i,0) = cmplx(ranr,rani,kind=Rkind)
      end do
      CALL ReNorm_CplxVec(Krylov_vectors(:,0))

      IF (para_CRP%LinSolv_type == 'matinv') THEN
        CALL sub_MatOp(tab_Op(1),print_Op=.FALSE.) ! H
        DO i=3,nb_Op ! for the CAP
          IF (i == para_CRP%iOp_CAP_Reactif .OR. i == para_CRP%iOp_CAP_Product) &
             CALL sub_MatOp(tab_Op(i),print_Op=.FALSE.)
        END DO

        CALL alloc_NParray(Ginv,shape(tab_Op(1)%Rmat),'Ginv',name_sub)
        CALL alloc_NParray(G,   shape(tab_Op(1)%Rmat),'G',   name_sub)
        CALL alloc_NParray(gGgG,shape(tab_Op(1)%Rmat),'gGgG',name_sub)

        Ginv(:,:) = -tab_Op(1)%Rmat + EYE*HALF * (tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat+ &
                                                  tab_Op(para_CRP%iOp_CAP_Product)%Rmat)

       DO i=1,tab_Op(1)%nb_tot
          Ginv(i,i) = Ginv(i,i) + Ene
       END DO

        CALL inv_m1_TO_m2_cplx(Ginv,G,tab_Op(1)%nb_tot,0,ZERO)

        gGgG(:,:) = matmul(tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat,               &
           matmul(G,matmul(tab_Op(para_CRP%iOp_CAP_Product)%Rmat,conjg(G))))

       CALL dealloc_NParray(Ginv,'Ginv',name_sub)
       CALL dealloc_NParray(G,   'G',   name_sub)
     ELSE IF (para_CRP%LinSolv_type == 'matlinsolv') THEN
       CALL sub_MatOp(tab_Op(1),print_Op=.FALSE.) ! H
       DO i=3,nb_Op ! for the CAP
         IF (i == para_CRP%iOp_CAP_Reactif .OR. i == para_CRP%iOp_CAP_Product) &
            CALL sub_MatOp(tab_Op(i),print_Op=.FALSE.)
       END DO

       CALL alloc_NParray(Ginv,shape(tab_Op(1)%Rmat),'Ginv',name_sub)
       CALL alloc_NParray(trav,[tab_Op(1)%nb_tot],'trav',name_sub)
       CALL alloc_NParray(indx,[tab_Op(1)%nb_tot],'indx',name_sub)

       Ginv(:,:) = -tab_Op(1)%Rmat + EYE*HALF * (tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat+ &
                                                 tab_Op(para_CRP%iOp_CAP_Product)%Rmat)

       DO i=1,tab_Op(1)%nb_tot
          Ginv(i,i) = Ginv(i,i) + Ene
       END DO

       CALL ludcmp_cplx(Ginv,tab_Op(1)%nb_tot,trav,indx,d)
       !now in Ginv we have its LU decomposition

       CALL dealloc_NParray(trav,'trav',name_sub)
     ELSE IF (para_CRP%LinSolv_type == 'qmr' .OR. para_CRP%LinSolv_type == 'gmres') THEN
       CALL alloc_NParray(M1,[tab_Op(1)%nb_tot],'M1',name_sub)

       IF (allocated(tab_Op(1)%BasisnD%EneH0)) THEN
         M1(:) = ONE/(Ene-tab_Op(1)%BasisnD%EneH0(:)) ! approximation of 1/(Ene-H(i,i))
         write(out_unitp,*) 'precon /= 1. DML'
       ELSE
         M1(:)        = CONE
         write(out_unitp,*) 'precon = 1. DML'
       END IF
       !M1(:)        = CONE
       !write(out_unitp,*) 'precon = 1. DML'
     END IF

      ! Begin Lanczos scheme
      oldcrp = ZERO
      do nks=1,para_CRP%KS_max_it

         IF (print_level > 1) then
           write(out_unitp,*) '######################'
           write(out_unitp,*) '# in KS iterations, n=',nks
           write(out_unitp,*) '# before p_multiply'
           call flush_perso(out_unitp)
         end if

         SELECT CASE ( para_CRP%LinSolv_type )

         CASE ( 'matinv' )

            Krylov_vectors(:,nks) = matmul(gGgG,Krylov_vectors(:,nks-1))

         CASE ( 'matlinsolv')

            call p_multiplyLU(Krylov_vectors(:,nks-1),Krylov_vectors(:,nks),    &
                              tab_Op,nb_Op,Ene,tab_Op(1)%nb_tot,Ginv,indx,      &
                              para_CRP%iOp_CAP_Reactif,para_CRP%iOp_CAP_Product)

         CASE ( 'qmr' )

            call p_multiplyQMR(Krylov_vectors(:,nks-1),Krylov_vectors(:,nks),   &
                         tab_Op,nb_Op,Ene,ncooked,M1,para_CRP%LinSolv_accuracy, &
                         para_CRP%iOp_CAP_Reactif,para_CRP%iOp_CAP_Product)

         CASE ( 'gmres' )
#if __CERFACS == 1
            call p_multiplyGMRES(Krylov_vectors(:,nks-1),Krylov_vectors(:,nks), &
                         tab_Op,nb_Op,Ene,ncooked,M1,para_CRP%LinSolv_accuracy, &
                         para_CRP%iOp_CAP_Reactif,para_CRP%iOp_CAP_Product)

#else
           write(out_unitp,*) ' ERROR in',name_sub
           write(out_unitp,*) '  CERFACS GMRES is not implemented.'
           write(out_unitp,*) '  You have to choose between: "MatInv" or "QMR".'
           STOP ' ERROR CERFACS GMRES is not implemented'
#endif
         CASE Default
           write(out_unitp,*) ' ERROR in',name_sub
           write(out_unitp,*) '  No Default for LinSolv_type:',para_CRP%LinSolv_type
           write(out_unitp,*) '  You have to choose between: "MatInv" or "QMR".'
           STOP ' ERROR No Default for LinSolv_type'
         END SELECT

         ! Calculate matrix
         do mks = 0, nks-1
            h(mks+1, nks) = dot_product(Krylov_vectors(:,mks),Krylov_vectors(:,nks))
            h(nks, mks+1) = conjg(h(mks+1,nks))
         end do
         IF (debug) write(out_unitp,*) '# in KS iterations, h'
         IF (debug) CALL Write_Mat(h(1:nks,1:nks),out_unitp,5)

         ! Orthogonalize vectors
         do mks = 0, nks-1
            y = dot_product(Krylov_vectors(:,mks),Krylov_vectors(:,nks))
            Krylov_vectors(:,nks) = Krylov_vectors(:,nks) - y * Krylov_vectors(:,mks)
         end do
         ! Normalize vector
         CALL ReNorm_CplxVec(Krylov_vectors(:,nks))

         if (nks > 1) then

            Eigvals(:) = ZERO
            CALL alloc_NParray(EVec,[nks,nks],'EVec',name_sub)
            CALL diagonalization_HerCplx(h(1:nks,1:nks),Eigvals,EVec,nks,3,.FALSE.,.FALSE.)
            IF (debug) write(out_unitp,*) '# in KS iterations, Eigvals',Eigvals(1:nks)
            write(out_unitp,*) '# in KS iterations, Eigvals',Eigvals(1:nks)
            CALL dealloc_NParray(EVec,'EVec',name_sub)

            crp = ZERO
            do i=1,nks
              crp = crp + h(i,i)
            end do

            DeltaCRP = oldcrp-crp
            if (abs(DeltaCRP) < para_CRP%KS_accuracy) EXIT
            oldcrp = crp
         endif

      !write(out_unitp,*) 'Krylov_vectors(:,0)',Krylov_vectors(:,0)
      !write(out_unitp,*) 'Krylov_vectors(:,1)',Krylov_vectors(:,1) ; stop

      end do

      !actual_iterations = nks
      write(out_unitp,*) '# in KS iterations, n=',nks
      write(out_unitp,*) 'CRP diago, minval: ',sum(Eigvals(1:nks)),minval(abs(Eigvals(1:nks)))
      write(out_unitp,*) 'accuracy: ',DeltaCRP
      IF (nks > para_CRP%KS_max_it .OR. DeltaCRP >= para_CRP%KS_accuracy) THEN
         write(out_unitp,*) 'WARNING: Lanczos did not converge'
      END IF

      RWU_E  = REAL_WU(Ene,'au','E')
      if (para_CRP%With_Eckart) then
        write(out_unitp,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                          CRP,CRP_Eckart(Ene,para_CRP%Eckart)
      else
        write(out_unitp,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                          CRP
      end if

    end if

    IF (allocated(gGgG))           CALL dealloc_NParray(gGgG,'gGgG',name_sub)

    IF (allocated(M1))             CALL dealloc_NParray(M1,'M1',name_sub)
    IF (allocated(Krylov_vectors)) CALL dealloc_NParray(Krylov_vectors,'Krylov_vectors',name_sub)
    IF (allocated(h))              CALL dealloc_NParray(h,'h',name_sub)
    IF (allocated(Eigvals))        CALL dealloc_NParray(Eigvals,'Eigvals',name_sub)

    IF (allocated(indx))           CALL dealloc_NParray(indx,'indx',name_sub)
    IF (allocated(Ginv))           CALL dealloc_NParray(Ginv,'Ginv',name_sub)


END SUBROUTINE calc_crp_p_lanczos_old
SUBROUTINE p_multiplyLU(Vin,Vut,tab_Op,nb_Op,Ene,N,Ginv_LU,indx,                &
                        iOp_CAP_Reactif,iOp_CAP_Product)
      use mod_system
      USE mod_Op
      implicit none

      integer,             intent(in)    :: N
      complex(kind=Rkind), intent(in)    :: Vin(N)
      complex(kind=Rkind), intent(inout) :: Vut(N)
!----- Operator variables ----------------------------------------------
      integer,             intent(in)    :: nb_Op,iOp_CAP_Reactif,iOp_CAP_Product
      TYPE (param_Op),     intent(in)    :: tab_Op(nb_Op)
      real (kind=Rkind),   intent(in)    :: Ene
      complex(kind=Rkind), intent(in)    :: Ginv_LU(N,N)
      integer,             intent(in)    :: indx(N)



      complex(kind=Rkind) :: b(N)


!----- for debuging --------------------------------------------------
      integer   :: err
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub ='p_multiplyLU'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Vin',Vin(:)
        write(out_unitp,*)
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

!     |b>=e_r|0>
      b(:)=Vin(:)
      call OpOnVec(b,tab_Op(iOp_CAP_Reactif),'NOC')
      IF (debug) write(out_unitp,*) 'e_r |Vin>',b(:)

!     |b>=1/(H-E-ie)|b>
      IF (print_level > 1) write(out_unitp,*) '# here before LU 1 '
      b(:) = conjg(b)
      CALL lubksb_cplx(Ginv_LU,N,indx,b)
      b(:) = conjg(b)
      IF (debug) write(out_unitp,*) '1/(H-E-ie)|b>',b(:)

!     |b>=e_p|x>
      call OpOnVec(b,tab_Op(iOp_CAP_Product),'NOC')
      IF (debug) write(out_unitp,*) 'e_p |b>',b(:)

!     |b>=1/(H-E+ie)|b>
      IF (print_level > 1) write(out_unitp,*) '# here before LU 2 '
      CALL lubksb_cplx(Ginv_LU,N,indx,b)

      Vut(:)=b(:)

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Vut',Vut(:)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

END SUBROUTINE p_multiplyLU

SUBROUTINE Gpsi(Vect,tab_Op,nb_Op,Ene,iOp_CAP_Reactif,iOp_CAP_Product,l_conjg)
!SUBROUTINE Gpsi(Vect,tab_Op,nb_Op,Ene,l_conjg)
      use mod_system
      USE mod_psi,     ONLY : param_psi,alloc_psi,dealloc_psi
      USE mod_Op
      implicit none

      integer,             intent(in)           :: nb_Op
      TYPE (param_Op)                           :: tab_Op(nb_Op)
      integer,             intent(in)           :: iOp_CAP_Reactif,iOp_CAP_Product
      complex(kind=Rkind), intent(inout)        :: Vect(tab_Op(1)%nb_tot)
      real(kind=Rkind),    intent(in)           :: Ene
      character(len=3),    intent(in)           :: l_conjg

      TYPE (param_psi)   :: Psi
      TYPE (param_psi)   :: OpPsi
      integer            :: i

      if (l_conjg == 'CJG') then
         Vect(:)=conjg(Vect(:))
      end if

      CALL init_psi(Psi,tab_Op(1),cplx=.TRUE.)
      CALL alloc_psi(Psi,BasisRep=.TRUE.,GridRep=.FALSE.)

      Psi%cvecB(:)=Vect(:)
      OpPsi = Psi

      call sub_OpPsi(Psi,OpPsi,tab_Op(1))

      Vect(:)= Vect(:)*Ene - OpPsi%cvecB(:)


      call sub_OpPsi(Psi,OpPsi,tab_Op(iOp_CAP_Reactif))

      Vect(:) = Vect(:)+EYE*HALF*OpPsi%cvecB(:)

      call sub_OpPsi(Psi,OpPsi,tab_Op(iOp_CAP_Product))

      Vect(:) = Vect(:)+EYE*HALF*OpPsi%cvecB(:)

      if (l_conjg == 'CJG') then
         Vect(:)=conjg(Vect(:))
      end if

      call dealloc_psi(Psi, .TRUE.)
      call dealloc_psi(OpPsi, .TRUE.)

END SUBROUTINE Gpsi

SUBROUTINE G_Mat(H,CAP_Reactif,CAP_Product,Ene,G)
      use mod_system
      USE mod_Op
      implicit none

      TYPE (param_Op),      intent(in)           :: H,CAP_Reactif,CAP_Product
      real(kind=Rkind),     intent(in)           :: Ene

      complex (kind=Rkind), intent(inout)        :: G(H%nb_tot,H%nb_tot)


      complex (kind=Rkind) :: Ginv(H%nb_tot,H%nb_tot)
      integer              :: i

!----- for debuging --------------------------------------------------
      integer   :: err
      logical, parameter :: debug=.FALSE.
!      logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'G_Mat'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Ene',Ene
        write(out_unitp,*)
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

      Ginv(:,:) = -H%Rmat + EYE*HALF * (CAP_Reactif%Rmat+CAP_Product%Rmat)

      DO i=1,H%nb_tot
        Ginv(i,i) = Ginv(i,i) + Ene
      END DO

      CALL inv_m1_TO_m2_cplx(Ginv,G,H%nb_tot,0,ZERO)

      IF (debug) THEN
        Ginv = matmul(Ginv,G)
        DO i=1,H%nb_tot
          Ginv(i,i) = Ginv(i,i) - CONE
        END DO
        write(out_unitp,*) 'id diff ?',maxval(abs(Ginv))
    END IF

    IF (debug) THEN
      write(out_unitp,*) 'END ',name_sub
      CALL flush_perso(out_unitp)
    END IF

END SUBROUTINE G_Mat

SUBROUTINE FluxOp_Mat(H,HStep_Op,FluxOp)
      use mod_system
      USE mod_Op
      implicit none

      TYPE (param_Op)                           :: H,HStep_Op
      real(kind=Rkind),    intent(inout)        :: FluxOp(H%nb_tot,H%nb_tot)



      complex(kind=Rkind)        :: Rvp(H%nb_tot,H%nb_tot)
      real(kind=Rkind)           :: Rdiag(H%nb_tot)

      integer :: i,nb_col


      FluxOp = EYE*(matmul(H%Rmat,HStep_Op%Rmat) - matmul(HStep_Op%Rmat,H%Rmat))


      CALL diagonalization_HerCplx(FluxOp,Rdiag,Rvp,H%nb_tot,3,2,.TRUE.)

      write(out_unitp,*) 'Eigenvalues'
      DO i=1,H%nb_tot
        write(out_unitp,*) i,Rdiag(i)
      END DO

      nb_col = 5
      write(out_unitp,*) 'Flux eigenvectors in column'
      write(out_unitp,*) nb_col,H%nb_tot,H%nb_tot
      CALL Write_Mat(Rvp,out_unitp,nb_col)


      write(out_unitp,*) 'Ortho ?'
      Rvp = matmul(transpose(Rvp),Rvp)
      CALL Write_Mat(Rvp,out_unitp,nb_col)


      !write(out_unitp,*) 'Diag',Rdiag

END SUBROUTINE FluxOp_Mat
SUBROUTINE FluxOp_Mat_old(H,HStep_Op,FluxOp)
      use mod_system
      USE mod_Op
      implicit none

      TYPE (param_Op)                           :: H,HStep_Op
      real(kind=Rkind),    intent(inout)        :: FluxOp(H%nb_tot,H%nb_tot)



      real(kind=Rkind)        :: Rdiag(H%nb_tot),Rvp(H%nb_tot,H%nb_tot)
      integer :: i,nb_col


      FluxOp = matmul(H%Rmat,HStep_Op%Rmat) - matmul(HStep_Op%Rmat,H%Rmat)
      !CALL  diagonalization(FluxOp,Rdiag,Rvp,H%nb_tot,4,1,.TRUE.)

      FluxOp = matmul(FluxOp,FluxOp)
      CALL  diagonalization(FluxOp,Rdiag,Rvp,H%nb_tot,3,1,.TRUE.)

      write(out_unitp,*) 'Eigenvalues'
      DO i=1,H%nb_tot
        write(out_unitp,*) i,Rdiag(i)
      END DO

      nb_col = 5
      write(out_unitp,*) 'Flux eigenvectors in column'
      write(out_unitp,*) nb_col,H%nb_tot,H%nb_tot
      CALL Write_Mat(Rvp,out_unitp,nb_col)


      !write(out_unitp,*) 'Ortho ?'
      !Rvp = matmul(transpose(Rvp),Rvp)
      !CALL Write_Mat(Rvp,out_unitp,nb_col)


      !write(out_unitp,*) 'Diag',Rdiag

END SUBROUTINE FluxOp_Mat_old
SUBROUTINE FluxOp_Mat_v0(H,HStep_Op,FluxOp)
      use mod_system
      USE mod_Op
      implicit none

      TYPE (param_Op)                           :: H,HStep_Op
      real(kind=Rkind),    intent(inout)        :: FluxOp(H%nb_tot,H%nb_tot)



      real(kind=Rkind)        :: Rdiag(H%nb_tot),Rvp(H%nb_tot,H%nb_tot)
      integer :: i,nb_col


      FluxOp = matmul(H%Rmat,HStep_Op%Rmat) - matmul(HStep_Op%Rmat,H%Rmat)
      CALL  diagonalization(FluxOp,Rdiag,Rvp,H%nb_tot,4,1,.TRUE.)

      ! WARNNING: since FluxOp is a skew matrix, they are pairs of  eigenvalues (i*wk, -iwk) ...
      ! and the eigenvectors are V(:,k) + i*V(:,k+1) and V(:,k) - i*V(:,k+1)
      ! its means the V(:,k) and V(:,k+1) are not normalized to one
      DO i=1,H%nb_tot
        Rvp(:,i) = Rvp(:,i)/sqrt(dot_product(Rvp(:,i),Rvp(:,i)))
      end do
      nb_col = 5
      write(out_unitp,*) 'Flux eigenvectors in column'
      write(out_unitp,*) nb_col,H%nb_tot,H%nb_tot
      CALL Write_Mat(Rvp,out_unitp,nb_col)


      !write(out_unitp,*) 'Ortho ?'
      !Rvp = matmul(transpose(Rvp),Rvp)
      !CALL Write_Mat(Rvp,out_unitp,nb_col)


      !write(out_unitp,*) 'Diag',Rdiag

END SUBROUTINE FluxOp_Mat_v0
SUBROUTINE OpOnVec(Vect,tab_Op,l_conjg)
      use mod_system
      USE mod_psi,     ONLY : param_psi,alloc_psi,dealloc_psi
      USE mod_Op
      implicit none

      TYPE (param_Op)     :: tab_Op
      logical             :: print_Op
      logical, parameter  :: cplx=.TRUE.
      TYPE (param_psi)    :: Psi
      TYPE (param_psi)    :: OpPsi
      character(len=3)    :: l_conjg
      complex(kind=Rkind), dimension(tab_Op%nb_tot) :: Vect
      integer         :: i

      if (l_conjg == 'CJG') then
         Vect(:)=conjg(Vect(:))
      end if

      CALL init_psi(Psi,tab_Op,cplx)
      CALL alloc_psi(Psi,BasisRep=.TRUE.,GridRep=.FALSE.)

      Psi%cvecB(:)=Vect(:)
      OpPsi = Psi

      call sub_OpPsi(Psi,OpPsi,tab_Op)

      Vect(:) = OpPsi%cvecB(:)

      if (l_conjg == 'CJG') then
         Vect(:)=conjg(Vect(:))
      end if

      call dealloc_psi(Psi, .TRUE.)
      call dealloc_psi(OpPsi, .TRUE.)
END SUBROUTINE OpOnVec

SUBROUTINE ReNorm_CplxVec(Vect)
      use mod_system
      implicit none

      complex (kind=Rkind), intent(inout) :: Vect(:)

      Vect(:) = Vect(:)/sqrt(dot_product(Vect,Vect))

END SUBROUTINE ReNorm_CplxVec
SUBROUTINE Random_CplxVec(Vect)
      use mod_system
      implicit none

      complex (kind=Rkind), intent(inout) :: Vect(:)

      integer           :: i
      real (kind=Rkind) :: ranr,rani

      ! Generate first Krylov vector randomly
      do i =1,size(Vect)
         CALL random_number(ranr)
         CALL random_number(rani)
         Vect(i) = cmplx(ranr,rani,kind=Rkind)
      end do

      Vect(:) = Vect(:)/sqrt(dot_product(Vect,Vect))

END SUBROUTINE Random_CplxVec
SUBROUTINE SchmidtProjectOut_CplxVec(Vect,tab_Vect)
      use mod_system
      implicit none

      complex (kind=Rkind), intent(inout) :: Vect(:)
      complex (kind=Rkind), intent(in)    :: tab_Vect(:,:)

      integer               :: i
      complex (kind=Rkind)  :: s


    ! Orthogonalize vectors
    do i = lbound(tab_Vect,dim=2),ubound(tab_Vect,dim=2)
      s = dot_product(Vect,tab_Vect(:,i))
      Vect(:) = Vect(:) - s * tab_Vect(:,i)
    end do
    Vect(:) = Vect(:)/sqrt(dot_product(Vect,Vect))

END SUBROUTINE SchmidtProjectOut_CplxVec
FUNCTION CRP_Eckart(E,Eckart)
  USE mod_system
  IMPLICIT NONE
  real (kind=Rkind)                 :: CRP_Eckart
  real (kind=Rkind),    intent(in)  :: E
  TYPE (CRP_Eckart_t),  intent(in)  :: Eckart


      real (kind=Rkind) :: b,c
      !real (kind=Rkind), parameter :: V0=0.0156_Rkind,m=1060._Rkind,L=ONE
      !real (kind=Rkind), parameter :: V0=0.015625_Rkind,m=1061._Rkind,L=ONE

       b = Eckart%L * Pi * sqrt(TWO*Eckart%m*E)
       c = (Pi/TWO) * sqrt(EIGHT * Eckart%V0*Eckart%m*Eckart%L**2 - ONE)

       CRP_Eckart = ONE / (ONE + (cosh(c)/sinh(b))**2)

END FUNCTION CRP_Eckart
FUNCTION combination(nv,ny)

  IMPLICIT NONE
  INTEGER nv, ny
  INTEGER combination
  INTEGER i


  combination=1

  DO i=nv+1,nv+ny-1
     combination = combination*i
  END DO

  DO i=2,ny-1
     combination=combination/i
  END DO

END FUNCTION combination

SUBROUTINE Read_Channel_AT_TS(Channel_AT_TS_var,ny)
USE mod_system
USE mod_RealWithUnit
USE mod_dnSVM
USE mod_nDindex
USE mod_Op
IMPLICIT NONE

  TYPE (CRP_Channel_AT_TS_t),     intent(inout)      :: Channel_AT_TS_var
  integer,                        intent(in)         :: ny

  integer                         :: option,err_unit
  TYPE (REAL_WU)                  :: w1,EneTS
  character (len=Name_len)        :: w_unit = 'cm-1'
  real (kind=Rkind)               :: conv
  integer                         :: nb_channels_added

  NAMELIST / Channel_AT_TS / option,w1,EneTS,w_unit,nb_channels_added


!----- for debuging --------------------------------------------------
      integer   :: err
      !logical, parameter :: debug=.FALSE.
      logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'Read_Channel_AT_TS'
!-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*) 'ny',ny
    CALL flush_perso(out_unitp)
  END IF
!-----------------------------------------------------------


  w_unit              = 'cm-1'
  EneTS               = REAL_WU(0.0105_Rkind,'au','E')
  w1                  = REAL_WU(0.015625_Rkind,'au','E')
  option              = 1
  nb_channels_added   = 1

  read(in_unitp,Channel_AT_TS)
  write(out_unitp,Channel_AT_TS)

  Channel_AT_TS_var%EneTS             = convRWU_TO_RWU(EneTS)
  Channel_AT_TS_var%w1                = convRWU_TO_RWU(w1)
  Channel_AT_TS_var%option            = option
  Channel_AT_TS_var%nb_channels_added = nb_channels_added

  IF (option == 2) Then
    conv = get_Conv_au_TO_unit(quantity='E',Unit=w_unit,err_unit=err_unit)
    IF (err_unit /= 0) STOP 'in Read_Channel_AT_TS: Wrong w_unit !'
    write(out_unitp,*) 'For w_unit= "',trim(w_unit),'", conv=',conv

    allocate(Channel_AT_TS_var%w(ny))
    read(in_unitp,*) Channel_AT_TS_var%w(:)
    Channel_AT_TS_var%w(:) = Channel_AT_TS_var%w(:)/conv
  END IF


  IF (debug) THEN
    CALL Write_Channel_AT_TS(Channel_AT_TS_var)
    write(out_unitp,*)
    write(out_unitp,*) 'END ',name_sub
    CALL flush_perso(out_unitp)
  END IF


END SUBROUTINE Read_Channel_AT_TS
SUBROUTINE Write_Channel_AT_TS(Channel_AT_TS)
USE mod_system
USE mod_Constant
USE mod_dnSVM
USE mod_nDindex
USE mod_Op
IMPLICIT NONE

  TYPE (CRP_Channel_AT_TS_t),     intent(in)      :: Channel_AT_TS

!----- for debuging --------------------------------------------------
      integer   :: err
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'Write_Channel_AT_TS'
!-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    CALL flush_perso(out_unitp)
  END IF
!-----------------------------------------------------------

    write(out_unitp,*) 'option            ',Channel_AT_TS%option
    write(out_unitp,*) 'EneTS (au)        ',Channel_AT_TS%EneTS
    write(out_unitp,*) 'nb_channels_added ',Channel_AT_TS%nb_channels_added
    SELECT CASE (Channel_AT_TS%option)
    CASE(1)
      write(out_unitp,*) 'w1 (au) ',Channel_AT_TS%w1
    CASE(2)
      IF (allocated(Channel_AT_TS%w)) THEN
        write(out_unitp,*) 'w(:) (au) ',Channel_AT_TS%w
      ELSE
        write(out_unitp,*) 'w(:) is not allocated !'
      END IF
    END SELECT

  IF (debug) THEN
    write(out_unitp,*)
    write(out_unitp,*) 'END ',name_sub
    CALL flush_perso(out_unitp)
  END IF


END SUBROUTINE Write_Channel_AT_TS
FUNCTION ChannelNumber_AT_TS(Ene,para_CRP,para_H) RESULT(nb_channels)
USE mod_system
USE mod_dnSVM
USE mod_nDindex
USE mod_Op
IMPLICIT NONE

  INTEGER                               :: nb_channels
  TYPE (param_CRP),     intent(in)      :: para_CRP
  real (kind=Rkind),    intent(in)      :: Ene
  TYPE (param_Op),      intent(in)      :: para_H


  integer :: option = 2 ! 1: Lucien Dupuy, one degenerate frequency
                        ! 2: ny frequencies at TS
                        ! 3: general, Energy levels calculation at the TS

  ! parameters for option=1 (Lucien)
  integer             :: i,ny,nv
  real (kind=Rkind)   :: wy,zpe,EneTS

  ! more general parameters for option=2
  TYPE (Type_nDindex)             :: nDindB_Channels
  TYPE (Type_IntVec), allocatable :: tab_i_TO_l(:)
  integer,            allocatable :: nbSize(:),tab_ib(:)
  integer                         :: LB,ib,nb,n
  real (kind=Rkind),  allocatable :: w(:)
  real (kind=Rkind)               :: EneChannel

  ! more general parameters for option=3
  TYPE (basis),        pointer    :: basisnD
  integer,            allocatable :: nDval(:)
  real (kind=Rkind)               :: E0_func_of_s


!----- for debuging --------------------------------------------------
      integer   :: err
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'ChannelNumber_AT_TS'
!-----------------------------------------------------------
  basisnD => para_H%para_AllBasis%BasisnD

  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*) 'Ene',Ene
    write(out_unitp,*) 'EneTS',para_CRP%Channel_AT_TS%EneTS
    write(out_unitp,*)
    IF (allocated(BasisnD%EneH0)) THEN
      write(out_unitp,*) 'size BasisnD%EneH0',size(BasisnD%EneH0)
      write(out_unitp,*) 'BasisnD%EneH0',BasisnD%EneH0
    END IF
    CALL flush_perso(out_unitp)
  END IF
!-----------------------------------------------------------


  ny = para_H%mole%nb_act-1

  SELECT CASE (para_CRP%Channel_AT_TS%option)
  CASE (1) ! Lucien Dupuy, one degenerate frequency
    wy    = para_CRP%Channel_AT_TS%w1
    zpe   = HALF*ny*wy
    EneTS = para_CRP%Channel_AT_TS%EneTS
    nv    = int( (Ene-zpe-EneTS)/wy )

    nb_channels = 1
    IF ( nv > 0 ) THEN
      DO i=1,nv
        nb_channels = nb_channels + combination(i,ny)
      END DO
    END IF

  CASE (2) ! ny frequencies at TS with basis set (SG4)
    allocate(nbSize(para_H%para_AllBasis%BasisnD%nb_basis-1))
    allocate(tab_ib(para_H%para_AllBasis%BasisnD%nb_basis-1))
    EneTS = para_CRP%Channel_AT_TS%EneTS

    LB = para_H%para_AllBasis%BasisnD%L_SparseBasis

    DO ib=2,para_H%para_AllBasis%BasisnD%nb_basis
      nbSize(ib-1) = para_H%para_AllBasis%BasisnD%tab_basisPrimSG(LB,ib)%nb
    END DO


    allocate(tab_i_TO_l(para_H%para_AllBasis%BasisnD%nb_basis-1))
    DO ib=2,para_H%para_AllBasis%BasisnD%nb_basis
      nb = para_H%para_AllBasis%BasisnD%tab_basisPrimSG(LB,ib)%nb
      CALL alloc_dnSVM(tab_i_TO_l(ib-1),nb)
      IF (para_H%para_AllBasis%BasisnD%tab_basisPrimSG(LB,ib)%nb_basis < 1) THEN
        tab_i_TO_l(ib-1)%vec(:) = para_H%para_AllBasis%BasisnD%tab_basisPrimSG(LB,ib)%nDindB%Tab_L(:)
      ELSE
        DO i=1,nb
          n = para_H%para_AllBasis%BasisnD%tab_basisPrimSG(LB,ib)%nDindB%Tab_L(i)
          tab_i_TO_l(ib-1)%vec(i) = get_L_FROM_Basis_L_TO_n(para_H%para_AllBasis%BasisnD%tab_basisPrimSG(LB,ib)%L_TO_nb,n)
        END DO
      END IF
    END DO


    nDindB_Channels%packed = .TRUE. ! with false the mapping is too long !!
    CALL init_nDindexPrim(nDindB_Channels,ny,nbSize,                            &
                          type_OF_nDindex=5,Lmax=LB,tab_i_TO_l=tab_i_TO_l)

    nb_channels = para_CRP%Channel_AT_TS%nb_channels_added

    CALL init_nDval_OF_nDindex(nDindB_Channels,tab_ib)
    DO ib=1,nDindB_Channels%Max_nDI
      CALL ADD_ONE_TO_nDindex(nDindB_Channels,tab_ib,iG=ib)
      tab_ib(:) = tab_ib(:)-1
      EneChannel = sum(tab_ib(:)*para_CRP%Channel_AT_TS%w(:)) +                 &
                  HALF*sum(para_CRP%Channel_AT_TS%w)
      IF (debug) write(out_unitp,*) 'ib,tab_ib(:)-1',ib,tab_ib,' :',            &
                                    EneChannel,EneTS + EneChannel
      IF (Ene >= EneTS + EneChannel) nb_channels = nb_channels + 1
    END DO


    DO ib=1,size(tab_i_TO_l)
      CALL dealloc_dnSVM(tab_i_TO_l(ib))
    END DO
    deallocate(tab_i_TO_l)

    CALL dealloc_nDindex(nDindB_Channels)

    deallocate(nbSize)
    deallocate(tab_ib)

  CASE (3) ! With the energy of the "inactive" basis functions

    !write(out_unitp,*) 'size BasisnD%EneH0',size(BasisnD%EneH0)
    !write(out_unitp,*) 'BasisnD%EneH0',BasisnD%EneH0(:)
    IF (.NOT. allocated(BasisnD%EneH0))                                         &
        STOP 'ERROR in ChannelNumber_AT_TS: EneH0 is not allocated'


    ! first the energy of BasisnD%tab
    SELECT CASE (BasisnD%SparseGrid_type)
    CASE (0) ! Direct product
      E0_func_of_s = BasisnD%tab_Pbasis(1)%Pbasis%EneH0(1)
      ! IF (debug) Then
      !   write(out_unitp,*) 'BasisnD%tab_Pbasis(1)%Pbasis%EneH0',BasisnD%tab_Pbasis(1)%Pbasis%EneH0
      !   write(out_unitp,*) 'BasisnD%tab_Pbasis(2)%Pbasis%EneH0',BasisnD%tab_Pbasis(2)%Pbasis%EneH0
      ! END IF
    CASE (1) ! Sparse basis
      E0_func_of_s = BasisnD%tab_basisPrimSG(1,BasisnD%L_SparseBasis)%EneH0(1)
    CASE (2,4) ! Sparse basis
      E0_func_of_s = BasisnD%tab_basisPrimSG(BasisnD%L_SparseBasis,1)%EneH0(1)
    END SELECT
    IF (debug) write(out_unitp,*) 'E0_func_of_s',E0_func_of_s

    allocate(nDval(BasisnD%nb_basis))

    EneTS = para_CRP%Channel_AT_TS%EneTS

    nb_channels = para_CRP%Channel_AT_TS%nb_channels_added

    IF (allocated(BasisnD%nDindB_contracted)) THEN
      CALL init_nDval_OF_nDindex(BasisnD%nDindB_contracted,nDval)
      DO ib=1,BasisnD%nDindB_contracted%Max_nDI
        CALL ADD_ONE_TO_nDindex(BasisnD%nDindB_contracted,nDval)
        IF (nDval(1) == 1) Then
          EneChannel = BasisnD%EneH0(ib) - E0_func_of_s
          IF (debug) THEN
            write(out_unitp,*) ib,'nDval',nDval
            write(out_unitp,*) ib,'EneChannel',EneChannel
          END IF
          write(out_unitp,*) 'Ene       ',Ene
          write(out_unitp,*) 'EneChannel',EneChannel

          IF (Ene >= EneChannel) nb_channels = nb_channels + 1
        END IF
      END DO
      nb_channels = min(nb_channels,BasisnD%nDindB_contracted%Max_nDI)
    ELSE
      CALL init_nDval_OF_nDindex(BasisnD%nDindB,nDval)
      DO ib=1,BasisnD%nDindB%Max_nDI
        CALL ADD_ONE_TO_nDindex(BasisnD%nDindB,nDval)
        IF (nDval(1) == 1) Then
          EneChannel = BasisnD%EneH0(ib) - E0_func_of_s
          IF (debug) THEN
            write(out_unitp,*) ib,'nDval',nDval
            write(out_unitp,*) ib,'EneChannel',EneChannel
          END IF
          write(out_unitp,*) 'Ene       ',Ene
          write(out_unitp,*) 'EneChannel',EneChannel

          IF (Ene >= EneChannel) nb_channels = nb_channels + 1
        END IF
      END DO
      nb_channels = min(nb_channels,BasisnD%nDindB%Max_nDI)
    END IF
    deallocate(nDval)

  END SELECT

  IF (print_level > 1 .OR. debug)  write(out_unitp,*) 'nb_channels',nb_channels
  IF (debug) THEN
    write(out_unitp,*)
    write(out_unitp,*) 'END ',name_sub
    CALL flush_perso(out_unitp)
  END IF

END FUNCTION ChannelNumber_AT_TS


END MODULE mod_CRP
