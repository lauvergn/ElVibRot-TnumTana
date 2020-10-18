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
    TYPE (CRP_Eckart_t)      :: Eckart

  END TYPE param_CRP

CONTAINS
SUBROUTINE read_CRP(para_CRP)
USE mod_system
USE mod_Constant
IMPLICIT NONE

  !----- variables pour la namelist analyse ----------------------------
  TYPE (param_CRP),     intent(inout)  :: para_CRP


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
                Eckart,With_Eckart

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

  IF (debug) write(out_unitp,*) 'E,DE,nb_E   : ',para_CRP%Ene,para_CRP%DEne,para_CRP%nb_Ene
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


      integer              :: i
      real (kind=Rkind)    :: Ene

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

      SELECT CASE (para_CRP%CRP_type)
      CASE ('withmat') ! old one
        CALL sub_CRP_BasisRep_WithMat(tab_Op,nb_Op,print_Op,para_CRP)

      CASE ('withmat_flux')
        CALL sub_CRP_BasisRep_WithMat_flux(tab_Op,nb_Op,print_Op,para_CRP)

      CASE ('lanczos') ! lanczos (Lucien Dupuy)
        DO i = 0, para_CRP%nb_Ene-1
          Ene = para_CRP%Ene+real(i,kind=Rkind)*para_CRP%DEne
          CALL calc_crp_P_lanczos(tab_Op, nb_Op,para_CRP,Ene)
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
      integer       ::    i,k,ie
      real (kind=Rkind), allocatable :: EneH(:),Vec(:,:) ! for debuging


      complex (kind=Rkind), allocatable :: G(:,:)
      complex (kind=Rkind), allocatable :: Ginv(:,:)
      complex (kind=Rkind), allocatable :: gGgG(:,:)
      complex :: CRP
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

      CALL sub_MatOp(tab_Op(1),print_Op) ! H
      DO i=3,nb_Op ! for the CAP
        IF (i == para_CRP%iOp_CAP_Reactif .OR. i == para_CRP%iOp_CAP_Product)   &
           CALL sub_MatOp(tab_Op(i),print_Op)
      END DO

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

        CALL inv_m1_TO_m2_cplx(Ginv,G,tab_Op(1)%nb_tot,0,ZERO)
        !Ginv = matmul(Ginv,G)
        !DO i=1,tab_Op(1)%nb_tot
        !  Ginv(i,i) = Ginv(i,i) - CONE
        !END DO
        !write(out_unitp,*) 'id diff ?',maxval(abs(Ginv))

        gGgG(:,:) = matmul(tab_Op(para_CRP%iOp_CAP_Reactif)%Rmat,               &
           matmul(G,matmul(tab_Op(para_CRP%iOp_CAP_Product)%Rmat,conjg(G))))

        CRP = ZERO
        DO i=1,tab_Op(1)%nb_tot
          CRP = CRP + gGgG(i,i)
        END DO

        RWU_E  = REAL_WU(Ene,'au','E')

        if (para_CRP%With_Eckart) then
          write(out_unitp,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                            real(CRP,kind=Rkind),aimag(CRP),CRP_Eckart(Ene,para_CRP%Eckart)
        else
          write(out_unitp,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                            real(CRP,kind=Rkind),aimag(CRP)
        end if

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
      complex :: CRP
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

      write(out_unitp,*) 'Op name: ',tab_Op(1)%name_Op
      CALL sub_MatOp(tab_Op(1),print_Op) ! H

      CALL flush_perso(out_unitp)

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
          CALL FluxOp_Mat(tab_Op(1),tab_Op(i),mEYE_FluxOpReactif_mat)
        END IF
        IF (i == para_CRP%iOp_Flux_Product) Then
          write(out_unitp,*) 'Op name: ',tab_Op(i)%name_Op
          CALL sub_MatOp(tab_Op(i),print_Op)
          CALL FluxOp_Mat(tab_Op(1),tab_Op(i),mEYE_FluxOpProduct_mat)
        END IF

        ! here the CAP
        IF (i == para_CRP%iOp_CAP_Reactif .OR. i == para_CRP%iOp_CAP_Product) THEN
          write(out_unitp,*) 'Op name: ',tab_Op(i)%name_Op
          CALL sub_MatOp(tab_Op(i),print_Op)
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
                            real(CRP,kind=Rkind),aimag(CRP),CRP_Eckart(Ene,para_CRP%Eckart)
        else
          write(out_unitp,*) 'CRP at ',RWU_Write(RWU_E,WithUnit=.TRUE.,WorkingUnit=.FALSE.),&
                            real(CRP,kind=Rkind),aimag(CRP)
        end if

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
SUBROUTINE calc_crp_p_lanczos(tab_Op,nb_Op,para_CRP,Ene)

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

      character (len=*), parameter :: name_sub = 'calc_crp_p_lanczos'


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

      CALL alloc_NParray(h,[para_CRP%KS_max_it,para_CRP%KS_max_it],'h',name_sub)

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

         ! Orthogonalize vectors
         do mks = 0, nks-1
            y = dot_product(Krylov_vectors(:,mks),Krylov_vectors(:,nks))
            Krylov_vectors(:,nks) = Krylov_vectors(:,nks) - y * Krylov_vectors(:,mks)
         end do
         ! Normalize vector
         CALL ReNorm_CplxVec(Krylov_vectors(:,nks))

         if (nks > 1) then
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

    IF (allocated(indx))           CALL dealloc_NParray(indx,'indx',name_sub)
    IF (allocated(Ginv))           CALL dealloc_NParray(Ginv,'Ginv',name_sub)


END SUBROUTINE calc_crp_p_lanczos

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



      real(kind=Rkind)        :: Rdiag(H%nb_tot),Rvp(H%nb_tot,H%nb_tot)


      FluxOp = matmul(H%Rmat,HStep_Op%Rmat) - matmul(HStep_Op%Rmat,H%Rmat)

      !FluxOp = ONE - HStep_Op%Rmat
      !FluxOp = HStep_Op%Rmat


      CALL  diagonalization(FluxOp,Rdiag,Rvp,H%nb_tot,4,1,.TRUE.)
      !write(out_unitp,*) 'Diag',Rdiag

END SUBROUTINE FluxOp_Mat

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

END MODULE mod_CRP
