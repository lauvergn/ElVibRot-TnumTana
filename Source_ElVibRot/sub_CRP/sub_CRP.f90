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

  TYPE param_CRP
    real (kind=Rkind) :: Ene    = ZERO            ! Total energy for CRP
    real (kind=Rkind) :: DEne   = ZERO            ! Energy increment for the CRP
    integer           :: nb_Ene = 1               ! Number of CRP calculation

    character (len=Name_len) :: CRP_Type         = 'lanczos'
    integer                  :: KS_max_it        = 100
    real (kind=Rkind)        :: KS_accuracy      = ONETENTH**5
    character (len=Name_len) :: LinSolv_Type     = 'MatInv'
    integer                  :: LinSolv_max_it   = 100
    real (kind=Rkind)        :: LinSolv_accuracy = ONETENTH**7
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

  character (len=Name_len) :: CRP_Type         = 'lanczos'
  integer                  :: KS_max_it        = 100
  real (kind=Rkind)        :: KS_accuracy      = ONETENTH**5
  character (len=Name_len) :: LinSolv_Type     = 'MatInv'
  integer                  :: LinSolv_max_it   = 100
  real (kind=Rkind)        :: LinSolv_accuracy = ONETENTH**7

  !----- for debuging --------------------------------------------------
  integer :: err_mem,memory
  character (len=*), parameter :: name_sub = "read_CRP"
  !logical, parameter :: debug=.FALSE.
  logical, parameter :: debug=.TRUE.
  !-----------------------------------------------------------

  NAMELIST /CRP/Ene,DEne,nb_Ene,CRP_Type,                               &
                KS_max_it,KS_accuracy,                                  &
                LinSolv_Type,LinSolv_max_it,LinSolv_accuracy

  Ene              = REAL_WU(ZERO,'cm-1','E')
  DEne             = REAL_WU(ZERO,'cm-1','E')
  nb_Ene           = 1

  CRP_Type         = 'lanczos'
  KS_max_it        = 100
  KS_accuracy      = ONETENTH**5
  LinSolv_Type     = 'QMR'
  LinSolv_max_it   = 100
  LinSolv_accuracy = ONETENTH**7

  read(in_unitp,CRP)
  write(out_unitp,CRP)

  CALL string_uppercase_TO_lowercase(CRP_Type)
  CALL string_uppercase_TO_lowercase(LinSolv_Type)

  IF (print_level > 0) write(out_unitp,CRP)
  write(out_unitp,*)

  para_CRP%Ene              = convRWU_TO_R(Ene)
  para_CRP%DEne             = convRWU_TO_R(DEne)
  para_CRP%nb_Ene           = nb_Ene
  para_CRP%CRP_Type         = CRP_Type
  para_CRP%KS_max_it        = KS_max_it
  para_CRP%KS_accuracy      = KS_accuracy
  para_CRP%LinSolv_Type     = LinSolv_Type
  para_CRP%LinSolv_max_it   = LinSolv_max_it
  para_CRP%LinSolv_accuracy = LinSolv_accuracy

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
      SELECT CASE (para_CRP%CRP_type)
      CASE ('withmat') ! old one
        CALL sub_CRP_BasisRep_WithMat(tab_Op,nb_Op,print_Op,para_CRP%Ene,para_CRP%DEne,para_CRP%nb_Ene)

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
      SUBROUTINE sub_CRP_BasisRep_WithMat(tab_Op,nb_Op,print_Op,CRP_Ene,CRP_DEne,nb_CRP_Ene)

      USE mod_system
      USE mod_Coord_KEO
      USE mod_basis
      USE mod_Op
      IMPLICIT NONE


!----- Operator variables ----------------------------------------------
      integer           :: nb_Op
      TYPE (param_Op)   :: tab_Op(nb_Op)
      logical           :: print_Op
      real (kind=Rkind) :: CRP_Ene,CRP_DEne
      integer           :: nb_CRP_Ene

!---- variable for the Z-matrix ----------------------------------------
      TYPE (CoordType), pointer  :: mole
      TYPE (Tnum), pointer       :: para_Tnum


!----- working variables -----------------------------
      integer       ::    i,k,ie

      complex (kind=Rkind), allocatable :: G(:,:)
      complex (kind=Rkind), allocatable :: Ginv(:,:)
      complex (kind=Rkind), allocatable :: gGgG(:,:)
      complex :: CRP
      real (kind=Rkind) :: Ene



!----- for debuging --------------------------------------------------
      integer   :: err
      logical, parameter :: debug=.FALSE.
!      logical, parameter :: debug=.TRUE.
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

      write(out_unitp,*) 'shape H',shape(tab_Op(1)%Rmat)
      write(out_unitp,*) 'nb_tot of H',tab_Op(1)%nb_tot
      CALL flush_perso(out_unitp)

      DO i=1,nb_Op
        IF (i == 2) CYCLE
        CALL sub_MatOp(tab_Op(i),print_Op)
      END DO

      CALL alloc_NParray(Ginv,shape(tab_Op(1)%Rmat),'Ginv',name_sub)
      CALL alloc_NParray(G,shape(tab_Op(1)%Rmat),'G',name_sub)
      CALL alloc_NParray(gGgG,shape(tab_Op(1)%Rmat),'gGgG',name_sub)

      write(out_unitp,*) 'Ginv calc'
      Ginv(:,:) = -tab_Op(1)%Rmat + EYE*HALF * (tab_Op(3)%Rmat+tab_Op(4)%Rmat)

      DO i=1,tab_Op(1)%nb_tot
          Ginv(i,i) = Ginv(i,i) + CRP_Ene-CRP_DEne
      END DO
      Ene = CRP_Ene-CRP_DEne

      DO ie=0,nb_CRP_Ene-1

        DO i=1,tab_Op(1)%nb_tot
          Ginv(i,i) = Ginv(i,i) + CRP_DEne
        END DO
        Ene = Ene + CRP_DEne

        CALL inv_m1_TO_m2_cplx(Ginv,G,tab_Op(1)%nb_tot,0,ZERO)
        !Ginv = matmul(Ginv,G)
        !DO i=1,tab_Op(1)%nb_tot
        !  Ginv(i,i) = Ginv(i,i) - CONE
        !END DO
        !write(out_unitp,*) 'id diff ?',maxval(abs(Ginv))


        gGgG(:,:) = matmul(tab_Op(3)%Rmat,matmul(G,matmul(tab_Op(4)%Rmat,conjg(G))))

        CRP = ZERO
        DO i=1,tab_Op(1)%nb_tot
          CRP = CRP + gGgG(i,i)
        END DO
        write(out_unitp,*) 'CRP at E (ua)',CRP_Ene+real(ie,kind=Rkind)*CRP_DEne,&
                                real(CRP,kind=Rkind),aimag(CRP),CRP_Eckart(Ene)
        !write(out_unitp,*) 'CRP at E (ua)',Ene,CRP
      END DO


      CALL flush_perso(out_unitp)
      CALL dealloc_NParray(Ginv,'Ginv',name_sub)
      CALL dealloc_NParray(gGgG,'gGgG',name_sub)
      CALL dealloc_NParray(G,'G',name_sub)
!----------------------------------------------------------
      IF (debug) THEN
      END IF
      write(out_unitp,*) 'END ',name_sub
!----------------------------------------------------------

END SUBROUTINE sub_CRP_BasisRep_WithMat
SUBROUTINE calc_crp_p_lanczos(tab_Op, nb_Op, para_CRP,Ene)
      USE mod_Op
      implicit none

!----- Operator variables ----------------------------------------------
      integer,            intent(in)      :: nb_Op
      TYPE (param_Op),    intent(inout)   :: tab_Op(nb_Op)

      TYPE (param_CRP),   intent(in)      :: para_CRP
      real(kind=Rkind),   intent(in)      :: Ene

!----- working variables -----------------------------

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


      integer       ::    k,ie
      complex (kind=Rkind), allocatable :: G(:,:)
      complex (kind=Rkind), allocatable :: Ginv(:,:)
      complex (kind=Rkind), allocatable :: gGgG(:,:)
      real (kind=Rkind) :: CRP, oldcrp,crp2
      complex(kind=Rkind):: M1(tab_Op(1)%nb_tot)

      character (len=*), parameter :: name_sub = 'calc_crp_p_lanczos'


      ncooked   = tab_Op(1)%nb_tot


    ! If need be, generate explicit matrix representation of operators
    if ( Inv ) then

         DO i=1,nb_Op
            IF (i == 2) CYCLE
            CALL sub_MatOp(tab_Op(i),print_Op=.FALSE.)
         END DO

         CALL alloc_NParray(Ginv,shape(tab_Op(1)%Rmat),'Ginv',name_sub)
         CALL alloc_NParray(G,   shape(tab_Op(1)%Rmat),'G',   name_sub)
         CALL alloc_NParray(gGgG,shape(tab_Op(1)%Rmat),'gGgG',name_sub)

         Ginv(:,:) = - tab_Op(1)%Rmat + EYE*HALF * (tab_Op(3)%Rmat+tab_Op(4)%Rmat)

         DO i=1,tab_Op(1)%nb_tot
            Ginv(i,i) = Ginv(i,i) + Ene
         END DO

         CALL inv_m1_TO_m2_cplx(Ginv,G,tab_Op(1)%nb_tot,0,ZERO)

         gGgG(:,:) = matmul(tab_Op(3)%Rmat,matmul(G,matmul(tab_Op(4)%Rmat,conjg(G))))


         crp2 = ZERO
         do i=1,tab_Op(1)%nb_tot
            crp2 = crp2 + real( gGgG(i,i), kind=Rkind)
         end do

         CALL dealloc_NParray(Ginv,'Ginv',name_sub)
         CALL dealloc_NParray(gGgG,'gGgG',name_sub)
         CALL dealloc_NParray(G,'G',name_sub)

         write(out_unitp,*) 'CRP at E (ua)', Ene, crp,'CRP with explicit inversion =', crp2
    ELSE

!      IF (allocated(tab_Op(1)%BasisnD%EneH0)) THEN
!        M1(:) = ONE/(Ene-tab_Op(1)%BasisnD%EneH0(:)) ! approximation of 1/(Ene-H(i,i))
!        !M1(:) = (Ene-tab_Op(1)%BasisnD%EneH0(:))
!        write(out_unitp,*) 'precon /= 1. DML'
!      ELSE
!        M1(:)        = CONE
!        write(out_unitp,*) 'precon = 1. DML'
!      END IF
      M1(:)        = CONE
      write(out_unitp,*) 'precon = 1. DML'

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
      !write(out_unitp,*) 'Krylov_vectors(:,0)',Krylov_vectors(:,0) ; stop

      IF (para_CRP%LinSolv_type == 'matinv') THEN
         DO i=1,nb_Op
            IF (i == 2) CYCLE
            CALL sub_MatOp(tab_Op(i),print_Op=.FALSE.)
         END DO

         CALL alloc_NParray(Ginv,shape(tab_Op(1)%Rmat),'Ginv',name_sub)
         CALL alloc_NParray(G,   shape(tab_Op(1)%Rmat),'G',   name_sub)
         CALL alloc_NParray(gGgG,shape(tab_Op(1)%Rmat),'gGgG',name_sub)

         Ginv(:,:) = - tab_Op(1)%Rmat + EYE*HALF * (tab_Op(3)%Rmat+tab_Op(4)%Rmat)

         DO i=1,tab_Op(1)%nb_tot
            Ginv(i,i) = Ginv(i,i) + Ene
         END DO

         CALL inv_m1_TO_m2_cplx(Ginv,G,tab_Op(1)%nb_tot,0,ZERO)

         gGgG(:,:) = matmul(tab_Op(3)%Rmat,matmul(G,matmul(tab_Op(4)%Rmat,conjg(G))))

         CALL dealloc_NParray(Ginv,'Ginv',name_sub)
         CALL dealloc_NParray(G,   'G',   name_sub)

      END IF

      ! Begin Lanczos scheme
      oldcrp = ZERO
      do nks=1,para_CRP%KS_max_it

         write(out_unitp,*) '######################'
         write(out_unitp,*) '# in KS iterations, n=',nks
         write(out_unitp,*) '# before p_multiply'
         call flush_perso(out_unitp)

         SELECT CASE ( para_CRP%LinSolv_type )

         CASE ( 'matinv' )

            Krylov_vectors(:,nks) = matmul(gGgG,Krylov_vectors(:,nks-1))

         CASE ( 'qmr' )

            call p_multiplyQMR(Krylov_vectors(:,nks-1),Krylov_vectors(:,nks), &
                               tab_Op,nb_Op,Ene,ncooked,M1,para_CRP%LinSolv_accuracy)

         CASE ( 'gmres' )
#if __CERFACS == 1
            call p_multiplyGMRES(Krylov_vectors(:,nks-1),Krylov_vectors(:,nks),&
                                 tab_Op,nb_Op,Ene,ncooked,M1,para_CRP%LinSolv_accuracy)
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

            if (dabs(oldcrp-crp) < para_CRP%KS_accuracy) EXIT
            oldcrp = crp
         endif

      end do

      !actual_iterations = nks
      write(out_unitp,*) 'CRP at E (ua)', Ene, crp

    end if

    IF (allocated(gGgG))           CALL dealloc_NParray(gGgG,'gGgG',name_sub)

    IF (allocated(Krylov_vectors)) CALL dealloc_NParray(Krylov_vectors,'Krylov_vectors',name_sub)
    IF (allocated(h))              CALL dealloc_NParray(h,'h',name_sub)


END SUBROUTINE calc_crp_p_lanczos
SUBROUTINE Gpsi(Vect,tab_Op,nb_Op,Ene,l_conjg)
      use mod_system
      USE mod_psi_set_alloc
      USE mod_Op
      implicit none

      real(kind=Rkind)  :: Ene
      integer           :: nb_Op
      TYPE (param_Op)   :: tab_Op(nb_Op)
      logical           :: print_Op
      !logical, parameter:: cplx=.TRUE.
      TYPE (param_psi)   :: Tab_Psi
      TYPE (param_psi)   :: Tab_OpPsi
      character(len=3) :: l_conjg
      complex(kind=Rkind), dimension(tab_Op(1)%nb_tot) :: Vect
      integer         :: i

      if (l_conjg == 'CJG') then
         Vect(:)=conjg(Vect(:))
      end if

      CALL init_psi(Tab_Psi,tab_Op(1),cplx=.TRUE.)
      CALL alloc_psi(Tab_Psi,BasisRep=.TRUE.,GridRep=.FALSE.)

      Tab_Psi%cvecB(:)=Vect(:)
      Tab_OpPsi = Tab_Psi

      call sub_OpPsi(Tab_Psi,Tab_OpPsi,tab_Op(1))

      Vect(:)= Vect(:)*Ene - Tab_OpPsi%cvecB(:)


      call sub_OpPsi(Tab_Psi,Tab_OpPsi,tab_Op(4))

      Vect(:) = Vect(:)+EYE*HALF*Tab_OpPsi%cvecB(:)

      call sub_OpPsi(Tab_Psi,Tab_OpPsi,tab_Op(3))

      Vect(:) = Vect(:)+EYE*HALF*Tab_OpPsi%cvecB(:)

      if (l_conjg == 'CJG') then
         Vect(:)=conjg(Vect(:))
      end if

      call dealloc_psi(Tab_Psi, .TRUE.)
      call dealloc_psi(Tab_OpPsi, .TRUE.)
END SUBROUTINE

SUBROUTINE OpOnVec(Vect,tab_Op,l_conjg)
      use mod_system
      USE mod_psi_set_alloc
      USE mod_Op
      implicit none

      TYPE (param_Op)   :: tab_Op
      logical           :: print_Op
      logical, parameter:: cplx=.TRUE.
      TYPE (param_psi)   :: Tab_Psi
      TYPE (param_psi)   :: Tab_OpPsi
      character(len=3) :: l_conjg
      complex(kind=Rkind), dimension(tab_Op%nb_tot) :: Vect
      integer         :: i

      if (l_conjg == 'CJG') then
         Vect(:)=conjg(Vect(:))
      end if

      CALL init_psi(Tab_Psi,tab_Op,cplx)
      CALL alloc_psi(Tab_Psi,BasisRep=.TRUE.,GridRep=.FALSE.)

      Tab_Psi%cvecB(:)=Vect(:)
      Tab_OpPsi = Tab_Psi

      call sub_OpPsi(Tab_Psi,Tab_OpPsi,tab_Op)

      Vect(:) = Tab_OpPsi%cvecB(:)

      if (l_conjg == 'CJG') then
         Vect(:)=conjg(Vect(:))
      end if

      call dealloc_psi(Tab_Psi, .TRUE.)
      call dealloc_psi(Tab_OpPsi, .TRUE.)
END SUBROUTINE OpOnVec

SUBROUTINE ReNorm_CplxVec(Vect)
      use mod_system
      implicit none

      complex (kind=Rkind), intent(inout) :: Vect(:)

      Vect(:) = Vect(:)/sqrt(dot_product(Vect,Vect))

END SUBROUTINE ReNorm_CplxVec
FUNCTION CRP_Eckart(E)
      USE mod_system
      IMPLICIT NONE
      real (kind=Rkind) :: CRP_Eckart
      real (kind=Rkind) :: E

      real (kind=Rkind) :: b,c
      real (kind=Rkind), parameter :: V0=0.0156_Rkind,m=1061._Rkind,a=ONE

       b = a * Pi * sqrt(TWO*m*E)
       c = (Pi/TWO) * sqrt(EIGHT * V0*m*a**2 - 1)

       CRP_Eckart = ONE / (ONE + (cosh(c)/sinh(b))**2)

END FUNCTION CRP_Eckart

END MODULE mod_CRP
