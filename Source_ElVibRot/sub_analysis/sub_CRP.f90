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

!================================================================
!     CRP
!================================================================
      SUBROUTINE sub_CRP(tab_Op,nb_Op,print_Op,CRPtype,CRP_Ene,CRP_DEne,nb_CRP_Ene)

      USE mod_system
      USE mod_Coord_KEO
      USE mod_basis
      USE mod_Op
      USE mod_analysis
      IMPLICIT NONE


!----- Operator variables ----------------------------------------------
      integer           :: nb_Op
      TYPE (param_Op)   :: tab_Op(nb_Op)
      logical           :: print_Op
      real (kind=Rkind) :: CRP_Ene,CRP_DEne
      integer           :: nb_CRP_Ene,CRPtype

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

      CALL sub_CRP_BasisRep_WithMat(tab_Op,nb_Op,print_Op,CRP_Ene,CRP_DEne,nb_CRP_Ene)

      !CALL sub_CRP_V0(tab_Op,nb_Op,print_Op,CRP_Ene)
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
      USE mod_analysis
      IMPLICIT NONE


!----- Operator variables ----------------------------------------------
      integer           :: nb_Op
      TYPE (param_Op)   :: tab_Op(nb_Op)
      logical           :: print_Op
      real (kind=Rkind) :: CRP_Ene,CRP_DEne
      integer           :: nb_CRP_Ene

!---- variable for the Z-matrix ----------------------------------------
      TYPE (zmatrix), pointer    :: mole
      TYPE (Tnum), pointer       :: para_Tnum


!----- working variables -----------------------------
      integer       ::    i,k,ie

      complex (kind=Rkind), allocatable :: G(:,:)
      complex (kind=Rkind), allocatable :: Ginv(:,:)
      complex (kind=Rkind), allocatable :: gGgG(:,:)
      complex :: CRP
      real (kind=Rkind) :: Ene
      real (kind=Rkind) :: CRP_Eckart ! function



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
        !write(6,*) 'id diff ?',maxval(abs(Ginv))


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
      FUNCTION CRP_Eckart(E)
      USE mod_system
      IMPLICIT NONE
      real (kind=Rkind) :: CRP_Eckart
      real (kind=Rkind) :: E

      real (kind=Rkind) :: b,c
      real (kind=Rkind), parameter :: V0=0.0156_Rkind,m=1060._Rkind,a=ONE

       b = a * Pi * sqrt(TWO*m*E)
       c = (Pi/TWO) * sqrt(EIGHT * V0*m*a**2 - 1)

       CRP_Eckart = ONE / (ONE + (cosh(c)/sinh(b))**2)

      END FUNCTION CRP_Eckart

      SUBROUTINE sub_CRP_V0(tab_Op,nb_Op,print_Op,Ene_CRP)

      USE mod_system
      USE mod_Coord_KEO
      USE mod_basis
      USE mod_Op
      USE mod_analysis
      IMPLICIT NONE


!----- Operator variables ----------------------------------------------
      integer           :: nb_Op
      TYPE (param_Op)   :: tab_Op(nb_Op)
      logical           :: print_Op
      real (kind=Rkind) :: Ene_CRP


!---- variable for the Z-matrix ----------------------------------------
      TYPE (zmatrix), pointer    :: mole
      TYPE (Tnum), pointer       :: para_Tnum


!----- working variables -----------------------------
      integer       ::    i,k,ie

      complex (kind=Rkind), allocatable :: G(:,:)
      complex (kind=Rkind), allocatable :: Ginv(:,:)
      complex (kind=Rkind), allocatable :: gGgG(:,:)
      complex :: CRP
      real (kind=Rkind) :: Ene0_CRP,Dene_CRP


!----- for debuging --------------------------------------------------
      integer   :: err
      logical, parameter :: debug=.FALSE.
!      logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'sub_CRP_V0'
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
          Ginv(i,i) = Ginv(i,i) + Ene_CRP
      END DO
      Dene_CRP = ONETENTH**4

      DO ie=0,100

        DO i=1,tab_Op(1)%nb_tot
          Ginv(i,i) = Ginv(i,i) + DEne_CRP
        END DO

        CALL inv_m1_TO_m2_cplx(Ginv,G,tab_Op(1)%nb_tot,0,ZERO)
        !Ginv = matmul(Ginv,G)
        !DO i=1,tab_Op(1)%nb_tot
        !  Ginv(i,i) = Ginv(i,i) - CONE
        !END DO
        !write(6,*) 'id diff ?',maxval(abs(Ginv))


        gGgG(:,:) = matmul(tab_Op(3)%Rmat,matmul(G,matmul(tab_Op(4)%Rmat,conjg(G))))

        CRP = ZERO
        DO i=1,tab_Op(1)%nb_tot
          CRP = CRP + gGgG(i,i)
        END DO
        write(out_unitp,*) 'CRP at E (ua)',Ene_CRP+real(ie,kind=Rkind)*DEne_CRP,CRP

        !gGgG(:,:) = matmul(tab_Op(3)%Rmat,matmul(G,matmul(tab_Op(3)%Rmat,conjg(G))))
        !CRP = ZERO
        !DO i=1,tab_Op(1)%nb_tot
        !  CRP = CRP + gGgG(i,i)
        !END DO
        !write(out_unitp,*) 'CRP at E (ua)',Ene_CRP,CRP
        !
        !gGgG(:,:) = matmul(tab_Op(4)%Rmat,matmul(G,matmul(tab_Op(4)%Rmat,conjg(G))))
        !CRP = ZERO
        !DO i=1,tab_Op(1)%nb_tot
        !  CRP = CRP + gGgG(i,i)
        !END DO
        !write(out_unitp,*) 'CRP at E (ua)',Ene_CRP,CRP
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

      end subroutine sub_CRP_V0

