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

MODULE mod_fullanalysis
USE mod_Constant
CONTAINS
!================================================================
!
!     write the energy ene(i) and the vectors i psi(.,i)
!
!
!================================================================
      SUBROUTINE sub_analyse(Tab_Psi,nb_psi_in,para_H,                  &
                          para_ana,para_intensity,para_AllOp,const_phys)
      USE mod_system
      USE mod_Coord_KEO
      USE mod_basis
      USE mod_ana_psi
      USE mod_psi_set_alloc
      USE mod_psi_B_TO_G
      USE mod_psi_io
      USE mod_type_ana_psi
      USE mod_Op
      USE mod_analysis
      USE mod_MPI
      IMPLICIT NONE

!----- variables for the construction of H ----------------------------
      TYPE (param_AllOp)  :: para_AllOp
      TYPE (param_Op)     :: para_H

!----- variables pour la namelist analyse ----------------------------
      TYPE (param_ana)           :: para_ana
      TYPE (param_intensity)     :: para_intensity

!----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix),pointer     :: mole
      TYPE (Tnum),pointer        :: para_Tnum

!----- physical and mathematical constants ---------------------------
      TYPE (constant)            :: const_phys

!----- for the basis set ----------------------------------------------
      TYPE (Basis), pointer      :: BasisnD

!----- variables for the WP propagation ----------------------------
      integer            :: nb_psi_in
      TYPE (param_psi)   :: Tab_Psi(nb_psi_in)

!------ working variables ---------------------------------

      real (kind=Rkind), allocatable :: ene(:)
      real (kind=Rkind), allocatable :: Mat_psi(:,:)
      character(len=:),  allocatable :: info

      logical                        :: cube = .FALSE.

      TYPE (param_ana_psi)           :: ana_psi

      integer                        :: i,nb_col,ib
      real (kind=Rkind)              :: Q,E,DE
      TYPE (param_file)              :: file_WPspectral
      integer                        :: nioWP
      character (len=Name_longlen)   :: lformat
      TYPE(REAL_WU)                  :: RWU_ZPE,RWU_E,RWU_DE

      real (kind=Rkind), allocatable :: AllPsi_max_RedDensity(:)

      integer  :: Version_File,nb_psi,nb_tot
      namelist / headerFile / Version_File,nb_psi,nb_tot


!----- FUNCTION --------------------------------------------------
      real (kind=Rkind) :: part_func
!---- FUNCTION ---------------------------------------------------

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_analyse'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------

      mole       => para_H%mole
      para_Tnum  => para_H%para_Tnum
      BasisnD    => para_H%para_AllBasis%BasisnD

      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_ba,nb_qa',para_H%nb_ba,para_H%nb_qa
        write(out_unitp,*) 'nb_bi',para_H%nb_bi
        write(out_unitp,*) 'nb_bie',para_H%nb_bie
        write(out_unitp,*) 'nb_act1',mole%nb_act1
        write(out_unitp,*) 'max_ana,max_ene',para_ana%max_ana,para_ana%max_ene
        write(out_unitp,*) 'nb_psi_in',nb_psi_in

      END IF
!-----------------------------------------------------------
      CALL alloc_NParray(ene,shape(Tab_psi),'ene',name_sub)

      write(*,*) 'ene check1:',ene(1),Tab_psi(1)%CAvOp,size(ene),size(Tab_psi),' from ',MPI_id
      ene(:) = real(Tab_psi(:)%CAvOp,kind=Rkind)
      write(*,*) 'ene check2:',ene(1),Tab_psi(1)%CAvOp,size(ene),size(Tab_psi),' from ',MPI_id
      CALL Set_ZPE_OF_ComOp(para_H%ComOp,Ene=ene)
      write(*,*) 'ene check3:',ene(1),Tab_psi(1)%CAvOp,size(ene),size(Tab_psi),' from ',MPI_id
      IF (count(ene(:)-para_H%ComOp%ZPE <= para_ana%max_ene) == 0) RETURN

      RWU_ZPE = REAL_WU(para_H%ComOp%ZPE,'au','E')
      RWU_E   = REAL_WU(sum(ene) / real(nb_psi_in,kind=Rkind),'au','E')
      write(out_unitp,*) 'BEGINNING ',name_sub
      write(out_unitp,*)
      write(out_unitp,*) 'ZPE        : ',RWU_Write(RWU_ZPE,WithUnit=.TRUE. ,WorkingUnit=.FALSE.)
      write(out_unitp,*) 'Average_ene: ',RWU_Write(RWU_E,  WithUnit=.TRUE. ,WorkingUnit=.FALSE.)
      CALL flush_perso(out_unitp)
      IF (para_ana%max_ana > nb_psi_in) para_ana%max_ana = nb_psi_in

      IF (para_ana%intensity .AND. para_intensity%l_IntVR) THEN
        CALL alloc_NParray(para_intensity%ABC,(/3,nb_psi_in /),            &
                          "para_intensity%ABC",name_sub)
      END IF

       write(out_unitp,*)
       Q =  part_func(ene,nb_psi_in,para_ana%Temp,const_phys)

       write(out_unitp,*) 'population at T, Q',para_ana%Temp,Q
       write(out_unitp,*) 'Energy level (',const_phys%ene_unit,') pop and means :'
       CALL flush_perso(out_unitp)

      file_WPspectral%name = make_FileName(para_ana%name_file_spectralWP)
      CALL file_open(file_WPspectral,nioWP,lformatted=para_ana%formatted_file_WP)

      ! For the header of the file
      nb_psi        = count((ene(:)-para_H%ComOp%ZPE) <= para_ana%max_ene)
      CALL Write_header_saveFile_psi(tab_Psi,nb_psi,file_WPspectral)

      ! write the energy level + save the psi
      DO i=1,nb_psi_in

        IF (ene(i)-para_H%ComOp%ZPE > para_ana%max_ene) CYCLE

        ana_psi%Ene = ene(i)
        ana_psi%num_psi = i

        RWU_E  = REAL_WU(ana_psi%Ene,'au','E')
        RWU_DE = REAL_WU(ana_psi%Ene-para_H%ComOp%ZPE,'au','E')
        E  = convRWU_TO_R(RWU_E ,WorkingUnit=.FALSE.)
        DE = convRWU_TO_R(RWU_DE,WorkingUnit=.FALSE.)


        IF (ana_psi%num_psi < 10000) THEN
          lformat = '("lev0: ",i4,i4,l3,3(1x,' // trim(adjustl(EneIO_format)) // '))'
        ELSE
          lformat = '("lev0: ",i0,i0,l3,3(1x,' // trim(adjustl(EneIO_format)) // '))'
        END IF

        write(out_unitp,lformat) ana_psi%num_psi,0,tab_Psi(i)%convAvOp,E,DE

        CALL Write_Psi_nDBasis(tab_Psi(i),nioWP,i,ZERO,file_WPspectral%formatted,FilePsiVersion)

      END DO
      close(nioWP)

      ana_psi%ZPE        = para_H%ComOp%ZPE
      ana_psi%Part_Func  = Q
      ana_psi%Temp       = para_ana%Temp

      write(out_unitp,*) 'population at T, Q',para_ana%Temp,Q
      write(out_unitp,*) 'Energy level (',const_phys%ene_unit,') pop and means :'
      CALL flush_perso(out_unitp)

      DO i=1,nb_psi_in

        IF (ene(i)-para_H%ComOp%ZPE > para_ana%max_ene) CYCLE

        IF (.NOT. tab_Psi(i)%BasisRep) THEN
          CALL sub_PsiGridRep_TO_BasisRep(tab_Psi(i))
        END IF
        ana_psi%Ene     = ene(i)
        ana_psi%num_psi = i

        info = String_TO_String( " " //                                 &
           real_TO_char( ene(i)*const_phys%auTOenergy,"f12.6" ) // " : ")

        CALL sub_analyze_psi(tab_Psi(i),ana_psi)

        IF (allocated(ana_psi%max_RedDensity)) THEN
          IF (.NOT. allocated(AllPsi_max_RedDensity)) THEN
            CALL alloc_NParray(AllPsi_max_RedDensity,shape(ana_psi%max_RedDensity), &
                              "AllPsi_max_RedDensity",name_sub)
            AllPsi_max_RedDensity(:) = ZERO
          END IF

          DO ib=1,size(AllPsi_max_RedDensity)
            AllPsi_max_RedDensity(ib) = max(AllPsi_max_RedDensity(ib),ana_psi%max_RedDensity(ib))
          END DO
        END IF

        IF (para_ana%intensity .AND. para_intensity%l_IntVR) THEN
          CALL sub_moyABC(tab_Psi(i),i,info,para_intensity%ABC(:,i),para_AllOp)
        ELSE IF (para_AllOp%tab_Op(1)%para_PES%nb_scalar_Op > 0 .AND. ana_psi%AvScalOp) THEN
          CALL sub_moyScalOp(tab_Psi(i),i,info,para_AllOp)
        END IF


        write(out_unitp,*)

        deallocate(info)

      END DO

      IF (allocated(AllPsi_max_RedDensity)) THEN

        CALL Write_Vec(AllPsi_max_RedDensity,out_unitp,6,Rformat='e10.3',name_info='For all psi max_RedDensity ')
        !write(out_unitp,*) 'For all psi max_RedDensity ',AllPsi_max_RedDensity(:)
        CALL dealloc_NParray(AllPsi_max_RedDensity,"AllPsi_max_RedDensity",name_sub)
      END IF

      CALL flush_perso(out_unitp)

!----------------------------------------------------------


!----------------------------------------------------------
!     writing the eigenvectors
      para_ana%print_psi = min(para_ana%print_psi,nb_psi_in)
      IF (debug) para_ana%print_psi = nb_psi_in
      write(out_unitp,*) 'para_ana%print_psi',para_ana%print_psi
      CALL flush_perso(out_unitp)

      IF (cube) CALL write_cube(Tab_Psi)


      IF (para_ana%print_psi > 0 .OR. debug) THEN


        CALL alloc_NParray(Mat_psi,(/ tab_Psi(1)%nb_tot,para_ana%print_psi /), &
                        "Mat_psi",name_sub)
        DO i=1,para_ana%print_psi
          Mat_psi(:,i) = tab_Psi(i)%RvecB(:)
        END DO
        IF (.NOT. para_H%ComOp%contrac_ba_ON_HAC .AND. mole%nb_act1 < 3) THEN
          CALL write_psi(Mat_psi,para_ana%psi2,para_ana%print_psi,      &
                          tab_Psi(1)%nb_tot,                            &
                          para_H%nb_ba,para_H%nb_qa,para_H%nb_bie,       &
                          para_Tnum,mole,BasisnD,para_H)
        !ELSE IF (para_H%ComOp%contrac_ba_ON_HAC .AND. mole%nb_act1 < 3) THEN
        !  CALL write_psi2_new(tab_Psi)
        END IF
        CALL flush_perso(out_unitp)
        nb_col = 5
        write(out_unitp,*) 'eigenvectors in column'
        write(out_unitp,*) nb_col,para_ana%print_psi,tab_Psi(1)%nb_tot
        CALL Write_Mat(Mat_psi(:,1:para_ana%print_psi),out_unitp,nb_col)
        CALL dealloc_NParray(Mat_psi,"Mat_psi",name_sub)

      END IF
!----------------------------------------------------------

      CALL dealloc_ana_psi(ana_psi)
      CALL dealloc_NParray(ene,'ene',name_sub)

!----------------------------------------------------------
      IF (debug) THEN
      END IF
      write(out_unitp,*) 'END ',name_sub
!----------------------------------------------------------


      END SUBROUTINE sub_analyse
!================================================================
!
!     calculation of <psi | Mhu | psi>
!
!================================================================

!================================================================
!
!     calculation of <psi | Mhu | psi>
!
!================================================================
       SUBROUTINE sub_moyABC(Psi,iPsi,info,ABC,para_AllOp)

      USE mod_system
      USE mod_Op
      USE mod_psi_set_alloc
      IMPLICIT NONE

      TYPE (param_psi)               :: Psi
      integer                        :: iPsi
      character (len=*)              :: info
      real (kind=Rkind)              :: ABC(3)
      TYPE (param_AllOp)             :: para_AllOp



!----- for the zmatrix and Tnum --------------------------------------
       TYPE (param_psi)          :: OpPsi
      TYPE (zmatrix),pointer     :: mole      ! true pointer
      TYPE (Tnum),pointer        :: para_Tnum ! true pointer

      real (kind=Rkind) :: avMhu(3,3),TensorI(3,3),mat(3,3)
      real (kind=Rkind) :: trav1(3),mat1(3,3),dummy
      integer           :: index(3)
      complex(kind=Rkind) :: avOp
!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING sub_moyABC'
         write(out_unitp,*) 'ipsi,info',iPsi,info
       END IF
!-----------------------------------------------------------
      mole       => para_AllOp%tab_Op(1)%mole
      para_Tnum  => para_AllOp%tab_Op(1)%para_Tnum


       OpPsi = psi  ! for the initialization

       CALL sub_PsiOpPsi(avOp,Psi,OpPsi,para_AllOp%tab_Op(6))
       TensorI(1,1) = real(avOp,kind=Rkind)
       CALL sub_PsiOpPsi(avOp,Psi,OpPsi,para_AllOp%tab_Op(7))
       TensorI(1,2) = real(avOp,kind=Rkind)
       TensorI(2,1) = real(avOp,kind=Rkind)
       CALL sub_PsiOpPsi(avOp,Psi,OpPsi,para_AllOp%tab_Op(8))
       TensorI(1,3) = real(avOp,kind=Rkind)
       TensorI(3,1) = real(avOp,kind=Rkind)
       CALL sub_PsiOpPsi(avOp,Psi,OpPsi,para_AllOp%tab_Op(10))
       TensorI(2,2) = real(avOp,kind=Rkind)
       CALL sub_PsiOpPsi(avOp,Psi,OpPsi,para_AllOp%tab_Op(11))
       TensorI(2,3) = real(avOp,kind=Rkind)
       TensorI(3,2) = real(avOp,kind=Rkind)
       CALL sub_PsiOpPsi(avOp,Psi,OpPsi,para_AllOp%tab_Op(14))
       TensorI(3,3) = real(avOp,kind=Rkind)

       IF (para_Tnum%Inertia) THEN
         write(out_unitp,*) iPsi,'TensorI',info
         write(out_unitp,"(3(3(f15.6,1x),/))") TensorI(:,:)
         CALL inversion(avMhu,TensorI,trav1,index,3)
         avMhu = avMhu * HALF
       ELSE
         avMhu = TensorI
         mat   = TensorI

         CALL inversion(TensorI,mat,trav1,index,3)
         write(out_unitp,*) iPsi,'TensorI (invers of G)',info
         write(out_unitp,"(3(3(f15.6,1x),/))") TensorI(:,:)
       END IF

       write(out_unitp,*) iPsi,'avMhu',info
       write(out_unitp,"(3(3(f15.9,1x),/))") avMhu(:,:)

       CALL diagonalization(avMhu,ABC(:),mat1,3,1,1,.FALSE.)
       dummy = ABC(1)
       ABC(1) = ABC(3)
       ABC(3) = dummy


       write(out_unitp,21) iPsi,' ABC (cm-1) at ',info,ABC(:) *         &
                                          get_Conv_au_TO_unit('E','cm-1')
       write(out_unitp,21) iPsi,' ABC (GHz) at ',info,ABC(:) *          &
                                          get_Conv_au_TO_unit('E','GHz')
 21    format(i4,2A,3f18.5)

       CALL dealloc_psi(OpPsi)


!----------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'END sub_moyABC'
        END IF
!----------------------------------------------------------


        end subroutine sub_moyABC

!================================================================
!
!     calculation of <psi | Mhu | psi>
!
!================================================================
       SUBROUTINE sub_moyScalOp(Psi,iPsi,info,para_AllOp)

      USE mod_system
      USE mod_Op
      USE mod_psi_set_alloc
      IMPLICIT NONE

      TYPE (param_psi)               :: Psi
      integer                        :: iPsi
      character (len=*)              :: info
      TYPE (param_AllOp)             :: para_AllOp



      !----- local variables --------------------------------------
      TYPE (param_psi)           :: OpPsi
      TYPE (param_Op), pointer   :: ScalOp(:) => null() ! true pointer


      real (kind=Rkind)   :: avScalOp(para_AllOp%tab_Op(1)%para_PES%nb_scalar_Op)
      complex(kind=Rkind) :: avOp
      integer             :: iOp,nb_scalar_Op
!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
       nb_scalar_Op = para_AllOp%tab_Op(1)%para_PES%nb_scalar_Op
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING sub_moyScalOp'
         write(out_unitp,*) 'ipsi,info',iPsi,info
         write(out_unitp,*) 'nb_scalar_Op',nb_scalar_Op
         DO iOp=1,size(para_AllOp%tab_Op)
           write(out_unitp,*) iOp,'Save_MemGrid_done', &
              para_AllOp%tab_Op(iOp)%para_ReadOp%para_FileGrid%Save_MemGrid_done
         END DO
       END IF
!-----------------------------------------------------------

      ScalOp(1:nb_scalar_Op) => para_AllOp%tab_Op(3:2+nb_scalar_Op)

       OpPsi = psi  ! for the initialization

       DO iOp=1,nb_scalar_Op

         CALL sub_PsiOpPsi(avOp,Psi,OpPsi,ScalOp(iOp))
         avScalOp(iOp) = real(avOp,kind=Rkind)

       END DO

       write(out_unitp,"(i0,2a,100(f15.9,1x))") iPsi,' avScalOp: ',info,avScalOp

       CALL dealloc_psi(OpPsi)
       nullify(ScalOp)


!----------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'END sub_moyScalOp'
        END IF
!----------------------------------------------------------


        end subroutine sub_moyScalOp

!================================================================
!
!     write psi or psi^2 on the grid point
!
!================================================================
      SUBROUTINE write_psi(Mat_psi,psi2,nb_psi,nb_tot,                  &
                           nb_ba,nb_qa,n_h,                             &
                           para_Tnum,mole,                              &
                           BasisnD,para_Op)

      USE mod_system
      USE mod_Coord_KEO
      USE mod_basis
      USE mod_Op
      IMPLICIT NONE


!----- for the zmatrix and Tnum --------------------------------------
      TYPE (zmatrix) :: mole
      TYPE (Tnum)    :: para_Tnum

!----- variables for the construction of H ----------------------------
      TYPE (param_Op)   :: para_Op


      integer           :: nb_ba,nb_qa,n_h
      integer           :: nb_psi,nb_tot
      real (kind=Rkind) :: Mat_psi(nb_tot,nb_psi)
      logical           :: psi2

!----- for the basis set ----------------------------------------------
      TYPE (Basis) :: BasisnD


!------ working variables ---------------------------------
      integer       :: i,k,l,i_q
      integer       :: ih,ihk
      real (kind=Rkind), allocatable :: Qact1(:)
      real (kind=Rkind), allocatable :: psi_q(:,:)
      integer       :: ih_print
      integer, parameter :: max_print = 200

      real (kind=Rkind), allocatable :: d0b(:)
      real (kind=Rkind), allocatable :: psid0b_k(:)

      real (kind=Rkind)              :: WrhonD
      real (kind=Rkind), save        :: Q1


!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING write_psi'
        write(out_unitp,*) 'nb_psi,psi2',nb_psi,psi2
        write(out_unitp,*) 'nb_ba',nb_ba
        write(out_unitp,*) 'nb_qa',nb_qa
        write(out_unitp,*) 'n_h',n_h
      END IF
!-----------------------------------------------------------

      IF (BasisnD%cplx) THEN
         write(out_unitp,*) ' ERROR in write_psi'
         write(out_unitp,*) ' the basis is complex'
         STOP
      END IF


!-----------------------------------------------------------
      write(out_unitp,*)
      write(out_unitp,*) 'eigenvectors on a grid',nb_psi

!     - initisalisation ----------------------------------
      CALL alloc_NParray(d0b,     (/ nb_ba /),       'd0b','write_psi')
      CALL alloc_NParray(psid0b_k,(/ nb_ba /),       'psid0b_k','write_psi')
      CALL alloc_NParray(Qact1,   (/ mole%nb_act1 /),'Qact1','write_psi')
      CALL alloc_NParray(psi_q,   (/ nb_psi,n_h /),  'psi_q','write_psi')


!      - check the phase of psi(:,i) ----------
       DO i=1,nb_psi
         IF (Mat_psi(1,i) < ZERO) Mat_psi(:,i) = -Mat_psi(:,i)
       END DO


       DO i_q=1,nb_qa

        CALL Rec_Qact(Qact1,BasisnD,i_q,mole)

        IF (i_q == 1) Q1 = Qact1(1)

!       - calculation of WrhonD ------------------------------
         WrhonD = Rec_WrhonD(BasisnD,i_q)
         CALL calc_d0b(d0b,BasisnD,i_q)

         DO i=1,nb_psi

           psi_q(i,:) = ZERO
           IF (psi2) THEN


             DO ih=1,n_h
               ihk = (ih-1)*nb_ba
               DO k=1,nb_ba
                 ihk = ihk + 1
                 psid0b_k(k) = Mat_psi(ihk,i) * d0b(k)
               END DO

               DO k=1,nb_ba
               DO l=1,nb_ba
                 psi_q(i,1) = psi_q(i,1) + psid0b_k(k) * psid0b_k(l)
               END DO
               END DO

             END DO


           ELSE


             DO ih=1,n_h
               ihk = (ih-1)*nb_ba
               psi_q(i,ih)=dot_product(Mat_psi(ihk+1:ihk+nb_ba,i),d0b(:))
             END DO

           END IF

         END DO

         ih_print = n_h
         IF (psi2) ih_print = 1

          IF (mole%nb_act1 == 2 .AND. abs(Q1-Qact1(1)) > ONETENTH**6) THEN
             write(out_unitp,*)
             Q1 = Qact1(1)
          END IF

         write(out_unitp,31) Qact1(:),WrhonD,psi_q(1:min(max_print,nb_psi),1:ih_print)
 31      format(3f20.10,200f20.10)

       END DO

      CALL dealloc_NParray(d0b,     'd0b',     'write_psi')
      CALL dealloc_NParray(psid0b_k,'psid0b_k','write_psi')
      CALL dealloc_NParray(Qact1,   'Qact1',   'write_psi')
      CALL dealloc_NParray(psi_q,   'psi_q',   'write_psi')
!----------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'END write_psi'
        END IF
!----------------------------------------------------------


      end subroutine write_psi

      SUBROUTINE write_psi2_new(Tab_Psi)

      USE mod_system
      USE mod_ana_psi
      USE mod_psi_set_alloc
      USE mod_psi_B_TO_G
      IMPLICIT NONE


      TYPE (param_psi)   :: Tab_Psi(:)


      integer           :: nb_psi


!------ working variables ---------------------------------
      integer       :: i,i_q,ie,ieq,nb_bi,nb_qa
      real (kind=Rkind), allocatable :: Qact1(:)
      real (kind=Rkind), allocatable :: psi_q(:,:)
      integer, parameter :: max_print = 200

      real (kind=Rkind)              :: WrhonD

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='write_psi2_new'
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      nb_psi = size(Tab_Psi)
      IF (nb_psi < 1) RETURN

      nb_bi  = Tab_Psi(1)%nb_bi
      nb_qa  = Tab_Psi(1)%nb_qa
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_psi',nb_psi
      END IF
!-----------------------------------------------------------
      write(out_unitp,*)
      write(out_unitp,*) 'eigenvectors on a grid',nb_psi


      CALL alloc_NParray(psi_q,(/nb_qa,nb_psi/),'psi_q',name_sub)
      CALL alloc_NParray(Qact1,(/Tab_Psi(1)%nb_act1/),'Qact1',name_sub)


!-----------------------------------------------------------


       DO i=1,nb_psi
         CALL sub_PsiBasisRep_TO_GridRep(Tab_Psi(i))

         DO i_q=1,nb_qa
           psi_q(i_q,i)  = ZERO
           DO ie=1,Tab_Psi(i)%nb_be*nb_bi
             ieq = (ie-1)*nb_qa + i_q
             IF (Tab_Psi(i)%cplx) THEN
               psi_q(i_q,i) = psi_q(i_q,i) + abs(Tab_Psi(i)%CvecG(ieq))**2
             ELSE
               psi_q(i_q,i) = psi_q(i_q,i) + Tab_Psi(i)%RvecG(ieq)**2
             END IF
           END DO

         END DO
       END DO


       DO i_q=1,nb_qa

         CALL Rec_x(Qact1,Tab_Psi(1)%BasisnD,i_q)

         !- calculation of WrhonD ------------------------------
         WrhonD = Rec_WrhonD(Tab_Psi(1)%BasisnD,i_q)

         write(out_unitp,31) Qact1(:),WrhonD,psi_q(i_q,1:nb_psi)
 31      format(3f20.10,200f20.10)

       END DO

      CALL dealloc_NParray(psi_q,'psi_q',name_sub)
      CALL dealloc_NParray(Qact1,'Qact1',name_sub)


!----------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'END ',name_sub
        END IF
!----------------------------------------------------------


      end subroutine write_psi2_new

      SUBROUTINE write_cube(Tab_Psi)

      USE mod_system
      USE mod_ana_psi
      USE mod_psi_set_alloc
      USE mod_psi_B_TO_G
      IMPLICIT NONE


      TYPE (param_psi)   :: Tab_Psi(:)


      integer           :: nb_psi


!------ working variables ---------------------------------
      integer       :: i,i_q,ie,ieq,nb_bi,nb_qa
      real (kind=Rkind), allocatable :: psi_q(:,:)
      integer, parameter :: max_print = 200


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='write_cube'
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      nb_psi = size(Tab_Psi)
      IF (nb_psi < 1) RETURN

      nb_bi  = Tab_Psi(1)%nb_bi
      nb_qa  = Tab_Psi(1)%nb_qa
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_psi',nb_psi
      END IF
!-----------------------------------------------------------
      write(out_unitp,*)
      write(out_unitp,*) 'eigenvectors on a cube (test)',nb_psi

!-----------------------------------------------------------

      CALL alloc_NParray(psi_q,(/nb_qa,nb_psi/),'psi_q',name_sub)


       DO i=1,nb_psi
         CALL sub_PsiBasisRep_TO_GridRep(Tab_Psi(i))

         DO i_q=1,nb_qa
           psi_q(i_q,i)  = ZERO
           DO ie=1,Tab_Psi(i)%nb_be*nb_bi
             ieq = (ie-1)*nb_qa + i_q
             IF (Tab_Psi(i)%cplx) THEN
               psi_q(i_q,i) = psi_q(i_q,i) + abs(Tab_Psi(i)%CvecG(ieq))**2
             ELSE
               psi_q(i_q,i) = psi_q(i_q,i) + Tab_Psi(i)%RvecG(ieq)**2
             END IF
           END DO

         END DO
         !write cube file one for each Tab_Psi(i)
         write(100+i,*)  psi_q(:,i)
       END DO

      CALL dealloc_NParray(psi_q,'psi_q',name_sub)

!----------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'END ',name_sub
        END IF
!----------------------------------------------------------


      end subroutine write_cube

END MODULE mod_fullanalysis
