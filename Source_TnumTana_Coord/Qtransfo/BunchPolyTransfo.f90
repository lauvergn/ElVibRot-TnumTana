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
      MODULE mod_BunchPolyTransfo
      use mod_system
      use mod_dnSVM ! only all
      use mod_constant,     only: table_atom, get_mass_tnum
      use mod_Lib_QTransfo, only: write_dnx, sub3_dnvec_toxf, func_ic
      USE mod_Tana_OpEl
      USE mod_Tana_Op1D
      USE mod_Tana_OpnD
      USE mod_Tana_sum_opnd
      USE mod_Tana_VecSumOpnD
      IMPLICIT NONE

      PRIVATE

      !!@description: TODO
      !!@param: TODO
      TYPE Type_BunchTransfo
          integer                    :: ncart=0,ncart_act=0
          integer                    :: nat0=0,nat=0,nat_act=0
          integer                    :: nb_var=0,nb_vect=0


          integer, pointer           :: ind_vect(:,:) => null()


          integer                    :: nb_X                   = 0       ! 0 (default) dummy atoms and centers of mass
          integer                    :: nb_G                   = 0       ! 0 (default) centers of mass

          real (kind=Rkind), pointer :: COM(:,:) => null()               ! COM(nat_act,nb_centers) defined all centers of mass as function active atom (not dummy)
          real (kind=Rkind), pointer :: Mat_At_TO_centers(:,:) => null() ! Mat_At_TO_centers(nat_act,nb_centers) (nb_centers=nat)
          real (kind=Rkind), pointer :: masses(:)              => null() ! masses (for the CC)
          real (kind=Rkind), pointer :: masses_OF_At(:)        => null() ! masses (for the atoms)

          real (kind=Rkind), pointer :: A(:,:)                 => null() ! Relations between the vectors and the X
          real (kind=Rkind), pointer :: A_inv(:,:)             => null() ! Relations between the vectors and the X
          real (kind=Rkind), pointer :: M_Tana(:,:)            => null() ! M for analytical KEO



          integer, pointer                  :: type_Qin(:)   => null() ! TRUE pointer
          character (len=Name_len), pointer :: name_Qin(:)   => null() ! TRUE pointer
          integer, pointer                  :: Z(:)          => null()
          character (len=Name_len),pointer  :: symbole(:)    => null()

      END TYPE Type_BunchTransfo

      !!@description: TODO
      !!@param: TODO
      TYPE Type_BFTransfo

          integer                  :: nb_var                 = 0
          integer                  :: nb_var_Rot             = 0

          integer                  :: nb_vect                = 0
          integer                  :: nb_vect_tot            = 0

          integer                  :: num_vect_in_Frame      = 0
          integer                  :: num_vect_in_BF         = 0

          integer                  :: num_Frame_in_BF        = 0
          integer                  :: num_Frame_in_Container = 0
          character (len=Name_len) :: name_Frame             = "F^(BF)"   ! BF
          integer, pointer         :: Tab_num_Frame(:)       => null()

          logical                    :: Frame                   = .FALSE.
          logical                    :: BF                      = .FALSE.
          integer                    :: Frame_type              = 0
          real (kind=Rkind), pointer :: Coef_Vect_FOR_xFrame(:) => null()
          real (kind=Rkind), pointer :: Coef_Vect_FOR_yFrame(:) => null()
          real (kind=Rkind), pointer :: Coef_Vect_FOR_zFrame(:) => null()
          real (kind=Rkind), pointer :: Coef_OF_Vect1(:)        => null()
          real (kind=Rkind), pointer :: Coef_OF_Vect2(:)        => null()
          character (len=3)          :: Type_Vect = 'zxy'
                 ! correspondance between Vec1 or vec2 and the BF axis
                 ! Vect1 => Type_Vect(1) axis (default z axis)
                 ! Vect2 => Type_Vect(2) axis (default x axis)
                 ! Vect1 ^ Vect2 => Type_Vect(3) axis (default y axis)


          logical                  :: cart                   = .FALSE.
          logical                  :: Li                     = .FALSE.
          logical                  :: Euler(3)               = (/ .FALSE., .FALSE., .FALSE. /)
                              ! F,F,F => for the true BF or F1
                              ! T,T,T => when a new BF is defined with 2 vectors
                              ! F,T,T => when a new BF is defined with ONE vector
                              ! Plus other cases
          logical                  :: cos_th                   = .TRUE.
          logical                  :: Def_cos_th               = .TRUE.

          integer                  :: iAtA=0,iAtB=0 ! enables to define the vector from 2 particles
          integer, pointer                  :: type_Qin(:) => null() ! TRUE pointer
          character (len=Name_len), pointer :: name_Qin(:) => null() ! TRUE pointer
          integer, pointer                  :: list_Qpoly_TO_Qprim(:) => null()
          integer, pointer                  :: list_Qprim_TO_Qpoly(:) => null()


          TYPE (Type_BFTransfo), pointer    :: tab_BFTransfo(:) => null() ! dim: nb_vect

          ! variables use for the calculation (Tana)
          ! They are defined here, because of the recursive structure
          type(opel)                     :: Qvec(3)   ! R,theta or u_theta, phi   or x,y,z
          type(opel)                     :: QEuler(3) ! alpha, beta or u_beta, gamma
          type(vec_sum_opnd)             :: Unit_Vector

          type(sum_opnd), pointer        :: M_mass(:,:)=>null()    ! mass matrix
          type(vec_sum_opnd)             :: J                      ! total angular momentum
          type(vec_sum_opnd)             :: Jdag                   ! adjoint of J
          integer,          pointer      :: listVFr(:) =>null()    ! index of the vector in the BF
          type(sum_opnd)                 :: KEO                    ! Output kEO
      END TYPE Type_BFTransfo

      INTERFACE alloc_array
        MODULE PROCEDURE alloc_array_OF_BFTransfodim1
      END INTERFACE
      INTERFACE dealloc_array
        MODULE PROCEDURE dealloc_array_OF_BFTransfodim1
      END INTERFACE

      PUBLIC :: Type_BunchTransfo, alloc_BunchTransfo, dealloc_BunchTransfo
      PUBLIC :: Read_BunchTransfo, Read2_BunchTransfo
      PUBLIC :: M_Tana_FROM_Bunch2Transfo, Write_BunchTransfo, calc_BunchTransfo
      PUBLIC :: BunchTransfo1TOBunchTransfo2

      PUBLIC :: Type_BFTransfo, RecRead_BFTransfo, RecWrite_BFTransfo
      PUBLIC :: dealloc_BFTransfo, alloc_array, dealloc_array
      PUBLIC :: calc_PolyTransfo, calc_PolyTransfo_outTOin, Rec_BFTransfo1TOBFTransfo2


      CONTAINS

      SUBROUTINE alloc_FrameType(BFTransfo,nb_vect)
      TYPE (Type_BFTransfo),intent(inout) :: BFTransfo

      integer, intent(in) :: nb_vect

      character (len=*), parameter :: name_sub='alloc_FrameType'

      !write(out_unitp,*) 'BEGINNING ',name_sub

      CALL dealloc_FrameType(BFTransfo)

      IF (nb_vect > 1) THEN

        CALL alloc_array(BFTransfo%Coef_Vect_FOR_xFrame,(/nb_vect/),    &
                        "BFTransfo%Coef_Vect_FOR_xFrame",name_sub)
        BFTransfo%Coef_Vect_FOR_xFrame(:) = 0

        CALL alloc_array(BFTransfo%Coef_Vect_FOR_yFrame,(/nb_vect/),    &
                        "BFTransfo%Coef_Vect_FOR_yFrame",name_sub)
        BFTransfo%Coef_Vect_FOR_yFrame(:) = 0

        CALL alloc_array(BFTransfo%Coef_Vect_FOR_zFrame,(/nb_vect/),    &
                        "BFTransfo%Coef_Vect_FOR_zFrame",name_sub)
        BFTransfo%Coef_Vect_FOR_zFrame(:) = 0

        CALL alloc_array(BFTransfo%Coef_OF_Vect1,(/nb_vect/),           &
                        "BFTransfo%Coef_OF_Vect1",name_sub)
        BFTransfo%Coef_OF_Vect1(:) = 0

        CALL alloc_array(BFTransfo%Coef_OF_Vect2,(/nb_vect/),           &
                        "BFTransfo%Coef_OF_Vect2",name_sub)
        BFTransfo%Coef_OF_Vect2(:) = 0
      END IF


      !write(out_unitp,*) 'END ',name_sub
      END SUBROUTINE alloc_FrameType

      SUBROUTINE dealloc_FrameType(BFTransfo)
      TYPE (Type_BFTransfo),intent(inout) :: BFTransfo

      character (len=*), parameter :: name_sub='dealloc_FrameType'

      !write(out_unitp,*) 'BEGINNING ',name_sub

      IF (associated(BFTransfo%Coef_Vect_FOR_xFrame)) THEN
        CALL dealloc_array(BFTransfo%Coef_Vect_FOR_xFrame,              &
                          "BFTransfo%Coef_Vect_FOR_xFrame",name_sub)
      END IF
      IF (associated(BFTransfo%Coef_Vect_FOR_yFrame)) THEN
        CALL dealloc_array(BFTransfo%Coef_Vect_FOR_yFrame,              &
                          "BFTransfo%Coef_Vect_FOR_yFrame",name_sub)

      END IF
      IF (associated(BFTransfo%Coef_Vect_FOR_zFrame)) THEN
        CALL dealloc_array(BFTransfo%Coef_Vect_FOR_zFrame,              &
                          "BFTransfo%Coef_Vect_FOR_zFrame",name_sub)
      END IF

      IF (associated(BFTransfo%Coef_OF_Vect1)) THEN
        CALL dealloc_array(BFTransfo%Coef_OF_Vect1,                     &
                          "BFTransfo%Coef_OF_Vect1",name_sub)
      END IF
      IF (associated(BFTransfo%Coef_OF_Vect2)) THEN
        CALL dealloc_array(BFTransfo%Coef_OF_Vect2,                     &
                          "BFTransfo%Coef_OF_Vect2",name_sub)
      END IF
      BFTransfo%Type_Vect = 'zxy'

      !write(out_unitp,*) 'END ',name_sub
      END SUBROUTINE dealloc_FrameType

      SUBROUTINE Write_FrameType(BFTransfo)
      TYPE (Type_BFTransfo),intent(in) :: BFTransfo


      character (len=*), parameter :: name_sub='Write_FrameType'

      !write(out_unitp,*) 'BEGINNING ',name_sub


      write(out_unitp,*) 'Frame,Frame_type',BFTransfo%Frame,BFTransfo%Frame_type
      IF (BFTransfo%Frame_type /= 0) THEN

        IF (associated(BFTransfo%Coef_Vect_FOR_xFrame)) THEN
            write(out_unitp,*) 'Coef_Vect_FOR_xFrame(:)',BFTransfo%Coef_Vect_FOR_xFrame(:)
        END IF
        IF (associated(BFTransfo%Coef_Vect_FOR_yFrame)) THEN
            write(out_unitp,*) 'Coef_Vect_FOR_yFrame(:)',BFTransfo%Coef_Vect_FOR_yFrame(:)
        END IF
        IF (associated(BFTransfo%Coef_Vect_FOR_zFrame)) THEN
            write(out_unitp,*) 'Coef_Vect_FOR_zFrame(:)',BFTransfo%Coef_Vect_FOR_zFrame(:)
        END IF

        IF (associated(BFTransfo%Coef_OF_Vect1)) THEN
          write(out_unitp,*) 'Coef_OF_Vect1(:)',BFTransfo%Coef_OF_Vect1(:)
        END IF
        IF (associated(BFTransfo%Coef_OF_Vect2)) THEN
          write(out_unitp,*) 'Coef_OF_Vect2(:)',BFTransfo%Coef_OF_Vect2(:)
        END IF
        write(out_unitp,*) 'Type_Vect: ',BFTransfo%Type_Vect
      END IF

      !write(out_unitp,*) 'END ',name_sub
      END SUBROUTINE Write_FrameType


      SUBROUTINE calc_Rot_Vect(tab_dnXVect,BFTransfo,nderiv)

      TYPE (Type_BFTransfo), intent(in)     :: BFTransfo
      TYPE (Type_dnVec),    intent(inout)   :: tab_dnXVect(:)
      integer, intent(in)                   :: nderiv

      TYPE (Type_dnS)   :: Rot(3,3),Vect1(3),Vect2(3)
      TYPE (Type_dnVec) :: dnVeczBF,dnVecxBF,dnVecyBF
      TYPE (Type_dnS)   :: dnSx,dnSy,dnSz


      integer           :: i,ixn,iyn,izn,iv,nb_var_deriv,nb_var_vec
      real (kind=Rkind) :: Coefx,Coefz
      logical           :: lx,ly,lz,l12
      logical           :: With_Rot = .FALSE.

!      -----------------------------------------------------------------
      integer :: nderiv_debug = 0
      logical, parameter :: debug = .FALSE.
!      logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='calc_Rot_Vect'
!      -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nderiv',nderiv
        CALL RecWrite_BFTransfo(BFTransfo,.FALSE.)

        DO iv=1,size(tab_dnXVect)
          write(out_unitp,*) 'tab_dnXVect(:)',iv
          CALL Write_dnVec(tab_dnXVect(iv))
        END DO

      END IF
!      -----------------------------------------------------------------

      lx = (sum(abs(BFTransfo%Coef_Vect_FOR_xFrame)) == 0)
      ly = (sum(abs(BFTransfo%Coef_Vect_FOR_yFrame)) == 0)
      lz = (sum(abs(BFTransfo%Coef_Vect_FOR_zFrame)) == 0)
      l12 = (sum(abs(BFTransfo%Coef_OF_Vect1)) == 0 .OR.                &
             sum(abs(BFTransfo%Coef_OF_Vect2)) == 0)
      IF (debug) write(out_unitp,*) 'lx,ly,lz,l12',lx,ly,lz,l12

      IF (.NOT. lx .AND. .NOT. lz .AND. ly .AND. l12) THEN ! with ZX
         CALL calc_Vect_zx(dnVecxBF,dnVecyBF,dnVeczBF,tab_dnXVect,BFTransfo,nderiv)
      ELSE IF (.NOT. lx .AND. lz .AND. ly .AND. .NOT. l12) THEN ! 12=>Z and X
         CALL calc_Vect_12zx(dnVecxBF,dnVecyBF,dnVeczBF,tab_dnXVect,BFTransfo,nderiv)
      ELSE ! not yet
        STOP 'calc_Rot_Vect not yet other cases!'
      END IF


      nb_var_deriv = tab_dnXVect(1)%nb_var_deriv
      nb_var_vec   = tab_dnXVect(1)%nb_var_vec

      CALL alloc_MatOFdnS(Rot,nb_var_deriv,nderiv)
      CALL alloc_VecOFdnS(Vect1,nb_var_deriv,nderiv)
      CALL alloc_VecOFdnS(Vect2,nb_var_deriv,nderiv)


      IF (With_Rot) THEN
        !Rotational matrix
        ixn=1
        iyn=2
        izn=3
        DO i=1,3
          CALL sub_dnVec_TO_dnS(dnVecxBF,Rot(ixn,i),i,nderiv)
          CALL sub_dnVec_TO_dnS(dnVecyBF,Rot(iyn,i),i,nderiv)
          CALL sub_dnVec_TO_dnS(dnVeczBF,Rot(izn,i),i,nderiv)
        END DO
        IF (debug) THEN
          write(out_unitp,*) 'Rotational matrix:'
          CALL Write_MatOFdnS(Rot,nderiv=0)
        END IF
        !Rotation to change the orientation of the frame
        DO iv=1,BFTransfo%nb_vect_tot

          IF (debug) write(out_unitp,*) 'Old Vector: ',iv
          IF (debug) CALL Write_dnVec(tab_dnXVect(iv),nderiv_debug)

          CALL sub_dnVec_TO_dnS(tab_dnXVect(iv),Vect1(1),1) ! x
          CALL sub_dnVec_TO_dnS(tab_dnXVect(iv),Vect1(2),2) ! y
          CALL sub_dnVec_TO_dnS(tab_dnXVect(iv),Vect1(3),3) ! z

          CALL Mat1OFdnS_MUL_Vec2OFdnS_TO_Vec3OFdnS(Rot,Vect1,Vect2,nderiv)

          CALL sub_dnS_TO_dnVec(Vect2(1),tab_dnXVect(iv),1) !x
          CALL sub_dnS_TO_dnVec(Vect2(2),tab_dnXVect(iv),2) !y
          CALL sub_dnS_TO_dnVec(Vect2(3),tab_dnXVect(iv),3) !z

          IF (debug) write(out_unitp,*) 'New Vector: ',iv
          IF (debug) CALL Write_dnVec(tab_dnXVect(iv),nderiv_debug)
        END DO
      ELSE
        DO iv=1,BFTransfo%nb_vect_tot
          CALL sub_dot_product_dnVec1_dnVec2_TO_dnS(dnVecxBF,tab_dnXVect(iv),dnSx,nderiv)
          CALL sub_dot_product_dnVec1_dnVec2_TO_dnS(dnVecyBF,tab_dnXVect(iv),dnSy,nderiv)
          CALL sub_dot_product_dnVec1_dnVec2_TO_dnS(dnVeczBF,tab_dnXVect(iv),dnSz,nderiv)
          IF (debug) write(out_unitp,*) 'New Vector: ',iv,dnSx%d0,dnSy%d0,dnSz%d0
          CALL sub_dnS_TO_dnVec(dnSx,tab_dnXVect(iv),1,nderiv)
          CALL sub_dnS_TO_dnVec(dnSy,tab_dnXVect(iv),2,nderiv)
          CALL sub_dnS_TO_dnVec(dnSz,tab_dnXVect(iv),3,nderiv)
      END DO
      END IF



      CALL dealloc_MatOFdnS(Rot)
      CALL dealloc_VecOFdnS(Vect1)
      CALL dealloc_VecOFdnS(Vect2)

      CALL dealloc_dnSVM(dnSx)
      CALL dealloc_dnSVM(dnSy)
      CALL dealloc_dnSVM(dnSz)

      CALL dealloc_dnSVM(dnVecxBF)
      CALL dealloc_dnSVM(dnVecyBF)
      CALL dealloc_dnSVM(dnVeczBF)


      !----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
        write(out_unitp,*)
      END IF


      END SUBROUTINE calc_Rot_Vect
      !  ez = V1 x V2
      !  ex ortho to Vx ez
      SUBROUTINE calc_Vect_12zx(dnVeczBF,dnVecxBF,dnVecyBF,tab_dnXVect,BFTransfo,nderiv)

      TYPE (Type_BFTransfo), intent(in)   :: BFTransfo
      TYPE (Type_dnVec),    intent(inout) :: tab_dnXVect(:)
      TYPE (Type_dnVec),    intent(inout) :: dnVeczBF,dnVecxBF,dnVecyBF

      integer, intent(in)                 :: nderiv


      TYPE (Type_dnVec) :: dnVec,dnVec1,dnVec2
      TYPE (Type_dnS)   :: dnNorm_inv
      integer           :: iv,nb_var_deriv,nb_var_vec
      real (kind=Rkind) :: Coefx,Coef1,Coef2

!      -----------------------------------------------------------------
      integer :: nderiv_debug = 0
      logical, parameter :: debug = .FALSE.
!      logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='calc_Vect_12zx'
!      -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nderiv',nderiv
        CALL RecWrite_BFTransfo(BFTransfo,.FALSE.)

        DO iv=1,size(tab_dnXVect)
          write(out_unitp,*) 'tab_dnXVect(:)',iv
          CALL Write_dnVec(tab_dnXVect(iv))
        END DO

      END IF
!      -----------------------------------------------------------------
      nb_var_deriv = tab_dnXVect(1)%nb_var_deriv
      nb_var_vec   = tab_dnXVect(1)%nb_var_vec

      CALL alloc_dnS(dnNorm_inv,nb_var_deriv,nderiv)
      CALL alloc_dnVec(dnVec,nb_var_vec,nb_var_deriv,nderiv)
      CALL alloc_dnVec(dnVec1,nb_var_vec,nb_var_deriv,nderiv)
      CALL alloc_dnVec(dnVec2,nb_var_vec,nb_var_deriv,nderiv)
      CALL alloc_dnVec(dnVecxBF,nb_var_vec,nb_var_deriv,nderiv)
      CALL alloc_dnVec(dnVecyBF,nb_var_vec,nb_var_deriv,nderiv)
      CALL alloc_dnVec(dnVeczBF,nb_var_vec,nb_var_deriv,nderiv)

      DO iv=1,size(tab_dnXVect)

         Coefx = BFTransfo%Coef_Vect_FOR_xFrame(iv)
         Coef1 = BFTransfo%Coef_OF_Vect1(iv)
         Coef2 = BFTransfo%Coef_OF_Vect2(iv)
         IF (debug) write(out_unitp,*) 'iv,Coef1,Coef2,Coefx',iv,Coef1,Coef2,Coefx

         IF (Coefx /= ZERO .OR. Coef1 /= ZERO .OR. Coef2 /= ZERO) THEN
           CALL sub_dnVec1_TO_dnVec2(tab_dnXVect(iv),dnVec,nderiv=nderiv)
           CALL sub_Normalize_dnVec(dnVec)

           IF (debug) CALL Write_dnVec(dnVec)

           IF (Coef1 /= ZERO) THEN
              CALL sub_dnVec1_wADDTO_dnVec2(dnVec,1,Coef1,dnVec1,1,ONE,3,nderiv=nderiv) ! add normalyzed vector
           END IF
           IF (Coef2 /= ZERO) THEN
              CALL sub_dnVec1_wADDTO_dnVec2(dnVec,1,Coef2,dnVec2,1,ONE,3,nderiv=nderiv) ! add normalyzed vector
           END IF
           IF (Coefx /= ZERO) THEN
              CALL sub_dnVec1_wADDTO_dnVec2(dnVec,1,Coefx,dnVecxBF,1,ONE,3,nderiv=nderiv) ! add normalyzed vector
           END IF
         END IF
         IF (debug) CALL Write_dnVec(dnVeczBF)
         IF (debug) CALL Write_dnVec(dnVecxBF)

      END DO

      ! new zBZ
      CALL Sub_crossproduct_dnVec1_dnVec2_TO_dnVec3(dnVec1,dnVec2,dnVeczBF,nderiv)
      CALL sub_Normalize_dnVec(dnVeczBF) ! normalyzed dnVeczBF

      IF (debug) write(out_unitp,*) 'New ez: '
      IF (debug) CALL Write_dnVec(dnVeczBF,nderiv=0)

      ! new xBZ
      ! orthogonalization of dnVecxBF
      CALL sub_dnVec1_TO_dnVec2(dnVeczBF,dnVec1,nderiv=nderiv)
      CALL sub_dot_product_dnVec1_dnVec2_TO_dnS(dnVecxBF,dnVec1,dnNorm_inv,nderiv)
      CALL sub_dnVec1_PROD_dnS2_TO_dnVec3(dnVec1,dnNorm_inv,dnVec1)
      CALL sub_dnVec1_wADDTO_dnVec2(dnVec1,1,-ONE,dnVecxBF,1,ONE,3,nderiv)
      CALL sub_Normalize_dnVec(dnVecxBF)

      IF (debug) write(out_unitp,*) 'New ex: '
      IF (debug) CALL Write_dnVec(dnVecxBF,nderiv=0)

      ! new yBZ
      CALL Sub_crossproduct_dnVec1_dnVec2_TO_dnVec3(dnVeczBF,dnVecxBF,dnVecyBF,nderiv)
      CALL sub_Normalize_dnVec(dnVecyBF)

      IF (debug) write(out_unitp,*) 'New ey: '
      IF (debug) CALL Write_dnVec(dnVecyBF,nderiv=0)


      CALL dealloc_dnSVM(dnNorm_inv)
      CALL dealloc_dnSVM(dnVec)
      CALL dealloc_dnSVM(dnVec1)
      CALL dealloc_dnSVM(dnVec2)


      !----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
        write(out_unitp,*)
      END IF


      END SUBROUTINE calc_Vect_12zx
      SUBROUTINE calc_Vect_zx(dnVecxBF,dnVecyBF,dnVeczBF,tab_dnXVect,BFTransfo,nderiv)

      TYPE (Type_BFTransfo), intent(in)   :: BFTransfo
      TYPE (Type_dnVec),    intent(inout) :: tab_dnXVect(:)
      TYPE (Type_dnVec),    intent(inout) :: dnVeczBF,dnVecxBF,dnVecyBF
      integer, intent(in)                 :: nderiv

      TYPE (Type_dnVec) :: dnVec1
      TYPE (Type_dnS)   :: dnNorm_inv
      integer           :: iv,nb_var_deriv,nb_var_vec
      real (kind=Rkind) :: Coefx,Coefz

!      -----------------------------------------------------------------
      integer :: nderiv_debug = 0
      logical, parameter :: debug = .FALSE.
!      logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='calc_Rot_Vect_zx'
!      -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nderiv',nderiv
        CALL RecWrite_BFTransfo(BFTransfo,.FALSE.)

        DO iv=1,size(tab_dnXVect)
          write(out_unitp,*) 'tab_dnXVect(:)',iv
          CALL Write_dnVec(tab_dnXVect(iv))
        END DO

      END IF
!      -----------------------------------------------------------------
      nb_var_deriv = tab_dnXVect(1)%nb_var_deriv
      nb_var_vec   = tab_dnXVect(1)%nb_var_vec


      CALL alloc_dnS(dnNorm_inv,nb_var_deriv,nderiv)
      CALL alloc_dnVec(dnVec1,nb_var_vec,nb_var_deriv,nderiv)
      CALL alloc_dnVec(dnVecxBF,nb_var_vec,nb_var_deriv,nderiv)
      CALL alloc_dnVec(dnVecyBF,nb_var_vec,nb_var_deriv,nderiv)
      CALL alloc_dnVec(dnVeczBF,nb_var_vec,nb_var_deriv,nderiv)

      DO iv=1,size(tab_dnXVect)
         Coefz = BFTransfo%Coef_Vect_FOR_zFrame(iv)
         Coefx = BFTransfo%Coef_Vect_FOR_xFrame(iv)
         IF (debug) write(out_unitp,*) 'iv,Coefz,Coefx',iv,Coefz,Coefx

         IF (Coefx /= ZERO .OR. Coefz /= ZERO) THEN
           CALL sub_dnVec1_TO_dnVec2(tab_dnXVect(iv),dnVec1,nderiv=nderiv)
           CALL sub_Normalize_dnVec(dnVec1)
           IF (debug) CALL Write_dnVec(dnVec1)

           IF (Coefz /= ZERO) THEN
              CALL sub_dnVec1_wADDTO_dnVec2(dnVec1,1,Coefz,dnVeczBF,1,ONE,3,nderiv=nderiv) ! add normalyzed vector
           END IF
           IF (Coefx /= ZERO) THEN
              CALL sub_dnVec1_wADDTO_dnVec2(dnVec1,1,Coefx,dnVecxBF,1,ONE,3,nderiv=nderiv) ! add normalyzed vector
           END IF
         END IF
         !CALL Write_dnVec(dnVeczBF)
         !CALL Write_dnVec(dnVecxBF)

      END DO

      ! new zBF
      CALL sub_Normalize_dnVec(dnVeczBF)

      IF (debug) write(out_unitp,*) 'New ez: ',dnVeczBF%d0
      IF (debug) CALL Write_dnVec(dnVeczBF,nderiv=0)

      ! new xBF
      ! orthogonalization of dnVecxBF
      CALL sub_dnVec1_TO_dnVec2(dnVeczBF,dnVec1,nderiv=nderiv)
      CALL sub_dot_product_dnVec1_dnVec2_TO_dnS(dnVecxBF,dnVec1,dnNorm_inv,nderiv)
      CALL sub_dnVec1_PROD_dnS2_TO_dnVec3(dnVec1,dnNorm_inv,dnVec1)
      CALL sub_dnVec1_wADDTO_dnVec2(dnVec1,1,-ONE,dnVecxBF,1,ONE,3,nderiv)
      CALL sub_Normalize_dnVec(dnVecxBF)

      IF (debug) write(out_unitp,*) 'New ex: ',dnVecxBF%d0
      IF (debug) CALL Write_dnVec(dnVecxBF,nderiv=0)

      ! new yBF
      CALL Sub_crossproduct_dnVec1_dnVec2_TO_dnVec3(dnVeczBF,dnVecxBF,dnVecyBF,nderiv)
      CALL sub_Normalize_dnVec(dnVecyBF)

      IF (debug) write(out_unitp,*) 'New ey: ',dnVecyBF%d0
      IF (debug) CALL Write_dnVec(dnVecyBF,nderiv=0)

      CALL dealloc_dnSVM(dnNorm_inv)
      CALL dealloc_dnSVM(dnVec1)


      !----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
        write(out_unitp,*)
      END IF


      END SUBROUTINE calc_Vect_zx

      SUBROUTINE FrameType1TOBFrameType2(BFTransfo1,BFTransfo2)
      TYPE (Type_BFTransfo),intent(in)    :: BFTransfo1
      TYPE (Type_BFTransfo),intent(inout) :: BFTransfo2

      integer :: nb_vect

!      -----------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!      logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='FrameType1TOBFrameType2'
!      -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL Write_FrameType(BFTransfo1)
        CALL flush_perso(out_unitp)
      END IF


      CALL dealloc_FrameType(BFTransfo2)

      BFTransfo2%Frame_type             = BFTransfo1%Frame_type

      IF (associated(BFTransfo1%Coef_Vect_FOR_xFrame)) THEN
        nb_vect = size(BFTransfo1%Coef_Vect_FOR_xFrame)
      ELSE
        nb_vect = 0
      END IF

      IF (nb_vect > 0) THEN
        CALL alloc_FrameType(BFTransfo2,nb_vect)
        BFTransfo2%Coef_Vect_FOR_xFrame = BFTransfo1%Coef_Vect_FOR_xFrame(:)
        BFTransfo2%Coef_Vect_FOR_yFrame = BFTransfo1%Coef_Vect_FOR_yFrame(:)
        BFTransfo2%Coef_Vect_FOR_zFrame = BFTransfo1%Coef_Vect_FOR_zFrame(:)

        BFTransfo2%Coef_OF_Vect1        = BFTransfo1%Coef_OF_Vect1(:)
        BFTransfo2%Coef_OF_Vect2        = BFTransfo1%Coef_OF_Vect2(:)
      END IF

      BFTransfo2%Type_Vect = BFTransfo1%Type_Vect

      IF (debug) THEN
        write(out_unitp,*)
        CALL Write_FrameType(BFTransfo2)
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF

      END SUBROUTINE FrameType1TOBFrameType2


      RECURSIVE SUBROUTINE dealloc_BFTransfo(BFTransfo)
      TYPE (Type_BFTransfo),intent(inout) :: BFTransfo

      integer :: iv,jv
      character (len=*), parameter :: name_sub='dealloc_BFTransfo'

      !write(out_unitp,*) 'BEGINNING ',name_sub

      ! the variables type_Qin and name_Qin are TRUE pointers.
      ! => they cannot be deallocated, but just nullified
      nullify(BFTransfo%type_Qin) ! true pointer
      nullify(BFTransfo%name_Qin) ! true pointer

      IF (associated(BFTransfo%tab_BFTransfo)) THEN
        DO iv=1,BFTransfo%nb_vect
          CALL dealloc_BFTransfo(BFTransfo%tab_BFTransfo(iv))
        END DO
        CALL dealloc_array(BFTransfo%tab_BFTransfo,                     &
                          "BFTransfo%tab_BFTransfo",name_sub)
      END IF

      IF (BFTransfo%BF) THEN
        IF (associated(BFTransfo%list_Qpoly_TO_Qprim)) THEN
          CALL dealloc_array(BFTransfo%list_Qpoly_TO_Qprim,             &
                            "BFTransfo%list_Qpoly_TO_Qprim",name_sub)
        END IF

        IF (associated(BFTransfo%list_Qprim_TO_Qpoly)) THEN
          CALL dealloc_array(BFTransfo%list_Qprim_TO_Qpoly,             &
                            "BFTransfo%list_Qprim_TO_Qpoly",name_sub)
        END IF
      END IF

      BFTransfo%nb_var                 = 0
      BFTransfo%nb_var_Rot             = 0
      BFTransfo%nb_vect                = 0
      BFTransfo%nb_vect_tot            = 0
      BFTransfo%num_Frame_in_BF        = 0
      BFTransfo%num_Frame_in_Container = 0
      BFTransfo%num_vect_in_Frame      = 0
      BFTransfo%num_vect_in_BF         = 0
      BFTransfo%name_Frame             = "F(BF)" ! BF
      BFTransfo%Frame                  = .FALSE.
      BFTransfo%BF                     = .FALSE.
      BFTransfo%Frame_type             = 0
      BFTransfo%cart                   = .FALSE.
      BFTransfo%Li                     = .FALSE.
      BFTransfo%Euler(:)               = (/ .FALSE., .FALSE., .FALSE. /)

      BFTransfo%cos_th                 = .TRUE.
      BFTransfo%Def_cos_th             = .TRUE.

      CALL dealloc_FrameType(BFTransfo)

      IF (associated(BFTransfo%tab_num_Frame))  THEN
        CALL dealloc_array(BFTransfo%tab_num_Frame,                     &
                          "BFTransfo%tab_num_Frame",name_sub)
      END IF

      CALL dealloc_TanaVar_FROM_BFTransfo(BFTransfo)

      !write(out_unitp,*) 'END ',name_sub
      END SUBROUTINE dealloc_BFTransfo

   RECURSIVE SUBROUTINE dealloc_TanaVar_FROM_BFTransfo(BFTransfo)
     type(Type_BFTransfo),            intent(inout)      :: BFTransfo

     integer                            :: i
     character (len=*), parameter       :: routine_name='dealloc_TanaVar_FROM_BFTransfo'


        IF (associated(BFTransfo%M_mass)) CALL dealloc_array(BFTransfo%M_mass,'BFTransfo%M_mass',routine_name)
        CALL delete_op(BFTransfo%J)
        CALL delete_op(BFTransfo%Jdag)
        IF (associated(BFTransfo%listVFr))                               &
                   CALL dealloc_array(BFTransfo%listVFr,'listVFr',routine_name)
        CALL delete_op(BFTransfo%KEO)
        CALL delete_op(BFTransfo%Unit_Vector)


       do i=1, BFTransfo%nb_vect
         CALL dealloc_TanaVar_FROM_BFTransfo(BFTransfo%tab_BFTransfo(i))
       end do

        CALL delete_op(BFTransfo%keo)


      BFTransfo%Qvec(1)   = czero
      BFTransfo%Qvec(2)   = czero
      BFTransfo%Qvec(3)   = czero
      BFTransfo%QEuler(1) = czero
      BFTransfo%QEuler(2) = czero
      BFTransfo%QEuler(3) = czero

   END SUBROUTINE dealloc_TanaVar_FROM_BFTransfo

      SUBROUTINE alloc_array_OF_BFTransfodim1(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      TYPE (Type_BFTransfo), pointer, intent(inout) :: tab(:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=1
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_BFTransfodim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (associated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       CALL sub_test_tab_ub(tab_ub,ndim,name_sub_alloc,name_var,name_sub)

       IF (present(tab_lb)) THEN
         CALL sub_test_tab_lb(tab_lb,ndim,name_sub_alloc,name_var,name_sub)

         memory = product(tab_ub(:)-tab_lb(:)+1)
         allocate(tab(tab_lb(1):tab_ub(1)),stat=err_mem)
       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Type_BFTransfo')

      END SUBROUTINE alloc_array_OF_BFTransfodim1
      SUBROUTINE dealloc_array_OF_BFTransfodim1(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (Type_BFTransfo), pointer, intent(inout) :: tab(:)
      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_BFTransfodim1'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN
       IF (.NOT. associated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Type_BFTransfo')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_BFTransfodim1

      !!@description: TODO
      !!@param: TODO
      RECURSIVE SUBROUTINE RecRead_BFTransfo(BFTransfo,BunchTransfo,    &
                                             num_Frame_in_Container)

      TYPE (Type_BFTransfo),intent(inout)    :: BFTransfo
      TYPE (Type_BunchTransfo),intent(inout) :: BunchTransfo

      integer, intent(inout)                 :: num_Frame_in_Container


      integer :: nb_vect,nb_var,iv,num_Frame_in_Container_rec,Frame_type
      logical :: Frame,cos_th,cart,zmat_order,Li
      character (len=Name_len) :: name_d,name_th,name_u,name_dih,        &
                                  name_x,name_y,name_z,                  &
                                  name_alpha,name_beta,name_gamma
      character (len=Name_len) :: name_Frame

      character(len=:), allocatable     :: name_F,name_v

      integer       :: iAtA,iAtB
      integer       :: i,nb_Qin,i_Qprim
      integer       :: iQalpha_TO_ivTot,iQbeta_TO_ivTot,iQgamma_TO_ivTot

      logical, save :: zmat_order_save        = .FALSE.
      integer, save :: i_Qpoly                = 0
      integer, save :: rec_level              = 0
      integer, save :: num_Frame_in_BF        = 0
      integer, save :: num_vect_in_BF         = 0
      integer, save :: iv_tot                 = 0

      integer, pointer :: tab_num_Frame(:)

      integer, parameter :: max_vect = 100
      real (kind=Rkind) :: Coef_Vect_FOR_xFrame(max_vect)
      real (kind=Rkind) :: Coef_Vect_FOR_yFrame(max_vect)
      real (kind=Rkind) :: Coef_Vect_FOR_zFrame(max_vect)
      real (kind=Rkind) :: Coef_OF_Vect1(max_vect)
      real (kind=Rkind) :: Coef_OF_Vect2(max_vect)
      character (len=3) :: Type_Vect

      NAMELIST /vector/ nb_vect,Frame,Frame_type,name_Frame,            &
                        zmat_order,                                     &
                        cos_th,name_d,name_th,name_u,name_dih,          &
                        cart,Li,name_x,name_y,name_z,                   &
                        name_alpha,name_beta,name_gamma,                &
                        iAtA,iAtB

      NAMELIST /Vect_FOR_AxisFrame / Coef_OF_Vect1,Coef_OF_Vect2,       &
                                     Type_Vect,Coef_Vect_FOR_xFrame,    &
                               Coef_Vect_FOR_yFrame,Coef_Vect_FOR_zFrame

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory,err_io
      character (len=*), parameter :: name_sub = "RecRead_BFTransfo"
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
       END IF
!-----------------------------------------------------------
      IF (BunchTransfo%nb_vect > max_vect) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) 'max_vect is too small',max_vect
        write(out_unitp,*) ' it MUST be larger than',BunchTransfo%nb_vect
        STOP
      END IF

      IF (iv_tot == 0) THEN
        nb_Qin = size(BFTransfo%type_Qin)
        CALL alloc_array(BFTransfo%list_Qpoly_TO_Qprim,(/nb_Qin/),      &
                        "BFTransfo%list_Qpoly_TO_Qprim",name_sub)
        BFTransfo%list_Qpoly_TO_Qprim(:) = 0

        CALL alloc_array(BFTransfo%list_Qprim_TO_Qpoly,(/nb_Qin/),      &
                        "BFTransfo%list_Qprim_TO_Qpoly",name_sub)
      END IF

      ! initialization for the namelist "vector"
      iv_tot        = iv_tot + 1
      iAtA          = 0
      iAtB          = 0
      nb_vect       = 0
      Frame         = .FALSE.
      Frame_type    = 0
      cos_th        = BFTransfo%Def_cos_th
      cart          = .FALSE.
      Li            = .FALSE.
      zmat_order    = .FALSE.
      name_x        = "x"
      name_y        = "y"
      name_z        = "z"
      name_d        = "R"
      name_th       = "th"
      name_u        = "u"
      name_dih      = "phi"
      name_alpha    = "alpha_"
      name_beta     = "beta_"
      name_gamma    = "gamma_"
      name_Frame    = ""

      read(in_unitp,vector,IOSTAT=err_io)
      IF (err_io < 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) '  while reading the namelist "vector"'
        write(out_unitp,*) ' end of file or end of record'
        write(out_unitp,*) ' Probably, nb_vect is to large ...'
        write(out_unitp,*) '   or you have forgotten the namelist.'
        write(out_unitp,*) ' Check your data !!'
        STOP
      END IF
      IF (err_io > 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) '  while reading the namelist "vector"'
        write(out_unitp,*) ' Probably, some arguments of namelist are wrong.'
        write(out_unitp,*) ' Check your data !!'
        STOP
      END IF

      IF (debug .OR. print_level > 1) write(out_unitp,vector)
      CALL flush_perso(out_unitp)

      IF (iv_tot == 1) zmat_order_save = zmat_order

      IF (BunchTransfo%ind_vect(3,iv_tot) == 0 .OR.                     &
          BunchTransfo%ind_vect(4,iv_tot) == 0) THEN
        BunchTransfo%ind_vect(3,iv_tot) = iAtA
        BunchTransfo%ind_vect(4,iv_tot) = iAtB
      END IF

      IF (len_trim(name_Frame) > 0) name_F = String_TO_string("  " // trim(name_Frame))
      BFTransfo%Frame             = Frame
      IF (.NOT. Frame .AND. Frame_type /= 0) THEN
        write(out_unitp,*) ' WARNING in ',name_sub
        write(out_unitp,*) ' Frame=F and Frame_type /= 0'
        write(out_unitp,*) ' Check your DATA !!'
      END IF
      IF (.NOT. Frame) Frame_type = 0
      IF (Frame_type /= 0) THEN
        Coef_Vect_FOR_xFrame(:) = ZERO
        Coef_Vect_FOR_yFrame(:) = ZERO
        Coef_Vect_FOR_zFrame(:) = ZERO
        Coef_OF_Vect1(:) = ZERO
        Coef_OF_Vect2(:) = ZERO
        Type_Vect        = 'zxy'

        read(in_unitp,Vect_FOR_AxisFrame,IOSTAT=err_io)
        IF (err_io < 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  while reading the namelist "Vect_FOR_AxisFrame"'
          write(out_unitp,*) ' end of file or end of record'
          write(out_unitp,*) ' Probably, you have forgotten the namelist.'
          write(out_unitp,*) ' Check your data !!'
          STOP
        END IF
        IF (err_io > 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  while reading the namelist "Vect_FOR_AxisFrame"'
          write(out_unitp,*) ' Probably, some arguments of namelist are wrong.'
          write(out_unitp,*) ' Check your data !!'
          STOP
        END IF

        CALL alloc_FrameType(BFTransfo,BunchTransfo%nb_vect)

        BFTransfo%Coef_Vect_FOR_xFrame = Coef_Vect_FOR_xFrame(1:BunchTransfo%nb_vect)
        BFTransfo%Coef_Vect_FOR_yFrame = Coef_Vect_FOR_yFrame(1:BunchTransfo%nb_vect)
        BFTransfo%Coef_Vect_FOR_zFrame = Coef_Vect_FOR_zFrame(1:BunchTransfo%nb_vect)

        BFTransfo%Coef_OF_Vect1 = Coef_OF_Vect1(1:BunchTransfo%nb_vect)
        BFTransfo%Coef_OF_Vect2 = Coef_OF_Vect2(1:BunchTransfo%nb_vect)

        CALL string_uppercase_TO_lowercase(Type_Vect)
        BFTransfo%Type_Vect = Type_Vect
        SELECT CASE(Type_Vect)
        CASE ("zxy","zyx","xyz","xzy","yzx","yxz")
          CONTINUE
        CASE default
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' Wrong "Type_Vect": ',Type_Vect
          write(out_unitp,*) ' The 6 possibilities are "zxy", "zyx", "xyz", "xzy", "yzx" and "yxz"'
          write(out_unitp,*) ' Check your data !!'
          STOP
        END SELECT

      END IF

      BFTransfo%cart              = cart
      BFTransfo%Li                = Li
      BFTransfo%nb_vect           = nb_vect
      num_vect_in_BF              = num_vect_in_BF + 1
      BFTransfo%num_vect_in_BF    = num_vect_in_BF
      BFTransfo%cos_th            = cos_th

      IF (nb_vect < 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' the number of vector is < 0!!',nb_vect
        STOP
      END IF

      IF (.NOT. BFTransfo%Frame .AND. BFTransfo%num_vect_in_BF == 1) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' The first vector MUST defined a Frame (BF)'
        write(out_unitp,*) ' Frame, num_vect_in_BF',Frame,BFTransfo%num_vect_in_BF
        STOP
      END IF
      IF (BFTransfo%Frame .AND. BFTransfo%num_vect_in_BF == 1) THEN
        BFTransfo%BF = .TRUE.
      ELSE
        BFTransfo%BF = .FALSE.
      END IF


      IF (BFTransfo%Frame) THEN
        BFTransfo%num_vect_in_Frame      = 1
        num_Frame_in_Container           = num_Frame_in_Container  + 1 ! recursive number of the frame
        rec_level                        = rec_level + 1               ! recursive index
        num_Frame_in_BF                  = num_Frame_in_BF + 1         ! recursive index
        BFTransfo%num_Frame_in_Container = num_Frame_in_Container
        BFTransfo%num_Frame_in_BF        = num_Frame_in_BF

        IF (.NOT. associated(BFTransfo%tab_num_Frame)) THEN ! rec_level MUST 1
          IF (rec_level /= 1) STOP 'WRONG rec_level'
          CALL alloc_array(BFTransfo%tab_num_Frame,(/rec_level/),       &
                          "BFTransfo%tab_num_Frame",name_sub)
          BFTransfo%tab_num_Frame(:) = 0
          BFTransfo%tab_num_Frame(rec_level) = 1
        ELSE ! the size of the table is not correct (rec_level-1))
          IF (size(BFTransfo%tab_num_Frame) /= rec_level-1) STOP 'PROBLEM with rec_level'

          nullify(tab_num_Frame)
          CALL alloc_array(tab_num_Frame,(/rec_level-1/),               &
                          "tab_num_Frame",name_sub)
          ! save BFTransfo%tab_num_Frame
          tab_num_Frame(:) = BFTransfo%tab_num_Frame(:)

          ! dealloc BFTransfo%tab_num_Frame
          CALL dealloc_array(BFTransfo%tab_num_Frame,                   &
                            "BFTransfo%tab_num_Frame",name_sub)

          ! alloc BFTransfo%tab_num_Frame with the new size
          CALL alloc_array(BFTransfo%tab_num_Frame,(/rec_level/),       &
                          "BFTransfo%tab_num_Frame",name_sub)

          ! set the new BFTransfo%tab_num_Frame
          BFTransfo%tab_num_Frame(:) = 0
          BFTransfo%tab_num_Frame(1:rec_level-1) = tab_num_Frame(:)
          BFTransfo%tab_num_Frame(rec_level)     = num_Frame_in_Container

          CALL dealloc_array(tab_num_Frame,"tab_num_Frame",name_sub)
        END IF

        ! Frame name
        IF (.NOT. allocated(name_F)) THEN
          name_F = String_TO_String("F(")
          DO i=rec_level,2,-1
            name_F = String_TO_String(name_F // int_TO_char(BFTransfo%tab_num_Frame(i)) // ",")
          END DO
          name_F = String_TO_String( name_F // "BF)" )
        END IF
        BFTransfo%name_Frame = name_F
        IF (allocated(name_F)) deallocate(name_F)

        name_v = String_TO_String( int_TO_char(BFTransfo%num_vect_in_Frame) // &
                                  "_" // trim(adjustl(BFTransfo%name_Frame)))


        IF (debug .OR. print_level > 1) THEN
          write(out_unitp,*) '.................................................'
          write(out_unitp,*) ' NEW Frame, rec_level: ',rec_level
          write(out_unitp,*) '   name_Frame,tab_num_Frame: ',trim(adjustl(BFTransfo%name_Frame)),' ',BFTransfo%tab_num_Frame(:)
          write(out_unitp,*) '   num_Frame_in_BF',BFTransfo%num_Frame_in_BF
          write(out_unitp,*) '   num_Frame_in_Container',BFTransfo%num_Frame_in_Container
          write(out_unitp,*)
          write(out_unitp,*) ' NEW vector (1st): '
          write(out_unitp,*) '   num_vect_in_Frame,num_vect_in_BF',BFTransfo%num_vect_in_Frame,BFTransfo%num_vect_in_BF
          write(out_unitp,*) '   name_v: ',trim(adjustl(name_v))
        END IF


        i_Qpoly = i_Qpoly + 1

        BFTransfo%Qvec(1) = set_opel(idf=1, idq=2, alfa=1, indexq=i_Qpoly, coeff=CONE) !R
        CALL get_unit_vector_Ei(BFTransfo%Unit_Vector, BFTransfo%Qvec, index_v=1, cart=cart)

        IF (zmat_order_save) THEN
          IF (iv_tot == 1) THEN
            BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly) = 1
          ELSE IF (iv_tot == 2) THEN
             BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly) = 2
          ELSE
            BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly) = (iv_tot-2)*3+1
          END IF
        ELSE
          BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly) = i_Qpoly
        END IF
        BFTransfo%list_Qprim_TO_Qpoly(BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)) = i_Qpoly

        i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)
        BFTransfo%type_Qin(i_Qprim) = 2 ! distance
        BFTransfo%name_Qin(i_Qprim) = "R" // trim(adjustl(name_v))
        BFTransfo%nb_vect_tot  = 1
        iQalpha_TO_ivTot  = iv_tot
        iQbeta_TO_ivTot   = iv_tot
        iQgamma_TO_ivTot  = iv_tot+1 ! If the vector "iv_tot+1" exist !!


        IF (allocated(name_v)) deallocate(name_v)

        IF (cart .OR. Li) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' The first vector of a frame CANNOT be defined'
          write(out_unitp,*) '   with Cartesian Coordinates or with Li'
          write(out_unitp,*) ' Check your data!'
          STOP
        END IF
        IF (nb_vect > 0) THEN
          CALL alloc_array(BFTransfo%tab_BFTransfo,(/nb_vect/),         &
                          "BFTransfo%tab_BFTransfo",name_sub)
          num_Frame_in_Container_rec = 0

          DO iv=1,nb_vect
            BFTransfo%tab_BFTransfo(iv)%num_vect_in_Frame      = iv+1
            BFTransfo%tab_BFTransfo(iv)%Frame                  = Frame
            BFTransfo%tab_BFTransfo(iv)%cart                   = cart
            BFTransfo%tab_BFTransfo(iv)%Li                     = Li
            BFTransfo%tab_BFTransfo(iv)%name_Frame             = BFTransfo%name_Frame
            BFTransfo%tab_BFTransfo(iv)%num_Frame_in_Container = num_Frame_in_Container
            BFTransfo%tab_BFTransfo(iv)%num_Frame_in_BF        = num_Frame_in_BF

            BFTransfo%tab_BFTransfo(iv)%euler(:)          = .TRUE.

            CALL alloc_array(BFTransfo%tab_BFTransfo(iv)%tab_num_Frame, &
                                                        (/rec_level/),  &
                            "BFTransfo%tab_BFTransfo(iv)%tab_num_Frame",name_sub)
            BFTransfo%tab_BFTransfo(iv)%tab_num_Frame(:) = BFTransfo%tab_num_Frame(:)

            IF (iv == 1) BFTransfo%tab_BFTransfo(iv)%euler(1) = .FALSE.

            BFTransfo%tab_BFTransfo(iv)%type_Qin => BFTransfo%type_Qin
            BFTransfo%tab_BFTransfo(iv)%name_Qin => BFTransfo%name_Qin

            BFTransfo%tab_BFTransfo(iv)%list_Qprim_TO_Qpoly => BFTransfo%list_Qprim_TO_Qpoly
            BFTransfo%tab_BFTransfo(iv)%list_Qpoly_TO_Qprim => BFTransfo%list_Qpoly_TO_Qprim

            BFTransfo%tab_BFTransfo(iv)%Def_cos_th = BFTransfo%Def_cos_th


            CALL RecRead_BFTransfo(BFTransfo%tab_BFTransfo(iv),         &
                                   BunchTransfo,                        &
                                   num_Frame_in_Container_rec)
            ! The next two lines HAVE to be after the "CALL RecRead_BFTransfo".
            BFTransfo%nb_vect_tot = BFTransfo%nb_vect_tot +             &
                                BFTransfo%tab_BFTransfo(iv)%nb_vect_tot

          END DO
        END IF
        IF (BFTransfo%nb_vect_tot == 1) BFTransfo%euler(3) = .FALSE.
        nb_var = max(1,3*BFTransfo%nb_vect_tot-3)
        BFTransfo%nb_var = nb_var
        !write(out_unitp,*) 'nb_var',nb_var

        BFTransfo%Frame_type     = 0
        IF (BFTransfo%euler(1) .AND. BFTransfo%euler(2) .AND. BFTransfo%euler(3)) THEN
          BFTransfo%Frame_type     = Frame_type
        END IF
        IF (.NOT.BFTransfo%euler(1) .AND. .NOT.BFTransfo%euler(2) .AND. .NOT.BFTransfo%euler(3)) THEN
          BFTransfo%Frame_type     = Frame_type
        END IF

        IF (BFTransfo%Frame_type == 0) CALL dealloc_FrameType(BFTransfo)

        ! Euler Angles
        IF (BFTransfo%euler(1)) THEN ! alpha
          i_Qpoly = i_Qpoly + 1
          BFTransfo%Qeuler(1) = set_opel(idf=1, idq=6, alfa=1, indexq=i_Qpoly, coeff=CONE)

          IF (zmat_order_save) THEN
            BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly) = (iQalpha_TO_ivTot-2)*3+3 ! like phi
          ELSE
            BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly) = i_Qpoly
          END IF
          BFTransfo%list_Qprim_TO_Qpoly(BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)) = i_Qpoly
          i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)
          BFTransfo%type_Qin(i_Qprim) = 4
          BFTransfo%name_Qin(i_Qprim) = trim(adjustl(name_alpha)) //         &
                 trim(adjustl(BFTransfo%name_Frame))
        END IF
        IF (BFTransfo%euler(2)) THEN ! beta
          i_Qpoly = i_Qpoly + 1
          IF (zmat_order_save) THEN
            IF (iQbeta_TO_ivTot == 2) THEN
              BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly) = 3
            ELSE
              BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly) = (iQbeta_TO_ivTot-2)*3+2! like theta
            END IF
          ELSE
            BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly) = i_Qpoly
          END IF
          BFTransfo%list_Qprim_TO_Qpoly(BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)) = i_Qpoly
          i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)
          IF (cos_th) THEN
            BFTransfo%Qeuler(2) = set_opel(idf=1, idq=-7, alfa=1, indexq=i_Qpoly, coeff=CONE)
            BFTransfo%type_Qin(i_Qprim) = -3
            BFTransfo%name_Qin(i_Qprim) = "u" // trim(adjustl(name_beta)) // &
                 trim(adjustl(BFTransfo%name_Frame))
          ELSE
            BFTransfo%Qeuler(2) = set_opel(idf=1, idq=7, alfa=1, indexq=i_Qpoly, coeff=CONE)
            BFTransfo%type_Qin(i_Qprim) = 3
            BFTransfo%name_Qin(i_Qprim) = trim(adjustl(name_beta)) //        &
                 trim(adjustl(BFTransfo%name_Frame))
          END IF
        END IF
        IF (BFTransfo%euler(3)) THEN ! gamma
          i_Qpoly = i_Qpoly + 1
          BFTransfo%Qeuler(3) = set_opel(idf=1, idq=8, alfa=1, indexq=i_Qpoly, coeff=CONE)

          IF (zmat_order_save) THEN
            BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly) = (iQgamma_TO_ivTot-2)*3+3 ! like phi
          ELSE
            BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly) = i_Qpoly
          END IF
          BFTransfo%list_Qprim_TO_Qpoly(BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)) = i_Qpoly
          i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)
          BFTransfo%type_Qin(i_Qprim) = 4
          BFTransfo%name_Qin(i_Qprim) = trim(adjustl(name_gamma)) //         &
                 trim(adjustl(BFTransfo%name_Frame))
        END IF

        IF (debug .OR. print_level > 1) THEN
          write(out_unitp,*) '    END new Frame'
          write(out_unitp,*) '.................................................'
        END IF
        rec_level = rec_level - 1
        IF (rec_level == 0) THEN
          write(out_unitp,*) 'list_Qpoly_TO_Qprim',BFTransfo%list_Qpoly_TO_Qprim(:)
          write(out_unitp,*) 'list_Qprim_TO_Qpoly',BFTransfo%list_Qprim_TO_Qpoly(:)
          !STOP
        END IF
      ELSE   ! new vector in the frame

        name_v = String_TO_String( int_TO_char(BFTransfo%num_vect_in_Frame) // &
                                  "_" // trim(adjustl(BFTransfo%name_Frame)) )


        IF (debug .OR. print_level > 1) THEN
          write(out_unitp,*) ' NEW vector: '
          write(out_unitp,*) '   num_vect_in_Frame,num_vect_in_BF',BFTransfo%num_vect_in_Frame,BFTransfo%num_vect_in_BF
          write(out_unitp,*) '   name_v: ',trim(adjustl(name_v))
        END IF

        BFTransfo%nb_vect_tot = BFTransfo%nb_vect_tot + 1
        IF ( (cart .OR. Li) .AND. BFTransfo%num_vect_in_Frame < 3) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' The second vector of a frame CANNOT be defined:'
          write(out_unitp,*) '   with Cartesian Coordinates or with Li'
          write(out_unitp,*) ' Check your data!'
          STOP
        END IF

        IF (cart) THEN
          i_Qpoly = i_Qpoly + 1
          BFTransfo%Qvec(1) = set_opel(idf=1, idq=1, alfa=1, indexq=i_Qpoly, coeff=CONE) !x
          IF (zmat_order_save) THEN
            BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly) = (iv_tot-2)*3+1
          ELSE
            BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly) = i_Qpoly
          END IF
          BFTransfo%list_Qprim_TO_Qpoly(BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)) = i_Qpoly
          i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)
          BFTransfo%type_Qin(i_Qprim) = 1 ! cart :x
          BFTransfo%name_Qin(i_Qprim) = trim(adjustl(name_x)) // trim(adjustl(name_v))

          i_Qpoly = i_Qpoly + 1
          BFTransfo%Qvec(2) = set_opel(idf=1, idq=1, alfa=1, indexq=i_Qpoly, coeff=CONE) ! y
          IF (zmat_order_save) THEN
            BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly) = (iv_tot-2)*3+2
          ELSE
            BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly) = i_Qpoly
          END IF
          BFTransfo%list_Qprim_TO_Qpoly(BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)) = i_Qpoly
          i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)
          BFTransfo%type_Qin(i_Qprim) = 1 ! cart :y
          BFTransfo%name_Qin(i_Qprim) = trim(adjustl(name_y)) // trim(adjustl(name_v))

          i_Qpoly = i_Qpoly + 1
          BFTransfo%Qvec(3) = set_opel(idf=1, idq=1, alfa=1, indexq=i_Qpoly, coeff=CONE) ! z
          IF (zmat_order_save) THEN
            BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly) = (iv_tot-2)*3+3
          ELSE
            BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly) = i_Qpoly
          END IF
          BFTransfo%list_Qprim_TO_Qpoly(BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)) = i_Qpoly
          i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)
          BFTransfo%type_Qin(i_Qprim) = 1 ! cart :z
          BFTransfo%name_Qin(i_Qprim) = trim(adjustl(name_z)) // trim(adjustl(name_v))

        ELSE
          i_Qpoly = i_Qpoly + 1
          BFTransfo%Qvec(1) = set_opel(idf=1, idq=2, alfa=1, indexq=i_Qpoly, coeff=CONE) !R

          IF (zmat_order_save) THEN
            IF (iv_tot == 2) THEN
              BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly) = 2
            ELSE
              BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly) = (iv_tot-2)*3+1
            END IF
          ELSE
            BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly) = i_Qpoly
          END IF
          BFTransfo%list_Qprim_TO_Qpoly(BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)) = i_Qpoly
          i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)
          BFTransfo%type_Qin(i_Qprim) = 2 ! distance
          BFTransfo%name_Qin(i_Qprim) = trim(adjustl(name_d)) // trim(adjustl(name_v))

          i_Qpoly = i_Qpoly + 1
          IF (zmat_order_save) THEN
            IF (iv_tot == 2) THEN
              BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly) = 3
            ELSE
              BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly) = (iv_tot-2)*3+2
            END IF
          ELSE
            BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly) = i_Qpoly
          END IF
          BFTransfo%list_Qprim_TO_Qpoly(BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)) = i_Qpoly
          i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)
          IF (cos_th) THEN
            BFTransfo%Qvec(2) = set_opel(idf=1, idq=-3, alfa=1, indexq=i_Qpoly, coeff=CONE) !cos(th)

            BFTransfo%type_Qin(i_Qprim) = -3 ! cos(th)
            BFTransfo%name_Qin(i_Qprim) = trim(adjustl(name_u)) // trim(adjustl(name_v))
          ELSE
            BFTransfo%Qvec(2) = set_opel(idf=1, idq=3, alfa=1, indexq=i_Qpoly, coeff=CONE) ! theta
            BFTransfo%type_Qin(i_Qprim) = 3 ! th
            BFTransfo%name_Qin(i_Qprim) = trim(adjustl(name_th)) // trim(adjustl(name_v))
          END IF

          IF (BFTransfo%num_vect_in_Frame > 2) THEN
            i_Qpoly = i_Qpoly + 1
            BFTransfo%Qvec(3) = set_opel(idf=1, idq=4, alfa=1, indexq=i_Qpoly, coeff=CONE) ! phi

            IF (zmat_order_save) THEN
              BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly) = (iv_tot-2)*3+3
            ELSE
              BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly) = i_Qpoly
            END IF
            BFTransfo%list_Qprim_TO_Qpoly(BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)) = i_Qpoly
            i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)
            BFTransfo%type_Qin(i_Qprim) = 4 ! dih
            BFTransfo%name_Qin(i_Qprim) = trim(adjustl(name_dih)) // trim(adjustl(name_v))
          END IF
        END IF
        CALL get_unit_vector_Ei(BFTransfo%Unit_Vector, BFTransfo%Qvec,  &
                     index_v = BFTransfo%num_vect_in_Frame, cart = cart)

      IF (allocated(name_v)) deallocate(name_v)

      END IF



      IF (rec_level == 0) THEN

        BFTransfo%QEuler(1) = set_opel(idf=1, idq=6, alfa=1, indexq=BFTransfo%nb_var+1, coeff=CONE) ! alpha (true Euler)
        IF (cos_th) THEN
          BFTransfo%QEuler(2) = set_opel(idf=1, idq=-7, alfa=1, indexq=BFTransfo%nb_var+2, coeff=CONE) ! beta or cos(beta) (true Euler)
        ELSE
          BFTransfo%QEuler(2) = set_opel(idf=1, idq=7, alfa=1, indexq=BFTransfo%nb_var+2, coeff=CONE) ! beta or cos(beta) (true Euler)
        END IF
        BFTransfo%QEuler(3) = set_opel(idf=1, idq=8, alfa=1, indexq=BFTransfo%nb_var+3, coeff=CONE) ! gamma (true Euler)

        !CALL write_op(BFTransfo%QEuler(1))
        !CALL write_op(BFTransfo%QEuler(2))
        !CALL write_op(BFTransfo%QEuler(3))
      END IF

      IF (debug) CALL RecWrite_BFTransfo(BFTransfo,recur=.FALSE.)
      IF (debug) write(out_unitp,*) 'END ',name_sub

      END SUBROUTINE RecRead_BFTransfo

      !!@description: TODO
      !!@param: TODO
      RECURSIVE SUBROUTINE RecWrite_BFTransfo(BFTransfo,recur)

      TYPE (Type_BFTransfo),intent(in)  :: BFTransfo
      logical,intent(in), optional      :: recur

      integer :: iv,iq
      logical :: recur_loc
      character (len=*), parameter :: name_sub='RecWrite_BFTransfo'

      !IF (.NOT. BFTransfo%Frame) RETURN
      recur_loc = .TRUE.
      IF (present(recur)) recur_loc = recur

      write(out_unitp,*) 'BEGINNING ',name_sub

      CALL Write_FrameType(BFTransfo)

      write(out_unitp,*) 'num_vect_in_Frame',BFTransfo%num_vect_in_Frame
      write(out_unitp,*) 'num_vect_in_BF',BFTransfo%num_vect_in_BF

      write(out_unitp,*) 'num_Frame_in_BF',BFTransfo%num_Frame_in_BF
      write(out_unitp,*) 'num_Frame_in_Container',BFTransfo%num_Frame_in_Container
      write(out_unitp,*) 'name_Frame: ',BFTransfo%name_Frame
      IF (associated(BFTransfo%tab_num_Frame))                          &
           write(out_unitp,*) 'tab_num_Frame',BFTransfo%tab_num_Frame(:)

      write(out_unitp,*) 'nb_vect,nb_vect_tot',                         &
                         BFTransfo%nb_vect,BFTransfo%nb_vect_tot

      write(out_unitp,*) 'nb_var,nb_var_Rot',                           &
                         BFTransfo%nb_var,BFTransfo%nb_var_Rot

      write(out_unitp,*) 'BF,local Frame,euler',BFTransfo%BF,BFTransfo%Frame,BFTransfo%euler(:)

      write(out_unitp,*) 'Cart',BFTransfo%cart
      write(out_unitp,*) 'Li',BFTransfo%Li

      write(out_unitp,*) 'Def_cos_th,cos_th',BFTransfo%Def_cos_th,BFTransfo%cos_th
      write(out_unitp,*) ' Elementary operators for Tana:'
      write(out_unitp,*) ' R, theta ( or u_theta), phi or x, y, z:'
      CALL write_op(BFTransfo%Qvec(1))
      CALL write_op(BFTransfo%Qvec(2))
      CALL write_op(BFTransfo%Qvec(3))
      write(out_unitp,*) ' alpha, beta ( or u_beta), gamma:'
      CALL write_op(BFTransfo%QEuler(1))
      CALL write_op(BFTransfo%QEuler(2))
      CALL write_op(BFTransfo%QEuler(3))
      write(out_unitp,*) ' Unit_vector:'
      CALL write_op(BFTransfo%Unit_Vector)

      IF (associated(BFTransfo%type_Qin) .AND.                          &
          associated(BFTransfo%name_Qin)  )      THEN
        write(out_unitp,*) 'type_Qin',BFTransfo%type_Qin(:)
        write(out_unitp,*) 'name_Qin: ',                                &
              (trim(BFTransfo%name_Qin(iq))," ",iq=1,BFTransfo%nb_var)
      END IF

      IF (associated(BFTransfo%list_Qpoly_TO_Qprim)) THEN
        write(out_unitp,*) 'list_Qpoly_TO_Qprim',BFTransfo%list_Qpoly_TO_Qprim(:)
      END IF

      IF (associated(BFTransfo%list_Qprim_TO_Qpoly)) THEN
        write(out_unitp,*) 'list_Qprim_TO_Qpoly',BFTransfo%list_Qprim_TO_Qpoly(:)
      END IF

      IF (recur_loc .AND. associated(BFTransfo%tab_BFTransfo)) THEN
        DO iv=1,ubound(BFTransfo%tab_BFTransfo,dim=1)
          CALL RecWrite_BFTransfo(BFTransfo%tab_BFTransfo(iv))
        END DO
      END IF

      write(out_unitp,*) 'END ',name_sub
      CALL flush_perso(out_unitp)

      END SUBROUTINE RecWrite_BFTransfo

      RECURSIVE SUBROUTINE calc_PolyTransfo(dnQin,i_Qpoly,dnQout,tab_dnXVect,iv_in,   &
                                            BFTransfo,nderiv)

      TYPE (Type_BFTransfo), intent(in)     :: BFTransfo
      TYPE (Type_dnVec),     intent(inout)  :: dnQin,dnQout
      TYPE (Type_dnVec),     intent(inout) :: tab_dnXVect(:)

      integer, intent (inout) :: i_Qpoly ! index for dnQin
      integer, intent (in) :: iv_in

      integer, intent(in) :: nderiv


      TYPE (Type_dnS)   :: dnd,dnQval,dnCval,dnSval,dnQdih,dnCdih,dnSdih ! one coordinate
      TYPE (Type_dnS)   :: dnf1,dnf2,dnf3,dnfx,dnfy,dnfz
      TYPE (Type_dnS)   :: dna,dnCa,dnSa

      integer :: i_q,iv,liv,uiv,ieuler,i_Qprim,dnErr
      logical :: check
      integer :: nb_var_deriv

!      -----------------------------------------------------------------
      integer :: nderiv_debug = 0
      logical, parameter :: debug = .FALSE.
!      logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='calc_PolyTransfo'
!      -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING RECURSIVE ',name_sub
        write(out_unitp,*) 'nderiv',nderiv
        write(out_unitp,*) 'i_Qpoly,iv_in',i_Qpoly,iv_in
        CALL RecWrite_BFTransfo(BFTransfo,.FALSE.)
      END IF
!      -----------------------------------------------------------------
      nb_var_deriv = dnQin%nb_var_deriv

      ! initialization : allocation....
      CALL alloc_dnSVM(dnd,nb_var_deriv,nderiv)
      CALL alloc_dnSVM(dnQval,nb_var_deriv,nderiv)
      CALL alloc_dnSVM(dnCval,nb_var_deriv,nderiv)
      CALL alloc_dnSVM(dnSval,nb_var_deriv,nderiv)
      CALL alloc_dnSVM(dnQdih,nb_var_deriv,nderiv)
      CALL alloc_dnSVM(dnCdih,nb_var_deriv,nderiv)
      CALL alloc_dnSVM(dnSdih,nb_var_deriv,nderiv)

      CALL alloc_dnSVM(dnf1,nb_var_deriv,nderiv)
      CALL alloc_dnSVM(dnf2,nb_var_deriv,nderiv)
      CALL alloc_dnSVM(dnf3,nb_var_deriv,nderiv)


      IF (BFTransfo%Frame) THEN
        !write(6,*) 'BFTransfo%euler : ',BFTransfo%euler(:)
        !==========================================
        ! 1st vector
        !write(out_unitp,*) 'd0,iv_in,name_frame',iv_in,BFTransfo%name_Frame
        iv = 0
        i_Qpoly = i_Qpoly + 1
        IF (.NOT. associated(BFTransfo%list_Qpoly_TO_Qprim)) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  "list_Qpoly_TO_Qprim" is not associated'
          CALL RecWrite_BFTransfo(BFTransfo,.FALSE.)
          write(out_unitp,*) ' Check the fortran source !!'
          STOP
        END IF
        i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)
        CALL Set_ZERO_TO_dnSVM(tab_dnXVect(iv+1))
        CALL sub_dnVec_TO_dnS(dnQin,dnd,i_Qprim)
        CALL sub_dnS_TO_dnVec(dnd,tab_dnXVect(iv+1),3,nderiv)

        !---------------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*)
          write(out_unitp,*) '-------------------------------------------------'
          write(out_unitp,*) 'Vector :',iv+1,' in ',BFTransfo%name_Frame
          CALL write_dnx(1,3,tab_dnXVect(iv+1),nderiv_debug)
        END IF
        !---------------------------------------------------------------

        liv = 2
        DO iv=1,BFTransfo%nb_vect
          uiv = ubound(tab_dnXVect(:),dim=1)
          !write(out_unitp,*) 'iv,BFTransfo%nb_vect,liv,uiv',iv,BFTransfo%tab_BFTransfo(iv)%nb_vect_tot,liv,uiv
          CALL calc_PolyTransfo(dnQin,i_Qpoly,dnQout,tab_dnXVect(liv:uiv),&
                                iv,BFTransfo%tab_BFTransfo(iv),nderiv)
          liv = liv + BFTransfo%tab_BFTransfo(iv)%nb_vect_tot

        END DO

        !Rotation to change the orientation of the frame
        IF (BFTransfo%Frame_type /= 0) CALL calc_Rot_Vect(tab_dnXVect,BFTransfo,nderiv)

        ! For the overall rotation in the local BF frame (not F0)
        !write(out_unitp,*) 'Overall rotation: ',BFTransfo%name_Frame
        CALL alloc_dnSVM(dna,nb_var_deriv,nderiv)
        CALL alloc_dnSVM(dnCa,nb_var_deriv,nderiv)
        CALL alloc_dnSVM(dnSa,nb_var_deriv,nderiv)
        CALL alloc_dnSVM(dnfx,nb_var_deriv,nderiv)
        CALL alloc_dnSVM(dnfy,nb_var_deriv,nderiv)
        CALL alloc_dnSVM(dnfz,nb_var_deriv,nderiv)
        ! rotation of gamma with respect to Z (not for the first vector)
        IF (BFTransfo%euler(3)) THEN
          ieuler = count(BFTransfo%euler(:))
          i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly+ieuler)

          !write(out_unitp,*) 'Rot gamma / Z , iQgamma',i_Qpoly+ieuler
          CALL sub_dnVec_TO_dnS(dnQin,dna,i_Qprim) ! gamma
          CALL sub_dnS1_TO_dntR2(dna,dnCa,2,nderiv)
          CALL sub_dnS1_TO_dntR2(dna,dnSa,3,nderiv)

          DO iv=2,BFTransfo%nb_vect_tot
            !write(out_unitp,*) 'Rot gamma / Z ',iv
            CALL sub_dnVec_TO_dnS(tab_dnXVect(iv),dnfx,1) ! x
            CALL sub_dnVec_TO_dnS(tab_dnXVect(iv),dnfy,2) ! y
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnfx,dnCa,dnf1,nderiv) ! Rx = cos*x
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnfy,dnSa,dnf3,nderiv) ! temp=sin*y

            CALL sub_dnS1_MINUS_dnS2_TO_dnS3(dnf1,dnf3,dnf1,nderiv)

   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnfx,dnSa,dnf3,nderiv) ! temp= sin*x
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnfy,dnCa,dnf2,nderiv) ! Ry=cos*y

            CALL sub_dnS1_PLUS_dnS2_TO_dnS3(dnf2,dnf3,dnf2,nderiv)! Ry = Ry+temp


            CALL sub_dnS_TO_dnVec(dnf1,tab_dnXVect(iv),1,nderiv) !x
            CALL sub_dnS_TO_dnVec(dnf2,tab_dnXVect(iv),2,nderiv) !y
          END DO
        END IF

        ! rotation of beta with respect to Y (all vectors)
        IF (BFTransfo%euler(2)) THEN
          ieuler = count(BFTransfo%euler(1:2))
          i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly+ieuler)

          !write(out_unitp,*) 'Rot beta / Y ,iQbeta',i_Qpoly+ieuler
          CALL sub_dnVec_TO_dnS(dnQin,dna,i_Qprim) ! beta
          IF (BFTransfo%type_Qin(i_Qprim) == 3) THEN
            CALL sub_dnS1_TO_dntR2(dna,dnCa,2,nderiv)
            CALL sub_dnS1_TO_dntR2(dna,dnSa,3,nderiv)
          ELSE ! type_Qin(i_Qprim) == -3 (using cos(beta) as coordinate)
            CALL sub_dnS1_TO_dnS2(dna,dnCa,nderiv)
            CALL sub_dnS1_TO_dntR2(dna,dnSa,4,nderiv,dnErr=dnErr)
            IF (dnErr /= 0) THEN
              write(out_unitp,*) ' ERROR in ',name_sub
              write(out_unitp,*) '   ERROR in the sub_dntf call for the coordinates, i_Qprim:',i_Qprim
              STOP 'ERROR in sub_dntf called from calc_PolyTransfo'
            END IF
          END IF

          DO iv=1,BFTransfo%nb_vect_tot
            !write(out_unitp,*) 'Rot beta / Y ',iv
            CALL sub_dnVec_TO_dnS(tab_dnXVect(iv),dnfx,1) ! x
            CALL sub_dnVec_TO_dnS(tab_dnXVect(iv),dnfz,3) ! z
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnfx,dnCa,dnf1,nderiv) ! Rx = cos*x
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnfz,dnSa,dnf2,nderiv) ! temp=sin*z

            CALL sub_dnS1_PLUS_dnS2_TO_dnS3(dnf1,dnf2,dnf1,nderiv)!  Rx = Rx+temp


   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnfx,dnSa,dnf2,nderiv) ! temp = sin*x
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnfz,dnCa,dnf3,nderiv) ! Rz=cos*z

            CALL sub_dnS1_MINUS_dnS2_TO_dnS3(dnf3,dnf2,dnf3,nderiv)! Rz = Rz-temp

            CALL sub_dnS_TO_dnVec(dnf1,tab_dnXVect(iv),1,nderiv) !x
            CALL sub_dnS_TO_dnVec(dnf3,tab_dnXVect(iv),3,nderiv) !z

          END DO
        END IF

        ! rotation of alpha with respect to Z (all vectors)
        IF (BFTransfo%euler(1)) THEN
          ieuler = 1
          i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly+ieuler)

          !write(out_unitp,*) 'Rot alpha / Z , iQalpha',i_Qpoly+ieuler
          CALL sub_dnVec_TO_dnS(dnQin,dna,i_Qprim) ! alpha
          CALL sub_dnS1_TO_dntR2(dna,dnCa,2,nderiv)
          CALL sub_dnS1_TO_dntR2(dna,dnSa,3,nderiv)
          DO iv=1,BFTransfo%nb_vect_tot
            !write(out_unitp,*) 'Rot alpha / Z ',iv
            CALL sub_dnVec_TO_dnS(tab_dnXVect(iv),dnfx,1) ! x
            CALL sub_dnVec_TO_dnS(tab_dnXVect(iv),dnfy,2) ! y
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnfx,dnCa,dnf1,nderiv) ! Rx = cos*x
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnfy,dnSa,dnf3,nderiv) ! temp=sin*y

            CALL sub_dnS1_MINUS_dnS2_TO_dnS3(dnf1,dnf3,dnf1,nderiv)! Rx = Rx-temp


   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnfx,dnSa,dnf3,nderiv) ! temp= sin*x
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnfy,dnCa,dnf2,nderiv) ! Ry=cos*y

            CALL sub_dnS1_PLUS_dnS2_TO_dnS3(dnf2,dnf3,dnf2,nderiv)! Ry = Rx+temp

            CALL sub_dnS_TO_dnVec(dnf1,tab_dnXVect(iv),1,nderiv) !x
            CALL sub_dnS_TO_dnVec(dnf2,tab_dnXVect(iv),2,nderiv) !y
          END DO
        END IF

        i_Qpoly = i_Qpoly + count(BFTransfo%euler(:))

        CALL dealloc_dnSVM(dna)
        CALL dealloc_dnSVM(dnCa)
        CALL dealloc_dnSVM(dnSa)
        CALL dealloc_dnSVM(dnfx)
        CALL dealloc_dnSVM(dnfy)
        CALL dealloc_dnSVM(dnfz)
      ELSE
        IF (iv_in == 1) THEN
          !write(out_unitp,*) 'd1,th1,iv_in,name_frame',iv_in,BFTransfo%name_Frame
          i_Qpoly = i_Qpoly + 1
          i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)
          CALL sub_dnVec_TO_dnS(dnQin,dnd,i_Qprim)
          !CALL Write_dnSVM(dnd,nderiv)

          i_Qpoly = i_Qpoly + 1
          i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)
          CALL sub_dnVec_TO_dnS(dnQin,dnQval,i_Qprim)
          !CALL Write_dnSVM(dnQval,nderiv)
          IF (BFTransfo%type_Qin(i_Qprim) == 3) THEN
            CALL sub_dnS1_TO_dntR2(dnQval,dnCval,2,nderiv)
            CALL sub_dnS1_TO_dntR2(dnQval,dnSval,3,nderiv)
          ELSE ! type_Qin(i_Qprim) == -3 (using Q=cos(val) as coordinate)
            CALL sub_dnS1_TO_dnS2(dnQval,dnCval,nderiv)
            CALL sub_dnS1_TO_dntR2(dnQval,dnSval,4,nderiv,dnErr=dnErr)
            IF (dnErr /= 0) THEN
              write(out_unitp,*) ' ERROR in ',name_sub
              write(out_unitp,*) '   ERROR in the sub_dntf call for the coordinates, i_Qprim:',i_Qprim
              STOP 'ERROR in sub_dntf called from calc_PolyTransfo'
            END IF
          END IF
          ! CALL Write_dnSVM(dnCval,nderiv)
          ! CALL Write_dnSVM(dnSval,nderiv)

          CALL Set_ZERO_TO_dnSVM(tab_dnXVect(1))

          !for x coordinates
          ! d0w = d0d * d0sin(Qval)
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnd,dnSval,dnf1,nderiv)
          CALL sub_dnS_TO_dnVec(dnf1,tab_dnXVect(1),1,nderiv)

          !for z coordinates
          ! d0w = d0d * d0cos(Qval)
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnd,dnCval,dnf2,nderiv)
          CALL sub_dnS_TO_dnVec(dnf2,tab_dnXVect(1),3,nderiv)

        ELSE IF (iv_in > 1) THEN
          i_Qpoly = i_Qpoly + 1
          i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)

          !write(out_unitp,*) 'd2,th2,iv_in,name_frame',iv_in,BFTransfo%name_Frame
          IF (BFTransfo%type_Qin(i_Qprim) == 1) THEN ! cartesian coordinates
            CALL sub_dnVec_TO_dnS(dnQin,dnf1,i_Qprim) !x
            i_Qpoly = i_Qpoly + 1
            i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)
            CALL sub_dnVec_TO_dnS(dnQin,dnf2,i_Qprim) !y
            i_Qpoly = i_Qpoly + 1
            i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)
            CALL sub_dnVec_TO_dnS(dnQin,dnf3,i_Qprim) !z
          ELSE ! spherical coordinates
            CALL sub_dnVec_TO_dnS(dnQin,dnd,i_Qprim)
            i_Qpoly = i_Qpoly + 1
            i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)
            CALL sub_dnVec_TO_dnS(dnQin,dnQval,i_Qprim)
            IF (BFTransfo%type_Qin(i_Qprim) == 3) THEN
              CALL sub_dnS1_TO_dntR2(dnQval,dnCval,2,nderiv)
              CALL sub_dnS1_TO_dntR2(dnQval,dnSval,3,nderiv)
            ELSE ! type_Qin(i_Qprim) == -3 (using Q=cos(val) as coordinate)
              CALL sub_dnS1_TO_dnS2(dnQval,dnCval,nderiv)
              CALL sub_dnS1_TO_dntR2(dnQval,dnSval,4,nderiv,dnErr=dnErr)
              IF (dnErr /= 0) THEN
                write(out_unitp,*) ' ERROR in ',name_sub
                write(out_unitp,*) '   ERROR in the sub_dntf call for the coordinates, i_Qprim:',i_Qprim
                STOP 'ERROR in sub_dntf called from calc_PolyTransfo'
              END IF
            END IF

            i_Qpoly = i_Qpoly + 1
            i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)
            CALL sub_dnVec_TO_dnS(dnQin,dnQdih,i_Qprim)
            CALL sub_dnS1_TO_dntR2(dnQdih,dnCdih,2,nderiv)
            CALL sub_dnS1_TO_dntR2(dnQdih,dnSdih,3,nderiv)


            !-----------------------------------------------------------
            !d0f3 = d0d * d0sval (tempory)
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnd,dnSval,dnf3,nderiv)

            !d0f1 = d0f3 * d0cdih  = (d0d * d0sval) * d0cdih
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnf3,dnCdih,dnf1,nderiv)

            !d0f2 = d0f3 * d0sdih  = (d0d * d0sval) * d0sdih
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnf3,dnSdih,dnf2,nderiv)
            !-----------------------------------------------------------

            !-----------------------------------------------------------
            ! d0f3 = d0d * d0cval
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnd,dnCval,dnf3,nderiv)

            !-----------------------------------------------------------
          END IF

          CALL sub_dnS_TO_dnVec(dnf1,tab_dnXVect(1),1,nderiv) !x
          CALL sub_dnS_TO_dnVec(dnf2,tab_dnXVect(1),2,nderiv) !y
          CALL sub_dnS_TO_dnVec(dnf3,tab_dnXVect(1),3,nderiv) !z

        ELSE
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  Wrong iv_in for Frame = F',iv_in,BFTransfo%Frame
          STOP
        END IF

        !-------------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*)
          write(out_unitp,*) '-----------------------------------------------'
          write(out_unitp,*) 'Vector :',iv_in,' in ',BFTransfo%name_Frame
          CALL write_dnx(1,3,tab_dnXVect(1),nderiv_debug)
        END IF
        !-------------------------------------------------------------

      END IF


      CALL dealloc_dnSVM(dnd)
      CALL dealloc_dnSVM(dnQval)
      CALL dealloc_dnSVM(dnCval)
      CALL dealloc_dnSVM(dnSval)
      CALL dealloc_dnSVM(dnQdih)
      CALL dealloc_dnSVM(dnCdih)
      CALL dealloc_dnSVM(dnSdih)
      CALL dealloc_dnSVM(dnf1)
      CALL dealloc_dnSVM(dnf2)
      CALL dealloc_dnSVM(dnf3)

      ! finalization : transfert tab_dnXVect => dnQout
      IF (BFTransfo%Frame .AND. size(BFTransfo%tab_num_Frame) == 1 ) THEN

        CALL Set_ZERO_TO_dnSVM(dnQout) ! initialization: useless !!!

        DO iv=1,BFTransfo%nb_vect_tot
          IF (debug) CALL write_dnx(1,3,tab_dnXVect(iv),nderiv_debug)
          CALL sub3_dnVec_TOxf(dnQout,3*iv-2,tab_dnXVect(iv),nderiv)
        END DO
      END IF

      !----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'dnQout'
        CALL Write_dnSVM(dnQout,nderiv_debug)
        write(out_unitp,*)
        write(out_unitp,*) 'END RECURSIVE ',name_sub
        write(out_unitp,*)
      END IF


      END SUBROUTINE calc_PolyTransfo

      SUBROUTINE calc_PolyTransfo_outTOin(dnQin,dnQout,BFTransfo,nderiv)

      TYPE (Type_BFTransfo), intent(in)    :: BFTransfo
      TYPE (Type_dnVec),     intent(inout) :: dnQin,dnQout

      integer, intent(in) :: nderiv


      TYPE (Type_dnVec) :: UnitVect_F0(3) ! unit vectors (ndim=3)
      TYPE (Type_dnVec) :: dnVect(BFTransfo%nb_vect_tot) ! table of vectors


      integer :: iv,iv_tot,iQin

!      -----------------------------------------------------------------
      integer :: nderiv_debug = 0
      integer :: err_mem,memory
      logical, parameter :: debug = .FALSE.
!      logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='calc_PolyTransfo_outTOin'
!      -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nderiv',nderiv
        write(out_unitp,*) 'dnQvect (dnQout): '
        CALL Write_dnSVM(dnQout,nderiv_debug)
        CALL RecWrite_BFTransfo(BFTransfo,.FALSE.)
      END IF
!      -----------------------------------------------------------------

      IF (nderiv /= 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) '  This subroutine cannot be use with nderiv > 0',nderiv
        write(out_unitp,*) '  Check the fortran source !!'
        STOP
      END IF


      CALL alloc_dnSVM(UnitVect_F0(1),nb_var_vec=3,nderiv=0)
      UnitVect_F0(1)%d0(1) = ONE

      CALL alloc_dnSVM(UnitVect_F0(2),nb_var_vec=3,nderiv=0)
      UnitVect_F0(2)%d0(3) = ONE

      CALL alloc_dnSVM(UnitVect_F0(3),nb_var_vec=3,nderiv=0)
      UnitVect_F0(3)%d0(3) = ONE

      DO iv=1,BFTransfo%nb_vect_tot
        CALL alloc_dnSVM(dnVect(iv),nb_var_vec=3,nderiv=0)
        dnVect(iv)%d0(1:3) = dnQout%d0(3*iv-2:3*iv)
      END DO

      iv_tot = 0
      iQin   = 0
      CALL RecGet_Vec_Fi_For_poly(dnVect,BFTransfo%nb_vect_tot,         &
                                  iv_tot,1,UnitVect_F0,                 &
                                  dnQin,iQin,BFTransfo)
      !write(out_unitp,*) 'iv_tot',iv_tot

      CALL dealloc_dnSVM(UnitVect_F0(1))
      CALL dealloc_dnSVM(UnitVect_F0(2))
      CALL dealloc_dnSVM(UnitVect_F0(3))


        !----------------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*)
          write(out_unitp,*) 'dnQin (dnQpoly):'
          CALL Write_dnSVM(dnQin,nderiv_debug)
          write(out_unitp,*)
          write(out_unitp,*) 'END ',name_sub
          write(out_unitp,*)
        END IF


      END SUBROUTINE calc_PolyTransfo_outTOin

      RECURSIVE SUBROUTINE Export_Fortran_PolyTransfo(dnQin,i_Qpoly,dnQout,&
                                               tab_dnXVect,iv_in,BFTransfo)

      TYPE (Type_BFTransfo), intent(in)     :: BFTransfo
      TYPE (Type_dnVec),     intent(inout)  :: dnQin,dnQout
      TYPE (Type_dnVec),     intent(inout) :: tab_dnXVect(:)

      integer, intent (inout) :: i_Qpoly ! index for dnQin
      integer, intent (in) :: iv_in


      TYPE (Type_dnS)   :: dnd,dnQval,dnCval,dnSval,dnQdih,dnCdih,dnSdih ! one coordinate
      TYPE (Type_dnS)   :: dnf1,dnf2,dnf3,dnfx,dnfy,dnfz
      TYPE (Type_dnS)   :: dna,dnCa,dnSa

      integer :: i_q,iv,liv,uiv,ieuler,i_Qprim,dnErr
      logical :: check
      integer :: nb_var_deriv,nderiv=0

!      -----------------------------------------------------------------
      integer :: nderiv_debug = 0
      !logical, parameter :: debug = .FALSE.
      logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='Export_Fortran_PolyTransfo'
!      -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING RECURSIVE ',name_sub
        write(out_unitp,*) 'nderiv',nderiv
        write(out_unitp,*) 'i_Qpoly,iv_in',i_Qpoly,iv_in
        CALL RecWrite_BFTransfo(BFTransfo,.FALSE.)
      END IF
!      -----------------------------------------------------------------
      nb_var_deriv = dnQin%nb_var_deriv

      ! initialization : allocation....
      CALL alloc_dnSVM(dnd,nb_var_deriv,nderiv)
      CALL alloc_dnSVM(dnQval,nb_var_deriv,nderiv)
      CALL alloc_dnSVM(dnCval,nb_var_deriv,nderiv)
      CALL alloc_dnSVM(dnSval,nb_var_deriv,nderiv)
      CALL alloc_dnSVM(dnQdih,nb_var_deriv,nderiv)
      CALL alloc_dnSVM(dnCdih,nb_var_deriv,nderiv)
      CALL alloc_dnSVM(dnSdih,nb_var_deriv,nderiv)

      CALL alloc_dnSVM(dnf1,nb_var_deriv,nderiv)
      CALL alloc_dnSVM(dnf2,nb_var_deriv,nderiv)
      CALL alloc_dnSVM(dnf3,nb_var_deriv,nderiv)


      IF (BFTransfo%Frame) THEN
        write(6,*) 'BFTransfo%euler : ',BFTransfo%euler(:)
        !==========================================
        ! 1st vector
        write(out_unitp,*) 'd0,iv_in,name_frame',iv_in,BFTransfo%name_Frame
        iv = 0
        i_Qpoly = i_Qpoly + 1
        IF (.NOT. associated(BFTransfo%list_Qpoly_TO_Qprim)) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  "list_Qpoly_TO_Qprim" is not associated'
          CALL RecWrite_BFTransfo(BFTransfo,.FALSE.)
          write(out_unitp,*) ' Check the fortran source !!'
          STOP
        END IF
        i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)
        CALL Set_ZERO_TO_dnSVM(tab_dnXVect(iv+1))
        CALL sub_dnVec_TO_dnS(dnQin,dnd,i_Qprim)
        CALL sub_dnS_TO_dnVec(dnd,tab_dnXVect(iv+1),3,nderiv)

        !---------------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*)
          write(out_unitp,*) '-------------------------------------------------'
          write(out_unitp,*) 'Vector :',iv+1,' in ',BFTransfo%name_Frame
          CALL write_dnx(1,3,tab_dnXVect(iv+1),nderiv_debug)
        END IF
        !---------------------------------------------------------------

        liv = 2
        DO iv=1,BFTransfo%nb_vect
          uiv = ubound(tab_dnXVect(:),dim=1)
          !write(out_unitp,*) 'iv,BFTransfo%nb_vect,liv,uiv',iv,BFTransfo%tab_BFTransfo(iv)%nb_vect_tot,liv,uiv
          CALL calc_PolyTransfo(dnQin,i_Qpoly,dnQout,tab_dnXVect(liv:uiv),&
                                iv,BFTransfo%tab_BFTransfo(iv),nderiv)
          liv = liv + BFTransfo%tab_BFTransfo(iv)%nb_vect_tot

        END DO

        !Rotation to change the orientation of the frame
        IF (BFTransfo%Frame_type /= 0) CALL calc_Rot_Vect(tab_dnXVect,BFTransfo,nderiv)

        ! For the overall rotation in the local BF frame (not F0)
        write(out_unitp,*) 'Overall rotation: ',BFTransfo%name_Frame
        CALL alloc_dnSVM(dna,nb_var_deriv,nderiv)
        CALL alloc_dnSVM(dnCa,nb_var_deriv,nderiv)
        CALL alloc_dnSVM(dnSa,nb_var_deriv,nderiv)
        CALL alloc_dnSVM(dnfx,nb_var_deriv,nderiv)
        CALL alloc_dnSVM(dnfy,nb_var_deriv,nderiv)
        CALL alloc_dnSVM(dnfz,nb_var_deriv,nderiv)
        ! rotation of gamma with respect to Z (not for the first vector)
        IF (BFTransfo%euler(3)) THEN
          ieuler = count(BFTransfo%euler(:))
          i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly+ieuler)

          !write(out_unitp,*) 'Rot gamma / Z , iQgamma',i_Qpoly+ieuler
          CALL sub_dnVec_TO_dnS(dnQin,dna,i_Qprim) ! gamma
          CALL sub_dnS1_TO_dntR2(dna,dnCa,2,nderiv)
          CALL sub_dnS1_TO_dntR2(dna,dnSa,3,nderiv)

          DO iv=2,BFTransfo%nb_vect_tot
            !write(out_unitp,*) 'Rot gamma / Z ',iv
            CALL sub_dnVec_TO_dnS(tab_dnXVect(iv),dnfx,1) ! x
            CALL sub_dnVec_TO_dnS(tab_dnXVect(iv),dnfy,2) ! y
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnfx,dnCa,dnf1,nderiv) ! Rx = cos*x
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnfy,dnSa,dnf3,nderiv) ! temp=sin*y

            CALL sub_dnS1_MINUS_dnS2_TO_dnS3(dnf1,dnf3,dnf1,nderiv)

   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnfx,dnSa,dnf3,nderiv) ! temp= sin*x
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnfy,dnCa,dnf2,nderiv) ! Ry=cos*y

            CALL sub_dnS1_PLUS_dnS2_TO_dnS3(dnf2,dnf3,dnf2,nderiv)! Ry = Ry+temp


            CALL sub_dnS_TO_dnVec(dnf1,tab_dnXVect(iv),1,nderiv) !x
            CALL sub_dnS_TO_dnVec(dnf2,tab_dnXVect(iv),2,nderiv) !y
          END DO
        END IF

        ! rotation of beta with respect to Y (all vectors)
        IF (BFTransfo%euler(2)) THEN
          ieuler = count(BFTransfo%euler(1:2))
          i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly+ieuler)

          !write(out_unitp,*) 'Rot beta / Y ,iQbeta',i_Qpoly+ieuler
          CALL sub_dnVec_TO_dnS(dnQin,dna,i_Qprim) ! beta
          IF (BFTransfo%type_Qin(i_Qprim) == 3) THEN
            CALL sub_dnS1_TO_dntR2(dna,dnCa,2,nderiv)
            CALL sub_dnS1_TO_dntR2(dna,dnSa,3,nderiv)
          ELSE ! type_Qin(i_Qprim) == -3 (using cos(beta) as coordinate)
            CALL sub_dnS1_TO_dnS2(dna,dnCa,nderiv)
            CALL sub_dnS1_TO_dntR2(dna,dnSa,4,nderiv,dnErr=dnErr)
            IF (dnErr /= 0) THEN
              write(out_unitp,*) ' ERROR in ',name_sub
              write(out_unitp,*) '   ERROR in the sub_dntf call for the coordinates, i_Qprim:',i_Qprim
              STOP 'ERROR in sub_dntf called from calc_PolyTransfo'
            END IF
          END IF

          DO iv=1,BFTransfo%nb_vect_tot
            !write(out_unitp,*) 'Rot beta / Y ',iv
            CALL sub_dnVec_TO_dnS(tab_dnXVect(iv),dnfx,1) ! x
            CALL sub_dnVec_TO_dnS(tab_dnXVect(iv),dnfz,3) ! z
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnfx,dnCa,dnf1,nderiv) ! Rx = cos*x
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnfz,dnSa,dnf2,nderiv) ! temp=sin*z

            CALL sub_dnS1_PLUS_dnS2_TO_dnS3(dnf1,dnf2,dnf1,nderiv)!  Rx = Rx+temp


   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnfx,dnSa,dnf2,nderiv) ! temp = sin*x
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnfz,dnCa,dnf3,nderiv) ! Rz=cos*z

            CALL sub_dnS1_MINUS_dnS2_TO_dnS3(dnf3,dnf2,dnf3,nderiv)! Rz = Rz-temp

            CALL sub_dnS_TO_dnVec(dnf1,tab_dnXVect(iv),1,nderiv) !x
            CALL sub_dnS_TO_dnVec(dnf3,tab_dnXVect(iv),3,nderiv) !z

          END DO
        END IF

        ! rotation of alpha with respect to Z (all vectors)
        IF (BFTransfo%euler(1)) THEN
          ieuler = 1
          i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly+ieuler)

          !write(out_unitp,*) 'Rot alpha / Z , iQalpha',i_Qpoly+ieuler
          CALL sub_dnVec_TO_dnS(dnQin,dna,i_Qprim) ! alpha
          CALL sub_dnS1_TO_dntR2(dna,dnCa,2,nderiv)
          CALL sub_dnS1_TO_dntR2(dna,dnSa,3,nderiv)
          DO iv=1,BFTransfo%nb_vect_tot
            !write(out_unitp,*) 'Rot alpha / Z ',iv
            CALL sub_dnVec_TO_dnS(tab_dnXVect(iv),dnfx,1) ! x
            CALL sub_dnVec_TO_dnS(tab_dnXVect(iv),dnfy,2) ! y
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnfx,dnCa,dnf1,nderiv) ! Rx = cos*x
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnfy,dnSa,dnf3,nderiv) ! temp=sin*y

            CALL sub_dnS1_MINUS_dnS2_TO_dnS3(dnf1,dnf3,dnf1,nderiv)! Rx = Rx-temp


   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnfx,dnSa,dnf3,nderiv) ! temp= sin*x
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnfy,dnCa,dnf2,nderiv) ! Ry=cos*y

            CALL sub_dnS1_PLUS_dnS2_TO_dnS3(dnf2,dnf3,dnf2,nderiv)! Ry = Rx+temp

            CALL sub_dnS_TO_dnVec(dnf1,tab_dnXVect(iv),1,nderiv) !x
            CALL sub_dnS_TO_dnVec(dnf2,tab_dnXVect(iv),2,nderiv) !y
          END DO
        END IF

        i_Qpoly = i_Qpoly + count(BFTransfo%euler(:))

        CALL dealloc_dnSVM(dna)
        CALL dealloc_dnSVM(dnCa)
        CALL dealloc_dnSVM(dnSa)
        CALL dealloc_dnSVM(dnfx)
        CALL dealloc_dnSVM(dnfy)
        CALL dealloc_dnSVM(dnfz)
      ELSE
        IF (iv_in == 1) THEN
          !write(out_unitp,*) 'd1,th1,iv_in,name_frame',iv_in,BFTransfo%name_Frame
          i_Qpoly = i_Qpoly + 1
          i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)
          CALL sub_dnVec_TO_dnS(dnQin,dnd,i_Qprim)
          !CALL Write_dnSVM(dnd,nderiv)

          i_Qpoly = i_Qpoly + 1
          i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)
          CALL sub_dnVec_TO_dnS(dnQin,dnQval,i_Qprim)
          !CALL Write_dnSVM(dnQval,nderiv)
          IF (BFTransfo%type_Qin(i_Qprim) == 3) THEN
            CALL sub_dnS1_TO_dntR2(dnQval,dnCval,2,nderiv)
            CALL sub_dnS1_TO_dntR2(dnQval,dnSval,3,nderiv)
          ELSE ! type_Qin(i_Qprim) == -3 (using Q=cos(val) as coordinate)
            CALL sub_dnS1_TO_dnS2(dnQval,dnCval,nderiv)
            CALL sub_dnS1_TO_dntR2(dnQval,dnSval,4,nderiv,dnErr=dnErr)
            IF (dnErr /= 0) THEN
              write(out_unitp,*) ' ERROR in ',name_sub
              write(out_unitp,*) '   ERROR in the sub_dntf call for the coordinates, i_Qprim:',i_Qprim
              STOP 'ERROR in sub_dntf called from calc_PolyTransfo'
            END IF
          END IF
          ! CALL Write_dnSVM(dnCval,nderiv)
          ! CALL Write_dnSVM(dnSval,nderiv)

          CALL Set_ZERO_TO_dnSVM(tab_dnXVect(1))

          !for x coordinates
          ! d0w = d0d * d0sin(Qval)
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnd,dnSval,dnf1,nderiv)
          CALL sub_dnS_TO_dnVec(dnf1,tab_dnXVect(1),1,nderiv)

          !for z coordinates
          ! d0w = d0d * d0cos(Qval)
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnd,dnCval,dnf2,nderiv)
          CALL sub_dnS_TO_dnVec(dnf2,tab_dnXVect(1),3,nderiv)

        ELSE IF (iv_in > 1) THEN
          i_Qpoly = i_Qpoly + 1
          i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)

          !write(out_unitp,*) 'd2,th2,iv_in,name_frame',iv_in,BFTransfo%name_Frame
          IF (BFTransfo%type_Qin(i_Qprim) == 1) THEN ! cartesian coordinates
            CALL sub_dnVec_TO_dnS(dnQin,dnf1,i_Qprim) !x
            i_Qpoly = i_Qpoly + 1
            i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)
            CALL sub_dnVec_TO_dnS(dnQin,dnf2,i_Qprim) !y
            i_Qpoly = i_Qpoly + 1
            i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)
            CALL sub_dnVec_TO_dnS(dnQin,dnf3,i_Qprim) !z
          ELSE ! spherical coordinates
            CALL sub_dnVec_TO_dnS(dnQin,dnd,i_Qprim)
            i_Qpoly = i_Qpoly + 1
            i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)
            CALL sub_dnVec_TO_dnS(dnQin,dnQval,i_Qprim)
            IF (BFTransfo%type_Qin(i_Qprim) == 3) THEN
              CALL sub_dnS1_TO_dntR2(dnQval,dnCval,2,nderiv)
              CALL sub_dnS1_TO_dntR2(dnQval,dnSval,3,nderiv)
            ELSE ! type_Qin(i_Qprim) == -3 (using Q=cos(val) as coordinate)
              CALL sub_dnS1_TO_dnS2(dnQval,dnCval,nderiv)
              CALL sub_dnS1_TO_dntR2(dnQval,dnSval,4,nderiv,dnErr=dnErr)
              IF (dnErr /= 0) THEN
                write(out_unitp,*) ' ERROR in ',name_sub
                write(out_unitp,*) '   ERROR in the sub_dntf call for the coordinates, i_Qprim:',i_Qprim
                STOP 'ERROR in sub_dntf called from calc_PolyTransfo'
              END IF
            END IF

            i_Qpoly = i_Qpoly + 1
            i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(i_Qpoly)
            CALL sub_dnVec_TO_dnS(dnQin,dnQdih,i_Qprim)
            CALL sub_dnS1_TO_dntR2(dnQdih,dnCdih,2,nderiv)
            CALL sub_dnS1_TO_dntR2(dnQdih,dnSdih,3,nderiv)


            !-----------------------------------------------------------
            !d0f3 = d0d * d0sval (tempory)
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnd,dnSval,dnf3,nderiv)

            !d0f1 = d0f3 * d0cdih  = (d0d * d0sval) * d0cdih
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnf3,dnCdih,dnf1,nderiv)

            !d0f2 = d0f3 * d0sdih  = (d0d * d0sval) * d0sdih
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnf3,dnSdih,dnf2,nderiv)
            !-----------------------------------------------------------

            !-----------------------------------------------------------
            ! d0f3 = d0d * d0cval
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnd,dnCval,dnf3,nderiv)

            !-----------------------------------------------------------
          END IF

          CALL sub_dnS_TO_dnVec(dnf1,tab_dnXVect(1),1,nderiv) !x
          CALL sub_dnS_TO_dnVec(dnf2,tab_dnXVect(1),2,nderiv) !y
          CALL sub_dnS_TO_dnVec(dnf3,tab_dnXVect(1),3,nderiv) !z

        ELSE
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  Wrong iv_in for Frame = F',iv_in,BFTransfo%Frame
          STOP
        END IF

        !-------------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*)
          write(out_unitp,*) '-----------------------------------------------'
          write(out_unitp,*) 'Vector :',iv_in,' in ',BFTransfo%name_Frame
          CALL write_dnx(1,3,tab_dnXVect(1),nderiv_debug)
        END IF
        !-------------------------------------------------------------

      END IF


      CALL dealloc_dnSVM(dnd)
      CALL dealloc_dnSVM(dnQval)
      CALL dealloc_dnSVM(dnCval)
      CALL dealloc_dnSVM(dnSval)
      CALL dealloc_dnSVM(dnQdih)
      CALL dealloc_dnSVM(dnCdih)
      CALL dealloc_dnSVM(dnSdih)
      CALL dealloc_dnSVM(dnf1)
      CALL dealloc_dnSVM(dnf2)
      CALL dealloc_dnSVM(dnf3)

      ! finalization : transfert tab_dnXVect => dnQout
      IF (BFTransfo%Frame .AND. size(BFTransfo%tab_num_Frame) == 1 ) THEN

        CALL Set_ZERO_TO_dnSVM(dnQout) ! initialization: useless !!!

        DO iv=1,BFTransfo%nb_vect_tot
          IF (debug) CALL write_dnx(1,3,tab_dnXVect(iv),nderiv_debug)
          CALL sub3_dnVec_TOxf(dnQout,3*iv-2,tab_dnXVect(iv),nderiv)
        END DO
      END IF

      !----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'dnQout'
        CALL Write_dnSVM(dnQout,nderiv_debug)
        write(out_unitp,*)
        write(out_unitp,*) 'END RECURSIVE ',name_sub
        write(out_unitp,*)
      END IF


      END SUBROUTINE Export_Fortran_PolyTransfo

      RECURSIVE SUBROUTINE RecGet_Vec_Fi_For_poly(tab_Vect_Fi,          &
                                              nb_vect_tot,iv_tot,iv_Fi, &
                                                  UnitVect_Fi,          &
                                               dnQpoly,iQpoly,BFTransfo)

      integer, intent(in)    :: nb_vect_tot,iv_Fi
      integer, intent(inout) :: iv_tot,iQpoly
      TYPE (Type_dnVec), intent(in) :: tab_Vect_Fi(nb_vect_tot) ! table of vectors (ndim=3)
      TYPE (Type_dnVec),     intent(inout) :: dnQpoly
      TYPE (Type_dnVec), intent(in) :: UnitVect_Fi(3) ! unit vectors (ndim=3)

      TYPE (Type_dnVec) :: UnitVect_Fij(3) ! unit vectors if frame=t (ndim=3)
      TYPE (Type_BFTransfo), intent(in) :: BFTransfo


      integer :: iv_Fij,i_Qprim

      integer :: nb_vect
      logical :: Frame,cart,cos_th
      character (len=Name_len) :: name_d,name_th,name_dih,              &
                                  name_x,name_y,name_z,                 &
                                  name_alpha,name_beta,name_gamma

      real (kind=Rkind) :: Riv,px,py,pz
      real (kind=Rkind) :: alphaiv,betaiv,ubetaiv,gammaiv,sgamma,cgamma
      real (kind=Rkind) :: uiv,thiv,phiv

      real (kind=Rkind), parameter :: radTOdeg = 180._Rkind / pi

!      -----------------------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug = .FALSE.
!      logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='RecGet_Vec_Fi_For_poly'
!      -----------------------------------------------------------------

      iv_tot = iv_tot + 1
      !write(out_unitp,*) 'RecGet_Vec_Fi: nb_vect_tot,iv_Fi,iv_tot',nb_vect_tot,iv_Fi,iv_tot
      !write(out_unitp,*) 'vect:',tab_Vect_Fi(iv_tot)%d0

      nb_vect = BFTransfo%nb_vect
      Frame   = BFTransfo%Frame

      IF (Frame) THEN
        IF (debug) write(out_unitp,*) '============================================='
        IF (debug) write(out_unitp,*) ' Frame = T, iv_tot',iv_tot
        ! norm of the vector (distance)
        Riv = sqrt(dot_product(tab_Vect_Fi(iv_tot)%d0,tab_Vect_Fi(iv_tot)%d0))

        ! The 3 unit vectors in the new frame Fij

        CALL alloc_dnSVM(UnitVect_Fij(1),nb_var_vec=3,nderiv=0)
        CALL alloc_dnSVM(UnitVect_Fij(2),nb_var_vec=3,nderiv=0)
        CALL alloc_dnSVM(UnitVect_Fij(3),nb_var_vec=3,nderiv=0)

        UnitVect_Fij(3)%d0 = tab_Vect_Fi(iv_tot)%d0 / Riv ! ez_Fij
        IF (nb_vect > 0) THEN
          pz = dot_product(tab_Vect_Fi(iv_tot+1)%d0,UnitVect_Fij(3)%d0)
          UnitVect_Fij(1)%d0 = tab_Vect_Fi(iv_tot+1)%d0 - pz * UnitVect_Fij(3)%d0
          UnitVect_Fij(1)%d0 = UnitVect_Fij(1)%d0 /                     &
               sqrt(dot_product(UnitVect_Fij(1)%d0,UnitVect_Fij(1)%d0))

          UnitVect_Fij(2)%d0(1) =                                       &
                       UnitVect_Fij(3)%d0(2) * UnitVect_Fij(1)%d0(3) -  &
                       UnitVect_Fij(3)%d0(3) * UnitVect_Fij(1)%d0(2)
          UnitVect_Fij(2)%d0(2) =                                       &
                       UnitVect_Fij(3)%d0(3) * UnitVect_Fij(1)%d0(1) -  &
                       UnitVect_Fij(3)%d0(1) * UnitVect_Fij(1)%d0(3)
          UnitVect_Fij(2)%d0(3) =                                       &
                       UnitVect_Fij(3)%d0(1) * UnitVect_Fij(1)%d0(2) -  &
                       UnitVect_Fij(3)%d0(2) * UnitVect_Fij(1)%d0(1)
          !write(out_unitp,*) 'ex_Fij',UnitVect_Fij(1)%d0
          !write(out_unitp,*) 'ey_Fij',UnitVect_Fij(2)%d0
          !write(out_unitp,*) 'ez_Fij',UnitVect_Fij(3)%d0
        END IF

        ! Riv = sqrt(dot_product(tab_Vect_Fi(iv_tot)%d0,tab_Vect_Fi(iv_tot)%d0)) !already calculated

        pz = dot_product(tab_Vect_Fi(iv_tot)%d0,UnitVect_Fi(3)%d0)
        px = dot_product(tab_Vect_Fi(iv_tot)%d0,UnitVect_Fi(1)%d0)
        py = dot_product(tab_Vect_Fi(iv_tot)%d0,UnitVect_Fi(2)%d0)

        ubetaiv = pz / Riv  ! cos(beta)
        betaiv = acos(ubetaiv)
        alphaiv = atan2(py,px)
        CALL dihedral_range(alphaiv,2) ! [0:2pi]

        IF (debug) write(out_unitp,*) 'R       : ',iv_Fi,':',Riv,Riv
        iQpoly  = iQpoly + 1
        i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(iQpoly)
        dnQpoly%d0(i_Qprim) = Riv

        ! loop on the vectors in frame Fij
        DO iv_Fij=2,nb_vect+1
!          write(out_unitp,*) 'iv_Fij,nb_vect+1',iv_Fij,nb_vect+1
!          write(out_unitp,*) 'shape tab_BFTransfo',shape(BFTransfo%tab_BFTransfo)
!          write(out_unitp,*) 'lbound tab_BFTransfo',lbound(BFTransfo%tab_BFTransfo)
!          write(out_unitp,*) 'ubound tab_BFTransfo',ubound(BFTransfo%tab_BFTransfo)

          CALL RecGet_Vec_Fi_For_poly(tab_Vect_Fi,nb_vect_tot,          &
                                      iv_tot,iv_Fij,UnitVect_Fij,       &
                                      dnQpoly,iQpoly,                   &
                                      BFTransfo%tab_BFTransfo(iv_Fij-1))
        END DO

        IF (iv_Fi == 2) THEN
          iQpoly  = iQpoly + 1
          IF (iQpoly <= dnQpoly%nb_var_vec) THEN
            i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(iQpoly)
            IF (BFTransfo%cos_th) THEN
              dnQpoly%d0(i_Qprim) = ubetaiv
            ELSE
              dnQpoly%d0(i_Qprim) = betaiv
            END IF
            IF (debug) write(out_unitp,*) 'beta (u): ',iv_Fi,':',betaiv*radTOdeg,ubetaiv
          END IF
        ELSE ! iv_Fi /= 2
          iQpoly  = iQpoly + 1
          IF (iQpoly <= dnQpoly%nb_var_vec) THEN
            i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(iQpoly)
            dnQpoly%d0(i_Qprim) = alphaiv
            IF (debug) write(out_unitp,*) 'alpha   : ',iv_Fi,':',alphaiv*radTOdeg,alphaiv
          END IF

          iQpoly = iQpoly + 1
          IF (iQpoly <= dnQpoly%nb_var_vec) THEN
            i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(iQpoly)
            IF (BFTransfo%cos_th) THEN
              dnQpoly%d0(i_Qprim) = ubetaiv
            ELSE
              dnQpoly%d0(i_Qprim) = betaiv
            END IF
            IF (debug) write(out_unitp,*) 'beta (u): ',iv_Fi,':',betaiv*radTOdeg,ubetaiv
          END IF
        END IF

        IF (nb_vect > 0) THEN
          ! for gamma we use ex_BF projected on the SF
          px = dot_product(UnitVect_Fij(1)%d0,UnitVect_Fi(3)%d0)
          py = dot_product(UnitVect_Fij(2)%d0,UnitVect_Fi(3)%d0)
          !write(out_unitp,*) 'for gamma, px,py',px,py
          gammaiv = atan2(py,-px)
          CALL dihedral_range(gammaiv,2) ! [0:2pi]

          iQpoly = iQpoly + 1
          IF (iQpoly <= dnQpoly%nb_var_vec) THEN
            i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(iQpoly)
            dnQpoly%d0(i_Qprim) = gammaiv
            IF (debug) write(out_unitp,*) 'gamma   : ',iv_Fi,':',gammaiv*radTOdeg,gammaiv
          END IF

        END IF

        CALL dealloc_dnSVM(UnitVect_Fij(1))
        CALL dealloc_dnSVM(UnitVect_Fij(2))
        CALL dealloc_dnSVM(UnitVect_Fij(3))

        IF (debug) write(out_unitp,*) '============================================='
      ELSE

!        write(out_unitp,*) 'ex_Fi',UnitVect_Fi(1)%d0
!        write(out_unitp,*) 'ey_Fi',UnitVect_Fi(2)%d0
!        write(out_unitp,*) 'ez_Fi',UnitVect_Fi(3)%d0

        IF (iv_Fi == 1) THEN
          write(out_unitp,*) ' ERROR in RecGet_Vec_Fi'
          write(out_unitp,*) ' iv_Fi=1 and frame=F is NOT possible'
          write(out_unitp,*) ' check the fortran!'
          STOP
        ELSE IF (iv_Fi == 2) THEN
          ! norm of the vector (distance)
          Riv = sqrt(dot_product(tab_Vect_Fi(iv_tot)%d0,tab_Vect_Fi(iv_tot)%d0))

          pz = dot_product(tab_Vect_Fi(iv_tot)%d0,UnitVect_Fi(3)%d0)
          uiv = pz / Riv  ! cos(thi)
          thiv = acos(uiv)

          iQpoly = iQpoly + 1
          i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(iQpoly)
          dnQpoly%d0(i_Qprim) = Riv
          IF (debug) write(out_unitp,*) 'R       : ',iv_Fi,':',Riv,Riv

          iQpoly = iQpoly + 1
          i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(iQpoly)

          IF (BFTransfo%cos_th) THEN
              dnQpoly%d0(i_Qprim) = uiv
          ELSE
              dnQpoly%d0(i_Qprim) = thiv
          END IF
          IF (debug) write(out_unitp,*) 'th (u)  : ',iv_Fi,':',thiv*radTOdeg,uiv

        ELSE
          IF (BFTransfo%cart) THEN
            iQpoly = iQpoly + 1
            i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(iQpoly)
            !dnQpoly%d0(i_Qprim) = tab_Vect_Fi(iv_tot)%d0(1) ! wrong, it is not projected in the BF
            dnQpoly%d0(i_Qprim) = dot_product(tab_Vect_Fi(iv_tot)%d0,UnitVect_Fi(1)%d0)

            iQpoly = iQpoly + 1
            i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(iQpoly)
            !dnQpoly%d0(i_Qprim) = tab_Vect_Fi(iv_tot)%d0(2)
            dnQpoly%d0(i_Qprim) = dot_product(tab_Vect_Fi(iv_tot)%d0,UnitVect_Fi(2)%d0)

            iQpoly = iQpoly + 1
            i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(iQpoly)
            !dnQpoly%d0(i_Qprim) = tab_Vect_Fi(iv_tot)%d0(3)
            dnQpoly%d0(i_Qprim) = dot_product(tab_Vect_Fi(iv_tot)%d0,UnitVect_Fi(3)%d0)

            IF (debug) write(out_unitp,*) 'x       : ',iv_Fi,':',tab_Vect_Fi(iv_tot)%d0(1)
            IF (debug) write(out_unitp,*) 'y       : ',iv_Fi,':',tab_Vect_Fi(iv_tot)%d0(2)
            IF (debug) write(out_unitp,*) 'z       : ',iv_Fi,':',tab_Vect_Fi(iv_tot)%d0(3)
          ELSE
            ! norm of the vector (distance)
            Riv = sqrt(dot_product(tab_Vect_Fi(iv_tot)%d0,tab_Vect_Fi(iv_tot)%d0))

            pz = dot_product(tab_Vect_Fi(iv_tot)%d0,UnitVect_Fi(3)%d0)
            uiv = pz / Riv  ! cos(thi)
            thiv = acos(uiv)

            px = dot_product(tab_Vect_Fi(iv_tot)%d0,UnitVect_Fi(1)%d0)
            py = dot_product(tab_Vect_Fi(iv_tot)%d0,UnitVect_Fi(2)%d0)
            phiv = atan2(py,px)
            CALL dihedral_range(phiv,2) ! [0:2pi]

            iQpoly = iQpoly + 1
            i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(iQpoly)
            dnQpoly%d0(i_Qprim) = Riv
            IF (debug) write(out_unitp,*) 'R       : ',iv_Fi,':',Riv,Riv

            iQpoly = iQpoly + 1
            i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(iQpoly)
            IF (BFTransfo%cos_th) THEN
              dnQpoly%d0(i_Qprim) = uiv
            ELSE
              dnQpoly%d0(i_Qprim) = thiv
            END IF
            IF (debug) write(out_unitp,*) 'th (u)  : ',iv_Fi,':',thiv*radTOdeg,uiv

            iQpoly = iQpoly + 1
            i_Qprim = BFTransfo%list_Qpoly_TO_Qprim(iQpoly)
            dnQpoly%d0(i_Qprim) = phiv
            IF (debug) write(out_unitp,*) 'phi     : ',iv_Fi,':',phiv*radTOdeg,phiv
          END IF

        END IF

      END IF

      END SUBROUTINE RecGet_Vec_Fi_For_poly



      !!@description: TODO
      !!@param: TODO
      RECURSIVE SUBROUTINE Rec_BFTransfo1TOBFTransfo2(BFTransfo1,BFTransfo2)
      TYPE (Type_BFTransfo),intent(in)    :: BFTransfo1
      TYPE (Type_BFTransfo),intent(inout) :: BFTransfo2

      integer :: i,n

!      -----------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!      logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='Rec_BFTransfo1TOBFTransfo2'
!      -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING RECURSIVE ',name_sub
        CALL RecWrite_BFTransfo(BFTransfo1,.FALSE.)
        CALL flush_perso(out_unitp)
      END IF


      CALL dealloc_BFTransfo(BFTransfo2)

      BFTransfo2%nb_vect                = BFTransfo1%nb_vect
      BFTransfo2%nb_vect_tot            = BFTransfo1%nb_vect_tot
      BFTransfo2%nb_var                 = BFTransfo1%nb_var
      BFTransfo2%nb_var_Rot             = BFTransfo1%nb_var_Rot

      BFTransfo2%num_vect_in_Frame      = BFTransfo1%num_vect_in_Frame
      BFTransfo2%num_vect_in_BF         = BFTransfo1%num_vect_in_BF



      BFTransfo2%num_Frame_in_BF        = BFTransfo1%num_Frame_in_BF
      BFTransfo2%num_Frame_in_Container = BFTransfo1%num_Frame_in_Container
      BFTransfo2%name_Frame             = BFTransfo1%name_Frame
      BFTransfo2%Frame                  = BFTransfo1%Frame
      BFTransfo2%BF                     = BFTransfo1%BF

      CALL FrameType1TOBFrameType2(BFTransfo1,BFTransfo2)

      BFTransfo2%cart                   = BFTransfo1%cart
      BFTransfo2%Li                     = BFTransfo1%Li
      BFTransfo2%cos_th                 = BFTransfo1%cos_th
      BFTransfo2%Def_cos_th             = BFTransfo1%Def_cos_th
      BFTransfo2%Euler(:)               = BFTransfo1%Euler(:)

      IF (associated(BFTransfo1%tab_num_Frame)) THEN
        n = size(BFTransfo1%tab_num_Frame)
        CALL alloc_array(BFTransfo2%tab_num_Frame,(/n/),                &
                        "BFTransfo2%tab_num_Frame",name_sub)
        BFTransfo2%tab_num_Frame(:) = BFTransfo1%tab_num_Frame(:)
      END IF

      IF (BFTransfo2%BF) THEN
        ! the lists are copied only for the BF frame.
        n = size(BFTransfo1%list_Qpoly_TO_Qprim)

        CALL alloc_array(BFTransfo2%list_Qpoly_TO_Qprim,(/n/),          &
                        "BFTransfo2%list_Qpoly_TO_Qprim",name_sub)
        BFTransfo2%list_Qpoly_TO_Qprim(:) = BFTransfo1%list_Qpoly_TO_Qprim(:)

        CALL alloc_array(BFTransfo2%list_Qprim_TO_Qpoly,(/n/),          &
                        "BFTransfo2%list_Qprim_TO_Qpoly",name_sub)
        BFTransfo2%list_Qprim_TO_Qpoly(:) = BFTransfo1%list_Qprim_TO_Qpoly(:)
      END IF

      ! they are true pointers, they are linked
      BFTransfo2%type_Qin => BFTransfo1%type_Qin
      BFTransfo2%name_Qin => BFTransfo1%name_Qin


      IF (BFTransfo2%Frame .AND. associated(BFTransfo1%tab_BFTransfo)) THEN
        n = size(BFTransfo1%tab_BFTransfo)
        CALL alloc_array(BFTransfo2%tab_BFTransfo,(/n/),                &
                        "BFTransfo2%tab_BFTransfo",name_sub)
        DO i=1,n
          !CALL RecWrite_BFTransfo(BFTransfo1%tab_BFTransfo(i),.FALSE.)
          BFTransfo2%tab_BFTransfo(i)%list_Qpoly_TO_Qprim => BFTransfo2%list_Qpoly_TO_Qprim
          BFTransfo2%tab_BFTransfo(i)%list_Qprim_TO_Qpoly => BFTransfo2%list_Qprim_TO_Qpoly

          CALL Rec_BFTransfo1TOBFTransfo2(BFTransfo1%tab_BFTransfo(i),  &
                                          BFTransfo2%tab_BFTransfo(i))
      END DO
      END IF

      IF (debug) THEN
        write(out_unitp,*)
        CALL RecWrite_BFTransfo(BFTransfo2,.FALSE.)
        write(out_unitp,*) 'END RECURSIVE ',name_sub
        CALL flush_perso(out_unitp)
      END IF

      END SUBROUTINE Rec_BFTransfo1TOBFTransfo2
!=======================================================================
!
!     Read Bunch transfo: Bunch of vectors to x
!
!=======================================================================
      !!@description: Read Bunch transfo: Bunch of vectors to x
      !!@param: TODO
      SUBROUTINE alloc_BunchTransfo(BunchTransfo)
      TYPE (Type_BunchTransfo), intent(inout) :: BunchTransfo

       character (len=*), parameter :: name_sub = 'alloc_BunchTransfo'


!      write(out_unitp,*) 'BEGINNING ',name_sub
!      write(out_unitp,*) 'nat,nb_var,ncart',BunchTransfo%nat,BunchTransfo%nb_var,BunchTransfo%ncart

       IF (BunchTransfo%nat < 3 .OR. BunchTransfo%nb_vect < 1 .OR.      &
           BunchTransfo%ncart < 9 .OR. BunchTransfo%nat_act < 2) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' wrong value of nat, nat_act, nb_vect or ncart',&
                                      BunchTransfo%nat,BunchTransfo%nat_act, &
                                      BunchTransfo%nb_vect,BunchTransfo%ncart
         write(out_unitp,*) ' CHECK the source !!'
         STOP
       END IF


       CALL alloc_array(BunchTransfo%ind_vect,(/5,BunchTransfo%nb_vect/),&
                       "BunchTransfo%ind_vect",name_sub)
       BunchTransfo%ind_vect(:,:) = 0

       CALL alloc_array(BunchTransfo%Mat_At_TO_centers,                 &
                           (/BunchTransfo%nat_act,BunchTransfo%nat/),   &
                       "BunchTransfo%Mat_At_TO_centers",name_sub)
       BunchTransfo%Mat_At_TO_centers(:,:) = 0

       CALL alloc_array(BunchTransfo%COM,                               &
                           (/BunchTransfo%nat_act,BunchTransfo%nat/),   &
                       "BunchTransfo%COM",name_sub)
       BunchTransfo%COM(:,:) = 0

       CALL alloc_array(BunchTransfo%A,                                 &
                     (/BunchTransfo%nb_vect+1,BunchTransfo%nb_vect+1/), &
                       "BunchTransfo%A",name_sub)
       BunchTransfo%A(:,:) = ZERO

       CALL alloc_array(BunchTransfo%A_inv,                             &
                   (/BunchTransfo%nb_vect+1,BunchTransfo%nb_vect+1/),   &
                       "BunchTransfo%A_inv",name_sub)
       BunchTransfo%A_inv(:,:) = ZERO

       CALL alloc_array(BunchTransfo%M_Tana,                            &
                       (/BunchTransfo%nb_vect,BunchTransfo%nb_vect/),   &
                       "BunchTransfo%M_Tana",name_sub)
       BunchTransfo%M_Tana(:,:) = ZERO

       CALL alloc_array(BunchTransfo%masses,(/BunchTransfo%ncart/),     &
                       "BunchTransfo%masses",name_sub)
       BunchTransfo%masses(:) = ZERO

       IF (associated(BunchTransfo%Z))  THEN
         CALL dealloc_array(BunchTransfo%Z,                             &
                           "BunchTransfo%Z",name_sub)
       END IF
       CALL alloc_array(BunchTransfo%Z,(/BunchTransfo%nat/),        &
                       "BunchTransfo%Z",name_sub)
       BunchTransfo%Z(:) = 0


       IF (associated(BunchTransfo%symbole))  THEN
         CALL dealloc_array(BunchTransfo%symbole,                       &
                           "BunchTransfo%symbole",name_sub)
       END IF
       CALL alloc_array(BunchTransfo%symbole,(/BunchTransfo%nat/),  &
                       "BunchTransfo%symbole",name_sub)
       BunchTransfo%symbole(:) = ""

       CALL alloc_array(BunchTransfo%masses_OF_At,(/BunchTransfo%nat_act/),&
                       "BunchTransfo%masses_OF_At",name_sub)
       BunchTransfo%masses_OF_At(:) = ZERO

!      write(out_unitp,*) 'END ',name_sub

      END SUBROUTINE alloc_BunchTransfo

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE dealloc_BunchTransfo(BunchTransfo)
       TYPE (Type_BunchTransfo),pointer, intent(inout) :: BunchTransfo

      integer :: err_mem,memory
       character (len=*), parameter :: name_sub='dealloc_BunchTransfo'

       !write(out_unitp,*) 'BEGINNING ',name_sub ; CALL flush_perso(out_unitp)

      IF (.NOT. associated(BunchTransfo)) RETURN


       IF (associated(BunchTransfo%ind_vect)) THEN
         CALL dealloc_array(BunchTransfo%ind_vect,                      &
                           "BunchTransfo%ind_vect",name_sub)
       END IF

       BunchTransfo%nb_X     = 0
       IF (associated(BunchTransfo%Mat_At_TO_centers))      THEN
         CALL dealloc_array(BunchTransfo%Mat_At_TO_centers,             &
                           "BunchTransfo%Mat_At_TO_centers",name_sub)
       END IF

       IF (associated(BunchTransfo%COM))      THEN
         CALL dealloc_array(BunchTransfo%COM,"BunchTransfo%COM",name_sub)
       END IF

       IF (associated(BunchTransfo%A))                THEN
         CALL dealloc_array(BunchTransfo%A,"BunchTransfo%A",name_sub)
       END IF

       IF (associated(BunchTransfo%A_inv))            THEN
         CALL dealloc_array(BunchTransfo%A_inv,"BunchTransfo%A_inv",name_sub)
       END IF
       IF (associated(BunchTransfo%M_Tana))           THEN
         CALL dealloc_array(BunchTransfo%M_Tana,"BunchTransfo%M_Tana",name_sub)
       END IF

       IF (associated(BunchTransfo%masses))           THEN
         CALL dealloc_array(BunchTransfo%masses,                        &
                           "BunchTransfo%masses",name_sub)
       END IF

       IF (associated(BunchTransfo%Z))  THEN
         CALL dealloc_array(BunchTransfo%Z,                             &
                           "BunchTransfo%Z",name_sub)
       END IF

       IF (associated(BunchTransfo%symbole))  THEN
         CALL dealloc_array(BunchTransfo%symbole,                       &
                           "BunchTransfo%symbole",name_sub)
       END IF

       IF (associated(BunchTransfo%masses_OF_At))     THEN
         CALL dealloc_array(BunchTransfo%masses_OF_At,                  &
                           "BunchTransfo%masses_OF_At",name_sub)
       END IF

       BunchTransfo%ncart     = 0
       BunchTransfo%ncart_act = 0
       BunchTransfo%nat0      = 0
       BunchTransfo%nat       = 0
       BunchTransfo%nat_act   = 0
       BunchTransfo%nb_var    = 0
       BunchTransfo%nb_vect   = 0


       nullify(BunchTransfo%type_Qin) ! TRUE pointer
       nullify(BunchTransfo%name_Qin) ! TRUE pointer


       deallocate(BunchTransfo,stat=err_mem)
       memory = 1
       CALL error_memo_allo(err_mem,-memory,'BunchTransfo',name_sub,'Type_BunchTransfo')
       nullify(BunchTransfo)

       !write(out_unitp,*) 'END dealloc_BunchTransfo'  ; CALL flush_perso(out_unitp)

      END SUBROUTINE dealloc_BunchTransfo

      SUBROUTINE Read_BunchTransfo(BunchTransfo,mendeleev)

      TYPE (Type_BunchTransfo),intent(inout) :: BunchTransfo
      TYPE (table_atom), intent(in)          :: mendeleev


      integer :: i,j,i_at,icf,icf1,iat1
      integer :: nb_vect,nat_dum,nb_at

      real (kind=Rkind)  :: at,Mtot
      character (len=Name_len) :: name1_at,name2_at

      integer, pointer :: tab_iAtTOiCart(:)
      real (kind=Rkind), pointer :: weight_vect(:)

      integer, parameter :: max_atG = 100
      integer :: tab_At_TO_G(max_atG),GAt
      integer :: tab_At_recenter(max_atG)
      NAMELIST /recenterG / tab_At_TO_G,tab_At_recenter,GAt

      !--------------------------------------------------------
      real (kind=Rkind), allocatable :: A(:,:),M_Tana(:,:) ! for F. Gatti and M. Ndong
      real (kind=Rkind), allocatable :: BB(:,:),B(:,:)     ! for F. Gatti and M. Ndong
      real (kind=Rkind), allocatable :: BG(:)              ! for F. Gatti and M. Ndong
      real (kind=Rkind), allocatable :: masses(:)          ! for F. Gatti and M. Ndong
      !--------------------------------------------------------


      integer :: err_io
      character (len=*), parameter :: name_sub = "Read_BunchTransfo"
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      !--------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
       END IF
      !---------------------- Tana -------------------------------------
      !--------------------------------------------------------
      write(out_unitp,*) '  Read the vector connectivities'
      write(out_unitp,*) '  nb_vect',BunchTransfo%nb_vect
      write(out_unitp,*) '  nat',BunchTransfo%nat


      BunchTransfo%type_Qin(:) = 1

      BunchTransfo%nat_act = BunchTransfo%nb_vect + 1
      CALL alloc_BunchTransfo(BunchTransfo)

      nullify(tab_iAtTOiCart)
      CALL alloc_array(tab_iAtTOiCart,(/BunchTransfo%nat/),             &
                      "tab_iAtTOiCart",name_sub)
      tab_iAtTOiCart(:) = 0

      nullify(weight_vect)
      CALL alloc_array(weight_vect,(/BunchTransfo%nb_vect/),            &
                      "weight_vect",name_sub)
      weight_vect(:) = ZERO

      !--------------------------------------------------------
      CALL alloc_NParray(M_Tana,                                        &
                     (/BunchTransfo%nb_vect+1,BunchTransfo%nb_vect+1/), &
                      "M_Tana",name_sub)
      CALL alloc_NParray(BB,                                            &
                  (/BunchTransfo%nb_vect+1,BunchTransfo%nb_vect+1/),    &
                      "BB",name_sub)
      CALL alloc_NParray(A,                                             &
                  (/BunchTransfo%nb_vect+1,BunchTransfo%nb_vect+1/),    &
                      "A",name_sub)
      CALL alloc_NParray(B,(/BunchTransfo%nat,BunchTransfo%nb_vect+1/), &
                      "B",name_sub)
      CALL alloc_NParray(BG,(/BunchTransfo%nb_vect+1/),"BG",name_sub)
      CALL alloc_NParray(masses,(/BunchTransfo%nat/),"masses",name_sub)
      B(:,:) = ZERO
      masses(:) = ZERO
      !---------------------- Tana -------------------------------------

      !--------------------------------------------------------
      !-- Read Bunch Vector -----------------------------------
      !--------------------------------------------------------

       nb_vect = BunchTransfo%nb_vect

       BunchTransfo%nat_act = 0
       ! The total center of mass will be at BunchTransfo%nat
       nat_dum = BunchTransfo%nat-1

       IF (nb_vect < 1) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' the number of vector is < 1',nb_vect
         STOP
       END IF

       ! for the first atom (mass = 0.)
       icf = func_ic(nat_dum)
       tab_iAtTOiCart(1)    = icf
       BunchTransfo%masses(icf+0:icf+2)  = ZERO

      !--------------------------------------------------------
       B(1,:)    = ZERO
       masses(1) = ZERO
      !--------------------------------------------------------

      write(out_unitp,*) '---------------------------------------------'
      write(out_unitp,*) '                                             '
      write(out_unitp,*) 'vector,V i: |----.-------->                  '
      write(out_unitp,*) '            A----X--------B                  '
      write(out_unitp,*) '                                             '
      write(out_unitp,*) '            | -w |          : A = X -w.V     '
      write(out_unitp,*) '                 |   1-w  | : B = X + (1-w).V'
      write(out_unitp,*) ' Remark: V = B-A                             '
      write(out_unitp,*) ' Read: #vect:i, #linked_at:X, weight:w, ',    &
                         'new_at:A and B'
      write(out_unitp,*) '---------------------------------------------'

      i_at = 1
      DO i=1,nb_vect
        CALL make_nameQ(BunchTransfo%name_Qin(3*i-2),"Xbunch",i)
        CALL make_nameQ(BunchTransfo%name_Qin(3*i-1),"Ybunch",i)
        CALL make_nameQ(BunchTransfo%name_Qin(3*i-0),"Zbunch",i)
        read(in_unitp,*,IOSTAT=err_io) BunchTransfo%ind_vect(1,i),iat1, &
                                       weight_vect(i),name1_at,name2_at

        IF (err_io /= 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  while the vector definition (bunch)'
          write(out_unitp,*) ' Check your data !!'
          STOP
        END IF

        IF (print_level > 1) THEN
          write(out_unitp,*) BunchTransfo%ind_vect(1,i),iat1,           &
                                        weight_vect(i),name1_at,name2_at
        END IF

        icf1 = tab_iAtTOiCart(iat1)
        BunchTransfo%ind_vect(2,i) = icf1
        IF (icf1 == 0) THEN
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) '  The vector (',                       &
             BunchTransfo%ind_vect(1,i),                              &
            ') is linked to an undefined atom (',iat1,').'

           write(out_unitp,*) i_at,                                   &
                 '  atoms have been of already defined'
           STOP
        END IF

        i_at                       = i_at +1
        BunchTransfo%Z(i_at)       = -1
        BunchTransfo%symbole(i_at) = name1_at
        at =get_mass_Tnum(mendeleev,Z=BunchTransfo%Z(i_at),name=name1_at)
        !write(out_unitp,*) 'atom1:',i,BunchTransfo%Z(i),at
        IF (at > ZERO) THEN
          BunchTransfo%nat_act = BunchTransfo%nat_act + 1
          icf = func_ic(BunchTransfo%nat_act)
        ELSE
          nat_dum = nat_dum - 1
          icf = func_ic(nat_dum)
        END IF

        !- Tana -----------------------------------
        B(i_at,:) = B(iat1,:)
        B(i_at,BunchTransfo%ind_vect(1,i)) =                          &
                                 B(i_at,BunchTransfo%ind_vect(1,i)) - &
                                                        weight_vect(i)
        masses(i_at) = at
        !- Tana -----------------------------------

        BunchTransfo%masses(icf+0:icf+2)  = at
        BunchTransfo%ind_vect(3,i)        = icf
        tab_iAtTOiCart(i_at) = icf
        write(out_unitp,*) '# vect,# at1,icf1,mass1',i,i_at,icf,at

        i_at                       = i_at +1
        BunchTransfo%Z(i_at)       = -1
        BunchTransfo%symbole(i_at) = name2_at
        at =get_mass_Tnum(mendeleev,Z=BunchTransfo%Z(i_at),name=name2_at)
        !write(out_unitp,*) 'atom2:',i,BunchTransfo%Z(i),at
        IF (at > ZERO) THEN
          BunchTransfo%nat_act = BunchTransfo%nat_act + 1
          icf = func_ic(BunchTransfo%nat_act)
        ELSE
          nat_dum = nat_dum - 1
          icf = func_ic(nat_dum)
        END IF

        !- Tana -----------------------------------
        B(i_at,:) = B(iat1,:)
        B(i_at,BunchTransfo%ind_vect(1,i)) =                          &
                                 B(i_at,BunchTransfo%ind_vect(1,i)) + &
                                                  ONE - weight_vect(i)
        masses(i_at) = at
        !- Tana -----------------------------------

        BunchTransfo%masses(icf+0:icf+2)  = at
        BunchTransfo%ind_vect(4,i)        = icf
        tab_iAtTOiCart(i_at)              = icf
        write(out_unitp,*) '# vect,# at2,icf2,mass2',i,i_at,icf,at

      END DO
      IF (BunchTransfo%nat_act /= BunchTransfo%nb_vect+1) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) '  The number of true atoms (not dummy) is ' // &
               'wrong ',BunchTransfo%nat_act
        write(out_unitp,*) '  It MUST be equal to ',BunchTransfo%nb_vect+1
        write(out_unitp,*) '  Probably, you have too many or not ' //   &
                                                  'enough dummy atoms'

        write(out_unitp,*) '  CHECK your data!!'
        STOP
      END IF

      !--------------------------------------------------------
      write(out_unitp,*) '---------------------------------------------'
      write(out_unitp,*) '                                             '
      write(out_unitp,*) '   ADD centers of masses'
      write(out_unitp,*) '     nb_G:',BunchTransfo%nb_G


      DO i=1,BunchTransfo%nb_G
        tab_At_TO_G(:)     = 0
        tab_At_recenter(:) = 0
        GAt                = 0

        read(in_unitp,recenterG,IOSTAT=err_io)
        IF (err_io /= 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  while the nemalist "recenterG"'
          write(out_unitp,*) ' Check your data !!'
          STOP
        END IF
         !write(out_unitp,recenterG)
         IF (count(tab_At_recenter > 0) == 0) tab_At_recenter(:) = tab_At_TO_G(:)
         write(out_unitp,*) '    G:',i,'Gat:',Gat

         nb_at = count(tab_At_TO_G > 0)
         IF (nb_at > BunchTransfo%nat) STOP 'wrong nat!!!'
         write(out_unitp,*) '    tab_At_TO_G     :',tab_At_TO_G(1:nb_at)

         nb_at = count(tab_At_recenter > 0)
         IF (nb_at > BunchTransfo%nat) STOP 'wrong nat!!!'

         write(out_unitp,*) '    tab_At_recenter :',tab_At_recenter(1:nb_at)

         IF (GAt < 1) THEN
            write(out_unitp,*) 'ERROR in ',name_sub
            write(out_unitp,*) 'GAt is < 0',Gat
            STOP
         END IF

         !---------------------- Tana -----------------------------------
         ! center of masses: i
         BG(:) = ZERO
         Mtot  = ZERO
         !write(out_unitp,*) 'i,G',i
         DO j=1,count(tab_At_TO_G(:) > 0)
            i_at = tab_At_TO_G(j)
            Mtot = Mtot + masses(i_at)
            !write(out_unitp,*) 'i,j,i_At',i,j,iAt,
            BG(:) = BG(:) + masses(i_at) * B(i_at,:)
         END DO

         ! Add the position of GAt
         !write(out_unitp,*) 'i,iAt of Gat',i,GAt
         BG(:) = B(GAt,:) - BG(:)/Mtot

         ! recenter the atom: i
         DO j=1,count(tab_At_recenter(:) > 0)
            i_at = tab_At_recenter(j)
            !write(out_unitp,*) 'i,j,iAt',i,j,iAt
            B(i_at,:) = B(i_at,:) + BG(:)
         END DO
        !---------------------- Tana -----------------------------------

      END DO
      write(out_unitp,*) '   END ADD centers of masses                 '
      write(out_unitp,*) '                                             '
      write(out_unitp,*) '---------------------------------------------'

      !ncart_act number of active cartesian coordinates (without dummy atom and G)
      BunchTransfo%ncart_act = 3 * BunchTransfo%nat_act
      IF (debug) write(out_unitp,*) 'nat_act,ncart_act',                &
                            BunchTransfo%nat_act,BunchTransfo%ncart_act
      !-------------------------------------------------------
      !-- end Read bunch of vectors  -------------------------
      !-------------------------------------------------------
      !-------------------------------------------------------

      !---------------------------------------------------------------
      !---------------------- Tana -----------------------------------
      !---------------------------------------------------------------
      !Here, we get the matrix M used in the analitycal calulation
      ! of the KEO with polysherical.
      ! see Gatti and Iung, Phys. Rep. v484, p1-69, 2009 : p10, just after the fig3

      !Center of mass (CM)
      B(BunchTransfo%nat,:)    = ZERO
      Mtot                     = ZERO
      DO i=1,BunchTransfo%nat-1
        B(BunchTransfo%nat,:)  = B(BunchTransfo%nat,:) + B(i,:)*masses(i)
        Mtot                   = Mtot + masses(i)
      END DO
      B(BunchTransfo%nat,:)    = B(BunchTransfo%nat,:) / Mtot
      masses(BunchTransfo%nat) = ZERO

      !All the B(:,i) are recentered with respect to the CM
      DO i=1,BunchTransfo%nat-1
        B(i,:) = B(i,:)-B(BunchTransfo%nat,:)
      END DO
      B(BunchTransfo%nat,:) = ZERO

      !Add the position of the CM with respect the LF
      B(:,BunchTransfo%nb_vect+1) = ONE

      !Use the relations and sort the dummy atoms (centers of mass)
      i_at = 0
      DO i=1,BunchTransfo%nat
        IF (masses(i) > 0) THEN
          i_at                       = i_at + 1
          BB(i_at,:)                 = B(i,:)
          BunchTransfo%A_inv(i_at,:) = B(i,:)
        END IF
      END DO
      IF (print_level > 1) THEN
        write(out_unitp,*) 'masses',masses(:)
        write(out_unitp,*) 'A_inv'
        CALL Write_Mat(BB,out_unitp,5)
      END IF

      !calculation of A (not mass weighted)
      CALL inv_m1_TO_m2(BB,A,BunchTransfo%nb_vect+1,0,ZERO) ! not SVD

      IF (print_level > 1) THEN
        write(out_unitp,*) 'A'
        CALL Write_Mat(A,out_unitp,5)
      END IF

      !Use the massweighted relation and sort the dummy atoms (centers of mass)
      i_at = 0
      DO i=1,BunchTransfo%nat
        IF (masses(i) > 0) THEN
          i_at                       = i_at + 1
          BB(i_at,:)                 = B(i,:) * sqrt(masses(i))
          BunchTransfo%A_inv(i_at,:) = B(i,:)
        END IF
      END DO

      !calculation of A (mass weighted)
      CALL inv_m1_TO_m2(BB,A,BunchTransfo%nb_vect+1,0,ZERO) ! not SVD

      !calculation of M (included the M of the center of mass)
      M_Tana = matmul(A,transpose(A))
      BunchTransfo%M_Tana(:,:) = M_Tana(1:nb_vect,1:nb_vect)
      write(out_unitp,*) 'M_Tana (without the center-of-mass contribution)'
      CALL Write_Mat(BunchTransfo%M_Tana,out_unitp,5)

      CALL dealloc_array(tab_iAtTOiCart,"tab_iAtTOiCart",name_sub)
      CALL dealloc_array(weight_vect,   "weight_vect",   name_sub)
      CALL dealloc_NParray(M_Tana,        "M_Tana",        name_sub)
      CALL dealloc_NParray(A,             "A",             name_sub)
      CALL dealloc_NParray(B,             "B",             name_sub)
      CALL dealloc_NParray(BB,            "BB",            name_sub)
      CALL dealloc_NParray(BG,            "BG",            name_sub)
      CALL dealloc_NParray(masses,        "masses",        name_sub)

      !--------------------------------------------------------
      IF (debug) write(out_unitp,*) 'END ',name_sub
      !--------------------------------------------------------

      END SUBROUTINE Read_BunchTransfo

      SUBROUTINE Read2_BunchTransfo(BunchTransfo,mendeleev,with_vect)

      TYPE (Type_BunchTransfo),intent(inout) :: BunchTransfo
      TYPE (table_atom), intent(in)          :: mendeleev


      logical :: with_vect
      integer :: i,j,i_at,icf,icf1,iat1,iv,iAtA,iAtB,jat,iAt0
      integer :: nb_vect,nat_dum,nb_at,nb_At_FOR_G

      real (kind=Rkind)  :: at,Mtot,MtotG,wX
      character (len=Name_len) :: name1_at,name2_at,type_dummyX
      character (len=Name_len), pointer :: name_at(:)


      integer, parameter :: max_atG = 1000
      integer :: tab_At_TO_G(max_atG)
      integer :: tab_At_TO_X(max_atG)
      real (kind=Rkind) :: Mat_At_TO_centers(BunchTransfo%nat_act)

      NAMELIST /recenterG / tab_At_TO_G
      NAMELIST /dummyX /    tab_At_TO_X,wX,type_dummyX


      !--------------------------------------------------------
      integer :: err_io
      character (len=*), parameter :: name_sub = "Read2_BunchTransfo"
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      !--------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
      END IF
      !--------------------------------------------------------
        write(out_unitp,*) '  Read the atoms'
        write(out_unitp,*) '  nb_vect',BunchTransfo%nb_vect
        write(out_unitp,*) '  nat',BunchTransfo%nat

        ! allocation of the variables:
        CALL alloc_BunchTransfo(BunchTransfo)

        BunchTransfo%type_Qin(:) = 1

        !--------------------------------------------------------
        ! first read the atoms (not dummy)
        !--------------------------------------------------------
        BunchTransfo%Mat_At_TO_centers(:,:) = ZERO

        nb_vect = BunchTransfo%nb_vect

        nb_at = BunchTransfo%nat_act + BunchTransfo%nb_G + BunchTransfo%nb_X
        write(out_unitp,*) '  nat_act',BunchTransfo%nat_act
        nullify(name_at)
        CALL alloc_array(name_at,(/nb_at/),"name_at",name_sub)

        read(in_unitp,*,IOSTAT=err_io) (name_at(i),i=1,BunchTransfo%nat_act)

        IF (err_io /= 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) '  while reading the list of atoms or masses'
          write(out_unitp,*) ' Check your data !!'
          STOP
        END IF

        Mtot = ZERO
        DO i_at=1,BunchTransfo%nat_act
          BunchTransfo%Z(i_at) = -1
          BunchTransfo%symbole(i_at) = name_at(i_at)
          at =get_mass_Tnum(mendeleev,Z=BunchTransfo%Z(i_at),name=name_at(i_at))
          icf = func_ic(i_at)

          BunchTransfo%Mat_At_TO_centers(i_at,i_at) = ONE
          BunchTransfo%COM(i_at,i_at)               = ONE

          BunchTransfo%masses_OF_At(i_at)  = at
          BunchTransfo%masses(icf+0:icf+2) = at
          !write(out_unitp,*) 'atom:',i_at,icf,BunchTransfo%Z(i_at),at
          Mtot = Mtot + at
          IF (at == ZERO) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' The readed atom cannot be dummy !'
            write(out_unitp,*) 'atom:',i_at,icf,BunchTransfo%Z(i_at),at
            STOP
          END IF

        END DO

        CALL dealloc_array(name_at,"name_at",name_sub)
        !write(out_unitp,*) 'Masses: ',BunchTransfo%masses(:)

        ! for the centers of mass
        IF (BunchTransfo%nb_G > 0)                                      &
              write(out_unitp,'(a,i0)') '  Read the centers of mass',BunchTransfo%nb_G

        DO i_at=BunchTransfo%nat_act+1,BunchTransfo%nat_act+BunchTransfo%nb_G

          tab_At_TO_G(:)     = 0

          read(in_unitp,recenterG,IOSTAT=err_io)
          IF (err_io /= 0) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) '  while the namelist "recenterG"'
            write(out_unitp,*) ' Check your data !!'
            STOP
          END IF
          IF (debug) write(out_unitp,recenterG)

          CALL Add_DummyG(BunchTransfo%Mat_At_TO_centers,               &
                          BunchTransfo%COM,i_at,                        &
                          BunchTransfo%masses_OF_At,MtotG,tab_At_TO_G)

        END DO

        IF (BunchTransfo%nb_X > 0)                                      &
                        write(out_unitp,*) '  Read the dummy atoms, X'
        DO i_at=BunchTransfo%nat_act+BunchTransfo%nb_G+1,BunchTransfo%nat_act+BunchTransfo%nb_G+BunchTransfo%nb_x
          tab_At_TO_X(:)     = 0
          wX                 = ZERO
          type_dummyX        = 'COM'

          read(in_unitp,dummyX,IOSTAT=err_io)
          IF (err_io < 0) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) '  while reading the namelist "dummyX"'
            write(out_unitp,*) ' end of file or end of record'
            write(out_unitp,*) ' Probably, nb_X is to large in the namelist "Coord_transfo"'
            write(out_unitp,*) '   or you have forgotten the namelist.'
            write(out_unitp,*) ' Check your data !!'
            STOP
          END IF
          IF (err_io > 0) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) '  while reading the namelist "dummyX"'
            write(out_unitp,*) ' Probably, some arguments of namelist are wrong.'
            write(out_unitp,*) ' Check your data !!'
            STOP
          END IF
          IF (debug) write(out_unitp,dummyX)

          CALL string_uppercase_TO_lowercase(type_dummyX)

          SELECT CASE (type_dummyX)

          CASE ("com","g")

            CALL Add_DummyG(BunchTransfo%Mat_At_TO_centers,             &
                            BunchTransfo%COM,i_at,                      &
                            BunchTransfo%masses_OF_At,MtotG,tab_At_TO_X)
          CASE ("wx")

            iAtA = tab_At_TO_X(1)
            iAtB = tab_At_TO_X(2)
            IF (iAtA < 1 .OR. iAtA > i_at-1 .OR. iAtB<1 .OR. iAtB > i_at-1) THEN
              write(out_unitp,*) ' ERROR in ',name_sub
              write(out_unitp,*) ' the atom indexes iAtA, iAtB are out of range: [1,',i_at-1,']'
              write(out_unitp,*) ' tab_At_TO_X(:) ',tab_At_TO_X(1:count(tab_At_TO_X > 0))
              STOP
            END IF
            BunchTransfo%Mat_At_TO_centers(:,i_at) =                    &
                          wX * BunchTransfo%Mat_At_TO_centers(:,iAtB) + &
                  (ONE - wX) * BunchTransfo%Mat_At_TO_centers(:,iAtA)

          CASE ("radau")

            CALL Add_DummyRadau(BunchTransfo%Mat_At_TO_centers,i_at,    &
                                BunchTransfo%masses_OF_At,tab_At_TO_X)

          CASE DEFAULT
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' Wrong type_dummyX: ',trim(type_dummyX)
            write(out_unitp,*) ' The possibilities are:'
            write(out_unitp,*) '  COM : center of mass'
            write(out_unitp,*) '  wX : weighted average of two centers'
            write(out_unitp,*) '  Radau : confocal point of Radau coordinates'
            STOP
          END SELECT

        END DO

        !     -------------------------------------------------------
        !--------------------------------------------------------------

        IF (with_vect) THEN
          write(out_unitp,*) '  Read the vector connectivities'

          !--------------------------------------------------------------
          !  read the vecteurs
          !--------------------------------------------------------------
          DO iv=1,BunchTransfo%nb_vect
            read(in_unitp,*,IOSTAT=err_io) BunchTransfo%ind_vect(3:4,iv)
            IF (err_io /= 0) THEN
              write(out_unitp,*) ' ERROR in ',name_sub
              write(out_unitp,*) '  while reading the vecors indices.'
              write(out_unitp,*) ' Check your data !!'
              STOP
            END IF
            IF (BunchTransfo%ind_vect(3,iv) < 1     .OR.                &
                BunchTransfo%ind_vect(3,iv) > nb_at .OR.                &
                BunchTransfo%ind_vect(4,iv) < 1     .OR.                &
                BunchTransfo%ind_vect(4,iv) > nb_at) THEN

               write(out_unitp,*) ' ERROR in ',name_sub
               write(out_unitp,*) '  The vector (',iv,') indexes are out-of-range'
               write(out_unitp,*) '   BunchTransfo%ind_vect(3:4,iv)',BunchTransfo%ind_vect(3:4,iv)
               write(out_unitp,*) ' Check your data !!'
               STOP
            END IF
          END DO
          !---------------------------------------------------------------
        END IF

        ! ncart_act number of active cartesian coordinates (without dummy atom and G)
        BunchTransfo%ncart_act = 3 * BunchTransfo%nat_act
        write(out_unitp,*) '  nat_act,ncart_act',                       &
                            BunchTransfo%nat_act,BunchTransfo%ncart_act

        IF (print_level > 1) THEN
          write(out_unitp,*) '  Mat_At_TO_centers'
          CALL Write_Mat(BunchTransfo%Mat_At_TO_centers,out_unitp,5)

          write(out_unitp,*) '  COM'
          DO i=1,ubound(BunchTransfo%COM,dim=2)
            write(out_unitp,*) '  Atom list for COM',i,': ',BunchTransfo%COM(:,i)
          END DO
        END IF

      !--------------------------------------------------------
      IF (debug) write(out_unitp,*) 'END ',name_sub
      !--------------------------------------------------------
      END SUBROUTINE Read2_BunchTransfo

      SUBROUTINE Add_DummyG(Mat_At_TO_centers,COM,iG,masses_OF_At,MtotG,tab_At_TO_G)

      real(kind=Rkind), intent(inout) :: Mat_At_TO_centers(:,:)
      real(kind=Rkind), intent(inout) :: COM(:,:)
      real(kind=Rkind), intent(inout) :: masses_OF_At(:)
      integer, intent(in)             :: tab_At_TO_G(:)
      integer, intent(in)             :: iG
      real(kind=Rkind), intent(inout) :: MtotG


      integer          :: i,j,jat,nat_act

      !--------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = "Add_DummyG"
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      !--------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
       END IF
      !--------------------------------------------------------

       nat_act = ubound(Mat_At_TO_centers,dim=1)
       Mat_At_TO_centers(:,iG) = ZERO
       COM(:,iG) = ZERO

       DO j=1,count(tab_At_TO_G > 0)
         jat = tab_At_TO_G(j)
         IF (jat < 1 .OR. jat > iG-1) THEN
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) ' The index of the center-of-mass is out of range: [1,',iG-1,']'
           write(out_unitp,*) ' tab_At_TO_G(:) ',tab_At_TO_G(1:count(tab_At_TO_G > 0))
           STOP
         END IF

         DO i=1,nat_act
           IF (Mat_At_TO_centers(i,iG) /= ZERO .AND. Mat_At_TO_centers(i,jat) /= ZERO ) THEN
             write(out_unitp,*) ' ERROR in ',name_sub
             write(out_unitp,*) ' Two partial centers of mass are developped on the same atoms!'
             write(out_unitp,*) ' tab_At_TO_G(:) ',tab_At_TO_G(1:count(tab_At_TO_G > 0))
             write(out_unitp,*) ' Check your data !!'
             STOP
           ELSE IF (Mat_At_TO_centers(i,jat) /= ZERO ) THEN
             Mat_At_TO_centers(i,iG) = masses_OF_At(i)
             COM(i,iG) = ONE
           END IF
         END DO
       END DO

       MtotG = sum(Mat_At_TO_centers(:,iG))

       IF (abs(MtotG) == ZERO) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' MtotG = 0'

         DO j=1,ubound(Mat_At_TO_centers,dim=2)
           write(out_unitp,*) 'MtotG,Mat_At_TO_centers(:,j):',Mat_At_TO_centers(:,j)
         END DO
         write(out_unitp,*) ' Check your data !!'
         STOP
       END IF
       Mat_At_TO_centers(:,iG) = Mat_At_TO_centers(:,iG) / MtotG

      !--------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'iG,MtotG,Mat_At_TO_centers(:,iG):',iG,MtotG,Mat_At_TO_centers(:,iG)
        write(out_unitp,*) 'END ',name_sub
      END IF
      !--------------------------------------------------------

      END SUBROUTINE Add_DummyG

      SUBROUTINE Add_DummyRadau(Mat_At_TO_centers,iX,masses_OF_At,tab_At_TO_X)

      real(kind=Rkind), intent(inout) :: Mat_At_TO_centers(:,:)
      real(kind=Rkind), intent(inout) :: masses_OF_At(:)
      integer, intent(in)            :: tab_At_TO_X(:)
      integer, intent(in)            :: iX

      integer           :: i,j,jat,iAt0,nb_At_FOR_G
      real (kind=Rkind) :: MtotG,M_iAt0,wX
      integer           :: tab_At_TO_G(size(tab_At_TO_X))
      real (kind=Rkind) :: Mat_At_TO_G(ubound(Mat_At_TO_centers,dim=1))
      real (kind=Rkind) :: COM(ubound(Mat_At_TO_centers,dim=1),ubound(Mat_At_TO_centers,dim=2))


      !--------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = "Add_DummyRadau"
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      !--------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
       END IF
      !--------------------------------------------------------

            iAt0 = tab_At_TO_X(1)
            IF (iAt0 < 1 .OR. iAt0 > iX-1) THEN
              write(out_unitp,*) ' ERROR in ',name_sub
              write(out_unitp,*) ' the atom index iAt0 is out of range: [1,',iX-1,']'
              write(out_unitp,*) ' tab_At_TO_X(:) ',tab_At_TO_X(1:count(tab_At_TO_X > 0))
              STOP
            END IF

            !M_iAt0 = masses_OF_At(iAt0)
            M_iAt0 = get_Mass(Mat_At_TO_centers,iAt0,masses_OF_At,tab_At_TO_G)
            IF (M_iAt0 == ZERO) THEN
              write(out_unitp,*) ' ERROR in ',name_sub
              write(out_unitp,*) ' the mass of the atom iAt0 is zero'
              STOP
            END IF
            IF (debug) write(out_unitp,*) 'iAt0,M_iAt0,Mat_At_TO_iAt0(:)', &
                                     iAt0,M_iAt0,Mat_At_TO_centers(:,iAt0)

            nb_At_FOR_G = count(tab_At_TO_X > 0) - 1
            IF (nb_At_FOR_G < 2) THEN
                write(out_unitp,*) ' ERROR in ',name_sub
                write(out_unitp,*) ' For Radau coordinates, you need at list three centers'
                write(out_unitp,*) ' tab_At_TO_X(:) ',tab_At_TO_X(1:count(tab_At_TO_X > 0))
                STOP
            END IF

            ! first the center of mass of centers : tab_At_TO_X(i) i > 1
            ! It is added in place of iX (the canonical point for the Radau)
            CALL Add_DummyG(Mat_At_TO_centers,COM,iX,masses_OF_At,MtotG,&
                            tab_At_TO_X(2:nb_At_FOR_G+1))
            ! .... then it is transfered in Mat_At_TO_G
            Mat_At_TO_G(:) = Mat_At_TO_centers(:,iX)
            IF (debug) write(out_unitp,*) 'MtotG,Mat_At_TO_G(:)',MtotG,Mat_At_TO_G(:)

            ! Then, wX and the dummy atom for the Radau coordinates
            ! from F. T. Smith, PRL, 1980, vol 45, p1157
            !   OD . CD = BD^2
            ! O : Mat_At_TO_centers(:,iAt0)
            ! D : Mat_At_TO_centers(:)  (center of mass partial)
            ! C : center of mass total
            ! B : canonical point : Mat_At_TO_centers(:,i_at)
            ! si OD = R alors CD = mass(O) / Mtot * R
            ! OC + CD = OD = R => OC = R - CD = R - sqrt(R.R. MO/Mtot) = R(1-sqrt(MO/Mtot))
            ! C = OC + O = O + wX*OD = O*(1-wX) + wX*D


            wX = ONE - sqrt(M_iAt0/(MtotG+M_iAt0))


            Mat_At_TO_centers(:,iX) = ZERO
            Mat_At_TO_centers(:,iX) = wX * Mat_At_TO_G(:) +             &
                                  (ONE - wX) * Mat_At_TO_centers(:,iAt0)

      !--------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'iX,xW',iX,wX
        write(out_unitp,*) 'Mat_At_TO_centers(:,iX)',Mat_At_TO_centers(:,iX)
        write(out_unitp,*) 'END ',name_sub
      END IF
      !--------------------------------------------------------

      END SUBROUTINE Add_DummyRadau

      FUNCTION get_Mass(Mat_At_TO_centers,jat,masses_OF_At,tab_At_TO_G)

      real(kind=Rkind)                :: get_Mass
      real(kind=Rkind), intent(inout) :: Mat_At_TO_centers(:,:)
      real(kind=Rkind), intent(inout) :: masses_OF_At(:)
      integer, intent(in)             :: tab_At_TO_G(:)
      integer                         :: jat



      real(kind=Rkind) :: MtotG
      real(kind=Rkind) :: Vec_At_TO_centers(ubound(Mat_At_TO_centers,dim=1))
      integer          :: i,j,nat_act,nat

      !--------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub = "get_Mass"
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      !--------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
       END IF
      !--------------------------------------------------------

       nat_act = ubound(Mat_At_TO_centers,dim=1)
       nat     = ubound(Mat_At_TO_centers,dim=2)
       Vec_At_TO_centers(:) = ZERO

       IF (jat < 1 .OR. jat > nat) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' The index of the center-of-mass is out of range: [1,',nat,']'
         STOP
       END IF

       DO i=1,nat_act
         IF (Mat_At_TO_centers(i,jat) /= ZERO ) THEN
           Vec_At_TO_centers(i) = masses_OF_At(i)
         END IF
       END DO

       MtotG = sum(Vec_At_TO_centers(:))

       IF (abs(MtotG) == ZERO) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         DO j=1,ubound(Mat_At_TO_centers,dim=2)
           write(out_unitp,*) 'MtotG,Mat_At_TO_centers(:,j):',Mat_At_TO_centers(:,j)
         END DO
         write(out_unitp,*) ' The mass of the center is ZERO!!'
         write(out_unitp,*) ' Check you data !!'

         STOP
       END IF
       Vec_At_TO_centers(:) = Vec_At_TO_centers(:) / MtotG

       get_Mass = MtotG

      !--------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'MtotG,Vec_At_TO_centers(:):',MtotG,Vec_At_TO_centers(:)
        write(out_unitp,*) 'END ',name_sub
      END IF
      !--------------------------------------------------------

      END FUNCTION get_Mass

      SUBROUTINE M_Tana_FROM_Bunch2Transfo(BunchTransfo)

      TYPE(type_bunchtransfo) :: BunchTransfo
      real (kind=Rkind)       :: Mtot


      !--------------------------------------------------------
      real (kind=Rkind), allocatable :: A(:,:)             ! for F. Gatti and M. Ndong
      real (kind=Rkind), allocatable :: Asave(:,:)         ! for F. Gatti and M. Ndong
      real (kind=Rkind), allocatable :: M_Tana(:,:)        ! for F. Gatti and M. Ndong
      real (kind=Rkind), allocatable :: BB(:,:)            ! for F. Gatti and M. Ndong

      integer :: i_at,iv
      integer :: iAtA,iAtB
      integer :: err_mem,memory

      !--------------------------------------------------------
      character (len=*), parameter :: name_sub='M_Tana_FROM_Bunch2Transfo'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      !--------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
       END IF
      !--------------------------------------------------------

      CALL alloc_NParray(A,(/BunchTransfo%nat_act,BunchTransfo%nat_act/), &
                        "A",name_sub)
      CALL alloc_NParray(Asave,                                         &
                         (/BunchTransfo%nat_act,BunchTransfo%nat_act/), &
                        "Asave",name_sub)
      CALL alloc_NParray(M_Tana,                                        &
                        (/BunchTransfo%nat_act,BunchTransfo%nat_act/),  &
                        "M_Tana",name_sub)
      CALL alloc_NParray(BB,(/BunchTransfo%nat_act,BunchTransfo%nat_act/),&
                        "BB",name_sub)

      Mtot = sum(BunchTransfo%masses(1:BunchTransfo%ncart_act:3))

      IF (debug) THEN
        DO iatA=1,ubound(BunchTransfo%Mat_At_TO_centers,dim=2)
          write(out_unitp,*) 'iat,Mat_At_TO_centers',iatA,              &
                                BunchTransfo%Mat_At_TO_centers(:,iAtA)
        END DO
      END IF

      ! The A matrix (not mass-weighted) to get the B one
      A(:,:) = ZERO
      DO iv=1,BunchTransfo%nb_vect
        iAtA = BunchTransfo%ind_vect(3,iv)
        iAtB = BunchTransfo%ind_vect(4,iv)
        IF (debug) write(out_unitp,*) 'iAtA,iAtB',iAtA,iAtB
        A(iv,:) = BunchTransfo%Mat_At_TO_centers(:,iAtB) -            &
                  BunchTransfo%Mat_At_TO_centers(:,iAtA)
      END DO
      DO i_at=1,BunchTransfo%nat_act
        A(BunchTransfo%nb_vect+1,i_at) = BunchTransfo%masses_OF_At(i_at) / Mtot
      END DO
      Asave(:,:) = A(:,:)
      BunchTransfo%A(:,:) = A(:,:)
      IF (debug) THEN
        write(out_unitp,*) 'A'
        CALL Write_Mat(A,out_unitp,5)
      END IF
      ! calculation of A_inv (B)
      CALL inv_m1_TO_m2(Asave,BB,BunchTransfo%nat_act,0,ZERO) ! not SVD

      BunchTransfo%A_inv(:,:) = ZERO
      DO i_at=1,BunchTransfo%nat_act
        BunchTransfo%A_inv(i_at,:) = BB(i_at,:)
      END DO
      IF (debug) THEN
        write(out_unitp,*) 'A_inv'
        CALL Write_Mat(BB,out_unitp,5)
      END IF


      ! for the M-matrix, we use mass-weighted A matrix
      DO iv=1,BunchTransfo%nb_vect+1
        A(iv,:) = A(iv,:) / sqrt(BunchTransfo%masses_OF_At(1:BunchTransfo%nat_act))
      END DO

      ! calculation of M (included the M of the center of mass)
      M_Tana = matmul(A,transpose(A))
      BunchTransfo%M_Tana(:,:) = M_Tana(1:BunchTransfo%nb_vect,1:BunchTransfo%nb_vect)

      write(out_unitp,*) '  M_Tana (without the center-of-mass contribution)'
      CALL Write_Mat(BunchTransfo%M_Tana,out_unitp,5)

      CALL dealloc_NParray(A,     "A",     name_sub)
      CALL dealloc_NParray(Asave, "Asave", name_sub)
      CALL dealloc_NParray(M_Tana,"M_Tana",name_sub)
      CALL dealloc_NParray(BB,    "BB",    name_sub)

      !--------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
      !--------------------------------------------------------
      END SUBROUTINE M_Tana_FROM_Bunch2Transfo

      SUBROUTINE Write_BunchTransfo(BunchTransfo)
      TYPE (Type_BunchTransfo), pointer, intent(in) :: BunchTransfo

      integer :: i,nb_at
      character (len=*), parameter :: name_sub='Write_BunchTransfo'


      IF (.NOT. associated(BunchTransfo)) RETURN

      write(out_unitp,*) 'BEGINNING ',name_sub

      write(out_unitp,*) 'ncart_act,ncart',                             &
                  BunchTransfo%ncart_act,BunchTransfo%ncart

      write(out_unitp,*) 'nat_act,nat0,nat,',                           &
                 BunchTransfo%nat_act,BunchTransfo%nat0,BunchTransfo%nat

      write(out_unitp,*) 'nb_var',BunchTransfo%nb_var
      write(out_unitp,*) 'nb_vect',BunchTransfo%nb_vect


     write(6,*) 'A(iv,iat): Xat=>Vect'
     CALL Write_Mat(BunchTransfo%A,out_unitp,5)

     write(6,*) 'A_inv(iat,iv): Vect=>Xat'
     CALL Write_Mat(BunchTransfo%A_inv,out_unitp,5)


      write(out_unitp,*) '---------------------------------------------'
      write(out_unitp,*) '------------ Tana ---------------------------'
      write(out_unitp,*) '---------------------------------------------'
      write(out_unitp,*) 'M_Tana (without center-of-mass contribution)'
      CALL Write_Mat(BunchTransfo%M_Tana,out_unitp,5)
      write(out_unitp,*) '---------------------------------------------'
      write(out_unitp,*) '---------------------------------------------'

      write(out_unitp,*) 'END ',name_sub

      END SUBROUTINE Write_BunchTransfo

      SUBROUTINE calc_BunchTransfo(dnQvect,dnx,BunchTransfo,nderiv,inTOout)

      TYPE (Type_dnVec), intent(inout)    :: dnQvect,dnx
      TYPE (Type_BunchTransfo),pointer, intent(in) :: BunchTransfo
      integer, intent(in)                 :: nderiv
      logical, intent(in)                 :: inTOout


      integer           ::iv,ivect_iv,i,i_at,ivect_at
      integer           :: iG,ixyz_G,ixyz_At,j,jat
      real(kind=Rkind)  :: weight,MtotG


      !--------------------------------------------------------
      integer :: nderiv_debug = 0
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='calc_BunchTransfo'
      !--------------------------------------------------------
      IF (.NOT. associated(BunchTransfo)) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'BunchTransfo is not associated!'
        write(out_unitp,*) ' CHECK the fortran!!'
        STOP
      END IF

      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nderiv',nderiv
        write(out_unitp,*)
        CALL Write_BunchTransfo(BunchTransfo)
        write(out_unitp,*) 'dnx: '
        CALL write_dnx(1,dnx%nb_var_vec,dnx,nderiv_debug)
        write(out_unitp,*) 'dnQvect: '
        CALL Write_dnSVM(dnQvect,nderiv_debug)
        CALL flush_perso(out_unitp)
      END IF
      !--------------------------------------------------------


      IF (inTOout) THEN
        ! Xvect => Xat
        CALL Set_ZERO_TO_dnSVM(dnx)

        DO iv=1,BunchTransfo%nb_vect
          ivect_iv = func_ic(iv)
          DO i_at=1,BunchTransfo%nat_act
            weight = BunchTransfo%A_inv(i_at,iv)

            IF (abs(weight) > ZERO) THEN
              ivect_at = func_ic(i_at)  ! 3*i-2
              CALL sub_dnVec1_wADDTO_dnVec2(dnQvect,ivect_iv,weight,    &
                                            dnx,ivect_at,ONE,3,nderiv)
            END IF
            !write(out_unitp,*) 'i_at,iv,A_inv',i_at,iv,weight
          END DO
        END DO

      ELSE

        ! Xat => Xvect
        CALL Set_ZERO_TO_dnSVM(dnQvect)

        !write(out_unitp,*) 'BunchTransfo%A',BunchTransfo%A
        DO iv=1,BunchTransfo%nb_vect
          ivect_iv = func_ic(iv)
          DO i_at=1,BunchTransfo%nat_act
            weight = BunchTransfo%A(iv,i_at)

            IF (abs(weight) > ZERO) THEN
              ivect_at = func_ic(i_at)  ! 3*i-2
              !write(out_unitp,*) 'i_at,weight',i_at,weight
              CALL sub_dnVec1_wADDTO_dnVec2(dnx,ivect_at,weight,        &
                                          dnQvect,ivect_iv,ONE,3,nderiv)
            END IF
          END DO
        END DO
      END IF
      !--------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'dnx: '
        CALL write_dnx(1,dnx%nb_var_vec,dnx,nderiv_debug)
        write(out_unitp,*) 'dnQvect: '
        CALL Write_dnSVM(dnQvect,nderiv_debug)
        write(out_unitp,*) 'END ',name_sub
        write(out_unitp,*)
        CALL flush_perso(out_unitp)
      END IF
      !--------------------------------------------------------

      END SUBROUTINE calc_BunchTransfo

      SUBROUTINE BunchTransfo1TOBunchTransfo2(BunchTransfo1,BunchTransfo2)
      TYPE (Type_BunchTransfo),pointer, intent(in)    :: BunchTransfo1
      TYPE (Type_BunchTransfo),pointer, intent(inout) :: BunchTransfo2


      !--------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='BunchTransfo1TOBunchTransfo2'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      !--------------------------------------------------------
      IF (.NOT. associated(BunchTransfo1)) RETURN

      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL Write_BunchTransfo(BunchTransfo1)
      END IF
      !--------------------------------------------------------

      CALL dealloc_BunchTransfo(BunchTransfo2)
      allocate(BunchTransfo2,stat=err_mem)
      memory = 1
      CALL error_memo_allo(err_mem,memory,'BunchTransfo2',name_sub,'Type_BunchTransfo')

      BunchTransfo2%ncart     = BunchTransfo1%ncart
      BunchTransfo2%ncart_act = BunchTransfo1%ncart_act

      BunchTransfo2%nat0      = BunchTransfo1%nat0
      BunchTransfo2%nat       = BunchTransfo1%nat
      BunchTransfo2%nat_act   = BunchTransfo1%nat_act

      BunchTransfo2%nb_var    = BunchTransfo1%nb_var

      BunchTransfo2%nb_vect   = BunchTransfo1%nb_vect
      BunchTransfo2%nb_G      = BunchTransfo1%nb_G
      BunchTransfo2%nb_X      = BunchTransfo1%nb_X

      CALL alloc_BunchTransfo(BunchTransfo2)

      IF (associated(BunchTransfo1%masses)) THEN
        BunchTransfo2%masses(:) = BunchTransfo1%masses(:)
      END IF
      IF (associated(BunchTransfo1%Z)) THEN
        BunchTransfo2%Z(:) = BunchTransfo1%Z(:)
      END IF
      IF (associated(BunchTransfo1%symbole)) THEN
        BunchTransfo2%symbole(:) = BunchTransfo1%symbole(:)
      END IF

      IF (associated(BunchTransfo1%ind_vect)) THEN
        BunchTransfo2%ind_vect(:,:) = BunchTransfo1%ind_vect(:,:)
      END IF
      IF (associated(BunchTransfo1%M_Tana)) THEN
        BunchTransfo2%M_Tana(:,:)       = BunchTransfo1%M_Tana(:,:)
      END IF
      IF (associated(BunchTransfo1%A_inv)) THEN
        BunchTransfo2%A_inv(:,:)        = BunchTransfo1%A_inv(:,:)
      END IF
      IF (associated(BunchTransfo1%A)) THEN
        BunchTransfo2%A(:,:)            = BunchTransfo1%A(:,:)
      END IF

      IF (associated(BunchTransfo1%Mat_At_TO_centers)) THEN
        BunchTransfo2%Mat_At_TO_centers(:,:) = BunchTransfo1%Mat_At_TO_centers(:,:)
      END IF

      !--------------------------------------------------------
      IF (debug) THEN
        CALL Write_BunchTransfo(BunchTransfo2)
        write(out_unitp,*) 'END ',name_sub
      END IF
      !-------------------------------------------------------
      END SUBROUTINE BunchTransfo1TOBunchTransfo2

   !! @description: Defines the spherical unit vectors ei.
   !!               Eq 10, from Ndong et al. J. Chem. Phys. V136, pp034107, 2012
   !! @param:       V      The unit vector (type: vec_sum_opnd)
   !! @param:       theta       an elementary op  which contains the needed
   !!                          information on \theta or u coordinate
   !! @param:       phi       an elementary op  which contains the needed
   !!                          information on \phi  coordinate
   SUBROUTINE get_unit_vector_Ei(V, Qvec, index_v, cart)
     type(vec_sum_opnd),      intent(inout)      :: V
     type(opel),              intent(in)         :: Qvec(3) ! R,th,phi
     integer,                 intent(in)         :: index_v
     logical,                 intent(in)         :: cart

     character (len = *), parameter :: routine_name='get_unit_vector_Ei'

     call allocate_op(V, 3)

     if (index_v == 1) then

       V%vec_sum(1) = czero
       V%vec_sum(2) = czero
       V%vec_sum(3) = cone

     else if (index_v == 2) then ! it deals automatically idq=3 or -3

       V%vec_sum(1) = get_sin(Qvec(2))
       V%vec_sum(2) = czero
       V%vec_sum(3) = get_cos(Qvec(2))

     else
       IF (cart) THEN
         V%vec_sum(1) = czero
         V%vec_sum(2) = czero
         V%vec_sum(3) = czero
       ELSE
         V%vec_sum(1) =  get_sin(Qvec(2)) * get_cos(Qvec(3))
         V%vec_sum(2) =  get_sin(Qvec(2)) * get_sin(Qvec(3))
         V%vec_sum(3) =  get_cos(Qvec(2))
       END IF

     end if

   END SUBROUTINE get_unit_vector_Ei

      END MODULE mod_BunchPolyTransfo

