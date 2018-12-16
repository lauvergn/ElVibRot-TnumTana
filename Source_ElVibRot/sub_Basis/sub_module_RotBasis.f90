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
      MODULE mod_RotBasis_Param
      USE mod_system
      use mod_nDindex, only: alloc_nparray, dealloc_nparray
      IMPLICIT NONE

      PRIVATE

        TYPE RotBasis_Param
          integer :: Jrot   = -1                    !  J value
          integer :: nb_Rot = 0                     !  size of the operator (2*Jrot+1)

          integer :: nb_term = 0
          integer :: tab_der_TO_iterm(-3:0,-3:0)    ! i1 or i2 =-3,-2,-1   => Jz, Jy, Jx
                                                    ! ex: -2,-1            => JyJx+JxJy operator
                                                    ! ex: -2, 0 or 0,-2    => Jy operator

          integer, allocatable :: tab_iterm_TO_der(:,:) !  ...(2,nb_term)

          real(kind=Rkind), allocatable :: tab_RotOp(:,:,:)  ! tab_RotOp(nb_Rot,nb_Rot,0:nb_term)
                                                    ! tab_RotOp(:,:,0) is not used but it needs when tab_der_TO_iterm(0,0)=0

        END TYPE RotBasis_Param

        INTERFACE assignment (=)
          MODULE PROCEDURE RotBasis_Param2TORotBasis_Param1
        END INTERFACE

      PUBLIC  RotBasis_Param, assignment (=), alloc_RotBasis_Param, &
              dealloc_RotBasis_Param, Init_RotBasis_Param, Write_RotBasis_Param, &
              RotBasis_Param2TORotBasis_Param1

      CONTAINS

      ! symmetrized version => JxJy+JyJx ...
      SUBROUTINE alloc_RotBasis_Param(RotBasis_Para,Jrot)

      TYPE (RotBasis_Param), intent(inout) :: RotBasis_Para
      integer, intent (in)        :: Jrot

      integer :: J1,J2,iterm
!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='alloc_RotBasis_Param'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (Jrot < 0 .OR. RotBasis_Para%nb_Rot > 0) RETURN
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Jrot',Jrot
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------

      RotBasis_Para%Jrot    = Jrot
      RotBasis_Para%nb_Rot  = Jrot+Jrot+1

      RotBasis_Para%nb_term = 0
      DO J1=-3,-1
      DO J2=J1,-1
         RotBasis_Para%nb_term   = RotBasis_Para%nb_term + 1

         RotBasis_Para%tab_der_TO_iterm(J1,J2) = RotBasis_Para%nb_term
         RotBasis_Para%tab_der_TO_iterm(J2,J1) = RotBasis_Para%nb_term
      END DO
      END DO
      DO J1=-3,-1
         RotBasis_Para%nb_term  = RotBasis_Para%nb_term + 1

         RotBasis_Para%tab_der_TO_iterm(J1,0) = RotBasis_Para%nb_term
         RotBasis_Para%tab_der_TO_iterm(0,J1) = RotBasis_Para%nb_term
      END DO
      RotBasis_Para%tab_der_TO_iterm(0,0) = 0


      CALL alloc_NParray(RotBasis_Para%tab_iterm_TO_der,                &
                                         (/ 2,RotBasis_Para%nb_term /), &
                        "RotBasis_Para%tab_iterm_TO_der",name_sub, (/ 1,0 /))

      DO J1=-3,0
      DO J2=J1,0
         iterm = RotBasis_Para%tab_der_TO_iterm(J1,J2)
         RotBasis_Para%tab_iterm_TO_der(:,iterm) = (/ J1,J2 /)
      END DO
      END DO

      CALL alloc_NParray(RotBasis_Para%tab_RotOp,                       &
 (/ RotBasis_Para%nb_Rot,RotBasis_Para%nb_Rot,RotBasis_Para%nb_term /), &
                        "RotBasis_Para%tab_RotOp",name_sub, (/ 1,1,0 /))
      RotBasis_Para%tab_RotOp(:,:,:) = ZERO



!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        CALL Write_RotBasis_Param(RotBasis_Para)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE alloc_RotBasis_Param

      SUBROUTINE alloc_RotBasis_Param_old(RotBasis_Para,Jrot)

      TYPE (RotBasis_Param), intent(inout) :: RotBasis_Para
      integer, intent (in)        :: Jrot

      integer :: J1,J2
!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='alloc_RotBasis_Param_old'
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (Jrot < 0 .OR. RotBasis_Para%nb_Rot > 0) RETURN
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Jrot',Jrot
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------

      RotBasis_Para%Jrot    = Jrot
      RotBasis_Para%nb_Rot  = Jrot+Jrot+1

      RotBasis_Para%nb_term = 0
      DO J1=-3,-1
      DO J2=-3,-1
         RotBasis_Para%nb_term   = RotBasis_Para%nb_term + 1

         RotBasis_Para%tab_der_TO_iterm(J1,J2) = RotBasis_Para%nb_term
      END DO
      END DO
      DO J1=-3,-1
         RotBasis_Para%nb_term  = RotBasis_Para%nb_term + 1

         RotBasis_Para%tab_der_TO_iterm(J1,0) = RotBasis_Para%nb_term
         RotBasis_Para%tab_der_TO_iterm(0,J1) = RotBasis_Para%nb_term
      END DO
      RotBasis_Para%tab_der_TO_iterm(0,0) = 0


      CALL alloc_NParray(RotBasis_Para%tab_RotOp,                       &
 (/ RotBasis_Para%nb_Rot,RotBasis_Para%nb_Rot,RotBasis_Para%nb_term /), &
                        "RotBasis_Para%tab_RotOp",name_sub, (/ 1,1,0 /))
      RotBasis_Para%tab_RotOp(:,:,:) = ZERO



!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        CALL Write_RotBasis_Param(RotBasis_Para)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE alloc_RotBasis_Param_old


      SUBROUTINE dealloc_RotBasis_Param(RotBasis_Para)

      TYPE (RotBasis_Param), intent(inout) :: RotBasis_Para

!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='dealloc_RotBasis_Param'
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------

      RotBasis_Para%Jrot    = -1
      RotBasis_Para%nb_Rot  = 0
      RotBasis_Para%nb_term = 0

      RotBasis_Para%tab_der_TO_iterm(:,:) = 0

      IF (allocated(RotBasis_Para%tab_RotOp)) THEN
        CALL dealloc_NParray(RotBasis_Para%tab_iterm_TO_der,            &
                            "RotBasis_Para%tab_iterm_TO_der",name_sub)
      END IF

      IF (allocated(RotBasis_Para%tab_RotOp)) THEN
        CALL dealloc_NParray(RotBasis_Para%tab_RotOp,                   &
                            "RotBasis_Para%tab_RotOp",name_sub)
      END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        CALL Write_RotBasis_Param(RotBasis_Para)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE dealloc_RotBasis_Param

      SUBROUTINE Init_RotBasis_Param(RotBasis_Para,Jrot)

      TYPE (RotBasis_Param), intent(inout) :: RotBasis_Para
      integer, intent (in)        :: Jrot

      integer :: J1,J2,K,iterm,iterm1,iterm2
      complex (kind=Rkind) :: s2
      complex (kind=Rkind) :: is2

      complex(kind=Rkind), allocatable :: PWang(:,:)     ! PWang(nb_Rot,nb_Rot)
      complex(kind=Rkind), allocatable :: PWangDag(:,:)  ! PWangDag(nb_Rot,nb_Rot)
      complex(kind=Rkind), allocatable :: test_inv(:,:)  ! PWangDag(nb_Rot,nb_Rot)

      complex(kind=Rkind), allocatable :: Jx(:,:)    ! Jx(-J:J,-J:J)
      complex(kind=Rkind), allocatable :: Jy(:,:)    ! Jx(-J:J,-J:J)
      complex(kind=Rkind), allocatable :: Jz(:,:)    ! Jx(-J:J,-J:J)
      complex(kind=Rkind), allocatable :: Ji(:,:,:)  ! Ji(-J:J,-J:J,-3:-1)

!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Init_RotBasis_Param'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (Jrot < 0) RETURN
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Jrot',Jrot
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------
      s2  = cmplx( sqrt(HALF), ZERO,kind=Rkind )
      is2 = cmplx( ZERO, sqrt(HALF),kind=Rkind )


      CALL alloc_RotBasis_Param(RotBasis_Para,Jrot)

!---------------------------------------------------------------------
      ! -1- First the Wang basis
      CALL alloc_NParray(PWang,(/ Jrot,RotBasis_Para%nb_Rot /),         &
                        'PWang',name_sub,(/-Jrot,1 /) )

      CALL alloc_NParray(PWangDag,(/ RotBasis_Para%nb_Rot,Jrot /),      &
                        'PWangDag',name_sub,(/1,-Jrot /) )


      PWangDag(:,:) = ZERO
      PWangDag(1,0) = ONE
      DO K=1,Jrot
        IF (mod(K,2) == 0) THEN
          PWangDag(K+1, K) =  s2
          PWangDag(K+1,-K) =  s2
        ELSE
          PWangDag(K+1, K) =  s2
          PWangDag(K+1,-K) = -s2
        END IF
      END DO

      DO K=1,Jrot
        IF (mod(K,2) == 0) THEN
          PWangDag(K+Jrot+1, K) =  is2
          PWangDag(K+Jrot+1,-K) = -is2
        ELSE
          PWangDag(K+Jrot+1, K) =  is2
          PWangDag(K+Jrot+1,-K) =  is2
        END IF
      END DO
      PWang(:,:) = Transpose(conjg(PWangDag))

      IF (debug) THEN
        write(out_unitp,*) 'PWangDag (in line)'
        CALL Write_Mat(PWangDag,out_unitp,5)

        write(out_unitp,*) 'PWang (in column)'
        CALL Write_Mat(PWang,out_unitp,5)
      END IF

      !CALL alloc_NParray(test_inv,(/ RotBasis_Para%nb_Rot,RotBasis_Para%nb_Rot /),         &
      !                  'test_inv',name_sub)
      !test_inv = matmul(PWang,PWangDag)
      !write(out_unitp,*) 'PWang . PWangDag'
      !CALL Write_Mat(test_inv,out_unitp,5)
      !test_inv = matmul(PWangDag,PWang)
      !write(out_unitp,*) 'PWangDag . PWang'
      !CALL Write_Mat(test_inv,out_unitp,5)
      !CALL dealloc_NParray(test_inv,'test_inv',name_sub)


!---------------------------------------------------------------------
      ! -2- Jx,Jy,Jz in the I JKM > basis
      CALL alloc_NParray(Ji,(/ Jrot,Jrot,-1 /),'Ji',name_sub,(/-Jrot,-Jrot,-3 /) )

      CALL alloc_NParray(Jx,(/ Jrot,Jrot /),'Jx',name_sub,(/-Jrot,-Jrot /) )
      Jx(:,:) = ZERO
      CALL alloc_NParray(Jy,(/ Jrot,Jrot /),'Jy',name_sub,(/-Jrot,-Jrot /) )
      Jy(:,:) = ZERO
      CALL alloc_NParray(Jz,(/ Jrot,Jrot /),'Jz',name_sub,(/-Jrot,-Jrot /) )
      Jz(:,:) = ZERO


      DO K=-Jrot,Jrot
        Jz(K,K) = cmplx(K,0,kind=Rkind)

        IF (K < Jrot) THEN
          Jx(K,K+1) =      cmplx(funcP_JK(Jrot,K),kind=Rkind)
          Jy(K,K+1) = -EYE*cmplx(funcP_JK(Jrot,K),kind=Rkind)
        END IF

        IF (K > -Jrot) THEN
          Jx(K,K-1) =      cmplx(funcM_JK(Jrot,K),kind=Rkind)
          Jy(K,K-1) =  EYE*cmplx(funcM_JK(Jrot,K),kind=Rkind)
        END IF

      END DO

      !write(out_unitp,*)
      !write(out_unitp,*) 'Re Jx'
      !CALL Write_Mat(real(Jx),out_unitp,5)
      !write(out_unitp,*) 'Jy'
      !CALL Write_Mat(aimag(Jy),out_unitp,5)
      !write(out_unitp,*) 'Jz'
      !CALL Write_Mat(real(Jz),out_unitp,5)

!---------------------------------------------------------------------
      ! -3- Jx,Jy,Jz in the Wang basis
      Ji(:,:,-1) = matmul(PWangDag,matmul(Jx,PWang))
      Ji(:,:,-2) = matmul(PWangDag,matmul(Jy,PWang))
      Ji(:,:,-3) = matmul(PWangDag,matmul(Jz,PWang))

      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'Im(Jx), SumRe(Jx)',sum(abs(real(Ji(:,:,-1))))
        CALL Write_Mat(aimag(Ji(:,:,-1)),out_unitp,5)
        write(out_unitp,*) 'Im(Jy), SumRe(Jy)',sum(abs(real(Ji(:,:,-2))))
        CALL Write_Mat(aimag(Ji(:,:,-2)),out_unitp,5)
        write(out_unitp,*) 'Im(Jz), SumRe(Jz)',sum(abs(real(Ji(:,:,-3))))
        CALL Write_Mat(aimag(Ji(:,:,-3)),out_unitp,5)
      END IF
      !CALL alloc_NParray(test_inv,(/ RotBasis_Para%nb_Rot,RotBasis_Para%nb_Rot /),         &
      !                  'test_inv',name_sub)
      !test_inv = matmul(Jx,Jx)+matmul(Jy,Jy)+matmul(Jz,Jz)
      !write(out_unitp,*) 'J^2'
      !CALL Write_Mat(test_inv,out_unitp,5)
      !CALL dealloc_NParray(test_inv,'test_inv',name_sub)

!---------------------------------------------------------------------
      ! -4- Transfert Jx,Jy .... (Jx*Jy+JyJx) ... in RotBasis_Para%tab_RotOp
      DO J1=-3,-1
        iterm = RotBasis_Para%tab_der_TO_iterm(J1,J1)
        RotBasis_Para%tab_RotOp(:,:,iterm) =                            &
                                     real(matmul(Ji(:,:,J1),Ji(:,:,J1)))
        DO J2=J1+1,-1
         iterm = RotBasis_Para%tab_der_TO_iterm(J1,J2)
         RotBasis_Para%tab_RotOp(:,:,iterm) =                           &
                                  real(matmul(Ji(:,:,J1),Ji(:,:,J2))) + &
                                  real(matmul(Ji(:,:,J2),Ji(:,:,J1)))
        END DO
      END DO
      ! A minus sign is added because the Hamitonian has i*Ji and the Wang basis are imaginary
      DO J1=-3,-1
         iterm = RotBasis_Para%tab_der_TO_iterm(J1,0)
         RotBasis_Para%tab_RotOp(:,:,iterm) = -aimag(Ji(:,:,J1))
      END DO

      ! Id
      iterm = RotBasis_Para%tab_der_TO_iterm(0,0)
      RotBasis_Para%tab_RotOp(:,:,iterm) = ZERO
      DO K=1,RotBasis_Para%nb_Rot
        RotBasis_Para%tab_RotOp(K,K,iterm) = ONE
      END DO



!---------------------------------------------------------------------
      CALL dealloc_NParray(PWangDag,'PWangDag',name_sub)
      CALL dealloc_NParray(PWang,'PWang',name_sub)
      CALL dealloc_NParray(Jx,'Jx',name_sub)
      CALL dealloc_NParray(Jy,'Jy',name_sub)
      CALL dealloc_NParray(Jz,'Jz',name_sub)
      CALL dealloc_NParray(Ji,'Ji',name_sub)


!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        CALL Write_RotBasis_Param(RotBasis_Para)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------

!STOP 'Init_RotBasis_Param'
!---------------------------------------------------------------------

      END SUBROUTINE Init_RotBasis_Param

      SUBROUTINE Init_RotBasis_Param_old(RotBasis_Para,Jrot)

      TYPE (RotBasis_Param), intent(inout) :: RotBasis_Para
      integer, intent (in)        :: Jrot

      integer :: J1,J2,K,iterm,iterm1,iterm2
      complex (kind=Rkind) :: s2
      complex (kind=Rkind) :: is2

      complex(kind=Rkind), allocatable :: PWang(:,:)     ! PWang(nb_Rot,nb_Rot)
      complex(kind=Rkind), allocatable :: PWangDag(:,:)  ! PWangDag(nb_Rot,nb_Rot)
      complex(kind=Rkind), allocatable :: test_inv(:,:)  ! PWangDag(nb_Rot,nb_Rot)

      complex(kind=Rkind), allocatable :: Jx(:,:)    ! Jx(-J:J,-J:J)
      complex(kind=Rkind), allocatable :: Jy(:,:)    ! Jx(-J:J,-J:J)
      complex(kind=Rkind), allocatable :: Jz(:,:)    ! Jx(-J:J,-J:J)
      complex(kind=Rkind), allocatable :: Ji(:,:,:)  ! Ji(-J:J,-J:J,-3:-1)

!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Init_RotBasis_Param_old'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (Jrot <= 0 .OR. RotBasis_Para%nb_Rot > 0) RETURN
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Jrot',Jrot
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------
      s2  = cmplx( sqrt(HALF), ZERO,kind=Rkind )
      is2 = cmplx( ZERO, sqrt(HALF),kind=Rkind )


      CALL alloc_RotBasis_Param(RotBasis_Para,Jrot)

!---------------------------------------------------------------------
      ! -1- First the Wang basis
      CALL alloc_NParray(PWang,(/ Jrot,RotBasis_Para%nb_Rot /),         &
                        'PWang',name_sub,(/-Jrot,1 /) )

      CALL alloc_NParray(PWangDag,(/ RotBasis_Para%nb_Rot,Jrot /),      &
                        'PWangDag',name_sub,(/1,-Jrot /) )


      PWangDag(:,:) = CZERO
      PWangDag(1,0) = CONE
      DO K=1,Jrot
        IF (mod(K,2) == 0) THEN
          PWangDag(K+1, K) =  s2
          PWangDag(K+1,-K) =  s2
        ELSE
          PWangDag(K+1, K) =  s2
          PWangDag(K+1,-K) = -s2
        END IF
      END DO

      DO K=1,Jrot
        IF (mod(K,2) == 0) THEN
          PWangDag(K+Jrot+1, K) =  is2
          PWangDag(K+Jrot+1,-K) = -is2
        ELSE
          PWangDag(K+Jrot+1, K) =  is2
          PWangDag(K+Jrot+1,-K) =  is2
        END IF
      END DO
      PWang(:,:) = Transpose(conjg(PWangDag))

      IF (debug) THEN
        write(out_unitp,*) 'PWangDag (in line)'
        CALL Write_Mat(PWangDag,out_unitp,5)

        write(out_unitp,*) 'PWang (in column)'
        CALL Write_Mat(PWang,out_unitp,5)
      END IF

      !CALL alloc_NParray(test_inv,(/ RotBasis_Para%nb_Rot,RotBasis_Para%nb_Rot /),         &
      !                  'test_inv',name_sub)
      !test_inv = matmul(PWang,PWangDag)
      !write(out_unitp,*) 'PWang . PWangDag'
      !CALL Write_Mat(test_inv,out_unitp,5)
      !test_inv = matmul(PWangDag,PWang)
      !write(out_unitp,*) 'PWangDag . PWang'
      !CALL Write_Mat(test_inv,out_unitp,5)
      !CALL dealloc_NParray(test_inv,'test_inv',name_sub)


!---------------------------------------------------------------------
      ! -2- Jx,Jy,Jz in the I JKM > basis
      CALL alloc_NParray(Ji,(/ Jrot,Jrot,-1 /),'Ji',name_sub,(/-Jrot,-Jrot,-3 /) )

      CALL alloc_NParray(Jx,(/ Jrot,Jrot /),'Jx',name_sub,(/-Jrot,-Jrot /) )
      Jx(:,:) = CZERO
      CALL alloc_NParray(Jy,(/ Jrot,Jrot /),'Jy',name_sub,(/-Jrot,-Jrot /) )
      Jy(:,:) = CZERO
      CALL alloc_NParray(Jz,(/ Jrot,Jrot /),'Jz',name_sub,(/-Jrot,-Jrot /) )
      Jz(:,:) = CZERO


      DO K=-Jrot,Jrot
        Jz(K,K) = cmplx(K,0,kind=Rkind)

        IF (K < Jrot) THEN
          Jx(K,K+1) =      cmplx(funcP_JK(Jrot,K),kind=Rkind)
          Jy(K,K+1) = -EYE*cmplx(funcP_JK(Jrot,K),kind=Rkind)
        END IF

        IF (K > -Jrot) THEN
          Jx(K,K-1) =      cmplx(funcM_JK(Jrot,K),kind=Rkind)
          Jy(K,K-1) =  EYE*cmplx(funcM_JK(Jrot,K),kind=Rkind)
        END IF

      END DO

      !write(out_unitp,*)
      !write(out_unitp,*) 'Re Jx'
      !CALL Write_Mat(real(Jx),out_unitp,5)
      !write(out_unitp,*) 'Jy'
      !CALL Write_Mat(aimag(Jy),out_unitp,5)
      !write(out_unitp,*) 'Jz'
      !CALL Write_Mat(real(Jz),out_unitp,5)

!---------------------------------------------------------------------
      ! -3- Jx,Jy,Jz in the Wang basis
      Ji(:,:,-1) = matmul(PWangDag,matmul(Jx,PWang))
      Ji(:,:,-2) = matmul(PWangDag,matmul(Jy,PWang))
      Ji(:,:,-3) = matmul(PWangDag,matmul(Jz,PWang))

      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'Im(Jx), SumRe(Jx)',sum(abs(real(Ji(:,:,-1))))
        CALL Write_Mat(aimag(Ji(:,:,-1)),out_unitp,5)
        write(out_unitp,*) 'Im(Jy), SumRe(Jy)',sum(abs(real(Ji(:,:,-2))))
        CALL Write_Mat(aimag(Ji(:,:,-2)),out_unitp,5)
        write(out_unitp,*) 'Im(Jz), SumRe(Jz)',sum(abs(real(Ji(:,:,-3))))
        CALL Write_Mat(aimag(Ji(:,:,-3)),out_unitp,5)
      END IF
      !CALL alloc_NParray(test_inv,(/ RotBasis_Para%nb_Rot,RotBasis_Para%nb_Rot /),         &
      !                  'test_inv',name_sub)
      !test_inv = matmul(Jx,Jx)+matmul(Jy,Jy)+matmul(Jz,Jz)
      !write(out_unitp,*) 'J^2'
      !CALL Write_Mat(test_inv,out_unitp,5)
      !CALL dealloc_NParray(test_inv,'test_inv',name_sub)

!---------------------------------------------------------------------
      ! -4- Transfert Jx,Jy .... Jx*Jy ... in RotBasis_Para%tab_RotOp
      DO J1=-3,-1
      DO J2=-3,-1
         iterm = RotBasis_Para%tab_der_TO_iterm(J1,J2)
         RotBasis_Para%tab_RotOp(:,:,iterm) = real(matmul(Ji(:,:,J1),Ji(:,:,J2)))
      END DO
      END DO
      DO J1=-3,-1
         iterm = RotBasis_Para%tab_der_TO_iterm(J1,0)
         RotBasis_Para%tab_RotOp(:,:,iterm) = aimag(Ji(:,:,J1))
      END DO

      ! J^2
      DO J1=-3,-1
        iterm = RotBasis_Para%tab_der_TO_iterm(J1,J1)
        RotBasis_Para%tab_RotOp(:,:,0) = RotBasis_Para%tab_RotOp(:,:,0) + &
            RotBasis_Para%tab_RotOp(:,:,iterm)
      END DO


      test_inv=matmul(Ji(:,:,-1),Ji(:,:,-2))-matmul(Ji(:,:,-2),Ji(:,:,-1))
      !write(out_unitp,*) '[Jx,Jy]'
      !CALL Write_Mat(test_inv,out_unitp,5)
      test_inv = test_inv + EYE*Ji(:,:,-3)
      write(out_unitp,*) '[Jx,Jy] + i Jz = 0 ?',sum(abs(test_inv))
      !CALL Write_Mat(test_inv,out_unitp,5)

      test_inv=matmul(Ji(:,:,-3),Ji(:,:,-1))-matmul(Ji(:,:,-1),Ji(:,:,-3))
      !write(out_unitp,*) '[Jz,Jx]'
      !CALL Write_Mat(test_inv,out_unitp,5)
      test_inv = test_inv+EYE*Ji(:,:,-2)
      write(out_unitp,*) '[Jz,Jx] + i Jy = 0 ?',sum(abs(test_inv))
      !CALL Write_Mat(test_inv,out_unitp,5)

      test_inv=matmul(Ji(:,:,-2),Ji(:,:,-3))-matmul(Ji(:,:,-3),Ji(:,:,-2))
      !write(out_unitp,*) '[Jy,Jz]'
      !CALL Write_Mat(test_inv,out_unitp,5)
      test_inv = test_inv+EYE*Ji(:,:,-1)
      write(out_unitp,*) '[Jy,Jz] + i Jx = 0 ?',sum(abs(test_inv))
      !CALL Write_Mat(test_inv,out_unitp,5)



!---------------------------------------------------------------------
      CALL dealloc_NParray(PWangDag,'PWangDag',name_sub)
      CALL dealloc_NParray(PWang,'PWang',name_sub)
      CALL dealloc_NParray(Jx,'Jx',name_sub)
      CALL dealloc_NParray(Jy,'Jy',name_sub)
      CALL dealloc_NParray(Jz,'Jz',name_sub)
      CALL dealloc_NParray(Ji,'Ji',name_sub)


!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        CALL Write_RotBasis_Param(RotBasis_Para)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------

!STOP 'Init_RotBasis_Param'
!---------------------------------------------------------------------

      END SUBROUTINE Init_RotBasis_Param_old

      FUNCTION funcP_JK(J,K)
        real (kind=Rkind) :: funcP_JK
        integer, intent(in) :: J,K

        funcP_JK = HALF*sqrt(real(J*J+J-K*K-K,kind=Rkind))

      END  FUNCTION funcP_JK
      FUNCTION funcM_JK(J,K)
        real (kind=Rkind) :: funcM_JK
        integer, intent(in) :: J,K

        funcM_JK = HALF*sqrt(real(J*J+J-K*K+K,kind=Rkind))

      END  FUNCTION funcM_JK
      SUBROUTINE Write_RotBasis_Param(RotBasis_Para)

      TYPE (RotBasis_Param), intent(in) :: RotBasis_Para


      integer :: J1,J2,iterm
!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Write_RotBasis_Param'
      !logical,parameter :: debug=.FALSE.
      logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------

      write(out_unitp,*) 'RotBasis_Para%Jrot   ',RotBasis_Para%Jrot
      write(out_unitp,*) 'RotBasis_Para%nb_Rot ',RotBasis_Para%nb_Rot
      write(out_unitp,*) 'RotBasis_Para%nb_term',RotBasis_Para%nb_term

      IF (allocated(RotBasis_Para%tab_RotOp)) THEN
        write(out_unitp,*) 'J(-1) => -i*Jx, J(-2) => -i*Jy, J(-3) => -i*Jz'
        DO J1=-3,-1
        DO J2=-3,-1
          iterm = RotBasis_Para%tab_der_TO_iterm(J1,J2)
          write(out_unitp,*) 'JiJj op.',J1,J2,' iterm: ',iterm
          CALL Write_Mat(RotBasis_Para%tab_RotOp(:,:,iterm),out_unitp,5)
        END DO
        END DO
        DO J1=-3,-1
          iterm = RotBasis_Para%tab_der_TO_iterm(J1,0)
          write(out_unitp,*) 'Im(Ji op.)',J1,' iterm: ',iterm
          CALL Write_Mat(RotBasis_Para%tab_RotOp(:,:,iterm),out_unitp,5)
        END DO

        write(out_unitp,*) 'Id'
        CALL Write_Mat(RotBasis_Para%tab_RotOp(:,:,0),out_unitp,5)
      END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE Write_RotBasis_Param

      SUBROUTINE RotBasis_Param2TORotBasis_Param1(RotBasis_Para1,       &
                                                         RotBasis_Para2)

      TYPE (RotBasis_Param), intent(inout) :: RotBasis_Para1
      TYPE (RotBasis_Param), intent(in)    :: RotBasis_Para2

!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='RotBasis_Param2TORotBasis_Param1'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------

      CALL alloc_RotBasis_Param(RotBasis_Para1,RotBasis_Para2%Jrot)

      IF (allocated(RotBasis_Para2%tab_RotOp)) THEN
        RotBasis_Para1%tab_RotOp(:,:,:) = RotBasis_Para2%tab_RotOp(:,:,:)
      END IF

      IF (allocated(RotBasis_Para2%tab_iterm_TO_der)) THEN
        RotBasis_Para1%tab_iterm_TO_der(:,:) = RotBasis_Para2%tab_iterm_TO_der(:,:)
      END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        CALL Write_RotBasis_Param(RotBasis_Para2)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE RotBasis_Param2TORotBasis_Param1



      END MODULE mod_RotBasis_Param
