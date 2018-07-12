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
!      with contributions of Mamadou Ndong
!
!===========================================================================
!===========================================================================
      MODULE mod_LinearNMTransfo
      USE mod_system
      USE mod_dnSVM
      USE mod_constant
      USE mod_file
      USE mod_string
      IMPLICIT NONE

      !! @description: TODO
      !! @param: TODO
      TYPE Type_LinearTransfo
        real (kind=Rkind), pointer  :: mat(:,:)=>null()
        real (kind=Rkind), pointer  :: mat_inv(:,:)=>null()
        logical :: inv    =.FALSE.
        logical :: transp =.FALSE.

        logical :: check_LinearTransfo=.TRUE.

      END TYPE Type_LinearTransfo


      !!@description: TODO
      !!@param: TODO
      TYPE Type_NMTransfo
        logical                    :: hessian_old      = .TRUE.
        logical                    :: hessian_cart     = .TRUE.
        logical                    :: hessian_onthefly = .FALSE.
        TYPE (param_file)          :: file_hessian
        integer                    :: nb_NM            = 0      ! nb_act

        logical                    :: d0c_read         = .FALSE.
        logical                    :: hessian_read     = .FALSE.
        logical                    :: k_read           = .FALSE.
        integer                    :: nb_read          = 0
        real (kind=Rkind), pointer :: d0h(:,:)         =>null()
        real (kind=Rkind), pointer :: d0k(:,:)         =>null()

        ! to set up automaticaly the HObasis
        real (kind=Rkind), pointer :: Q0_HObasis(:)        =>null()
        real (kind=Rkind), pointer :: scaleQ_HObasis(:)    =>null()


        real (kind=Rkind), pointer :: d0c_inv(:,:)     =>null()
        real (kind=Rkind), pointer :: d0c(:,:)         =>null()
        real (kind=Rkind), pointer :: phase(:)         =>null() ! to change the sign of d0c
        real (kind=Rkind), pointer :: d0eh(:)          =>null()
        logical                    :: purify_hess      = .FALSE.! if .TRUE., we use an hessian purified (default : .FALSE.)
        logical                    :: eq_hess          = .FALSE.! if .TRUE., we use an hessian purified (default : .FALSE.)
        logical                    :: k_Half           = .FALSE.! if .TRUE., make a transfo, to get -1/2D2./dx2
        integer, pointer           :: Qact1_sym(:)     =>null() ! Qact1_sym(nb_NM)  : for the symmetrized H0
        integer, pointer           :: Qact1_eq(:,:)    =>null() ! Qact1_eq(nb_NM,nb_NM) : for the symmetrized H0
        integer                    :: nb_equi          = 0      ! number of set of equivalent variables
        integer, pointer           :: dim_equi(:)      =>null() ! dimension for each set of equivalence
        integer, pointer           :: tab_equi(:,:)    =>null() ! list of variables for each set of equivalence
      END TYPE Type_NMTransfo

      INTERFACE alloc_array
        ! for RPHTransfo
        MODULE PROCEDURE alloc_array_OF_NMTransfodim0
      END INTERFACE
      INTERFACE dealloc_array
        ! for RPHTransfo
        MODULE PROCEDURE dealloc_array_OF_NMTransfodim0
      END INTERFACE

      CONTAINS

!================================================================
!      Subroutines for the linear Transfo:
!       alloc_LinearTransfo
!       dealloc_LinearTransfo
!       Read_linearTransfo
!       Check_linearTransfo
!       calc_lineartransfo
!================================================================
      !!@description: ubroutines for the linear Transfo:
      !!       alloc_LinearTransfo
      !!       dealloc_LinearTransfo
      !!       Read_linearTransfo
      !!       Check_linearTransfo
      !!       calc_lineartransfo
      !!@param: TODO
      SUBROUTINE alloc_LinearTransfo(LinearTransfo,nb_Qin)

      TYPE (Type_LinearTransfo), intent(inout) :: LinearTransfo
      integer, intent(in) :: nb_Qin

      character (len=*), parameter :: name_sub='alloc_LinearTransfo'


      IF (associated(LinearTransfo%mat))  THEN
        CALL dealloc_array(LinearTransfo%mat,                           &
                          "LinearTransfo%mat",name_sub)
      END IF
      IF (associated(LinearTransfo%mat_inv))  THEN
        CALL dealloc_array(LinearTransfo%mat_inv,                       &
                          "LinearTransfo%mat_inv",name_sub)
      END IF

      CALL alloc_array(LinearTransfo%mat,(/nb_Qin,nb_Qin/),             &
                      "LinearTransfo%mat",name_sub)
      LinearTransfo%mat(:,:) = ZERO
      CALL alloc_array(LinearTransfo%mat_inv,(/nb_Qin,nb_Qin/),         &
                      "LinearTransfo%mat_inv",name_sub)
      LinearTransfo%mat_inv(:,:) = ZERO

      END SUBROUTINE alloc_LinearTransfo
      !-----------------------------------------------------------------------
      !!@description: TODO
      !!@param: TODO
      SUBROUTINE dealloc_LinearTransfo(LinearTransfo)

      TYPE (Type_LinearTransfo), intent(inout) :: LinearTransfo

      character (len=*), parameter :: name_sub='dealloc_LinearTransfo'

      IF (associated(LinearTransfo%mat))  THEN
        CALL dealloc_array(LinearTransfo%mat,                           &
                          "LinearTransfo%mat",name_sub)
      END IF
      IF (associated(LinearTransfo%mat_inv))  THEN
        CALL dealloc_array(LinearTransfo%mat_inv,                       &
                          "LinearTransfo%mat_inv",name_sub)
      END IF

      LinearTransfo%inv                 = .FALSE.
      LinearTransfo%transp              = .FALSE.

      LinearTransfo%check_LinearTransfo = .TRUE.

      END SUBROUTINE dealloc_linearTransfo

    SUBROUTINE alloc_array_OF_NMTransfodim0(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (Type_NMTransfo), pointer, intent(out) :: tab

      character (len=*), intent(in) :: name_var,name_sub

      integer, parameter :: ndim=0
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_NMTransfodim0'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


       IF (associated(tab))                                             &
             CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

       memory = 1
       allocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Type_NMTransfo')

      END SUBROUTINE alloc_array_OF_NMTransfodim0
      SUBROUTINE dealloc_array_OF_NMTransfodim0(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (Type_NMTransfo), pointer, intent(inout) :: tab
      character (len=*), intent(in) :: name_var,name_sub

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_NMTransfodim0'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN
       IF (.NOT. associated(tab))                                       &
             CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       memory = 1
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Type_NMTransfo')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_NMTransfodim0

      SUBROUTINE Read_linearTransfo(LinearTransfo,nb_Qin)

      TYPE (Type_LinearTransfo), intent(inout) :: LinearTransfo
      integer, intent(in) :: nb_Qin

      integer :: i,it,err,nbcol

      integer :: err_mem,memory
      !logical, parameter :: debug=.TRUE.
      logical, parameter :: debug=.FALSE.
      character (len=*), parameter :: name_sub='Read_LinearTransfo'

      CALL alloc_LinearTransfo(LinearTransfo,nb_Qin)


      read(in_unitp,*,IOSTAT=err)
      IF (err /= 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' "End of file", while reading an empty line.'
        write(out_unitp,*) ' Check your data !!'
        STOP
      END IF
      read(in_unitp,*,IOSTAT=err) nbcol
      IF (err /= 0) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' "End of file", while reading nbcol'
        write(out_unitp,*) ' Check your data !!'
        STOP
      END IF

      IF (print_level > 1) write(out_unitp,*)'nbcol=',nbcol

      IF (LinearTransfo%inv) THEN

        CALL Read_RMat(LinearTransfo%mat_inv,in_unitp,nbcol,err)
        IF (LinearTransfo%transp) THEN
           LinearTransfo%mat_inv = transpose(LinearTransfo%mat_inv)
        END IF
        IF (err /= 0) THEN
          write(out_unitp,*) 'ERROR ',name_sub
          write(out_unitp,*) ' while reading the matrix "LinearTransfo%mat_inv"'
          write(out_unitp,*) ' Check your data !!'
          STOP
        END IF
        write(out_unitp,*) 'mat_inv of LinearTransfo has been read'

        CALL inv_m1_TO_m2(LinearTransfo%mat_inv,LinearTransfo%mat,nb_Qin,1,ONETENTH**10) ! SVD

      ELSE

        CALL Read_RMat(LinearTransfo%mat,in_unitp,nbcol,err)
        IF (LinearTransfo%transp) THEN
           LinearTransfo%mat = transpose(LinearTransfo%mat)
        END IF
        IF (err /= 0) THEN
          write(out_unitp,*) 'ERROR ',name_sub
          write(out_unitp,*) ' while reading the matrix "LinearTransfo%mat"'
          write(out_unitp,*) ' Check your data !!'
          STOP
        END IF
        write(out_unitp,*) 'mat of LinearTransfo has been read'

        CALL inv_m1_TO_m2(LinearTransfo%mat,LinearTransfo%mat_inv,nb_Qin,1,ONETENTH**10) ! SVD

      END IF

      IF (print_level > 1) THEN
        write(out_unitp,*)  'mat of LinearTransfo: '
        CALL Write_Mat(LinearTransfo%mat,out_unitp,4)
        write(out_unitp,*)  'mat_inv of LinearTransfo: '
        CALL Write_Mat(LinearTransfo%mat_inv,out_unitp,4)
      END IF
      CALL flush_perso(out_unitp)

      END SUBROUTINE Read_LinearTransfo
!-----------------------------------------------------------------------

      SUBROUTINE Read_LC_projectionTransfo(LinearTransfo,               &
                                  nb_transfo,opt_transfo,not_all,nb_Qin)

      TYPE (Type_LinearTransfo), intent(inout) :: LinearTransfo
      integer, intent(in) :: nb_transfo,opt_transfo,nb_Qin
      logical, intent(in) :: not_all

      real (kind=Rkind) :: mat_inv(nb_Qin,nb_Qin)
      real (kind=Rkind) :: LC(nb_transfo,nb_Qin),LC_read(nb_transfo,nb_Qin)
      real (kind=Rkind) :: S_LC(nb_transfo,nb_transfo)
      real (kind=Rkind) :: Sinv_LC(nb_transfo,nb_transfo)
      real (kind=Rkind) :: V(nb_Qin)
      logical           :: Tab_Qi_in_Mat(nb_Qin)

      !real (kind=Rkind) :: S(nb_Qin,nb_Qin)
      !real (kind=Rkind) :: SS


      integer           :: iLC,i,ii,j,it,err,nbcol,k,ic,nb_coef
      real (kind=Rkind) :: norm_ii,norm_jj,norm_ij,coef

      integer :: err_mem,memory
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='Read_LC_projectionTransfo'

      CALL alloc_LinearTransfo(LinearTransfo,nb_Qin)


      LC_read(:,:) = ZERO
      IF (not_all) THEN
        DO iLC=1,nb_transfo
          read(in_unitp,*,IOSTAT=err) nb_coef
          DO i=1,nb_coef
            read(in_unitp,*,IOSTAT=err) ic,coef
            LC_read(iLC,ic) = coef
          END DO
          IF (err /= 0) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' "End of file" or "end of record", while reading the linear combination'
            write(out_unitp,*) ' Check your data !!'
            STOP
          END IF
          write(out_unitp,*) "iLC,norm",iLC,dot_product(LC_read(iLC,:),LC_read(iLC,:)),'Read vect:',LC_read(iLC,:)
        END DO
      ELSE
        DO iLC=1,nb_transfo
          read(in_unitp,*,IOSTAT=err) i,LC_read(iLC,:)
          IF (err /= 0) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' "End of file" or "end of record", while reading the linear combination'
            write(out_unitp,*) ' Check your data !!'
            STOP
          END IF
          write(out_unitp,*) "iLC,norm",iLC,dot_product(LC_read(iLC,:),LC_read(iLC,:)),'Read vect:',LC_read(iLC,:)
        END DO
      END IF
      LC(:,:) = LC_read(:,:)

      ! here the projection
      ! first Schmidt orthogonalization of the linear combinations
      DO i=1,nb_transfo
        DO j=1,i-1
          norm_jj = dot_product(LC(j,:),LC(j,:))
          norm_ij = dot_product(LC(i,:),LC(j,:))
          LC(i,:) = LC(i,:)*norm_jj - LC(j,:)*norm_ij
        END DO
        norm_ii = dot_product(LC(i,:),LC(i,:))
        IF (debug) write(out_unitp,*) "iLC,norm",i,norm_ii,'vect:',LC(i,:)

        IF (norm_ii < ONETENTH**5) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' Your linear combinations are not independent!!'
          write(out_unitp,*) ' CHECK your data'
          STOP
        END IF
        LC(i,:) = LC(i,:) / sqrt(norm_ii)
      END DO


      ! then the projection of the LC vectors on the mat_inv
      CALL mat_id(mat_inv,nb_Qin,nb_Qin) ! initializationwith the identity matrix

      DO i=1,nb_transfo
        norm_ii = dot_product(LC(i,:),LC(i,:))
        !write(out_unitp,*) 'LC i norm:',i,norm_ii
        DO j=1,nb_Qin
          norm_ij = dot_product(LC(i,:),mat_inv(j,:))
          !write(out_unitp,*) 'over i,j:',i,j,norm_ij
          mat_inv(j,:) = mat_inv(j,:)*norm_ii - LC(i,:)*norm_ij
        END DO
      END DO

      ! Finnally Schmidt orthogonalization of mat_inv
      DO i=1,nb_Qin
        DO j=1,i-1
          norm_jj = dot_product(mat_inv(j,:),mat_inv(j,:))
          IF (norm_jj < ONETENTH**5) CYCLE
          norm_ij = dot_product(mat_inv(i,:),mat_inv(j,:))
          mat_inv(i,:) = mat_inv(i,:)*norm_jj - mat_inv(j,:)*norm_ij
        END DO
        norm_ii = dot_product(mat_inv(i,:),mat_inv(i,:))
        IF (norm_ii > ONETENTH**5) mat_inv(i,:) = mat_inv(i,:) / sqrt(norm_ii)
      END DO

      ! Set LinearTransfo%mat_inv. Here, just from the LC_read (several options)
      SELECT CASE (opt_transfo)
      CASE (1,11)
        LinearTransfo%mat_inv(1:nb_transfo,:) = LC_read(1:nb_transfo,:)
      CASE (2,12)
        LinearTransfo%mat_inv(1:nb_transfo,:) = LC(1:nb_transfo,:)
      CASE (3,13)
        ! Overlap matrix of the LC_read and its invers
        S_LC(:,:) = matmul(LC_read,transpose(LC_read))
        CALL inv_m1_TO_m2(S_LC,Sinv_LC,nb_transfo,1,ONETENTH**10) ! SVD

        IF (debug) THEN
          write(out_unitp,*) 'S_LC'
          CALL Write_Mat(S_LC,out_unitp,4)
          write(out_unitp,*) 'Sinv_LC'
          CALL Write_Mat(Sinv_LC,out_unitp,4)
        END IF

        ! the contribution of mat_inv
        !LinearTransfo%mat_inv(1:nb_transfo,:) = matmul(Sinv_LC,LC_read)
        LC = matmul(Sinv_LC,LC_read)

      CASE DEFAULT
        ! Overlap matrix of the LC_read and its invers
        S_LC(:,:) = matmul(LC_read,transpose(LC_read))
        DO i=1,nb_transfo
        DO j=1,nb_transfo
          S_LC(i,j) = S_LC(i,j)/ sqrt(S_LC(i,i)*S_LC(j,j))
        END DO
        END DO
        CALL inv_m1_TO_m2(S_LC,Sinv_LC,nb_transfo,1,ONETENTH**10) ! SVD

        IF (debug) THEN
          write(out_unitp,*) 'S_LC'
          CALL Write_Mat(S_LC,out_unitp,4)

          write(out_unitp,*) 'Sinv_LC'
          CALL Write_Mat(Sinv_LC,out_unitp,4)
        END IF

        ! the contribution of mat_inv
        LC = matmul(Sinv_LC,LC_read)

      END SELECT

      IF (debug) THEN
        DO i=1,nb_transfo
          write(out_unitp,*) 'LC_read',i,LC_read(i,:)
          write(out_unitp,*) 'LC',i,LC(i,:)
        END DO
      END IF

      ! With this procedure a vector (LC(i,:) or mat_inv(i,:) is transfer close to
      !   the first none-zero coeficient of the vector (if opt_transfo > 10)
      ! otherwise, the vectors are transfer one after each orther.
      Tab_Qi_in_Mat(:) = .FALSE.
      DO i=1,nb_transfo
        k = kNoneZero(LC(i,:),Tab_Qi_in_Mat,opt_transfo)
        IF (debug) write(out_unitp,*) 'LC',i,k,LC(i,:)
        LinearTransfo%mat_inv(k,:) = LC(i,:)
      END DO
      DO i=1,nb_Qin
        IF (count(Tab_Qi_in_Mat) == nb_Qin) EXIT
        norm_ii = dot_product(mat_inv(i,:),mat_inv(i,:))
        IF (norm_ii < ONETENTH**5) CYCLE
        k = kNoneZero(mat_inv(i,:),Tab_Qi_in_Mat,opt_transfo)
        IF (debug) write(out_unitp,*) 'mat_inv',i,k,mat_inv(i,:)
        LinearTransfo%mat_inv(k,:) = mat_inv(i,:)
      END DO

!      DO i=1,nb_Qin
!      DO j=1,nb_Qin
!          S(i,j) = dot_product(LinearTransfo%mat_inv(i,:),LinearTransfo%mat_inv(j,:))
!          SS = dot_product(LinearTransfo%mat_inv(i,:),LinearTransfo%mat_inv(j,:))
!          IF (i == j .AND. abs(SS-ONE) > ONETENTH**5) THEN
!            write(out_unitp,*) 'Overlap of mat_inv,i,j',i,j,SS
!          END IF
!          IF (i /= j .AND. abs(SS) > ONETENTH**5) THEN
!            write(out_unitp,*) 'Overlap of mat_inv,i,j',i,j,SS
!          END IF
!      END DO
!      END DO
      !write(out_unitp,*)  'Overlap of mat_inv: '
      !CALL Write_Mat(S,out_unitp,4)
      write(out_unitp,*) 'mat_inv is set-up'


      CALL inv_m1_TO_m2(LinearTransfo%mat_inv,LinearTransfo%mat,nb_Qin,1,ONETENTH**10) ! SVD

      DO i=1,nb_transfo
         V(:) = matmul(LinearTransfo%mat_inv,LC_read(i,:))
         !write(6,*) 'mat_inv*LC_read,',i,V
         V(i) = V(i)-ONE
         write(out_unitp,*) 'Norm of mat_inv*LC_read-LC_read,',i,dot_product(V,V)
      END DO



      IF (.NOT. LinearTransfo%inv) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' inv=.TRUE. should not be possible'
        write(out_unitp,*) ' CHECK the fortran!!'
        STOP
      END IF


      IF (print_level > 1) THEN
        write(out_unitp,*)  'mat of LinearTransfo: '
        CALL Write_Mat(LinearTransfo%mat,out_unitp,4)
        write(out_unitp,*)  'mat_inv of LinearTransfo: '
        CALL Write_Mat(LinearTransfo%mat_inv,out_unitp,4)
      END IF
      CALL flush_perso(out_unitp)

      END SUBROUTINE Read_LC_projectionTransfo

      integer FUNCTION kNoneZero(V,Tab_Qi_in_Mat,opt)

      real (kind=Rkind), intent(in) :: V(:)
      logical, intent(inout)        :: Tab_Qi_in_Mat(:)
      integer, intent(in)           :: opt


      integer           :: k

      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='kNoneZero'

      IF (opt > 10) THEN

        kNoneZero = 0
        DO k=1,size(V)
          IF (V(k) /= ZERO .AND. .NOT. Tab_Qi_in_Mat(k)) THEN
            kNoneZero = k
            Tab_Qi_in_Mat(k) = .TRUE.
            EXIT
          END IF
        END DO
        !write(out_unitp,*) 'Vec',kNoneZero,V(:)
        IF (kNoneZero == 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' Problem with the projection !!!!!'
          write(out_unitp,*) ' CHECK the fortran'
        END IF
      ELSE
        k = count(Tab_Qi_in_Mat) + 1
        kNoneZero = k
        Tab_Qi_in_Mat(k) = .TRUE.
      END IF

      END FUNCTION kNoneZero

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE calc_LinearTransfo(dnQin,dnQout,LinearTransfo,nderiv,inTOout)

        TYPE (Type_dnVec), intent(inout)      :: dnQin,dnQout
        TYPE (Type_LinearTransfo), intent(in) :: LinearTransfo

        integer, intent(in)                   :: nderiv
        logical, intent(in)                   :: inTOout


        integer :: i,j,k
        character (len=*), parameter :: name_sub='calc_LinearTransfo'



        CALL check_alloc_dnVec(dnQin,'dnQin',name_sub)
        CALL check_alloc_dnVec(dnQout,'dnQout',name_sub)

        IF (inTOout) THEN
          IF (nderiv == 0) THEN
            dnQout%d0 = matmul(LinearTransfo%mat,dnQin%d0)
          ELSE IF (nderiv == 1) THEN
            dnQout%d0 = matmul(LinearTransfo%mat,dnQin%d0)
            DO i=1,dnQin%nb_var_deriv ! mole%nb_act
              dnQout%d1(:,i) = matmul(LinearTransfo%mat,dnQin%d1(:,i))
            END DO
          ELSE IF (nderiv == 2) THEN
            dnQout%d0 = matmul(LinearTransfo%mat,dnQin%d0)
            DO i=1,dnQin%nb_var_deriv ! mole%nb_act
              dnQout%d1(:,i) = matmul(LinearTransfo%mat,dnQin%d1(:,i))
            END DO
            DO i=1,dnQin%nb_var_deriv ! mole%nb_act
            DO j=1,dnQin%nb_var_deriv ! mole%nb_act
              dnQout%d2(:,i,j) = matmul(LinearTransfo%mat,dnQin%d2(:,i,j))
            END DO
            END DO
          ELSE IF (nderiv == 3) THEN
            dnQout%d0 = matmul(LinearTransfo%mat,dnQin%d0)
            DO i=1,dnQin%nb_var_deriv ! mole%nb_act
              dnQout%d1(:,i) = matmul(LinearTransfo%mat,dnQin%d1(:,i))
            END DO
            DO i=1,dnQin%nb_var_deriv ! mole%nb_act
            DO j=1,dnQin%nb_var_deriv ! mole%nb_act
              dnQout%d2(:,i,j) = matmul(LinearTransfo%mat,dnQin%d2(:,i,j))
            END DO
            END DO
            DO i=1,dnQin%nb_var_deriv ! mole%nb_act
            DO j=1,dnQin%nb_var_deriv ! mole%nb_act
            DO k=1,dnQin%nb_var_deriv ! mole%nb_act
              dnQout%d3(:,i,j,k) = matmul(LinearTransfo%mat,dnQin%d3(:,i,j,k))
            END DO
            END DO
            END DO
          ELSE
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' nderiv > 4 is NOT possible',nderiv
            write(out_unitp,*) 'It should never append! Check the source'
            STOP
          END IF
        ELSE
          IF (nderiv == 0) THEN
            dnQin%d0 = matmul(LinearTransfo%mat_inv,dnQout%d0)
          ELSE IF (nderiv == 1) THEN
            dnQin%d0 = matmul(LinearTransfo%mat_inv,dnQout%d0)
            DO i=1,dnQout%nb_var_deriv ! mole%nb_act
              dnQin%d1(:,i) = matmul(LinearTransfo%mat_inv,dnQout%d1(:,i))
            END DO
          ELSE IF (nderiv == 2) THEN
            dnQin%d0 = matmul(LinearTransfo%mat_inv,dnQout%d0)
            DO i=1,dnQout%nb_var_deriv ! mole%nb_act
              dnQin%d1(:,i) = matmul(LinearTransfo%mat_inv,dnQout%d1(:,i))
            END DO
            DO i=1,dnQout%nb_var_deriv ! mole%nb_act
            DO j=1,dnQout%nb_var_deriv ! mole%nb_act
              dnQin%d2(:,i,j) = matmul(LinearTransfo%mat_inv,dnQout%d2(:,i,j))
            END DO
            END DO
          ELSE IF (nderiv == 3) THEN
            dnQin%d0 = matmul(LinearTransfo%mat_inv,dnQout%d0)
            DO i=1,dnQout%nb_var_deriv ! mole%nb_act
              dnQin%d1(:,i) = matmul(LinearTransfo%mat_inv,dnQout%d1(:,i))
            END DO
            DO i=1,dnQout%nb_var_deriv ! mole%nb_act
            DO j=1,dnQout%nb_var_deriv ! mole%nb_act
              dnQin%d2(:,i,j) = matmul(LinearTransfo%mat_inv,dnQout%d2(:,i,j))
            END DO
            END DO
            DO i=1,dnQout%nb_var_deriv ! mole%nb_act
            DO j=1,dnQout%nb_var_deriv ! mole%nb_act
            DO k=1,dnQout%nb_var_deriv ! mole%nb_act
              dnQin%d3(:,i,j,k) = matmul(LinearTransfo%mat_inv,dnQout%d3(:,i,j,k))
            END DO
            END DO
            END DO
          ELSE
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' nderiv > 4 is NOT possible',nderiv
            write(out_unitp,*) 'It should never append! Check the source'
            STOP
          END IF
        END IF


      END SUBROUTINE calc_LinearTransfo
!================================================================
!      Subroutines for the NM (Normal Modes) Transfo:
!       alloc_NMTransfo
!       dealloc_NMTransfo
!       Read_NMTransfo
!       calc_NMTransfo
!================================================================
      !!@description: Subroutines for the NM (Normal Modes) Transfo:
      !!       alloc_NMTransfo
      !!       dealloc_NMTransfo
      !!       Read_NMTransfo
      !!       calc_NMTransfo
      !!@param: TODO
      SUBROUTINE dealloc_NMTransfo(NMTransfo)

      TYPE (Type_NMTransfo), intent(inout) :: NMTransfo

      character (len=*), parameter :: name_sub='dealloc_NMTransfo'

      IF (associated(NMTransfo%d0c_inv))  THEN
        CALL dealloc_array(NMTransfo%d0c_inv,"NMTransfo%d0c_inv",name_sub)
      END IF
      IF (associated(NMTransfo%d0c))  THEN
        CALL dealloc_array(NMTransfo%d0c,"NMTransfo%d0c",name_sub)
      END IF
      IF (associated(NMTransfo%d0eh))  THEN
        CALL dealloc_array(NMTransfo%d0eh,"NMTransfo%d0eh",name_sub)
      END IF
      IF (associated(NMTransfo%dim_equi))  THEN
        CALL dealloc_array(NMTransfo%dim_equi,"NMTransfo%dim_equi",name_sub)
      END IF
      IF (associated(NMTransfo%tab_equi))  THEN
        CALL dealloc_array(NMTransfo%tab_equi,"NMTransfo%tab_equi",name_sub)
      END IF
      IF (associated(NMTransfo%Qact1_sym))  THEN
        CALL dealloc_array(NMTransfo%Qact1_sym,"NMTransfo%Qact1_sym",name_sub)
      END IF
      IF (associated(NMTransfo%Qact1_eq))  THEN
        CALL dealloc_array(NMTransfo%Qact1_eq,"NMTransfo%Qact1_eq",name_sub)
      END IF

      NMTransfo%hessian_old      = .TRUE.
      NMTransfo%hessian_cart     = .TRUE.
      NMTransfo%hessian_onthefly = .FALSE.

      NMTransfo%file_hessian%name      = "file_hessian"
      NMTransfo%file_hessian%unit      = 0
      NMTransfo%file_hessian%formatted = .TRUE.
      NMTransfo%file_hessian%append    = .FALSE.
      NMTransfo%file_hessian%old       = NMTransfo%hessian_old

      NMTransfo%hessian_read     = .FALSE.
      NMTransfo%k_read           = .FALSE.
      NMTransfo%d0c_read         = .FALSE.

      IF (associated(NMTransfo%d0h)) THEN
        CALL dealloc_array(NMTransfo%d0h,"NMTransfo%d0h",name_sub)
      END IF
      IF (associated(NMTransfo%d0k)) THEN
        CALL dealloc_array(NMTransfo%d0k,"NMTransfo%d0k",name_sub)
      END IF

      IF (associated(NMTransfo%Q0_HObasis)) THEN
        CALL dealloc_array(NMTransfo%Q0_HObasis,"NMTransfo%Q0_HObasis",name_sub)
      END IF
      IF (associated(NMTransfo%scaleQ_HObasis)) THEN
        CALL dealloc_array(NMTransfo%scaleQ_HObasis,"NMTransfo%scaleQ_HObasis",name_sub)
      END IF

      NMTransfo%nb_NM       = 0

      NMTransfo%purify_hess = .FALSE.
      NMTransfo%eq_hess     = .FALSE.
      NMTransfo%k_Half      = .FALSE.
      NMTransfo%nb_equi     = 0

      END SUBROUTINE dealloc_NMTransfo

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE Read_NMTransfo(NMTransfo,nb_Qin)

      integer, intent(in) :: nb_Qin
      TYPE (Type_NMTransfo), intent(inout) :: NMTransfo

      integer                  :: i,k,it,nb_col,nb_NM
      character (len=Name_len) :: name0
      real(kind=Rkind), allocatable :: mat(:,:)

      logical, parameter :: debug=.TRUE.
      !logical, parameter :: debug=.FALSE.
      integer            :: err_read
      character (len=*), parameter :: name_sub='Read_NMTransfo'

      IF (NMTransfo%d0c_read) THEN

        read(in_unitp,*,IOSTAT=err_read)     ! for read a title (like d0c)
        IF (err_read /= 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' "End of file", while reading an empty ', &
             ' line or title of d0c matrix.'
          write(out_unitp,*) ' NMTransfo%d0c_read: ',NMTransfo%d0c_read
          write(out_unitp,*) ' Check your data !!'
          STOP
        END IF
        read(in_unitp,*,IOSTAT=err_read) nb_NM,nb_col
        IF (err_read /= 0) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' "End of file", while reading nb_NM,nb_col',&
                          ' d0c matrix.'
          write(out_unitp,*) ' Check your data !!'
          STOP
        END IF

        NMTransfo%nb_NM = nb_NM
        IF (.NOT. associated(NMTransfo%d0c)) THEN
          CALL alloc_array(NMTransfo%d0c,(/nb_NM,nb_NM/),"NMTransfo%d0c",name_sub)
          NMTransfo%d0c(:,:) = ZERO
        END IF

        CALL Read_RMat(NMTransfo%d0c(:,:),in_unitp,nb_col,err_read)
        IF (err_read /= 0) THEN
          write(out_unitp,*) 'ERROR ',name_sub
          write(out_unitp,*) ' while reading the matrix "NMTransfo%d0c"'
          write(out_unitp,*) ' Check your data !!'
          STOP
        END IF
        IF (debug) CALL Write_RMat(NMTransfo%d0c(:,:),out_unitp,nb_col,Rformat='f10.6')

        IF (debug) THEN
          IF (allocated(mat)) CALL dealloc_NParray(mat,"mat",name_sub)
          CALL alloc_NParray(mat,(/nb_NM,nb_NM/),"mat",name_sub)
          mat = matmul(transpose(NMTransfo%d0c),NMTransfo%d0c)
          write(out_unitp,*) ' td0c.d0c'
          CALL Write_RMat(mat,out_unitp,nb_col,Rformat='f10.6')
          CALL dealloc_NParray(mat,"mat",name_sub)
        END IF

      END IF
      CALL flush_perso(out_unitp)


      IF (NMTransfo%hessian_read .NEQV. NMTransfo%k_read) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' You MUST read both hessian and k matrix'
        write(out_unitp,*) ' hessian_read: ',NMTransfo%hessian_read
        write(out_unitp,*) ' k_read:       ',NMTransfo%k_read
        write(out_unitp,*) ' Check your data !!'
        STOP
      END IF

      IF (NMTransfo%hessian_read .AND. NMTransfo%nb_read < 1) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' You WANT to read both hessian and k matrix'
        write(out_unitp,*) ' but nb_read is < 1: ',NMTransfo%nb_read
        write(out_unitp,*) ' Check your data !!'
        STOP
      END IF
      IF (NMTransfo%hessian_read) THEN

        DO i=1,NMTransfo%nb_read

          read(in_unitp,*,IOSTAT=err_read)     ! for read a title (like d0h)
          IF (err_read /= 0) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' "End of file", while reading an empty ', &
               ' line or title of d0h matrices.'
            write(out_unitp,*) ' i_read,NMTransfo%nb_read: ',i,NMTransfo%nb_read
            write(out_unitp,*) ' NMTransfo%hessian_read: ',NMTransfo%hessian_read
            write(out_unitp,*) ' Check your data !!'
            STOP
          END IF
          read(in_unitp,*,IOSTAT=err_read) nb_NM,nb_col
          IF (err_read /= 0) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' "End of file", while reading nb_NM,nb_col ',&
                            ' of d0h matrices.'
            write(out_unitp,*) ' i_read,NMTransfo%nb_read: ',i,NMTransfo%nb_read
            write(out_unitp,*) ' NMTransfo%hessian_read: ',NMTransfo%hessian_read
            write(out_unitp,*) ' Check your data !!'
            STOP
          END IF

          NMTransfo%nb_NM = nb_NM
          IF (.NOT. associated(NMTransfo%d0h)) THEN
            CALL alloc_array(NMTransfo%d0h,(/nb_NM,nb_NM/),"NMTransfo%d0h",name_sub)
            NMTransfo%d0h(:,:) = ZERO
            CALL alloc_NParray(mat,(/nb_NM,nb_NM/),"mat",name_sub)
          END IF


          CALL Read_RMat(mat,in_unitp,nb_col,err_read)
          IF (err_read /= 0) THEN
            write(out_unitp,*) 'ERROR ',name_sub
            write(out_unitp,*) ' reading the matrix "NMTransfo%d0h"'
            write(out_unitp,*) ' i_read,NMTransfo%nb_read: ',i,NMTransfo%nb_read
            write(out_unitp,*) ' NMTransfo%hessian_read: ',NMTransfo%hessian_read
            write(out_unitp,*) ' Check your data !!'
            STOP
          END IF
          IF (debug) CALL Write_RMat(mat,out_unitp,nb_col)

          NMTransfo%d0h(:,:) = NMTransfo%d0h(:,:) + mat(:,:)

        END DO
        NMTransfo%d0h(:,:) = NMTransfo%d0h(:,:)/real(NMTransfo%nb_read,kind=Rkind)

        IF (allocated(mat)) CALL dealloc_NParray(mat,"mat",name_sub)

      END IF
      CALL flush_perso(out_unitp)

      IF (NMTransfo%k_read) THEN

        DO i=1,NMTransfo%nb_read

          read(in_unitp,*,IOSTAT=err_read)     ! for read a title (like d0h)
          IF (err_read /= 0) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' "End of file", while reading an empty ', &
               ' line or title of d0k matrices.'
            write(out_unitp,*) ' i_read,NMTransfo%nb_read: ',i,NMTransfo%nb_read
            write(out_unitp,*) ' NMTransfo%k_read: ',NMTransfo%k_read
            write(out_unitp,*) ' Check your data !!'
            STOP
          END IF
          read(in_unitp,*,IOSTAT=err_read) nb_NM,nb_col
          IF (err_read /= 0) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' "End of file", while reading nb_NM,nb_col',&
                            ' of d0k matrices.'
            write(out_unitp,*) ' i_read,NMTransfo%nb_read: ',i,NMTransfo%nb_read
            write(out_unitp,*) ' NMTransfo%k_read: ',NMTransfo%k_read
            write(out_unitp,*) ' Check your data !!'
            STOP
          END IF

          IF (.NOT. associated(NMTransfo%d0k)) THEN
            CALL alloc_array(NMTransfo%d0k,(/nb_NM,nb_NM/),"NMTransfo%d0k",name_sub)
            NMTransfo%d0k(:,:) = ZERO
            CALL alloc_NParray(mat,(/nb_NM,nb_NM/),"mat",name_sub)
          END IF


          CALL Read_RMat(mat,in_unitp,nb_col,err_read)
          IF (err_read /= 0) THEN
            write(out_unitp,*) 'ERROR ',name_sub
            write(out_unitp,*) ' reading the matrix "NMTransfo%d0k"'
            write(out_unitp,*) ' i_read,NMTransfo%nb_read: ',i,NMTransfo%nb_read
            write(out_unitp,*) ' NMTransfo%k_read: ',NMTransfo%k_read
            write(out_unitp,*) ' Check your data !!'
            STOP
          END IF
          IF (debug) CALL Write_RMat(mat,out_unitp,nb_col)

          NMTransfo%d0k(:,:) = NMTransfo%d0k(:,:) + mat(:,:)

        END DO
        NMTransfo%d0k(:,:) = NMTransfo%d0k(:,:)/real(NMTransfo%nb_read,kind=Rkind)

        IF (allocated(mat)) CALL dealloc_NParray(mat,"mat",name_sub)

      END IF
      CALL flush_perso(out_unitp)

      IF (NMTransfo%purify_hess) THEN

        IF (.NOT. associated(NMTransfo%Qact1_sym)) THEN
          CALL alloc_array(NMTransfo%Qact1_sym,(/nb_Qin/),              &
                          "NMTransfo%Qact1_sym",name_sub)
        END IF
        IF (.NOT. associated(NMTransfo%Qact1_eq)) THEN
          CALL alloc_array(NMTransfo%Qact1_eq,(/nb_Qin,nb_Qin/),        &
                          "NMTransfo%Qact1_eq",name_sub)
        END IF
        IF (.NOT. associated(NMTransfo%dim_equi)) THEN
          CALL alloc_array(NMTransfo%dim_equi,(/nb_Qin/),               &
                          "NMTransfo%dim_equi",name_sub)
        END IF
        IF (.NOT. associated(NMTransfo%tab_equi)) THEN
          CALL alloc_array(NMTransfo%tab_equi,(/nb_Qin,nb_Qin/),        &
                          "NMTransfo%tab_equi",name_sub)
        END IF

        write(out_unitp,*)
        write(out_unitp,*) "========================================"
        write(out_unitp,*) 'Hessian purification',NMTransfo%purify_hess

        write(out_unitp,*) 'nb_Qin',nb_Qin
        NMTransfo%Qact1_sym(:)  = 0
        NMTransfo%Qact1_eq(:,:) = 0

        read(in_unitp,*,IOSTAT=err_read) name0,NMTransfo%Qact1_sym(:)
        IF (err_read /= 0) THEN
          write(out_unitp,*) 'ERROR ',name_sub
          write(out_unitp,*) ' while reading Qact1_sym',name0,NMTransfo%Qact1_sym(:)
          write(out_unitp,*) ' Check your data !!'
          STOP
        END IF

        IF (NMTransfo%eq_hess) THEN
          DO i=1,nb_Qin
            read(in_unitp,*,IOSTAT=err_read) name0,NMTransfo%Qact1_eq(i,:)
            IF (err_read /=0) THEN
              write(out_unitp,*) ' while reading Qact1_eq',name0,NMTransfo%Qact1_eq(i,:)
              EXIT
            END IF
          END DO
          IF (err_read /= 0) THEN
            write(out_unitp,*) 'WARNING ',name_sub
            write(out_unitp,*) ' while reading Qact1_eq'
            write(out_unitp,*) ' The matrix, Qact1_eq, is assumed to be zero'
            NMTransfo%Qact1_eq(:,:) = 0
          END IF
        END IF

        NMTransfo%tab_equi(:,:) = 0
        DO i=1,nb_Qin
          NMTransfo%dim_equi(i) = 1
          NMTransfo%tab_equi(i,NMTransfo%dim_equi(i)) = i
          DO k=1,nb_Qin
            IF (NMTransfo%Qact1_eq(i,k) == 1) THEN
              NMTransfo%dim_equi(i) = NMTransfo%dim_equi(i) + 1
              NMTransfo%tab_equi(i,NMTransfo%dim_equi(i)) = k
            END IF
          END DO
        END DO


        write(out_unitp,*) 'Hessian purification parameters'
        write(out_unitp,*) 'Qact1_sym',NMTransfo%Qact1_sym(:)
        DO i=1,nb_Qin
          write(out_unitp,*) 'Qact1_eq',i,NMTransfo%Qact1_eq(i,:)
        END DO
        write(out_unitp,*) 'dim_equi',NMTransfo%dim_equi(:)

        DO i=1,nb_Qin
          write(out_unitp,*) 'tab_equi:',i,NMTransfo%dim_equi(i),':',     &
                       NMTransfo%tab_equi(i,NMTransfo%dim_equi(i))
        END DO
        write(out_unitp,*) 'END Hessian purification'
        write(out_unitp,*) "========================================"
        write(out_unitp,*)
      END IF
      END SUBROUTINE Read_NMTransfo

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE Write_NMTransfo(NMTransfo)

      TYPE (Type_NMTransfo), intent(inout) :: NMTransfo

      integer           :: i,nb_NM,nb_Qin
      character (len=Name_len) :: name0

      character (len=*), parameter :: name_sub='Write_NMTransfo'

      write(out_unitp,*) 'BEGINNING ',name_sub

      write(out_unitp,*) 'hessian_old,hessian_cart,hessian_onthefly',   &
                          NMTransfo%hessian_old,NMTransfo%hessian_cart, &
                          NMTransfo%hessian_onthefly
      IF (NMTransfo%hessian_old) write(out_unitp,*) 'file_hessian',     &
                                            NMTransfo%file_hessian%name
      write(out_unitp,*) 'k_Half',NMTransfo%k_Half



      IF (associated(NMTransfo%Q0_HObasis)) THEN
        write(out_unitp,*) 'Q0_HObasis'
        CALL Write_VecMat(NMTransfo%Q0_HObasis,out_unitp,5)
      END IF
      IF (associated(NMTransfo%scaleQ_HObasis)) THEN
        write(out_unitp,*) 'scaleQ_HObasis'
        CALL Write_VecMat(NMTransfo%scaleQ_HObasis,out_unitp,5)
      END IF


      write(out_unitp,*) 'hessian_read,k_read',NMTransfo%hessian_read,NMTransfo%k_read
      write(out_unitp,*) 'nb_read',NMTransfo%nb_read

      write(out_unitp,*) 'd0c_read',NMTransfo%d0c_read


      IF (associated(NMTransfo%d0h)) THEN
        write(out_unitp,*) 'd0h'
        CALL Write_Mat(NMTransfo%d0h,out_unitp,5)
      END IF
      IF (associated(NMTransfo%d0k)) THEN
        write(out_unitp,*) 'd0k'
        CALL Write_Mat(NMTransfo%d0k,out_unitp,5)
      END IF


      nb_NM = NMTransfo%nb_NM
      write(out_unitp,*) 'nb_NM',nb_NM
      CALL flush_perso(out_unitp)

      IF (nb_NM > 0) THEN
        IF (associated(NMTransfo%d0c_inv)) THEN
          write(out_unitp,*)  'd0c_inv: '
          CALL Write_Mat(NMTransfo%d0c_inv,out_unitp,4)
        END IF
        CALL flush_perso(out_unitp)

        IF (associated(NMTransfo%d0c)) THEN
          write(out_unitp,*)  'd0c: '
          CALL Write_Mat(NMTransfo%d0c,out_unitp,4)
        END IF
        CALL flush_perso(out_unitp)

        IF (associated(NMTransfo%d0eh)) THEN
          write(out_unitp,*)  'd0eh: ',NMTransfo%d0eh(:)
        END IF
        CALL flush_perso(out_unitp)

      END IF

      IF (NMTransfo%purify_hess) THEN

        write(out_unitp,*)
        write(out_unitp,*) "========================================"
        write(out_unitp,*) 'Hessian purification',NMTransfo%purify_hess
        nb_Qin = size(NMTransfo%Qact1_sym(:))
        write(out_unitp,*)  'Qact1_sym: ',NMTransfo%Qact1_sym(:)

        IF (NMTransfo%eq_hess) THEN
          DO i=1,nb_Qin
            write(out_unitp,*) 'Qact1_eq',i,NMTransfo%Qact1_eq(i,:)
          END DO
        END IF

        write(out_unitp,*) 'dim_equi',NMTransfo%dim_equi(:)
        write(out_unitp,*)  'tab_equi: '
        DO i=1,nb_Qin
          write(out_unitp,*) 'tab_equi(i,:)',i,':',                     &
                           NMTransfo%tab_equi(i,1:NMTransfo%dim_equi(i))
        END DO

        write(out_unitp,*) 'END Hessian purification'
        write(out_unitp,*) "========================================"
        write(out_unitp,*)

      write(out_unitp,*) 'END ',name_sub

      END IF
      CALL flush_perso(out_unitp)
      END SUBROUTINE Write_NMTransfo

      !!@description: TODO
      !!@param: TODO
      SUBROUTINE NMTransfo1TONMTransfo2(NMTransfo1,NMTransfo2)

      TYPE (Type_NMTransfo), intent(in)  :: NMTransfo1
      TYPE (Type_NMTransfo), intent(inout) :: NMTransfo2
      integer :: nb_Qin,n1
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='NMTransfo1TONMTransfo2'

      CALL dealloc_NMTransfo(NMTransfo2)

      NMTransfo2%hessian_old      = NMTransfo1%hessian_old
      NMTransfo2%hessian_cart     = NMTransfo1%hessian_cart
      NMTransfo2%hessian_onthefly = NMTransfo1%hessian_onthefly
      NMTransfo2%file_hessian     = NMTransfo1%file_hessian

      NMTransfo2%hessian_read     = NMTransfo1%hessian_read
      NMTransfo2%k_read           = NMTransfo1%k_read
      NMTransfo2%nb_read          = NMTransfo1%nb_read
      NMTransfo2%d0c_read         = NMTransfo1%d0c_read


      IF (associated(NMTransfo1%d0h)) THEN
        CALL alloc_array(NMTransfo2%d0h,shape(NMTransfo1%d0h),          &
                        "NMTransfo2%d0h",name_sub)
        NMTransfo2%d0h(:,:) = NMTransfo1%d0h(:,:)
      END IF

      IF (associated(NMTransfo1%d0k)) THEN
        CALL alloc_array(NMTransfo2%d0k,shape(NMTransfo1%d0k),          &
                        "NMTransfo2%d0k",name_sub)
        NMTransfo2%d0k(:,:) = NMTransfo1%d0k(:,:)
      END IF


      IF (associated(NMTransfo1%Q0_HObasis)) THEN
        CALL alloc_array(NMTransfo2%Q0_HObasis,shape(NMTransfo1%Q0_HObasis),&
                        "NMTransfo2%Q0_HObasis",name_sub)
        NMTransfo2%Q0_HObasis(:) = NMTransfo1%Q0_HObasis(:)
      END IF
      IF (associated(NMTransfo1%scaleQ_HObasis)) THEN
        CALL alloc_array(NMTransfo2%scaleQ_HObasis,shape(NMTransfo1%scaleQ_HObasis),&
                        "NMTransfo2%scaleQ_HObasis",name_sub)
        NMTransfo2%scaleQ_HObasis(:) = NMTransfo1%scaleQ_HObasis(:)
      END IF



      IF (NMTransfo2%nb_NM > 0) THEN
        NMTransfo2%nb_NM            = NMTransfo1%nb_NM

        CALL alloc_array(NMTransfo2%d0c_inv,shape(NMTransfo1%d0c_inv),  &
                        "NMTransfo2%d0c_inv",name_sub)
        NMTransfo2%d0c_inv = NMTransfo1%d0c_inv

        CALL alloc_array(NMTransfo2%d0c,shape(NMTransfo1%d0c),          &
                        "NMTransfo2%d0c",name_sub)
        NMTransfo2%d0c = NMTransfo1%d0c

        CALL alloc_array(NMTransfo2%d0eh,shape(NMTransfo1%d0eh),        &
                        "NMTransfo2%d0eh",name_sub)
        NMTransfo2%d0eh = NMTransfo1%d0eh

      END IF

      NMTransfo2%k_Half      = NMTransfo1%k_Half
      NMTransfo2%purify_hess = NMTransfo1%purify_hess
      NMTransfo2%eq_hess     = NMTransfo1%eq_hess

      IF (NMTransfo2%purify_hess) THEN
        nb_Qin = size(NMTransfo1%Qact1_sym)

        CALL alloc_array(NMTransfo2%Qact1_sym,shape(NMTransfo1%Qact1_sym),&
                        "NMTransfo2%Qact1_sym",name_sub)
        NMTransfo2%Qact1_sym = NMTransfo1%Qact1_sym

        CALL alloc_array(NMTransfo2%Qact1_eq,shape(NMTransfo1%Qact1_eq),&
                        "NMTransfo2%Qact1_eq",name_sub)
        NMTransfo2%Qact1_eq  = NMTransfo1%Qact1_eq

        NMTransfo2%nb_equi = NMTransfo1%nb_equi

        CALL alloc_array(NMTransfo2%dim_equi,shape(NMTransfo1%dim_equi),&
                        "NMTransfo2%dim_equi",name_sub)
        NMTransfo2%dim_equi = NMTransfo1%dim_equi

        CALL alloc_array(NMTransfo2%tab_equi,shape(NMTransfo1%tab_equi),&
                        "NMTransfo2%tab_equi",name_sub)
        NMTransfo2%tab_equi = NMTransfo1%tab_equi
      END IF

      END SUBROUTINE NMTransfo1TONMTransfo2


      END MODULE mod_LinearNMTransfo

