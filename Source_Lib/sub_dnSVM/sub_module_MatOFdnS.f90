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
      MODULE mod_MatOFdnS
      use mod_system, only: out_unitp, write_error_not_null, sub_test_tab_ub, &
                            sub_test_tab_lb, error_memo_allo, write_error_null, &
                            rkind, alloc_array, write_vecmat, dealloc_array, &
                            flush_perso, one, zero, ten, onetenth, half, pi, two
      use mod_dnS, only: type_dns, alloc_dns, dealloc_dns, check_alloc_dns, &
                         write_dns, alloc_array, dealloc_array, sub_dns1_to_dns2, &
                         sub_dns1_prod_dns2_to_dns3, sub_dns1_wplus_dns2_to_dns3, &
                         sub_dns1_prod_w_to_dns2, sub_zero_to_dns, sub_dns1_to_dntr2, &
                         sub_dns1_plus_dns2_to_dns3, sub_weight_dns
      use mod_VecOFdnS, only: alloc_array, dealloc_array, normalization_of_vecofdns, &
                              alloc_vecofdns, vec1ofdns_dotproduct_vec2ofdns_to_dns3, &
                              dealloc_vecofdns, check_alloc_vecofdns
      IMPLICIT NONE

      PRIVATE

      INTERFACE alloc_array
        MODULE PROCEDURE alloc_array_OF_dnSdim2
      END INTERFACE
      INTERFACE dealloc_array
        MODULE PROCEDURE dealloc_array_OF_dnSdim2
      END INTERFACE

      PUBLIC :: alloc_array, dealloc_array
      PUBLIC :: alloc_MatOFdnS, dealloc_MatOFdnS, check_alloc_MatOFdnS, Write_MatOFdnS
      PUBLIC :: sub_Mat1OFdnS_TO_Mat2OFdnS
      PUBLIC :: DET_Mat3x3OFdnS_TO_dnS
      PUBLIC :: TRANS_Mat1OFdnS_TO_Mat2OFdnS
      PUBLIC :: DIAG_MatOFdnS
      PUBLIC :: Mat1OFdnS_MUL_Mat2OFdnS_TO_Mat3OFdnS

      PUBLIC :: Mat1OFdnS_wPLUS_Mat2OFdnS_TO_Mat3OFdnS
      PUBLIC :: Mat1OFdnS_MUL_Vec2OFdnS_TO_Vec3OFdnS
      PUBLIC :: MatOFdnS_TO_VecOFdnS
      PUBLIC :: sub_ZERO_TO_MatOFdnS
      PUBLIC :: sub_Weight_MatOFdnS
      PUBLIC :: sub_Id_TO_MatOFdnS

      CONTAINS
!
!================================================================
!
!     allocation
!
!================================================================

      SUBROUTINE alloc_MatOFdnS(MatOFdnS,nb_var_deriv,nderiv)

        TYPE (Type_dnS) :: MatOFdnS(:,:)
        integer, optional :: nb_var_deriv,nderiv

        integer :: i,j

        IF (size(MatOFdnS) < 1) THEN
          write(out_unitp,*) 'ERROR in alloc_MatOFdnS'
          write(out_unitp,*) 'the matrix MatOFdnS(:,:) is not allocated'
          write(out_unitp,*) ' each element cannot be allocated!'
          write(out_unitp,*) ' Check the fortran!'
          STOP
        END IF
        IF (present(nderiv)) MatOFdnS(:,:)%nderiv = nderiv
        IF (present(nb_var_deriv)) MatOFdnS(:,:)%nb_var_deriv = nb_var_deriv

        IF (MatOFdnS(1,1)%nb_var_deriv == 0) MatOFdnS(:,:)%nderiv = 0

        DO i=lbound(MatOFdnS,dim=1),ubound(MatOFdnS,dim=1)
        DO j=lbound(MatOFdnS,dim=2),ubound(MatOFdnS,dim=2)
          CALL alloc_dnS(MatOFdnS(i,j))
        END DO
        END DO

      END SUBROUTINE alloc_MatOFdnS
      SUBROUTINE dealloc_MatOFdnS(MatOFdnS)

        TYPE (Type_dnS) :: MatOFdnS(:,:)

        integer :: i,j

        DO i=lbound(MatOFdnS,dim=1),ubound(MatOFdnS,dim=1)
        DO j=lbound(MatOFdnS,dim=2),ubound(MatOFdnS,dim=2)
          CALL dealloc_dnS(MatOFdnS(i,j))
        END DO
        END DO

      END SUBROUTINE dealloc_MatOFdnS

      SUBROUTINE alloc_array_OF_dnSdim2(tab,tab_ub,name_var,name_sub,tab_lb)
      IMPLICIT NONE

      TYPE (Type_dnS), pointer, intent(inout) :: tab(:,:)
      integer, intent(in) :: tab_ub(:)
      integer, intent(in), optional :: tab_lb(:)

      character (len=*), intent(in) :: name_var,name_sub


      integer, parameter :: ndim=2
      logical :: memory_test

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_dnSdim2'
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
         allocate(tab(tab_lb(1):tab_ub(1),                              &
                      tab_lb(2):tab_ub(2)),stat=err_mem)
       ELSE
         memory = product(tab_ub(:))
         allocate(tab(tab_ub(1),tab_ub(2)),stat=err_mem)
       END IF
       CALL error_memo_allo(err_mem,memory,name_var,name_sub,'Type_dnS')

      END SUBROUTINE alloc_array_OF_dnSdim2
      SUBROUTINE dealloc_array_OF_dnSdim2(tab,name_var,name_sub)
      IMPLICIT NONE

      TYPE (Type_dnS), pointer, intent(inout) :: tab(:,:)
      character (len=*), intent(in) :: name_var,name_sub

      integer :: i1,i2
!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_dnSdim2'
      integer :: err_mem,memory
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------

       !IF (.NOT. associated(tab)) RETURN

       IF (.NOT. associated(tab))                                       &
                 CALL Write_error_null(name_sub_alloc,name_var,name_sub)

       DO i1=ubound(tab,dim=1),lbound(tab,dim=1)
       DO i2=ubound(tab,dim=2),lbound(tab,dim=2)
         CALL dealloc_dnS(tab(i1,i2))
       END DO
       END DO


       memory = size(tab)
       deallocate(tab,stat=err_mem)
       CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'Type_dnS')
       nullify(tab)

      END SUBROUTINE dealloc_array_OF_dnSdim2


!================================================================
!
!     check if alloc has been done
!
!================================================================

      SUBROUTINE check_alloc_MatOFdnS(A,name_A,name_sub)
        TYPE (Type_dnS), intent(in) :: A(:,:)
        character (len=*), intent(in) :: name_A
        character (len=*), intent(in) :: name_sub

        integer :: i,j

        DO i=lbound(A,dim=1),ubound(A,dim=1)
        DO j=lbound(A,dim=2),ubound(A,dim=2)
          CALL check_alloc_dnS(A(i,j),name_A,name_sub)
        END DO
        END DO

      END SUBROUTINE check_alloc_MatOFdnS

!================================================================
!        write the derived type
!================================================================
      SUBROUTINE Write_MatOFdnS(MatOFdnS,nderiv)
        TYPE (Type_dnS) :: MatOFdnS(:,:)
        integer, optional :: nderiv
        integer :: i,j,id,jd,kd,nb_var_deriv,nderiv_loc
        real (kind=Rkind), pointer :: mat(:,:)
        logical :: old = .FALSE.
        !logical :: old = .TRUE.

        CALL check_alloc_MatOFdnS(MatOFdnS,'MatOFdnS','Write_MatOFdnS')

        nderiv_loc = minval(MatOFdnS%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)


        write(out_unitp,*) 'MatOFdnS'
        nb_var_deriv = MatOFdnS(lbound(MatOFdnS,dim=1),lbound(MatOFdnS,dim=2))%nb_var_deriv
        IF (old) THEN
          DO i=lbound(MatOFdnS,dim=1),ubound(MatOFdnS,dim=1)
          DO j=lbound(MatOFdnS,dim=2),ubound(MatOFdnS,dim=2)
            write(out_unitp,*) 'MatOFdnS',i,j
            CALL Write_dnS(MatOFdnS(i,j),nderiv_loc)
          END DO
          END DO
        ELSE
          nullify(mat)
          CALL alloc_array(mat,ubound(MatOFdnS),'mat','Write_MatOFdnS',lbound(MatOFdnS))

          IF (nderiv_loc >= 0) THEN
            mat(:,:) = MatOFdnS(:,:)%d0
            write(out_unitp,*) 'd0'
            CALL Write_VecMat(mat,out_unitp,5)
          END IF
          IF (nderiv_loc >= 1) THEN
            DO id=1,nb_var_deriv
              DO i=lbound(MatOFdnS,dim=1),ubound(MatOFdnS,dim=1)
              DO j=lbound(MatOFdnS,dim=2),ubound(MatOFdnS,dim=2)
                 mat(i,j) = MatOFdnS(i,j)%d1(id)
              END DO
              END DO
              write(out_unitp,*) 'd1',id
              CALL Write_VecMat(mat,out_unitp,5)
            END DO
          END IF
          IF (nderiv_loc >= 2) THEN
            DO id=1,nb_var_deriv
            DO jd=1,nb_var_deriv
              DO i=lbound(MatOFdnS,dim=1),ubound(MatOFdnS,dim=1)
              DO j=lbound(MatOFdnS,dim=2),ubound(MatOFdnS,dim=2)
                 mat(i,j) = MatOFdnS(i,j)%d2(id,jd)
              END DO
              END DO
              write(out_unitp,*) 'd2',id,jd
              CALL Write_VecMat(mat,out_unitp,5)
            END DO
            END DO
          END IF
          IF (nderiv_loc >= 3) THEN
            DO id=1,nb_var_deriv
            DO jd=1,nb_var_deriv
            DO kd=1,nb_var_deriv
              DO i=lbound(MatOFdnS,dim=1),ubound(MatOFdnS,dim=1)
              DO j=lbound(MatOFdnS,dim=2),ubound(MatOFdnS,dim=2)
                 mat(i,j) = MatOFdnS(i,j)%d3(id,jd,kd)
              END DO
              END DO
              write(out_unitp,*) 'd3',id,jd,kd
              CALL Write_VecMat(mat,out_unitp,5)
            END DO
            END DO
            END DO
          END IF
          CALL dealloc_array(mat,'mat','Write_MatOFdnS')
        END IF
        CALL flush_perso(out_unitp)
      END SUBROUTINE Write_MatOFdnS
!================================================================
!        dnS2 = dnS1 , dnVec2 = dnVec1 ...
!        transfer Vec(iVec) => R or R => Vec(iVec)
!================================================================


      SUBROUTINE sub_Mat1OFdnS_TO_Mat2OFdnS(Mat1OFdnS,Mat2OFdnS,nderiv)
        TYPE (Type_dnS) :: Mat1OFdnS(:,:),Mat2OFdnS(:,:)
        integer, optional :: nderiv

         integer :: i,j,nderiv_loc
        character (len=*), parameter :: name_sub='sub_Mat1OFdnS_TO_Mat2OFdnS'


        CALL check_alloc_MatOFdnS(Mat1OFdnS,'Mat1OFdnS',name_sub)
        CALL check_alloc_MatOFdnS(Mat2OFdnS,'Mat2OFdnS',name_sub)

        nderiv_loc = min(minval(Mat1OFdnS%nderiv),minval(Mat2OFdnS%nderiv))
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (Mat1OFdnS(1,1)%nb_var_deriv /= Mat2OFdnS(1,1)%nb_var_deriv) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nb_var_deriv in Mat1OFdnS and Mat2OFdnS are different!',&
                    Mat1OFdnS(1,1)%nb_var_deriv,Mat2OFdnS(1,1)%nb_var_deriv
          STOP
        END IF

        DO i=lbound(Mat1OFdnS,dim=1),ubound(Mat1OFdnS,dim=1)
        DO j=lbound(Mat1OFdnS,dim=2),ubound(Mat1OFdnS,dim=2)
           CALL sub_dnS1_TO_dnS2(Mat1OFdnS(i,j),Mat2OFdnS(i,j),nderiv_loc)
        END DO
        END DO

      END SUBROUTINE sub_Mat1OFdnS_TO_Mat2OFdnS

      SUBROUTINE DET_Mat3x3OFdnS_TO_dnS(MatOFdnS,dnDet,nderiv)
        TYPE (Type_dnS), intent (in)     :: MatOFdnS(3,3)
        TYPE (Type_dnS), intent (inout)  :: dnDet
        integer, optional   :: nderiv

        integer :: nderiv_loc
        TYPE (Type_dnS)  :: dnW1,dnW2,dnW3

        character (len=*), parameter :: name_sub='DET_Mat3x3OFdnS_TO_dnS'

        CALL check_alloc_MatOFdnS(MatOFdnS,'MatOFdnS',name_sub)
        IF (.NOT. dnDet%alloc) THEN
          CALL alloc_dnS(dnDet,MatOFdnS(1,1)%nb_var_deriv,MatOFdnS(1,1)%nderiv)
        END IF

        nderiv_loc = min(minval(MatOFdnS%nderiv),dnDet%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (MatOFdnS(1,1)%nb_var_deriv /= dnDet%nb_var_deriv) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nb_var_deriv in MatOFdnS and dnDet are different!',&
                    MatOFdnS(1,1)%nb_var_deriv,dnDet%nb_var_deriv
          STOP
        END IF

        ! check the determniant of MatOFdnS (and its derivatives)
        !first line
        CALL sub_dnS1_PROD_dnS2_TO_dnS3(MatOFdnS(2,2),MatOFdnS(3,3),dnW1)
        CALL sub_dnS1_PROD_dnS2_TO_dnS3(MatOFdnS(2,3),MatOFdnS(3,2),dnW2)
        CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnW1,ONE,dnW2,-ONE,dnW3)
        CALL sub_dnS1_PROD_dnS2_TO_dnS3(MatOFdnS(1,1),dnW3,dnDet)

        !second line
        CALL sub_dnS1_PROD_dnS2_TO_dnS3(MatOFdnS(1,2),MatOFdnS(3,3),dnW1)
        CALL sub_dnS1_PROD_dnS2_TO_dnS3(MatOFdnS(1,3),MatOFdnS(3,2),dnW2)
        CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnW1,ONE,dnW2,-ONE,dnW3)
        CALL sub_dnS1_PROD_dnS2_TO_dnS3(MatOFdnS(2,1),dnW3,dnW2)
        CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnDet,ONE,dnW2,-ONE,dnDet)

        !third line
        CALL sub_dnS1_PROD_dnS2_TO_dnS3(MatOFdnS(1,2),MatOFdnS(2,3),dnW1)
        CALL sub_dnS1_PROD_dnS2_TO_dnS3(MatOFdnS(1,3),MatOFdnS(2,2),dnW2)
        CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnW1,ONE,dnW2,-ONE,dnW3)
        CALL sub_dnS1_PROD_dnS2_TO_dnS3(MatOFdnS(3,1),dnW3,dnW2)
        CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnDet,ONE,dnW2, ONE,dnDet)

        !write(out_unitp,*) 'det(MatOFdnS)'
        !CALL Write_dnS(dnDet)

        CALL dealloc_dnS(dnW1)
        CALL dealloc_dnS(dnW2)
        CALL dealloc_dnS(dnW3)


      END SUBROUTINE DET_Mat3x3OFdnS_TO_dnS


      SUBROUTINE TRANS_Mat1OFdnS_TO_Mat2OFdnS(Mat1OFdnS,Mat2OFdnS,nderiv)
        TYPE (Type_dnS)     :: Mat1OFdnS(:,:),Mat2OFdnS(:,:)
        integer, optional   :: nderiv

        integer :: i,j


        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='TRANS_Mat1OFdnS_TO_Mat2OFdnS'

        CALL check_alloc_MatOFdnS(Mat1OFdnS,'Mat1OFdnS',name_sub)
        IF (.NOT. Mat2OFdnS(1,1)%alloc) THEN
          CALL alloc_MatOFdnS(Mat2OFdnS,Mat1OFdnS(1,1)%nb_var_deriv,Mat1OFdnS(1,1)%nderiv)
        END IF

        nderiv_loc = min(minval(Mat1OFdnS%nderiv),minval(Mat2OFdnS%nderiv))
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (Mat1OFdnS(1,1)%nb_var_deriv /= Mat2OFdnS(1,1)%nb_var_deriv) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nb_var_deriv in Mat1OFdnS and Mat2OFdnS are different!',&
                    Mat1OFdnS(1,1)%nb_var_deriv,Mat2OFdnS(1,1)%nb_var_deriv
          STOP
        END IF

        DO i=lbound(Mat1OFdnS,dim=2),ubound(Mat1OFdnS,dim=2)
        DO j=lbound(Mat1OFdnS,dim=1),ubound(Mat1OFdnS,dim=1)
           CALL sub_dnS1_TO_dnS2(Mat1OFdnS(j,i),Mat2OFdnS(i,j),nderiv_loc)
        END DO
        END DO

      END SUBROUTINE TRANS_Mat1OFdnS_TO_Mat2OFdnS

      SUBROUTINE DIAG_MatOFdnS(MatdnS,EigVecdnS,nderiv,sort,type_diago,phase,type_cs)
        TYPE (Type_dnS)     :: MatdnS(:,:),EigVecdnS(:,:)
        integer, optional   :: nderiv,sort,type_diago,type_cs
        logical, optional   :: phase


        integer             :: sort_loc,nderiv_loc,type_diago_loc,type_cs_loc
        logical             :: phase_loc

        integer             :: N

        TYPE (Type_dnS), pointer :: MatPdnS(:,:),MatWorkdnS(:,:)
        TYPE (Type_dnS), pointer :: MatdnS_save(:,:),Vec1dnS(:),Vec2dnS(:)

        integer, parameter           :: max_it = 200
        real (kind=Rkind), parameter :: tresh  = tiny(ONE)
        !real (kind=Rkind), parameter :: tresh  = TEN**(-30)
        real (kind=Rkind) :: d0,d1,d2,d3
        integer :: i,j

        !logical, parameter           :: check  = .TRUE.
        logical, parameter           :: check  = .FALSE.

!------ for debuging -------------------------------------------
        integer :: err_mem,memory
        logical, parameter :: debug=.FALSE.
!        logical, parameter :: debug=.TRUE.
        character (len=*), parameter :: name_sub='DIAG_MatOFdnS'
!---------------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
        END IF

        CALL check_alloc_MatOFdnS(MatdnS,'MatdnS',name_sub)
        IF (.NOT. EigVecdnS(1,1)%alloc) THEN
          CALL alloc_MatOFdnS(EigVecdnS,MatdnS(1,1)%nb_var_deriv,MatdnS(1,1)%nderiv)
        END IF
        IF (MatdnS(1,1)%nb_var_deriv /= EigVecdnS(1,1)%nb_var_deriv) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nb_var_deriv in MatdnS and EigVecdnS are different!',&
                    MatdnS(1,1)%nb_var_deriv,EigVecdnS(1,1)%nb_var_deriv
          STOP
        END IF

        nderiv_loc = min(MatdnS(1,1)%nderiv,EigVecdnS(1,1)%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        sort_loc = 0
        IF (present(sort)) sort_loc = sort

        type_diago_loc = -1
        IF (present(type_diago)) type_diago_loc = type_diago

        phase_loc = .FALSE.
        IF (present(phase)) phase_loc = phase

        IF (present(type_cs)) THEN
          type_cs_loc = type_cs
        ELSE
          type_cs_loc = 0
        END IF

        IF (check) THEN
          N = ubound(MatdnS,dim=2)
          nullify(MatdnS_save)
          CALL alloc_array(MatdnS_save,(/N,N/),'MatdnS_save',name_sub)
          CALL alloc_MatOFdnS(MatdnS_save,MatdnS(1,1)%nb_var_deriv,nderiv_loc)
          CALL sub_Mat1OFdnS_TO_Mat2OFdnS(MatdnS,MatdnS_save)
        END IF



        SELECT CASE(type_diago_loc)
        CASE(0)
          CALL DIAG01_MatOFdnS(MatdnS,EigVecdnS,nderiv_loc,max_it,tresh,type_cs_loc)

        CASE(1)
          CALL DIAG1_MatOFdnS(MatdnS,EigVecdnS,nderiv_loc,max_it,tresh,type_cs_loc)

        CASE(2)
          CALL DIAG2_MatOFdnS(MatdnS,EigVecdnS,nderiv_loc,max_it,tresh,type_cs_loc)

        CASE(3)
          CALL DIAG3_MatOFdnS(MatdnS,EigVecdnS,nderiv_loc,max_it,tresh,type_cs_loc)

        CASE(4)
          CALL DIAG4_MatOFdnS(MatdnS,EigVecdnS,nderiv_loc,max_it,tresh,type_cs_loc)

        CASE DEFAULT
          CALL DIAG3_MatOFdnS(MatdnS,EigVecdnS,nderiv_loc,max_it,tresh,type_cs_loc)

        END SELECT

        ! normalization of the vecors
        DO i=1,ubound(MatdnS,dim=2)
          CALL NORMALIZATION_OF_VecOFdnS(EigVecdnS(:,i),nderiv_loc)
        END DO


        IF (sort_loc == 1) CALL SORT_EigVectMat_OF_dnS(MatdnS,EigVecdnS,nderiv_loc)

        IF (phase_loc) CALL PHASE_EigVectMat_OF_dnS(EigVecdnS,nderiv_loc)

        IF (check) THEN
          ! allocation
          nullify(MatPdnS)
          CALL alloc_array(MatPdnS,(/N,N/),'MatPdnS',name_sub)
          nullify(MatWorkdnS)
          CALL alloc_array(MatWorkdnS,(/N,N/),'MatWorkdnS',name_sub)

          write(out_unitp,*) '=========================================='
          write(out_unitp,*) '======= CHECK DIAGO ======================'
          write(out_unitp,*) '=========================================='

          ! check line or column
          nullify(Vec1dnS)
          CALL alloc_array(Vec1dnS,(/N/),'Vec1dnS',name_sub)
          CALL alloc_VecOFdnS(Vec1dnS,MatdnS(1,1)%nb_var_deriv,nderiv_loc)
          nullify(Vec2dnS)
          CALL alloc_array(Vec2dnS,(/N/),'Vec2dnS',name_sub)
          CALL alloc_VecOFdnS(Vec2dnS,MatdnS(1,1)%nb_var_deriv,nderiv_loc)


          CALL MatOFdnS_TO_VecOFdnS(EigVecdnS,Vec1dnS,1,.TRUE.,nderiv_loc)
          CALL Mat1OFdnS_MUL_Vec2OFdnS_TO_Vec3OFdnS(MatdnS_save,        &
                                             Vec1dnS,Vec2dnS,nderiv_loc)
          CALL Vec1OFdnS_DOTPRODUCT_Vec2OFdnS_TO_dnS3(Vec1dnS,Vec2dnS,  &
                                              MatWorkdnS(1,1),nderiv_loc)
          write(out_unitp,*) 'vec #1 in line?',MatWorkdnS(1,1)%d0,MatdnS(1,1)%d0
          CALL flush_perso(out_unitp)

          CALL MatOFdnS_TO_VecOFdnS(EigVecdnS,Vec1dnS,1,.FALSE.,nderiv_loc)
          CALL Mat1OFdnS_MUL_Vec2OFdnS_TO_Vec3OFdnS(MatdnS_save,        &
                                             Vec1dnS,Vec2dnS,nderiv_loc)
          CALL Vec1OFdnS_DOTPRODUCT_Vec2OFdnS_TO_dnS3(Vec1dnS,Vec2dnS,  &
                                              MatWorkdnS(1,1),nderiv_loc)
          write(out_unitp,*) 'vec #1 in column?',MatWorkdnS(1,1)%d0,MatdnS(1,1)%d0
          CALL flush_perso(out_unitp)



          write(out_unitp,*) 'diagonal Matrix of dnS'
          CALL Write_MatOFdnS(MatdnS,nderiv)
          !write(out_unitp,*) 'Eigenvector (in line) Matrix of dnS'
          !CALL Write_MatOFdnS(EigVecdnS,nderiv)


          CALL TRANS_Mat1OFdnS_TO_Mat2OFdnS(EigVecdnS,MatPdnS,nderiv_loc)
          CALL Mat1OFdnS_MUL_Mat2OFdnS_TO_Mat3OFdnS(MatdnS,             &
                                              MatPdnS,MatWorkdnS,nderiv_loc)
          CALL Mat1OFdnS_MUL_Mat2OFdnS_TO_Mat3OFdnS(EigVecdnS,          &
                                               MatWorkdnS,MatdnS,nderiv_loc)

          CALL Mat1OFdnS_wPLUS_Mat2OFdnS_TO_Mat3OFdnS(MatdnS,ONE,       &
                                     MatdnS_save,-ONE,MatWorkdnS,nderiv_loc)

          !write(out_unitp,*) 'initial Matrix of dnS'
          !CALL Write_MatOFdnS(MatdnS,nderiv)

          d0 = ZERO
          d1 = ZERO
          d2 = ZERO
          d3 = ZERO
          DO i=1,ubound(MatWorkdnS,dim=1)
          DO j=1,ubound(MatWorkdnS,dim=2)
            IF (abs(MatWorkdnS(i,j)%d0) > d0) d0 = abs(MatWorkdnS(i,j)%d0)
            IF (nderiv_loc == 0) CYCLE

            IF (maxval(abs(MatWorkdnS(i,j)%d1)) > d1) d1 = maxval(abs(MatWorkdnS(i,j)%d1))
            IF (nderiv_loc == 1) CYCLE

            IF (maxval(abs(MatWorkdnS(i,j)%d2)) > d2) d2 = maxval(abs(MatWorkdnS(i,j)%d2))
            IF (nderiv_loc == 2) CYCLE

            IF (maxval(abs(MatWorkdnS(i,j)%d3)) > d3) d3 = maxval(abs(MatWorkdnS(i,j)%d3))

          END DO
          END DO
          write(out_unitp,*) 'diff mat',d0,d1,d2,d3
          !CALL Write_MatOFdnS(MatWorkdnS,nderiv_loc)


!         write(out_unitp,*) 'diff mat',MatdnS(:,:)%d0-MatdnS_save(:,:)%d0
          CALL dealloc_MatOFdnS(MatdnS_save)
          CALL dealloc_array(MatdnS_save,'MatdnS_save',name_sub)
          CALL dealloc_MatOFdnS(MatPdnS)
          CALL dealloc_array(MatPdnS,'MatPdnS',name_sub)
          CALL dealloc_MatOFdnS(MatWorkdnS)
          CALL dealloc_array(MatWorkdnS,'MatWorkdnS',name_sub)
          CALL dealloc_VecOFdnS(Vec1dnS)
          CALL dealloc_array(Vec1dnS,'Vec1dnS',name_sub)
          CALL dealloc_VecOFdnS(Vec2dnS)
          CALL dealloc_array(Vec2dnS,'Vec2dnS',name_sub)

          !write(out_unitp,*) 'STOP: check DIAG ', name_sub
          !STOP
          write(out_unitp,*) '=========================================='
          write(out_unitp,*) '=========================================='

        END IF




        IF (debug) THEN
          write(out_unitp,*) 'diagonal Matrix of dnS'
          CALL Write_MatOFdnS(MatdnS,nderiv)
          write(out_unitp,*) 'Eigenvector (in line) Matrix of dnS'
          CALL Write_MatOFdnS(EigVecdnS,nderiv)

          write(out_unitp,*) 'END ',name_sub
        END IF

      END SUBROUTINE DIAG_MatOFdnS

      SUBROUTINE DIAG4_MatOFdnS(MatdnS,EigVecdnS,nderiv,max_it,tresh,type_cs)
        TYPE (Type_dnS)                :: MatdnS(:,:),EigVecdnS(:,:)
        integer, intent (in)           :: nderiv,max_it
        real (kind=Rkind), intent (in) :: tresh
        integer, intent(in)            :: type_cs

        TYPE (Type_dnS), pointer :: MatPdnS(:,:),MatWorkdnS(:,:)
        TYPE (Type_dnS), pointer :: MatdnS_save(:,:)


        integer             :: it,i,j,ip,iq,max_ip,max_iq,N

        logical             :: conv
        real (kind=Rkind)   :: max_matpq,th,cth,sth,x,y,max_diag,max_dnSmatpq
        TYPE (Type_dnS)     :: dnTh,dnCos,dnSin,dnWork1,dnWork2,dnWork3



!------ for debuging -------------------------------------------
        integer :: err_mem,memory
        logical, parameter :: debug=.FALSE.
!        logical, parameter :: debug=.TRUE.
        character (len=*), parameter :: name_sub='DIAG4_MatOFdnS'
!---------------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
        END IF

        ! allocation
        N = ubound(MatdnS,dim=2)
        nullify(MatPdnS)
        CALL alloc_array(MatPdnS,(/N,N/),'MatPdnS',name_sub)
        CALL alloc_MatOFdnS(MatPdnS,MatdnS(1,1)%nb_var_deriv,nderiv)
        nullify(MatWorkdnS)
        CALL alloc_array(MatWorkdnS,(/N,N/),'MatWorkdnS',name_sub)
        CALL alloc_MatOFdnS(MatWorkdnS,MatdnS(1,1)%nb_var_deriv,nderiv)

        ! initialization
        CALL sub_Id_TO_MatOFdnS(EigVecdnS)

        IF (debug) THEN
           write(out_unitp,*) 'Initial Matrix of dnS'
           CALL Write_MatOFdnS(MatdnS,nderiv)
         END IF

        DO it=1,max_it
          conv      = .FALSE.
          max_matpq = ZERO
          DO ip=1,N
          DO iq=ip+1,N
            CALL Calc_dnCos_AND_dnSin(dnCos,dnSin,                      &
                                      MatdnS(ip,ip),                    &
                                      MatdnS(iq,iq),                    &
                                      MatdnS(ip,iq),                    &
                                      dnTh,dnWork1,dnWork2,dnWork3,type_cs)
            CALL Calc_MaxVal_OF_dnS(dnTh,max_dnSmatpq)
            !write(6,*) 'it,ip,iq,max_dnSmatpq',it,ip,iq,max_dnSmatpq
            IF (max_dnSmatpq < tresh) CYCLE
            IF (max_dnSmatpq > max_matpq)  max_matpq = max_dnSmatpq

            !-----------------------------------------------------------
            ! rotational matrix
            CALL sub_Id_TO_MatOFdnS(MatPdnS)

            !MatPdnS(max_ip,max_ip)%d0 =  cth
            CALL sub_dnS1_TO_dnS2(dnCos,MatPdnS(ip,ip),nderiv)
            !MatPdnS(max_iq,max_iq)%d0 =  cth
            CALL sub_dnS1_TO_dnS2(dnCos,MatPdnS(iq,iq),nderiv)
            !MatPdnS(max_ip,max_iq)%d0 = -sth
            CALL sub_dnS1_PROD_w_TO_dnS2(dnSin,-ONE,MatPdnS(ip,iq),nderiv)
            !MatPdnS(max_iq,max_ip)%d0 =  sth
            CALL sub_dnS1_TO_dnS2(dnSin,MatPdnS(iq,ip),nderiv)

            IF (debug) THEN
              !write(out_unitp,*) 'Rotational matrix',it
              !CALL Write_MatOFdnS(MatPdnS,nderiv=0)
            END IF
            !-----------------------------------------------------------

            !-----------------------------------------------------------
            ! eigenvectors
            CALL Mat1OFdnS_MUL_Mat2OFdnS_TO_Mat3OFdnS(EigVecdnS,MatPdnS,MatWorkdnS)
            CALL sub_Mat1OFdnS_TO_Mat2OFdnS(MatWorkdnS,EigVecdnS,nderiv)
            !-----------------------------------------------------------

            !-----------------------------------------------------------
            ! matrix
            CALL Mat1OFdnS_MUL_Mat2OFdnS_TO_Mat3OFdnS(MatdnS,MatPdnS,MatWorkdnS)
            CALL sub_Mat1OFdnS_TO_Mat2OFdnS(MatWorkdnS,MatdnS,nderiv)

            !MatPdnS(max_iq,max_ip)%d0 = -sth
            CALL sub_dnS1_PROD_w_TO_dnS2(dnSin,-ONE,MatPdnS(iq,ip),nderiv)
            !MatPdnS(max_ip,max_iq)%d0 =  sth
            CALL sub_dnS1_TO_dnS2(dnSin,MatPdnS(ip,iq),nderiv)

            CALL Mat1OFdnS_MUL_Mat2OFdnS_TO_Mat3OFdnS(MatPdnS,MatdnS,MatWorkdnS)
            CALL sub_Mat1OFdnS_TO_Mat2OFdnS(MatWorkdnS,MatdnS,nderiv)
            CALL sub_ZERO_TO_dnS(MatdnS(ip,iq),nderiv)
            CALL sub_ZERO_TO_dnS(MatdnS(iq,ip),nderiv)

            IF (debug) THEN
              write(out_unitp,*) 'Intermediate Matrix of dnS',it
              !CALL Write_MatOFdnS(MatdnS,nderiv=0)
              CALL Write_MatOFdnS(MatdnS)
            END IF
            !-----------------------------------------------------------

          END DO
          END DO
          conv = (max_matpq < tresh)
          IF (conv) EXIT

        END DO
        !write(out_unitp,*) 'it',it
        IF (.NOT. conv) THEN
          write(out_unitp,*) 'diagonal Matrix of dnS (not converged)'
          CALL Write_MatOFdnS(MatdnS,nderiv)
          write(out_unitp,*) 'Eigenvector (in line) Matrix of dnS (not converged)'
          CALL Write_MatOFdnS(EigVecdnS,nderiv)
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) max_it,' iterations should never happen'
          STOP
        END IF

        ! deallocation
        CALL dealloc_dnS(dnTh)
        CALL dealloc_dnS(dnCos)
        CALL dealloc_dnS(dnSin)
        CALL dealloc_dnS(dnWork1)
        CALL dealloc_dnS(dnWork2)
        CALL dealloc_dnS(dnWork3)
        CALL dealloc_MatOFdnS(MatPdnS)
        CALL dealloc_MatOFdnS(MatWorkdnS)
        CALL dealloc_array(MatPdnS,'MatPdnS',name_sub)
        CALL dealloc_array(MatWorkdnS,'MatWorkdnS',name_sub)

        IF (debug) THEN
          write(out_unitp,*) 'diagonal Matrix of dnS'
          CALL Write_MatOFdnS(MatdnS,nderiv)
          write(out_unitp,*) 'Eigenvector (in line) Matrix of dnS'
          CALL Write_MatOFdnS(EigVecdnS,nderiv)
          write(out_unitp,*) 'END ',name_sub
        END IF


      END SUBROUTINE DIAG4_MatOFdnS


      SUBROUTINE DIAG3_MatOFdnS(MatdnS,EigVecdnS,nderiv,max_it,tresh,type_cs)
        TYPE (Type_dnS)     :: MatdnS(:,:),EigVecdnS(:,:)
        integer, intent (in)           :: nderiv,max_it
        real (kind=Rkind), intent (in) :: tresh
        integer, intent(in)            :: type_cs

        TYPE (Type_dnS), pointer :: MatPdnS(:,:),MatWorkdnS(:,:)
        TYPE (Type_dnS), pointer :: MatdnS_save(:,:)


        integer             :: it,i,j,ip,iq,max_ip,max_iq,N

        logical             :: conv
        real (kind=Rkind)   :: max_matpq,th,cth,sth,x,y,max_diag,max_dnSmatpq
        TYPE (Type_dnS)     :: dnTh,dnCos,dnSin,dnWork1,dnWork2,dnWork3



!------ for debuging -------------------------------------------
        integer :: err_mem,memory
        logical, parameter :: debug=.FALSE.
!        logical, parameter :: debug=.TRUE.
        character (len=*), parameter :: name_sub='DIAG3_MatOFdnS'
!---------------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
        END IF

        ! allocation
        N = ubound(MatdnS,dim=2)
        nullify(MatPdnS)
        CALL alloc_array(MatPdnS,(/N,N/),'MatPdnS',name_sub)
        CALL alloc_MatOFdnS(MatPdnS,MatdnS(1,1)%nb_var_deriv,nderiv)
        nullify(MatWorkdnS)
        CALL alloc_array(MatWorkdnS,(/N,N/),'MatWorkdnS',name_sub)
        CALL alloc_MatOFdnS(MatWorkdnS,MatdnS(1,1)%nb_var_deriv,nderiv)


        ! initialization
        CALL sub_Id_TO_MatOFdnS(EigVecdnS)

        IF (debug) THEN
           write(out_unitp,*) 'Initial Matrix of dnS'
           CALL Write_MatOFdnS(MatdnS,nderiv)
         END IF

        DO it=1,max_it
          conv      = .FALSE.
          max_matpq = ZERO
          DO ip=1,N
          DO iq=ip+1,N
            CALL Calc_MaxVal_OF_dnS(MatdnS(ip,iq),max_dnSmatpq)
            IF (max_dnSmatpq < TEN*tresh) CYCLE
            CALL Calc_dnCos_AND_dnSin(dnCos,dnSin,                      &
                                      MatdnS(ip,ip),                    &
                                      MatdnS(iq,iq),                    &
                                      MatdnS(ip,iq),                    &
                                      dnTh,dnWork1,dnWork2,dnWork3,type_cs)
            CALL Calc_MaxVal_OF_dnS(dnTh,max_dnSmatpq)
            !max_dnSmatpq = abs(dnTh%d0)
            IF (max_dnSmatpq > max_matpq) THEN
              max_matpq = max_dnSmatpq
              max_ip    = ip
              max_iq    = iq
            END IF
          END DO
          END DO
          conv = (max_matpq < tresh)
          IF (conv) EXIT

          ip = max_ip
          iq = max_iq

          CALL Calc_dnCos_AND_dnSin(dnCos,dnSin,                        &
                                    MatdnS(ip,ip),                      &
                                    MatdnS(iq,iq),                      &
                                    MatdnS(ip,iq),                      &
                                    dnTh,dnWork1,dnWork2,dnWork3,type_cs)
          IF (debug) THEN
            write(out_unitp,*) 'iteration',it
            write(out_unitp,*) 'max_matpq',max_ip,max_iq,max_matpq
            write(out_unitp,*) 'th',dnTh%d0
            CALL Write_dnS(dnth)
          END IF

          ! rotational matrix
          CALL sub_Id_TO_MatOFdnS(MatPdnS)

          !MatPdnS(max_ip,max_ip)%d0 =  cth
          CALL sub_dnS1_TO_dnS2(dnCos,MatPdnS(max_ip,max_ip),nderiv)
          !MatPdnS(max_iq,max_iq)%d0 =  cth
          CALL sub_dnS1_TO_dnS2(dnCos,MatPdnS(max_iq,max_iq),nderiv)
          !MatPdnS(max_ip,max_iq)%d0 = -sth
          CALL sub_dnS1_PROD_w_TO_dnS2(dnSin,-ONE,MatPdnS(max_ip,max_iq),nderiv)
          !MatPdnS(max_iq,max_ip)%d0 =  sth
          CALL sub_dnS1_TO_dnS2(dnSin,MatPdnS(max_iq,max_ip),nderiv)

          IF (debug) THEN
            !write(out_unitp,*) 'Rotational matrix',it
            !CALL Write_MatOFdnS(MatPdnS,nderiv=0)
          END IF

          ! eigenvectors
          CALL Mat1OFdnS_MUL_Mat2OFdnS_TO_Mat3OFdnS(EigVecdnS,MatPdnS,MatWorkdnS)
          CALL sub_Mat1OFdnS_TO_Mat2OFdnS(MatWorkdnS,EigVecdnS,nderiv)

          ! matrix
          CALL Mat1OFdnS_MUL_Mat2OFdnS_TO_Mat3OFdnS(MatdnS,MatPdnS,MatWorkdnS)
          CALL sub_Mat1OFdnS_TO_Mat2OFdnS(MatWorkdnS,MatdnS,nderiv)

          !MatPdnS(max_iq,max_ip)%d0 = -sth
          CALL sub_dnS1_PROD_w_TO_dnS2(dnSin,-ONE,MatPdnS(max_iq,max_ip),nderiv)
          !MatPdnS(max_ip,max_iq)%d0 =  sth
          CALL sub_dnS1_TO_dnS2(dnSin,MatPdnS(ip,iq),nderiv)

          CALL Mat1OFdnS_MUL_Mat2OFdnS_TO_Mat3OFdnS(MatPdnS,MatdnS,MatWorkdnS)
          CALL sub_Mat1OFdnS_TO_Mat2OFdnS(MatWorkdnS,MatdnS,nderiv)
          CALL sub_ZERO_TO_dnS(MatdnS(ip,iq),nderiv)
          CALL sub_ZERO_TO_dnS(MatdnS(iq,ip),nderiv)

          IF (debug) THEN
            write(out_unitp,*) 'Intermediate Matrix of dnS',it
            !CALL Write_MatOFdnS(MatdnS,nderiv=0)
            CALL Write_MatOFdnS(MatdnS)
          END IF

        END DO
        !write(out_unitp,*) 'it',it

        IF (.NOT. conv) THEN
          write(out_unitp,*) 'diagonal Matrix of dnS (not converged)'
          CALL Write_MatOFdnS(MatdnS,nderiv)
          write(out_unitp,*) 'Eigenvector (in line) Matrix of dnS (not converged)'
          CALL Write_MatOFdnS(EigVecdnS,nderiv)
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) max_it,' iterations should never happen'
          STOP
        END IF

        IF (debug) THEN
            write(out_unitp,*) 'iteration',it,'in ',name_sub
            write(out_unitp,*) 'max_matpq',max_ip,max_iq,max_matpq
            write(out_unitp,*) 'th',dnTh%d0
        END IF


        ! deallocation
        CALL dealloc_dnS(dnTh)
        CALL dealloc_dnS(dnCos)
        CALL dealloc_dnS(dnSin)
        CALL dealloc_dnS(dnWork1)
        CALL dealloc_dnS(dnWork2)
        CALL dealloc_dnS(dnWork3)
        CALL dealloc_MatOFdnS(MatPdnS)
        CALL dealloc_MatOFdnS(MatWorkdnS)
        CALL dealloc_array(MatPdnS,'MatPdnS',name_sub)
        CALL dealloc_array(MatWorkdnS,'MatWorkdnS',name_sub)

        IF (debug) THEN
          write(out_unitp,*) 'diagonal Matrix of dnS'
          CALL Write_MatOFdnS(MatdnS,nderiv)
          write(out_unitp,*) 'Eigenvector (in line) Matrix of dnS'
          CALL Write_MatOFdnS(EigVecdnS,nderiv)
          write(out_unitp,*) 'END ',name_sub
        END IF


      END SUBROUTINE DIAG3_MatOFdnS

      SUBROUTINE DIAG2_MatOFdnS(MatdnS,EigVecdnS,nderiv,max_it,tresh,type_cs)
        TYPE (Type_dnS)     :: MatdnS(:,:),EigVecdnS(:,:)
        integer, intent (in)           :: nderiv,max_it
        real (kind=Rkind), intent (in) :: tresh
        integer, intent(in)            :: type_cs

        TYPE (Type_dnS), pointer :: MatPdnS(:,:),MatWorkdnS(:,:)
        TYPE (Type_dnS), pointer :: MatdnS_save(:,:)


        integer             :: it,i,j,ip,iq,max_ip,max_iq,N

        logical             :: conv
        real (kind=Rkind)   :: max_matpq,th,cth,sth,x,y,max_diag,max_dnSmatpq,diff_ppqq
        TYPE (Type_dnS)     :: dnTh,dnCos,dnSin,dnWork1,dnWork2,dnWork3

!------ for debuging -------------------------------------------
        integer :: err_mem,memory
        logical, parameter :: debug=.FALSE.
!        logical, parameter :: debug=.TRUE.
        character (len=*), parameter :: name_sub='DIAG2_MatOFdnS'
!---------------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
        END IF

        ! allocation
        N = ubound(MatdnS,dim=2)
        CALL alloc_array(MatPdnS,(/N,N/),'MatPdnS',name_sub)
        CALL alloc_MatOFdnS(MatPdnS,MatdnS(1,1)%nb_var_deriv,nderiv)
        CALL alloc_array(MatWorkdnS,(/N,N/),'MatWorkdnS',name_sub)
        CALL alloc_MatOFdnS(MatWorkdnS,MatdnS(1,1)%nb_var_deriv,nderiv)


        ! initialization
        CALL sub_Id_TO_MatOFdnS(EigVecdnS)

        IF (debug) THEN
           write(out_unitp,*) 'Initial Matrix of dnS'
           CALL Write_MatOFdnS(MatdnS,nderiv)
         END IF

        DO it=1,max_it
          conv      = .FALSE.
          max_matpq = ZERO
          max_diag  = ZERO
          DO ip=1,N

            IF (abs(MatdnS(ip,ip)%d0-MatdnS(1,1)%d0) > max_diag)        &
                         max_diag = abs(MatdnS(ip,ip)%d0-MatdnS(1,1)%d0)

            DO iq=ip+1,N
              diff_ppqq = abs(MatdnS(ip,ip)%d0-MatdnS(iq,iq)%d0)

              IF (diff_ppqq > ONETENTH**6) THEN
                max_dnSmatpq = abs(MatdnS(ip,iq)%d0)/diff_ppqq
              ELSE
                max_dnSmatpq = TEN**6
              END IF
              !CALL Calc_MaxVal_OF_dnS(MatdnS(ip,iq),max_dnSmatpq)
              IF (max_dnSmatpq > max_matpq) THEN
                max_matpq = max_dnSmatpq
                max_ip    = ip
                max_iq    = iq
              END IF
            END DO
          END DO
          max_diag = max(max_diag,max_matpq)
          conv = (max_matpq/max_diag) < tresh
          conv = max_matpq < tresh
          IF (conv) EXIT

          ip = max_ip
          iq = max_iq

          IF (abs(MatdnS(ip,iq)%d0/(MatdnS(ip,ip)%d0-MatdnS(iq,iq)%d0)) &
              < ONETENTH) THEN
            CALL Calc2_dnCos_AND_dnSin(dnCos,dnSin,                     &
                                       MatdnS(ip,ip),                   &
                                       MatdnS(iq,iq),                   &
                                       MatdnS(ip,iq),                   &
                                       dnTh,dnWork1,dnWork2,dnWork3)
          ELSE
            CALL Calc1_dnCos_AND_dnSin(dnCos,dnSin,                     &
                                       MatdnS(ip,ip),                   &
                                       MatdnS(iq,iq),                   &
                                       MatdnS(ip,iq),                   &
                                       dnTh,dnWork1,dnWork2,dnWork3)
          END IF
          IF (debug) THEN
            write(out_unitp,*) 'max_matpq',max_ip,max_iq,max_matpq
            write(out_unitp,*) 'th',dnTh%d0
          END IF

          ! rotational matrix
          CALL sub_Id_TO_MatOFdnS(MatPdnS)

          !MatPdnS(max_ip,max_ip)%d0 =  cth
          CALL sub_dnS1_TO_dnS2(dnCos,MatPdnS(max_ip,max_ip),nderiv)
          !MatPdnS(max_iq,max_iq)%d0 =  cth
          CALL sub_dnS1_TO_dnS2(dnCos,MatPdnS(max_iq,max_iq),nderiv)
          !MatPdnS(max_ip,max_iq)%d0 = -sth
          CALL sub_dnS1_PROD_w_TO_dnS2(dnSin,-ONE,MatPdnS(max_ip,max_iq),nderiv)
          !MatPdnS(max_iq,max_ip)%d0 =  sth
          CALL sub_dnS1_TO_dnS2(dnSin,MatPdnS(max_iq,max_ip),nderiv)

          IF (debug) THEN
            write(out_unitp,*) 'Rotational matrix',it
            CALL Write_MatOFdnS(MatPdnS,nderiv=0)
          END IF

          ! eigenvectors
          CALL Mat1OFdnS_MUL_Mat2OFdnS_TO_Mat3OFdnS(EigVecdnS,MatPdnS,MatWorkdnS)
          CALL sub_Mat1OFdnS_TO_Mat2OFdnS(MatWorkdnS,EigVecdnS,nderiv)

          ! matrix
          CALL Mat1OFdnS_MUL_Mat2OFdnS_TO_Mat3OFdnS(MatdnS,MatPdnS,MatWorkdnS)
          CALL sub_Mat1OFdnS_TO_Mat2OFdnS(MatWorkdnS,MatdnS,nderiv)

          !MatPdnS(max_iq,max_ip)%d0 = -sth
          CALL sub_dnS1_PROD_w_TO_dnS2(dnSin,-ONE,MatPdnS(max_iq,max_ip),nderiv)
          !MatPdnS(max_ip,max_iq)%d0 =  sth
          CALL sub_dnS1_TO_dnS2(dnSin,MatPdnS(max_ip,max_iq),nderiv)

          CALL Mat1OFdnS_MUL_Mat2OFdnS_TO_Mat3OFdnS(MatPdnS,MatdnS,MatWorkdnS)
          CALL sub_Mat1OFdnS_TO_Mat2OFdnS(MatWorkdnS,MatdnS,nderiv)

          IF (debug) THEN
            write(out_unitp,*) 'Intermediate Matrix of dnS',it
            CALL Write_MatOFdnS(MatdnS,nderiv=0)
          END IF

        END DO

        IF (.NOT. conv) THEN
          write(out_unitp,*) 'diagonal Matrix of dnS (not converged)'
          CALL Write_MatOFdnS(MatdnS,nderiv)
          write(out_unitp,*) 'Eigenvector (in line) Matrix of dnS (not converged)'
          CALL Write_MatOFdnS(EigVecdnS,nderiv)
          write(out_unitp,*) ' ERROR in ',name_sub
          write(6,*) max_it,' iterations should never happen'
          STOP
        END IF


        ! deallocation
        CALL dealloc_dnS(dnTh)
        CALL dealloc_dnS(dnCos)
        CALL dealloc_dnS(dnSin)
        CALL dealloc_dnS(dnWork1)
        CALL dealloc_dnS(dnWork2)
        CALL dealloc_dnS(dnWork3)
        CALL dealloc_MatOFdnS(MatPdnS)
        CALL dealloc_MatOFdnS(MatWorkdnS)
        CALL dealloc_array(MatPdnS,'MatPdnS',name_sub)
        CALL dealloc_array(MatWorkdnS,'MatWorkdnS',name_sub)

      END SUBROUTINE DIAG2_MatOFdnS


      SUBROUTINE DIAG1_MatOFdnS(MatdnS,EigVecdnS,nderiv,max_it,tresh,type_cs)
        TYPE (Type_dnS)     :: MatdnS(:,:),EigVecdnS(:,:)
        integer, intent (in)           :: nderiv,max_it
        real (kind=Rkind), intent (in) :: tresh
        integer, intent(in)            :: type_cs

        TYPE (Type_dnS), pointer :: MatPdnS(:,:),MatWorkdnS(:,:)
        TYPE (Type_dnS), pointer :: MatdnS_save(:,:)


        integer             :: it,i,j,ip,iq,max_ip,max_iq,N

        logical             :: conv
        real (kind=Rkind)   :: max_matpq,th,cth,sth,x,y,max_diag,max_dnSmatpq
        TYPE (Type_dnS)     :: dnTh,dnCos,dnSin,dnWork1,dnWork2,dnWork3



!------ for debuging -------------------------------------------
        integer :: err_mem,memory
        logical, parameter :: debug=.FALSE.
!        logical, parameter :: debug=.TRUE.
        character (len=*), parameter :: name_sub='DIAG1_MatOFdnS'
!---------------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
        END IF


        ! allocation
        N = ubound(MatdnS,dim=2)
        IF (N==1) THEN
          write(out_unitp,*) 'WARNNING : The dimensin of the matrix is 1'
          write(out_unitp,*) 'No diagonalization!!'
          CALL sub_ZERO_TO_MatOFdnS(EigVecdnS)
          EigVecdnS(1,1)%d0 = ONE
          RETURN
        END IF

        nullify(MatPdnS)
        CALL alloc_array(MatPdnS,(/N,N/),'MatPdnS',name_sub)
        CALL alloc_MatOFdnS(MatPdnS,MatdnS(1,1)%nb_var_deriv,nderiv)
        nullify(MatWorkdnS)
        CALL alloc_array(MatWorkdnS,(/N,N/),'MatWorkdnS',name_sub)
        CALL alloc_MatOFdnS(MatWorkdnS,MatdnS(1,1)%nb_var_deriv,nderiv)


        ! initialization
        CALL sub_Id_TO_MatOFdnS(EigVecdnS)

        IF (debug) THEN
           write(out_unitp,*) 'Initial Matrix of dnS'
           CALL Write_MatOFdnS(MatdnS,nderiv)
         END IF

        DO it=1,max_it
          conv      = .FALSE.
          max_matpq = ZERO
          max_diag  = ZERO
          DO ip=1,N
            IF (abs(MatdnS(ip,ip)%d0-MatdnS(1,1)%d0) > max_diag)        &
                         max_diag = abs(MatdnS(ip,ip)%d0-MatdnS(1,1)%d0)
            DO iq=ip+1,N
              CALL Calc_MaxVal_OF_dnS(MatdnS(ip,iq),max_dnSmatpq)
              IF (max_dnSmatpq > max_matpq) THEN
                max_matpq = max_dnSmatpq
                max_ip    = ip
                max_iq    = iq
              END IF
            END DO
          END DO
          max_diag = max(max_diag,max_matpq)
          conv = (max_matpq/max_diag) < tresh
          IF (conv) EXIT

          ip = max_ip
          iq = max_iq

          CALL Calc_dnCos_AND_dnSin(dnCos,dnSin,                        &
                                    MatdnS(ip,ip),                      &
                                    MatdnS(iq,iq),                      &
                                    MatdnS(ip,iq),                      &
                                    dnTh,dnWork1,dnWork2,dnWork3,type_cs)
          IF (debug) THEN
            write(out_unitp,*) 'max_matpq',max_ip,max_iq,max_matpq
            write(out_unitp,*) 'th',dnTh%d0
            CALL Write_dnS(dnSin)
          END IF

          ! rotational matrix
          CALL sub_Id_TO_MatOFdnS(MatPdnS)

          !MatPdnS(max_ip,max_ip)%d0 =  cth
          CALL sub_dnS1_TO_dnS2(dnCos,MatPdnS(max_ip,max_ip),nderiv)
          !MatPdnS(max_iq,max_iq)%d0 =  cth
          CALL sub_dnS1_TO_dnS2(dnCos,MatPdnS(max_iq,max_iq),nderiv)
          !MatPdnS(max_ip,max_iq)%d0 = -sth
          CALL sub_dnS1_PROD_w_TO_dnS2(dnSin,-ONE,MatPdnS(max_ip,max_iq),nderiv)
          !MatPdnS(max_iq,max_ip)%d0 =  sth
          CALL sub_dnS1_TO_dnS2(dnSin,MatPdnS(max_iq,max_ip),nderiv)

          IF (debug) THEN
            !write(out_unitp,*) 'Rotational matrix',it
            !CALL Write_MatOFdnS(MatPdnS,nderiv=0)
          END IF

          ! eigenvectors
          CALL Mat1OFdnS_MUL_Mat2OFdnS_TO_Mat3OFdnS(EigVecdnS,MatPdnS,MatWorkdnS)
          CALL sub_Mat1OFdnS_TO_Mat2OFdnS(MatWorkdnS,EigVecdnS,nderiv)

          ! matrix
          CALL Mat1OFdnS_MUL_Mat2OFdnS_TO_Mat3OFdnS(MatdnS,MatPdnS,MatWorkdnS)
          CALL sub_Mat1OFdnS_TO_Mat2OFdnS(MatWorkdnS,MatdnS,nderiv)

          !MatPdnS(max_iq,max_ip)%d0 = -sth
          CALL sub_dnS1_PROD_w_TO_dnS2(dnSin,-ONE,MatPdnS(max_iq,max_ip),nderiv)
          !MatPdnS(max_ip,max_iq)%d0 =  sth
          CALL sub_dnS1_TO_dnS2(dnSin,MatPdnS(max_ip,max_iq),nderiv)

          CALL Mat1OFdnS_MUL_Mat2OFdnS_TO_Mat3OFdnS(MatPdnS,MatdnS,MatWorkdnS)
          CALL sub_Mat1OFdnS_TO_Mat2OFdnS(MatWorkdnS,MatdnS,nderiv)

          IF (debug) THEN
            !write(out_unitp,*) 'Intermediate Matrix of dnS',it
            !CALL Write_MatOFdnS(MatdnS,nderiv=0)
          END IF

        END DO

        IF (.NOT. conv) THEN
          write(out_unitp,*) 'diagonal Matrix of dnS (not converged)'
          CALL Write_MatOFdnS(MatdnS,nderiv)
          write(out_unitp,*) 'Eigenvector (in line) Matrix of dnS (not converged)'
          CALL Write_MatOFdnS(EigVecdnS,nderiv)
          write(out_unitp,*) ' ERROR in ',name_sub
          write(6,*) max_it,' iterations should never happen'
          STOP
        END IF

        ! deallocation
        CALL dealloc_dnS(dnTh)
        CALL dealloc_dnS(dnCos)
        CALL dealloc_dnS(dnSin)
        CALL dealloc_dnS(dnWork1)
        CALL dealloc_dnS(dnWork2)
        CALL dealloc_dnS(dnWork3)
        CALL dealloc_MatOFdnS(MatPdnS)
        CALL dealloc_MatOFdnS(MatWorkdnS)
        CALL dealloc_array(MatPdnS,'MatPdnS',name_sub)
        CALL dealloc_array(MatWorkdnS,'MatWorkdnS',name_sub)

      END SUBROUTINE DIAG1_MatOFdnS

      SUBROUTINE DIAG01_MatOFdnS(MatdnS,EigVecdnS,nderiv,max_it,tresh,type_cs)
        TYPE (Type_dnS)     :: MatdnS(:,:),EigVecdnS(:,:)
        integer, intent (in)           :: nderiv,max_it
        real (kind=Rkind), intent (in) :: tresh
        integer, intent(in)            :: type_cs

        TYPE (Type_dnS), pointer :: MatPdnS(:,:),MatWorkdnS(:,:)
        TYPE (Type_dnS), pointer :: MatdnS_save(:,:)


        integer             :: it,i,j,ip,iq,max_ip,max_iq,N

        logical             :: conv
        real (kind=Rkind)   :: max_matpq,th,cth,sth,x,y,max_diag,max_dnSmatpq
        TYPE (Type_dnS)     :: dnTh,dnCos,dnSin,dnWork1,dnWork2,dnWork3

!------ for debuging -------------------------------------------
        integer :: err_mem,memory
        logical, parameter :: debug=.FALSE.
!        logical, parameter :: debug=.TRUE.
        character (len=*), parameter :: name_sub='DIAG01_MatOFdnS'
!---------------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
        END IF

        ! allocation
        N = ubound(MatdnS,dim=2)
        nullify(MatPdnS)
        CALL alloc_array(MatPdnS,(/N,N/),'MatPdnS',name_sub)
        CALL alloc_MatOFdnS(MatPdnS,MatdnS(1,1)%nb_var_deriv,nderiv)
        nullify(MatWorkdnS)
        CALL alloc_array(MatWorkdnS,(/N,N/),'MatWorkdnS',name_sub)
        CALL alloc_MatOFdnS(MatWorkdnS,MatdnS(1,1)%nb_var_deriv,nderiv)

        ! initialization
        CALL sub_Id_TO_MatOFdnS(EigVecdnS)

        IF (debug) THEN
           write(out_unitp,*) 'Initial Matrix of dnS'
           CALL Write_MatOFdnS(MatdnS,nderiv)
         END IF

        DO it=1,max_it
          conv      = .FALSE.
          max_matpq = ZERO
          max_diag  = ZERO
          max_ip    = 0
          max_iq    = 0
          DO ip=1,N
            IF (abs(MatdnS(ip,ip)%d0-MatdnS(1,1)%d0) > max_diag)        &
                         max_diag = abs(MatdnS(ip,ip)%d0-MatdnS(1,1)%d0)
            DO iq=ip+1,N
              CALL Calc_MaxVal_OF_dnS(MatdnS(ip,iq),max_dnSmatpq)
              IF (max_dnSmatpq > max_matpq) THEN
                max_matpq = max_dnSmatpq
                max_ip    = ip
                max_iq    = iq
              END IF
            END DO
          END DO
          IF (debug) write(out_unitp,*) 'max_matpq,max_diag',max_ip,max_iq,max_matpq,max_diag
          IF (max_diag < tresh) max_diag = ONE
          !max_diag = max(max_diag,max_matpq)
          conv = (max_matpq/max_diag) < tresh
          IF (conv) EXIT

          ip = max_ip
          iq = max_iq

          CALL Calc_dnCos_AND_dnSin(dnCos,dnSin,                        &
                                    MatdnS(ip,ip),                      &
                                    MatdnS(iq,iq),                      &
                                    MatdnS(ip,iq),                      &
                                    dnTh,dnWork1,dnWork2,dnWork3,type_cs)
          IF (debug) THEN
            write(out_unitp,*) 'max_matpq',max_ip,max_iq,max_matpq
            write(out_unitp,*) 'th',dnTh%d0
            CALL Write_dnS(dnSin)
          END IF

          ! rotational matrix
          CALL sub_Id_TO_MatOFdnS(MatPdnS)

          !MatPdnS(max_ip,max_ip)%d0 =  cth
          CALL sub_dnS1_TO_dnS2(dnCos,MatPdnS(max_ip,max_ip),nderiv)
          !MatPdnS(max_iq,max_iq)%d0 =  cth
          CALL sub_dnS1_TO_dnS2(dnCos,MatPdnS(max_iq,max_iq),nderiv)
          !MatPdnS(max_ip,max_iq)%d0 = -sth
          CALL sub_dnS1_PROD_w_TO_dnS2(dnSin,-ONE,MatPdnS(max_ip,max_iq),nderiv)
          !MatPdnS(max_iq,max_ip)%d0 =  sth
          CALL sub_dnS1_TO_dnS2(dnSin,MatPdnS(max_iq,max_ip),nderiv)

          IF (debug) THEN
            !write(out_unitp,*) 'Rotational matrix',it
            !CALL Write_MatOFdnS(MatPdnS,nderiv=0)
          END IF

          ! eigenvectors
          CALL Mat1OFdnS_MUL_Mat2OFdnS_TO_Mat3OFdnS(EigVecdnS,MatPdnS,MatWorkdnS)
          CALL sub_Mat1OFdnS_TO_Mat2OFdnS(MatWorkdnS,EigVecdnS,nderiv)

          ! matrix
          CALL Mat1OFdnS_MUL_Mat2OFdnS_TO_Mat3OFdnS(MatdnS,MatPdnS,MatWorkdnS)
          CALL sub_Mat1OFdnS_TO_Mat2OFdnS(MatWorkdnS,MatdnS,nderiv)

          !MatPdnS(max_iq,max_ip)%d0 = -sth
          CALL sub_dnS1_PROD_w_TO_dnS2(dnSin,-ONE,MatPdnS(max_iq,max_ip),nderiv)
          !MatPdnS(max_ip,max_iq)%d0 =  sth
          CALL sub_dnS1_TO_dnS2(dnSin,MatPdnS(max_ip,max_iq),nderiv)

          CALL Mat1OFdnS_MUL_Mat2OFdnS_TO_Mat3OFdnS(MatPdnS,MatdnS,MatWorkdnS)
          CALL sub_Mat1OFdnS_TO_Mat2OFdnS(MatWorkdnS,MatdnS,nderiv)

          IF (debug) THEN
            !write(out_unitp,*) 'Intermediate Matrix of dnS',it
            !CALL Write_MatOFdnS(MatdnS,nderiv=0)
          END IF

        END DO

        IF (.NOT. conv) THEN
          write(out_unitp,*) 'diagonal Matrix of dnS (not converged)'
          CALL Write_MatOFdnS(MatdnS,nderiv)
          write(out_unitp,*) 'Eigenvector (in line) Matrix of dnS (not converged)'
          CALL Write_MatOFdnS(EigVecdnS,nderiv)
          write(out_unitp,*) ' ERROR in ',name_sub
          write(6,*) max_it,' iterations should never happen'
          STOP
        END IF


        ! deallocation
        CALL dealloc_dnS(dnTh)
        CALL dealloc_dnS(dnCos)
        CALL dealloc_dnS(dnSin)
        CALL dealloc_dnS(dnWork1)
        CALL dealloc_dnS(dnWork2)
        CALL dealloc_dnS(dnWork3)
        CALL dealloc_MatOFdnS(MatPdnS)
        CALL dealloc_MatOFdnS(MatWorkdnS)
        CALL dealloc_array(MatPdnS,'MatPdnS',name_sub)
        CALL dealloc_array(MatWorkdnS,'MatWorkdnS',name_sub)

      END SUBROUTINE DIAG01_MatOFdnS

      SUBROUTINE SORT_EigVectMat_OF_dnS(EigMatdnS,EigVecdnS,nderiv)
        TYPE (Type_dnS)     :: EigMatdnS(:,:),EigVecdnS(:,:)
        integer, optional   :: nderiv

        integer             :: i,j,k,N,nderiv_loc
        TYPE (Type_dnS)     :: dnWork

!------ for debuging -------------------------------------------
        integer :: err_mem,memory
        logical, parameter :: debug=.FALSE.
!        logical, parameter :: debug=.TRUE.
        character (len=*), parameter :: name_sub='SORT_EigVectMat_OF_dnS'
!---------------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
        END IF

        nderiv_loc = min(EigMatdnS(1,1)%nderiv,EigVecdnS(1,1)%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)


        N = ubound(EigMatdnS,dim=2)
        DO i=1,N
        DO j=i+1,N
          IF (EigMatdnS(i,i)%d0 > EigMatdnS(j,j)%d0) THEN
            IF (debug) write(out_unitp,*) 'Switch:',i,j
            !permutation

            CALL sub_dnS1_TO_dnS2(EigMatdnS(i,i),dnWork,nderiv_loc)
            CALL sub_dnS1_TO_dnS2(EigMatdnS(j,j),EigMatdnS(i,i),nderiv_loc)
            CALL sub_dnS1_TO_dnS2(dnWork,EigMatdnS(j,j),nderiv_loc)

            DO k=1,N
              CALL sub_dnS1_TO_dnS2(EigVecdnS(k,i),dnWork,nderiv_loc)
              CALL sub_dnS1_TO_dnS2(EigVecdnS(k,j),EigVecdnS(k,i),nderiv_loc)
              CALL sub_dnS1_TO_dnS2(dnWork,EigVecdnS(k,j),nderiv_loc)
            END DO

          END IF
        END DO
        END DO

        CALL dealloc_dnS(dnWork)
        IF (debug) THEN
          write(out_unitp,*) 'END ',name_sub
        END IF

      END SUBROUTINE SORT_EigVectMat_OF_dnS

      SUBROUTINE PHASE_EigVectMat_OF_dnS(EigVecdnS,nderiv)
        TYPE (Type_dnS)     :: EigVecdnS(:,:)
        integer, optional   :: nderiv

        integer             :: N,nderiv_loc
        real(kind=Rkind)    :: max_val
        integer             :: i,j,max_j

!------ for debuging -------------------------------------------
        integer :: err_mem,memory
        logical, parameter :: debug=.FALSE.
!        logical, parameter :: debug=.TRUE.
        character (len=*), parameter :: name_sub='PHASE_EigVectMat_OF_dnS'
!---------------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
        END IF

        nderiv_loc = EigVecdnS(1,1)%nderiv
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)


        N = ubound(EigVecdnS,dim=2)
        DO i=1,N
          max_val = abs(EigVecdnS(1,i)%d0)
          max_j   = 1
          DO j=2,N
            IF ( abs(EigVecdnS(j,i)%d0) > max_val ) THEN
              max_val = abs(EigVecdnS(j,i)%d0)
              max_j   = j
            END IF
          END DO
          IF (EigVecdnS(max_j,i)%d0 < ZERO) THEN
            DO j=1,N
              CALL sub_dnS1_PROD_w_TO_dnS2(EigVecdnS(j,i),-ONE,EigVecdnS(j,i),nderiv_loc)
            END DO
          END IF
        END DO

        IF (debug) THEN
          write(out_unitp,*) 'END ',name_sub
        END IF

      END SUBROUTINE PHASE_EigVectMat_OF_dnS

      SUBROUTINE Calc_MaxVal_OF_dnS(dnS,MaxValdnS)
        TYPE (Type_dnS), intent(in) :: dnS
        real (kind=Rkind) :: MaxValdnS

        MaxValdnS = abs(dnS%d0)
        IF (dnS%nderiv > 0) MaxValdnS = max(MaxValdnS,maxval(abs(dnS%d1)))
        IF (dnS%nderiv > 1) MaxValdnS = max(MaxValdnS,maxval(abs(dnS%d2)))
        IF (dnS%nderiv > 2) MaxValdnS = max(MaxValdnS,maxval(abs(dnS%d3)))

      END SUBROUTINE Calc_MaxVal_OF_dnS

      SUBROUTINE Calc_dnCos_AND_dnSin(dnCos,dnSin,dnMatpp,dnMatqq,dnMatpq, &
                                      dnTh,dnWork1,dnWork2,dnWork3,type_cs)
        TYPE (Type_dnS), intent(in)    :: dnMatpp,dnMatqq,dnMatpq
        TYPE (Type_dnS), intent(inout) :: dnCos,dnSin
        integer, intent(in)            :: type_cs


        TYPE (Type_dnS)     :: dnTh,dnWork1,dnWork2,dnWork3
        real (kind=Rkind)   :: delta,c,s,norm

!----- for debuging --------------------------------------------------
        integer :: err_mem,memory
        logical, parameter :: debug=.FALSE.
!        logical, parameter :: debug=.TRUE.
        character (len=*), parameter :: name_sub='Calc_dnCos_AND_dnSin'
!-----------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
        END IF


        SELECT CASE (type_cs)
        CASE (2)
          CALL CalcType2_dnCos_AND_dnSin(dnCos,dnSin,                   &
                                       dnMatpp,dnMatqq,dnMatpq,         &
                                       dnTh,dnWork1,dnWork2,dnWork3)

        CASE (1) ! direct calculation of cos and sin
          CALL CalcType1_dnCos_AND_dnSin(dnCos,dnSin,                   &
                                       dnMatpp,dnMatqq,dnMatpq,         &
                                       dnTh,dnWork1,dnWork2,dnWork3)
        CASE (0)
          IF (abs(dnMatpq%d0/(dnMatpp%d0-dnMatqq%d0)) < ONE) THEN
            CALL Calc2_dnCos_AND_dnSin(dnCos,dnSin,                     &
                                       dnMatpp,dnMatqq,dnMatpq,         &
                                       dnTh,dnWork1,dnWork2,dnWork3)
          ELSE
            CALL Calc1_dnCos_AND_dnSin(dnCos,dnSin,                     &
                                       dnMatpp,dnMatqq,dnMatpq,         &
                                       dnTh,dnWork1,dnWork2,dnWork3)
          END IF
        CASE DEFAULT
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) 'This option is not possible:',type_cs
          STOP
        END SELECT

        !write(out_unitp,*) 'for tested',dnTh%d0,(abs(dnTh%d0)<pi/FOUR)
        !write(out_unitp,*) 'cos(th)'
        !CALL Write_dnS(dnCos,nderiv=1)
        !write(out_unitp,*) 'sin(th)'
        !CALL Write_dnS(dnSin,nderiv=1)

        IF (debug) THEN
          write(out_unitp,*) 'type of cos/sin, type_cs:',type_cs
          write(out_unitp,*) 'cos(th),sin(th)',dnCos%d0,dnSin%d0
          write(out_unitp,*) 'END ',name_sub
        END IF

      END SUBROUTINE Calc_dnCos_AND_dnSin

      SUBROUTINE CalcType2_dnCos_AND_dnSin(dnCos,dnSin,dnMatpp,dnMatqq,dnMatpq, &
                                         dnTh,dnWork1,dnWork2,dnWork3)

        TYPE (Type_dnS), intent(in) :: dnMatpp,dnMatqq,dnMatpq
        TYPE (Type_dnS), intent(inout) :: dnCos,dnSin

        TYPE (Type_dnS)     :: dnTh,dnWork1,dnWork2,dnWork3


        TYPE (Type_dnS)     :: dnDelta,dnBigDelta
        TYPE (Type_dnS)     :: dnC1,dnS1,dnNorm1,dnC2,dnS2,dnNorm2


        real (kind=Rkind) :: norm1,norm2


!----- for debuging --------------------------------------------------
        integer :: err_mem,memory
        logical, parameter :: debug=.FALSE.
        !logical, parameter :: debug=.TRUE.
        character (len=*), parameter :: name_sub='CalcType2_dnCos_AND_dnSin'
!-----------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
        END IF

        CALL check_alloc_dnS(dnMatpp,'dnMatpp',name_sub)
        CALL check_alloc_dnS(dnMatqq,'dnMatqq',name_sub)
        CALL check_alloc_dnS(dnMatpq,'dnMatqp',name_sub)


        ! here we want: -pi/4 < theta < pi/4

        !delta = (dnMatpp%d0-dnMatqq%d0)*HALF
        ! delta = (Mat(p,p)-Mat(q,q))/2
        CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnMatpp,HALF,dnMatqq,-HALF,dnDelta)


        IF (abs(dnMatpq%d0) <= abs(dnDelta%d0) ) THEN
          ! we use 2theta = atan( Mat(p,q)/delta )   => -pi/4 < theta < pi/4

          ! 1/delta => dnWork1
          CALL sub_dnS1_TO_dntR2(dnDelta,dnWork1,90)
          ! Mat(p,q)/delta => dnWork2
          CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnMatpq,dnWork1,dnWork2)

          ! atan(Mat(p,q)/delta)=atan(dnWork2) => dnWork1  (2*th)
          CALL sub_dnS1_TO_dntR2(dnWork2,dnWork1,70)
          !write(6,*) '2th',dnWork1%d0

          ! theta
          CALL sub_dnS1_PROD_w_TO_dnS2(dnWork1,HALF,dnTh)


        ELSE ! use the relation between atan(x)+atan(1/x) = pi/2 (x>0) or -pi/2 (x<0)
          ! 1/Mat(p,q) => dnWork1
          CALL sub_dnS1_TO_dntR2(dnMatpq,dnWork1,90)
          ! delta/Mat(p,q) => dnWork2
          CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnDelta,dnWork1,dnWork2)

          ! atan(delta/Mat(p,q))=atan(dnWork2) => dnWork1  (2*thbis)
          CALL sub_dnS1_TO_dntR2(dnWork2,dnWork1,70)

          !write(6,*) '2thbis',dnWork1%d0


          IF (dnWork1%d0 > 0) THEN
            dnWork1%d0 = dnWork1%d0 - PI/TWO
          ELSE
            dnWork1%d0 = dnWork1%d0 + PI/TWO
          END IF

          ! theta
          CALL sub_dnS1_PROD_w_TO_dnS2(dnWork1,-HALF,dnTh)


        END IF

        ! dnTh => cos(dnTh)
        CALL sub_dnS1_TO_dntR2(dnTh,dnCos,2)
        ! dnTh => sin(dnTh)
        CALL sub_dnS1_TO_dntR2(dnTh,dnSin,3)


        !write(6,*) 'type_cs=2',dnTh%d0, (abs(dnTh%d0)< PI/FOUR)


        CALL dealloc_dnS(dnDelta)


        IF (debug) THEN
          write(out_unitp,*) 'for tested',dnTh%d0
          write(out_unitp,*) 'New cos(th)'
          CALL Write_dnS(dnCos)
          write(out_unitp,*) 'New sin(th)'
          CALL Write_dnS(dnSin)
          write(out_unitp,*) 'END ',name_sub
        END IF

      END SUBROUTINE CalcType2_dnCos_AND_dnSin

      SUBROUTINE CalcType1_dnCos_AND_dnSin(dnCos,dnSin,dnMatpp,dnMatqq,dnMatpq, &
                                           dnTh,dnWork1,dnWork2,dnWork3)

        TYPE (Type_dnS), intent(in) :: dnMatpp,dnMatqq,dnMatpq
        TYPE (Type_dnS), intent(inout) :: dnCos,dnSin

        TYPE (Type_dnS)     :: dnTh,dnWork1,dnWork2,dnWork3


        TYPE (Type_dnS)     :: dnDelta,dnBigDelta
        TYPE (Type_dnS)     :: dnC1,dnS1,dnNorm1,dnC2,dnS2,dnNorm2


        real (kind=Rkind) :: norm1,norm2


!----- for debuging --------------------------------------------------
        integer :: err_mem,memory
        logical, parameter :: debug=.FALSE.
        !logical, parameter :: debug=.TRUE.
        character (len=*), parameter :: name_sub='CalcType1_dnCos_AND_dnSin'
!-----------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
        END IF

        CALL check_alloc_dnS(dnMatpp,'dnMatpp',name_sub)
        CALL check_alloc_dnS(dnMatqq,'dnMatqq',name_sub)
        CALL check_alloc_dnS(dnMatpq,'dnMatqp',name_sub)

        ! new formula
        ! delta = (dnMatpp%d0-dnMatqq%d0)*HALF
        ! delta = (Mat(p,p)-Mat(q,q))/2
        CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnMatpp,HALF,dnMatqq,-HALF,dnDelta)

        !BigDelta = sqrt(delta**2+dnMatpq%d0**2)
        CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnDelta,dnDelta,dnWork1)
        CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnMatpq,dnMatpq,dnWork2)
        CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnWork1,ONE,dnWork2,ONE,dnWork3)
        CALL sub_dnS1_TO_dntR2(dnWork3,dnBigDelta,91) ! sqrt

        !Cos and sin from the first formula
        ! c = delta+sqrt(delta**2+dnMatpq%d0**2) = delta + BigDelta
        ! s = dnMatpq%d0
        CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnDelta,ONE,dnBigDelta,ONE,dnC1)
        CALL sub_dnS1_TO_dnS2(dnMatpq,dnS1)
        ! norm : Cos^2 + Sin^2
        CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnC1,dnC1,dnWork1)
        CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnS1,dnS1,dnWork2)
        CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnWork1,ONE,dnWork2,ONE,dnNorm1)
        CALL Calc_MaxVal_OF_dnS(dnNorm1,norm1)
        !write(6,*) 'dnNorm1',norm1

        !Cos and sin from the second formula
        ! s = -delta+sqrt(delta**2+dnMatpq%d0**2) = -delta + BigDelta
        ! c = dnMatpq%d0
        CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnDelta,ONE,dnBigDelta,-ONE,dnC2)
        CALL sub_dnS1_TO_dnS2(dnMatpq,dnS2)
        ! norm : Cos^2 + Sin^2
        CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnC2,dnC2,dnWork1)
        CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnS2,dnS2,dnWork2)
        CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnWork1,ONE,dnWork2,ONE,dnNorm2)
        CALL Calc_MaxVal_OF_dnS(dnNorm2,norm2)
        !write(6,*) 'dnNorm2',norm2
        IF (abs(dnNorm1%d0) < ONETENTH**10 .AND. abs(dnNorm2%d0) < ONETENTH**10) THEN
          write(out_unitp,*) 'Warning: degenerated matrix'
        END IF
        !write(6,*) 'type_cs=1, norm1,norm2',norm1,norm2

        IF (norm1 >= norm2) THEN
          CALL sub_dnS1_TO_dntR2(dnNorm1,dnWork1,91) ! sqrt
          CALL sub_dnS1_TO_dntR2(dnWork1,dnNorm1,90) ! 1/x

          CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnC1,dnNorm1,dnCos)
          CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnS1,dnNorm1,dnSin)

          !write(6,*) 'type_cs=1, 1st formula',dnCos%d0,dnSin%d0

        ELSE
          CALL sub_dnS1_TO_dntR2(dnNorm2,dnWork2,91) ! sqrt
          CALL sub_dnS1_TO_dntR2(dnWork2,dnNorm2,90) ! 1/x

          CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnC2,dnNorm2,dnCos)
          CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnS2,dnNorm2,dnSin)

          !write(6,*) 'type_cs=1, 2d formula',dnCos%d0,dnSin%d0

        END IF
        CALL sub_dnS1_TO_dnS2(dnSin,dnTh)  ! to be able to test ...

        ! deallocation
        CALL dealloc_dnS(dnC1)
        CALL dealloc_dnS(dnS1)
        CALL dealloc_dnS(dnNorm1)
        CALL dealloc_dnS(dnC2)
        CALL dealloc_dnS(dnS2)
        CALL dealloc_dnS(dnNorm2)

        CALL dealloc_dnS(dnDelta)
        CALL dealloc_dnS(dnBigDelta)


        IF (debug) THEN
          write(out_unitp,*) 'for tested',dnTh%d0
          write(out_unitp,*) 'New cos(th)'
          CALL Write_dnS(dnCos)
          write(out_unitp,*) 'New sin(th)'
          CALL Write_dnS(dnSin)
          write(out_unitp,*) 'END ',name_sub
        END IF

      END SUBROUTINE CalcType1_dnCos_AND_dnSin


      SUBROUTINE Calc2_dnCos_AND_dnSin(dnCos,dnSin,dnMatpp,dnMatqq,dnMatpq, &
                                      dnTh,dnWork1,dnWork2,dnWork3)
        TYPE (Type_dnS), intent(in) :: dnMatpp,dnMatqq,dnMatpq
        TYPE (Type_dnS), intent(inout) :: dnCos,dnSin

        TYPE (Type_dnS)     :: dnTh,dnWork1,dnWork2,dnWork3


        real (kind=Rkind) :: th,x,y


!----- for debuging --------------------------------------------------
        integer :: err_mem,memory
        logical, parameter :: debug=.FALSE.
!        logical, parameter :: debug=.TRUE.
        character (len=*), parameter :: name_sub='Calc2_dnCos_AND_dnSin'
!-----------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
          write(out_unitp,*) 'calc atan2(y,x) from atan(y/x)'
        END IF

        CALL check_alloc_dnS(dnMatpp,'dnMatpp',name_sub)
        CALL check_alloc_dnS(dnMatqq,'dnMatqq',name_sub)
        CALL check_alloc_dnS(dnMatpq,'dnMatqp',name_sub)


        ! numerateur : dnWork2
        CALL sub_dnS1_PROD_w_TO_dnS2(dnMatpq,TWO,dnWork2)
        y = dnWork2%d0
        ! denominateur : dnWork3
        CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnMatpp,ONE,dnMatqq,-ONE,dnWork3)
        x = dnWork3%d0
        ! 1/dnWork3 => dnWork1
        CALL sub_dnS1_TO_dntR2(dnWork3,dnWork1,90)

        ! dnWork2*dnWork1 => dnWork3 = 2*Mat(p,q)/(Mat(p,p)-Mat(q,q))
        CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnWork2,dnWork1,dnWork3)

        ! atan(dnWork3) => dnWork1  (2*th')
        CALL sub_dnS1_TO_dntR2(dnWork3,dnWork1,70)

        IF (x > ZERO) THEN
          CONTINUE ! th is unchanged
        ELSE IF (x < ZERO) THEN
          IF (y < ZERO) THEN
            dnWork1%d0 = dnWork1%d0 - pi
          ELSE
            dnWork1%d0 = dnWork1%d0 + pi
          END IF
        ELSE
          STOP 'pb with x=0'
        END IF

        ! dnWork1*HALF => dnTh
        CALL sub_dnS1_PROD_w_TO_dnS2(dnWork1,HALF,dnTh)

        th = HALF * atan2(y,x)

        ! dnTh => cos(dnTh)
        CALL sub_dnS1_TO_dntR2(dnTh,dnCos,2)
        ! dnTh => sin(dnTh)
        CALL sub_dnS1_TO_dntR2(dnTh,dnSin,3)

        IF (debug) THEN
          write(out_unitp,*) 'y: num',y
          write(out_unitp,*) 'x: deno',x
          write(out_unitp,*) 'dnTh%d0,th',dnTh%d0,th
          write(out_unitp,*) 'cos(th),sin(th)',dnCos%d0,dnSin%d0
          write(out_unitp,*) 'END ',name_sub
        END IF

      END SUBROUTINE Calc2_dnCos_AND_dnSin
      SUBROUTINE Calc1_dnCos_AND_dnSin(dnCos,dnSin,dnMatpp,dnMatqq,dnMatpq, &
                                      dnTh,dnWork1,dnWork2,dnWork3)
        TYPE (Type_dnS), intent(in) :: dnMatpp,dnMatqq,dnMatpq
        TYPE (Type_dnS), intent(inout) :: dnCos,dnSin

        TYPE (Type_dnS)     :: dnTh,dnWork1,dnWork2,dnWork3


        real (kind=Rkind) :: th,x,y


!----- for debuging --------------------------------------------------
        integer :: err_mem,memory
        logical, parameter :: debug=.FALSE.
!        logical, parameter :: debug=.TRUE.
        character (len=*), parameter :: name_sub='Calc1_dnCos_AND_dnSin'
!-----------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
          write(out_unitp,*) 'calc atan2(y,x) from atan(x/y)'
        END IF

        CALL check_alloc_dnS(dnMatpp,'dnMatpp',name_sub)
        CALL check_alloc_dnS(dnMatqq,'dnMatqq',name_sub)
        CALL check_alloc_dnS(dnMatpq,'dnMatqp',name_sub)


        ! calc atan2(y,x) from atan(x/y)
        ! Here y cannot be zero (y is proportional to the off diagonal term)

        ! y => dnWork2
        CALL sub_dnS1_PROD_w_TO_dnS2(dnMatpq,TWO,dnWork2)
        y = dnWork2%d0
        ! 1/dnWork2 => dnWork1
        CALL sub_dnS1_TO_dntR2(dnWork2,dnWork1,90)

        ! x => dnWork3
        CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnMatpp,ONE,dnMatqq,-ONE,dnWork3)
        x = dnWork3%d0


        ! dnWork3*dnWork1 => dnWork2
        CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnWork3,dnWork1,dnWork2)

        ! atan(dnWork2) => dnWork1  (2*th')
        CALL sub_dnS1_TO_dntR2(dnWork2,dnWork1,70)

        IF (y > ZERO) THEN
          dnWork1%d0 = dnWork1%d0 - pi*HALF
        ELSE IF (y < ZERO) THEN
          dnWork1%d0 = dnWork1%d0 + pi*HALF
        ELSE ! y=0
            STOP 'pb with y=0'
        END IF


        ! dnWork1*HALF => dnTh
        CALL sub_dnS1_PROD_w_TO_dnS2(dnWork1,-HALF,dnTh)

        th = HALF * atan2(y,x)

        ! dnTh => cos(dnTh)
        CALL sub_dnS1_TO_dntR2(dnTh,dnCos,2)
        ! dnTh => sin(dnTh)
        CALL sub_dnS1_TO_dntR2(dnTh,dnSin,3)

        IF (debug) THEN
          write(out_unitp,*) 'y: deno',y
          write(out_unitp,*) 'x: num',x
          write(out_unitp,*) 'dnTh%d0,th',dnTh%d0,th
          write(out_unitp,*) 'cos(th),sin(th)',dnCos%d0,dnSin%d0
          write(out_unitp,*) 'END ',name_sub
        END IF

      END SUBROUTINE Calc1_dnCos_AND_dnSin

      SUBROUTINE Mat1OFdnS_MUL_Mat2OFdnS_TO_Mat3OFdnS(Mat1OFdnS,&
                                         Mat2OFdnS,Mat3OFdnS,nderiv)
        TYPE (Type_dnS)   :: Mat1OFdnS(:,:),Mat2OFdnS(:,:),Mat3OFdnS(:,:)
        integer, optional :: nderiv

        TYPE (Type_dnS)   :: dnWork
        integer :: i,j,k


        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='Mat1OFdnS_MUL_Mat2OFdnS_TO_Mat3OFdnS'

        CALL check_alloc_MatOFdnS(Mat1OFdnS,'Mat1OFdnS',name_sub)
        CALL check_alloc_MatOFdnS(Mat2OFdnS,'Mat2OFdnS',name_sub)
        CALL check_alloc_MatOFdnS(Mat3OFdnS,'Mat3OFdnS',name_sub)

        nderiv_loc = min(Mat1OFdnS(1,1)%nderiv,Mat2OFdnS(1,1)%nderiv,Mat3OFdnS(1,1)%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (Mat1OFdnS(1,1)%nb_var_deriv /= Mat2OFdnS(1,1)%nb_var_deriv .OR. &
            Mat1OFdnS(1,1)%nb_var_deriv /= Mat3OFdnS(1,1)%nb_var_deriv) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nb_var_deriv in Mat1OFdnS or Mat2OFdnS or Mat3OFdnS are different!',&
                    Mat1OFdnS(1,1)%nb_var_deriv,Mat2OFdnS(1,1)%nb_var_deriv,Mat3OFdnS(1,1)%nb_var_deriv
          STOP
        END IF

        CALL alloc_dnS(dnWork,Mat3OFdnS(1,1)%nb_var_deriv,nderiv_loc)


        !  MatdnS3(i,j) = sum_k MatdnS1(i,k) . MatdnS2(k,j)
        DO i=lbound(Mat3OFdnS,dim=1),ubound(Mat3OFdnS,dim=1)
        DO j=lbound(Mat3OFdnS,dim=2),ubound(Mat3OFdnS,dim=2)
          CALL sub_ZERO_TO_dnS(Mat3OFdnS(i,j),nderiv_loc)
          DO k=lbound(Mat1OFdnS,dim=2),ubound(Mat1OFdnS,dim=2)
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(Mat1OFdnS(i,k),Mat2OFdnS(k,j),dnWork,nderiv_loc)
            CALL sub_dnS1_PLUS_dnS2_TO_dnS3(Mat3OFdnS(i,j),dnWork,Mat3OFdnS(i,j),nderiv)

          END DO
        END DO
        END DO

        CALL dealloc_dnS(dnWork)

      END SUBROUTINE Mat1OFdnS_MUL_Mat2OFdnS_TO_Mat3OFdnS

      SUBROUTINE Mat1OFdnS_wPLUS_Mat2OFdnS_TO_Mat3OFdnS(Mat1OFdnS,w1,&
                                         Mat2OFdnS,w2,Mat3OFdnS,nderiv)
        TYPE (Type_dnS)   :: Mat1OFdnS(:,:),Mat2OFdnS(:,:),Mat3OFdnS(:,:)
        real (kind=Rkind) :: w1,w2
        integer, optional :: nderiv

        integer :: i,j


        integer :: nderiv_loc
        character (len=*), parameter :: name_sub='Mat1OFdnS_wPLUS_Mat2OFdnS_TO_Mat3OFdnS'

        CALL check_alloc_MatOFdnS(Mat1OFdnS,'Mat1OFdnS',name_sub)
        CALL check_alloc_MatOFdnS(Mat2OFdnS,'Mat2OFdnS',name_sub)
        CALL check_alloc_MatOFdnS(Mat3OFdnS,'Mat3OFdnS',name_sub)

        nderiv_loc = min(Mat1OFdnS(1,1)%nderiv,Mat2OFdnS(1,1)%nderiv,Mat3OFdnS(1,1)%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (Mat1OFdnS(1,1)%nb_var_deriv /= Mat2OFdnS(1,1)%nb_var_deriv .OR. &
            Mat1OFdnS(1,1)%nb_var_deriv /= Mat3OFdnS(1,1)%nb_var_deriv) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nb_var_deriv in Mat1OFdnS or Mat2OFdnS or Mat3OFdnS are different!',&
                    Mat1OFdnS(1,1)%nb_var_deriv,Mat2OFdnS(1,1)%nb_var_deriv,Mat3OFdnS(1,1)%nb_var_deriv
          STOP
        END IF


        !  MatdnS3(i,j) = sum_k MatdnS1(i,k) . MatdnS2(k,j)
        DO i=lbound(Mat3OFdnS,dim=1),ubound(Mat3OFdnS,dim=1)
        DO j=lbound(Mat3OFdnS,dim=2),ubound(Mat3OFdnS,dim=2)
          CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(Mat1OFdnS(i,j),w1,           &
                            Mat2OFdnS(i,j),w2,Mat3OFdnS(i,j),nderiv_loc)
        END DO
        END DO

      END SUBROUTINE Mat1OFdnS_wPLUS_Mat2OFdnS_TO_Mat3OFdnS

      SUBROUTINE Mat1OFdnS_MUL_Vec2OFdnS_TO_Vec3OFdnS(Mat1OFdnS,&
                                          Vec2OFdnS,Vec3OFdnS,nderiv)
        TYPE (Type_dnS)   :: Mat1OFdnS(:,:),Vec2OFdnS(:),Vec3OFdnS(:)
        integer, optional :: nderiv

        TYPE (Type_dnS)   :: dnWork
        integer :: i,k


        integer :: nderiv_loc,nb_var_deriv_loc
        character (len=*), parameter ::                                 &
                         name_sub='Mat1OFdnS_MUL_Vec2OFdnS_TO_Vec3OFdnS'

        nderiv_loc = min(minval(Mat1OFdnS%nderiv),                      &
                  minval(Vec2OFdnS%nderiv),minval(Vec3OFdnS%nderiv))
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (minval(Mat1OFdnS%nb_var_deriv) /= minval(Vec2OFdnS%nb_var_deriv) &
       .OR. minval(Mat1OFdnS%nb_var_deriv) /= minval(Vec3OFdnS%nb_var_deriv)) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) ' nb_var_deriv in Mat1OFdnS or Vec2OFdnS',  &
                            ' or Vec3OFdnS are different!',             &
                             minval(Mat1OFdnS%nb_var_deriv),            &
                             minval(Vec2OFdnS%nb_var_deriv),            &
                             minval(Vec3OFdnS%nb_var_deriv)
          STOP
        END IF
        nb_var_deriv_loc = minval(Vec2OFdnS%nb_var_deriv)

        CALL alloc_dnS(dnWork,nb_var_deriv_loc,nderiv_loc)


        ! Vec3OFdnS(i) = sum_k Mat1OFdnS(i,k) . Vec2OFdnS(k)
        DO i=lbound(Vec3OFdnS,dim=1),ubound(Vec3OFdnS,dim=1)
          CALL sub_ZERO_TO_dnS(Vec3OFdnS(i),nderiv_loc)
          DO k=lbound(Mat1OFdnS,dim=2),ubound(Mat1OFdnS,dim=2)
            CALL sub_dnS1_PROD_dnS2_TO_dnS3(Mat1OFdnS(i,k),Vec2OFdnS(k),dnWork,nderiv_loc)
            CALL sub_dnS1_PLUS_dnS2_TO_dnS3(Vec3OFdnS(i),dnWork,Vec3OFdnS(i),nderiv)

          END DO
        END DO

        CALL dealloc_dnS(dnWork)

      END SUBROUTINE Mat1OFdnS_MUL_Vec2OFdnS_TO_Vec3OFdnS

      SUBROUTINE MatOFdnS_TO_VecOFdnS(MatOFdnS,VecOFdnS,ilc,line,nderiv)
        TYPE (Type_dnS), intent(in)    :: MatOFdnS(:,:)
        TYPE (Type_dnS), intent(inout) :: vecOFdnS(:)
        integer, intent(in)            :: ilc
        logical, intent(in)            :: line

        integer, optional :: nderiv


        integer :: i
        integer :: nderiv_loc

        character (len=*), parameter :: name_sub='MatOFdnS_TO_VecOFdnS'

        CALL check_alloc_MatOFdnS(MatOFdnS,'MatOFdnS',name_sub)
        CALL check_alloc_VecOFdnS(VecOFdnS,'VecOFdnS',name_sub)

        nderiv_loc = min(MatOFdnS(1,1)%nderiv,VecOFdnS(1)%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        IF (MatOFdnS(1,1)%nb_var_deriv /= VecOFdnS(1)%nb_var_deriv) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nb_var_deriv in MatOFdnS and VecOFdnS are different!',&
                    MatOFdnS(1,1)%nb_var_deriv,VecOFdnS(1)%nb_var_deriv
          STOP
        END IF


        IF (line) THEN
          DO i=lbound(MatOFdnS,dim=1),ubound(MatOFdnS,dim=1)
            CALL sub_dnS1_TO_dnS2(MatOFdnS(ilc,i),VecOFdnS(i),nderiv_loc)
          END DO
        ELSE
          DO i=lbound(MatOFdnS,dim=1),ubound(MatOFdnS,dim=1)
            CALL sub_dnS1_TO_dnS2(MatOFdnS(i,ilc),VecOFdnS(i),nderiv_loc)
          END DO
        END IF

        END SUBROUTINE MatOFdnS_TO_VecOFdnS
!================================================================
!
!     dnMat = 0
!
!================================================================

      SUBROUTINE sub_ZERO_TO_MatOFdnS(MatOFdnS,nderiv)
        TYPE (Type_dnS) :: MatOFdnS(:,:)
        integer, optional :: nderiv
        integer :: nderiv_loc,i,j

        nderiv_loc = minval(MatOFdnS%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        DO i=lbound(MatOFdnS,dim=1),ubound(MatOFdnS,dim=1)
        DO j=lbound(MatOFdnS,dim=2),ubound(MatOFdnS,dim=2)
           CALL sub_ZERO_TO_dnS(MatOFdnS(i,j),nderiv=nderiv_loc)
        END DO
        END DO

      END SUBROUTINE sub_ZERO_TO_MatOFdnS

      SUBROUTINE sub_Weight_MatOFdnS(MatOFdnS,w,nderiv)
        TYPE (Type_dnS) :: MatOFdnS(:,:)
        real (kind=Rkind) :: w
        integer, optional :: nderiv

        integer :: nderiv_loc,i,j

        nderiv_loc = minval(MatOFdnS%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        DO j=lbound(MatOFdnS,dim=2),ubound(MatOFdnS,dim=2)
        DO i=lbound(MatOFdnS,dim=1),ubound(MatOFdnS,dim=1)
           CALL sub_Weight_dnS(MatOFdnS(i,j),w,nderiv=nderiv_loc)
        END DO
        END DO

      END SUBROUTINE sub_Weight_MatOFdnS


      SUBROUTINE sub_Id_TO_MatOFdnS(MatOFdnS,nderiv)
        TYPE (Type_dnS) :: MatOFdnS(:,:)
        integer, optional :: nderiv
        integer :: nderiv_loc,i

        nderiv_loc = minval(MatOFdnS%nderiv)
        IF (present(nderiv)) nderiv_loc = min(nderiv_loc,nderiv)

        CALL sub_ZERO_TO_MatOFdnS(MatOFdnS,nderiv=nderiv_loc)

        DO i=lbound(MatOFdnS,dim=1),ubound(MatOFdnS,dim=1)
           MatOFdnS(i,i)%d0 = ONE
        END DO

      END SUBROUTINE sub_Id_TO_MatOFdnS

      END MODULE mod_MatOFdnS

