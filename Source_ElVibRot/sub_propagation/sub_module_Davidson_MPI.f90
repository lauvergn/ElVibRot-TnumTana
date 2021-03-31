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
!         Université de Montpellier, France
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
MODULE mod_Davidson_MPI
  USE mod_Constant
  USE mod_MPI 
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: Schmidt_process_MPI
  PUBLIC :: MakeResidual_Davidson_MPI3,MakeResidual_Davidson_MPI4
  PUBLIC :: MakeResidual_Davidson_j_MPI3
  PUBLIC :: exit_Davidson_external_MPI

  CONTAINS

!=======================================================================================
!> @brief Schmidt process with MPI
!> <psi|psi> calculated in the meantime 
!
! note, new psi(ndim+1) are distributed to the other threads here
!=======================================================================================
  SUBROUTINE Schmidt_process_MPI(S_Overlap1D,psi,ndim,isym,With_Grid) 
    USE mod_system
    USE mod_ana_psi_MPI,ONLY:norm_psi_MPI
    USE mod_psi,        ONLY:param_psi
    USE mod_psi_Op_MPI, ONLY:Set_symab_OF_psiBasisRep_MPI,calculate_overlap_MPI,       &
                             calculate_overlap1D_MPI,distribute_psi_MPI
    USE mod_propa,      ONLY:param_Davidson
    USE mod_MPI_aux
    IMPLICIT NONE

    Real(kind=Rkind),intent(inout)             :: S_Overlap1D(:)
    TYPE(param_psi), intent(inout)             :: psi(:)
    Integer,         intent(in)                :: ndim
    Integer,         intent(in)                :: isym
    Logical,optional,intent(in)                :: With_Grid
    
    Real(kind=Rkind)                           :: RS,RS2
    Integer                                    :: twice
    Integer                                    :: ii
    Integer                                    :: d1
    Integer                                    :: d2
    Logical                                    :: With_Grid_loc
    Character(len=*),parameter                 :: name_sub='Schmidt_process_MPI'

#if(run_MPI)

    With_Grid_loc=.FALSE.
    IF(present(With_Grid)) With_Grid_loc=With_Grid
    
    IF(MPI_scheme/=1) CALL distribute_psi_MPI(psi,ndim+1,ndim+1,With_Grid_loc)

    d1=bounds_MPI(1,MPI_id)
    d2=bounds_MPI(2,MPI_id)

    DO twice=0,1
      IF(With_Grid_loc) THEN
        ! normalization
        ! nb_per_MPI=psi(ndim+1)%nb_qaie/MPI_np
        RS=real(dot_product(psi(ndim+1)%RvecG(d1:d2),psi(ndim+1)%RvecG(d1:d2)))
        CALL MPI_Reduce_sum_Bcast(RS)
        psi(ndim+1)%RvecG=psi(ndim+1)%RvecG/sqrt(RS)

        CALL calculate_overlap1D_MPI(psi,ndim+1,With_Grid=With_Grid_loc,               &
                                     S_overlap1D=S_overlap1D)

        DO ii=1,ndim
          RS=S_overlap1D(ii)
          IF(RS==ZERO) CYCLE
          psi(ndim+1)%RvecG=(psi(ndim+1)%RvecG-psi(ii)%RvecG*RS)/sqrt(ONE-RS**2)
        ENDDO ! ii=1,ndim

      ELSE !----------------------------------------------------------------------------

        ! normalization
        ! nb_per_MPI=psi(ndim+1)%nb_tot/MPI_np
        RS=dot_product(psi(ndim+1)%RvecB(d1:d2),psi(ndim+1)%RvecB(d1:d2))
        CALL MPI_Reduce_sum_Bcast(RS)
        psi(ndim+1)%RvecB=psi(ndim+1)%RvecB/sqrt(RS)

        ! calculate S_overlap_add
        CALL calculate_overlap1D_MPI(psi,ndim+1,With_Grid=With_Grid_loc,               &
                                     S_overlap1D=S_overlap1D)

        DO ii=1,ndim
          RS=S_overlap1D(ii)
          IF(RS==ZERO) CYCLE
          ! be careful on the difference of vec on different threads
          psi(ndim+1)%RvecB=psi(ndim+1)%RvecB-psi(ii)%RvecB*RS
        ENDDO ! ii=1,ndim
      ENDIF
    ENDDO ! do Schmidt process twice

#endif
  END SUBROUTINE Schmidt_process_MPI
!=======================================================================================


!=======================================================================================
! MPI for calculating residual in the main Davidson procedure
! V3
!=======================================================================================
  SUBROUTINE MakeResidual_Davidson_MPI3(ndim,g,psi,Hpsi,Ene,Vec,conv,converge,         &
                                        VecToBeIncluded,tab_norm2g,norm2g,convergeResi,&
                                        convergeEne,fresidu,iresidu,nb_diago,epsi)
    USE mod_system
    USE mod_psi,     ONLY : param_psi,norm2_psi,Set_symab_OF_psiBasisRep
    USE mod_MPI_aux
    IMPLICIT NONE

    TYPE(param_psi), intent(inout)              :: g
    TYPE(param_psi), intent(in)                 :: psi(:)
    TYPE(param_psi), intent(in)                 :: Hpsi(:)
    Real(kind=Rkind),intent(in)                 :: Ene(:)
    Real(kind=Rkind),intent(in)                 :: Vec(:,:)
    Real(kind=Rkind),intent(inout)              :: tab_norm2g(:)
    Real(kind=Rkind),intent(inout)              :: norm2g
    Real(kind=Rkind),intent(in)                 :: epsi
    Integer,         intent(in)                 :: ndim
    Integer,         intent(in)                 :: nb_diago
    Integer,         intent(inout)              :: fresidu
    Integer,         intent(inout)              :: iresidu
    Logical,         intent(in)                 :: VecToBeIncluded(:)
    Logical,         intent(in)                 :: convergeEne(:)
    Logical,         intent(inout)              :: convergeResi(:)
    Logical,         intent(inout)              :: converge(:)
    Logical,         intent(inout)              :: conv

    Integer                                     :: case_vec
    Integer                                     :: size_vec
    Integer                                     :: isym
    Integer                                     :: d1
    Integer                                     :: d2
    Integer                                     :: ii
    Integer                                     :: jj

#if(run_MPI)

    IF(MPI_id==0) THEN
      IF(allocated(g%RvecB)) THEN   
        case_vec=1 
        size_vec=size(g%RvecB)
      ELSEIF(allocated(g%CvecB)) THEN
        case_vec=2
        size_vec=size(g%CvecB)
      ELSEIF(allocated(g%RvecG)) THEN
        case_vec=3
        size_vec=size(g%RvecG)
      ELSEIF(allocated(g%CvecG)) THEN
        case_vec=4
        size_vec=size(g%CvecG)
      ELSE                       
        case_vec=0
        STOP 'ERROR in g%vec of MakeResidual_Davidson_MPI'
      ENDIF
    ENDIF
    CALL MPI_BCAST(case_vec,size1_MPI,Int_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
    CALL MPI_BCAST(size_vec,size1_MPI,Int_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)

    IF(.NOT. allocated(bounds_MPI)) allocate(bounds_MPI(1:2,0:MPI_np-1))
    nb_per_MPI=size_vec/MPI_np
    nb_rem_MPI=mod(size_vec,MPI_np) 
    DO i_MPI=0,MPI_np-1
      bounds_MPI(1,i_MPI)=i_MPI*nb_per_MPI+1+MIN(i_MPI,nb_rem_MPI)
      bounds_MPI(2,i_MPI)=(i_MPI+1)*nb_per_MPI+MIN(i_MPI,nb_rem_MPI)                   &
                                              +merge(1,0,nb_rem_MPI>i_MPI)
    ENDDO
    
    d1=bounds_MPI(1,MPI_id)
    d2=bounds_MPI(2,MPI_id)

  !  SELECT CASE (case_vec)
  !  CASE(1,3)
  !    IF(MPI_id==0) allocate(Rvec(size_vec))
  !    IF(MPI_id/=0) allocate(Rvec(bounds_MPI(1,MPI_id):bounds_MPI(2,MPI_id)))
  !  CASE(2,4)
  !    IF(MPI_id==0) allocate(Cvec(size_vec))
  !    IF(MPI_id/=0) allocate(Cvec(bounds_MPI(1,MPI_id):bounds_MPI(2,MPI_id)))
  !  END SELECT
    
    IF(MPI_id/=0) THEN
      SELECT CASE (case_vec)
      CASE(1)
        CALL allocate_array(g%RvecB,d1,d2)
      CASE(2)
        CALL allocate_array(g%CvecB,d1,d2)
      CASE(3)
        CALL allocate_array(g%RvecG,d1,d2)
      CASE(4)
        CALL allocate_array(g%CvecG,d1,d2)
      END SELECT
    ENDIF

    ! converge,VecToBeIncluded,fresidu,norm2g,epsi,convergeEne synchronized 
    !CALL MPI_BCAST(converge,ndim,MPI_LOGICAL,root_MPI,MPI_COMM_WORLD,MPI_err)
    DO jj=1,ndim
      IF (.NOT. converge(jj) .AND. VecToBeIncluded(jj)) THEN
        isym = maxloc(abs(Vec(:,jj)),dim=1) ! to find the rigth symmetry
        
        CALL Residual_Davidson_sum_MPI(g,Hpsi,psi,Vec,Ene,ndim,case_vec,size_vec,jj)

        IF(MPI_id==0) THEN
          CALL Set_symab_OF_psiBasisRep(g,symab=psi(isym)%symab)
          CALL norm2_psi(g)
          tab_norm2g(jj) = sqrt(g%norm2)
        ENDIF
        
      ENDIF ! for .NOT. converge(jj) .AND. VecToBeIncluded(jj)  
    ENDDO ! for jj=1,ndim

    CALL MPI_BCAST(tab_norm2g,ndim,MPI_Real8,root_MPI,MPI_COMM_WORLD,MPI_err)
    
    DO jj=1,ndim
      IF (.NOT. converge(jj) .AND. VecToBeIncluded(jj)) THEN
        IF(fresidu==0) fresidu=jj
        IF(tab_norm2g(jj)>norm2g) THEN
          iresidu=jj
          norm2g=tab_norm2g(iresidu)
        ENDIF
        convergeResi(jj)=tab_norm2g(jj)<epsi
      ENDIF ! for .NOT. converge(jj) .AND. VecToBeIncluded(jj)
      converge(jj) = (convergeEne(jj) .AND. convergeResi(jj))
    ENDDO ! for jj=1,ndim
    conv = all(converge(1:nb_diago))

    IF(MPI_id/=0) THEN
      SELECT CASE (case_vec)
      CASE(1)
        deallocate(g%RvecB)
      CASE(2)
        deallocate(g%CvecB)
      CASE(3)
        deallocate(g%RvecG)
      CASE(4)
        deallocate(g%CvecG)
      END SELECT
    ENDIF

#endif
  END SUBROUTINE MakeResidual_Davidson_MPI3
!=======================================================================================

!=======================================================================================
! MPI for calculating residual in the main Davidson procedure
! V4, current
!
! note: bounds_MPI is calculated here
!=======================================================================================
  SUBROUTINE MakeResidual_Davidson_MPI4(ndim,g,psi,Hpsi,Ene,Vec,conv,converge,         &
                                        VecToBeIncluded,tab_norm2g,norm2g,convergeResi,&
                                        convergeEne,fresidu,iresidu,nb_diago,epsi)
    USE mod_system
    USE mod_psi,        ONLY:param_psi
    USE mod_ana_psi_MPI,ONLY:norm_psi_MPI
    USE mod_psi_Op_MPI, ONLY:Set_symab_OF_psiBasisRep_MPI
    USE mod_MPI_aux
    IMPLICIT NONE

    TYPE(param_psi), intent(inout)              :: g
    TYPE(param_psi), intent(in)                 :: psi(:)
    TYPE(param_psi), intent(in)                 :: Hpsi(:)
    Real(kind=Rkind),intent(in)                 :: Ene(:)
    Real(kind=Rkind),intent(in)                 :: Vec(:,:)
    Real(kind=Rkind),intent(inout)              :: tab_norm2g(:)
    Real(kind=Rkind),intent(inout)              :: norm2g
    Real(kind=Rkind),intent(in)                 :: epsi
    Integer,         intent(in)                 :: ndim
    Integer,         intent(in)                 :: nb_diago
    Integer,         intent(inout)              :: fresidu
    Integer,         intent(inout)              :: iresidu
    Logical,         intent(in)                 :: VecToBeIncluded(:)
    Logical,         intent(in)                 :: convergeEne(:)
    Logical,         intent(inout)              :: convergeResi(:)
    Logical,         intent(inout)              :: converge(:)
    Logical,         intent(inout)              :: conv

    Integer                                     :: case_vec
    Integer                                     :: size_vec
    Integer                                     :: isym
    Integer                                     :: d1
    Integer                                     :: d2
    Integer                                     :: ii
    Integer                                     :: jj

#if(run_MPI)

    IF(keep_MPI) THEN
      IF(allocated(g%RvecB)) THEN   
        case_vec=1 
        size_vec=size(g%RvecB)
      ELSEIF(allocated(g%CvecB)) THEN
        case_vec=2
        size_vec=size(g%CvecB)
      ELSEIF(allocated(g%RvecG)) THEN
        case_vec=3
        size_vec=size(g%RvecG)
      ELSEIF(allocated(g%CvecG)) THEN
        case_vec=4
        size_vec=size(g%CvecG)
      ELSE                       
        case_vec=0
        STOP 'ERROR in g%vec of MakeResidual_Davidson_MPI'
      ENDIF
    ENDIF

    IF(MPI_scheme/=1) THEN
      CALL MPI_BCAST(case_vec,size1_MPI,Int_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
      CALL MPI_BCAST(size_vec,size1_MPI,Int_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
    ENDIF

    IF(.NOT. allocated(bounds_MPI)) allocate(bounds_MPI(1:2,0:MPI_np-1))
    nb_per_MPI=size_vec/MPI_np
    nb_rem_MPI=mod(size_vec,MPI_np) 
    DO i_MPI=0,MPI_np-1
      bounds_MPI(1,i_MPI)=i_MPI*nb_per_MPI+1+MIN(i_MPI,nb_rem_MPI)
      bounds_MPI(2,i_MPI)=(i_MPI+1)*nb_per_MPI+MIN(i_MPI,nb_rem_MPI)                   &
                                              +merge(1,0,nb_rem_MPI>i_MPI)
    ENDDO

    d1=bounds_MPI(1,MPI_id)
    d2=bounds_MPI(2,MPI_id)

    ! IF(MPI_id/=0 .AND. MPI_scheme/=1) THEN
    IF(.NOT. keep_MPI) THEN
      SELECT CASE (case_vec)
      CASE(1)
        CALL allocate_array(g%RvecB,d1,d2)
      CASE(2)
        CALL allocate_array(g%CvecB,d1,d2)
      CASE(3)
        CALL allocate_array(g%RvecG,d1,d2)
      CASE(4)
        CALL allocate_array(g%CvecG,d1,d2)
      END SELECT
    ENDIF

    !-----------------------------------------------------------------------------------
    ! converge,VecToBeIncluded,fresidu,norm2g,epsi,convergeEne synchronized 
    !CALL MPI_BCAST(converge,ndim,MPI_LOGICAL,root_MPI,MPI_COMM_WORLD,MPI_err)
    DO jj=1,ndim
      IF (.NOT. converge(jj) .AND. VecToBeIncluded(jj)) THEN
        isym=maxloc(abs(Vec(:,jj)),dim=1) ! to find the rigth symmetry
        
        CALL Residual_Davidson_sum_MPI(g,Hpsi,psi,Vec,Ene,ndim,case_vec,size_vec,jj)

        CALL Set_symab_OF_psiBasisRep_MPI(g,symab=psi(isym)%symab)
        CALL norm_psi_MPI(g,1)
        tab_norm2g(jj)=sqrt(g%norm2)
      ENDIF ! for .NOT. converge(jj) .AND. VecToBeIncluded(jj)  
    ENDDO ! for jj=1,ndim

    DO jj=1,ndim
      IF (.NOT. converge(jj) .AND. VecToBeIncluded(jj)) THEN
        IF(fresidu==0) fresidu=jj
        IF(tab_norm2g(jj)>norm2g) THEN
          iresidu=jj
          norm2g=tab_norm2g(iresidu)
        ENDIF
        convergeResi(jj)=tab_norm2g(jj)<epsi
      ENDIF ! for .NOT. converge(jj) .AND. VecToBeIncluded(jj)
      converge(jj)=(convergeEne(jj) .AND. convergeResi(jj))
    ENDDO ! for jj=1,ndim
    conv=all(converge(1:nb_diago))
    !-----------------------------------------------------------------------------------

    ! IF(MPI_id/=0 .AND. MPI_scheme/=1) THEN
    IF(.NOT. keep_MPI) THEN
      SELECT CASE (case_vec)
      CASE(1)
        deallocate(g%RvecB)
      CASE(2)
        deallocate(g%CvecB)
      CASE(3)
        deallocate(g%RvecG)
      CASE(4)
        deallocate(g%CvecG)
      END SELECT
    ENDIF

#endif
  END SUBROUTINE MakeResidual_Davidson_MPI4
!=======================================================================================

!=======================================================================================
!> @brief MPI for calculating residual at jth
!=======================================================================================
  SUBROUTINE MakeResidual_Davidson_j_MPI3(jj,g,psi,Hpsi,Ene,Vec)
    USE mod_system
    USE mod_psi,       ONLY:param_psi
    USE mod_psi_Op_MPI,ONLY:Set_symab_OF_psiBasisRep_MPI
    USE mod_MPI_aux
    IMPLICIT NONE
    
    TYPE(param_psi), intent(inout)              :: g
    TYPE(param_psi), intent(in)                 :: psi(:)
    TYPE(param_psi), intent(in)                 :: Hpsi(:)
    Real(kind=Rkind),intent(in)                 :: Ene(:)
    Real(kind=Rkind),intent(in)                 :: Vec(:,:)
    Integer         ,intent(in)                 :: jj

    Integer                                     :: case_vec
    Integer                                     :: size_vec
    Integer                                     :: isym
    Integer                                     :: ndim
    Integer                                     :: d1
    Integer                                     :: d2
    Integer                                     :: ii

#if(run_MPI)

    ! IF(MPI_id==0 .OR. MPI_scheme==1) THEN
    IF(keep_MPI) THEN
      IF(allocated(g%RvecB)) THEN   
        case_vec=1 
        size_vec=size(g%RvecB)
      ELSEIF(allocated(g%CvecB)) THEN
        case_vec=2
        size_vec=size(g%CvecB)
      ELSEIF(allocated(g%RvecG)) THEN
        case_vec=3
        size_vec=size(g%RvecG)
      ELSEIF(allocated(g%CvecG)) THEN
        case_vec=4
        size_vec=size(g%CvecG)
      ELSE                       
        case_vec=0
        STOP 'ERROR in g%vec of MakeResidual_Davidson_MPI'
      ENDIF
    ENDIF

    IF(MPI_scheme/=1) THEN
      CALL MPI_BCAST(case_vec,size1_MPI,Int_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
      CALL MPI_BCAST(size_vec,size1_MPI,Int_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
    ENDIF

    ndim=size(Vec,dim=1)
    
    d1=bounds_MPI(1,MPI_id)
    d2=bounds_MPI(2,MPI_id)

    !> allocate g%RvecB for slave threads 
    ! IF(MPI_id/=0 .AND. MPI_scheme/=1 ) THEN
    IF(.NOT. keep_MPI) THEN
      SELECT CASE (case_vec)
      CASE(1)
        CALL allocate_array(g%RvecB,d1,d2)
      CASE(2)
        CALL allocate_array(g%CvecB,d1,d2)
      CASE(3)
        CALL allocate_array(g%RvecG,d1,d2)
      CASE(4)
        CALL allocate_array(g%CvecG,d1,d2)
      END SELECT
    ENDIF

    isym = maxloc(abs(Vec(:,jj)),dim=1) ! to find the rigth symmetry

    CALL Residual_Davidson_sum_MPI(g,Hpsi,psi,Vec,Ene,ndim,case_vec,size_vec,jj)
    
    ! IF(MPI_id==0) CALL Set_symab_OF_psiBasisRep(g,symab=psi(isym)%symab)
    CALL Set_symab_OF_psiBasisRep_MPI(g,symab=psi(isym)%symab)

    !> deallocate g%RvecB
    ! IF(MPI_id/=0 .AND. MPI_scheme/=1) THEN
    IF(.NOT. keep_MPI) THEN
      SELECT CASE (case_vec)
      CASE(1)
        deallocate(g%RvecB)
      CASE(2)
        deallocate(g%CvecB)
      CASE(3)
        deallocate(g%RvecG)
      CASE(4)
        deallocate(g%CvecG)
      END SELECT
    ENDIF
    !-----------------------------------------------------------------------------------
    
    ! call MakeResidual_Davidson_core instead, but with repecct allcoation of Rvec
    !CALL MakeResidual_Davidson_core(jj,g,psi,Hpsi,Ene,Vec,case_vec,size_vec,ndim)

#endif
  END SUBROUTINE MakeResidual_Davidson_j_MPI3
!=======================================================================================


!=======================================================================================
!> @brief MPI for calculating residual at jth
! not used any more, delete next update
!=======================================================================================
  SUBROUTINE MakeResidual_Davidson_core(jj,g,psi,Hpsi,Ene,Vec,case_vec,size_vec,ndim)
    USE mod_system
    USE mod_psi,     ONLY : param_psi,Set_symab_OF_psiBasisRep
    IMPLICIT NONE
    
    TYPE(param_psi), intent(inout)              :: g
    TYPE(param_psi), intent(in)                 :: psi(:)
    TYPE(param_psi), intent(in)                 :: Hpsi(:)
    Real(kind=Rkind),intent(in)                 :: Ene(:)
    Real(kind=Rkind),intent(in)                 :: Vec(:,:)
    Integer         ,intent(in)                 :: case_vec
    Integer         ,intent(in)                 :: size_vec
    Integer         ,intent(in)                 :: ndim
    Integer         ,intent(in)                 :: jj

    Real(kind=Rkind),allocatable                :: Rvec(:)
    Complex(kind=Rkind),allocatable             :: Cvec(:)
    Integer                                     :: isym
    Integer                                     :: ii

#if(run_MPI)

    nb_per_MPI=size_vec/MPI_np
    nb_rem_MPI=mod(size_vec,MPI_np) 
    
    bound1_MPI=MPI_id*nb_per_MPI+1+MIN(MPI_id,nb_rem_MPI)
    bound2_MPI=(MPI_id+1)*nb_per_MPI+MIN(MPI_id,nb_rem_MPI)+merge(1,0,nb_rem_MPI>MPI_id)
    
    SELECT CASE (case_vec)
    CASE(1,3)
      allocate(Rvec(size_vec))
    CASE(2,4)
      allocate(Cvec(size_vec))
    END SELECT
    
    isym = maxloc(abs(Vec(:,jj)),dim=1) ! to find the rigth symmetry
    
    SELECT CASE (case_vec) 
    CASE(1)
      Rvec=ZERO
      DO ii=1,ndim
        Rvec(bound1_MPI:bound2_MPI)=Rvec(bound1_MPI:bound2_MPI)                        &
                         +Hpsi(ii)%RvecB(bound1_MPI:bound2_MPI)*Vec(ii,jj)             &
                          -psi(ii)%RvecB(bound1_MPI:bound2_MPI)*(Ene(jj)*Vec(ii,jj))
      ENDDO
      CALL MPI_Reduce(Rvec,g%RvecB,size_vec,MPI_Real8,MPI_SUM,root_MPI,                &
                      MPI_COMM_WORLD,MPI_err)
    CASE(2)
      Cvec=ZERO
      DO ii=1,ndim
        Cvec(bound1_MPI:bound2_MPI)=Cvec(bound1_MPI:bound2_MPI)                        &
                         +Hpsi(ii)%CvecB(bound1_MPI:bound2_MPI)*Vec(ii,jj)             &
                          -psi(ii)%CvecB(bound1_MPI:bound2_MPI)*(Ene(jj)*Vec(ii,jj))
      ENDDO
      CALL MPI_Reduce(Cvec,g%CvecB,size_vec,MPI_Complex8,MPI_SUM,root_MPI,             &
                      MPI_COMM_WORLD,MPI_err)
    CASE(3)
      Rvec=ZERO
      DO ii=1,ndim
        Rvec(bound1_MPI:bound2_MPI)=Rvec(bound1_MPI:bound2_MPI)                        &
                      +Hpsi(ii)%RvecG(bound1_MPI:bound2_MPI)*Vec(ii,jj)                &
                       -psi(ii)%RvecG(bound1_MPI:bound2_MPI)*(Ene(jj)*Vec(ii,jj))
      ENDDO
      CALL MPI_Reduce(Rvec,g%RvecG,size_vec,MPI_Real8,MPI_SUM,root_MPI,                &
                      MPI_COMM_WORLD,MPI_err)
    CASE(4)
      Cvec=ZERO
      DO ii=1,ndim
        Cvec(bound1_MPI:bound2_MPI)=Cvec(bound1_MPI:bound2_MPI)                        &
                         +Hpsi(ii)%CvecG(bound1_MPI:bound2_MPI)*Vec(ii,jj)             &
                          -psi(ii)%CvecG(bound1_MPI:bound2_MPI)*(Ene(jj)*Vec(ii,jj))
      ENDDO
      CALL MPI_Reduce(Cvec,g%CvecG,size_vec,MPI_Complex8,MPI_SUM,root_MPI,             &
                      MPI_COMM_WORLD,MPI_err)
    END SELECT
    
    IF(MPI_id==0) CALL Set_symab_OF_psiBasisRep(g,symab=psi(isym)%symab)
    
    SELECT CASE (case_vec)
    CASE(1,3)
      deallocate(Rvec)
    CASE(2,4)
      deallocate(Cvec)
    END SELECT

#endif
  ENDSUBROUTINE MakeResidual_Davidson_core
!=======================================================================================


!=======================================================================================
!> @brief: summary part for MakeResidual_Davidson_MPI
!=======================================================================================
  SUBROUTINE Residual_Davidson_sum_MPI(g,Hpsi,psi,Vec,Ene,ndim,case_vec,size_vec,jj)
    USE mod_system
    USE mod_psi,ONLY:param_psi
    USE mod_MPI_aux
    IMPLICIT NONE
    
    TYPE(param_psi),              intent(inout) :: g
    TYPE(param_psi),              intent(in)    :: psi(:)
    TYPE(param_psi),              intent(in)    :: Hpsi(:)
    Real(kind=Rkind),             intent(in)    :: Ene(:)
    Real(kind=Rkind),             intent(in)    :: Vec(:,:)
    Integer,                      intent(in)    :: ndim
    Integer,                      intent(in)    :: case_vec
    Integer,                      intent(in)    :: size_vec
    Integer,                      intent(in)    :: jj

    Integer                                     :: d1
    Integer                                     :: d2
    Integer                                     :: ii

#if(run_MPI)

    d1=bounds_MPI(1,MPI_id)
    d2=bounds_MPI(2,MPI_id)

    SELECT CASE (case_vec) 
    CASE(1)
      g%RvecB=ZERO
      DO ii=1,ndim
        g%RvecB(d1:d2)=g%RvecB(d1:d2)+Hpsi(ii)%RvecB(d1:d2)*Vec(ii,jj)                 &
                                      -psi(ii)%RvecB(d1:d2)*(Ene(jj)*Vec(ii,jj))
      ENDDO
      ! MPI_Gatherv or MPI_Reduce alternate, avoid creaction of extra memory
      CALL MPI_combine_array(g%RvecB,MS=MPI_scheme)

    CASE(2)
      g%CvecB=ZERO
      DO ii=1,ndim
        g%CvecB(d1:d2)=g%CvecB(d1:d2)+Hpsi(ii)%CvecB(d1:d2)*Vec(ii,jj)                 &
                                      -psi(ii)%CvecB(d1:d2)*(Ene(jj)*Vec(ii,jj))
      ENDDO
      CALL MPI_combine_array(g%CvecB,MS=MPI_scheme)

    CASE(3)
      g%RvecG=ZERO
      DO ii=1,ndim
        g%RvecG(d1:d2)=g%RvecG(d1:d2)+Hpsi(ii)%RvecG(d1:d2)*Vec(ii,jj)                 &
                                      -psi(ii)%RvecG(d1:d2)*(Ene(jj)*Vec(ii,jj))
      ENDDO
      CALL MPI_combine_array(g%RvecG,MS=MPI_scheme)
    CASE(4)
      g%CvecG=ZERO
      DO ii=1,ndim
        g%CvecG(d1:d2)=g%CvecG(d1:d2)+Hpsi(ii)%CvecG(d1:d2)*Vec(ii,jj)                 &
                                      -psi(ii)%CvecG(d1:d2)*(Ene(jj)*Vec(ii,jj))
      ENDDO
      CALL MPI_combine_array(g%CvecG,MS=MPI_scheme)

    END SELECT  

#endif
  ENDSUBROUTINE Residual_Davidson_sum_MPI
!=======================================================================================


!=======================================================================================
! MPI for calculating residual in the main Davidson procedure
!=======================================================================================
  SUBROUTINE MakeResidual_Davidson_MPI(ndim,g,psi,Hpsi,Ene,Vec,conv,converge,          &
                                       VecToBeIncluded,tab_norm2g,norm2g,convergeResi, &
                                       convergeEne,fresidu,iresidu,nb_diago,epsi)
    USE mod_system
    USE mod_psi,     ONLY : param_psi,norm2_psi,Set_symab_OF_psiBasisRep
    USE mod_propa,   ONLY : param_Davidson
    IMPLICIT NONE
    
    TYPE(param_psi), intent(inout)              :: g
    TYPE(param_psi), intent(in)                 :: psi(:)
    TYPE(param_psi), intent(in)                 :: Hpsi(:)
    Real(kind=Rkind),intent(in)                 :: Ene(:)
    Real(kind=Rkind),intent(in)                 :: Vec(:,:)
    Real(kind=Rkind),intent(inout)              :: tab_norm2g(:)
    Real(kind=Rkind),intent(inout)              :: norm2g
    Real(kind=Rkind),intent(in)                 :: epsi
    Integer,         intent(in)                 :: ndim
    Integer,         intent(in)                 :: nb_diago
    Integer,         intent(inout)              :: fresidu
    Integer,         intent(inout)              :: iresidu
    Logical,         intent(in)                 :: VecToBeIncluded(:)
    Logical,         intent(in)                 :: convergeEne(:)
    Logical,         intent(inout)              :: convergeResi(:)
    Logical,         intent(inout)              :: converge(:)
    Logical,         intent(inout)              :: conv

    Real(kind=Rkind),allocatable                :: Rvec(:)
    Complex(kind=Rkind),allocatable             :: Cvec(:)
    Integer                                     :: case_vec
    Integer                                     :: size_vec
    Integer                                     :: isym
    Integer                                     :: ii
    Integer                                     :: jj

#if(run_MPI)

    IF(MPI_id==0) THEN
      IF(allocated(g%RvecB)) THEN   
        g%RvecB=ZERO
        case_vec=1 
        size_vec=size(g%RvecB)
      ELSEIF(allocated(g%CvecB)) THEN
        g%CvecB=ZERO
        case_vec=2
        size_vec=size(g%CvecB)
      ELSEIF(allocated(g%RvecG)) THEN
        g%RvecG=ZERO
        case_vec=3
        size_vec=size(g%RvecG)
      ELSEIF(allocated(g%CvecG)) THEN
        g%CvecG=ZERO
        case_vec=4
        size_vec=size(g%CvecG)
      ELSE                       
        case_vec=0
        STOP 'ERROR in g%vec of MakeResidual_Davidson_MPI'
      ENDIF
    ENDIF
    CALL MPI_BCAST(case_vec,size1_MPI,Int_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
    CALL MPI_BCAST(size_vec,size1_MPI,Int_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
    
    SELECT CASE (case_vec)
    CASE(1,3)
      allocate(Rvec(size_vec))
    CASE(2,4)
      allocate(Cvec(size_vec))
    END SELECT
    
    nb_per_MPI=ndim/MPI_np
    nb_rem_MPI=mod(ndim,MPI_np) !remainder jobs
    
    bound1_MPI=MPI_id*nb_per_MPI+1+MIN(MPI_id,nb_rem_MPI)
    bound2_MPI=(MPI_id+1)*nb_per_MPI+MIN(MPI_id,nb_rem_MPI)+merge(1,0,nb_rem_MPI>MPI_id)
    
    ! converge,VecToBeIncluded,fresidu,norm2g,epsi,convergeEne synchronized 
    !CALL MPI_BCAST(converge,ndim,MPI_LOGICAL,root_MPI,MPI_COMM_WORLD,MPI_err)
    DO jj=1,ndim
      IF (.NOT. converge(jj) .AND. VecToBeIncluded(jj)) THEN
        isym = maxloc(abs(Vec(:,jj)),dim=1) ! to find the rigth symmetry

        SELECT CASE (case_vec) 
        CASE(1)
          Rvec=ZERO
          DO ii=bound1_MPI,bound2_MPI
            Rvec=Rvec+Hpsi(ii)%RvecB*Vec(ii,jj)-psi(ii)%RvecB*(Ene(jj)*Vec(ii,jj))
          ENDDO
          CALL MPI_Reduce(Rvec,g%RvecB,size_vec,MPI_Real8,MPI_SUM,root_MPI,            &
                          MPI_COMM_WORLD,MPI_err)
        CASE(2)
          Cvec=ZERO
          DO ii=bound1_MPI,bound2_MPI
            Cvec=Cvec+Hpsi(ii)%CvecB*Vec(ii,jj)-psi(ii)%CvecB*(Ene(jj)*Vec(ii,jj))
          ENDDO
          CALL MPI_Reduce(Cvec,g%CvecB,size_vec,MPI_Complex8,MPI_SUM,root_MPI,         &
                          MPI_COMM_WORLD,MPI_err)
        CASE(3)
          Rvec=ZERO
          DO ii=bound1_MPI,bound2_MPI
            Rvec=Rvec+Hpsi(ii)%RvecG*Vec(ii,jj)-psi(ii)%RvecG*(Ene(jj)*Vec(ii,jj))
          ENDDO
          CALL MPI_Reduce(Rvec,g%RvecG,size_vec,MPI_Real8,MPI_SUM,root_MPI,            &
                          MPI_COMM_WORLD,MPI_err)
        CASE(4)
          Cvec=ZERO
          DO ii=bound1_MPI,bound2_MPI
            Cvec=Cvec+Hpsi(ii)%CvecG*Vec(ii,jj)-psi(ii)%CvecG*(Ene(jj)*Vec(ii,jj))
          ENDDO
          CALL MPI_Reduce(Cvec,g%CvecG,size_vec,MPI_Complex8,MPI_SUM,root_MPI,         &
                          MPI_COMM_WORLD,MPI_err)
        END SELECT
        
        IF(MPI_id==0) THEN
          CALL Set_symab_OF_psiBasisRep(g,symab=psi(isym)%symab)
          CALL norm2_psi(g)
          tab_norm2g(jj) = sqrt(g%norm2)
        ENDIF
      ENDIF ! for .NOT. converge(jj) .AND. VecToBeIncluded(jj)  
        
    ENDDO ! for jj=1,ndim

    CALL MPI_BCAST(tab_norm2g,ndim,MPI_Real8,root_MPI,MPI_COMM_WORLD,MPI_err)
    
    DO jj=1,ndim
      IF (.NOT. converge(jj) .AND. VecToBeIncluded(jj)) THEN
        IF(fresidu==0) fresidu=jj
        IF(tab_norm2g(jj)>norm2g) THEN
          iresidu=jj
          norm2g=tab_norm2g(iresidu)
        ENDIF
        convergeResi(jj)=tab_norm2g(jj)<epsi
      ENDIF ! for .NOT. converge(jj) .AND. VecToBeIncluded(jj)
      converge(jj) = (convergeEne(jj) .AND. convergeResi(jj))
    ENDDO ! for jj=1,ndim
    conv = all(converge(1:nb_diago))

    SELECT CASE (case_vec)
    CASE(1,3)
      deallocate(Rvec)
    CASE(2,4)
      deallocate(Cvec)
    END SELECT

#endif
  END SUBROUTINE MakeResidual_Davidson_MPI
!=======================================================================================


!=======================================================================================
! MPI for calculating residual at jth
!=======================================================================================
  SUBROUTINE MakeResidual_Davidson_j_MPI(jj,g,psi,Hpsi,Ene,Vec)
    USE mod_system
    USE mod_psi,     ONLY : param_psi,Set_symab_OF_psiBasisRep
    USE mod_propa,   ONLY : param_Davidson
    IMPLICIT NONE
    
    TYPE(param_psi), intent(inout)               :: g
    TYPE(param_psi), intent(in)                  :: psi(:)
    TYPE(param_psi), intent(in)                  :: Hpsi(:)
    Real(kind=Rkind),intent(in)                  :: Ene(:)
    Real(kind=Rkind),intent(in)                  :: Vec(:,:)

    Real(kind=Rkind),allocatable                 :: Rvec(:)
    Complex(kind=Rkind),allocatable              :: Cvec(:)
    Integer                                      :: case_vec
    Integer                                      :: size_vec
    Integer                                      :: isym
    Integer                                      :: ndim
    Integer                                      :: ii
    Integer                                      :: jj

#if(run_MPI)

    IF(MPI_id==0) THEN
      IF(allocated(g%RvecB)) THEN   
        g%RvecB=ZERO
        case_vec=1 
        size_vec=size(g%RvecB)
      ELSEIF(allocated(g%CvecB)) THEN
        g%CvecB=ZERO
        case_vec=2
        size_vec=size(g%CvecB)
      ELSEIF(allocated(g%RvecG)) THEN
        g%RvecG=ZERO
        case_vec=3
        size_vec=size(g%RvecG)
      ELSEIF(allocated(g%CvecG)) THEN
        g%CvecG=ZERO
        case_vec=4
        size_vec=size(g%CvecG)
      ELSE                       
        case_vec=0
        STOP 'ERROR in g%vec of MakeResidual_Davidson_MPI'
      ENDIF
    ENDIF
    CALL MPI_BCAST(case_vec,size1_MPI,Int_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
    CALL MPI_BCAST(size_vec,size1_MPI,Int_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
    
    SELECT CASE (case_vec)
    CASE(1,3)
      allocate(Rvec(size_vec))
    CASE(2,4)
      allocate(Cvec(size_vec))
    END SELECT
    
    ndim=size(Vec,dim=1)
    nb_per_MPI=ndim/MPI_np
    nb_rem_MPI=mod(ndim,MPI_np) !remainder jobs
    
    bound1_MPI=MPI_id*nb_per_MPI+1+MIN(MPI_id,nb_rem_MPI)
    bound2_MPI=(MPI_id+1)*nb_per_MPI+MIN(MPI_id,nb_rem_MPI)+merge(1,0,nb_rem_MPI>MPI_id)
    
    ! converge,VecToBeIncluded,fresidu,norm2g,epsi,convergeEne synchronized 
    isym = maxloc(abs(Vec(:,jj)),dim=1) ! to find the rigth symmetry

    SELECT CASE (case_vec) 
    CASE(1)
      Rvec=ZERO
      DO ii=bound1_MPI,bound2_MPI
        Rvec=Rvec+Hpsi(ii)%RvecB*Vec(ii,jj)-psi(ii)%RvecB*(Ene(jj)*Vec(ii,jj))
      ENDDO
      CALL MPI_Reduce(Rvec,g%RvecB,size_vec,MPI_Real8,MPI_SUM,root_MPI,                  &
                      MPI_COMM_WORLD,MPI_err)
    CASE(2)
      Cvec=ZERO
      DO ii=bound1_MPI,bound2_MPI
        Cvec=Cvec+Hpsi(ii)%CvecB*Vec(ii,jj)-psi(ii)%CvecB*(Ene(jj)*Vec(ii,jj))
      ENDDO
      CALL MPI_Reduce(Cvec,g%CvecB,size_vec,MPI_Complex8,MPI_SUM,root_MPI,               &
                      MPI_COMM_WORLD,MPI_err)
    CASE(3)
      Rvec=ZERO
      DO ii=bound1_MPI,bound2_MPI
        Rvec=Rvec+Hpsi(ii)%RvecG*Vec(ii,jj)-psi(ii)%RvecG*(Ene(jj)*Vec(ii,jj))
      ENDDO
      CALL MPI_Reduce(Rvec,g%RvecG,size_vec,MPI_Real8,MPI_SUM,root_MPI,                  &
                      MPI_COMM_WORLD,MPI_err)
    CASE(4)
      Cvec=ZERO
      DO ii=bound1_MPI,bound2_MPI
        Cvec=Cvec+Hpsi(ii)%CvecG*Vec(ii,jj)-psi(ii)%CvecG*(Ene(jj)*Vec(ii,jj))
      ENDDO
      CALL MPI_Reduce(Cvec,g%CvecG,size_vec,MPI_Complex8,MPI_SUM,root_MPI,               &
                      MPI_COMM_WORLD,MPI_err)
    END SELECT
    
    IF(MPI_id==0) CALL Set_symab_OF_psiBasisRep(g,symab=psi(isym)%symab)

    SELECT CASE (case_vec)
    CASE(1,3)
      deallocate(Rvec)
    CASE(2,4)
      deallocate(Cvec)
    END SELECT

#endif
  END SUBROUTINE MakeResidual_Davidson_j_MPI
!=======================================================================================


!=======================================================================================
! make save_WP=.true. with external control
!=======================================================================================
  SUBROUTINE exit_Davidson_external_MPI(exit_Davidson,save_WP,it)
    Logical,                       intent(inout) :: exit_Davidson
    Logical,                       intent(inout) :: save_WP
    Integer,                       intent(in)    :: it
    Logical                                      :: exist
    Integer                                      :: stat

!!#if(run_MPI)

    IF(it>2) THEN
      INQUIRE(FILE='Davidson_exit',EXIST=exist) 
      IF(exist) THEN
        save_WP=.TRUE.
        exit_Davidson=.TRUE.
        CALL RENAME('Davidson_exit','done_exit')
      ENDIF
    ENDIF

!!#endif
  END SUBROUTINE exit_Davidson_external_MPI
!=======================================================================================

ENDMODULE mod_Davidson_MPI
