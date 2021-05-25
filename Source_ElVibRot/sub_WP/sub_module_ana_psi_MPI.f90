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
MODULE mod_ana_psi_MPI
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: norm_psi_MPI,share_psi_nodes_MPI

CONTAINS
!=======================================================================================      
!> MPI norm^2 of psi (BasisRep or GridRep)
!> NOTE: ReNorm controls if nromalizing the the vector or
!> just keep the nromalization constrant
!> ReNorm: 1 calculate normlization constant only; 
!!         2 normalize with existing normlization constant
!!         3 calculate and normalize
!
! Warning: the relevant part of RvecB should be pre-ready on the each threads
!          so this is not a normlization subroutine for general case
!=======================================================================================      
  SUBROUTINE norm_psi_MPI(psi,ReNorm,GridRep,BasisRep)
    USE mod_system
    USE mod_psi_set_alloc
    IMPLICIT NONE

    TYPE(param_psi),intent(inout)           :: psi
    Integer,intent(in)                      :: ReNorm
    Logical,optional,intent(in)             :: GridRep
    Logical,optional,intent(in)             :: BasisRep

    Real(kind=Rkind),allocatable            :: tab_WeightChannels(:,:)
    Real(kind=Rkind)                        :: WrhonD
    Real(kind=Rkind)                        :: temp
    Integer                                 :: i_qa
    Integer                                 :: i_qaie
    Integer                                 :: i_be
    Integer                                 :: i_bi
    Integer                                 :: i_ba
    Integer                                 :: i_baie
    Integer                                 :: ii_baie
    Integer                                 :: if_baie
    Integer                                 :: i_max_w
    Logical                                 :: norm2GridRep
    Logical                                 :: norm2BasisRep

#if(run_MPI)

    IF(present(GridRep)) THEN
      IF(present(BasisRep)) THEN
        norm2GridRep =GridRep
        norm2BasisRep=BasisRep
      ELSE
        norm2GridRep =GridRep
        norm2BasisRep=.FALSE.
      ENDIF
    ELSE
      IF(present(BasisRep)) THEN
        norm2BasisRep = BasisRep
        norm2GridRep  = .FALSE.
      ELSE
        IF(psi%BasisRep .AND. psi%GridRep) THEN
          norm2BasisRep = .TRUE.
          norm2GridRep  = .FALSE.
        ELSE
          norm2BasisRep = psi%BasisRep
          norm2GridRep  = psi%GridRep
        ENDIF
      ENDIF
    ENDIF

    IF (norm2GridRep .AND. norm2BasisRep) THEN
      write(out_unitp,*) ' ERROR in norm2_psi'
      write(out_unitp,*) ' norm2GridRep=t and norm2BasisRep=t !'
      write(out_unitp,*) ' BasisRep,GridRep',psi%BasisRep,psi%GridRep
      STOP
    ENDIF
    
    IF(ReNorm==1 .OR. ReNorm==3) THEN
      ! NOTE: tab_WeightChannels obtained correctly only on master
      CALL Channel_weight_MPI(tab_WeightChannels,psi,norm2GridRep,norm2BasisRep)

      IF(MPI_id==0) psi%norm2=sum(tab_WeightChannels)
      CALL MPI_Bcast(psi%norm2,size1_MPI,Real_MPI,root_MPI,MPI_COMM_WORLD,MPI_err)
      
      IF (allocated(tab_WeightChannels)) THEN
        CALL dealloc_NParray(tab_WeightChannels,"tab_WeightChannels","Channel_weight")
      END IF
    ENDIF
    
    IF(ReNorm==2 .OR. ReNorm==3) THEN
      IF (psi%norm2 .EQ. ZERO ) THEN
        write(out_unitp,*) ' ERROR in norm2_psi'
        write(out_unitp,*) ' the norm2 is zero !',psi%norm2
        STOP
      END IF
      
      temp=sqrt(ONE/psi%norm2)

      !-normalization---------------------------------------------------------------------
      IF(norm2GridRep) THEN
        IF(psi%cplx) THEN
          psi%CvecG=psi%CvecG*cmplx(temp,ZERO,kind=Rkind)
        ELSE
          psi%RvecG=psi%RvecG*temp
        ENDIF
      ELSE
        IF(psi%cplx) THEN
          psi%CvecB=psi%CvecB*cmplx(temp,ZERO,kind=Rkind)
        ELSE
          psi%RvecB=psi%RvecB*temp
        ENDIF
      ENDIF
      psi%norm2=ONE
    ENDIF
      
    IF(.NOT. (ReNorm==1 .OR. ReNorm==2 .OR. ReNorm==3)) THEN
      STOP 'error in norm_psi_MPI, ReNorm not provided '
    ENDIF ! for ReNorm

#endif
  END SUBROUTINE norm_psi_MPI
!=======================================================================================


!=======================================================================================  
!> calculate weight with MPI
!> be careful with the way the vector distributed 
!=======================================================================================  
  SUBROUTINE Channel_weight_MPI(tab_WeightChannels,psi,                                &
                                GridRep,BasisRep,Dominant_Channel)
    USE mod_system
    USE mod_psi_set_alloc
    USE mod_basis,        ONLY:Rec_WrhonD
    USE mod_param_SGType2,ONLY:OldParam
    USE mod_ana_psi,      ONLY:Channel_weight_contracHADA
    USE mod_MPI_aux
    IMPLICIT NONE

    TYPE(param_psi), intent(in)                   :: psi
    Real(kind=Rkind),intent(inout),allocatable    :: tab_WeightChannels(:,:)
    Integer,optional,intent(inout)                :: Dominant_Channel(2)
    Logical,         intent(in)                   :: GridRep,BasisRep

    TYPE(OldParam)                                :: OldPara
    Real(kind=Rkind)                              :: WeightSG
    Real(kind=Rkind)                              :: WrhonD
    Real(kind=Rkind)                              :: max_w
    Real(kind=Rkind)                              :: temp
    Integer                                       :: i_qa
    Integer                                       :: i_qaie
    Integer                                       :: i_be
    Integer                                       :: i_bi
    Integer                                       :: i_ba
    Integer                                       :: i_baie
    Integer                                       :: ii_baie
    Integer                                       :: if_baie
    Integer                                       :: nb_be
    Integer                                       :: nb_bi
    Integer                                       :: iSG
    Integer                                       :: iqSG

#if(run_MPI)

    nb_be=get_nb_be_FROM_psi(psi)
    nb_bi=get_nb_bi_FROM_psi(psi)

    IF(GridRep .AND. BasisRep) THEN
      write(out_unitp,*) ' ERROR in Channel_weight'
      write(out_unitp,*) ' GridRep=t and BasisRep=t !'
      STOP
    END IF

    IF(.NOT. allocated(tab_WeightChannels) .AND. nb_bi > 0 .AND. nb_be > 0) THEN
      CALL alloc_NParray(tab_WeightChannels,(/nb_bi,nb_be/),                             &
                        "tab_WeightChannels","Channel_weight")
      tab_WeightChannels(:,:) = ZERO
    END IF

    IF(psi%para_AllBasis%basis_ext2n%contrac_ba_ON_HAC) THEN
      IF(MPI_id==0) CALL Channel_weight_contracHADA(tab_WeightChannels(:,1),psi)
    ELSE IF (psi%nb_baie==psi%nb_tot) THEN
      IF(BasisRep .AND. (allocated(psi%CvecB) .OR. allocated(psi%RvecB)) ) THEN
        nb_per_MPI=psi%nb_tot/MPI_np
        nb_rem_MPI=mod(psi%nb_tot,MPI_np) 
        bound1_MPI=MPI_id*nb_per_MPI+1+MIN(MPI_id,nb_rem_MPI)
        bound2_MPI=(MPI_id+1)*nb_per_MPI+MIN(MPI_id,nb_rem_MPI)                          &
                                        +merge(1,0,nb_rem_MPI>MPI_id)

        DO i_be=1,psi%nb_be
        DO i_bi=1,psi%nb_bi
          ii_baie=1+((i_bi-1)+(i_be-1)*psi%nb_bi)*psi%nb_ba
          if_baie=ii_baie-1+psi%nb_ba

          IF(bound2_MPI>ii_baie .AND. bound1_MPI<if_baie) THEN
            if_baie=MIN(if_baie,bound2_MPI)
            ii_baie=MAX(ii_baie,bound1_MPI)
            IF(psi%cplx) THEN
              tab_WeightChannels(i_bi,i_be)=real(dot_product(psi%CvecB(ii_baie:if_baie), &
                                                 psi%CvecB(ii_baie:if_baie)),kind=Rkind)
            ELSE
              tab_WeightChannels(i_bi,i_be)=dot_product(psi%RvecB(ii_baie:if_baie),      &
                                                        psi%RvecB(ii_baie:if_baie))
            ENDIF
          ENDIF
        END DO
        END DO

        CALL MPI_Reduce_sum_matrix(tab_WeightChannels,ONE_1,nb_bi,ONE_1,nb_be,root_MPI)

      ELSE IF(GridRep .AND. (allocated(psi%CvecG) .OR. allocated(psi%RvecG))) THEN
        tab_WeightChannels(:,:) = ZERO
        
        nb_per_MPI=psi%nb_baie/MPI_np
        nb_rem_MPI=mod(psi%nb_baie,MPI_np) 
        bound1_MPI=MPI_id*nb_per_MPI+1+MIN(MPI_id,nb_rem_MPI)
        bound2_MPI=(MPI_id+1)*nb_per_MPI+MIN(MPI_id,nb_rem_MPI)                          &
                                        +merge(1,0,nb_rem_MPI>MPI_id)
        
        DO i_qa=1,psi%nb_qa
          WrhonD=Rec_WrhonD(psi%BasisnD,i_qa)
          DO i_be=1,nb_be
          DO i_bi=1,nb_bi
            i_qaie=i_qa+((i_bi-1)+(i_be-1)*nb_bi )*psi%nb_qa

            IF(bound2_MPI>i_qaie .AND. bound1_MPI<i_qaie) THEN            
              IF(psi%cplx) THEN
                temp=abs(psi%CvecG(i_qaie))
              ELSE
                temp=psi%RvecG(i_qaie)
              END IF
              tab_WeightChannels(i_bi,i_be)=tab_WeightChannels(i_bi,i_be)+WrhonD*temp**2
            ENDIF
          ENDDO
          ENDDO
        ENDDO
       
        CALL MPI_Reduce_sum_matrix(tab_WeightChannels,ONE_1,nb_bi,ONE_1,nb_be,root_MPI)
        
      ELSE
        write(out_unitp,*) ' ERROR in Channel_weight',' from ',MPI_id
        IF (GridRep)  write(out_unitp,*) ' impossible to calculate the weights with the Grid'
        IF (BasisRep) write(out_unitp,*) ' impossible to calculate the weights with the Basis'
        STOP
      END IF
    ELSE
      !-To deal with a spectral representation--------------------------------------------
      tab_WeightChannels(:,:)=ZERO
      
      nb_per_MPI=psi%nb_tot/MPI_np
      nb_rem_MPI=mod(psi%nb_tot,MPI_np)
      bound1_MPI=MPI_id*nb_per_MPI+1+MIN(MPI_id,nb_rem_MPI)
      bound2_MPI=(MPI_id+1)*nb_per_MPI+MIN(MPI_id,nb_rem_MPI)+merge(1,0,nb_rem_MPI>MPI_id)

      IF(psi%cplx) THEN
        temp_cplx=dot_product(psi%CvecB(bound1_MPI:bound2_MPI),                          &
                              psi%CvecB(bound1_MPI:bound2_MPI))
        CALL MPI_Reduce(temp_cplx,temp_cplx1,size1_MPI,MPI_Complex8,MPI_SUM,root_MPI,    &
                        MPI_COMM_WORLD,MPI_err)
        tab_WeightChannels(1,1)=real(temp_real1)
        !tab_WeightChannels(1,1)=real(dot_product(psi%CvecB,psi%CvecB),kind=Rkind)
      ELSE
        temp_real=dot_product(psi%RvecB(bound1_MPI:bound2_MPI),                          &
                              psi%RvecB(bound1_MPI:bound2_MPI))
        CALL MPI_Reduce(temp_real,temp_real1,size1_MPI,MPI_Real8,MPI_SUM,root_MPI,       &
                        MPI_COMM_WORLD,MPI_err)
        tab_WeightChannels(1,1)=temp_real1
        !tab_WeightChannels(1,1)=dot_product(psi%RvecB,psi%RvecB)
      ENDIF
    ENDIF

    IF(present(Dominant_Channel)) THEN
      Dominant_Channel(:)=1
      max_w              =ZERO
      DO i_be=1,nb_be
      DO i_bi=1,nb_bi
        IF(tab_WeightChannels(i_bi,i_be)>max_w) THEN
          max_w=tab_WeightChannels(i_bi,i_be)
          Dominant_Channel(:)=(/ i_be,i_bi /)
        END IF
      END DO
      END DO
    END IF

#endif

  END SUBROUTINE Channel_weight_MPI
!=======================================================================================  

!=======================================================================================
! share psi%RvecB with the other processors in MPI scheme 3
!=======================================================================================
  SUBROUTINE share_psi_nodes_MPI(psi,psi_size)
    USE mod_system
    USE mod_psi_set_alloc
    USE mod_MPI_aux

    TYPE(param_psi),               intent(inout) :: psi(:)
    Integer,                       intent(in)    :: psi_size

    Real(kind=Rkind),allocatable                 :: RvecB_all(:)
    Complex(kind=Rkind),allocatable              :: CvecB_all(:)
    Integer                                      :: size_vecB
    Integer                                      :: size_psi
    Integer                                      :: d1,d2
    Integer                                      :: ii

#if(run_MPI)

    If(MPI_id==0) THEN
      IF(psi(1)%cplx) THEN
        size_vecB=size(psi(1)%CvecB)
      ELSE
        size_vecB=size(psi(1)%RvecB)
      ENDIF
      size_psi=psi_size
    ENDIF
    CALL MPI_Bcast_(size_psi,  size1_MPI,root_MPI)
    CALL MPI_Bcast_(size_vecB,size1_MPI,root_MPI)

    ! IF(MPI_id==0) THEN
    !   CALL allocate_array(RvecB_all,1,size_vecB*size_psi)
    !   DO ii=1,size_psi
    !     d1=size_vecB*(ii-1)+1
    !     d2=size_vecB* ii
    !     RvecB_all(d1:d2)=psi(ii)%RvecB
    !   ENDDO

    !   DO ii=1,MPI_nodes_num-1
    !     CALL MPI_Send(RvecB_all,size_vecB*size_psi,Real_MPI,MPI_nodes_p00(ii),        &
    !                   MPI_nodes_p00(ii),MPI_COMM_WORLD,MPI_err)
    !   ENDDO

    ! ELSEIF(MPI_nodes_p0) THEN !---------------------------------------------------------
    !   CALL allocate_array(RvecB_all,1,size_vecB*size_psi)

    !   CALL MPI_Recv(RvecB_all,size_vecB*size_psi,Real_MPI,root_MPI,MPI_id,            &
    !                 MPI_COMM_WORLD,MPI_stat,MPI_err)

    !   DO ii=1,size_psi
    !     allocate(psi(ii)%RvecB(size_vecB))
    !     d1=size_vecB*(ii-1)+1
    !     d2=size_vecB* ii
    !     psi(ii)%RvecB=RvecB_all(d1:d2)
    !   ENDDO

    ! ENDIF
    ! IF(allocated(RvecB_all)) deallocate(RvecB_all)

    IF(psi(1)%cplx) THEN
      CALL allocate_array(CvecB_all,1,size_vecB*size_psi)
    ELSE
      CALL allocate_array(RvecB_all,1,size_vecB*size_psi)
    ENDIF
    IF(MPI_id==0) THEN
      DO ii=1,size_psi
        d1=size_vecB*(ii-1)+1
        d2=size_vecB* ii
        IF(psi(1)%cplx) THEN
          CvecB_all(d1:d2)=psi(ii)%CvecB
        ELSE
          RvecB_all(d1:d2)=psi(ii)%RvecB
        ENDIF
      ENDDO
    ENDIF

    IF(MPI_nodes_p0) THEN
      IF(psi(1)%cplx) THEN
        CALL MPI_BCAST(CvecB_all,size_vecB*size_psi,Cplx_MPI,root_MPI,                &
                       MPI_NODE_0_COMM,MPI_err)
      ELSE
        CALL MPI_BCAST(RvecB_all,size_vecB*size_psi,Real_MPI,root_MPI,                &
                       MPI_NODE_0_COMM,MPI_err)
      ENDIF
    ENDIF

    IF(MPI_id/=0) THEN
      DO ii=1,size_psi
        IF(psi(1)%cplx) THEN
          allocate(psi(ii)%CvecB(size_vecB))
        ELSE
          allocate(psi(ii)%RvecB(size_vecB))
        ENDIF

        d1=size_vecB*(ii-1)+1
        d2=size_vecB* ii

        IF(psi(1)%cplx) THEN
          psi(ii)%RvecB=CvecB_all(d1:d2)
        ELSE
          psi(ii)%RvecB=RvecB_all(d1:d2)
        ENDIF
      ENDDO
    ENDIF

    IF(allocated(RvecB_all)) deallocate(RvecB_all)
    IF(allocated(CvecB_all)) deallocate(CvecB_all)

    DO ii=1,size_psi
      IF(MPI_nodes_p0) THEN
        CALL MPI_BCAST(psi(ii)%norm2,size1_MPI,Real_MPI,root_MPI,MPI_NODE_0_COMM,MPI_err)
        CALL MPI_BCAST(psi(ii)%symab,size1_MPI,Int_MPI,root_MPI,MPI_NODE_0_COMM,MPI_err)
      ENDIF
    ENDDO

#endif
  END SUBROUTINE share_psi_nodes_MPI
!=======================================================================================

END MODULE mod_ana_psi_MPI

