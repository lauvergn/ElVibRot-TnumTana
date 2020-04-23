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
      MODULE mod_basis_BtoG_GtoB
      USE mod_system
      USE mod_basis_set_alloc
      USE mod_basis_BtoG_GtoB_SGType2
      IMPLICIT NONE

      CONTAINS
!      ==========================================================
!       VecGridRep <=> VecBasisRep
!      ==========================================================
      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
       RECURSIVE SUBROUTINE RecRVecG_TO_RvecB(RVecG,RvecB,nq,nb,basis_set)
       USE mod_basis_BtoG_GtoB_SGType4
       USE mod_basis_RCVec_SGType4
        IMPLICIT NONE


        TYPE (basis), intent(in)         :: basis_set
        integer, intent(in)              :: nq,nb
        real (kind=Rkind), intent(inout) :: RVecG(:)
        real (kind=Rkind), intent(inout) :: RvecB(:)

        real (kind=Rkind), allocatable   :: RTempB(:,:)
        real (kind=Rkind), allocatable   :: RTempG(:,:,:)

        integer                          :: nbb,ibb1,ibb2

        real (kind=Rkind), allocatable   :: RB(:)
        real (kind=Rkind), allocatable   :: RG(:)

        !real (kind=Rkind), allocatable   :: d0b(:,:)
        !real (kind=Rkind), allocatable   :: w(:)

        integer                          :: nnb,nb2,ib,ib2,newnb2

        integer                          :: nnq,nq2,iq,iq2,nqi,nbi
        integer                          :: ibasis
        integer                          :: i_SG,iq0_SG,iq1_SG,nq_SG
        integer                          :: nb_thread
        integer                          :: der00(2) = (/0,0/)

        TYPE(Type_SmolyakRep)            :: SRep ! smolyak rep for SparseGrid_type=4
        integer                          :: itabR,iG,nR

!----- for debuging --------------------------------------------------
        integer :: err_mem,memory
        character (len=*), parameter :: name_sub='RecRVecG_TO_RvecB'
        logical, parameter :: debug=.FALSE.
        !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
          write(out_unitp,*) 'RVecG(:)',RVecG(:)
        END IF
        IF (basis_set%cplx) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' the basis is complex but the'
          write(out_unitp,*) ' the vector is real !!'
          STOP
        END IF


        IF (basis_set%packed_done) THEN
          IF (NewBasisEl .AND. basis_set%ndim == 0) THEN
             RvecB(:) = RVecG
          ELSE

            nb_mult_GTOB = nb_mult_GTOB + int(nb,kind=ILkind)*size(RVecG,kind=ILkind)
          !write(out_unitp,*) 'nb_mult_GTOB,nb,size(RVecG)',nb_mult_BTOG,nb,size(RVecG)

!          write(out_unitp,*) 'nq,nb',nq,nb
!
!          write(out_unitp,*) 'asso d0b',associated(basis_set%dnRBGwrho%d0)
!          write(out_unitp,*) 'shape d0b',shape(basis_set%dnRBGwrho%d0)
!          write(out_unitp,*) 'shape RVecG',shape(RVecG)
!          write(out_unitp,*) 'shape RVecB',shape(RVecB)
!          flush(out_unitp)


             RvecB(:) = matmul( basis_set%dnRBGwrho%d0,RVecG)
           END IF

        ELSE ! basis_set%nb_basis MUST BE > 0
          IF (basis_set%nb_basis == 0 ) STOP ' ERROR with packed!!!'

          SELECT CASE (basis_set%SparseGrid_type)
          CASE (0) ! Direct product

            nnq  = nq
            nnb  = 1
            nb2  = 1
            nbb  = 1
            nq2  = 1

            CALL alloc_NParray(RTempB,(/ nnq,nbb /),"RTempB",name_sub)
            RTempB(:,:) = reshape(RVecG,shape=(/ nnq,nbb /))


            DO ibasis=1,basis_set%nb_basis

              nb2 = basis_set%tab_Pbasis(ibasis)%Pbasis%nb
              nq2 = get_nq_FROM_basis(basis_set%tab_Pbasis(ibasis)%Pbasis)
              nnq = nnq / nq2

              CALL alloc_NParray(RTempG,(/ nnq,nq2,nnb /),"RTempG",name_sub)
              RTempG(:,:,:) = reshape(RTempB,shape=(/ nnq,nq2,nnb /))
              !write(out_unitp,*) 'ibasis shape G',ibasis,shape(RTempG)
              nbb = sum(basis_set%Tab_OF_Tabnb2(ibasis)%vec)

!              write(out_unitp,*) 'G=>B ibasis,nnb,Tab_OF_Tabnb2',ibasis,nnb,   &
!                                        basis_set%Tab_OF_Tabnb2(ibasis)%vec

              CALL dealloc_NParray(RTempB,"RTempB",name_sub)
              CALL alloc_NParray(RTempB,(/ nnq,nbb /),"RTempB",name_sub)
              !write(out_unitp,*) 'ibasis shape B',ibasis,shape(RTempB)


              CALL alloc_NParray(RG,(/ nq2 /),"RG",name_sub)
              CALL alloc_NParray(RB,(/ nb2 /),"RB",name_sub)

              IF (basis_set%tab_Pbasis(ibasis)%Pbasis%packed) THEN

                ibb1 = 1
                ibb2 = 0
                DO ib=1,nnb
                  newnb2 = basis_set%Tab_OF_Tabnb2(ibasis)%vec(ib)
                  ibb2 = ibb2 + newnb2

                  DO iq=1,nnq

                    RTempB(iq,ibb1:ibb2) = matmul(      &
                      basis_set%tab_Pbasis(ibasis)%Pbasis%dnRBGwrho%d0(1:newnb2,:), &
                                                  RTempG(iq,:,ib))

                    nb_mult_GTOB = nb_mult_GTOB + int(newnb2,kind=ILkind)*size(RTempG(iq,:,ib),kind=ILkind)
                    !write(out_unitp,*) 'nb_mult_GTOB,nb,size(RVecG)',nb_mult_GTOB,newnb2,size(RG)
                  END DO

                  ibb1 = ibb1 + newnb2
                END DO

              ELSE

                ibb1 = 1
                ibb2 = 0
                DO ib=1,nnb
                  newnb2 = basis_set%Tab_OF_Tabnb2(ibasis)%vec(ib)
                  ibb2 = ibb2 + newnb2

                  DO iq=1,nnq

                  RG(:) = RTempG(iq,:,ib)
                  CALL RecRVecG_TO_RvecB(RG,RB,nq2,newnb2,            &
                                         basis_set%tab_Pbasis(ibasis)%Pbasis)

                  RTempB(iq,ibb1:ibb2) = RB(1:newnb2)

                  END DO

                  ibb1 = ibb1 + newnb2
                END DO
              END IF

              CALL dealloc_NParray(RB,"RB",name_sub)
              CALL dealloc_NParray(RG,"RG",name_sub)

              CALL dealloc_NParray(RTempG,"RTempG",name_sub)

              nnb = nbb

            END DO

            RvecB(:) = reshape(RTempB, shape=(/ nnq*nnb /) )
            CALL dealloc_NParray(RTempB,"RTempB",name_sub)

          CASE (1) ! Sparse basis (Smolyak 1st implementation)
            RvecB(:) = ZERO
            CALL alloc_NParray(RB,(/ nb /),"RB",name_sub)
            iq0_SG = 1
            iq1_SG = 0
            DO i_SG=1,basis_set%nb_SG
              nq_SG = get_nq_FROM_basis(basis_set%tab_PbasisSG(i_SG)%Pbasis)
              iq1_SG = iq1_SG + nq_SG
              CALL RecRVecG_TO_RvecB(RVecG(iq0_SG:iq1_SG),RB,nq_SG,nb,  &
                                     basis_set%tab_PbasisSG(i_SG)%Pbasis)
              iq0_SG = iq0_SG + nq_SG
              RvecB(:) = RvecB(:) + RB(:) * basis_set%WeightSG(i_SG)

            END DO
            CALL dealloc_NParray(RB,"RB",name_sub)

          CASE (2) ! Sparse basis (Smolyak 2d implementation)

            CALL sub_G_TO_B(RVecG,RVecB,basis_set%WeightSG,             &
                      basis_set%para_SGType2%nDind_SmolyakRep%Tab_DInd, &
                                     basis_set%nDindB%Tab_DInd,         &
                                     basis_set%tab_basisPrimSG,         &
                                     D=basis_set%nb_basis,              &
                                     LG=basis_set%L_SparseGrid,         &
                                     LB=basis_set%L_SparseBasis)

            !write(out_unitp,*) ' ERROR in ',name_sub
            !STOP 'SparseGrid_type=2'

          CASE (4) ! Sparse basis (Smolyak 4th implementation)

            !!  RVecG TO SRep
            CALL alloc2_SmolyakRep(SRep,basis_set%para_SGType2%nDind_SmolyakRep,&
                                   basis_set%tab_basisPrimSG,grid=.TRUE.,       &
                                   nb0=basis_set%para_SGType2%nb0)
            CALL tabR2bis_TO_SmolyakRep1(SRep,RVecG) ! on the grid

            IF (allocated(basis_set%para_SGType2%nDind_SmolyakRep%Tab_nDval)) THEN
              CALL GSmolyakRep_TO_BSmolyakRep(SRep,                     &
                    basis_set%para_SGType2%nDind_SmolyakRep%Tab_nDval,  &
                    basis_set%tab_basisPrimSG,basis_set%para_SGType2%nb0)
            ELSE
              CALL GSmolyakRep_TO3_BSmolyakRep(SRep,basis_set%para_SGType2, &
                                               basis_set%tab_basisPrimSG)
            END IF

            CALL SmolyakRepBasis_TO_tabPackedBasis(SRep,RVecB,          &
                             basis_set%tab_basisPrimSG,basis_set%nDindB,&
                             basis_set%para_SGType2,basis_set%WeightSG)

            CALL dealloc_SmolyakRep(SRep)

          CASE DEFAULT
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' WRONG SparseGrid_type',basis_set%SparseGrid_type
            write(out_unitp,*) ' The possibilities are: 0, 1, 2, 4'
            STOP
          END SELECT

        END IF

       IF (debug) THEN
         write(out_unitp,*) 'RvecB(:)',RvecB(:)
         write(out_unitp,*) 'END ',name_sub
       END IF
      END SUBROUTINE RecRVecG_TO_RvecB


      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      RECURSIVE SUBROUTINE RecCVecG_TO_CVecB(CVecG,CVecB,nq,nb,basis_set)
       USE mod_basis_BtoG_GtoB_SGType4
        IMPLICIT NONE

        TYPE (basis), intent(in)            :: basis_set
        integer, intent(in)                 :: nq,nb
        complex (kind=Rkind), intent(inout) :: CVecG(:)
        complex (kind=Rkind), intent(inout) :: CVecB(:)

        complex (kind=Rkind), allocatable   :: CTempG(:,:,:)
        complex (kind=Rkind), allocatable   :: CTempB(:,:)

        real (kind=Rkind), allocatable :: RVecG(:)
        real (kind=Rkind), allocatable :: RVecB(:)

        TYPE(Type_SmolyakRep)            :: SRep ! smolyak rep for SparseGrid_type=4
        integer                          :: itabR,iG,nR

        integer                             :: nbb,ibb1,ibb2

        complex (kind=Rkind), allocatable       :: CB(:)
        complex (kind=Rkind), allocatable       :: CG(:)
        complex (kind=Rkind), allocatable       :: d0cb(:,:)
        real    (kind=Rkind), allocatable       :: d0b(:,:)
        real    (kind=Rkind), allocatable       :: w(:)

        integer                             :: nnb,nb2,ib,ib2,newnb2
        integer                             :: nnq,nq2,iq,iq2

        integer                             :: ibasis
        integer                             :: i_SG,iq0_SG,iq1_SG,nq_SG
        integer                             :: nb_thread
        integer                             :: der00(2) = (/0,0/)

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='RecCVecG_TO_CVecB'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) ' BEGINNING in ',name_sub
        write(out_unitp,*) 'nq,nb',nq,nb
        write(out_unitp,*) 'shape CVecG',shape(CVecG)
        write(out_unitp,*) 'shape CVecB',shape(CVecB)
        CALL flush_perso(out_unitp)
      END IF

        IF (basis_set%cplx .AND. basis_set%packed_done) THEN
          IF (NewBasisEl .AND. basis_set%ndim == 0) THEN
            CVecB(:) = CVecG(:)
          ELSE
            CVecG(:) = CVecG(:) * cmplx(get_wrho_OF_basis(basis_set),kind=Rkind)
            DO ib=1,nb
              CVecB(ib) = sum( CVecG(:) * basis_set%dnCGB%d0(:,ib) )
            END DO
            !CVecB(1:nb) = matmul(CVecG(:),basis_set%dnCGB%d0(:,1:nb))
          END IF

        ELSE IF (.NOT. basis_set%cplx .AND. basis_set%packed_done) THEN
          IF (NewBasisEl .AND. basis_set%ndim == 0) THEN
            CVecB(:) = CVecG(:)
          ELSE
            CVecG(:) = CVecG(:) * get_wrho_OF_basis(basis_set)
            DO ib=1,nb
              CVecB(ib) = sum( CVecG(:) * cmplx(basis_set%dnRGB%d0(:,ib),kind=Rkind) )
            END DO
            !CVecB(1:nb) = matmul(CVecG(:),basis_set%dnRGB%d0(:,1:nb))
          END IF
        ELSE
          IF (basis_set%nb_basis == 0 ) STOP ' ERROR with packed!!!'

          SELECT CASE (basis_set%SparseGrid_type)
          CASE (0) ! Direct product
            nnq  = nq
            nnb  = 1
            nb2  = 1
            nbb  = 1
            nq2  = 1

            CALL alloc_NParray(CTempB,(/ nnq,nbb /),"CTempB",name_sub)
            CTempB(:,:) = reshape(CVecG,shape=(/ nnq,nbb /))

            DO ibasis=1,basis_set%nb_basis
              nb2 = basis_set%tab_Pbasis(ibasis)%Pbasis%nb
              nq2 = get_nq_FROM_basis(basis_set%tab_Pbasis(ibasis)%Pbasis)

              nnq = nnq / nq2

              CALL alloc_NParray(CTempG,(/ nnq,nq2,nnb /),"CTempG",name_sub)

              CTempG(:,:,:) = reshape(CTempB,shape=(/ nnq,nq2,nnb /))

              nbb = sum(basis_set%Tab_OF_Tabnb2(ibasis)%vec)

              !write(out_unitp,*) 'G=>B ibasis,nnb,Tab_OF_Tabnb2',ibasis,nnb,      &
              !                          basis_set%Tab_OF_Tabnb2(ibasis)%vec

              CALL dealloc_NParray(CTempB,"CTempB",name_sub)
              CALL alloc_NParray(CTempB,(/ nnq,nbb /),"CTempB",name_sub)


              CALL alloc_NParray(CG,(/ nq2 /),"CG",name_sub)
              CALL alloc_NParray(CB,(/ nb2 /),"CB",name_sub)

              IF (basis_set%tab_Pbasis(ibasis)%Pbasis%packed) THEN

                IF (basis_set%tab_Pbasis(ibasis)%Pbasis%cplx) THEN
                  CALL alloc_NParray(d0cb,(/ nq2,nb2 /),"d0cb",name_sub)
                  CALL Get_MatdnCGB(basis_set%tab_Pbasis(ibasis)%Pbasis,d0cb,der00)

                  w = get_wrho_OF_basis(basis_set%tab_Pbasis(ibasis)%Pbasis)


                  ibb1 = 1
                  ibb2 = 0
                  DO ib=1,nnb
                    newnb2 = basis_set%Tab_OF_Tabnb2(ibasis)%vec(ib)
                    ibb2 = ibb2 + newnb2

                    DO iq=1,nnq

                      CG(:) = CTempG(iq,:,ib)*cmplx(w,kind=Rkind)
                      CTempB(iq,ibb1:ibb2) = matmul(CG(:),d0cb(:,1:newnb2))

                    END DO

                    ibb1 = ibb1 + newnb2
                  END DO

                  CALL dealloc_NParray(d0cb,"d0cb",name_sub)
                  CALL dealloc_NParray(w,"w",name_sub)
                ELSE
                  CALL alloc_NParray(d0b,(/ nq2,nb2 /),"d0cb",name_sub)
                  CALL Get_MatdnRGB(basis_set%tab_Pbasis(ibasis)%Pbasis,d0b,der00)

                  w = get_wrho_OF_basis(basis_set%tab_Pbasis(ibasis)%Pbasis)

                  ibb1 = 1
                  ibb2 = 0
                  DO ib=1,nnb
                    newnb2 = basis_set%Tab_OF_Tabnb2(ibasis)%vec(ib)
                    ibb2 = ibb2 + newnb2

                    DO iq=1,nnq

                      CG(:) = CTempG(iq,:,ib)*cmplx(w,kind=Rkind)
                      CTempB(iq,ibb1:ibb2) = matmul(CG(:),d0b(:,1:newnb2))

                    END DO

                    ibb1 = ibb1 + newnb2
                  END DO

                  CALL dealloc_NParray(d0b,"d0cb",name_sub)
                  CALL dealloc_NParray(w,"w",name_sub)
                END IF

              ELSE

                ibb1 = 1
                ibb2 = 0
                DO ib=1,nnb
                  newnb2 = basis_set%Tab_OF_Tabnb2(ibasis)%vec(ib)
                  ibb2 = ibb2 + newnb2

                  DO iq=1,nnq

                    CG(:) = CTempG(iq,:,ib)
                    CALL RecCVecG_TO_CVecB(CG,CB,nq2,newnb2,          &
                                          basis_set%tab_Pbasis(ibasis)%Pbasis)

                    CTempB(iq,ibb1:ibb2) = CB(1:newnb2)

                  END DO

                  ibb1 = ibb1 + newnb2
                END DO
              END IF

              CALL dealloc_NParray(CB,"CB",name_sub)
              CALL dealloc_NParray(CG,"CG",name_sub)
              CALL dealloc_NParray(CTempG,"CTempG",name_sub)

              nnb = nbb
            END DO

            CVecB(:) = reshape(CTempB, shape=(/ nnq*nnb /) )
            CALL dealloc_NParray(CTempB,"CTempB",name_sub)
          CASE (1) ! Sparse basis (Smolyak 1st implementation)
            CVecB(:) = CZERO
            CALL alloc_NParray(CB,(/ nb /),"CB",name_sub)
            iq0_SG = 1
            iq1_SG = 0
            DO i_SG=1,basis_set%nb_SG
              nq_SG = get_nq_FROM_basis(basis_set%tab_PbasisSG(i_SG)%Pbasis)
              iq1_SG = iq1_SG + nq_SG
              CALL RecCVecG_TO_CVecB(CVecG(iq0_SG:iq1_SG),CB,           &
                                           nq_SG,nb,                    &
                                           basis_set%tab_PbasisSG(i_SG)%Pbasis)
              iq0_SG = iq0_SG + nq_SG
              CVecB(:) = CVecB(:) + CB(:) * cmplx(basis_set%WeightSG(i_SG),kind=Rkind)
            END DO
            CALL dealloc_NParray(CB,"CB",name_sub)

          CASE (2) ! Sparse basis (Smolyak 2d  implementation)

            CALL alloc_NParray(RVecB,shape(CVecB),'RVecB',name_sub)
            CALL alloc_NParray(RVecG,shape(CVecG),'RVecG',name_sub)

            RVecG(:) = real(CVecG,kind=Rkind)

            CALL sub_G_TO_B(RVecG,RVecB,basis_set%WeightSG,             &
                    basis_set%para_SGType2%nDind_SmolyakRep%Tab_DInd, &
                                     basis_set%nDindB%Tab_DInd,         &
                                     basis_set%tab_basisPrimSG,         &
                                     D=basis_set%nb_basis,              &
                                     LG=basis_set%L_SparseGrid,         &
                                     LB=basis_set%L_SparseBasis)
            CVecB(:) = cmplx(RVecB,kind=Rkind)

            RVecG(:) = aimag(CVecG)
            CALL sub_G_TO_B(RVecG,RVecB,basis_set%WeightSG,             &
                    basis_set%para_SGType2%nDind_SmolyakRep%Tab_DInd, &
                                     basis_set%nDindB%Tab_DInd,         &
                                     basis_set%tab_basisPrimSG,         &
                                     D=basis_set%nb_basis,              &
                                     LG=basis_set%L_SparseGrid,         &
                                     LB=basis_set%L_SparseBasis)
            CVecB(:) = CVecB + EYE*cmplx(RVecB,kind=Rkind)

            CALL dealloc_NParray(RVecB,'RVecB',name_sub)
            CALL dealloc_NParray(RVecG,'RVecG',name_sub)


            !write(out_unitp,*) ' ERROR in ',name_sub
            !STOP 'SparseGrid_type=2'

          CASE (4) ! Sparse basis (Smolyak 4th implementation)
            CALL alloc_NParray(RVecB,shape(CVecB),'RVecB',name_sub)

            !!  RVecG TO SRep
            CALL alloc2_SmolyakRep(SRep,basis_set%para_SGType2%nDind_SmolyakRep,&
                                   basis_set%tab_basisPrimSG,grid=.TRUE.,       &
                                   nb0=basis_set%para_SGType2%nb0)
            CALL tabR2bis_TO_SmolyakRep1(SRep,real(CVecG,kind=Rkind)) ! on the grid

            IF (allocated(basis_set%para_SGType2%nDind_SmolyakRep%Tab_nDval)) THEN
              CALL GSmolyakRep_TO_BSmolyakRep(SRep,                        &
                        basis_set%para_SGType2%nDind_SmolyakRep%Tab_nDval, &
                        basis_set%tab_basisPrimSG,basis_set%para_SGType2%nb0)
            ELSE
              CALL GSmolyakRep_TO3_BSmolyakRep(SRep,basis_set%para_SGType2,&
                                               basis_set%tab_basisPrimSG)
            END IF

            CALL SmolyakRepBasis_TO_tabPackedBasis(SRep,RVecB,          &
                             basis_set%tab_basisPrimSG,basis_set%nDindB,&
                             basis_set%para_SGType2,basis_set%WeightSG)

            CALL dealloc_SmolyakRep(SRep)

            CVecB(:) = cmplx(RVecB,kind=Rkind)


            !!  RVecG TO SRep
            CALL alloc2_SmolyakRep(SRep,basis_set%para_SGType2%nDind_SmolyakRep,&
                                   basis_set%tab_basisPrimSG,grid=.TRUE.,       &
                                   nb0=basis_set%para_SGType2%nb0)
            CALL tabR2bis_TO_SmolyakRep1(SRep,aimag(CVecG)) ! on the grid

            IF (allocated(basis_set%para_SGType2%nDind_SmolyakRep%Tab_nDval)) THEN
              CALL GSmolyakRep_TO_BSmolyakRep(SRep,                     &
                    basis_set%para_SGType2%nDind_SmolyakRep%Tab_nDval,  &
                    basis_set%tab_basisPrimSG,basis_set%para_SGType2%nb0)
            ELSE
              CALL GSmolyakRep_TO3_BSmolyakRep(SRep,basis_set%para_SGType2,&
                                               basis_set%tab_basisPrimSG)
            END IF

            CALL SmolyakRepBasis_TO_tabPackedBasis(SRep,RVecB,          &
                             basis_set%tab_basisPrimSG,basis_set%nDindB,&
                             basis_set%para_SGType2,basis_set%WeightSG)

            CALL dealloc_SmolyakRep(SRep)


            CVecB(:) = CVecB + EYE*cmplx(RVecB,kind=Rkind)

            CALL dealloc_NParray(RVecB,'RVecB',name_sub)

          CASE DEFAULT
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' WRONG SparseGrid_type',basis_set%SparseGrid_type
            write(out_unitp,*) ' The possibilities are: 0, 1, 2, 4'
            STOP
          END SELECT
        END IF

      IF (debug) write(out_unitp,*) ' END in ',name_sub


      END SUBROUTINE RecCVecG_TO_CVecB

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      RECURSIVE SUBROUTINE RecRvecB_TO_RVecG(RvecB,RVecG,nb,nq,basis_set,tab_der)
        USE mod_basis_BtoG_GtoB_SGType4
        IMPLICIT NONE
        TYPE (basis), intent(in)         :: basis_set
        integer, intent(in)              :: nq,nb
        real (kind=Rkind), intent(inout) :: RVecG(:)
        real (kind=Rkind), intent(in)    :: RvecB(:)
        integer, optional                :: tab_der(2)

        integer :: tab_der_loc(2),dnba_ind(2),iQact_ba,k
        real (kind=Rkind), allocatable       :: RTempG(:,:,:)
        real (kind=Rkind), allocatable       :: RTempB(:,:)
        TYPE(Type_SmolyakRep)                :: SRep ! smolyak rep for SparseGrid_type=4

        real (kind=Rkind), allocatable       :: RG(:)
        real (kind=Rkind), allocatable       :: RB(:)
        real (kind=Rkind), allocatable       :: dnb(:,:)

        integer                          :: nbb,ibb1,ibb2
        integer                          :: nnb,nb2,ib,ib2,newnb2
        integer                          :: nnq,nq2,iq,iq2
        integer                          :: ibasis
        integer                          :: i_SG,iq0_SG,iq1_SG,nq_SG,ii_SG,nq0_SG

        integer                          :: itabR,iG,nR


        integer                          :: nb_thread


!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='RecRvecB_TO_RVecG'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!---------------------------------------------------------------------

        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
          write(out_unitp,*) 'RvecB(:)',RvecB(:)
          !CALL RecWrite_basis(basis_set)
          CALL flush_perso(out_unitp)
        END IF

        IF (present(tab_der)) THEN
          tab_der_loc(:) = tab_der(:)
        ELSE
          tab_der_loc(:) = 0
        END IF
        WHERE (tab_der_loc < 0) tab_der_loc = 0

        IF (basis_set%cplx) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' the basis is complex but the'
          write(out_unitp,*) ' the vector is real !!'
          STOP
        END IF

        IF (BasisTOGrid_omp == 0) THEN
          nb_thread = 1
        ELSE
          nb_thread = BasisTOGrid_maxth
        END IF
        !write(out_unitp,*) 'nb_thread in ',name_sub,' : ',nb_thread

        IF (basis_set%packed_done) THEN ! Packed basis
          nb_mult_BTOG = nb_mult_BTOG + int(nb,kind=ILkind)*size(RVecG,kind=ILkind)

          dnba_ind(:) = basis_set%Tabder_Qdyn_TO_Qbasis(tab_der_loc(:))
          CALL flush_perso(out_unitp)
          IF (dnba_ind(1) == 0 .AND. dnba_ind(2) == 0) THEN ! dnba_ind(:)=0 => no derivative

            RVecG(:) = matmul(basis_set%dnRGB%d0,RvecB)
          ELSE IF (dnba_ind(1) == 0) THEN ! first derivative
              IF (basis_set%dnBBRep_done) THEN
                CALL alloc_NParray(RB, (/ nb /),"RB",name_sub)
                RB(:) = matmul(basis_set%dnRBB%d1(1:nb,1:nb,dnba_ind(2)),RvecB(1:nb))
                DO iq=1,nq
                  RVecG(iq) = sum(basis_set%dnRGB%d0(iq,1:nb)*RB)
                END DO
!                RVecG(:) = matmul(basis_set%dnRGB%d0,RB)
                CALL dealloc_NParray(RB,"RB",name_sub)
              ELSE
                DO iq=1,nq
                  RVecG(iq) = sum(basis_set%dnRGB%d1(iq,1:nb,dnba_ind(2))*RvecB(1:nb))
                END DO
!                RVecG(:) = matmul(basis_set%dnRGB%d1(:,1:nb,dnba_ind(2)),RvecB(1:nb))
              END IF
          ELSE IF (dnba_ind(2) == 0) THEN ! first derivative
            IF (basis_set%dnBBRep_done) THEN
              CALL alloc_NParray(RB,(/ nb /),"RB",name_sub)
              RB(:) = matmul(basis_set%dnRBB%d1(1:nb,1:nb,dnba_ind(1)),RvecB(1:nb))
              DO iq=1,nq
                RVecG(iq) = sum(basis_set%dnRGB%d0(iq,1:nb)*RB)
              END DO
!              RVecG(:) = matmul(basis_set%dnRGB%d0,RB)
              CALL dealloc_NParray(RB,"RB",name_sub)
            ELSE
              DO iq=1,nq
                RVecG(iq) = sum(basis_set%dnRGB%d1(iq,1:nb,dnba_ind(1))*RvecB(1:nb))
              END DO
!              RVecG(:) = matmul(basis_set%dnRGB%d1(:,1:nb,dnba_ind(1)),RvecB(1:nb))
            END IF
          ELSE
            IF (basis_set%dnBBRep_done) THEN
              CALL alloc_NParray(RB,(/ nb /),"RB",name_sub)
              RB(:) = matmul(basis_set%dnRBB%d2(1:nb,1:nb,dnba_ind(1),dnba_ind(2)),RvecB(1:nb))
              DO iq=1,nq
                RVecG(iq) = sum(basis_set%dnRGB%d0(iq,1:nb)*RB)
              END DO
!              RVecG(:) = matmul(basis_set%dnRGB%d0,RB)
              CALL dealloc_NParray(RB,"RB",name_sub)
            ELSE
              DO iq=1,nq
                RVecG(iq) = sum(basis_set%dnRGB%d2(iq,1:nb,dnba_ind(1),dnba_ind(2))*RvecB(1:nb))
              END DO
!              RVecG(:) = matmul(basis_set%dnRGB%d2(:,1:nb,dnba_ind(1),dnba_ind(2)),RvecB(1:nb))
            END IF
          END IF
        ELSE ! basis_set%nb_basis MUST BE > 0
          IF (basis_set%nb_basis == 0 ) STOP ' ERROR with packed_done !!!'

          SELECT CASE (basis_set%SparseGrid_type)
          CASE (0) ! Direct product
            nnb = nb
            nbb = nb
            nnq = 1
            nb2 = 1
            nq2 = 1

            CALL alloc_NParray(RTempG,(/ nnq,nq2,nnb /),"RTempG",name_sub)

            RTempG(:,:,:) = reshape(RvecB,shape=(/ nnq,nq2,nnb /))

            DO ibasis=basis_set%nb_basis,1,-1

              nb2 = basis_set%tab_Pbasis(ibasis)%Pbasis%nb
              nq2 = get_nq_FROM_basis(basis_set%tab_Pbasis(ibasis)%Pbasis)
              nnb = basis_set%Tab_OF_Tabnb2(ibasis)%nb_var_vec

              CALL alloc_NParray(RTempB, (/nnq,nbb /),"RTempB",name_sub)

              RTempB(:,:) = reshape(RTempG,shape=(/ nnq,nbb /))

              CALL dealloc_NParray(RTempG,"RTempG",name_sub)
              CALL alloc_NParray(RTempG,(/ nnq,nq2,nnb /),"RTempG",name_sub)

              IF (basis_set%tab_Pbasis(ibasis)%Pbasis%packed) THEN
                dnba_ind(:) = basis_set%tab_Pbasis(ibasis)%Pbasis%Tabder_Qdyn_TO_Qbasis(tab_der_loc(:))

                CALL alloc_NParray(dnb,(/ nq2,nb2 /),"dnb",name_sub)
                CALL Get_MatdnRGB(basis_set%tab_Pbasis(ibasis)%Pbasis,dnb,dnba_ind)

                ibb1 = 1
                ibb2 = 0
                DO ib=1,nnb
                  newnb2 = basis_set%Tab_OF_Tabnb2(ibasis)%vec(ib)
                  ibb2 = ibb2 + newnb2
                  DO iq=1,nnq
                    RTempG(iq,:,ib) = matmul(dnb(:,1:newnb2),RTempB(iq,ibb1:ibb2))

                    nb_mult_BTOG = nb_mult_BTOG + int(newnb2,kind=ILkind)*size(dnb(:,1),kind=ILkind)
                  END DO
                  ibb1 = ibb1 + newnb2
                END DO

                CALL dealloc_NParray(dnb,"d0b",name_sub)
              ELSE
                ibb1 = 1
                ibb2 = 0
                DO ib=1,nnb
                  newnb2 = basis_set%Tab_OF_Tabnb2(ibasis)%vec(ib)
                  ibb2 = ibb2 + newnb2

                  DO iq=1,nnq
                    CALL RecRvecB_TO_RVecG(RTempB(iq,ibb1:ibb2),RTempG(iq,:,ib),&
                                           newnb2,nq2,                          &
                                           basis_set%tab_Pbasis(ibasis)%Pbasis, &
                                           tab_der_loc)
                  END DO
                  ibb1 = ibb1 + newnb2
                END DO

              END IF

              nbb = nnb
              nnq = nnq * nq2
              CALL dealloc_NParray(RTempB,"RTempB",name_sub)
            END DO

            RvecG(:) = reshape(RTempG, shape=(/ nnq*nnb /) )
            CALL dealloc_NParray(RTempG,"RTempG",name_sub)

          CASE (1) ! Sparse basis (Smolyak 1st implementation)
            iq0_SG = 1
            iq1_SG = 0
            DO i_SG=1,basis_set%nb_SG
               nq_SG = get_nq_FROM_basis(basis_set%tab_PbasisSG(i_SG)%Pbasis)
               !write(out_unitp,*) 'i_SG,nq_SG',i_SG,nq_SG
!               iq1_SG = iq1_SG + nq_SG
               iq0_SG = 1
               DO ii_SG=1,i_SG-1
                 nq0_SG = get_nq_FROM_basis(basis_set%tab_PbasisSG(ii_SG)%Pbasis)
                 iq0_SG = iq0_SG + nq0_SG
               END DO
               iq1_SG = iq0_SG-1 + nq_SG
               !write(out_unitp,*) 'iq0_SG,iq1_SG',iq0_SG,iq1_SG
               CALL RecRvecB_TO_RVecG(RvecB,RVecG(iq0_SG:iq1_SG),       &
                                      nb,nq_SG,                         &
                                     basis_set%tab_PbasisSG(i_SG)%Pbasis,&
                                      tab_der_loc)
!               iq0_SG = iq0_SG + nq_SG
            END DO

          CASE (2) ! Sparse basis (Smolyak 2d  implementation)

            CALL sub_B_TO_G(RVecB,RVecG,                                &
                    basis_set%para_SGType2%nDind_SmolyakRep%Tab_DInd, &
                                     basis_set%nDindB%Tab_DInd,         &
                                     basis_set%tab_basisPrimSG,         &
                                     D=basis_set%nb_basis,              &
                                     LG=basis_set%L_SparseGrid,         &
                                     LB=basis_set%L_SparseBasis)

          CASE (4) ! Sparse basis (Smolyak 4th implementation)

            CALL tabPackedBasis_TO_SmolyakRepBasis(SRep,RVecB,basis_set%tab_basisPrimSG,basis_set%nDindB,basis_set%para_SGType2)

            IF (allocated(basis_set%para_SGType2%nDind_SmolyakRep%Tab_nDval)) THEN
              CALL BSmolyakRep_TO_GSmolyakRep(SRep,                         &
                 basis_set%para_SGType2%nDind_SmolyakRep%Tab_nDval,         &
                 basis_set%tab_basisPrimSG,basis_set%para_SGType2%nb0)
            ELSE
              CALL BSmolyakRep_TO3_GSmolyakRep(SRep,basis_set%para_SGType2, &
                                               basis_set%tab_basisPrimSG)
            END IF

            CALL SmolyakRep2_TO_tabR1bis(RVecG,SRep)
            CALL dealloc_SmolyakRep(SRep)


          CASE DEFAULT
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' WRONG SparseGrid_type',basis_set%SparseGrid_type
            write(out_unitp,*) ' The possibilities are: 0, 1, 2, 4'
            STOP
          END SELECT

        END IF

       IF (debug) THEN
         write(out_unitp,*) 'RvecG(:)',RvecG(:)
         write(out_unitp,*) 'END ',name_sub
       END IF

      END SUBROUTINE RecRvecB_TO_RVecG

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      RECURSIVE SUBROUTINE RecCVecB_TO_CVecG(CVecB,CVecG,nb,nq,basis_set,tab_der)
        USE mod_basis_BtoG_GtoB_SGType4
        IMPLICIT NONE
        TYPE (basis), intent(in)            :: basis_set
        integer, intent(in)                 :: nq,nb
        complex (kind=Rkind), intent(inout) :: CVecG(:)
        complex (kind=Rkind), intent(in)    :: CVecB(:)
        integer, optional                   :: tab_der(2)

        integer :: tab_der_loc(2),dnba_ind(2),iQact_ba,k
        complex (kind=Rkind), allocatable   :: CTempG(:,:,:)
        complex (kind=Rkind), allocatable   :: CTempB(:,:)
        TYPE(Type_SmolyakRep)               :: SRep ! smolyak rep for SparseGrid_type=4

        complex (kind=Rkind), allocatable   :: CG(:)
        complex (kind=Rkind), allocatable   :: CB(:)
        real (kind=Rkind), allocatable      :: dnb(:,:)
        complex (kind=Rkind), allocatable   :: dncb(:,:)


        real (kind=Rkind), allocatable :: RVecG(:)
        real (kind=Rkind), allocatable :: RVecB(:)

        integer                             :: nbb,ibb1,ibb2
        integer                             :: nnb,nb2,ib,ib2,newnb2
        integer                             :: nnq,nq2,iq,iq2
        integer                             :: ibasis
        integer                             :: i_SG,iq0_SG,iq1_SG,nq_SG

        integer                             :: itabR,iG,nR

        integer                             :: nb_thread

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='RecCVecB_TO_CVecG'
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------

        IF (present(tab_der)) THEN
          tab_der_loc(:) = tab_der(:)
        ELSE
          tab_der_loc(:) = 0
        END IF
        WHERE (tab_der_loc < 0) tab_der_loc = 0

        IF (basis_set%cplx) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' the basis is complex but the'
          write(out_unitp,*) ' the vector is real !!'
          STOP
        END IF


        !IF (BasisTOGrid_omp == 0) THEN
        !  nb_thread = 1
        !ELSE
        !  nb_thread = BasisTOGrid_maxth
        !END IF
        !write(out_unitp,*) 'nb_thread in ',name_sub,' : ',nb_thread

        IF (basis_set%cplx .AND. basis_set%packed_done) THEN

          dnba_ind(:) = basis_set%Tabder_Qdyn_TO_Qbasis(tab_der_loc(:))

          IF (dnba_ind(1) == 0 .AND. dnba_ind(2) == 0) THEN ! dnba_ind(:)=0 => no derivative
            CVecG(:) = matmul(basis_set%dnCGB%d0,CVecB)
          ELSE IF (dnba_ind(1) == 0) THEN ! first derivative
            IF (basis_set%dnBBRep_done) THEN
              CALL alloc_NParray(CB,(/ nb /),"CB",name_sub)
              CB(:) = matmul(basis_set%dnCBB%d1(:,1:nb,dnba_ind(2)),CVecB(1:nb))
              CVecG(:) = matmul(basis_set%dnCGB%d0(:,1:nb),CB)
              CALL dealloc_NParray(CB,"CB",name_sub)
            ELSE
              CVecG(:) = matmul(basis_set%dnCGB%d1(:,1:nb,dnba_ind(2)),CVecB(1:nb))
            END IF
          ELSE IF (dnba_ind(2) == 0) THEN ! first derivative
            IF (basis_set%dnBBRep_done) THEN
              CALL alloc_NParray(CB,(/ nb /),"CB",name_sub)
              CB(:) = matmul(basis_set%dnCBB%d1(:,1:nb,dnba_ind(1)),CVecB(1:nb))
              CVecG(:) = matmul(basis_set%dnCGB%d0(:,1:nb),CB)
              CALL dealloc_NParray(CB,"CB",name_sub)
            ELSE
              CVecG(:) = matmul(basis_set%dnCGB%d1(:,1:nb,dnba_ind(1)),CVecB(1:nb))
            END IF
          ELSE ! second derivative
            IF (basis_set%dnBBRep_done) THEN
              CALL alloc_NParray(CB,(/ nb /),"CB",name_sub)
              CB(:) = matmul(basis_set%dnCBB%d2(:,1:nb,dnba_ind(1),dnba_ind(2)),CVecB(1:nb))
              CVecG(:) = matmul(basis_set%dnCGB%d0(:,1:nb),CB)
              CALL dealloc_NParray(CB,"CB",name_sub)
            ELSE
              CVecG(:) = matmul(basis_set%dnCGB%d2(:,1:nb,dnba_ind(1),dnba_ind(2)),CVecB(1:nb))
            END IF
          END IF

        ELSE IF (.NOT. basis_set%cplx .AND. basis_set%packed_done) THEN

          dnba_ind(:) = basis_set%Tabder_Qdyn_TO_Qbasis(tab_der_loc(:))

          IF (dnba_ind(1) == 0 .AND. dnba_ind(2) == 0) THEN ! dnba_ind(:)=0 => no derivative
            CVecG(:) = matmul(basis_set%dnRGB%d0,CVecB)
          ELSE IF (dnba_ind(1) == 0) THEN ! first derivative
            IF (basis_set%dnBBRep_done) THEN
              CALL alloc_NParray(CB,(/ nb /),"CB",name_sub)
              CB(:) = matmul(basis_set%dnRBB%d1(:,1:nb,dnba_ind(2)),CVecB(1:nb))
              CVecG(:) = matmul(basis_set%dnRGB%d0(:,1:nb),CB)
              CALL dealloc_NParray(CB,"CB",name_sub)
            ELSE
              CVecG(:) = matmul(basis_set%dnRGB%d1(:,1:nb,dnba_ind(2)),CVecB(1:nb))
            END IF
          ELSE IF (dnba_ind(2) == 0) THEN ! first derivative
            IF (basis_set%dnBBRep_done) THEN
              CALL alloc_NParray(CB,(/ nb /),"CB",name_sub)
              CB(:) = matmul(basis_set%dnRBB%d1(:,1:nb,dnba_ind(1)),CVecB(1:nb))
              CVecG(:) = matmul(basis_set%dnRGB%d0(:,1:nb),CB)
              CALL dealloc_NParray(CB,"CB",name_sub)
            ELSE
              CVecG(:) = matmul(basis_set%dnRGB%d1(:,1:nb,dnba_ind(1)),CVecB(1:nb))
            END IF
          ELSE ! second derivative
            IF (basis_set%dnBBRep_done) THEN
              CALL alloc_NParray(CB,(/ nb /),"CB",name_sub)
              CB(:) = matmul(basis_set%dnRBB%d2(:,1:nb,dnba_ind(1),dnba_ind(2)),CVecB(1:nb))
              CVecG(:) = matmul(basis_set%dnRGB%d0(:,1:nb),CB)
              CALL dealloc_NParray(CB,"CB",name_sub)
            ELSE
              CVecG(:) = matmul(basis_set%dnRGB%d2(:,1:nb,dnba_ind(1),dnba_ind(2)),CVecB(1:nb))
            END IF
          END IF

        ELSE ! basis_set%nb_basis MUST BE > 0
          IF (basis_set%nb_basis == 0 ) STOP ' ERROR with packed!!!'

          SELECT CASE (basis_set%SparseGrid_type)
          CASE (0) ! Direct product

            nnb = nb
            nbb = nb
            nnq = 1
            nb2 = 1
            nq2 = 1

            CALL alloc_NParray(CTempG,(/ nnq,nq2,nnb /),"CTempG",name_sub)
            CTempG(:,:,:) = reshape(CVecB,shape=(/ nnq,nq2,nnb /))

            DO ibasis=basis_set%nb_basis,1,-1

              nb2 = basis_set%tab_Pbasis(ibasis)%Pbasis%nb
              nq2 = get_nq_FROM_basis(basis_set%tab_Pbasis(ibasis)%Pbasis)
              nnb = basis_set%Tab_OF_Tabnb2(ibasis)%nb_var_vec

!             write(out_unitp,*) 'B=>G ibasis,Tab_OF_Tabnb2',ibasis,              &
!                                        basis_set%Tab_OF_Tabnb2(ibasis)%vec

              CALL alloc_NParray(CTempB,(/ nnq,nbb /),"CTempB",name_sub)

              CTempB(:,:) = reshape(CTempG,shape=(/ nnq,nbb /))
              CALL dealloc_NParray(CTempG,"CTempG",name_sub)
              CALL alloc_NParray(CTempG,(/ nnq,nq2,nnb /),"CTempG",name_sub)

              IF (basis_set%tab_Pbasis(ibasis)%Pbasis%packed) THEN
                IF (basis_set%tab_Pbasis(ibasis)%Pbasis%cplx) THEN
                  dnba_ind(:) = basis_set%tab_Pbasis(ibasis)%Pbasis%Tabder_Qdyn_TO_Qbasis(tab_der_loc(:))

                  CALL alloc_NParray(dncb,(/ nq2,nb2 /),"dncb",name_sub)
                  CALL Get_MatdnCGB(basis_set%tab_Pbasis(ibasis)%Pbasis,dncb,dnba_ind)

                  ibb1 = 1
                  ibb2 = 0
                  DO ib=1,nnb
                    newnb2 = basis_set%Tab_OF_Tabnb2(ibasis)%vec(ib)
                    ibb2 = ibb2 + newnb2
                    DO iq=1,nnq
                      CTempG(iq,:,ib) = matmul(dnb(:,1:newnb2),CTempB(iq,ibb1:ibb2))
                    END DO
                    ibb1 = ibb1 + newnb2
                  END DO

                  CALL dealloc_NParray(dnb,"d0b",name_sub)
                ELSE
                  dnba_ind(:) = basis_set%tab_Pbasis(ibasis)%Pbasis%Tabder_Qdyn_TO_Qbasis(tab_der_loc(:))

                  CALL alloc_NParray(dnb,(/ nq2,nb2 /),"dnb",name_sub)
                  CALL Get_MatdnRGB(basis_set%tab_Pbasis(ibasis)%Pbasis,dnb,dnba_ind)

                  ibb1 = 1
                  ibb2 = 0
                  DO ib=1,nnb
                    newnb2 = basis_set%Tab_OF_Tabnb2(ibasis)%vec(ib)
                    ibb2 = ibb2 + newnb2
                    DO iq=1,nnq
                      CTempG(iq,:,ib) = matmul(dnb(:,1:newnb2),CTempB(iq,ibb1:ibb2))
                    END DO
                    ibb1 = ibb1 + newnb2
                  END DO

                  CALL dealloc_NParray(dnb,"d0b",name_sub)
                END IF
              ELSE
                ibb1 = 1
                ibb2 = 0
                DO ib=1,nnb
                  newnb2 = basis_set%Tab_OF_Tabnb2(ibasis)%vec(ib)
                  ibb2 = ibb2 + newnb2
                  DO iq=1,nnq
                    CALL RecCVecB_TO_CVecG(CTempB(iq,ibb1:ibb2),        &
                                           CTempG(iq,:,ib),newnb2,nq2,  &
                                   basis_set%tab_Pbasis(ibasis)%Pbasis, &
                                                             tab_der_loc)
                  END DO
                  ibb1 = ibb1 + newnb2
                END DO
              END IF

              nbb = nnb
              nnq = nnq * nq2
              CALL dealloc_NParray(CTempB,'CTempB',name_sub)
            END DO

            CvecG(:) = reshape(CTempG, shape=(/ nnq*nnb /) )
            CALL dealloc_NParray(CTempG,"CTempG",name_sub)

          CASE (1) ! Sparse basis (Smolyak 1st implementation)
            iq0_SG = 1
            iq1_SG = 0
            DO i_SG=1,basis_set%nb_SG
               nq_SG  = get_nq_FROM_basis(basis_set%tab_PbasisSG(i_SG)%Pbasis)
               iq1_SG = iq1_SG + nq_SG
               CALL RecCVecB_TO_CVecG(CVecB,CVecG(iq0_SG:iq1_SG),       &
                                      nb,nq_SG,                         &
                                   basis_set%tab_PbasisSG(i_SG)%Pbasis, &
                                      tab_der_loc)
               iq0_SG = iq0_SG + nq_SG
            END DO

          CASE (2) ! Sparse basis (Smolyak 2d  implementation)
            CALL alloc_NParray(RVecB,shape(CVecB),'RVecB',name_sub)
            CALL alloc_NParray(RVecG,shape(CVecG),'RVecG',name_sub)

            RVecB(:) = real(CVecB,kind=Rkind)
            CALL sub_B_TO_G(RVecB,RVecG,                                &
                    basis_set%para_SGType2%nDind_SmolyakRep%Tab_DInd, &
                                     basis_set%nDindB%Tab_DInd,         &
                                     basis_set%tab_basisPrimSG,         &
                                     D=basis_set%nb_basis,              &
                                     LG=basis_set%L_SparseGrid,         &
                                     LB=basis_set%L_SparseBasis)
            CVecG(:) = cmplx(RVecG,kind=Rkind)

            RVecB(:) = aimag(CVecB)
            CALL sub_B_TO_G(RVecB,RVecG,                                &
                    basis_set%para_SGType2%nDind_SmolyakRep%Tab_DInd, &
                                     basis_set%nDindB%Tab_DInd,         &
                                     basis_set%tab_basisPrimSG,         &
                                     D=basis_set%nb_basis,              &
                                     LG=basis_set%L_SparseGrid,         &
                                     LB=basis_set%L_SparseBasis)
            CVecG(:) = CVecG + EYE*cmplx(RVecG,kind=Rkind)

            CALL dealloc_NParray(RVecB,'RVecB',name_sub)
            CALL dealloc_NParray(RVecG,'RVecG',name_sub)
            !write(out_unitp,*) ' ERROR in ',name_sub
            !STOP 'SparseGrid_type=2'

          CASE (4) ! Sparse basis (Smolyak 4th implementation)
            CALL alloc_NParray(RVecB,shape(CVecB),'RVecB',name_sub)
            !CALL alloc_NParray(RVecG,shape(CVecG),'RVecG',name_sub)

            RVecB(:) = real(CVecB,kind=Rkind)
            CALL tabPackedBasis_TO_SmolyakRepBasis(SRep,RVecB,basis_set%tab_basisPrimSG,basis_set%nDindB,basis_set%para_SGType2)

            IF (allocated(basis_set%para_SGType2%nDind_SmolyakRep%Tab_nDval)) THEN
              CALL BSmolyakRep_TO_GSmolyakRep(SRep,                         &
                 basis_set%para_SGType2%nDind_SmolyakRep%Tab_nDval,         &
                 basis_set%tab_basisPrimSG,basis_set%para_SGType2%nb0)
            ELSE
              CALL BSmolyakRep_TO3_GSmolyakRep(SRep,basis_set%para_SGType2, &
                                               basis_set%tab_basisPrimSG)
            END IF

            itabR = 0
            DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)
              nR = size(SRep%SmolyakRep(iG)%V)
              CVecG(itabR+1:itabR+nR) = cmplx(SRep%SmolyakRep(iG)%V,kind=Rkind)
              itabR = itabR + nR
            END DO
            CALL dealloc_SmolyakRep(SRep)


            RVecB(:) = aimag(CVecB)
            CALL tabPackedBasis_TO_SmolyakRepBasis(SRep,RVecB,basis_set%tab_basisPrimSG,basis_set%nDindB,basis_set%para_SGType2)

            IF (allocated(basis_set%para_SGType2%nDind_SmolyakRep%Tab_nDval)) THEN
              CALL BSmolyakRep_TO_GSmolyakRep(SRep,                         &
                 basis_set%para_SGType2%nDind_SmolyakRep%Tab_nDval,         &
                 basis_set%tab_basisPrimSG,basis_set%para_SGType2%nb0)
            ELSE
              CALL BSmolyakRep_TO3_GSmolyakRep(SRep,basis_set%para_SGType2, &
                                               basis_set%tab_basisPrimSG)
            END IF

            itabR = 0
            DO iG=lbound(SRep%SmolyakRep,dim=1),ubound(SRep%SmolyakRep,dim=1)
              nR = size(SRep%SmolyakRep(iG)%V)
              CVecG(itabR+1:itabR+nR) = CVecG(itabR+1:itabR+nR) +       &
                             EYE*cmplx(SRep%SmolyakRep(iG)%V,kind=Rkind)
              itabR = itabR + nR
            END DO
            CALL dealloc_SmolyakRep(SRep)

            CALL dealloc_NParray(RVecB,'RVecB',name_sub)
            !CALL dealloc_NParray(RVecG,'RVecG',name_sub)

          CASE DEFAULT
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' WRONG SparseGrid_type',basis_set%SparseGrid_type
            write(out_unitp,*) ' The possibilities are: 0, 1, 2, 4'
            STOP
          END SELECT


          IF (basis_set%nb_SG > 0) THEN ! for Sparse Grid

          ELSE ! for direct-product grid

          END IF
        END IF

      END SUBROUTINE RecCVecB_TO_CVecG

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO

      RECURSIVE SUBROUTINE DerivOp_TO_RVecG(RVecG,nq,basis_set,tab_der)
        USE mod_basis_BtoG_GtoB_SGType4
        IMPLICIT NONE
        TYPE (basis), intent(in)         :: basis_set
        integer, intent(in)              :: nq
        real (kind=Rkind), intent(inout) :: RVecG(:)
        integer, optional                :: tab_der(2)

        integer :: tab_der_loc(2),dnba_ind(2),iQact_ba,k
        real (kind=Rkind), allocatable       :: RG1(:,:,:)
        real (kind=Rkind), allocatable       :: RG2(:,:,:)
        real (kind=Rkind), allocatable       :: BGG(:,:)

        integer                          :: nnq1,nnq3,nq2,iq1,iq2,iq3
        integer                          :: iqi,iqe
        logical                          :: skip
        integer                          :: ibasis
        integer                          :: i_SG,iq0_SG,iq1_SG,nq_SG,ii_SG,nq0_SG
        integer                          :: nb_thread
        TYPE(Type_SmolyakRep)            :: SRep ! smolyak rep for SparseGrid_type=4
        integer                          :: itabR,iG,nR

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='DerivOp_TO_RVecG'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!---------------------------------------------------------------------

        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
          write(out_unitp,*) 'RvecG(:)',RvecG(:)
          write(out_unitp,*) 'packed_done',basis_set%packed_done
          write(out_unitp,*) 'nb_basis',basis_set%nb_basis
          write(out_unitp,*) 'SparseGrid_type',basis_set%SparseGrid_type

          !CALL RecWrite_basis(basis_set)
          CALL flush_perso(out_unitp)
        END IF

        IF (present(tab_der)) THEN
          tab_der_loc(:) = tab_der(:)
        ELSE
          tab_der_loc(:) = 0
        END IF
        IF (debug) write(out_unitp,*) 'tab_der_loc',tab_der_loc
        WHERE (tab_der_loc < 0) tab_der_loc = 0

        IF (basis_set%cplx) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' the basis is complex but the'
          write(out_unitp,*) ' the vector is real !!'
          STOP
        END IF

        IF (BasisTOGrid_omp == 0) THEN
          nb_thread = 1
        ELSE
          nb_thread = BasisTOGrid_maxth
        END IF
        !write(out_unitp,*) 'nb_thread in ',name_sub,' : ',nb_thread

        IF (basis_set%packed_done) THEN ! Packed basis

          dnba_ind(:) = basis_set%Tabder_Qdyn_TO_Qbasis(tab_der_loc(:))

          CALL flush_perso(out_unitp)
          IF (dnba_ind(1) == 0 .AND. dnba_ind(2) == 0) THEN ! dnba_ind(:)=0 => no derivative
            CONTINUE ! notthing to do
          ELSE
            nq2 = get_nq_FROM_basis(basis_set)
            CALL alloc_NParray(BGG,(/ nq2,nq2 /),"BGG",name_sub)
            CALL Get_MatdnRGG(basis_set,BGG,dnba_ind)
            RVecG(:) = matmul(BGG,RVecG(:))
            CALL dealloc_NParray(BGG,"BGG",name_sub)
          END IF

        ELSE ! basis_set%nb_basis MUST BE > 0
          IF (basis_set%nb_basis == 0 ) STOP ' ERROR with packed_done !!!'

          SELECT CASE (basis_set%SparseGrid_type)
          CASE (0) ! Direct product
            nnq3 = nq
            nnq1 = 1
            nq2  = 1

            CALL alloc_NParray(RG1,(/ nnq1,nq2,nnq3 /),"RG1",name_sub)
            RG1(:,:,:) = reshape(RvecG,shape=(/ nnq1,nq2,nnq3 /))

            DO ibasis=basis_set%nb_basis,1,-1

              nq2  = get_nq_FROM_basis(basis_set%tab_Pbasis(ibasis)%Pbasis)
              nnq3 = nnq3 / nq2



              dnba_ind(:) = basis_set%tab_Pbasis(ibasis)%Pbasis%Tabder_Qdyn_TO_Qbasis(tab_der_loc(:))

              IF (dnba_ind(1) /= 0 .OR. dnba_ind(2) /= 0) THEN
                CALL alloc_NParray(RG2,(/ nnq1,nq2,nnq3 /),"RG2",name_sub)
                RG2(:,:,:) = reshape(RG1,shape=(/ nnq1,nq2,nnq3 /))

                IF (basis_set%tab_Pbasis(ibasis)%Pbasis%packed) THEN

                  CALL alloc_NParray(BGG,(/ nq2,nq2 /),"BGG",name_sub)

                  CALL Get_MatdnRGG(basis_set%tab_Pbasis(ibasis)%Pbasis,BGG,dnba_ind)

                 !$OMP parallel do default(none)                        &
                 !$OMP shared(BGG,RG2,nnq3,nnq1)                        &
                 !$OMP private(iq1,iq3)                                 &
                 !$OMP num_threads(nb_thread)
                  DO iq3=1,nnq3
                  DO iq1=1,nnq1
                     RG2(iq1,:,iq3) = matmul(BGG,RG2(iq1,:,iq3))
                  END DO
                  END DO
                 !$OMP end parallel do

                  CALL dealloc_NParray(BGG,"BGG",name_sub)

                ELSE
                  DO iq3=1,nnq3
                  DO iq1=1,nnq1
                    CALL DerivOp_TO_RVecG(RG2(iq1,:,iq3),nq2,           &
                                   basis_set%tab_Pbasis(ibasis)%Pbasis, &
                                                            tab_der_loc)
                  END DO
                  END DO
                END IF

                CALL dealloc_NParray(RG1,"RG1",name_sub)
                CALL alloc_NParray(RG1,(/ nnq1,nq2,nnq3 /),"RG1",name_sub)
                RG1 = RG2
                CALL dealloc_NParray(RG2,"RG2",name_sub)
              END IF

              nnq1 = nnq1 * nq2

            END DO

            RvecG(:) = reshape(RG1, shape=(/ nq /) )
            CALL dealloc_NParray(RG1,"RG1",name_sub)

          CASE (1) ! Sparse basis (Smolyak 1st implementation)
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' it does not work with SparseGrid_type=1'
            STOP 'SparseGrid_type=1'

            iq0_SG = 1
            iq1_SG = 0
            DO i_SG=1,basis_set%nb_SG
               nq_SG = get_nq_FROM_basis(basis_set%tab_PbasisSG(i_SG)%Pbasis)
               !write(out_unitp,*) 'i_SG,nq_SG',i_SG,nq_SG
!               iq1_SG = iq1_SG + nq_SG
               iq0_SG = 1
               DO ii_SG=1,i_SG-1
                 nq0_SG = get_nq_FROM_basis(basis_set%tab_PbasisSG(ii_SG)%Pbasis)
                 iq0_SG = iq0_SG + nq0_SG
               END DO
               iq1_SG = iq0_SG-1 + nq_SG
               !write(out_unitp,*) 'iq0_SG,iq1_SG',iq0_SG,iq1_SG
               CALL DerivOp_TO_RVecG(RVecG(iq0_SG:iq1_SG),nq_SG,        &
                                    basis_set%tab_PbasisSG(i_SG)%Pbasis,&
                                                            tab_der_loc)
!               iq0_SG = iq0_SG + nq_SG
            END DO

          CASE (2) ! Sparse basis (Smolyak 2d implementation)
            dnba_ind(:) = basis_set%Tabder_Qdyn_TO_Qbasis(tab_der_loc(:))

            IF (dnba_ind(1) == 0 .AND. dnba_ind(2) == 0) THEN
              CONTINUE ! nothing to do RvecG is unchanged
            ELSE
              CALL DerivOp_TO_RVecG_SGType2(RVecG,nq,                   &
                    basis_set%para_SGType2%nDind_SmolyakRep%Tab_DInd, &
                                     basis_set%nDindB%Tab_DInd,         &
                                     basis_set%tab_basisPrimSG,         &
                                     D=basis_set%nb_basis,              &
                                     LG=basis_set%L_SparseGrid,         &
                                     LB=basis_set%L_SparseBasis,        &
                                     tab_der=tab_der_loc)
            END IF

          CASE (4) ! Sparse basis (Smolyak 4th implementation)
            dnba_ind(:) = basis_set%Tabder_Qdyn_TO_Qbasis(tab_der_loc(:))

            IF (dnba_ind(1) == 0 .AND. dnba_ind(2) == 0) THEN
              CONTINUE ! nothing to do RvecG is unchanged
            ELSE
              !!  RVecG TO SRep
              CALL alloc2_SmolyakRep(SRep,basis_set%para_SGType2%nDind_SmolyakRep,&
                                   basis_set%tab_basisPrimSG,grid=.TRUE.,         &
                                   nb0=basis_set%para_SGType2%nb0)
              CALL tabR2bis_TO_SmolyakRep1(SRep,RVecG)

              IF (allocated(basis_set%para_SGType2%nDind_SmolyakRep%Tab_nDval)) THEN
                CALL DerivOp_TO3_GSmolyakRep(SRep,basis_set%para_SGType2,&
                           basis_set%tab_basisPrimSG,tab_der=tab_der_loc)
              ELSE
                CALL DerivOp_TO3_GSmolyakRep(SRep,basis_set%para_SGType2,&
                      basis_set%tab_basisPrimSG,tab_der=tab_der_loc)
              END IF

              CALL SmolyakRep2_TO_tabR1bis(RVecG,SRep)

              CALL dealloc_SmolyakRep(SRep)

            END IF

          CASE DEFAULT
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' WRONG SparseGrid_type',basis_set%SparseGrid_type
            write(out_unitp,*) ' The possibilities are: 0, 1, 2, 4'
            STOP
          END SELECT

        END IF

       IF (debug) THEN
         write(out_unitp,*) 'RvecG(:)',RvecG(:)
         write(out_unitp,*) 'END ',name_sub
       END IF

      END SUBROUTINE DerivOp_TO_RVecG


      RECURSIVE SUBROUTINE DerivOp_TO_CVecG(CVecG,nq,basis_set,tab_der)
        USE mod_basis_BtoG_GtoB_SGType4
        IMPLICIT NONE
        TYPE (basis), intent(in)         :: basis_set
        integer, intent(in)              :: nq
        complex (kind=Rkind), intent(inout) :: CVecG(:)
        integer, optional                :: tab_der(2)

        integer :: tab_der_loc(2),dnba_ind(2),iQact_ba,k
        complex (kind=Rkind), allocatable       :: CG1(:,:,:)
        complex (kind=Rkind), allocatable       :: CG2(:,:,:)
        real (kind=Rkind), allocatable          :: BGG(:,:)

        real (kind=Rkind), allocatable :: RVecG(:)
        !TYPE(Type_SmolyakRep)               :: SRep ! smolyak rep for SparseGrid_type=4
        TYPE(Type_SmolyakRepC)               :: SRep ! smolyak rep for SparseGrid_type=4


        integer                          :: nnq1,nnq3,nq2,iq1,iq2,iq3
        logical                          :: skip

        integer                          :: ibasis
        integer                          :: i,i_SG,iq0_SG,iq1_SG,nq_SG,ii_SG,nq0_SG
        integer                          :: nb_thread
        real (kind=Rkind) :: a

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='DerivOp_TO_CVecG'
      logical, parameter :: debug=.FALSE.
!      logical, parameter :: debug=.TRUE.
!---------------------------------------------------------------------

        IF (debug) THEN
          write(out_unitp,*) 'BEGINNING ',name_sub
          write(out_unitp,*) 'CvecG(:)',CvecG(:)
          !CALL RecWrite_basis(basis_set)
          CALL flush_perso(out_unitp)
        END IF

        IF (present(tab_der)) THEN
          tab_der_loc(:) = tab_der(:)
        ELSE
          tab_der_loc(:) = 0
        END IF
        WHERE (tab_der_loc < 0) tab_der_loc = 0

        IF (basis_set%cplx) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' the basis is complex but the'
          write(out_unitp,*) ' the vector is real !!'
          STOP
        END IF

        IF (BasisTOGrid_omp == 0) THEN
          nb_thread = 1
        ELSE
          nb_thread = BasisTOGrid_maxth
        END IF
        !write(out_unitp,*) 'nb_thread in ',name_sub,' : ',nb_thread

        IF (basis_set%packed_done) THEN ! Packed basis

          dnba_ind(:) = basis_set%Tabder_Qdyn_TO_Qbasis(tab_der_loc(:))

          IF (dnba_ind(1) == 0 .AND. dnba_ind(2) == 0) THEN ! dnba_ind(:)=0 => no derivative
            CONTINUE ! notthing to do
          ELSE
            CALL alloc_NParray(BGG,(/ nq,nq /),"BGG",name_sub)
            CALL Get_MatdnRGG(basis_set,BGG,dnba_ind)
            CVecG(:) = matmul(BGG,CVecG)
            CALL dealloc_NParray(BGG,"BGG",name_sub)

          END IF

        ELSE ! basis_set%nb_basis MUST BE > 0
          IF (basis_set%nb_basis == 0 ) STOP ' ERROR with packed_done !!!'

          SELECT CASE (basis_set%SparseGrid_type)
          CASE (0) ! Direct product

            nnq3 = nq
            nnq1 = 1
            nq2  = 1

            CALL alloc_NParray(CG1,(/ nnq1,nq2,nnq3 /),"CG1",name_sub)
            CG1(:,:,:) = reshape(CvecG,shape=(/ nnq1,nq2,nnq3 /))

            DO ibasis=basis_set%nb_basis,1,-1

              nq2  = get_nq_FROM_basis(basis_set%tab_Pbasis(ibasis)%Pbasis)
              nnq3 = nnq3 / nq2


              dnba_ind(:) = basis_set%tab_Pbasis(ibasis)%Pbasis%Tabder_Qdyn_TO_Qbasis(tab_der_loc(:))

              IF (dnba_ind(1) /= 0 .OR. dnba_ind(2) /= 0) THEN
                CALL alloc_NParray(CG2,(/ nnq1,nq2,nnq3 /),"CG2",name_sub)
                CG2(:,:,:) = reshape(CG1,shape=(/ nnq1,nq2,nnq3 /))

                IF (basis_set%tab_Pbasis(ibasis)%Pbasis%packed) THEN

                  CALL alloc_NParray(BGG,(/ nq2,nq2 /),"BGG",name_sub)
                  CALL Get_MatdnRGG(basis_set%tab_Pbasis(ibasis)%Pbasis,BGG,dnba_ind)

                 !$OMP parallel do default(none)                        &
                 !$OMP shared(BGG,CG2,nnq3,nnq1)                        &
                 !$OMP private(iq1,iq3)                                 &
                 !$OMP num_threads(nb_thread)
                  DO iq3=1,nnq3
                  DO iq1=1,nnq1
                     CG2(iq1,:,iq3) = matmul(BGG,CG2(iq1,:,iq3))
                  END DO
                  END DO
                 !$OMP end parallel do

                  CALL dealloc_NParray(BGG,"BGG",name_sub)

                ELSE
                  DO iq3=1,nnq3
                  DO iq1=1,nnq1
                    CALL DerivOp_TO_CVecG(CG2(iq1,:,iq3),nq2,           &
                                   basis_set%tab_Pbasis(ibasis)%Pbasis, &
                                                            tab_der_loc)
                  END DO
                  END DO
                END IF

                CALL dealloc_NParray(CG1,"CG1",name_sub)
                CALL alloc_NParray(CG1,(/ nnq1,nq2,nnq3 /),"CG1",name_sub)
                CG1 = CG2
                CALL dealloc_NParray(CG2,"CG2",name_sub)
              END IF

              nnq1 = nnq1 * nq2

            END DO

            CvecG(:) = reshape(CG1, shape=(/ nq /) )
            CALL dealloc_NParray(CG1,"CG1",name_sub)

          CASE (1) ! Sparse basis (Smolyak 1st implementation)
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' it does not work with SparseGrid_type=1'
            STOP 'SparseGrid_type=1'

            iq0_SG = 1
            iq1_SG = 0
            DO i_SG=1,basis_set%nb_SG
               nq_SG = get_nq_FROM_basis(basis_set%tab_PbasisSG(i_SG)%Pbasis)
               !write(out_unitp,*) 'i_SG,nq_SG',i_SG,nq_SG
!               iq1_SG = iq1_SG + nq_SG
               iq0_SG = 1
               DO ii_SG=1,i_SG-1
                 nq0_SG = get_nq_FROM_basis(basis_set%tab_PbasisSG(ii_SG)%Pbasis)
                 iq0_SG = iq0_SG + nq0_SG
               END DO
               iq1_SG = iq0_SG-1 + nq_SG
               !write(out_unitp,*) 'iq0_SG,iq1_SG',iq0_SG,iq1_SG
               CALL DerivOp_TO_CVecG(CVecG(iq0_SG:iq1_SG),nq_SG,        &
                                    basis_set%tab_PbasisSG(i_SG)%Pbasis,&
                                                            tab_der_loc)
!               iq0_SG = iq0_SG + nq_SG
            END DO

          CASE (2) ! Sparse basis (Smolyak 2d and 4th implementation)

            dnba_ind(:) = basis_set%Tabder_Qdyn_TO_Qbasis(tab_der_loc(:))

            IF (dnba_ind(1) == 0 .AND. dnba_ind(2) == 0) THEN
              CONTINUE ! nothing to do CvecG is unchanged
            ELSE

              CALL alloc_NParray(RVecG,shape(CVecG),'RVecG',name_sub)

              RVecG(:) = real(CVecG,kind=Rkind)
              CALL DerivOp_TO_RVecG_SGType2(RVecG,nq,                   &
                    basis_set%para_SGType2%nDind_SmolyakRep%Tab_DInd, &
                                     basis_set%nDindB%Tab_DInd,         &
                                     basis_set%tab_basisPrimSG,         &
                                     D=basis_set%nb_basis,              &
                                     LG=basis_set%L_SparseGrid,         &
                                     LB=basis_set%L_SparseBasis,        &
                                     tab_der=tab_der_loc)

              DO i=1,size(CVecG)
                a = RVecG(i)
                RVecG(i) = aimag(CVecG(i))
                CVecG(i) = cmplx(a,kind=Rkind)
              END DO
              CALL DerivOp_TO_RVecG_SGType2(RVecG,nq,                   &
                    basis_set%para_SGType2%nDind_SmolyakRep%Tab_DInd, &
                                     basis_set%nDindB%Tab_DInd,         &
                                     basis_set%tab_basisPrimSG,         &
                                     D=basis_set%nb_basis,              &
                                     LG=basis_set%L_SparseGrid,         &
                                     LB=basis_set%L_SparseBasis,        &
                                     tab_der=tab_der_loc)
              CVecG(:) = CVecG(:) + EYE*cmplx(RVecG(:),kind=Rkind)

              CALL dealloc_NParray(RVecG,'RVecG',name_sub)

            END IF

          CASE (4) ! Sparse basis (Smolyak 4th implementation)

            dnba_ind(:) = basis_set%Tabder_Qdyn_TO_Qbasis(tab_der_loc(:))

            IF (all(dnba_ind == 0)) THEN
              CONTINUE ! nothing to do RvecG is unchanged
            ELSE
              !!  RVecG TO SRep
              CALL alloc2_SmolyakRepC(SRep,basis_set%para_SGType2%nDind_SmolyakRep,&
                                     basis_set%tab_basisPrimSG,grid=.TRUE.,       &
                                     nb0=basis_set%para_SGType2%nb0)
              CALL tabC2bis_TO_SmolyakRepC1(SRep,CVecG)

              IF (allocated(basis_set%para_SGType2%nDind_SmolyakRep%Tab_nDval)) THEN
                CALL DerivOp_TO3_GSmolyakRepC(SRep,basis_set%para_SGType2,&
                      basis_set%tab_basisPrimSG,tab_der=tab_der_loc)
              ELSE
                CALL DerivOp_TO3_GSmolyakRepC(SRep,basis_set%para_SGType2,&
                      basis_set%tab_basisPrimSG,tab_der=tab_der_loc)
              END IF

              CALL SmolyakRepC2_TO_tabC1bis(CVecG,SRep)

              CALL dealloc_SmolyakRepC(SRep)

            END IF

          CASE DEFAULT
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' WRONG SparseGrid_type',basis_set%SparseGrid_type
            write(out_unitp,*) ' The possibilities are: 0, 1, 2, 4'
            STOP
          END SELECT

        END IF

       IF (debug) THEN
         write(out_unitp,*) 'CvecG(:)',CvecG(:)
         write(out_unitp,*) 'END ',name_sub
       END IF

      END SUBROUTINE DerivOp_TO_CVecG

      END MODULE mod_basis_BtoG_GtoB
