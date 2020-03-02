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
 MODULE mod_psi_Op
   IMPLICIT NONE

   CONTAINS


!=======================================================================================
!
!     Symmetrization (with abelian group) of psi in BasisRep
!
!=======================================================================================
      SUBROUTINE Set_symab_OF_psiBasisRep(psi,symab)
      USE mod_system
      USE mod_psi_set_alloc
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_psi)   :: psi

      integer, intent(in), optional :: symab

      integer :: loc_symab,ib
      integer :: Get_symabOFSymAbelianOFBasis_AT_ib ! function

!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='Set_symab_OF_psiBasisRep'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_ba',psi%nb_ba
        write(out_unitp,*) 'present(symab)? ',present(symab)
        IF (present(symab)) write(out_unitp,*) 'symab',symab
        write(out_unitp,*) 'psi BasisRep'
        !CALL ecri_psi(ZERO,psi)
      END IF

      IF (psi%BasisRep) THEN
        IF (psi%nb_bi == 1 .AND. psi%nb_be == 1) THEN
          IF (present(symab)) THEN
            loc_symab = symab
          ELSE
            ! find the symmtry (symab of the largest coef)
            IF (psi%cplx) THEN
              ib = maxloc(abs(psi%CvecB),dim=1)
            ELSE
              ib = maxloc(abs(psi%RvecB),dim=1)
            END IF
            loc_symab = Get_symabOFSymAbelianOFBasis_AT_ib(psi%BasisnD,ib)
            IF (debug) write(out_unitp,*) 'maxloc,symab',ib,loc_symab
          END IF
        ELSE
          loc_symab = -1
        END IF

        psi%symab = loc_symab
      ELSE
        psi%symab = -1
      END IF

      IF (debug) THEN
        write(out_unitp,*) 'symab, bits(symab)',WriteTOstring_symab(psi%symab)
      END IF

!-----------------------------------------------------------
      IF (psi%symab >= 0 .AND. psi%symab <= 7) THEN
        IF (psi%cplx .AND. allocated(psi%CvecB)) THEN
          DO ib=1,psi%nb_tot
            IF (psi%symab /= Get_symabOFSymAbelianOFBasis_AT_ib(psi%BasisnD,ib) ) &
                                                 psi%CvecB(ib) = CZERO
          END DO
        ELSE IF (.NOT. psi%cplx .AND. allocated(psi%RvecB)) THEN
          DO ib=1,psi%nb_tot
             IF (psi%symab /= Get_symabOFSymAbelianOFBasis_AT_ib(psi%BasisnD,ib) ) &
                                            psi%RvecB(ib) = ZERO
          END DO
        END IF
      END IF

!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'symab psi BasisRep',psi%symab
        !CALL ecri_psi(ZERO,psi)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------

      END SUBROUTINE Set_symab_OF_psiBasisRep
!=======================================================================================


!=======================================================================================
! MPI version, Symmetrization (with abelian group) of psi in BasisRep
!=======================================================================================
SUBROUTINE Set_symab_OF_psiBasisRep_MPI(psi,symab,changes)
  USE mod_system
  USE mod_psi_set_alloc
  USE mod_MPI
  IMPLICIT NONE

  TYPE(param_psi), intent(inout)      :: psi
  Integer,optional,intent(in)         :: symab
  Logical,optional,intent(inout)      :: changes

  Integer                             :: loc_symab
  Integer                             :: ib
  Integer                             :: Get_symabOFSymAbelianOFBasis_AT_ib ! function
  Logical                             :: broadcast

  broadcast=.FALSE.
  changes=.FALSE.
  
  IF(psi%BasisRep) THEN
    IF(psi%nb_bi==1 .AND. psi%nb_be==1) THEN
      IF(present(symab)) THEN
        loc_symab=symab
      ELSE
        broadcast=.TRUE.
        ! find the symmtry (symab of the largest coef)
        IF(MPI_id==0) THEN
          IF(psi%cplx) THEN
            ib=maxloc(abs(psi%CvecB),dim=1)
          ELSE
            ib=maxloc(abs(psi%RvecB),dim=1)
          ENDIF
          loc_symab=Get_symabOFSymAbelianOFBasis_AT_ib(psi%BasisnD,ib)
        ENDIF ! for MPI_id==0
      ENDIF
    ELSE
      loc_symab=-1
    ENDIF
    psi%symab=loc_symab
  ELSE
    psi%symab=-1
  ENDIF

  IF(broadcast) CALL MPI_Bcast(psi%symab,size1_MPI,MPI_Int_fortran,root_MPI,           &
                               MPI_COMM_WORLD,MPI_err)

  IF(psi%symab >= 0 .AND. psi%symab <= 7) THEN
    nb_per_MPI=psi%nb_tot/MPI_np
    nb_rem_MPI=mod(psi%nb_tot,MPI_np) 
    bound1_MPI=MPI_id*nb_per_MPI+1+MIN(MPI_id,nb_rem_MPI)
    bound2_MPI=(MPI_id+1)*nb_per_MPI+MIN(MPI_id,nb_rem_MPI)+merge(1,0,nb_rem_MPI>MPI_id)
    IF(MPI_id==0) THEN
      bound1_MPI=1
      bound2_MPI=psi%nb_tot
    ENDIF
    IF(psi%cplx .AND. allocated(psi%CvecB)) THEN      
      DO ib=bound1_MPI,bound2_MPI
        IF(psi%symab/=Get_symabOFSymAbelianOFBasis_AT_ib(psi%BasisnD,ib)) THEN
          psi%CvecB(ib)=CZERO
          changes=.TRUE.
        ENDIF
      ENDDO
    ELSE IF(.NOT. psi%cplx .AND. allocated(psi%RvecB)) THEN
      DO ib=bound1_MPI,bound2_MPI
        IF(psi%symab/=Get_symabOFSymAbelianOFBasis_AT_ib(psi%BasisnD,ib)) THEN
          psi%RvecB(ib)=ZERO
          changes=.TRUE.
        ENDIF
      ENDDO
    ENDIF
  ENDIF

END SUBROUTINE Set_symab_OF_psiBasisRep_MPI
!=======================================================================================

!================================================================
!
!     Overlap : <psi1 I psi2>
!
!================================================================
      !!@description: Overlap : $\langle \psi_1 | \psi_2\rangle
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      SUBROUTINE Overlap_psi1_psi2(Overlap,psi1,psi2,With_Grid,Channel_ie)
      USE mod_system
      USE mod_psi_set_alloc
      USE mod_MPI
      IMPLICIT NONE

!----- variables for the WP ----------------------------------------
      TYPE (param_psi), intent(in)    :: psi1,psi2
      complex (kind=Rkind)            :: Overlap
      logical, optional, intent(in)   :: With_Grid
      integer, optional, intent(in)   :: Channel_ie

!------ working variables ---------------------------------
      logical              :: With_Grid_loc
      integer              :: locChannel_ie
      integer              :: i_qa,i_qaie
      integer              :: i_be,i_bi,i_ba
      integer              :: i_baie,f_baie
      integer              :: i_modif_q
      real (kind=Rkind)    :: WrhonD
      complex (kind=Rkind) :: temp
      real (kind=Rkind)    :: Roverlap,Rtemp
      integer              :: iie,fie
      real (kind=Rkind), allocatable :: wrho(:)

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Overlap_psi1_psi2'
      logical,parameter :: debug = .FALSE.
!     logical,parameter :: debug = .TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'psi1'
        CALL ecri_psi(psi=psi1)

        write(out_unitp,*) 'psi2'
        CALL ecri_psi(psi=psi2)
        write(out_unitp,*) 'GridRep,BasisRep ?'
        IF (present(With_Grid)) write(out_unitp,*) 'With_Grid',With_Grid
        IF (present(Channel_ie)) write(out_unitp,*) 'Channel_ie',Channel_ie
      END IF
!-----------------------------------------------------------

      With_Grid_loc = .FALSE.

      IF (present(With_Grid)) With_Grid_loc = With_Grid

      locChannel_ie = 0
      IF (present(Channel_ie)) locChannel_ie = Channel_ie

      IF (psi1%nb_baie > psi1%nb_tot) THEN
        With_Grid_loc = .FALSE.
      END IF

      ! With_Grid_loc: F
      IF(MPI_id==0) THEN
        IF (With_Grid_loc) THEN
          IF (psi1%cplx .AND.                                             &
           allocated(psi1%CvecG) .AND. allocated(psi2%CvecG)) THEN
          ELSE IF (.NOT. psi1%cplx .AND.                                  &
           allocated(psi1%RvecG) .AND. allocated(psi2%RvecG)) THEN
          ELSE
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' impossible to calculate the GridRep overlap'
            write(out_unitp,*) ' With_Grid_loc=t but problem with the allocation GridRep'
            write(out_unitp,*) 'allocated(psi1%CvecG)',allocated(psi1%CvecG)
            write(out_unitp,*) 'allocated(psi2%CvecG)',allocated(psi2%CvecG)
            write(out_unitp,*) 'allocated(psi1%RvecG)',allocated(psi1%RvecG)
            write(out_unitp,*) 'allocated(psi2%RvecG)',allocated(psi2%RvecG)
            write(out_unitp,*) ' psi1'
            CALL ecri_psi(psi=psi1,ecri_GridRep=.TRUE.)
            write(out_unitp,*) ' psi2'
            CALL ecri_psi(psi=psi2,ecri_GridRep=.TRUE.)
            STOP
          END IF
        ELSE
          IF (psi1%cplx .AND.                                             &
           allocated(psi1%CvecB) .AND. allocated(psi2%CvecB)) THEN
          ELSE IF (.NOT. psi1%cplx .AND.                                  &
           allocated(psi1%RvecB) .AND. allocated(psi2%RvecB)) THEN
          ELSE
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' impossible to calculate the BasisRep overlap'
            write(out_unitp,*) ' With_Grid_loc=f (on basis) but problem with the allocation of BasisRep'
            write(out_unitp,*) 'allocated(psi1%CvecB)',allocated(psi1%CvecB)
            write(out_unitp,*) 'allocated(psi2%CvecB)',allocated(psi2%CvecB)
            write(out_unitp,*) 'allocated(psi1%RvecB)',allocated(psi1%RvecB)
            write(out_unitp,*) 'allocated(psi2%RvecB)',allocated(psi2%RvecB)
            write(out_unitp,*) ' psi1'
            CALL ecri_psi(psi=psi1,ecri_BasisRep=.TRUE.)
            write(out_unitp,*) ' psi2'
            CALL ecri_psi(psi=psi2,ecri_BasisRep=.TRUE.)
            STOP
          END IF
        END IF

        IF (.NOT. With_Grid_loc) THEN
          i_baie=1
          f_baie=psi1%nb_tot
          IF (psi1%nb_tot == psi1%nb_baie .AND.  locChannel_ie > 0 .AND.  &
                                  locChannel_ie <= psi1%ComOp%nb_bie) THEN
            i_baie = 1 + (locChannel_ie-1)*psi1%nb_ba
            f_baie = i_baie-1 + psi1%nb_ba
          END IF
          IF (psi1%symab > -1 .AND. psi2%symab > -1 .AND. psi1%symab /= psi2%symab) THEN
            Overlap = cmplx(ZERO,ZERO,kind=Rkind)
          ELSE
            IF (psi1%cplx) THEN
              Overlap = dot_product( psi1%CvecB(i_baie:f_baie) ,          &
                                     psi2%CvecB(i_baie:f_baie) )
            ELSE
              ROverlap = dot_product( psi1%RvecB(i_baie:f_baie) ,         &
                                      psi2%RvecB(i_baie:f_baie) )
              Overlap = cmplx(ROverlap,ZERO,kind=Rkind)
            END IF
          END IF

        ELSE

  !       - initialization ----------------------------------
          Overlap = cmplx(ZERO,ZERO,kind=Rkind)

          CALL alloc_NParray(wrho,(/ psi1%nb_qa/),"wrho",name_sub)
          DO i_qa=1,psi1%nb_qa
            wrho(i_qa) = Rec_WrhonD(psi1%BasisnD,i_qa)
          END DO

          IF (psi1%cplx) THEN
            iie = 1
            fie = psi1%nb_qa
            DO i_be=1,psi1%nb_be
            DO i_bi=1,psi1%nb_bi
              Overlap = Overlap + dot_product(                            &
                psi1%CvecG(iie:fie),wrho*psi2%CvecG(iie:fie))
              iie = iie + psi1%nb_qa
              fie = fie + psi1%nb_qa
            END DO
            END DO
          ELSE
            iie = 1
            fie = psi1%nb_qa
            DO i_be=1,psi1%nb_be
            DO i_bi=1,psi1%nb_bi
              Overlap = Overlap + cmplx(dot_product(                      &
                psi1%RvecG(iie:fie),wrho*psi2%RvecG(iie:fie)) ,kind=Rkind)
              iie = iie + psi1%nb_qa
              fie = fie + psi1%nb_qa
            END DO
            END DO
          END IF

          CALL dealloc_NParray(wrho,"wrho",name_sub)

        END IF
      ENDIF ! for MPI_id==0

!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Overlap : ',Overlap
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------

      END SUBROUTINE Overlap_psi1_psi2

!=======================================================================================
! subroutine for calculation of matrix H_overlap(i,j) for <psi(i)|Hpsi(j)> 
! MPI version, takes too much memory
!=======================================================================================
#if(run_MPI)
SUBROUTINE Overlap_psi_Hpsi_matrix_MPI(H_overlap,S_overlap,psi,Hpsi,ndim,With_Grid)
  USE mod_system
  USE mod_psi_set_alloc
  USE mod_MPI
  USE mod_MPI_Aid
  IMPLICIT NONE
  
  TYPE(param_psi), intent(inout)              :: psi(:)  !< inout only non-root threats
  TYPE(param_psi), intent(inout)              :: Hpsi(:) !< inout only non-root threats
  Logical,optional,intent(in)                 :: With_Grid
  Integer,         intent(in)                 :: ndim
  
  Real(kind=Rkind),allocatable,intent(inout)  :: H_overlap(:,:)
  Real(kind=Rkind),allocatable,intent(inout)  :: S_overlap(:,:)

  Character(len=*),parameter                  :: name_sub='Overlap_psi_Hpsi_matrix_MPI'
  !Real(kind=Rkind),allocatable                :: H_flat(:)
  !Real(kind=Rkind),allocatable                :: S_flat(:)
  Complex(kind=Rkind)                         :: Overlap
  Integer                                     :: num
  Integer                                     :: ii
  Integer                                     :: jj

  IF(allocated(H_overlap)) deallocate(H_overlap)
  IF(allocated(S_overlap)) deallocate(S_overlap)
  CALL alloc_NParray(H_overlap,(/ ndim,ndim /),"H",name_sub)
  CALL alloc_NParray(S_overlap,(/ ndim,ndim /),"S",name_sub)

  !allocate(H_flat(ndim*ndim))
  !allocate(S_flat(ndim*ndim))
  
  ! calculate on master without MPI
  IF(MPI_id==0) THEN
    num=0
    DO ii=1,ndim
      DO jj=1,ndim
        num=num+1
        CALL Overlap_psi1_psi2(Overlap,psi(jj),Hpsi(ii),With_Grid=With_Grid)
        H_overlap(jj,ii)=real(Overlap,kind=Rkind)
        
        CALL Overlap_psi1_psi2(Overlap,psi(jj), psi(ii),With_Grid=With_Grid)
        S_overlap(jj,ii)=real(Overlap,kind=Rkind)
        
        !H_flat(num)=H_overlap(jj,ii)
        !S_flat(num)=S_overlap(jj,ii)
      ENDDO
    ENDDO
  ENDIF ! for MPI_id==0
 
  CALL MPI_Bcast_matrix(H_overlap,1,ndim,1,ndim,root_MPI)
  CALL MPI_Bcast_matrix(S_overlap,1,ndim,1,ndim,root_MPI)
  
!  CALL MPI_Bcast(H_flat,ndim*ndim,MPI_Real8,root_MPI,MPI_COMM_WORLD,MPI_err)
!  CALL MPI_Bcast(S_flat,ndim*ndim,MPI_Real8,root_MPI,MPI_COMM_WORLD,MPI_err)
!
!  IF(MPI_id/=0) THEN
!    num=0
!    DO ii=1,ndim
!      DO jj=1,ndim
!        num=num+1
!        H_overlap(jj,ii)=H_flat(num)
!        S_overlap(jj,ii)=S_flat(num)
!      ENDDO
!    ENDDO
!  ENDIF ! for MPI_id/=0
  
  ! Using MPI, takes too much memory'
  ! psi(1:ndim) are share on all threads, Hpsi(:) are distributed to each threads
  !CALL Overlap_psi1_psi2_MPI(H_overlap,S_overlap,psi,Hpsi,ndim,With_Grid)
  
END SUBROUTINE Overlap_psi_Hpsi_matrix_MPI
#endif
!=======================================================================================


!=======================================================================================
! subroutine for calculation of matrix H_overlap(i,j) for <psi(i)|Hpsi(j)> 
! MPI version, less memory
!=======================================================================================
#if(run_MPI)
SUBROUTINE Overlap_psi_Hpsi_matrix_MPI2(H_overlap,S_overlap,psi,Hpsi,ndim,With_Grid)
  USE mod_system
  USE mod_psi_set_alloc
  USE mod_MPI
  USE mod_MPI_Aid
  IMPLICIT NONE
  
  TYPE(param_psi), intent(inout)              :: psi(:)  !< inout only non-root threats
  TYPE(param_psi), intent(inout)              :: Hpsi(:) !< inout only non-root threats
  Logical,optional,intent(in)                 :: With_Grid
  Integer,         intent(in)                 :: ndim
  
  Real(kind=Rkind),allocatable,intent(inout)  :: H_overlap(:,:)
  Real(kind=Rkind),allocatable,intent(inout)  :: S_overlap(:,:)

  Character(len=*),parameter                  :: name_sub='Overlap_psi_Hpsi_matrix_MPI'
  Complex(kind=Rkind)                         :: Overlap
  Integer                                     :: ii
  Integer                                     :: jj

  IF(allocated(H_overlap)) deallocate(H_overlap)
  IF(allocated(S_overlap)) deallocate(S_overlap)
  CALL alloc_NParray(H_overlap,(/ ndim,ndim /),"H",name_sub)
  CALL alloc_NParray(S_overlap,(/ ndim,ndim /),"S",name_sub)

  !-------------------------------------------------------------------------------------
  ! calculate on master without MPI
!  IF(MPI_id==0) THEN
!    DO ii=1,ndim
!      DO jj=1,ndim
!        CALL Overlap_psi1_psi2(Overlap,psi(jj),Hpsi(ii),With_Grid=With_Grid)
!        H_overlap(jj,ii)=real(Overlap,kind=Rkind)
!
!        CALL Overlap_psi1_psi2(Overlap,psi(jj), psi(ii),With_Grid=With_Grid)
!        S_overlap(jj,ii)=real(Overlap,kind=Rkind)
!
!      ENDDO
!    ENDDO
!  ENDIF ! for MPI_id==0
!  
!  CALL MPI_Bcast_matrix(H_overlap,1,ndim,1,ndim,root_MPI)
!  CALL MPI_Bcast_matrix(S_overlap,1,ndim,1,ndim,root_MPI)

  !-------------------------------------------------------------------------------------
  ! calculate with MPI
  CALL Overlap_psi1_psi2_MPI2(H_overlap,S_overlap,psi,Hpsi,ndim,With_Grid)

END SUBROUTINE Overlap_psi_Hpsi_matrix_MPI2
#endif
!=======================================================================================

!=======================================================================================
!> Subroutine for the calculation of matrix H_overlap(i,j) for <psi(i)|Hpsi(j)> 
!>                                      and S_overlap(i,j) for <psi(i)| psi(j)> 
!> MPI V3

!> overlap is cut into submatrix for MPI parallel   
!> note, this requre new  Overlap_psipsi_MPI3 routine
!
!> the submatrix of psi and Hpsi are keeped also for the calculation 
!> of Residual g in MakeResidual_Davidson_MPI2
!=======================================================================================
#if(run_MPI)
SUBROUTINE Overlap_HS_matrix_MPI3(H_overlap,S_overlap,psi,Hpsi,ndim0,ndim,With_Grid)
  USE mod_system
  USE mod_psi_set_alloc
  USE mod_MPI
  USE mod_MPI_Aid
  IMPLICIT NONE
  
  TYPE(param_psi), intent(inout)              :: psi(:)  !< inout only non-root threats
  TYPE(param_psi), intent(inout)              :: Hpsi(:) !< inout only non-root threats
  Logical,optional,intent(in)                 :: With_Grid
  Integer,         intent(in)                 :: ndim0
  Integer,         intent(in)                 :: ndim
  
  Real(kind=Rkind),allocatable,intent(inout)  :: H_overlap(:,:)
  Real(kind=Rkind),allocatable,intent(inout)  :: S_overlap(:,:)
  
  Real(kind=Rkind),allocatable                :: H0_overlap(:,:)
  Real(kind=Rkind),allocatable                :: S0_overlap(:,:)
  
  Real(kind=Rkind),allocatable                :: H_overlapp(:,:)
  Real(kind=Rkind),allocatable                :: S_overlapp(:,:)

  Character(len=*),parameter                  :: name_sub='Overlap_psi_Hpsi_matrix_MPI'
  Complex(kind=Rkind)                         :: Overlap
  Integer                                     :: ii
  Integer                                     :: jj

  
  CALL increase_martix(H_overlap,name_sub,ndim0,ndim)
  CALL increase_martix(S_overlap,name_sub,ndim0,ndim)
  
!  IF(allocated(H_overlap)) THEN
!!    CALL alloc_NParray(H0_overlap,(/ ndim,ndim /),"H0",name_sub)
!!    H0_overlap(1:ndim0,1:ndim0)=H_overlap(1:ndim0,1:ndim0)
!!    deallocate(H_overlap)
!!    CALL alloc_NParray(H_overlap,(/ ndim,ndim /),"H0",name_sub)
!!    H_overlap(1:ndim0,1:ndim0)=H0_overlap(1:ndim0,1:ndim0)
!    CALL alloc_NParray(H0_overlap,(/ ndim,ndim /),"H0",name_sub)
!    H0_overlap(1:ndim0,1:ndim0)=H_overlap(1:ndim0,1:ndim0)
!    CALL move_alloc(H0_overlap,H_overlap) ! moves the allocation from H0_overlap to ...
!  ELSE
!    CALL alloc_NParray(H_overlap,(/ ndim,ndim /),"H",name_sub)
!  ENDIF
!  
!  IF(allocated(S_overlap)) THEN
!!    CALL alloc_NParray(S0_overlap,(/ ndim,ndim /),"S0",name_sub)
!!    S0_overlap(1:ndim0,1:ndim0)=S_overlap(1:ndim0,1:ndim0)
!!    deallocate(S_overlap)
!!    CALL alloc_NParray(S_overlap,(/ ndim,ndim /),"S",name_sub)
!!    S_overlap(1:ndim0,1:ndim0)=S0_overlap(1:ndim0,1:ndim0)
!    CALL alloc_NParray(S0_overlap,(/ ndim,ndim /),"S0",name_sub)
!    S0_overlap(1:ndim0,1:ndim0)=S_overlap(1:ndim0,1:ndim0)
!    CALL move_alloc(S0_overlap,S_overlap) ! moves the allocation from S0_overlap to ...
!  ELSE
!    CALL alloc_NParray(S_overlap,(/ ndim,ndim /),"H",name_sub)
!  ENDIF

  !-------------------------------------------------------------------------------------
  ! calculate on master without MPI
  CALL alloc_NParray(H_overlapp,(/ ndim,ndim /),"H",name_sub)
  CALL alloc_NParray(S_overlapp,(/ ndim,ndim /),"S",name_sub)
  IF(MPI_id==0) THEN
    DO ii=1,ndim
      DO jj=1,ndim
        CALL Overlap_psi1_psi2(Overlap,psi(jj),Hpsi(ii),With_Grid=With_Grid)
        H_overlapp(jj,ii)=real(Overlap,kind=Rkind)

        CALL Overlap_psi1_psi2(Overlap,psi(jj), psi(ii),With_Grid=With_Grid)
        S_overlapp(jj,ii)=real(Overlap,kind=Rkind)

      ENDDO
    ENDDO
  ENDIF ! for MPI_id==0
  
  CALL MPI_Bcast_matrix(H_overlapp,1,ndim,1,ndim,root_MPI)
  CALL MPI_Bcast_matrix(S_overlapp,1,ndim,1,ndim,root_MPI)

  !-------------------------------------------------------------------------------------
  ! calculate with MPI
  CALL Overlap_psi1_psi2_MPI5(H_overlap,S_overlap,psi,Hpsi,ndim0,ndim,With_Grid)
  
!  write(*,*) 'overlapp check',MAXVAL(ABS(H_overlapp-H_overlap)),                       &
!                              MAXVAL(ABS(S_overlapp-S_overlap)),                       &
!                              MAXVAL(ABS(H_overlap)),MAXVAL(ABS(S_overlap)),           &
!                              MAXVAL(ABS(H_overlapp)),MAXVAL(ABS(S_overlapp))
!  Stop

END SUBROUTINE Overlap_HS_matrix_MPI3
#endif
!=======================================================================================

!=======================================================================================
!> MPI Subroutine for the calculation of matrix H_overlap(i,j) for <psi(i)|Hpsi(j)> 

!> overlap is cut into submatrix for MPI parallel   
!> note, this requre new Overlap_psipsi_MPI3 routine
!
!> the submatrix of psi have been distribed in sub_NewVec_Davidson
!> Hpsi are distribed and keeped also for the calculation
!> of Residual g in MakeResidual_Davidson_MPI3
!=======================================================================================
#if(run_MPI)
SUBROUTINE Overlap_H_matrix_MPI4(H_overlap,psi,Hpsi,ndim0,ndim,With_Grid)
  USE mod_system
  USE mod_psi_set_alloc
  USE mod_MPI
  USE mod_MPI_Aid
  IMPLICIT NONE
  
  TYPE(param_psi), intent(inout)              :: psi(:)  !< inout only non-root threats
  TYPE(param_psi), intent(inout)              :: Hpsi(:) !< inout only non-root threats
  Logical,optional,intent(in)                 :: With_Grid
  Integer,         intent(in)                 :: ndim0
  Integer,         intent(in)                 :: ndim
  
  Real(kind=Rkind),allocatable,intent(inout)  :: H_overlap(:,:)
  Real(kind=Rkind),allocatable                :: H0_overlap(:,:)
!  Real(kind=Rkind),allocatable                :: H_overlapp(:,:)

  Character(len=*),parameter                  :: name_sub='Overlap_H_matrix_MPI4'
  Complex(kind=Rkind)                         :: Overlap
  Integer                                     :: ii
  Integer                                     :: jj


  CALL increase_martix(H_overlap,name_sub,ndim0,ndim)
  
!  IF(allocated(H_overlap)) THEN
!!    CALL alloc_NParray(H0_overlap,(/ ndim,ndim /),"H0",name_sub)
!!    H0_overlap(1:ndim0,1:ndim0)=H_overlap(1:ndim0,1:ndim0)
!!    deallocate(H_overlap)
!!    CALL alloc_NParray(H_overlap,(/ ndim,ndim /),"H0",name_sub)
!!    H_overlap(1:ndim0,1:ndim0)=H0_overlap(1:ndim0,1:ndim0)
!    CALL alloc_NParray(H0_overlap,(/ ndim,ndim /),"H0",name_sub)
!    H0_overlap(1:ndim0,1:ndim0)=H_overlap(1:ndim0,1:ndim0)
!    CALL move_alloc(H0_overlap,H_overlap) ! moves the allocation from H0_overlap to ...
!  ELSE
!    CALL alloc_NParray(H_overlap,(/ ndim,ndim /),"H",name_sub)
!  ENDIF

  !-------------------------------------------------------------------------------------
  ! calculate on master without MPI
!  CALL alloc_NParray(H_overlapp,(/ ndim,ndim /),"H",name_sub)
!  IF(MPI_id==0) THEN
!    DO ii=1,ndim
!      DO jj=1,ndim
!        CALL Overlap_psi1_psi2(Overlap,psi(jj),Hpsi(ii),With_Grid=With_Grid)
!        H_overlapp(jj,ii)=real(Overlap,kind=Rkind)
!      ENDDO
!    ENDDO
!  ENDIF ! for MPI_id==0
!  
!  CALL MPI_Bcast_matrix(H_overlapp,1,ndim,1,ndim,root_MPI)

  !-------------------------------------------------------------------------------------
  ! calculate with MPI
  CALL Overlap_psi_Hpsi_MPI(H_overlap,psi,Hpsi,ndim0,ndim,With_Grid)
  
!  write(*,*) 'H_overlapp check',MAXVAL(ABS(H_overlapp-H_overlap)),    &
!                                MAXVAL(ABS(H_overlap)),MAXVAL(ABS(H_overlapp))

END SUBROUTINE Overlap_H_matrix_MPI4
#endif
!=======================================================================================

!=======================================================================================
!> MPI Subroutine for the calculation of matrix S_overlap(i,j) for <psi(i)|psi(j)> 

!> overlap is cut into submatrix for MPI parallel   
!> note, this requre new Overlap_psipsi_MPI3 routine
!
!> psi are distribed and keeped also for the calculation
!> of H_overlap in sub_MakeHPsi_Davidson 
!> and Residual g in MakeResidual_Davidson_MPI3
!=======================================================================================
#if(run_MPI)
SUBROUTINE Overlap_S_matrix_MPI4(S_overlap,psi,ndim0,ndim,With_Grid)
  USE mod_system
  USE mod_psi_set_alloc
  USE mod_MPI
  USE mod_MPI_Aid
  IMPLICIT NONE
  
  TYPE(param_psi), intent(inout)              :: psi(:)  !< inout only non-root threats
  Logical,optional,intent(in)                 :: With_Grid
  Integer,         intent(in)                 :: ndim0
  Integer,         intent(in)                 :: ndim
  
  Real(kind=Rkind),allocatable,intent(inout)  :: S_overlap(:,:)
  Real(kind=Rkind),allocatable                :: S0_overlap(:,:)
  !Real(kind=Rkind),allocatable                :: S_overlapp(:,:)

  Character(len=*),parameter                  :: name_sub='Overlap_S_matrix_MPI4'
  Complex(kind=Rkind)                         :: Overlap
  Integer                                     :: ii
  Integer                                     :: jj


  CALL increase_martix(S_overlap,name_sub,ndim0,ndim)
  
!  IF(allocated(S_overlap)) THEN
!!    CALL alloc_NParray(S0_overlap,(/ ndim,ndim /),"S0",name_sub)
!!    S0_overlap(1:ndim0,1:ndim0)=S_overlap(1:ndim0,1:ndim0)
!!    deallocate(S_overlap)
!!    CALL alloc_NParray(S_overlap,(/ ndim,ndim /),"S",name_sub)
!!    S_overlap(1:ndim0,1:ndim0)=S0_overlap(1:ndim0,1:ndim0)
!    CALL alloc_NParray(S0_overlap,(/ ndim,ndim /),"S0",name_sub)
!    S0_overlap(1:ndim0,1:ndim0)=S_overlap(1:ndim0,1:ndim0)
!    CALL move_alloc(S0_overlap,S_overlap) ! moves the allocation from S0_overlap to ...
!  ELSE
!    CALL alloc_NParray(S_overlap,(/ ndim,ndim /),"H",name_sub)
!  ENDIF
  
  !-------------------------------------------------------------------------------------
  ! calculate on master without MPI
!  CALL alloc_NParray(S_overlapp,(/ ndim,ndim /),"S",name_sub)
!  IF(MPI_id==0) THEN
!    DO ii=1,ndim
!      DO jj=1,ndim
!        CALL Overlap_psi1_psi2(Overlap,psi(jj), psi(ii),With_Grid=With_Grid)
!        S_overlapp(jj,ii)=real(Overlap,kind=Rkind)
!      ENDDO
!    ENDDO
!  ENDIF ! for MPI_id==0
!  
!  CALL MPI_Bcast_matrix(S_overlapp,1,ndim,1,ndim,root_MPI)

  !-------------------------------------------------------------------------------------
  ! calculate with MPI
  !CALL Overlap_psi_psi_MPI(S_overlap,psi,ndim0,ndim,With_Grid)
  ! no distribution  of psi
  CALL calculate_overlap_MPI(psi,ndim0+1,ndim,With_Grid=With_Grid,S_overlap=S_overlap)


!  write(*,*) 'S_overlapp check',MAXVAL(ABS(S_overlapp-S_overlap)),    &
!                                MAXVAL(ABS(S_overlap)),MAXVAL(ABS(S_overlapp))

END SUBROUTINE Overlap_S_matrix_MPI4
#endif
!=======================================================================================

!=======================================================================================
!> Overlap_S_matrix_MPI without the distribution of psi
!=======================================================================================
#if(run_MPI)
SUBROUTINE Overlap_S0_matrix_MPI4(S_overlap,psi,ndim0,ndim,With_Grid)
  USE mod_system
  USE mod_psi_set_alloc
  USE mod_MPI
  USE mod_MPI_Aid
  IMPLICIT NONE
  
  TYPE(param_psi), intent(inout)              :: psi(:)  !< inout only non-root threats
  Real(kind=Rkind),allocatable,intent(inout)  :: S_overlap(:,:)
  Logical,optional,intent(in)                 :: With_Grid
  Integer,         intent(in)                 :: ndim0
  Integer,         intent(in)                 :: ndim
  
  Character(len=*),parameter                  :: name_sub='Overlap_S_matrix_MPI4'
  Complex(kind=Rkind)                         :: Overlap
  Integer                                     :: ii
  Integer                                     :: jj


  CALL increase_martix(S_overlap,name_sub,ndim0,ndim)
  
  ! calculate without the distribution of psi
  CALL calculate_overlap_MPI(psi,ndim0+1,ndim,With_Grid=With_Grid,S_overlap=S_overlap)

END SUBROUTINE Overlap_S0_matrix_MPI4
#endif
!=======================================================================================

!=======================================================================================
!> subroutine for the calculation of Overlap_psi_Hpsi with MPI 
!> overlap is cut into submatrix for MPI parallel   
!> note, this requre new Overlap_psi1_psi2 routine
!
!> the submatrix of Hpsi are keeped also for the calculation 
!> of Residual g in MakeResidual_Davidson_MPI2
!=======================================================================================
#if(run_MPI)
SUBROUTINE Overlap_psi_Hpsi_MPI(H_overlap,psi,Hpsi,ndim0,ndim,With_Grid)
  USE mod_system
  USE mod_psi_set_alloc
  USE mod_MPI
  USE mod_MPI_Aid
  IMPLICIT NONE
  
  TYPE(param_psi), intent(inout)          :: psi(:)  !< on non-root threads allocated
  TYPE(param_psi), intent(inout)          :: Hpsi(:) !< on non-root threads allocated
  Real(kind=Rkind),intent(inout)          :: H_overlap(:,:)
  Logical,optional,intent(in)             :: With_Grid
  Integer,         intent(in)             :: ndim0
  Integer,         intent(in)             :: ndim

  complex(kind=Rkind)                     :: Overlap
  logical                                 :: With_Grid_loc
  Integer                                 :: i
  Integer                                 :: j
  Logical                                 :: root_jobs
  
  With_Grid_loc=.FALSE.
  IF (present(With_Grid)) With_Grid_loc = With_Grid

  ! only on master: psi%CvecG,psi%RvecG,psi%CvecB,psi%RvecB
  CALL distribute_psi_pack_MPI(Hpsi,ndim0+1,ndim,With_Grid=With_Grid_loc)
  CALL calculate_overlap_MPI(psi,ndim0+1,ndim,With_Grid=With_Grid_loc,Hpsi=Hpsi,       &
                             H_overlap=H_overlap)
                             
END SUBROUTINE Overlap_psi_Hpsi_MPI
#endif
!=======================================================================================

!=======================================================================================
!> subroutine for the calculation of S_overlap(Overlap_psi_psi) with MPI 
!> overlap is cut into submatrix for MPI parallel   
!> note, this requre new Overlap_psi1_psi2 routine
!
!> the submatrix of psi are keeped also for the calculation 
!> of H_overlap and Residual g in MakeResidual_Davidson_MPI3
!=======================================================================================
#if(run_MPI)
SUBROUTINE Overlap_psi_psi_MPI(S_overlap,psi,ndim0,ndim,With_Grid)
  USE mod_system
  USE mod_psi_set_alloc
  USE mod_MPI
  USE mod_MPI_Aid
  IMPLICIT NONE
  
  TYPE(param_psi), intent(inout)          :: psi(:)  !< on non-root threads allocated
  Real(kind=Rkind),intent(inout)          :: S_overlap(:,:)
  Logical,optional,intent(in)             :: With_Grid
  Integer,         intent(in)             :: ndim0
  Integer,         intent(in)             :: ndim

  logical                                 :: With_Grid_loc
  Integer                                 :: i,j
  Logical                                 :: root_jobs
  
  With_Grid_loc=.FALSE.
  IF(present(With_Grid)) With_Grid_loc=With_Grid

  ! only on master: psi%CvecG,psi%RvecG,psi%CvecB,psi%RvecB
  CALL distribute_psi_pack_MPI(psi,ndim0+1,ndim,With_Grid=With_Grid_loc)
  CALL calculate_overlap_MPI(psi,ndim0+1,ndim,With_Grid=With_Grid_loc,                 &
                             S_overlap=S_overlap)
                             
END SUBROUTINE Overlap_psi_psi_MPI
#endif
!=======================================================================================

!=======================================================================================
!> subroutine for distribute psi(ndim1:ndim2) to different threads 
!> the distribution is depends on the length of vec (RvecB, CvecB, RvecG, or CvecG)
!=======================================================================================
SUBROUTINE distribute_psi_MPI(psi,ndim1,ndim2,With_Grid)
  USE mod_system
  USE mod_psi_set_alloc
  USE mod_MPI
  USE mod_MPI_Aid
  IMPLICIT NONE
  
  TYPE(param_psi), intent(inout)          :: psi(:) 
  Integer,         intent(in)             :: ndim1
  Integer,         intent(in)             :: ndim2
  Logical,optional,intent(in)             :: With_Grid
 
  Integer                                 :: i,j
  Logical                                 :: With_Grid_loc

  With_Grid_loc=.FALSE.
  IF (present(With_Grid)) With_Grid_loc=With_Grid

  IF(MPI_id==0) THEN
    DO i=ndim1,ndim2
      IF(.NOT. With_Grid_loc) THEN
        nb_per_MPI=psi(i)%nb_tot/MPI_np
        nb_rem_MPI=mod(psi(i)%nb_tot,MPI_np) !remainder jobs
      ELSE
        nb_per_MPI=psi(i)%nb_qaie/MPI_np
        nb_rem_MPI=mod(psi(i)%nb_qaie,MPI_np) !remainder jobs
      ENDIF
      
      ! Send array
      DO i_MPI=1,MPI_np-1
        bound1_MPI=i_MPI*nb_per_MPI+1+MIN(i_MPI,nb_rem_MPI)
        bound2_MPI=(i_MPI+1)*nb_per_MPI+MIN(i_MPI,nb_rem_MPI)                          &
                                       +merge(1,0,nb_rem_MPI>i_MPI)
        IF (.NOT. With_Grid_loc) THEN
          IF (psi(i)%cplx) THEN
            CALL MPI_Send(psi(i)%CvecB(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1, &
                          MPI_Complex8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
          ELSE
            CALL MPI_Send(psi(i)%RvecB(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1, &
                          MPI_Real8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
          ENDIF
        ELSE
          IF (psi(i)%cplx) THEN
            CALL MPI_Send(psi(i)%CvecG(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1, &
                          MPI_Complex8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)        
          ELSE
            CALL MPI_Send(psi(i)%RvecG(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1, &
                          MPI_Real8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
          ENDIF
        ENDIF ! for .NOT. With_Grid_loc 
      ENDDO ! for i_MPI=1,MPI_np-1
    ENDDO ! for i=ndim1,ndim2
  ENDIF ! for MPI_id==0
  
  ! MPI/=0------------------------------------------------------------------------------
  IF(MPI_id/=0) THEN
    DO i=ndim1,ndim2
      IF(.NOT. With_Grid_loc) THEN
        nb_per_MPI=psi(i)%nb_tot/MPI_np
        nb_rem_MPI=mod(psi(i)%nb_tot,MPI_np) !remainder jobs
      ELSE
        nb_per_MPI=psi(i)%nb_qaie/MPI_np
        nb_rem_MPI=mod(psi(i)%nb_qaie,MPI_np) !remainder jobs
      ENDIF
      bound1_MPI=MPI_id*nb_per_MPI+1+MIN(MPI_id,nb_rem_MPI)
      bound2_MPI=(MPI_id+1)*nb_per_MPI+MIN(MPI_id,nb_rem_MPI)                          &
                                      +merge(1,0,nb_rem_MPI>MPI_id)
      
      ! receive array
      IF (.NOT. With_Grid_loc) THEN
        IF (psi(i)%cplx) THEN
          IF(.NOT. allocated(psi(i)%CvecB))                                            &
               allocate(psi(i)%CvecB(bound1_MPI:bound2_MPI))
          CALL MPI_Recv(psi(i)%CvecB(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,   &
                        MPI_Complex8,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
        ELSE
          IF(.NOT. allocated(psi(i)%RvecB))                                            &
               allocate(psi(i)%RvecB(bound1_MPI:bound2_MPI))
          CALL MPI_Recv( psi(i)%RvecB(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,  &
                        MPI_Real8,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
        ENDIF
      ELSE
        IF (psi(i)%cplx) THEN
          IF(.NOT. allocated(psi(i)%CvecG))                                            &
               allocate(psi(i)%CvecG(bound1_MPI:bound2_MPI))
          CALL MPI_Recv(psi(i)%CvecG(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,   &
                        MPI_Complex8,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
        ELSE
          IF(.NOT. allocated(psi(i)%RvecG))                                            &
               allocate(psi(i)%RvecG(bound1_MPI:bound2_MPI))
          CALL MPI_Recv(psi(i)%RvecG(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,   &
                        MPI_Real8,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
        ENDIF
      ENDIF
    ENDDO ! for i=ndim1,ndim2
  ENDIF ! for MPI_id/=0  

END SUBROUTINE distribute_psi_MPI
!=======================================================================================

!=======================================================================================
!> subroutine for distribute psi(ndim1:ndim2) to different threads 
!> the distribution is depends on the length of vec (RvecB, CvecB, RvecG, or CvecG)
!> vector is packed on master and unpacked on theards to reduce comm. time.
!=======================================================================================
SUBROUTINE distribute_psi_pack_MPI(psi,ndim1,ndim2,With_Grid)
  USE mod_system
  USE mod_psi_set_alloc
  USE mod_MPI
  USE mod_MPI_Aid
  IMPLICIT NONE
  
  TYPE(param_psi), intent(inout)          :: psi(:) 
  Integer,         intent(in)             :: ndim1
  Integer,         intent(in)             :: ndim2
  Logical,optional,intent(in)             :: With_Grid
 
  Complex(kind=Rkind),allocatable         :: Cvec(:)
  Real(kind=Rkind),allocatable            :: Rvec(:)
  Integer                                 :: length(MPI_np)
  Integer                                 :: count
  Integer                                 :: i,j
  Logical                                 :: With_Grid_loc

  With_Grid_loc=.FALSE.
  IF (present(With_Grid)) With_Grid_loc=With_Grid

  DO i_MPI=1,MPI_np-1
    length(i_MPI)=0
    DO i=ndim1,ndim2
      IF(.NOT. With_Grid_loc) THEN
        nb_per_MPI=psi(i)%nb_tot/MPI_np
        nb_rem_MPI=mod(psi(i)%nb_tot,MPI_np) 
      ELSE
        nb_per_MPI=psi(i)%nb_qaie/MPI_np
        nb_rem_MPI=mod(psi(i)%nb_qaie,MPI_np) 
      ENDIF
      bound1_MPI=i_MPI*nb_per_MPI+1+MIN(i_MPI,nb_rem_MPI)
      bound2_MPI=(i_MPI+1)*nb_per_MPI+MIN(i_MPI,nb_rem_MPI)+merge(1,0,nb_rem_MPI>i_MPI)
      
      length(i_MPI)=length(i_MPI)+(bound2_MPI-bound1_MPI+1)
    ENDDO
  ENDDO

  IF(MPI_id==0) THEN
    DO i_MPI=1,MPI_np-1
      IF(psi(1)%cplx) THEN
        allocate(Cvec(length(i_MPI)))
      ELSE
        allocate(Rvec(length(i_MPI)))
      ENDIF
      count=1

      DO i=ndim1,ndim2
        IF(.NOT. With_Grid_loc) THEN
          nb_per_MPI=psi(i)%nb_tot/MPI_np
          nb_rem_MPI=mod(psi(i)%nb_tot,MPI_np) !remainder jobs
        ELSE
          nb_per_MPI=psi(i)%nb_qaie/MPI_np
          nb_rem_MPI=mod(psi(i)%nb_qaie,MPI_np) !remainder jobs
        ENDIF
        bound1_MPI=i_MPI*nb_per_MPI+1+MIN(i_MPI,nb_rem_MPI)
        bound2_MPI=(i_MPI+1)*nb_per_MPI+MIN(i_MPI,nb_rem_MPI)                          &
                                       +merge(1,0,nb_rem_MPI>i_MPI)

        ! pack array
        IF (.NOT. With_Grid_loc) THEN
          IF (psi(i)%cplx) THEN
            Cvec(count:bound2_MPI-bound1_MPI+count)=psi(i)%CvecB(bound1_MPI:bound2_MPI)
          ELSE
            Rvec(count:bound2_MPI-bound1_MPI+count)=psi(i)%RvecB(bound1_MPI:bound2_MPI)
          ENDIF
        ELSE
          IF (psi(i)%cplx) THEN
            Cvec(count:bound2_MPI-bound1_MPI+count)=psi(i)%CvecG(bound1_MPI:bound2_MPI)
          ELSE
            Rvec(count:bound2_MPI-bound1_MPI+count)=psi(i)%RvecG(bound1_MPI:bound2_MPI)
          ENDIF
        ENDIF ! for .NOT. With_Grid_loc 
        count=bound2_MPI-bound1_MPI+count+1
        
      ENDDO ! for i=ndim1,ndim2 

      ! send array
      IF(psi(1)%cplx) THEN
        CALL MPI_Send(Cvec,count-1,MPI_Complex8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
        deallocate(Cvec)
      ELSE
        CALL MPI_Send(Rvec,count-1,MPI_Real8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
        deallocate(Rvec)
      ENDIF
    ENDDO ! for i_MPI=1,MPI_np-1
  ENDIF ! for MPI_id==0
  
  ! MPI/=0------------------------------------------------------------------------------
  IF(MPI_id/=0) THEN
    ! receive array
    IF(psi(1)%cplx) THEN
      allocate(Cvec(length(MPI_id)))
      CALL MPI_Recv(Cvec,length(MPI_id),MPI_Complex8,root_MPI,MPI_id,                  &
                    MPI_COMM_WORLD,MPI_stat,MPI_err)
    ELSE
      allocate(Rvec(length(MPI_id)))
      CALL MPI_Recv(Rvec,length(MPI_id),MPI_Real8,root_MPI,MPI_id,                     &
                    MPI_COMM_WORLD,MPI_stat,MPI_err)
    ENDIF

    count=1
    DO i=ndim1,ndim2
      IF(.NOT. With_Grid_loc) THEN
        nb_per_MPI=psi(i)%nb_tot/MPI_np
        nb_rem_MPI=mod(psi(i)%nb_tot,MPI_np) !remainder jobs
      ELSE
        nb_per_MPI=psi(i)%nb_qaie/MPI_np
        nb_rem_MPI=mod(psi(i)%nb_qaie,MPI_np) !remainder jobs
      ENDIF
      bound1_MPI=MPI_id*nb_per_MPI+1+MIN(MPI_id,nb_rem_MPI)
      bound2_MPI=(MPI_id+1)*nb_per_MPI+MIN(MPI_id,nb_rem_MPI)                          &
                                      +merge(1,0,nb_rem_MPI>MPI_id)

      ! unpack
      IF (.NOT. With_Grid_loc) THEN
        IF (psi(i)%cplx) THEN
          IF(.NOT. allocated(psi(i)%CvecB))                                            &
                    allocate(psi(i)%CvecB(bound1_MPI:bound2_MPI))
          psi(i)%CvecB(bound1_MPI:bound2_MPI)=Cvec(count:count+bound2_MPI-bound1_MPI)
        ELSE
          IF(.NOT. allocated(psi(i)%RvecB))                                            &
                    allocate(psi(i)%RvecB(bound1_MPI:bound2_MPI))
          psi(i)%RvecB(bound1_MPI:bound2_MPI)=Rvec(count:count+bound2_MPI-bound1_MPI)
        ENDIF
      ELSE
        IF (psi(i)%cplx) THEN
          IF(.NOT. allocated(psi(i)%CvecG))                                            &
                    allocate(psi(i)%CvecG(bound1_MPI:bound2_MPI))
          psi(i)%CvecG(bound1_MPI:bound2_MPI)=Cvec(count:count+bound2_MPI-bound1_MPI)
        ELSE
          IF(.NOT. allocated(psi(i)%RvecG))                                            &
                    allocate(psi(i)%RvecG(bound1_MPI:bound2_MPI))
          psi(i)%RvecG(bound1_MPI:bound2_MPI)=Rvec(count:count+bound2_MPI-bound1_MPI)
        ENDIF
      ENDIF
      count=bound2_MPI-bound1_MPI+count+1
    ENDDO ! for i=ndim1,ndim2
    
    IF(psi(1)%cplx) THEN
      deallocate(Cvec)
    ELSE
      deallocate(Rvec)
    ENDIF
  ENDIF ! for MPI_id/=0  

END SUBROUTINE distribute_psi_pack_MPI
!=======================================================================================

!=======================================================================================
!> for calculating the overlap of <psi|psi> or <psi|Hpsi> with MPI
!> working for ndim1=<i<=ndim2 
!
!>   NOTE: the pre-distribution of vectors (RvecB, CvecB, RvecG, or CvecG) are required
!>         i.e. subroutine distribute_psi_MPI
!>   NOTE: works according to the vectors stored on different threads 
!>         limited by bound1_MPI and bound2_MPI, see "Overlap_psipsi_MPI3"
!=======================================================================================
SUBROUTINE calculate_overlap_MPI(psi,ndim1,ndim2,With_Grid,Hpsi,S_overlap,H_overlap)
  USE mod_system
  USE mod_psi_set_alloc
  USE mod_MPI
  USE mod_MPI_Aid
  IMPLICIT NONE
  
  TYPE(param_psi),          intent(in)    :: psi(:) 
  TYPE(param_psi),optional, intent(in)    :: Hpsi(:) 
  Integer,                  intent(in)    :: ndim1
  Integer,                  intent(in)    :: ndim2
  Logical,optional,         intent(in)    :: With_Grid
  Real(kind=Rkind),optional,intent(inout) :: H_overlap(:,:)
  Real(kind=Rkind),optional,intent(inout) :: S_overlap(:,:)

  Complex(kind=Rkind)                     :: Overlap
  Integer                                 :: i,j
  Logical                                 :: With_Grid_loc

  With_Grid_loc=.FALSE.
  IF(present(With_Grid)) With_Grid_loc=With_Grid

  IF((present(Hpsi) .AND. (.NOT. present(H_overlap))) .OR.                             &
     ((.NOT. present(H_overlap)) .AND. (.NOT. present(S_overlap)))) THEN
    write(out_unitp,*) 'variable presented error in calculate_overlap_MPI'  
    STOP
  ENDIF
  
  !-overlap <psi|Hpsi>------------------------------------------------------------------
  IF(present(Hpsi) .AND. present(H_overlap)) THEN 
    DO i=1,ndim1-1
      DO j=ndim1,ndim2
        CALL Overlap_psipsi_MPI3(Overlap,psi(j),Hpsi(i),With_Grid=With_Grid)
        H_overlap(j,i)=real(Overlap,kind=Rkind) 
      ENDDO 
    END DO
    
    DO i=ndim1,ndim2
      DO j=1,ndim2
        CALL Overlap_psipsi_MPI3(Overlap,psi(j),Hpsi(i),With_Grid=With_Grid)
        H_overlap(j,i)=real(Overlap,kind=Rkind) 
      ENDDO 
    ENDDO
    
    ! collect and broadcast result
    CALL MPI_Reduce_sum_matrix(H_overlap,ndim1,ndim2,1,ndim1-1,root_MPI)
    CALL MPI_Bcast_matrix     (H_overlap,ndim1,ndim2,1,ndim1-1,root_MPI)
    CALL MPI_Reduce_sum_matrix(H_overlap,1,ndim2,ndim1,ndim2,root_MPI)
    CALL MPI_Bcast_matrix     (H_overlap,1,ndim2,ndim1,ndim2,root_MPI)
  ENDIF 
    
  !-overlap <psi|psi>-------------------------------------------------------------------  
  IF(present(S_overlap)) THEN 
    DO i=1,ndim1-1
      DO j=ndim1,ndim2
        IF(j<=i) THEN
          CALL Overlap_psipsi_MPI3(Overlap,psi(j),psi(i),With_Grid=With_Grid)
          S_overlap(j,i)=real(Overlap,kind=Rkind)
          S_overlap(i,j)=S_overlap(j,i)
        ENDIF
      ENDDO 
    ENDDO 
    
    DO i=ndim1,ndim2
      DO j=1,ndim2
        IF(j<=i) THEN
          CALL Overlap_psipsi_MPI3(Overlap,psi(j),psi(i),With_Grid=With_Grid)
          S_overlap(j,i)=real(Overlap,kind=Rkind)
          S_overlap(i,j)=S_overlap(j,i)
        ENDIF
      ENDDO 
    ENDDO
    
    ! collect and broadcast result 
    CALL MPI_Reduce_sum_matrix(S_overlap,ndim1,ndim2,1,ndim1-1,root_MPI)
    CALL MPI_Bcast_matrix     (S_overlap,ndim1,ndim2,1,ndim1-1,root_MPI)
    CALL MPI_Reduce_sum_matrix(S_overlap,1,ndim2,ndim1,ndim2,root_MPI)
    CALL MPI_Bcast_matrix     (S_overlap,1,ndim2,ndim1,ndim2,root_MPI)
  ENDIF

END SUBROUTINE calculate_overlap_MPI
!=======================================================================================

!=======================================================================================
! calculate overlap for <psi(n+1)|psi(i)> or/and <Hpsi(n+1)|psi(i)> i=1...n
!=======================================================================================
SUBROUTINE calculate_overlap1D_MPI(psi,ndim,With_Grid,Hpsi,S_overlap1D,H_overlap1D)
  USE mod_system
  USE mod_psi_set_alloc
  USE mod_MPI
  USE mod_MPI_Aid
  IMPLICIT NONE
  
  TYPE(param_psi),          intent(in)    :: psi(:) 
  TYPE(param_psi),optional, intent(in)    :: Hpsi(:) 
  Integer,                  intent(in)    :: ndim
  Logical,optional,         intent(in)    :: With_Grid
  Real(kind=Rkind),optional,intent(inout) :: H_overlap1D(:)
  Real(kind=Rkind),optional,intent(inout) :: S_overlap1D(:)
  
  Complex(kind=Rkind)                     :: Overlap
  Real(kind=Rkind)                        :: Overlap1D_temp(ndim)
  Integer                                 :: i,j
  Logical                                 :: With_Grid_loc

  With_Grid_loc=.FALSE.
  IF(present(With_Grid)) With_Grid_loc=With_Grid

  IF((present(Hpsi) .AND. (.NOT. present(H_overlap1D))) .OR.                             &
     ((.NOT. present(S_overlap1D)) .AND. (.NOT. present(H_overlap1D)))) THEN
    write(out_unitp,*) 'variable presented error in calculate_overlap_MPI'  
    STOP
  ENDIF
  
  !-overlap <psi|Hpsi>------------------------------------------------------------------
  IF(present(Hpsi) .AND. present(H_overlap1D)) THEN 
    DO i=1,ndim
      CALL Overlap_psipsi_MPI3(Overlap,psi(i),Hpsi(ndim),With_Grid=With_Grid)
      H_overlap1D(i)=real(Overlap,kind=Rkind) 
    END DO
    
    ! collect and broadcast result
    CALL MPI_Reduce(H_overlap1D,Overlap1D_temp,ndim,MPI_Real8,MPI_SUM,root_MPI,        &
                    MPI_COMM_WORLD,MPI_err)
    IF(MPI_id==0) H_overlap1D=Overlap1D_temp
    CALL MPI_BCAST(H_overlap1D,ndim,MPI_Real8,root_MPI,MPI_COMM_WORLD,MPI_err)
  ENDIF 
    
  !-overlap <psi|psi>-------------------------------------------------------------------  
  IF(present(S_overlap1D)) THEN 
    DO i=1,ndim
      CALL Overlap_psipsi_MPI3(Overlap,psi(i),psi(ndim),With_Grid=With_Grid)
      S_overlap1D(i)=real(Overlap,kind=Rkind)
    ENDDO 
    
    ! collect and broadcast result
    CALL MPI_Reduce(S_overlap1D,Overlap1D_temp,ndim,MPI_Real8,MPI_SUM,root_MPI,        &
                    MPI_COMM_WORLD,MPI_err)
    IF(MPI_id==0) S_overlap1D=Overlap1D_temp
    CALL MPI_BCAST(S_overlap1D,ndim,MPI_Real8,root_MPI,MPI_COMM_WORLD,MPI_err)
  ENDIF

END SUBROUTINE calculate_overlap1D_MPI
!=======================================================================================

!=======================================================================================
!> subroutine for the calculation of Overlap_psi1_psi2 with MPI 
!> overlap is cut into submatrix for MPI parallel   
!> note, this requre new  Overlap_psi1_psi2 routine
!
!> the submatrix of psi and Hpsi are keeped also for the calculation 
!> of Residual g in MakeResidual_Davidson_MPI2
!=======================================================================================
#if(run_MPI)
SUBROUTINE Overlap_psi1_psi2_MPI5(H_overlap,S_overlap,psi,Hpsi,ndim0,ndim,With_Grid)
  USE mod_system
  USE mod_psi_set_alloc
  USE mod_MPI
  USE mod_MPI_Aid
  IMPLICIT NONE
  
  TYPE(param_psi), intent(inout)          :: psi(:)  !< on non-root threads allocated
  TYPE(param_psi), intent(inout)          :: Hpsi(:) !< on non-root threads allocated
  Real(kind=Rkind),intent(inout)          :: H_overlap(:,:)
  Real(kind=Rkind),intent(inout)          :: S_overlap(:,:)
  Logical,optional,intent(in)             :: With_Grid
  Integer,         intent(in)             :: ndim0
  Integer,         intent(in)             :: ndim

  complex(kind=Rkind)                     :: Overlap
  logical                                 :: With_Grid_loc
  Integer                                 :: i
  Integer                                 :: j
  Logical                                 :: root_jobs
  
  With_Grid_loc=.FALSE.
  IF (present(With_Grid)) With_Grid_loc = With_Grid

  ! only on master: psi%CvecG,psi%RvecG,psi%CvecB,psi%RvecB
  CALL distribute_psi_pack_MPI( psi,ndim0+1,ndim,With_Grid)
  CALL distribute_psi_pack_MPI(Hpsi,ndim0+1,ndim,With_Grid)

  CALL calculate_overlap_MPI(psi,ndim0+1,ndim,With_Grid,Hpsi=Hpsi,                     &
                             S_overlap=S_overlap,H_overlap=H_overlap)
  
END SUBROUTINE Overlap_psi1_psi2_MPI5
#endif
!=======================================================================================

!=======================================================================================
!> subroutine for the calculation of Overlap_psi1_psi2 with MPI 
!> overlap is cut into submatrix for MPI parallel   
!> note, this requre new  Overlap_psi1_psi2 routine
!
!> the submatrix of psi and Hpsi are keeped also for the calculation 
!> of Residual g in MakeResidual_Davidson_MPI2
!=======================================================================================
#if(run_MPI)
SUBROUTINE Overlap_psi1_psi2_MPI3(H_overlap,S_overlap,psi,Hpsi,ndim0,ndim,With_Grid)
  USE mod_system
  USE mod_psi_set_alloc
  USE mod_MPI
  USE mod_MPI_Aid
  IMPLICIT NONE
  
  TYPE(param_psi), intent(inout)          :: psi(:)  !< on non-root threads allocated
  TYPE(param_psi), intent(inout)          :: Hpsi(:) !< on non-root threads allocated
  Real(kind=Rkind),intent(inout)          :: H_overlap(:,:)
  Real(kind=Rkind),intent(inout)          :: S_overlap(:,:)
  Logical,optional,intent(in)             :: With_Grid
  Integer,         intent(in)             :: ndim0
  Integer,         intent(in)             :: ndim

  complex(kind=Rkind)                     :: Overlap
  logical                                 :: With_Grid_loc
  Integer                                 :: i
  Integer                                 :: j
  Logical                                 :: root_jobs
  
  With_Grid_loc=.FALSE.
  IF (present(With_Grid)) With_Grid_loc = With_Grid

  ! only on master: psi%CvecG,psi%RvecG,psi%CvecB,psi%RvecB
  IF(MPI_id==0) THEN
    DO i=ndim0+1,ndim
      
      IF(.NOT. With_Grid_loc) THEN
        nb_per_MPI=Hpsi(i)%nb_tot/MPI_np
        nb_rem_MPI=mod(Hpsi(i)%nb_tot,MPI_np) !remainder jobs
      ELSE
        nb_per_MPI=Hpsi(i)%nb_qaie/MPI_np
        nb_rem_MPI=mod(Hpsi(i)%nb_qaie,MPI_np) !remainder jobs
      ENDIF
      
      ! Send array
      DO i_MPI=1,MPI_np-1
        bound1_MPI=i_MPI*nb_per_MPI+1+MIN(i_MPI,nb_rem_MPI)
        bound2_MPI=(i_MPI+1)*nb_per_MPI+MIN(i_MPI,nb_rem_MPI)+merge(1,0,nb_rem_MPI>i_MPI)
        IF (.NOT. With_Grid_loc) THEN
          IF (psi(i)%cplx) THEN
            CALL MPI_Send(Hpsi(i)%CvecB(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,&
                          MPI_complex8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
            CALL MPI_Send(psi(i)%CvecB(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1, &
                          MPI_Complex8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
          ELSE
            CALL MPI_Send(Hpsi(i)%RvecB(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,&
                          MPI_Real8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
            CALL MPI_Send(psi(i)%RvecB(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1, &
                          MPI_Real8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
          
          ENDIF
        ELSE
          IF (psi(i)%cplx) THEN
            CALL MPI_Send(Hpsi(i)%CvecG(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,&
                          MPI_complex8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
            CALL MPI_Send(psi(i)%CvecG(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1, &
                          MPI_Complex8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)        
          ELSE
            CALL MPI_Send(Hpsi(i)%RvecG(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,&
                          MPI_Real8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
            CALL MPI_Send(psi(i)%RvecG(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1, &
                          MPI_Real8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
          
          ENDIF
        ENDIF ! for .NOT. With_Grid_loc 
      ENDDO ! for i_MPI=1,MPI_np-1
    ENDDO ! for i=ndim0+1,ndim

    !> calcuation on master
    DO i=1,ndim0
      DO j=ndim0+1,ndim
        CALL Overlap_psipsi_MPI3(Overlap,psi(j),Hpsi(i),With_Grid=With_Grid)
        H_overlap(j,i)=real(Overlap,kind=Rkind) 
        CALL Overlap_psipsi_MPI3(Overlap,psi(j), psi(i),With_Grid=With_Grid)
        S_overlap(j,i)=real(Overlap,kind=Rkind)
      ENDDO ! for j=ndim0+1,ndim
    END DO ! for i=ndim0+1,ndim
    
    DO i=ndim0+1,ndim
      DO j=1,ndim
        CALL Overlap_psipsi_MPI3(Overlap,psi(j),Hpsi(i),With_Grid=With_Grid)
        H_overlap(j,i)=real(Overlap,kind=Rkind) 
        CALL Overlap_psipsi_MPI3(Overlap,psi(j), psi(i),With_Grid=With_Grid)
        S_overlap(j,i)=real(Overlap,kind=Rkind)
      ENDDO ! for j=ndim0+1,ndim
    END DO ! for i=ndim0+1,ndim

  ENDIF ! for MPI_id==0

  ! MPI/=0------------------------------------------------------------------------------
  IF(MPI_id/=0) THEN
    DO i=ndim0+1,ndim
      IF(.NOT. With_Grid_loc) THEN
        nb_per_MPI=Hpsi(i)%nb_tot/MPI_np
        nb_rem_MPI=mod(Hpsi(i)%nb_tot,MPI_np) !remainder jobs
      ELSE
        nb_per_MPI=Hpsi(i)%nb_qaie/MPI_np
        nb_rem_MPI=mod(Hpsi(i)%nb_qaie,MPI_np) !remainder jobs
      ENDIF
      bound1_MPI=MPI_id*nb_per_MPI+1+MIN(MPI_id,nb_rem_MPI)
      bound2_MPI=(MPI_id+1)*nb_per_MPI+MIN(MPI_id,nb_rem_MPI)+merge(1,0,nb_rem_MPI>MPI_id)
      
      ! receive array
      IF (.NOT. With_Grid_loc) THEN
        IF (psi(i)%cplx) THEN
          IF(.NOT. allocated(Hpsi(i)%CvecB))                                           &
              allocate(Hpsi(i)%CvecB(bound1_MPI:bound2_MPI))
          CALL MPI_Recv(Hpsi(i)%CvecB(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,  &
                        MPI_complex8,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
          IF(.NOT. allocated(psi(i)%CvecB))                                            &
               allocate(psi(i)%CvecB(bound1_MPI:bound2_MPI))
          CALL MPI_Recv( psi(i)%CvecB(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,  &
                        MPI_Complex8,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
        ELSE
          IF(.NOT. allocated(Hpsi(i)%RvecB))                                           &
               allocate(Hpsi(i)%RvecB(bound1_MPI:bound2_MPI))
          CALL MPI_Recv(Hpsi(i)%RvecB(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,  &
                        MPI_Real8,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
        
          IF(.NOT. allocated(psi(i)%RvecB))                                            &
               allocate(psi(i)%RvecB(bound1_MPI:bound2_MPI))
          CALL MPI_Recv( psi(i)%RvecB(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,  &
                        MPI_Real8,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
          
        ENDIF
      ELSE
        IF (psi(i)%cplx) THEN
          IF(.NOT. allocated(Hpsi(i)%CvecG))                                           &
               allocate(Hpsi(i)%CvecG(bound1_MPI:bound2_MPI))
          CALL MPI_Recv(Hpsi(i)%CvecG(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,  &
                        MPI_complex8,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
          IF(.NOT. allocated(psi(i)%CvecG))                                            &
               allocate(psi(i)%CvecG(bound1_MPI:bound2_MPI))
          CALL MPI_Recv( psi(i)%CvecG(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,  &
                        MPI_Complex8,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
        ELSE
          IF(.NOT. allocated(Hpsi(i)%RvecB))                                           &
               allocate(Hpsi(i)%RvecG(bound1_MPI:bound2_MPI))
          CALL MPI_Recv(Hpsi(i)%RvecB(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,  &
                        MPI_Real8,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
          IF(.NOT. allocated(psi(i)%RvecB))                                            &
               allocate(psi(i)%RvecG(bound1_MPI:bound2_MPI))
          CALL MPI_Recv( psi(i)%RvecB(bound1_MPI:bound2_MPI),bound2_MPI-bound1_MPI+1,  &
                        MPI_Real8,root_MPI,MPI_id,MPI_COMM_WORLD,MPI_stat,MPI_err)
        ENDIF
      ENDIF
    ENDDO ! for i=ndim0+1,ndim

    !> calcuation on each thread
    DO i=1,ndim0
      DO j=ndim0+1,ndim
        CALL Overlap_psipsi_MPI3(Overlap,psi(j),Hpsi(i),With_Grid=With_Grid)
        H_overlap(j,i)=real(Overlap,kind=Rkind) 
        CALL Overlap_psipsi_MPI3(Overlap,psi(j), psi(i),With_Grid=With_Grid)
        S_overlap(j,i)=real(Overlap,kind=Rkind)
      ENDDO ! for j=ndim0+1,ndim
    END DO ! for i=ndim0+1,ndim
    
    DO i=ndim0+1,ndim
      DO j=1,ndim
        CALL Overlap_psipsi_MPI3(Overlap,psi(j),Hpsi(i),With_Grid=With_Grid)
        H_overlap(j,i)=real(Overlap,kind=Rkind) 
        CALL Overlap_psipsi_MPI3(Overlap,psi(j), psi(i),With_Grid=With_Grid)
        S_overlap(j,i)=real(Overlap,kind=Rkind)
      ENDDO ! for j=ndim0+1,ndim
    END DO ! for i=ndim0+1,ndim

  ENDIF ! for MPI_id/=0  

  CALL MPI_Reduce_sum_matrix(H_overlap,ndim0+1,ndim,1,ndim0,root_MPI)
  CALL MPI_Bcast_matrix     (H_overlap,ndim0+1,ndim,1,ndim0,root_MPI)
  CALL MPI_Reduce_sum_matrix(H_overlap,1,ndim,ndim0+1,ndim,root_MPI)
  CALL MPI_Bcast_matrix     (H_overlap,1,ndim,ndim0+1,ndim,root_MPI)

  CALL MPI_Reduce_sum_matrix(S_overlap,ndim0+1,ndim,1,ndim0,root_MPI)
  CALL MPI_Bcast_matrix     (S_overlap,ndim0+1,ndim,1,ndim0,root_MPI)
  CALL MPI_Reduce_sum_matrix(S_overlap,1,ndim,ndim0+1,ndim,root_MPI)
  CALL MPI_Bcast_matrix     (S_overlap,1,ndim,ndim0+1,ndim,root_MPI)

END SUBROUTINE Overlap_psi1_psi2_MPI3
#endif
!=======================================================================================

!=======================================================================================
! subroutine for the calculation of Overlap_psi1_psi2 with MPI 
! ONLY the overlap of i=j are calculated on the other threads to reduce memory 
! psi and Hpsi are keeped also for the calculation 
! of Residual g in MakeResidual_Davidson_MPI2
!=======================================================================================
#if(run_MPI)
SUBROUTINE Overlap_psi1_psi2_MPI2(H_overlap,S_overlap,psi,Hpsi,ndim,With_Grid)
  USE mod_system
  USE mod_psi_set_alloc
  USE mod_MPI
  USE mod_MPI_Aid
  IMPLICIT NONE
  
  TYPE(param_psi), intent(inout)          :: psi(:)  !< on non-root threads allocated
  TYPE(param_psi), intent(inout)          :: Hpsi(:) !< on non-root threads allocated
  Real(kind=Rkind),intent(inout)          :: H_overlap(:,:)
  Real(kind=Rkind),intent(inout)          :: S_overlap(:,:)
  Logical,optional,intent(in)             :: With_Grid
  Integer,         intent(in)             :: ndim

  complex(kind=Rkind)                     :: Overlap
  logical                                 :: With_Grid_loc
  Integer                                 :: i
  Integer                                 :: j
  Logical                                 :: send_once(MPI_np-1)
  Logical                                 :: root_jobs
  
  With_Grid_loc=.FALSE.
  IF (present(With_Grid)) With_Grid_loc = With_Grid

  ! only on master: psi%CvecG,psi%RvecG,psi%CvecB,psi%RvecB
  nb_per_MPI=(ndim)/MPI_np
  nb_rem_MPI=mod(ndim,MPI_np) !remainder jobs
  IF(MPI_id==0) THEN
    !> send Hpsi and psi
    DO i_MPI=1,MPI_np-1
      bound1_MPI=i_MPI*nb_per_MPI+1+MIN(i_MPI,nb_rem_MPI)
      bound2_MPI=(i_MPI+1)*nb_per_MPI+MIN(i_MPI,nb_rem_MPI)+merge(1,0,nb_rem_MPI>i_MPI)
      DO i=1,ndim
        IF(i>=bound1_MPI .AND. i<=bound2_MPI) THEN
          IF (.NOT. With_Grid_loc) THEN
            IF (psi(i)%cplx) THEN
              CALL MPI_Send(Hpsi(i)%CvecB,Hpsi(i)%nb_tot,MPI_complex8,i_MPI,           &
                            i_MPI,MPI_COMM_WORLD,MPI_err)
              CALL MPI_Send(psi(i)%CvecB,psi(i)%nb_tot,MPI_Complex8,i_MPI,             &
                            i_MPI,MPI_COMM_WORLD,MPI_err)

            ELSE
              CALL MPI_Send(Hpsi(i)%RvecB,Hpsi(i)%nb_tot,MPI_Real8,i_MPI,              &
                            i_MPI,MPI_COMM_WORLD,MPI_err)
              CALL MPI_Send(psi(i)%RvecB,psi(i)%nb_tot,MPI_Real8,i_MPI,                &
                            i_MPI,MPI_COMM_WORLD,MPI_err)
            ENDIF
          ELSE
            IF (psi(i)%cplx) THEN
              CALL MPI_Send(Hpsi(i)%CvecG,Hpsi(i)%nb_qaie,MPI_complex8,i_MPI,          &
                            i_MPI,MPI_COMM_WORLD,MPI_err)
              CALL MPI_Send(psi(i)%CvecG,psi(i)%nb_qaie,MPI_Complex8,i_MPI,            &
                            i_MPI,MPI_COMM_WORLD,MPI_err)
            ELSE
              CALL MPI_Send(Hpsi(i)%RvecG,Hpsi(i)%nb_qaie,MPI_Real8,i_MPI,             &
                            i_MPI,MPI_COMM_WORLD,MPI_err)
              CALL MPI_Send(psi(i)%RvecG,psi(i)%nb_qaie,MPI_Real8,i_MPI,               &
                            i_MPI,MPI_COMM_WORLD,MPI_err)
            ENDIF
          ENDIF ! for .NOT. With_Grid_loc 
        ENDIF ! for i>=bound1_MPI .AND. i<=bound2_MPI
      ENDDO ! for i=1,ndim
    ENDDO ! for i_MPI=1,MPI_np-1

    !> calcuation on master
    DO i=1,ndim
      DO j=1,ndim
        root_jobs=.TRUE.
        DO i_MPI=1,MPI_np-1
          bound1_MPI=i_MPI*nb_per_MPI+1+MIN(i_MPI,nb_rem_MPI)
          bound2_MPI=(i_MPI+1)*nb_per_MPI+MIN(i_MPI,nb_rem_MPI)                        &
                                         +merge(1,0,nb_rem_MPI>i_MPI)
          IF((i>=bound1_MPI .AND. i<=bound2_MPI) .AND.                                 &
             (j>=bound1_MPI .AND. j<=bound2_MPI)) THEN
             root_jobs=.FALSE.
          ENDIF
        ENDDO ! i_MPI=1,MPI_np-1
        
        IF(root_jobs) THEN
          CALL Overlap_psipsi_MPI(Overlap,psi(j),Hpsi(i),With_Grid=With_Grid)
          H_overlap(j,i)=real(Overlap,kind=Rkind) 
          CALL Overlap_psipsi_MPI(Overlap,psi(j), psi(i),With_Grid=With_Grid)
          S_overlap(j,i)=real(Overlap,kind=Rkind)
        ENDIF
      ENDDO ! for j=1,ndim
    END DO ! for i=1,ndim

    DO i_MPI=1,MPI_np-1
      bound1_MPI=i_MPI*nb_per_MPI+1+MIN(i_MPI,nb_rem_MPI)
      bound2_MPI=(i_MPI+1)*nb_per_MPI+MIN(i_MPI,nb_rem_MPI)+merge(1,0,nb_rem_MPI>i_MPI)    
      CALL MPI_Recv_matrix(H_overlap,bound1_MPI,bound2_MPI,bound1_MPI,bound2_MPI,      &  
                           i_MPI,i_MPI)
      CALL MPI_Recv_matrix(S_overlap,bound1_MPI,bound2_MPI,bound1_MPI,bound2_MPI,      &  
                           i_MPI,i_MPI)
    ENDDO 

  ENDIF ! for MPI_id==0

  !-------------------------------------------------------------------------------------
  bound1_MPI=MPI_id*nb_per_MPI+1+MIN(MPI_id,nb_rem_MPI)
  bound2_MPI=(MPI_id+1)*nb_per_MPI+MIN(MPI_id,nb_rem_MPI)+merge(1,0,nb_rem_MPI>MPI_id)
  IF(MPI_id/=0) THEN
    DO i=1,ndim
      IF(i>=bound1_MPI .AND. i<=bound2_MPI) THEN
        IF (.NOT. With_Grid_loc) THEN
          IF (psi(i)%cplx) THEN
            IF(.NOT. allocated(Hpsi(i)%CvecB)) CALL alloc_NParray(Hpsi(i)%CvecB,       &
                                            (/Hpsi(i)%nb_tot/),'Hpsi%CvecB','alloc_psi')
            CALL MPI_Recv(Hpsi(i)%CvecB,Hpsi(i)%nb_tot,MPI_complex8,root_MPI,MPI_id,   &
                          MPI_COMM_WORLD,MPI_stat,MPI_err)
            IF(.NOT. allocated(psi(i)%CvecB))  CALL alloc_NParray(psi(i)%CvecB,        &
                                              (/psi(i)%nb_tot/),'psi%CvecB','alloc_psi')
            CALL MPI_Recv( psi(i)%CvecB, psi(i)%nb_tot,MPI_Complex8,root_MPI,MPI_id,   &
                          MPI_COMM_WORLD,MPI_stat,MPI_err)
          ELSE
            IF(.NOT. allocated(Hpsi(i)%RvecB)) CALL alloc_NParray(Hpsi(i)%RvecB,       &
                                            (/Hpsi(i)%nb_tot/),'Hpsi%RvecB','alloc_psi')
            CALL MPI_Recv(Hpsi(i)%RvecB,Hpsi(i)%nb_tot,MPI_Real8,root_MPI,MPI_id,      &
                          MPI_COMM_WORLD,MPI_stat,MPI_err)
          
            IF(.NOT. allocated(psi(i)%RvecB))  CALL alloc_NParray(psi(i)%RvecB,        &
                                              (/psi(i)%nb_tot/),'psi%RvecB','alloc_psi')
            CALL MPI_Recv( psi(i)%RvecB, psi(i)%nb_tot,MPI_Real8,root_MPI,MPI_id,      &
                          MPI_COMM_WORLD,MPI_stat,MPI_err)
          ENDIF
        ELSE
          IF (psi(i)%cplx) THEN
            IF(.NOT. allocated(Hpsi(i)%CvecG)) CALL alloc_NParray(Hpsi(i)%CvecG,       &
                                           (/Hpsi(i)%nb_qaie/),'Hpsi%CvecG','alloc_psi')
            CALL MPI_Recv(Hpsi(i)%CvecG,Hpsi(i)%nb_qaie,MPI_complex8,root_MPI,MPI_id,  &
                          MPI_COMM_WORLD,MPI_stat,MPI_err)
            IF(.NOT. allocated(psi(i)%CvecG))  CALL alloc_NParray(psi(i)%CvecG,        &
                                             (/psi(i)%nb_qaie/),'psi%CvecG','alloc_psi')
            CALL MPI_Recv( psi(i)%CvecG, psi(i)%nb_qaie,MPI_Complex8,root_MPI,MPI_id,  &
                        MPI_COMM_WORLD,MPI_stat,MPI_err)
          ELSE
            IF(.NOT. allocated(Hpsi(i)%RvecB)) CALL alloc_NParray(Hpsi(i)%RvecG,       &
                                            (/Hpsi(i)%nb_tot/),'Hpsi%RvecG','alloc_psi')
            CALL MPI_Recv(Hpsi(i)%RvecB,Hpsi(i)%nb_tot,MPI_Real8,root_MPI,MPI_id,      &
                          MPI_COMM_WORLD,MPI_stat,MPI_err)
            IF(.NOT. allocated(psi(i)%RvecB))  CALL alloc_NParray(psi(i)%RvecG,        &
                                              (/psi(i)%nb_tot/),'psi%RvecG','alloc_psi')
            CALL MPI_Recv( psi(i)%RvecB, psi(i)%nb_tot,MPI_Real8,root_MPI,MPI_id,      &
                        MPI_COMM_WORLD,MPI_stat,MPI_err)
          ENDIF
        ENDIF ! for .NOT. With_Grid_loc
      ENDIF ! for i>=bound1_MPI .AND. i<=bound2_MPI
    ENDDO ! for i=1,ndim

    ! main calculation
    DO i=1,ndim
      DO j=1,ndim
        IF((i>=bound1_MPI .AND. i<=bound2_MPI) .AND.                                   &
           (j>=bound1_MPI .AND. j<=bound2_MPI)) THEN
          IF (.NOT. With_Grid_loc) THEN
            IF (psi(j)%symab > -1 .AND. Hpsi(i)%symab > -1                             &
                                  .AND. Hpsi(i)%symab /= Hpsi(i)%symab) THEN
            ELSE
              CALL Overlap_psipsi_MPI(Overlap,psi(j),Hpsi(i),With_Grid=With_Grid)
              H_overlap(j,i)=real(Overlap,kind=Rkind)
              CALL Overlap_psipsi_MPI(Overlap,psi(j), psi(i),With_Grid=With_Grid)
              S_overlap(j,i)=real(Overlap,kind=Rkind)
            ENDIF
          ENDIF
        ENDIF    
      ENDDO
    ENDDO

    CALL MPI_Send_matrix(H_overlap,bound1_MPI,bound2_MPI,bound1_MPI,bound2_MPI,        &  
                         root_MPI,MPI_id)
    CALL MPI_Send_matrix(S_overlap,bound1_MPI,bound2_MPI,bound1_MPI,bound2_MPI,        &  
                         root_MPI,MPI_id)
                         write(*,*) 'send done',MPI_id

  ENDIF ! for MPI_id/=0  

!  DO i_MPI=1,MPI_np-1
!    bound1_MPI=i_MPI*nb_per_MPI+1+MIN(i_MPI,nb_rem_MPI)
!    bound2_MPI=(i_MPI+1)*nb_per_MPI+MIN(i_MPI,nb_rem_MPI)+merge(1,0,nb_rem_MPI>i_MPI)
!    CALL MPI_Bcast_matrix(H_overlap,bound1_MPI,bound2_MPI,bound1_MPI,bound2_MPI,i_MPI)
!    CALL MPI_Bcast_matrix(S_overlap,bound1_MPI,bound2_MPI,bound1_MPI,bound2_MPI,i_MPI)
!  ENDDO
  ! a bit waste of comm. time
  CALL MPI_Bcast_matrix(H_overlap,1,ndim,1,ndim,root_MPI)
  CALL MPI_Bcast_matrix(S_overlap,1,ndim,1,ndim,root_MPI)

END SUBROUTINE Overlap_psi1_psi2_MPI2
#endif
!=======================================================================================

!=======================================================================================
! subroutine for the calculation of Overlap_psi1_psi2 with MPI 
! in loop: 
! DO i=i_l,i_u
! DO j=j_l,j_u
!=======================================================================================
#if(run_MPI)
SUBROUTINE Overlap_psi1_psi2_MPI(H_overlap,S_overlap,psi,Hpsi,ndim,With_Grid)
  USE mod_system
  USE mod_psi_set_alloc
  USE mod_MPI
  USE mod_MPI_Aid
  IMPLICIT NONE

  TYPE(param_psi), intent(inout)          :: psi(:)  !< on non-root threads allocated
  TYPE(param_psi), intent(inout)          :: Hpsi(:) !< on non-root threads allocated
  Real(kind=Rkind),intent(inout)          :: H_overlap(:,:)
  Real(kind=Rkind),intent(inout)          :: S_overlap(:,:)
  Logical,optional,intent(in)             :: With_Grid
  Integer,intent(in)                      :: ndim

  complex(kind=Rkind)                     :: Overlap
  logical                                 :: With_Grid_loc
  Integer                                 :: i
  Integer                                 :: j
  Logical                                 :: send_once(MPI_np-1)

  With_Grid_loc=.FALSE.
  IF (present(With_Grid)) With_Grid_loc = With_Grid

  ! only on master: psi%CvecG,psi%RvecG,psi%CvecB,psi%RvecB
  nb_per_MPI=(ndim)/MPI_np
  nb_rem_MPI=mod(ndim,MPI_np) !remainder jobs
  
  IF(MPI_id==0) THEN
    !> send Hpsi and psi
    DO i_MPI=1,MPI_np-1
      bound1_MPI=i_MPI*nb_per_MPI+1+MIN(i_MPI,nb_rem_MPI)
      bound2_MPI=(i_MPI+1)*nb_per_MPI+MIN(i_MPI,nb_rem_MPI)                            &
                                     +merge(1,0,nb_rem_MPI>i_MPI)
      DO i=1,ndim
        IF(.NOT. With_Grid_loc) THEN
          IF(psi(i)%cplx) THEN
            IF(i>=bound1_MPI .AND. i<=bound2_MPI) CALL MPI_Send(Hpsi(i)%CvecB,         &
                         Hpsi(i)%nb_tot,MPI_complex8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
            CALL MPI_Send(psi(i)%CvecB,psi(i)%nb_tot,MPI_Complex8,i_MPI,               &
                          i_MPI,MPI_COMM_WORLD,MPI_err)
          ELSE
            IF(i>=bound1_MPI .AND. i<=bound2_MPI) CALL MPI_Send(Hpsi(i)%RvecB,         &
                            Hpsi(i)%nb_tot,MPI_Real8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
            CALL MPI_Send(psi(i)%RvecB,psi(i)%nb_tot,MPI_Real8,i_MPI,                  &
                          i_MPI,MPI_COMM_WORLD,MPI_err)
          ENDIF
        ELSE
          IF(psi(i)%cplx) THEN
            IF(i>=bound1_MPI .AND. i<=bound2_MPI) CALL MPI_Send(Hpsi(i)%CvecG,         &
                        Hpsi(i)%nb_qaie,MPI_complex8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
            CALL MPI_Send(psi(i)%CvecG,psi(i)%nb_qaie,MPI_Complex8,i_MPI,              &
                          i_MPI,MPI_COMM_WORLD,MPI_err)
          ELSE
            IF(i>=bound1_MPI .AND. i<=bound2_MPI) CALL MPI_Send(Hpsi(i)%RvecG,         &
                           Hpsi(i)%nb_qaie,MPI_Real8,i_MPI,i_MPI,MPI_COMM_WORLD,MPI_err)
            CALL MPI_Send(psi(i)%RvecG,psi(i)%nb_qaie,MPI_Real8,i_MPI,                 &
                          i_MPI,MPI_COMM_WORLD,MPI_err)
          ENDIF
        ENDIF ! for .NOT. With_Grid_loc  
      ENDDO ! for i=1,ndim
    ENDDO ! for i_MPI=1,MPI_np-1
    
    !> calcuation on master
    i_MPI=0
    DO i=1,ndim
      bound1_MPI=i_MPI*nb_per_MPI+1+MIN(i_MPI,nb_rem_MPI)
      bound2_MPI=(i_MPI+1)*nb_per_MPI+MIN(i_MPI,nb_rem_MPI)+merge(1,0,nb_rem_MPI>i_MPI)
      IF(i>=bound1_MPI .AND. i<=bound2_MPI) THEN
        DO j=1,ndim
          ! main claculation on master
          CALL Overlap_psipsi_MPI(Overlap,psi(j),Hpsi(i),With_Grid=With_Grid)
          H_overlap(j,i)=real(Overlap,kind=Rkind) 
          CALL Overlap_psipsi_MPI(Overlap,psi(j), psi(i),With_Grid=With_Grid)
          S_overlap(j,i)=real(Overlap,kind=Rkind)
        ENDDO ! for j=1,ndim
      ENDIF ! for i>=bound1_MPI .AND. i<=bound2_MPI
    END DO ! for i=1,ndim
    
  ENDIF ! for MPI_id==0
    
  !-------------------------------------------------------------------------------------
  bound1_MPI=MPI_id*nb_per_MPI+1+MIN(MPI_id,nb_rem_MPI)
  bound2_MPI=(MPI_id+1)*nb_per_MPI+MIN(MPI_id,nb_rem_MPI)+merge(1,0,nb_rem_MPI>MPI_id)
  IF(MPI_id/=0) THEN
    DO i=1,ndim
      IF (.NOT. With_Grid_loc) THEN
        IF (psi(i)%cplx) THEN
          IF(i>=bound1_MPI .AND. i<=bound2_MPI) THEN
            IF(.NOT. allocated(Hpsi(i)%CvecB)) CALL alloc_NParray(Hpsi(i)%CvecB,       &
                                            (/Hpsi(i)%nb_tot/),'Hpsi%CvecB','alloc_psi')
            CALL MPI_Recv(Hpsi(i)%CvecB,Hpsi(i)%nb_tot,MPI_complex8,root_MPI,MPI_id,   &
                          MPI_COMM_WORLD,MPI_stat,MPI_err)
          ENDIF          
          IF(.NOT. allocated(psi(i)%CvecB)) CALL alloc_NParray(psi(i)%CvecB,           &
                                              (/psi(i)%nb_tot/),'psi%CvecB','alloc_psi')
          CALL MPI_Recv( psi(i)%CvecB, psi(i)%nb_tot,MPI_Complex8,root_MPI,MPI_id,     &
                        MPI_COMM_WORLD,MPI_stat,MPI_err)
        ELSE
          IF(i>=bound1_MPI .AND. i<=bound2_MPI) THEN
            IF(.NOT. allocated(Hpsi(i)%RvecB)) CALL alloc_NParray(Hpsi(i)%RvecB,       &
                                            (/Hpsi(i)%nb_tot/),'Hpsi%RvecB','alloc_psi')
            CALL MPI_Recv(Hpsi(i)%RvecB,Hpsi(i)%nb_tot,MPI_Real8,root_MPI,MPI_id,      &
                          MPI_COMM_WORLD,MPI_stat,MPI_err)
          ENDIF
          IF(.NOT. allocated(psi(i)%RvecB))  CALL alloc_NParray(psi(i)%RvecB,          &
                                              (/psi(i)%nb_tot/),'psi%RvecB','alloc_psi')
          CALL MPI_Recv( psi(i)%RvecB, psi(i)%nb_tot,MPI_Real8,root_MPI,MPI_id,        &
                        MPI_COMM_WORLD,MPI_stat,MPI_err)
        ENDIF
      ELSE
        IF (psi(i)%cplx) THEN
          IF(i>=bound1_MPI .AND. i<=bound2_MPI) THEN
            IF(.NOT. allocated(Hpsi(i)%CvecG)) CALL alloc_NParray(Hpsi(i)%CvecG,       &
                                           (/Hpsi(i)%nb_qaie/),'Hpsi%CvecG','alloc_psi')
            CALL MPI_Recv(Hpsi(i)%CvecG,Hpsi(i)%nb_qaie,MPI_complex8,root_MPI,MPI_id,  &
                          MPI_COMM_WORLD,MPI_stat,MPI_err)
          ENDIF
          IF(.NOT. allocated(psi(i)%CvecG))  CALL alloc_NParray(psi(i)%CvecG,          &
                                             (/psi(i)%nb_qaie/),'psi%CvecG','alloc_psi')
          CALL MPI_Recv( psi(i)%CvecG, psi(i)%nb_qaie,MPI_Complex8,root_MPI,MPI_id,    &
                        MPI_COMM_WORLD,MPI_stat,MPI_err)
        ELSE
          IF(i>=bound1_MPI .AND. i<=bound2_MPI) THEN
            IF(.NOT. allocated(Hpsi(i)%RvecB)) CALL alloc_NParray(Hpsi(i)%RvecB,       &
                                            (/Hpsi(i)%nb_tot/),'Hpsi%RvecB','alloc_psi')
            CALL MPI_Recv(Hpsi(i)%RvecB,Hpsi(i)%nb_tot,MPI_Real8,root_MPI,MPI_id,      &
                          MPI_COMM_WORLD,MPI_stat,MPI_err)
          ENDIF
          IF(.NOT. allocated(psi(i)%RvecB))  CALL alloc_NParray(psi(i)%RvecB,          &
                                              (/psi(i)%nb_tot/),'psi%RvecB','alloc_psi')
          CALL MPI_Recv( psi(i)%RvecB, psi(i)%nb_tot,MPI_Real8,root_MPI,MPI_id,        &
                        MPI_COMM_WORLD,MPI_stat,MPI_err)
        ENDIF
      ENDIF ! for .NOT. With_Grid_loc
    ENDDO ! for i=1,ndim
    
    ! main calculation
    DO i=1,ndim
      IF(i>=bound1_MPI .AND. i<=bound2_MPI) THEN
        DO j=1,ndim
          IF (.NOT. With_Grid_loc) THEN
            IF (psi(j)%symab > -1 .AND. Hpsi(i)%symab > -1                             &
                                  .AND. Hpsi(i)%symab /= Hpsi(i)%symab) THEN
            ELSE
              CALL Overlap_psipsi_MPI(Overlap,psi(j),Hpsi(i),With_Grid=With_Grid)
              H_overlap(j,i)=real(Overlap,kind=Rkind)
              CALL Overlap_psipsi_MPI(Overlap,psi(j), psi(i),With_Grid=With_Grid)
              S_overlap(j,i)=real(Overlap,kind=Rkind)
            ENDIF
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDIF ! for MPI_id/=0  

  DO i_MPI=0,MPI_np-1
    bound1_MPI=i_MPI*nb_per_MPI+1+MIN(i_MPI,nb_rem_MPI)
    bound2_MPI=(i_MPI+1)*nb_per_MPI+MIN(i_MPI,nb_rem_MPI)+merge(1,0,nb_rem_MPI>i_MPI)
    CALL MPI_Bcast_matrix(H_overlap,1,ndim,bound1_MPI,bound2_MPI,i_MPI)
    CALL MPI_Bcast_matrix(S_overlap,1,ndim,bound1_MPI,bound2_MPI,i_MPI)
  ENDDO

END SUBROUTINE Overlap_psi1_psi2_MPI
#endif
!=======================================================================================

!=======================================================================================
!> subroutine calculating overlap of psi1 and psi2 with MPI
!>  NOTE: work with Overlap_psi1_psi2_MPI3
!>  be careful with the way distributing array, 
!>  which is ready in Overlap_psi1_psi2_MPI3
!=======================================================================================
#if(run_MPI)
SUBROUTINE Overlap_psipsi_MPI3(Overlap,psi1,psi2,With_Grid,Channel_ie)
  USE mod_system
  USE mod_psi_set_alloc
  USE mod_MPI
  IMPLICIT NONE

  !-variables for the WP----------------------------------------------------------------
  TYPE(param_psi),intent(in)                :: psi1
  TYPE(param_psi),intent(in)                :: psi2
  Complex(kind=Rkind)                       :: Overlap
  Logical,optional,intent(in)               :: With_Grid
  Integer,optional,intent(in)               :: Channel_ie

  !-working variables-------------------------------------------------------------------
  Logical                                   :: With_Grid_loc
  Complex(kind=Rkind)                       :: temp
  Real(kind=Rkind)                          :: WrhonD
  Real(kind=Rkind)                          :: Roverlap
  Real(kind=Rkind)                          :: Rtemp
  Real(kind=Rkind),allocatable              :: wrho(:)
  Integer                                   :: locChannel_ie
  Integer                                   :: i_qa
  Integer                                   :: i_qaie
  Integer                                   :: i_be
  Integer                                   :: i_bi
  Integer                                   :: i_ba
  Integer                                   :: i_baie
  Integer                                   :: f_baie
  Integer                                   :: iie
  Integer                                   :: fie
  Integer                                   :: iie_MPI
  Integer                                   :: fie_MPI
  
  !-for debuging------------------------------------------------------------------------
  character (len=*), parameter :: name_sub='Overlap_psipsi_MPI3'
  logical,parameter :: debug = .FALSE.
  ! logical,parameter :: debug = .TRUE.

  !-------------------------------------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*) 'psi1'
    CALL ecri_psi(psi=psi1)

    write(out_unitp,*) 'psi2'
    CALL ecri_psi(psi=psi2)
    write(out_unitp,*) 'GridRep,BasisRep ?'
    IF (present(With_Grid)) write(out_unitp,*) 'With_Grid',With_Grid
    IF (present(Channel_ie)) write(out_unitp,*) 'Channel_ie',Channel_ie
  ENDIF
  !-------------------------------------------------------------------------------------

  With_Grid_loc = .FALSE.

  IF(present(With_Grid)) With_Grid_loc=With_Grid

  locChannel_ie = 0
  IF(present(Channel_ie)) locChannel_ie=Channel_ie

  IF (psi1%nb_baie>psi1%nb_tot) THEN
    With_Grid_loc = .FALSE.
  ENDIF

  ! get bound1_MPI and bound2_MPI
  IF(.NOT. With_Grid_loc) THEN
    nb_per_MPI=psi1%nb_tot/MPI_np
    nb_rem_MPI=mod(psi1%nb_tot,MPI_np) 
  ELSE
    nb_per_MPI=psi1%nb_qaie/MPI_np
    nb_rem_MPI=mod(psi1%nb_qaie,MPI_np) 
  ENDIF
  bound1_MPI=MPI_id*nb_per_MPI+1+MIN(MPI_id,nb_rem_MPI)
  bound2_MPI=(MPI_id+1)*nb_per_MPI+MIN(MPI_id,nb_rem_MPI)+merge(1,0,nb_rem_MPI>MPI_id)

  ! With_Grid_loc: F
  IF(With_Grid_loc) THEN
    IF(psi1%cplx .AND. allocated(psi1%CvecG) .AND. allocated(psi2%CvecG)) THEN
    ELSE IF(.NOT. psi1%cplx .AND. allocated(psi1%RvecG) .AND.                          &
             allocated(psi2%RvecG)) THEN
    ELSE
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' impossible to calculate the GridRep overlap'
      write(out_unitp,*) ' With_Grid_loc=t but problem with the allocation GridRep'
      write(out_unitp,*) 'allocated(psi1%CvecG)',allocated(psi1%CvecG)
      write(out_unitp,*) 'allocated(psi2%CvecG)',allocated(psi2%CvecG)
      write(out_unitp,*) 'allocated(psi1%RvecG)',allocated(psi1%RvecG)
      write(out_unitp,*) 'allocated(psi2%RvecG)',allocated(psi2%RvecG)
      write(out_unitp,*) ' psi1'
      CALL ecri_psi(psi=psi1,ecri_GridRep=.TRUE.)
      write(out_unitp,*) ' psi2'
      CALL ecri_psi(psi=psi2,ecri_GridRep=.TRUE.)
      STOP
    ENDIF
  ELSE
    IF(psi1%cplx .AND.allocated(psi1%CvecB) .AND. allocated(psi2%CvecB)) THEN
    ELSE IF(.NOT. psi1%cplx .AND. allocated(psi1%RvecB) .AND.                          &
            allocated(psi2%RvecB)) THEN
    ELSE
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' impossible to calculate the BasisRep overlap'
      write(out_unitp,*) ' With_Grid_loc=f (on basis) but problem with the allocation of BasisRep'
      write(out_unitp,*) 'allocated(psi1%CvecB)',allocated(psi1%CvecB)
      write(out_unitp,*) 'allocated(psi2%CvecB)',allocated(psi2%CvecB)
      write(out_unitp,*) 'allocated(psi1%RvecB)',allocated(psi1%RvecB)
      write(out_unitp,*) 'allocated(psi2%RvecB)',allocated(psi2%RvecB)
      write(out_unitp,*) ' psi1'
      CALL ecri_psi(psi=psi1,ecri_BasisRep=.TRUE.)
      write(out_unitp,*) ' psi2'
      CALL ecri_psi(psi=psi2,ecri_BasisRep=.TRUE.)
      STOP
    ENDIF
  ENDIF

  Overlap = cmplx(ZERO,ZERO,kind=Rkind)
  IF (.NOT. With_Grid_loc) THEN
    i_baie=1
    f_baie=psi1%nb_tot
    IF (psi1%nb_tot == psi1%nb_baie .AND.  locChannel_ie > 0 .AND.  &
                            locChannel_ie <= psi1%ComOp%nb_bie) THEN
      i_baie = 1 + (locChannel_ie-1)*psi1%nb_ba
      f_baie = i_baie-1 + psi1%nb_ba
    END IF

    IF(bound2_MPI>i_baie .AND. bound1_MPI<f_baie) THEN
      f_baie=MIN(f_baie,bound2_MPI)
      i_baie=MAX(i_baie,bound1_MPI)
      
      IF(psi1%symab>-1 .AND. psi2%symab>-1 .AND. psi1%symab/=psi2%symab) THEN
        !Overlap = cmplx(ZERO,ZERO,kind=Rkind)
      ELSE
        IF(psi1%cplx) THEN
          Overlap=dot_product(psi1%CvecB(i_baie:f_baie),psi2%CvecB(i_baie:f_baie))
        ELSE
          ROverlap=dot_product(psi1%RvecB(i_baie:f_baie),psi2%RvecB(i_baie:f_baie))
          Overlap=cmplx(ROverlap,ZERO,kind=Rkind)
        ENDIF
      ENDIF
    ENDIF ! bound2_MPI>i_baie .AND. bound1_MPI<f_baie

  ELSE ! With_Grid_loc
  
    CALL alloc_NParray(wrho,(/psi1%nb_qa/),"wrho",name_sub)
    DO i_qa=1,psi1%nb_qa
      wrho(i_qa) = Rec_WrhonD(psi1%BasisnD,i_qa)
    ENDDO

    IF(psi1%cplx) THEN
      iie=1
      fie=psi1%nb_qa
      DO i_be=1,psi1%nb_be
      DO i_bi=1,psi1%nb_bi
        IF(bound2_MPI>iie .AND. bound1_MPI<fie) THEN
          fie_MPI=MIN(fie,bound2_MPI)
          iie_MPI=MAX(iie,bound1_MPI)
          Overlap=Overlap+dot_product(psi1%CvecG(iie_MPI:fie_MPI),                     &
                          wrho(Mod(iie_MPI,psi1%nb_qa):Mod(fie_MPI,psi1%nb_qa))        &
                          *psi2%CvecG(iie_MPI:fie_MPI))
        ENDIF
        iie=iie+psi1%nb_qa
        fie=fie+psi1%nb_qa
      ENDDO
      ENDDO
    ELSE
      iie=1
      fie=psi1%nb_qa
      DO i_be=1,psi1%nb_be
      DO i_bi=1,psi1%nb_bi
        IF(bound2_MPI>iie .AND. bound1_MPI<fie) THEN
          fie_MPI=MIN(fie,bound2_MPI)
          iie_MPI=MAX(iie,bound1_MPI)
          ROverlap=ROverlap+dot_product(psi1%RvecG(iie_MPI:fie_MPI),                   &
                            wrho(Mod(iie_MPI,psi1%nb_qa):Mod(fie_MPI,psi1%nb_qa))      &
                            *psi2%RvecG(iie_MPI:fie_MPI))
        ENDIF
        iie=iie + psi1%nb_qa
        fie=fie + psi1%nb_qa
      ENDDO
      ENDDO
      Overlap=cmplx(ROverlap,ZERO,kind=Rkind)
    ENDIF

    CALL dealloc_NParray(wrho,"wrho",name_sub)

  ENDIF

  !-------------------------------------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'Overlap : ',Overlap
    write(out_unitp,*) 'END ',name_sub
  ENDIF
  !-------------------------------------------------------------------------------------
  
END SUBROUTINE Overlap_psipsi_MPI3
#endif
!=======================================================================================

!=======================================================================================
! this is a temp subroutine for running Overlap_psi1_psi2 on all threads with MPI
! it will be replaced by the original 'Overlap_psi1_psi2' 
! and thus completely removed later
!=======================================================================================
#if(run_MPI)
      SUBROUTINE Overlap_psipsi_MPI(Overlap,psi1,psi2,With_Grid,Channel_ie)
      USE mod_system
      USE mod_psi_set_alloc
      USE mod_MPI
      IMPLICIT NONE

!----- variables for the WP ----------------------------------------
      TYPE (param_psi), intent(in)    :: psi1,psi2
      complex (kind=Rkind)            :: Overlap
      logical, optional, intent(in)   :: With_Grid
      integer, optional, intent(in)   :: Channel_ie

!------ working variables ---------------------------------
      logical              :: With_Grid_loc
      integer              :: locChannel_ie
      integer              :: i_qa,i_qaie
      integer              :: i_be,i_bi,i_ba
      integer              :: i_baie,f_baie
      integer              :: i_modif_q
      real (kind=Rkind)    :: WrhonD
      complex (kind=Rkind) :: temp
      real (kind=Rkind)    :: Roverlap,Rtemp
      integer              :: iie,fie
      real (kind=Rkind), allocatable :: wrho(:)

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Overlap_psipsi_MPI'
      logical,parameter :: debug = .FALSE.
!     logical,parameter :: debug = .TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'psi1'
        CALL ecri_psi(psi=psi1)

        write(out_unitp,*) 'psi2'
        CALL ecri_psi(psi=psi2)
        write(out_unitp,*) 'GridRep,BasisRep ?'
        IF (present(With_Grid)) write(out_unitp,*) 'With_Grid',With_Grid
        IF (present(Channel_ie)) write(out_unitp,*) 'Channel_ie',Channel_ie
      END IF
!-----------------------------------------------------------

      With_Grid_loc = .FALSE.

      IF (present(With_Grid)) With_Grid_loc = With_Grid

      locChannel_ie = 0
      IF (present(Channel_ie)) locChannel_ie = Channel_ie

      IF (psi1%nb_baie > psi1%nb_tot) THEN
        With_Grid_loc = .FALSE.
      END IF

      ! With_Grid_loc: F
      IF (With_Grid_loc) THEN
        IF (psi1%cplx .AND.                                             &
         allocated(psi1%CvecG) .AND. allocated(psi2%CvecG)) THEN
        ELSE IF (.NOT. psi1%cplx .AND.                                  &
         allocated(psi1%RvecG) .AND. allocated(psi2%RvecG)) THEN
        ELSE
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' impossible to calculate the GridRep overlap'
          write(out_unitp,*) ' With_Grid_loc=t but problem with the allocation GridRep'
          write(out_unitp,*) 'allocated(psi1%CvecG)',allocated(psi1%CvecG)
          write(out_unitp,*) 'allocated(psi2%CvecG)',allocated(psi2%CvecG)
          write(out_unitp,*) 'allocated(psi1%RvecG)',allocated(psi1%RvecG)
          write(out_unitp,*) 'allocated(psi2%RvecG)',allocated(psi2%RvecG)
          write(out_unitp,*) ' psi1'
          CALL ecri_psi(psi=psi1,ecri_GridRep=.TRUE.)
          write(out_unitp,*) ' psi2'
          CALL ecri_psi(psi=psi2,ecri_GridRep=.TRUE.)
          STOP
        END IF
      ELSE
        IF (psi1%cplx .AND.                                             &
         allocated(psi1%CvecB) .AND. allocated(psi2%CvecB)) THEN
        ELSE IF (.NOT. psi1%cplx .AND.                                  &
         allocated(psi1%RvecB) .AND. allocated(psi2%RvecB)) THEN
        ELSE
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' impossible to calculate the BasisRep overlap'
          write(out_unitp,*) ' With_Grid_loc=f (on basis) but problem with the allocation of BasisRep'
          write(out_unitp,*) 'allocated(psi1%CvecB)',allocated(psi1%CvecB)
          write(out_unitp,*) 'allocated(psi2%CvecB)',allocated(psi2%CvecB)
          write(out_unitp,*) 'allocated(psi1%RvecB)',allocated(psi1%RvecB)
          write(out_unitp,*) 'allocated(psi2%RvecB)',allocated(psi2%RvecB)
          write(out_unitp,*) ' psi1'
          CALL ecri_psi(psi=psi1,ecri_BasisRep=.TRUE.)
          write(out_unitp,*) ' psi2'
          CALL ecri_psi(psi=psi2,ecri_BasisRep=.TRUE.)
          STOP
        END IF
      END IF

      IF (.NOT. With_Grid_loc) THEN
        i_baie=1
        f_baie=psi1%nb_tot
        IF (psi1%nb_tot == psi1%nb_baie .AND.  locChannel_ie > 0 .AND.  &
                                locChannel_ie <= psi1%ComOp%nb_bie) THEN
          i_baie = 1 + (locChannel_ie-1)*psi1%nb_ba
          f_baie = i_baie-1 + psi1%nb_ba
        END IF
        IF (psi1%symab > -1 .AND. psi2%symab > -1 .AND. psi1%symab /= psi2%symab) THEN
          Overlap = cmplx(ZERO,ZERO,kind=Rkind)
        ELSE
          IF (psi1%cplx) THEN
            Overlap = dot_product( psi1%CvecB(i_baie:f_baie) ,          &
                                   psi2%CvecB(i_baie:f_baie) )
          ELSE
            ROverlap = dot_product( psi1%RvecB(i_baie:f_baie) ,         &
                                    psi2%RvecB(i_baie:f_baie) )
            Overlap = cmplx(ROverlap,ZERO,kind=Rkind)
          END IF
        END IF

      ELSE

!       - initialization ----------------------------------
        Overlap = cmplx(ZERO,ZERO,kind=Rkind)

        CALL alloc_NParray(wrho,(/ psi1%nb_qa/),"wrho",name_sub)
        DO i_qa=1,psi1%nb_qa
          wrho(i_qa) = Rec_WrhonD(psi1%BasisnD,i_qa)
        END DO

        IF (psi1%cplx) THEN
          iie = 1
          fie = psi1%nb_qa
          DO i_be=1,psi1%nb_be
          DO i_bi=1,psi1%nb_bi
            Overlap = Overlap + dot_product(                            &
              psi1%CvecG(iie:fie),wrho*psi2%CvecG(iie:fie))
            iie = iie + psi1%nb_qa
            fie = fie + psi1%nb_qa
          END DO
          END DO
        ELSE
          iie = 1
          fie = psi1%nb_qa
          DO i_be=1,psi1%nb_be
          DO i_bi=1,psi1%nb_bi
            Overlap = Overlap + cmplx(dot_product(                      &
              psi1%RvecG(iie:fie),wrho*psi2%RvecG(iie:fie)) ,kind=Rkind)
            iie = iie + psi1%nb_qa
            fie = fie + psi1%nb_qa
          END DO
          END DO
        END IF

        CALL dealloc_NParray(wrho,"wrho",name_sub)

      END IF

!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Overlap : ',Overlap
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------

      END SUBROUTINE Overlap_psipsi_MPI
#endif
!=======================================================================================

!============================================================
!
!   trie des vecteur dans l'ordre croissant
!   le vecteur i est psi(i) with ene(i)
!
!============================================================
!
      SUBROUTINE trie_psi(psi,ene,nb_wp)
      USE mod_system
      USE mod_psi_set_alloc
      IMPLICIT NONE

!----- variables for the WP ----------------------------------------
      integer :: nb_wp
      TYPE (param_psi), intent(inout) :: psi(nb_wp)
      TYPE (param_psi)                :: psi_temp

      real (kind=Rkind) :: ene(nb_wp)
      real (kind=Rkind) :: a


      integer       :: i,j,k


      DO i=1,nb_wp
        DO j=i+1,nb_wp
          IF (ene(i) > ene(j)) THEN
            !permutation
            a=ene(i)
            ene(i)=ene(j)
            ene(j)=a

            psi_temp = psi(i)
            psi(i) = psi(j)
            psi(j) = psi_temp
          END IF
        END DO
      END DO
      CALL dealloc_psi(psi_temp)

      END SUBROUTINE trie_psi


!================================================================
!
!     Save vectors
!
!================================================================
      SUBROUTINE sub_LCpsi_TO_psi(psi,Vec,ndim,nb_save)
      USE mod_system
      USE mod_psi_set_alloc
      USE mod_psi_SimpleOp
      IMPLICIT NONE


!----- variables for the WP propagation ----------------------------
      integer            :: ndim,nb_save
      TYPE (param_psi)   :: psi(ndim)
      real (kind=Rkind)  :: Vec(ndim,ndim)



!------ working parameters --------------------------------
      integer            :: i,k,isym
      real (kind=Rkind), allocatable  :: PsiRk(:)
      integer            :: symab_psi_old(ndim)



!----- for debuging --------------------------------------------------
      character (len=*), parameter ::name_sub='sub_LCpsi_TO_psi'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) ' nb_save,ndim',nb_save,ndim
        CALL Write_Mat(Vec,out_unitp,5)
        write(out_unitp,*)
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

      IF (debug) write(out_unitp,*) 'sum(abs(Vec))',sum(abs(Vec))
      IF (sum(abs(Vec)) > ONETENTH**10) THEN

        IF (psi(1)%BasisRep) THEN

          DO i=1,ndim
            symab_psi_old(i) = psi(i)%symab
          END DO

          !$OMP parallel default(none) &
          !$OMP shared(ndim,psi,Vec) &
          !$OMP private(i,k,PsiRk)

          CALL alloc_NParray(PsiRk,(/ndim/),'PsiRk',name_sub)

          !$OMP do
          DO k=1,size(psi(1)%RvecB)
            DO i=1,ndim
              PsiRk(i) = psi(i)%RvecB(k)
            END DO


            PsiRk(:) = matmul(PsiRk,Vec)

            DO i=1,ndim
              psi(i)%RvecB(k) = PsiRk(i)
            END DO
          END DO
          !$OMP end do
          CALL dealloc_NParray(PsiRk,'PsiRk',name_sub)
          !$OMP end parallel

          DO i=1,ndim
            isym = maxloc(abs(Vec(:,i)),dim=1)
            CALL Set_symab_OF_psiBasisRep(psi(i),symab_psi_old(isym))
          END DO
        ELSE

          !$OMP parallel default(none) &
          !$OMP shared(ndim,psi,Vec) &
          !$OMP private(i,k,PsiRk)

          CALL alloc_NParray(PsiRk,(/ndim/),'PsiRk',name_sub)

          !$OMP do
          DO k=1,size(psi(1)%RvecG)
            DO i=1,ndim
              PsiRk(i) = psi(i)%RvecG(k)
            END DO

            PsiRk = matmul(PsiRk,Vec)

            DO i=1,ndim
              psi(i)%RvecG(k) = PsiRk(i)
            END DO
          END DO
          !$OMP end do
          CALL dealloc_NParray(PsiRk,'PsiRk',name_sub)
          !$OMP end parallel
        END IF


      ELSE
        CONTINUE ! nothing!!!
      END IF

!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
!----------------------------------------------------------

      END SUBROUTINE sub_LCpsi_TO_psi
!================================================================
!
!     Schmidt ortho
!
!================================================================
      SUBROUTINE sub_Schmidt(psi,nb_psi)
      USE mod_system
      USE mod_psi_set_alloc
      USE mod_psi_SimpleOp
      USE mod_ana_psi
      IMPLICIT NONE

      integer            :: nb_psi
      TYPE (param_psi)   :: psi(nb_psi)


!------ working parameters --------------------------------
      complex (kind=Rkind) :: Overlap
      real    (kind=Rkind) :: ROverlap
      integer       :: i,j,sym

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_Schmidt'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_psi',nb_psi
       END IF
!-----------------------------------------------------------

      DO i=1,nb_psi
        sym = psi(i)%symab
        DO j=1,i-1
          CALL Overlap_psi1_psi2(Overlap,psi(i),psi(j))
          IF (abs(Overlap) == ZERO) CYCLE
          IF (psi(j)%cplx) THEN
            psi(i) = psi(i) - psi(j) * Overlap
          ELSE
            ROverlap = real(Overlap,kind=Rkind)
            psi(i) = psi(i) - psi(j) * ROverlap
          END IF
!         write(out_unitp,*) 'j,i,S',j,i,Overlap
!         CALL flush_perso(out_unitp)
        END DO
        CALL Set_symab_OF_psiBasisRep(psi(i),sym)

!       CALL norm2_psi(psi(i))
!       write(out_unitp,*) ' Ortho: norme',i,psi(i)%norme
        CALL renorm_psi(psi(i))
        !write(out_unitp,*) 'symab, bits(symab)',WriteTOstring_symab(psi(i)%symab)

      END DO

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
       END IF
!-----------------------------------------------------------
      END SUBROUTINE sub_Schmidt

      SUBROUTINE sub_Lowdin(psi,nb_psi)
      USE mod_system
      USE mod_psi_set_alloc
      USE mod_psi_SimpleOp
      USE mod_ana_psi
      IMPLICIT NONE

      integer            :: nb_psi
      TYPE (param_psi)   :: psi(nb_psi)


!------ working parameters --------------------------------
      TYPE (param_psi)   :: TempPsi(nb_psi)

      real    (kind=Rkind) :: RS(nb_psi,nb_psi)
      real    (kind=Rkind) :: Vec(nb_psi,nb_psi)
      real    (kind=Rkind) :: Eig(nb_psi)

      complex (kind=Rkind) :: CS(nb_psi,nb_psi)

      complex (kind=Rkind) :: Overlap
      real    (kind=Rkind) :: ROverlap
      integer       :: i,j

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_Lowdin'
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_psi',nb_psi
       END IF
!-----------------------------------------------------------

      IF (psi(1)%cplx) THEN
        DO i=1,nb_psi
        CALL renorm_psi(psi(i))
        DO j=1,i-1
          CALL Overlap_psi1_psi2(Overlap,psi(i),psi(j))
          CS(i,j) =  Overlap
          CS(j,i) =  Overlap
        END DO
        END DO
        STOP 'complex not yet'
      ELSE
        ! first the overlap matrix
        DO i=1,nb_psi
        CALL renorm_psi(psi(i))
        DO j=1,i
          CALL Overlap_psi1_psi2(Overlap,psi(i),psi(j))
          RS(i,j) =  Real(Overlap,kind=Rkind)
          RS(j,i) =  Real(Overlap,kind=Rkind)
        END DO
        END DO
        IF (debug) CALL Write_Mat(RS,out_unitp,5)

        CALL diagonalization(RS,Eig,Vec,nb_psi,1,-1,.FALSE.)
        IF (debug) THEN
          write(out_unitp,*) 'Eig S ',Eig(:)
          write(out_unitp,*) 'nb large vp ',count(Eig>ONETENTH**6)
        END IF

        DO i=1,nb_psi
          TempPsi(i) = psi(i)
        END DO

        DO i=1,nb_psi
          psi(i) = ZERO
          DO j=1,nb_psi
            psi(i) = psi(i) + Vec(j,i) * TempPsi(j)
          END DO
          CALL renorm_psi(psi(i))
          psi(i)%CAvOp    = cmplx(Eig(i),kind=Rkind)
        END DO

        DO i=1,nb_psi
          CALL dealloc_psi(TempPsi(i),delete_all=.TRUE.)
        END DO

      END IF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Eig S ',Eig(:)
        write(out_unitp,*) 'END ',name_sub
       END IF
!-----------------------------------------------------------
      END SUBROUTINE sub_Lowdin


      END MODULE mod_psi_Op

