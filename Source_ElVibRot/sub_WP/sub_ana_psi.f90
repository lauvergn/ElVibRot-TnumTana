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

!==============================================================
!
!      calc_DM
!
!==============================================================
      SUBROUTINE calc_DM(psi,max_1D,T,info,print_DM)
      USE mod_system
      USE mod_psi_set_alloc
      IMPLICIT NONE


!----- variables for the WP propagation ----------------------------
      TYPE (param_psi)   :: psi


      complex (kind=Rkind),allocatable :: DM(:,:),DM2(:,:)
      complex (kind=Rkind),allocatable :: DM_v1(:,:)
      complex (kind=Rkind),allocatable :: DM2_v1(:,:)
      real (kind=Rkind),allocatable    :: weight(:)
      complex (kind=Rkind) :: tr,tr2

      real (kind=Rkind)    :: T ! time
      logical          :: print_DM
      character (len=*)    :: info

      integer          :: nb_ba
      integer          :: ib,ibie
      integer          :: max_1D,max_print

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='calc_DM'
      logical, parameter :: debug =.FALSE.
!     logical, parameter :: debug =.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_tot',psi%nb_tot
        write(out_unitp,*) 'nb_bi,nb_be',psi%nb_bi,psi%nb_be
        write(out_unitp,*) 'nb_bie',psi%nb_bi*psi%nb_be
      END IF
!-----------------------------------------------------------

      IF (print_DM) THEN
        max_print = min(50,psi%nb_tot)
        CALL alloc_NParray(weight,(/psi%nb_tot/),"weight",name_sub)
        weight(:) = real(conjg(psi%CvecB)*psi%CvecB,kind=Rkind)
!       write(out_unitp,*) 'T weight_tot',T,weight(1:max_print)
        write(out_unitp,'(a,f12.1,1x,a,50(1x,f8.6))')                           &
                 'T weight_tot ',T,trim(info),weight(1:max_print)
        CALL dealloc_NParray(weight,"weight",name_sub)
      END IF
      IF (psi%nb_baie /= psi%nb_tot) THEN
        CALL alloc_NParray(DM ,(/psi%nb_tot,psi%nb_tot/),"DM" ,name_sub)
        CALL alloc_NParray(DM2,(/psi%nb_tot,psi%nb_tot/),"DM2",name_sub)

        DM(:,:) = CZERO
        DM(:,:) = matmul(                                             &
                conjg(reshape(psi%CvecB,(/psi%nb_tot,1/)) ),          &
                      reshape(psi%CvecB,(/1,psi%nb_tot/))  )

        DM2 = matmul(DM,DM)

        tr  = CZERO
        tr2 = CZERO
        DO  ibie=1,psi%nb_tot
          tr = tr + DM(ibie,ibie)
          tr2 = tr2 + DM2(ibie,ibie)
        END DO
        write(out_unitp,*) 'T tr(DM),tr(DM2)',T,abs(tr),abs(tr2)
        CALL dealloc_NParray(DM ,"DM" ,name_sub)
        CALL dealloc_NParray(DM2,"DM2",name_sub)
      END IF

      IF (psi%nb_baie == psi%nb_tot .AND. psi%nb_bi*psi%nb_be > 1) THEN
        nb_ba = psi%nb_ba
        CALL alloc_NParray(DM_v1 ,(/nb_ba,nb_ba/),"DM_v1",name_sub)
        CALL alloc_NParray(DM2_v1,(/nb_ba,nb_ba/),"DM2_v1",name_sub)

        DM_v1(:,:) = matmul(                                            &
                conjg(reshape(psi%CvecB,(/nb_ba,1/)) ),                 &
                      reshape(psi%CvecB,(/1,nb_ba/))  )
!       DM_v1(:,:) = DM(1:psi%nb_ba,1:psi%nb_ba)
        DM2_v1 = matmul(DM_v1,DM_v1)
        tr  = CZERO
        tr2 = CZERO
        DO  ib=1,psi%nb_ba
          tr = tr + DM_v1(ib,ib)
          tr2 = tr2 + DM2_v1(ib,ib)
        END DO
        write(out_unitp,*) 'T v1 tr(DM),tr(DM2)',T,abs(tr),abs(tr2)

        CALL dealloc_NParray(DM_v1 ,"DM_v1" ,name_sub)
        CALL dealloc_NParray(DM2_v1,"DM2_v1",name_sub)
      END IF

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------



      END SUBROUTINE calc_DM
!==============================================================
!
!      calc_1Dweight(psi)
!
!==============================================================
      SUBROUTINE calc_1Dweight(psi,ana_psi,max_1D,T,info,print_w)
      USE mod_system
      USE mod_dnSVM
      USE mod_psi_set_alloc
      USE mod_psi_Op
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_psi)     :: psi
      TYPE (param_ana_psi) :: ana_psi

      integer              :: max_1D
      real (kind=Rkind)    :: T ! time
      character (len=*)    :: info
      logical              :: print_w


!----- for debuging --------------------------------------------------
      integer :: err_mem
      character (len=*), parameter :: name_sub='calc_1Dweight'
      logical, parameter :: debug =.FALSE.
!      logical, parameter :: debug =.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL RecWrite_basis(psi%BasisnD,write_all=.TRUE.)
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

      CALL Set_symab_OF_psiBasisRep(psi)

      IF (psi%nb_baie*psi%nb_bRot == psi%nb_tot .AND. ana_psi%ana) THEN
        CALL calc_1Dweight_act1(psi,ana_psi,max_1D,T,info,print_w)
      END IF

      CALL calc_1Dweight_inact2n_elec(psi,ana_psi,max_1D,T,info,print_w)

      CALL calc_MaxCoef_psi(psi,T,info)

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------


      END SUBROUTINE calc_1Dweight

      SUBROUTINE calc_1Dweight_inact2n_elec(psi,ana_psi,max_1D,T,info,print_w)
      USE mod_system
      USE mod_psi_set_alloc
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_psi)     :: psi
      TYPE (param_ana_psi) :: ana_psi

      integer          :: max_herm
      real (kind=Rkind),allocatable :: weight1D(:,:)
      integer          :: tab(psi%Basis2n%nb_basis+1)

      real (kind=Rkind),allocatable :: weight1Dact(:,:)
      real (kind=Rkind)    :: a

      real (kind=Rkind)    :: T ! time
      character (len=*)    :: info
      logical          :: print_w

      integer              :: i,j,j_herm,beg_e
      integer              :: ie,ii,ib,ibie,iq,ibiq
      integer              :: max_dim,max_1D
      integer, allocatable :: nDval(:)

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='calc_1Dweight_inact2n_elec'
      integer :: err_mem,memory
      logical, parameter :: debug =.FALSE.
!      logical, parameter :: debug =.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_inact2n',psi%Basis2n%nb_basis
        write(out_unitp,*) 'nb_bi,nb_be',psi%nb_bi,psi%nb_be
        write(out_unitp,*) 'nb_bie',psi%nb_bi*psi%nb_be
        write(out_unitp,*) 'tab_WeightChannels',ana_psi%tab_WeightChannels
      END IF
!-----------------------------------------------------------
      CALL flush_perso(out_unitp)

      IF (psi%Basis2n%nb_basis == 0) THEN
        max_herm = psi%nb_be-1
      ELSE
        CALL alloc_NParray(nDval,(/psi%Basis2n%nb_basis/),'nDval',name_sub)
        max_herm = 0
        DO i=1,psi%Basis2n%nb_basis
          IF (psi%Basis2n%tab_Pbasis(i)%Pbasis%nb > max_herm)           &
                      max_herm = psi%Basis2n%tab_Pbasis(i)%Pbasis%nb
        END DO
        max_herm = max_herm - 1
      END IF

      CALL alloc_NParray(weight1D,(/psi%Basis2n%nb_basis+1,max_herm/),    &
                        "weight1D",name_sub,(/1,0/))

      IF (psi%nb_bi > 1 .OR. psi%nb_be > 1) THEN
        weight1D(:,:) = ZERO
        i = 0
        DO ie=1,psi%nb_be
        DO ii=1,psi%nb_bi
          i = i + 1
          IF (psi%Basis2n%nb_basis > 0) THEN
            CALL calc_nDindex(psi%Basis2n%nDindB,i,nDval)
            tab(1:psi%Basis2n%nb_basis) = nDval(:)-1
            tab(psi%Basis2n%nb_basis+1) = ie
          ELSE
            tab(psi%Basis2n%nb_basis+1) = ie-1
          END IF


          DO j=1,psi%Basis2n%nb_basis+1
            j_herm = tab(j)
            weight1D(j,j_herm) = weight1D(j,j_herm) +                   &
                                        ana_psi%tab_WeightChannels(ii,ie)
          END DO
        END DO
        END DO

        IF (print_w .OR. debug) THEN
          DO i=1,psi%Basis2n%nb_basis
            write(out_unitp,11) 'harm T W ',trim(info),i,T,             &
                    weight1D(i,0:psi%Basis2n%tab_Pbasis(i)%Pbasis%nb-1)
          END DO
          beg_e = 1
          IF (psi%Basis2n%nb_basis == 0) beg_e = 0
          i = psi%Basis2n%nb_basis+1
          write(out_unitp,11) 'elec T W ',trim(info),i,T,               &
                weight1D(i,beg_e:beg_e+psi%nb_be-1)

 11       format(a,a,1x,i3,1x,f15.2,20(1x,f6.3))
        END IF
      END IF

      CALL dealloc_NParray(weight1D,"weight1D",name_sub)
      IF (allocated(nDval)) CALL dealloc_NParray(nDval,'nDval',name_sub)

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------



      END SUBROUTINE calc_1Dweight_inact2n_elec

      SUBROUTINE calc_MaxCoef_psi(psi,T,info)
      USE mod_system
      USE mod_psi_set_alloc
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_psi)     :: psi

      real (kind=Rkind)    :: T ! time
      character (len=*)    :: info

      real (kind=Rkind)  :: C,maxC1,maxC2
      integer            :: i_bhe,i_b,i_e,i_h,i_R
      integer            :: i_bhe1,i_bhe2,i_b_maxC1,i_b_maxC2
      integer            :: i_R_maxC1,i_R_maxC2
      integer            :: i_e_maxC1,i_e_maxC2
      integer            :: i_h_maxC1,i_h_maxC2
      integer            :: ndim_index(psi%BasisnD%ndim)


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='calc_MaxCoef_psi'
      logical, parameter :: debug =.FALSE.
!      logical, parameter :: debug =.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'asso psi%BasisnD',associated(psi%BasisnD)
      END IF
!-----------------------------------------------------------
      IF (psi%nb_baie*psi%nb_bRot /= psi%nb_tot .AND.                   &
          .NOT. psi%ComOp%contrac_ba_ON_HAC) RETURN   ! should be spectral WP

      maxC1 = ZERO
      maxC2 = ZERO
      i_b_maxC1 = 0
      i_b_maxC2 = 0
      i_h_maxC1 = 0
      i_h_maxC2 = 0
      i_e_maxC1 = 0
      i_e_maxC2 = 0
      i_R_maxC1 = 0
      i_R_maxC2 = 0
      i_bhe1    = 0
      i_bhe2    = 0

      i_bhe = 0
      DO i_R=1,psi%nb_bRot
      DO i_e=1,psi%nb_be
      DO i_h=1,psi%nb_bi
      !DO i_b=1,psi%ComOp%nb_ba_ON_HAC(i_h)
      DO i_b=1,psi%nb_ba

        i_bhe = i_bhe + 1

        IF (psi%cplx) THEN
          C = abs(psi%CvecB(i_bhe))
        ELSE
          C = abs(psi%RvecB(i_bhe))
        END IF

        IF ( C > maxC1) THEN
          ! test if maxC1 > maxC2 (the old maxC1)
          IF ( maxC1 > maxC2) THEN
            maxC2     = maxC1
            i_bhe2    = i_bhe1
            i_b_maxC2 = i_b_maxC1
            i_h_maxC2 = i_h_maxC1
            i_e_maxC2 = i_e_maxC1
            i_R_maxC2 = i_R_maxC1
          END IF
          maxC1     = C
          i_bhe1    = i_bhe
          i_b_maxC1 = i_b
          i_h_maxC1 = i_h
          i_e_maxC1 = i_e
          i_R_maxC1 = i_R
        ELSE IF ( C > maxC2) THEN
          maxC2     = C
          i_bhe2    = i_bhe
          i_b_maxC2 = i_b
          i_h_maxC2 = i_h
          i_e_maxC2 = i_e
          i_R_maxC2 = i_R
        END IF
        !write(6,*) i_bhe,C,'i_b_maxC1,i_b_maxC2',i_b_maxC1,i_b_maxC2
      END DO
      END DO
      END DO
      END DO

      IF (psi%ComOp%contrac_ba_ON_HAC) THEN

        IF (psi%cplx) THEN
          write(out_unitp,*) 'max1 psi%vecBasisRep ',T,trim(info),      &
                         i_b_maxC1,i_h_maxC1,i_e_maxC1,psi%CvecB(i_bhe1)
          IF (i_b_maxC2 > 0)                                            &
               write(out_unitp,*) 'max2 psi%vecBasisRep ',T,trim(info), &
                          i_b_maxC2,i_h_maxC2,i_e_maxC2,psi%CvecB(i_bhe2)
        ELSE
          write(out_unitp,*) 'max1 psi%vecBasisRep ',T,trim(info),      &
                         i_b_maxC1,i_h_maxC1,i_e_maxC1,psi%RvecB(i_bhe1)
          IF (i_b_maxC2 > 0)                                            &
               write(out_unitp,*) 'max2 psi%vecBasisRep ',T,trim(info), &
                          i_b_maxC2,i_h_maxC2,i_e_maxC2,psi%RvecB(i_bhe2)
        END IF
      ELSE
        IF (i_b_maxC1 < 1 .OR. i_b_maxC1 > psi%nb_ba) THEN
          write(out_unitp,*) '   WARNING, i_b_maxC1 is out-of-range !! ',i_b_maxC1
          write(out_unitp,*) '    it should > 0 and <= nb_ba'
        ELSE
          CALL Rec_ndim_index(psi%BasisnD,ndim_index,i_b_maxC1)
          IF (psi%cplx) THEN
            write(out_unitp,*) 'max1 psi%vecBasisRep ',T,trim(info),      &
                        ndim_index(:),i_h_maxC1,i_e_maxC1,psi%CvecB(i_bhe1)

            IF (i_b_maxC2 > 0) THEN
              CALL Rec_ndim_index(psi%BasisnD,ndim_index,i_b_maxC2)
              write(out_unitp,*) 'max2 psi%vecBasisRep ',T,trim(info),    &
                        ndim_index(:),i_h_maxC2,i_e_maxC2,psi%CvecB(i_bhe2)
            END IF
          ELSE
            write(out_unitp,*) 'max1 psi%vecBasisRep ',T,trim(info),      &
                        ndim_index(:),i_h_maxC1,i_e_maxC1,psi%RvecB(i_bhe1)

            IF (i_b_maxC2 > 0) THEN
              CALL Rec_ndim_index(psi%BasisnD,ndim_index,i_b_maxC2)
              write(out_unitp,*) 'max2 psi%vecBasisRep ',T,trim(info),    &
                        ndim_index(:),i_h_maxC2,i_e_maxC2,psi%RvecB(i_bhe2)
            END IF

          END IF
        END IF


      END IF
      write(out_unitp,*) 'Abelian symmetry (symab):',psi%symab
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      end subroutine calc_MaxCoef_psi

      SUBROUTINE calc_1Dweight_act1(psi,ana_psi,max_1D,T,info,print_w)
      USE mod_system
      USE mod_psi_set_alloc
      USE mod_ana_psi
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_psi)     :: psi
      TYPE (param_ana_psi) :: ana_psi

      real (kind=Rkind), allocatable :: weight1Dact(:,:)
      real (kind=Rkind)    :: a

      real (kind=Rkind)    :: T ! time
      character (len=*)    :: info
      logical          :: print_w

      integer          :: i,ie,ii,ib,ibie,iq,ibiq
      integer          :: max_dim,max_1D
      integer          :: max_indGr(psi%BasisnD%nDindB%ndim)
      integer          :: ndim_AT_ib(psi%BasisnD%nDindB%ndim)
      integer          :: nDval(psi%BasisnD%nDindB%ndim)

      real (kind=Rkind):: r2

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='calc_1Dweight_act1'
      logical, parameter :: debug =.FALSE.
!     logical, parameter :: debug =.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_inact2n',psi%Basis2n%nb_basis
        write(out_unitp,*) 'nb_bi,nb_be',psi%nb_bi,psi%nb_be
      END IF
!-----------------------------------------------------------
      CALL flush_perso(out_unitp)
      IF (psi%nb_baie /= psi%nb_tot) RETURN

      IF (allocated(Psi%BasisnD%nDindB%Tab_nDval)) THEN
        max_dim = maxval(Psi%BasisnD%nDindB%Tab_nDval)
      ELSE
        max_dim = maxval(psi%BasisnD%nDindB%nDsize(1:psi%BasisnD%nDindB%ndim))
      END IF
      CALL alloc_NParray(weight1Dact,(/psi%BasisnD%nDindB%ndim,max_dim/), &
                        "weight1Dact",name_sub)

      weight1Dact(:,:) = ZERO
      ndim_AT_ib(:)    = 0

      IF (psi%cplx) THEN
        ibie = 0
        ie=1
        ii=1
        DO ib=1,psi%nb_ba
          ibie = ibie + 1
          a = abs(psi%CvecB(ibie))**2

          CALL calc_nDindex(psi%BasisnD%nDindB,ib,nDval)

          DO iq=1,psi%BasisnD%nDindB%ndim
            ibiq = nDval(iq)
            ndim_AT_ib(iq) = max(ibiq,ndim_AT_ib(iq))
            !write(out_unitp,*) 'calc_1Dweight',iq,ibiq,ib
            weight1Dact(iq,ibiq) = weight1Dact(iq,ibiq) + a
          END DO
        END DO
      ELSE
        ibie = 0
        ie=1
        ii=1
        DO ib=1,psi%nb_ba
          ibie = ibie + 1
          a = psi%RvecB(ibie)**2

          CALL calc_nDindex(psi%BasisnD%nDindB,ib,nDval)

          DO iq=1,psi%BasisnD%nDindB%ndim
            ibiq = nDval(iq)
            ndim_AT_ib(iq) = max(ibiq,ndim_AT_ib(iq))
            !write(out_unitp,*) 'calc_1Dweight',iq,ibiq,ib
            weight1Dact(iq,ibiq) = weight1Dact(iq,ibiq) + a
          END DO
        END DO
      END IF

      !write(6,*) 'ndim_AT_ib',ndim_AT_ib(:)
      IF (print_w .OR. debug) THEN
        DO iq=1,psi%BasisnD%nDindB%ndim
          IF (sum(weight1Dact(iq,1:ndim_AT_ib(iq)))-ONE > ONETENTH**7)  &
                write(out_unitp,21) 'Grd Channel Sum(RD)/=1',trim(info),&
                             iq,T,sum(weight1Dact(iq,1:ndim_AT_ib(iq)))
          write(out_unitp,21) 'Grd Channel ',trim(info),iq,T,           &
                           weight1Dact(iq,1:min(max_1D,ndim_AT_ib(iq)))
 21       format(a,a,i3,1x,f17.4,300(1x,e10.3))
          max_indGr(iq) = sum(maxloc(weight1Dact(iq,1:ndim_AT_ib(iq))))
        END DO

        !write(out_unitp,'(a,a,20i4)') 'max_indGr at ',trim(info),max_indGr(:)
        write(out_unitp,'(a,a)',advance='no') 'max_indGr at ',trim(info)
        DO i=1,size(max_indGr)-1
          write(out_unitp,'(X,i0)',advance='no') max_indGr(i)
        END DO
        write(out_unitp,'(X,i0)') max_indGr(size(max_indGr))

        CALL flush_perso(out_unitp)
      END IF

      !write(6,*) 'max_RedDensity
      IF (.NOT. allocated(ana_psi%max_RedDensity)) THEN
        CALL alloc_NParray(ana_psi%max_RedDensity,(/psi%BasisnD%nDindB%ndim/), &
                          "ana_psi%max_RedDensity",name_sub)
      END IF
      DO iq=1,psi%BasisnD%nDindB%ndim
        ana_psi%max_RedDensity(iq) = ZERO
        DO ib=1,ndim_AT_ib(iq)
          r2 = real(ndim_AT_ib(iq)-ib,kind=Rkind)**2
          ana_psi%max_RedDensity(iq) = ana_psi%max_RedDensity(iq) + weight1Dact(iq,ib)*exp(-r2)
        END DO
      END DO
      !write(out_unitp,*) 'max_RedDensity ',ana_psi%max_RedDensity(:)
      CALL Write_Vec(ana_psi%max_RedDensity,out_unitp,6,Rformat='e9.2',name_info='max_RedDensity ')
      CALL flush_perso(out_unitp)

      CALL dealloc_NParray(weight1Dact,"weight1Dact",name_sub)
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------

      END SUBROUTINE calc_1Dweight_act1
!==============================================================
!
!      calc_nDTk(psi)
!
!==============================================================
      SUBROUTINE calc_nDTk(psi,T)
      USE mod_system
      USE mod_psi_set_alloc
      IMPLICIT NONE

!----- variables for the WP propagation ----------------------------
      TYPE (param_psi)   :: psi


      real (kind=Rkind)    :: T ! time

      integer          :: k
      integer          :: iei,ib,ibie,nb_bie

      real (kind=Rkind),allocatable :: w_ei(:)

      complex (kind=Rkind),allocatable :: Ek_ei(:)
      complex (kind=Rkind),allocatable :: Ck_ei(:),Sk_ei(:)
      complex (kind=Rkind)    :: a

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='calc_nDTk'
      integer :: err_mem,memory
      logical, parameter :: debug =.FALSE.
!     logical, parameter :: debug =.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb_bi,nb_be',psi%nb_bi,psi%nb_be
      END IF
!-----------------------------------------------------------

      IF (psi%nb_baie /= psi%nb_tot) RETURN


      nb_bie = psi%nb_bi*psi%nb_be
      CALL alloc_NParray(w_ei,(/nb_bie/),"w_ei",name_sub)
      CALL alloc_NParray(Ek_ei,(/nb_bie/),"Ek_ei",name_sub)
      CALL alloc_NParray(Ck_ei,(/nb_bie/),"Ck_ei",name_sub)
      CALL alloc_NParray(Sk_ei,(/nb_bie/),"Sk_ei",name_sub)

      DO ib=1,psi%nb_ba
        k = int(ib/2)
        DO iei=1,nb_bie
          ibie = ib + (iei-1) * psi%nb_ba

          IF (psi%cplx) a = psi%CvecB(ibie)
          IF (.NOT. psi%cplx) a=cmplx(psi%RvecB(ibie),ZERO,kind=Rkind)
          IF (mod(ib-1,2) == 0) THEN
            Ck_ei(iei) = a
          ELSE
            Sk_ei(iei) = a
          END IF
        END DO
        IF (mod(ib-1,2) == 0) THEN

          Ek_ei = CHALF*(Ck_ei-EYE*Sk_ei)
          DO iei=1,nb_bie
            w_ei(iei) = abs(Ek_ei(iei))**2
          END DO
          write(998,*) 'T,+k',T,k,w_ei

          Ek_ei = CHALF*(Ck_ei+EYE*Sk_ei)
          DO iei=1,nb_bie
            w_ei(iei) = abs(Ek_ei(iei))**2
          END DO
          write(999,*) 'T,-k',T,k,w_ei

        END IF
      END DO
      write(998,*)
      write(999,*)
      CALL flush_perso(998)
      CALL flush_perso(999)

      CALL dealloc_NParray(w_ei,"w_ei",name_sub)
      CALL dealloc_NParray(Ek_ei,"Ek_ei",name_sub)
      CALL dealloc_NParray(Ck_ei,"Ck_ei",name_sub)
      CALL dealloc_NParray(Sk_ei,"Sk_ei",name_sub)

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------



      END SUBROUTINE calc_nDTk
