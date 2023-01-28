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
 MODULE mod_psi_Op
  USE mod_basis
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
        CALL ecri_psi(ZERO,psi)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------

      END SUBROUTINE Set_symab_OF_psiBasisRep
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
      IF(keep_MPI) THEN
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
          IF (psi1%nb_tot == psi1%nb_baie .AND.  locChannel_ie > 0 .AND.&
                            locChannel_ie <= psi1%nb_bi*psi1%nb_be) THEN
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

          CALL alloc_NParray(wrho,[psi1%nb_qa],"wrho",name_sub)
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
      ENDIF ! for keep_MPI

!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Overlap : ',Overlap
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------

      END SUBROUTINE Overlap_psi1_psi2

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
      USE mod_MPI_aux
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
        flush(out_unitp)
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

          CALL alloc_NParray(PsiRk,[ndim],'PsiRk',name_sub)

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

          CALL alloc_NParray(PsiRk,[ndim],'PsiRk',name_sub)

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
        flush(out_unitp)
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
!         flush(out_unitp)
        END DO
        IF (i > 1) THEN
          CALL Set_symab_OF_psiBasisRep(psi(i),sym)

          !CALL norm2_psi(psi(i))
          !write(out_unitp,*) ' Ortho: norm2',i,psi(i)%norm2
          CALL renorm_psi(psi(i))
          !write(out_unitp,*) 'symab, bits(symab)',WriteTOstring_symab(psi(i)%symab)
        END IF
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
