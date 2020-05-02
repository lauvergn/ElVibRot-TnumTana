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
      PROGRAM Tnum_f90
      USE mod_system
      USE mod_Coord_KEO
      IMPLICIT NONE

!     - parameters for para_Tnum -----------------------
      TYPE (zmatrix) :: mole
      TYPE (Tnum)    :: para_Tnum

      real (kind=Rkind) :: vep,rho
      real (kind=Rkind), pointer :: Tdef2(:,:),Tdef1(:)
      real (kind=Rkind), pointer :: Tcor2(:,:),Tcor1(:),Trot(:,:)

      integer :: nderiv,ia,ib,ic,id
      real (kind=Rkind) :: DV1opt,DV2opt
      real (kind=Rkind), pointer :: trav(:)
      real (kind=Rkind), pointer :: hess(:,:),hess_inv(:,:),hessG36(:,:),grad(:),DQinact2n(:)
      integer, pointer           :: symhessG36(:,:),trav_index(:)

      real (kind=Rkind), pointer :: d0c(:,:)
      real (kind=Rkind), pointer :: d0k(:,:),d0eh(:)
      real (kind=Rkind), pointer :: A(:,:)

      TYPE(Type_dnMat) :: dng,dnGG
      TYPE(Type_dnVec) :: dnx
      TYPE(Type_dnS)   :: dnE,dnMu(3)


!     - for the coordinate values ----------------------------------
      TYPE (param_Q)    :: para_Q
      TYPE (PrimOp_t)   :: PrimOp


!     - working parameters ------------------------------------------
      integer :: i,j,i1,i1i,i1f,i2,i2i,i2f,i12
      integer :: err_mem,memory

!=======================================================================
!=======================================================================
      CALL versionEVRT(.TRUE.)

      !-----------------------------------------------------------------
      !     - read the coordinate tansformations :
      !     -   zmatrix, polysperical, bunch...
      !     ------------------------------------------------------------
      CALL Read_mole(mole,para_Tnum)
      !     ------------------------------------------------------------
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !     - read coordinate values -----------------------------------
      !     ------------------------------------------------------------
      CALL read_RefGeom(para_Q,mole,para_Tnum)
      !     ------------------------------------------------------------
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      !     ---- TO finalize the coordinates (NM) and the KEO ----------
      !     ------------------------------------------------------------
      CALL Finalize_TnumTana_Coord_PrimOp(para_Q%Qact,para_Tnum,mole,PrimOp)
      !-----------------------------------------------------------------
!=======================================================================
!=======================================================================

!-------------------------------------------------
!     - Cartesian coordinates --------------------
!     --------------------------------------------

       write(out_unitp,*) "======================================"
       write(out_unitp,*) "======================================"
       write(out_unitp,*) "======================================"
       write(out_unitp,*) "======================================"
       nderiv = 0
       CALL alloc_dnSVM(dnx,mole%ncart,mole%nb_act,nderiv)
       write(out_unitp,*) "======================================"
       CALL time_perso('dnx')
       write(out_unitp,*) "======================================"
       CALL sub3_QTOdnx(dnx,para_Q,mole,nderiv,.FALSE.)
       write(out_unitp,*) 'dnx: ',mole%ncart
       CALL write_dnx(1,mole%ncart,dnx,nderiv)
       CALL dealloc_dnSVM(dnx)
       write(out_unitp,*) "======================================"
       write(out_unitp,*) "======================================"
       write(out_unitp,*) "======================================"
       write(out_unitp,*) "======================================"


!===========================================================
!===========================================================

!-------------------------------------------------
!     - Potential : on the fly calculation -------
!     --------------------------------------------
        nderiv = 2
        PrimOp%nb_elec = 1
        PrimOp%OnTheFly = .TRUE.
        PrimOp%Read_OnTheFly_only = .TRUE.
        write(out_unitp,*) 'on the fly calculation'
        CALL alloc_dnSVM(dnE,mole%nb_act,nderiv)
        CALL alloc_dnSVM(dnMu(1),mole%nb_act,nderiv-1)
        CALL alloc_dnSVM(dnMu(2),mole%nb_act,nderiv-1)
        CALL alloc_dnSVM(dnMu(3),mole%nb_act,nderiv-1)


        memory = product( (/ mole%nb_inact2n,mole%nb_inact2n /) )
        allocate(hess(mole%nb_inact2n,mole%nb_inact2n),stat=err_mem) ! change alloc done
        CALL error_memo_allo(err_mem,memory,"hess","main")
        memory = product( (/ mole%nb_inact2n,mole%nb_inact2n /) )
        allocate(hess_inv(mole%nb_inact2n,mole%nb_inact2n),stat=err_mem) ! change alloc done
        CALL error_memo_allo(err_mem,memory,"hess_inv","main")
        memory = product( (/ mole%nb_inact2n /) )
        allocate(trav(mole%nb_inact2n),stat=err_mem) ! change alloc done
        CALL error_memo_allo(err_mem,memory,"trav","main")
        memory = product( (/ mole%nb_inact2n /) )
        allocate(grad(mole%nb_inact2n),stat=err_mem) ! change alloc done
        CALL error_memo_allo(err_mem,memory,"grad","main")
        memory = product( (/ mole%nb_inact2n /) )
        allocate(dQinact2n(mole%nb_inact2n),stat=err_mem) ! change alloc done
        CALL error_memo_allo(err_mem,memory,"dQinact2n","main")
        memory = product( (/ mole%nb_inact2n /) )
        allocate(trav_index(mole%nb_inact2n),stat=err_mem) ! change alloc done
        CALL error_memo_allo(err_mem,memory,"trav_index","main")


        memory = product( (/ mole%nb_act,mole%nb_act /) )
        allocate(hessG36(mole%nb_act,mole%nb_act),stat=err_mem) ! change alloc done
        CALL error_memo_allo(err_mem,memory,"hessG36","main")
        memory = product( (/ mole%nb_act,mole%nb_act /) )
        allocate(symhessG36(mole%nb_act,mole%nb_act),stat=err_mem) ! change alloc done
        CALL error_memo_allo(err_mem,memory,"symhessG36","main")

        CALL NdnOp_grid(para_Q%Qact,dnE,nderiv,dnMu,nderiv-1,           &
                        mole,para_Tnum,PrimOp)

        dnE%d2(:,:) = (dnE%d2+transpose(dnE%d2)) * HALF ! symetrization


        write(out_unitp,*) 'dnE%d0',mole%ActiveTransfo%Qdyn0,dnE%d0
        write(out_unitp,*) 'Dipole Moment',dnMu%d0

        IF (nderiv > 0) THEN
          DO i=1,mole%nb_act
            write(out_unitp,*) 'dnE%d1',mole%ActiveTransfo%Qact0(i),dnE%d1(i)
          END DO
        END IF
        !write(out_unitp,*) 'Grad',dnE%d1(mole%nb_act1+1:mole%nb_act)

!        IF (nderiv == 2) CALL Write_Mat(dnE%d2,out_unitp,5)


        hess(:,:) = dnE%d2(mole%nb_act1+1:mole%nb_act,mole%nb_act1+1:mole%nb_act)
        CALL inversion(hess_inv,hess,trav,trav_index,mole%nb_inact2n)
        hess(:,:) = dnE%d2(mole%nb_act1+1:mole%nb_act,mole%nb_act1+1:mole%nb_act)
        grad(:)   = dnE%d1(mole%nb_act1+1:mole%nb_act)

        dQinact2n(:) = - matmul(hess_inv,grad)

!        DO i=1,mole%nb_inact2n
!          write(out_unitp,*) 'dQinact2n',i,dQinact2n(i)
!        END DO

!        DV1opt = dot_product(dQinact2n,grad) +                          &
!          HALF * dot_product(dQinact2n,matmul(hess,dQinact2n))
!        grad = grad + matmul(hess,dQinact2n)
!        DO i=1,mole%nb_inact2n
!          write(out_unitp,*) 'new grad',i,para_Q%Qact(mole%nb_act1+i),grad(i)
!        END DO

        DV1opt = HALF* dot_product(grad,dQinact2n) ! it's equivalent

        para_Q%Qact(mole%nb_act1+1:mole%nb_act) =                       &
             para_Q%Qact0(mole%nb_act1+1:mole%nb_act) + dQinact2n(:)

         CALL symetrization_G36(dnE%d2,hessG36,symhessG36,mole%nb_act)
         i1i=8
         i1f=15
         i2i=8
         i2f=15

         i12 = 2
         DO i1=i1i,i1f
         DO i2=i2i,i2f
           i12 = i12 +1
           write(out_unitp,*) 'dat_fit_v5G36_HessG ',i12,' 0,1,',       &
                         int_TO_char(symhessG36(i2,i1)), hessG36(i2,i1)
         END DO
         END DO

        write(out_unitp,*) 'Vopt',dnE%d0,DV1opt*mole%const_phys%auTOcm_inv,dnE%d0+DV1opt
        write(out_unitp,*) 'Qopt',para_Q%Qact(mole%nb_act1+1:mole%nb_act)
        !write(out_unitp,*) 'hessG36 ',hessG36(i2i:i2f,i1i:i1f)
!       CALL Write_Mat(hessG36,out_unitp,5)


      CALL dealloc_zmat(mole)
      CALL dealloc_param_Q(para_Q)


      write(out_unitp,*) 'END Tnum'

      end program Tnum_f90
!=======================================================================
!     symetrization hess => hessG36
!     The coordinates has to be the symetrized ones and the must be sorted as followed:
!     each components of degenerated symetry must be be together:
!     Eia,Eib and Ga, Gb, Gc, Gd
!
!     sym:  A1,A2,A3,A4, E1a,E1b, E2a,E2b, E3a,E3b, E4a,E4b, Ga,Gb,Gc,Gd
!     #sym:  1  2  3  4   5   6    7   8    9  10   11  12   13 14 15 16
!=======================================================================
      SUBROUTINE symetrization_G36(hess,hessG36,symhessG36,nb_coord)
      USE mod_system
      IMPLICIT NONE

      integer :: nb_coord
      real (kind=Rkind), intent(in)  :: hess(nb_coord,nb_coord)
      real (kind=Rkind), intent(inout) :: hessG36(nb_coord,nb_coord)
      integer, intent(inout)           :: symhessG36(nb_coord,nb_coord)


      integer :: sym_coord(nb_coord)

      real (kind=Rkind),allocatable :: blocG36_FROM_blocS1xS2(:,:,:),blocG36(:)
      integer,allocatable :: symblocG36_FROM_blocS1xS2(:)

      integer :: degen_FROM_symG36(0:16) = (/0, 1,1,1,1,  2,0,2,0,2,0,2,0,  4,0,0,0 /)
      integer :: sym1,sym2,n1,n2,i1,i2,k12,i
      integer :: err_mem,memory

!     read the symetry of each coordinate
      read(in_unitp,*) sym_coord(:)

!     check that coordinates are correctly sorted
      DO i=1,nb_coord

        IF (sym_coord(i) < 0 .OR. sym_coord(i) > 16) THEN
          write(out_unitp,*) ' ERROR in symetrization_G36'
          write(out_unitp,*) ' sym_coord(i) is out of range',i,sym_coord(i)
          write(out_unitp,*) ' sym_coord(i) sould be in the range [1...16]'
          STOP
        END IF

        IF (sym_coord(i) == 5 .OR. sym_coord(i) == 7 .OR.          &
                 sym_coord(i) == 9 .OR. sym_coord(i) == 11 ) THEN   ! E1,E2,E3,E4
          IF (sym_coord(i+1) /= sym_coord(i)+1 ) THEN
            write(out_unitp,*) ' ERROR in symetrization_G36'
            write(out_unitp,*) ' sym_coord(:) are not sorted (E symetry)'
            write(out_unitp,*) ' i',i,'sym_coord(i:i+1)',sym_coord(i:i+1)
            STOP
          END IF
        ELSE IF (sym_coord(i) == 13) THEN  !G
          IF (sym_coord(i+1) /= sym_coord(i)+1 .AND.                    &
              sym_coord(i+2) /= sym_coord(i)+2 .AND.                    &
              sym_coord(i+3) /= sym_coord(i)+3            ) THEN
            write(out_unitp,*) ' ERROR in symetrization_G36'
            write(out_unitp,*) ' sym_coord(:) are not sorted (G symetry)'
            write(out_unitp,*) ' i',i,'sym_coord(i:i+1)',sym_coord(i:i+1)
            STOP
          END IF
        ELSE ! the other symetries (6,8,10,12 and 14,15,16) have been already tested
          CONTINUE
        END IF
      END DO

!     symetrization of hess => hessG36
      hessG36(:,:) = ZERO
       DO i1=1,nb_coord
       !DO i1=18,18

        sym1 = sym_coord(i1)
        n1 = degen_FROM_symG36(sym1)
        IF (n1 == 0) CYCLE

        DO i2=1,nb_coord

          sym2 = sym_coord(i2)
          n2 = degen_FROM_symG36(sym2)
          IF (n2 == 0) CYCLE
          !write(out_unitp,*) i1,i2,sym1,sym2,n1,n2
          memory = product( (/ n1*n2 /) )
          allocate(symblocG36_FROM_blocS1xS2(n1*n2),stat=err_mem) ! change alloc done
          CALL error_memo_allo(err_mem,memory,                          &
                        "symblocG36_FROM_blocS1xS2","symetrization_G36")
          memory = product( (/ n1*n2,n1,n2 /) )
          allocate(blocG36_FROM_blocS1xS2(n1*n2,n1,n2),stat=err_mem) ! change alloc done
          CALL error_memo_allo(err_mem,memory,"blocG36_FROM_blocS1xS2", &
                                                    "symetrization_G36")
          memory = product( (/ n1*n2 /) )
          allocate(blocG36(n1*n2),stat=err_mem) ! change alloc done
          CALL error_memo_allo(err_mem,memory,"blocG36",                &
                                                    "symetrization_G36")
          CALL init_blocG36_FROM_blocS1xS2(blocG36_FROM_blocS1xS2,      &
                                           symblocG36_FROM_blocS1xS2,   &
                                           sym1,sym2,n1,n2)

          DO k12=1,n1*n2
            blocG36(k12) = sum(hess(i1:i1-1+n1,i2:i2-1+n2)*blocG36_FROM_blocS1xS2(k12,:,:))
          END DO

          hessG36(i1:i1-1+n1,i2:i2-1+n2)    = reshape(blocG36,(/ n1,n2 /) )
          symhessG36(i1:i1-1+n1,i2:i2-1+n2) = reshape(symblocG36_FROM_blocS1xS2,(/ n1,n2 /) )
          !write(out_unitp,*) 'blocG36',blocG36(:)
          memory = size(symblocG36_FROM_blocS1xS2)
          deallocate(symblocG36_FROM_blocS1xS2,stat=err_mem) ! change dealloc done
          CALL error_memo_allo(err_mem,-memory,                         &
                        "symblocG36_FROM_blocS1xS2","symetrization_G36")
          memory = size(blocG36_FROM_blocS1xS2)
          deallocate(blocG36_FROM_blocS1xS2,stat=err_mem) ! change dealloc done
          CALL error_memo_allo(err_mem,-memory,"blocG36_FROM_blocS1xS2", &
                                                    "symetrization_G36")
          memory = size(blocG36)
          deallocate(blocG36,stat=err_mem) ! change dealloc done
          CALL error_memo_allo(err_mem,-memory,"blocG36",               &
                                                    "symetrization_G36")
        END DO
      END DO

      end subroutine symetrization_G36
      SUBROUTINE init_blocG36_FROM_blocS1xS2(blocG36_FROM_blocS1xS2,    &
                                             symblocG36_FROM_blocS1xS2, &
                                             sym1,sym2,n1,n2)
      USE mod_system
      IMPLICIT NONE

      integer, intent(in) :: sym1,sym2,n1,n2

      real (kind=Rkind), intent(inout)  :: blocG36_FROM_blocS1xS2(n1*n2,n1,n2)
      integer, intent(inout)            :: symblocG36_FROM_blocS1xS2(n1*n2)

      integer :: i,is,ia,ib,ic,id

      IF ( (sym1 >= 1  .AND. sym1 <= 4  .AND. n1 == 1) .OR.             &
           (sym1 >= 5  .AND. sym1 <= 12 .AND. n1 == 2) .OR.             &
           (sym1 >= 13 .AND. sym1 <= 16 .AND. n1 == 4) ) THEN
        CONTINUE
      ELSE
        write(out_unitp,*) ' ERROR in init_blocG36_FROM_blocS1xS2'
        write(out_unitp,*) ' sym1 or n1 are out of range! ',sym1,n1
        write(out_unitp,*) ' sym1 E [1... 4] (sym Ai)  => n1 = 1 '
        write(out_unitp,*) ' sym1 E [5... 12] (sym Ei) => n1 = 2 '
        write(out_unitp,*) ' sym1 E [13...16] (sym Gi) => n1 = 4 '
        STOP
      END IF

      IF ( (sym2 >= 1  .AND. sym2 <= 4  .AND. n2 == 1) .OR.             &
           (sym2 >= 5  .AND. sym2 <= 12 .AND. n2 == 2) .OR.             &
           (sym2 >= 13 .AND. sym2 <= 16 .AND. n2 == 4) ) THEN
        CONTINUE
      ELSE
        write(out_unitp,*) ' ERROR in init_blocG36_FROM_blocS1xS2'
        write(out_unitp,*) ' sym2 or n2 are out of range! ',sym2,n2
        write(out_unitp,*) ' sym2 E [1... 4] (sym Ai)  => n2 = 1 '
        write(out_unitp,*) ' sym2 E [5... 12] (sym Ei) => n2 = 2 '
        write(out_unitp,*) ' sym2 E [13...16] (sym Gi) => n2 = 4 '
        STOP
      END IF

      blocG36_FROM_blocS1xS2(:,:,:) = ZERO
      symblocG36_FROM_blocS1xS2(:)  = 0
      IF (sym1 == 1) THEN       ! A1

        IF (sym2 == 1) THEN       ! A1xA1
          blocG36_FROM_blocS1xS2(1,1,1) = ONE
          symblocG36_FROM_blocS1xS2(1)  = 1 ! A1
        ELSE IF (sym2 == 2) THEN  ! A1xA2
          blocG36_FROM_blocS1xS2(1,1,1) = ONE
          symblocG36_FROM_blocS1xS2(1)  = 2 ! A2
        ELSE IF (sym2 == 3) THEN  ! A1xA3
          blocG36_FROM_blocS1xS2(1,1,1) = ONE
          symblocG36_FROM_blocS1xS2(1)  = 3 ! A3
        ELSE IF (sym2 == 4) THEN  ! A1xA4
          blocG36_FROM_blocS1xS2(1,1,1) = ONE
          symblocG36_FROM_blocS1xS2(1)  = 1 ! A4
        ELSE IF (sym2 == 13) THEN ! A1xG
          blocG36_FROM_blocS1xS2(1,1,1)  = ONE
          blocG36_FROM_blocS1xS2(2,1,2)  = ONE
          blocG36_FROM_blocS1xS2(3,1,3)  = ONE
          blocG36_FROM_blocS1xS2(4,1,4)  = ONE
          symblocG36_FROM_blocS1xS2(1:4) = (/ 13,14,15,16 /) ! Ga,Gb,Gc,Gd
        ELSE
          STOP 'E sym are not implemented !'
        END IF

      ELSE IF (sym1 == 2) THEN  ! A2
          STOP 'A2 sym are not implemented !'

      ELSE IF (sym1 == 3) THEN  ! A3
          STOP 'A3 sym are not implemented !'

      ELSE IF (sym1 == 4) THEN  ! A4
        IF (sym2 == 1) THEN       ! A4xA1
          blocG36_FROM_blocS1xS2(1,1,1) = ONE
          symblocG36_FROM_blocS1xS2(1)  = 4 ! A4
        ELSE IF (sym2 == 2) THEN  ! A4xA2
          STOP 'A2 sym are not implemented !'
        ELSE IF (sym2 == 3) THEN  ! A4xA3
          STOP 'A2 sym are not implemented !'
        ELSE IF (sym2 == 4) THEN  ! A4xA4
          blocG36_FROM_blocS1xS2(1,1,1) = ONE
          symblocG36_FROM_blocS1xS2(1)  = 1 ! A1
        ELSE IF (sym2 == 13) THEN ! A4xG
          blocG36_FROM_blocS1xS2(1,1,2)  = ONE
          blocG36_FROM_blocS1xS2(2,1,1)  = ONE
          blocG36_FROM_blocS1xS2(3,1,4)  = ONE
          blocG36_FROM_blocS1xS2(4,1,3)  = ONE
          symblocG36_FROM_blocS1xS2(1:4) = (/ 13,14,15,16 /) ! Ga,Gb,Gc,Gd
        ELSE
          STOP 'E sym are not implemented !'
        END IF

      ELSE IF (sym1 == 13) THEN ! G
        IF (sym2 == 1) THEN       ! A1
          blocG36_FROM_blocS1xS2(1,1,1)  = ONE
          blocG36_FROM_blocS1xS2(2,2,1)  = ONE
          blocG36_FROM_blocS1xS2(3,3,1)  = ONE
          blocG36_FROM_blocS1xS2(4,4,1)  = ONE
          symblocG36_FROM_blocS1xS2(1:4) = (/ 13,14,15,16 /) ! Ga,Gb,Gc,Gd
        ELSE IF (sym2 == 2) THEN  ! A2
          STOP 'A2 sym are not implemented !'
        ELSE IF (sym2 == 3) THEN  ! A3
          STOP 'A2 sym are not implemented !'
        ELSE IF (sym2 == 4) THEN  ! A4
          blocG36_FROM_blocS1xS2(1,2,1)  = ONE
          blocG36_FROM_blocS1xS2(2,1,1)  = ONE
          blocG36_FROM_blocS1xS2(3,4,1)  = ONE
          blocG36_FROM_blocS1xS2(4,3,1)  = ONE
          symblocG36_FROM_blocS1xS2(1:4) = (/ 13,14,15,16 /) ! Ga,Gb,Gc,Gd
        ELSE IF (sym2 == 13) THEN ! G
         ia=1
         ib=2
         ic=3
         id=4

         is=1 !GxG => sym A1
         blocG36_FROM_blocS1xS2(is,ia,ia) = HALF
         blocG36_FROM_blocS1xS2(is,ib,ib) = HALF
         blocG36_FROM_blocS1xS2(is,ic,ic) = HALF
         blocG36_FROM_blocS1xS2(is,id,id) = HALF

         is=2 !GxG => sym A2
         blocG36_FROM_blocS1xS2(is,ia,ib) = HALF
         blocG36_FROM_blocS1xS2(is,ib,ia) =-HALF
         blocG36_FROM_blocS1xS2(is,ic,id) =-HALF
         blocG36_FROM_blocS1xS2(is,id,ic) = HALF

         is=3 !GxG => sym A3
         blocG36_FROM_blocS1xS2(is,ia,ib) = HALF
         blocG36_FROM_blocS1xS2(is,ib,ia) =-HALF
         blocG36_FROM_blocS1xS2(is,ic,id) = HALF
         blocG36_FROM_blocS1xS2(is,id,ic) =-HALF

         is=4 !GxG => sym A4
         blocG36_FROM_blocS1xS2(is,ia,ia) = HALF
         blocG36_FROM_blocS1xS2(is,ib,ib) = HALF
         blocG36_FROM_blocS1xS2(is,ic,ic) =-HALF
         blocG36_FROM_blocS1xS2(is,id,id) =-HALF

         is=5 !GxG => sym E1a
         blocG36_FROM_blocS1xS2(is,ia,ic) = HALF
         blocG36_FROM_blocS1xS2(is,ic,ia) = HALF
         blocG36_FROM_blocS1xS2(is,ib,id) = HALF
         blocG36_FROM_blocS1xS2(is,id,ib) = HALF

         is=6 ! GxG => sym E1b
         blocG36_FROM_blocS1xS2(is,ia,id) = HALF
         blocG36_FROM_blocS1xS2(is,id,ia) = HALF
         blocG36_FROM_blocS1xS2(is,ib,ic) =-HALF
         blocG36_FROM_blocS1xS2(is,ic,ib) =-HALF

         is=7 !GxG => sym E2a
         blocG36_FROM_blocS1xS2(is,ia,ic) = HALF
         blocG36_FROM_blocS1xS2(is,ic,ia) =-HALF
         blocG36_FROM_blocS1xS2(is,ib,id) = HALF
         blocG36_FROM_blocS1xS2(is,id,ib) =-HALF

         is=8 !GxG => sym E2b
         blocG36_FROM_blocS1xS2(is,ia,id) = HALF
         blocG36_FROM_blocS1xS2(is,id,ia) =-HALF
         blocG36_FROM_blocS1xS2(is,ib,ic) =-HALF
         blocG36_FROM_blocS1xS2(is,ic,ib) = HALF


         is=9 !GxG => sym E3a
         blocG36_FROM_blocS1xS2(is,ia,ic) = HALF
         blocG36_FROM_blocS1xS2(is,ic,ia) = HALF
         blocG36_FROM_blocS1xS2(is,ib,id) =-HALF
         blocG36_FROM_blocS1xS2(is,id,ib) =-HALF

         is=10 !GxG => sym E3b
         blocG36_FROM_blocS1xS2(is,ia,id) = HALF
         blocG36_FROM_blocS1xS2(is,id,ia) = HALF
         blocG36_FROM_blocS1xS2(is,ib,ic) = HALF
         blocG36_FROM_blocS1xS2(is,ic,ib) = HALF

         is=11 !GxG => sym E4a
         blocG36_FROM_blocS1xS2(is,ia,ic) = HALF
         blocG36_FROM_blocS1xS2(is,ic,ia) =-HALF
         blocG36_FROM_blocS1xS2(is,ib,id) =-HALF
         blocG36_FROM_blocS1xS2(is,id,ib) = HALF

         is=12 !GxG => sym E4b
         blocG36_FROM_blocS1xS2(is,ia,id) = HALF
         blocG36_FROM_blocS1xS2(is,id,ia) =-HALF
         blocG36_FROM_blocS1xS2(is,ib,ic) = HALF
         blocG36_FROM_blocS1xS2(is,ic,ib) =-HALF

         is=13 !GxG => sym Ga
         blocG36_FROM_blocS1xS2(is,ib,ib) = ONE/sqrt(TWO)
         blocG36_FROM_blocS1xS2(is,ia,ia) =-ONE/sqrt(TWO)

         is=14 ! GxG => sym Gb
         blocG36_FROM_blocS1xS2(is,ia,ib) = ONE/sqrt(TWO)
         blocG36_FROM_blocS1xS2(is,ib,ia) = ONE/sqrt(TWO)

         is=15 !GxG => sym Gc
         blocG36_FROM_blocS1xS2(is,id,id) = ONE/sqrt(TWO)
         blocG36_FROM_blocS1xS2(is,ic,ic) =-ONE/sqrt(TWO)

         is=16 !GxG => sym Gd
         blocG36_FROM_blocS1xS2(is,ic,id) = ONE/sqrt(TWO)
         blocG36_FROM_blocS1xS2(is,id,ic) = ONE/sqrt(TWO)

         symblocG36_FROM_blocS1xS2(1:16) = (/ (i,i=1,16) /) ! A1...., Ga,Gb,Gc,Gd

        ELSE
          STOP 'E sym are not implemented !'
        END IF


      ELSE
        STOP 'E sym are not implemented !'
      END IF

!      DO is=1,n1*n2
!        write(out_unitp,*) 'sym:',symblocG36_FROM_blocS1xS2(is)
!        CALL Write_Mat(blocG36_FROM_blocS1xS2(is,:,:),out_unitp,5)
!      END DO


      end subroutine init_blocG36_FROM_blocS1xS2
      SUBROUTINE symetrization_G36_old(hess,hessG36,mat_symG36,nb_coord)
      USE mod_system
      IMPLICIT NONE

      integer :: nb_coord
      real (kind=Rkind), intent(in)  :: hess(nb_coord,nb_coord)
      real (kind=Rkind), intent(inout) :: hessG36(nb_coord,nb_coord)
      integer, intent(inout)           :: mat_symG36(nb_coord,nb_coord)


      integer :: sym_coord(nb_coord),degen_coord(nb_coord)

      real (kind=Rkind) :: blocG36_FROM_blocGxG(16,4,4)
      integer :: symblocG36_FROM_blocGxG(16)

      real (kind=Rkind) :: blocG36_FROM_blocA1xG(4,1,4)
      real (kind=Rkind) :: blocG36_FROM_blocGxA1(4,4,1)
      real (kind=Rkind) :: blocG36_FROM_blocA4xG(4,1,4)
      real (kind=Rkind) :: blocG36_FROM_blocGxA4(4,4,1)
      integer :: symblocG36_FROM_blocA1xG(4)
      integer :: symblocG36_FROM_blocGxA1(4)
      integer :: symblocG36_FROM_blocA4xG(4)
      integer :: symblocG36_FROM_blocGxA4(4)

      real (kind=Rkind) :: blocG36_FROM_blocA1xA1(1,1,1)
      real (kind=Rkind) :: blocG36_FROM_blocA1xA4(1,1,1)
      real (kind=Rkind) :: blocG36_FROM_blocA4xA1(1,1,1)
      real (kind=Rkind) :: blocG36_FROM_blocA4xA4(1,1,1)
      integer :: symblocG36_FROM_blocA1xA1(1)
      integer :: symblocG36_FROM_blocA1xA4(1)
      integer :: symblocG36_FROM_blocA4xA1(1)
      integer :: symblocG36_FROM_blocA4xA4(1)

      integer :: i,is,ia,ib,ic,id

!     read the symetry of each coordinate
      read(in_unitp,*) sym_coord(:)
      degen_coord(:) = 0

!     check that coordinates are correctly sorted
!     affect the degeneracy for each coordinates
      DO i=1,nb_coord

        IF (sym_coord(i) < 1 .OR. sym_coord(i) > 16) THEN
          write(out_unitp,*) ' ERROR in symetrization_G36'
          write(out_unitp,*) ' sym_coord(i) is out of range',i,sym_coord(i)
          write(out_unitp,*) ' sym_coord(i) sould be in the range [1...16]'
          STOP
        END IF

        IF (sym_coord(i) >= 1 .OR. sym_coord(i) <=4) THEN  ! A1,A2,A3,A4
          degen_coord(i) = 1
        ELSE IF (sym_coord(i) >= 5 .OR. sym_coord(i) >= 7 .OR.          &
                 sym_coord(i) >= 9 .OR. sym_coord(i) >= 11 ) THEN   ! E1,E2,E3,E4
          degen_coord(i) = 2
          IF (sym_coord(i+1) /= sym_coord(i)+1 ) THEN
            write(out_unitp,*) ' ERROR in symetrization_G36'
            write(out_unitp,*) ' sym_coord(:) are not sorted (E symetry)'
            write(out_unitp,*) ' i',i,'sym_coord(i:i+1)',sym_coord(i:i+1)
            STOP
          END IF
        ELSE IF (sym_coord(i) >= 13) THEN  !G
          degen_coord(i) = 4
          IF (sym_coord(i+1) /= sym_coord(i)+1 .AND.                    &
              sym_coord(i+2) /= sym_coord(i)+2 .AND.                    &
              sym_coord(i+3) /= sym_coord(i)+3            ) THEN
            write(out_unitp,*) ' ERROR in symetrization_G36'
            write(out_unitp,*) ' sym_coord(:) are not sorted (G symetry)'
            write(out_unitp,*) ' i',i,'sym_coord(i:i+1)',sym_coord(i:i+1)
            STOP
          END IF
        ELSE ! the other symetries (6,8,10,12 and 14,15,16) have been already tested
          CONTINUE
        END IF

      END DO

!     initialization of blocS1xS2_TO_blocG36
!     A1xA1 A1xA4, A4xA1, A4xA4
      blocG36_FROM_blocA1xA1(1,1,1) = ONE
      symblocG36_FROM_blocA1xA1(1) = 1 ! A1

      blocG36_FROM_blocA1xA4(1,1,1) = ONE
      symblocG36_FROM_blocA1xA4(1) = 4 ! A4

      blocG36_FROM_blocA4xA1(1,1,1) = ONE
      symblocG36_FROM_blocA4xA1(1) = 4 ! A4

      blocG36_FROM_blocA4xA4(1,1,1) = ONE
      symblocG36_FROM_blocA4xA4(1) = 1 ! A1

!     A1xG
      blocG36_FROM_blocA1xG(1,1,:) = (/ ONE, ZERO, ZERO, ZERO /)
      blocG36_FROM_blocA1xG(2,1,:) = (/ ZERO, ONE, ZERO, ZERO /)
      blocG36_FROM_blocA1xG(3,1,:) = (/ ZERO, ZERO, ONE, ZERO /)
      blocG36_FROM_blocA1xG(4,1,:) = (/ ZERO, ZERO, ZERO, ONE /)
      symblocG36_FROM_blocA1xG(1:4) = (/ 13,14,15,16 /) ! Ga,Gb,Gc,Gd

!     GxA1
      blocG36_FROM_blocGxA1(1,:,1) = (/ ONE, ZERO, ZERO, ZERO /)
      blocG36_FROM_blocGxA1(2,:,1) = (/ ZERO, ONE, ZERO, ZERO /)
      blocG36_FROM_blocGxA1(3,:,1) = (/ ZERO, ZERO, ONE, ZERO /)
      blocG36_FROM_blocGxA1(4,:,1) = (/ ZERO, ZERO, ZERO, ONE /)
      symblocG36_FROM_blocGxA1(1:4) = (/ 13,14,15,16 /) ! Ga,Gb,Gc,Gd

!     A4xG
      blocG36_FROM_blocA4xG(1,1,:) = (/ ZERO, ONE, ZERO, ZERO /)
      blocG36_FROM_blocA4xG(2,1,:) = (/ ONE, ZERO, ZERO, ZERO /)
      blocG36_FROM_blocA4xG(3,1,:) = (/ ZERO, ZERO, ZERO, ONE /)
      blocG36_FROM_blocA4xG(4,1,:) = (/ ZERO, ZERO, ONE, ZERO /)
      symblocG36_FROM_blocA4xG(1:4) = (/ 13,14,15,16 /) ! Ga,Gb,Gc,Gd

!     GxA4
      blocG36_FROM_blocGxA4(1,:,1) = (/ ZERO, ONE, ZERO, ZERO /)
      blocG36_FROM_blocGxA4(2,:,1) = (/ ONE, ZERO, ZERO, ZERO /)
      blocG36_FROM_blocGxA4(3,:,1) = (/ ZERO, ZERO, ZERO, ONE /)
      blocG36_FROM_blocGxA4(4,:,1) = (/ ZERO, ZERO, ONE, ZERO /)
      symblocG36_FROM_blocGxA1(1:4) = (/ 13,14,15,16 /) ! Ga,Gb,Gc,Gd

!     GxG
      blocG36_FROM_blocGxG(:,:,:) = ZERO
      ia=1
      ib=2
      ic=3
      id=4

      write(out_unitp,*) ' GxG => sym A1'
      is=1
      blocG36_FROM_blocGxG(is,ia,ia) = HALF
      blocG36_FROM_blocGxG(is,ib,ib) = HALF
      blocG36_FROM_blocGxG(is,ic,ic) = HALF
      blocG36_FROM_blocGxG(is,id,id) = HALF
      CALL Write_Mat(blocG36_FROM_blocGxG(is,:,:),out_unitp,5)

      write(out_unitp,*) ' GxG => sym A2'
      is=2
      blocG36_FROM_blocGxG(is,ia,ib) = HALF
      blocG36_FROM_blocGxG(is,ib,ia) =-HALF
      blocG36_FROM_blocGxG(is,ic,id) =-HALF
      blocG36_FROM_blocGxG(is,id,ic) = HALF
      CALL Write_Mat(blocG36_FROM_blocGxG(is,:,:),out_unitp,5)

      write(out_unitp,*) ' GxG => sym A3'
      is=3
      blocG36_FROM_blocGxG(is,ia,ib) = HALF
      blocG36_FROM_blocGxG(is,ib,ia) =-HALF
      blocG36_FROM_blocGxG(is,ic,id) = HALF
      blocG36_FROM_blocGxG(is,id,ic) =-HALF
      CALL Write_Mat(blocG36_FROM_blocGxG(is,:,:),out_unitp,5)

      write(out_unitp,*) ' GxG => sym A4'
      is=4
      blocG36_FROM_blocGxG(is,ia,ia) = HALF
      blocG36_FROM_blocGxG(is,ib,ib) = HALF
      blocG36_FROM_blocGxG(is,ic,ic) =-HALF
      blocG36_FROM_blocGxG(is,id,id) =-HALF
      CALL Write_Mat(blocG36_FROM_blocGxG(is,:,:),out_unitp,5)

      write(out_unitp,*) ' GxG => sym E1a'
      is=5
      blocG36_FROM_blocGxG(is,ia,ic) = HALF
      blocG36_FROM_blocGxG(is,ic,ia) = HALF
      blocG36_FROM_blocGxG(is,ib,id) = HALF
      blocG36_FROM_blocGxG(is,id,ib) = HALF
      CALL Write_Mat(blocG36_FROM_blocGxG(is,:,:),out_unitp,5)
      write(out_unitp,*) ' GxG => sym E1b'
      is=6
      blocG36_FROM_blocGxG(is,ia,id) = HALF
      blocG36_FROM_blocGxG(is,id,ia) = HALF
      blocG36_FROM_blocGxG(is,ib,ic) =-HALF
      blocG36_FROM_blocGxG(is,ic,ib) =-HALF
      CALL Write_Mat(blocG36_FROM_blocGxG(is,:,:),out_unitp,5)

      write(out_unitp,*) ' GxG => sym E2a'
      is=7
      blocG36_FROM_blocGxG(is,ia,ic) = HALF
      blocG36_FROM_blocGxG(is,ic,ia) =-HALF
      blocG36_FROM_blocGxG(is,ib,id) = HALF
      blocG36_FROM_blocGxG(is,id,ib) =-HALF
      CALL Write_Mat(blocG36_FROM_blocGxG(is,:,:),out_unitp,5)
      write(out_unitp,*) ' GxG => sym E2b'
      is=8
      blocG36_FROM_blocGxG(is,ia,id) = HALF
      blocG36_FROM_blocGxG(is,id,ia) =-HALF
      blocG36_FROM_blocGxG(is,ib,ic) =-HALF
      blocG36_FROM_blocGxG(is,ic,ib) = HALF
      CALL Write_Mat(blocG36_FROM_blocGxG(is,:,:),out_unitp,5)

      write(out_unitp,*) ' GxG => sym E3a'
      is=9
      blocG36_FROM_blocGxG(is,ia,ic) = HALF
      blocG36_FROM_blocGxG(is,ic,ia) = HALF
      blocG36_FROM_blocGxG(is,ib,id) =-HALF
      blocG36_FROM_blocGxG(is,id,ib) =-HALF
      CALL Write_Mat(blocG36_FROM_blocGxG(is,:,:),out_unitp,5)
      write(out_unitp,*) ' GxG => sym E3b'
      is=10
      blocG36_FROM_blocGxG(is,ia,id) = HALF
      blocG36_FROM_blocGxG(is,id,ia) = HALF
      blocG36_FROM_blocGxG(is,ib,ic) = HALF
      blocG36_FROM_blocGxG(is,ic,ib) = HALF
      CALL Write_Mat(blocG36_FROM_blocGxG(is,:,:),out_unitp,5)

      write(out_unitp,*) ' GxG => sym E4a'
      is=11
      blocG36_FROM_blocGxG(is,ia,ic) = HALF
      blocG36_FROM_blocGxG(is,ic,ia) =-HALF
      blocG36_FROM_blocGxG(is,ib,id) =-HALF
      blocG36_FROM_blocGxG(is,id,ib) = HALF
      CALL Write_Mat(blocG36_FROM_blocGxG(is,:,:),out_unitp,5)
      write(out_unitp,*) ' GxG => sym E4b'
      is=12
      blocG36_FROM_blocGxG(is,ia,id) = HALF
      blocG36_FROM_blocGxG(is,id,ia) =-HALF
      blocG36_FROM_blocGxG(is,ib,ic) = HALF
      blocG36_FROM_blocGxG(is,ic,ib) =-HALF
      CALL Write_Mat(blocG36_FROM_blocGxG(is,:,:),out_unitp,5)

      write(out_unitp,*) ' GxG => sym Ga'
      is=13
      blocG36_FROM_blocGxG(is,ib,ib) = ONE/sqrt(TWO)
      blocG36_FROM_blocGxG(is,ia,ia) =-ONE/sqrt(TWO)
      CALL Write_Mat(blocG36_FROM_blocGxG(is,:,:),out_unitp,5)
      write(out_unitp,*) ' GxG => sym Gb'
      is=14
      blocG36_FROM_blocGxG(is,ia,ib) = ONE/sqrt(TWO)
      blocG36_FROM_blocGxG(is,ib,ia) = ONE/sqrt(TWO)
      CALL Write_Mat(blocG36_FROM_blocGxG(is,:,:),out_unitp,5)
      write(out_unitp,*) ' GxG => sym Gc'
      is=15
      blocG36_FROM_blocGxG(is,id,id) = ONE/sqrt(TWO)
      blocG36_FROM_blocGxG(is,ic,ic) =-ONE/sqrt(TWO)
      CALL Write_Mat(blocG36_FROM_blocGxG(is,:,:),out_unitp,5)
      write(out_unitp,*) ' GxG => sym Gd'
      is=16
      blocG36_FROM_blocGxG(is,ic,id) = ONE/sqrt(TWO)
      blocG36_FROM_blocGxG(is,id,ic) = ONE/sqrt(TWO)
      CALL Write_Mat(blocG36_FROM_blocGxG(is,:,:),out_unitp,5)


stop
      symblocG36_FROM_blocGxG(1:16) = (/ (i,i=1,16) /) ! Ga,Gb,Gc,Gd



      end subroutine symetrization_G36_old
