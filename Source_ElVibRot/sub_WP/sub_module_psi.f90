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

      MODULE mod_psi
      USE mod_system
      USE mod_nDindex
      USE mod_Constant
      IMPLICIT NONE

      CONTAINS

      SUBROUTINE write1D2D_psi(psi,ana_psi)
      USE mod_system
      USE mod_dnSVM
      USE mod_psi_set_alloc
      USE mod_psi_B_TO_G
      IMPLICIT NONE

!----- variables for the WP ----------------------------------------
      TYPE (param_psi),     intent(inout)  :: psi
      TYPE (param_ana_psi), intent(inout)  :: ana_psi



!----- working variables ------------------------------------------
      real (kind=Rkind) :: val,val_mini,val_psi2,val_Rpsi,val_Cpsi

      integer           :: i_qa,i_qaie,iq,jq
      integer           :: ib,jb,i_bie


      TYPE (Type_dnVec) :: Tab1D_Qact(psi%BasisnD%nb_basis)
      TYPE (Type_dnVec) :: Tab1D_psi2(psi%BasisnD%nb_basis)
      TYPE (Type_dnVec) :: Tab1D_Rpsi(psi%BasisnD%nb_basis)
      TYPE (Type_dnVec) :: Tab1D_Cpsi(psi%BasisnD%nb_basis)

      TYPE (Type_dnMat) :: Tab2D_psi2(psi%BasisnD%nb_basis,psi%BasisnD%nb_basis)
      TYPE (Type_dnMat) :: Tab2D_Rpsi(psi%BasisnD%nb_basis,psi%BasisnD%nb_basis)
      TYPE (Type_dnMat) :: Tab2D_Cpsi(psi%BasisnD%nb_basis,psi%BasisnD%nb_basis)

      integer :: nDval0(psi%BasisnD%nb_basis)
      integer :: nDval(psi%BasisnD%nb_basis)
      integer :: nDvalib(psi%BasisnD%nb_basis)
      integer :: nDval0ib(psi%BasisnD%nb_basis)

!------ working variables ---------------------------------
      TYPE (param_file)              :: file_psi
      integer                        :: nio,nqi,nqj


!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='write1D2D_psi'
      logical,parameter :: debug = .FALSE.
      !logical,parameter :: debug = .TRUE.
!-----------------------------------------------------------

      IF (.NOT.ana_psi%psi1D_Q0 .AND. .NOT.ana_psi%psi2D_Q0) RETURN


      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
      END IF
!-----------------------------------------------------------

      IF (psi%nb_act1 /= psi%BasisnD%nb_basis) THEN
        write(out_unitp,*) 'WARNNING in',name_sub
        write(out_unitp,*) 'nb_act1 /= BasisnD%nb_basis',psi%nb_act1,psi%BasisnD%nb_basis
        RETURN
      END IF


      CALL sub_PsiBasisRep_TO_GridRep(psi)

      nDval0(:) = 0
      DO ib=1,psi%BasisnD%nb_basis
        nqi = get_nq_FROM_basis(psi%BasisnD%tab_Pbasis(ib)%Pbasis)
        CALL alloc_dnVec(Tab1D_Qact(ib),nqi)
        CALL alloc_dnVec(Tab1D_psi2(ib),nqi)
        CALL alloc_dnVec(Tab1D_Rpsi(ib),nqi)
        CALL alloc_dnVec(Tab1D_Cpsi(ib),nqi)

        Tab1D_Qact(ib)%d0(:) = psi%BasisnD%tab_Pbasis(ib)%Pbasis%x(1,:)

        ! set the table nDval0 from Qana and the Tab1D_Qact(ib)
        val_mini  = huge(ONE)
        DO iq=1,Tab1D_Qact(ib)%nb_var_vec
          val = abs(ana_psi%Qana(ib)-Tab1D_Qact(ib)%d0(iq))
          IF (val < val_mini) THEN
            val_mini = val
            nDval0(ib) = iq
          END IF
        END DO

        DO jb=1,psi%BasisnD%nb_basis
          nqj = get_nq_FROM_basis(psi%BasisnD%tab_Pbasis(jb)%Pbasis)
          CALL alloc_dnMat(Tab2D_psi2(jb,ib),nb_var_Matl=nqj,nb_var_Matc=nqi)
          CALL alloc_dnMat(Tab2D_Rpsi(jb,ib),nb_var_Matl=nqj,nb_var_Matc=nqi)
          CALL alloc_dnMat(Tab2D_Cpsi(jb,ib),nb_var_Matl=nqj,nb_var_Matc=nqi)
        END DO

      END DO
      !write(6,*) 'nDval0,Qana',nDval0,ana_psi%Qana

      i_bie = 1

      IF (allocated(psi%CvecG)) THEN

        DO i_qa=1,psi%nb_qa
          i_qaie = i_qa + (i_bie-1) * psi%nb_qa
          val_psi2 = abs(psi%CvecG(i_qaie))**2
          val_Rpsi = Real(psi%CvecG(i_qaie),kind=Rkind)
          val_Cpsi = aimag(psi%CvecG(i_qaie))

          CALL calc_nDindex(psi%BasisnD%nDindG,i_qa,nDval)

          DO ib=1,psi%BasisnD%nb_basis
            nDvalib(:)   = nDval(:)
            nDvalib(ib)  = 0

            nDval0ib(:)  = nDval0(:)
            nDval0ib(ib) = 0

            IF (compare_tab(nDval0ib,nDvalib)) THEN
              Tab1D_psi2(ib)%d0(nDval(ib)) = val_psi2
              Tab1D_Rpsi(ib)%d0(nDval(ib)) = val_Rpsi
              Tab1D_Cpsi(ib)%d0(nDval(ib)) = val_Cpsi
            END IF

            DO jb=ib+1,psi%BasisnD%nb_basis

              nDvalib(jb)  = 0
              nDval0ib(jb) = 0

              IF (compare_tab(nDval0ib,nDvalib)) THEN
                Tab2D_psi2(ib,jb)%d0(nDval(ib),nDval(jb)) = val_psi2
                Tab2D_Rpsi(ib,jb)%d0(nDval(ib),nDval(jb)) = val_Rpsi
                Tab2D_Cpsi(ib,jb)%d0(nDval(ib),nDval(jb)) = val_Cpsi
              END IF

            END DO

          END DO

        END DO
      ELSE IF (allocated(psi%RvecG)) THEN

        DO i_qa=1,psi%nb_qa
          i_qaie = i_qa + (i_bie-1) * psi%nb_qa
          val_psi2 = psi%RvecG(i_qaie)**2
          val_Rpsi = psi%RvecG(i_qaie)
          val_Cpsi = ZERO

          CALL calc_nDindex(psi%BasisnD%nDindG,i_qa,nDval)

          DO ib=1,psi%BasisnD%nb_basis
            nDvalib(:)   = nDval(:)
            nDvalib(ib)  = 0

            nDval0ib(:)  = nDval0(:)
            nDval0ib(ib) = 0

            IF (compare_tab(nDval0ib,nDvalib)) THEN
              Tab1D_psi2(ib)%d0(nDval(ib)) = val_psi2
              Tab1D_Rpsi(ib)%d0(nDval(ib)) = val_Rpsi
              Tab1D_Cpsi(ib)%d0(nDval(ib)) = val_Cpsi
            END IF

            DO jb=ib+1,psi%BasisnD%nb_basis

              nDvalib(jb)  = 0
              nDval0ib(jb) = 0

              IF (compare_tab(nDval0ib,nDvalib)) THEN
                Tab2D_psi2(ib,jb)%d0(nDval(ib),nDval(jb)) = val_psi2
                Tab2D_Rpsi(ib,jb)%d0(nDval(ib),nDval(jb)) = val_Rpsi
                Tab2D_Cpsi(ib,jb)%d0(nDval(ib),nDval(jb)) = val_Cpsi
              END IF

            END DO
          END DO


        END DO

      ELSE
        write(out_unitp,*) 'WARNING in ',name_sub
        write(out_unitp,*) ' Impossible to write CvecG or RvecG'
        RETURN
      END IF


      IF (ana_psi%psi1D_Q0) THEN

        DO ib=1,psi%BasisnD%nb_basis

          file_psi%name = 'psi1D_' // int_TO_char(ana_psi%iwp) // '-' // int_TO_char(ib)

          CALL file_open(file_psi,nio)

          DO iq=1,Tab1D_Qact(ib)%nb_var_vec
            write(nio,*) iq,Tab1D_Qact(ib)%d0(iq),Tab1D_Rpsi(ib)%d0(iq),  &
                             Tab1D_Cpsi(ib)%d0(iq),Tab1D_psi2(ib)%d0(iq)
          END DO
          close(nio)
        END DO
      END IF

      IF (ana_psi%psi2D_Q0) THEN

        DO ib=1,psi%BasisnD%nb_basis
        DO jb=ib+1,psi%BasisnD%nb_basis
          file_psi%name = 'psi2D_' // int_TO_char(ana_psi%iwp)
          file_psi%name = trim(file_psi%name) // '-' // int_TO_char(ib)
          file_psi%name = trim(file_psi%name) // '-' // int_TO_char(jb)
          CALL file_open(file_psi,nio)

          DO iq=1,Tab1D_Qact(ib)%nb_var_vec
          DO jq=1,Tab1D_Qact(jb)%nb_var_vec

            write(nio,*) iq,jq,Tab1D_Qact(ib)%d0(iq),Tab1D_Qact(jb)%d0(jq),&
              Tab2D_Rpsi(ib,jb)%d0(iq,jq),Tab2D_Cpsi(ib,jb)%d0(iq,jq),  &
                                             Tab2D_psi2(ib,jb)%d0(iq,jq)
          END DO
          write(nio,*)
          END DO
          close(nio)

        END DO
        END DO
      END IF


      DO ib=1,psi%BasisnD%nb_basis
        CALL dealloc_dnVec(Tab1D_Qact(ib))
        CALL dealloc_dnVec(Tab1D_psi2(ib))
        CALL dealloc_dnVec(Tab1D_Rpsi(ib))
        CALL dealloc_dnVec(Tab1D_Cpsi(ib))
        DO jb=1,psi%BasisnD%nb_basis
          CALL dealloc_dnMat(Tab2D_psi2(jb,ib))
          CALL dealloc_dnMat(Tab2D_Rpsi(jb,ib))
          CALL dealloc_dnMat(Tab2D_Cpsi(jb,ib))
        END DO
      END DO


!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'END ',name_sub
       END IF
!-----------------------------------------------------------


      END SUBROUTINE write1D2D_psi

!================================================================
!
!     sub_analyze_psi :
!       ....
!       Rho1D, Rho2D
!
!================================================================
      SUBROUTINE sub_analyze_psi(psi,info,ana_psi)
      USE mod_system
      USE mod_psi_set_alloc
      USE mod_psi_Op
      IMPLICIT NONE

      !----- variables for the WP -------------------------------
      TYPE (param_psi),     intent(inout) :: psi
      character (len=*),    intent(in)    :: info
      TYPE (param_ana_psi), intent(inout) :: ana_psi


      integer :: ibi,i_max_w,ih
      real (kind=Rkind) :: max_w,pop,Etemp,T
      character(len=:), allocatable     :: lformat
      TYPE(REAL_WU)                     :: RWU_Temp,RWU_E,RWU_DE
      real (kind=Rkind)                 :: E,DE

      real (kind=Rkind), allocatable :: tab_WeightChannels(:,:)


      !----- dynamic allocation memory ------------------------------------
      real (kind=Rkind), allocatable :: moy_Qba(:)

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_analyze_psi'
      integer :: err_mem,memory
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------

      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
      END IF

      IF (psi%ComOp%contrac_ba_ON_HAC) THEN
        ana_psi%AvQ = .FALSE.
      END IF


      CALL Channel_weight(ana_psi%tab_WeightChannels,psi,GridRep=.FALSE.,BasisRep=.TRUE.)

      i_max_w   = 1
      max_w     = ZERO
      DO ibi=1,psi%nb_bi
        IF (ana_psi%tab_WeightChannels(ibi,1) > max_w) THEN
          max_w = ana_psi%tab_WeightChannels(ibi,1)
          i_max_w = ibi
        END IF
      END DO

      RWU_Temp = REAL_WU(ana_psi%Temp,'°K','E')
      !Etemp   = RWU_Temp  ! Temperature convertion in Hartree
      Etemp    = convRWU_TO_R(RWU_Temp ,WorkingUnit=.TRUE.)

      pop = exp(-(ana_psi%Ene-ana_psi%ZPE)/Etemp) / ana_psi%Part_func

      RWU_E  = REAL_WU(ana_psi%Ene,'au','E')
      RWU_DE = REAL_WU(ana_psi%Ene-ana_psi%ZPE,'au','E')
      E  = convRWU_TO_R(RWU_E ,WorkingUnit=.FALSE.)
      DE = convRWU_TO_R(RWU_DE,WorkingUnit=.FALSE.)

      IF (ana_psi%AvQ) THEN
        CALL alloc_NParray(moy_Qba,(/2*Psi%BasisnD%ndim/),"moy_Qba",name_sub)
        CALL sub_Qmoy(psi,moy_Qba,ana_psi)

        IF (ana_psi%iwp < 10000 .AND. i_max_w < 10000) THEN
          lformat = String_TO_String('("lev: ",i4,i4,l3,' //            &
                                       int_TO_char(3+size(moy_Qba)) //  &
                           "(1x," // trim(adjustl(EneIO_format)) // "))")
        ELSE
          lformat = String_TO_String('("lev: ",i0,i0,l3,' //            &
                                       int_TO_char(3+size(moy_Qba)) //  &
                           "(1x," // trim(adjustl(EneIO_format)) // "))")
        END IF

        write(out_unitp,lformat) ana_psi%iwp,i_max_w,psi%convAvOp,  &
                                 E,DE,pop,moy_Qba(:)

        CALL dealloc_NParray(moy_Qba,"moy_Qba",name_sub)
      ELSE
        IF (ana_psi%iwp < 10000 .AND. i_max_w < 10000) THEN
          lformat = String_TO_String( '("lev: ",i4,i4,l3,3(1x,' //      &
                                   trim(adjustl(EneIO_format)) // '))' )
        ELSE
          lformat = String_TO_String( '("lev: ",i0,i0,l3,3(1x,' //      &
                                   trim(adjustl(EneIO_format)) // '))' )
        END IF

        write(out_unitp,lformat) ana_psi%iwp,i_max_w,psi%convAvOp,E,DE,pop

      END IF
      CALL flush_perso(out_unitp)

      IF (psi%nb_bi > 1) THEN

        lformat = String_TO_String( '("% HAC: ",' //                    &
                                int_TO_char(psi%nb_bi) // "(1x,f4.0) )" )

        write(out_unitp,lformat) (ana_psi%tab_WeightChannels(ih,1)*TEN**2,ih=1,psi%nb_bi)
      END IF

      deallocate(lformat)

      CALL calc_1Dweight(psi,ana_psi,20,real(ana_psi%iwp,kind=Rkind),info,.TRUE.)

      CALL Rho1D_Rho2D_psi(psi,ana_psi)

      CALL write1D2D_psi(psi,ana_psi)


      psi%GridRep = .FALSE.
      CALL alloc_psi(psi)

      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
      CALL flush_perso(out_unitp)

      END SUBROUTINE sub_analyze_psi

      SUBROUTINE sub_analyze_WP_forPropa(T,WP,nb_WP,field)
      USE mod_system
      USE mod_psi_set_alloc
      USE mod_psi_Op
      IMPLICIT NONE

      real (kind=Rkind) :: T      ! time
      real (kind=Rkind), intent(in), optional ::field(3)

      integer            :: nb_WP
      TYPE (param_psi)   :: WP(:)

      TYPE (param_ana_psi) :: ana_psi

!------ working parameters --------------------------------
      complex (kind=Rkind)   :: C12
      real (kind=Rkind)      :: Psi_norm2

      integer                       :: i
      integer                       :: max_ecri
      character(len=:), allocatable :: info

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='sub_analyze_WP_forPropa'
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
        write(out_unitp,*) 'nb_ba,nb_qa',WP(1)%nb_ba,WP(1)%nb_qa
        write(out_unitp,*) 'nb_bi',WP(1)%nb_bi
        write(out_unitp,*)

        DO i=1,nb_WP

          write(out_unitp,*) 'WP(i)%BasisRep',i
          CALL ecri_psi(T=ZERO,psi=WP(1),                               &
                        ecri_GridRep=.FALSE.,ecri_BasisRep=.TRUE.)
          write(out_unitp,*) 'WP(i)%GridRep',i
          CALL ecri_psi(T=ZERO,psi=WP(1),                               &
                        ecri_GridRep=.TRUE.,ecri_BasisRep=.FALSE.)
        END DO
       END IF
!-----------------------------------------------------------

        max_ecri = min(25,WP(1)%nb_tot)

        DO i=1,nb_WP

         CALL Channel_weight(ana_psi%tab_WeightChannels,WP(i),          &
                             GridRep=.FALSE.,BasisRep=.TRUE.)

          Psi_norm2 = sum(ana_psi%tab_WeightChannels)

          info = String_TO_String( '#WP ' // int_TO_char(i) )
          IF (present(field)) THEN
            write(out_unitp,21) info,T,real(WP(i)%CAvOp,kind=Rkind),    &
                           field(:),Psi_norm2,ana_psi%tab_WeightChannels
          ELSE
            write(out_unitp,21) info,T,real(WP(i)%CAvOp,kind=Rkind),    &
                                 Psi_norm2,ana_psi%tab_WeightChannels
          END IF
 21       format('norm^2-WP ',a,1x,f12.2,2(1x,f8.5),30(1x,f10.7))

          CALL calc_1Dweight(WP(i),ana_psi,max_ecri,T,info,.TRUE.)

          CALL psi_Qba_ie_psi(T,WP(i),ana_psi,info)

!         CALL calc_DM(WP(i),max_ecri,T,info,.TRUE.)
!         C12 = WP(i)%CvecB(1)*conjg(WP(i)%CvecB(2))
!         write(out_unitp,31) T,i,C12,abs(C12)
!31       format('C12',f12.1,1x,i2,3(1x,f12.6))

          !CALL calc_nDTk(WP(i),T)

          CALL write1D2D_psi(WP(i),ana_psi)

          deallocate(info)
        END DO


        CALL flush_perso(out_unitp)


      CALL dealloc_ana_psi(ana_psi)

!----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'END ',name_sub
       END IF
!----------------------------------------------------------
      END SUBROUTINE sub_analyze_WP_forPropa

!================================================================
!
!     calculation of <psi | Qact(j) | psi>
!
!================================================================
      SUBROUTINE sub_Qmoy(Psi,moy_Qba,ana_psi)
      USE mod_system
      USE mod_dnSVM
      USE mod_psi_set_alloc
      USE mod_psi_B_TO_G
      IMPLICIT NONE

      TYPE (param_psi)   :: Psi
      TYPE (param_ana_psi), intent(inout) :: ana_psi

!------ means value ---------------------------------------
      real (kind=Rkind)  ::    moy_Qba(2*Psi%nb_act1)


!------ working variables ---------------------------------
      integer           :: i,i_q,i_ie,iqie,dnErr
      real (kind=Rkind) :: psi2,WrhonD
      real (kind=Rkind) :: x(Psi%BasisnD%ndim)
      TYPE (Type_dnS)   :: dnt
      real (kind=Rkind) :: cte(20)

!----- for debuging --------------------------------------------------
      logical, parameter :: debug =.FALSE.
      !logical, parameter :: debug =.TRUE.
      !-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING sub_Qmoy'
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

      cte(:) = ZERO
      IF (allocated(ana_psi%Qtransfo_type)) CALL alloc_dnSVM(dnt,1,3)

      CALL sub_PsiBasisRep_TO_GridRep(Psi)

      moy_Qba(:)  = ZERO

      DO i_q=1,Psi%nb_qa
        psi2 = ZERO
        DO i_ie=1,Psi%nb_bi*Psi%nb_be

          iqie = (i_ie-1)*Psi%nb_qa + i_q

          IF (Psi%cplx) THEN
            psi2 = psi2 + abs(Psi%CvecG(iqie))**2
          ELSE
            psi2 = psi2 + abs(Psi%RvecG(iqie))**2
          END IF
        END DO

        !- calculation of x -------------------------------
        CALL Rec_x(x,psi%BasisnD,i_q)

        IF (allocated(ana_psi%Qtransfo_type)) THEN
          DO i=1,size(x)
            CALL sub_dntf(ana_psi%Qtransfo_type(i),dnt,x(i),cte,dnErr)
            IF (dnErr /= 0) THEN
              write(out_unitp,*) ' ERROR in sub_Qmoy'
              write(out_unitp,*) '   ERROR in the sub_dntf call for the coordinates, iQdyn:',psi%BasisnD%iQdyn(i)
              STOP 'ERROR in sub_dntf called from sub_Qmoy'
            END IF
            x(i) = dnt%d0
          END DO
        END IF
        !write(out_unitp,*) 'x',x

        !- calculation of WrhonD ------------------------------
        WrhonD = Rec_WrhonD(Psi%BasisnD,i_q)
        ! write(out_unitp,*) i,'WrhonD',WrhonD

        psi2 = psi2 * WrhonD

        moy_Qba(1:Psi%nb_act1) = moy_Qba(1:Psi%nb_act1) +             &
                                                         psi2 * x(:)
        moy_Qba(1+Psi%nb_act1:2*Psi%nb_act1) =                         &
                              moy_Qba(1+Psi%nb_act1:2*Psi%nb_act1) +   &
                                                      psi2 * x(:)**2
      END DO

      moy_Qba(1+Psi%nb_act1:2*Psi%nb_act1) =                           &
                              moy_Qba(1+Psi%nb_act1:2*Psi%nb_act1) -   &
                                            moy_Qba(1:Psi%nb_act1)**2

      Psi%GridRep= .FALSE.
      CALL alloc_psi(Psi)
      IF (allocated(ana_psi%Qtransfo_type)) CALL dealloc_dnSVM(dnt)

!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END sub_Qmoy'
      END IF
!----------------------------------------------------------

      END SUBROUTINE sub_Qmoy

!================================================================
!
!     means of Qact1(i) psi  : <psi|Qact(i)|psi>
!
!================================================================
      SUBROUTINE psi_Qba_ie_psi(T,psi,ana_psi,info)
      USE mod_system
      USE mod_psi_set_alloc
      USE mod_psi_Op
      USE mod_psi_B_TO_G
      IMPLICIT NONE

!----- variables for the WP ----------------------------------------
      TYPE (param_psi),     intent(inout) :: psi
      TYPE (param_ana_psi), intent(in)    :: ana_psi

      character (len=*) :: info

      real (kind=Rkind) :: T


      real (kind=Rkind) :: Qmean(psi%nb_act1)
      real (kind=Rkind) :: Qmean_ie(psi%nb_act1,psi%nb_bi,psi%nb_be)
      real (kind=Rkind) :: x(Psi%BasisnD%ndim)

!------ working variables ---------------------------------
      integer       :: i_qa,i_qaie
      integer       :: i_be,i_bi,i_ba,i_baie
      integer       :: ii_baie,if_baie
      integer       :: i
      real (kind=Rkind) :: WrhonD,temp

      logical       :: psiN,normeGridRep,normeBasisRep

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub = 'psi_Qba_ie_psi'
      logical,parameter :: debug = .FALSE.
      !logical,parameter :: debug = .TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'psi'
        CALL ecri_psi(Psi=psi)
      END IF
!-----------------------------------------------------------

      IF (psi%nb_baie > psi%nb_tot) RETURN


      IF (.NOT. allocated(ana_psi%tab_WeightChannels)) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'ana_psi%tab_WeightChannels is not allocated!!'
        write(out_unitp,*) ' It should be done in sub_analyze_WP_forPropa'
        STOP
      END IF

      Qmean(:) = ZERO
      Qmean_ie(:,:,:) = ZERO


      psi%GridRep = .TRUE.
      CALL alloc_psi(psi)
      CALL sub_PsiBasisRep_TO_GridRep(psi)

      DO i_qa=1,psi%nb_qa

        !- calculation of WrhonD ------------------------------
        WrhonD = Rec_WrhonD(psi%BasisnD,i_qa)

        !- calculation of x -------------------------------
        CALL Rec_x(x,psi%BasisnD,i_qa)

        DO i_be=1,psi%nb_be
        DO i_bi=1,psi%nb_bi
          i_qaie = i_qa + ( (i_bi-1)+(i_be-1)*psi%nb_bi ) * psi%nb_qa

          IF (psi%cplx) THEN
            temp = abs( psi%CvecG(i_qaie) )
          ELSE
            temp = psi%RvecG(i_qaie)
          END IF
          temp = WrhonD * temp**2

          Qmean(:) = Qmean(:) + x(:) * temp
          Qmean_ie(:,i_bi,i_be) = Qmean_ie(:,i_bi,i_be) + x(:) * temp
        END DO
        END DO
      END DO


      Qmean(:) = Qmean(:) / sum(ana_psi%tab_WeightChannels)

      DO i_be=1,psi%nb_be
      DO i_bi=1,psi%nb_bi
        IF (ana_psi%tab_WeightChannels(i_bi,i_be) > ONETENTH**7) THEN
          Qmean_ie(:,i_bi,i_be) = Qmean_ie(:,i_bi,i_be) /               &
                                  ana_psi%tab_WeightChannels(i_bi,i_be)
        ELSE
          Qmean_ie(:,i_bi,i_be) = ZERO
        END IF
      END DO
      END DO

      DO i=1,psi%nb_act1
          write(out_unitp,11) 'T iQbasis Qmean_ie ',info,T,i,Qmean_ie(i,:,:)
          write(out_unitp,11) 'T iQbasis Qmean    ',info,T,i,Qmean(i)
      END DO
 11   format(2a,' ',f12.4,' ',i4,' ',100(' ',f6.3))
      CALL flush_perso(out_unitp)

      psi%GridRep = .FALSE.
      CALL alloc_psi(psi)

!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'Qmean',Qmean
        DO i=1,psi%nb_act1
          write(out_unitp,*) 'Qmean_ie',i,Qmean_ie(i,:,:)
        END DO
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------


      END SUBROUTINE psi_Qba_ie_psi

      SUBROUTINE Rho1D_Rho2D_psi(psi,ana_psi)
      USE mod_system
      USE mod_psi_set_alloc
      USE mod_psi_B_TO_G
      IMPLICIT NONE
!----- variables for the WP ----------------------------------------
      TYPE (param_psi),     intent(inout) :: psi
      TYPE (param_ana_psi), intent(inout) :: ana_psi


!------ working variables ---------------------------------
      TYPE (param_file)              :: file_Rho
      integer                        :: i_basis_act1,j_basis_act1,nioRho
      real (kind=Rkind), allocatable :: rho2D(:,:)
      real (kind=Rkind), allocatable :: rho1D(:)


      integer               :: i_qa,i_qaie
      integer               :: i_be,i_bi,i_ba
      integer               :: i_baie,f_baie
      real (kind=Rkind)     :: WrhonD,Weight
      complex (kind=Rkind)  :: temp
      real (kind=Rkind)     :: Rtemp,Qix

      integer :: iq_act1,jq_act1,i,j,iq,ix,ix_TO_iQdyn
      integer :: nDval(psi%BasisnD%nDindG%ndim)

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Rho1D_Rho2D_psi'
      logical,parameter :: debug = .FALSE.
!     logical,parameter :: debug = .TRUE.
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'psi'
        CALL ecri_psi(psi=psi)
      END IF
!-----------------------------------------------------------
      IF (ana_psi%ana .AND. (ana_psi%Rho1D .OR. ana_psi%Rho2D)) THEN
        CALL sub_PsiBasisRep_TO_GridRep(psi)


        IF (ana_psi%Rho1D) THEN
          !- loop on coordinates ----------------------------------
          DO i_basis_act1=1,psi%BasisnD%nb_basis

            file_Rho%name = 'Rho1D_' // int_TO_char(ana_psi%iwp) // &
                                    '-' // int_TO_char(i_basis_act1)

            CALL file_open(file_Rho,nioRho)

            CALL alloc_NParray(rho1D,                                   &
                           (/psi%BasisnD%nDindG%nDsize(i_basis_act1)/), &
                              'rho1D',name_sub)
            rho1D(:) = ZERO

            DO i_qa=1,psi%nb_qa

              CALL calc_nDindex(psi%BasisnD%nDindG,i_qa,nDval)

              !- calculation of WrhonD ------------------------------
              WrhonD   = ONE
              DO i=1,psi%BasisnD%nb_basis
                iq  = nDval(i)

                IF (allocated(ana_psi%Weight_Rho) .AND.            &
                              allocated(ana_psi%Qana_weight) ) THEN

                  DO ix=1,psi%BasisnD%tab_Pbasis(i)%Pbasis%ndim
                    ix_TO_iQdyn = psi%BasisnD%tab_Pbasis(i)%Pbasis%iQdyn(ix)
                    IF (ana_psi%Weight_Rho(ix_TO_iQdyn) == 1) THEN
                      Qix = psi%BasisnD%tab_Pbasis(i)%Pbasis%x(ix,iq)
                      IF (Qix > ana_psi%Qana_weight(ix_TO_iQdyn)) THEN
                        WrhonD   = ZERO
                      END IF
                    ELSE IF (ana_psi%Weight_Rho(ix_TO_iQdyn) == -1) THEN
                      Qix = psi%BasisnD%tab_Pbasis(i)%Pbasis%x(ix,iq)
                      IF (Qix < ana_psi%Qana_weight(ix_TO_iQdyn)) THEN
                        WrhonD   = ZERO
                      END IF
                    END IF
                  END DO

                END IF

                IF (WrhonD == ZERO) CYCLE

                IF (i == i_basis_act1) THEN
                  WrhonD = WrhonD * Rec_rhonD(psi%BasisnD%tab_Pbasis(i)%Pbasis,iq)
                ELSE
                  WrhonD = WrhonD * Rec_WrhonD(psi%BasisnD%tab_Pbasis(i)%Pbasis,iq)
                END IF
              END DO

              DO i_be=1,psi%nb_be
              DO i_bi=1,psi%nb_bi
                i_qaie=i_qa + ( (i_bi-1)+(i_be-1)*psi%nb_bi ) * psi%nb_qa

                IF (psi%cplx) THEN
                  temp = conjg(psi%CvecG(i_qaie))*psi%CvecG(i_qaie)
                ELSE
                  temp = cmplx(psi%RvecG(i_qaie)*psi%RvecG(i_qaie),ZERO,kind=Rkind)
                END IF

                iq_act1 = nDval(i_basis_act1)
                rho1D(iq_act1) = rho1D(iq_act1) + WrhonD * real(temp,kind=Rkind)
              END DO
              END DO
            END DO

            write(out_unitp,*) 'rho1D of ',ana_psi%iwp,i_basis_act1
            DO i=1,psi%BasisnD%nDindG%nDsize(i_basis_act1)
             write(nioRho,*) ana_psi%iwp,i,                         &
                psi%BasisnD%tab_Pbasis(i_basis_act1)%Pbasis%x(:,i),rho1D(i)
            END DO

            CALL dealloc_NParray(rho1D,'rho1D',name_sub)
            close(nioRho)
          END DO
        END IF

        IF (ana_psi%Rho2D) THEN

          !- loop on coordinates ----------------------------------
          DO i_basis_act1=1,psi%BasisnD%nb_basis
          DO j_basis_act1=i_basis_act1+1,psi%BasisnD%nb_basis
            file_Rho%name = 'Rho2D_' // int_TO_char(ana_psi%iwp)
            file_Rho%name = trim(file_Rho%name) // '-' // int_TO_char(i_basis_act1)
            file_Rho%name = trim(file_Rho%name) // '-' // int_TO_char(j_basis_act1)
            CALL file_open(file_Rho,nioRho)

            CALL alloc_NParray(rho2D,                                   &
                            (/psi%BasisnD%nDindG%nDsize(i_basis_act1),  &
                              psi%BasisnD%nDindG%nDsize(j_basis_act1)/),&
                              'rho2D',name_sub)
            rho2D(:,:) = ZERO

            DO i_qa=1,psi%nb_qa

              CALL calc_nDindex(psi%BasisnD%nDindG,i_qa,nDval)

              !- calculation of WrhonD ------------------------------
              WrhonD   = ONE
              DO i=1,psi%BasisnD%nb_basis
                iq  = nDval(i)

                IF (allocated(ana_psi%Weight_Rho) .AND.            &
                              allocated(ana_psi%Qana_weight) ) THEN
                  DO ix=1,psi%BasisnD%tab_Pbasis(i)%Pbasis%ndim
                    ix_TO_iQdyn = psi%BasisnD%tab_Pbasis(i)%Pbasis%iQdyn(ix)

                    IF (ana_psi%Weight_Rho(ix_TO_iQdyn) == 1) THEN
                      Qix = psi%BasisnD%tab_Pbasis(i)%Pbasis%x(ix,iq)
                      IF (Qix > ana_psi%Qana_weight(ix_TO_iQdyn)) THEN
                        WrhonD   = ZERO
                      END IF
                    ELSE IF (ana_psi%Weight_Rho(ix_TO_iQdyn) == -1) THEN
                      Qix = psi%BasisnD%tab_Pbasis(i)%Pbasis%x(ix,iq)
                      IF (Qix < ana_psi%Qana_weight(ix_TO_iQdyn)) THEN
                        WrhonD   = ZERO
                      END IF
                    END IF
                  END DO

                  IF (WrhonD == ZERO) CYCLE
                END IF


                IF (i == i_basis_act1 .OR. i == j_basis_act1) THEN
                  WrhonD = WrhonD * Rec_rhonD(psi%BasisnD%tab_Pbasis(i)%Pbasis,iq)
                ELSE
                  WrhonD = WrhonD * Rec_WrhonD(psi%BasisnD%tab_Pbasis(i)%Pbasis,iq)
                END IF
              END DO

              DO i_be=1,psi%nb_be
              DO i_bi=1,psi%nb_bi
                i_qaie=i_qa + ( (i_bi-1)+(i_be-1)*psi%nb_bi ) * psi%nb_qa

                IF (psi%cplx) THEN
                  temp = conjg(psi%CvecG(i_qaie))*psi%CvecG(i_qaie)
                ELSE
                  temp = cmplx(psi%RvecG(i_qaie)*psi%RvecG(i_qaie),ZERO,kind=Rkind)
                END IF

                iq_act1 = nDval(i_basis_act1)
                jq_act1 = nDval(j_basis_act1)
                rho2D(iq_act1,jq_act1) = rho2D(iq_act1,jq_act1) +       &
                                       WrhonD * real(temp,kind=Rkind)
              END DO
              END DO
            END DO

            write(out_unitp,*) 'rho2D of ',ana_psi%iwp,i_basis_act1,j_basis_act1
            DO i=1,psi%BasisnD%nDindG%nDsize(i_basis_act1)
            DO j=1,psi%BasisnD%nDindG%nDsize(j_basis_act1)
               write(nioRho,*) ana_psi%iwp,i,j,                     &
                 psi%BasisnD%tab_Pbasis(i_basis_act1)%Pbasis%x(:,i),    &
                 psi%BasisnD%tab_Pbasis(j_basis_act1)%Pbasis%x(:,j),    &
                 rho2D(i,j)
            END DO
            write(nioRho,*)
            END DO

            CALL dealloc_NParray(rho2D,'rho2D',name_sub)
            close(nioRho)
          END DO
          END DO
        END IF

        ! deallocate psi on the grid
        psi%GridRep = .FALSE.
        CALL alloc_psi(psi)


      END IF

!----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!----------------------------------------------------------


      END SUBROUTINE Rho1D_Rho2D_psi


      END MODULE mod_psi

