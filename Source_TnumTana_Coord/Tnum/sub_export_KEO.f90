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
MODULE mod_export_KEO
  USE mod_system
  use mod_dnSVM,    only: type_dnmat, alloc_dnsvm, dealloc_dnsvm,       &
                          alloc_array, set_zero_to_dnsvm, dealloc_array
  use mod_Tnum,     only: CoordType, tnum
  use mod_dnGG_dng, only: get_dng_dngg
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: export3_MCTDH_T,export_Taylor_dnG

  CONTAINS

!===========================================================
!     Export T for MCTDH
!===========================================================
      SUBROUTINE export3_MCTDH_T(Qact,para_Tnum,mole)
      IMPLICIT NONE

!     - for Tnum -------------------------------------------
      TYPE (CoordType) :: mole
      TYPE (Tnum)      :: para_Tnum
      real (kind=Rkind), intent(inout) :: Qact(mole%nb_var)

!     - G,g ... --------------------------------------------
      TYPE(Type_dnMat) :: dnGG
      integer :: ndimG



!     - divers ------------------------------------------
      integer           :: i
      real (kind=Rkind) :: epsi_MCTDH
      logical :: label,grid1D,perio,EVRT
      logical :: periodic(mole%nb_var)

      NAMELIST /MCTDH/epsi_MCTDH,grid1D,perio,EVRT

      epsi_MCTDH  = ONETENTH**8
      grid1D      = .FALSE.
      perio       = .FALSE.
      periodic(:) = .FALSE.
      EVRT        = .FALSE.
      read(in_unitp,MCTDH,end=999,err=999)
      IF (perio) read(in_unitp,*) periodic(:)

 999  CONTINUE

      ndimG = mole%ndimG
      CALL alloc_dnSVM(dnGG,ndimG,ndimG,mole%nb_act,2)

      para_Tnum%WriteT  = .FALSE.
      CALL get_dng_dnGG(Qact,para_Tnum,mole,dnGG=dnGG,nderiv=2)


      IF (EVRT) CALL export_Taylor_dnG(dnGG,Qact,epsi_MCTDH,option=1)

      write(out_unitp,'(a)') '# =============================================='
      write(out_unitp,'(a)') '# =============================================='
      write(out_unitp,'(a)') '# === Taylor expansion of KEO for MCTDH ========'
      write(out_unitp,'(a)') '# === !! The coordinates are the Qi ============'
      write(out_unitp,'(a)') '# === But we define NQi for the KEO ============'
      write(out_unitp,'(a)') '# === NQi=q[p] : NQi = ( Qi - p ) =============='
      write(out_unitp,'(a)') '# = or NQi=sin[1.,p] : NQi = sin( 1.(Qi - p)) =='
      write(out_unitp,'(a)') '# =============================================='
      write(out_unitp,'(a)') '# =============================================='

      write(out_unitp,'(a)') '-----------------------------------------------'
      write(out_unitp,*) 'LABELS-SECTION'
      write(out_unitp,'(a)') '-----------------------------------------------'
      DO i=1,mole%nb_act
        IF (periodic(i)) THEN
          write(out_unitp,*) 'NQ',int_TO_char(i),'=sin[1., ',           &
                                        real_TO_char_MCTDH(Qact(i)),' ]'
       ELSE
          write(out_unitp,*) 'NQ',int_TO_char(i),'=q[ ',                &
                                        real_TO_char_MCTDH(Qact(i)),' ]'
        END IF
      END DO
      IF (grid1D) THEN
        label = .TRUE.
        STOP 'grid1D does not work yet'
        CALL export3_d0G_grid1D(Qact,para_Tnum,mole,dnGG,label,epsi_MCTDH)
      ENDIF
      write(out_unitp,'(a)') '-----------------------------------------------'
      write(out_unitp,*) 'END-LABELS-SECTION'
      write(out_unitp,'(a)') '-----------------------------------------------'

      write(out_unitp,*)
      write(out_unitp,*)
      write(out_unitp,'(a)') '-----------------------------------------------'
      write(out_unitp,'(a)') 'HAMILTONIAN-SECTION'
      write(out_unitp,'(a)',advance='no') 'modes  '
      DO i=1,mole%nb_act
        write(out_unitp,'(2a)',advance='no') ' | ',trim(mole%name_Qact(i))
      END DO
      write(out_unitp,'(a)',advance='yes')
      write(out_unitp,'(a)') '-----------------------------------------------'

      IF (grid1D) THEN
        label = .FALSE.
        STOP 'grid1D does not work yet'
        CALL export3_d0G_grid1D(Qact,para_Tnum,mole,dnGG,label,epsi_MCTDH)
      END IF

      CALL export3_MCTDH_dnG(dnGG,grid1D,epsi_MCTDH)

      write(out_unitp,'(a)') '-----------------------------------------------'
      write(out_unitp,'(a)') 'END-HAMILTONIAN-SECTION'
      write(out_unitp,'(a)') '-----------------------------------------------'

      CALL dealloc_dnSVM(dnGG)

      write(out_unitp,'(a)') '# =============================================='
      write(out_unitp,'(a)') '# =============================================='
      write(out_unitp,'(a)') '# === END: Taylor expansion of T for MCTDH ====='
      write(out_unitp,'(a)') '# =============================================='
      write(out_unitp,'(a)') '# =============================================='

      end subroutine export3_MCTDH_T
!===========================================================
!     Export G in 1D grid for MCTDH
!===========================================================
     SUBROUTINE export3_d0G_grid1D(Qact,para_Tnum,mole,dnGG,label,epsi_MCTDH)
      IMPLICIT NONE

!----- for the CoordType and Tnum --------------------------------------
      TYPE (CoordType) :: mole
      TYPE (Tnum)      :: para_Tnum

      real (kind=Rkind), intent(inout) :: Qact(mole%nb_var)


      TYPE(Type_dnMat) :: dnGG

      TYPE(Type_dnMat), pointer :: dnGG_grid(:)

      !character (len=Name_len) :: name_i,name_j,name_iQact
      !character (len=Name_len) :: name_Opi,name_Opj
!     ------------------------------------------------------

!     ------------------------------------------------------
      type (param_file) :: file_Ggrid
      character (len=Name_len) ::name_file
      integer :: no

!     - divers ------------------------------------------
      integer :: i,j,iq,iQact,iiQact
      integer :: nb_nq
      real (kind=Rkind), pointer :: q(:)
      real (kind=Rkind) :: Qact_ref(mole%nb_var)
      real (kind=Rkind) :: val
      real (kind=Rkind) :: epsi_MCTDH
      logical :: label
!     - function ----------------------------------------


!===========================================================
!===========================================================
      nullify(dnGG_grid)
      nullify(q)

      Qact_ref(:) = Qact(:)

!-------------------------------------------------
!     - calculation of G,g and d1G,d1g.... -------
!     --------------------------------------------
      IF (.NOT. label) THEN
        write(out_unitp,'(a)') '-----------------------------------------------'
        write(out_unitp,*)
        write(out_unitp,'(a)') '-----------------------------------------------'
        write(out_unitp,'(a)') '# Zero order part: -1/2*G^ij(Qref)'
        write(out_unitp,'(a)') '-----------------------------------------------'
        DO i=1,mole%nb_act1
          IF (abs(dnGG%d0(i,i)) < epsi_MCTDH) cycle

          write(out_unitp,*) -HALF*dnGG%d0(i,i),' |' // int_TO_char(i),' dq^2'

        END DO
        DO i=1,mole%nb_act1
        DO j=i+1,mole%nb_act1
          IF (abs(dnGG%d0(i,j)) < epsi_MCTDH) cycle

          write(out_unitp,*) -dnGG%d0(i,j),' |' // int_TO_char(i),' dq',&
                              ' |' // int_TO_char(j),' dq'

        END DO
        END DO
        write(out_unitp,'(a)') '-----------------------------------------------'
      END IF
!==============================================
!     loop on active coordinates
!==============================================
      DO iQact=1,mole%nb_act1
        read(in_unitp,*) nb_nq,iiQact
        IF (iQact /= iiQact) THEN
          write(out_unitp,*) 'ERROR iQact /= iiQact',iQact,iiQact
          STOP
        END IF

        CALL alloc_array(q,[nb_nq],"q","export3_d0G_grid1D")

        CALL alloc_array(dnGG_grid,[nb_nq],                           &
                        "dnGG_grid","export3_d0G_grid1D")
        DO iq=1,nb_nq
          CALL alloc_dnSVM(dnGG_grid(iq),                               &
                           mole%ndimG,mole%ndimG,mole%nb_act,0)
          CALL Set_ZERO_TO_dnSVM(dnGG_grid(iq))
        END DO

        DO iq=1,nb_nq
          read(in_unitp,*) iiQact,q(iq)
        END DO
        read(in_unitp,*)
        Qact(:) = Qact_ref(:)

!       - grid calculation
        DO iq=1,nb_nq

          Qact(iQact) = q(iq)


          CALL get_dng_dnGG(Qact,para_Tnum,mole,dnGG=dnGG_grid(iq),nderiv=0)

          dnGG_grid(iq)%d0(:,:) = dnGG_grid(iq)%d0(:,:) - dnGG%d0(:,:)

        END DO




        IF (label) THEN
        write(out_unitp,'(a)') '-------------------------------------------------'
          write(out_unitp,'(2a)') '# G : 1D-grid along ',                       &
                            trim(adjustl(mole%name_Qact(iQact)))
        write(out_unitp,'(a)') '-------------------------------------------------'
          DO i=1,mole%nb_act1
          DO j=i,mole%nb_act1

            val = ZERO
            DO iq=1,nb_nq
              val = val + abs(dnGG_grid(iq)%d0(i,j))
            END DO
            val = val / real(nb_nq,kind=Rkind)

            IF (val > epsi_MCTDH*ONETENTH) THEN
!           write(out_unitp,'(a,2(i4,x),e10.4)') '# crit: ',i,j,val

              name_file = 'G' // trim(adjustl(mole%name_Qact(iQact))) //&
                      '_' // int_TO_char(i) // '_' // int_TO_char(j)
              write(out_unitp,*) trim(adjustl(name_file)),              &
                      ' = read1d{',trim(adjustl(name_file)),' ascii}'

              file_Ggrid%name = name_file
              CALL file_open(file_Ggrid,no)
              DO iq=1,nb_nq
                write(no,*) dnGG_grid(iq)%d0(i,j)
              END DO
              close(no)
            END IF
          END DO
          END DO
        ELSE


          write(out_unitp,'(a)') '-------------------------------------------------'
          write(out_unitp,'(2a)') '# G : 1D-grid along ',                       &
                            trim(adjustl(mole%name_Qact(iQact)))
          write(out_unitp,'(a)') '-------------------------------------------------'
          DO i=1,mole%nb_act1
          DO j=i,mole%nb_act1

            val = ZERO
            DO iq=1,nb_nq
              val = val + abs(dnGG_grid(iq)%d0(i,j))
            END DO
            val = val / real(nb_nq,kind=Rkind)

            IF (val > epsi_MCTDH*ONETENTH) THEN

              name_file = 'G' // trim(adjustl(mole%name_Qact(iQact))) //&
                          '_' // int_TO_char(i) // '_' // int_TO_char(j)

              IF (i == j) THEN
                IF (i == iQact) THEN
                  write(out_unitp,*) '-0.5 |',int_TO_char(i),' dq*',    &
                                          trim(adjustl(name_file)),'*dq'
                ELSE
                  write(out_unitp,*) '-0.5 |',int_TO_char(i),' dq^2 |',   &
                         int_TO_char(iQact),' ',trim(adjustl(name_file))
                END IF
              ELSE ! i/=j
                IF (i == iQact .AND. j/=iQact) THEN
                  write(out_unitp,*) '-0.5 |',int_TO_char(i),           &
                   '  dq*',trim(adjustl(name_file)),' |',int_TO_char(j),' dq'

                  write(out_unitp,*) '-0.5 |',int_TO_char(i),           &
                    ' ',trim(adjustl(name_file)),'*dq |',int_TO_char(j),' dq'

                ELSE IF (i /= iQact .AND. j==iQact) THEN
                  write(out_unitp,*) '-0.5 |',int_TO_char(j),           &
                    ' dq*',trim(adjustl(name_file)),' |',int_TO_char(i),' dq'

                  write(out_unitp,*) '-0.5 |',int_TO_char(j),           &
                    ' ',trim(adjustl(name_file)),'*dq |',int_TO_char(i),' dq'

                ELSE IF (i /= iQact .AND. j/=iQact) THEN
                  write(out_unitp,*) '-1.0 |',int_TO_char(iQact),' ',   &
                    trim(adjustl(name_file)),                           &
                    ' |' // int_TO_char(i),' dq',' |' // int_TO_char(j),' dq'
                ELSE
                  write(out_unitp,*) 'Should never append !',i,j,iQact
                  STOP
                END IF
              END IF
            END IF
        END DO
        END DO
      END IF


      DO iq=1,nb_nq
        CALL dealloc_dnSVM(dnGG_grid(iq))
      END DO

      CALL dealloc_array(q,"q","export3_d0G_grid1D")

      CALL dealloc_array(dnGG_grid,"dnGG_grid","export3_d0G_grid1D")
      END DO
!==============================================
!     END loop on active coordinates
!==============================================

      end subroutine export3_d0G_grid1D
!================================================================
!      Export second order Taylor expansion of G for MCTDH
!================================================================
      SUBROUTINE export3_MCTDH_dnG(dnGG,grid1D,epsi_MCTDH)
      IMPLICIT NONE

!     - G,g ... --------------------------------------------
      TYPE(Type_dnMat) :: dnGG
      real (kind=Rkind) :: epsi_MCTDH
      logical :: grid1D

!     - divers ------------------------------------------
      integer :: i,j,k,l,n0,n1,n2


      IF (.NOT. grid1D) THEN

        write(out_unitp,'(a)') '-----------------------------------------------'
        write(out_unitp,'(a)') '# Zero order part: -1/2*G^ij(Qref)'
        write(out_unitp,'(a)') '-----------------------------------------------'
        n0 = 0
        k=0
        l=0
        DO i=1,dnGG%nb_var_deriv
        DO j=1,dnGG%nb_var_deriv

          IF (abs(dnGG%d0(i,j)) < epsi_MCTDH) CYCLE

          !CALL T_Operator_MCTDH(name_Op,i,j,k,l,dnGG%nb_var_deriv)

          write(out_unitp,*) real_TO_char_MCTDH( -HALF*dnGG%d0(i,j) ),  &
                         get_T_Operator_MCTDH(i,j,k,l,dnGG%nb_var_deriv)
          n0 = n0+1

        END DO
        END DO

        write(out_unitp,'(a)') '-----------------------------------------------'
        write(out_unitp,'(a)') '# First order part:'
        write(out_unitp,'(a)') '#         -1/2*dG^ij/dNQk d./dQi * NQk * d./dQj'
        write(out_unitp,'(a)') '-----------------------------------------------'
        n1 = 0
        l=0
        DO i=1,dnGG%nb_var_deriv
        DO j=1,dnGG%nb_var_deriv
        DO k=1,dnGG%nb_var_deriv

          IF (abs(dnGG%d1(i,j,k)) < epsi_MCTDH) CYCLE

          !CALL T_Operator_MCTDH(name_Op,i,j,k,l,dnGG%nb_var_deriv)

          write(out_unitp,*) real_TO_char_MCTDH( -HALF*dnGG%d1(i,j,k) ),&
                        get_T_Operator_MCTDH(i,j,k,l,dnGG%nb_var_deriv)

          n1 = n1+1

        END DO
        END DO
        END DO
      END IF
      write(out_unitp,*)
      write(out_unitp,'(a)') '-----------------------------------------------'
      write(out_unitp,'(a)') '# Second order part:'
      write(out_unitp,'(a)') '# -1/4*d^2G^ij/dNQk^2   d./dQi * NQk*NQl * d./dQj'
      write(out_unitp,'(a)') '-----------------------------------------------'

      n2=0
      DO i=1,dnGG%nb_var_deriv
      DO j=1,dnGG%nb_var_deriv
      DO k=1,dnGG%nb_var_deriv
      DO l=1,dnGG%nb_var_deriv
        IF (abs(dnGG%d2(i,j,k,l)) < epsi_MCTDH) CYCLE
        IF (k == l .AND. grid1D) CYCLE

        !CALL T_Operator_MCTDH(name_Op,i,j,k,l,dnGG%nb_var_deriv)

        write(out_unitp,*) real_TO_char_MCTDH( -QUARTER*dnGG%d2(i,j,k,l) ),&
                           get_T_Operator_MCTDH(i,j,k,l,dnGG%nb_var_deriv)

        n2 = n2+1

      END DO
      END DO
      END DO
      END DO

      end subroutine export3_MCTDH_dnG
!================================================================
!      Export second order Taylor expansion of G for EVR
!================================================================
      SUBROUTINE export_Taylor_dnG(dnGG,Qact,epsi_G,file_name,option)
      IMPLICIT NONE

!     - G,g ... --------------------------------------------
      TYPE(Type_dnMat),  intent(in)           :: dnGG
      real (kind=Rkind), intent(in)           :: Qact(dnGG%nb_var_deriv)
      real (kind=Rkind), intent(in)           :: epsi_G
      character(len=*),  intent(in), optional :: file_name
      integer,           intent(in), optional :: option

!     - local variables -------------------------------------
      integer :: i,j,k,l,n0,n1,n2,nb_act,nio,option_loc


!     -----------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='export_Taylor_dnG'
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'Qact',Qact
        write(out_unitp,*) 'epsi_G',epsi_G

        IF (present(file_name)) THEN 
          write(out_unitp,*) 'file_name',file_name
        ELSE
          write(out_unitp,*) 'file_name is not present'
        END IF

        IF (present(option)) THEN 
          write(out_unitp,*) 'option',option
        ELSE
          write(out_unitp,*) 'option is not present'
        END IF
      END IF
      flush(out_unitp)
!     -----------------------------------------------------------------

      nb_act = dnGG%nb_var_deriv
      IF (size(Qact) < dnGG%nb_var_deriv) THEN
         write(out_unitp,*) ' ERROR in export3_dnG'
         write(out_unitp,*) ' size(Qact) < dnGG%nb_var_deriv'
         write(out_unitp,*) ' size(Qact),dnGG%nb_var_deriv',            &
                              size(Qact),dnGG%nb_var_deriv
        STOP 'ERROR in export_Taylor_dnG: wrong Qact size'
      END IF
      IF (debug) write(out_unitp,*) ' nb_act',nb_act

      IF (present(option)) THEN
        option_loc = option
      ELSE
        option_loc = 0
      END IF
      IF (debug) write(out_unitp,*) ' option_loc',option_loc

      IF (present(file_name)) THEN
        CALL file_open2(name_file = file_name, iunit=nio)
      ELSE
        nio = out_unitp
      END IF
      IF (debug) write(out_unitp,*) ' nio',nio

!     For ElVibRot (Tnum.op)
      write(nio,'(a)') '-------------------------------------------------'
      write(nio,'(a)') '  Taylor expansion of G at Qact (Taylor_G.keo)'
      write(nio,'(a)') '-------------------------------------------------'
      write(nio,'(" ",a)',advance='no') 'Qact: '
      DO i=1,nb_act
        write(nio,'(" ",f12.6)',advance='no') Qact(i)
      END DO
      write(nio,'(a)',advance='yes')

      SELECT CASE (option_loc)
      CASE (1)
         ! Zero order (cte)
        n0 = 0
        DO i=1,nb_act
        DO j=i,nb_act
          IF (abs(dnGG%d0(i,j)) > epsi_G) n0 = n0+1
        END DO
        END DO
        write(nio,'(a,i0,a)') 'Zero order. ',n0,' terms: '
        write(nio,'(a)') 'i j G(i,j)'
        DO i=1,nb_act
        DO j=i,nb_act
          IF (abs(dnGG%d0(i,j)) > epsi_G) write(nio,*) i,j,dnGG%d0(i,j)
        END DO
        END DO
        !     First order
        n1 = 0
        DO i=1,nb_act
        DO j=i,nb_act
        DO k=1,nb_act
          IF (abs(dnGG%d1(i,j,k)) > epsi_G) n1 = n1+1
        END DO
        END DO
        END DO
        write(nio,'(a,i0,a)') 'First order. ',n1,' terms: '
        write(nio,'(a)') 'i j k d G(i,j)/dQ(k)'
        write(nio,*) n1
        DO k=1,nb_act
        DO i=1,nb_act
        DO j=i,nb_act
          IF (abs(dnGG%d1(i,j,k)) > epsi_G) write(nio,*) i,j,k,dnGG%d1(i,j,k)
        END DO
        END DO
        END DO
        ! second order
        n2 = 0
        DO i=1,nb_act
        DO j=i,nb_act
        DO k=1,nb_act
        DO l=k,nb_act
          IF (abs(dnGG%d2(i,j,k,l)) > epsi_G) n2 = n2+1
        END DO
        END DO
        END DO
        END DO
        write(nio,'(a,i0,a)') 'Second order. ',n2,' terms: '
        write(nio,'(a)') 'i j k l d^2G(i,j)/dQ(k)dQ(l)'
        write(nio,*) n2
        DO k=1,nb_act
        DO l=k,nb_act
        DO i=1,nb_act
        DO j=i,nb_act
          IF (abs(dnGG%d2(i,j,k,l)) > epsi_G)                            &
                              write(nio,*) i,j,k,l,dnGG%d2(i,j,k,l)
        END DO
        END DO
        END DO
        END DO
      CASE DEFAULT
       ! Zero order (cte)
        write(nio,'(a,i0,a)') 'Zero order'
        write(nio,*) dnGG%d0(:,:)
        ! First order
        write(nio,'(a,i0,a)') 'First order'
        DO k=1,nb_act
          write(nio,*) k
          write(nio,*) dnGG%d1(:,:,k)
        END DO
        ! second order
        write(nio,'(a,i0,a)') 'Second order'
        DO k=1,nb_act
        DO l=k,nb_act
          write(nio,*) k,l
          write(nio,*) dnGG%d2(:,:,k,l)
        END DO
        END DO
      END SELECT
      write(nio,'(a)') '-------------------------------------------------'
      write(nio,'(a)') ' END Taylor expansion of G at Qact (Taylor_G.keo)'
      write(nio,'(a)') '-------------------------------------------------'
      flush(nio)
      IF (present(file_name)) close(nio)

!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
      flush(out_unitp)
!     -----------------------------------------------------------------

      end subroutine export_Taylor_dnG

      ! The operartor is d./dQi^i * d./dQj^j Qk^k Ql^l
      SUBROUTINE T_Operator_MCTDH(name_Op,i,j,k,l,n)
      IMPLICIT NONE

      integer :: n
      integer :: i,j,k,l
      character (len=Name_longlen) :: name_Op

      character (len=1) :: sepa
      integer :: iQ,dQi(n),dQj(n),Q(n)
      logical :: First_op,T(n)


      dQi(:) = 0
      dQj(:) = 0
      Q(:) = 0
      T(:) = .FALSE.

      IF (i>0) THEN
        dQi(i) = 1
        T(i)   = .TRUE.
      END IF
      IF (j>0) THEN
        dQj(j) = 1
        T(j)   = .TRUE.
      END IF
      IF (k>0) THEN
        Q(k) = 1
        T(k)   = .TRUE.
      END IF
      IF (l>0) THEN
        Q(l) = Q(l) + 1
        T(l)   = .TRUE.
      END IF

!     write(out_unitp,*) 'dQi',dQi
!     write(out_unitp,*) 'dQj',dQj
!     write(out_unitp,*) ' Q ',Q
!     write(out_unitp,*) 'T: ',T

      name_Op = ' '
      DO iQ=1,n
        IF ( T(iQ) ) THEN

          name_Op = trim(name_Op) // ' |' // int_TO_char(iQ)

          IF ( dQi(iQ) == 1 .AND. dQj(iQ) == 1 .AND. Q(iQ) == 0) THEN
            name_Op = trim(name_Op) // ' dq^2'
            CYCLE
          END IF

          First_Op = .FALSE.
          IF ( dQi(iQ) == 1 ) THEN
            name_Op = trim(name_Op) // ' dq'
            First_op = .TRUE.
          END IF
          IF (First_op) THEN
            sepa = '*'
          ELSE
            sepa = ' '
          END IF
          First_Op = .FALSE.

          IF ( Q(iQ) == 1 ) THEN

            name_Op = trim(name_Op) // sepa // 'NQ' // int_TO_char(iQ)
            First_op = .TRUE.
          ELSE IF ( Q(iQ) == 2 ) THEN
            name_Op = trim(name_Op) // sepa // 'NQ' // int_TO_char(iQ) // '^2'
            First_op = .TRUE.
          END IF
          IF (First_op) THEN
            sepa = '*'
          ELSE
            sepa = ' '
          END IF
          IF ( dQj(iQ) == 1 ) name_Op = trim(name_Op) // sepa // 'dq'
!         write(out_unitp,*) 'iQ, Op',iQ,name_Op
        END IF
      END DO

      end subroutine T_Operator_MCTDH

      FUNCTION get_T_Operator_MCTDH(i,j,k,l,n) RESULT(name_Op)
      IMPLICIT NONE

      integer :: n
      integer :: i,j,k,l

      character (len=:), allocatable           :: name_Op

      character (len=1) :: sepa
      integer :: iQ,dQi(n),dQj(n),Q(n)
      logical :: First_op,T(n)


      dQi(:) = 0
      dQj(:) = 0
      Q(:) = 0
      T(:) = .FALSE.

      IF (i>0) THEN
        dQi(i) = 1
        T(i)   = .TRUE.
      END IF
      IF (j>0) THEN
        dQj(j) = 1
        T(j)   = .TRUE.
      END IF
      IF (k>0) THEN
        Q(k) = 1
        T(k)   = .TRUE.
      END IF
      IF (l>0) THEN
        Q(l) = Q(l) + 1
        T(l)   = .TRUE.
      END IF

!     write(out_unitp,*) 'dQi',dQi
!     write(out_unitp,*) 'dQj',dQj
!     write(out_unitp,*) ' Q ',Q
!     write(out_unitp,*) 'T: ',T

      name_Op = String_TO_String(' ',ltrim=.FALSE.)

      DO iQ=1,n
        IF ( T(iQ) ) THEN

          name_Op = String_TO_String( name_Op // ' |' // int_TO_char(iQ) )

          IF ( dQi(iQ) == 1 .AND. dQj(iQ) == 1 .AND. Q(iQ) == 0) THEN
            name_Op = String_TO_String( name_Op // ' dq^2' )
            CYCLE
          END IF

          First_Op = .FALSE.
          IF ( dQi(iQ) == 1 ) THEN
            name_Op = String_TO_String( name_Op // ' dq' )
            First_op = .TRUE.
          END IF
          IF (First_op) THEN
            sepa = '*'
          ELSE
            sepa = ' '
          END IF
          First_Op = .FALSE.

          IF ( Q(iQ) == 1 ) THEN
            name_Op = String_TO_String( name_Op // sepa // 'NQ' // int_TO_char(iQ) )
            First_op = .TRUE.
          ELSE IF ( Q(iQ) == 2 ) THEN
            name_Op = String_TO_String( name_Op // sepa // 'NQ' // int_TO_char(iQ) // '^2')
            First_op = .TRUE.
          END IF

          IF (First_op) THEN
            sepa = '*'
          ELSE
            sepa = ' '
          END IF

          IF ( dQj(iQ) == 1 ) name_Op = String_TO_String( name_Op // sepa // 'dq')

        END IF
      END DO

      END FUNCTION get_T_Operator_MCTDH

  FUNCTION real_TO_char_MCTDH(r) RESULT(string)
    IMPLICIT NONE

    character (len=:), allocatable           :: string
    real (kind=Rkind), intent(in)            :: r

    integer                                  :: ie

    string = real_TO_char(r)

    ! change the D or E or e character into a "d" character for MCTDH
    ie = scan(string,'DEe')
    IF (ie > 0) string(ie:ie) = 'd'


  END FUNCTION real_TO_char_MCTDH

END MODULE mod_export_KEO

