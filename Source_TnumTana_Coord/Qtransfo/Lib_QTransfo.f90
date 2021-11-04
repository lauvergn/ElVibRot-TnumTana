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
MODULE mod_Lib_QTransfo
      use mod_system
      use mod_dnSVM, only: type_dnvec, type_dns, write_dnsvm,        &
                           sub_dnvec1_prod_dns2_to_dnvec3,           &
                           sub_dns1_prod_dns2_to_dns3,               &
                           sub_dnvec1_plus_dnvec2_to_dnvec3,         &
                           sub_crossproduct_dnvec1_dnvec2_to_dnvec3, &
                           sub_normalize_dnvec, sub_dns_to_dnvec,    &
                           check_alloc_dnvec, dealloc_dnS,           &
                           sub_dnS1_wPLUS_dnS2_TO_dnS3
      IMPLICIT NONE

      PRIVATE
      PUBLIC :: sub3_dnx_AT1, sub3_dnx_AT2_new, sub3_dnx_AT3_new
      PUBLIC :: sub3_dnx_AT4, sub3_dnx_AT4_cart, sub3_dnx_AT4_poly, sub3_dnx_AT4_cart_new
      PUBLIC :: sub3_dnx_TO_dnVec, sub3_dnVec_PLUS_x1TOxf, sub3_dnVec_TOxf
      PUBLIC :: Write_Cart, Write_dnx
      PUBLIC :: calc_vector, calc_vector2, calc_cross_product
      PUBLIC :: calc_angle, calc_angle_d, calc_OutOfPlane
      PUBLIC :: check_Valence, func_ic, func_iat

      CONTAINS

!================================================================
!       atom 1 : d0x(0 , 0, 0)
!================================================================
      SUBROUTINE sub3_dnx_AT1(dnx,icf,nderiv)
      IMPLICIT NONE


        integer :: nb_act,ncart

        integer :: nderiv

        integer :: icf
        TYPE (Type_dnVec) :: dnx


!      -----------------------------------------------------------------
!      logical, parameter :: debug = .TRUE.
       logical, parameter :: debug = .FALSE.
       integer, parameter :: nderiv_debug = 3
      character (len=*), parameter :: name_sub = 'sub3_dnx_AT1'
!      -----------------------------------------------------------------

!      -----------------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'nderiv',nderiv
         write(out_unitp,*) 'icf',icf
       END IF
!      -----------------------------------------------------------------


       IF (nderiv >= 0) dnx%d0(icf:icf+2)       = ZERO
       IF (nderiv >= 1) dnx%d1(icf:icf+2,:)     = ZERO
       IF (nderiv >= 2) dnx%d2(icf:icf+2,:,:)   = ZERO
       IF (nderiv == 3) dnx%d3(icf:icf+2,:,:,:) = ZERO

!      -----------------------------------------------------------------
       IF (debug) THEN
         CALL write_dnx(icf,3,dnx,nderiv_debug)
         write(out_unitp,*) 'END ',name_sub
       END IF
!      -----------------------------------------------------------------

        RETURN
        end subroutine sub3_dnx_AT1
!================================================================
!       atom 2 : d0x(0 , 0, d0d)
!================================================================
      SUBROUTINE sub3_dnx_AT2_new(dnx,icf,ic1,dnd,dnz2,nderiv,check)
      IMPLICIT NONE

      TYPE (Type_dnVec) :: dnx
      integer :: icf,ic1

      TYPE (Type_dnS)   :: dnd
      TYPE (Type_dnVec) :: dnz2

      integer :: nderiv
      logical :: check

!     -----------------------------------------------------------------
      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      integer, parameter :: nderiv_debug = 1
      character (len=*), parameter :: name_sub = 'sub3_dnx_AT2_new'
!     -----------------------------------------------------------------

!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nderiv',nderiv
        write(out_unitp,*) 'icf',icf
        write(out_unitp,*) 'dnd'
        CALL Write_dnSVM(dnd)
      END IF
!     -----------------------------------------------------------------

      IF (dnd%d0 < ONETENTH**10 .AND. check) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'The distance is too small. dnd:'
        CALL Write_dnSVM(dnd)
        STOP
      END IF

      CALL sub_dnVec1_PROD_dnS2_TO_dnVec3(dnz2,dnd,dnz2)

      CALL sub3_dnVec_PLUS_x1TOxf(dnx,icf,ic1,dnz2,nderiv)

!     -----------------------------------------------------------------
      IF (debug) THEN
        CALL write_dnx(icf,3,dnx,nderiv_debug)
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -----------------------------------------------------------------

      end subroutine sub3_dnx_AT2_new
!================================================================
!       atom 3 : with d0d and d0val
!=======================================================================
!
!      ------------
!
!          atf                        ex3
!          /                           ^
!      d  /                            |
!        /                             |
!       /  val                         |
!      at(ic1)---------------at(ic2)   --->ez3
!
!      ex3 and ez3 are the unit vector in the at3=atf frame
!=======================================================================
      SUBROUTINE sub3_dnx_AT3_new(dnx,icf,ic1,check,                    &
                                 dnd,dncval,dnsval,dnz3,dnx3,dnw,nderiv)
      IMPLICIT NONE



        TYPE (Type_dnVec) :: dnx
        integer :: icf,ic1
        logical :: check
        TYPE (Type_dnS) :: dnd,dncval,dnsval
        TYPE (Type_dnVec) :: dnz3,dnx3
        TYPE (Type_dnS) :: dnw
        integer :: nderiv

!     -----------------------------------------------------------------
      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      integer, parameter :: nderiv_debug = 1
      character (len=*), parameter :: name_sub = 'sub3_dnx_AT3_new'
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nderiv',nderiv
        write(out_unitp,*) 'icf ic1',icf,ic1
        write(out_unitp,*) 'dnd'
        CALL Write_dnSVM(dnd)
        write(out_unitp,*) 'dnsin(val)'
        CALL Write_dnSVM(dnsval)
        write(out_unitp,*) 'dncos(val)'
        CALL Write_dnSVM(dncval)
        write(out_unitp,*) 'dnx'
        CALL write_dnx(1,dnx%nb_var_vec,dnx,nderiv_debug)
      END IF
!      -----------------------------------------------------------------


       IF (dnd%d0 < ONETENTH**10  .AND. check) THEN
         write(out_unitp,*) 'ERROR in ',name_sub
         write(out_unitp,*) 'The distance is too small',dnd%d0
         STOP
       END IF


!      -----------------------------------------------------------------
!      for x coordinates
!      d0w = d0d * d0sin * d0x3
       CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnd,dnsval,dnw,nderiv)

       CALL sub_dnVec1_PROD_dnS2_TO_dnVec3(dnx3,dnw,dnx3)

!      -----------------------------------------------------------------


!      -----------------------------------------------------------------
!      for y coordinates
!      icy = icf+1
!      d0x(icf+1)=0   has been done in the initialization
!      -----------------------------------------------------------------

!      -----------------------------------------------------------------
!      for z coordinates
!      d0w = d0d * d0cos * d0z3
       CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnd,dncval,dnw,nderiv)

       CALL sub_dnVec1_PROD_dnS2_TO_dnVec3(dnz3,dnw,dnz3)

!      -----------------------------------------------------------------
!      add the x and z contribution to d0x
       CALL sub_dnVec1_PLUS_dnVec2_TO_dnVec3(dnx3,dnz3,dnz3,nderiv)

       CALL sub3_dnVec_PLUS_x1TOxf(dnx,icf,ic1,dnz3,nderiv)
!      -----------------------------------------------------------------

!      -----------------------------------------------------------------
       IF (debug) THEN
         CALL write_dnx(icf,3,dnx,nderiv_debug)
         write(out_unitp,*) 'END ',name_sub
       END IF
!      -----------------------------------------------------------------

        end subroutine sub3_dnx_AT3_new
!================================================================
!       atom 4 : with d0d, (d0cval,d0sval) and (d0cdih,d0sdih)
!================================================================
!       genere un atome fictif suivant une distance, un angle de valence
!          et un angle dihedre
!       nx n1 d n2 val n3 dih
!       avec   d   = distance(nx,n1)
!              val = angle(nx,n1,n2)
!              dih = angle_d(nx,n1,n2,n3)
!================================================================
!
      SUBROUTINE sub3_dnx_AT4(dnx,icf,ic1,ic2,ic3,check,                &
                               dnd,dncval,dnsval,dncdih,dnsdih,         &
                               dnv1,dnv2,dnv3,dnf1,dnf2,dnf3,           &
                               nderiv)
      IMPLICIT NONE

      TYPE (Type_dnVec) :: dnx
      integer :: icf,ic1,ic2,ic3
      logical :: check

      TYPE (Type_dnS) :: dnd,dncval,dnsval,dncdih,dnsdih

        TYPE (Type_dnVec) :: dnv1,dnv2,dnv3
        TYPE (Type_dnS) :: dnf1,dnf2,dnf3
        integer :: nderiv

!     -----------------------------------------------------------------
!     logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      integer, parameter :: nderiv_debug = 3
      character (len=*), parameter :: name_sub = 'sub3_dnx_AT4'
!     -----------------------------------------------------------------

!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nderiv',nderiv
        write(out_unitp,*) 'icf,ic1,ic2,ic3',icf,ic1,ic2,ic3

        write(out_unitp,*) 'dnd'
        CALL Write_dnSVM(dnd)
        write(out_unitp,*) 'dnsin(val)'
        CALL Write_dnSVM(dnsval)
        write(out_unitp,*) 'dncos(val)'
        CALL Write_dnSVM(dncval)
        write(out_unitp,*) 'dnsin(dih)'
        CALL Write_dnSVM(dnsdih)
        write(out_unitp,*) 'dncos(dih)'
        CALL Write_dnSVM(dncdih)

        write(out_unitp,*) 'dnx'
        CALL write_dnx(1,dnx%nb_var_vec,dnx,nderiv_debug)
      END IF
!     -----------------------------------------------------------------

!================================================

      IF (dnd%d0 < ONETENTH**10 .AND. check) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'The distance is too small',dnd%d0
        STOP
      END IF


      CALL sub3_dnx_TO_dnVec(dnx,ic2,ic3,dnv1,nderiv)
      CALL sub3_dnx_TO_dnVec(dnx,ic2,ic1,dnv2,nderiv)

      CALL Sub_crossproduct_dnVec1_dnVec2_TO_dnVec3(dnv2,dnv1,dnv3,nderiv)

!     d0v1 is reused instead of d0v4
      CALL Sub_crossproduct_dnVec1_dnVec2_TO_dnVec3(dnv3,dnv2,dnv1,nderiv)

!     ----------------------------------------------------------
!     d0v2 d0v3 d0v1 normalization
      CALL sub_Normalize_dnVec(dnv2)
      CALL sub_Normalize_dnVec(dnv3)


!     d0v1 (d0v4) norm is only the product of d0norm2 and d0norm3
!     since d0v1 = d0v2*d0v3 and d0v2 is perpendicular to d0v3
!     simplification not used yet !!!!
      CALL sub_Normalize_dnVec(dnv1)

!     -----------------------------------------------------------------
!     d0f2 = d0d * d0sval (tempory)
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnd,dnsval,dnf2,nderiv)

!     d0f1 = d0f2 * d0cdih  = (d0d * d0sval) * d0cdih
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnf2,dncdih,dnf1,nderiv)

!     d0f3 = d0f2 * d0sdih  = (d0d * d0sval) * d0sdih
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnf2,dnsdih,dnf3,nderiv)
!      -----------------------------------------------------------------


!     -----------------------------------------------------------------
!     d0f2 = d0d * d0cval
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnd,dncval,dnf2,nderiv)
!     -----------------------------------------------------------------


!     -----------------------------------------------------------------
!     -----------------------------------------------------------------
!     NOW :
!     d0f1 = (d0d * d0sval) * d0cdih    : associated with d0v1
!     d0f2 =  d0d * d0cval              : associated with d0v2
!     d0f3 = (d0d * d0sval) * d0sdih    : associated with d0v3
!     -----------------------------------------------------------------
!     -----------------------------------------------------------------

!     ---------------------------------------------------------
!     d0v1 = d0v1 * d0f1
!     d0v2 = d0v2 * d0f2
!     d0v3 = d0v3 * d0f3

      CALL sub_dnVec1_PROD_dnS2_TO_dnVec3(dnv1,dnf1,dnv1)
      CALL sub_dnVec1_PROD_dnS2_TO_dnVec3(dnv2,dnf2,dnv2)
      CALL sub_dnVec1_PROD_dnS2_TO_dnVec3(dnv3,dnf3,dnv3)

!     ----------------------------------------------------------

!     ---------------------------------------------------------
!     d0x(icf) = d0x(ic1) - d0v2 + d0v1 + d0v3
      CALL sub3_dnVec123_TOx(dnx,icf,ic1,dnv2,dnv1,dnv3,nderiv)
!      ---------------------------------------------------------


!     -----------------------------------------------------------------
      IF (debug) THEN
        CALL write_dnx(1,dnx%nb_var_vec,dnx,nderiv_debug)
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -----------------------------------------------------------------

      end subroutine sub3_dnx_AT4
!================================================================
!       atom 4 :  en Cartesienne
!================================================================
!
      SUBROUTINE sub3_dnx_AT4_cart(dnx,icf,dna,dnb,dnc,nderiv)
      IMPLICIT NONE

      TYPE (Type_dnVec) :: dnx
      integer :: icf

      TYPE (Type_dnS) :: dna,dnb,dnc
      integer :: nderiv


!     -----------------------------------------------------------------
!      logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      integer, parameter :: nderiv_debug = 1
      character (len=*), parameter :: name_sub = 'sub3_dnx_AT4_cart'
!     -----------------------------------------------------------------

!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nderiv',nderiv
        write(out_unitp,*) 'icf',icf
        write(out_unitp,*) 'dna'
        CALL Write_dnSVM(dna,nderiv=nderiv_debug)
        write(out_unitp,*) 'dnb'
        CALL Write_dnSVM(dnb,nderiv=nderiv_debug)
        write(out_unitp,*) 'dnb'
        CALL Write_dnSVM(dnb,nderiv=nderiv_debug)
      END IF
!     -----------------------------------------------------------------


      CALL sub_dnS_TO_dnVec(dna,dnx,icf,nderiv)
      icf=icf+1
      CALL sub_dnS_TO_dnVec(dnb,dnx,icf,nderiv)
      icf=icf+1
      CALL sub_dnS_TO_dnVec(dnc,dnx,icf,nderiv)

!     -----------------------------------------------------------------
      IF (debug) THEN
        icf=icf-2
        write(out_unitp,*) 'icf',icf
        CALL write_dnx(1,dnx%nb_var_vec,dnx,nderiv_debug)
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -----------------------------------------------------------------

      end subroutine sub3_dnx_AT4_cart
!================================================================
!       atom 4 :  en Cartesienne
!================================================================
!
      SUBROUTINE sub3_dnx_AT4_cart_new(dnx,icf,dna,dnb,dnc,dnx3,dny3,dnz3,nderiv)
      IMPLICIT NONE

      TYPE (Type_dnVec) :: dnx
      integer           :: icf

      TYPE (Type_dnS)   :: dna,dnb,dnc
      TYPE (Type_dnVec) :: dnx3,dny3,dnz3
      integer           :: nderiv


      TYPE (Type_dnS)   :: dnS1,dnS2
!     -----------------------------------------------------------------
      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      integer, parameter :: nderiv_debug = 1
      character (len=*), parameter :: name_sub = 'sub3_dnx_AT4_cart_new'
!     -----------------------------------------------------------------

!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nderiv',nderiv
        write(out_unitp,*) 'icf',icf
        write(out_unitp,*) 'dna'
        CALL Write_dnSVM(dna,nderiv=nderiv_debug)
        write(out_unitp,*) 'dnb'
        CALL Write_dnSVM(dnb,nderiv=nderiv_debug)
        write(out_unitp,*) 'dnc'
        CALL Write_dnSVM(dnc,nderiv=nderiv_debug)

        write(out_unitp,*) 'dnx3'
        CALL Write_dnSVM(dnx3,nderiv=nderiv_debug)
        write(out_unitp,*) 'dny3'
        CALL Write_dnSVM(dny3,nderiv=nderiv_debug)
        write(out_unitp,*) 'dnz3'
        CALL Write_dnSVM(dnz3,nderiv=nderiv_debug)
      END IF
!     -----------------------------------------------------------------

      ! x component
      !dnS = dna * dnx3%d0(1) + dnb * dny3%d0(1) + dnc * dnz3%d0(1)
      CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dna,dnx3%d0(1),dnb,dny3%d0(1),dnS1,nderiv)
      CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnS1,ONE,dnc,dnz3%d0(1),dnS2,nderiv)

      CALL sub_dnS_TO_dnVec(dnS2,dnx,icf,nderiv)

      ! y component
      !dnS = dna * dnx3%d0(2) + dnb * dny3%d0(2) + dnc * dnz3%d0(2)
      CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dna,dnx3%d0(2),dnb,dny3%d0(2),dnS1,nderiv)
      CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnS1,ONE,dnc,dnz3%d0(2),dnS2,nderiv)

      icf=icf+1
      CALL sub_dnS_TO_dnVec(dnS2,dnx,icf,nderiv)

      ! z component
      !dnS = dna * dnx3%d0(3) + dnb * dny3%d0(3) + dnc * dnz3%d0(3)
      CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dna,dnx3%d0(3),dnb,dny3%d0(3),dnS1,nderiv)
      CALL sub_dnS1_wPLUS_dnS2_TO_dnS3(dnS1,ONE,dnc,dnz3%d0(3),dnS2,nderiv)

      icf=icf+1
      CALL sub_dnS_TO_dnVec(dnS2,dnx,icf,nderiv)


      ! deallocate
      CALL dealloc_dnS(dnS1)
      CALL dealloc_dnS(dnS2)


!     -----------------------------------------------------------------
      IF (debug) THEN
        icf=icf-2
        write(out_unitp,*) 'icf',icf
        CALL write_dnx(1,dnx%nb_var_vec,dnx,nderiv_debug)
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -----------------------------------------------------------------

      end subroutine sub3_dnx_AT4_cart_new

!================================================================
!       atom 4 : with d0d, (d0cval,d0sval) and (d0cdih,d0sdih)
!       with polyspherical coordinates
!
!
!          atf
!          /                      x
!      d  /                       |
!        /                        |
!       /  val                    |
!      at1---------------at2      --->z
!
!
!      at(icf) : at(ic1) + ( d*sin(val)*cos(dih) ; d*sin(val)*sin(dih) ; d*cos(val) )
!
!================================================================
!
      SUBROUTINE sub3_dnx_AT4_poly(dnx,icf,ic1,ic2,ic3,check,           &
                                   dnd,dncval,dnsval,dncdih,dnsdih,     &
                                   dnv1,dnf1,dnf2,dnf3,                 &
                                   nderiv)
      IMPLICIT NONE

      TYPE (Type_dnVec) :: dnx
      integer :: icf,ic1,ic2,ic3
      logical :: check

      TYPE (Type_dnS) :: dnd,dncval,dnsval,dncdih,dnsdih

      TYPE (Type_dnS) :: dnf1,dnf2,dnf3
      TYPE (Type_dnVec) :: dnv1
      integer :: nderiv

!     -----------------------------------------------------------------
      !logical, parameter :: debug = .TRUE.
      logical, parameter :: debug = .FALSE.
      integer, parameter :: nderiv_debug = 0
      character (len=*), parameter :: name_sub = 'sub3_dnx_AT4_poly'
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nderiv',nderiv
        write(out_unitp,*) 'icf,ic1,ic2,ic3',icf,ic1,ic2,ic3

        write(out_unitp,*) 'dnd'
        CALL Write_dnSVM(dnd)
        write(out_unitp,*) 'dnsin(val)'
        CALL Write_dnSVM(dnsval)
        write(out_unitp,*) 'dncos(val)'
        CALL Write_dnSVM(dncval)
        write(out_unitp,*) 'dnsin(dih)'
        CALL Write_dnSVM(dnsdih)
        write(out_unitp,*) 'dncos(dih)'
        CALL Write_dnSVM(dncdih)

        write(out_unitp,*) 'dnx'
        CALL write_dnx(1,dnx%nb_var_vec,dnx,nderiv_debug)
      END IF
!     -----------------------------------------------------------------


!================================================

      IF (dnd%d0 < ONETENTH**10 .AND. check) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) 'The distance is too small',dnd%d0
        STOP
      END IF


!     -----------------------------------------------------------------
!     d0f3 = d0d * d0sval (tempory)
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnd,dnsval,dnf3,nderiv)

!     d0f1 = d0f3 * d0cdih  = (d0d * d0sval) * d0cdih
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnf3,dncdih,dnf1,nderiv)

!     d0f2 = d0f3 * d0sdih  = (d0d * d0sval) * d0sdih
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnf3,dnsdih,dnf2,nderiv)
!     -----------------------------------------------------------------


!     -----------------------------------------------------------------
!     d0f3 = d0d * d0cval
   CALL sub_dnS1_PROD_dnS2_TO_dnS3(dnd,dncval,dnf3,nderiv)
!     -----------------------------------------------------------------


!      ----------------------------------------------------------------
!      ----------------------------------------------------------------
!      NOW :
!      d0f1 = (d0d * d0sval) * d0cdih
!      d0f2 = (d0d * d0sval) * d0sdih
!      d0f3 =  d0d * d0cval
!      -----------------------------------------------------------------
      CALL sub_dnS_TO_dnVec(dnf1,dnv1,1,nderiv)
      CALL sub_dnS_TO_dnVec(dnf2,dnv1,2,nderiv)
      CALL sub_dnS_TO_dnVec(dnf3,dnv1,3,nderiv)
!     -----------------------------------------------------------------
      IF (debug) write(out_unitp,*) 'ic1,dnx',ic1,dnx%d0(ic1:ic1+2)
      IF (debug) write(out_unitp,*) 'dnv1'
      IF (debug) CALL Write_dnSVM(dnv1)

!     ---------------------------------------------------------
!     d0x(icf) = d0x(ic1) + (d0f1,d0f2,d0f3)
!     d0x(icf) = d0x(ic1) + d0v1

      CALL sub3_dnVec_PLUS_x1TOxf(dnx,icf,ic1,dnv1,nderiv)

!     ---------------------------------------------------------


!     -----------------------------------------------------------------
      IF (debug) THEN
        CALL write_dnx(1,dnx%nb_var_vec,dnx,nderiv_debug)
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -----------------------------------------------------------------

      end subroutine sub3_dnx_AT4_poly


!================================================================
!       vector d0v = d0x2 - d0x1
!================================================================
      SUBROUTINE sub3_dnx_TO_dnVec(dnx,ic1,ic2,dnVec,nderiv)
      IMPLICIT NONE


      TYPE (Type_dnVec) :: dnx
      integer           :: ic1,ic2
      TYPE (Type_dnVec) :: dnVec
      integer           :: nderiv

       integer          :: icx1,icx2
       integer          :: icz1,icz2

!      -----------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='sub3_dnx_TO_dnVec'
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nderiv',nderiv
        write(out_unitp,*) 'ic1 ic2',ic1,ic2
        write(out_unitp,*) 'dnx'
        CALL write_dnx(1,dnx%nb_var_vec,dnx,nderiv)
      END IF
!     -----------------------------------------------------------------

      icx1 = ic1+0
      icz1 = ic1+2

      icx2 = ic2+0
      icz2 = ic2+2
!     -----------------------------------------------------------------
!     vector
!     -----------------------------------------------------------------
      IF (nderiv >= 0) THEN
        dnVec%d0(:) = dnx%d0(icx2:icz2) - dnx%d0(icx1:icz1)
      END IF
!     -----------------------------------------------------------------
      IF (nderiv >= 1) THEN
        dnVec%d1(:,:) = dnx%d1(icx2:icz2,:) - dnx%d1(icx1:icz1,:)
      END IF
!     -----------------------------------------------------------------
      IF (nderiv >= 2) THEN
        dnVec%d2(:,:,:) = dnx%d2(icx2:icz2,:,:) - dnx%d2(icx1:icz1,:,:)
      END IF
!     -----------------------------------------------------------------
      IF (nderiv >= 3) THEN
         dnVec%d3(:,:,:,:) =                                            &
                    dnx%d3(icx2:icz2,:,:,:) - dnx%d3(icx1:icz1,:,:,:)
      END IF
!     -----------------------------------------------------------------

!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'dnVec'
        CALL Write_dnSVM(dnVec)
        write(out_unitp,*) 'END ',name_sub
      END IF
!      -----------------------------------------------------------------

       end subroutine sub3_dnx_TO_dnVec
!================================================================
!       dnx%d0(0   or  1  or 2 ) = d0f
!    ic =   icx or icy or icz
!
!     SUBROUTINE d0d1d2d3fTOx(nb_act,ncart,
! remplacer par sub_dnS_TO_dnVec(dnR,dnVec,iVec)
!================================================================

!================================================================
!      dnx%d0(icf) = dnx%d0(ic1) - d0v1 + d0v2 + d0v3
!================================================================
      SUBROUTINE sub3_dnVec123_TOx(dnx,icf,ic1,                         &
                                   dnVec1,dnVec2,dnVec3,                &
                                   nderiv)
      IMPLICIT NONE


      TYPE (Type_dnVec) :: dnx
      integer           :: icf,ic1
      TYPE (Type_dnVec) :: dnVec1,dnVec2,dnVec3
      integer           :: nderiv

       integer i,j,k
       integer ic1x,ic1y,ic1z
       integer icfx,icfy,icfz
!      -----------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='sub3_dnVec123_TOx'
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nderiv',nderiv
        write(out_unitp,*) 'icf,ic1 ',icf,ic1
        write(out_unitp,*) 'dnx'
        CALL write_dnx(1,dnx%nb_var_vec,dnx,nderiv)
      END IF
!     -----------------------------------------------------------------

!      -----------------------------------------------------------------

       icfx = icf + 0
       icfz = icf + 2

       ic1x = ic1 + 0
       ic1z = ic1 + 2

!     -----------------------------------------------------------------
      dnx%d0(icf:icfz) = dnx%d0(ic1:ic1z) -                             &
                 dnVec1%d0(:) + dnVec2%d0(:) + dnVec3%d0(:)
!     -----------------------------------------------------------------
      IF (nderiv .GE. 1) THEN
        dnx%d1(icf:icfz,:) = dnx%d1(ic1:ic1z,:) -                       &
                   dnVec1%d1(:,:) + dnVec2%d1(:,:) + dnVec3%d1(:,:)
      END IF
!     -----------------------------------------------------------------
      IF (nderiv .GE. 2) THEN
        dnx%d2(icf:icfz,:,:) = dnx%d2(ic1:ic1z,:,:) -                   &
                 dnVec1%d2(:,:,:) + dnVec2%d2(:,:,:) + dnVec3%d2(:,:,:)
      END IF
!     -----------------------------------------------------------------
      IF (nderiv .GE. 3) THEN
        dnx%d3(icf:icfz,:,:,:) = dnx%d3(ic1:ic1z,:,:,:) -               &
            dnVec1%d3(:,:,:,:) + dnVec2%d3(:,:,:,:) + dnVec3%d3(:,:,:,:)
      END IF

!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'dnx'
        CALL write_dnx(icf,3,dnx,nderiv)
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -----------------------------------------------------------------

      end subroutine sub3_dnVec123_TOx
!================================================================
!      dnx%d0(icf) = dnx%d0(ic1) + d0v1
!================================================================
      SUBROUTINE sub3_dnVec_PLUS_x1TOxf(dnx,icf,ic1,dnVec,nderiv)
      IMPLICIT NONE


      TYPE (Type_dnVec) :: dnx
      integer           :: ic1,icf
      TYPE (Type_dnVec) :: dnVec
      integer           :: nderiv

!      -----------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='sub3_dnVec_PLUS_x1TOxf'
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nderiv',nderiv
        write(out_unitp,*) 'ic1 icf',ic1,icf
        write(out_unitp,*) 'dnx,ic1'
        CALL write_dnx(ic1,3,dnx,nderiv)
      END IF
!     -----------------------------------------------------------------

!      -----------------------------------------------------------------
       dnx%d0(icf:icf+2) = dnx%d0(ic1:ic1+2) + dnVec%d0(:)
!      -----------------------------------------------------------------
       IF (nderiv .GE. 1) THEN
         dnx%d1(icf:icf+2,:) = dnx%d1(ic1:ic1+2,:) + dnVec%d1(:,:)
       END IF
!      -----------------------------------------------------------------
       IF (nderiv .GE. 2) THEN
         dnx%d2(icf:icf+2,:,:) = dnx%d2(ic1:ic1+2,:,:) + dnVec%d2(:,:,:)
       END IF
!      -----------------------------------------------------------------
       IF (nderiv .GE. 3) THEN
         dnx%d3(icf:icf+2,:,:,:) = dnx%d3(ic1:ic1+2,:,:,:) +            &
                                                       dnVec%d3(:,:,:,:)
       END IF

!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'dnx,icf'
        CALL write_dnx(icf,3,dnx,nderiv)
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -----------------------------------------------------------------

       end subroutine sub3_dnVec_PLUS_x1TOxf
!================================================================
!      dnx%d0(icf:icf+2) = d0Vect(1:3)
!================================================================
      SUBROUTINE sub3_dnVec_TOxf(dnx,icf,dnVec,nderiv)
      IMPLICIT NONE


      TYPE (Type_dnVec) :: dnx
      integer           :: icf
      TYPE (Type_dnVec) :: dnVec
      integer           :: nderiv

!      -----------------------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub='sub3_dnVec_TOxf'
!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nderiv',nderiv
        write(out_unitp,*) 'icf ',icf
        write(out_unitp,*) 'dnx'
        CALL write_dnx(1,dnx%nb_var_vec,dnx,nderiv)
      END IF
!     -----------------------------------------------------------------

!      -----------------------------------------------------------------
       dnx%d0(icf:icf+2) = dnVec%d0(:)
!      -----------------------------------------------------------------
       IF (nderiv .GE. 1) THEN
         dnx%d1(icf:icf+2,:) = dnVec%d1(:,:)
       END IF
!      -----------------------------------------------------------------
       IF (nderiv .GE. 2) THEN
         dnx%d2(icf:icf+2,:,:) = dnVec%d2(:,:,:)
       END IF
!      -----------------------------------------------------------------
       IF (nderiv .GE. 3) THEN
         dnx%d3(icf:icf+2,:,:,:) = dnVec%d3(:,:,:,:)
       END IF

!     -----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'dnx'
        CALL write_dnx(icf,3,dnx,nderiv)
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -----------------------------------------------------------------

       end subroutine sub3_dnVec_TOxf

!================================================================
!       Write Cartesian coordinates
!================================================================
      SUBROUTINE Write_Cart(ncart,d0x)
      IMPLICIT NONE

      integer           :: ncart
      real (kind=Rkind) :: d0x(ncart)

      integer       :: i


      IF (sqrt(dot_product(d0x,d0x)) < ONETENTH**6) THEN

        write(out_unitp,*) 'coordinates : 0'

      ELSE

        DO i=1,ncart,3
          write(out_unitp,111) d0x(i+0),d0x(i+1),d0x(i+2)
 111      format(1x,3(2x,f20.9))
        END DO

      END IF

      END SUBROUTINE Write_Cart
      SUBROUTINE Write_dnx(ic,ncart_e,dnx,nderiv)
      IMPLICIT NONE

        integer :: ic,ncart_e,ncart,nb_act,nderiv,nderiv_loc
        TYPE (Type_dnVec) :: dnx

        integer :: i,j,k
        character (len=*), parameter :: name_sub = 'write_dnx'

        CALL check_alloc_dnVec(dnx,'dnx',name_sub)

        nderiv_loc = min(nderiv,dnx%nderiv)
        nb_act = dnx%nb_var_deriv

        IF (ncart_e > dnx%nb_var_vec) THEN
          write(out_unitp,*) ' ERROR in write_dnx'
          write(out_unitp,*) ' ncart_e > nb_var_vec',ncart_e,dnx%nb_var_vec
          write(out_unitp,*) ' Check the fortran !!!!'
          STOP
        END IF

         write(out_unitp,*) 'd0x ='
         CALL Write_Cart(ncart_e,dnx%d0(ic:ic+ncart_e-1))

         IF (nderiv_loc >= 1) THEN
           DO i=1,nb_act
             write(out_unitp,*) 'd1x =',i
             CALL Write_Cart(ncart_e,dnx%d1(ic:ic+ncart_e-1,i))
           END DO
         END IF

         IF (nderiv_loc >= 2 ) THEN
           DO i=1,nb_act
           DO j=1,nb_act
             write(out_unitp,*) 'd2x =',i,j
             CALL Write_Cart(ncart_e,dnx%d2(ic:ic+ncart_e-1,i,j))
           END DO
           END DO
         END IF

         IF (nderiv_loc >= 3) THEN
           DO i=1,nb_act
           DO j=1,nb_act
           DO k=1,nb_act
             write(out_unitp,*) 'd3x =',i,j,k
             CALL Write_Cart(ncart_e,dnx%d3(ic:ic+ncart_e-1,i,j,k))
           END DO
           END DO
           END DO
         END IF

      END SUBROUTINE write_dnx

!================================================================
!       vector n1-n2
!================================================================
      SUBROUTINE calc_vector(v,norm,n1,n2,x,ndim)
      IMPLICIT NONE

      integer           :: n1,n2
      integer           :: ndim
      real (kind=Rkind) :: x(ndim),v(3),norm

      integer           :: nc1,nc2

      nc1 = 3*n1-2
      nc2 = 3*n2-2

      v(:) = x(nc2+0:nc2+2) - x(nc1+0:nc1+2)

      norm = sqrt( dot_product(v,v))

!     write(out_unitp,*) 'v,n1,n2,norm',v,n1,n2,norm
!      write(out_unitp,*) 'v,nc1,nc2,norm',v,nc1,nc2,norm


      end subroutine calc_vector
      SUBROUTINE calc_vector2(v,norm,nc1,nc2,x,ndim)
      IMPLICIT NONE

      integer           :: nc1,nc2
      integer           :: ndim
      real (kind=Rkind) :: x(ndim),v(3),norm

      !write(out_unitp,*) 'X1',nc1,x(nc1+0:nc1+2)
      !write(out_unitp,*) 'X2',nc2,x(nc2+0:nc2+2)

      v(:) = x(nc2+0:nc2+2) - x(nc1+0:nc1+2)

      norm = sqrt( dot_product(v,v))

      !write(out_unitp,*) 'v,nc1,nc2,norm',v,nc1,nc2,norm


      end subroutine calc_vector2
!================================================================
!       angle
!================================================================
      SUBROUTINE calc_angle(angle,v1,norm1,v2,norm2)
      IMPLICIT NONE

      real (kind=Rkind) :: v1(3),norm1
      real (kind=Rkind) :: v2(3),norm2
      real (kind=Rkind) :: angle

      real (kind=Rkind) :: c

!     write(out_unitp,*) 'v1,norm1',v1,norm1
!     write(out_unitp,*) 'v2,norm2',v2,norm2

      IF (abs(norm1) < ONETENTH**10 .OR. abs(norm2) < ONETENTH**10) THEN
        write(out_unitp,*) 'ERROR in angle'
        write(out_unitp,*)  ' norm1 = 0 or norm2 = 0'
        write(out_unitp,*) 'norm1,norm2',norm1,norm2
        STOP
      END IF

      c = dot_product(v1,v2)/(norm1*norm2)

      IF (c < -ONE) THEN
        angle = pi
      ELSE IF (c > ONE) THEN
        angle = ZERO
      ELSE
        angle = acos(c)
      END IF
      !write(out_unitp,*) 'angle',angle

      end subroutine calc_angle
!================================================================
!       produit vectoriel
!================================================================
      SUBROUTINE calc_cross_product(v1,norm1,v2,norm2,v3,norm3)
      IMPLICIT NONE

      real (kind=Rkind) :: v1(3),norm1
      real (kind=Rkind) :: v2(3),norm2
      real (kind=Rkind) :: v3(3),norm3

      v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
      v3(2) =-v1(1)*v2(3) + v1(3)*v2(1)
      v3(3) = v1(1)*v2(2) - v1(2)*v2(1)

      norm3 = sqrt(dot_product(v3,v3))


      end subroutine calc_cross_product
!================================================================
!       angle oriente (dihedre)
!================================================================
      SUBROUTINE calc_angle_d(angle_d,v1,norm1,v2,norm2,v3,norm3)
      IMPLICIT NONE

      real (kind=Rkind) :: v1(3),norm1
      real (kind=Rkind) :: v2(3),norm2
      real (kind=Rkind) :: v3(3),norm3
      real (kind=Rkind) :: angle_d
      real (kind=Rkind) :: ca,sa

!     write(out_unitp,*) 'v1,norm1',v1,norm1
!     write(out_unitp,*) 'v2,norm2',v2,norm2
!     write(out_unitp,*) 'v3,norm3',v3,norm3

      IF (abs(norm1) < ONETENTH**10 .OR.                                &
          abs(norm2) < ONETENTH**10 .OR.                                &
          abs(norm3) < ONETENTH**10) THEN
        write(out_unitp,*) 'ERROR in angle'
        write(out_unitp,*)  ' norm1 = 0 or norm2 = 0 or norm3 = 0'
        write(out_unitp,*) 'norm1,norm2,norm3',norm1,norm2,norm3
        STOP
      END IF

      ca = dot_product(v1,v2)/(norm1*norm2)

      sa = (v1(1)*(v2(2)*v3(3)-v2(3)*v3(2))                             &
           -v1(2)*(v2(1)*v3(3)-v2(3)*v3(1))                             &
           +v1(3)*(v2(1)*v3(2)-v2(2)*v3(1)))                            &
           /(norm1*norm2*norm3)

      angle_d = atan2(sa,ca)
!     write(out_unitp,*) 'angle_d',angle_d

      end subroutine calc_angle_d
!================================================================
!       out-of-plane angle : theta
!  See Wilson, Decius and Cross, Molecular Vibrations pp58-59
!
!                   3
!                  /
!                 /
!    1-----------4  ) phi1
!                 \
!                  \
!                   \
!                    2
!
!    sin(theta) = (e42 x e43).e41 / sin(phi1)
!               = (v42 x v43).v41 / sin(phi1) / (norm41*norm42*norm43)
!       OR
!    sin(theta) = (v42 x v43).v41 / ( norm(v42 x v43)*norm41 )
!
!    v41 = v1 ; v42 = v2 ; v43 = v3
!================================================================
      SUBROUTINE calc_OutOfPlane(theta,v1,norm1,v2,norm2,v3,norm3)
      IMPLICIT NONE

      real (kind=Rkind) :: v1(3),norm1
      real (kind=Rkind) :: v2(3),norm2
      real (kind=Rkind) :: v3(3),norm3
      real (kind=Rkind) :: vcp(3),normcp
      real (kind=Rkind) :: theta,stheta

!     write(out_unitp,*) 'v1,norm1',v1,norm1
!     write(out_unitp,*) 'v2,norm2',v2,norm2
!     write(out_unitp,*) 'v3,norm3',v3,norm3

      IF (abs(norm1) < ONETENTH**10 .OR.                                &
          abs(norm2) < ONETENTH**10 .OR.                                &
          abs(norm3) < ONETENTH**10) THEN
        write(out_unitp,*) 'ERROR in angle'
        write(out_unitp,*)  ' norm1 = 0 or norm2 = 0 or norm3 = 0'
        write(out_unitp,*) 'norm1,norm2,norm3',norm1,norm2,norm3
        STOP
      END IF

      CALL calc_cross_product(v2,norm2,v3,norm3,vcp,normcp)

      stheta = dot_product(vcp,v1)/(norm1*normcp)

      theta = asin(stheta)


!     write(out_unitp,*) 'OutOfPlane',theta

      end subroutine calc_OutOfPlane


!================================================================
!     Check the range of the valence angle
!================================================================
      SUBROUTINE check_Valence(iz,Q,type_Q)
      IMPLICIT NONE

      integer :: iz,type_Q
      real (kind=Rkind) :: Q

      integer, parameter :: max_cc = 500
      logical, save :: liste_pb(max_cc) = .FALSE.

      IF (print_level < 0) RETURN
!     write(out_unitp,*) 'iz,Q,type_Q',iz,Q,type_Q
      IF (iz > max_cc) RETURN

!$OMP CRITICAL (check_Valence_CRIT)

      IF ( type_Q == 3 .AND. (Q <= ZERO .OR. Q >= Pi)                   &
                                         .AND. .NOT. liste_pb(iz) ) THEN
        write(out_unitp,*) ' WARNNING in check_Valence'
        write(out_unitp,*) ' The value of the valence angle (',iz,') : ',Q
        write(out_unitp,*) ' is out of range ]0,Pi(',pi,')['
        write(out_unitp,*) ' type_Q',type_Q
        write(out_unitp,*) ' ONLY one warning for this coordinate'
        write(out_unitp,*) ' END WARNNING in check_Valence'

        liste_pb(iz) = .TRUE.
      ELSE IF ( type_Q == -3 .AND. abs(Q) >= ONE                        &
                               .AND. .NOT. liste_pb(iz) ) THEN
        write(out_unitp,*) ' WARNNING in check_Valence'
        write(out_unitp,*) ' The value of the cos(valence angle) (',iz,') :',Q
        write(out_unitp,*) ' is out of range ]-1,1['
        write(out_unitp,*) ' type_Q',type_Q
        write(out_unitp,*) ' ONLY one warning for this coordinate'
        write(out_unitp,*) ' END WARNNING in check_Valence'
        liste_pb(iz) = .TRUE.
      ELSE IF ( type_Q == 2 .AND. Q <= ZERO                             &
                               .AND. .NOT. liste_pb(iz) ) THEN
        write(out_unitp,*) ' WARNNING in check_Valence'
        write(out_unitp,*) ' The value of the distance (',iz,') :',Q
        write(out_unitp,*) ' is out of range ]0,+inf['
        write(out_unitp,*) ' type_Q',type_Q
        write(out_unitp,*) ' ONLY one warning for this coordinate'
        write(out_unitp,*) ' END WARNNING in check_Valence'
        liste_pb(iz) = .TRUE.
      END IF
!$OMP END CRITICAL (check_Valence_CRIT)

      end subroutine check_Valence
!================================================================
!     Calculation of the index ic in d0x(..) with the atomic number i
!================================================================
      FUNCTION func_ic(i)
      IMPLICIT NONE
      integer :: func_ic

      integer :: i

!     write(out_unitp,*) 'func_ic',i,3*i-2
      IF (i > 0) THEN
        func_ic = 3*i-2
      ELSE
        func_ic = 3*i+2
      END IF

      end function func_ic
      FUNCTION func_iat(if)
      IMPLICIT NONE
      integer :: func_iat

      integer :: if

!     write(out_unitp,*) 'func_iat',i,3*i-2
      IF (if > 0) THEN
        func_iat = (if+2)/3
      ELSE
        STOP 'ERROR in func_iat, if<1'
      END IF

      end function func_iat

END MODULE mod_Lib_QTransfo
