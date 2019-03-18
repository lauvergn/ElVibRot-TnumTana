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
MODULE mod_basis
      USE mod_nDindex
      USE mod_Coord_KEO

      USE mod_RotBasis_Param
      USE mod_Basis_Grid_Param
      USE mod_Basis_L_TO_n
      USE mod_SymAbelian
      USE mod_param_SGType2

      USE mod_basis_set_alloc
      USE mod_basis_BtoG_GtoB
      IMPLICIT NONE

      CONTAINS

      RECURSIVE SUBROUTINE Set_basis_para_FOR_optimization(basis_set,Set_Val)
      USE mod_system
      IMPLICIT NONE


!----- for the zmatrix and Tnum --------------------------------------
      TYPE (basis)          :: basis_set
      integer, intent(in)   :: Set_Val

      integer :: nopt,ib,i,i1,i2
!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'Set_basis_para_FOR_optimization'
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
      END IF
!---------------------------------------------------------------------

      IF (basis_set%opt_param <= 0) RETURN
      IF (basis_set%nb_basis == 0) THEN

        ! for basis_set%A
        nopt = count(basis_set%opt_A /= 0)
        IF (debug) write(out_unitp,*) 'nopt A',nopt
        IF (nopt > 0) THEN
          i1 = para_FOR_optimization%i_OptParam+1
          i2 = para_FOR_optimization%i_OptParam+nopt
          para_FOR_optimization%nb_OptParam =                           &
                                para_FOR_optimization%nb_OptParam + nopt
          IF (Set_Val == -1) THEN
            DO i=1,nopt
              IF (basis_set%opt_A(i) /= 0)                              &
                para_FOR_optimization%Val_RVec(i1+i-1) = basis_set%A(i)
            END DO
          ELSE IF (Set_Val == 1) THEN
            DO i=1,nopt
              IF (basis_set%opt_A(i) /= 0)                              &
                basis_set%A(i) = para_FOR_optimization%Val_RVec(i1+i-1)
            END DO
          END IF
          para_FOR_optimization%i_OptParam = i2
        END IF

         ! for basis_set%B
        nopt = count(basis_set%opt_B /= 0)
        IF (debug) write(out_unitp,*) 'nopt B',nopt

        IF (nopt > 0) THEN
          i1 = para_FOR_optimization%i_OptParam+1
          i2 = para_FOR_optimization%i_OptParam+nopt
          para_FOR_optimization%nb_OptParam =                           &
                                para_FOR_optimization%nb_OptParam + nopt

          IF (Set_Val == -1) THEN
            DO i=1,nopt
              IF (basis_set%opt_B(i) /= 0)                              &
                para_FOR_optimization%Val_RVec(i1+i-1) = basis_set%B(i)
            END DO
          ELSE IF (Set_Val == 1) THEN
            DO i=1,nopt
              IF (basis_set%opt_B(i) /= 0)                              &
                basis_set%B(i) = para_FOR_optimization%Val_RVec(i1+i-1)
            END DO
          END IF
          para_FOR_optimization%i_OptParam = i2
        END IF

        ! for basis_set%Q0
        nopt = count(basis_set%opt_Q0 /= 0)
        IF (debug) write(out_unitp,*) 'nopt Q0',nopt
        IF (nopt > 0) THEN
          i1 = para_FOR_optimization%i_OptParam+1
          i2 = para_FOR_optimization%i_OptParam+nopt
          para_FOR_optimization%nb_OptParam =                           &
                                para_FOR_optimization%nb_OptParam + nopt

          IF (Set_Val == -1) THEN
            DO i=1,nopt
              IF (basis_set%opt_Q0(i) /= 0)                              &
                para_FOR_optimization%Val_RVec(i1+i-1) = basis_set%Q0(i)
            END DO
          ELSE IF (Set_Val == 1) THEN
            DO i=1,nopt
              IF (basis_set%opt_Q0(i) /= 0)                              &
                basis_set%Q0(i) = para_FOR_optimization%Val_RVec(i1+i-1)
            END DO
          END IF
          para_FOR_optimization%i_OptParam = i2
        END IF

        ! for basis_set%scaleQ
        nopt = count(basis_set%opt_scaleQ /= 0)
        IF (debug) write(out_unitp,*) 'nopt scaleQ',nopt

        IF (nopt > 0) THEN
          i1 = para_FOR_optimization%i_OptParam+1
          i2 = para_FOR_optimization%i_OptParam+nopt
          para_FOR_optimization%nb_OptParam =                           &
                                para_FOR_optimization%nb_OptParam + nopt

          IF (Set_Val == -1) THEN
            DO i=1,nopt
              IF (basis_set%opt_scaleQ(i) /= 0)                              &
                para_FOR_optimization%Val_RVec(i1+i-1) = basis_set%scaleQ(i)
            END DO
          ELSE IF (Set_Val == 1) THEN
            DO i=1,nopt
              IF (basis_set%opt_scaleQ(i) /= 0)                              &
                basis_set%scaleQ(i) = para_FOR_optimization%Val_RVec(i1+i-1)
            END DO
          END IF
          para_FOR_optimization%i_OptParam = i2

        END IF

      ELSE
        DO ib=1,size(basis_set%tab_Pbasis)
          CALL Set_basis_para_FOR_optimization(basis_set%tab_Pbasis(ib)%Pbasis,Set_Val)
        END DO
      END IF

      !-----------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'nb_OptParam ',para_FOR_optimization%nb_OptParam
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF

      END SUBROUTINE Set_basis_para_FOR_optimization

!================================================================
! ++    Construct a primitive basis set
!================================================================
      SUBROUTINE construct_primitive_basis(basis_primi)
      use mod_nDindex
      IMPLICIT NONE

!----- for the active basis set ---------------------------------------
      TYPE (basis)  :: basis_primi
      logical       :: deriv,nosym

      integer :: nb_shift,tab_shift(10),nstep,ic,nb_word,nq,nb,k
      character (len=Name_len) :: word

      TYPE (basis)  :: basis_temp


      real (kind=Rkind) :: weight

      logical :: err_grid

      integer :: L,nbi,nb0

!---------------------------------------------------------------------
      character (len=*), parameter :: name_sub='construct_primitive_basis'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL RecWrite_basis(basis_primi)
      END IF
!---------------------------------------------------------------------


      IF (.NOT. basis_primi%active) RETURN


!     - analyze the name -------------------------------
      CALL analysis_name(basis_primi%name,word,1,nb_word)
      IF (nb_word == 1) THEN
        nstep = 1
        nb_shift = 1
        tab_shift(1) = 0
        word = basis_primi%name
      ELSE IF (nb_word >= 3) THEN
        nb_shift = nb_word -2
        read(basis_primi%name,*) word,nstep,(tab_shift(ic),ic=1,nb_shift)
      ELSE
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' the analysis of the basis_primi%name has failed'
        write(out_unitp,*) ' basis:',basis_primi%name
        write(out_unitp,*) ' nb_word=',nb_word
        CALL RecWrite_basis(basis_primi)
        STOP
      END IF
      CALL string_uppercase_TO_lowercase(word)


      CALL basis2TObasis1(basis_temp,basis_primi)



      IF (basis_temp%L_SparseGrid > -1) THEN
        nq = Get_nq_FROM_l_OF_PrimBasis(basis_temp%L_SparseGrid,basis_temp)
        CALL Set_nq_OF_basis(basis_temp,nq)
      END IF

      SELECT CASE (word)
      CASE ("ylm","ylm_0a","ylm_1a")
        CONTINUE ! will be done in the sub_quadra... subroutine
      CASE ("cuba_ho","cuba_hm","cuba_hermite")
        CONTINUE ! will be done in the sub_quadra... subroutine
      CASE DEFAULT
        basis_temp%nb = Get_nb_FROM_l_OF_PrimBasis(basis_temp%L_SparseBasis,basis_temp)
      END SELECT



      SELECT CASE (word)


      CASE ("0")
        basis_temp%type = 0
        write(out_unitp,*) ' ERROR : in ',name_sub
        write(out_unitp,*) 'Basis "0" CANNOT be a primitive basis'
        STOP


      CASE ("direct_prod","sb","sg","sparse")
        basis_temp%type = 1
        write(out_unitp,*) ' ERROR : in ',name_sub
        write(out_unitp,*) 'Basis 1 (',trim(word),') CANNOT be a primitive basis'
        STOP

      CASE ("el")
        basis_temp%type = 2
        write(out_unitp,*) 'Basis "El": for several diabatic PES'
        CALL sub_basis_El(basis_temp)

      CASE ("pl0")
        basis_temp%type = 10
        CALL sub_quadra_legendre(basis_temp,-1)
      CASE ("pl0_0")
        basis_temp%type = 100
        CALL sub_quadra_legendre(basis_temp,0)
      CASE ("pl0_1")
        basis_temp%type = 101
        CALL sub_quadra_legendre(basis_temp,1)

      CASE ("pl0_a")
        basis_temp%type = 11
        CALL sub_quadra_legendre(basis_temp,-1)
        !CALL transfo_cosTOangle(basis_temp)
        CALL transfo_Q_TO_tQ(basis_temp)

      CASE ("pl0_a_0")
        basis_temp%type = 110
        CALL sub_quadra_legendre(basis_temp,0)
        CALL transfo_cosTOangle(basis_temp)
      CASE ("pl0_a_1")
        basis_temp%type = 111
        CALL sub_quadra_legendre(basis_temp,1)
        CALL transfo_cosTOangle(basis_temp)


      CASE ("ho","hm","hermite")
        basis_temp%type = 20

        IF (basis_temp%Nested == 2) THEN
          CALL sub_quadra_HermiteNested2(basis_temp,-1)
        ELSE IF (basis_temp%Nested <= 1) THEN
          CALL sub_quadra_hermite(basis_temp,-1)
        ELSE
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' hermite with the parameter value of Nested is not possible',basis_temp%Nested
          write(out_unitp,*) ' Nested',basis_temp%Nested
          STOP
        END IF
      CASE ("ho_0","hm_0")
        basis_temp%type = 200
        CALL sub_quadra_hermite(basis_temp,0)
      CASE ("ho_1","hm_1")
        basis_temp%type = 201
        CALL sub_quadra_hermite(basis_temp,1)
      CASE ("hobox","hmbox")
        basis_temp%type = 21
        CALL sub_quadra_hermitebox(basis_temp,-1)
      CASE ("cuba_ho","cuba_hm","cuba_hermite")
        basis_temp%type = 2000
        err_grid = .TRUE.
        CALL sub_quadra_hermite_cuba_DML(basis_temp,err_grid)
        IF (err_grid) CALL sub_quadra_hermite_cuba(basis_temp)

      CASE ("lm","laguerre")
        basis_temp%type = 90
        CALL sub_quadra_laguerre(basis_temp)

      CASE ("cos","sin","fourier")
        basis_temp%type = 30
        nosym = .FALSE.
        CALL sub_quadra_fourier(basis_temp,nosym,nstep,nb_shift,tab_shift)
      CASE ("cosab","cosabnosym" )
        basis_temp%type = 40
        nosym = word .EQ. "cosabnosym"
        CALL sub_quadra_fourier(basis_temp,nosym,nstep,nb_shift,tab_shift)

        IF (.NOT. basis_temp%xPOGridRep_done) THEN ! because the scaling factors are already calculated
          basis_temp%Q0     = (basis_temp%B+basis_temp%A)*HALF
          basis_temp%scaleQ = (pi+pi)/(basis_temp%B-basis_temp%A)
          write(out_unitp,*) 'A,B,Q0,scaleQ',basis_temp%A,basis_temp%B,         &
                                     basis_temp%Q0,basis_temp%scaleQ
        END IF
      CASE ("ft")
        basis_temp%type = 60
        CALL sub_quadra_FT(basis_temp,nosym,nstep,nb_shift,tab_shift)

        IF (.NOT. basis_temp%xPOGridRep_done) THEN ! because the scaling factors are already calculated
          basis_temp%Q0     = (basis_temp%B+basis_temp%A)*HALF
          basis_temp%scaleQ = (pi+pi)/(basis_temp%B-basis_temp%A)
          write(out_unitp,*) 'A,B,Q0,scaleQ',basis_temp%A,basis_temp%B,         &
                                     basis_temp%Q0,basis_temp%scaleQ
        END IF


      CASE ("boxab","boxabnosym")
        basis_temp%type = 50
        nosym = word .EQ. "boxabnosym"
        CALL sub_quadra_box(basis_temp,nosym)

        IF (.NOT. basis_temp%xPOGridRep_done) THEN ! because the scaling factors are already calculated
          basis_temp%Q0     = basis_temp%A
          basis_temp%scaleQ = pi/(basis_temp%B-basis_temp%A)
          write(out_unitp,*) 'A,B,Q0,scaleQ',basis_temp%A,basis_temp%B,         &
                                     basis_temp%Q0,basis_temp%scaleQ
        END IF


      CASE ("ylm")
        basis_temp%type = 60
        CALL sub_quadra_Ylm(basis_temp,-1,-1)
      CASE ("ylm_0a")
        basis_temp%type = 600
        CALL sub_quadra_Ylm(basis_temp,0,-1)
      CASE ("ylm_1a")
        basis_temp%type = 601
        CALL sub_quadra_Ylm(basis_temp,1,-1)
      CASE ("coll_abpluscd")
        basis_temp%type = 70
        CALL sub_quadra_ABplusCD(basis_temp)


      CASE DEFAULT
        write(out_unitp,*) ' This basis is unknown: ',word
        write(out_unitp,*) ' The possibilities :'
        write(out_unitp,*) '  0 : No basis!                           : 0'
        write(out_unitp,*)
        write(out_unitp,*) '  2 : Diabatic Electronic states          : El'
        write(out_unitp,*)
        write(out_unitp,*) '  1 : Direct product basis                : direct_prod'
        write(out_unitp,*) '  1 : Sparse basis and Smolyak Grids      : sparse or SB or SG'
        write(out_unitp,*)
        write(out_unitp,*) ' 10 : Legendre Poly.                      : Pl0'
        write(out_unitp,*) ' 100: Legendre Poly. (even)               : Pl0_0'
        write(out_unitp,*) ' 101: Legendre Poly. (odd)                : Pl0_1'
        write(out_unitp,*)
        write(out_unitp,*) ' 11 : Legendre Poly. theta                : Pl0_a'
        write(out_unitp,*) ' 110: Legendre Poly. theta (even)         : Pl0_a_0'
        write(out_unitp,*) ' 111: Legendre Poly. theta (odd)          : Pl0_a_1'
        write(out_unitp,*)
        write(out_unitp,*) ' 20 : Hermite Poly.                       : HO or Hm or hermite'
        write(out_unitp,*) ' 200: Hermite Poly. (even)                : HO_0 or Hm_0'
        write(out_unitp,*) ' 201: Hermite Poly. (odd)                 : HO_1 or Hm_1'
        write(out_unitp,*) ' 21 : PHermite Poly. + points from boxAB  : HObox or Hmbox'
        write(out_unitp,*) ' 2000:Hermite Poly. + cubature (nD)       : cuba_HO or cuba_Hm or cuba_hermite'

        write(out_unitp,*)
        write(out_unitp,*) ' 90 : Laguerre Poly.                      : Lm or Laguerre'

        write(out_unitp,*)
        write(out_unitp,*) ' 30 : Fourier Series                      : cos ou sin ou fourier'
        write(out_unitp,*) ' 40 : Fourier Series [A B]                : cosAB'
        write(out_unitp,*) ' 40 : Fourier Series [A B]                : cosABnosym'
        write(out_unitp,*) ' 60 : Fourier Transform [A B]             : FT'
        write(out_unitp,*) ' 50 : Particle-in-a-box[A B]              : boxAB'
        write(out_unitp,*) ' 50 : Particle-in-a-box[A B]              : boxABnosym'
        write(out_unitp,*)

        write(out_unitp,*) ' 60 : Ylm                                 : Ylm'
        write(out_unitp,*) ' 600: Ylm_0a                              : Ylm, even in l'
        write(out_unitp,*) ' 601: Ylm_1a                              : Ylm,  odd in l'
        write(out_unitp,*)

        write(out_unitp,*) ' 70 : Pl1m*Pl2m*Fourier                   : coll_ABplusCD'

        STOP
      END SELECT

      CALL basis2TObasis1(basis_primi,basis_temp)

      CALL dealloc_basis(basis_temp)

      nq = get_nq_FROM_basis(basis_primi)
      !write(6,*) 'basis_primi nq',nq
      IF (nq < 1) RETURN

      basis_primi%primitive      = .TRUE.
      basis_primi%primitive_done = .TRUE.

      CALL dealloc_nDindex(basis_primi%nDindG)
      CALL init_nDindexPrim(basis_primi%nDindG,1,(/ nq /))

      IF (basis_primi%type == 2000) THEN
        CONTINUE ! nD-HO with cubature: already done
      ELSE IF (basis_primi%type == 60 .OR. basis_primi%type == 600 .OR. basis_primi%type == 601) THEN
        CONTINUE ! Ylm: already done
      ELSE IF (basis_primi%type == 2) THEN
        CONTINUE ! El basis set: already done
      ELSE
        CALL dealloc_nDindex(basis_primi%nDindB)
        IF (basis_primi%With_L) THEN
          basis_primi%nDindB%packed = .TRUE.
          CALL init_nDindexPrim(basis_primi%nDindB,1,(/ basis_primi%nb /))
          basis_primi%nDindB%With_L = .TRUE.
          basis_primi%nDindB%Tab_L(:)    = -1
          basis_primi%nDindB%Tab_Norm(:) = -ONE

          DO L=basis_primi%L_SparseBasis,0,-1
            nb = Get_nb_FROM_l_OF_PrimBasis(L,basis_primi)
            basis_primi%nDindB%Tab_L(1:nb)    = L
            basis_primi%nDindB%Tab_Norm(1:nb) = real(L,kind=Rkind)
          END DO

          !CALL Write_nDindex(basis_primi%nDindB)

          !STOP 'With_L'
        ELSE
          weight = basis_primi%weight_OF_nDindB
          basis_primi%nDindB%packed = .TRUE.
          CALL init_nDindexPrim(basis_primi%nDindB,ndim=1,               &
                            Type_OF_nDindex=basis_primi%Type_OF_nDindB,  &
                            MaxNorm=basis_primi%Norm_OF_nDindB,          &
                       nb_OF_MinNorm=basis_primi%nb_OF_MinNorm_OF_nDindB,&
                     Div_nb_TO_Norm=basis_primi%Div_nb_TO_Norm_OF_nDindB,&
                            nDinit=(/ basis_primi%nDinit_OF_nDindB /),   &
                            nDsize=(/ basis_primi%nb /),                 &
                            nDweight=(/ weight /)      )
        END IF
      END IF

      CALL Basis_Grid_ParamTOBasis_Grid_Param_init(basis_primi%Basis_Grid_Para)
      basis_primi%nb_init = basis_primi%nb

      CALL pack_nDindex(basis_primi%nDindB)

    !- symmetry of the basis ---------------------------------
    IF (basis_primi%type == 60 .OR. basis_primi%type == 600 .OR. basis_primi%type == 601) THEN
      CONTINUE ! Ylm: already done
      ELSE IF (basis_primi%type == 2) THEN
        CONTINUE ! El basis set: already done
    ELSE
      CALL Set_tab_SymAbelian(basis_primi%P_SymAbelian,basis_primi%nb)
    END IF

      !- scaling of the basis ---------------------------------
      CALL sub_scale_basis(basis_primi)


      IF (basis_primi%xPOGridRep_done) RETURN
      CALL flush_perso(out_unitp)
!     - d1b => d1BasisRep and  d2b => d2BasisRep ------------
      CALL sub_dnGB_TO_dnBB(basis_primi)
      CALL flush_perso(out_unitp)
!     - d1b => dnBGG%d1 and  d2b => dnBGG%d2 ------------
      CALL sub_dnGB_TO_dnGG(basis_primi)
      !- d0b => transpose(d0b) ... transpose(d0bwrho) ---
      CALL sub_dnGB_TO_dnBG(basis_primi)

!     - check the overlap matrix -----------------------------
      CALL check_ortho_basis(basis_primi)

      IF (print_level > 2 .AND. allocated(basis_primi%x) ) THEN
        write(out_unitp,*) '---------------------------------------'
        write(out_unitp,*) 'x Grid:',basis_primi%x
        write(out_unitp,*) '---------------------------------------'
      END IF

!    DO k=1,size(basis_primi%x,dim=2)
!      write(666,'(101(x,f10.5))') basis_primi%x(:,k),basis_primi%dnRGB%d0(k,1:min(10,basis_primi%nb))
!    END DO

!---------------------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(basis_primi)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------
      END SUBROUTINE construct_primitive_basis

      FUNCTION d0b_OF_primitive_basis_AT_Q(basis_temp,ib,Qbasis) result(d0b)
      IMPLICIT NONE

!----- for the active basis set ---------------------------------------
      real (kind=Rkind)             :: d0b ! result of the function
      TYPE (basis), intent(inout)   :: basis_temp
      integer, intent(in)           :: ib
      real (kind=Rkind), intent(in) :: Qbasis(basis_temp%ndim)



      character (len=Name_len) :: word
      real (kind=Rkind)        :: poly_Hermite_exp ! function
      real (kind=Rkind)        :: x,xx


!---------------------------------------------------------------------
      character (len=*), parameter :: name_sub='d0b_OF_primitive_basis_AT_Q'
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL RecWrite_basis(basis_temp)
      END IF
!---------------------------------------------------------------------

      x = basis_temp%scaleQ(1)*(Qbasis(1)-basis_temp%Q0(1))


      word = basis_temp%name
      CALL string_uppercase_TO_lowercase(word)

      SELECT CASE (word)

      CASE ("pl0")
        STOP
      CASE ("pl0_0")
        STOP

      CASE ("pl0_1")
        STOP

      CASE ("pl0_a")
        STOP

      CASE ("pl0_a_0")
        STOP

      CASE ("pl0_a_1")
        STOP


      CASE ("ho","hm","hermite")

        d0b = sqrt(basis_temp%scaleQ(1))*poly_Hermite_exp(x,ib-1)
        !write(6,*) 'ib,x,d0b',ib,x,d0b
      CASE ("ho_0","hm_0")
        STOP

      CASE ("ho_1","hm_1")
         STOP

      CASE ("hobox","hmbox")
        STOP

      CASE ("cuba_ho","cuba_hm","cuba_hermite")
        STOP

      CASE ("lm","laguerre")
        STOP

      CASE ("cos","sin","fourier")
        STOP

      CASE ("cosab","cosabnosym" )
        STOP

      CASE ("ft")
        STOP

      CASE ("boxab","boxabnosym")

       xx = mod(x*real(ib,kind=Rkind),pi+pi)
       d0b = sin(xx) / sqrt(pi*HALF)

      CASE ("ylm")
        STOP

      CASE ("ylm_0a")
        STOP

      CASE ("ylm_1a")
        STOP

      CASE ("coll_abpluscd")
        STOP



      CASE DEFAULT

        STOP
      END SELECT


!---------------------------------------------------------------------
      IF (debug) THEN
        CALL RecWrite_basis(basis_temp)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------
      END FUNCTION d0b_OF_primitive_basis_AT_Q


!=============================================================
!
!      contraction of the basis set basis_set
!
!      1) read the new set of vector v10(nb_bc,nb_b)
!   or 2) use basis_set%Rvec (without_read = T)
!
!=============================================================
      !!@description: contraction of the basis set
      !!@param: basis_set          Basis set to be contracted
      !!@param: without_read  if this parameter is .TRUE., the contraction coeficients are already in basis_set,
      !!                      otherwise, the coeficients are read from a file
      SUBROUTINE sub_contraction_basis(basis_set,without_read)
      USE mod_system
      USE mod_dnSVM
      use mod_nDindex
      IMPLICIT NONE
!---------------------------------------------------------------------
!---------- variables passees en argument ----------------------------
      TYPE (basis), intent(inout)  :: basis_set
      logical,      intent(in)     :: without_read


!     - working parameters --------------------------------------------
      integer          :: i,j,k,l,nb_col
      integer          :: ic,inc,isym,jsym,iloc,jloc
      real(kind=Rkind) :: norm,Sij
      integer          :: nb_q,nb_qc,nb_b,nb_bc,ib,ibc
      integer          :: nq
      integer          :: nb_b1,nb_bc1
      integer          :: nio
      character (len=Name_longlen) :: nom

      real (kind=Rkind), allocatable :: Mat_read(:,:)  ! Mat_read(nb_b1,nb_bc1)
      TYPE (Type_dnMat) :: dnRGBuncontract
      TYPE (Type_dnMat) :: dnRBBuncontract

      integer, allocatable           :: tab_contract_symab(:)
      integer, allocatable           :: Tab_nDval(:,:)

      integer                    :: nderiv

      real (kind=Rkind)          :: weight_AT_ib

      integer :: err,err_sub
      integer :: Get_symabOFSymAbelianOFBasis_AT_ib ! function

!---------------------------------------------------------------------

!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sub_contraction_basis'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        CALL RecWrite_basis(basis_set,.TRUE.)
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------


      basis_set%packed            = .TRUE.
      basis_set%packed_done       = .TRUE.


      IF (basis_set%cplx) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' I cannot contract a COMPLEX basis set !'
        write(out_unitp,*) ' Not yet implemented'
        STOP
      END IF


      nb_bc = basis_set%nbc
      nb_b  = basis_set%nb

      IF (debug .OR. print_level > -1) write(out_unitp,*) 'nb',basis_set%nb
      IF (debug .OR. print_level > -1) write(out_unitp,*) 'nbc',basis_set%nbc

      IF (.NOT. without_read) THEN
        IF (allocated(basis_set%Rvec))  THEN
          CALL dealloc_NParray(basis_set%Rvec,"basis_set%Rvec",name_sub)
        END IF
        CALL alloc_NParray(basis_set%Rvec,(/ nb_b,nb_b /),"basis_set%Rvec",name_sub)
        basis_set%Rvec(:,:) = ZERO
!       - test and read a matrix ------------------------
        nio = 5
        IF (basis_set%read_contrac_file) CALL file_open(basis_set%file_contrac,nio)

        read(nio,*)
        read(nio,*) nb_col,nb_bc1,nb_b1
        IF (debug .OR. print_level > -1) THEN
          write(out_unitp,*) 'nb_col,nb_bc1,nb_b1',nb_col,nb_bc1,nb_b1
          write(out_unitp,*) '    The basis is contracted'
          write(out_unitp,*) '      nbc nbc_lect: ',nb_bc,nb_bc1
          write(out_unitp,*) '      nb: ',nb_b1
        END IF

        IF (nb_b1 /= nb_b .OR. nb_bc1 < nb_bc) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' nb_bc =',nb_bc,' and the matrice v10 has',       &
                       nb_bc1,'columns'
          write(out_unitp,*) ' nb_b =',nb_b,' and the matrice v10 has',         &
                       nb_b1,'lignes'
          write(out_unitp,*) ' Rq: You MUST have nb_b1=nb_b and nb_bc1=>nb_bc'
          STOP
        END IF

!       - read the matrix for the contraction ----------
!       ------------------ Mat_read --------------------
!
        CALL alloc_NParray(Mat_read,(/ nb_b1,nb_bc1 /),"Mat_read",name_sub)
        CALL Read_Mat(Mat_read,nio,nb_col,err)
        IF (err /= 0) THEN
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) ' reading the matrix "Mat_read"'
          STOP
        END IF
        IF (debug) THEN
          write(out_unitp,*) ' contraction basis'
          write(out_unitp,*)
          CALL Write_Mat(Mat_read,out_unitp,nb_col)
          CALL flush_perso(out_unitp)
        END IF
!       -------------------------------------------------
        IF (basis_set%read_contrac_file) close(nio)
        basis_set%Rvec(:,1:nb_bc1) = Mat_read(:,:)
        CALL dealloc_NParray(Mat_read,"Mat_read",name_sub)
      ELSE
        IF (.NOT. allocated(basis_set%Rvec)) THEN
          write(out_unitp,*) ' ERROR in ',name_sub
          write(out_unitp,*) ' without_read=t and basis_set%Rvec is not allocated!'
          write(out_unitp,*) ' CHECK the source!!'
          STOP
        END IF
      END IF

      !-------------------------------------------------
      !- first the symmetry because Rvec can be modified
      CALL alloc_NParray(tab_contract_symab,(/ basis_set%nbc /),          &
                        "tab_contract_symab",name_sub)
      !analyze the symmetry of each Rvec(:,i)
      !write(out_unitp,*) 'basis_set%tab_symab',basis_set%tab_symab
      DO i=1,nb_bc
        ! symmetry of the basis function with the largest coefficient
        iloc = sum(maxloc(abs(basis_set%Rvec(:,i))))
        isym = Get_symabOFSymAbelianOFBasis_AT_ib(basis_set,iloc)
        tab_contract_symab(i) = isym ! symmetry of Rvec(:,i)

        IF (isym == -1) CYCLE ! the symetry is not used

        DO j=1,nb_b
          jsym = Get_symabOFSymAbelianOFBasis_AT_ib(basis_set,j)
          IF (jsym /= isym) basis_set%Rvec(j,i) = ZERO
        END DO

        ! schmidt orthogonalization
        DO j=1,i-1
          ! symmetry of the basis function with the largest coefficient
          jloc = sum(maxloc(abs(basis_set%Rvec(:,j))))
          jsym = Get_symabOFSymAbelianOFBasis_AT_ib(basis_set,jloc)
          IF (isym == jsym) THEN
             Sij = dot_product(basis_set%Rvec(:,i),basis_set%Rvec(:,j))
             !write(out_unitp,*) 'i,j,Sij',i,isym,j,jsym,Sij
             basis_set%Rvec(:,i) = basis_set%Rvec(:,i) - Sij*basis_set%Rvec(:,j)
          END IF
        END DO
        norm = dot_product(basis_set%Rvec(:,i),basis_set%Rvec(:,i))
        IF ((ONE-norm) > ONETENTH**6)                                    &
          write(out_unitp,*) 'WARNING the contracted coefficients break the symmetry',norm
        basis_set%Rvec(:,i) = basis_set%Rvec(:,i)/sqrt(norm)

        ! Once again schmidt orthogonalization
        DO j=1,i-1
          ! symmetry of the basis function with le largest coeficient
          jloc = sum(maxloc(abs(basis_set%Rvec(:,j))))
          jsym = Get_symabOFSymAbelianOFBasis_AT_ib(basis_set,jloc)
          IF (isym == jsym) THEN
             Sij = dot_product(basis_set%Rvec(:,i),basis_set%Rvec(:,j))
             !write(out_unitp,*) 'i,j,Sij',i,j,Sij
             basis_set%Rvec(:,i) = basis_set%Rvec(:,i) - Sij*basis_set%Rvec(:,j)
          END IF
        END DO
        norm = dot_product(basis_set%Rvec(:,i),basis_set%Rvec(:,i))
        basis_set%Rvec(:,i) = basis_set%Rvec(:,i)/sqrt(norm)
      END DO


      CALL Set_tab_symabOFSymAbelian_WITH_tab(basis_set%P_SymAbelian,     &
                                              tab_contract_symab)

      CALL dealloc_NParray(tab_contract_symab,"tab_contract_symab",name_sub)
      !-------------------------------------------------

      !-------------------------------------------------
      !-- secondly, the basis d0b,d1b,d2b or d0b, dnRBB dnCBB
      !write(out_unitp,*) ' matmul contra'
      nq = get_nq_FROM_basis(basis_set)

      nderiv = 0
      IF (associated(basis_set%dnRGB%d1)) nderiv = 1
      IF (associated(basis_set%dnRGB%d2)) nderiv = 2

      dnRGBuncontract = basis_set%dnRGB
      CALL dealloc_dnMat(basis_set%dnRGB)


      CALL alloc_dnMat(basis_set%dnRGB,nq,basis_set%nbc,basis_set%ndim,nderiv)

      basis_set%dnRGB%d0(:,:) =  matmul(dnRGBuncontract%d0(:,:),        &
                                             basis_set%Rvec(:,1:nb_bc))

      IF (nderiv > 0) THEN
        DO i=1,basis_set%ndim
          basis_set%dnRGB%d1(:,:,i) =                             &
                                matmul(dnRGBuncontract%d1(:,:,i), &
                                             basis_set%Rvec(:,1:nb_bc))
        END DO
      END IF

      IF (nderiv > 1) THEN
        DO i=1,basis_set%ndim
        DO j=1,basis_set%ndim
          basis_set%dnRGB%d2(:,:,i,j) =                           &
                              matmul(dnRGBuncontract%d2(:,:,i,j), &
                                             basis_set%Rvec(:,1:nb_bc))
        END DO
        END DO
      END IF

      CALL dealloc_dnMat(dnRGBuncontract)


      IF (basis_set%dnBBRep_done) THEN

        dnRBBuncontract = basis_set%dnRBB
        CALL dealloc_dnMat(basis_set%dnRBB)
        CALL alloc_dnMat(basis_set%dnRBB,basis_set%nbc,basis_set%nbc,basis_set%ndim,nderiv=2)


        DO k=1,basis_set%nbc
        DO l=1,basis_set%nbc

          basis_set%dnRBB%d1(k,l,:)   = ZERO
          basis_set%dnRBB%d2(k,l,:,:) = ZERO

          DO i=1,basis_set%nb
          DO j=1,basis_set%nb

            basis_set%dnRBB%d1(k,l,:)   = basis_set%dnRBB%d1(k,l,:) +   &
              basis_set%Rvec(i,k)*basis_set%Rvec(j,l) * dnRBBuncontract%d1(i,j,:)
            basis_set%dnRBB%d2(k,l,:,:) = basis_set%dnRBB%d2(k,l,:,:) + &
              basis_set%Rvec(i,k)*basis_set%Rvec(j,l) * dnRBBuncontract%d2(i,j,:,:)

          END DO
          END DO
        END DO
        END DO

        IF (print_level > -1)  THEN
          CALL Write_dnMat(basis_set%dnRBB,nderiv=1)
        END IF

      END IF
      !-------------------------------------------------

      !-------------------------------------------------
      !- Save the uncontracted nDindB
      IF (associated(basis_set%nDindB_uncontracted)) THEN
        CALL dealloc_nDindex(basis_set%nDindB_uncontracted)
        CALL dealloc_array(basis_set%nDindB_uncontracted,               &
                          'basis_set%nDindB_uncontracted',name_sub)
      END IF
      CALL alloc_array(basis_set%nDindB_uncontracted,                   &
                      'basis_set%nDindB_uncontracted',name_sub)
      CALL nDindex2TOnDindex1(basis_set%nDindB_uncontracted,basis_set%nDindB)

      !- Finally, nDindB + the new Tab_Norm
      IF (basis_set%contrac_WITH_nDindB) THEN

        CALL dealloc_nDindex(basis_set%nDindB)
        basis_set%nDindB%packed      = .TRUE.
        CALL init_nDindexPrim(basis_set%nDindB,ndim=1,                 &
                            Type_OF_nDindex=basis_set%Type_OF_nDindB,  &
                            MaxNorm=basis_set%Norm_OF_nDindB,          &
                       nb_OF_MinNorm=basis_set%nb_OF_MinNorm_OF_nDindB,&
                     Div_nb_TO_Norm=basis_set%Div_nb_TO_Norm_OF_nDindB,&
                            nDinit=(/ basis_set%nDinit_OF_nDindB /),   &
                            nDsize=(/ basis_set%nb /),                 &
                            nDweight=(/ basis_set%weight_OF_nDindB /) )
        CALL pack_nDindex(basis_set%nDindB)

      ELSE

        CALL alloc_NParray(Tab_nDval,(/ 1,basis_set%nbc /),"Tab_nDval",name_sub)

        CALL dealloc_nDindex(basis_set%nDindB)
        CALL alloc_NParray(basis_set%nDindB%Tab_Norm,(/basis_set%nbc/),   &
                          "basis_set%nDindB%Tab_Norm",name_sub)
        ! calculation of the Norm of each contracted vector
        DO ibc=1,basis_set%nbc
          Norm = ZERO
          DO ib=1,basis_set%nb
            Norm = Norm + abs(basis_set%Rvec(ib,ibc))**2 *              &
                      calc_Norm_OF_nDI(basis_set%nDindB_uncontracted,ib)
          END DO
          basis_set%nDindB%Tab_Norm(ibc)    = Norm
          Tab_nDval(1,ibc) = ibc
        END DO


        basis_set%nDindB%Tab_Norm(:) = basis_set%nDindB%Tab_Norm(:) -   &
                                       minval(basis_set%nDindB%Tab_Norm)


        CALL init_nDindex_typeTAB(basis_set%nDindB,1,Tab_nDval,basis_set%nbc,err_sub)
        IF (err_sub /= 0) THEN
          write(*,*) ' ERROR in ',name_sub
          STOP ' from init_nDindex_typeTAB'
        END IF
!        basis_set%nDindB%packed      = .TRUE.
!        CALL pack_nDindex(basis_set%nDindB)

        CALL dealloc_NParray(Tab_nDval,"Tab_nDval",name_sub)

      END IF
      !-------------------------------------------------

      IF (allocated(basis_set%tab_ndim_index))  THEN
        CALL dealloc_NParray(basis_set%tab_ndim_index,                    &
                            "basis_set%tab_ndim_index",name_sub)
      END IF
      CALL alloc_NParray(basis_set%tab_ndim_index,                        &
                                    (/ basis_set%ndim,basis_set%nbc /), &
                      "basis_set%tab_ndim_index",name_sub)
      basis_set%tab_ndim_index(:,:) = 0
      basis_set%tab_ndim_index(1,:) = (/ (ibc,ibc=1,basis_set%nbc) /)

      basis_set%nb  = basis_set%nbc

      IF (.NOT. basis_set%contrac_WITH_nDindB) CALL sort_basis(basis_set)
      IF (print_level > -1) THEN
        write(out_unitp,*) 'tab_ndim_index contr',basis_set%tab_ndim_index(1,:)

        write(out_unitp,*) 'Tab_Norm(:) of Contracted basis',           &
                                            basis_set%nDindB%Tab_Norm(:)
      END IF

      !-------------------------------------------------
      !- d1b => dnGG%d1 and  d2b => dnGG%d2 ----------
      CALL sub_dnGB_TO_dnGG(basis_set)
      !-------------------------------------------------

      !-------------------------------------------------
      !- d0b => transpose(d0b) ... traspose(d0bwrho) ---
      !CALL sub_dnGB_TO_dnBG(basis_set)
      !-------------------------------------------------

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        CALL RecWrite_basis(basis_set)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE sub_contraction_basis

      !!@description: contraction of the basis set
      !!@param: basis_set          Basis set to be contracted
      !!@param: without_read  if this parameter is .TRUE., the contraction coeficients are already in basis_set,
      !!                      otherwise, the coeficients are read from a file
      SUBROUTINE sort_basis(basis_set)
      use mod_nDindex
      IMPLICIT NONE
      TYPE (basis), intent(inout) :: basis_set


!     - working parameters --------------------------------------------
      integer                          :: ib,jb,isym,jsym,i,j,nb_uncc
      integer                          :: nq,Li

      integer, allocatable             :: nDval(:)
      real(kind=Rkind)                 :: Normi,Normj
      real(kind=Rkind), allocatable    :: bi(:),biF(:),bcci(:)
      complex(kind=Rkind), allocatable :: cbi(:),cbiF(:)
      integer, allocatable             :: ni(:)

      integer :: Get_symabOFSymAbelianOFBasis_AT_ib ! function


!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='sort_basis'
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        !CALL RecWrite_basis(basis_set,.TRUE.)
        write(out_unitp,*)
        write(out_unitp,*) 'sort basis set, nb,iQdyn',basis_set%nb,basis_set%iQdyn(:)
        write(out_unitp,*) 'type_OF_nDindex',basis_set%nDindB%type_OF_nDindex
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------

     IF (.NOT. basis_set%nDindB%packed_done .OR. .NOT. basis_set%packed_done) THEN
        CALL RecWrite_basis(basis_set,.TRUE.)
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' I cannot sort the basis set: '
        write(out_unitp,*) ' The multidimensional index (nDindB) is unpacked'
        write(out_unitp,*) ' or the basis_set is unpacked'

        write(out_unitp,*) ' CHECK the fortran source !!'
        STOP
      END IF
      nq = get_nq_FROM_basis(basis_set)

      CALL alloc_NParray(nDval,(/ basis_set%nDindB%ndim /),'nDval',name_sub)
      CALL alloc_NParray(bi,   (/ nq /),                   'bi',name_sub)
      CALL alloc_NParray(cbi,  (/ nq /),                   'cbi',name_sub)
      CALL alloc_NParray(biF,  (/ basis_set%nb /),         'biF',name_sub)
      CALL alloc_NParray(cbiF, (/ basis_set%nb /),         'cbiF',name_sub)
      CALL alloc_NParray(ni,   (/ basis_set%ndim /),       'ni',name_sub)

      IF (debug) THEN
        write(out_unitp,*) ' Unsorted basis set with Norm of nDindB'
        DO ib=1,basis_set%nb
          Normi = calc_Norm_OF_nDI(basis_set%nDindB,ib)
          CALL calc_nDindex(basis_set%nDindB,ib,nDval)
          write(out_unitp,*) 'ib,nDval,Normi',ib,':',nDval,Normi
        END DO
      END IF

      !write(out_unitp,*) ' coucou sorted basis'

      !write(out_unitp,*) 'Permutation'
      DO ib=1,basis_set%nb
        DO jb=ib+1,basis_set%nb
          Normi = calc_Norm_OF_nDI(basis_set%nDindB,ib)
          Normj = calc_Norm_OF_nDI(basis_set%nDindB,jb)


          IF (Normi > Normj) THEN   ! permutation: d0b, d1b, d2b ...

            IF (allocated(basis_set%Rvec)) THEN
              nb_uncc = ubound(basis_set%Rvec,dim=1)
              CALL alloc_NParray(bcci,(/ nb_uncc /),'bcci',name_sub)

              bcci(:)              = basis_set%Rvec(:,ib)
              basis_set%Rvec(:,ib) = basis_set%Rvec(:,jb)
              basis_set%Rvec(:,jb) = bcci(:)

              CALL dealloc_NParray(bcci,'bcci',name_sub)
            END IF

            !write(out_unitp,*) 'permutation',ib,jb
            ! first the basis coef : d0b, d1b, d2b
            IF (allocated(basis_set%tab_ndim_index)) THEN
              ni(:)                          = basis_set%tab_ndim_index(:,ib)
              basis_set%tab_ndim_index(:,ib) = basis_set%tab_ndim_index(:,jb)
              basis_set%tab_ndim_index(:,jb) = ni(:)
            END IF
            IF (associated(basis_set%dnRGB%d0)) THEN
              bi(:)                    = basis_set%dnRGB%d0(:,ib)
              basis_set%dnRGB%d0(:,ib) = basis_set%dnRGB%d0(:,jb)
              basis_set%dnRGB%d0(:,jb) = bi(:)
            END IF
            IF (associated(basis_set%dnCGB%d0)) THEN
              cbi(:)                   = basis_set%dnCGB%d0(:,ib)
              basis_set%dnCGB%d0(:,ib) = basis_set%dnCGB%d0(:,jb)
              basis_set%dnCGB%d0(:,jb) = cbi(:)
            END IF

            IF (associated(basis_set%dnRGB%d1)) THEN
              DO i=1,basis_set%ndim
                bi(:)                      = basis_set%dnRGB%d1(:,ib,i)
                basis_set%dnRGB%d1(:,ib,i) = basis_set%dnRGB%d1(:,jb,i)
                basis_set%dnRGB%d1(:,jb,i) = bi(:)
              END DO
            END IF
           IF (basis_set%dnRBB%alloc) THEN
              DO i=1,basis_set%ndim
                biF(:)                     = basis_set%dnRBB%d1(:,ib,i)
                basis_set%dnRBB%d1(:,ib,i) = basis_set%dnRBB%d1(:,jb,i)
                basis_set%dnRBB%d1(:,jb,i) = biF(:)
                biF(:)                     = basis_set%dnRBB%d1(ib,:,i)
                basis_set%dnRBB%d1(ib,:,i) = basis_set%dnRBB%d1(jb,:,i)
                basis_set%dnRBB%d1(jb,:,i) = biF(:)
              END DO
            END IF
            IF (associated(basis_set%dnCGB%d1)) THEN
              DO i=1,basis_set%ndim
                cbi(:)                     = basis_set%dnCGB%d1(:,ib,i)
                basis_set%dnCGB%d1(:,ib,i) = basis_set%dnCGB%d1(:,jb,i)
                basis_set%dnCGB%d1(:,jb,i) = cbi(:)
              END DO
            END IF
           IF (basis_set%dnCBB%alloc) THEN
              DO i=1,basis_set%ndim
                cbiF(:)                    = basis_set%dnCBB%d1(:,ib,i)
                basis_set%dnCBB%d1(:,ib,i) = basis_set%dnCBB%d1(:,jb,i)
                basis_set%dnCBB%d1(:,jb,i) = cbiF(:)
                cbiF(:)                    = basis_set%dnCBB%d1(ib,:,i)
                basis_set%dnCBB%d1(ib,:,i) = basis_set%dnCBB%d1(jb,:,i)
                basis_set%dnCBB%d1(jb,:,i) = cbiF(:)
              END DO
            END IF

            IF (associated(basis_set%dnRGB%d2)) THEN
              DO i=1,basis_set%ndim
              DO j=1,basis_set%ndim
                bi(:)                        = basis_set%dnRGB%d2(:,ib,i,j)
                basis_set%dnRGB%d2(:,ib,i,j) = basis_set%dnRGB%d2(:,jb,i,j)
                basis_set%dnRGB%d2(:,jb,i,j) = bi(:)
              END DO
              END DO
            END IF
           IF (basis_set%dnRBB%alloc) THEN
              DO i=1,basis_set%ndim
              DO j=1,basis_set%ndim
                biF(:)                       = basis_set%dnRBB%d2(:,ib,i,j)
                basis_set%dnRBB%d2(:,ib,i,j) = basis_set%dnRBB%d2(:,jb,i,j)
                basis_set%dnRBB%d2(:,jb,i,j) = biF(:)
                biF(:)                       = basis_set%dnRBB%d2(ib,:,i,j)
                basis_set%dnRBB%d2(ib,:,i,j) = basis_set%dnRBB%d2(jb,:,i,j)
                basis_set%dnRBB%d2(jb,:,i,j) = biF(:)
              END DO
              END DO
            END IF
            IF (associated(basis_set%dnCGB%d2)) THEN
              DO i=1,basis_set%ndim
              DO j=1,basis_set%ndim
                cbi(:)                       = basis_set%dnCGB%d2(:,ib,i,j)
                basis_set%dnCGB%d2(:,ib,i,j) = basis_set%dnCGB%d2(:,jb,i,j)
                basis_set%dnCGB%d2(:,jb,i,j) = cbi(:)
              END DO
              END DO
            END IF
           IF (basis_set%dnCBB%alloc) THEN
              DO i=1,basis_set%ndim
              DO j=1,basis_set%ndim
                cbiF(:)                      = basis_set%dnCBB%d2(:,ib,i,j)
                basis_set%dnCBB%d2(:,ib,i,j) = basis_set%dnCBB%d2(:,jb,i,j)
                basis_set%dnCBB%d2(:,jb,i,j) = cbiF(:)
                cbiF(:)                      = basis_set%dnCBB%d2(ib,:,i,j)
                basis_set%dnCBB%d2(ib,:,i,j) = basis_set%dnCBB%d2(jb,:,i,j)
                basis_set%dnCBB%d2(jb,:,i,j) = cbiF(:)
              END DO
              END DO
            END IF

            ! permutation of the elements of tab_symab
            IF (SymAbelian_IS_initialized(basis_set%P_SymAbelian)) THEN
              isym = Get_symabOFSymAbelianOFBasis_AT_ib(basis_set,ib)
              jsym = Get_symabOFSymAbelianOFBasis_AT_ib(basis_set,jb)

              CALL Set_symabOFSymAbelian_AT_ib(basis_set%P_SymAbelian,ib,jsym)
              CALL Set_symabOFSymAbelian_AT_ib(basis_set%P_SymAbelian,jb,isym)
            END IF

            ! Permutation of tab_val of nDindB
            IF (basis_set%nDindB%packed) THEN
              nDval(:)                         = basis_set%nDindB%Tab_nDval(:,ib)
              basis_set%nDindB%Tab_nDval(:,ib) = basis_set%nDindB%Tab_nDval(:,jb)
              basis_set%nDindB%Tab_nDval(:,jb) = nDval(:)

              Normi                            = basis_set%nDindB%Tab_Norm(ib)
              basis_set%nDindB%Tab_Norm(ib)    = basis_set%nDindB%Tab_Norm(jb)
              basis_set%nDindB%Tab_Norm(jb)    = Normi

              Li                               = basis_set%nDindB%Tab_L(ib)
              basis_set%nDindB%Tab_L(ib)       = basis_set%nDindB%Tab_L(jb)
              basis_set%nDindB%Tab_L(jb)       = Li
            END IF


          END IF


        END DO
      END DO
      !write(out_unitp,*) 'END Permutation'

      CALL sub_dnGB_TO_dnBG(basis_set)

      IF (debug) THEN
        write(out_unitp,*) ' Sorted basis set with Norm of nDindB'
        DO ib=1,basis_set%nb
          Normi = calc_Norm_OF_nDI(basis_set%nDindB,ib)
          CALL calc_nDindex(basis_set%nDindB,ib,nDval)
          write(out_unitp,*) 'ib,nDval,Normi',ib,':',nDval,Normi
        END DO
      END IF


      CALL dealloc_NParray(nDval,'nDval',name_sub)
      CALL dealloc_NParray(bi,   'bi',   name_sub)
      CALL dealloc_NParray(cbi,  'cbi',  name_sub)
      CALL dealloc_NParray(biF,  'biF',  name_sub)
      CALL dealloc_NParray(cbiF, 'cbiF', name_sub)
      CALL dealloc_NParray(ni,   'ni',   name_sub)

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        !CALL RecWrite_basis(basis_set)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE sort_basis

!=============================================================
!
!      scale de la basis_set et des points (poids) de la quadrature
!
!      Qnew = Q0new + Qquadra / scale
!
!=============================================================
      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      SUBROUTINE sub_scale_basis(basis_sc)

      TYPE (basis)  :: basis_sc

!---------- working variables ----------------------------------------
      logical       ::  not_scaled

      real (kind=Rkind) ::  scaleQ,scale_inv
      real (kind=Rkind) ::  scale_d0b
      real (kind=Rkind) ::  scale_d1b(basis_sc%ndim)
      real (kind=Rkind) ::  scale_d2b(basis_sc%ndim,basis_sc%ndim)

      integer           :: i,j

!---------------------------------------------------------------------
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING sub_scale_basis'
        !CALL RecWrite_basis(basis_sc)
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------


      not_scaled = (basis_sc%type == 1 .OR. basis_sc%ndim == 0) ! direct_product
      DO i=1,basis_sc%ndim
        not_scaled  =  not_scaled .AND. basis_sc%Q0(i)     == ZERO
        not_scaled  =  not_scaled .AND. basis_sc%scaleQ(i) == ONE
      END DO


      IF ( .NOT. not_scaled ) THEN

        IF (basis_sc%print_info_OF_basisDP .AND. print_level > -1) THEN
          write(out_unitp,*) '    The basis is scaled'
          write(out_unitp,*) '          Q0(:) ',basis_sc%Q0(:)
          write(out_unitp,*) '      scaleQ(:) ',basis_sc%scaleQ(:)
        END IF

        scaleQ = product(basis_sc%scaleQ)

        scale_inv = ONE / scaleQ
        scale_d0b    = sqrt(scaleQ)
        DO i=1,basis_sc%ndim
          scale_d1b(i) = sqrt(scaleQ) * basis_sc%scaleQ(i)
        END DO
        DO i=1,basis_sc%ndim
        DO j=1,basis_sc%ndim
          scale_d2b(i,j) = sqrt(scaleQ) * basis_sc%scaleQ(i) * basis_sc%scaleQ(j)
        END DO
        END DO

        !write(6,*) 'alloc + shape x',allocated(basis_sc%x),shape(basis_sc%x)
        !write(6,*) 'nq',get_nq_FROM_basis(basis_sc)
        !flush(6)
        DO i=1,basis_sc%ndim
          basis_sc%x(i,:) = basis_sc%Q0(i) + basis_sc%x(i,:) / basis_sc%scaleQ(i)
        END DO

        basis_sc%w(:)            = basis_sc%w(:)          * scale_inv
        basis_sc%wrho(:)         = basis_sc%wrho(:)       * scale_inv

        IF (basis_sc%cplx) THEN
          basis_sc%dnCGB%d0(:,:)       = basis_sc%dnCGB%d0(:,:)     * cmplx(scale_d0b,kind=Rkind)
          DO i=1,basis_sc%ndim
            basis_sc%dnCGB%d1(:,:,i)   = basis_sc%dnCGB%d1(:,:,i)   * cmplx(scale_d1b(i),kind=Rkind)
          END DO
          DO i=1,basis_sc%ndim
          DO j=1,basis_sc%ndim
            basis_sc%dnCGB%d2(:,:,i,j) = basis_sc%dnCGB%d2(:,:,i,j) * cmplx(scale_d2b(i,j),kind=Rkind)
          END DO
          END DO
        ELSE
          basis_sc%dnRGB%d0(:,:)        = basis_sc%dnRGB%d0(:,:)      * scale_d0b
          DO i=1,basis_sc%ndim
            basis_sc%dnRGB%d1(:,:,i)    = basis_sc%dnRGB%d1(:,:,i)    * scale_d1b(i)
          END DO
          DO i=1,basis_sc%ndim
          DO j=1,basis_sc%ndim
            basis_sc%dnRGB%d2(:,:,i,j)  = basis_sc%dnRGB%d2(:,:,i,j)  * scale_d2b(i,j)
          END DO
          END DO
        END IF

      END IF



!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        !write(out_unitp,*) 'basis_sc%x',basis_sc%x
        CALL RecWrite_basis(basis_sc)
        CALL check_ortho_basis(basis_sc)
        write(out_unitp,*) 'END sub_scale_basis'
      END IF
!---------------------------------------------------------------------


      END SUBROUTINE sub_scale_basis

      SUBROUTINE sub_dnGB_TO_dnBB(basis_set)
      use mod_dnSVM
      IMPLICIT NONE

      TYPE (basis)  :: basis_set
!---------------------------------------------------------------------
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      integer           :: i,j,iq,ib,jb,nq
      real(kind=Rkind)  :: wrho
      real (kind=Rkind) :: d0,d1(basis_set%ndim)
      real (kind=Rkind) :: d2(basis_set%ndim,basis_set%ndim)
      real(kind=Rkind), allocatable :: d0bwrho(:,:)

      complex (kind=Rkind) :: d0c,d1c(basis_set%ndim)
      complex (kind=Rkind) :: d2c(basis_set%ndim,basis_set%ndim)
      complex (kind=Rkind), allocatable :: d0cbwrho(:,:)

!---------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug=.FALSE.
!      logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='sub_dnGB_TO_dnBB'
!-----------------------------------------------------------
      IF (basis_set%ndim == 0) RETURN

      nq = get_nq_FROM_basis(basis_set)
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb,nq',basis_set%nb,nq
      END IF
!-----------------------------------------------------------

      IF (.NOT. basis_set%packed_done .OR. .NOT. basis_set%dnBBRep) RETURN

      IF (basis_set%cplx) THEN
         CALL alloc_dnCplxMat(basis_set%dnCBB,basis_set%nb,basis_set%nb,nb_var_deriv=basis_set%ndim,nderiv=2)

        CALL alloc_NParray(d0cbwrho,(/ nq,basis_set%nb /),'d0cbwrho',name_sub)

        DO iq=1,nq
          wrho = Rec_WrhonD(basis_set,iq)
          DO ib=1,basis_set%nb
            d0cbwrho(iq,ib) = Rec_d0cbnD(basis_set,iq,ib) * wrho
          END DO
        END DO

        basis_set%dnCBB%d1(:,:,:)   = CZERO
        basis_set%dnCBB%d2(:,:,:,:) = CZERO

        DO ib=1,basis_set%nb
        DO iq=1,nq

          CALL Rec_d0d1d2cbnD(d0c,d1c(:),d2c(:,:),basis_set,iq,ib)
          DO jb=1,basis_set%nb

            basis_set%dnCBB%d1(jb,ib,:)   = basis_set%dnCBB%d1(jb,ib,:) + &
                                        d0cbwrho(iq,jb) * d1c(:)
            basis_set%dnCBB%d2(jb,ib,:,:) = basis_set%dnCBB%d2(jb,ib,:,:)+&
                                      d0cbwrho(iq,jb) * d2c(:,:)

          END DO
        END DO
        END DO
        CALL dealloc_NParray(d0cbwrho,'d0cbwrho',name_sub)

      ELSE
        CALL alloc_dnMat(basis_set%dnRBB,basis_set%nb,basis_set%nb,nb_var_deriv=basis_set%ndim,nderiv=2)

        CALL alloc_NParray(d0bwrho,(/ nq,basis_set%nb /),'d0bwrho',name_sub)


        DO iq=1,nq
          wrho = Rec_WrhonD(basis_set,iq)
          DO ib=1,basis_set%nb
            d0bwrho(iq,ib) = Rec_d0bnD(basis_set,iq,ib) * wrho
          END DO
        END DO

        basis_set%dnRBB%d1(:,:,:)   = ZERO
        basis_set%dnRBB%d2(:,:,:,:) = ZERO

        DO ib=1,basis_set%nb
        DO iq=1,nq

          CALL Rec_d0d1d2bnD(d0,d1(:),d2(:,:),basis_set,iq,ib)
          DO jb=1,basis_set%nb

            basis_set%dnRBB%d1(jb,ib,:)   = basis_set%dnRBB%d1(jb,ib,:)   + &
                                               d0bwrho(iq,jb) * d1(:)
            basis_set%dnRBB%d2(jb,ib,:,:) = basis_set%dnRBB%d2(jb,ib,:,:) + &
                                               d0bwrho(iq,jb) * d2(:,:)

          END DO
        END DO
        END DO
        CALL dealloc_NParray(d0bwrho,'d0bwrho',name_sub)

      END IF
      basis_set%dnBBRep_done = .TRUE.


!-----------------------------------------------------------
      IF (debug) THEN
        IF (basis_set%cplx) THEN
          CALL Write_dnCplxMat(basis_set%dnCBB)
        ELSE
          CALL Write_dnMat(basis_set%dnRBB)
        END IF

        !CALL RecWrite_basis(basis_set)
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------
      END SUBROUTINE sub_dnGB_TO_dnBB

      SUBROUTINE sub_dnGB_TO_dnBG(basis_set)
      use mod_dnSVM
      IMPLICIT NONE

      TYPE (basis), intent(inout) :: basis_set
!---------------------------------------------------------------------
!---------------------------------------------------------------------

      integer           :: nq,ib


!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='sub_dnGB_TO_dnBG'
!-----------------------------------------------------------

      nq = get_nq_FROM_basis(basis_set)
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb,nq',basis_set%nb,nq
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

  IF (.NOT. basis_set%packed_done) RETURN

  IF (NewBasisEl .AND. basis_set%ndim == 0) THEN
    IF (basis_set%cplx) THEN

      CALL dealloc_dnCplxMat(basis_set%dnCBG)
      CALL alloc_dnCplxMat(basis_set%dnCBG,basis_set%nb,basis_set%nb,nb_var_deriv=basis_set%ndim,nderiv=0)
      CALL Cplx_mat_id(basis_set%dnCBG%d0,basis_set%nb,basis_set%nb)

      CALL dealloc_dnCplxMat(basis_set%dnCBGwrho)
      CALL alloc_dnCplxMat(basis_set%dnCBGwrho,basis_set%nb,basis_set%nb,nb_var_deriv=basis_set%ndim,nderiv=0)
      CALL Cplx_mat_id(basis_set%dnCBGwrho%d0,basis_set%nb,basis_set%nb)

    ELSE

      CALL dealloc_dnMat(basis_set%dnRBG)
      CALL alloc_dnMat(basis_set%dnRBG,basis_set%nb,basis_set%nb,nb_var_deriv=basis_set%ndim,nderiv=0)

      CALL mat_id(basis_set%dnRBG%d0,basis_set%nb,basis_set%nb)

      CALL dealloc_dnMat(basis_set%dnRBGwrho)
      CALL alloc_dnMat(basis_set%dnRBGwrho,basis_set%nb,basis_set%nb,nb_var_deriv=basis_set%ndim,nderiv=0)

      CALL mat_id(basis_set%dnRBGwrho%d0,basis_set%nb,basis_set%nb)

    END IF
  ELSE


      IF (basis_set%cplx) THEN

        CALL dealloc_dnCplxMat(basis_set%dnCBG)
        CALL alloc_dnCplxMat(basis_set%dnCBG,basis_set%nb,nq,nb_var_deriv=basis_set%ndim,nderiv=0)

        basis_set%dnCBG%d0     = transpose(basis_set%dnCGB%d0)

        CALL dealloc_dnCplxMat(basis_set%dnCBGwrho)
        CALL alloc_dnCplxMat(basis_set%dnCBGwrho,basis_set%nb,nq,nb_var_deriv=basis_set%ndim,nderiv=0)

        DO ib=1,get_nb_FROM_basis(basis_set)
          basis_set%dnCBGwrho%d0(ib,:) = basis_set%dnCGB%d0(:,ib) * cmplx(basis_set%wrho,kind=Rkind)
        END DO

      ELSE

        CALL dealloc_dnMat(basis_set%dnRBG)

        CALL alloc_dnMat(basis_set%dnRBG,basis_set%nb,nq,nb_var_deriv=basis_set%ndim,nderiv=0)

        basis_set%dnRBG%d0 = transpose(basis_set%dnRGB%d0)

        CALL dealloc_dnMat(basis_set%dnRBGwrho)

        CALL alloc_dnMat(basis_set%dnRBGwrho,basis_set%nb,nq,nb_var_deriv=basis_set%ndim,nderiv=0)

        DO ib=1,get_nb_FROM_basis(basis_set)
          basis_set%dnRBGwrho%d0(ib,:) = basis_set%dnRGB%d0(:,ib) * basis_set%wrho
        END DO

      END IF
  END IF

!-----------------------------------------------------------
      IF (debug) THEN
        IF (basis_set%cplx) THEN
          CALL Write_dnCplxMat(basis_set%dnCGB)
        ELSE
          CALL Write_dnMat(basis_set%dnRGB)
        END IF
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------
      END SUBROUTINE sub_dnGB_TO_dnBG

      RECURSIVE SUBROUTINE Set_dnGGRep(basis_set,With_GG)

      TYPE (basis), intent(inout)  :: basis_set
      logical, intent(in)          :: With_GG
      !---------------------------------------------------------------------

      integer           :: ib


!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='Set_dnGGRep'
!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'With_GG',With_GG
      END IF
!-----------------------------------------------------------

      basis_set%dnGGRep = With_GG

      DO ib=1,basis_set%nb_basis
       CALL Set_dnGGRep(basis_set%tab_Pbasis(ib)%Pbasis,With_GG)
      END DO

!-----------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!-----------------------------------------------------------
      END SUBROUTINE Set_dnGGRep
      SUBROUTINE sub_dnGB_TO_dnGG(basis_set)
      USE mod_dnSVM
      IMPLICIT NONE

      TYPE (basis)  :: basis_set
!---------------------------------------------------------------------
!---------------------------------------------------------------------

!---------------------------------------------------------------------
      integer           :: i,j,nq,nb
      real(kind=Rkind), allocatable :: d0bxd0bT(:,:)      ! (qb,nq)
      real(kind=Rkind), allocatable :: vecp(:,:),valp(:)      ! (nq,nq)

      real(kind=Rkind), allocatable :: d0bxd0bT_inv(:,:)  ! (qb,nq)
      real(kind=Rkind), allocatable :: d0b_pseudoInv(:,:) ! (nb,nq)
      real(kind=Rkind), allocatable :: Check_bGB(:,:)     ! (nb,nq)
      real(kind=Rkind) :: Max_err_Check_bGB


      !logical :: SVD = .TRUE.
      logical :: SVD = .FALSE.
      !logical :: PseudoInverse = .FALSE.
      logical :: PseudoInverse = .TRUE.
      real(kind=Rkind) :: epsi,wmin
      logical :: pseudo_inv_OK
      logical :: Check = .TRUE.
      !logical :: Check = .FALSE.

!---------------------------------------------------------------------

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub='sub_dnGB_TO_dnGG'
!-----------------------------------------------------------
      IF (basis_set%ndim == 0) RETURN
      IF (.NOT. basis_set%dnGGRep) RETURN
      nq = get_nq_FROM_basis(basis_set)
      nb = basis_set%nb


      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'nb,nq',nb,nq

        write(out_unitp,*) 'dnRGB'
        CALL write_dnSVM(basis_set%dnRGB)
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------

      IF (.NOT. basis_set%packed_done) RETURN

      IF (basis_set%cplx) THEN
        STOP ' STOP in sub_dnGB_TO_dnGG, cplx impossible'
      ELSE


        CALL alloc_dnMat(basis_set%dnRGG,nq,nq,nb_var_deriv=basis_set%ndim,nderiv=2)

        CALL alloc_NParray(d0bxd0bT,     (/ nq,nq /),'d0bxd0bT',     name_sub)

        CALL alloc_NParray(d0bxd0bT_inv, (/ nq,nq /),'d0bxd0bT_inv', name_sub)
        CALL alloc_NParray(d0b_pseudoInv,(/ nb,nq /),'d0b_pseudoInv',name_sub)

        CALL alloc_NParray(Check_bGB,    (/ nq,nb /),'Check_bGB',    name_sub)

        CALL sub_dnGB_TO_dnBG(basis_set)


        IF (PseudoInverse) THEN
          d0bxd0bT = matmul(basis_set%dnRGB%d0,basis_set%dnRBG%d0)


          IF (debug) THEN
            write(out_unitp,*) 'd0bxd0bT'
            CALL Write_VecMat(d0bxd0bT,out_unitp,5)
            CALL flush_perso(out_unitp)
          END IF

          epsi = ONETENTH**10
          IF (SVD) THEN
            CALL inv_m1_TO_m2(d0bxd0bT,d0bxd0bT_inv,nq,1,epsi) ! SVD
            ! The SVD procedure does not work all the time !! Why ????

          IF (debug) THEN
            write(out_unitp,*) ' Computation with pseudonverse, SVD'
          END IF

          ELSE
            CALL alloc_NParray(vecP,(/ nq,nq /),'vecP',name_sub)
            CALL alloc_NParray(valP,(/ nq /),'valP',name_sub)

            CALL diagonalization(d0bxd0bT,valP,VecP,nq,3,-1,.FALSE.)

            IF (debug) THEN
              write(out_unitp,*) 'diago : ValP(:)',ValP(:)
              CALL flush_perso(out_unitp)
            END IF
            d0bxd0bT(:,:) = ZERO
            DO i=1,nb
              d0bxd0bT(:,i) = VecP(:,i)/valP(i)
            END DO

!           pseudo_inv_OK = .FALSE.
!           DO
!             IF (pseudo_inv_OK) EXIT
!             IF (debug) write(out_unitp,*) 'diago : count non zero',count(valP > epsi)
!             IF (count(valP > epsi) > nb) THEN
!               epsi = epsi * TWO
!             ELSE
!               pseudo_inv_OK = .TRUE.
!             END IF
!           END DO

!           d0bxd0bT(:,:) = ZERO
!           DO i=1,nq
!             IF (abs(valP(i)) > epsi) d0bxd0bT(:,i) = VecP(:,i)/valP(i)
!           END DO
            d0bxd0bT_inv = matmul(d0bxd0bT,transpose(VecP))
            !CALL Write_VecMat(matmul(d0bxd0bT,d0bxd0bT_inv),out_unitp,5)
            !CALL Write_VecMat(matmul(d0bxd0bT_inv,d0bxd0bT),out_unitp,5)

            CALL dealloc_NParray(valP,'valP',name_sub)
            CALL dealloc_NParray(vecP,'vecP',name_sub)

            IF (debug) THEN
              write(out_unitp,*) ' Computation with pseudonverse, diago'
            END IF

          END IF

          IF (debug) THEN
            write(out_unitp,*) 'd0bxd0bT_inv'
            CALL Write_VecMat(d0bxd0bT_inv,out_unitp,5)
            CALL flush_perso(out_unitp)
          END IF



          d0b_pseudoInv =  matmul(basis_set%dnRBG%d0,d0bxd0bT_inv)
          IF (debug) THEN
            write(out_unitp,*) 'True d0b_pseudoInv'
            CALL Write_VecMat(d0b_pseudoInv,out_unitp,5)
            CALL flush_perso(out_unitp)
          END IF
        ELSE

          IF (debug) THEN
            write(out_unitp,*) ' Standard computation'
          END IF

          DO i=1,nb
            d0b_pseudoInv(i,:) =  basis_set%dnRBG%d0(i,:) * basis_set%wrho(:)
          END DO

          IF (debug) THEN
            write(out_unitp,*) 'Standard d0b_pseudoInv'
            CALL Write_VecMat(d0b_pseudoInv,out_unitp,5)
            CALL flush_perso(out_unitp)
          END IF
        END IF

        basis_set%dnRGG%d0 =  matmul(basis_set%dnRGB%d0(:,1:nb),d0b_pseudoInv)

        IF (debug) THEN
          !write(out_unitp,*) 'dnRGG%d0'
          !CALL Write_VecMat(basis_set%dnRGG%d0,out_unitp,5)
          write(out_unitp,*) 'dnRGG%d0 done'
          CALL flush_perso(out_unitp)
        END IF

        IF (.NOT. associated(basis_set%dnRGB%d1) .OR.                   &
                              .NOT. associated(basis_set%dnRGB%d2)) THEN
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) '  dnRGB%d1 or dnRGB%d2',                  &
            associated(basis_set%dnRGB%d1),associated(basis_set%dnRGB%d2)
          write(out_unitp,*) '  dnRGB%d1 or dnRGB%d2 are not associated!!'
          STOP
        END IF

        DO i=1,basis_set%ndim
          basis_set%dnRGG%d1(:,:,i) =  matmul(basis_set%dnRGB%d1(:,1:nb,i),d0b_pseudoInv)
        END DO
        IF (debug) THEN
          write(out_unitp,*) 'dnRGG%d1'
          DO i=1,basis_set%ndim
          !  CALL Write_VecMat(basis_set%dnRGG%d1(:,:,i),out_unitp,5)
          END DO
          write(out_unitp,*) 'dnRGG%d1 done'
          CALL flush_perso(out_unitp)
        END IF

        DO i=1,basis_set%ndim
        DO j=1,basis_set%ndim
          basis_set%dnRGG%d2(:,:,i,j) =  matmul(basis_set%dnRGB%d2(:,1:nb,i,j),d0b_pseudoInv)
        END DO
        END DO
        IF (debug) THEN
          write(out_unitp,*) 'dnRGG%d2 done'
          CALL flush_perso(out_unitp)
        END IF

        IF (Check) THEN
          IF (debug) THEN
            write(out_unitp,*) '-------------------------------------------'
            write(out_unitp,*) 'Check dnRGG, nb,nq',nb,nq
            CALL flush_perso(out_unitp)
          END IF
          Check_bGB = basis_set%dnRGB%d0(:,1:nb)-matmul(basis_set%dnRGG%d0,basis_set%dnRGB%d0(:,1:nb))
          Max_err_Check_bGB = maxval(abs(Check_bGB))
          IF (debug .OR. maxval(abs(Check_bGB)) > ONETENTH**8) THEN
            write(out_unitp,'(a,4x,e9.2)') 'Check_bGB%d0',maxval(abs(Check_bGB))
            !CALL Write_VecMat(Check_bGB,out_unitp,5)
          END IF
          CALL flush_perso(out_unitp)

          DO i=1,basis_set%ndim
            Check_bGB = basis_set%dnRGB%d1(:,1:nb,i)-matmul(basis_set%dnRGG%d1(:,:,i),basis_set%dnRGB%d0(:,1:nb))
            Max_err_Check_bGB = max(Max_err_Check_bGB,maxval(abs(Check_bGB)))
            IF (debug .OR. maxval(abs(Check_bGB)) > ONETENTH**8) THEN
              write(out_unitp,'(a,x,i0,2x,e9.2)') 'Check_bGB%d1',i,maxval(abs(Check_bGB))
              !CALL Write_VecMat(Check_bGB,out_unitp,5)
            END IF
          END DO
          CALL flush_perso(out_unitp)

          DO i=1,basis_set%ndim
          DO j=1,basis_set%ndim
            Check_bGB = basis_set%dnRGB%d2(:,1:nb,i,j)-matmul(basis_set%dnRGG%d2(:,:,i,j),basis_set%dnRGB%d0(:,1:nb))
            Max_err_Check_bGB = max(Max_err_Check_bGB,maxval(abs(Check_bGB)))
            IF (debug .OR. maxval(abs(Check_bGB)) > ONETENTH**8) THEN
              write(out_unitp,'(a,x,i0x,i0,e9.2)') 'Check_bGB%d2',i,j,maxval(abs(Check_bGB))
              !CALL Write_VecMat(Check_bGB,out_unitp,5)
            END IF
          END DO
          END DO
          IF (debug) write(out_unitp,*) '-------------------------------------------'
          CALL flush_perso(out_unitp)

          IF (Max_err_Check_bGB > ONETENTH**6) THEN
            write(out_unitp,*) 'ERROR in ',name_sub
            write(out_unitp,'(a,e9.2)') 'Max_err_Check_bGB is too large! ',Max_err_Check_bGB
            STOP
          END IF
        END IF

        CALL dealloc_NParray(d0bxd0bT,'d0bxd0bT',name_sub)
        CALL dealloc_NParray(d0bxd0bT_inv,'d0bxd0bT_inv',name_sub)

        CALL dealloc_NParray(d0b_pseudoInv,'d0b_pseudoInv',name_sub)
        CALL dealloc_NParray(Check_bGB,'Check_bGB',name_sub)

      END IF

      basis_set%dnGGRep_done = .TRUE.

!-----------------------------------------------------------
      IF (debug) THEN
        CALL alloc_dnMat(basis_set%dnRGG)
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
!-----------------------------------------------------------
      END SUBROUTINE sub_dnGB_TO_dnGG

      SUBROUTINE pack_basis(basis_set,sortX)

      TYPE (basis), intent(inout) :: basis_set
      logical, optional :: sortX

      integer              :: ib,iq,iq1,iq0,i,j,nq,i_SG,ndimi,iib,idim1,idim2
      integer              :: nqo
      integer              :: nDvalB(basis_set%nDindB%ndim)
      logical              :: sortX_loc
      integer, allocatable :: tab_iqXmin(:)
      real (kind=Rkind)    :: Xmin(basis_set%ndim)
      real (kind=Rkind)    :: X(basis_set%ndim)
      TYPE (basis)         :: basis_loc
      real (kind=Rkind), allocatable    :: RvecB(:),RvecG(:)


      logical :: inferior ! function

!----- for debuging --------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='pack_basis'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!-----------------------------------------------------------
      IF (basis_set%ndim == 0) RETURN
      IF (.NOT. basis_set%packed) RETURN
      IF (basis_set%packed_done) RETURN

      nqo = get_nq_FROM_basis(basis_set)
      !write(out_unitp,*) 'BEGINNING ',name_sub
      !write(out_unitp,*) 'nb,nq',basis_set%nb,nqo

      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) '--------------------------'
        CALL RecWrite_basis(basis_set,write_all=.TRUE.)
        write(out_unitp,*) '--------------------------'

        write(out_unitp,*) 'unpacked Grid: ',nqo
        DO iq=1,nqo
          CALL Rec_x(X,basis_set,iq)
          write(out_unitp,*) 'iq,x',iq,X(:)
        END DO
        !STOP
      END IF
!-----------------------------------------------------------------------

   sortX_loc = .FALSE.
   IF (present(sortX)) sortX_loc = sortX
   IF (basis_set%SparseGrid_type == 2) sortX_loc = .FALSE.

   IF (basis_set%nb_basis > 0 .AND. basis_set%SparseGrid_type == 2) THEN
     IF (sortX_loc) THEN
STOP 'pack and SG2 does not work!!!'
       ! a temp basis with all grid points
       CALL Set_nq_OF_basis(basis_loc,nqo)
       basis_loc%ndim = basis_set%ndim
       basis_loc%nb   = basis_set%nb

       ! first the basis without sortX
       CALL alloc_dnb_OF_basis(basis_loc)
       CALL alloc_NParray(RvecB,(/basis_loc%nb/),'RvecB',name_sub)

       DO ib=1,basis_loc%nb
         RvecB = ZERO
         RvecB(ib) = ONE
         CALL RecRvecB_TO_RVecG(RvecB,basis_loc%dnRGB%d0(:,ib),basis_set%nb,nqo,basis_set)
         DO i=1,basis_loc%ndim
           basis_loc%dnRGB%d1(:,ib,i) = basis_loc%dnRGB%d0(:,ib)
           CALL DerivOp_TO_RVecG(basis_loc%dnRGB%d1(:,ib,i),nqo,basis_set,(/basis_set%iQdyn(i),0/))
         END DO

         DO i=1,basis_loc%ndim
         DO j=1,basis_loc%ndim
           basis_loc%dnRGB%d2(:,ib,i,j) = basis_loc%dnRGB%d0(:,ib)
           CALL DerivOp_TO_RVecG(basis_loc%dnRGB%d2(:,ib,i,j),nqo,basis_set, &
             (/basis_set%iQdyn(i),basis_set%iQdyn(j)/))
         END DO
         END DO

       END DO

       CALL dealloc_NParray(RvecB,'RvecB',name_sub)


       !then with the new points
       CALL alloc_NParray(tab_iqXmin,(/ nqo /),'tab_iqXmin',name_sub)
       tab_iqXmin(:) = (/ (i,i=1,nqo) /)

       ! first the new number of points
       DO iq=1,nqo

         CALL Rec_x(Xmin(:),basis_set,tab_iqXmin(iq))
         !write(out_unitp,*) 'tab_iqXmin(iq),Xmin(:)',iq,tab_iqXmin(iq),Xmin(:)

         ! 3d: find iqm
         DO iq1=iq+1,nqo
           CALL Rec_x(X,basis_set,tab_iqXmin(iq1))

           IF (inferior_tab(X,Xmin)) THEN
             iq0             = tab_iqXmin(iq)
             tab_iqXmin(iq)  = tab_iqXmin(iq1)
             tab_iqXmin(iq1) = iq0
             Xmin(:) = X(:)
           END IF
         END DO
         !write(out_unitp,*) 'tab_iqXmin(iq),Xmin(:)',iq,tab_iqXmin(iq),Xmin(:)
       END DO

       ! first the new number of points
       CALL Rec_x(Xmin(:),basis_set,tab_iqXmin(1))
       nq = 1
       DO iq=2,nqo

         CALL Rec_x(X(:),basis_set,tab_iqXmin(iq))
         IF (compare_tab(X,Xmin)) THEN
           tab_iqXmin(iq) = -tab_iqXmin(iq)
         ELSE
           Xmin(:) = X(:)
           nq = nq + 1
         END IF

       END DO
       !IF (print_level > -1) write(out_unitp,*) 'old/new nq',nqo,nq
       write(out_unitp,*) 'old/new nq',nqo,nq
!nq=nqo
!tab_iqXmin = abs(tab_iqXmin)
       ! set up the basis set with the right number of points
       CALL Set_nq_OF_basis(basis_set,nq)

       CALL alloc_xw_OF_basis(basis_set)
       CALL alloc_dnb_OF_basis(basis_set)

       iq0 = 0
       DO iq=1,nqo ! old nq

         IF (tab_iqXmin(iq) > 0) THEN
           iq1 = tab_iqXmin(iq)
           iq0 = iq0 + 1

           CALL Rec_x(basis_set%x(:,iq0),basis_set,iq1)
           basis_set%rho(iq0)  = Rec_rhonD(basis_set,iq1)
           basis_set%wrho(iq0) = Rec_WrhonD(basis_set,iq1)
           basis_set%w(iq0)    = Rec_WnD(basis_set,iq1)

           basis_set%dnRGB%d0(iq0,:)     = basis_loc%dnRGB%d0(iq1,:)
           basis_set%dnRGB%d1(iq0,:,:)   = basis_loc%dnRGB%d1(iq1,:,:)
           basis_set%dnRGB%d2(iq0,:,:,:) = basis_loc%dnRGB%d2(iq1,:,:,:)

         ELSE  ! two (or more) identical points, just one is consider. But the weight must be changed
           iq1 = abs(tab_iqXmin(iq))

           write(6,*) 'd0RGB,wrho',iq0,basis_set%dnRGB%d0(iq0,:),basis_set%wrho(iq0)
           write(6,*) 'd0RGB,wrho',iq1,basis_loc%dnRGB%d0(iq1,:),Rec_WrhonD(basis_set,iq1)


           ! be carrefull, rho remains unchanged. It just function of Xmin
           ! the sum has to be done only on w(:) and wrho(:)
           basis_set%wrho(iq0) = basis_set%wrho(iq0) + Rec_WrhonD(basis_set,iq1)
           basis_set%w(iq0)    = basis_set%w(iq0)    + Rec_WnD(basis_set,iq1)

         END IF

       END DO
       IF (iq0 /= nq) STOP 'Wrong nq'

       CALL dealloc_NParray(tab_iqXmin,'tab_iqXmin',name_sub)
       CALL dealloc_basis(basis_loc)

     ELSE

       CALL Set_nq_OF_basis(basis_set,nqo)

       ! first the basis without sortX
       CALL alloc_xw_OF_basis(basis_set)
       CALL alloc_dnb_OF_basis(basis_set)
       CALL alloc_NParray(RvecB,(/basis_set%nb/),'RvecB',name_sub)

       DO ib=1,basis_set%nb
         RvecB = ZERO
         RvecB(ib) = ONE
         CALL RecRvecB_TO_RVecG(RvecB,basis_set%dnRGB%d0(:,ib),basis_set%nb,nqo,basis_set)
         DO i=1,basis_set%ndim
           basis_set%dnRGB%d1(:,ib,i) = basis_set%dnRGB%d0(:,ib)
           CALL DerivOp_TO_RVecG(basis_set%dnRGB%d1(:,ib,i),nqo,basis_set,(/basis_set%iQdyn(i),0/))
         END DO

         DO i=1,basis_set%ndim
         DO j=1,basis_set%ndim
           basis_set%dnRGB%d2(:,ib,i,j) = basis_set%dnRGB%d0(:,ib)
           CALL DerivOp_TO_RVecG(basis_set%dnRGB%d2(:,ib,i,j),nqo,basis_set, &
             (/basis_set%iQdyn(i),basis_set%iQdyn(j)/))
         END DO
         END DO

       END DO

       CALL dealloc_NParray(RvecB,'RvecB',name_sub)

       DO iq=1,nqo
         CALL Rec_x(basis_set%x(:,iq),basis_set,iq)
         basis_set%rho(iq)  = Rec_rhonD(basis_set,iq)
         basis_set%wrho(iq) = Rec_WrhonD(basis_set,iq)
         basis_set%w(iq)    = Rec_WnD(basis_set,iq)
       END DO
     END IF

   ELSE

     IF (sortX_loc) THEN

       CALL alloc_NParray(tab_iqXmin,(/ nqo /),'tab_iqXmin',name_sub)
       tab_iqXmin(:) = (/ (i,i=1,nqo) /)

       ! first the new number of points
       DO iq=1,nqo

         CALL Rec_x(Xmin(:),basis_set,tab_iqXmin(iq))
         !write(out_unitp,*) 'tab_iqXmin(iq),Xmin(:)',iq,tab_iqXmin(iq),Xmin(:)

         ! 3d: find iqm
         DO iq1=iq+1,nqo
           CALL Rec_x(X,basis_set,tab_iqXmin(iq1))

           IF (inferior_tab(X,Xmin)) THEN
             iq0             = tab_iqXmin(iq)
             tab_iqXmin(iq)  = tab_iqXmin(iq1)
             tab_iqXmin(iq1) = iq0
             Xmin(:) = X(:)
           END IF
         END DO
         !write(out_unitp,*) 'tab_iqXmin(iq),Xmin(:)',iq,tab_iqXmin(iq),Xmin(:)
       END DO

       ! first the new number of points
       CALL Rec_x(Xmin(:),basis_set,tab_iqXmin(1))
       nq = 1
       DO iq=2,nqo

         CALL Rec_x(X(:),basis_set,tab_iqXmin(iq))
         IF (compare_tab(X,Xmin)) THEN
           tab_iqXmin(iq) = -tab_iqXmin(iq)
         ELSE
           Xmin(:) = X(:)
           nq = nq + 1
         END IF

       END DO
       IF (print_level > -1) write(out_unitp,*) 'old/new nq',nqo,nq

       ! set up the basis set with the right number of points
       CALL Set_nq_OF_basis(basis_set,nq)

       CALL alloc_xw_OF_basis(basis_set)
       CALL alloc_dnb_OF_basis(basis_set)

       iq0 = 0
       DO iq=1,nqo ! old nq

         IF (tab_iqXmin(iq) > 0) THEN
           iq1 = tab_iqXmin(iq)
           iq0 = iq0 + 1

           CALL Rec_x(basis_set%x(:,iq0),basis_set,iq1)
           basis_set%rho(iq0)  = Rec_rhonD(basis_set,iq1)
           basis_set%wrho(iq0) = Rec_WrhonD(basis_set,iq1)
           basis_set%w(iq0)    = Rec_WnD(basis_set,iq1)

           DO ib=1,basis_set%nb

             CALL Rec_d0d1d2bnD(basis_set%dnRGB%d0(iq0,ib),               &
                                basis_set%dnRGB%d1(iq0,ib,:),             &
                                basis_set%dnRGB%d2(iq0,ib,:,:),           &
                                basis_set,iq1,ib)
           END DO

         ELSE  ! two (or more) identical points, just one is consider. But the weight must be changed
           iq1 = abs(tab_iqXmin(iq))

           ! be carrefull, rho remains unchanged. It just function of Xmin
           ! the sum has to be done only on w(:) and wrho(:)
           basis_set%wrho(iq0) = basis_set%wrho(iq0) + Rec_WrhonD(basis_set,iq1)
           basis_set%w(iq0)    = basis_set%w(iq0)    + Rec_WnD(basis_set,iq1)

         END IF

       END DO
       IF (iq0 /= nq) STOP 'Wrong nq'

       CALL dealloc_NParray(tab_iqXmin,'tab_iqXmin',name_sub)

     ELSE

       CALL alloc_xw_OF_basis(basis_set)
       CALL alloc_dnb_OF_basis(basis_set)

       DO iq=1,nqo

         CALL Rec_x(basis_set%x(:,iq),basis_set,iq)

         basis_set%rho(iq)  = Rec_rhonD(basis_set,iq)
         basis_set%wrho(iq) = Rec_WrhonD(basis_set,iq)
         basis_set%w(iq)    = Rec_WnD(basis_set,iq)

         DO ib=1,basis_set%nb
           CALL Rec_d0d1d2bnD(basis_set%dnRGB%d0(iq,ib),                 &
                              basis_set%dnRGB%d1(iq,ib,:),               &
                              basis_set%dnRGB%d2(iq,ib,:,:),             &
                              basis_set,iq,ib)
         END DO
       END DO

     END IF
   END IF

   DO ib=1,basis_set%nb
     CALL Rec_ndim_index(basis_set,basis_set%tab_ndim_index(:,ib),ib)
     IF (debug) write(out_unitp,*) 'ib,ndim_index',ib,basis_set%tab_ndim_index(:,ib)
   END DO

   IF (associated(basis_set%tab_PbasisSG)) THEN
     CALL dealloc_array(basis_set%tab_PbasisSG,                    &
                       'basis_set%tab_PbasisSG',name_sub)
     basis_set%nb_SG = 0
     !basis_set%SparseGrid = .FALSE.
   END IF

   basis_set%packed_done = .TRUE.

   CALL check_ortho_basis(basis_set)

  !-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'pack done '
    CALL RecWrite_basis(basis_set,write_all=.TRUE.)

    nqo = get_nq_FROM_basis(basis_set)
    write(out_unitp,*) 'packed Grid: ',nqo
    DO iq=1,nqo
      write(out_unitp,*) 'iq,x',iq,basis_set%x(:,iq)
    END DO

    write(out_unitp,*) 'END ',name_sub


  END IF
  !-----------------------------------------------------------

END SUBROUTINE pack_basis
!=============================================================
!
!      check the overlap matrix of a basis
!
!=============================================================
      !!@description: check the overlap matrix of a basis
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      SUBROUTINE check_ortho_basis(basis_temp,test_stop)

      TYPE (basis)  :: basis_temp
      logical, optional :: test_stop

!---------- working variables ----------------------------------------
      real (kind=Rkind), allocatable :: tbasiswrho(:,:)
      real (kind=Rkind), allocatable :: matS(:,:)

      real (kind=Rkind) ::  Sij,Sii
      integer           ::  i,j,nq

      real (kind=Rkind) :: max_Sii,max_Sij
      logical           :: loc_test_stop,Print_basis

!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='check_ortho_basis'
      logical,parameter :: debug= .FALSE.
      !logical,parameter :: debug= .TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'check_basis: ',basis_temp%check_basis
        write(out_unitp,*) 'packed_done: ',basis_temp%packed_done
        CALL RecWrite_basis(basis_temp,write_all=.TRUE.)
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
      IF (basis_temp%ndim == 0) RETURN
      IF (.NOT. basis_temp%check_basis .OR. .NOT. basis_temp%packed_done) RETURN
      !IF (.NOT. basis_temp%packed_done) RETURN
      nq = get_nq_FROM_basis(basis_temp)
      IF ( (basis_temp%nb*nq) > 1000000000) RETURN

      Print_basis = basis_temp%print_info_OF_basisDP .AND. print_level > -1 .OR. debug


      IF (present(test_stop)) THEN
        loc_test_stop = test_stop
      ELSE
        loc_test_stop = .TRUE.
      END IF
      max_Sii = ZERO
      max_Sij = ZERO
      IF (basis_temp%cplx) THEN
        DO i=1,basis_temp%nb
          Sii = dot_product(basis_temp%dnCGB%d0(:,i)*basis_temp%wrho(:),    &
                            basis_temp%dnCGB%d0(:,i) ) -CONE
          IF (abs(Sii) .GT. max_Sii) max_Sii = abs(Sii)
          DO j=i+1,basis_temp%nb
            Sij = dot_product(basis_temp%dnCGB%d0(:,i)*basis_temp%wrho(:),  &
                              basis_temp%dnCGB%d0(:,j))
            IF (abs(Sij) .GT. max_Sij) max_Sij = abs(Sij)
          END DO
        END DO
      ELSE
        CALL alloc_NParray(tbasiswrho,(/ basis_temp%nb,nq /),             &
                          'tbasiswrho',name_sub)

        CALL alloc_NParray(matS,(/ basis_temp%nb,basis_temp%nb /),        &
                          'matS',name_sub)

        DO i=1,basis_temp%nb
          tbasiswrho(i,:) = basis_temp%dnRGB%d0(:,i) * basis_temp%wrho(:)
        END DO
        matS(:,:) = matmul(tbasiswrho,basis_temp%dnRGB%d0)

        CALL sub_ana_S(matS,basis_temp%nb,max_Sii,max_Sij,Print_basis)


        IF (max_Sii > ONETENTH**5 .OR. max_Sij > ONETENTH**5) THEN
          write(out_unitp,*) ' Overlap Matrix:'
          CALL Write_nDindex(basis_temp%nDindB)
          CALL Write_Mat(matS,out_unitp,5)
        END IF


        CALL dealloc_NParray(tbasiswrho,'tbasiswrho',name_sub)
        CALL dealloc_NParray(matS,'matS',name_sub)

      END IF


      IF (max_Sii > ONETENTH**5 .OR. max_Sij > ONETENTH**5) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' the basis is not orthonormal !!'
        IF (loc_test_stop) THEN
          CALL RecWrite_basis(basis_temp)
          write(out_unitp,*) ' The basis is not orthonormal !!'
          STOP 'ERROR: None orthonormal Basis'
        END IF
      END IF


!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


      END SUBROUTINE check_ortho_basis

      integer FUNCTION Get_nb_FROM_l_OF_PrimBasis(l,PrimBasis)
      IMPLICIT NONE
      TYPE (Basis), intent(inout)  :: PrimBasis
      integer, intent(in)          :: l
      integer :: u

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Get_nb_FROM_l_OF_PrimBasis'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) '  With_L',PrimBasis%With_L
         write(out_unitp,*) '  Norm_OF_nDindB',PrimBasis%Norm_OF_nDindB
         write(out_unitp,*) '  weight_OF_nDindB',PrimBasis%weight_OF_nDindB
       END IF
!-----------------------------------------------------------
       IF (PrimBasis%With_L) THEN
         IF (PrimBasis%contrac) THEN
           PrimBasis%nbc = get_n_FROM_Basis_L_TO_n(PrimBasis%L_TO_nb,L)
           Get_nb_FROM_l_OF_PrimBasis = PrimBasis%nb
         ELSE
           Get_nb_FROM_l_OF_PrimBasis = get_n_FROM_Basis_L_TO_n(PrimBasis%L_TO_nb,L)
         END IF
       ELSE
          IF (PrimBasis%nb < 1) THEN
            Get_nb_FROM_l_OF_PrimBasis = 1 +                            &
                int(PrimBasis%Norm_OF_nDindB/PrimBasis%weight_OF_nDindB)
          ELSE
            Get_nb_FROM_l_OF_PrimBasis = PrimBasis%nb
          END IF
       END IF

       IF (Get_nb_FROM_l_OF_PrimBasis < 1) THEN
         write(out_unitp,*) 'ERROR in ',name_sub
         write(out_unitp,*) '  L',L
         write(out_unitp,*) '  contrac',PrimBasis%contrac
         write(out_unitp,*) '  nb',PrimBasis%nb
         write(out_unitp,*) '  With_L',PrimBasis%With_L
         write(out_unitp,*) '  Norm_OF_nDindB',PrimBasis%Norm_OF_nDindB
         write(out_unitp,*) '  weight_OF_nDindB',PrimBasis%weight_OF_nDindB
         STOP
       END IF

!     -------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'nb',Get_nb_FROM_l_OF_PrimBasis
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -------------------------------------------------------

      END FUNCTION Get_nb_FROM_l_OF_PrimBasis

      integer FUNCTION Get_nq_FROM_l_OF_PrimBasis(L,PrimBasis)
      IMPLICIT NONE
      TYPE (Basis), intent(in)  :: PrimBasis
      integer, intent(in)       :: L

      integer :: u

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='Get_nq_FROM_l_OF_PrimBasis'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'L',L
         write(out_unitp,*) 'PrimBasis%L_SparseBasis',PrimBasis%L_SparseBasis
       END IF
!-----------------------------------------------------------

       IF (PrimBasis%L_SparseBasis < 0) THEN
         Get_nq_FROM_l_OF_PrimBasis =                                   &
                            get_n_FROM_Basis_L_TO_n(PrimBasis%L_TO_nq,L)
       ELSE
         Get_nq_FROM_l_OF_PrimBasis =                                   &
                get_n_FROM_Basis_L_TO_n(PrimBasis%L_TO_nq,L,PrimBasis%L_SparseBasis)
       END IF
!     -------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'nq',Get_nq_FROM_l_OF_PrimBasis
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -------------------------------------------------------

      END FUNCTION Get_nq_FROM_l_OF_PrimBasis

      RECURSIVE SUBROUTINE Rec_ndim_index(BasisnD,ndim_index,ib)
      USE mod_system
      implicit none

!----- variables for the Basis and quadrature points -----------------
      TYPE (Basis)         :: BasisnD
      integer              :: ndim_index(:)
      integer, intent(in)  :: ib

!------ working variables ---------------------------------
      integer       :: i,idim1,idim2,iib,i_SG,L
      integer       :: nDvalB(BasisnD%nDindB%ndim)


!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
       character (len=*), parameter :: name_sub = 'Rec_ndim_index'

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
       END IF
!-----------------------------------------------------------
       IF (BasisnD%packed_done) THEN
         ndim_index(:) = BasisnD%tab_ndim_index(:,ib)
       ELSE ! BasisnD%nb_basis MUST BE > 0
         IF (BasisnD%nb_basis == 0 ) THEN
            CALL RecWriteMini_basis(BasisnD)
            write(6,*) ' ERROR in ',name_sub
            write(6,*) ' This basis should be packed'

            STOP
         END IF

         SELECT CASE (BasisnD%SparseGrid_type)
         CASE (0) ! Direct product
           CALL calc_nDindex(BasisnD%nDindB,ib,nDvalB)
           idim1 = 1
           idim2 = 0
           DO i=1,BasisnD%nb_basis
             iib = nDvalB(i)
             idim2 = idim2 + BasisnD%tab_Pbasis(i)%Pbasis%ndim
             CALL Rec_ndim_index(BasisnD%tab_Pbasis(i)%Pbasis,ndim_index(idim1:idim2),iib)
             idim1 = idim1 + BasisnD%tab_Pbasis(i)%Pbasis%ndim
           END DO

         CASE (1) ! Sparse basis (Smolyak 1st implementation)
           CALL Rec_ndim_index(BasisnD%tab_PbasisSG(1)%Pbasis,ndim_index,ib)

         CASE (2,4) ! Sparse basis (Smolyak 2d  implementation)

           CALL calc_nDindex(BasisnD%nDindB,ib,nDvalB)
           idim1 = 1
           idim2 = 0
           DO i=1,BasisnD%nb_basis
             iib = nDvalB(i)
             L   = BasisnD%L_SparseGrid
             idim2 = idim2 + BasisnD%tab_basisPrimSG(L,i)%ndim
             CALL Rec_ndim_index(BasisnD%tab_basisPrimSG(L,i),ndim_index(idim1:idim2),iib)
             idim1 = idim1 + BasisnD%tab_basisPrimSG(L,i)%ndim
           END DO

         CASE DEFAULT
           write(out_unitp,*) ' ERROR in',name_sub
           write(out_unitp,*) ' WRONG SparseGrid_type',BasisnD%SparseGrid_type
           write(out_unitp,*) ' The possibilities are: 0, 1, 2, 4'
           STOP
         END SELECT

       END IF

      !-------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'ib,ndim_index',ib,ndim_index
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
      !-------------------------------------------------------
      END SUBROUTINE Rec_ndim_index
      SUBROUTINE calc_InD_FROM_ndim_index(BasisnD,ndim_index,InD,ibGuess)
      USE mod_system
      implicit none

!----- variables for the Basis and quadrature points -----------------
      TYPE (Basis)           :: BasisnD
      integer,           intent(in)    :: ndim_index(:)
      integer,           intent(inout) :: InD
      integer, optional, intent(in)    :: ibGuess

!------ working variables ---------------------------------
      integer              :: ib,ibGuess_loc
      integer              :: ndim_index_loc(BasisnD%ndim)


!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING calc_InD_FROM_ndim_index'
       END IF
!-----------------------------------------------------------
       IF (present(ibGuess)) THEN
         ibGuess_loc = ibGuess
       ELSE
         ibGuess_loc = 1
       END IF

       DO ib=ibGuess_loc,BasisnD%nb
         CALL Rec_ndim_index(BasisnD,ndim_index_loc,ib)
         IF (compare_tab(ndim_index_loc,ndim_index)) EXIT
       END DO
       InD = ib

       IF (present(ibGuess) .AND. InD > BasisnD%nb) THEN
         DO ib=1,ibGuess_loc-1
           CALL Rec_ndim_index(BasisnD,ndim_index_loc,ib)
           IF (compare_tab(ndim_index_loc,ndim_index)) EXIT
         END DO
         IF (ib < ibGuess_loc) InD = ib
       END IF

       !write(out_unitp,*) 'ndim_index',ndim_index,':',ibGuess_loc,ib,BasisnD%nb


!     -------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'InD,ndim_index',InD,ndim_index
        write(out_unitp,*)
        write(out_unitp,*) 'END calc_InD_FROM_ndim_index'
      END IF
!     -------------------------------------------------------

      END SUBROUTINE calc_InD_FROM_ndim_index
!=====================================================================
!
!  calculation of WrhonD, rhonD, WnD
!                 x
!
!=====================================================================
      !!@description: calculation of WrhonD, rhonD, Wn
      !!@param: TODO
      !!@param: TODO
      RECURSIVE FUNCTION Rec_WrhonD(BasisnD,iq,OldPara) result(WrhonD)
      USE mod_system
      USE mod_basis_BtoG_GtoB_SGType4, ONLY : getbis_tab_nq
      implicit none

!----- variables for the Basis and quadrature points -----------------
      real (kind=Rkind)   :: WrhonD
      TYPE (Basis)        :: BasisnD
      integer, intent(in) :: iq
      TYPE (OldParam), intent(inout), optional :: OldPara

!------ working variables ---------------------------------
      integer           :: i,iq_SG,nq,L,ib,iq_ib
      integer           :: nDval(BasisnD%nDindG%ndim)
      integer           :: nDval_SG2(BasisnD%nb_basis)
      integer           :: nDl_SG2(BasisnD%nb_basis)

      integer             :: tab_nq(BasisnD%nb_basis)
      integer             :: tab_l(BasisnD%nb_basis)
      integer             :: tab_iq(BasisnD%nb_basis)

      integer           :: i_SG
      integer           :: err_sub


!----- for debuging --------------------------------------------------
       logical, parameter :: debug=.FALSE.
       !logical, parameter :: debug=.TRUE.
       character (len=*), parameter :: name_sub = 'Rec_WrhonD'
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'iq',iq
         write(out_unitp,*) 'BasisnD%nb_basis',BasisnD%nb_basis
         write(out_unitp,*) 'BasisnD%packed_done',BasisnD%packed_done
       END IF
!-----------------------------------------------------------
!
       IF (BasisnD%ndim == 0) THEN
         WrhonD = ONE ! no grid point
       ELSE IF (BasisnD%packed_done) THEN
         !CALL RecWrite_basis(BasisnD,write_all=.TRUE.)
         !IF (.NOT. allocated(BasisnD%wrho)) STOP 'problem with primitive basis. wrho is not allocated !!'
         !write(6,*) 'shape wrho, iq',shape(BasisnD%wrho),iq
         WrhonD = BasisnD%wrho(iq)
       ELSE ! BasisnD%nb_basis MUST BE > 0
         IF (BasisnD%nb_basis == 0 ) STOP ' ERROR with packed in Rec_WrhonD !!!'

         SELECT CASE (BasisnD%SparseGrid_type)
         CASE (0) ! Direct product
           CALL calc_nDindex(BasisnD%nDindG,iq,nDval)
           IF (debug) write(out_unitp,*) 'DP: nDval',nDval(:)
           WrhonD = ONE
           DO i=1,BasisnD%nb_basis
             WrhonD = WrhonD * Rec_WrhonD(BasisnD%tab_Pbasis(i)%Pbasis,nDval(i))
           END DO

         CASE (1) ! Sparse basis (Smolyak 1st implementation)
           iq_SG = iq
           DO i_SG=1,BasisnD%nb_SG
             nq = get_nq_FROM_basis(BasisnD%tab_PbasisSG(i_SG)%Pbasis)
             IF (iq_SG <= nq) EXIT
             iq_SG = iq_SG - nq
           END DO
           IF (debug) write(out_unitp,*) 'i_SG,iq_SG',i_SG,iq_SG
           WrhonD = BasisnD%WeightSG(i_SG) *                            &
                    Rec_WrhonD(BasisnD%tab_PbasisSG(i_SG)%Pbasis,iq_SG)

         CASE (2) ! Sparse basis (Smolyak 2d  implementation)
           CALL get_Tabiq_Tabil_FROM_iq_old(nDval_SG2,nDl_SG2,          &
                                     i_SG,iq_SG,iq,BasisnD%para_SGType2)
!           IF (present(OldPara)) THEN
!             !write(6,*) 'OldPara in Rec_WrhonD:',OldPara
!             CALL get_Tabiq_Tabil_FROM_iq(nDval_SG2,nDl_SG2,          &
!                                     i_SG,iq_SG,iq,BasisnD%para_SGType2,&
!                                     OldPara,err_sub=err_sub)
!           ELSE
!             CALL get_Tabiq_Tabil_FROM_iq(nDval_SG2,nDl_SG2,          &
!                                     i_SG,iq_SG,iq,BasisnD%para_SGType2,&
!                                     err_sub=err_sub)
!           END IF
!           IF (err_sub /= 0) STOP 'Rec_WrhonD'

           WrhonD = BasisnD%WeightSG(i_SG)
           DO ib=1,BasisnD%nb_basis
             L     = nDl_SG2(ib)
             iq_ib = nDval_SG2(ib)
             WrhonD = WrhonD * Rec_WrhonD(BasisnD%tab_basisPrimSG(L,ib),iq_ib)
           END DO
           IF (debug) write(out_unitp,*) 'iq,i_SG,iq_SG',iq,i_SG,iq_SG
           IF (debug) write(out_unitp,*) 'DP: nDval',nDval_SG2(:)
           IF (debug) write(out_unitp,*) 'WrhonD',iq,WrhonD

         CASE (4) ! Sparse basis (Smolyak 4th implementation)

           IF (present(OldPara)) THEN
             CALL get_Tabiq_Tabil_FROM_iq(tab_iq,tab_l,i_SG,iq_SG,iq,BasisnD%para_SGType2,OldPara,err_sub)
           ELSE
             CALL get_Tabiq_Tabil_FROM_iq(tab_iq,tab_l,i_SG,iq_SG,iq,BasisnD%para_SGType2,err_sub=err_sub)
           END IF
           IF (err_sub /= 0) STOP 'Error in get_Tabiq_Tabil_FROM_iq in Rec_rhonD'

           tab_nq(:) = getbis_tab_nq(tab_l,BasisnD%tab_basisPrimSG)

           CALL calc_nDval_m1(tab_iq,iq_SG,tab_nq,BasisnD%nb_basis)

           WrhonD = BasisnD%WeightSG(i_SG)
           DO ib=1,BasisnD%nb_basis
             WrhonD = WrhonD * Rec_WrhonD(BasisnD%tab_basisPrimSG(tab_l(ib),ib),tab_iq(ib))
           END DO

         CASE DEFAULT
           write(out_unitp,*) ' ERROR in',name_sub
           write(out_unitp,*) ' WRONG SparseGrid_type',BasisnD%SparseGrid_type
           write(out_unitp,*) ' The possibilities are: 0, 1, 2, 4'
           STOP
         END SELECT
       END IF

!     -------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) ' WrhonD',WrhonD
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -------------------------------------------------------

      END FUNCTION Rec_WrhonD


      RECURSIVE SUBROUTINE Rec_tab_iq(tab_iq,BasisnD,iq,iqi)
      USE mod_system
      implicit none

!----- variables for the Basis and quadrature points -----------------
      TYPE (Basis)           :: BasisnD
      integer, intent(in)    :: iq
      integer, intent(inout) :: tab_iq(BasisnD%ndim)
      integer, optional      :: iqi

!------ working variables ---------------------------------
      integer           :: idim,ndim,iqi_loc,Lm
      integer           :: i,iq_SG,nq,L,ib,iq_ib
      real (kind=Rkind) :: WrhonD
      integer           :: nDval(BasisnD%nDindG%ndim)
      integer           :: nDval_SG2(BasisnD%nb_basis)
      integer           :: nDl_SG2(BasisnD%nb_basis)

      integer           :: i_SG
      integer           :: err_sub

!----- for debuging --------------------------------------------------
       logical, parameter :: debug=.FALSE.
!       logical, parameter :: debug=.TRUE.
       character (len=*), parameter :: name_sub = 'Rec_tab_iq'
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'iq',iq
       END IF
!-----------------------------------------------------------
!
       IF (BasisnD%ndim == 0) THEN
         CONTINUE ! no grid point
       ELSE IF (BasisnD%primitive) THEN
         IF (present(iqi)) THEN
           tab_iq(:) = iq + iqi
         ELSE
           tab_iq(:) = iq
         END IF
       ELSE ! BasisnD%nb_basis MUST BE > 0
         IF (BasisnD%nb_basis == 0 ) STOP ' ERROR with packed in Rec_tab_iq!!!'

         SELECT CASE (BasisnD%SparseGrid_type)
         CASE (0) ! Direct product
           CALL calc_nDindex(BasisnD%nDindG,iq,nDval)
           IF (debug) write(out_unitp,*) 'DP: nDval',nDval(:)
           idim = 0
           DO i=1,BasisnD%nb_basis
             ndim = BasisnD%tab_Pbasis(i)%Pbasis%ndim
             CALL Rec_tab_iq(tab_iq(idim+1:idim+ndim),                  &
                                  BasisnD%tab_Pbasis(i)%Pbasis,nDval(i))
             idim = idim+ndim
           END DO

         CASE (1) ! Sparse basis (Smolyak 1st implementation)

           iq_SG = iq
           DO i_SG=1,BasisnD%nb_SG
             nq = get_nq_FROM_basis(BasisnD%tab_PbasisSG(i_SG)%Pbasis)
             IF (iq_SG <= nq) EXIT
             iq_SG = iq_SG - nq
           END DO
           IF (debug) write(out_unitp,*) 'i_SG,iq_SG',i_SG,iq_SG
           CALL Rec_tab_iq(tab_iq,BasisnD%tab_PbasisSG(i_SG)%Pbasis,iq_SG)

         CASE (2) ! Sparse basis (Smolyak 2d  implementation)

           CALL get_Tabiq_Tabil_FROM_iq(nDval_SG2,nDl_SG2,            &
                      i_SG,iq_SG,iq,BasisnD%para_SGType2,err_sub=err_sub)
           IF (err_sub /= 0) STOP 'Rec_tab_iq'

           idim  = 0
           DO ib=1,BasisnD%nb_basis
             L     = nDl_SG2(ib)
             iq_ib = nDval_SG2(ib)
             iqi_loc = 0
             DO Lm=0,L-1
               iqi_loc = iqi_loc + get_nq_FROM_basis(BasisnD%tab_basisPrimSG(Lm,ib))
             END DO
             ndim  = BasisnD%tab_basisPrimSG(L,ib)%ndim
             write(6,*) 'ib,L,idim,ndim',ib,L,idim,ndim
             CALL Rec_tab_iq(tab_iq(idim+1:idim+ndim),                  &
                         BasisnD%tab_basisPrimSG(L,ib),iq_ib,iqi=iqi_loc)
             idim = idim+ndim
           END DO
           IF (debug) write(out_unitp,*) 'iq,i_SG,iq_SG',iq,i_SG,iq_SG
           IF (debug) write(out_unitp,*) 'DP: nDval',nDval_SG2(:)
           IF (debug) write(out_unitp,*) 'WrhonD',iq,WrhonD

         CASE (4) ! Sparse basis (Smolyak 4th implementation)
           write(out_unitp,*) ' ERROR in',name_sub
           STOP 'Rec_tab_iq: SparseGrid_type=4'

         CASE DEFAULT
           write(out_unitp,*) ' ERROR in',name_sub
           write(out_unitp,*) ' WRONG SparseGrid_type',BasisnD%SparseGrid_type
           write(out_unitp,*) ' The possibilities are: 0, 1, 2, 4'
           STOP
         END SELECT
       END IF

!     -------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) ' tab_iq',tab_iq(:)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -------------------------------------------------------

      END SUBROUTINE Rec_tab_iq

!=====================================================================
!
!  calculation of rhoD  as a function of nDval and the basis
!
!=====================================================================
      RECURSIVE FUNCTION Rec_rhonD(BasisnD,iq,OldPara) result(rhonD)
      USE mod_system
      USE mod_basis_BtoG_GtoB_SGType4, ONLY : getbis_tab_nq
      implicit none

!----- variables for the Basis and quadrature points -----------------
      TYPE (Basis)        :: BasisnD
      integer, intent(in) :: iq
      real (kind=Rkind)   :: rhonD
      TYPE (OldParam), intent(inout), optional :: OldPara

!------ working variables ---------------------------------
      integer           :: i,iqi,iq_SG,nq,L,ib,iq_ib
      integer           :: nDval(BasisnD%nDindG%ndim)
      integer           :: nDval_SG2(BasisnD%nb_basis)
      integer           :: nDl_SG2(BasisnD%nb_basis)

      integer             :: tab_nq(BasisnD%nb_basis)
      integer             :: tab_l(BasisnD%nb_basis)
      integer             :: tab_iq(BasisnD%nb_basis)

      integer           :: i_SG
      integer           :: err_sub

!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
      !logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'Rec_rhonD'

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'BasisnD%nb_basis',BasisnD%nb_basis
         write(out_unitp,*) 'BasisnD%SparseGrid_type',BasisnD%SparseGrid_type
         write(out_unitp,*) 'BasisnD%packed_done',BasisnD%packed_done
         write(out_unitp,*) 'iq',iq
         CALL flush_perso(out_unitp)
       END IF
!-----------------------------------------------------------

       IF (BasisnD%ndim == 0) THEN
         rhonD = ONE ! no grid point
       ELSE IF (BasisnD%packed_done) THEN
         rhonD = BasisnD%rho(iq)
       ELSE ! BasisnD%nb_basis MUST BE > 0
         IF (BasisnD%nb_basis == 0 ) STOP ' ERROR with packed in Rec_rhonD!!!'

         SELECT CASE (BasisnD%SparseGrid_type)
         CASE (0) ! Direct product
           CALL calc_nDindex(BasisnD%nDindG,iq,nDval)
           IF (debug) write(out_unitp,*) 'iq',iq,'nDval:',nDval
           rhonD = ONE
           DO i=1,BasisnD%nb_basis
             rhonD = rhonD * Rec_rhonD(BasisnD%tab_Pbasis(i)%Pbasis,nDval(i))
           END DO

         CASE (1) ! Sparse basis (Smolyak 1st implementation)
           iq_SG = iq
           DO i_SG=1,BasisnD%nb_SG
             nq = get_nq_FROM_basis(BasisnD%tab_PbasisSG(i_SG)%Pbasis)
             IF (iq_SG <= nq) EXIT
             iq_SG = iq_SG - nq
           END DO
           rhonD = Rec_rhonD(BasisnD%tab_PbasisSG(i_SG)%Pbasis,iq_SG)

         CASE (2) ! Sparse basis (Smolyak 2d)
           CALL get_Tabiq_Tabil_FROM_iq_old(nDval_SG2,nDl_SG2,          &
                                     i_SG,iq_SG,iq,BasisnD%para_SGType2)

           rhonD = ONE
           DO ib=1,BasisnD%nb_basis
             L     = nDl_SG2(ib)
             iq_ib = nDval_SG2(ib)
             rhonD = rhonD * Rec_rhonD(BasisnD%tab_basisPrimSG(L,ib),iq_ib)
           END DO

         CASE (4) ! Sparse basis (Smolyak 4th implementation)

           IF (present(OldPara)) THEN
             CALL get_Tabiq_Tabil_FROM_iq(tab_iq,tab_l,i_SG,iq_SG,iq,BasisnD%para_SGType2,OldPara,err_sub)
           ELSE
             CALL get_Tabiq_Tabil_FROM_iq(tab_iq,tab_l,i_SG,iq_SG,iq,BasisnD%para_SGType2,err_sub=err_sub)
           END IF
           IF (err_sub /= 0) STOP 'Error in get_Tabiq_Tabil_FROM_iq in Rec_rhonD'

           tab_nq(:) = getbis_tab_nq(tab_l,BasisnD%tab_basisPrimSG)

           CALL calc_nDval_m1(tab_iq,iq_SG,tab_nq,BasisnD%nb_basis)

           rhonD = ONE
           DO ib=1,BasisnD%nb_basis
             rhonD = rhonD * Rec_rhonD(BasisnD%tab_basisPrimSG(tab_l(ib),ib),tab_iq(ib))
           END DO

         CASE DEFAULT
           write(out_unitp,*) ' ERROR in',name_sub
           write(out_unitp,*) ' WRONG SparseGrid_type',BasisnD%SparseGrid_type
           write(out_unitp,*) ' The possibilities are: 0, 1, 2, 4'
           STOP
         END SELECT
       END IF

!     -------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) ' rhonD',rhonD
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
!     -------------------------------------------------------

      END FUNCTION Rec_rhonD

      RECURSIVE FUNCTION Rec_WnD(BasisnD,iq) result(WnD)
      USE mod_system
      implicit none

!----- variables for the Basis and quadrature points -----------------
      TYPE (Basis)        :: BasisnD
      integer, intent(in) :: iq
      real (kind=Rkind)   :: WnD

!------ working variables ---------------------------------
      integer           :: i,iqi,iq_SG,nq,L,ib,iq_ib
      integer           :: nDval(BasisnD%nDindG%ndim)
      integer           :: nDval_SG2(BasisnD%nb_basis)
      integer           :: nDl_SG2(BasisnD%nb_basis)

      integer           :: i_SG
      integer           :: err_sub

!----- for debuging --------------------------------------------------
      logical, parameter :: debug=.FALSE.
!     logical, parameter :: debug=.TRUE.
      character (len=*), parameter :: name_sub = 'Rec_WnD'
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'BasisnD%nb_basis',BasisnD%nb_basis
       END IF
!-----------------------------------------------------------
!
       IF (BasisnD%ndim == 0) THEN
         WnD = ONE ! no grid point
       ELSE IF (BasisnD%packed_done) THEN
         WnD = BasisnD%w(iq)
       ELSE ! BasisnD%nb_basis MUST BE > 0
         IF (BasisnD%nb_basis == 0 ) STOP ' ERROR with packed in Rec_WnD!!!'

         SELECT CASE (BasisnD%SparseGrid_type)
         CASE (0) ! Direct product
           CALL calc_nDindex(BasisnD%nDindG,iq,nDval)
           WnD = ONE
           DO i=1,BasisnD%nb_basis
             WnD = WnD * Rec_WnD(BasisnD%tab_Pbasis(i)%Pbasis,nDval(i))
           END DO

         CASE (1) ! Sparse basis (Smolyak 1st implementation)
           iq_SG = iq
           DO i_SG=1,BasisnD%nb_SG
             nq = get_nq_FROM_basis(BasisnD%tab_PbasisSG(i_SG)%Pbasis)
             IF (iq_SG <= nq) EXIT
             iq_SG = iq_SG - nq
           END DO
           WnD = BasisnD%WeightSG(i_SG) *                               &
                 Rec_WnD(BasisnD%tab_PbasisSG(i_SG)%Pbasis,iq_SG)

         CASE (2) ! Sparse basis (Smolyak 2d and 4th implementation)
           CALL get_Tabiq_Tabil_FROM_iq_old(nDval_SG2,nDl_SG2,          &
                                     i_SG,iq_SG,iq,BasisnD%para_SGType2)

!           CALL get_Tabiq_Tabil_FROM_iq(nDval_SG2,nDl_SG2,            &
!                      i_SG,iq_SG,iq,BasisnD%para_SGType2,err_sub=err_sub)
!           IF (err_sub /= 0) STOP 'Rec_WnD'

           WnD = BasisnD%WeightSG(i_SG)
           DO ib=1,BasisnD%nb_basis
             L     = nDl_SG2(ib)
             iq_ib = nDval_SG2(ib)
             WnD = WnD * Rec_WnD(BasisnD%tab_basisPrimSG(L,ib),iq_ib)
           END DO

         CASE (4) ! Sparse basis (Smolyak 2d and 4th implementation)

           CALL get_Tabiq_Tabil_FROM_iq(nDval_SG2,nDl_SG2,            &
                      i_SG,iq_SG,iq,BasisnD%para_SGType2,err_sub=err_sub)
           IF (err_sub /= 0) STOP 'Rec_WnD'

           WnD = BasisnD%WeightSG(i_SG)
           DO ib=1,BasisnD%nb_basis
             L     = nDl_SG2(ib)
             iq_ib = nDval_SG2(ib)
             WnD = WnD * Rec_WnD(BasisnD%tab_basisPrimSG(L,ib),iq_ib)
           END DO

         CASE DEFAULT
           write(out_unitp,*) ' ERROR in',name_sub
           write(out_unitp,*) ' WRONG SparseGrid_type',BasisnD%SparseGrid_type
           write(out_unitp,*) ' The possibilities are: 0, 1, 2, 4'
           STOP
         END SELECT

       END IF

!     -------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) ' WnD',WnD
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -------------------------------------------------------

      END FUNCTION Rec_WnD

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      RECURSIVE SUBROUTINE Rec_x(x,BasisnD,iq,OldPara)
      USE mod_system
      USE mod_basis_BtoG_GtoB_SGType4, ONLY : getbis_tab_nq
      USE mod_nDindex
      IMPLICIT NONE


!----- variables for the Basis and quadrature points -----------------
      TYPE (Basis)         :: BasisnD
      integer              :: iq
      real (kind=Rkind)    :: x(:)
      TYPE (OldParam), intent(inout), optional :: OldPara

!------ working variables ---------------------------------
      integer             :: iq_SG,nq
      integer             :: i,j,iqi,i0,i1,ib,L,iq_ib
      integer             :: nDval(BasisnD%nDindG%ndim)
      integer             :: nDval_SG2(BasisnD%nb_basis)
      integer             :: nDl_SG2(BasisnD%nb_basis)
      integer             :: tab_nq(BasisnD%nb_basis)
      integer             :: tab_l(BasisnD%nb_basis)
      integer             :: tab_iq(BasisnD%nb_basis)

      integer           :: i_SG
      integer           :: err_sub


!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'Rec_x'

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         CALL flush_perso(out_unitp)
       END IF
!-----------------------------------------------------------
      !write(6,*) ' in ',name_sub,'iq',iq


       IF (BasisnD%ndim == 0) THEN
         CONTINUE ! no grid point
       ELSE IF (BasisnD%packed_done) THEN
         x(:) = BasisnD%x(:,iq)
       ELSE ! BasisnD%nb_basis MUST BE > 0
         IF (BasisnD%nb_basis == 0 ) STOP ' ERROR with packed in Rec_x!!!'

         SELECT CASE (BasisnD%SparseGrid_type)
         CASE (0,3) ! Direct product and SG for type 21
           CALL calc_nDindex(BasisnD%nDindG,iq,nDval)
           i0 = 0
           DO i=1,BasisnD%nb_basis
             i1 = i0 + BasisnD%tab_Pbasis(i)%Pbasis%ndim
             CALL Rec_x(x(i0+1:i1),BasisnD%tab_Pbasis(i)%Pbasis,nDval(i))
             i0 = i1
           END DO

         CASE (1) ! Sparse basis (Smolyak 1st implementation)
           iq_SG = iq
           DO i_SG=1,BasisnD%nb_SG
             nq = get_nq_FROM_basis(BasisnD%tab_PbasisSG(i_SG)%Pbasis)
             IF (iq_SG <= nq) EXIT
             iq_SG = iq_SG - nq
           END DO
           CALL Rec_x(x,BasisnD%tab_PbasisSG(i_SG)%Pbasis,iq_SG)

!         CASE (2) ! Sparse basis (Smolyak 2d implementation)
!
!           IF (present(OldPara)) THEN
!             !write(6,*) 'OldPara in Rec_x:',OldPara
!             CALL get_Tabiq_Tabil_FROM_iq(nDval_SG2,nDl_SG2,          &
!                                     i_SG,iq_SG,iq,BasisnD%para_SGType2,&
!                                     OldPara,err_sub=err_sub)
!           ELSE
!             CALL get_Tabiq_Tabil_FROM_iq(nDval_SG2,nDl_SG2,          &
!                                     i_SG,iq_SG,iq,BasisnD%para_SGType2,&
!                                     err_sub=err_sub)
!           END IF
!           IF (err_sub /= 0) STOP 'Rec_x'
!
!           i0 = 0
!           DO ib=1,BasisnD%nb_basis
!             L     = nDl_SG2(ib)
!             iq_ib = nDval_SG2(ib)
!
!             i1 = i0 + BasisnD%tab_basisPrimSG(L,ib)%ndim
!             CALL Rec_x(x(i0+1:i1),BasisnD%tab_basisPrimSG(L,ib),iq_ib)
!             i0 = i1
!           END DO
!           !write(out_unitp,*) ' Qbasis',iq,x

         CASE (2,4) ! Sparse basis (Smolyak 4th implementation)

           IF (present(OldPara)) THEN
             CALL get_Tabiq_Tabil_FROM_iq(tab_iq,tab_l,i_SG,iq_SG,iq,BasisnD%para_SGType2,OldPara,err_sub)
           ELSE
             CALL get_Tabiq_Tabil_FROM_iq(tab_iq,tab_l,i_SG,iq_SG,iq,BasisnD%para_SGType2,err_sub=err_sub)
           END IF
           IF (err_sub /= 0) STOP 'Error in get_Tabiq_Tabil_FROM_iq'

           tab_nq(:) = getbis_tab_nq(tab_l,BasisnD%tab_basisPrimSG)

           CALL calc_nDval_m1(tab_iq,iq_SG,tab_nq,BasisnD%nb_basis)

           i0 = 0
           DO ib=1,BasisnD%nb_basis
             i1 = i0 +BasisnD%tab_basisPrimSG(tab_l(ib),ib)%ndim
             CALL Rec_x(x(i0+1:i1),BasisnD%tab_basisPrimSG(tab_l(ib),ib),tab_iq(ib))
             i0 = i1
           END DO


         CASE DEFAULT
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) ' WRONG SparseGrid_type',BasisnD%SparseGrid_type
           write(out_unitp,*) ' The possibilities are: 0, 1, 2, 4, 3'
           STOP
         END SELECT

       END IF

      !write(6,*) ' out ',name_sub,'iq,x',iq,x

!     -------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) ' Qbasis, x',x
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
!     -------------------------------------------------------

      END SUBROUTINE Rec_x
  SUBROUTINE Rec_x_SG4(x,tab_ba,tab_l,nDind_DPG,iq,err_sub)
  USE mod_system
  USE mod_nDindex
  IMPLICIT NONE

  real (kind=Rkind),               intent(inout)          :: x(:)

  TYPE(basis),                     intent(in)             :: tab_ba(0:,:) ! tab_ba(0:L,D)
  integer,                         intent(in)             :: tab_l(:)
  TYPE (Type_nDindex),             intent(in)             :: nDind_DPG    ! multidimensional DP index (nb_SG)

  integer,                         intent(in)             :: iq
  integer,                         intent(inout)          :: err_sub


  !------ working variables ---------------------------------
  integer           :: ib,i0,i1,L,iq_ib
  integer           :: tab_iq(size(tab_l))

  !----- for debuging --------------------------------------------------
  logical, parameter :: debug = .FALSE.
  !logical, parameter :: debug = .TRUE.
  character (len=*), parameter :: name_sub = 'Rec_x_SG4'

  !-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    CALL Write_nDindex(nDind_DPG,' in ' // name_sub // ': ')
    CALL flush_perso(out_unitp)
  END IF
  !-----------------------------------------------------------

  CALL calc_nDindex(nDind_DPG,iq,tab_iq,err_sub)
  IF (err_sub /= 0) THEN
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) '  from nDind_DPG'
    STOP 'calc_nDindex'
  END IF

  i0 = 0
  DO ib=1,size(tab_l)
    L     = tab_l(ib)
    iq_ib = tab_iq(ib)

    i1 = i0 + tab_ba(L,ib)%ndim
    CALL Rec_x(x(i0+1:i1),tab_ba(L,ib),iq_ib)
    i0 = i1
  END DO


! -------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*)
    write(out_unitp,*) ' Qbasis, x',x
    write(out_unitp,*)
    write(out_unitp,*) 'END ',name_sub
    CALL flush_perso(out_unitp)
  END IF
! -------------------------------------------------------

  END SUBROUTINE Rec_x_SG4

  SUBROUTINE Rec_Qact(Qact,BasisnD,iq,mole,OldPara)
  USE mod_system
  USE mod_Coord_KEO

  IMPLICIT NONE


  !- variables for the Basis and quadrature points -----------------
  TYPE (Basis) :: BasisnD
  integer      :: iq
  TYPE (OldParam), intent(inout), optional :: OldPara


  !- for the zmatrix  --------------------------------------
  TYPE (zmatrix) :: mole

  real (kind=Rkind) :: Qact(BasisnD%ndim)

  !-- working variables ---------------------------------
  integer           :: j_act,j
  real (kind=Rkind) :: x(BasisnD%ndim)

  !- for debuging --------------------------------------------------
  logical, parameter :: debug = .FALSE.
  !logical, parameter :: debug = .TRUE.
  character (len=*), parameter :: name_sub = 'Rec_Qact'
  !-------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*) 'nb_act1',mole%nb_act1
    write(out_unitp,*) 'iq',iq
    write(out_unitp,*) 'BasisnD%ndim',BasisnD%ndim
    CALL flush_perso(out_unitp)
   END IF
  !-------------------------------------------------------


  IF (present(OldPara)) THEN
    CALL Rec_x(x,BasisnD,iq,OldPara)
  ELSE
    CALL Rec_x(x,BasisnD,iq)
  END IF

  DO j=1,BasisnD%ndim
    j_act = mole%ActiveTransfo%list_QdynTOQact(BasisnD%iQdyn(j))
    Qact(j_act) = x(j)
  END DO

  ! -------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*)
    write(out_unitp,*) ' Qact',Qact
    write(out_unitp,*)
    write(out_unitp,*) 'END ',name_sub
    CALL flush_perso(out_unitp)
  END IF
  ! -------------------------------------------------------

  END SUBROUTINE Rec_Qact

  SUBROUTINE Rec_Qact_SG4(Qact,tab_ba,tab_l,nDind_DPG,iq,mole,err_sub)
  USE mod_system
  USE mod_nDindex
  USE mod_Coord_KEO
  IMPLICIT NONE

  real (kind=Rkind),               intent(inout)          :: Qact(:)
  TYPE(basis),                     intent(in)             :: tab_ba(0:,:) ! tab_ba(0:L,D)
  integer,                         intent(in)             :: tab_l(:)
  TYPE (Type_nDindex),             intent(in)             :: nDind_DPG    ! multidimensional DP index
  integer,                         intent(in)             :: iq
  TYPE (zmatrix),                  intent(in)             :: mole
  integer,                         intent(inout)          :: err_sub

!------ working variables ---------------------------------
      integer                        :: ndim,i,i_act,ib,ix
      real (kind=Rkind), allocatable :: x(:)

!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
      ! logical, parameter :: debug = .TRUE.
   character (len=*), parameter :: name_sub = 'Rec_Qact_SG4'
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'nb_act1',mole%nb_act1
         write(out_unitp,*) 'iq',iq
         write(out_unitp,*) 'tab_l',tab_l
         CALL Write_nDindex(nDind_DPG,' in ' // name_sub // ': ')
         CALL flush_perso(out_unitp)
       END IF
!-----------------------------------------------------------
       ndim = 0
       DO i=1,size(tab_l)
         ndim = ndim + tab_ba(tab_l(i),i)%ndim
       END DO

       CALL alloc_NParray(x,(/ndim/),'x',name_sub)

       CALL Rec_x_SG4(x,tab_ba,tab_l,nDind_DPG,iq,err_sub)
       IF (err_sub /= 0) THEN
         write(out_unitp,*) ' ERROR in ',name_sub
         write(out_unitp,*) '  from nDind_DPG'
         STOP 'Rec_x_SG4'
       END IF


       ix = 0
       DO ib=1,size(tab_l)
         DO i=1,tab_ba(tab_l(ib),ib)%ndim
           i_act = mole%ActiveTransfo%list_QdynTOQact(tab_ba(tab_l(ib),ib)%iQdyn(i))
           ix = ix + 1
           Qact(i_act) = x(ix)
         END DO
       END DO

       CALL dealloc_NParray(x,'x',name_sub)


!     -------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) ' Qact',Qact
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
!     -------------------------------------------------------

      END SUBROUTINE Rec_Qact_SG4

  SUBROUTINE Rec_Qact_SG4_with_Tab_iq(Qact,tab_ba,tab_l,tab_iq,mole,err_sub)
  USE mod_system
  USE mod_Coord_KEO
  IMPLICIT NONE

  real (kind=Rkind),               intent(inout)          :: Qact(:)
  TYPE(basis),                     intent(in)             :: tab_ba(0:,:) ! tab_ba(0:L,D)
  integer,                         intent(in)             :: tab_l(:)
  integer,                         intent(in)             :: tab_iq(:)
  TYPE (zmatrix),                  intent(in)             :: mole
  integer,                         intent(inout)          :: err_sub

!------ working variables ---------------------------------
  integer                        :: ndim,i,i_act,ib
  real (kind=Rkind), allocatable :: x(:)

!----- for debuging --------------------------------------------------
  logical, parameter :: debug = .FALSE.
  !logical, parameter :: debug = .TRUE.
  character (len=*), parameter :: name_sub = 'Rec_Qact_SG4_with_Tab_iq'
!-----------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
    write(out_unitp,*) 'nb_act1',mole%nb_act1
    write(out_unitp,*) 'tab_l',tab_l
    write(out_unitp,*) 'tab_iq',tab_iq
    CALL flush_perso(out_unitp)
  END IF
!-----------------------------------------------------------

  DO ib=1,size(tab_l)
    ndim = tab_ba(tab_l(ib),ib)%ndim
    CALL alloc_NParray(x,(/ndim/),'x',name_sub)

    CALL Rec_x(x,tab_ba(tab_l(ib),ib),tab_iq(ib))

    DO i=1,tab_ba(tab_l(ib),ib)%ndim
      i_act = mole%ActiveTransfo%list_QdynTOQact(tab_ba(tab_l(ib),ib)%iQdyn(i))
      Qact(i_act) = x(i)
    END DO

    CALL dealloc_NParray(x,'x',name_sub)

  END DO


! -------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*)
    write(out_unitp,*) ' Qact',Qact
    write(out_unitp,*)
    write(out_unitp,*) 'END ',name_sub
    CALL flush_perso(out_unitp)
  END IF
! -------------------------------------------------------

  END SUBROUTINE Rec_Qact_SG4_with_Tab_iq

  SUBROUTINE Rec_Qact_SG4_with_Tab_iq_old(Qact,tab_ba,tab_l,tab_iq,mole,err_sub)
  USE mod_system
  USE mod_Coord_KEO
  IMPLICIT NONE

  real (kind=Rkind),               intent(inout)          :: Qact(:)
  TYPE(basis),                     intent(in)             :: tab_ba(0:,:) ! tab_ba(0:L,D)
  integer,                         intent(in)             :: tab_l(:)
  integer,                         intent(in)             :: tab_iq(:)
  TYPE (zmatrix),                  intent(in)             :: mole
  integer,                         intent(inout)          :: err_sub

!------ working variables ---------------------------------
      integer                        :: ndim,i,i_act,ib,ix,i0,i1
      real (kind=Rkind), allocatable :: x(:)

!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
      ! logical, parameter :: debug = .TRUE.
   character (len=*), parameter :: name_sub = 'Rec_Qact_SG4_with_Tab_iq_old'
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'nb_act1',mole%nb_act1
         write(out_unitp,*) 'tab_l',tab_l
         write(out_unitp,*) 'tab_iq',tab_iq
         CALL flush_perso(out_unitp)
       END IF
!-----------------------------------------------------------
       ndim = 0
       DO i=1,size(tab_l)
         ndim = ndim + tab_ba(tab_l(i),i)%ndim
       END DO

       CALL alloc_NParray(x,(/ndim/),'x',name_sub)

       i0 = 0
       DO ib=1,size(tab_l)
         i1 = i0 + tab_ba(tab_l(ib),ib)%ndim
         CALL Rec_x(x(i0+1:i1),tab_ba(tab_l(ib),ib),tab_iq(ib))
         i0 = i1

       END DO

       ix = 0
       DO ib=1,size(tab_l)
         DO i=1,tab_ba(tab_l(ib),ib)%ndim
           i_act = mole%ActiveTransfo%list_QdynTOQact(tab_ba(tab_l(ib),ib)%iQdyn(i))
           ix = ix + 1
           Qact(i_act) = x(ix)
         END DO
       END DO

       CALL dealloc_NParray(x,'x',name_sub)


!     -------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) ' Qact',Qact
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
        CALL flush_perso(out_unitp)
      END IF
!     -------------------------------------------------------

      END SUBROUTINE Rec_Qact_SG4_with_Tab_iq_old

!=====================================================================
!
!  calculation of d0bnD
!  calculation of d0cbnD (in complex)
!
!=====================================================================
      !!@description: calculation of d0bnD, calculation of d0cbnD (in complex)
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      RECURSIVE FUNCTION Rec_d0bnD(BasisnD,iq,ib) result(d0bnD)
      USE mod_system
      implicit none

!----- variables for the Basis and quadrature points -----------------
      TYPE (Basis)         :: BasisnD
      integer, intent(in)  :: iq
      integer, intent(in)  :: ib
      real (kind=Rkind)    :: d0bnD

!------ working variables ---------------------------------
      integer       :: i_SG,iq_SG,nq
      integer       :: i
      integer :: nDvalG(BasisnD%nDindG%ndim)
      integer :: nDvalB(BasisnD%nDindB%ndim)


!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'Rec_d0bnD'
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
       END IF
!-----------------------------------------------------------

       IF (BasisnD%packed_done) THEN
         d0bnD = BasisnD%dnRGB%d0(iq,ib)
       ELSE ! BasisnD%nb_basis MUST BE > 0
         IF (BasisnD%nb_basis == 0 ) STOP ' ERROR with packed in Rec_d0bnD!!!'

         SELECT CASE (BasisnD%SparseGrid_type)
         CASE (0,3) ! Direct product
           CALL calc_nDindex(BasisnD%nDindG,iq,nDvalG)
           CALL calc_nDindex(BasisnD%nDindB,ib,nDvalB)
           d0bnD   = ONE
           DO i=1,BasisnD%nb_basis
             d0bnD = d0bnD * Rec_d0bnD(BasisnD%tab_Pbasis(i)%Pbasis,    &
                                       nDvalG(i),nDvalB(i))
           END DO

         CASE (1) ! Sparse basis (Smolyak 1st implementation)
           iq_SG = iq
           DO i_SG=1,BasisnD%nb_SG
             nq = get_nq_FROM_basis(BasisnD%tab_PbasisSG(i_SG)%Pbasis)
             IF (iq_SG <= nq) EXIT
             iq_SG = iq_SG - nq
           END DO
           d0bnD = Rec_d0bnD(BasisnD%tab_PbasisSG(i_SG)%Pbasis,iq_SG,ib)

         CASE (2) ! Sparse basis (Smolyak 2d  implementation)
           write(out_unitp,*) ' ERROR in',name_sub
           STOP 'Rec_d0bnD: SparseGrid_type=2'

         CASE (4) ! Sparse basis (Smolyak 4th implementation)
           write(out_unitp,*) ' ERROR in',name_sub
           STOP 'Rec_d0bnD: SparseGrid_type=4'

         CASE DEFAULT
           write(out_unitp,*) ' ERROR in ',name_sub
           write(out_unitp,*) ' WRONG SparseGrid_type',BasisnD%SparseGrid_type
           write(out_unitp,*) ' The possibilities are: 0, 1, 2, 4, 3 (type 21)'
           STOP
         END SELECT

       END IF

!     -------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) ' d0bnD',d0bnD
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -------------------------------------------------------

      END FUNCTION Rec_d0bnD
      RECURSIVE FUNCTION Rec_d0bnD_AT_Q(BasisnD,ib,Qbasis) result(d0bnD)
      USE mod_system
      implicit none

!----- variables for the Basis and quadrature points -----------------
      TYPE (Basis), intent(inout)   :: BasisnD
      real (kind=Rkind), intent(in) :: Qbasis(BasisnD%ndim)
      integer, intent(in)           :: ib
      real (kind=Rkind)             :: d0bnD


!------ working variables ---------------------------------
      integer :: i,i0,i1
      integer :: nDvalB(BasisnD%nDindB%ndim)


!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'Rec_d0bnD_AT_Q'
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
       END IF
!-----------------------------------------------------------

       IF (BasisnD%primitive) THEN
         d0bnD = d0b_OF_primitive_basis_AT_Q(BasisnD,ib,Qbasis)
       ELSE ! BasisnD%nb_basis MUST BE > 0
         IF (BasisnD%nb_basis == 0 ) STOP ' ERROR with packed in Rec_d0bnD_AT_Q!!!'
         CALL calc_nDindex(BasisnD%nDindB,ib,nDvalB)
         d0bnD   = ONE
         i0 = 0
         DO i=1,BasisnD%nb_basis
           i1 = i0 + BasisnD%tab_Pbasis(i)%Pbasis%ndim
           d0bnD = d0bnD * Rec_d0bnD_AT_Q(BasisnD%tab_Pbasis(i)%Pbasis, &
                                          nDvalB(i),Qbasis(i0+1:i1))
           i0 = i1
         END DO

       END IF

!     -------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) ' d0bnD',d0bnD
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -------------------------------------------------------

      END FUNCTION Rec_d0bnD_AT_Q
      RECURSIVE FUNCTION Rec_d0cbnD(BasisnD,iq,ib) result(d0cbnD)
      USE mod_system
      implicit none

!----- variables for the Basis and quadrature points -----------------
      TYPE (Basis)              :: BasisnD
      integer, intent(in)       :: iq
      integer, intent(in)       :: ib
      complex (kind=Rkind)      :: d0cbnD

!------ working variables ---------------------------------
      integer       :: i_SG,iq_SG,nq
      integer       :: i
      integer :: nDvalG(BasisnD%nDindG%ndim)
      integer :: nDvalB(BasisnD%nDindB%ndim)

!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'Rec_d0cbnD'
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
       END IF
!-----------------------------------------------------------

       IF (BasisnD%packed_done) THEN
         d0cbnD = BasisnD%dnCGB%d0(iq,ib)
       ELSE ! BasisnD%nb_basis MUST BE > 0
         IF (BasisnD%nb_basis == 0 ) STOP ' ERROR with packed in Rec_d0cbnD!!!'

         SELECT CASE (BasisnD%SparseGrid_type)
         CASE (0) ! Direct product
           CALL calc_nDindex(BasisnD%nDindG,iq,nDvalG)
           CALL calc_nDindex(BasisnD%nDindB,ib,nDvalB)
           d0cbnD   = CONE
           DO i=1,BasisnD%nb_basis
             d0cbnD = d0cbnD * Rec_d0cbnD(BasisnD%tab_Pbasis(i)%Pbasis, &
                                          nDvalG(i),nDvalB(i))
           END DO

         CASE (1) ! Sparse basis (Smolyak 1st implementation)
           iq_SG = iq
           DO i_SG=1,BasisnD%nb_SG
             nq = get_nq_FROM_basis(BasisnD%tab_PbasisSG(i_SG)%Pbasis)
             IF (iq_SG <= nq) EXIT
             iq_SG = iq_SG - nq
           END DO
           d0cbnD = Rec_d0cbnD(BasisnD%tab_PbasisSG(i_SG)%Pbasis,iq_SG,ib)

         CASE (2) ! Sparse basis (Smolyak 2d  implementation)
           write(out_unitp,*) ' ERROR in',name_sub
           STOP 'Rec_d0cbnD: SparseGrid_type=2'

         CASE (4) ! Sparse basis (Smolyak 4th implementation)
           write(out_unitp,*) ' ERROR in',name_sub
           STOP 'Rec_d0cbnD: SparseGrid_type=4'

         CASE DEFAULT
           write(out_unitp,*) ' ERROR in',name_sub
           write(out_unitp,*) ' WRONG SparseGrid_type',BasisnD%SparseGrid_type
           write(out_unitp,*) ' The possibilities are: 0, 1, 2, 4'
           STOP
         END SELECT

       END IF

!     -------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) ' d0bnD',d0cbnD
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -------------------------------------------------------

      END FUNCTION Rec_d0cbnD
!=====================================================================
!
!  calculation of d0bnD d1bnD(:) and d2bnD(:,:)
!
!=====================================================================
      !!@description: calculation of d0bnD d1bnD(:) and d2bnD(:,:)
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      RECURSIVE SUBROUTINE Rec_d0d1d2bnD(d0b,d1b,d2b,BasisnD,iq,ib)

      USE mod_system
      implicit none

!----- variables for the Basis and quadrature points -----------------
      TYPE (Basis) :: BasisnD

      integer           :: iq,ib,ndim
      real (kind=Rkind) :: d0b
      real (kind=Rkind) :: d1b(BasisnD%ndim)
      real (kind=Rkind) :: d2b(BasisnD%ndim,BasisnD%ndim)

!------ working variables ---------------------------------
      integer           :: L,iq_SG,nq
      integer           :: i,j,i0,i1
      integer           :: iqi,ibi,ndimi
      real (kind=Rkind) :: d0bi
      real (kind=Rkind), allocatable :: d1bi(:)
      real (kind=Rkind), allocatable :: d2bi(:,:)
      integer :: nDvalG(BasisnD%nDindG%ndim)
      integer :: nDvalB(BasisnD%nDindB%ndim)

      integer           :: nDval_SG2(BasisnD%nb_basis)
      integer           :: nDl_SG2(BasisnD%nb_basis)

      integer           :: i_SG = 0
      integer           :: err_sub

!----- for debuging --------------------------------------------------
       integer :: err_mem,memory
       logical, parameter :: debug = .FALSE.
       !logical, parameter :: debug = .TRUE.
       character (len=*), parameter :: name_sub = 'Rec_d0d1d2bnD'
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'iq,ib',iq,ib
         !write(out_unitp,*) '--------------------------'
         !CALL RecWrite_basis(BasisnD,write_all=.TRUE.)
         !write(out_unitp,*) '--------------------------'
       END IF
!-----------------------------------------------------------

       IF (ib < 0 .OR. ib > BasisnD%nb) THEN
         write(out_unitp,*) ' ERROR in',name_sub
         write(out_unitp,*) ' Wrong range of ib',ib
         write(out_unitp,*) ' nb',BasisnD%nb
         STOP
       END IF

       IF (BasisnD%packed_done) THEN
         d0b      = BasisnD%dnRGB%d0(iq,ib)

         IF (associated(BasisnD%dnRGB%d1) .AND. associated(BasisnD%dnRGB%d2)) THEN
           d1b(:)   = BasisnD%dnRGB%d1(iq,ib,:)
           d2b(:,:) = BasisnD%dnRGB%d2(iq,ib,:,:)
         ELSE
           DO i=1,BasisnD%ndim
             d1b(i) = dot_product(BasisnD%dnRGB%d0(iq,:),BasisnD%dnRBB%d1(:,ib,i))
           END DO
           DO i=1,BasisnD%ndim
           DO j=1,BasisnD%ndim
             d2b(i,j) = dot_product(BasisnD%dnRGB%d0(iq,:),BasisnD%dnRBB%d2(:,ib,i,j))
           END DO
           END DO
         END IF
       ELSE ! BasisnD%nb_basis MUST BE > 0
         IF (BasisnD%nb_basis == 0 ) STOP ' ERROR with packed Rec_d0d1d2bnD!!!'

         SELECT CASE (BasisnD%SparseGrid_type)
         CASE (0) ! Direct product
           CALL calc_nDindex(BasisnD%nDindG,iq,nDvalG)
           CALL calc_nDindex(BasisnD%nDindB,ib,nDvalB)
           d0b      = ONE
           d1b(:)   = ONE
           d2b(:,:) = ONE
           i0 = 0
           ndim = BasisnD%ndim
           DO i=1,BasisnD%nb_basis
             i1 = i0 + BasisnD%tab_Pbasis(i)%Pbasis%ndim
             iqi   = nDvalG(i)
             ibi   = nDvalB(i)
             ndimi = BasisnD%tab_Pbasis(i)%Pbasis%ndim

             CALL alloc_NParray(d1bi,(/ ndimi /),      'd1bi','name_sub')
             CALL alloc_NParray(d2bi,(/ ndimi,ndimi /),'d2bi','name_sub')

             CALL Rec_d0d1d2bnD(d0bi,d1bi,d2bi,BasisnD%tab_Pbasis(i)%Pbasis,iqi,ibi)
             ! no derivative
             d0b = d0b * d0bi
             ! first derivatives
             d1b(1:i0)      = d1b(1:i0)      * d0bi
             d1b(i0+1:i1)   = d1b(i0+1:i1)   * d1bi(:)
             d1b(i1+1:ndim) = d1b(i1+1:ndim) * d0bi
             ! second derivatives
             d2b(1:i0,1:i0)      = d2b(1:i0,1:i0)      * d0bi
             d2b(1:i0,i1+1:ndim) = d2b(1:i0,i1+1:ndim) * d0bi
             DO j=1,i0
               d2b(j,i0+1:i1)    = d2b(j,i0+1:i1)      * d1bi(:)
             END DO

             DO j=1,ndimi
               d2b(i0+j,1:i0)      = d2b(i0+j,1:i0)       * d1bi(j)
               d2b(i0+j,i1+1:ndim) = d2b(i0+j,i1+1:ndim)  * d1bi(j)
             END DO
             d2b(i0+1:i1,i0+1:i1)  = d2b(i0+1:i1,i0+1:i1) * d2bi(:,:)

             d2b(i1+1:ndim,1:i0)      = d2b(i1+1:ndim,1:i0)      * d0bi
             DO j=i1+1,ndim
               d2b(j,i0+1:i1)         = d2b(j,i0+1:i1)           * d1bi(:)
             END DO
             d2b(i1+1:ndim,i1+1:ndim) = d2b(i1+1:ndim,i1+1:ndim) * d0bi

             CALL dealloc_NParray(d1bi,'d1bi','name_sub')
             CALL dealloc_NParray(d2bi,'d2bi','name_sub')

             i0 = i1
           END DO

         CASE (1) ! Sparse basis (Smolyak 1st implementation)
           iq_SG = iq
           DO i_SG=1,BasisnD%nb_SG
             nq = get_nq_FROM_basis(BasisnD%tab_PbasisSG(i_SG)%Pbasis)
             IF (iq_SG <= nq) EXIT
             iq_SG = iq_SG - nq
           END DO
           CALL Rec_d0d1d2bnD(d0b,d1b,d2b,                              &
                              BasisnD%tab_PbasisSG(i_SG)%Pbasis,iq_SG,ib)

         CASE (2) ! Sparse basis (Smolyak 2d  implementation)
           CALL calc_nDindex(BasisnD%nDindB,ib,nDvalB)
           write(6,*) 'ib,nDvalB',ib,':',nDvalB

           CALL get_Tabiq_Tabil_FROM_iq_old(nDval_SG2,nDl_SG2,          &
                      i_SG,iq_SG,iq,BasisnD%para_SGType2)

           write(6,*) 'iq,i_SG',iq,i_SG
           write(6,*) 'tab_l',nDl_SG2
           write(6,*) 'tab_iq',nDval_SG2

           d0b      = ONE
           d1b(:)   = ONE
           d2b(:,:) = ONE
           i0 = 0
           ndim = BasisnD%ndim
           DO i=1,BasisnD%nb_basis
             L     = nDl_SG2(i)
             iqi   = nDval_SG2(i)
             ibi   = nDvalB(i)

             write(6,*) 'ibasis,L,iqi,ibi',i,L,iqi,ibi

             ndimi = BasisnD%tab_basisPrimSG(L,ib)%ndim
             i1    = i0 + ndimi

             CALL alloc_NParray(d1bi,(/ ndimi /),      'd1bi',name_sub)
             CALL alloc_NParray(d2bi,(/ ndimi,ndimi /),'d2bi',name_sub)

             CALL Rec_d0d1d2bnD(d0bi,d1bi,d2bi,BasisnD%tab_basisPrimSG(L,ib),iqi,ibi)

             ! no derivative
             d0b = d0b * d0bi
             ! first derivatives
             d1b(1:i0)      = d1b(1:i0)      * d0bi
             d1b(i0+1:i1)   = d1b(i0+1:i1)   * d1bi(:)
             d1b(i1+1:ndim) = d1b(i1+1:ndim) * d0bi
             ! second derivatives
             d2b(1:i0,1:i0)      = d2b(1:i0,1:i0)      * d0bi
             d2b(1:i0,i1+1:ndim) = d2b(1:i0,i1+1:ndim) * d0bi
             DO j=1,i0
               d2b(j,i0+1:i1)    = d2b(j,i0+1:i1) * d1bi(:)
             END DO

             DO j=1,ndimi
               d2b(i0+j,1:i0)      = d2b(i0+j,1:i0)      * d1bi(j)
               d2b(i0+j,i1+1:ndim) = d2b(i0+j,i1+1:ndim) * d1bi(j)
             END DO
             d2b(i0+1:i1,i0+1:i1)  = d2b(i0+1:i1,i0+1:i1) * d2bi(:,:)

             d2b(i1+1:ndim,1:i0)      = d2b(i1+1:ndim,1:i0) * d0bi
             DO j=i1+1,ndim
               d2b(j,i0+1:i1)         = d2b(j,i0+1:i1) * d1bi(:)
             END DO
             d2b(i1+1:ndim,i1+1:ndim) = d2b(i1+1:ndim,i1+1:ndim) * d0bi

             CALL dealloc_NParray(d1bi,'d1bi',name_sub)
             CALL dealloc_NParray(d2bi,'d2bi',name_sub)

             i0 = i1
           END DO


           !write(out_unitp,*) ' ERROR in',name_sub
           !STOP 'Rec_d0d1d2bnD: SparseGrid_type=2'

         CASE (4) ! Sparse basis (Smolyak 4th implementation)
           write(out_unitp,*) ' ERROR in',name_sub
           STOP 'Rec_d0d1d2bnD: SparseGrid_type=4'

         CASE DEFAULT
           write(out_unitp,*) ' ERROR in',name_sub
           write(out_unitp,*) ' WRONG SparseGrid_type',BasisnD%SparseGrid_type
           write(out_unitp,*) ' The possibilities are: 0, 1, 2, 4'
           STOP
         END SELECT

       END IF

!     -------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) ' d0b',d0b
        write(out_unitp,*) ' d1b',d1b
        write(out_unitp,*) ' d2b',d2b
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -------------------------------------------------------


      END SUBROUTINE Rec_d0d1d2bnD

      RECURSIVE SUBROUTINE Rec_d0d1d2cbnD(d0cb,d1cb,d2cb,BasisnD,iq,ib)
      USE mod_system
      implicit none

!----- variables for the Basis and quadrature points -----------------
      TYPE (Basis) :: BasisnD

      integer              :: iq,ib,ndim
      complex (kind=Rkind) :: d0cb
      complex (kind=Rkind) :: d1cb(BasisnD%ndim)
      complex (kind=Rkind) :: d2cb(BasisnD%ndim,BasisnD%ndim)

!------ working variables ---------------------------------
      integer           :: i_SG,iq_SG,nq
      integer           :: i,j,i0,i1
      integer           :: iqi,ibi,ndimi
      complex (kind=Rkind) :: d0cbi
      complex (kind=Rkind), allocatable :: d1cbi(:)
      complex (kind=Rkind), allocatable :: d2cbi(:,:)
      real (kind=Rkind) :: d0bi
      real (kind=Rkind), allocatable :: d1bi(:)
      real (kind=Rkind), allocatable :: d2bi(:,:)
      integer :: nDvalG(BasisnD%nDindG%ndim)
      integer :: nDvalB(BasisnD%nDindB%ndim)

!----- for debuging --------------------------------------------------
       integer :: err_mem,memory
       logical, parameter :: debug = .FALSE.
!       logical, parameter :: debug = .TRUE.
      character (len=*), parameter :: name_sub = 'Rec_d0d1d2cbnD'

!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'iq,ib',iq,ib
       END IF
!-----------------------------------------------------------

       IF (BasisnD%packed_done) THEN
         d0cb      = BasisnD%dnCGB%d0(iq,ib)
         IF (associated(BasisnD%dnCGB%d1) .AND. associated(BasisnD%dnCGB%d2)) THEN
           d1cb(:)   = BasisnD%dnCGB%d1(iq,ib,:)
           d2cb(:,:) = BasisnD%dnCGB%d2(iq,ib,:,:)
         ELSE
           DO i=1,BasisnD%ndim
             d1cb(i) = dot_product(BasisnD%dnCGB%d0(iq,:),BasisnD%dnCBB%d1(:,ib,i))
           END DO
           DO i=1,BasisnD%ndim
           DO j=1,BasisnD%ndim
             d2cb(i,j) = dot_product(BasisnD%dnCGB%d0(iq,:),BasisnD%dnCBB%d2(:,ib,i,j))
           END DO
           END DO
         END IF
       ELSE ! BasisnD%nb_basis MUST BE > 0
         IF (BasisnD%nb_basis == 0 ) STOP ' ERROR with packed in Rec_d0d1d2cbnD!!!'

         SELECT CASE (BasisnD%SparseGrid_type)
         CASE (0) ! Direct product
           CALL calc_nDindex(BasisnD%nDindG,iq,nDvalG)
           CALL calc_nDindex(BasisnD%nDindB,ib,nDvalB)
           d0cb      = CONE
           d1cb(:)   = CONE
           d2cb(:,:) = CONE
           i0 = 0
           ndim = BasisnD%ndim
           DO i=1,BasisnD%nb_basis
             i1 = i0 + BasisnD%tab_Pbasis(i)%Pbasis%ndim
             iqi   = nDvalG(i)
             ibi   = nDvalB(i)
             ndimi = BasisnD%tab_Pbasis(i)%Pbasis%ndim
             IF (BasisnD%tab_Pbasis(i)%Pbasis%cplx) THEN
               CALL alloc_NParray(d1cbi,(/ ndimi /),'d1cbi',name_sub)
               CALL alloc_NParray(d2cbi,(/ ndimi,ndimi /),'d2cbi',name_sub)

               CALL Rec_d0d1d2cbnD(d0cbi,d1cbi,d2cbi,BasisnD%tab_Pbasis(i)%Pbasis,iqi,ibi)
               ! no derivative
               d0cb = d0cb * d0cbi
               ! first derivatives
               d1cb(1:i0)      = d1cb(1:i0) * d0cbi
               d1cb(i0+1:i1)   = d1cb(i0+1:i1) * d1cbi(:)
               d1cb(i1+1:ndim) = d1cb(i1+1:ndim) * d0cbi
               ! second derivatives
               d2cb(1:i0,1:i0)      = d2cb(1:i0,1:i0) * d0cbi
               d2cb(1:i0,i1+1:ndim) = d2cb(1:i0,i1+1:ndim) * d0cbi
               DO j=1,i0
                 d2cb(j,i0+1:i1)    = d2cb(j,i0+1:i1) * d1cbi(:)
               END DO

               DO j=1,ndimi
                 d2cb(i0+j,1:i0)      = d2cb(i0+j,1:i0) * d1cbi(j)
                 d2cb(i0+j,i1+1:ndim) = d2cb(i0+j,i1+1:ndim) * d1cbi(j)
               END DO
               d2cb(i0+1:i1,i0+1:i1)  = d2cb(i0+1:i1,i0+1:i1) * d2cbi(:,:)

               d2cb(i1+1:ndim,1:i0)      = d2cb(i1+1:ndim,1:i0) * d0cbi
               DO j=i1+1,ndim
                 d2cb(j,i0+1:i1)         = d2cb(j,i0+1:i1) * d1cbi(:)
               END DO
               d2cb(i1+1:ndim,i1+1:ndim) = d2cb(i1+1:ndim,i1+1:ndim) * d0cbi

               CALL dealloc_NParray(d1cbi,'d1cbi',name_sub)
               CALL dealloc_NParray(d2cbi,'d2cbi',name_sub)
             ELSE
               CALL alloc_NParray(d1bi,(/ ndimi /),      'd1bi',name_sub)
               CALL alloc_NParray(d2bi,(/ ndimi,ndimi /),'d2bi',name_sub)

               CALL Rec_d0d1d2bnD(d0bi,d1bi,d2bi,BasisnD%tab_Pbasis(i)%Pbasis,iqi,ibi)
               ! no derivative
               d0cb = d0cb * cmplx(d0bi,kind=Rkind)
               ! first derivatives
               d1cb(1:i0)      = d1cb(1:i0) * cmplx(d0bi,kind=Rkind)
               d1cb(i0+1:i1)   = d1cb(i0+1:i1) * cmplx(d1bi(:),kind=Rkind)
               d1cb(i1+1:ndim) = d1cb(i1+1:ndim) * cmplx(d0bi,kind=Rkind)
               ! second derivatives
               d2cb(1:i0,1:i0)      = d2cb(1:i0,1:i0) * cmplx(d0bi,kind=Rkind)
               d2cb(1:i0,i1+1:ndim) = d2cb(1:i0,i1+1:ndim) * cmplx(d0bi,kind=Rkind)
               DO j=1,i0
                 d2cb(j,i0+1:i1)    = d2cb(j,i0+1:i1) * cmplx(d1bi(:),kind=Rkind)
               END DO

               DO j=1,ndimi
                 d2cb(i0+j,1:i0)      = d2cb(i0+j,1:i0) * cmplx(d1bi(j),kind=Rkind)
                 d2cb(i0+j,i1+1:ndim) = d2cb(i0+j,i1+1:ndim) * cmplx(d1bi(j),kind=Rkind)
               END DO
               d2cb(i0+1:i1,i0+1:i1)  = d2cb(i0+1:i1,i0+1:i1) * cmplx(d2bi(:,:),kind=Rkind)

               d2cb(i1+1:ndim,1:i0)      = d2cb(i1+1:ndim,1:i0) * cmplx(d0bi,kind=Rkind)
               DO j=i1+1,ndim
                 d2cb(j,i0+1:i1)         = d2cb(j,i0+1:i1) * cmplx(d1bi(:),kind=Rkind)
               END DO
               d2cb(i1+1:ndim,i1+1:ndim) = d2cb(i1+1:ndim,i1+1:ndim) * cmplx(d0bi,kind=Rkind)

               CALL dealloc_NParray(d1bi,'d1bi',name_sub)
               CALL dealloc_NParray(d2bi,'d2bi',name_sub)
             END IF
             i0 = i1
           END DO

         CASE (1) ! Sparse basis (Smolyak 1st implementation)
           iq_SG = iq
           DO i_SG=1,BasisnD%nb_SG
             nq = get_nq_FROM_basis(BasisnD%tab_PbasisSG(i_SG)%Pbasis)
             IF (iq_SG <= nq) EXIT
             iq_SG = iq_SG - nq
           END DO
           CALL Rec_d0d1d2cbnD(d0cb,d1cb,d2cb,                          &
                              BasisnD%tab_PbasisSG(i_SG)%Pbasis,iq_SG,ib)

         CASE (2) ! Sparse basis (Smolyak 2d  implementation)
           write(out_unitp,*) ' ERROR in',name_sub
           STOP 'Rec_d0d1d2cbnD: SparseGrid_type=2'

         CASE (4) ! Sparse basis (Smolyak 4th implementation)
           write(out_unitp,*) ' ERROR in',name_sub
           STOP 'Rec_d0d1d2cbnD: SparseGrid_type=4'

         CASE DEFAULT
           write(out_unitp,*) ' ERROR in',name_sub
           write(out_unitp,*) ' WRONG SparseGrid_type',BasisnD%SparseGrid_type
           write(out_unitp,*) ' The possibilities are: 0, 1, 2, 4'
           STOP
         END SELECT

       END IF

!     -------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) ' d0cb',d0cb
        write(out_unitp,*) ' d1cb',d1cb
        write(out_unitp,*) ' d2cb',d2cb
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!     -------------------------------------------------------


      END SUBROUTINE Rec_d0d1d2cbnD
      SUBROUTINE calc_d0b(d0b,BasisnD,iq)
      USE mod_system
      IMPLICIT NONE



!----- for the basis set ----------------------------------------------
      TYPE (Basis) :: BasisnD

!------ nD basis set for the k point of the nD quadrature points --
      real (kind=Rkind) :: d0b(BasisnD%nb)

!------ working variables ---------------------------------
      integer       :: ib,iq


!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING calc_d0b'
         write(out_unitp,*) 'nb_ba',BasisnD%nb
       END IF
!-----------------------------------------------------------


!-----------------------------------------------------------
!     test if the basis is real
!     ------------------------------------------------------
      IF (BasisnD%cplx) THEN
         write(out_unitp,*) ' ERROR in calc_d0b'
         write(out_unitp,*) ' You are using the REAL subroutine '
         write(out_unitp,*) ' for the d0b calculation, but the basis is COMPLEX'
         STOP
      END IF
!-----------------------------------------------------------

!-----------------------------------------------------------
!     -----------------------------------------------------
!     transfert de d0b_k
!
      DO ib=1,BasisnD%nb
       d0b(ib) = Rec_d0bnD(BasisnD,iq,ib)
      END DO
!       -----------------------------------------------------

!=====================================================================
!=====================================================================

!     -------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) ' d0b',BasisnD%nb
        CALL Write_Vec(d0b,out_unitp,8)
        write(out_unitp,*)
        write(out_unitp,*) 'END calc_d0b'
      END IF
!     -------------------------------------------------------

      END SUBROUTINE calc_d0b
!=====================================================================
!
!  and d0b(Q,ib) calculations for the nD quadrature point
!  in complex
!
!=====================================================================
      SUBROUTINE calc_d0cb(d0cb,BasisnD,iq)
      USE mod_system
      IMPLICIT NONE

!----- for the basis set ----------------------------------------------
      TYPE (Basis) :: BasisnD

!------ nD basis set for the k point of the nD quadrature points --
      complex (kind=Rkind) :: d0cb(BasisnD%nb)

!------ working variables ---------------------------------
      integer       :: ib,iq

!----- for debuging --------------------------------------------------
      logical, parameter :: debug = .FALSE.
!     logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING calc_d0cb'
         write(out_unitp,*) 'nb_ba',BasisnD%nb
       END IF
!-----------------------------------------------------------

!-----------------------------------------------------------
!     test if the basis is complex
!     ------------------------------------------------------
      IF ( .NOT. BasisnD%cplx) THEN
         write(out_unitp,*) ' ERROR in calc_d0cb'
         write(out_unitp,*) ' You are using the COMPLEX subroutine '
         write(out_unitp,*) ' for the d0b calculation, but the basis is REAL'
         STOP
      END IF
!-----------------------------------------------------------

      DO ib=1,BasisnD%nb
       d0cb(ib) = Rec_d0cbnD(BasisnD,iq,ib)
      END DO
!     -----------------------------------------------------
!-----------------------------------------------------------


!-----------------------------------------------------------
!     -------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) ' d0cb',BasisnD%nb
        CALL Write_Vec(d0cb,out_unitp,8)
        write(out_unitp,*)
        write(out_unitp,*) 'END calc_d0cb'
      END IF
!     -------------------------------------------------------

      END SUBROUTINE calc_d0cb

      !!@description: TODO
      !!@param: TODO
      !!@param: TODO
      !!@param: TODO
      SUBROUTINE nrho_Basis_TO_nhro_Tnum(para_AllBasis,mole)
         TYPE (param_AllBasis), intent(in) :: para_AllBasis
         TYPE (zmatrix), intent(inout)     :: mole


         integer      :: iQbasis,iQact,iQdyn,nrho


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='nrho_Basis_TO_nhro_Tnum'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------
       IF (debug) THEN
         write(out_unitp,*) 'BEGINNING ',name_sub
         write(out_unitp,*) 'asso: nrho of basisnD',                            &
                                 allocated(para_AllBasis%BasisnD%nrho)
         IF (allocated(para_AllBasis%BasisnD%nrho))                    &
                write(out_unitp,*) 'nrho of basisnD',para_AllBasis%BasisnD%nrho
       END IF
!-----------------------------------------------------------


         IF (allocated(para_AllBasis%BasisnD%nrho)) THEN
           DO iQbasis=1,para_AllBasis%BasisnD%ndim
             iQdyn = para_AllBasis%BasisnD%iQdyn(iQbasis)
             iQact = mole%ActiveTransfo%list_QdynTOQact(iQdyn)
             nrho = para_AllBasis%BasisnD%nrho(iQbasis)
             mole%nrho_OF_Qdyn(iQdyn) = nrho
             mole%nrho_OF_Qact(iQact) = nrho
           END DO
         END IF


         IF (allocated(para_AllBasis%Basis2n%nrho)) THEN
           DO iQbasis=1,para_AllBasis%Basis2n%ndim
             iQdyn = para_AllBasis%Basis2n%iQdyn(iQbasis)
             iQact = mole%ActiveTransfo%list_QdynTOQact(iQdyn)
             nrho = para_AllBasis%Basis2n%nrho(iQbasis)
             mole%nrho_OF_Qdyn(iQdyn) = nrho
             mole%nrho_OF_Qact(iQact) = nrho
           END DO
         END IF

         CALL Write_rho(mole)

         IF (debug) THEN
           write(out_unitp,*) 'mole%nrho_OF_Qdyn',mole%nrho_OF_Qdyn
           write(out_unitp,*) 'END  ',name_sub
         END IF

      END SUBROUTINE nrho_Basis_TO_nhro_Tnum


END MODULE mod_basis
