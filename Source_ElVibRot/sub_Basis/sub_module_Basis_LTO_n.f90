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
      MODULE mod_Basis_L_TO_n
      USE mod_system
      IMPLICIT NONE

        PRIVATE

        TYPE Basis_L_TO_n
          integer                     :: L_TO_n_type     = 0        ! for the initalization (after we use the Tab_L_TO_n)

          integer                     :: Lmax            = -1      ! parameter for the number of points of SparseGrid
          !!! relation between L and n: n(L) = A + B * L**expo or n(L) = Tab_L_TO_n(L)

          integer                     :: A                = 1
          integer                     :: B                = 1
          integer                     :: expo             = 1
          integer                     :: C                = 0       !  use nq(L) with basis parameters. ???? with Lexpo_TO_nq > 1
                                                                    ! Then when L>LB add (L-LB)*L_TO_nq_C
          integer                     :: max_n            = huge(1) ! value such n(L)<= max_n
          integer, allocatable        :: Tab_L_TO_n(:)
          integer, allocatable        :: Tab_n_TO_L(:)

        END TYPE Basis_L_TO_n

        INTERFACE assignment (=)
          MODULE PROCEDURE L_TO_n_para2_TO_L_TO_n_para1
        END INTERFACE

        PUBLIC :: Basis_L_TO_n, assignment (=), init_Basis_L_TO_n,      &
                  Write_Basis_L_TO_n, Set_Basis_L_TO_n,                 &
                  alloc_Basis_L_TO_n, dealloc_Basis_L_TO_n,             &
                  Get_n_FROM_Basis_L_TO_n, Get_L_FROM_Basis_L_TO_n,     &
                  check_Basis_L_TO_n

      CONTAINS

      SUBROUTINE alloc_Basis_L_TO_n(L_TO_n_para,Lmax)

      TYPE (Basis_L_TO_n), intent(inout) :: L_TO_n_para
      integer,             intent (in)   :: Lmax


!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='alloc_Basis_L_TO_n'
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------

      CALL dealloc_Basis_L_TO_n(L_TO_n_para)

      IF (Lmax < 0 ) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' Wrong Lmax value',Lmax
        STOP
      END IF

      L_TO_n_para%Lmax           = Lmax

      CALL alloc_NParray(L_TO_n_para%tab_L_TO_n,(/Lmax/),               &
                        "L_TO_n_para%tab_L_TO_n",name_sub,(/0/) )

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        CALL Write_Basis_L_TO_n(L_TO_n_para)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE alloc_Basis_L_TO_n

      SUBROUTINE dealloc_Basis_L_TO_n(L_TO_n_para)

      TYPE (Basis_L_TO_n), intent(inout) :: L_TO_n_para


!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='dealloc_Basis_L_TO_n'
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------

      L_TO_n_para%L_TO_n_type = 0
      L_TO_n_para%Lmax        = -1
      L_TO_n_para%A           = 1
      L_TO_n_para%B           = 1
      L_TO_n_para%C           = 0
      L_TO_n_para%expo        = 1


      IF (allocated(L_TO_n_para%tab_L_TO_n)) THEN
        CALL dealloc_NParray(L_TO_n_para%tab_L_TO_n,                    &
                            "L_TO_n_para%tab_L_TO_n",name_sub)
      END IF

      IF (allocated(L_TO_n_para%tab_n_TO_L)) THEN
        CALL dealloc_NParray(L_TO_n_para%tab_n_TO_L,                    &
                            "L_TO_n_para%tab_n_TO_L",name_sub)
      END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        CALL Write_Basis_L_TO_n(L_TO_n_para)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE dealloc_Basis_L_TO_n

      SUBROUTINE Write_Basis_L_TO_n(L_TO_n_para,Rec_line)
      TYPE (Basis_L_TO_n), intent(in) :: L_TO_n_para
      character (len=*), optional     :: Rec_line

      integer :: L
!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Write_Basis_L_TO_n'
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------


      IF (present(Rec_line)) THEN
        write(out_unitp,*) trim(Rec_line),'L_TO_n_type     ',L_TO_n_para%L_TO_n_type
        write(out_unitp,*) trim(Rec_line),'Lmax            ',L_TO_n_para%Lmax
        write(out_unitp,*) trim(Rec_line),'A               ',L_TO_n_para%A
        write(out_unitp,*) trim(Rec_line),'B               ',L_TO_n_para%B
        write(out_unitp,*) trim(Rec_line),'C               ',L_TO_n_para%C
        write(out_unitp,*) trim(Rec_line),'expo            ',L_TO_n_para%expo
        write(out_unitp,*) trim(Rec_line),'max_n           ',L_TO_n_para%max_n

        write(out_unitp,*) trim(Rec_line),'alloc tab_L_TO_n',allocated(L_TO_n_para%tab_L_TO_n)

        IF (allocated(L_TO_n_para%tab_L_TO_n)) THEN
          write(out_unitp,*) trim(Rec_line),'tab_L_TO_n'
          write(out_unitp,*) trim(Rec_line),':  ',L_TO_n_para%tab_L_TO_n
        END IF

        write(out_unitp,*) trim(Rec_line),'alloc tab_n_TO_L',allocated(L_TO_n_para%tab_n_TO_L)

        IF (allocated(L_TO_n_para%tab_n_TO_L)) THEN
          write(out_unitp,*) trim(Rec_line),'tab_n_TO_L'
          write(out_unitp,*) trim(Rec_line),':  ',L_TO_n_para%tab_n_TO_L
        END IF

      ELSE
        write(out_unitp,*) 'L_TO_n_type     ',L_TO_n_para%L_TO_n_type
        write(out_unitp,*) 'Lmax            ',L_TO_n_para%Lmax
        write(out_unitp,*) 'A               ',L_TO_n_para%A
        write(out_unitp,*) 'B               ',L_TO_n_para%B
        write(out_unitp,*) 'C               ',L_TO_n_para%C
        write(out_unitp,*) 'expo            ',L_TO_n_para%expo
        write(out_unitp,*) 'max_n           ',L_TO_n_para%max_n
        write(out_unitp,*) 'alloc tab_L_TO_n',allocated(L_TO_n_para%tab_L_TO_n)

        IF (allocated(L_TO_n_para%tab_L_TO_n)) THEN
          write(out_unitp,*) 'tab_L_TO_n'
          write(out_unitp,*) '   ',L_TO_n_para%tab_L_TO_n
        END IF

        write(out_unitp,*) 'alloc tab_n_TO_L',allocated(L_TO_n_para%tab_n_TO_L)

        IF (allocated(L_TO_n_para%tab_n_TO_L)) THEN
          write(out_unitp,*) 'tab_n_TO_L'
          write(out_unitp,*) '   ',L_TO_n_para%tab_n_TO_L
        END IF

      END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE Write_Basis_L_TO_n

      SUBROUTINE init_Basis_L_TO_n(L_TO_n_para,Lmax)

      TYPE (Basis_L_TO_n), intent(inout) :: L_TO_n_para
      integer,             intent (in)   :: Lmax

      integer  :: L,errBasis_L_TO_n,nmin,nmax,n


!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='init_Basis_L_TO_n'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
        CALL Write_Basis_L_TO_n(L_TO_n_para)

      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------

    IF (L_TO_n_para%L_TO_n_type /=2) THEN
      IF (allocated(L_TO_n_para%tab_L_TO_n))                            &
           CALL dealloc_NParray(L_TO_n_para%tab_L_TO_n,"L_TO_n_para%tab_L_TO_n",name_sub)

      IF (Lmax < 0 ) THEN
        write(out_unitp,*) ' ERROR in ',name_sub
        write(out_unitp,*) ' Wrong Lmax value',Lmax
        STOP
      END IF

      L_TO_n_para%Lmax = Lmax
      CALL alloc_NParray(L_TO_n_para%tab_L_TO_n,(/Lmax/),               &
                        "L_TO_n_para%tab_L_TO_n",name_sub,(/0/) )
    END IF

      SELECT CASE (L_TO_n_para%L_TO_n_type)
      CASE (0) ! n(L) = A + B * L**expo
        DO L=0,Lmax
          L_TO_n_para%tab_L_TO_n(L) = L_TO_n_para%A + L_TO_n_para%B * L**L_TO_n_para%expo
        END DO

      CASE (1) ! Delta_n(L) = B * L**(expo-1)  => n(L) = A + B * L*(L+1)/2 (with expo=2)
        L_TO_n_para%tab_L_TO_n(0) = L_TO_n_para%A
        DO L=1,Lmax
          L_TO_n_para%tab_L_TO_n(L) = L_TO_n_para%tab_L_TO_n(L-1) + L_TO_n_para%B * L**(L_TO_n_para%expo-1)
        END DO

      CASE (2) ! The tab_L_TO_n is already present
        CONTINUE !

      CASE (3) ! Delta_n(L) = B * L**(expo-1)  => n(L) = A + B * L*(L-1)/2 (with expo=2)
        L_TO_n_para%tab_L_TO_n(0) = L_TO_n_para%A
        DO L=1,Lmax
          L_TO_n_para%tab_L_TO_n(L) = L_TO_n_para%tab_L_TO_n(L-1) + L_TO_n_para%B * int(ONE+log(real(L,kind=Rkind)))
        END DO

      CASE (4) ! n(L) = A+ B * tanh(L/4)
        DO L=0,Lmax
          L_TO_n_para%tab_L_TO_n(L) = L_TO_n_para%A +       &
            int( real(L_TO_n_para%B,kind=Rkind) * tanh(real(L,kind=Rkind)/FOUR) )
        END DO

      CASE Default
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) '  WRONG L_TO_n_type',L_TO_n_para%L_TO_n_type
        STOP
      END SELECT

      WHERE (L_TO_n_para%tab_L_TO_n > L_TO_n_para%max_n)
        L_TO_n_para%tab_L_TO_n = L_TO_n_para%max_n
      END WHERE

      CALL check_Basis_L_TO_n(L_TO_n_para%tab_L_TO_n,errBasis_L_TO_n)

      IF (errBasis_L_TO_n /= 0) THEN
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) '  Problem with initialization'
        CALL Write_Basis_L_TO_n(L_TO_n_para)
        write(out_unitp,*)
        STOP
      END IF

      nmax = get_n_FROM_Basis_L_TO_n(L_TO_n_para,Lmax)
      nmin = get_L_FROM_Basis_L_TO_n(L_TO_n_para,nmax) ! to set up the table

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        CALL Write_Basis_L_TO_n(L_TO_n_para)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE init_Basis_L_TO_n

      FUNCTION Get_n_FROM_Basis_L_TO_n(L_TO_n_para,L,L2) RESULT(n)
      integer  :: n
      TYPE (Basis_L_TO_n), intent(in)    :: L_TO_n_para
      integer,             intent (in)   :: L
      integer, optional,   intent (in)   :: L2

      integer :: Ltmp,u,L2_loc,n2B,n2C

!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Get_n_FROM_Basis_L_TO_n'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'L',L
        write(out_unitp,*) 'L2?',present(L2)
        IF (present(L2)) write(out_unitp,*) 'L2',L2
        CALL Write_Basis_L_TO_n(L_TO_n_para)
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------



    IF (L < 0) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) 'L',L
      write(out_unitp,*) 'L2?',present(L2)
      IF (present(L2)) write(out_unitp,*) 'L2',L2
      CALL Write_Basis_L_TO_n(L_TO_n_para)
      write(out_unitp,*) '  L < 0',L
      STOP
    END IF

    L2_loc = L
    IF (present(L2)) L2_loc = L2

    IF (L2_loc < 0) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) 'L',L
      write(out_unitp,*) 'L2?',present(L2)
      IF (present(L2)) write(out_unitp,*) 'L2',L2
      CALL Write_Basis_L_TO_n(L_TO_n_para)
      write(out_unitp,*) '  L2 < 0',L2_loc
      STOP
    END IF

    IF (allocated(L_TO_n_para%tab_L_TO_n)) THEN
      u = ubound(L_TO_n_para%Tab_L_TO_n,dim=1)
      IF (L <= u) THEN
        n = L_TO_n_para%tab_L_TO_n(L)
      ELSE
        n = L_TO_n_para%Tab_L_TO_n(u) +                                 &
           (l-u)*(L_TO_n_para%Tab_L_TO_n(u)-L_TO_n_para%Tab_L_TO_n(u-1))
      END IF
    ELSE

      SELECT CASE (L_TO_n_para%L_TO_n_type)
      CASE (0) ! A + B * L**expo
        IF (L <= L2_loc) THEN
          n   = L_TO_n_para%A + L_TO_n_para%B * L     **L_TO_n_para%expo
        ELSE
          n2B = L_TO_n_para%A + L_TO_n_para%B * L2_loc**L_TO_n_para%expo
          n2C = L_TO_n_para%A + L_TO_n_para%C * L2_loc**L_TO_n_para%expo
          n   = L_TO_n_para%A + (n2B-n2C)     + L_TO_n_para%C * L     **L_TO_n_para%expo
        END IF

      CASE (1) ! Delta_n(L) = B * L**(expo-1)  => n(L) = A + B * L*(L-1)/2 (with expo=2)
        n = L_TO_n_para%A
        DO Ltmp=1,L
          n = n + L_TO_n_para%B * Ltmp**(L_TO_n_para%expo-1)
        END DO

      CASE (2) ! read tab_L_TO_n
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) '   L_TO_n_type=2 and tab_L_TO_n is not allocated !'
        STOP

      CASE (3)
        n = L_TO_n_para%A
        DO Ltmp=1,L
          n = n + L_TO_n_para%B * int(ONE+log(real(Ltmp,kind=Rkind)))
        END DO

      CASE (4)
        n = L_TO_n_para%A +       &
          int( real(L_TO_n_para%B,kind=Rkind) * tanh(real(L,kind=Rkind)/FOUR) )

      CASE Default
        write(out_unitp,*) 'ERROR in ',name_sub
        write(out_unitp,*) '  WRONG L_TO_n_type',L_TO_n_para%L_TO_n_type
        STOP
      END SELECT

    END IF

    n = min(n,L_TO_n_para%max_n)

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'L,max_n,n',L,L_TO_n_para%max_n,n
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------
      END FUNCTION Get_n_FROM_Basis_L_TO_n

      FUNCTION Get_L_FROM_Basis_L_TO_n(L_TO_n_para,n) RESULT(L)
      integer  :: L
      TYPE (Basis_L_TO_n), intent(inout)    :: L_TO_n_para
      integer,             intent (in)   :: n


      integer :: nni,nn0,nn,nmin,nmax

!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Get_L_FROM_Basis_L_TO_n'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'n',n
        CALL Write_Basis_L_TO_n(L_TO_n_para)
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------

    nmax = get_n_FROM_Basis_L_TO_n(L_TO_n_para,L_TO_n_para%Lmax)
    nmin = get_n_FROM_Basis_L_TO_n(L_TO_n_para,0)

    IF (n < nmin .OR. n > nmax) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      CALL Write_Basis_L_TO_n(L_TO_n_para)
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) 'n is out of range. n:',n
      write(out_unitp,*) '   range:',nmin,nmax
      STOP
    END IF

    IF (.NOT. allocated(L_TO_n_para%tab_n_TO_L)) THEN
      CALL alloc_NParray(L_TO_n_para%tab_n_TO_L,(/nmax/),'tab_n_TO_L',name_sub,(/nmin/))

      DO L=L_TO_n_para%Lmax,0,-1
        nn = get_n_FROM_Basis_L_TO_n(L_TO_n_para,L)
        L_TO_n_para%tab_n_TO_L(nmin:nn) = L
      END DO
    END IF

    L = L_TO_n_para%tab_n_TO_L(n)



!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'L',L
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------
      END FUNCTION Get_L_FROM_Basis_L_TO_n

      SUBROUTINE Set_Basis_L_TO_n(L_TO_n_para,A,B,C,expo,Tab_L_TO_n,max_n,L_TO_n_type)

      TYPE (Basis_L_TO_n), intent(inout) :: L_TO_n_para
      integer, optional,   intent (in)   :: A,B,C,expo,L_TO_n_type,max_n
      integer, optional,   intent (in), allocatable   :: Tab_L_TO_n(:)


      integer  :: L,Aloc,Bloc,Cloc,expoloc,L_TO_n_type_loc,max_n_loc
      integer  :: errBasis_L_TO_n


!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='Set_Basis_L_TO_n'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------

      Aloc = 1
      IF (present(A)) Aloc=A

      Bloc = 1
      IF (present(B)) Bloc=B

      Cloc = 0
      IF (present(C)) Cloc=C

      expoloc = 1
      IF (present(expo)) expoloc=expo

      L_TO_n_type_loc = 0
      IF (present(L_TO_n_type)) L_TO_n_type_loc=L_TO_n_type

      max_n_loc = huge(1)
      IF (present(max_n)) max_n_loc=max_n

      CALL dealloc_Basis_L_TO_n(L_TO_n_para)

      L_TO_n_para%A           = Aloc
      L_TO_n_para%B           = Bloc
      L_TO_n_para%C           = Cloc
      L_TO_n_para%expo        = expoloc
      L_TO_n_para%max_n       = max_n_loc
      L_TO_n_para%L_TO_n_type = L_TO_n_type_loc

      IF (present(Tab_L_TO_n)) THEN
        IF (.NOT. allocated(Tab_L_TO_n)) STOP ' in Set_Basis_L_TO_n: Tab_L_TO_n is not allocated'
        L_TO_n_para%L_TO_n_type = 2
        L_TO_n_para%Tab_L_TO_n  = Tab_L_TO_n
      END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        CALL Write_Basis_L_TO_n(L_TO_n_para)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE Set_Basis_L_TO_n

      SUBROUTINE check_Basis_L_TO_n(tab_L_TO_n,errBasis_L_TO_n)

      integer, allocatable,     intent (in)    :: tab_L_TO_n(:)
      integer,                  intent (inout) :: errBasis_L_TO_n


      integer :: L
!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='check_Basis_L_TO_n'
      logical,parameter :: debug=.FALSE.
!      logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*)
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------
      errBasis_L_TO_n = 0

      IF (.NOT. allocated(tab_L_TO_n)) errBasis_L_TO_n = 1
      IF (errBasis_L_TO_n /= 0) RETURN

      IF (tab_L_TO_n( lbound(tab_L_TO_n,dim=1) ) < 0)  errBasis_L_TO_n = -1 ! tab_L_TO_n(0) < 1
      IF (errBasis_L_TO_n /= 0) RETURN


      DO L=lbound(tab_L_TO_n,dim=1)+1,ubound(tab_L_TO_n,dim=1) ! tab_L_TO_n(L) < tab_L_TO_n(L-1)
        IF ( tab_L_TO_n( L) < tab_L_TO_n( L-1) ) THEN
          errBasis_L_TO_n = -1
          RETURN
        END IF
      END DO



!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------


!---------------------------------------------------------------------

      END SUBROUTINE check_Basis_L_TO_n

      SUBROUTINE L_TO_n_para2_TO_L_TO_n_para1(L_TO_n_para1,L_TO_n_para2)
      TYPE (Basis_L_TO_n), intent(inout) :: L_TO_n_para1
      TYPE (Basis_L_TO_n), intent(in)    :: L_TO_n_para2


!---------------------------------------------------------------------
      integer :: err_mem,memory
      character (len=*), parameter :: name_sub='L_TO_n_para2_TO_L_TO_n_para1'
      logical,parameter :: debug=.FALSE.
      !logical,parameter :: debug=.TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'BEGINNING ',name_sub
        write(out_unitp,*) 'L_TO_n_para2'
        CALL Write_Basis_L_TO_n(L_TO_n_para2)
        CALL flush_perso(out_unitp)
      END IF
!---------------------------------------------------------------------
!---------------------------------------------------------------------


      L_TO_n_para1%L_TO_n_type = L_TO_n_para2%L_TO_n_type
      L_TO_n_para1%Lmax        = L_TO_n_para2%Lmax
      L_TO_n_para1%A           = L_TO_n_para2%A
      L_TO_n_para1%B           = L_TO_n_para2%B
      L_TO_n_para1%C           = L_TO_n_para2%C
      L_TO_n_para1%expo        = L_TO_n_para2%expo
      L_TO_n_para1%max_n       = L_TO_n_para2%max_n

      IF (allocated(L_TO_n_para2%tab_L_TO_n)) L_TO_n_para1%tab_L_TO_n  = L_TO_n_para2%tab_L_TO_n
      IF (allocated(L_TO_n_para2%tab_n_TO_L)) L_TO_n_para1%tab_n_TO_L  = L_TO_n_para2%tab_n_TO_L

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*)
        write(out_unitp,*) 'L_TO_n_para1'
        CALL Write_Basis_L_TO_n(L_TO_n_para1)
        write(out_unitp,*)
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------

      END SUBROUTINE L_TO_n_para2_TO_L_TO_n_para1

      END MODULE mod_Basis_L_TO_n
