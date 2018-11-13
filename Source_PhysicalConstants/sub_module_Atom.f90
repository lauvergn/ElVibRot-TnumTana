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
!      with contributions of Mamadou Ndong
!
!===========================================================================
!===========================================================================

!> @brief Module which enables to use isotopic masses.
!!
!> @author David Lauvergnat
!! @date 16/12/2016
!!
!> @brief This module has two options:
!! @li Masses from Handbook of Chemistry and Physics 70th edition (B-228) with the
!! "construct_old_table_at" subroutine.
!! @li Masses from <a href="https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses">NIST</a>.
!!   with the "construct_table_at" subroutine.
!!  This subroutine uses an internal file "internal_data/IsotopicMass.txt" download in 2012 from NIST.
      MODULE mod_Atom
      use mod_system, only: name_len, rkind, one, zero, param_file, int_to_char,   &
                            ten, out_unitp, evrt_path, file_open, error_memo_allo, &
                            alloc_array, file_close, print_level,                  &
                            string_uppercase_to_lowercase, dealloc_array
      IMPLICIT NONE

      PRIVATE

!> @brief Derived type in which the isotopic mass will be set-up.
!!
!> @author David Lauvergnat
!! @date 17/12/2016
!!
!! @param isotope        string: the atomic symbole with A (ex "12C")
!! @param symbol         string: the atomic symbole
!! @param Z              integer: number of electrons
!! @param A              integer: number of nucleons
!! @param mass           real: the mass
!! @param mass_unit      string: the mass unit
!! @param spin           real: the nucleus spin (in atomic unit, *h/2Pi)
!! @param average_mass   real: the average mass
!! @param abundance      real: the isotopic abundance
  TYPE atom
    character (len=Name_len) :: isotope      = ''      !<  the atomic symbole with A (ex "12C")
    character (len=Name_len) :: symbol       = ''      !<  the atomic symbole
    integer                  :: Z            = 0       !<  number of electrons
    integer                  :: A            = 0       !<  number of nucleons
    real (kind=Rkind)        :: mass         = -ONE    !< the mass
    character (len=Name_len) :: mass_unit    = 'g/mol' !< the mass unit
    real (kind=Rkind)        :: spin         = ZERO    !< the nucleus spin (in atomic unit, *h/2Pi)
    real (kind=Rkind)        :: average_mass = ZERO    !< the average mass
    real (kind=Rkind)        :: abundance    = ZERO    !< the isotopic abundance
  END TYPE atom

!> @brief Derived type in which all isotopes are set-up.
!!
!> @author David Lauvergnat
!! @date 17/12/2016
!!
!! @param at             table of derived type "atom". 0:max_Z,0:max_A)
!! @param Def_isotope    table of integer: to define the most abundant isotope
!! @param max_Z          integer: upper bound of the 1st dimension of at(:,:)
!! @param max_A          integer: upper bound of the 2d dimension of at(:,:)
!! @param construct      logical: flag to check if the tables have been build
!! @param List_Isotope   derived type file: file where the list of the isotopes are read.

  TYPE table_atom
    TYPE (atom), dimension(:,:), pointer            :: at          => null() !< table of derived type "atom"
    integer, dimension(:), pointer                  :: Def_isotope => null() !< table to define the most abundant isotope
    character (len=Name_len), dimension(:), pointer :: name_at_Z   => null() !< ????
    integer                                         :: max_A       = 0 !< upper bound of the 1st dimension of at(:,:)
    integer                                         :: max_Z       = 0 !< upper bound of the 2d dimension of at(:,:)
    logical                                         :: construct   = .FALSE. !< flag to check if the tables have been build

    TYPE(param_file)                                :: List_Isotope !< file where the list of the isotopes are read


  END TYPE table_atom

  PUBLIC :: table_atom, construct_table_at, construct_old_table_at, dealloc_table_at
  PUBLIC :: get_mass_Tnum

  CONTAINS

!> @brief Subroutine to build the derived type "atom".
!!
!> @author David Lauvergnat
!! @date 17/12/2016
!!
!! @param at            derived type "atom" in which the isotopic mass will be set-up.
!! @param isotope       optional string: the atomic symbole
!! @param Z             optional integer: number of electrons
!! @param A             optional integer: number of nucleons
!! @param mass          optional real: the mass
!! @param mass_unit     optional string: the mass unit
!! @param spin          optional real: the nucleus spin
!! @param average_mass  optional real: the average mass
!! @param abundance     optional real: the isotopic abundance

  SUBROUTINE construct_atom(at,isotope,Z,A,mass,mass_unit,spin,   &
                            average_mass,abundance)

    TYPE (atom), intent(inout)                   :: at

    character (len=*) , intent(in), optional     :: isotope
    integer           , intent(in), optional     :: Z,A
    real (kind=Rkind)     , intent(in), optional :: mass
    character (len=*), intent(in), optional      :: mass_unit
    real (kind=Rkind)     , intent(in), optional :: spin
    real (kind=Rkind)     , intent(in), optional :: average_mass
    real (kind=Rkind)     , intent(in), optional :: abundance

    IF ( present(isotope) )      at%isotope      = isotope
    IF ( present(Z) )            at%Z            = Z
    IF ( present(A) )            at%A            = A
    IF ( present(mass) )         at%mass         = mass
    IF ( present(mass_unit) )    at%mass_unit    = mass_unit
    IF ( present(spin) )         at%spin         = spin
    IF ( present(average_mass) ) at%average_mass = average_mass
    IF ( present(abundance) )    at%abundance    = abundance

  END SUBROUTINE construct_atom

!> @brief Subroutine which reads isotpic masses from the
!!        an internal file "internal_data/IsotopicMass.txt" download in 2012 from NIST.
!!
!> @author David Lauvergnat
!! @date 16/12/2016
!!
!! @param at         derived type "atom" in which the isotopic mass will be set-up.
!! @param nio        file unit of the internal file.
!! @param err        integer to handle err

  SUBROUTINE Read_atom(at,nio,err)

    TYPE (atom), intent(inout)    :: at
    integer, intent(in)           :: nio
    integer, intent(inout)        :: err

    err = 0
    read(nio,*,IOSTAT=err) at%Z,at%A,at%symbol,at%mass,at%abundance
    IF (err /= 0) RETURN

    at%isotope = int_TO_char(at%A) // trim(adjustl(at%symbol))

    IF (at%Z == 1 .AND. at%A == 1) at%isotope='H'
    IF (at%Z == 1 .AND. at%A == 2) at%isotope='D'
    IF (at%Z == 1 .AND. at%A == 3) at%isotope='T'

    IF (abs(real(at%A,kind=Rkind)/at%mass-ONE) > TEN**(-2)) THEN
      write(out_unitp,*) '  WARNNING in Read_atom'
      write(out_unitp,*) '  The atomic mass in g/mol is probably too different from A (> 1%)'
      write(out_unitp,*) '  mass, A: ',at%mass,at%A
      write(out_unitp,*) '  atom: ',at
      !STOP
    END IF

    !write(out_unitp,*) 'at: ',at

  END SUBROUTINE Read_atom

!> @brief Subroutine which initializes the isotpic masses from
!!        from <a href="https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses">NIST</a>.
!!
!> @author David Lauvergnat
!! @date 16/12/2016
!!
!! @param mendeleev     table of "atoms" with known isotopes
!! @param auTOgPmol     conversion factor: au to g.mol-1
!! @param mass_unit     name of the mass unit (eg gPmol)

  SUBROUTINE construct_table_at(mendeleev,auTOgPmol,mass_unit)

    TYPE (table_atom) :: mendeleev
    character (len=*) :: mass_unit
    real (kind=Rkind) :: auTOgPmol

    integer           :: nio,err
    integer           :: memory
    TYPE (atom)       :: at
    integer           :: Z,A,err_mem
    real (kind=Rkind) :: Max_abundance

    character (len=*), parameter ::                               &
       alphab1="ABCBEFGHIJKLMNOPQRSTUVWXYZ"
    character (len=*), parameter ::                               &
       alphab2="abcdefghijklmnopqrstuvwxyz"

    integer :: i,pos1,pos2,pos3

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='construct_table_at'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
      END IF
!---------------------------------------------------------------------

    mendeleev%List_Isotope%name = trim(EVRT_path) // '/' //       &
                                  'Internal_data/IsotopicMass.txt'

    IF (.NOT. mendeleev%construct) THEN

      CALL file_open(mendeleev%List_Isotope,nio,old=.TRUE.,err_file=err)
      IF (err /= 0) THEN
        write(out_unitp,*) ' WARNNING in construct_table_at'
        write(out_unitp,*) ' The file "Internal_data/IsotopicMass.txt" can not be open'
        write(out_unitp,*) ' Two main raisons:'
        write(out_unitp,*) '  -1- the ElVibRot directory has been moved after the compilation'
        write(out_unitp,*) '  -2- the "EVRT_path" is badly defined in the namelist "system"'
        write(out_unitp,*) '      EVRT_path: ',trim(EVRT_path)

        write(out_unitp,*) ' Two solutions:'
        write(out_unitp,*) '  -1- recompile ElVibRot: "make clean ; make"'
        write(out_unitp,*) '  -2- set-up the "EVRT_path" in the namelist "system"'

        write(out_unitp,*) ' WARNNING the old construct_table_at is used'
        RETURN
      END IF

      mendeleev%construct = .TRUE.

      mendeleev%max_Z = 110
      mendeleev%max_A = 300
      memory = product( (/ 1+mendeleev%max_Z,1+mendeleev%max_A /) &
                                                                 )
      allocate(mendeleev%at(0:mendeleev%max_Z,0:mendeleev%max_A), &
                                                     stat=err_mem) ! change alloc done
      CALL error_memo_allo(err_mem,memory,"mendeleev%at",         &
                                      "construct_table_at","atom")

      CALL alloc_array(mendeleev%Def_isotope,(/mendeleev%max_Z/), &
                      "mendeleev%Def_isotope","construct_table_at",(/0/))
      mendeleev%Def_isotope(:) = 0

      CALL alloc_array(mendeleev%name_at_Z,(/mendeleev%max_Z+3/),Name_len, &
                      "mendeleev%name_at_Z","construct_table_at",(/0/))
      mendeleev%name_at_Z(:) = ' '


      DO Z=0,mendeleev%max_Z
      DO A=0,mendeleev%max_A
        mendeleev%at(Z,A) =                                       &
          atom(' ',' ',Z,A,-ONE,'g/mol',ZERO,ZERO,ZERO)
      END DO
      END DO

      read(nio,*)   ! first line
      DO
        CALL Read_atom(at,nio,err)
        IF (err /= 0) EXIT
        Z = at%Z
        A = at%A
        IF (Z < lbound(mendeleev%at,dim=1) .OR. Z > ubound(mendeleev%at,dim=1)) CYCLE
        IF (A < lbound(mendeleev%at,dim=2) .OR. A > ubound(mendeleev%at,dim=2)) CYCLE

        mendeleev%at(Z,A) = at
        !IF (debug) write(out_unitp,*) 'Z,A,at',Z,A,mendeleev%at(Z,A)
        IF (debug) write(out_unitp,*) 'at: ',mendeleev%at(Z,A)

      END DO
      CALL file_close(mendeleev%List_Isotope)


      DO Z=0,mendeleev%max_Z
        mendeleev%Def_isotope(Z) = 0
        Max_abundance = ZERO
        DO A=0,mendeleev%max_A
          mendeleev%at(Z,A)%mass = mendeleev%at(Z,A)%mass/auTOgPmol
          mendeleev%at(Z,A)%mass_unit = mass_unit
          IF (mendeleev%at(Z,A)%abundance > Max_abundance) THEN
            Max_abundance            = mendeleev%at(Z,A)%abundance
            mendeleev%Def_isotope(Z) = A
            !mendeleev%name_at_Z(Z)   = mendeleev%at(Z,A)%symbol
          END IF
        END DO
        IF (debug) write(out_unitp,*) 'Z,Def_isotope,symbol',Z,mendeleev%Def_isotope(Z),mendeleev%name_at_Z(Z)
      END DO


      DO Z=0,mendeleev%max_Z+3
          mendeleev%name_at_Z(Z) = ' '
      END DO
      i = 0
      !- atom with 2 letters ------------------
      DO Z=0,mendeleev%max_Z
        A = mendeleev%Def_isotope(Z)

        !write(6,*) 'Z,A,mendeleev%at(Z,A)%isotope ',Z,A, " ",mendeleev%at(Z,A)%isotope
        pos1=scan(mendeleev%at(Z,A)%isotope,alphab2)
        pos2=len(mendeleev%at(Z,A)%isotope)
        IF (pos1 > 1) THEN
          i = i+1
          mendeleev%name_at_Z(i) =                                &
          trim(mendeleev%at(Z,A)%isotope(pos1-1:pos2)) // '   ' // int_TO_char(Z)
        END IF
      END DO
      !- atom with 1 letter ------------------
      DO Z=0,mendeleev%max_Z
        A = mendeleev%Def_isotope(Z)
        !write(6,*) 'Z,A,mendeleev%at(Z,A)%isotope ',Z,A, " ",mendeleev%at(Z,A)%isotope

        pos3=scan(mendeleev%at(Z,A)%isotope,alphab2)
        IF (pos3 == 0) THEN
          pos1=scan(mendeleev%at(Z,A)%isotope,alphab1)
          pos2=len(mendeleev%at(Z,A)%isotope)
          IF (pos1 > 0) THEN
            i = i+1
            mendeleev%name_at_Z(i) =                              &
         trim(mendeleev%at(Z,A)%isotope(pos1:pos2)) // '   ' // int_TO_char(Z)
          END IF
        END IF
      END DO
      i = i+1
      mendeleev%name_at_Z(i) = "D   2"
      i = i+1
      mendeleev%name_at_Z(i) = "T   3"



      IF (debug .OR. print_level > 1) THEN
        write(out_unitp,*) '  Atomic list:'
        DO Z=0,mendeleev%max_Z
          IF (len_trim(mendeleev%name_at_Z(Z)) > 0) THEN
            IF (Z /= 1) write(out_unitp,*) Z,mendeleev%name_at_Z(Z)
            IF (Z == 1) write(out_unitp,*) Z,'H'

          END IF
        END DO
        write(out_unitp,*) '1 D'
        write(out_unitp,*) '1 T'

        write(out_unitp,*) '  END Atomic list:'
      END IF

    END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------

  END SUBROUTINE construct_table_at
!> @brief Subroutine which initializes the isotpic masses from
!!        the Handbook of Chemistry and Physics 70th edition (B-228).
!!
!> @author David Lauvergnat
!! @date 16/12/2016
!!
!! @param mendeleev     table of "atoms" with known isotopes
!! @param auTOgPmol     conversion factor: au to g.mol-1
!! @param mass_unit     name of the mass unit (eg gPmol)
  SUBROUTINE construct_old_table_at(mendeleev,auTOgPmol,mass_unit)


    TYPE (table_atom) :: mendeleev
    character (len=*) :: mass_unit
    real (kind=Rkind)     :: auTOgPmol

    character (len=*), parameter ::                               &
       alphab1="ABCBEFGHIJKLMNOPQRSTUVWXYZ"
    character (len=*), parameter ::                               &
       alphab2="abcdefghijklmnopqrstuvwxyz"

    integer :: i,Z,A,pos1,pos2,pos3
    integer :: err_mem,memory

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='construct_old_table_at'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
      END IF
!---------------------------------------------------------------------

    IF (.NOT. mendeleev%construct) THEN
      mendeleev%construct = .TRUE.

      mendeleev%max_Z = 110
      mendeleev%max_A = 300
      memory = product( (/ 1+mendeleev%max_Z,1+mendeleev%max_A /) &
                                                                 )
      allocate(mendeleev%at(0:mendeleev%max_Z,0:mendeleev%max_A), &
                                                     stat=err_mem) ! change alloc done
      CALL error_memo_allo(err_mem,memory,"mendeleev%at",         &
                                      "construct_old_table_at","atom")
      CALL alloc_array(mendeleev%Def_isotope,(/mendeleev%max_Z/), &
                      "mendeleev%Def_isotope","construct_old_table_at",(/0/))
      mendeleev%Def_isotope(:) = 0

      CALL alloc_array(mendeleev%name_at_Z,(/mendeleev%max_Z+3/),Name_len, &
                      "mendeleev%name_at_Z","construct_old_table_at")


      DO Z=0,mendeleev%max_Z
      DO A=0,mendeleev%max_A
        mendeleev%at(Z,A) =                                       &
                    atom(' ',' ',Z,A,-ONE,'g/mol',ZERO,ZERO,ZERO)
      END DO
      END DO

      CALL construct_atom(                                        &
         mendeleev%at(  0,  0),isotope='X',mass=ZERO)

      CALL construct_atom(                                        &
         mendeleev%at(  1,  1),isotope='H',mass=1.007825_Rkind)
      CALL construct_atom(                                        &
         mendeleev%at(  1,  2),isotope='D',mass=2.014000_Rkind)
      CALL construct_atom(                                        &
         mendeleev%at(  1,  3),isotope='T',mass=3.016050_Rkind)

      CALL construct_atom(                                        &
         mendeleev%at(  2,  3),isotope='3He',mass=3.016030_Rkind)
      CALL construct_atom(                                        &
         mendeleev%at(  2,  4),isotope='4He',mass=4.002600_Rkind)
      CALL construct_atom(                                        &
         mendeleev%at(  2,  5),isotope='5He',mass=5.012220_Rkind)
      CALL construct_atom(                                        &
         mendeleev%at(  2,  6),isotope='6He',mass=6.018886_Rkind)
      CALL construct_atom(                                        &
         mendeleev%at(  2,  8),isotope='8He',mass=8.033920_Rkind)

      mendeleev%Def_isotope(0:2) = (/0,1,4/) ! most abundant isotopes

      CALL construct_atom(                                        &
         mendeleev%at(  3,  5),isotope='5Li',mass=5.012540_Rkind)
      CALL construct_atom(                                        &
         mendeleev%at(  3,  6),isotope='6Li',mass=6.015210_Rkind)
      CALL construct_atom(                                        &
         mendeleev%at(  3,  7),isotope='7Li',mass=7.016003_Rkind)
      CALL construct_atom(                                        &
         mendeleev%at(  3,  8),isotope='8Li',mass=8.022485_Rkind)
      CALL construct_atom(                                        &
         mendeleev%at(  3,  9),isotope='9Li',mass=9.026789_Rkind)
      CALL construct_atom(                                        &
         mendeleev%at(  4,  9),isotope='9Be',mass=9.012182_Rkind)
      CALL construct_atom(                                        &
         mendeleev%at(  5, 10),isotope='10B',mass=10.012937_Rkind)
      CALL construct_atom(                                        &
         mendeleev%at(  5, 11),isotope='11B',mass=11.009305_Rkind)
      CALL construct_atom(                                        &
         mendeleev%at(  6, 12),isotope='12C',mass=12.000000_Rkind)
      CALL construct_atom(                                        &
         mendeleev%at(  7, 14),isotope='14N',mass=14.003074_Rkind)
      CALL construct_atom(                                        &
         mendeleev%at(  8, 16),isotope='16O',mass=15.994915_Rkind)
      CALL construct_atom(                                        &
         mendeleev%at(  9, 19),isotope='19F',mass=18.998403_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 10, 20),isotope='20Ne',mass=19.992435_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 10, 22),isotope='22Ne',mass=21.991383_Rkind)

      mendeleev%Def_isotope(3:10) = (/7,9,11,12,14,16,19,20/) ! most abundant isotopes

      CALL construct_atom(                                        &
        mendeleev%at( 11, 23),isotope='23Na',mass=22.989767_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 12, 24),isotope='24Mg',mass=23.985042_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 12, 25),isotope='25Mg',mass=24.985837_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 12, 26),isotope='26Mg',mass=25.982593_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 13, 27),isotope='27Al',mass=26.981540_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 14, 28),isotope='28Si',mass=27.976927_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 14, 29),isotope='29Si',mass=28.976495_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 14, 30),isotope='30Si',mass=29.97370_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 15, 31),isotope='31P',mass=30.973762_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 16, 32),isotope='32S',mass=31.972070_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 17, 35),isotope='35Cl',mass=34.968852_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 17, 37),isotope='37Cl',mass=36.965903_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 18, 40),isotope='40Ar',mass=39.962384_Rkind)

      mendeleev%Def_isotope(11:18) = (/23,24,27,28,31,32,35,40/) ! most abundant isotopes

      CALL construct_atom(                                        &
        mendeleev%at( 19, 39),isotope='39K',mass=38.963707_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 20, 40),isotope='40Ca',mass=39.962591_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 21, 45),isotope='45Sc',mass=44.955910_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 22, 48),isotope='48Ti',mass=47.947947_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 23, 51),isotope='51V',mass=50.943962_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 24, 52),isotope='52Cr',mass=51.940509_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 24, 52),isotope='53Cr',mass=52.940651_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 25, 55),isotope='55Mn',mass=54.938047_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 26, 56),isotope='56Fe',mass=55.934939_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 27, 59),isotope='59Co',mass=58.933198_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 28, 58),isotope='58Ni',mass=57.935346_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 28, 60),isotope='60Ni',mass=59.930788_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 29, 63),isotope='63Cu',mass=62.939598_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 29, 65),isotope='65Cu',mass=64.927793_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 30, 64),isotope='64Zn',mass=63.929145_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 30, 66),isotope='66Zn',mass=65.926034_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 30, 68),isotope='68Zn',mass=67.924846_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 31, 69),isotope='69Ga',mass=68.925580_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 31, 71),isotope='71Ga',mass=70.924700_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 32, 70),isotope='70Ge',mass=69.924250_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 32, 72),isotope='72Ge',mass=71.922079_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 32, 74),isotope='74Ge',mass=73.921177_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 33, 75),isotope='75As',mass=74.921594_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 34, 78),isotope='78Se',mass=78.000000_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 34, 80),isotope='80Se',mass=79.916520_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 35, 79),isotope='79Br',mass=78.918336_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 35, 81),isotope='81Br',mass=80.916289_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 36, 82),isotope='82Kr',mass=81.913482_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 36, 83),isotope='83Kr',mass=82.914135_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 36, 84),isotope='84Kr',mass=83.911507_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 36, 86),isotope='86Kr',mass=85.910616_Rkind)

      mendeleev%Def_isotope(19:36) = (/39,40,45,48,51,52,55,56,   &
                                   59,58,63,64,69,74,75,80,79,84/) ! most abundant isotopes

      CALL construct_atom(                                        &
        mendeleev%at( 78, 194),isotope='194Pt',mass=193.962655_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 78, 195),isotope='195Pt',mass=194.964766_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 78, 196),isotope='196Pt',mass=195.964926_Rkind)
      mendeleev%Def_isotope(78) = 195 ! most abundant isotopes



      DO Z=1,mendeleev%max_Z+3
          mendeleev%name_at_Z(Z) = ' '
      END DO
      i = 0
!     - atom with 2 letters ------------------
      DO Z=0,mendeleev%max_Z
        A = mendeleev%Def_isotope(Z)

        !write(6,*) 'Z,A,mendeleev%at(Z,A)%isotope ',Z,A, " ",mendeleev%at(Z,A)%isotope
        pos1=scan(mendeleev%at(Z,A)%isotope,alphab2)
        pos2=len(mendeleev%at(Z,A)%isotope)
        IF (pos1 > 1) THEN
          i = i+1
          mendeleev%name_at_Z(i) =                                &
          trim(mendeleev%at(Z,A)%isotope(pos1-1:pos2)) // '   ' // int_TO_char(Z)
        END IF
      END DO
!     - atom with 1 letter ------------------
      DO Z=0,mendeleev%max_Z
        A = mendeleev%Def_isotope(Z)

        pos3=scan(mendeleev%at(Z,A)%isotope,alphab2)
        IF (pos3 == 0) THEN
          pos1=scan(mendeleev%at(Z,A)%isotope,alphab1)
          pos2=len(mendeleev%at(Z,A)%isotope)
          IF (pos1 > 0) THEN
            i = i+1
            mendeleev%name_at_Z(i) =                              &
         trim(mendeleev%at(Z,A)%isotope(pos1:pos2)) // '   ' // int_TO_char(Z)
          END IF
        END IF
      END DO
!      i = i+1
!      mendeleev%name_at_Z(i) = "D 2"
!      i = i+1
!      mendeleev%name_at_Z(i) = "T 3"

      IF (debug .OR. print_level > 1) THEN
        write(out_unitp,*) '  Atomic list:'
        DO i=1,mendeleev%max_Z+3
          IF (len_trim(mendeleev%name_at_Z(i)) > 0) THEN
            write(out_unitp,*) i,mendeleev%name_at_Z(i)
          END IF
        END DO
        write(out_unitp,*) '  END Atomic list:'
      END IF

      DO Z=0,mendeleev%max_Z
      DO A=0,mendeleev%max_A
        mendeleev%at(Z,A)%mass = mendeleev%at(Z,A)%mass/auTOgPmol
        mendeleev%at(Z,A)%mass_unit = mass_unit
      END DO
      END DO


    END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------

  END SUBROUTINE construct_old_table_at
!> @brief Function: enables to get an isotopic mass
!!
!> @author David Lauvergnat
!! @date 16/12/2016
!! @param mendeleev     table of derived type "atom" with the known isotopes
!! @param Z             optional integer: number of electrons
!! @param A             optional integer: number of nucleons
!! @param name          optional string which contains the mass (a real) or the atomic symbol or "Z_A" or "A"symbol (12C)
!! @param err_mass      optional integer to handle error (err_mass=0 => no error)
!! @return              the isotopic mass in au (atomic unit)
!!
!> @li   At least "name" or "Z" must be given.
!! @li   When "A" is unknown, you get the most abundant isotope.
!! @li   When the atomic symbol is given (or "Z_A" or A"symbol") and when Z and A are present, then
!!       the function gives also the correponding Z and A values.
!! @li   Dummy atom is defined with the symbol "X" and/or Z=0.
!!
        FUNCTION get_mass_Tnum(mendeleev,Z,A,name,err_mass)
          real (kind=Rkind)                          :: get_mass_Tnum !< the isotopic mass in au (atomic unit)
          integer,           intent(inout), optional :: Z !< number of electrons
          integer,           intent(inout), optional :: A !< number of nucleons
          integer,           intent(out), optional   :: err_mass !< to handle error (err_mass=0 => no error)

          character (len=*), intent(in), optional    :: name !< string which contains the mass (a real), the atomic symbol, or "Z_A" or "A"symbol

          TYPE (table_atom),  intent(in)             :: mendeleev !< table of derived type "atom" with the known isotopes

          integer                         :: err_mass_loc
          integer                         :: AA,ZZ
          integer                         :: pos,i,clen
          real (kind=Rkind)               :: mass
          character (len=Name_len)        :: at_Z,name_Z,symb,name2
          character (len=:), allocatable  :: name_lower
          logical                         :: find

          err_mass_loc = 0
          get_mass_Tnum = -ONE
         !IF (present(A)) write(out_unitp,*) 'A present',A
         !IF (present(Z)) write(out_unitp,*) 'Z present',Z
         !IF (present(name)) write(out_unitp,*) 'name present',name

          AA = -1
          IF (present(A)) AA = A
          ZZ = -1
          IF (present(Z)) ZZ = Z


          IF ( AA > -1 .AND. ZZ >-1 ) THEN
            CONTINUE
          ELSE IF ( AA < 0 .AND. ZZ >-1 ) THEN
            AA = mendeleev%Def_isotope(Z)
          ELSE IF ( present(name) ) THEN
            ! name has the value of the mass....
            IF (index(name,".") .NE. 0) THEN
              read(name,*,IOSTAT=err_mass_loc) mass
              IF (err_mass_loc /= 0) THEN
                write(out_unitp,*) ' ERROR in : get_mass_Tnum'
                write(out_unitp,*) ' I CANNOT read the mass in "',trim(adjustl(name)),'"'
                IF (present(err_mass)) THEN
                  mass = -ONE
                  err_mass = err_mass_loc
                ELSE
                  STOP ' ERROR in : get_mass_Tnum'
                END IF
              END IF
              get_mass_Tnum = mass
              ZZ = 0
              AA = 0
              IF (present(Z)) Z = ZZ
              IF (present(A)) A = AA
              IF (print_level > 1)                                      &
                  write(out_unitp,*) 'get_mass_Tnum (read): ',ZZ,AA,mass
              RETURN

            ELSE IF (index(name,"_") .NE. 0) THEN
              ! isotope defined as : Z_A
              pos = index(name,"_")
              read(name(1:pos-1),*) ZZ

              IF (len_trim(name) < pos+1) THEN
                AA = mendeleev%Def_isotope(ZZ)
              ELSE
                read(name(pos+1:len(name)),*) AA
              END IF
            ELSE

              find = .FALSE.

              clen = len_trim(adjustl(name))
              allocate(character(len=clen) :: name_lower)
              name_lower = trim(adjustl(name))

              CALL string_uppercase_TO_lowercase(name_lower)
              find = .FALSE.
              DO i=1,mendeleev%max_Z
                !- extract the symbol and Z of name_at_Z(i)
                at_Z = mendeleev%name_at_Z(i)
                CALL string_uppercase_TO_lowercase(at_Z)
                pos = index(at_Z," ")
                IF (pos < 2) exit
                symb   = at_Z(1:pos-1)
                name_Z = at_Z(pos+1:len(at_Z))

                read(name_Z,*) ZZ

                !- find the right atom in name
                pos = index(name_lower,trim(symb))
                IF (pos .NE. 0) THEN
                  name2 = name(pos:len(name_lower))
                  IF (len_trim(name2) == len_trim(symb)) THEN
                    find = .TRUE.
                    IF (pos-1 < 1) THEN
                      AA = mendeleev%Def_isotope(ZZ)
                    ELSE
                      read(name_lower(1:pos-1),*) AA
                    END IF
                  END IF
                END IF
                IF (find) exit
              END DO
              !----------------------------------------------
              !  the symbol can be D or T
              IF (.NOT. find) THEN
                IF (index(name_lower,"D") .NE. 0 .OR. index(name_lower,"d") .NE. 0) THEN
                  name2 = "D"
                  ZZ = 1
                  AA = 2
                  find = .TRUE.
                END IF
                IF (index(name_lower,"T") .NE. 0 .OR. index(name_lower,"t") .NE. 0) THEN
                  name2 = "T"
                  ZZ = 1
                  AA = 3
                  find = .TRUE.
                END IF
              END IF
              !----------------------------------------------
              IF (.NOT. find) THEN
                write(out_unitp,*) ' ERROR : get_mass_Tnum'
                write(out_unitp,*) ' I CANNOT get the right atom'
                write(out_unitp,*) ' Your atom is NOT in my list !!'
                write(out_unitp,*) ' Your atom :',name
                write(out_unitp,*) ' My list:'
                DO i=1,mendeleev%max_Z+3
                  IF (len_trim(mendeleev%name_at_Z(i)) > 0) THEN
                    write(out_unitp,*) i,mendeleev%name_at_Z(i)
                  END IF
                END DO
                IF (present(err_mass)) THEN
                  err_mass = -1
                ELSE
                  STOP ' ERROR in : get_mass_Tnum'
                END IF
                RETURN
              END IF

            END IF
          ELSE
            write(out_unitp,*) ' ERROR : get_mass_Tnum'
            write(out_unitp,*) ' Only A is present !!!'
            IF (present(err_mass)) THEN
                  err_mass = -1
            ELSE
                  STOP ' ERROR in : get_mass_Tnum'
            END IF
            RETURN
          END IF


          IF (mendeleev%at(ZZ,AA)%mass <= ZERO .AND. ZZ /=0) THEN
            write(out_unitp,*) ' ERROR in : get_mass_Tnum'
            write(out_unitp,*) 'the mass of this isotope is NOT defined !'
            write(out_unitp,*) 'Z,A:',ZZ,AA
            IF (present(err_mass)) THEN
                  err_mass = -1
            ELSE
                  STOP ' ERROR in : get_mass_Tnum'
            END IF
            RETURN
          END IF

          IF (present(Z)) Z = ZZ
          IF (present(A)) A = AA


!         write(out_unitp,*) 'atom: ',mendeleev%at(ZZ,AA)
          IF (print_level > 1)                                          &
             write(out_unitp,*) 'get_mass_Tnum: ',ZZ,AA,mendeleev%at(ZZ,AA)%mass

          get_mass_Tnum = mendeleev%at(ZZ,AA)%mass

        END FUNCTION get_mass_Tnum
!> @brief Subroutine which deallocates the derived type table_atom
!!
!> @author David Lauvergnat
!! @date 16/12/2016
!! @param mendeleev     table of "atoms" with known isotopes
        SUBROUTINE dealloc_table_at(mendeleev)


          TYPE (table_atom) :: mendeleev

          integer :: err_mem,memory

          memory = size(mendeleev%at)
          deallocate(mendeleev%at,stat=err_mem) ! change alloc done
          CALL error_memo_allo(err_mem,-memory,"mendeleev%at",          &
                                              "dealloc_table_at","atom")

          CALL dealloc_array(mendeleev%Def_isotope,                     &
                            "mendeleev%Def_isotope","dealloc_table_at")

          CALL dealloc_array(mendeleev%name_at_Z,                       &
                            "mendeleev%name_at_Z","dealloc_table_at")

        END SUBROUTINE dealloc_table_at

      END MODULE mod_Atom
