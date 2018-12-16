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
!! @date 29/11/2018
!!
!> @brief This module has two options:
!! @li Masses from Handbook of Chemistry and Physics 70th edition (B-228) with the
!! "construct_old_table_at" subroutine.
!! @li Masses from <a href="https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses">NIST</a>.
!!   with the "construct_table_at" subroutine.
!!  This subroutine uses an internal file "internal_data/IsotopicMass.txt" download in 2012 from NIST.
MODULE mod_Atom
  use mod_system
  IMPLICIT NONE

PRIVATE

!> @brief Derived type in which the isotopic mass will be set-up.
!!
!> @author David Lauvergnat
!! @date 29/11/2018
!!
!! @param isotope        string:  the atomic symbole with A (ex "12C")
!! @param symbol         string:  the atomic symbole
!! @param Z              integer: number of electrons
!! @param A              integer: number of nucleons
!! @param mass           real:    the mass
!! @param mass_unit      string:  the mass unit
!! @param spin           real:    the nucleus spin (in atomic unit, *h/2Pi)
!! @param average_mass   real:    the average mass
!! @param abundance      real:    the isotopic abundance
!! @param MainIsotope    logical: flag to define the most abundant isotope
!! @param SetIsotope     logical: flag to define when an isotope is set
  TYPE atom
    character (len=Name_len) :: isotope      = ''      !<  the atomic symbole with A (ex "12C")
    character (len=Name_len) :: symbol       = ''      !<  the atomic symbole
    integer                  :: Z            = 0       !<  number of electrons
    integer                  :: A            = 0       !<  number of nucleons
    real (kind=Rkind)        :: mass         = -ONE    !< the mass
    character (len=Name_len) :: mass_unit    = 'g/mol' !< the mass unit
    real (kind=Rkind)        :: spin         = ZERO    !< the nucleus spin (in atomic unit, *h/2Pi)
    real (kind=Rkind)        :: average_mass = ZERO    !< the average mass
    real (kind=Rkind)        :: abundance    = -ONE    !< the isotopic abundance
    logical                  :: MainIsotope  = .FALSE. !< flag to define the most abundant isotope
    logical                  :: SetIsotope   = .FALSE. !< flag to define when an isotope is set

  END TYPE atom

!> @brief Derived type in which all isotopes are set-up.
!!
!> @author David Lauvergnat
!! @date 29/11/2018
!!
!! @param at             table of derived type "atom". 0:max_Z,0:max_A)
!! @param max_Z          integer: upper bound of the 1st dimension of at(:,:)
!! @param max_A          integer: upper bound of the 2d dimension of at(:,:)
!! @param construct      logical: flag to check if the tables have been build
!! @param List_Isotope   derived type file: file where the list of the isotopes are read.
  TYPE table_atom
    TYPE (atom), dimension(:,:), pointer            :: at          => null() !< table of derived type "atom"
    integer                                         :: max_A       = 0 !< upper bound of the 1st dimension of at(:,:)
    integer                                         :: max_Z       = 0 !< upper bound of the 2d dimension of at(:,:)
    logical                                         :: construct   = .FALSE. !< flag to check if the tables have been build

    TYPE(param_file)                                :: List_Isotope !< file where the list of the isotopes are read

  END TYPE table_atom

  PUBLIC :: table_atom, dealloc_table_at
  PUBLIC :: construct_table_at_NIST2018, construct_table_at_NIST2012, construct_table_at_HandBook70ed
  PUBLIC :: get_mass_Tnum

  CONTAINS
!> @brief Subroutine to build the derived type "atom".
!!
!> @author David Lauvergnat
!! @date 29/11/2018
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
                            average_mass,abundance,MainIsotope)

    TYPE (atom), intent(inout)                   :: at

    character (len=*),      intent(in), optional :: isotope
    integer,                intent(in), optional :: Z,A
    real (kind=Rkind),      intent(in), optional :: mass
    character (len=*),      intent(in), optional :: mass_unit
    real (kind=Rkind),      intent(in), optional :: spin
    real (kind=Rkind),      intent(in), optional :: average_mass
    real (kind=Rkind),      intent(in), optional :: abundance
    logical,                intent(in), optional :: MainIsotope

    integer :: pos,i,err

    IF ( present(isotope) )      at%isotope      = isotope
    IF ( present(Z) )            at%Z            = Z
    IF ( present(A) )            at%A            = A
    IF ( present(mass) )         at%mass         = mass
    IF ( present(mass_unit) )    at%mass_unit    = mass_unit
    IF ( present(spin) )         at%spin         = spin
    IF ( present(average_mass) ) at%average_mass = average_mass
    IF ( present(abundance) )    at%abundance    = abundance
    IF ( present(MainIsotope) )  at%MainIsotope  = MainIsotope

    at%SetIsotope = .TRUE.

    IF ( present(isotope) ) THEN ! defined symbol from isotope
      pos = 0
      DO
        read(at%isotope(pos+1:pos+1),'(i1)',IOSTAT=err) i
        IF (err /= 0) EXIT
        pos = pos + 1
      END DO

      at%symbol = at%isotope(pos+1:len_trim(at%isotope))

    END IF


  END SUBROUTINE construct_atom

!> @brief Subroutine which reads isotpic masses from the
!!        an internal file "internal_data/IsotopicMass.txt" download in 2012 from NIST.
!!
!> @author David Lauvergnat
!! @date 29/11/2018
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

    at%isotope    = int_TO_char(at%A) // trim(adjustl(at%symbol))
    at%SetIsotope = .TRUE.

    IF (at%Z == 1 .AND. at%A == 1) at%isotope='H'
    IF (at%Z == 1 .AND. at%A == 2) at%isotope='D'
    IF (at%Z == 1 .AND. at%A == 3) at%isotope='T'

    IF (at%A > 0 .AND. abs(real(at%A,kind=Rkind)/at%mass-ONE) > TEN**(-2)) THEN
      write(out_unitp,*) '  WARNNING in Read_atom'
      write(out_unitp,*) '  The atomic mass in g/mol is probably too different from A (> 1%)'
      write(out_unitp,*) '  mass, A: ',at%mass,at%A
      write(out_unitp,*) '  atom: ',at
      !STOP
    END IF

    !write(out_unitp,*) 'at: ',at

  END SUBROUTINE Read_atom
  SUBROUTINE Read2018_atom(at,nio,err)

    TYPE (atom), intent(inout)    :: at
    integer, intent(in)           :: nio
    integer, intent(inout)        :: err

    err = 0
    read(nio,*,IOSTAT=err) at%Z,at%symbol,at%A,at%mass,at%abundance
    !write(out_unitp,*) err,'Z,symbol,A,mass,abundance',at%Z,at%symbol,at%A,at%mass,at%abundance
    IF (err /= 0) RETURN

    at%isotope    = int_TO_char(at%A) // trim(adjustl(at%symbol))
    at%SetIsotope = .TRUE.

    IF (at%Z == 1 .AND. at%A == 1) at%isotope='H'
    IF (at%Z == 1 .AND. at%A == 2) at%isotope='D'
    IF (at%Z == 1 .AND. at%A == 3) at%isotope='T'

    IF (at%Z > 0 .AND. abs(real(at%A,kind=Rkind)/at%mass-ONE) > TEN**(-2)) THEN
      write(out_unitp,*) '  WARNNING in Read_atom'
      write(out_unitp,*) '  The atomic mass in g/mol is probably too different from A (> 1%)'
      write(out_unitp,*) '  mass, A: ',at%mass,at%A
      write(out_unitp,*) '  atom: ',at
      !STOP
    END IF
    IF (at%abundance > ONE .AND. at%abundance < ZERO) THEN
      write(out_unitp,*) '  ERROR in Read_atom'
      write(out_unitp,*) '  The abundance is greater than ONE or lower than ZERO!!'
      write(out_unitp,*) ' Z,symbol,A,mass,abundance',at%Z,at%symbol,at%A,at%mass,at%abundance
      STOP
    END IF
    !write(out_unitp,*) 'at: ',at

  END SUBROUTINE Read2018_atom

!> @brief Subroutine which initializes the isotpic masses from
!!        from <a href="https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses">NIST</a>.
!!
!> @author David Lauvergnat
!! @date 29/11/2018
!!
!! @param mendeleev     table of "atoms" with known isotopes
!! @param auTOgPmol     conversion factor: au to g.mol-1
!! @param mass_unit     name of the mass unit (eg gPmol)
  SUBROUTINE construct_table_at_NIST2012(mendeleev,auTOgPmol,mass_unit)

    TYPE (table_atom), intent(inout) :: mendeleev
    character (len=*), intent(in)    :: mass_unit
    real (kind=Rkind), intent(in)    :: auTOgPmol

    integer           :: nio,err
    integer           :: memory
    TYPE (atom)       :: at
    integer           :: Z,A,err_mem,A_OF_Max_abundance
    real (kind=Rkind) :: Max_abundance

    character (len=*), parameter ::                               &
       alphab1="ABCBEFGHIJKLMNOPQRSTUVWXYZ"
    character (len=*), parameter ::                               &
       alphab2="abcdefghijklmnopqrstuvwxyz"

    integer :: i,pos1,pos2,pos3,isot

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='construct_table_at_NIST2012'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
      END IF
!---------------------------------------------------------------------

    mendeleev%List_Isotope%name = trim(EVRT_path) // '/' //       &
                                  'Internal_data/IsotopicMass_2012.txt'

    IF (.NOT. mendeleev%construct) THEN

      CALL file_open(mendeleev%List_Isotope,nio,old=.TRUE.,err_file=err)
      IF (err /= 0) THEN
        write(out_unitp,*) ' WARNNING in ',name_sub
        write(out_unitp,*) ' The file "Internal_data/IsotopicMass_2012.txt" can not be open'
        write(out_unitp,*) ' Two main raisons:'
        write(out_unitp,*) '  -1- the ElVibRot directory has been moved after the compilation'
        write(out_unitp,*) '  -2- the "EVRT_path" is badly defined in the namelist "system"'
        write(out_unitp,*) '      EVRT_path: ',trim(EVRT_path)

        write(out_unitp,*) ' Two solutions:'
        write(out_unitp,*) '  -1- recompile ElVibRot: "make clean ; make"'
        write(out_unitp,*) '  -2- set-up the "EVRT_path" in the namelist "system"'

        write(out_unitp,*) ' WARNNING the old table (construct_table_at_HandBook70ed) is used'
        CALL construct_table_at_HandBook70ed(mendeleev,auTOgPmol,mass_unit)
        RETURN
      END IF

      mendeleev%construct = .TRUE.

      mendeleev%max_Z = 110
      mendeleev%max_A = 300
      memory = product( (/ 1+mendeleev%max_Z,1+mendeleev%max_A /) &
                                                                 )
      allocate(mendeleev%at(0:mendeleev%max_Z,0:mendeleev%max_A), &
                                                     stat=err_mem) ! change alloc done
      CALL error_memo_allo(err_mem,memory,"mendeleev%at",name_sub,"atom")


      DO Z=0,mendeleev%max_Z
      DO A=0,mendeleev%max_A
        mendeleev%at(Z,A) = atom(' ',' ',Z,A,-ONE,'g/mol',ZERO,ZERO,-ONE,.FALSE.,.FALSE.)
      END DO
      END DO

      read(nio,*)   ! first line
      isot = 0
      DO
        CALL Read_atom(at,nio,err)
        IF (err /= 0) EXIT
        isot = isot + 1
        Z    = at%Z
        A    = at%A
        IF (Z < lbound(mendeleev%at,dim=1) .OR. Z > ubound(mendeleev%at,dim=1)) CYCLE
        IF (A < lbound(mendeleev%at,dim=2) .OR. A > ubound(mendeleev%at,dim=2)) CYCLE

        mendeleev%at(Z,A) = at
        IF (debug) write(out_unitp,*) 'at: ',mendeleev%at(Z,A)

      END DO
      write(out_unitp,*) 'The number of read isotopes:',isot
      CALL file_close(mendeleev%List_Isotope)

      DO Z=0,mendeleev%max_Z
        Max_abundance            = -ONE
        A_OF_Max_abundance       = -1
        DO A=0,mendeleev%max_A
          mendeleev%at(Z,A)%mass = mendeleev%at(Z,A)%mass/auTOgPmol
          mendeleev%at(Z,A)%mass_unit = mass_unit
          IF (mendeleev%at(Z,A)%abundance > Max_abundance) THEN
            Max_abundance            = mendeleev%at(Z,A)%abundance
            A_OF_Max_abundance       = A
          END IF
        END DO
        IF (A_OF_Max_abundance /= -1) mendeleev%at(Z,A_OF_Max_abundance)%MainIsotope = .TRUE.
      END DO


      IF (debug .OR. print_level > 1) CALL List_OF_table_at(mendeleev)
    END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------

  END SUBROUTINE construct_table_at_NIST2012
!> @brief Subroutine which initializes the isotpic masses download (24/11/2018) from
!!        from <a href="https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses">NIST</a>.
!!
!> @author David Lauvergnat
!! @date 29/11/2018
!!
!! @param mendeleev     table of "atoms" with known isotopes
!! @param auTOgPmol     conversion factor: au to g.mol-1
!! @param mass_unit     name of the mass unit (eg gPmol)
  SUBROUTINE construct_table_at_NIST2018(mendeleev,auTOgPmol,mass_unit)

    TYPE (table_atom), intent(inout) :: mendeleev
    character (len=*), intent(in)    :: mass_unit
    real (kind=Rkind), intent(in)    :: auTOgPmol

    integer           :: nio,err
    integer           :: memory
    TYPE (atom)       :: at
    integer           :: Z,A,err_mem
    real (kind=Rkind) :: Max_abundance

    character (len=*), parameter ::                               &
       alphab1="ABCBEFGHIJKLMNOPQRSTUVWXYZ"
    character (len=*), parameter ::                               &
       alphab2="abcdefghijklmnopqrstuvwxyz"

    integer :: i,pos1,pos2,pos3,isot,A_OF_Max_abundance

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='construct_table_at_NIST2018'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
      END IF
!---------------------------------------------------------------------

    mendeleev%List_Isotope%name = trim(EVRT_path) // '/' //       &
                                  'Internal_data/IsotopicMass_2018.txt'

    IF (.NOT. mendeleev%construct) THEN

      CALL file_open(mendeleev%List_Isotope,nio,old=.TRUE.,err_file=err)
      IF (err /= 0) THEN
        write(out_unitp,*) ' WARNNING in ',name_sub
        write(out_unitp,*) ' The file "Internal_data/IsotopicMass_2018.txt" can not be open'
        write(out_unitp,*) ' Two main raisons:'
        write(out_unitp,*) '  -1- the ElVibRot directory has been moved after the compilation'
        write(out_unitp,*) '  -2- the "EVRT_path" is badly defined in the namelist "system"'
        write(out_unitp,*) '      EVRT_path: ',trim(EVRT_path)

        write(out_unitp,*) ' Two solutions:'
        write(out_unitp,*) '  -1- recompile ElVibRot: "make clean ; make"'
        write(out_unitp,*) '  -2- set-up the "EVRT_path" in the namelist "system"'

        write(out_unitp,*) ' WARNNING the old table (construct_table_at_HandBook70ed) is used'
        CALL construct_table_at_HandBook70ed(mendeleev,auTOgPmol,mass_unit)
        RETURN
      END IF

      mendeleev%construct = .TRUE.

      mendeleev%max_Z = 118
      mendeleev%max_A = 400
      memory = product( (/ 1+mendeleev%max_Z,1+mendeleev%max_A /) )
      allocate(mendeleev%at(0:mendeleev%max_Z,0:mendeleev%max_A),stat=err_mem) ! change alloc done
      CALL error_memo_allo(err_mem,memory,"mendeleev%at",name_sub,"atom")

      DO Z=0,mendeleev%max_Z
      DO A=0,mendeleev%max_A
        mendeleev%at(Z,A) = atom(' ',' ',Z,A,-ONE,'g/mol',ZERO,ZERO,-ONE,.FALSE.,.FALSE.)
      END DO
      END DO

      IF (debug) write(out_unitp,*) '-----------------------------------------------------'
      read(nio,*)   ! first line
      isot = 0
      DO
        CALL Read2018_atom(at,nio,err)
        IF (debug)  write(out_unitp,*) 'Read2018_atom, err',err
        IF (err /= 0) EXIT
        isot = isot + 1
        Z = at%Z
        A = at%A
        IF (Z < lbound(mendeleev%at,dim=1) .OR. Z > ubound(mendeleev%at,dim=1)) CYCLE
        IF (A < lbound(mendeleev%at,dim=2) .OR. A > ubound(mendeleev%at,dim=2)) CYCLE

        mendeleev%at(Z,A) = at
        IF (debug) write(out_unitp,*) 'at: ',mendeleev%at(Z,A)
        IF (debug) write(out_unitp,*) '-----------------------------------------------------'

      END DO
      write(out_unitp,*) 'The number of read isotopes:',isot
      CALL file_close(mendeleev%List_Isotope)


      DO Z=0,mendeleev%max_Z
        Max_abundance            = -ONE
        A_OF_Max_abundance       = -1
        DO A=0,mendeleev%max_A
          mendeleev%at(Z,A)%mass = mendeleev%at(Z,A)%mass/auTOgPmol
          mendeleev%at(Z,A)%mass_unit = mass_unit
          IF (mendeleev%at(Z,A)%abundance > Max_abundance) THEN
            Max_abundance            = mendeleev%at(Z,A)%abundance
            A_OF_Max_abundance       = A
          END IF
        END DO
        IF (A_OF_Max_abundance /= -1) mendeleev%at(Z,A_OF_Max_abundance)%MainIsotope = .TRUE.
      END DO

      IF (debug .OR. print_level > 1) CALL List_OF_table_at(mendeleev)

    END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------

  END SUBROUTINE construct_table_at_NIST2018
!> @brief Subroutine which initializes the isotpic masses from
!!        the Handbook of Chemistry and Physics 70th edition (B-228).
!!
!> @author David Lauvergnat
!! @date 29/11/2018
!!
!! @param mendeleev     table of "atoms" with known isotopes
!! @param auTOgPmol     conversion factor: au to g.mol-1
!! @param mass_unit     name of the mass unit (eg gPmol)
  SUBROUTINE construct_table_at_HandBook70ed(mendeleev,auTOgPmol,mass_unit)

    TYPE (table_atom), intent(inout) :: mendeleev
    character (len=*), intent(in)    :: mass_unit
    real (kind=Rkind), intent(in)    :: auTOgPmol

    integer :: Z,A
    integer :: err_mem,memory

!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='construct_table_at_HandBook70ed'
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
      CALL error_memo_allo(err_mem,memory,"mendeleev%at",name_sub,"atom")


      DO Z=0,mendeleev%max_Z
      DO A=0,mendeleev%max_A
        mendeleev%at(Z,A) = atom(' ',' ',Z,A,-ONE,'g/mol',ZERO,ZERO,-ONE,.FALSE.,.FALSE.)
      END DO
      END DO

      CALL construct_atom(                                        &
         mendeleev%at(  0,  0),isotope='X',mass=ZERO,MainIsotope=.TRUE.)

      CALL construct_atom(                                        &
         mendeleev%at(  1,  1),isotope='H',mass=1.007825_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
         mendeleev%at(  1,  2),isotope='D',mass=2.014000_Rkind)
      CALL construct_atom(                                        &
         mendeleev%at(  1,  3),isotope='T',mass=3.016050_Rkind)

      CALL construct_atom(                                        &
         mendeleev%at(  2,  3),isotope='3He',mass=3.016030_Rkind)
      CALL construct_atom(                                        &
         mendeleev%at(  2,  4),isotope='4He',mass=4.002600_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
         mendeleev%at(  2,  5),isotope='5He',mass=5.012220_Rkind)
      CALL construct_atom(                                        &
         mendeleev%at(  2,  6),isotope='6He',mass=6.018886_Rkind)
      CALL construct_atom(                                        &
         mendeleev%at(  2,  8),isotope='8He',mass=8.033920_Rkind)

      CALL construct_atom(                                        &
         mendeleev%at(  3,  5),isotope='5Li',mass=5.012540_Rkind)
      CALL construct_atom(                                        &
         mendeleev%at(  3,  6),isotope='6Li',mass=6.015210_Rkind)
      CALL construct_atom(                                        &
         mendeleev%at(  3,  7),isotope='7Li',mass=7.016003_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
         mendeleev%at(  3,  8),isotope='8Li',mass=8.022485_Rkind)
      CALL construct_atom(                                        &
         mendeleev%at(  3,  9),isotope='9Li',mass=9.026789_Rkind)
      CALL construct_atom(                                        &
         mendeleev%at(  4,  9),isotope='9Be',mass=9.012182_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
         mendeleev%at(  5, 10),isotope='10B',mass=10.012937_Rkind)
      CALL construct_atom(                                        &
         mendeleev%at(  5, 11),isotope='11B',mass=11.009305_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
         mendeleev%at(  6, 12),isotope='12C',mass=12.000000_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
         mendeleev%at(  7, 14),isotope='14N',mass=14.003074_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
         mendeleev%at(  8, 16),isotope='16O',mass=15.994915_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
         mendeleev%at(  9, 19),isotope='19F',mass=18.998403_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
        mendeleev%at( 10, 20),isotope='20Ne',mass=19.992435_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
        mendeleev%at( 10, 22),isotope='22Ne',mass=21.991383_Rkind)

      CALL construct_atom(                                        &
        mendeleev%at( 11, 23),isotope='23Na',mass=22.989767_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
        mendeleev%at( 12, 24),isotope='24Mg',mass=23.985042_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
        mendeleev%at( 12, 25),isotope='25Mg',mass=24.985837_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 12, 26),isotope='26Mg',mass=25.982593_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 13, 27),isotope='27Al',mass=26.981540_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
        mendeleev%at( 14, 28),isotope='28Si',mass=27.976927_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
        mendeleev%at( 14, 29),isotope='29Si',mass=28.976495_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 14, 30),isotope='30Si',mass=29.97370_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 15, 31),isotope='31P',mass=30.973762_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
        mendeleev%at( 16, 32),isotope='32S',mass=31.972070_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
        mendeleev%at( 17, 35),isotope='35Cl',mass=34.968852_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
        mendeleev%at( 17, 37),isotope='37Cl',mass=36.965903_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 18, 40),isotope='40Ar',mass=39.962384_Rkind,MainIsotope=.TRUE.)

      CALL construct_atom(                                        &
        mendeleev%at( 19, 39),isotope='39K',mass=38.963707_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
        mendeleev%at( 20, 40),isotope='40Ca',mass=39.962591_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
        mendeleev%at( 21, 45),isotope='45Sc',mass=44.955910_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
        mendeleev%at( 22, 48),isotope='48Ti',mass=47.947947_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
        mendeleev%at( 23, 51),isotope='51V',mass=50.943962_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
        mendeleev%at( 24, 52),isotope='52Cr',mass=51.940509_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
        mendeleev%at( 24, 52),isotope='53Cr',mass=52.940651_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 25, 55),isotope='55Mn',mass=54.938047_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
        mendeleev%at( 26, 56),isotope='56Fe',mass=55.934939_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
        mendeleev%at( 27, 59),isotope='59Co',mass=58.933198_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
        mendeleev%at( 28, 58),isotope='58Ni',mass=57.935346_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
        mendeleev%at( 28, 60),isotope='60Ni',mass=59.930788_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 29, 63),isotope='63Cu',mass=62.939598_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
        mendeleev%at( 29, 65),isotope='65Cu',mass=64.927793_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 30, 64),isotope='64Zn',mass=63.929145_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
        mendeleev%at( 30, 66),isotope='66Zn',mass=65.926034_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 30, 68),isotope='68Zn',mass=67.924846_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 31, 69),isotope='69Ga',mass=68.925580_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
        mendeleev%at( 31, 71),isotope='71Ga',mass=70.924700_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 32, 70),isotope='70Ge',mass=69.924250_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 32, 72),isotope='72Ge',mass=71.922079_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 32, 74),isotope='74Ge',mass=73.921177_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
        mendeleev%at( 33, 75),isotope='75As',mass=74.921594_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
        mendeleev%at( 34, 78),isotope='78Se',mass=78.000000_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 34, 80),isotope='80Se',mass=79.916520_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
        mendeleev%at( 35, 79),isotope='79Br',mass=78.918336_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
        mendeleev%at( 35, 81),isotope='81Br',mass=80.916289_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 36, 82),isotope='82Kr',mass=81.913482_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 36, 83),isotope='83Kr',mass=82.914135_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 36, 84),isotope='84Kr',mass=83.911507_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
        mendeleev%at( 36, 86),isotope='86Kr',mass=85.910616_Rkind)

      CALL construct_atom(                                        &
        mendeleev%at( 78, 194),isotope='194Pt',mass=193.962655_Rkind)
      CALL construct_atom(                                        &
        mendeleev%at( 78, 195),isotope='195Pt',mass=194.964766_Rkind,MainIsotope=.TRUE.)
      CALL construct_atom(                                        &
        mendeleev%at( 78, 196),isotope='196Pt',mass=195.964926_Rkind)


      DO Z=0,mendeleev%max_Z
      DO A=0,mendeleev%max_A
        IF (.NOT. mendeleev%at(Z,A)%SetIsotope) CYCLE
        mendeleev%at(Z,A)%mass = mendeleev%at(Z,A)%mass/auTOgPmol
        mendeleev%at(Z,A)%mass_unit = mass_unit
        IF (debug) write(out_unitp,*) mendeleev%at(Z,A)
      END DO
      END DO

      IF (debug .OR. print_level > 1) CALL List_OF_table_at(mendeleev)

    END IF

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------

  END SUBROUTINE construct_table_at_HandBook70ed

!> @brief Function: enables to get an isotopic mass
!!
!> @author David Lauvergnat
!! @date 29/11/2018
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
    integer,           intent(inout), optional :: err_mass !< to handle error (err_mass=0 => no error)
    character (len=*), intent(in),    optional :: name !< string which contains the mass (a real), the atomic symbol, or "Z_A" or "A"symbol
    TYPE (table_atom), intent(in)              :: mendeleev !< table of derived type "atom" with the known isotopes

    integer                         :: err_mass_loc
    integer                         :: AA,ZZ
    integer                         :: pos,i,clen,err
    real (kind=Rkind)               :: mass
    character (len=Name_len)        :: symb,symb_at
    character (len=:), allocatable  :: name2
    logical                         :: find,MainIsotope

    err_mass_loc = 0
    IF (present(err_mass)) err_mass = 0

    get_mass_Tnum = -ONE
    !IF (present(A)) write(out_unitp,*) 'A present',A
    !IF (present(Z)) write(out_unitp,*) 'Z present',Z
    !IF (present(name)) write(out_unitp,*) 'name present: ',name
    !flush(out_unitp)

    AA = -1
    IF (present(A)) AA = A
    ZZ = -1
    IF (present(Z)) ZZ = Z

    MainIsotope = (AA == -1)

    IF (present(A) .AND. .NOT. present(Z) .AND. .NOT. present(name)) THEN
      write(out_unitp,*) ' ERROR : get_mass_Tnum'
      write(out_unitp,*) ' Only A is present !!!'
      IF (present(err_mass)) THEN
            err_mass = -1
      ELSE
            STOP ' ERROR in : get_mass_Tnum'
      END IF
      RETURN
    END IF



    IF ( present(name) ) THEN
       name2 = String_TO_String(name)

      ! name has the value of the mass....
      IF (index(name2,".") .NE. 0) THEN
        !============================================================
        ! isotope mass read as real number
        !============================================================
        read(name2,*,IOSTAT=err_mass_loc) mass
        IF (err_mass_loc /= 0) THEN
          write(out_unitp,*) ' ERROR in : get_mass_Tnum'
          write(out_unitp,*) ' I CANNOT read the mass in "',trim(adjustl(name2)),'"'
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

      ELSE IF (index(name2,"_") > 0) THEN
        !============================================================
        ! isotope defined as : Z_A
        !============================================================
        MainIsotope = .FALSE.
        pos = index(name2,"_")
        name2(pos:pos) = ' ' ! change '_' in ' '
        read(name2,*,IOSTAT=err_mass_loc) ZZ,AA

        IF (err_mass_loc /= 0) THEN
          mass = -ONE
          ZZ   = -1
          AA   = -1
        END IF
      ELSE
        !============================================================
        ! isotope defined as:
        !       AX (2H or 13C or 28Si)
        ! or
        !      just as X (most abundant isotope)
        !
        ! where X is the symbol
        !============================================================
        ! find the end of the position of "A" in name
        !change name for "D" and "T"
        CALL string_uppercase_TO_lowercase(name2)
        IF (name2 == 'd') name2 = '2h'
        IF (name2 == 't') name2 = '3h'

        pos = 0
        DO
          read(name2(pos+1:pos+1),'(i1)',IOSTAT=err) i
          IF (err /= 0) EXIT
          pos = pos + 1
        END DO
        MainIsotope = (pos == 0)

        ! first the symbol, then ZZ
        symb = name2(pos+1:len(name2))
        !write(6,*) 'pos,MainIsotope,ZZ,AA,symb',pos,MainIsotope,ZZ,AA,symb

        DO ZZ=0,mendeleev%max_Z
          DO AA=0,mendeleev%max_A ! find the first isotope with Z
            IF (mendeleev%at(ZZ,AA)%SetIsotope) EXIT
          END DO
          IF (AA > mendeleev%max_A) CYCLE

          symb_at = mendeleev%at(ZZ,AA)%symbol
          CALL string_uppercase_TO_lowercase(symb_at)
          IF (symb == symb_at) EXIT
        END DO
        AA = -1

        ! Second AA, if need (pos > 0)
        IF (pos > 0) THEN ! isotope defined as: AX (2H or 13C or 28Si)
          read(name2(1:pos),*,IOSTAT=err_mass_loc) AA
        END IF
        !write(6,*) 'pos,MainIsotope,ZZ,AA,symb',pos,MainIsotope,ZZ,AA,symb

        deallocate(name2)
      END IF
    END IF

    IF (ZZ < 0 .OR. ZZ > ubound(mendeleev%at,dim=1)) THEN
       write(out_unitp,*) ' ERROR in : get_mass_Tnum'
       write(out_unitp,*) 'ZZ,AA',ZZ,AA
       write(out_unitp,*) ' I CANNOT find Z in "',trim(adjustl(name)),'"'
       IF (present(Z)) THEN
         write(out_unitp,*) ' ... or Z (from the argument) is out of range'
       END IF
       err_mass_loc = -1
    END IF

    IF (err_mass_loc == 0) THEN
      IF (MainIsotope) THEN
        DO AA=0,mendeleev%max_A ! find the first isotope with Z
          IF (mendeleev%at(ZZ,AA)%MainIsotope) EXIT
        END DO
        !write(6,*) 'pos,MainIsotope,ZZ,AA,symb',pos,MainIsotope,ZZ,AA,symb
        !flush(6)

      ELSE IF (AA < 0 .OR. AA > ubound(mendeleev%at,dim=2)) THEN
        write(out_unitp,*) ' ERROR in : get_mass_Tnum'
        write(out_unitp,*) 'ZZ,AA',ZZ,AA
        write(out_unitp,*) ' I CANNOT read A in "',trim(adjustl(name)),'"'
        IF (present(A)) THEN
          write(out_unitp,*) ' ... or A (from the argument) is out of range'
        END IF
        err_mass_loc = -1
      ELSE IF (.NOT. mendeleev%at(ZZ,AA)%SetIsotope) THEN
        write(out_unitp,*) 'ZZ,AA',ZZ,AA
        write(out_unitp,*) ' ERROR in : get_mass_Tnum'
        write(out_unitp,*) ' This isotope is not defined in "',trim(adjustl(name)),'"'
        err_mass_loc = -1
      END IF
    END IF


    IF (present(Z)) Z = ZZ
    IF (present(A)) A = AA

    IF (err_mass_loc == 0) THEN
      !write(out_unitp,*) 'atom: ',mendeleev%at(ZZ,AA)
      IF (print_level > 1)                                        &
        write(out_unitp,*) 'get_mass_Tnum: ',ZZ,AA,mendeleev%at(ZZ,AA)%mass

      get_mass_Tnum = mendeleev%at(ZZ,AA)%mass
    ELSE
          write(out_unitp,*) ' ERROR : get_mass_Tnum'
          write(out_unitp,*) ' I CANNOT get the right isotope'
          write(out_unitp,*) ' Your atom is NOT in my list !!'
          write(out_unitp,*) ' Your atom: "',trim(name),'"'
          CALL List_OF_table_at(mendeleev)
      IF (present(err_mass)) THEN
        mass = -ONE
        err_mass = -1
      ELSE
       STOP ' ERROR in : get_mass_Tnum'
      END IF
    END IF

  END FUNCTION get_mass_Tnum
!> @brief Subroutine which deallocates the derived type table_atom
!!
!> @author David Lauvergnat
!! @date 29/11/2018
!! @param mendeleev     table of "atoms" with known isotopes
  SUBROUTINE dealloc_table_at(mendeleev)
    TYPE (table_atom), intent(inout) :: mendeleev

    integer :: err_mem,memory

    memory = size(mendeleev%at)
    deallocate(mendeleev%at,stat=err_mem) ! change alloc done
    CALL error_memo_allo(err_mem,-memory,"mendeleev%at","dealloc_table_at","atom")
    nullify(mendeleev%at)

  END SUBROUTINE dealloc_table_at
!> @brief Subroutine which write the list of isotopes from table_atom
!!
!> @author David Lauvergnat
!! @date 29/11/2018
!! @param mendeleev     table of "atoms" with known isotopes
  SUBROUTINE List_OF_table_at(mendeleev)

    TYPE (table_atom), intent(in) :: mendeleev


    integer           :: Z,A
    logical           :: SetIsotope


!----- for debuging --------------------------------------------------
      character (len=*), parameter :: name_sub='List_OF_table_at'
      logical, parameter :: debug = .FALSE.
      !logical, parameter :: debug = .TRUE.
!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING ',name_sub
      END IF
!---------------------------------------------------------------------

    IF (.NOT. mendeleev%construct) RETURN

      write(out_unitp,*) '---------------------------------------------'
      write(out_unitp,*) 'Isotopes list:'

      DO Z=0,mendeleev%max_Z

        !check if at least one isotope is defined
        SetIsotope = .FALSE.
        DO A=0,mendeleev%max_A
          SetIsotope = mendeleev%at(Z,A)%SetIsotope
          IF (SetIsotope) EXIT
        END DO

        IF (SetIsotope) THEN
          write(out_unitp,'(i4,":")',advance='no') Z

          DO A=0,mendeleev%max_A
            IF (.NOT. mendeleev%at(Z,A)%SetIsotope) CYCLE
            IF (Z == 0) THEN
              write(out_unitp,'(x,a)',advance='no') trim(mendeleev%at(Z,A)%symbol)
            ELSE IF (Z == 1 .AND. (A == 2 .OR. A == 3)) THEN
              write(out_unitp,'(x,a)',advance='no') trim(mendeleev%at(Z,A)%symbol)
              write(out_unitp,'(" (or ",i0,a,")")',advance='no') A,trim(mendeleev%at(Z,1)%symbol)
            ELSE
              write(out_unitp,'(x,i0,a)',advance='no') A,trim(mendeleev%at(Z,A)%symbol)
            END IF
          END DO
          write(out_unitp,*)
        END IF
      END DO
      write(out_unitp,*) 'END Isotopes list'
      write(out_unitp,*) '---------------------------------------------'

!---------------------------------------------------------------------
      IF (debug) THEN
        write(out_unitp,*) 'END ',name_sub
      END IF
!---------------------------------------------------------------------

  END SUBROUTINE List_OF_table_at

END MODULE mod_Atom
