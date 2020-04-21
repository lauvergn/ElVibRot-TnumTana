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
!> @mainpage Physical Constants
!! This module enables to use fundamental physical constants
!! (speed of light in vacuum, Planck constant ... and isotopic masses).
!! Three versions can be selected:
!!
!! @li The CODATA 2014 ones, downloaded from
!! <a href="https://physics.nist.gov/cuu/Constants/index.html">NIST</a>
!! and the NIST <a href="https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses">masses</a>
!!  downloaded in 2018 (not the default).
!!
!! @li The CODATA 2006 ones, downloaded from
!! <a href="https://physics.nist.gov/cuu/Constants/archive2006.html">NIST</a>
!! and the NIST <a href="https://www.nist.gov/pml/atomic-weights-and-isotopic-compositions-relative-atomic-masses">masses</a>
!!  downloaded in 2012 (default).
!!
!! @li Constants and masses from the 70th edition of the Handbook of Chemistry and Physics.
!!
!! From these fundamental constants, some conversion factors are calculated automatically
!! and can be used easily.
!!
!! @remark The actual mass values of the NIST web page differ slightly from the module ones.
!!
!! @author David Lauvergnat
!! @date 29/11/2018
!!
!! @licence GNU Lesser General Public License
!!
!! @section install_sec Installation
!!
!! Dependencies: this module needs the fortran modules in the @e Source_Lib/sub_system directory.
!!
!! Build the module (with dependencies):
!!
!!     make PhysConst
!!
!! Build the module documentation (with doxygen):
!!
!!     make doxy
!!
!! @section test_sec Tests
!!
!! Example data/script files:
!!     Examples/exa_PhysicalConstants/dat_PhysConst_NIST2018
!!     Examples/exa_PhysicalConstants/dat_PhysConst_NIST2012
!!     Examples/exa_PhysicalConstants/dat_PhysConst_HandBook70ed
!!
!! To test the installation, you can run the test examples.
!!
!!     cd Examples/exa_PhysicalConstants ; ./run_tests
!!
!! The results will be compared to previous results in Examples/exa_PhysicalConstants/output_29nov2018
!!
!! @section Dev_sec Developer use
!!     The main program "PhysicalConstants" shows:
!!     @li how to initialize physical constants and masses
!!     @li examples to use physical constants and conversion factors
!!     @li examples to use and read isotopic masses
!!
!! @subsection Init_Dev_sec Initialization
!! Declaration of \a constant derived type
!!
!!      TYPE(constant) :: const_phys
!!
!! where @a const_phys is a derived type which contains all the constants and the masses.
!!
!! Call the subroutine:
!!
!!      CALL sub_constantes(const_phys,Read_Namelist=.TRUE.)
!!
!! When @a Read_Namelist is set to .TRUE., the namelist "constantes" will be read.
!!
!!
!! @subsection Const_Dev_sec Use of some constants
!!  @li  const_phys\%c          :    Speed of light (in m.s-1)
!!  @li  const_phys\%h          :    Planck constant (in J.s)
!!  @li  const_phys\%e          :    Electron charge (in C)
!!  @li  const_phys\%me         :    Electron mass (in kg)
!!  @li  const_phys\%a0         :    Bohr radius (in m)
!!  @li  const_phys\%Eh         :    Hartree constant (in J)
!!  @li  const_phys\%auTOcm_inv :    Conversion factor: au to cm-1
!!  @li  and many others ...
!!
!! @subsection Mass_Dev_sec Get mass
!! @li with the atomic symbol
!!
!!     mass = get_mass_Tnum(const_phys%mendeleev,name='H')
!!     mass = get_mass_Tnum(const_phys%mendeleev,name='D')
!! @li with the atomic symbol plus the number of nucleons
!!
!!     mass = get_mass_Tnum(const_phys%mendeleev,name='12C')
!! @li with the number of electrons (get the most abundant isotope)
!!
!!     mass = get_mass_Tnum(const_phys%mendeleev,Z=6)
!! @li with the numbers of electrons and nucleons
!!
!!     mass = get_mass_Tnum(const_phys%mendeleev,Z=6,A=12)

  PROGRAM PhysicalConstants
    USE mod_system
    USE mod_constant
    USE mod_RealWithUnit
    IMPLICIT NONE

    !- parameters for para_Tnum -----------------------
    TYPE (constant)  :: const_phys

    real(kind=Rkind) :: mass
    integer          :: Z,A,err_mass
    character (len=Name_len) :: mass_name
    !- working parameters ------------------------------------------
    integer :: err_read
    character (len=*), parameter :: name_sub='PhysicalConstants'

    !=======================================================================
    !=======================================================================
    write(out_unitp,*) 'BEGINNING ',name_sub
    print_level=0

    write(out_unitp,*) '==========================================='
    write(out_unitp,*) '= Module test: RealWithUnit(RWU) =========='
    CALL Test_RWU()
    write(out_unitp,*) '==========================================='


    write(out_unitp,*) '==========================================='
    write(out_unitp,*) '= Usefull conversion factors or constants ='
    CALL sub_constantes(const_phys,Read_Namelist=.TRUE.)

 11 format (a,e17.10)
 21 format (a,f18.6)
    write(out_unitp,*)
    write(out_unitp,*) 'pi =                ',const_phys%pi
    write(out_unitp,*) 'cos(pi) =           ',cos(const_phys%pi)
    write(out_unitp,*)
    write(out_unitp,11) 'au => m            ',const_phys%a0
    write(out_unitp,11) 'au => Angstrom     ',const_phys%a0 * TEN**10
    write(out_unitp,*)

    write(out_unitp,11) ' au => s           ',const_phys%Ta
    write(out_unitp,21) ' au => fs          ',const_phys%Ta*TEN**15
    write(out_unitp,*)

    write(out_unitp,11) ' au => J           ',const_phys%Eh
    write(out_unitp,21) ' au => cm-1        ',const_phys%auTOcm_inv
    write(out_unitp,21) ' au => eV          ',const_phys%auTOeV
    write(out_unitp,*)

    write(out_unitp,21) ' g.mol-1 => au     ',const_phys%inv_Name/TEN**3
    write(out_unitp,11) ' Debye => au       ',const_phys%convDebyeTOau
    write(out_unitp,*)

    write(out_unitp,11) ' au => V.cm-1 (E0) ',const_phys%E0
    write(out_unitp,11) ' au => W.cm-2 (I0) ',const_phys%I0
    write(out_unitp,*)

    write(out_unitp,*) '==========================================='
    CALL flush_perso(out_unitp)

    write(out_unitp,*) '==========================================='
    write(out_unitp,*) '====== TEST to get masses ================='
    CALL flush_perso(out_unitp)

    write(out_unitp,*) ' -------------------------------------- '
    mass = get_mass_Tnum(const_phys%mendeleev,name='X')
    write(out_unitp,*) 'mass of X in au',mass
    CALL flush_perso(out_unitp)

    write(out_unitp,*) ' -------------------------------------- '
    mass = get_mass_Tnum(const_phys%mendeleev,name='100.')
    write(out_unitp,*) 'mass of "100." in au',mass
    CALL flush_perso(out_unitp)

    write(out_unitp,*) ' -------------------------------------- '
    mass = get_mass_Tnum(const_phys%mendeleev,name='H')
    write(out_unitp,*) 'mass of H in au',mass
    CALL flush_perso(out_unitp)

    write(out_unitp,*) ' -------------------------------------- '
    mass = get_mass_Tnum(const_phys%mendeleev,name='1_2')
    write(out_unitp,*) 'mass of "1_2" in au',mass
    CALL flush_perso(out_unitp)

    write(out_unitp,*) ' -------------------------------------- '
    mass = get_mass_Tnum(const_phys%mendeleev,name='D')
    write(out_unitp,*) 'mass of D in au',mass
    CALL flush_perso(out_unitp)

    write(out_unitp,*) ' -------------------------------------- '
    mass = get_mass_Tnum(const_phys%mendeleev,name='T')
    write(out_unitp,*) 'mass of T in au',mass
    CALL flush_perso(out_unitp)

    write(out_unitp,*) ' -------------------------------------- '
    mass = get_mass_Tnum(const_phys%mendeleev,name='1_4',err_mass=err_mass)
    write(out_unitp,*) 'mass of "1_4" in au',mass
    CALL flush_perso(out_unitp)

    write(out_unitp,*) ' -------------------------------------- '
    mass = get_mass_Tnum(const_phys%mendeleev,name='He')
    write(out_unitp,*) 'mass of He in au',mass
    CALL flush_perso(out_unitp)

    write(out_unitp,*) ' -------------------------------------- '
    mass = get_mass_Tnum(const_phys%mendeleev,name='HE')
    write(out_unitp,*) 'mass of HE in au',mass
    CALL flush_perso(out_unitp)

    write(out_unitp,*) ' -------------------------------------- '
    mass = get_mass_Tnum(const_phys%mendeleev,name='he')
    write(out_unitp,*) 'mass of he in au',mass
    CALL flush_perso(out_unitp)

    write(out_unitp,*) ' -------------------------------------- '
    Z=-1
    A=-1
    mass = get_mass_Tnum(const_phys%mendeleev,Z,A,name='He')
    write(out_unitp,*) 'mass of He in au',mass,' and then Z and A',Z,A
    CALL flush_perso(out_unitp)

    write(out_unitp,*) ' -------------------------------------- '
    Z=6
    A=-1
    mass = get_mass_Tnum(const_phys%mendeleev,Z,A)
    write(out_unitp,*) 'mass of Z=6 in au',mass,' and then A',A
    CALL flush_perso(out_unitp)

    write(out_unitp,*) ' -------------------------------------- '
    Z=6
    A=13
    mass = get_mass_Tnum(const_phys%mendeleev,Z,A,err_mass=err_mass)
    write(out_unitp,*) 'mass of Z=6,A=13 in au',mass
    CALL flush_perso(out_unitp)

    write(out_unitp,*) '==========================================='
    DO
      read(in_unitp,*,IOSTAT=err_read) mass_name
      IF (err_read /= 0) EXIT
      mass = get_mass_Tnum(const_phys%mendeleev,name=mass_name,err_mass=err_mass)
      write(out_unitp,*) 'mass of "',trim(adjustl(mass_name)),'" in au',mass

      write(out_unitp,*) ' -------------------------------------- '
      CALL flush_perso(out_unitp)

    END DO


    write(out_unitp,*) '==========================================='

    write(out_unitp,*) 'END ',name_sub

  END PROGRAM PhysicalConstants
