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
!> This fortran module enables to initialyze and use physical constants and masses.
!> Only the "sub_constantes", "get_mass_Tnum" (from the module mod_Atom) subroutines can used.
!> \author David Lauvergnat
!> \date 29/11/2018
 MODULE mod_Constant
 use mod_system
 use mod_Atom, only: table_atom, dealloc_table_at, get_mass_Tnum,        &
                     construct_table_at_NIST2012, construct_table_at_NIST2018,&
                     construct_table_at_HandBook70ed

 USE mod_RealWithUnit
 IMPLICIT NONE

 PRIVATE


!> This derived type contains some fundamental physical constants, conversion factors and isotopic masses (in the mendeleev variable)
 TYPE constant

   logical           :: constant_done = .FALSE. !< Flag to check the initialization
   real (kind=Rkind) :: pi                      !<  pi=3.14159....

   real (kind=Rkind) :: auTOcm_inv              !< Conversion factor: au to cm-1
   real (kind=Rkind) :: auTOeV                  !< Conversion factor: au to eV
   real (kind=Rkind) :: auTOGHz                 !< Conversion factor: au to GHz
   real (kind=Rkind) :: auTOkJmol_inv           !< Conversion factor: au to kJ.mol-1
   real (kind=Rkind) :: auTOkcalmol_inv         !< Conversion factor: au to kcal.mol-1

   real (kind=Rkind) :: auTOenergy              !< Conversion factor: au to selected energy unit (default cm-1)
   character(len=Name_len) :: ene_unit          !< The name of selected energy unit (default cm-1)

   real (kind=Rkind) :: inv_Name         !<  Conversion factor: kg/mol => au

   real (kind=Rkind) :: c                !<  Speed of light (exact) (in m s-1)
   real (kind=Rkind) :: mhu0             !<  Vacuum permeability (exact) (in N.A-2)
   real (kind=Rkind) :: epsi0            !<  Vacuum permittivity (exact) (in F.m-1)
   real (kind=Rkind) :: G                !<  Gravitational constant (in m3.kg-1.s-2)
   real (kind=Rkind) :: h,hb             !<  Planck constant (h and hb=h/2pi) (in J.s)
   real (kind=Rkind) :: e                !<  Electron charge (in C)
   real (kind=Rkind) :: me               !<  Electron mass (in kg)
   real (kind=Rkind) :: mp               !<  Proton mass (in kg)
   real (kind=Rkind) :: alpha            !<  Fine-structure constant: alpha (without unit)
   real (kind=Rkind) :: Na               !<  Avogadro number (in mol-1)
   real (kind=Rkind) :: R                !<  Ideal gas constant: R (in J.mol-1.K-1)
   real (kind=Rkind) :: k                !<  Boltzmann constant (in J.K-1)

   real (kind=Rkind) :: mhu              !<  Reduce mass of the hydrogen atom (in kg)
   real (kind=Rkind) :: a0               !<  Bohr radius (in m)
   real (kind=Rkind) :: Eh               !<  Hartree constant (in J) or conversion factor: au to J
   real (kind=Rkind) :: Ta               !<  au of time (in s) or conversion factor: au to s
   real (kind=Rkind) :: E0               !<  Conversion factor of an electric field in au to N m-1
   real (kind=Rkind) :: I0               !<  Conversion factor a plane wave intensity in au to W m-2

   real (kind=Rkind) :: convAif          !<  Conversion factor (Einstein coefficient, A) for Spontaneous emission: au to s-1 facteur de convertion pour le coef d'Einstein Aif
   real (kind=Rkind) :: convIDif         !<  Conversion factor for the electric-dipolar intenisties: au to m.mol-1
   real (kind=Rkind) :: convIQif         !<  Conversion factor for the electric-quadripolar intenisties: au to ???


   real (kind=Rkind) :: convDebyeTOau    !< Conversion factor: Debye TO au
   real (kind=Rkind) :: conv_auTOCm      !< Conversion factor: au TO C.m


   TYPE (table_atom)        :: mendeleev !< Mendeleev table with isotopic masses
   character(len=Name_len)  :: mass_unit !< Mass unit (default : au)
   real (kind=Rkind)        :: auTOmass  !< Converstion factor:  au TO g/mol (1/1822...)

 END TYPE constant
 !==============================================

 INTERFACE alloc_array
   MODULE PROCEDURE alloc_array_OF_Constantdim0
 END INTERFACE
 INTERFACE dealloc_array
   MODULE PROCEDURE dealloc_array_OF_Constantdim0
 END INTERFACE

 PUBLIC  :: table_atom,get_mass_Tnum,dealloc_table_at

 PUBLIC  :: constant,sub_constantes

 PUBLIC :: REAL_WU, RWU_Write, RWU_WriteUnit
 PUBLIC :: convRWU_TO_R_WITH_WorkingUnit,convRWU_TO_R_WITH_WritingUnit
 PUBLIC :: get_Conv_au_TO_WriteUnit,get_Conv_au_TO_Unit

 CONTAINS


  SUBROUTINE alloc_array_OF_Constantdim0(tab,name_var,name_sub)
    IMPLICIT NONE

    TYPE (constant),   pointer, intent(inout) :: tab
    character (len=*),          intent(in)    :: name_var,name_sub

    integer, parameter :: ndim=0
    logical :: memory_test

!----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub_alloc = 'alloc_array_OF_Constantdim0'
    integer :: err_mem,memory
    logical,parameter :: debug=.FALSE.
!    logical,parameter :: debug=.TRUE.
!----- for debuging --------------------------------------------------


    IF (associated(tab))                                             &
         CALL Write_error_NOT_null(name_sub_alloc,name_var,name_sub)

    memory = 1
    allocate(tab,stat=err_mem)
    CALL error_memo_allo(err_mem,memory,name_var,name_sub,'constant')

  END SUBROUTINE alloc_array_OF_Constantdim0
  SUBROUTINE dealloc_array_OF_Constantdim0(tab,name_var,name_sub)
  IMPLICIT NONE

  TYPE (constant), pointer, intent(inout) :: tab
  character (len=*), intent(in) :: name_var,name_sub

!- for debuging --------------------------------------------------
  character (len=*), parameter :: name_sub_alloc = 'dealloc_array_OF_Constantdim0'
  integer :: err_mem,memory
  logical,parameter :: debug=.FALSE.
!  logical,parameter :: debug=.TRUE.
!- for debuging --------------------------------------------------

   !IF (.NOT. associated(tab)) RETURN
   IF (.NOT. associated(tab))                                       &
         CALL Write_error_null(name_sub_alloc,name_var,name_sub)

   memory = 1
   deallocate(tab,stat=err_mem)
   CALL error_memo_allo(err_mem,-memory,name_var,name_sub,'constant')
   nullify(tab)

  END SUBROUTINE dealloc_array_OF_Constantdim0

!> Subroutine: Initializes the physical constants, conversion factors and masses
!> \author David Lauvergnt
!> \date 29/11/2018
!! \param const_phys  derived type with the physical constants, conversion factors and masses
!! \param Read_Namelist an optional logical flag to be able to read the namelist "constantes"
!
  SUBROUTINE sub_constantes(const_phys,Read_Namelist,version,mass_version)
  IMPLICIT NONE

  TYPE (constant),         intent(inout)            :: const_phys
  logical,                 intent(in),    optional  :: Read_Namelist
  character(len=*),        intent(in),    optional  :: version,mass_version


  logical :: Read_Namelist_loc
  integer :: err_read,err_unit

  real (kind=Rkind) :: c
  real (kind=Rkind) :: mhu0
  real (kind=Rkind) :: epsi0
  real (kind=Rkind) :: G
  real (kind=Rkind) :: h,hb
  real (kind=Rkind) :: e,me
  real (kind=Rkind) :: mp
  real (kind=Rkind) :: alpha
  real (kind=Rkind) :: Na
  real (kind=Rkind) :: R,k
  real (kind=Rkind) :: mhu
  real (kind=Rkind) :: a0,Eh,Ta,E0,I0

  real (kind=Rkind) :: convAif,convIDif,convCrSecif,convIQif
  real (kind=Rkind) :: convDebyeTOau,conv_auTOCm
  real (kind=Rkind) :: inv_Name

  real (kind=Rkind) :: auTOcm_inv,auTOenergy,auTOeV,auTOGHz,auTOkJmol_inv,auTOkcalmol_inv
  character(len=Name_len)  :: ene_unit
  character(len=Name_len)  :: Time_unit

  character(len=Name_len)  :: mass_unit  !  the energy unit (default : au)
  real (kind=Rkind)        :: auTOmass          !  au => g/mol (1/1822...)
  real (kind=Rkind)        :: gPmolTOmass       !  au => g/mol (1/1822...)
  character(len=Name_len)  :: version_loc,mass_version_loc


  integer :: i


!- for debuging --------------------------------------------------
  character (len=*), parameter :: name_sub='sub_constantes'
  logical, parameter :: debug = .FALSE.
  !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
  END IF
!-----------------------------------------------------------------
  const_phys%constant_done = .TRUE.


  IF (present(Read_Namelist)) THEN
    Read_Namelist_loc = Read_Namelist
  ELSE
    Read_Namelist_loc = .FALSE.
  END IF

  ! initialization with or without the namelist
  IF (present(version)) THEN
    version_loc = trim(adjustl(version))
  ELSE
    version_loc = 'CODATA2006'
  END IF
  IF (present(mass_version)) THEN
    mass_version_loc = trim(adjustl(mass_version))
  ELSE
    mass_version_loc = 'NIST2012'
  END IF

  ene_unit   = 'cm-1'
  auTOenergy = -ONE
  auTOcm_inv = -ONE

  mass_unit  = "au"
  auTOmass   = -ONE
  inv_Name   = -ONE

  time_unit  = "au"


  IF (Read_Namelist_loc) THEN
    CALL sub_ReadNMLconstantes(auTOcm_inv,ene_unit,auTOenergy,          &
                               inv_Name,mass_unit,auTOmass,             &
                               time_unit,                               &
                               mass_version_loc,version_loc)
  END IF

  write(out_unitp,*) 'EVRT_path: ',EVRT_path


  CALL string_uppercase_TO_lowercase(version_loc,lower=.FALSE.) ! conversion in capital letters
  SELECT CASE (version_loc)
  CASE ('CODATA2006')
    CALL constantes_CODATA2006(c,mhu0,G,h,e,me,mp,Na,R)
  CASE ('CODATA2014')
    CALL constantes_CODATA2014(c,mhu0,G,h,e,me,mp,Na,R)
  CASE ('CODATA2018')
    CALL constantes_CODATA2018(c,mhu0,G,h,e,me,mp,Na,R)
  CASE ('HANDBOOK70ED','HANDBOOK')
    CALL constantes_HandBook70ed(c,mhu0,G,h,e,me,mp,Na,R)
  CASE ('PUBLI2001')
    CALL constantes_HandBook70ed_2001(c,mhu0,G,h,e,me,mp,Na,R)
  CASE Default
    CALL constantes_CODATA2006(c,mhu0,G,h,e,me,mp,Na,R)
  END SELECT

  ! Vacuum permittivity (exact) (in F.m-1) :: OK CODATA2006
  epsi0 = ONE/(mhu0*c*c)
  ! Planck constant (in J.s)
  hb = h/(TWO*pi)
  ! Fine structure constant (without unit)
  alpha = mhu0*c*e*e/(TWO*h)
  ! Ideal gas constant: R (in J.mol-1.K-1)
  k = R/Na


  IF (print_level > 1)                                              &
      write(out_unitp,*) 'c,mhu0,epsi0,G,h,hb,e,me,mp,alpha,Na,R,k',&
                          c,mhu0,epsi0,G,h,hb,e,me,mp,alpha,Na,R,k

!   Atomic constants

!   Reduced mass of the Hydrogen atom
    mhu = me*(mp/(mp+me))
    IF (print_level > 0) write(out_unitp,*) 'mhu (g)',mhu
!   Bohr radius (in m)
    a0 = (hb/e)**2 /me * FOUR*pi*epsi0
    IF (print_level > 0) write(out_unitp,*) 'a0 (m)',a0
!   Atomic energy unit (hartree) (en J)
    Eh = (e*e/a0)/(FOUR*pi*epsi0)
    IF (print_level > 0) write(out_unitp,*) 'Eh (J)',Eh
!   temps en unite atomique (en s)
    Ta = hb/Eh
    IF (print_level > 0) write(out_unitp,*) 'Ta (s)',Ta
!   Mass conversion factor: kg.mol-1 => au
    IF (inv_Name < ZERO) inv_Name = ONE/(Na*me)

!   Energy conversion factor: Hartree (au) => cm-1
    IF (auTOcm_inv < ZERO) auTOcm_inv = ONE/(Ta*TWO*pi*c*HUNDRED)
!   Energy conversion factor: au => eV
    auTOeV     = Eh/e
!   Energy conversion factor: GHz => eV
    auTOGHz    = auTOcm_inv * ONETENTH**7 * c

!   Energy conversion factor: au => kJ.mol-1
    auTOkJmol_inv    = Eh * ONETENTH**3 * Na
!   Energy conversion factor: au => kcal.mol-1
    auTOkcalmol_inv    = Eh * ONETENTH**3 * Na /4.184_Rkind

!   Electric field conversion factor: au => V.cm-1
    E0 = Eh / (e * a0) / HUNDRED
    E0 = hb**2 / (e * a0**3 * me) / HUNDRED
!   Electric field intensity (plane wave) conversion factor: au => W.cm-2
    I0 = epsi0 * c * E0**2 / TWO

!   Conversion for the Einstein coefficient (Aif)
    convAif = (TWO**4 * pi**3)/(THREE*epsi0*h) * (e*e*a0*a0) * (Eh/(c*h))**3

!   Convertion for the dipolar intensity
    convIDif = (EIGHT*pi**3*Na)/(FOUR*pi*epsi0*THREE*h*c) * (Eh/(c*h)) * (e*a0)**2

!   Convertion for the dipolar intensity
    convCrSecif = pi/(hb*epsi0*c) * (Eh/hb) * (e*a0)**2

!   Convertion for the quadripolar intensity
    convIQif = (FOUR*pi**5*Na)/(FOUR*pi*epsi0*FIVE*h*c) * (Eh/(c*h))**3 *  (e*a0*a0)**2

!   Convertion for the dipole moment (Debye TO au)
    IF (print_level > 0) write(out_unitp,*) 'Debye TO C.m',TEN/(HUNDRED*c)/TEN**20

    convDebyeTOau = TEN/(HUNDRED*c)/TEN**20 / (e*a0)
    IF (print_level > 0) write(out_unitp,*) 'Debye TO au',convDebyeTOau

    conv_auTOCm = e*a0
!--------------------------------------------------------------

!--------------------------------------------------------------


  !------------------------------------------------------------------
  !     atomic mass of isotopes
  IF (mass_unit .EQ. "au") THEN
    IF (auTOmass < ZERO) const_phys%auTOmass  = ONE
    const_phys%mass_unit = "au"
  ELSE IF (mass_unit .EQ. "g/mol" .OR. mass_unit .EQ. "gPmol" .OR. mass_unit .EQ. "g.mol-1") THEN
    IF (auTOmass < ZERO) const_phys%auTOmass  = TEN**3 / inv_Name
    const_phys%mass_unit = "g.mol-1"
  ELSE
    const_phys%auTOmass  = auTOmass
    const_phys%mass_unit = mass_unit
  END IF

  ! in this subroutine the mass are read in g/mol
  gPmolTOmass          = (TEN**3/inv_Name) / const_phys%auTOmass
  CALL string_uppercase_TO_lowercase(mass_version_loc,lower=.FALSE.) ! conversion in capital letters
  SELECT CASE (mass_version_loc)
  CASE ('NIST2012')
    write(out_unitp,*) 'MASSES, version: ',mass_version_loc
    CALL construct_table_at_NIST2012(const_phys%mendeleev,gPmolTOmass,mass_unit)
  CASE ('NIST2018')
    write(out_unitp,*) 'MASSES, version: ',mass_version_loc
    CALL construct_table_at_NIST2018(const_phys%mendeleev,gPmolTOmass,mass_unit)
  CASE ('HANDBOOK70ED','HANDBOOK')
    write(out_unitp,*) 'MASSES, version: ','HandBook70ed'
    CALL construct_table_at_HandBook70ed(const_phys%mendeleev,gPmolTOmass,mass_unit)
  CASE Default
    write(out_unitp,*) 'MASSES, version: ',mass_version_loc
    CALL construct_table_at_NIST2012(const_phys%mendeleev,gPmolTOmass,mass_unit)
  END SELECT

  !------------------------------------------------------------------

  ! for the automatic energy (E) conversion
  CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(ONE,'au','E'),Work_unit=.TRUE.)
  CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(ONE,'hartree','E'))
  CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(ONE/auTOeV,'eV','E'))
  CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(k/Eh,'°K','E')) ! Kelvin
  CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(ONE/auTOcm_inv,'cm-1','E'))
  CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(ONE/auTOGHz,'GHz','E'))
  CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(ONE/auTOkcalmol_inv,'kcal.mol-1','E'))
  CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(ONE/auTOkJmol_inv,'kJ.mol-1','E'))

  IF (ene_unit == "") THEN
    const_phys%auTOenergy = get_Conv_au_TO_unit('E','cm-1',err_unit=err_unit)
    const_phys%ene_unit   = 'cm-1'
    IF (err_unit /= 0) THEN
      write(out_unitp,*) 'ERROR in ',name_sub
      write(out_unitp,*) ' Problem to get the conversion factor for "cm-1"'
      write(out_unitp,*) '  It should never append.'
      write(out_unitp,*) '  Check the fortran!!'
      write(out_unitp,*) 'List of available units:'
      CALL Write_TabConvRWU_dim1(Tab_conv_FOR_quantity)
      STOP
    END IF
  ELSE
    IF (auTOenergy > ZERO) THEN
      const_phys%auTOenergy = auTOenergy
      const_phys%ene_unit   = ene_unit
      err_unit              = 0
    ELSE
      const_phys%auTOenergy = get_Conv_au_TO_unit('E',ene_unit,err_unit=err_unit)
      const_phys%ene_unit   = trim(adjustl(ene_unit))
      IF (err_unit /= 0) THEN
        IF(MPI_id==0) THEN
          write(out_unitp,*) 'ERROR in ',name_sub
          write(out_unitp,*) ' Problem with "ene_unit" and/or "auTOenergy"'
          write(out_unitp,*) '   energy unit: ',ene_unit
          write(out_unitp,*) '   auTOenergy:  ',auTOenergy
          IF (auTOenergy < ZERO) write(out_unitp,*) ' The value of "auTOenergy" is wrong'
          write(out_unitp,*) 'List of available units:'
        ENDIF
        CALL Write_TabConvRWU_dim1(Tab_conv_FOR_quantity)
        STOP
      END IF
    END IF
  END IF
  ! Now we can add the writing unit for the energy (E).
  CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(ONE/const_phys%auTOenergy,const_phys%ene_unit,'E'),Write_unit=.TRUE.)


  ! for the automatic time (t) conversion
  CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(ONE,'ua','t'),Work_unit=.TRUE.)
  CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(TEN**(-15)/Ta,'fs','t'),Write_unit=.TRUE.) ! fs => atmic unit
  CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(TEN**(-12)/Ta,'ps','t')) ! fs => atmic unit

  ! for the automatic lenght (L) conversion
  CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(ONE,'au','L'),Work_unit=.TRUE.)
  CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(ONE,'bohr','L'))
  CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(TEN**(-10)/a0,'Angs','L'),Write_unit=.TRUE.)
  CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(TEN**(-9)/a0,'nm','L'))
  CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(TEN**(-12)/a0,'pm','L'))


  ! for the automatic angle (angle) conversion
  CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(ONE,'Rad','angle'),Work_unit=.TRUE.)
  CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(pi/180_Rkind,'°','angle'),Write_unit=.TRUE.)

  ! for the automatic electric dipole moment (QL) conversion
  CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(ONE,'au','QL'),Work_unit=.TRUE.)
  CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(convDebyeTOau,'D','QL'),Write_unit=.TRUE.)
  CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(ONE/conv_auTOCm,'C.m','QL'))

  ! for the automatic Electric field (???) conversion
  CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(ONE,'au','Electric field'),Work_unit=.TRUE.)
  CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(ONE/E0,'V.cm-1','Electric field'))

  ! for the automatic Electric field Intensity (???) conversion
  CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(ONE,'au','EF intensity'),Work_unit=.TRUE.)
  CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(ONE/I0,'W.cm-2','EF intensity'))

  ! When no dimension
  CALL ADD_RWU_TO_Tab_conv_FOR_quantity(REAL_WU(ONE,'','No_Dim'),Work_unit=.TRUE.,Write_unit=.TRUE.)

  IF (print_level > 1) THEN
    write(out_unitp,*) 'List of available units:'
    CALL Write_TabConvRWU_dim1(Tab_conv_FOR_quantity)
  END IF

  const_phys%pi         = pi

  const_phys%auTOcm_inv      = auTOcm_inv
  const_phys%auTOeV          = auTOeV
  const_phys%auTOGHz         = auTOGHz
  const_phys%auTOkJmol_inv   = auTOkJmol_inv
  const_phys%auTOkcalmol_inv = auTOkcalmol_inv

  const_phys%inv_Name   = inv_Name

  const_phys%c          = c
  const_phys%mhu0       = mhu0
  const_phys%epsi0      = epsi0
  const_phys%G          = G
  const_phys%h          = h
  const_phys%hb         = hb
  const_phys%e          = e
  const_phys%me         = me
  const_phys%mp         = mp
  const_phys%alpha      = alpha
  const_phys%Na         = Na
  const_phys%R          = R
  const_phys%k          = k

  const_phys%mhu        = mhu
  const_phys%a0         = a0
  const_phys%Eh         = Eh
  const_phys%Ta         = Ta
  const_phys%E0         = E0
  const_phys%I0         = I0


  const_phys%convAif    = convAif
  const_phys%convIDif   = convIDif
  const_phys%convIQif   = convIQif

  const_phys%convDebyeTOau = convDebyeTOau
  const_phys%conv_auTOCm   = conv_auTOCm

!-----------------------------------------------------------------

!-- Write some constantes ------------------------------------
  IF(MPI_id==0) THEN
                         write(out_unitp,*)  ' Writing energy unit : ',RWU_WriteUnit('E',WorkingUnit=.FALSE.)
                         write(out_unitp,*)  ' Working energy unit : ',RWU_WriteUnit('E',WorkingUnit=.TRUE.)
                         write(out_unitp,21) ' auTOenergy      = ',const_phys%auTOenergy
                         write(out_unitp,21) ' auTOcm_inv      = ',const_phys%auTOcm_inv
                         write(out_unitp,21) ' auTOeV          = ',const_phys%auTOeV
                         write(out_unitp,21) ' auTOGHz         = ',const_phys%auTOGHz
                         write(out_unitp,21) ' auTOkcalmol_inv = ',const_phys%auTOkcalmol_inv
                         write(out_unitp,21) ' auTOkJmol_inv   = ',const_phys%auTOkJmol_inv

    IF (print_level > 0) write(out_unitp,*)  ' pi              = ',const_phys%pi
    IF (print_level > 0) write(out_unitp,*)  ' cos(pi)         = ',cos(const_phys%pi)
    IF (print_level > 0) write(out_unitp,11) ' a0 (m-1)        = ',const_phys%a0
    IF (print_level > 0) write(out_unitp,11) ' a0 (Angs)       = ',const_phys%a0*TEN**10
    IF (print_level > 0) write(out_unitp,11) ' Eh (J)          = ',const_phys%Eh
    IF (print_level > 0) write(out_unitp,11) ' Ta (s)          = ',const_phys%Ta
                         write(out_unitp,21) ' Ta (fs)         = ',const_phys%Ta*TEN**15

    IF (print_level > 0) write(out_unitp,21) ' inv_Name        = ',const_phys%inv_Name
    IF (print_level > 0) write(out_unitp,11) ' E0 (V cm-1)     = ',const_phys%E0
    IF (print_level > 0) write(out_unitp,11) ' I0 (W cm-2)     = ',const_phys%I0
    IF (print_level > 0) write(out_unitp,11) ' convAif         = ',const_phys%convAif
    IF (print_level > 0) write(out_unitp,11) ' convIDif        = ',const_phys%convIDif
    IF (print_level > 0) write(out_unitp,11) ' convIQif        = ',const_phys%convIQif
    IF (print_level > 0) write(out_unitp,11) ' convDebyeTOau   = ',convDebyeTOau
  ENDIF

11 format (a,e17.10)
21 format (a,f18.6)

  IF (debug) THEN
    write(out_unitp,*) 'END ',name_sub
  END IF
!-----------------------------------------------------------------
  END SUBROUTINE sub_constantes
  SUBROUTINE sub_ReadNMLconstantes(auTOcm_inv1,ene_unit1,auTOenergy1,   &
                                   inv_Name1,mass_unit1,auTOmass1,      &
                                   time_unit1,                          &
                                   mass_version1,version1)
  IMPLICIT NONE

  real (kind=Rkind),       intent(inout)   :: auTOcm_inv1,auTOenergy1,auTOmass1,inv_Name1
  character(len=Name_len), intent(inout)   :: ene_unit1,mass_unit1,time_unit1
  character(len=Name_len), intent(inout)   :: version1,mass_version1


  integer :: err_read,err_unit


  real (kind=Rkind)        :: inv_Name,auTOmass
  real (kind=Rkind)        :: auTOcm_inv,auTOenergy
  character(len=Name_len)  :: ene_unit
  character(len=Name_len)  :: time_unit
  character(len=Name_len)  :: mass_unit  !  the energy unit (default : au)
  character(len=Name_len)  :: version,mass_version


  NAMELIST /constantes/ auTOcm_inv,inv_Name,mass_unit,auTOmass,     &
                        mass_version,version,EVRT_path,             &
                        ene_unit,auTOenergy,time_unit

!- for debuging --------------------------------------------------
  character (len=*), parameter :: name_sub='sub_ReadNMLconstantes'
  logical, parameter :: debug = .FALSE.
  !logical, parameter :: debug = .TRUE.
!-----------------------------------------------------------------
  IF (debug) THEN
    write(out_unitp,*) 'BEGINNING ',name_sub
  END IF
!-----------------------------------------------------------------

!-- read the namelist ------------------------------------
  version      = version1
  mass_version = mass_version1

  ene_unit     = ene_unit1
  auTOenergy   = auTOenergy1
  auTOcm_inv   = auTOcm_inv1

  mass_unit    = mass_unit1
  auTOmass     = auTOmass1
  inv_Name     = inv_Name1

  time_unit    = time_unit1

  read(in_unitp,constantes,IOSTAT=err_read)
  IF (err_read < 0) THEN
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) ' End-of-file or End-of-record'
    write(out_unitp,*) ' The namelist "constantes" is probably absent'
    write(out_unitp,*) ' check your data!'
    write(out_unitp,*) ' ERROR in ',name_sub
    STOP
  ELSE IF (err_read > 0) THEN
    write(out_unitp,*) ' ERROR in ',name_sub
    write(out_unitp,*) ' Some parameter name of the namelist "constantes" are probaly wrong'
    write(out_unitp,*) ' check your data!'
    write(out_unitp,constantes)
    write(out_unitp,*) ' ERROR in ',name_sub
    STOP
  END IF

  version1      = version
  mass_version1 = mass_version

  ene_unit1     = ene_unit
  auTOenergy1   = auTOenergy
  auTOcm_inv1   = auTOcm_inv

  mass_unit1    = mass_unit
  auTOmass1     = auTOmass
  inv_Name1     = inv_Name

  time_unit1    = time_unit

  IF (debug) THEN
    write(out_unitp,*) 'END ',name_sub
  END IF

END SUBROUTINE sub_ReadNMLconstantes
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!> @brief fundamental physical constants:
!!
!! The CODATA 2014 one in 2018. They can be download from
!! the <a href="http://physics.nist.gov/constants">NIST</a>
  SUBROUTINE constantes_CODATA2018(c,mhu0,G,h,e,me,mp,Na,R)
  IMPLICIT NONE

  !----- physical constants ---------------------------
  real (kind=Rkind) :: c      !< speed of light in vacuum (exact) (in m s-1)
  real (kind=Rkind) :: mhu0   !< Magnetic Constant (exact) (in N A-2)
  real (kind=Rkind) :: G      !< Newtonian constant of gravitation (in m^3 kg-1 s-2)
  real (kind=Rkind) :: h      !< Planck Constant (h et hb) (in J s)
  real (kind=Rkind) :: e      !< Atomic unit of charge (in C)
  real (kind=Rkind) :: me     !< Atomic unit of mass (in kg) (Electron mass)
  real (kind=Rkind) :: mp     !< Proton mass (in kg)
  real (kind=Rkind) :: Na     !< Avogadro constant (in mol-1)
  real (kind=Rkind) :: R      !< Molar gas constant (in J mol−1 K−1)
  real (kind=Rkind) :: k      !< Boltzmann constant (in J K^-1)

  character (len=*), parameter :: version='CODATA 2018'
  !---------------------------------------------------------------------
  write(out_unitp,*) 'PHYSICAL CONSTANTS, version: ',version
  !---------------------------------------------------------------------

  !------ Physical constant of CODATA2018 ---------------------------
  ! http://physics.nist.gov/constants
  ! in 2021: https://physics.nist.gov/cuu/Constants/Table/allascii.txt
  ! in 2021: https://physics.nist.gov/cuu/pdf/wall_2018.pdf

  ! Speed of light in vacuum (exact) (in m s-1)
  ! speed of light in vacuum                                    299792458              (exact)                  m s^-1
  c = 299792458._Rkind
  ! Magnetic Constant (exact) (in N A-2)
  mhu0 = pi*4.e-7_Rkind
  ! Newtonian constant of gravitation (in m^3 kg-1 s-2)
  ! Newtonian constant of gravitation                           6.67430 e-11            0.000 15 e-11            m^3 kg^-1 s^-2
  G = 6.67430e-11_Rkind
  ! Planck Constant (exact) (h et hb) (in J s)
  !Planck constant                                             6.62607015 e-34        (exact)                  J Hz^-1
  h = 6.62607015e-34_Rkind
  ! Elementary Charge -electron- (en C) (exact)
  ! atomic unit of charge                                       1.602176634 e-19       (exact)                  C
  e = 1.602176634e-19_Rkind
  ! Electron mass (en kg)
  ! atomic unit of mass                                         9.1093837015 e-31      0.000 000 0028 e-31      kg
  me = 9.1093837015e-31_Rkind
  ! Proton mass (en kg)
  ! proton mass                                                 1.67262192369 e-27    0.000 000 000 51 e-27    kg
  mp = 1.67262192369e-27_Rkind
  ! Avogadro constant (mol-1) (exact)
  ! Avogadro constant                                           6.02214076 e23         (exact)                  mol^-1
  Na = 6.02214076e23_Rkind
  ! Boltzmann constant                                          1.380649 e-23           (exact)                  J K^-1
  k = 1.380649e-23_Rkind
  ! Molar gas constant R (J mol-1 K-1): R=k*Na
  R = k * Na

END SUBROUTINE constantes_CODATA2018
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!> @brief fundamental physical constants:
!!
!! The CODATA 2014 one in 2018. They can be download from
!! the <a href="http://physics.nist.gov/constants">NIST</a>
  SUBROUTINE constantes_CODATA2014(c,mhu0,G,h,e,me,mp,Na,R)
  IMPLICIT NONE

  !----- physical constants ---------------------------
  real (kind=Rkind) :: c      !< speed of light in vacuum (exact) (in m s-1)
  real (kind=Rkind) :: mhu0   !< Magnetic Constant (exact) (in N A-2)
  real (kind=Rkind) :: G      !< Newtonian constant of gravitation (in m^3 kg-1 s-2)
  real (kind=Rkind) :: h      !< Planck Constant (h et hb) (in J s)
  real (kind=Rkind) :: e      !< Atomic unit of charge (in C)
  real (kind=Rkind) :: me     !< Atomic unit of mass (in kg) (Electron mass)
  real (kind=Rkind) :: mp     !< Proton mass (in kg)
  real (kind=Rkind) :: Na     !< Avogadro constant (in mol-1)
  real (kind=Rkind) :: R      !< Molar gas constant (in J mol−1 K−1)
  real (kind=Rkind) :: k      !< Boltzmann constant (in J K^-1)

  character (len=*), parameter :: version='CODATA 2014'
  !---------------------------------------------------------------------
  write(out_unitp,*) 'PHYSICAL CONSTANTS, version: ',version
  !---------------------------------------------------------------------

  !------ Physical constant of CODATA2014 ---------------------------
  ! http://physics.nist.gov/constants
  ! in 2018: https://ws680.nist.gov/publication/get_pdf.cfm?pub_id=920687
  !  DOI: 10.1103/RevModPhys.88.035009

  ! Speed of light in vacuum (exact) (in m s-1)
  ! speed of light in vacuum                                    299792458              (exact)                  m s^-1
  c = 299792458._Rkind
  ! Magnetic Constant (exact) (in N A-2)
  mhu0 = pi*4.e-7_Rkind
  ! Newtonian constant of gravitation (in m^3 kg-1 s-2)
  ! Newtonian constant of gravitation                           6.67408 e-11            0.000 31 e-11            m^3 kg^-1 s^-2
  G = 6.67408e-11_Rkind
  ! Planck Constant (h et hb) (in J s)
  ! Planck constant                                             6.626070040 e-34       0.000 000 081 e-34       J s
  h = 6.626070040e-34_Rkind
  ! Elementary Charge -electron- (en C)
  ! atomic unit of charge                                       1.6021766208 e-19      0.000 000 0098 e-19      C
  e = 1.6021766208e-19_Rkind
  ! Electron mass (en kg)
  ! Atomic unit of mass                                         9.10938356 e-31        0.000 000 11 e-31        kg
  me = 9.10938356e-31_Rkind
  ! Proton mass (en kg)
  ! Proton mass                                                 1.672621898 e-27       0.000 000 021 e-27       kg
  mp = 1.672621898e-27_Rkind
  ! Avogadro constant (mol-1)
  ! Avogadro constant                                           6.022140857 e23        0.000 000 074 e23        mol^-1
  Na = 6.022140857e23_Rkind
  ! Boltzmann constant                                          1.38064852e-23        0.000 000 79 e-23        J K^-1
  k = 1.38064852e-23_Rkind
  ! Molar gas constant R (J mol-1 K-1): R=k*Na
  R = k * Na

  END SUBROUTINE constantes_CODATA2014
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------

!> @brief fundamental physical constants:
!!
!! The CODATA 2006 one in 2006. There are not available anymore, but they can be download from
!! the <a href="http://physics.nist.gov/cuu/Constants/archive2006.html">NIST</a>
  SUBROUTINE constantes_CODATA2006(c,mhu0,G,h,e,me,mp,Na,R)
  IMPLICIT NONE

  !----- physical constants ---------------------------
  real (kind=Rkind) :: c      !< Speed of light (exact) (in m s-1)
  real (kind=Rkind) :: mhu0   !< Magnetic Constant (exact) (in N A-2)
  real (kind=Rkind) :: G      !< Gravitational Constant (in m^3 kg-1 s-2)
  real (kind=Rkind) :: h      !< Planck Constant (h et hb) (in J s)
  real (kind=Rkind) :: e      !< Elementary Charge -electron- (in C)
  real (kind=Rkind) :: me     !< Electron mass (in kg)
  real (kind=Rkind) :: mp     !< Proton mass (in kg)
  real (kind=Rkind) :: Na     !< Avogadro constant (in mol-1)
  real (kind=Rkind) :: R      !< Molar gas constant (in J mol−1 K−1)

  character (len=*), parameter :: version='CODATA 2006'
  !---------------------------------------------------------------------
  write(out_unitp,*) 'PHYSICAL CONSTANTS, version: ',version
  !---------------------------------------------------------------------

  !------ Physical constant of CODATA2006 ---------------------------
  ! http://www.codata.org/resources/databases/index.html (from NIST now)
  ! http://physics.nist.gov/cuu/Constants/archive2006.html
  ! https://physics.nist.gov/cuu/Constants/RevModPhys_80_000633acc.pdf
  !    DOI: 10.1103/RevModPhys.80.633
  !
  ! Speed of light (exact) (in m s-1)
  c = 299792458._Rkind
  ! Magnetic Constant (exact) (in N A-2)
  mhu0 = pi*4.e-7_Rkind
  ! Gravitational Constant (in m^3 kg-1 s-2)
  G = 6.67428e-11_Rkind
  ! Planck Constant (h et hb) (in J s)
  h = 6.62606896e-34_Rkind
  ! Elementary Charge -electron- (en C)
  e = 1.602176487e-19_Rkind
  ! Electron mass (en kg)
  me = 9.10938215e-31_Rkind
  ! Proton mass (en kg)
  mp = 1.672621637e-27_Rkind
  ! Avogadro constant (mol-1)
  Na = 6.02214179e23_Rkind
  ! constante des gaz parfait R (J mol-1 K-1)
  R = 8.314472_Rkind

  END SUBROUTINE constantes_CODATA2006


  SUBROUTINE constantes_HandBook70ed(c,mhu0,G,h,e,me,mp,Na,R)
  IMPLICIT NONE

  !----- physical constants ---------------------------
  real (kind=Rkind) :: c      !< Speed of light (exact) (in m s-1)
  real (kind=Rkind) :: mhu0   !< Magnetic Constant (exact) (in N A-2)
  real (kind=Rkind) :: G      !< Gravitational Constant (in m^3 kg-1 s-2)
  real (kind=Rkind) :: h      !< Planck Constant (h et hb) (in J s)
  real (kind=Rkind) :: e      !< Elementary Charge -electron- (in C)
  real (kind=Rkind) :: me     !< Electron mass (in kg)
  real (kind=Rkind) :: mp     !< Proton mass (in kg)
  real (kind=Rkind) :: Na     !< Avogadro constant (in mol-1)
  real (kind=Rkind) :: R      !< Molar gas constant (in J mol−1 K−1)

  character (len=*), parameter :: version='HandBook70ed'
  !---------------------------------------------------------------------
  write(out_unitp,*) 'PHYSICAL CONSTANTS, version: ',version
  !---------------------------------------------------------------------

  !------ affectation des constantes avec ---------------------------
  !       Constante du Handbook of Chemisry and Physics (70th ed)
  !       pp F-215 F-219 pour les constantes physiques

  ! vitesse de la lumiere (exacte) (en m s-1)
  c = 299792458._Rkind
  write(out_unitp,*) 'c = ',c
  ! permeabilite du vide (exacte) (en N A-2)
  mhu0 = pi*4.e-7_Rkind
  write(out_unitp,*) 'mhu0 = ',mhu0
  ! permitivite du vide (exacte) (en F m-1)
  G = 6.6725985e-11_Rkind
  ! constante de Planck (h et hb) (en J s)
  h = 6.626075540e-34_Rkind
  write(out_unitp,*) 'h = ',h
  ! charge de l electron (en C)
  e = 1.6021773349e-19_Rkind
  write(out_unitp,*) 'e = ',e
  ! masse de l electron (en kg)
  me = 9.109389754e-31_Rkind

  write(out_unitp,*) 'me = ',me
  ! masse du proton (en kg)
  mp = 1.672623110e-27_Rkind
  ! constante d'Avogadro (mol-1)
  Na = 6.022136736e23_Rkind
  write(out_unitp,*) 'Na = ',Na
  ! constante des gaz parfait R (J mol-1 K-1)
  R = 8.31451070_Rkind
  write(out_unitp,*) 'R = ',R

  END SUBROUTINE constantes_HandBook70ed


  SUBROUTINE constantes_HandBook70ed_2001(c,mhu0,G,h,e,me,mp,Na,R)
  IMPLICIT NONE

  !----- physical constants ---------------------------
  real (kind=Rkind) :: c      !< Speed of light (exact) (in m s-1)
  real (kind=Rkind) :: mhu0   !< Magnetic Constant (exact) (in N A-2)
  real (kind=Rkind) :: G      !< Gravitational Constant (in m^3 kg-1 s-2)
  real (kind=Rkind) :: h      !< Planck Constant (h et hb) (in J s)
  real (kind=Rkind) :: e      !< Elementary Charge -electron- (in C)
  real (kind=Rkind) :: me     !< Electron mass (in kg)
  real (kind=Rkind) :: mp     !< Proton mass (in kg)
  real (kind=Rkind) :: Na     !< Avogadro constant (in mol-1)
  real (kind=Rkind) :: R      !< Molar gas constant (in J mol−1 K−1)

  character (len=*), parameter :: version='constantes_HandBook70ed_2001'
  !---------------------------------------------------------------------
  write(out_unitp,*) 'PHYSICAL CONSTANTS, version: ',version
  !---------------------------------------------------------------------

  !------ affectation des constantes avec ---------------------------
  !       Constante du Handbook of Chemisry and Physics (70th ed)
  !       pp F-215 F-219 pour les constantes physiques
  ! several caonstants were defined in single precision (e, me, mp)

  ! vitesse de la lumiere (exacte) (en m s-1)
  c = 299792458._Rkind
  write(out_unitp,*) 'c = ',c
  ! permeabilite du vide (exacte) (en N A-2)
  mhu0 = pi*4.e-7_Rkind
  write(out_unitp,*) 'mhu0 = ',mhu0
  ! permitivite du vide (exacte) (en F m-1)
  G = 6.6725985e-11_Rkind
  ! constante de Planck (h et hb) (en J s)
  h = 6.626075540e-34_Rkind
  write(out_unitp,*) 'h = ',h
  ! charge de l electron (en C)
  !e = 1.6021773349e-19_Rkind
  e = 1.6021773349e-19
  write(out_unitp,*) 'e = ',e
  ! masse de l electron (en kg)
  !me = 9.109389754e-31_Rkind
  me = 9.109389754e-31

  write(out_unitp,*) 'me = ',me
  ! masse du proton (en kg)
  mp = 1.672623110e-27_Rkind
  mp = 1.672623110e-27

  ! constante d'Avogadro (mol-1)
  Na = 6.022136736e23_Rkind
  write(out_unitp,*) 'Na = ',Na
  ! constante des gaz parfait R (J mol-1 K-1)
  R = 8.31451070_Rkind
  write(out_unitp,*) 'R = ',R

 END SUBROUTINE constantes_HandBook70ed_2001
 END MODULE mod_constant
