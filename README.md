# ElVibRot-TnumTana

General quantum dynamics code using curvilinear coordinates:
1. Vibrational levels, intensities for floppy molecular systems
2. Wave-packet propagation with or without time dependant Hamiltonian
3. Quantum gate and optimal control.
4. Optimization with the given set of curvilinear coordiantes

## 1) Introduction

Originalities of this code:
 * No built-in limitation in terms of number of degrees of freedom.
 * the use of a numerical but exact kinetic energy operator with Tnum (Automatic differentiation), which enables a large flexibility in the choice of the curvilinear coordinates.
 * the use of Smolyak scheme, which enables to avoid a direct-product basis set and grids.

## 2) Input file
 The input file has four mains sections:
* *SYSTEM* and *CONSTANTS*, which define general parameters for parallelization, printing levels, energy unit, physical constants...
* *COORDINATES*, which defines the curvilinear coordinates, the coordinates transformations and some aspects of the physical models (constraints....). This section is part of Tnum.
* *OPERATORS* and *BASIS SETS*, which define parameters of scalar operators (potential, dipole moments...) and the active and inactive basis set (contracted).
* *ANALYSIS*, which defines parameters for time dependent (including optimal control) or independent calculations, intensities.


### 2.1) &system namelist

This namelist can be used without parameter. The default values should work for most of the situations. However, the following parameters can be used in this namelist:

* **printlevel** (default **-1**): The value **-1** gives the minimal output. The value **0** gives more output and **1** gives even more.


* **RMatFormat** (default **"f18.10"**): This parameter controls the format to print a REAL matrix.
* **CMatFormat** (default **"f15.7"**): This parameter controls the format to print a COMPLEX matrix.
* **EneFormat** (default **"f18.10"**): This parameter controls the format to print the energy.


* **EVR** (default **T**): The value **T** enables to use **ElVibRot**. The value **F** disables the **ElVibRot** calculation.
* **intensity_only** (default **F**): The value **T** enables to calculate only the intensities with ElVibRot (for time independent calculation) with the restart file (*restart.int*). It is useful, when one want to change temperature.
* **analysis_only** (default **F**): The value **T** enables to analyse the wave functions from the "file_spectralWP" file obtained from a previous ElVibRot calculation with the same bassis set.
* **optimization** (default **0**):  The value **0** disables the optimization. The value **1** enables the optimization (geometry or other ElVibRot parameters).
* **nDGrid** (default **F**): The value **T** enables to generate a multidimensional grid (for the nDFit) using a set of curvilinear coordinates
* **nDFit** (default **F**): The value **T** enables to fit the multidimensional grid to an analytical form. Then it can be used with ElVibRot.
* **cart** (default **F**): The value **T** enables to use **cart**. The value **F** disables the **cart** calculation.

* **GridTOBasis_test** (default **F**): The value **T** enables to make a grid-to-basis and a basis-to-grid transformations. Developpers only.
* **OpPsi_test** (default **F**): The value **T** enables to calculate one Hamiltonian action. Developpers only.
* **main_test** (default **F**): The value **T** enables to use some testing units. Developpers only.



* **Popenmp** (default **T**): The value **T** enables OPENMP parallelization in some parts of the code. The value **F** disables OpenMP parallelization. When its value is **T**, the following parameters enable to control more precisely the OpenMP parallelization.

* **PGrid_omp** (default **1**): The value **1** enables the parallelization of the grid calculation of operators. The value **0** disables this parallelization. Efficient parallelization.
* **PGrid_maxth** (default **$OMP_NUM_THREADS value**): The maximum number of threads when **PGrid_omp** > 0.

* **PSG4_omp** (default **1**): The value **1** enables parallelization when Smolyak scheme (the 4th one) is used. The value **0** disables this parallelization. Relatively efficient parallelization.
* **PSG4_maxth** (default **$OMP_NUM_THREADS value**): The maximum number of threads when **PSG4_omp** > 0.

* **PMatOp_omp** (default **0**): The value **1** enables the parallelization of the matrix calculation of an operator. The value **0** disables this parallelization. Do not use, it is not efficient.
* **PMatOp_maxth** (default **$OMP_NUM_THREADS value**): The maximum number of threads when **PMatOp_omp** > 0.

* **POpPsi_omp** (default **0**): The value **1** enables the parallelization of the action of an operator on a wave function. The value **0** disables this parallelization. Do not use, it is not efficient.
* **POpPsi_maxth** (default **$OMP_NUM_THREADS value**): The maximum number of threads when **POpPsi_omp** > 0.

* **PBasisTOGrid_omp** (default **0**): The value**1** enables the parallelization of the grid <-> basis transformation. The value **0** disables this parallelization. Do not use, it is not efficient.
* **PBasisTOGrid_maxth** (default **$OMP_NUM_THREADS value**): The maximum number of threads when **POpPsi_omp** > 0.


* **EVRT_path** (default **directory of compilation**): it enables to read the isotopic masses from an internal file (this parameter can be change also in the namelist "&constantes") and other files used in ElVibRot (Hermite cubature and Lebedev grids).

* **File_path**: it enables to define the path were most of the files will be save.


### 2.2) &constantes namelist (physical constants and units)

This namelist can be used without parameter, it enables to change energy unit and selects different version of physical constants or to redefine some of them.

In ElVibRot, most of the quantities (energy, time, length, angle...) can be read with an unit. Then, they are converted to the working unit (atomic unit).

The following parameters can be used in this namelist:

* **version** (default **"CODATA2006"**): The value **"CODATA2006"** enables to use physical constants defined in codata2006. The other possibilities are: **"CODATA2014"**, **"HANDBOOK70ED"** (70ed of the Physical chemistry Handbook) and **"PUBLI2001"** (modification of 70ed of the Physical chemistry Handbook). The default is used when the "version" is not correct.

* **mass_version** (default **"NIST2012"**): The value **"NIST2012"** enables to use masses published by the NIST in 2012. The other possibilities are: **"NIST2018"**, **"HANDBOOK70ED"** (70ed of the Physical chemistry Handbook). The default is used when the "mass_version" is not correct.


* **ene_unit** (default **"cm-1"**): This value stores the name of the printed energy unit and it enables to select some energy conversion (see **auTOenergy**). If the value of **ene_unit** is different from the ones in following list, a personal unit is used and **ene_unit** is just the name of this new unit. The conversion from atomic unit to the personal the unit is selected with **auTOenergy** parameter. The possible values are of **ene_unit** are:

    * **"cm-1"**: wave number unit (default).
    * **"au"** or **"hartree"**: atomic unit or Hartree unit (the working energy unit).
    * **"ev"** or **"eV"**: electron volt unit.
    * **"GHz"**: energy expressed as GHz.
    * **"°K"**: energy expressed as Kelvin.
    * **"kcal.mol-1"**
    * **"kJ.mol-1"**

* **auTOenergy**: This parameter enables to change the printed energy unit conversion.


* **EVRT_path** (default **directory of compilation**): it enables to read the isotopic masses from an internal file.
If this parameter is changed  in the namelist "&system", the value from the &constantes namelist will be used.


The following parameters can be use to modify some physical constants (to reproduce calculations with wrong ones).

* **auTOcm_inv** : This parameter enables to modify energy conversion factor to **"cm-1"**. It has an effect only if **ene_unit="cm-1"**.
* **inv_Name** : This parameter enables to modify mass conversion factor (au <=> g.mol-1).

### 2.3) &variables namelist (Tnum)

This namelist is part of **Tnum** (or **Tana**) and it enables to define the curvilinear coordinates used in the dynamics. It enables to set reduced dimensionality models (rigid or flexible constraints, some part of adiabatic approximation...).
See the **Tnum** manual.

### 2.4) &minimun namelist (Tnum)

The parameters of this namelist define two mains features:

#### a) The reference geometry (Tnum).

This geometry is used to define the rigid constraints in **Tnum**. The relevant parameters are the following:

* **read_Qdyn0** (default **T**): When the value is **T**, the reference geometry is defined with the dynamical coordinates (see **Tnum**). Remark, in previous version, this parameter was defined as **read_Qsym0**.
* **read_Qact0** (default **F**): When the value is **T**, the reference geometry is defined with the active coordinates (see **Tnum**).
* **read_xyz0** (default **F**): When the value is **T**, the reference geometry is defined with the Cartesian coordinates. Be careful, the transformation from Cartesian to dynamical coordinates is not always possible with the Tnum. See **read_xyz0_with_dummy**.
* **read_itQ0transfo** (default **-1**): The value defines the coordinate transformation number for which the reference geometry is read. Any numbers from **0** (for Cartesian coordinates), **nb_Qtransfo** (for active coordinates) are possible (see **Tnum** for more details).
* **read_xyz0_with_dummy** (default **T**): When the value is **T**, the Cartesian coordinates are read with the dummy atoms. If the value is **F**, the transformation from Cartesian coordinates to the active ones are not possible. **Be carrefull, the dummy atoms are at the end in reverse order.**
* **unit** (default **"au"**): When the value is **"au"** (for atomic unit), the coordinate units are "bohr" and "radian". Otherwise, the units are the "Angstrom" and the "degree". This is relevant for Cartesian coordinates (**read_xyz0=T**) and when the unit are not read with the coordinates values.


Old parameters (not used anymore)
* **read_nameQ** (default **F**): When the value is **T**, the coordinate names are read.
* **read_Qsym0**: see **read_Qdyn0**

#### b) Options for the scalar operators (part of **ElVibRot**).

The relevant parameters are the following:

* **nb_elec** (default **1**): This parameter defines the number of diabatic electronic surface used in the dynamic. The potential PES) is defined as a matrix (nb_elec x nb_elec).
* **pot_cplx** (default **F**): When the value is **T**, the PES are complex and the imaginary part is defined as a matrix (nb_elec x nb_elec).
* **pot0** (default **0.**): This value is the energy reference of the PES. This value is not automatically defined from the PES grid, therefore it has to be set-up. Otherwise, the printed energy with respect to the **pot0** value might be not printable. Nevertheless, the quantum dynamics calculation will be correct.
* **nb_scalar_Op** (default **0**): This parameter defines the number of scalar operators (without the PES) such as the dipole moment. It is used when intensities are calculated. Remark, when  **nb_scalar_Op=0** and when **intensity=t** is set up in the **&analysis** namelist, the **nb_scalar_Op** is modified to **3**. All these operators are defined as matrices (nb_elec x nb_elec x nb_scalar_Op).

By default, the operators (PES and the scalar ones) are obtained in *calcN_op* fortran subroutine from the *sub_system.f* or *sub_system.f90* files (default *sub_system.f*. This default can be changed in the *makfile*). However, they are several other posibilities:

* **OnTheFly** (default **F**): When the value is **T**, the scalar operators and the potential energy surfaces are obtained "on-the-fly" with quantum chemistry codes (only Gaussian and Gamess-US are defined). Only the 3 components of the dipole moment can be calculated (**nb_scalar_Op=3**).
* **QMLib**  (default **F**): When the value is **T**, the potential energy surfaces are obtained from the Quantum Model library (https://github.com/lauvergn/QuantumModelLib). When the metric tensor, **G**, is assumed to be constant (**Gcte=t** in &variables namelist), the metric tensor is used also from the **QMLib**.
* **nDfit_Op**  (default **F**): When the value is **T**, the potential energy surfaces are obtained from the the nDfit procedure. The data needed for the fits are stored in files with names starting with "BaseName_nDfit_file".
* **BaseName_nDfit_file**: base name of the files used for the nDfit procedure.

In link with the Tnum coordinates, one can select with the different sets of coordinates which will be used for the operator calculations as follows:

* **pot_act** (default **T**): When the value is **T**, the scalar operators and the potential energy surfaces are evaluated with active coordinates (in active order).
* **pot_cart** (default **F**): When the value is **T**, the scalar operators and the potential energy surfaces are evaluated with Cartesian coordinates.
* **pot_itQtransfo** (default **-1**): The value defines the coordinate transformation number for which the scalar operators and the potential energy surfaces are evaluated. Any numbers from **0** (for Cartesian coordinates), **nb_Qtransfo** (for active coordinates) are possible (see **Tnum** for more details).
Remark: when **pot_act=F**, **pot_cart=F** and **pot_itQtransfo=-1**, the scalar operators and the potential energy surfaces are evaluated with dynamical coordinates (not in active order).


Other parameters:

* **HarD** (default **T**): This parameter is relevant only when Harmonic adiabatic separation is used (HADA or cHAC, coordinate types, 21 or 22 in Tnum). When the value is **T**, the scalar operators and the potential energy surfaces are evaluated only along the "active" coordinates (coordinate type, 1) and a linear and a quadratic contributions are added with *d0d1d2_g* *d0d1d2_h* subroutines (*sub_system.f* or *sub_system.f90* files). When the value is **F**, the scalar operators and the potential energy surfaces are evaluated only along all coordinates (coordinate types, 1, 21 or 22). It enables to recover exact calculation with coupled Harmonic Adiabatic Channels (cHAC).

* **deriv_WITH_FiniteDiff** (default **F**): when the value is **T**, it forces the numerical hessian calculation with finite differences.

* **opt** (default **F**): when the value is **T** and when the **optimization > 0** in the *&system* namelist, the active coordinates could be optimized (see Optimization procedure).


## 3) Installation

The installation is simple. However, we do not have an fully automatic procedure (like configure ...).
However, the program uses some fortran 2003 features. Therefore, the compilers gfortran or ifort or pgf90 can be used (probably others as well).

You have to select or define your compiler in the "makefile" (the default is gfortran).
Then, you have to compile the program with the unix command "make".

Currently, the program can be compiled on:
* linux platform with
-   gfortran (version= 6.3, 8.3, 9.1)
-   ifort    (version= 18.0.3) Pb with 17.0.1
-   pgf90    (version= 17.10-0) Pb with 16.4-0
* osx platform with gfortran


There are severals options with can be modified in the makefile:
* F90: the compiler (see above)
* OPT: the compiler optimization. 0 => no optimization, 1 => optimization
* OMP: compilation with (OMP=1) or without (OMP=0) OpenMP.
* INT: This enables to change the integer king default (INT=4) during the compilation to a "long integer" INT=8. This is usefull for large calculations (Smolyak).
* ARPACK: this enables the use of ARPACK diagonalisation (ARPACK=1). It needs Lapack as well. 0 => without ARPACK (default), 1 => with ARPACK.
* LAPACK: this enables the use of BLAS and LAPACK libraries. 0 => without LAPACK, 1 => with LAPACK (default)
* QML: Quantum Model Lib.  0 => without QML (default), 1 => with QML.
* extf: it enables to change the "sub_system" file extention. Possible values: f or f90. Alternativaly, one can use an external unix variable: EXTFextern.

To check that the program has be compile correctly, you can run some tests from the directory: **exa_hcn-dist**. Each input data starts with dat....
To run the **dat_Davidson** example, just the type the **dat_Davidson** command. The output will be in the **res** file.
