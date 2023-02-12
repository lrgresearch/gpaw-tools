---
layout: default
nav_order: 2
title: Input File Keywords
parent: development
---
# Input File Keywords
There are two scripts that use input files in gpaw-tools. One of them is the main DFT script `gpawsolve.py`, and other is MD script `asapsolve.py`. You can find the keyword lists of each script below:
* [**gpawsolve.py Keyword List**](inputfilekeywords.md#gpawsolvepy-keyword-list)
* [**asapsolve.py Keyword List**](inputfilekeywords.md#asapsolvepy-keyword-list)

# `gpawsolve.py` Keyword List

[**General Keywords:**](inputfilekeywords.md#general-keywords) [Mode](inputfilekeywords.md#mode), [Geo_optim](inputfilekeywords.md#geo_optim), [Elastic_calc](inputfilekeywords.md#elastic_calc), [DOS_calc](inputfilekeywords.md#dos_calc), [Band_calc](inputfilekeywords.md#band_calc), [Density_calc](inputfilekeywords.md#density_calc), [Optical_calc](inputfilekeywords.md#optical_calc), [MPI_cores](inputfilekeywords.md#mpi_cores)

[**Geometric Optimization Keywords:**](inputfilekeywords.md#geometric-optimization-keywords) [Optimizer](inputfilekeywords.md#optimizer), [Max_F_tolerance](inputfilekeywords.md#max_f_tolerance), [Max_step](inputfilekeywords.md#max_step), [Alpha](inputfilekeywords.md#alpha), [Damping](inputfilekeywords.md#damping), [Fix_symmetry](inputfilekeywords.md#fix_symmetry), [Relax_cell](inputfilekeywords.md#relax_cell)

[**Electronic Calculations Keywords:**](inputfilekeywords.md#electronic-calculations-keywords) [Cut_off_energy](inputfilekeywords.md#cut_off_energy), [Ground_kpts_density](inputfilekeywords.md#ground_kpts_density), [Ground_kpts_x](inputfilekeywords.md#ground_kpts_x), [Ground_kpts_y](inputfilekeywords.md#ground_kpts_y), [Ground_kpts_z](inputfilekeywords.md#ground_kpts_z), [Ground_gpts_density](inputfilekeywords.md#ground_gpts_density), [Ground_gpts_x](inputfilekeywords.md#ground_gpts_x), [Ground_gpts_y](inputfilekeywords.md#ground_gpts_y), [ground_gpts_z](inputfilekeywords.md#ground_gpts_z), [Gamma](inputfilekeywords.md#gamma), [Band_path](inputfilekeywords.md#band_path), [Band_npoints](inputfilekeywords.md#band_npoints), [Energy_max](inputfilekeywords.md#energy_max), [Setup_params](inputfilekeywords.md#setup_params), [XC_calc](inputfilekeywords.md#xc_calc), [Ground_convergence](inputfilekeywords.md#ground_convergence), [Band_convergence](inputfilekeywords.md#band_convergence), [DOS_convergence](inputfilekeywords.md#dos_convergence), [Occupations](inputfilekeywords.md#occupations), [Mixer_type](inputfilekeywords.md#mixer_type) [DOS_npoints](inputfilekeywords.md#dos_npoints), [DOS_width](inputfilekeywords.md#dos_width), [Spin_calc](inputfilekeywords.md#spin_calc), [Magmom_per_atom](inputfilekeywords.md#magmom_per_atom), [Refine_grid](inputfilekeywords.md#refine_grid)

[**GW Calculations Keywords:**](inputfilekeywords.md#gw-calculations-keywords) [GW_calc_type](inputfilekeywords.md#gw_calc_type), [GW_kpoints_list](inputfilekeywords.md#gw_kpoints_list), [GW_truncation](inputfilekeywords.md#gw_truncation), [GW_cut_off_energy](inputfilekeywords.md#gw_cut_off_energy), [GW_valence_band_no](inputfilekeywords.md#gw_valence_band_no), [GW_conduction_band_no](inputfilekeywords.md#gw_conduction_band_no), [GW_PPA](inputfilekeywords.md#GW_PPA), [GW_q0_correction](inputfilekeywords.md#gw_q0_correction), [GW_nblocks_max](inputfilekeywords.md#gw_nblocks_max), [GW_interpolate_band](inputfilekeywords.md#gw_interpolate_band)

[**Optical Calculations Keywords:**](inputfilekeywords.md#optical-calculations-keywords) [Opt_calc_type](inputfilekeywords.md#opt_calc_type), [Opt_shift_en](inputfilekeywords.md#opt_shift_en), [Opt_BSE_valence](inputfilekeywords.md#opt_bse_valence), [Opt_BSE_conduction](inputfilekeywords.md#opt_bse_conduction), [Opt_BSE_min_en](inputfilekeywords.md#opt_bse_min_en), [Opt_BSE_max_en](inputfilekeywords.md#opt_bse_max_en), [Opt_BSE_num_of_data](inputfilekeywords.md#opt_bse_num_of_data), [Opt_num_of_bands](inputfilekeywords.md#opt_num_of_bands), [Opt_FD_smearing](inputfilekeywords.md#opt_fd_smearing), [Opt_eta](inputfilekeywords.md#opt_eta), [Opt_domega0](inputfilekeywords.md#opt_domega0), [Opt_omega2](inputfilekeywords.md#opt_omega2), [Opt_cut_of_energy](inputfilekeywords.md#opt_cut_of_energy), [Opt_nblocks](inputfilekeywords.md#opt_nblocks)

## General Keywords

---

### Mode
#### Keyword type
String

#### Description
This keyword controls the running mode of the GPAW. Available options are:

* PW
* PW-GW
* EXX
* LCAO
* FD

#### Default
PW

#### Example
Mode = 'PW'

---

### Geo_optim
#### Keyword type
Logical

#### Description
This keyword controls the execution of geometric optimization. Available options are:

* True
* False

User can implement a filter for optimization of supercell and atoms with keyword `whichstrain`. More information about [whichstrain](inputfilekeywords.md#whichstrain).

#### Default
True

#### Example
Geo_optim = False

---

### Elastic_calc
#### Keyword type
Logical

#### Description
This keyword controls the performing of Elastic calculations or not. Available options are:

* True
* False

#### Default
False

#### Example
Elastic_calc = True

---

### DOS_calc
#### Keyword type
Logical

#### Description
This keyword controls the performing of DOS calculations or not. Available options are:

* True
* False

#### Default
False

#### Example
DOS_calc = True

---

### Band_calc
#### Keyword type
Logical

#### Description
This keyword controls the performing of Band calculations or not. Available options are:

* True
* False

#### Default
False

#### Example
Band_calc = False

---

### Density_calc
#### Keyword type
Logical

#### Description
This keyword controls the performing of electron density calculations or not. Available options are:

* True
* False

#### Default
False

#### Example
Density_calc = True

---

### Optical_calc
#### Keyword type
Logical

#### Description
This keyword controls the performing of optical calculations or not. Must be used independently from DOS_calc, Band_calc and Density_calc. Please visit examples directory for the example usage. Available options are:

* True
* False

#### Default
False

#### Example
Optical_calc = False

---

### MPI_cores
#### Keyword type
Integer

#### Description
This keyword controls the number of cores used in calculation. This parameter is not used with `gpawsolve.py`. It is only needed for the `gg.py`.
NOTE: `gg.py` can run `gpawsolve.py` with only `mpirun -np <corenumber>` command. Therefore, for CPUs with hyperthreading support, you can run only the half of your thread number. In the future, `gg.py` will have an option for threads. Please control this variable in future.

#### Default
4

#### Example
MPI_cores = 4

## Geometric Optimization Keywords
### Optimizer
#### Keyword type
String

#### Description
This keyword controls the energy minimization algorithm for the geometry optimization. Available options are:

* LBFGS
* FIRE

#### Default
LBFGS

#### Example
Optimizer = 'FIRE'

---
### Max_F_tolerance
#### Keyword type
Float

#### Description
This keyword controls the maximum force tolerance in BFGS type geometry optimization. Unit is eV/Ang.

#### Default
0.05

#### Example
Max_F_tolerance = 0.05 # eV/Ang

---
### Max_step
#### Keyword type
Float

#### Description
This keyword controls how far a single atom allowed to move. Default is 0.2 Ang.

#### Default
0.2

#### Example
Max_step = 0.2 # Ang

---
### Alpha
#### Keyword type
Float

#### Description
Initial guess for the Hessian (curvature of energy surface)

#### Default
70.0

#### Example
Alpha = 70.0

---
### Damping
#### Keyword type
Float

#### Description
The calculated step is multiplied with this number before added to the positions

#### Default
1.0

#### Example
Damping = 1.0

---
### Fix_symmetry
#### Keyword type
Logical

#### Description
This keyword controls the preserving the spacegroup symmetry during optimisation. Available options are:

* True
* False

#### Default
False

#### Example
Fix_symmetry = True

---

### Relax_cell
#### Keyword type
Python List of Logical values

#### Description
This keyword controls the hich components of strain will be relaxed. There are six independent components indicating the strain are relaxed or not. Here:

* True = relax to zero
* False = fixed

And these six independent components are in order:

* EpsilonX
* EpsilonY
* EpsilonZ
* ShearYZ
* ShearXZ
* ShearXY

**IMPORTANT**: This keyword is only working when `Geo_optim = True` and under PW mode. This feature is not implemented in LCAO mode.


#### Default
[False, False, False, False, False, False]

#### Example
Relax_cell = [True, True, False, False, False, False] #For a x-y 2D nanosheet only first 2 component will be true

---
## Electronic Calculations Keywords
### Cut_off_energy
#### Keyword type
Integer

#### Description
This keyword controls the plane wave cut off energy value. Unit is eV. Can be used in PW mode.

#### Default
340 eV

#### Example
Cut_off_energy = 500 # eV

---

### Ground_kpts_density
#### Keyword type
Float

#### Description
This keyword controls kpoint density. It is deactivated normally. Monkhorst-Pack mesh is used with `Ground_kpts_x`, `Ground_kpts_y` and `Ground_kpts_z` variables. If `Ground_kpts_density` is included in an input file, the `Ground_kpts_x`, `Ground_kpts_y` and `Ground_kpts_z` variables will be ignored automatically. Unit is pts per Å^-1.

#### Default
Not used in default.

#### Example
Ground_kpts_density = 2.5     # pts per Å^-1

---

### Ground_kpts_x
#### Keyword type
Integer

#### Description
This keyword controls the number of kpoints in x direction. If `Ground_kpts_density` is included in an input file, the `Ground_kpts_x` variable will be ignored automatically. Unit is number of points.

#### Default
5

#### Example
Ground_kpts_x = 5

---

### Ground_kpts_y
#### Keyword type
Integer

#### Description
This keyword controls the number of kpoints in y direction. If `Ground_kpts_density` is included in an input file, the `Ground_kpts_y` variable will be ignored automatically. Unit is number of points.

#### Default
5

#### Example
Ground_kpts_y = 5

---

### Ground_kpts_z
#### Keyword type
Integer

#### Description
This keyword controls the number of kpoints in z direction. If `Ground_kpts_density` is included in an input file, the `Ground_kpts_z` variable will be ignored automatically. Unit is number of points.

#### Default
5

#### Example
Ground_kpts_z = 5

---

### Ground_gpts_density
#### Keyword type
Float

#### Description
This keyword controls gpoint density (LCAO only). If `Ground_gpts_density` is included in an input file, the `Ground_gpts_x`, `Ground_gpts_y` and `Ground_gpts_z` variables will be ignored automatically. Unit is pts per Å.

#### Default
0.2.

#### Example
Ground_gpts_density = 0.2     # pts per Å^-1

---

### Ground_gpts_x
#### Keyword type
Integer

#### Description
This keyword controls the number of grid points in x direction (LCAO only). If `Ground_gpts_density` is included in an input file, the `Ground_gpts_x` variable will be ignored automatically. Unit is number of points.

#### Default
8

#### Example
Ground_gpts_x = 8

---

### Ground_gpts_y
#### Keyword type
Integer

#### Description
This keyword controls the number of grid points in y direction (LCAO only). If `Ground_gpts_density` is included in an input file, the `Ground_gpts_y` variable will be ignored automatically. Unit is number of points.

#### Default
8

#### Example
Ground_gpts_y = 8

---

### Ground_gpts_z
#### Keyword type
Integer

#### Description
This keyword controls the number of grid points in z direction (LCAO only). If `Ground_gpts_density` is included in an input file, the `Ground_gpts_z` variable will be ignored automatically. Unit is number of points.

#### Default
8

#### Example
Ground_gpts_z = 8

---

### Gamma
#### Keyword type
Logical

#### Description
This keyword controls the inclusion of Gamma point in band calculations. Available options are:

* True
* False

#### Default
True

#### Example
Gamma = False

---

### Band_path
#### Keyword type
String

#### Description
This keyword controls the path of high-symmetry points in band structure diagram. Use 'G' for Gamma point. 

#### Default
'LGL'

#### Example
Band_path = 'GMKG'

---

### Band_npoints
#### Keyword type
Integer

#### Description
This keyword controls the number of points between the first and the last high symmetry points. 

#### Default
60

#### Example
Band_npoints = 50

---

### Energy_max
#### Keyword type
Integer

#### Description
This keyword controls the maximum energy value when the software is used with -d (draw) argument. number of points between the first and the last high symmetry points. Unit is eV.

#### Default
15

#### Example
Energy_max = 10 # eV

---

### Setup_params
#### Keyword type
Python dictionary

#### Description
This keyword controls the implementation of setup parameter on the related orbitals of related elements. For none use {}. More information can be found [here](https://wiki.fysik.dtu.dk/gpaw/documentation/basic.html#manual-setups).Unit is eV.

#### Default
{}

#### Example
Setup_params = {'N': ':p,6.0'} # eV

---

### XC_calc
#### Keyword type
String

#### Description
This keyword controls the which exchange-correlation functional is used in the calculation.Available options are:

* LDA
* PBE
* GLLBSC (-)
* revPBE
* RPBE
* HSE03 (-)
* HSE06 (-)
* B3LYP (can be used only with PW-EXX)
* PBE0  (can be used only with PW-EXX)

(-): whichstrain keyword must be [False, False, False, False, False, False]

Because GPAW is using libxc, there are many exchange-correlation functionals available to use. However, the above functionals are used and tested successfully with gpaw-tools. Please try other possible functionals, make us know, send us input files.

#### Default
LDA

#### Example
XC_calc = 'PBE'

---
### Ground_convergence
#### Keyword type
Python dictionary

#### Description
This keyword controls the convergence parameters for the ground-state calculations. For default use {}.

#### Default
{'energy': 0.0005,  # eV / electron
 'density': 1.0e-4,  # electrons / electron
 'eigenstates': 4.0e-8,  # eV^2 / electron
 'forces': np.inf,
 'bands': None,
 'maximum iterations': None}

#### Example
Ground_convergence = {'energy': 0.005} # eV

---
### Band_convergence
#### Keyword type
Python dictionary

#### Description
This keyword controls the convergence parameters for the band calculations.

#### Default
{'bands':8} 

#### Example
Band_convergence = {'bands':8, 'eigenstates': 1.0e-8} 

---

### DOS_convergence
#### Keyword type
Python dictionary

#### Description
This keyword controls the convergence parameters for the DOS calculations.

#### Default
{} 

#### Example
DOS_convergence = {'maximum iterations': 100} 

---
### Occupations
#### Keyword type
Python dictionary

#### Description
This keyword controls the smearing of the occupation numbers. You can use 4 types:

* improved-tetrahedron-method
* tetrahedron-method
* fermi-dirac
* marzari-vanderbilt

#### Default
{'name': 'fermi-dirac', 'width': 0.05}

#### Example
Occupations = {'name': 'marzari-vanderbilt', 'width': 0.2}

---

### Mixer_type
#### Keyword type
Python import

#### Description
This keyword controls a number of density mixing posibilities. Detailed information can be found on [GPAW's webpage about density mixing](https://wiki.fysik.dtu.dk/gpaw/documentation/densitymix/densitymix.html).

You can use
* Mixer()
* MixerSum()
* MixerDif()

You need to import these modules in the input file like:

    from gpaw import Mixer
    
or 

    from gpaw import MixerSum
    
or

    from gpaw import MixerDif
    
The values of mixer modules corresponds (beta, nmaxold, weight). If you have convergence problems, you can try  (0.02, 5, 100) and (0.05, 5, 50)

#### Default
MixerSum(0.1,3,50)

#### Example
Mixer_type = Mixer(0.02, 5, 100)

---

### DOS_npoints
#### Keyword type
Integer

#### Description
This keyword controls the number of data points for DOS data:

#### Default
501

#### Example
DOS_npoints = 1001

---

### DOS_width
#### Keyword type
Float

#### Description
This keyword controls the width of Gaussian smearing in DOS calculation. Use 0.0 for linear tetrahedron interpolation.

#### Default
0.1

#### Example
DOS_width = 0.0 #Using tetrahedron interpolation

---

### Spin_calc
#### Keyword type
Logical

#### Description
This keyword controls the inclusion of spin based calculations. Please do not forget to set `Magmom_per_atom` variable. Available options are:

* True
* False

Because GPAW is using libxc, there are many exchange-correlation functionals available to use. However, the above functionals are used and tested successfully with gpaw-tools. Please try other possible functionals, make us know, send us input files.

#### Default
False

#### Example
Spin_calc = True

---

### Magmom_per_atom
#### Keyword type
Float

#### Description
This keyword controls the value of magnetic moment of each atom. Please do not forget to set `Spin_calc` variable to `True`. Unit is μB.

#### Default
1.0

#### Example
Magmom_per_atom = 1.0

## GW Calculations Keywords
### GW_calc_type
#### Keyword type
String

#### Description
This keyword controls the type GW calculation. Available options are:

* GW0
* G0W0

#### Default
GW0

#### Example
GW_calc_type = 'GW0'

---

### GW_kpoints_list
#### Keyword type
NumPy Array

#### Description
This keyword represents the kpoint coordinates for the GW calculation.

#### Default
np.array([[0.0, 0.0, 0.0], [1 / 3, 1 / 3, 0], [0.0, 0.0, 0.0]])

#### Example
GW_kpoints_list = np.array([[0.0, 0.0, 0.0], [1 / 3, 1 / 3, 0], [0.0, 0.0, 0.0]])

---

### GW_truncation 
#### Keyword type
NumPy Array

#### Description
This keyword controls the truncation of Coulomb potential for the GW calculations. Available options are:

* None
* 2D
* 1D
* 0D
* wigner-seitz

#### Default
None

#### Example
GW_truncation = '2D'

---

### GW_cut_off_energy
#### Keyword type
Integer

#### Description
This keyword controls the cut off energy value for the GW calculations. Unit is eV.

#### Default
50 eV

#### Example
GW_cut_off_energy = 50

---

### GW_valence_band_no
#### Keyword type
Integer

#### Description
This keyword controls the number of the band for the valence band for GW calculations.

#### Default
8 (Default value is not a general value. Please find correct band for your calculation.)

#### Example
GW_valence_band_no = 8

---

### GW_conduction_band_no
#### Keyword type
Integer

#### Description
This keyword controls the number of the band for the conduction band for GW calculations.

#### Default
18 (Default value is not a general value. Please find correct band for your calculation.)

#### Example
GW_conduction_band_no = 18

---

### GW_PPA
#### Keyword type
Logical

#### Description
This keyword controls the usage of Plasmon Pole Approximation (PPA) for GW calculations.

#### Default
True

#### Example
GW_PPA = True

---

### GW_q0_correction
#### Keyword type
Logical

#### Description
This keyword controls the usage of analytic correction to the q=0 contribution applicable to 2D systems.

#### Default
True

#### Example
GW_q0_correction = True

---

### GW_nblocks_max
#### Keyword type
Logical

#### Description
This keyword controls the behaviour of cutting chi0 into as many blocks to reduce memory requirement as much as possible.

#### Default
True

#### Example
GW_nblocks_max = True

---

### GW_interpolate_band
#### Keyword type
Logical

#### Description
This keyword controls the behaviour of drawing the band with interpolating the values between the points.

#### Default
True

#### Example
GW_interpolate_band = True

## Optical Calculations Keywords
### Opt_calc_type
#### Keyword type
String

#### Description
This keyword controls the optical calculation type: random phase approximation (RPA) or Bethe-Salpeter Equation (BSE).

#### Default
BSE

#### Example
Opt_calc_type = 'BSE'

---

### Opt_shift_en
#### Keyword type
Float

#### Description
This keyword add a shifting to energy values. Unit is eV. Works on BSE calculations only!

#### Default
0.0

#### Example
Opt_shift_en = 1.0 #eV

---

### Opt_BSE_valence
#### Keyword type
Sequence of integers

#### Description
This keyword shows the valence bands that will be used in BSE calculation.

#### Default
range(0,3)

#### Example
Opt_BSE_valence = range(120,124)

---

### Opt_BSE_conduction
#### Keyword type
Sequence of integers

#### Description
This keyword shows the conduction bands that will be used in BSE calculation.

#### Default
range(4,7)

#### Example
Opt_BSE_conduction = range(124,128)

---

### Opt_BSE_min_en
#### Keyword type
Float

#### Description
This keyword shows the starting energy value of result data that will be used in BSE calculation.

#### Default
0.0

#### Example
Opt_BSE_min_en = 0.0

---

### Opt_BSE_max_en
#### Keyword type
Float

#### Description
This keyword shows the ending energy value of result data that will be used in BSE calculation.

#### Default
20.0

#### Example
Opt_BSE_max_en = 10.0

---

### Opt_BSE_num_of_data
#### Keyword type
Integer

#### Description
This keyword shows the number of data points in BSE calculation.

#### Default
1001

#### Example
Opt_BSE_num_of_data = 401

---

### Opt_num_of_bands
#### Keyword type
Integer

#### Description
This keyword controls the number of bands used in optical calculations.

#### Default
16

#### Example
Opt_num_of_bands = 8

---

### Opt_FD_smearing
#### Keyword type
Float

#### Description
This keyword controls the Fermi Dirac smearing for optical calculations.

#### Default
0.05

#### Example
Opt_FD_smearing = 0.02

---

### Opt_eta
#### Keyword type
Float

#### Description
This keyword controls the broadening parameter -eta- used in dielectric function calculations.

#### Default
0.2

#### Example
Opt_eta = 0.1

---

### Opt_domega0
#### Keyword type
Float

#### Description
This keyword controls the  [Δω0 parameter](https://wiki.fysik.dtu.dk/gpaw/tutorialsexercises/opticalresponse/dielectric_response/dielectric_response.html) for non-linear frequency grid used in dielectric function calculations. Unit is eV.

#### Default
0.1 eV

#### Example
Opt_domega0 = 0.05 # eV

---

### Opt_omega2
#### Keyword type
Float

#### Description
This keyword controls the  [ω2 parameter](https://wiki.fysik.dtu.dk/gpaw/tutorialsexercises/opticalresponse/dielectric_response/dielectric_response.html) for non-linear frequency grid used in dielectric function calculations. Unit is eV.

#### Default
10.0 eV

#### Example
Opt_omega2 = 2.0 # eV

---

### Opt_cut_of_energy
#### Keyword type
Float

#### Description
This keyword controls the planewave energy cutoff in dielectric function calculations. Determines the size of dielectric matrix. Unit is eV.

#### Default
10.0 eV

#### Example
Opt_cut_of_energy = 20.0 # eV

---

### Opt_nblocks
#### Keyword type
Integer

#### Description
This keyword controls the split matrices in nblocks blocks and distribute them G-vectors or frequencies over processes.

#### Default
4

#### Example
Opt_nblocks = 4

# `asapsolve.py` Keyword List

[**MD Keywords:**](inputfilekeywords.md#md-keywords) [OpenKIM_potential](inputfilekeywords.md#openkim_potential), [Temperature](inputfilekeywords.md#temperature), [Time](inputfilekeywords.md#time), [Friction](inputfilekeywords.md#friction), [Scaled](inputfilekeywords.md#scaled), [Manual_PBC](inputfilekeywords.md#manual_pbc), [pbc_constraints](inputfilekeywords.md#pbc_constraints), [Solve_double_element_problem](inputfilekeywords.md#solve_double_element_problem)

## MD Keywords

---

### OpenKIM_potential
#### Keyword type
String

#### Description
This keyword controls the interatomic potential used in the calculation.

#### Default
'LJ_ElliottAkerson_2015_Universal__MO_959249795837_003'

#### Example
OpenKIM_potential = 'LJ_ElliottAkerson_2015_Universal__MO_959249795837_003'

---

### Temperature
#### Keyword type
Integer

#### Description
This keyword controls the temperature used in the calculation. Unit is Kelvin.

#### Default
1

#### Example
Temperature = 300 #K

---

### Time
#### Keyword type
Float

#### Description
This keyword controls the timestep used in the calculation. Unit is femtosecond

#### Default
5

#### Example
Time = 10 #fs

---

### Friction
#### Keyword type
Float

#### Description
This keyword controls the friction used in the calculation.

#### Default
0.05

#### Example
Friction = 0.1

---

### Scaled
#### Keyword type
Boolean

#### Description
This keyword controls the usage of scaled or cartesian coordinates in the calculation.

#### Default
False

#### Example
Scaled = True

---

### Manual_PBC
#### Keyword type
Boolean

#### Description
This keyword controls the usage manual constraint axis in the calculation. If it used as True, [PBC_constraints](inputfilekeywords.md#pbc_constraints) keyword must be used.

#### Default
False

#### Example
Manual_PBC = True

---

### PBC_constraints
#### Keyword type
Python List of Logical values

#### Description
This keyword controls the which components of axes will be constrained. Here:

* True = constrained
* False = not constrained

And these 3 independent components are in order:

* X
* Y
* Z

**IMPORTANT**: This keyword is only working when `Manual_PBC = True`

#### Default
[True, True, False]

#### Example
PBC_constraints = [True, False, False]

---

### Solve_double_element_problem
#### Keyword type
Boolean

#### Description
This keyword is used for a possible problem. If you have double number of elements in your final file, please use this keyword as `True`.

#### Default
True

#### Example
Solve_double_element_problem = False

---
