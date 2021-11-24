---
layout: default
navigation_weight: 6
title: Input File Keywords
---

# Input File Keywords

[**General Keywords:**](inputfilekeywords.md#general-keywords) [Mode](inputfilekeywords.md#mode), [DOS_calc](inputfilekeywords.md#dos_calc), [Band_calc](inputfilekeywords.md#band_calc), [Density_calc](inputfilekeywords.md#density_calc), [Optical_calc](inputfilekeywords.md#optical_calc), [whichstrain](inputfilekeywords.md#whichstrain), [MPIcores](inputfilekeywords.md#mpicores)

[**Electronic Calculations Keywords:**](inputfilekeywords.md#electronic-calculations-keywords) [fmaxval](inputfilekeywords.md#fmaxval), [cut_off_energy](inputfilekeywords.md#cut_off_energy), [kpts_density](inputfilekeywords.md#kpts_density), [kpts_x](inputfilekeywords.md#kpts_x), [kpts_y](inputfilekeywords.md#kpts_y), [kpts_z](inputfilekeywords.md#kpts_z), [Gamma](inputfilekeywords.md#gamma), [band_path](inputfilekeywords.md#band_path), [band_npoints](inputfilekeywords.md#band_npoints), [energy_max](inputfilekeywords.md#energy_max), [Hubbard](inputfilekeywords.md#hubbard), [XC_calc](inputfilekeywords.md#xc_calc), [Spin_calc](inputfilekeywords.md#spin_calc), [Magmom_per_atom](inputfilekeywords.md#magmom_per_atom), [gridref](inputfilekeywords.md#gridref)

[**GW Calculations Keywords:**](inputfilekeywords.md#gw-calculations-keywords) [GWtype](inputfilekeywords.md#gwtype), [GWkpoints](inputfilekeywords.md#gwkpoints), [GWtruncation](inputfilekeywords.md#gwtruncation), [GWcut_off_energy](inputfilekeywords.md#gwcut_off_energy), [GWbandVB](inputfilekeywords.md#gwbandvb), [GWbandCB](inputfilekeywords.md#gwbandcb), [GWppa](inputfilekeywords.md#gwppa), [GWq0correction](inputfilekeywords.md#gwq0correction), [GWnblock](inputfilekeywords.md#gwnblock)

[**Optical Calculations Keywords:**](inputfilekeywords.md#optical-calculations-keywords)
[num_of_bands](inputfilekeywords.md#num_of_bands), [optFDsmear](inputfilekeywords.md#optfdsmear), [opteta](inputfilekeywords.md#opteta), [optdomega0](inputfilekeywords.md#optdomega0), [optomega2](inputfilekeywords.md#optomega2), [optecut](inputfilekeywords.md#optecut), [optnblocks](inputfilekeywords.md#optnblocks)

# All Keywords

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
* FD.

#### Default
PW

#### Example
Mode = 'PW'

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

### whichstrain
#### Keyword type
Python List of Logical values

#### Description
This keyword controls the hich components of strain will be relaxed. There are six independent components indicating the strain are relaxed or not. Here:

* True = relax to zero
* False = fixed

And these six independent components are in order:

* EpsilonX
* EpsilonY, 
* EpsilonZ
* ShearYZ
* ShearXZ
* ShearXY


#### Default
[False, False, False, False, False, False]

#### Example
whichstrain=[True, True, False, False, False, False] #For a x-y 2D nanosheet only first 2 component will be true

---

### MPIcores
#### Keyword type
Integer

#### Description
This keyword controls the number of cores used in calculation. This parameter is not used with `gpawsolve.py`. It is only needed for the `gg.py`.
NOTE: `gg.py` can run `gpawsolve.py` with only `mpirun -np <corenumber>` command. Therefore, for CPUs with hyperthreading support, you can run only the half of your thread number. In the future, `gg.py` will have an option for threads. Please control this variable in future.

#### Default
4

#### Example
MPIcores = 4

## Electronic Calculations Keywords
### fmaxval
#### Keyword type
Float

#### Description
This keyword controls the maximum force tolerance in BFGS type geometry optimization. Unit is eV/Ang.

#### Default
0.05

#### Example
fmaxval = 0.05 # eV/Ang

---

### cut_off_energy
#### Keyword type
Integer

#### Description
This keyword controls the plane wave cut off energy value. Unit is eV. Can be used in PW mode.

#### Default
340 eV

#### Example
cut_off_energy = 500 # eV

---

### kpts_density
#### Keyword type
Float

#### Description
This keyword controls kpoint density. It is deactivated normally. Monkhorst-Pack mesh is used with `kpts_x`, `kpts_y` and `kpts_z` variables. If `kpts_density` is included in an input file, the `kpts_x`, `kpts_y` and `kpts_z` variables will be ignored automatically. Unit is pts per Å^-1.

#### Default
Not used in default.

#### Example
kpts_density = 2.5     # pts per Å^-1

---

### kpts_x
#### Keyword type
Integer

#### Description
This keyword controls the number of kpoints in x direction. If `kpts_density` is included in an input file, the `kpts_x` variable will be ignored automatically. Unit is number of points.

#### Default
5

#### Example
kpts_x = 5

---

### kpts_y
#### Keyword type
Integer

#### Description
This keyword controls the number of kpoints in y direction. If `kpts_density` is included in an input file, the `kpts_y` variable will be ignored automatically. Unit is number of points.

#### Default
5

#### Example
kpts_y = 5

---

### kpts_z
#### Keyword type
Integer

#### Description
This keyword controls the number of kpoints in z direction. If `kpts_density` is included in an input file, the `kpts_z` variable will be ignored automatically. Unit is number of points.

#### Default
5

#### Example
kpts_z = 5

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

### band_path
#### Keyword type
String

#### Description
This keyword controls the path of high-symmetry points in band structure diagram. Use 'G' for Gamma point. 

#### Default
'LGL'

#### Example
band_path = 'GMKG'

---

### band_npoints
#### Keyword type
Integer

#### Description
This keyword controls the number of points between the first and the last high symmetry points. 

#### Default
60

#### Example
band_npoints = 50

---

### energy_max
#### Keyword type
Integer

#### Description
This keyword controls the maximum energy value when the software is used with -d (draw) argument. number of points between the first and the last high symmetry points. Unit is eV.

#### Default
15

#### Example
energy_max = 10 # eV

---

### energy_max
#### Keyword type
Integer

#### Description
This keyword controls the maximum energy value when the software is used with -d (draw) argument. number of points between the first and the last high symmetry points. Unit is eV.

#### Default
15

#### Example
energy_max = 10 # eV

---

### Hubbard
#### Keyword type
Python dictionary

#### Description
This keyword controls the implementation of Hubbard parameter on the related orbitals of related elements. For none use {}. Unit is eV.

#### Default
{}

#### Example
energy_max = {'N': ':p,6.0'} # eV

---

### XC_calc
#### Keyword type
String

#### Description
This keyword controls the which exchange-correlation functional is used in the calculation.Available options are:

* LDA
* PBE
* GLLBSC
* revPBE
* RPBE
* B3LYP
* PBE0  (can be used only with PW-EXX)
* HSE06 (can be used only with PW-EXX)

Because GPAW is using libxc, there are many exchange-correlation functionals available to use. However, the above functionals are used and tested successfully with gpaw-tools. Please try other possible functionals, make us know, send us input files.

#### Default
LDA

#### Example
XC_calc = 'PBE'

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
### GWtype
#### Keyword type
String

#### Description
This keyword controls the type GW calculation. Available options are:

* GW0
* G0W0

#### Default
GW0

#### Example
GWtype = 'GW0'

---

### GWkpoints
#### Keyword type
NumPy Array

#### Description
This keyword represents the kpoint coordinates for the GW calculation.

#### Default
np.array([[0.0, 0.0, 0.0], [1 / 3, 1 / 3, 0], [0.0, 0.0, 0.0]])

#### Example
GWkpoints = np.array([[0.0, 0.0, 0.0], [1 / 3, 1 / 3, 0], [0.0, 0.0, 0.0]])

---

### GWtruncation 
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
GWtruncation = '2D'

---

### GWcut_off_energy
#### Keyword type
Integer

#### Description
This keyword controls the cut off energy value for the GW calculations. Unit is eV.

#### Default
50 eV

#### Example
GWcut_off_energy = 50

---

### GWbandVB
#### Keyword type
Integer

#### Description
This keyword controls the number of the band for the valence band for GW calculations.

#### Default
8 (Default value is not a general value. Please find correct band for your calculation.)

#### Example
GWbandVB = 8

---

### GWbandCB
#### Keyword type
Integer

#### Description
This keyword controls the number of the band for the conduction band for GW calculations.

#### Default
18 (Default value is not a general value. Please find correct band for your calculation.)

#### Example
GWbandCB = 18

---

### GWppa
#### Keyword type
Logical

#### Description
This keyword controls the usage of Plasmon Pole Approximation (PPA) for GW calculations.

#### Default
True

#### Example
GWppa = True

---

### GWq0correction
#### Keyword type
Logical

#### Description
This keyword controls the usage of analytic correction to the q=0 contribution applicable to 2D systems.

#### Default
True

#### Example
GWq0correction = True

---

### GWnblock
#### Keyword type
Logical

#### Description
This keyword controls the behaviour of cutting chi0 into as many blocks to reduce memory requirement as much as possible.

#### Default
True

#### Example
GWnblock = True

## Optical Calculations Keywords
### num_of_bands
#### Keyword type
Integer

#### Description
This keyword controls the number of bands used in optical calculations.

#### Default
16

#### Example
num_of_bands = 8

---

### optFDsmear
#### Keyword type
Float

#### Description
This keyword controls the Fermi Dirac smearing for optical calculations.

#### Default
0.05

#### Example
optFDsmear = 0.02

---

### opteta
#### Keyword type
Float

#### Description
This keyword controls the broadening parameter -eta- used in dielectric function calculations.

#### Default
0.2

#### Example
optFDsmear = 0.05

---

### optdomega0
#### Keyword type
Float

#### Description
This keyword controls the  [Δω0 parameter](https://wiki.fysik.dtu.dk/gpaw/tutorialsexercises/opticalresponse/dielectric_response/dielectric_response.html) for non-linear frequency grid used in dielectric function calculations. Unit is eV.

#### Default
0.1 eV

#### Example
optdomega0 = 0.05 # eV

---

### optomega2
#### Keyword type
Float

#### Description
This keyword controls the  [ω2 parameter](https://wiki.fysik.dtu.dk/gpaw/tutorialsexercises/opticalresponse/dielectric_response/dielectric_response.html) for non-linear frequency grid used in dielectric function calculations. Unit is eV.

#### Default
10.0 eV

#### Example
optomega2 = 2.0 # eV

---

### optecut
#### Keyword type
Float

#### Description
This keyword controls the planewave energy cutoff in dielectric function calculations. Determines the size of dielectric matrix. Unit is eV.

#### Default
10.0 eV

#### Example
optecut = 20.0 # eV

---

### optnblocks
#### Keyword type
Integer

#### Description
This keyword controls the Split matrices in nblocks blocks and distribute them G-vectors or frequencies over processes.

#### Default
4

#### Example
optnblocks = 4

