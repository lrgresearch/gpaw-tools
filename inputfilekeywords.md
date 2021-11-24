---
layout: default
navigation_weight: 6
title: Input File Keywords
---

# Input File Keywords

```
NOTE: This page is under construction.
```
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
## Electronic Properties Calculation Variables
### fmaxval
#### Keyword type
Floating point

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
Floating point

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
