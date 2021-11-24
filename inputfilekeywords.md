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
False

#### Example
fmaxval = 0.05 # eV/Ang
