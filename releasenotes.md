---
layout: default
navigation_weight: 5
title: Release Notes
---

## Release notes

### Development version

* In GW calculations, calculation could not be done because of interpolation in drawing the figure when the data did not have a minimum of 3 points. Now there is a variable to use interpolation or not.
* In GW calculations, `gpawsolve.py` can write quasiparticle energies to a file separately.
* DFT+U calculation ability is added for PW and LCAO modes.
* `gg.py` is better now. It is compatible with new arguments and removed variables. It can run in any directory and handling of `GWkpoints` and `GWtruncation` variables are correct now.
* A new argument '-d' is added. This argument makes script draw the DOS and band calculation results. In the past, it was a varible in the config file.
* `WantedCIFexport` variable is removed. 
* A new argument '-r' is added. This argument makes script pass the ground state calculations and continue with the next calculation.

### Version 21.10.1

* An important bug made it impossible to work with existing examples with `gg.py`. It is now resolved.

### Version 21.10.0

* 5 different examples are added to show simply different usage cases.
* Initializing magnetic moment problem is solved.
* Version argument is added.
* GW parameters are also added to `gg.py`
* Add some optical parameters to config files and `gg.py`.
* Major change: `gpawsolve.py` and `gg.py` are now working as commands. Config files can be used as general input files, you can put ASE Atom object plus every parameter that `gpawsolve.py` accepts. If you want, you can provide a CIF file instead of using an atom object. You can run `gpawsolve.py` and `gg.py` from any folder.
* Optical: Refractive index, extinction index, absorption, and reflectivity calculations.
* `gg.py` is now opening ase gui when the user clicks the structure image.
* New argument parsing scheme for better future usages.
* Very basic PW-EXX mode with HSE06 and PBE06. (Only some ground-state calculations.)
* Adding GW0 and G0W0-GW0 selector.
* Adding GW approximation to `gpawsolve.py` (only bands).
* Many other small corrections.

### Version 21.9.0

* Corrected `quickoptimize.py` behaviour.
* Many code quality and folder structure improvements.
* Comment additions to code.
* Better README.md.
* `gg.py` which is a GUI for gpaw-tools is added to project. It can do all `gpawsolve.py`'s features in a graphical way!
* `gpawsolve.py` can be run solely as a command now (This is needed for a GUI project).
* All three scripts`PW-Electronic.py`, `LCAO-Electronic.py` and `PW-Optical-SingleCoreOnly.py` scripts becomes a single for-all-case script: `gpawsolve.py`.
* `PW-Electronic-changename.py` script becomes `PW-Electronic.py`.
* Spin-polarized results in `PW-Electronic-changename.py` script.
* All-electron density calculations in `PW-Electronic-changename.py`.
* CIF Export in `PW-Electronic-changename.py` script.
* Better parallel computation.
* Several XCs available for PW.
* `LCAO-Electronic.py` script.
* Strain minimization in PW only. 
* BFGS to LBFGS, Small many changes have been done.

### Preversion
* `PW-Optical-SingleCoreOnly.py` script for optical calculations.
* `PW-Electronic-changename.py` script for electronic calculations. 
* First scripts for personal usage.
