# gpaw-tools
## Introduction
gpaw-tools is a bunch for python scripts for easy performing of GPAW calculations:
1. A force-field quick optimization script `quickoptimization.py` for preliminary calculations using ASAP3/OpenKIM potentials. 
2. `ciftoase.py` script for transform CIF files to ASE's own Atoms object.
3. To choose better cut off energy, lattice parameter and k points, there are 3 scripts called `Optimize-CutOff.py`, `Optimize-Lattice.py` and `Optimize-KPoints.py`.
4. And, the main solver script `gpawsolver.py` which can be run in PW or LCAO mode. It can do strain minimization, can use several different XCs, can do spin-polarized calculations, can calculate, draw and save tidily DOS and band structures, can calculate and save all-electron densities and can calculate optical properties in a very simple and organized way.

## Usage
We will add more information soon.

## Release notes
Because this is a bunch of scripts, there will be no versioning, no releases. Please try to use the latest github repo zip.

#### August 2021
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

# July 2021 
* `PW-Optical-SingleCoreOnly.py` script for optical calculations.
* `PW-Electronic-changename.py` script for electronic calculations.

# March 2020 
* First scripts for personal usage.
