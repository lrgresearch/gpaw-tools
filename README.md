# gpaw-tools
## Introduction
gpaw-tools is a bunch for python scripts for easy performing of GPAW calculations:
1. A force-field quick optimization script `quickoptimization.py` for preliminary calculations using ASAP3/OpenKIM potentials. 
2. `ciftoase.py` script for transform CIF files to ASE's own Atoms object.
3. To choose better cut off energy, lattice parameter and k points, there are 3 scripts called `Optimize-CutOff.py`, `Optimize-Lattice.py` and `Optimize-KPoints.py`.
4. And, the main solver script `gpawsolver.py` which can be run in PW or LCAO mode. It can do strain minimization, can use several different XCs, can do spin-polarized calculations, can calculate, draw and save tidily DOS and band structures, can calculate and save all-electron densities and can calculate optical properties in a very simple and organized way.

## Usage
When you download `gpaw-tools` from GitHub and extract it to a folder you will have a folder structure as:

```
gpaw-tools/
└── Benchmarks/
│   └── GPAWSimpleBenchmark2021.py
├── CIF-to-ASE/
│   ├── ciftoase.py
│   └── visualize.py
├── Cutoff-Lattice-Kpoint-optimizations/
│   ├── Optimize-CutOff.py
│   ├── Optimize-KPoints.py
│   └── Optimize-Lattice.py
├── QuickOptimize/
|   └── quickoptimize.py
├── gui_files/
└── gpawsolve.py
└── gg.py
└── config.py
```

### gpawsolve.py
This is the main script for easy and ordered PW/LCAO Calculations with ASE/GPAW

Command line usage: `gpawsolve.py -ochi <inputfile.cif>`

Argument list:
```
'-i, --Input'  : Use input CIF file
-c, --Config : Use configuration file in the main directory for parameters (config.py)
-o, --Outdir : Save everything to a output directory with naming /inputfile. 
               If there is no input file given and Atoms object is used in gpawsolve.py file 
               then the directory name will be /gpawsolve. If you change gpawsolve.py name to 
               anyname.py then the directory name will be /anyname
-h --Help    : Help
 ```
 
 #### General command structure:
 Change <core_number> with core numbers/threads to use. For getting a maximum performance from your PC you can use `total number of cores(or threads) - 1`. or `total RAM/2Gb` as a <core_number>

Usage:
`$ gpaw -P8 python gpawsolve.py`

For AMD CPUs or using Intel CPUs without hyperthreading: (Example CPU is intel here, 4 cores or 8 threads)
`$ mpirun -n 4 gpaw python gpawsolve.py`

For using all threads provided by Intel Hyperthreading technology
`$ mpirun --use-hwthread-cpus -n 8 gpaw python gpawsolve.py`

#### Calculation selector

| Method | Strain_minimization | Several XCs | Spin polarized | DOS | Band | Electron Density | Optical |
| ------ | ------------------- | ----------- | -------------- | --- | ---- | ---------------- | ------- |
|   PW   | Yes                 | Yes         | Yes            | Yes | Yes  | Yes              | Yes     |
|  LCAO  | No                  | No          | No             | Yes | Yes  | Yes              | No      |

## Release notes
Because this is a bunch of scripts, there will be no strict versioning, rolling releases. Please try to use the latest github repo zip.

#### September 2021
* Comment additions to code.
* Better README.md

#### August 2021
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

#### July 2021 
* `PW-Optical-SingleCoreOnly.py` script for optical calculations.
* `PW-Electronic-changename.py` script for electronic calculations.

#### March 2020 
* First scripts for personal usage.

## Licensing
This project is licensed under the terms of the MIT license.
