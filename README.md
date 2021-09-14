# gpaw-tools
## Introduction
*gpaw-tools* is a bunch for python scripts for easy performing of GPAW calculations. It is mostly written for new DFT users who are running codes in their own PCs or on small group clusters.

`gpaw-tools` have:
1. A force-field quick optimization script `quickoptimization.py` for preliminary calculations using ASAP3/OpenKIM potentials. 
2. `ciftoase.py` script for transform CIF files to ASE's own Atoms object.
3. To choose better cut off energy, lattice parameter and k points, there are 3 scripts called `optimize_cutoff.py`, `optimize_latticeparam.py` and `optimize_kpoints.py`.
4. And, the main solver script `gpawsolver.py` which can be run in PW or LCAO mode. It can do strain minimization, can use several different XCs, can do spin-polarized calculations, can calculate, draw and save tidily DOS and band structures, can calculate and save all-electron densities and can calculate optical properties in a very simple and organized way.

## Usage
When you download `gpaw-tools` from GitHub and extract it to a folder you will have a folder structure as:

```
gpaw-tools/
└── benchmarks/
│   └── simple_benchmark_2021.py
├── optimizations/
│   ├── ciftoase.py
│   ├── optimize_cutoff.py
│   ├── optimize_kpoints.py
│   └── optimize_latticeparam.py
├── quick_optimization/
|   └── quickoptimize.py
├── gui_files/
└── gpawsolve.py
└── gg.py
└── config.py
```

### gpawsolve.py
This is the main script for easy and ordered PW/LCAO Calculations with ASE/GPAW. It can run as a stand-alone script or as a command.

#### As a command:
Command line usage: `gpawsolve.py -ochi <inputfile.cif>`

Argument list:
```
-i, --Input  : Use input CIF file
-c, --Config : Use configuration file in the main directory for parameters (config.py) If you do not
               use this argument, parameters will be taken from the Lines 41-86 of gpawsolve.py
-o, --Outdir : Save everything to a output directory with naming /inputfile. 
               If there is no input file given and Atoms object is used in gpawsolve.py file 
               then the directory name will be /gpawsolve. If you change gpawsolve.py name to 
               anyname.py then the directory name will be /anyname
-h --Help    : Help
 ```
 
 #### As a stand alone script
 * Change the parameters from Line 41 to 86 for each simulation OR change `config.py` once and use `-c` argument.
 * If you want to use CIF files for structure, use `-i` argument like `gpawsolve.py -i structurefile.cif`.
 * If you want to use ASE atoms method for structure, just copy/paste your `Atoms` info into the part mentioned with "Bulk Structure".
 * If you have CIF file but want to use Atoms method you can use `CIF-to-ASE/ciftoase.py` to convert your CIF files to ASE Atoms.
 * If you use Atoms method, change the name of `gpawsolve.py` to your simulation name like `graphene7x7-Fe-onsite32.py`. The naming will be used for naming of all output/result files.
 * If you use CIF file as an input, name of the input file will be used for naming of all output/result files.
 
 #### How to run?
 Change `<core_number>` with core numbers/threads to use. For getting a maximum performance from your PC you can use `total number of cores(or threads) - 1`. or `total RAM/2Gb` as a `<core_number>`

Usage:
`$ gpaw -P<core_number> python gpawsolve.py`

For AMD CPUs or using Intel CPUs without hyperthreading:
`$ mpirun -n <core_number> gpaw python gpawsolve.py`

For using all threads provided by Intel Hyperthreading technology
`$ mpirun --use-hwthread-cpus -n <core_number> gpaw python gpawsolve.py`

#### Calculation selector

| Method | Strain_minimization | Several XCs | Spin polarized | DOS | Band | Electron Density | Optical |
| ------ | ------------------- | ----------- | -------------- | --- | ---- | ---------------- | ------- |
|   PW   | Yes                 | Yes         | Yes            | Yes | Yes  | Yes              | Yes     |
|  LCAO  | No                  | No          | No             | Yes | Yes  | Yes              | No      |

### gg.py
Basic DFT calculations can be done graphically with the script `gg.py`. This script is behaving as a GUI to run `gpawsolve.py` script. To execute the GUI, type simply:
  python gg.py

### quick_optimize/quickoptimize.py
Inter-atomic potentials are useful tool to perform a quick geometric optimization of the studied system before starting a precise DFT calculation. The `quickoptimize.py` script is written for geometric optimizations with inter-atomic potentials. The bulk configuration of atoms can be provided by the user in the script as an ASE Atoms object or given as an argument for the CIF file. A general potential is given for any calculation. However, user can provide the necessary OpenKIM potentialby changing the related line in the script.

Mainly, quickoptimize.py is not related to GPAW. However it is dependent to ASAP3/OpenKIM and Kimpy. Therefore, the user must install necessary libraries before using the script:

    pip install --upgrade --user ase asap3
    sudo add-apt-repository ppa:openkim/latest
    sudo apt-get update
    sudo apt-get install libkim-api-dev openkim-models libkim-api2 pkg-config
    pip3 install kimpy

The script can be called as: from the command line  in the script itself:

    python quickoptimize.py                   (if the user wants to provide structure as ASE Atoms object)
    python quickoptimize.py <inputfile.cif>   (if the user wants to provide structure as a CIF file


### optimizations/ciftoase.py
For `quickoptimize.py` or other optimization scripts, user may need to give ASE Atoms object instead of using a CIF file. This script changes a CIF file information to ASE Atoms object. Because there is a problem in the read method of ASE.io, sometimes it can give a double number of atoms. If the user lives this kind of problem, there is a setting inside the script. User can run the script like:

    python ciftoase.py <inputfile.cif>

Result will be printed to screen and will be saved as `inputfile.py` in the same folder.

### optimizations/optimize_cutoff (and kpoints)(and latticeparam).py
Users must provide ASE Atoms object and simply insert the object inside these scripts. With the scripts, the user can do convergence tests for cut-off energy, k-points and can calculate the energy dependent lattice parameter values. These codes are mainly based on Prof. J. Kortus, R. Wirnata's Electr. Structure & Properties of Solids course notes and GPAW's tutorials. Scripts can easily called with MPI as:

    gpaw -P <core_number> python optimize_cutoff.py
    gpaw -P <core_number> python optimize_kpoints.py
    gpaw -P <core_number> python optimize_latticeparam.py

### benchmarks/
GPAW has many test scripts for many cases. However, new users may need something easy to run and compare. Some very easy single file test scripts will be listed [here](https://github.com/lrgresearch/gpaw-tools/tree/main/benchmarks) with some hardware benchmark information. Your timings are always welcomed.

## Release notes

### Version 21.9.0
#### September 2021
* Corrected `quickoptimize.py` behaviour
* Many code quality and folder structure improvements.
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
