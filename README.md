# gpaw-tools
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Issues:](https://img.shields.io/github/issues/lrgresearch/gpaw-tools)](https://github.com/lrgresearch/gpaw-tools/issues)
[![Pull requests:](https://img.shields.io/github/issues-pr/lrgresearch/gpaw-tools)](https://github.com/lrgresearch/gpaw-tools/pulls)
[![Latest version:](https://img.shields.io/github/v/release/lrgresearch/gpaw-tools)](https://github.com/lrgresearch/gpaw-tools/releases/)
![Release date:](https://img.shields.io/github/release-date/lrgresearch/gpaw-tools)
[![Commits:](https://img.shields.io/github/commit-activity/m/lrgresearch/gpaw-tools)](https://github.com/lrgresearch/gpaw-tools/commits/main)
[![Last Commit:](https://img.shields.io/github/last-commit/lrgresearch/gpaw-tools)](https://github.com/lrgresearch/gpaw-tools/commits/main)
[![Total alerts](https://img.shields.io/lgtm/alerts/g/lrgresearch/gpaw-tools.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/lrgresearch/gpaw-tools/alerts/)
[![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/lrgresearch/gpaw-tools.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/lrgresearch/gpaw-tools/context:python)
## Introduction
*gpaw-tools* is a collection of python scripts that use ASE and GPAW for performing Density Functional Theory (DFT) calculations. Its aim is to lower the entry barrier and to provide an easy-to-use command line and graphical user interfaces for GPAW. It is mostly written for new DFT users who are running codes on their own PCs or on small group clusters.

`gpaw-tools` have:
1. A force-field quick optimization script `quickoptimization.py` for preliminary calculations using ASAP3/OpenKIM potentials. 
2. `ciftoase.py` script for transform CIF files to ASE's own Atoms object.
3. To choose better cut off energy, lattice parameter and k points, there are 3 scripts called `optimize_cutoff.py`, `optimize_latticeparam.py` and `optimize_kpoints.py`.
4. The main solver script `gpawsolver.py` which can be run in PW (also with GW and EXX) or LCAO mode. It can do structure optimization, Equation of State and elastic tensor calculations, can use several different XCs, can do spin-polarized calculations, can calculate, draw and save tidily DOS and band structures, can calculate and save all-electron densities and can calculate optical properties in a very simple and organized way.
5. A simple Graphical User Interface (GUI) for `gpawsolve.py` (and also you may say that GUI for GPAW) which is called `gg.py`. 

## Usage
### Installation
When you download `gpaw-tools-main` from GitHub and extract it to a folder you will have a folder structure as:

```
gpaw-tools-main/
└── benchmarks/
│   └── simple_benchmark_2021.py
├── examples/
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
```
To make the `gpawsolve.py` and `gg.py` as system-wide commands, user must include the `gpaw-tools-main` folder to the $PATH variable in the `.bashrc` file. In case of user  downloaded and extracted the `gpaw-tools-main` file to user's home directory, and to make the change permanent, user must need to define the $PATH variable in the shell configuration file `.bashrc` as

    export PATH="/home/username/gpaw-tools-main:$PATH"
    
also you may need to give execute rights to `gpawsolve.py` and `gg.py` to execute these scripts as a command

    cd /home/username/gpaw-tools-main
    chmod +x gpawsolve.py gg.py

### gpawsolve.py
This is the main script for easy and ordered PW/LCAO Calculations with ASE/GPAW. It can run as a stand-alone script or as a command.

#### As a command:
Command line usage: `gpawsolve.py -v -o -r -d -i <inputfile.py> -h -g <geometryfile.cif>`

Argument list:
```
-g, --geometry   : Use a CIF file for geometry
-i, --input      : Use an input file for variables (input.py) If you do not use this argument, parameters 
                   will be taken from the related lines of gpawsolve.py. Visit "Input File Keywords" webpage for more. 
-o, --outdir     : Save everything to a output directory with naming /inputfile. 
                   If there is no input file given and Atoms object is used in gpawsolve.py file 
                   then the directory name will be /gpawsolve. If you change gpawsolve.py name to 
                   anyname.py then the directory name will be /anyname
-h --help        : Help
-d --drawfigures : Draws DOS and band structure figures at the end of calculation.
-r --restart     : Passing ground calculations and continue with the next required calculation.
-v --version     : Version information of running code and the latest stable code. Also gives download link.
 ```
  You can put ASE Atoms object in to your config file and therefore can use it like an input file. As an example please note the example at: `examples\Bulk-aluminum` folder.
  
 #### As a stand alone script
 * Change the parameters from related lines for each simulation OR create an input file (as you can see in examples) once and use `-i` argument.
 * If you want to use CIF files for structure, use `-g` argument like `gpawsolve.py -g geometryfile.cif`.
 * If you want to use ASE atoms method for structure, just copy/paste your `Atoms` info into the part mentioned with "Bulk Structure".
 * If you have CIF file but want to use Atoms method you can use `CIF-to-ASE/ciftoase.py` to convert your CIF files to ASE Atoms.
 * If you use Atoms method, change the name of `gpawsolve.py` to your simulation name like `graphene7x7-Fe-onsite32.py`. The naming will be used for naming of all output/result files.
 * If you use CIF file as an geometry, name of the input file will be used for naming of all output/result files.
 * **Performance note:** When you want to use `gpawsolve.py` as a script, you can copy `gpawsolve.py` to your working folder where your config file and input file are ready. You must rename `gpawsolve.py` to something else like `gpawsolve1.py` or `gs-graphene.py`, something you like and then you can now run `gpaw -P<core> python gpawsolve1.py <args>` type command. Initializing with gpaw command in your system will give you better parallel computing, therefore shorter computation times. Initialization with gpaw can not be done when `gpawsolve.py` is used as command, because of the structure of initialization of Gpaw, as we know. If you know a solution from the point of view of gpaw-tools, please use issues to discuss or pull request for a solution.
 
 #### How to run?
 Change `<core_number>` with core numbers to use. For getting a maximum performance from your PC you can use `total number of cores - 1` or `total RAM/2Gb` as a `<core_number>`. For CPUs supporting hyperthreading, users can use more than one instance of `gpawsolve.py` to achive maximum efficiency. 

Usage:
`$ mpirun -np <core_number> gpawsolve.py <args>`

#### Calculation selector

| Method | Structure optim.    | Different XCs | Spin polarized | Elastic | DOS | DFT+U | Band | Electron Density | Optical |
| ------ | ------------------- | ------------- | -------------- | ------- | --- | ----- | ---- | ---------------- | ------- |
|   PW   | Yes (No for GLLBSC) | Yes           | Yes            | Yes     | Yes | Yes   | Yes  | Yes              | Yes     |
| PW-G0W0| Yes                 | Yes           | No             | No      | No  | No    | Yes  | No               | No      |
| PW-EXX*| Yes (with PBE)      | No            | No             | No      | No  | No    | No   | No               | No      |
|  LCAO  | No                  | No            | No             | No      | Yes | Yes   | Yes  | Yes              | No      |

*: Just some ground state energy calculations for PBE0 and HSE06.

### gg.py
Basic DFT calculations can be done graphically with the script `gg.py`. This script is behaving as a GUI to run `gpawsolve.py` script. To execute the GUI, type simply:
  gg.py

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
    python quickoptimize.py <geometryfile.cif>   (if the user wants to provide structure as a CIF file


### optimizations/ciftoase.py
For `quickoptimize.py` or other optimization scripts, user may need to give ASE Atoms object instead of using a CIF file. This script changes a CIF file information to ASE Atoms object. Because there is a problem in the read method of ASE.io, sometimes it can give a double number of atoms. If the user lives this kind of problem, there is a setting inside the script. User can run the script like:

    python ciftoase.py <geometryfile.cif>

Result will be printed to screen and will be saved as `geometryfile.py` in the same folder.

### optimizations/optimize_cutoff (and kpoints)(and latticeparam).py
Users must provide ASE Atoms object and simply insert the object inside these scripts. With the scripts, the user can do convergence tests for cut-off energy, k-points and can calculate the energy dependent lattice parameter values. These codes are mainly based on Prof. J. Kortus, R. Wirnata's Electr. Structure & Properties of Solids course notes and GPAW's tutorials. Scripts can easily called with MPI as:

    gpaw -P <core_number> python optimize_cutoff.py
    gpaw -P <core_number> python optimize_kpoints.py
    gpaw -P <core_number> python optimize_latticeparam.py

### benchmarks/
GPAW has many test scripts for many cases. However, new users may need something easy to run and compare. Some very easy single file test scripts will be listed [here](https://github.com/lrgresearch/gpaw-tools/tree/main/benchmarks) with some hardware benchmark information. Your timings are always welcomed.

## examples/
There are [some example calculations](https://github.com/lrgresearch/gpaw-tools/tree/main/examples) given with different usage scenarios. Please send us more calculations to include in this folder.

## Input File Keywords
There are many keywords can be used in input files. You can find more at [here](https://www.lrgresearch.org/gpaw-tools/inputfilekeywords/)

## Release notes
Release notes are listed at [here](https://www.lrgresearch.org/gpaw-tools/releasenotes/).

## Citing
Please do not forget that, gpaw-tools is a UI/GUI software. For the main DFT calculations, it uses ASE and GPAW. It also uses Elastic python package for elastic tensor solutions and ASAP with KIM database for interatomic interaction calculations. Therefore, you must know what you use, and cite them properly. Here, the basic citation information of each packages are given. There are many other packages needed to be cited. With GPAW, you may needed to cite LibXC, LCAO  Please visit their pages for many other citation possibilities. 

* **ASE** : Ask Hjorth Larsen et al. "[The Atomic Simulation Environment—A Python library for working with atoms](https://doi.org/10.1088/1361-648X/aa680e)" J. Phys.: Condens. Matter Vol. 29 273002, 2017.
* **GPAW**: J. J. Mortensen, L. B. Hansen, and K. W. Jacobsen "[Real-space grid implementation of the projector augmented wave method](https://doi.org/10.1103/PhysRevB.71.035109)" Phys. Rev. B 71, 035109 (2005) and J. Enkovaara, C. Rostgaard, J. J. Mortensen et al. "[Electronic structure calculations with GPAW: a real-space implementation of the projector augmented-wave method](https://doi.org/10.1088/0953-8984/22/25/253202)" J. Phys.: Condens. Matter 22, 253202 (2010) [OTHER POSSIBLE CITATION](https://wiki.fysik.dtu.dk/gpaw/faq.html#citation-how-should-i-cite-gpaw)
* **KIM** : E. B. Tadmor, R. S. Elliott, J. P. Sethna, R. E. Miller and C. A. Becker, "The Potential of Atomistic Simulations and the Knowledgebase of Interatomic Models" JOM, 63, 17 (2011). doi:10.1007/s11837-011-0102-6. [OTHER POSSIBLE CITATION](https://openkim.org/how-to-cite/)

For `gpaw-tools` usage you can use the following citation:

* S.B. Lisesivdin, B. Sarikavak-Lisesivdin "[gpaw-tools – higher-level user interaction scripts for GPAW calculations and interatomic potential based structure optimization](https://doi.org/10.1016/j.commatsci.2022.111201)" Comput. Mater. Sci. 204, 111201 (2022).

## Licensing
This project is licensed under the terms of the [MIT license](https://opensource.org/licenses/MIT).
