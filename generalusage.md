---
layout: default
nav_order: 5
title: General Usage
---

# General Usage
When you download `gpaw-tools` from GitHub and extract it to a folder you will have a folder structure as:

```
gpaw-tools/
└── benchmarks/
│   └── simple_benchmark_2021.py
├── examples/
├── optimizations/
│   ├── ciftoase.py
│   ├── optimize_cutoff.py
│   ├── optimize_kpoints.py
│   ├── optimize_kptsdensity.py
│   └── optimize_latticeparam.py
├── gui_files/
└── asapsolve.py
└── gpawsolve.py
└── gg.py
└── shrinkgpw.py
```

## gpawsolve.py
This is the main script for easy and ordered PW/LCAO Calculations with ASE/GPAW. It can run as a stand-alone script or as a command.

### As a command:
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
 
### As a stand alone script
* Change the parameters from related lines for each simulation OR create an input file (as you can see in examples) once and use `-i` argument.
 * If you want to use CIF files for structure, use `-g` argument like `gpawsolve.py -g geometryfile.cif`.
 * If you want to use ASE atoms method for structure, just copy/paste your `Atoms` info into the part mentioned with "Bulk Structure".
 * If you have CIF file but want to use Atoms method you can use `CIF-to-ASE/ciftoase.py` to convert your CIF files to ASE Atoms.
 * If you use Atoms method, change the name of `gpawsolve.py` to your simulation name like `graphene7x7-Fe-onsite32.py`. The naming will be used for naming of all output/result files.
 * If you use CIF file as an geometry, name of the input file will be used for naming of all output/result files.
 * **Performance note:** When you want to use `gpawsolve.py` as a script, you can copy `gpawsolve.py` to your working folder where your config file and input file are ready. You must rename `gpawsolve.py` to something else like `gpawsolve1.py` or `gs-graphene.py`, something you like and then you can now run `gpaw -P<core> python gpawsolve1.py <args>` type command. Initializing with gpaw command in your system will give you better parallel computing, therefore shorter computation times. Initialization with gpaw can not be done when `gpawsolve.py` is used as command, because of the structure of initialization of Gpaw, as we know. If you know a solution from the point of view of gpaw-tools, please use issues to discuss or pull request for a solution.
 
### How to run?
Change `<core_number>` with core numbers to use. For getting a maximum performance from your PC you can use `total number of cores - 1` or `total RAM/2Gb` as a `<core_number>`. For CPUs supporting hyperthreading, users can use more than one instance of `gpawsolve.py` to achive maximum efficiency. 

Usage:
`$ mpirun -np <core_number> gpawsolve.py <args>`

### Calculation selector (Not complete, not up-to-date information)

 | Method | XCs                 | Structure optim. | Spin polarized | Ground | Elastic | DOS | DFT+U | Band | Electron Density | Optical |
 | ------ | ------------------- | ---------------- | -------------- | ------ | ------- | --- | ----- | ---- | ---------------- | ------- |
 |   PW   | Local and LibXC     | Yes              | Yes            | Yes    | Yes     | Yes | Yes   | Yes  | Yes              | Yes     |
 |   PW   | GLLBSC / M          | No               | Yes            | Yes    | Yes     | Yes | No    | Yes  | Yes              | Yes     |
 |   PW   | HSE03, HSE06        | No               | Yes            | Yes    | n/a     | Yes | No    | No   | No               | No      |
 | PW-G0W0| Local and LibXC     | No               | No             | Yes    | No      | No  | No    | Some | No               | No      |
 | PW-EXX*| B3LYP, PBE0         | Yes (with PBE)   | No             | Yes    | No      | No  | No    | No   | No               | No      |
 |  LCAO  | Local and LibXC     | Yes              | Yes            | Yes    | Yes     | Yes | Yes   | Yes  | Yes              | No      |

*: Just some ground state energy calculations for PBE0 and HSE06.

## gg.py
Basic DFT calculations can be done graphically with the script `gg.py`. This script is behaving as a GUI to run `gpawsolve.py` script. To execute the GUI, type simply:
  gg.py

## asapsolve.py
Inter-atomic potentials are useful tool to perform a quick geometric optimization of the studied system before starting a precise DFT calculation. The `asapsolve.py` script is written for geometric optimizations with inter-atomic potentials. The bulk configuration of atoms can be provided by the user given as CIF file. A general potential is given for any calculation. However, user can provide the necessary OpenKIM potential by changing the related line in the input file.

Mainly, `asapsolve.py` is not related to GPAW. However it is dependent to ASAP3/OpenKIM and Kimpy. Therefore, the user must install necessary libraries before using the script:

    pip install --upgrade --user ase asap3
    sudo add-apt-repository ppa:openkim/latest
    sudo apt-get update
    sudo apt-get install libkim-api-dev openkim-models libkim-api2 pkg-config
    pip3 install kimpy

Main usage is:

`$ asapsolve.py <args>`

### Arguments

`asapsolve.py -v -i <inputfile.py> -h -g <geometryfile.cif>`

Argument list:
```
-g, --geometry   : Use a CIF file for geometry
-i, --input      : Use an input file for variables (input.py) 

-h --help        : Help
-v --version     : Version information of running code and the latest stable code. Also gives download link.
 ```
 
## optimizations/ciftoase.py
For `quickoptimize.py` or other optimization scripts, user may need to give ASE Atoms object instead of using a CIF file. This script changes a CIF file information to ASE Atoms object. Because there is a problem in the read method of ASE.io, sometimes it can give a double number of atoms. If the user lives this kind of problem, there is a setting inside the script. User can run the script like:

    python ciftoase.py <geometryfile.cif>

Result will be printed to screen and will be saved as `geometryfile.py` in the same folder.

## optimizations/optimize_cutoff (and kpoints)(and latticeparam).py
Users must provide ASE Atoms object and simply insert the object inside these scripts. With the scripts, the user can do convergence tests for cut-off energy, k-points and can calculate the energy dependent lattice parameter values. These codes are mainly based on Prof. J. Kortus, R. Wirnata's Electr. Structure & Properties of Solids course notes and GPAW's tutorials. Scripts can easily called with MPI as:

    gpaw -P <core_number> python optimize_cutoff.py
    gpaw -P <core_number> python optimize_kpoints.py
    gpaw -P <core_number> python optimize_latticeparam.py

`optimize_latticeparam.py` can perform simultaneous calculation for lattice parameters a and c. And can also draw 3D contour graph for Energy versus lattice parameters (a and c).

## benchmarks/
GPAW has many test scripts for many cases. However, new users may need something easy to run and compare. Some very easy single file test scripts will be listed [here](https://github.com/lrgresearch/gpaw-tools/tree/main/benchmarks) with some hardware benchmark information. Your timings are always welcomed.

## examples/
There are some example calculations given with different usage scenarios in the code. Please send us more calculations to include.

### GPAW Example List

| Name              | Notes  | 
| ----------------- | ------ |
| Bulk-Al-noCIF     | Ground, DOS and Band calculations of Bulk Aluminum with PW. Positions are given with Atom object.          |
| Cr2O-spin         |Spin-dependent electronic properties of CrO2 |
| Graphene-LCAO     | Pristine graphene and graphene with defect with LCAO. Uses single config for two calculations. |
| MoS2-GW           | GW Approximation calculation for MoS2 |
| Si-2atoms-optical | Three step calculation. First step ground, DOS and Band calculations. Second and third steps for RPA and BSE optical calculations, respectively. Structure is given with CIF file. |
| ZnO with DFT+U    | Wurtzite ZnO calculation with DFT+U. Positions are given with Bulk object. Hubbard params are: O-p: 7eV, Zn-d: 10eV|
| TiC-elastic-electronic | Elastic (EoS and Elastic Tensor) and Electronic Properties of Rocksalt TiC |
| Si-with-HSE | Ground state, DOS and band structure of Si with HSE06 Hybrid XC | 

### ASAP3 Example List

| Name              | Notes  | 
| ----------------- | ------ |
| ASAP3-Example     | Germanene nanosheet example with a general potential.          |
