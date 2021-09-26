---
layout: default
navigation_weight: 4
title: Usage
---

# Usage
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

## gpawsolve.py
This is the main script for easy and ordered PW/LCAO Calculations with ASE/GPAW. It can run as a stand-alone script or as a command.

### As a command:
Command line usage: `gpawsolve.py -o -c <configfile.py> -h -i <inputfile.cif>`

Argument list:
```
-i, --Input  : Use input CIF file
-c, --Config : Use configuration file in the main directory for parameters (config.py) If you do not
               use this argument, parameters will be taken from the related lines of gpawsolve.py
-o, --Outdir : Save everything to a output directory with naming /inputfile. 
               If there is no input file given and Atoms object is used in gpawsolve.py file 
               then the directory name will be /gpawsolve. If you change gpawsolve.py name to 
               anyname.py then the directory name will be /anyname
-h --Help    : Help
 ```
 
### As a stand alone script
 * Change the parameters for each simulation OR change `config.py` once and use `-c` argument.
 * If you want to use CIF files for structure, use `-i` argument like `gpawsolve.py -i structurefile.cif`.
 * If you want to use ASE atoms method for structure, just copy/paste your `Atoms` info into the part mentioned with "Bulk Structure".
 * If you have CIF file but want to use Atoms method you can use `CIF-to-ASE/ciftoase.py` to convert your CIF files to ASE Atoms.
 * If you use Atoms method, change the name of `gpawsolve.py` to your simulation name like `graphene7x7-Fe-onsite32.py`. The naming will be used for naming of all output/result files.
 * If you use CIF file as an input, name of the input file will be used for naming of all output/result files.
 
### How to run?
 Change `<core_number>` with core numbers/threads to use. For getting a maximum performance from your PC you can use `total number of cores(or threads) - 1`. or `total RAM/2Gb` as a `<core_number>`

Usage:
For AMD CPUs or using Intel CPUs without hyperthreading:
`$ mpirun -np <core_number> python gpawsolve.py <args>`

For using all threads provided by Intel Hyperthreading technology
`$ mpirun --use-hwthread-cpus -np <core_number> python gpawsolve.py <args>`

### Calculation selector

| Method | Strain_minimization | Several XCs | Spin polarized | DOS | Band | Electron Density | Optical |
| ------ | ------------------- | ----------- | -------------- | --- | ---- | ---------------- | ------- |
|   PW   | Yes                 | Yes         | Yes            | Yes | Yes  | Yes              | Yes     |
| PW-G0W0| Yes                 | Yes         | No             | No  | Yes  | No               | No      |
| PW-EXX*| Yes (with PBE)      | No          | No             | No  | No   | No               | No      |
|  LCAO  | No                  | No          | No             | Yes | Yes  | Yes              | No      |

*: Just some ground state energy calculations for PBE0 and HSE06.

## gg.py
Basic DFT calculations can be done graphically with the script `gg.py`. This script is behaving as a GUI to run `gpawsolve.py` script. To execute the GUI, type simply:
  python gg.py

## quick_optimize/quickoptimize.py
Inter-atomic potentials are useful tool to perform a quick geometric optimization of the studied system before starting a precise DFT calculation. The `quickoptimize.py` script is written for geometric optimizations with inter-atomic potentials. The bulk configuration of atoms can be provided by the user in the script as an ASE Atoms object or given as an argument for the CIF file. A general potential is given for any calculation. However, user can provide the necessary OpenKIM potentialby changing the related line in the script.

Mainly, quickoptimize.py is not related to GPAW. However it is dependent to ASAP3/OpenKIM and Kimpy. Therefore, the user must install necessary libraries before using the script.

The script can be called as: from the command line  in the script itself:

    python quickoptimize.py                   (if the user wants to provide structure as ASE Atoms object)
    python quickoptimize.py <inputfile.cif>   (if the user wants to provide structure as a CIF file


## optimizations/ciftoase.py
For `quickoptimize.py` or other optimization scripts, user may need to give ASE Atoms object instead of using a CIF file. This script changes a CIF file information to ASE Atoms object. Because there is a problem in the read method of ASE.io, sometimes it can give a double number of atoms. If the user lives this kind of problem, there is a setting inside the script. User can run the script like:

    python ciftoase.py <inputfile.cif>

Result will be printed to screen and will be saved as `inputfile.py` in the same folder.

## optimizations/optimize_cutoff (and kpoints)(and latticeparam).py
Users must provide ASE Atoms object and simply insert the object inside these scripts. With the scripts, the user can do convergence tests for cut-off energy, k-points and can calculate the energy dependent lattice parameter values. These codes are mainly based on Prof. J. Kortus, R. Wirnata's Electr. Structure & Properties of Solids course notes and GPAW's tutorials. Scripts can easily called with MPI as:

    gpaw -P <core_number> python optimize_cutoff.py
    gpaw -P <core_number> python optimize_kpoints.py
    gpaw -P <core_number> python optimize_latticeparam.py

## benchmarks/
GPAW has many test scripts for many cases. However, new users may need something easy to run and compare. Some very easy single file test scripts will be listed [here](https://github.com/lrgresearch/gpaw-tools/tree/main/benchmarks) with some hardware benchmark information. Your timings are always welcomed.
