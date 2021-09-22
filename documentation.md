# Documentation
## Installation
 You must have `ase` and `gpaw` codes on your computer. You can find more information about installation of [ASE](https://wiki.fysik.dtu.dk/ase/install.html) and [GPAW](https://wiki.fysik.dtu.dk/gpaw/install.html) from their related sites.
 However, for a simple installation on a Windows 10 computer with WSL, a small installation information is given below:

### Full Installation of WSL, ASE, GPAW and gpaw-tools on a Windows 10 system

#### Windows Subsystem for Linux (WSL)
The Windows Subsystem for Linux allows developers to run a wide range of Linux-based apps and utilities on Windows without the need for a traditional virtual machine or a dualboot setup. More information about installation can be found [here](https://docs.microsoft.com/en-us/windows/wsl/install).

We are suggesting Ubuntu 20.04 LTS version for gpaw-tools studies.

#### WSL Linux Distribution Version Control
After installing WSL, we need to control the WSL version of Linux distribution installed in your system. Open Powershell and type:

    wsl --list --verbose 

Result will be something like:

    NAME            STATE           VERSION
    Ubuntu-20.04    Running         1

Your WSL distribution version must be "1", not "2". If it is not, change the WSL distribution version with command:

    wsl --set-version Ubuntu-20.04 1

#### ASE/GPAW Installation
To use a X server on windows, firstly you need to install few files (not necessary if you do not use matplotlib features.)

    sudo apt install libglu1-mesa-dev freeglut3 freeglut3-dev mesa-common-dev

You also need to install [XMing](https://sourceforge.net/projects/xming/) on Windows. In new WSL installations, the following line must be added to `.bashrc` file

    export LIBGL_ALWAYS_INDIRECT=0
    export DISPLAY=localhost:0.0
    
also you need Tk library

    sudo apt install python3-tk

To further package installations, we need PIP installer

    sudo apt install python3-pip

Install ASE and other math, parallel, dev libraries

    pip3 install --upgrade --user ase
    
At this point, PIP can give some warnings as:

    WARNING: The scripts f2py, f2py3 and f2py3.8 are installed in '/home/YOURUSERNAME/.local/bin' which is not on PATH.
    Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.
    WARNING: The scripts ase, ase-build, ase-db, ase-gui, ase-info and ase-run are installed in '/home/YOURUSERNAME/.local/bin' 
    which is not on PATH.
    Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.

You must add this folder to your ``~/.bashrc`` file.

    export PATH=/home/YOURUSERNAME/.local/bin:$PATH
    
Then continue,

    sudo apt install python3-dev libopenblas-dev libxc-dev libscalapack-mpi-dev libfftw3-dev

Create a siteconfig.py file:

```
$ mkdir -p ~/.gpaw
$ cat > ~/.gpaw/siteconfig.py
fftw = True
scalapack = True
libraries = ['xc', 'blas', 'fftw3', 'scalapack-openmpi']
^D
```
Then install gpaw

    pip3 install --upgrade --user gpaw

Use `gpaw info` to see information about installation. However, PAW-datasets are not installed yet. To install it, firstly create a directory under `~/.gpaw` then install PAW datasets

    mkdir ~/.gpaw/gpaw-setups
    gpaw install-data ~/.gpaw/gpaw-setups/

Now, you can test your GPAW with https://github.com/lrgresearch/gpaw-tools/blob/main/benchmarks/simple_benchmark_2021.py file

#### Installation of ASAP and KIM for Quick Optimization

For quick optimization, we need simple interatomic modelling. For this, we need ASAP3 for ASE, then we must install KIM with OpenKIM models and kimpy libraries.

    pip install --upgrade --user ase asap3
    sudo add-apt-repository ppa:openkim/latest
    sudo apt-get update
    sudo apt-get install libkim-api-dev openkim-models libkim-api2 pkg-config
    pip3 install kimpy

Then you can use files in https://github.com/lrgresearch/gpaw-tools/tree/main/QuickOptimize

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
Command line usage: `gpawsolve.py -o -c <configfile.py> -h -i <inputfile.cif>`

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
 * Change the parameters for each simulation OR change `config.py` once and use `-c` argument.
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
| PW-G0W0| Yes                 | Yes         | No             | No  | Yes  | No               | No      |
| PW-EXX*| Yes (with PBE)      | No          | No             | No  | No   | No               | No      |
|  LCAO  | No                  | No          | No             | Yes | Yes  | Yes              | No      |

*: Just some ground state energy calculations for PBE0 and HSE06.

### gg.py
Basic DFT calculations can be done graphically with the script `gg.py`. This script is behaving as a GUI to run `gpawsolve.py` script. To execute the GUI, type simply:
  python gg.py

### quick_optimize/quickoptimize.py
Inter-atomic potentials are useful tool to perform a quick geometric optimization of the studied system before starting a precise DFT calculation. The `quickoptimize.py` script is written for geometric optimizations with inter-atomic potentials. The bulk configuration of atoms can be provided by the user in the script as an ASE Atoms object or given as an argument for the CIF file. A general potential is given for any calculation. However, user can provide the necessary OpenKIM potentialby changing the related line in the script.

Mainly, quickoptimize.py is not related to GPAW. However it is dependent to ASAP3/OpenKIM and Kimpy. Therefore, the user must install necessary libraries before using the script.

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
