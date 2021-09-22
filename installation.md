---
layout: default
navigation_weight: 3
title: Installation
---

# Installation
 You must have `ase` and `gpaw` codes on your computer. You can find more information about installation of [ASE](https://wiki.fysik.dtu.dk/ase/install.html) and [GPAW](https://wiki.fysik.dtu.dk/gpaw/install.html) from their related sites.
 However, for a simple installation on a Windows 10 computer with WSL, a small installation information is given below:

## Full Installation of WSL, ASE, GPAW and gpaw-tools on a Windows 10 system

### Windows Subsystem for Linux (WSL)
The Windows Subsystem for Linux allows developers to run a wide range of Linux-based apps and utilities on Windows without the need for a traditional virtual machine or a dualboot setup. More information about installation can be found [here](https://docs.microsoft.com/en-us/windows/wsl/install).

We are suggesting Ubuntu 20.04 LTS version for gpaw-tools studies.

### WSL Linux Distribution Version Control
After installing WSL, we need to control the WSL version of Linux distribution installed in your system. Open Powershell and type:

    wsl --list --verbose 

Result will be something like:

    NAME            STATE           VERSION
    Ubuntu-20.04    Running         1

Your WSL distribution version must be "1", not "2". If it is not, change the WSL distribution version with command:

    wsl --set-version Ubuntu-20.04 1

### ASE/GPAW Installation
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

### Installation of ASAP and KIM for Quick Optimization

For quick optimization, we need simple interatomic modelling. For this, we need ASAP3 for ASE, then we must install KIM with OpenKIM models and kimpy libraries.

    pip install --upgrade --user ase asap3
    sudo add-apt-repository ppa:openkim/latest
    sudo apt-get update
    sudo apt-get install libkim-api-dev openkim-models libkim-api2 pkg-config
    pip3 install kimpy

Then you can use files in https://github.com/lrgresearch/gpaw-tools/tree/main/QuickOptimize
