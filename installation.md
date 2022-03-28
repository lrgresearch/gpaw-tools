---
layout: default
navigation_weight: 3
title: Installation
---

# Installation

This installation notes are based on installations on Linux operating system. Here, we are using Ubuntu distro as Linux system. However, other linux distros will also work with similar commands. If you are new to Linux, please continue with Ubuntu Linux and our commands given below,

Required softwares for successfully running *gpaw-tools* software:

* Linux Operating system (can be independently installed, virtually installed or installed as WSL distro under Windows Operating system)
* Python 3.6 or above (with PIP and TK packages)
* Atomic Simulation Environment (ASE) (will automatically install to NumPy, SciPy, Matplotlib packages)
* GPAW
* BLAS, LibXC, MPI, LibKIM, OpenKIM and ScaLAPACK packages
* Some other Python packages like Elastic, spglib, kimpy, ASAP3.

There are two-possible installation methods:

1. conda installation
2. Manual installation

## Conda installation

The best and the easiest way to install ASE/GPAW/Elastic system with gpaw-tools is a conda installation. Download and install the miniconda. You can say ‘yes’ or ‘no’ to initialization after installing it:

    $ wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
    $ chmod +x Miniconda3-latest-Linux-x86_64.sh
    $ ./Miniconda3-latest-Linux-x86_64.sh

then you can update miniconda:

    $ eval "$(/home/$USER/miniconda3/bin/conda shell.bash hook)"
    $ conda update conda

Now, we can create an environment (here ‘gpaw-env’ name is used. You can use any name) and activate it:

    $ conda create --name gpaw-env
    $ conda activate gpaw-env

Then, install GPAW and Elastic packages

    $ conda install -c conda-forge gpaw elastic requests

Lastly, [download and install gpaw-tools](https://www.lrgresearch.org/gpaw-tools/installation/#4-installation-of-gpaw-tools).

## Manual installation

1. Choosing and preparing your Linux System,
2. Install ASE and GPAW,
3. Install ASAP3 and KIM,
4. Install Elastic package,
5. Install *gpaw-tools*.


### 1. Choosing and preparing your Linux System

Windows Subsystem for Linux on Windows, Virtual machine or a Linux machine?

You can use *gpaw-tools* on many different systems that are supporting required libraries. In this installation notes, we are using Linux systems. We would like to give some information on installation on Linux on Windows machines, virtual machines or a pure Linux machine. You can also use a Mac system, however because we do not have a Mac system to control it we couldn't give the installation notes for a Mac system. Please do not hesitate to send us if you have installation notes for a Mac system.

 * [Installation on a Windows 10 system](installation.md#installation-on-a-windows-10-system-with-wsl1)
 * [Installation on a Windows 11 system](installation.md#installation-on-a-windows-11-system-with-wslg)
 * [Installation on an independent Linux system](installation.md#installation-on-an-independent-linux-system)
 * [Installation on a Virtual Machine Linux system](installation.md#installation-on-a-virtual-machine-linux-system)

#### Running on Windows machines with WSL
The Windows Subsystem for Linux (WSL) allows developers to run a wide range of Linux-based apps and utilities on Windows without the need for a traditional virtual machine or a dualboot setup. More information about installation can be found [here](https://docs.microsoft.com/en-us/windows/wsl/install). There are two versions of WSL which are called WSL1 and WSL2. For GPAW calculations, WSL1 is giving better computation times than WSL2. However, with the announcement of WSLg on April 2021 at the Microsoft Build 2021, WSL2 is seemed to be the future of WSL.

##### Installation on a Windows 10 system (with WSL1)
This installation note is explaining how to install WSL1 and other required tools to study gpaw-tools on a Windows 10 system. WSLg is coming on default at Windows 11 and installation note for Windows 11 is explaining similar things but using WSLg. If you are using Windows 11, please continue from here.
We are suggesting Ubuntu 20.04 LTS version for *gpaw-tools* studies.

##### WSL Version Control
[After activating the WSL](https://www.windowscentral.com/install-windows-subsystem-linux-windows-10), and installing an Ubuntu system for Microsoft Store application, we need to control the WSL version of Linux distribution installed in your system. Open Powershell and type:

    wsl --list --verbose 

Result will be something like:

    NAME            STATE           VERSION
    Ubuntu-20.04    Running         1

Your WSL distribution version must be "1", not "2". If it is not, change the WSL distribution version with command:

    wsl --set-version Ubuntu-20.04 1

##### Update, upgrade and XWindows
Open Ubuntu, finish the installation of it, then update your Linux system with:

    sudo apt update
    sudo apt upgrade
    
To use a X server on windows, firstly you need to install few files (not necessary if you do not use matplotlib features.)

    sudo apt install libglu1-mesa-dev freeglut3 freeglut3-dev mesa-common-dev

You also need to install [XMing](https://sourceforge.net/projects/xming/) on Windows. In new WSL installations, the following line must be added to `~/.bashrc` file. You can use nano editor to do this:

    nano ~/.bashrc

and add these lines at the end of the file

    export LIBGL_ALWAYS_INDIRECT=0
    export DISPLAY=localhost:0.0

After editing ~/.bashrc file quit the current shell session and start a new one (or you can use `source ~/.bashrc` command).

##### Installation on a Windows 11 system (with WSLg)
This installation note is explaining how to install WSLg and other required tools to study gpaw-tools on a Windows 11 system. WSLg is coming on default at Windows 11.
We are suggesting Ubuntu 20.04 LTS version for *gpaw-tools* studies.

Open Ubuntu, finish the installation of it, then update your Linux system with:

    sudo apt update
    sudo apt upgrade
    
You do not need to install X server on your Windows to use with WSLg. 

#### Running on Linux machines
Instead of using a Windows system as a host, you can use Linux system to do your calculations. It is rather simple and it will give you more performance.

##### Installation on an independent Linux system
After [installing Ubuntu on your PC](https://ubuntu.com/tutorials/install-ubuntu-desktop#1-overview), you just need to update and upgrade, before continuing:

    sudo apt update
    sudo apt upgrade
    
##### Installation on a Virtual Machine Linux system
After [installing your Linux system inside Windows or other Linux using a virtualization software like VirtualBox](https://itsfoss.com/install-linux-in-virtualbox/), you just need to update and upgrade, before continuing:

    sudo apt update
    sudo apt upgrade
    
### 2. Installation of ASE and GPAW
After preparing your Linux system, you must have `ase` and `gpaw` codes on your computer. You can find more information about installation of [ASE](https://wiki.fysik.dtu.dk/ase/install.html) and [GPAW](https://wiki.fysik.dtu.dk/gpaw/install.html) from their related sites.

You need Tk library for GUI, unzip for file unzipping and for further package installations, we need PIP installer

    sudo apt install python3-tk python3-pip unzip python-is-python3

Install ASE and other math, parallel, dev libraries

    pip3 install --upgrade --user ase
    
At this point, PIP can give some warnings as:

    WARNING: The scripts f2py, f2py3 and f2py3.8 are installed in '/home/YOURUSERNAME/.local/bin' which is not on PATH.
    Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.
    WARNING: The scripts ase, ase-build, ase-db, ase-gui, ase-info and ase-run are installed in '/home/YOURUSERNAME/.local/bin' 
    which is not on PATH.
    Consider adding this directory to PATH or, if you prefer to suppress this warning, use --no-warn-script-location.

Add the following line at the end of your ``~/.bashrc`` file.

    export PATH=/home/YOURUSERNAME/.local/bin:$PATH
    

After editing ~/.bashrc file quit the current shell session and start a new one (or you can use `source ~/.bashrc` command). Then continue,

    sudo apt install python3-dev libopenblas-dev libxc-dev libscalapack-mpi-dev libfftw3-dev

Create a `siteconfig.py` file:

```
$ mkdir -p ~/.gpaw
$ cat > ~/.gpaw/siteconfig.py
fftw = True
scalapack = True
libraries = ['xc', 'blas', 'fftw3', 'scalapack-openmpi']
^D
```

NOTE: If the user wants to use exchange correlations listed in [libxc library](https://www.tddft.org/programs/libxc/), 'xc' must be listed in the libraries line as shown above.
{: .text-red-200 }


Then install gpaw

    pip3 install --upgrade --user gpaw

Use `gpaw info` to see information about installation. However, PAW-datasets are not installed yet. To install it, firstly create a directory under `~/.gpaw` then install PAW datasets

    mkdir ~/.gpaw/gpaw-setups
    gpaw install-data ~/.gpaw/gpaw-setups/

### 3. Installation of ASAP and KIM for Quick Optimization

For quick optimization, we need simple interatomic modelling. For this, we need ASAP3 for ASE, then we must install KIM with OpenKIM models and kimpy libraries.

    pip install --upgrade --user ase asap3
    sudo add-apt-repository ppa:openkim/latest
    sudo apt-get update
    sudo apt-get install libkim-api-dev openkim-models libkim-api2 pkg-config
    pip3 install kimpy

Then you can use files in https://github.com/lrgresearch/gpaw-tools/tree/main/QuickOptimize

### 4. Installation of gpaw-tools

Before, starting to installation of `gpaw-tools`, we need to install `spglib` 

    pip3 install spglib

and some packages to run elastic package, which is needed to run gpaw-tools: 

    pip3 install setuptools_scm
    pip3 install docutils
    pip3 install elastic

Now, all needed packages are installed and we can continue with installation of `gpaw-tools`. In your home folder (~), let's download the latest development release (you can prefer stable release also, please visit https://www.lrgresearch.org/gpaw-tools/ to get the latest URL)

    cd ~
    wget https://github.com/lrgresearch/gpaw-tools/archive/refs/heads/main.zip
    unzip main.zip

All files will be extracted to a folder called `gpaw-tools-main`. We need to make some files executable, and add this folder to `~/.bashrc` file to system-wide reach.

    cd gpaw-tools-main/
    chmod +x gg.py gpawsolve.py
    nano ~/.bashrc

Add the following line at the end of your ``~/.bashrc`` file.

    export PATH=/home/YOURUSERNAME/gpaw-tools-main:$PATH


After editing ~/.bashrc file quit the current shell session and start a new one (or you can use `source ~/.bashrc` command). 

Congratulations! You installed all necessary files to run *gpaw-tools*. You can continue with our [usage](usage.md) page, or continue with the `examples` folder in your `gpaw-tools-main` folder. All examples have README.md files.
