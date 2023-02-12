---
layout: default
nav_order: 4
title: Installation on WSL
parent: installation
---

# Installation on WSL

## Installation on a Windows 10 system (with WSL1)
This installation note is explaining how to install WSL1 and other required tools to study gpaw-tools on a Windows 10 system. WSLg is coming on default at Windows 11 and installation note for Windows 11 is explaining similar things but using WSLg. If you are using Windows 11, please continue from here.
We are suggesting Ubuntu 20.04 LTS version for *gpaw-tools* studies.

### WSL Version Control
[After activating the WSL](https://www.windowscentral.com/install-windows-subsystem-linux-windows-10), and installing an Ubuntu system for Microsoft Store application, we need to control the WSL version of Linux distribution installed in your system. Open Powershell and type:

    wsl --list --verbose 

Result will be something like:

    NAME            STATE           VERSION
    Ubuntu-20.04    Running         1

Your WSL distribution version must be "1", not "2". If it is not, change the WSL distribution version with command:

    wsl --set-version Ubuntu-20.04 1

### Update, upgrade and XWindows
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

## Installation on a Windows 11 system (with WSL2 and WSLg)
This installation note is explaining how to install WSLg and other required tools to study gpaw-tools on a Windows 11 system. WSL2 is coming on default at Windows 11. With WSLg you will not need to install xming to your windows system and Ubuntu packages libglu1-mesa-dev freeglut3 freeglut3-dev mesa-common-dev to your system.
We are suggesting Ubuntu 20.04 LTS version for *gpaw-tools* studies.

Open Ubuntu, finish the installation of it, then update your Linux system with:

    sudo apt update
    sudo apt upgrade
    
You do not need to install X server on your Windows to use with WSLg. 

### WSL2 memory problem
By default the WSL2 will consume up to 50% of the total system memory upto 8GB at max. However, it is possible to configure an upper limit for the memory and swap usage. Firstly, you must create a .wslconfig file in your Windows related home directory (C:\Users\<user>). And for 14GB RAM and 32 GB Swap, add following information to that file:

    [wsl2]
    memory=14GB
    swap=32GB

Then, open powershell or terminal and run:

    wsl --shutdown

Then, you can continue with new RAM-Swap settings.
    
## Installation of ASE and GPAW
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

## Installation of ASAP and KIM for Quick Optimization

For quick optimization, we need simple interatomic modelling. For this, we need ASAP3 for ASE, then we must install KIM with OpenKIM models and kimpy libraries.

    pip install --upgrade --user ase asap3
    sudo add-apt-repository ppa:openkim/latest
    sudo apt-get update
    sudo apt-get install libkim-api-dev openkim-models libkim-api2 pkg-config
    pip3 install kimpy

Then you can continue on [installation of gpaw-tools](installationofgpawtools.md).

