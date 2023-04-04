---
layout: default
nav_order: 6
title: Installation of gpaw-tools
parent: installation
---

# Installation of gpaw-tools

After installing ASE, GPAW, ASAP, KIM and other necessary packages, you can proceed with the installation of *gpaw-tools*.

Before, starting to installation of `gpaw-tools`, we need to install `spglib` and `elastic` packages and their dependencies. If you used conda to install previous packages, you do not need to run the following commands (If you try to run these two command as a single command you may have receive an error).

    pip3 install setuptools_scm
    pip3 install spglib docutils elastic
    pip3 install requests

If you want to use energy consumption measurement feauture, install:

    pip3 install pyrapl pymondo pandas

Also, lastly, it is good to use a job queue system when you have many inputs to run. GPAW / gpaw-tools can be run with task managers like SLURM. However, if you use your GPAW / gpaw-tools system on your local server/workstation, using Task Spooler is a good idea. It only works on one server for one user. It makes a queue, and run your commands in order. To install tsp command to your Ubuntu system, use:

    sudo apt install task-spooler

Now, all needed packages are installed and we can continue with installation of `gpaw-tools`. In your home folder (~), let's download the latest development release (you can prefer stable release also, please visit https://www.lrgresearch.org/gpaw-tools/ to get the latest URL)

    cd ~
    wget https://github.com/lrgresearch/gpaw-tools/archive/refs/heads/main.zip
    unzip main.zip

All files will be extracted to a folder called `gpaw-tools-main`. We need to make this folder to `~/.bashrc` file to system-wide reach.

    nano ~/.bashrc

Add the following line at the end of your ``~/.bashrc`` file.

    export PATH=/home/YOURUSERNAME/gpaw-tools-main:$PATH


After editing ~/.bashrc file quit the current shell session and start a new one (or you can use `source ~/.bashrc` command). 

Congratulations! You installed all necessary files to run *gpaw-tools*. You can continue with our [usage](generalusage.md) page, or continue with the `examples` folder in your `gpaw-tools-main` folder. All examples have README.md files.

