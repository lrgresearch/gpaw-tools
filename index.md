---
layout: default
navigation_weight: 1
title: Home
---

# Welcome to *gpaw-tools*
{: .fs-9 }

*gpaw-tools* is a collection of python scripts that use ASE and GPAW for performing Density Functional Theory (DFT) calculations. Its aim is to lower the entry barrier and to provide an easy-to-use command line and graphical user interfaces for GPAW. It is mostly written for new DFT users who are running codes on their own PCs or on small group clusters.
{: .fs-6 .fw-300 }

[Download now](#download){: .btn .btn-primary .fs-5 .mb-4 .mb-md-0 .mr-2 } [View it on GitHub](https://github.com/lrgresearch/gpaw-tools){: .btn .fs-5 .mb-4 .mb-md-0 }

`gpaw-tools` have:
1. A force-field quick optimization script `quickoptimization.py` for preliminary calculations using ASAP3/OpenKIM potentials. 
2. `ciftoase.py` script for transform CIF files to ASE's own Atoms object.
3. To choose better cut off energy, lattice parameter and k points, there are 3 scripts called `optimize_cutoff.py`, `optimize_latticeparam.py` and `optimize_kpoints.py`.
4. The main solver script `gpawsolver.py` which can be run in PW (also with GW and EXX) or LCAO mode. It can do structure optimization, can use several different XCs, can do spin-polarized calculations, can calculate, draw and save tidily DOS and band structures, can calculate and save all-electron densities and can calculate optical properties in a very simple and organized way.
5. A simple Graphical User Interface (GUI) for `gpawsolve.py` (and also you may say that GUI for GPAW) which is called `gg.py`.

More information about [gpaw-tools idea](about.md), [installation](installation.md), [usage](usage.md) and [release notes](releasenotes.md) can be found at related pages.

## Download

**Latest stable release: v21.12.0 [download (tar.gz)](https://github.com/lrgresearch/gpaw-tools/archive/refs/tags/v21.12.0.tar.gz), [download (zip)](https://github.com/lrgresearch/gpaw-tools/archive/refs/tags/v21.12.0.zip)**

Latest development release: [download (tar.gz)](https://github.com/lrgresearch/gpaw-tools/archive/refs/heads/main.tar.gz), [download (zip)](https://github.com/lrgresearch/gpaw-tools/archive/refs/heads/main.zip)

## News
* Our paper about *gpaw-tools* is published in Computational Material Science. [Download without any registration or fees](https://authors.elsevier.com/a/1ePsf3In-uvQ8R) before March 5th 2022.
* **[gpaw-tools](releasenotes.md#version-21120)** version 21.12.0 released (Dec 2, 2021).
* **[gpaw-tools](releasenotes.md#version-21110)** version 21.11.0 released (Nov 2, 2021).
* **[gpaw-tools](releasenotes.md#version-21101)** version 21.10.1 released (Oct 1, 2021).
* **[gpaw-tools](releasenotes.md#version-21100)** version 21.10.0 released (Oct 1, 2021).
* **[gpaw-tools](releasenotes.md#version-2190)** version 21.9.0 released (Sep 14, 2021).

## Licensing
This project is licensed under the terms of the [MIT license](https://opensource.org/licenses/MIT).
