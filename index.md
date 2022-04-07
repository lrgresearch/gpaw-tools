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
3. To choose better cut off energy, lattice parameter and k points, there are 4 scripts called `optimize_cutoff.py`, `optimize_latticeparam.py` and `optimize_kpoints.py` and `optimize_kptsdensity.py`.
4. The main solver script `gpawsolver.py` which can be run in PW (also with GW and EXX) or LCAO mode. It can do structure optimization, Equation of State and elastic tensor calculations, can use several different XCs, can do spin-polarized calculations, can calculate, draw and save tidily DOS and band structures, can calculate and save all-electron densities and can calculate optical properties in a very simple and organized way.
5. A simple Graphical User Interface (GUI) for `gpawsolve.py` (and also you may say that GUI for GPAW) which is called `gg.py`.

More information about [gpaw-tools idea](about.md), [installation](installation.md), [usage](usage.md) and [release notes](releasenotes.md) can be found at related pages.

## Download

**Latest stable release: v22.4.0 [download (tar.gz)](https://github.com/lrgresearch/gpaw-tools/archive/refs/tags/v22.4.0.tar.gz), [download (zip)](https://github.com/lrgresearch/gpaw-tools/archive/refs/tags/v22.4.0.zip)**

Latest development release: [download (tar.gz)](https://github.com/lrgresearch/gpaw-tools/archive/refs/heads/main.tar.gz), [download (zip)](https://github.com/lrgresearch/gpaw-tools/archive/refs/heads/main.zip)

## News
* **[gpaw-tools](releasenotes.md#version-2240)** version 22.4.0 released (Apr 7, 2022).
* **[gpaw-tools](releasenotes.md#version-2230)** version 22.3.0 released (Mar 4, 2022).
* Our paper about *gpaw-tools* is published in Computational Material Science. [Download without any registration or fees](https://authors.elsevier.com/a/1ePsf3In-uvQ8R) before March 5th 2022.
* **[gpaw-tools](releasenotes.md#version-21120)** version 21.12.0 released (Dec 2, 2021).
* **[gpaw-tools](releasenotes.md#version-21110)** version 21.11.0 released (Nov 2, 2021).
* **[gpaw-tools](releasenotes.md#version-21101)** version 21.10.1 released (Oct 1, 2021).
* **[gpaw-tools](releasenotes.md#version-21100)** version 21.10.0 released (Oct 1, 2021).
* **[gpaw-tools](releasenotes.md#version-2190)** version 21.9.0 released (Sep 14, 2021).

## Citing
Please do not forget that, gpaw-tools is a UI/GUI software. For the main DFT calculations, it uses ASE and GPAW. It also uses Elastic python package for elastic tensor solutions and ASAP with KIM database for interatomic interaction calculations. Therefore, you must know what you use, and cite them properly. Here, the basic citation information of each packages are given. There are many other packages needed to be cited. With GPAW, you may needed to cite LibXC or cite for LCAO, TDDFT, lineer-response calculations. Please visit their pages for many other citation possibilities. 

* **ASE** : Ask Hjorth Larsen et al. "[The Atomic Simulation Environment—A Python library for working with atoms](https://doi.org/10.1088/1361-648X/aa680e)" J. Phys.: Condens. Matter Vol. 29 273002, 2017.
* **GPAW**: J. J. Mortensen, L. B. Hansen, and K. W. Jacobsen "[Real-space grid implementation of the projector augmented wave method](https://doi.org/10.1103/PhysRevB.71.035109)" Phys. Rev. B 71, 035109 (2005) and J. Enkovaara, C. Rostgaard, J. J. Mortensen et al. "[Electronic structure calculations with GPAW: a real-space implementation of the projector augmented-wave method](https://doi.org/10.1088/0953-8984/22/25/253202)" J. Phys.: Condens. Matter 22, 253202 (2010) [OTHER POSSIBLE CITATION](https://wiki.fysik.dtu.dk/gpaw/faq.html#citation-how-should-i-cite-gpaw)
* **KIM** : E. B. Tadmor, R. S. Elliott, J. P. Sethna, R. E. Miller and C. A. Becker "The Potential of Atomistic Simulations and the Knowledgebase of Interatomic Models" JOM, 63, 17 (2011). doi:10.1007/s11837-011-0102-6. [OTHER POSSIBLE CITATION](https://openkim.org/how-to-cite/)
* **Elastic**: P.T. Jochym, K. Parlinski and M. Sternik "[TiC lattice dynamics from ab initio calculations](https://doi.org/10.1007/s100510050823)", European Physical Journal B; 10, 9 (1999).

And for `gpaw-tools` usage, please use the following citation:

* S.B. Lisesivdin, B. Sarikavak-Lisesivdin "[gpaw-tools – higher-level user interaction scripts for GPAW calculations and interatomic potential based structure optimization](https://doi.org/10.1016/j.commatsci.2022.111201)" Comput. Mater. Sci. 204, 111201 (2022).

## Licensing
This project is licensed under the terms of the [MIT license](https://opensource.org/licenses/MIT).
