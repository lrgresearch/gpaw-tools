## Welcome to *gpaw-tools*

| Note: |
| ----- |
| These pages will be the future website of gpaw-tools. Please visit later. |

gpaw-tools is a bunch for python scripts for easy performing of GPAW calculations. It is mostly written for new DFT users who are running codes in their own PCs or on small group clusters.

gpaw-tools have:

A force-field quick optimization script quickoptimization.py for preliminary calculations using ASAP3/OpenKIM potentials.
ciftoase.py script for transform CIF files to ASE's own Atoms object.
To choose better cut off energy, lattice parameter and k points, there are 3 scripts called optimize_cutoff.py, optimize_latticeparam.py and optimize_kpoints.py.
And, the main solver script gpawsolver.py which can be run in PW or LCAO mode. It can do strain minimization, can use several different XCs, can do spin-polarized calculations, can calculate, draw and save tidily DOS and band structures, can calculate and save all-electron densities and can calculate optical properties in a very simple and organized way.
