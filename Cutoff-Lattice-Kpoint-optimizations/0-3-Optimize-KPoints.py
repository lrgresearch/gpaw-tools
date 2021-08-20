import numpy as np
from ase import *
from gpaw import GPAW, PW
from ase.parallel import paropen, world, parprint

# Sample K-point Optimization GPAW Input for LRG Studies
# by Sefer Bora Lisesivdin
# August 2021 - Better parallel computation and its handling
# July 2021 - Corrected version
# March 2020 - First Version 
# Usage: $ gpaw -P<core_number> python Optimize_KPoints.py
# 
# -------------------------------------------------------------
# ENTER PARAMETERS
# -------------------------------------------------------------
cutoffenergy = 200
kpoint_min = 4
kpoint_max = 11


# -------------------------------------------------------------
# Bulk Configuration
# -------------------------------------------------------------
bulk_configuration = Atoms(
    [
    Atom('Si', ( 0.0, 0.0, 0.0 )),
    Atom('Si', ( 1.35765, 1.35765, 1.35765 )),
    ],
    cell=[(0.0, 2.7153, 2.7153), (2.7153, 0.0, 2.7153), (2.7153, 2.7153, 0.0)],
    pbc=True,
    )
# -------------------------------------------------------------
# CONVERGE KPOINTS
# -------------------------------------------------------------
cell0 = bulk_configuration.cell
f = paropen('Optimize-KPoints_Table-KPoints.txt', 'a')
f.write('K-point  Total_Energy\n')
for k in range(kpoint_min, kpoint_max):
    bulk_configuration.calc = GPAW(mode=PW(cutoffenergy),
                              xc='PBE',
                              kpts=(k, k, k),
                              parallel={'band': 1},
                              basis='dzp',
                              txt='Optimize-KPoints_KPoints-%02d.txt' % k)
    parprint("K-point:"+str(k)+"   Potential Energy:"+str(bulk_configuration.get_potential_energy()))
    f.write(str(k)+'  '+str(bulk_configuration.get_potential_energy())+'\n')
f.close()

