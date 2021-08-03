import numpy as np
from ase import *
from gpaw import GPAW, PW
# Sample K-point Optimization GPAW Input for LRG Studies
# by Sefer Bora Lisesivdin
# July 2021 - Corrected version
# March 2020 - First Version 
# Usage: $ gpaw python script.py
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
f = open('Table-K_points.txt', 'a')
f.write('K-point  Total_Energy\n')
for k in range(kpoint_min, kpoint_max):
    bulk_configuration.calc = GPAW(mode=PW(cutoffenergy),
                              xc='PBE',
                              kpts=(k, k, k),
                              parallel={'band': 1},
                              basis='dzp',
                              txt='Optimize-0-kpoints-%02d.txt' % k)
    print ("K-point:"+str(k))
#    for eps in np.linspace(-0.02, 0.02, 5):
#        bulk_configuration.cell = (1 + eps) * cell0
#        bulk_configuration.get_potential_energy()
#        print ("eps:"+str(eps))
    f.write(str(k)+'  '+str(bulk_configuration.get_potential_energy())+'\n')
f.close()

