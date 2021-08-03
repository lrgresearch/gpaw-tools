import numpy as np
from ase import *
from gpaw import GPAW, PW
# Sample Lattice Parameter Optimization GPAW Input for LRG Studies
# by Sefer Bora Lisesivdin
# July 2021 - Corrected version
# March 2020 - First Version 
# Usage: $ gpaw python script.py
# 
# -------------------------------------------------------------
# ENTER PARAMETERS
# -------------------------------------------------------------
cutoff = 200

kpts_x = 7
kpts_y = 7
kpts_z = 7

a0 = 2.7153 # Creating a lattice parameter list
percent = 0.05 # changing the lattice parameter how much percent?

a_list = a0 * (1 + np.linspace(-percent , +percent , 11)) # do not change this line
latt_a = a0
# -------------------------------------------------------------
# Bulk Configuration
# -------------------------------------------------------------
## USE latt_a for lattice parameter 
bulk_configuration = Atoms(
    [
    Atom('Si', ( 0.0, 0.0, 0.0 )),
    Atom('Si', ( 1.35765, 1.35765, 1.35765 )),
    ],
    cell=[(0.0, latt_a, latt_a), (latt_a, 0.0, latt_a), (latt_a, latt_a, 0.0)],
    pbc=True,
    )
# -------------------------------------------------------------
# CONVERGE CUTOFF
# -------------------------------------------------------------
cell0 = bulk_configuration.cell
f = open('Table-LatticeParam.txt', 'a')
f.write('Lattice_Parameter  Total_Energy\n')
for latt_a in a_list:
    bulk_configuration.calc = GPAW(mode=PW(cutoff),
                              xc='PBE',
                              kpts=(kpts_x, kpts_y, kpts_z),
                              parallel={'band': 1},
                              basis='dzp',
                              txt='Optimize-0-Lattice-%.2f.txt' % latt_a)
    print ("Lattice Parameter:"+str(latt_a))
#    for eps in np.linspace(-0.02, 0.02, 5):
#        bulk_configuration.cell = (1 + eps) * cell0
#        bulk_configuration.get_potential_energy()
#        print ("eps:"+str(eps))
    f.write(str(latt_a)+'  '+str(bulk_configuration.get_potential_energy())+'\n')
f.close()



