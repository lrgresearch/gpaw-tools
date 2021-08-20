import numpy as np
from ase import *
from gpaw import GPAW, PW
from ase.io import Trajectory
from ase.parallel import paropen, world, parprint

# Sample Lattice Parameter Optimization GPAW Input for LRG Studies
# by Sefer Bora Lisesivdin
# August 2021 - Rewritten with using Prof. J. Kortus, R. Wirnata - WS 2019 Course notes
# July 2021 - Corrected version
# March 2020 - First Version 
# Usage: $ gpaw python 0-2-Optimize-0-Lattice.py
# 
# -------------------------------------------------------------
# ENTER PARAMETERS
# -------------------------------------------------------------
cutoff = 400

kpts_x = 3
kpts_y = 3
kpts_z = 1

a0 = 3.7 # Creating a lattice parameter list
percent = 0.05 # changing the lattice parameter how much percent?

a_list = a0 * (1 + np.linspace(-percent , +percent , 11)) # do not change this line

# -------------------------------------------------------------
# Bulk Configuration
# -------------------------------------------------------------
## USE latt_a for lattice parameter 

# prepare arrays for total energies and correspoding volumes 
etots = []
vols = []
# prepare trajectory object where all cells are stored in
traj = Trajectory("0-2-Optimize-Lattice-Trajectory.traj", "w")

# print nice header
parprint(" a      volume total energy")
parprint("-----------------------------")

f = paropen('Table-LatticeParam.txt', 'w')

f.write(' a(Ang)   volume total energy\n')
f.write('-----------------------------\n')
for latt_a in a_list:
# -------------------------------------------------------------
# DO NOT FORGET TO INSERT BULK CONFIGURATION UNDER HERE
# WITH PROPER INDENTATION!!!
# -------------------------------------------------------------
    bulk_configuration = Atoms(
        [
        Atom('Sn', ( 0.0, 0.0, 5.02 )),
        Atom('C', ( 0.0, 2.3094022314590417*latt_a/4.0, 4.98 )),
        ],
        cell=[(4.0, 0.0, 0.0), (-1.9999999999999991*latt_a/4.0, 3.464101615137755*latt_a/4.0, 0.0), (0.0, 0.0, 20.0)],
        pbc=True,
        )

# --------------------------------------------------------------
    
    # create the cell, then find, store and print the total energy for a in a_list:
    bulk_configuration.calc = GPAW(mode=PW(cutoff), 
                                parallel={'domain': world.size},
                                xc='PBE',
                                kpts=(kpts_x, kpts_y, kpts_z),
                                basis='dzp',
                                txt='Optimize-0-Lattice-%.2f.txt' % latt_a)
    # we need volume and energy for E(V)−curve;
    # use corresponding getter functions and append values to lists 
    vol = bulk_configuration.get_volume()
    vols.append(vol)
    etot = bulk_configuration.get_potential_energy()
    etots.append(etot)
    parprint("{:1.3f}      {:2.3f} {:2.6f}".format(latt_a, vol, etot))
    f.write("{:1.3f}      {:2.3f} {:2.6f}\n".format(latt_a, vol, etot))

    # write the entire configuration (cell + energy) as trajectory
    traj.write(bulk_configuration)

f.close()
