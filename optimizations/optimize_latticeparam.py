'''
optimize_latticeparam.py: Sample Lattice Paramater Optimization with GPAW

Usage: $ gpaw -P<core_number> python optimize_latticeparam.py
'''

import numpy as np
from ase import *
from gpaw import GPAW, PW
from ase.io import Trajectory
from ase.parallel import paropen, world, parprint

# -------------------------------------------------------------
# ENTER PARAMETERS
# -------------------------------------------------------------
cutoff = 700

kpts_x = 5
kpts_y = 5
kpts_z = 5

xc_used = 'GLLBSCM'

kpts_density =2.5

use_density = True # Change to "True", if you want to use k-point density instead of kpoints.

a0 = 2.53 # Creating a lattice parameter list
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
traj = Trajectory("Optimize-Lattice_Trajectory.traj", "w")

# print nice header
parprint(" a      volume total energy")
parprint("-----------------------------")

with paropen('Optimize-Lattice_Table-LatticeParam.txt', 'w') as f:
    f.write('LattParam volume total_energy\n')
    f.write('-----------------------------\n')
    for latt_a in a_list:
        # -------------------------------------------------------------
        # DO NOT FORGET TO INSERT BULK CONFIGURATION UNDER HERE
        # WITH PROPER INDENTATION!!!
        # -------------------------------------------------------------           
        bulk_configuration = Atoms(['C' for i in range(2)],
              [(   3.790492*latt_a/2.526995,    2.188441*latt_a/2.526995,    1.547462*latt_a/2.526995),
               (   2.526995*latt_a/2.526995,    1.458961*latt_a/2.526995,    1.031641*latt_a/2.526995)],
              pbc = (True,True,True))
        bulk_configuration.set_cell([[    2.526995*latt_a/2.526995,     0.000000,     0.000000],
                [    1.263497*latt_a/2.526995,     2.188441*latt_a/2.526995,     0.000000],
                [    1.263497*latt_a/2.526995,     0.729480*latt_a/2.526995,     2.063282*latt_a/2.526995]],
                scale_atoms = False)

        # --------------------------------------------------------------
        # create the cell, then find, store and print the total energy for a in a_list:
        if use_density == True:
            bulk_configuration.calc = GPAW(mode=PW(cutoff),
                                        parallel={'domain': world.size},
                                        xc=xc_used,
                                        kpts={'density': kpts_density, 'gamma': True},
                                        basis='dzp',
                                        txt='Optimize-Lattice_%.2f.txt' % latt_a)
        else:
            bulk_configuration.calc = GPAW(mode=PW(cutoff),
                                        parallel={'domain': world.size},
                                        xc=xc_used,
                                        kpts=(kpts_x, kpts_y, kpts_z),
                                        basis='dzp',
                                        txt='Optimize-Lattice_%.2f.txt' % latt_a)
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
