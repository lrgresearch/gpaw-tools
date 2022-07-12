'''
optimize_latticeparam.py: Sample Lattice Paramater Optimization with GPAW

Usage: $ gpaw -P<core_number> python optimize_latticeparam.py geometry_file.cif

'''

import numpy as np
from ase import *
from ase.io import read
from ase.geometry.cell import cell_to_cellpar, cellpar_to_cell
from gpaw import GPAW, PW
import sys, getopt
from ase.io import Trajectory
from ase.parallel import paropen, world, parprint, broadcast
from argparse import ArgumentParser
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D
# -------------------------------------------------------------
# ENTER PARAMETERS
# -------------------------------------------------------------
mode = 'LCAO' # PW or LCAO
cut_off_energy = 340 # only for PW

kpts_x = 5
kpts_y = 5
kpts_z = 5

gpts_density = 0.2

xc_used = 'PBE'
Spin_calc = True
Hubbard = {}      # Example: {'O': ':p,7.0','Zn': ':d,12.0'}
Gamma = True


a0 = 3.1 # Creating a lattice parameter list
c0 = 5.05 # Needed for hexagonal structures
percent_a = 0.10 # changing the lattice parameter how much percent? (for lattice param a)
percent_c = 0.10 # changing the lattice parameter how much percent? (for lattice param c)

# If you want to optimize only on lattice param a or c, just make the other datapoint_count_x = 1
datapoint_count_a = 11 # how many data points for lattice param a?
datapoint_count_c = 11 # how many data points for lattice param c?

Draw_figure = False # Draw 3D figure for latt_a - latt_c optimization.

# --------------------------------------
# DO NOT CHANGE ANYTHING BELOW THIS LINE
# --------------------------------------
# Control draw figure is used correctly or not.
if (datapoint_count_a == 1 or datapoint_count_c == 1) and Draw_figure == True:
    parprint("ERROR: -d arg can not be used for single lattice param optimization.")
    quit()

# Populize lists
a_list = a0 * (1 + np.linspace(-percent_a , +percent_a , datapoint_count_a)) # do not change this line
c_list = c0 * (1 + np.linspace(-percent_c , +percent_c , datapoint_count_c)) # do not change this line

# -------------------------------------------------------------
# Bulk Configuration
# -------------------------------------------------------------
# Read bulk structure from CIF
if len(sys.argv) > 1:
    inFile = sys.argv[1]
    bulk_configuration = read(inFile, index='-1')

a, b, c, alpha, beta, gamma = cell_to_cellpar(bulk_configuration.get_cell(), radians=False)

# prepare arrays for total energies and correspoding volumes
etots = []
vols = []
# prepare trajectory object where all cells are stored in
traj = Trajectory("Optimize-Lattice_Trajectory.traj", "w")

# print nice header
parprint(" a      c      volume total energy")
parprint("----------------------------------")

with paropen('Optimize-Lattice_Table-LatticeParam.txt', 'w') as f:
    for latt_a in a_list:
        for latt_c in c_list:
            if datapoint_count_c == 1:
                latt_c = c
                bulk_configuration.set_cell(cellpar_to_cell([a*latt_a/a, b*latt_a/a, c, alpha, beta, gamma]), scale_atoms = True)
            else:
                bulk_configuration.set_cell(cellpar_to_cell([a*latt_a/a, b*latt_a/a, c*latt_c/c, alpha, beta, gamma]), scale_atoms = True)

            # --------------------------------------------------------------
            # create the cell, then find, store and print the total energy for a in a_list:
            if mode == "PW":
                bulk_configuration.calc = GPAW(mode=PW(cut_off_energy), xc=xc_used, nbands='200%', 
                            setups= Hubbard, parallel={'domain': world.size}, h = gpts_density,
                            spinpol=Spin_calc, kpts={'size': (kpts_x, kpts_y, kpts_z),'gamma': Gamma},
                            txt='Optimize-Lattice_%.2f_%.2f.txt' % (latt_a, latt_c))
            else:
                bulk_configuration.calc = GPAW(mode='lcao', parallel={'domain': world.size}, xc=xc_used,
                            setups= Hubbard, h= gpts_density, kpts={'size': (kpts_x, kpts_y, kpts_z),'gamma': Gamma},
                            basis='dzp', txt='Optimize-Lattice_%.2f_%.2f.txt' % (latt_a, latt_c))
            # we need volume and energy for E(V)−curve;
            # use corresponding getter functions and append values to lists 
            vol = bulk_configuration.get_volume()
            vols.append(vol)
            etot = bulk_configuration.get_potential_energy()
            etots.append(etot)
            # Writing to screen
            parprint("{:1.3f} {:1.3f}     {:2.3f} {:2.6f}".format(latt_a, latt_c, vol, etot))
            # Writing to text file
            f.write("{:1.3f} {:1.3f}     {:2.3f} {:2.6f}\n".format(latt_a, latt_c, vol, etot))

            # write the entire configuration (cell + energy) as trajectory
            traj.write(bulk_configuration)

if Draw_figure == True:
    if world.rank == 0:
        chunk = np.loadtxt('./Optimize-Lattice_Table-LatticeParam.txt')
        data=np.array(chunk)

        Xs = data[:,0]
        Ys = data[:,1]
        Zs = data[:,3]

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        surf = ax.plot_trisurf(Xs, Ys, Zs, cmap=cm.jet, linewidth=0)
        fig.colorbar(surf)

        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.yaxis.set_major_locator(MaxNLocator(6))
        ax.zaxis.set_major_locator(MaxNLocator(5))

        fig.tight_layout()
        plt.show()
        
