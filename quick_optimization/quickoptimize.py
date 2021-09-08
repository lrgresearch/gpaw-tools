'''
quickoptimize.py: Quick Geometric Optimization
                  using interatomic potentials.
Usage: $ gpaw python quickoptimize.py ciffilename.cif
'''
import sys
from ase import *
from ase.io import read
from ase.io.cif import write_cif
from asap3 import Atoms, units
from asap3.md.verlet import VelocityVerlet
from asap3.md.langevin import Langevin
from ase.calculators.kim import KIM


# -------------------------------------------------------------
# Parameters
# -------------------------------------------------------------

Filename = 'Optimized'

Scaled = False #Scaled or cartessian coordinates
Manualpbc = False # If you need manual constraint axis

# IF Manualpbc is true then change following:
pbcmanual=[True,True,False]

WantCIFexport = True

# -------------------------------------------------------------
# Bulk Configuration
# -------------------------------------------------------------

bulk_configuration = Atoms(
    [
    Atom('Ge', ( 1.222474e-31, 4.094533468076675e-32, 5.02 )),
    Atom('Ge', ( -1.9999999993913775e-06, 2.3094022314590417, 4.98 )),
    ],
    cell=[(4.0, 0.0, 0.0), (-1.9999999999999991, 3.464101615137755, 0.0), (0.0, 0.0, 20.0)],
    pbc=True,
    )


# -------------------------------------------------------------
# ///////   YOU DO NOT NEED TO CHANGE ANYTHING BELOW    \\\\\\\
# -------------------------------------------------------------

bulk_configuration.set_calculator(KIM('LJ_ElliottAkerson_2015_Universal__MO_959249795837_003', options={"ase_neigh": False}))

# dyn = VelocityVerlet(bulk_configuration, 5.0 * units.fs, trajectory='md.traj', logfile=Filename+'-Log-GPAW-QuickOptim.txt')
dyn = Langevin(bulk_configuration, timestep=5*units.fs, trajectory='md.traj', logfile=Filename+'-Log-GPAW-QuickOptim.txt', temperature=1*units.kB, friction=0.05)


print("")
print("Energy per atom:")
print("  %15s %15s %15s" % ("Pot. energy", "Kin. energy", "Total energy"))

for i in range(25):
    dyn.run(10)
    epot = bulk_configuration.get_potential_energy()/len(bulk_configuration)
    ekin = bulk_configuration.get_kinetic_energy()/len(bulk_configuration)
    print("%15.5f %15.5f %15.5f" % (epot, ekin, epot+ekin))

# PRINT TO FILE PART -----------------------------------
with open(Filename+'.py', 'w') as f:
    f.write("bulk_configuration = Atoms(\n")
    f.write("    [\n")
    if Scaled == True:
        positions = bulk_configuration.get_scaled_positions()
    else:
        positions = bulk_configuration.get_positions()
    nn=-1
    mm=-1
    for n in bulk_configuration.get_chemical_symbols():
        nn=nn+1
        for m in positions:
            mm=mm+1
            if mm == nn:
                if (mm % 2) != 0:
                    f.write("    Atom('"+n+"', ( "+str(m[0])+", "+str(m[1])+", "+str(m[2])+" )),\n")
        mm=0
    f.write("    ],\n")
    f.write("    cell=[("+str(bulk_configuration.cell[0,0])+", "+str(bulk_configuration.cell[0,1])+", "+str(bulk_configuration.cell[0,2])+"), ("+str(bulk_configuration.cell[1,0])+", "+str(bulk_configuration.cell[1,1])+", "+str(bulk_configuration.cell[1,2])+"), ("+str(bulk_configuration.cell[2,0])+", "+str(bulk_configuration.cell[2,1])+", "+str(bulk_configuration.cell[2,2])+")],\n")
    if Manualpbc == False:
        f.write("    pbc=True,\n")
    else:
        f.write("    pbc=["+str(pbcmanual[0])+","+str(pbcmanual[1])+","+str(pbcmanual[2])+"],\n")
    f.write("    )\n")


# PRINT SCREEEN PART -----------------------------------
print("bulk_configuration = Atoms(")
print("    [")
if Scaled == True:
    positions = bulk_configuration.get_scaled_positions()
else:
    positions = bulk_configuration.get_positions()
nn=-1
mm=-1
for n in bulk_configuration.get_chemical_symbols():
    nn=nn+1
    for m in positions:
        mm=mm+1
        if mm == nn:
            if (mm % 2) != 0:
                print("    Atom('"+n+"', ( "+str(m[0])+", "+str(m[1])+", "+str(m[2])+" )),")
    mm=0
print("    ],")
print("    cell=[("+str(bulk_configuration.cell[0,0])+", "+str(bulk_configuration.cell[0,1])+", "+str(bulk_configuration.cell[0,2])+"), ("+str(bulk_configuration.cell[1,0])+", "+str(bulk_configuration.cell[1,1])+", "+str(bulk_configuration.cell[1,2])+"), ("+str(bulk_configuration.cell[2,0])+", "+str(bulk_configuration.cell[2,1])+", "+str(bulk_configuration.cell[2,2])+")],")
if Manualpbc == False:
    print("    pbc=True,")
else:
    print("    pbc=["+str(pbcmanual[0])+","+str(pbcmanual[1])+","+str(pbcmanual[2])+"],")
print("    )")

# CIF export
if WantCIFexport == True:
    write_cif(Filename+'.cif', bulk_configuration)
