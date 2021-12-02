'''
quickoptimize.py: Quick Geometric Optimization
                  using interatomic potentials.
Usage: $ python quickoptimize.py geometryfile.cif
'''
import sys, os
from pathlib import Path
from ase import *
from ase.io import read
from ase.io.cif import write_cif
from asap3 import Atoms, units
from asap3.md.langevin import Langevin
from ase.calculators.kim import KIM


# -------------------------------------------------------------
# Parameters
# -------------------------------------------------------------

# Simulation parameters

#More potantial can be found from https://openkim.org
PotentialUsed = 'LJ_ElliottAkerson_2015_Universal__MO_959249795837_003'
Temperature = 1 #Kelvin
Time = 5 # fs
Friction = 0.05

Scaled = False #Scaled or cartessian coordinates
Manualpbc = False # If you need manual constraint axis

# IF Manualpbc is true then change following:
pbcmanual=[True,True,False]

# If you have double number of elements in your final file
SolveDoubleElementProblem = True

# If you do not want to use a CIF file for geometry, please provide
# ASE Atoms object information below. You can use ciftoase.py to
# make your own ASE Atoms object from a CIF file.
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
asestruct = bulk_configuration

# Basic argument triage
if len(sys.argv) == 1:
    Filename = 'Optimized'
    print("Number of atoms provided in Atoms object:"+str(asestruct.get_global_number_of_atoms()))
elif len(sys.argv) == 2:
    inFile = sys.argv[1]
    Filename = os.path.splitext(inFile)[0]
    struct = Path(inFile).stem
    asestruct = read(inFile, index='-1')
    print("Number of atoms imported from CIF file:"+str(asestruct.get_global_number_of_atoms()))
else:
    print("Wrong number of argument. Usage: python quickoptimize.py inputfile.cif")
    quit()

asestruct.set_calculator(KIM(PotentialUsed, options={"ase_neigh": False}))

dyn = Langevin(asestruct, timestep=Time*units.fs, trajectory='quickoptim.traj', logfile=Filename+'-Log-GPAW-QuickOptim.txt', temperature_K=Temperature, friction=Friction)


print("")
print("Energy per atom:")
print("  %15s %15s %15s" % ("Pot. energy", "Kin. energy", "Total energy"))

for i in range(25):
    dyn.run(10)
    epot = asestruct.get_potential_energy()/len(asestruct)
    ekin = asestruct.get_kinetic_energy()/len(asestruct)
    print("%15.5f %15.5f %15.5f" % (epot, ekin, epot+ekin))

# PRINT TO FILE PART -----------------------------------
with open(Filename+'.py', 'w') as f:
    f.write("bulk_configuration = Atoms(\n")
    f.write("    [\n")
    if Scaled == True:
        positions = asestruct.get_scaled_positions()
    else:
        positions = asestruct.get_positions()
    nn=0
    mm=0

    if SolveDoubleElementProblem == True:
        for n in asestruct.get_chemical_symbols():
            nn=nn+1
            for m in positions:
                mm=mm+1
                if mm == nn:
                    f.write("    Atom('"+n+"', ( "+str(m[0])+", "+str(m[1])+", "+str(m[2])+" )),\n")
            mm=0
    else:
        for n in asestruct.get_chemical_symbols():
            for m in positions:
                f.write("    Atom('"+n+"', ( "+str(m[0])+", "+str(m[1])+", "+str(m[2])+" )),\n")
    f.write("    ],\n")
    f.write("    cell=[("+str(asestruct.cell[0,0])+", "+str(asestruct.cell[0,1])+", "+str(asestruct.cell[0,2])+"), ("+str(asestruct.cell[1,0])+", "+str(asestruct.cell[1,1])+", "+str(asestruct.cell[1,2])+"), ("+str(asestruct.cell[2,0])+", "+str(asestruct.cell[2,1])+", "+str(asestruct.cell[2,2])+")],\n")
    if Manualpbc == False:
        f.write("    pbc=True,\n")
    else:
        f.write("    pbc=["+str(pbcmanual[0])+","+str(pbcmanual[1])+","+str(pbcmanual[2])+"],\n")
    f.write("    )\n")

# PRINT SCREEEN PART -----------------------------------
print("Result:")
with open(Filename+'.py', 'r') as f:
    lines = f.readline()
    while lines:
        print(lines)
        lines = f.readline()



print("    )")
print("ATTENTION: If you have double number of atoms, it may be caused by ")
print("           repeating ASE bug https://gitlab.com/ase/ase/-/issues/169 ")
print("           Please assign SolveDoubleElementProblem variable as True in this script if necessary.")


# CIF export
write_cif(Filename+'-Optimized.cif', bulk_configuration)
