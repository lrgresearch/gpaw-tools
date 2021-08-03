from ase import *
from ase.build import bulk
from ase.io import read, write
from ase.io.cif import write_cif
from gpaw import GPAW, PW
from ase.optimize.lbfgs import LBFGS
from ase.constraints import UnitCellFilter
import sys, os

# Quick Geometric Optimization for CIF Files for LRG Studies
# by Sefer Bora Lisesivdin
# August 2021 - First version
# Usage: $ gpaw python quickoptimize.py ciffilename.cif
# -------------------------------------------------------------
# Parameters
# -------------------------------------------------------------

Filename = 'Optimized'

Scaled = False #Scaled or cartessian coordinates
Manualpbc = False # If you need manual constraint axis

# IF Manualpbc is true then change following:
pbcmanual=[True,True,False]
# Which components of strain will be relaxed
# EpsX, EpsY, EpsZ, ShearYZ, ShearXZ, ShearXY
# Example: For a x-y 2D nanosheet only first 2 component will be true
whichstrain=[True, True, False, False, False, False]

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

f = open(Filename+'.py', 'w')

calc = GPAW(mode=PW(340), kpts=[2, 2, 1], xc='PBE', txt=Filename+'-Log-GPAW-QuickOptim.txt')
bulk_configuration.calc = calc
uf = UnitCellFilter(bulk_configuration, mask=whichstrain)
relax = LBFGS(uf, trajectory=Filename+'-Trajectory-GPAW-QuickOptim.traj')
relax.run(fmax=0.005)  # Consider much tighter fmax!

# PRINT TO FILE PART -----------------------------------
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
f.close()

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