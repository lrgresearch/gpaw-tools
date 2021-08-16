from ase import Atoms
from ase.io import read, write
from ase.io.cif import read_cif
import sys, os
#from ase.io import read, write
# Converting CIF files to ASE Atoms object to use in GPAW Input for LRG Studies
# by Sefer Bora Lisesivdin
# August 2021 - First version
# Usage: $ python ciftoase.py file.cif
# -------------------------------------------------------------
# Parameters
# -------------------------------------------------------------
Scaled = False #Scaled or cartessian coordinates
Manualpbc = False # If you need manual constraint axis

# IF Manualpbc is true then change following:
pbcmanual=[True,True,False]

# If you have double number of elements in your final file
SolveDoubleElementProblem = True
# -------------------------------------------------------------
# ///////   YOU DO NOT NEED TO CHANGE ANYTHING BELOW    \\\\\\\ 
# -------------------------------------------------------------
inFile = sys.argv[1]
asestruct = read(inFile, index='-1')
print(inFile)
print(os.path.splitext(inFile)[0])
f = open(os.path.splitext(inFile)[0]+'.py', 'w')

# PRINT TO FILE PART -----------------------------------
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
f.close()

# PRINT SCREEEN PART -----------------------------------
print("Result:")
f = open(os.path.splitext(inFile)[0]+'.py', 'r')
lines = f.readline()
while lines:
    print(lines)
    lines = f.readline()
f.close()
