'''
gpaw-tools: Simple Benchmark Calculation for GPAW
Usage: 
  $ gpaw -P8 python simple_benchmark_2023.py
For AMD CPUs or using Intel CPUs without hyperthreading: (Example CPU is intel here, 4 cores or 8 threads)
  $ mpirun -n 4 gpaw python simple_benchmark_2023.py
For using all threads provided by Intel Hyperthreading technology:
  $ mpirun --use-hwthread-cpus -n 8 gpaw python simple_benchmark_2023.py
'''

from ase import *
from gpaw import GPAW, PW
from ase.optimize.lbfgs import LBFGS
from ase.dft.dos import DOS
from pathlib import Path

# -------------------------------------------------------------
# Parameters
# -------------------------------------------------------------
fmaxval = 0.05 			#
cut_off_energy = 240 	# eV
kpts_x = 1 			# kpoints in x direction
kpts_y = 1				# kpoints in y direction
kpts_z = 20				# kpoints in z direction
band_path = 'GZ'	# Brillouin zone high symmetry points
band_npoints = 40		# Number of points between high symmetry points 
energy_max = 15 		# eV. It is the maximum energy value for band structure figure.
num_of_bands = 40		#
# -------------------------------------------------------------
# Bulk Configuration
# -------------------------------------------------------------
bulk_configuration = Atoms(
    [
    Atom('C', ( 5.0, 4.999995597, 2.84172142086 )),
    Atom('C', ( 5.0, 6.230496855, 3.55214857914 )),
    Atom('C', ( 5.0, 4.999995597, 1.42085857914 )),
    Atom('C', ( 5.0, 7.460998113, 2.84172142086 )),
    Atom('C', ( 5.0, 8.691499370999999, 3.55214857914 )),
    Atom('C', ( 5.0, 6.230496855, 0.71043142086 )),
    Atom('C', ( 5.0, 7.460998113, 1.42085857914 )),
    Atom('C', ( 5.0, 9.922000629, 2.84172142086 )),
    Atom('C', ( 5.0, 11.152501886999998, 3.55214857914 )),
    Atom('C', ( 5.0, 8.691499370999999, 0.71043142086 )),
    Atom('C', ( 5.0, 9.922000629, 1.42085857914 )),
    Atom('C', ( 5.0, 12.383003145, 2.84172142086 )),
    Atom('C', ( 5.0, 13.613504402999999, 3.55214857914 )),
    Atom('C', ( 5.0, 11.152501886999998, 0.71043142086 )),
    Atom('C', ( 5.0, 12.383003145, 1.42085857914 )),
    Atom('C', ( 5.0, 13.613504402999999, 0.71043142086 )),
    Atom('H', ( 5.0, 4.056030558, 3.3867178493399996 )),
    Atom('H', ( 5.0, 4.056030558, 0.87586215066 )),
    Atom('H', ( 5.0, 14.557469441999999, 3.00715215066 )),
    Atom('H', ( 5.0, 14.557469441999999, 1.25542784934 )),
    ],
    cell=[(10.0, 0.0, 0.0), (0.0, 18.6135, 0.0), (0.0, 0.0, 4.26258)],
    pbc=True,
    )

# -------------------------------------------------------------
# ///////   YOU DO NOT NEED TO CHANGE ANYTHING BELOW    \\\\\\\ 
# -------------------------------------------------------------
struct = Path(__file__).stem # All files will get their names from this file
# -------------------------------------------------------------
# Step 1 - GROUND STATE
# -------------------------------------------------------------
calc = GPAW(mode=PW(cut_off_energy), kpts=[kpts_x, kpts_y, kpts_z], txt=struct+'-GROUND-Log.txt')
bulk_configuration.calc = calc

relax = LBFGS(bulk_configuration, trajectory=struct+'-GROUND-Result.traj')
relax.run(fmax=fmaxval)  # Consider much tighter fmax!

bulk_configuration.get_potential_energy()
calc.write(struct+'-GROUND-Result.gpw')

# -------------------------------------------------------------
# Step 2 - DOS CALCULATION
# -------------------------------------------------------------
calc = GPAW(struct+'-GROUND-Result.gpw', txt=struct+'-DOS-Log.txt')
#energies, weights = calc.get_dos(npts=800, width=0)
dos = DOS(calc, npts=500, width=0)
energies = dos.get_energies()
weights = dos.get_dos()

fd = open(struct+'-DOS-Result.txt', "w")
for x in zip(energies, weights):
    print(*x, sep=", ", file=fd)
fd.close()

# -------------------------------------------------------------
# Step 3 - BAND STRUCTURE CALCULATION
# -------------------------------------------------------------
calc = GPAW(struct+'-GROUND-Result.gpw',
	    txt=struct+'-BAND-Log.txt',
	    nbands=num_of_bands,
	    fixdensity=True,
	    symmetry='off',
	    kpts={'path': band_path, 'npoints': band_npoints},
	    convergence={'bands': 8})

calc.get_potential_energy()
bs = calc.band_structure()
ef = calc.get_fermi_level()
#bs.write(struct+'-BAND-Result.json')
calc.write(struct+'-BAND-Result.gpw')

# Extract eigenenergies into a file for plotting with some external package

import numpy as np
calc = GPAW(struct+'-BAND-Result.gpw', txt=None)
eps_skn = np.array([[calc.get_eigenvalues(k,s)
                     for k in range(band_npoints)]
                    for s in range(1)]) - ef

f = open(struct+'-BAND-Result.dat', 'w')
for n in range(num_of_bands):
    for k in range(band_npoints):
        print(k, eps_skn[0, k, n], end="\n", file=f)
    print (end="\n", file=f)
f.close()


# - - - - - - - - - COLUMNED OUTPUT - - - - - - - - - - 
# create a  matrix of zeroes
arr = [[0 for col in range(2*num_of_bands+1)] for row in range(band_npoints+1)]
f = open(struct+'-BAND-Result.dat', 'r')
lines = f.readlines()
f.close()
a = 0 

for i in range(0, num_of_bands, 1):
   b = 0
   for a in range (a, a+band_npoints, 1):
      fields = lines[a].split()
      arr[b][2*i] = fields[0]
      arr[b][2*i+1] = fields[1]
      b = b + 1
   a = a + 2

# writing to output file
f = open(struct+'-BAND-Result-withColumns.dat', 'w')
stringline = ""

for i in range(0, band_npoints, 1):
    stringline = stringline + arr[i][0] + " " + arr[i][1] + " "
    for j in range(1, num_of_bands, 1):
        stringline = stringline + arr[i][2*j+1] + " "
    f.write(stringline + "\n")
    stringline = ""

f.close()
