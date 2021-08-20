from ase import *
from ase.parallel import paropen, world, parprint
from gpaw import GPAW, PW
from ase.optimize.lbfgs import LBFGS
from ase.io import read, write
import matplotlib.pyplot as plt
from ase.dft.dos import DOS
from ase.constraints import UnitCellFilter
from ase.io.cif import write_cif
from pathlib import Path
import numpy as np

# Sample Electronic Calculation GPAW Input for LRG Studies
# by Sefer Bora Lisesivdin
# August 2021 - BFGS to LBFGS, Small many changes , Strain, CIF Export, Spin polarized results, 
#               Several XC, better parallel computation, all-electron density
# July 2021 - Corrected version
# March 2020 - First Version 
# Usage: Change number with core numbers/threads to use. I am suggesting to use total number of cores(or threads) - 1
# Usage: $ gpaw -P8 python GPAWSimpleBenchmark2021.py
# For AMD CPUs or using Intel CPUs without hyperthreading: (Example CPU is intel here, 4 cores or 8 threads)
#        $ mpirun -n 4 gpaw python GPAWSimpleBenchmark2021.py
# For using all threads provided by Intel Hyperthreading technology:
#        $ mpirun --use-hwthread-cpus -n 8 gpaw python GPAWSimpleBenchmark2021.py 
# -------------------------------------------------------------
# Parameters
# -------------------------------------------------------------
fmaxval = 0.05 			#
cut_off_energy = 340 	# eV
kpts_x = 5 			    # kpoints in x direction
kpts_y = 5				# kpoints in y direction
kpts_z = 1				# kpoints in z direction
band_path = 'GMKG'	    # Brillouin zone high symmetry points
band_npoints = 40		# Number of points between high symmetry points 
energy_max = 15 		# eV. It is the maximum energy value for band structure figure.
#Exchange-Correlation, choose one:
#XC_calc = 'LDA'
XC_calc = 'PBE'
#XC_calc = 'revPBE'
#XC_calc = 'RPBE'
Spin_calc = False        # Spin polarized calculation?
Electron_density = False  # Calculate the all-electron density for 3D isosurface graphing? (Huge GPW file size)
gridref = 4             # refine grid for all electron density (1, 2 [=default] and 4)
draw_graphs = False		# Draw DOS and band structure on screen (yes for draw, small letters)

# Which components of strain will be relaxed
# EpsX, EpsY, EpsZ, ShearYZ, ShearXZ, ShearXY
# Example: For a x-y 2D nanosheet only first 2 component will be true
whichstrain=[False, False, False, False, False, False]

WantCIFexport = False
# -------------------------------------------------------------
# Bulk Configuration
# -------------------------------------------------------------
bulk_configuration = Atoms(
    [
    Atom('C', ( 0.0, 0.0, 5.0 )),
    Atom('C', ( -1.2339999999999995, 2.1373506965399947, 5.0 )),
    Atom('C', ( 2.4679999999999995, 0.0, 5.0 )),
    Atom('C', ( 1.234, 2.1373506965399947, 5.0 )),
    Atom('C', ( 2.468000000230841e-06, 1.424899039459532, 5.0 )),
    Atom('C', ( -1.2339975319999992, 3.5622497359995267, 5.0 )),
    Atom('C', ( 2.4680024680000003, 1.424899039459532, 5.0 )),
    Atom('C', ( 1.234002468000001, 3.5622497359995267, 5.0 )),
    ],
    cell=[(4.936, 0.0, 0.0), (-2.467999999999999, 4.274701393079989, 0.0), (0.0, 0.0, 20.0)],
    pbc=True,
    )

# -------------------------------------------------------------
# ///////   YOU DO NOT NEED TO CHANGE ANYTHING BELOW    \\\\\\\ 
# -------------------------------------------------------------
struct = Path(__file__).stem # All files will get their names from this file
# -------------------------------------------------------------
# Step 1 - GROUND STATE
# -------------------------------------------------------------
calc = GPAW(mode=PW(cut_off_energy), xc=XC_calc, parallel={'domain': world.size}, spinpol=Spin_calc, kpts=[kpts_x, kpts_y, kpts_z], txt=struct+'-1-Log-Ground.txt')
bulk_configuration.calc = calc

uf = UnitCellFilter(bulk_configuration, mask=whichstrain)
relax = LBFGS(uf, trajectory=struct+'-1-Result-Ground.traj')
relax.run(fmax=fmaxval)  # Consider tighter fmax!

bulk_configuration.get_potential_energy()
if Electron_density == True:
    #This line makes huge GPW files. Therefore it is better to use this if else
    calc.write(struct+'-1-Result-Ground.gpw', mode="all")
else:
    calc.write(struct+'-1-Result-Ground.gpw')

if WantCIFexport == True:
    write_cif(struct+'-Final.cif', bulk_configuration)
    
# -------------------------------------------------------------
# Step 2 - DOS CALCULATION
# -------------------------------------------------------------
calc = GPAW(struct+'-1-Result-Ground.gpw', txt=struct+'-2-Log-DOS.txt')
#energies, weights = calc.get_dos(npts=800, width=0)
dos = DOS(calc, npts=500, width=0)
if Spin_calc == True:
    energies = dos.get_energies()
    weights = dos.get_dos(spin=0)
    weightsup = dos.get_dos(spin=1)
else:
    energies = dos.get_energies()
    weights = dos.get_dos()

fd = open(struct+'-2-Result-DOS.txt', "w")
if Spin_calc == True:
    for x in zip(energies, weights, weightsup):
        print(*x, sep=", ", file=fd)
else:
    for x in zip(energies, weights):
        print(*x, sep=", ", file=fd)
fd.close()

# -------------------------------------------------------------
# Step 3 - BAND STRUCTURE CALCULATION
# -------------------------------------------------------------
calc = GPAW(struct+'-1-Result-Ground.gpw',
	    txt=struct+'-3-Log-Band.txt',
	    fixdensity=True,
	    symmetry='off',
	    kpts={'path': band_path, 'npoints': band_npoints},
	    convergence={'bands': 8})

calc.get_potential_energy()
bs = calc.band_structure()
ef = calc.get_fermi_level()
num_of_bands = calc.get_number_of_bands()
parprint('Num of bands:'+str(num_of_bands))

#bs.write(struct+'-3-Result-Band.json')
calc.write(struct+'-3-Result-Band.gpw')

if Spin_calc == True:
    eps_skn = np.array([[calc.get_eigenvalues(k,s)
                         for k in range(band_npoints)]
                        for s in range(2)]) - ef
    parprint(eps_skn.shape)
    f1 = open(struct+'-3-Result-Band-Down.dat', 'w')
    for n1 in range(num_of_bands):
        for k1 in range(band_npoints):
            print(k1, eps_skn[0, k1, n1], end="\n", file=f1)
        print (end="\n", file=f1)
    f1.close()
    f2 = open(struct+'-3-Result-Band-Up.dat', 'w')
    for n2 in range(num_of_bands):
        for k2 in range(band_npoints):
            print(k2, eps_skn[1, k2, n2], end="\n", file=f2)
        print (end="\n", file=f2)
    f2.close()
else:
    eps_skn = np.array([[calc.get_eigenvalues(k,s)
                         for k in range(band_npoints)]
                        for s in range(1)]) - ef
    f = open(struct+'-3-Result-Band.dat', 'w')
    for n in range(num_of_bands):
        for k in range(band_npoints):
            print(k, eps_skn[0, k, n], end="\n", file=f)
        print (end="\n", file=f)
    f.close()

# -------------------------------------------------------------
# Step 4 - ALL-ELECTRON DENSITY
# -------------------------------------------------------------
calc = GPAW(struct+'-1-Result-Ground.gpw', txt=struct+'-4-Log-ElectronDensity.txt')
bulk_configuration.calc = calc
np = calc.get_pseudo_density()
n = calc.get_all_electron_density(gridrefinement=gridref)

write(struct+'-4-Result-All-electron_n.cube', bulk_configuration, data=n)
write(struct+'-4-Result-All-electron_np.cube', bulk_configuration, data=np)

# -------------------------------------------------------------
# Step 5 - DRAWING BAND STRUCTURE AND DOS
# -------------------------------------------------------------
if draw_graphs == True:
    # Draw graphs only on master node
    if world.rank == 0:
        # DOS
        if Spin_calc == True:
            ax = plt.gca()
            ax.plot(energies, -1.0*weights, 'r')
            ax.plot(energies, weightsup, 'b')
            ax.set_xlabel('Energy [eV]')
            ax.set_ylabel('DOS [1/eV]')
        else:
            ax = plt.gca()
            ax.plot(energies, weights, 'b')
            ax.set_xlabel('Energy [eV]')
            ax.set_ylabel('DOS [1/eV]')
        plt.savefig(struct+'-2-Graph-DOS.png')
        #plt.show()
        # Band Structure
        bs.plot(filename=struct+'-3-Graph-Band.png', show=True, emax=energy_max)
