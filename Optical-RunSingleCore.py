from ase import *
from ase.parallel import paropen, world
from gpaw import GPAW, PW, FermiDirac
from gpaw.response.df import DielectricFunction
from ase.io import read, write
from pathlib import Path

# Sample Optical Calculation GPAW Input for LRG Studies
# by Sefer Bora Lisesivdin
# July 2021 - First version
# Usage: $ gpaw python script.py
# -------------------------------------------------------------
# Parameters
# -------------------------------------------------------------

struct = 'deneme'       # Name of the electronic calculating python file without .py
num_of_bands = 16		#
optFDsmear = 0.05      # Fermi Dirac smearing for optical calculations
opteta=0.05             # Eta for Optical calculations
optdomega0=0.02         # Domega0 for Optical calculations
optecut=150             # Ecut for Optical calculations


# -------------------------------------------------------------
# Step 4 - OPTICAL CALCULATION
#      /// RUN SINGLE CORE ONLY!!! \\\
# -------------------------------------------------------------
calc = GPAW(struct+'-1-Result-Ground.gpw',
	    txt=struct+'-4-Log-Optical.txt',
	    nbands=num_of_bands,
	    fixdensity=True,
	    symmetry='off',
        occupations=FermiDirac(optFDsmear))

calc.get_potential_energy()

calc.diagonalize_full_hamiltonian(nbands=num_of_bands)  # diagonalize Hamiltonian
calc.write(struct+'-4-Result-Optical.gpw', 'all')  # write wavefunctions

# Getting absorption spectrum
df = DielectricFunction(calc=struct+'-4-Result-Optical.gpw',
                        eta=opteta,
                        domega0=optdomega0,
                        ecut=optecut)
df.get_dielectric_function( direction='x', filename=struct+'-4-Result-Optical_abs_xdirection.csv')
df.get_dielectric_function( direction='y', filename=struct+'-4-Result-Optical_abs_ydirection.csv')
df.get_dielectric_function( direction='z', filename=struct+'-4-Result-Optical_abs_zdirection.csv')