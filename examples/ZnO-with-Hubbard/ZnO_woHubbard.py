from ase.build import bulk
import numpy as np

Outdirname = 'ZnO-results'

bulk_configuration = bulk('ZnO', 'wurtzite', a=3.25, c=5.2)

# -------------------------------------------------------------
Mode = 'PW'             # Use PW, PW-GW, PW-EXX, LCAO, FD  (PW is more accurate, LCAO is quicker mostly.)
# -------------------------------------------------------------
Elastic_calc = False    # Elastic calculation
DOS_calc = True         # DOS calculation
Band_calc = True        # Band structure calculation
Density_calc = False    # Calculate the all-electron density?
Optical_calc = False     # Calculate the optical properties

# -------------------------------------------------------------
# Parameters
# -------------------------------------------------------------
# ELECTRONIC
fmaxval = 0.05 			#
Fix_symmetry = False    # True for preserving the spacegroup symmetry during optimisation
cut_off_energy = 340 	# eV
#kpts_density = 2.5     # pts per Å^-1  If the user prefers to use this, kpts_x,y,z will not be used automatically.
kpts_x = 5 			    # kpoints in x direction
kpts_y = 5				# kpoints in y direction
kpts_z = 5				# kpoints in z direction
Gamma = True
band_path = 'ALMGAHKG'	    # Brillouin zone high symmetry points
band_npoints = 40		# Number of points between high symmetry points
energy_max = 15 		# eV. It is the maximum energy value for band structure figure.
Hubbard = {}  # Can be used like {'N': ':p,6.0,0'}, for none use {}
#Exchange-Correlation, choose one:
#XC_calc = 'LDA'
XC_calc = 'PBE'
#XC_calc = 'GLLBSC'
#XC_calc = 'revPBE'
#XC_calc = 'RPBE'
#Choose one for PW-EXX (Ground state calculations will be done with PBE):
#XC_calc = 'PBE0'
#XC_calc = 'HSE06'

Ground_convergence = {}   # Convergence items for ground state calculations
Band_convergence = {'bands':8}   # Convergence items for band calculations
Occupation = {'name': 'fermi-dirac', 'width': 0.05}  # Refer to GPAW docs: https://wiki.fysik.dtu.dk/gpaw/documentation/basic.html#occupation-numbers

DOS_npoints = 501        # Number of points
DOS_width = 0.1          # Width of Gaussian smearing. Use 0.0 for linear tetrahedron interpolation

Spin_calc = False        # Spin polarized calculation?
Magmom_per_atom = 1.0    # Magnetic moment per atom
gridref = 4             # refine grid for all electron density (1, 2 [=default] and 4)

#GW Parameters
GWtype = 'GW0'          # GW0 or G0W0
GWkpoints = np.array([[0.0, 0.0, 0.0], [1 / 3, 1 / 3, 0], [0.0, 0.0, 0.0]]) #Kpoints list
GWtruncation = '2D'     # Can be None, '2D', '1D', '0D' or 'wigner-seitz'
GWcut_off_energy = 50   # Cut-off energy
GWbandVB = 8            # Valence band number
GWbandCB = 18           # Conduction band number
GWppa = True            # Plasmon Pole Approximation
GWq0correction = True   # Analytic correction to the q=0 contribution applicable to 2D systems.
GWnblock = True         # Cuts chi0 into as many blocks to reduce mem. req. as much as possible.

# OPTICAL
opttype = 'BSE'         # BSE or RPA
optshift = 0.0          # Shifting of the energy
optBSEvb = range(0,3)  # Valence bands that will be used in BSE calculation
optBSEcb = range(4,7) # Conduction bands that will be used in BSE calculation
optBSEminEn = 0.0       # Results will be started from this energy (BSE only)
optBSEmaxEn = 20.0      # Results will be ended at this energy (BSE only)
optBSEnumdata = 1001   # Number of data points in BSE  calculation
num_of_bands = 8	# Number of bands
optFDsmear = 0.05       # Fermi Dirac smearing for optical calculations
opteta=0.05             # Eta for Optical calculations
optdomega0=0.05         # Domega0 for Optical calculations
optomega2=5.0           # Frequency at which the non-lin freq grid has doubled the spacing
optecut=100             # Cut-off energy for optical calculations
optnblocks=4            # Split matrices in nblocks blocks and distribute them G-vectors
                        # or frequencies over processes

#GENERAL
# Which components of strain will be relaxed
# EpsX, EpsY, EpsZ, ShearYZ, ShearXZ, ShearXY
# Example: For a x-y 2D nanosheet only first 2 component will be true
whichstrain=[True, True, True, False, False, False]
MPIcores = 4            # Number of cores in calculation.
