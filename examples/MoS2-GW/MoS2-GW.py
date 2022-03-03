import numpy as np

# -------------------------------------------------------------
Mode = 'PW-GW'             # Use PW, PW-GW, PW-EXX, LCAO, FD  (PW is more accurate, LCAO is quicker mostly.)
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
cut_off_energy = 300 	# eV
#kpts_density = 2.5     # pts per Å^-1  If the user prefers to use this, kpts_x,y,z will not be used automatically.
kpts_x = 3 			    # kpoints in x direction
kpts_y = 3				# kpoints in y direction
kpts_z = 1				# kpoints in z direction
Gamma = True
band_path = 'GKG'	    # Brillouin zone high symmetry points
band_npoints = 40		# Number of points between high symmetry points
energy_max = 10		# eV. It is the maximum energy value for band structure figure.
Hubbard = {}            # Can be used like {'N': ':p,6.0'}, for none use {}
#Exchange-Correlation, choose one:
XC_calc = 'LDA'
#XC_calc = 'PBE'
#XC_calc = 'GLLBSC'
#XC_calc = 'revPBE'
#XC_calc = 'RPBE'
#Choose one for PW-EXX (Ground state calculations will be done with PBE):
#XC_calc = 'PBE0'
#XC_calc = 'HSE06'
DOS_npoints = 501        # Number of points
DOS_width = 0.1          # Width of Gaussian smearing. Use 0.0 for linear tetrahedron interpolation

Spin_calc = False        # Spin polarized calculation?
Magmom_per_atom = 1.0    # Magnetic moment per atom
gridref = 4             # refine grid for all electron density (1, 2 [=default] and 4)

#GW Parameters
GWtype = 'G0W0'          # GW0 or G0W0
GWkpoints = np.array([[0.0, 0.0, 0.0], [1 / 3, 1 / 3, 0], [0.0, 0.0, 0.0]]) #Kpoints list
GWtruncation = '2D'     # Can be None, '2D', '1D', '0D' or 'wigner-seitz'
GWcut_off_energy = 60   # Cut-off energy
GWbandVB = 8            # Valence band number
GWbandCB = 18           # Conduction band number
GWppa = True            # Plasmon Pole Approximation
GWq0correction = True   # Analytic correction to the q=0 contribution applicable to 2D systems.
GWnblock = True         # Cuts chi0 into as many blocks to reduce mem. req. as much as possible.
GWbandinterpolation = False

# OPTICAL
num_of_bands = 16		#
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
whichstrain=[False, False, False, False, False, False]
MPIcores = 4            # Number of cores in calculation.
