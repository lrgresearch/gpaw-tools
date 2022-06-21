Outdirname = 'Si-results'

# -------------------------------------------------------------
Mode = 'PW'             # Use PW, PW-GW, PW-EXX, LCAO, FD  (PW is more accurate, LCAO is quicker mostly.)
# -------------------------------------------------------------
Geo_optim = False       # Geometric optimization with LFBGS
Elastic_calc = False    # Elastic calculation
DOS_calc = False         # DOS calculation
Band_calc = False        # Band structure calculation
Density_calc = False    # Calculate the all-electron density?
Optical_calc = True     # Calculate the optical properties

# -------------------------------------------------------------
# Parameters
# -------------------------------------------------------------
# ELECTRONIC
fmaxval = 0.05 			#
Fix_symmetry = False    # True for preserving the spacegroup symmetry during optimisation
cut_off_energy = 340 	# eV
#kpts_density = 2.5     # pts per Å^-1  If the user prefers to use this, kpts_x,y,z will not be used automatically.
kpts_x = 4 			    # kpoints in x direction
kpts_y = 4				# kpoints in y direction
kpts_z = 4				# kpoints in z direction
Gamma = True
band_path = 'GXWKL'	    # Brillouin zone high symmetry points
band_npoints = 40		# Number of points between high symmetry points
energy_max = 15 		# eV. It is the maximum energy value for band structure figure.
Hubbard = {}            # Can be used like {'N': ':p,6.0'}, for none use {}

XC_calc = 'PBE'         # Exchange-Correlation, choose one: LDA, PBE, GLLBSCM, HSE06, HSE03, revPBE, RPBE, PBE0(for PW-EXX)

Ground_convergence = {}   # Convergence items for ground state calculations
Band_convergence = {'bands':8}   # Convergence items for band calculations
Occupation = {'name': 'fermi-dirac', 'width': 0.05}  # Refer to GPAW docs: https://wiki.fysik.dtu.dk/gpaw/documentation/basic.html#occupation-numbers

DOS_npoints = 501        # Number of points
DOS_width = 0.1          # Width of Gaussian smearing. Use 0.0 for linear tetrahedron interpolation

Spin_calc = False        # Spin polarized calculation?
Magmom_per_atom = 1.0    # Magnetic moment per atom
gridref = 4             # refine grid for all electron density (1, 2 [=default] and 4)

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
MPIcores = 4            # Number of cores in calculation.
