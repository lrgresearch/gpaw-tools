from ase import Atoms

Outdirname = 'bulk-aluminum-results'

bulk_configuration = Atoms('Al',
                            cell=[[0, 2.025, 2.025],
                                 [2.025, 0, 2.025],
                                 [2.025, 2.025, 0]],
                                 pbc=True)

# -------------------------------------------------------------
Mode = 'PW'             # Use PW, PW-GW, PW-EXX, LCAO, FD  (PW is more accurate, LCAO is quicker mostly.)
# -------------------------------------------------------------
Geo_optim = True       # Geometric optimization with LFBGS
Elastic_calc = False    # Elastic calculation
DOS_calc = True         # DOS calculation
Band_calc = True        # Band structure calculation
Density_calc = False    # Calculate the all-electron density?
Optical_calc = False     # Calculate the optical properties

# -------------------------------------------------------------
# Parameters
# -------------------------------------------------------------
# GEOMETRY
Minimizer = 'FIRE'      # QuasiNewton, GPMin, LBFGS or FIRE
fmaxval = 0.05 			# Maximum force tolerance in LBFGS geometry optimization. Unit is eV/Ang.
Max_step = 0.2    # How far is a single atom allowed to move. Default is 0.2 Ang.
Alpha = 70.0      # LBFGS only: Initial guess for the Hessian (curvature of energy surface)
Damping = 1.0     # LBFGS only: The calculated step is multiplied with this number before added to the positions
Fix_symmetry = False    # True for preserving the spacegroup symmetry during optimisation
# Which components of strain will be relaxed: EpsX, EpsY, EpsZ, ShearYZ, ShearXZ, ShearXY
# Example: For a x-y 2D nanosheet only first 2 component will be true
whichstrain=[True, True, True, False, False, False]

# ELECTRONIC
cut_off_energy = 340 	# eV
kpts_density = 2.5     # pts per Å^-1  If the user prefers to use this, kpts_x,y,z will not be used automatically.
kpts_x = 4 			    # kpoints in x direction
kpts_y = 4				# kpoints in y direction
kpts_z = 4				# kpoints in z direction
Gamma = True
band_path = 'LGXG'	    # Brillouin zone high symmetry points
band_npoints = 40		# Number of points between high symmetry points
energy_max = 15 		# eV. It is the maximum energy value for band structure figure.
Hubbard = {}            # Can be used like {'N': ':p,6.0'}, for none use {}

XC_calc = 'LDA'         # Exchange-Correlation, choose one: LDA, PBE, GLLBSCM, HSE06, HSE03, revPBE, RPBE, PBE0(for PW-EXX)

Ground_convergence = {}   # Convergence items for ground state calculations
Band_convergence = {'bands':8}   # Convergence items for band calculations
Occupation = {'name': 'fermi-dirac', 'width': 0.05}  # Refer to GPAW docs: https://wiki.fysik.dtu.dk/gpaw/documentation/basic.html#occupation-numbers

DOS_npoints = 501        # Number of points
DOS_width = 0.1          # Width of Gaussian smearing.  Use 0.0 for linear tetrahedron interpolation

Spin_calc = False        # Spin polarized calculation?
Magmom_per_atom = 1.0    # Magnetic moment per atom
gridref = 4             # refine grid for all electron density (1, 2 [=default] and 4)

#GENERAL
MPIcores = 4            # Number of cores in calculation.
