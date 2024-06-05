# -------------------------------------------------------------
Mode = 'PW'             # Use PW, PW-GW, LCAO, FD  (PW is more accurate, LCAO is quicker mostly.)
# -------------------------------------------------------------
Ground_calc = True     # Ground state calculations
Geo_optim = False       # Geometric optimization with LFBGS
Elastic_calc = False    # Elastic calculation
DOS_calc = True         # DOS calculation
Band_calc = True        # Band structure calculation
Density_calc = False    # Calculate the all-electron density?
Optical_calc = False     # Calculate the optical properties

# -------------------------------------------------------------
# Parameters
# -------------------------------------------------------------
# GEOMETRY
Optimizer = 'LBFGS'     # QuasiNewton, GPMin, LBFGS or FIRE
Max_F_tolerance = 0.05 			# Maximum force tolerance in LBFGS geometry optimization. Unit is eV/Ang.
Max_step = 0.1          # How far is a single atom allowed to move. Default is 0.2 Ang.
Alpha = 60.0            # LBFGS only: Initial guess for the Hessian (curvature of energy surface)
Damping = 1.0           # LBFGS only: The calculated step is multiplied with this number before added to the positions
Fix_symmetry = True    # True for preserving the spacegroup symmetry during optimisation
# Which components of strain will be relaxed: EpsX, EpsY, EpsZ, ShearYZ, ShearXZ, ShearXY
# Example: For a x-y 2D nanosheet only first 2 component will be true
Relax_cell = [False, False, False, False, False, False]
Hydrostatic_pressure=0.0 #GPa

# ELECTRONIC
Cut_off_energy = 400 	# eV
#Ground_kpts_density = 3.0     # pts per Å^-1  If the user prefers to use this, kpts_x,y,z will not be used automatically.
Ground_kpts_x = 3 			    # kpoints in x direction
Ground_kpts_y = 3				# kpoints in y direction
Ground_kpts_z = 3				# kpoints in z direction
Gamma = True
Band_path = 'GXWKGLUWLK'	    # Brillouin zone high symmetry points
Band_npoints = 401		# Number of points between high symmetry points
Setup_params = {}            # Can be used like {'N': ':p,6.0'}, for none use {}

XC_calc = 'HSE06'       # Exchange-Correlation, choose one: LDA, PBE, GLLBSCM, HSE06, HSE03, revPBE, RPBE, PBE0, EXX, B3LYP

# These convergence values listed below kept low for this example to finish the calculation quicker.
# Otherwise HSE calculations are much more slower than standard PBE calculations (Sometimes few thousand times slower).
# Please use proper convergence values and always use HPC for your HSE calculations :)

Ground_convergence = {'energy':1e-1, 'eigenstates':1e-1, 'density':1e-1}   # Convergence items for ground state calculations
Band_convergence = {'bands':8, 'eigenstates':1e-1, 'density':1e-1}   # Convergence items for band calculations
Occupation = {'name': 'marzari-vanderbilt', 'width': 0.2}  # Refer to GPAW docs: https://wiki.fysik.dtu.dk/gpaw/documentation/basic.html#occupation-numbers

DOS_npoints = 301        # Number of points
DOS_width = 0.2          # Width of Gaussian smearing.  Use 0.0 for linear tetrahedron interpolation
DOS_convergence = {}  # Convergence items for DOS calculations

Spin_calc = False        # Spin polarized calculation?
Magmom_per_atom = 1.0    # Magnetic moment per atom
Refine_grid = 4             # refine grid for all electron density (1, 2 [=default] and 4)
Total_charge = 0.0       # Total charge. Normally 0.0 for a neutral system.

#GENERAL
MPI_cores = 4            # Number of cores in calculation.
Energy_min = -5 		# eV. It is the minimum energy value for band structure and DOS figures.
Energy_max = 5  		# eV. It is the maximum energy value for band structure and DOS figures.
Localisation = "en_UK"  # Localisation setting for figures. en_UK is default.
