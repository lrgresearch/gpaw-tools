import numpy as np

# -------------------------------------------------------------
Mode = 'PW-GW'             # Use PW, PW-GW, PW-EXX, LCAO, FD  (PW is more accurate, LCAO is quicker mostly.)
# -------------------------------------------------------------
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
Max_F_tolerance = 0.05 	# Maximum force tolerance in LBFGS geometry optimization. Unit is eV/Ang.
Max_step = 0.1          # How far is a single atom allowed to move. Default is 0.2 Ang.
Alpha = 60.0            # LBFGS only: Initial guess for the Hessian (curvature of energy surface)
Damping = 1.0           # LBFGS only: The calculated step is multiplied with this number before added to the positions
Fix_symmetry = True    # True for preserving the spacegroup symmetry during optimisation
# Which components of strain will be relaxed: EpsX, EpsY, EpsZ, ShearYZ, ShearXZ, ShearXY
# Example: For a x-y 2D nanosheet only first 2 component will be true
Relax_cell = [False, False, False, False, False, False]

# ELECTRONIC
Cut_off_energy = 300 	# eV
#Ground_kpts_density = 2.5     # pts per Å^-1  If the user prefers to use this, kpts_x,y,z will not be used automatically.
Ground_kpts_x = 3 			    # kpoints in x direction
Ground_kpts_y = 3				# kpoints in y direction
Ground_kpts_z = 1				# kpoints in z direction
Gamma = True
Band_path = 'GKG'	    # Brillouin zone high symmetry points
Band_npoints = 40		# Number of points between high symmetry points
Setup_params = {}            # Can be used like {'N': ':p,6.0'}, for none use {}

XC_calc = 'LDA'         # Exchange-Correlation, choose one: LDA, PBE, GLLBSCM, HSE06, HSE03, revPBE, RPBE, PBE0(for PW-EXX)


Ground_convergence = {}   # Convergence items for ground state calculations
Band_convergence = {'bands':8}   # Convergence items for band calculations
Occupation = {'name': 'fermi-dirac', 'width': 0.001}  # Refer to GPAW docs: https://wiki.fysik.dtu.dk/gpaw/documentation/basic.html#occupation-numbers

DOS_npoints = 501        # Number of points
DOS_width = 0.1          # Width of Gaussian smearing. Use 0.0 for linear tetrahedron interpolation
DOS_convergence = {}  # Convergence items for DOS calculations

Spin_calc = False        # Spin polarized calculation?
Magmom_per_atom = 1.0    # Magnetic moment per atom
Refine_grid = 4             # refine grid for all electron density (1, 2 [=default] and 4)

#GW Parameters
GW_calc_type = 'G0W0'          # GW0 or G0W0
GW_kpoints_list = np.array([[0.0, 0.0, 0.0], [1 / 3, 1 / 3, 0], [0.0, 0.0, 0.0]]) #Kpoints list
GW_truncation = '2D'     # Can be None, '2D', '1D', '0D' or 'wigner-seitz'
GW_cut_off_energy = 60   # Cut-off energy
GW_valence_band_no = 8            # Valence band number
GW_conduction_band_no = 18           # Conduction band number
GW_PPA = True            # Plasmon Pole Approximation
GW_q0_correction = True   # Analytic correction to the q=0 contribution applicable to 2D systems.
GW_nblocks_max = True         # Cuts chi0 into as many blocks to reduce mem. req. as much as possible.
GW_interpolate_band = False

#GENERAL
MPI_cores = 4            # Number of cores in calculation.
Energy_min = -5 		# eV. It is the minimum energy value for band structure and DOS figures.
Energy_max = 5  		# eV. It is the maximum energy value for band structure and DOS figures.
