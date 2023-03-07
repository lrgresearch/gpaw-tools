
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

Spin_calc = False        # Spin polarized calculation?

# OPTICAL
Opt_calc_type = 'RPA'      # BSE or RPA
Opt_shift_en = 0.0         # Shifting of the energy
Opt_BSE_valence = range(0,3)    # Valence bands that will be used in BSE calculation
Opt_BSE_conduction = range(4,7) # Conduction bands that will be used in BSE calculation
Opt_BSE_min_en = 0.0       # Results will be started from this energy (BSE only)
Opt_BSE_max_en = 20.0      # Results will be ended at this energy (BSE only)
Opt_BSE_num_of_data = 1001 # Number of data points in BSE  calculation
Opt_num_of_bands = 8	   # Number of bands
Opt_FD_smearing = 0.05     # Fermi Dirac smearing for optical calculations
Opt_eta = 0.05             # Eta for Optical calculations
Opt_domega0 = 0.05         # Domega0 for Optical calculations
Opt_omega2 = 5.0           # Frequency at which the non-lin freq grid has doubled the spacing
Opt_cut_of_energy = 100    # Cut-off energy for optical calculations
Opt_nblocks = 4            # Split matrices in nblocks blocks and distribute them G-vectors
                           # or frequencies over processes or can use world.size

#GENERAL
MPI_cores = 4            # Number of cores in calculation.
