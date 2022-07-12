#!/usr/bin/env python

'''
gpawsolve.py: High-level Interaction Script for GPAW
More information: $ gpawsolve.py -h
'''

Description = f''' 
 Usage: 
 $ mpirun -np <corenumbers> gpawsolve.py <args>
 -------------------------------------------------------------
 Calculation selector
 -------------------------------------------------------------
 | Method | XCs             | Structure optim. | Spin | Ground | Elastic | DOS | DFT+U | Band | Density | Optical |
 | ------ | --------------- | ---------------- | ---- | ------ | ------- | --- | ----- | ---- | ------- | ------- |
 |   PW   | Local and LibXC | Yes              | Yes  | Yes    | Yes     | Yes | Yes   | Yes  | Yes     | Yes     |
 |   PW   | GLLBSC / M      | No               | Yes  | Yes    | Yes     | Yes | No    | Yes  | Yes     | Yes     |
 |   PW   | HSE03, HSE06    | No               | Yes  | Yes    | n/a     | Yes | No    | No   | No      | No      |
 | PW-G0W0| Local and LibXC | No               | No   | Yes    | No      | No  | No    | Some | No      | No      |
 | PW-EXX*| B3LYP, PBE0     | Yes (with PBE)   | No   | Yes    | No      | No  | No    | No   | No      | No      |
 |  LCAO  | Local and LibXC | Yes              | Yes  | Yes    | Yes     | Yes | Yes   | Yes  | Yes     | No      |
 *: Just some ground state energy calculations.
'''

import getopt, sys, os, time
import textwrap
import requests
import pickle
import spglib as spg
from argparse import ArgumentParser, HelpFormatter
from ase import *
from ase.dft.kpoints import get_special_points
from ase.parallel import paropen, world, parprint, broadcast
from gpaw import GPAW, PW, Davidson, FermiDirac, MixerSum, MixerDif, Mixer
from ase.optimize import QuasiNewton
from ase.io import read, write
from ase.eos import calculate_eos
from ase.units import Bohr, GPa, kJ
import matplotlib.pyplot as plt
from ase.dft.dos import DOS
from ase.constraints import ExpCellFilter
from ase.spacegroup.symmetrize import FixSymmetry
from ase.io.cif import write_cif
from pathlib import Path
from gpaw.response.df import DielectricFunction
from gpaw.response.bse import BSE
from gpaw.response.g0w0 import G0W0
from gpaw.response.gw_bands import GWBands
from gpaw.xc.exx import EXX
from gpaw.dos import DOSCalculator
import numpy as np
from numpy import genfromtxt
from elastic import get_elastic_tensor, get_elementary_deformations

# DEFAULT VALUES
# These values (with bulk configuration) can be used to run this script without using inputfile (py file)
# and configuration file (cif file). 
# -------------------------------------------------------------
Mode = 'PW'             # Use PW, PW-GW, PW-EXX, LCAO, FD  (PW is more accurate, LCAO is quicker mostly.)
# -------------------------------------------------------------
Geo_optim = True       # Geometric optimization with LFBGS
Elastic_calc = False    # Elastic calculation
DOS_calc = False         # DOS calculation
Band_calc = False        # Band structure calculation
Density_calc = False    # Calculate the all-electron density?
Optical_calc = False     # Calculate the optical properties

# -------------------------------------------------------------
# Parameters
# -------------------------------------------------------------
# GEOMETRY
Optimizer = 'QuasiNewton' # QuasiNewton, GPMin, LBFGS or FIRE
fmaxval = 0.05 			# Maximum force tolerance in LBFGS geometry optimization. Unit is eV/Ang.
Max_step = 0.1          # How far is a single atom allowed to move. Default is 0.2 Ang.
Alpha = 60.0            # LBFGS only: Initial guess for the Hessian (curvature of energy surface)
Damping = 1.0           # LBFGS only: The calculated step is multiplied with this number before added to the positions
Fix_symmetry = False    # True for preserving the spacegroup symmetry during optimisation
# Which components of strain will be relaxed: EpsX, EpsY, EpsZ, ShearYZ, ShearXZ, ShearXY
# Example: For a x-y 2D nanosheet only first 2 component will be true
whichstrain=[False, False, False, False, False, False]

# ELECTRONIC
cut_off_energy = 340 	# eV
#kpts_density = 2.5     # pts per Å^-1  If the user prefers to use this, kpts_x,y,z will not be used automatically.
kpts_x = 5 			    # kpoints in x direction
kpts_y = 5				# kpoints in y direction
kpts_z = 5				# kpoints in z direction
gpts_density = 0.2      # (for LCAO) Unit is Å. If the user prefers to use this, gpts_x,y,z will not be used automatically.
gpts_x = 8              # grid points in x direction (for LCAO)
gpts_y = 8              # grid points in y direction (for LCAO)
gpts_z = 8              # grid points in z direction (for LCAO)

Gamma = True
band_path = 'LGL'	    # Brillouin zone high symmetry points
band_npoints = 60		# Number of points between high symmetry points
energy_max = 15 		# eV. It is the maximum energy value for band structure figure.
Hubbard = {}            # Can be used like {'N': ':p,6.0'}, for none use {}

XC_calc = 'LDA'         # Exchange-Correlation, choose one: LDA, PBE, GLLBSCM, HSE06, HSE03, revPBE, RPBE, PBE0(for PW-EXX)

Ground_convergence = {}   # Convergence items for ground state calculations
Band_convergence = {'bands':8}   # Convergence items for band calculations
Occupation = {'name': 'fermi-dirac', 'width': 0.05}  # Refer to GPAW docs: https://wiki.fysik.dtu.dk/gpaw/documentation/basic.html#occupation-numbers
Mixer_type = MixerSum(0.1, 3, 50) # MixerSum(beta,nmaxold, weight) default:(0.1,3,50), you can try (0.02, 5, 100) and (0.05, 5, 50)

DOS_npoints = 501        # Number of points
DOS_width = 0.1          # Width of Gaussian smearing. Use 0.0 for linear tetrahedron interpolation

Spin_calc = False        # Spin polarized calculation?
Magmom_per_atom = 1.0    # Magnetic moment per atom
gridref = 4             # refine grid for all electron density (1, 2 [=default] and 4)

#GW Parameters
GWtype = 'GW0'          # GW0 or G0W0
GWkpoints = np.array([[0.0, 0.0, 0.0], [1 / 3, 1 / 3, 0], [0.0, 0.0, 0.0]]) #Kpoints list
GWtruncation = 'None'     # Can be None, '2D', '1D', '0D' or 'wigner-seitz'
GWcut_off_energy = 50   # Cut-off energy
GWbandVB = 8            # Valence band number
GWbandCB = 18           # Conduction band number
GWppa = True            # Plasmon Pole Approximation
GWq0correction = True   # Analytic correction to the q=0 contribution applicable to 2D systems.
GWnblock = True         # Cuts chi0 into as many blocks to reduce mem. req. as much as possible.
GWbandinterpolation = True # Interpolate band

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
MPIcores = 4            # This is for gg.py. Not used in this script.

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
# Version
__version__ = "v22.7.0"

# Start time
time0 = time.strftime('%Y-%m-%d %H:%M:%S', time.gmtime(time.time()))
# To print Description variable with argparse
class RawFormatter(HelpFormatter):
    def _fill_text(self, text, width, indent):
        return "\n".join([textwrap.fill(line, width) for line in textwrap.indent(textwrap.dedent(text), indent).splitlines()])

# Arguments parsing
parser = ArgumentParser(prog ='gpawtools.py', description=Description, formatter_class=RawFormatter)


parser.add_argument("-o", "--outdir", dest = "outdir", action='store_true', 
                    help="""Save everything to a output directory with naming /inputfile. If there is no input file given and 
                    Atoms object is used in gpawsolve.py file then the directory name will be /gpawsolve. 
                    If you change gpawsolve.py name to anyname.py then the directory name will be /anyname.""")
parser.add_argument("-i", "--input", dest = "inputfile", help="Use input file for calculation variables (also you can insert geometry)")
parser.add_argument("-g", "--geometry",dest ="geometryfile", help="Use CIF file for geometry")
parser.add_argument("-v", "--version", dest="version", action='store_true')
parser.add_argument("-r", "--restart", dest="restart", action='store_true')
parser.add_argument("-d", "--drawfigures", dest="drawfigs", action='store_true', help="Draws DOS and band structure figures at the end of calculation.")

args = None

try:
    if world.rank == 0:
        args = parser.parse_args()
finally:
    args = broadcast(args, root=0, comm=world)

if args is None:
    exit(0)

outdir = False
restart = False
inFile = None
drawfigs = False
configpath = None
Outdirname = ''

try:
    if args.inputfile is not None:
        configpath = os.path.join(os.getcwd(),args.inputfile)
        sys.path.append(os.getcwd())
        # Works like from FILE import *
        conf = __import__(Path(configpath).stem, globals(), locals(), ['*'])
        for k in dir(conf):
            locals()[k] = getattr(conf, k)

    if args.geometryfile :
        inFile = os.path.join(os.getcwd(),args.geometryfile)

    if args.outdir == True:
        outdir = True

    if args.drawfigs == True:
        drawfigs = True
    
    if args.version == True:
        import gpaw
        import ase
        try:
            response = requests.get("https://api.github.com/repos/lrgresearch/gpaw-tools/releases/latest", timeout=5)
            parprint('-------------------------------------------------------------------------------------------------------')
            parprint('\033[95mgpaw-tools:\033[0m This is '+str(__version__)+' uses GPAW '+gpaw.__version__+', and ASE '+ase.__version__)
            parprint('-------------------------------------------------------------------------------------------------------')
            parprint('The latest STABLE release was '+response.json()["tag_name"]+', which is published at '+response.json()["published_at"])
            parprint('Download the latest STABLE tarball release at: '+response.json()["tarball_url"])
            parprint('Download the latest STABLE zipball release at: '+response.json()["zipball_url"])
            parprint('Download the latest DEV zipball release at: https://github.com/lrgresearch/gpaw-tools/archive/refs/heads/main.zip')
        except (requests.ConnectionError, requests.Timeout) as exception:
            parprint('-------------------------------------------------------------------------------------------------------')
            parprint('gpaw-tools: This is '+str(__version__)+' uses GPAW '+gpaw.__version__+', ASE '+ase.__version__)
            parprint('-------------------------------------------------------------------------------------------------------')
            parprint('No internet connection available.')
        quit()
        
    if args.restart == True:
        restart = True        

except getopt.error as err:
    # output error, and return with an error code
    parprint (str(err))

# If there is a CIF input, use it. Otherwise use the bulk configuration provided above.
if inFile is None:
    if Outdirname !='':
        struct = Outdirname
    else:
        struct = 'results' # All files will get their names from this file
    parprint("Number of atoms provided in Atoms object:"+str(bulk_configuration.get_global_number_of_atoms()))
else:
    struct = Path(inFile).stem
    bulk_configuration = read(inFile, index='-1')
    parprint("Number of atoms imported from CIF file:"+str(bulk_configuration.get_global_number_of_atoms()))
    parprint("Spacegroup of CIF file (SPGlib):",spg.get_spacegroup(bulk_configuration))
    parprint("Special Points usable for this spacegroup:",get_special_points(bulk_configuration.get_cell()))

# Control if outdir is set or not
if outdir is False:
    #No change is necessary
    parprint("Output directory is the main directory")
else:
    if Outdirname != '':
        structpath = os.path.join(os.getcwd(),Outdirname)
    else:
        structpath = os.path.join(os.getcwd(),struct)

    if not os.path.isdir(structpath):
        os.makedirs(structpath, exist_ok=True)
    struct = os.path.join(structpath,struct)

if Optical_calc == False:
    # -------------------------------------------------------------
    # Step 1 - GROUND STATE
    # -------------------------------------------------------------
    if Mode == 'PW':
        if XC_calc in ['B3LYP', 'PBE0']:
            parprint('\033[91mERROR:\033[0m'+XC_calc+' can be used only in PW-EXX mode...')
            quit()
        if Spin_calc == True:
           numm = [Magmom_per_atom]*bulk_configuration.get_global_number_of_atoms()
           bulk_configuration.set_initial_magnetic_moments(numm)
        
        if restart == False:
            # Start ground state
            time11 = time.time()
            # PW Ground State Calculations
            parprint("Starting PW ground state calculation...")
            if True in whichstrain:
                if XC_calc in ['GLLBSC', 'GLLBSCM', 'HSE06', 'HSE03']:
                    parprint("\033[91mERROR:\033[0m Structure optimization LBFGS can not be used with "+XC_calc+" xc.")
                    parprint("Do manual structure optimization, or do with PBE, then use its final CIF as input.")
                    parprint("Quiting...")
                    quit()
            if XC_calc in ['HSE06', 'HSE03']:
                parprint('Starting Hybrid XC calculations...')
                if 'kpts_density' in globals():
                    calc = GPAW(mode=PW(cut_off_energy), xc={'name': XC_calc, 'backend': 'pw'}, nbands='200%', parallel={'band': 1, 'kpt': 1},
                            eigensolver=Davidson(niter=1), mixer=Mixer_type,
                            spinpol=Spin_calc, kpts={'density': kpts_density, 'gamma': Gamma}, txt=struct+'-1-Log-Ground.txt',
                            convergence = Ground_convergence, occupations = Occupation)
                else:
                    calc = GPAW(mode=PW(cut_off_energy), xc={'name': XC_calc, 'backend': 'pw'}, nbands='200%', parallel={'band': 1, 'kpt': 1},
                            eigensolver=Davidson(niter=1), mixer=Mixer_type,
                            spinpol=Spin_calc, kpts={'size': (kpts_x, kpts_y, kpts_z), 'gamma': Gamma}, txt=struct+'-1-Log-Ground.txt',
                            convergence = Ground_convergence, occupations = Occupation)
            else:
                parprint('Starting calculations with '+XC_calc+'...')
                # Fix the spacegroup in the geometric optimization if wanted
                if Fix_symmetry == True:
                    bulk_configuration.set_constraint(FixSymmetry(bulk_configuration))
                if 'kpts_density' in globals():
                    calc = GPAW(mode=PW(cut_off_energy), xc=XC_calc, nbands='200%', setups= Hubbard, parallel={'domain': world.size}, 
                            spinpol=Spin_calc, kpts={'density': kpts_density, 'gamma': Gamma}, 
                            mixer=Mixer_type, txt=struct+'-1-Log-Ground.txt',
                            convergence = Ground_convergence, occupations = Occupation)
                else:
                    calc = GPAW(mode=PW(cut_off_energy), xc=XC_calc, nbands='200%', setups= Hubbard, parallel={'domain': world.size}, 
                            spinpol=Spin_calc, kpts={'size': (kpts_x, kpts_y, kpts_z), 'gamma': Gamma}, 
                            mixer=Mixer_type, txt=struct+'-1-Log-Ground.txt',
                            convergence = Ground_convergence, occupations = Occupation)
            bulk_configuration.calc = calc
            if Geo_optim == True:
                if True in whichstrain:
                    uf = ExpCellFilter(bulk_configuration, mask=whichstrain)
                    # Optimizer Selection
                    if Optimizer == 'FIRE':
                        from ase.optimize.fire import FIRE
                        relax = FIRE(uf, maxstep=Max_step, trajectory=struct+'-1-Result-Ground.traj')
                    elif  Optimizer == 'LBFGS':
                        from ase.optimize.lbfgs import LBFGS
                        relax = LBFGS(uf, maxstep=Max_step, alpha=Alpha, damping=Damping, trajectory=struct+'-1-Result-Ground.traj')
                    elif  Optimizer == 'GPMin':
                        from ase.optimize import GPMin
                        relax = GPMin(uf, trajectory=struct+'-1-Result-Ground.traj')
                    else:
                        relax = QuasiNewton(uf, maxstep=Max_step, trajectory=struct+'-1-Result-Ground.traj')       
                else:
                    # Optimizer Selection
                    if Optimizer == 'FIRE':
                        from ase.optimize.fire import FIRE
                        relax = FIRE(bulk_configuration, maxstep=Max_step, trajectory=struct+'-1-Result-Ground.traj')
                    elif  Optimizer == 'LBFGS':
                        from ase.optimize.lbfgs import LBFGS
                        relax = LBFGS(bulk_configuration, maxstep=Max_step, alpha=Alpha, damping=Damping, trajectory=struct+'-1-Result-Ground.traj')
                    elif  Optimizer == 'GPMin':
                        from ase.optimize import GPMin
                        relax = GPMin(bulk_configuration, trajectory=struct+'-1-Result-Ground.traj')
                    else:
                        relax = QuasiNewton(bulk_configuration, maxstep=Max_step, trajectory=struct+'-1-Result-Ground.traj')
                relax.run(fmax=fmaxval)  # Consider tighter fmax!
            else:
                bulk_configuration.set_calculator(calc)
                bulk_configuration.get_potential_energy()
            if Density_calc == True:
                #This line makes huge GPW files. Therefore it is better to use this if else
                calc.write(struct+'-1-Result-Ground.gpw', mode="all")
            else:
                calc.write(struct+'-1-Result-Ground.gpw')
            # Writes final configuration as CIF file
            write_cif(struct+'-Final.cif', bulk_configuration)
        else:
            parprint("Passing PW ground state calculation...")
            # Control the ground state GPW file
            if not os.path.exists(struct+'-1-Result-Ground.gpw'):
                parprint('\033[91mERROR:\033[0m'+struct+'-1-Result-Ground.gpw file can not be found. It is needed in restart mode. Quiting.')
                quit()

    elif Mode == 'PW-EXX':
        if restart == False:
            # PW Ground State Calculations
            parprint("Starting PW ground state calculation with PBE...")
            # Fix the spacegroup in the geometric optimization if wanted
            if Fix_symmetry == True:
                    bulk_configuration.set_constraint(FixSymmetry(bulk_configuration))
            if 'kpts_density' in globals():
                calc = GPAW(mode=PW(cut_off_energy), xc='PBE', parallel={'domain': world.size}, kpts={'density': kpts_density, 'gamma': Gamma},
                        convergence = Ground_convergence, mixer=Mixer_type, occupations = Occupation, txt=struct+'-1-Log-Ground.txt')
            else:
                calc = GPAW(mode=PW(cut_off_energy), xc='PBE', parallel={'domain': world.size}, kpts={'size': (kpts_x, kpts_y, kpts_z), 'gamma': Gamma},
                        convergence = Ground_convergence, mixer=Mixer_type, occupations = Occupation, txt=struct+'-1-Log-Ground.txt')
            bulk_configuration.calc = calc
            uf = ExpCellFilter(bulk_configuration, mask=whichstrain)
            # Optimizer Selection
            if Optimizer == 'FIRE':
                from ase.optimize.fire import FIRE
                relax = FIRE(uf, maxstep=Max_step, trajectory=struct+'-1-Result-Ground.traj')
            elif  Optimizer == 'LBFGS':
                from ase.optimize.lbfgs import LBFGS
                relax = LBFGS(uf, maxstep=Max_step, alpha=Alpha, damping=Damping, trajectory=struct+'-1-Result-Ground.traj')
            elif  Optimizer == 'GPMin':
                from ase.optimize import GPMin
                relax = GPMin(uf, trajectory=struct+'-1-Result-Ground.traj')
            else:
                relax = QuasiNewton(uf, maxstep=Max_step, trajectory=struct+'-1-Result-Ground.traj')       
            relax.run(fmax=fmaxval)  # Consider tighter fmax!
            calc.write(struct+'-1-Result-Ground.gpw', mode="all")
            # Writes final configuration as CIF file
            write_cif(struct+'-Final.cif', bulk_configuration)
            # Print final spacegroup information
            parprint("Final Spacegroup (SPGlib):",spg.get_spacegroup(bulk_configuration))
        else:
            parprint("Passing PW ground state calculation...")
            # Control the ground state GPW file
            if not os.path.exists(struct+'-1-Result-Ground.gpw'):
                parprint('\033[91mERROR:\033[0m'+struct+'-1-Result-Ground.gpw file can not be found. It is needed in restart mode. Quiting.')
                quit()

        if XC_calc in ['B3LYP', 'PBE0']:
            parprint('Starting PW EXX ground state calculation with '+XC_calc+' ...')
            calc_exx = EXX(struct+'-1-Result-Ground.gpw', xc=XC_calc, txt=struct+'-1-Log-EXX_mode.txt')
            bulk_configuration.calc_exx = calc_exx
            with paropen(struct+'-1-Result-Ground-EXX_mode.txt', "w") as fd:
                print('Eigenvalue contributions: ',calc_exx.get_eigenvalue_contributions() , file=fd)
                if np.isnan(calc_exx.get_exx_energy()):
                    print ('The EXX and therefore total energy is not be calculated, because we are only', file=fd)
                    print ('interested in a few eigenvalues for a few k-points.', file=fd)
                else:
                    print('EXX Energy: ',calc_exx.get_exx_energy() , file=fd)
                    print('Total Energy: ',calc_exx.get_total_energy() , file=fd)
            parprint('\033[94mINFORMATION:\033[0mEXX mode results are only listed at: '+struct+'-1-Result-Ground-EXX_mode.txt')
            parprint('          Other files (DOS, band, etc...) are the results calculated with PBE.')

    elif Mode == 'PW-GW':
        if restart == False:
            # PW Ground State Calculations
            parprint("Starting PW only ground state calculation for GW calculation...")
            # Fix the spacegroup in the geometric optimization if wanted
            if Fix_symmetry == True:
                    bulk_configuration.set_constraint(FixSymmetry(bulk_configuration))
            if 'kpts_density' in globals():
                calc = GPAW(mode=PW(cut_off_energy), xc=XC_calc, parallel={'domain': 1}, kpts={'density': kpts_density, 'gamma': Gamma},
                        convergence = Ground_convergence, 
                        mixer=Mixer_type, occupations = Occupation, txt=struct+'-1-Log-Ground.txt')
            else:
                calc = GPAW(mode=PW(cut_off_energy), xc=XC_calc, parallel={'domain': 1}, kpts={'size':(kpts_x, kpts_y, kpts_z), 'gamma': Gamma}, 
                        convergence = Ground_convergence, 
                        mixer=Mixer_type, occupations = Occupation, txt=struct+'-1-Log-Ground.txt')
            bulk_configuration.calc = calc
            uf = ExpCellFilter(bulk_configuration, mask=whichstrain)
            # Optimizer Selection
            if Optimizer == 'FIRE':
                from ase.optimize.fire import FIRE
                relax = FIRE(uf, maxstep=Max_step, trajectory=struct+'-1-Result-Ground.traj')
            elif  Optimizer == 'LBFGS':
                from ase.optimize.lbfgs import LBFGS
                relax = LBFGS(uf, maxstep=Max_step, alpha=Alpha, damping=Damping, trajectory=struct+'-1-Result-Ground.traj')
            elif  Optimizer == 'GPMin':
                from ase.optimize import GPMin
                relax = GPMin(uf, trajectory=struct+'-1-Result-Ground.traj')
            else:
                relax = QuasiNewton(uf, maxstep=Max_step, trajectory=struct+'-1-Result-Ground.traj')       
            relax.run(fmax=fmaxval)  # Consider tighter fmax!
            bulk_configuration.get_potential_energy()
            calc.diagonalize_full_hamiltonian()
            calc.write(struct+'-1-Result-Ground.gpw', mode="all")
            # Writes final configuration as CIF file
            write_cif(struct+'-Final.cif', bulk_configuration)
            # Print final spacegroup information
            parprint("Final Spacegroup (SPGlib):",spg.get_spacegroup(bulk_configuration))
        else:
            parprint("Passing ground state calculation for GW calculation...")
            # Control the ground state GPW file
            if not os.path.exists(struct+'-1-Result-Ground.gpw'):
                parprint('\033[91mERROR:\033[0m'+struct+'-1-Result-Ground.gpw file can not be found. It is needed in restart mode. Quiting.')
                quit()

        # We start by setting up a G0W0 calculator object
        gw = G0W0(struct+'-1-Result-Ground.gpw', filename=struct+'-1-', bands=(GWbandVB, GWbandCB), 
                  method=GWtype,truncation=GWtruncation, nblocksmax=GWnblock,
                  maxiter=5, q0_correction=GWq0correction,
                  mixing=0.5,savepckl=True,
                  ecut=GWcut_off_energy, ppa=GWppa)
        parprint("Starting PW ground state calculation with G0W0 approximation...")
        gw.calculate()
        results = pickle.load(open(struct+'-1-_results.pckl', 'rb'))
        with paropen(struct+'-1-BandGap.txt', "w") as fd:
            print('Quasi particle (QP) energies in eV. Take CB-VB for the bandgap', file=fd)
            print('To see other energy contributions, use python -mpickle <picklefile>', file=fd)
            for x in zip(results['qp']):
                    print(*x, sep=", ", file=fd)

    elif Mode == 'LCAO':
        if Spin_calc == True:
           numm = [Magmom_per_atom]*bulk_configuration.get_global_number_of_atoms()
           bulk_configuration.set_initial_magnetic_moments(numm)
        if restart == False:
            parprint("Starting LCAO ground state calculation...")
            # Fix the spacegroup in the geometric optimization if wanted
            if Fix_symmetry == True:
                    bulk_configuration.set_constraint(FixSymmetry(bulk_configuration))
            if 'gpts_density' in globals():
                if 'kpts_density' in globals():
                    calc = GPAW(mode='lcao', basis='dzp', setups= Hubbard, kpts={'density': kpts_density, 'gamma': Gamma},
                            convergence = Ground_convergence, h=gpts_density, spinpol=Spin_calc, txt=struct+'-1-Log-Ground.txt',
                            mixer=Mixer_type, occupations = Occupation, parallel={'domain': world.size})
                else:
                    calc = GPAW(mode='lcao', basis='dzp', setups= Hubbard, kpts={'size':(kpts_x, kpts_y, kpts_z), 'gamma': Gamma},
                            convergence = Ground_convergence, h=gpts_density, spinpol=Spin_calc, txt=struct+'-1-Log-Ground.txt',
                            mixer=Mixer_type, occupations = Occupation, parallel={'domain': world.size})
            else:
                if 'kpts_density' in globals():
                    calc = GPAW(mode='lcao', basis='dzp', setups= Hubbard, kpts={'density': kpts_density, 'gamma': Gamma},
                            convergence = Ground_convergence, gpts=(gpts_x, gpts_y, gpts_z), spinpol=Spin_calc, txt=struct+'-1-Log-Ground.txt',
                            mixer=Mixer_type, occupations = Occupation, parallel={'domain': world.size})
                else:
                    calc = GPAW(mode='lcao', basis='dzp', setups= Hubbard, kpts={'size':(kpts_x, kpts_y, kpts_z), 'gamma': Gamma},
                            convergence = Ground_convergence, gpts=(gpts_x, gpts_y, gpts_z), spinpol=Spin_calc, txt=struct+'-1-Log-Ground.txt',
                            mixer=Mixer_type, occupations = Occupation, parallel={'domain': world.size})
            bulk_configuration.calc = calc
            if Geo_optim == True:
                if True in whichstrain:
                    #uf = ExpCellFilter(bulk_configuration, mask=whichstrain)
                    #relax = LBFGS(uf, maxstep=Max_step, alpha=Alpha, damping=Damping, trajectory=struct+'-1-Result-Ground.traj')
                    parprint('\033[91mERROR:\033[0mModifying supercell and atom positions with a filter (whichstrain keyword) is not implemented in LCAO mode.')
                    quit()
                else:
                    # Optimizer Selection
                    if Optimizer == 'FIRE':
                        from ase.optimize.fire import FIRE
                        relax = FIRE(bulk_configuration, maxstep=Max_step, trajectory=struct+'-1-Result-Ground.traj')
                    elif  Optimizer == 'LBFGS':
                        from ase.optimize.lbfgs import LBFGS
                        relax = LBFGS(bulk_configuration, maxstep=Max_step, alpha=Alpha, damping=Damping, trajectory=struct+'-1-Result-Ground.traj')
                    elif  Optimizer == 'GPMin':
                        from ase.optimize import GPMin
                        relax = GPMin(bulk_configuration, trajectory=struct+'-1-Result-Ground.traj')
                    else:
                        relax = QuasiNewton(bulk_configuration, maxstep=Max_step, trajectory=struct+'-1-Result-Ground.traj')       
                relax.run(fmax=fmaxval)  # Consider tighter fmax!
            else:
                bulk_configuration.set_calculator(calc)
                bulk_configuration.get_potential_energy()
            #relax = LBFGS(bulk_configuration, maxstep=Max_step, alpha=Alpha, damping=Damping, trajectory=struct+'-1-Result-Ground.traj')
            #relax.run(fmax=fmaxval)  # Consider much tighter fmax!
            #bulk_configuration.get_potential_energy()
            if Density_calc == True:
                #This line makes huge GPW files. Therefore it is better to use this if else
                calc.write(struct+'-1-Result-Ground.gpw', mode="all")
            else:
                calc.write(struct+'-1-Result-Ground.gpw')
            # Writes final configuration as CIF file
            write_cif(struct+'-Final.cif', bulk_configuration)
            # Print final spacegroup information
            parprint("Final Spacegroup (SPGlib):",spg.get_spacegroup(bulk_configuration))
        else:
            parprint("Passing LCAO ground state calculation...")
            # Control the ground state GPW file
            if not os.path.exists(struct+'-1-Result-Ground.gpw'):
                parprint('\033[91mERROR:\033[0m'+struct+'-1-Result-Ground.gpw file can not be found. It is needed in restart mode. Quiting.')
                quit()

    elif Mode == 'FD':
        parprint("\033[91mERROR:\033[0mFD mode is not implemented in gpaw-tools yet...")
        quit()
    else:
        parprint("\033[91mERROR:\033[0mPlease enter correct mode information.")
        quit()
    # Finish ground state
    time12 = time.time()
    # -------------------------------------------------------------
    # Step 1.5 - ELASTIC CALCULATION
    # -------------------------------------------------------------
    if Elastic_calc == True:
        # Start elastic calc
        time151 = time.time()
        parprint('Starting elastic tensor calculations (\033[93mWARNING:\033[0mNOT TESTED FEATURE, PLEASE CONTROL THE RESULTS)...')
        calc = GPAW(struct+'-1-Result-Ground.gpw', fixdensity=True, txt=struct+'-1.5-Log-Elastic.txt')
        # Getting space group from SPGlib
        parprint('Spacegroup:',spg.get_spacegroup(bulk_configuration))
        # Calculating equation of state
        parprint('Calculating equation of state...')
        eos = calculate_eos(bulk_configuration, trajectory=struct+'-1-Result-Ground.traj')
        v, e, B = eos.fit()
        # Calculating elastic tensor
        parprint('Calculating elastic tensor...')
        Cij, Bij=get_elastic_tensor(bulk_configuration,get_elementary_deformations(bulk_configuration, n=5, d=2))
        with paropen(struct+'-1.5-Result-Elastic.txt', "w") as fd:
            print("Elastic calculation results (NOT TESTED FEATURE, PLEASE CONTROL THE RESULTS):", file=fd)
            print("EoS: Stabilized jellium equation of state (SJEOS)", file=fd)
            print("Refs: Phys.Rev.B 63, 224115 (2001) and Phys.Rev.B 67, 026103 (2003)", file=fd)
            print("Elastic constants: Standart elasticity theory calculated by -Elastic- library", file=fd)
            print("Ref: European Physical Journal B; 15, 2 (2000) 265-268", file=fd)
            print("-----------------------------------------------------------------------------", file=fd)
            print("Spacegroup: "+str(spg.get_spacegroup(bulk_configuration)), file=fd)
            print("B (GPa): "+str(B / kJ * 1.0e24), file=fd)
            print("e (eV): "+str(e), file=fd)
            print("v (Ang^3): "+str(v), file=fd)
            print("Cij (GPa): ",Cij/GPa, file=fd)
        # Finish elastic calc
        time152 = time.time()
    # -------------------------------------------------------------
    # Step 2 - DOS CALCULATION
    # -------------------------------------------------------------
    if DOS_calc == True:
        # Start DOS calc
        time21 = time.time()
        parprint("Starting DOS calculation...")
        calc = GPAW(struct+'-1-Result-Ground.gpw', fixdensity=True, txt=struct+'-2-Log-DOS.txt',
                convergence = Ground_convergence, occupations = Occupation)
        #energies, weights = calc.get_dos(npts=800, width=0)
        dos = DOS(calc, npts=DOS_npoints, width=DOS_width)
        
        if Spin_calc == True:
            energies = dos.get_energies()
            weights = dos.get_dos(spin=0)
            weightsup = dos.get_dos(spin=1)
        else:
            energies = dos.get_energies()
            weights = dos.get_dos()

        with paropen(struct+'-2-Result-DOS.csv', "w") as fd:
            if Spin_calc == True:
                for x in zip(energies, weights, weightsup):
                    print(*x, sep=", ", file=fd)
            else:
                for x in zip(energies, weights):
                    print(*x, sep=", ", file=fd)
        #PDOS
        parprint("Calculating and saving PDOS...")
        chem_sym = bulk_configuration.get_chemical_symbols()
        ef = calc.get_fermi_level()
        
        if Spin_calc == True:
            #Spin down
            with paropen(struct+'-2-Result-PDOS-Down.csv', "w") as fd:
                print("Energy, s-orbital, p-orbital, d-orbital, f-orbital", file=fd)
                for j in range(0, bulk_configuration.get_global_number_of_atoms()):
                    print("Atom no: "+str(j+1)+", Atom Symbol: "+chem_sym[j]+" --------------------", file=fd)
                    en, pdossd = calc.get_orbital_ldos(a=j, spin=0, angular='s', npts=DOS_npoints, width=DOS_width)
                    en, pdospd = calc.get_orbital_ldos(a=j, spin=0, angular='p', npts=DOS_npoints, width=DOS_width)
                    en, pdosdd = calc.get_orbital_ldos(a=j, spin=0, angular='d', npts=DOS_npoints, width=DOS_width)
                    en, pdosfd = calc.get_orbital_ldos(a=j, spin=0, angular='f', npts=DOS_npoints, width=DOS_width)
                    for x in zip(en-ef, pdossd, pdospd, pdosdd, pdosfd):
                        print(*x, sep=", ", file=fd)
                print("---------------------------------------------------- --------------------", file=fd)
                
            # RAW PDOS for spin down
            parprint("Calculating and saving Raw PDOS for spin down...")
            rawdos = DOSCalculator.from_calculator(struct+'-1-Result-Ground.gpw',soc=False, theta=0.0, phi=0.0, shift_fermi_level=True)
            energies = rawdos.get_energies(npoints=DOS_npoints)
              
            with paropen(struct+'-2-Result-RawPDOS-Down.csv', "w") as fd:
                print("Energy, s-total, p-total, px, py, pz, d-total, dxy, dyz, d3z2_r2, dzx, dx2_y2, f-total", file=fd)
                for j in range(0, bulk_configuration.get_global_number_of_atoms()):
                    print("Atom no: "+str(j+1)+", Atom Symbol: "+chem_sym[j]+" --------------------", file=fd)
                    pdoss = rawdos.raw_pdos(energies, a=j, l=0, m=None, spin=0, width=DOS_width)
                    pdosp = rawdos.raw_pdos(energies, a=j, l=1, m=None, spin=0, width=DOS_width)
                    pdospx = rawdos.raw_pdos(energies, a=j, l=1, m=2, spin=0, width=DOS_width)
                    pdospy = rawdos.raw_pdos(energies, a=j, l=1, m=0, spin=0, width=DOS_width)
                    pdospz = rawdos.raw_pdos(energies, a=j, l=1, m=1, spin=0, width=DOS_width)
                    pdosd = rawdos.raw_pdos(energies, a=j, l=2, m=None, spin=0, width=DOS_width)
                    pdosdxy = rawdos.raw_pdos(energies, a=j, l=2, m=0, spin=0, width=DOS_width)
                    pdosdyz = rawdos.raw_pdos(energies, a=j, l=2, m=1, spin=0, width=DOS_width)
                    pdosd3z2_r2 = rawdos.raw_pdos(energies, a=j, l=2, m=2, spin=0, width=DOS_width)
                    pdosdzx = rawdos.raw_pdos(energies, a=j, l=2, m=3, spin=0, width=DOS_width)
                    pdosdx2_y2 = rawdos.raw_pdos(energies, a=j, l=2, m=4, spin=0, width=DOS_width)
                    pdosf = rawdos.raw_pdos(energies, a=j, l=3, m=None, spin=0, width=DOS_width)
                    for x in zip(en-ef, pdoss, pdosp, pdospx, pdospy, pdospz, pdosd, pdosdxy, pdosdyz, pdosd3z2_r2, pdosdzx, pdosdx2_y2, pdosf):
                        print(*x, sep=", ", file=fd)

            #Spin up
            with paropen(struct+'-2-Result-PDOS-Up.csv', "w") as fd:
                print("Energy, s-orbital, p-orbital, d-orbital, f-orbital", file=fd)
                for j in range(0, bulk_configuration.get_global_number_of_atoms()):
                    print("Atom no: "+str(j+1)+", Atom Symbol: "+chem_sym[j]+" --------------------", file=fd)
                    en, pdossu = calc.get_orbital_ldos(a=j, spin=1, angular='s', npts=DOS_npoints, width=DOS_width)
                    en, pdospu = calc.get_orbital_ldos(a=j, spin=1, angular='p', npts=DOS_npoints, width=DOS_width)
                    en, pdosdu = calc.get_orbital_ldos(a=j, spin=1, angular='d', npts=DOS_npoints, width=DOS_width)
                    en, pdosfu = calc.get_orbital_ldos(a=j, spin=1, angular='f', npts=DOS_npoints, width=DOS_width)
                    for x in zip(en-ef, pdossu, pdospu, pdosdu, pdosfu):
                        print(*x, sep=", ", file=fd)
                print("---------------------------------------------------- --------------------", file=fd)

            # RAW PDOS for spin up
            parprint("Calculating and saving Raw PDOS for spin up...")
            rawdos = DOSCalculator.from_calculator(struct+'-1-Result-Ground.gpw',soc=False, theta=0.0, phi=0.0, shift_fermi_level=True)
            energies = rawdos.get_energies(npoints=DOS_npoints)

            with paropen(struct+'-2-Result-RawPDOS-Up.csv', "w") as fd:
                print("Energy, s-total, p-total, px, py, pz, d-total, dxy, dyz, d3z2_r2, dzx, dx2_y2, f-total", file=fd)
                for j in range(0, bulk_configuration.get_global_number_of_atoms()):
                    print("Atom no: "+str(j+1)+", Atom Symbol: "+chem_sym[j]+" --------------------", file=fd)
                    pdoss = rawdos.raw_pdos(energies, a=j, l=0, m=None, spin=1, width=DOS_width)
                    pdosp = rawdos.raw_pdos(energies, a=j, l=1, m=None, spin=1, width=DOS_width)
                    pdospx = rawdos.raw_pdos(energies, a=j, l=1, m=2, spin=1, width=DOS_width)
                    pdospy = rawdos.raw_pdos(energies, a=j, l=1, m=0, spin=1, width=DOS_width)
                    pdospz = rawdos.raw_pdos(energies, a=j, l=1, m=1, spin=1, width=DOS_width)
                    pdosd = rawdos.raw_pdos(energies, a=j, l=2, m=None, spin=1, width=DOS_width)
                    pdosdxy = rawdos.raw_pdos(energies, a=j, l=2, m=0, spin=1, width=DOS_width)
                    pdosdyz = rawdos.raw_pdos(energies, a=j, l=2, m=1, spin=1, width=DOS_width)
                    pdosd3z2_r2 = rawdos.raw_pdos(energies, a=j, l=2, m=2, spin=1, width=DOS_width)
                    pdosdzx = rawdos.raw_pdos(energies, a=j, l=2, m=3, spin=1, width=DOS_width)
                    pdosdx2_y2 = rawdos.raw_pdos(energies, a=j, l=2, m=4, spin=1, width=DOS_width)
                    pdosf = rawdos.raw_pdos(energies, a=j, l=3, m=None, spin=1, width=DOS_width)
                    for x in zip(en-ef, pdoss, pdosp, pdospx, pdospy, pdospz, pdosd, pdosdxy, pdosdyz, pdosd3z2_r2, pdosdzx, pdosdx2_y2, pdosf):
                        print(*x, sep=", ", file=fd)

        else:
            with paropen(struct+'-2-Result-PDOS.csv', "w") as fd:
                print("Energy, s-orbital, p-orbital, d-orbital, f-orbital", file=fd)
                for j in range(0, bulk_configuration.get_global_number_of_atoms()):
                    print("Atom no: "+str(j+1)+", Atom Symbol: "+chem_sym[j]+" --------------------", file=fd)
                    en, pdoss = calc.get_orbital_ldos(a=j, spin=0, angular='s', npts=DOS_npoints, width=DOS_width)
                    en, pdosp = calc.get_orbital_ldos(a=j, spin=0, angular='p', npts=DOS_npoints, width=DOS_width)
                    en, pdosd = calc.get_orbital_ldos(a=j, spin=0, angular='d', npts=DOS_npoints, width=DOS_width)
                    en, pdosf = calc.get_orbital_ldos(a=j, spin=0, angular='f', npts=DOS_npoints, width=DOS_width)
                    for x in zip(en-ef, pdoss, pdosp, pdosd, pdosf):
                        print(*x, sep=", ", file=fd)
                print("---------------------------------------------------- --------------------", file=fd)
            # RAW PDOS  
            parprint("Calculating and saving Raw PDOS...")
            rawdos = DOSCalculator.from_calculator(struct+'-1-Result-Ground.gpw',soc=False, theta=0.0, phi=0.0, shift_fermi_level=True)
            energies = rawdos.get_energies(npoints=DOS_npoints)
              
            with paropen(struct+'-2-Result-RawPDOS.csv', "w") as fd:
                print("Energy, s-total, p-total, px, py, pz, d-total, dxy, dyz, d3z2_r2, dzx, dx2_y2, f-total", file=fd)
                for j in range(0, bulk_configuration.get_global_number_of_atoms()):
                    print("Atom no: "+str(j+1)+", Atom Symbol: "+chem_sym[j]+" --------------------", file=fd)
                    pdoss = rawdos.raw_pdos(energies, a=j, l=0, m=None, spin=None, width=DOS_width)
                    pdosp = rawdos.raw_pdos(energies, a=j, l=1, m=None, spin=None, width=DOS_width)
                    pdospx = rawdos.raw_pdos(energies, a=j, l=1, m=2, spin=None, width=DOS_width)
                    pdospy = rawdos.raw_pdos(energies, a=j, l=1, m=0, spin=None, width=DOS_width)
                    pdospz = rawdos.raw_pdos(energies, a=j, l=1, m=1, spin=None, width=DOS_width)
                    pdosd = rawdos.raw_pdos(energies, a=j, l=2, m=None, spin=None, width=DOS_width)
                    pdosdxy = rawdos.raw_pdos(energies, a=j, l=2, m=0, spin=None, width=DOS_width)
                    pdosdyz = rawdos.raw_pdos(energies, a=j, l=2, m=1, spin=None, width=DOS_width)
                    pdosd3z2_r2 = rawdos.raw_pdos(energies, a=j, l=2, m=2, spin=None, width=DOS_width)
                    pdosdzx = rawdos.raw_pdos(energies, a=j, l=2, m=3, spin=None, width=DOS_width)
                    pdosdx2_y2 = rawdos.raw_pdos(energies, a=j, l=2, m=4, spin=None, width=DOS_width)
                    pdosf = rawdos.raw_pdos(energies, a=j, l=3, m=None, spin=None, width=DOS_width)
                    for x in zip(en-ef, pdoss, pdosp, pdospx, pdospy, pdospz, pdosd, pdosdxy, pdosdyz, pdosd3z2_r2, pdosdzx, pdosdx2_y2, pdosf):
                        print(*x, sep=", ", file=fd)
                print("---------------------------------------------------- --------------------", file=fd)
        # Finish DOS calc
        time22 = time.time()
    # -------------------------------------------------------------
    # Step 3 - BAND STRUCTURE CALCULATION
    # -------------------------------------------------------------
    if Band_calc == True:
        # Start Band calc
        time31 = time.time()
        parprint("Starting band structure calculation...")
        if Mode == 'PW-GW':      
            GW = GWBands(calc=struct+'-1-Result-Ground.gpw', fixdensity=True,
                 gw_file=struct+'-1-_results.pckl',kpoints=GWkpoints)

            # Getting results without spin-orbit
            results = GW.get_gw_bands(SO=False, interpolate=GWbandinterpolation, vac=True)

            # Extracting data
            X = results['X']
            ef = results['ef']
            xdata = results['x_k']
            banddata = results['e_kn']

            np.savetxt(struct+'-3-Result-Band.dat', np.c_[xdata,banddata])

            with open(struct+'-3-Result-Band.dat', 'a') as f:
                print ('Symmetry points: ', X, end="\n", file=f)
                print ('Fermi Level: ', ef, end="\n", file=f)

        else:
            if XC_calc in ['HSE06', 'HSE03']:
                calc = GPAW(struct+'-1-Result-Ground.gpw',
                        parallel={'band': 1, 'kpt': 1}, 
                        txt=struct+'-3-Log-Band.txt',
                        symmetry='off', occupations = Occupation,
                        kpts={'path': band_path, 'npoints': band_npoints}, convergence=Band_convergence)
                
            else:
                calc = GPAW(struct+'-1-Result-Ground.gpw',
                        txt=struct+'-3-Log-Band.txt',
                        fixdensity=True,
                        symmetry='off', occupations = Occupation,
                        kpts={'path': band_path, 'npoints': band_npoints},
                        convergence=Band_convergence)

            calc.get_potential_energy()
            bs = calc.band_structure()
            ef = calc.get_fermi_level()
            num_of_bands = calc.get_number_of_bands()
            parprint('Num of bands:'+str(num_of_bands))

            # No need to write an additional gpaw file. Use json file to use with ase band-structure command
            #calc.write(struct+'-3-Result-Band.gpw')
            bs.write(struct+'-3-Result-Band.json')

            if Spin_calc == True:
                eps_skn = np.array([[calc.get_eigenvalues(k,s)
                                    for k in range(band_npoints)]
                                    for s in range(2)]) - ef
                parprint(eps_skn.shape)
                with paropen(struct+'-3-Result-Band-Down.dat', 'w') as f1:
                    for n1 in range(num_of_bands):
                        for k1 in range(band_npoints):
                            print(k1, eps_skn[0, k1, n1], end="\n", file=f1)
                        print (end="\n", file=f1)

                with paropen(struct+'-3-Result-Band-Up.dat', 'w') as f2:
                    for n2 in range(num_of_bands):
                        for k2 in range(band_npoints):
                            print(k2, eps_skn[1, k2, n2], end="\n", file=f2)
                        print (end="\n", file=f2)

            else:
                eps_skn = np.array([[calc.get_eigenvalues(k,s)
                                    for k in range(band_npoints)]
                                    for s in range(1)]) - ef
                with paropen(struct+'-3-Result-Band.dat', 'w') as f:
                    for n in range(num_of_bands):
                        for k in range(band_npoints):
                            print(k, eps_skn[0, k, n], end="\n", file=f)
                        print (end="\n", file=f)
        # Finish Band calc
        time32 = time.time()
    # -------------------------------------------------------------
    # Step 4 - ALL-ELECTRON DENSITY
    # -------------------------------------------------------------
    if Density_calc == True:
        #Start Density calc
        time41 = time.time()
        parprint("Starting All-electron density calculation...")
        calc = GPAW(struct+'-1-Result-Ground.gpw', txt=struct+'-4-Log-ElectronDensity.txt')
        bulk_configuration.calc = calc
        np = calc.get_pseudo_density()
        n = calc.get_all_electron_density(gridrefinement=gridref)
        
        # Writing pseudo and all electron densities to cube file with Bohr unit
        write(struct+'-4-Result-All-electron_nall.cube', bulk_configuration, data=n * Bohr**3)
        write(struct+'-4-Result-All-electron_npseudo.cube', bulk_configuration, data=np * Bohr**3)
        # Finish Density calc
        time42 = time.time()

# -------------------------------------------------------------
# Step 5 - OPTICAL CALCULATION
# -------------------------------------------------------------
if Optical_calc == True:
    #Start Optical calc
    time51 = time.time()
    if Mode == 'PW':
        parprint("Starting optical calculation...")
        try:
            calc = GPAW(struct+'-1-Result-Ground.gpw',
                    txt=struct+'-5-Log-Optical.txt',
                    nbands=num_of_bands,
                    fixdensity=True,
                    symmetry='off',
                    occupations=FermiDirac(optFDsmear))
        except FileNotFoundError as err:
            # output error, and return with an error code
            parprint('\033[91mERROR:\033[0mOptical computations must be done separately. Please do ground calculations first.')
            quit()
    
        calc.get_potential_energy()

        calc.diagonalize_full_hamiltonian(nbands=num_of_bands)  # diagonalize Hamiltonian
        calc.write(struct+'-5-Result-Optical.gpw', 'all')  # write wavefunctions

        #from mpi4py import MPI
        if opttype == 'BSE':
            if Spin_calc == True:
               parprint('\033[91mERROR:\033[0mBSE calculations can not run with spin dependent data.')
               quit()
            parprint('Starting BSE calculations')
            bse = BSE(calc= struct+'-5-Result-Optical.gpw', ecut=optecut,
                         valence_bands=optBSEvb,
                         conduction_bands=optBSEcb,
                         nbands=num_of_bands,
                         eshift=optshift,
                         mode='BSE',
                         write_v=True,
                         integrate_gamma=0,
                         txt=struct+'-5-Log-Optical-BSE.txt')
            
            # Getting dielectric function spectrum
            parprint("Starting dielectric function calculation...")
            # Writing to files
            bse.get_dielectric_function(direction=0, q_c = [0.0, 0.0, 0.0], eta=opteta,
                                        w_w=np.linspace(optBSEminEn, optBSEmaxEn, optBSEnumdata),
                                        filename=struct+'-5-Result-Optical-BSE_dielec_xdirection.csv',
                                        write_eig=struct+'-5-Result-Optical-BSE_eig_xdirection.dat')
            bse.get_dielectric_function(direction=1, q_c = [0.0, 0.0, 0.0], eta=opteta,
                                        w_w=np.linspace(optBSEminEn, optBSEmaxEn, optBSEnumdata),
                                        filename=struct+'-5-Result-Optical-BSE_dielec_ydirection.csv',
                                        write_eig=struct+'-5-Result-Optical-BSE_eig_ydirection.dat')
            bse.get_dielectric_function(direction=2, q_c = [0.0, 0.0, 0.0], eta=opteta,
                                        w_w=np.linspace(optBSEminEn, optBSEmaxEn, optBSEnumdata),
                                        filename=struct+'-5-Result-Optical-BSE_dielec_zdirection.csv',
                                        write_eig=struct+'-5-Result-Optical-BSE_eig_zdirection.dat')
                                        
            # Loading dielectric function spectrum to numpy
            dielec_x = genfromtxt(struct+'-5-Result-Optical-BSE_dielec_xdirection.csv', delimiter=',')
            dielec_y = genfromtxt(struct+'-5-Result-Optical-BSE_dielec_ydirection.csv', delimiter=',')
            dielec_z = genfromtxt(struct+'-5-Result-Optical-BSE_dielec_zdirection.csv', delimiter=',')
            # dielec_x.shape[0] will give us the length of data.
            # c and h
            c_opt = 29979245800
            h_opt = 6.58E-16
            #c_opt = 1
            #h_opt = 1

            # Initialize arrays
            opt_n_bse_x = np.array ([1e-6,]*dielec_x.shape[0])
            opt_k_bse_x = np.array ([1e-6,]*dielec_x.shape[0])
            opt_abs_bse_x = np.array([1e-6,]*dielec_x.shape[0])
            opt_ref_bse_x = np.array([1e-6,]*dielec_x.shape[0])
            opt_n_bse_y = np.array ([1e-6,]*dielec_y.shape[0])
            opt_k_bse_y = np.array ([1e-6,]*dielec_y.shape[0])
            opt_abs_bse_y = np.array([1e-6,]*dielec_y.shape[0])
            opt_ref_bse_y = np.array([1e-6,]*dielec_y.shape[0])
            opt_n_bse_z = np.array ([1e-6,]*dielec_z.shape[0])
            opt_k_bse_z = np.array ([1e-6,]*dielec_z.shape[0])
            opt_abs_bse_z = np.array([1e-6,]*dielec_z.shape[0])
            opt_ref_bse_z = np.array([1e-6,]*dielec_z.shape[0])

            # Calculation of other optical data
            for n in range(dielec_x.shape[0]):
                # x-direction
                opt_n_bse_x[n] = np.sqrt((np.sqrt(np.square(dielec_x[n][1])+np.square(dielec_x[n][2]))+dielec_x[n][1])/2.0)
                opt_k_bse_x[n] = np.sqrt((np.sqrt(np.square(dielec_x[n][1])+np.square(dielec_x[n][2]))-dielec_x[n][1])/2.0)
                opt_abs_bse_x[n] = 2*dielec_x[n][0]*opt_k_bse_x[n]/(h_opt*c_opt)
                opt_ref_bse_x[n] = (np.square(1-opt_n_bse_x[n])+np.square(opt_k_bse_x[n]))/(np.square(1+opt_n_bse_x[n])+np.square(opt_k_bse_x[n]))
                # y-direction
                opt_n_bse_y[n] = np.sqrt((np.sqrt(np.square(dielec_y[n][1])+np.square(dielec_y[n][2]))+dielec_y[n][1])/2.0)
                opt_k_bse_y[n] = np.sqrt((np.sqrt(np.square(dielec_y[n][1])+np.square(dielec_y[n][2]))-dielec_y[n][1])/2.0)
                opt_abs_bse_y[n] = 2*dielec_y[n][0]*opt_k_bse_y[n]/(h_opt*c_opt)
                opt_ref_bse_y[n] = (np.square(1-opt_n_bse_y[n])+np.square(opt_k_bse_y[n]))/(np.square(1+opt_n_bse_y[n])+np.square(opt_k_bse_y[n]))
                # z-direction
                opt_n_bse_z[n] = np.sqrt((np.sqrt(np.square(dielec_z[n][1])+np.square(dielec_z[n][2]))+dielec_z[n][1])/2.0)
                opt_k_bse_z[n] = np.sqrt((np.sqrt(np.square(dielec_z[n][1])+np.square(dielec_z[n][2]))-dielec_z[n][1])/2.0)
                opt_abs_bse_z[n] = 2*dielec_z[n][0]*opt_k_bse_z[n]/(h_opt*c_opt)
                opt_ref_bse_z[n] = (np.square(1-opt_n_bse_z[n])+np.square(opt_k_bse_z[n]))/(np.square(1+opt_n_bse_z[n])+np.square(opt_k_bse_z[n]))

            # Saving other data for x-direction
            with paropen(struct+'-5-Result-Optical-BSE-AllData_xdirection.dat', 'w') as f1:
                print("Energy(eV) Eps_real Eps_img Refractive_Index Extinction_Index Absorption(1/cm) Reflectivity", end="\n", file=f1)
                for n in range(dielec_x.shape[0]):
                    print(dielec_x[n][0], dielec_x[n][1], dielec_x[n][2], opt_n_bse_x[n], opt_k_bse_x[n], opt_abs_bse_x[n], opt_ref_bse_x[n], end="\n", file=f1)
                print (end="\n", file=f1)

            # Saving other data for y-direction
            with paropen(struct+'-5-Result-Optical-BSE-AllData_ydirection.dat', 'w') as f1:
                print("Energy(eV) Eps_real Eps_img Refractive_Index Extinction_Index Absorption(1/cm) Reflectivity", end="\n", file=f1)
                for n in range(dielec_y.shape[0]):
                    print(dielec_y[n][0], dielec_y[n][1], dielec_y[n][2], opt_n_bse_y[n], opt_k_bse_y[n], opt_abs_bse_y[n], opt_ref_bse_y[n], end="\n", file=f1)
                print (end="\n", file=f1)

            # Saving other data for z-direction
            with paropen(struct+'-5-Result-Optical-BSE-AllData_zdirection.dat', 'w') as f1:
                print("Energy(eV) Eps_real Eps_img Refractive_Index Extinction_Index Absorption(1/cm) Reflectivity", end="\n", file=f1)
                for n in range(dielec_z.shape[0]):
                    print(dielec_z[n][0], dielec_z[n][1], dielec_z[n][2], opt_n_bse_z[n], opt_k_bse_z[n], opt_abs_bse_z[n], opt_ref_bse_z[n], end="\n", file=f1)
                print (end="\n", file=f1)

        elif opttype == 'RPA':
            parprint('Starting RPA calculations')
            df = DielectricFunction(calc=struct+'-5-Result-Optical.gpw',
                                eta=opteta, nblocks=world.size,
                                omega2=optomega2,
                                domega0=optdomega0,
                                ecut=optecut)
            # Writing to files as: omega, nlfc.real, nlfc.imag, lfc.real, lfc.imag 
            # Here lfc is local field correction
            # Getting dielectric function spectrum
            parprint("Starting dielectric function calculation...")
            df.get_dielectric_function( direction='x', 
                                        filename=struct+'-5-Result-Optical-RPA_dielec_xdirection.csv')
            df.get_dielectric_function( direction='y',
                                        filename=struct+'-5-Result-Optical-RPA_dielec_ydirection.csv')
            df.get_dielectric_function( direction='z',
                                        filename=struct+'-5-Result-Optical-RPA_dielec_zdirection.csv')

            # Loading dielectric function spectrum to numpy
            dielec_x = genfromtxt(struct+'-5-Result-Optical-RPA_dielec_xdirection.csv', delimiter=',')
            dielec_y = genfromtxt(struct+'-5-Result-Optical-RPA_dielec_ydirection.csv', delimiter=',')
            dielec_z = genfromtxt(struct+'-5-Result-Optical-RPA_dielec_zdirection.csv', delimiter=',')
            # dielec_x.shape[0] will give us the length of data.
            # c and h
            c_opt = 29979245800
            h_opt = 6.58E-16
            #c_opt = 1
            #h_opt = 1
            # ---- NLFC ----
            # Initialize arrays for NLFC
            opt_n_nlfc_x = np.array ([1e-6,]*dielec_x.shape[0])
            opt_k_nlfc_x = np.array ([1e-6,]*dielec_x.shape[0])
            opt_abs_nlfc_x = np.array([1e-6,]*dielec_x.shape[0])
            opt_ref_nlfc_x = np.array([1e-6,]*dielec_x.shape[0])
            opt_n_nlfc_y = np.array ([1e-6,]*dielec_y.shape[0])
            opt_k_nlfc_y = np.array ([1e-6,]*dielec_y.shape[0])
            opt_abs_nlfc_y = np.array([1e-6,]*dielec_y.shape[0])
            opt_ref_nlfc_y = np.array([1e-6,]*dielec_y.shape[0])
            opt_n_nlfc_z = np.array ([1e-6,]*dielec_z.shape[0])
            opt_k_nlfc_z = np.array ([1e-6,]*dielec_z.shape[0])
            opt_abs_nlfc_z = np.array([1e-6,]*dielec_z.shape[0])
            opt_ref_nlfc_z = np.array([1e-6,]*dielec_z.shape[0])

            # Calculation of other optical spectrum for NLFC
            for n in range(dielec_x.shape[0]):
                # NLFC-x
                opt_n_nlfc_x[n] = np.sqrt((np.sqrt(np.square(dielec_x[n][1])+np.square(dielec_x[n][2]))+dielec_x[n][1])/2.0)
                opt_k_nlfc_x[n] = np.sqrt((np.sqrt(np.square(dielec_x[n][1])+np.square(dielec_x[n][2]))-dielec_x[n][1])/2.0)
                opt_abs_nlfc_x[n] = 2*dielec_x[n][0]*opt_k_nlfc_x[n]/(h_opt*c_opt)
                opt_ref_nlfc_x[n] = (np.square(1-opt_n_nlfc_x[n])+np.square(opt_k_nlfc_x[n]))/(np.square(1+opt_n_nlfc_x[n])+np.square(opt_k_nlfc_x[n]))
                # NLFC-y
                opt_n_nlfc_y[n] = np.sqrt((np.sqrt(np.square(dielec_y[n][1])+np.square(dielec_y[n][2]))+dielec_y[n][1])/2.0)
                opt_k_nlfc_y[n] = np.sqrt((np.sqrt(np.square(dielec_y[n][1])+np.square(dielec_y[n][2]))-dielec_y[n][1])/2.0)
                opt_abs_nlfc_y[n] = 2*dielec_y[n][0]*opt_k_nlfc_y[n]/(h_opt*c_opt)
                opt_ref_nlfc_y[n] = (np.square(1-opt_n_nlfc_y[n])+np.square(opt_k_nlfc_y[n]))/(np.square(1+opt_n_nlfc_y[n])+np.square(opt_k_nlfc_y[n]))
                # NLFC-z
                opt_n_nlfc_z[n] = np.sqrt((np.sqrt(np.square(dielec_z[n][1])+np.square(dielec_z[n][2]))+dielec_z[n][1])/2.0)
                opt_k_nlfc_z[n] = np.sqrt((np.sqrt(np.square(dielec_z[n][1])+np.square(dielec_z[n][2]))-dielec_z[n][1])/2.0)
                opt_abs_nlfc_z[n] = 2*dielec_z[n][0]*opt_k_nlfc_z[n]/(h_opt*c_opt)
                opt_ref_nlfc_z[n] = (np.square(1-opt_n_nlfc_z[n])+np.square(opt_k_nlfc_z[n]))/(np.square(1+opt_n_nlfc_z[n])+np.square(opt_k_nlfc_z[n]))

            # Saving NLFC other optical spectrum for x-direction
            with paropen(struct+'-5-Result-Optical-RPA-NLFC-AllData_xdirection.dat', 'w') as f1:
                print("Energy(eV) Eps_real Eps_img Refractive_Index Extinction_Index Absorption(1/cm) Reflectivity", end="\n", file=f1)
                for n in range(dielec_x.shape[0]):
                    print(dielec_x[n][0], dielec_x[n][1], dielec_x[n][2], opt_n_nlfc_x[n], opt_k_nlfc_x[n], opt_abs_nlfc_x[n], opt_ref_nlfc_x[n], end="\n", file=f1)
                print (end="\n", file=f1)

            # Saving NLFC other optical spectrum for y-direction
            with paropen(struct+'-5-Result-Optical-RPA-NLFC-AllData_ydirection.dat', 'w') as f1:
                print("Energy(eV) Eps_real Eps_img Refractive_Index Extinction_Index Absorption(1/cm) Reflectivity", end="\n", file=f1)
                for n in range(dielec_y.shape[0]):
                    print(dielec_y[n][0], dielec_y[n][1], dielec_y[n][2], opt_n_nlfc_y[n], opt_k_nlfc_y[n], opt_abs_nlfc_y[n], opt_ref_nlfc_y[n], end="\n", file=f1)
                print (end="\n", file=f1)

            # Saving NLFC other optical spectrum for z-direction
            with paropen(struct+'-5-Result-Optical-RPA-NLFC-AllData_zdirection.dat', 'w') as f1:
                print("Energy(eV) Eps_real Eps_img Refractive_Index Extinction_Index Absorption(1/cm) Reflectivity", end="\n", file=f1)
                for n in range(dielec_z.shape[0]):
                    print(dielec_z[n][0], dielec_z[n][1], dielec_z[n][2], opt_n_nlfc_z[n], opt_k_nlfc_z[n], opt_abs_nlfc_z[n], opt_ref_nlfc_z[n], end="\n", file=f1)
                print (end="\n", file=f1)

            # ---- LFC ----
            # Initialize arrays for LFC
            opt_n_lfc_x = np.array ([1e-6,]*dielec_x.shape[0])
            opt_k_lfc_x = np.array ([1e-6,]*dielec_x.shape[0])
            opt_abs_lfc_x = np.array([1e-6,]*dielec_x.shape[0])
            opt_ref_lfc_x = np.array([1e-6,]*dielec_x.shape[0])
            opt_n_lfc_y = np.array ([1e-6,]*dielec_y.shape[0])
            opt_k_lfc_y = np.array ([1e-6,]*dielec_y.shape[0])
            opt_abs_lfc_y = np.array([1e-6,]*dielec_y.shape[0])
            opt_ref_lfc_y = np.array([1e-6,]*dielec_y.shape[0])
            opt_n_lfc_z = np.array ([1e-6,]*dielec_z.shape[0])
            opt_k_lfc_z = np.array ([1e-6,]*dielec_z.shape[0])
            opt_abs_lfc_z = np.array([1e-6,]*dielec_z.shape[0])
            opt_ref_lfc_z = np.array([1e-6,]*dielec_z.shape[0])

            # Calculation of other optical spectrum for LFC
            for n in range(dielec_x.shape[0]):
                # LFC-x
                opt_n_lfc_x[n] = np.sqrt((np.sqrt(np.square(dielec_x[n][3])+np.square(dielec_x[n][4]))+dielec_x[n][3])/2.0)
                opt_k_lfc_x[n] = np.sqrt((np.sqrt(np.square(dielec_x[n][3])+np.square(dielec_x[n][4]))-dielec_x[n][3])/2.0)
                opt_abs_lfc_x[n] = 2*dielec_x[n][0]*opt_k_nlfc_x[n]/(h_opt*c_opt)
                opt_ref_lfc_x[n] = (np.square(1-opt_n_lfc_x[n])+np.square(opt_k_lfc_x[n]))/(np.square(1+opt_n_lfc_x[n])+np.square(opt_k_lfc_x[n]))
                # LFC-y
                opt_n_lfc_y[n] = np.sqrt((np.sqrt(np.square(dielec_y[n][3])+np.square(dielec_y[n][4]))+dielec_y[n][3])/2.0)
                opt_k_lfc_y[n] = np.sqrt((np.sqrt(np.square(dielec_y[n][3])+np.square(dielec_y[n][4]))-dielec_y[n][3])/2.0)
                opt_abs_lfc_y[n] = 2*dielec_y[n][0]*opt_k_lfc_y[n]/(h_opt*c_opt)
                opt_ref_lfc_y[n] = (np.square(1-opt_n_lfc_y[n])+np.square(opt_k_lfc_y[n]))/(np.square(1+opt_n_lfc_y[n])+np.square(opt_k_lfc_y[n]))
                # LFC-z
                opt_n_lfc_z[n] = np.sqrt((np.sqrt(np.square(dielec_z[n][3])+np.square(dielec_z[n][4]))+dielec_z[n][3])/2.0)
                opt_k_lfc_z[n] = np.sqrt((np.sqrt(np.square(dielec_z[n][3])+np.square(dielec_z[n][4]))-dielec_z[n][3])/2.0)
                opt_abs_lfc_z[n] = 2*dielec_z[n][0]*opt_k_lfc_z[n]/(h_opt*c_opt)
                opt_ref_lfc_z[n] = (np.square(1-opt_n_lfc_z[n])+np.square(opt_k_lfc_z[n]))/(np.square(1+opt_n_lfc_z[n])+np.square(opt_k_lfc_z[n]))

            # Saving LFC other optical spectrum for x-direction
            with paropen(struct+'-5-Result-Optical-RPA-LFC-AllData_xdirection.dat', 'w') as f1:
                print("Energy(eV) Eps_real Eps_img Refractive_Index Extinction_Index Absorption(1/cm) Reflectivity", end="\n", file=f1)
                for n in range(dielec_x.shape[0]):
                    print(dielec_x[n][0], dielec_x[n][3], dielec_x[n][4], opt_n_lfc_x[n], opt_k_lfc_x[n], opt_abs_lfc_x[n], opt_ref_lfc_x[n], end="\n", file=f1)
                print (end="\n", file=f1)

            # Saving LFC other optical spectrum for y-direction
            with paropen(struct+'-5-Result-Optical-RPA-LFC-AllData_ydirection.dat', 'w') as f1:
                print("Energy(eV) Eps_real Eps_img Refractive_Index Extinction_Index Absorption(1/cm) Reflectivity", end="\n", file=f1)
                for n in range(dielec_y.shape[0]):
                    print(dielec_y[n][0], dielec_y[n][3], dielec_y[n][4], opt_n_lfc_y[n], opt_k_lfc_y[n], opt_abs_lfc_y[n], opt_ref_lfc_y[n], end="\n", file=f1)
                print (end="\n", file=f1)

            # Saving LFC other optical spectrum for z-direction
            with paropen(struct+'-5-Result-Optical-RPA-LFC-AllData_zdirection.dat', 'w') as f1:
                print("Energy(eV) Eps_real Eps_img Refractive_Index Extinction_Index Absorption(1/cm) Reflectivity", end="\n", file=f1)
                for n in range(dielec_z.shape[0]):
                    print(dielec_z[n][0], dielec_z[n][3], dielec_z[n][4], opt_n_lfc_z[n], opt_k_lfc_z[n], opt_abs_lfc_z[n], opt_ref_lfc_z[n], end="\n", file=f1)
                print (end="\n", file=f1)
            
        else:
            parprint('\033[91mERROR:\033[0mUnknown optical calculation type.')
            quit()

    elif Mode == 'LCAO':
        parprint('\033[91mERROR:\033[0mNot implemented in LCAO mode yet.')
    else:
        parprint('\033[91mERROR:\033[0mNot implemented in FD mode yet.')
    # Finish Optical calc
    time52 = time.time()

# -------------------------------------------------------------
# Step 6 - TIME
# -------------------------------------------------------------

with paropen(struct+'-6-Result-Log-Timings.txt', 'a') as f1:
    print("gpawsolve.py execution timings (seconds):", end="\n", file=f1)
    print("Execution started:", time0, end="\n", file=f1)
    if 'time11' in globals():
        print('Ground state: ', round((time12-time11),2), end="\n", file=f1)
    if Elastic_calc == True:
        print('Elastic calculation: ', round((time152-time151),2), end="\n", file=f1)
    if DOS_calc == True:
        print('DOS calculation: ', round((time22-time21),2), end="\n", file=f1)
    if Band_calc == True:
        print('Band calculation: ', round((time32-time31),2), end="\n", file=f1)
    if Density_calc == True:
        print('Density calculation: ', round((time42-time41),2), end="\n", file=f1)
    if Optical_calc == True:
        print('Optical calculation: ', round((time52-time51),2), end="\n", file=f1)
    print("---------------------------------", end="\n", file=f1)

# -------------------------------------------------------------
# Step Last - DRAWING BAND STRUCTURE AND DOS
# -------------------------------------------------------------
if drawfigs == True:
    # Draw graphs only on master node
    if world.rank == 0:
        # Elastic
        if Elastic_calc == True:
            eos.plot(struct+'-1.5-Graph-Elastic-EOS.png', show=True)
        # DOS
        if DOS_calc == True:
            if Spin_calc == True:
                ax = plt.gca()
                ax.plot(energies, -1.0*weights, 'y')
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
        if Band_calc == True:
            # Band Structure
            if Mode == 'PW-GW':
                f = plt.figure()
                plt.plot(xdata, banddata, '-b', '-r', linewidth=1)
                plt.xticks(X, GWkpoints, fontsize=8)
                plt.ylabel('Energy with respect to vacuum (eV)', fontsize=14)
                plt.tight_layout()
                plt.savefig(struct+'-3-Graph-Band.png')
                plt.show()
            else:
                bs.plot(filename=struct+'-3-Graph-Band.png', show=True, emax=energy_max)
else:
    # Draw graphs only on master node
    if world.rank == 0:
        # Elastic
        if Elastic_calc == True:
            eos.plot(struct+'-1.5-Result-Elastic-EOS.png')
        # DOS
        if DOS_calc == True:
            if Spin_calc == True:
                ax = plt.gca()
                ax.plot(energies, -1.0*weights, 'y')
                ax.plot(energies, weightsup, 'b')
                ax.set_xlabel('Energy [eV]')
                ax.set_ylabel('DOS [1/eV]')
            else:
                ax = plt.gca()
                ax.plot(energies, weights, 'b')
                ax.set_xlabel('Energy [eV]')
                ax.set_ylabel('DOS [1/eV]')
            plt.savefig(struct+'-2-Graph-DOS.png')
        if Band_calc == True:
            # Band Structure
            if Mode == 'PW-GW':
                f = plt.figure()
                plt.plot(xdata, banddata, '-b', '-r', linewidth=1)
                plt.xticks(X, GWkpoints, fontsize=8)
                plt.ylabel('Energy with respect to vacuum (eV)', fontsize=14)
                plt.tight_layout()
                plt.savefig(struct+'-3-Graph-Band.png')
                #plt.show()
            else:
                bs.plot(filename=struct+'-3-Graph-Band.png', show=False, emax=energy_max)
