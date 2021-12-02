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
 | Method | Structure optim.    | Different XCs | Spin polarized | DOS | DFT+U | Band | Electron Density | Optical |
 | ------ | ------------------- | ------------- | -------------- | --- | ----- | ---- | ---------------- | ------- |
 |   PW   | Yes (No for GLLBSC) | Yes           | Yes            | Yes | Yes   | Yes  | Yes              | Yes     |
 | PW-G0W0| Yes                 | Yes           | No             | No  | No    | Yes  | No               | No      |
 | PW-EXX*| Yes (with PBE)      | No            | No             | No  | No    | No   | No               | No      |
 |  LCAO  | No                  | No            | No             | Yes | Yes   | Yes  | Yes              | No      |
 *: Just some ground state energy calculations for PBE0 and HSE06.
'''

import getopt, sys, os
import textwrap
import requests
import pickle
from argparse import ArgumentParser, HelpFormatter
from ase import *
from ase.parallel import paropen, world, parprint, broadcast
from gpaw import GPAW, PW, FermiDirac
from ase.optimize.lbfgs import LBFGS
from ase.io import read, write
from ase.units import Bohr
import matplotlib.pyplot as plt
from ase.dft.dos import DOS
from ase.constraints import UnitCellFilter
from ase.io.cif import write_cif
from pathlib import Path
from gpaw.response.df import DielectricFunction
from gpaw.response.g0w0 import G0W0
from gpaw.response.gw_bands import GWBands
from gpaw.xc.exx import EXX
import numpy as np
from numpy import genfromtxt

# DEFAULT VALUES
# These values (with bulk configuration) can be used to run this script without using inputfile (py file)
# and configuration file (cif file). 
# -------------------------------------------------------------
Mode = 'PW'             # Use PW, PW-GW, PW-EXX, LCAO, FD  (PW is more accurate, LCAO is quicker mostly.)
# -------------------------------------------------------------
DOS_calc = False         # DOS calculation
Band_calc = False        # Band structure calculation
Density_calc = False    # Calculate the all-electron density?
Optical_calc = False     # Calculate the optical properties

# -------------------------------------------------------------
# Parameters
# -------------------------------------------------------------
# ELECTRONIC
fmaxval = 0.05 			#
cut_off_energy = 340 	# eV
#kpts_density = 2.5     # pts per Å^-1  If the user prefers to use this, kpts_x,y,z will not be used automatically.
kpts_x = 5 			    # kpoints in x direction
kpts_y = 5				# kpoints in y direction
kpts_z = 5				# kpoints in z direction
Gamma = True
band_path = 'LGL'	    # Brillouin zone high symmetry points
band_npoints = 60		# Number of points between high symmetry points
energy_max = 15 		# eV. It is the maximum energy value for band structure figure.
Hubbard = {}            # Can be used like {'N': ':p,6.0'}, for none use {}
#Exchange-Correlation, choose one:
XC_calc = 'LDA'
#XC_calc = 'PBE'
#XC_calc = 'GLLBSC'
#XC_calc = 'revPBE'
#XC_calc = 'RPBE'
#XC_calc = 'B3LYP'
#Choose one for PW-EXX (Ground state calculations will be done with PBE):
#XC_calc = 'PBE0'
#XC_calc = 'HSE06'
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
num_of_bands = 16		#
optFDsmear = 0.05       # Fermi Dirac smearing for optical calculations
opteta=0.2             # Eta for Optical calculations
optdomega0=0.1         # Domega0 for Optical calculations
optomega2=10.0          # Frequency at which the non-lin freq grid has doubled the spacing
optecut=10             # Cut-off energy for optical calculations
optnblocks=4            # Split matrices in nblocks blocks and distribute them G-vectors
                        # or frequencies over processes

#GENERAL
# Which components of strain will be relaxed
# EpsX, EpsY, EpsZ, ShearYZ, ShearXZ, ShearXY
# Example: For a x-y 2D nanosheet only first 2 component will be true
whichstrain=[False, False, False, False, False, False]
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
__version__ = "v21.12.1b1"

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
        response = requests.get("https://api.github.com/repos/lrgresearch/gpaw-tools/releases/latest")
        parprint('gpawtools: This is the version '+str(__version__))
        parprint('The latest stable release was '+response.json()["tag_name"]+', which is published at '+response.json()["published_at"])
        parprint('Download the latest stable release at: '+response.json()["tarball_url"])
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
        if XC_calc in ['HSE06', 'PBE0']:
            parprint('Error: '+XC_calc+' can be used only in PW-EXX mode...')
            quit()
        if Spin_calc == True:
           numm = [Magmom_per_atom]*bulk_configuration.get_global_number_of_atoms()
           bulk_configuration.set_initial_magnetic_moments(numm)
        
        if restart == False:
            # PW Ground State Calculations
            parprint("Starting PW ground state calculation...")
            if 'kpts_density' in globals():
                calc = GPAW(mode=PW(cut_off_energy), xc=XC_calc, nbands='200%', setups= Hubbard, parallel={'domain': world.size}, 
                        spinpol=Spin_calc, kpts={'density': kpts_density, 'gamma': Gamma}, txt=struct+'-1-Log-Ground.txt')
            else:
                calc = GPAW(mode=PW(cut_off_energy), xc=XC_calc, nbands='200%', setups= Hubbard, parallel={'domain': world.size}, 
                        spinpol=Spin_calc, kpts={'size': (kpts_x, kpts_y, kpts_z), 'gamma': Gamma}, txt=struct+'-1-Log-Ground.txt')
            bulk_configuration.calc = calc
            if True in whichstrain:
                if XC_calc == 'GLLBSC':
                    parprint("Structure optimization LBFGS can not be used with GLLBSC xc.")
                    parprint("Do structure optimization with PBE, then use its final CIF as input.")
                    parprint("Quiting...")
                    quit()
                uf = UnitCellFilter(bulk_configuration, mask=whichstrain)
                relax = LBFGS(uf, trajectory=struct+'-1-Result-Ground.traj')
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

    elif Mode == 'PW-EXX':
        if restart == False:
            # PW Ground State Calculations
            parprint("Starting PW ground state calculation with PBE...")
            if 'kpts_density' in globals():
                calc = GPAW(mode=PW(cut_off_energy), xc='PBE', parallel={'domain': world.size}, kpts={'density': kpts_density, 'gamma': Gamma}, txt=struct+'-1-Log-Ground.txt')
            else:
                calc = GPAW(mode=PW(cut_off_energy), xc='PBE', parallel={'domain': world.size}, kpts={'size': (kpts_x, kpts_y, kpts_z), 'gamma': Gamma}, txt=struct+'-1-Log-Ground.txt')
            bulk_configuration.calc = calc
            uf = UnitCellFilter(bulk_configuration, mask=whichstrain)
            relax = LBFGS(uf, trajectory=struct+'-1-Result-Ground.traj')
            relax.run(fmax=fmaxval)  # Consider tighter fmax!
            calc.write(struct+'-1-Result-Ground.gpw', mode="all")
            # Writes final configuration as CIF file
            write_cif(struct+'-Final.cif', bulk_configuration)
        else:
            parprint("Passing PW ground state calculation...")

        if XC_calc in ['HSE06', 'PBE0']:
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
            parprint('ATTENTION:EXX mode results are only listed at: '+struct+'-1-Result-Ground-EXX_mode.txt')
            parprint('          Other files (DOS, band, etc...) are the results calculated with PBE.')

    elif Mode == 'PW-GW':
        if restart == False:
            # PW Ground State Calculations
            parprint("Starting PW only ground state calculation for GW calculation...")
            if 'kpts_density' in globals():
                calc = GPAW(mode=PW(cut_off_energy), xc=XC_calc, parallel={'domain': 1}, occupations=FermiDirac(0.001), kpts={'density': kpts_density, 'gamma': Gamma}, txt=struct+'-1-Log-Ground.txt')
            else:
                calc = GPAW(mode=PW(cut_off_energy), xc=XC_calc, parallel={'domain': 1}, occupations=FermiDirac(0.001), kpts={'size':(kpts_x, kpts_y, kpts_z), 'gamma': Gamma}, txt=struct+'-1-Log-Ground.txt')
            bulk_configuration.calc = calc
            uf = UnitCellFilter(bulk_configuration, mask=whichstrain)
            relax = LBFGS(uf, trajectory=struct+'-1-Result-Ground.traj')
            relax.run(fmax=fmaxval)  # Consider tighter fmax!
            bulk_configuration.get_potential_energy()
            calc.diagonalize_full_hamiltonian()
            calc.write(struct+'-1-Result-Ground.gpw', mode="all")
            # Writes final configuration as CIF file
            write_cif(struct+'-Final.cif', bulk_configuration)
        else:
            parprint("Passing ground state calculation for GW calculation...")

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
        if restart == False:
            parprint("Starting LCAO ground state calculation...")
            if 'kpts_density' in globals():
                calc = GPAW(mode='lcao', basis='dzp', setups= Hubbard, kpts={'density': kpts_density, 'gamma': Gamma}, parallel={'domain': world.size})
            else:
                calc = GPAW(mode='lcao', basis='dzp', setups= Hubbard, kpts={'size':(kpts_x, kpts_y, kpts_z), 'gamma': Gamma}, parallel={'domain': world.size})
            bulk_configuration.calc = calc
            relax = LBFGS(bulk_configuration, trajectory=struct+'-1-Result-Ground.traj')
            relax.run(fmax=fmaxval)  # Consider much tighter fmax!
            bulk_configuration.get_potential_energy()
            calc.write(struct+'-1-Result-Ground.gpw', mode='all')
            # Writes final configuration as CIF file
            write_cif(struct+'-Final.cif', bulk_configuration)
        else:
            parprint("Passing LCAO ground state calculation...")

    elif Mode == 'FD':
        parprint("FD mode is not implemented in gpaw-tools yet...")
        quit()
    else:
        parprint("Please enter correct mode information.")
        quit()

    # -------------------------------------------------------------
    # Step 2 - DOS CALCULATION
    # -------------------------------------------------------------
    if DOS_calc == True:
        parprint("Starting DOS calculation...")
        calc = GPAW(struct+'-1-Result-Ground.gpw', fixdensity=True, txt=struct+'-2-Log-DOS.txt')
        #energies, weights = calc.get_dos(npts=800, width=0)
        dos = DOS(calc, npts=501, width=0.3)
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
                    print("Atom no: "+str(j)+", Atom Symbol: "+chem_sym[j]+" --------------------", file=fd)
                    en, pdossd = calc.get_orbital_ldos(a=j, spin=0, angular='s', npts=501, width=0.3)
                    en, pdospd = calc.get_orbital_ldos(a=j, spin=0, angular='p', npts=501, width=0.3)
                    en, pdosdd = calc.get_orbital_ldos(a=j, spin=0, angular='d', npts=501, width=0.3)
                    en, pdosfd = calc.get_orbital_ldos(a=j, spin=0, angular='f', npts=501, width=0.3)
                    for x in zip(en-ef, pdossd, pdospd, pdosdd, pdosfd):
                        print(*x, sep=", ", file=fd)
                print("---------------------------------------------------- --------------------", file=fd)
            #Spin up
            with paropen(struct+'-2-Result-PDOS-Up.csv', "w") as fd:
                print("Energy, s-orbital, p-orbital, d-orbital, f-orbital", file=fd)
                for j in range(0, bulk_configuration.get_global_number_of_atoms()):
                    print("Atom no: "+str(j)+", Atom Symbol: "+chem_sym[j]+" --------------------", file=fd)
                    en, pdossu = calc.get_orbital_ldos(a=j, spin=1, angular='s', npts=501, width=0.3)
                    en, pdospu = calc.get_orbital_ldos(a=j, spin=1, angular='p', npts=501, width=0.3)
                    en, pdosdu = calc.get_orbital_ldos(a=j, spin=1, angular='d', npts=501, width=0.3)
                    en, pdosfu = calc.get_orbital_ldos(a=j, spin=1, angular='f', npts=501, width=0.3)
                    for x in zip(en-ef, pdossu, pdospu, pdosdu, pdosfu):
                        print(*x, sep=", ", file=fd)
                print("---------------------------------------------------- --------------------", file=fd)
        else:
            with paropen(struct+'-2-Result-PDOS.csv', "w") as fd:
                print("Energy, s-orbital, p-orbital, d-orbital, f-orbital", file=fd)
                for j in range(0, bulk_configuration.get_global_number_of_atoms()):
                    print("Atom no: "+str(j)+", Atom Symbol: "+chem_sym[j]+" --------------------", file=fd)
                    en, pdoss = calc.get_orbital_ldos(a=j, spin=0, angular='s', npts=501, width=0.3)
                    en, pdosp = calc.get_orbital_ldos(a=j, spin=0, angular='p', npts=501, width=0.3)
                    en, pdosd = calc.get_orbital_ldos(a=j, spin=0, angular='d', npts=501, width=0.3)
                    en, pdosf = calc.get_orbital_ldos(a=j, spin=0, angular='f', npts=501, width=0.3)
                    for x in zip(en-ef, pdoss, pdosp, pdosd, pdosf):
                        print(*x, sep=", ", file=fd)
                print("---------------------------------------------------- --------------------", file=fd)
    # -------------------------------------------------------------
    # Step 3 - BAND STRUCTURE CALCULATION
    # -------------------------------------------------------------
    if Band_calc == True:
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

    # -------------------------------------------------------------
    # Step 4 - ALL-ELECTRON DENSITY
    # -------------------------------------------------------------
    if Density_calc == True:
        parprint("Starting All-electron density calculation...")
        calc = GPAW(struct+'-1-Result-Ground.gpw', txt=struct+'-4-Log-ElectronDensity.txt')
        bulk_configuration.calc = calc
        np = calc.get_pseudo_density()
        n = calc.get_all_electron_density(gridrefinement=gridref)
        
        # Writing pseudo and all electron densities to cube file with Bohr unit
        write(struct+'-4-Result-All-electron_nall.cube', bulk_configuration, data=n * Bohr**3)
        write(struct+'-4-Result-All-electron_npseudo.cube', bulk_configuration, data=np * Bohr**3)

# -------------------------------------------------------------
# Step 5 - OPTICAL CALCULATION
# -------------------------------------------------------------
if Optical_calc == True:
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
            parprint('ERROR: Optical computations must be done separately. Please do ground calculations first.')
            quit()
    
        calc.get_potential_energy()

        calc.diagonalize_full_hamiltonian(nbands=num_of_bands)  # diagonalize Hamiltonian
        calc.write(struct+'-5-Result-Optical.gpw', 'all')  # write wavefunctions

        # Getting dielectric function spectrum
        parprint("Starting dielectric function calculation...")
        #from mpi4py import MPI
        df = DielectricFunction(calc=struct+'-5-Result-Optical.gpw',
                                eta=opteta, nblocks=world.size,
                                omega2=optomega2,
                                domega0=optdomega0,
                                ecut=optecut)
        # Writing to files as: omega, nlfc.real, nlfc.imag, lfc.real, lfc.imag 
        # Here lfc is local field correction
        df.get_dielectric_function( direction='x', 
                                    filename=struct+'-5-Result-Optical_dielec_xdirection.csv')
        df.get_dielectric_function( direction='y',
                                    filename=struct+'-5-Result-Optical_dielec_ydirection.csv')
        df.get_dielectric_function( direction='z',
                                    filename=struct+'-5-Result-Optical_dielec_zdirection.csv')
        # Loading dielectric function spectrum to numpy
        dielec_x = genfromtxt(struct+'-5-Result-Optical_dielec_xdirection.csv', delimiter=',')
        dielec_y = genfromtxt(struct+'-5-Result-Optical_dielec_ydirection.csv', delimiter=',')
        dielec_z = genfromtxt(struct+'-5-Result-Optical_dielec_zdirection.csv', delimiter=',')
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
        with paropen(struct+'-5-Result-Optical-NLFC-AllData_xdirection.dat', 'w') as f1:
            print("Energy(eV) Eps_real Eps_img Refractive_Index Extinction_Index Absorption(1/cm) Reflectivity", end="\n", file=f1)
            for n in range(dielec_x.shape[0]):
                print(dielec_x[n][0], dielec_x[n][1], dielec_x[n][2], opt_n_nlfc_x[n], opt_k_nlfc_x[n], opt_abs_nlfc_x[n], opt_ref_nlfc_x[n], end="\n", file=f1)
            print (end="\n", file=f1)

        # Saving NLFC other optical spectrum for y-direction
        with paropen(struct+'-5-Result-Optical-NLFC-AllData_ydirection.dat', 'w') as f1:
            print("Energy(eV) Eps_real Eps_img Refractive_Index Extinction_Index Absorption(1/cm) Reflectivity", end="\n", file=f1)
            for n in range(dielec_y.shape[0]):
                print(dielec_y[n][0], dielec_y[n][1], dielec_y[n][2], opt_n_nlfc_y[n], opt_k_nlfc_y[n], opt_abs_nlfc_y[n], opt_ref_nlfc_y[n], end="\n", file=f1)
            print (end="\n", file=f1)

        # Saving NLFC other optical spectrum for z-direction
        with paropen(struct+'-5-Result-Optical-NLFC-AllData_zdirection.dat', 'w') as f1:
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
        with paropen(struct+'-5-Result-Optical-LFC-AllData_xdirection.dat', 'w') as f1:
            print("Energy(eV) Eps_real Eps_img Refractive_Index Extinction_Index Absorption(1/cm) Reflectivity", end="\n", file=f1)
            for n in range(dielec_x.shape[0]):
                print(dielec_x[n][0], dielec_x[n][3], dielec_x[n][4], opt_n_lfc_x[n], opt_k_lfc_x[n], opt_abs_lfc_x[n], opt_ref_lfc_x[n], end="\n", file=f1)
            print (end="\n", file=f1)

        # Saving LFC other optical spectrum for y-direction
        with paropen(struct+'-5-Result-Optical-LFC-AllData_ydirection.dat', 'w') as f1:
            print("Energy(eV) Eps_real Eps_img Refractive_Index Extinction_Index Absorption(1/cm) Reflectivity", end="\n", file=f1)
            for n in range(dielec_y.shape[0]):
                print(dielec_y[n][0], dielec_y[n][3], dielec_y[n][4], opt_n_lfc_y[n], opt_k_lfc_y[n], opt_abs_lfc_y[n], opt_ref_lfc_y[n], end="\n", file=f1)
            print (end="\n", file=f1)

        # Saving LFC other optical spectrum for z-direction
        with paropen(struct+'-5-Result-Optical-LFC-AllData_zdirection.dat', 'w') as f1:
            print("Energy(eV) Eps_real Eps_img Refractive_Index Extinction_Index Absorption(1/cm) Reflectivity", end="\n", file=f1)
            for n in range(dielec_z.shape[0]):
                print(dielec_z[n][0], dielec_z[n][3], dielec_z[n][4], opt_n_lfc_z[n], opt_k_lfc_z[n], opt_abs_lfc_z[n], opt_ref_lfc_z[n], end="\n", file=f1)
            print (end="\n", file=f1)

    elif Mode == 'LCAO':
        parprint('Not implemented in LCAO mode yet.')
    else:
        parprint('Not implemented in FD mode yet.')

# -------------------------------------------------------------
# Step Last - DRAWING BAND STRUCTURE AND DOS
# -------------------------------------------------------------
if drawfigs == True:
    # Draw graphs only on master node
    if world.rank == 0:
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
