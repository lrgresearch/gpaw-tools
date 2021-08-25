from ase import *
from ase.parallel import paropen, world, parprint
from gpaw import GPAW, PW, FermiDirac
from ase.optimize.lbfgs import LBFGS
from ase.io import read, write
import matplotlib.pyplot as plt
from ase.dft.dos import DOS
from ase.constraints import UnitCellFilter
from ase.io.cif import write_cif
from pathlib import Path
from gpaw.response.df import DielectricFunction
import numpy as np
import getopt, sys, os

# gpawsolve.py: Easy PW/LCAO Calculation Script for GPAW
# --------------------------------------------------------
HelpText = """ 
 Command line usage: gpawsolve.py -i <inputfile.cif> -c -o -h
 Argument list:
                   -i, --Input  : Use input CIF file
                   -c, --Config : Use configuration file in the main directory for parameters (config.py)
                   -o, --Outdir : Save everything to a output directory with naming /inputfile. If there is no input file given and Atoms object is used in
                      gpawsolve.py file then the directory name will be /gpawsolve. If you change gpawsolve.py name to anyname.py
                      then the directory name will be /anyname
                   -h --Help    : Help
 Usage: Change number with core numbers/threads to use. I am suggesting to use total number of cores(or threads) - 1
 Usage: $ gpaw -P8 python gpawsolve.py
 For AMD CPUs or using Intel CPUs without hyperthreading: (Example CPU is intel here, 4 cores or 8 threads)
        $ mpirun -n 4 gpaw python gpawsolve.py
 For using all threads provided by Intel Hyperthreading technology:
        $ mpirun --use-hwthread-cpus -n 8 gpaw python gpawsolve.py 
 -------------------------------------------------------------
 Calculation selector
 -------------------------------------------------------------

 | Method | Strain_minimization | Several XCs | Spin polarized | DOS | Band | Electron Density | Optical |
 | ------ | ------------------- | ----------- | -------------- | --- | ---- | ---------------- | ------- |
 |   PW   | Yes                 | Yes         | Yes            | Yes | Yes  | Yes              | Yes     |
 |  LCAO  | No                  | No          | No             | Yes | Yes  | Yes              | No     |
"""
# IF YOU WANT TO USE CONFIG FILE, PLEASE COPY/PASTE FROM HERE:>>>>>>>
# -------------------------------------------------------------
Use_PW = True          # Use PW or LCAO? (PW is more accurate, LCAO is quicker mostly.)
# -------------------------------------------------------------
DOS_calc = True         # DOS calculation
Band_calc = True        # Band structure calculation
Density_calc = True    # Calculate the all-electron density?
Optical_calc = False     # Calculate the optical properties

# -------------------------------------------------------------
# Parameters
# -------------------------------------------------------------
# ELECTRONIC
fmaxval = 0.05 			#
cut_off_energy = 340 	# eV
kpts_x = 5 			    # kpoints in x direction
kpts_y = 5				# kpoints in y direction
kpts_z = 1				# kpoints in z direction
band_path = 'GMKG'	    # Brillouin zone high symmetry points
band_npoints = 40		# Number of points between high symmetry points 
energy_max = 15 		# eV. It is the maximum energy value for band structure figure.
#Exchange-Correlation, choose one:
#XC_calc = 'LDA'
XC_calc = 'PBE'
#XC_calc = 'revPBE'
#XC_calc = 'RPBE'
Spin_calc = False        # Spin polarized calculation?
gridref = 4             # refine grid for all electron density (1, 2 [=default] and 4)

# OPTICAL
num_of_bands = 16		#
optFDsmear = 0.05       # Fermi Dirac smearing for optical calculations
opteta=0.05             # Eta for Optical calculations
optdomega0=0.02         # Domega0 for Optical calculations
optnblocks=4            # Split matrices in nblocks blocks and distribute them G-vectors or frequencies over processes

#GENERAL
draw_graphs = True		# Draw DOS and band structure on screen (yes for draw, small letters)

# Which components of strain will be relaxed
# EpsX, EpsY, EpsZ, ShearYZ, ShearXZ, ShearXY
# Example: For a x-y 2D nanosheet only first 2 component will be true
whichstrain=[False, False, False, False, False, False]

WantCIFexport = False
# <<<<<<< TO HERE TO FILE config.py IN SAME DIRECTORY AND USE -c FLAG WITH COMMAND
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
# Remove 1st argument from the
# list of command line arguments
argumentList = sys.argv[1:]
 
# Options
options = "ohci:"
 
# Long options
long_options = ["Outdir", "Help", "Config", " Input ="]
 
try:
    # Parsing argument
    arguments, values = getopt.getopt(argumentList, options, long_options)
     
    # checking each argument
    for currentArgument, currentValue in arguments:
 
        if currentArgument in ("-o", "--Outdir"):
            outdir = True
        
        elif currentArgument in ("-h", "--Help"):
            parprint (HelpText)
             
        elif currentArgument in ("-c", "--Config"):
            import config
            # There must be some elegant way to do this.
            Use_PW = config.Use_PW
            DOS_calc = config.DOS_calc
            Band_calc = config.Band_calc
            Density_calc = config.Density_calc
            Optical_calc = config.Optical_calc
            fmaxval = config.fmaxval
            cut_off_energy = config.cut_off_energy
            kpts_x = config.kpts_x
            kpts_y = config.kpts_y
            kpts_z = config.kpts_z
            band_path = config.band_path
            band_npoints = config.band_npoints
            energy_max = config.energy_max
            XC_calc = config.XC_calc
            Spin_calc = config.Spin_calc
            gridref = config.gridref
            num_of_bands = config.num_of_bands
            optFDsmear =config.optFDsmear
            opteta=config.opteta
            optdomega0=config.optdomega0
            optnblocks=config.optnblocks
            draw_graphs = config.draw_graphs
            whichstrain=config.whichstrain
            WantCIFexport = config.WantCIFexport
             
        elif currentArgument in ("-i", "--Input"):
            inFile = currentValue
    
except getopt.error as err:
    # output error, and return with an error code
    parprint (str(err))

# If there is a CIF input, use it. Otherwise use the bulk configuration provided above.
if inFile is None:
    struct = Path(__file__).stem # All files will get their names from this file
    parprint("Number of atoms provided in Atoms object:"+str(bulk_configuration.get_global_number_of_atoms()))
else:
    struct = Path(inFile).stem
    bulk_configuration = read(inFile, index='-1')
    parprint("Number of atoms imported from CIF file:"+str(bulk_configuration.get_global_number_of_atoms()))
    
# Control if outdir is set or not
if outdir is None:
    #No change is necessary
    parprint("Output directory is the main directory")
else:
    if not os.path.isdir(struct):
        os.makedirs(struct, exist_ok=True)
    struct = os.path.join(struct,struct)

# -------------------------------------------------------------
# Step 1 - GROUND STATE
# -------------------------------------------------------------
if Use_PW == True:
    parprint("Starting PW ground state calculation...")
    calc = GPAW(mode=PW(cut_off_energy), xc=XC_calc, parallel={'domain': world.size}, spinpol=Spin_calc, kpts=[kpts_x, kpts_y, kpts_z], txt=struct+'-1-Log-Ground.txt')
    bulk_configuration.calc = calc

    uf = UnitCellFilter(bulk_configuration, mask=whichstrain)
    relax = LBFGS(uf, trajectory=struct+'-1-Result-Ground.traj')
    relax.run(fmax=fmaxval)  # Consider tighter fmax!

    bulk_configuration.get_potential_energy()
    if Density_calc == True:
        #This line makes huge GPW files. Therefore it is better to use this if else
        calc.write(struct+'-1-Result-Ground.gpw', mode="all")
    else:
        calc.write(struct+'-1-Result-Ground.gpw')
else:
    parprint("Starting LCAO ground state calculation...")
    calc = GPAW(mode='lcao', basis='dzp', kpts=(kpts_x, kpts_y, kpts_z), parallel={'domain': world.size})
    bulk_configuration.calc = calc

    relax = LBFGS(bulk_configuration, trajectory=struct+'-1-Result-Ground.traj')
    relax.run(fmax=fmaxval)  # Consider much tighter fmax!

    bulk_configuration.get_potential_energy()
    calc.write(struct+'-1-Result-Ground.gpw', mode='all')

    
if WantCIFexport == True:
    write_cif(struct+'-Final.cif', bulk_configuration)
    
# -------------------------------------------------------------
# Step 2 - DOS CALCULATION
# -------------------------------------------------------------
if DOS_calc== True:
    parprint("Starting DOS calculation...")
    calc = GPAW(struct+'-1-Result-Ground.gpw', txt=struct+'-2-Log-DOS.txt')
    #energies, weights = calc.get_dos(npts=800, width=0)
    dos = DOS(calc, npts=500, width=0)
    if Spin_calc == True:
        energies = dos.get_energies()
        weights = dos.get_dos(spin=0)
        weightsup = dos.get_dos(spin=1)
    else:
        energies = dos.get_energies()
        weights = dos.get_dos()

    fd = open(struct+'-2-Result-DOS.txt', "w")
    if Spin_calc == True:
        for x in zip(energies, weights, weightsup):
            print(*x, sep=", ", file=fd)
    else:
        for x in zip(energies, weights):
            print(*x, sep=", ", file=fd)
    fd.close()

# -------------------------------------------------------------
# Step 3 - BAND STRUCTURE CALCULATION
# -------------------------------------------------------------
if Band_calc == True:
    parprint("Starting band structure calculation...")
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
        f1 = open(struct+'-3-Result-Band-Down.dat', 'w')
        for n1 in range(num_of_bands):
            for k1 in range(band_npoints):
                print(k1, eps_skn[0, k1, n1], end="\n", file=f1)
            print (end="\n", file=f1)
        f1.close()
        f2 = open(struct+'-3-Result-Band-Up.dat', 'w')
        for n2 in range(num_of_bands):
            for k2 in range(band_npoints):
                print(k2, eps_skn[1, k2, n2], end="\n", file=f2)
            print (end="\n", file=f2)
        f2.close()
    else:
        eps_skn = np.array([[calc.get_eigenvalues(k,s)
                            for k in range(band_npoints)]
                            for s in range(1)]) - ef
        f = open(struct+'-3-Result-Band.dat', 'w')
        for n in range(num_of_bands):
            for k in range(band_npoints):
                print(k, eps_skn[0, k, n], end="\n", file=f)
            print (end="\n", file=f)
        f.close()

# -------------------------------------------------------------
# Step 4 - ALL-ELECTRON DENSITY
# -------------------------------------------------------------
if Density_calc == True:
    parprint("Starting All-electron density calculation...")
    calc = GPAW(struct+'-1-Result-Ground.gpw', txt=struct+'-4-Log-ElectronDensity.txt')
    bulk_configuration.calc = calc
    np = calc.get_pseudo_density()
    n = calc.get_all_electron_density(gridrefinement=gridref)

    write(struct+'-4-Result-All-electron_n.cube', bulk_configuration, data=n)
    write(struct+'-4-Result-All-electron_np.cube', bulk_configuration, data=np)


# -------------------------------------------------------------
# Step 5 - OPTICAL CALCULATION
# -------------------------------------------------------------
if Optical_calc == True:
    if Use_PW == True:
        parprint("Starting optical calculation...")
        calc = GPAW(struct+'-1-Result-Ground.gpw',
                parallel={'domain': 1},
                txt=struct+'-5-Log-Optical.txt',
                nbands=num_of_bands,
                fixdensity=True,
                symmetry='off',
                occupations=FermiDirac(optFDsmear))

        calc.get_potential_energy()

        calc.diagonalize_full_hamiltonian(nbands=num_of_bands)  # diagonalize Hamiltonian
        calc.write(struct+'-5-Result-Optical.gpw', 'all')  # write wavefunctions

        # Getting absorption spectrum
        parprint("Starting dielectric function calculation...")
        df = DielectricFunction(calc=struct+'-5-Result-Optical.gpw',
                                eta=opteta,
                                nblocks=world.size,
                                domega0=optdomega0,
                                ecut=cut_off_energy)
        df.get_dielectric_function( direction='x', filename=struct+'-5-Result-Optical_abs_xdirection.csv')
        df.get_dielectric_function( direction='y', filename=struct+'-5-Result-Optical_abs_ydirection.csv')
        df.get_dielectric_function( direction='z', filename=struct+'-5-Result-Optical_abs_zdirection.csv')
    else:
        parprint('Not implemented in LCAO mode yet.')

# -------------------------------------------------------------
# Step Last - DRAWING BAND STRUCTURE AND DOS
# -------------------------------------------------------------
if draw_graphs == True:
    # Draw graphs only on master node
    if world.rank == 0:
        # DOS
        if DOS_calc == True:
            if Spin_calc == True:
                ax = plt.gca()
                ax.plot(energies, -1.0*weights, 'r')
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
            bs.plot(filename=struct+'-3-Graph-Band.png', show=True, emax=energy_max)
