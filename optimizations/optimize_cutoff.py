'''
optimize_cutoff.py: Sample Cut-off Energy Optimization with GPAW

Usage: $ gpaw -P<core_number> python optimize_cutoff.py <file.cif>
'''
from ase import *
from ase.io import read
from ase.geometry.cell import cell_to_cellpar
from gpaw import GPAW, PW
import sys
from ase.parallel import paropen, world, parprint

#-------------------------------------------------------------
# ENTER PARAMETERS
# -------------------------------------------------------------
cutoff_min = 200
cutoff_max = 400
cutoff_step = 50

kpts_x = 7
kpts_y = 7
kpts_z = 7

xc_used = 'GLLBSCM'

# -------------------------------------------------------------
# DO NOT NEED TO CHANGE ANYTHING UNDER THIS POINT
# -------------------------------------------------------------
# Read bulk structure from CIF
if len(sys.argv) > 1:
    inFile = sys.argv[1]
    bulk_configuration = read(inFile, index='-1')

a, b, c, alpha, beta, gamma = cell_to_cellpar(bulk_configuration.get_cell(), radians=False)
parprint("a:"+str(a)+" , b:"+str(b)+" , c:"+str(c)+" , Alpha:"+str(alpha)+" , Beta:"+str(beta)+" , Gamma:"+str(gamma))
# Start trying all cut-off energies
with paropen('OptimizeCutOff_Table-CutOff.txt', 'a') as f:
    f.write('CutOff_Energy  Total_Energy\n')
    for ecut in range(cutoff_min, cutoff_max+1, cutoff_step):
        bulk_configuration.calc = GPAW(mode=PW(ecut),
                                    xc=xc_used,
                                    kpts=(kpts_x, kpts_y, kpts_z),
                                    parallel={'domain': world.size},
                                    basis='dzp',
                                    txt='OptimizeCutOff_%d.txt' % ecut)
        parprint ("CutOff_energy:"+str(ecut)+"  Potential_Energy:"+str(bulk_configuration.get_potential_energy()))
        f.write(str(ecut)+'  '+str(bulk_configuration.get_potential_energy())+'\n')
