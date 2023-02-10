'''
optimize_kptsdensity.py: Sample K-point density optimization with GPAW

Usage: $ gpaw -P<core_number> python optimize_kptsdensity.py <file.cif>
'''

from ase import *
from ase.io import read
from ase.geometry.cell import cell_to_cellpar
from gpaw import GPAW, PW
import sys
from ase.parallel import paropen, parprint
import numpy as np

# -------------------------------------------------------------
# ENTER PARAMETERS
# -------------------------------------------------------------
cutoffenergy = 700

kptsdensity_min = 1.5
kptsdensity_max = 3.0
kptsdensity_step = 0.1

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
# Start trying all k-density values
with paropen('OptimizeKPoints_Table-KPointDensity.txt', 'a') as f:
    f.write('K-density  Total_Energy\n')
    for k in np.arange(kptsdensity_min, kptsdensity_max, kptsdensity_step):
        bulk_configuration.calc = GPAW(mode=PW(cutoffenergy),
                                    xc=xc_used,
                                    kpts={'density': k, 'gamma': True},
                                    parallel={'band': 1},
                                    basis='dzp',
                                    txt='OptimizeKPointDensity_KDensity-%02d.txt' % k)
        parprint("K-density:"+str(k)+"   Potential_Energy:"+str(bulk_configuration.get_potential_energy()))
        f.write(str(k)+'  '+str(bulk_configuration.get_potential_energy())+'\n')
