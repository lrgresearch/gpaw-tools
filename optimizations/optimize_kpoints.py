'''
optimize_kpoints.py: Sample K-point Optimization with GPAW

Usage: $ gpaw -P<core_number> python optimize_kpoints.py <file.cif>
'''

from ase import *
from ase.io import read
from ase.geometry.cell import cell_to_cellpar
import sys
from gpaw import GPAW, PW
from ase.parallel import paropen, parprint

# -------------------------------------------------------------
# ENTER PARAMETERS
# -------------------------------------------------------------
cutoffenergy = 200

kpoint_min = 4
kpoint_max = 11

xc_used = 'GLLBSCM'

# -------------------------------------------------------------
# DO NOT NEED TO CHANGE ANYTHING UNDER THIS POINT
# -------------------------------------------------------------
# Read bulk structure from CIF
if len(sys.argv) > 1:
    inFile = sys.argv[1]
    bulk_configuration = read(inFile, index='-1')

a, b, c, alpha, beta, gamma = cell_to_cellpar(bulk_configuration.get_cell(), radians=False)

# Start trying all k-points
with paropen('OptimizeKPoints_Table-KPoints.txt', 'a') as f:
    f.write('K-point  Total_Energy\n')
    for k in range(kpoint_min, kpoint_max):
        bulk_configuration.calc = GPAW(mode=PW(cutoffenergy),
                                    xc=xc_used,
                                    kpts=(k, k, k),
                                    parallel={'band': 1},
                                    basis='dzp',
                                    txt='OptimizeKPoints_KPoints-%02d.txt' % k)
        parprint("K-point:"+str(k)+"   Potential_Energy:"+str(bulk_configuration.get_potential_energy()))
        f.write(str(k)+'  '+str(bulk_configuration.get_potential_energy())+'\n')
