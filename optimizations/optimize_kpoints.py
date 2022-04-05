'''
optimize_kpoints.py: Sample K-point Optimization with GPAW

Usage: $ gpaw -P<core_number> python optimize_kpoints.py <file.cif>
If the atom positions are not given with an external file like <file.cif>
it must be provided as a ASE object below part "Bulk Configuration"
Related ASE object can be produced with ciftoase.py script.
'''

from ase import *
from ase.io import read
import sys
from gpaw import GPAW, PW
from ase.parallel import paropen, parprint

# -------------------------------------------------------------
# ENTER PARAMETERS
# -------------------------------------------------------------
cutoffenergy = 200
kpoint_min = 4
kpoint_max = 11


# -------------------------------------------------------------
# Bulk Configuration - YOU DON'T NEED IF YOU USE CIF FILE
# -------------------------------------------------------------
bulk_configuration = Atoms(['C' for i in range(2)],
              [(   3.790492,    2.188441,    1.547462),
               (   2.526995,    1.458961,    1.031641)],
              pbc = (True,True,True))
bulk_configuration.set_cell([[    2.526995,     0.000000,     0.000000],
                [    1.263497,     2.188441,     0.000000],
                [    1.263497,     0.729480,     2.063282]],
                scale_atoms = False)
# -------------------------------------------------------------
# CONVERGE KPOINTS
# -------------------------------------------------------------
# if bulk structure is given with CIF, XYZ, etc... file
if len(sys.argv) > 1:
    inFile = sys.argv[1]
    bulk_configuration = read(inFile, index='-1')

cell0 = bulk_configuration.cell
with paropen('OptimizeKPoints_Table-KPoints.txt', 'a') as f:
    f.write('K-point  Total_Energy\n')
    for k in range(kpoint_min, kpoint_max):
        bulk_configuration.calc = GPAW(mode=PW(cutoffenergy),
                                    xc='GLLBSCM',
                                    kpts=(k, k, k),
                                    parallel={'band': 1},
                                    basis='dzp',
                                    txt='OptimizeKPoints_KPoints-%02d.txt' % k)
        parprint("K-point:"+str(k)+"   Potential_Energy:"+str(bulk_configuration.get_potential_energy()))
        f.write(str(k)+'  '+str(bulk_configuration.get_potential_energy())+'\n')
