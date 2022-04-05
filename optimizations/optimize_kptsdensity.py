'''
optimize_kptsdensity.py: Sample K-point density optimization with GPAW

Usage: $ gpaw -P<core_number> python optimize_kptsdensity.py <file.cif>
If the atom positions are not given with an external file like <file.cif>
it must be provided as a ASE object below part "Bulk Configuration"
Related ASE object can be produced with ciftoase.py script.
'''

from ase import *
from gpaw import GPAW, PW
from ase.parallel import paropen, parprint
import numpy as np

# -------------------------------------------------------------
# ENTER PARAMETERS
# -------------------------------------------------------------
cutoffenergy = 700
kptsdensity_min = 1.5
kptsdensity_max = 3.0
kptsdensity_step = 0.1


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
with paropen('OptimizeKPoints_Table-KPointDensity.txt', 'a') as f:
    f.write('K-density  Total_Energy\n')
    for k in np.arange(kptsdensity_min, kptsdensity_max, kptsdensity_step):
        bulk_configuration.calc = GPAW(mode=PW(cutoffenergy),
                                    xc='GLLBSCM',
                                    kpts={'density': k, 'gamma': True},
                                    parallel={'band': 1},
                                    basis='dzp',
                                    txt='OptimizeKPointDensity_KDensity-%02d.txt' % k)
        parprint("K-density:"+str(k)+"   Potential_Energy:"+str(bulk_configuration.get_potential_energy()))
        f.write(str(k)+'  '+str(bulk_configuration.get_potential_energy())+'\n')
