'''
optimize_kpoints.py: Sample K-point Optimization with GPAW

Usage: $ gpaw -P<core_number> python optimize_kpoints.py <file.cif>
If the atom positions are not given with an external file like <file.cif>
it must be provided as a ASE object below part "Bulk Configuration"
Related ASE object can be produced with ciftoase.py script.
'''

from ase import *
import sys, os
from gpaw import GPAW, PW
from ase.parallel import paropen, parprint

# -------------------------------------------------------------
# ENTER PARAMETERS
# -------------------------------------------------------------
cutoffenergy = 200
kpoint_min = 4
kpoint_max = 11


# -------------------------------------------------------------
# Bulk Configuration
# -------------------------------------------------------------
bulk_configuration = Atoms(
    [
    Atom('Si', ( 0.0, 0.0, 0.0 )),
    Atom('Si', ( 1.35765, 1.35765, 1.35765 )),
    ],
    cell=[(0.0, 2.7153, 2.7153), (2.7153, 0.0, 2.7153), (2.7153, 2.7153, 0.0)],
    pbc=True,
    )
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
                                    xc='PBE',
                                    kpts=(k, k, k),
                                    parallel={'band': 1},
                                    basis='dzp',
                                    txt='OptimizeKPoints_KPoints-%02d.txt' % k)
        parprint("K-point:"+str(k)+"   Potential_Energy:"+str(bulk_configuration.get_potential_energy()))
        f.write(str(k)+'  '+str(bulk_configuration.get_potential_energy())+'\n')
