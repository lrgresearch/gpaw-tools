'''
optimize_cutoff.py: Sample Cut-off Energy Optimization with GPAW

Usage: $ gpaw -P<core_number> python optimize_cutoff.py <file.cif>
If the atom positions are not given with an external file like <file.cif>
it must be provided as a ASE object below part "Bulk Configuration"
Related ASE object can be produced with ciftoase.py script.
'''
from ase import *
from ase.io import read
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
# DO NOT NEED TO CHANGE ANYTHING UNDER THIS POINT
# -------------------------------------------------------------
# if bulk structure is given with CIF, XYZ, etc... file
if len(sys.argv) > 1:
    inFile = sys.argv[1]
    bulk_configuration = read(inFile, index='-1')

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
