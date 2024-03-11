# Example: Phonon dispersion of Bulk Aluminum

Phonon dispersion calculation of Bulk Aluminum. Ground state calculations will be done with PW, 700 eV cutoff, 5x5x5 kpoints. To run the calculation with MPI on 4 cores please execute the following command in this folder.

    mpirun -np 4 gpawsolve.py -i Al-phonon.py -g Al_mp-134_primitive.cif
    
or

    gpaw -P4 python ~/path-to-gpawtools/gpawsolve.py -- -i Al-phonon.py -g Al_mp-134_primitive.cif

Here, ~/path-to-gpawtools shows a full path your gpaw-tools folder.
	
WARNING: Phonon calculations are done with Phonopy package plus GPAW and ASE. However, Phonopy integration to gpaw-tools is far from mature. Please use this function with care.
