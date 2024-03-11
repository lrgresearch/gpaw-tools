# Example: Elastic and Electronic Properties of Rocksalt TiC

Ground, Elastic, DOS and Band calculations of Rocksalt TiC. PW with 600 eV cutoff, 7x7x7 kpoints as done in. To run the calculation with MPI on 4 cores please execute the following command in this folder.

    mpirun -np 4 gpawsolve.py -i TiC.py -g TiC_mp-631_primitive-Final.cif
    
or

    gpaw -P4 python ~/path-to-gpawtools/gpawsolve.py -- -i TiC.py -g TiC_mp-631_primitive-Final.cif

Here, ~/path-to-gpawtools shows a full path your gpaw-tools folder.
	
