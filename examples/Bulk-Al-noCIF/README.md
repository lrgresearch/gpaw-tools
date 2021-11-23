# Example: Electronic Properties of Bulk Aluminum

**Note:** You can not use this example with `gg.py`. It is only working with CIF files.

Ground, DOS and Band calculations of Bulk Aluminum. PW with 340 eV cutoff, 4x4x4 kpoints. The important thing is that the ositions are given with Atom object. To run the calculation with MPI on 4 cores please execute the following command in this folder.

    mpirun -np 4 gpawsolve.py -o -i bulk_aluminum.py
	
When you use Atoms object inside configuration file, please note that you must add

    from ase import Atoms

and, if you want to use -o argument like above you must give the name of the directory inside the configuration file like

    Outdirname = 'bulk-aluminum-results'
