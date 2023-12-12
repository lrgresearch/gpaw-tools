# Example: Electronic Properties of Bulk GaAs

**Note:** You can not use this example with `gg.py`. It is only working with CIF files.

Ground, DOS and Band calculations of Bulk GaAs. PW with 300 eV cutoff, 2.5 kpoint per Angstrom k-point density. The important thing is that the positions are given with Atom object. To run the calculation with MPI on 4 cores please execute the following command in this folder.

    mpirun -np 4 gpawsolve.py -i bulk_gaas.py

or

    gpaw -P4 python ~/path-to-gpawtools/gpawsolve.py -- -i bulk_gaas.py

Here, ~/path-to-gpawtools shows a full path your gpaw-tools folder.
	
When you use Atoms object inside configuration file, please note that you must add

    from ase import Atoms,Atom

and, if you want to use -o argument like above you must give the name of the directory inside the configuration file like

    Outdirname = 'bulk-gaas-results'
