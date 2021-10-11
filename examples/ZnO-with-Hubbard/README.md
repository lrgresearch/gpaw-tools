# Example: Electronic Properties of Wurtzite ZnO with DFT+U

**Note:** You can not use this example with `gg.py`. It is only working with CIF files.

Ground, DOS and Band calculations of Wurtzite bulk ZnO. PW with 340 eV cutoff, 5x5x5 kpoints. Hubbard parameters were chosen as 7 eV and 10 eV for O-p and Zn-d orbitals, respectively. The important thing is that the positions are given with Bulk object. Therefore the line `from ase.build import bulk` must be included at the beginning of the configuration file.

To run the calculation with MPI on 4 cores please execute the following command in this folder.

    mpirun -np 4 gpawsolve.py -o -c ZnO_withHubbard.py

or calculation with drawing band and DOS at the end:

	mpirun -np 4 gpawsolve.py -o -d -c ZnO_withHubbard.py

and, if you want to use -o argument like above you must give the name of the directory inside the configuration file like

    Outdirname = 'ZnO-withHubbard-results'

There is also a without Hubbard configuration file to compare:

	mpirun -np 4 gpawsolve.py -o -c ZnO_woHubbard.py
