# Example: Electronic Properties of Wurtzite ZnO with DFT+U

**Note:** You can not use this example with `gg.py`. It is only working with CIF files.

Ground, DOS and Band calculations of Wurtzite bulk ZnO. PW with 340 eV cutoff, 5x5x5 kpoints. Hubbard parameters were chosen as 7 eV and 10 eV for O-p and Zn-d orbitals, respectively. Hubbard parameter can be changed with changing Setup_params variable as {'Element': ':orbital,value,0'} if the user wants the Hubbard correction as normalized. In the example,  {'Element': ':orbital,value'} type not normalized input is used. The important thing is that the positions are given with Bulk object. Therefore the line `from ase.build import bulk` must be included at the beginning of the configuration file.

To run the calculation with MPI on 4 cores please execute the following command in this folder.

    mpirun -np 4 gpawsolve.py -i ZnO_withHubbard.py

or calculation with drawing band and DOS at the end:

	mpirun -np 4 gpawsolve.py -d -i ZnO_withHubbard.py

and,

    Outdirname = 'ZnO-withHubbard-results'

There is also a without Hubbard configuration file to compare:

	mpirun -np 4 gpawsolve.py -i ZnO_woHubbard.py
