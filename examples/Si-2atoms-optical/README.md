# Example: 2 Atoms Silicon Calculations (2 Step)

This example has two steps. First step is the ground state, DOS and band structure calculations. And the second step is the calculation of optical properties. Please do not forget to run optical calculation seperately.

To run the first step of calculation with MPI please execute the following command in this folder.

    mpirun -np 5 gpawsolve.py -o -c Si-Step1-ground_dos_band.py -i Si_mp-149_primitive_Example.cif
	
And then for the second step, execute the second command as

    gpawsolve.py -o -c Si-Step2-optical.py -i Si_mp-149_primitive_Example.cif