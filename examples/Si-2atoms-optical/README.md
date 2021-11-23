# Example: 2 Atoms Silicon Calculations (2 Step)

This example has two steps. First step is the ground state, DOS and band structure calculations. And the second step is the calculation of optical properties. Please do not forget to run optical calculation seperately.

To run the first step of calculation with MPI please execute the following command in this folder.

    mpirun -np 5 gpawsolve.py -o -i Si-Step1-ground_dos_band.py -g Si_mp-149_primitive_Example.cif
	
And then for the second step, execute the second command as

    gpawsolve.py -o -i Si-Step2-optical.py -g Si_mp-149_primitive_Example.cif

If the users wants to see results at the end of the calculation, they can use "-d" argument for the first script. Because of optical calculations does not support drawing results at the end, users can not use "-d" argument.