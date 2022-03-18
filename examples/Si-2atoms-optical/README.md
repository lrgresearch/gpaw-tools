# Example: 2 Atoms Silicon Calculations (3 Steps)

This example has three steps. First step is the ground state, DOS and band structure calculations. And second and third steps are the calculations of optical properties. Please do not forget to run optical calculations seperately.

To run the first step of calculation with MPI please execute the following command in this folder.

    mpirun -np 4 gpawsolve.py -o -i Si-Step1-ground_dos_band.py -g Si_mp-149_primitive_Example.cif
	
And then for the second step, there are two possibilities. Real and imaginary parts of dielectric function are usually calculated with random phase approximation (RPA). With GPAW, we can go beyond the RPA using the Bethe-Salpeter equation (BSE).

In this example, the second step is the optical properties calculation with RPA method. It can be executed as:

    gpawsolve.py -o -i Si-Step2-optical-RPA.py -g Si_mp-149_primitive_Example.cif

and, the third step is the optical properties calculation with BSE method. It can be executed as:

    gpawsolve.py -o -i Si-Step2-optical-BSE.py -g Si_mp-149_primitive_Example.cif

If the users wants to see results at the end of the calculation, they can use "-d" argument for the first script. Because of optical calculations does not support drawing results at the end, users can not use "-d" argument.

Output files of optical calculations will be named accordingly. There will be no data loss.