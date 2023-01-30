# Example: HSE calculations of Silicon

This example uses hybrid XC HSE06 for the calculations. It can do ground state, DOS and band structure calculations. Because HSE calculations are much more slower than standard PBE calculations (Sometimes few thousand times slower with conda), convergence values listed in the input file, kept low to finish the calculation quicker. 

Please use proper convergence values and always use HPC for your HSE calculations :)

You can run this example with:

    mpirun -np 4 gpawsolve.py -o -i Si-with-HSE.py -g Si_mp-149_primitive.cif
	
Normally, prior to HSE calculations, you can prefer to do PBE calculations with structure optimization. Then you can continue to use HSE.
