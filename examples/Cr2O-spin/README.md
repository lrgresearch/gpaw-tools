# Example: Spin-dependent electronic properties of Cr2O

This example is for showing spin dependent calculations. The input files is obtained from (https://materialsproject.org/materials/mp-1206821/) .

To run with MPI for 4 cores please execute the following command.

    mpirun -np 4 gpawsolve.py -o -i Cr2O.py -g Cr2O_mp-1206821_primitive.cif
	
To visualize the electron-densities, you can use the result cube files `Cr2O_mp-1206821_primitive-4-Result-All-electron_n.cube` and `Cr2O_mp-1206821_primitive-4-Result-All-electron_np.cube` with Vesta.
