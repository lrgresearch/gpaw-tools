# Example: GW Aproximation calculation for MoS2

To run the calculation with MPI on 4 cores please execute the following as:

    mpirun -np 4 gpawsolve.py -o -c MoS2-GW.py -i MoS2-structure.cif
	
or if you want to run it on a single core

    gpawsolve.py -o -c MoS2-GW.py -i MoS2-structure.cif