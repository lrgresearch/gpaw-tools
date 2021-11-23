# Example: Pristine graphene and graphene with defect with LCAO

This example has input files. Because Outdirname variable is not used in config file, the name of the result folder is coming from the input file. For calculating the pristine graphene with MPI 4 cores:

    mpirun -np 4 gpawsolve.py -o -i graphene.py -g graphene4x4.cif
	
And then for the calculation of graphene with defect, firstly dod not forget to close graphs of the first calculation,then execute the second command as

    mpirun -np 4 gpawsolve.py -o -i graphene.py -g graphene4x4withdefect.cif