# Example: Pristine graphene and graphene with defect with LCAO

This example has input files. For calculating the pristine graphene with MPI 4 cores:

    mpirun -np 4 gpawsolve.py -i graphene.py -g graphene4x4.cif
	
And then for the calculation of graphene with defect, firstly dod not forget to close graphs of the first calculation,then execute the second command as

    mpirun -np 4 gpawsolve.py -i graphene.py -g graphene4x4withdefect.cif
