# Example: Pristine graphene and graphene with defect with LCAO

This example has input files. For calculating the pristine graphene with MPI 4 cores:

    mpirun -np 4 gpawsolve.py -i graphene.py -g graphene4x4.cif
    
or

    gpaw -P4 python ~/path-to-gpawtools/gpawsolve.py -- -i graphene.py -g graphene4x4.cif

Here, ~/path-to-gpawtools shows a full path your gpaw-tools folder.
	
And then for the calculation of graphene with defect, firstly dod not forget to close graphs of the first calculation,then execute the second command as

    mpirun -np 4 gpawsolve.py -i graphene.py -g graphene4x4withdefect.cif
    
or

    gpaw -P4 python ~/path-to-gpawtools/gpawsolve.py -- -i graphene.py -g graphene4x4withdefect.cif
