# Example: Effect of charge in graphene with defect with LCAO

This example has two input files. For calculating the neutral defected graphene with MPI 4 cores:

    mpirun -np 4 gpawsolve.py -i graphene-neutral.py -g graphene4x4withdefect.cif
    
or

    gpaw -P4 python ~/path-to-gpawtools/gpawsolve.py -- -i graphene-neutral.py -g graphene4x4withdefect.cif

Here, ~/path-to-gpawtools shows a full path your gpaw-tools folder. Results will be saved to "Neutral" folder.
	
And then for the calculation of charged defected graphene, execute the second command as

    mpirun -np 4 gpawsolve.py -i graphene-charged.py -g graphene4x4withdefect.cif
    
or

    gpaw -P4 python ~/path-to-gpawtools/gpawsolve.py -- -i graphene-charged.py -g graphene4x4withdefect.cif
    
Results will be saved to "Charged" folder.
