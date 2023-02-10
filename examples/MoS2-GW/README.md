# Example: GW Aproximation calculation for MoS2

To run the calculation, do not usi mpirun and use `gpawsolve.py` as

    gpawsolve.py -i MoS2-GW.py -g MoS2-structure.cif

After first run, you can by pass ground state calculation with `-r` argument as:

	gpawsolve.py -r -i MoS2-GW.py -g MoS2-structure.cif

Also you can use `-d` argument to see the band graph at the end of the calculation.
