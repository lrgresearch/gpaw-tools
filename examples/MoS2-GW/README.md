# Example: GW Aproximation calculation for MoS2

To run the calculation, do not usi mpirun and use `gpawsolve.py` as

    gpawsolve.py -i MoS2-GW.py -g MoS2-structure.cif

After first run, you can by pass ground state calculation with changing the line `Ground_calc` as:

	Ground_calc = True     # Ground state calculations

Also you can use `-d` argument to see the band graph at the end of the calculation.
