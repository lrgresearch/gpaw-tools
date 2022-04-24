#!/usr/bin/env bash
echo "gpaw-tools: "
echo "Calculating all examples BASH script..."
CORENUMBER=4
SCRIPT=`realpath $0`
SCRIPTPATH=`dirname $SCRIPT`
# Examples
# Bulk-Al-noCIF -------------------
echo "Calculating: Bulk-Al-noCIF"
cd ./Bulk-Al-noCIF
time mpirun -np $CORENUMBER gpawsolve.py -o -i bulk_aluminum.py

# Cr2O-spin -------------------
echo "Calculating: Cr2O-spin"
cd ../Cr2O-spin
time mpirun -np $CORENUMBER gpawsolve.py -o -i Cr2O.py -g Cr2O_mp-1206821_primitive.cif

# Graphene-LCAO -------------
echo "Calculating: Graphene-LCAO"
cd ../Graphene-LCAO
echo "Step 1: Pristine graphene"
time mpirun -np $CORENUMBER gpawsolve.py -o -i graphene.py -g graphene4x4.cif
echo "Step 2: Graphene with defect"
time mpirun -np $CORENUMBER gpawsolve.py -o -i graphene.py -g graphene4x4withdefect.cif

# MoS2-GW -------------------
echo "Calculating: MoS2-GW"
cd ../MoS2-GW
time gpawsolve.py -o -i MoS2-GW.py -g MoS2-structure.cif

# Si-2atoms-optical ----------------
echo "Calculating: Si-2atoms-optical"
cd ../Si-2atoms-optical
echo "Step 1: Ground, DOS and Band"
time mpirun -np $CORENUMBER gpawsolve.py -o -i Si-Step1-ground_dos_band.py -g Si_mp-149_primitive_Example.cif
echo "Step 2: Optical - RPA"
time gpawsolve.py -o -i Si-Step2-optical-RPA.py -g Si_mp-149_primitive_Example.cif
echo "Step 3: Optical - BSE"
time gpawsolve.py -o -i Si-Step2-optical-BSE.py -g Si_mp-149_primitive_Example.cif

# Wurtzite ZnO with DFT+U
echo "Calculating: ZnO with DFT+U"
cd ../ZnO-with-Hubbard
echo "Calculating: Ground, DOS and Band with DFT+U"
time mpirun -np $CORENUMBER gpawsolve.py -o -i ZnO_withHubbard.py
echo "Calculating: Ground, DOS and Band without DFT+U"
time mpirun -np $CORENUMBER gpawsolve.py -o -i ZnO_woHubbard.py

# Rocksalt TiC with Elastic Calculations
echo "Calculating: Rocksalt TiC"
cd ../TiC-elastic-electronic
time mpirun -np $CORENUMBER gpawsolve.py -o -i TiC.py -g TiC_mp-631_primitive-Final.cif

# Finish
echo "All calculations are finished."
