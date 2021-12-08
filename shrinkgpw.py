#!/usr/bin/env python
 
'''
shrinkgpw.py: Extracting wave functions from gpw file to get smaller files

Usage: $ shrinkgpw.py <file.gpw>
'''
from ase import *
from ase.io import read
from gpaw import GPAW
import sys, os
from pathlib import Path

if len(sys.argv) > 1:
    inFile = os.path.join(os.getcwd(),sys.argv[1])
else:
    print("Please provide a GPW file to continue.")
    exit()

calc = GPAW(inFile)

# Writing result
calc.write(os.path.join(os.getcwd(),Path(inFile).stem+'_withoutwf.gpw'))
