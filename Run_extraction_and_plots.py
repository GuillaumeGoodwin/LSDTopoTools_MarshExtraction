"""
Motherscript.py

This Python script runs all the scripts necessary to:
1. Extract marsh platforms and outlines from an input DEM
2. Save these features as .bil files
3. Plot the figures of our paper, including a comparison with a reference DEM
"""

# First import the mecessary modules
import os
import sys

# Then run the scripts
os.system("python Marsh_topo_ID.py")
os.system("python Plots.py")