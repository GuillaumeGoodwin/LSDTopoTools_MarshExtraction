""""
This Python script converts multiple .bil DEMs into .xyz files
Those files are OK for use in almost all applications
To use with Incidenze you must open all the files with Excel and save them, then change the extension to .xyz
"""
#----------------------------------------------------------------
#1. Load useful Python packages

import os
import sys


Gauges=["BOU", "FEL", "CRO", "SHE", "WOR", "HEY", "HIN"]
#------------------------------------------------------------------
#STEP1
"""for gauge in Gauges:
    SRC_file = 'Input/Topography/%s/%s_DEM_WFILT.bil' % (gauge, gauge)
    DST_file = 'Input/Topography/%s_DEM.xyz' % (gauge)
    os.system('gdal_translate -of XYZ ' + SRC_file + ' ' + DST_file)

    crs = open('Input/Topography/%s_DEM.xyz' % (gauge), "r")
    dst = open('Input/Topography/%s_DEM.csv' % (gauge), "w")
    for columns in (raw.strip().split() for raw in crs ):
        A = columns[0] + ' ' + columns[1] + ' ' + columns[2] + '\n'
        dst.write(A)
    dst.close"""

#------------------------------------------------------------------
#STEP2 : comment STEP 1 before running
for gauge in Gauges:
    SRC_file = 'Input/Topography/%s_DEM.csv' % (gauge)
    DST_file = 'Input/Topography/%s_DEM.xyz' % (gauge)
    os.system('cp ' + SRC_file + ' ' + DST_file)
