"""
This is a Python script to prepare DEM files for the automated topographic detection of saltmarshes. This file is particularly important and performs the following actions:
1. Convert your DEM into the ENVI format used by LSDTopoTools
2. Apply a Wiener filter (optional)
3. Calculate basic topographic analysis (hillshade, slope and curvature)
"""

#------------------------------------------------------------------
#0. Set up display environment in putty if you are working on a terminal with no graphical interface.
import matplotlib
matplotlib.use('Agg')


#------------------------------------------------------------------
#1. Load useful Python packages
import os
import sys


#------------------------------------------------------------------
#2. Set up the important variables

#Select the file names (or partial file names) you are working with
Gauges=["HIN5", "HIN10", "HIN15", "HIN20", "HIN25", "HIN30", "HIN40", "HIN50", "HIN75", "HIN100"] #resolution in dm

# Set up value for empty DEM cells
Nodata_value = -9999

#------------------------------------------------------------------
#3. Run the preparation script. 

"""NOTE: in the test version make sure you only need to run this script once"""

for gauge in Gauges:

    #### Run this from the following directory: (this is where your files are)
    # /home/s1563094/Datastore/Software/LSDTopoTools/LSDTopoTools_MarshExtraction/Input/Topography/HIN_degraded
    
    #print "Converting DEM files to ENVI format"
    #sourcefile = "%sdm_DEM.asc" % (gauge)
    #destinationfile = "%sdm_DEM.bil" % (gauge)
    #os.system("gdal_translate -of ENVI -a_srs EPSG:27700 " + sourcefile + " " +  destinationfile)
    #Useful page for further use of GDAL: https://pcjericks.github.io/py-gdalogr-cookbook/
    

    #### Run this from the following directory:
    # /home/s1563094/Datastore/Software/LSDTopoTools/LSDTopoTools_ChannelExtraction/driver_functions_ChannelExtraction
    
    #print "Applying wiener filter"
    #sourcedir = "/home/s1563094/Datastore/Software/LSDTopoTools/LSDTopoTools_MarshExtraction/Input/Topography/HIN_degraded/ "
    #sourcefile = "%sdm_DEM " % (gauge)
    #sourceformat = "bil"
    #os.system("./Wiener_filter.out " + sourcedir + sourcefile +  sourceformat)


    #### Run this from the following directory:
    # /home/s1563094/Datastore/Software/LSDTopoTools/Git_projects/LSDTopoTools_AnalysisDriver/Analysis_driver
    
    print "Calculating slopes and curvatures"
    sourcedir = "/home/s1563094/Datastore/Software/LSDTopoTools/LSDTopoTools_MarshExtraction/Input/Topography/HIN_degraded/ "
    driverfile = "%s_slope_curv.LSDTT_driver " % (gauge)
    os.system("./LSDTT_analysis_from_paramfile.out " + sourcedir + driverfile)

