"""
This is a Python script that identifies marsh platforms from a DEM and properties (slopes, curvature)
"""

#Set up display environment in putty
import matplotlib
matplotlib.use('Agg')

#----------------------------------------------------------------
#1. Load useful Python packages
import os
import sys
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import numpy as np
import functools
import math as mt
import cmath
import scipy as sp
import scipy.stats as stats
from datetime import datetime
import cPickle
from matplotlib import cm
from pylab import *
import functools
import matplotlib.ticker as tk
from matplotlib import rcParams
from mpl_toolkits.axes_grid1.inset_locator import *
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import itertools as itt
from osgeo import gdal, osr
import matplotlib.ticker as tk
from matplotlib import rcParams
from mpl_toolkits.axes_grid1.inset_locator import *
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from osgeo import gdal, gdalconst
from osgeo.gdalconst import *

#------------------------------------------------------------------
# Import the cosutm-made marsh-finding functions
from DEM_functions import Open_tide_stats
from DEM_functions import MARSH_ID
from DEM_functions import Confusion
from DEM_functions import METRIC1_DEM
from DEM_functions import METRIC2_DEM
from DEM_functions import METRIC3_DEM
from DEM_functions import ENVI_raster_binary_to_2d_array
from DEM_functions import ENVI_raster_binary_from_2d_array


#------------------------------------------------------------------
#Use the cookbook
#https://pcjericks.github.io/py-gdalogr-cookbook/

#------------------------------------------------------------------
#These are the resolutions we work on

Gauges=["HIN5", "HIN10", "HIN15", "HIN20", "HIN25", "HIN30", "HIN40", "HIN50", "HIN75"]#, "HIN100"] # Sorted by resolution in dm
#Gauges=["HIN30", "HIN40", "HIN50", "HIN75", "HIN100"] # Sorted by resolution in dm

Resolutions = [0.5,1,1.5,2,2.5,3,4,5,7.5,10]


Nodata_value = -9999 # This is the value for empty DEM cells




##################


# FOR THIS YOU NEED TO COMPARE SHAPEFILES....

# OR DO YOU? YOU COULD JUST UPGRADE THE RESOLUTION OF THE OUTPUT RASTER

# Something's wrong

##################



# THIS IS THE NO-SALTINGS VERSION
print "Preparing reference data at a high resolution"

print " Converting reference marsh outline to a raster"
sourcefile = "Input/Reference/HIN_degraded/HIN05_marsh_nosaltings_DEM.shp"
destinationfile = "Input/Reference/HIN_degraded/HIN05_marsh_nosaltings_DEM.bil"
os.system("gdal_rasterize -a Raster_val -of ENVI -a_srs EPSG:27700 -tr 0.5 0.5 " + sourcefile + " " +  destinationfile)

print " Clipping reference marsh outline raster"
sourcefile = "Input/Reference/HIN_degraded/HIN05_marsh_nosaltings_DEM.bil"
cutfile = "Input/Reference/HIN_degraded/HIN_domain.shp"
destinationfile = "Input/Reference/HIN_degraded/HIN05_marsh_nosaltings_DEM_clip.bil"
os.system("gdalwarp -overwrite -of ENVI -cutline " + cutfile + " -crop_to_cutline " + sourcefile + " " +  destinationfile)

print " Loading reference marsh outline raster"
Reference, post_Reference, envidata_Reference =  ENVI_raster_binary_to_2d_array ("Input/Reference/HIN_degraded/HIN05_marsh_nosaltings_DEM_clip.bil", "HIN05")


for gauge in Gauges:
    print "Let's process site %s" % (gauge)
    print

    print "Preparing input data"
    cutfile = "Input/Reference/HIN_degraded/HIN_domain.shp"
    
    print " Clipping DEM raster"
    sourcefile = "Input/Topography/HIN_degraded/%sdm_DEM_WFILT.bil" % (gauge)
    destinationfile = "Input/Topography/HIN_degraded/%sdm_DEM_WFILT_clip.bil" % (gauge)
    os.system("gdalwarp -overwrite -of ENVI -cutline " + cutfile + " -crop_to_cutline " + sourcefile + " " +  destinationfile)

    print " Clipping slope raster"
    sourcefile = "Input/Topography/HIN_degraded/%sdm_slope.bil" % (gauge)
    destinationfile = "Input/Topography/HIN_degraded/%sdm_slope_clip.bil" % (gauge)
    os.system("gdalwarp -overwrite -of ENVI -cutline " + cutfile + " -crop_to_cutline " + sourcefile + " " +  destinationfile)
    
    print " Clipping curvature raster"
    sourcefile = "Input/Topography/HIN_degraded/%sdm_curvature.bil" % (gauge)
    destinationfile = "Input/Topography/HIN_degraded/%sdm_curvature_clip.bil" % (gauge)
    os.system("gdalwarp -overwrite -of ENVI -cutline " + cutfile + " -crop_to_cutline " + sourcefile + " " +  destinationfile)
    
    print " Clipping hillshade raster"
    sourcefile = "Input/Topography/HIN_degraded/%sdm_hs.bil" % (gauge)
    destinationfile = "Input/Topography/HIN_degraded/%sdm_hs_clip.bil" % (gauge)
    os.system("gdalwarp -overwrite -of ENVI -cutline " + cutfile + " -crop_to_cutline " + sourcefile + " " +  destinationfile)


    print "Loading input data"
    print " Loading tidalstatistix"
    Metric1_tide, Metric2_tide, Metric3_tide, Subsample = Open_tide_stats ("Input/Tide/HIN/HIN_", "HIN")
    print " Loading DEM"
    DEM, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array ("Input/Topography/HIN_degraded/%sdm_DEM_WFILT_clip.bil" % (gauge), gauge)
    print " Loading Slopes"
    Slope, post_Slope, envidata_Slope =  ENVI_raster_binary_to_2d_array ("Input/Topography/HIN_degraded/%sdm_slope_clip.bil" % (gauge), gauge)
    print " Loading Curvature"
    Curvature, post_Curvature, envidata_Curvature =  ENVI_raster_binary_to_2d_array ("Input/Topography/HIN_degraded/%sdm_curvature_clip.bil" % (gauge), gauge)
    #print " Loading Channels"
    #Channels, post_Channels, envidata_Channels =  ENVI_raster_binary_to_2d_array ("LiDAR_DTM_1m/%s/%s_Channels_SO_wiener.bil" % (gauge,gauge), gauge)



    print "Identifying the platform and scarps"
    DEM_work = np.copy(DEM)/1000
    Slope_work = np.copy(Slope)/1000
    Search_space, Scarps, Platform = MARSH_ID(DEM_work, Slope_work, Curvature, Metric2_tide, Nodata_value)
    Scarps[Scarps == 0] = Nodata_value
    
        #------------------------------------------------------------------------------------------------------
    #Save the results
    print "Saving marsh features"
    new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, "Output/HIN_degraded/%s_Search_space.bil" % (gauge), post_DEM, Search_space) 
    new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, "Output/HIN_degraded/%s_Scarps.bil" % (gauge), post_DEM, Scarps)
    new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, "Output/HIN_degraded/%s_Marsh.bil" % (gauge), post_DEM, Platform)
    
    
    # If the resolution is different from the reference, upgrade raster resolution
    if gauge != "HIN5":
        print "Upgrading output raster resolution to match the reference"
        sourcefile = "Output/HIN_degraded/%s_Marsh.bil" % (gauge)
        destinationfile = "Output/HIN_degraded/%s_Marsh.bil" % (gauge)
        os.system("gdal_translate -of ENVI -a_srs EPSG:27700 -tr 0.5 0.5 " + sourcefile + " " +  destinationfile)


    print "Measuring performances"
    
    print " Loading Marsh"
    Platform_work, post_Platform, envidata_Platform =  ENVI_raster_binary_to_2d_array ("Output/HIN_degraded/%s_Marsh.bil" % (gauge), gauge)
    
    Confusion_matrix, Performance, Metrix = Confusion (Platform_work, Reference, Nodata_value)

    new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, "Output/HIN_degraded/%s_Confusion_nosaltings_DEM.bil" % (gauge), post_DEM, Confusion_matrix)

    cPickle.dump(Performance,open("Output/HIN_degraded/%s_nosaltings_Performance.pkl" % (gauge), "wb"))
    cPickle.dump(Metrix,open("Output/HIN_degraded/%s_nosaltings_Metrix.pkl" % (gauge), "wb"))

    print
    print
