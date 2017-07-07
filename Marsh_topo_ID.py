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
#These are the tide gauges next to the marshes we work on
Gauges=["BOU", "FEL", "CRO", "SHE", "WOR", "HEY", "HIN"] # Sorted by Tidal Range
Gauges=["HIN"] # Sorted by Tidal Range

Nodata_value = -9999 # This is the value for empty DEM cells


for gauge in Gauges:
    print "Let's process site %s" % (gauge)
    print

    print "Preparing reference data"
    print " Converting reference marsh outline to a raster"
    sourcefile = "Input/Reference/%s/%s_marsh_DEM.shp" % (gauge,gauge)
    destinationfile = "Input/Reference/%s/%s_marsh_DEM.bil"  % (gauge,gauge)
    os.system("gdal_rasterize -a Raster_val -of ENVI -a_srs EPSG:27700 -tr 1 1 " + sourcefile + " " +  destinationfile)
    print " Clipping reference marsh outline raster"
    sourcefile = "Input/Reference/%s/%s_marsh_DEM.bil" % (gauge,gauge)
    cutfile = "Input/Reference/%s/%s_domain.shp" % (gauge,gauge)
    destinationfile = "Input/Reference/%s/%s_marsh_DEM_clip.bil"  % (gauge,gauge)
    os.system("gdalwarp -overwrite -of ENVI -cutline " + cutfile + " -crop_to_cutline " + sourcefile + " " +  destinationfile)
    print " Loading reference marsh outline raster"
    Reference, post_Reference, envidata_Reference =  ENVI_raster_binary_to_2d_array ("Input/Reference/%s/%s_marsh_DEM_clip.bil" % (gauge,gauge), gauge)


    """print "Preparing input data"
    print " Clipping DEM raster"
    sourcefile = "Input/Topography/%s/%s_DEM_WFILT.bil" % (gauge,gauge)
    cutfile = "Input/Reference/%s/%s_domain.shp" % (gauge,gauge)
    destinationfile = "Input/Topography/%s/%s_DEM_WFILT.bil" % (gauge,gauge)
    os.system("gdalwarp -overwrite -of ENVI -cutline " + cutfile + " -crop_to_cutline " + sourcefile + " " +  destinationfile)
    print " Clipping slope raster"
    sourcefile = "Input/Topography/%s/%s_slope.bil" % (gauge,gauge)
    cutfile = "Input/Reference/%s/%s_domain.shp" % (gauge,gauge)
    destinationfile = "Input/Topography/%s/%s_slope.bil" % (gauge,gauge)
    os.system("gdalwarp -overwrite -of ENVI -cutline " + cutfile + " -crop_to_cutline " + sourcefile + " " +  destinationfile)
    print " Clipping curvature raster"
    sourcefile = "Input/Topography/%s/%s_curvature.bil" % (gauge,gauge)
    cutfile = "Input/Reference/%s/%s_domain.shp" % (gauge,gauge)
    destinationfile = "Input/Topography/%s/%s_curvature.bil" % (gauge,gauge)
    os.system("gdalwarp -overwrite -of ENVI -cutline " + cutfile + " -crop_to_cutline " + sourcefile + " " +  destinationfile)
    print " Clipping hillshade raster"
    sourcefile = "Input/Topography/%s/%s_hs.bil" % (gauge,gauge)
    cutfile = "Input/Reference/%s/%s_domain.shp" % (gauge,gauge)
    destinationfile = "Input/Topography/%s/%s_hs.bil" % (gauge,gauge)
    os.system("gdalwarp -overwrite -of ENVI -cutline " + cutfile + " -crop_to_cutline " + sourcefile + " " +  destinationfile)"""


    print "Loading input data"
    print " Loading tidalstatistix"
    Metric1_tide, Metric2_tide, Metric3_tide, Subsample = Open_tide_stats ("Input/Tide/%s/%s_" % (gauge,gauge), gauge)
    print " Loading DEM"
    DEM, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array ("Input/Topography/%s/%s_DEM_WFILT.bil" % (gauge,gauge), gauge)
    print " Loading Slopes"
    Slope, post_Slope, envidata_Slope =  ENVI_raster_binary_to_2d_array ("Input/Topography/%s/%s_slope.bil" % (gauge,gauge), gauge)
    print " Loading Curvature"
    Curvature, post_Curvature, envidata_Curvature =  ENVI_raster_binary_to_2d_array ("Input/Topography/%s/%s_curvature.bil" % (gauge,gauge), gauge)
    #print " Loading Channels"
    #Channels, post_Channels, envidata_Channels =  ENVI_raster_binary_to_2d_array ("LiDAR_DTM_1m/%s/%s_Channels_SO_wiener.bil" % (gauge,gauge), gauge)



    print "Identifying the platform and scarps"
    DEM_work = np.copy(DEM)
    Search_space, Scarps, Platform = MARSH_ID(DEM_work, Slope, Curvature, Metric2_tide, Nodata_value)
    Platform_work = np.copy(Platform)

    print "Measuring performances"
    Confusion_matrix, Performance, Metrix = Confusion (Platform_work, Reference, Nodata_value)

    Scarps[Scarps == 0] = Nodata_value
    #------------------------------------------------------------------------------------------------------
    #Save the results
    print "Saving marsh features"
    new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, "Output/%s/%s_Search_spaceB.bil" % (gauge, gauge), post_DEM, Search_space)
    new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, "Output/%s/%s_ScarpsB.bil" % (gauge, gauge), post_DEM, Scarps)
    new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, "Output/%s/%s_MarshB.bil" % (gauge, gauge), post_DEM, Platform)
    new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, "Output/%s/%s_Confusion_DEMB.bil" % (gauge, gauge), post_DEM, Confusion_matrix)

    cPickle.dump(Performance,open("Output/%s/%s_PerformanceA.pkl" % (gauge,gauge), "wb"))
    cPickle.dump(Metrix,open("Output/%s/%s_MetrixA.pkl" % (gauge,gauge), "wb"))

    print
    print
