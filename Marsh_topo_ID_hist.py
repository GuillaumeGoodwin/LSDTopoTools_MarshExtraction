"""
This is a Python script. It does quite a lot of useful things:

1. If you give it a DEM, it will prepare it for analysis

2. It will then find the marsh platforms and scarps inside your DEM

"""

#Set up display environment in putty
import matplotlib
matplotlib.use('Agg')

#----------------------------------------------------------------
# Load useful Python packages
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
import fileinput, sys

#------------------------------------------------------------------
# Import the costum-made marsh-finding functions
from DEM_functions import Open_tide_stats
from DEM_functions import MARSH_ID
from DEM_functions import Confusion
from DEM_functions import METRIC1_DEM
from DEM_functions import METRIC2_DEM
from DEM_functions import METRIC3_DEM
from DEM_functions import ENVI_raster_binary_to_2d_array
from DEM_functions import ENVI_raster_binary_from_2d_array


#-------------------------------------------------------------------
# This is a function to make a driver file
def make_driver(templatedir,templatefile, gauge, Wfilter):
    src = open(templatedir+templatefile, "r")
    dst = open(templatedir+"AAA.txt", "w")
    counter = 0
    for columns in (raw.strip().split() for raw in src):
        line = ''
        for i in range(len(columns)):
            line = line + ' ' + columns[i]
        A = line + '\n'
        if counter == 8:
            A = "read path: /exports/csce/datastore/geos/users/s1563094/Software/LSDTopoTools/LSDTopoTools_MarshExtraction/Input/Topography/HIN_hist" + '\n'
        if counter == 9:
            A = "write path: /exports/csce/datastore/geos/users/s1563094/Software/LSDTopoTools/LSDTopoTools_MarshExtraction/Input/Topography/HIN_hist" + '\n'
        if counter == 10:
            if Wfilter == 1:
                A = "read fname: %s_DEM_clip_WFILT " % (gauge)  + '\n'
            elif Wfilter == 0:
                A = "read fname: %s_DEM_clip " % (gauge)  + '\n'
        if counter == 11:
            A = "write fname: %s" % (gauge) + '\n'
        dst.write(A)
        counter = counter + 1
    dst.close


#------------------------------------------------------------------
#Use the cookbook
#https://pcjericks.github.io/py-gdalogr-cookbook/

#------------------------------------------------------------------
# This is important
Nodata_value = -9999 # This is the value for empty DEM cells


#------------------------------------------------------------------
#These are the tide gauges next to the marshes we work on, sorted by Tidal Range
Gauges=["HIN_200703","HIN_200710", "HIN_200909", "HIN_201103", "HIN_201205", "HIN_201302", "HIN_201402","HIN_201410"] 

#------------------------------------------------------------------
# First, you need to prepare the data


for gauge in Gauges:
    
    # Start with your reference data, lovingly digitised in your favourite GIS software (i.e. not Arc)
    print "Preparing reference data"
    directory = "Input/Reference/HIN_hist/"
    
    print " Converting reference marsh outline to a raster"
    srcfile = "%s_marsh_DEM.shp" % (gauge)
    dstfile = "%s_marsh_DEM.bil"  % (gauge)
    os.system("gdal_rasterize -a Raster_val -of ENVI -a_srs EPSG:27700 -tr 1 1 " + directory+srcfile + " " +  directory+dstfile)
    
    print " Clipping reference marsh outline raster"
    srcfile = "%s_marsh_DEM.bil" % (gauge)
    cutfile = "HIN_domain.shp"
    dstfile = "%s_ref_DEM_clip.bil"  % (gauge)
    os.system("gdalwarp -overwrite -of ENVI -cutline " + directory+cutfile + " -crop_to_cutline " + directory+srcfile + " " +  directory+dstfile)
    
    print " Loading reference marsh outline raster"
    Reference, post_Reference, envidata_Reference =  ENVI_raster_binary_to_2d_array (directory+dstfile, gauge)

    
    
    # Then take care of your input data
    print "Converting Input DEM files to ENVI format of the right domain"
    directory = "Input/Topography/HIN_hist/"
    cutfile = "Input/Reference/HIN_hist/HIN_domain.shp"
    srcfile = "%s.tif" % (gauge)
    dstfile = "%s_DEM_clip.bil" % (gauge)
    os.system("gdalwarp -of ENVI -t_srs EPSG:27700 -cutline " + cutfile + " -crop_to_cutline " + directory+srcfile + " " +  directory+dstfile)

    
    """print "Applying wiener filter"
    workdir = "/home/s1563094/Datastore/Software/LSDTopoTools/LSDTopoTools_ChannelExtraction/driver_functions_ChannelExtraction"
    srcdir = "/home/s1563094/Datastore/Software/LSDTopoTools/LSDTopoTools_MarshExtraction/Input/Topography/HIN_hist/ "
    srcfile = "%s_DEM_clip " % (gauge)
    srcformat = "bil"
    os.system("(cd "+ workdir +" && ./Wiener_filter.out " + srcdir + srcfile +  srcformat + ")")"""


    print "Calculating slopes and curvatures"
    workdir = "/home/s1563094/Datastore/Software/LSDTopoTools/Git_projects/LSDTopoTools_AnalysisDriver/Analysis_driver"
    srcdir = "/home/s1563094/Datastore/Software/LSDTopoTools/LSDTopoTools_MarshExtraction/Input/Topography/HIN_hist/"
    srcdir_short = "Input/Topography/HIN_hist/"
    templatedir = "/home/s1563094/Datastore/Software/LSDTopoTools/LSDTopoTools_MarshExtraction/"
    templatefile = "Template_driver.LSDTT_driver"
    driverfile = "%s_driver.LSDTT_driver " % (gauge)

    print " Making the driver file"
    make_driver(templatedir,templatefile, gauge,0)
    os.system("mv " + templatedir + "AAA.txt " + srcdir+driverfile)
    os.system("(cd "+ workdir +" && ./LSDTT_analysis_from_paramfile.out " + srcdir + " " + driverfile + ")")


    print "Loading input data"
    print " Loading tidalstatistix"
    Metric1_tide, Metric2_tide, Metric3_tide, Subsample = Open_tide_stats ("Input/Tide/HIN/HIN_",gauge)
    print " Loading DEM"
    DEM, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array ("Input/Topography/HIN_hist/%s_DEM_clip.bil" % (gauge), gauge)
    print " Loading Slopes"
    Slope, post_Slope, envidata_Slope =  ENVI_raster_binary_to_2d_array ("Input/Topography/HIN_hist/%s_slope.bil" % (gauge), gauge)

    #This is just to sort out our crappy DEM in WOR
    # Change the NaNs into nodata
    DEM[np.isnan(DEM)] = Nodata_value
    Slope[np.isnan(Slope)] = Nodata_value
    Slope[Slope>1] = 0

    print " Loading Curvature"
    Curvature, post_Curvature, envidata_Curvature =  ENVI_raster_binary_to_2d_array ("Input/Topography/HIN_hist/%s_curvature.bil" % (gauge), gauge)
    #print " Loading Channels"
    #Channels, post_Channels, envidata_Channels =  ENVI_raster_binary_to_2d_array ("LiDAR_DTM_1m/%s/%s_%s_Channels_SO_wiener.bil" % (gauge,gauge,res), gauge)


    

    # This is where the magic happens

    
    print "Identifying the platform and scarps"
    DEM_work = np.copy(DEM)
    Search_space, Scarps, Platform = MARSH_ID(DEM_work, Slope, Curvature, Metric2_tide, Nodata_value)
    Platform_work = np.copy(Platform)
    Scarps[Scarps == 0] = Nodata_value


    #------------------------------------------------------------------------------------------------------
    #Save the results
    print "Saving marsh features"
    new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, "Output/HIN_hist/%s_Search_space.bil" % (gauge), post_DEM, Search_space)
    new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, "Output/HIN_hist/%s_Scarps.bil" % (gauge), post_DEM, Scarps)
    new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, "Output/HIN_hist/%s_Marsh.bil" % (gauge), post_DEM, Platform)

    #---------------------------------------------------------------
    # Calculate the performances of your algorithm

    print " Loading Marsh"
    Platform_work, post_Platform, envidata_Platform =  ENVI_raster_binary_to_2d_array ("Output/HIN_hist/%s_Marsh.bil" % (gauge), gauge)

    print "Measuring performances"
    Confusion_matrix, Performance, Metrix = Confusion (Platform_work, Reference, Nodata_value)
    new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_Platform, "Output/HIN_hist/%s_Confusion_DEM.bil" % (gauge), post_Platform, Confusion_matrix)

    cPickle.dump(Performance,open("Output/HIN_hist/%s_Performance.pkl" % (gauge), "wb"))
    cPickle.dump(Metrix,open("Output/HIN_hist/%s_Metrix.pkl" % (gauge), "wb"))


