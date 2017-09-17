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
def make_driver(templatedir,templatefile, gauge, res, Wfilter):
    src = open(templatedir+templatefile, "r")
    dst = open(templatedir+"AAA.txt", "w")
    counter = 0
    for columns in (raw.strip().split() for raw in src):
        line = ''
        for i in range(len(columns)):
            line = line + ' ' + columns[i]
        A = line + '\n'
        if counter == 8:
            A = "read path: /exports/csce/datastore/geos/users/s1563094/Software/LSDTopoTools/LSDTopoTools_MarshExtraction/Input/Topography/%s" % (gauge)  + '\n'
        if counter == 9:
            A = "write path: /exports/csce/datastore/geos/users/s1563094/Software/LSDTopoTools/LSDTopoTools_MarshExtraction/Input/Topography/%s" % (gauge)  + '\n'
        if counter == 10:
            if Wfilter == 1:
                A = "read fname: %s_%s_DEM_clip_WFILT " % (gauge,res)  + '\n'
            elif Wfilter == 0:
                A = "read fname: %s_%s_DEM_clip " % (gauge,res)  + '\n'
        if counter == 11:
            A = "write fname: %s_%s" % (gauge,res) + '\n'
        dst.write(A)
        counter = counter + 1
    dst.close



#------------------------------------------------------------------
# This is important
Nodata_value = -9999 # This is the value for empty DEM cells


#------------------------------------------------------------------
#These are the tide gauges next to the marshes we work on, sorted by Tidal Range
Gauges=["BOU","FEL", "CRO", "SHE", "HEY", "HIN"] 
Gauges=["BOU","FEL", "CRO"] 



# Here are the 3 optimisation points
Opt1=["-2.0", "-1.8", "-1.6", "-1.4", "-1.2", "-1.0", "-0.8", "-0.6", "-0.4", "-0.2"] 
Opt1=["-2.0"] # The optimal value, equals with -1.8
Opt1_num = np.asarray(Opt1).astype(float)

Opt2=["0.50", "0.55", "0.60", "0.65", "0.70", "0.75", "0.80", "0.85", "0.90", "0.95"]  
Opt2=["0.85"] # The optimal value, equals with 0.80
Opt2_num = np.asarray(Opt2).astype(float)


Opt3=["2.0", "4.0", "6.0", "8.0", "10.0", "12.0", "14.0", "16.0", "18.0", "20.0"] # This is necessarily a round number, but please keep the float format
Opt3_num = np.asarray(Opt3).astype(float)


#------------------------------------------------------------------
# First, you need to prepare the data


for gauge in Gauges:
    
    # Start with your reference data, lovingly digitised in your favourite GIS software (i.e. not Arc)
    print "Preparing reference data"
    directory = "Input/Reference/%s/" % (gauge)
    
    print " Converting reference marsh outline to a raster"
    srcfile = "%s_marsh_DEM2.shp" % (gauge)
    dstfile = "%s_marsh_DEM2.bil"  % (gauge)
    os.system("gdal_rasterize -a Raster_val -of ENVI -a_srs EPSG:27700 -tr 1 1 " + directory+srcfile + " " +  directory+dstfile)
    
    print " Clipping reference marsh outline raster"
    srcfile = "%s_marsh_DEM2.bil" % (gauge)
    cutfile = "%s_domain2.shp" % (gauge)
    dstfile = "%s_ref_DEM_clip.bil"  % (gauge)
    os.system("gdalwarp -overwrite -of ENVI -cutline " + directory+cutfile + " -crop_to_cutline " + directory+srcfile + " " +  directory+dstfile)
    
    print " Loading reference marsh outline raster"
    Reference, post_Reference, envidata_Reference =  ENVI_raster_binary_to_2d_array (directory+dstfile, gauge)

    
    """# Then take care of your input DEM
    print "Converting Input DEM files to ENVI format of the right domain"
    directory = "Input/Topography/%s/" % (gauge)
    cutfile = "Input/Reference/%s/%s_domain.shp" % (gauge,gauge)
    srcfile = "%s_trans.asc" % (gauge)
    dstfile = "%s_DEM_clip.bil" % (gauge)
    os.system("gdalwarp -of ENVI -t_srs EPSG:27700 -cutline " + cutfile + " -crop_to_cutline " + directory+srcfile + " " +  directory+dstfile)
    
    # Then take care of your input Slope
    print "Converting Input DEM files to ENVI format of the right domain"
    directory = "Input/Topography/%s/" % (gauge)
    cutfile = "Input/Reference/%s/%s_domain.shp" % (gauge,gauge)
    srcfile = "%s_slope_big.bil" % (gauge)
    dstfile = "%s_slope.bil" % (gauge)
    os.system("gdalwarp -of ENVI -t_srs EPSG:27700 -cutline " + cutfile + " -crop_to_cutline " + directory+srcfile + " " +  directory+dstfile)
    
    # Then take care of your input Curvature
    print "Converting Input DEM files to ENVI format of the right domain"
    directory = "Input/Topography/%s/" % (gauge)
    cutfile = "Input/Reference/%s/%s_domain.shp" % (gauge,gauge)
    srcfile = "%s_curvature_big.bil" % (gauge)
    dstfile = "%s_curvature.bil" % (gauge)
    os.system("gdalwarp -of ENVI -t_srs EPSG:27700 -cutline " + cutfile + " -crop_to_cutline " + directory+srcfile + " " +  directory+dstfile)"""

    
    
    print "Loading input data"
    print " Loading tidalstatistix"
    Metric1_tide, Metric2_tide, Metric3_tide, Subsample = Open_tide_stats ("Input/Tide/%s/%s_" % (gauge,gauge), gauge)
    print " Loading DEM"
    DEM, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array ("Input/Topography/%s/%s_DEM_clip.bil" % (gauge,gauge), gauge)
    #DEM, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array ("Input/Topography/%s/%s_%s_DEM_clip_WFILT.bil" % (gauge,gauge,res), gauge)
    print " Loading Slopes"
    Slope, post_Slope, envidata_Slope =  ENVI_raster_binary_to_2d_array ("Input/Topography/%s/%s_slope.bil" % (gauge,gauge), gauge)

    #This is just to sort out our crappy DEM in WOR
    Slope[Slope>1] = 0

    print " Loading Curvature"
    Curvature, post_Curvature, envidata_Curvature =  ENVI_raster_binary_to_2d_array ("Input/Topography/%s/%s_curvature.bil" % (gauge,gauge), gauge)
    #print " Loading Channels"
    #Channels, post_Channels, envidata_Channels =  ENVI_raster_binary_to_2d_array ("LiDAR_DTM_1m/%s/%s_%s_Channels_SO_wiener.bil" % (gauge,gauge,res), gauge)

    
    i=0
    for opt1 in Opt1:
        opt1_num = float(opt1)
        
        for opt2 in Opt2:
            opt2_num = float(opt2)
            
            for opt3 in Opt3:
                opt3_num = float(opt3)
    
                print "Identifying the platform and scarps"
                DEM_work = np.copy(DEM)
                Search_space, Scarps, Platform = MARSH_ID(DEM_work, Slope, Curvature, Metric2_tide, Nodata_value, opt1_num, opt2_num, opt3_num)
                Platform_work = np.copy(Platform)
                Scarps[Scarps == 0] = Nodata_value


                #------------------------------------------------------------------------------------------------------
                #Save the results
                print "Saving marsh features"
                new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, "Output/%s/%s_O1_%s_O2_%s_O3_%s_Search_space_nofilter.bil" % (gauge, gauge, opt1, opt2, opt3), post_DEM, Search_space)

                new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, "Output/%s/%s_O1_%s_O2_%s_O3_%s_Scarps_nofilter.bil" % (gauge, gauge,opt1, opt2, opt3), post_DEM, Scarps)

                new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, "Output/%s/%s_O1_%s_O2_%s_O3_%s_Marsh_nofilter.bil" % (gauge, gauge,opt1, opt2, opt3), post_DEM, Platform)

                #---------------------------------------------------------------
                # Calculate the performances of your algorithm


                print " Loading Marsh"
                Platform_work, post_Platform, envidata_Platform =  ENVI_raster_binary_to_2d_array ("Output/%s/%s_O1_%s_O2_%s_O3_%s_Marsh_nofilter.bil" % (gauge,gauge, opt1, opt2, opt3), gauge)

                print "Measuring performances"
                Confusion_matrix, Performance, Metrix = Confusion (Platform_work, Reference, Nodata_value)
                new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_Platform, "Output/%s/%s_O1_%s_O2_%s_O3_%s_Confusion_DEM_nofilter.bil" % (gauge, gauge, opt1, opt2, opt3), post_Platform, Confusion_matrix)

                cPickle.dump(Performance,open("Output/%s/%s_O1_%s_O2_%s_O3_%s_Performance_nofilter.pkl" % (gauge,gauge,opt1, opt2, opt3), "wb"))
                cPickle.dump(Metrix,open("Output/%s/%s_O1_%s_O2_%s_O3_%s_Metrix_nofilter.pkl" % (gauge,gauge, opt1, opt2, opt3), "wb"))



        i=i+1




