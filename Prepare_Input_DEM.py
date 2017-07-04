"""
This is a Python script to prepare files for analysis"""

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
#Use the cookbook
#https://pcjericks.github.io/py-gdalogr-cookbook/

#------------------------------------------------------------------
#These are the resolutions we work on
Gauges=["HIN5", "HIN10", "HIN15", "HIN20", "HIN25", "HIN30", "HIN40", "HIN50", "HIN75", "HIN100"] # Sorted by resolution in dm
Gauges=["HIN20"] # Sorted by resolution in dm

Nodata_value = -9999 # This is the value for empty DEM cells


for gauge in Gauges:

    #### Run this from the following directory: (this is where your files are)
    # /home/s1563094/Datastore/Software/LSDTopoTools/LSDTopoTools_MarshExtraction/Input/Topography/HIN_degraded
    
    #print "Converting DEM files to ENVI format"
    #sourcefile = "%sdm_DEM.asc" % (gauge)
    #destinationfile = "%sdm_DEM.bil" % (gauge)
    #os.system("gdal_translate -of ENVI -a_srs EPSG:27700 " + sourcefile + " " +  destinationfile)



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























