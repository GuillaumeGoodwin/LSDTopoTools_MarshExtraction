"""
This is a Python script plots the results of the MarshFinder


"""


#Set up display environment in putty
import matplotlib
matplotlib.use('Agg')

#----------------------------------------------------------------
#1. Load useful Python packages

import os
import sys


import numpy as np
import functools
import math as mt
import cmath
import scipy as sp
import scipy.stats as stats
from datetime import datetime
import cPickle

from pylab import *
import functools

import itertools as itt
from osgeo import gdal, osr
from osgeo import gdal, gdalconst
from osgeo.gdalconst import *
from copy import copy


from matplotlib import cm
import matplotlib.colors as colors
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import matplotlib.ticker as tk
from matplotlib import rcParams
from mpl_toolkits.axes_grid1.inset_locator import *
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from mpl_toolkits.basemap import Basemap, cm
from matplotlib.patches import Rectangle
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
#------------------------------------------------------------------
# Import the marsh-finding functions
from DEM_functions import Open_tide_stats
from DEM_functions import ENVI_raster_binary_to_2d_array
from DEM_functions import ENVI_raster_binary_from_2d_array
from DEM_functions import add_subplot_axes
from DEM_functions import Distribution

from DEM_functions import define_search_space



#------------------------------------------------------------------
#Use the cookbook
#https://pcjericks.github.io/py-gdalogr-cookbook/







# Straight out of the Python cookbook: https://stackoverflow.com/questions/26911898/matplotlib-curve-with-arrow-ticks
def add_arrow_to_line2D(
    axes, line, arrow_locs=[0.2, 0.4, 0.6, 0.8],
    arrowstyle='-|>', arrowsize=1, transform=None):
    """
    Add arrows to a matplotlib.lines.Line2D at selected locations.

    Parameters:
    -----------
    axes: 
    line: Line2D object as returned by plot command
    arrow_locs: list of locations where to insert arrows, % of total length
    arrowstyle: style of the arrow
    arrowsize: size of the arrow
    transform: a matplotlib transform instance, default to data coordinates

    Returns:
    --------
    arrows: list of arrows
    """
    if not isinstance(line, mlines.Line2D):
        raise ValueError("expected a matplotlib.lines.Line2D object")
    x, y = line.get_xdata(), line.get_ydata()

    arrow_kw = {
        "arrowstyle": arrowstyle,
        "mutation_scale": 10 * arrowsize,
    }

    color = line.get_color()
    use_multicolor_lines = isinstance(color, np.ndarray)
    if use_multicolor_lines:
        raise NotImplementedError("multicolor lines not supported")
    else:
        arrow_kw['color'] = color

    linewidth = line.get_linewidth()
    if isinstance(linewidth, np.ndarray):
        raise NotImplementedError("multiwidth lines not supported")
    else:
        arrow_kw['linewidth'] = linewidth

    if transform is None:
        transform = axes.transData

    arrows = []
    for loc in arrow_locs:
        s = np.cumsum(np.sqrt(np.diff(x) ** 2 + np.diff(y) ** 2))
        n = np.searchsorted(s, s[-1] * loc)
        arrow_tail = (x[n], y[n])
        arrow_head = (np.mean(x[n:n + 2]), np.mean(y[n:n + 2]))
        p = mpatches.FancyArrowPatch(
            arrow_tail, arrow_head, transform=transform,
            **arrow_kw)
        axes.add_patch(p)
        arrows.append(p)
    return arrows







#------------------------------------------------------------------
#These are the tide gauges next to the marshes we work on
#Gauges=["BOU", "FEL", "CRO", "SHE", "WOR", "HEY", "HIN"] # ALL gauges by tidal range
Gauges=["BOU", "FEL", "CRO", "SHE", "HEY", "HIN"] # ALL gauges by tidal range
# And these are the resolutions we work on, expressed in dm
Resolutions=["1.0", "1.5", "2.0", "2.5", "3.0", "4.0", "5.0", "7.5", "10.0"]  
Resolutions_num = np.asarray(Resolutions).astype(float)

#These are the gauges for evolution in time
HGauges=["HIN_200703","HIN_200710", "HIN_200909", "HIN_201103", "HIN_201205", "HIN_201302", "HIN_201402","HIN_201410"] 


#------------------------------------------------------------------
# These are some variables and arrays that need to be set up at the beginning
Nodata_value = -9999



#------------------------------------------------------------------
# This is where we load all the data we want to plot
"""i=0

for gauge in Gauges:
    print "Loading datasets for %s" % (gauge)
    print " Loading DEM"
    DEM, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array ("Input/Topography/%s/%s_DEM_WFILT.bil" % (gauge,gauge), gauge)
    print " Loading Slopes"
    Slope, post_Slope, envidata_Slope =  ENVI_raster_binary_to_2d_array ("Input/Topography/%s/%s_slope.bil" % (gauge,gauge), gauge)
    #print " Loading Curvature"
    #Curvature, post_Curvature, envidata_Curvature =  ENVI_raster_binary_to_2d_array ("Input/Topography/%s/%s_curvature.bil" % (gauge,gauge), gauge)
    #print " Loading Reference marsh"
    #Reference, post_Reference, envidata_Reference =  ENVI_raster_binary_to_2d_array ("Input/Reference/%s/%s_marsh_GE_clip.bil" % (gauge,gauge), gauge)
    #print " Loading Channels"
    #Channels, post_Channels, envidata_Channels =  ENVI_raster_binary_to_2d_array ("LiDAR_DTM_1m/%s/%s_Channels_SO_wiener.bil" % (gauge,gauge), gauge)


    print "Loading results for %s" % (gauge)
    print " Loading tidalstatistix"
    Metric1_tide, Metric2_tide, Metric3_tide, Subsample = Open_tide_stats ("Input/Tide/%s/%s_" % (gauge,gauge), gauge)
    print " Loading Search Space"
    Search_space, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array ("Output/%s/%s_Search_space.bil" % (gauge,gauge), gauge)
    print " Loading Scarps"
    Scarps, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array ("Output/%s/%s_Scarps.bil" % (gauge,gauge), gauge)
    print " Loading Platforms"
    
    """
    #Platform, post_Slope, envidata_Slope =  ENVI_raster_binary_to_2d_array ("Output/%s/%s_Marsh.bil" % (gauge,gauge), gauge)
    #print " Loading Confusion"
    #Confusion_matrix, post_Curvature, envidata_Curvature =  ENVI_raster_binary_to_2d_array ("Output/%s/%s_Confusion_DEM.bil" % (gauge,gauge), gauge)


    #print "Loading performance files"
    #with open ("Output/%s/%s_Performance.pkl" % (gauge,gauge), 'rb') as input_file:
        #Performance = cPickle.load(input_file)
    #with open ("Output/%s/%s_Metrix.pkl" % (gauge,gauge), 'rb') as input_file:
        #Metrix = cPickle.load(input_file)

        
        
    # Here you classify by tidal range
    #Metrix_gauges[i,0] = np.mean (Metric2_tide[3])-np.mean (Metric2_tide[0])
    
    # Here you classify by relief
    #Metrix_gauges[i,0] = np.amax(DEM) - np.amin(DEM[DEM>Nodata_value])
    
    #for j in np.arange(1,5,1):
        #Metrix_gauges[i,j] = Metrix[j-1]
    

    #i = i + 1
        
    

#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
#---------------------------------------------------------------------------
# Here are the paper plots

#---------------------------------------------------------------------------
# Figure 1 [1col]: This is a map of where our sites are from (a), complete with tidal range and distribution of elevations (b)
fig=plt.figure(1, facecolor='White',figsize=[3.2,4.9])
matplotlib.rc('xtick', labelsize=9) 

# First map the map
ax1 = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=1, axisbg='white')

# create Basemap
m = Basemap(llcrnrlon=-6.5,llcrnrlat=49.5,urcrnrlon=3.5,urcrnrlat=59, resolution='i', projection='cass', lon_0=-4.36, lat_0=54.7)
m.drawcoastlines()
m.drawparallels(np.arange(-40,61.,2.), labels=[1,0,0,0], fontsize = 9)
m.drawmeridians(np.arange(-20.,21.,2.), labels=[0,0,1,0], fontsize = 9)

# Plot the points on it
lats = [50.6, 51.8, 52.9, 51.5, 54.1, 51.15]
lons = [-1.9, 1.13, 0.8, 0.5, -2.8, -3.1]
Metrix_gauges = np.zeros((len(Gauges),5), dtype=np.float)


# Load the tidal data
i=0
for gauge in Gauges:
    Metric1_tide, Metric2_tide, Metric3_tide, Subsample = Open_tide_stats ("Input/Tide/%s/%s_" % (gauge,gauge), gauge)
    Metrix_gauges[i,0] = np.mean (Metric2_tide[3])-np.mean (Metric2_tide[0])
    if lons[i]>0:
        x, y = m(lons[i]-0.40,lats[i]+0.35)
    else:
        x, y = m(lons[i]+0.25,lats[i]-0.30)
    ax1.annotate('S%g' % (i+1), xy=(x,y), xycoords='data', fontsize=rcParams['font.size']-2, color='k')  
    i = i+1

# Plot the points
x, y = m(lons,lats)
Scatt = ax1.scatter(x,y, s = 50, color=plt.cm.winter(0.1*Metrix_gauges[:,0]), alpha = 0.9, linewidth = 5)

# Make a colourbar
ax2 = fig.add_axes([0.83, 0.10, 0.03, 0.8])
scheme = plt.cm.winter; norm = mpl.colors.Normalize(vmin=0, vmax=12)
bounds = [0,2,4,6,8,10,12]
cb = mpl.colorbar.ColorbarBase(ax2, cmap=scheme, norm=norm, ticks = bounds, orientation='vertical')
cb.set_label('Spring tidal range (m)', fontsize = 9)  
    


plt.savefig('Output/Paper/0_Main_Fig1.png')

 
# ADD COAST TYPE?




#---------------------------------------------------------------------------
# Figure 2 [1col]: This is a figure of the definition of the search space. It has examples of DEMxSlope (a) and the resulting search space (b) for an example
"""fig=plt.figure(2, facecolor='White',figsize=[3.2,6.5])

matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8) 


i=0
for gauge in Gauges:
    ax1 = plt.subplot2grid((12,5),(i,0),colspan=5, rowspan=1,axisbg='white')
    ax1.annotate('a%s.' % (i+1) , xy=(0.85,0.65), xycoords='axes fraction', fontsize=rcParams['font.size']-2, color='k') 
    
    majorLocatorx = MultipleLocator(0.02)
    majorFormatterx = FormatStrFormatter('%d')
    ax1.xaxis.set_major_locator(majorLocatorx)
    majorLocatory = MultipleLocator(0.1)
    majorFormattery = FormatStrFormatter('%d')
    ax1.yaxis.set_major_locator(majorLocatory)
    
    if i == 2:
        ax1.set_ylabel('Distribution', fontsize = 9)
    if i == 5:
        ax1.set_xlabel('X-value (m)', fontsize = 9)
    if i != 5:
        ax1.set_xticklabels([])
    ax1.set_xlim (0,0.08)
    ax1.set_ylim (0, 0.28)
    

    # Load the data
    DEM, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array ("Input/Topography/%s/%s_DEM_WFILT.bil" % (gauge,gauge), gauge)
    Slope, post_Slope, envidata_Slope =  ENVI_raster_binary_to_2d_array ("Input/Topography/%s/%s_slope.bil" % (gauge,gauge), gauge)
    Search_space, Crossover, bins, hist, Inflexion_point = define_search_space (DEM, Slope, Nodata_value)
    Search_space[Search_space==1] = Crossover[Search_space==1]

    ax1.plot(bins, hist, color=plt.cm.winter(0.1*Metrix_gauges[i,0]))
    Inflexion_index = np.where(hist < Inflexion_point)
    Inflexion_index = Inflexion_index[0]
    Inflexion_bin = bins[min(Inflexion_index)-1]
    ax1.fill_betweenx(bins, 0, Inflexion_bin, color='black', alpha = 0.5)
    
    if i == 0:

        #---------------------------
        # Maps now
        ax2 = plt.subplot2grid((12,5),(6,0), rowspan=6, colspan=4, axisbg='white')
        ax2.annotate('b.' , xy=(0.05,0.9), xycoords='axes fraction', fontsize=rcParams['font.size']-2, color='w') 
        
        majorLocator = MultipleLocator(30)
        majorFormatter = FormatStrFormatter('%d')
        ax2.xaxis.set_major_locator(majorLocator)
        ax2.yaxis.set_major_locator(majorLocator)
        
        ax2.set_ylabel('Distance (m)', fontsize = 9)
        ax2.set_xlabel('Distance (m)', fontsize = 9)
        
        ax2.set_ylim (250,100)
        ax2.set_xlim (100,250)
        
        # Set up the colour scheme
        Min_value = np.amin (Crossover[Crossover>Nodata_value])
        Search_space_mask = np.ma.masked_where(Search_space == 0, Search_space)

        # Plot the map
        Vmin = -0.1
        Vmax = 0.3
        Map_SS = ax2.imshow(Crossover, interpolation='None', cmap=plt.cm.gist_earth, vmin=Vmin, vmax=Vmax)
        Map_SS = ax2.imshow(Search_space_mask[:,0:181], interpolation='None', cmap=plt.cm.copper, vmin=Vmin, vmax=Vmax)
        
        # Make the colourbars
        ax12 = fig.add_axes([0.743, 0.115, 0.03, 0.295])
        ax22 = fig.add_axes([0.773, 0.115, 0.03, 0.295])

        X_scheme = plt.cm.copper; X_norm = mpl.colors.Normalize(vmin=Vmin, vmax=Vmax)
        SS_scheme = plt.cm.gist_earth; SS_norm = mpl.colors.Normalize(vmin=Vmin, vmax=Vmax)
        bounds = [-0.1, 0, 0.1, 0.2,0.3, 0.4, 0.5]
        cb1 = mpl.colorbar.ColorbarBase(ax12, cmap=X_scheme, norm=X_norm, ticks = [], orientation='vertical')
        cb2 = mpl.colorbar.ColorbarBase(ax22, cmap=SS_scheme, norm=SS_norm, ticks = bounds, orientation='vertical')
        cb2.set_label('X-value (m)', fontsize = 9)

    i=i+1



left  = 0.15  # the left side of the subplots of the figure
right = 0.2    # the right side of the subplots of the figure
bottom = 0.05   # the bottom of the subplots of the figure
top = 0.0      # the top of the subplots of the figure
wspace = 0.0   # the amount of width reserved for blank space between subplots
hspace = 0.0   # the amount of height reserved for white space between subplots
    
subplots_adjust(left=left, bottom=bottom, right=None, top=None, wspace=None, hspace=hspace)


plt.savefig('Output/Paper/0_Main_Fig2.png')"""





#---------------------------------------------------------------------------
# Figure 3 [2col]: This one shows the construction of scarps. It has a diagram showing how we proceed (a), and a slopes array with local max and scarp order (b)



fig=plt.figure(3, facecolor='White',figsize=[4.7,8])



"""matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8) 

ax1 = plt.subplot2grid((4,2),(0,0),colspan=2, rowspan=2,axisbg='white')
ax2 = plt.subplot2grid((4,2),(2,0),colspan=1, rowspan=1,axisbg='white')
ax3 = plt.subplot2grid((4,2),(2,1),colspan=1, rowspan=1,axisbg='white')
ax4 = plt.subplot2grid((4,1),(3,0),colspan=2, rowspan=1,axisbg='white')


ax1.set_ylabel('Distance (m)', fontsize = 9)
ax1.set_xlabel('Distance (m)', fontsize = 9)
ax1.annotate('a.', xy=(0.01,0.98), xycoords='axes fraction',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']-2)  
bbox_props = dict(boxstyle="square,pad=0.1", alpha = 0, fc="white", ec="b", lw=0)
ax2.text(365, 182, 'b.', ha="left", va="top", rotation=0, size=9, bbox=bbox_props, color='white')
ax3.text(345.5, 211.5, 'c.', ha="left", va="top", rotation=0, size=9, bbox=bbox_props, color='white')


ax1.xaxis.tick_top()
ax1.xaxis.set_label_position('top') 


i=0
for gauge in Gauges:
    
    # Choose your site
    if i == 4:
        # Load the data
        DEM, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array ("Input/Topography/%s/%s_DEM_WFILT.bil" % (gauge,gauge), gauge)
        Slope, post_Slope, envidata_Slope =  ENVI_raster_binary_to_2d_array ("Input/Topography/%s/%s_slope.bil" % (gauge,gauge), gauge)
        Search_space, Crossover, bins, hist, Inflexion_point = define_search_space (DEM, Slope, Nodata_value)
        Search_space[Search_space==1] = Slope[Search_space==1]
        Scarps, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array ("Output/%s/%s_Scarps.bil" % (gauge,gauge), gauge)

        SS_masked = np.ma.masked_where(Search_space == 0, Search_space)
        Scarps_masked = np.ma.masked_where(Scarps < 1, Scarps)
        
        
        Map_DEM = ax1.imshow(DEM, interpolation='None', cmap=plt.cm.Greys, vmin=0, vmax=6.5)
        Map_SS = ax1.imshow(SS_masked, interpolation='None', cmap=plt.cm.copper, vmin=0, vmax=0.4)
        Map_Scarps = ax1.imshow(Scarps_masked, interpolation='None', cmap=plt.cm.jet, vmin=0, vmax=12)
        
        Map_DEM = ax2.imshow(DEM, interpolation='None', cmap=plt.cm.Greys, vmin=0, vmax=6.5)
        Map_SS = ax2.imshow(SS_masked, interpolation='None', cmap=plt.cm.copper, vmin=0, vmax=0.4)
        Map_Scarps = ax2.imshow(Scarps_masked, interpolation='None', cmap=plt.cm.jet, vmin=0, vmax=12)
        
        Scarps_masked = np.ma.masked_where(Scarps_masked >= 2, Scarps_masked)
        
        Map_DEM = ax3.imshow(DEM, interpolation='None', cmap=plt.cm.Greys, vmin=0, vmax=6.5)
        Map_SS = ax3.imshow(SS_masked, interpolation='None', cmap=plt.cm.copper, vmin=0, vmax=0.4)
        Map_Scarps = ax3.imshow(Scarps_masked, interpolation='None', cmap=plt.cm.jet, vmin=0, vmax=12)
    
    
        
    
    
    i=i+1

ax1.add_patch(Rectangle((310, 180), 60, 60, fill=None, color = 'r', alpha=1))
ax2.add_patch(Rectangle((345, 211), 13, 13, fill=None, color = 'r', alpha=1))
    
ax2.set_ylim (240,180)
ax2.set_xlim (310,370)

ax3.set_ylim (224,211)
ax3.set_xlim (345,358)


    
# Make the colourbars
#ax2 = fig.add_axes([0.575, 0.25, 0.3, 0.02])
#TF_scheme = plt.cm.jet; TF_norm = mpl.colors.Normalize(vmin=0, vmax=10)
#cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=TF_scheme, norm=TF_norm, orientation='horizontal')
    

    
X = [352,352,352,351,351,352,352,352]
Y = [213,214,215,216,217,218,219,220]
line1, = ax3.plot(X,Y,'-r')
add_arrow_to_line2D(ax3, line1, arrow_locs=np.linspace(0., 1., 30), arrowstyle='->', arrowsize=0.8)    
X = [352,351,350]
Y = [213,212,211]
line1, = ax3.plot(X,Y,'-w')
add_arrow_to_line2D(ax3, line1, arrow_locs=np.linspace(0., 1., 30), arrowstyle='->', arrowsize=0.8)
      
X = [354,354,353,352,351,350,349,348,347,346,347,348,349,350,351]
Y = [219,220,221,220,221,221,221,221,221,220,219,219,219,219,219]
line2, = ax3.plot(X,Y,'-w')
add_arrow_to_line2D(ax3, line2, arrow_locs=np.linspace(0., 1., 30), arrowstyle='->', arrowsize=0.8)       
X = [354,354,355,356,356,356]
Y = [219,218,217,217,218,218]
line2, = ax3.plot(X,Y,'-r')
add_arrow_to_line2D(ax3, line2, arrow_locs=np.linspace(0., 1., 30), arrowstyle='->', arrowsize=0.8)
               
X = [358,357,356,356,356,356]
Y = [222,222,221,220,219,218]
line3, = ax3.plot(X,Y,'-w')
add_arrow_to_line2D(ax3, line3, arrow_locs=np.linspace(0., 1., 30), arrowstyle='->', arrowsize=0.8)       
      
X = [357,357,357,356]
Y = [214,215,216,217]
line4, = ax3.plot(X,Y,'-w')
add_arrow_to_line2D(ax3, line4, arrow_locs=np.linspace(0., 1., 30), arrowstyle='->', arrowsize=0.8)
X = [357,358,359]
Y = [214,213,212]
line4, = ax3.plot(X,Y,'-r')
add_arrow_to_line2D(ax3, line4, arrow_locs=np.linspace(0., 1., 30), arrowstyle='->', arrowsize=0.8)
 
     
Map_Scarps = ax3.imshow(Scarps_masked*6, interpolation='None', cmap = plt.cm.jet, vmin=0, vmax=12)   
  
    """




plt.savefig('Output/Paper/0_Main_Fig3.png')
    

#---------------------------------------------------------------------------
# Figure 4: This one shows a diagram of the process (a) and an array with filled-in platform (b) and cleaned-platform(c)




#---------------------------------------------------------------------------
# Figure 6 [2col]: This one showcases the results in a combined way
fig=plt.figure(5, facecolor='White',figsize=[4.7,7.5])

# Set up the fonts and stuff
matplotlib.rc('xtick', labelsize=9) 
matplotlib.rc('ytick', labelsize=9)

gs = gridspec.GridSpec(3, 2)
gs.update(bottom = 0.15,hspace=0.03,wspace = 0.03)

# Set up annotations
Annotations = ['a.','b.','c.','d.','e.','f.','g.']

i = 0
for gauge in Gauges: 
    # Set up the plot space 
    if (-1)**i > 0:
        Plot_col = 0
    else:
        Plot_col = 1
    Plot_row = int(np.floor(float(i)/2))

    ax1 = plt.subplot(gs[Plot_row,Plot_col])
    ax1.annotate(Annotations[i], xy=(0.05,0.90), xycoords='axes fraction', fontsize=rcParams['font.size']-2, color='w') 
    
    ax1.set_xticklabels([])
    ax1.set_yticklabels([])
    ax1.tick_params(axis='x', colors='white')
    ax1.tick_params(axis='y', colors='white')
    
    
    # Load the relevant data
    DEM, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array ("Input/Topography/%s/%s_DEM_WFILT.bil" % (gauge,gauge), gauge)
    Platform, post_Slope, envidata_Slope =  ENVI_raster_binary_to_2d_array ("Output/%s/%s_Marsh.bil" % (gauge,gauge), gauge)
    Confusion_matrix, post_Curvature, envidata_Curvature =  ENVI_raster_binary_to_2d_array ("Output/%s/%s_Confusion_DEM.bil" % (gauge,gauge), gauge)

    # Set up the colour scheme
    Min_value = np.amin (DEM[DEM>Nodata_value])
    DEM_mask = np.ma.masked_where(Platform == 0, DEM)
    Platform[Platform > 0] = DEM [Platform > 0]

    Confusion_matrix1 = np.ma.masked_where(Confusion_matrix !=-1, Confusion_matrix)
    Confusion_matrix1 [Confusion_matrix1 == -1] = DEM [Confusion_matrix1 == -1]
    Confusion_matrix2 = np.ma.masked_where(Confusion_matrix !=-2, Confusion_matrix) 
    Confusion_matrix2 [Confusion_matrix2 == -2] = DEM [Confusion_matrix2 == -2]

    # Plot the things
    Map_TF = ax1.imshow(DEM, interpolation='None', cmap=plt.cm.Greys, vmin=-5, vmax=5)
    Map_Marsh = ax1.imshow(DEM_mask, interpolation='None', cmap=plt.cm.gist_earth, vmin=0, vmax=8)
    Map_FP = ax1.imshow(Confusion_matrix1, interpolation='None', cmap=plt.cm.Oranges, vmin=0, vmax=5)
    Map_FN = ax1.imshow(Confusion_matrix2, interpolation='None', cmap=plt.cm.Purples, vmin=0, vmax=5)

    
    # Only display a hind-picked (square) area
    if i == 0:
        xmin = 0; xmax = 325
        ymin = 325; ymax = ymin - (xmax-xmin)
    elif i == 1:
        xmin = 0; xmax = 225
        ymin = 225; ymax = ymin - (xmax-xmin)
    elif i == 2:
        xmin = 140; xmax = 490
        ymin = 380; ymax = ymin - (xmax-xmin)
    elif i == 3:
        xmin = 100; xmax = 550
        ymin = 550; ymax = ymin - (xmax-xmin)
    elif i == 4:
        xmin = 0; xmax = 420
        ymin = 420; ymax = ymin - (xmax-xmin)
    elif i == 5:
        xmin = 0; xmax = 210
        ymin = 210; ymax = ymin - (xmax-xmin)
    
    ax1.set_xlim(xmin,xmax)
    ax1.set_ylim(ymin,ymax)
          
    i=i+1
        
        
        
# Make the colourbars
ax2 = fig.add_axes([0.175, 0.12, 0.3, 0.02])
ax3 = fig.add_axes([0.175, 0.10, 0.3, 0.02])
ax4 = fig.add_axes([0.575, 0.12, 0.3, 0.02])
ax5 = fig.add_axes([0.575, 0.10, 0.3, 0.02])

TF_scheme = plt.cm.Greys; TF_norm = mpl.colors.Normalize(vmin=-4, vmax=8)
Marsh_scheme = plt.cm.gist_earth; Marsh_norm = mpl.colors.Normalize(vmin=-4, vmax=8)
FP_scheme = plt.cm.Oranges; FP_norm = mpl.colors.Normalize(vmin=-0, vmax=5)
FN_scheme = plt.cm.Purples; FN_norm = mpl.colors.Normalize(vmin=-0, vmax=5)

bounds_1 = [-4,-2,0,2,4,6,8]
bounds_2 = [0,1,2,3,4,5]

cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=TF_scheme, norm=TF_norm, ticks = [], orientation='horizontal')
cb2 = mpl.colorbar.ColorbarBase(ax3, cmap=Marsh_scheme, norm=Marsh_norm, ticks = bounds_1, orientation='horizontal')
cb3 = mpl.colorbar.ColorbarBase(ax4, cmap=FP_scheme, norm=FP_norm, ticks = [], orientation='horizontal')
cb4 = mpl.colorbar.ColorbarBase(ax5, cmap=FN_scheme, norm=FN_norm, ticks = bounds_2, orientation='horizontal')
cb2.set_label('Elevation (m)', fontsize = 9)
cb4.set_label('Elevation (m)', fontsize = 9)
      
ax2.annotate('True negatives', xy=(0.05,0.85), xycoords='axes fraction',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']-3)  
ax3.annotate('True positives', xy=(0.35,0.85), xycoords='axes fraction',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']-3)  
ax4.annotate('False positives', xy=(0.05,0.85), xycoords='axes fraction',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']-3)  
ax5.annotate('False negatives', xy=(0.05,0.85), xycoords='axes fraction',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']-3)  
    

    


plt.savefig('Output/Paper/0_Main_Fig5_noresample.png')






# TICKS ARE EVERY 50m



#---------------------------------------------------------------------------
# Figure 6 [1col]: Performance results for different resolutions
fig=plt.figure(6, facecolor='White',figsize=[4.7,4])

# Set up the fonts and stuff
matplotlib.rc('xtick', labelsize=9) 
matplotlib.rc('ytick', labelsize=9)

ax1 = plt.subplot2grid((3,1),(0,0),colspan=1, rowspan=1)
ax3 = plt.subplot2grid((3,1),(1,0),colspan=1, rowspan=1)
ax4 = plt.subplot2grid((3,1),(2,0),colspan=1, rowspan=1)

ax1.set_ylim(0.36,1.05)
ax3.set_ylim(0.36,1.05)
ax4.set_ylim(0.36,1.05)
ax1.set_xlim(0.5,11)
ax3.set_xlim(0.5,11)
ax4.set_xlim(0.5,11)



ax1.annotate('a.', xy=(0.95,0.85), xycoords='axes fraction', fontsize=rcParams['font.size']-2, color='k') 
ax3.annotate('b.', xy=(0.95,0.85), xycoords='axes fraction', fontsize=rcParams['font.size']-2, color='k') 
ax4.annotate('c.', xy=(0.95,0.85), xycoords='axes fraction', fontsize=rcParams['font.size']-2, color='k') 



ax1.set_ylabel('Accuracy', fontsize = 9)
ax3.set_ylabel('Reliability', fontsize = 9)
ax4.set_ylabel('Sensitivity', fontsize = 9)

ax4.set_xlabel('x,y resolution (m)', fontsize = 9)

ax1.set_xticklabels([])
ax3.set_xticklabels([])

width = 0.07

i = 0
for gauge in Gauges:
    Accuracy = np.zeros(len(Resolutions),dtype=np.float)
    Reliability = np.zeros(len(Resolutions),dtype=np.float)
    Sensitivity = np.zeros(len(Resolutions),dtype=np.float)
      
    j=0
    for res in Resolutions:
        with open ("Output/%s/%s_%s_Metrix.pkl" % (gauge,gauge,res), 'rb') as input_file:
            Metrix = cPickle.load(input_file)
            Accuracy [j]= Metrix[0]
            Reliability [j]= Metrix[1]
            Sensitivity [j]= Metrix[2]

        j = j+1
    
    with open ("Output/%s/%s_Metrix.pkl" % (gauge,gauge), 'rb') as input_file:
        Metrix = cPickle.load(input_file)
        Accuracy [0]= Metrix[0]
        Reliability [0]= Metrix[1]
        Sensitivity [0]= Metrix[2]
    
    ax1.bar(Resolutions_num + i*width, Accuracy, width, color=plt.cm.winter(0.1*Metrix_gauges[i,0]), linewidth = 0)    
    ax3.scatter (Resolutions_num, Reliability, color=plt.cm.winter(0.1*Metrix_gauges[i,0]),edgecolors='none')
    ax4.scatter (Resolutions_num, Sensitivity, facecolors='none', edgecolors=plt.cm.winter(0.1*Metrix_gauges[i,0]))
     
    i = i+1

left  = 0.15  # the left side of the subplots of the figure
right = 0.2    # the right side of the subplots of the figure
bottom = 0.10   # the bottom of the subplots of the figure
top = 0.0      # the top of the subplots of the figure
wspace = 0.0   # the amount of width reserved for blank space between subplots
hspace = 0.02   # the amount of height reserved for white space between subplots
    
subplots_adjust(left=left, bottom=bottom, right=None, top=None, wspace=None, hspace=hspace)
   
    



plt.savefig('Output/Paper/0_Main_Fig6_noresample.png')



STOP

    



#---------------------------------------------------------------------------
# Figure 7: This one shows the different performances classified by tidal range (a), type of coast (b)
fig=plt.figure(7, facecolor='White',figsize=[4.7,5])
matplotlib.rc('xtick', labelsize=9) 

# First map the map
ax1 = plt.subplot2grid((7,2),(0,0),colspan=1, rowspan=7, axisbg='white')

# create Basemap
m = Basemap(llcrnrlon=-6.5,llcrnrlat=49.5,urcrnrlon=3.5,urcrnrlat=59, resolution='i', projection='cass', lon_0=-4.36, lat_0=54.7)
m.drawcoastlines()
m.drawparallels(np.arange(-40,61.,2.), labels=[1,0,0,0], fontsize = 9)
m.drawmeridians(np.arange(-20.,21.,2.), labels=[0,0,1,0], fontsize = 9)

# Plot the points on it
lats = [50.6, 51.8, 52.9, 51.5, 54.85, 54.1, 51.15]
lons = [-1.9, 1.13, 0.8, 0.5, -3.3, -2.8, -3.1]
Metrix_gauges = np.zeros((len(Gauges),5), dtype=np.float)

#this is the annotation of the plot number
ax1.annotate('a.', xy=(0.90,0.98), xycoords='axes fraction',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']-1) 


# Load the tidal data
i=0
for gauge in Gauges:
    Metric1_tide, Metric2_tide, Metric3_tide, Subsample = Open_tide_stats ("Input/Tide/%s/%s_" % (gauge,gauge), gauge)
    Metrix_gauges[i,0] = np.mean (Metric2_tide[3])-np.mean (Metric2_tide[0])
    i = i+1

# Plot the points
x, y = m(lons,lats)
Scatt = ax1.scatter(x,y, s = 50, color=plt.cm.winter(0.1*Metrix_gauges[:,0]), alpha = 0.9, linewidth = 5)

# Make a colourbar
ax2 = fig.add_axes([0.125, 0.10, 0.352, 0.02])
scheme = plt.cm.winter; norm = mpl.colors.Normalize(vmin=0, vmax=12)
bounds = [0,2,4,6,8,10,12]
cb = mpl.colorbar.ColorbarBase(ax2, cmap=scheme, norm=norm, ticks = bounds, orientation='horizontal')
cb.set_label('Spring tidal range (m)', fontsize = 9)  
    






    
# Then make the other plot
matplotlib.rc('xtick', labelsize=9)
matplotlib.rc('ytick', labelsize=9)
Annotations = ['b.1.','b.2.','b.3.','b.4.','b.5.','b.6.','b.7.']

i=0
for gauge in Gauges:
    ax3 = plt.subplot2grid((7,2),(i,1),colspan=1, rowspan=1,axisbg='white')
    
    ax3.set_xlim (0,8)
    ax3.set_ylim (0, 0.1)
    plt.yticks([0,0.04, 0.08], fontsize = 9)
    if i == 6:
        ax3.set_xlabel('Elevation (m)', fontsize = 9)
    else:
        ax3.set_xticklabels([])     
    if i == 3:
        ax3.set_ylabel('Frequency', fontsize = 9)
    ax3.yaxis.tick_right()


    # Load the data
    DEM, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array ("Input/Topography/%s/%s_DEM_WFILT.bil" % (gauge,gauge), gauge)
    
    Platform, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array ("Output/%s/%s_Marsh.bil" % (gauge,gauge), gauge)
    Platform [Platform>0] = DEM [Platform>0]
    
    
    Metric1_tide, Metric2_tide, Metric3_tide, Subsample = Open_tide_stats ("Input/Tide/%s/%s_" % (gauge,gauge), gauge)
    Metrix_gauges[i,0] = np.mean (Metric2_tide[3])-np.mean (Metric2_tide[0])
    
    
    

    bins, hist = Distribution (DEM, Nodata_value)
    Elevation_range = np.arange(0, max(bins), 0.1)
    
    bins_M, hist_M = Distribution (Platform, Nodata_value)
    Elevation_range_M = np.arange(0, max(bins_M), 0.1)
    
    
    ax3.fill_between(Elevation_range, 0, 1, where=Elevation_range<=max(Metric2_tide[3]), facecolor='blue', alpha = 0.2, linewidth = 0)
    ax3.fill_between(Elevation_range, 0, 1, where=Elevation_range<=np.mean(Metric2_tide[3]), facecolor='blue', alpha = 0.2, linewidth = 0)
    ax3.fill_between(Elevation_range, 0, 1, where=Elevation_range<=np.mean(Metric2_tide[2]), facecolor='blue', alpha = 0.2, linewidth = 0)
    
    ax3.plot( bins, hist, color=plt.cm.winter(0.1*Metrix_gauges[i,0]))
    ax3.plot( bins_M, hist_M, 'r') #color=plt.cm.winter(0.1*Metrix_gauges[i,0]))

    
    ax3.annotate(Annotations[i], xy=(0.80,0.925), xycoords='axes fraction',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']-2) 


    i=i+1



plt.savefig('Output/Paper/0_Fig7.png')



STOP


#---------------------------------------------------------------------------
# Figure 8: This one shows evolution of performance for degraded resolution. Look at relief in the thing too, or position within the tidal frame. Still in design






for i in [0,1]:
    













    
    
    
    
    
    #Paper Plots
    fig=plt.figure(11, facecolor='White',figsize=[30,15])
    ax = plt.subplot2grid((1,2),(0,0),colspan=1, rowspan=2)
    #ax.set_title('Scarp relief = %g' % (np.amax(Scarp_DEM)-np.amin(Scarp_DEM)), fontsize = 22)
    #ax.set_title('Slope (from polynomial fit)', fontsize = 22)
    ax.set_title('Platform', fontsize = 22)
    ax.set_xlabel('X-distance from origin (m)', fontsize = 18)
    ax.set_ylabel('Y-distance from origin (m)', fontsize = 18)
    

    #Map = ax.plot(Scarps_bins, Scarps_hist)
    
    #Map = ax.imshow(Rel_relief, interpolation='None', cmap=palette_Slope, vmin=0, vmax=1)#, norm=colors.Normalize(vmin=Zmin, vmax=Zmax))
    #Map = ax.imshow(Crossover, interpolation='None', cmap=palette_Slope, vmin=0, vmax=1)#, norm=colors.Normalize(vmin=Zmin, vmax=Zmax))
    #Map = ax.imshow(Search_space, interpolation='None', cmap=palette_Slope, vmin=0, vmax=1)#, norm=colors.Normalize(vmin=Zmin, vmax=Zmax))

    Map = ax.imshow(Platform, interpolation='None', vmin = 0)#, cmap=palette_Rel_relief)#,  vmin=0, vmax=1)#, norm=colors.Normalize(vmin=Smin, vmax=Smax))
    cbar = fig.colorbar(Map, extend='both', shrink=0.95, ax=ax)
    cbar.set_label('Platform elevation (m)', fontsize = 20)
    #cbar.set_label('Relative relief (m)', fontsize = 20)

    
    
    
    
    
    
    
    #This is where you define the cutoff spot!

    Platform_bins, Platform_hist = Distribution(Platform,0)

    #1. Find the highest and biggest local maximum of frequency distribution
    for j in range(1,len(Platform_hist)-1):
        if Platform_hist[j]>0.9*max(Platform_hist) and Platform_hist[j]>Platform_hist[j-1] and Platform_hist[j]>Platform_hist[j+1]:
            Index  = j

    #2. Now run a loop from there toward lower elevations.
    Counter = 0
    for j in range(Index,0,-1):
        # See if you cross the mean value. Count for how many indices you are under.
        if Platform_hist[j] < mean(Platform_hist):
            Counter = Counter + 1
        # Reset the counter value if you go above average again
        else:
            Counter = 0 
            
        #If you stay long enough under (10 is arbitrary for now), initiate cutoff
        if Counter > 10:
            Cutoff = j
            break
            

       
    
    
    
            
    ax = plt.subplot2grid((1,2),(0,1),colspan=1, rowspan=2)
    ax.set_title('Detected scarps', fontsize = 22)
    #ax.set_xlabel('X-distance from origin (m)', fontsize = 18)
    ax.set_xlabel('Elevation (m)', fontsize = 18)
    ax.set_ylabel('Y-distance from origin (m)', fontsize = 18)

    Map = ax.plot(Platform_bins, Platform_hist)
    Scatt = ax.scatter(Platform_bins[Index], Platform_hist[Index], c = 'red', alpha = 0.5)
    ax.fill_between(Platform_bins, 0, mean(Platform_hist), alpha = 0.5)
    
    #plt.axvline(x=Platform_bins[Cutoff], ymin=0, linewidth=0.1)
    
    #Map = ax.imshow(Scarps, interpolation='None', cmap=plt.cm.gist_heat)
    #Map = ax.imshow(Rel_slope, interpolation='None', cmap=palette_Rel_slope,  vmin=0, vmax=1)#, norm=colors.Normalize(vmin=Zmin, vmax=Zmax))
    #Map = ax.imshow(Scarps, interpolation='None', cmap=palette_Scarps)#, norm=colors.Normalize(vmin=Smin, vmax=Smax))
    #cbar = fig.colorbar(Map, extend='both', shrink=0.95, ax=ax)
    #cbar.set_label('Scarp slope (m/m)', fontsize = 20)

    plt.savefig('Output/Paper/%s_Paper_Scarps7.png' % (gauge))



    
    
    
    fig=plt.figure(13, facecolor='White',figsize=[60,30])
    ax = plt.subplot2grid((1,2),(0,0),colspan=1, rowspan=2)
    ax.set_title('Slope (from polynomial fit)', fontsize = 22)
    ax.set_xlabel('X-distance from origin (m)', fontsize = 18)
    ax.set_ylabel('Y-distance from origin (m)', fontsize = 18)
    Map = ax.imshow(Platform, interpolation='None', cmap=plt.cm.gist_heat)#, vmin=0, vmax=3)
    cbar = fig.colorbar(Map, extend='both', shrink=0.95, ax=ax)
    cbar.set_label('XXXX', fontsize = 20)

    ax = plt.subplot2grid((1,2),(0,1),colspan=1, rowspan=2)
    ax.set_title('Detected scarps', fontsize = 22)
    ax.set_xlabel('X-distance from origin (m)', fontsize = 18)
    ax.set_ylabel('Y-distance from origin (m)', fontsize = 18)
    Map = ax.imshow(Scarps, interpolation='None', cmap=palette_Scarps)#, vmin= 0, vmax=1)
    cbar = fig.colorbar(Map, extend='both', shrink=0.95, ax=ax)
    cbar.set_label('Scarp slope (m/m)', fontsize = 20)

    #plt.savefig('Output/Paper/%s_Paper_fig15.png' % (gauge))


    
    fig=plt.figure(14, facecolor='White',figsize=[15,15])
    ax = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=2)
    ax.set_title('Confusion map', fontsize = 22)
    ax.set_xlabel('X-distance from origin (m)', fontsize = 18)
    ax.set_ylabel('Y-distance from origin (m)', fontsize = 18)
    Map = ax.imshow(Confusion_matrix, interpolation='None', cmap=palette_Confusion, norm=colors.Normalize(vmin=-2, vmax=2))
    cbar = fig.colorbar(Map, ticks=[-2,-1,1,2], shrink=0.95, ax=ax)
    #cbar.ax.set_yticklabels(['False Negative', 'False Positive', 'True Positive', 'True Negative'], rotation = 90, fontsize = 16)
    cbar.ax.set_yticklabels(['FN', 'FP', 'TP', 'TN'], rotation = 90, fontsize = 16)
    cbar.set_label('Confusion value', fontsize = 20)

    #rect = [0.55,0.25,0.5,0.25]
    #axbis = add_subplot_axes(ax,rect)
    #TP=Performance[0]; TN=Performance[1]; FP=Performance[2]; FN=Performance[3]
    #axbis.set_title("Correct marsh points: %d " % (100*(TP+TN)/(TP+TN+FP+FN)) + "%", fontsize = 18)
    #sizes = Performance
    #colormap = [plt.cm.RdYlGn(200), plt.cm.RdYlGn(256), plt.cm.RdYlGn(64), plt.cm.RdYlGn(0)]
    #axbis.pie(sizes, autopct='%1.0f%%',shadow=False, startangle=90, colors = colormap)
    #from matplotlib import font_manager as fm
    #proptease = fm.FontProperties()
    #proptease.set_size('xx-small')
    #axbis.axis('equal')

    plt.savefig('Output/Paper/%s_Confusion_nopie2.png' % (gauge))

    
    
    
    i = i + 1

#-----------------------------------------------------------------------------------------
#Global figure for performance
Gauge_label= ["Shell Bay","Stour Estuary","Chalksock Lake","Campfield Marsh", "Dee Estuary", "Morecambe Bay", "Steart Point", "Severn Estuary"]


Gauge_label=["BOU", "FEL", "CRO", "SHE", "WOR", "HEY", "HIN"]

fig=plt.figure(15, facecolor='White',figsize=[20,8])
ax = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=2)
ax.set_xlabel('Spring tidal range (m)', fontsize = 18)
ax.set_ylabel('Performance values', fontsize = 18)
ax.grid(True)

for i in range(len(Gauge_label)):
    
    print i
    
    #ax.axvline(Metrix_gauges[i,0], ymin=75, ymax=0.95, color='black', lw=5, alpha=0.6)
    ax.annotate(Gauge_label[i], xy=(Metrix_gauges[i,0]-0.024, 0.62), xycoords='data',
        horizontalalignment='left', verticalalignment='bottom', fontsize=rcParams['font.size']-0.5, color='black', rotation = 90)

yerr_lower=np.zeros(len(Metrix_gauges[:,1]), dtype=np.float)
yerr_upper=np.zeros(len(Metrix_gauges[:,1]), dtype=np.float)
for i in range(len(Metrix_gauges[:,1])):
    if Metrix_gauges[i,2] > Metrix_gauges[i,3]:
        yerr_lower[i] =  Metrix_gauges[i,4] - Metrix_gauges[i,3]
        yerr_upper[i] =  Metrix_gauges[i,2] - Metrix_gauges[i,4]
    else:
        yerr_lower[i] =  Metrix_gauges[i,4] - Metrix_gauges[i,2]
        yerr_upper[i] =  Metrix_gauges[i,3] - Metrix_gauges[i,4]



errorbar(Metrix_gauges[:,0], Metrix_gauges[:,4], yerr=[yerr_lower, yerr_upper], fmt='o', ecolor='k', capthick=2, elinewidth = 5, alpha=0.6)


performance = Metrix_gauges[:,1]; performance_colour = performance**3
ax.bar(Metrix_gauges[:,0]-0.09, Metrix_gauges[:,1], 0.18, alpha=0.6, color=plt.cm.RdYlGn(performance_colour),label='Accuracy')
line, = ax.plot(Metrix_gauges[:,0], Metrix_gauges[:,2], 'oy', label='Reliability')
line, = ax.plot(Metrix_gauges[:,0], Metrix_gauges[:,3], 'or', label='Sensitivity')
line, = ax.plot(Metrix_gauges[:,0], Metrix_gauges[:,4], 'ow', label='F1')

ax.set_xlim(2,13.5)
ax.set_ylim(0.6,1.)
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), fancybox=True, shadow=False, ncol=4)



plt.savefig('Output/Paper/Global_Perf.png') 
    
    
    
    
    
    
    
    
    
    
    

for i in [0,1]:
    STOP

   
 
    
    
    

"""







    #Platform[Platform > 0] = DEM[Platform > 0]

    Search_x = Crossover
    Z_data = Search_x.ravel(); Z_data = Z_data[Z_data>0]

    step = (max(Z_data) - min(Z_data)) / 50
    #step = 0.05
    Z_value = np.arange(min(Z_data), max(Z_data), step)

    hist_gauge, bins_gauge = np.histogram (Z_data, Z_value, density=True)
    hist_gauge=hist_gauge/sum(hist_gauge); bins_gauge=bins_gauge[:-1]
    hist.append(hist_gauge)
    bins.append(bins_gauge)


    #Derivative
    hist_der = np.zeros(len(hist_gauge), dtype = np.float)
    hist_curv = np.zeros(len(hist_gauge), dtype = np.float)
    for j in range(1, len(hist_gauge), 1):
        hist_der[j] = (hist_gauge[j]-hist_gauge[j-1])/step

    for j in range(1, len(hist_gauge)-1, 1):
        if hist_der[j] < -1 and hist_der[j+1] >= -1:
            Inflexion_point[i] = bins_gauge[j]

    #Search_s = Scarps
    #S_data = Search_s.ravel(); S_data = S_data[S_data>0]; S_data = S_data[S_data<1]
    #step = (max(S_data) - min(S_data)) / 200
    #S_value = np.arange(min(S_data), max(S_data), step)

    #hist_gauge, bins_gauge = np.histogram (S_data, S_value, density=True); hist_gauge=hist_gauge/sum(hist_gauge); bins_gauge=bins_gauge[:-1]
    #histS.append(hist_gauge)
    #binsS.append(bins_gauge)


    i = i+1

cmap = mpl.cm.gist_heat
fig=plt.figure(12, facecolor='White',figsize=[20,10])
ax = plt.subplot2grid((1,2),(0,0),colspan=1, rowspan=1)
ax2 = plt.subplot2grid((1,2),(0,1),colspan=1, rowspan=1)
for i in range(len(Gauges)):
    ax.plot(bins[i], hist[i], color = cmap(i*40), alpha = 0.8)#, linewidth = 0)
    ax.axvline(Inflexion_point[i], color=cmap(i*40), lw=2, alpha=0.5)
    #ax2.plot(binsS[i], HIST[i], color = cmap(i*40))
    #ax2.axvline(Inflexion_point[i], color=cmap(i*40), lw=2, alpha=0.5)
    #ax2.plot(binsS[i], histS[i], color = cmap(i*40))

ax.set_xlim(xmax = 0.4)
ax.set_ylim(ymax = 0.4)
#ax2.set_xlim(xmax = 0.2)
#ax.set_ylim(ymax = 0.2)

#plt.savefig('Output/Paper/Crossover.png')






    #---------------------------------------------------------------------------
    #Poster Plots

    fig=plt.figure(10, facecolor='White',figsize=[15,15])
    ax = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=2)
    ax.set_title('Digital Elevation Model (from LiDAR)', fontsize = 22)
    ax.set_xlabel('X-distance from origin (m)', fontsize = 18)
    ax.set_ylabel('Y-distance from origin (m)', fontsize = 18)
    Map = ax.imshow(DEM, interpolation='None', cmap=palette_DEM, norm=colors.Normalize(vmin=Zmin, vmax=Zmax))
    cbar = fig.colorbar(Map, extend='both', shrink=0.95, ax=ax)
    cbar.set_label('Elevation (m)', fontsize = 20)

    plt.savefig('Output/Poster/%s_DEM.png' % (gauge))



    fig=plt.figure(11, facecolor='White',figsize=[15,15])
    ax = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=2)
    ax.set_title('Slope (from polynomial fit)', fontsize = 22)
    ax.set_xlabel('X-distance from origin (m)', fontsize = 18)
    ax.set_ylabel('Y-distance from origin (m)', fontsize = 18)
    Map = ax.imshow(Slope, interpolation='None', cmap=palette_Slope)#, norm=colors.Normalize(vmin=Zmin, vmax=Zmax))
    cbar = fig.colorbar(Map, extend='both', shrink=0.95, ax=ax)
    cbar.set_label('Slope (m/m)', fontsize = 20)

    plt.savefig('Output/Poster/%s_Slope.png' % (gauge))




    fig=plt.figure(12, facecolor='White',figsize=[15,15])
    ax = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=2)
    ax.set_title('Detected scarps', fontsize = 22)
    ax.set_xlabel('X-distance from origin (m)', fontsize = 18)
    ax.set_ylabel('Y-distance from origin (m)', fontsize = 18)
    #Map = ax.imshow(Scarps, interpolation='None', cmap=plt.cm.gist_heat)
    Map = ax.imshow(Scarps, interpolation='None', cmap=palette_Scarps, norm=colors.Normalize(vmin=Smin, vmax=Smax))
    cbar = fig.colorbar(Map, extend='both', shrink=0.95, ax=ax)
    cbar.set_label('Scarp slope (m/m)', fontsize = 20)

    plt.savefig('Output/Poster/%s_Scarps.png' % (gauge))




    fig=plt.figure(13, facecolor='White',figsize=[15,15])
    ax = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=2)
    ax.set_title('Detected platforms', fontsize = 22)
    ax.set_xlabel('X-distance from origin (m)', fontsize = 18)
    ax.set_ylabel('Y-distance from origin (m)', fontsize = 18)
    Map = ax.imshow(Platform, interpolation='None', cmap=palette_platform, norm=colors.Normalize(vmin=Zmin, vmax=Zmax))
    cbar = fig.colorbar(Map, extend='both', shrink=0.95, ax=ax)
    cbar.set_label('Elevation (m)', fontsize = 20)

    plt.savefig('Output/Poster/%s_Platform.png' % (gauge))




    fig=plt.figure(14, facecolor='White',figsize=[15,15])
    ax = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=2)
    ax.set_title('Confusion map', fontsize = 22)
    ax.set_xlabel('X-distance from origin (m)', fontsize = 18)
    ax.set_ylabel('Y-distance from origin (m)', fontsize = 18)
    Map = ax.imshow(Confusion_matrix, interpolation='None', cmap=palette_Confusion, norm=colors.Normalize(vmin=-2, vmax=2))
    cbar = fig.colorbar(Map, ticks=[-2,-1,1,2], shrink=0.95, ax=ax)
    #cbar.ax.set_yticklabels(['False Negative', 'False Positive', 'True Positive', 'True Negative'], rotation = 90, fontsize = 16)
    cbar.ax.set_yticklabels(['FN', 'FP', 'TP', 'TN'], rotation = 90, fontsize = 16)
    cbar.set_label('Confusion value', fontsize = 20)

    #rect = [0.55,0.25,0.5,0.25]
    #axbis = add_subplot_axes(ax,rect)
    #TP=Performance[0]; TN=Performance[1]; FP=Performance[2]; FN=Performance[3]
    #axbis.set_title("Correct marsh points: %d " % (100*(TP+TN)/(TP+TN+FP+FN)) + "%", fontsize = 18)
    #sizes = Performance
    #colormap = [plt.cm.RdYlGn(200), plt.cm.RdYlGn(256), plt.cm.RdYlGn(64), plt.cm.RdYlGn(0)]
    #axbis.pie(sizes, autopct='%1.0f%%',shadow=False, startangle=90, colors = colormap)
    #from matplotlib import font_manager as fm
    #proptease = fm.FontProperties()
    #proptease.set_size('xx-small')
    #axbis.axis('equal')

    plt.savefig('Output/Poster/%s_Confusion_nopie.png' % (gauge))






#-----------------------------------------------------------------------------------------
#Global figure for performance
Gauge_label= ["Shell Bay","Stour Estuary","Chalksock Lake","Campfield Marsh", "Dee Estuary", "Morecambe Bay", "Steart Point", "Severn Estuary"]


fig=plt.figure(15, facecolor='White',figsize=[20,8])
ax = plt.subplot2grid((1,1),(0,0),colspan=1, rowspan=2)
ax.set_xlabel('Spring tidal range (m)', fontsize = 18)
ax.set_ylabel('Performance values', fontsize = 18)
ax.grid(True)

for i in range(len(Gauge_label)):
    #ax.axvline(Metrix_gauges[i,0], ymin=75, ymax=0.95, color='black', lw=5, alpha=0.6)
    ax.annotate(Gauge_label[i], xy=(Metrix_gauges[i,0]-0.065, 0.62), xycoords='data',
        horizontalalignment='left', verticalalignment='bottom', fontsize=rcParams['font.size']-0.5, color='black', rotation = 90)

yerr_lower=np.zeros(len(Metrix_gauges[:,1]), dtype=np.float)
yerr_upper=np.zeros(len(Metrix_gauges[:,1]), dtype=np.float)
for i in range(len(Metrix_gauges[:,1])):
    if Metrix_gauges[i,2] > Metrix_gauges[i,3]:
        yerr_lower[i] =  Metrix_gauges[i,4] - Metrix_gauges[i,3]
        yerr_upper[i] =  Metrix_gauges[i,2] - Metrix_gauges[i,4]
    else:
        yerr_lower[i] =  Metrix_gauges[i,4] - Metrix_gauges[i,2]
        yerr_upper[i] =  Metrix_gauges[i,3] - Metrix_gauges[i,4]



errorbar(Metrix_gauges[:,0], Metrix_gauges[:,4], yerr=[yerr_lower, yerr_upper], fmt='o', ecolor='k', capthick=2, elinewidth = 5, alpha=0.6)


performance = Metrix_gauges[:,1]; performance_colour = performance**3
ax.bar(Metrix_gauges[:,0]-0.09, Metrix_gauges[:,1], 0.18, alpha=0.6, color=plt.cm.RdYlGn(performance_colour),label='Accuracy')
line, = ax.plot(Metrix_gauges[:,0], Metrix_gauges[:,2], 'oy', label='Reliability')
line, = ax.plot(Metrix_gauges[:,0], Metrix_gauges[:,3], 'or', label='Sensitivity')
line, = ax.plot(Metrix_gauges[:,0], Metrix_gauges[:,4], 'ow', label='F1')

ax.set_xlim(2,13.5)
ax.set_ylim(0.6,1.)
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
ax.legend(loc='upper center', bbox_to_anchor=(0.5, 1.1), fancybox=True, shadow=False, ncol=4)



plt.savefig('Output/Poster/Global_Perf.png')

for i in [0,1]:
    STOP





    #-------------------------------------------------------------------------------------------------------------------
    # Paper plots








fig = plt.figure(1, facecolor='White',figsize=[9,9])
ax5 = plt.subplot2grid((1,1),(0,0),colspan=1)
#ax6=ax5.twinx()


i=0

TR = Metrix_gauges[:,0]
Max_slope =  np.zeros(len(TR), dtype = np.float)
Med_slope =  np.zeros(len(TR), dtype = np.float)
Min_slope =  np.zeros(len(TR), dtype = np.float)

cmap = mpl.cm.rainbow

for gauge in Gauges:
    print " Loading Slopes"
    Slope, post_Slope, envidata_Slope =  ENVI_raster_binary_to_2d_array ("Input/Topography/%s/%s_slope.bil" % (gauge,gauge), gauge)
    #Z_data = Slope.ravel(); Z_data = Z_data[Z_data>0]
    #Z_value = np.arange(min(Z_data), max(Z_data), 0.005)

    #hist,bins = np.histogram (Z_data, Z_value, density=True); hist=hist/sum(hist); bins=bins[:-1]

    #Derivative?
    #hist_der = np.zeros(len(hist), dtype = np.float)
    #for j in range(1, len(hist), 1):
        #hist_der[j] = (hist[j]-hist[j-1])/0.005



    #ax6.plot(bins, hist_der, '--',color = cmap(TR[i]/max(TR)))
    #ax5.plot(bins, hist, color = cmap(TR[i]/max(TR)))


    Slope_filtered = Slope[Slope<1]
    Max_slope[i] = np.percentile(Slope_filtered, 95)
    Med_slope[i] = np.percentile(Slope_filtered, 50)
    Min_slope[i] = np.percentile(Slope_filtered, 10)

    i = i+1



ax5.scatter(TR, Max_slope, c='green')#, color = cmap(TR[i]/max(TR)))
ax5.scatter(TR, Med_slope)#, color = cmap(TR[i]/max(TR)))
ax5.scatter(TR, Min_slope, c='red')#, color = cmap(TR[i]/max(TR)))

ax5.set_xlim(xmin=0)
ax5.set_ylim(ymin=0, ymax=0.4)
ax5.grid(True)


#ax1 = fig.add_axes([0.91, 0.10, 0.02, 0.8])
#norm = mpl.colors.Normalize(vmin=min(TR), vmax=max(TR))
#cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='vertical')
#cb1.set_label('Tidal range (m)', fontsize = 10)

#plt.savefig('Output/Poster/00_Slope_distribution.png')

for i in [0,1]:
    STOP

"""



#----------------------------------------------------------------
# Figure: Displays map, Pie chart and R/S values
#fig=plt.figure(1, facecolor='White',figsize=[50,50])


# Platform Map
"""ax = plt.subplot2grid((4,2),(0,0),colspan=1, rowspan=2)
ax.set_title('a. Platform DEM')
Map = ax.imshow(Platform, interpolation='None', cmap=palette_platform, norm=colors.Normalize(vmin=Zmin, vmax=Zmax))
cbar = fig.colorbar(Map, extend='both', shrink=0.95, ax=ax)
cbar.set_label('elevation (m)')"""







# Confusion Map
"""ax = plt.subplot2grid((4,2),(2,0),colspan=1, rowspan=2)
ax.set_title('b. Confusion map')

Map = ax.imshow(Reference, interpolation = 'None', cmap = palette_platform, alpha = 0.8)
Map = ax.imshow(Confusion_matrix, interpolation = 'None', cmap = palette_Confusion, alpha = 0.8, norm=colors.Normalize(vmin=-2.0, vmax=2.0))
cbar = fig.colorbar(Map, ticks=[-2,-1,1,2], shrink=0.95, ax=ax)
cbar.ax.set_yticklabels(['FN', 'FP', 'TP', 'TN'])
cbar.set_label('Confusion value')"""


"""#Original DEM
ax = plt.subplot2grid((4,2),(0,0),colspan=1, rowspan=2)
ax.set_title('Original DEM')

DEM[DEM==Nodata_value]=0
Map = ax.imshow(DEM, interpolation='None', cmap=plt.cm.gist_earth)

cbar = fig.colorbar(Map, extend='both', shrink=0.95, ax=ax)
cbar.set_label('elevation (m)')



#Digitised platform
ax = plt.subplot2grid((4,2),(2,0),colspan=1, rowspan=2)
ax.set_title('Digitised Platform DEM')

Reference[Reference==1] = DEM[Reference==1]
Map = ax.imshow(Reference, interpolation='None', cmap=plt.cm.gist_earth)

cbar = fig.colorbar(Map, extend='both', shrink=0.95, ax=ax)
cbar.set_label('elevation (m)')




#Scarps
Scarps2 = np.copy(Scarps)
Scarps2[Scarps!=0] = Slope[Scarps!=0] * DEM[Scarps!=0]

ax = plt.subplot2grid((4,2),(0,1),colspan=1, rowspan=2)
ax.set_title('Scarps Slope')
Map = ax.imshow(Scarps2, interpolation='None', cmap=plt.cm.RdYlGn)
cbar = fig.colorbar(Map, extend='both', shrink=0.95, ax=ax)
cbar.set_label('Slope')"""


#Scarps2
"""Scarps[Scarps!=0] = Slope[Scarps!=0] * DEM[Scarps!=0] / (np.mean(Metric2_tide[3])-np.mean(Metric2_tide[0]))

ax = plt.subplot2grid((4,2),(2,1),colspan=1, rowspan=2)
ax.set_title('Scarps Elevation')
Map = ax.imshow(Scarps, interpolation='None', cmap=plt.cm.RdYlGn)
cbar = fig.colorbar(Map, extend='both', shrink=0.95, ax=ax)
cbar.set_label('Elevation')"""


# Platform Map

#Platform[Platform>0] = DEM[Platform>0]


"""ax = plt.subplot2grid((4,2),(2,1),colspan=1, rowspan=2)
ax.set_title('a. Platform DEM')
Map = ax.imshow(Platform, interpolation='None', cmap=plt.cm.gist_earth, norm=colors.Normalize(vmin=0.))
cbar = fig.colorbar(Map, extend='both', shrink=0.95, ax=ax)
cbar.set_label('elevation (m)')"""



#Search Space
"""ax = plt.subplot2grid((4,2),(2,1),colspan=1, rowspan=2)
ax.set_title('Search Space')
Map = ax.imshow(Search_space, interpolation = 'None', cmap = plt.cm.gist_earth, alpha = 0.8)
cbar = fig.colorbar(Map, shrink=0.95, ax=ax)
cbar.set_label('Slope')"""







# Pie
""" ax = plt.subplot2grid((4,2),(0,1),colspan=1, rowspan=2)
ax.set_title('c. Proportional confusion distribution')

labels = 'TP', 'TN', 'FP', 'FN'
sizes = Performance
colormap = [plt.cm.RdYlGn(200), plt.cm.RdYlGn(256), plt.cm.RdYlGn(64), plt.cm.RdYlGn(0)]
ax.pie(sizes, labels=labels, autopct='%1.0f%%',shadow=False, startangle=90, colors = colormap)
ax.axis('equal')


from matplotlib import font_manager as fm
proptease = fm.FontProperties()
proptease.set_size('xx-small')"""

# Metrix
"""ax = plt.subplot2grid((4,3),(2,1),colspan=1, rowspan=2)
ax.set_title('d. Metrics')

metrix = ('Accuracy', 'Precision', 'Sensitivity', 'F1')
y_pos = np.arange(len(metrix))
performance = Metrix
performance_colour = performance**3

ax.barh(y_pos, performance, align='center', color=plt.cm.RdYlGn(performance_colour))
ax.set_yticks(y_pos); ax.set_yticklabels(metrix)
ax.invert_yaxis(); ax.invert_xaxis()
ax.set_xlabel('non-dimensional value')
plt.yticks(y_pos, metrix, rotation=90, fontsize=11)
ax.set_xlim(0,1)
plt.margins(0.1)"""




#plt.savefig('Output/%s_Confusion.png' % (gauge))
#plt.savefig('Output/%s/%s_Confusion.png' % (gauge,gauge))




"""# Setting up the spectral analysis
print " Loading red"
sourcefile = "Input/Ortho_BOU.tif"
destinationfile = "Input/OUTPUT.bil"
os.system("gdal_translate -b 3 -of ENVI " + sourcefile + " " +  destinationfile)
Band1, post_Band1, envidata_Band1 =  ENVI_raster_binary_to_2d_array (destinationfile, 'Band1')

#print " Loading green"
#sourcefile = "Input/Ortho_HIN.tif"
#destinationfile = "Input/OUTPUT.bil"
#os.system("gdal_translate -b 2 -of ENVI " + sourcefile + " " +  destinationfile)
#Band2, post_Band2, envidata_Band2 =  ENVI_raster_binary_to_2d_array (destinationfile, 'Band2')

#print " Loading blue"
#sourcefile = "Input/Ortho_HIN.tif"
#destinationfile = "Input/OUTPUT.bil"
#os.system("gdal_translate -b 3 -of ENVI " + sourcefile + " " +  destinationfile)
#Band3, post_Band3, envidata_Band3 =  ENVI_raster_binary_to_2d_array (destinationfile, 'Band3')

print " Loading NIR"
sourcefile = "Input/Ortho_BOU.tif"
destinationfile = "Input/OUTPUT.bil"
os.system("gdal_translate -b 4 -of ENVI " + sourcefile + " " +  destinationfile)
Band4, post_Band4, envidata_Band4 =  ENVI_raster_binary_to_2d_array (destinationfile, 'Band4')


Red = np.copy(Band1)
NIR = np.copy(Band4)

#print 'Red'; print Red[0:4, 0:4]
#print 'NIR'; print NIR[0:4, 0:4]

NDVI_num =  NIR-Red
NDVI_denom =  NIR+Red
NDVI =  np.divide(NDVI_num.astype(float),NDVI_denom.astype(float))

#print 'Red2'; print Red[0:4, 0:4]
#print 'NIR2'; print NIR[0:4, 0:4]
#print 'NDVI top'; print NDVI_num[0:4, 0:4]
#print 'NDVI bottom'; print NDVI_denom[0:4, 0:4]
print 'NDVI'; print NDVI[0:4, 0:4]


#NDVI2 = (Band4-Band2)/(Band4+Band2)
#NDVI3 = (Band4-Band3)/(Band4+Band3)



fig=plt.figure(1, facecolor='White',figsize=[15,15])

ax = plt.subplot2grid((2,2),(0,0),colspan=1, rowspan=1)
ax.set_title('Band4', fontsize = 22)
ax.set_xlabel('X-distance from origin (m)', fontsize = 18)
ax.set_ylabel('Y-distance from origin (m)', fontsize = 18)
Map = ax.imshow(Band4, interpolation='None', cmap=plt.cm.Greys)
cbar = fig.colorbar(Map, extend='both', shrink=0.95, ax=ax)
cbar.set_label('value', fontsize = 20)

ax = plt.subplot2grid((2,2),(0,1),colspan=1, rowspan=1)
ax.set_title('NDVI', fontsize = 22)
ax.set_xlabel('X-distance from origin (m)', fontsize = 18)
ax.set_ylabel('Y-distance from origin (m)', fontsize = 18)
Map = ax.imshow(NDVI, interpolation='None', cmap=plt.cm.RdYlGn, vmin=-1, vmax=1)
cbar = fig.colorbar(Map, extend='both', shrink=0.95, ax=ax)
cbar.set_label('value', fontsize = 20)

ax = plt.subplot2grid((2,2),(1,0),colspan=1, rowspan=1)
ax.set_title('Band1', fontsize = 22)
ax.set_xlabel('X-distance from origin (m)', fontsize = 18)
ax.set_ylabel('Y-distance from origin (m)', fontsize = 18)
Map = ax.imshow(Band1, interpolation='None', cmap=plt.cm.Greys)
cbar = fig.colorbar(Map, extend='both', shrink=0.95, ax=ax)
cbar.set_label('value', fontsize = 20)

#ax = plt.subplot2grid((2,2),(1,1),colspan=1, rowspan=1)
#ax.set_title('NDVI3', fontsize = 22)
#ax.set_xlabel('X-distance from origin (m)', fontsize = 18)
#ax.set_ylabel('Y-distance from origin (m)', fontsize = 18)
#Map = ax.imshow(NDVI3, interpolation='None', cmap=plt.cm.gist_earth)
#cbar = fig.colorbar(Map, extend='both', shrink=0.95, ax=ax)
#cbar.set_label('value', fontsize = 20)

plt.savefig('Input/Bandspectrum_BOU.png')


#new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_Band4, "Input/NDVI.bil", post_Band4, NDVI)


STOP"""
