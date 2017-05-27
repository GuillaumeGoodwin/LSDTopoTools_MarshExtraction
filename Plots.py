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
import numpy as np
from osgeo import gdal, osr
import matplotlib.ticker as tk
from matplotlib import rcParams
from mpl_toolkits.axes_grid1.inset_locator import *
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
from osgeo import gdal, gdalconst
from osgeo.gdalconst import *
from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab

#------------------------------------------------------------------
# Import the marsh-finding functions
from DEM_functions import Open_tide_stats
from DEM_functions import ENVI_raster_binary_to_2d_array
from DEM_functions import ENVI_raster_binary_from_2d_array
from DEM_functions import add_subplot_axes


#------------------------------------------------------------------
#Use the cookbook
#https://pcjericks.github.io/py-gdalogr-cookbook/

#------------------------------------------------------------------
#These are the tide gauges next to the marshes we work on

Gauges=["BOU", "FEL", "CRO", "SHE", "WOR", "HEY", "HIN"] # ALL gauges by tidal range
#Gauges=["CRO", "SHE", "WOR", "HEY", "HIN"] # ALL gauges by tidal range
#Gauges=["WOR", "HEY", "HIN"] # ALL gauges by tidal range

#Gauges=["BOU", "HEY", "HIN"] # ALL gauges by tidal range



hist = []
bins = []


HIST=[]
BINS=[]

Inflexion = np.zeros(len(Gauges), dtype = np.float)
Inflexion_point = np.zeros(len(Gauges), dtype = np.float)


histS = []
binsS = []

Nodata_value = -9999




Metrix_gauges = np.zeros((len(Gauges),5), dtype=np.float)
i=0

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
    Platform, post_Slope, envidata_Slope =  ENVI_raster_binary_to_2d_array ("Output/%s/%s_Marsh.bil" % (gauge,gauge), gauge)
    #print " Loading Confusion"
    #Confusion_matrix, post_Curvature, envidata_Curvature =  ENVI_raster_binary_to_2d_array ("Output/%s/%s_Confusion_DEM.bil" % (gauge,gauge), gauge)


    print "Loading performance files"
    with open ("Output/%s/%s_Performance.pkl" % (gauge,gauge), 'rb') as input_file:
        Performance = cPickle.load(input_file)
    with open ("Output/%s/%s_Metrix.pkl" % (gauge,gauge), 'rb') as input_file:
        Metrix = cPickle.load(input_file)

    Metrix_gauges[i,0] = np.mean (Metric2_tide[3])-np.mean (Metric2_tide[0])
    for j in np.arange(1,5,1):
        Metrix_gauges[i,j] = Metrix[j-1]




    Relief = DEM-np.amin(DEM[DEM > Nodata_value])
    Rel_relief = Relief/np.amax(Relief)
    Rel_relief[DEM == Nodata_value] = Nodata_value

    #Scarps_relief = DEM-np.amin(DEM[DEM > Nodata_value])
    #Scarps_rel_relief = Scarps_relief/np.amax(Scarps_relief)
    #Scarps_rel_relief[Scarps == 0] = Nodata_value
    #Scarps_rel_relief[Scarps_rel_relief > 100] = Nodata_value




    Rel_slope = Slope/np.amax(Slope)
    Rel_slope[Slope == Nodata_value] = Nodata_value

    #Scarps_rel_slope = Slope/np.amax(Slope)
    #Scarps_rel_slope[Slope == Nodata_value] = Nodata_value
    #Scarps_rel_slope[Scarps == 0] = 0


    Crossover = Rel_relief * Rel_slope
    Crossover_copy = np.copy(Crossover)
    Crossover[Search_space == 0] = 0
    Crossover[DEM == Nodata_value] = Nodata_value
    Crossover_copy[DEM == Nodata_value] = Nodata_value


    #Scarps[Scarps>1]=0

    #Scarps_crossover = Scarps * Scarps_rel_slope
    #Scarps_crossover_copy = np.copy(Scarps_crossover)
    #Scarps_crossover[Scarps == 0] = Nodata_value
    #Scarps_crossover[DEM == Nodata_value] = Nodata_value
    #Scarps_crossover_copy[DEM == Nodata_value] = Nodata_value



    #--------------------------------------------------------------------------------------
    # Make the colour palettes


    #DEM colourmap
    palette_DEM = copy(plt.cm.gist_earth)
    palette_DEM.set_bad(alpha = 0.00)
    DEM = np.ma.masked_where(DEM <= Nodata_value, DEM)


    # Slopes colourmap
    palette_Slope = copy(plt.cm.gist_heat)
    palette_Slope.set_bad(alpha = 0.00)
    Slope = np.ma.masked_where(Slope <= Nodata_value, Slope)


    #Rel_relief colourmap
    #palette_Rel_relief = copy(plt.cm.gist_earth)
    #palette_Rel_relief.set_bad(alpha = 0.00)
    #Rel_relief = np.ma.masked_where(Rel_relief <= Nodata_value, Rel_relief)


    # Rel_slope colourmap
    #palette_Rel_slope = copy(plt.cm.gist_heat)
    #palette_Rel_slope.set_bad(alpha = 0.00)
    #Rel_slope = np.ma.masked_where(Rel_slope <= Nodata_value, Rel_slope)


    # Crossover colourmap
    #palette_Crossover = copy(plt.cm.seismic)
    #palette_Crossover.set_bad(alpha = 0.00)
    #Crossover = np.ma.masked_where(Crossover <= Nodata_value, Crossover)







    # Scarps colourmap
    #Scarps[Scarps > 0] =  Slope [Scarps > 0]
    #Scarps[Scarps == 0] = 0.00001
    Void = np.where (Slope == Nodata_value)
    Scarps[Void] = Nodata_value
    Smax = np.amax(Scarps); Smin = 0
    palette_Scarps = copy(plt.cm.jet)
    palette_Scarps.set_bad(alpha = 0.00)
    Scarps = np.ma.masked_where(Scarps <= Nodata_value, Scarps)





    """# Platform map colourmap
    #Platform[Platform > 0] = DEM [Platform > 0]
    Void = np.where (DEM == Nodata_value)
    Platform[Void] = Nodata_value

    Zmax = np.amax(Platform); Zmin = np.amin(Platform[Platform != Nodata_value])
    #Platform[Platform == 0] = Zmin-1
    Marsh = np.where(np.logical_and(Platform != 0, Platform != Zmin-1))
    Zmax = np.amax(Platform); Zmin = np.amin(Platform[Platform > Zmin-1])
    palette_platform = copy(plt.cm.gist_earth)
    palette_platform.set_under('k', Zmin)
    palette_platform.set_bad(alpha = 0.00)
    Platform = np.ma.masked_where(Platform <= Nodata_value, Platform)"""


    # Confusion map colourmap
    #palette_Confusion = copy(plt.cm.RdYlGn)
    #palette_Confusion.set_bad(alpha = 0.0)
    #Confusion_matrix = np.ma.masked_where(Confusion_matrix <= Nodata_value, Confusion_matrix)


    #---------------------------------------------------------------------------
    #Paper Plots
    fig=plt.figure(11, facecolor='White',figsize=[30,15])
    ax = plt.subplot2grid((1,2),(0,0),colspan=1, rowspan=2)
    ax.set_title('Slope (from polynomial fit)', fontsize = 22)
    ax.set_xlabel('X-distance from origin (m)', fontsize = 18)
    ax.set_ylabel('Y-distance from origin (m)', fontsize = 18)
    Map = ax.imshow(Rel_relief, interpolation='None', cmap=palette_Slope, vmin=0, vmax=1)#, norm=colors.Normalize(vmin=Zmin, vmax=Zmax))
    #Map = ax.imshow(Crossover, interpolation='None', cmap=palette_Slope, vmin=0, vmax=1)#, norm=colors.Normalize(vmin=Zmin, vmax=Zmax))
    #Map = ax.imshow(Search_space, interpolation='None', cmap=palette_Slope, vmin=0, vmax=1)#, norm=colors.Normalize(vmin=Zmin, vmax=Zmax))

    #Map = ax.imshow(Rel_relief, interpolation='None', cmap=palette_Rel_relief,  vmin=0, vmax=1)#, norm=colors.Normalize(vmin=Smin, vmax=Smax))
    #cbar = fig.colorbar(Map, extend='both', shrink=0.95, ax=ax)
    #cbar.set_label('Slope (m/m)', fontsize = 20)


    ax = plt.subplot2grid((1,2),(0,1),colspan=1, rowspan=2)
    ax.set_title('Detected scarps', fontsize = 22)
    ax.set_xlabel('X-distance from origin (m)', fontsize = 18)
    ax.set_ylabel('Y-distance from origin (m)', fontsize = 18)
    #Map = ax.imshow(Scarps, interpolation='None', cmap=plt.cm.gist_heat)
    #Map = ax.imshow(Rel_slope, interpolation='None', cmap=palette_Rel_slope,  vmin=0, vmax=1)#, norm=colors.Normalize(vmin=Zmin, vmax=Zmax))
    #Map = ax.imshow(Scarps, interpolation='None', cmap=palette_Scarps)#, norm=colors.Normalize(vmin=Smin, vmax=Smax))
    #cbar = fig.colorbar(Map, extend='both', shrink=0.95, ax=ax)
    #cbar.set_label('Scarp slope (m/m)', fontsize = 20)

    plt.savefig('Output/Paper/%s_Paper_fig11.png' % (gauge))





    fig=plt.figure(13, facecolor='White',figsize=[60,30])
    ax = plt.subplot2grid((1,2),(0,0),colspan=1, rowspan=2)
    ax.set_title('Slope (from polynomial fit)', fontsize = 22)
    ax.set_xlabel('X-distance from origin (m)', fontsize = 18)
    ax.set_ylabel('Y-distance from origin (m)', fontsize = 18)
    Map = ax.imshow(Platform, interpolation='None', cmap=plt.cm.gist_heat, vmin=0, vmax=3)
    cbar = fig.colorbar(Map, extend='both', shrink=0.95, ax=ax)
    cbar.set_label('XXXX', fontsize = 20)

    ax = plt.subplot2grid((1,2),(0,1),colspan=1, rowspan=2)
    ax.set_title('Detected scarps', fontsize = 22)
    ax.set_xlabel('X-distance from origin (m)', fontsize = 18)
    ax.set_ylabel('Y-distance from origin (m)', fontsize = 18)
    Map = ax.imshow(Scarps, interpolation='None', cmap=palette_Scarps)#, vmin= 0, vmax=1)
    cbar = fig.colorbar(Map, extend='both', shrink=0.95, ax=ax)
    cbar.set_label('Scarp slope (m/m)', fontsize = 20)

    plt.savefig('Output/Paper/%s_Paper_fig15.png' % (gauge))








    #Platform[Platform > 0] = DEM[Platform > 0]

    Search_x = Crossover
    Z_data = Search_x.ravel(); Z_data = Z_data[Z_data>0]


    print Z_data
    print max(Z_data)
    print min(Z_data)


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

plt.savefig('Output/Paper/Crossover.png')





for i in [0,1]:
    STOP






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




# Setting up the spectral analysis
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


STOP
