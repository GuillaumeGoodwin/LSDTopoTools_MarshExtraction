"""This is a Python script that calculates the past frequency and interarrival of tidal flooding, and generates a flooding probability (or frequency, tbd) map.

NB: we could make this script into a series of functions
"""

#---------------------------------------------------------------- 
#1. Load useful Python packages
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


import os
import sys

from osgeo import gdal, gdalconst 
from osgeo.gdalconst import * 


from DEM_functions import Open_tide_stats

from DEM_functions import METRIC1_DEM
from DEM_functions import METRIC2_DEM
from DEM_functions import METRIC3_DEM

from DEM_functions import ENVI_raster_binary_to_2d_array
from DEM_functions import ENVI_raster_binary_from_2d_array


#------------------------------------------------------------------
#Use the cookbook
#https://pcjericks.github.io/py-gdalogr-cookbook/


#------------------------------------------------------------------
#
#These are the tide gauges associated to the marshes we work on
Gauges=["BOU", "HIN", "CRO"] # test sample of gauges
#Gauges=["AVO","BOU","CRO","DEV","FEL","HEY","HIN","IMM","LIV","PTM","SHE"] # All the gauges

i=0
Zstep=0.1

# This is for figures
# Create array to store the values of stats for Metric2
# 0=Tidal range; 1; Average elevation of the "marsh"; 2=Standard dev of "marsh" elevation
Metric2_tide_array = np.zeros((3, len(Gauges)), dtype=np.float)
Metric2_array = np.zeros((6, len(Gauges)), dtype=np.float)

#start calculation
start=datetime.datetime.now()

for gauge in Gauges:

    Metric1_tide, Metric2_tide, Metric3_tide, Subsample = Open_tide_stats (gauge) # Load tidal stats
    
    DEM, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array ("LiDAR_DTM_1m/%s/%s_marsh1.bil" % (gauge,gauge), gauge) # Open the DEM file
    
    Channels, post_DEM, envidata_DEM =  ENVI_raster_binary_to_2d_array ("LiDAR_DTM_1m/%s/%s_marsh1_DArea.bil" % (gauge,gauge), gauge) # Open the DArea file

    #Slope, post_Slope, envidata_Slope =  ENVI_raster_binary_to_2d_array ("LiDAR_DTM_1m/%s/%s_marsh1_topo_slope.bil" % (gauge,gauge), gauge) # Open the Slope file
    #Curvature, post_Curvature, envidata_Curvature =  ENVI_raster_binary_to_2d_array ("LiDAR_DTM_1m/%s/%s_marsh1_topo_curvature.bil" % (gauge,gauge), gauge) # Open the Curvature file
    
    #Metric1_DEM = METRIC1_DEM (DEM, Zstep)

    DEM_work = np.copy(DEM)
    DEM_new, DEM_MarshChannels, DEM_TFChannels, Average_Z, Stdev_Z, Area, Perimeter, Channel_length = METRIC2_DEM (DEM_work, Metric2_tide, Channels)

    #print "Saving Metric2 %s" % (gauge)
    #new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, "LiDAR_DTM_1m/%s/%s_Metric2_DEM.bil" % (gauge,gauge), post_DEM, Metric2_DEM)
    
    #new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, "LiDAR_DTM_1m/%s/%s_MARSH_CH_DEM.bil" % (gauge,gauge), post_DEM, DEM_MarshChannels)
    
    #new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, "LiDAR_DTM_1m/%s/%s_TEST_TFCH_DEM.bil" % (gauge,gauge), post_DEM, DEM_TFChannels)
    

    Metric2_array[0,i]=np.mean(Metric2_tide[3])-np.mean(Metric2_tide[0])
    Metric2_array[1,i]=Average_Z[3]
    Metric2_array[2,i]=Stdev_Z[3]
    Metric2_array[3,i]=Area[3]
    Metric2_array[4,i]=Perimeter[3]
    Metric2_array[5,i]=Channel_length[1]
    
    Metric2_tide_array[0,i]=np.mean(Metric2_tide[3])-np.mean(Metric2_tide[0])
    Metric2_tide_array[1,i]=np.mean(Metric2_tide[3])
    Metric2_tide_array[2,i]=np.mean(Metric2_tide[2])


    DEM_work = np.copy(DEM)
    Metric3_DEM = METRIC3_DEM (DEM_work, Metric3_tide)

    print "Saving Metric3 %s" % (gauge)
    new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata_DEM, "LiDAR_DTM_1m/%s/%s_Metric3_DEM.bil" % (gauge,gauge), post_DEM, Metric3_DEM)
    
    i=i+1
    
        
#fig=plt.figure('TEST', facecolor='White',figsize=[7.08,14])
fig=plt.figure('TEST', facecolor='White',figsize=[10,20])
ax2 = plt.subplot2grid((3,1),(1,0),colspan=1)
i=0
for gauge in Gauges:
    Tidal_range=Metric2_tide_array[1, i] - Metric2_tide_array[2, i]
    ax2.plot(Metric2_array[0,i], (Metric2_array[1,i]-Metric2_tide_array[2,i])/Tidal_range,'og')
    ax2.errorbar(Metric2_array[0,i], (Metric2_array[1,i]-Metric2_tide_array[2,i])/Tidal_range, yerr=Metric2_array[2,i]/Tidal_range, color=plt.cm.Greens(200))
    plt.annotate(gauge, xy=(Metric2_array[0,i], (Metric2_array[1,i]-Metric2_tide_array[2,i])/Tidal_range), xytext=((Metric2_array[0,i]+0.2), (((Metric2_array[1,i]-Metric2_tide_array[2,i])/Tidal_range)-0.07)), fontsize=12 )
    i=i+1

ax1 = plt.subplot2grid((3,1),(0,0),colspan=1)
i=0
for gauge in Gauges:
    Tidal_range=np.mean(Metric2_tide[3])-np.mean(Metric2_tide[0])
    ax1.errorbar(Metric2_array[0,i], Metric2_array[1,i], yerr=Metric2_array[2,i], color=plt.cm.Greens(200))
    ax1.plot(Metric2_array[0,i], Metric2_array[1,i],'og')
    plt.annotate(gauge, xy=(Metric2_array[0,i], Metric2_array[1,i]), xytext=((Metric2_array[0,i]+0.2) , ((Metric2_array[1,i])-0.5)), fontsize=10 )
    i=i+1
    
Metric2_tide_array[0]= sorted(Metric2_tide_array[0])
Metric2_tide_array[1]= sorted(Metric2_tide_array[1])
Metric2_tide_array[2]= sorted(Metric2_tide_array[2])
ax1.fill_between(Metric2_tide_array[0], Metric2_tide_array[2], Metric2_tide_array[1], alpha=0.1)

ax1.set_ylabel('Platform elevation (m)', fontsize=14)
ax1.set_xlabel('Spring tidal range (m)', fontsize=14)# ax1.set_xticklabels([])
ax1.grid(True)
ax1.set_xlim(2,14)

ax2.set_ylabel('Elevation relative to the MHWN-MHWS range', fontsize=14)
ax2.set_xlabel('Spring tidal range (m)', fontsize=14) #ax2.set_xticklabels([])
ax2.grid(True)
ax2.set_ylim(0,1)
ax2.set_xlim(2,14)


ax2 = plt.subplot2grid((3,1),(2,0),colspan=1)
i=0
for gauge in Gauges:
    Tidal_range=Metric2_tide_array[1, i] - Metric2_tide_array[2, i]
    ax2.plot(Metric2_array[0,i], Metric2_array[3,i]/Metric2_array[4,i],'og')
    plt.annotate(gauge, xy=(Metric2_array[0,i], Metric2_array[3,i]/Metric2_array[4,i]), xytext=((Metric2_array[0,i]+0.2), ((Metric2_array[3,i]/Metric2_array[4,i])-0.07)), fontsize=10 )
    i=i+1

ax2.set_ylabel('Area/Perimeter of MHWN-MHWS zone (m-1)', fontsize=14)
ax2.set_xlabel('Spring tidal range (m)', fontsize=14)
ax2.grid(True)
ax2.set_xlim(2,14)

#fig.tight_layout()


plt.savefig('TEST.png')





fig=plt.figure('TEST2', facecolor='White',figsize=[10,6])
ax1 = plt.subplot2grid((1,1),(0,0),colspan=1)
i=0
for gauge in Gauges:
    Tidal_range=Metric2_tide_array[1, i] - Metric2_tide_array[2, i]
    ax1.plot(Metric2_array[0,i], Metric2_array[5,i]/Metric2_array[3,i],'og')
    
    plt.annotate (gauge, xy=(Metric2_array[0,i], Metric2_array[5,i]/Metric2_array[3,i]), xytext=((Metric2_array[0,i]+0.2), (Metric2_array[5,i]/Metric2_array[3,i]-0.001)), fontsize=12 )
    i=i+1


ax1.set_ylabel('Channel length/Area of MHWN-MHWS zone (m)', fontsize=14)
ax1.set_xlabel('Tidal range (m)', fontsize=14)

ax1.grid(True)
ax1.set_xlim(2,14)


plt.savefig('TEST2.png')









STOP
for i in Gauges:
    
    
    

    Metric3_DEM = METRIC3_DEM (DEM, Metric3_tide)
    #print "Saving Metric2 %s" % (gauge)
    #new_geotransform, new_projection, file_out = ENVI_raster_binary_from_2d_array (envidata, "LiDAR_DTM_1m/%s/%s_Metric3_DEM.bil" % (gauge,gauge), post, Metric3_DEM)
    
    
    i=i+1


    STOP   

    
    
    
    
   
     
    
    
    
    #fig=plt.figure('PDF_%s' % (gauge), facecolor='White',figsize=[7.08,8])
    #ax1 = plt.subplot2grid((1,1),(0,0),colspan=1)
    #ax1.plot(Metric3_tide[0], Metric3_tide[1],'-')
    
    #plt.show()
    
    
    
   
    """Metric1_DEM = METRIC1_DEM (DEM, Zstep)
    fig=plt.figure('PDF_%s' % (gauge), facecolor='White',figsize=[7.08,8])
    ax1 = plt.subplot2grid((1,1),(0,0),colspan=1)
    
    ax1.plot(Metric1_DEM[0], Metric1_DEM[1],'-')
    ax1.plot(Metric1_tide[0], Metric1_tide[1],'--')
    
    ax1.fill_betweenx(Metric1_DEM[0], np.mean(Metric2_tide[0]), np.mean(Metric2_tide[3]), alpha=0.1)
    ax1.fill_betweenx(Metric1_DEM[0], np.mean(Metric2_tide[1]), np.mean(Metric2_tide[2]), alpha=0.1)
    ax1.set_xlim(xmin=-8, xmax=8)
    ax1.set_ylim(ymin=0, ymax=max(Metric1_DEM[1]))
    ax1.grid(True)

    plt.savefig('PDF_%s.png' %(gauge))"""

    
"""A=DEM.shape
x = np.arange(A[0])
y = np.arange (A[1])
X, Y = np.meshgrid(y,x)

Fig=plt.figure(0, facecolor='White',figsize=[7.08,14])
ax = plt.subplot2grid((1,2),(0,0),colspan=1)

MAP=ax.contourf (Y,X,DEM, cmap=plt.cm.jet,alpha=0.4,interpolate=True)
#ax.set_xlim([7950, 8600])
#ax.set_ylim([2180, 2900])

cbar = plt.colorbar(MAP)
cbar.set_label('Elevation (m)', rotation=270)


ax = plt.subplot2grid((1,2),(0,1),colspan=1)
MAP2=ax.contourf (Y,X,DEM_new, cmap=plt.cm.jet, alpha=0.4, interpolate=True)
#ax.set_xlim([7950, 8600])
#ax.set_ylim([2180, 2900])

cbar = plt.colorbar(MAP2)
cbar.set_label('Submersion time (hrs)', rotation=270)
"""

    

DEM_resolution=1


ELEV= np.arange(np.amin(DEM), np.amax(DEM), 0.05)
Flooded_surface= np.zeros(len(ELEV), dtype=np.float)
Flooded_volume= np.zeros(len(ELEV), dtype=np.float)

print ELEV

for i in range(len(ELEV)):
    print ELEV[i]
    find=np.where(np.logical_and(DEM<ELEV[i], DEM>-2.)) 
    find=np.asarray(find)

    Flooded_surface[i]= len(find[0])*DEM_resolution**2
    Flooded_volume[i]= len(find[0])


    
    #os.system("gdal_translate -a_srs EPSG:27700 -of ENVI " + "name.tif" + " " +  "nametest.bil")
 
    
    
    STOP
    
    

    Z_data=DEM.ravel()
    
    if gauge == "SHE":
        print 'clipping boundaries from', len(Z_data)
        Z_data=np.delete(Z_data, np.where(Z_data==0.))
        print 'to', len(Z_data)
    
    Z_value=np.arange(min(Z_data),max(Z_data),Zstep)
    

    
    

STOP    
    
    
fig=plt.figure(0, facecolor='White',figsize=[7.08,8])

for i in range(len(Gauges)):
    ax = plt.subplot2grid((5,1),(i,0),colspan=1)   
    ax.plot(Zbins[i], Zhist[i],color=plt.cm.jet(60*i))
    #ax.plot(Tidebins[i], Tidehist[i],'--',color=plt.cm.jet(60*i))

    #ax.set_xlim(-7,9)    
    ax.set_xlim(0,15)    
    ax.set_ylim(ymin=0,ymax=0.1) 
    #ax.set_ylim(ymin=0,ymax=max(Tidehist[i])) 
    
    ax.set_xlabel ('DArea')
    
#plt.show() 

plt.savefig('Marsh_Tides_DArea.png')    
 
STOP






















    
#------------------------------------------
# Apply the transformation to the raster
# And throw in the submerged volume as well
ds = gdal.Open("Solway/LiDAR/DTM_1m/CLIPTEST.bil")
DEM = np.array(ds.GetRasterBand(1).ReadAsArray())
DEM[DEM==-1000.]=0.

A=DEM.shape
x = np.arange(A[0])
y = np.arange (A[1])

Sub= np.zeros((len(x),len(y)), dtype=np.float)
Sub_vol= np.zeros((len(x),len(y)), dtype=np.float)

for i in range(len(Th)):
    print Th[i]
    #print Th[i], i,'/',len(Th)
    find=np.where(np.logical_and(DEM>Th[i]-delta_Th, DEM<Th[i]+delta_Th))
    
    
    #this is for the volume bit
    find=np.where(DEM<Th[i])  ##Blah blah, get that done on monday
    
    
    print Sub[find]
    
    
    Sub[find]=Sub_AvTime[0,i]
    #print Sub[find]
    
    
    print Sub_AvTime[0,i]
    

#Diff=DEM-Sub

#print DEM

#print Sub

#print np.where(Diff<>0.)
    

    

    
#STOP
        

    
    

    
    
#----------------------------------------------------
#figures


X, Y = np.meshgrid(y,x)
         
#divider = make_axes_locatable(plt.gca())
                  
                  
                  
Fig=plt.figure(0, facecolor='White',figsize=[7.08,14])
ax = plt.subplot2grid((1,2),(0,0),colspan=1)

MAP=ax.contourf (Y,X,DEM, cmap=plt.cm.jet,alpha=0.4,interpolate=True)
ax.set_xlim([7950, 8600])
ax.set_ylim([2180, 2900])

cbar = plt.colorbar(MAP)
cbar.set_label('Elevation (m)', rotation=270)


ax = plt.subplot2grid((1,2),(0,1),colspan=1)
MAP2=ax.contourf (Y,X,Sub, cmap=plt.cm.jet, alpha=0.4, interpolate=True)
ax.set_xlim([7950, 8600])
ax.set_ylim([2180, 2900])

cbar = plt.colorbar(MAP2)
cbar.set_label('Submersion time (hrs)', rotation=270)



plt.show()
stop


#----------------------------------------------------
# Plot tide height
#Looks like -2.81 is where the harbour gets dry...
Fig=plt.figure(1, facecolor='White',figsize=[7.08,8])
ax = plt.subplot2grid((1,1),(0,0),colspan=1)
ax.plot(Time[0:1000],Tide_data[0:1000], '-') 
ax.fill_between(Time[0:1000],0,Tide_data[0:1000],alpha=0.5)
ax.set_xlabel ('hours')

#-----------------------------------------------------
#Plot the tidal stats
fig=plt.figure(2, facecolor='White',figsize=[7.08,8])
ax = plt.subplot2grid((1,1),(0,0),colspan=1)
ax.plot(HWlevels, 'b-')
ax.plot(LWlevels, 'r-')

#----------------------------------------------------
fig=plt.figure(3, facecolor='White',figsize=[7.08,8])
ax = plt.subplot2grid((1,1),(0,0),colspan=1)

ax.plot (Th,Sub_AvTime[0], 'b-')
ax.plot (Th,Inter_AvTime[0], 'r-')
ax.fill_between (Th,Sub_AvTime[0]-Sub_AvTime[1]/5,Sub_AvTime[0]+Sub_AvTime[1]/5, color= 'blue', alpha=0.2)
ax.fill_between (Th,Inter_AvTime[0]-Inter_AvTime[1]/5,Inter_AvTime[0]+Inter_AvTime[1]/5, color='red', alpha=0.2)

ax.set_yscale('log')

#------------------------------------------------------
fig=plt.figure(4, facecolor='White',figsize=[7.08,8])
ax = plt.subplot2grid((1,1),(0,0),colspan=1)
ax2 = ax.twinx()
ax.plot(Th[:len(Th)-1],hist, '-')
ax2.plot(Th[:len(Th)-1],cdf, 'g-')
  

plt.show()

STOP

"""YOU HAVE STUFF TO SORT OUT!
1. Get that time conversion going. (Fer F***'s sake) -- DONE
2. Figure out what metrics are useful. For now, I think a/flooding duration, b/ interarrival time are good. We can work out flooding probability in a more comprehensive way later.
2.1 Compute area under curve for submerged periods. 
2.2 Compute useful tidal statistics.
2.3 Compute tidal prism and use it to find channels! You can then measure their migration as a measure of dynaimcs for different tidal ranges
3. Figure out how to apply this to a raster. This is basically a static tidal flooding model. Could we later make it dynamic?

4. Once all of this is done, we can figure out how to ralte this to other metrics like slope and curvature.

Question: Do marshes close to a major and flooding river sit higher on the tidal scale?


"""



