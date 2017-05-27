"""This is a Python script that makes a confusion matrix by comparing two arrays.


NB: It displays the results in map and pie chart form.

Later on, we will make it so that it displays multiple pie chart results on a world map

NB: This is actually a function

"""


#Set up display environment in putty
import matplotlib
matplotlib.use('Agg')

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
from DEM_functions import MARSH_ID
from DEM_functions import METRIC1_DEM
from DEM_functions import METRIC2_DEM
from DEM_functions import METRIC3_DEM
from DEM_functions import ENVI_raster_binary_to_2d_array
from DEM_functions import ENVI_raster_binary_from_2d_array




#-------------------------------------------------------------------------------------
# Input arrays: these might have to be rasterized from polygons
Nodata_value = -9999

# Input array 1: the subject
Subject = np.array([[1,0,1,1,1],[1,1,0,1,1],[1,1,1,0,1],[1,1,0,1,1],[Nodata_value,1,1,0,1],[Nodata_value,1,1,1,0]])

# Input array 2: the reference
Reference = np.array([[0,1,1,1,1],[1,0,1,1,1],[1,1,0,1,1],[1,1,0,1,1],[Nodata_value,1,1,0,1],[Nodata_value,1,1,0,0]])

Height = len(Subject[:,0]); Width = len(Subject[0,:])


#-------------------------------------------------------------------------------------
# Make the confusion matrix
Confusion_matrix = Nodata_value*np.ones((Height, Width), dtype = np.float)

for i in range (Height):
    for j in range (Width):
        if Subject[i,j] == 1 and Reference[i,j] == 1: # TRUE POSITIVE
            Confusion_matrix[i,j] = 1
        elif Subject[i,j] == 0 and Reference[i,j] == 0: # TRUE NEGATIVE
            Confusion_matrix[i,j] = 2
        elif Subject[i,j] == 1 and Reference[i,j] == 0: # FALSE POSITIVE
            Confusion_matrix[i,j] = -1
        elif Subject[i,j] == 0 and Reference[i,j] == 1: # FALSE NEGATIVE
            Confusion_matrix[i,j] = -2
            
            
True_positive = np.sum(Confusion_matrix[Confusion_matrix == 1])
print 'TP', True_positive

True_negative = np.sum(Confusion_matrix[Confusion_matrix == 2])/2
print 'TN',True_negative

False_positive = -np.sum(Confusion_matrix[Confusion_matrix == -1])
print 'FP',False_positive

False_negative = -np.sum(Confusion_matrix[Confusion_matrix == -2])/2
print 'FN',False_negative


Reliability = True_positive / (True_positive+False_positive)
Sensitivity = True_positive / (True_positive+False_negative)
print 'R', Reliability
print 'S', Sensitivity


#--------------------------------------------------------------------------------------
# Display the results
from copy import copy
import matplotlib.colors as colors
import matplotlib.mlab as mlab

# Set up a colormap:
palette = copy(plt.cm.RdYlGn)
palette.set_over('r', 3.0)
palette.set_under('g', -3.0)
palette.set_bad(alpha = 0.0)
Confusion_matrix = np.ma.masked_where(Confusion_matrix <= Nodata_value, Confusion_matrix)




#----------------------------------------------------------------
# Figure: Displays map, Pie chart and R/S values
fig=plt.figure(1, facecolor='White',figsize=[10,6])

# a/ Map
ax = plt.subplot2grid((3,2),(0,0),colspan=1, rowspan=3)
Map = ax.imshow(Confusion_matrix, interpolation = 'None', cmap = palette, alpha = 0.8, norm=colors.Normalize(vmin=-2.0, vmax=2.0))
ax.set_title('a. Map of the confusion matrix')


# b/ Pie
ax = plt.subplot2grid((3,2),(0,1),colspan=1, rowspan=2)
labels = 'TP', 'TN', 'FP', 'FN'
sizes = [True_positive, True_negative, False_positive, False_negative]
colormap = [plt.cm.RdYlGn(200), plt.cm.RdYlGn(256), plt.cm.RdYlGn(64), plt.cm.RdYlGn(0)]
ax.pie(sizes, labels=labels, autopct='%1.0f%%',shadow=False, startangle=90, colors = colormap)
ax.axis('equal')

from matplotlib import font_manager as fm
proptease = fm.FontProperties()
proptease.set_size('xx-small')

ax.set_title('b. Performance of the algorithm')

plt.margins(0.3)





# c/ Metrix
ax = plt.subplot2grid((3,2),(2,1),colspan=1, rowspan=1)
metrix = ('R', 'S')
y_pos = np.arange(len(metrix))
performance = np.array([Reliability, Sensitivity])
performance_colour = performance**2

ax.barh(y_pos, performance, align='center', color=plt.cm.RdYlGn(performance_colour))
ax.set_yticks(y_pos); ax.set_yticklabels(metrix)
ax.invert_yaxis(); ax.invert_xaxis()  
ax.set_xlabel('non-dimensional value')
plt.yticks(y_pos, metrix, rotation=0)
ax.set_xlim(0,1)
plt.margins(0.3)
ax.set_title('c. Reliability and Sensitivity of the algorithm')





plt.savefig('Confusion.png')





















