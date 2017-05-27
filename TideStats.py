"""This is a Python script that calculates flooding statistics from tide gauge data."""


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
from osgeo import gdal
import matplotlib.ticker as tk
from matplotlib import rcParams
from mpl_toolkits.axes_grid1.inset_locator import *
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable


from Tide_functions import Open_tide_data
from Tide_functions import METRIC1
from Tide_functions import METRIC2
from Tide_functions import SUBSAMPLE
from Tide_functions import METRIC3

#------------------------------------------------------------------
#2. Read data from those files
# NO DATA value is -100

Gauges=["BOU"]
#Gauges=["BOU","FEL","SHE","HEY","AVO"]
#Gauges=["AVO","BOU","CRO","DEV","FEL","HEY","HIN","IMM","LEI","LIV","PTM","SHE","WOR"]


#------------------------------------------------------------------
# Variable names
Tidestep=0.05

#------------------------------------------------------------------
#Calculate metrics

for gauge in Gauges:
    Tide_data, Tide_astrodata, Tide_value, Tide_astrovalue = Open_tide_data ('Tide_gauges_BODC/%s_TideData.pkl'% (gauge), Tidestep)
    
    Metric2_time, Metric2 = METRIC2 (Tide_data, Tide_value)
    
    Subsample, bandwidth_pos, bandwidth_neg=SUBSAMPLE(Tide_data, Metric2)
    #Subsample_astro, bandwidth_pos, bandwidth_neg=SUBSAMPLE(Tide_astrodata, Metric2)

    
    
    
    

    #Metric1 = METRIC1 (Subsample, Tide_value)
    Metric3 = METRIC3 (Subsample, Tide_value)
    
    
    
    
    
    Index1 = np.where (Metric3[1]>12)
    Threshold1 = min(Metric3[0, min(Index1)])
    
    Index2 = np.where (Metric3[1]>24)
    Threshold2 = min(Metric3[0, min(Index2)])
    
    Index3 = np.where (Metric3[1]>72)
    Threshold3 = min(Metric3[0, min(Index3)])
    
    print Threshold1, Threshold2, Threshold3 
    
    
    
    
    
    

    
    
    fig=plt.figure('1_%s' % (gauge), facecolor='White',figsize=[7.08,16])
    
    ax = plt.subplot2grid((1,1),(0,0),colspan=1)
    ax2=ax.twinx()
    
    ax2.plot (Tide_data[:,0],Tide_data[:,1],'r-')
    ax.plot (Tide_astrodata[:,0],Tide_astrodata[:,1],'b-')
    ax2.plot (Tide_data[:,0],Tide_data[:,1]-Tide_astrodata[:,1],'r--')

    ax.fill_between(Tide_data[:,0], np.mean(Metric2[2]),np.mean(Metric2[3]) , alpha=0.1)
    ax.fill_between(Tide_data[:,0], np.mean(Metric2[2]),np.mean(Metric2[1]) , alpha=0.2)
    ax.fill_between(Tide_data[:,0], np.mean(Metric2[1]),np.mean(Metric2[0]) , alpha=0.3)
    
    ax.plot(Tide_data[:,0], Threshold1* np.ones(len(Tide_data[:,0]), dtype=np.float) , 'r')
    ax.plot(Tide_data[:,0], Threshold2* np.ones(len(Tide_data[:,0]), dtype=np.float) , 'r')
    ax.plot(Tide_data[:,0], Threshold3* np.ones(len(Tide_data[:,0]), dtype=np.float) , 'r')
  
    
    
    

    fig.autofmt_xdate()
    ax.set_xlim([datetime.date(2009, 3, 3), datetime.date(2009, 3, 5)]) 
    ax.set_ylim(np.mean(Metric2[0])-0.1, np.mean(Metric2[3]+0.1)) 
    ax2.set_ylim(np.mean(Metric2[0])-0.1, np.mean(Metric2[3]+0.1)) 
    
    ax.annotate('MHWS', xy=(0.01,0.95), xycoords='axes fraction',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+0) 
    ax.annotate('MHWN', xy=(0.01,0.635), xycoords='axes fraction',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+0) 
    ax.annotate('MLWN', xy=(0.01,0.415), xycoords='axes fraction',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+0) 
    ax.annotate('MLWS', xy=(0.01,0.03), xycoords='axes fraction',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+0) 
    
    ax.annotate('3 days without flooding', xy=(0.01,0.82), xycoords='axes fraction',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+0) 
    ax.annotate('1 day without flooding', xy=(0.01,0.74), xycoords='axes fraction',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+0) 
    ax.annotate('Flooded twice a day', xy=(0.01,0.685), xycoords='axes fraction',horizontalalignment='left', verticalalignment='top', fontsize=rcParams['font.size']+0) 
    
    
    ax.set_ylabel(r'Astronomic tidal elevation ($m$)', fontsize=12, color='b')
    ax2.set_ylabel(r'Observed tidal elevation ($m$) and storm surge (--)', fontsize=12, color='r')
    ax.set_xlabel('Time (March 3 to 4, 2009)', fontsize=12)
    
    ax.grid(True)
    
    
    
    #plt.show()
    
    plt.savefig('POSTER.jpg')
    
    
    STOP
    
    print 'Saving stats!'
    cPickle.dump(Metric1,open("Tide_gauges_BODC/%s_Metric1.pkl" % (gauge), "wb"))
    cPickle.dump(Metric2_time,open("Tide_gauges_BODC/%s_Metric2_time.pkl" % (gauge), "wb"))
    cPickle.dump(Metric2,open("Tide_gauges_BODC/%s_Metric2.pkl" % (gauge), "wb"))
    cPickle.dump(Metric3,open("Tide_gauges_BODC/%s_Metric3.pkl" % (gauge), "wb"))
    cPickle.dump(Subsample,open("Tide_gauges_BODC/%s_Subsample.pkl" % (gauge), "wb"))
    

STOP    
    
    
    
    
    
    
    
"""# TEST DATA
Tide_data=np.zeros((24*4,2),dtype=datetime.datetime)
for j in range (len(Tide_data[:,0])):
    Tide_data[j,0]=datetime.datetime(2000,1,1,0,0,0)+datetime.timedelta(minutes=15*j)
    Tide_data[j,1] = j/12 * mt.sin(3.14 * j/12)
        
Tide_value=np.arange(min(Tide_data[:,1]),max(Tide_data[:,1]),Tidestep)
#print min(Tide_data[:,1]),max(Tide_data[:,1])"""

    
    
    
fig=plt.figure('1_%s' % (gauge), facecolor='White',figsize=[7.08,8])
ax = plt.subplot2grid((1,1),(0,0),colspan=1)
ax.plot (Tide_data[:,0],Tide_data[:,1],'b-')
ax.plot (Subsample[:,0],Subsample[:,1],'r--')

ax.fill_between(Subsample[:,0], np.mean(Metric2[0]),np.mean(Metric2[3]) , alpha=0.1)
ax.fill_between(Subsample[:,0], np.mean(Metric2[1]),np.mean(Metric2[2]) , alpha=0.1)

fig.autofmt_xdate()
ax.grid(True)


fig=plt.figure('2_%s' % (gauge), facecolor='White',figsize=[7.08,8])
ax1 = plt.subplot2grid((2,2),(0,0),colspan=1)
ax2 = plt.subplot2grid((2,2),(0,1),colspan=1)
ax3 = plt.subplot2grid((2,2),(1,0),colspan=1)
ax4 = plt.subplot2grid((2,2),(1,1),colspan=1)

ax1.plot(Metric1[0], Metric1[1],'-',color=plt.cm.jet(5+100*i))

ax2.plot(Metric1[0], Metric1[2],'-',color=plt.cm.jet(5+100*i))

ax3.plot(Metric3[0], Metric3[1],'-',color=plt.cm.jet(5+100*i))
ax3.plot(Metric3[0], Metric3[3],'--',color=plt.cm.jet(5+100*i))

ax4.plot(Metric3[0], Metric3[5],'-',color=plt.cm.jet(5+100*i))

ax.grid(True)

i=i+1



fig=plt.figure('2_%s' % (gauge), facecolor='White',figsize=[7.08,8])
ax = plt.subplot2grid((1,1),(0,0),colspan=1)

ax.plot (Metric2_time, Metric2[0])
ax.plot (Metric2_time, Metric2[1])
ax.plot (Metric2_time, Metric2[2])
ax.plot (Metric2_time, Metric2[3])

fig.autofmt_xdate()
ax.set_xlim([datetime.date(2010, 2, 25), datetime.date(2010, 3, 10)]) 
ax.grid(True)

#plt.show()
#STOP
    
    
    #Metric1 = METRIC1 (Tide_data, Tide_value)
    
    #ax2.set_xlim([datetime.date(2010, 2, 25), datetime.date(2010, 3, 10)])    
    #ax2.set_ylim(-1.5,1.5) 
    
    #ax = plt.subplot2grid((2,2),(0,0),colspan=1)
    #ax2 = plt.subplot2grid((2,2),(1,0),colspan=1)
    #ax3 = plt.subplot2grid((1,2),(0,1),colspan=1)



    #ax3 = plt.subplot2grid((2,2),(0,1),rowspan=4)
    #ax4 = plt.subplot2grid((2,2),(0,0),colspan=1)
    #ax5 = plt.subplot2grid((2,2),(1,0),colspan=1)

    
    #ax3.plot (Tide_data[:,0],Tide_data[:,1])
    

    #fig.autofmt_xdate()
    #ax2.set_xlim([datetime.date(2010, 2, 25), datetime.date(2010, 3, 10)])    
    #ax2.set_ylim(-1.5,1.5) 
 
    

        
    #ax.plot(Metric1[0], Metric1[1],'-',color=plt.cm.jet(5+100*i))
    #ax2.plot(Metric1[0], Metric1[2],'-',color=plt.cm.jet(5+100*i))
    
    #ax3.plot (Tide_data[:,0],Tide_data[:,1])
    
    #ax3.plot (Metric2_time,Metric2[0],'-o',color=plt.cm.jet(5+100*i))
    #ax3.plot (Metric2_time,Metric2[1],'-o',color=plt.cm.jet(5+100*i))
    #ax3.plot (Metric2_time,Metric2[2],'-o',color=plt.cm.jet(5+100*i))
    #ax3.plot (Metric2_time,Metric2[3],'-o',color=plt.cm.jet(5+100*i))
    
    #ax4.plot(Metric3[0], Metric3[1],'-',color=plt.cm.jet(5+100*i))
    #ax4.plot(Metric3[0], Metric3[3],'-',color=plt.cm.jet(5+100*i))
    
    #ax5.plot(Metric3[0], Metric3[5],'-',color=plt.cm.jet(5+100*i))
    
    #i=i+1

#ax.set_xlim(min(Tide_value),max(Tide_value))
#ax2.set_xlim(min(Tide_value),max(Tide_value))  
#ax4.set_xlim(min(Tide_value),max(Tide_value)); ax4.set_ylabel('Time (hours)') #ax4.set_yscale('log'); 
#ax5.set_xlim(min(Tide_value),max(Tide_value)); ax5.set_ylabel('height (m)')

#fig.autofmt_xdate()
#ax3.set_xlim([datetime.date(2000, 1, 1), datetime.date(2000, 1, 2)])    

#plt.show()

    


#----------------------------------------------------
# Save metrics













      
#--------------------------------------------------------------
fig=plt.figure(1, facecolor='White',figsize=[7.08,8])
ax = plt.subplot2grid((1,2),(0,0),colspan=1)





for i in range(len(Gauges)):

    print i
    
    #print METRIC3[0][i]
    #print METRIC3[1][i]
    
    #ax.fill_betweenx(-4,4)
    
    ax.plot (METRIC3[0][i],METRIC3[1][i],'b')#,color=plt.cm.jet(5+100*i))
    ax.plot (METRIC3[0][i],METRIC3[3][i],'or')#,color=plt.cm.jet(5+100*i))

ax.set_yscale('log')

#fig.autofmt_xdate()   
#ax.set_ylim([datetime.date(1950, 1, 1), datetime.date(2015, 1, 1)])

ax.set_xlim(-3,3)







ax = plt.subplot2grid((1,2),(0,1),colspan=1)
for i in range(len(Gauges)):
    #ax.plot(High_tides_time,High_tides,'--b')
    #ax.plot(Low_tides_time,Low_tides,'--b')
    #ax.plot(Tide_data[:,0],Tide_data[:,2],'r')
    print i
    
    ax.plot (METRIC2[0][i],METRIC2[1][i],color=plt.cm.jet(5+100*i))
    ax.plot (METRIC2[0][i],METRIC2[2][i],color=plt.cm.jet(5+100*i))
    ax.plot (METRIC2[0][i],METRIC2[3][i],color=plt.cm.jet(5+100*i))
    ax.plot (METRIC2[0][i],METRIC2[4][i],color=plt.cm.jet(5+100*i))
    
    #ax.plot (METRIC2_astro[0][i],METRIC2_astro[1][i],'--',color=plt.cm.jet(5+100*i))
    #ax.plot (METRIC2_astro[0][i],METRIC2_astro[2][i],'--',color=plt.cm.jet(5+100*i))
    #ax.plot (METRIC2_astro[0][i],METRIC2_astro[3][i],'--',color=plt.cm.jet(5+100*i))
    #ax.plot (METRIC2_astro[0][i],METRIC2_astro[4][i],'--',color=plt.cm.jet(5+100*i))
    
    


fig.autofmt_xdate()
ax.set_xlim([datetime.date(1950, 1, 1), datetime.date(2015, 1, 1)])    


plt.show()




plt.show()        
                

    
    
STOP
                    
                 
                    
                    
            
            
   
        
        
        
        
        
        
        
        
        
    
#--------------------------------------------------------------
fig=plt.figure(1, facecolor='White',figsize=[7.08,8])
ax = plt.subplot2grid((1,1),(0,0),colspan=1)
for i in range(len(Gauges)):
    #ax.plot(High_tides_time,High_tides,'--b')
    #ax.plot(Low_tides_time,Low_tides,'--b')
    #ax.plot(Tide_data[:,0],Tide_data[:,2],'r')
    print i
    
    ax.plot (METRIC2[0][i],METRIC2[1][i],color=plt.cm.jet(5+100*i))
    ax.plot (METRIC2[0][i],METRIC2[2][i],color=plt.cm.jet(5+100*i))
    ax.plot (METRIC2[0][i],METRIC2[3][i],color=plt.cm.jet(5+100*i))
    ax.plot (METRIC2[0][i],METRIC2[4][i],color=plt.cm.jet(5+100*i))
    
    #ax.plot (METRIC2_astro[0][i],METRIC2_astro[1][i],'--',color=plt.cm.jet(5+100*i))
    #ax.plot (METRIC2_astro[0][i],METRIC2_astro[2][i],'--',color=plt.cm.jet(5+100*i))
    #ax.plot (METRIC2_astro[0][i],METRIC2_astro[3][i],'--',color=plt.cm.jet(5+100*i))
    #ax.plot (METRIC2_astro[0][i],METRIC2_astro[4][i],'--',color=plt.cm.jet(5+100*i))
    
    


fig.autofmt_xdate()
ax.set_xlim([datetime.date(1950, 1, 1), datetime.date(2015, 1, 1)])    


plt.show()
        
#plt.savefig('TESTING_%s.png' %(gauge))

        
STOP
        



#----------------------------------------------------
# Save metrics
print 'Saving Tidal statistics!'
cPickle.dump(METRIC1,open("METRIC1.pkl", "wb"))
cPickle.dump(METRIC1_astro,open("METRIC1_astro.pkl", "wb"))


#-----------------------------------------------------------------------------
# Print figures
  
fig=plt.figure(1, facecolor='White',figsize=[7.08,8])
ax = plt.subplot2grid((2,1),(0,0),colspan=1)
for i in range(len(Gauges)):
    ax.plot(METRIC1_astro[0][i], METRIC1_astro[1][i], color=plt.cm.jet(60*i))
    ax.plot(METRIC1[0][i], METRIC1[1][i],'--',color=plt.cm.jet(60*i))
ax.set_xlim(-7,9)

ax = plt.subplot2grid((2,1),(1,0),colspan=1)
for i in range(len(Gauges)):
    ax.plot(METRIC1_astro[0][i], METRIC1_astro[2][i], color=plt.cm.jet(60*i))
    ax.plot(METRIC1[0][i], METRIC1[2][i],'--',color=plt.cm.jet(60*i))
    
ax.set_xlim(-7,9) 


plt.show()

STOP     


        
        
        
        
        
        
        
  
#--------------------------------------------------------------
fig=plt.figure(1, facecolor='White',figsize=[7.08,8])
ax = plt.subplot2grid((1,1),(0,0),colspan=1)
for i in range(len(Gauges)):
    ax.plot(High_tides_time,High_tides,'--b')
    ax.plot(Low_tides_time,Low_tides,'--b')
    ax.plot(Tide_data[:,0],Tide_data[:,2],'r')

    ax.fill_between(Tide_data[:,0],Spring_LT,Spring_HT, alpha=0.1)
    ax.fill_between(Tide_data[:,0],Neap_LT,Neap_HT, alpha=0.1)

fig.autofmt_xdate()
ax.set_xlim([datetime.date(2006, 1, 1), datetime.date(2006, 2, 1)])    
ax.set_ylim(Min,Max)
            

plt.savefig('GRAPH_%s.png' %(gauge))
        
print 'DONE'

        
        
        
        
    
