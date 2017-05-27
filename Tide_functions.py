"""This file contains useful functions to extract flooding statistics from tidal data"""

#---------------------------------------------------------------- 
#1. Load useful Python packages
import numpy as np
import functools
import math as mt
import cmath
import scipy as sp
import scipy.stats as stats
from datetime import datetime
import cPickle
from pylab import *


#----------------------------------------------------------
# This function opens the cPickle file and returns the data
# Call this function with Tide_data, Tide_astrodata, Tide_value, Tide_astrovalue = Open_tide_data ('Tide_gauges_BODC/%s_TideData.pkl'% (gauge), Tidestep)

def Open_tide_data (Inputfile, Tidestep):
    print 'Opening', Inputfile
    with open (Inputfile, 'rb') as input_file:
        Data = cPickle.load(input_file)
        Tide_data = np.delete(Data,(2), axis=1);Tide_astrodata = np.delete(Data,(1), axis=1)
        Max=Tide_data[:,1]; Max=np.max(Max[Max<15]); Min=Tide_data[:,1]; Min=np.min(Min[Min>-15])
        Max_astro=Tide_astrodata[:,1]; Max_astro=np.max(Max_astro[Max_astro<15]); Min_astro=Tide_astrodata[:,1]; Min_astro=np.min(Min_astro[Min_astro>-15])
        Tide_value= np.arange(Min,Max+Tidestep,Tidestep)
        Tide_astrovalue= np.arange(Min_astro,Max_astro+Tidestep,Tidestep)
        
    return Tide_data, Tide_astrodata, Tide_value, Tide_astrovalue





#----------------------------------------------------------
# This function calculates the METRIC1 array for one gauge
# Call this function with METRIC1 = METRIC1 (Tide_data, Tide_value) or with the 'astro' versions
# METRIC1 dimensions are 3*length of Tide_value (0 = elevation bins; 1 = pdf of bins; 2 = cdf of each bin)

def METRIC1 (Tide_data, Tide_value):
    print 'Calculating pdf and cdf'
    METRIC1=np.zeros((3, len(Tide_value)-1),dtype=np.float)
    hist,bins = np.histogram (Tide_data[:,1],Tide_value,density=True); hist=hist/sum(hist); bins=bins[:-1]; cdf=[0] 
    for i in range(1,len(hist)-1):
        if hist[i]>=5*hist[i-1] and hist[i]>=5*hist[i+1]:          
            hist[i]=0
        cdf.append(hist[i]+cdf[i-1])
    cdf.append(hist[-1]+cdf[-2])
    METRIC1[0]=bins; METRIC1[1]=hist; METRIC1[2]=cdf

    return METRIC1





#----------------------------------------------------------
# This function calculates the METRIC2 array for one gauge
# Call this function with METRIC2_time, METRIC2 = METRIC2 (Tide_data, Tide_value) or with the 'astro' versions 
# METRIC2 dimensions are 5*Nyears of data for each gauge (0=year; 1 = spring low tide; 2 = neap low tide; 3 = median value; 4 = neap high tide; 5 = spring high tide)

def METRIC2 (Tide_data, Tide_value):
    print 'Calculating yearly tide stage percentiles (95,5)'
    
    Metric2=[];Metric2.append([]);Metric2.append([]);Metric2.append([]);Metric2.append([]);Metric2.append([]); Metric2.append([])    
    High_tides=[]; Low_tides=[]; Mid_tides=[]; High_tides_time=[]; Low_tides_time=[]; Mid_tides_time=[]

    for i in range (1,len(Tide_data)-1):
        if Tide_data[i,1]>Tide_data[i-1,1] and Tide_data[i,1]>Tide_data[i+1,1] and Tide_data[i,1]>0 and Tide_data[i,1]<10:
            High_tides.append(Tide_data[i,1])
            High_tides_time.append(Tide_data[i,0]) 
        elif Tide_data[i,1]<Tide_data[i-1,1] and Tide_data[i,1]<Tide_data[i+1,1] and Tide_data[i,1]<0 and Tide_data[i,1]>-10:
            Low_tides.append(Tide_data[i,1])
            Low_tides_time.append(Tide_data[i,0])

        if Tide_data[i,0].year > Tide_data[i-1,0].year: 
            YEAR= datetime.datetime (Tide_data[i,0].year,1,1,1,0,0,0)
            Metric2[0].append(YEAR)
            Metric2[1].append(np.percentile(Low_tides,5)); Metric2[2].append(np.percentile(Low_tides,95)); Metric2[3].append(np.percentile(High_tides,5)); Metric2[4].append(np.percentile(High_tides,95))
            
            High_tides=[]; High_tides_time=[]; Low_tides=[]; Low_tides_time=[]; Mid_tides=[]; Mid_tides_time=[]

    METRIC2_time=Metric2[0]; METRIC20=Metric2[1]; METRIC21=Metric2[2]; METRIC22=Metric2[3]; METRIC23=Metric2[4]
    METRIC2=np.row_stack((METRIC20, METRIC21, METRIC22, METRIC23))

    return METRIC2_time, METRIC2





#-----------------------------------
# This function calculates the METRIC3 array for one gauge
# Call this function with Subsample, bandwidth_pos, bandwidth_neg = SUBSAMPLE (Tide_data, Metric2) or with the 'astro' versions 
def SUBSAMPLE (Tide_data, Metric2):
    print 'Sampling the most recent year of complete data'
    bandwidth_pos=1.5 * np.mean (Metric2[3]); bandwidth_neg=1.5 * np.mean (Metric2[0])
    
    for i in range(len(Tide_data)-1,0,-1):
        t1=Tide_data[i,0]; t2=Tide_data[i-1,0]; z1=Tide_data[i,1]; z2=Tide_data[i-1,1]
        if z1*z2 <= 0 and z2 < bandwidth_pos and z2 > bandwidth_neg:
            Start_time=t2; Start_index=i-1
            End_time = Start_time - datetime.timedelta(days=365); End_index = np.where (Tide_data[:,0]==End_time); End_index=End_index[0]

            if max(Tide_data[End_index:Start_index,1])/bandwidth_pos > 1 or min(Tide_data[End_index:Start_index,1])/bandwidth_neg > 1:
                print max(Tide_data[End_index:Start_index,1])/bandwidth_pos, min(Tide_data[End_index:Start_index,1])/bandwidth_neg
                print 'Initialising new subsample'

            else:
                Subsample = Tide_data[End_index:Start_index]
                print 'Done', Start_time, End_time, Start_time-End_time
                break

    return Subsample, bandwidth_pos, bandwidth_neg





#-----------------------------------
# This function calculates the METRIC3 array for one gauge
# Call this function with METRIC3 = METRIC3 (Tide_data, Tide_value) or with the 'astro' versions 
# METRIC3 dimensions are 8*length of Tide_values (0 = elevations bins; 1 = Average Interarrival Duration; 2 = StDev of Interarrival Duration; 3 = Average Flooding Duration; 4 = StDev of Flooding Duration; 5 = Average Flooded Volume; 6 = Total Flooded Volume; 7 = Flooding Fraction)


"""THE HEIGHT CALCULATION BIT IS STILL DODGY BUT WE CAN DO WITHOUT IT FOR NOW"""

def METRIC3 (Tide_data, Tide_value):
    print 'Calculating flooding parameters'
    Metric3=[]; Metric3.append([]); Metric3.append([]); Metric3.append([]); Metric3.append([]); Metric3.append([]); Metric3.append([]); Metric3.append([]); Metric3.append([])

    for i in range (1,len(Tide_value)):
        print i,'/',len(Tide_value)
        #Reset DRY_TIME, FLOODED_TIME and FLOODED_HEIGHT
        Dry_time=0; Flooded_time=0; Flooded_height=0
        DRY_TIME=[]; FLOODED_TIME=[]; FLOODED_HEIGHT=[]

        for j in range (1,len(Tide_data)): #this is the height value checked against threshold
            t1=Tide_data[j-1,0]; t2=Tide_data[j,0]; Timestep= (t2-t1).total_seconds()
            a1=Tide_data[j-1,1]-Tide_value[i]; a2=Tide_data[j,1]-Tide_value[i]; step_fraction = abs(a2)/(abs(a1)*(1+abs(a2)/abs(a1)))
            
            if a1<-10 or a2<-10: #this is the safety case against nodata values
                pass

            elif a1*a2<0 and a2>0: #Meaning you cross the threshold from under when actual data is recorded=> you are starting to get flooded!

                #Initialise Flooded_time
                Flooded_time = Timestep * step_fraction
                Flooded_height= a2/2

                #Record and reset Dry_time
                Dry_time = Dry_time + Timestep * (1-step_fraction)
                DRY_TIME.append(Dry_time); Dry_time=0

            elif a1*a2<0 and a1>0: #Meaning you cross the threshold from above when actual data is recorded=> you are drying up!

                #Initialise Dry_time
                Dry_time = Timestep * (step_fraction)

                #Record and reset Flooded_time  
                Flooded_time = Flooded_time + Timestep * (1-step_fraction)
                Flooded_height= Flooded_height + a1/2

                FLOODED_TIME.append(Flooded_time); Flooded_time=0
                FLOODED_HEIGHT.append(Flooded_height); Flooded_height=0

            elif a1*a2>0 and a1>0: #This is the flooded case
                Flooded_time = Flooded_time + Timestep
                Flooded_height = Flooded_height + (a1+a2)/2

            elif a1*a2>0 and a1<0: #This is the dry case
                Dry_time = Dry_time + Timestep
                

        DRY_TIME_Mean=np.mean(DRY_TIME)/(3600)
        DRY_TIME_Stdev=np.std(DRY_TIME)/(3600)
        FLOODED_TIME_Mean=np.mean(FLOODED_TIME)/(3600)
        FLOODED_TIME_Stdev=np.std(FLOODED_TIME)/(3600)
        FLOODED_HEIGHT=np.mean(FLOODED_HEIGHT)/FLOODED_TIME # this is the flooded height divided by the flooded time

        Metric3[1].append(DRY_TIME_Mean)
        Metric3[2].append(DRY_TIME_Stdev)
        Metric3[3].append(FLOODED_TIME_Mean)
        Metric3[4].append(FLOODED_TIME_Stdev)
        Metric3[5].append(FLOODED_HEIGHT)
        Metric3[6].append(FLOODED_HEIGHT)

    METRIC30=Tide_value[:-1]
    METRIC31=Metric3[1]
    METRIC32=Metric3[2]
    METRIC33=Metric3[3]
    METRIC34=Metric3[4]
    METRIC35=Metric3[5]
    METRIC36=Metric3[6]
    
    #print len(METRIC30), len(METRIC31), len(METRIC32), len(METRIC33), len(METRIC34), len(METRIC35),len(METRIC36)
    
    METRIC3=np.row_stack((METRIC30, METRIC31, METRIC32, METRIC33, METRIC34, METRIC35, METRIC36))
    
    
    
    return METRIC3

