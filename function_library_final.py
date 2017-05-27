"""
Created 19/01/2016
This file defines a series of functions that are useful to run Andrea's code.
@author Guillaume Goodwin
"""
import numpy as np
from scipy import optimize
import cmath
import math
from scipy.optimize import curve_fit# Import the curve fitting module
from scipy.optimize import leastsq# Import the curve fitting module
import simple_distributions as sd
import peak_detect as pd
import scipy.stats as stats
import matplotlib.pyplot as plt

#-------------------------------------------------------------------
#1.
#The Newton-Raphson method

def NewtonRaphson (k, deltak, deltak_min, depth, omega, g):
    #the Newton-Raphson method is uded to iteratively determine an
    #intercept: https://en.wikipedia.org/wiki/Newton%27s_method
    #Here it is written in a general form

    #k is the x-coordinate of the intercept: That's what we want to know
    #deltak is our indicator of how close we are to the solution
    #Fk is the the function for which we are looking for an intercept
    #dFk is its derivative with respect to k
    #deltak_min is our accuracy objective for deltak
    
    while abs(deltak)>deltak_min:
                Fk=np.sqrt(g*k*np.tanh(k*depth))-omega #Lamb's equation expressed as f(k)
                dFk=0.5*np.sqrt(g)*(np.tanh(k*depth)-k*depth*(np.tanh(k*depth))**2+k*depth)/(np.sqrt(k*np.tanh(k*depth))) #f'(k)
                deltak=-Fk/dFk
                k=k+deltak
                deltak=deltak/k   
    k_final=k
    return k_final

#-------------------------------------------------------------------
#1.5. Bisection method
def Bisection_k (kmin,kmax,deltak_min, depth, omega, g):
    iteration=0
    knext=0
    while abs(kmin-kmax)>deltak_min and iteration<1000:
        iteration=iteration+1
        Fk_0=np.sqrt(g*kmin*np.tanh(kmin*depth))-omega #Lamb's equation expressed as f(k)
        Fk_1=np.sqrt(g*kmax*np.tanh(kmax*depth))-omega #Lamb's equation expressed as f(k)

        if Fk_0*Fk_1<0: 
            knext=(kmin+kmax)/2
            Fk_0=np.sqrt(g*knext*np.tanh(knext*depth))-omega

            if Fk_0*Fk_1<0:
                kmin=knext

            else:
                kmax=knext

        k_final=knext

    return k_final





#-------------------------------------------------------------------
#2.
#Accurate and iterative calculation of the terms of wave motion.
#This function requires the existence of k defined by function 1.

def wave_energy(Emin,Emax,acc,windspeed,depth,k,g,Twave,omega,RHOw,RHOs,gammaPM,cd,cbf):   
    #The objective of this function is to express wave energy E as
    #E(n+1)=f(E(n)) and then x=f(x) and to find its limit/solution L
    #This is done using the Bisection method.
    
    #Note that we use the equations described in Andrea's code, with a twist.

    #1. Definition of initial parameters
    alpha=80*((RHOs*cd*windspeed**2)/(g*k))**2*omega**2 #alpha coefficient for the wind-generation term
    beta=5*RHOs/Twave*((windspeed*k/omega)-0.9) #beta coefficient for the wind-generation term
    HWmax=0.78*depth #Maximum wave height
    
    #2. Expression of the wave motion conservation equation
    #We use a function to make it displaceable
    
 
    def Conservation(E,alpha,beta,k,cbf,omega,gammaPM,Twave,HWmax,g,RHOw):
        
        #Relationship between wave action N and wave energy E
        N=E/omega
        H=np.sqrt(8*N*omega/(RHOw*g))
        
        #Wind Wave Generation source term
        #Swind=alpha+beta*E
        
        #Wave Bottom friction
        #Sbf=Wbf*E
        Wbf=-4*cbf*H*np.pi*k/(Twave*np.sinh(k*depth)*np.sinh(2*k*depth))
        
        
        #Wave White capping
        #Swc=Wwc*E
        gammawc=H**2*omega**4/(8*g**2) #Integral wave steepness parameter
        Wwc=-3.33*10**(-5)*omega*(gammawc/gammaPM)**2 
        
        
        #Wave Breaking  
        #First we calculate Qbr, the breaking probability
        if (H/HWmax)>=1:
            Qbr=1
        elif (H/HWmax)>=0.4:
            Qbr=1.8*(H/HWmax-0.4)**3+1.7*(H/HWmax-0.4)**2  
        else:
            Qbr=0
        
        #Then we calculate the corresponding source term
        #Sbr=Qbr*HWmax**2*RHOw*g/(np.sqrt(2)*Twave)
        #Or Sbr=Wbr*E
        Wbr=-2*Qbr/Twave*(H/HWmax)**2

        #Conservation equation
        #Cons=Sbf+Swc+Sbr-Swind
        Cons=alpha+(beta+Wbf+Wwc+Wbr)*N
        
        
        return Cons
    
    
    #3. Initialisation of our iteration
    iteration=0
    

    #4. Iterative calculation of E using the Bisection method
    while abs(Emin-Emax)>10**(-acc) and iteration<1000:
        iteration=iteration+1
        Cons_0=Conservation(Emin,alpha,beta,k,cbf,omega,gammaPM,Twave,HWmax,g,RHOw)
        Cons_1=Conservation(Emax,alpha,beta,k,cbf,omega,gammaPM,Twave,HWmax,g,RHOw)
        
        if Cons_0*Cons_1<0:      
            Enext=(Emin+Emax)/2
            Cons_0=Conservation(Enext,alpha,beta,k,cbf,omega,gammaPM,Twave,HWmax,g,RHOw)
            
            if Cons_0*Cons_1<0:
                Emin=Enext
            #elif Cons_0*Cons_1==0:
                #print 'found it', Enext, Cons_0
            else:
                Emax=Enext
                
        #else:
            #print windspeed, depth,'There is no solution?'
            #break

    E_final=Enext
   
    return E_final



#-------------------------------------------------------------------
#2.5 Calculation of Wave Energy Source Terms (WEST) after (Mariotti & Fagherazzi, 2010)

def WEST(E,windspeed,depth,k,g,Twave,sigma,RHOw,RHOs,gammaPM,cd,cbf):  
    
    H=np.sqrt((8*E)/(RHOw*g))#Wave height
    HWmax=0.78*depth #Maximum wave height

    #1. Wind generation term (generation)
    alpha=80*((RHOs*cd*windspeed**2)/(g*k))**2*sigma**2 
    
    if E>0:
        beta=5*RHOs/Twave*((windspeed*k/sigma)-0.9) 
        Sw=alpha+beta*E

        #2. Bottom friction term (decay)
        Wbf=-4*cbf*H*np.pi*k/(Twave*np.sinh(k*depth)*np.sinh(2*k*depth))
        Sbf=Wbf*E



        #3. Whitecapping term (decay)
        gammawc=H**2*sigma**4/(8*g**2) #Integral wave steepness parameter
        Wwc=-3.33*10**(-5)*sigma*(gammawc/gammaPM)**2 
        Swc=Wwc*E

        #4. Breaking term (decay)
        #First we calculate Qbr, the breaking probability
        if (H/HWmax)>=1:
            Qbr=1
        elif (H/HWmax)>=0.4:
            Qbr=1.8*(H/HWmax-0.4)**3+1.7*(H/HWmax-0.4)**2  
        else:
            Qbr=0

        #Then we calculate the corresponding source term
        Wbr=-2*Qbr/Twave*(H/HWmax)**2
        Sbr=Wbr*E

    else:
        Sw=alpha
        Sbf=0
        Swc=0
        Sbr=0
    
    return Sw,Sbf,Swc,Sbr




#-------------------------------------------------------------------
#3.
#Calculation of the bed shear stress tau in function of water depth (this is for a given wind speed)

def shear_stress(H,depth,k,Twave,RHOw,ni,D50):
# depth is the water depth [m]
# k is the wave number

#these variables are in fact model parameters, and no not need to be defined
# Twave is the the wave period in [s]
# RHOw is the volumetric density ofwater in [kg/m**3]
# ni is the kinematic viscosity of water
#D50 is the average grain size

    #Maximum horizontal orbital velocity at the bottom     
    um=np.pi*H/(Twave*np.sinh(k*depth))
            
    if um <> 0:
        #Wave friction factor
        #Rough wave
        fwave_s=1.39*(um*Twave/(2*np.pi*D50/12))**(-0.52)
        #Smooth wave
        Rewaves=um**2*Twave/(2*np.pi*ni) #wave Reynolds number
        if (Rewaves<=500000):
            fwave_l=2*Rewaves**(-0.5)
            regime=0
        else:
            fwave_l=0.0521*Rewaves**(-0.187)

        fwave=max(fwave_s,fwave_l)
        tau=0.5*RHOw*fwave*um**2
        tau_ST=0.5*RHOw*fwave_l*um**2

    else:
        tau=0
        tau_ST=0
        fwave=0
            
    return tau,tau_ST,fwave

 
    
#-------------------------------------------------------------------
#4.  
# Interpolation

#This function calculates the values of the bottom shear stress over a tide cycle
#Because the step chosen for D is smaller than the step of depth categories
#We must interpolate the values of tau linearly

def interpolate(X, depth, d_depth, Y):
    
    #These are the values of depth neighbouring D
    #This extra step compensates the fact that python seems to have trouble dividing by 0.1...
    A=1/d_depth
    Xlow=d_depth*np.floor(X*A)
    Xhigh=d_depth*np.ceil(X*A)
    
    #x-distance
    #This is the step between Dhigh and Dlow
    Xstep=Xhigh-Xlow
    
    #x-index-distance
    #these are the indices of Dlow/Dhigh in depth
    XlowIndex=depth.searchsorted(Xlow)
    XhighIndex=depth.searchsorted(Xhigh)
    
    #these are the associated values of Y
    Ylow=Y[XlowIndex]
    Yhigh=Y[XhighIndex]

    #And here is the difference between the two
    Ystep=Yhigh-Ylow

    #Now we calculate the slope, only where Dstep is not nil(how?)
    Slope=Ystep/Xstep
    Ynew=Ylow+Slope*(X-Xlow)

    #now we deal with the impossible values
    A=np.where(np.isnan(Ynew))
    Ynew[A]=Ylow[A]


    return Ynew

    
#-------------------------------------------------------------------
#5. Fitting an exponantial to wind data
def fit_exponential(hist, bins, threshold_velocity):
    # This takes the left-side value of all bins except the first one.
    # We create this variable to conceptualize the pdf x-axis:
    # it starts at 1 because the probability of a number to be smaller than 0 is useless.
    x_hist = bins[1:]
    #print "bin edges are: " + str(x_hist)
    
    
    #We now transform the histogram (already a density) into a probability density function by...
    #Dividing it by the max of the cumulative sum, which is 1. But, you know, just in case...
    y_pdf = hist/hist.cumsum().max()
    
    # now get the cdf
    y_cdf = hist.cumsum()/hist.cumsum().max()  # Normalise the cumulative sum
    
    # now get a truncated version of the data, cutting off any values less
    # than the threshold   
    truncated_x = []
    truncated_pdf = []
    truncated_cdf = []    
    for x,pdf_value,cdf_value in zip(x_hist,y_pdf,y_cdf):
        if x >= threshold_velocity:
            truncated_x.append(x)
            truncated_pdf.append(pdf_value)
            truncated_cdf.append(cdf_value)
            
    # turn the data into numpy arrays, necessary for the fitting routine below
    trunc_x = np.asarray(truncated_x)
    #print 'trunc_x', trunc_x
    trunc_pdf = np.asarray(truncated_pdf)
    #print 'trunc_pdf', trunc_pdf

    # fit the pdf for the thresholded data. 
    popt_exp_pdf, pcov_exp_pdf = curve_fit(sd.exponential_fxn, trunc_x, trunc_pdf)
     
    # get the fitted pdf
    #print "The fitted decay coefficient is: "
    #print popt_exp_pdf  
    y_exp_fit = sd.exponential_fxn(trunc_x,popt_exp_pdf)        
    #print y_exp_fit
    
    RMSE =  sd.RMSE(truncated_pdf,y_exp_fit)
    #print "The RMSE is: " + str(RMSE)
    
    return x_hist,y_pdf,trunc_x,trunc_pdf,popt_exp_pdf,y_exp_fit,RMSE


#-------------------------------------------------------------------
#6. Fitting a double gaussian to erosion results

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.))) 

def double_gaussian(x,midpoint,sig,amplitude,spacing):
    mu = midpoint-spacing/2
    first_gauss = gaussian(x,mu,sig)
    mu2 = midpoint+spacing/2
    second_gauss = gaussian(x,mu2,sig)
    double_gauss = np.multiply(amplitude,(np.add(first_gauss,second_gauss)))
    return double_gauss


def fit_double_gaussian(depth,erate):
    # State the filename
    #filename = 'c:\\Users\\smudd\\Documents\\Papers\\Tidal_paper_padova\\Wind_Data\\s_andrea_ok.txt'    
    
    # First load the erosion data file. This data just has an erosion rate 
    # and a depth on every line
    
    
    # Find some peaks
    _max, _min = pd.peakdetect(erate, depth,11)# 20, 0.10)
    #print _max
    #print _min
    
    xm = [p[0] for p in _max]
    ym = [p[1] for p in _max]
    
    #print 'xm', xm
    #print 'ym', ym

    # default initial guess
    midpoint = 0
    sig = 0.1
    amplitude = 1
    spacing = 1
    
    # now construct the initial guess from these data
    n_peaks = len(xm)
    
    if n_peaks==0:#If we don't see any peaks, we continue with a null gaussian.
        midpoint = 0
        sig = 0.0
        amplitude = 0
        spacing = 0
    
    elif n_peaks == 1:
        amplitude = ym[0]
        spacing = 0
        midpoint = xm[0]
    else:
        amplitude = max(ym)
        spacing = xm[1]-xm[0]
        midpoint = 0.5*(xm[1]+xm[0])
    
    # new initial guess    
    initial_guess = [midpoint,sig,amplitude,spacing]
    
    #print "Initial guess: "
    #print "midpoint: " + str(midpoint)
    #print "sigma: " + str(sig)
    #print "amplitude: " + str(amplitude)
    #print "spacing: " + str(spacing)
        

    # try to get some reasonable inital guesses based on the data

    # now fit the erosion data to a double gaussian
    popt_dg, pcov_dg = curve_fit(sd.double_gaussian, depth, erate,initial_guess)
   
   
    # get the fitted pdf
    #print "The fitted components are: "
    #print popt_dg 
    y_dg_fit = sd.double_gaussian(depth,popt_dg[0],popt_dg[1],popt_dg[2],popt_dg[3])        
    
    
    RMSE =  sd.RMSE(erate,y_dg_fit)
    #print "The RMSE is: " + str(RMSE)
    
    return depth,erate,popt_dg,y_dg_fit,RMSE

    
    
    
#-------------------------------------------------------------------
#6.  
# Fitting a Weibull distribution

def weibull_pdf(x,k,lambda_wb,loc):
    #We manually determine the shape coefficient

    Wb=stats.exponweib.pdf(x, k, lambda_wb) #(Could be x-loc, better RMSE but lousy at the end)
    return Wb

def Weibull_fit(hist, bins):
    
    x_hist = bins[1:]
    print "bin edges are: " + str(x_hist)
    
    #We now transform the histogram (already a density) into a probability density function by...
    #Dividing it by the max of the cumulative sum, which is 1. But, you know, just in case...
    y_pdf = hist/hist.cumsum().max()
    
    # now get the cdf
    y_cdf = hist.cumsum()/hist.cumsum().max()  # Normalise the cumulative sum
    
    # fit the pdf to the data. 
    popt_w_pdf, pcov_w_pdf = curve_fit(weibull_pdf, x_hist, y_pdf)
     
    # get the fitted pdf
    print "The fitted Coefficients are: "
    print popt_w_pdf  
    y_w_fit = weibull_pdf(x_hist,popt_w_pdf[0],popt_w_pdf[1],popt_w_pdf[2])         
    
    RMSE =  sd.RMSE(y_pdf,y_w_fit)
    print "The RMSE is: " + str(RMSE)
    
    return x_hist,popt_w_pdf,y_w_fit,RMSE 
    
     
#-------------------------------------------------------------------
#7.  
# Some subaxes
#http://stackoverflow.com/questions/17458580/embedding-small-plots-inside-subplots-in-matplotlib
def add_subplot_axes(ax,rect,axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax
    
#-------------------------------------------------------------------
#8. Fitting a power or linear law to erosion results

def power(x, a, n, b):
    return a*np.power(x,n)+b

def fit_power(U,my_array):
    guess=[1,1,1]
    popt, pcov = curve_fit(power, U, my_array,guess)
    Y=power(U,popt[0],popt[1],popt[2])
    RMSE =  sd.RMSE(my_array,Y)

    
    return popt,Y,RMSE  
    
    
    
    
