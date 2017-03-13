#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 08:58:57 2017

@author: mglez
"""

import scipy
import scipy.constants
import numpy as np
import matplotlib.pyplot as plt
import lmfit

import appendix2

#------------------PHYSICAL CONSTANTS------------------
# Electron rest mass in cgs units
m = scipy.constants.m_e*1e3     # g 

# Elementary charge in cgs units
e = scipy.constants.e*10*scipy.constants.c       # statC

# Light speed in cgs units
c = scipy.constants.c*1e2       # cm/s

#------------------------------------------------------

#--------------RADIATIVE LOSSES MODELS-----------------

#Constants and values
c1 = appendix2.c1
c2 = appendix2.c2
c3 = appendix2.c3
c7 = appendix2.c7
#Theta==pitch 

#------------------------EXAMPLE------------------------
"""
Fig. 3 MurgiaAA-99. Fit the data shown at the graphical with the two models
(low and high frequencies).
"""
# Frequency data \nu [Ghz]
freq = [0.033000347911253
,0.071097094323124
,0.159846710643432
,0.296634883917773
,0.316227766016838
,0.375036963792269
,0.591864405186301
,0.681292069057961
,0.90658265844882
,1.31938143966415
,2.61015721568254
,4.84378864954104
,4.64158883361278
,7.74263682681127
,10.4356262977539
,14.6779926762207]

# Flux density data [mJy]
flux = [44478.2967612763
,24484.3674682223
,13768.5716485276
,9182.54283565628
,8079.92644838296
,7909.48393241786
,6528.52114112785
,5446.41977741133
,4641.58883361278
,2903.77500727351
,1598.46710643432
,938.041866639814
,807.992644838296
,435.400465365665
,330.003479112528
,234.622884814226
]

print('----------------EXPERIMENTAL DATA------------------')
plt.loglog(freq,flux, 'bo')
plt.xlabel(r'$\nu$ [GHz]')
plt.ylabel('Flux density [mJy]')
plt.show()


#------------------LOW FREQUENCY ----------------------
print('----------------LOW FREQUENCIES MODEL------------------')
#
def I_low(nu,N_0,s,p,H):
    """

    Source brightness or intensity per unit frequency I_$nu$(t) 
    in the case of an optcally thin medium for frequencies 
    much lower nu_T (Pacholzcyk 1970, eq. 6.27)
    
    :param nu: radiation frequency.
    :param N_0: initial number of electrons.
    :param s: extent of the source along the line of sight.
    :param p: exponent of the power-law.
    :param H: magnetic field.
    """

    Il = 2*np.pi*c3*c1**((p-1)/2)*nu**((1-p)/2)*appendix2.c8(9)*appendix2.c9(p)*H**((p-1)/2)*N_0*s
    
    
    return Il

lowFreqIntensity = lmfit.Model(I_low)
print(lowFreqIntensity.independent_vars)
print(lowFreqIntensity.param_names)

lowFreqResult = lowFreqIntensity.fit(flux, nu=freq, N_0=1e11, s=0.1, p=2.44, H=100)
print(lowFreqResult.fit_report())

try:
    plt.loglog(freq, flux,         'bo')
    plt.loglog(freq, lowFreqResult.init_fit, 'k--')
    plt.loglog(freq, lowFreqResult.best_fit, 'r-')
    plt.xlabel(r'$\nu$ [GHz]')
    plt.ylabel('Flux density [mJy]')
    plt.show()
except:
    pass




lowFreqParams = lowFreqResult.best_values
N_0_low = lowFreqParams['N_0']
s_low = lowFreqParams['s']
p_low = lowFreqParams['p']
H_low = lowFreqParams['H']




#--------------------------HIGH FREQUENCY----------------------
print('----------------HIGH FREQUENCIES MODEL------------------')
    
def I_high(nu,N_0,s,p,H,t):
    """
    Source brightness or intensity per unit frequency I_$nu$(t) 
    in the case of an optcally thin medium for frequencies 
    much greater than nu_T (Pacholzcyk 1970, eq. 6.28)
    
    :param nu: radiation frequency.
    :param N_0: initial number of electrons.
    :param s: extent of the source along the line of sight.
    :param p: exponent of the power-law.
    :param H: magnetic field.
    :param t: lifetime of the source since it began to behave as a power-law distribution of electrons.
    """
    Ih = 2*np.pi*c1**((2*p+1)/3)*c2**(-(p+5)/3)*nu**(-(2*p+1)/3)*t**(-(p+5)/3)*H**(-2)*N_0*s

    
    return Ih
    

highFreqIntensity = lmfit.Model(I_high)
print(highFreqIntensity.independent_vars)
print(highFreqIntensity.param_names)

highFreqResult = highFreqIntensity.fit(flux, nu=freq, N_0=1e8, s=0.1, p=1, H=100, t=1e8)
print(highFreqResult.fit_report())

#Calculation of the break frequency (Pacholzcyk 1970 eq. 6.26)
highFreqParams = highFreqResult.best_values
N_0_high = highFreqParams['N_0']
s_high = highFreqParams['s']
p_high = highFreqParams['p']
H_high = highFreqParams['H']
t_high = highFreqParams['t']

freq_break = c7*H_high**(-3)*t_high**(-2)
alpha_inj = (p_high-1)/2


# Plot of the data and the fitting
try:
    plt.loglog(freq, flux,         'bo')
    plt.loglog(freq, highFreqResult.init_fit, 'k--')
    plt.loglog(freq, highFreqResult.best_fit, 'g-')
    plt.xlabel(r'$\nu$ [GHz]')
    plt.ylabel('Flux density [mJy]')
    plt.show()
except:
    pass

print('The break frequency is ',freq_break,' GHz')
print('The injection exponent alpha is ',alpha_inj)
print('The lifetime of the source is ',t_high,' yrs')