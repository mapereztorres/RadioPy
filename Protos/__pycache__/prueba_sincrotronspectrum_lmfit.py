#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 06:06:37 2017

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
c2 = appendix2.c2
c3 = appendix2.c3



# function
def KP(t, E, H, pitch, p, x):
    E_T = 1/(c2*H**2*(sin(pitch))**2*E*t)
    angular_integral = scipy.integrate.quad(1/((sin(pitch))**2),pitch,np.pi/2.)
    energy_integral = scipy.integrate.quad(1/(E**(-p)*appendix2.function_F(x)), 0, E_T)
    
    
    
    KP = 4*np.pi*c3 *angular_integral*energy_integral*(1-c2*H**2*sin(pitch)**2*E*t)**(p-2)
    return KP
    
modelKP = lmfit.Model(KP)
modelKP.param_names



# Model of original function
#from lmfit import Model
#>>> gmod = Model(gaussian)
#>>> gmod.param_names
#set(['amp', 'wid', 'cen'])
#>>> gmod.independent_vars)
#['x']

