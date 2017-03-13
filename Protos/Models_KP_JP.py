#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  6 18:05:41 2016

TITLE
-----
??


DESCRIPTION
-----------
?


METADATA
--------
:author: Miguel Glez San Emeterio
:organization: Universidad de Zaragoza
:contact: mglez.sanemeterio@gmail.com
:license: yet to be determined

"""
#Imported packages
import scipy
import scipy.constants
import numpy as np
import matplotlib.pyplot as plt

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
#Theta==pitch angle
#ET = 1/(c2*H**2*(np.sin(theta))**2*E*t)

#Kardashev-Pacholczyk radiative losses for a frequency \nu
#def model_KP(t):
    #inside_Itheta = lambda th: sin(th)**2
    #I_theta = scipy.integrate.quad(lambda th: sin(th)**2, 0, +np.inf)
    #I_energy = scipy.integrate.quad(insideintI, 0, +np.inf)
    #constants = 4*np.pi*c3*N_0*s*B
    
    #L_nu = constants*I_theta*I_energy
