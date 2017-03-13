#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 10:04:46 2017

@author: mglez
"""

import scipy
import scipy.constants
import numpy as np
import matplotlib.pyplot as plt
import lmfit

import appendix2
from appendix2 import function_F


"""
A little simulation to calculate the Riemann integral of function F(x) from
1e-4 to 50.
This result will reproduce the values shown at A2.Table 1 (Pacholzcyk 1970)
and serves to proof that c10(1) gives a wrong result on A2.Table 7.


"""
x = np.linspace(1e-4,50,1e3)
suma = 0

for i in range(0,(len(x)-1)):
    #Base times the average height between two consecutive points
    area = (x[i+1]-x[i])*(function_F(x[i+1])+function_F(x[i]))/2
    suma = suma + area

print(suma)


"""
For more evidence we calculate from 10-6 to 1e-4

"""
y = np.linspace(1e-6,1e-4,1e3)
suma2 = 0

for j in range(0,(len(y)-1)):
    #Base times the average height between two consecutive points
    area2 = (y[j+1]-y[j])*(function_F(y[j+1])+function_F(y[j]))/2
    suma2 = suma2 + area2

print(suma2)

"""
Below x = 1e-04 the influence of the integral of F(x) is neglectable
"""


    
