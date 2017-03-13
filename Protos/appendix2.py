#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov  6 18:05:41 2016

TITLE
-----
Appendix 2: Tables of special 


DESCRIPTION
-----------
The appendix2 module assembles a group of constants and functions 
used in the theory for better comprehension and use of the mathematical 
expresions derived from the theoretical formalism 
of the synchrotron radiation of relativistic electrons.


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



#------------------PHYSICAL CONSTANTS------------------
# Electron rest mass in cgs units
m = scipy.constants.m_e*1e3     # g 

# Elementary charge in cgs units
e = scipy.constants.e*10*scipy.constants.c       # statC

# Light speed in cgs units
c = scipy.constants.c*1e2       # cm/s


#-------------------SPECIAL FUNCTIONS-----------------

#Function K
def function_K(x):
    """
    Special function K(x) = K_5/3 with K_5/3. 
    The second order modified Bessel function with imaginary arguments.
    
    :param x: the ratio between frequency nu and critical frequency nu_c.
    """
    K = scipy.special.kv(5/3,x)
    return K

    
#Function F
def function_F(x):
    """
    Special function F(x) = x int^inf_x[K_5/3(z) dz].
    Function F(x) derives from the modified bessel function and appears
    at expression 3.38 in Pacholczyk (1970), pp.89. It's purpose is
    to facilitate the understanding of the mathematical development.
    F(x) attains maximum at x=0.29.
    
    :param x: the ratio between frequency nu and critical frequency nu_c.
    """
    integralF = scipy.integrate.quad(lambda x: scipy.special.kv(5/3,x),
                                     x, +np.inf)
    
    F = x * integralF[0]
    return F
    

#Function S
def function_S(x):
    """
    Special function S(x) = F(x)/K(x).
    It appears at expression 3.46 (source function) and its purpose is to comprehend a part of
    the equation in order to ease the mathematical understanding.
    
    :param x: the ratio between frequency nu and critical frequency nu_c.
    """
    S = function_F(x)/function_K(x)
    
    return S


#Function Jx
def function_Jx(tauM,x):
    """
    Function J(tauM,x) = S(x)/S(0.29)*[1-exp(-tauM*(K(x)/K(0.29))).]
    It appears at expression 3.48 (transfer equation) and its purpose is to comprehend a part of
    the equation in order to ease the mathematical understanding.
    
    :param tauM:	The optical depth for frequency nu_M. nu_M=0.29*nu_c(E_0).
    :param x: the ratio between frequency nu and critical frequency nu_c.
    """
    Jx = function_S(x)/function_S(0.29) * (1-np.e**(-tauM*function_K(x)/function_K(0.29)))

    return Jx

#Function Jz
def function_Jz(z,p):
    """
    Function J(z,p) = z^5/2 * [1-exp(-z^(-(p+4)/2))].
    It appears at expression 3.53 (transfer equation of a power-law
    distribution) and its purpose is to comprehend a part of
    the equation in order to ease the mathematical understanding.
    
    Just a variable of the function. At expression 3.53, z=nu/nu_1 
    where nu_1 is defined by the condition tau(nu_1)=1.
    """
    Jz = z**(5/2) * (1-np.exp(-z**(-(p+4)/2)))

    return Jz
    
    
#Function I
def function_I(xM):
    """
    Function I(xM) = 1/xM * int_0^inf[z^2*exp(-z)F(xM/z^2)dz].    
    It appears at expressiona 3.54, 3.55 (emission and absorption coefficients
    for a relativistic Maxwellian distribution of electrons) and its
    purpose is to comprehend a part of the equation in order to ease 
    the mathematical understanding.
    
    For small values of xM (<1e-3) I(xM) can be expressed
    by an asymptotic formula.
    
    :param xM: the ratio between frequency nu and critical frequency nu_M. 
    """
    if xM > 1e-3:
        insideintI = lambda z: (z**2) * np.exp(-z) * function_F(xM/(z**2))
        integral_I=scipy.integrate.quad(insideintI, 0, +np.inf)
        I = integral_I[0]/xM
    else:
        I = 2.56*xM**(-2/3)
    
    return I

#Funcion G
def function_G(x):
    """
    Function G(x) = x * K_2/3(x).
    Function G(x) derives from the modified bessel function and appears
    at expression (3.38) for a better understanding of the mathematical
    development.
    
    :param x: the ratio between frequency nu and critical frequency nu_c.
    """
    G = x * scipy.special.kv(2/3,x)
    
    return G
    

#Function PI
def function_PI(x):
    """
    Function PI(x) = G(x)/F(x).
    Defined by expression 3.80 PI(x) states the degree of polarization for a
    totally thin homogeneous source.
    
    :param x: the ratio between frequency nu and critical frequency nu_c.
    """
    PI = function_G(x)/function_F(x)

    return PI    
    

#--------------------CONSTANTS----------------------------

# Constant c1
c1 = 3*e/(4*np.pi*(m**3)*(c**5))     # cgs: statC*s^5/(g^3*cm^5)


# Constant c2
c2 = 2*(e**4)/(3*(m**4)*(c**7))      # cgs: statC^4*s^7/(g^4*cm^7)
        

# Constant c3
c3 = np.sqrt(3)*e**3/(4*np.pi*m*c**2)    # cgs: statC^3*s^2/(g*cm^2)


# Constant c4
c4 = c1**(3/2) * c3 * (c**2)   # cgs: statC^9/2*s^(15/2)/(g^5*cm^(15/2))


# Constant c7
c7 = c1/(c2**2)        # cgs: 


#-------------------- c FUNCTIONS --------------------------

# c-function c5
def c5(p):
    """
    c-function c5 dependant on the power p.
    First appearance at 
    
    :param p: exponent of the function.
    """
    c5 = np.sqrt(3)/(16*np.pi) * (e**3)/(m*(c**2)) * (p+7/3)/(p+1) * scipy.special.gamma((3*p-1)/12) * scipy.special.gamma((3*p+7)/12)

    return c5

# c-function c6
def c6(p):
    """
    c-function c6 dependant on the power p.
    
    :param p: exponent of the function.
    """
    c6 = np.sqrt(3)*np.pi/72 *e *(m**5) *(c**10) *(p+10/3) * scipy.special.gamma((3*p+2)/12) * scipy.special.gamma((3*p+10)/12)
    
    return c6

    
# c-function c8
def c8(p):
    """
    c-function c8(p).
    
    :param p: exponent of the function.
    """
    insidec8 = lambda x: x**((p-3)/2) *function_F(x)
    c8 = scipy.integrate.quad(insidec8,0, +np.inf)
    
    return c8[0]
#HAY UN ERROR NO DESDEÑABLE (crece proporcionalmente a p)!!

# c-function c9
def c9(p):
    """
    c-function c9(p)
    
    :param p: exponent of the function.
    """
    c9 = np.sqrt(np.pi)/2 *scipy.special.gamma((p+5)/4) /scipy.special.gamma((p+7)/4)
    
    return c9
    

# c-function c10
def c10(p):
    """
    c-function c10(p).
    
    :param p: exponent of the function.
    """
    insidec10 = lambda x: x**(2/3*(p-1)) *function_F(x)
    c10 = scipy.integrate.quad(insidec10,0, +np.inf)
    
    return c10[0]
    #NO SALE 

# c-function c11
def c11(p):
    """
    c-function c11(p).
    
    :param p: exponent of the function.
    """

    insidec11 = lambda x: (1-x**(3/2))**(p-2) *x**((p+3)/2)
    c11 = scipy.integrate.quad(insidec11,0,1)
    
    return c11[0]

# c-function c14
def c14(p):
    """
    c-function c14(p)
    
    :param p: exponent of the function.
    """
    c14 = np.sqrt(6/np.pi)/4 *np.sqrt(e/c/m**3) *(p+1)*(3*p+10)/(3*p+7) *scipy.special.gamma((3*p+2)/12) *scipy.special.gamma((3*p+10)/12) /scipy.special.gamma((3*p-1)/12) /scipy.special.gamma((3*p+7)/12)
    
    return c14
#NO SALEE!!

#------------- ++ c-FUNCTIONS-----------------------------

# c-function c12
def c12(alpha,nu1,nu2):
    """
    c-function c12(p,nu)
    
    :param alpha:	Spectral index alpha=(p-1)/2. Exceptions for alpha=0.5, 1.
    :param nu1: frequency nu1.
    :param nu2: frequency nu2.
    """
    c12 = np.sqrt(c1)/c2 *(2*alpha-2)/(2*alpha-1) *(nu1**((1-2*alpha)/2)-nu2**((1-2*alpha)/2)) /(nu1**(1-alpha)-nu2**(1-alpha))

    return c12


#c-function c13
def c13(alpha, nu1, nu2):
    """
    c-function c13(alpha, nu1, nu2)
    
    :param alpha:	Spectral index alpha=(p-1)/2. Exceptions for alpha=0.5, 1.
    :param nu1: frequency nu1.
    :param nu2: frequency nu2.    
    """
    c13 = 0.921 *c12(alpha,nu1,nu2)**(4/7)
    
    return c13



#-----------------------TABLES----------------------------

#Table 1
# Values of variable x in Table 1
def table_1():
    """
    From Pacholczyk (1970), Appendix 2 - extended Table 1 (pp. 221-225).
    """
    t1_xvalues = np.concatenate((np.linspace(1e-4,1e-1,50),np.linspace(1.01e-1,9.99e-1,100),np.linspace(1,10,200),np.linspace(1e1,5e1,50)))
    col_F = np.vectorize(function_F)
    col_K = np.vectorize(function_K)
    col_S = np.vectorize(function_S)
    table1 = np.column_stack((t1_xvalues, col_F(t1_xvalues), col_K(t1_xvalues), col_S(t1_xvalues)))
        
    return table1
    #NO SALE
    
def table_2():
    """
    From Pacholczyk (1970), Appendix 2 - extended Table 2 (pp. 226).
    """
    t2_xvalues = np.concatenate((np.linspace(1e-4,1e-1,50),np.linspace(1.01e-1,9.99e-1,80),np.linspace(1,4.5,50)))
    col_J1 = np.vectorize(function_Jx)
    table2 = np.column_stack((t2_xvalues,col_J1(t2_xvalues,0.01), col_J1(t2_xvalues,0.1), col_J1(t2_xvalues,1),col_J1(t2_xvalues,10),col_J1(t2_xvalues,100)))
    
    return table2
    

#-----------------------PLOTS------------------------------

def figure_3p4():
    """
    Graphical representation of F(x) and G(x). Pacholczyk (1970), pp. 90.
    """
    x = np.concatenate((np.linspace(1e-4,1e-1,50),np.linspace(1.01e-1,9.99e-1,100),np.linspace(1,10,200),np.linspace(1e1,5e1,50)))
    F = np.vectorize(function_F)
    G = np.vectorize(function_G)
    
    plt.loglog(x,F(x),'ro',x,G(x),'bs', markersize=4)
    plt.ylim(1e-4,1e2)
    plt.xlim(1e-3,1e1)
    plt.text(4e-2,2,'F(x)')
    plt.text(1e-1,1e-1,'G(x)')
    plt.xlabel('$x$')
    
    
def figure_3p6():
    """
    Graphical representation of K(x), F(x) and S(x). Pacholczyk (1970), pp. 94.
    """
    x = np.concatenate((np.linspace(1e-4,1e-1,50),np.linspace(1.01e-1,9.99e-1,100),np.linspace(1,10,200),np.linspace(1e1,5e1,50)))
    F = np.vectorize(function_F)
    K = np.vectorize(function_K)
    S = np.vectorize(function_S)
    
    plt.loglog(x,F(x),'ro',x,K(x),'bs',x,S(x),'g^',markersize=4)
    plt.ylim(1e-5,1e4)
    plt.xlim(0.01,10)
    plt.text(0.1,2e2,'K(x)')
    plt.text(2e-2,2,'F(x)')
    plt.text(0.2,8e-3,'S(x)')
    plt.xlabel('$x$')
    
    
def figure_3p7():
    """
    Graphical representation of J(x,tauM) for different optical depths.Pacholczyk (1970), pp. 96.
    """
    x = np.concatenate((np.linspace(1e-4,1e-1,100),np.linspace(1.01e-1,9.99e-1,250),np.linspace(1,5e1,250)))
    Jx = np.vectorize(function_Jx)
    
    plt.loglog(x,Jx(0.01,x),'ko',x,Jx(0.1,x),'rs',x,Jx(1,x),'bo',markersize=4)
    plt.loglog(x,Jx(10,x),'gs',x,Jx(100,x),'c^',markersize=4)
    plt.ylim(1e-10,1e2)
    plt.xlim(1e-4,1e2)
    plt.text(6e-3,3e-3,'$10^{-2}$')
    plt.text(2e-2,3e-2,'$10^{-1}$')
    plt.text(1e-1,5e-1,'$10^0$')
    plt.text(4e-1,5,'$10^1$')
    plt.text(3.2e0,2e1,r'$\tau_M = 10^2$')
    plt.xlabel('$x$')
    
    
def figure_3p8():
    """
    Graphical representation of J(z,p). Pacholczyk (1970), pp. 97.
    """
    z = np.concatenate((np.linspace(1e-2,1e-1,10),np.linspace(1.01e-1,9.99e-1,200),np.linspace(1,10,250),np.linspace(10,1e3,200)))
    Jz = np.vectorize(function_Jz)
    
    plt.loglog(z,Jz(z,1),'ko',z,Jz(z,1.5),'rs',z,Jz(z,2),'b^', markersize=4)
    plt.loglog(z,Jz(z,2.5),'k^',z,Jz(z,3),'bs',z,Jz(z,3),'ro',markersize=4)
    plt.loglog(z,Jz(z,4),'ks',z,Jz(z,5),'bo',z,Jz(z,6),'r^',markersize=4)
    plt.ylim(1e-3,1)
    plt.xlim(1e-1,1e2)
    plt.text(0.9e1,1.8e-3,'p=6.0')
    plt.text(1e1,6e-3,'5.0')
    plt.text(1e1,2e-2,'4.0')
    plt.text(1.1e1,7e-2,'3.0')
    plt.text(1.05e1,1.25e-1,'2.5')
    plt.text(1.05e1,2.5e-1,'2.0')
    plt.text(1.05e1,5e-1,'1.5')
    plt.text(7e1,6e-1,'1.0')
    plt.ylabel('$J(z,p)$')
    plt.xlabel('$z$')
    
    
def figure_3p9():
    """
    Graphical representation of I(x_M). Pacholczyk (1970), pp. 99. 
    """
    xM = np.linspace(1e-6,5e3,300)
    I = np.vectorize(function_I)
    
    plt .loglog(xM,I(xM),'bo',markersize=4)
    plt.ylim(1e-10,5e2)
    plt.xlim(5e-5,1e4)
    plt.xlabel('$x_M$')
    plt.ylabel('$I(x_M)$')
    
    
def figure_3p14():
    """
    Graphical representation of the degree of polarization -PI(x)- for
    an optically thin source of synchrotron radiation. Figure 3.14 in 
    Pacholczyk (1970), pp. 106.
    """
    x = np.linspace(0.001,100,500)
    PI = np.vectorize(function_PI)
    
    plt.semilogx(x,PI(x),'ro',markersize=4)
    plt.ylim(0.5,1)
    plt.xlim(1e-3,1e2)
    plt.ylabel('$\Pi(x)$')
    plt.xlabel('$x$')
    

    

    
    
    
    

