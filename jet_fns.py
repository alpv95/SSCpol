from __future__ import division
import numpy as np

def freqtoeV(freq):
    '''Fn to convert frequency of light to energy in eV'''
    return 6.63E-34*freq/1.6E-19

def R_0(E_j, gamma, B0, A_eq = 1):
    '''Calculate starting radius of jet given, B_0, W_j, gamma_bulk and A_eq'''
    sqrt = (2.0*E_j*A_eq*4.0*np.pi*10.0**-7)/(gamma*gamma*(np.pi*3.0E8*B0*B0)*(1.0+A_eq))
    ans = sqrt**0.5
    return ans
