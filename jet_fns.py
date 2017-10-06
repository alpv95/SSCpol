from __future__ import division
import numpy as np

def freqtoeV(freq):
    '''Fn to convert frequency of light to energy in eV'''
    return 6.63E-34*freq/1.6E-19
