import numpy as np
import ctypes
import os
from numpy.ctypeslib import ndpointer

_model = ctypes.CDLL(os.path.join(os.path.dirname(__file__),
                                  [f for f in os.listdir(os.path.dirname(__file__))
                                   if f.startswith("func_model")][0]))
_model.jetmodel.argtypes = [ctypes.c_char_p,
                            ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                            ndpointer(ctypes.c_int, flags="C_CONTIGUOUS"),
                            ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                            ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                            ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                            ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]


def freqtoeV(freq):
    '''Fn to convert frequency of light to energy in eV'''
    return 6.63E-34*freq/1.6E-19


def run_ssc(params, nblocks=1, seed=42, rand_gamma=0):
    '''
    Input: [W_j, E_max, alpha, theta_open_p, gamma_bulk, B0, theta_obs, A_eq, E_min]
    Output: IC_stokes, S_stokes,
    '''
    params = np.squeeze(np.array(params))
    IC_stokes, S_stokes, fIC, fS = np.zeros(
        50*3), np.zeros(50*3), np.empty(50), np.empty(50)
    nrings = [0, 1, 2, 3, 4, 5, 6][[1, 7, 19, 37, 61, 91, 127].index(nblocks)]
    blocks = np.array([nblocks, nrings, seed, rand_gamma], dtype=np.int32)
    fg_file = bytes(os.path.join(os.path.dirname(
        __file__), "..", "src", "FG.txt"), "utf-8")

    _model.jetmodel(fg_file, params, blocks, IC_stokes, S_stokes, fIC, fS)
    S_stokes = S_stokes.reshape((3, 50))
    IC_stokes = IC_stokes.reshape((3, 50))

    beta_bulk = (1.0-(params[4]**(-2.0)))**(0.5)
    doppler_factor = 1.0 / \
        (params[4]*(1.0-beta_bulk*np.cos(np.deg2rad(params[6]))))
    fIC *= doppler_factor
    fS *= doppler_factor

    return np.vstack([freqtoeV(fS), S_stokes]), np.vstack([freqtoeV(fIC), IC_stokes])
