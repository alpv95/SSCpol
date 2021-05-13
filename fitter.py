import numpy as np
import ctypes
from jet_fns import *
import multiprocessing as mp
from scipydirect import minimize
import os
import sys
import argparse
# from ipopt import minimize_ipopt
import pyswarms as ps
from numpy.ctypeslib import ndpointer


parser = argparse.ArgumentParser()
parser.add_argument('--method', type=str, choices=["ip","direct", "ps"], default="ip",
                    help='Whether to resume from a different run')
parser.add_argument('--blazar', type=str, choices=["J2011","TXS", "S5low", "S5flare"], default="J2011",
                    help='Which blazar to fit')
args = parser.parse_args()


_model = ctypes.CDLL('func_model.so')
_model.jetmodel.argtypes = [ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                          ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"),
                         ndpointer(ctypes.c_double, flags="C_CONTIGUOUS")]

if args.blazar == "J2011":
    d_Blazar = 1000E6*3.08E18 
    z = 0.2
    data = np.loadtxt('data/new_data_sed_CGRaBSJ0211+1051_XMM.txt')
    data2 = np.loadtxt('data/new_data_sed_CGRaBSJ0211+1051.txt')
    data = np.concatenate([data[6:],data2[12:]], axis=0)
elif args.blazar == "S5low":
    d_Blazar = 1627E6*3.08E18 
    z = 0.31
    data = np.loadtxt('data/S5low.txt')
    data = np.insert(data, 1, 0, axis=1)
elif args.blazar == "S5flare":
    d_Blazar = 1627E6*3.08E18 
    z = 0.31
    data = np.loadtxt('data/S5flare.txt')
else:
    d_Blazar = 1789E6*3.08E18 
    z = 0.3365
    data = np.loadtxt('data/txs.txt')
    data = np.insert(data, 1, 0, axis=1)

def run_ssc(params,):
    '''
    Input: [W_j, E_max, alpha, theta_open_p, gamma_bulk, B0, theta_obs, A_eq]
    Output: IC_stokes, S_stokes,
    '''
    params = np.squeeze(np.array(params))     
    IC_stokes, S_stokes, fIC, fS =  np.zeros(50*3), np.zeros(50*3), np.empty(50), np.empty(50)
    _model.jetmodel(params, IC_stokes, S_stokes, fIC, fS)
    S_stokes = S_stokes.reshape((3,50))
    IC_stokes = IC_stokes.reshape((3,50))

    beta_bulk=(1.0-(params[4]**(-2.0)))**(0.5)
    doppler_factor = 1.0/(params[4]*(1.0-beta_bulk*np.cos(np.deg2rad(params[6]))))
    fIC *= doppler_factor
    fS *= doppler_factor

    return np.vstack([freqtoeV(fS), S_stokes]), np.vstack([freqtoeV(fIC), IC_stokes])



def loglike(params):
    params = np.squeeze(np.array(params))
    params[0] = 10**params[0]
    params[1] = 10**params[1]
    params[5] = 10**params[5]
    IC_stokes, S_stokes, fIC, fS =  np.zeros(50*3), np.zeros(50*3), np.empty(50), np.empty(50)
    _model.jetmodel(params, IC_stokes, S_stokes, fIC, fS)
    S_stokes = S_stokes.reshape((3,50))
    IC_stokes = IC_stokes.reshape((3,50))

    beta_bulk=(1.0-(params[4]**(-2.0)))**(0.5)
    doppler_factor = 1.0/(params[4]*(1.0-beta_bulk*np.cos(np.deg2rad(params[6]))))
    fIC *= doppler_factor
    fS *= doppler_factor

    dist_factor = 1.0E7*(1.0/((4.0*np.pi*d_Blazar**2.0)*(1.0+z)**2.0)) 
    loss = np.sum((data[:,2] - np.log10(np.interp(10**data[:,0], freqtoeV(fIC), IC_stokes[0,:] * dist_factor) 
                                        + np.interp(10**data[:,0], freqtoeV(fS), S_stokes[0,:] * dist_factor)) )**2 )
    if args.method == 'ps':
        print('LOSS: ', loss, "   ", params[2], '\n')
        return np.array([loss])
    print('LOSS: ', loss, '\n')
    return loss
    

def main():
    print("BLAZAR:  ", args.blazar)
    bounds = [(36.3, 39),(10,13),(1.65,2.05),(15,75),(7,20),(-4.7,-3.1),(0.1,6),(0.9,1.4)]
    print("bounds: ", bounds)
    
    x0 = [1.3E37,1.7E11,1.85,40.0,16,1E-4,1.5,1]

    if args.method == "ip":
        res = minimize_ipopt(lambda cube: loglike(cube,active=active), x0=x0, bounds=bounds, tol=5e-2)
    elif args.method == "direct":
        res = minimize(loglike, bounds, fglobal=1.0, fglper=(0.2*100), disp=True, maxt=10000)
    elif args.method == "ps":
        options = {'c1': 0.5, 'c2': 0.3, 'w':0.9}
        # Call instance of GlobalBestPSO
        optimizer = ps.single.GlobalBestPSO(n_particles=15, dimensions=len(bounds), bounds=(np.array([b[0] for b in bounds]),np.array([b[1] for b in bounds])),
                                    options=options)
        # Perform optimization
        res = optimizer.optimize(loglike, iters=10000, n_processes=15)
    print(res)



if __name__ == "__main__":
    main()
