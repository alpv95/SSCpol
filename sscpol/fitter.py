import numpy as np
from sscpol.jet_fns import *
import multiprocessing as mp
from scipydirect import minimize as minimize_direct
from scipy.optimize import minimize
# import argparse
import pyswarms as ps


class SSC_Fitter(object):
    def __init__(self, blazar, method="standard", nblocks=1, seed=42,):
        if blazar not in ["J2011", "TXS", "S5low", "S5flare"]:
            raise ValueError("blazar must be one of J2011, TXS, S5low, S5flare")
        if blazar == "J2011":
            self.d_Blazar = 1000E6*3.08E18 
            self.z = 0.2
            data = np.loadtxt('data/new_data_sed_CGRaBSJ0211+1051_XMM.txt')
            data2 = np.loadtxt('data/new_data_sed_CGRaBSJ0211+1051.txt')
            self.data = np.concatenate([data[6:],data2[12:]], axis=0)
        elif blazar == "S5low":
            self.d_Blazar = 1627E6*3.08E18 
            self.z = 0.31
            data = np.loadtxt('data/S5low.txt')
            self.data = np.insert(data, 1, 0, axis=1)
        elif blazar == "S5flare":
            self.d_Blazar = 1627E6*3.08E18 
            self.z = 0.31
            self.data = np.loadtxt('data/S5flare.txt')
        else:
            self.d_Blazar = 1789E6*3.08E18 
            self.z = 0.3365
            data = np.loadtxt('data/txs.txt')
            self.data = np.insert(data, 1, 0, axis=1)
        
        if method not in ["standard", "direct", "ps"]:
            raise ValueError("method must be one of ip, direct, ps")
        self.method = method
        self.nblocks = nblocks
        self.seed = seed

    def run(self,):
        bounds = [(36.3, 39),(10,13),(1.65,2.05),(15,75),(7,20),(-4.7,-3.1),(0.1,6),(0.9,1.4),(6,6.8)]
        print("bounds: ", bounds)
        x0 = [1.3E37,1.7E11,1.85,40.0,16,1E-4,1.5,1,6.7]
        
        if self.method == 'standard':
            res = minimize(self.loglike, x0, method='L-BFGS-B', bounds=bounds, options={'disp': True})
        elif self.method == 'ps':
            options = {'c1': 0.5, 'c2': 0.3, 'w':0.9}
            # Call instance of GlobalBestPSO
            optimizer = ps.single.GlobalBestPSO(n_particles=15, dimensions=len(bounds), bounds=(np.array([b[0] for b in bounds]),np.array([b[1] for b in bounds])),
                                        options=options)
            # Perform optimization
            res = optimizer.optimize(self.loglike, iters=10000, n_processes=15)
        elif self.method == "direct":
            res = minimize_direct(self.loglike, bounds, fglobal=1.0, fglper=(0.2*100), disp=True, maxt=10000)

        print("res: ", res)
        return res

    def loglike(self, params):
        params = np.squeeze(np.array(params))
        params[0] = 10**params[0]
        params[1] = 10**params[1]
        params[5] = 10**params[5]
        params[8] = 10**params[8]

        sync, ic = run_ssc(params, nblocks=self.nblocks, seed=self.seed)

        dist_factor = 1.0E7*(1.0/((4.0*np.pi*self.d_Blazar**2.0)*(1.0+self.z)**2.0)) 
        loss = np.sum((self.data[:,2] - np.log10(np.interp(10**self.data[:,0], ic[0,:], ic[1,:] * dist_factor) 
                                            + np.interp(10**self.data[:,0], sync[0,:], sync[1,:] * dist_factor)) )**2 )
        if self.method == 'ps':
            print('LOSS: ', loss, "   ", params[2], '\n')
            return np.array([loss])
        print('LOSS: ', loss, '\n')
        return loss
