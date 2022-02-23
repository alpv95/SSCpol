import numpy as np
from sscpol.jet_fns import *
import multiprocessing as mp
from scipydirect import minimize as minimize_direct
from scipy.optimize import minimize
import pyswarms as ps
import multiprocessing as mp
import pickle
from os.path import exists


class SSC_Fitter(object):
    def __init__(self, blazar, method="standard", nblocks=1, seed=100, nprocs=20, rand_gamma=0):
        if blazar not in ["J2011", "TXS", "S5low", "S5flare", "J2011flare"]:
            raise ValueError("blazar must be one of J2011, TXS, S5low, S5flare")
        if blazar == "J2011":
            self.d_Blazar = 1000E6*3.08E18 
            self.z = 0.2
            self.pol = 0.231 #+- 6.93%
            # self.x0 = [37.4,9.8,1.894,40.0,8.7,-4.6,1.5,1,5.73]
            self.x0 = [37.50040529 ,10.01467913 , 1.89450649, 41.99447311,  8.31244801 ,-4.58830858, 2.16303875 , 0.97451764 , 5.91011898]
            #[4.46684e+37, 1.42666e+10, 1.96358e+00, 1.80864e+01, 1.14136e+01, 3.02740e-05, 2.99335e+00, 1.13302e+00]
            data = np.loadtxt('data/new_data_sed_CGRaBSJ0211+1051_XMM.txt')
            data2 = np.loadtxt('data/new_data_sed_CGRaBSJ0211+1051.txt')
            self.data = np.concatenate([data[6:],data2[12:]], axis=0)
        elif blazar == "J2011flare":
            self.d_Blazar = 1000E6*3.08E18 
            self.z = 0.2
            self.pol = 0.231 #+- 6.93%
            # self.x0 = [37.4,10.3,2.03,40.0,18.7,-4.2,1.53,1,5.73]
            self.x0 = [37.68246069 ,10.15999064 , 2.04750948, 42.89329777, 18.78745097, -4.19312741, 2.67427829 , 0.87634521 , 5.79508428]
            data = np.loadtxt('data/CGRaBSJ02111051-builderX.txt')
            data2 = np.loadtxt('data/new_data_sed_CGRaBSJ0211+1051.txt')
            self.data = np.concatenate([data,data2[12:]], axis=0)
        elif blazar == "S5low":
            self.d_Blazar = 1627E6*3.08E18 
            self.z = 0.31
            self.pol =  0.109 #+-4.89%
            # self.x0 = [np.log10(4.52530e+38), np.log10(1.92333e+10), 1.79025, 58.1863,	15.1244, np.log10(3.36434e-05), 3.38444, 1.13037, 5.73]
            self.x0 = [38.58034712 ,10.10411188 , 1.72720727 ,40.17342854 ,15.32114366 , -4.46809398, 3.98353423 , 1.08008377 , 5.86056596]
            data = np.loadtxt('data/S5low.txt')
            self.data = np.insert(data, 1, 0, axis=1)
        elif blazar == "S5flare":
            self.d_Blazar = 1627E6*3.08E18 
            self.z = 0.31
            # self.x0 = [np.log10(4.52530e+38), np.log10(1.92333e+10), 1.79025, 40.1863,	15.1244, np.log10(3.36434e-05), 3.38444, 1.13037, 5.73]
            self.x0 = [38.68526327, 10.30926079,  1.67885216, 40.9791205,  15.02985019, -4.54369785, 3.31040824,  1.13898073,  5.82027554,]
            self.pol =  0.109 #+-4.89%
            self.data = np.loadtxt('data/S5flare.txt')
        else:
            self.d_Blazar = 1789E6*3.08E18
            self.pol = 0.09
            #1.8 +/- 0.8%, 6.2 +/- 1.2%, 6.6 +/- 2.3%, 14%, 7.9+/-1.1%.
            self.x0 = [37.81, 10.2, 1.97865, 32.7639, 17.3750, -4.1, 1.58989, 1.20151, 5.73]
            # self.x0 = [37.4426634 ,  9.70125571 , 1.67955008, 34.60545266, 16.70775311 , -3.1, 2.85939946 , 1.08650853 , 5.93446008,]
            self.z = 0.3365
            data = np.loadtxt('data/txs.txt')
            self.data = np.insert(data, 1, 0, axis=1)
        
        if method not in ["standard", "direct", "ps", 'cross_entropy']:
            raise ValueError("method must be one of ip, direct, ps")
        self.method = method
        self.nblocks = nblocks
        self.seed_n = seed
        self.blazar = blazar
        self.n_procs = nprocs
        self.rand_gamma = rand_gamma
        self.bounds = [(36.3, 39),(9.7,13),(1.65,2.05),(10,45),(7,30),(-4.7,-3.1),(0.1,6),(0.8,1.4),(5.72,6.8)]
        self.filename = 'data/checkpoints/'+str(self.blazar)+'_'+str(self.nblocks)+'_g'+str(self.rand_gamma)+'.p'

    def run(self,):
        print("bounds: ", self.bounds)
        
        if self.method == 'standard':
            res = minimize(self.loglike, self.x0, method='Nelder-Mead', bounds=self.bounds, options={'disp': True})
        elif self.method == 'cross_entropy':
            res = self.cross_entropy_method()
        elif self.method == 'ps':
            options = {'c1': 0.5, 'c2': 0.3, 'w':0.9}
            # Call instance of GlobalBestPSO
            optimizer = ps.single.GlobalBestPSO(n_particles=15, dimensions=len(self.bounds), 
                                                bounds=(np.array([b[0] for b in self.bounds]),np.array([b[1] for b in self.bounds])),
                                                options=options)
            # Perform optimization
            res = optimizer.optimize(self.loglike, iters=10000, n_processes=15)
        elif self.method == "direct":
            res = minimize_direct(self.loglike, self.bounds, fglobal=1.0, fglper=(0.2*100), disp=True, maxt=10000)

        print("res: ", res)
        return res

    def loglike(self, params):
        params = np.squeeze(np.array(params))
        params[0] = 10**params[0]
        params[1] = 10**params[1]
        params[5] = 10**params[5]
        params[8] = 10**params[8]

        seed = np.random.RandomState().randint(self.seed_n)
        sync, ic = run_ssc(params, nblocks=self.nblocks, seed=seed, rand_gamma=self.rand_gamma)

        dist_factor = 1.0E7*(1.0/((4.0*np.pi*self.d_Blazar**2.0)*(1.0+self.z)**2.0)) 
        loss = np.sum((self.data[:,2] - np.log10(np.interp(10**self.data[:,0], ic[0,:], ic[1,:] * dist_factor) 
                                            + np.interp(10**self.data[:,0], sync[0,:], sync[1,:] * dist_factor)) )**2 )
        if self.method == 'ps':
            print('LOSS: ', loss, "   ", params[2], '\n')
            return np.array([loss])
        print('LOSS: ', loss, '\n')
        return loss, sync, ic

    def get_checkpoint(self):
        if exists(self.filename):
            print("Checkpoint found, loading...")
            return pickle.load(open(self.filename, 'rb'))
        else:
            return self.x0, np.diag((np.array(self.bounds)[:,1] - np.array(self.bounds)[:,0]) / 50)
    
    def save_checkpoint(self, mu, Sigma):
        pickle.dump((mu, Sigma), open(self.filename, 'wb'))

    def cross_entropy_method(self, n_samples=60, n_elite=15, max_k=12,):
        bounds = np.array(self.bounds)
        current_mean_likelihood = 1000
        mu, Sigma = self.get_checkpoint()

        sample_buffer = []
        likelihood_buffer = []
        best_estimate = mu
        best_likelihood = current_mean_likelihood
        for k in range(max_k):
            #Sample from mutlivariate gaussian with mean mu and covariance Sigma
            samples = np.random.multivariate_normal(mu, Sigma, size=n_samples)
            #Check if samples are within bounds
            samples = np.array([np.clip(s, bounds[:,0], bounds[:,1]) for s in samples])

            #Evaluate the likelihood of each sample in parallel with multiprocessing using the loglikelihood function
            with mp.Pool(processes=self.n_procs) as pool:
                results = pool.map(self.loglike, samples)

            log_likelihoods, syncs, ics = zip(*results)
            mask = np.isfinite(log_likelihoods)
            log_likelihoods = [l for i,l in enumerate(log_likelihoods) if mask[i]]
            samples = samples[mask]

            likelihood_buffer = log_likelihoods + likelihood_buffer
            sample_buffer = [samples] + sample_buffer

            samples = np.concatenate(sample_buffer, axis=0)
            log_likelihoods = np.array(likelihood_buffer)

            #If number of nan values is less than 2*n_elite, then skip back to the previous iteration
            if len(likelihood_buffer) < 4*n_elite:
                print('skipped', k)
                continue
            sample_buffer = []
            likelihood_buffer = []

            #Sort the samples according to the log likelihoods
            sorted_samples = np.array(samples)[np.argsort(log_likelihoods)]
            sorted_likelihoods = log_likelihoods[np.argsort(log_likelihoods)]

            #Select the best n_elite samples
            elite_samples = sorted_samples[:n_elite]
            elite_likelihoods = sorted_likelihoods[:n_elite]

            #Update the mean and covariance of the gaussian
            mu = np.mean(elite_samples, axis=0)
            Sigma = np.cov(elite_samples.T)
            current_mean_likelihood = np.mean(elite_likelihoods)
            if current_mean_likelihood < best_likelihood:
                best_estimate = mu
                best_likelihood = current_mean_likelihood
                print('best', best_likelihood, best_estimate)

            print("current_mean_likelihood: ", current_mean_likelihood)
            print("current mean and Sigma: ", repr(mu), repr(Sigma))

            self.save_checkpoint(mu, Sigma)

            if k == max_k - 1: #last iteration
                #save syncs and ics as a joint numpy array
                syncs = np.array(syncs)
                ics = np.array(ics)
                np.save('data/fits/syncs_'+str(self.blazar)+'_'+str(self.nblocks)+'_g'+str(self.rand_gamma)+'.npy', syncs)
                np.save('data/fits/ics_'+str(self.blazar)+'_'+str(self.nblocks)+'_g'+str(self.rand_gamma)+'.npy', ics)

        return current_mean_likelihood, mu, Sigma



