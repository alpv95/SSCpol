import time
import sys
import os
import utils.junk as junk
import subprocess
import argparse

#Command Line argument Parsing: ###########################
parser = argparse.ArgumentParser()
parser.add_argument('--nblocks', type=int, default=1, choices=[1,7,19,37,64,91,127],
                    help='number of blocks in jet model, default=1')
parser.add_argument('--nsteps', type=int, default=100,
                    help='number of jet calculation steps, default=100')
parser.add_argument('--nworkers', type=int, default=1,
                    help='number of multithread workers, default=1')
parser.add_argument('theta_obs', type=float,
                    help='theta_obs in the lab frame')
parser.add_argument('-ssc','--SSC', action='store_true',
                    help='include SSC calculation (long)')
parser.add_argument('-d','--dir', default="workdir",
                    help='results directory, default=workdir')
args = parser.parse_args()

n_rings =[0,1,2,3,4,5,6][[1,7,19,37,64,91,127].index(args.nblocks)] #rings fixed by number of blocks
task_id = 0 #for multiprocessing
############################################################

subprocess.call(['gcc-9','-fopenmp','src/jet_model.c','src/mtwister.c','include/mtwister.h', 'src/jet_fns.c', 'include/jet_fns.h','-lm','-I','include/']) #Compile C
inpts = [(0,i,1,args.theta_obs, args.nblocks, n_rings, task_id, args.nsteps, args.dir, int(args.SSC)) for i in range(args.nworkers)]

try:
    # Create results Directory
    os.mkdir(args.dir)
    print("Directory " , args.dir,  " Created ")
except FileExistsError:
    print("Directory ", args.dir,  " already exists")

for i in range(args.nworkers):
    inpt = inpts[i]
    subprocess.check_call(['./a.out', str(inpt[0]),str(inpt[1]),str(inpt[2]),str(inpt[3]),str(inpt[4]),str(inpt[5]),str(inpt[6]),str(inpt[7]),str(inpt[8]),str(inpt[9]) ])

#plot results
junk.plot_SED(args.dir + "/pi0_0.txt",args.dir + "/keyparams0_0.txt",args.dir + "/freqrange0_0.txt")

#then run this once all workers have finished, useful if nworkers > 1
subprocess.check_call('cat ' + args.dir + '/pi' + str(task_id) + '_* > ' + args.dir + '/pi' + str(args.theta_obs) + '.txt', shell=True) #combine all the output fil    es in order of inputs
subprocess.check_call('cat ' + args.dir + '/basicdata' + str(task_id) + '_* > ' + args.dir + '/basicdata'+ str(args.theta_obs) +'.txt',shell=True)
subprocess.check_call('cat ' + args.dir + '/freqrange' + str(task_id) + '_* > ' + args.dir + '/freqrange'+ str(args.theta_obs) +'.txt',shell=True)
subprocess.check_call('cat ' + args.dir + '/keyparams' + str(task_id) + '_* > ' + args.dir + '/keyparams'+ str(args.theta_obs) +'.txt',shell=True) 
subprocess.check_call('cat ' + args.dir + '/IC_Z' + str(task_id) + '_* > ' + args.dir + '/IC_Z'+ str(args.theta_obs) +'.txt',shell=True) 
subprocess.check_call('cat ' + args.dir + '/S_Z' + str(task_id) + '_* > ' + args.dir + '/S_Z'+ str(args.theta_obs) +'.txt',shell=True) 

subprocess.check_call('rm ' + args.dir + '/pi' + str(task_id) + '_*',shell=True) #delete files after concat
subprocess.check_call('rm ' + args.dir + '/basicdata' + str(task_id) + '_*',shell=True)
subprocess.check_call('rm ' + args.dir + '/freqrange' + str(task_id) + '_*',shell=True)
subprocess.check_call('rm ' + args.dir + '/keyparams' + str(task_id) + '_*',shell=True)
subprocess.check_call('rm ' + args.dir + '/IC_Z' + str(task_id) + '_*',shell=True)
subprocess.check_call('rm ' + args.dir + '/S_Z' + str(task_id) + '_*',shell=True)

