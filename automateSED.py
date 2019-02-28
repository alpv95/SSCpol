'''Program to automate production of Synchrotron SED from all files'''
'''this is one task, ntasks>1 runs this whole file multiple times on different cores - need .mpi file to handle multiple tasks properly'''
'''what we want is ntasks = 1, cpu per task > 1, ie a multithreaded program, can use import threading instead of import multiprocessing if necessary'''
'''use task_id + thread id to seed random numbers'''
import multiprocessing
import subprocess
import sys
import time
import sys
import os
task_id = sys.argv[1]

results_dir = "results"
try:
    # Create results Directory
    os.mkdir(results_dir)
    print("Directory " , results_dir,  " Created ") 
except FileExistsError:
    print("Directory ", results_dir,  " already exists")

subprocess.call(['gcc', 'jet_synconlyALPWFlux.c','mtwister.c','mtwister.h', 'jet_fns.c', 'jet_fns.h','nrutil.c','nrutil.h','-lm'])

q = multiprocessing.Queue()
inputs = [(0,i,1,1.5,1,0,task_id) for i in range(1)] #these are saved in keyparams along with more

for inpt in inputs:
  q.put(inpt)

def worker():
  while True:
    inpt = q.get()
    if inpt is None:  # EOF?
      return
    print("task_id: ",task_id," thread_id: " ,inpt[1])
    #time.sleep(0.1)
    p = subprocess.check_call(['./a.out',str(inpt[0]),str(inpt[1]),str(inpt[2]),str(inpt[3]),str(inpt[4]),str(inpt[5]),str(inpt[6])])
    #checksum = collect_md5_result_for(fileName)
    #result[fileName] = checksum  # store it

threads = [multiprocessing.Process(target=worker) for _i in range(1) ]
for thread in threads:
  #time.sleep(0.3) #make sure C program gets different random seed
  thread.start()
  q.put(None)  # one EOF marker for each thread

for thread in threads: #wait for all workers to finish
  thread.join()

#then run this once all workers have finished
subprocess.check_call('cat results/TESTFIL' + str(task_id) + '_* > results/TESTFIL' + str(task_id) + '.txt', shell=True) #combine all the output files in order of inputs
subprocess.check_call('cat results/basicdata' + str(task_id) + '_* > results/basicdata'+ str(task_id) +'.txt',shell=True)
subprocess.check_call('cat results/freqrange' + str(task_id) + '_* > results/freqrange'+ str(task_id) +'.txt',shell=True)
subprocess.check_call('cat results/keyparams' + str(task_id) + '_* > results/keyparams'+ str(task_id) +'.txt',shell=True)

subprocess.check_call('rm results/TESTFIL' + str(task_id) + '_*',shell=True) #delete files after concat
subprocess.check_call('rm results/basicdata' + str(task_id) + '_*',shell=True)
subprocess.check_call('rm results/freqrange' + str(task_id) + '_*',shell=True)
subprocess.check_call('rm results/keyparams' + str(task_id) + '_*',shell=True)



