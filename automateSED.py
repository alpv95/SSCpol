'''Program to automate production of Synchrotron SED from all files'''
import multiprocessing
import subprocess
import sys
#import junk

subprocess.call(['gcc', 'jet_synconlyALPWFlux.c','mtwister.c','mtwister.h', 'jet_fns.c', 'jet_fns.h','-lm'])

q = multiprocessing.Queue()
inputs = [(0,i,1,0.5,7,1) for i in range(3)]
for inpt in inputs:
  q.put(inpt)

def worker():
  while True:
    inpt = q.get()
    if inpt is None:  # EOF?
      return
    print(inpt[1])
    p = subprocess.check_call(['./a.out',str(inpt[0]),str(inpt[1]),str(inpt[2]),str(inpt[3]),str(inpt[4]),str(inpt[5])])
    #checksum = collect_md5_result_for(fileName)
    #result[fileName] = checksum  # store it

threads = [multiprocessing.Process(target=worker) for _i in range(1) ]
for thread in threads:
  thread.start()
  q.put(None)  # one EOF marker for each thread

