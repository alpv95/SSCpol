'''Program to automate production of Sychrotron SED from all files'''

import subprocess

subprocess.call(['gcc', 'jet_synconlyALPWFlux.c', 'jet_fns.c', 'jet_fns.h'])
'''for i in range(100):'''
subprocess.call(['./a.out'])
subprocess.call(['python','process_syncALP.py'])
