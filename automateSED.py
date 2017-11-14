'''Program to automate production of Sychrotron SED from all files'''

import subprocess
import sys

subprocess.call(['gcc', 'jet_synconlyALPWFlux.c', 'jet_fns.c', 'jet_fns.h'])
for i in range(30):
#i=1
    #if (i>10) and (i<21):
    EVPA_rotation = 1
    DD = 1
    #else:
        #EVPA_rotation = 0
    subprocess.call(['./a.out',str(EVPA_rotation),str(i),str(DD)])
    subprocess.call(['python','process_syncALP.py'])
