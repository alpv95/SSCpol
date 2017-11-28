'''Program to automate production of Sychrotron SED from all files'''

import subprocess
import sys

subprocess.call(['gcc', 'jet_synconlyALPWFlux.c', 'jet_fns.c', 'jet_fns.h'])
for i in range(1):
    #if (i>10) and (i<21):
    EVPA_rotation = 0
    DD = 1
    #else:
        #EVPA_rotation = 0
    subprocess.call(['./a.out',str(EVPA_rotation),str(i),str(DD)])
    subprocess.call(['python','process_syncALP.py'])

plot = 0
id = 489
if plot == 1:
    import math as math
    import numpy as np
    import matplotlib.pyplot as plt
    # used to plot the polarisation and EVPA
    Pol = np.loadtxt('EVPA+Pol.txt')
    t = [i for i in range(len(Pol))]
    plt.figure(12)
    plt.plot(t[id-1:(id-1+17)],Pol[(id-1):(id-1+17),5],'g--')
    plt.plot(t[(id-1):(id-1+17)],Pol[(id-1):(id-1+17),1],'g-')
    plt.plot(t[(id-1):(id-1+17)],Pol[(id-1):(id-1+17),6],'b--')
    plt.plot(t[(id-1):(id-1+17)],Pol[(id-1):(id-1+17),2],'b-')
    plt.plot(t[(id-1):(id-1+17)],Pol[(id-1):(id-1+17),7],'r--')
    plt.plot(t[(id-1):(id-1+17)],Pol[(id-1):(id-1+17),3],'r-')
    plt.xlabel('Time [d]')
    plt.ylabel('EVPA [rad], $\Pi$')
    plt.text(id+3,1.4,'$\Pi_{meanGam} =$ %.5f' % np.mean(Pol[(id-1):(id-1+17),7]), fontsize=10)
    plt.text(id+3,1.25,'$\sigma_{\Pi} =$ %.5f' % np.std(Pol[(id-1):(id-1+17),7]), fontsize=10)
    plt.text(id+3,1.1,'$\Pi_{meanRad} =$ %.5f' % np.mean(Pol[(id-1):(id-1+17),6]), fontsize=10)
    plt.text(id+3,0.95,'$\sigma_{\Pi} =$ %.5f' % np.std(Pol[(id-1):(id-1+17),6]), fontsize=10)
    plt.text(id+3,0.8,'$\Pi_{meanOpt} =$ %.5f' % np.mean(Pol[(id-1):(id-1+17),5]), fontsize=10)
    plt.text(id+3,0.65,'$\sigma_{\Pi} =$ %.5f' % np.std(Pol[(id-1):(id-1+17),5]), fontsize=10)
    plt.show()