'''Program to automate production of Sychrotron SED from all files'''

import subprocess
import sys

subprocess.call(['gcc', 'jet_synconlyALPWFlux.c', 'jet_fns.c', 'jet_fns.h'])
for i in range(100):
    print(i)
    EVPA_rotation = 0
    #if (i>=20) and (i<60):
        #EVPA_rotation = 1
    DD = 0
    #else:
        #EVPA_rotation = 0
    subprocess.call(['./a.out',str(EVPA_rotation),str(i),str(DD)])
    subprocess.call(['python','process_syncALP.py'])

#1737
#1837
'''
import numpy as np
keydat = np.loadtxt('keyparams.txt')
theta_obs=keydat[2]
gamma_bulk = keydat[1]
plot = 0
id = 1557
data_rotate = 1
rot_angle = -40*np.pi/180
if data_rotate == 1:
    Pol = np.loadtxt('EVPA+Pol.txt')
    for i in range(80):
        for j in range(4):
            A = Pol[(id-1+i),j]
            if (A + rot_angle) < -np.pi/2:
                Pol[(id-1+i),j] = np.pi/2 + ((Pol[(id-1+i),j] + rot_angle) + np.pi/2)
            elif (A + rot_angle) > np.pi/2:
                Pol[(id-1+i),j] = -np.pi/2 + ((Pol[(id-1+i),j] + rot_angle) - np.pi/2)
            else:
                Pol[(id-1+i),j] += rot_angle
    np.savetxt('EVPA+PolTEST.txt', Pol,'%.4lf','\t')

Rotation_unwrapper = 0
if Rotation_unwrapper ==1:
    Pol = np.loadtxt('EVPA+Pol.txt')
    for j in range(4):
        Pol[(id-1+28):,j] += np.pi
        Pol[(id-1+44):,j] += np.pi
        Pol[(id-1+52):,j] += np.pi
        Pol[(id-1+61):,j] += np.pi
    Pol[(id-1+39):,3] += np.pi
    Pol[(id-1+37):,1] += np.pi
    Pol[(id-1+36):,2] += np.pi
    #Pol[(id-1+61):,3] += np.pi
    #Pol[(id-1+61):,1] += np.pi
    #Pol[(id-1+61):,2] += np.pi
    np.savetxt('EVPA+PolTEST.txt', Pol,'%.4lf','\t')

if plot == 1:
    import math as math
    import numpy as np
    import matplotlib.pyplot as plt
    # used to plot the polarisation and EVPA
    Pol = np.loadtxt('EVPA+PolTEST.txt')
    t = [i for i in range(len(Pol))]
    fig = plt.figure(15)
    #plt.plot(t[id-1:(id-1+80)],Pol[(id-1):(id-1+80),5],'g--')
    #plt.plot(t[(id-1):(id-1+80)],Pol[(id-1):(id-1+80),1],'g-')
    #plt.plot(t[(id-1):(id-1+80)],Pol[(id-1):(id-1+80),6],'b--')
   # plt.plot(t[(id-1):(id-1+80)],Pol[(id-1):(id-1+80),2],'b-')
    plt.plot(t[(id-1):(id-1+80)],Pol[(id-1):(id-1+80),7],'r--')
    plt.plot(t[(id-1):(id-1+80)],Pol[(id-1):(id-1+80),3],'r-')
    plt.xlabel('Time [d]')
    plt.ylabel('EVPA [rad], $\Pi$')
    '''
'''
    plt.text(id+3,1.4,'$\Pi_{meanGam} =$ %.5f' % np.mean(Pol[(id-1):(id-1+100),7]), fontsize=10)
    plt.text(id+3,1.25,'$\sigma_{\Pi} =$ %.5f' % np.std(Pol[(id-1):(id-1+100),7]), fontsize=10)
    plt.text(id+3,1.1,'$\Pi_{meanRad} =$ %.5f' % np.mean(Pol[(id-1):(id-1+100),6]), fontsize=10)
    plt.text(id+3,0.95,'$\sigma_{\Pi} =$ %.5f' % np.std(Pol[(id-1):(id-1+100),6]), fontsize=10)
    plt.text(id+3,0.8,'$\Pi_{meanOpt} =$ %.5f' % np.mean(Pol[(id-1):(id-1+100),5]), fontsize=10)
    plt.text(id+3,0.65,'$\sigma_{\Pi} =$ %.5f' % np.std(Pol[(id-1):(id-1+100),5]), fontsize=10)
    plt.text(id+3,0.5,'$\sigma_{EVPA_{Gam}} =$ %.5f' % np.std(Pol[(id-1):(id-1+100),3]), fontsize=10)
    plt.text(id+3,0.35,'$\sigma_{EVPA_{Rad}} =$ %.5f' % np.std(Pol[(id-1):(id-1+100),2]), fontsize=10)
    plt.text(id+3,0.2,'$\sigma_{EVPA_{Opt}} =$ %.5f' % np.std(Pol[(id-1):(id-1+100),1]), fontsize=10)
    '''
'''
    plt.text(id+3,10,'$\Pi_{meanGam} =$ %.5f' % np.mean(Pol[(id-1):(id-1+20),7]), fontsize=10)
    plt.text(id+3,9.25,'$\sigma_{\Pi} =$ %.5f' % np.std(Pol[(id-1):(id-1+20),7]), fontsize=10)
    plt.text(id+3,8.5,'$\Pi_{meanRad} =$ %.5f' % np.mean(Pol[(id-1):(id-1+20),6]), fontsize=10)
    plt.text(id+3,7.75,'$\sigma_{\Pi} =$ %.5f' % np.std(Pol[(id-1):(id-1+20),6]), fontsize=10)
    plt.text(id+3,7.0,'$\Pi_{meanOpt} =$ %.5f' % np.mean(Pol[(id-1):(id-1+20),5]), fontsize=10)
    plt.text(id+3,6.25,'$\sigma_{\Pi} =$ %.5f' % np.std(Pol[(id-1):(id-1+20),5]), fontsize=10)
    plt.text(id+3,5.5,'$\sigma_{EVPA_{Gam}} =$ %.5f' % np.std(Pol[(id-1):(id-1+20),3]), fontsize=10)
    plt.text(id+3,4.75,'$\sigma_{EVPA_{Rad}} =$ %.5f' % np.std(Pol[(id-1):(id-1+20),2]), fontsize=10)
    plt.text(id+3,4.0,'$\sigma_{EVPA_{Opt}} =$ %.5f' % np.std(Pol[(id-1):(id-1+20),1]), fontsize=10)
    plt.text(id+55,10,'$EVPA_{Gammean} =$ %.5f' % np.mean(Pol[(id-1):(id-1+20),3]), fontsize=10)
    plt.text(id+55,9.25,'$EVPA_{Radmean} =$ %.5f' % np.mean(Pol[(id-1):(id-1+20),2]), fontsize=10)
    plt.text(id+55,8.5,'$EVPA_{Optmean} =$ %.5f' % np.mean(Pol[(id-1):(id-1+20),1]), fontsize=10)
    plt.show()

    fig.savefig('/Users/ALP/Desktop/DD'+str(DD)+'_Rot'+str(EVPA_rotation)+'_theta'+str(theta_obs)+'_gamma'+str(gamma_bulk)+'.png') '''