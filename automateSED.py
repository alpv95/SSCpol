'''Program to automate production of Synchrotron SED from all files'''

import subprocess
import sys

subprocess.call(['gcc', 'jet_synconlyALPWFlux.c','mtwister.c','mtwister.h', 'jet_fns.c', 'jet_fns.h'])
for i in range(1):
    EVPA_rotation = 1
    theta_obs=3.3
    DD = 1
    #if i>=150:
        #theta_obs = 1.5
    #if (i>=40) and (i<110):
        #EVPA_rotation = 1
    #elif (i>=190) and (i<260):
       # EVPA_rotation = 1
    #if i>=300:
        #DD = 0
    #if i>=340:
        #theta_obs = 4.0
    #if i<10:
        ##ting = (i+3)*500-1
    #else:
        #ting = i*(i-4)+4

    #else:
        #EVPA_rotation = 0
    subprocess.call(['./a.out',str(EVPA_rotation),str(i),str(DD),str(theta_obs)])
    subprocess.call(['python','process_syncALP.py'])
    print(i)


import numpy as np
import math as math
import matplotlib.pyplot as plt
keydat = np.loadtxt('keyparams.txt')
#theta_obs=keydat[2]
#gamma_bulk = keydat[1]
DD=1
EVPA_rotation = 1
plot = 0
id = 4197  #2161-2192 => 127blocks,NoRot; 2192-2224 => 127blocks,Rot;2225-2256 => 61blocks,NoRot;2257-2288 => 61blocks,Rot;
#2289-
#2322 norot 127
#2352 rot 127
#2385 rot 127
#2417 rot 127 DD0
#2449 norot 127 DD0
#2481 rot 127 theta1
#2513 rot 127 theta3
#2545 100step 127 theta3 gamma13
#2625 rot 61 theta3 gamma13
#2657 rot 61 theta2 gamma13
#2688 rot 127 theta2 gamma13 64steps
#2753 rot 127 theta1 gamma9.5
#2785 rot 127 theta1 gamma11.5
#2817 rot 127 theta1 gamma13.5
#2848 rot 127 theta1 gamma15.5
#2880 rot 127 theta1 gamma19
#2912 rot 127 theta1 gamma14 100steps
#3012 rot 127 theta3 gamma6
#3044 rot 127 theta3 gamma6 100steps
#3144 rot 127 theta4 gamma5 good one!
#3176 rot 127 theta4 gamma5 only 2 inner circles are helical here
#3208 rot 127 theta4 gamma5 only 2 inner circles helical 100steps in paper! (Using first 20 of this as rndom exposition)
#3308 rot 127 theta1.5 gamma5 only 2 inner circles helical 100steps in paper!
#3408 norot 127 theta4 gamma5 noDD for paper
#3428 ...... 40 norot, 70rot, 40 norot theta4, same again for theta1.5, 40 norot DD0 theta1.5, 40 norot DD0 theta4
#3808 ...... only m 1st inner circle is helical, same as above but cut short at 4118
#4118 norot 127 theta1.5 gamma5 sameseed
#4158 norot 127 theta1.5 gamma5 DD0 samessed as above
#4198 norot 127 theta1.5 gamma5 sameseed
#4230 norot 127 theta1.5 gamma5 DD0
#4262
data_rotate = 0
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
    #Pol[(id-1+19):,3] += np.pi
    Pol[(id-1+29):,3] += np.pi
    Pol[(id-1+45):,3] += np.pi
    Pol[(id-1+61):,3] += np.pi
    Pol[(id-1+77):,3] += np.pi
    Pol[(id-1+7):,2] += np.pi
    Pol[(id-1+23):,2] += np.pi
    Pol[(id-1+7):,1] += np.pi
    Pol[(id-1+23):,1] += np.pi
    #Pol[(id-1+57):,j] += np.pi
    #Pol[(id-1+61):,j] += np.pi
    #Pol[(id-1+39):,3] += np.pi
    #Pol[(id-1+37):,1] += np.pi
    #Pol[(id-1+36):,2] += np.pi
    #Pol[(id-1+61):,3] += np.pi
    #Pol[(id-1+61):,1] += np.pi
    #Pol[(id-1+61):,2] += np.pi
    np.savetxt('EVPA+PolTEST.txt', Pol,'%.4lf','\t')

nplot = 32
if plot == 1:
    # used to plot the polarisation and EVPA
    Pol = np.loadtxt('EVPA+Pol.txt')
    t = [i for i in range(len(Pol))]
    #fig = plt.figure(15)
    #plt.plot(t[id-1:(id-1+nplot)],Pol[(id-1):(id-1+nplot),5],'g--')
    #plt.plot(t[(id-1):(id-1+nplot)],Pol[(id-1):(id-1+nplot),1],'g-')
    #plt.plot(t[(id-1):(id-1+nplot)],Pol[(id-1):(id-1+nplot),6],'r--')
    #plt.plot(t[(id-1):(id-1+nplot)],Pol[(id-1):(id-1+nplot),2],'r-')
    f, (ax1,ax2) = plt.subplots(2,sharex=True,sharey=True)
    ax1.plot(t[(id-1):(id-1+nplot)],Pol[(id-1):(id-1+nplot),2]*180/np.pi,'r-',label='Radio EVPA') #10
    ax1.plot(t[(id-1):(id-1+nplot)],Pol[(id-1):(id-1+nplot),1]*180/np.pi,'g-', label='Optical EVPA') #11
    ax1.plot(t[(id-1):(id-1+nplot)],Pol[(id-1):(id-1+nplot),3]*180/np.pi,'b-',label='X-ray EVPA') #11
    ax1.set_ylabel('EVPA [deg]')
    ax1.set_yticks([-90,90,-45,45,0])

    ax1.tick_params(axis='y')
    ax15 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    ax15.set_ylabel('$\Pi$')  # we already handled the x-label with ax1

    ax15.plot(t[(id-1):(id-1+nplot)],Pol[(id-1):(id-1+nplot),7],'--',color='b')
    ax15.plot(t[(id-1):(id-1+nplot)],Pol[(id-1):(id-1+nplot),6],'r--')
    ax15.plot(t[id-1:(id-1+nplot)],Pol[(id-1):(id-1+nplot),5],'g--')
    ax15.set_yticks([1,0.75,0.5,0.25,0])

    ax2.plot(t[(id-1):(id-1+nplot)],Pol[(id-1+nplot):(id-1+2*nplot),2]*180/np.pi,'r-',label='Radio EVPA') #10
    ax2.plot(t[(id-1):(id-1+nplot)],Pol[(id-1+nplot):(id-1+2*nplot),1]*180/np.pi,'g-', label='Optical EVPA') #11
    ax2.plot(t[(id-1):(id-1+nplot)],Pol[(id-1+nplot):(id-1+2*nplot),3]*180/np.pi,'b-',label='X-ray EVPA') #11
    ax2.tick_params(axis='y')
    ax2.set_xlabel('Time [d]')
    ax25 = ax2.twinx()  # instantiate a second axes that shares the same x-axis
    #ax15.get_shared_y_axes().join(ax15, ax25)
    ax25.plot(t[(id-1):(id-1+nplot)],Pol[(id-1+nplot):(id-1+2*nplot),7],'--',color='b')
    ax25.plot(t[(id-1):(id-1+nplot)],Pol[(id-1+nplot):(id-1+2*nplot),6],'r--')
    ax25.plot(t[id-1:(id-1+nplot)],Pol[(id-1+nplot):(id-1+2*nplot),5],'g--')
    ax25.set_yticks([1,0.75,0.5,0.25,0])

    # Fine-tune figure; make subplots close to each other and hide x ticks for
    # all but bottom plot.
    f.subplots_adjust(hspace=0)
    plt.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    plt.show()
    '''
    fig, ax1 = plt.subplots()

    color = 'tab:red'
    ax1.set_xlabel('Time [d]')
    ax1.set_ylabel('EVPA [deg]')
    PAR, = ax1.plot(t[(id-1):(id-1+nplot)],Pol[(id-1):(id-1+nplot),2]*180/np.pi,'r-',label='Radio EVPA') #10
    PAO, = ax1.plot(t[(id-1):(id-1+nplot)],Pol[(id-1):(id-1+nplot),1]*180/np.pi,'g-', label='Optical EVPA') #11
    PA, = ax1.plot(t[(id-1):(id-1+nplot)],Pol[(id-1):(id-1+nplot),3]*180/np.pi,'b-',label='X-ray EVPA') #11
    '''
    '''
    ax1.plot(t[(id-1+10):(id-1+26)],Pol[(id-1+10):(id-1+26),2]*180/np.pi,'r-',label='Radio EVPA')
    ax1.plot(t[(id-1+11):(id-1+27)],Pol[(id-1+11):(id-1+27),1]*180/np.pi,'g-', label='Optical EVPA')
    ax1.plot(t[(id-1+11):(id-1+27)],Pol[(id-1+11):(id-1+27),3]*180/np.pi,'b-',label='X-ray EVPA')
    ax1.plot(t[(id-1+26):(id-1+42)],Pol[(id-1+26):(id-1+42),2]*180/np.pi,'r-',label='Radio EVPA')
    ax1.plot(t[(id-1+27):(id-1+43)],Pol[(id-1+27):(id-1+43),1]*180/np.pi,'g-', label='Optical EVPA')
    ax1.plot(t[(id-1+27):(id-1+43)],Pol[(id-1+27):(id-1+43),3]*180/np.pi,'b-',label='X-ray EVPA')
    ax1.plot(t[(id-1+42):(id-1+58)],Pol[(id-1+42):(id-1+58),2]*180/np.pi,'r-',label='Radio EVPA')
    ax1.plot(t[(id-1+43):(id-1+59)],Pol[(id-1+43):(id-1+59),1]*180/np.pi,'g-', label='Optical EVPA')
    ax1.plot(t[(id-1+43):(id-1+59)],Pol[(id-1+43):(id-1+59),3]*180/np.pi,'b-',label='X-ray EVPA')
    ax1.plot(t[(id-1+58):(id-1+70)],Pol[(id-1+58):(id-1+70),2]*180/np.pi,'r-',label='Radio EVPA')
    ax1.plot(t[(id-1+59):(id-1+70)],Pol[(id-1+59):(id-1+70),1]*180/np.pi,'g-', label='Optical EVPA')
    ax1.plot(t[(id-1+59):(id-1+70)],Pol[(id-1+59):(id-1+70),3]*180/np.pi,'b-',label='X-ray EVPA')
    '''
    '''
    #ax1.plot([3616,3686],[-25,767],color='k')
    ax1.tick_params(axis='y')

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    color = 'tab:blue'
    ax2.set_ylabel('$\Pi$')  # we already handled the x-label with ax1
    Pi, = ax2.plot(t[(id-1):(id-1+nplot)],Pol[(id-1):(id-1+nplot),7],'--',color='b')
    PiR, = ax2.plot(t[(id-1):(id-1+nplot)],Pol[(id-1):(id-1+nplot),6],'r--')
    PiO, = ax2.plot(t[id-1:(id-1+nplot)],Pol[(id-1):(id-1+nplot),5],'g--')
    ax2.set_yticks([-1,-0.5,0,0.5,1])
    #fig.legend([PA, Pi], ['EVPA', 'Polarization Fraction'], loc='upper left')
    #fig.tight_layout()  # otherwise the right y-label is slightly clipped
    #ax1.legend([PA,PAR,PAO,Pi,PiR,PiO],['X-ray EVPA','Radio EVPA','Optical EVPA',r'X-ray $\Pi$',r'Radio $\Pi$',r'Optical $\Pi$'])
    #ax1.legend([PA,Pi],['EVPA',r'$\Pi$'])
    #plt.show()
    '''
    '''
    plt.text(id+53,6,'$\Pi_{meanGam} =$ %.5f' % np.mean(Pol[(id-1+20):(id-1+60),7]), fontsize=10)
    plt.text(id+53,5.25,'$\sigma_{\Pi} =$ %.5f' % np.std(Pol[(id-1+20):(id-1+60),7]), fontsize=10)
    plt.text(id+53,4.5,'$\Pi_{meanRad} =$ %.5f' % np.mean(Pol[(id-1+20):(id-1+60),6]), fontsize=10)
    plt.text(id+53,3.75,'$\sigma_{\Pi} =$ %.5f' % np.std(Pol[(id-1+20):(id-1+60),6]), fontsize=10)
    plt.text(id+53,3.0,'$\Pi_{meanOpt} =$ %.5f' % np.mean(Pol[(id-1+20):(id-1+60),5]), fontsize=10)
    plt.text(id+53,2.25,'$\sigma_{\Pi} =$ %.5f' % np.std(Pol[(id-1+20):(id-1+60),5]), fontsize=10)
    plt.text(id+53,1.5,'$\sigma_{EVPA_{Gam}} =$ %.5f' % np.std(Pol[(id-1+20):(id-1+60),3]), fontsize=10)
    plt.text(id+53,0.75,'$\sigma_{EVPA_{Rad}} =$ %.5f' % np.std(Pol[(id-1+20):(id-1+60),2]), fontsize=10)
    plt.text(id+53,0.2,'$\sigma_{EVPA_{Opt}} =$ %.5f' % np.std(Pol[(id-1+20):(id-1+60),1]), fontsize=10)
    '''
    '''
    ax1.text(id+3,-85,'$\Pi_{meanGam} =$ %.5f' % np.mean(Pol[(id-1):(id-1+nplot),7]), fontsize=10)
    ax1.text(id+3,-75,'$\sigma_{\Pi} =$ %.5f' % np.std(Pol[(id-1):(id-1+nplot),7]), fontsize=10)
    ax1.text(id+3,-65,'$\Pi_{meanRad} =$ %.5f' % np.mean(Pol[(id-1):(id-1+nplot),6]), fontsize=10)
    ax1.text(id+3,-55,'$\sigma_{\Pi} =$ %.5f' % np.std(Pol[(id-1):(id-1+nplot),6]), fontsize=10)
    ax1.text(id+3,-45,'$\Pi_{meanOpt} =$ %.5f' % np.mean(Pol[(id-1):(id-1+nplot),5]), fontsize=10)
    ax1.text(id+3,-35,'$\sigma_{\Pi} =$ %.5f' % np.std(Pol[(id-1):(id-1+nplot),5]), fontsize=10)
    ax1.text(id+3,-25,'$\sigma_{EVPA_{Gam}} =$ %.5f' % np.std(Pol[(id-1):(id-1+nplot),3]*180/np.pi), fontsize=10)
    ax1.text(id+3,-15,'$\sigma_{EVPA_{Rad}} =$ %.5f' % np.std(Pol[(id-1):(id-1+nplot),2]*180/np.pi), fontsize=10)
    ax1.text(id+3,-5,'$\sigma_{EVPA_{Opt}} =$ %.5f' % np.std(Pol[(id-1):(id-1+nplot),1]*180/np.pi), fontsize=10)
    ax1.text(id+3,5,'$EVPA_{Gammean} =$ %.5f' % np.mean(Pol[(id-1):(id-1+nplot),3]*180/np.pi), fontsize=10)
    ax1.text(id+3,15,'$EVPA_{Radmean} =$ %.5f' % np.mean(Pol[(id-1):(id-1+nplot),2]*180/np.pi), fontsize=10)
    ax1.text(id+3,25,'$EVPA_{Optmean} =$ %.5f' % np.mean(Pol[(id-1):(id-1+nplot),1]*180/np.pi), fontsize=10)
    plt.show()
    '''
'''
    fig.savefig('/Users/ALP/Desktop/DD'+str(DD)+'_Rot'+str(EVPA_rotation)+'_theta'+str(theta_obs)+'_gamma'+str(gamma_bulk)+'_id100'+str(id)+'.png')
'''