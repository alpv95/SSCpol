import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import math as math
from jet_fns import *
import copy as copy
import csv
reader = csv.reader(open("MRK501DATA.csv", "r"), delimiter=",")
x = list(reader)

Mkn501_pts = np.loadtxt('MKN501.txt')
ph_energy_MK501 = Mkn501_pts[:,0]
flux_MK501 = Mkn501_pts[:,1]
mkn501IC = np.array(x).astype("float")
ph_energy_MK501 = np.append(ph_energy_MK501,mkn501IC[:,1])
flux_MK501 = np.append(flux_MK501,mkn501IC[:,0])
deleters = []
for i,ph in enumerate(ph_energy_MK501):
    if ph >= 2e9 or (ph>1e5 and ph<2e5):
        deleters.append(i)
ph_energy_MK501 = np.delete(ph_energy_MK501,deleters)
flux_MK501 = np.delete(flux_MK501,deleters)


Pol = np.loadtxt('POL4.txt')
EVPA = np.loadtxt('EVPA4.txt')
P_detected = np.loadtxt('PSYNC.txt')
P_detected_IC = np.loadtxt('PIC.txt')
fq_mids = np.loadtxt('FREQSYNC.txt')
fq_mids_IC = np.loadtxt('FREQIC.txt')

P_detected3 = np.loadtxt('PSYNC3.txt')
P_detected_IC3 = np.loadtxt('PIC3.txt')
fq_mids3 = np.loadtxt('FREQSYNC3.txt')
fq_mids_IC3 = np.loadtxt('FREQIC3.txt')

Pol2 = np.loadtxt('POL3.txt')
EVPA2 = np.loadtxt('EVPA3.txt')
#fq_mids2 = np.loadtxt('FREQSYNC2.txt')



EVPA = np.mean(np.array(EVPA),axis=0)
EVPA2 = np.mean(np.array(EVPA2),axis=0)

data_rotate = 1
rot_angle = -40*np.pi/180
if data_rotate == 1:
    for i in range(len(EVPA)):
        A = EVPA[i]
        if (A + rot_angle) < -np.pi/2:
            EVPA[i] = np.pi/2 + ((EVPA[i] + rot_angle) + np.pi/2)
        elif (A + rot_angle) > np.pi/2:
            EVPA[i] = -np.pi/2 + ((EVPA[i] + rot_angle) - np.pi/2)
        else:
            EVPA[i] += rot_angle
    for i in range(len(EVPA2)):
        A = EVPA2[i]
        if (A + rot_angle) < -np.pi/2:
            EVPA2[i] = np.pi/2 + ((EVPA2[i] + rot_angle) + np.pi/2)
        elif (A + rot_angle) > np.pi/2:
            EVPA2[i] = -np.pi/2 + ((EVPA2[i] + rot_angle) - np.pi/2)
        else:
            EVPA2[i] += rot_angle

fig = plt.figure(2)
# set height ratios for sublots
gs = gridspec.GridSpec(3, 1, height_ratios=[1,1,2])
# the fisrt subplot
ax0 = plt.subplot(gs[0])
# log scale for axis X of the first subplot
#ax0.set_title('$W_j=3*10^{37}W$, $L_{jet} = 5*10^{20}m$, $B_0=8*10^{-5}T$, $E_{max}=50*10^9eV$,\n' + r'$\alpha=1.95$, $\theta_{opening} = 9^{\circ}$, $\theta_{obs}=4.0^{\circ}$, $\gamma_{bulk}=7.5$')
ax0.set_xscale("log")
ax0.set_xlim([1E-6, 1E13])
ax0.set_ylim([0, 0.2])
ax0.set_ylabel(r'$\Pi(\omega)$', size='13')
line0 = ax0.plot(freqtoeV(fq_mids[1:39]), np.mean(Pol[:,1:39],axis=0),'m',label='Non Rotating')
line011 = ax0.plot(freqtoeV(fq_mids3[1:39]), np.mean(Pol2[:,1:39],axis=0),'m--',label='Rotating')
#line1 = ax0.plot(freqtoeV(fq_mids), Pol_Init,'r',label='Initial Population')
#line2 = ax0.plot(freqtoeV(fq_mids), Pol_single, color='g',label='Single e')
#ax0.text(1E7,0.3,'$\Pi_{optical} =$ %.4f' % Pol_tot_opt, fontsize=10)
#ax0.text(1E7,0.1,'$\Pi_{radio} =$ %.4f' % Pol_tot_rad, fontsize=10)
#ax0.text(1E7,0.45,'$\Pi_{x-ray} =$ %.4f' % Pol_tot_gam, fontsize=10)
#ax0.text(100,0.2,'$\Pi^{Initial}_{tot} =$ %.5f' % Pol_Init_tot, fontsize=10)

#ax0.legend()
#subplot for the EVPA
ax05 = plt.subplot(gs[1], sharex = ax0)
line05 = ax05.plot(freqtoeV(fq_mids[1:39]), EVPA[1:39]*180/np.pi, color='g',label='EVPA Non Rotating')
line055 = ax05.plot(freqtoeV(fq_mids3[1:39]), EVPA2[1:39]*180/np.pi,'--' ,color='g',label='EVPA Rotating')
#ax05.set_yscale("log")
ax05.set_xlim([1E-6, 1E13])
#ax05.set_ylim([-1.58, 1.58])
ax05.set_ylim([-90, 90])
yticks = ax05.yaxis.get_major_ticks()
#ax05.text(1E7,-0.65*180/np.pi,'$EVPA_{optical} =$ %.4f' % (EVPA_tot_opt*180/np.pi), fontsize=10)
##ax05.text(1E7,-1.2*180/np.pi,'$EVPA_{radio} =$ %.4f' % (EVPA_tot_rad*180/np.pi), fontsize=10)
#ax05.text(1E7,-0.1*180/np.pi,'$EVPA_{x-ray} =$ %.4f' % (EVPA_tot_gam*180/np.pi), fontsize=10)
ax05.set_ylabel(r'$\theta_{sky}[deg]$', size='13')
#ax05.legend()
#the second subplot
# shared axis X
ax1 = plt.subplot(gs[2], sharex = ax0)
line3 = ax1.plot(freqtoeV(fq_mids_IC), P_detected_IC, 'r-', label='SSC')#Inverse Compton
line3 = ax1.plot(freqtoeV(fq_mids_IC3), P_detected_IC3, 'r--', label='SSC')#Inverse Compton
line4 = ax1.plot(freqtoeV(fq_mids), P_detected, 'b-', label='Synchrotron') #synchrotron
line4 = ax1.plot(freqtoeV(fq_mids3), P_detected3,'b--', label='Synchrotron') #synchrotron
line5 = ax1.plot(ph_energy_MK501, flux_MK501, 'k.', label='Data 2008-2010')
#line6 = ax1.plot(ph_energy_MK421[0:30], flux_MK421[0:30], 'g.', label='data 2008-2009')
#line7 = ax1.plot(ph_energy[0:30], flux_BL[0:30], 'r.', label='data 2008-2009')
ax1.set_yscale('log')
ax1.set_ylim([1E-14, 9.99E-10])
ax1.set_xlim([1E-5, 1E10])
ax1.set_ylabel(r'$\nu F_{\nu}$ [erg s$^{-1}$ cm$^{-2}$]', size='13')
ax1.set_xlabel('eV', size='13')
#ax1.legend()
plt.setp(ax0.get_xticklabels(), visible=False)
plt.setp(ax05.get_xticklabels(), visible=False)
# remove last tick label for the second subplot
yticks = ax1.yaxis.get_major_ticks()
yticks[-2].label1.set_visible(False)

plt.subplots_adjust(hspace=.0)
#fig.savefig('/Users/ALP/Desktop/DD'+str(10)+'.png')
plt.show()