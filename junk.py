import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import math as math
from jet_fns import *
import copy as copy


keydat = np.loadtxt('keyparams.txt')
opdat = np.loadtxt('EBLOpacity.txt') #high energy gamma opacities for different z due to EBL

#loading up data points for different blazars (need to change distances and redshift as well!!!)
d_Blazar = 144E6*3.08E18 # pc in cm #276E6pc BL_Lac, 144E6 Mkn501, 131E6 Mkn421, 891E6 J2143
z = 0.034 #0.0686 BL-Lac, 0.034 Mkn501, 0.031 Mkn421, 0.211 J2143
theta_obs=keydat[2]#2.8
BLLac_pts = np.loadtxt('BLLacpointsfinal.data')
Mkn501_pts = np.loadtxt('MKN501.txt')
Mkn421_pts = np.loadtxt('MKN421.txt')
ph_energy = BLLac_pts[:,0]
flux_BL = BLLac_pts[:,1]
ph_energy_MK501 = Mkn501_pts[:,0]
flux_MK501 = Mkn501_pts[:,1]
ph_energy_MK421 = Mkn421_pts[:,0]
flux_MK421 = Mkn421_pts[:,1]
L_jet=keydat[0]#5.0E20
gamma_bulk=keydat[1]#7.5#12.0
beta_bulk=(1.0-(gamma_bulk**(-2.0)))**(0.5)
doppler_factor = 1.0/(gamma_bulk*(1.0-beta_bulk*np.cos(np.deg2rad(theta_obs))))


#define the bins
frdata = np.loadtxt('freqrange.txt')*doppler_factor#*doppler_factor#*gamma_bulk*2.8
fq_mins = frdata[:,0]
fq_maxs = frdata[:,1]
fq_mids = frdata[:,2]
fq_mins_IC = frdata[:,4]
fq_maxs_IC = frdata[:,5]
fq_mids_IC = frdata[:,6]

pi = np.loadtxt("TESTFIL2.txt")
Pol = pi[:,0]
EVPA = pi[:,1]
P_detected = pi[:,2]*1.0E7*(1.0/((4.0*np.pi*d_Blazar**2.0)*(1.0+z)**2.0)) #convert flux to ergs per cm and incluse blazar distance


fig = plt.figure(2)
# set height ratios for sublots
gs = gridspec.GridSpec(3, 1, height_ratios=[1,1,2])
# the fisrt subplot
ax0 = plt.subplot(gs[0])
# log scale for axis X of the first subplot
ax0.set_title('$W_j=3*10^{37}W$, $L_{jet} = 5*10^{20}m$, $B_0=8*10^{-5}T$, $E_{max}=50*10^9eV$,\n' + r'$\alpha=1.95$, $\theta_{opening} = 9^{\circ}$, $\theta_{obs}=4.0^{\circ}$, $\gamma_{bulk}=7.5$')
ax0.set_xscale("log")
ax0.set_xlim([1E-6, 1E13])
ax0.set_ylim([0, 1.0])
ax0.set_ylabel(r'$\Pi(\omega)$', size='13')
line0 = ax0.plot(freqtoeV(fq_mids), Pol,'m',label='Pol Fraction')
#line01 = ax0.plot(freqtoeV(fq_mids_IC), PolIC,'r',label='Pol Fraction IC')
#line1 = ax0.plot(freqtoeV(fq_mids), Pol_Init,'r',label='Initial Population')
#line2 = ax0.plot(freqtoeV(fq_mids), Pol_single, color='g',label='Single e')
#ax0.text(1E7,0.3,'$\Pi_{optical} =$ %.4f' % Pol_tot_opt, fontsize=10)
#ax0.text(1E7,0.1,'$\Pi_{radio} =$ %.4f' % Pol_tot_rad, fontsize=10)
#ax0.text(1E7,0.45,'$\Pi_{x-ray} =$ %.4f' % Pol_tot_gam, fontsize=10)
#ax0.text(100,0.2,'$\Pi^{Initial}_{tot} =$ %.5f' % Pol_Init_tot, fontsize=10)

ax0.legend()
#subplot for the EVPA
ax05 = plt.subplot(gs[1], sharex = ax0)
line05 = ax05.plot(freqtoeV(fq_mids), np.array(EVPA)*180/np.pi, color='g',label='EVPA')
#line051 = ax05.plot(freqtoeV(fq_mids_IC), np.array(EVPAIC)*180/np.pi, color='r',label='EVPAIC')
#ax05.set_yscale("log")
ax05.set_xlim([1E-6, 1E13])
#ax05.set_ylim([-1.58, 1.58])
ax05.set_ylim([-90, 90])
yticks = ax05.yaxis.get_major_ticks()
#ax05.text(1E7,-0.65*180/np.pi,'$EVPA_{optical} =$ %.4f' % (EVPA_tot_opt*180/np.pi), fontsize=10)
#ax05.text(1E7,-1.2*180/np.pi,'$EVPA_{radio} =$ %.4f' % (EVPA_tot_rad*180/np.pi), fontsize=10)
#ax05.text(1E7,-0.1*180/np.pi,'$EVPA_{x-ray} =$ %.4f' % (EVPA_tot_gam*180/np.pi), fontsize=10)
ax05.set_ylabel(r'$\theta_{sky}[deg]$', size='13')
ax05.legend()
#the second subplot
# shared axis X
ax1 = plt.subplot(gs[2], sharex = ax0)
#line3 = ax1.plot(freqtoeV(fq_mids_IC), P_detected_rawIC, 'r-.', label='ICRAW')#Inverse Compton
line4 = ax1.plot(freqtoeV(fq_mids), P_detected, 'b-', label='synchrotron') #synchrotron
#line5 = ax1.plot(freqtoeV(fq_mids), P_detected_raw, 'b-.', label='synchrotronRAW') #synchrotron
line6 = ax1.plot(ph_energy_MK501[0:30], flux_MK501[0:30], 'k.', label='data 2008-2009')
#line6 = ax1.plot(ph_energy_MK421[0:30], flux_MK421[0:30], 'g.', label='data 2008-2009')
#line7 = ax1.plot(ph_energy[0:30], flux_BL[0:30], 'r.', label='data 2008-2009')
ax1.set_yscale('log')
ax1.set_ylim([1E-14, 9.99E-10])
ax1.set_xlim([1E-6, 1E13])
ax1.set_ylabel(r'$\nu F_{\nu}$ [erg s$^{-1}$ cm$^{-2}$]', size='13')
ax1.set_xlabel('eV', size='13')
ax1.legend()
plt.setp(ax0.get_xticklabels(), visible=False)
# remove last tick label for the second subplot
yticks = ax1.yaxis.get_major_ticks()
yticks[-2].label1.set_visible(False)

plt.subplots_adjust(hspace=.0)
#fig.savefig('/Users/ALP/Desktop/DD'+str(10)+'.png')
plt.show()

