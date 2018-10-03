import numpy as np
import matplotlib.pyplot as plt
from jet_fns import *

p_perpIC = np.loadtxt('PperpfileIC.txt')
p_paraIC = np.loadtxt('PparafileIC.txt')
p_perp = np.loadtxt('Pperpfile.txt')
p_para = np.loadtxt('Pparafile.txt')
elex = np.loadtxt('elecpop.txt')
ph = np.loadtxt('photonpop.txt')


keydat = np.loadtxt('keyparams.txt')
opdat = np.loadtxt('EBLOpacity.txt') #high energy gamma opacities for different z due to EBL

#loading up data points for different blazars (need to change distances and redshift as well!!!)
d_Blazar = 0.276E9*3.08E18 # pc in cm #276E6pc BL_Lac, 144E6 Mkn501, 131E6 Mkn421, 891E6 J2143
z = 0.0686 #0.0686 BL-Lac, 0.034 Mkn501, 0.031 Mkn421, 0.211 J2143
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


plt.subplot(411)
plt.plot(freqtoeV(fq_mids_IC),(p_perpIC-p_paraIC)/(p_perpIC+p_paraIC),label='IC')
plt.plot(freqtoeV(fq_mids),(p_perp-p_para)/(p_perp+p_para),label='sync')
plt.legend()
plt.xscale('log')
plt.title('Polarization (Sync + IC)')

plt.subplot(412)
plt.plot(elex[:,0],elex[:,4])
plt.xscale('log')
plt.yscale('log')
plt.title('Electron N(E)')

plt.subplot(413)
plt.plot(ph[:,0],ph[:,1])
plt.xscale('log')
plt.yscale('log')
plt.title('Synchrotron Photon number density')

plt.subplot(414)
plt.plot(ph[:,0],ph[:,2])
plt.xscale('log')
plt.yscale('log')
plt.title('Synchrotron Photon energy density')
