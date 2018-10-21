import numpy as np
import matplotlib.pyplot as plt

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

pi = np.loadtxt("TESTFIL2.txt")
plt.figure(1)
plt.subplot(311)
plt.xscale("log")
plt.title(r"$\Pi$")
plt.plot(fq_mids,pi[:50,0])
plt.plot(fq_mids_IC,pi[:50,3])
# plt.plot(fq_mids,pi[50:100,0],'r')
# plt.plot(fq_mids_IC,pi[50:100,3],'r')
plt.subplot(312)
plt.xscale("log")
plt.title("EVPA")
plt.plot(fq_mids,pi[:50,1])
plt.plot(fq_mids_IC,pi[:50,4])
# plt.plot(fq_mids,pi[50:100,1],'r')
# plt.plot(fq_mids_IC,pi[50:100,4],'r')
plt.subplot(313)
plt.xscale("log")
plt.yscale("log")
plt.title(r"$P v \omega$")
plt.plot(fq_mids,pi[:50,2])
plt.plot(fq_mids_IC,pi[:50,5])
# plt.plot(fq_mids,pi[50:100,2],'r')
# plt.plot(fq_mids_IC,pi[50:100,5],'r')


