__author__ = 'ALP'
# This interpretation assumes that the jet material moves at ~ c so you only see one electron population evolving along the conical jet.
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import math as math
from jet_fns import *
import copy as copy
import sys

keydat = np.loadtxt('keyparams.txt')
Basics = np.loadtxt('basicdata.txt')
#-------------------Load in Polarisation---------------------------------
P_paraArray = np.loadtxt('Pparafile.txt')
P_perpArray = np.loadtxt('Pperpfile.txt')
Pol_single = np.loadtxt('Pol_single.txt') #for comparison
### First calculate initial electron PL population polarisation to compare with full dynamic pop.
Pol_Init = np.zeros(50)
for i,ting in enumerate(P_perpArray[0,:]):
    if (ting > 0 ) or (P_paraArray[0,i]>0):
        Pol_Init[i] = (ting - P_paraArray[0,i])/(ting + P_paraArray[0,i])
    else:
        Pol_Init[i] = 1.0
### Then polarisation for each frequency interval for whole jet emission
P_perp = P_perpArray.sum(axis=0)
P_para = P_paraArray.sum(axis=0)
Pol = np.zeros(50)
for i,ting in enumerate(P_perp):
    if (ting > 0 ) or (P_para[i]>0):
        Pol[i] = (ting - P_para[i])/(ting + P_para[i])
    else:
        Pol[i] = 1.0
###Now Total polarisation over all frequencies for both init and whole jet
Pol_Init_tot = (P_perpArray[0,8:40].sum() - P_paraArray[0,8:40].sum())/(P_perpArray[0,8:40].sum() + P_paraArray[0,8:40].sum()) #should = alpha+1/alpha+7/3
Pol_tot = (P_perp.sum() - P_para.sum())/(P_perp.sum() + P_para.sum())

# initialise EVPA
EVPA = []

#find time taken to traverse first section dx in days, this will be Poisson increment:
t0 = (Basics[0,0] / (3E8)) / (60*60*24)
lamb = 7
Poiss = [np.random.poisson(lamb) for i in range(31)]
T = [t0*i for i in range(sum(Poiss))]

phi = 0
R = Basics[0,3]
theta_obs=keydat[2]#2.8
c =1000  #Basics[0,0]/R
print(c)

k = 0
for i,Pois in enumerate(Poiss):
    theta_B = math.atan2(math.cos(k*t0/50 + phi),((c/R)*math.sin(theta_obs)-math.cos(theta_obs)*math.sin(k*t0/50 +phi))) #maybe use arctan2
    for j in range(Pois):
        EVPA.append(theta_B * 180/math.pi)
    k += Pois



plt.figure()
plt.plot(T, EVPA, 'r--', label='EVPA')#Inverse Compton
#plt.plot([1E10, 1E10], [1E-14, 1E-10], 'k-', lw=2)
#plt.plot([3E11, 3E11], [1E-14, 1E-10], 'k-', lw=1)
#plt.xscale('log')
#plt.yscale('log')
plt.ylabel(r'$EVPA [deg]$', size='14')
plt.xlabel('Days', size='14')
#plt.ylim(1E-14, 1E-9)
#plt.xlim(1E-6, 1E14)#plt.xlim(1E-6, 1E12)
#plt.legend()
plt.yticks(size='12')
plt.xticks(size='12')
plt.show()