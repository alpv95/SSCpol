import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec
from jet_fns import *

X = np.loadtxt('Edens0.5_19_100/Xfile.txt')
Y = np.loadtxt('Edens0.5_19_100/Yfile.txt')
f = np.loadtxt('Edens0.5_19_100/freqrange1.txt')
keydat = np.loadtxt('Edens0.5_19_100/keyparams1.txt')

d_Blazar = 1000E6*3.08E18 # pc in cm #276E6pc BL_Lac, 144E6 Mkn501, 131E6 Mkn421, 891E6 J2143, 1000E6 J0721, 1800E6 3C279
z = 0.3 #0.0686 BL-Lac, 0.034 Mkn501, 0.031 Mkn421, 0.211 J2143, 0.3 J0721, 0.536 3C279

theta_obs=keydat[2]
W_j=keydat[0]#5.0E20
gamma_bulk=keydat[1]#7.5#12.0

theta_open_p = keydat[3]
alpha = keydat[4]
B0 = keydat[5]
E_max = keydat[6]
n_blocks = keydat[7]
array_size = int(keydat[8]) #50, can change this here by hand

beta_bulk=(1.0-(gamma_bulk**(-2.0)))**(0.5)
doppler_factor = 1.0/(gamma_bulk*(1.0-beta_bulk*np.cos(np.deg2rad(theta_obs))))
f *= doppler_factor
f = freqtoeV(f)

plt.figure(1)
gs = gridspec.GridSpec(1, 1, height_ratios=[1])
# the fisrt subplot
ax0 = plt.subplot(gs[0])
ax0.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False,
                     bottom=True, top=False, left=True, right=True)
ax0.tick_params(axis="x", labelsize=12)
ax0.tick_params(axis="y", labelsize=12)
# log scale for axis X of the first subplot
#ax0.set_title('$W_j= %.2G W$, $B_0= %.2G T$, $E_{max}= %.2G eV$,\n' % (W_j, B0, E_max)
   #           + r'$\alpha= %.2f $, $\theta_{opening} = %.1f ^{\circ}$, $\theta_{obs}= %.1f ^{\circ}$,'
    #            r'$\gamma_{bulk}= %.1f $, $n_{blocks}= %d $' % (
     #         alpha, 180 / np.pi * np.arctan(np.tan(np.pi * theta_open_p / 180) / gamma_bulk), theta_obs,
      #        gamma_bulk, n_blocks))
ax0.set_xscale("log")
#ax0.set_xlim([1E-5, 1E17])
#ax0.set_ylim([0, 1.0])
ax0.set_ylabel(r'$P_{jet} [W]$', size='15')
ax0.set_xlabel(r'$h\nu [eV]$', size='15')
ax0.plot(f[:,2],np.sum(X,axis=0),'r')
ax0.plot(f[:,2],np.sum(Y,axis=0),'b--')

ax1 = ax0.twiny()
ax1.plot((1.6*10**-19 * f[:,2]) / (6.63*10**-34), np.sum(X,axis=0),'r')
ax1.tick_params(axis="x", labelsize=12)
ax1.set_xlabel(r'$\nu [Hz]$', size='15')
#ax05 = plt.subplot(gs[1], sharex=ax0)
#ax05.set_xlim([1E-5, 1E17])
#ax05.set_ylim([-90, 90])
#yticks = ax05.yaxis.get_major_ticks()
#ax05.set_ylabel(r'$\theta_{sky}[deg]$', size='13')
#ax1 = plt.subplot(gs[1], sharex=ax0)
ax0.set_yscale('log')
ax1.set_xscale('log')
ax1.set_yscale('log')
#ax1.set_ylim([0, 1.1])
#ax1.set_xlim([1E-5, 1E17])
#ax1.set_ylabel(r'$\rho_{self} / \rho_{total}$', size='13')
#ax1.set_xlabel(r'$x$ [m]', size='13')
#ax1.plot(x[:,1],e[:,6],'k')
#ax1.plot(x[:,1],e[:,7],'k--')
#plt.setp(ax0.get_xticklabels(), visible=False)
# remove last tick label for the second subplot
#yticks = ax1.yaxis.get_major_ticks()
#yticks[-2].label1.set_visible(False)

#plt.subplots_adjust(hspace=.0)
plt.show()
