import numpy as np
import matplotlib
#matplotlib.use("TkAgg")
from matplotlib import pyplot as plt
from matplotlib import gridspec
from scipy.signal import savgol_filter
import math as math
from jet_fns import *
import matplotlib.cm as cm
import copy as copy

def Norm_factor(alpha,flux,Erange=(2,10)): #Erange must be a tuple (keV) eg, (2,10), alpha is spectral index across range,
    #assumes alpha is positive!!!
    #flux in ergs cm^-2 s^-1
    #calculates Norm for ixpeobssim
    A = flux * (alpha-1) / (Erange[0]**(1-alpha) - Erange[1]**(1-alpha))
    return A / (1.6E-16 * 1E7) #convert to counts norm

def Dfactor(theta, gamma): #doppler factor
    beta = np.sqrt(1 - 1/gamma**2)
    return 1 / (gamma * (1 - beta * np.cos(np.deg2rad(theta))))

def Theta_lab(theta_op_jet, gamma): #theta_opening in lab frame
    return np.rad2deg(np.arctan(np.tan(np.deg2rad(theta_op_jet))/gamma))

def Theta_jet(theta_op_lab, gamma): #theta_opening in jet frame
    return np.rad2deg(np.arctan(np.tan(np.deg2rad(theta_op_lab))*gamma))

def plot_SED(filename,keyfile,freqfile,IC=True, save=False): #plots SED with polarisation fraction and EVPA as a function of energy
    #can plot
    #filename can be single string or list of strings to plot simultaneously
    if isinstance(filename, str):
        filename = [filename]
        legends = 1

    keydat = np.loadtxt(keyfile)
    if keydat.ndim == 1:
        keydat = np.expand_dims(keydat,axis=0)
    #opdat = np.loadtxt('EBLOpacity.txt') #high energy gamma opacities for different z due to EBL

    #loading up data points for different blazars (need to change distances and redshift as well!!!)
    d_Blazar = 1000E6 * 3.08E18 # pc in cm #1000E6 J0211, 1627E6 S5, 1789E6 TXS
    z = 0.2 #0.2 J0211, 0.3365 TXS, 0.31 S5
    theta_obs=keydat[0,2]#2.8
    J0211_pts = np.loadtxt('data/new_data_sed_CGRaBSJ0211+1051_XMM.txt')
    J0211_pts2 = np.loadtxt('data/new_data_sed_CGRaBSJ0211+1051.txt')
    J0211_pts = np.concatenate([J0211_pts,J0211_pts2])
    S5_pts = np.loadtxt('data/tofit_flare_0716+016.txt')
    S5_pts2 = np.loadtxt('data/tofit_historical_0716+016.txt')
    S5_pts2 = np.insert(S5_pts2, 1, 0,axis=1)
    S5_pts = np.concatenate([S5_pts,S5_pts2])
    TXS_pts = np.loadtxt('data/to_fit_SED_txs0506+056.txt')
    TXS_pts = np.insert(TXS_pts, 1, 0, axis=1)
    TXS_pts2 = np.loadtxt('data/txs_xray_data.txt')
    TXS_pts = np.concatenate([TXS_pts,TXS_pts2])
    # BLLac_pts = np.loadtxt('BLLacpointsfinal.data')
    #Mkn501_pts = np.loadtxt('MKN501.txt')
    #Mkn421_pts = np.loadtxt('MKN421.txt')
    #C279_pts = np.loadtxt("3C279.txt")
    #C279_pts = C279_pts[:,[0,2]]
    #S50716_pts = np.loadtxt("S50716714.txt")
    #S50716_pts = S50716_pts[:,[0,2]]
    #ph_energy = BLLac_pts[:,0]
    #flux_BL = BLLac_pts[:,1]
    #ph_energy_MK501 = Mkn501_pts[:,0]
    #flux_MK501 = Mkn501_pts[:,1]
    #ph_energy_MK421 = Mkn421_pts[:,0]
    #flux_MK421 = Mkn421_pts[:,1]

    W_j=keydat[0,0]#5.0E20
    gamma_bulk=keydat[0,1]#7.5#12.0

    theta_open_p = keydat[0,3]
    alpha = keydat[0,4]
    B0 = keydat[0,5]
    E_max = keydat[0,6]
    n_blocks = keydat[0,7]
    array_size = int(keydat[0,8]) #50, can change this here by hand

    beta_bulk=(1.0-(gamma_bulk**(-2.0)))**(0.5)
    doppler_factor = 1.0/(gamma_bulk*(1.0-beta_bulk*np.cos(np.deg2rad(theta_obs))))


    #define the bins
    frdata = np.loadtxt(freqfile)*doppler_factor#*doppler_factor#*gamma_bulk*2.8
    fq_mins = frdata[:array_size,0]
    fq_maxs = frdata[:array_size,1]
    fq_mids = frdata[:array_size,2]
    fq_mins_IC = frdata[:array_size,4]
    fq_maxs_IC = frdata[:array_size,5]
    fq_mids_IC = frdata[:array_size,6]

    fig = plt.figure(figsize=(7,7))
    # set height ratios for sublots
    gs = gridspec.GridSpec(3, 1, height_ratios=[1, 1, 2])
    # the fisrt subplot
    ax0 = plt.subplot(gs[0])
    # log scale for axis X of the first subplot
    #ax0.set_title('$W_j= %.2G W$, $B_0= %.2G T$, $E_{max}= %.2G eV$,\n' % (W_j, B0, E_max)
    #              + r'$\alpha= %.2f $, $\theta_{opening} = %.1f ^{\circ}$, $\theta_{obs}= %.1f ^{\circ}$,'
    #                r'$\gamma_{bulk}= %.1f $, $n_{blocks}= %d $' % (
    #              alpha, 180 / np.pi * np.arctan(np.tan(np.pi * theta_open_p / 180) / gamma_bulk), theta_obs,
    #              gamma_bulk, n_blocks))
    ax0.set_xscale("log")
    ax0.set_xlim([1E-7, 1E18])
    ax0.set_ylim([0, 1.0])
    ax0.set_ylabel(r'$\Pi(\omega)$', size='15')
    ax0.tick_params(labelbottom=False, labeltop=False, labelleft=True, labelright=False,
                     bottom=True, top=True, left=True, right=True)
    ax0.tick_params(axis="y", labelsize=12)

    ax05 = plt.subplot(gs[1], sharex=ax0)
    ax05.set_xlim([1E-7, 1E18])
    ax05.set_ylim([-90, 90])
    ax05.tick_params(labelbottom=False, labeltop=False, labelleft=True, labelright=False,
                     bottom=True, top=True, left=True, right=True)
    ax05.tick_params(axis="y", labelsize=12)
    yticks = ax05.yaxis.get_major_ticks()
    ax05.set_ylabel(r'$\theta_{sky}[deg]$', size='15')
    ax1 = plt.subplot(gs[2], sharex=ax0)
    ax1.set_yscale('log')
    ax1.set_ylim([1E-16, 9.99E-10])
    ax1.set_xlim([1E-7, 1E18])
    ax1.set_ylabel(r'$\nu F_{\nu}$ [erg s$^{-1}$ cm$^{-2}$]', size='15')
    ax1.set_xlabel(r'$h\nu$ [eV]', size='15')
    ax1.tick_params(axis="x", labelsize=12)
    ax1.tick_params(axis="y", labelsize=12)
    ax1.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False,
                     bottom=True, top=True, left=True, right=True)
    plt.setp(ax0.get_xticklabels(), visible=False)
    # remove last tick label for the second subplot
    yticks = ax1.yaxis.get_major_ticks()
    yticks[-2].label1.set_visible(False)

    plt.subplots_adjust(hspace=.0)

    for j,f in enumerate(filename):
        fullpi = np.loadtxt(f)
        #fullpi = fullpi[(29*array_size):(29*array_size)+array_size,:]
        # check0 = fullpi==0 #replacing 0s with nans for smart binning
        # fullpi[check0] = np.nan
        pi = np.zeros((array_size,3))
        stdpi = np.zeros((array_size,3))
        pi_IC = np.zeros((array_size, 3))
        stdpi_IC = np.zeros((array_size, 3))
        pi_ICS = np.zeros((array_size, 3))
        stdpi_ICS = np.zeros((array_size, 3))
        n_examples = int(fullpi.shape[0]/array_size)

        for i in range(n_examples): #calculates average pi/evpa/power over many jet realisations
            #i=j*20
            pi[:,0] += fullpi[(i*array_size):(i*array_size)+array_size,0]
            pi[:,1] += fullpi[(i*array_size):(i*array_size)+array_size,1]
            pi[:,2] += fullpi[(i*array_size):(i*array_size)+array_size,2]

            pi_IC[:, 0] += fullpi[(i * array_size):(i * array_size) + array_size, 3]
            pi_IC[:, 1] += fullpi[(i * array_size):(i * array_size) + array_size, 4]
            pi_IC[:, 2] += fullpi[(i * array_size):(i * array_size) + array_size, 5]

            pi_ICS[:, 0] += fullpi[(i * array_size):(i * array_size) + array_size, 6]
            pi_ICS[:, 1] += fullpi[(i * array_size):(i * array_size) + array_size, 7]
            pi_ICS[:, 2] += fullpi[(i * array_size):(i * array_size) + array_size, 8]

        pi = pi / n_examples
        pi_IC = pi_IC / n_examples
        pi_ICS = pi_ICS / n_examples

        for i in range(n_examples): #calculates std of pi/evpa/power over many jet realisations
            #i=j*20
            stdpi[:,0] += (fullpi[(i*array_size):(i*array_size)+array_size,0] - pi[:,0])**2
            stdpi[:,1] += (fullpi[(i*array_size):(i*array_size)+array_size,1] - pi[:,1])**2
            stdpi[:,2] += ((fullpi[(i*array_size):(i*array_size)+array_size,2] - pi[:,2])*1.0E7*(1.0/((4.0*np.pi*d_Blazar**2.0)*(1.0+z)**2.0)))**2

            stdpi_IC[:, 0] += (fullpi[(i * array_size):(i * array_size) + array_size, 3] - pi_IC[:, 0]) ** 2
            stdpi_IC[:, 1] += (fullpi[(i * array_size):(i * array_size) + array_size, 4] - pi_IC[:, 1]) ** 2
            stdpi_IC[:, 2] += ((fullpi[(i * array_size):(i * array_size) + array_size, 5] - pi_IC[:, 2]) * 1.0E7 * (1.0 / ((4.0 * np.pi * d_Blazar ** 2.0) * (1.0 + z) ** 2.0)))**2

            stdpi_ICS[:, 0] += (fullpi[(i * array_size):(i * array_size) + array_size, 6] - pi_ICS[:, 0]) ** 2
            stdpi_ICS[:, 1] += (fullpi[(i * array_size):(i * array_size) + array_size, 7] - pi_ICS[:, 1]) ** 2
            stdpi_ICS[:, 2] += ((fullpi[(i * array_size):(i * array_size) + array_size, 8] - pi_ICS[:, 2]) * 1.0E7 * (
                         1.0 / ((4.0 * np.pi * d_Blazar ** 2.0) * (1.0 + z) ** 2.0))) ** 2

        stdpi = np.sqrt(stdpi / (n_examples))
        stdpi_IC = np.sqrt(stdpi_IC / (n_examples))
        stdpi_ICS = np.sqrt(stdpi_ICS / (n_examples))


        Pol = pi[:,0]
        EVPA = pi[:,1]
        P_detected = pi[:,2]*1.0E7*(1.0/((4.0*np.pi*d_Blazar**2.0)*(1.0+z)**2.0)) #convert flux to ergs per cm and incluse blazar distance
        Pol_IC = pi_IC[:, 0]
        EVPA_IC = pi_IC[:, 1]
        P_detected_IC = pi_IC[:, 2] * 1.0E7 * (1.0 / ((4.0 * np.pi * d_Blazar ** 2.0) * (1.0 + z) ** 2.0))
        Pol_ICS = pi_ICS[:, 0]
        EVPA_ICS = pi_ICS[:, 1]
        P_detected_ICS = pi_ICS[:, 2] * 1.0E7 * (1.0 / ((4.0 * np.pi * d_Blazar ** 2.0) * (1.0 + z) ** 2.0))

        line0 = ax0.plot(freqtoeV(fq_mids), Pol,'b',label='Pol Fraction')
        if (n_examples > 1):
            ax0.fill_between(freqtoeV(fq_mids), Pol - stdpi[:,0], Pol + stdpi[:,0], alpha=0.2, facecolor='b')
            ax0.fill_between(freqtoeV(fq_mids_IC), Pol_IC - stdpi_IC[:,0], Pol_IC + stdpi_IC[:,0], alpha=0.2, facecolor='r')
            #line001 = ax0.plot(freqtoeV(fq_mids), stdpi[:,0],'b--',label='Pol Fraction std')
            #line01 = ax0.plot(freqtoeV(fq_mids_IC), stdpi_IC[:,0], 'r--',label='Pol IC std')
           
        line0143 = ax0.plot(freqtoeV(fq_mids_IC[P_detected_IC!=0.0]), Pol_IC[P_detected_IC!=0.0],'r',label='Pol Fraction IC',ls='-.')
        line02 = ax0.plot(freqtoeV(fq_mids), Pol_ICS, 'k', ls='-',label='Pol Fraction ICS')
        #line1 = ax0.plot(freqtoeV(fq_mids), Pol_Init,'r',label='Initial Population')
        #line2 = ax0.plot(freqtoeV(fq_mids), Pol_single, color='g',label='Single e')
        #ax0.text(1E7,0.3,'$\Pi_{optical} =$ %.4f' % Pol_tot_opt, fontsize=10)
        #ax0.text(1E7,0.1,'$\Pi_{radio} =$ %.4f' % Pol_tot_rad, fontsize=10)
        #ax0.text(1E7,0.45,'$\Pi_{x-ray} =$ %.4f' % Pol_tot_gam, fontsize=10)
        #ax0.text(100,0.2,'$\Pi^{Initial}_{tot} =$ %.5f' % Pol_Init_tot, fontsize=10)

        #subplot for the EVPA

        line05 = ax05.plot(freqtoeV(fq_mids), EVPA*180/np.pi, color='b',label='EVPA')
        if (n_examples > 1):
            ax05.fill_between(freqtoeV(fq_mids), EVPA*180/np.pi - stdpi[:,1]*180/np.pi, EVPA*180/np.pi + stdpi[:,1]*180/np.pi, alpha=0.2, facecolor='b')
            ax05.fill_between(freqtoeV(fq_mids_IC), EVPA_IC*180/np.pi - stdpi_IC[:,1]*180/np.pi, EVPA_IC*180/np.pi + stdpi_IC[:,1]*180/np.pi, alpha=0.2, facecolor='r')
            #line005 = ax05.plot(freqtoeV(fq_mids), stdpi[:,1]*180/np.pi, color='b',linestyle='--',label='EVPA std')
            #line0053 = ax05.plot(freqtoeV(fq_mids_IC), stdpi_IC[:, 1] * 180 / np.pi, color='r', linestyle='--', label='EVPAIC std')
        line051 = ax05.plot(freqtoeV(fq_mids_IC[P_detected_IC!=0.0]), EVPA_IC[P_detected_IC!=0.0]*180/np.pi, color='r',label='EVPAIC', ls='-.')
        line052 = ax05.plot(freqtoeV(fq_mids), EVPA_ICS * 180 / np.pi, color='k',ls='-', label='EVPAICS')
        #the second subplot
        # shared axis X

        #savgol filter to smooth bumpy ICS from rebinning larger IC bins in synchrotron ones.
        #np.savetxt("singleSED.txt", np.array([freqtoeV(fq_mids),P_detected]))
        line3 = ax1.plot(freqtoeV(fq_mids_IC[P_detected_IC!=0.0]), P_detected_IC[P_detected_IC!=0.0],linestyle='-.', color='r')#Inverse Compton color=(j/17,0.2,0.2)
        line355 = ax1.plot(freqtoeV(fq_mids), savgol_filter(P_detected_ICS,9,3), 'k')  # Inverse Compton + Synchrotron
        line4 = ax1.plot(freqtoeV(fq_mids), P_detected,'-', color='b') #synchrotron color=(0.1,j/17,1)
        if (n_examples > 1 ):
            ax1.fill_between(freqtoeV(fq_mids), P_detected - stdpi[:,2], P_detected + stdpi[:,2], alpha=0.2, facecolor='b')
            ax1.fill_between(freqtoeV(fq_mids_IC[P_detected_IC!=0.0]), P_detected_IC[P_detected_IC!=0.0] - stdpi_IC[:,2][P_detected_IC!=0.0], 
            P_detected_IC[P_detected_IC!=0.0] + stdpi_IC[:,2][P_detected_IC!=0.0], alpha=0.2, facecolor='r')
            # line004 = ax1.plot(freqtoeV(fq_mids), stdpi[:,2], 'b--', label='synchrotron std')
            # line006 = ax1.plot(freqtoeV(fq_mids_IC), stdpi_IC[:, 2], 'r--', label='IC std')
        #line5 = ax1.plot(freqtoeV(fq_mids), P_detected_raw, 'b-.', label='synchrotronRAW') #synchrotron
        #line6 = ax1.errorbar(10**S5_pts[:,0], 10**S5_pts[:,2], yerr=np.stack([np.abs(10**S5_pts[:,2] - 10**(S5_pts[:,2] - S5_pts[:,3])),np.abs(10**S5_pts[:,2] - 10**(S5_pts[:,2] + S5_pts[:,3]))]), color='k', marker='o', label='S50716',ls="", markersize=1.2)
        #line6 = ax1.errorbar(10**TXS_pts[:,0], 10**TXS_pts[:,2], yerr=np.stack([np.abs(10**TXS_pts[:,2] - 10**(TXS_pts[:,2] - TXS_pts[:,3])),np.abs(10**TXS_pts[:,2] - 10**(TXS_pts[:,2] + TXS_pts[:,3]))]), color='g', marker='o', label='TXS0506',ls="", markersize=1.2)
        line6 = ax1.errorbar(10**J0211_pts[:,0], 10**J0211_pts[:,2], yerr=np.stack([np.abs(10**J0211_pts[:,2] - 10**(J0211_pts[:,2] - J0211_pts[:,3])),np.abs(10**J0211_pts[:,2] - 10**(J0211_pts[:,2] + J0211_pts[:,3]))]), color='r', marker='o', label='J0211',ls="", markersize=1.2)
        ax1.axvline(10**3,ls=':',color='k')
        ax1.axvline(10**4,ls=':',color='k')
        ax05.axvline(10**3,ls=':',color='k')
        ax05.axvline(10**4,ls=':',color='k')
        ax0.axvline(10**3,ls=':',color='k')
        ax0.axvline(10**4,ls=':',color='k')


        #line6 = ax1.plot(10**(S50716_pts[:,0]), 10**(S50716_pts[:,1]), 'k.', label='observation')
        #line6 = ax1.plot(ph_energy_MK421, flux_MK421, 'g.', label='data 2008-2009')
        #line7 = ax1.plot(ph_energy[0:30], flux_BL[0:30], 'r.', label='data 2008-2009')

        if (legends):
            # ax0.legend()
            # ax05.legend()
            ax1.legend()
    if save:
       fig.savefig("SED.pdf", bbox_inches='tight')
    plt.show()

def PlotAndSaveB(count): #Visualises jet zones with their B-fields given Proj_Bfile, for making movies

    keydat = np.loadtxt('keyparams.txt')
    Zones = np.loadtxt("Proj_Bfile.txt")
    zone_loc = Zones[:,:3]
    B_thetas = Zones[:,3:]
    n_blocks = Zones.shape[0]
    theta_obs = keydat[2] * np.pi/180
    gamma = keydat[1]
    theta_open =  np.arctan(np.tan(keydat[3]*np.pi/180)/gamma)

    X = np.zeros((n_blocks,2)) #position of blocks on diagram
    Y = np.zeros((n_blocks,2))
    U = np.zeros((n_blocks,2)) #B-field direction of blocks on diagram
    V = np.zeros((n_blocks,2))

    fig = plt.figure(1)
    ax = fig.add_subplot(111)
    ax.plot(-theta_obs,0,'gx')

    circle1 = plt.Circle((-theta_obs, 0), 1/gamma, color='g',fill=False,linestyle='--')
    ax.add_artist(circle1)
    ax.set_ylabel("Y [rad]")
    ax.set_xlabel("X [rad]")
    ax.set_aspect('equal', adjustable='box')
    colors = np.array([0.7,0.1])
    colormap = cm.inferno
    for i in range(n_blocks): #position of blocks on diagram
        X[i,0] = zone_loc[i,0] * np.cos(zone_loc[i,1])
        Y[i,0] = zone_loc[i,0] * np.sin(zone_loc[i,1])
        U[i,0] = np.cos(B_thetas[i,3]) #first random set
        V[i,0] = np.sin(B_thetas[i,3])
        # X[i,1] = 100 * zone_loc[i,0] * np.cos(zone_loc[i,1])
        # Y[i,1] = 100 * zone_loc[i,0] * np.sin(zone_loc[i,1])
        # U[i,1] = np.cos(B_thetas[i,1]) #first random set
        # V[i,1] = np.sin(B_thetas[i,1])
        X[i,1] = zone_loc[i,0] * np.cos(zone_loc[i,1])
        Y[i,1] = zone_loc[i,0] * np.sin(zone_loc[i,1])
        U[i,1] = np.cos(B_thetas[i,4]) #first random set
        V[i,1] = np.sin(B_thetas[i,4])
        ax.quiver(X[i], Y[i],U[i],V[i],color=colormap(colors))
        ax.plot(X[i,0],Y[i,0],'k.')

    if count <10:
        plt.savefig("movie00"+str(count)+".png")
    elif (count <100):
        plt.savefig("movie0"+str(count)+".png")
    else:
        plt.savefig("movie"+str(count)+".png")

    plt.close(fig)

def Save_movie(count): #joins IXPE pol_ang & pol_frac plots with movie
    period = 1/5.79E-8 #period of viewing time in s
    xp = np.loadtxt("IXPEobsim.txt")
    model = np.loadtxt("TESTFIL2.txt")
    model = np.rad2deg(model[:,5])+90
    model = np.array([h-50 if h>=50 else h+130 for h in model]) #align with jet axis along the sky

    keydat = np.loadtxt('keyparams.txt')
    Zones = np.loadtxt("Proj_Bfile.txt")
    zone_loc = Zones[:,:3]
    B_thetas = Zones[:,3:]
    n_blocks = Zones.shape[0]
    theta_obs = keydat[2] * np.pi/180
    gamma = keydat[1]
    theta_open =  np.arctan(np.tan(keydat[3]*np.pi/180)/gamma)

    X = np.zeros((n_blocks,2)) #position of blocks on diagram
    Y = np.zeros((n_blocks,2))
    U = np.zeros((n_blocks,2)) #B-field direction of blocks on diagram
    V = np.zeros((n_blocks,2))

    fig = plt.figure(figsize=(12, 6))
    gs = gridspec.GridSpec(2,1)
    gs.update(top=0.93,bottom=0.07,right=0.999,left=0.36, hspace=0.)
    # the fisrt subplot
    ax1 = plt.subplot(gs[0])
    ax1.set_ylim(0,0.4)
    ax1.set_ylabel(r"$\Pi$", fontsize="8")
    ax1.errorbar(xp[:,0]*period/(24*60*60), xp[:,2], xp[:,3], xp[:,1]*period/(24*60*60), fmt='o',
                 label="IXPEobs", markersize=3.2)
    ax1.plot([(count/200)*period/(24*60*60),(count/200)*period/(24*60*60)],[0,50],'r')
    ax1.plot([-10,-10],[-20,-21],'m--',label="Jet Projection") #just for legend
    ax1.plot([-10,-10],[-20,-21],'r.-',label="Helix Only",alpha=0.48) #just for legend
    ax1.legend(prop={'size': 7})
    for tick in ax1.yaxis.get_major_ticks():
                tick.label.set_fontsize(8)

    ax2 = plt.subplot(gs[1],sharex=ax1)
    ax2.set_ylim(-10,210)
    ax2.set_xlim(-5,205)
    ax2.set_xlabel(r"Time [Days]", fontsize="8")
    ax2.set_ylabel(r"PA [$^{\circ}$]", fontsize="8")
    ax2.errorbar(xp[:,0]*period/(24*60*60), xp[:,4], xp[:,5], xp[:,1]*period/(24*60*60), fmt='o',
                  markersize=3.2)
    ax2.plot([(count/200)*period/(24*60*60),(count/200)*period/(24*60*60)],[-10,210],'r')
    ax2.plot([-5,205],[50,50],'m--',label="Jet Projection")
    ax2.plot(np.array([i+40 for i in range(64)])*period/(24*60*60*200),model,'r.-',label="Jet Projection",alpha=0.48)
    ax2.plot(np.array([i+151 for i in range(32)])*period/(24*60*60*200),np.append(model[46:],model[:14]),'r.-',label="Jet Projection",alpha=0.48)
    #ax2.legend()
    for tick in ax2.yaxis.get_major_ticks():
                tick.label.set_fontsize(8)
    for tick in ax2.xaxis.get_major_ticks():
                tick.label.set_fontsize(8)

    gs2 = gridspec.GridSpec(1,1)
    gs2.update(top=0.8,bottom=0.2,right=0.3,left=0.07,hspace=0.)
    ax = plt.subplot(gs2[0])
    ax.plot(-theta_obs*np.cos(np.deg2rad(140)),-theta_obs*np.sin(np.deg2rad(140)),'gx')

    circle1 = plt.Circle((-theta_obs*np.cos(np.deg2rad(140)), -theta_obs*np.sin(np.deg2rad(140))), 1/gamma, color='g',fill=False,linestyle='--')
    ax.add_artist(circle1)
    ax.set_ylabel("Y [rad]", fontsize="8")
    ax.set_xlabel("X [rad]", fontsize="8")
    ax.set_aspect('equal', adjustable='box')
    colors = np.array([0.7,0.1])
    colormap = cm.inferno
    for i in range(n_blocks): #position of blocks on diagram
        X[i,0] = zone_loc[i,0] * np.cos(zone_loc[i,1]) * np.cos(np.deg2rad(140)) - zone_loc[i,0] * np.sin(zone_loc[i,1]) * np.sin(np.deg2rad(140))
        Y[i,0] = zone_loc[i,0] * np.sin(zone_loc[i,1]) * np.cos(np.deg2rad(140)) + zone_loc[i,0] * np.cos(zone_loc[i,1]) * np.sin(np.deg2rad(140))
        U[i,0] = np.cos(B_thetas[i,3]) * np.cos(np.deg2rad(140)) - np.sin(B_thetas[i,3]) * np.sin(np.deg2rad(140)) #No DD
        V[i,0] = np.sin(B_thetas[i,3]) * np.cos(np.deg2rad(140)) + np.sin(np.deg2rad(140)) * np.cos(B_thetas[i,3])
        X[i,1] = zone_loc[i,0] * np.cos(zone_loc[i,1]) * np.cos(np.deg2rad(140)) - zone_loc[i,0] * np.sin(zone_loc[i,1]) * np.sin(np.deg2rad(140))
        Y[i,1] = zone_loc[i,0] * np.sin(zone_loc[i,1]) * np.cos(np.deg2rad(140)) + zone_loc[i,0] * np.cos(zone_loc[i,1]) * np.sin(np.deg2rad(140))
        U[i,1] = np.cos(B_thetas[i,4]) * np.cos(np.deg2rad(140)) - np.sin(B_thetas[i,4]) * np.sin(np.deg2rad(140)) #DD
        V[i,1] = np.sin(B_thetas[i,4]) * np.cos(np.deg2rad(140)) + np.sin(np.deg2rad(140)) * np.cos(B_thetas[i,4])
        ax.quiver(X[i], Y[i],U[i],V[i],color=colormap(colors))
        ax.plot(X[i,0],Y[i,0],'k.',markersize="3.75")

    for tick in ax.yaxis.get_major_ticks():
                tick.label.set_fontsize(8)
    for tick in ax.xaxis.get_major_ticks():
                tick.label.set_fontsize(8)
    plt.setp(ax1.get_xticklabels(), visible=False)
    # plt.show()
    if count <10:
        plt.savefig("movie00"+str(count)+".png")
    elif (count <100):
        plt.savefig("movie0"+str(count)+".png")
    else:
        plt.savefig("movie"+str(count)+".png")

    plt.close(fig)

if __name__ == "__main__":
    #Save_movie(1)
    plot_SED(["TESTFIL2.txt","TESTFIL4.txt","TESTFIL3.txt","TESTFIL5.txt","TESTFIL1.txt"],IC=True)



