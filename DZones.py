#Functions to work out how many zones dominate at each frequency along the SED

import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import math as math
from jet_fns import *
from matplotlib import gridspec
from scipy.signal import savgol_filter
from scipy.stats import truncnorm

def get_truncated_normal(mean=0, sd=1, low=0, upp=10):
    return truncnorm(
        (low - mean) / sd, (upp - mean) / sd, loc=mean, scale=sd)

sed = np.loadtxt("singleSED2.txt") #load template SED
sed[:, 0] = freqtoeV(sed[:, 0])
def real1ZoneSED(sed, nu, numax, D):
    #numerical single zone SED function from model
    if isinstance(nu,np.ndarray):
        flux = np.zeros(len(nu))
        for i,n in enumerate(nu):
            flux[i] = D**4 * (6.1E11**0.15 / (655914 * numax)**0.15) * sed[np.argmin(np.abs(D*sed[:,0]-n)),1] * np.exp(-n/(D*numax))
        return flux
    else:
        return D**4 * (6.1E11**0.15 / (655914 * numax)**0.15) * sed[np.argmin(np.abs(D*sed[:,0]-nu)),1] * np.exp(-nu/(D*numax))

def singleZoneSED(nu, numax, p, D, const,n):
    #analytical function for a single zone's SED, frequency nu, frequency where exp starts to kick in numax,
    #spectral index p, doppler factor D and const to match flux with real SED
    #numax > 100
    sumexp = 0
    for i in range(n):
        if i < n/10:
            sumexp += const * np.exp(-nu/(D * numax/(10*np.log(0.1*i+1)+1))) / ((i+1)**(-1/2)) #cooling + B field diminishing (probably should be logarithmic instead of linear)
            prev = (i+1)**(-1/2)
            previ = i
        else:
            sumexp += const * np.exp(-nu / (D * numax / (10 * np.log(0.1*i + 1) + 1))) / (((i + (prev)**(2)/1 - previ - 1/1) + 1) ** (1/2))
    return D**4 * nu ** -p * sumexp

def Ddist(gamma, theta_obs, N_0, theta_op, plot=False):
    #gives distribution of Doppler factors of zones, outputs a list with each zone and its D
    #assumes theta_op = 0.2/gamma atm
    #assume theta_obs >= theta_op for now
    beta = np.sqrt(1-1/gamma**2)
    #theta_op = 0.3 / gamma #0.2 from papers yannis
    r_1 = theta_op / np.sqrt(N_0)
    N_D = {}
    theta = 0
    sum_counter = 0

    for i in range(int(np.sqrt(N_0)) + 3):
        if (theta_obs > theta_op):
            theta = 2 * r_1 * i + theta_obs - theta_op
        else:
            theta = 2 * r_1 * i
        D = 1 / (gamma * (1 - beta * np.cos(theta)))
        if i == 0:
            N_D[D] = 1
            sum_counter += 1
        else:
            N_D[D] = np.ceil(circle_intersection(theta, theta_op, theta_obs)[0] / (2*r_1)) #how much of circle theta is in circle theta_op
            sum_counter += np.ceil(circle_intersection(theta, theta_op, theta_obs)[0] / (2*r_1))

    if plot:
        plt.plot(N_D.keys(), N_D.values())
        plt.show()

    return N_D, sum_counter

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def Ddist_gauss(gamma, theta_obs, N_0, plot=False):
    #gives distribution of Doppler factors of zones, outputs a list with each zone and its D
    #assumes theta_op = 0.2/gamma atm
    #assume theta_obs >= theta_op for now
    beta = np.sqrt(1-1/gamma**2)
    Doppler = 1/(gamma*(1-beta*np.cos(theta_obs)))
    N_D = {}
    sum_counter = 0

    for i in np.linspace(0.2,1.9,13)*Doppler:
        N_D[i] = np.ceil(gaussian(i,Doppler,Doppler/3) * 8)
        sum_counter += np.ceil(gaussian(i,Doppler,Doppler/3) * 8)

    if plot:
        plt.plot(N_D.keys(), N_D.values())

    return N_D, sum_counter


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def zoneDominator(nup,N_0,gamma,theta_obs, theta_op, B, Emax=np.array([None]), Donly=False, Dlist=None): #returns how many zones dominate at a given nu + expected polarization
    beta = np.sqrt(1 - 1/gamma**2)
    if Dlist == None:
        N_D, _ = Ddist(gamma, theta_obs, N_0, theta_op)
        Dlist = np.array([[k for i in range(int(v))] for k,v in N_D.items()])
        Dlist = np.array([item for sublist in Dlist for item in sublist])

    nu = 10**np.linspace(-3.5,8,100) #have to shift these forward with respect to SED freqs due to D factor
    Fnu = np.zeros((len(Dlist),100))
    psi = np.zeros(len(Dlist)) #polarization angles
    Stokes = np.zeros(3)
    # plt.figure()
    # D_prev = -1
    for i,D in enumerate(Dlist):
        psi[i] = np.arctan2(B[i,1],B[i,0])
        if Donly:
            numax = 10**3.12
        elif (Emax.any() != None):
            numax = Emax[i] * 10**3.45
        else:
            #projection along normal of transverse jet shock
            numax = ((B[i,2]*np.cos( np.arccos((np.cos(theta_obs) - beta)/(1-beta*np.cos(theta_obs))) ) + B[i,0]*np.sin( np.arccos((np.cos(theta_obs) - beta)/(1-beta*np.cos(theta_obs))) ))**2 / (B[i,1]**2 + B[i,2]**2 + B[i,0]**2)) * 10**3.45 #( Bpara / B )**2 controls gamma_max / numax -> marscher
        if numax < 10: #gamma_min
            numax = 10
        #Fnu[i,:] = singleZoneSED(nu, numax, -0.5, D, 1, 100)
        Fnu[i, :] = real1ZoneSED(sed, nu, numax, D) #* pow(sin(acos(B_effectives[g][g][2])),(effective_alpha[l]+1)/2) should add q_theta
        #if D != D_prev:
        #    plt.plot(nu, Fnu[i, :]*N_D[D], '--',label=str(int(N_D[D])) + ' at D = ' +str("%.2f" % D))
        # D_prev = D

    Fnu_tot = np.sum(Fnu,axis=0)
    nupeak = np.argmax(Fnu_tot)
    target = find_nearest(nu, nup * nu[nupeak])

    Stokes[0] = Fnu_tot[target]
    for i,D in enumerate(Dlist):
        Stokes[1] += Fnu[i,target] * np.cos(2*psi[i]) * 0.75 * min(1.333,max(np.log(1E11*nu[target] / (D * numax))/np.log(1E11),1)) #increase polarization on rolloff
        Stokes[2] += Fnu[i,target] * np.sin(2*psi[i]) * 0.75 * min(1.333,max(np.log(1E11*nu[target] / (D * numax))/np.log(1E11),1)) #deviation from power law

    # # print(nupeak)
    # # print(Fnu_tot[nupeak])
    # plt.plot(nu,Fnu_tot,'-b')
    # for i in range(len(Dlist)):
    #      plt.plot(nu,Fnu[i,:],'--')
    # plt.xscale("log")
    # plt.yscale("log")
    # plt.xlabel(r"$h\nu [eV]$")
    # plt.ylabel(r"Observed Flux [arbitrary]")
    # plt.legend()
    # plt.show()

    flux_counter = 0
    num_counter = 0
    fluxes = np.sort(Fnu[:,target])[::-1]
    for i,f in enumerate(fluxes):
        flux_counter += f
        num_counter += 1
        if flux_counter >= 0.5 * Fnu_tot[target]:
            return num_counter, Stokes


def N_zone_plot(N_0,gamma,theta_obs, theta_op):
    _, N0_real = Ddist(gamma, theta_obs, N_0, theta_op)
    print("N0_real = ", N0_real)
    num = 100
    nup = 10 ** np.linspace(-1, 5, 60) #peak is at idx 30, nup is nu/nupeak
    N_total = np.zeros(60)
    Pi_total = np.zeros(60)
    PA_total = np.zeros(60)
    for i in range(num):
        B =np.random.normal(size=(int(N0_real),3))
        Emax = np.array([None])
        for j,nu in enumerate(nup):
            numz, stoke = zoneDominator(nu, N_0, gamma, theta_obs, theta_op, B, Emax, Donly=True)
            N_total[j] += numz
            Pi_total[j] += np.sqrt(stoke[1] ** 2 + stoke[2] ** 2) / stoke[0]
            PA_total[j] += 0.5 * np.arctan2(stoke[2], stoke[1])

    N_av = N_total / num
    Pi_av = Pi_total / num
    PA_av = PA_total / num
    #np.savetxt('Donly1.txt',np.array(N)/N[30])
    #plt.plot(nup, np.array(N)/N[30],'r')
    #plt.xscale("log")
    #plt.yscale("log")
    #plt.show()
    return N_av / N_av[10], Pi_av, PA_av


def N_zone_av(N_0,gamma,theta_obs, theta_op, Marsch=True):
    _, N0_real = Ddist(gamma, theta_obs, N_0, theta_op)
    print("N0_real = ", N0_real)
    num = 100
    nup = 10 ** np.linspace(-1, 5, 60)
    N_total = np.zeros(60)
    Pi_total = np.zeros(60)
    PA_total = np.zeros(60)
    random_variable = get_truncated_normal(mean=1, sd=1, low=-1, upp=1)

    for i in range(num):
        B = np.random.normal(size=(int(N0_real), 3))
        if Marsch:
            Emax = np.array([None])
        else:
            Emax = np.random.random(size=int(N0_real))**2 #np.random.lognormal(sigma=2.5, size=int(N0_real)) #random_variable.rvs(int(N0_real))**2

        for j,nu in enumerate(nup):
            numz, stoke = zoneDominator(nu, N_0, gamma, theta_obs,theta_op, B, Emax, Donly=False)
            N_total[j] += numz
            Pi_total[j] += np.sqrt(stoke[1]**2 + stoke[2]**2) / stoke[0]
            PA_total[j] += 0.5 * np.arctan2(stoke[2],stoke[1])


    N_av = N_total / num
    Pi_av = Pi_total / num
    PA_av = PA_total / num
    #np.savetxt('Monly1.txt',N_av/N_av[30])
    #plt.plot(nup, N_av/N_av[30],'b')
    #plt.xlabel(r"$\nu / \nu_{peak}$")
    #plt.ylabel(r"$N / N_{peak}$")
    #plt.xscale("log")
    #plt.yscale("log")
    #plt.show()

    return N_av/N_av[10], Pi_av, PA_av





def circle_intersection(r_1,r_2,d): #find length of overlapping arc between two circles, d < r_1 + r_2 and d > |r_1 - r_2|
    #need r_1 > r_2 ! Gives length of arc in bigger circle if smaller one is inside
    if (d > r_1 + r_2):
        return (0,0) #left is length of r_1 in r_1, right is vice-versa
    elif (d < abs(r_1 - r_2)):
        if (r_1 < r_2):
            return (2*np.pi * r_1,0)
        else:
            return (0,2*np.pi * r_2)
    else:
        return r_1 * 2*np.arccos((r_1**2 + d**2 - r_2**2) / (2 * r_1 * d)), r_2 * 2 * np.arccos((r_2 ** 2 + d ** 2 - r_1 ** 2) / (2 * r_2 * d))

def main():
    #sed = np.loadtxt("singleSED.txt") #load template SED
    #sed[:, 0] = freqtoeV(sed[:, 0])
    nup = 10 ** np.linspace(-1, 5, 60)

    Donly, Donly_pi, Donly_pa = N_zone_plot(150, 10, 0.15/10, theta_op=0.2/10)
    Donly2, Donly_pi2, Donly_pa2 = N_zone_plot(150, 10, 0.5/10, theta_op=0.5/10)
    Marsch, Marsch_pi, Marsch_pa = N_zone_av(150, 10, 0.3/10, theta_op=0.3/10)
    #Marsch2, Marsch_pi2, Marsch_pa2 = N_zone_av(100, 10, 0.6/10, theta_op=0.6/10)
    Romani, Romani_pi, Romani_pa = N_zone_av(150, 10, 0.3/10, theta_op=0.3/10, Marsch=False)

    np.savez("Pol_Zones", nup = nup, Donly=Donly, Donly_pi=Donly_pi, Donly2=Donly2, Donly_pi2=Donly_pi2, Marsch=Marsch, Marsch_pi=Marsch_pi, Romani=Romani, Romani_pi=Romani_pi)

    fig = plt.figure(1)
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 2])
    ax0 = plt.subplot(gs[0])
    ax0.set_xscale("log")
    #ax0.set_xlim([1E-5, 1E17])
    ax0.set_ylim([0, 1.0])
    ax0.set_ylabel(r'$\Pi$', size='13')
    #ax05 = plt.subplot(gs[1], sharex=ax0)
    #ax05.set_xlim([1E-5, 1E17])
    #ax05.set_ylim([-90, 90])
    #yticks = ax05.yaxis.get_major_ticks()
    #ax05.set_ylabel(r'$\theta_{sky}[deg]$', size='13')
    ax1 = plt.subplot(gs[1], sharex=ax0)
    ax1.set_yscale('log')
    #ax1.set_ylim([1E-17, 9.99E-11])
    #ax1.set_xlim([1E-5, 1E17])
    ax1.set_ylabel(r"$N / N_{peak}$", size='13')
    ax1.set_xlabel(r'$\nu / \nu_{peak}$', size='13')
    plt.setp(ax0.get_xticklabels(), visible=False)
    # remove last tick label for the second subplot
    yticks = ax1.yaxis.get_major_ticks()
    yticks[-2].label1.set_visible(False)

    plt.subplots_adjust(hspace=.0)

    line0 = ax0.plot(nup, Donly_pi, 'b', label='Pol Fraction')
    #line1 = ax05.plot(nup, np.array(Donly_pa) * 180 / np.pi, 'b', label='EVPA')
    line2 = ax1.plot(nup, Donly, 'b', label='NZones')

    line0 = ax0.plot(nup, Donly_pi2, 'b--', label='Pol Fraction')
    # line1 = ax05.plot(nup, np.array(Donly_pa) * 180 / np.pi, 'b', label='EVPA')
    line2 = ax1.plot(nup, Donly2, 'b--', label='NZones')

    line3 = ax0.plot(nup, Marsch_pi, 'r', label='Pol Fraction')
    #line4 = ax05.plot(nup, np.array(Marsch_pa) * 180 / np.pi, 'r', label='EVPA')
    line5 = ax1.plot(nup, Marsch, 'r', label='NZones')

    #line3 = ax0.plot(nup, Marsch_pi2, 'r--', label='Pol Fraction')
    # line4 = ax05.plot(nup, np.array(Marsch_pa) * 180 / np.pi, 'r', label='EVPA')
    #line5 = ax1.plot(nup, Marsch2, 'r--', label='NZones')

    line6 = ax0.plot(nup, Romani_pi, 'g', label='Pol Fraction')
    #line7 = ax05.plot(nup, np.array(Romani_pa) * 180 / np.pi, 'g', label='EVPA')
    line8 = ax1.plot(nup, Romani, 'g', label='NZones')

    plt.show()


    # #plt.plot(nup, Marsch, 'r')
    # plt.plot(nup, Donly, 'b')
    # plt.xlabel(r"$\nu / \nu_{peak}$")
    # plt.ylabel(r"$N / N_{peak}$")
    # plt.xscale("log")
    # plt.yscale("log")
    # plt.show()



if __name__ == '__main__':
    main()