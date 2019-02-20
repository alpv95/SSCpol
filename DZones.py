#Functions to work out how many zones dominate at each frequency along the SED

import numpy as np
import matplotlib
from matplotlib import pyplot as plt
import math as math
from jet_fns import *
from scipy.signal import savgol_filter



sed = np.loadtxt("singleSED.txt")
def real1ZoneSED(sed, nu, numax, D):
    #numerical single zone SED function from model
    if isinstance(nu,np.ndarray):
        flux = np.zeros(len(nu))
        for i,n in enumerate(nu):
            flux[i] = D**4 * sed[1,np.argmin(np.abs(D*sed[0,:]-n))] * np.exp(-n/(D*numax))
        return flux
    else:
        return D**4 * sed[1,np.argmin(np.abs(D*sed[0,:]-nu))] * np.exp(-nu/(D*numax))

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

def Ddist(gamma, theta_obs, N_0, plot=False):
    #gives distribution of Doppler factors of zones, outputs a list with each zone and its D
    #assumes theta_op = 0.2/gamma atm
    #assume theta_obs >= theta_op for now
    beta = np.sqrt(1-1/gamma**2)
    theta_op = 0.3 / gamma #0.2 from papers yannis
    r_1 = theta_op / np.sqrt(N_0)
    N_D = {}
    theta = 0
    sum_counter = 0

    for i in range(int(np.sqrt(N_0)) + 1):
        theta = 2 * r_1 * i + theta_obs - theta_op
        D = 1 / (gamma * (1 - beta * np.cos(theta)))
        if i == 0:
            N_D[D] = 1
            sum_counter += 1
        else:
            N_D[D] = np.ceil(circle_intersection(theta, theta_op, theta_obs) / (2*r_1))
            sum_counter += np.ceil(circle_intersection(theta, theta_op, theta_obs) / (2*r_1))

    if plot:
        plt.plot(N_D.keys(), N_D.values())

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

def zoneDominator(nup,N_0,gamma,theta_obs, Dlist=None): #returns how many zones dominate at a given nu
    if Dlist == None:
        N_D, _ = Ddist(gamma, theta_obs, N_0)
        Dlist = np.array([[k for i in range(int(v))] for k,v in N_D.items()])
        Dlist = np.array([item for sublist in Dlist for item in sublist])

    nu = 10**np.linspace(-3,7,100)
    Fnu = np.zeros((len(Dlist),100))
    # plt.figure()
    # D_prev = -1
    for i,D in enumerate(Dlist):
        numax = 0.334*10**4.0
        #B = np.random.random(3) #is this really isotropic?
        #numax = (B[2]**2 / (B[1]**2 + B[2]**2 + B[0]**2)) * 10**4.0 #( Bpara / B )**2 controls gamma_max / numax -> marscher
        if numax < 10: #gamma_min
            numax = 10
        #Fnu[i,:] = singleZoneSED(nu, numax, -0.5, D, 1, 100)
        Fnu[i, :] = real1ZoneSED(sed, nu, numax, D)
        # if D != D_prev:
        #     plt.plot(nu, Fnu[i, :]*N_D[D], '--',label=str(int(N_D[D])) + ' at D = ' +str("%.2f" % D))
        # D_prev = D

    Fnu_tot = np.sum(Fnu,axis=0)
    nupeak = np.argmax(Fnu_tot)
    target = find_nearest(nu, nup * nu[nupeak])
    # print(nupeak)
    # print(Fnu_tot[nupeak])
    # plt.plot(nu,Fnu_tot,'-b')
    # for i in range(len(Dlist)):
    #     plt.plot(nu,Fnu[i,:],'--')
    # plt.xscale("log")
    # plt.yscale("log")
    # plt.xlabel(r"$\nu [eV]$")
    # plt.ylabel(r"Flux [arbitrary]")
    # #plt.legend()
    # plt.show()

    flux_counter = 0
    num_counter = 0
    fluxes = np.sort(Fnu[:,target])[::-1]
    for i,f in enumerate(fluxes):
        flux_counter += f
        num_counter += 1
        if flux_counter >= 0.5 * Fnu_tot[target]:
            return num_counter


def N_zone_plot(N_0,gamma,theta_obs):
    _, N0_real = Ddist(gamma, theta_obs, N_0)
    print("N0_real = ", N0_real)
    nup = 10 ** np.linspace(-3, 5, 80) #peak is at idx 30
    N = []
    for nu in nup:
        N.append(zoneDominator(nu, N_0, gamma, theta_obs,))
    np.savetxt('Donly.txt',np.array(N)/N[30])
    plt.plot(nup, np.array(N)/N[30],'r')
    plt.xscale("log")
    plt.yscale("log")
    plt.show()


def N_zone_av(N_0,gamma,theta_obs):
    _, N0_real = Ddist(gamma, theta_obs, N_0)
    print("N0_real = ", N0_real)
    nup = 10 ** np.linspace(-3, 4, 70)
    N_total = np.zeros(70)
    for i in range(70):
        for j,nu in enumerate(nup):
            N_total[j] += zoneDominator(nu, N_0, gamma, theta_obs,)

    N_av = N_total / 70
    np.savetxt('Monly.txt',N_av/N_av[30])
    plt.plot(nup, N_av/N_av[30],'b')
    plt.xlabel(r"$\nu / \nu_{peak}$")
    plt.ylabel(r"$N / N_{peak}$")
    plt.xscale("log")
    plt.yscale("log")

    return N_av/N0_real





def circle_intersection(r_1,r_2,d): #find length of overlapping arc between two circles, d < r_1 + r_2 and d > |r_1 - r_2|
    #need r_1 > r_2 ! Gives length of arc in bigger circle if smaller one is inside
    return r_1 * 2*np.arccos((r_1**2 + d**2 - r_2**2) / (2 * r_1 * d))


# if __name__ == "__main__":
#