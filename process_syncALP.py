
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec
import math as math
from jet_fns import *
import copy as copy

#ROWS (x) are values for each electron bin
#columns (y) are values for each jet section

#N_ICbins = 50

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


fcritdata = np.loadtxt('critfreqs.txt')*doppler_factor#*frequency gets doppler boosted from jet rest frame to observers frame!!!
opdata = np.loadtxt('kvalues.txt')
powdata = np.loadtxt('secpow.txt')*1.0E7*(1.0/((4.0*np.pi*d_Blazar**2.0)*(1.0+z)**2.0))#convert to flux in ergs
basics = np.loadtxt('basicdata.txt')
""" -------------------------------------"""
ICpowdata = np.loadtxt('ICfile.txt')*1.0E7*(1.0/((4.0*np.pi*d_Blazar**2.0)*(1.0+z)**2.0))#convert to flux in ergs
ICfreqdata = np.loadtxt('ICfreqfile.txt')*doppler_factor
'''--------------------------------------'''
##### Loading in Polarisation powers + theta for each block
P_paraArray = np.loadtxt('Pparafile.txt')
P_perpArray = np.loadtxt('Pperpfile.txt')
block_theta = np.loadtxt('block_thetafile.txt')
Pol_single = np.loadtxt('Pol_single.txt') #for comparison
### First calculate initial electron PL population polarisation to compare with full dynamic pop.
Pol_Init = np.zeros(50)
for i,ting in enumerate(P_perpArray[0,:]):
    if (ting > 0 ) or (P_paraArray[0,i]>0):
        Pol_Init[i] = (ting - P_paraArray[0,i])/(ting + P_paraArray[0,i])
    else:
        Pol_Init[i] = 1.0
### Then polarisation for each frequency interval for whole jet emission
ProjB_theta = np.loadtxt('Proj_Bfile.txt')

# find the doppler factor for each block:
doppler = 1.0/(gamma_bulk*(1.0-beta_bulk*np.cos(block_theta)))
'''
P_perp = P_perpArray.sum(axis=0) * (1/np.shape(ProjB_theta)[1]) * (doppler**4).sum() #have to include this average doppler weighting for denominator
P_para = P_paraArray.sum(axis=0) * (1/np.shape(ProjB_theta)[1]) * (doppler**4).sum()
'''

# form stokes vectors for each parallel and perp component of light in each section (same for all energies)
# list of tuples which represent 2 middle components of stokes vectors
Stokes_para = [[np.array([math.cos(2*projB),math.sin(2*projB)])for projB in projB_row]for projB_row in ProjB_theta]
Stokes_perp = [[np.array([math.cos(2*(projB-math.pi/2)),math.sin(2*(projB-math.pi/2))]) if projB >= -math.pi/2 else np.array([math.cos(2*(projB+3*math.pi/2)),math.sin(2*(projB+3*math.pi/2))]) for projB in ProjB_row] for ProjB_row in ProjB_theta]
Stokes_total = [np.zeros(2) for i in P_perpArray[0,:]]
Stokes_totaltot = np.zeros(2) #Stokes for whole spectrum
Stokes_total_radio = np.zeros(2) # Stokes for a band in the radio range
Stokes_total_radio_denominator = 0
Stokes_total_optical = np.zeros(2) #stokes for a band in the optical-xray range
Stokes_total_optical_denominator = 0
Stokes_total_gamma = np.zeros(2) #stokes for a band in the gamma range
Stokes_total_gamma_denominator = 0
Stokes_total_denominator = np.zeros(len(P_perpArray[0,:]))
'''
for j, item in enumerate(P_perpArray[0,:]): #loop over number of frequency bins
    if (freqtoeV(fq_mids[j])>1.5) and (freqtoeV(fq_mids[j])<3.4):
        Stokes_total_optical_denominator += (P_perp[j]+P_para[j])
    elif (freqtoeV(fq_mids[j])<1E-3) and (freqtoeV(fq_mids[j])>1E-4):
        Stokes_total_radio_denominator += (P_perp[j]+P_para[j])
    elif (freqtoeV(fq_mids[j])<1E4) and (freqtoeV(fq_mids[j])>1E3):
        Stokes_total_gamma_denominator += (P_perp[j]+P_para[j])

    for i, item2 in enumerate(P_perpArray[:,0]): #loop over number of jet segments
        for k,item3 in enumerate(ProjB_theta[0,:]): #loop over number of B-field blocks in jet segment

            Stokes_total[j] += doppler[k]**4 * (P_perpArray[i,j]*Stokes_perp[i][k] + P_paraArray[i,j]*Stokes_para[i][k]) / ((P_perp[j]+P_para[j])*np.shape(ProjB_theta)[1]) #dividing by number of blocks at the end here (assumes unifrm distr of electrons through jet segement)
            Stokes_totaltot += doppler[k]**4 * (P_perpArray[i,j]*Stokes_perp[i][k] + P_paraArray[i,j]*Stokes_para[i][k]) / ((P_perp.sum()+P_para.sum())*np.shape(ProjB_theta)[1])
            #now include Polarisation in 3 different energy bands: radio,optical,gamma
            if (freqtoeV(fq_mids[j])>1.5) and (freqtoeV(fq_mids[j])<3.4):
                Stokes_total_optical += doppler[k]**4 * (P_perpArray[i,j]*Stokes_perp[i][k] + P_paraArray[i,j]*Stokes_para[i][k]) / (np.shape(ProjB_theta)[1])
            elif (freqtoeV(fq_mids[j])<1E-3) and (freqtoeV(fq_mids[j])>1E-4):
                Stokes_total_radio += doppler[k]**4 * (P_perpArray[i,j]*Stokes_perp[i][k] + P_paraArray[i,j]*Stokes_para[i][k]) / (np.shape(ProjB_theta)[1])
            elif (freqtoeV(fq_mids[j])<1E4) and (freqtoeV(fq_mids[j])>1E3):
                Stokes_total_gamma += doppler[k]**4 * (P_perpArray[i,j]*Stokes_perp[i][k] + P_paraArray[i,j]*Stokes_para[i][k]) / (np.shape(ProjB_theta)[1])

'''
for i, item2 in enumerate(P_perpArray[:,0]): #loop over number of jet segments
    for k,item3 in enumerate(ProjB_theta[0,:]): #loop over number of B-field blocks in jet segment
        #if fq_mids[j]*(1 - doppler[k]/doppler_factor) > (fq_mids[j]-fq_mids[j-1])/2:
        for j, item in enumerate(P_perpArray[0,:]): #loop over number of frequency bins
            j_newd = (np.abs(fq_mids-fq_mids[j]*doppler[k]/doppler_factor)).argmin() #adjusting for different doppler factor of each block, might need to use log nearest?? if too  grainy

            if (freqtoeV(fq_mids[j_newd])>1.5) and (freqtoeV(fq_mids[j_newd])<3.4):
                Stokes_total_optical_denominator += (P_perpArray[i,j_newd] + P_paraArray[i,j_newd]) * (doppler[k]**4)
            elif (freqtoeV(fq_mids[j_newd])<1E-3) and (freqtoeV(fq_mids[j_newd])>1E-4):
                Stokes_total_radio_denominator += (P_perpArray[i,j_newd] + P_paraArray[i,j_newd]) * (doppler[k]**4)
            elif (freqtoeV(fq_mids[j_newd])<1E4) and (freqtoeV(fq_mids[j_newd])>1E3):
                Stokes_total_gamma_denominator += (P_perpArray[i,j_newd] + P_paraArray[i,j_newd]) * (doppler[k]**4)
            Stokes_total_denominator[j_newd] += (P_perpArray[i,j_newd] + P_paraArray[i,j_newd]) * (doppler[k]**4)


            Stokes_total[j_newd] += doppler[k]**4 * (P_perpArray[i,j_newd]*Stokes_perp[i][k] + P_paraArray[i,j_newd]*Stokes_para[i][k]) #dividing by number of blocks at the end here (assumes unifrm distr of electrons through jet segement)

            #now include Polarisation in 3 different energy bands: radio,optical,gamma
            if (freqtoeV(fq_mids[j_newd])>1.5) and (freqtoeV(fq_mids[j_newd])<3.4):
                Stokes_total_optical += doppler[k]**4 * (P_perpArray[i,j_newd]*Stokes_perp[i][k] + P_paraArray[i,j_newd]*Stokes_para[i][k])
            elif (freqtoeV(fq_mids[j_newd])<1E-3) and (freqtoeV(fq_mids[j_newd])>1E-4):
                Stokes_total_radio += doppler[k]**4 * (P_perpArray[i,j_newd]*Stokes_perp[i][k] + P_paraArray[i,j_newd]*Stokes_para[i][k])
            elif (freqtoeV(fq_mids[j_newd])<=1E4) and (freqtoeV(fq_mids[j_newd])>=1E3):
                Stokes_total_gamma += doppler[k]**4 * (P_perpArray[i,j_newd]*Stokes_perp[i][k] + P_paraArray[i,j_newd]*Stokes_para[i][k])


Stokes_total_optical = Stokes_total_optical / Stokes_total_optical_denominator
Stokes_total_radio = Stokes_total_radio / Stokes_total_radio_denominator
Stokes_total_gamma = Stokes_total_gamma / Stokes_total_gamma_denominator
for i,item in enumerate(Stokes_total):
    Stokes_total[i] = Stokes_total[i] / Stokes_total_denominator[i]

Pol = [math.sqrt(item[0]**2 + item[1]**2) for item in Stokes_total]
EVPA = [0.5*math.atan2(item[1],item[0]) for item in Stokes_total]

Pol_tot_opt = math.sqrt(Stokes_total_optical[0]**2 + Stokes_total_optical[1]**2)
EVPA_tot_opt = 0.5*math.atan2(Stokes_total_optical[1],Stokes_total_optical[0])
Pol_tot_rad = math.sqrt(Stokes_total_radio[0]**2 + Stokes_total_radio[1]**2)
EVPA_tot_rad = 0.5*math.atan2(Stokes_total_radio[1],Stokes_total_radio[0])
Pol_tot_gam = math.sqrt(Stokes_total_gamma[0]**2 + Stokes_total_gamma[1]**2)
EVPA_tot_gam = 0.5*math.atan2(Stokes_total_gamma[1],Stokes_total_gamma[0])


###Now Total polarisation over all frequencies for both init and whole jet
Pol_Init_tot = (P_perpArray[0,:].sum() - P_paraArray[0,:].sum())/(P_perpArray[0,:].sum() + P_paraArray[0,:].sum()) #should = alpha+1/alpha+7/3
'''Pol_tot = (P_perp.sum() - P_para.sum())/(P_perp.sum() + P_para.sum())'''

'''-----------------------------------------------'''
print(len(opdata[:,0]), len(opdata[0]))

dx_steps = basics[:,0]
x_cumu = basics[:,1]
B_steps = basics[:,2]
R_steps = basics[:,3]

print(len(dx_steps), len(x_cumu), len(B_steps), len(R_steps))

nfreqs = len(powdata[0]) #gives the number of frequency bins
nSecs = len(powdata[:,0]) #gives the number of jet sections

print('Computing ', nfreqs, ' frequencies...')
print('There are ', nSecs, ' jet sections.')

#create an array to store optical depths along the LoS for every frequency bin from every jet section
opt_depths = np.zeros(nSecs*nfreqs)
opt_depths = np.reshape(opt_depths, (nSecs, nfreqs)) #nSecs cols, nfreqs rows

print('new array shape is:', np.shape(opt_depths))
print('required shape is:', np.shape(fcritdata))

#print dx_steps

P_detected = np.zeros(len(fq_mins)) #store all emitted power here
P_detected_raw = np.zeros(len(fq_mins)) #as above for non opacity data



for i in range(nSecs):#nSecs): #loop over jet sections no of ROWS (column elements)

#    print 'working'

    B = B_steps[i]
    R = R_steps[i]
    dx = dx_steps[i]

    power_array = powdata[i]#emitted power in jet rest frame

    #just try to do the first section of the jet
    for j in range(nfreqs): #loop over all frequencies emitted by a jet section (row elements/no cols)
    #Power_array = powdata[i] #probably need this here for whole jet
#         print i, j
         freq_emit = fcritdata[i, j] #frequency at which power is emitted
         #print freq_emit #gives the expected values          

         index=0
         f_emit_max=0
         f_emit_min=0

         #PSUEDOCODE:
         #find the corresponding frequency bin for this emission
         for a in range(len(fq_mins)): #number of energy bins
             if freq_emit > fq_mins[a] and freq_emit <= fq_maxs[a]:
                  f_emit_max = fq_maxs[a]
                  f_emit_min = fq_mins[a]
                  index=a #changes index in detected power array  
                  
#         print "%.3e" % f_emit_min, "%.3e" % f_emit_max, "%.3e" %  freq_emit

         power_emit = powdata[i,j] #power at which this frequency was emitted
                  #optical_depth = 0 # initialise
         for k in range(nSecs-i): #loop over remaining jet sections 
             k+=i #ensure loops over correct section jet to end not from beginning
             if fcritdata[k, j] <= f_emit_max and fcritdata[k, j] > f_emit_min:
                 opt_depths[i,j] += opdata[k,j]*dx_steps[k]


P_rad_sec1 = np.array([])
#P_emit_sec1 = np.array([]) #jua=st the first row of fcrits

for i in range(nfreqs):
     seen_power = powdata[0,i]*np.exp(-opt_depths[0, i]*gamma_bulk**2*(1/np.cos(np.deg2rad(theta_obs))-beta_bulk))
     #seen_power = powdata[0,i]*np.exp(-opt_depths[0, i])*(1/np.cos(np.deg2rad(theta_obs)))
     P_rad_sec1 = np.append(P_rad_sec1, seen_power)

#plt.figure()
#plt.plot(B_steps[0:-1], fcritdata[:,0])
#plt.xscale('log')
#plt.yscale('log')
#plt.show()

flux = P_rad_sec1
flux_l = powdata[0]
max = np.max(flux)


obs_power = copy.deepcopy(fcritdata)

#some code to plot an entire synchrotron spectrum
#first get the power seen from each section
for i in range(nSecs):
    for j in range(nfreqs):
        #obs_power[i,j] = powdata[i, j]*np.exp(-opt_depths[i, j])*(1/(np.cos(np.deg2rad(theta_obs))))
        obs_power[i,j] = powdata[i,j]*np.exp(-opt_depths[i, j]*gamma_bulk**2*((1/np.cos(np.deg2rad(theta_obs)))-beta_bulk))

counts = 0
#now sum up the total power
for i in range(nSecs):
    for j in range(nfreqs):
        for k in range(len(fq_mins)):
            freq_emit = fcritdata[i, j] #frequency at which power is emitted                                      
            counts = 0
            if freq_emit > fq_mins[k] and freq_emit <= fq_maxs[k]:
                counts +=1 
#                print counts #tests for double counting
                P_detected[k] += obs_power[i,j]
                P_detected_raw[k] += powdata[i,j]

#Ptot_xerr=[np.log10(freqtoeV(fq_mids-fq_mins)), np.log10(freqtoeV(fq_maxs-fq_mids))]
Ptot_xerr=[(freqtoeV(fq_mids-fq_mins)), (freqtoeV(fq_maxs-fq_mids))]

'''__________________________________'''
P_detected_IC = np.zeros(len(fq_mins_IC))

for i in range(nSecs):
    for j in range(nfreqs):
        for k in range(len(fq_mins_IC)):
            freq_emit = ICfreqdata[i, j] #frequency at which power is emitted
            counts = 0
            if freq_emit > fq_mins_IC[k] and freq_emit <= fq_maxs_IC[k]:
                counts +=1
                #                print counts #tests for double counting
                P_detected_IC[k] += ICpowdata[i,j]

'''
                                                       # & AXIONS AXIONS ALPS!!!!
#THIS PART FOR WORKING OUT L_GAMMAS WITH EBL______________________________________________#
#Now also total radiant flux between 0.1 and 100GeV -> F_gamma:
F_gammaAx = np.zeros(10)
L_gammaAx = np.zeros(10) #Luminostiy between 0.1 and 100GeV
for i in range(len(P_detected_IC)):
    if (freqtoeV(fq_mids_IC[i]) >= 2E12 and freqtoeV(fq_mids_IC[i]) <= 1E13):
        for j in range(len(opdat)-1):
            if (freqtoeV(fq_mids_IC[i]) >= opdat[j,0] and freqtoeV(fq_mids_IC[i]) < opdat[j+1,0]):
                for k in range(len(F_gammaAx)):
                    F_gammaAx[k] += P_detected_IC[i] * (1/16 + (9/16)*math.exp(-opdat[j,k+1])) #need to add function that accounts for opacity due to EBL at particular redshift here, have to use freqtoeV(fq_mids_IC) to form this
            elif (freqtoeV(fq_mids_IC[i])>= opdat[-1,0]):
                 for k in range(len(F_gammaAx)):
                    F_gammaAx[k] += P_detected_IC[i] * (1/16 + (9/16)*math.exp(-opdat[-1,k+1]))

#convert flux between 0.1 and 100GeV back to Luminosity in ERGS
L_gammaAx = F_gammaAx *((4.0*np.pi*d_Blazar**2.0)*(1.0+z)**2.0)
print(F_gammaAx,L_gammaAx)
'''
'''
with open('RepAx>2<10Tev_theta0.txt', 'a') as f:
    f.write('\n%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf' % (L_gammaAx[0], L_gammaAx[1], L_gammaAx[2], L_gammaAx[3], L_gammaAx[4], L_gammaAx[5], L_gammaAx[6], L_gammaAx[7], L_gammaAx[8], L_gammaAx[9]))
'''
#___________________________________________________________________________________________#
'''
#THIS PART FOR WORKING OUT L_GAMMAS WITH EBL______________________________________________#
#Now also total radiant flux between 0.1 and 100GeV -> F_gamma:
F_gammaEBL = np.zeros(10)
L_gammaEBL = np.zeros(10) #Luminostiy between 0.1 and 100GeV
for i in range(len(P_detected_IC)):
    if (freqtoeV(fq_mids_IC[i]) >= 2E12 and freqtoeV(fq_mids_IC[i]) <= 1E13):
        for j in range(len(opdat)-1):
            if (freqtoeV(fq_mids_IC[i]) >= opdat[j,0] and freqtoeV(fq_mids_IC[i]) < opdat[j+1,0]):
                for k in range(len(F_gammaEBL)):
                    F_gammaEBL[k] += P_detected_IC[i] * math.exp(-opdat[j,k+1]) #need to add function that accounts for opacity due to EBL at particular redshift here, have to use freqtoeV(fq_mids_IC) to form this
            elif (freqtoeV(fq_mids_IC[i])>= opdat[-1,0]):
                 for k in range(len(F_gammaEBL)):
                    F_gammaEBL[k] += P_detected_IC[i] * math.exp(-opdat[-1,k+1])

#convert flux between 0.1 and 100GeV back to Luminosity in ERGS
L_gammaEBL = F_gammaEBL *((4.0*np.pi*d_Blazar**2.0)*(1.0+z)**2.0)
print(F_gammaEBL,L_gammaEBL)
'''
'''
with open('RepEBL>2<10Tev_theta0.txt', 'a') as f:
    f.write('\n%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf' % (L_gammaEBL[0], L_gammaEBL[1], L_gammaEBL[2], L_gammaEBL[3], L_gammaEBL[4], L_gammaEBL[5], L_gammaEBL[6], L_gammaEBL[7], L_gammaEBL[8], L_gammaEBL[9]))
'''
#___________________________________________________________________________________________#

'''
#THIS PART IS WITHOUT EBL______________________________________________________#
#Now also total radiant flux between 0.1 and 100GeV -> F_gamma:
F_gamma = 0
L_gamma = 0  #Luminostiy between 0.1 and 100GeV
for i in range(len(P_detected_IC)):
    if (freqtoeV(fq_mids_IC[i]) >= 2E12 and freqtoeV(fq_mids_IC[i]) <= 1E13):
        F_gamma += P_detected_IC[i]

#convert flux between 0.1 and 100GeV back to Luminosity in ERGS
L_gamma = F_gamma *((4.0*np.pi*d_Blazar**2.0)*(1.0+z)**2.0)
print(F_gamma,L_gamma)
'''
#with open('RepInt>2<10Tev_theta0.txt', 'a') as f:
    #f.write('\n%lf' % L_gamma)


'''
#__________________________________________________#
#WITHOUT EBL, PHOTONS s^-1
PhF_gamma = 0
PhL_gamma = 0  #Luminostiy between 0.1 and 100GeV
for i in range(len(P_detected_IC)):
    if (freqtoeV(fq_mids_IC[i]) >= 1E9 and freqtoeV(fq_mids_IC[i]) <= 1E11):
        PhF_gamma += (P_detected_IC[i] * 1.0E-7)/(6.62607004E-34*fq_mids_IC[i])

#convert photon flux to number of photons per second back to Luminosity IN ERGS
PhL_gamma = PhF_gamma *((4.0*np.pi*d_Blazar**2.0)*(1.0+z)**2.0)
print(PhF_gamma,PhL_gamma)

with open('O3Gamma_Photons100GeV_theta5.txt', 'a') as f:
    f.write('\n%lf' % PhL_gamma)

'''
EVPA_tot = 0
Pol_tot = 0
#____________________________________________#
'''
with open('EVPA+Pol.txt', 'a') as f:
    f.write('%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\t%.4lf\n' % (EVPA_tot,EVPA_tot_opt,EVPA_tot_rad,EVPA_tot_gam,Pol_tot,Pol_tot_opt,Pol_tot_rad,Pol_tot_gam))
print("SUCCESS")
'''
'''_____________________________________________'''


'''
    #now add some new code to plot the IC emission
ICdat = np.loadtxt('ICdat.txt')
IC_x = ICdat[:,0]*doppler_factor
IC_y = ICdat[:,1]*doppler_factor**4*1E7*(1/(4*np.pi*d_Blazar**2))#*5E4                           


#save data
#syncdat = np.c_[freqtoeV(fq_mids), P_detected]
#np.savetxt('synchrotron.txt', syncdat)


#add thin and thick contributions

thin = np.loadtxt('thin.txt')
#thin_flux = thin[:,0]*doppler_factor**4*1E7*(1/(4*np.pi*d_Blazar**2))
thin_flux = thin*doppler_factor**4*1E7*(1/(4*np.pi*d_Blazar**2))

thick = np.loadtxt('thick.txt')
#thick_flux = thick[:,0]*doppler_factor**4*1E7*(1/(4*np.pi*d_Blazar**2))
thick_flux = thick*doppler_factor**4*1E7*(1/(4*np.pi*d_Blazar**2))
'''

'''-----Original plot of just SED--------
#np.savetxt('RepSSC.txt',(freqtoeV(fq_mids_IC),P_detected_IC))
#plot the total synchrotron emission + IC
plt.figure()
plt.plot(freqtoeV(fq_mids_IC), P_detected_IC, 'r-.', label='IC')#Inverse Compton
plt.plot(freqtoeV(fq_mids), P_detected, 'b-', label='synchrotron') #synchrotron
plt.plot(ph_energy[0:30], flux_BL[0:30], 'k.', label='data 2008-2009')
plt.plot(ph_energy_MK, flux_MK, 'gx', label='data MK')
#plt.plot([1E10, 1E10], [1E-14, 1E-10], 'k-', lw=2)
#plt.plot([3E11, 3E11], [1E-14, 1E-10], 'k-', lw=1)
plt.xscale('log')
plt.yscale('log')
plt.ylabel(r'$\nu F_{\nu}$ (erg s$^{-1}$ cm$^{-2}$)', size='14')
plt.xlabel('eV', size='14')
plt.ylim(1E-14, 1E-9)
plt.xlim(1E-6, 1E14)#plt.xlim(1E-6, 1E12)
#plt.legend()
plt.yticks(size='12')
plt.xticks(size='12')
#plt.savefig('SED.png')

plt.show()'''

'''-----Plot with SED and polarisations----'''
'''
for i,f in enumerate(fq_mids):
    if freqtoeV(f) > 1.5E2 and freqtoeV(f) < 2.5E2:
        P_detected[i] = (P_detected[i-1] + P_detected[i+1])/2
for i,f in enumerate(fq_mids_IC):
    if freqtoeV(f) > 2.5E11 and freqtoeV(f) < 3.5E11:
        P_detected_IC[i] = (P_detected_IC[i-1] + P_detected_IC[i+1])/2
'''

#np.savetxt('POL3.txt',Pol)
#np.savetxt('EVPA3.txt',EVPA)
'''
with open('EVPA3.txt', 'a') as f:
    for i in range(len(EVPA)):
        if i == len(EVPA)-1:
            f.write('%.4lf\n' % (EVPA[i]))
        else:
            f.write('%.4lf\t' % (EVPA[i]))

with open('POL3.txt', 'a') as f:
    for i in range(len(Pol)):
        if i == len(Pol)-1:
            f.write('%.4lf\n' % (Pol[i]))
        else:
            f.write('%.4lf\t' % (Pol[i]))

print("SUCCESS")
'''


#np.savetxt('PSYNC3.txt',P_detected)
#np.savetxt('PIC3.txt',P_detected_IC)
#np.savetxt('FREQSYNC3.txt',fq_mids)
#np.savetxt('FREQIC3.txt',fq_mids_IC)


fig = plt.figure(1)
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
#line1 = ax0.plot(freqtoeV(fq_mids), Pol_Init,'r',label='Initial Population')
#line2 = ax0.plot(freqtoeV(fq_mids), Pol_single, color='g',label='Single e')
ax0.text(1E7,0.3,'$\Pi_{optical} =$ %.4f' % Pol_tot_opt, fontsize=10)
ax0.text(1E7,0.1,'$\Pi_{radio} =$ %.4f' % Pol_tot_rad, fontsize=10)
ax0.text(1E7,0.45,'$\Pi_{x-ray} =$ %.4f' % Pol_tot_gam, fontsize=10)
#ax0.text(100,0.2,'$\Pi^{Initial}_{tot} =$ %.5f' % Pol_Init_tot, fontsize=10)

ax0.legend()
#subplot for the EVPA
ax05 = plt.subplot(gs[1], sharex = ax0)
line05 = ax05.plot(freqtoeV(fq_mids), np.array(EVPA)*180/np.pi, color='g',label='EVPA')
#ax05.set_yscale("log")
ax05.set_xlim([1E-6, 1E13])
#ax05.set_ylim([-1.58, 1.58])
ax05.set_ylim([-90, 90])
yticks = ax05.yaxis.get_major_ticks()
ax05.text(1E7,-0.65*180/np.pi,'$EVPA_{optical} =$ %.4f' % (EVPA_tot_opt*180/np.pi), fontsize=10)
ax05.text(1E7,-1.2*180/np.pi,'$EVPA_{radio} =$ %.4f' % (EVPA_tot_rad*180/np.pi), fontsize=10)
ax05.text(1E7,-0.1*180/np.pi,'$EVPA_{x-ray} =$ %.4f' % (EVPA_tot_gam*180/np.pi), fontsize=10)
ax05.set_ylabel(r'$\theta_{sky}[deg]$', size='13')
ax05.legend()
#the second subplot
# shared axis X
ax1 = plt.subplot(gs[2], sharex = ax0)
line3 = ax1.plot(freqtoeV(fq_mids_IC), P_detected_IC, 'r-.', label='IC')#Inverse Compton
line4 = ax1.plot(freqtoeV(fq_mids), P_detected, 'b-', label='synchrotron') #synchrotron
line5 = ax1.plot(ph_energy_MK501[0:30], flux_MK501[0:30], 'k.', label='data 2008-2009')
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

'''
fig = plt.figure(figsize=(4,2))
ax1 = fig.add_axes([0.1, 0.5, 0.8, 0.4])
#ax1.set_yticks([0.1, 1, 10, 100])
#ax1.set_xticks([1E9, 1E10, 1E11, 1E12, 1E13])
line0, = ax1.plot(freqtoeV(fq_mids), Pol, color='b',label='Pol')
line1 = ax1.plot(freqtoeV(fq_mids), Pol_Init, color='r',label='Pol_Init')
line2 = ax1.plot(freqtoeV(fq_mids), Pol_single, color='g',label='Pol_single')
ax1.set_xscale('log')
ax1.set_ylim([0,1.0])
ax1.set_xlim([1E-6,1E13])
ax1.text(1E5,0.5,'$\Pi_{tot} =$ %.5f' % Pol_tot, fontsize=10)
ax1.text(1E10,0.7,'$\Pi^{Initial}_{tot} =$ %.5f' % Pol_Init_tot, fontsize=15)

ax2 = fig.add_axes([0.1, 0.5, 0.4, 0.4])
ax2.plot(NG[0,:], NG[1,:], 'm-',label='Intrinsic (No EBL)')
ax2.plot(NG[0,:], NG[2,:], 'b--', label='EBL')
ax2.plot(NG[0,:], NG[3,:], 'r-', label='EBL+ALP')
ax2.set_yscale('log')
ax2.set_xscale('log')
ax2.set_xticklabels('')
ax2.set_ylim([1.0001E-4,999.999])
ax2.set_xlim([1.0001E-16,9.9999E-12])
ax2.set_ylabel('N($>F_{\gamma}$)$[deg^{-2}]$',fontsize=18)
ax2.text(5E-15,100,'$20GeV<E_{\gamma}<200GeV$', fontsize=15)
ax2.plot([2.5E-20, 2.5E-20], [1E-5, 100], 'c-.', lw=2, label='CTA South (5$\sigma$,50hr)')
ax2.plot([5.8E-20, 5.8E-20], [1E-5, 100], 'k-.', lw=2, label='CTA North (5$\sigma$,50hr)')
ax2.plot([5.8E-20, 5.8E-20], [1E-5, 100], '-.',color=(1.0,0.41,0.7), lw=2, label='HAWC (5$\sigma$,1yr)')
ax2.plot([5.8E-20, 5.8E-20], [1E-5, 100], '-.',color=(1,0.79,0), lw=2, label='MAGIC (5$\sigma$,50hr)')
ax2.legend(fontsize=10, loc=3)
'''