__author__ = 'ALP'
import math as math
import numpy as np
def cross(a, b):
    c = [a[1]*b[2] - a[2]*b[1],
         a[2]*b[0] - a[0]*b[2],
         a[0]*b[1] - a[1]*b[0]]

    return c

def dot(a, b):
    c = 0
    for i,item in enumerate(a):
        c += a[i]*b[i]

    return c

def Dop_Dep(B,v, Gamma, n = np.array([0,0,1])): #usually have observer looking down z axis for the moment
    B = B / np.sqrt(np.dot(B,B))
    n = n / np.sqrt(np.dot(n,n))
    v = np.sqrt(1 - 1/Gamma**2)* v / np.sqrt(np.dot(v,v))
    q = B + np.cross(n,np.cross(v,B)) - (Gamma/(1+Gamma))*(np.dot(B,v))*v
    e = np.cross(n,q) / np.sqrt(np.dot(q,q) - (np.dot(n,q))**2)
    e_stationary = np.cross(n,B)
    if (e_stationary[0]!=0) or (e_stationary[1]!=0) or (e_stationary[2]!=0):
        e_stationary = e_stationary / np.sqrt(np.dot(e_stationary,e_stationary))
    return e , e_stationary

def Perp_B_extender(B,v,Gamma,n = np.array([0,0,1])): #NOT the same effect as Dop_Dep above, only the same in special case
    B = B / np.sqrt(np.dot(B,B))
    n = n / np.sqrt(np.dot(n,n))
    v = v / np.sqrt(np.dot(v,v))
    B_perp = np.cross(v,np.cross(B,v))     #B × (A×B / |B|) / |B| perpendicular to v
    B_para = np.dot(B,v)*v       #A•B * B/|B|2
    B_perp = B_perp*Gamma
    B_lab = B_perp + B_para
    B_lab = B_lab / np.sqrt(np.dot(B_lab,B_lab))
    e_stationary = np.cross(n,B)
    e = np.cross(n,B_lab)
    return e,e_stationary

#now want to check if two e_stationarys are perp, are the corresponding es also perp -> not necessarily the case


#visualisation
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.cm as cm
from matplotlib.colors import Normalize

def depol_plot(B,v,Gamma,n = [0,0,1],plot = True):
    B = [np.array(item) for item in B] # B input is given as list of lists now]
    v = np.array(v)
    n = np.array(n)


    U = [[] for item in B]
    V = [[] for item in B]
    W = [[] for item in B]
    X = [[] for item in B]
    Y = [[] for item in B]
    #Z = [[]]

    for i,item in enumerate(U):
        U[i], V[i], W[i] = zip(*Dop_Dep(B[i],v,Gamma,n))

        if i<=20:
            X[i] = (i*20,i*20)
            Y[i] = (0,0)
        elif i<=30:
            X[i] = ((i-11)*0.2,(i-11)*0.2)
            Y[i] = (0.5,0.5)
        elif i<=40:
            X[i] = ((i-21)*0.2,(i-21)*0.2)
            Y[i] = (1.0,1.0)
        elif i<=50:
            X[i] = ((i-31)*0.2,(i-31)*0.2)
            Y[i] = (1.5,1.5)
        elif i<= 60:
            X[i] = ((i-41)*0.2,(i-41)*0.2)
            Y[i] = (2.0,2.0)
        else:
            X[i] = ((i-51)*0.2,(i-51)*0.2)
            Y[i] = (2.5,2.5)

    #add jet axis arrow
    Xv = (-50,-50)
    Yv = (0,0)
    Uv = (v[0],-v[0])
    Vv = (v[1],-v[1])

    #for double ended arrows:
    Ud = []; Vd = []; Xd = []; Yd=[]
    for i,ting in enumerate(U):
        Ud.extend([U[i],(-U[i][0],-U[i][1])])
        Vd.extend([V[i],(-V[i][0],-V[i][1])])
        Xd.extend([X[i],X[i]])
        Yd.extend([Y[i],Y[i]])
    U = Ud; V=Vd; X=Xd; Y = Yd;
    '''
    EVPAs = [np.array([U[i][0],V[i][0],W[i][0]]) for i,u in enumerate(U)]
    EVPAorigs = [np.array([U[i][1],V[i][1],W[i][1]]) for i,u in enumerate(U)]
    EVPAsky = [np.cross(n,np.cross(EVPA,n)) for EVPA in EVPAs]
    EVPAorigsky = [np.cross(n,np.cross(EVPA,n)) for EVPA in EVPAorigs]
    vsky = np.cross(n,np.cross(v,n))
    angles = [np.arccos(np.dot(vsky,EVPAsky[i])/(np.sqrt(np.dot(vsky,vsky))*np.sqrt(np.dot(EVPAsky[i],EVPAsky[i]))))*(180/np.pi) for i, ting in enumerate(EVPAsky)]
    angles_orig = [np.arccos(np.dot(vsky,EVPAorigsky[i])/(np.sqrt(np.dot(vsky,vsky))*np.sqrt(np.dot(EVPAorigsky[i],EVPAorigsky[i]))))*(180/np.pi) for i, ting in enumerate(EVPAorigsky)]
    '''


    #Now find mean and std of arrows:
    #print(np.mean(np.array(U)[:,0]))
    #print(np.mean(np.array(U)[:,1]))
    #print(np.mean(np.array(V)[:,0]))
    #print(np.mean(np.array(V)[:,1]))
    '''
    xe_mean = np.mean(np.array(U)[:,0])
    ye_mean = np.mean(np.array(V)[:,0])
    xestat_mean = np.mean(np.array(U)[:,1])
    yestat_mean = np.mean(np.array(V)[:,1])

    return (xe_mean,ye_mean),(xestat_mean,yestat_mean)
    '''
    '''
    colormap = cm.inferno

    fig = plt.figure(figsize=(6, 6)) #plot means with jet direction also
    ax = plt.subplot()##fig.add_subplot(figsize=(4,4))
    colors = np.array([0.7])
    ax.quiver(0, 0, xe_mean, ye_mean,color=colormap(colors),angles='xy', scale_units='xy', scale=1)
    colors = np.array([0.1])
    ax.quiver(0, 0, xestat_mean, yestat_mean,color=colormap(colors),angles='xy', scale_units='xy', scale=1)
    colors = np.array([0.4]) #gotta scale v to be roughly the same size as the two mean arrows just for comparison:
    ax.quiver(0, 0, (np.sqrt(xe_mean**2 + ye_mean**2)/np.sqrt(v[0]**2+v[1]**2))*v[0], (np.sqrt(xe_mean**2 + ye_mean**2)/np.sqrt(v[0]**2+v[1]**2))*v[1],color=colormap(colors),angles='xy', scale_units='xy', scale=1)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_xlim([-3*np.sqrt(xe_mean**2 + ye_mean**2), 3*np.sqrt(xe_mean**2 + ye_mean**2)])
    ax.set_ylim([-3*np.sqrt(xe_mean**2 + ye_mean**2), 3*np.sqrt(xe_mean**2 + ye_mean**2)])
    plt.show()
    '''
    if plot:
        colors = np.array([0.7,0.1]) #orange is after doppler dep, black is no doppler dep (B field stationary)
        vcolors = np.array([0.5,0.5])

        #norm = Normalize()
        #norm.autoscale(colors)
        # we need to normalize our colors array to match it colormap domain
        # which is [0, 1]

        colormap = cm.inferno

        fig = plt.figure()
        ax = fig.add_subplot(111)
        for i, item in enumerate(U):
            if not i:
                ax.quiver(X[i][0], Y[i][0],U[i][0],V[i][0],color=colormap(colors[0]),label='EVPA')
                ax.quiver(X[i][1], Y[i][1],U[i][1],V[i][1],color=colormap(colors[1]),label='EVPA_orig')
            else:
                ax.quiver(X[i], Y[i],U[i],V[i],color=colormap(colors))
        colormap = cm.Greens
        ax.quiver(Xv, Yv,Uv,Vv,color=colormap(vcolors),label='Jet Projection')
        ax.set_xlabel('$\phi$ [deg]')
        ax.set_ylabel('y')
        ax.set_xlim([-70, 360])
        ax.set_ylim([-20, 20])
        plt.text(0, 35, '$\Gamma =$ {}'.format(Gamma), fontsize=12)
        #plt.legend(fontsize=10, loc=4)
        plt.show()

#Some different sets of B - fields here to test out
import random
random_B_set = [[random.randint(-50,50),random.randint(-50,50),random.randint(-50,50)] for i in range(50000)]
random_B_set2 = [[1,0,0],[0,1,0],[1,1,0],[-1,1,0],[-1,-1,0],[1,1,1],[2,-1,1],[3,-4,-2],[-5,1,1],[0,0,1],[0,1,1],[1,0,1],[1,1,5],[5,1,5]]
theta_obs = 4.0 * np.pi/180
v = [np.sin(theta_obs),0,np.cos(theta_obs)]
c_helix = 1
R = 1
helical_set = [[c_helix*np.sin(theta_obs)-R*np.cos(theta_obs)*np.sin(i*np.pi/8),R*np.cos(i*np.pi/8),R*np.sin(theta_obs)*np.sin(i*np.pi/8)+c_helix*np.cos(theta_obs)] for i in range(17)]
#helical_setUP = [[c_helix*np.sin(theta_obs)-R*np.cos(theta_obs)*np.sin(i*np.pi/8),R*np.cos(i*np.pi/8),R*np.sin(theta_obs)*np.sin(i*np.pi/8)+c_helix*np.cos(theta_obs)] for i in range(17)]


def depol_mean(v,Gamma,n = [0,0,1],plot = False):
    U_mean = [[] for item in range(20)]
    V_mean = [[] for item in range(20)]
    for i in range(20):
        random_B_set = [[random.randint(-50,50),random.randint(-50,50),random.randint(-50,50)] for i in range(10000)]
        U_mean[i], V_mean[i] = zip(*depol_plot(random_B_set,v,Gamma,n))

    v = np.array(v) / np.sqrt(np.dot(v,v))
    ev_mean = np.mean(abs(np.array(U_mean)[:,0])) #mean of e component along v
    estatv_mean = np.mean(abs(np.array(U_mean)[:,1])) #mean of e_stat along v

    evperp_mean = np.mean(abs(np.array(V_mean)[:,0]))
    evstatperp_mean = np.mean(abs(np.array(V_mean)[:,1]))
    return ev_mean,estatv_mean,evperp_mean,evstatperp_mean

#plotting gamma vs angle of EVPA to v
def depol_graph(B,v,Gamma,n = [0,0,1]):
    B = np.array(B) # B input is given as list of lists now]
    v = np.array(v)
    n = np.array(n)

    U = [[] for item in Gamma]
    V = [[] for item in Gamma]
    W = [[] for item in Gamma]
    for i,item in enumerate(U):
        U[i], V[i], W[i] = zip(*Dop_Dep(B,v,Gamma[i],n))
    return U,V,W

Gamma = [i for i in range(1,202)]
theta_obs = 4.0 * np.pi/180
v = [np.sin(theta_obs),0,np.cos(theta_obs)]
ths = [0.1,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2*0.99,(np.pi-0.1),5*np.pi/8,3*np.pi/4,7*np.pi/8]
#B = [[np.cos(th),np.sin(th),1.5] for th in ths]
#B = [[random.randint(-50,50),random.randint(-50,50),random.randint(-50,50)] for i in range(10)]
c_helix = 1
R = 1
B = [[c_helix*np.sin(theta_obs)-R*np.cos(theta_obs)*np.sin(i*np.pi/9),R*np.cos(i*np.pi/9),R*np.sin(theta_obs)*np.sin(i*np.pi/9)+c_helix*np.cos(theta_obs)] for i in range(18)]
fig = plt.figure(1)
plt.xlabel('$\Gamma$')
plt.ylabel('EVPA [deg] (cyan $\phi = 0$, blue $\phi = 120$)')
plt.title('$\Theta_{obs} = 4.0\circ$, B helical')
axes = plt.gca()
axes.set_xlim([-0.1, 100])
axes.set_ylim([-90.1, 90.1])
for i,b in enumerate(B):
    U,V,_ = depol_graph(b,v,Gamma,n = [0,0,1])
    vectors = [np.array([U[i][0],V[i][0]]) for i, ting in enumerate(U)]
    angles = [np.arctan(vector[1]/vector[0])*(180/np.pi) for vector in vectors]
    #angles = [np.arccos(np.dot(np.array([1,0]),vector)/(np.sqrt(np.dot(vector,vector))))*(180/np.pi) for vector in vectors]
    if not i:
        plt.plot(Gamma, angles,color='c',linewidth=1.2)
    elif i == 6:
        plt.plot(Gamma, angles,color='b',linewidth=1.2)
    else:
        plt.plot(Gamma, angles,color=(0.98,0.43,0),linewidth=1.2)
    plt.plot([-1, 1], [np.arctan(V[0][1]/U[0][0])*(180/np.pi) , np.arctan(V[0][1]/U[0][0])*(180/np.pi) ], color='k', linestyle='-', linewidth=2)
plt.plot([1/theta_obs,1/theta_obs],[-89.5,89.5],color='r',linestyle='-.',linewidth=2,label='$\Gamma = 1/\Theta_{obs}$')
#plt.plot([1,1],[-89.5,89.5],color='g',linestyle='--',linewidth=1,label='r$\Gamma = 1$')
plt.plot([4,4],[-89.5,89.5],color='b',linestyle='--',linewidth=1.3,label='$\Gamma = 4$')
plt.plot([10,10],[-89.5,89.5],color='g',linestyle='--',linewidth=1.3,label='$\Gamma = 10$')
plt.plot([90,90],[-89.5,89.5],color='m',linestyle='--',linewidth=1.3,label='$\Gamma = 90$')
plt.legend(loc=1)
    #can +90 to angles if want jet axis to be at theta=90


Gamma = [i for i in range(1,101)]
theta_obs = 4.0 * np.pi/180
v = [1,0,0]
ths = [0.1,np.pi/8,np.pi/4,3*np.pi/8,np.pi/2*0.99,(np.pi-0.1),5*np.pi/8,3*np.pi/4,7*np.pi/8]
#B = [[np.cos(th),np.sin(th),1.5] for th in ths]
#B = [[random.randint(-50,50),random.randint(-50,50),random.randint(-50,50)] for i in range(10)]
c_helix = 1
R = 1
#B = [[c_helix,R*np.cos(i*20 * np.pi/180),R*np.sin(i*20*np.pi/180)] for i in range(19)]
B = [[c_helix,R*np.cos(i * np.pi/8),R*np.sin(i*np.pi/8)] for i in range(17)]
phis = [i*180/8 for i in range(17)]
n = [np.cos(theta_obs),0,np.sin(theta_obs)]
fig = plt.figure(2)
plt.xlabel('$\Gamma$')
plt.ylabel('EVPA [deg]')
plt.title('$\Theta_{obs} = 4.0\circ$, B helical')
for j,b in enumerate(B):
    U,V,W = depol_graph(b,v,Gamma,n = [np.cos(theta_obs),0,np.sin(theta_obs)])
    EVPAs = [np.array([U[i][0],V[i][0],W[i][0]]) for i,u in enumerate(U)]
    EVPAorigs = [np.array([U[i][1],V[i][1],W[i][1]]) for i,u in enumerate(U)]
    EVPAsky = [np.cross(n,np.cross(EVPA,n)) for EVPA in EVPAs]
    EVPAorigsky = [np.cross(n,np.cross(EVPA,n)) for EVPA in EVPAorigs]
    vsky = np.cross(n,np.cross(v,n))
    angles = [np.arccos(np.dot(vsky,EVPAsky[i])/(np.sqrt(np.dot(vsky,vsky))*np.sqrt(np.dot(EVPAsky[i],EVPAsky[i]))))*(180/np.pi) for i, ting in enumerate(EVPAsky)]
    angles_orig = [np.arccos(np.dot(vsky,EVPAorigsky[0])/(np.sqrt(np.dot(vsky,vsky))*np.sqrt(np.dot(EVPAorigsky[0],EVPAorigsky[0]))))*(180/np.pi)]
    print(vsky)
    #print(angles[0])
    #vectors = [np.array([U[i][0],V[i][0]]) for i, ting in enumerate(U)]
    #angles = [-np.arctan(vector[1]/vector[0])*(180/np.pi) for vector in vectors]
    plt.plot(Gamma, angles,'b')
    plt.plot([1, 3], [angles_orig[0], angles_orig[0]], color='r', linestyle='-', linewidth=1)
    plt.plot([1/theta_obs,1/theta_obs],[0.5,179.5],color='g',linestyle='--',linewidth=1)

#need to evaluate angles in plane where n is the norm (this is the new plane of the sky)

jet_bias_obs = []
jet_bias_orig = []

for Gamma in range(2,100):
    random_B_set = [[random.randint(-50,50),random.randint(-50,50),random.randint(-50,50)] for i in range(50000)]
    print(Gamma)
    PA_obs = []
    PA_orig = []
    for b in random_B_set:
        e, e_orig = (Dop_Dep(np.array(b),np.array(v), Gamma, n = np.array([0,0,1])))
        if not math.isnan(e[0]) and not math.isnan(e[1]) and not math.isnan(e_orig[0]) and not math.isnan(e_orig[1]):
            PA_orig.append(e_orig)
            PA_obs.append(e)

    jet_bias_obs.append(np.std(np.array(PA_obs)[:,0]) - np.std(np.array(PA_obs)[:,1]))
    jet_bias_orig.append(np.std(np.array(PA_orig)[:,0]) - np.std(np.array(PA_orig)[:,1]))

plt.figure(1)
plt.plot(range(2,100),jet_bias_obs,color='orange',label='DPAR')
plt.plot(range(2,100),jet_bias_orig,color='k',label='Original')
plt.plot([1/theta_obs,1/theta_obs],[-0.005,0.055],color='r',linestyle='-.',linewidth=1.6,label='$\Gamma = 1/\Theta_{obs}$')
plt.ylabel('EVPA bias along jet axis')
plt.xlabel('$\Gamma$')
plt.title('$\Theta_{obs} = 4.0^{\circ}$')
plt.legend()


#print(sum(PA_obs)/50000)
#print(sum(PA_orig)/50000)
#print(np.std(np.array(PA_obs)[:,0]) - np.std(np.array(PA_obs)[:,1]))
#print(np.std(np.array(PA_orig)[:,0]) - np.std(np.array(PA_orig)[:,1]))
