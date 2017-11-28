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

def depol_plot(B,v,Gamma,n = [0,0,1],plot = False):
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

        if i<=10:
            X[i] = (i*0.2,i*0.2)
            Y[i] = (0,0)
        elif i<=20:
            X[i] = ((i-11)*0.2,(i-11)*0.2)
            Y[i] = (0.5,0.5)
        elif i<=30:
            X[i] = ((i-21)*0.2,(i-21)*0.2)
            Y[i] = (1.0,1.0)
        elif i<=40:
            X[i] = ((i-31)*0.2,(i-31)*0.2)
            Y[i] = (1.5,1.5)
        elif i<= 50:
            X[i] = ((i-41)*0.2,(i-41)*0.2)
            Y[i] = (2.0,2.0)
        else:
            X[i] = ((i-51)*0.2,(i-51)*0.2)
            Y[i] = (2.5,2.5)



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

    #return (xe_mean,ye_mean),(xestat_mean,yestat_mean)

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

        #norm = Normalize()
        #norm.autoscale(colors)
        # we need to normalize our colors array to match it colormap domain
        # which is [0, 1]

        colormap = cm.inferno

        fig = plt.figure()
        ax = fig.add_subplot(111)
        for i, item in enumerate(U):
            ax.quiver(X[i], Y[i], U[i], V[i],color=colormap(colors))
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_xlim([-0.1, 2.1])
        ax.set_ylim([-0.5, 3.5])
        plt.show()

#Some different sets of B - fields here to test out
import random
random_B_set = [[random.randint(-50,50),random.randint(-50,50),random.randint(-50,50)] for i in range(50000)]
random_B_set2 = [[1,0,0],[0,1,0],[1,1,0],[-1,1,0],[-1,-1,0],[1,1,1],[2,-1,1],[3,-4,-2],[-5,1,1],[0,0,1],[0,1,1],[1,0,1],[1,1,5],[5,1,5]]
theta_obs = 2.8 * np.pi/180
c_helix = 1
R = 1
helical_set = [[c_helix*np.sin(theta_obs)-R*np.cos(theta_obs)*np.sin(i*np.pi/8),R*np.cos(i*np.pi/8),R*np.sin(theta_obs)*np.sin(i*np.pi/8)+c_helix*np.cos(theta_obs)] for i in range(17)]

def depol_mean(v,Gamma,n = [0,0,1],plot = False):
    U_mean = [[] for item in range(20)]
    V_mean = [[] for item in range(20)]
    for i in range(20):
        random_B_set = [[random.randint(-50,50),random.randint(-50,50),random.randint(-50,50)] for i in range(10000)]
        U_mean[i], V_mean[i] = zip(*depol_plot(random_B_set,v,Gamma,n))

    v = np.array(v) / np.sqrt(np.dot(v,v))
    ev_mean = np.mean(abs(np.array(U_mean)[:,0]*v[0]+np.array(V_mean)[:,0]*v[1])) #mean of e component along v
    estatv_mean = np.mean(abs(np.array(U_mean)[:,1]*v[0]+np.array(V_mean)[:,1]*v[1])) #mean of e_stat along v

    evperp_mean = np.mean(abs(np.array(V_mean)[:,0]))
    evstatperp_mean = np.mean(abs(np.array(V_mean)[:,1]))
    return ev_mean,estatv_mean,evperp_mean,evstatperp_mean