'''
Visualises jet zones with their B-fields given Proj_Bfile
'''

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

def PlotAndSave(count):

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




