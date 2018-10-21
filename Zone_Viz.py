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




# colors = np.array([0.7,0.1]) #orange is after doppler dep, black is no doppler dep (B field stationary)
# vcolors = np.array([0.5,0.5])
#
#
#
#
# colormap = cm.inferno
# fig = plt.figure(1)
# ax = fig.add_subplot(111)
# for i, item in enumerate(U):
#     if not i:
#         ax.quiver(X[i][0], Y[i][0],U[i][0],V[i][0],color=colormap(colors[0]),label='EVPA')
#         ax.quiver(X[i][1], Y[i][1],U[i][1],V[i][1],color=colormap(colors[1]),label='EVPA_orig')
#     else:
#         ax.quiver(X[i], Y[i],U[i],V[i],color=colormap(colors))
# colormap = cm.Greens
# ax.quiver(Xv, Yv,Uv,Vv,color=colormap(vcolors),label='Jet Projection')
# ax.set_xlabel('$\phi$ [deg]')
# ax.set_ylabel('y')
# ax.set_xlim([-70, 360])
# ax.set_ylim([-20, 20])
# plt.text(0, 35, '$\Gamma =$ {}'.format(Gamma), fontsize=12)