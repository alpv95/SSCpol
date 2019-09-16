import math as math
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as sc
x = np.random.uniform(-5,5,15)
y = np.random.uniform(-5,5,15)
D = np.array([x*x,x*y,y*y,x,y,np.ones(len(x))]).T
S = np.dot(D.T,D)
C = np.zeros((6,6))
C[0,2] = 2
C[1,1] = -1
C[2,0] = 2
w,v = np.linalg.eig(np.dot(np.linalg.inv(S),C))
g = None
for i,ting in enumerate(w):
    if ting > 0:
        g = ting
        u = v[:,i]
        break
#print(w)
#print(g)
if g is not None:
    mu = np.sqrt(g/(np.dot(np.expand_dims(u,0),np.dot(S,np.expand_dims(u,1)))))
    #print(mu)
    a = mu*u
    a = np.squeeze(a)
    print(a)
    def ellipse(x):
        return a[0]*x[0]**2 + a[1]*x[0]*x[1] + a[2]*x[1]**2 + a[3]*x[0] + a[4]*x[1] + a[5]

    X = []
    Y1 = []
    Y2 =[]

    for i in range(1000):
        XX = 15/500 * (i-500)
        try:
            YY1 = -(np.sqrt((a[1]*XX+a[4])**2 - 4*a[2]*(XX*(a[0]*XX+a[3])+a[5]))+a[1]*XX+a[4])/(2*a[2])
            YY2 = -(-np.sqrt((a[1]*XX+a[4])**2 - 4*a[2]*(XX*(a[0]*XX+a[3])+a[5]))+a[1]*XX+a[4])/(2*a[2])
        except:
            continue
        if not math.isnan(YY1) and not math.isnan(YY2):
            X.append(XX)
            Y1.append(YY1)
            Y2.append(YY2)

    eta = np.linalg.det(np.array([[a[0],a[1]/2,a[3]/2],[a[1]/2,a[2],a[4]/2],[a[3]/2,a[4]/2,a[5]]]))
    if eta > 0:
        E = np.sqrt((2*np.sqrt((a[0]-a[2])**2 + a[1]**2))/(-(a[0]+a[2]) + np.sqrt((a[0]-a[2])**2 + a[1]**2)))
    elif eta < 0:
        E = np.sqrt((2*np.sqrt((a[0]-a[2])**2 + a[1]**2))/((a[0]+a[2]) + np.sqrt((a[0]-a[2])**2 + a[1]**2)))

    plt.figure(1)
    plt.plot(x,y,'.')
    plt.plot(X,Y1,'r')
    plt.plot(X,Y2,'r')
    plt.text(0,0,'Eccentricity = {:.3f}'.format(E))


