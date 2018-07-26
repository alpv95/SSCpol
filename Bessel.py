__author__ = 'ALP'
from scipy import integrate
import numpy as np
import math as math
import matplotlib.pyplot as plt
'''Everything needed to calculate modified Bessel functions of the second kind and from those the synchrotron emission function F(x),G(x)'''

butts = lambda t: (math.gamma(13/6)*(2*5)**(13/6)/math.sqrt(math.pi))*math.cos(t) / ((t**2 + 5**2)**(13/6))
x = integrate.quad(butts, 0, np.inf)
print(x)

def BesselK(b,x):
    '''Calculates Bessel (K) function'''
    butts = lambda t: math.cos(t) / ((t**2 + x**2)**(b+1/2))
    int1 = integrate.quad(butts, 0, np.inf)
    return((math.gamma(b+1/2)*(2*x)**b/math.sqrt(math.pi)) * int1[0])

def F(x):
    #calculates synchrotron spectral function F(x), where x=f/f_c
    func = lambda z: x*BesselK(5/3,z)
    int = integrate.quad(func, x, 10)
    return int[0]

def G(x):
    #function used to find parallel and perpendicular pol power for electron emiiting synchrotron, Longair pg207
    return x*BesselK(2/3,x)

'''To get total polarisation for a single electron emitting synchrotron tot_pol = integral0->inf (F(x) + G(x) / F(x) - G(x) )'''
def pol1(x):
    func = lambda z: F(z)+G(z)
    int = integrate.quad(func, 0, 10)
    return int[0]

def pol2(x):
    func = lambda z: F(z)-G(z)
    int = integrate.quad(func, 0, 10)
    return int[0]

def pol3(x):
    func = lambda z: F(z)*z**(-1/2)
    int = integrate.quad(func, 0, 1000000)
    return int[0]

def pol4(x):
    func = lambda z: G(z)*z**(-1/2)
    int = integrate.quad(func, 0, 1000000)
    return int[0]



'''def INT3(f, a, b, N):
    x = np.linspace(a+(b-a)/(2*N), b-(b-a)/(2*N), N)
    fx=np.zeros(len(x))
    for i in range(len(x)):
        fx[i] = f(x[i])
    area = np.sum(fx)*(b-a)/N
    return area

def F(x):
    func = lambda z: x*BesselK(5/3,z)
    int = INT3(func, x, 5,1000)
    return int'''

x = np.logspace(-6,1,num=75)
x = x.tolist()
P_perp = []
P_para = []
for i,ting in enumerate(x):
    P_perp.append(F(ting)+G(ting))
    P_para.append(F(ting)-G(ting))

print(len(P_perp))
plt.figure()
plt.plot(x,P_perp,'b.')
plt.plot(x,P_para,'r.')
plt.show()

for i,ting in enumerate(P_perp):
    if (ting < 0):
        P_perp[i] = 0.0

xmid = []
for i,ting in enumerate(x):
    if (i<len(x)-1):
        xmid.append((x[i+1] - ting)/2 + ting)
    else:
        xmid.append(0.0)

thefile = open('FG.txt', 'w')
for i,ting in enumerate(xmid):
    y = "%.8f \t %.8f \t %.8f \n" % (ting, P_perp[i], P_para[i])
    print(y)
    thefile.write(y)
thefile.close()