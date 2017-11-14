__author__ = 'ALP'
import math as math
import numpy as np
import matplotlib.pyplot as plt
# used to plot the polarisation and EVPA
Pol = np .loadtxt('EVPA+Pol.txt')
t = [i for i in range(len(Pol))]
EVPAband = 1
Polband = 5
'''
plt.figure(11)
plt.plot(t[185:],Pol[185:,0],'--')
plt.plot(t[155:185],Pol[155:185,0],'--')
plt.plot(t[125:155],Pol[125:155,0],'-.')
plt.plot(t[95:125],Pol[95:125,0],'-.')
plt.plot(t[:30],Pol[:30,0],'-')
plt.plot(t[62:92],Pol[62:92,0],'--')
plt.plot(t[:],Pol[:,1],'-')
plt.xlabel('Time [d]')
plt.ylabel('EVPA [rad], $\Pi$')


plt.text(3,1.4,'$\Pi_{mean} =$ %.5f' % np.mean(Pol[:30,1]), fontsize=10)
plt.text(3,1.25,'$\sigma_{\Pi} =$ %.5f' % np.std(Pol[:30,1]), fontsize=10)
plt.text(3,1.1,'$EVPA_{mean} =$ %.5f' % np.mean(Pol[:30,0]), fontsize=10)
plt.text(3,0.95,'$\sigma_{EVPA} =$ %.5f' % np.std(Pol[:30,0]), fontsize=10)

plt.text(59,1.4,'$\Pi_{mean} =$ %.5f' % np.mean(Pol[62:92,1]), fontsize=10)
plt.text(59,1.25,'$\sigma_{\Pi} =$ %.5f' % np.std(Pol[62:92,1]), fontsize=10)
plt.text(59,1.1,'$EVPA_{mean} =$ %.5f' % np.mean(Pol[62:92,0]), fontsize=10)
plt.text(59,0.95,'$\sigma_{EVPA} =$ %.5f' % np.std(Pol[62:92,0]), fontsize=10)

plt.text(92,1.4,'$\Pi_{mean} =$ %.5f' % np.mean(Pol[95:125,1]), fontsize=10)
plt.text(92,1.25,'$\sigma_{\Pi} =$ %.5f' % np.std(Pol[95:125,1]), fontsize=10)
plt.text(92,1.1,'$EVPA_{mean} =$ %.5f' % np.mean(Pol[95:125,0]), fontsize=10)
plt.text(92,0.95,'$\sigma_{EVPA} =$ %.5f' % np.std(Pol[95:125,0]), fontsize=10)

plt.text(125,1.4,'$\Pi_{mean} =$ %.5f' % np.mean(Pol[125:155,1]), fontsize=10)
plt.text(125,1.25,'$\sigma_{\Pi} =$ %.5f' % np.std(Pol[125:155,1]), fontsize=10)
plt.text(125,1.1,'$EVPA_{mean} =$ %.5f' % np.mean(Pol[125:155,0]), fontsize=10)
plt.text(125,0.95,'$\sigma_{EVPA} =$ %.5f' % np.std(Pol[125:155,0]), fontsize=10)

plt.text(155,1.4,'$\Pi_{mean} =$ %.5f' % np.mean(Pol[155:185,1]), fontsize=10)
plt.text(155,1.25,'$\sigma_{\Pi} =$ %.5f' % np.std(Pol[155:185,1]), fontsize=10)
plt.text(155,1.1,'$EVPA_{mean} =$ %.5f' % np.mean(Pol[155:185,0]), fontsize=10)
plt.text(155,0.95,'$\sigma_{EVPA} =$ %.5f' % np.std(Pol[155:185,0]), fontsize=10)

plt.text(185,1.4,'$\Pi_{mean} =$ %.5f' % np.mean(Pol[185:,1]), fontsize=10)
plt.text(185,1.25,'$\sigma_{\Pi} =$ %.5f' % np.std(Pol[185:,1]), fontsize=10)
plt.text(185,1.1,'$EVPA_{mean} =$ %.5f' % np.mean(Pol[185:,0]), fontsize=10)
plt.text(185,0.95,'$\sigma_{EVPA} =$ %.5f' % np.std(Pol[185:,0]), fontsize=10)
'''
plt.figure(12)
plt.plot(t[2:32],Pol[2:32,Polband],'--')
plt.plot(t[2:32],Pol[2:32,EVPAband],'-')
plt.plot(t[32:62],Pol[32:62,Polband],'--')
plt.plot(t[32:62],Pol[32:62,EVPAband],'-')
plt.plot(t[62:92],Pol[62:92,Polband],'--')
plt.plot(t[62:92],Pol[62:92,EVPAband],'-')
plt.plot(t[92:122],Pol[92:122,Polband],'--')
plt.plot(t[92:122],Pol[92:122,EVPAband],'-')
plt.plot(t[124:154],Pol[124:154,Polband],'--')
plt.plot(t[124:154],Pol[124:154,EVPAband],'-')
plt.plot(t[154:],Pol[154:,Polband],'--')
plt.plot(t[154:],Pol[154:,EVPAband],'-')
plt.xlabel('Time [d]')
plt.ylabel('EVPA [rad], $\Pi$')

plt.text(2,1.4,'$\Pi_{mean} =$ %.5f' % np.mean(Pol[2:32,Polband]), fontsize=10)
plt.text(2,1.25,'$\sigma_{\Pi} =$ %.5f' % np.std(Pol[2:32,Polband]), fontsize=10)
plt.text(2,1.1,'$EVPA_{mean} =$ %.5f' % np.mean(Pol[2:32,EVPAband]), fontsize=10)
plt.text(2,0.95,'$\sigma_{EVPA} =$ %.5f' % np.std(Pol[2:32,EVPAband]), fontsize=10)

plt.text(32,1.4,'$\Pi_{mean} =$ %.5f' % np.mean(Pol[32:62,Polband]), fontsize=10)
plt.text(32,1.25,'$\sigma_{\Pi} =$ %.5f' % np.std(Pol[32:62,Polband]), fontsize=10)
plt.text(32,1.1,'$EVPA_{mean} =$ %.5f' % np.mean(Pol[32:62,EVPAband]), fontsize=10)
plt.text(32,0.95,'$\sigma_{EVPA} =$ %.5f' % np.std(Pol[32:62,EVPAband]), fontsize=10)

plt.text(62,1.4,'$\Pi_{mean} =$ %.5f' % np.mean(Pol[62:92,Polband]), fontsize=10)
plt.text(62,1.25,'$\sigma_{\Pi} =$ %.5f' % np.std(Pol[62:92,Polband]), fontsize=10)
plt.text(62,1.1,'$EVPA_{mean} =$ %.5f' % np.mean(Pol[62:92,EVPAband]), fontsize=10)
plt.text(62,0.95,'$\sigma_{EVPA} =$ %.5f' % np.std(Pol[62:92,EVPAband]), fontsize=10)

plt.text(92,1.4,'$\Pi_{mean} =$ %.5f' % np.mean(Pol[92:122,Polband]), fontsize=10)
plt.text(92,1.25,'$\sigma_{\Pi} =$ %.5f' % np.std(Pol[92:122,Polband]), fontsize=10)
plt.text(92,1.1,'$EVPA_{mean} =$ %.5f' % np.mean(Pol[92:122,EVPAband]), fontsize=10)
plt.text(92,0.95,'$\sigma_{EVPA} =$ %.5f' % np.std(Pol[92:122,EVPAband]), fontsize=10)

plt.text(122,1.4,'$\Pi_{mean} =$ %.5f' % np.mean(Pol[124:156,Polband]), fontsize=10)
plt.text(122,1.25,'$\sigma_{\Pi} =$ %.5f' % np.std(Pol[124:156,Polband]), fontsize=10)
plt.text(122,1.1,'$EVPA_{mean} =$ %.5f' % np.mean(Pol[124:156,EVPAband]), fontsize=10)
plt.text(122,0.95,'$\sigma_{EVPA} =$ %.5f' % np.std(Pol[124:156,EVPAband]), fontsize=10)


plt.text(152,1.4,'$\Pi_{mean} =$ %.5f' % np.mean(Pol[156:186,Polband]), fontsize=10)
plt.text(152,1.25,'$\sigma_{\Pi} =$ %.5f' % np.std(Pol[156:186,Polband]), fontsize=10)
plt.text(152,1.1,'$EVPA_{mean} =$ %.5f' % np.mean(Pol[156:186,EVPAband]), fontsize=10)
plt.text(152,0.95,'$\sigma_{EVPA} =$ %.5f' % np.std(Pol[156:186,EVPAband]), fontsize=10)


plt.show()