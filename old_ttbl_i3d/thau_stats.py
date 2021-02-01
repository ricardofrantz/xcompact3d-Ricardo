import numpy as np
from numpy import zeros, float32, float64, concatenate, sqrt
import scipy as sp
import scipy.integrate
from scipy.integrate import simps
from scipy.special import erf
import matplotlib as mpl
import matplotlib.pyplot as plt
import csv

transp = False
formt = 'png'
qual = 1500

#width, height in inches
#single -> 90mm  -> 3.54
#1.5    -> 140mm -> 5.50
#double -> 190mm -> 7.48
siz1=3.54
siz2=2.54 #minimum 1.18

nx=216
ny=289
nz=216
ns=1

xlx=50.
yly=30.
zlz=30.

nclx=2
ncly=2
nclz=0

re=2000.
nu = 1./re
uinf = 1.

dt=0.016
imodulo=500
dsnapshot=imodulo*dt

x = zeros((nx), dtype=float64)
y = zeros((ny), dtype=float64)
z = zeros((nz), dtype=float64)

def simpson_nonuniform(f, x):
    """
    Simpson rule for irregularly spaced data.

        Parameters
        ----------
        x : list or np.array of floats
                Sampling points for the function values
        f : list or np.array of floats
                Function values at the sampling points

        Returns
        -------
        float : approximation for the integral
    """
    N = len(x) - 1
    h = np.diff(x)

    result = 0.0
    for i in range(1, N, 2):
        hph = h[i] + h[i - 1]
        result += f[i] * ( h[i]**3 + h[i - 1]**3
                           + 3. * h[i] * h[i - 1] * hph )\
                     / ( 6 * h[i] * h[i - 1] )
        result += f[i - 1] * ( 2. * h[i - 1]**3 - h[i]**3
                              + 3. * h[i] * h[i - 1]**2)\
                     / ( 6 * h[i - 1] * hph)
        result += f[i + 1] * ( 2. * h[i]**3 - h[i - 1]**3
                              + 3. * h[i - 1] * h[i]**2)\
                     / ( 6 * h[i] * hph )

    if (N + 1) % 2 == 0:
        result += f[N] * ( 2 * h[N - 1]**2
                          + 3. * h[N - 2] * h[N - 1])\
                     / ( 6 * ( h[N - 2] + h[N - 1] ) )
        result += f[N - 1] * ( h[N - 1]**2
                           + 3*h[N - 1]* h[N - 2] )\
                     / ( 6 * h[N - 2] )
        result -= f[N - 2] * h[N - 1]**3\
                     / ( 6 * h[N - 2] * ( h[N - 2] + h[N - 1] ) )
    return result

if nclx == 0:
    dx = xlx / nx
if nclx == 1 or nclx == 2:
    dx = xlx / (nx - 1.)
if ncly == 0:
    dy = yly / ny
if ncly == 1 or ncly == 2:
    dy = yly / (ny - 1.)
if nclz == 0:
    dz = zlz / nz
if nclz == 1 or nclz == 2:
    dz = zlz / (nz - 1.)

for i in range(nx):
    x[i] = float(i) * dx
for j in range(ny):
    y[j] = float(j) * dy
for k in range(nz):
    z[k] = float(k) * dz

class Theta(object):
    def __init__(self, filename):
        data = np.transpose(np.genfromtxt(filename))
        self.t  = data[0]
        self.td = data[1]
        self.r  = data[2]

f = Theta('../theta')
fig = plt.figure()
fig.set_size_inches(siz1, siz2)
plt.plot(f.t,f.td, lw=0.5, ls='-', color='k',label='thetad')
plt.plot(f.t,f.r*1000., lw=0.5, ls='--', color='r',label='res')
#plt.ylim(0.95, 1.05); plt.ylabel(r'$\theta$')
#plt.xlim(5e-1, 1300.); plt.xscale('log');
plt.xlabel(r'$t$')
plt.legend(loc='upper left',fontsize=6)
filename = 'ttbl_theta_dot'+ '.' + formt
print('Saving '+filename)
plt.savefig(filename,format=formt,dpi=qual,bbox_inches='tight',transparent=transp); plt.close(); plt.close('All')

# Integral quantities:
# Re_{\theta}   =       2000.703
# Re_{\delta^*} =       2827.938
# Re_{\tau}     =       671.1240
# H_{12}        =       1.413472
# c_f           =    0.003539148
# y/\delta_{99}       y+          U+          urms+       vrms+       wrms+       uv+         prms+       pu+         pv+         S(u)        F(u)        dU+/dy+     V+        omxrms^+    omyrms^+    omzrms^+ 

kth = np.fromfile('re_theta2000_kth', dtype=float, count=-1, sep=' ').reshape((513, 17)) 

#cf = 0.003539148
#re = 450
#utau = np.sqrt(cf/2.)
#lstar = re*utau
#upkth = np.asarray(kth[:,2])
ypk = np.asarray((kth[:,1]))
upk = np.asarray(kth[:,2])
#urms= np.asarray(kth[:,3])*utau
#vrms= np.asarray(kth[:,4])*utau
#wrms= np.asarray(kth[:,5])*utau

#fig = plt.figure()
#fig.set_size_inches(siz1, siz2)
#plt.scatter(up,yp,s=0.1)
#plt.scatter(erf(0.14*yp),yp,s=0.1)
#plt.scatter(erf(0.14*yp),yp,s=0.1)
#plt.scatter(erf(0.24*yp),yp,s=0.45)
#plt.scatter(urms,yp,s=0.1,color='r')
#plt.scatter(vrms,yp,s=0.1,color='g')
#plt.scatter(wrms,yp,s=0.1,color='b')
#plt.scatter(((1-erf(0.1*yp))*0.15),yp,s=1)
#plt.scatter((-0.5*(yp)**2),yp,s=1,color='k')
#plt.xlim(0,1)
#plt.ylim(0,30)
#plt.axhline(y=30)
#plt.savefig('umean'+ '.' + formt, format=formt, dpi=qual, bbox_inches='tight', transparent=transp); plt.close(); plt.close('All')

isk = 25000
itf = 67600
itf2 = 675

#t_uf = np.fromfile('0_t_ufric', dtype=float, count=-1, sep=' ').reshape((itf2, 2)) 
#fig = plt.figure()
#fig.set_size_inches(siz1, siz2)
#plt.scatter(2.*t_uf[:,0]**2.,t_uf[:,1],s=0.0005,color='b')
#plt.xlim(0,1)
#plt.ylim(0.04,0.05)
#plt.xlabel(r'$t$')
#plt.ylabel(r'$cf$')
#plt.axhline(y=0.003539148,lw=0.1)
#plt.savefig('ttbl_cf'+ '.' + formt, format=formt, dpi=qual, bbox_inches='tight', transparent=transp); plt.close(); plt.close('All')

yp = np.fromfile('../yp.dat', dtype=float, count=-1, sep=' ').reshape((ny)) 

#um = erf(0.2350*yp)
#print('erf 0.2350 theta =',simps(um * (1. - um), dx=dy))
#print('erf 0.2350 theta =',simps(um * (1. - um), yp))
#print('erf 0.2350 theta =',simpson_nonuniform(um * (1. - um), yp))

#um = erf(0.1*yp)
#print('erf 0.1 theta =',simps(um * (1. - um), yp))


#um = erf(0.15*yp)
#print('erf 0.15 theta =',simps(um * (1. - um), yp))

umean = np.fromfile('03_ux', dtype=float, count=-1, sep=' ').reshape((itf, ny))
#print('yp ',yp.shape)
#print('umean ',umean.shape)

fig = plt.figure()
fig.set_size_inches(siz1, siz2)
plt.scatter(kth[:,1],kth[:,2],s=0.2,color='k',label='KTH')

#delta = np.zeros(len(x))
#theta = np.zeros(len(x))
#for i in range(0, len(x)):
#    delta[i] = scipy.integrate.simps(1. - uxm[:, i], dx=dy, axis=0)
#    theta[i] = scipy.integrate.simps(uxm[:, i] * (1. - uxm[:, i]), dx=dy, axis=0)
        
col=int(itf*0.999)
um = umean[col,:]
ut = np.sqrt((um[1]/yp[1])*(1./re))
theta = simps(um * (1. - um), yp)
print(theta)
plt.scatter(yp*ut*re/theta,um/ut,s=0.2,color='b',label=r'TTBL $t=$'+str(round((itf+isk)*dt,2)))

um = 0.
iti=int(itf*0.5)
for i in range(iti,itf):
    um = um + umean[i,:]
um=um/(itf-iti)
print(um.shape)
ut = np.sqrt((um[1]/yp[1])*(1./re))
theta = simps(um * (1. - um), yp)
print(theta)
plt.scatter(yp*ut*re/theta,um/ut,s=0.2,color='r',label='mean')

yy = np.linspace(1.e-3, 2.9e3, 3.e5, endpoint=True); dyy = 2.9e3 / len(yy)
plt.plot(yy[600:int(200. / dyy)], (1. / 0.41) * np.log(yy[600:int(200. / dyy)]) + 5., linewidth=0.5,  linestyle='dashed', color='black', label=r'$log\ law$')
plt.plot(yy[0:775], yy[0:775], linewidth=0.5, linestyle='dotted', color='black', label=r'$U^{+}=y^{+}$')
plt.ylim(0., 25.); plt.ylabel(r'$u_1^{+}$')
plt.xlim(5e-1, 1300.); plt.xscale('log'); plt.ylabel(r'$u_1^{+}$')
plt.legend(loc='upper left',fontsize=6)
filename = 'ttbl_u'+ '.' + formt
print('Saving '+filename)
plt.savefig(filename,format=formt,dpi=qual,bbox_inches='tight',transparent=transp); plt.close(); plt.close('All')


#fig = plt.figure()
#fig.set_size_inches(siz1, siz2)
#umean = np.fromfile('15_epsilon', dtype=float, count=-1, sep=' ').reshape((itf, ny)) 
#plt.scatter(umean[-1,:],yp,s=0.0005,color='r')
#plt.scatter(umean[-2,:],yp,s=0.0005,color='g')
#plt.scatter(umean[-3,:],yp,s=0.0005,color='b')
#plt.scatter(umean[-4,:],yp,s=0.0005,color='m')
##plt.xlim(0,1)
#plt.ylim(0,5)
#plt.xlabel(r'$t$')
#plt.ylabel(r'$cf$')
#plt.savefig('ttbl_eps'+ '.' + formt, format=formt, dpi=qual, bbox_inches='tight', transparent=transp); plt.close(); plt.close('All')


'''

fig = plt.figure()
fig.set_size_inches(siz1, siz2)




umean = np.fromfile('07_uxp2', dtype=float, count=-1, sep=' ').reshape((itf, ny)) 
plt.scatter(yp*utau*re,umean[-1,:]/utau,s=0.5,color='r',label=r'TTBL $u_{1}^{\prime,2}$')

umean = np.fromfile('08_uyp2', dtype=float, count=-1, sep=' ').reshape((itf, ny)) 
plt.scatter(yp*utau*re,umean[-1,:]/utau,s=0.5,color='g',label=r'TTBL $u_{2}^{\prime,2}$')

umean = np.fromfile('09_uzp2', dtype=float, count=-1, sep=' ').reshape((itf, ny)) 
plt.scatter(yp*utau*re,umean[-1,:]/utau,s=0.5,color='b',label=r'TTBL $u_{3}^{\prime,2}$')


plt.ylim(0., 0.4)
plt.xlim(5e-1, 1300.)

plt.ylabel(r'$u^{\prime +}$')
plt.xlabel(r'$x_2^{+}$')
plt.xscale('log')
plt.legend(loc='best',fontsize=6)

plt.savefig('ttbl_uxp2'+ '.' + formt, format=formt, dpi=qual, bbox_inches='tight', transparent=transp); plt.close(); plt.close('All')


fig = plt.figure()
fig.set_size_inches(siz1, siz2)

umean = np.fromfile('10_uxp_uyp', dtype=float, count=-1, sep=' ').reshape((itf, ny)) 
plt.scatter(yp*utau*re,umean[-1,:]/utau,s=0.5,color='r',label=r'TTBL $u_{1}^{\prime}u_{2}^{\prime}$')

umean = np.fromfile('11_uyp_uzp', dtype=float, count=-1, sep=' ').reshape((itf, ny)) 
plt.scatter(yp*utau*re,umean[-1,:]/utau,s=0.5,color='g',label=r'TTBL $u_{2}^{\prime}u_{3}^{\prime}$')

umean = np.fromfile('12_uxp_uzp', dtype=float, count=-1, sep=' ').reshape((itf, ny)) 
plt.scatter(yp*utau*re,umean[-1,:]/utau,s=0.5,color='b',label=r'TTBL $u_{1}^{\prime}u_{3}^{\prime}$')

#plt.ylim(0., 0.05)
plt.xlim(5e-1, 1300.)

plt.ylabel(r'$u_1^{\prime +}u_2^{\prime +}$')
plt.xlabel(r'$x_2^{+}$')
plt.xscale('log')
plt.legend(loc='best',fontsize=6)
plt.savefig('ttbl_uxp_uyp'+ '.' + formt, format=formt, dpi=qual, bbox_inches='tight', transparent=transp); plt.close(); plt.close('All')



'''
class Theta(object):
    def __init__(self, filename):
        data = np.transpose(np.genfromtxt(filename))
        self.t  = data[0]
        self.th = data[1]
        self.r  = data[2]

f = Theta('../theta')
fig = plt.figure()
fig.set_size_inches(siz1, siz2)
plt.plot(f.t,f.th, lw=0.5, ls='-', color='black')
plt.ylim(0.95, 1.05); plt.ylabel(r'$\theta$')
#plt.xlim(5e-1, 1300.); plt.xscale('log');
plt.xlabel(r'$t$')
#plt.legend(loc='upper left',fontsize=6)
filename = 'ttbl_theta'+ '.' + formt
print('Saving '+filename)
plt.savefig(filename,format=formt,dpi=qual,bbox_inches='tight',transparent=transp); plt.close(); plt.close('All')

fig = plt.figure()
fig.set_size_inches(siz1, siz2)
plt.plot(f.t,f.r, lw=0.5, ls='-', color='black')
#plt.ylim(0., 25.); plt.ylabel(r'$u_1^{+}$')
#plt.xlim(5e-1, 1300.); plt.xscale('log'); plt.ylabel(r'$u_1^{+}$')
#plt.legend(loc='upper left',fontsize=6)
filename = 'ttbl_theta_res'+ '.' + formt
print('Saving '+filename)
plt.savefig(filename,format=formt,dpi=qual,bbox_inches='tight',transparent=transp); plt.close(); plt.close('All')


f = Theta('../thetad')
fig = plt.figure()
fig.set_size_inches(siz1, siz2)
plt.plot(f.t,f.th, lw=0.5, ls='-', color='black')
plt.plot(f.t,f.r, lw=0.5, ls='-', color='black')

#plt.ylim(0., 25.); plt.ylabel(r'$u_1^{+}$')
#plt.xlim(5e-1, 1300.); plt.xscale('log'); plt.ylabel(r'$u_1^{+}$')
#plt.legend(loc='upper left',fontsize=6)
filename = 'ttbl_thetd'+ '.' + formt
print('Saving '+filename)
plt.savefig(filename,format=formt,dpi=qual,bbox_inches='tight',transparent=transp); plt.close(); plt.close('All')


