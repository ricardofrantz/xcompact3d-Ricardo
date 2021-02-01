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
qual = 500

#width, height in inches
#single -> 90mm  -> 3.54
#1.5    -> 140mm -> 5.50
#double -> 190mm -> 7.48
siz1=3.54
siz2=2.54 #minimum 1.18

nx=128
ny=257
nz=128
ns=1

xlx=50.
yly=30.
zlz=30.

nclx=2
ncly=2
nclz=0

re=2000.

dt=0.0010
imodulo=500
dsnapshot=imodulo*dt

x = zeros((nx), dtype=float64)
y = zeros((ny), dtype=float64)
z = zeros((nz), dtype=float64)

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

f = Theta('theta')
fig = plt.figure()
fig.set_size_inches(siz1, siz2)
plt.scatter(f.t,f.td, s=0.05, color='k',label=r'$\dot{\theta}$')
#plt.plot(f.t,f.r*1000., lw=0.5, ls='--', color='r',label='res')
#plt.ylim(-0.01, 0.01); plt.ylabel(r'$\theta$')
plt.xlim(0., 500.)#; plt.xscale('log');
plt.xlabel(r'$t$')
plt.legend(loc='upper left',fontsize=6)
filename = 'ttbl_theta_dot'+ '.' + formt
print('Saving '+filename)
plt.savefig(filename,format=formt,dpi=qual,bbox_inches='tight',transparent=transp); plt.close(); plt.close('All')

fig = plt.figure()
fig.set_size_inches(siz1, siz2)
plt.scatter(f.t,1.0+f.r, s=0.051, color='k',label=r'$\theta$')
#plt.ylim(-0.01, 0.01); plt.ylabel(r'$\theta$')
plt.xlim(0., 500.)#; plt.xscale('log');
plt.xlabel(r'$t$')
plt.legend(loc='upper left',fontsize=6)
filename = 'ttbl_theta_res'+ '.' + formt
print('Saving '+filename)
plt.savefig(filename,format=formt,dpi=qual,bbox_inches='tight',transparent=transp); plt.close(); plt.close('All')




'''

itf = 0
yp = np.fromfile('yp.dat', dtype=float, count=-1, sep=' ').reshape((ny)) 
ux2m = np.fromfile('ux2m', dtype=float, count=-1, sep=' ').reshape((itf, ny))
ux2mt1 = np.fromfile('ux2m_theta1', dtype=float, count=-1, sep=' ').reshape((itf, ny))
ux2mt2 = np.fromfile('ux2m_theta2', dtype=float, count=-1, sep=' ').reshape((itf, ny))
du2dy22m = np.fromfile('du2dy22m', dtype=float, count=-1, sep=' ').reshape((itf, ny))
dudy2m = np.fromfile('dudy2m', dtype=float, count=-1, sep=' ').reshape((itf, ny))
duxuy2pm = np.fromfile('duxuy2pm', dtype=float, count=-1, sep=' ').reshape((itf, ny))


fig = plt.figure()
fig.set_size_inches(siz1, siz2)
plt.plot(ux2m[-1,:],yp,color='k',label='ux2m',lw=3)
plt.plot(ux2mt1[-1,:],yp,color='b',label='ux2mt1',lw=1.5)
plt.plot(ux2mt2[-1,:],yp,color='g',label='ux2mt2',lw=0.6)
#plt.ylim(0., 30.)#; plt.ylabel(r'$u_1^{+}$')
#plt.xlim(5e-1, 1300.); plt.xscale('log'); plt.ylabel(r'$u_1^{+}$')
plt.legend(loc='upper left',fontsize=6)
filename = 'ux2m'+ '.' + formt
print('Saving '+filename)
plt.savefig(filename,format=formt,dpi=qual,bbox_inches='tight',transparent=transp); plt.close(); plt.close('All')

fig = plt.figure()
fig.set_size_inches(siz1, siz2)
#plt.plot(ux2m[-1,:],yp,color='k',label='ux2m')
plt.plot(dudy2m[-1,:],yp,color='r',label='dudy2m')
plt.plot(du2dy22m[-1,:],yp,color='g',label='du2dy22m')
plt.ylim(0., 10.)#; plt.ylabel(r'$u_1^{+}$')
#plt.xlim(5e-1, 1300.); plt.xscale('log'); plt.ylabel(r'$u_1^{+}$')
plt.legend(loc='upper left',fontsize=6)
filename = 'du2'+ '.' + formt
print('Saving '+filename)
plt.savefig(filename,format=formt,dpi=qual,bbox_inches='tight',transparent=transp); plt.close(); plt.close('All')

fig = plt.figure()
fig.set_size_inches(siz1, siz2)
plt.plot(duxuy2pm[-1,:],yp,color='b',label='duxuy2pm')
plt.ylim(0., 10.)#; plt.ylabel(r'$u_1^{+}$')
#plt.xlim(5e-1, 1300.); plt.xscale('log'); plt.ylabel(r'$u_1^{+}$')
plt.legend(loc='upper left',fontsize=6)
filename = 'duxuy2pm'+ '.' + formt
print('Saving '+filename)
plt.savefig(filename,format=formt,dpi=qual,bbox_inches='tight',transparent=transp); plt.close(); plt.close('All')
'''
