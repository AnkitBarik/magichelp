#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from magic import *
from koreviz import *

radratio = 0.35
ri = radratio/(1-radratio)
ro = 1/(1-radratio)

nr     = 100
ntheta = 256
nphi   = 512

eps = 1e-8

r     = np.linspace(ri,ro,nr)
theta = np.linspace(eps,np.pi-eps,ntheta)
phi   = np.linspace(0,2*np.pi,nphi)

def get_grid(r, theta, phi, nr, ntheta, nphi):

    r3D  = np.zeros([nphi, ntheta, nr])
    th3D = np.zeros([nphi, ntheta, nr])
    p3D  = np.zeros([nphi, ntheta, nr])

    for i in range(nr):
        r3D[..., i] = r[i]
    for j in range(ntheta):
        th3D[:, j, :] = theta[j]
    for k in range(nphi):
        p3D[k, ...] = phi[k]

    return r3D, th3D, p3D

r3D, th3D, p3D = get_grid(r, theta, phi, nr, ntheta, nphi)



# Test1
########
# Bz = 1

# Br     = Bz * np.cos(th3D)
# Btheta = -Bz * np.sin(th3D)
# Bphi   = np.zeros_like(r3D)

# vr     = np.zeros_like(r3D)
# vtheta = np.zeros_like(r3D)
# vphi   = -0.25*(np.sqrt(7/np.pi) * (3*np.sin(th3D) - 15*np.cos(th3D)**2 * np.sin(th3D)))/r3D
# bgradu_p_analytical = (3*Bz*np.sqrt(7/np.pi)*(2*np.sin(2*th3D) - 5*np.sin(4*th3D)))/(8.*r3D**2)

# Test 2
########
# Bz = 1

# Br     = Bz * np.cos(th3D)
# Btheta = -Bz * np.sin(th3D)
# Bphi   = np.zeros_like(r3D)

# vr     = np.zeros_like(r3D)
# vtheta = -((np.sqrt(105/(2.*np.pi))*np.cos(p3D)*np.cos(th3D)*np.sin(p3D)*np.sin(th3D))/r3D)
# vphi   = (np.sqrt(105/(2.*np.pi))*np.cos(2*p3D)*(np.sin(th3D) - 3*np.sin(3*th3D)))/(16.*r3D)

# bgradu_r_a = -((Bz*np.sqrt(105/(2.*np.pi))*np.cos(p3D)*np.cos(th3D)*np.sin(p3D)*np.sin(th3D)**2)/r3D**2)
# bgradu_t_a = (Bz*np.sqrt(105/(2.*np.pi))*(1 + 3*np.cos(2*th3D))*np.sin(2*p3D)*np.sin(th3D))/(4.*r3D**2)
# bgradu_p_a = (Bz*np.sqrt(105/(2.*np.pi))*np.cos(2*p3D)*(-2*np.sin(2*th3D) + 3*np.sin(4*th3D)))/(8.*r3D**2)

# Test 3
########

vr     = (-3*np.sqrt(21/np.pi)*np.cos(p3D)*(3 + 5*np.cos(2*th3D))*np.sin(th3D))/(4.*r3D**2)
vtheta = -((np.sqrt(105/(2.*np.pi))*np.cos(p3D)*np.cos(th3D)*np.sin(p3D)*np.sin(th3D))/r3D)
vphi   = (np.sqrt(105/(2.*np.pi))*np.cos(2*p3D)*(np.sin(th3D) - 3*np.sin(3*th3D)))/(16.*r3D)

Br     = (15*np.sqrt(5/(2.*np.pi))*np.cos(2*p3D)*(5 + 7*np.cos(2*th3D))*np.sin(th3D)**2)/(4.*r3D**2)
Btheta = (3*np.sqrt(5/np.pi)*np.cos(th3D)*(1 + 7*np.cos(2*th3D))*np.sin(p3D))/(16.*r3D)
Bphi   = (3*np.sqrt(5/np.pi)*np.cos(p3D)*(np.cos(2*th3D) + 7*np.cos(4*th3D)))/(16.*r3D)

bgradu_r_a = (-72*np.sqrt(105)*r3D*(116*np.cos(2*th3D) + 35*(3 + np.cos(4*th3D)))*np.sin(2*p3D)*np.sin(th3D)**2 + 1440*np.sqrt(210)*(np.cos(p3D) + np.cos(3*p3D))*(65 + 92*np.cos(2*th3D) + 35*np.cos(4*th3D))*np.sin(th3D)**3 -
          15*np.sqrt(42)*r3D**2*(np.cos(p3D)*(13*np.sin(th3D) - 24*np.sin(3*th3D) + 4*(-7 + 9*np.cos(2*p3D))*np.sin(5*th3D) - 21*np.sin(7*th3D)) + 3*np.cos(3*p3D)*(7*np.sin(th3D) + 4*np.sin(3*th3D) - 7*np.sin(7*th3D))))/(2048.*np.pi*r3D**5)

bgradu_t_a = (3*np.sqrt(21)*(1600*np.cos(th3D)*(5 + 7*np.cos(2*th3D))*np.sin(4*p3D)*np.sin(th3D)**3 - 10*np.sqrt(2)*r3D*(np.cos(p3D)*np.cos(th3D)*(25 + np.cos(2*th3D) + 123*np.cos(4*th3D) - 21*np.cos(6*th3D)) + 2*np.cos(3*p3D)*(-81*np.cos(th3D) - 4*np.cos(3*th3D) + 21*np.cos(5*th3D))*np.sin(th3D)**2) -
            6*np.sqrt(5)*np.sin(2*p3D)*(47*np.sin(2*th3D) + 52*np.sin(4*th3D) + 35*np.sin(6*th3D))))/(2048.*np.pi*r3D**4)

bgradu_p_a = (3*np.sqrt(21)*(5*np.sqrt(2)*r3D*(70 + 177*np.cos(2*th3D) + 90*np.cos(4*th3D) + 175*np.cos(6*th3D))*np.sin(p3D) + 62*(25 + 12*np.sqrt(5))*np.sin(th3D) - 20*np.sqrt(2)*r3D*(39 + 40*np.cos(2*th3D) + 49*np.cos(4*th3D))*np.sin(3*p3D)*np.sin(th3D)**2 +
           800*np.cos(4*p3D)*(1 + 3*np.cos(2*th3D))*(5 + 7*np.cos(2*th3D))*np.sin(th3D)**3 + 2*(175 + 72*np.sqrt(5))*np.sin(3*th3D) + (950 - 288*np.sqrt(5))*np.sin(5*th3D) + 24*np.sqrt(5)*np.cos(2*p3D)*(31*np.sin(th3D) + 6*np.sin(3*th3D) - 12*np.sin(5*th3D) - 35*np.sin(7*th3D)) -
           210*(5 + 4*np.sqrt(5))*np.sin(7*th3D)))/(4096.*np.pi*r3D**4)

ugradb_r_a = (15*np.sqrt(21)*np.sin(th3D)*(48*np.sqrt(10)*(np.cos(p3D) + np.cos(3*p3D))*(65 + 92*np.cos(2*th3D) + 35*np.cos(4*th3D))*np.sin(th3D)**2 + 160*r3D*(9 + 7*np.cos(2*th3D))*np.sin(4*p3D)*np.sin(th3D)**3 +
            np.sqrt(2)*r3D**2*(np.cos(p3D)*(21 + 55*np.cos(2*th3D) + 31*np.cos(4*th3D) + 21*np.cos(6*th3D)) - 12*np.cos(3*p3D)*(10 + 15*np.cos(2*th3D) + 7*np.cos(4*th3D))*np.sin(th3D)**2)))/(1024.*np.pi*r3D**5)

ugradb_t_a = (3*np.sqrt(21)*(20*np.cos(th3D)*np.sin(th3D)**2*(np.sqrt(2)*r3D*(-4*np.cos(3*p3D)*np.cos(th3D)**2*(19 + 21*np.cos(2*th3D)) + np.cos(p3D)*(1 + 4*np.cos(2*th3D) - 21*np.cos(4*th3D))) - 40*(5 + 7*np.cos(2*th3D))*np.sin(4*p3D)*np.sin(th3D)) +
            3*np.sqrt(5)*np.sin(2*p3D)*(47*np.sin(2*th3D) + 52*np.sin(4*th3D) + 35*np.sin(6*th3D))))/(1024.*np.pi*r3D**4)

ugradb_p_a = (3*np.sqrt(21)*((-186*np.sqrt(5)*np.cos(2*p3D) + (-125 + 492*np.sqrt(5))*np.cos(2*th3D) + (50 + 564*np.sqrt(5))*np.cos(4*th3D) + 15*(-30 + 4*np.sqrt(5) + 7*(5 + 4*np.sqrt(5))*np.cos(6*th3D)))*np.sin(th3D) +
           5*np.sqrt(2)*r3D*((333 + 596*np.cos(2*th3D) + 287*np.cos(4*th3D))*np.sin(p3D) + (147 + 332*np.cos(2*th3D) + 161*np.cos(4*th3D))*np.sin(3*p3D))*np.sin(th3D)**2 - 200*np.cos(4*p3D)*(1 + 3*np.cos(2*th3D))*(5 + 7*np.cos(2*th3D))*np.sin(th3D)**3 +
           6*np.sqrt(5)*np.cos(2*p3D)*(-6*np.sin(3*th3D) + 12*np.sin(5*th3D) + 35*np.sin(7*th3D))))/(1024.*np.pi*r3D**4)



# r
####

durdr = rderavg(vr, rad=r, exclude=False)
durdt = (1.0 / r3D) * thetaderavg(vr, order=4)
durdp = (1.0 / (r3D * np.sin(th3D))) * phideravg(vr, order=4)

# theta
#######


dutdr = rderavg(vtheta, rad=r, exclude=False)
dutdt = (1.0 / r3D) * thetaderavg(vtheta, order=4)
dutdp = (1.0 / (r3D * np.sin(th3D))) * phideravg(vtheta, order=4)

# phi
#####

dupdr = rderavg(vphi, rad=r, exclude=False)
dupdt = (1.0 / r3D) * thetaderavg(vphi, order=4)
dupdp = (1.0 / (r3D * np.sin(th3D))) * phideravg(vphi, order=4)


# r
###

dbrdr = rderavg(Br, rad=r, exclude=False)
dbrdt = (1.0 / r3D) * thetaderavg(Br, order=4)
dbrdp = (1.0 / (r3D * np.sin(th3D))) * phideravg(Br, order=4)

# theta
#######

dbtdr = rderavg(Btheta, rad=r, exclude=False)
dbtdt = (1.0 / r3D) * thetaderavg(Btheta, order=4)
dbtdp = (1.0 / (r3D * np.sin(th3D))) * phideravg(Btheta, order=4)

# phi
#####

dbpdr = rderavg(Bphi, rad=r, exclude=False)
dbpdt = (1.0 / r3D) * thetaderavg(Bphi, order=4)
dbpdp = (1.0 / (r3D * np.sin(th3D))) * phideravg(Bphi, order=4)


bgradu_r = (
    Br * durdr
    + Btheta * durdt
    + Bphi * durdp
    - ( (Btheta * vtheta + Bphi * vphi)/r3D )
)


bgradu_t = (
    Br * dutdr
    + Btheta * dutdt
    + Bphi * dutdp
    + (Btheta * vr)/r3D
    - ((Bphi * vphi) / (r3D * np.tan(th3D)))
)


bgradu_p = (
    Br * dupdr
    + Btheta * dupdt
    + Bphi * dupdp
    + (Bphi * vr) / r3D
    + (Bphi * vtheta) / (r3D * np.tan(th3D))
)


ugradb_r = (
    vr * dbrdr
    + vtheta * dbrdt
    + vphi * dbrdp
    - ( (vtheta * Btheta + vphi * Bphi)/r3D )
)


ugradb_t = (
    vr * dbtdr
    + vtheta * dbtdt
    + vphi * dbtdp
    + (vtheta * Br)/r3D
    - ((vphi * Bphi) / (r3D * np.tan(th3D)))
)


ugradb_p = (
    vr * dbpdr
    + vtheta * dbpdt
    + vphi * dbpdp
    + (vphi * Br) / r3D
    + (vphi * Btheta) / (r3D * np.tan(th3D))
)

def compare_fields(ax,xx,yy,field1, field2, phiIdx,r,theta):

    diff   = field1 - field2

    for k in range(len(ax)):
        if k == 0:
            dat = field1
        elif k == 1:
            dat = field2
        elif k == 2:
            dat = diff

        datMax = dat.max()
        datMin = dat.min()
        datCenter = 0

        divnorm = colors.TwoSlopeNorm(vmin=datMin, vcenter=datCenter, vmax=datMax)

        cont = ax[k].contourf(xx,yy,dat[phiIdx,...],norm=divnorm,cmap='seismic',levels=100)
        # cont.set_clim(vmin=datMin,vmax=datMax)
        cbar = plt.colorbar(cont)
        ax[k].plot(r[0]*np.sin(theta),r[0]*np.cos(theta),'k',lw=1)
        ax[k].plot(r[-1]*np.sin(theta),r[-1]*np.cos(theta),'k',lw=1)
        ax[k].plot([0,0],[r[0],r[-1]],'k',lw=1)
        ax[k].plot([0,0],[-r[0],-r[-1]],'k',lw=1)

        ax[k].set_aspect('equal')
        ax[k].set_axis_off()

rr,tth = np.meshgrid(r,theta)

xx = rr*np.sin(tth)
yy = rr*np.cos(tth)

phiPlot = 45
phiIdx = np.argmin(np.abs(phi*180/np.pi - phiPlot))

nrows = 6
ncols = 3

fig,axs = plt.subplots(figsize=(16,16),nrows=nrows,ncols=ncols)

for nrow in range(nrows):
    if nrow == 0:
        field1 = bgradu_r
        field2 = bgradu_r_a
    elif nrow == 1:
        field1 = bgradu_t
        field2 = bgradu_t_a
    elif nrow == 2:
        field1 = bgradu_p
        field2 = bgradu_p_a
    if nrow == 3:
        field1 = ugradb_r
        field2 = ugradb_r_a
    elif nrow == 4:
        field1 = ugradb_t
        field2 = ugradb_t_a
    elif nrow == 5:
        field1 = ugradb_p
        field2 = ugradb_p_a

    compare_fields(axs[nrow,:],xx,yy,field1,field2,phiIdx,r,theta)

plt.tight_layout()
plt.show()