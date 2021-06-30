#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

from pylab import *
from magic import *
from evtk.hl import gridToVTK
from get_Prod import getPro
from dat_filt import filt
from copy import deepcopy

def get_grid(r,theta,phi,nr,nt,np):

    r3D  = zeros([np,nt,nr])
    th3D = zeros([np,nt,nr])
    p3D  = zeros([np,nt,nr])

    for i in range(nr):
        r3D[...,i] = r[i]
    for j in range(nt):
        th3D[:,j,:] = theta[j]
    for k in range(np):
        p3D[k,...] = phi[k]

    s = r3D * sin(th3D)
    x = s * cos(p3D)
    y = s * sin(p3D)
    z = r3D * cos(th3D)

    return r3D,th3D,p3D, x,y,z, s

def get_cart(vr,vt,vp,r3D,th3D,p3D):
    
    vs = vr * sin(th3D) + vt * cos(th3D)
    vz = vr * cos(th3D) - vt * sin(th3D)

    vx = vs * cos(p3D) - vp * sin(p3D)
    vy = vs * sin(p3D) + vp * cos(p3D)

    return vx,vy,vz


s1 = getPro(ave=True,tag='Graph',term='ugradu',comp='all')

r      = deepcopy(s1.gr.radius)
theta  = deepcopy(s1.gr.colatitude)
nr     = deepcopy(s1.gr.nr)
ntheta = deepcopy(s1.gr.ntheta)
nphi   = deepcopy(s1.gr.npI)

phi = linspace(0.,2*pi,nphi)

r3D,th3D,p3D, x,y,z, s = get_grid(r,theta,phi,nr,ntheta,nphi)

udu_r, udu_t, udu_p = filt(s1.ugradu_r,s1.ugradu_t,s1.ugradu_p)

#udu_x,udu_y,udu_z = get_cart(s3.ugradu_r,s4.ugradu_t,s5.ugradu_p,r3D,th3D,p3D)
udu_x,udu_y,udu_z = get_cart(udu_r, udu_t, udu_p,r3D,th3D,p3D)

gridToVTK("ugradu",x,y,z,pointData= {"radius":r3D,
                                   "cyl_radius":s,
                                   "ugradu":(udu_x, udu_y,udu_z)})

del s1

s2 = getPro(ave=True,tag='Graph',term='Coriolis',comp='all')
cor_r,cor_t,cor_p = filt(s2.cor_r,s2.cor_t,s2.cor_p)

cor_x,cor_y,cor_z = get_cart(cor_r,cor_t,cor_p,r3D,th3D,p3D)

gridToVTK("Coriolis",x,y,z,pointData= {"radius":r3D,
                                   "cyl_radius":s,
                                   "coriolis force":(cor_x, cor_y,cor_z)})

del s2

s3 = getPro(ave=True,tag='Graph',term='pressure',comp='all')

dp_r,dp_t,dp_p = filt(s3.dpdr,s3.dpdt,s3.dpdp)

dp_x,dp_y,dp_z = get_cart(dp_r,dp_t,dp_p,r3D,th3D,p3D)

gridToVTK("pressure",x,y,z,pointData= {"radius":r3D,
                                   "cyl_radius":s,
                                   "pressure gradient":(dp_x, dp_y,dp_z)})
