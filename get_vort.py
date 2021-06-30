#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

from pylab import *
from magic import *
from copy import deepcopy
from scipy.io import savemat

def getRms(u,r,th,th2D,phi):

    if u.ndim == 3:    
        vol = 4./3. * pi * (r.max()**3 - r.min()**3)
        uPhiInt = trapz(u,phi,axis=0)
    elif u.ndim == 2:
        vol = 2./3. * (r.max()**3 - r.min()**3)
        uPhiInt = u
    
    uThetaInt = trapz(uPhiInt*sin(th2D),th,axis=0)
    uRInt     =-trapz(uThetaInt*r**2,r)

    return sqrt(uRInt/vol)

def getVort(vr,vt,vp,r,th,phi):

    eta = r.min()/r.max()

    if vr.ndim == 3:
        th3D = zeros_like(vr)
        r3D  = zeros_like(vr)
        th2D = zeros_like(vr[0,...])
        r2D  = zeros_like(vr[0,...])
    
        for k in range(G.ntheta):
            th3D[:,k,:] = th[k]
            th2D[k,:] = th[k]

        for k in range(len(r)):
            r3D[...,k] = r[k]
            r2D[...,k] = r[k]

        curl_r = 1./(r3D * sin(th3D)) * ( thetaderavg(vp * sin(th3D)) - phideravg(vt) )
        curl_t = 1./r3D * ( 1./(sin(th3D)) * phideravg(vr) - rderavg(r3D * vp, eta = eta) )
        curl_p = 1./r3D * ( rderavg(r3D * vt, eta=eta) - thetaderavg(vr)  )

    elif vr.ndim == 2:
        th2D = zeros_like(vr)
        r2D  = zeros_like(vr)
      
        for k in range(G.ntheta):
            th2D[k,:] = th[k]

        for k in range(len(r)):
            r2D[...,k] = r[k]
     
        curl_r = 1./(r2D * sin(th2D)) * ( thetaderavg(vp * sin(th2D)) )
        curl_t = 1./r2D * ( - rderavg(r2D * vp, eta = eta) )
        curl_p = 1./r2D * ( rderavg(r2D * vt, eta=eta) - thetaderavg(vr)  )



    vort2 = curl_r**2 + curl_t**2 + curl_p**2

    return th2D, vort2


dirs = sort(os.walk('.').next()[1])[:-1]

delDirs = ['2.15','2.25','2.35']

for i in range(len(delDirs)):
    Idx = where(dirs==delDirs[i])
    dirs = delete(dirs,Idx)

ro = -float32(dirs)

rM = [1.6,1.0,0.7,0.6,0.58]
l = zeros([len(rM),len(ro)])
#rM = [1.5,1.4,1.3,1.0,0.7,0.5]
#rM = [1.6]

for rIt in range(len(rM)):
    print("r < ",rM[rIt])
    print("------------")
    for i in progressbar(list(range(len(dirs)))):

        os.chdir(dirs[i] + '/Graph')

        G = MagicGraph(ave=True)

        vr = deepcopy(G.vr)
        vp = deepcopy(G.vphi)
        vt = deepcopy(G.vtheta)
        r = deepcopy(G.radius)
        th = deepcopy(G.colatitude)
        phi = linspace(0.,2*pi,G.npI)


        mask = r < rM[rIt]

        r = r[mask]
        vr = vr[...,mask]
        vt = vt[...,mask]
        vp = vp[...,mask]

#        vp = vp.mean(axis=0)
#        vr = vr.mean(axis=0)
#        vt = vt.mean(axis=0)
        vp -= vp.mean(axis=0)
        vr -= vr.mean(axis=0)
        vt -= vt.mean(axis=0)

        th2D, vort2 = getVort(vr,vt,vp,r,th,phi)

        vortRms = getRms(vort2,r,th,th2D,phi)

        u2 = vr**2 + vt**2 + vp**2

        uRms = getRms(u2,r,th,th2D,phi)

        l[rIt,i] = uRms/vortRms
        
        del G

        os.chdir('../..')

#    plot(ro,l,'-o',label='s < %.1f' %rM[rIt])

#axvline(x=-2.3,ls='--',color='k')
#legend()
#show()
#savemat('lOmegaNA.mat',{'NA':False,'ro':ro,'rM':rM,'l':l})

