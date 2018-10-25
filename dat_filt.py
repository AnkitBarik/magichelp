#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

from pylab import *
from magic import *
from get_Prod import getPro
import shtns


def filt(dat_r,dat_t,dat_p,lmax = None, lCut = 10, filt_type='high'):

    dat_r = float64(dat_r)
    dat_t = float64(dat_t)
    dat_p = float64(dat_p)

    nphi = dat_r.shape[0]
    nlat = dat_r.shape[1]
    nr   = dat_r.shape[2]

    # Convert from MagIC to shtns format: nphi,ntheta,nr -> ntheta,nphi,nr

    dat_r = transpose(dat_r,(1,0,2))
    dat_t = transpose(dat_t,(1,0,2))
    dat_p = transpose(dat_p,(1,0,2))

    if lmax is None:
        lmax = nphi/3

    sh = shtns.sht(lmax = lmax,mmax=lmax,norm=shtns.sht_schmidt|shtns.SHT_NO_CS_PHASE)
    polar_opt=1e-10                                                                 
    nlat,nphi = sh.set_grid(nlat=nlat,nphi=nphi,polar_opt=polar_opt)  

    Q  = zeros([nr,sh.nlm],dtype='complex128')
    S  = zeros([nr,sh.nlm],dtype='complex128')
    T  = zeros([nr,sh.nlm],dtype='complex128')

    Q1 = zeros([nr,sh.nlm],dtype='complex128')
    S1 = zeros([nr,sh.nlm],dtype='complex128')
    T1 = zeros([nr,sh.nlm],dtype='complex128')

    # Forward transform: data -> SH
    
    for i in range(nr):
        Q[i,:],S[i,:],T[i,:] = sh.analys(dat_r[...,i],dat_t[...,i],dat_p[...,i])


    # Apply filter
    
    if filt_type == 'high':
        mask = sh.l >= lCut
    elif filt_type == 'low':
        mask = sh.l <= lCut
    else:
        print 'Filter type must be "high" or "low"'
        exit()

    Q1[:,mask] = Q[:,mask]
    S1[:,mask] = S[:,mask]
    T1[:,mask] = T[:,mask]

    # Inverse transform: SH -> data

    for k in range(nr):
         dat_r[...,k], dat_t[...,k], dat_p[...,k] = sh.synth(Q1[k,:],S1[k,:],T1[k,:])

    # Convert from shtns to MagIC format: ntheta,nphi,nr -> nphi,ntheta,nr
    
    dat_r = transpose(dat_r,(1,0,2))
    dat_t = transpose(dat_t,(1,0,2))
    dat_p = transpose(dat_p,(1,0,2))

    return dat_r, dat_t, dat_p


if __name__ == '__main__':

    G = MagicGraph(ave=True,tag='Graph')

    vr,vt,vp = filt(G.vr, G.vtheta, G.vphi)

    print vr
    print vt
    print vp

