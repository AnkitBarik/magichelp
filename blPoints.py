#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

from pylab import *
import math

def get_grid(a,b,nr,mapping='default'):

    bma = 0.5*(b-a)
    bpa = 0.5*(b+a)

    k = arange(1.,nr + 2)

    y = cos(pi * (k-1)/nr)

    a1 = 0.8
    a2 = 0.

    K=arctan(a1*(1.+a2))/arctan(a1*(1.-a2))
    x0=(K-1.)/(K+1.)
    lambd=arctan(a1*(1.-a2))/(1.-x0)

    if mapping == 'bay':
        x = 0.5 * (a2+tan(lambd*(y-x0))/a1) + bpa
    if mapping == 'ktl':
        x = 0.5 * arcsin(a1*y)/arcsin(a1) + bpa
    if mapping == 'default':
        x = bma * y + bpa

    return x

def get_blPts(ek,ri,ro,mapping='default'):

    r = get_grid(ri,ro,nr,mapping=mapping)

    d = sqrt(ek)

    nrBl = len(r[r < ri + d]) - 1

    return nrBl

def nextpow2(x):
    return 2**(math.ceil(math.log(x, 2)))

a = 0.35

ri = a/(1.-a)
ro = 1./(1.-a)

ek = 1e-15

for i in range(10000):
    nr = 4*i + 1

    nrBl = get_blPts(ek,ri,ro)

    if nrBl > 2:
        break

m = 4e4

np = nextpow2( 2*m )

nt = np/2.

N = nr*np*nt

print(("Nr,Nt,Np = %d,%d,%d" %(nr,nt,np)))

nPB = N*8./2**50

print(("Total memory required: %.2f PB" %nPB))

minc = 1
                                                                                
l_max = (20*nt)/30                                                                                                                                                                       
m_max = (l_max/minc)*minc                                        
lm_max = m_max*(l_max+1)/minc-m_max*(m_max-minc)/(2*minc)+(l_max+1-m_max)

nPBspec = nr * lm_max * 8./2**50

print(("Total memory required in spectral space: %.2f PB" %nPBspec))

