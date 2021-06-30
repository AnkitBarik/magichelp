#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

from pylab import *
import math

def get_grid(a,b,nr,mapping='default'):

    bma = 0.5*(b-a)
    bpa = 0.5*(b+a)

    k = arange(1.,nr + 2)

    y = cos(pi * (k-1)/nr)

    a1 = 0.5
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

def get_blPts(ek,ri,ro,nr,mapping='default'):

    r = get_grid(ri,ro,nr,mapping=mapping)

    d = sqrt(ek)

    nrBl = len(r[r < ri + d]) - 1

    return nrBl


a = 0.35

ek = 1e-5

ri = a/(1.-a)
ro = 1./(1.-a)

nr = 101

figure(figsize=(16,9))

subplot(3,1,1)

mapping='default'

x = get_grid(ri,ro,nr,mapping=mapping)

print((get_blPts(ek,ri,ro,nr,mapping=mapping)))


plot(x,ones_like(x),'o')

d = sqrt(ek)

axvline(x=ri+d,color='k')
axvline(x=ro-d,color='k')

title(mapping,fontsize=30)

subplot(3,1,2)

mapping='bay'

x = get_grid(ri,ro,nr,mapping=mapping)

print((get_blPts(ek,ri,ro,nr,mapping=mapping)))


plot(x,ones_like(x),'o')

d = sqrt(ek)

axvline(x=ri+d,color='k')
axvline(x=ro-d,color='k')
title(mapping,fontsize=30)


subplot(3,1,3)

mapping='ktl'

x = get_grid(ri,ro,nr,mapping=mapping)

print((get_blPts(ek,ri,ro,nr,mapping=mapping)))


plot(x,ones_like(x),'o')

d = sqrt(ek)

axvline(x=ri+d,color='k')
axvline(x=ro-d,color='k')

title(mapping,fontsize=30)

tight_layout()
show()
