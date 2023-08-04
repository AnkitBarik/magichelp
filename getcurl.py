#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

from magic import *
import numpy as np

def getRms(u, r, th, th2D, phi):

    if u.ndim == 3:
        vol = 4.0 / 3.0 * np.pi * (r.max() ** 3 - r.min() ** 3)
        uPhiInt = np.trapz(u, phi, axis=0)
    elif u.ndim == 2:
        vol = 2.0 / 3.0 * (r.max() ** 3 - r.min() ** 3)
        uPhiInt = u

    uThetaInt = np.trapz(uPhiInt * np.sin(th2D), th, axis=0)
    uRInt = -np.trapz(uThetaInt * r ** 2, r)

    return np.sqrt(uRInt / vol)


def curlVec(vr,vtheta,vphi,r3D,th3D,r):

    dupdt = 1./r3D * thetaderavg(vphi)
    dutdp = 1./(r3D * np.sin(th3D)) * phideravg(vtheta)

    dupdr = rderavg(vphi,rad=r)
    durdp = 1./(r3D * np.sin(th3D)) * phideravg(vr)

    durdt = 1/r3D * thetaderavg(vr)
    dutdr = rderavg(vtheta,rad=r)

    curl_r = dupdt - dutdp + (vphi*np.cos(th3D))/(r3D * np.sin(th3D))
    curl_t = -dupdr + durdp - vphi/r3D
    curl_p = -durdt + dutdr + vtheta/r3D


    return curl_r,curl_t,curl_p
