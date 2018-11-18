# -*- coding: utf-8 -*-
from magic import MagicSetup, scanDir
from .setup import labTex, defaultCm, defaultLevels, labTex, buildSo
from .libmagic import *
import os, re, sys, time
import numpy as np
import matplotlib.pyplot as plt
from .npfile import *

#if buildSo:
#    if sys.version_info.major == 3:
#        import magic.checkReader3 as chk
#    elif  sys.version_info.major == 2:
#        import magic.checkReader2 as chk

import magic.checkReader as chk


def getEndianness(filename):
    """
    This function determines the endianness of the potential files

    :param filename: input of the filename
    :type filename: str
    :returns: the endianness of the file ('B'='big_endian' or 'l'='little_endian')
    :rtype: str
    """
    f = npfile(filename, endian='B')
    try:
        f.fort_read('i8')
        endian = 'B'
    except TypeError:
        endian = 'l'
    f.close()

    return endian


def getlm(l_max,m_max):
    l = np.arange(l_max)
    m = np.arange(m_max)
    idx = np.zeros([l_max,m_max])
    count = 0
    for i in range(l_max):
        for j in range(i):
            idx[i,j] = count
            count += 1

    idx = np.int32(idx)
    return l,m,idx

class MagicCheckpoint(MagicSetup):
    """
    This class allows to load and display the content of the checkpoint files

    """

    def __init__(self, datadir='.', tag=None, ave=False):
        """
        :param datadir: the working directory
        :type datadir: str
        :param tag: if you specify a pattern, it tries to read the corresponding files
        :type tag: str
        :param ave: read a time-averaged checkpoint file
        :type ave: bool
        """

        if ave:
            self.name = 'checkpoint_ave' 
        else:
            self.name = 'checkpoint_end' 

        if tag is not None:
            pattern = os.path.join(datadir, '%s*%s' % (self.name, tag))
            files = scanDir(pattern)
            if len(files) != 0:
                filename = files[-1]
            else:
                print('No such tag... try again')
                return

            if os.path.exists(os.path.join(datadir, 'log.%s' % tag)):
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % tag)
        else:
            pattern = os.path.join(datadir, '%s*' % self.name)
            files = scanDir(pattern)
            filename = files[-1]
            # Determine the setup
            mask = re.compile(r'.*\.(.*)')
            ending = mask.search(files[-1]).groups(0)[0]
            if os.path.exists(os.path.join(datadir, 'log.%s' % ending)):
                MagicSetup.__init__(self, datadir=datadir, quiet=True,
                                    nml='log.%s' % ending)


        if (self.radial_scheme=='CHEB'):
            l_cheb = True
        else:
            l_cheb = False

        # Determine file endianness
#        endian = getEndianness(filename)
        endian = 'B'
        Prd = chk.rstreader
        Prd.readrst(filename, endian, l_cheb, self.lm_max)
        self.minc = Prd.minc
        self.m_max = int((self.l_max/self.minc)*self.minc)
        self.n_m_max = int(self.m_max/self.minc+1)
        self.l,self.m,self.idx = getlm(self.l_max,self.m_max)
        self.l_cond_ic = Prd.l_cond_ic
        self.l_heat = Prd.l_heat
        self.l_chemical_conv = Prd.l_chemical_conv
        self.l_mag = Prd.l_mag

        self.wlm = Prd.wlm 
        self.zlm = Prd.zlm 
        self.plm = Prd.plm 
        self.dwdt = Prd.dwdt
        self.dzdt = Prd.dzdt
        self.dpdt = Prd.dpdt

        if self.l_heat:
            self.slm = Prd.slm
            self.dsdt = Prd.dsdt
        if self.l_chemical_conv:
            self.xilm = Prd.xilm
            self.dxidt = Prd.dxidt
        if self.l_mag:
            self.blm = Prd.blm
            self.jlm = Prd.jlm
            self.dbdt = Prd.dbdt
            self.djdt = Prd.djdt
            if self.l_cond_ic:
                self.blm_ic = Prd.blm_ic 
                self.jlm_ic = Prd.jlm_ic 
                self.dbdt_ic = Prd.dbdt_ic 
                self.djdt_ic = Prd.djdt_ic 
