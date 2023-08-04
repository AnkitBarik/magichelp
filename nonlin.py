#!/usr/bin/env python3

from magic import *
import numpy as np
from getcurl import curlVec
import os


class getnl(MagicSetup):
    def __init__(
        self,
        ivar=None,
        tag=None,
        datadir=".",
        term="ugradb",
        comp="r",
        ave=False,
        precision=np.float32,
        fluct=False,
    ):

        if os.path.exists("log.%s" % tag):
            MagicSetup.__init__(self, datadir=datadir, quiet=True, nml="log.%s" % tag)

        self.gr = MagicGraph(
            ivar=ivar, tag=tag, datadir=datadir, ave=ave, precision=precision
        )

        self.term = term
        self.comp = comp

        r3D = np.zeros_like(self.gr.vr)
        th3D = np.zeros_like(self.gr.vr)

        for i in range(len(self.gr.radius)):
            r3D[:, :, i] = self.gr.radius[i]

        for i in range(len(self.gr.colatitude)):
            th3D[:, i, :] = self.gr.colatitude[i]

        if fluct and term != "ugradu":
            self.gr.vr -= self.gr.vr.mean(axis=0)
            self.gr.vtheta -= self.gr.vtheta.mean(axis=0)
            self.gr.vphi -= self.gr.vphi.mean(axis=0)

            self.gr.Br -= self.gr.Br.mean(axis=0)
            self.gr.Btheta -= self.gr.Btheta.mean(axis=0)
            self.gr.Bphi -= self.gr.Bphi.mean(axis=0)
        elif fluct and term == "ugradu":
            self.gr.vr -= self.gr.vr.mean(axis=0)
            self.gr.vtheta -= self.gr.vtheta.mean(axis=0)
            self.gr.vphi -= self.gr.vphi.mean(axis=0)

        ########################
        # Velocity derivatives #
        ########################

        if self.term == "bgradu" or self.term == "curluxb" or self.term == "ugradu":

            # r
            ###

            if self.comp == "r" or self.comp == "tot" or self.comp == "all":

                self.durdr = rderavg(self.gr.vr, rad=self.gr.radius, exclude=False)
                self.durdt = (1.0 / r3D) * thetaderavg(self.gr.vr, order=4)
                self.durdp = (1.0 / (r3D * np.sin(th3D))) * phideravg(self.gr.vr, order=4)

            # theta
            #######

            if self.comp == "theta" or self.comp == "tot" or self.comp == "all":

                self.dutdr = rderavg(self.gr.vtheta, rad=self.gr.radius, exclude=False)
                self.dutdt = (1.0 / r3D) * thetaderavg(self.gr.vtheta, order=4)
                self.dutdp = (1.0 / (r3D * np.sin(th3D))) * phideravg(self.gr.vtheta, order=4)

            # phi
            #####

            if self.comp == "phi" or self.comp == "tot" or self.comp == "all":

                self.dupdr = rderavg(self.gr.vphi, rad=self.gr.radius, exclude=False)
                self.dupdt = (1.0 / r3D) * thetaderavg(self.gr.vphi, order=4)
                self.dupdp = (1.0 / (r3D * np.sin(th3D))) * phideravg(self.gr.vphi, order=4)

        ##########################################################

        ##############################
        # Magnetic field derivatives #
        ##############################

        if self.term == "ugradb" or self.term == "curluxb":

            # r
            ###

            if self.comp == "r" or self.comp == "all":

                self.dbrdr = rderavg(self.gr.Br, rad=self.gr.radius, exclude=False)
                self.dbrdt = (1.0 / r3D) * thetaderavg(self.gr.Br, order=4)
                self.dbrdp = (1.0 / (r3D * np.sin(th3D))) * phideravg(self.gr.Br, order=4)

            # theta
            #######

            if self.comp == "theta" or self.comp == "all":

                self.dbtdr = rderavg(self.gr.Btheta, rad=self.gr.radius, exclude=False)
                self.dbtdt = (1.0 / r3D) * thetaderavg(self.gr.Btheta, order=4)
                self.dbtdp = (1.0 / (r3D * np.sin(th3D))) * phideravg(self.gr.Btheta, order=4)

            # phi
            #####

            if self.comp == "phi" or self.comp == "all":

                self.dbpdr = rderavg(self.gr.Bphi, rad=self.gr.radius, exclude=False)
                self.dbpdt = (1.0 / r3D) * thetaderavg(self.gr.Bphi, order=4)
                self.dbpdp = (1.0 / (r3D * np.sin(th3D))) * phideravg(self.gr.Bphi, order=4)

        #################
        #   B.grad U    #
        #################

        if self.term == "bgradu":

            if self.comp == "r" or self.comp == "all":
                self.bgradu_r = (
                    self.gr.Br * self.durdr
                    + self.gr.Btheta * self.durdt
                    + self.gr.Bphi * self.durdp
                    - ( (self.gr.Btheta * self.gr.vtheta + self.gr.Bphi * self.gr.vphi)/r3D )
                )

            if self.comp == "theta" or self.comp == "all":
                self.bgradu_t = (
                    self.gr.Br * self.dutdr
                    + self.gr.Btheta * self.dutdt
                    + self.gr.Bphi * self.dutdp
                    + (self.gr.Btheta * self.gr.vr)/r3D
                    - ((self.gr.Bphi * self.gr.vphi) / (r3D * np.tan(th3D)))
                )

            if self.comp == "phi" or self.comp == "all":
                self.bgradu_p = (
                    self.gr.Br * self.dupdr
                    + self.gr.Btheta * self.dupdt
                    + self.gr.Bphi * self.dupdp
                    + (self.gr.Bphi * self.gr.vr) / r3D
                    + (self.gr.Bphi * self.gr.vtheta) / (r3D * np.tan(th3D))
                )

        #################
        #   U.grad B    #
        #################

        if self.term == "ugradb":

            if self.comp == "r" or self.comp == "all":
                self.ugradb_r = (
                    self.gr.vr * self.dbrdr
                    + self.gr.vtheta * self.dbrdt
                    + self.gr.vphi * self.dbrdp
                    - ( (self.gr.vtheta * self.gr.Btheta + self.gr.vphi * self.gr.Bphi)/r3D )
                )

            if self.comp == "theta" or self.comp == "all":
                self.ugradb_t = (
                    self.gr.vr * self.dbtdr
                    + self.gr.vtheta * self.dbtdt
                    + self.gr.vphi * self.dbtdp
                    + (self.gr.vtheta * self.gr.Br)/r3D
                    - ((self.gr.vphi * self.gr.Bphi) / (r3D * np.tan(th3D)))
                )

            if self.comp == "phi" or self.comp == "all":
                self.ugradb_p = (
                    self.gr.vr * self.dbpdr
                    + self.gr.vtheta * self.dbpdt
                    + self.gr.vphi * self.dbpdp
                    + (self.gr.vphi * self.gr.Br) / r3D
                    + (self.gr.vphi * self.gr.Btheta) / (r3D * np.tan(th3D))
                )

        #################
        #   curl(UxB)   #
        #################

        if self.term == "curluxb":
            uxb_r = -self.gr.Btheta * self.gr.vphi + self.gr.Bphi * self.gr.vtheta
            uxb_t =  self.gr.Br     * self.gr.vphi - self.gr.Bphi * self.gr.vr
            uxb_p =  self.gr.Btheta * self.gr.vr   - self.gr.Br   * self.gr.vtheta

            self.curluxb_r, self.curluxb_t, self.curluxb_p = curlVec(uxb_r,uxb_t,uxb_p,
                                                                     r3D, th3D, self.gr.radius)

        #################
        #   U.grad U    #
        #################

        if self.term == "ugradu":

            if self.comp == "r" or self.comp == "all":
                self.ugradu_r = (
                    self.gr.vr * self.durdr
                    + self.gr.vtheta * self.durdt
                    + self.gr.vphi * self.durdp
                    - ((self.gr.vtheta * self.gr.vtheta + self.gr.vphi * self.gr.vphi)/r3D
                    )
                )

            if self.comp == "theta" or self.comp == "all":
                self.ugradu_t = (
                    self.gr.vr * self.dutdr
                    + self.gr.vtheta * self.dutdt
                    + self.gr.vphi * self.dutdp
                    - ((self.gr.vphi * self.gr.vphi) / (r3D * np.tan(th3D)))
                )

            if self.comp == "phi" or self.comp == "all":
                self.ugradu_p = (
                    self.gr.vr * self.dupdr
                    + self.gr.vtheta * self.dupdt
                    + self.gr.vphi * self.dupdp
                    + (self.gr.vphi * self.gr.vr) / r3D
                    + (self.gr.vphi * self.gr.vtheta) / (r3D * np.tan(th3D))
                )

            if self.comp == "tot":
                self.ugradu_r = (
                    self.gr.vr * self.durdr
                    + self.gr.vtheta * self.durdt
                    + self.gr.vphi * self.durdp
                    - (
                        (self.gr.vtheta * self.gr.vtheta + self.gr.vphi * self.gr.vphi)
                        / r3D
                    )
                )
                self.ugradu_t = (
                    self.gr.vr * self.dutdr
                    + self.gr.vtheta * self.dutdt
                    + self.gr.vphi * self.dutdp
                    - ((self.gr.vphi * self.gr.vphi) / (r3D * np.tan(th3D)))
                )
                self.ugradu_p = (
                    self.gr.vr * self.dupdr
                    + self.gr.vtheta * self.dupdt
                    + self.gr.vphi * self.dupdp
                    + (self.gr.vphi * self.gr.vr) / r3D
                    + (self.gr.vphi * self.gr.vtheta) / (r3D * np.tan(th3D))
                )
                self.ugradu_tot = np.sqrt(
                    self.ugradu_r ** 2 + self.ugradu_t ** 2 + self.ugradu_p ** 2
                )

        #####################
        #   Coriolis Force  #
        #####################

        if self.term == "Coriolis":

            fac = 2.0 / self.ek

            if self.comp == "r" or self.comp == "all":
                self.cor_r = fac * self.gr.vphi * np.sin(th3D)
            if self.comp == "theta" or self.comp == "all":
                self.cor_t = fac * self.gr.vphi * np.cos(th3D)
            if self.comp == "phi" or self.comp == "all":
                self.cor_p = -fac * (
                    self.gr.vtheta * np.cos(th3D) + self.gr.vr * np.sin(th3D)
                )
            if self.comp == "tot":
                self.cor_r = fac * self.gr.vphi * np.sin(th3D)
                self.cor_t = fac * self.gr.vphi * np.cos(th3D)
                self.cor_p = -fac * (
                    self.gr.vtheta * np.cos(th3D) + self.gr.vr * np.sin(th3D)
                )
                self.cor_tot = np.sqrt(self.cor_r ** 2 + self.cor_t ** 2 + self.cor_p ** 2)

        #####################
        #  Pressure gradient#
        #####################

        if self.term == "pressure":

            if self.comp == "r" or self.comp == "all":
                self.dpdr = -rderavg(
                    self.gr.pre, rad=self.gr.radius, exclude=False
                )
            if self.comp == "theta" or self.comp == "all":
                self.dpdt = -(1.0 / r3D) * thetaderavg(self.gr.pre, order=4)
            if self.comp == "phi" or self.comp == "all":
                self.dpdp = -(1.0 / (r3D * np.sin(th3D))) * phideravg(self.gr.pre, order=4)

    # def surf(
    #     self, r=0.5, cm="seismic", tit=True, cbar=True, vmax=None, vmin=None, levels=50
    # ):

    #     figure(figsize=(9, 4.5))

    #     r /= 1 - self.gr.radratio
    #     ind = nonzero(
    #         where(abs(self.gr.radius - r) == min(abs(self.gr.radius - r)), 1, 0)
    #     )
    #     indPlot = ind[0][0]
    #     rad = self.gr.radius[indPlot] * (1.0 - self.gr.radratio)

    #     datPlot = symmetrize(self.dat[:, :, indPlot], self.gr.minc)

    #     phi = linspace(-pi, pi, self.gr.nphi)
    #     theta = linspace(pi / 2, -pi / 2, self.gr.ntheta)
    #     pphi, ttheta = mgrid[
    #         -pi : pi : self.gr.nphi * 1j, pi / 2.0 : -pi / 2.0 : self.gr.ntheta * 1j
    #     ]

    #     xxout, yyout = hammer2cart(theta, -pi)
    #     xxin, yyin = hammer2cart(theta, pi)
    #     x, y = hammer2cart(ttheta, pphi)

    #     plot(xxin, yyin, "k-")
    #     plot(xxout, yyout, "k-")
    #     axis("off")

    #     if vmax is not None or vmin is not None:
    #         cs = N.linspace(vmin, vmax, levels)
    #         contourf(x, y, datPlot, cs, cmap=cm, extend="both")
    #     else:
    #         contourf(x, y, datPlot, levels, cmap=cm)

    #     if cbar:
    #         colorbar()

    #     if tit:
    #         if self.term == "bgradu":

    #             if self.comp == "r":
    #                 title(
    #                     r"$(\mathbf{B}\cdot\nabla \mathbf{U})_r$ $r/r_o = %.2f$" % rad
    #                 )
    #             if self.comp == "theta":
    #                 title(
    #                     r"$(\mathbf{B}\cdot\nabla \mathbf{U})_{\theta}$ $r/r_o = %.2f$"
    #                     % rad
    #                 )
    #             if self.comp == "phi":
    #                 title(
    #                     r"$(\mathbf{B}\cdot\nabla \mathbf{U})_{\phi}$ $r/r_o = %.2f$"
    #                     % rad
    #                 )

    #         if self.term == "ugradb":

    #             if self.comp == "r":
    #                 title(
    #                     r"$(\mathbf{U}\cdot\nabla \mathbf{B})_r$ $r/r_o = %.2f$" % rad
    #                 )
    #             if self.comp == "theta":
    #                 title(
    #                     r"$(\mathbf{U}\cdot\nabla \mathbf{B})_{\theta}$ $r/r_o = %.2f$"
    #                     % rad
    #                 )
    #             if self.comp == "phi":
    #                 title(
    #                     r"$(\mathbf{U}\cdot\nabla \mathbf{B})_{\phi}$ $r/r_o = %.2f$"
    #                     % rad
    #                 )

    #         if self.term == "curluxb":

    #             if self.comp == "r":
    #                 title(
    #                     r"$\nabla \times (\mathbf{U} \times \mathbf{B})_r$ $r/r_o = %.2f$"
    #                     % rad
    #                 )
    #             if self.comp == "theta":
    #                 title(
    #                     r"$\nabla \times (\mathbf{U} \times \mathbf{B})_{\theta}$ $r/r_o = %.2f$"
    #                     % rad
    #                 )
    #             if self.comp == "phi":
    #                 title(
    #                     r"$\nabla \times (\mathbf{U} \times \mathbf{B})_{\phi}$ $r/r_o = %.2f$"
    #                     % rad
    #                 )

    #         if self.term == "ugradu":

    #             if self.comp == "r":
    #                 title(
    #                     r"$(\mathbf{U}\cdot\nabla \mathbf{U})_r$ $r/r_o = %.2f$" % rad
    #                 )
    #             if self.comp == "theta":
    #                 title(
    #                     r"$(\mathbf{U}\cdot\nabla \mathbf{U})_{\theta}$ $r/r_o = %.2f$"
    #                     % rad
    #                 )
    #             if self.comp == "phi":
    #                 title(
    #                     r"$(\mathbf{U}\cdot\nabla \mathbf{U})_{\phi}$ $r/r_o = %.2f$"
    #                     % rad
    #                 )

    #     datMax = max(abs(datPlot.max()), abs(datPlot.min()))

    #     clim(-datMax, datMax)

    #     tight_layout()
    #     show()

    # def slice(
    #     self, lon=0, cm="seismic", tit=True, cbar=True, vmax=None, vmin=None, levels=50
    # ):

    #     figure(figsize=(6, 10))
    #     phi = linspace(0.0, 360, self.gr.nphi)
    #     indPlot = where(abs(lon - phi) == min(abs(lon - phi)))[0][0]

    #     datPlot = self.dat[indPlot, :, :]

    #     rr, tth = meshgrid(self.gr.radius, self.gr.colatitude)

    #     xx = rr * np.sin(tth)
    #     yy = rr * np.cos(tth)

    #     plot(
    #         self.gr.radius[0] * np.sin(self.gr.colatitude),
    #         self.gr.radius[0] * np.cos(self.gr.colatitude),
    #         "k-",
    #     )
    #     plot(
    #         self.gr.radius[-1] * np.sin(self.gr.colatitude),
    #         self.gr.radius[-1] * np.cos(self.gr.colatitude),
    #         "k-",
    #     )
    #     plot([0.0, 0], [self.gr.radius[-1], self.gr.radius[0]], "k-")
    #     plot([0.0, 0], [-self.gr.radius[-1], -self.gr.radius[0]], "k-")

    #     if vmax is not None or vmin is not None:
    #         cs = linspace(vmin, vmax, levels)
    #         contourf(xx, yy, datPlot, cs, cmap=cm, extend="both")
    #     else:
    #         contourf(xx, yy, datPlot, levels, cmap=cm)

    #     if cbar:
    #         colorbar()

    #     if tit:
    #         if self.term == "bgradu":

    #             if self.comp == "r":
    #                 title(
    #                     r"$(\mathbf{B}\cdot\nabla \mathbf{U})_r$ $\phi = %d$"
    #                     % phi[indPlot]
    #                 )
    #             if self.comp == "theta":
    #                 title(
    #                     r"$(\mathbf{B}\cdot\nabla \mathbf{U})_{\theta}$ $\phi = %d$"
    #                     % phi[indPlot]
    #                 )
    #             if self.comp == "phi":
    #                 title(
    #                     r"$(\mathbf{B}\cdot\nabla \mathbf{U})_{\phi}$ $\phi = %d$"
    #                     % phi[indPlot]
    #                 )

    #         if self.term == "ugradb":

    #             if self.comp == "r":
    #                 title(
    #                     r"$(\mathbf{U}\cdot\nabla \mathbf{B})_r$ $\phi = %d$"
    #                     % phi[indPlot]
    #                 )
    #             if self.comp == "theta":
    #                 title(
    #                     r"$(\mathbf{U}\cdot\nabla \mathbf{B})_{\theta}$ $\phi = %d$"
    #                     % phi[indPlot]
    #                 )
    #             if self.comp == "phi":
    #                 title(
    #                     r"$(\mathbf{U}\cdot\nabla \mathbf{B})_{\phi}$ $\phi = %d$"
    #                     % phi[indPlot]
    #                 )

    #         if self.term == "curluxb":

    #             if self.comp == "r":
    #                 title(
    #                     r"$\nabla \times (\mathbf{U} \times \mathbf{B})_r$ $\phi = %d$"
    #                     % phi[indPlot]
    #                 )
    #             if self.comp == "theta":
    #                 title(
    #                     r"$\nabla \times (\mathbf{U} \times \mathbf{B})_{\theta}$ $\phi = %d$"
    #                     % phi[indPlot]
    #                 )
    #             if self.comp == "phi":
    #                 title(
    #                     r"$\nabla \times (\mathbf{U} \times \mathbf{B})_{\phi}$ $\phi = %d$"
    #                     % phi[indPlot]
    #                 )

    #         if self.term == "ugradu":

    #             if self.comp == "r":
    #                 title(
    #                     r"$(\mathbf{U}\cdot\nabla \mathbf{U})_r$ $\phi = %d$"
    #                     % phi[indPlot]
    #                 )
    #             if self.comp == "theta":
    #                 title(
    #                     r"$(\mathbf{U}\cdot\nabla \mathbf{U})_{\theta}$ $\phi = %d$"
    #                     % phi[indPlot]
    #                 )
    #             if self.comp == "phi":
    #                 title(
    #                     r"$(\mathbf{U}\cdot\nabla \mathbf{U})_{\phi}$ $\phi = %d$"
    #                     % phi[indPlot]
    #                 )

    #     datMax = max(abs(datPlot.max()), abs(datPlot.min()))

    #     clim(-datMax, datMax)

    #     axis("off")
    #     tight_layout()
    #     show()
