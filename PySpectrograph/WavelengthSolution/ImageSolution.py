"""ImageSolution is a task describing the functional form for fitting an
spectrograph model to an 2-D spectroscopic data.  The model will output
the best fit spectrograph to the observed data

In the Spectrograph model, the three things that we can alter are the
x position, the y position, and focus.

HISTORY
20100614  SMC Initially Written by SM Crawford

LIMITATIONSa
20101001  Add the ability to alter n--the index of refraction

"""

import math
import numpy as np

import fit as ft
from CreateImage import CreateImage

import pylab as pl


class fitimage:

    """Fit legendre polynomials to a function"""

    def __init__(self, order=3):
        self.order = order
        self.setcoef()

    def setcoef(self, x=1):
        self.coef = []
        for i in range(self.order):
            try:
                self.coef.append(ft.Parameter(x[i]))
            except TypeError:
                self.coef.append(ft.Parameter(x))

    def legendre(self, x):
        v = 0
        for i in range(self.order):
            v += self.coef[i]() * legendre(i)(x)
        return v

    def fit(self, data, x=None, var=None):
        return ft.fit(self.legendre, self.coef, data, x=x, var=var)


class ImageSolution:

    """ImageSolution is a task describing the functional form for transforming
       pixel position to wavelength for a 2-D image using a model for a
       spectrograph.

       data--2-D array to be fit
    """

    def __init__(self, data, sgraph, spectrum):

        # set up the variables
        self.data = data
        self.sgraph = sgraph
        self.spectrum = spectrum

        # set up the sizes
        self.ylen = len(self.data)
        self.xlen = len(self.data[0])

        # set the parameters
        self.xpos = ft.Parameter(self.sgraph.detector.xpos)
        self.ypos = ft.Parameter(self.sgraph.detector.ypos)
        self.fcam = ft.Parameter(self.sgraph.camera.focallength)
        self.nd0 = ft.Parameter(1)
        self.nd1 = ft.Parameter(0.00)
        self.nd2 = ft.Parameter(0.00)
        # self.xpos=ft.Parameter(2*0.015*8.169)
        # self.fcam=ft.Parameter(327.85)
        #self.coef=[self.xpos, self.ypos, self.fcam]
        self.coef = [self.nd0, self.nd1, self.nd2]
        # print self.sgraph.gratingequation

        j = int(self.ylen / 2)
        xarr = np.arange(self.xlen)
        yarr = np.arange(self.ylen)
        self.data1 = self.data[j, :]
        self.err1 = abs(self.data1) + self.data1.mean()
        # print self.makewave(500)
        # print self.xpos(), self.ypos(), self.fcam(), self.nd0()
        # print self.nd0(), self.nd1(), self.nd2()
        #print (self.data1-self.makeflux(xarr)).sum()
        self.fit(self.makeflux, self.coef, self.data1, var=self.err1)

        warr = self.makewave(xarr)
        # print self.makewave(500)
        # print self.xpos(), self.ypos(), self.fcam()
        # print self.nd0(), self.nd1(), self.nd2()
        #print (self.data1-self.makeflux(xarr)).sum()

        # check the results
        pl.figure()
        pl.plot(self.spectrum.wavelength, self.spectrum.flux * self.data.max() / self.spectrum.flux.max())
        pl.plot(warr, self.data[j, :])
        pl.show()

    def makend(self, x):
        return self.nd0() + self.nd1() * x + self.nd2() * x ** 2

    def makewave(self, x):
        dx = self.xpos()
        dy = self.ypos()
        alpha = self.sgraph.gratang
        beta = self.sgraph.gratang - self.sgraph.camang
        dw = 0.5 * 1e7 * self.sgraph.calc_resolelement(alpha, beta)
        # dx=self.sgraph.detector.xpos/(self.sgraph.detector.xbin*self.sgraph.detector.pix_scale)
        # dy=self.sgraph.detector.ypos/(self.sgraph.detector.ybin*self.sgraph.detector.pix_scale)
        dbeta = np.degrees(np.arctan(self.sgraph.detector.xbin *
                                     self.sgraph.detector.pix_size *
                                     (x -
                                      0.5 *
                                      self.xlen +
                                      dx) /
                                     self.fcam()))
        gamma = np.degrees(np.arctan(self.sgraph.detector.ybin * self.sgraph.detector.pix_size * (dy) / self.fcam()))
        return 1e7 * self.sgraph.calc_wavelength(alpha, beta - dbeta, gamma=gamma) * self.makend(x)

    def makeflux(self, x):
        return self.data.max() / self.spectrum.flux.max() * self.spectrum.get_flux(self.makewave(x))

    def fit(self, func, coef, data, var):
        return ft.fit(func, coef, data, var=var)
