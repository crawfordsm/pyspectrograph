"""ModelSolution is a task describing the functional form for fitting an
spectrograph model to an 1-D spectroscopic data.  The model will output
the best fit spectrograph to the observed data

In the Spectrograph model, the three things that we can alter are the
x position, the y position, and focus.

HISTORY
20100614  SMC Initially Written by SM Crawford

LIMITATIONSa
20101001  Add the ability to alter n--the index of refraction

"""

import numpy as np
from PySpectrograph.Utilities import fit as ft


class ModelSolution:

    """ModelSolution is a task describing the functional form for transforming
       pixel position to wavelength for a 2-D image using a model for a
       spectrograph.

       xarr--pixel position of the spectrum
       farr--the flux values for the spectrum
       sgraph--model of a spectrograph
       spectrum--model spectrum
       yval--offset on the ccd in the y-direction
       var--array with the variance values
    """

    def __init__(self, x, y, sgraph, xlen=3162, yval=0, order=3):

        # set up the variables
        self.x = x
        self.y = y
        self.sgraph = sgraph
        self.yval = yval
        self.xlen = xlen
        self.order = order
        self.function = 'model'

        # set the parameters
        self.reset_coef()

    def reset_coef(self):
        xpos = self.sgraph.detector.xpos
        ypos = self.yval + self.sgraph.detector.ypos
        fcam = self.sgraph.camera.focallength
        self.set_spcoef(fcam, xpos, ypos)
        self.spcoef = [self.fcam, self.xpos, self.ypos]
        self.set_ndcoef()

        # the total coefficient
        self.set_coef()

    def set_spcoef(self, fcam, xpos, ypos):
        self.fcam = ft.Parameter(fcam)
        self.xpos = ft.Parameter(xpos)
        self.ypos = ft.Parameter(ypos)
        self.spcoef = [self.fcam, self.xpos, self.ypos]

    def set_coef(self):
        self.coef = self.spcoef + self.ndcoef

    def set_ndcoef(self, x=[1.00, 0.00, 0.00]):
        self.ndcoef = []
        for i in range(self.order):
            try:
                self.ndcoef.append(ft.Parameter(x[i]))
            except:
                self.ndcoef.append(ft.Parameter(0.0))

    def set_xarr(self, xarr):
        self.xarr = xarr
        if self.xarr is None:
            npix = self.sgraph.detector.width / (self.sgraph.detector.pix_size * self.sgraph.detector.xbin)
            self.xarr = np.arange(npix)
            self.xlen = len(self.xarr)
        else:
            self.xlen = len(self.xarr)

    def nd(self, x):
        v = 0
        for i in range(self.order):
            v += self.ndcoef[i]() * x ** i
        return v

    def value(self, x):
        """Return the wavelength value at x due to the model and current values for the model"""
        # these are just left here for the record
        # dx=self.sgraph.detector.xpos/(self.sgraph.detector.xbin*self.sgraph.detector.pix_scale)
        # dy=self.sgraph.detector.ypos/(self.sgraph.detector.ybin*self.sgraph.detector.pix_scale)
        fcam = self.spcoef[0]()
        dx = self.spcoef[1]()
        dy = self.spcoef[2]()
        alpha = self.sgraph.gratang
        beta = self.sgraph.gratang - self.sgraph.camang
        dw = 0.5 * 1e7 * self.sgraph.calc_resolelement(alpha, beta)
        dbeta = np.degrees(
            np.arctan(self.sgraph.detector.xbin * self.sgraph.detector.pix_size * (x - 0.5 * self.xlen + dx) / fcam))
        gamma = np.degrees(np.arctan(self.sgraph.detector.ybin * self.sgraph.detector.pix_size * (dy) / fcam))
        return 1e7 * self.sgraph.calc_wavelength(alpha, beta - dbeta, gamma=gamma) * self.nd(x)

    def erf(self, x, y):
        retrun(y - self.value(x))

    def fit_coef(self, coef):
        self.result = ft.fit(self.value, coef, self.y, x=self.x)

    def fit(self, cfit='all'):
        """Fit can be set to be several different possibilites:

           pscoef--only fit xpos, ypos, fcam
           ndcoef--only fit the index of refraction
           all--fit coef and then fit ndcoef
           both--fit coef and ndcoef
        """
        if cfit in ['pscoef', 'all']:
            self.fit_coef(self.spcoef)
        if cfit in ['ndcoef', 'all']:
            self.fit_coef(self.ndcoef)
        if cfit in ['both']:
            self.fit_coef(self.coef)
        # update the total coefficient
        self.set_coef()

    def sigma(self, x, y):
        """Return the RMS of the fit """
        return (((y - self.value(x)) ** 2).mean()) ** 0.5

    def chisq(self, x, y, err):
        """Return the chi^2 of the fit"""
        return (((y - self.value(x)) / err) ** 2).sum()
