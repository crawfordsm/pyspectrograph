"""LineFit is a task describing the functional form for transforming
pixel position to wavelength by fitting a model spectrum.  The inputs for this task
are an observed spectrum and a calibrated spectrum
and the corresponding wavelength.  The user selects an input functional form and
order for that form.  The task then calculates the coefficients for that form.
Possible options for the wavelength solution include polynomial, legendre, spline.

HISTORY
20090915  SMC Initially Written by SM Crawford
20101018  SMC Updated to use interfit

LIMITATIONS

"""

import math
import numpy as np
from scipy import optimize, interpolate


from PySpectrograph.Utilities.fit import fit, interfit


class LineFit(interfit):

    """LineSolution is a task describing the functional form for transforming
       pixel position to wavelength.

    * obsspec --observed spectrum
    * calspec --calibrated spectrum
    * function - function to be fit to the data:
                 options include polynomial, legendre, chebyshev, or spline
    * order - order of the function that is fit

    """

    def __init__(self, obsspec, calspec, function='poly', order=3, coef=None):

        # set up the spectrum
        self.obs_spec = obsspec
        self.cal_spec = calspec

        # set up the function
        self.order = order
        self.set_func(function)

        # set up the coef
        self.set_coef(coef)

    def flux(self, x):
        """Return the calibrated flux at the transformed position x"""
        return self.cal_spec.get_flux(self(x))

    def value(self, x):
        """Calculate the value of the array
        """
        return self(x)

    def errfit(self, coef, x, y, err=1):
        self.set_coe(coef)
        return (y - self.flux(x))

    def lfit(self, task=0, s=None, t=None, full_output=1, warn=False):
        self.results = optimize.leastsq(self.errfit, self.coef,
                                        args=(self.obs_spec.wavelength, self.obs_spec.flux), full_output=full_output)
        self.set_coef(self.results[0])
