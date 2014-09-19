"""Wavelength Solution is a task describing the functional form for transforming
pixel position to wavelength.  The inputs for this task are the given pixel position
and the corresponding wavelength.  The user selects an input functional form and
order for that form.  The task then calculates the coefficients for that form.
Possible options for the wavelength solution include polynomial, legendre, spline.

HISTORY
20090915  SMC Initially Written by SM Crawford

LIMITATIONS
20090915 SMC Need to add legendre, spline functions

"""

import numpy as np


from .LineSolution import LineSolution
from .ModelSolution import ModelSolution


class WavelengthSolution:

    """Wavelength Solution is a task describing the functional form for transforming
       pixel position to wavelength.
    """
    func_options = ['poly', 'polynomial', 'spline', 'legendre', 'chebyshev', 'model']

    def __init__(self, x, w, function='poly', order=3, niter=5, thresh=3,
                 sgraph=None, cfit='both', xlen=3162, yval=0):
        self.sgraph = sgraph
        self.function = function
        self.order = order
        self.niter = niter
        self.thresh = thresh
        self.cfit = cfit
        self.xlen = xlen
        self.yval = yval

        self.set_array(x, w)
        self.set_func()

    def set_array(self, x, w):
        self.x_arr = x
        self.w_arr = w

    def set_thresh(self, thresh):
        self.thresh = thresh

    def set_niter(self, niter):
        self.niter = niter

    def set_func(self):
        if self.function in ['poly', 'polynomial', 'spline', 'legendre', 'chebyshev']:
            self.func = LineSolution(self.x_arr, self.w_arr, function=self.function,
                                     order=self.order, niter=self.niter, thresh=self.thresh)
        if self.function == 'model':
            self.func = ModelSolution(self.x_arr, self.w_arr, sgraph=self.sgraph,
                                      xlen=self.xlen, yval=self.yval, order=self.order)

    def fit(self):
        if self.function in ['poly', 'polynomial', 'spline', 'legendre', 'chebyshev']:
            self.func.interfit()
            self.coef = self.func.coef
        if self.function in ['model']:
            self.func.fit(cfit=self.cfit)
            self.coef = np.array([c() for c in self.func.coef])
        # self.set_coef(coef)

    def set_coef(self, coef):
        if self.function in ['poly', 'polynomial', 'spline', 'legendre', 'chebyshev']:
            self.func.coef = coef
            self.coef = self.func.coef
        if self.function in ['model']:
            for i in range(len(self.func.coef)):
                self.func.coef[i].set(coef[i])
            self.coef = np.array([c() for c in self.func.coef])

    def value(self, x):
        return self.func.value(x)

    def invvalue(self, w):
        """Given a wavelength, return the pixel position

        """
        return w

    def sigma(self, x, y):
        """Return the RMS of the fit """
        return (((y - self.value(x)) ** 2).mean()) ** 0.5

    def chisq(self, x, y, err):
        """Return the chi^2 of the fit"""
        return (((y - self.value(x)) / err) ** 2).sum()
