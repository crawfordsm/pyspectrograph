
"""FIT.PY--General fitting routines

fit--fit is a general fitting routine using the scipy.optimize.leastsq
technique to fit a function.  This example has been taken from the scipy
cookbook.

"""

import numpy as np
from scipy import optimize, interpolate
from scipy.special import legendre, chebyt


class power:

    """A class to produce a polynomial term of power n.

       This has similar behavior to scipy.special.legendre

    """

    def __init__(self, n):
        self.n = n

    def __call__(self, x):
        return x ** self.n


class Parameter:

    def __init__(self, value):
        self.value = value

    def set(self, value):
        self.value = value

    def __call__(self):
        return self.value


def fit(function, parameters, y, x=None, var=1, warn=False):
    def f(params):
        i = 0
        for p in parameters:
            p.set(params[i])
            i += 1
        return (y - function(x)) / var

    if x is None:
        x = np.arange(y.shape[0])
    p = [param() for param in parameters]
    return optimize.leastsq(f, p, full_output=1, warning=warn)


class curfit:

    """Given an x and y data arrays,  find the best fitting curve

    * x - list or array of x data
    * y - list or array of y data
    * yerr - error on y data
    * coef - Initial coefficients for fit
    * function - function to be fit to the data:
                 options include polynomial, legendre, chebyshev, or spline
    * order - order of the function that is fit


    """

    def __init__(self, x, y, yerr=None, coef=None, function='poly', order=3):

        # set up the variables
        self.x = x
        self.y = y
        if yerr is None:
            self.yerr = 1
        else:
            self.yerr = yerr
        self.order = order

        self.set_func(function)
        self.set_coef(coef)

    def set_coef(self, coef=None):
        """Set the coefficients for the fits for poly, legendre, and chebyshev"""
        if coef is None:
            coef = np.ones(self.order + 1)
        if isinstance(coef, np.ndarray):
            self.coef = coef
        elif isinstance(coef, list):
            self.coef = np.array(coef)
        else:
            self.coef = np.array([coef])

    def set_func(self, function):
        """Set the function that will be used.
        * function - name of function to be used

        It will throw an error if an inappropriate function is given
        """
        self.function = function
        if self.function == 'poly' or self.function == 'polynomial' or self.function == 'power':
            self.func = power
        elif self.function == 'legendre':
            self.func = legendre
        elif self.function == 'chebyshev':
            self.func = chebyt
        elif self.function == 'spline':
            self.func = None
        else:
            msg = '%s is not a valid function' % self.function
            raise Exception(msg)

    def set_weight(self, err):
        """Set the weighting for spline fitting """
        if isinstance(err, np.ndarray):
            if err.any() != 0:
                self.weight = 1 / err
                return
        self.weight = None

    def __call__(self, x):
        """Return the value of the function evaluated at x"""
        if self.function == 'spline':
            return interpolate.splev(x, self.coef, der=0)
        v = x * 0.0
        for i in range(self.order + 1):
            v += self.coef[i] * self.func(i)(x)
        return v

    def erf(self, coef, x, y, v):
        """Error function to be minimized in least-squares fit"""
        self.set_coef(coef)
        return (y - self.__call__(x)) / v

    def sigma(self, x, y):
        """Return the RMS of the fit """
        return (((y - self(x)) ** 2).mean()) ** 0.5

    def chisq(self, x, y, err):
        """Return the chi^2 of the fit"""
        return (((y - self(x)) / err) ** 2).sum()

    def fit(self, task=0, s=None, t=None, full_output=1, warn=False):
        """Fit the function to the data"""
        if self.function == 'spline':
            self.set_weight(self.yerr)
            self.results = interpolate.splrep(
                self.x,
                self.y,
                w=self.weight,
                task=0,
                s=None,
                t=None,
                k=self.order,
                full_output=full_output)
            # w=None, k=self.order, s=s, t=t, task=task,
            # full_output=full_output)
            self.set_coef(self.results[0])

        else:
            self.results = optimize.leastsq(self.erf, self.coef,
                                            args=(self.x, self.y, self.yerr),
                                            full_output=full_output)
            self.set_coef(self.results[0])


class interfit(curfit):

    """Given an x and y data arrays,  find the best fitting curve.
        After the initial fit, iterate on the solution to reject any
        points which are away from the solution

    * x - list or array of x data
    * y - list or array of y data
    * yerr - error on y data
    * coef - Initial coefficients for fit
    * function - function to be fit to the data:
                 options include polynomial, legendre, chebyshev, or spline
    * order - order of the function that is fit
    * thresh - threshold for rejection
    * niter - number of times to iterate


    """

    def __init__(self, x, y, yerr=None, coef=None, function='poly', order=3,
                 thresh=3, niter=5):
        # set up the variables
        self.x_orig = x
        self.y_orig = y
        self.npts = len(self.x_orig)
        if yerr is None:
            self.yerr_orig = np.ones(self.npts)
        else:
            self.yerr_orig = yerr

        self.order = order
        self.thresh = thresh
        self.niter = niter

        self.set_func(function)
        self.set_coef(coef)
        self.set_mask(init=True)
        self.set_arrays(self.x_orig, self.y_orig, self.mask, err=self.yerr_orig)

    def set_mask(self, init=False):
        """Set the mask according to the values for rejecting points"""
        self.mask = np.ones(self.npts, dtype=bool)
        if init:
            return

        # difference the arrays
        diff = self.y_orig - self(self.x_orig)
        sigma = self.sigma(self.x, self.y)
        self.mask = (abs(diff) < self.thresh * sigma)

    def set_arrays(self, x, y, mask, err=None):
        """set the arrays using a mask"""
        self.x = x[mask]
        self.y = y[mask]
        if err is not None:
            self.yerr = err[mask]

    def interfit(self):
        """Fit a function and then iterate it to reject possible outlyiers"""
        self.fit()
        for i in range(self.niter):
            self.set_mask()
            self.set_arrays(self.x_orig, self.y_orig, self.mask, err=self.yerr_orig)
            self.fit()
