"""LineSolution is a task describing the functional form for transforming
pixel position to wavelength.  The inputs for this task are the given pixel position
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
from PySpectrograph.Utilities.fit import interfit


class LineSolution(interfit):

    """LineSolution is a task describing the functional form for transforming
       pixel position to wavelength.

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

    def value(self, x):
        """Calculate the value of the array
        """
        return self(x)
