"""Functions supports the creation of general 1 and n-dim profiles including
normal and sersic profiles

TODO:

Limitations:
1. The sersic profile assumes radial symmetry

"""
import numpy as np


class FunctionsError(Exception):

    """Exception Raised for Function errors"""
    pass


def Normal(arr, mean, sigma, scale=1):
    """Return the normal distribution of N dimensions

       arr--The input array of N dimensions
       mean--the mean of the distribution--either a scalor or N-dimensional array
       sigma--the std deviation of the distribution--either a scalor or N-dimensial array
       scale--the scale of the distrubion

    """
    arr = arr
    mean = mean
    sigma = sigma
    scale = scale
    dim = np.ndim(arr) - 1

    # check to make sure that mean and sigma have dimensions that match dim
    # and create the arrays to calculate the distribution
    if isinstance(mean, np.ndarray):
        if len(mean) != dim:
            raise FunctionsError('Mean and input array are different number of dimensions')
        mean = np.reshape(mean, (dim, 1, 1))
    if isinstance(sigma, np.ndarray):
        if len(sigma) != dim:
            raise FunctionsError('Sigma and input array are different number of dimensions')
        sigma = np.reshape(sigma, (dim, 1, 1))

    # calculate the gaussian
    z = scale * np.exp(-0.5 * (arr - mean) ** 2 / sigma)
    return z


def sersic(arr, deg, r_e, ell=1, scale=1):
    """Produce a light distribution given by a sersic profile
        I=I_e exp(-b [ (R/R_E)^(1/n)-1])

        where:
        Gamma(2n)=gamma(2n, b)
    """

    dim = np.ndim(arr) - 1
    # assume radial symmetry
    if dim == 0:
        r = abs(arr / r_e)
    else:
        pass

    z = scale * exp(-b * (r ** (1 / deg) - 1))
    return z
