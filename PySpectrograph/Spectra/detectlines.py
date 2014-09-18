"""detectlines includes tasks and tools for handling 1-d spectra


"""

import numpy as np
from PySpectrograph import SpectrographError

default_kernal = [0, -1, -2, -3, -2, -1, 0, 1, 2, 3, 2, 1, 0]


def centroid(xarr, yarr, kern=default_kernal, mask=None, mode='same'):
    """Find the centroid of a line following a similar algorithm as
        the center1d algorithm in IRAF.   xarr and yarr should be an area
        around the desired feature to be centroided.  The default kernal
        is used if the user does not specific one.

        The algorithm solves for the solution to the equation

        ..math:: \int (I-I_0) f(x-x_0) dx = 0

        returns xc
    """
    if len(yarr) < len(kern):
        raise SpectrographError('Array has to be larger than kernal')

    if mask is not None:
        # catch the fact that at the edges it
        if mask.sum() < len(default_kernal):
            warr = np.convolve(yarr, kern, mode=mode)
            xc = np.interp(0, warr[mask], xarr[mask])
            return xc
        else:
            yarr = yarr[mask]
            xarr = xarr[mask]

    # convle the input array with the default kernal
    warr = np.convolve(yarr, kern, mode=mode)

    # interpolate the results
    xc = np.interp(0, warr, xarr)

    return xc


def detectlines(w_arr, f_arr, sigma=3, bsigma=None, niter=5, mask=None, kern=default_kernal, center=False):
    """Detect lines goes through a 1-D spectra and detect peaks

       w_arr--xaxis array (pixels, wavelength, etc)
       f_arr--yaxis array (flux, counts, etc)
       sigma--Threshold for detecting sources
       bsigma--Threshold for determining background statistics
       niter--iterations to determine background
       center--return centroids and not pixels
       mask--Pixels not to use
    """
    # set up the variables
    if bsigma is None:
        bsigma = sigma

    if mask:
        f_arr = f_arr[mask]
        w_arr = w_arr[mask]

    # find all peaks
    peaks = find_peaks(f_arr, sigma, niter, bsigma=bsigma)

    # set the output values
    xp = w_arr[peaks]

    if center:
        xdiff = int(0.5 * len(kern) + 1)
        x_arr = np.arange(len(w_arr))
        xp = xp * 1.0
        for i in range(len(peaks)):
            cmask = (abs(x_arr - peaks[i]) < xdiff)
            xp[i] = centroid(w_arr, f_arr, kern=kern, mask=cmask)

    return xp


def find_peaks(f_arr, sigma, niter, bsigma=None):
    """Go through an ordered array and find any element which is a peak"""
    # set up the variables
    if bsigma is None:
        bsigma = sigma

    # determine the background statistics
    back_ave, back_std = find_backstats(f_arr, sigma, niter)

    # calculate the differences between the pixels
    dfh = f_arr[1:-1] - f_arr[:-2]
    dfl = f_arr[1:-1] - f_arr[2:]

    # find the objects
    mask = (dfh > 0) * (dfl > 0) * (abs(f_arr[1:-1] - back_ave) > back_std * sigma)
    t = np.where(mask)[0]
    return t + 1


def find_backstats(f_arr, sigma, niter):
    """Iteratively calculate the statistics of an array"""
    ave = f_arr.mean()
    std = f_arr.std()
    for i in range(niter):
        mask = (abs(f_arr - ave) < sigma * std)
        ave = f_arr[mask].mean()
        std = f_arr[mask].std()
    return ave, std
