"""
SPECTOOLS contains useful functions for handling spectroscopic data


Author                 Version      Date
-----------------------------------------------
S. M. Crawford (SAAO)    1.0        8 Nov 2009

TODO
----


LIMITATIONS
-----------

"""
import pyfits
import numpy as np
from scipy import interpolate as scint
from pyraf import iraf
import saltprint
import saltkey
import saltio
import salttime
import saltstat
import saltsafekey
import saltsafeio
import slottool
from salterror import SaltError


from PySpectrograph import WavelengthSolution
from PySpectrograph.detectlines import detectlines, centroid
from PySpectrograph.Spectrum import Spectrum
from PySpectrograph.apext import apext


class SALTSpecError(SaltError):

    """Errors involving Spec package should cause this exception to be raised."""
    pass


default_kernal = [0, -1, -2, -3, -2, -1, 0, 1, 2, 3, 2, 1, 0]


def mcentroid(xarr, yarr, kern=default_kernal, xc=None, xdiff=None, mode='same'):
    """Find the centroid of a line following a similar algorithm as
       the centroid algorithm in IRAF.   xarr and yarr should be an area
       around the desired feature to be centroided.  The default kernal
       is used if the user does not specific one.

       The algorithm solves for the solution to the equation

       ..math:: \int (I-I_0) f(x-x_0) dx = 0

       returns xc
    """
    if xdiff < len(kern):
        xdiff = len(kern)
    if xc and xdiff:
        mask = (abs(xarr - xc) < xdiff)
    else:
        mask = np.ones(len(xarr))

    return centroid(xarr, yarr, kern=kern, mask=mask, mode=mode)


def interpolate(x, x_arr, y_arr, type='interp', order=3, left=None, right=None):
    """Perform interpolation on value x using arrays x_arr
       and y_arr.  The type of interpolate is defined by interp

       type:
       interp--use numpy.interp
       spline--use scipy.splrep and splev

       return
    """
    if type == 'interp':
        y = np.interp(x, x_arr, y_arr, left=left, right=right)
    if type == 'spline':
        if not left is None:
            y_arr[0] = left
        if not right is None:
            y_arr[-1] = right

        tk = scint.splrep(x_arr, y_arr, k=order)
        y = scint.splev(x, tk, der=0)

    return y


def clipstats(yarr, thresh, iter):
    """Return sigma-clipped mean of yarr"""
    mean = yarr.mean()
    std = yarr.std()
    for i in range(iter):
        mask = (abs(yarr - mean) < thresh * std)
        mean = yarr[mask].mean()
        std = yarr[mask].std()

    return mean, std


def findpoints(xarr, farr, sigma, niter):
    """Find all the peaks and the peak flux in a spectrum"""
    xp = detectlines(xarr, farr, sigma=sigma, niter=niter)
    print(xp)
    mask = [(xp == k).any() for k in xarr]
    xf = np.compress(mask, farr)
    # repeat the second time, but get the centroids for the points
    xp = detectlines(xarr, farr, sigma=sigma, niter=niter, center=True)
    print(xp)
    return xp, xf


def flatspectrum(xarr, yarr, mode='mean', thresh=3, iter=5, order=3):
    """Remove the continuum from a spectrum either by masking it or fitting and subtracting it.
       xarr= input x-vales (pixels or wavelength)
       yarr= flux or counts for the spectrum
       mode=None--no subtraction
            mean--subtract off the mean
            poly--subtact off a fit
            mask--return a spectra with continuum set to zero
    """
    if mode == 'mean':
        # subtract off the mean value
        sarr = yarr - clipstats(yarr, thresh, iter)[0]
    elif mode == 'poly':
        # calculate the statistics and mask all of the mask with values above these
        mean, std = clipstats(yarr, thresh, iter)
        mask = (yarr < mean + thresh * std)
        coef = np.polyfit(xarr[mask], yarr[mask], order)
        sarr = yarr - np.polyval(coef, xarr)
    elif mode == 'mask':
        # mask the values
        mean, std = clipstats(yarr, thresh, iter)
        mask = (yarr < mean + thresh * std)
        sarr = yarr.copy()
        sarr[mask] = 0
    else:
        sarr = yarr.copy()
    return sarr


def findwavelengthsolution(xarr, farr, sl, sf, ws, sigma=5, niter=5):
    """Calculates the wavelength solution given a spectra and a set of lines.  Hopefully
       an accurate first guess (ws) is provided and relative fluxes are provided as well,
       but if not, then the program is still designed to attempt to handle it.

       returns ws
    """

    # detect lines in the input spectrum and identify the peaks and peak values
    xp, xf = findpoints(xarr, farr, sigma, niter)

    # return no solution if no peaks were found
    if len(xp) == 0:
        return None

    # find the best match to the lines
    try:
        wp = findmatch(xarr, farr, xp, xf, sl, sf, ws)
        for i in range(len(xp)):
            if wp[i] > -1:
                print(xp[i], wp[i])
    except Exception as e:
        message = 'Unable to match line lists because %s' % e
        raise SALTSpecError(message)

    # if wavelength solution does not already exist, create it3dd:
    if not isinstance(ws, WavelengthSolution.WavelengthSolution):
        message = 'Wavelength solution does not exist'
        raise SALTSpecError(message)

    # find the solution to the best fit
    mask = (wp > 0)
    if mask.sum():
        nws = WavelengthSolution.WavelengthSolution(xp[mask], wp[mask], order=ws.order, function=ws.function)
        nws.fit()
    else:
        nws = None

    return nws


def findfeatures(xarr, farr, sl, sf, ws, sigma=5, niter=5):
    """Given a spectra, detect lines in the spectra, and find lines in
       the line list that correspond to those lines
    """

    # detect lines in the input spectrum and identify the peaks and peak values
    xp, xf = findpoints(xarr, farr, sigma, niter)

    # return no solution if no peaks were found
    if len(xp) == 0:
        return None

    # find the best match to the lines
    try:
        wp = findmatch(xarr, farr, xp, xf, sl, sf, ws)
        for i in range(len(xp)):
            if wp[i] > -1:
                print(xp[i], wp[i])
    except Exception as e:
        message = 'Unable to match line lists because %s' % e
        raise SALTSpecError(message)
    return xp, wp


def findmatch(xarr, farr, xp, xf, sl, sf, ws, xlimit=5, wlimit=2):
    """Find the best match between the observed arc lines and the spectral line list.  If available,
       use the line fluxes and the wavelength solution.  Returns a an array that is a wavelength
       for each peak wavelength

       returns wp
    """
    wp = xp * 0.0 - 1

    # calculate it using only xp and sl
    if sf is None and not ws:
        print('no')

    # calculate it without any wavelength solution
    elif not ws:
        print(ws)

    # calculate it without any flux information
    elif sf is None and ws:
        for i in xf.argsort()[::-1]:
            cx = mcentroid(xarr, farr, xc=xp[i], xdiff=4)
            if abs(cx - xp[i]) < xlimit:
                w = wavematch(ws.value(cx), wp, sl)
                wp[i] = w

    # calculate it using all of the information
    else:
        dcoef = ws.coef * 0.0
        dcoef[-1] = 10
        dcoef[-2] = dcoef[-2] * 0.2
        nstep = 20
        nws = spectramatch(xarr, farr, sl, sf, ws, dcoef, nstep=nstep, res=2, dres=0.1)
        print('nws:', nws.coef)
        for i in xf.argsort()[::-1]:
            cx = mcentroid(xarr, farr, xc=xp[i], xdiff=4)
            if abs(cx - xp[i]) < xlimit:
                w = wavematch(nws.value(cx), wp, sl, wlimit=4 * dcoef[-1] / nstep)
                wp[i] = w

    return wp


def spectramatch(xarr, farr, sw, sf, ws, dcoef, nstep, res=2, dres=0.1):
    """Using all the information which is available, cross correlate the observed spectra
       and the wavelength spectra to find the best coefficients and match the data
    """
    # create an artificial spectrum of the lines
    lmax = farr.max()
    swarr, sfarr = makeartificial(sw, sf, lmax, res, dres)

    # Now find the best fitting coefficients for the wavelength solution
    nws = WavelengthSolution.WavelengthSolution(ws.x_arr, ws.w_arr)
    nws.coef = ws.coef

    # create the range of coefficents
    dlist = mod_coef(ws.coef, dcoef, 0, nstep)

    # loop through them and deteremine the best cofficient
    cc_arr = np.zeros(len(dlist), dtype=float)
    for i in range(len(dlist)):
        # set the coeficient
        nws.coef = dlist[i]

        # set the wavelegnth coverage
        warr = nws.value(xarr)

        # resample the artificial spectrum at the same wavelengths as the
        asfarr = np.interp(warr, swarr, sfarr, left=0.0, right=0.0)

        # calculate the correlation value
        cc_arr[i] = ncor(farr, asfarr)

    nws.coef = dlist[cc_arr.argmax()]

    # return the best fit solution
    return nws


def mod_coef(coef, dcoef, index, nstep):
    """For a given index, return
    """
    dlist = []
    if index >= len(coef):
        return dlist
    if dcoef[index] == 0:
        if index < len(coef):
            dlist.extend((mod_coef(coef, dcoef, index + 1, nstep)))
        else:
            dlist.append(coef)
        return dlist

    if index < len(coef) - 1:
        for x in np.arange(-dcoef[index], dcoef[index], 2 * dcoef[index] / nstep):
            ncoef = coef.copy()
            ncoef[index] = coef[index] + x
            dlist.extend(mod_coef(ncoef, dcoef, index + 1, nstep))
    else:
        for x in np.arange(-dcoef[index], dcoef[index], 2 * dcoef[index] / nstep):
            ncoef = coef.copy()
            ncoef[index] = coef[index] + x
            dlist.append(ncoef)
    return dlist


def makeartificial(sw, sf, fmax, res, dw, pad=10, nkern=200, wrange=None):
    """For a given line list with fluxes, create an artifical spectrum"""
    if wrange is None:
        wrange = [sw.min() - pad, sw.max() + pad]
    spec = Spectrum(sw, sf, wrange=wrange, dw=dw, stype='line', sigma=res)
    spec.flux = spec.flux * fmax / spec.flux.max()

    return spec.wavelength, spec.flux


def ncor(x, y):
    """Calculate the normalized correlation of two arrays"""
    return np.correlate(x, y, old_behavior=False) / \
        (np.correlate(x, x, old_behavior=False) * np.correlate(y, y, old_behavior=False)) ** 0.5


def wavematch(w, wp, sl, wlimit=10):
    """Compare a wavelength to an observed list and see if it matches up.  Skip
       if the lines is already in the wp list

    """

    # first remove anything already in the self.wp from the sl list
    lines = []
    for x in sl:
        if x not in wp:
            lines.append(x)
    lines = np.array(lines)

    # find the best match
    dist = abs(lines - w)
    if dist.min() < wlimit:
        i = dist.argmin()
    else:
        return -1

    # return the values
    return lines[i]


def findfit(xp, wp, order=3, function='poly'):
    """Find the fit using just the matched points of xp and wp"""
    ws = WavelengthSolution.WavelengthSolution(xp, wp, order=order, function=function)
    ws.fit()
    return ws


def findzeropoint(xarr, farr, swarr, sfarr, ws, dc=10, nstep=20, inttype='interp'):
    """Uses cross-correlation to find the best fitting zeropoint"""

    # if an initial solution, then cut the template lines to just be the length of the spectrum
    if ws is None:
        return ws

    # set up the the dc coefficient
    dcoef = ws.coef * 0.0
    dcoef[-1] = dc

    ws = findxcor(xarr, farr, swarr, sfarr, ws, dcoef=dcoef, nstep=nstep, inttype=inttype)
    return ws


def findxcor(xarr, farr, swarr, sfarr, ws, dcoef=None, nstep=20, inttype='interp'):
    """Find the solution using crosscorrelation of the wavelength solution.  An initial
       guess needs to be supplied along with the variation in each coefficient and the
       number of steps to calculate the correlation.  The input wavelength and flux
       for the known spectral features should be in the format where they have already
       been convolved with the response function of the spectrograph

       xarr--Pixel coordinates of the image

       farr--Flux values for each pixel

       swarr--Input wavelengths of known spectral features

       sfarr--fluxes of known spectral features

       ws--current wavelength solution

       dcoef--Variation over each coefficient for correlation

       nstep--number of steps to sample over

       inttype--type of interpolation

    """

    # cross-correlate the spectral lines and the observed fluxes in order to refine the solution
    nws = WavelengthSolution.WavelengthSolution(ws.x_arr, ws.w_arr, order=ws.order, function=ws.function)
    nws.setcoef(ws.coef)

    # create the range of coefficents
    if dcoef is None:
        dcoef = ws.coef * 0.0 + 1.0
    dlist = mod_coef(ws.coef, dcoef, 0, nstep)

    # loop through them and deteremine the best cofficient
    cc_arr = np.zeros(len(dlist), dtype=float)
    zp_arr = np.zeros(len(dlist), dtype=float)
    dp_arr = np.zeros(len(dlist), dtype=float)
    for i in range(len(dlist)):
        # set the coeficient
        nws.setcoef(dlist[i])

        # set the wavelegnth coverage
        warr = nws.value(xarr)

        # resample the artificial spectrum at the same wavelengths as the
        asfarr = interpolate(warr, swarr, sfarr, type=inttype, left=0.0, right=0.0)

        # calculate the correlation value
        cc_arr[i] = ncor(farr, asfarr)

    # now set the best coefficients
    i = cc_arr.argmax()
    bcoef = dlist[i]
    nws.setcoef(bcoef)
    darr = np.array(dlist)
    for j in range(len(nws.coef)):
        if dcoef[j] != 0.0:
            tk = np.polyfit(darr[:, j], cc_arr, 2)
            bval = -0.5 * tk[1] / tk[0]
            if abs(bval - bcoef[j]) < dcoef[j]:
                bcoef[j] = bval

            #coef=np.polyfit(dlist[:][j], cc_arr, 2)
            # nws.coef[j]=-0.5*coef[1]/coef[0]

    return nws

#------------------------------------------------------------------
# Read in the line list file


def readlinelist(linelist):
    """Read in the line lists.  Determine what type of file it is.  The default is
        an ascii file with line and relative intensity.  The other types are just line,
        or a wavelenght calibrated fits file

       return lines, fluxes, and status
    """
    slines = []
    sfluxes = []
    status = 0

    # Check to see if it is a fits file
    # if not, then read in the ascii file
    if linelist[-4:] == 'fits':
        try:
            slines, sfluxes = readfitslinelist(linelist)
        except Exception as e:
            message = 'Unable to read in the line list %s because %s' % (linelist, e)
            raise SALTSpecError(message)
    else:
        try:
            slines, sfluxes = readasciilinelist(linelist)
        except Exception as e:
            message = 'Unable to read in the line list %s because %s' % (linelist, e)
            raise SALTSpecError(message)

    # conver to numpy arrays
    try:
        slines = np.asarray(slines)
        sfluxes = np.asarray(sfluxes)
    except Exception as e:
        message = 'Unable to create numpy arrays because %s' % (e)
        raise SALTSpecError(logfile, message)

    return slines, sfluxes

#------------------------------------------------------------------
# Read in the line list file


def readfitslinelist(linelist):
    """Read in the line lists from an fits file.  If it is a 2-D array
       it will assume that it is an image and select the central wavlength

       return lines, fluxes, and status
    """
    slines = []
    sfluxes = []

    # open the image
    shdu = pyfits.open(linelist)
    nhdu = len(shdu)
    # determine if it is a one or two-d image
    # if ndhu=0 then assume that it is in the zeroth image
    # otherwise assume the data is in the first extension
    # assumes the x-axis is the wavelength axis
    if nhdu == 1:
        ctype1 = shdu[0].header['CTYPE1']
        crval1 = shdu[0].header['CRVAL1']
        cdelt1 = shdu[0].header['CDELT1']
        if shdu[0].data.ndim == 1:
            data = shdu[0].data
            wave = crval1 + cdelt1 * np.arange(len(shdu[0].data))

    # detect lines in the input spectrum and identify the peaks and peak values
    slines, sfluxes = findpoints(wave, data, 3, 5)
    """
   figure(figsize=(8,8), dpi=72)
   axes([0.1, 0.1, 0.8, 0.8])
   plot(wave, data, ls='-')
   plot(slines, sfluxes, ls='', marker='o')
   xlim(4220,4900)
   show()
   """

    return slines, sfluxes

#------------------------------------------------------------------
# Read in the line list file


def readasciilinelist(linelist):
    """Read in the line lists from an ascii file.  It can either be a
        file with one or two columns.  Only read in lines that are not
        commented out.

       return lines, fluxes, and status
    """
    slines = []
    sfluxes = []

    # read in the file
    f = open(linelist)
    lines = f.readlines()
    f.close()

    # for each line,
    for l in lines:
        l = l.strip()
        if not (l and l.startswith('#')):
            l = l.split()
            slines.append(float(l[0]))
            try:
                sfluxes.append(float(l[1]))
            except IndexError:
                sfluxes.append(-1)
    return slines, sfluxes
