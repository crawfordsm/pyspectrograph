
"""Unit Test for LineSolution.  LineSolution will calculate the
   best spectrograph design given an input spectrograph and spectrum

"""
import pyfits
import numpy as np

from LineSolution import LineSolution
import RSSModel
import Spectrum

inimage = 'fmbxpP200610180009.fits'
inspectra = 'Xe.dat'

import pylab as pl

xp = np.array([1134.87, 1239.11, 1498.22, 1687.21, 1904.49, 1997.15, 2025.77, 2202.59, 2559.72, 2673.74, 3124.34])
wp = np.array([4500.9772,
               4524.6805,
               4582.7474,
               4624.2757,
               4671.226,
               4690.9711,
               4697.02,
               4734.1524,
               4807.019,
               4829.709,
               4916.51])


def test_imagesolution():
    # load the image and determine its spectrograph parameters
    hdu = pyfits.open(inimage)

    # create the data arra
    data = hdu[1].data

    # create the header information
    instrume = hdu[1].header['INSTRUME'].strip()
    grating = hdu[1].header['GRATING'].strip()
    grang = hdu[1].header['GR-ANGLE']
    arang = hdu[1].header['AR-ANGLE']
    filter = hdu[1].header['FILTER'].strip()
    slit = float(hdu[1].header['MASKID'])
    xbin, ybin = hdu[1].header['CCDSUM'].strip().split()

    print(instrume, grating, grang, arang, filter)
    print(xbin, ybin)

    # create the RSS Model
    rssmodel = RSSModel.RSSModel(grating_name=grating, gratang=grang,
                                 camang=arang, slit=slit, xbin=int(xbin),
                                 ybin=int(ybin))
    alpha = rssmodel.rss.gratang
    beta = rssmodel.rss.gratang - rssmodel.rss.camang

    sigma = 1e7 * rssmodel.rss.calc_resolelement(alpha, beta)

    # create the spectrum
    stype = 'line'
    w, s = np.loadtxt(inspectra, usecols=(0, 1), unpack=True)
    spec = Spectrum.Spectrum(w, s, wrange=[4000, 5000], dw=0.1, stype='line', sigma=sigma)
    # spec.flux=spec.set_dispersion(sigma=sigma)

    # Now having the model and the data, set up the variables
    j = int(len(data) / 2)
    xlen = len(data[0])
    xarr = np.arange(len(data[0]))
    farr = data[j, :]
    var = abs(farr) + farr.mean()
    var = var * (farr > 1000) + 1

    if 1:
        imsol = LineSolution(rssmodel.rss, xarr=xarr, farr=farr, spectrum=spec, yval=0, var=var, order=2)
        output = imsol.fit(imsol.makeflux, imsol.coef, imsol.xarr, imsol.farr, imsol.var)
        #imsol.fit(imsol.makeflux, imsol.ndcoef, imsol.xarr, imsol.farr, imsol.var)
        #imsol.fit(imsol.makeflux, imsol.coef, imsol.xarr, imsol.farr, imsol.var)
        # for i in range(len(imsol.coef)):
        #    print imsol.coef[i]()
        # for i in range(len(imsol.ndcoef)):
        #    print imsol.ndcoef[i]()

        print(output.beta)
        print(imsol.value(output.beta, 500))

        # check the results
        warr = imsol.value(output.beta, imsol.xarr)
        print((wp - imsol.value(output.beta, xp)).mean(), (wp - imsol.value(output.beta, xp)).std())
        pl.figure()
        pl.plot(imsol.spectrum.wavelength, imsol.spectrum.flux * imsol.farr.max() / imsol.spectrum.flux.max())
        pl.plot(warr, imsol.farr)
        #pl.plot(xp, wp-imsol.value(xp), ls='', marker='o')
        pl.show()

    # okay now test the results for a purely matched lines
    fp = xp * 0.0 + 2.0
    var = xp * 0.0 + 1.0
    xspec = Spectrum.Spectrum(xp, fp, dw=0.1, stype='line', sigma=sigma)
    wspec = Spectrum.Spectrum(wp, fp, wrange=[4000, 5000], dw=0.1, stype='line', sigma=sigma)

    imsol = LineSolution(rssmodel.rss, xarr=xspec.wavelength, farr=xspec.flux, spectrum=wspec, yval=0, var=var, order=3)
    imsol.xlen = xlen
    output = imsol.fit(imsol.value, imsol.coef, xp, wp, var)
    print(output.beta)
    #imsol.fit(imsol.value, imsol.ndcoef, xp, wp, var)
    # for i in range(len(imsol.coef)):
    # print imsol.coef[i])
    # for i in range(len(imsol.ndcoef)):
    #    print imsol.ndcoef[i]()

    print(imsol.value(output.beta, 500))
    print((wp - imsol.value(output.beta, xp)).mean(), (wp - imsol.value(output.beta, xp)).std())

    # check the results
    warr = imsol.value(output.beta, xarr)
    pl.figure()
    pl.plot(spec.wavelength, spec.flux * farr.max() / spec.flux.max())
    pl.plot(warr, farr)
    #pl.plot(xp, wp-imsol.value(output.beta, xp), ls='', marker='o')
    pl.show()


test_imagesolution()
