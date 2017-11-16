import numpy as np
import pylab as pl

from astropy.io import fits

from PySpectrograph.WavelengthSolution import WavelengthSolution as WS
from PySpectrograph.Models import RSSModel
from PySpectrograph.Spectra import Spectrum

inimage = 'fmbxpP200610180009.fits'
inspectra = 'Xe.dat'


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
ep = xp * 0.0 + 0.1


def test_LineSolution():
    ws = WS.WavelengthSolution(xp, wp, function='poly')
    ws.fit()
    print(ws.func.coef)
    print(ws.value(2000))
    print(ws.sigma(ws.func.x, ws.func.y))
    print(ws.chisq(ws.func.x, ws.func.y, ws.func.yerr))

    pl.figure()
    pl.plot(xp, wp - ws.value(xp), ls='', marker='o')
    # pl.plot(ls.x, ls.y-ls(ls.x), ls='', marker='o')
    pl.show()
    return


def test_ModelSolution():

    hdu = fits.open(inimage)

    # create the data arra
    data = hdu[1].data

    # create the header information
    grating = hdu[1].header['GRATING'].strip()
    grang = hdu[1].header['GR-ANGLE']
    arang = hdu[1].header['AR-ANGLE']
    slit = float(hdu[1].header['MASKID'])
    xbin, ybin = hdu[1].header['CCDSUM'].strip().split()

    # print instrume, grating, grang, arang, filter
    # print xbin, ybin
    # print len(data), len(data[0])

    # create the RSS Model
    rssmodel = RSSModel.RSSModel(grating_name=grating, gratang=grang,
                                 camang=arang, slit=slit, xbin=int(xbin),
                                 ybin=int(ybin))

    xarr = np.arange(len(data[0]), dtype='int64')
    rss = rssmodel.rss
    alpha = rss.gratang
    beta = rss.gratang - rss.camang
    d = rss.detector.xbin * rss.detector.pix_size * (xarr - 0.5 * len(xarr))
    dbeta = np.degrees(np.arctan(d / rss.camera.focallength))
    y = 1e7 * rss.calc_wavelength(alpha, beta - dbeta)

    # ws=WS.WavelengthSolution(xp, wp, function='model', sgraph=rssmodel.rss, xlen=len(data[0]), order=4, cfit='all')
    ws = WS.WavelengthSolution(
        xarr,
        y,
        function='model',
        sgraph=rssmodel.rss,
        xlen=len(
            data[0]),
        order=4,
        cfit='ndcoef')

    # ws=WS.WavelengthSolution(xarr, y, function='poly', order=3)
    # ws=WS.WavelengthSolution(xp, wp, function='poly', order=3)

    ws.fit()
    print(ws.coef)
    print(ws.func.result)
    for c in ws.func.spcoef:
        print(c())
    for c in ws.func.ndcoef:
        print(c())
    print(ws.value(2000))
    print(ws.sigma(xp, wp))
    print(ws.chisq(xp, wp, ep))

    pl.figure()
    # pl.plot(xarr,y)
    # pl.plot(xarr, ws.value(xarr))
    # pl.plot(xp, wp-ws.value(xp), ls='', marker='o')
    pl.plot(xarr, y - ws.value(xarr))
    # pl.plot(ls.x, ls.y-ls(ls.x), ls='', marker='o')
    pl.show()


def test_Linefit():

    hdu = fits.open(inimage)

    # create the header information
    grating = hdu[1].header['GRATING'].strip()
    grang = hdu[1].header['GR-ANGLE']
    arang = hdu[1].header['AR-ANGLE']
    slit = float(hdu[1].header['MASKID'])
    xbin, ybin = hdu[1].header['CCDSUM'].strip().split()

    # print instrume, grating, grang, arang, filter
    # print xbin, ybin
    # print len(data), len(data[0])

    # create the RSS Model
    rssmodel = RSSModel.RSSModel(grating_name=grating, gratang=grang,
                                 camang=arang, slit=slit, xbin=int(xbin),
                                 ybin=int(ybin))
    alpha = rssmodel.rss.gratang
    beta = rssmodel.rss.gratang - rssmodel.rss.camang

    sigma = 1e7 * rssmodel.rss.calc_resolelement(alpha, beta)

    # create artificial spectrum
    # create the spectrum
    stype = 'line'
    w, s = np.loadtxt(inspectra, usecols=(0, 1), unpack=True)
    spec = Spectrum(w, s, wrange=[4000, 5000], dw=0.1, stype=stype, sigma=sigma)
    spec.flux = spec.set_dispersion(sigma=sigma)


# test_LineSolution()
# test_ModelSolution()
test_Linefit()
