import pyfits
import numpy as np
import pylab as pl

from PySpectrograph.WavelengthSolution import LineFit as LF
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


def test_Linefit():

    hdu = pyfits.open(inimage)

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
    alpha = rssmodel.rss.gratang
    beta = rssmodel.rss.gratang - rssmodel.rss.camang

    sigma = 1e7 * rssmodel.rss.calc_resolelement(alpha, beta)

    # create the observed spectrum
    midpoint = int(0.5 * len(data))
    xarr = np.arange(len(data[0]), dtype='float')
    farr = data[midpoint, :]
    obs_spec = Spectrum(xarr, farr, stype='continuum')

    # create artificial spectrum
    stype = 'line'
    w, s = np.loadtxt(inspectra, usecols=(0, 1), unpack=True)
    cal_spec = Spectrum(w, s, wrange=[4000, 5000], dw=0.1, stype=stype, sigma=sigma)
    cal_spec.flux = cal_spec.set_dispersion(sigma=sigma)
    cal_spec.flux = cal_spec.flux * obs_spec.flux.max() / cal_spec.flux.max() + 1

    lf = LF.LineFit(obs_spec, cal_spec, function='legendre', order=3)
    lf.set_coef([4.23180070e+03, 2.45517852e-01, -4.46931562e-06, -2.22067766e-10])
    print(lf(2000))
    print(lf.obs_spec.get_flux(2000), lf.flux(2000))
    print('chisq ', (lf.errfit(lf.coef, xarr, farr) ** 2).sum() / 1e7)
    lf.set_coef([4.23280070e+03, 2.45517852e-01, -4.46931562e-06, -2.22067766e-10])
    print(lf(2000))
    print(lf.obs_spec.get_flux(2000), lf.flux(2000))
    print('chisq ', (lf.errfit(lf.coef, xarr, farr) ** 2).sum() / 1e7)
    # print lf.lfit(xarr)
    # print lf.coef
    # print lf(2000)
    # print lf.results

    pl.figure()

    pl.plot(lf(lf.obs_spec.wavelength), lf.obs_spec.get_flux(xarr))
    pl.plot(lf.cal_spec.wavelength, lf.cal_spec.flux)
    pl.show()


# test_LineSolution()
# test_ModelSolution()
test_Linefit()
