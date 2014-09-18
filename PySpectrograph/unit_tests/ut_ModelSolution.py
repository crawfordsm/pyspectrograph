
import pyfits
import numpy as np
import pylab as pl
from PySpectrograph.Models import RSSModel
from PySpectrograph.WavelengthSolution import ModelSolution as MS

inimage = 'fmbxpP200610180009.fits'

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


def test_ModelSolution():

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

    print instrume, grating, grang, arang, filter
    print xbin, ybin
    print len(data), len(data[0])

    # create the RSS Model
    rssmodel = RSSModel.RSSModel(grating_name=grating, gratang=grang,
                                 camang=arang, slit=slit, xbin=int(xbin),
                                 ybin=int(ybin))

    err = xp * 0.0 + 0.1
    ls = MS.ModelSolution(xp, wp, rssmodel.rss, xlen=len(data[0]), order=4)
    ls.fit(cfit='all')
    ls.fit(ls.ndcoef)
    for c in ls.coef:
        print c(),
    print
    print ls.sigma(ls.x, ls.y)
    print ls.chisq(ls.x, ls.y, err)
    print ls.value(2000)

    pl.figure()
    pl.plot(xp, wp - ls.value(xp), ls='', marker='o')
    #pl.plot(ls.x, ls.y-ls(ls.x), ls='', marker='o')
    pl.show()
    return

test_ModelSolution()
