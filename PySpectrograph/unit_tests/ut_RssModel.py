
"""Unit Test for LineSolution.  LineSolution will calculate the
   best spectrograph design given an input spectrograph and spectrum

"""
import pyfits
import numpy as np

from PySpectrograph.Models import RSSModel

inimage = 'fmbxpP200610180009.fits'
inspectra = 'Xe.dat'


def test_rssmodel():
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

    # create the RSS Model
    rss = RSSModel.RSSModel(grating_name=grating, gratang=grang,
                            camang=arang, slit=slit, xbin=int(xbin),
                            ybin=int(ybin))
    alpha = rss.alpha()
    beta = rss.beta()

    sigma = 1e7 * rss.calc_resolelement(alpha, beta)
    print "SIGMA: ", sigma

    # test to see if giving the wrong name it will raise an error
    try:
        rss = RSSModel.RSSModel(grating_name="not a grating", gratang=grang,
                                camang=arang, slit=slit, xbin=int(xbin),
                                ybin=int(ybin))
    except RSSModel.RSSError as e:
        pass


test_rssmodel()
