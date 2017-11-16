
"""Unit Test for ImageSolution.  ImageSolution will calculate the
   best spectrograph design given an input spectrograph and data array

"""
import numpy as np
from astropy.io import fits

from ImageSolution import ImageSolution
import RSSModel
import Spectrum

inimage = 'fmbxpP200610180009.fits'
inspectra = 'Xe.dat'


def test_imagesolution():
    # load the image and determine its spectrograph parameters
    hdu = fits.open(inimage)

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

    # create the spectrum
    stype = 'line'
    w, s = np.loadtxt(inspectra, usecols=(0, 1), unpack=True)
    spec = Spectrum.Spectrum(w, s, wrange=[4000, 5000], dw=0.1, stype=stype)

    # Now having the model and the data, set up the variables
    imsol = ImageSolution(data, rssmodel.rss, spec)
    # imsol.fit()
    # for i in imsol.coef:
    #    print imsol.coef[i]()


test_imagesolution()
