import pylab as pl

from PySpectrograph.Models import RSSModel
from PySpectrograph.Spectra import Spectrum

# create the spectrograph model
rss = RSSModel.RSSModel(grating_name="PG0900", gratang=15.875, camang=31.76496,
                        slit=1.50, xbin=2, ybin=2)


# print out some basic statistics
print(1e7 * rss.calc_bluewavelength(), 1e7 * rss.calc_centralwavelength(), 1e7 * rss.calc_redwavelength())
R = rss.calc_resolution(rss.calc_centralwavelength(), rss.alpha(), -rss.beta())
res = 1e7 * rss.calc_resolelement(rss.alpha(), -rss.beta())
print(R, res)

# set up the detector
ycen = rss.detector.get_ypixcenter()
d_arr = rss.detector.make_detector()[ycen, :]
w = 1e7 * rss.get_wavelength(d_arr)

# set up the artificial spectrum
sw, sf = pl.loadtxt('Ne.txt', usecols=(0, 1), unpack=True)
wrange = [1e7 * rss.calc_bluewavelength(), 1e7 * rss.calc_redwavelength()]
spec = Spectrum.Spectrum(sw, sf, wrange=wrange, dw=res / 10, stype='line', sigma=res)

# interpolate it over the same range as the detector
spec.interp(w)


# plot it
pl.figure()
pl.plot(spec.wavelength, d_arr * ((spec.flux) / spec.flux.max()))
pl.show()
