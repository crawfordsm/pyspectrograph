"""Spectrum is a class to describe and generate a spectrum.

HISTORY
20100601 SMC  First written by SM Crawford

Limitations:

"""

import numpy as np
from PySpectrograph.Utilities.Functions import Normal


class SpectrumError(Exception):

    """Exception Raised for Spectrograph errors"""
    pass


class Spectrum:

    """Spectrum is a class for handling and creating spectra.  It can either
       be created from a line list or a continuum flux.  If no inputs are given,
       it creates an empty object that has a wavelength range from 3000-9000 and
       a sampling of 0.1.

       Parameters
       ----------
       wavelength:  An array of wavelengths for values in flux

       flux:  An array of fluxes.  It can either be lines or continuum

       wrange:  wavelenght range of flux

       dw: sampling of the spectrum

       stype: line--the input spectrum is a list of lines
              continuum--the input spectrum is continuum values

       wavelength_units:  Units of the wavelength array

       flux_units:  Units of the wavelength array
    """

    def __init__(self, wavelength=None, flux=None, var=None, wrange=None, dw=0.1,
                 sigma=1e-1, stype='line', wavelength_units=None, flux_units=None):

        if wavelength is None and wrange is None:
            raise SpectrumError('Please specify either wavelength or wrange')

        # set the variables
        self.wrange = wrange
        self.dw = dw
        self.wavelength_units = wavelength_units
        self.stype = stype
        self.flux_units = flux_units
        self.sigma = sigma
        self.var = var

        # set the wavelength
        self.set_wavelength(wavelength)

        # set the flux
        self.set_flux(wavelength, flux)

        return

    def set_wavelength(self, wavelength, new=False):
        """Set the wavelength scale

           If new is True, then it will create a new wavelength
           using wrange and dw

        """

        if self.wrange is None and wavelength is not None:
            self.wrange = [wavelength.min(), wavelength.max()]

        if self.stype == 'line' or new:
            self.wavelength = np.arange(self.wrange[0], self.wrange[1], self.dw)
        elif self.stype == 'continuum':
            self.wavelength = wavelength
        else:
            self.wavelength = wavelength

        self.nwave = len(self.wavelength)

    def set_flux(self, wavelength, flux):
        """Set the flux levels"""

        if flux is not None and len(flux) == self.nwave:
            self.flux = flux

        elif flux is not None and wavelength is not None:
            if self.stype == 'line':
                self.flux = np.zeros(self.nwave, dtype=float)
                for w, f in zip(wavelength, flux):
                    self.flux += Normal(self.wavelength, w, self.sigma * self.dw, f)
            elif self.stype == 'continuum':
                self.flux = np.interp(self.wavelength, wavelength, flux)
            else:
                raise SpectrumError('%s is not an acceptable stype option' % stype)
        else:
            self.flux = np.zeros(self.nwave, dtype=float)

    def set_dispersion(self, sigma=1.0, nkern=20, ktype='Gaussian', func=None):
        """Applies an dispersion to the spectrum.

           sigma:  Dispersion to apply if assuming a Gaussian in units of wavelength

           nkern:  sampling of kernal

           ktype:  'Gaussian'--Apply a Gaussian function
                   'User'  -- Use a user supplied kernal

           func:   User supplied kernal

        """
        if ktype == 'Gaussian':
            xkern = np.arange(nkern)
            kern = np.exp(-(xkern - 0.5 * nkern) ** 2 / (sigma / self.dw) ** 2)
        elif ktype == 'User':
            kern = func
        else:
            raise SpectrumError('%s is not an acceptable kernal option' % ktype)

        return np.convolve(self.flux, kern, mode='same')

    def get_flux(self, w):
        """Given a wavelength w, return what the flux value is"""
        return np.interp(w, self.wavelength, self.flux)

    def interp(self, warr):
        """Re-interpolate the spectrum such that the wavelength sampling is given by warr"""
        self.flux = np.interp(warr, self.wavelength, self.flux)
        self.wavelength = warr

"""
TODO:
Here's the full text from Morton 1991 ApJS 77, 119 and it is a
good question whether we should adopt the data from Peck and
Reeder as IR data will be important for us.

The IAU standard for conversion between air and vacuum wavelengths is,
according to Oosterhoff (1957) and Edlen (1953),

\begin{equation}
\frac{\lambda_{vac}-\lambda_{air}}{\lambda_{air}}=(n-1) =
   6.4328\times10^{-5}+ \frac{2.94981 times10^{-2}}{146-\sigma^2}
   +\frac{2.5540\times10^{-4}}{41-\sigma^2}
\end{equation}
where $\sigma=10^4/$, wiht $\lambda$ in angstroms.  Edlen (1966) and Peck
& Reeder (1972) have proposed improvmeents that primarily affect infrared
wavelengths, but equation (3) was used here.

More recently, this has been updated in IDLASTRO with a
formula from Ciddor 1996, Applied Optics 62, 958
See http://idlastro.gsfc.nasa.gov/ftp/pro/astro/airtovac.pro
"""


def air2vac(w_air, mode='Morton'):
    """Given an wavelength in units of

    w--wavelength in air in units of angstrom

    mode--method to use for conversion
       Morton--Morton 1991 ApJS 77, 119
       Ciddor--Ciddor 1996, Applied Optics 62, 958

    """
    if mode == 'Morton':
        sigmasq = (1e4 / w_air) ** 2
        w_vac = w_air * (1 + 6.4328e-5 + 2.94981e-2 / (146.0 - sigmasq) + 2.5540e-4 / (41.0 - sigmasq))
    elif mode == 'Ciddor':
        sigmasq = (1e4 / w_air) ** 2
        w_vac = w_air * (1 + 5.792105e-2 / (238.0185 - sigmasq) + 1.67917e-3 / (57.362 - sigmasq))
    else:
        raise SpectrumError('%s is an invalid mode' % mode)

    return w_vac


def vac2air(w_vac, mode='Morton'):
    """Given an wavelength in units of

    w--wavelength in air in units of angstrom

    mode--method to use for conversion
       Morton--Morton 1991 ApJS 77, 119
       Ciddor--Ciddor 1996, Applied Optics 62, 958

    """
    if mode == 'Morton':
        sigmasq = (1e4 / w_vac) ** 2
        w_air = w_vac / (1 + 6.4328e-5 + 2.94981e-2 / (146.0 - sigmasq) + 2.5540e-4 / (41.0 - sigmasq))
    elif mode == 'Ciddor':
        sigmasq = (1e4 / w_vac) ** 2
        w_air = w_vac / (1 + 5.792105e-2 / (238.0185 - sigmasq) + 1.67917e-3 / (57.362 - sigmasq))
    else:
        raise SpectrumError('%s is an invalid mode' % mode)

    return w_air


def fnutofwave(warr, farr):
    """Converts farr in ergs/s/cm2/Hz to ergs/s/cm2/A"""
    c = 2.99792458e18  # spped of light in Angstroms/s
    return farr * c / warr ** 2


def magtoflux(marr, fzero):
    """Convert from magnitude to flux.
     marr--input array in mags
     fzero--zero point for the conversion
    """
    return fzero * 10 ** (-0.4 * marr)


def fluxtomag(farr, fzero):
    """"Convert from flux to magnitudes
        farr--input array in flux units
        fzero--zero point for the converion
    """
    return -2.5 * np.log10(farr / fzero)


if __name__ == '__main__':
    import sys
    from pylab import *

    infile = sys.argv[1]
    stype = sys.argv[2]
    w1 = float(sys.argv[3])
    w2 = float(sys.argv[4])
    dw = float(sys.argv[5])

    w, s = np.loadtxt(infile, usecols=(0, 1), unpack=True)

    spec = Spectrum(w, s, wrange=[w1, w2], dw=dw, stype=stype)
    spec.set_dispersion(2, nkern=200)

    figure()
    plot(spec.wavelength, spec.flux)
    show()
