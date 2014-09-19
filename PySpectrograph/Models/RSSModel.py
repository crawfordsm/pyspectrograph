"""RSSmodel is a class that describes the optical arm of the Robert Stobie
Spectrograph.  This model inherents from Spectragraph.  The RSS is currently
described by the grating, grating angle, and camera angle.  All other components
are fixed.

20090913  SMC   First version
20120325  SMC   -Update to include deviations in the as built spectrograph in the grating angle
                -include the most up to date chip geometry
20120511  SMC   -Updated to inherit directly from Spectrograph
                -Included get_wavelength

Limitations:
-Set the geometry of the CCDs using the current geometry file instead of the
 default given here

"""

import numpy as np

from PySpectrograph.Spectrograph import Spectrograph, Grating, Optics, CCD, Detector, Slit


class RSSError(Exception):
    pass


class RSSModel (Spectrograph):

    """A model describing the RSS spectrograph"""

    def __init__(self, grating_name='None', gratang=45, camang=45, slit=1.0,
                 xbin=2, ybin=2, xpos=-0.2666, ypos=0.0117, wavelength=None):

        # set up the parts of the grating
        self.grating_name = grating_name
        self.slitang = slit

        # set the telescope
        self.set_telescope('RSS')

        # set the collimator
        self.set_collimator('RSS')

        # set the camera
        self.set_camera('RSS', wavelength=wavelength)

        # set the detector
        self.set_detector('RSS', xbin=xbin, ybin=ybin, xpos=xpos, ypos=ypos)

        # set up the grating
        self.set_grating(self.grating_name)

        # set up the slit
        self.set_slit(self.slitang)

        # set up the grating angle
        self.gratang = gratang
        self.camang = camang

    def alpha(self, da=0.23):
        """Return the value of alpha for the spectrograph"""
        return self.gratang + da

    def beta(self, mF=4.2e-5, db=0.00):
        """Return the value of beta for the spectrograph

           Beta_o=(1+fA)*(camang)-gratang+beta_o
        """
        return (1 + mF) * self.camang - self.alpha() + db

    def focallength(self, focallength=328.0, wavelength=None):
        """The camera focal lenght set oaccording to the following
           Wavelength should be provided in angstroms

        """
        if wavelength is None:
            return focallength

        L = (wavelength - 4000.0) / 1000.0
        return 327.66 + -0.1861 * L + 0.5061 * L ** 2 + -0.2100 * L ** 3 + 0.0365 * L ** 4 + -0.0023 * L ** 5

    def get_wavelength(self, xarr, gamma=0.0):
        """For a given spectrograph configuration, return the wavelength coordinate
           associated with a pixel coordinate.

           xarr: 1-D Array of pixel coordinates
           gamma: Value of gamma for the row being analyzed

           returns an array of wavelengths in mm
        """
        d = self.detector.xbin * self.detector.pix_size * (xarr - self.detector.get_xpixcenter())
        dbeta = -np.degrees(np.arctan(d / self.camera.focallength))
        return self.calc_wavelength(self.alpha(), -self.beta() + dbeta, gamma=gamma)

    def set_telescope(self, name='RSS'):
        if name == 'RSS':
            self.telescope = Optics(name=name, focallength=46200.0)
        elif name == 'SALT':
            self.telescope = Optics(name=name, focallength=46200.0)
        else:
            raise RSSError('%s is not a supported Telescope' % name)

    def set_collimator(self, name='RSS', focallength=630.0):
        if name == 'RSS':
            self.collimator = Optics(name=name, focallength=focallength)
        else:
            raise RSSError('%s is not a supported collimator' % name)

    def set_camera(self, name='RSS', wavelength=None):
        if name == 'RSS':
            self.camera = Optics(name=name, focallength=self.focallength())
        else:
            raise RSSError('%s is not a supported camera' % name)

    def set_detector(self, name='RSS', geom=None, xbin=2, ybin=2, xpos=0, ypos=0):
        if name == 'RSS':
            if geom:
                pass
            else:
                # version 1
                #ccd1=Spectrograph.CCD(name='CCD1',  xpix=2032, ypix=4102, pix_size=0.015, xpos=-32.19, ypos=0)
                # updated on 20120325
                ccd1 = CCD(name='CCD1', xpix=2032, ypix=4102, pix_size=0.015, xpos=-32.40, ypos=0.0486)
                ccd2 = CCD(name='CCD1', xpix=2032, ypix=4102, pix_size=0.015, xpos=0, ypos=0)
                ccd3 = CCD(name='CCD1', xpix=2032, ypix=4102, pix_size=0.015, xpos=32.218, ypos=0.0196)
            self.detector = Detector(name=name, ccd=[ccd1, ccd2, ccd3], xbin=xbin, ybin=ybin,
                                     xpos=xpos, ypos=ypos, plate_scale=0.224)
        else:
            raise RSSError('%s is not a supported detector' % name)

    def set_grating(self, name=None):
        if name == 'PG0300':
            self.grating = Grating(name='PG0300', spacing=300)
        elif name == 'PG0900':
            self.grating = Grating(name='PG0900', spacing=903.20)
        elif name == 'PG1300':
            self.grating = Grating(name='PG1300', spacing=1301.20)
        elif name == 'PG1800':
            self.grating = Grating(name='PG1800', spacing=1801.65)
        elif name == 'PG2300':
            self.grating = Grating(name='PG2300', spacing=2302.15)
        elif name == 'PG3000':
            self.grating = Grating(name='PG3000', spacing=2999.98)
        else:
            raise RSSError('%s is not a supported RSS grating' % name)

    def set_slit(self, slitang=1.0):
        self.slit = Slit(name='LongSlit', phi=slitang)
        self.slit.width = self.slit.calc_width(self.telescope.focallength)
