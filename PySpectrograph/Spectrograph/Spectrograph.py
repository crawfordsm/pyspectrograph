"""Spectrograph is a class that general describes a spectrograph.  This includes
describing the telescope, slit,  collimator, grating, camera, and detector.

HISTORY
20090912 SMC  First written by SM Crawford

Limitations:
-Still need to verify how alpha, grating angle, beta, and camera angle
 to see if I can hardwire some of the tasks

"""


import math

from .SpectrographEquations import *
from .Grating import Grating
from .Optics import Optics
from .Slit import Slit
from .Detector import Detector
from .CCD import CCD


class Spectrograph(Grating, Optics, Slit, Detector):

    """A class describing a spectrograph and functions
       related to a spectrograph.  All angles are in degrees.
    """

    def __init__(self, camang=45, gratang=45, grating=Grating(), camera=Optics(),
                 collimator=Optics(), telescope=Optics(), slit=Slit(),
                 detector=Detector()):

        # initiate the grating
        self.grating = grating

        # initiate the telescope
        self.telescope = telescope

        # initiate the collimator
        self.collimator = collimator

        # initiate the camera
        self.camera = camera

        # initiate the slit
        self.slit = slit

        # initiate the detector
        self.detector = detector

        # set up the angles in the system
        self.gratang = gratang
        self.camang = camang

        return

    def alpha(self):
        return self.gratang

    def beta(self):
        return self.camang - self.gratang

    def gamma(self):
        return self.gamma

    def calc_wavelength(self, alpha, beta, gamma=0.0, nd=n_index):
        """Apply the grating equation to determine the wavelength
           returns wavelength in mm
        """
        w = gratingequation(self.grating.sigma, self.grating.order, self.grating.sign, alpha, beta, gamma=gamma, nd=nd)
        return w

    def calc_angdisp(self, beta):
        """Calculate the angular dispersion according to m/sigma/cos beta

           returns angular dispersion in 1/mm
        """
        A = calc_angdisp(self.grating.sigma, self.grating.order, beta)
        return A

    def calc_lindisp(self, beta):
        """Calculate the linear dispersion according to f_cam * A

           return linear dispersion in mm/mm

        """
        return calc_lindisp(self.camera.focallength, self.grating.sigma, self.grating.order, beta)

    def calc_demagspatial(self):
        """Calculate the spatial demagnification

           returns the spatial demagnification
        """
        return calc_demagspatial(self.collimator.focallength, self.camera.focallength)

    def calc_demagspectral(self, alpha, beta):
        """Calculate the spectral demagnification

           returns the spectral demagnification
        """
        return self.calc_demagspatial() / se.calc_anamorph(alpha, beta)

    def calc_spatslitimage(self):
        """Calculate the spatial extant of the slit image

           return in mm
        """
        return self.slit.width / self.calc_demagspatial()

    def calc_specslitimage(self, beta):
        """Calculate the spectral extant of the slit image

           return in mm
        """
        return self.slit.width * self.calc_lindisp(beta)

    def calc_resolelement(self, alpha, beta):
        """Calculate the resolution of a single element for a filled slit

           return the wavelength resolution in mm
        """
        dw = calc_resolelement(self.slit.width, self.collimator.focallength,
                               self.grating.sigma, self.grating.order,
                               alpha, beta)
        return dw

    def calc_resolution(self, w, alpha, beta):
        """Calculate the resolution at a given wavelength.  w/dw

           returns resolution
        """
        return w / self.calc_resolelement(alpha, beta)

    def calc_centralwavelength(self):
        """Calculate the central wavlength

           return waveleng in mm
        """
        return self.calc_wavelength(self.alpha(), -self.beta())

    def calc_redwavelength(self):
        """For the detector, calculate the maximum red wavelength
           Assume just the width of the detector

           return waveleng in mm
        """
        dbeta = math.degrees(math.atan(0.5 * self.detector.width / self.camera.focallength))
        return self.calc_wavelength(self.alpha(), -self.beta() - dbeta)

    def calc_bluewavelength(self):
        """For the detector, calculate the maximum blue wavelength
           Assume just the width of the detector

           return waveleng in mm
        """
        dbeta = math.degrees(math.atan(0.5 * self.detector.width / self.camera.focallength))
        return self.calc_wavelength(self.alpha(), -self.beta() + dbeta)
