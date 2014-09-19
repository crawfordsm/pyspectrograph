"""Spectrograph package is part of PySpectrograph and provides
a model for a grating spectrograph.

This task includes:

Spectrograph--The basic equations and components of a spectrograph
Models--models for different spectrographs
Spectra--classes for the handling of spectra
WavelengthSolutions--Methods for calculating the wavelength solution
Identify--Tasks to measure the wavelength solution
Utilities--General useful tasks

"""

from .Spectrograph import Spectrograph
from .Grating import Grating
from .Optics import Optics
from .Slit import Slit
from .Detector import Detector
from .CCD import CCD
