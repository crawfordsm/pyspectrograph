"""The PySpectrograph package holds the tasks to model and analyze spectroscopic
data from grating spectrographs.

These tasks include:

Spectrograph--The basic equations and components of a spectrograph
Models--models for different spectrographs
Spectra--classes for the handling of spectra
WavelengthSolutions--Methods for calculating the wavelength solution
Identify--Tasks to measure the wavelength solution 
Utilities--General useful tasks

"""
class SpectrographError(Exception):
   """Exception Raised for Spectrograph errors"""
   pass

class PySpectrographError(Exception):
   """Exception Raised for PySpectrograph errors"""
   pass

import Utilities
#import Spectrograph
from Spectrograph import *
import Models
import Spectra
from Spectra import *
import WavelengthSolution
#import Identify

__all__=['Identify','Models','Spectra','Spectrograph','Utilities','WavelengthSolution']
__version__=0.23

#general class for errors



