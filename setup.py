#!/usr/bin/env python

from distutils.core import setup
DISTUTILS_DEBUG=True


setup(name='PySpectrograph',
      version='0.3',
      description='Spectrograph Modelling Software',
      author='Steve Crawford',
      author_email='crawfordsm@gmail.com',
      url='http://code.google.com/p/pyspectrograph/',
      packages=['PySpectrograph', 'PySpectrograph/Spectrograph', 'PySpectrograph/Utilities', 'PySpectrograph/WavelengthSolution', 'PySpectrograph/Spectra', 'PySpectrograph/Models'],
     )

