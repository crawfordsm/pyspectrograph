#!/usr/bin/env python
################################# LICENSE ##################################
# Copyright (c) 2009, South African Astronomical Observatory (SAAO)        #
# All rights reserved.                                                     #
#                                                                          #
# Redistribution and use in source and binary forms, with or without       #
# modification, are permitted provided that the following conditions       #
# are met:                                                                 #
#                                                                          #
#     * Redistributions of source code must retain the above copyright     #
#       notice, this list of conditions and the following disclaimer.      #
#     * Redistributions in binary form must reproduce the above copyright  #
#       notice, this list of conditions and the following disclaimer       #
#       in the documentation and/or other materials provided with the      #
#       distribution.                                                      #
#     * Neither the name of the South African Astronomical Observatory     #
#       (SAAO) nor the names of its contributors may be used to endorse    #
#       or promote products derived from this software without specific    #
#       prior written permission.                                          #
#                                                                          #
# THIS SOFTWARE IS PROVIDED BY THE SAAO ''AS IS'' AND ANY EXPRESS OR       #
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED           #
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE   #
# DISCLAIMED. IN NO EVENT SHALL THE SAAO BE LIABLE FOR ANY                 #
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL       #
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS  #
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)    #
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,      #
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN #
# ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE          #
# POSSIBILITY OF SUCH DAMAGE.                                              #
############################################################################
"""
SPECIDENTIFY  is a program to read in SALT RSS spectroscopic arc lamps and
determine the wavelength solution for that data.  The input data should be
a SALT arc lamp and a line list or another arc lamp image with high quality
wavelength solution. The lamp list can be either wavelengths, wavelengths
and fluxes, or an arc image with a high quality solution.  The line lamp
can also be left unspecified and the user will manually enter the data.

From there, the user has several different possible choices.  They can provide
a first guess of the coefficients for the wavelength solution or transformation,
indicate the model for the spectragraph for the first guess, or provide an
image with a solution already as the first guess.

Author                 Version      Date
-----------------------------------------------
S. M. Crawford (SAAO)    1.0       10 Oct 2009

TODO
----


LIMITATIONS
-----------
1. Currently assumes that the linelist is of the form of an ascii file with
   either appropriate information in either one or two columns

"""
# Ensure python 2.5 compatibility
from __future__ import with_statement

import os
import sys
import time
import numpy as np

from pyraf import iraf
import saltprint
import saltio
import saltkey
import saltsafekey
import saltsafeio
from saltsafelog import logging
from salterror import SaltError, SaltIOError


from PySpectrograph import RSSModel
from PySpectrograph import apext
from PySpectrograph import WavelengthSolution
from PySpectrograph import LineSolution
from PySpectrograph.detectlines import detectlines


from . import spectools as st
from .spectools import SALTSpecError
from .InterIdentify import InterIdentify
from .AutoIdentify import AutoIdentify

debug = True


# -----------------------------------------------------------
# core routine

def specidentify(images, linelist, outfile, guesstype, guessfile, function,
                 order, rstep, interact, clobber, logfile, verbose, status):

    with logging(logfile, debug) as log:

        # set up the variables
        infiles = []
        outfiles = []
        status = 0

        # Check the input images
        infiles = saltsafeio.argunpack('Input', images)
        print infiles

        # create list of output files
        outfiles = saltsafeio.argunpack('Input', outfile)

        # if outfiles is a single image, turn it into a list
        # of the same length as infiles
        if len(outfiles) != len(infiles):
            if len(outfiles) == 1:
                outfiles = outfiles * len(infiles)
            elif len(outfiles) == 0:
                outfiles = [None] * len(infiles)
            else:
                msg = 'Please enter an appropriate number of outfiles'
                raise SALTSpecError(msg)

        # open the line lists
        slines, sfluxes = st.readlinelist(linelist)

        # Identify the lines in each file
        for img, oimg in zip(infiles, outfiles):
            log.message('Proccessing image %s' % img)
            identify(img, oimg, slines, sfluxes, guesstype, guessfile, function,
                     order, rstep, interact, clobber, log, verbose)


#------------------------------------------------------------------
# Find the solution for lines in a file

def identify(img, oimg, slines, sfluxes, guesstype, guessfile, function, order,
             rstep, interact, clobber, log, verbose):
    """For a given image, find the solution for each row in the file.  Use the appropriate first guess and
       guess type along with the appropriate function and order for the fit.

       Write out the new image with the solution in the headers and/or as a table in the multi-extension
       fits file

        returns the status
    """
    status = 0
    ImageSolution = {}
    dcstep = 3
    nstep = 50
    res = 2.0
    dres = 0.1
    centerrow = None
    nrows = 1
    sigma = 3
    niter = 5
    xdiff = 2 * res
    method = 'MatchZero'

    # Open up the image
    hdu = saltsafeio.openfits(img)

    # Read in important keywords

    # determine the central row and read it in
    try:
        data = hdu[1].data
        midline = int(0.5 * len(data))
        xarr = np.arange(len(data[midline]))
        specarr = data
    except Exception as e:
        message = 'Unable to read in data array in %s because %s' % (img, e)
        raise SALTSpecError(message)

    # determine the type of first guess.  Assumes none
    if guesstype == 'user':
        pass
    elif guesstype == 'rss':
        dateobs = saltsafekey.get('DATE-OBS', hdu[0], img)
        utctime = saltsafekey.get('UTC-OBS', hdu[0], img)
        instrume = saltsafekey.get('INSTRUME', hdu[0], img)
        grating = saltsafekey.get('GRATING', hdu[0], img)
        grang = saltsafekey.get('GR-ANGLE', hdu[0], img)
        arang = saltsafekey.get('AR-ANGLE', hdu[0], img)
        filter = saltsafekey.get('FILTER', hdu[0], img)
        slit = float(saltsafekey.get('MASKID', hdu[0], img))
        xbin, ybin = saltsafekey.ccdbin(hdu[0], img)
        # set up the rss model
        rssmodel = RSSModel.RSSModel(grating_name=grating.strip(), gratang=grang,
                                     camang=arang, slit=slit, xbin=xbin, ybin=ybin)
        rss = rssmodel.rss
        res = 1e7 * rss.calc_resolelement(rss.gratang, rss.gratang - rss.camang)

        if not instrume in ['PFIS', 'RSS']:
            msg = '%s is not a currently supported instrument' % instrume
            raise SALTSpecError(msg)
        ws = useRSSModel(xarr, rss, function=function, order=order)
    elif guesstype == 'image':
        pass
    else:
        ws = None

    # run in either interactive or non-interactive mode
    if interact:
        ImageSolution = InterIdentify(xarr, specarr, slines, sfluxes, ws, xdiff=xdiff, function=function,
                                      order=order, verbose=True)
    else:
        ImageSolution = AutoIdentify(xarr, specarr, slines, sfluxes, ws,
                                     rstep=rstep, method=method, icenter=centerrow, nrows=nrows,
                                     res=res, dres=dres, dc=dcstep, nstep=nstep, sigma=sigma, niter=niter,
                                     verbose=verbose)

    # set up the list of solutions to into an array
    key_arr = np.array(ImageSolution.keys())
    arg_arr = key_arr.argsort()
    ws_arr = np.zeros((len(arg_arr), len(ws.coef) + 1), dtype=float)

    # write the solution to an array
    for j, i in enumerate(arg_arr):
        if isinstance(ImageSolution[key_arr[i]], WavelengthSolution.WavelengthSolution):
            ws_arr[j, 0] = key_arr[i]
            ws_arr[j, 1:] = ImageSolution[key_arr[i]].coef

    # write the solution as an file
    if outfile:
        # write header to the file that should include the order and function
        if os.path.isfile(outfile) and not clobber:
            dout = open(outfile, 'a')
        else:
            dout = open(outfile, 'w')

        msg = '#WS: Wavelength solution for image %s\n' % img
        msg += '#The following parameters were used in determining the solution:\n'
        msg += '#name=%s\n' % img
        msg += '#time-obs=%s %s\n' % (dateobs, utctime)
        msg += '#instrument=%s\n' % instrume
        msg += '#grating=%s\n' % grating.strip()
        msg += '#graang=%s\n' % grang
        msg += '#arang=%s\n' % arang
        msg += '#filter=%s\n' % filter.strip()
        msg += '#Function=%s\n' % function
        msg += '#Order=%s\n' % order
        msg += '#Starting Data\n'
        dout.write(msg)

        for i in range(len(ws_arr)):
            if ws_arr[i, 0]:
                msg = '%5.2f ' % ws_arr[i, 0]
                msg += ' '.join(['%e' % k for k in ws_arr[i, 1:]])
                dout.write(msg + '\n')
        dout.write('\n')
        dout.close()

    # write the solution as an extension
    hdu.close()

    return


def useRSSModel(xarr, rss, function='poly', order=3):
    """Returns the wavelength solution using the RSS model for the spectrograph


    """

    # now for each position on the detector, calculate the wavelength at that position
    if function == 'poly' or function == 'legendre':
        d = rss.detector.xbin * rss.detector.pix_size * (xarr - 0.5 * len(xarr))
        alpha = rss.gratang
        beta = rss.gratang - rss.camang
        dbeta = -np.degrees(np.arctan(d / rss.camera.focallength))
        y = 1e7 * rss.calc_wavelength(alpha, beta + dbeta)

        # for these models, calculate the wavelength solution
        ws = WavelengthSolution.WavelengthSolution(xarr, y, order=order, function=function)
        ws.fit()
    elif function == 'line':
        ws = LineSolution.LineSolution(rss)
    else:
        message = '%s is not an acceptable form for the function' % function
        raise SALTSpecError(message)

    return ws


# main code

parfile = iraf.osfn("saltspec$specidentify.par")
t = iraf.IrafTaskFactory(taskname="specidentify", value=parfile, function=specidentify, pkgname='saltspec')
