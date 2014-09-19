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
AutoIDENTIFY  is a program to automatically identify spectral lines in
an arc image.

Author                 Version      Date
-----------------------------------------------
S. M. Crawford (SAAO)    1.0       21 Aug 2010

TODO
----


LIMITATIONS
-----------

"""
# Ensure python 2.5 compatibility

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
from PySpectrograph.detectlines import detectlines


from . import spectools as st
from .spectools import SALTSpecError

debug = True


def AutoIdentify(xarr, specarr, slines, sfluxes, ws, method='Zeropoint',
                 rstep=1, icenter=None, nrows=1, res=2, dres=0.1,
                 sigma=5, niter=5,
                 dc=20, nstep=20,
                 verbose=True):
    """Automatically find the wavlength solution for the entire image.  The following
       methods are used:

       zeropoint--Assume that the form for the initial guess of the wavelength solution
                  is correct
    """
    ImageSolution = {}

    # run it if only the zeropoint needs to be calculated
    if method == 'Zeropoint':
        func = st.findzeropoint
        ImageSolution = runsolution(xarr, specarr, slines, sfluxes, ws, func, fline=False, oneline=False,
                                    rstep=rstep, icenter=icenter, nrows=nrows, res=res, dres=dres,
                                    dsigma=sigma, dniter=niter, verbose=verbose, dc=dc, nstep=nstep)
        print method

    # use a line matching algorithm to match the lines
    # in the image with those in the line list
    if method == 'Matchlines':
        func = st.findwavelengthsolution
        ImageSolution = runsolution(xarr, specarr, slines, sfluxes, ws, func, fline=True, oneline=False,
                                    rstep=rstep, icenter=icenter, nrows=nrows, res=res, dres=dres,
                                    dsigma=sigma, dniter=niter, verbose=verbose, sigma=sigma, niter=niter)

    # first fit a zeropoint, then match the lines, and then
    # find the rest of the points by using only the zeropoint
    if method == 'MatchZero':
        print ws.coef
        func = st.findzeropoint
        ws = runsolution(xarr, specarr, slines, sfluxes, ws, func, fline=False, oneline=True,
                         rstep=rstep, icenter=icenter, nrows=nrows, res=res, dres=dres,
                         dsigma=sigma, dniter=niter, verbose=verbose, dc=10, nstep=20)

        func = st.findwavelengthsolution
        ws = runsolution(xarr, specarr, slines, sfluxes, ws, func, fline=True, oneline=True,
                         rstep=rstep, icenter=icenter, nrows=nrows, res=res, dres=dres,
                         dsigma=sigma, dniter=niter, verbose=verbose, sigma=sigma, niter=niter)
        print 'Running zero now'
        func = st.findzeropoint
        ImageSolution = runsolution(xarr, specarr, slines, sfluxes, ws, func, fline=False, oneline=False,
                                    rstep=rstep, icenter=icenter, nrows=nrows, res=res, dres=dres,
                                    dsigma=sigma, dniter=niter, verbose=verbose, dc=dc, nstep=nstep)

        print method

    if method == 'FullXcor':
        func = st.findxcor
        dcoef = ws.coef * 0.1
        dcoef[-1] = dc
        print dcoef
        ws = runsolution(xarr, specarr, slines, sfluxes, ws, func, fline=True, oneline=True,
                         rstep=rstep, icenter=icenter, nrows=nrows, res=res, dres=dres,
                         dsigma=sigma, dniter=niter, verbose=verbose, dcoef=dcoef)
        print 'Running zero now'
        func = st.findzeropoint
        ImageSolution = runsolution(xarr, specarr, slines, sfluxes, ws, func, fline=False, oneline=False,
                                    rstep=rstep, icenter=icenter, nrows=nrows, res=res, dres=dres,
                                    dsigma=sigma, dniter=niter, verbose=verbose, dc=dc, nstep=nstep)

        print method

    return ImageSolution


def runsolution(xarr, specarr, slines, sfluxes, ws, func, ivar=None,
                fline=True, oneline=False,
                rstep=20,
                icenter=None, nrows=1, dsigma=5, dniter=5, res=2.0, dres=0.1, verbose=True, **kwargs):
    """Starting in the middle of the image, it will determine the solution
       by working its way out to either edge and compiling all the results into
       ImageSolution

       xarr--Full range in x of pixels to solve for

       specarr--Input 2D flux

       func--function to use for the solution


    """
    # set up the variables
    ImageSolution = {}

    # Setup the central line if it isn't specified
    if icenter is None:
        icenter = int(0.5 * len(specarr))

    # set up the flux from the central line (or the line specified by the user in icenter)
    specext = apext.apext(xarr, specarr, ivar=ivar)
    farr = apext.makeflat(specarr, icenter, icenter + nrows)
    farr = st.flatspectrum(xarr, farr, mode='poly', order=2)
    cxp = st.detectlines(xarr, farr, dsigma, dniter)
    nlines = len(cxp)

    # first set up the artificial spectrum
    swarr, sfarr = st.makeartificial(slines, sfluxes, farr.max(), res, dres)

    # find the solution for the central wavelegnth
    k = icenter
    min_lines = 0.1 * len(cxp)
    if fline:
        mws = solution(xarr, specarr, slines, sfluxes, ws, func, k, k + nrows,
                       min_lines=min_lines, dsigma=dsigma, dniter=dniter, **kwargs)
    else:
        mws = solution(xarr, specarr, swarr, sfarr, ws, func, k, k + nrows,
                       min_lines=min_lines, dsigma=dsigma, dniter=dniter, **kwargs)

    print 'runsolution:', mws.coef
    if oneline:
        return mws

    ImageSolution[k] = mws

    # now loop through each step, and calculate the wavelengths for the given
    for i in range(rstep, int(0.5 * len(specarr)), rstep):
        for k in [icenter - i, icenter + i]:
            lws = getwsfromIS(k, ImageSolution)
            if fline:
                fws = solution(xarr, specarr, slines, sfluxes, lws, func, k, k + nrows,
                               min_lines=min_lines, dsigma=dsigma, dniter=dniter, **kwargs)
            else:
                fws = solution(xarr, specarr, swarr, sfarr, lws, func, k, k + nrows,
                               min_lines=min_lines, dsigma=dsigma, dniter=dniter, **kwargs)
            ImageSolution[k] = fws

            if verbose:
                p_new = i * 100.0 / (0.5 * len(specarr))
                ctext = 'Percentage Complete: %d %d %f\r' % (i, p_new, time.clock())  # p_new
                sys.stdout.write(ctext)
                sys.stdout.flush()

    return ImageSolution


def solution(xarr, specarr, sl, sf, ws, func, y1, y2, min_lines=2, dsigma=5, dniter=3, pad=50, **kwargs):
    """Extract a single line and calculate the wavelneght solution"""

    # set up the flux from the set of lines
    farr = apext.makeflat(specarr, y1, y2)
    farr = st.flatspectrum(xarr, farr, mode='poly', order=2)

    # check to see if there are any points
    xp = st.detectlines(xarr, farr, dsigma, dniter)

    if len(xp) > min_lines and ws:
        # make the artificial list
        wmin = ws.value(xarr.min())
        wmax = ws.value(xarr.max())
        smask = (sl > wmin - pad) * (sl < wmax + pad)

        # fit the function
        fws = func(xarr, farr, sl[smask], sf[smask], ws, **kwargs)
        return fws

    return None


def getwsfromIS(k, ImageSolution):
    """From the imageSolution dictionary, find the ws which is nearest to the value k

    """
    ISkeys = np.array(ImageSolution.keys())
    ws = ImageSolution[ISkeys[abs(ISkeys - k).argmin()]]
    if ws is None:
        dist = abs(ISkeys[0] - k)
        ws = ImageSolution[ISkeys[0]]
        for i in ISkeys:
            if ImageSolution[i] and abs(i - k) < dist:
                dist = abs(i - k)
                ws = ImageSolution[i]
    return ws
