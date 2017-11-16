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

"""INTERIDENTIFY provides an interactive method for identifying
lines in an arc image.  The tasks displays the full image, a
line extracted from the image, and residuals to the fit of that line.
The task will display the total image so the user can extract the lines
to be fit.  Or the user can automate the process so only certain lines are
fit by the user.  On the next tab, the task displays the arc line
and the fit to the line including what lines have been detected and
are being used for the fit.  Finally the task displays the residual in
the fit and the user can select different options to be displayed.

Author                 Version      Date
-----------------------------------------------
S. M. Crawford (SAAO)    1.0       10 Oct 2009

TODO
----


LIMITATIONS
-----------


"""

# Ensure Python 2.5 compatibility


# General imports
import sys
import os
import numpy as np
import pyfits
from pyraf import iraf
from pyraf.iraf import pysalt

# Gui library imports
from PyQt4 import QtGui, QtCore
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QTAgg

# Salt imports
import saltsafeio
from saltgui import ImageDisplay, MplCanvas
from salterror import SaltIOError

from PySpectrograph import apext
from PySpectrograph import Spectrum
#from PySpectrograph import RSSModel
#from PySpectrograph import WavelengthSolution
#from PySpectrograph.detectlines import detectlines

from . import spectools as st
from .spectools import SALTSpecError


class InterIdentifyWindow(QtGui.QMainWindow):

    """Main application window."""

    def __init__(self, xarr, specarr, slines, sfluxes, ws, hmin=150, wmin=400,
                 filename=None, res=2.0, dres=0.1, sigma=5, niter=5,
                 ivar=None, nrows=1, nsteps=100, cmap='gray', scale='zscale', contrast=1.0):
        """Default constructor."""

        # set up the variables
        self.y1 = int(0.5 * len(specarr))
        self.y2 = self.y1 + nrows
        self.specarr = specarr
        self.xarr = xarr
        self.ivar = ivar
        self.slines = slines
        self.sfluxes = sfluxes
        self.hmin = hmin
        self.wmin = wmin
        self.ws = ws
        self.res = res
        self.dres = dres
        self.sigma = sigma
        self.niter = niter
        self.nrows = nrows
        self.cmap = cmap
        self.scale = scale
        self.contrast = contrast
        self.filename = filename
        self.ImageSolution = {}

        # Setup widget
        QtGui.QMainWindow.__init__(self)

        # Set main widget
        self.main = QtGui.QWidget(self)

        # Set window title
        self.setWindowTitle("InterIdentify")

        # create the Image page
        self.imagePage = imageWidget(self.specarr, y1=self.y1, y2=self.y2, hmin=self.hmin, wmin=self.wmin, cmap=self.cmap,
                                     name=self.filename, scale=self.scale, contrast=self.contrast)

        # set up the arc page
        self.farr = apext.makeflat(self.specarr, self.y1, self.y2)

        # set up variables
        self.arcdisplay = ArcDisplay(
            xarr,
            self.farr,
            slines,
            sfluxes,
            self.ws,
            res=self.res,
            dres=self.dres,
            niter=self.niter,
            sigma=self.sigma)
        self.arcPage = arcWidget(self.arcdisplay, hmin=hmin, wmin=wmin, y1=self.y1, y2=self.y2)
        # set up the residual page
        self.errPage = errWidget(self.arcdisplay, hmin=hmin, wmin=wmin)

        # create the tabs
        self.tabWidget = QtGui.QTabWidget()
        self.tabWidget.addTab(self.imagePage, 'Image')
        self.tabWidget.addTab(self.arcPage, 'Arc')
        self.tabWidget.addTab(self.errPage, 'Residual')

        # layout the widgets
        mainLayout = QtGui.QVBoxLayout(self.main)
        mainLayout.addWidget(self.tabWidget)
        # self.setLayout(mainLayout)

        # Set focus to main widget
        # self.main.setFocus()

        # Set the main widget as the central widget
        self.setCentralWidget(self.main)

        # Destroy widget on close
        self.setAttribute(QtCore.Qt.WA_DeleteOnClose)

        # Close when config dialog is closed
        #self.connect(self.conf, QtCore.SIGNAL('destroyed()'), self, QtCore.SLOT('close()'))
        self.connect(self.tabWidget, QtCore.SIGNAL('currentChanged(int)'), self.currentChanged)
        self.connect(self.imagePage, QtCore.SIGNAL('regionChange(int,int)'), self.regionChange)

    def keyPressEvent(self, event):
        print("Key Pressed:", event.key)

    def currentChanged(self, event):
        print(event)

    def regionChange(self, y1, y2):
        self.saveWS()
        print("RegionChange:", y1, y2)
        self.y1 = y1
        self.y2 = y2
        self.farr = apext.makeflat(self.specarr, self.y1, self.y2)
        # set up variables
        print("FINAL WS:", self.ws.coef)
        self.ws = self.newWS(0.5 * (self.y1 + self.y2))
        self.arcdisplay = ArcDisplay(
            self.xarr,
            self.farr,
            self.slines,
            self.sfluxes,
            self.ws,
            res=self.res,
            dres=self.dres,
            niter=self.niter,
            sigma=self.sigma,
            xp=[],
            wp=[])
        self.arcPage = arcWidget(self.arcdisplay, hmin=self.hmin, wmin=self.wmin, y1=self.y1, y2=self.y2)
        # set up the residual page
        self.errPage = errWidget(self.arcdisplay, hmin=self.hmin, wmin=self.wmin)
        # reset the pages
        self.tabWidget.removeTab(2)
        self.tabWidget.removeTab(1)
        self.tabWidget.insertTab(1, self.arcPage, 'Arc')
        self.tabWidget.insertTab(2, self.errPage, 'Residual')

    def saveWS(self):
        self.ws = self.arcdisplay.ws
        k = 0.5 * (self.y1 + self.y2)
        self.ImageSolution[k] = self.ws

    def newWS(self, y):
        """Determine the WS closest to the values given by y1 and y2"""
        keys = np.array(list(self.ImageSolution.keys()))
        print("keys:", keys)
        i = abs(keys - y).argmin()
        print(i, keys[i])
        return self.ImageSolution[keys[i]]


class imageWidget(QtGui.QWidget):

    def __init__(self, imarr, y1=None, y2=None, nrows=1, nsteps=100, hmin=150, wmin=400,
                 name=None, cmap='Gray', scale='zscale', contrast=0.1, parent=None):
        super(imageWidget, self).__init__(parent)

        self.y1 = y1
        self.y2 = y2
        self.x1 = 0
        self.x2 = len(imarr[0])
        self.nrows = nrows
        self.nsteps = nsteps

        # Add FITS display widget with mouse interaction and overplotting
        self.imdisplay = ImageDisplay()
        self.imdisplay.setMinimumHeight(hmin)
        self.imdisplay.setMinimumWidth(wmin)

        # Set colormap
        self.imdisplay.setColormap(cmap)

        # Set scale mode for dynamic range
        self.imdisplay.scale = scale
        self.imdisplay.contrast = contrast
        self.imdisplay.aspect = 'auto'
        self.imdisplay.loadImage(imarr)
        self.imdisplay.drawImage()
        self.y1line, = self.imdisplay.axes.plot([self.x1, self.x2], [self.y1, self.y1], ls='-', color='#00FF00')
        self.y2line, = self.imdisplay.axes.plot([self.x1, self.x2], [self.y2, self.y2], ls='-', color='#00FF00')

        # Add navigation toolbars for each widget to enable zooming
        self.toolbar = NavigationToolbar2QTAgg(self.imdisplay, self)

        # set up the information panel
        self.infopanel = QtGui.QWidget()

        # add the name of the file
        self.NameLabel = QtGui.QLabel("Filename:")
        self.NameLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
        self.NameValueLabel = QtGui.QLabel("%s" % name)
        self.NameValueLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Sunken)

        # add the rows that are extracted
        self.y1Label = QtGui.QLabel("Y1:")
        self.y1Label.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
        self.y1ValueEdit = QtGui.QLineEdit("%6i" % self.y1)
        self.y2Label = QtGui.QLabel("Y2:")
        self.y2Label.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
        self.y2ValueEdit = QtGui.QLineEdit("%6i" % self.y2)
        self.updateButton = QtGui.QPushButton("Update")
        self.updateButton.clicked.connect(self.updatesection)

        # add the update for automatically updating it
        self.nrLabel = QtGui.QLabel("nrows:")
        self.nrLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
        self.nrValueEdit = QtGui.QLineEdit("%5i" % self.nrows)
        self.nsLabel = QtGui.QLabel("nsteps:")
        self.nsLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
        self.nsValueEdit = QtGui.QLineEdit("%6i" % self.nsteps)
        self.nextButton = QtGui.QPushButton("Next")
        self.nextButton.clicked.connect(self.nextsection)

        # set up the info panel layout
        infoLayout = QtGui.QGridLayout(self.infopanel)
        infoLayout.addWidget(self.NameLabel, 0, 0, 1, 1)
        infoLayout.addWidget(self.NameValueLabel, 0, 1, 1, 5)
        infoLayout.addWidget(self.y1Label, 1, 0, 1, 1)
        infoLayout.addWidget(self.y1ValueEdit, 1, 1, 1, 1)
        infoLayout.addWidget(self.y2Label, 1, 2, 1, 1)
        infoLayout.addWidget(self.y2ValueEdit, 1, 3, 1, 1)
        infoLayout.addWidget(self.updateButton, 1, 4, 1, 1)
        infoLayout.addWidget(self.nrLabel, 2, 0, 1, 1)
        infoLayout.addWidget(self.nrValueEdit, 2, 1, 1, 1)
        infoLayout.addWidget(self.nsLabel, 2, 2, 1, 1)
        infoLayout.addWidget(self.nsValueEdit, 2, 3, 1, 1)
        infoLayout.addWidget(self.nextButton, 2, 4, 1, 1)

        # Set up the layout
        mainLayout = QtGui.QVBoxLayout()
        mainLayout.addWidget(self.imdisplay)
        mainLayout.addWidget(self.toolbar)
        mainLayout.addWidget(self.infopanel)
        self.setLayout(mainLayout)

    def updatesection(self):
        self.y1 = int(self.y1ValueEdit.text())
        self.y2 = int(self.y2ValueEdit.text())
        self.y1line.set_ydata([self.y1, self.y1])
        self.y2line.set_ydata([self.y2, self.y2])
        self.imdisplay.draw()
        self.emit(QtCore.SIGNAL("regionChange(int,int)"), self.y1, self.y2)

    def nextsection(self):
        self.nrows = int(self.nrValueEdit.text())
        self.nsteps = int(self.nsValueEdit.text())
        self.y1 = self.y1 + self.nsteps
        self.y2 = self.y1 + self.nrows
        self.y1ValueEdit.setText('%6i' % self.y1)
        self.y2ValueEdit.setText('%6i' % self.y2)
        self.updatesection()


class arcWidget(QtGui.QWidget):

    def __init__(self, arcdisplay, hmin=150, wmin=450, name=None, x1=0, w1=0, y1=None, y2=None, parent=None):
        super(arcWidget, self).__init__(parent)

        # Add FITS display widget with mouse interaction and overplotting
        self.arcdisplay = arcdisplay
        self.arcdisplay.arcfigure.setMinimumHeight(hmin)
        self.arcdisplay.arcfigure.setMinimumWidth(wmin)
        self.arcdisplay.plotArc()

        # Add navigation toolbars for each widget to enable zooming
        self.toolbar = NavigationToolbar2QTAgg(self.arcdisplay.arcfigure, self)

        # set up the information panel
        self.infopanel = QtGui.QWidget()

        # add the name of the file
        self.NameLabel = QtGui.QLabel("Filename:")
        self.NameLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
        self.NameValueLabel = QtGui.QLabel("%s" % name)
        self.NameValueLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Sunken)

        # add the rows that are extracted
        self.y1Label = QtGui.QLabel("Y1:")
        self.y1Label.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
        self.y1ValueLabel = QtGui.QLabel("%6i" % y1)
        self.y1ValueLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Sunken)
        self.y2Label = QtGui.QLabel("Y2:")
        self.y2Label.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
        self.y2ValueLabel = QtGui.QLabel("%6i" % y2)
        self.y2ValueLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Sunken)

        # add in what the value is for a x and w position
        self.x1Label = QtGui.QLabel("X1:")
        self.x1Label.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
        self.w1Label = QtGui.QLabel("w1:")
        self.w1Label.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
        self.x1ValueLabel = QtGui.QLabel("%6.2f" % x1)
        self.x1ValueLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Sunken)
        w1 = self.arcdisplay.ws.value(x1)
        self.w1ValueEdit = QtGui.QLineEdit("%6i" % w1)
        self.addButton = QtGui.QPushButton("Add")
        self.addButton.clicked.connect(self.addpoints)

        # add in radio buttons for pixel or wavelength
        self.pixelradio = QtGui.QRadioButton("Pixel")
        self.wavelengthradio = QtGui.QRadioButton("Wavelength")
        self.pixelradio.setChecked(True)

        # add in information about the order and type of solution
        self.funcLabel = QtGui.QLabel("Function:")
        self.funcLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
        self.funcValueLabel = QtGui.QLabel("%s" % self.arcdisplay.ws.function)
        self.funcValueLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Sunken)
        self.orderLabel = QtGui.QLabel("Order:")
        self.orderLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
        self.orderValueLabel = QtGui.QLabel("%2i" % self.arcdisplay.ws.order)
        self.orderValueLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Sunken)

        # provide the full layout of the information panel
        infoLayout = QtGui.QGridLayout(self.infopanel)
        infoLayout.addWidget(self.NameLabel, 0, 0, 1, 1)
        infoLayout.addWidget(self.NameValueLabel, 0, 1, 1, 5)
        infoLayout.addWidget(self.y1Label, 1, 0, 1, 1)
        infoLayout.addWidget(self.y1ValueLabel, 1, 1, 1, 1)
        infoLayout.addWidget(self.y2Label, 1, 2, 1, 1)
        infoLayout.addWidget(self.y2ValueLabel, 1, 3, 1, 1)
        infoLayout.addWidget(self.x1Label, 2, 0, 1, 1)
        infoLayout.addWidget(self.x1ValueLabel, 2, 1, 1, 1)
        infoLayout.addWidget(self.w1Label, 2, 2, 1, 1)
        infoLayout.addWidget(self.w1ValueEdit, 2, 3)
        infoLayout.addWidget(self.addButton, 2, 4, 1, 1)
        infoLayout.addWidget(self.pixelradio, 3, 0, 1, 2)
        infoLayout.addWidget(self.wavelengthradio, 3, 2, 1, 2)
        infoLayout.addWidget(self.funcLabel, 4, 0, 1, 1)
        infoLayout.addWidget(self.funcValueLabel, 4, 1, 1, 1)
        infoLayout.addWidget(self.orderLabel, 4, 2, 1, 1)
        infoLayout.addWidget(self.orderValueLabel, 4, 3, 1, 1)

        # Set up the layout
        mainLayout = QtGui.QVBoxLayout()
        mainLayout.addWidget(self.arcdisplay.arcfigure)
        mainLayout.addWidget(self.toolbar)
        mainLayout.addWidget(self.infopanel)
        self.setLayout(mainLayout)

        self.connect(self.arcdisplay, QtCore.SIGNAL('keyPressEvent(string)'), self.keyPressEvent)
        self.connect(self.arcdisplay, QtCore.SIGNAL('updatex(float)'), self.updatexlabel)

    def keyPressEvent(self, event):
        print("Arc Widget, keyPress:", event)

    def updatexlabel(self, value):
        try:
            self.x1ValueLabel.setText("%6.2f" % value)
            self.w1ValueEdit.setText("%6.2f" % self.arcdisplay.ws.value(value))
        except TypeError:
            pass

    def addpoints(self):
        """Add the x and w points to the list of matched points"""
        x = float(self.x1ValueLabel.text())
        w = float(self.w1ValueEdit.text())
        #x=[1904.5, 1687.22, 3124.3499999999999, 632.57000000000005]
        #w=[4671.2259999999997, 4624.2757000000001, 4916.5100000000002, 4383.9092000000001]
        self.arcdisplay.addpoints(x, w)


class errWidget(QtGui.QWidget):

    def __init__(self, arcdisplay, hmin=150, wmin=450, name=None, parent=None):
        super(errWidget, self).__init__(parent)

        # Add FITS display widget with mouse interaction and overplotting
        self.arcdisplay = arcdisplay
        self.arcdisplay.errfigure.setMinimumHeight(hmin)
        self.arcdisplay.errfigure.setMinimumWidth(wmin)
        self.arcdisplay.plotErr()

        # Add navigation toolbars for each widget to enable zooming
        self.toolbar = NavigationToolbar2QTAgg(self.arcdisplay.errfigure, self)

        # set up the information panel
        self.infopanel = QtGui.QWidget()

        # add the name of the file
        self.NameLabel = QtGui.QLabel("Filename:")
        self.NameLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
        self.NameValueLabel = QtGui.QLabel("%s" % name)
        self.NameValueLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Sunken)

        # add the rows that are extracted
        self.aveLabel = QtGui.QLabel("Average:")
        self.aveLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
        self.aveValueLabel = QtGui.QLabel("")
        self.aveValueLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Sunken)
        self.stdLabel = QtGui.QLabel("Std:")
        self.stdLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Raised)
        self.stdValueLabel = QtGui.QLabel("")
        self.stdValueLabel.setFrameStyle(QtGui.QFrame.Panel | QtGui.QFrame.Sunken)

        # provide the full layout of the information panel
        infoLayout = QtGui.QGridLayout(self.infopanel)
        infoLayout.addWidget(self.NameLabel, 0, 0, 1, 1)
        infoLayout.addWidget(self.NameValueLabel, 0, 1, 1, 5)
        infoLayout.addWidget(self.aveLabel, 1, 0)
        infoLayout.addWidget(self.aveValueLabel, 1, 1)
        infoLayout.addWidget(self.stdLabel, 1, 2)
        infoLayout.addWidget(self.stdValueLabel, 1, 3)

        # Set up the layout
        mainLayout = QtGui.QVBoxLayout()
        mainLayout.addWidget(self.arcdisplay.errfigure)
        mainLayout.addWidget(self.toolbar)
        mainLayout.addWidget(self.infopanel)
        self.setLayout(mainLayout)

        self.connect(self.arcdisplay, QtCore.SIGNAL('fitUpdate()'), self.fitUpdate)

    def fitUpdate(self):
        try:
            xp = np.array(self.arcdisplay.xp)
            wp = np.array(self.arcdisplay.wp)
            w = self.arcdisplay.ws.value(xp)
            value = (wp - w).mean()
            print(value)
            self.aveValueLabel.setText("%4.2g" % value)
            value = (wp - w).std()
            print(value)
            self.stdValueLabel.setText("%4.2g" % value)
        except Exception as e:
            print(e)


class ArcDisplay(QtGui.QWidget):

    """Class for displaying Arc Spectra using matplotlib and embedded in a Qt 4 GUI.
    """

    def __init__(self, xarr, farr, slines, sfluxes, ws, xp=[], wp=[], xdiff=20,
                 res=2.0, dres=0.1, sigma=5, niter=5):
        """Default constructor."""
        QtGui.QWidget.__init__(self)

        # Initialize base class
        self.arcfigure = MplCanvas()
        self.errfigure = MplCanvas()

        # Add central axes instance
        self.axes = self.arcfigure.figure.add_subplot(111)
        self.erraxes = self.errfigure.figure.add_subplot(111)

        # Connect mouse events
        self.arcfigure.connectMatplotlibMouseMotion()
        self.arcfigure.mpl_connect('button_press_event', self.onButtonPress)
        self.arcfigure.mpl_connect('key_press_event', self.onKeyPress)

        self.errfigure.connectMatplotlibMouseMotion()
        self.errfigure.mpl_connect('button_press_event', self.onButtonPress)
        self.errfigure.mpl_connect('key_press_event', self.onKeyPress)

        # load the data
        self.xarr = xarr
        self.farr = farr
        self.slines = slines
        self.sfluxes = sfluxes
        self.ws = ws
        self.xdiff = xdiff
        self.sigma = sigma
        self.niter = niter
        self.res = res
        self.dres = dres

        self.xp = xp
        self.wp = wp

        # set up the artificial spectra
        self.spectrum = Spectrum.Spectrum(self.slines, self.sfluxes, dw=self.dres, stype='line', sigma=self.res)
        self.swarr = self.spectrum.wavelength
        self.sfarr = self.spectrum.flux * self.farr.max() / self.spectrum.flux.max()

        # set up the wavelength solution
        if self.ws.function == 'line':
            self.ws.set_xarr(self.xarr)
            self.ws.farr = self.farr
            self.ws.spectrum = self.spectrum

        # set up the list of deleted points
        self.dxp = []
        self.dwp = []

        # set up other variables
        self.isArt = False
        self.isFeature = False

        # Set display parameters
        self.xmin = self.xarr.min()
        self.xmax = self.xarr.max()
        self.ymin = self.farr.min()
        self.ymax = self.farr.max()

    def onKeyPress(self, event):
        """Emit signal on key press"""
        if event.key == '?':
            # return the help file
            print('?:', event.key)
        elif event.key == 'c':
            # return the centroid
            if event.xdata:
                cx = st.mcentroid(self.xarr, self.farr, xc=event.xdata, xdiff=self.xdiff)
                self.emit(QtCore.SIGNAL("updatex(float)"), cx)
        elif event.key == 'x':
            # return the x position
            if event.xdata:
                self.emit(QtCore.SIGNAL("updatex(float)"), event.xdata)
        elif event.key == 'f':
            # find the fit
            self.findfit()
            self.emit(QtCore.SIGNAL("fitUpdate()"))
        elif event.key == 'b':
            # auto-idenitfy features
            self.findfeatures()
        elif event.key == 'z':
            # Assume the solution is correct and find the zeropoint
            # that best matches it from cross correllation
            self.findzp()
        elif event.key == 'e':
            # find closest feature from existing fit and line list
            # and match it
            pass
        elif event.key == 'l':
            # plot the features from existing list
            self.plotFeatures()
            self.redraw_canvas()
        elif event.key == 'i':
            # reset identified features
            pass
        elif event.key == 'r':
            # redraw graph
            self.redraw_canvas()
        elif event.key == 'a':
            # draw artificial spectrum
            self.isArt = not self.isArt
            self.redraw_canvas()
        elif event.key == 'd':
            # Delete feature
            save = False
            if event.canvas == self.errfigure:
                save = True
            self.deletepoints(event.xdata, save=save)
            self.redraw_canvas(keepzoom=True)
        elif event.key:
            self.emit(QtCore.SIGNAL("keyPressEvent(string)"), event.key)

    def onButtonPress(self, event):
        """Emit signal on selecting valid image position."""

        if event.xdata and event.ydata:
            print(event.xdata, event.ydata, event.canvas, event.name)
            self.emit(QtCore.SIGNAL("positionSelected(float, float)"),
                      float(event.xdata), float(event.ydata))

    def plotArc(self):
        """Draw image to canvas."""

        # plot the spectra
        self.spcurve, = self.axes.plot(self.xarr, self.farr, linewidth=0.5, linestyle='-', marker=None, color='b')

    def plotArt(self):
        """Plot the artificial spectrum"""
        self.isArt = True
        warr = self.ws.value(self.xarr)
        asfarr = st.interpolate(warr, self.swarr, self.sfarr, left=0.0, right=0.0)
        self.fpcurve, = self.axes.plot(self.xarr, asfarr, linewidth=0.5, linestyle='-',
                                       marker=None, color='r')

    def plotFeatures(self):
        """Plot features identified in the line list"""
        self.isFeature = True
        fl = np.array(self.xp) * 0.0 + 0.25 * self.farr.max()
        self.splines = self.axes.plot(self.xp, fl, ls='', marker='|', ms=20, color='#00FF00')
        # set up the text position
        tsize = 0.83
        self.ymin, self.ymax = self.axes.get_ylim()
        ppp = (self.ymax - self.ymin) / (self.arcfigure.figure.get_figheight() * self.arcfigure.figure.get_dpi())
        f = self.ymax - 10 * tsize * ppp
        for x, w in zip(self.xp, self.wp):
            w = '%6.2f' % float(w)
            self.axes.text(x, f, w, size='small', rotation='vertical', color='#00FF00')

    def plotErr(self):
        """Draw image to canvas."""
        if self.xp and self.wp:
            # plot the spectra
            w = self.ws.value(self.xp)
            self.errcurve, = self.erraxes.plot(self.xp, self.wp - w, linewidth=0.5, linestyle='', marker='o', color='b')
        if self.dxp and self.dwp:
            # plot the spectra
            dw = self.ws.value(self.dxp)
            self.delerrcurve, = self.erraxes.plot(
                self.dxp, self.dwp - dw, linewidth=0.5, linestyle='', marker='x', color='b')

    def findfeatures(self):
        """Given a set of features, find other features that might
           correspond to those features
        """
        xp, wp = st.findfeatures(self.xarr, self.farr, self.slines, self.sfluxes,
                                 self.ws, sigma=self.sigma, niter=self.niter)
        for x, w in zip(xp, wp):
            if w not in self.wp and w > -1:
                self.xp.append(x)
                self.wp.append(w)
        self.plotFeatures()
        self.redraw_canvas()

    def findzp(self):
        """Find the zeropoint for the source and plot of the new value
        """
        self.ws = st.findzeropoint(self.xarr, self.farr, self.swarr, self.sfarr,
                                   self.ws, dc=10, nstep=20, inttype='interp')
        self.plotArt()
        self.redraw_canvas()

    def findfit(self):
        print(self.xp)
        print(self.wp)
        self.ws = st.findfit(self.xp, self.wp, function=self.ws.function, order=self.ws.order)
        self.err_redraw_canvas()

    def addpoints(self, x, w):
        """Add points to the line list
        """
        if isinstance(x, list) and isinstance(w, list):
            self.xp.extend(x)
            self.wp.extend(w)
        else:
            self.xp.append(x)
            self.wp.append(w)
        print(self.xp)

    def deletepoints(self, x, save=False):
        """ Delete points from the line list
        """
        in_minw = np.abs(np.array(self.xp) - x).argmin()
        if save:
            self.dxp.append(self.xp[in_minw])
            self.dwp.append(self.wp[in_minw])
        self.xp.__delitem__(in_minw)
        self.wp.__delitem__(in_minw)

    def redraw_canvas(self, keepzoom=False):
        if keepzoom:
            # Store current zoom level
            xmin, xmax = self.axes.get_xlim()
            ymin, ymax = self.axes.get_ylim()

        # Clear plot
        self.axes.clear()

        # Draw image
        self.plotArc()

        # if necessary, redraw the features
        if self.isFeature:
            self.plotFeatures()

        # if necessary, draw the artificial spectrum
        if self.isArt:
            self.plotArt()

        # Restore zoom level
        if keepzoom:
            self.axes.set_xlim((self.xmin, self.xmax))
            self.axes.set_ylim((self.ymin, self.ymax))

        # Force redraw
        self.arcfigure.draw()

        self.err_redraw_canvas()

    def err_redraw_canvas(self, keepzoom=False):
        if keepzoom:
            # Store current zoom level
            xmin, xmax = self.erraxes.get_xlim()
            ymin, ymax = self.erraxes.get_ylim()
        else:
            self.xmin, self.xmax = self.axes.get_xlim()

        # Clear plot
        self.erraxes.clear()

        # Draw image
        self.plotErr()

        # Restore zoom level
        if keepzoom:
            self.erraxes.set_xlim((xmin, xmax))
            self.erraxes.set_ylim((ymin, ymax))
        else:
            self.erraxes.set_xlim((self.xmin, self.xmax))

        self.errfigure.draw()

        self.emit(QtCore.SIGNAL("fitUpdate()"))


def InterIdentify(xarr, specarr, slines, sfluxes, ws, xdiff=20, function='poly',
                  filename=None, order=3, scale='zscale', cmap='gray', contrast=1.0, verbose=True):

    # Create GUI
    App = QtGui.QApplication(sys.argv)
    aw = InterIdentifyWindow(
        xarr,
        specarr,
        slines,
        sfluxes,
        ws,
        cmap=cmap,
        scale=scale,
        contrast=contrast,
        filename=filename)
    aw.show()

    # Start application event loop
    exit = App.exec_()

    # Check if GUI was executed succesfully
    if exit != 0:
        raise SALTSpecError('InterIdentify GUI has unexpected exit status ' + str(exit))

    return aw.ImageSolution
