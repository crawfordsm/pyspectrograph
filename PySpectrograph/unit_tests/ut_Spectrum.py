#import Spectrum
from PySpectrograph import Spectrum
import numpy as np
import pylab as pl

infile='Xe.10.spec'
inlist='Xe.dat'
oneline='one.line'

def test_spectrum():
   print 'testing spectrum...'

   #create a spectrum from a single line
   w,f=np.loadtxt(oneline, usecols=(0,1),unpack=True)
   sp1=Spectrum.Spectrum([w], [f], wrange=[4500, 4900], sigma=1)
   pl.figure(figsize=(10,10))
   pl.axes([0.1, 0.1, 0.8, 0.8])
   pl.plot(sp1.wavelength,sp1.flux)

   #create a spectrum from a line list
   w,f=np.loadtxt(inlist, usecols=(0,1),unpack=True)
   spl=Spectrum.Spectrum(w, f, sigma=5)
   pl.plot(spl.wavelength,spl.flux)

   #create a spectrum from a spectrum
   w,f=np.loadtxt(infile, usecols=(0,1),unpack=True)
   spf=Spectrum.Spectrum(w, f, sigma=5, stype='continuum')
   pl.plot(spf.wavelength,spf.flux)

   #test error that mean is the wrong number of dimensions
   pl.savefig('out.png')
   pl.show()

def test_vacuum():
    #test code to change vacuum to air wavlengths and back again
    #test case form IDL script
    w_air=6056.125
    w_vac=6057.8019
    if abs(Spectrum.air2vac(w_air, mode='Morton')-w_vac)>0.01: 
       print "ERROR in MORTON AIR2VAC calculation"
    if abs(Spectrum.air2vac(w_air, mode='Ciddor')-w_vac)>0.01: 
       print "ERROR in CIDDOR AIR2VAC calculation"
    print w_air, w_vac
    if abs(Spectrum.vac2air(w_vac, mode='Ciddor')-w_air)>0.01:
       print "ERROR in CIDDOR VAC2AIR calculation"
    if abs(Spectrum.vac2air(w_vac, mode='Morton')-w_air)>0.01: 
       print "ERROR in MORTON VAC2AIR calculation"

    #check that it works with a 

test_spectrum()
test_vacuum()
