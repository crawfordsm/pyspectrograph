import numpy as np
from detectlines import *
import pylab as pl

infile = 'Xe.01.spec'
inlist = 'Xe.dat'
plist = [4500.98500001, 4624.28500002, 4671.23500002, 4734.16500002, 4807.02500002, 4916.51500002, 4923.16500002]
wlist = [4500.9772, 4624.2757, 4671.226, 4734.1524, 4807.019, 4916.51, 4923.152]


def centroid2(warr, farr, mask):
    return (warr[mask] * farr[mask]).sum() / farr[mask].sum()


def test_centroid():
    xc = 30.2
    gaussian = lambda x: 3 * np.exp(-0.5 * (xc - x) ** 2 / (3.16 ** 2.))
    x = np.arange(-50, 50)
    f = gaussian(x)

    xc1 = centroid(x, f)

    xdiff = 10
    mask = (abs(x - xc) < xdiff)

    xcm = centroid(x, f, mask=mask)

    if abs(xc - xc1) > 0.01:
        print("FAIL", xc, xc1)
    if abs(xc - xcm) > 0.01:
        print("FAIL", xc, xcm)


def test_find_backstats():
    w, f = np.loadtxt(infile, unpack=True)
    ave, std = find_backstats(f, sigma=3, niter=5)
    if abs(ave - 50) > 1 or (std - 10) > 1:
        print("FAIL")


def test_find_peaks():
    warr, farr = np.loadtxt(infile, unpack=True)
    xp = find_peaks(farr, 10, 5)
    wdiff = 10
    for x, w in zip(xp, plist):
        if warr[x] != w:
            print("FAIL", warr[x], w)


def test_detectlines():
    wl, fl = np.loadtxt(inlist, usecols=[0, 1], unpack=True)
    warr, farr = np.loadtxt(infile, unpack=True)
    wp = detectlines(warr, farr, sigma=10, niter=5, bsigma=3, mask=None, kern=default_kernal, center=False)
    cwp = detectlines(warr, farr, sigma=10, niter=5, bsigma=3, mask=None, kern=default_kernal, center=True)
    if (np.array(wlist, dtype=float) - wp).mean() > 0.01:
        print('FAIL center=False')
    if (np.array(wlist, dtype=float) - cwp).mean() > 0.01:
        print('FAIL center=True')

test_centroid()
test_find_backstats()
test_find_peaks()
test_detectlines()
