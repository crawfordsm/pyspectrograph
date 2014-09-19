#
# MAKEPLOTS--A library for making plots for demaniacs
#
#
#

from pylab import *
import numpy


def plotframe(data):
    """Plot the entire data array
       returns a figure
    """
    nimg = 10
    ywidth = 0.08
    xlen = len(data[0]) / nimg
    for i in range(nimg):
        yax = 0.90 - ywidth * 1.1 * i
        x1 = xlen * i
        x2 = xlen * (1 + i)
        f = axes([0.1, yax, 0.8, ywidth])
        f.imshow(data[:, x1:x2], cmap=cm.gray, aspect='auto', vmin=-5, vmax=50)
        f.axis('off')
    return f


def plotfeature(f, wave, data, w1, w2, z):
    """Plot a section of the data array
       as indicated by w1 and w2
    """
    w1 = w1 * (1 + z)
    w2 = w2 * (1 + z)
    if w1 > wave.max():
        return f
    mask = (w1 < wave) * (wave < w2)
    mdata = data[:, mask]
    f.imshow(mdata, cmap=cm.gray, aspect='auto', vmin=-5, vmax=50)
    # set up the axis labels
    x1 = wave[mask][0]
    x2 = wave[mask][-1]
    dw = (x2 - x1) / 5
    xtarr = arange(x1, x2 + dw, dw)
    xtlab = []
    for x in xticks()[0]:
        if x >= 0 and x < len(wave[mask]):
            x = wave[mask][x]
            xtlab.append('%4.2f' % x)
        else:
            xtlab.append('0')
    f.set_yticklabels([])
    f.set_xticklabels([])
    return f


def plotlinefeature(f, wave, flux, w1, w2, z):
    w1 = w1 * (1 + z)
    w2 = w2 * (1 + z)
    mask = (w1 < wave) * (wave < w2)
    f = plotline(f, wave[mask], flux[mask])
    return f


def plotline(f, wave, flux, color=None):
    if color:
        f.plot(wave, flux, ls='-', color=color, lw=1.55)
    else:
        f.plot(wave, flux, ls='-', lw=1.55)
    f.set_xlim((wave[0], wave[-1]))
    # f.set_yticklabels([])
    return f
