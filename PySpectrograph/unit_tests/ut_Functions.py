
from Functions import Normal, FunctionsError
import numpy as np
import pylab as pl


def test_normal():

    xarr = np.arange(-25, 25, 0.1)
    arr = np.indices([50, 50]) - 25
    d3arr = np.indices([50, 50, 50]) - 25
    # create a 1-D gaussian
    g1d = Normal(xarr, 10, 5, 15)
    pl.figure(figsize=(10, 10))
    pl.axes([0.1, 0.6, 0.8, 0.3])
    pl.plot(xarr, g1d)

    # create a 2-D gaussian
    g2d = Normal(arr, np.array([10, 0]), np.array([5, 10]), 15)
    pl.axes([0.1, 0.1, 0.3, 0.3])
    zmap = np.product(g2d, axis=0)
    pl.imshow(zmap)

    # create a 3-D gaussian
    g3d = Normal(d3arr, 10, 5, 15)
    z3map = np.product(g3d, axis=0)
    m3d = z3map[35].max()
    pl.axes([0.5, 0.1, 0.2, 0.2])
    pl.imshow(z3map[10], vmin=0, vmax=m3d)
    pl.axes([0.5, 0.35, 0.2, 0.2])
    pl.imshow(z3map[33], vmin=0, vmax=m3d)
    pl.axes([0.75, 0.1, 0.2, 0.2])
    pl.imshow(z3map[35], vmin=0, vmax=m3d)
    pl.axes([0.75, 0.35, 0.2, 0.2])
    pl.imshow(z3map[37], vmin=0, vmax=m3d)
    # test error that mean is the wrong number of dimensions
    try:
        g2d = Normal(arr, np.array([10, 10, 0]), np.array([5, 10]), 15)
    except FunctionsError:
        pass

    pl.show()


test_normal()
