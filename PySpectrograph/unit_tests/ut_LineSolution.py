
import numpy as np
import pylab as pl

from PySpectrograph.WavelengthSolution import LineSolution as LS

xp=np.array([1134.87, 1239.11, 1498.22, 1687.21, 1904.49, 1997.15, 2025.77, 2202.59, 2559.72, 2673.74, 3124.34])
wp=np.array([4500.9772, 4524.6805, 4582.7474, 4624.2757, 4671.226, 4690.9711, 4697.02, 4734.1524, 4807.019, 4829.709,  4916.51] )


def test_LineSolution():
   ls=LS.LineSolution(xp, wp, function='legendre')
   ls.interfit()
   print ls.coef
   print ls.sigma(ls.x, ls.y)
   print ls.value(2000)

   pl.figure()
   pl.plot(xp, wp-ls(xp), ls='', marker='o')
   pl.plot(ls.x, ls.y-ls(ls.x), ls='', marker='o')
   pl.show()
   return

test_LineSolution()
