#
# APEXT--Extract a 1-D spectra from a 2-D image
#
#
#
import math
import numpy as np


def makeflat(arr, y1, y2):
    """For a 2-D array, compress it along the y-dimension to make a one dimensional
       array
    """
    y1 = max(0, y1)
    y2 = min(len(arr), y2)
    # return the 1-D array if only a single line is requested
    if abs(y1 - y2) <= 1:
        return arr[y1, :]

    # sum up the array along for the full value
    return arr[y1:y2, :].sum(axis=0)


class apext:

    """A class for extracting a 1-D spectra from a 2-D image"""

    def __init__(self, wave, data, ivar=None):
        self.data = data
        self.wave = wave
        self.nwave = len(wave)
        if ivar is None:
            self.ivar = None
        else:
            self.ivar = ivar

    def flatten(self, y1, y2):
        """Compress the 2-D array down to 1-D. This is just either done by
           a straight summation or by a weighted summation if IVAR is present
        """
        if self.ivar is None:
            self.ldata = makeflat(self.data, y1, y2)
            # self.lvar=self.data[y1:y2,:].std(axis=0)
            self.lvar = self.ldata
        else:
            self.lvar = ((self.ivar[y1:y2, :] ** 2).sum(axis=0))
            # self.ldata=(self.data[y1:y2,:]*self.ivar[y1:y2,:]).sum(axis=0)/self.lvar
            self.ldata = makeflat(self.data, y1, y2)
        return
