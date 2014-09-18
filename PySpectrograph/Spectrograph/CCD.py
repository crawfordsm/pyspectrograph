import numpy as np


class CCD:

    """Defines a CCD by x and y position, size, and pixel size.  The x and y position are
       set such that they are zero relative to the detector position.  This assumes that
       the x and y positions are in the center of the pixels and that the ccd is symmetric.

       pix_size is in mm
    """

    def __init__(self, name='', height=0, width=0, xpos=0, ypos=0,
                 pix_size=0.015, xpix=2048, ypix=2048):
        # set the variables
        self.xpos = xpos
        self.ypos = ypos
        self.pix_size = pix_size
        self.xpix = xpix
        self.ypix = ypix
        self.height = self.set_height(height)
        self.width = self.set_width(width)

    def set_width(self, w):
        """If the width  is less than the number of pixels, then the width is
           given by the number of pixels
        """
        min = self.xpix * self.pix_size
        return max(w, min)

    def set_height(self, h):
        """If the height is less than the number of pixels, then the height is
           given by the number of pixels
        """
        min = self.ypix * self.pix_size
        return max(h, min)

    def find_corners(self):
        """Return the corners of the ccd"""
        x1 = self.xpos - 0.5 * self.width
        x2 = self.xpos + 0.5 * self.width
        y1 = self.ypos - 0.5 * self.height
        y2 = self.ypos + 0.5 * self.height
        return x1, x2, y1, y2
