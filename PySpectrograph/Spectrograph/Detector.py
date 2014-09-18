import numpy as np
from .CCD import CCD


class Detector(CCD):

    """A class that describing the Detector.  It inherets from the CCD class as there could be
       multiple ccds at each position.

       name--Name of the detector
       ccd--a CCD class or list describing the CCDs in the detecfor
       xpos--Offset of the x center of the ccd from the central ray in mm
       ypos--Offset of the y center of the ccd from the central ray in mm
       zpos--Offset of the z center of the ccd from the central ray in mm
       xbin--ccd binning in x-direction
       ybin--ccd binning in y-direction
       plate_scale--plate scale in mm/"
    """

    def __init__(self, name='', ccd=CCD(), zpos=0, xpos=0, ypos=0, xbin=2, ybin=2, plate_scale=0.224):

        # Set the detector up as a list of CCDs.
        self.detector = []
        self.pix_size = None
        if isinstance(ccd, CCD):
            self.detector = [ccd]
            self.pix_size = ccd.pix_size
        elif isinstance(ccd, list):
            for c in ccd:
                if isinstance(c, CCD):
                    self.detector.append(c)
                    if self.pix_size:
                        self.pix_size = min(self.pix_size, c.pix_size)
                    else:
                        self.pix_size = c.pix_size
        else:
            return

        self.nccd = len(self.detector)

        # set up the zero points for the detector
        self.name = name
        self.zpos = zpos
        self.xpos = xpos
        self.ypos = ypos
        self.xbin = xbin
        self.ybin = ybin
        self.plate_scale = plate_scale
        self.pix_scale = self.plate_scale / self.pix_size

        # check to make sure that the ccds don't overlap
        self.real = self.check_ccds()

        # determine the max width and height for the detector
        self.width = self.find_width()

        # determine the max width and height for the detector
        self.height = self.find_height()

    def check_ccds(self):
        """Check to make sure none of the ccds overlap"""
        if self.nccd <= 1:
            return True

        # loop over each ccd and check to see if any of the ccd
        # overlaps with the coordinates of another ccd
        for i in range(self.nccd):
            ax1, ax2, ay1, ay2 = self.detector[i].find_corners()
            for j in range(i + 1, self.nccd):
                bx1, bx2, by1, by2 = self.detector[j].find_corners()
                if ax1 <= bx1 < ax2 or ax1 < bx2 < ax2:
                    if ay1 <= by1 < ay2 or ay1 < by2 < ay2:
                        return False

        return True

    def get_xpixcenter(self):
        """Return the xpixel center based on the x and y position"""
        return int((0.5 * self.find_width() - self.xpos) / self.pix_size / self.xbin)

    def get_ypixcenter(self):
        """Return the xpixel center based on the x and y position"""
        return int((0.5 * self.find_height() - self.ypos) / self.pix_size / self.ybin)

    def find_width(self):
        """Loop over all the ccds in detector and find the width"""
        width = 0
        # return zero if no detector
        if self.nccd < 1:
            return width

        # handle a single detector
        width = self.detector[0].width
        if self.nccd == 1:
            return width

        # Loop over multipe CCDs to find the width
        ax1, ax2, ay1, ay2 = self.detector[0].find_corners()
        xmin = min(ax1, ax2)
        xmax = max(ax1, ax2)
        for ccd in self.detector[1:]:
            ax1, ax2, ay1, ay2 = ccd.find_corners()
            xmin = min(xmin, ax1, ax2)
            xmax = max(xmax, ax1, ax2)
        return xmax - xmin

    def find_height(self):
        """Loop over all the ccds in detector and find the height"""
        height = 0
        # return zero if no detector
        if self.nccd < 1:
            return height

        # handle a single detector
        height = self.detector[0].height
        if self.nccd == 1:
            return height

        # Loop over multipe CCDs to find the width
        ax1, ax2, ay1, ay2 = self.detector[0].find_corners()
        ymin = min(ay1, ay2)
        ymax = max(ay1, ay2)
        for ccd in self.detector[1:]:
            ax1, ax2, ay1, ay2 = ccd.find_corners()
            ymin = min(ymin, ay1, ay2)
            ymax = max(ymax, ay1, ay2)
        height = ymax - ymin
        return height

    def make_detector(self):
        """Given the information about the detector, return an array with values of
           either 1 or 0 for where the CCDs are
        """

        # find the minimum pixel scale and set the number of pixels
        xps = self.xbin * self.pix_size
        yps = self.ybin * self.pix_size
        pw = round(self.width / xps)
        ph = round(self.height / yps)

        # create the array
        arr = np.zeros((ph, pw), dtype=float)
        y, x = np.indices((ph, pw))

        # set up where the detectors are
        for ccd in self.detector:
            x1, x2, y1, y2 = ccd.find_corners()
            x1 = (x1 + 0.5 * self.width) / xps
            x2 = (x2 + 0.5 * self.width) / xps
            y1 = (y1 + 0.5 * self.height) / yps
            y2 = (y2 + 0.5 * self.height) / yps
            mask = (x1 <= x) * (x < x2) * (y1 <= y) * (y < y2)
            arr[mask] = 1

        return arr
