import math


class Slit:

    """A class that describing the slit.  Only assuming a single slit.  All sizes are in mm.
       All positions assume the center of the slit. Phi is in arcseconds
    """

    def __init__(self, name='', height=100, width=100, zpos=0, xpos=0, ypos=0, phi=1):

        # define the variables that describe the grating
        self.height = height
        self.width = width
        self.name = name
        self.zpos = zpos
        self.xpos = xpos
        self.ypos = ypos
        self.phi = phi

    def set_phi(self, phi):
        self.phi = phi

    def calc_phi(self, ftel):
        """Calculate phi(angle on sky) assuming w/ftel

           returns phi in arcseconds
        """
        return 3600.0 * math.degrees(self.width / ftel)

    def calc_width(self, ftel):
        """Calculate the width assuming ftel*phi(rad).

           returns the slit width in mm
        """
        return ftel * math.radians(self.phi / 3600.0)
