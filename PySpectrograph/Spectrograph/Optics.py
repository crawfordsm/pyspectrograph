

class Optics:

    """A class that describing optics.  All dimensions should in mm.  This assumes all optics
       can be desribed by a diameter and focal length. zpos is the distance in mm that the center
       of the optic is from the primary mirror.   focas is the offset from that position
    """

    def __init__(self, name='', diameter=100, focallength=100, width=100,
                 zpos=0, focus=0):
        # define the variables that describe the opticsg
        self.diameter = diameter
        self.focallength = focallength
        self.width = width
        self.name = name
        # set distances of the optics
        self.zpos = zpos
        self.focus = focus

    def speed(self):
        """Camera Speed f/ = d/f
        """
        return self.diameter / self.focallength

    def platescale(self):
        """Plate scale for a given optic in arcsec/mm
        """
        return 206265 / self.focallength

    def position(self):
        """Current position along the optical path of the element
        """
        return self.zpos + self.focus
