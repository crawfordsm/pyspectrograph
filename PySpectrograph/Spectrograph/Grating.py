

class Grating:

    """A class that describing gratings.  Sigma should be in lines/mm and the
       units of the dimensions should be mm.
    """

    def __init__(self, name='', spacing=600, order=1, height=100, width=100,
                 thickness=100, blaze=0, type='transmission'):
        # define the variables that describe the grating
        self.order = order
        self.height = height
        self.width = width
        self.thickness = thickness
        self.sigma = 1.0 / spacing
        self.blaze = blaze
        self.name = name
        self.type = type
        # set the sign for the grating equation
        self.sign = 1
        if self.type == 'transmission':
            self.sign = -1
