
from PySpectrograph.Utilities.degreemath import cosd, sind


def n_index():
    return 1.0000


def gratingequation(sigma, order, sign, alpha, beta, gamma=0, nd=n_index):
    """Apply the grating equation to determine the wavelength
       w = sigma/m cos (gamma) * n_ind *(sin alpha +- sin beta)

       returns wavelength in mm
    """
    angle = cosd(gamma) * nd() * (sind(alpha) + sign * sind(beta))
    return sigma / order * angle


def calc_angdisp(sigma, order, beta, gamma=0):
    """Calculate the angular dispersion according to m/sigma/cos beta

       returns angular dispersion in 1/mm
    """
    return order / sigma / cosd(beta) / cosd(gamma)


def calc_lindisp(f, sigma, order, beta, gamma=0.0):
    """Calculate the linear dispersion according to f_cam * A

       return linear dispersion in mm/mm

    """
    return f * calc_angdisp(sigma, order, beta, gamma=gamma)


def calc_anamorph(alpha, beta):
    """Calculates the anamorphic magnification

       returns the anamorpic magnification
    """
    return cosd(alpha) / cosd(beta)


def calc_demagspatial(fcol, fcam):
    """Calculate the spatial demagnification

           returns the spatial demagnification
    """
    return fcol / fcam


def calc_demagspectral(fcol, fcam, alpha, beta):
    """Calculate the spectral demagnification

       returns the spectral demagnification
    """
    return calc_demagspatial(fcol, fcam) / calc_anamorph(alpha, beta)


def calc_resolelement(w, fcol, sigma, order, alpha, beta, gamma=0.0):
    """Calculate the resolution element using dw=r*w/A/fcol

       returns resolution element in mm
    """
    r = calc_anamorph(alpha, beta)
    A = calc_angdisp(sigma, order, beta, gamma=gamma)
    return _calc_resolelement(w, fcol, r, A)


def _calc_resolelement(w, fcol, r, A):
    """Calculate the resolution element using dw=r*w/A/fcol

       returns resolution element in mm
    """
    return r * w / A / fcol


def calc_resolution(w, dw):
    """Calcualte the resolution R=w/dw

       returns the resolution
    """
    return w / dw
