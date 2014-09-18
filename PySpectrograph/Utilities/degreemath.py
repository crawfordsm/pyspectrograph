import math
import numpy


def sind(x):
    """Return the sin of x where x is in degrees"""
    if isinstance(x, numpy.ndarray):
        return numpy.sin(math.pi * x / 180.0)
    return math.sin(math.radians(x))


def cosd(x):
    """Return the cos of x where x is in degrees"""
    if isinstance(x, numpy.ndarray):
        return numpy.cos(math.pi * x / 180.0)
    return math.cos(math.radians(x))


def tand(x):
    """Return the cos of x where x is in degrees"""
    if isinstance(x, numpy.ndarray):
        return numpy.tan(math.pi * x / 180.0)
    return math.cos(math.radians(x))
