# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np

import pytest

from ..Optics import Optics


def test_make_optic():
    op = Optics(diameter=100, focallength=100, width=100, zpos=0, focus=0)
    
    assert op.diameter == 100
    assert op.focallength == 100
    assert op.width == 100
    assert op.focus == 0
    assert op.zpos == 0
