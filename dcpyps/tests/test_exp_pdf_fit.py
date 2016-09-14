#! /usr/bin/python

__author__ = "Remis"
__date__ = "$14-Sep-2016 19:59:09$"

import numpy as np
from nose.tools import assert_almost_equal

def test_load_patches():
    from dcpyps.dcfits.expPDF import load_patches, extract_intervals
    concentration = 0.1
    dir = '../samples/etc/EKDIST_patches'
    recs = load_patches(dir, concentration)
    lim = 0.1
    shints = extract_intervals(recs, lim, 'shut')
    assert len(recs) == 2
    assert recs[1].filenames == ['../samples/etc/EKDIST_patches/15071504.scn']
    assert_almost_equal( shints[0][0], 0.00050565, delta = 1e-8)

