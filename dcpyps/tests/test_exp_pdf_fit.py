#! /usr/bin/python

__author__ = "Remis"
__date__ = "$14-Sep-2016 19:59:09$"

import numpy as np
from nose.tools import assert_almost_equal

def test_load_patches():
    from dcpyps.dcfits.expPDF import load_data, extract_intervals
    concentration = 0.1
    dir = './dcpyps/samples/etc/EKDIST_patches/EKDIST_patches.csv'
    lim = 0.1
    recs, intervals, tres = load_data(dir, concentration, is_open=False, limit=lim)
    #shints = extract_intervals(recs, lim, 'shut')
    assert len(recs) == 2
    assert recs[1].filenames == ['./dcpyps/samples/etc/EKDIST_patches/15071504.scn']
    assert_almost_equal( intervals[0][0], 0.00050565, delta = 1e-8)

def test_ExponentialPDF_class():
    from dcpyps.dcfits.equations import ExponentialPDF