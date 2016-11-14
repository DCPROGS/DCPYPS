#! /usr/bin/python

__author__ = "Remis"
__date__ = "$14-Sep-2016 19:59:09$"

import numpy as np
from scipy.optimize import minimize
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

def test_MultiExponentialPDF_class():
    from dcpyps.dcfits.equations import MultiExponentialPDF
    from dcpyps import dcio
    # LOAD DATA.
    intervals = np.array(dcio.txt_load_one_col("./dcpyps/samples/etc/intervals.txt"))
    taus, areas = [0.036, 1.1], [0.20]
    expPDF = MultiExponentialPDF(intervals, taus, areas)
    assert len(expPDF.areas) == 2
    assert len(expPDF.fixed) == 4
    
    theta = [0.036, 1.1, 0.20]
    start_lik = expPDF.loglik(theta)
    assert_almost_equal( start_lik, 87.31806715582867, delta = 1e-8)
    res = minimize(expPDF.loglik, theta, method='Nelder-Mead')
    assert_almost_equal(res.fun, 87.288287733277258, delta = 1e-8)
    assert len(res.x) == 3
    expPDF.theta = res.x
    
    from dcpyps.dcfits.stats import ObservedInformation
    taus, areas = res.x[:2], [res.x[-1]]
    #expPDFlik = MultiExponentialPDFLogLik(intervals, taus, areas)
    errs = ObservedInformation(res.x, expPDF, expPDF.loglik)
    assert_almost_equal(errs.approximateSD[0], 0.014062302913235354, delta = 1e-6)
    assert_almost_equal(errs.approximateSD[1], 0.15042035132630516, delta = 1e-6)
    assert_almost_equal(errs.approximateSD[2], 0.05795821682555078, delta = 1e-6)
    
    #from dcpyps.dcfits.stats import LikelihoodIntervals
    #ll = LikelihoodIntervals(res.x, expPDF, errs.approximateSD)
    #assert_almost_equal(ll.limits[0][1], 16.001557908087079, delta = 1e-9)
    #assert_almost_equal(ll.limits[1][0], 1.4776049294397993, delta = 1e-9)
    #assert_almost_equal(ll.limits[1][1], 8.9233105381283568, delta = 1e-9)
    