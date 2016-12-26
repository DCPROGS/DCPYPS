#! /usr/bin/python

__author__ = "Remis"
__date__ = "$22-Mar-2016 14:27:44$"

import numpy as np
from nose.tools import assert_almost_equal

from dcpyps.dcfits.dataIO import XYDataSet
from dcpyps.dcfits.equations import SSR
from dcpyps.dcfits.simplex import Simplex
from dcpyps.dcfits.equations import Linear
from dcpyps.dcfits.stats import ObservedInformation
from dcpyps.dcfits.stats import LikelihoodIntervals

def test_regression_stats():
    
    X = np.array([1., 1., 2., 2., 3., 3., 4., 4., 5., 5.])
    Y = np.array([3.17, 13.25, 19.8, 14.18, 11.43, 25.85, 13.81, 25.49,
                  26.94, 38.86])
    W = np.ones((10, ))
    dataset = XYDataSet()
    dataset.from_columns(X, Y)
    equation = Linear('Linear')
    pars = np.array([3., 5.])
    ssr = SSR(equation, dataset)
    simp = Simplex(ssr.equation, pars)
    result = simp.run()
    assert_almost_equal( result.rd['x'][0], 3.66600016, delta = 1e-8)
    assert_almost_equal( result.rd['x'][1], 5.20399993, delta = 1e-8)
    assert_almost_equal( result.rd['fval'], 395.86544, delta = 1e-8)
    
    maxlik = ssr.loglik(result.rd['x'])
    assert_almost_equal(maxlik, 32.581831644727835, delta = 1e-12)
    
    errs = ObservedInformation(result.rd['x'], ssr, ssr.to_fit)
    assert_almost_equal(errs.approximateSD[0], 5.2168715720984222, delta = 1e-9)
    assert_almost_equal(errs.approximateSD[1], 1.5729459621937478, delta = 1e-9)

    ssr.pars = result.rd['x']
    #ll = LikelihoodIntervals(ssr, errs.approximateSD)
    ll = LikelihoodIntervals(result.rd['x'], ssr, errs.approximateSD)
    assert_almost_equal(ll.limits[0][1], 16.001557908087079, delta = 1e-9)
    assert_almost_equal(ll.limits[1][0], 1.4776049294397993, delta = 1e-9)
    assert_almost_equal(ll.limits[1][1], 8.9233105381283568, delta = 1e-9)
       
def test_stats_tvalue():
    from dcpyps.dcfits.stats import tvalue
    assert tvalue(2) == 4.303
    assert_almost_equal(tvalue(100), 1.9866666, delta = 1e-6)
        
#def test_SSR():
#    from dcpyps.dcfits.equations import SSR
#    def func(X, pars):
#        return np.array([2.0, 2.0])
#    X, Y, W = np.array([1.0, 1.0]), np.array([1.0, 3.0]), np.array([1.0, 1.0])
#    ssr = SSR(func, (X, Y, W))
#    res = ssr.residuals(None)
#    assert (res[0] == -1.0) and (res[1] == 1.0)
#    sr = ssr.equation(None)
#    assert sr == 2.0
