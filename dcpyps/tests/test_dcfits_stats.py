#! /usr/bin/python

__author__ = "Remis"
__date__ = "$22-Mar-2016 14:27:44$"

import numpy as np
from nose.tools import assert_almost_equal

def test_regression_stats():
    from dcpyps.dcfits.stats import approximateSD, lik_intervals, SSR
    from dcpyps.dcfits.simplex import Simplex
    X = np.array([1., 1., 2., 2., 3., 3., 4., 4., 5., 5.])
    Y = np.array([3.17, 13.25, 19.8, 14.18, 11.43, 25.85, 13.81, 25.49,
                  26.94, 38.86])
    W = np.ones((10, ))
    from dcpyps.dcfits.equations import Linear
    equation = Linear('Linear')
    pars = np.array([3., 5.])
    ssr = SSR(equation.to_fit, (X, Y, W))
    simp = Simplex(ssr.SSR, pars)
    result = simp.run()
    assert_almost_equal( result.rd['x'][0], 3.66600016, delta = 1e-8)
    assert_almost_equal( result.rd['x'][1], 5.20399993, delta = 1e-8)
    assert_almost_equal( result.rd['fval'], 395.86544, delta = 1e-8)
    
    maxlik = ssr.loglik(result.rd['x'])
    assert_almost_equal(maxlik, 32.581831644727835, delta = 1e-12)
    
    SD = approximateSD(result.rd['x'], equation.to_fit, (X, Y, W))
    assert_almost_equal(SD[0], 5.2168715720984222, delta = 1e-12)
    assert_almost_equal(SD[1], 1.5729459621937478, delta = 1e-12)

    Llimits = lik_intervals(result.rd['x'], SD, equation, (X, Y, W))
    assert_almost_equal(Llimits[0][1], 16.001557908087079, delta = 1e-12)
    assert_almost_equal(Llimits[1][0], 1.4776049294397993, delta = 1e-12)
    assert_almost_equal(Llimits[1][1], 8.9233105381283568, delta = 1e-12)
    
def test_stats_tvalue():
    from dcpyps.dcfits.stats import tvalue
    assert tvalue(2) == 4.303
    assert_almost_equal(tvalue(100), 1.9866666, delta = 1e-6)
    
def test_EstimateErrors():
    X = np.array([1., 1., 2., 2., 3., 3., 4., 4., 5., 5.])
    Y = np.array([3.17, 13.25, 19.8, 14.18, 11.43, 25.85, 13.81, 25.49,
                  26.94, 38.86])
    from dcpyps.dcfits.data import XYDataSet
    dataset = XYDataSet()
    dataset.from_columns(X, Y)
    from dcpyps.dcfits.equations import Linear
    equation = Linear('Linear')
    pars = np.array([3., 5.])
    from dcpyps.dcfits.simplex import Simplex
    from dcpyps.dcfits.stats import SSR
    ssr = SSR(equation.to_fit, (dataset.X, dataset.Y, dataset.W))
    simp = Simplex(ssr.SSR, pars)
    result = simp.run()
    equation.theta = result.rd['x']
    from dcpyps.dcfits.stats import EstimateErrors
    errs = EstimateErrors(equation, dataset)
    
    assert_almost_equal(errs.approximateSD[0], 5.2168715720984222, delta = 1e-12)
    assert_almost_equal(errs.approximateSD[1], 1.5729459621937478, delta = 1e-12)
    
    assert_almost_equal(errs.Llimits[0][1], 16.001557908087079, delta = 1e-6)
    assert_almost_equal(errs.Llimits[1][0], 1.4776049294397993, delta = 1e-6)
    assert_almost_equal(errs.Llimits[1][1], 8.9233105381283568, delta = 1e-6)
    
def test_SSR():
    from dcpyps.dcfits.stats import SSR
    def func(X, pars):
        return np.array([2.0, 2.0])
    X, Y, W = np.array([1.0, 1.0]), np.array([1.0, 3.0]), np.array([1.0, 1.0])
    ssr = SSR(func, (X, Y, W))
    res = ssr.residuals(None)
    assert (res[0] == -1.0) and (res[1] == 1.0)
    sr = ssr.SSR(None)
    assert sr == 2.0
