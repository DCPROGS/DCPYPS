#! /usr/bin/python

__author__ = "Remis"
__date__ = "$22-Mar-2016 14:27:44$"

import numpy as np
from nose.tools import assert_almost_equal

def test_compare_scipy_to_dc_simplex():
    # scipy Nelder-Mead
    from dcpyps.dcfits.simplex import rosen, rosen1, Simplex
    from scipy.optimize import minimize
    x0 = np.array([1.3, 0.7, 0.8, 1.9, 1.2])
    res = minimize(rosen, x0, method='nelder-mead',
                options={'xtol': 1e-8, 'disp': True})
    # dc simplex            
    kwargs = {'parerror' : 1e-8, 'funerror' : 1e-8}
    x0 = np.array([1.3, 0.7, 0.8, 1.9, 1.2])
    simp = Simplex(rosen1, x0, **kwargs)
    result = simp.run()
    np.testing.assert_almost_equal(result.rd['x'], res.x, 8) 
    assert result.rd['fval'] < res.fun
        
def test_simplex_iteration_number():
    from dcpyps.dcfits.simplex import Simplex, rosen1
    x0 = np.array([1.3, 0.7, 0.8, 1.9, 1.2])
    kwargs = {'maxiter' : 10, 'maxevaluatons' : 100}
    simp = Simplex(rosen1, x0, **kwargs)
    assert kwargs['maxiter'] == simp.settings['maxiter']
    result = simp.run()
    assert result.rd['nit'] == kwargs['maxiter']
    
def test_simplex_kwargs():
    kwargs = {'step' : -1, 'reflect' : -1, 'extend' : -1}
    from dcpyps.dcfits.simplex import Simplex, rosen1
    x0 = np.array([1.3, 0.7, 0.8, 1.9, 1.2])
    simp = Simplex(rosen1, x0, **kwargs)
    for o in kwargs.keys():
        assert kwargs[o] == simp.settings[o]

def test_regression_stats():
    from dcpyps.dcfits.stats import SSD, SSDlik, approximateSD, lik_intervals
    from dcpyps.dcfits.simplex import Simplex
    X = np.array([1., 1., 2., 2., 3., 3., 4., 4., 5., 5.])
    Y = np.array([3.17, 13.25, 19.8, 14.18, 11.43, 25.85, 13.81, 25.49,
                  26.94, 38.86])
    W = np.ones((10, ))
    from dcpyps.dcfits.equations import Linear
    equation = Linear('Linear')
    pars = np.array([3., 5.])
    simp = Simplex(SSD, pars, *(equation.to_fit, (X, Y, W)))
    result = simp.run()
    assert_almost_equal( result.rd['x'][0], 3.66600016, delta = 1e-8)
    assert_almost_equal( result.rd['x'][1], 5.20399993, delta = 1e-8)
    assert_almost_equal( result.rd['fval'], 395.86544, delta = 1e-8)
    
    maxlik = SSDlik(result.rd['x'], equation.to_fit, (X, Y, W))
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
    
def test_SingleFitSession():
    X = np.array([1., 1., 2., 2., 3., 3., 4., 4., 5., 5.])
    Y = np.array([3.17, 13.25, 19.8, 14.18, 11.43, 25.85, 13.81, 25.49,
                  26.94, 38.86])
    from dcpyps.dcfits.data import XYDataSet
    set0 = XYDataSet()
    set0.from_columns(X, Y)
    from dcpyps.dcfits.equations import Linear
    equation = Linear('Linear')
    from dcpyps.dcfits.fitting import SingleFitSession
    fsession = SingleFitSession(set0, equation)
    fsession.fit()
    fsession.calculate_errors()
    assert_almost_equal(fsession.Llimits[0][1], 16.001557908087079, delta = 1e-6)
    assert_almost_equal(fsession.Llimits[1][0], 1.4776049294397993, delta = 1e-6)
    assert_almost_equal(fsession.Llimits[1][1], 8.9233105381283568, delta = 1e-6)