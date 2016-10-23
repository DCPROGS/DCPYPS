#! /usr/bin/python

__author__ = "Remis"
__date__ = "$22-Oct-2016 21:08:31$"

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

def test_dc_simplex_function_returns():
    """
    Test simplex with functions which return value or value and argument.
    """
    from dcpyps.dcfits.simplex import rosen, rosen1, Simplex
    x = np.array([1.3, 0.7, 0.8, 1.9, 1.2])
    kwargs = {'parerror' : 1e-8, 'funerror' : 1e-8}
    # dc simplex with function returning single value            
    simp1 = Simplex(rosen, x, **kwargs)
    result1 = simp1.run()
    simp2 = Simplex(rosen1, x, **kwargs)
    result2 = simp2.run()
    np.testing.assert_almost_equal(result1.rd['x'], result2.rd['x'], 8) 
    assert result1.rd['fval'] == result2.rd['fval']
        
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


