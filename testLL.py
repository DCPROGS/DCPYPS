#! /usr/bin/python

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

__author__ = "User"
__date__ = "$25-Oct-2016 20:05:02$"

import numpy as np
from nose.tools import assert_almost_equal

if __name__ == "__main__":
    from dcpyps.dcfits.stats import lik_intervals
    from dcpyps.dcfits.equations import SSR
    from dcpyps.dcfits.simplex import Simplex
    X = np.array([1., 1., 2., 2., 3., 3., 4., 4., 5., 5.])
    Y = np.array([3.17, 13.25, 19.8, 14.18, 11.43, 25.85, 13.81, 25.49,
                  26.94, 38.86])
    W = np.ones((10, ))
    from dcpyps.dcfits.equations import Linear
    equation = Linear('Linear')
    pars = np.array([3., 5.])
    ssr = SSR(equation, (X, Y, W))
    simp = Simplex(ssr.equation, pars)
    result = simp.run()
    assert_almost_equal( result.rd['x'][0], 3.66600016, delta = 1e-8)
    assert_almost_equal( result.rd['x'][1], 5.20399993, delta = 1e-8)
    assert_almost_equal( result.rd['fval'], 395.86544, delta = 1e-8)
    
    maxlik = ssr.loglik(result.rd['x'])
    assert_almost_equal(maxlik, 32.581831644727835, delta = 1e-12)
    
    from dcpyps.dcfits.stats import ObservedInformation
    errs = ObservedInformation(result.rd['x'], ssr.equation, X.size)
    assert_almost_equal(errs.approximateSD[0], 5.2168715720984222, delta = 1e-9)
    assert_almost_equal(errs.approximateSD[1], 1.5729459621937478, delta = 1e-9)

    from dcpyps.dcfits.stats import LikelihoodIntervals
    li = LikelihoodIntervals(result.rd['x'], ssr, errs.approximateSD, X.size)
    print(li.limits)
#    assert_almost_equal(li.limits[0][1], 16.001557908087079, delta = 1e-9)
#    assert_almost_equal(li.limits[1][0], 1.4776049294397993, delta = 1e-9)
#    assert_almost_equal(li.limits[1][1], 8.9233105381283568, delta = 1e-9)
    
    Llimits = lik_intervals(result.rd['x'], errs.approximateSD, equation, (X, Y, W))
    assert_almost_equal(Llimits[0][1], 16.001557908087079, delta = 1e-9)
    assert_almost_equal(Llimits[1][0], 1.4776049294397993, delta = 1e-9)
    assert_almost_equal(Llimits[1][1], 8.9233105381283568, delta = 1e-9)
