import copy
from math import fabs, sqrt
from scipy import optimize
import numpy as np
from numpy import linalg as nplin

class ObservedInformation(object):
    """
    Approximate standard deviations of estimated values are calculated using
    the observed information matrix. The observed information, or observed
    Fisher information, is the negative of the second derivative (the Hessian
    matrix) of the "log-likelihood".
    """
    def __init__(self, theta, equation, func, *args):
        self.theta = theta
        self.eq = equation
        self.func = func
        self.args = args
        self.__calculate_variances()

    def __calculate_variances(self):
        """ """
        self.hessian = self.__calculate_hessian()
        weight = self.func(self.theta) / (self.eq.nfit - self.eq.kfit)
        self.covariance = nplin.inv(0.5 * self.hessian) * weight 
        self.approximateSD = np.sqrt(self.covariance.diagonal())
        self.correlations = self.__calculate_correlations(self.covariance)
        
    def __calculate_correlations(self, covar):
        """
        Compute correlation matrix.
        """
        correl = np.zeros((len(covar),len(covar)))
        for i1 in range(len(covar)):
            for j1 in range(len(covar)):
                correl[i1,j1] = (covar[i1,j1] /
                    np.sqrt(np.multiply(covar[i1,i1],covar[j1,j1])))
        return correl

    def __calculate_hessian(self):
        """
        Compute Hessian matrix.
        """
        hess = np.zeros((self.eq.kfit, self.eq.kfit))
        deltas = self.__optimal_deltas()
        # Diagonal elements of Hessian
        coe11 = np.array([self.theta.copy(), ] * self.eq.kfit) + np.diag(deltas)
        coe33 = np.array([self.theta.copy(), ] * self.eq.kfit) - np.diag(deltas)
        for i in range(self.eq.kfit):
            hess[i, i] = ((self.func(coe11[i]) - 
                2.0 * self.func(self.theta) +
                self.func(coe33[i])) / (deltas[i]  ** 2))
        # Non diagonal elements of Hessian
        for i in range(self.eq.kfit):
            for j in range(self.eq.kfit):
                coe1, coe2 = self.theta.copy(), self.theta.copy()
                coe3, coe4 = self.theta.copy(), self.theta.copy()
                if i != j:                
                    coe1[i] += deltas[i]
                    coe1[j] += deltas[j]
                    coe2[i] += deltas[i]
                    coe2[j] -= deltas[j]
                    coe3[i] -= deltas[i]
                    coe3[j] += deltas[j]
                    coe4[i] -= deltas[i]
                    coe4[j] -= deltas[j]
                    hess[i, j] = ((
                        self.func(coe1) - self.func(coe2) -
                        self.func(coe3) + self.func(coe4)) /
                        (4 * deltas[i] * deltas[j]))
        return hess

    def __optimal_deltas(self):
        """ """
        Lcrit = (self.func(self.theta) + fabs(self.func(self.theta)*0.005))
        deltas = 0.001 * self.theta
        L = self.func(self.theta + deltas)
        if L < Lcrit:
            count = 0
            while L < Lcrit and count < 100:
                L, deltas, count = self.__change_L(2, deltas, count)
        elif L > Lcrit:
            count = 0
            while L > Lcrit and count < 100:
                L, deltas, count = self.__change_L(0.5, deltas, count)
        return deltas

    def __change_L(self, factor, deltas, count):
        """ """
        new_deltas = deltas * factor
        L = self.func(self.theta + new_deltas)
        return L, new_deltas, count+1
    
    def print_errs(self):
        """ """
        for par, sd in zip(self.theta, self.approximateSD):
            print('Par: {0:.6f}; approximate SD = {1:.6f}'.format(par, sd))


class LikelihoodIntervals(object):
    """The most likely values of the parameters are those that correspond to
    maximal log-likelihood, Lmax. The likelihood intervals for for a particular
    parameter are estimated by holding that parameter at a fixed value away
    from Lmax and repeating the fit with all other parameters free to vary. 
    Log-likelihood value so achieved, Lcrit, is less than Lmax. Lcrit is 
    obtained for two values of the parameter, one above and one below its 
    maximum likelihood estimate. These upper and lower limits are usually 
    asymmetrical to the best fit value. 
    """
    def __init__(self, theta, equation, SD):
        self.theta = theta
        self.eq = equation
        self.SD = SD
        self.limits = self.lik_intervals()   
                 
    def lik_intervals(self):
        """
        Calculate upper and lower limits for all free parameters.
        """
        self.__max_crit_likelihoods()
        Llimits = []
        i = 0
        for j in range(len(self.theta)):
            if not self.eq.fixed[j]:
                xhi1, xlo1, xlo2, xhi2 = self.__initial_guesses(i)
                xlowlim = self.__bisect_limit(j, xlo1, xhi1)
                xhighlim = self.__bisect_limit(j, xhi2, xlo2)
                Llimits.append([xlowlim, xhighlim])
                i += 1
        return Llimits

    def __max_crit_likelihoods(self):
        """Set critical likelihood at 0.5-unit interval (corresponding to one
        standard error).
        """
        #TODO: take input to how specify m
        self.m = tvalue(self.eq.nfit - self.eq.kfit)**2 / 2.0
        self.Lmax = -self.eq.loglik(self.theta)
        self.clim = sqrt(2. * self.m)    # equivalent to SD for Gaussian
        self.Lcrit = self.Lmax - self.m

    def __initial_guesses(self, i):
        """
        Get initial limits for bisection. 
        """
        xhi1 = self.theta[i]
        #TODO: if parameter constrained to be positive- ensure that xlow is positive
        xlo1 = self.theta[i] - 2 * self.clim * self.SD[i]
        xlo2 = xhi1
        xhi2 = self.theta[i] + 5 * self.clim * self.SD[i]
        return xhi1, xlo1, xlo2, xhi2

    def __bisect_limit(self, j, lo1, hi1):
        """ """
        xlo1, xhi1 = lo1, hi1
        xlowlim = None
        found = False
        iter = 0
        while not found and iter < 100: 
            xav1 = (xlo1 + xhi1) / 2
            L = -self.__lik_contour(xav1, j) 
            if fabs(self.Lcrit - L) > 0.01:
                if L < self.Lcrit:
                    xlo1 = xav1
                else:
                    xhi1 = xav1
            else:
                found = True
                xlowlim = xav1
                if xlowlim < 0:
                    xlowlim = None
            iter += 1
        return xlowlim

    def __lik_contour(self, x, num):
        """ """
        functemp = copy.deepcopy(self.eq)
        functemp.fixed[num] = True
        functemp.pars[num] = x
        theta = functemp.theta
        result = optimize.minimize(functemp.loglik, theta, method='Nelder-Mead')
        return result.fun
    
    def print_intervals(self):
        """ """
        for par, sd, lims in zip(self.theta, self.SD, self.limits):
            print('Par: {0:.6f}; approximate SD = {1:.6f}: '.format(par, sd) +
                'lower limit = {0:.6f}; upper limit = {1:.6f}'.
                format(lims[0], lims[1]))
         

def tvalue(ndf):
    """
    Return P=0.95 value of Student's t, with ndf degrees of freedom.
    """
    ttable = [12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365,
        2.306, 2.262, 2.228, 2.201, 2.179, 2.160, 2.145, 2.131,
        2.120, 2.210, 2.101, 2.093, 2.086, 2.080, 2.074, 2.069,
        2.064, 2.060, 2.056, 2.052, 2.048, 2.045, 2.042]
    if ndf in range(1, 31):
        tval = ttable[ndf-1]
    elif ndf in range(31, 41):
        frac = (ndf - 30.0) / 10.0 #=0.1 to 1
        tval = 2.042 - frac * (2.042 - 2.021)
    elif ndf in range(41, 61):
        frac = (ndf - 40.0) / 20.0 #=0.05 to 1
        tval = 2.021 - frac * (2.021 - 2.000)
    elif ndf in range(61, 121):
        frac = (ndf - 60.0) / 60.0 #=1/60 to 1
        tval = 2.000 - frac * (2.000 - 1.980)
    elif ndf > 120:
        tval = 1.96
    else:
        print(' ERROR IN TVALUE ')
    return tval           
