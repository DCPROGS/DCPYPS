import copy
import math
from scipy import optimize
import numpy as np
from numpy import linalg as nplin

class Variation():
    def __init__(self, XYW, func, pars):
        '''
        Calculate measures of variation.

        Parameters
        ----------
        X, Y, W : ndarrays, shape (n, )
            x, y and weights of data points.
        func : callable func(x, args)
            The objective function to be minimized.
        args : array_like, shape (k, )
            Function parameters.
        '''
        self.func = func
        self.X, self.Y, self.W = XYW
        self.args = pars
        
    def residuals(self):
        '''
        Calculate the weighted residuals.

        Returns
        -------
        residuals : ndarray, shape (n, )
            Weighted sum of squared deviations.
        '''
        return self.W * (self.Y - self.func(self.X, self.args))
    
    def SSD(self):
        """
        Calculate weighted sum of squared deviations.

        Returns
        -------
        SSD : float
            Weighted sum of squared deviations.
        """

        return np.sum(self.residuals()**2)

def residuals(pars, func, X, Y, W):
    '''
    Calculate the weighted residuals.

    Parameters
    ----------
    pars : array_like, shape (k, )
        Function parameters.
    func : callable func(x, args)
        The objective function to be minimized.
    X, Y, W : ndarrays, shape (n, )
        x, y and weights of data points.

    Returns
    -------
    residuals : ndarray, shape (n, )
        Weighted sum of squared deviations.
    '''
    return W * (Y - func(pars, X))

def SSD(pars, args):
    """
    Calculate weighted sum of squared deviations.

    Parameters
    ----------
    pars : array_like, shape (k, )
        Function parameters.
    func : callable func(x, args)
        The objective function to be minimized.
    X, Y, W : ndarrays, shape (n, )
        x, y and weights of data points.

    Returns
    -------
    SSD : float
        Weighted sum of squared deviations.
    """
    func = args[0]
    X, Y, W = args[1]
    return np.sum(W * (Y - func(X, pars))**2), pars

def SSDlik(theta, func, args):
    """
    Calculate log likelihood coresponding to the sum of the squared deviations
    assuming that errors follow Gausian distribution.

    Parameters
    ----------
    theta : array_like, shape (k, )
        Initial guess.
    func : callable func(x, args)
        The objective function to be minimized.
    set : cvfit.data.XYDataSet type object
        Extra argument passed to func.

    Returns
    -------
    SSDlik : float
        Log likelihood of SSD function.
    """
    X, Y, W = args[0], args[1], args[2]
    S, theta = SSD(theta, (func, (X, Y, W)))
    Sres = math.sqrt(S / float(X.size)) # - len(func.fixed)))
    return X.size * math.log(math.sqrt(2 * math.pi) * Sres) + S / (2 * Sres**2)

def SSDlik_contour(x, num, theta, func, set):
    functemp = copy.deepcopy(func)
    functemp.fixed[num] = True
    functemp.pars[num] = x
    theta = functemp.theta
    result = optimize.minimize(SSDlik, theta, args=(functemp.to_fit, set), 
        method='Nelder-Mead', jac=None, hess=None)
    return -result.fun

def optimal_deltas(theta, func, args):
    """ """

    Lcrit = 1.005 * -0.5 * SSD(theta, (func, args))[0]
    deltas = 0.1 * theta
    L = -0.5 * SSD(theta + deltas, (func, args))[0]
    if L > Lcrit:
        count = 0
        while L > Lcrit and count < 100:
            deltas *= 2
            L = -0.5 * SSD(theta + deltas, (func, args))[0]
            count += 1
    elif L < Lcrit:
        count = 0
        while L < Lcrit and count < 100:
            deltas *= 0.5
            L = -0.5 * SSD(theta + deltas, (func, args))[0]
            count += 1
    return deltas

def hessian(theta, func, args):
    """
    """
    hessian = np.zeros((theta.size, theta.size))
    deltas = optimal_deltas(theta, func, args)
    # Diagonal elements of Hessian
    coe11 = np.array([theta.copy(), ] * theta.size) + np.diag(deltas)
    coe33 = np.array([theta.copy(), ] * theta.size) - np.diag(deltas)
    for i in range(theta.size):
        hessian[i, i] = ((SSD(coe11[i], (func, args))[0] - 
            2.0 * SSD(theta, (func, args))[0] +
            SSD(coe33[i], (func, args))[0]) / (deltas[i]  ** 2))
    # Non diagonal elements of Hessian
    for i in range(theta.size):
        for j in range(theta.size):
            coe1, coe2, coe3, coe4 = theta.copy(), theta.copy(), theta.copy(), theta.copy()
            if i != j:                
                coe1[i] += deltas[i]
                coe1[j] += deltas[j]
                coe2[i] += deltas[i]
                coe2[j] -= deltas[j]
                coe3[i] -= deltas[i]
                coe3[j] += deltas[j]
                coe4[i] -= deltas[i]
                coe4[j] -= deltas[j]
                hessian[i, j] = ((
                    SSD(coe1, (func, args))[0] -
                    SSD(coe2, (func, args))[0] -
                    SSD(coe3, (func, args))[0] +
                    SSD(coe4, (func, args))[0]) /
                    (4 * deltas[i] * deltas[j]))
    return 0.5 * hessian

def covariance_matrix(theta, func, args, weightmode=1):
    """ """
    cov = nplin.inv(hessian(theta, func, args))
    if weightmode == 1:
        errvar = SSD(theta, (func, args))[0] / (args[0].size - theta.size)
    else:
        errvar = 1.0
    return cov * errvar

def approximateSD(theta, func, args):
    """
    Calculate approximate standard deviation of the estimates from the inverse
    of the Hessian matrix ('observed information matrix').

    Parameters
    ----------
    theta : array_like, shape (k, )
        Initial guess.
    func : callable func(x, args)
        The objective function to be minimized.
    args : object
        Extra arguments passed to func.

    Returns
    -------
    approximateSD : ndarray, shape (k,)
        Approximate SD.
    """
    cov = covariance_matrix(theta, func, args)
    return np.sqrt(cov.diagonal())
    
def correlation_matrix(covar):
    correl = np.zeros((len(covar),len(covar)))
    for i1 in range(len(covar)):
        for j1 in range(len(covar)):
            correl[i1,j1] = (covar[i1,j1] / 
                np.sqrt(np.multiply(covar[i1,i1],covar[j1,j1])))
    return correl
            
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

def lik_intervals(theta, SD, func, set):
    
    ndf = set[0].size - theta.size
    m = tvalue(ndf)**2 / 2.0
    Lmax = -SSDlik(theta, func.to_fit, set)
    clim = math.sqrt(2. * m)
    Lcrit = Lmax - m
    Llimits = []
    
    i = 0
    for j in range(len(func.pars)):
        if not func.fixed[j]:
            #print('\nCalculating Lik limits for parameter- {0} = {1:.3f}'.format(func.names[j], theta[i]))
            xhigh1 = theta[i]
            #TODO: if parameter constrained to be positive- ensure that xlow is positive
            xlow1 = theta[i] - 2 * clim * SD[i]
            xlow2 = xhigh1
            xhigh2 = theta[i] + 5 * clim * SD[i]
            #print('\tInitial guesses for lower limit: {0:.3f} and {1:.3f}'.format(xlow1, xhigh1))
            #print('\tInitial guesses for higher limit: {0:.3f} and {1:.3f}'.format(xlow2, xhigh2))

            found = False
            iter = 0
            xlowlim, xhighlim = None, None
            while not found and iter < 100: 
                L = SSDlik_contour(((xlow1 + xhigh1) / 2), j, theta,
                    func, set) 
                if math.fabs(Lcrit - L) > 0.01:
                    if L < Lcrit:
                        xlow1 = (xlow1 + xhigh1) / 2
                    else:
                        xhigh1 = (xlow1 + xhigh1) / 2
                else:
                    found = True
                    xlowlim = (xlow1 + xhigh1) / 2
                    if xlowlim < 0:
                        xlowlim = None
                    #print 'lower limit found: ', xlowlim
                iter += 1
            found = False
            iter = 0   
            while not found and iter < 100: 
                #L1, L2, L3 = SSDlik_bisect(xlow2, xhigh2, j, theta, notfixed, hill_equation, dataset)
                L = SSDlik_contour(((xlow2 + xhigh2) / 2), j, theta,
                    func, set) 
                if math.fabs(Lcrit - L) > 0.01:
                    if L > Lcrit:
                        xlow2 = (xlow2 + xhigh2) / 2
                    else:
                        xhigh2 = (xlow2 + xhigh2) / 2
                else:
                    found = True
                    xhighlim = (xlow2 + xhigh2) / 2
                    if xhighlim < 0:
                        xhighlim = None
                    #print 'higher limit found: ', xhighlim
                iter += 1
            Llimits.append([xlowlim, xhighlim])
            i += 1
    return Llimits