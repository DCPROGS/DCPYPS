import sys
from math import sqrt#, fabs
import numpy as np

from dcpyps.dcfits.simplex import Simplex
from dcpyps.dcfits.stats import SSD, SSDlik, lik_intervals
from dcpyps.dcfits.stats import covariance_matrix
from dcpyps.dcfits.stats import correlation_matrix
from dcpyps.dcfits.stats import approximateSD
from dcpyps.dcfits.stats import tvalue

class SingleFitSession(object):
    def __init__(self, dataset, equation):
        """
        """
        self.dataset = dataset
        self.XYW = (self.dataset.X, self.dataset.Y, self.dataset.W)
        self.equation = equation
        self.equation.propose_guesses(self.dataset)
        
    def fit(self):
        simp = Simplex(SSD, self.equation.theta, *(self.equation.to_fit, self.XYW))
        result = simp.run()
        self.equation.theta = result.rd['x']
       
    def calculate_errors(self):        
        ASD = ApproximateSD(self.equation, self.dataset)
        
        self.Smin, self.equation.theta = SSD(self.equation.theta, (self.equation.to_fit,
            self.XYW))
        self.kfit = len(np.nonzero(np.invert(self.equation.fixed))[0])
        self.ndf = self.dataset.size() - self.kfit
        self.var = self.Smin / self.ndf
        self.Sres, self.Lmax = sqrt(self.var), -SSDlik(self.equation.theta, self.equation.to_fit, self.XYW)


        self.m = tvalue(self.ndf) ** 2 / 2.0
        self.clim = sqrt(2. * self.m)
        self.Lcrit = self.Lmax - self.m
        
        self.Llimits = lik_intervals(self.equation.theta, ASD.approximateSD, self.equation, self.XYW)
        
        #self.output.write(self.string_estimates())
        #self.output.write(self.string_liklimits())
        
    def string_estimates(self):
        j = 0
        str = 'Number of point fitted = {0:d}'.format(self.dataset.size())
        
        str += '\nNumber of parameters estimated = {0:d}'.format(self.kfit)
        str += '\nDegrees of freedom = {0:d}'.format(self.ndf)
        
        str += ('\nResidual error SD = {0:.3f}      (variance = {1:.3f})'.
            format(self.Sres, self.var))
        
        for i in range(len(self.eq.names)):
            str += '\nParameter {0:d}: {1}  \t= {2:.6g}  \t'.format(i+1, self.eq.names[i], self.eq.pars[i])
            if not self.eq.fixed[i]:
                str += '  Approx SD = {0:.6g}\t'.format(self.aproxSD[j])
                str += '  CV = {0:.1f}'.format(self.CVs[j])
                j += 1
            else:
                str += '  (fixed)'

        str += ('\nMinimum SSD = {0:.3f}; \nMax log-likelihood = {1:.3f}'.
            format(self.Smin, self.Lmax))
        str += ('\nCorrelation matrix = ' + 
            '[!!!! PRINTOUT OF CORRELATION MATRIX NOT IMPLEMENTED YET. SORRY.\n')
#        self.output.write(correl)
        if np.any(np.absolute(self.correl - np.identity(self.kfit)) > 0.9):
            str += ("\nWARNING: SOME PARAMETERS ARE STRONGLY CORRELATED (coeff > 0.9); try different guesses")

        if np.any(self.CVs > 33):
            str += "\nWARNING: SOME PARAMETERS POORLY DEFINED (CV > 33%); try different guesses"
        return str

        
    def string_liklimits(self):
        j = 0
        str = '\nLIKELIHOOD INTERVALS\n'
        str += ('{0:.3g}-unit Likelihood Intervals'.format(self.m) +
            '  (equivalent SD for Gaussian- {0:.3g})'.format(self.clim))
        str += '\nLmax= {0:.6g};   Lcrit= {1:.6g}'.format(self.Lmax, self.Lcrit)
        for i in range(len(self.eq.names)):
            str += '\nParameter {0:d}:   {1}\t= {2:.6g}'.format(i+1, self.eq.names[i], self.eq.pars[i])
            if not self.eq.fixed[i]:
                try:
                    str += '\t  LOWER = {0:.6g}'.format(self.Llimits[j][0])
                except:
                    str += '\t  LOWER limit not found'
                try:
                    str += '\t  UPPER = {0:.6g}'.format(self.Llimits[j][1])
                except:
                    str += '\t  UPPER limit not found'
                j += 1
            else:
                str += '\t  (fixed)'
        return str
    
class ApproximateSD(object):
    def __init__(self, equation, dataset):
        self.theta = equation.theta
        self.equation = equation
        self.func = equation.to_fit
        self.dataset = dataset
        self.XYW = (dataset.X, dataset.Y, dataset.W)
        self.__calculate_all()
        
    def __calculate_all(self):
        self.covariance = covariance_matrix(self.theta, self.func, self.XYW)
        self.correlations = correlation_matrix(self.covariance)
        self.approximateSD = approximateSD(self.theta, self.func, self.XYW)
        self.CVs = 100.0 * self.approximateSD / self.theta
        
    def __optimal_deltas(self):
        """ """

        Lcrit = 1.005 * -0.5 * SSD(self.theta, (self.func, self.XYW))[0]
        deltas = 0.1 * self.theta
        L = -0.5 * SSD(self.theta + deltas, (self.func, self.XYW))[0]
        if L > Lcrit:
            count = 0
            while L > Lcrit and count < 100:
                deltas *= 2
                L = -0.5 * SSD(self.theta + deltas, (self.func, self.XYW))[0]
                count += 1
        elif L < Lcrit:
            count = 0
            while L < Lcrit and count < 100:
                deltas *= 0.5
                L = -0.5 * SSD(self.theta + deltas, (self.func, self.XYW))[0]
                count += 1
        return deltas

    def __hessian(self):
        """
        """
        hessian = np.zeros((self.theta.size, self.theta.size))
        deltas = self.__optimal_deltas()
        # Diagonal elements of Hessian
        coe11 = np.array([self.theta.copy(), ] * self.theta.size) + np.diag(deltas)
        coe33 = np.array([self.theta.copy(), ] * self.theta.size) - np.diag(deltas)
        for i in range(self.theta.size):
            hessian[i, i] = ((SSD(coe11[i], (self.func, self.XYW))[0] - 
                2.0 * SSD(self.theta, (self.func, self.XYW))[0] +
                SSD(coe33[i], (self.func, self.XYW))[0]) / (deltas[i]  ** 2))
        # Non diagonal elements of Hessian
        for i in range(self.theta.size):
            for j in range(self.theta.size):
                coe1, coe2, coe3, coe4 = self.theta.copy(), self.theta.copy(), self.theta.copy(), self.theta.copy()
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
                        SSD(coe1, (self.func, self.XYW))[0] -
                        SSD(coe2, (self.func, self.XYW))[0] -
                        SSD(coe3, (self.func, self.XYW))[0] +
                        SSD(coe4, (self.func, self.XYW))[0]) /
                        (4 * deltas[i] * deltas[j]))
        return 0.5 * hessian

    def covariance_matrix(self):
        """ """
        cov = nplin.inv(self.__hessian(self.theta, self.func, self.XYW))
        if self.dataset.weightmode == 1:
            errvar = SSD(self.theta, (self.func, self.XYW))[0] / (self.XYW[0].size - self.theta.size)
        else:
            errvar = 1.0
        return cov * errvar

    def approximateSD(self):
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
        cov = self.covariance_matrix(self.theta, self.func, self.XYW)
        return np.sqrt(cov.diagonal())

    def correlation_matrix(self, covar):
        correl = np.zeros((len(covar),len(covar)))
        for i1 in range(len(covar)):
            for j1 in range(len(covar)):
                correl[i1,j1] = (covar[i1,j1] / 
                    np.sqrt(np.multiply(covar[i1,i1],covar[j1,j1])))
        return correl

class LikelihoodIntervals(object):
    def __init__(self):
        pass

class MultipleFitSession(object):
    def __init__(self, output=sys.stdout):
        self.output = output
        self.list = []
        self.pooled = None
        
    def add(self, fitset):
        self.list.append(fitset)
        
    def pool(self, norm=False, output=sys.stdout):
        dataset = data.XYDataSet()
        for session in self.list:
            if norm:
                dataset.pool(session.data.X, session.data.normY, session.data.normS)
            else:
                dataset.pool(session.data.X, session.data.Y, session.data.S)
        dataset.weightmode = 2
        self.pooled = SingleFitSession(dataset, self.list[0].eq, output)
        
    def string_average_estimates(self):
        str = 'Average of estimates of {0:d} sets:'.format(len(self.list))
       
        for i in range(len(self.list[0].eq.names)):
            pars = []
            for j in range(len(self.list)):
                pars.append(self.list[j].eq.pars[i])
            str += ('\nParameter {0:d}: {1}  \t= {2:.6g} +/- {3:.6g}'.
                format(i+1, self.list[0].eq.names[i], np.mean(pars), 
                np.std(pars)/sqrt(len(pars))))
        return str
    
    def prepare_fplot(self, type):
        plot_list = []
        if type == 'pooled':
            X, Y = self.pooled.eq.calculate_plot(self.pooled.data.X, self.pooled.eq.pars)
            plot_list.append([X,Y])
        else:
            for session in self.list:
                if type == 'fit':
                    pars = session.eq.pars
                elif type == 'guess':
                    pars = session.eq.guess
                elif type == 'norm':
                    pars = session.eq.normpars
                X, Y = session.eq.calculate_plot(session.data.X, pars)
                plot_list.append([X,Y])
        return plot_list
    




