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
    def __init__(self, dataset, equation, output=sys.stdout):
        """
        """
        self.output = output
        self.data = dataset
        self.eq = equation
        self.eq.propose_guesses(self.data)
        self.text = '\n\tFitting session for ' + self.data.title + ' initialised...'
        self.output.write(self.text)
        
    def fit(self):
        # Least square fitting
        #coeffs, cov, dict, mesg, ier = optimize.leastsq(residuals, self.eq.theta,
        #    args=(self.eq, self.data.X, self.data.Y, self.data.W), full_output=1)
        #coeffs, smin = simplex(SSD, self.eq.theta, self.eq, self.data.X, self.data.Y, self.data.W)
        
        simp = Simplex(SSD, self.eq.theta, *(self.eq.to_fit, (self.data.X, self.data.Y, self.data.W)))
        result = simp.run()
        self.eq.theta = result.rd['x']
       
    def calculate_errors(self):

        self.Smin, self.eq.theta = SSD(self.eq.theta, (self.eq.to_fit,
            (self.data.X, self.data.Y, self.data.W)))
        
        #print '\n SSD \n', Smin
        #hes = errors.hessian(coeffs, eq, set)
        #print '\n Observed information matrix = \n', hes
        self.covar = covariance_matrix(self.eq.theta, self.eq.to_fit, (self.data.X, self.data.Y, self.data.W))
        #print '\n Covariance matrix = \n', covar
        self.correl = correlation_matrix(self.covar)
        self.aproxSD = approximateSD(self.eq.theta, self.eq.to_fit, (self.data.X, self.data.Y, self.data.W))
        self.CVs = 100.0 * self.aproxSD / self.eq.theta
        self.kfit = len(np.nonzero(np.invert(self.eq.fixed))[0])
        self.ndf = self.data.size() - self.kfit
        self.var = self.Smin / self.ndf
        self.Sres, self.Lmax = sqrt(self.var), -SSDlik(self.eq.theta, self.eq.to_fit, (self.data.X, self.data.Y, self.data.W))

        tval = tvalue(self.ndf)
        self.m = tval * tval / 2.0
        self.clim = sqrt(2. * self.m)
        self.Lcrit = self.Lmax - self.m
        self.Llimits = lik_intervals(self.eq.theta, self.aproxSD, self.eq, (self.data.X, self.data.Y, self.data.W))
        
        #self.output.write(self.string_estimates())
        #self.output.write(self.string_liklimits())
        
    def string_estimates(self):
        j = 0
        str = 'Number of point fitted = {0:d}'.format(self.data.size())
        
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
    




