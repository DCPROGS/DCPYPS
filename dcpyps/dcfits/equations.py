from math import sqrt, log, pi, ceil
from scipy import stats
import numpy as np

class Equation(object):
    def __init__(self):
        """        """
        self.eqname = None
        self.ncomp = 1
        self.pars = None
        self.fixed = []
        self.names = []
        self.data = None
        self.guess = None
        self._theta = None
        self.normalised = False
        
    def equation(self, x, coeff):
        ''' '''
        pass
    
    def calculate_random(self, x, sd):
        """ """
        if isinstance(x, float):
            return np.random.normal(self.equation(x, self.pars), sd, 1)[0]
        elif isinstance(x, list) or isinstance(x, np.ndarray):
            resp = []
            for each in x:
                resp.append(np.random.normal(self.equation(each, self.pars), sd, 1)[0])
            return np.array(resp)
        
    def to_fit(self, x, theta):
        self._set_theta(theta)
        return self.equation(x, self.pars)
    
    def normalise(self, data):
        pass
    
    def _set_theta(self, theta):
        for each in np.nonzero(self.fixed)[0]:   
            theta = np.insert(theta, each, self.pars[each])
        self.pars = theta
    def _get_theta(self):
        theta = self.pars[np.nonzero(np.invert(self.fixed))[0]]
        if isinstance(theta, float):
            theta = np.array([theta])
        return theta
    theta = property(_get_theta, _set_theta)
    
    def __repr__(self):
        txt = "equation " + self.eqname + "\n"
        for name, par in zip(self.names, self.pars):
            txt += "{} = {}\n".format(name, par)
        return txt
    
class SSR(Equation):
    def __init__(self, equation, dataset, pars=None):
        '''
        Calculate the sum of the squares of residuals (deviations predicted
        from actual empirical values of data). It is a measure of the 
        discrepancy between the data and an estimation model. 

        Parameters
        ----------
        func : callable func(x, args)
            The objective function which calculates predicted values of Y.
        X, Y, W : ndarrays, shape (n, )
            x, y and weights of data points.
        '''
        self.func = equation
        self.X, self.Y, self.W = dataset.X, dataset.Y, dataset.W
        self.nfit = self.X.size
        self.fixed = equation.fixed
        self.kfit = len(self.fixed) - np.sum(self.fixed)
        self.pars = pars
        

        
    def residuals(self, pars):
        '''
        Calculate the weighted residuals.

        Parameters
        ----------
        pars : array_like, shape (k, )
            Function parameters.

        Returns
        -------
        residuals : ndarray, shape (n, )
            Weighted sum of squared deviations.
        '''
        
        return self.W * (self.Y - self.func.to_fit(self.X, pars))
    
    def equation(self, pars):
        """
        Calculate weighted sum of squared residuals.
        
        Parameters
        ----------
        pars : array_like, shape (k, )
            Function parameters.

        Returns
        -------
        SSR : float
            Weighted sum of squared residuals.
        """
        return np.sum(self.residuals(pars)**2)
    
    def to_fit(self, theta):
        self._set_theta(theta)
        return self.equation(self.pars)

    def loglik(self, pars):
        """
        Calculate log likelihood coresponding to the sum of the squared 
        residuals assuming that errors follow Gausian distribution.

        Parameters
        ----------
        pars : array_like, shape (k, )
            Function parameters.

        Returns
        -------
        loglik : float
            Log likelihood of SSR function.
        """
        #n = self.X.size
        S = self.equation(pars)
        Sres = sqrt(S / float(self.nfit)) # - len(func.fixed)))
        return self.nfit * log(sqrt(2 * pi) * Sres) + S / (2 * Sres**2)
    

class GHK(Equation):
    """Goldman-Hodgkin-Katz equation for bi-ionic condition"""
    def __init__(self, eqname, pars=None):
        """
        pars = [r]
        """
        self.eqname = eqname
        self.pars = pars
        self.fixed = [False, True, True, True]
        self.names = ['r', 'totOut', 'In1', 'In2']
        
    def equation(self, x, pars):
        '''
        The GHK equation.
        '''
        return 25 * np.log( (x + pars[0] * (pars[1] - x)) / (pars[2] + pars[0] * pars[3]) )
   
    def propose_guesses(self, data=None):
        '''
        '''
        self.guess = self.pars.copy()
        
    def calculate_plot(self, X, coeff):
        plX = np.linspace(np.floor(np.amin(X)), np.ceil(np.amax(X)), 100)
        plY = self.equation(plX, coeff)
        return plX, plY

    
class Linear(Equation):
    def __init__(self, eqname, pars=None):
        """
        pars = [a, b]
        """
        self.eqname = eqname
        self.ncomp = 1
        self.pars = pars
        self.fixed = [False, False]
        self.names = ['a', 'b']
        
    def equation(self, x, coeff):
        '''
        The linear equation.
        '''
        return coeff[0] + coeff[1] * x 
   
    def propose_guesses(self, data):
        '''
        Calculate the initial guesses for fitting with Linear equation.
        '''
        #if self.Component == 1:
        slope, intercept, r, p, stderr = stats.linregress(data.X, data.Y)
        self.guess = np.array([intercept, slope])
        self.pars = self.guess.copy()
        
    def calculate_plot(self, X, coeff):
        plotX = np.linspace(np.floor(np.amin(X) - 1),
            np.ceil(np.amax(X) + 1), 100)
        plotY = self.equation(plotX, coeff)
        return plotX, plotY


class MultiExponentialPDF(Equation):
    def __init__(self, dataset, taus=None, areas=None, eqname='ExpPDF'):
        self.X = dataset # numpy array
        self.nfit = self.X.size
        self.eqname = eqname
        self.fixed = None
        self.set_pars(taus, areas)
        
    def set_pars(self, taus, areas, fixed=None):
        self.taus = taus
        self.areas = areas
        if taus is not None:
            self.ncomp = len(taus)
            if areas is None:
                pass
            elif len(areas) == len(taus) - 1:
                self.areas = np.append(self.areas, 1 - np.sum(areas))
            if fixed is None and self.fixed is None:
                self.fixed = [False, False] * self.ncomp
                self.fixed[-1] = True
        #self.names = ['tau', 'area']
        self.kfit = len(self.fixed) - np.sum(self.fixed)
        self.pars = np.append(self.taus, self.areas)
        
    def to_fit(self, theta):
        self._set_theta(theta)
        return self.equation(self.taus, self.areas)
    
    def equation(self, taus, areas):
        y = np.array([])
        for t in np.nditer(self.X):
            y = np.append(y, np.sum((areas / taus) * np.exp(-t / taus)))
        return y
    
    def to_plot(self, t):
        y = np.array([])
        for ti in np.nditer(t):
            y = np.append(y, np.sum((self.areas / self.taus) * np.exp(-ti / self.taus)))
        return y
    
    def loglik(self, theta):
        self._set_theta(theta)
        self.taus[self.taus < 1.0e-30] = 1e-8
        self.areas[self.areas > 1.0] = 0.99999
        self.areas[self.areas < 0.0] = 1e-6
        if np.sum(self.areas[:-1]) >= 1: 
            self.areas[:-1] = 0.99 * self.areas[:-1] / np.sum(self.areas[:-1])
        self.areas[-1] = 1 - np.sum(self.areas[:-1])
        d = np.sum( self.areas * (np.exp(-min(self.X) / self.taus) - np.exp(-max(self.X)/ self.taus)))
        if d < 1.e-37:
            print (' ERROR in EXPLIK: d = ', d)
        s = 0.0
        for t in np.nditer(self.X):
            s -= log(np.sum((self.areas / self.taus) * np.exp(-t / self.taus)))
        return s + len(self.X) * log(d) #, theta

    def _expand_pars(self):
        self.taus, self.areas = np.split(np.asarray(self.pars), 2)
        self.areas[-1] = 1 - np.sum(self.areas[:-1])

    def _set_theta(self, theta):
        for each in np.nonzero(self.fixed)[0]:   
            theta = np.insert(theta, each, self.pars[each])
        self.pars = theta
        self._expand_pars()
    def _get_theta(self):
        theta = self.pars[np.nonzero(np.invert(self.fixed))[0]]
        if isinstance(theta, float):
            theta = np.array([theta])
        self._expand_pars()
        return theta
    theta = property(_get_theta, _set_theta)
    
    def __number_per_comp(self):
        f1 = np.sum(self.areas * np.exp(-min(self.X) / self.taus))  #Prob(obs>ylow)
        f2 = np.sum(self.areas * np.exp(-max(self.X) / self.taus))  #Prob(obs>yhigh)
        antrue = len(self.X) / (f1 - f2)
        en = antrue * self.areas
        enout = [antrue * (1. - f1), antrue * f2]
        return en, enout
    
    def __repr__(self):
        repstring = ''
        numb, numout = self.__number_per_comp()
        for ta, ar, nu in zip(self.taus, self.areas, numb):
            repstring += ('Tau = {0:.6f}; lambda (1/s)= {1:.6f}\n'.
                format(ta, 1.0 / ta))
            repstring += ('Area= {0:.6f}; number = {1:.3f};'.format(ar, nu) +
                'amplitude (1/s) = {0:.3f}\n'.format(ar / ta))
        mean = np.sum(self.areas * self.taus)
        repstring += '\nOverall mean = {0:.6f}\n'.format(mean)
        repstring += 'Predicted true number of events = {0:d}\n'.format(int(np.sum(numb)))
        repstring += 'Number of fitted = {0:d}\n'.format(len(self.X))
        repstring += ('Number below Ylow = {0:.3f}; '.format(numout[0]) + 
            'number above Yhigh = {0:.3f}\n'.format(numout[1]))
        return repstring


class Exponential(Equation):
    def __init__(self, eqname, pars=None):
        """
        pars = [a, b]
        """
        self.eqname = eqname
        self.ncomp = 1
        self.pars = pars
        self.fixed = [False, False]
        self.names = ['area', 'tau']
        
    def equation(self, x, coeff):
        '''
        The exponential equation.
        '''
        return coeff[0] / coeff[1] + math.exp(-x /coeff[1]) 

    def calculate_plot(self, X, coeff):
        plotX = np.linspace(np.floor(np.amin(X) - 1),
            np.ceil(np.amax(X) + 1), 100)
        plotY = self.equation(plotX, coeff)
        return plotX, plotY


class Hill(Equation):
    def __init__(self, eqname, pars=None):
        """
        pars = [Ymin, Ymax, EC50, nH]
        """
        self.eqname = eqname
        self.ncomp = 1
        if pars:
            self.pars = pars
        else:
            self.pars = [0.0, 100.0, 10.0, 1.0]
        if eqname == 'Hill':
            self.fixed = [True, False, False, False]
            self.names = ['Ymin', 'Ymax', 'EC50', 'nH  ']
        if eqname == 'Langmuir':
            self.fixed = [True, False, False, True]
            self.names = ['Ymin', 'Ymax', 'EC50']
        self.normpars = None
        self.normalised = False
        
    def equation(self, conc, coeff):
        '''
        The hill equation.
        '''
        return (coeff[0] + ((coeff[1] - coeff[0]) * (conc / coeff[2]) ** coeff[3]) / 
            (1 + (conc / coeff[2]) ** coeff[3]))
            
    def normalise(self, data):
        '''
        Nomalise Y to the fitted maximum.
        '''
        # Nomalise the coefficients by fixing the Y(0) and Ymax
        self.normpars = self.pars.copy()
        self.normpars[0], self.normpars[1] = 0, 1
        #data.normY = (data.Y - self.pars[0]) / self.pars[1]
        data.normY = data.Y / self.pars[1]
        data.normS = data.S / self.pars[1]
        self.normalised = True
    
    def propose_guesses(self, data):
        '''
        Calculate the initial guesses for fitting with Hill equation.
        '''
        #if self.Component == 1:
        self.guess = np.empty(4)
        if data.increase: # Response increases with concentration
            # Determine Y(0)
            if self.fixed[0]:
                self.guess[0] = 0
            else:
                self.guess[0] = np.mean(data.Y[data.X == data.X[0]])
            if self.fixed[1]:
                self.guess[1] = 1
            else:
                # Determine Ymax
                self.guess[1] = np.mean(data.Y[data.X == data.X[-1]])# - self.guess[0]
        else: # Response decreases with concentration
            # Determine Y(0)
            if self.fixed[0]:
                self.guess[0] = 0
            else:
                self.guess[0] = np.mean(data.Y[data.X == data.X[-1]])
            # Determine Ymin
            self.guess[1] = np.mean(data.Y[data.X == data.X[0]])# - self.guess[0]
        # Determine Kr
        self.guess[2] = 10 ** ((np.log10(data.X[0]) + np.log10(data.X[-1])) / 2)
        # Determine nH  
        LinRegressX = np.log10(data.X[data.Y < np.amax(data.Y)]) - np.log10(self.guess[2])
        ratio = data.Y[data.Y < np.amax(data.Y)] / np.amax(data.Y)
        LinRegressY = np.log10(ratio / (1 - ratio))
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            LinRegressX, LinRegressY)
        if math.isnan(slope):
            self.guess[3] = 1.0 if data.increase else -1.0
        else:
            self.guess[3] = slope
        if self.eqname == 'Langmuir':
            self.guess[3] = slope / math.fabs(slope)
#        elif self.Component == 2:
#            print 'Two Components fitting is not completed.'
#            sys.exit(0)
        self.pars = self.guess.copy()
        
    def calculate_plot(self, X, coeff):
        plotX = 10 ** np.linspace(np.floor(np.amin(np.log10(X)) - 1),
            np.ceil(np.amax(np.log10(X)) + 1), 100)
        plotY = self.equation(plotX, coeff)
        return plotX, plotY

