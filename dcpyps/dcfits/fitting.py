import os
import sys
import numpy as np

import dcpyps.dcfits.dataIO as dataIO
from dcpyps.dcfits.simplex import Simplex

#import cvfit
#from cvfit import data
#from cvfit import errors
#from cvfit.errors import residuals
#from cvfit.errors import SSD
#from cvfit.errors import SSDlik

class MultipleFitSession(object):
    def __init__(self, output=sys.stdout):
        self.output = output
        self.allsessions = []
        self.pooled = None
        
    def add(self, fitset):
        self.allsessions.append(fitset)
        
    def pool(self, norm=False, output=sys.stdout):
        dataset = dataIO.XYDataSet()
        for session in self.allsessions:
            if norm:
                dataset.pool(session.data.X, session.data.normY, session.data.normS)
            else:
                dataset.pool(session.data.X, session.data.Y, session.data.S)
        dataset.weightmode = 2
        self.pooled = SingleFitSession(dataset, self.allsessions[0].eq, output)
        
    def string_average_estimates(self):
        str = 'Average of estimates of {0:d} sets:'.format(len(self.allsessions))
       
        for i in range(len(self.allsessions[0].eq.names)):
            pars = []
            for j in range(len(self.allsessions)):
                pars.append(self.allsessions[j].eq.pars[i])
            str += ('\nParameter {0:d}: {1}  \t= {2:.6g} +/- {3:.6g}'.
                format(i+1, self.allsessions[0].eq.names[i], np.mean(pars), 
                np.std(pars)/sqrt(len(pars))))
        return str
    
    def prepare_fplot(self, type):
        plot_list = []
        if type == 'pooled':
            X, Y = self.pooled.eq.calculate_plot(self.pooled.data.X, self.pooled.eq.pars)
            plot_list.append([X,Y])
        else:
            for session in self.allsessions:
                if type == 'fit':
                    pars = session.eq.pars
                elif type == 'guess':
                    pars = session.eq.guess
                elif type == 'norm':
                    pars = session.eq.normpars
                X, Y = session.eq.calculate_plot(session.data.X, pars)
                plot_list.append([X,Y])
        return plot_list
    

class SingleFitSession(object):
    def __init__(self, dataset, equation, output=sys.stdout):
        """
        """
        self.output = output
        self.data = dataset
        self.eq = equation
        self.output.write('\n\tFitting session for ' + self.data.title + ' initialised!')
        
    def propose_guesses(self):
        try:
            self.eq.propose_guesses(self.data)
            print("Guesses proposed")
        except:
            print("{} equation might not be able to propose guesses".
                format(self.eq.eqname))
                
    def fit(self, func):
        simp = Simplex(func, self.eq.theta)
        result = simp.run()
        self.eq.theta = result.rd['x']
        print("Single fit finished")
       
    def print_estimates(self, errs, liklimits):
        j = 0
        txt = 'Number of point fitted = {0:d}'.format(self.data.size())
        
        txt += '\nNumber of parameters estimated = {0:d}'.format(self.kfit)
        txt += '\nDegrees of freedom = {0:d}'.format(self.ndf)
        
        txt += ('\nResidual error SD = {0:.3f}      (variance = {1:.3f})'.
            format(self.Sres, self.var))
        
        for i in range(len(self.eq.names)):
            txt += '\nParameter {0:d}: {1}  \t= {2:.6g}  \t'.format(i+1, self.eq.names[i], self.eq.pars[i])
            if not self.eq.fixed[i]:
                txt += '  Approx SD = {0:.6g}\t'.format(self.aproxSD[j])
                txt += '  CV = {0:.1f}'.format(self.CVs[j])
                j += 1
            else:
                txt += '  (fixed)'

        txt += ('\nMinimum SSD = {0:.3f}; \nMax log-likelihood = {1:.3f}'.
            format(self.Smin, self.Lmax))
        txt += ('\nCorrelation matrix = ' + 
            '[!!!! PRINTOUT OF CORRELATION MATRIX NOT IMPLEMENTED YET. SORRY.\n')
#        self.output.write(correl)
        if np.any(np.absolute(self.correl - np.identity(self.kfit)) > 0.9):
            txt += ("\nWARNING: SOME PARAMETERS ARE STRONGLY CORRELATED (coeff > 0.9); try different guesses")

        if np.any(self.CVs > 33):
            txt += "\nWARNING: SOME PARAMETERS POORLY DEFINED (CV > 33%); try different guesses"
        return txt

        
    def string_liklimits(self):
        j = 0
        txt = '\nLIKELIHOOD INTERVALS\n'
        txt += ('{0:.3g}-unit Likelihood Intervals'.format(self.m) +
            '  (equivalent SD for Gaussian- {0:.3g})'.format(self.clim))
        txt += '\nLmax= {0:.6g};   Lcrit= {1:.6g}'.format(self.Lmax, self.Lcrit)
        for i in range(len(self.eq.names)):
            txt += '\nParameter {0:d}:   {1}\t= {2:.6g}'.format(i+1, self.eq.names[i], self.eq.pars[i])
            if not self.eq.fixed[i]:
                try:
                    txt += '\t  LOWER = {0:.6g}'.format(self.Llimits[j][0])
                except:
                    txt += '\t  LOWER limit not found'
                try:
                    txt += '\t  UPPER = {0:.6g}'.format(self.Llimits[j][1])
                except:
                    txt += '\t  UPPER limit not found'
                j += 1
            else:
                txt += '\t  (fixed)'
        return txt


def load_data(example=False):

    if example:
        filename = (os.path.dirname(os.path.dirname(cvfit.__file__)) +
            "./Example/Example.xlsx")
    else:
        filename = data.ask_for_file()
    try:
        #allsets = data.read_sets_from_csv(filename, 'csv', col=2, header=0, namesin=False, weight=1)
        allsets = data.read_sets_from_Excel(filename, 2, 0, 0) #, namesin=False, weight=1)
    except ValueError:
        print('fitting.py: WARNING: Oops! File did not load properly...')
    return allsets, filename

