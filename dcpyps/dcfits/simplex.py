from math import sqrt
import numpy as np

def rosen(x, args=None):
    """The Rosenbrock function"""
    return sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0)

def rosen1(x, args=None):
    """The Rosenbrock function"""
    return sum(100.0*(x[1:]-x[:-1]**2.0)**2.0 + (1-x[:-1])**2.0), x

class SimplexResult():
    def __init__(self, rd):
        self.rd = rd
        
    def __repr__(self):
        text = ''
        for o in self.rd:
            text += str(o).rjust(8) + '  :  ' + str(self.rd[o]) + '\n'
        return text
        
class Simplex():
    """Minimize a function using the Nelder-Mead simplex algorithm to find the
    minimum of function."""
    def __init__(self, func, theta, *args, **kwargs):
        self.func = func
        self.theta = theta
        self.args = args
        self.k = np.size(theta)
        self._set_settings(kwargs)
        
        self._can_run = True
        self.nit = 1
        self.nfev = 0
        
    def _func(self, x, *args):
        """
        Note that function to be minimised should return function value and
        parameters. There could be functions which modify parameters inside 
        (eg parameters out of limits, etc).
        """
        self.nfev += 1
        ans = self.func(x, *args)
        if isinstance(ans, tuple):
            return ans
        elif isinstance(ans, float):
            return ans, x
        
    def _set_settings(self, kwargs):
        """
        """
        self._settings = {'step' : 0.02, 'reflect' : 1.0, 'extend' : 2.0, 
            'contract' : 0.5, 'shrink' : 0.5, 'restart' : 10.0, 'perturb' : 0.1,
            'parerror' : 1e-8, 'funerror' : 1e-8, 'printiter' : False, 
            'maxiter' : 10000, 'maxevaluations' : 100000, 'maxrestarts' : 3} 
        intersect = set(self._settings.keys()).intersection(set(kwargs.keys()))
        for o in intersect:
            self._settings[o] = kwargs[o]
        
        self._stepfac = self._settings['step']
        self._reflectfac = self._settings['reflect']
        self._extendfac = self._settings['extend']
        self._contractfac = self._settings['contract']
        self._shrinkfac = self._settings['shrink']
        self._restartfac = self._settings['restart']
        self._perturbfac = self._settings['perturb']
        self._parerror = self._settings['parerror']
        self._funerror = self._settings['funerror']
        self._printiter = self._settings['printiter']
        self._maxiter = self._settings['maxiter']
        self._maxfeval = self._settings['maxevaluations']
        self._maxrestart = self._settings['maxrestarts']
        
    def _get_settings(self):
        return self._settings
    settings = property(_get_settings, _set_settings)
                
    def run(self):
        while self._can_run:
            self._prepare()
            while (self.nfev < self._maxfeval and self.nit < self._maxiter):
                centre = np.sum(self.simp[:-1,:], axis=0) / float(self.k)
                fxr, xr = self._reflect(centre) # Reflect
                if fxr < self.fval[0]:
                    fxe, xe = self._extend(centre, xr) # Extend
                    if fxe < fxr:
                        self.simp[-1], self.fval[-1] = xe, fxe
                    else:
                        self.simp[-1], self.fval[-1] = xr, fxr
                else: # fval[0] <= fxr reflected vertex is not better than the best
                    if fxr < self.fval[-2]:
                        self.simp[-1], self.fval[-1] = xr, fxr
                    else: # fxr >= fval[-2]
                        if fxr < self.fval[-1]:
                            self.simp[-1], self.fval[-1] = xr, fxr
                        fxc, xc = self._contract(centre) # Perform contraction
                        if fxc <= self.fval[-1]:
                            self.simp[-1], self.fval[-1] = xc, fxc
                        else:
                            self._shrink() # Shrink
                self._sort()
                self.nit += 1
                if self._printiter: self._print_iteration()
                
                if self._did_converge(): # Check for convergence.
                    self._message = "Simplex terminated successfully."
                    return self._result()
                if self.nfev >= self._maxfeval:
                    self._message = "Maximum number of evaluations has been exceeded."
                    return self._result()
                if self.nit >= self._maxiter:
                    self._message = "Maximum number of iterations has been exceeded."
                    return self._result()
         
    def _result(self):
        return SimplexResult({'nit' : self.nit, 'nfev' : self.nfev, 'message' : self._message,
                'fval' : self.fval[0], 'x' : self.simp[0]})
            
    def _reflect(self, centre):
        xr = (centre + self._reflectfac * (centre - self.simp[-1]))
        return self._func(xr)
    
    def _extend(self, centre, xr):
        xe = centre + self._extendfac * (xr - centre)
        return self._func(xe)

    def _contract(self, centre):
        xc = centre + self._contractfac * (self.simp[-1] - centre)
        return self._func(xc)

    def _shrink(self): #fval, simp, shrfac, func, args):
        for j in range(1, np.size(self.fval)):
            ar = self.simp[0] + self._shrinkfac * (self.simp[j] - self.simp[0])
            self.fval[j], self.simp[j] = self._func(ar)

    def _prepare(self):
        self.simp = np.zeros((self.k + 1, self.k), dtype=self.theta.dtype)
        self.fval = np.zeros((self.k + 1,))
        self.fval[0], self.theta = self._func(self.theta)
        self.simp[0] = self.theta
        step = np.ones((self.k)) * self._stepfac
        fac = (sqrt(self.k + 1) - 1.) / (self.k * sqrt(2.))
        for i in range(self.k):
            ar = self.theta + step * fac
            ar[i] = self.theta[i] + step[i] * (fac + 1. / sqrt(2))
            self.fval[i+1], self.simp[i+1] = self._func(ar)
        self._sort()
        
    def _sort(self):
        ind = np.argsort(self.fval)
        self.fval = np.take(self.fval,ind,0)
        # sort so simp[0,:] has the lowest function value
        self.simp = np.take(self.simp,ind,0)

    def _did_converge(self):
        if (max(np.ravel(abs(self.simp[1:] - self.simp[0]))) <= self._parerror \
            and max(abs(self.fval[0] - self.fval[1:])) <= self._funerror):
            #do local search
            floc, xloc = self._local_search()
            if floc < self.fval[0]:
                self.theta = xloc
                self._stepfac = self._restartfac * self._parerror
                self._can_run = True
            else:
                self._can_run = False
            return True
        else: return False

    def _local_search(self):
        """
        Do local search with each param +/- crtstp at current best vertex.
        """
        fxl1, xl1 = self._func(self.theta + self._parerror)
        fxl2, xl2 = self._func(self.theta - self._parerror)
        if fxl1 < fxl2: return fxl1, xl1
        else: return fxl2, xl2

    def _print_iteration(self):
        if (self.nit % 10) == 0:
            print ("\niter# {0:d}\tlik= {1:f}\ttheta=\n".format(self.nit, -self.fval[0]))
            for th in self.simp[0]: print("{0:6e}\t".format(th))


if __name__ == "__main__":
    x0 = np.array([1.3, 0.7, 0.8, 1.9, 1.2])
    from scipy.optimize import minimize
    res = minimize(rosen, x0, method='nelder-mead',
                options={'xtol': 1e-8, 'disp': True})
    print(res)
    print("=================================")
    print("Testing dcfits Simplex class")
    kwargs = {'maxiter' : 10000, 'maxevaluations' : 100000,
              'parerror' : 1e-8, 'funerror' : 1e-8}
    x0 = np.array([1.3, 0.7, 0.8, 1.9, 1.2])
    simp = Simplex(rosen1, x0, **kwargs)
    result = simp.run()
    print(result)
