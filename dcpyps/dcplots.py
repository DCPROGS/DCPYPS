"""
Collection of functions useful for plotting in DCPYPS.
"""
__author__ = "Remis"
__date__ = "$22-Feb-2016 09:13:22$"

import math
import numpy as np

def xlog_hist_data(ax, X, tres, shut=True, unit='s'):
    """
    Plot dwell time histogram in log x and square root y.
    """
    
    xout, yout, dx = prepare_xlog_hist(X, tres)
    ax.semilogx(xout, np.sqrt(yout))
    ax.set_xlabel('Apparent {0} periods ({1})'.
        format('shut' if shut else 'open', unit))
    ax.set_ylabel('Square root of frequency density')
    
def xlog_hist_HJC_fit(ax, tres, X=None, pdf=None, ipdf=None, iscale=None, 
                       shut=True, tcrit=None, legend=True, unit='s'):
    """
    Plot dwell time histogram and predicted pdfs (ideal and corrected for the 
    missed events) in log x and square root y.
    """

    scale = 1.0
    if X is not None:
        tout, yout, dt = prepare_xlog_hist(X, tres)
        ax.semilogx(tout, np.sqrt(yout))
        scale = len(X) * math.log10(dt) * math.log(10)
        
    if pdf and ipdf:
        if shut:
            t = np.logspace(math.log10(tres), math.log10(tcrit), 512)
            ax.semilogx(t, np.sqrt(t * ipdf(t) * scale * iscale), '--r', 
                label='Ideal distribution')
            ax.semilogx(t, np.sqrt(t * pdf(t) * scale), '-b', 
                label='Corrected distribution')
            t = np.logspace(math.log10(tcrit), math.log10(max(X) *2), 512)
            ax.semilogx(t, np.sqrt(t * pdf(t) * scale), '--b')
            ax.semilogx(t, np.sqrt(t * ipdf(t) * scale * iscale), '--r')
            ax.axvline(x=tcrit, color='g')
        else:
            t = np.logspace(math.log10(tres), math.log10(2 * max(X)), 512)
            ax.semilogx(t, np.sqrt(t * ipdf(t) * scale * iscale), '--r',
                label='Ideal distribution')
            ax.semilogx(t, np.sqrt(t * pdf(t) * scale), '-b', 
                label='Corrected distribution')
                
    ax.set_xlabel('Apparent {0} periods ({1})'.
        format('shut' if shut else 'open', unit))
    ax.set_ylabel('Square root of frequency density')
    if legend: ax.legend(loc=(1 if shut else 3))
    
def xlog_hist_EXP_fit(ax, tres, X=None, pdf=None, pars=None, shut=True, 
                      tcrit=None, unit='s'):
    """
    Plot dwell time histogram and multi-exponential pdf with single 
    components in log x and square root y.
    """
    
    theta = np.asarray(pars)
    tau, area = np.split(theta, [int(math.ceil(len(theta) / 2))])
    area = np.append(area, 1 - np.sum(area))

    scale = 1.0
    if X is not None:
        tout, yout, dt = prepare_xlog_hist(X, tres)
        ax.semilogx(tout, np.sqrt(yout))
        scale = (len(X) * math.log10(dt) * math.log(10) *
            (1 / np.sum(area * np.exp(-tres / tau))))
        
    t = np.logspace(math.log10(tres), math.log10(2 * max(X)), 512)
    ax.plot(t, np.sqrt(scale * t * pdf(pars, t)), '-b')
    for ta, ar in zip(tau, area):
        ax.plot(t, np.sqrt(scale * t * (ar / ta) * np.exp(-t / ta)), '--b')
        
    if tcrit is not None:
        tcrit = np.asarray(tcrit)
        for tc in tcrit:
            ax.axvline(x=tc, color='g')
        
    ax.set_xlabel('Apparent {0} periods ({1})'.
        format('shut' if shut else 'open', unit))
    ax.set_ylabel('Square root of frequency density')
    
def prepare_xlog_hist(X, tres):
    """
    Prepare data points for x-log histogram to plot in matplotlib
    
    eg. 
    xout, yout, dx = dcplots.prepare_xlog_hist(intervals, resolution)
    plt.plot(xout, yout)
    plt.xscale('log')
    
    Parameters
    ----------
    X :  1-D array or sequence of scalar
    tres : float
        Temporal resolution, shortest resolvable time interval. It is
        histogram's starting point.

    Returns
    -------
    xout, yout :  list of scalar
        x and y values to plot histogram.
    dx : float
        Histogram bin width.

    """

    # Defines bin width and number of bins.
    # Number of bins/decade
    n = len(X)
    if (n <= 300): nbdec = 5
    if (n > 300) and (n <= 1000): nbdec = 8
    if (n > 1000) and (n <= 3000): nbdec = 10
    if (n > 3000): nbdec = 12
    dx = math.exp(math.log(10.0) / float(nbdec))
    xstart = tres    # histogramm starts at
    xmax = max(X)
    # round up maximum value, so get Xmax for distribution
    xend = math.exp(math.ceil(math.log(xmax)))
    nbin = int(math.log(xend / xstart) / math.log(dx))
    
#    xaxis = np.arange(xstart, xend, dx)

    # Make bins.
    xaxis = np.zeros(nbin+1)
    xaxis[0] = xstart
    # For log scale.
    for i in range(1, nbin+1):
        xaxis[i] = xstart * (dx**i)

    # Sorts data into bins.
    freq = np.zeros(nbin)
    for i in range(n):
        for j in range(nbin):
            if X[i] >= xaxis[j] and X[i] < xaxis[j+1]:
                freq[j] = freq[j] + 1

    xout = np.zeros((nbin + 1) * 2)
    yout = np.zeros((nbin + 1) * 2)

    xout[0] = xaxis[0]
    yout[0] = 0
    for i in range(0, nbin):
        xout[2*i+1] = xaxis[i]
        xout[2*i+2] = xaxis[i+1]
        yout[2*i+1] = freq[i]
        yout[2*i+2] = freq[i]
    xout[-1] = xaxis[-1]
    yout[-1] = 0

    return xout, yout, dx


def moving_average(x, n):
    """
    Compute an n period moving average.
    """
    x = np.asarray(x)
    weights = np.ones(n)
    weights /= weights.sum()
    a =  np.convolve(x, weights, mode='full')[:len(x)]
    a[:n] = a[n]
    return a