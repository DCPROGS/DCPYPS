
__author__ = "Remis"
__date__ = "$22-Feb-2016 09:13:22$"

import math
import numpy as np
import matplotlib.pyplot as plt

def xlog_hist_data(X, tres, shut=True, unit='s'):
    """
    Plot dwell time histogram in log x and square root y.
    """
    xout, yout, dx = prepare_xlog_hist(X, tres)
    plt.plot(xout, np.sqrt(yout))
    plt.xlabel('Apparent {0} periods ({1})'.
        format('shut' if shut else 'open', unit))
    plt.ylabel('Square root of frequency density')
    plt.xscale('log')
    
def xlog_hist_data_fit(tres, X=None, pdf=None, ipdf=None, iscale=None, 
                       shut=True, tcrit=None, unit='s'):
    """
    Plot dwell time histogram and predicted pdfs (ideal and corrected for the 
    missed events) in log x and square root y.
    """

    # Plot apparent shut period histogram
    if X:
        tout, yout, dt = prepare_xlog_hist(X, tres)
        plt.plot(tout, np.sqrt(yout))
        scale = len(X) * math.log10(dt) * 2.30259
        
    if pdf and ipdf:
        if shut:
            t = np.logspace(math.log10(tres), math.log10(tcrit), 512)
            plt.plot(t, np.sqrt(t * ipdf(t) * scale * iscale), '--r', 
                label='Ideal distribution')
            plt.plot(t, np.sqrt(t * pdf(t) * scale), '-b', 
                label='Corrected distribution')
            t = np.logspace(math.log10(tcrit), math.log10(max(X) *2), 512)
            plt.plot(t, np.sqrt(t * pdf(t) * scale), '--b')
            plt.plot(t, np.sqrt(t * ipdf(t) * scale * iscale), '--r')
            plt.axvline(x=tcrit, color='g')
        else:
            t = np.logspace(math.log10(tres), math.log10(2 * max(X)), 512)
            plt.plot(t, np.sqrt(t * ipdf(t) * scale * iscale), '--r',
                label='Ideal distribution')
            plt.plot(t, np.sqrt(t * pdf(t) * scale), '-b', 
                label='Corrected distribution')
                
    plt.xlabel('Apparent {0} periods ({1})'.
        format('shut' if shut else 'open', unit))
    plt.ylabel('Square root of frequency density')
    plt.legend(loc=(1 if shut else 3))
    plt.xscale('log')

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