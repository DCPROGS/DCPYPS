import math
import sys

import numpy as np
from scipy.optimize import bisect

def expPDF(t, tau, area):
    """
    Calculate exponential probabolity density function.

    Parameters
    ----------
    t : float or ndarray.
        Time.
    tau : ndarray, shape(k, 1)
        Time constants.
    area : ndarray, shape(k, 1)
        Component relative area.

    Returns
    -------
    f : float or ndarray.
    """

    if type(tau) == type(np.array(())):
        f = np.zeros(t.shape)
        for i in range(tau.shape[0]):
            f += (area[i] / tau[i]) * np.exp(-t / tau[i])
    else:
        f = (area / tau) * np.exp(-t / tau)
    return f

def expPDF_mean_sd(tau, area):
    """
    Calculate mean and standard deviation for exponential PDF.

    Parameters
    ----------
    tau : ndarray, shape(k, 1)
        Time constants.
    area : ndarray, shape(k, 1)
        Component relative area.

    Returns
    -------
    m : float
        Mean.
    sd : float
        Standard deviation.
    """

    m = np.sum(area * tau)
    var = np.sum(area * tau * tau)
    sd = math.sqrt(2 * var - m * m)

    return m, sd

def expPDF_printout(eigs, ampl):
    """
    """

    str = ('term\tw\trate (1/sec)\ttau (ms)\tarea (%)\n')
    for i in range(eigs.shape[0]):
        str += ('{0:d}'.format(i+1) +
            '\t{0:.5g}'.format(ampl[i]) +
            '\t{0:.5g}'.format(eigs[i]) +
            '\t{0:.5g}'.format(1000 / eigs[i]) +
            '\t{0:.5g}\n'.format(100 * ampl[i] / eigs[i]))

    mean, sd = expPDF_mean_sd(1 / eigs, ampl / eigs)
    str += ('Mean (ms) =\t {0:.5g}'.format(mean * 1000) +
        '\tSD =\t {0:.5g}'.format(sd * 1000) +
        '\tSD/mean =\t {0:.5g}\n'.format(sd / mean))
        
    return str

def expPDF_misclassified(tcrit, tau, area, comp):
    """
    Calculate number and fraction of misclassified events after division into
    bursts by critical time, tcrit.
    """

    tfast = tau[:comp]
    tslow = tau[comp:]
    afast = area[:comp]
    aslow = area[comp:]

    # Number of misclassified.
    enf = np.sum(afast * np.exp(-tcrit / tfast))
    ens = np.sum(aslow * (1 - np.exp(-tcrit / tslow)))

    # Fraction misclassified.
    pf = enf / np.sum(afast)
    ps = ens / np.sum(aslow)

    return enf, ens, pf, ps

def expPDF_misclassified_printout(tcrit, enf, ens, pf, ps):
    """
    """

    return ('tcrit = {0:.5g} ms\n'.format(tcrit * 1000) +
        '% misclassified: short = {0:.5g};'.format(pf * 100) +
        ' long = {0:.5g}\n'.format(ps * 100) +
        '# misclassified (out of 100): short = {0:.5g};'.format(enf * 100) +
        ' long = {0:.5g}\n'.format(ens * 100) +
        'Total # misclassified (out of 100) = {0:.5g}\n\n'
        .format((enf + ens) * 100))

def theta_unsqueeze(theta):
    theta = np.asarray(theta)
    tau, area = np.split(theta, [int(math.ceil(len(theta) / 2))])
    area = np.append(area, 1 - np.sum(area))
    return tau, area

def calculate_all_tcrits(theta):
    tau, area = theta_unsqueeze(theta)
    tcrits = np.empty((3, len(tau)-1))
    misclassified = np.empty((3, len(tau)-1, 4))
    for i in range(len(tau)-1):
        try:
            tcrit = bisect(expPDF_tcrit_DC, tau[i], tau[i+1], args=(tau, area, i+1))
            enf, ens, pf, ps = expPDF_misclassified(tcrit, tau, area, i+1)
        except:
            print('Bisection with DC criterion failed.\n')
            tcrit = None
            enf, ens, pf, ps = None, None, None, None
        tcrits[0, i] = tcrit
        misclassified[0, i] = np.array([enf, ens, pf, ps])

        try:
            tcrit = bisect(expPDF_tcrit_CN, tau[i], tau[i+1], args=(tau, area, i+1))
            enf, ens, pf, ps = expPDF_misclassified(tcrit, tau, area, i+1)
        except:
            print('Bisection with Clapham & Neher criterion failed.\n')
            tcrit = None
            enf, ens, pf, ps = None, None, None, None
        tcrits[1, i] = tcrit
        misclassified[1, i] = np.array([enf, ens, pf, ps])

        try:
            tcrit = bisect(expPDF_tcrit_Jackson, tau[i], tau[i+1], args=(tau, area, i+1))
            enf, ens, pf, ps = expPDF_misclassified(tcrit, tau, area, i+1)
        except:
            print('Bisection with Jackson criterion failed.\n')
            tcrit = None
            enf, ens, pf, ps = None, None, None, None
        tcrits[2, i] = tcrit
        misclassified[2, i] = np.array([enf, ens, pf, ps])
    return tcrits, misclassified

def printout_all_tcrits(tcrits, misclassified):
    for i in range(len(tcrits)-1):
        print('\nCritical time between components {0:d} and {1:d}\n'.
                format(i+1, i+2) + '\nEqual % misclassified (DC criterion)')
        tcrit = tcrits[0, i]
        if tcrit is not None:
            miscl = misclassified[0, i] 
            print(expPDF_misclassified_printout(tcrit, miscl[0], miscl[1], miscl[2], miscl[3]))

        print('Equal # misclassified (Clapham & Neher criterion)')
        tcrit = tcrits[1, i] 
        if tcrit is not None:
            miscl = misclassified[1, i] 
            print(expPDF_misclassified_printout(tcrit, miscl[0], miscl[1], miscl[2], miscl[3]))

        print('Minimum total # misclassified (Jackson et al criterion)')
        tcrit = tcrits[2, i] 
        if tcrit is not None:
            miscl = misclassified[2, i] 
            print(expPDF_misclassified_printout(tcrit, miscl[0], miscl[1], miscl[2], miscl[3]))

    print('\n\nSUMMARY of tcrit values:\nComponents\t\tDC\t\tC&N\t\tJackson\n')
    for i in range(len(tcrits)-1):
        print('{0:d} to {1:d} '.format(i+1, i+2) +
                '\t\t\t{0:.5g}'.format(tcrits[0, i] * 1000) +
                '\t\t{0:.5g}'.format(tcrits[1, i] * 1000) +
                '\t\t{0:.5g}\n'.format(tcrits[2, i] * 1000))


def expPDF_tcrit_DC(tcrit, tau, area, comp):
    """
    """

    enf, ens, pf, ps = expPDF_misclassified(tcrit, tau, area, comp)
    return ps - pf

def expPDF_tcrit_CN(tcrit, tau, area, comp):
    """
    """

    enf, ens, pf, ps = expPDF_misclassified(tcrit, tau, area, comp)
    return ens - enf

def expPDF_tcrit_Jackson(tcrit, tau, area, comp):
    """
    """

    tfast = tau[:comp]
    tslow = tau[comp:]
    afast = area[:comp]
    aslow = area[comp:]

    # Number of misclassified.
    enf = np.sum((afast / tfast) * np.exp(-tcrit / tfast))
    ens = np.sum((aslow / tslow) * np.exp(-tcrit / tslow))

    return enf - ens

def geometricPDF_mean_sd(rho, w):
    """
    Calculate mean and standard deviation for geometric PDF.

    Parameters
    ----------
    rho : ndarray, shape(k, 1)
        Probabilities.
    w : ndarray, shape(k, 1)
        Component amplitudes.

    Returns
    -------
    m : float
        Mean.
    sd : float
        Standard deviation.
    """

    k = rho.shape[0]
    m = np.sum(w / np.power(np.ones((k)) - rho, 2))
    var = np.sum(w * (np.ones((k)) + rho) / np.power(np.ones((k)) - rho, 3))
    sd = math.sqrt(var - m * m)

    return m, sd

def geometricPDF_printout(rho, w):
    """
    """

    norm = 1 / (np.ones((rho.shape[0])) - rho)
    str = ('term\tw\trho\tarea(%)\tNorm mean')
    for i in range(rho.shape[0]):
        str += ('{0:d}'.format(i+1) +
            '\t{0:.5g}'.format(w[i]) +
            '\t{0:.5g}'.format(rho[i]) +
            '\t{0:.5g}'.format(w[i] * norm[i] * 100) +
            '\t{0:.5g}\n'.format(norm[i]))

    mean, sd = geometricPDF_mean_sd(rho, w)
    str += ('Mean number of openings per burst =\t {0:.5g}'.format(mean) +
        '\n\tSD =\t {0:.5g}'.format(sd) +
        '\tSD/mean =\t {0:.5g}\n'.format(sd / mean))
    return str
