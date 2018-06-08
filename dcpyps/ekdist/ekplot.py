import numpy as np
import matplotlib.pyplot as plt

from dcpyps.ekdist import eklib
from dcpyps import dcplots

def plot_stability_intervals(rec, open=True, shut=True, popen=True):
    opma, shma, poma = eklib.moving_average_open_shut_Popen(rec.opint[:-1], rec.shint, window=50)
    x = np.linspace(0, np.prod(opma.shape), num=np.prod(opma.shape), endpoint=True)
    fig = plt.figure(figsize=(6,3))
    ax = fig.add_subplot(111)
    if open:
        ax.semilogy(x, opma, 'r', label='Open periods')
    if shut:
        ax.semilogy(x, poma, 'b', label='Popen')
    if popen:
        ax.semilogy(x, shma, 'g', label='Shut periods')
    #plt.legend()
    print('RED- Open periods\nGREEN- Shut intervals\nBLUE- Popen')
    
def plot_stability_amplitudes(rec, fc, n=2):

    long_opamp = eklib.amplitudes_openings_longer_Tr(rec, fc, n)
    fig = plt.figure(figsize=(6,3))
    ax = fig.add_subplot(111)
    ax.plot(long_opamp, '.b')
    ax.set_ylim([0, 1.2 * max(long_opamp)])
    print('Average open amplitude = ', np.average(long_opamp))

def plot_fitted_amplitude_histogram(rec, fc, n=2, nbins=20):
    long_opamp = eklib.amplitudes_openings_longer_Tr(rec, fc, n)
    fig = plt.figure(figsize=(6,3))
    ax = fig.add_subplot(111)
    ax.hist(long_opamp, nbins, normed=True); #, 50
    ax.set_xlim([0, 1.2 * max(long_opamp)])
    print('Range of amplitudes: {0:.3f} - {1:.3f}'.
          format(min(long_opamp), max(long_opamp)))
          
def plot_xlog_interval_histogram(intervals, tres, shut=False):
    #oxout, oyout, odx = dcplots.prepare_xlog_hist(intervals, rec.tres)
    fig = plt.figure(figsize=(8,3))
    ax = fig.add_subplot(121)
    dcplots.xlog_hist_data(ax, intervals, tres, shut)
    print('Mean and SD of {0:d} open periods = {1:.6g} +/- {2:.6g} ms'.format
          (len(intervals), np.average(intervals)*1000, np.std(intervals)*1000))
    print('\tRange from {0:.6g} to {1:.6g} ms'.format(min(intervals)*1000,
                                                  max(intervals)*1000))
                                                  
def plot_xlog_interval_histogram_fit(intervals, tres, func, pars, shut=False):
    fig = plt.figure(figsize=(4,3))
    ax = fig.add_subplot(111)
    dcplots.xlog_hist_EXP_fit(ax, tres, intervals, pdf=func, pars=pars, shut=shut) 
