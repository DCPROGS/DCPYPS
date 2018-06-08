import numpy as np

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

def moving_average_open_shut_Popen(opints, shints, window=50):
    # window : moving average interval
    opma = moving_average(opints, window) # Moving average for open periods
    shma = moving_average(shints, window) # Moving average for shut periods
    poma = opma / (opma + shma) # Moving average for Popen
    return opma, shma, poma

def filter_risetime(fc):
    return 0.3321 / fc

def amplitudes_openings_longer_Tr(rec, fc, n=2):
    all_resolved_ops = np.array(rec.rtint)[np.where( np.fabs(np.asarray(rec.rampl)) > 0.0)]
    all_resolved_opamp = np.array(rec.rampl)[np.where( np.fabs(np.asarray(rec.rampl)) > 0.0)]
    #long_ops = all_resolved_ops[np.where( all_resolved_ops > filter_risetime(fc))]
    long_opamp = all_resolved_opamp[np.where( all_resolved_ops > n * filter_risetime(fc))]
    return np.absolute(long_opamp)