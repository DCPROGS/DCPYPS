#! /usr/bin/python

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

__author__ = "Remis"
__date__ = "$14-Sep-2016 10:15:40$"

import os
import math
from scipy.optimize import minimize
import numpy as np
from numpy import linalg as nplin
from pylab import *

from dcpyps import dataset
from dcpyps import dcplots
#from dcpyps import samples

def hessian(theta, LLfunc, args):
    """
    """
    hess = np.zeros((theta.size, theta.size))
    deltas = optimal_deltas(theta, LLfunc, args)
    # Diagonal elements of Hessian
    coe11 = np.array([theta.copy(), ] * theta.size) + np.diag(deltas)
    coe33 = np.array([theta.copy(), ] * theta.size) - np.diag(deltas)
    for i in range(theta.size):
        hess[i, i] = ((LLfunc(coe11[i], args) - 
            2.0 * LLfunc(theta, args) +
            LLfunc(coe33[i], args)) / (deltas[i]  ** 2))
    # Non diagonal elements of Hessian
    for i in range(theta.size):
        for j in range(theta.size):
            coe1, coe2, coe3, coe4 = theta.copy(), theta.copy(), theta.copy(), theta.copy()
            if i != j:                
                coe1[i] += deltas[i]
                coe1[j] += deltas[j]
                coe2[i] += deltas[i]
                coe2[j] -= deltas[j]
                coe3[i] -= deltas[i]
                coe3[j] += deltas[j]
                coe4[i] -= deltas[i]
                coe4[j] -= deltas[j]
                hess[i, j] = ((
                    LLfunc(coe1, args) -
                    LLfunc(coe2, args) -
                    LLfunc(coe3, args) +
                    LLfunc(coe4, args)) /
                    (4 * deltas[i] * deltas[j]))
    return hess

def optimal_deltas(theta, LLfunc, args):
    """ """
          
    Lcrit = LLfunc(theta, args) + math.fabs(LLfunc(theta, args)*0.005)
    deltas = 0.001 * theta
    L = LLfunc(theta + deltas, args)
    if L < Lcrit:
        count = 0
        while L < Lcrit and count < 100:
            deltas *= 2
            L = LLfunc(theta + deltas, args)
            count += 1
    elif L > Lcrit:
        count = 0
        while L > Lcrit and count < 100:
            deltas *= 0.5
            L = LLfunc(theta + deltas, args)
            count += 1
    return deltas


def covariance_matrix(theta, func, args, weightmode=1):
    """ """
    cov = nplin.inv(hessian(theta, func, args))
#    if weightmode == 1:
#        errvar = SSD(theta, (func, args))[0] / (args[0].size - theta.size)
#    else:
#        errvar = 1.0
    return cov #* errvar

def correlation_matrix(covar):
    correl = np.zeros((len(covar),len(covar)))
    for i1 in range(len(covar)):
        for j1 in range(len(covar)):
            correl[i1,j1] = (covar[i1,j1] / 
                np.sqrt(np.multiply(covar[i1,i1],covar[j1,j1])))
    return correl

def theta_unsqueeze_simult(theta, npatch):
    '''
    length is the number of patches.
    returns the list of thetas.
    '''
    num = (len(theta) + npatch) // (1 + npatch)
    tau = theta[:num]
    thetas = []
    for i in range(npatch):
        new_theta = tau[:]
        new_theta = np.append(new_theta,theta[num+i*(num-1): 2*num+i*(num-1)-1])
        thetas.append(new_theta)

    return thetas

def theta_unsqueeze(theta):
    theta = np.asarray(theta)
    tau, area = np.split(theta, [int(math.ceil(len(theta) / 2))])
    area = np.append(area, 1 - np.sum(area))

    return tau, area

def myexp(theta, X):
    tau, area = theta_unsqueeze(theta)
    X = np.asarray(X)
    y = np.array([])
    for t in np.nditer(X):
        y = np.append(y, np.sum((area / tau) * np.exp(-t / tau)))
    return y

def LL_simult(theta, X):
    num_set = len(X)
    thetas = theta_unsqueeze_simult(theta,num_set)
    s = 0.0
    for i in range(num_set):
        s += LL(thetas[i], X[i])
    return s
    
def LL(theta, X):
    tau, area = theta_unsqueeze(theta)
    tau[tau < 1.0e-30] = 1e-8
    area[area > 1.0] = 0.99999
    area[area < 0.0] = 1e-6
    if np.sum(area[:-1]) >= 1: 
        area[:-1] = 0.99 * area[:-1] / np.sum(area[:-1])
    area[-1] = 1 - np.sum(area[:-1])
    
    d = np.sum( area * (np.exp(-min(X) / tau) - np.exp(-max(X)/ tau)))
    if d < 1.e-37:
        print (' ERROR in EXPLIK: d = ', d)
    X = np.asarray(X)
    s = 0.0
    for t in np.nditer(X):
        s -= math.log(np.sum((area / tau) * np.exp(-t / tau)))
    #theta = np.append(tau, area[:-1])
    return s + len(X) * math.log(d) #, theta

def number_per_comp(theta, X):
    tau, area = theta_unsqueeze(theta)
    f1 = np.sum(area * np.exp(-min(X) / tau))  #Prob(obs>ylow)
    f2 = np.sum(area * np.exp(-max(X) / tau))  #Prob(obs>yhigh)
    antrue = len(X) / (f1 - f2)
    en = antrue * area
    enout = [antrue * (1. - f1), antrue * f2]
    return en, enout
          
def load_patches(dir, conc):
    """
    Load 
    """
    
    data = np.loadtxt(dir + '.csv', 
                  delimiter=',', dtype={
                  'names': ('filename', 'concentration', 'res', 'tcrit'),
                  'formats': ('S8', 'S3', 'f4', 'f4')})
    patches = []
    for patch in data:
        if float(patch['concentration']) == conc:
            filename = "{}/{}.scn".format(dir, patch['filename'].decode('utf8'))
            rec = dataset.SCRecord([filename], float(patch['concentration'])*1e-3, 
                                       patch['res']*1e-6, patch['tcrit']*1e-6)
            patches.append(rec)
    return patches

def extract_intervals(recs, lim, interval_type='open'):
    intervals = []
    for rec in recs:
        if interval_type == 'open':
            temp = np.array(rec.opint)
        elif interval_type == 'shut':
            temp = np.array(rec.shint)
        temp = temp[temp < lim]
        intervals.append(temp)
    return intervals
          
def extract_taus_areas_from_theta(theta, ncomp):
    taus = theta[:ncomp]
    areas = np.reshape(theta[ncomp:], 
                       (len(theta[ncomp:]) // (ncomp - 1), (ncomp - 1)))
    return taus, areas

def are_taus_areas_positive(taus, areas, verbose=False):
    
    if verbose: print('tau=', taus)
    if verbose: print('areas=', areas)
    if verbose: print('all(tau > 0)\t all(area < 1)\t ' +
                      'all(area > 0)\t all(sum_areas < 1)')
    if verbose: print(str(all(taus > 0)) + '\t' +
                      str(all(areas < 1)) + '\t' +
                      str(all(areas > 0)) + '\t' + 
                      str(all(np.sum(areas, axis=1) < 1)))

    return (all(taus > 0) and all(areas < 1) and 
            all(areas > 0) and all(np.sum(areas, axis=1) < 1))

def get_new_theta(taus, areas, ncomp):
            
    if are_taus_areas_positive(taus, areas):
        mean_area = np.mean(areas,axis = 0)
        theta = np.append(taus, np.tile(mean_area, len(areas)))
    else:
        if any(taus < 0):
            taus = np.abs(taus)
        satisified = []
        for area in areas:
            if all(area > 0) and sum(area) < 1:
                satisified.append(area)
        if satisified:
            area = np.vstack(satisified)
            mean_area = np.mean(area,axis = 0)
        else:
            mean_area = np.random.uniform(low=0.0, high=1.0, 
                                          size=len(taus)-1)/len(taus)
        theta = np.append(taus, np.tile(mean_area, ncomp))
    return theta

def generate_random_theta_multiExpPDF(ncomp):
    suggest_tau = np.array([1 / math.pow(10, 3-i) for i in range(ncomp)])
    return np.hstack(
        (np.random.random_sample(ncomp) * suggest_tau, #np.array([1e-4,1e-3,1e-2,1e-1]),
        np.random.random_sample(ncomp-1)/(ncomp-1)))

def fit_histo_expPDF_simult(data, ig, ncomp, verbose=False):
    theta = ig
    likelihood = 0.0
    repeat = True
    while repeat:
        res = minimize(LL_simult, theta, args=data, method='Nelder-Mead')
        taus, areas = extract_taus_areas_from_theta(res.x, ncomp)

        if are_taus_areas_positive(taus, areas, verbose) and res.fun >= likelihood:   
            repeat = False
            if verbose: print('Fitting finished.')
            if verbose: print('Final likelihood = ', likelihood)
        else:
            if are_taus_areas_positive(taus, areas, verbose) and res.fun < likelihood:
                likelihood = res.fun
                if verbose: print('New likelihood = ', likelihood)
            best = res.x
            theta = get_new_theta(taus, areas, ncomp)    
    return best
    
def fit_histo_expPDF_batch(data, ig, ncomp, verbose=False):
    
    thetas = ig
    liks = np.empty(len(data))
    liks.fill(-float('inf'))
    repeat = True
    while repeat:
        if verbose: print('liks=', liks) 
        repeat = False
        theta = np.average(thetas, axis = 0)
        for i in range(len(data)):    
            res = minimize(LL, theta, args=np.array(data[i]), method='Nelder-Mead')
            while any(res.x <0):
                if verbose: print('negative area encountered')
                #if verbose: print(res.x)
                theta = generate_random_theta_multiExpPDF(ncomp)
                #theta = np.abs(theta) * np.random.random_sample(len(theta))
                res = minimize(LL, theta, args=np.array(data[i]), method='Nelder-Mead')
            if res.fun > liks[i]:
                repeat = True                
                thetas[i] = res.x
                liks[i] = res.fun
    return thetas, liks
            
def display_fits(thetas, ints, tres, is_shut=False):
    
    fig, ax  = subplots(1, len(ints), figsize=(4 * len(ints), 4))
    for i in range(len(ints)):
        dcplots.xlog_hist_EXP_fit(ax[i], tres[i], ints[i], 
            pdf=myexp, pars=thetas[i], shut=is_shut) 
#        print_exps(thetas[i], recs[i].shint)
    show()

def print_simult_fit_results(thetas, ints):
    pars = theta_unsqueeze_simult(thetas, len(ints))
    for i in range(len(ints)):
        print_exps(pars[i], ints[i])
    print_averages(pars)

def print_averages(thetas):
    taus, areas = [], []
    for theta in thetas:
        tau, area = theta_unsqueeze(theta)
        taus.append(tau)
        areas.append(area)
    tau_av = np.average(np.transpose(np.array(taus)), axis=1)
    tau_se = np.std(np.transpose(np.array(taus)), axis=1) / math.sqrt(len(thetas))
    area_av = np.average(np.transpose(np.array(areas)), axis=1)
    area_se = np.std(np.transpose(np.array(areas)), axis=1) / math.sqrt(len(thetas))
    print('\nAverages of {0} patches:'.format(len(thetas)))
    for t, dt, a, da in zip(tau_av, tau_se, area_av, area_se):
        print('tau (ms)= {0:.6f} +/- {1:.6f} ; area (%)= {2:.3f} +/- {3:.3f}'.
            format(t*1000, dt*1000, a*100, da*100))
    return tau_av, tau_se, area_av, area_se

def print_exps(theta, X):
    tau, area = theta_unsqueeze(theta)
    numb, numout = number_per_comp(theta, X)
    for ta, ar, nu in zip(tau, area, numb):
        print('\nTau = {0:.6f} ms; lambda (1/s)= {1:.6f}'.format(ta*1000, 1.0 / ta))
        print('Area= {0:.6f}; number = {1:.3f}; amplitude (1/s) = {2:.3f}'.format(ar, nu, ar / ta))
    mean = np.sum(area * tau)
    print('\nOverall mean = {0:.6f}'.format(mean * 1000))
    print('Predicted true number of events = ', np.sum(numb))
    print('Number of fitted = ', len(X))
    print('Number below Ylow = {0:.3f}; number above Yhigh = {1:.3f}'.
          format(numout[0], numout[1]))


def load_shut_intervals(datadir, conc, lim):
    recs = load_patches(datadir, conc)
    shints = extract_intervals(recs, lim, interval_type='shut')
    tres = [recs[i].tres for i in range(len(recs))]
    return shints, tres

def fit_expPDF_batch_independent(data, ncomp):
    ig = np.tile(generate_random_theta_multiExpPDF(ncomp), (len(data), 1))
    return fit_histo_expPDF_batch(data, ig, ncomp, verbose=True)

def load_data(csvfile, concentration, is_open=True, limit=None):
    data = np.loadtxt(csvfile, delimiter=',', dtype={
                  'names': ('filename', 'concentration', 'res', 'tcrit'),
                  'formats': ('S8', 'S3', 'f4', 'f4')})
    recs, intervals = [], []
    for patch in data:
        if float(patch['concentration']) == concentration:
            filename = "{0}/{1}.scn".format(os.path.dirname(csvfile), 
                                            patch['filename'].decode('utf8'))
            print(filename)
            rec = dataset.SCRecord([filename], float(patch['concentration'])*1e-3, 
                                       patch['res']*1e-6, patch['tcrit']*1e-6)
            recs.append(rec)
            if is_open:
                tempint = np.array(rec.opint)
            else:
                tempint = np.array(rec.shint)
            if limit is not None:
                tempint = tempint[tempint<0.1]
            intervals.append(tempint)
    tres = [recs[i].tres for i in range(len(recs))]
    return recs, intervals, tres
    
def test_fit_10_shuts_independent(conc, datadir):
    lim = 1 #if conc > 0.3 else 0.1
    shints, tres = load_shut_intervals(datadir, conc, lim)
    ncomp = 3
    thetas, liks = fit_expPDF_batch_independent(shints, ncomp)
    print_results(thetas, shints)
    print_averages(thetas)
    display_fits(thetas, shints, tres, is_shut=True)
    
def test_fit_10_open_simult(conc, csvfile, tau_ig):
    recs, data, tres = load_data(csvfile, conc)
    ig = np.append(tau_ig[:], [1 / len(tau_ig)] * (len(tau_ig) - 1) * len(recs))
    restheta = fit_histo_expPDF_simult(data, ig, len(tau_ig), verbose=True)
    print_simult_fit_results(restheta, data)
    theta_list = theta_unsqueeze_simult(restheta, len(data))
    display_fits(theta_list, data, tres)
    cov = covariance_matrix(restheta, LL_simult, data)
    apprSD = np.sqrt(cov.diagonal())
    for ta, td in zip(restheta, apprSD):
        print('par = {0:.9f}; approximate SD = {1:.9f}'.format(ta, td))
    
if __name__ == "__main__":
    
    datadir = '../samples/etc/EKDIST_patches'
    csvfile = '../samples/etc/EKDIST_patches/EKDIST_patches.csv'
    conc = 10
    #Initial guess of Tau
    tau_ig1 = [1.38E-05,9.85E-05,3.35E-04,2.48E-03]
    tau_ig1a = [1.38E-05,10.9E-05,2.85E-04,2.7E-03]
    tau_ig2 = [14E-05,3.9E-04,2.1E-03]
    #test_fit_10_open_simult(conc, csvfile, tau_ig2)
    test_fit_10_shuts_independent(conc, datadir)
    
#    tau_ig1 = [2.31469875e-05 ,  7.79896578e-04  , 6.73495738e-03  , 1.18752898e-01]
#    theta = np.append(tau_ig1[:], [1 / ncomp] * (ncomp - 1) * len(recs))
#    theta = fit_histo_expPDF_simult(shints, theta, ncomp, verbose=True)
#    thetas = theta_unsqueeze_simult(theta, len(recs))
    