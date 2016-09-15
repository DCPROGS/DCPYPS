#! /usr/bin/python

# To change this license header, choose License Headers in Project Properties.
# To change this template file, choose Tools | Templates
# and open the template in the editor.

__author__ = "Remis"
__date__ = "$14-Sep-2016 10:15:40$"

import math
from scipy.optimize import minimize
import numpy as np
from pylab import *

from dcpyps import dataset
from dcpyps import dcplots
#from dcpyps import samples
from dcpyps.dcfits.equations import ExponentialPDF

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
          
#def extract_shut_intervals(data):
#    patch_list = []
#    data_list = []
#
#    for patch in data:
#        if float(patch['concentration']) == concentration:
#            filename = "../samples/etc/EKDIST_scn/{}.scn".format(patch['filename'].decode('utf8'))
#            print(filename)
#            rec = dataset.SCRecord([filename], float(patch['concentration'])*1e-3, 
#                                       patch['res']*1e-6, patch['tcrit']*1e-6)
#            patch_list.append(rec)
#            shint = np.array(rec.shint)
#            if concentration >= 0.3:
#                shint = shint[shint<0.1]
#            else:
#                shint = shint[shint<1]
#            data_list.append(shint)
#    return patch_list, data_list

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

def get_new_theta(taus, areas):
            
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
        theta = np.append(taus, np.tile(mean_area, len(data)))
    return theta

def generate_random_theta_multiExpPDF(ncomp):
    return np.hstack(
          (np.random.random_sample(ncomp)* np.array([1e-4,1e-3,1e-2,1e-1]),
           np.random.random_sample(ncomp-1)/(ncomp-1)))


def fit_histo_expPDF_simult(data, theta, ncomp, verbose=False):
    likelihood = 0.0
    repeat = True
    while repeat:
        res = minimize(LL_simult, theta, args=data, method='Nelder-Mead')
        if verbose: print('Fit success: '+str(res.success))
        taus, areas = extract_taus_areas_from_theta(res.x, ncomp)

        if are_taus_areas_positive(taus, areas, verbose) and res.fun >= likelihood:   
            repeat = False
            if verbose: print('Fitting finished.')
            if verbose: print('Final likelihood = ', res.fun)
        else:
            if are_taus_areas_positive(taus, areas, verbose) and res.fun < likelihood:
                likelihood = res.fun
            theta = get_new_theta(taus, areas)    
    return res.x
    
def fit_histo_expPDF(data, theta, ncomp, verbose=False):
    
    thetas = np.empty([len(data), (2*ncomp-1)])
    thetas.fill(-float('inf'))
    likehood_list = np.empty(len(data))
    likehood_list.fill(-float('inf'))
    
    for i in range(len(data)):
        if verbose: print('Patch number: ', i+1)
        repeat = True
        while repeat:
            res = minimize(LL, theta, args=np.array(data[i]), method='Nelder-Mead')
            if verbose: print('Fit success: '+str(res.success))
            repeat = False
            if (res.x < 0).any():
                if verbose: print('negative tau or area encountered')
                theta = generate_random_theta_multiExpPDF(ncomp)
                if verbose: print(theta)
                repeat = True
            else:
                theta = np.average(thetas, axis = 0)
                print(theta)
                res = minimize(LL, theta, args=np.array(data[i]), method='Nelder-Mead')
                print(res.x)
                if res.fun > likehood_list[i]:
                    thetas[i] = res.x
                    likehood_list[i] = res.fun
                    repeat = True
    return thetas
            
def display_fits(recs, thetas):
    
    k = len(recs)
    fig, ax  = subplots(1, k, figsize=(4 * k, 4))
    for i in range(k):
        dcplots.xlog_hist_EXP_fit(ax[i], recs[i].tres, recs[i].shint, 
            pdf=myexp, pars=thetas[i], shut=True) 
        print_exps(thetas[i], recs[i].shint)
    show()

    
if __name__ == "__main__":
    
    dir = '../samples/etc/EKDIST_patches'
    concentration = 0.1
    lim = 0.1
    recs = load_patches(dir, concentration)
    shints = extract_intervals(recs, lim, interval_type='shut')
    print('Concentration: ', concentration) 
    
    ncomp = 4
    tau_ig1 = [2.31469875e-05 ,  7.79896578e-04  , 6.73495738e-03  , 1.18752898e-01]
    theta = np.append(tau_ig1[:], [1 / ncomp] * (ncomp - 1) * len(recs))
   
    theta = fit_histo_expPDF_simult(shints, theta, ncomp, verbose=True)
    thetas = theta_unsqueeze_simult(theta, len(recs))
    
    
#    theta = generate_random_theta_multiExpPDF(ncomp)
#    print ('random theta:', theta)
#    thetas = fit_histo_expPDF(shints, theta, ncomp, verbose=True)
    
    display_fits(recs, thetas)
