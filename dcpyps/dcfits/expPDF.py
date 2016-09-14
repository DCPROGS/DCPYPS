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
from dcpyps.dcfits import equations

def mega_unsqueeze(theta, length):
    '''
    length is the number of patches.
    returns the list of thetas.
    '''
    num = (len(theta) + length) // (1 + length)
    tau = theta[:num]
    theta_list = []
    for i in range(length):
        new_theta = tau[:]
        new_theta = np.append(new_theta,theta[num+i*(num-1): 2*num+i*(num-1)-1])
        theta_list.append(new_theta)

    return theta_list


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

def sum_LL(theta, X):
    num_set = len(X)
    theta_list = mega_unsqueeze(theta,num_set)
    s = 0.0
    for i in range(num_set):
        s += LL(theta_list[i], X[i])
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
          
def extract_shut_intervals(data):
    patch_list = []
    data_list = []

    for patch in data:
        if float(patch['concentration']) == concentration:
            filename = "../samples/etc/EKDIST_scn/{}.scn".format(patch['filename'].decode('utf8'))
            print(filename)
            rec = dataset.SCRecord([filename], float(patch['concentration'])*1e-3, 
                                       patch['res']*1e-6, patch['tcrit']*1e-6)
            patch_list.append(rec)
            shint = np.array(rec.shint)
            if concentration >= 0.3:
                shint = shint[shint<0.1]
            else:
                shint = shint[shint<1]
            data_list.append(shint)
    return patch_list, data_list

def fit_histo_expPDF_simult(data_list, tau_ig, verbose=False):
    
    ncomp = len(tau_ig)
    patchnum = len(data_list)

    theta = tau_ig[:]
    theta = np.append(theta,[1/ncomp]*(ncomp-1)*patchnum)
    if verbose: print('theta= ', theta)
        #theta1 = tau_ig[:]
        #theta1 = np.append(theta1,[1/compnum]*(compnum)*patchnum)
        #if verbose: print('theta= ', theta1)
    
    likelihood = 0.0
    repeat = True
    while repeat:
        res = minimize(sum_LL, theta, args=data_list, method='Nelder-Mead')
        if verbose: print('Fit success: '+str(res.success))

        tau = res.x[:ncomp]
        if verbose: print('tau=', tau)
        area = res.x[ncomp:]
        areas = np.reshape(area, (len(area) // (ncomp - 1), (ncomp - 1)))
        if verbose: print('areas=', areas)
        sum_areas = np.sum(areas, axis=1)

        if verbose: print('all(tau > 0)\t all(area < 1)\t all(area > 0)\t all(sum_areas < 1)')
        if verbose: print(str(all(tau > 0))+'\t'+str(all(area < 1))+
                          '\t'+str(all(area > 0))+'\t'+str(all(sum_areas < 1)))

        if all(tau > 0) and all(area < 1) and all(area > 0) and all(sum_areas < 1):
            new = sum_LL(res.x,data_list)
            print('New likelihood = ', new)
            if new < likelihood:
                #best = res.x
                likelihood = new
                mean_area = np.mean(areas,axis = 0)
                theta = np.append(tau,np.tile(mean_area, patchnum))
            else:
                repeat = False
                print('Fitting finished.')
        else:
            if any(tau < 0):
                tau = np.abs(tau)
            satisified = []
            for area in areas:
                if all(area > 0) and sum(area) < 1:
                    satisified.append(area)
            if satisified:
                area = np.vstack(satisified)
                mean_area = np.mean(area,axis = 0)
            else:
                mean_area = np.random.uniform(low=0.0, high=1.0, size=len(tau)-1)/len(tau)
            theta = np.append(tau,np.tile(mean_area, patchnum))
            
    return tau, areas
    
def fit_histo_expPDF(data_list, verbose=False):
    
    
    npatch = len(data_list)
    theta_list = np.empty([npatch, 7])
    theta_list.fill(-float('inf'))

    likehood_list = np.empty(npatch)
    likehood_list.fill(-float('inf'))
    repeat = True
    while repeat:
        repeat = False
        for i in range(npatch):
            shut_ints = np.array(data_list[i])
            shut_ints = shut_ints[shut_ints < 1]

            if (theta_list < 0).any():
                while any(theta_list[i] < 0):
                    if verbose: print('negative area encountered')
                    theta = np.hstack((np.random.random_sample(4) * np.array([1e-4,1e-3,1e-2,1e-1]),
                                      np.random.random_sample(3)/3))
                    if verbose: print('Patch number: ', i+1)
                    if verbose: print(theta)
                    res = minimize(LL, theta, args=shut_ints, method='Nelder-Mead')
                    if verbose: print(res.x)    
                    theta_list[i] = res.x
                    likehood_list[i] = LL(res.x, shut_ints)
                    repeat = True
            else:
                theta = np.average(theta_list, axis = 0)
                if verbose: print('Patch number: ', i+1)
                print(theta)
                res = minimize(LL, theta, args=shut_ints, method='Nelder-Mead')
                print(res.x)
                likehood = LL(res.x, shut_ints)
                if likehood > likehood_list[i]:
                    theta_list[i] = res.x
                    likehood_list[i] = LL(res.x, shut_ints)
                    repeat = True
    return theta_list
            
def display_fits(recs, tau, areas):
    
    k = len(recs)
    fig, ax  = subplots(1, k, figsize=(4 * k, 4))
    for i in range(k):
        dcplots.xlog_hist_EXP_fit(ax[i], recs[i].tres, recs[i].shint, 
            pdf=myexp, pars=np.append(tau, areas[i]), shut=True) 
        print_exps(np.append(tau, areas[i]), recs[i].shint)
    show()

def display_fits2(recs, thetas):
    
    k = len(recs)
    fig, ax  = subplots(1, k, figsize=(4 * k, 4))
    for i in range(k):
        dcplots.xlog_hist_EXP_fit(ax[i], recs[i].tres, recs[i].shint, 
            pdf=myexp, pars=thetas[i], shut=True) 
#        print_exps(np.append(tau, areas[i]), recs[i].shint)
    show()

    
if __name__ == "__main__":
    data = np.loadtxt('../samples/etc/EKDIST_patches.csv', 
                  delimiter=',', dtype={
                  'names': ('filename', 'concentration', 'res', 'tcrit'),
                  'formats': ('S8', 'S3', 'f4', 'f4')})
    concentration = 0.1
    print('Concentration: ', concentration) 

    recs, shut_list = extract_shut_intervals(data)
            
    tau_ig1 = [2.31469875e-05 ,  7.79896578e-04  , 6.73495738e-03  , 1.18752898e-01]
    #tau1, areas1 = fit_histo_expPDF_simult(shut_list, tau_ig1, verbose=True)
    
    
    thetas = fit_histo_expPDF(shut_list, verbose=True)
    
    #display_fits(recs, tau1, areas1)
    display_fits2(recs, thetas)
    