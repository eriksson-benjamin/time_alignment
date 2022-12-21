#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 12:05:12 2022

@author: beriksso
"""

import useful_defs as udfs
import numpy as np
import matplotlib.pyplot as plt
from calculate_coincidences import get_g_loc
import scipy as sp
import sys
sys.path.insert(0, '/home/beriksso/TOFu/analysis/benjamin/github/TOFu/functions')
import tofu_functions as dfs

def gauss(x, amp, mu, sigma):
    """Gaussian function."""
    return amp * np.exp(-(x - mu)**2 / (2 * sigma**2))
    
    
def fit_function(par, x_obs, y_obs):
    """Function to minimize in fitting routine."""
    # Calculate Gaussian
    y_mod = gauss(x_obs, *par)

    return np.sum((y_obs - y_mod)**2 / y_mod)


def S1_shift(coincidences, S1_detector):
    """Find shift of S1 detectors by fitting a Gaussian to the peak."""
    # Get approximate location
    g_loc = get_g_loc(S1_detector)
    
    # Make histogram of coincidences
    bins = np.arange(g_loc - 5, g_loc + 5, 0.1)
    bin_centres = udfs.get_bin_centres(bins)
    counts, _ = np.histogram(coincidences, bins=bins)

    # Minimize chi2
    init_guess = (counts.max(), g_loc, 0.4) 
    x_lim = ((bin_centres > g_loc - 0.7) & (bin_centres < g_loc + 0.7))
    popt = sp.optimize.minimize(fit_function, init_guess, 
                                args=(bin_centres[x_lim], counts[x_lim]))
    
    return popt.x, bin_centres, counts


def S2_shift(coincidences, S2_detector):
    """Find shift of S2 detectors by taking the mode."""
    # Get approximate location
    g_loc = get_g_loc(S2_detector)
    
    # Make histogram of coincidences
    bins = np.arange(g_loc - 10, g_loc + 10, 0.4)
    bin_centres = udfs.get_bin_centres(bins)
    counts, _ = np.histogram(coincidences, bins=bins)

    # Find mode
    argm = np.argmax(counts)
    
    return bin_centres[argm], bin_centres, counts


def plot_S1(par, bin_centres, counts, det_name):
    """Plot S1 data with Gaussian fit."""
    # Plot data    
    plt.figure(det_name)
    plt.errorbar(bin_centres, counts, yerr=np.sqrt(counts), color='k', 
                 linestyle='None', marker='.', markersize=2)
    
    # Plot fit
    x_fit = np.arange(bin_centres[0], bin_centres[-1], 
                      np.diff(bin_centres)[0]/5)
    y_fit = gauss(x_fit, *par)
    plt.plot(x_fit, y_fit, 'r-')
    
    # Configure plot
    plt.xlabel('$t_{TOF}$ (ns)')
    plt.ylabel('counts')
    

def plot_S2(mode, bin_centres, counts, det_name):
    """Plot S2 data with mode of distribution."""
    # Plot data
    plt.figure(det_name)
    plt.errorbar(bin_centres, counts, yerr=np.sqrt(counts), color='k', 
                 linestyle='None', marker='.', markersize=2)
    
    # Plot mode
    plt.axvline(mode, color='r', linestyle='--')
    
    # Configure plot
    plt.xlabel('$t_{TOF}$ (ns)')
    plt.ylabel('counts')
    
    g_loc = get_g_loc(det_name)
    plt.xlim(g_loc - 10, g_loc + 10)
    
if __name__ == '__main__':
    # Read coincidence data
    dat = udfs.unpickle('output/gamma_events.pickle')
    S1_dat = dat['S1']
    S2_dat = dat['S2']
    
    shifts = {}
    
    # S1 shifts
    # ---------
    S1_dict = dfs.get_dictionaries('S1')
    S1_dict.pop('S1_05')
    for S1 in S1_dict.keys():
        det_int = int(S1.replace('S1_', '')) - 1 
        
        # Fit Gaussian
        pars, bin_centres, counts = S1_shift(S1_dat[det_int], S1)
        
        # Plot fit
        plot_S1(pars, bin_centres, counts, S1)
        
        # Add to shift file
        shifts[S1] = pars[1]
        
    # S2 shifts
    # ---------
    S2_dict = dfs.get_dictionaries('S2')
    for S2 in S2_dict.keys():
        det_int = int(S2.replace('S2_', '')) - 1 
        
        # Find mode
        mode, bin_centres, counts = S2_shift(S2_dat[det_int], S2)
    
        # Plot gamma peak + mode
        plot_S2(mode, bin_centres, counts, S2)
        
        # Add to shift file
        shifts[S2] = mode
        
    # Write to file
    with open('shifts.txt', 'w') as handle:
        for sx, shift in shifts.items():
            handle.write(f'{sx} {shift}\n')
        
            
            
            
            
        
        