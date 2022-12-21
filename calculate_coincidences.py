#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 10:03:09 2022

@author: beriksso
"""

"""
Calculate coincidences between S1-5 and the other S1 detectors, as well as S1-5
and the S2 detectors.
"""

import sys
sys.path.insert(0, '/home/beriksso/TOFu/analysis/benjamin/github/TOFu/functions')
import tofu_functions as dfs
import useful_defs as udfs
import numpy as np
import matplotlib.pyplot as plt
udfs.set_nes_plot_style()
import os


def plot_tof_hist(tof, bins):
    """Plot histogram of times-of-flight."""
    plt.figure()
    bin_centres = udfs.get_bin_centres(bins)
    counts, _ = np.histogram(tof, bins)
    
    plt.plot(bin_centres, counts, 'k.')
    plt.errorbar(bin_centres, counts, yerr=np.sqrt(counts), color='k', 
                 linestyle='None')
    plt.xlabel('$t_{TOF}$ (ns)')
    plt.ylabel('counts')

    return bin_centres, counts


def S2_coincidences(data, select=True):
    """
    Calculate coincidences between S1-5 and S2's. Save gamma peak events.

    Parameters
    ----------
    data : dict,
        Dictionary containing time data for S1 and S2 events.
    select : bool, optional
        If True, select events from gamma peak only.
        If False, return all TOFs. Default is True.

    Returns
    -------
    tof_list : list,
        List of arrays of TOFs between S1-5 and S2 events.
    ref_list : ndarray,
        Numpy array of unique indices of S1-5 events that have caused a 
        coincidence with an S2.
    """
    tof_list = []  # for all S1-5 vs. S1 or S2 times-of-flight
    ref_list = np.array([])  # for saving S1-5 indices for coincidences with an S2
    
    s2_dict = dfs.get_dictionaries('S2')
    for s2 in s2_dict:
        g_loc = get_g_loc(s2)

        # Select gamma events
        if select:
            # Set search window around gamma peak
            tofs, inds = dfs.sTOF4(data['times_S1']['S1_05'], 
                                   data['times_S2'][s2], t_back=g_loc + 10, 
                                   t_forward=-(g_loc - 10), 
                                   return_indices=True)

        else:
            tofs, inds = dfs.sTOF4(data['times_S1']['S1_05'], 
                                   data['times_S2'][s2], t_back=160, 
                                   t_forward=0, return_indices=True)
        
        tof_list.append(tofs)
        
        # Save unique S1-5 indices which have caused a coincidence with an S2
        ref_list = np.unique(np.append(ref_list, inds))
    
    # Cast ref_inds to int type
    ref_list = np.array(ref_list, dtype='int')

    return tof_list, ref_list

def S1_coincidences(data, ref_list):
    """Calculate coincidences between S1-5 and S1's."""
    t_window = 20
    tof_list = []
    
    # Select only S1-5 events contributing to gamma peak
    selected = data['times_S1']['S1_05'][ref_list]
    for s1 in ['S1_01', 'S1_02', 'S1_03', 'S1_04']:
        tofs, inds = dfs.sTOF4(data['times_S1'][s1], selected, t_back=t_window, 
                               t_forward=t_window, return_indices=True)
        tof_list.append(tofs)
    
    return tof_list


def add_to_array(big_array, small_array):
    """Add events to big S1/S2 arrays of events."""
    if len(big_array) == 0:
        return small_array
    else:
        for i in range(len(big_array)):
            big_array[i] = np.append(big_array[i], small_array[i])
            
        return big_array

def get_g_loc(detector):
    """Return the approximate gamma location for the given detector."""
    g_locs = udfs.json_read_dictionary('input/gamma_locations.json')
    
    return g_locs[detector]
    
def main(data_paths, select=True):
    S1_events = []
    S2_events = []
    for path in data_paths:
        # Import data
        data = udfs.unpickle(path)
        
        # Calculate S1-5 vs. S2 gamma peak cgoincidences (indices + tof)
        S2_tof, ref_list = S2_coincidences(data, select)
    
        # Find coincidences between selected S1-5 events and other S1's
        S1_tof = S1_coincidences(data, ref_list)
    
        # Add events to big array
        S1_events = add_to_array(S1_events, S1_tof)
        S2_events = add_to_array(S2_events, S2_tof)
        
    return S1_events, S2_events


if __name__ == '__main__':
    directory = '/common/scratch/beriksso/TOFu/data/time_alignment/dataset_3'
    files = os.listdir(directory)
    data_paths = [f'{directory}/{file}' for file in files if '.pickle' in file]
    S1_tof, S2_tof = main(data_paths, select=True)
    
    for i in range(4):
        _, counts = plot_tof_hist(S1_tof[i], bins=np.arange(-10, 10, 0.1))
    
    for i in range(32):
        _, counts = plot_tof_hist(S2_tof[i], bins=np.arange(-199.8, 200, 0.4))
    
    
    # Write to file
    udfs.pickler('gamma_events.pickle', {'S1':S1_tof, 'S2':S2_tof})
    
    
    
    
    