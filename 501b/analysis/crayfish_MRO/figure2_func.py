import matplotlib
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import pyabf
import pandas as pd
from scipy.signal import medfilt
from scipy.signal import find_peaks
from sklearn.cluster import KMeans
from scipy.optimize import curve_fit

def_params = {'height':100,'distance':25}
def_timewindow = [0,1]
plt.style.use('ggplot')

get_ifr = lambda peak_times: 1/np.diff(peak_times)

def decay_func(t, r0, tau, r_inf):
    # R_t = R_{0}\cdot  e^{-t/\tau} +R_{\infty}
    return r0 * np.exp(-t/tau) + r_inf

def extract_paired(data_path):
    '''
    input: path string
    output: raw data (see dict)
    '''
    data_tmp = pyabf.ABF(data_path)
    
    #nerve
    data_tmp.setSweep(sweepNumber=0, channel=0) 
    ch1 = data_tmp.sweepY
    #ch1_label = data_tmp.sweepLabelY

    times = data_tmp.sweepX

    # muscle
    data_tmp.setSweep(sweepNumber=0, channel=1) 
    ch2 = data_tmp.sweepY
    #muscle_label = data_tmp.sweepLabelY
    
    
    return {'ch1':ch1,
            'ch2':ch2,
            'times':times,
            'misc':data_tmp,
           }


def filter_data(data_array, kernel_size=1501):
    '''
    
    performs median filtering (low pass filtering) 
    
    input: data array (np array 1d), kernel size (int)
    output: data array (np array 1d, same shape at input), correction (np array 1d, same shape at input)
    '''
    correction = medfilt(data_array, kernel_size=kernel_size)
    filtered_data = data_array - correction + np.median(data_array)
    
    return filtered_data, correction

        
def peak_finding(data_array, active_window, params=def_params):
    '''
    input: data array (np array 1d), 
        active_window [2,1] start and end of the window
        params
    
    '''
    
    peaks_ids,_ = find_peaks(data_array, 
                         height=params['height'], 
                         distance=params['distance'])
    
    
    non_edge_peaks = peaks_ids[np.bitwise_and(peaks_ids>active_window[0], peaks_ids<active_window[1])]
    peak_heights = data_array[non_edge_peaks]

    return non_edge_peaks, peak_heights
    
    

def fit_decay_curve(non_edge_peaks, times, remove_outliers=False, thresh=0):
    
    # get firing freq & timestamps
    ifr_val = get_ifr(times[non_edge_peaks])
    ifr_times = np.convolve(times[non_edge_peaks], np.ones(2), 'valid') / 2

    
    #decay curve data, only from peak forward
    ifr_peak = np.argmax(ifr_val)
    decay_data = ifr_val[ifr_peak:]
    decay_times = ifr_times[ifr_peak:]
    
    #remove edge outliers
    ok_decay_ids = np.where(np.diff(decay_data)<thresh)[0]
    decay_data = decay_data[ok_decay_ids]
    decay_times = decay_times[ok_decay_ids]
    
    # shift decay times to start w 0
    decay_times_shift = decay_times - decay_times[0]

    _, smoothed_decay_data = filter_data(decay_data, kernel_size=7)
    
    smooth_peak = np.argmax(smoothed_decay_data)
    smoothed_decay_data = smoothed_decay_data[smooth_peak:]
    decay_times_shift = decay_times_shift[smooth_peak:]
    
    #remove low outliers
    if remove_outliers:
        #remove edge outliers
        ok_decay_ids = np.where(np.diff(smoothed_decay_data)<thresh)[0]
        smoothed_decay_data = smoothed_decay_data[ok_decay_ids]
        decay_times_shift = decay_times_shift[ok_decay_ids]
        if ok_decay_ids.shape[0]>5:
            not_too_low_ids = np.where(smoothed_decay_data>smoothed_decay_data[-5])
            smoothed_decay_data=smoothed_decay_data[not_too_low_ids]
            decay_times_shift=decay_times_shift[not_too_low_ids]

    
    # fit
    if smoothed_decay_data.shape[0]>5:
        fit_param, pcov = curve_fit(decay_func, decay_times_shift, smoothed_decay_data, (smoothed_decay_data[0],.5,smoothed_decay_data[-1]))
        perr = np.sqrt(np.diag(pcov))

        curve_data = [decay_times_shift, smoothed_decay_data]
        return fit_param, curve_data, perr 
    return False



    

    
    