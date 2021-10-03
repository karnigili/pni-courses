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

def_params = {'height':58,'distance':25}
def_timewindow = [0,1]
plt.style.use('ggplot')


def extract_paired(data_path):
    '''
    input: path string
    output: raw data (see dict)
    '''
    data_tmp = pyabf.ABF(data_path)
    
    #nerve
    data_tmp.setSweep(sweepNumber=0, channel=0) 
    nerve_v = data_tmp.sweepY
    nerve_label = data_tmp.sweepLabelY

    times = data_tmp.sweepX

    # muscle
    data_tmp.setSweep(sweepNumber=0, channel=1) 
    muscle_v = data_tmp.sweepY
    muscle_label = data_tmp.sweepLabelY
    
    
    return {'nerve_data':nerve_v,
            'muscle_data':muscle_v,
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

def _get_waveform(data_array, peaks, edge=2000, window=30):
    '''
    get waveform of symm window (int), removing edge peaks (int)
    '''
    waveforms=[]
    non_edge_peaks=[]
    for peak in peaks:
        if(peak>=edge and peak<=
len(data_array)-edge): #edge cases ignored
            waveforms.append(data_array[peak-window:peak+window])
            non_edge_peaks.append(peak)
        
    return non_edge_peaks, waveforms
        
        
def peak_finding(data_array, params=def_params):
    '''
    input: data array (np array 1d)
    '''
    
    peaks,_ = find_peaks(data_array, 
                         height=params['height'], 
                         distance=params['distance'])
    
    
    non_edge_peaks, _ = _get_waveform(data_array, peaks)
    peak_heights = data_array[non_edge_peaks]

    return non_edge_peaks, peak_heights
    
    
def svd_2d(data_array, peaks):
    
    _, waveforms = _get_waveform(data_array, peaks)
    
    M = np.array(waveforms)
    covariance_matrix=np.cov(M.T)
    u, s, vh = np.linalg.svd(covariance_matrix)

    Mo = M-np.mean(M, axis=0) 
    projection_on_PC1=np.dot(Mo,u[:,0])
    projection_on_PC2=np.dot(Mo,u[:,1])
    projection_on_PC3=np.dot(Mo,u[:,2])

    pc_X = np.array([projection_on_PC1, projection_on_PC2, projection_on_PC3]).T
    
    return pc_X

def kmeans_spike_sorting(X, k):
    
    km = KMeans(n_clusters=k).fit(X)
    
    return km.labels_

def plot_mean_waveforms(clusters, waveforms, peaks, peak_heights):
    
    fig, axes = plt.subplots(1,2,figsize=(10, 3),)
    plt.set_cmap('Set3')
    colors = plt.cm.Dark2(np.linspace(0, 1, max(clusters)+2))

    axes = axes.ravel()
    for cluster,color in zip(set(clusters),colors):

        
        wave = [waveforms[i] for i in range(len(peaks)) 
                               if clusters[i]==cluster]
        
        amp = [peak_heights[i] for i in range(len(peaks)) 
                               if clusters[i]==cluster]
        
        std_wave = np.std(wave,axis=0)
        mean_wave = np.mean(wave,axis=0)
        mean_amp = np.mean(amp,axis=0)
        std_amp = np.std(amp,axis=0)
        
        half_len_wave = int(len(mean_wave)/2)
        timeline = range(-1*half_len_wave, half_len_wave)
        

        axes[0].plot(timeline, mean_wave, label='Cluster {}'.format(cluster+1), color=color)
        axes[0].fill_between(timeline, mean_wave-std_wave, mean_wave+std_wave, alpha=.4, color=color)
    
        axes[1].bar(cluster,mean_amp, yerr= std_amp, label='Cluster {}'.format(cluster+1), color=color)
        
    axes[0].set_xlim(-10,10)   
    axes[0].set_label('mean waveform')
    axes[0].set_xlabel('Times (s)')
    axes[0].set_ylabel('Amplitude (mV)')
    
    axes[1].set_label('mean amplitude')
    axes[1].set_xlabel('Cluster')
    axes[1].set_ylabel('Amplitude (mV)')
    
    plt.tight_layout()
    plt.legend(loc="best")
    plt.show()
    

    
    
def plot_two_panels(array_nerve, array_muscle, times, x_lim=False):
    '''
    simple plot of nerve and muscle side-by-side
    '''
    fig, axes = plt.subplots(2, 1, figsize=(10, 5), sharey=False)

    axes = axes.ravel()

    axes[0].plot(times, array_nerve)
    axes[0].set_xlabel('times (s)')
    axes[0].set_label('uV')
    axes[1].scatter(times, array_muscle)
    axes[1].set_xlabel('times (s)')
    axes[1].set_ylabel('mV')
    
    if x_lim:
        axes[0].set_xlim(x_lim[0],x_lim[1])
        axes[1].set_xlim(x_lim[0],x_lim[1])

    
    plt.tight_layout()
    plt.show()
    
def plot_highlighted_peaks(nerve_data, times, peaks, x_lim=False):
    
    fig, axes = plt.subplots(2, 1, figsize=(10, 5), sharey=False)

    axes[0].plot(times, nerve_data, alpha=0.5)
    axes[0].scatter(times[peaks], nerve_data[peaks], color='red')

    axes[1].hist(nerve_data[peaks], bins=100)

    
    if x_lim:
        axes[0].set_xlim(x_lim[0],x_lim[1])

    
    plt.tight_layout()
    plt.show()


    
def plot_highlighted_AP_JP_old(nerve_data, muscle_data, times, peaks, k_labels, x_lim=False,y_lim=False, splash_onset=False):
    
    colors = plt.cm.Dark2(np.linspace(0, 1, max(k_labels)+2))
    plt.set_cmap('Dark2')
    
    fig, axes = plt.subplots(2, 1, figsize=(10, 5), sharex = True, sharey=False)

    axes[0].plot(times, nerve_data, alpha=0.8)
    axes[0].scatter(times[peaks], nerve_data[peaks], color=[colors[c] for c in k_labels])

    axes[1].plot(times, muscle_data)
    for x,c in zip(times[peaks], k_labels):
        axes[1].axvline(x=x, linestyle='dashed', color=colors[c], alpha=0.6)
    

    if x_lim:
        axes[0].set_xlim(x_lim[0],x_lim[1])
        axes[1].set_xlim(x_lim[0],x_lim[1])
    if y_lim:
        axes[1].set_ylim(y_lim[0],y_lim[1])
        
    if splash_onset:
        axes[0].axvline(splash_onset, linestyle='solid', color='grey', label='splash onset')
        axes[1].axvline(splash_onset, linestyle='solid', color='grey',label='splash onset')
        plt.legend(loc="best")
        
    axes[0].set_title('Action Potentials')
    axes[0].set_xlabel('times (s)')
    axes[0].set_ylabel('uV')
    
    axes[1].set_title('Junction Potentials')
    axes[1].set_xlabel('times (s)')
    axes[1].set_ylabel('mV')
    
    
    plt.tight_layout()
    plt.show()

    
def plot_highlighted_AP_JP(nerve_data, muscle_data, times, peaks, k_labels, x_lim=False,y_lim=False, splash_onset=False):
    
    colors = plt.cm.Dark2(np.linspace(0, 1, max(med_amp_clusters)+2))


    plt.set_cmap('Dark2')

    fig, axes0 = plt.subplots(figsize=(10, 3), sharex = True, sharey=False)

    axes0.plot(times, med_filter_n, alpha=0.8)
    axes0.scatter(times[peaks], nerve_data[peaks], color=[colors[c] for c in k_labels])
    axes1 = axes0.twinx()
    axes1.plot(times, muscle_data, c='grey')
    for x,c in zip(times[peaks], k_labels):
        axes1.axvline(x=x, linestyle='dashed', color=colors[c], alpha=0.2)


    if x_lim:
        axes[0].set_xlim(x_lim[0],x_lim[1])
        axes[1].set_xlim(x_lim[0],x_lim[1])
    if y_lim:
        axes[1].set_ylim(y_lim[0],y_lim[1])

    if splash_onset:
        axes[0].axvline(splash_onset, linestyle='solid', color='grey', label='splash onset')
        axes[1].axvline(splash_onset, linestyle='solid', color='grey',label='splash onset')
        plt.legend(loc="best")

    axes0.set_title('Paired recording: Action Potentials/ Junction Potentials')
    axes0.set_xlabel('times (s)')
    axes0.set_ylabel('AP Amplitude (uV)')

    axes1.set_ylabel('JP Amplitude (mV)')

    plt.tight_layout()
    plt.show()


def run_spike_sorting_analysis(file_path, params={'k':3, 'peak_finding':def_params}):

   
    #extract data
    file_data = extract_paired(file_path)
    #filter
    filter_n,_ = filter_data(file_data['nerve_data'])
    filter_m,_ = filter_data(file_data['muscle_data'])
    ## correct offest 
    filter_m -= 70
    times = file_data['times']
    #spike sorting
    non_edge_peaks, peaks_height = peak_finding(filter_n, params['peak_finding'])
    
    amp_clusters = kmeans_spike_sorting(filter_n[non_edge_peaks].reshape([-1,1]), params['k'])
          

    return non_edge_peaks, peaks_height, amp_clusters, times, filter_n,filter_m  


def plot_summary(non_edge_peaks, peaks_height, amp_clusters, times, filter_n,filter_m, params={'splash_onset':False, 'x_lim':def_timewindow, 'y_lim':False}):
    
    _, waveforms = _get_waveform(filter_n, non_edge_peaks) 
    
    plot_highlighted_AP_JP(filter_n, filter_m, times, 
                                   non_edge_peaks, amp_clusters, x_lim=params['x_lim'],y_lim=params['y_lim'],splash_onset=params['splash_onset'])
        
    plot_mean_waveforms(amp_clusters, waveforms, non_edge_peaks, peaks_height)


    
    