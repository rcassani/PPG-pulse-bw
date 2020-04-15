#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to perform the analysis of the PPG signal to derive an optimal fitler
for PPG peak detection.

This script is an implementation of the methods described in:

"Optimal filter characterization for photoplethysmography-based pulse
rate and pulse power spectrum estimation" 
Raymundo Cassani, Abhishek Tiwari, and Tiago H. Falk"  

It uses data from CapnoDatabase
http://www.capnobase.org/database/pulse-oximeter-ieee-tbme-benchmark/

Requires the 
Amplitude Modulation Analysis Module
https://github.com/MuSAELab/amplitude-modulation-analysis-module
"""

import numpy as np
from am_analysis import am_analysis as ama
import scipy.signal
import matplotlib.pyplot as plt
import os
import h5py
import scipy.io as sio

# %% data folder
data_dir  = 'D:\\data\\ppg databases\\capnobase\\TBME2013-PPGRR-Benchmark_R3\\data\\'

# %% parameters
age_limits = [16, 180] # keeps data from subjects in this range (inclusive)

# odd number of have the same number of samples before and after the pulse peak
samples_pls = 301  #  1s   @300Hz
samples_ppg = 521  #  1.7s @300Hz
diff_samples = samples_ppg - samples_pls

# create windows
win_pls = scipy.signal.get_window(('tukey', 0.8), samples_pls, fftbins=False)
win_ppg = scipy.signal.get_window(('tukey', 0.8), samples_ppg, fftbins=False) 
win_pls_rms = np.mean(win_pls**2)
win_ppg_rms = np.mean(win_ppg**2)


# %% processing
# search for .mat files
files_all = os.listdir(data_dir)
filename_bases = [i[0:-4] for i in files_all if i.endswith('.mat')]

pls_peaks_all = []
ppg_signals_all =[]
age_all = []

for filename_base in filename_bases:    
    print(filename_base)

    # load mat file
    with h5py.File(data_dir + filename_base + '.mat') as file:
        # get age
        age = file['meta']['subject']['age'][0,0]
        # get fs
        fs = file['param']['samplingrate']['pleth'][0,0]       
        # get PPG signal
        x_p = np.squeeze(file['signal']['pleth']['y'][:])
        # get pulses location
        ixs = np.squeeze(file['labels']['pleth']['peak']['x'][:]).astype(int)

    # age of one child ('0032_8min.mat') is missing, children average is used 
    if np.isnan(age):
        age = 8.7

    ixs = ixs[2:-2] # ignore two first and two last beats
    
    pls_peaks = np.zeros((samples_ppg, len(ixs)))  
    ppg_signals = np.zeros((samples_ppg, len(ixs)))  

    for ipls, ix in enumerate(ixs):
        # get PPG centered at the pulse 
        ppg_t = x_p[ix-(samples_ppg//2) : ix+(samples_ppg//2) + 1] 
        ppg_t = ppg_t - np.mean(ppg_t)
        ppg_t_win = (ppg_t * win_ppg) / win_ppg_rms
        ppg_signals[:, ipls] = ppg_t_win
        # get pulse from this segment
        pls_t = ppg_t[(samples_ppg//2)+1 - (samples_pls//2) : (samples_ppg//2) + 1 + (samples_pls//2) + 1]
        # Tukey window for Pulse peak
        pls_t_win = (pls_t * win_pls) / win_pls_rms
        # place at center
        pls_peaks[(diff_samples//2) : (diff_samples//2) + len(pls_t_win), ipls] = pls_t_win   # _ _ * _ _              
            
    if  age_limits[0] <= age <= age_limits[1]:                
        # appends pulse peaks from one file to ALL files list
        pls_peaks_all.append(pls_peaks)
        ppg_signals_all.append(ppg_signals)
        age_all.append(age)


# %% analysis

# concatenate all PPG signals and Pulse peaks
pls_peaks_all_cat = np.concatenate(pls_peaks_all, axis=1)   
ppg_signals_all_cat = np.concatenate(ppg_signals_all, axis=1)   

# PSD
pls_peaks_psd = ama.rfft_psd(pls_peaks_all_cat, fs, 1024)
ppg_signals_psd = ama.rfft_psd(ppg_signals_all_cat, fs, 1024)

# average across instances 
pls_peaks_psd['PSD'] = np.mean(pls_peaks_psd['PSD'], axis=1)[:, None]
ppg_signals_psd['PSD'] = np.mean(ppg_signals_psd['PSD'], axis=1)[:, None]

# compute transfer function
n_fft = 2048
pls_peaks_fft = np.fft.fft(pls_peaks_all_cat, n_fft, axis=0)
ppg_signals_fft = np.fft.fft(ppg_signals_all_cat, n_fft, axis=0)
# adjust for even and odd number of elements
if n_fft % 2 != 0:
    # odd case
    n_freqs = int((n_fft + 1) / 2)
else:
    # even case 
    n_freqs = int((n_fft / 2) + 1)
f_axis = np.arange(0,n_freqs) * fs / n_fft

# compute autospectra
sxx = pls_peaks_fft * np.conjugate(pls_peaks_fft)
syy = ppg_signals_fft * np.conjugate(ppg_signals_fft)
syx = ppg_signals_fft * np.conjugate(pls_peaks_fft)

# average accross instances
sxx_avg = sxx.mean(axis= 1)
syy_avg = syy.mean(axis= 1)
syx_avg = syx.mean(axis= 1)

# magnitude square coherence
msc = np.abs(syx_avg) ** 2 / ((sxx_avg * syy_avg) + np.finfo(float).eps) 
# SNR
snr = msc / (1-msc)

# scale SNR, Syy and Sxx
snr_n = snr / snr.max()
syy_n = syy_avg / syy_avg.max()
sxx_n = sxx_avg / syy_avg.max()

# %% plots
# Figure: Sxx, Syy, MSC and SNR as functions of frequency
fig = plt.figure()
ax1 = fig.add_subplot(111)

asxx = ax1.plot(f_axis, sxx_n[0:n_freqs])
asyy = ax1.plot(f_axis, syy_n[0:n_freqs])
amsc = ax1.plot(f_axis, msc[0:n_freqs])
plt.xlim([0,10])
plt.xlabel('frequency (Hz)')
plt.ylabel('normalized')

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
asnr = ax2.plot(f_axis, 10*np.log10(snr[0:n_freqs]), 'k')

ax1.legend(asxx + asyy + amsc + asnr, ['Sxx', 'Syy', 'MSC', 'SNR'])


#%% Mean PPG signal and mean Pulse peak
x_mean = pls_peaks_all_cat.mean(axis=1) 
x_max  = x_mean.max()
x_std  = pls_peaks_all_cat.std(axis=1)
y_mean = ppg_signals_all_cat.mean(axis=1) 
y_max  = y_mean.max()
y_std  = ppg_signals_all_cat.std(axis=1)

plt.figure()
plt.subplot(2,1,1)
ym = plt.plot(y_mean / y_max,  color='black')
ys = plt.plot( (y_mean + y_std) / y_max, color=(0.5,0.5,0.5))
plt.plot( (y_mean - y_std) / y_max, color=(0.5,0.5,0.5))
plt.ylim(-1.2, 1.7)
plt.xlabel('time (s)')
plt.legend(ym+ys, ['mean', '+- 1 std'])
plt.title('PPG signal')

plt.subplot(2,1,2)
plt.plot(x_mean / x_max, color='black' )
plt.plot( (x_mean + x_std) / x_max, color=(0.5,0.5,0.5))
plt.plot( (x_mean - x_std) / x_max, color=(0.5,0.5,0.5))
plt.ylim(-1.2, 1.7)
plt.xlabel('time (s)')
plt.legend(ym+ys, ['mean', '+- 1 std'])
plt.title('Pulse peak')


#%% save signals as mat file
savemat_data = {
                # average pulse and PPG signals
                'x_mean' : x_mean,
                'x_max'  : x_max,
                'x_std'  : x_std, 
                'y_mean' : y_mean,
                'y_max'  : y_max,
                'y_std'  : y_std, 
                # frequency curves
                'f_axis' : f_axis,
                'msc'    : msc[0:n_freqs],
                'sxx_n'  : sxx_n[0:n_freqs],
                'syy_n'  : syy_n[0:n_freqs],
                'snr_db' : 10*np.log10(snr[0:n_freqs])
        }

sio.savemat('matfile.mat', savemat_data)






