# Characterization of the bandwith of PPG signal and its pulse peak

This script implements the methods presented in the article: 
**"Optimal filter characterization for photoplethysmography-based pulse rate and pulse power spectrum estimation"**
by Raymundo Cassani, Abhishek Tiwari, and Tiago H. Falk

In recent years, PPG-based heart rate measurement has gained significant attention due to its popularity in wearable devices. Studies comparing the dynamics of ECG- and PPG-based heart rate measures have found small differences between these two modalities; differences related to the physiological processes behind each technique. We analyzed the spectral coherence and the SNR between isolated PPG pulses and the raw PPG signal to: 

1. Determine the optimal filter to enhance pulse detection from raw PPG for improved heart rate estimation, and
2. Characterize the spectral content of the PPG pulse. 

The methods evaluated in the [Capnobase IEEE TBME Respiratory Rate Benchmark dataset](http://www.capnobase.org/database/pulse-oximeter-ieee-tbme-benchmark/), which contained over 27,000 pulses from PPG recordints from 13 adults and 29 children.

We hope that the results presented herein serve as a baseline for pulse detection algorithms and assist with the development of more sophisticated PPG processing algorithms.

## Results
![alt](https://user-images.githubusercontent.com/8238803/79368838-094cdc80-7f1e-11ea-9001-3cc9fa076e7f.png)  
Figure 1. Average pulse signals for both groups: adults (blue) and children (red). Signals have been scaled for sake of comparison, and the shadowed area presents one standard deviation.

![alt](https://user-images.githubusercontent.com/8238803/79368846-0baf3680-7f1e-11ea-8ade-ce9964b144f7.png)  
Figure 2. Power spectra, coherence and SNR for: (A) Adult, and (B) Children. The dotted lines indicate the BW of the optima filter, and the pulse signal.
