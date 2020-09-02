# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 08:26:23 2020

@author: R. Jordan Lear
"""

import os

import numpy as np
import matplotlib.pyplot as plt

from scipy.io import wavfile
from scipy.signal import hilbert, butter, sosfilt, correlate

def get_wav(file):
    if not os.path.isfile(file):
        raise ValueError("%s is not  a file." % file)
    return wavfile.read(file)

if __name__=="__main__":
    file = r".\data\the_supremes_where_did_our_love_go.wav"
    
    samplerate, data = get_wav(file)
    duration = data.shape[0]/samplerate
 
    start = 0
    dur = 10
    
    stop = np.floor(samplerate*(start+dur)).astype(int)+1
    start = np.floor(samplerate*start).astype(int)
    t_sample = np.arange(start, stop)/samplerate
    sample = data[start:stop,0]
    
    hilb = hilbert(sample)
    env = np.abs(hilb)
    
    sos = butter(10, 100, 'lp', fs=samplerate, output='sos')
    filtered = sosfilt(sos, env-env[0])+env[0]
    
    plt.plot(t_sample, sample)
    plt.plot(t_sample, env, 'r')
    plt.plot(t_sample, filtered, 'k--')
    plt.show()
    
    corr = correlate(filtered, filtered, mode='same')
    plt.plot(t_sample,corr)
    plt.show()