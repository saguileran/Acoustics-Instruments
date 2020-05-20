#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import librosa, librosa.display
import numpy as np
import matplotlib.pyplot as plt
from os import listdir
from os.path import isfile, join
from scipy.io import wavfile
from scipy.signal import find_peaks
import pandas as pd

#mypath = '/home/sebas/Documents/Recorders/WavFiles/'
mypath = '/home/sebas/Documents/Recorders/Felipe/MedicionesFelipe/NotasFelipe/'
#mypath = input('Enter location of folder: ')
#file = 'Flute-D.wav'

def PlotSpectrogram(file):
    plt.clf()
    y, sr = librosa.load(mypath + file)
    librosa.feature.melspectrogram(y=y, sr=sr)
    
    D = np.abs(librosa.stft(y))**2
    S = librosa.feature.melspectrogram(S=D, sr=sr)
    
    S = librosa.feature.melspectrogram(y=y, sr=sr, n_mels=128, fmax=sr/2)
    
    plt.figure(figsize=(10, 4))
    S_dB = librosa.power_to_db(S, ref=np.max)
    librosa.display.specshow(S_dB, x_axis='time',y_axis='mel', sr=sr, fmax=sr/2)
    plt.colorbar(format='%+2.0f dB')
    plt.title(file[:-4]+' Spectrogram')
    plt.savefig(mypath+'Images/'+file[:-4]+'-spectro')
    plt.pause(1)
    
def Plot(file):
    plt.clf()
    samplingFrequency, signalData = wavfile.read(mypath+file)
    dB, dim = 20*np.log10(signalData), len(signalData)
    times = np.arange(0, dim/samplingFrequency, 1/samplingFrequency)
    plt.plot(times, dB)
    plt.xlabel('Time (s)')
    plt.yticks(np.arange(0, (np.nanmax(dB)//10+2)*10, step=10))
    plt.ylabel('dB')
    plt.title(file[:-4] + ' dB')
    plt.savefig(mypath+'Images/'+file[:-4])
    plt.pause(1)

def PlotDbvsHz(file):
    plt.clf()
    data = np.loadtxt(mypath+file,skiprows=1)
    Y, X = [dat[1] for dat in data], [dat[0] for dat in data]
    Y, X = np.array(Y), np.array(X)
    
    plt.fill_between(X,np.min(Y),Y)
    
    plt.ylabel('dB')
    plt.xlabel('Frequency (Hz)')
    plt.grid()
    
    plt.xticks(np.arange(0, 10000, step=1000))
    plt.xlim(0, 10000)
    plt.ylim(np.min(Y), np.max(Y)+5)
    
    peaks, _ = find_peaks(Y, prominence=1) 
    #print(peaks, X[peaks])
    plt.plot(X[peaks], Y[peaks], "ob")
    
    peakscoord = np.array([X[peaks], Y[peaks]])
    #np.savetxt(mypath+file[:-4]+"-peaks.csv", peakscoord, delimiter=' ' )
    
    #exporting with pandas
    df = pd.DataFrame(np.transpose(peakscoord), columns=['Hz', 'dB'])
    df.to_csv(path_or_buf=mypath+'Images/'+file+'-peaks.csv', index=False, sep=' ')
    plt.title(file[:-4]+' Decibels Vs Frequency')
    plt.savefig(mypath+'Images/'+file[:-4]+'-FreVsDb')
    plt.pause(1)

onlyfiles = [file for file in listdir(mypath) if isfile(join(mypath, file))]
formats, names = [formato[-3:] for formato in onlyfiles], [name[:-4] for name in onlyfiles]

for file in onlyfiles:
    if 'wav' in file and 'peak' not in file:
        Plot(file)
        PlotSpectrogram(file)
        print(file)
        print(' ')
    elif 'txt' in file and 'peak' not in file: PlotDbvsHz(file)
