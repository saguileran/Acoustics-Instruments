# -*- coding: utf-8 -*-
import scipy.io.wavfile, matplotlib
import numpy as np
import matplotlib.pyplot as plt
from numpy import fft as fft
import librosa, librosa.display
from scipy.signal import find_peaks
from os import listdir
from os.path import isfile, join
import pandas as pd

#Necesary functions
def ReduceVectors(x, y, first = 9):
    divisores, i = x.shape[0]/np.array(range(2,1000))%1, first #primer divisor
    divisores = divisores == 0
    while divisores[i] == False: i += 1
    return np.mean(x.reshape(-1, i+2), axis=1), np.mean(y.reshape(-1, i+2), axis=1)
    
def dft_shift(X):
    N = len(X)
    if (N % 2 == 0):
        # even-length: return N+1 values
        return np.arange(-int(N/2), int(N/2) + 1), np.concatenate((X[int(N/2):], X[:int(N/2)+1]))
    else:
        # odd-length: return N values
        return np.arange(-int((N-1)/2), int((N-1)/2) + 1), np.concatenate((X[int((N+1)/2):], X[:int((N+1)/2)]))

def dft_map(X, Fs, shift=True):
    resolution = float(Fs) / len(X)
    if shift:
        n, Y = dft_shift(X)
    else:
        Y, n = X, np.arange(0, len(X))
    f = n * resolution
    return f, Y

def Peaks(coefficients_weight, frequencies, distance=150, height=-20):
     peaks, _ = find_peaks(coefficients_weight[:samples//2], distance=distance, height=height)
     return frequencies[peaks], coefficients_weight[peaks]

def Smooth(x,window_len=11,window='hanning'):
    dim = len(x)
    #possibl windows ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
    s = np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]

    if window == 'flat': w = np.ones(window_len,'d')#moving average
    else: w=eval('np.'+window+'(window_len)')

    y = np.convolve(w/w.sum(),s,mode='valid')
    smooted = y[(window_len//2-1):-(window_len//2)]
    if dim == len(smooted):return smooted
    else: return smooted[:dim]
    
#------------Main Program----------------
if __name__ == "__main__":
    #a temp folder for downloads
    
    #temp_folder = '/home/sebas/Documents/Acoustics-Instruments/Experiment/Measurements/Audacity/Flute/Test/' #'C.wav'
    temp_folder = '/home/sebas/Documents/Recorders/Test/'
    
    onlyfiles = [file for file in listdir(temp_folder ) if isfile(join(temp_folder , file))]
    formats, names = [formato[-3:] for formato in onlyfiles], [name[:-4] for name in onlyfiles]
    
    wav_files = [wave for wave in onlyfiles if "wav" in wave]
    
    for file in wav_files:
    
    #file = wav_files[0]
        song = temp_folder + file 
        
        rate, audData = scipy.io.wavfile.read(song) #read wav file
        samples, channels = audData.shape[0], len(audData.shape)
        #if chanels==2: audData = np.sum(audData.astype(float), axis=1, keepdims= True)/2
        #else: audData = audData.reshape(samples, 1)
        
        time = np.arange(0, float(samples), 1) / rate #time
        fourier = fft.fft(audData)  #FFT to data
        fourier_scaled = fourier / float(samples) #scaling dat to not depend of length
        frequencies, coefficients_weight = dft_map(fourier_scaled , rate, shift=False) #Fourier space, counting frequencies
        coefficients_weight = abs(coefficients_weight) #Convert to real data, weights are complex
        Maxfretoplot = np.where(frequencies//1 == 5000)[0][0] #Maximum frequency to plot
        
        fourier1 = fourier[0:(samples//2)] / float(samples) # scale by the number of points so that the magnitude does not depend on the length, half of points bcause fouri space is symtric
        
        freqArray = np.arange(0, (samples//2), 1.0) * (rate*1.0/samples) #calculate the frequency at each point in Hz
        
        SPL = 20*np.log10(audData) #Sound Pressure Level
        p0, c = 1, 20
        pressure = p0*10**(SPL/20)
        
        #4,2,2
        x_pure, y_pure = freqArray/1000,  20*np.log10(fourier1)
        
        #4,2,4
        X_pro, Y_pro = ReduceVectors(x_pure, y_pure)
        Y_pro = Smooth(Y_pro, 30)
        X_peaks, Y_peaks = Peaks(Y_pro, X_pro, distance=10, height=5)
        
        #4,2,6
        T, P = ReduceVectors(time, pressure/c, 1)
        
        
        fig = plt.figure(figsize=(40, 24))
        fig.suptitle(file[:-4]) # or plt.suptitle('Main title')
        plt.subplots_adjust(hspace = 0.4, wspace = 0.2)
        
        fig.add_subplot(4,2,1)
        plt.title('Waveform')
        plt.plot(time, audData, linewidth=0.1, color='red')
        plt.xticks(np.arange(0, time[-1], step=0.25))
        plt.xlabel('Time (s)'); plt.ylabel('Amplitude')
        
        fig.add_subplot(4,2,2)
        plt.title('Pure Tone')
        plt.plot(x_pure, y_pure, color='#ff7f00', linewidth=0.5)
        plt.xticks(np.arange(0, 10, step=0.25)); #plt.yticks(np.arange(-50, 50, step=10))
        plt.xlabel('Frequency (kHz)'); plt.ylabel('Power (dB)')
        plt.xlim(0, 10/2)
        
        fig.add_subplot(4,2,3)
        plt.title('dB')
        plt.plot(time, Smooth(SPL, 10), color='red')
        #plt.plot(time, SPL, color='red')
        plt.xlabel('time (s)'); plt.ylabel('dB')
        #plt.xlim([2, 3])
        
        fig.add_subplot(4,2,4)
        plt.title('Smooth Tone and Peaks')
        plt.plot(X_pro , Y_pro.real, color='#ff7f00', linewidth=1.0)
        plt.scatter(X_peaks , Y_peaks.real, s=100)
        for x,y in zip(X_peaks, Y_peaks.real):
            #label = "({:.2f}, {:.2f})".format(1000*x,y)
            label = "{:.2f}".format(1000*x)
            plt.annotate(label, # this is the text
                         (x,y), # this is the point to label
                         textcoords="offset points", # how to position the text
                         xytext=(0,15), # distance from text to points (x,y)
                         ha='center') # horizontal alignment can be left, right or center
        #plt.yticks(np.arange(-50, 50, step=10))
        plt.xlabel('Frequency (kHz)'); plt.ylabel('Power (dB)'); plt.grid()
        plt.xticks(np.arange(0, 10, step=0.25))
        plt.xlim(0, 10/2); plt.ylim(np.min(Y_pro.real), np.max(Y_pro.real)+10)
        
        fig.add_subplot(4,2,5)
        plt.title('Amplitude Coefficients of FFT')
        plt.plot(frequencies, coefficients_weight)
        plt.xlabel('Frequencies (Hz)'); plt.ylabel('Weight')
        plt.xticks(np.arange(0, 2500, step=250))
        plt.xlim(0, 2500)
        
        fig.add_subplot(4,2,6)
        plt.title('Pressure')
        plt.loglog(T, P, color='red')
        plt.xlabel('time (s)'); plt.ylabel('Pressure (Pa)')
        plt.xlim(10**-2, 10)
        
        fig.add_subplot(4,2,7)
        plt.title('Spectrogram')
        y, sr = librosa.load(song)
        librosa.feature.melspectrogram(y=y, sr=sr)
        D = np.abs(librosa.stft(y))**2
        S = librosa.feature.melspectrogram(S=D, sr=sr)
        S = librosa.feature.melspectrogram(y=y, sr=sr, n_mels=128, fmax=sr/2)
        S_dB = librosa.power_to_db(S, ref=np.max)
        librosa.display.specshow(S_dB, x_axis='time',y_axis='mel', sr=sr, fmax=sr/2)
        #plt.yticks(np.arange(0, 10000, step=1000))
        plt.colorbar(format='%+2.0f dB')
        
        plt.savefig(temp_folder + "Plots - " + file[:-4])
        plt.show()
        plt.close()
        
        Data = pd.DataFrame({'freq': X_peaks*1000, 'dB': Y_peaks.real})
        Data.to_csv(temp_folder + "Peaks - " + file[:-4])
        
        #----------------_Harmonics analysis----------
        #2,2,2
        Y_norm = Y_peaks.real - np.min(Y_peaks.real)
        
        #2,2,3
        index_main_frequecy = np.where(np.max(Y_peaks.real) == Y_peaks.real)
        X_norm = X_peaks/X_peaks[index_main_frequecy]
        
        #2,2,4
        Y_peaks_diff = [abs(Y_peaks[i+1] - Y_peaks[i]) for i in range(len(Y_peaks)-1)]
        
        
        fig = plt.figure(figsize=(20, 20))
        plt.subplots_adjust(hspace = 0.2, wspace = 0.2)
        
        
        fig.add_subplot(2,2,1)
        plt.plot(X_peaks, 'o')
        plt.xticks(np.arange(0, len(X_peaks), step=1)); plt.yticks(np.arange(0, np.max(X_peaks)+0.1, step=0.25))
        plt.xlabel("Frequency number (f_i)"); plt.ylabel("Frequency (KHz)"); plt.grid()
        
        fig.add_subplot(2,2,2)
        plt.plot(Y_norm , 'o')
        plt.xticks(np.arange(0, len(Y_norm), step=1)); #plt.yticks(np.arange(0, 5.1, step=0.25))
        plt.xlabel("Frequency number (f_i)"); plt.ylabel("Manitud dB"); plt.grid()
        
        fig.add_subplot(2,2,3)
        plt.plot(X_norm , 'o')
        plt.xticks(np.arange(0, len(X_norm ), step=1)); plt.yticks(np.arange(0, np.max(X_norm ), step=0.5))
        plt.xlabel("Frequency number (f_i)"); plt.ylabel("Frequency/F_0"); plt.grid()
        
        fig.add_subplot(2,2,4)
        plt.plot(Y_peaks_diff, 'o')
        plt.xticks(np.arange(0, len(X_peaks), step=1)); #plt.yticks(np.arange(0, 5.1, step=0.25))
        plt.xlabel("Frequency number (f_i)"); plt.ylabel("Differences (dB)"); plt.grid()
        
        plt.savefig(temp_folder + "Harmonics Analysis - " + file[:-4])
        plt.show()
        plt.close()
        