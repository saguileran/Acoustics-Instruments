# -*- coding: utf-8 -*-
import scipy.io.wavfile, matplotlib
import numpy as np
import matplotlib.pyplot as plt
from numpy import fft as fft
import librosa, librosa.display
from scipy.signal import find_peaks

#Necesary functions
def ReduceVectors(x, y, first = 25):
    divisores, i = x.shape[0]/np.array(range(2,100))%1, first #primer divisor
    while divisores[i] != 0: i += 1
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
        Y = X
        n = np.arange(0, len(Y))
    f = n * resolution
    return f, Y

def Peaks(coefficients_weight, frequencies, distance=150, height=-20):
     peaks, _ = find_peaks(coefficients_weight[:samples//2], distance=distance, height=height)
     #return peaks, frequencies[peaks]
     return frequencies[peaks], coefficients_weight[peaks]

def Smooth(x,window_len=11,window='hanning'):
    dim = len(x)
    if x.ndim != 1: raise ValueError( "smooth only accepts 1 dimension arrays.")
    if x.size < window_len: raise ValueError( "Input vector needs to be bigger than window size.")

    if window_len<3: return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError ("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s = np.r_[x[window_len-1:0:-1],x,x[-2:-window_len-1:-1]]
    
    if window == 'flat': w = np.ones(window_len,'d')#moving average
    else: w=eval('np.'+window+'(window_len)')

    y = np.convolve(w/w.sum(),s,mode='valid')
    smooted = y[(window_len//2-1):-(window_len//2)]
    if dim == len(smooted):return smooted
    else: return smooted[:dim]
    
#a temp folder for downloads
temp_folder = '/home/sebas/Documents/Acoustics-Instruments/Experiment/Measurements/Audacity/Flute/Lina/C.wav'
#file = '/C.wav'
#temp_folder = "/Users/sebas/Downloads/BlackClover.wav"

#read wav file
rate, audData = scipy.io.wavfile.read(temp_folder)
samples = audData.shape[0]
chanels = len(audData.shape)
#audData = 
#if chanels==2: audData = np.sum(audData.astype(float), axis=1, keepdims= True)/2
#else: audData = audData.reshape(samples, 1)

#energy = np.sum(audData.astype(float)**2)
#power - energy per unit of time
#1.0/(2*(audData.size)+1)*np.sum(audData.astype(float)**2)/rate

#create a time variables to plot
time = np.arange(0, float(samples), 1) / rate
fourier = fft.fft(audData) 
fourier_scaled = fourier / float(samples)
frequencies, coefficients_weight = dft_map(fourier_scaled , rate, shift=False)
coefficients_weight = abs(coefficients_weight)
Maxfretoplot = np.where(frequencies//1 == 5000)[0][0]


# scale by the number of points so that the magnitude does not depend on the length
fourier1 = fourier[0:(samples//2)] / float(samples)

#calculate the frequency at each point in Hz
freqArray = np.arange(0, (samples//2), 1.0) * (rate*1.0/samples);

#Calculate the SPL for each point
SPL = 20*np.log10(audData)
p0, c = 1, 20
pressure = p0*10**(SPL/20)

#plt.plot(freqArray/1000, 10*np.log10(fourier1), color='#ff7f00', linewidth=0.02)
#plt.xlabel('Frequency (kHz)')
#plt.ylabel('Power (dB)')

fig = plt.figure(figsize=(40, 24))
plt.subplots_adjust(hspace = 0.4, wspace = 0.2)

fig.add_subplot(4,2,1)
plt.title('Waveform')
plt.plot(time, audData, linewidth=0.1, alpha=0.7, color='red')
plt.xticks(np.arange(0, time[-1], step=0.25))
plt.xlabel('Time (s)')
plt.ylabel('Amplitude')

fig.add_subplot(4,2,5)
plt.title('Amplitude Coefficients of FFT')
plt.xlabel('Frequencies (Hz)')
plt.ylabel('Weight')
plt.xticks(np.arange(0, 2500, step=250))
plt.xlim(0, 2500)
plt.plot(frequencies, coefficients_weight)
#plt.plot(audData, np.abs(fourier[0:(samples//2)]), color='red')
#plt.xlabel('k')
#plt.ylabel('Amplitude')

fig.add_subplot(4,2,3)
plt.title('dB')
#Time, SPLR = ReduceVectors(time, SPL, 10)
plt.plot(time, Smooth(SPL, 25), color='red')
plt.xlim([2, 3])
plt.xlabel('time (s)')
plt.ylabel('dB')

fig.add_subplot(4,2,2)
plt.title('Pure Tone')
x, y = freqArray/1000,  20*np.log10(fourier1)
X_pro, Y_pro = ReduceVectors(x, y)
#plt.plot(X_pro, abs(Y_pro), color='#ff7f00', linewidth=0.5)
plt.plot(x, y, color='#ff7f00', linewidth=0.5)
plt.xticks(np.arange(0, 10, step=0.25))
#plt.yticks(np.arange(-50, 50, step=10))
plt.xlabel('Frequency (kHz)')
plt.ylabel('Power (dB)')
plt.xlim(0, 10/2)


fig.add_subplot(4,2,4)
plt.title('Smooth Tone and Peaks')

Y_pro = Smooth(Y_pro, 30)
x = (1/1000)*np.arange(0, (len(Y_pro)), 1.0) * (rate*1.0/len(Y_pro));
X_peaks, Y_peaks = Peaks(Y_pro, X_pro, distance=1)

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
plt.xlabel('Frequency (kHz)')
plt.ylabel('Power (dB)')
plt.grid()
plt.xticks(np.arange(0, 10, step=0.25))
plt.xlim(0, 10/2)
plt.ylim(np.min(Y_pro.real), np.max(Y_pro.real)+10)




'''
ax = fig.add_subplot(4,2,6)
plt.title('Peaks')
X_peaks, Y_peaks = Peaks(Y_pro, X_pro, distance=1)
#plt.xticks(np.arange(0, 10, step=0.25))
#plt.grid()
plt.xticks(np.arange(0, 10, step=0.25))
plt.xlim(0, 10/2)

plt.scatter(X_peaks , Y_peaks, s=100)

for x,y in zip(X_peaks, Y_peaks):
    #label = "({:.2f}, {:.2f})".format(1000*x,y)
    label = "{:.2f}".format(1000*x)
    plt.annotate(label, # this is the text
                 (x,y), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,15), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center
#for xy in zip(X_peaks , Y_peaks): 
#    print(xy)                                      # <--
#    ax.annotate('(%s, %s)' % truncate(xy, 2), xy=xy, textcoords='data') # <--

plt.xlabel('Frequency (kHs)')
plt.ylabel('Power (dB)')
'''

fig.add_subplot(4,2,6)
plt.title('Pressure')
T, P = ReduceVectors(time, pressure/c, 1)
plt.loglog(T, P, color='red')
plt.xlim([10**-3, 10])
plt.xlabel('time (s)')
plt.ylabel('Pressure (Pa)')
plt.xlim(10**-2, 10)


fig.add_subplot(4,2,7)
'''
Pxx, freqs, bins, im = plt.specgram(audData, Fs=rate, NFFT=1024, cmap=plt.get_cmap('autumn_r'))
cbar=plt.colorbar(im)
plt.xlabel('Time (s)')
plt.ylabel('Frequency (Hz)')
cbar.set_label('Intensity dB')
'''
y, sr = librosa.load(temp_folder)
librosa.feature.melspectrogram(y=y, sr=sr)
D = np.abs(librosa.stft(y))**2
S = librosa.feature.melspectrogram(S=D, sr=sr)
S = librosa.feature.melspectrogram(y=y, sr=sr, n_mels=128, fmax=sr/2)
S_dB = librosa.power_to_db(S, ref=np.max)
librosa.display.specshow(S_dB, x_axis='time',y_axis='mel', sr=sr, fmax=sr/2)
#plt.yticks(np.arange(0, 10000, step=1000))
plt.colorbar(format='%+2.0f dB')
plt.title('Spectrogram')



plt.show()
plt.close()


fig = plt.figure(figsize=(20, 20))
plt.subplots_adjust(hspace = 0.2, wspace = 0.2)


fig.add_subplot(2,2,1)
plt.plot(X_peaks, 'o')
plt.xticks(np.arange(0, len(X_peaks), step=1))
plt.yticks(np.arange(0, np.max(X_peaks)+0.1, step=0.25))
plt.xlabel("Frequency number (f_i)")
plt.ylabel("Frequency (KHz)")
plt.grid()

fig.add_subplot(2,2,2)
Y_norm = Y_peaks.real - np.min(Y_peaks.real)
plt.plot(Y_norm , 'o')
plt.xticks(np.arange(0, len(Y_norm), step=1))
#plt.yticks(np.arange(0, 5.1, step=0.25))
plt.xlabel("Frequency number (f_i)")
plt.ylabel("Manitud dB")
plt.grid()

fig.add_subplot(2,2,3)
index_main_frequecy = np.where(np.max(Y_peaks.real) == Y_peaks.real)
X_norm = X_peaks/X_peaks[index_main_frequecy]
plt.plot(X_norm , 'o')
plt.xticks(np.arange(0, len(X_norm ), step=1))
plt.yticks(np.arange(0, np.max(X_norm ), step=0.5))
plt.xlabel("Frequency number (f_i)")
plt.ylabel("Frequency/F_0")
plt.grid()

fig.add_subplot(2,2,4)
Y_peaks_diff = [abs(Y_peaks[i+1] - Y_peaks[i]) for i in range(len(Y_peaks)-1)]
#Y_norm_fraction = Y_peaks/Y_peaks[index_main_frequecy]
#Y_norm_fraction= Y_norm_fraction**-1
plt.plot(Y_peaks_diff, 'o')
plt.xticks(np.arange(0, len(X_peaks), step=1))
#plt.yticks(np.arange(0, 5.1, step=0.25))
plt.xlabel("Frequency number (f_i)")
plt.ylabel("Differences (dB)")
plt.grid()


