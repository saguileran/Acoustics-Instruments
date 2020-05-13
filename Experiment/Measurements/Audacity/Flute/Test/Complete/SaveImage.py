# -*- coding: utf-8 -*-
#import the pyplot and wavfile modules 
import matplotlib.pyplot as plot
from scipy.io import wavfile
import numpy as np
# Read the wav file (mono)
#path = '/home/sebas/Documents/Recorders'
file = "Flute-C#.wav"
samplingFrequency, signalData = wavfile.read(file)
dB, dim = 20*np.log10(signalData), len(signalData)
times = np.arange(0, dim/samplingFrequency, 1/samplingFrequency)
# Plot the signal read from wav file
#plot.subplot(211)
#plot.subplots(nrows=2, ncols=1, sharex=False)
plot.title( file)

plot.plot(times, dB)
plot.xlabel('Time (s)')
plot.ylabel('Amplitude (dB)')
plot.pause(1)

#plot.subplot(212)
plot.specgram(signalData,Fs=samplingFrequency)  
plot.xlabel('Time')
plot.ylabel('Frequency (Hz)')


plot.show()