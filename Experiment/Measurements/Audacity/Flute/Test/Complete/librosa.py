#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import librosa
import librosa.display
import numpy as np

path = '/home/sebas/Documents/Recorders/WavFiles/'
file = 'Flute-F.wav'

y, sr = librosa.load(path + file)
librosa.feature.melspectrogram(y=y, sr=sr)


D = np.abs(librosa.stft(y))**2
S = librosa.feature.melspectrogram(S=D, sr=sr)

# Display of mel-frequency spectrogram coefficients, with custom
# arguments for mel filterbank construction (default is fmax=sr/2):

# Passing through arguments to the Mel filters
S = librosa.feature.melspectrogram(y=y, sr=sr, n_mels=128, fmax=sr/2)

import matplotlib.pyplot as plt
plt.figure(figsize=(10, 4))
S_dB = librosa.power_to_db(S, ref=np.max)
librosa.display.specshow(S_dB, x_axis='time',
                         y_axis='mel', sr=sr,
                         fmax=sr/2)
plt.colorbar(format='%+2.0f dB')
plt.title(file)
#plt.tight_layout()
#plt.show()
plt.savefig(path+file[:-4]+'-spectro')