#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
author: root
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from numpy import fft as fft
import matplotlib

rootfolder = input('Enter folder direction: ')
filename = input('Enter file name: ')
#filename = '/home/sebas/Documents/Acoustics-Instruments/Simulation/Scripts/Examples/60mm.dat'
data = pd.read_csv(rootfolder + filename,
                   sep = ',',
                   names=['t', 'rho'])
p0 = 343**2
#data['rho'] = (data['rho'] -1) #* p0#[:8000]
data['fft'] = fft.fft(data['rho'])


#----------Plot settings---------------
fig = plt.figure(figsize=(22, 10))
#fig.suptitle('Z = ' + str(AT[j]))
matplotlib.rcParams.update({'font.size': 26})
plt.subplots_adjust(hspace = 1.2)
    
fig.add_subplot(1,2,1)
plt.plot(data['t'], data['rho'])
plt.xlabel('t'); plt.ylabel('rho')
plt.title('density')
#plt.pause(1)

fig.add_subplot(1,2,2)
plt.plot(data['t'], np.log10(abs(data['fft'])))
plt.title('Fourier Transform')
plt.xlabel('t');
plt.xlim(0,5000)

plt.savefig(rootfolder + filename[:-4])
#Vrm = []

#for i in range(0,len(data['rho'])):
#    Vrm.append(np.sqrt(np.sum(data['rho'][:i]**2)/i))
    
SL = 20 * np.log10(data['rho'])

#for i in range(len(data['rho'])):
#    if data['rho'][i] == 0: data['rho'][i] = 1

#data['dB']  = -20 * np.log10(data['rho'])
#plt.plot(data['t'], Vrm)
#plt.pause(1)

#plt.plot(data['t'], SL)
#plt.pause(1)

#plt.show()

#plt.savefig(filename[:-4])
