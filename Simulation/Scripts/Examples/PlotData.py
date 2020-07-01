#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
author: root
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from numpy import fft as fft

filename = input('Enter file name: ')
#filename = 'Pulse10k-0mm.dat' #'0mm.dat'
data = pd.read_csv('./' + filename, 
                   sep = ',',
                   names=['t', 'rho'])
p0 = 343**2
data['rho'] = (data['rho'] -1) * p0#[:8000]
data['fft'] = fft.fft(data['rho'])

plt.plot(data['t'], data['rho'])
plt.xlabel('t'); plt.ylabel('rho')
plt.title('density')
plt.pause(1)

plt.plot(data['t'], data['fft'])
plt.pause(1)

Vrm = []

for i in range(0,len(data['rho'])):
    Vrm.append(np.sqrt(np.sum(data['rho'][:i]**2)/i))
    
SL = 20 * np.log10(Vrm)

#for i in range(len(data['rho'])):
#    if data['rho'][i] == 0: data['rho'][i] = 1

#data['dB']  = -20 * np.log10(data['rho'])
plt.plot(data['t'], Vrm)
plt.pause(1)

plt.plot(data['t'], SL)
plt.pause(1)

plt.show()

#plt.savefig(filename[:-4])
