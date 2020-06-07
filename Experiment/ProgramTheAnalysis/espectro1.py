import random as rn
import scipy.stats as st
import math
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
#::::::::::::::::::::::::::::::::::::
#:::::::: ESPECTRO EXPERIMENTAL::::::
#::::::::::::::::::::::::::::::::::::

from scipy.optimize import curve_fit
from scipy.signal import find_peaks
from astropy.io import ascii



datos=ascii.read('espectro1.txt', data_start=0) #Acá se debe poner el nombre 
#del archivo del cual se desean extraer los datos. La divición de decimales 
#debe estar con "." no con ",". 

g=len(datos)

frecuencia=np.zeros(g)
decibeles=np.zeros(g)
for i in range(0,g):
    frecuencia[i]=datos[i][0]
    decibeles[i]=datos[i][1]
    

fig, axs=plt.subplots(1,1,sharey=False)


fig.suptitle('ESPECTRO', fontsize=18)
plt.xlabel('FRECUENCIAS', fontsize=14)
plt.ylabel('DECIBELES', fontsize=14)
axs.plot(frecuencia,decibeles,color='purple')
plt.show()