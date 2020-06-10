This program generate acoustic grpahs of all wav files of a folder, it create Sound Pressure Level (SPL), Fourier Transform, pure data, Pressure, Spectrogram and Spectral plots for every wav file.
In addition, it create a plot of harmonics making visible lineatity. 

In this folder you will find three scripts:

*  **SpectrogramPlot.py:**
Generate just three graphs, spectrogram, spectral and peaks plots. It create these three plots for all wav files in a folder and save in the same.

*  **PlottingFeatures.py:**
Generate all plots describe before wut is make just for one file.

*  **AllPlots.py:**
  Generate all plots for all wav files of a folfer. It is the final version.
  
  
We suggest that run AllPlots script to get all the information and analysis, to run it you need the folowing libraries

* matplotlib
* librosa
* numpy
* scipy
* pandas
 
If you are using linux machine it can be done executing: **sudo apt-get install python-scipy python-matplotlib python-numpy python-librosa python-pandas**.
On windows you run **pip install librosa scipy matplotlib numpy librosa pandas** in any anaconda (spyder) terminal. 

The most easy way is install anaconda and make it from graphical interface.
To run you enter **python  AllPlots.py** in terminal, if you can measure the time spent add time at begining. 
After enter this code you have to give the folder direction, it has to end in / in other case the script doesn't work.

It is tasted in linux SO, windows will be tested soon.
