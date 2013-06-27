# siganalysis

Provide various analysis routines required for analyzing signals in
Python. Functions include:

* Calculating Short-Term Fourier Transform
* Smoothing a signal
* Plotting an STFT's spectrogram
* Calculating the peak hold of an STFT in the freq domain
* Plotting the peak hold of an STFT

The above functions are handy when analyzing signals recorded in the
time domain, such as using a TEAC LX-10 data recorder, and seeing the
frequency spectrum. This is usefull for Electromagnetic Compatibiliity
(EMC) analyses.

## Requirements

* [numpy][]
* [scipy][]
* [matplotlib][]

[numpy]: http://www.numpy.org
[scipy]: http://www.scipy.org
[matplotlib]: http://matplotlib.org
