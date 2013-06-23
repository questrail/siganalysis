# -*- coding: utf-8 -*-
"""
siganalysis.py

Provide signal analysis routines.

"""

# TODO: Add unit tests for this code.

# Try to future proof code so that it's Python 3.x ready
from __future__ import print_function
# Importing unicode_literals broke the convoluation on line 132
# window='hanning'
#from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import

# Numerical analysis related imports
import numpy as np
import scipy
import matplotlib.pyplot as plt
# TODO(mdr): Should I split this into a package so that matplotlib
# is only imported if the plotting functions are needed?


def time_slice_zip(number_of_samples, samples_per_time_slice):
    """Create a zipped list of tuples for time slicing a numpy array

    When dealing with large numpy arrays containing time series data, it is
    often desirable to time slice the data on a fixed duration, such as one
    minute. This function creates a list of tuples (similar to the Python zip
    function) to iterate through a numpy array using slices.

    Args:
        number_of_samples: Number of samples in the time series numpy array
        samples_per_time_slice: Desired number of samples per time slice not
            including the last time slice which will be limited to the length
            of the time series

    Returns:
        A list of tuples that can be used to time slice the data.

    """
    current_index = 0
    zipped = []
    while current_index < (number_of_samples - samples_per_time_slice):
        this_tuple = current_index, current_index + samples_per_time_slice
        zipped.append(this_tuple)
        current_index += samples_per_time_slice
    zipped.append((current_index, number_of_samples))
    return zipped


def stft(input_data, sampling_frequency_hz, frame_size_sec, hop_size_sec,
         use_hamming_window=False):
    """Calculates the Short Time Fourier Transform

    Using code based on http://stackoverflow.com/a/6891772/95592 calculate
    the STFT.

    Args:
        input_data: A 1D numpy ndarray containing the signal in the time
            domain that will be converted to the freq domain via STFT.
        sampling_frequency_hz: Sampling frequency originally used to capture
            the input_data
        frame_size_sec: Frame size given in seconds. The frame size determines
            how long each FFT will be in the time domain.
        hop_size_sec: Hop size given in seconds. The hop size is the time
            by which the frame should be shifted forward for the next
            FFT. It is not uncommon for this to be less than the frame
            size so that there is some amount of overlap.
        use_hamming_window: A Boolean indicating if the Hamming window
            should be used when performing the FFT. Using a Hamming window
            helps.

    Returns:
        A tuple containing:
            1. A 2D numpy ndarray providing the amplitude of the STFT with
                respect to the frequency and time having a shape of
                (time, freq). This array is trimmed to be single-sided instead
                of returning the double-sided FFT, and it is normalized by
                2/N where N is the length of the frequency domain info. The
                DC component is not multiplied by 2 though, it is just
                normalized by 1/N.
            2. A 1D numpy ndarray [shape = (time,)] containing the time in
                seconds for each value in the stft_data along the time axes.
            3. A 1D numpy ndarray [shape = (freq,)] containing the freq in
                Hz for each value in the stft_data along the frequency axes.
            4. A float indicating the frequency bin size in Hz or what is
                also referred to as the frequency domain step size (not
                to be confused with or equal to the sampling frequency).

    """
    # FIXME(mdr) The windowing functionality is currently broken!!! Using a
    # Hamming window will result in amplitudes which are approximately 1/2
    # what they should be. Therefore, regardless of what is passed in this
    # function, I'm turning off the windowing.
    use_hamming_window = False

    # TODO(mdr): The Agilent 35670A uses a Hann (aka Hanning) window, which
    # is slightly different from a Hamming window. I should change this code
    # so that the user can select which type of window is used.

    num_frame_samples = int(frame_size_sec * sampling_frequency_hz)
    num_hop_samples = int(hop_size_sec * sampling_frequency_hz)
    # TODO(mdr): Add a log entry here showing the frame size in seconds
    # and samples, as well as the hop size in seconds and samples

    # TODO(mdr): Add a log entry that lists which window was used
    if (use_hamming_window):
        X = np.array([
            scipy.fft(
                scipy.hamming(num_frame_samples) *
                input_data[i:i+num_frame_samples])
            for i in range(
                0,
                len(input_data)-num_frame_samples,
                num_hop_samples)])
    else:
        X = np.array([
            scipy.fft(input_data[i:i+num_frame_samples])
            for i in range(
                0,
                len(input_data)-num_frame_samples,
                num_hop_samples)])

    # Normalize the FFT results
    # See "Description and Application of Fourier Transforms and Fourier
    # Series" rev A05 by Matthew Rankin for a description on why the
    # normalization is 2 / N except for the DC component which is 1 / N
    # Only deal with the single-sided FFT, so cut it in half
    X = X[:, :num_frame_samples//2]
    # Convert from complex to absolute values
    X = np.abs(X)
    # Divide all components by the num_frame_samples
    # Multiply all but the DC component by 2
    non_dc_normalization = 2 / num_frame_samples
    X[:, 1:] = X[:, 1:] * non_dc_normalization
    X[:, 0] = X[:, 0] / num_frame_samples

    # Create the time vector
    time_vector_stft = np.linspace(
        frame_size_sec / 2,
        (X.shape[0] - 1) * hop_size_sec + frame_size_sec / 2,
        X.shape[0])

    # Calculate the width of each frequency bin
    hz_per_freq_bin = sampling_frequency_hz / num_frame_samples

    # Create the frequency vector
    freq_vector_stft = np.arange(X.shape[1]) * hz_per_freq_bin

    return (X, time_vector_stft, freq_vector_stft, hz_per_freq_bin)


def hz2khz(frequency_in_hz):
    """Convert from Hz to kHz

    Args:
        frequency_in_hz: A float containing the frequency value in Hz
            that is to be converted.

    Return:
        The frequency in kHz.

    """
    return frequency_in_hz / 1000


def smooth(x, window_len=11, window='hanning'):
    """smooth the data using a window with requested size.

    cookb_signalsmooth.py

    from: http://scipy.org/Cookbook/SignalSmooth

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    Args:
        x: The input signal to be smoothed
        window_len: the dimension of the smoothing window
        window: The type of window from 'flat', 'hanning', 'hamming',
            'bartlett', 'blackman' flat window will produce a moving
            average smoothing.

    Returns:
        the smoothed signal

    example:

    import numpy as np
    t = np.linspace(-2,2,0.1)
    x = np.sin(t)+np.random.randn(len(t))*0.1
    y = smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman,
    numpy.convolve, scipy.signal.lfilter

    """
    # TODO(mdr): Convert the above example to a doctest or remove.
    # TODO(mdr): The window parameter could be an array containing the window
    # itself, instead of being a string that calls out the window to use.

    if x.ndim != 1:
        raise ValueError('Function smooth only accepts 1D arrays.')

    if x.size < window_len:
        raise ValueError('Input vector needs to be bigger than window size.')

    if window_len < 3:
        return x

    if window_len & 1:
        pass
    else:
        window_len += 1

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window must be one of: 'flat', 'hanning', "
                         "'hamming', 'bartlett', 'blackman'")

    s = np.r_[x[window_len-1:0:-1], x, x[-1:-window_len:-1]]

    if window == 'flat':
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')
    y = np.convolve(w/w.sum(), s, mode='valid')
    samples_to_strip = (window_len - 1) / 2
    return y[samples_to_strip:len(y)-samples_to_strip]


def smooth2(x, beta=3, window_len=11):
    """Smooth function using Kaiser window

    Args:
        x: ndarray containing the signal to be smoothed
        beta: beta to use as part of the Kaiser smoothing
        window_len: Integer length of window to be used in Kaiser
            smoothing, which must be odd or it will be made odd.

    Returns:
        An ndarrary containing the smoothed signal.

    """
    # If window_len is not odd, add one so that it is odd
    if window_len & 1:
        pass
    else:
        window_len += 1
    s = np.r_[x[window_len-1:0:-1], x, x[-1:-window_len:-1]]
    w = np.kaiser(window_len, beta)
    y = np.convolve(w/w.sum(), s, mode='valid')
    print("Window length =", window_len)
    print("y length =", len(y))
    samples_to_strip = (window_len - 1) / 2
    return y[samples_to_strip:len(y)-samples_to_strip]


def plot_spectrogram(stft_data, start_plot_freq, stop_plot_freq, hz_per_freq_bin,
                      axes, time_vector, title, xlabel, ylabel,
                      colorbar_label=False, colorbar_size=8):
    """Create a spectrogram plot

    Take a numpy ndarray containing amplitude vs. frequency vs. time info and
    create a spectrogram. Currently, this assumes that the stft_data starts at
    0 Hz and uses the given hz_per_freq_bin. It would be better if I passed in
    a freq array similar to the time_array that is passed.

    Args:
        stft_data: A 2D numpy ndarray of shape (time, freq) containing the
            amplitude over both freq and time.
        start_plot_freq: Starting spectrogam plot frequency in Hz
        stop_plot_freq: Ending spectrogram plot frequency in Hz
        hz_per_freq_bin: Width of each freq bin in Hz
        axes: matplotlip axes that this plot should be added to
        time_vector: 1D np.ndarray containing the time in sec of the stft_data
        title: Title to be added to the plot
        xlabel: matplotlib x-axis label to be added to the plot
        ylabel: matplotlib y-axis label to be added to the plot
        colorbar_label: matplotlib label to be added to the colorbar. If
            excluded then the colorbar is not plotted.
        colorbar_size: Integer of the colorbar font size.

    Returns:
        matplolib handle to the spectrogram
    """
    # FIXME(mdr): I need to refactor this code to use a freq array as an input
    # instead of using hz_per_freq_bin, since that inherently assumes the
    # stft_data starts at 0 Hz, which might not always be the case.
    # FIXME(mdr): Need to add a start and stop plot for the time domain as well.
    # TODO(mdr): I could convert the start/stop freq and time args to be tuples
    # with one for time and another for freq.

    # Determine the frequency bins for the start and stop freqs
    start_freq_bin = int(start_plot_freq / hz_per_freq_bin)
    stop_freq_bin = int(stop_plot_freq / hz_per_freq_bin)
    # Create the spectrogram
    spectrogram = axes.imshow(stft_data[:,start_freq_bin:stop_freq_bin].T,
            origin='lower', aspect='auto', interpolation='nearest')
    if colorbar_label:
        cb = plt.colorbar(spectrogram, ax=axes)
        cb.ax.tick_params(labelsize=colorbar_size)
        cb.set_label(colorbar_label)
    spectrogram.set_extent([time_vector[0], time_vector[-1],
        start_plot_freq, stop_plot_freq])
    axes.set_title(title)
    axes.set_xlabel(xlabel)
    axes.set_ylabel(ylabel)
    return spectrogram


def calculate_peak_hold(stft_data, frequency_array):
    """Calculate the peak hold for a given STFT dataset.

    Args:
        stft_data: A 2D numpy ndarray with shape (time, freq) containing
            the amplitude vs freq vs time.
        frequency_array: A 1d numpy ndarray containing the frequencies
            for the stft_data.

    Returns:
        peak_hold: A 2D numpy structured array containing the frequency
            and amplitude as the dtype.
    """
    data_type = np.dtype([('frequency', 'f8'), ('amplitude', 'f8')])
    peak_hold = np.zeros(frequency_array.size, dtype=data_type)
    peak_hold['frequency'] = frequency_array
    peak_hold['amplitude'] = np.amax(stft_data, axis=0)
    return peak_hold

# Plot peak hold
