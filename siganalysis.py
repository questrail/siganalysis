# -*- coding: utf-8 -*-
# Copyright (c) 2013-2016 The siganalysis developers. All rights reserved.
# Project site: https://github.com/questrail/siganalysis
# Use of this source code is governed by a MIT-style license that
# can be found in the LICENSE.txt file for the project.
"""Provide Python routines for signal analysis

Provide various analysis routines required for analyzing signals in Python,
such as calculating a Short-Time Fourier Transform, plotting an STFT's
spectrogram, calculating the peak hold values for an STFT, etc.
"""

# Try to future proof code so that it's Python 3.x ready
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import division
from __future__ import absolute_import

# Numerical analysis related imports
import numpy as np
import scipy
import matplotlib.pyplot as plt

__version__ = '0.4.0'


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
         use_hamming_window=True):
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
                seconds for each value in the stft_data along the time axis.
            3. A 1D numpy ndarray [shape = (freq,)] containing the freq in
                Hz for each value in the stft_data along the frequency axis.
            4. A float indicating the frequency bin size in Hz or what is
                also referred to as the frequency domain step size (not
                to be confused with or equal to the sampling frequency).

    """
    num_frame_samples = int(frame_size_sec * sampling_frequency_hz)
    num_hop_samples = int(hop_size_sec * sampling_frequency_hz)
    if (use_hamming_window):
        x = np.array([
            scipy.fft(
                2 * scipy.hamming(num_frame_samples) *
                input_data[i:i+num_frame_samples])
            for i in range(
                0,
                len(input_data)-num_frame_samples,
                num_hop_samples)])
    else:
        x = np.array([
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
    x = x[:, :num_frame_samples//2]
    # Convert from complex to absolute values
    x = np.abs(x)
    # Divide all components by the num_frame_samples
    # Multiply all but the DC component by 2
    non_dc_normalization = 2 / num_frame_samples
    x[:, 1:] = x[:, 1:] * non_dc_normalization
    x[:, 0] = x[:, 0] / num_frame_samples

    # Create the time vector
    # FIXME(mdr): Need to add test to make sure this is correctly calculated.
    # Might want to refactor into separate function.
    time_vector_stft = np.linspace(
        0,
        (x.shape[0] - 1) * hop_size_sec,
        x.shape[0])

    # Calculate the width of each frequency bin
    hz_per_freq_bin = sampling_frequency_hz / num_frame_samples

    # Create the frequency vector
    freq_vector_stft = np.arange(x.shape[1]) * hz_per_freq_bin

    return (x, time_vector_stft, freq_vector_stft, hz_per_freq_bin)


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
    if x.ndim != 1:
        raise ValueError('Function smooth only accepts 1D arrays.')

    if x.size < window_len:
        raise IndexError('Input vector needs to be bigger than window size.')

    if window_len < 3:
        return x

    if window_len & 1:
        pass
    else:
        window_len += 1

    if window not in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
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
    samples_to_strip = (window_len - 1) / 2
    return y[samples_to_strip:len(y)-samples_to_strip]


def calculate_peak_hold(stft_data, frequency_array):
    """Calculate the peak hold for a given STFT dataset.

    Args:
        stft_data: A 2D numpy ndarray with shape (time, freq) containing
            the amplitude vs freq vs time.
        frequency_array: A 1d numpy ndarray containing the frequencies
            for the stft_data.

    Returns:
        peak_hold: A 1D numpy structured array containing the frequency
            and amplitude with the dtype [(freq, amp)]

    Raises:
        ValueError: The frequency_array and stft_data[1] are not the same
            length.
    """
    if frequency_array.size != stft_data.shape[1]:
        raise IndexError('The size of the frequency_array does not match '
                         'the STFT data.')
    data_type = np.dtype([('frequency', 'f8'), ('amplitude', 'f8')])
    peak_hold = np.zeros(frequency_array.size, dtype=data_type)
    peak_hold['frequency'] = frequency_array
    peak_hold['amplitude'] = np.amax(stft_data, axis=0)
    return peak_hold


def plot_spectrogram(stft_data,
                     time_vector,
                     freq_vector,
                     plot_axis,
                     freq_plot_range=False,
                     time_plot_range=False,
                     plot_title=False,
                     plot_xlabel=False,
                     plot_ylabel=False,
                     colorbar_label=False,
                     colorbar_fontsize=8):
    """Create a spectrogram plot

    Take a numpy ndarray containing amplitude vs. frequency vs. time info and
    create a spectrogram. Currently, this assumes that the stft_data starts at
    0 Hz and uses the given hz_per_freq_bin. It would be better if I passed in
    a freq array similar to the time_array that is passed.

    Args:
        stft_data: A 2D numpy ndarray of shape (time, freq) containing the
            amplitude over both freq and time.
        time_vector: A 1d numpy ndarray containing the time in seconds for each
            value in the stft_data along the time axis. time_vector is assumed
            to be sorted and to contain equal time steps.
        freq_vector: A 1d numpy ndarray containing the freq in Hz for each
            value in the stft_data along the frequency axis. freq_vector is
            assumed to be sorted and to contain equal frequency steps.
        plot_axis: matplotlip axis that this plot should be added to
        freq_plot_range: A tuple containing the start and stop frequency in Hz
            for the spectrogram plot (frequencies are inclusive)
        time_plot_range: A tuple containing the start and stop time in seconds
            for the spectrogram plot (time are inclusive)
        plot_title: An optional string with the plot title
        plot_xlabel: An optional string with the x-axis label
        plot_ylabel: An optional string with the y-axis label
        colorbar_label: An optional string with the label to be added to the
            colorbar. If excluded then the colorbar is not plotted.
        colorbar_fontsize: Integer of the colorbar font size.

    Returns:
        matplolib handle to the spectrogram
    """
    if freq_plot_range is False:
        start_freq_plot = freq_vector[0]
        stop_freq_plot = freq_vector[-1]
    else:
        start_freq_plot, stop_freq_plot = freq_plot_range

    # FIXME: Is there an error in the time plot range or the calculation of the
    # start and stop time bins?
    if time_plot_range is False:
        start_time_plot = time_vector[0]
        stop_time_plot = time_vector[-1]
    else:
        start_time_plot, stop_time_plot = time_plot_range
    # Calculate the hz_per_freq_bin assuming that the frequency steps are
    # equal.
    hz_per_freq_bin = freq_vector[1] - freq_vector[0]
    sec_per_time_bin = time_vector[1] - time_vector[0]
    # Determine the frequency bins for the start and stop freqs
    start_freq_bin = int((start_freq_plot - freq_vector[0]) / hz_per_freq_bin)
    stop_freq_bin = int((stop_freq_plot - freq_vector[0]) / hz_per_freq_bin)
    start_time_bin = int((start_time_plot - time_vector[0]) / sec_per_time_bin)
    stop_time_bin = int((stop_time_plot - time_vector[0]) / sec_per_time_bin)
    # Create the spectrogram
    spectrogram = plot_axis.imshow(
        stft_data[start_time_bin:stop_time_bin,
                  start_freq_bin:stop_freq_bin].T,
        origin='lower',
        aspect='auto',
        interpolation='nearest')
    if colorbar_label:
        cb = plt.colorbar(spectrogram, ax=plot_axis)
        cb.ax.tick_params(labelsize=colorbar_fontsize)
        cb.set_label(colorbar_label)
    spectrogram.set_extent([start_time_plot, stop_time_plot,
                           start_freq_plot, stop_freq_plot])
    if plot_title:
        plot_axis.set_title(plot_title)
    if plot_xlabel:
        plot_axis.set_xlabel(plot_xlabel)
    if plot_ylabel:
        plot_axis.set_ylabel(plot_ylabel)
    return spectrogram


def plot_peak_hold(axis,
                   stft_data,
                   frequency_array,
                   title=False,
                   xlabel=False,
                   ylabel=False,
                   plot_freq_limits=False,
                   plot_amp_limits=False,
                   limit_array=False,
                   trace_label=False):
    """Plot the peak hold for a 2D STFT array

    Args:
        axis: matplotlip axis that this plot should be added to
        stft_data: A 2D numpy ndarray of shape (time, freq) containing the
            amplitude over both freq and time.
        frequency_array: A 1D numpy ndarray containing hte frequencies in
            Hz of the stft_data.
        title: An optional title to be added to the plot
        xlabel: An optional x-axis label to be added to the plot
        ylabel: An optional y-axis label to be added to the plot
        plot_freq_limits: An optional tuple containing the starting and ending
            frequencies to be used in the plot
        limit_array: An optional 1D numpy ndarray containing the limits for the
            plotted data of dtype = [('frequency', 'f8'), ('amplitude', 'f8')]

    Returns:
        matplolib handle to the axis

    Raises:
    """
    peak_hold = calculate_peak_hold(stft_data, frequency_array)
    if trace_label is not False:
        axis.loglog(peak_hold['frequency'],
                    peak_hold['amplitude'],
                    label=trace_label)
    else:
        axis.loglog(peak_hold['frequency'],
                    peak_hold['amplitude'])
    if limit_array is not False:
        axis.loglog(limit_array['frequency'],
                    limit_array['amplitude'])
    if plot_freq_limits is not False:
        axis.set_xlim(plot_freq_limits)
    if plot_amp_limits is not False:
        axis.set_ylim(plot_amp_limits)
    if title is not False:
        axis.set_title(title)
    if xlabel is not False:
        axis.set_xlabel(xlabel)
    if ylabel is not False:
        axis.set_ylabel(ylabel)
    axis.xaxis.set_major_formatter(plt.FormatStrFormatter('%g'))
    axis.yaxis.set_major_formatter(plt.FormatStrFormatter('%g'))
    axis.grid(b=True, which='major', color='0.25', linestyle='-')
    axis.grid(b=True, which='minor', color='0.75', linestyle='-')
    axis.set_axisbelow(True)


def single_frequency_over_time(stft_data,
                               freq_array,
                               time_array,
                               frequency):
    """Determine the amplitude vs. time for a particular frequency

    Given an STFT data array and its supporting frequency and time arrays, as
    well as a desired frequency, determine the amplitude for just that
    frequency.

    Args:
        stft_data: A 2D numpy ndarray containing the amplitude vs. frequency
            vs. time from a Short-Time Fourier Transform.
        freq_array: A 1D numpy ndarray containing the frequencies in Hz for the
            given STFT data.
        time_array: A 1D numpy ndarray containing the time values in seconds
            for the given STFT data.
        frequency: A float or int of the desired frequency

    Returns:
        A 1D numpy structured array of dtype
            [('time', 'f8'), ('amplitude', 'f8')]

    Raises:
        IndexError: The size of the STFT does not match the given frequency
            and/or time arrays.
    """
    # Check that the arrays are the correct size
    if freq_array.size != stft_data.shape[1]:
        raise IndexError('The size of the freq_array does not match '
                         'the STFT data.')
    if time_array.size != stft_data.shape[0]:
        raise IndexError('The size of the time_array does not match '
                         'the STFT data.')
    # Create the array to return the time and amplitude
    data_type = np.dtype([('time', 'f8'), ('amplitude', 'f8')])
    stft_at_frequency = np.zeros(time_array.size, dtype=data_type)
    stft_at_frequency['time'] = time_array
    freq_bin = int(frequency / (freq_array[1] - freq_array[0]))
    stft_at_frequency['amplitude'] = stft_data[:, freq_bin]
    return stft_at_frequency
