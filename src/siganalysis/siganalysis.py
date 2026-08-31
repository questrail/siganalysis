# Copyright (c) 2013-2026 The siganalysis developers. All rights reserved.
# Project site: https://github.com/questrail/siganalysis
# Use of this source code is governed by a MIT-style license that
# can be found in the LICENSE.txt file for the project.
"""Provide Python (3.12+) routines for signal analysis.

Provide various analysis routines required for analyzing signals in Python,
such as calculating a Short-Time Fourier Transform, plotting an STFT's
spectrogram, calculating the peak hold values for an STFT, etc.
"""

# Standard module imports
import math
from importlib.metadata import version

import matplotlib.axes
import matplotlib.image
import matplotlib.pyplot as plt
import matplotlib.ticker

# Numerical analysis related imports
import numpy as np
import numpy.typing as npt
from scipy import signal
from scipy.fft import fft

__version__ = version("siganalysis")


def time_slice_zip(
    number_of_samples: int,
    samples_per_time_slice: int,
    minimum_samples_in_last_slice: int | None = None,
) -> list[tuple[int, int]]:
    """Create a zipped list of tuples for time slicing a numpy array.

    When dealing with large numpy arrays containing time series data, it is
    often desirable to time slice the data on a fixed duration, such as one
    minute. This function creates a list of tuples (similar to the Python zip
    function) to iterate through a numpy array using slices.

    Args:
        number_of_samples: Number of samples in the time series numpy array.
        samples_per_time_slice: Desired number of samples per time slice not
            including the last time slice which will be limited to the length
            of the time series.
        minimum_samples_in_last_slice: An optional minimum number of samples
            for the last time slice. A last slice shorter than this is folded
            into the slice before it, giving a longer last slice rather than a
            short one. Pass None, the default, to leave the last slice however
            short it falls. This cannot exceed samples_per_time_slice, since
            folding in one short slice cannot make up any more than that.

    Returns:
        A list of tuples that can be used to time slice the data.

    Raises:
        ValueError: There is less than one sample per time slice, or the
            minimum for the last slice exceeds samples_per_time_slice.

    Example:
        A number of samples just past a multiple of the slice size leaves a
        very short last slice, which is sometimes too short to process:

        >>> time_slice_zip(25, 10)
        [(0, 10), (10, 20), (20, 25)]

        Asking for at least 8 samples in the last slice folds those 5 samples
        into the slice before them:

        >>> time_slice_zip(25, 10, 8)
        [(0, 10), (10, 25)]

    """
    if samples_per_time_slice < 1:
        raise ValueError(
            f"There must be at least one sample per time slice, but "
            f"{samples_per_time_slice} were asked for."
        )
    if (
        minimum_samples_in_last_slice is not None
        and minimum_samples_in_last_slice > samples_per_time_slice
    ):
        raise ValueError(
            f"The last time slice cannot be held to a minimum of "
            f"{minimum_samples_in_last_slice} samples, which is longer than "
            f"the {samples_per_time_slice} samples in a time slice."
        )
    current_index = 0
    zipped = []
    while current_index < (number_of_samples - samples_per_time_slice):
        this_tuple = current_index, current_index + samples_per_time_slice
        zipped.append(this_tuple)
        current_index += samples_per_time_slice
    zipped.append((current_index, number_of_samples))
    if minimum_samples_in_last_slice is not None and len(zipped) > 1:
        last_start, last_stop = zipped[-1]
        if last_stop - last_start < minimum_samples_in_last_slice:
            # Fold the short slice into the one before it rather than dropping
            # it, so that the samples in it are still covered.
            zipped[-2:] = [(zipped[-2][0], last_stop)]
    return zipped


# The windows that stft() accepts, mapped to the function creating each one.
# Each function takes the number of samples and returns the window.
#
# A window trades frequency resolution, the width of the main lobe, against
# the suppression of spectral leakage, the height of the side lobes. hamming
# is what siganalysis has always applied, and hann is what the Agilent 35670A
# applies. blackman and blackmanharris widen the main lobe further to push the
# side lobes down, which helps when a small tone sits near a large one.
# flattop has the widest main lobe and the flattest top, which makes it the
# most accurate for the amplitude of a tone that falls between two bins, at
# the cost of resolving those bins.
STFT_WINDOWS = {
    "hamming": signal.windows.hamming,
    "hann": signal.windows.hann,
    # hann is also known as hanning, which is what smooth() calls it.
    "hanning": signal.windows.hann,
    "blackman": signal.windows.blackman,
    "blackmanharris": signal.windows.blackmanharris,
    "flattop": signal.windows.flattop,
}


def stft(
    input_data: npt.NDArray,
    sampling_frequency_hz: float,
    frame_size_sec: float,
    hop_size_sec: float,
    window: str | None = "hamming",
) -> tuple[npt.NDArray, npt.NDArray, npt.NDArray, float]:
    """Calculate the Short Time Fourier Transform.

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
            size so that there is some amount of overlap. Both the frame
            size and the hop size are truncated to a whole number of
            samples at the given sampling frequency.
        window: The name of the window to apply to each frame before
            transforming it, one of the keys of STFT_WINDOWS, or None to
            apply no window at all. A window reduces the spectral leakage
            from a tone that does not fall on the center of a frequency bin.
            The returned amplitudes are corrected for the gain of the window,
            so they are comparable to the amplitudes returned for an
            unwindowed frame whichever window is chosen.

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

    Raises:
        ValueError: The frame size or the hop size is shorter than one
            sample at the given sampling frequency, or the window is not
            one that STFT_WINDOWS names.
        IndexError: The input_data is not longer than one frame, so there
            is nothing to transform.

    """
    num_frame_samples = int(frame_size_sec * sampling_frequency_hz)
    num_hop_samples = int(hop_size_sec * sampling_frequency_hz)
    if num_frame_samples < 1:
        raise ValueError(
            f"A {frame_size_sec} sec frame is shorter than one sample at "
            f"{sampling_frequency_hz} Hz."
        )
    if num_hop_samples < 1:
        raise ValueError(
            f"A {hop_size_sec} sec hop is shorter than one sample at "
            f"{sampling_frequency_hz} Hz."
        )
    if len(input_data) <= num_frame_samples:
        raise IndexError(
            f"The input_data needs to be longer than the {num_frame_samples} "
            f"samples in a {frame_size_sec} sec frame, but it contains "
            f"{len(input_data)} samples."
        )
    if window is not None and window not in STFT_WINDOWS:
        raise ValueError(
            f"Window must be None or one of: "
            f"{', '.join(repr(name) for name in STFT_WINDOWS)}"
        )
    if window is None:
        # No window is a rectangular window, which leaves a frame untouched.
        window_samples = np.ones(num_frame_samples)
    else:
        # Applying a window scales every amplitude by the mean of the window,
        # its coherent gain, so divide the window through by that mean to get
        # back the amplitudes an unwindowed frame would have given.
        window_samples = STFT_WINDOWS[window](num_frame_samples)
        window_samples = window_samples / window_samples.mean()
    x = np.array(
        [
            fft(window_samples * input_data[i : i + num_frame_samples])
            for i in range(0, len(input_data) - num_frame_samples, num_hop_samples)
        ]
    )

    # Normalize the FFT results
    # See "Description and Application of Fourier Transforms and Fourier
    # Series" rev A05 by Matthew Rankin for a description on why the
    # normalization is 2 / N except for the DC component which is 1 / N
    # Only deal with the single-sided FFT, so cut it in half
    x = x[:, : num_frame_samples // 2]
    # Convert from complex to absolute values
    x = np.abs(x)
    # Divide all components by the num_frame_samples
    # Multiply all but the DC component by 2
    non_dc_normalization = 2 / num_frame_samples
    x[:, 1:] = x[:, 1:] * non_dc_normalization
    x[:, 0] = x[:, 0] / num_frame_samples

    # Create the time vector. Each frame advances by num_hop_samples, which is
    # the requested hop truncated to a whole number of samples, so the time
    # step has to be derived from that rather than from hop_size_sec.
    sec_per_hop = num_hop_samples / sampling_frequency_hz
    time_vector_stft = np.arange(x.shape[0]) * sec_per_hop

    # Calculate the width of each frequency bin
    hz_per_freq_bin = sampling_frequency_hz / num_frame_samples

    # Create the frequency vector
    freq_vector_stft = np.arange(x.shape[1]) * hz_per_freq_bin

    return (x, time_vector_stft, freq_vector_stft, hz_per_freq_bin)


def hz2khz(frequency_in_hz: float) -> float:
    """Convert a value from Hz to kHz.

    Args:
        frequency_in_hz: An interger or floating point number containing the
            frequency value in Hz that is to be converted.

    Return:
        The frequency in kHz.

    """
    return frequency_in_hz / 1000


# The windows that smooth() accepts, mapped to the function creating each one.
# A flat window is all ones, which produces a moving average.
WINDOW_FUNCTIONS = {
    "flat": np.ones,
    "hanning": np.hanning,
    "hamming": np.hamming,
    "bartlett": np.bartlett,
    "blackman": np.blackman,
}


def smooth(
    x: npt.NDArray, window_len: int = 11, window: str = "hanning"
) -> npt.NDArray:
    """Smooth the data using a window with requested size.

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

    Raises:
        ValueError: The input is not one dimensional, or the window is not
            one that WINDOW_FUNCTIONS names.
        IndexError: The input is shorter than the window.

    Example:
        Smoothing an alternating signal with a flat window, which is to say
        with a moving average, pulls each sample toward the average of the
        window around it:

        >>> import numpy as np
        >>> x = np.array([1., 5., 1., 5., 1., 5., 1., 5., 1., 5., 1.])
        >>> np.round(smooth(x, window_len=5, window='flat'), 2)
        array([2.6, 3.4, 2.6, 3.4, 2.6, 3.4, 2.6, 3.4, 2.6, 2.6, 2.6])

        The smoothed signal is as long as the signal given:

        >>> smooth(x, window_len=5).size == x.size
        True

        A window shorter than three samples has nothing to average, so the
        signal comes back untouched:

        >>> smooth(x, window_len=2) is x
        True

    See also:
        numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman,
        numpy.convolve, scipy.signal.lfilter

    """
    if x.ndim != 1:
        raise ValueError("Function smooth only accepts 1D arrays.")

    if window_len < 3:
        return x

    # If window_len is not odd, add one so that it is odd. This happens before
    # the length is checked below, since it is the window that gets used that
    # has to fit within the signal.
    if window_len & 1:
        pass
    else:
        window_len += 1

    if x.size < window_len:
        raise IndexError("Input vector needs to be bigger than window size.")

    if window not in WINDOW_FUNCTIONS:
        raise ValueError(
            f"Window must be one of: {', '.join(repr(w) for w in WINDOW_FUNCTIONS)}"
        )

    s = np.r_[x[window_len - 1 : 0 : -1], x, x[-1:-window_len:-1]]

    w = WINDOW_FUNCTIONS[window](window_len)
    y = np.convolve(w / w.sum(), s, mode="valid")
    samples_to_strip = (window_len - 1) // 2
    return y[samples_to_strip : len(y) - samples_to_strip]


def smooth2(x: npt.NDArray, beta: int = 3, window_len: int = 11) -> npt.NDArray:
    """Smooth function using Kaiser window

    Args:
        x: ndarray containing the signal to be smoothed
        beta: beta to use as part of the Kaiser smoothing
        window_len: Integer length of window to be used in Kaiser
            smoothing, which must be odd or it will be made odd.

    Returns:
        An ndarrary containing the smoothed signal.

    Raises:
        ValueError: The input is not one dimensional.
        IndexError: The input is shorter than the window.

    """
    if x.ndim != 1:
        raise ValueError("Function smooth2 only accepts 1D arrays.")

    # If window_len is not odd, add one so that it is odd. This happens before
    # the length is checked below, since it is the window that gets used that
    # has to fit within the signal.
    if window_len & 1:
        pass
    else:
        window_len += 1

    if x.size < window_len:
        raise IndexError("Input vector needs to be bigger than window size.")

    s = np.r_[x[window_len - 1 : 0 : -1], x, x[-1:-window_len:-1]]
    w = np.kaiser(window_len, beta)
    y = np.convolve(w / w.sum(), s, mode="valid")
    samples_to_strip = (window_len - 1) // 2
    return y[samples_to_strip : len(y) - samples_to_strip]


def _check_stft_data(
    stft_data: npt.NDArray,
    time_vector: npt.NDArray | None = None,
    freq_vector: npt.NDArray | None = None,
) -> None:
    """Check STFT data against the vectors describing its axes.

    Args:
        stft_data: The 2D array of shape (time, freq) to check.
        time_vector: An optional vector that has to match the time axis.
        freq_vector: An optional vector that has to match the frequency axis.

    Raises:
        ValueError: The stft_data is not two dimensional.
        IndexError: A vector does not describe the axis it is given for.

    """
    if stft_data.ndim != 2:
        raise ValueError(
            f"The STFT data needs to be 2D with a shape of (time, freq), but "
            f"it has {stft_data.ndim} dimension(s)."
        )
    if time_vector is not None and time_vector.size != stft_data.shape[0]:
        raise IndexError(
            f"The size of the time vector, {time_vector.size}, does not match "
            f"the {stft_data.shape[0]} time bins in the STFT data."
        )
    if freq_vector is not None and freq_vector.size != stft_data.shape[1]:
        raise IndexError(
            f"The size of the freq vector, {freq_vector.size}, does not match "
            f"the {stft_data.shape[1]} frequency bins in the STFT data."
        )


def calculate_peak_hold(
    stft_data: npt.NDArray, frequency_array: npt.NDArray
) -> npt.NDArray:
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
        ValueError: The stft_data is not two dimensional.
        IndexError: The frequency_array does not describe the frequency axis
            of the stft_data.

    """
    _check_stft_data(stft_data, freq_vector=frequency_array)
    data_type = np.dtype([("frequency", "f8"), ("amplitude", "f8")])
    peak_hold = np.zeros(frequency_array.size, dtype=data_type)
    peak_hold["frequency"] = frequency_array
    peak_hold["amplitude"] = np.amax(stft_data, axis=0)
    return peak_hold


def _bin_holding(value: float, vector: npt.NDArray, step: float) -> int:
    """Determine the bin of vector holding value.

    A bin covers half a step either side of its own value, so the bin holding
    a value is the one whose value is nearest to it. That is the rule
    freq_bin() applies to a frequency, applied here to any vector.

    Taking the nearest bin at each end of a range also gives exactly the bins
    that overlap that range. The lowest bin reaching up into the range is the
    nearest one to its start, and the highest bin reaching down into it is the
    nearest one to its stop, so no bin holding any part of the range is left
    out and no bin outside it is drawn in.

    The result is clamped to the vector, so a value outside the range the
    vector covers gives the first or the last bin rather than an out of range
    index.
    """
    return min(max(freq_bin(value, vector[0], step), 0), vector.size - 1)


def plot_spectrogram(
    stft_data: npt.NDArray,
    time_vector: npt.NDArray,
    freq_vector: npt.NDArray,
    plot_axis: matplotlib.axes.Axes,
    freq_plot_range: tuple[int, int] | None = None,
    time_plot_range: tuple[int, int] | None = None,
    plot_title: str | None = None,
    plot_xlabel: str | None = None,
    plot_ylabel: str | None = None,
    colorbar_label: str | None = None,
    colorbar_fontsize: int = 8,
) -> matplotlib.image.AxesImage:
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
        plot_axis: matplotlib axis to which this plot should be added.
        freq_plot_range: An optional tuple containing the start and stop
            frequency in Hz for the spectrogram plot. Both ends are
            inclusive, and every bin holding any part of the range is
            plotted, so a range need not land on a bin.
        time_plot_range: An optional tuple containing the start and stop time
            in seconds for the spectrogram plot. Both ends are inclusive, on
            the same terms as freq_plot_range.
        plot_title: An optional string containing the plot title.
        plot_xlabel: An optional string containing the x-axis label.
        plot_ylabel: An optional string containing the y-axis label.
        colorbar_label: An optional string with the label to be added to the
            colorbar. If excluded then the colorbar is not plotted.
        colorbar_fontsize: An optional integer providing the colorbar font
            size.

    Returns:
        matplotlib handle to the spectrogram

    Raises:
        ValueError: The stft_data is not two dimensional.
        IndexError: A vector does not describe the axis it is given for, or
            holds fewer than the two values needed to give a step size.

    """
    _check_stft_data(stft_data, time_vector=time_vector, freq_vector=freq_vector)
    for name, vector in (("time_vector", time_vector), ("freq_vector", freq_vector)):
        if vector.size < 2:
            raise IndexError(
                f"The {name} needs at least two values to give a step size, "
                f"but it holds {vector.size}."
            )
    if freq_plot_range is None:
        start_freq_plot = freq_vector[0]
        stop_freq_plot = freq_vector[-1]
    else:
        start_freq_plot, stop_freq_plot = freq_plot_range

    if time_plot_range is None:
        start_time_plot = time_vector[0]
        stop_time_plot = time_vector[-1]
    else:
        start_time_plot, stop_time_plot = time_plot_range
    # Calculate the hz_per_freq_bin assuming that the frequency steps are
    # equal.
    hz_per_freq_bin = freq_vector[1] - freq_vector[0]
    sec_per_time_bin = time_vector[1] - time_vector[0]
    # Determine the bins holding the requested start and stop values. Both
    # ends are inclusive, so the stop bins are included in the slices.
    start_freq_bin = _bin_holding(start_freq_plot, freq_vector, hz_per_freq_bin)
    stop_freq_bin = _bin_holding(stop_freq_plot, freq_vector, hz_per_freq_bin)
    start_time_bin = _bin_holding(start_time_plot, time_vector, sec_per_time_bin)
    stop_time_bin = _bin_holding(stop_time_plot, time_vector, sec_per_time_bin)
    # Create the spectrogram
    spectrogram = plot_axis.imshow(
        stft_data[
            start_time_bin : stop_time_bin + 1,
            start_freq_bin : stop_freq_bin + 1,
        ].T,
        origin="lower",
        aspect="auto",
        interpolation="nearest",
    )
    if colorbar_label is not None:
        cb = plt.colorbar(spectrogram, ax=plot_axis)
        cb.ax.tick_params(labelsize=colorbar_fontsize)
        cb.set_label(colorbar_label)
    # imshow puts the extent at the outer edges of the drawn bins, so extend
    # by half a bin on each side to place each bin center at its own time and
    # frequency. Use the bins that were actually drawn rather than the
    # requested values, which need not line up with a bin.
    spectrogram.set_extent(
        (
            time_vector[start_time_bin] - sec_per_time_bin / 2,
            time_vector[stop_time_bin] + sec_per_time_bin / 2,
            freq_vector[start_freq_bin] - hz_per_freq_bin / 2,
            freq_vector[stop_freq_bin] + hz_per_freq_bin / 2,
        )
    )
    if plot_title is not None:
        plot_axis.set_title(plot_title)
    if plot_xlabel is not None:
        plot_axis.set_xlabel(plot_xlabel)
    if plot_ylabel is not None:
        plot_axis.set_ylabel(plot_ylabel)
    return spectrogram


def plot_peak_hold(
    axis: matplotlib.axes.Axes,
    stft_data: npt.NDArray,
    frequency_array: npt.NDArray,
    title: str | None = None,
    xlabel: str | None = None,
    ylabel: str | None = None,
    plot_freq_limits: tuple[int, int] | None = None,
    plot_amp_limits: tuple[float, float] | None = None,
    limit_array: npt.NDArray | None = None,
    trace_label: str | None = None,
):
    """Plot the peak hold for a 2D STFT array

    Args:
        axis: matplotlib axis to which this plot should be added.
        stft_data: A 2D numpy ndarray of shape (time, freq) containing the
            amplitude over both freq and time.
        frequency_array: A 1D numpy ndarray containing the frequencies in
            Hz of the stft_data.
        title: An optional title to be added to the plot.
        xlabel: An optional x-axis label to be added to the plot.
        ylabel: An optional y-axis label to be added to the plot.
        plot_freq_limits: An optional tuple containing the starting and ending
            frequencies to be used in the plot.
        plot_amp_limits: An optional tuple containing the minimum and maximum
            amplitude values.
        limit_array: An optional 1D numpy ndarray containing the limits for the
            plotted data of dtype = [('frequency', 'f8'), ('amplitude', 'f8')]

    The peak hold is drawn onto the given axis. Nothing is returned.

    Raises:
        ValueError: The stft_data is not two dimensional, or the limit_array
            does not carry the fields it is read by.
        IndexError: The frequency_array does not describe the frequency axis
            of the stft_data.

    """
    if limit_array is not None:
        fields = limit_array.dtype.names
        if fields is None or not {"frequency", "amplitude"}.issubset(fields):
            raise ValueError(
                "The limit_array needs a structured dtype carrying "
                "'frequency' and 'amplitude', such as the array that "
                "calculate_peak_hold() returns."
            )
    peak_hold = calculate_peak_hold(stft_data, frequency_array)
    if trace_label is not None:
        axis.loglog(peak_hold["frequency"], peak_hold["amplitude"], label=trace_label)
    else:
        axis.loglog(peak_hold["frequency"], peak_hold["amplitude"])
    if limit_array is not None:
        axis.loglog(limit_array["frequency"], limit_array["amplitude"])
    if plot_freq_limits is not None:
        axis.set_xlim(plot_freq_limits)
    if plot_amp_limits is not None:
        axis.set_ylim(plot_amp_limits)
    if title is not None:
        axis.set_title(title)
    if xlabel is not None:
        axis.set_xlabel(xlabel)
    if ylabel is not None:
        axis.set_ylabel(ylabel)
    axis.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
    axis.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter("%g"))
    axis.grid(visible=True, which="major", color="0.25", linestyle="-")
    axis.grid(visible=True, which="minor", color="0.75", linestyle="-")
    axis.set_axisbelow(True)


def single_frequency_over_time(
    stft_data: npt.NDArray,
    freq_array: npt.NDArray,
    time_array: npt.NDArray,
    frequency: float,
) -> npt.NDArray:
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
        ValueError: The stft_data is not two dimensional.
        IndexError: The size of the STFT does not match the given frequency
            and/or time arrays.

    """
    _check_stft_data(stft_data, time_vector=time_array, freq_vector=freq_array)
    # Create the array to return the time and amplitude
    data_type = np.dtype([("time", "f8"), ("amplitude", "f8")])
    stft_at_frequency = np.zeros(time_array.size, dtype=data_type)
    stft_at_frequency["time"] = time_array
    bin_number = freq_bin(frequency, freq_array[0], freq_array[1] - freq_array[0])
    stft_at_frequency["amplitude"] = stft_data[:, bin_number]
    return stft_at_frequency


def freq_bin(desired_freq: float, first_freq: float, hz_per_freq_bin: float) -> int:
    """Determine the frequency bin for the desired frequency.

    Given the width of each frequency bin (Hz) and the frequency of the first
    bin, determine the bin number for the given desired frequency (Hz).

    Frequencies falling exactly halfway between two bins are rounded up to the
    higher bin. Note that this differs from the builtin round(), which rounds
    a halfway value to the nearest even bin.

    Args:
        desired_freq: The frequency in Hz for which the bin is desired.
        first_freq: The frequency in Hz of the first bin.
        hz_per_freq_bin: The width of each frequency bin in Hz.

    Returns:
        The bin number containing the desired frequency.

    """
    return math.floor((desired_freq - first_freq) / hz_per_freq_bin + 0.5)
