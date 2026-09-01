# Copyright (c) 2013-2026 The siganalysis developers. All rights reserved.
# Project site: https://github.com/questrail/siganalysis
# Use of this source code is governed by a MIT-style license that
# can be found in the LICENSE.txt file for the project.
"""Plot the results of the signal analysis routines.

These live apart from the analysis routines so that matplotlib is imported
only by someone who actually plots. matplotlib is an optional dependency,
installed with the plotting extra:

    pip install siganalysis[plotting]

The functions here are also reachable from the package itself, as
siganalysis.plot_spectrogram and siganalysis.plot_peak_hold, which imports
this module on first use.
"""

try:
    import matplotlib.axes
    import matplotlib.image
    import matplotlib.pyplot as plt
    import matplotlib.ticker
except ImportError as error:
    raise ImportError(
        "siganalysis.plotting needs matplotlib, which is not installed. "
        "Install it with the plotting extra: pip install siganalysis[plotting]"
    ) from error

import numpy.typing as npt

from .siganalysis import _check_stft_data, calculate_peak_hold, freq_bin


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


def plot_spectrogram(  # noqa: PLR0913, PLR0917
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
        if vector.size < 2:  # noqa: PLR2004
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


def plot_peak_hold(  # noqa: PLR0913, PLR0917
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
