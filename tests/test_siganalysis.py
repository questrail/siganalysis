import itertools
import subprocess
import sys
import textwrap

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest

import siganalysis
from siganalysis.plotting import _bin_holding

# Draw into a buffer rather than opening a window when testing the plots.
matplotlib.use("Agg")

# Signals are sampled at 96 kHz (common for the LX-10) for 10 seconds.
SAMPLING_RATE_HZ = 96000
NUM_SAMPLES = 960000
SIGNAL_1_FREQUENCY_HZ = 11000
SIGNAL_1_AMPLITUDE = 1
SIGNAL_2_FREQUENCY_HZ = 3000
SIGNAL_2_AMPLITUDE = 0.4

FRAME_SIZE_SEC = 1
HOP_SIZE_SEC = 0.5


def sine_wave(frequency_hz: float, amplitude: float) -> np.ndarray:
    """Create a time domain sine wave at the given frequency and amplitude."""
    time_array = (1 / SAMPLING_RATE_HZ) * np.arange(NUM_SAMPLES)
    return amplitude * np.sin(2 * np.pi * frequency_hz * time_array)


@pytest.fixture(scope="module")
def two_tone_signal():
    """Create a time-domain signal comprised of two frequencies."""
    return sine_wave(SIGNAL_1_FREQUENCY_HZ, SIGNAL_1_AMPLITUDE) + sine_wave(
        SIGNAL_2_FREQUENCY_HZ, SIGNAL_2_AMPLITUDE
    )


@pytest.fixture(scope="module")
def single_tone_stft():
    """Calculate the STFT of a single tone without using a Hamming window.

    Returns the (stft_data, time_array, freq_array, freq_bin_size) tuple
    returned by stft().
    """
    return siganalysis.stft(
        sine_wave(SIGNAL_1_FREQUENCY_HZ, SIGNAL_1_AMPLITUDE),
        SAMPLING_RATE_HZ,
        FRAME_SIZE_SEC,
        HOP_SIZE_SEC,
        window=None,
    )


class TestTimeSlicing:
    def test_small_sample(self):
        expected_zip = [(0, 10), (10, 20), (20, 30), (30, 40), (40, 50), (50, 53)]
        assert siganalysis.time_slice_zip(53, 10) == expected_zip

    def test_small_sample_multiple_of_sampling_rate(self):
        expected_zip = [(0, 10), (10, 20), (20, 30), (30, 40), (40, 50)]
        assert siganalysis.time_slice_zip(50, 10) == expected_zip

    def test_large_sample(self):
        expected_zip = [
            (0, 96000),
            (96000, 192000),
            (192000, 288000),
            (288000, 384000),
            (384000, 480000),
            (480000, 576000),
            (576000, 672000),
            (672000, 768000),
            (768000, 864000),
            (864000, 960000),
            (960000, 960101),
        ]
        assert siganalysis.time_slice_zip(960101, 96000) == expected_zip

    def test_large_sample_multiple_of_sampling_rate(self):
        expected_zip = [
            (0, 96000),
            (96000, 192000),
            (192000, 288000),
            (288000, 384000),
            (384000, 480000),
            (480000, 576000),
            (576000, 672000),
            (672000, 768000),
            (768000, 864000),
            (864000, 960000),
        ]
        assert siganalysis.time_slice_zip(960000, 96000) == expected_zip

    def test_minimum_folds_a_short_last_slice_into_the_one_before(self):
        # 25 samples in slices of 10 leaves a last slice of 5.
        assert siganalysis.time_slice_zip(25, 10) == [(0, 10), (10, 20), (20, 25)]
        assert siganalysis.time_slice_zip(25, 10, 8) == [(0, 10), (10, 25)]

    def test_minimum_leaves_a_long_enough_last_slice_alone(self):
        expected_zip = [(0, 10), (10, 20), (20, 29)]
        assert siganalysis.time_slice_zip(29, 10, 8) == expected_zip

    def test_minimum_leaves_an_exact_multiple_alone(self):
        expected_zip = [(0, 10), (10, 20), (20, 30), (30, 40), (40, 50)]
        assert siganalysis.time_slice_zip(50, 10, 10) == expected_zip

    def test_minimum_covers_every_sample_exactly_once(self):
        # Folding must not drop or repeat a sample, whatever the remainder.
        for number_of_samples in range(1, 60):
            zipped = siganalysis.time_slice_zip(number_of_samples, 10, 7)
            assert zipped[0][0] == 0
            assert zipped[-1][1] == number_of_samples
            for (_, stop), (start, _) in itertools.pairwise(zipped):
                assert stop == start

    def test_minimum_is_met_whenever_there_is_more_than_one_slice(self):
        for number_of_samples in range(1, 60):
            zipped = siganalysis.time_slice_zip(number_of_samples, 10, 7)
            if len(zipped) > 1:
                start, stop = zipped[-1]
                assert stop - start >= 7

    def test_minimum_with_a_single_short_slice(self):
        # There is no earlier slice to fold into, so the minimum cannot be
        # met and the samples are returned as they are rather than dropped.
        assert siganalysis.time_slice_zip(5, 10, 8) == [(0, 5)]

    def test_the_issue_examples(self):
        # The examples given in issue #20.
        assert siganalysis.time_slice_zip(960101, 96000)[-1] == (960000, 960101)
        assert siganalysis.time_slice_zip(960101, 96000, 12000)[-1] == (864000, 960101)

    def test_minimum_longer_than_a_slice(self):
        with pytest.raises(ValueError, match="longer than"):
            siganalysis.time_slice_zip(25, 10, 11)

    def test_no_samples_per_slice(self):
        with pytest.raises(ValueError, match="at least one sample"):
            siganalysis.time_slice_zip(25, 0)


class TestShortTimeFourierTransform:
    def test_stft_shape_size(self, two_tone_signal):
        data_stft, time_array_stft, freq_array_stft, _ = siganalysis.stft(
            two_tone_signal, SAMPLING_RATE_HZ, FRAME_SIZE_SEC, HOP_SIZE_SEC
        )
        known_freq_size = int(FRAME_SIZE_SEC * SAMPLING_RATE_HZ / 2)
        assert data_stft.shape[1] == known_freq_size
        assert data_stft.shape[1] == freq_array_stft.size
        assert data_stft.shape[0] == time_array_stft.size

    def test_stft_freq_bin_size(self, two_tone_signal):
        *_, freq_bin_size = siganalysis.stft(
            two_tone_signal, SAMPLING_RATE_HZ, FRAME_SIZE_SEC, HOP_SIZE_SEC
        )
        assert freq_bin_size == 1 / FRAME_SIZE_SEC

    def test_stft_time_vector_starts_at_zero(self, two_tone_signal):
        _, time_array_stft, *_ = siganalysis.stft(
            two_tone_signal, SAMPLING_RATE_HZ, FRAME_SIZE_SEC, HOP_SIZE_SEC
        )
        assert time_array_stft[0] == 0

    def test_stft_time_vector_steps_by_the_hop_size(self, two_tone_signal):
        data_stft, time_array_stft, *_ = siganalysis.stft(
            two_tone_signal, SAMPLING_RATE_HZ, FRAME_SIZE_SEC, HOP_SIZE_SEC
        )
        assert time_array_stft.size == data_stft.shape[0]
        assert np.diff(time_array_stft) == pytest.approx(HOP_SIZE_SEC)
        assert time_array_stft[-1] == pytest.approx(
            (data_stft.shape[0] - 1) * HOP_SIZE_SEC
        )

    def test_stft_time_vector_uses_the_truncated_hop(self):
        # A hop of 0.0015 sec is 1.5 samples at 1 kHz, and a frame is taken
        # every whole sample, so the frames really advance by 0.001 sec.
        sampling_rate_hz = 1000
        hop_size_sec = 0.0015
        truncated_hop_sec = int(hop_size_sec * sampling_rate_hz) / sampling_rate_hz
        assert truncated_hop_sec != hop_size_sec

        _, time_array_stft, *_ = siganalysis.stft(
            np.zeros(10000), sampling_rate_hz, 0.5, hop_size_sec
        )
        assert np.diff(time_array_stft) == pytest.approx(truncated_hop_sec)

    def test_stft_time_vector_matches_the_frame_start_times(self):
        # Each frame starts at a known sample, so the reported time for a
        # frame has to be that sample divided by the sampling frequency.
        sampling_rate_hz = 1000
        frame_size_sec = 0.1
        hop_size_sec = 0.03
        num_hop_samples = int(hop_size_sec * sampling_rate_hz)

        data_stft, time_array_stft, *_ = siganalysis.stft(
            np.zeros(1000), sampling_rate_hz, frame_size_sec, hop_size_sec
        )
        expected = [
            frame * num_hop_samples / sampling_rate_hz
            for frame in range(data_stft.shape[0])
        ]
        assert list(time_array_stft) == pytest.approx(expected)

    def test_input_shorter_than_one_frame(self):
        with pytest.raises(IndexError, match="longer than"):
            siganalysis.stft(
                np.zeros(100), SAMPLING_RATE_HZ, FRAME_SIZE_SEC, HOP_SIZE_SEC
            )

    def test_input_exactly_one_frame(self):
        # A frame is taken only while the start index is short of the end of
        # the signal, so a signal exactly one frame long yields no frames.
        with pytest.raises(IndexError, match="longer than"):
            siganalysis.stft(
                np.zeros(SAMPLING_RATE_HZ), SAMPLING_RATE_HZ, 1, HOP_SIZE_SEC
            )

    def test_frame_shorter_than_one_sample(self, two_tone_signal):
        with pytest.raises(ValueError, match="frame is shorter than one sample"):
            siganalysis.stft(two_tone_signal, SAMPLING_RATE_HZ, 0, HOP_SIZE_SEC)

    def test_hop_shorter_than_one_sample(self, two_tone_signal):
        with pytest.raises(ValueError, match="hop is shorter than one sample"):
            siganalysis.stft(two_tone_signal, SAMPLING_RATE_HZ, FRAME_SIZE_SEC, 0)

    def test_stft_another_freq_bin_size(self, two_tone_signal):
        frame_size_sec = 0.5
        *_, freq_bin_size = siganalysis.stft(
            two_tone_signal, SAMPLING_RATE_HZ, frame_size_sec, HOP_SIZE_SEC
        )
        assert freq_bin_size == 1 / frame_size_sec

    def test_stft_amplitude_with_no_window(self, two_tone_signal):
        # Both tones fall on the center of a bin for a 1 sec frame, so an
        # unwindowed frame reports each amplitude exactly.
        data_stft, _, freq_array_stft, freq_bin_size = siganalysis.stft(
            two_tone_signal,
            SAMPLING_RATE_HZ,
            FRAME_SIZE_SEC,
            HOP_SIZE_SEC,
            window=None,
        )
        for frequency_hz, amplitude in (
            (SIGNAL_1_FREQUENCY_HZ, SIGNAL_1_AMPLITUDE),
            (SIGNAL_2_FREQUENCY_HZ, SIGNAL_2_AMPLITUDE),
        ):
            bin_number = siganalysis.freq_bin(
                frequency_hz, freq_array_stft[0], freq_bin_size
            )
            assert data_stft[0, bin_number] == pytest.approx(amplitude)

    def test_stft_amplitude_is_corrected_for_the_window_gain(self, two_tone_signal):
        # A window scales every amplitude by its coherent gain, so without a
        # correction a windowed frame reports the wrong amplitude. The mean of
        # a Hamming window is about 0.54, so an uncorrected frame would report
        # about half of each amplitude.
        data_stft, _, freq_array_stft, freq_bin_size = siganalysis.stft(
            two_tone_signal,
            SAMPLING_RATE_HZ,
            FRAME_SIZE_SEC,
            HOP_SIZE_SEC,
            window="hamming",
        )
        for frequency_hz, amplitude in (
            (SIGNAL_1_FREQUENCY_HZ, SIGNAL_1_AMPLITUDE),
            (SIGNAL_2_FREQUENCY_HZ, SIGNAL_2_AMPLITUDE),
        ):
            bin_number = siganalysis.freq_bin(
                frequency_hz, freq_array_stft[0], freq_bin_size
            )
            assert data_stft[0, bin_number] == pytest.approx(amplitude, rel=1e-3)

    def test_stft_amplitude_matches_with_and_without_a_window(self, two_tone_signal):
        # A tone on the center of a bin has no leakage to suppress, so the
        # window may not change the amplitude reported for it.
        windowed, _, freq_array_stft, freq_bin_size = siganalysis.stft(
            two_tone_signal,
            SAMPLING_RATE_HZ,
            FRAME_SIZE_SEC,
            HOP_SIZE_SEC,
            window="hamming",
        )
        unwindowed, *_ = siganalysis.stft(
            two_tone_signal,
            SAMPLING_RATE_HZ,
            FRAME_SIZE_SEC,
            HOP_SIZE_SEC,
            window=None,
        )
        bin_number = siganalysis.freq_bin(
            SIGNAL_1_FREQUENCY_HZ, freq_array_stft[0], freq_bin_size
        )
        assert windowed[0, bin_number] == pytest.approx(
            unwindowed[0, bin_number], rel=1e-3
        )

    @pytest.mark.parametrize("window_name", sorted(siganalysis.STFT_WINDOWS))
    def test_stft_amplitude_is_corrected_for_every_window(
        self, two_tone_signal, window_name
    ):
        # The gain correction is the mean of the window, so it has to give the
        # right amplitude for every window offered, not just for the Hamming
        # window it was written against.
        data_stft, _, freq_array_stft, freq_bin_size = siganalysis.stft(
            two_tone_signal,
            SAMPLING_RATE_HZ,
            FRAME_SIZE_SEC,
            HOP_SIZE_SEC,
            window=window_name,
        )
        for frequency_hz, amplitude in (
            (SIGNAL_1_FREQUENCY_HZ, SIGNAL_1_AMPLITUDE),
            (SIGNAL_2_FREQUENCY_HZ, SIGNAL_2_AMPLITUDE),
        ):
            bin_number = siganalysis.freq_bin(
                frequency_hz, freq_array_stft[0], freq_bin_size
            )
            assert data_stft[0, bin_number] == pytest.approx(amplitude, rel=1e-3)

    def test_stft_defaults_to_the_hamming_window(self, two_tone_signal):
        without_argument, *_ = siganalysis.stft(
            two_tone_signal, SAMPLING_RATE_HZ, FRAME_SIZE_SEC, HOP_SIZE_SEC
        )
        hamming, *_ = siganalysis.stft(
            two_tone_signal,
            SAMPLING_RATE_HZ,
            FRAME_SIZE_SEC,
            HOP_SIZE_SEC,
            window="hamming",
        )
        assert np.array_equal(without_argument, hamming)

    def test_stft_hanning_is_an_alias_for_hann(self, two_tone_signal):
        hann, *_ = siganalysis.stft(
            two_tone_signal, SAMPLING_RATE_HZ, FRAME_SIZE_SEC, HOP_SIZE_SEC, "hann"
        )
        hanning, *_ = siganalysis.stft(
            two_tone_signal, SAMPLING_RATE_HZ, FRAME_SIZE_SEC, HOP_SIZE_SEC, "hanning"
        )
        assert np.array_equal(hann, hanning)

    def test_stft_windows_differ_from_one_another(self, two_tone_signal):
        # Guard against a registry entry silently pointing at the wrong
        # window, which the amplitude tests above would not catch.
        results = {
            name: siganalysis.stft(
                two_tone_signal,
                SAMPLING_RATE_HZ,
                FRAME_SIZE_SEC,
                HOP_SIZE_SEC,
                window=name,
            )[0]
            for name in ("hamming", "hann", "blackman", "blackmanharris", "flattop")
        }
        for name, data_stft in results.items():
            for other_name, other_stft in results.items():
                if name < other_name:
                    assert not np.array_equal(data_stft, other_stft), (
                        f"{name} and {other_name} gave identical results"
                    )

    def test_stft_accepts_integer_input(self, two_tone_signal):
        # A recorder gives integer samples, so integers have to work. This is
        # why the arguments are not checked for a floating point dtype.
        as_integers = (1000 * two_tone_signal).astype(np.int16)
        data_stft, *_ = siganalysis.stft(
            as_integers, SAMPLING_RATE_HZ, FRAME_SIZE_SEC, HOP_SIZE_SEC
        )
        assert data_stft.dtype == np.float64
        assert data_stft.max() > 0

    def test_stft_unknown_window(self, two_tone_signal):
        with pytest.raises(ValueError, match="Window must be None or one of"):
            siganalysis.stft(
                two_tone_signal,
                SAMPLING_RATE_HZ,
                FRAME_SIZE_SEC,
                HOP_SIZE_SEC,
                window="kaiser",
            )


class TestPeakHold:
    def test_peak_hold_size(self, single_tone_stft):
        data_stft, time_array_stft, freq_array_stft, _ = single_tone_stft
        peak_hold = siganalysis.calculate_peak_hold(data_stft, freq_array_stft)
        assert data_stft.shape[1] == freq_array_stft.size
        assert data_stft.shape[0] == time_array_stft.size
        assert peak_hold.size == freq_array_stft.size

    def test_peak_hold_amplitude(self, single_tone_stft):
        data_stft, _, freq_array_stft, _ = single_tone_stft
        peak_hold = siganalysis.calculate_peak_hold(data_stft, freq_array_stft)
        assert np.amax(peak_hold["amplitude"]) == np.amax(data_stft)

    def test_peak_hold_amplitude_calculation(self, single_tone_stft):
        data_stft, _, freq_array_stft, _ = single_tone_stft
        peak_hold = siganalysis.calculate_peak_hold(data_stft, freq_array_stft)
        assert np.amax(peak_hold["amplitude"]) == pytest.approx(SIGNAL_1_AMPLITUDE)

    def test_peak_hold_frequency(self, single_tone_stft):
        data_stft, _, freq_array_stft, _ = single_tone_stft
        peak_hold = siganalysis.calculate_peak_hold(data_stft, freq_array_stft)
        assert peak_hold["frequency"][-1] == freq_array_stft[-1]

    def test_peak_hold_size_error(self, single_tone_stft):
        data_stft, _, freq_array_stft, _ = single_tone_stft
        with pytest.raises(IndexError):
            siganalysis.calculate_peak_hold(data_stft, freq_array_stft[:-1])

    def test_peak_hold_one_dimensional_stft_data(self, single_tone_stft):
        _, _, freq_array_stft, _ = single_tone_stft
        with pytest.raises(ValueError, match="needs to be 2D"):
            siganalysis.calculate_peak_hold(np.arange(3, dtype=float), freq_array_stft)


class TestFrequencyConversion:
    @pytest.mark.parametrize(
        ("given_freq_hz", "expected_khz"),
        [
            (1, 0.001),
            (20, 0.02),
            (300, 0.3),
            (4000, 4),
            (50000, 50),
            (600000, 600),
            (7000000, 7000),
        ],
    )
    def test_convert_hz_to_khz(self, given_freq_hz, expected_khz):
        assert siganalysis.hz2khz(given_freq_hz) == expected_khz


class TestCalcFreqBin:
    @pytest.mark.parametrize(
        ("given_freq_hz", "expected_bin"),
        [
            (0, 0),
            (0.9, 0),
            (1, 1),
            (2, 1),
            (2.9, 1),
            (3, 2),
            (4, 2),
            (398, 199),
            (399, 200),
            (400, 200),
            (401, 201),
        ],
    )
    def test_calc_freq_bin_0_hz_starting_freq(self, given_freq_hz, expected_bin):
        assert siganalysis.freq_bin(given_freq_hz, 0, 2) == expected_bin

    @pytest.mark.parametrize(
        ("given_freq_hz", "expected_bin"),
        [(10, 0), (12, 1), (14, 2), (399, 195), (400, 195), (401, 196)],
    )
    def test_calc_freq_bin_10_hz_starting_freq(self, given_freq_hz, expected_bin):
        assert siganalysis.freq_bin(given_freq_hz, 10, 2) == expected_bin


class TestSingleFrequencyOverTime:
    @pytest.fixture
    def stft_data(self):
        """Three time steps by ten frequency bins.

        The amplitude in each bin is the bin number so that the selected bin
        is easy to identify.
        """
        return np.tile(np.arange(10, dtype=float), (3, 1))

    @pytest.fixture
    def time_array(self):
        return np.arange(3, dtype=float)

    def test_freq_array_starting_at_0_hz(self, stft_data, time_array):
        freq_array = np.arange(10) * 2.0
        amplitude_at_freq = siganalysis.single_frequency_over_time(
            stft_data, freq_array, time_array, 4
        )
        assert list(amplitude_at_freq["amplitude"]) == [2, 2, 2]

    def test_freq_array_not_starting_at_0_hz(self, stft_data, time_array):
        freq_array = np.arange(10) * 2.0 + 100
        amplitude_at_freq = siganalysis.single_frequency_over_time(
            stft_data, freq_array, time_array, 104
        )
        assert list(amplitude_at_freq["amplitude"]) == [2, 2, 2]

    def test_time_array_is_returned(self, stft_data, time_array):
        freq_array = np.arange(10) * 2.0
        amplitude_at_freq = siganalysis.single_frequency_over_time(
            stft_data, freq_array, time_array, 4
        )
        assert list(amplitude_at_freq["time"]) == list(time_array)

    def test_freq_array_size_error(self, stft_data, time_array):
        with pytest.raises(IndexError):
            siganalysis.single_frequency_over_time(
                stft_data, np.arange(9) * 2.0, time_array, 4
            )

    def test_one_dimensional_stft_data(self, time_array):
        with pytest.raises(ValueError, match="needs to be 2D"):
            siganalysis.single_frequency_over_time(
                np.arange(3, dtype=float), np.arange(3, dtype=float), time_array, 1
            )

    def test_time_array_size_error(self, stft_data):
        with pytest.raises(IndexError):
            siganalysis.single_frequency_over_time(
                stft_data, np.arange(10) * 2.0, np.arange(2, dtype=float), 4
            )


class TestPlotSpectrogram:
    @pytest.fixture
    def spectrogram_data(self):
        """Ten time bins by ten frequency bins.

        Each value is 10 * time bin + frequency bin, so every cell identifies
        which bin it came from.
        """
        return np.arange(100, dtype=float).reshape(10, 10)

    @pytest.fixture
    def time_vector(self):
        return np.arange(10) * 1.0

    @pytest.fixture
    def freq_vector(self):
        return np.arange(10) * 2.0

    @pytest.fixture
    def axis(self):
        figure, axis = plt.subplots()
        yield axis
        plt.close(figure)

    def test_plots_every_bin_by_default(
        self, spectrogram_data, time_vector, freq_vector, axis
    ):
        spectrogram = siganalysis.plot_spectrogram(
            spectrogram_data, time_vector, freq_vector, axis
        )
        assert spectrogram.get_array().shape == (freq_vector.size, time_vector.size)

    def test_plot_ranges_include_both_ends(
        self, spectrogram_data, time_vector, freq_vector, axis
    ):
        # Time bins 2 through 5 and frequency bins 2 through 5, which is four
        # of each because both ends are inclusive.
        spectrogram = siganalysis.plot_spectrogram(
            spectrogram_data,
            time_vector,
            freq_vector,
            axis,
            time_plot_range=(2, 5),
            freq_plot_range=(4, 10),
        )
        plotted = spectrogram.get_array()
        assert plotted.shape == (4, 4)
        # The array is transposed for plotting, so it is (freq, time).
        assert plotted[0][0] == 10 * 2 + 2
        assert plotted[-1][-1] == 10 * 5 + 5

    def test_plot_ranges_are_clamped_to_the_data(
        self, spectrogram_data, time_vector, freq_vector, axis
    ):
        spectrogram = siganalysis.plot_spectrogram(
            spectrogram_data,
            time_vector,
            freq_vector,
            axis,
            time_plot_range=(-5, 99),
            freq_plot_range=(-5, 99),
        )
        assert spectrogram.get_array().shape == (freq_vector.size, time_vector.size)

    def test_extent_puts_each_bin_center_on_its_own_value(
        self, spectrogram_data, time_vector, freq_vector, axis
    ):
        spectrogram = siganalysis.plot_spectrogram(
            spectrogram_data,
            time_vector,
            freq_vector,
            axis,
            time_plot_range=(2, 5),
            freq_plot_range=(4, 10),
        )
        left, right, bottom, top = spectrogram.get_extent()
        num_freq_bins, num_time_bins = spectrogram.get_array().shape

        time_centers = [
            left + (bin + 0.5) * (right - left) / num_time_bins
            for bin in range(num_time_bins)
        ]
        freq_centers = [
            bottom + (bin + 0.5) * (top - bottom) / num_freq_bins
            for bin in range(num_freq_bins)
        ]
        assert time_centers == pytest.approx(list(time_vector[2:6]))
        assert freq_centers == pytest.approx(list(freq_vector[2:6]))

    def test_mismatched_time_vector(
        self, spectrogram_data, time_vector, freq_vector, axis
    ):
        # This used to plot happily, mislabelling the axes with a time vector
        # that did not describe the data.
        with pytest.raises(IndexError, match="time vector"):
            siganalysis.plot_spectrogram(
                spectrogram_data, time_vector[:-1], freq_vector, axis
            )

    def test_mismatched_freq_vector(
        self, spectrogram_data, time_vector, freq_vector, axis
    ):
        with pytest.raises(IndexError, match="freq vector"):
            siganalysis.plot_spectrogram(
                spectrogram_data, time_vector, freq_vector[:-1], axis
            )

    def test_one_dimensional_stft_data(self, time_vector, freq_vector, axis):
        with pytest.raises(ValueError, match="needs to be 2D"):
            siganalysis.plot_spectrogram(
                np.arange(10, dtype=float), time_vector, freq_vector, axis
            )

    def test_a_vector_with_a_single_value(self, axis):
        # The step size comes from the first two values, so one is not enough.
        with pytest.raises(IndexError, match="at least two values"):
            siganalysis.plot_spectrogram(
                np.ones((1, 1)), np.zeros(1), np.zeros(1), axis
            )


class TestBinHolding:
    """The bin holding a value, which is issue #11.

    A bin covers half a step either side of its own value, so the bin holding
    a value is the one whose value is nearest to it. freq_bin() applies that
    rule to a frequency and _bin_holding() applies it to any vector; the two
    have to agree.
    """

    @pytest.fixture
    def freq_vector(self):
        # Ten bins, 10 Hz wide, centered on 0, 10, ... 90 Hz.
        return np.arange(10) * 10.0

    def test_the_issue_example(self, freq_vector):
        # For 10 Hz bins, 12 Hz falls in the second bin, the one covering
        # 10 Hz +/- 5 Hz, rather than the first.
        assert _bin_holding(12, freq_vector, 10.0) == 1

    def test_a_value_in_the_upper_half_of_a_bin(self, freq_vector):
        # 6 Hz is nearer to the 10 Hz bin than to the 0 Hz bin. Truncating,
        # as this used to, put it in the 0 Hz bin.
        assert _bin_holding(6, freq_vector, 10.0) == 1

    def test_agrees_with_freq_bin(self, freq_vector):
        for value in np.arange(-20, 120, 0.5):
            expected = siganalysis.freq_bin(value, freq_vector[0], 10.0)
            expected = min(max(expected, 0), freq_vector.size - 1)
            assert _bin_holding(value, freq_vector, 10.0) == expected

    def test_covers_every_bin_overlapping_a_range(self, freq_vector):
        # Taking the nearest bin at each end of a range gives exactly the
        # bins overlapping that range. Checked against a brute force search
        # over which bins reach into the range, away from the exact bin edges
        # where a value touches two bins and freq_bin() rounds up by rule.
        step = 10.0
        half = step / 2
        for start in np.arange(0, 90, 2.5):
            for stop in np.arange(start, 90, 2.5):
                on_an_edge = any(
                    abs(((value - freq_vector[0]) / step + 0.5) % 1.0) < 1e-9
                    for value in (start, stop)
                )
                if on_an_edge:
                    continue
                overlapping = [
                    index
                    for index, center in enumerate(freq_vector)
                    if center + half >= start and center - half <= stop
                ]
                assert _bin_holding(start, freq_vector, step) == overlapping[0]
                assert _bin_holding(stop, freq_vector, step) == overlapping[-1]

    def test_clamped_to_the_vector(self, freq_vector):
        assert _bin_holding(-1000, freq_vector, 10.0) == 0
        assert _bin_holding(1000, freq_vector, 10.0) == freq_vector.size - 1


class TestSmooth:
    @pytest.fixture
    def alternating(self):
        return np.array([1.0, 5.0] * 10 + [1.0])

    @pytest.mark.parametrize("window", sorted(siganalysis.WINDOW_FUNCTIONS))
    def test_output_is_as_long_as_the_input(self, alternating, window):
        assert siganalysis.smooth(alternating, 5, window).size == alternating.size

    @pytest.mark.parametrize("window_len", range(3, 16))
    def test_output_length_for_every_window_length(self, window_len):
        x = np.arange(40, dtype=float)
        assert siganalysis.smooth(x, window_len).size == x.size

    def test_flat_window_is_a_moving_average(self):
        x = np.array([0.0, 0.0, 3.0, 0.0, 0.0])
        # The middle sample averages with its two neighbors either side.
        assert siganalysis.smooth(x, 3, "flat")[2] == pytest.approx(1.0)

    def test_smoothing_reduces_variation(self, alternating):
        smoothed = siganalysis.smooth(alternating, 5, "hanning")
        assert smoothed.std() < alternating.std()

    def test_an_even_window_is_made_odd(self):
        x = np.arange(40, dtype=float)
        assert np.array_equal(
            siganalysis.smooth(x, 4, "flat"), siganalysis.smooth(x, 5, "flat")
        )

    def test_a_window_shorter_than_three_returns_the_input(self, alternating):
        assert siganalysis.smooth(alternating, 2) is alternating

    def test_two_dimensional_input(self):
        with pytest.raises(ValueError, match="only accepts 1D arrays"):
            siganalysis.smooth(np.zeros((4, 4)))

    def test_input_shorter_than_the_window(self):
        with pytest.raises(IndexError, match="bigger than window size"):
            siganalysis.smooth(np.arange(4, dtype=float), 11)

    def test_input_shorter_than_the_window_once_made_odd(self):
        # A window of 4 becomes 5, which does not fit in 4 samples.
        with pytest.raises(IndexError, match="bigger than window size"):
            siganalysis.smooth(np.arange(4, dtype=float), 4)

    def test_unknown_window(self, alternating):
        with pytest.raises(ValueError, match="Window must be one of"):
            siganalysis.smooth(alternating, 5, "kaiser")


class TestSmooth2:
    @pytest.fixture
    def alternating(self):
        return np.array([1.0, 5.0] * 10 + [1.0])

    @pytest.mark.parametrize("window_len", range(3, 16))
    def test_output_is_as_long_as_the_input(self, window_len):
        x = np.arange(40, dtype=float)
        assert siganalysis.smooth2(x, 3, window_len).size == x.size

    def test_smoothing_reduces_variation(self, alternating):
        assert siganalysis.smooth2(alternating, 3, 5).std() < alternating.std()

    def test_an_even_window_is_made_odd(self):
        x = np.arange(40, dtype=float)
        assert np.array_equal(
            siganalysis.smooth2(x, 3, 4), siganalysis.smooth2(x, 3, 5)
        )

    def test_beta_changes_the_result(self, alternating):
        assert not np.array_equal(
            siganalysis.smooth2(alternating, 1, 5),
            siganalysis.smooth2(alternating, 10, 5),
        )

    def test_two_dimensional_input(self):
        with pytest.raises(ValueError, match="only accepts 1D arrays"):
            siganalysis.smooth2(np.zeros((4, 4)))

    def test_input_shorter_than_the_window(self):
        # This used to return an empty array rather than saying anything.
        with pytest.raises(IndexError, match="bigger than window size"):
            siganalysis.smooth2(np.arange(3, dtype=float), 3, 11)

    def test_input_shorter_than_the_window_once_made_odd(self):
        with pytest.raises(IndexError, match="bigger than window size"):
            siganalysis.smooth2(np.arange(4, dtype=float), 3, 4)


class TestPlotPeakHold:
    @pytest.fixture
    def axis(self):
        figure, axis = plt.subplots()
        yield axis
        plt.close(figure)

    @pytest.fixture
    def stft_data(self):
        return np.array([[1.0, 2.0, 3.0], [4.0, 1.0, 1.0]])

    @pytest.fixture
    def freq_array(self):
        return np.array([10.0, 20.0, 30.0])

    def test_plots_the_peak_of_each_frequency(self, axis, stft_data, freq_array):
        siganalysis.plot_peak_hold(axis, stft_data, freq_array)
        assert len(axis.lines) == 1
        assert list(axis.lines[0].get_xdata()) == list(freq_array)
        assert list(axis.lines[0].get_ydata()) == [4.0, 2.0, 3.0]

    def test_returns_nothing(self, axis, stft_data, freq_array):
        assert siganalysis.plot_peak_hold(axis, stft_data, freq_array) is None

    def test_both_axes_are_logarithmic(self, axis, stft_data, freq_array):
        siganalysis.plot_peak_hold(axis, stft_data, freq_array)
        assert axis.get_xscale() == "log"
        assert axis.get_yscale() == "log"

    def test_a_limit_array_is_drawn_as_a_second_line(self, axis, stft_data, freq_array):
        limit_array = siganalysis.calculate_peak_hold(stft_data, freq_array)
        limit_array["amplitude"] = 10.0
        siganalysis.plot_peak_hold(axis, stft_data, freq_array, limit_array=limit_array)
        assert len(axis.lines) == 2
        assert list(axis.lines[1].get_ydata()) == [10.0, 10.0, 10.0]

    def test_the_trace_label_is_applied(self, axis, stft_data, freq_array):
        siganalysis.plot_peak_hold(axis, stft_data, freq_array, trace_label="peak hold")
        assert axis.lines[0].get_label() == "peak hold"

    def test_the_titles_are_applied(self, axis, stft_data, freq_array):
        siganalysis.plot_peak_hold(
            axis,
            stft_data,
            freq_array,
            title="A title",
            xlabel="An x label",
            ylabel="A y label",
        )
        assert axis.get_title() == "A title"
        assert axis.get_xlabel() == "An x label"
        assert axis.get_ylabel() == "A y label"

    def test_the_limits_are_applied(self, axis, stft_data, freq_array):
        siganalysis.plot_peak_hold(
            axis,
            stft_data,
            freq_array,
            plot_freq_limits=(10, 30),
            plot_amp_limits=(0.5, 5),
        )
        assert axis.get_xlim() == (10, 30)
        assert axis.get_ylim() == (0.5, 5)

    def test_mismatched_frequency_array(self, axis, stft_data, freq_array):
        with pytest.raises(IndexError):
            siganalysis.plot_peak_hold(axis, stft_data, freq_array[:-1])

    def test_a_limit_array_without_the_fields(self, axis, stft_data, freq_array):
        # A plain ndarray is indexed by field name further down, which used to
        # fail with an error about how numpy indexing works.
        with pytest.raises(ValueError, match="structured dtype"):
            siganalysis.plot_peak_hold(
                axis, stft_data, freq_array, limit_array=np.arange(3, dtype=float)
            )

    def test_one_dimensional_stft_data(self, axis, freq_array):
        with pytest.raises(ValueError, match="needs to be 2D"):
            siganalysis.plot_peak_hold(axis, np.arange(3, dtype=float), freq_array)


class TestPackageLayout:
    """The plotting lives apart from the analysis, which is issue #5.

    matplotlib is an optional dependency, so importing the package must not
    import it, while the plotting functions stay reachable from the package.
    """

    def _run(self, code: str) -> str:
        # A subprocess, because this test session has already imported
        # matplotlib and sys.modules cannot be un-imported reliably.
        result = subprocess.run(
            [sys.executable, "-c", textwrap.dedent(code)],
            capture_output=True,
            text=True,
            check=True,
        )
        return result.stdout.strip()

    def test_importing_the_package_does_not_import_matplotlib(self):
        assert (
            self._run("""
            import sys
            import siganalysis
            print("matplotlib" in sys.modules)
            """)
            == "False"
        )

    def test_touching_a_plotting_name_imports_matplotlib(self):
        assert (
            self._run("""
            import sys
            import siganalysis
            siganalysis.plot_spectrogram
            print("matplotlib" in sys.modules)
            """)
            == "True"
        )

    def test_the_analysis_runs_without_matplotlib_installed(self):
        # Block matplotlib outright, standing in for it not being installed.
        assert (
            self._run("""
            import sys
            from importlib.abc import MetaPathFinder

            class Block(MetaPathFinder):
                def find_spec(self, fullname, path=None, target=None):
                    if fullname.split(".")[0] == "matplotlib":
                        raise ModuleNotFoundError(fullname)
                    return None

            sys.meta_path.insert(0, Block())
            import numpy as np
            import siganalysis
            data, *_ = siganalysis.stft(np.zeros(1000), 1000, 0.1, 0.1)
            print(data.shape[0])
            """)
            == "9"
        )

    def test_plotting_without_matplotlib_says_what_to_install(self):
        assert (
            self._run("""
            import sys
            from importlib.abc import MetaPathFinder

            class Block(MetaPathFinder):
                def find_spec(self, fullname, path=None, target=None):
                    if fullname.split(".")[0] == "matplotlib":
                        raise ModuleNotFoundError(fullname)
                    return None

            sys.meta_path.insert(0, Block())
            import siganalysis
            try:
                siganalysis.plot_peak_hold
            except ImportError as error:
                print("siganalysis[plotting]" in str(error))
            """)
            == "True"
        )

    @pytest.mark.parametrize("name", ["plot_spectrogram", "plot_peak_hold"])
    def test_the_plotting_names_resolve_from_the_package(self, name):
        assert callable(getattr(siganalysis, name))

    @pytest.mark.parametrize("name", ["plot_spectrogram", "plot_peak_hold"])
    def test_the_plotting_names_are_listed(self, name):
        # __getattr__ would otherwise hide them from dir() and tab completion.
        assert name in dir(siganalysis)

    def test_an_unknown_name_still_raises_attribute_error(self):
        missing = "not_a_real_function"
        with pytest.raises(AttributeError, match="has no attribute"):
            getattr(siganalysis, missing)

    def test_everything_in_all_is_reachable(self):
        for name in siganalysis.__all__:
            assert getattr(siganalysis, name) is not None
