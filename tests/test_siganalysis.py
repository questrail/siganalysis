import itertools

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pytest

import siganalysis

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
