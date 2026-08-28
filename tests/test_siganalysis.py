import numpy as np
import pytest

import siganalysis

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
        use_hamming_window=False,
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
