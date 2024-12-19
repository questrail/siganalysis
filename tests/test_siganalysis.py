import unittest

import numpy as np

import siganalysis


class TestTimeSlicing(unittest.TestCase):
    def test_small_sample(self):
        sampling_rate = 10
        sample_size = 53
        low_small_zip = [(0, 10), (10, 20), (20, 30), (30, 40), (40, 50), (50, 53)]
        self.assertEqual(
            low_small_zip, siganalysis.time_slice_zip(sample_size, sampling_rate)
        )

    def test_small_sample_multiple_of_sampling_rate(self):
        sampling_rate = 10
        sample_size = 50
        low_small_zip = [(0, 10), (10, 20), (20, 30), (30, 40), (40, 50)]
        self.assertEqual(
            low_small_zip, siganalysis.time_slice_zip(sample_size, sampling_rate)
        )

    def test_large_sample(self):
        sampling_rate = 96000
        sample_size = 960101
        zip = [
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
        self.assertEqual(zip, siganalysis.time_slice_zip(sample_size, sampling_rate))

    def test_large_sample_multiple_of_sampling_rate(self):
        sampling_rate = 96000
        sample_size = 960000
        zip = [
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
        self.assertEqual(zip, siganalysis.time_slice_zip(sample_size, sampling_rate))


class TestShortTimeFourierTransform(unittest.TestCase):
    def setUp(self):
        """Create a time-domain signal comprised of two frequencies."""
        # Setup a 11kHz 1V and 3kHz 0.4V signal in time domain
        # Signals are sampled at 96 kHz (common for LX-10) for 10 seconds
        self.sampling_rate = 96000
        num_samples = 960000
        time_array = (1 / self.sampling_rate) * np.arange(num_samples)
        self.signal_1_frequency = 11000
        self.signal_1_amplitude = 1
        self.signal_2_frequency = 3000
        self.signal_2_amplitude = 0.4
        signal_1 = self.signal_1_amplitude * np.sin(
            2 * np.pi * self.signal_1_frequency * time_array
        )
        signal_2 = self.signal_2_amplitude * np.sin(
            2 * np.pi * self.signal_2_frequency * time_array
        )
        self.test_signal = signal_1 + signal_2

        # Setup variables needed for the STFT
        self.frame_size_sec = 1
        self.hop_size_sec = 0.5
        self.use_hamming_window = True

    def test_stft_shape_size(self):
        data_stft, time_array_stft, freq_array_stft, freq_bin_size = siganalysis.stft(
            self.test_signal,
            self.sampling_rate,
            self.frame_size_sec,
            self.hop_size_sec,
            self.use_hamming_window,
        )
        known_freq_size = int(self.frame_size_sec * self.sampling_rate / 2)
        self.assertEqual(data_stft.shape[1], known_freq_size)
        self.assertEqual(data_stft.shape[1], freq_array_stft.size)
        self.assertEqual(data_stft.shape[0], time_array_stft.size)

    def test_stft_freq_bin_size(self):
        data_stft, time_array_stft, freq_array_stft, freq_bin_size = siganalysis.stft(
            self.test_signal,
            self.sampling_rate,
            self.frame_size_sec,
            self.hop_size_sec,
            self.use_hamming_window,
        )
        self.assertEqual(freq_bin_size, 1 / self.frame_size_sec)

    def test_stft_another_freq_bin_size(self):
        frame_size_sec = 0.5
        correct_freq_per_bin = 1 / frame_size_sec
        data_stft, time_array_stft, freq_array_stft, freq_bin_size = siganalysis.stft(
            self.test_signal,
            self.sampling_rate,
            frame_size_sec,
            self.hop_size_sec,
            self.use_hamming_window,
        )
        self.assertEqual(freq_bin_size, correct_freq_per_bin)


class TestPeakHold(unittest.TestCase):
    def setUp(self):
        """Create a time-domain signal comprised of two frequencies."""
        # Setup a 11kHz 1V signal in time domain
        # Signals are sampled at 96 kHz (common for LX-10) for 10 seconds
        sampling_rate = 96000
        num_samples = 960000
        time_array = (1.0 / sampling_rate) * np.arange(num_samples)
        signal_1_frequency = 11000
        self.signal_1_amplitude = 1
        signal_1 = self.signal_1_amplitude * np.sin(
            2 * np.pi * signal_1_frequency * time_array
        )
        test_signal = signal_1

        # Setup variables needed for the STFT
        frame_size_sec = 1
        hop_size_sec = 0.5
        use_hamming_window = False
        self.data_stft, self.time_array_stft, self.freq_array_stft, freq_bin_size = (
            siganalysis.stft(
                test_signal,
                sampling_rate,
                frame_size_sec,
                hop_size_sec,
                use_hamming_window,
            )
        )

    def test_peak_hold_size(self):
        peak_hold = siganalysis.calculate_peak_hold(
            self.data_stft, self.freq_array_stft
        )
        self.assertEqual(self.data_stft.shape[1], self.freq_array_stft.size)
        self.assertEqual(self.data_stft.shape[0], self.time_array_stft.size)
        self.assertEqual(peak_hold.size, self.freq_array_stft.size)

    def test_peak_hold_amplitude(self):
        peak_hold = siganalysis.calculate_peak_hold(
            self.data_stft, self.freq_array_stft
        )
        self.assertEqual(np.amax(peak_hold["amplitude"]), np.amax(self.data_stft))

    def test_peak_hold_amplitude_calculation(self):
        peak_hold = siganalysis.calculate_peak_hold(
            self.data_stft, self.freq_array_stft
        )
        self.assertAlmostEqual(np.amax(peak_hold["amplitude"]), self.signal_1_amplitude)

    def test_peak_hold_frequency(self):
        peak_hold = siganalysis.calculate_peak_hold(
            self.data_stft, self.freq_array_stft
        )
        self.assertEqual(peak_hold["frequency"][-1], self.freq_array_stft[-1])

    def test_peak_hold_size_error(self):
        self.assertRaises(
            IndexError,
            siganalysis.calculate_peak_hold,
            self.data_stft,
            self.freq_array_stft[:-1],
        )


class TestFrequencyConversion(unittest.TestCase):
    def test_convert_hz_to_khz(self):
        given_frequencies_hz = [1, 20, 300, 4000, 50000, 600000, 7000000]
        expected_khz = [0.001, 0.02, 0.3, 4, 50, 600, 7000]
        calculated_khz = [
            siganalysis.hz2khz(input_frequency_hz)
            for input_frequency_hz in given_frequencies_hz
        ]
        self.assertEqual(calculated_khz, expected_khz)


class TestCalcFreqBin(unittest.TestCase):
    def test_calc_freq_bin_0_hz_starting_freq(self):
        given_freqs_hz = [0, 0.9, 1, 2, 2.9, 3, 4, 398, 399, 400, 401]
        expected_bins = [0, 0, 1, 1, 1, 2, 2, 199, 200, 200, 201]
        calculated_bins = [siganalysis.freq_bin(x, 0, 2) for x in given_freqs_hz]
        self.assertEqual(calculated_bins, expected_bins)

    def test_calc_freq_bin_10_hz_starting_freq(self):
        given_freqs_hz = [10, 12, 14, 399, 400, 401]
        expected_bins = [0, 1, 2, 195, 195, 196]
        calculated_bins = [siganalysis.freq_bin(x, 10, 2) for x in given_freqs_hz]
        self.assertEqual(calculated_bins, expected_bins)
