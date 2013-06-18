import os.path
import unittest

import numpy as np

import siganalysis

class TestTimeSlicing(unittest.TestCase):

    def test_small_sample(self):
        sampling_rate = 10
        sample_size = 53
        low_small_zip = [(0,10),(10,20),(20,30),(30,40),(40,50),(50,53)]
        self.assertEqual(low_small_zip,
                        siganalysis.time_slice_zip(sample_size,
                                                   sampling_rate))

    def test_small_sample_multiple_of_sampling_rate(self):
        sampling_rate = 10
        sample_size = 50
        low_small_zip = [(0,10),(10,20),(20,30),(30,40),(40,50)]
        self.assertEqual(low_small_zip,
                        siganalysis.time_slice_zip(sample_size,
                                                   sampling_rate))

    def test_large_sample(self):
        sampling_rate = 96000
        sample_size = 960101
        zip = [(0,96000),(96000,192000),(192000,288000),(288000,384000),
               (384000,480000),(480000,576000),(576000,672000),
               (672000,768000),(768000,864000),(864000,960000),
               (960000,960101)]
        self.assertEqual(zip,
                         siganalysis.time_slice_zip(sample_size,
                                                    sampling_rate))

    def test_large_sample_multiple_of_sampling_rate(self):
        sampling_rate = 96000
        sample_size = 960000
        zip = [(0,96000),(96000,192000),(192000,288000),(288000,384000),
               (384000,480000),(480000,576000),(576000,672000),
               (672000,768000),(768000,864000),(864000,960000)]
        self.assertEqual(zip,
                         siganalysis.time_slice_zip(sample_size,
                                                    sampling_rate))
