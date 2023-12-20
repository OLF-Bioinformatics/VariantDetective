"""
This module contains a class for describing fragment length distributions
(described by the gamma distribution) and related functions.

Copyright (C) 2024 Phil Charron (phil.charron@inspection.gc.ca)
https://github.com/OLF-Bioinformatics/VariantDetective
Portions Copyright (C) 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Badread
"""

import numpy as np
import sys

### NEEDED
class FragmentLengths(object):
    def __init__(self, mean, stdev, output=sys.stderr):
        self.mean = mean
        self.stdev = stdev
        if self.stdev == 0:
            self.gamma_k, self.gamma_t = None, None
        else:  # gamma distribution
            gamma_a, gamma_b, self.gamma_k, self.gamma_t = gamma_parameters(mean, stdev)

    def get_fragment_length(self):
        if self.stdev == 0:
            return int(round(self.mean))
        else:  # gamma distribution
            fragment_length = int(round(np.random.gamma(self.gamma_k, self.gamma_t)))
            return max(fragment_length, 1)

def gamma_parameters(gamma_mean, gamma_stdev):
    # Shape and rate parametrisation:
    gamma_a = (gamma_mean ** 2) / (gamma_stdev ** 2)
    gamma_b = gamma_mean / (gamma_stdev ** 2)

    # Shape and scale parametrisation:
    gamma_k = (gamma_mean ** 2) / (gamma_stdev ** 2)
    gamma_t = (gamma_stdev ** 2) / gamma_mean

    return gamma_a, gamma_b, gamma_k, gamma_t