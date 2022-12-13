"""
This module contains  class for describing fragment length distributions
(described by the gamma distribution) and related functions.

Portions Copyright (C) 2022 Phil Charron (phil.charron@inspection.gc.ca)
https://github.com/philcharron-cfia/VariantDetective
Portions Copyright (C) 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Badread
"""

import numpy as np
import scipy.special
import scipy.stats
import sys
from .simulate_tools import float_to_str, print_in_two_columns

### NEEDED
class FragmentLengths(object):
    def __init__(self, mean, stdev, output=sys.stderr):
        self.mean = mean
        self.stdev = stdev
        #print('', file=output)
        if self.stdev == 0:
            self.gamma_k, self.gamma_t = None, None
            #print(f'Using a constant fragment length of {mean} bp', file=output)
        else:  # gamma distribution
            #print('Generating fragment lengths from a gamma distribution:', file=output)
            gamma_a, gamma_b, self.gamma_k, self.gamma_t = gamma_parameters(mean, stdev)
            n50 = int(round(find_n_value(gamma_a, gamma_b, 50)))
            #print_in_two_columns(f'  mean  = {float_to_str(mean):>6} bp',
            #                     f'  stdev = {float_to_str(stdev):>6} bp',
            #                     f'  N50   = {n50:>6} bp',
            #                     'parameters:',
            #                     f'  k (shape)     = {self.gamma_k:.4e}',
            #                     f'  theta (scale) = {self.gamma_t:.4e}',
            #                     output=output)

    def get_fragment_length(self):
        if self.stdev == 0:
            return int(round(self.mean))
        else:  # gamma distribution
            fragment_length = int(round(np.random.gamma(self.gamma_k, self.gamma_t)))
            return max(fragment_length, 1)

def base_distribution_integral(a, b, x):
    # Uses logs to avoid overflow:
    integral = 1.0 - np.exp(inc_gamma_ln(a+1, b*x) - scipy.special.gammaln(a+1))
    return integral

def find_n_value(a, b, n):
    """
    Uses the integral of the base distribution function to binary search for the N50 (or N-whatever)
    value.
    """
    target = 1.0 - (n / 100.0)
    bottom_range = 0.0
    top_range = 1.0
    while base_distribution_integral(a, b, top_range) < target:
        bottom_range = top_range
        top_range *= 2
    guess = (bottom_range + top_range) / 2.0
    while True:
        integral = base_distribution_integral(a, b, guess)
        if top_range - bottom_range < 0.01:
            return guess
        elif integral < target:
            bottom_range = guess
            guess = (bottom_range + top_range) / 2.0
        else:  # integral > target:
            top_range = guess
            guess = (bottom_range + top_range) / 2.0

def gamma_parameters(gamma_mean, gamma_stdev):
    # Shape and rate parametrisation:
    gamma_a = (gamma_mean ** 2) / (gamma_stdev ** 2)
    gamma_b = gamma_mean / (gamma_stdev ** 2)

    # Shape and scale parametrisation:
    gamma_k = (gamma_mean ** 2) / (gamma_stdev ** 2)
    gamma_t = (gamma_stdev ** 2) / gamma_mean

    return gamma_a, gamma_b, gamma_k, gamma_t

def inc_gamma_ln(a, b):
    """
    Natural log of the inc_gamma function.
    """
    return scipy.special.gammaln(a) + np.log(1-scipy.stats.gamma.cdf(b, a))