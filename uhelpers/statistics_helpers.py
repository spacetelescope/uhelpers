"""Helper functions for statistical analyses

Authors
-------

    Johannes Sahlmann

Use
---

"""

import numpy as np
import scipy
try:
    from scipy.stats import betai
except ImportError:
    from scipy.special import betainc as betai


def f_test_probability(N, p1, Chi2_1, p2, Chi2_2):
    """Return F-Test probability that the simpler model is correct.

      e.g. p1 = 5.; //number of PPM parameters
      e.g. p2 = p1 + 7.; // number of PPM + orbital parameters

    :param N: int
        Number of data points
    :param p1: int
        Number of parameters of the simpler model
    :param Chi2_1: float
        chi^2 corresponding to the simpler model
    :param p2: int
        Number of parameters of the model with more parameters
        p2 > p1
    :param Chi2_2: float
        chi^2 corresponding to the model with more parameters
    :return:
        prob: float
        probability

    """

    nu1 = p2 - p1
    nu2 = N - p2  # degrees of freedom

    if (Chi2_1 < Chi2_2):
        raise RuntimeWarning('Solution better with less parameters')

    # F test
    F0 = nu2 / nu1 * (Chi2_1 - Chi2_2) / Chi2_2

    # probability
    prob = betai(0.5 * nu2, 0.5 * nu1, nu2 / (nu2 + F0 * nu1))

    return prob


def sigma_to_fraction(sigma):
    """
    Convert between a sigma value and a fractional probability
    e.g. 1-sigma ~= 67 %

    two-sided

    :param sigma:
    :return:
    """

    return 1. - scipy.stats.norm.cdf(-1. * sigma) * 2


def fraction_to_sigma(fraction):
    """
    Inverse of sigma_to_fraction

    :param fraction:
    :return:
    """
    return scipy.stats.norm.ppf(1. - (1 - fraction) / 2.)


def binomial_fraction_uncertainty(N_events, N_reference):
    """

    :param N_events:
    :param N_reference:
    :return:
    """
    c = 0.683 #    k = 3;    n = 20;
    k = N_events
    n = N_reference
    if n == 0:
        p_lower = 0
        p_upper = 0
        frac = 0
    else:
        p_lower = dist.beta.ppf((1-c)/2.,k+1,n-k+1)
        p_upper = dist.beta.ppf(1-(1-c)/2.,k+1,n-k+1)
        frac = k/np.float(n)

    return np.array([frac,frac-p_lower,p_upper-frac])
