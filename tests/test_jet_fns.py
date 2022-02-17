import numpy as np
import unittest
from sscpol.jet_fns import *

def test_run_ssc():
    """
    run jet model
    Input: [W_j, E_max, alpha, theta_open_p, gamma_bulk, B0, theta_obs, A_eq, E_min]
    """
    x0 = [10**38.436129, 10**10.12309571, 1.7230574,  39.97781039,
               15.28831413, 10**-4.39991182, 3.4409881, 1.07342862, 10**5.76709691]
    sync, ic = run_ssc(x0, nblocks=1, seed=42)
    assert not np.isnan(sync.sum())
    assert not np.isnan(ic.sum())
    np.testing.assert_approx_equal(np.log2(sync[2,30]), np.log2(4.141691912633195e+38), significant=4)
    np.testing.assert_approx_equal(np.log2(ic[2,20]), 123.98748246880184, significant=4)


def test_run_ssc_block7():
    """
    run jet model
    Input: [W_j, E_max, alpha, theta_open_p, gamma_bulk, B0, theta_obs, A_eq, E_min]
    """
    x0 = [10**38.436129, 10**10.12309571, 1.7230574,  39.97781039,
               15.28831413, 10**-4.39991182, 3.4409881, 1.07342862, 10**5.76709691]
    sync, ic = run_ssc(x0, nblocks=7, seed=41)
    assert not np.isnan(sync.sum())
    assert not np.isnan(ic.sum())
    np.testing.assert_approx_equal(np.log2(-sync[2,30]), 131.83795208670585, significant=4)
    np.testing.assert_approx_equal(np.log2(-ic[2,20]), 125.76012219483809, significant=3)
