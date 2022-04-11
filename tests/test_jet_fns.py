import numpy as np
import unittest
from sscpol.jet_fns import *

# def test_run_ssc_full():
#     """
#     run jet model
#     Input: [W_j, E_max, alpha, theta_open_p, gamma_bulk, B0, theta_obs, A_eq, E_min]
#     """
#     x0 = [10**38.436129, 10**10.12309571, 1.7230574,  39.97781039,
#                15.28831413, 10**-4.39991182, 3.4409881, 1.07342862, 10**5.76709691]
#     sync, ic = run_ssc_full(x0, nblocks=1, seed=42)
#     #Print all values being tested
#     print(np.log2(sync[1,25]),np.log2(ic[1,19]),
#             np.sqrt(sync[2,25]**2 + sync[3,25]**2) / sync[1,25],np.sqrt(ic[2,19]**2 + ic[3,19]**2) / ic[1,19],
#             np.log2(sync[2,25]),np.log2(ic[2,19]), np.log2(-sync[3,25]),np.log2(-ic[3,19]),)
#     assert not np.isnan(sync.sum())
#     assert not np.isnan(ic.sum())
#     #Test flux density
#     np.testing.assert_approx_equal(np.log2(sync[1,25]), 129.9947468025604, significant=3)
#     np.testing.assert_approx_equal(np.log2(ic[1,19]), 126.52255838032376, significant=3)
#     #Test polarization fraction
#     np.testing.assert_approx_equal(np.sqrt(sync[2,25]**2 + sync[3,25]**2) / sync[1,25], 0.6709805756019441, significant=3)
#     np.testing.assert_approx_equal(np.sqrt(ic[2,19]**2 + ic[3,19]**2) / ic[1,19], 0.16587311755562797, significant=3)
#     #Test Stokes parameters
#     np.testing.assert_approx_equal(np.log2(sync[2,25]), 128.99177001456783, significant=3)
#     np.testing.assert_approx_equal(np.log2(ic[2,19]), 123.503390683244, significant=3)
#     np.testing.assert_approx_equal(np.log2(-sync[3,25]), 128.83825737222233, significant=3)
#     np.testing.assert_approx_equal(np.log2(-ic[3,19]), 123.34987804089852, significant=3)

def test_run_ssc():
    """
    run jet model
    Input: [W_j, E_max, alpha, theta_open_p, gamma_bulk, B0, theta_obs, A_eq, E_min]
    """
    x0 = [10**38.436129, 10**10.12309571, 1.7230574,  39.97781039,
               15.28831413, 10**-4.39991182, 3.4409881, 1.07342862, 10**5.76709691]
    sync, ic = run_ssc(x0, nblocks=1, seed=42)
    #Print all values being tested
    print(np.log2(sync[1,25]),np.log2(ic[1,19]),
            np.sqrt(sync[2,25]**2 + sync[3,25]**2) / sync[1,25],np.sqrt(ic[2,19]**2 + ic[3,19]**2) / ic[1,19],
            np.log2(sync[2,25]),np.log2(ic[2,19]), np.log2(-sync[3,25]),np.log2(-ic[3,19]),)
    assert not np.isnan(sync.sum())
    assert not np.isnan(ic.sum())
    #Test flux density
    np.testing.assert_approx_equal(np.log2(sync[1,25]), 129.9947468025604, significant=3)
    np.testing.assert_approx_equal(np.log2(ic[1,19]), 126.52255838032376, significant=3)
    #Test polarization fraction
    np.testing.assert_approx_equal(np.sqrt(sync[2,25]**2 + sync[3,25]**2) / sync[1,25], 0.6709805756019441, significant=3)
    np.testing.assert_approx_equal(np.sqrt(ic[2,19]**2 + ic[3,19]**2) / ic[1,19], 0.16587311755562797, significant=3)
    #Test Stokes parameters
    np.testing.assert_approx_equal(np.log2(sync[2,25]), 128.99177001456783, significant=3)
    np.testing.assert_approx_equal(np.log2(ic[2,19]), 123.503390683244, significant=3)
    np.testing.assert_approx_equal(np.log2(-sync[3,25]), 128.83825737222233, significant=3)
    np.testing.assert_approx_equal(np.log2(-ic[3,19]), 123.34987804089852, significant=3)

# def test_run_ssc_block7():
#     """
#     run jet model
#     Input: [W_j, E_max, alpha, theta_open_p, gamma_bulk, B0, theta_obs, A_eq, E_min]
#     """
#     x0 = [10**38.436129, 10**10.12309571, 1.7230574,  39.97781039,
#                15.28831413, 10**-4.39991182, 3.4409881, 1.07342862, 10**5.76709691]
#     sync, ic = run_ssc(x0, nblocks=7, seed=41)
# #     #Print all values being tested
#     print(np.log2(sync[1,25]),np.log2(ic[1,19]),
#             np.sqrt(sync[2,25]**2 + sync[3,25]**2) / sync[1,25],np.sqrt(ic[2,19]**2 + ic[3,19]**2) / ic[1,19],
#             np.log2(-sync[2,25]),np.log2(-ic[2,19]), np.log2(sync[3,25]),np.log2(ic[3,19]),)
#     assert not np.isnan(sync.sum())
#     assert not np.isnan(ic.sum())
#     np.testing.assert_approx_equal(np.log2(-sync[2,30]), 131.83795208670585, significant=4)
#     np.testing.assert_approx_equal(np.log2(-ic[2,20]), 125.76012219483809, significant=3)
