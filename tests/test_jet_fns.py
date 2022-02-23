import numpy as np
import unittest
from sscpol.jet_fns import *

def test_run_ssc_full():
    """
    run jet model
    Input: [W_j, E_max, alpha, theta_open_p, gamma_bulk, B0, theta_obs, A_eq, E_min]
    """
    x0 = [10**38.436129, 10**10.12309571, 1.7230574,  39.97781039,
               15.28831413, 10**-4.39991182, 3.4409881, 1.07342862, 10**5.76709691]
    sync, ic = run_ssc_full(x0, nblocks=1, seed=42)
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

# 0.08015695791501892
# 0.7314559253458528
#129.90801369869504 126.28267348553015 0.7313059237185293 0.09700675895222936 129.02924095682144 122.48958287126139 128.87572831447594 122.33607022891591

#parallel one step, full
#120.85059163350853 106.57853522119169 0.6991753194106786 0.09294730066563177 119.90699810267498 102.72377230389276 119.75348546032949 102.57025966154727
# Jet Parameters: 2.72979e+38     1.32769e+10     1.72306e+00     3.99778e+01     1.52883e+01     3.98188e-05     3.44099e+00     1.07343e+00
# Radius at jet base: 8.73321e+14 
# marker_list init = 15Beginning jet analysis... 
# nStep: 1
# x: 2.04636e+11
# B: 3.98110e-05
# R: 8.73492e+14


#parallel one step, fast
#120.74869135451321 106.37440868299436 0.6991753194106785 0.16501099985762077 119.80509782367966 103.34772309263032 119.65158518133417 103.19421045028483
# Jet Parameters: 2.72979e+38     1.32769e+10     1.72306e+00     3.99778e+01     1.52883e+01     3.98188e-05     3.44099e+00     1.07343e+00     5.84921e+05
# Radius at jet base: 8.73321e+14 
# Equipartition fraction U_B / U_e = 1.07343e+00 
# nStep: 1
# x: 2.04682e+11
# B: 3.98110e-05
# R: 8.73492e+14

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
