import numpy as np
from sscpol.jet_fns import *


def test_run_ssc():
    """
    run jet model
    Input: [W_j, E_max, alpha, theta_open_p, gamma_bulk, B0, theta_obs, A_eq, E_min]
    """
    x0 = [10**38.436129, 10**10.12309571, 1.7230574,  39.97781039,
          15.28831413, 10**-4.39991182, 3.4409881, 1.07342862, 10**5.76709691]
    sync, ic = run_ssc(x0, nblocks=1, seed=42)
    # Print all values being tested
    print(np.log2(sync[1, 25]), np.log2(ic[1, 19]),
          np.sqrt(sync[2, 25]**2 + sync[3, 25]**2) / sync[1,
                                                          25], np.sqrt(ic[2, 19]**2 + ic[3, 19]**2) / ic[1, 19],
          np.log2(sync[2, 25]), np.log2(ic[2, 19]), np.log2(-sync[3, 25]), np.log2(-ic[3, 19]),)
    assert not np.isnan(sync.sum())
    assert not np.isnan(ic.sum())
    # Test flux density
    np.testing.assert_approx_equal(
        np.log2(sync[1, 25]), 129.9947468025604, significant=3)
    np.testing.assert_approx_equal(
        np.log2(ic[1, 19]), 126.52255838032376, significant=3)
    # Test polarization fraction
    np.testing.assert_approx_equal(np.sqrt(
        sync[2, 25]**2 + sync[3, 25]**2) / sync[1, 25], 0.731455925, significant=3)  # 0.67
    np.testing.assert_approx_equal(np.sqrt(
        ic[2, 19]**2 + ic[3, 19]**2) / ic[1, 19], 0.09583, significant=3)  # 0.165
