import numpy as np
import unittest
from sscpol.jet_fns import *

def test_run_ssc():
    """
    run jet model
    Input: [W_j, E_max, alpha, theta_open_p, gamma_bulk, B0, theta_obs, A_eq, E_min]
    """
    x0 = [1.3E37, 1.7E11, 1.85, 40.0, 16, 1E-4, 1.5, 1, 5.11E6]
    sync, ic = run_ssc(x0, nblocks=1, seed=42)
    sync2, ic2 = run_ssc(x0, nblocks=1, seed=42)
    assert not np.isnan(sync.sum())
    assert not np.isnan(ic.sum())
    np.testing.assert_almost_equal(sync2, sync, decimal=6)
    np.testing.assert_almost_equal(ic2, ic, decimal=6)
