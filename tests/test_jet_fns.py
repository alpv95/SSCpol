import sys
import numpy as np
import unittest
from jet_fns import *

def test_run_ssc():
    """
    run jet model
    """
    x0 = [1.3E37,1.7E11,1.85,40.0,16,1E-4,1.5,1]
    sync, ic = run_ssc(x0, nblocks=7, seed=13)
    assert not np.isnan(sync.sum())
    assert not np.isnan(ic.sum())
    #np.testing.assert_almost_equal(result, 3.6868957, decimal=3)
