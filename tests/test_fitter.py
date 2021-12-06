import numpy as np
from sscpol.fitter import SSC_Fitter
import math

def test_fitter():
    """
    run jet model
    Input: [W_j, E_max, alpha, theta_open_p, gamma_bulk, B0, theta_obs, A_eq, E_min]
    """
    x0 = [37.1,11.1,1.85,40.0,16,-4,1.5,1,6.7]
    fitter = SSC_Fitter("J2011")
    loss = fitter.loglike(x0)
    print('LOSS: ', loss)
    assert not np.isnan(loss)
    assert not math.isinf(loss)
    assert loss > 0