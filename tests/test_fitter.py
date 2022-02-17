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

def test_log_likelihood():
    """
    run jet model
    Input: [W_j, E_max, alpha, theta_open_p, gamma_bulk, B0, theta_obs, A_eq, E_min]
    """
    bounds = [(36.3, 39),(9.7,13),(1.65,2.05),(15,75),(7,20),(-4.7,-3.1),(0.1,6),(0.8,1.4),(5.72,6.8)]
    fitter = SSC_Fitter("J2011flare")
    loss = fitter.loglike([37.1,11.1,1.85,40.0,16,-4,1.5,1,6.7])
    loss2 = fitter.loglike([37.1,11.1,1.85,40.0,16,-4,1.5,1,6.7])
    assert not np.isnan(loss)
    assert loss != loss2


# def test_cross_entropy():
#     """
#     run jet model
#     Input: [W_j, E_max, alpha, theta_open_p, gamma_bulk, B0, theta_obs, A_eq, E_min]
#     """
#     bounds = [(36.3, 39),(9.7,13),(1.65,2.05),(15,75),(7,20),(-4.7,-3.1),(0.1,6),(0.8,1.4),(5.72,6.8)]
#     fitter = SSC_Fitter("J2011flare")
#     loss, mu, Sigma,  = fitter.cross_entropy_method(bounds, n_elite=4, n_samples=10, max_k=2)
#     print('LOSS: ',loss, mu, Sigma, )