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
    loss, _, _ = fitter.loglike(x0)
    print('LOSS: ', loss)
    assert not np.isnan(loss)
    assert not math.isinf(loss)
    assert loss > 0

def test_log_likelihood():
    """
    run jet model
    Input: [W_j, E_max, alpha, theta_open_p, gamma_bulk, B0, theta_obs, A_eq, E_min]
    """
    fitter = SSC_Fitter("J2011flare")
    loss, _, _ = fitter.loglike([37.1,11.1,1.85,40.0,16,-4,1.5,1,6.7])
    loss2, _, _ = fitter.loglike([37.1,11.1,1.85,40.0,16,-4,1.5,1,6.7])
    assert not np.isnan(loss)
    assert loss != loss2

def test_cross_entropy():
    """
    run jet model
    Input: [W_j, E_max, alpha, theta_open_p, gamma_bulk, B0, theta_obs, A_eq, E_min]
    """
    fitter = SSC_Fitter("J2011flare", nprocs=4)
    mean_likelihood, mu, Sigma = fitter.cross_entropy_method(n_samples=8, n_elite=3, max_k=1,)
    assert not np.isnan(mean_likelihood)
    assert mu.shape == (9,)
    assert Sigma.shape == (9,9)
    assert mean_likelihood > 0

def test_cross_entropy_gamma_rand():
    """
    run jet model
    Input: [W_j, E_max, alpha, theta_open_p, gamma_bulk, B0, theta_obs, A_eq, E_min]
    """
    fitter = SSC_Fitter("S5flare", nprocs=4, rand_gamma=1)
    mean_likelihood, mu, Sigma = fitter.cross_entropy_method(n_samples=8, n_elite=3, max_k=1,)
    assert not np.isnan(mean_likelihood)
    assert mu.shape == (9,)
    assert Sigma.shape == (9,9)
    assert mean_likelihood > 0
