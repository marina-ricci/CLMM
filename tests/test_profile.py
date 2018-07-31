
from __future__ import absolute_import, print_function
import numpy as np
import scipy
#import matplotlib.pyplot as plt
from clmm import Model, Parameter, Profile1D

def nfwprofile(r, params):
    rs, rho_0 = params
    rho = rho_0/((r/rs)*(1+r/rs)**2)
    return rho

def nfwsigma_analytic(r, rs, rho_0):
    k = 2*rho_0

    if isinstance(r,np.ndarray):
        f = np.piecewise(r, [r>rs, r<rs,r==rs], [lambda r: (rs*k/(((r/rs)**2) - 1))*(1 - 2*(np.arctan((((r/rs) - 1)/((r/rs) + 1))**0.5)/(((r/rs)**2) -1)**0.5)), lambda r:(rs*k/(((r/rs)**2) - 1))*(1 - 2*(np.arctanh((((-r/rs) + 1)/((r/rs) + 1))**0.5)/((1 - ((r/rs)**2))**0.5))), lambda r: rs*k/3.])
    else:
        if r>rs:
            f = (rs*k/(((r/rs)**2) - 1))*(1 - 2*(np.arctan((((r/rs) - 1)/((r/rs) + 1))**0.5)/(((r/rs)**2) -1)**0.5))
        elif r<rs:
            f = (rs*k/(((r/rs)**2) - 1))*(1 - 2*(np.arctanh((((-r/rs) + 1)/((r/rs) + 1))**0.5)/((1 - ((r/rs)**2))**0.5)))
        else:
            f = rs*k/3.
    return f



def test_profile():
    r = np.linspace(.05, 3., 15)

    sample_profile = profile.Profile1D(z_lens=0.25, mass_definition='m200c', cosmology='WMAP7', func=nfwprofile, params=[.5, 1.E15])
    sigma = sample_profile.surface_density(r)
    sigma_analytic = nfwsigma_analytic(r, rs=.5, rho_0=1.E15)

    print(sigma/sigma_analytic-1.)
    for i in range(len(sigma)):
        np.testing.assert_almost_equal(sigma[i]/sigma_analytic[i]-1.,0.,decimal=5)

    
