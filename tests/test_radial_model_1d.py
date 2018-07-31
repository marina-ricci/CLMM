
from __future__ import absolute_import, print_function
import numpy as np
from clmm.models.radial_models.radial_model_1d import RadialModel1D
from colossus.halo.mass_so import M_to_R
from colossus.halo.mass_so import densityThreshold

def charOverdensity(Delta, c):
        sigma_c = (Delta/3.)*(c**3.)/(np.log(1. + c) - c/(1. + c))
        return sigma_c #unitless

def nfwprofile(r, params, z_lens, mass_definition):
    Delta = int(mass_definition[:-1])
    M_mass_definition = params['M'] #[M] = M_dot
    c = params['c']
    h = 0.7040000000000001
    r_mass_definition = M_to_R(M_mass_definition*h, z_lens, mass_definition)/1E3/h #Mpc from kpc/h
    rs = r_mass_definition/c
    rho_mass_definition = (densityThreshold(z_lens, mass_definition) * ((1E3)**3.) *(h)**2.)/Delta
    rho_0 = charOverdensity(Delta, c)*rho_mass_definition
    rho = rho_0/((r/rs)*(1+r/rs)**2)
    return rho

def nfwsigma_analytic(r, params, z_lens, mass_definition):
    Delta = int(mass_definition[:-1])
    M_mass_definition = params['M'] #[M] = M_dot
    c = params['c']
    h = 0.7040000000000001
    r_mass_definition = M_to_R(M_mass_definition*h, z_lens, mass_definition)/1E3/h #Mpc from kpc/h
    rs = r_mass_definition/c
    rho_mass_definition = (densityThreshold(z_lens, mass_definition) * ((1E3)**3.) *(h)**2.)/Delta
    rho_0 = charOverdensity(Delta, c)*rho_mass_definition
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



def test_radial_model_1d():
    r = np.linspace(.05, 3., 15)

    sample_profile = RadialModel1D(z_lens=0.25, func=nfwprofile, 
                                       params={'M': 1E15, 'c': 4}, 
                                       z_source = 1.)
    sigma = sample_profile.surface_density(r)
    sigma_analytic = nfwsigma_analytic(r, params={'M': 1E15, 'c': 4}, 
                                       z_lens=0.25, mass_definition='200c')

    print(sigma/sigma_analytic-1.)
    for i in range(len(sigma)):
        np.testing.assert_almost_equal(sigma[i]/sigma_analytic[i]-1.,0.,decimal=5)
