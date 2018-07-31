"""
Collection of weak lensing functions 

"""

from astropy import constants
import numpy as np


        
def sigma_crit(d_a_lens, d_a_source, d_a_lens_source) :
    """ 
    Parameters
    ----------

    returns sigma crit in units of Msun / Mpc^2
    """ 
    
    return constants.c.to('Mpc/s') * constants.c.to('Mpc/s') / \
        (4. * np.pi *  constants.G.to('Mpc3 / (Msun  s2)') * \
         d_a_lens * beta_geometric_lensing_efficiency(d_a_source, d_a_lens_source) )

def beta_geometric_lensing_efficiency(d_a_source, d_a_lens_source) :
    """
    Schrabback2016 arXiv:1611.03866
    """
    
    if d_a_lens_source < 0:
        return 0.
    else:
        return d_a_lens_source / d_a_source 

def beta_discrete_mean(d_a_lens, d_a_source, d_a_lens_source, weight) :
    """
    Arithmetic mean for discrete lensing efficiency <beta> for given source redshift or redshifts
    Schrabback2016 arXiv:1611.03866
    
    Parameters
    ----------
    weight : shape weight

    z_source : float or array-like of floats
        Source redshift or redshifts
    
    """

    if np.iterable(d_a_source) and all(isinstance(d_a_s, float) for d_a_s in d_a_source) :
        betaxweight = beta_geometric_lensing_efficiency(d_a_source, d_a_lens_source)*weight
        return sum(betaxweight)/sum(weight)
    
        #return np.mean(beta_geometric_lensing_efficiency(d_a_source, d_a_lens_source))

    elif isinstance(d_a_source, float) :
        return beta_geometric_lensing_efficiency(d_a_source, d_a_lens_source)

    else :
        raise TypeError("d_a_source must be float or array-like of floats ")
        


def beta_square_discrete_mean(d_a_lens, d_a_source, d_a_lens_source, weight) :
    """
    Mean of the square of the discrete lensing efficiency <beta^2>
    Schrabback2016 arXiv:1611.03866
    
    Parameters
    ----------
    weight : shape weight

    z_source : float or array-like of floats
        Source redshift or redshifts

    """

    if np.iterable(d_a_source) and all(isinstance(d_a_z, float) for d_a_z in d_a_source) :
        beta2xweight = (beta_geometric_lensing_efficiency(d_a_source, d_a_lens_source)**2.)*weight
        return sum(beta2xweight)/sum(weight)
        #return np.mean(beta_geometric_lensing_efficiency(d_a_source, d_a_lens_source)**2.)

    elif isinstance(d_a_source, float) :
        return beta_geometric_lensing_efficiency(d_a_source, d_a_lens_source)**2.

    else :
        raise TypeError("d_a_source must be float or array-like of floats ")

