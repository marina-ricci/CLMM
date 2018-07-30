'''profile utilities including sigma crit functionality and get cosmology function'''

import astropy.cosmology as apycosmo
from astropy import constants
import numpy as np


class sigma_crit
    """
    This calculates a sigma crit for a given cosmology, mass definition, source and lens redshifts.
    
    Attributes
    ----------

    z_lens : float
        Redshift of the lensing object

    z_source : 

    mass_definition : str
        Definition of mass, e.g. 'm200c'
 
    cosmology : str
        Label of cosmology

    """
    

    def __init__(self, z_lens, z_source, mass_definition, cosmology) :

        """ 
        Parameters
        ----------

        z_lens : float
            Redshift of the lensing object

        z_source : 

        mass_definition : str
            Definition of mass, e.g. 'm200c'

        cosmology : str
            Label of cosmology

        Methods
        -------

        calculate_sigma_crit

        beta_function

        beta_discrete_mean
        
        beta_square_discrete_mean


        """ 
        self.cosmology = _get_cosmology( cosmology ) 

    def calculate_sigma_crit(self) :
        
        return constants.c.to('Mpc/s') * constants.c.to('Mpc/s') / \
            (4. * np.pi *  constants.G.to('Mpc3 / (Msun  s2)') * \
             self.cosmology.Dl * self.beta_function(z_source) )



    def _beta_function(self, z_source) :
        
        return self.cosmo.angular_diameter_distance( self.z_lens ) / self.cosmo.angular_diameter_distance( self.z_source )

        pass

    def _beta_discrete_mean( z_source ) :
        pass

    def _beta_square_discrete_mean( z_source ) :
        pass


def _angular_diameter_distance_two_objects(z_lens, z_source, cosmology) :
    """
    Hogg+00 Eqn (19)

    """
    pass


def _get_cosmology( cosmology ) :
    """
    Currently only can parse for astropy cosmology options
    
    cosmology : str
        Label of cosmology


    """

    pass

    
