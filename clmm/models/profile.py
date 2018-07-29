'''General profile class that inherits from Models'''

import Models

class Profile1D(Models) :
    """
    Generalized superclass for 1D model profiles. It inherits from Models.

    Attributes
    -----------------------------------------------
    z_lens: float
        Redshift of the lensing object

    mass_definition: string
        Definition of mass, e.g. 'm200c'

    cosmology: string
        Label of cosmology, e.g. 'WMAP7'

    func: callable
        Functional form of the model, should be wrapped by the class

    params: list
        List of Parameter objects (default to None)
    
    """
    def __init__(self, z_lens, mass_definition, cosmology, func, params=None) :
        pass
    
    def density_3d(self,r):
        """
        3D volume density profile, :math:'\rho'. It is a function of radius.

        Abstract function which must be overwritten by child classes.

        Parameters
        -----------------------------------------------
        r: ndarray
            The radius in units of Mpc.

        Returns
        -----------------------------------------------
        rho: ndarray
            :math:'\rho', the 3D density in units of :math:'\mathrm{M}_{\odot}/\mathrm{Mpc}^3'.
            It has the same dimensions as r. 

        """
        pass
    
    def surface_density(self, r):
        """
        Surface density profile, :math:'\Sigma'. It is a function of radius.

        Abstract function which must be overwritten by child classes.

        Parameters
        -----------------------------------------------
        r: ndarray
            The radius in units of Mpc.

        Returns
        -----------------------------------------------
        sigma: ndarray
            :math:'\Sigma', the surface density in units of :math:'\mathrm{M}_{\odot}/\mathrm{Mpc}^2'. 
            It has the same dimensions as r.
        
        """

        def integrand(r, R):
            ret = 2.0 * r * self.density(r) / numpy.sqrt(r**2 - R**2)
            return ret

        r_use = r.tolist()
        surfaceDensity = (0.0 * r).tolist()
        for i in range(len(r_use)): 
            ReturnResult = scipy.integrate.quad(integrand, r_use[i] + 0.0000001, self.rmax, args = r_use[i], epsrel = 1E-6, limit = 100000)   
            surfaceDensity[i] = ReturnResult[0]
        surfaceDensity = np.array(surfaceDensity)

        return surfaceDensity
    
    def mean_surface_density(self, r):
        """
        Mean enclosed surface density, :math:'\bar{\Sigma}'. It is a function of radius.

        Abstract function which must be overwritten by child classes.

        Parameters
        -----------------------------------------------
        r: ndarray
            The radius in units of Mpc.

        Returns
        -----------------------------------------------
        mean_sigma: ndarray
            :math:'\bar{\Sigma}', the mean enclosed surface density in units of :math:'\mathrm{M}_{\odot}/\mathrm{Mpc}^2'
            It has the same dimensions as r.
        
        """
        pass
    
    def delta_sigma(self, r):
        """
        Difference in mean_surface density and surface density, :math:'\Delta\Sigma = \bar{\Sigma} - \Sigma'. It is a function of radius.

        Abstract function which must be overwritten by child classes.

        Parameters
        -----------------------------------------------
        r: ndarray
            The radius in units of Mpc.

        Returns
        -----------------------------------------------
        delta_sigma: ndarray
            :math:'\Delta\Sigma = \bar{\Sigma} - \Sigma' in units of :math:'\mathrm{M}_{\odot}/\mathrm{Mpc}^2'
            It has the same dimensions as r.
        
        """
        pass
    
    def convergence(self, r, z_source):
        """
        Convergence, or dimensionless surface mass density, :math:'\kappa=\Sigma/\Sigma_{crit}'. It is a function of radius.

        Parameters
        -----------------------------------------------
        r: ndarray
            The radius in units of Mpc.
        z_source: float
            Mean effective redshift of the background galaxies.

        Returns
        -----------------------------------------------
        kappa: ndarray
            :math:'\kappa', the convergence, which is unitless. It has the same dimensions of r.
        
        """
        pass
    
    def mean_convergence(self, r, z_source):
        """
        Mean enclosed convergence, :math:'\bar{\kappa}=\bar{\Sigma}/\Sigma_{crit}'. It is a function of radius.

        Parameters
        -----------------------------------------------
        r: ndarray
            The radius in units of Mpc.
        z_source: float
            Mean effective redshift of the background galaxies.

        Returns
        -----------------------------------------------
        mean_kappa: ndarray
            :math:'\bar{\kappa}', the mean convergence enclosed by r, which is unitless. It has the same dimensions as r.
        
        """
        pass
    
    def shear(self, r, z_source):
        """
        Tangential shear, :math:'\gamma_{t}=\Delta\Sigma / \Sigma_{crit}'. It is a function of radius.

        Parameters
        -----------------------------------------------
        r: ndarray
            The radius in units of Mpc.
        z_source: float
            Mean effective redshift of the background galaxies.

        Returns
        -----------------------------------------------
        gamma: ndarray
            :math:'\gamma_t', the tangential shear, which is unitless. It has the same dimensions as r.
        
        """
        pass
    
    def reduced_shear(self, r, z_source):
        """
        Tangential reduced shear, :math:'g_t'. It is a function of radius.

        Parameters
        -----------------------------------------------
        r: ndarray
            The radius in units of Mpc.
        z_source: float
            Mean effective redshift of the background galaxies.

        Returns
        -----------------------------------------------
        g: ndarray
            :math:'g_t', the tangential reduced shear, which is unitless. It has the same dimensions as r.
        
        """
        pass
    
    