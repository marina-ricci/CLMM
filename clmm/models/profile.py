'''General profile class that inherits from Models'''

from models.models import Model
from scipy import integrate
from _profile_utils import sigma_crit

class Profile1D(Model) :
    """
    Generalized superclass for 1D model profiles. It inherits from Models.

    Attributes
    ----------
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
        
        if isinstance(z_lens, float) :
            self.z_lens = zlens
        else :
            raise TypeError('z_lens should be a float')

        if isinstance(mass_definition, str) :
            self.mass_definition = mass_definition
        else :
            # Need to implement test of valid mass_definitions
            raise TypeError('mass_definition should be a string')

        if isinstance(cosmology, string) :
            self.cosmology = cosmology
        else :
            raise TypeError('cosmology should be a string')

        
        # Need to implement r as the independent_var for all Profile1D models
        super().__init__(func, params=params)
    
    def density_3d(self,r):
        """
        3D volume density profile, :math:'\rho'. It is a function of radius.

        Abstract function which must be overwritten by child classes.

        Parameters
        ----------
        r: ndarray
            The radius in units of Mpc.

        Returns
        -------
        rho: ndarray
            :math:'\rho', the 3D density in units of :math:'\mathrm{M}_{\odot}/\mathrm{Mpc}^3'.
            It has the same dimensions as r. 

        """
        pass
    
    def surface_density(self, r):
        """
        Projected surface density profile, :math:'\Sigma'. It is a function of radius.

        Abstract function which must be overwritten by child classes.

        Parameters
        ----------
        r: ndarray
            The radius in units of Mpc.

        Returns
        -------
        sigma: ndarray
            :math:'\Sigma', the surface density in units of :math:'\mathrm{M}_{\odot}/\mathrm{Mpc}^2'. 
            It has the same dimensions as r.
        
        """

        def integrand(r, R):
            ret = 2.0 * r * self.density_3d(r) / numpy.sqrt(r**2 - R**2)
            return ret

        r_use = r.tolist()
        surface_density = (0.0 * r).tolist()
        for i in range(len(r_use)): 
            ReturnResult = scipy.integrate.quad(integrand, r_use[i] + 0.0000001, self.rmax, args = r_use[i], epsrel = 1E-6, limit = 100000)   
            surface_density[i] = ReturnResult[0]
        surface_density = np.array(surface_density)

        return surface_density
    
    def mean_surface_density(self, r, return_sigma=False):
        """
        Mean enclosed surface density, :math:'\bar{\Sigma}'. It is a function of radius.

        Abstract function which must be overwritten by child classes.

        Parameters
        ----------
        r: ndarray
            The radius in units of Mpc.
        return_sigma: bool
            A flag for whether or not to also return the surface density, since it's calculated anyways in this method

        Returns
        -------
        sigmaMean: ndarray
            :math:'\bar{\Sigma}', the mean enclosed surface density in units of :math:'\mathrm{M}_{\odot}/\mathrm{Mpc}^2'
            It has the same dimensions as r.
        sigma
            The surface density. Optional return, since surface density is calculated in the process
        
        """
        R = np.asarray(r)
        add = np.arange(0.0001, R[0], 0.001)
        r = []
        r.extend(add)
        r.extend(R)
        r = np.asarray(r)
            
        sigma = self.surface_density(r)
        sigma_int = integrate.cumtrapz(sigma*r, r, initial = 0)
        sigma_mean = 2.*sigma_int/(r**2.)
            
        sigma = sigma[len(add):]
        sigma_mean = sigma_mean[len(add):]

        if return_sigma:
            return (sigma_mean, sigma)
        else:
            return sigma_mean

    
    def delta_sigma(self, r):
        """
        Difference in mean_surface density and surface density, :math:'\Delta\Sigma = \bar{\Sigma} - \Sigma'. It is a function of radius.

        Abstract function which must be overwritten by child classes.

        Parameters
        ----------
        r: ndarray
            The radius in units of Mpc.

        Returns
        -------
        delta_sigma: ndarray
            :math:'\Delta\Sigma = \bar{\Sigma} - \Sigma' in units of :math:'\mathrm{M}_{\odot}/\mathrm{Mpc}^2'
            It has the same dimensions as r.
        
        """
        sigma_mean, sigma = self.mean_surface_density(r, return_sigma=True)
        delta_sigma = sigma_mean - sigma
        return delta_sigma
    
    def convergence(self, r, z_source):
        """
        Convergence, or dimensionless surface mass density, :math:'\kappa=\Sigma/\Sigma_{crit}'. It is a function of radius.

        Parameters
        ----------
        r: ndarray
            The radius in units of Mpc.
        z_source: float
            Mean effective redshift of the background galaxies.

        Returns
        -------
        kappa: ndarray
            :math:'\kappa', the convergence, which is unitless. It has the same dimensions of r.
        
        """

        sigma_crit = sigma_crit(self.z_lens, z_source, self.mass_definition, self.cosmology)
        sigma = self.surface_density(r)
        kappa = sigma/sigmaC
        return kappa
    
    def mean_convergence(self, r, z_source):
        """
        Mean enclosed convergence, :math:'\bar{\kappa}=\bar{\Sigma}/\Sigma_{crit}'. It is a function of radius.

        Parameters
        ----------
        r: ndarray
            The radius in units of Mpc.
        z_source: float
            Mean effective redshift of the background galaxies.

        Returns
        -------
        mean_kappa: ndarray
            :math:'\bar{\kappa}', the mean convergence enclosed by r, which is unitless. It has the same dimensions as r.
        
        """
        mean_sigma = self.mean_surface_density(r)
        mean_kappa = mean_sigma/sigmaC
        return mean_kappa
    
    def shear(self, r, z_source):
        """
        Tangential shear, :math:'\gamma_{t}=\Delta\Sigma / \Sigma_{crit}'. It is a function of radius.

        Parameters
        ----------
        r: ndarray
            The radius in units of Mpc.
        z_source: float
            Mean effective redshift of the background galaxies.

        Returns
        -------
        gamma: ndarray
            :math:'\gamma_t', the tangential shear, which is unitless. It has the same dimensions as r.
        
        """
        delta_sigma = self.delta_sigma(r)
        gamma = delta_sigma/sigmaC
        return gamma
    
    def reduced_shear(self, r, z_source):
        """
        Tangential reduced shear, :math:'g_t'. It is a function of radius.

        Parameters
        ----------
        r: ndarray
            The radius in units of Mpc.
        z_source: float
            Mean effective redshift of the background galaxies.

        Returns
        -------
        redg: ndarray
            :math:'g_t', the tangential reduced shear, which is unitless. It has the same dimensions as r.
        
        """
        mean_sigma, sigma = self.mean_surface_density(r, return_sigma=True)

        redg = (mean_sigma-sigma)/(sigmaC-sigma)
        return redg    
    
