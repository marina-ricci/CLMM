# -*- coding: utf-8 -*-
"""
Created on Thu Nov 19 15:05:13 2015

@author: Matthew Fong
"""


import colossus.cosmology.cosmology as Cosmology
import colossus.halo.profile_dk14 as profile_dk14
import colossus.halo as Halo
import colossus.halo.concentration as hc
from clmm.models import Profile1D


import numpy as np


'''
############################################################################
                                  NFW
############################################################################
'''
class nfwProfile(Profile1D):
    ##### We're going to swap out ``profile'' with ``Profile1D''
    def __init__(self, parameters, zL, mdef, chooseCosmology, esp = None):
        Profile1D.__init__(self, zL, mdef, chooseCosmology)
        
        self.parameters = parameters
        self.M_mdef = parameters['M'].value #[M] = M_dot
        self.c = parameters['c'].value
        self.zL = zL
        self.mdef = mdef
        self.chooseCosmology = chooseCosmology
        #input [M] = M_dot/h
        self.r_mdef = Halo.mass_so.M_to_R(self.M_mdef*self.cosmo.h, self.zL, self.mdef)/1E3/self.cosmo.h #Mpc from kpc/h
        self.Delta = int(mdef[:-1])

        #[rho_mdef] = M_dot/Mpc^3 from M_{\odot}h^2/kpc^3
        self.rho_mdef = (Halo.mass_so.densityThreshold(self.zL, self.mdef) * ((1E3)**3.) *(self.cosmo.h)**2.)/self.Delta

        self.rs = self.r_mdef/self.c #Mpc
        self.profile = 'nfw'
        if esp == None:
            self.esp = 1E-5
        else:
            self.esp = esp
        return
        
    '''
    ############################################################################
                               Analytic NFW profile
    ############################################################################
    '''
    def charOverdensity(self):
        sigma_c = (self.Delta/3.)*(self.c**3.)/(np.log(1. + self.c) - self.c/(1. + self.c))
        return sigma_c #unitless
    
    def nfwrho(self, R):
        #R in Mpc
        #[sigma_c] = unitless
        #[rho_mdef] = M_dot / Mpc^3
        const =  self.rho_mdef * self.charOverdensity() 
        rhoForm = 1./( (R/self.rs) * (1. + R/self.rs)**2.)
        return (const * rhoForm)
        
    def nfwSigma(self, r):
        #[r] = Mpc
        
        #[rs] = Mpc        
        rs = self.rs
        expSig = np.empty(len(r))
        for i in range(len(r)):
            if r[i]/self.rs < 1.0 - self.esp:
                expSig[i] = (1./((r[i]/rs)**2. - 1.))*(1. - 2.*np.arctanh(np.sqrt(\
                (1. - (r[i]/rs))/(1. + (r[i]/rs))))/np.sqrt(1. - (r[i]/rs)**2.))
            #if r[i]/rs == 1.0:
            #    expSig[i] = 1./3.
            if r[i]/rs >= 1.0 - self.esp and r[i]/rs <= 1.0 + self.esp:
                expSig[i] = 1./3.
            if r[i]/self.rs > 1.0 + self.esp:
                expSig[i] = (1./((r[i]/rs)**2. - 1.))*(1. - 2.*np.arctan(np.sqrt(\
                ((r[i]/rs) - 1.)/(1. + (r[i]/rs))))/np.sqrt((r[i]/rs)**2. - 1.))
        const = 2.*(self.rs)*self.charOverdensity()*(self.rho_mdef)
        #[Sigma] = M_dot / Mpc^2 
        return (expSig * const)
    
    def nfwSigmaMean(self, r):
        #[r] = Mpc
        x = r/self.rs
        const = 4.*self.rs*self.charOverdensity()*self.rho_mdef
        if type(x) is np.ndarray:
            #print x
            fnfw = np.piecewise(x,[x<1., x==1., x>1.], \
                        [lambda x:(2./np.sqrt(1.-x*x)*np.arctanh(np.sqrt((1.-x)/(1.+x)))+np.log(x/2.))/(x*x), \
                         lambda x:(1.+np.log(0.5)), \
                         lambda x:(2./np.sqrt(x*x-1.)*np.arctan2(np.sqrt(x-1.),np.sqrt(1.+x))+np.log(x/2.))/(x*x)])
            return const*fnfw
    
        else:
            if x<1:
                return const*(2./np.sqrt(1.-x*x)*np.arctanh(np.sqrt((1.-x)/(1.+x)))+np.log(x/2.))/(x*x)
            elif x==1:
                return const*(1.+np.log(0.5))
            else:
                return const*(2./np.sqrt(x*x-1.)*np.arctan2(np.sqrt(x-1.),np.sqrt(1.+x))+np.log(x/2.))/(x*x)
        
        
    
    def nfwDeltaSigma(self, r):
        #[r] = Mpc
        rs = self.rs
        
        x = r/rs
        expG = np.empty(len(x))
        for i in range(len(x)):
            if x[i] < 1.0 - self.esp:
                expG[i] = 8.0*np.arctanh(np.sqrt((1.0-(x[i]))/(1.0+(x[i]))))/(((x[i])**2.0)*np.sqrt(1.0-(x[i])**2.0))\
                + (4.0/((x[i])**2.0))*np.log((x[i])/2.0) \
                - 2.0/((x[i])**2.0 - 1.0)\
                + 4.0*np.arctanh(np.sqrt((1.0-(x[i]))/(1.0+(x[i]))))/(((x[i])**2.0 - 1.0)*np.sqrt(1.0-(x[i])**2.0))
            if x[i] > 1.0 - self.esp and x[i] < 1.0 + self.esp:
                expG[i] = (10.0 / 3.0) + (4.0 * np.log(1.0/2.0))
            if x[i] > 1.0 + self.esp:
                expG[i] = 8.0*np.arctan(np.sqrt(((x[i])-1.0)/(1.0+(x[i]))))/(((x[i])**2.0)*np.sqrt((x[i])**2.0 - 1.0))\
                + (4.0/((x[i])**2.0))*np.log((x[i])/2.0) \
                - 2.0/((x[i])**2.0 - 1.0)\
                + 4.0*np.arctan(np.sqrt(((x[i])-1.0)/(1.0+(x[i]))))/(((x[i])**2.0 - 1.0)**(3.0/2.0))
        
        
        #[rho_mdef] = M_dot/Mpc^3 
        #[charOverdensity] = unitless
        charOverdensity = self.charOverdensity()
        Const = (rs) * (self.rho_mdef) * (charOverdensity)
        # [Const] = 1 
        Const = Const
        return Const*expG
    
        
        
    
'''
############################################################################
                               Truncated NFW
############################################################################
'''
class nfwBMOProfile(Profile1D):
    
    def __init__(self, parameters, zL, n, mdef, chooseCosmology, Tau = None, cM_relation = None, esp = None):
        Profile1D.__init__(self, zL, mdef, chooseCosmology)
        
        cosmo = Cosmology.setCosmology(chooseCosmology)
        self.parameters = parameters
        self.M_mdef = parameters['M'].value #M200 input in M_dot/h
        if cM_relation == True:
            self.c = hc.concentration(self.M_mdef*cosmo.h, self.mdef, self.zL)
        else:
            self.c = parameters['c'].value
        self.zL = zL
        self.n = n #sharpness of truncation (n = 1 or 2)
        if Tau == None: #dimensionless truncation radius (T = rt/rvir => fit =2.6)
            self.T = 2.6
        else:
            self.T = Tau
        self.r_mdef = Halo.mass_so.M_to_R(self.M, self.zL, self.mdef)/1E3 #Mpc/h
        self.rs = self.r_mdef/self.c
        self.rt = self.T*self.rs
        self.mdef = mdef
        self.chooseCosmology = chooseCosmology
        self.G, self.v_c, self.H2, self.cosmo = self.calcConstants()
        self.Delta = int(mdef[:-1])
        #[rho_mdef] = M_dot Mpc^3 from M_{\odot}h^2/kpc^3
        self.rho_mdef = (Halo.mass_so.densityThreshold(self.zL, self.mdef) * 1E9 *(self.cosmo.h)**2.)/self.Delta
        
        self.profile = 'nfwBMO'
        if esp == None:
            self.esp = 1E-5
        else:
            self.esp = esp
        
        return
        #super(nfwProfile, self).__init__()

        
    '''
    ############################################################################
                               Analytic BMO profile
    ############################################################################
    '''
    def mnfw(self):
        m = (np.log(1. + self.c) - self.c/(1. + self.c))
        return m

    def densityParameter(self):
        rho_s = self.M/(4.*np.pi*(self.rs**3.)*self.mnfw())
        #[rho_s] = M_dot / Mpc^3 from M_dot h^2 / Mpc^3
        rho_s = rho_s * (self.cosmo.h**2.)
        #[rho_s] = M_dot / Mpc^3
        return rho_s
        
    def bmorho(self, R):
        #R in Mpc/h
        #[rho_s] = M_dot / Mpc^3
        rho_s = self.densityParameter()
        rho_nfw = (rho_s/((R/self.rs)*(1.+R/self.rs)**2.))*((self.rt**2.)/(R**2. + self.rt**2.))**self.n
        #[rho_nfw] = M_dot / Mpc^3
        return rho_nfw
    
    def F(self, r):
        #[r] = Mpc/h
        func = np.empty(len(r))
        x = r/self.rs
        for i in range(len(x)):
            if x[i] < 1.0:
                func[i] = (1./(np.sqrt(1. - x[i]**2.)))*np.arctanh(np.sqrt(1.-x[i]**2.))
            if x[i] > 1.0:
                func[i] = (1./(np.sqrt(x[i]**2. - 1.)))*np.arctan(np.sqrt(x[i]**2.-1.))
        return func
    
    def L(self, r):
        #[r] = Mpc/h
        x = r/self.rs
        T = self.T
        func = np.log(x/(np.sqrt(T**2. + x**2.) + T))
        return func
        
    def bmoSigma(self, r):
        #[r] = Mpc/h
        x = r/self.rs
        T = self.T
        #[densityParameter] = M_dot / Mpc^3
        #[Const] = M_dot/Mpc^2 from M_dot / h / Mpc^2
        Const = (4.*self.densityParameter()*self.rs)*((T**2.)/(2.*(T**2. + 1.)**2.))
        Const = Const/self.cosmo.h
        if self.n == 1:
            func = ((T**2. + 1.)/(x**2. - 1.))*(1.-self.F(r)) + 2.*self.F(r)\
                    - np.pi/(np.sqrt(T**2. + x**2.))\
                    +((T**2. - 1.)/(T*(np.sqrt(T**2. + x**2.))))*self.L(r)
        if self.n == 2:
            Const = Const*(T**2./(2.*(T**2.+1.)))
            func = ((2.*(T**2. + 1.))/(x**2. - 1.))*(1. - self.F(r))\
                    + 8.*self.F(r) + (T**4. - 1.)/((T**2.)*(T**2. + x**2.))\
                    - np.pi*((4.*(T**2. + x**2.)\
                    + T**2. + 1.)/((T**2. + x**2.)**(3./2.)))\
                    + (((T**2.)*(T**4. - 1.) + (T**2. + x**2.)*(3.*T**4. - 6.*T**2. - 1.))\
                    / ((T**3.)*(T**2. + x**2.)**(3./2.)))*self.L(r)
        #[Sigma] = M_dot/Mpc^2
        return Const*func
        
        
    def bmoSigmaMean(self, r):
        x = r/self.rs
        T = self.T
        Const = (4.*self.densityParameter()*self.rs)
        #[Const] = M_dot / Mpc^2 from M_dot / h / Mpc^2
        Const = Const / self.cosmo.h
        if self.n == 1:
            func = ((T**2.)/((x**2.)*(T**2. + 1.)**2.))\
                    *(\
                    (T**2. + 2.*x**2. + 1.)*self.F(r)\
                    + T*np.pi + (T**2. - 1.)*np.log(T)\
                    + np.sqrt(T**2. + x**2.)*(-np.pi + ((T**2. - 1.)/T)*self.L(r))\
                    )
        if self.n == 2:
            func = ((T**4.)/(2.*(x**2.)*(T**2. + 1.)**3.))\
                    *(2.*(T**2. + 4.*x**2. - 3.)*self.F(r)\
                    +(1./T)*(\
                    np.pi*(3.*T**2. - 1.)\
                    + 2.*T*(T**2. - 3.)*np.log(T)\
                    )\
                    + (1./((T**3.)*(np.sqrt(T**2. + x**2.))))\
                    *(\
                    (-(T**3.))*np.pi*(4.*x**2. + 3.*T**2. - 1.)\
                    + ((2.*T**4.)*(T**2. - 3.)\
                    + (x**2.)*(3.*T**4. - 6.*T**2. - 1.))*self.L(r)\
                    ))
        #[SigmaMean] = M_dot / Mpc^2
        return Const*func   



'''
############################################################################
                             Numerical dk14
############################################################################
'''

class dkProfile(Profile1D):
    
    def __init__(self, parameters, zL, mdef, chooseCosmology, part = None, \
                 se = None, be = None, cM_relation = None):
        Profile1D.__init__(self, zL, mdef, chooseCosmology)
        
        
        self.M_mdef = parameters['M'].value #M200
        
        if cM_relation == True:
            self.c = hc.concentration(self.M*self.cosmo.h, self.mdef, self.zL)
            #self.c = 3.614*((1+self.zL)**(-0.424))*(self.M/self.cosmo.h/1E14)**(-0.105)
        else:
            self.c = parameters['c'].value
        
        if se is None:
            self.se = 1.5
        else:
            self.se = se
        if be is None:
            self.be = 1.0
        else:
            self.be = be
        
        
        self.zL = zL
        self.mdef = mdef
        if part is not None:
            self.part = part
        else:
            self.part = 'both'
        #[rs] = Mpc/h
        self.r_mdef = Halo.mass_so.M_to_R(self.M_mdef*self.cosmo.h, self.zL, self.mdef)/1E3/self.cosmo.h #Mpc
        self.rs = self.r_mdef/self.c #Mpc
        self.Delta = int(mdef[:-1])
        #[rho_mdef] = M_dot Mpc^3 from M_{\odot}h^2/kpc^3
        self.rho_mdef = (Halo.mass_so.densityThreshold(self.zL, self.mdef) * 1E9 *(self.cosmo.h)**2.)/self.Delta
        
        '''
        self.dk14Prof = HaloDensityProfile.DK14Profile(M = self.M_mdef, c = self.c, z = \
                                                       self.zL, mdef = self.mdef, \
                                                       be = self.be, se = self.se, \
                                                       part = self.part)
        '''
        self.dk14Prof = profile_dk14.getDK14ProfileWithOuterTerms(M = self.M_mdef*self.cosmo.h, c = self.c, z = self.zL, 
                                     mdef = self.mdef, 
                                     outer_term_names = ['pl'])
        #self.dk14Prof.par.se = self.se
        #self.dk14Prof.par.be = self.be
        self.rmaxMult = 2.
        
        #self.dk14Prof.par.rs = self.rs*1E3 #[rs] = kpc/h from Mpc/h
        #self.dk14Prof.selected = 'by_accretion_rate' #beta = 6 gamma = 4; more accurate results
        self.dk14Prof.selected = 'by_mass' 
        self.profile = 'dk'
        #super(dkProfile, self).__init__()
        
        return
    
    def dkrho(self,r):
        #input [R] = Mpc
        R = r*1E3*self.cosmo.h #[R] = kpc/h from Mpc for Diemer input
        rho = self.dk14Prof.density(R) *1E9 #[rho] = M_dot h^2 / Mpc^3 from M_{\odot} h^2/ kpc^3
        #[rho] = M_dot / Mpc^3 from M_dot h^2 / Mpc^3
        rho = rho * (self.cosmo.h**2.)
        return rho 
    
    def surface_density(self,r):
        """
        DK14 projected surface density profile, :math:'\Sigma'. It is a function of radius. See Diemer & Kravtsov (2014)

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
        #input [R] = Mpc
        #[R] = kpc/h from Mpc for Diemer input
        R = r*1E3*self.cosmo.h
        #self.dk14Prof.rmax = R[-1]*self.rmaxMult
        #[surfaceDensity] = M_dot h/Mpc^2 from M_{\odot} h/kpc^2
        SigmaDiemer = self.dk14Prof.surfaceDensity(R) * 1E6
        #[Sigma] = M_dot / Mpc^2 from `M_{\odot} h/Mpc^2`
        SigmaDiemer = SigmaDiemer * self.cosmo.h
        return SigmaDiemer #[Sigma] = M_dot/Mpc^2
    
