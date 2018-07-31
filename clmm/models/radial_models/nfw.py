#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 22:32:04 2018

@author: MattFong
"""
import colossus

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
        super(nfwProfile, self).__init__(zL, mdef, chooseCosmology)
        #Profile1D.__init__(self, zL, mdef, chooseCosmology)
        
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
    
        