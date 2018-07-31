#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 30 22:33:55 2018

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
                             Numerical dk14
############################################################################
'''

class dkProfile(Profile1D):
    
    def __init__(self, parameters, zL, mdef, chooseCosmology, part = None, \
                 se = None, be = None, cM_relation = None):
        Profile1D.__init__(self, zL, mdef, chooseCosmology)
        
        
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
        #self.rmaxMult = 2.
        
        #self.dk14Prof.par.rs = self.rs*1E3 #[rs] = kpc/h from Mpc/h
        #self.dk14Prof.selected = 'by_accretion_rate' #beta = 6 gamma = 4; more accurate results
        #self.dk14Prof.selected = 'by_mass' 
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
    