"""Tests for the models.cosmo_dependent_models.sigma_crit module"""
from __future__ import absolute_import, print_function

import numpy as np
from numpy.testing import assert_raises
import six
#from clmm import Model, Parameter
from clmm import Model


def assert_block(test_sigma_crit):
    """Block of asserts for type checks in models.cosmo_dependent_models.sigma_crit

    Parameters
    ----------
    test_sigma_crit : SigmaCrit instance
        Instance of the sigma_crit class to run asserts on
    """
