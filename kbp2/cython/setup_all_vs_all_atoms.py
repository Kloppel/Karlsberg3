# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 12:44:23 2012

@author: tmeyer
"""

# call this file with: python setup_all_vs_all_atoms.py build_ext --inplace

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("all_vs_all_atoms", ["all_vs_all_atoms.pyx"])]
)
