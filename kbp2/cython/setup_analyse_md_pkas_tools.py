# -*- coding: utf-8 -*-

# call this file with: python setup_analyse_md_pkas_tools.py build_ext --inplace

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("analyse_md_pkas_tools", ["analyse_md_pkas_tools.pyx"])]
)
