# -*- coding: utf-8 -*-

# call this file with: python setup_combine_md_tools.py build_ext --inplace

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = [Extension("combine_md_tools", ["combine_md_tools.pyx"])]
)
