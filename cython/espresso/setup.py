
#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy as np

ext_params = {}
ext_params['include_dirs'] = [np.get_include(), "../../src/", "../test_bin/src"]
ext_params['extra_compile_args'] = ["-O2"]
ext_params['extra_link_args'] = ["-Wl", "-Wl"]  # TODO: ad-neeeded ignored
ext_params['library_dirs'] = [ "/Users/stefankesselheim/Physik/espresso/cython" ]

ext_modules=[
    Extension("espresso", ["espresso.pyx"], libraries=['espresso_main','tcl'], **ext_params),
    Extension("particle_data", ["particle_data.pyx"], **ext_params),
    Extension("global_variables", ["global_variables.pyx"], **ext_params),
]

setup(
#  name = 'BLAS and LAPACK wrapper',
  name = 'espresso',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules,
)

