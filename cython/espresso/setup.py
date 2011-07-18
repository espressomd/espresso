
#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy as np

ext_params = {}
ext_params['include_dirs'] = [np.get_include(), "../../src/", "../espresso_bin/src"]
ext_params['extra_compile_args'] = ["-fPIC"]
ext_params['extra_link_args'] = ["-Wl", "-Wl"]  # TODO: ad-neeeded ignored
ext_params['library_dirs'] = [ "..", "/usr/lib64" ]

ext_modules=[
    Extension("espresso", ["espresso.pyx"], libraries=['espresso_main','tcl8.5', 'mpi', 'fftw3'], **ext_params),
    Extension("particle_data", ["particle_data.pyx"], libraries=['espresso_main','tcl8.5', 'mpi', 'fftw3'], **ext_params),
    Extension("interaction_data", ["interaction_data.pyx"], libraries=['espresso_main','tcl8.5', 'mpi', 'fftw3'], **ext_params),
    Extension("global_variables", ["global_variables.pyx"], libraries=['espresso_main','tcl8.5', 'mpi', 'fftw3'], **ext_params),
]

setup(
#  name = 'BLAS and LAPACK wrapper',
  name = 'espresso',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules,
)

