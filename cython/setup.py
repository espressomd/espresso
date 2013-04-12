#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os

import numpy as np

ext_params = {}
ext_params['include_dirs'] = [np.get_include(), "../src"]
ext_params['extra_compile_args'] = ["-fPIC"]
ext_params["depends"] = ["myconfig.pxi"]
#ext_params['extra_link_args'] = ["-Wl", "-Wl"]  # TODO: ad-neeeded ignored

# Convert the cython myconfig.pxi to a python dictionary (we can't just include it)
config = {}
with open("myconfig.pxi") as f:
    myconfig = f.readlines()
    for line in myconfig:
        line = line.split(" = ")
        try:
            config[line[0][4:]] = int(line[1][:-1])
        except ValueError:
            config[line[0][4:]] = line[1][:-1]


libs = ['Espresso','tcl8.5', 'mpi', 'fftw3']
lib_dirs = ['']

ext_modules=[
    Extension("espresso", ["espresso.pyx"], libraries=libs, **ext_params),
    Extension("particle_data", ["particle_data.pyx"], libraries=libs, **ext_params),
    Extension("interaction_data", ["interaction_data.pyx"], libraries=libs, **ext_params),
    Extension("global_variables", ["global_variables.pyx"], libraries=libs, **ext_params),
    Extension("debye_hueckel", ["debye_hueckel.pyx"], libraries=libs, **ext_params),
    Extension("integrate", ["integrate.pyx"], libraries=libs, **ext_params),
    Extension("thermostat", ["thermostat.pyx"], libraries=libs, **ext_params),
    Extension("changeVolume", ["changeVolume.pyx"], libraries=libs, **ext_params),
    Extension("cuda_init", ["cuda_init.pyx"], libraries=libs, **ext_params),
    Extension("invalidateSystem", ["invalidateSystem.pyx"], libraries=libs, **ext_params),
    Extension("code_info", ["code_info.pyx"], libraries=libs, **ext_params),
    Extension("cellsystem", ["cellsystem.pyx"], libraries=libs, **ext_params),
    Extension("analyze", ["analyze.pyx"], libraries=libs, **ext_params),
    Extension("utils", ["utils.pyx"], libraries=libs, **ext_params),
]

if config["LB"] == 1:
    ext_modules.append(Extension("lb", ["lb.pyx"], libraries=libs, **ext_params))

setup(
#  name = 'BLAS and LAPACK wrapper',
    name = 'espresso',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules,
)
