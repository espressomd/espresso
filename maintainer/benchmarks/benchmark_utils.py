#
# Copyright (C) 2018-2026 The ESPResSo project
#
# This file is part of ESPResSo.
#
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#


CONFIGS = [{"config": "maxset", "mpi": True}, {"config": "default", "mpi": True}, {"config": "empty", "mpi": True}]

# Defintion of all individual benchmark tests
BENCHMARKS = [
    {"file": "lj.py", "args": ["--particles_per_core=1000", "--volume_fraction=0.50"]},
    {"file": "lj.py", "args": ["--particles_per_core=1000", "--volume_fraction=0.02"]},
    {"file": "lj.py", "args": ["--particles_per_core=10000", "--volume_fraction=0.50"]},
    {"file": "lj.py", "args": ["--particles_per_core=10000", "--volume_fraction=0.02"]},
    {"file": "mc_acid_base_reservoir.py", "args": ["--particles_per_core=500"]},
    {"file": "lj.py", "args": ["--particles_per_core=1000", "--volume_fraction=0.10", "--bonds"]},
    {"file": "lj.py", "args": ["--particles_per_core=10000", "--volume_fraction=0.10", "--bonds"]},
    {"file": "p3m.py", "args": ["--particles_per_core=1000", "--volume_fraction=0.25", "--prefactor=4"]},
    {"file": "p3m.py", "args": ["--particles_per_core=10000", "--volume_fraction=0.25", "--prefactor=4"]},
    {"file": "ferrofluid.py", "args": ["--particles_per_core=400"]},
    {"file": "lb.py", "args": ["--particles_per_core=125", "--volume_fraction=0.03", "--lb_sites_per_particle=28"]},
    {"file": "lb.py", "args": ["--box_l=32", "--particles_per_core=0", "--single_precision"]},
    {"file": "lb.py", "args": ["--box_l=32", "--particles_per_core=0"]},
    {"file": "lb.py", "args": ["--box_l=64", "--particles_per_core=0", "--single_precision"]},
    {"file": "lb.py", "args": ["--box_l=64", "--particles_per_core=0"]},
    {"file": "lb.py", "args": ["--box_l=128", "--particles_per_core=0", "--single_precision"]},
    {"file": "lb.py", "args": ["--box_l=128", "--particles_per_core=0"]},
    {"file": "lb.py", "args": ["--box_l=196", "--particles_per_core=0", "--single_precision"]},
    {"file": "lb.py", "args": ["--box_l=196", "--particles_per_core=0"]},
]

CORES_LIST = [1, 2, 4, 8, 12]


def generate_build_parameters():
    "Generate all build parameter configurations."


def generate_test_parameters():
    """Generates all test parameter configurations."""
    params = []
    for benchmark in BENCHMARKS:
        # Replicate the CMake scaling logic
        cores_list = CORES_LIST if benchmark.get('mpi', True) else [1]
        for cores in cores_list:
            # Reframe parameters must be hashable, so lists become tuples
            params.append((benchmark['file'], (benchmark['args']), cores))
    return params