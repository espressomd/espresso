# ESPResSo

[![GitLab CI](https://gitlab.icp.uni-stuttgart.de/espressomd/espresso/badges/python/pipeline.svg)](https://gitlab.icp.uni-stuttgart.de/espressomd/espresso/commits/python)
[![codecov](https://codecov.io/gh/espressomd/espresso/branch/python/graph/badge.svg)](https://codecov.io/gh/espressomd/espresso)

This is the Molecular Dynamics software ESPResSo ("Extensible
Simulation Package for Research on Soft Matter Systems").

ESPResSo is a highly versatile software package for performing and
analyzing scientific Molecular Dynamics many-particle simulations of
"coarse-grained" bead-spring models as they are used in soft-matter
research in physics, chemistry and molecular biology. It can be used
to simulate systems as for example polymers, liquid crystals,
colloids, ferrofluids and biological systems such as DNA and lipid
membranes.

In "coarse-grained" models, a whole group of atoms or molecules are
treated as a single bead.  Although many details are coarse-grained
away in these models, they can often predict qualitative properties,
such as for example the scaling behavior of a system, and can give
insight into theoretical models.  Due to the drastic reduction of
degrees of freedom, coarse-graining allows to investigate systems
which would be out of reach of the commonly used atom-based
simulations, due to the large time- and length scales of the studied
processes.

ESPResSo is capable of doing classical Molecular dynamics simulations
of many types of systems in different statistical ensembles (NVE, NVT,
NPT) and non-equilibrium situations, using standard potentials such as
the Lennard-Jones or Morse potential. It contains many advanced
simulation algorithms, which take into account hydrodynamic (lattice
Boltzmann) and electrostatic interactions (P3M, ELC, MMMxD). Rigid
bodies can be modelled by virtual site interactions, and it can
integrate rotationally non-invariant particles.

ESPResSo is free, open-source software (GPL). It is parallelized and
can be employed on desktop machines, convenience clusters as well as
on supercomputers with hundreds of CPUs. The parallel code is
controlled via the scripting language Python, which gives the software
its great flexibility and allows for many unconventional simulation
protocols, as are often required when studying coarse-grained models.

ESPResSo is used in scientific working groups all over the world both
as a production platform as well as a research platform for developing
new algorithms and methods and designing new models for coarse-grained
simulations.  It is mainly developed at the Institute for
Computational Physics of the University of Stuttgart, but has
contributors from all over the world.


## Documentation

You can find documentation on how to compile, use and develop ESPResSo
on the homepage at [http://espressomd.org/html/doc/index.html](http://espressomd.org/html/doc/index.html).

ESPResSo is intended to be used by people that have proper knowledge
of simulation techniques and know how to use them. We do not take
responsibility if you use ESPResSo to create and publish bogus
results. You have been warned!

## PLEASE CITE US!

If you use ESPResSo and obtain scientific results that you publish, we
would ask you to acknowledge the usage of ESPResSo by citing the relevant papers:
http://espressomd.org/wordpress/about/please-cite-us/

A number of algorithms in ESPResSo are fairly advanced and unique to
ESPResSo. The authors of these contributions kindly ask you to cite the
relevant publications, as indicated in the ESPResSo User's Guide.

## License

Copyright (C) 2010-2018 The ESPResSo project
Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
Max-Planck-Institute for Polymer Research, Theory Group
  
This file is part of ESPResSo.
  
ESPResSo is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.
 
ESPResSo is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.
 
You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

