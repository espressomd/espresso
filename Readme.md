# ESPResSo

[![GitLab CI](https://gitlab.icp.uni-stuttgart.de/espressomd/espresso/badges/4.1/pipeline.svg)](https://gitlab.icp.uni-stuttgart.de/espressomd/espresso/commits/4.1)
[![codecov](https://codecov.io/gh/espressomd/espresso/branch/4.1/graph/badge.svg)](https://codecov.io/gh/espressomd/espresso/branch/4.1)

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
simulation algorithms, which take into account hydrodynamic
(lattice-Boltzmann) and electrostatic interactions (P3M, ELC, MMMxD).
Rigid bodies can be modelled by virtual site interactions, and it can
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

The [user guide](http://espressomd.org/html/doc/index.html) will
walk you through the basic usage of ESPResSo. Advanced simulation
methods are extensively documented, with examples and links to the
relevant literature. Additional resources can be found on the
homepage at http://espressomd.org/wordpress/documentation/, such
as tutorials and doxygen documentation.

## Installation

Detailed installation instructions for Ubuntu and macOS can be found in the
documentation, section [Installation](http://espressomd.org/html/doc/installation.html).

For most users, we recommend downloading the latest release version of ESPResSo. You
can find it in the [release page](https://github.com/espressomd/espresso/releases),
together with past releases until 4.0. When choosing a release, we recommend that
you get the latest bugfix release in that line. For example, for 4.1 you would like
to use 4.1.4.

### Join the community

Please consider subscribing to our
[mailing list](http://espressomd.org/wordpress/community-and-support/mailing-lists/)
if you're actively using ESPResSo, as we occasionally need community
feedback when making decisions on the future of specific features in
upcoming releases. You'll also get notifications on bugfix releases.

### Please cite us!

If you use ESPResSo to publish scientific results, we would ask you to
acknowledge this usage by mentioning the software with its version number and
[citing the relevant papers](http://espressomd.org/wordpress/about/please-cite-us/).
A number of algorithms in ESPResSo are fairly advanced and unique to ESPResSo.
The authors of these contributions kindly ask you to cite the relevant
publications, as indicated in the documentation. For detailed instructions, see
[How to cite ESPResSo](http://espressomd.org/html/doc/introduction.html#how-to-cite-espresso).

## License

Copyright (C) 2010-2019 The ESPResSo project

Copyright (C) 2002-2010 Max-Planck-Institute for Polymer Research, Theory Group

ESPResSo is free software: you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the
Free Software Foundation, either version 3 of the License, or (at your
option) any later version.

You should have received a [copy](COPYING) of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
