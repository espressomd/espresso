# ESPResSo

[![GitLab CI](https://gitlab.icp.uni-stuttgart.de/espressomd/espresso/badges/doc/pipeline.svg)](https://gitlab.icp.uni-stuttgart.de/espressomd/espresso/commits/doc)
[![codecov](https://codecov.io/gh/espressomd/espresso/branch/python/graph/badge.svg)](https://codecov.io/gh/espressomd/espresso)

This is the Molecular Dynamics software ESPResSo ("Extensible
Simulation Package for the Research on Soft Matter Systems").

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
would ask you to acknowledge the usage of ESPResSo by referencing

      Hans-Jörg Limbach, Axel Arnold, Bernward A. Mann and Christian
      Holm. 
      "ESPResSo - An Extensible Simulation Package for Research on
      Soft Matter Systems". 
      Comput. Phys. Commun. 174(9) (704-727), 2006.

For a BibTeX record, please refer to the key "espresso" in
doc/ug/citations.bib.

A number of algorithms in ESPResSo are fairly advanced and unique to
ESPResSo. The authors of these contributions kindly ask you to cite the
relevant publications, as indicated in the ESPResSo User's Guide
(see http://espressomd.org or doc/ug/ug.pdf). For your convenience, the
BibTeX records are compiled in doc/ug/citations.bib.

## Licence

The following files are licensed and copyrighted as indicated below:

```
doc/dg/figs/bond_angle.fig
doc/dg/figs/datastorage.gif
doc/dg/figs/dihedral_angle.fig
doc/dg/figs/directions.fig
doc/dg/figs/elc_errordist.gif
doc/dg/figs/ghost_cells.fig
doc/dg/figs/ghost_communication.fig
doc/dg/figs/linked_cells.fig
doc/dg/figs/logo.png
doc/dg/figs/move_to_p_buf.fig
doc/dg/figs/particles.fig
doc/ug/figures/diamond.png
doc/ug/figures/dihedral-angle.fig
doc/ug/figures/dihedral-angle.pdf
doc/ug/figures/fullerene.png
doc/ug/figures/hbond.fig
doc/ug/figures/hbond.pdf
doc/ug/figures/logo.png
doc/ug/figures/maggs-charge-assignment.pdf
doc/ug/figures/maggs-initial-scheme.pdf
doc/ug/figures/maggs-rotation.pdf
doc/ug/figures/nacl-rdf.pdf
doc/ug/figures/salt.png
doc/tutorials/tut2/figures/data/neutral-rho.data
doc/tutorials/tut2/figures/data/nonneutral-rho.data
doc/tutorials/tut2/figures/data/rdf_from_melt_00.data
doc/tutorials/tut2/figures/data/rdf_from_melt_10.data
doc/tutorials/tut2/figures/data/rdf_from_melt_11.data
doc/tutorials/tut2/figures/data/rdf_lj_00.data
doc/tutorials/tut2/figures/nacl-rdf.pdf
doc/tutorials/tut2/figures/nacl.plot
doc/tutorials/tut2/figures/neutral-rho.pdf
doc/tutorials/tut2/figures/nonneutral-rho.pdf
doc/tutorials/tut2/figures/salt.png
packages/mbtools/doc/colloid_model.jpg
packages/mbtools/doc/cylinder_membrane.jpg
packages/mbtools/doc/flat_membrane.jpg
packages/mbtools/doc/protein_model.jpg
packages/mbtools/doc/sphere_membrane.jpg
packages/mbtools/doc/torus_membrane.jpg
packages/mbtools/doc/wrapped_colloid_densitymap.jpg
packages/mbtools/examples/forcetables/9_095_11.tab
packages/mbtools/examples/forcetables/n9_c140_22.tab
packages/mbtools/examples/forcetables/n9_c160_22.tab
packages/mbtools/examples/forcetables/sr_e10_c25.tab
samples/adress/bond_tetra.tab
samples/adress/cg_ic_tetra.tab
samples/adress/cg_tetra.tab
samples/adress/thermo_force.tab
samples/maggs_correct_rdf.dat
samples/object-in-fluid/object-in-fluidUG/figures/arealocal.eps
samples/object-in-fluid/object-in-fluidUG/figures/bending.eps
samples/object-in-fluid/object-in-fluidUG/figures/bloodCell.eps
samples/object-in-fluid/object-in-fluidUG/figures/stretching.eps
samples/object-in-fluid/object-in-fluidUG/figures/volume.eps
testsuite/analysis_system.data.00.gz
testsuite/analysis_system.data.01.gz
testsuite/analysis_system.data.02.gz
testsuite/analysis_system.data.03.gz
testsuite/analysis_system.data.04.gz
testsuite/analysis_system.data.05.gz
testsuite/analysis_system.data.06.gz
testsuite/analysis_system.data.07.gz
testsuite/analysis_system.data.08.gz
testsuite/analysis_system.data.09.gz
testsuite/analysis_system.data.10.gz
testsuite/analysis_system.data.chk
testsuite/angle_cosine.data
testsuite/angle_cossquare.data
testsuite/angle_harmonic.data
testsuite/comfixed_system.data
testsuite/comforce_system.data
testsuite/constraints_system.data
testsuite/dh_system.data
testsuite/el2d_system.data
testsuite/el2d_system_die.data
testsuite/gb_system.data
testsuite/lb_system.data
testsuite/lj-cos_system.data
testsuite/lj_system.data
testsuite/maggs_correct_rdf.data
testsuite/mass_system.data
testsuite/mdlc_expected_energy.data
testsuite/mdlc_expected_force_torque.data
testsuite/mdlc_system.data
testsuite/mmm1d_system.data
testsuite/npt_lj_system.data
testsuite/object_in_fluid_system-final.data
testsuite/object_in_fluid_system-init.data
testsuite/object_in_fluid_system-nodes.data
testsuite/object_in_fluid_system-triangles.data
testsuite/p3m_magnetostatics2_expected.data
testsuite/p3m_magnetostatics2_system.data
testsuite/p3m_magnetostatics.data
testsuite/p3m_system.data
testsuite/p3m_wall_system.data
testsuite/pe_micelle.data
testsuite/tabulated_system.data
testsuite/thermostat.data
testsuite/thermostat_rot.data
testsuite/uwerr.data
```

> "Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
> Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
>  Max-Planck-Institute for Polymer Research, Theory Group
>  
> This file is part of ESPResSo.
>  
> ESPResSo is free software: you can redistribute it and/or modify it
> under the terms of the GNU General Public License as published by the
> Free Software Foundation, either version 3 of the License, or (at your
> option) any later version.
>  
> ESPResSo is distributed in the hope that it will be useful, but
> WITHOUT ANY WARRANTY; without even the implied warranty of
> MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
> General Public License for more details.
>  
> You should have received a copy of the GNU General Public License
> along with this program.  If not, see <http://www.gnu.org/licenses/>."

This file is licensed as follows:

> Copyright (C) 2010,2011,2012,2013,2014 The ESPResSo project
> 
> Copying and distribution of this file, with or without modification,
> are permitted in any medium without royalty provided the copyright
> notice and this notice are preserved.  This file is offered as-is,
>without any warranty.
