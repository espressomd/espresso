.. _Electrokinetics:

Electrokinetics
===============

The electrokinetics setup in |es| allows for the description of
electro-hydrodynamic systems on the level of ion density distributions
coupled to a lattice-Boltzmann (LB) fluid. The ion density distributions
may also interact with explicit charged particles, which are
interpolated on the LB grid. In the following paragraph we briefly
explain the electrokinetic model implemented in |es|, before we come to the
description of the interface.

.. note::
    Please cite :cite:t:`godenschwager13a` and :cite:t:`bauer21a` (BibTeX keys
    ``godenschwager13a`` and ``bauer21a`` in :file:`doc/bibliography.bib`) if
    you use the LB fluid. When generating your own kernels with pystencils and
    lbmpy, please also cite :cite:t:`bauer19a` and :cite:t:`bauer21b` (BibTeX
    key ``bauer19a`` resp. ``bauer21b`` in :file:`doc/bibliography.bib`).

.. note::

    Requires external features ``WALBERLA`` and optionally ``WALBERLA_FFT``
    (for the FFT-based Poisson solver), enabled with the CMake options
    ``-D ESPRESSO_BUILD_WITH_WALBERLA=ON -D ESPRESSO_BUILD_WITH_WALBERLA_FFT=ON``.

.. _Electrokinetic equations:

Electrokinetic equations
------------------------

In the electrokinetics code we solve the following system of coupled
continuity, diffusion-advection, Poisson, and Navier-Stokes equations:

.. math::

   \begin{aligned}
   \label{eq:ek-model-continuity} \frac{\partial n_k}{\partial t} & = & -\, \nabla \cdot \vec{j}_k \vphantom{\left(\frac{\partial}{\partial}\right)} ; \\
   \label{eq:ek-model-fluxes} \vec{j}_{k} & = & -D_k \nabla n_k - \nu_k \, q_k n_k\, \nabla \Phi + n_k \vec{v}_{\mathrm{fl}} \vphantom{\left(\frac{\partial}{\partial}\right)} + \sqrt{n_k}\vec{\mathcal{W}}_k; \\
   \label{eq:ek-model-poisson} \Delta \Phi & = & -4 \pi \, {l_\mathrm{B}}\, {k_\mathrm{B}T}\sum_k q_k n_k \vphantom{\left(\frac{\partial}{\partial}\right)}; \\
   \nonumber \left(\frac{\partial \vec{v}_{\mathrm{fl}}}{\partial t} + \vec{v}_{\mathrm{fl}} \cdot \nabla \vec{v}_{\mathrm{fl}} \right) \rho_\mathrm{fl} & = & -{k_\mathrm{B}T}\, \nabla \rho_\mathrm{fl} - q_k n_k \nabla \Phi \\
   \label{eq:ek-model-velocity} & & +\, \eta \Delta \vec{v}_{\mathrm{fl}} + (\eta / 3 + \eta_{\text{b}}) \nabla (\nabla \cdot \vec{v}_{\mathrm{fl}}) \vphantom{\left(\frac{\partial}{\partial}\right)} ; \\
   \label{eq:ek-model-continuity-fl} \frac{\partial \rho_\mathrm{fl}}{\partial t} & = & -\,\nabla\cdot\left( \rho_\mathrm{fl} \vec{v}_{\mathrm{fl}} \right) \vphantom{\left(\frac{\partial}{\partial}\right)} , \end{aligned}

which define relations between the following observables

:math:`n_k`
    the number density of the particles of species :math:`k`,

:math:`\vec{j}_k`
    the number density flux of the particles of species :math:`k`,

:math:`\Phi`
    the electrostatic potential,

:math:`\rho_{\mathrm{fl}}`
    the mass density of the fluid,

:math:`\vec{v}_{\mathrm{fl}}`
    the advective velocity of the fluid,

and input parameters

:math:`D_k`
    the diffusion constant of species :math:`k`,

:math:`\nu_k`
    the mobility of species :math:`k`,

:math:`\vec{\mathcal{W}}_k`
    the white-noise term for the fluctuations of species :math:`k`,

:math:`q_k`
    the charge of a single particle of species :math:`k`,

:math:`{l_\mathrm{B}}`
    the Bjerrum length,

:math:`{k_\mathrm{B}T}`
    | the thermal energy given by the product of Boltzmann's constant
      :math:`k_\text{B}`
    | and the temperature :math:`T`,

:math:`\eta`
    the dynamic viscosity of the fluid,

:math:`\eta_{\text{b}}`
    the bulk viscosity of the fluid.

The temperature :math:`T`, and diffusion constants :math:`D_k` and
mobilities :math:`\nu_k` of individual species are linked through the
Einstein-Smoluchowski relation :math:`D_k /
\nu_k = {k_\mathrm{B}T}`. This system of equations
combining diffusion-advection, electrostatics, and hydrodynamics is
conventionally referred to as the *Electrokinetic Equations*.

The electrokinetic equations have the following properties:

-  On the coarse time and length scale of the model, the dynamics of the
   particle species can be described in terms of smooth density
   distributions and potentials as opposed to the microscale where
   highly localized densities cause singularities in the potential.

   In most situations, this restricts the application of the model to
   species of monovalent ions, since ions of higher valency typically
   show strong condensation and correlation effects – the localization
   of individual ions in local potential minima and the subsequent
   correlated motion with the charges causing this minima.

-  Only the entropy of an ideal gas and electrostatic interactions are
   accounted for. In particular, there is no excluded volume.

   This restricts the application of the model to monovalent ions and
   moderate charge densities. At higher valencies or densities,
   overcharging and layering effects can occur, which lead to
   non-monotonic charge densities and potentials, that can not be
   covered by a mean-field model such as Poisson--Boltzmann or this one.

   Even in salt free systems containing only counter ions, the
   counter-ion densities close to highly charged objects can be
   overestimated when neglecting excluded volume effects. Decades of the
   application of Poisson--Boltzmann theory to systems of electrolytic
   solutions, however, show that those conditions are fulfilled for
   monovalent salt ions (such as sodium chloride or potassium chloride)
   at experimentally realizable concentrations.

-  Electrodynamic and magnetic effects play no role. Electrolytic
   solutions fulfill those conditions as long as they don't contain
   magnetic particles.

-  The diffusion coefficient is a scalar, which means there can not be
   any cross-diffusion. Additionally, the diffusive behavior has been
   deduced using a formalism relying on the notion of a local
   equilibrium. The resulting diffusion equation, however, is known to
   be valid also far from equilibrium.

-  The temperature is constant throughout the system.

-  The density fluxes instantaneously relax to their local equilibrium
   values. Obviously one can not extract information about processes on
   length and time scales not covered by the model, such as dielectric
   spectra at frequencies, high enough that they correspond to times
   faster than the diffusive time scales of the charged species.

.. _EK Setup:

Setup
-----

.. _EK Initialization:

Initialization
^^^^^^^^^^^^^^

Here is a minimal working example::

    import espressomd
    import espressomd.electrokinetics

    system = espressomd.System(box_l=3 * [6.0])
    system.time_step = 0.01
    system.cell_system.skin = 1.0

    ek_lattice = espressomd.electrokinetics.LatticeWalberla(agrid=0.5, n_ghost_layers=1)
    ek_solver = espressomd.electrokinetics.EKNone(lattice=ek_lattice)
    system.ekcontainer.solver = ek_solver
    system.ekcontainer.tau = system.time_step

where ``system.ekcontainer`` is the EK system, ``ek_solver`` is the Poisson
solver (here ``EKNone`` doesn't actually solve the electrostatic field, but
instead imposes a zero field), and ``ek_lattice`` contains the grid parameters.
In this setup, the EK system doesn't contain any species. The following
sections will show how to add species that can diffuse, advect, react and/or
electrostatically interact. An EK system can be set up at the same time as a
LB system.

.. _Diffusive species:

Diffusive species
^^^^^^^^^^^^^^^^^
::

    ek_species = espressomd.electrokinetics.EKSpecies(
        lattice=ek_lattice,
        single_precision=False,
        kT=1.0,
        density=0.85,
        valency=0.0,
        diffusion=0.1,
        advection=False,
        friction_coupling=False,
        ext_efield=[0., 0., 0.]
    )

:class:`~espressomd.electrokinetics.EKSpecies` is used to initialize a diffusive
species. Here the options specify: the electrokinetic *number densities*
``density`` (independent from the LB ``density``), the diffusion coefficient
``diffusion``, the valency of the particles of that species ``valency``,
the optional external (electric) force ``ext_efield`` which is applied to
the diffusive species, the thermal energy ``kT`` for thermal fluctuations,
``friction_coupling`` to enable coupling of the diffusive species to the
LB fluid force and ``advection`` to add an advective contribution to the
diffusive species' fluxes from the LB fluid.
Multiple species can be added to the EK system.

To add species to the EK system::

    system.ekcontainer.add(ek_species)

To remove species from the EK system::

    system.ekcontainer.remove(ek_species)

Individual nodes and slices of the species lattice can be accessed and
modified using the syntax outlined in :ref:`Reading and setting properties
of single lattice nodes`.

As mentioned before, the LB density is completely decoupled from the
electrokinetic densities. This has the advantage that greater freedom can
be achieved in matching the internal parameters to an experimental system.
Moreover, it is possible to choose parameters for which the LB is more stable.

Performance considerations
^^^^^^^^^^^^^^^^^^^^^^^^^^

The CPU implementation of the EK has an extra flag ``single_precision`` to
use single-precision floating point values. These are approximately 10%
faster than double-precision, at the cost of a small loss in precision.

.. _Checkpointing EK:

Checkpointing
-------------

::

    ek.save_checkpoint(path, binary)
    ek.load_checkpoint(path, binary)

The first command saves all of the EK nodes' properties to an ASCII
(``binary=False``) or binary (``binary=True``) format respectively.
The second command loads the EK nodes' properties.
In both cases ``path`` specifies the location of the
checkpoint file. This is useful for restarting a simulation either on the same
machine or a different machine. Some care should be taken when using the binary
format as the format of doubles can depend on both the computer being used as
well as the compiler.

.. _EK VTK output:

VTK output
----------

The waLBerla library implements a globally-accessible VTK registry.
A VTK stream can be attached to an EK actor to periodically write
one or multiple fluid field data into a single file using
:class:`~espressomd.electrokinetics.VTKOutput`::

    vtk_obs = ["density"]
    # create a VTK callback that automatically writes every 10 EK steps
    ek_vtk = espressomd.electrokinetics.VTKOutput(
        identifier="ek_vtk_automatic", observables=vtk_obs, delta_N=10)
    ek.add_vtk_writer(vtk=ek_vtk)
    system.integrator.run(100)
    # can be deactivated
    ek_vtk.disable()
    system.integrator.run(10)
    ek_vtk.enable()
    # create a VTK callback that writes only when explicitly called
    ek_vtk_on_demand = espressomd.electrokinetics.VTKOutput(
        identifier="ek_vtk_now", observables=vtk_obs)
    ek.add_vtk_writer(vtk=ek_vtk_on_demand)
    ek_vtk_on_demand.write()

Currently only supports the species density.
By default, the properties of the current state
of the species are written to disk on demand. To add a stream that writes
to disk continuously, use the optional argument ``delta_N`` to indicate
the level of subsampling. Such a stream can be deactivated.

The VTK format is readable by visualization software such as ParaView [5]_
or Mayavi2 [6]_, as well as in |es| (see :ref:`Reading VTK files`).
If you plan to use ParaView for visualization, note that also the particle
positions can be exported using the VTK format
(see :meth:`~espressomd.particle_data.ParticleList.writevtk`).

Important: these VTK files are written in multi-piece format, i.e. each MPI
rank writes its local domain to a new piece in the VTK uniform grid to avoid
a MPI reduction. ParaView can handle the topology reconstruction natively.
However, when reading the multi-piece file with the Python ``vtk`` package,
the topology must be manually reconstructed. In particular, calling the XML
reader ``GetOutput()`` method directly after the update step will erase all
topology information. While this is not an issue for VTK files obtained from
simulations that ran with 1 MPI rank, for parallel simulations this will lead
to 3D grids with incorrectly ordered data. Automatic topology reconstruction
is available through :class:`~espressomd.io.vtk.VTKReader`::

    import pathlib
    import tempfile
    import numpy as np
    import espressomd
    import espressomd.electrokinetics
    import espressomd.io.vtk

    system = espressomd.System(box_l=[12., 14., 10.])
    system.cell_system.skin = 0.4
    system.time_step = 0.1

    lattice = espressomd.electrokinetics.LatticeWalberla(agrid=1.)
    species = espressomd.electrokinetics.EKSpecies(
            lattice=lattice, density=1., kT=1., diffusion=0.1, valency=0.,
            advection=False, friction_coupling=False, tau=system.time_step)
    system.ekcontainer.tau = species.tau
    system.ekcontainer.add(species)
    system.integrator.run(10)

    vtk_reader = espressomd.io.vtk.VTKReader()
    label_density = "density"

    with tempfile.TemporaryDirectory() as tmp_directory:
        path_vtk_root = pathlib.Path(tmp_directory)
        label_vtk = "ek_vtk"
        path_vtk = path_vtk_root / label_vtk / "simulation_step_0.vtu"

        # write VTK file
        ek_vtk = espressomd.electrokinetics.VTKOutput(
            identifier=label_vtk, delta_N=0,
            observables=["density"],
            base_folder=str(path_vtk_root))
        species.add_vtk_writer(vtk=ek_vtk)
        ek_vtk.write()

        # read VTK file
        vtk_grids = vtk_reader.parse(path_vtk)
        vtk_density = vtk_grids[label_density]

        # check VTK values match node values
        ek_density = np.copy(lbf[:, :, :].density)
        np.testing.assert_allclose(vtk_density, ek_density, rtol=1e-10, atol=0.)

.. _Setting up EK boundary conditions:

Setting up boundary conditions
------------------------------

It is possible to impose a fixed density and a fixed flux on EK species.

Under the hood, a boundary field is added to the blockforest, which contains
pre-calculated information for the streaming operations.

.. _Per-node EK boundary conditions:

Per-node boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

One can set (or update) the boundary conditions of individual nodes::

    import espressomd
    import espressomd.electrokinetics
    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    system.cell_system.skin = 0.1
    system.time_step = 0.01
    lattice = espressomd.electrokinetics.LatticeWalberla(agrid=0.5, n_ghost_layers=1)
    ek_species = espressomd.electrokinetics.EKSpecies(
        kT=1.5, lattice=self.lattice, density=0.85, valency=0., diffusion=0.1,
        advection=False, friction_coupling=False, tau=system.time_step)
    system.ekcontainer.tau = species.tau
    system.ekcontainer.add(ek_species)
    # set node fixed density boundary conditions
    lbf[0, 0, 0].boundary = espressomd.electrokinetics.DensityBoundary(1.)
    # update node fixed density boundary conditions
    lbf[0, 0, 0].boundary = espressomd.electrokinetics.DensityBoundary(2.)
    # remove node boundary conditions
    lbf[0, 0, 0].boundary = None

.. _Shape-based EK boundary conditions:

Shape-based boundary conditions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Adding a shape-based boundary is straightforward::

    import espressomd
    import espressomd.electrokinetics
    import espressomd.shapes
    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    system.cell_system.skin = 0.1
    system.time_step = 0.01
    lattice = espressomd.electrokinetics.LatticeWalberla(agrid=0.5, n_ghost_layers=1)
    ek_species = espressomd.electrokinetics.EKSpecies(
        kT=1.5, lattice=self.lattice, density=0.85, valency=0.0, diffusion=0.1,
        advection=False, friction_coupling=False, tau=system.time_step)
    system.ekcontainer.tau = species.tau
    system.ekcontainer.add(ek_species)
    # set fixed density boundary conditions
    wall = espressomd.shapes.Wall(normal=[1., 0., 0.], dist=2.5)
    ek_species.add_boundary_from_shape(
        shape=wall, value=1., boundary_type=espressomd.electrokinetics.DensityBoundary)
    # clear fixed density boundary conditions
    ek_species.clear_density_boundaries()

For a position-dependent flux, the argument to ``value`` must be a 4D grid
(the first three dimensions must match the EK grid shape, the fourth
dimension has size 3 for the flux).

For a complete description of all available shapes, refer to
:mod:`espressomd.shapes`.

.. _Prototyping new EK methods:

Prototyping new EK methods
--------------------------

Start by installing the code generator dependencies:

.. code-block:: bash

    python3 -m pip install --user -c requirements.txt numpy sympy lbmpy pystencils islpy

Next, edit the code generator script to configure new kernels, then execute it:

.. code-block:: bash

    python3 maintainer/walberla_kernels/generate_lb_kernels.py

The script takes optional arguments to control the CPU or GPU architecture,
as well as the floating-point precision. The generated source code files need
to be written to :file:`src/walberla_bridge/src/electrokinetics/generated_kernels/`
and :file:`src/walberla_bridge/src/electrokinetics/reactions/generated_kernels/`.
These steps can be automated with the convenience shell functions documented in
:file:`maintainer/walberla_kernels/Readme.md`.
Edit the :file:`CMakeLists.txt` file in the destination folders to include the
new kernels in the build system.
Then, adapt :file:`src/walberla_bridge/src/electrokinetics/EKinWalberlaImpl.hpp`
to use the new EK kernels.


.. [5]
   https://www.paraview.org/
.. [6]
   http://code.enthought.com/projects/mayavi/
