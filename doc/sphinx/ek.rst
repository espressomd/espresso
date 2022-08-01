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
~~~~~~~~~~~~~~

Here is a minimal working example::

    import espressomd
    import espressomd.lb
    import espressomd.EKSpecies

    system = espressomd.System(box_l=3 * [6.0])
    system.time_step = 0.01
    system.cell_system.skin = 1.0

    ek_lattice = espressomd.lb.LatticeWalberla(agrid=0.5, n_ghost_layers=1)
    ek_solver = espressomd.EKSpecies.EKNone(lattice=ek_lattice)
    system.ekcontainer.tau = system.time_step
    system.ekcontainer.solver = ek_solver

.. note::

    Requires external feature ``LB_WALBERLA``, enabled with the CMake option
    ``-D WITH_WALBERLA=ON``.

An EK system can be set up at the same time as a LB system. The EK ``density``
represents the electrokinetic *number densities* and is independent of the
LB ``density``. The thermal energy ``kT`` controls thermal fluctuations,
``friction_coupling`` controls coupling of the diffusive species with the
LB fluid force, ``advection`` controls whether there should be an advective
contribution to the diffusive species' fluxes from the LB fluid.

.. _Diffusive species:

Diffusive species
~~~~~~~~~~~~~~~~~
::

    ek_species = espressomd.EKSpecies.EKSpecies(
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

:class:`~espressomd.EKSpecies.EKSpecies` is used to initialize a diffusive
species. Here the options specify: the number density ``density``,
the diffusion coefficient ``diffusion``, the valency of the particles
of that species ``valency``, and an optional external (electric) force
``ext_efield`` which is applied to the diffusive species. As mentioned
before, the LB density is completely decoupled from the electrokinetic
densities. This has the advantage that greater freedom can be achieved
in matching the internal parameters to an experimental system. Moreover,
it is possible to choose parameters for which the LB is more stable.

To add species to the EK solver::

    system.ekcontainer.add(ek_species)

To remove species from the EK solver::

    system.ekcontainer.remove(ekspecies)

..  .. _EK boundaries:

    EK boundaries
    ~~~~~~~~~~~~~

    :class:`~espressomd.ekboundaries.EKBoundary` is used to set up
    internal (or external) boundaries for the electrokinetics algorithm in much
    the same way as the :class:`~espressomd.lbboundaries.LBBoundary` class is
    used for the LB fluid::

        ek_boundary = espressomd.ekboundaries.EKBoundary(charge_density=1.0, shape=my_shape)
        system.ekboundaries.add(ek_boundary)

    .. note:: Feature ``EK_BOUNDARIES`` required

    The major difference with the LB class is the option ``charge_density``,
    with which a boundary can be endowed with a volume charge density.
    To create a surface charge density, a combination of two
    oppositely charged boundaries, one inside the other, can be used. However,
    care should be taken to maintain the surface charge density when the value of ``agrid``
    is changed. Examples for possible shapes are wall, sphere, ellipsoid, cylinder,
    rhomboid and hollow conical frustum. We refer to the documentation of the
    :class:`espressomd.shapes` module for more possible shapes and information on
    the options associated to these shapes. In order to properly set up the
    boundaries, the ``charge_density`` and ``shape`` must be specified.

    .. _Output:

    Output
    ~~~~~~

    .. _Fields:

    Fields
    """"""

    ::

        ek.write_vtk_boundary(path)
        ek.write_vtk_density(path)
        ek.write_vtk_velocity(path)
        ek.write_vtk_potential(path)

    A property of the fluid field can be exported into a file in one go.
    Currently supported fields are: density, velocity, potential and boundary,
    which give the LB fluid density, the LB fluid velocity,
    the electrostatic potential, and the location and type of the
    boundaries, respectively. The boundaries can only be printed when the
    ``EK_BOUNDARIES`` is compiled in. The output is a vtk-file, which is readable by
    visualization software such as ParaView [5]_ and Mayavi2 [6]_.

    ::

        species.write_vtk_flux(path)
        species.write_vtk_density(path)

    These commands are similar to the above. They enable the
    export of diffusive species properties, namely: ``density`` and ``flux``, which specify the
    number density and flux of species ``species``, respectively.

    .. _Local quantities:

    Local quantities
    """"""""""""""""

    Local quantities like velocity or fluid density for single nodes can be accessed in the same way
    as for an LB fluid, see :ref:`Lattice-Boltzmann`. The only EK-specific quantity is the potential.

    ::

        ek[0, 0, 0].potential
        ek[0, 0, 0].velocity
        ek[0, 0, 0].boundary

    The local ``density`` and ``flux`` of a species can be obtained in the same fashion:

    ::

        species[0, 0, 0].density
        species[0, 0, 0].flux

    .. [5]
       https://www.paraview.org/
    .. [6]
       http://code.enthought.com/projects/mayavi/
