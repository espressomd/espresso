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
   \nonumber \left(\frac{\partial \vec{v}_{\mathrm{fl}}}{\partial t} + \vec{v}_{\mathrm{fl}} \cdot \vec{\nabla} \vec{v}_{\mathrm{fl}} \right) \rho_\mathrm{fl} & = & -{k_\mathrm{B}T}\, \nabla \rho_\mathrm{fl} - q_k n_k \nabla \Phi \\
   \label{eq:ek-model-velocity} & & +\, \eta \vec{\Delta} \vec{v}_{\mathrm{fl}} + (\eta / 3 + \eta_{\text{b}}) \nabla (\nabla \cdot \vec{v}_{\mathrm{fl}}) \vphantom{\left(\frac{\partial}{\partial}\right)} ; \\
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

.. _Setup:

Setup
-----

.. _Initialization:

Initialization
~~~~~~~~~~~~~~
::

    import espressomd
    import espressomd.electrokinetics
    system = espressomd.System(box_l=[10.0, 10.0, 10.0])
    system.time_step = 0.0
    system.cell_system.skin = 0.4
    ek = espressomd.electrokinetics.Electrokinetics(agrid=1.0, lb_density=1.0,
        viscosity=1.0, ext_force_density = [1,0,0], friction=1.0, T=1.0, prefactor=1.0,
        stencil='linkcentered', advection=True, fluid_coupling='friction')
    system.actors.add(ek)

.. note::

    Requires external feature ``LB_WALBERLA``, enabled with the CMake option
    ``-D WITH_WALBERLA=ON``.

The above is a minimal example how to initialize the LB fluid, and
it is very similar to the lattice-Boltzmann command in set-up. We
therefore refer the reader to chapter :ref:`Lattice-Boltzmann` for details on the
implementation of LB in |es| and describe only the major differences here.

To set up a proper LB fluid using this command one has to specify at
least the following options: ``agrid``, ``lb_density``, ``viscosity``,
``friction``, ``T``, and ``prefactor``. The other options can be
used to modify the behavior of the LB fluid. Note that the command does
not allow the user to set the time step parameter as is the case for the
lattice-Boltzmann command, this parameter is instead taken directly from the value set for
:attr:`espressomd.system.System.time_step`. The LB *mass density* is set independently from the
electrokinetic *number densities*, since the LB fluid serves only as a
medium through which hydrodynamic interactions are propagated, as will
be explained further in the next paragraph. If no ``lb_density`` is specified, then our
algorithm assumes ``lb_density= 1.0``. The two 'new' parameters are the temperature ``T`` at
which the diffusive species are simulated and the ``prefactor``
associated with the electrostatic properties of the medium. See the
above description of the electrokinetic equations for an explanation of
the introduction of a temperature, which does not come in directly via a
thermostat that produces thermal fluctuations.

``advection`` can be set to ``True`` or ``False``. It controls whether there should be an
advective contribution to the diffusive species' fluxes. Default is
``True``.

``fluid_coupling`` can be set to ``"friction"`` or ``"estatics"``. This option determines the force
term acting on the fluid. The former specifies the force term to be the
sum of the species fluxes divided by their respective mobilities while
the latter simply uses the electrostatic force density acting on all
species. Note that this switching is only possible for the ``"linkcentered"``
stencil. For all other stencils, this choice is hardcoded. The default
is ``"friction"``.

``es_coupling`` enables the action of the electrostatic potential due to the
electrokinetics species and charged boundaries on the MD particles. The
forces on the particles are calculated by interpolation from the
electric field which is in turn calculated from the potential via finite
differences. This only includes interactions between the species and
boundaries and MD particles, not between MD particles and MD particles.
To get complete electrostatic interactions a particles Coulomb method
like Ewald or P3M has to be activated too.

The fluctuation of the EK species can be turned on by the flag ``fluctuations``.
This adds a white-noise term to the fluxes. The amplitude of this noise term
can be controlled by ``fluctuation_amplitude``. To circumvent that these fluctuations
lead to negative densities, they are modified by a smoothed Heaviside function,
which decreases the magnitude of the fluctuation for densities close to 0.
By default the fluctuations are turned off.

..
    .. _Diffusive species:

    Diffusive species
    ~~~~~~~~~~~~~~~~~
    ::

        species = electrokinetics.Species(density=density, D=D, valency=valency,
            ext_force_density=ext_force)

    :class:`espressomd.electrokinetics.Species` is used to initialize a diffusive species. Here the
    options specify: the number density ``density``, the diffusion coefficient ``D``, the
    valency of the particles of that species ``valency``, and an optional external
    (electric) force which is applied to the diffusive species. As mentioned
    before, the LB density is completely decoupled from the electrokinetic
    densities. This has the advantage that greater freedom can be achieved
    in matching the internal parameters to an experimental system. Moreover,
    it is possible to choose parameters for which the LB is more stable.
    The species can be added to a LB fluid::

        ek.add_species(species)

    One can also add the species during the initialization step of the
    :class:`~espressomd.electrokinetics.Electrokinetics` class by defining
    the list variable ``species``::

        ek = espressomd.electrokinetics.Electrokinetics(species=[species], ...)

    The variables ``density``, ``D``, and
    ``valency`` must be set to properly initialize the diffusive species; the
    ``ext_force_density`` is optional.

    .. _EK boundaries:

    EK boundaries
    ~~~~~~~~~~~~~
    ::

        ek_boundary = espressomd.ekboundaries.EKBoundary(charge_density=1.0, shape=my_shape)
        system.ekboundaries.add(ek_boundary)

    .. note:: Feature ``EK_BOUNDARIES`` required

    The :class:`~espressomd.ekboundaries.EKBoundary` class is used to set up
    internal or external) boundaries for the electrokinetics algorithm in much
    the same way as the :class:`~espressomd.lbboundaries.LBBoundary` class is
    used for the LB fluid. The major difference with the LB class is the option
    ``charge_density``, with which a boundary can be endowed with a volume
    charge density. To create a surface charge density, a combination of two
    oppositely charged boundaries, one inside the other, can be used. However,
    care should be taken to maintain the surface charge density when the value of ``agrid``
    is changed. Examples for possible shapes are wall, sphere, ellipsoid, cylinder,
    rhomboid and hollow conical frustum. We refer to the documentation of the
    :class:`espressomd.shapes` module for more possible shapes and information on
    the options associated to these shapes. In order to properly set up the
    boundaries, the ``charge_density`` and ``shape`` must be specified.

.. _Checkpointing EK:

Checkpointing
~~~~~~~~~~~~~
::

    ek.save_checkpoint(path)

Checkpointing in the EK works quite similar to checkpointing in the LB,
because the density is not saved within the :class:`espressomd.checkpointing`
object. However one should keep in mind, that the EK not only saves the density
of the species but also saves the population of the LB fluid in a separate file.
To load a checkpoint the ``espressomd.electrokinetics.Electrokinetics``
should have the same name as in the script it was saved, but to use the species
one need to extract them from the ``espressomd.electrokinetics.Electrokinetics``
via ``species``::

    checkpoint.load(cpt_path)
    species = ek.get_params()['species']
    ek.load_checkpoint(path)

..
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

    A property of the fluid field can be exported into a
    file in one go. Currently supported
    are: density, velocity, potential and boundary, which give the LB fluid density, the LB fluid velocity,
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

    .. _Local Quantities:

    Local Quantities
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
