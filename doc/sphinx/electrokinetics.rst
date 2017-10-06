Electrokinetics
===============

The electrokinetics setup in allows for the description of
electro-hydrodynamic systems on the level of ion density distributions
coupled to a Lattice-Boltzmann (LB) fluid. The ion density distributions
may also interact with explicit charged particles, which are
interpolated on the LB grid. In the following paragraph we briefly
explain the electrokinetic model implemented in , before we come to the
description of the interface.

If you are interested in using the electrokinetic implementation in for
scientific purposes, please contact G. Rempfer before you start your
project.

Electrokinetic Equations
------------------------

In the electrokinetics code we solve the following system of coupled
continuity, diffusion-advection, Poisson, and Navier-Stokes equations:

.. math::

   \begin{aligned}
   \label{eq:ek-model-continuity} \frac{\partial n_k}{\partial t} & = & -\, \nabla \cdot \vec{j}_k \vphantom{\left(\frac{\partial}{\partial}\right)} ; \\
   \label{eq:ek-model-fluxes} \vec{j}_{k} & = & -D_k \nabla n_k - \nu_k \, q_k n_k\, \nabla \Phi + n_k \vec{v}_{\mathrm{fl}} \vphantom{\left(\frac{\partial}{\partial}\right)} ; \\
   \label{eq:ek-model-poisson} \Delta \Phi & = & -4 \pi \, {l_\mathrm{B}}\, {k_\mathrm{B}T}\sum_k q_k n_k \vphantom{\left(\frac{\partial}{\partial}\right)}; \\
   \nonumber \left(\frac{\partial \vec{v}_{\mathrm{fl}}}{\partial t} + \vec{v}_{\mathrm{fl}} \cdot \vec{\nabla} \vec{v}_{\mathrm{fl}} \right) \rho_\mathrm{fl} & = & -{k_\mathrm{B}T}\, \nabla \rho_\mathrm{fl} - q_k n_k \nabla \Phi \\
   \label{eq:ek-model-velocity} & & +\, \eta \vec{\Delta} \vec{v}_{\mathrm{fl}} + (\eta / 3 + \eta_{\text{\,b}}) \nabla (\nabla \cdot \vec{v}_{\mathrm{fl}}) \vphantom{\left(\frac{\partial}{\partial}\right)} ; \\
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

:math:`q_k`
    the charge of a single particle of species :math:`k`,

:math:`{l_\mathrm{B}}`
    the Bjerrum length,

:math:`{k_\mathrm{B}T}`
    | the thermal energy given by the product of Boltzmann’s constant
      :math:`k_\text{B}`
    | and the temperature :math:`T`,

:math:`\eta`
    the dynamic viscosity of the fluid,

:math:`\eta_{\text{\,b}}`
    the bulk viscosity of the fluid.

The temperature :math:`T`, and diffusion constants :math:`D_k` and
mobilities :math:`\nu_k` of individual species are linked through the
Einstein-Smoluchowski relation :math:`D_k /
\nu_k = {k_\mathrm{B}T}`. The system of equations described in Eqs. -,
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
   covered by a mean-field model such as Poisson-Boltzmann or this one.

   Even in salt free systems containing only counter ions, the
   counter-ion densities close to highly charged objects can be
   overestimated when neglecting excluded volume effects. Decades of the
   application of Poisson-Boltzmann theory to systems of electrolytic
   solutions, however, show that those conditions are fulfilled for
   monovalent salt ions (such as sodium chloride or potassium chloride)
   at experimentally realizable concentrations.

-  Electrodynamic and magnetic effects play no role. Electrolytic
   solutions fulfill those conditions as long as they don’t contain
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

Setup
-----

[ssec:ek-init]Initialization
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The command initializes the LB fluid with a given set of parameters, and
it is very similar to the Lattice-Boltzmann command in set-up. We
therefore refer the reader to Chapter [sec:lb] for details on the
implementation of LB in and describe only the major differences here.

The first major difference with the LB implementation is that the
electrokinetics set-up is a Graphics Processing Unit (GPU) only
implementation. There is no Central Processing Unit (CPU) version, and
at this time there are no plans to make a CPU version available in the
future. To use the electrokinetics features it is therefore imperative
that your computer contains a CUDA capable GPU which is sufficiently
modern.

To set up a proper LB fluid using the command one has to specify at
least the following options: , , , , , and . The other options can be
used to modify the behavior of the LB fluid. Note that the command does
not allow the user to set the time step parameter as is the case for the
command, this parameter is instead taken directly from the input of the
``t_step`` command. The LB *mass density* is set independently from the
electrokinetic *number densities*, since the LB fluid serves only as a
medium through which hydrodynamic interactions are propagated, as will
be explained further in the next paragraph. If no is specified, then our
algorithm assumes = 1.0. The two ‘new’ parameters are the temperature at
which the diffusive species are simulated and the Bjerrum length
associated with the electrostatic properties of the medium. See the
above description of the electrokinetic equations for an explanation of
the introduction of a temperature, which does not come in directly via a
thermostat that produces thermal fluctuations.

can be set to *on* or *off*. It controls whether there should be an
advective contribution to the diffusive species’ fluxes. Default is
*on*.

can be set to *friction* or *estatics*. This option determines the force
term acting on the fluid. The former specifies the force term to be the
sum of the species fluxes divided by their respective mobilities while
the latter simply uses the electrostatic force density acting on all
species. Note that this switching is only possible for the linkcentered
stencil. For all other stencils, this choice is hardcoded. The default
is *friction*.

enables the action of the electrostatic potential due to the
electrokinetics species and charged boundaries on the MD particles. The
forces on the particles are calculated by interpolation from the
electric field which is in turn calculated from the potential via finite
differences. This only includes interactions between the species and
boundaries and MD particles, not between MD particles and MD particles.
To get complete electrostatic interactions a particles Coulomb method
like Ewald or P3M has to be activate too.

[ssec:ek-diff-species]Diffusive Species
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The command followed by an integer (in the range 0 to 10) and several
options can be used to initialize the diffusive species. Here the
options specify: the number density , the diffusion coefficient , the
valency of the particles of that species , and an optional external
(electric) force which is applied to the diffusive species. As mentioned
before, the LB density is completely decoupled from the electrokinetic
densities. This has the advantage that greater freedom can be achieved
in matching the internal parameters to an experimental system. Moreover,
it is possible to choose parameters for which the LB is more stable. The
LB fluid must already be (partially) set up using the ... command,
before the diffusive species can be initialized. The variables , , and
must be set to properly initialize the diffusive species; the is
optional.

[ssec:ek-boundaries]Boundaries
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The command allows one to set up (internal or external) boundaries for
the electrokinetics algorithm in much the same way as the command is
used for the LB fluid. The major difference with the LB command is given
by the option , with which a boundary can be endowed with a volume
charge density. To create a surface charge density, a combination of two
oppositely charged boundaries, one inside the other, can be used.
However, care should be taken to maintain the surface charge density
when the value of is changed. Currently, the following s are available:
wall, sphere, cylinder, rhomboid, pore, stomatocyte, hollow\_cone, and
spherocylinder. We refer to the documentation of the command
(Chapter [sec:lb]) for information on the options associated to these
shapes. In order to properly set up the boundaries, the and relevant
must be specified.

[ssec:ek-output]Output
----------------------

[ssec:ek-output-fields]Fields
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print

The print parameter of the command enables simple visualization of
simulation data. A property of the fluid field can be exported into a
file with name in one go. Currently, supported values of the parameter
are: , , , and , which give the LB fluid density, the LB fluid velocity,
the electrostatic potential, and the location and type of the
boundaries, respectively. The boundaries can only be printed when the
``EK_BOUNDARIES`` is compiled in. The additional option can be used to
directly export in the vtk format. The vtk format is readable by
visualization software such as paraview [1]_ and mayavi2 [2]_. If the
option is not specified, a gnuplot readable data file will be exported.

print

This print statement is similar to the above command. It enables the
export of diffusive species properties, namely: and , which specify the
number density and flux of species , respectively.

[ssec:ek-local-quantities]Local Quantities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

node velocity

The option of the command allows one to output the value of a quantity
on a single LB node. The node is addressed using three integer values
which run from 0 to /, /, and /, respectively. Thus far, only the
velocity of the LB fluid can be printed in the standard electrokinetics
implementation. For other quantities the command may be used.

node density

This command can be used to output the number density of the -th
diffusive species on a single LB node.

[ssec:ek-checkpointing]Checkpointing
------------------------------------

checkpoint save checkpoint load

Variant writes the species density fields as well as all necessary LB
fields into two files called and . Variant reads these files back and
restores the fields. All algorithm parameters must be set via the
simulation script, as they are not part of the checkpointed data.

The format of the checkpoint is binary and no special care is taken with
respect to the specific binary layout of the machine.

.. _Catalytic Reaction:

Catalytic Reaction
------------------

Concept
~~~~~~~

The electrokinetics solver implemented in can be used to simulate a
system, for which in addition to the electrokinetic equations, there is
a (local) catalytic reaction which converts one species into another.

If you are interested in using this implementation in for scientific
purposes, please contact J. de Graaf before you start your project.

Currently, a linear reaction is implemented which converts one species
into two others, in order to model the catalytic decomposition of
hydrogen peroxide in the presence of a platinum catalyst:
:math:`2 \mathrm{H}_{2}\mathrm{O}_{2} \rightarrow 
2 \mathrm{H}_{2}\mathrm{O} + \mathrm{O}_{2}`. The decomposition of
:math:`\mathrm{H}_{2}\mathrm{O}_{2}` is in reality more complicated than
the linear reaction introduced here, since it is assumed to proceed via
several intermediate complexed-states, but our model can be thought of
as modeling the rate-limiting step. If we assume that there are three
non-ionic species with number densities :math:`n_{k}`, where
:math:`n_{0} = [ \mathrm{H}_{2}\mathrm{O}_{2} ]`,
:math:`n_{1} = [ \mathrm{H}_{2}\mathrm{O} ]`, and
:math:`n_{2} = [ \mathrm{O}_{2} ]`, then we can write the
(electro)kinetic equations for this system as

.. math::

   \begin{aligned}
   \label{eq:ek-reaction-continuity} \frac{\partial n_k}{\partial t} & = & -\, \nabla \cdot \vec{j}_k +\, f_{k} c n_{k} \vphantom{\left(\frac{\partial}{\partial}\right)} ; \\
   \label{eq:ek-reaction-fluxes} \vec{j}_{k} & = & -D_k \nabla n_k + n_k \vec{v}_{\mathrm{fl}} \vphantom{\left(\frac{\partial}{\partial}\right)} ; \\
   \nonumber \left(\frac{\partial \vec{v}_{\mathrm{fl}}}{\partial t} + \vec{v}_{\mathrm{fl}} \cdot \vec{\nabla} \vec{v}_{\mathrm{fl}} \right) \rho_\mathrm{fl} & = & -{k_\mathrm{B}T}\, \sum_{k} \nabla n_k   \\
   \label{eq:ek-reaction-velocity} & & +\, \eta \vec{\Delta} \vec{v}_{\mathrm{fl}} + (\eta / 3 + \eta_{\text{\,b}}) \nabla (\nabla \cdot \vec{v}_{\mathrm{fl}}) \vphantom{\left(\frac{\partial}{\partial}\right)} ; \\
   \label{eq:ek-reaction-continuity-fl} \frac{\partial \rho_\mathrm{fl}}{\partial t} & = & -\,\nabla\cdot\left( \rho_\mathrm{fl} \vec{v}_{\mathrm{fl}} \right) \vphantom{\left(\frac{\partial}{\partial}\right)} ,\end{aligned}

which define relations between the following observables

:math:`n_k`
    the number density of the particles of species :math:`k`,

:math:`\vec{j}_k`
    the number density flux of the particles of species :math:`k`,

:math:`\rho_{\mathrm{fl}}`
    the mass density of the fluid,

:math:`\vec{v}_{\mathrm{fl}}`
    the advective velocity of the fluid,

and input parameters

:math:`D_k`
    the diffusion constant of species :math:`k`,

:math:`{k_\mathrm{B}T}`
    | the thermal energy given by the product of Boltzmann’s constant
      :math:`k_\text{B}`
    | and the temperature :math:`T`,

:math:`\eta`
    the dynamic viscosity of the fluid,

:math:`\eta_{\text{\,b}}`
    the bulk viscosity of the fluid,

:math:`f_{k}`
    the reaction constant :math:`f_{0} \equiv -1`, :math:`f_{1} = 1` and
    :math:`f_{2} = 0.5` for the above reaction,

:math:`c`
    the reaction rate.

In this set of equations we have fully decoupled the number densities
and the fluid mass density. N.B. We have set the initial fluid mass
density is not necessarily equal to the sum of the initial species
number densities. This means that some care needs to be taken in the
interpretation of the results obtained using this feature. In
particular, the solution of the Navier-Stokes equation exclusively
models the momentum transport through the (multicomponent) fluid, while
the diffusive properties of the individual chemical species are handled
by Eqs.  and .

It is important to note that to ensure mass conservation the reaction
must satisfy:

.. math:: \label{eq:ek-mass-balance} \sum_{k} f_{k} m_{k} = 0 ,

where :math:`m_{k}` is the molecular mass of a reactive species.
Unfortunately, the current electrokinetic implementation does not
conserve mass flux locally. That is to say, the LB fluid is compressible
and the sum of the fluxes of the three species is not equal to zero in
the frame co-moving with the advective fluid velocity. It is therefore
debatable whether it is necessary to impose Eq. , since the EK algorithm
itself does not conserve mass density. However, we strived to be as
accurate as possible and in future versions of the EK algorithm the lack
of incompressiblity will be addressed.

The reaction is specified by the second term on the right-hand side of
Eq. . It is important to note that this term can be set locally, as
opposed to the other terms in the equation system Eqs. -, in our
implementation, as will become clear in the following. This has the
advantage that catalytic surfaces may be modeled.

[ssec:ek-reac-init]Initialization and Geometry Definition
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The command is used to set up the catalytic reaction between three
previously defined the diffusive species, of which the i identifiers are
given by , , and , respectively. In the 1:2 reaction, these fulfill the
role of the reactant and the two products, as indicated by the naming
convention. For each species a reservoir (number) density must be set,
given by the variables , , and , respectively. These reservoir densities
correspond to the initial number densities associated with the reactive
species. The reservoir densities, in tandem with reservoir nodes, see
below, can be used to keep the reaction from depleting all the reactant
in the simulation box. The variable specifies the speed at which the
reaction proceeds. The three masses (typically given in the atomic
weight equivalent) are used to determine the total mass flux provided by
the reaction, as described above, and are also used to check whether the
reaction ratios that are given satisfy the chemical requirement of mass
conservation. Finally, the parameters and specify what fractions of the
product are generated when a given quantity of reactant is catalytically
converted. To use a chemical reaction, all options for the command must
be specified.

The option of the command allows one to set up regions in which the
reaction takes place with the help of the constraints that are available
to set up boundaries. The integer value can be used to select the
reaction: 0 no reaction takes place for this region, 1 the catalytic
reaction takes place in this region, and 2 the region functions as a
reservoir, wherein the species densities are reset to their initial (or
reservoir) concentrations. The rest of the command follows the same
format of the command. Currently, the following s are available: box,
wall, sphere, cylinder, rhomboid, pore, stomatocyte, hollow\_cone, and
spherocylinder. The box shape is a specific command, which can be used
to set the entire simulation box to a specific reaction value. To use
the command, one must first set up a reaction, as described above. To
successfully specify a region all the relevant arguments that go with
the shape constraints must be provided.

[sssec:ek-pdb-parse]Parsing PDB Files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The feature allows the user to parse simple PDB files, a file format
introduced by the protein database to encode molecular structures.
Together with a topology file (here ) the structure gets interpolated to
the grid. For the input you will need to prepare a PDB file with a force
field to generate the topology file. Normally the PDB file extension is
, the topology file extension is . Obviously the PDB file is placed
instead of and the topology file instead of .

[ssec:ek-reac-output]Reaction-Specific Output
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print

The print parameter of the command can be used in combination with the
``EK_REACTION`` feature to give advanced output options. Currently,
supported values of the parameter are: and , which give the location and
type of the reactive regions and the ideal-gas pressure coming from the
diffusive species, respectively. To use this command a reaction must be
set up.

.. [1]
   http://www.paraview.org/

.. [2]
   http://code.enthought.com/projects/mayavi/
