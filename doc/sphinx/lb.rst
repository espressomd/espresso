Lattice-Boltzmann
=================

For an implicit treatment of a solvent, |es| allows to couple the molecular
dynamics simulation to a Lattice-Boltzmann fluid. The Lattice-Boltzmann-Method (LBM) is a fast, lattice based method that, in its
"pure" form, allows to calculate fluid flow in different boundary
conditions of arbitrarily complex geometries. Coupled to molecular
dynamics, it allows for the computationally efficient inclusion of
hydrodynamic interactions into the simulation. The focus of the |es| implementation
of the LBM is, of course, the coupling to MD and therefore available
geometries and boundary conditions are somewhat limited in comparison to
"pure" LB codes.

Here we restrict the documentation to the interface. For a more detailed
description of the method, please refer to the literature.

.. note:: Please cite :cite:`espresso2` (Bibtex key espresso2 in doc/sphinx/zref.bib) if you use the LB fluid and :cite:`lbgpu` (Bibtex key lbgpu in doc/sphinx/zref.bib) if you use the GPU implementation.

Setting up a LB fluid
---------------------

The following minimal example illustrates how to use the LBM in |es|::

    import espressomd
    sys = espressomd.System()
    sys.box_l = [10, 20, 30]
    sys.time_step = 0.01
    lb = espressomd.lb.LBFluid(agrid=1.0, dens=1.0, visc=1.0, tau=0.01)
    sys.actors.add(lb)
    sys.integrator.run(100)

.. note:: `Feature LB or LB_GPU required`

To use the (much faster) GPU implementation of the LBM, use
:class:`espressomd.lb.LBFluid_GPU` in place of :class:`espressomd.lb.LBFluid`.

The command initializes the fluid with a given set of parameters. It is
also possible to change parameters on the fly, but this will only rarely
be done in practice. Before being able to use the LBM, it is necessary
to set up a box of a desired size. The parameter is used to set the
lattice constant of the fluid, so the size of the box in every direction
must be a multiple of ``agrid``.

In the following, we discuss the parameters that can be supplied to the LBM in |es|. The detailed interface definition is available at :class:`espressomd.lb.LBFluid`.

In the LB scheme and the MD scheme are not synchronized: In one LB time
step typically several MD steps are performed. This allows to speed up
the simulations and is adjusted with the parameter ``tau``, the LB timestep.
The parameters ``dens`` and ``visc`` set up the density and (kinematic) viscosity of the
LB fluid in (usual) MD units. Internally the LB implementation works
with a different set of units: all lengths are expressed in ``agrid``, all times
in ``tau`` and so on.
LB nodes are located at 0.5, 1.5, 2.5, etc.
(in terms of ``agrid``). This has important implications for the location of
hydrodynamic boundaries which are generally considered to be halfway
between two nodes for flat, axis-aligned walls. For more complex boundary geometries, the hydrodynamic boundary location deviates from this midpoint and the deviation decays like to first order in ``agrid``. 
The LBM should
*not be used as a black box*, but only after a careful check of all
parameters that were applied.

In the following, we describe a number of optional parameters.
The parameter ``ext_force`` takes a three dimensional vector as an `array_like`, representing a homogeneous external body force density in MD units to be applied to the fluid. The
parameter allows to tune the bulk viscosity of the fluid and is given in
MD units. In the limit of low Mach number, the flow does not compress the fluid and the resulting flow field is therefore independent of the bulk viscosity. It is however known that the values of the viscosity does affect
the quality of the implemented link-bounce-back method. and are the
relaxation parameters for the kinetic modes. Due to their somewhat
obscure nature they are to be given directly in LB units.

Before running a simulation at least the following parameters must be
set up: , , , , . For the other parameters, the following are taken: =0,
=0, =0, = = = 0.

If the feature is activated, the Lattice Boltzmann code (so far GPU
version only) is extended to a two-component Shan-Chen (SC) method. The
command requires in this case to supply two values, for the respective
fluid components, to each of the options , , , , and , when they are
used, otherwise they are set to the default values. The three elements
of the coupling matrix can be supplied with the option , and the
mobility coefficient can be specified with the option . By default no
copuling is activated, and the relaxation parameter associated to the
mobility is zero, corresponding to an infinite value for . Additional
details are given in [sec:shanchen] and [sec:scmd-coupling].

lbfluid print\_interpolated\_velocity

This variant returns the velocity at point in countinous space. This can
make it easier to calculate flow profiles independent of the lattice
constant.

| lbfluid save\_ascii\_checkpoint
| lbfluid save\_binary\_checkpoint
| lbfluid load\_ascii\_checkpoint
| lbfluid load\_binary\_checkpoint

The first two save commands save all of the LB fluid nodes’ populations
to in ascii or binary format respectively. The two load commands load
the populations from . This is useful for restarting a simulation either
on the same machine or a different machine. Some care should be taken
when using the binary format as the format of doubles can depend on both
the computer being used as well as the compiler. One thing that one
needs to be aware of is that loading the checkpoint also requires the
used to reuse the old forces. This is necessary since the coupling force
between the paricles and the fluid has already been applied to the
fluid. Failing to reuse the old forces breaks momentum conservation,
which is in general a problem. It is particularly problematic for bulk
simulations as the system as a whole acquires a drift of the center of
mass, causing errors in the calculation of velocities and diffusion
coefficients. The correct way to restart an LB simulation is to first
load in the particles with the correct forces, and use “integrate
*steps* reuse\_forces” upon the first call to integrate. This causes the
old forces to be reused and thus conserves momentum.

LB as a thermostat
------------------

thermostat

The LBM implementation in uses Ahlrichs and Dünweg’s point coupling
method to couple MD particles the LB fluid. This coupling consists in a
frictional force and a random force:

.. math:: \vec{F} = -\gamma \left(\vec{v}-\vec{u}\right) + \vec{F}_R.

The momentum acquired by the particles is then transferred back to the
fluid using a linear interpolation scheme, to preserve total momentum.
In the GPU implementation the force can alternatively be interpolated
using a three point scheme which couples the particles to the nearest 27
LB nodes. This can be called using “lbfluid 3pt” and is described in
Dünweg and Ladd by equation 301 :cite:`duenweg08a`. Note that
the three point coupling scheme is incompatible with the Shan Chen
Lattice Boltmann. The frictional force tends to decrease the relative
velocity between the fluid and the particle whereas the random forces
are chosen so large that the average kinetic energy per particle
corresponds to the given temperature, according to a fluctuation
dissipation theorem. No other thermostatting mechanism is necessary
then. Please switch off any other thermostat before starting the LB
thermostatting mechanism.

The LBM implementation provides a fully thermalized LB fluid, all
nonconserved modes, including the pressure tensor, fluctuate correctly
according to the given temperature and the relaxation parameters. All
fluctuations can be switched off by setting the temperature to 0.

Regarind the unit of the temperature, please refer to
Section [sec:units].

The Shan Chen bicomponent fluid
-------------------------------

Please cite  if you use the Shan Chen implementation described below.

The Lattice Boltzmann variant of Shan and
Chan :cite:`shan93a` is widely used as it is simple and yet
very effective in reproducing the most important traits of
multicomponent or multiphase fluids. The version of the Shan-Chen method
implemented in is an extension to bi-component fluids of the
multi-relaxation-times Lattice Boltzmann with fluctuations applied to
all modes, that is already present in . It features, in addition,
coupling with particles :cite:`sega13c` and
component-dependent particle interactions (see sections
[sec:scmd-coupling] and[sec:scmd-affinity]).

The Shan-Chen fluid is set up using the command, supplying two values
(one per component) to the option. Optionally, two values can be set for
each of the usual transport coefficients (shear and bulk viscosity), and
for the ghost modes. It is possible to set a relaxation time also for
the momentum modes, since they are not conserved quantities in the
Shan-Chen method, by using the option . The mobility transport
coefficient expresses the propensity of the two components to mutually
diffuse, and, differently from other transport coefficients, only one
value is needed, as it carachterizes the mixture as a whole. When
thermal fluctuations are switched on, a random noise is added, in
addition, also to the momentum modes. Differently from the other modes,
a correlated noise is added to the momentum ones, in order to preserve
the *total* momentum.

The fluctuating hydrodynamic equations that are simulated using the
Shan-Chen approach are

.. math::

   \label{eq:shanchen-NS}
   \rho \left(\frac{\partial }{\partial  t} {\vec {u}} + ({\vec {u}}\cdot {\vec {\nabla}})  {\vec {u}} \right)=-{\vec {\nabla}} p+{\vec {\nabla}} \cdot ({\vec {\Pi}}+\hat{{\vec {\sigma}}})+\sum_{\zeta} {\vec {g}}_{\zeta},

.. math::

   \label{eq:shanchen-cont}
   \frac{\partial }{\partial  t} \rho_{\zeta}+{\vec {\nabla}} \cdot (\rho_{\zeta} {\vec {u}}) = {\vec {\nabla}} \cdot  ({\vec {D}}_{\zeta}+\hat{{\vec {\xi}}}_{\zeta}),

.. math::

   \label{eq:shanchen-globalcont}
   \partial_t \rho+{\vec {\nabla}} \cdot (\rho {\vec {u}}) = 0,

where the index :math:`\zeta=1,2` specifies the component,
:math:`\vec{u}` is the fluid (baricentric) velocity,
:math:`\rho=\sum_\zeta\rho_\zeta` is the total density, and
:math:`p=\sum_{\zeta} p_{\zeta}=\sum_{\zeta} c_s^2
\rho_{\zeta}` is the internal pressure of the mixture (:math:`c_s` being
the sound speed). Two fluctuating terms :math:`\hat{{\vec{\sigma}}}` and
:math:`\hat{{\vec{\xi}}}_{\zeta}` are associated, respectivelu, to the
diffusive current :math:`{\vec{D}}_{\zeta}` and to the viscous stress
tensor :math:`{\vec{\Pi}}`.

The coupling between the fluid components is realized by the force

.. math::

   \vec{g}_{\zeta}(\vec{r}) =  - \rho_{\zeta}(\vec{r})
    \sum_{\vec{r}'}\sum_{\zeta'}  g_{\zeta \zeta'} \rho_{\zeta'}
    (\vec{r}') (\vec{r}'-\vec{r}),

that acts on the component :math:`\zeta` at node position
:math:`\vec{r}`, and depends on the densities on the neighboring nodes
located at :math:`\vec{r}'`. The width of the interfacial regions
between two components, that can be obtained with the Shan-Chen method
is usually 5-10 lattice units. The coupling matrix
:math:`g_{\zeta \zeta'}` is in general symmetric, so in the present
implementation only three real values need to be specified with the
option . The command sets the density of the two components to the
values specified by the option , and these can be modified with the
command. Note that the number of active fluid components can be accessed
through the global variable .

SC as a thermostat
------------------

The coupling of particle dynamics to the Shan-Chen fluid has been
conceived as an extension of the Ahlrichs and Dünweg’s point coupling,
with the force acting on a particle given by

.. math:: \vec{F} = -\frac{\sum_\zeta \gamma_\zeta \rho_\zeta(\vec{r})}{\sum_\zeta \rho_\zeta(\vec{r}_\zeta)} \left(\vec{v}-\vec{u}\right) + \vec{F}_R + \vec{F}^{ps},

where :math:`\zeta` identifies the component,
:math:`\rho_\zeta(\vec{r})` is a linear interpolation of the component
density on the nodes surrounding the particle, :math:`\gamma_\zeta` is
the component-dependent friction coefficient, :math:`\vec{F}_R` is the
usual random force, and

.. math:: \vec{F}^{\mathrm{ps}}= -  \sum_{\zeta} \kappa_{\zeta} \nabla \rho_{\zeta}(\vec{r}).

This is an effective solvation force, that can drive the particle
towards density maxima or minima of each component, depending on the
sign of the constant :math:`\kappa_\zeta`. Note that by setting the
coupling constant to the same negative value for both components will,
in absence of other forces, push the particle to the interfacial region.

In addition to the solvation force acting on particles, another one that
acts on the fluid components is present, representing the solvation
force of particles on the fluid.

.. math:: \vec{F}_{\zeta}^{\mathrm{fs}}(\vec{r}) = -\lambda_{\zeta} \rho_{\zeta}(\vec{r}) \sum_i \sum_{\vec{r}'} \Theta \left[\frac{(\vec{r}_i-\vec{r})}{\|\vec{r}_i-\vec{r}\|} \cdot \frac{(\vec{r}'-\vec{r})}{\|\vec{r}'-\vec{r}\|} \right] \frac{\vec{r}'-\vec{r}}{\|\vec{r}'-\vec{r}\|^2},

where :math:`\Theta(x)=1` if :math:`0<x<1`, and 0 otherwise, the sum
over lattice nodes is performed on the neighboring sites of
:math:`\vec{r}` and the index :math:`i` runs over all particles. Note
that a dependence on the particle index :math:`i` is assumed for
:math:`\kappa_\zeta` and :math:`\lambda_\zeta`. This force has the
effect of raising or lowering (depending on the sign of the coupling
constant :math:`\lambda_\zeta`) the density in the eight nodes around a
particle. The particle property (Chap. [chap:part]) sets the coupling
constants :math:`\lambda_A`,\ :math:`\kappa_A`,\ :math:`\lambda_B` and
:math:`\kappa_B`, where :math:`A` and :math:`B` denote the first and
second fluid component, respectively. A complete description of the
copuling scheme can be found in :cite:`sega13c`.

SC component-dependent interactions between particles
-----------------------------------------------------

Often particle properties depend on the type of solvent in which they
are. For example, a polymer chain swells in a good solvent, and
collapses in a bad one. One of the possible ways to model the good or
bad solvent condition in coarse-grained models is to employ a WCA or a
LJ (attractive) potential, respectively. If one wants to model the two
components of the SC fluid as good/bad solvent, it is possible to do it
using the argument of the command. This non-bonded interaction type acts
as a modifier to other interactions. So far only the Lennard-Jones
interaction is changed by the , so that it switches in a continuous way
(after the potential minimum) from the full interaction to the WCA one.
For more information see [sec:LennardJones] and [sec:affinity].

Reading and setting single lattice nodes
----------------------------------------

lbnode

| The command allows to inspect () and modify () single LB nodes. Note
  that the indexing in every direction starts with 0. For both commands
  you have to specify what quantity should be printed or modified. Print
  allows the following arguments:

+------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|                              | the density (one scalar\ :math:`^{1,2}` or two scalars\ :math:`^3`).                                                                                                                                |
+------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|                              | the fluid velocity (three floats: :math:`u_x`, :math:`u_y`, :math:`u_z`)                                                                                                                            |
+------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|                              | the fluid velocity (six floats: :math:`\Pi_{xx}`, :math:`\Pi_{xy}`, :math:`\Pi_{yy}`, :math:`\Pi_{xz}`, :math:`\Pi_{yz}`, :math:`\Pi_{zz}`)                                                         |
+------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|                              | the nonequilbrium part of the pressure tensor, components as above.                                                                                                                                 |
+------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|                              | the 19 population (check the order from the source code please).                                                                                                                                    |
+------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
|                              | the flag indicating whether the node is a fluid node (:math:`\lit{boundary}=0`) or a boundary node (:math:`\lit{boundary}\ne 0`). Does not support . Refer to the command for this functionality.   |
+------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
| :math:`^1` or ; :math:`^2`   |                                                                                                                                                                                                     |
+------------------------------+-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+

Example: The line

puts [ lbnode 0 0 0 print u ]

prints the fluid velocity in node 0 0 0 to the screen. The command
allows to change the density or fluid velocity in a single node. Setting
the other quantities can easily be implemented. Example:

puts [ lbnode 0 0 0 set u 0.01 0. 0.]

Removing total fluid momentum
-----------------------------

lbfluid remove\_momentum

In some cases, such as free energy profile calculations, it might be
useful to prevent interface motion. This can be achieved using the
command , that removes the total momentum of the fluid.

Visualization
-------------

lbfluid print lbfluid print vtk velocity

The print parameter of the command is a feature to simplify
visualization. It allows for the export of the whole fluid field data
into a file with name at once. Currently supported values for the
parameter are boundary and velocity when using or and density and
velocity when using . The additional option enables export in the vtk
format which is readable by visualization software such as paraview [1]_
or mayavi2 [2]_. Otherwise gnuplot readable data will be exported. If
you plan to use paraview for visualization, note that also the particle
positions can be exported in the VTK format [sec:writevtk]. allows you
to only output part of the flow field by specifiying an axis aligned
bounding box through the coordinates of two of its corners. This
bounding box can be used to output a slice of the flow field. As an
example, executing ``lbfluid print vtk velocity 0 0 5 10 10 5 filename``
will output the cross-section of the velocity field in a plane
perpendicular to the :math:`z`-axis at :math:`z = 5` (assuming the box
size is 10 in the :math:`x`- and :math:`y`-direction). If the
bicomponent fluid is used, two filenames have to be supplied when
exporting the density field, to save both components.

Setting up boundary conditions
------------------------------

lbboundary lbboundary force

If nothing else is specified, periodic boundary conditions are assumed
for the LB fluid. Variant allows to set up other (internal or external)
boundaries.

The command syntax is very close to the syntax, as usually one wants the
hydrodynamic boundary conditions to be shaped similarily to the MD
boundaries. Currently the shapes mentioned above are available and their
syntax exactly follows the syntax of the constraint command. For example

lbboundary wall dist 1.5 normal 1. 0. 0.

creates a planar boundary condition at distance 1.5 from the origin of
the coordinate system where the half space :math:`x>1.5` is treated as
normal LB fluid, and the other half space is filled with boundary nodes.

Intersecting boundaries are in principle possible but must be treated
with care. In the current, only partly satisfactory, all nodes that are
within at least one boundary are treated as boundary nodes. Improving
this is nontrivial, and suggestions are very welcome.

Currently, only the so called “link-bounce-back” algorithm for wall
nodes is available. This creates a boundary that is located
approximately midway between the lattice nodes, so in the above example
this corresponds indeed to a boundary at :math:`x=1.5`. Note that the
location of the boundary is unfortunately not entirely independent of
the viscosity. This can be seen when using the sample script with a high
viscosity.

The bounce back boundary conditions allow to set velocity at a boundary
to a nonzero value. This allows to create shear flow and boundaries
moving relative to each other. This could be a fixed sphere in a channel
moving at a finite speed – corresponding to the galilei-transform of a
moving sphere in a fixed channel. The velocity boundary conditions are
implemented according to :cite:`succi01a` eq. 12.58. Using
this implementation as a blueprint for the boundary treatment an
implementation of the Ladd-Coupling should be relatively
straightforward.

Variant prints out the force on boundary number .

Choosing between the GPU and CPU implementations
------------------------------------------------

lbfluid cpu lbfluid gpu

A very recent development is an implementation of the LBM for NVIDIA
GPUs using the CUDA framework. On CUDA-supporting machines this can be
activated by configuring with and activating the feature . Within the
-Tcl-script, the command can be used to choose between the CPU and GPU
implementations of the Lattice-Boltzmann algorithm, for further
information on CUDA support see section [sec:cuda].

Variant is the default and turns on the standard CPU implementation of
the Lattice-Boltzmann fluid, while variant turns on the GPU
implementation, implying that all following LB-related commands are
executed on the GPU.

Currently only a subset of the CPU commands are available for the GPU
implementation. For boundary conditions analogous to the CPU
implementation, the feature has to be activated.

Electrohydrodynamics
--------------------

setmd mu\_E

If the feature is activated, the (non-GPU) Lattice Boltzmann Code can be
used to implicitely model surrounding salt ions in an external electric
field by having the charged particles create flow.

For that to work, you need to set the electrophoretic mobility
(multiplied by the external :math:`E`-field) :math:`\mu E` in all 3
dimensions for your system. The three given parameters are float values
and should, for a meaningful system, be less than :math:`1.0`.

For more information on this method and how it works, read the
publication :cite:`hickey10a`.

.. [1]
   http://www.paraview.org/

.. [2]
   http://code.enthought.com/projects/mayavi/
