====================
Setting up particles
====================

Creating single particles
=========================


Defining particle properties
----------------------------

The first step when writing a simulation script is to import the :mod:`espressomd`
module and to create a :class:`espressomd.system.System` instance::

    >>> import espressomd

    >>> system = espressomd.System()

In order to add particles to the system, call
:meth:`espressomd.particle_data.ParticleList.add`::

    >>> system.part.add(pos=[1.0, 1.0, 1.0], id=0, type=0)

This command adds a single particle to the system with properties given
as arguments. All available particle properties are members of
:class:`espressomd.particle_data.ParticleHandle` and are listed below.

    - :attr:`espressomd.particle_data.ParticleHandle.bonds`
    - :attr:`espressomd.particle_data.ParticleHandle.dip`

        ..  note::

            `Feature DIPOLES required.`

    - :attr:`espressomd.particle_data.ParticleHandle.dipm`

        ..  note::

            `Feature DIPOLES required.`

    - :attr:`espressomd.particle_data.ParticleHandle.director`

        ..  note::

            `Feature ROTATION required.`

    - :attr:`espressomd.particle_data.ParticleHandle.exclude`

        ..  note::

            `Feature EXCLUSIONS required.`

    - :attr:`espressomd.particle_data.ParticleHandle.ext_force`

        ..  note::
            
            `Feature EXTERNAL_FORCE required.`

    - :attr:`espressomd.particle_data.ParticleHandle.ext_torque`

        ..  note::
            
            `Feature ROTATION and EXTERNAL_FORCE required.`

    - :attr:`espressomd.particle_data.ParticleHandle.f`
    - :attr:`espressomd.particle_data.ParticleHandle.fix`

        ..  note::
            
            `Feature EXTERNAL_FORCE required.`

    - :attr:`espressomd.particle_data.ParticleHandle.gamma`

        ..  note::
            
            `Feature LANGEVIN_PER_PARTICLE required.`

    - :attr:`espressomd.particle_data.ParticleHandle.gamma_rot`

        ..  note::
            
            `Feature LANGEVIN_PER_PARTICLE, ROTATION and ROTATIONAL_INERTIA required.`

    - :attr:`espressomd.particle_data.ParticleHandle.type`
    - :attr:`espressomd.particle_data.ParticleHandle.pos`
    - :attr:`espressomd.particle_data.ParticleHandle.pos_folded`
    - :attr:`espressomd.particle_data.ParticleHandle.mass`

        ..  note::

            `Feature MASS required.`

    - :attr:`espressomd.particle_data.ParticleHandle.omega_body`

        ..  note::

            `Feature ROTATION required.`

    - :attr:`espressomd.particle_data.ParticleHandle.omega_lab`

        ..  note::

            `Feature ROTATION required.`

    - :attr:`espressomd.particle_data.ParticleHandle.q`

        ..  note::

            `Feature ELECTROSTATICS required.`

    - :attr:`espressomd.particle_data.ParticleHandle.quat`

        ..  note::

            `Feature ROTATION required.`

    - :attr:`espressomd.particle_data.ParticleHandle.rotation`

        ..  note::

            `Feature ROTATION required.`

    - :attr:`espressomd.particle_data.ParticleHandle.rinertia`

        ..  note::

            `Feature ROTATIONAL_INERTIA required.`

    - :attr:`espressomd.particle_data.ParticleHandle.smaller_timestep`

        ..  note::

            `Feature MULTI_TIMESTEP required.`

    - :attr:`espressomd.particle_data.ParticleHandle.swimming`

        ..  note::

            `Feature ENGINE required.`

    - :attr:`espressomd.particle_data.ParticleHandle.temp`

        ..  note::
            
            `Feature LANGEVIN_PER_PARTICLE required.`

    - :attr:`espressomd.particle_data.ParticleHandle.torque_lab`

        ..  note::

            `Feature ROTATION required.`

    - :attr:`espressomd.particle_data.ParticleHandle.v`
    - :attr:`espressomd.particle_data.ParticleHandle.virtual`

        ..  note::

            `Feature VIRTUAL_SITES required.`

    - :attr:`espressomd.particle_data.ParticleHandle.vs_relative`

        ..  note::

            `Feature VIRTUAL_SITES required.`

Properties of already existing particles can be set using::

    >>> system.part[<ID>].<PROPERTY> = <SOME_VALUE>

This sets the property ``PROPERTY`` for the particle with id ``ID`` to 
``SOME_VALUE``.

Getting particle properties
~~~~~~~~~~~~~~~~~~~~~~~~~~~

If a certain particle property is set, it can be accessed like a class
member. To access property ``PROPERTY`` of the particle with id ``ID``::

    system.part[<ID>].<PROPERTY>

For example, to print the current position of all particles in the system, call::

    >>> print(system.part[:].pos)
    [[ 3.77651228  2.74802277  7.8614655 ]
     [ 3.16587857  2.88714253  3.0495119 ]
     [ 3.13657179  6.26879849  1.1182947 ]
     ..., 
     [ 1.42691672  8.39444662  7.61834009]
     [ 3.63801655  3.05804937  5.30344039]
     [ 8.13603676  3.91630721  2.70701524]]

Deleting particles
~~~~~~~~~~~~~~~~~~

Particles can be easily deleted in Python using particle ids or ranges of particle ids. For example, to delete all particles with particle id greater than 10, run::

    >>> system.part[10:].remove()

TODO: Check following text
--------------------------

Enables the particle to be self-propelled in the direction determined by
its quaternion. For setting the quaternion of the particle see . The
self-propulsion speed will relax to a constant velocity, that is
specified by . Alternatively it is possible to achieve a constant
velocity by imposing a constant force term that is balanced by friction
of a (Langevin) thermostat. The way the velocity of the particle decays
to the constant terminal velocity in either of these methods is
completely determined by the friction coefficient. You may only set one
of the possibilities *or* as you cannot relax to constant force *and*
constant velocity at the same time. The option (re)sets and both to
:math:`0.0` and thus disables swimming. This option applies to all
non-lattice-Boltzmann thermostats. Note that there is no real difference
between and , since the latter may aways be chosen such that the same
terminal velocity is achieved for a given friction coefficient.

For an explanation of the parameters , and see the previous item. In
lattice-Boltzmann self-propulsion is less trivial than for normal MD,
because the self-propulsion is achieved by a force-free mechanism, which
has strong implications for the far-field hydrodynamic flow field
induced by the self-propelled particle. In only the dipolar component of
the flow field of an active particle is taken into account. This flow
field can be generated by a *pushing* or a *pulling* mechanism, leading
to change in the sign of the dipolar flow field with respect to the
direction of motion. You can specify the nature of the particle’s flow
field by using the keywords or . You will also need to specify a which
determines the distance of the source of propulsion from the particle’s
center. Note that you should not put this distance to zero; (currently)
does not support mathematical dipole flow fields. The key can be used to
set the friction that causes the orientation of the particle to change
in shear flow. The torque on the particle is determined by taking the
cross product of the difference between the fluid velocity at the center
of the particle and at the source point and the vector connecting the
center and source.

You may ask: “Why are there two methods and for the self-propulsion
using the lattice-Bolzmann algorithm?” The answer is straightforward.
When a particle is accelerating, it has a monopolar flow-field
contribution which vanishes when it reaches its terminal velocity (for
which there will only be a dipolar flow field). The major difference
between the above two methods is that with the flow field *only* has a
monopolar moment and *only* while the particle is accelerating. As soon
as the particle reaches a constant speed (given by ) this monopolar
moment is gone and the flow field is zero! In contrast, always, i.e.,
while accelerating *and* while swimming at constant force possesses a
dipolar flow field.

Variant will return a list of the specified properties of particle , or
all properties, if no keyword is specified. Variant will return a list
of all properties of all particles.

Note that there is a difference between the and . The first prints the
variable in the co-rotating frame, whereas the second gives the variable
in the stationary frame, the body and laboratory frames, respectively.
One would typically want to output the variable in the laboratory frame,
since it is the frame of interest. However for some tests involving
reading and writing the variable it may be desireable to know it in the
body frame as well. Be careful with reading and writing, if you write in
the lab frame, then read in the lab frame. If you are setting the
variable in the lab frame, the orientation of the particle’s must be set
before, otherwise the conversion from lab to body frame will not be
handled properly. Also be careful about the order in which you write and
read in data from a blockfile, for instance if you output the variable
in both frames!

The command is a print-only command that gives the velocity in the body
frame, which can be useful for determining the translational diffusion
tensor of an anisotropic particle via the velocity auto-correlation
(Green-Kubo) method.

part 40 print id pos q bonds

will return a list like

40 8.849 1.8172 1.4677 1.0

This routine is primarily intended for effective use in Tcl scripts.

When the keyword is specified, it returns the connectivity of the
particle up to (defaults to 1). For particle 5 in a linear chain the
result up to = 3 would look like:

 4 6 4 3 6 7 4 3 2 6 7 8

The function is useful when you want to create bonded interactions to
all other particles a certain particle is connected to. Note that this
output can not be used as input to the part command. Check results if
you use them in ring structures.

If none of the options is specified, it returns all properties of the
particle, if it exists, in the form

0 pos 2.1 6.4 3.1 type 0 q -1.0 v 0.0 0.0 0.0 f 0.0 0.0 0.0 bonds 0 480
0 368 ...

which may be used as an input to this function later on. The first
integer is the particle number.

Variant returns the properties of all stored particles in a tcl-list
with the same format as specified above:

0 pos 2.1 6.4 3.1 type 0 q -1.0 v 0.0 0.0 0.0 f 0.0 0.0 0.0 bonds0 4800
368... 1 pos 1.0 2.0 3.0 type 0 q 1.0 v 0.0 0.0 0.0 f 0.0 0.0 0.0 bonds0
3400 83... 2......... 3......... ...

When using ``pos``, the particle position returned is **unfolded**, for
convenience in diffusion calculations etc. Note that therefore
blockfiles will contain imaged positions, but un-imaged velocities,
which should not be interpreted together. However, that is fine for
restoring the simulation, since the particled data is loaded the same
way.

Exclusions
~~~~~~~~~~

part auto\_exclusions part delete\_exclusions

Variant will create exclusions for all particles pairs connected by not
more than bonds ( defaults to 2). This is typically used in atomistic
simulations, where nearest and next nearest neighbour interactions along
the chain have to be omitted since they are included in the bonding
potentials. For example, if the system contains particles :math:`0`
…\ :math:`100`, where particle :math:`n` is bonded to particle
:math:`n-1` for :math:`1 \leq n \leq 100`, then it will result in the
exclusions:

-  particle 1 does not interact with particles 2 and 3

-  particle 2 does not interact with particles 1, 3 and 4

-  particle 3 does not interact with particles 1, 2, 4 and 5

-  ...

Variant deletes all exclusions currently present in the system.

Creating groups of particle
---------------------------

``polymer``: Setting up polymer chains
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

:: 

    from espressomd.polymer import create_polymer

A function that allow to create a number of polymers and polyelectrolytes. 
See :attr:`espressomd.polymer.create_polymer()` for a detailed list of 
arguments. 

The distance between adjacent monomers
during the course of the simulation depends on the applied potentials.
For fixed bond length please refer to the Rattle Shake
algorithm:raw-latex:`\cite{andersen83a}`. The algorithm is based on
Verlet algorithm and satisfy internal constraints for molecular models
with internal constrains, using Lagrange multipliers.

The polymer can be created using several different random walk modes:

 (Random walk)
    mode = 1 The monomers are randomly placed by a random walk with a
    steps size of ``bond_length``.

 (Pruned self-avoiding walk)
    mode = 2 The position of a monomer is randomly chosen in a distance
    of to the previous monomer. If the position is closer to another
    particle than ``shield``, the attempt is repeated up to ``max_tries`` times. Note, that this
    is not a real self-avoiding random walk, as the particle
    distribution is not the same. If you want a real self-avoiding walk, use
    the mode 0. However, this mode is several orders of magnitude faster than a
    true self-avoiding random walk, especially for long chains.

 (Self-avoiding random walk)
    mode = 0 The positions of the monomers are chosen as in the plain
    random walk. However, if this results in a chain that has a monomer
    that is closer to another particle than ``shield``, a new attempt of setting
    up the whole chain is done, up to ``max_tries`` times.


``counterions``: Setting up counterions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

counterions

This command will create counterions in the simulation box.

Sets the particle id of the first counterion. It defaults to the current
number of particles, counterions are placed after all previously defined
particles.

Specifies the setup method to place the counterions. It defaults to .
See the command for a detailed description.

Specifies the charge of the counterions. If not set, it defaults to
:math:`-1.0`.

Specifies the particle type of the counterions. It defaults to
:math:`2`.

``salt``: Setting up salt ions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

salt

Create positively and negatively charged salt ions of charge and within
the simulation box.

Sets the particle id of the first (positively charged) salt ion. It
defaults to the current number of particles.

Specifies the setup method to place the counterions. It defaults to .
See the command for a detailed description.

Sets the charge of the positive salt ions to and the one of the
negatively charged salt ions to . If not set, the values default to
:math:`1.0` and :math:`-1.0`, respectively.

Specifies the particle type of the salt ions. It defaults to :math:`3`
respectively :math:`4`.

The salt ions are only placed in a sphere with radius around the origin.

``diamond``: Setting up diamond polymer networks
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
::

    from espressomd import Diamond

Creates a diamond-shaped polymer network with 8 tetra-functional nodes
connected by :math:`2*8` polymer chains of length (MPC) in a unit cell
of length :math:`a`. Chain monomers are placed at a mutual distance along the
vector connecting network nodes. The polymer is created starting from
particle ID 0. Nodes are assigned type 0, monomers (both charged and
uncharged) are type 1 and counterions type 2. For inter-particle bonds
interaction :math:`0` is taken which must be a two-particle bond.


.. figure:: figures/diamond.png
   :alt: Diamond-like polymer network with MPC=15.
   :align: center
   :height: 6.00000cm

   Diamond-like polymer network with MPC=15.

See :meth:`espressomd.diamond.Diamond` for more details. 

``icosaeder``: Setting up an icosaeder
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

icosaeder

Creates a modified icosaeder to model a fullerene (or soccer ball). The
edges are modeled by polymer chains connected at the corners of the
icosaeder. For inter-particle bonds interaction :math:`0` is taken which
must be a two-particle bond. Two particle types are used for the
pentagons and the interconnecting links. For an example, see figure
[fig:fullerene].

.. figure:: figures/fullerene.png
   :alt: Icosaeder with =15.
   :align: center
   :height: 6.00000cm

   Icosaeder with =15.

Length of the links. Defines the size of the icosaeder.

Specifies the number of chain monomers along one edge.

Specifies the number of counterions to be placed into the system.

Set the charges of the monomers to and the charges of the counterions to
.

Specifies the distance between two charged monomer along the edge. If
:math:`\var{d_\mathrm{charged}} > 1` the remaining monomers are
uncharged.

``crosslink``: Cross-linking polymers
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

crosslink

Attempts to end-crosslink the current configuration of equally long
polymers with monomers each, returning how many ends are successfully
connected.

specifies the first monomer of the chains to be linked. It has to be
specified if the polymers do not start at id 0.

Set the radius around each monomer which is searched for possible new
monomers to connect to. defaults to :math:`1.9`.

The minimal distance of two interconnecting links. It defaults to
:math:`2`.

The minimal distance for an interconnection along the same chain. It
defaults to :math:`0`. If set to , no interchain connections are
created.

Sets the bond type for the connections to .

If not specified, defaults to :math:`30000`.

``constraint``: Setting up constraints
--------------------------------------

constraint wall normal dist type

constraint sphere center radius direction type

constraint cylinder center axis radius length direction type

| constraint rhomboid corner a b
| c direction type

constraint maze nsphere dim sphrad cylrad type

constraint pore center axis radius length type

constraint stomatocyte center orientation outer\_radius inner\_radius
layer\_width direction type

constraint slitpore pore\_mouth channel\_width pore\_width pore\_length
upper\_smoothing\_radius lower\_smoothing\_radius

constraint mindist\_position

constraint hollow\_cone center orientation outer\_radius inner\_radius
width opening\_angle direction type

constraint spherocylinder center axis radius length direction type

The command offers a variety of surfaces that can be defined to interact
with desired particles. Variants to create interactions via a non-bonded
interaction potential, where the distance between the two particles is
replaced by the distance of the center of the particle to the surface.

.. warning:: When using shapes with concave edges and corners, the fact that a particle only interacts with the closest point on the constraint surface leads to discontinuous force fields acting on the particles. This breaks energy conservation in otherwise symplectic integrators. Often, the total energy of the system increases exponentially.

The constraints are identified like a particle via its type for the
non-bonded interaction. After a type is defined for each constraint one
has to define the interaction of all different particle types with the
constraint using the command. In variants to , constraints are able to
be penetrated if is set to 1. Otherwise, when the penetrable option is
ignored or is set to 0, the constraint cannot be violated, i.e. no
particle can go through the constraint surface. In variants to and it is
also possible to specify a flag indicating if the constraints should be
reflecting. The flags can equal 1 or 2. The flag 1 corresponds to a
reflection process where the normal component of the velocity is
reflected and the tangential component remains unchanged. If the flag is
2, also the tangential component is turned around, so that a bounce back
motion is performed. The second variant is useful for boundaries of DPD.
The reflection property is only activated if an interaction is defined
between a particular particle and the constraint! This will usually be a
lennard-jones interaction with :math:`\epsilon=0`, but finite
interaction range.

In variant if the flag is set to 1, interactions are only calculated if
the particle is on the side of the wall in which the normal vector is
pointing. This has only an effect for penetrable walls. If the flag is
set to 1, then slip boundary interactions apply that are essential for
microchannel flows like the Plane Poiseuille or Plane Couette Flow. You
also need to use the tunable\_slip interaction (see [sec:tunableSlip])
for this too work.

Variants and create interactions based on electrostatic interactions.
The corresponding force acts in direction of the normal vector of the
surface and applies to all charged particles. For the normal vector
which is used in the implementation lies in z-direction.

Variant does not define a surface but is based on magnetic dipolar
interaction with an external magnetic field. It applies to all particles
with a dipole moment.

The resulting surface in variant is a plane defined by the normal vector
and the distance from the origin (in the direction of the normal
vector). The force acts in direction of the normal. Note that the
describes the distance from the origin in units of the normal vector so
that the product of :math:`d` and :math:`n` is a point on the surface.
Therefore negative distances are quite common!

The resulting surface in variant is a sphere with center and radius .
The determines the force direction, -1 or for inward and +1 or for
outward.

The resulting surface in variant is a cylinder with center and radius .
The parameter is **half** of the cylinder length. The is a vector along
the cylinder axis, which is normalized in the program. The is defined
the same way as for the spherical constraint.

The resulting surface in variant is a rhomboid, defined by one corner
located at and three adjacent edges, defined by the three vectors
connecting the corner p with it’s three neighboring corners, a ( ), b (
) and c ( ).

The resulting surface in variant is spheres of radius along each
dimension, connected by cylinders of radius . The spheres have simple
cubic symmetry. The spheres are distributed evenly by dividing the by .
Dimension of the maze can be controlled by : 0 for one dimensional, 1
for two dimensional and 2 for three dimensional maze.

Variant sets up a cylindrical pore similar to variant with a center and
radius . The parameter is **half** of the cylinder length. The is a
vector along the cylinder axis, which is normalized in the program.
Optionally the outer radius of the pore can be specified. By default
this is (numerical) infinity and thus results in an infinite wall with
one pore. The argument can be replaced by the argument to obtain a pore
with a conical shape and corresponding opening radii. The first radius
is in the direction opposite to the axis vector. The same applies for
which can be replaced with . Per default sharp edges are replaced by
circles of unit radius. The radius of this smoothing can be set with the
optional keyword .

Variant creates a stomatocyte shaped boundary. This command should be
used with care. The position can be any point in the simulation box, and
the orientation of the (cylindrically symmetric) stomatocyte is given by
a vector, which points in the direction of the symmetry axis, it does
not need to be normalized. The parameters: outer\_radius , inner\_radius
, and layer\_width , specify the shape of the stomatocyte. Here
inappropriate choices of these parameters can yield undersired results.
The width is used as a scaling parameter. That is, a stomatocyte given
by :: = 7:3:1 is half the size of the stomatocyte given by 7:3:2. Not
all choices of the parameters give reasonable values for the shape of
the stomatocyte, but the combination 7:3:1 is a good point to start from
when trying to modify the shape.

In variant , a slit-shaped pore in a T-orientation to a flat channel is
created. The geometry is depicted in Fig. [fig:slitpore]. It
translationally invariant in y direction. The pore (lower vertical part)
extends in z-direction, and the channel (upper horizontal part). The
pore mouth is defined as the z-coordinate, where the lower plane of the
channel and the slit pore intersect. It is always centered in the
x-direction. A corresponding command decorates the surface with surface
charges that can be calculated with the ICC\ :math:`\star` algorithm.

[fig:slitpore]

.. figure:: figures/slitpore.pdf
   :alt: The slitpore created by the
   :align: center
   :height: 6.00000cm

   The slitpore created by the 

Variant specifies an electrostatic interaction between the charged
particles in the system to an infinitely long rod with a line charge of
which is alinge along the z-axis and centered at and .

Variant specifies the electrostatic interactinos between the charged
particles in the system and an inifinitely large plate in the x-y-plane
at height . The plate carries a charge density of .

Variant specifies the dipolar coupling of particles with a dipolar
moment to an external field .

Variant creates an infinite plane at a fixed position. For
non-initializing a direction of the constraint values of the positions
have to be negative. For the tunable-slip boundary interactions you have
to set *two* constraints.

Variant calculates the smallest distance to all non-penetrable
constraints, that can be repulsive (wall, cylinder, sphere, rhomboid,
maze, pore, slitpore). Negative distances mean that the position is
“within” the area that particles should not access. Helpful to find
initial configurations.)

Variant creates a hollow-cone shaped boundary. The position can be any
point in the simulation box, and the orientation of the (cylindrically
symmetric) cone is given by a vector, which points in the direction of
the symmetry axis, it does not need to be normalized. The parameters:
outer\_radius , inner\_radius , width , and opening\_angle , specify the
shape of the object. The inner radius gives the narrow end opening size,
the outer radius the length of the shaft, and the width the layer width,
i.e., the thickness of the cone. The opening\_angle (between 0 and
:math:`\pi/2`) specifies the angle of the cone.

Variant creates a spherocylinder, that is, a cylinder capped by a
hemisphere on either side. The parameter length specifies the length of
the shaft, excluding the two hemispherical caps.

To create an infinite plane in :math:`z`-direction at :math:`z=20.0` of
type id 1 which is directed inwards to the origin (0,0,0), use:

constraint wall normal 0 0 -1 dist -20 type 1

Python Syntax::

    import espressomd from espressomd.shapes import <SHAPE>
    system=espressomd.System()

``<SHAPE>`` can be any of:

* ::

     wall = Wall(normal=[n_x, n_y, n_z], dist=d)

creates a wall shape object with normal vector at distance from the
origin in the direction of the normal vector.

* ::

    sphere = Sphere(pos=[x, y, z], rad=R, direction=D)

    
creates a sphere object with its center at posistion and radius ``R``.

* ::

    cylinder = Cylinder(pos=[x, y, z], axis=[a\_x, a\_y, a\_z], rad=R,
    length=L, direction=D)

creates a cylinder object at position with its axis pointing torwards
   , radius and length .

* ::

    rhomboid = Rhomboid(pos=[x, y, z], a=[a\_x, a\_y, a\_z], b=[b\_x,
    b\_y, b\_z], c=[c\_x, c\_y, c\_z], direction=D)

creates a rhomboid defined by one corner located at and three
adjacent edges, defined by the three vectors connecting the corner p
with it’s three neighboring corners, a , b and c .

* ::

    maze = Maze(nsphere=NS, dim=Dim, cylrad=CR, sphrad=SR)

creates a -dimensional maze by spheres with radius per dimension. The
spheres build a grid of simple cubic symmetry and are connected by
cylinders with radius .

* ::

    pore = Pore(pos=[x, y, z], axis=[a\_x, a\_y, a\_z], length=L,
    smoothing\_radius = SR, rad\_left=LR, rad\_right=RR,
    outer\_rad\_left=ORL, outer\_rad\_right=ORR)

creates a pore at position , with its axis pointing in the direction
of and length .

* ::

    stomatocyte = Stomatocyte(position\_x=x, position\_y=y,
    position\_z=z, orientation\_x = o\_x, orientation\_y = o\_y,
    orientation\_z = o\_z, outer\_radius = OR, inner\_radius = IR,
    layer\_width = LW, direction = D)

creates a Stomatocyte at position with orientation whose outer radius
is , its inner radius is and the layer width will be .

* ::

    slitpore = Slitpore(pore\_mouth = z, channel\_width = CW, pore\_width
    = PW, pore\_length = PL, upper\_smoothing\_radius = USR,
    lower\_smoothing\_radius = LSR)

creates a Slitpore, the meaning of the geometrical parameters can be
inferred from fig. [fig:slitpore].

* ::

    spherocylinder=SpheroCylinder(pos=[x, y, z], axis=[a\_x, a\_y, a\_z],
    length=L, rad=R)

creates a Sphero-Cylinder at position , whose cylindrical element is
aligned in direction . The cylinder will have a length of and a
radius of . The spherical cap will have the same radius.

* ::

    hollowCone = HollowCone(position\_x = x, position\_y = y, position\_z
    = z, orientation\_x = o\_x, orientation\_y = o\_y, orientation\_z =
    o\_z, outer\_radius = OR, inner\_radius = IR, width = W,
    opening\_angle = a, direction = D)

creates a hollow cone whose axis is aligned to located at position ,
which inner and outer radii are and , respectively. The width and
opening angle are given by and .

The direction paramter for the shapes specifies wheter it will act
torwards the outside or inside . All those shapes can be used as
constraints and added to the systems constraints by passing a
initialized shape object to

::
    system.constraints.add(shape = shape\_object, particle\_type=p\_type)

The extra argument specifies the nonbonded interaction to be used with
that constraint. There are two further optional parameters and that can
be used to fine tune the behavior of the constraint. If penetrable is
set to then particles can move through the constraint in this case the
other option controls wheter the particle is subject to the interaction
potential of the wall. If set to then the constraint will only act in
the direction of the normal vector.

Deleting a constraint
~~~~~~~~~~~~~~~~~~~~~

constraint delete

This command will delete constraints. If is specified only this
constraint will deleted, otherwise all constraints will be removed from
the system.

Getting the force on a constraint
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

constraint force

Returns the force acting on the th constraint. Note, however, that this
are only forces due to interactions with particles, not with other
constraints. Also, these forces still do not mean that the constraints
move, they are just the negative of the sum of forces acting on all
particles due to this constraint. Similarly, the total energy does not
containt constraint-constraint contributions.

Getting the currently defined constraints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

constraint

Prints out all constraint information. If is specified only this
constraint is displayed, otherwise all constraints will be printed.

``harmonic_well``: Creating a harmonic trap
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

harmonic\_well { }

Calculates a spring force for all particles, where the equilibrium
position of the spring is at and it’s force constant is . A more
flexible trap can be constructed with constraints, but this one runs on
the GPU.

Virtual sites
-------------

Virtual sites are particles, the positions and velocities of which are
not obtained by integrating an equation of motion. Rather, their
coordinates are obtained from the position (and orientation) of one or
more other particles. In this way, rigid arrangements of particles can
be constructed and a particle can be placed in the center of mass of a
set of other particles. Virtual sites can interact with other particles
in the system by means of interactions. Forces are added to them
according to their respective particle type. Before the next integration
step, the forces accumulated on a virtual site are distributed back to
those particles, from which the virtual site was derived.

There are two distinct types of virtual sites, described in the
following.

Virtual sites in the center of mass of a molecule
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To activate this implementation, enable the feature (sec.
[sec:myconfig]). Virtual sites are then placed in the center of mass of
a set of particles (as defined below). Their velocity will also be that
of the center of mass. Forces accumulating on the virtual sites are
distributed back to the particles which form the molecule. To place a
virtual site at the center of a molecule, perform the following steps in
that order

#. Create a particle of the desired type for each molecule. It should be
   placed at least roughly in the center of the molecule to make sure,
   it’s on the same node as the other particles forming the molecule, in
   a simulation with more than one cpu.

#. Make it a virtual site using

   part virtual 1

#. Declare the list of molecules and the particles they consist of:

   analyze set { ...} ...

   The lists of particles in a molecule comprise the non-virtual
   particles as well as the virtual site. The id of this molecule is its
   index in this list. For example,

   analyze set {0 1 2 3 4} {0 5 6 7 8} {1 9 10 11}

   declares three molecules, of which the first two consist of three
   particles and a virtual site each (particles 1–4 and 5–8,
   respectively). The third molecule has type 1 and consists of two
   particles and a virtual site. The virtual sites were determined
   before by setting the flag. You can choose freely one out of each
   molecule, for example particles 1, 5, and 9.

#. Assign to all particles that belong to the same molecule the
   molecule’s id

   part mol

   The molid is the index of the particle in the above list, so you
   would assign 0 to particles 1-4, 1 to particles 5-8 and 2 to
   particles 9-11. Alternatively, you can call

   analyze set topo\_part\_sync

   to set the s from the molecule declarations.

#. Update the position of all virtual particles (optional)

   integrate 0

Please note that the use of virtual sites requires that the particles
are numbered consecutively. I.e., the particle ids should go from zero
to :math:`N-1`, where :math:`N` is the number of particles.

The type of the molecule you can choose freely, it is only used in
certain analysis functions, namely ``energy_kinetic_mol``,
``pressure_mol`` and ``dipmom_mol``, which compute kinetic energy,
pressure and dipole moment per molecule type, respectively.

Rigid arrangements of particles
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The “relative” implementation of virtual sites allows for the simulation
of rigid arrangements of particles. It can be used, , for extended
dipoles and raspberry-particles, but also for more complex
configurations. Position and velocity of a virtual site are obtained
from the position and orientation of exactly one non-virtual particle,
which has to be placed in the center of mass of the rigid body. Several
virtual sites can be related to one and the same non-virtual particle.
The position of the virtual site is given by

.. math:: \vec{x_v} =\vec{x_n} +O_n (O_v \vec{E_z}) d,

where :math:`\vec{x_n}` is the position of the non-virtual particle,
:math:`O_n` is the orientation of the non-virtual particle, :math:`O_v`
denotes the orientation of the vector :math:`\vec{x_v}-\vec{x_n}` with
respect to the non-virtual particle’s body fixed frame and :math:`d` the
distance between virtual and non-virtual particle. In words: The virtual
site is placed at a fixed distance from the non-virtual particle. When
the non-virtual particle rotates, the virtual sites rotates on an orbit
around the non-virtual particle’s center.

To use this implementation of virtual sites, activate the feature (see
sec. [sec:myconfig]). To set up a virtual site,

#. Place the particle to which the virtual site should be related. It
   needs to be in the center of mass of the rigid arrangement of
   particles you create. Let its particle id be n.

#. Place a particle at the desired relative position, make it virtual
   and relate it to the first particle

   part pos virtual 1 vs\_auto\_relate

#. Repeat the previous step with more virtual sites, if desired.

#. To update the positions of all virtual sites, call

   integrate 0

Please note:

-  The relative position of the virtual site is defined by its distance
   from the non-virtual particle, the id of the non-virtual particle and
   a quaternion which defines the vector from non-virtual particle to
   virtual site in the non-virtual particle’s body-fixed frame. This
   information is saved in the virtual site’s vs\_relative-attribute.
   Take care, not to overwrite these after using vs\_auto\_relate.

-  Virtual sites can not be placed relative to other virtual sites, as
   the order in which the positions of virtual sites are updated is not
   guaranteed. Always relate a virtual site to a non-virtual particle
   placed in the center of mass of the rigid arrangement of particles.

-  Don’t forget to declare the particle virtual in addition to calling
   vs\_auto\_relate

-  In case you know the correct quaternions, you can also setup a
   virtual site using

   part virtual 1 vs\_relative

   where n is the id of the non-virtual particle, d is its distance from
   the virtual site, and q are the quaternions.

-  In a simulation on more than one CPU, the effective cell size needs
   to be larger than the largest distance between a non-virtual particle
   and its associated virtual sites. To this aim, you need to set the
   global variable to this largest distance. issues a warning when
   creating a virtual site with and the cutoff is insufficient.

-  If the virtual sites represent actual particles carrying a mass, the
   inertia tensor of the non-virtual particle in the center of mass
   needs to be adapted.

-  The presence of rigid bodies constructed by means of virtual sites
   adds a contribution to the pressure and stress tensor.

-  The use of virtual sites requires that the particles are numbered
   consecutively, , the particle ids should go from zero to :math:`N-1`,
   where :math:`N` is the number of particles.

Additional features
~~~~~~~~~~~~~~~~~~~

The behaviour of virtual sites can be fine-tuned with the following
switches in ``myconfig.hpp`` (sec. [sec:myconfig])

-  specifies that the velocity of virtual sites is not computed

-  specifies that the Langevin thermostat should also act on virtual
   sites

-  specifies that the thermostat does not act on non-virtual particles

Grand canonical feature
-----------------------

For using conveniently for simulations in the grand canonical ensemble,
or other purposes, when particles of certain types are created and
deleted frequently. Particle ids can be stored in lists for each
individual type and so random ids of particles of a certain type can be
drawn.

from espressomd import grand\_canonical grand\_canonical.setup([\_type])
grand\_canonical.delete\_particles(\_type)
grand\_canonical.find\_particle(\_type)
grand\_canonical.number\_of\_particles(\_type)

If you want to keep track of particle ids of a certain type you have to
initialize the method by calling

part gc

grand\_canonical.setup([\_type])

After that will keep track of particle ids of that type. When using the
keyword ``find`` and a particle type, the command will return a randomly
chosen particle id, for a particle of the given type. The keyword
``status`` will return a list with all particles with the given type,
similarly giving ``number`` as argument will return the number of
particles which share the given type.
