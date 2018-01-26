.. _Advanced Methods:

Advanced Methods
================

.. todo:: Write short introduction

.. _Creating bonds when particles collide:

Creating bonds when particles collide
-------------------------------------

Please cite :cite:`espresso2` when using dynamic bonding.

With the help of this feature, bonds between particles can be created
automatically during the simulation, every time two particles collide.
This is useful for simulations of chemical reactions and irreversible
adhesion processes. Both, sliding and non-sliding contacts can be created.

The collision detection is controlled via the :attr:`espressomd.system.System.collision_detection` attribute, which is an instance of the class :class:`espressomd.collision_detection.CollisionDetection`.

Several modes are available for different types of binding.

* "bind_centers": adds a pair-bond between two particles at their first collision. By making the bonded interaction `stiff` enough, the particles can be held together after the collision. Note that the particles can still slide on each others' surface, as the pair bond is not directional. This mode is set up as follows::
    import espressomd
    from espressomd.interactions import HarmonicBond
    
    system=espressomd.System()
    bond_centers=HarmonicBond(k=1000,r_0=<CUTOFF>)
    system.bonded_inter.add(bond_centers)
    system.collision_detection.set_params(mode="bind_centers",distance=<CUTOFF>, bond_centers=bond_centers)
  
  The parameters are as follows:
  
    * `distance` is the distance between two particles at which the binding is triggered. This cutoff distance, `<CUTOFF>` in the example above, is typically chosen slightly larger than the particle diameter. It is also a good choice for the equilibrium length of the bond.
    * `bond_centers` is the bonded interaction (an instance of :class:espressomd.interactions.BondedInteraction`) to be created between the particles. No guarantees are made regarding which of the two colliding particles gets the bond. Once there is a bond of this type on any of the colliding particles, no further binding occurs for this pair of particles.

* "bind_at_point_of_collision": this mode prevents sliding of the colliding particles at the contact. This is achieved by
  creating two virtual sites at the point of collision. They are
  rigidly connected to the colliding particles, respectively. A bond is
  then created between the virtual sites, or an angular bond between
  the two colliding particles and the virtual particles. In the latter case,
  the virtual particles are the centers of the angle potentials
  (particle 2 in the description of the angle potential (see :ref:`Bond-angle interactions`).
  Due to the rigid connection between each of the
  particles in the collision and its respective virtual site, a sliding
  at the contact point is no longer possible. See the documentation on
  :ref:`Rigid arrangements of particles` for details. In addition to the bond between the virtual
  sites, the bond between the colliding particles is also created, i.e., the "bind_at_point_of_collision" mode implicitly includes the "bind_centers" mode. You
  can either use a real bonded interaction to prevent wobbling around
  the point of contact or you can use :class:`espressomd.interactions.Virtual` which acts as a marker, only.
  The method is setup as follows::
     
     system.collision_detection.set_params(mode="bind_at_point_of_collision", distance=<CUTOFF>, bond_centers=<BOND_CENTERS>, bond_vs=<BOND_VS>, part_type_vs=<PART_TYPE_VS>, vs_placement=<VS_PLACEMENT>)

  
  The parameters `distance` and `bond_centers` have the same meaning as in the `bind_centers` mode. The remaining parameters are as follows:
    
    * `bond_vs` is the bond to be added between the two virtual sites created on collision. This is either a pair-bond with an equilibrium length matching the distance between the virtual sites, or an angle bond fully stretched in its equilibrium configuration.
    * `part_type_vs` is the particle type assigned to the virtual sites created on collision. In nearly all cases, no non-bonded interactions should be defined for this particle type.
    * `vs_placement` controls, where on the line connecting the centers of the colliding particles, the virtual sites are placed. A value of 0 means that the virtual sites are placed at the same position as the colliding particles on which they are based. A value of 0.5 will result in the virtual sites being placed ad the mid-point between the two colliding particles. A value of 1 will result the virtual site associated to the first colliding particle to be placed at the position of the second colliding particle. In most cases, 0.5, is a good choice. Then, the bond connecting the virtual sites should have an equilibrium length of zero.

* "glue_to_surface": This mode is used to irreversibly attach small particles to the surface of a big particle. It is asymmetric in that several small particles can be bound to a big particle but not vice versa. The small particles can change type after collision to make them `inert`. On collision, a single virtual site is placed and related to the big particle. Then, a bond (`bond_centers`) connects the big and the small particle. A second bond (`bond_vs`) connects the virtual site and the small particle. Further required parameters are:
  
  * `part_type_to_attach_vs_to`: Type of the particle to which the virtual site is attached, i.e., the `big` particle.
  * `part_type_to_be_glued`: Type of the particle bound to the virtual site (the `small` particle).
  * `part_type_after_glueing`: The type assigned to the particle bound to the virtual site (`small` particle) after the collision.
  * `part_type_vs`: Particle type assigned to the virtual site created during the collision. 
  * `distance_glued_particle_to_vs`: Distance of the virtual site to the particle being bound to it (`small` particle).




- "bind_three_particles" allows for the creation of agglomerates which maintain their shape
  similarly to those create by the mode "bind_at_point_of_collision". The present approach works
  without virtual sites. Instead, for each two-particle collision, the
  surrounding is searched for a third particle. If one is found,
  angular bonds are placed to maintain the local shape.
  If all three particles are within the cutoff distance, an angle bond is added
  on each of the three particles in addition
  to the distance based bonds between the particle centers. 
  If two particles are within the cutoff of a centrla particle (e.g., chain of three particles)
  an angle bond is placed on the central particle.
  The angular bonds being added are determined from the angle between the particles.
  This method does not depend on the particles’ rotational
  degrees of freedom being integrated. Virtual sites are also not
  required.
  The method, along with the corresponding bonds are setup as follows::
        
        n_anlge_bonds=181 # 0 to 180 degrees in one degree steps
        for i in range(0,res,1):
           self.s.bonded_inter[i]=Angle_Harmonic(bend=1,phi0=float(i)/(res-1)*np.pi)
        
        # Create the bond passed to bond_centers here and add it to the system
        
        self.s.collision_detection.set_params(mode="bind_three_particles",bond_centers=<BOND_CENTERS>,bond_three_particles=0,three_particle_binding_angle_resolution=res,distance=<CUTOFF>)

  Important: The bonds for the angles are mapped via their numerical bond ids. In this example, ids from 0 to 180 are used. All other bonds required for the simulation need to be added to the system after those bonds. In particular, this applies to the bonded interaction passed via `bond_centers` 


The following limitations currently apply for the collision detection:
* No distinction is currently made between different particle types for the `bind_centers` method.
* The “bind at point of collision” and "glue to surface"  approaches require the feature `VIRTUAL_SITES_RELATIVE` to be activated in `myconfig.hpp`.

* The “bind at point of collision” approach cannot handle collisions
  between virtual sites

.. _Catalytic Reactions:

Catalytic Reactions
-------------------


:ref:`Catalytic Reaction theory`

With the help of the feature ``CATALYTIC_REACTIONS``, one can define three particle types to act as reactant (e.g. :math:`\mathrm{H_2 O_2}`), catalyzer (e.g. platinum), and product (e.g. :math:`\mathrm{O_2}` and :math:`\mathrm{H_2 O}`). The current setup allows one to simulate active swimmers and their chemical propulsion.

For a Janus swimmer consisting of platinum on one hemisphere and gold on the other hemisphere, both surfaces catalytically induce a reaction. We assume an initial abundance of hydrogen peroxide and absence of products, so that back (recombination) reactions seldomly occur at the surface. A typical model for the propulsion of such a particle assumes

.. math::

    \begin{aligned}
      \mathrm{H_2 O_2} &\xrightarrow{\text{Pt}} \mathrm{2 H^{+} + 2 e^{-} + O_2} \\
      \mathrm{2 H^{+} + 2 e^{-} + H_2 O_2} &\xrightarrow{\text{Au}} \mathrm{2 H_2 O}
    \end{aligned}

That is, catalytic surfaces induce a reactions that produce charged species by consuming hydrogen peroxide. It is the change in distribution of charged species that leads to motion of the swimmer, a process referred to as self-electrophoresis. A minimal model for this would be

.. math::

    \begin{aligned}
      A &\xrightarrow{C^{+}} B \\
      B &\xrightarrow{C^{-}} A
    \end{aligned}

where on the upper half of the catalyst :math:`C^{+}` a species :math:`A` is converted into :math:`B`, and on the lower half :math:`C^{-}` the opposite reaction takes place. Note that when :math:`A` and :math:`B` are charged, this reaction conserves charge, provided the rates are equal.

In |es| the orientation of a catalyzer particle is used to define hemispheres; half spaces going through the particle's center. The reaction region is bounded by the *reaction range*: :math:`r`. Inside the reaction range, we react only rectant-product pairs. The particles in a pair are swapped from hemisphere to another with a rate prescribed by

.. math::

    P_{\text{move}} = 1 - \mathrm{e}^{-k_{\mathrm{ct}}\,\Delta t} ,

with the reaction rate :math:`k_{\mathrm{ct}}` and the simulation time step :math:`\Delta t`. A pair may be swapped only once per MD time step, to avoid a no-net-effect situation. That is, we allow an exchange move only when the following conditions are met:

1. Both partners of the reactant-product pair have to reside within the reaction range.
2. The product has to reside in the upper half-space of the reaction range.
3. The reactant has to reside in the lower half-space of the reaction range.

Self-propulsion is achieved by imposing an interaction asymmetry between the partners of a swapped pair. That is, the heterogeneous distribution of chemical species induced by the swapping leads to a net force on the particle, counter balanced by friction.

To set up the system for catalytic reactions the class :class:`espressomd.reaction.Reaction`
can be used.::

    from espressomd.reaction import Reaction

    system = espressomd.System()

    # setting up particles etc

    r = Reaction(product_type=1, reactant_type=2, catalyzer_type=0, ct_range=2, ct_rate=0.2, eq_rate=0)
    r.start()
    r.stop()

    print r

* the first invocation of ``Reaction``, in the above example,  defines a
  reaction with particles of type number 2 as reactant, type 0 as catalyzer and
  type 1 as product [#1]_. The catalytic reaction rate constant is given by :math:`\mathrm{ct\_rate}`
  [#2]_ and to override the default rate constant for the equilibrium reaction
  ( = 0), one can specify it by as ``eq_rata``.  By default each reactant particle is checked
  against each catalyst particle (``react_once=False``). However, when creating
  smooth surfaces using many catalyst particles, it can be desirable to let the
  reaction rate be independent of the surface density of these particles. That
  is, each particle has a likelihood of reacting in the vicinity of the surface
  (distance is less than :math:`r`) as specified by the rate constant, i.e.,
  *not* according to :math:`P_{\text{cvt}} = 1 - \exp \left( - n k\Delta t
  \right)`, with :math:`n` the number of local catalysts. To accomplish this,
  each reactant is considered only once each time step by using the option
  ``react_once=True`` . The reaction command is set up such that the different
  properties may be influenced individually.

*  ``r.stop()`` disables the reaction. Note that at the moment, there can
   only be one reaction in the simulation.

*  ``print r``  returns the current reaction parameters.

In future versions of |es| the capabilities of the ``CATALYTIC_REACTIONS`` feature may be generalized
to handle multiple reactant, catalyzer, and product types, as well as
more general reaction schemes. Other changes may involve merging the
current implementation with the ``COLLISION_DETECTION`` feature.

.. rubric:: Footnotes

.. [#1]
   Only one type of particle can be assigned to each of these three
   reaction species and no particle type may be assigned to multiple
   species. That is, currently does not support particles of type 1 and
   2 both to be reactants, nor can particles of type 1 be a reactant as
   well as a catalyst. Moreover, only one of these reactions can be
   implemented in a single Tcl script. If, for instance, there is a
   reaction involving particle types 1, 2, and 4, there cannot be a
   second reaction involving particles of type 5, 6, and 8. It is
   however possible to modify the reaction properties for a given set of
   types during the simulation.

.. [#2]
   Currently only strictly positive values of the catalytic conversion
   rate constant are allowed. Setting the value to zero is equivalent to
   ``r.stop()``.

..
    .. _\`\`nemd\`\`\: Setting up non-equilibirum MD:

    ``nemd``: Setting up non-equilibrium MD
    ---------------------------------------

    .. todo::
        This is not implemented for the python interface yet

    nemd exchange nemd shearrate nemd off nemd nemd profile nemd viscosity

    Use NEMD (Non Equilibrium Molecular Dynamics) to simulate a system under
    shear with help of an unphysical momentum change in two slabs in the
    system.

    Variants and will initialise NEMD. Two distinct methods exist. Both
    methods divide the simulation box into slabs that lie parallel to the
    x-y-plane and apply a shear in x direction. The shear is applied in the
    top and the middle slabs. Note, that the methods should be used with a
    DPD thermostat or in an NVE ensemble. Furthermore, you should not use
    other special features like or inside the top and middle slabs. For
    further reference on how NEMD is implemented into see
    :cite:`soddeman01a`.

    Variant chooses the momentum exchange method. In this method, in each
    step the largest positive x-components of the velocity in the middle
    slab are selected and exchanged with the largest negative x-components
    of the velocity in the top slab.

    Variant chooses the shear-rate method. In this method, the targeted
    x-component of the mean velocity in the top and middle slabs are given
    by

    .. math:: {target\_velocity} = \pm {shearrate}\,\frac{L_z}{4}

    where :math:`L_z` is the simulation box size in z-direction. During the
    integration, the x-component of the mean velocities of the top and
    middle slabs are measured. Then, the difference between the mean
    x-velocities and the target x-velocities are added to the x-component of
    the velocities of the particles in the respective slabs.

    Variant will turn off NEMD, variant will print usage information of the
    parameters of NEMD. Variant will return the velocity profile of the
    system in x-direction (mean velocity per slab).

    Variant will return the viscosity of the system, that is computed via

    .. math:: \eta = \frac{F}{\dot{\gamma} L_x L_y}

    where :math:`F` is the mean force (momentum transfer per unit time)
    acting on the slab, :math:`L_x L_y` is the area of the slab and
    :math:`\dot{\gamma}` is the shearrate.

    NEMD as implemented generates a Poiseuille flow, with shear flow rate
    varying over a finite wavelength determined by the box. For a planar
    Couette flow (constant shear, infinite wavelength), consider using
    Lees-Edwards boundary conditions (see ) to drive the shear.

.. _Lees-Edwards boundary conditions:

Lees-Edwards boundary conditions
--------------------------------

To use the Lees-Edwards boundary conditions, the feature ``LEES_EDWARDS`` is required.

Lees-Edwards boundary conditions can be used to introduce a shear flow to the MD simulation. An introduction can be found in :cite:`lees72`. Compared to NEMD simulations they have two big advantages: First, the bulk behavior of the system remains unchanged. Second, the image boxes are moved, whereas the flow within the primary simulation box has to develop on its own. Hence, this allows two additional phenomena: Shear banding can occur as well as non-linear shear profiles can be observed. This makes Lees-Edwards boundary conditions suitable for comparison with rheological experiments. 

Lees-Edwards boundary conditions impose a shear flow of speed :math:`\dot\gamma` by moving the periodic image boxes along the x-direction according to:

.. math:: v_{\text{x, unfolded}} = v_{\text{x, folded}} + \dot\gamma \cdot y_{\text{imagecount}}

:math:`v_{\text{x, unfolded}}` refers to the velocity of a particle outside the main simulation box, :math:`y_{\text{imagecount}}` is the amount of periodic boundaries crossed in the  :math:`y`-direction. 

The absolute offset of the periodic images can be set via

* :py:attr:`~espressomd.System().lees_edwards_offset`

The following example introduces the usage::
    
    import espressomd
    system = espressomd.System()
    absolute_offset = 0.2
    system.lees_edwards_offset = absolute_offset

Lees-Edwards boundary conditions can be used to obtain the shear modulus :math:`G = \frac{\tau}{\gamma}` or the shear viscosity :math:`\eta = \frac{\tau}{\dot\gamma}` outside the linear regime, where Green-Kubo relations are not valid anymore. For this purpose a lees_edwards_offset is set followed by one integration step for multiple times. Strain, strain rate and the shear stress need to be recorded for the calculation. Alternatively a sinusoidal lees_edwards_offset series can be used to carry out oscillatory experiments to calculate viscoelastic moduli (:math:`G', G''`). Furthermore a lees_edwards_offset can be set followed by many integration steps obtain the relaxation behaviour of a system. 

When applying a constant shear rate :math:`\dot\gamma` the velocity of the particles changes from :math:`-\frac{\dot\gamma}{2}` at the bottom of the box to :math:`\frac{\dot\gamma}{2}` at the top of the box. 

Physical meaningful values for systems where hydrodynamics play a major role, can only be obtained by including hydrodynamic interactions. Lees-Edwards boundary conditions are implemented in the :ref:`Lattice-Boltzmann` algorithms. For this algorithm the feature ``LB_GPU`` is required. Please refer to chapter :ref:`Lattice-Boltzmann` for more information. 

Lees-Edwards boundary conditions work with the DPD thermostat. In order to correctly observe transport properties, symmetry-breaking or entropy production in relation to shear flow is probably better to use the DPD thermostat (:ref:`Dissipative Particle Dynamics (DPD)`) once the initial heat-up has been carried out. The DPD thermostat removes kinetic energy from the system based on a frictional term defined relative to a local reference frame of a given particle-pair, without enforcing any specific flow pattern apriori. At high rates of dissipation, this can however lead to an artefactual shear-banding type effect at the periodic boundaries, such that the bulk fluid is nearly stationary. y. This effect is removed using the modification proposed to the DPD thermostat by Chatterjee :cite:`chatterjee2007` to allow treatment of systems with high dissipation rates, which is applied automatically if ``LEES_EDWARDS`` is compiled in. Chatterjee’s modification is just to skip calculation of DPD forces (both dissipative and random) for particle pairs which cross a boundary in y.

The command::

  print(system.lees_edwards_offset)

returns the current value of the offset. If ``LEES_EDWARDS`` is compiled in, then coordinates are folded into the primary simulation box as the integration progresses, to prevent a numerical overflow.

.. _Immersed Boundary Method for soft elastic objects:

Immersed Boundary Method for soft elastic objects
-------------------------------------------------

Please contact the Biofluid Simulation and Modeling Group at the
University of Bayreuth if you plan to use this feature.

This section describes an alternative way to include soft elastic
objects somewhat different from the previous chapter. In the Immersed
Boundary Method (IBM), soft particles are considered as an infinitely
thin shell filled with liquid (see e.g.
 :cite:`Peskin2002,Crowl2010,KruegerThesis`). When the
shell is deformed by an external flow it responds by elastic restoring
forces which are transmitted into the fluid. In the present case, the
inner and outer liquid are of the same type and are simulated using
Lattice-Boltzmann.

Numerically, the shell is discretized by a set of marker points
connected by triangles. The marker points are advected with *exactly*
the local fluid velocity, i.e., they do not possess a mass nor a
friction coefficient (this is different from the Object-in-Fluid method
of the previous chapter). We implement these marker points as virtual
particles in which are not integrated using the usual velocity-verlet
scheme, but instead are propagated using a simple Euler algorithm with
the local fluid velocity (if the ``IMMERSED_BOUNDARY`` feature is turned
on).

To compute the elastic forces, three new bonded interactions are defined
ibm\_triel, ibm\_tribend and ibm\_volCons:

-  ibm\_triel is a discretized elastic force with the following syntax

   inter ibm\_triel

   where , and represent the indices of the three marker points making
   up the triangle. The parameter specifies the maximum stretch above
   which the bond is considered broken. The final parameter can be
   either

   ::

       NeoHookean <k>

   or

   ::

       Skalak <k1> <k2>

   which specifies the elastic law and its corresponding parameters (see
   e.g. :cite:`KruegerThesis`).

-  ibm\_tribend is a discretized bending potential with the following
   syntax

   inter ibm\_tribend

   where , , and are four marker points corresponding to two neighboring
   triangles. The indices and contain the shared edge. Note that the
   marker points within a triangle must be labelled such that the normal
   vector
   :math:`\vec{n} = (\vec{r}_\text{ind2} - \vec{r}_\text{ind1}) \times (\vec{r}_\text{ind3} - \vec{r}_\text{ind1})`
   points outward of the elastic object.

   The parameter allows to specify different numerical ways of computing
   the bending interaction. Currently, two methods are implemented,
   where the first one () follows :cite:`KruegerThesis` and
   the second one () follows :cite:`Gompper1996`. In both
   cases, is the bending modulus. The options or specify whether the
   reference shape is a flat configuration or whether the initial
   configuration is taken as reference shape, this option is only
   available for the method.

-  ibm\_volCons is a volume-conservation force. Without this correction,
   the volume of the soft object tends to shrink over time due to
   numerical inaccuracies. Therefore, this implements an artificial
   force intended to keep the volume constant. If volume conservation is
   to be used for a given soft particle, the interaction must be added
   to every marker point belonging to that object. The syntax is

   inter ibm\_volCons

   where identifies the soft particle and is a volumetric spring
   constant :cite:`KruegerThesis`.


.. _Object-in-fluid:

Object-in-fluid
---------------

Please cite  if you use the object-in-fluid implementation described
below. For more details also see the documentation at
http://cell-in-fluid.fri.uniza.sk/oif-documentation or contact the
Cell-in-fluid Research Group at University of Žilina.

Simulations using work mostly with objects (molecules, atoms, polymers,
colloids, crystals, …) that are physically composed of points linked
together with bonds. These objects are like skeletons, without inner or
outer volume.

The idea behind this module, is to use for objects that do have inner
volume, for example blood cells, magnetic beads, capsules, …The boundary
of an object is covered with triangular mesh. The vertices of the mesh
are declared in as particles. The edges of the mesh define elastic
forces keeping the shape of the object. The movement of object is
achieved by adding forces to the mesh points.

Modelled elastic or rigid objects are immersed in the LB fluid flow. The
fluid interacts with an elastic object resulting in its deformation;
this immediately generates forces acting back on the fluid. The aim is
to describe the immersed object using the notion of particles, and to
create bonds between these particles representing elastic or rigid
forces.

The objects are composed of a membrane encapsulating the fluid inside
the object. For now, the inside fluid must have the same density and
viscosity as the outside fluid. The object is represented by its
membrane (boundary), that is discretized using a triangulation. Such
triangulation defines interacting particles distributed on the surface
of the immersed object :cite:`dupin07`:

-  between two particles, corresponding to the edges in the
   triangulation (modelling the stretching of the membrane),

-  between three particles, corresponding to the triangles of the
   triangulation (local area, or local surface preservation of the
   membrane),

-  between four particles, corresponding to two triangles from the
   triangulation sharing a common edge (bending of the membrane).

The object immersed in the fluid moves under the influence of the
deforming forces, defined through the bonds, and under the influence of
the fluid motion. This interaction is based on the frictional force
between the fluid and the surface particles. Therefore the object moves
in the flow only if there is a nonzero difference between the fluid
velocity and the particle velocity. In other words, there has to be at
least small flow through the membrane, which is in most cases
unphysical. However, this unphysical flow through the membrane is
probably negligible in larger scales.

.. _Membranes:

Membranes
~~~~~~~~~

With this approach, it is easy to model also elastic sheets, or free
membranes that do not necessarily enclose a 3D object. In this case,
area\_force\_global and volume\_force interactions are not needed, since
these two interactions are meant for closed immersed objects.

.. _Parameters:

Parameters
~~~~~~~~~~

There are several parameters involved in this model. All of them should
be calibrated according to the intended application.

-  Mass of the particles. Every particle has its mass, which influences
   the dynamics.

-  Friction coefficient. The main parameter describing the
   fluid-particle interaction is the ``riction \ parameter ``\ rom the
   command ``bf``\ uid .

-  Parameters of elastic moduli. Elastic behaviour can be described by
   five different elastic moduli: hyperelastic stretching, linear
   stretching, bending, local and global area preservation and volume
   preservation. Each of them has its own scaling parameter:
   :math:`k_s, ks_{lin}, k_b, k_{al}, k_{ag}, k_v`. Their mathematical
   formulations have been taken from :cite:`dupin07`.

The mass of the particles and the friction coefficient can be calibrated
using the drag coefficients of the ellipsoidal objects. These drag
coefficients have known analytical values and the mass and friction can
be calibrated to fit this values. More details about the calibration can
be found in :cite:`cimrak`.

The elastic parameters are specific to the immersed objects. They
correspond to their physical values. More details about their mechanical
and biological meaning is presented in :cite:`dao03`
specifically for red blood cells. However, the proper calibration to fit
the experimental data has been performed in :cite:`cimrak`.

.. _Geometry:

Geometry
~~~~~~~~

The membrane of the immersed object is triangulated. In
doc/tutorials/03-object\_in\_fluid you can find an example using
deformable objects in the fluid.

|image|

Triangulation can be obtained using various software tools. Two files
are needed for mesh input:
``mesh-nodes.dat`` and ``mesh-triangles.dat``. The parameters of
the mesh are the number of particles on the surface of the immersed
object, denoted by ``mesh_nnode``, and the number of triangular faces
in the triangulation, denoted by ``mesh_ntriangle``. These parameters
are obtained automatically from ``mesh-nodes.dat`` and ``mesh-triangles.dat`` 
by counting the number of lines in respective files.

The ``mesh-nodes.dat`` thus contains ``mesh_nnode`` lines with three
real numbers separated by blank space, representing three coordinates of
the corresponding particle. The membrane is thus discretized into
``mesh_nnode`` particles with IDs starting from 0 to ``mesh_nnode-1``.
The IDs are assigned in the same order as in the ``mesh-nodes.dat``
file.

The ``mesh-triangles.dat`` contains ``mesh_ntriangle`` lines with three
nonnegative integers separated by blank space. Each line represents one
triangle in the triangulation. For algorithmic purposes it is crucial to
have defined a correct orientation of the triangle. The orientation is
defined using the normal vector associated with the triangle. The
important rule is that the normal vector of the triangle must point
inside the immersed object.

As an example, let us have one line in the file ``mesh-triangles.dat``
with numbers 4, 0 and 7. This means that particles with IDs 4, 0 and 7
form one triangular face of the triangulation. The orientation is
defined as follows: create two vectors :math:`v_1` and :math:`v_2`, such
that :math:`v_1` is pointing from particle 4 to particle 0, and
:math:`v_2` is pointing from particle 4 to particle 7. Be careful, the
order of vectors and particles matters!

The normal vector :math:`n` is computed as a vector product
:math:`v_1 \times v_2`. The direction of :math:`n` can be determined by
the rule of right hand: the thumb points in the :math:`v_1` direction,
the index finger in the :math:`v_2` direction and the middle finger in
the :math:`n` direction. Following this principle, all the lines in the
``mesh-triangles.dat`` files must be such that the normal vectors of the
corresponding triangles points inside the immersed object.

These two files are sufficient to describe the geometry and topology of
the triangulation. The following geometric entities are necessary for
the definition of bonded interactions: position of the particles, edges,
lengths of the edges, triangles, areas of triangles, angles between two
triangles sharing a common edge, surface of the immersed object, volume
of the immersed object. All these geometrical entities can be computed
using the information from the files ``mesh-nodes.dat`` and
``mesh-triangles.dat`` and the computation is done in the script
``scripts/object_in_fluid.tcl`` .

The script ``scripts/object_in_fluid.tcl`` reads both mesh files,
generates list of edges, and computes all geometrical entities needed
for definition of bonded interactions. It then executes commands
creating the particles, interactions and bonds.An example of ``part``
command is as follows:

::

     part 0 pos 3.0 3.0 6.0 type 1 mol 1 mass 1 

Note, the is feature ``mol`` that used for the particles. We use this
feature we distinguish between different objects. The upper limit for
the number of objects is 10000. However it can be increased by changing
the ``MAX_OBJECTS_IN_FLUID`` constant.

The following example shows an interaction.

::

     inter 106 oif_local_force 1.0 0.5 0.0 1.7 0.6 0.2 0.3 1.1 

This command (“invisible” for the user who executes the
``cript``/object\_in\_fluid.tcl  script) takes care of stretching,
bending and local area conservation all in one interaction with ID 106.
Detailed description of the available types of interactions is presented
in Section [sec:inter-bonded-oif].

.. _Available commands:

Available commands
~~~~~~~~~~~~~~~~~~

In order to use the object-in-fluid (OIF) commands and work with
immersed objects, the following features need to be compiled in:
``ASS, \ \verb EXTERNAL_FORCES . We do not specifically require \verb LB, \ \verb LB_BOUNDARIES, \ \verb CONSTRAINTS, \ \verb SOFT_SPHERE, \ \verb ``\ EMBRANE\_COLLISION,
 ``IF_L``\ CAL\_FORCES,  ``IF_GL``\ BAL\_FORCES.  They are most likely
to be used (for objects immersed in fluid and interacting with
boundaries and each other), but they are not necessary for the following
commands. For up-to-date overview of available oif commands see the OIF
user guide at cell-in-fluid.fri.uniza.sk/oif-documentation.

.. _Initialisation:

Initialisation
^^^^^^^^^^^^^^

oif\_init

Must be used before any other OIF command, initializes all global
variables and lists, does not take any arguments.

.. _Information about object-in-fluid structures:

Information about object-in-fluid structures
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

oif\_info

Prints information about whole framework, all global variables,
currently available templates and objects, etc. Does not take any
arguments.

.. _Templates for objects:

Templates for objects
^^^^^^^^^^^^^^^^^^^^^

template-id nodes-file triangles-file

This command creates a template that will be used for all objects that
share the same elastic properties and have the same triangulation.

specifies a unique ID for each template. The first template has the ID
0. The following ones need to be be numbered consecutively.

input file, each line contains three real numbers. These are the
:math:`x, y, z` coordinates of individual surface mesh nodes of the
objects.

input file, each line contains three integers. These are the ID numbers
of the mesh nodes as they appear in . Note, the first node has ID 0.

coefficients by which the coordinates stored in will be stretched in the
:math:`x, y, z` direction. The default values are 1.0 1.0 1.0.

| whether the respective coordinate will be flipped around 0.
  Coefficients :math:`x, y, z` must be either 0 or 1. The reflection of
  only one coordinate is allowed so at most one number is 1, others are
  0. For example ``mirror`` *0 1 0* results in flipping the spatial
  point :math:`(x,y,z)` to :math:`(x,-y,z)`. The default value is 0 0 0.

elastic modulus for hyperelastic stretching forces

elastic modulus for linear stretching forces

elastic modulus for bending forces

elastic modulus for local area forces

elastic modulus for global area forces

elastic modulus for volume forces

switch to turn on the computation of local outward normal vectors

The four switches ``ks``, ``kslin``, ``kb`` and ``kal`` set elastic
parameters for local interactions - ``ks`` for hyperelastic edge
stiffness, ``kslin`` for linear edge stiffness, ``kb`` for angle
preservation stiffness and ``kal`` for triangle surface preservation
stiffness. This stiffness can be either uniform over the whole object,
or non-uniform. In case of stretching modulus, we can have spring
stiffness the same for all edges in the whole object, or we can choose
the value for every edge of the object separately. Analogically, for
``kslin``, for ``kal`` and ``kb``. Therefore, there are two options for
setting ``ks``, ``kslin``, ``kal`` and ``kb`` stiffness. Here is the
explanation for ``ks``:

-  **Uniform stiffness:** To set uniform hyperelastic stiffness for all
   edges in the object, use ``ks``

-  **Non-uniform stiffness:** To set non-uniform hyperelastic stiffness,
   prepare a file with number of lines equal to the number of edges of
   the triangulation. Each line should contain a real number between 0
   and 1, so called “weight”. Then call ``ks`` This command reads the
   weights :math:`weight_i` for each edge and the stiffness for that
   edge is set to

   .. math:: ks_i = ksMin * (1 - weight_i) + ksMax*(weight_i)

   For bending stiffness, must contain the same number of lines as there
   are edges in the object. However, for local area preservation, the
   stiffness constant is linked to triangles. Therefore, must contain
   the same number of lines as there are triangles in the object.

.. _Elastic objects:

Elastic objects
^^^^^^^^^^^^^^^

object-id template-id origin part-type

Using a previously defined template , this command creates a new object.
Features ``OIF_LOCAL_FORCES``, ``OIF_GLOBAL_FORCES``,
``OIF_MEMBRANE_COLLISION`` are needed, if the template used the
corresponding elastic moduli.

unique ID for each object, the first object has the ID 0. The following
ones should be numbered consecutively.

object will be created using nodes, triangle incidences, elasticity
parameters and initial stretching saved in this template.

center of the object will be at this point.

can be any integer starting at 0. All particles of one object have the
same ``part-type``. One can have more objects with the same type of
particles, but this is not recommended, because the interactions between
objects are set up using these types.

angles in radians, by which the object is initially rotated around the
:math:`x, y, z` axis. Default values are 0.0 0.0 0.0.

this parameter refers to the mass of one particle (one mesh point of the
triangulation). For the proper setting, the mass of the whole membrane
must be distributed to all mesh points. Default value is 1.0.

.. _Mesh analysis:

Mesh analysis
^^^^^^^^^^^^^

oif\_mesh\_analyze nodes-file triangles-file

This command is useful for some preparatory work with mesh before it is
used for creating elastic objects.

- file with coordinates of the mesh nodes. The center of the object
  should be as close to (0,0,0) as possible.

- file with incidences for all triangles. Each line of this file
  contains three integer IDs (starting from 0) with indices of three
  vertices forming one triangle.

checks whether all triangles of the surface mesh are properly oriented.
For now, only works for convex (or almost convex) objects.

outputs the corrected file into . For now, only works for convex (or
almost convex) objects. needs to be set to 1.

subtracts 1 from all numbers in and saves a new file . This is useful,
if the mesh generating software starts numbering the particles from 1
instead of 0.

.. _Output information about specific object:

Output information about specific object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

object-id

This command is used to output information about the object that can be
used for visualisation or as input for other simulations.

- the id of the object

outputs the mesh of the object to the desired . Paraview can directly
visualize this file.

the same as the previous option, however the whole object is shift such
that it is visualized within the simulation box. This option is useful
for simulating periodical processes when objects flowing out on one side
of simulation box are transferred to the opposite side.

outputs affinity bonds that are currently activated. If no bonds are
present, the file will be generated anyway with no bonds to visualize.
Paraview can directly visualize this file.

outputs the positions of the mesh nodes to . In fact, this command
creates a new file that can be used by ``oif_object_set``. The center of
the object is located at point (0,0,0). This command is aimed to store
the deformed shape in order to be loaded later.

.. _Descriptive information about specific object:

Descriptive information about specific object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

object-id

This command is used to output information about the properties of the
object. Some of these properties can also be visualized.

- the id of the object

- outputs the location of the center of the object

computes six extremal coordinates of the object. More precisely, runs
through the all mesh points and remembers the minimal and maximal
:math:`x`-coordinate, :math:`y`-coordinate and :math:`z`-coordinate. If
is one of these: *z-min, z-max, x-min, x-max, y-min, y-max* then the
procedure returns one number according to the value of . If is , then
the procedure returns a list of six numbers, namely *x-min, x-max,
y-min, y-max, z-min, z-max*.

- outputs the approximate location of the center of the object. It is
  computed as average of 6 mesh points that have extremal :math:`x`,
  :math:`y` and :math:`z` coordinates at the time of object loading.

- outputs the minimum, average and maximum edge length of the object and
  corresponding standard deviation

- outputs the current volume of the object

- outputs the current surface of the object

- outputs the current average velocity of the object. Runs over all mesh
  points and calculates their average velocity.

.. _Setting properties for specific object:

Setting properties for specific object
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

object-id

This command sets some properties of the object.

- the id of the object

- sets the force vector () to all mesh nodes of the object. Setting is
  done using command ``part $i set ext_force`` . Note, that this command
  sets the external force in each integrate step. So if you want to use
  the external force only in one iteration, you need to set zero external
  force in the following integrate step

- moves the object so that the origin has coordinates

- deforms the object such that its origin stays unchanged, however the
  relative positions of the mesh points are taken from file . The file
  should contain the coordinates of the mesh points with the origin’s
  location at (0,0,0). The procedure also checks whether number of lines
  in the file is the same as the number of triangulation nodes of the
  object.

- stops all the particles in the object (analogue to the ``part ``
  ``fix 1 1 1`` command for single particles).

- releases the particles in the object (analogue to the ``part ``
  `` unfix`` command for single particles).

.. |image| image:: figures/oif.png


.. _Object-in-fluid interactions:

Object-in-fluid interactions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Please cite :cite:`cimrak` when using the interactions in this section in order to
simulate extended objects embedded in a LB fluid. For more details also
see the documentation at http://cell-in-fluid.fri.uniza.sk/oif-documentation.

The following interactions are implemented in order to mimic the
mechanics of elastic or rigid objects immersed in the LB fluid flow.
Their mathematical formulations were inspired by
:cite:`dupin07`. Details on how the bonds can be used for
modeling objects are described in section :ref:`Object-in-fluid`.

.. _OIF local forces:

OIF local forces
^^^^^^^^^^^^^^^^

OIF local forces are available through the :class:`espressomd.interactions.OifLocalForces` class.

This type of interaction is available for closed 3D immersed objects as
well as for 2D sheet flowing in the 3D flow.

This interaction comprises three different concepts. The local
elasticity of biological membranes can be captured by three different
elastic moduli. Stretching of the membrane, bending of the membrane and
local preservation of the surface area. Parameters
:math:`{L^0_{AB}},\ {k_s},\ {k_{slin}}` define the stretching,
parameters :math:`\phi,\ k_b` define the bending, and
:math:`A_1,\ A_2,\ k_{al}` define the preservation of local area. They
can be used all together, or, by setting any of
:math:`k_s, k_{slin}, k_b, k_{al}` to zero, the corresponding modulus
can be turned off.

.. _Stretching:

Stretching
""""""""""

For each edge of the mesh, :math:`L_{AB}` is the current distance between point :math:`A` and
point :math:`B`. :math:`L^0_{AB}` is the distance between these points in the relaxed state, that
is if the current edge has the length exactly , then no forces are
added. :math:`\Delta L_{AB}` is the deviation from the relaxed
state, that is :math:`\Delta L_{AB} = L_{AB} - L_{AB}^0`. The
stretching force between :math:`A` and :math:`B` is calculated using

.. math:: F_s(A,B) = (k_s\kappa(\lambda_{AB}) + k_{s,\mathrm{lin}})\Delta L_{AB}n_{AB}.

Here, :math:`n_{AB}` is the unit vector pointing from :math:`A` to :math:`B`, `k_s` is the
constant for nonlinear stretching, :math:`k_{s,\mathrm{lin}}` is the constant for 
linear stretching, :math:`\lambda_{AB} = L_{AB}/L_{AB}^0`, and :math:`\kappa`
is a nonlinear function that resembles neo-Hookean behavior

.. math::

   \kappa(\lambda_{AB}) = \frac{\lambda_{AB}^{0.5} + \lambda_{AB}^{-2.5}}
   {\lambda_{AB} + \lambda_{AB}^{-3}}.

Typically, one wants either nonlinear or linear behavior and therefore
one of :math:`k_s, k_{s,\mathrm{lin}}` is zero. Nonetheless the interaction will work if
both constants are non-zero.

|image_oif_streching|

.. _Bending:

Bending
"""""""

The tendency of an elastic object to maintain the resting shape is
achieved by prescribing the preferred angles between neighboring
triangles of the mesh.

Denote the angle between two triangles in the resting shape by
:math:`\theta^0`. For closed immersed objects, one always has to set the
inner angle. The deviation of this angle
:math:`\Delta \theta = \theta - \theta^0` defines two bending forces for
two triangles :math:`A_1BC` and :math:`A_2BC`

.. math:: F_{bi}(A_iBC) = k_b\frac{\Delta \theta}{\theta^0} n_{A_iBC}

Here, :math:`n_{A_iBC}` is the unit normal vector to the triangle :math:`A_iBC`.
The force :math:`F_{bi}(A_iBC)` is assigned
to the vertex not belonging to the common edge. The opposite force
divided by two is assigned to the two vertices lying on the common edge.
This procedure is done twice, for :math:`i=1` and for
:math:`i=2`.

|image_oif_bending|

.. _Local area conservation:

Local area conservation
"""""""""""""""""""""""

This interaction conserves the area of the triangles in the
triangulation.

The deviation of the triangle surface :math:`S_{ABC}` is computed from the triangle
surface in the resting shape
:math:`\Delta S_{ABC} = S_{ABC} - S_{ABC}^0`. The area
constraint assigns the following shrinking/expanding force to every
vertex

.. math:: F_{al}(A) = -k_{al}\frac{\Delta S_{ABC}}{\sqrt{S_{ABC}}}w_{A}

where :math:`k_{al}` is the area constraint coefficient, and :math:`w_{A}` is the unit vector
pointing from the centroid of triangle :math:`ABC` to the vertex :math:`A`. Similarly the
analogical forces are assigned to :math:`B` and :math:`C`.

.. todo:: Rest of this section is still Tcl syntax

OIF local force is asymmetric. After creating the interaction

::

    inter 33 oif_local_force 1.0 0.5 0.0 1.7 0.6 0.2 0.3 1.1

it is important how the bond is created. Particles need to be mentioned
in the correct order. Command

::

    part 0 bond 33 1 2 3

creates a bond related to the triangles 012 and 123. The particle 0
corresponds to point A1, particle 1 to C, particle 2 to B and particle 3
to A2. There are two rules that need to be fulfilled:

-  there has to be an edge between particles 1 and 2

-  orientation of the triangle 012, that is the normal vector defined as
   a vector product :math:`01 \times 02`, must point to the inside of
   the immersed object.

Then the stretching force is applied to particles 1 and 2, with the
relaxed length being 1.0. The bending force is applied to preserve the
angle between triangles 012 and 123 with relaxed angle 1.7 and finally,
local area force is applied to both triangles 012 and 123 with relaxed
area of triangle 012 being 0.2 and relaxed area of triangle 123 being
0.3.

Notice that also concave objects can be defined. If :math:`\theta_0` is
larger than :math:`\pi`, then the inner angle is concave.

.. _OIF global forces:

OIF global forces
^^^^^^^^^^^^^^^^^

OIF global forces are available through the
:class:`espressomd.interactions.OifGlobalForces` class.

This type of interaction is available solely for closed 3D immersed
objects.

It comprises two concepts: preservation of global surface
and of volume of the object. The parameters :math:`S^0, k_{ag}`
define preservation of the surface while parameters
:math:`V^0, k_{v}` define volume preservation. They can be
used together, or, by setting either :math:`k_{ag}` or :math:`k_{v}` to
zero, the corresponding modulus can be turned off.

.. _Global area conservation:

Global area conservation
""""""""""""""""""""""""

The global area conservation force is defined as

.. math:: F_{ag}(A) = - k_{ag}\frac{\Delta S}{S}w_{A},

where :math:`S` denotes the current surface of the immersed object, :math:`S_0` the surface in
the relaxed state and :math:`\Delta S = S - S_0`.

Here, the above mentioned force divided by 3 is added to all three
particles.

|image_oif_area|

.. _Volume conservation:

Volume conservation
"""""""""""""""""""

The deviation of the objects volume :math:`V` is computed from the volume in the
resting shape :math:`\Delta V = V - V^0`. For each
triangle the following force is computed

.. math:: F_v(ABC) = -k_v\frac{\Delta V}{V^0} S_{ABC} n_{ABC}

where :math:`S_{ABC}` is the area of triangle :math:`ABC`, :math:`n_{ABC}` is the
normal unit vector of the plane spanned by :math:`ABC`, and :math:`k_v`
is the volume constraint coefficient. The volume of one immersed object
is computed from

.. math:: V = \sum_{ABC}S_{ABC}\ n_{ABC}\cdot h_{ABC},

where the sum is computed over all triangles of the mesh and :math:`h_{ABC}` is the
normal vector from the centroid of triangle :math:`ABC` to any plane which does not
cross the cell. The force :math:`F_v(ABC)` is equally distributed to all three vertices
:math:`A, B, C.`

|image_oif_volume|

.. todo:: Rest of section still Tcl syntax

This interaction is symmetric. After the definition of the interaction
by

::

    inter 22 oif_global_force 65.3 3.0 57.0 2.0

the order of vertices is crucial. By the following command the bonds are
defined

::

    part 0 bond 22 1 2

Triangle 012 must have correct orientation, that is the normal vector
defined by a vector product :math:`01\times02`. The orientation must
point inside the immersed object.

.. _Out direction:

Out direction
^^^^^^^^^^^^^

inter oif_out_direction

This type of interaction is primarily for closed 3D immersed objects to
compute the input for membrane collision. After creating the interaction

::

    inter 66 oif_out_direction

it is important how the bond is created. Particles need to be mentioned
in the correct order. Command

::

    part 0 bond 66 1 2 3

calculates the outward normal vector of triangle defined by particles 1,
2, 3 (these should be selected in such a way that particle 0 lies
approximately at its centroid - for OIF objects, this is automatically
handled by oif_create_template command, see Section
[ssec:oif-create-template]). In order for the direction to be outward
with respect to the underlying object, the triangle 123 needs to be
properly oriented (as explained in the section on volume in
oif_global_forces interaction).

.. |image_oif_streching| image:: figures/stretching.png
.. |image_oif_bending| image:: figures/bending.png
.. |image_oif_area| image:: figures/arealocal.png
.. |image_oif_volume| image:: figures/volume.png
