Features
========

This chapter describes the features that can be activated in . Even if
possible, it is not recommended to activate all features, because this
will negatively effect ’s performance.

Features can be activated in the configuration header ``myconfig.hpp``
(see section ). Too activate ``FEATURE``, add the following line to the
header file:

::

    #define FEATURE

General features
----------------

-  By default, all coordinates in are periodic. With
   ``PARTIAL_PERIODIC`` turned on, the global variable ``periodic`` (see
   section ) controls the periodicity of the individual coordinates.
   Note that this slows the integrator down by around :math:`10-30\%`.

-  This switches on the various electrostatics algorithms, such as P3M.
   See section for details on these algorithms.

-  This activates the dipole-moment property of particles; In addition,
   the various magnetostatics algorithms, such as P3M are switched on.
   See section for details on these algorithms.

-  Switch on rotational degrees of freedom for the particles, as well as
   the corresponding quaternion integrator. See section for details.
   Note, that when the feature is activated, every particle has three
   additional degrees of freedom, which for example means that the
   kinetic energy changes at constant temperature is twice as large.

-  Allows to set whether a particle has rotational degrees of freedom.

-  Allows to choose the Langevin temperature and friction coefficient
   per particle.

-  
-  Allows to define an arbitrary constant force for each particle
   individually. Also allows to fix individual coordinates of particles,
   keep them at a fixed position or within a plane.

-  Turns on various spatial constraints such as spherical compartments
   or walls. This constraints interact with the particles through
   regular short ranged potentials such as the Lennard–Jones potential.
   See section for possible constraint forms.

-  Switch on tunable slip conditions for planar wall boundary
   conditions. See section for details.

-  Allows particles to have individual masses. Note that some analysis
   procedures have not yet been adapted to take the masses into account
   correctly.

-  Allows to exclude specific short ranged interactions within
   molecules.

-  Allows to pull apart groups of particles

-  Allows to fix the center of mass of all particles of a certain type.

-  
-  Turns on the RATTLE integrator which allows for fixed lengths bonds
   between particles.

-  Virtual sites are particles, the position and velocity of which is
   not obtained by integrating equations of motion. Rather, they are
   placed using the position (and orientation) of other particles. The
   feature allows to place a virtual particle into the center of mass of
   a set of other particles. See section [sec:virtual] for details.

-  Virtual sites are particles, the position and velocity of which is
   not obtained by integrating equations of motion. Rather, they are
   placed using the position (and orientation) of other particles. The
   feature allows for rigid arrangements of particles. See section
   [sec:virtual] for details.

-  
-  
-  
-  
-  
-  
-  
-  Allows to define the temperature and friction coefficient for
   individual particles. See [ssec:particleproperties] for details.

-  Allows the user to define three particle types to be reactant,
   catalyzer, and product. Reactants get converted into products in the
   vicinity of a catalyst according to a used-defined reaction rate
   constant. It is also possible to set up a chemical equilibrium
   reaction between the reactants and products, with another rate
   constant. See section [sec:Reactions] for details.

-  
-  Allows particles to be bound on collision. See section
   [sec:collision]

-  This switches back to the old, *wrong* random walk code of the
   polymer. Only use this if you rely on the old behaviour and *know
   what you are doing*.

-  Allows to write data to H5MD formatted hdf5 files.

In addition, there are switches that enable additional features in the
integrator or thermostat:

-  Enables the non-equilbrium (shear) MD support (see section ).

-  Enables an on–the–fly NPT integration scheme (see section ).

-  Enables the dissipative particle dynamics thermostat (see section ).

-  Enables the transversal dissipative particle dynamics thermostat (see
   section ).

-  Enables the dissipative particle dynamics thermostat implemented as
   an interaction, allowing to choose different parameters between
   different particle types (see section ).

-  
-  Enables masses in DPD using reduced, dimensionless mass units.

-  Enables masses in DPD using absolute mass units.

-  Enables the lattice-Boltzmann fluid code (see section ).

-  Enables the lattice-Boltzmann fluid code support for GPU (see section
   ).

-  Enables the Shan Chen bicomponent fluid code on the GPU (see section
   ).

-  Enables the implicit calculation of electro-hydrodynamics for charged
   particles and salt ions in an electric field.

Interactions
------------

The following switches turn on various short ranged interactions (see
section ):

-  Enable support for user–defined interactions.

-  Enable the Lennard–Jones potential.

-  Enable the generic Lennard–Jones potential with configurable
   exponents and individual prefactors for the two terms.

-  Enable the Lennard–Jones potential with a cosine–tail.

-  Same as LJCOS, but using a slightly different way of smoothing the
   connection to 0.

-  Enable the directional Lennard–Jones potential.

-  
-  
-  
-  
-  Enable the Morse potential.

-  Enable the Buckingham potential.

-  Enable the soft sphere potential.

-  Enable the smooth step potential, a step potential with two length
   scales.

-  Enable the Born-Meyer-Huggins-Tosi-Fumi potential, which can be used
   to model salt melts.

Some of the short range interactions have additional features:

-  This adds an additional check to the Lennard–Jones potentials that
   prints a warning if particles come too close so that the simulation
   becomes unphysical.

-  Switch the interface of the dihedral potential to its old, less
   flexible form. Use this for older scripts that are not yet adapted to
   the new interface of the dihedral potential.

If you want to use bond-angle potentials (see section ), you need the
followig features.

-  
-  
-  

Debug messages
--------------

Finally, there are a number of flags for debugging. The most important
one are

-  Enables numerous additional checks which can detect inconsistencies
   especially in the cell systems. This checks are however too slow to
   be enabled in production runs.

-  Enables an internal memory allocation checking system. This produces
   output for each allocation and freeing of a memory chunk, and
   therefore allows to track down memory leaks. This works by internally
   replacing ``malloc``, ``realloc`` and ``free``.

The following flags control the debug output of various sections of
Espresso. You will however understand the output very often only by
looking directly at the code.

-  Output from the asynchronous communication code.

-  Notifications for event calls, i. e. the ``on_?`` functions in
   ``initialize.c``. Useful if some module does not correctly respond to
   changes of e. g. global variables.

-  Integrator output.

-  Cellsystem output.

-  Cellsystem output specific to the handling of ghost cells and the
   ghost cell communication.

-  
-  Debugging of the Verlet list code of the domain decomposition cell
   system.

-  Universal lattice structure debugging.

-  
-  
-  Output from the particle handling code.

-  
-  debugging of P\ :math:`^3`\ Ms real space part.

-  debugging of P\ :math:`^3`\ Ms :math:`k`–space part.

-  
-  Output from the unified FFT code.

-  
-  
-  Output from the force calculation loops.

-  Output from the pressure tensor calculation loops.

-  Output from the thermostats.

-  Output from the Lennard–Jones code.

-  Output from the Morse code.

-  
-  Define to a number of a particle to obtain output on the forces
   calculated for this particle.

-  
-  
-  
-  Output from the lattice–Boltzmann code.

-  
-  Introduce a barrier after each asynchronous command completion. Helps
   in detection of mismatching communication.

-  Causes to try to provoke a core dump when exiting unexpectedly.

-  Causes to try this even with MPI errors.
