Features
========

This chapter describes the features that can be activated in |es|. Even if
possible, it is not recommended to activate all features, because this
will negatively effect |es| ’s performance.

Features can be activated in the configuration header ``myconfig.hpp`` (see
section :ref:`myconifg.hpp\: Activating and deactivating features`). To
activate ``FEATURE``, add the following line to the header file:

::

    #define FEATURE

General features
----------------

-  ``PARTIAL_PERIODIC`` By default, all coordinates in |es| are periodic. With
   ``PARTIAL_PERIODIC`` turned on, the |es| global variable ``periodic``
   controls the periodicity of the individual coordinates.

   .. note:: This slows the integrator down by around :math:`10-30\%`.

   .. seealso:: :ref:`Setting global variables in Python`

-  ``ELECTROSTATICS`` This enables the use of the various electrostatics algorithms, such as P3M.

   .. seealso:: :ref:`Coulomb interaction`

-  ``INTER_RF``

-  ``MMM1D_GPU``

-  ``EWALD_GPU``

-  ``_P3M_GPU_FLOAT``

-  ``_P3M_GPU_DOUBLE``

-  ``DIPOLES`` This activates the dipole-moment property of particles; In addition,
   the various magnetostatics algorithms, such as P3M are switched on.

   .. seealso:: :ref:`Coulomb interaction`

-  ``SCAFACOS_DIPOLES``

-  ``ROTATION`` Switch on rotational degrees of freedom for the particles, as well as
   the corresponding quaternion integrator. 
   
   .. seealso:: :ref:`Setting up particles`

   .. note:: 
      Note, that when the feature is activated, every particle has three
      additional degrees of freedom, which for example means that the
      kinetic energy changes at constant temperature is twice as large.

-  ``LANGEVIN_PER_PARTICLE`` Allows to choose the Langevin temperature and friction coefficient
   per particle.

-  ``ROTATIONAL_INERTIA``

-  ``EXTERNAL_FORCES`` Allows to define an arbitrary constant force for each particle
   individually. Also allows to fix individual coordinates of particles,
   keep them at a fixed position or within a plane.

-  ``CONSTRAINTS`` Turns on various spatial constraints such as spherical compartments
   or walls. This constraints interact with the particles through
   regular short ranged potentials such as the Lennard–Jones potential.
   See section for possible constraint forms.

-  ``MASS`` Allows particles to have individual masses. Note that some analysis
   procedures have not yet been adapted to take the masses into account
   correctly.

   .. seealso:: :attr:`espressomd.particle_data.ParticleHandle.mass`

-  ``EXCLUSIONS`` Allows to exclude specific short ranged interactions within
   molecules.

   .. seealso:: :attr:`espressomd.particle_data.ParticleHandle.exclude`

-  ``COMFORCE`` Allows to pull apart groups of particles

-  ``COMFIXED`` Allows to fix the center of mass of all particles of a certain type.

-  ``MOLFORCES``

-  ``BOND_CONSTRAINT`` Turns on the RATTLE integrator which allows for fixed lengths bonds
   between particles.

-  ``VIRTUAL_SITES_COM`` Virtual sites are particles, the position and velocity of which is
   not obtained by integrating equations of motion. Rather, they are
   placed using the position (and orientation) of other particles. The
   feature allows to place a virtual particle into the center of mass of
   a set of other particles.
   
   .. seealso:: :ref:`Virtual sites` 

-  ``VIRTUAL_SITES_RELATIVE`` Virtual sites are particles, the position and velocity of which is
   not obtained by integrating equations of motion. Rather, they are
   placed using the position (and orientation) of other particles. The
   feature allows for rigid arrangements of particles.

   .. seealso:: :ref:`Virtual sites` 

-  ``VIRTUAL_SITES_NO_VELOCITY``

   .. seealso:: :ref:`Virtual sites` 

-  ``VIRTUAL_SITES_THERMOSTAT``

   .. seealso:: :ref:`Virtual sites` 

-  ``THERMOSTAT_IGNORE_NON_VIRTUAL``

-  ``MODES``

-  ``METADYNAMICS``

-  ``CATALYTIC_REACTIONS`` Allows the user to define three particle types to be reactant,
   catalyzer, and product. Reactants get converted into products in the
   vicinity of a catalyst according to a used-defined reaction rate
   constant. It is also possible to set up a chemical equilibrium
   reaction between the reactants and products, with another rate
   constant. 
   
   .. seealso:: :ref:`Catalytic reaction`

-  ``OVERLAPPED``

-  ``COLLISION_DETECTION`` Allows particles to be bound on collision.

-  ``H5MD`` Allows to write data to H5MD formatted hdf5 files.

   .. seealso:: :ref:`Writing H5MD-Files`

In addition, there are switches that enable additional features in the
integrator or thermostat:


-  ``NEMD`` Enables the non-equilbrium (shear) MD support.

   .. seealso:: :ref:`\`\`nemd\`\`\: Setting up non-equilibirum MD`

-  ``NPT`` Enables an on–the–fly NPT integration scheme.
   
   .. seealso:: :ref:`\`\`thermostat\`\`\: Setting up the thermostat`


-  ``MEMBRANE_COLLISION``

-  ``LEES_EDWARDS``

-  ``REACTION_ENSEMBLE``

-  ``GHMC``

-  ``MULTI_TIMESTEP``

-  ``ENGINE``

-  ``PARTICLE_ANISOTROPY``


Fluid dynamics and fluid structure interaction
----------------------------------------------

-  ``DPD`` Enables the dissipative particle dynamics thermostat and interaction.

   .. seealso:: :ref:`DPD interaction`

-  ``DPD_MASS_RED``

-  ``DPD_MASS_LIN``

-  ``LB`` Enables the lattice-Boltzmann fluid code.

   .. seealso:: :attr:`espressomd.lb`, :ref:`Lattice-Boltzmann`

-  ``LB_GPU`` Enables the lattice-Boltzmann fluid code support for GPU.

   .. seealso:: :attr:`espressomd.lb`, :ref:`Lattice-Boltzmann`

-  ``LB_BOUNDARIES``

-  ``LB_BOUNDARIES_GPU``

-  ``SHANCHEN`` Enables the Shan Chen bicomponent fluid code on the GPU.

-  ``AFFINITY``

-  ``LB_ELECTROHYDRODYNAMICS`` Enables the implicit calculation of electro-hydrodynamics for charged
   particles and salt ions in an electric field.

-  ``ELECTROKINETICS``

-  ``EK_BOUNDARIES``

-  ``EK_ELECTROSTATIC_COUPLING``

-  ``EK_DEBUG``

-  ``EK_DOUBLE_PREC``

-  ``IMMERSED_BOUNDARY`` Immersed-Boundary Bayreuth version.

-  ``OIF_LOCAL_FORCES``

-  ``OIF_GLOBAL_FORCES``


Interactions
------------

The following switches turn on various short ranged interactions (see
section :ref:`Isotropic non-bonded interactions`):

-  ``TABULATED`` Enable support for user–defined interactions.

-  ``LENNARD_JONES`` Enable the Lennard–Jones potential.

-  ``LENNARD_JONES_GENERIC`` Enable the generic Lennard–Jones potential with configurable
   exponents and individual prefactors for the two terms.

-  ``LJCOS`` Enable the Lennard–Jones potential with a cosine–tail.

-  ``LJCOS2`` Same as LJCOS, but using a slightly different way of smoothing the
   connection to 0.

-  ``LJ_ANGLE`` Enable the directional Lennard–Jones potential.

-  ``GAY_BERNE``

-  ``HERTZIAN``

-  ``NO_INTRA_NB``

-  ``MORSE`` Enable the Morse potential.

-  ``BUCKINGHAM`` Enable the Buckingham potential.

-  ``SOFT_SPHERE`` Enable the soft sphere potential.

-  ``SMOOTH_STEP`` Enable the smooth step potential, a step potential with two length
   scales.

-  ``BMHTF_NACL`` Enable the Born-Meyer-Huggins-Tosi-Fumi potential, which can be used
   to model salt melts.

Some of the short range interactions have additional features:

-  ``LJ_WARN_WHEN_CLOSE`` This adds an additional check to the Lennard–Jones potentials that
   prints a warning if particles come too close so that the simulation
   becomes unphysical.

-  ``OLD_DIHEDRAL`` Switch the interface of the dihedral potential to its old, less
   flexible form. Use this for older scripts that are not yet adapted to
   the new interface of the dihedral potential.

If you want to use bond-angle potentials (see section :ref:`Bond-angle interactions`), you need the
following features.

-  ``BOND_ANGLE``

-  ``BOND_ANGLEDIST``

-  ``BOND_ENDANGLEDIST``

-  ``BOND_ANGLEDIST_HARMONIC``

-  ``BOND_ENDANGLEDIST_HARMONIC``

-  ``LJGEN_SOFTCORE``

-  ``COS2``

-  ``GAUSSIAN``

-  ``HAT``

-  ``UMBRELLA``


DNA Model (Fyta DNA)
--------------------

-  ``CG_DNA``

-  ``TWIST_STACK``

-  ``HYDROGEN_BOND``

-  ``COULOMB_DEBYE_HUECKEL``

Dipolar interactions on the gpu:

-  ``DIPOLAR_BARNES_HUT`` Enable Barnes-Hut octree sum on GPU algorithm for a dipole-dipole interaction calculation.

-  ``DIPOLAR_DIRECT_SUM`` Enables calculation of dipole-dipole interactions via direct summation on the gpu




Miscellaneous
-------------

-  ``FLATNOISE`` Shape of the noise in ther (Langevin) thermostat.

-  ``GAUSSRANDOM`` Shape of the noise in ther (Langevin) thermostat.

-  ``GAUSSRANDOMCUT`` Shape of the noise in ther (Langevin) thermostat.

-  ``GHOSTS_HAVE_BONDS`` Ghost particles also have the bond information.


Debug messages
--------------

Finally, there are a number of flags for debugging. The most important
one are

-  ``ADDITIONAL_CHECKS`` Enables numerous additional checks which can detect inconsistencies
   especially in the cell systems. This checks are however too slow to
   be enabled in production runs.

The following flags control the debug output of various sections of
|es|. You will however understand the output very often only by
looking directly at the code.

-  ``COMM_DEBUG`` Output from the asynchronous communication code.

-  ``EVENT_DEBUG`` Notifications for event calls, i. e. the ``on_...`` functions in
   ``initialize.c``. Useful if some module does not correctly respond to
   changes of e. g. global variables.

-  ``INTEG_DEBUG`` Integrator output.

-  ``CELL_DEBUG`` Cellsystem output.

-  ``GHOST_DEBUG`` Cellsystem output specific to the handling of ghost cells and the
   ghost cell communication.

-  ``GHOST_FORCE_DEBUG``

-  ``VERLET_DEBUG`` Debugging of the Verlet list code of the domain decomposition cell
   system.

-  ``LATTICE_DEBUG`` Universal lattice structure debugging.

-  ``HALO_DEBUG``

-  ``GRID_DEBUG``

-  ``PARTICLE_DEBUG`` Output from the particle handling code.

-  ``P3M_DEBUG``

-  ``ESR_DEBUG`` debugging of P\ :math:`^3`\ Ms real space part.

-  ``ESK_DEBUG`` debugging of P\ :math:`^3`\ Ms :math:`k` –space part.

-  ``FFT_DEBUG`` Output from the unified FFT code.

-  ``MAGGS_DEBUG``

-  ``RANDOM_DEBUG``

-  ``FORCE_DEBUG`` Output from the force calculation loops.

-  ``PTENSOR_DEBUG`` Output from the pressure tensor calculation loops.

-  ``THERMO_DEBUG`` Output from the thermostats.

-  ``LJ_DEBUG`` Output from the Lennard–Jones code.

-  ``MORSE_DEBUG`` Output from the Morse code.

-  ``FENE_DEBUG``

-  ``ONEPART_DEBUG`` Define to a number of a particle to obtain output on the forces
   calculated for this particle.

-  ``STAT_DEBUG``

-  ``POLY_DEBUG``

-  ``MOLFORCES_DEBUG``

-  ``LB_DEBUG`` Output from the lattice–Boltzmann code.

-  ``VIRTUAL_SITES_DEBUG``

-  ``ASYNC_BARRIER`` Introduce a barrier after each asynchronous command completion. Helps
   in detection of mismatching communication.

-  ``FORCE_CORE`` Causes |es| to try to provoke a core dump when exiting unexpectedly.

-  ``MPI_CORE`` Causes |es| to try this even with MPI errors.

-  ``ESIF_DEBUG``

-  ``LE_DEBUG``

-  ``SD_DEBUG``

-  ``CUDA_DEBUG``

-  ``H5MD_DEBUG``

-  ``ONEPART_DEBUG_ID`` Use this define to supply a particle ID for which to output debug messages. For example: ``#define ONEPART_DEBUG_ID 13``

