.. _Input and Output:

Input and Output
================

.. _No generic checkpointing:

Checkpointing and restoring a simulation
----------------------------------------

One of the most asked-for feature that seems to be missing is
*checkpointing*, a simple way to store and restore the current
state of the simulation, and to be able to write this state to or read
it from a file. This would be most useful to be able to restart a
simulation from a specific point in time.

Unfortunately, it is impossible to provide a simple command
(``checkpoint``), out of two reasons. The main reason is that it has no
way to determine what information constitutes the actual state of the
simulation. Scripts sometimes use variables that
contain essential information about a simulation: the stored values of
an observable that was computed in previous time steps, counters, etc.
These would have to be contained in a checkpoint. However, not all
variables are of interest.

Another problem with a generic checkpoint would be the control flow of
the script. In principle, the checkpoint would have to store where in
the script the checkpointing function was called to be able to return
there. All this is even further complicated by the fact that |es| is
running in parallel.

Having said that, |es| does provide functionality which aims to store the state of the simulation engine.
In addition, variables declared in the simulation script can be added to the checkpoint.
The checkpoint data can then later be restored by calling one
load function that will automatically process the checkpoint data by
setting the user variables and restore the components of the simulation.
Furthermore, the checkpointing can be triggered by system signals that
are invoked for example when the simulation is aborted by the user or by
a timeout.

The checkpointing functionality is difficult to test for all possible simulation setups. Therefore, it is to be used with care.
It is strongly recommended to keep track of the times in the simulation run where a checkpoint was written and restored and manually verify that the observables of interest do not jump or drift after restoring the checkpoint.
Moreover, please carefully read the limitations mentioned below.

Checkpointing is implemented by the :class:`espressomd.checkpointing.Checkpoint` class. It is instanced as follows::

    from espressomd import checkpointing
    checkpoint = checkpointing.Checkpoint(checkpoint_id="mycheckpoint", checkpoint_path=".")

Here, ``checkpoint_id`` denotes the identifier for a checkpoint. Legal characters for an id
are "0-9", "a-zA-Z", "-", "_".
The parameter ``checkpoint_path``, specifies the relative or absolute path where the checkpoints are
stored. The current working directory is assumed, when this parameter is skipped.

After the simulation system and user variables are set up, they can be
registered for checkpointing.
Name the string of the object or user variable that should be registered for
checkpointing.

To give an example::

    my_var = "some variable value"
    system = espressomd.System(box_l=[100.0, 100.0, 100.0])
    # ... set system properties like time_step here ...
    checkpoint.register("system")
    checkpoint.register("my_var")
    # ...

will register the user variable ``my_var`` and the instance of the simulation system. The checkpoint can be saved via::


    checkpoint.save()

To trigger the checkpoint when Ctrl+C is pressed during a running simulation, the corresponding signal has to be registered::


    import signal
    # signal.SIGINT: signal 2, is sent when ctrl+c is pressed
    checkpoint.register_signal(signal.SIGINT)

In the above example checkpointing is triggered, when the user interrupts by
pressing Ctrl+C. In this case a new checkpoint is written and the simulation
quits.

An existing checkpoint can be loaded with::

    import espressomd
    from espressomd import checkpointing
    import signal

    checkpoint = checkpointing.Checkpoint(checkpoint_id="mycheckpoint")
    checkpoint.load()

This will restore the state of the objects registered for checkpointing.
The checkpointing instance itself will also be restored. I.e., the same variables will be registered for the next checkpoint and the same system signals will be caught as in the initial setup of the checkpointing.

Be aware of the following limitations:

  * Checkpointing makes use of the ``pickle`` python package. Objects will only be restored as far as they support pickling. This is the case for Python's basic data types, ``numpy`` arrays and many other objects. Still, pickling support cannot be taken for granted.

  * Pickling support of the :class:`espressomd.system.System` instance and contained objects such as bonded and non-bonded interactions and electrostatics methods. However, there are many more combinations of active interactions and algorithms than can be tested.

  * The active actors, i.e., the content of ``system.actors``, are checkpointed. For lattice-Boltzmann fluids, this only includes the parameters such as the lattice constant (``agrid``). The actual flow field has to be saved separately with the lattice-Boltzmann specific methods
    :meth:`espressomd.lb.HydrodynamicInteraction.save_checkpoint`
    and loaded via :meth:`espressomd.lb.HydrodynamicInteraction.load_checkpoint` after restoring the checkpoint.

  * References between Python objects are not maintained during checkpointing. For example, if an instance of a shape and an instance of a constraint containing the shape are checkpointed, these two objects are equal before checkpointing but independent copies which have the same parameters after restoring the checkpoint. Changing one will no longer affect the other.

  * The state of the cell system as well as the MPI node grid are checkpointed. Therefore, checkpoints can only be loaded, when the script runs on the same number of MPI ranks.

  * Checkpoints are not compatible between different |es| versions.

  * Checkpoints may depend on the presence of other Python modules at specific versions. It may therefore not be possible to load a checkpoint in a different environment than where it was loaded.

For additional methods of the checkpointing class, see :class:`espressomd.checkpointing.Checkpoint`.

.. _Writing H5MD-files:

Writing H5MD-files
------------------

.. note::

    Requires ``H5MD`` external feature, enabled with ``-DWITH_HDF5=ON``. Also
    requires a parallel version of HDF5. On Ubuntu, this can be installed via
    either ``libhdf5-openmpi-dev`` for OpenMPI or ``libhdf5-mpich-dev`` for
    MPICH, but not ``libhdf5-dev`` which is the serial version.

For large amounts of data it's a good idea to store it in the hdf5 (H5MD
is based on hdf5) file format (see https://www.hdfgroup.org/ for
details). Currently |es| supports some basic functions for writing simulation
data to H5MD files. The implementation is MPI-parallelized and is capable
of dealing with varying numbers of particles.

To write data in a hdf5-file according to the H5MD proposal (https://nongnu.org/h5md/), first an object of the class
:class:`espressomd.io.writer.h5md.H5md` has to be created and linked to the
respective hdf5-file. This may, for example, look like:

.. code:: python

    from espressomd.io.writer import h5md
    system = espressomd.System(box_l=[100.0, 100.0, 100.0])
    # ... add particles here
    h5 = h5md.H5md(filename="trajectory.h5", write_pos=True, write_vel=True)

If a file with the given filename exists and has a valid H5MD structures,
it will be backed up to a file with suffix ".bak". This backup file will be
deleted when the new file is closed at the end of the simulation with
``h5.close()``.

The current implementation allows to write the following properties: positions,
velocities, forces, species (|es| types), and masses of the particles. In order
to write any property, you have to set the respective boolean flag as an option
to the :class:`~espressomd.io.writer.h5md.H5md` class. Currently available:

    - ``write_pos``: particle positions

    - ``write_vel``: particle velocities

    - ``write_force``: particle forces

    - ``write_species``: particle types

    - ``write_mass``: particle masses

    - ``write_ordered``: if particles should be written ordered according to their
      id (implies serial write).

In simulations with varying numbers of particles (MC or reactions), the
size of the dataset will be adapted if the maximum number of particles
increases but will not be decreased. Instead a negative fill value will
be written to the trajectory for the id. If you have a parallel
simulation, please keep in mind that the sequence of particles in general
changes from timestep to timestep. Therefore you have to always use the
dataset for the ids to track which position/velocity/force/type/mass
entry belongs to which particle. To write data to the hdf5 file, simply
call the H5md object :meth:`~espressomd.io.writer.h5md.H5md.write` method without any arguments.

.. code:: python

    h5.write()


After the last write call, you have to call the
:meth:`~espressomd.io.writer.h5md.H5md.close` method to remove
the backup file, close the datasets, etc.

H5MD files can be read and modified with the python module h5py (for
documentation see `h5py <https://docs.h5py.org/en/stable/>`_). For example,
all positions stored in the file called "h5mdfile.h5" can be read using:

.. code:: python

    import h5py
    h5file = h5py.File("h5mdfile.h5", 'r')
    positions = h5file['particles/atoms/position/value']

Furthermore, the files can be inspected with the GUI tool hdfview or visually with the
H5MD VMD plugin (see `H5MD plugin <https://github.com/h5md/VMD-h5mdplugin>`_).

For other examples, see :file:`/samples/h5md.py`


.. _Writing MPI-IO binary files:

Writing MPI-IO binary files
---------------------------

This method outputs binary data in parallel and is, thus, also suitable for
large-scale simulations. Generally, H5MD is the preferred method because the
data is easily accessible. In contrast to H5MD, the MPI-IO functionality
outputs data in a *machine-dependent format*, but has write and read
capabilities. The usage is quite simple:

.. code:: python

    from espressomd.io.mppiio import mpiio
    system = espressomd.System()
    # ... add particles here
    mpiio.write("/tmp/mydata", positions=True, velocities=True, types=True, bonds=True)

Here, :file:`/tmp/mydata` is the prefix used for several files. The call will output
particle positions, velocities, types and their bonds to the following files in
folder :file:`/tmp`:

    - :file:`mydata.head`
    - :file:`mydata.id`
    - :file:`mydata.pos`
    - :file:`mydata.pref`
    - :file:`mydata.type`
    - :file:`mydata.vel`
    - :file:`mydata.boff`
    - :file:`mydata.bond`

Depending on the chosen output, not all of these files might be created.
To read these in again, simply call :meth:`espressomd.io.mpiio.Mpiio.read`. It has the same signature as
:meth:`espressomd.io.mpiio.Mpiio.write`.

*WARNING*: Do not attempt to read these binary files on a machine with a different
architecture!

.. _Writing VTF files:

Writing VTF files
-----------------

The formats VTF (**V**\ TF **T**\ rajectory **F**\ ormat), VSF
(**V**\ TF **S**\ tructure **F**\ ormat) and VCF (**V**\ TF
**C**\ oordinate **F**\ ormat) are formats for the visualization
software VMD: :cite:`humphrey96a`. They are intended to
be human-readable and easy to produce automatically and modify.

The format distinguishes between *structure blocks* that contain the
topological information of the system (the system size, particle names,
types, radii and bonding information, amongst others), while *coordinate
blocks* (a.k.a. as *timestep blocks*) contain the coordinates for the
particles at a single timestep. For a visualization with VMD, one
structure block and at least one coordinate block is required.

Files in the VSF format contain a single structure block, files in the
VCF format contain at least one coordinate block, while files in the VTF
format contain a single structure block (usually as a header) and an arbitrary number of
coordinate blocks (time frames) afterwards, thus allowing to store all information for
a whole simulation in a single file. For more details on the format,
refer to the VTF homepage (https://github.com/olenz/vtfplugin/wiki).

Creating files in these formats from within is supported by the commands :meth:`espressomd.io.writer.vtf.writevsf`
and :meth:`espressomd.io.writer.vtf.writevcf`, that write a structure and coordinate block (respectively) to the
given file. To create a standalone VTF file, first use ``writevsf`` at the beginning of
the simulation to write the particle definitions as a header, and then ``writevcf``
to generate a timeframe of the simulation state. For example:

A standalone VTF file can simply be

.. code:: python

    import espressomd
    from espressomd.io.writer import vtf
    system = espressomd.System(box_l=[100.0, 100.0, 100.0])
    fp = open('trajectory.vtf', mode='w+t')

    # ... add particles here

    # write structure block as header
    vtf.writevsf(system, fp)
    # write initial positions as coordinate block
    vtf.writevcf(system, fp)

    # integrate and write the frame
    for n in num_steps:
        system.integrator.run(100)
        vtf.writevcf(system, fp)
    fp.close()

The structure definitions in the VTF/VSF formats are incremental, the user
can easily add further structure lines to the VTF/VSF file after a
structure block has been written to specify further particle properties
for visualization.

Note that the ``ids`` of the particles in |es| and VMD may differ. VMD requires
the particle ids to be enumerated continuously without any holes, while
this is not required in |es|. When using ``writevsf``
and ``writevcf``, the particle ids are
automatically translated into VMD particle ids. The function allows the
user to get the VMD particle id for a given |es| particle id.

One can specify the coordinates of which particles should be written using ``types``.
If ``types='all'`` is used, all coordinates will be written (in the ordered timestep format).
Otherwise, has to be a list specifying the pids of the particles.

Also note, that these formats can not be used to write trajectories
where the number of particles or their types varies between the
timesteps. This is a restriction of VMD itself, not of the format.

.. _writevsf\: Writing the topology:

``writevsf``: Writing the topology
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:meth:`espressomd.io.writer.vtf.writevsf`

Writes a structure block describing the system's structure to the given channel, for example:

.. code:: python

    import espressomd
    from espressomd.io.writer import vtf
    system = espressomd.System(box_l=[100.0, 100.0, 100.0])
    # ... add particles here
    fp = open('trajectory.vsf', mode='w+t')
    vtf.writevsf(system, fp, types='all')

The output of this command can be
used for a standalone VSF file, or at the beginning of a VTF file that
contains a trajectory of a whole simulation.

.. _writevcf\: Writing the coordinates:

``writevcf``: Writing the coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:meth:`espressomd.io.writer.vtf.writevcf`

Writes a coordinate (or timestep) block that contains all coordinates of
the system's particles.

.. code:: python

    import espressomd
    from espressomd.io.writer import vtf
    system = espressomd.System(box_l=[100.0, 100.0, 100.0])
    # ... add particles here
    fp = open('trajectory.vcf', mode='w+t')
    vtf.writevcf(system, fp, types='all')

.. _vtf_pid_map\: Going back and forth between |es| and VTF indexing:

:meth:`espressomd.io.writer.vtf.vtf_pid_map`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Generates a dictionary which maps |es| particle ``id`` to VTF indices.
This is motivated by the fact that the list of |es| particle ``id`` is allowed to contain *holes* but VMD
requires increasing and continuous indexing. The |es| ``id`` can be used as *key* to obtain the VTF index as the *value*, for example:

.. code:: python

    import espressomd
    from espressomd.io.writer import vtf
    system = espressomd.System(box_l=[100.0, 100.0, 100.0])
    system.part.add(id=5, pos=[0, 0, 0])
    system.part.add(id=3, pos=[0, 0, 0])
    vtf_index = vtf.vtf_pid_map(system)
    vtf_index[3]

Note that the |es| particles are ordered in increasing order, thus ``id=3`` corresponds to the zeroth VTF index.

.. _Writing various formats using MDAnalysis:

Writing various formats using MDAnalysis
----------------------------------------

If the MDAnalysis package (https://mdanalysis.org) is installed, it
is possible to use it to convert frames to any of the supported
configuration/trajectory formats, including PDB, GROMACS, GROMOS,
CHARMM/NAMD, AMBER, LAMMPS, ...

To use MDAnalysis to write in any of these formats, one has first to prepare a stream from
the |es| particle data using the class :class:`espressomd.MDA_ESP`, and then read from it
using MDAnalysis. A simple example is the following:

.. code:: python

    import espressomd
    import MDAnalysis as mda
    from espressomd import MDA_ESP
    system = espressomd.System(box_l=[100.0, 100.0, 100.0])
    # ... add particles here
    eos = MDA_ESP.Stream(system)  # create the stream
    u = mda.Universe(eos.topology, eos.trajectory)  # create the MDA universe

    # example: write a single frame to PDB
    u.atoms.write("system.pdb")

    # example: save the trajectory to GROMACS format
    from MDAnalysis.coordinates.TRR import TRRWriter
    W = TRRWriter("traj.trr", n_atoms=len(system.part))  # open the trajectory file
    for i in range(100):
        system.integrator.run(1)
        u.load_new(eos.trajectory)  # load the frame to the MDA universe
        W.write_next_timestep(u.trajectory.ts)  # append it to the trajectory

For other examples, see :file:`/samples/MDAnalysisIntegration.py`

.. _Reading various formats using MDAnalysis:

Reading various formats using MDAnalysis
----------------------------------------

MDAnalysis can read various formats, including MD topologies and trajectories.
To read a PDB file containing a single frame::

    import MDAnalysis
    import numpy as np
    import espressomd
    from espressomd.interactions import HarmonicBond

    # parse protein structure
    universe = MDAnalysis.Universe("protein.pdb")
    # extract only the C-alpha atoms of chain A
    chainA = universe.select_atoms("name CA and segid A")
    # use the unit cell as box
    box_l = np.ceil(universe.dimensions[0:3])
    # setup system
    system = espressomd.System(box_l=box_l)
    system.time_step = 0.001
    system.cell_system.skin = 0.4
    # configure sphere size sigma and create a harmonic bond
    system.non_bonded_inter[0, 0].lennard_jones.set_params(
        epsilon=1, sigma=1.5, cutoff=2, shift="auto")
    system.bonded_inter[0] = HarmonicBond(k=0.5, r_0=1.5)
    # create particles and add bonds between them
    system.part.add(pos=np.array(chainA.positions, dtype=float))
    for i in range(0, len(chainA) - 1):
        system.part[i].add_bond((system.bonded_inter[0], system.part[i + 1].id))
    # visualize protein in 3D
    from espressomd import visualization
    visualizer = visualization.openGLLive(system, bond_type_radius=[0.2])
    visualizer.run(0)
