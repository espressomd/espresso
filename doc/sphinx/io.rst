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

    import espressomd
    import espressomd.checkpointing
    checkpoint = espressomd.checkpointing.Checkpoint(checkpoint_id="mycheckpoint", checkpoint_path=".")

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
    import espressomd.checkpointing
    import signal

    checkpoint = espressomd.checkpointing.Checkpoint(checkpoint_id="mycheckpoint")
    checkpoint.load()

This will restore the state of the objects registered for checkpointing.
The checkpointing instance itself will also be restored. I.e., the same
variables will be registered for the next checkpoint and the same system
signals will be caught as in the initial setup of the checkpointing.

Be aware of the following limitations:

* Checkpointing makes use of the ``pickle`` python package. Objects will only
  be restored as far as they support pickling. This is the case for Python's
  basic data types, ``numpy`` arrays and many other objects. Still, pickling
  support cannot be taken for granted.

* Pickling of the :class:`espressomd.system.System` instance and
  contained objects such as bonded and non-bonded interactions and
  electrostatics methods is covered by basic tests. However, not all
  combinations of algorithms can be tested. If you encounter an issue
  for a specific combination of features, please share your findings
  with the |es| community.

* The active actors, i.e., the content of ``system.actors`` resp.
  ``system.ekcontainers``, are checkpointed. For lattice-based methods like
  lattice-Boltzmann fluids and advection-diffusion-reaction models, this only
  includes the parameters such as the lattice constant (``agrid``) and initial
  densities.
  The actual fields have to be saved separately with the lattice-specific
  methods :meth:`espressomd.lb.LBFluidWalberla.save_checkpoint
  <espressomd.detail.walberla.LatticeModel.save_checkpoint>` resp.
  :meth:`espressomd.electrokinetics.EKSpecies.save_checkpoint
  <espressomd.detail.walberla.LatticeModel.save_checkpoint>`
  and loaded via :meth:`espressomd.lb.LBFluidWalberla.load_checkpoint
  <espressomd.detail.walberla.LatticeModel.load_checkpoint>` resp.
  :meth:`espressomd.electrokinetics.EKSpecies.load_checkpoint
  <espressomd.detail.walberla.LatticeModel.load_checkpoint>`
  after restoring the checkpoint. See :ref:`LB checkpointing <Checkpointing LB>`
  resp. :ref:`EK checkpointing <Checkpointing EK>` for more details.

* References between Python objects are not maintained during checkpointing.
  For example, if an instance of a shape and an instance of a constraint
  containing the shape are checkpointed, these two objects are equal before
  checkpointing but independent copies which have the same parameters after
  restoring the checkpoint. Changing one will no longer affect the other.

* The state of the cell system as well as the MPI node grid are checkpointed.
  Therefore, checkpoints can only be loaded, when the script runs on the same
  number of MPI ranks.

* Checkpoints are not compatible between different |es| versions.

* Checkpoints may depend on the presence of other Python modules at specific
  versions. It may therefore not be possible to load a checkpoint in a
  different environment than where it was written.
  In particular, all |es| modules whose classes have been used to
  instantiate objects in the checkpoint file also need to be imported
  in the script that loads the checkpoint (because importing an |es|
  module does more than just making classes visibles: it also registers
  them as script interface classes).
  Loading a checkpoint without the proper |es| module imports will generally
  raise an exception indicating which module is missing.

* It is only possible to checkpoint objects at global scope.
  When wrapping the checkpointing logic in a function, objects local to
  that function won't be checkpointed, since their origin cannot be safely
  stored in the checkpoint file. To circumvent this limitation, make any
  local object explicitly global, so that it belongs to the global scope::

      import espressomd
      import espressomd.checkpointing

      def setup_system():
          global system  # attach 'system' to global scope for checkpointing
          checkpoint = espressomd.checkpointing.Checkpoint(checkpoint_id="mycheckpoint")
          if not checkpoint.has_checkpoints():
              system = espressomd.System(box_l=[1., 1., 1.])
              system.part.add(pos=[0.1, 0.2, 0.3])
              checkpoint.register("system")
              checkpoint.save()
          else:
              checkpoint.load()
              print(system.part.by_id(0).pos)
          return system

      system = setup_system()

* To be fully deterministic when loading from a checkpoint with an active
  thermostat, the first step of the integration should be called with the flag
  ``reuse_forces=True``, e.g. ``system.integrator.run(2, reuse_forces=True)``.
  This is because loading a checkpoint reinitializes the system and enforces
  a recalculation of the forces. However, this computes the forces from the
  velocities at the current time step and not at the previous half time step.
  Please note that long-range actors can make trajectories non-reproducible.
  For example, lattice-Boltzmann introduces errors of the order of 1e-15 with
  binary checkpoint files, or 1e-7 with ASCII checkpoint files. In addition,
  several electrostatic and magnetostatic actors automatically introduce
  a deviation of the order of 1e-7, either due to floating-point rounding
  errors (:class:`~espressomd.electrostatics.P3MGPU`), or due to re-tuning
  using the most recent system state (:class:`~espressomd.electrostatics.MMM1D`,
  :class:`~espressomd.electrostatics.MMM1DGPU`).
  When in doubt, you can easily verify the absence of a "force jump" when
  loading from a checkpoint by replacing the electrostatics actor with your
  combination of features in files :file:`samples/save_checkpoint.py` and
  :file:`samples/load_checkpoint.py` and running them sequentially.

For additional methods of the checkpointing class, see
:class:`espressomd.checkpointing.Checkpoint`.

.. _Writing H5MD-files:

Writing H5MD-files
------------------

.. note::

    Requires ``H5MD`` external feature, enabled with ``-D ESPRESSO_BUILD_WITH_HDF5=ON``.
    Also requires a parallel version of HDF5. On Ubuntu, this can be installed
    via either ``libhdf5-openmpi-dev`` for OpenMPI or ``libhdf5-mpich-dev`` for
    MPICH, but not ``libhdf5-dev`` which is the serial version.

For long simulations, it's a good idea to store data in the hdf5 file format
(see https://www.hdfgroup.org for details, H5MD is based on hdf5).
Currently |es| supports some basic functions for writing simulation
data to H5MD files. The implementation is MPI-parallelized and is capable
of dealing with a varying number of particles.

To write data in a hdf5-file according to the H5MD proposal
(https://nongnu.org/h5md), first an object of the class
:class:`espressomd.io.writer.h5md.H5md` has to be created and linked to the
respective hdf5-file. This may, for example, look like:

.. code-block:: python

    import espressomd.io.writer.h5md
    system = espressomd.System(box_l=[100.0, 100.0, 100.0])
    # ... add particles here
    h5 = espressomd.io.writer.h5md.H5md(file_path="trajectory.h5")

An optional argument to the constructor of :class:`espressomd.io.writer.h5md.H5md` is
an instance of :class:`espressomd.io.writer.h5md.UnitSystem` which encapsulates
physical units for time, mass, length and electrical charge.

If a file at the given filepath exists and has a valid H5MD structure,
it will be backed up to a file with suffix ".bak" and loaded into
a new file. Therefore H5MD can be used together with checkpointing.
The backup file will be deleted when the new file is closed at the end of the
simulation with :meth:`~espressomd.io.writer.h5md.H5md.close()`. The backup
file is not be erased if the simulation terminates unexpectedly.

To write data to the HDF5 file, simply call the method
:meth:`~espressomd.io.writer.h5md.H5md.write` without any arguments.
After the last write, call :meth:`~espressomd.io.writer.h5md.H5md.flush()`
and then :meth:`~espressomd.io.writer.h5md.H5md.close()`
to close the datasets and remove the backup file.

The current implementation writes the following properties by default: folded
positions, periodic image count, velocities, forces, species (|es| types),
charges and masses of the particles. While folded positions are written
to disk, the unfolded coordinates can be reconstructed from the image count.
The time-dependent box size and Lees-Edwards parameters are also stored.
Some of these properties can be opted out by specifying in argument
``fields`` the subset of fields to write to the trajectory file;
call method :meth:`~espressomd.io.writer.h5md.H5md.valid_fields()`
to find out which string corresponds to which field.

In simulations with a varying number of particles (Monte-Carlo reactions), the
size of the dataset will be adapted if the maximum number of particles
increases but will not be decreased. Instead a negative fill value will
be written to the trajectory for the id.

If you have a parallel
simulation, please keep in mind that the sequence of particles in general
changes from timestep to timestep. Therefore you have to always use the
dataset for the ids to track which position/velocity/force/type/mass
entry belongs to which particle.

For an example involving physical units, see :file:`/samples/h5md.py`.

.. _Reading H5MD-files:

Reading H5MD-files
------------------

H5MD files can be read and sometimes modified by many tools. If the data was
stored with `physical units <https://nongnu.org/h5md/modules/units.html>`__,
they can be accessed by reading the group attributes. Since the data is
written in parallel, the particles are unsorted; if particles were created
with increasing particle id and no particle deletion occurred during the
simulation, the coordinates can be sorted with a simply numpy operation.

To read with the python module ``h5py`` (documentation:
`HDF5 for Python <https://docs.h5py.org/en/stable>`__)::

    import h5py
    with h5py.File("sample.h5", mode='r') as h5file:
        positions = h5file['particles/atoms/position/value']
        positions.attrs['unit']
        forces = h5file['particles/atoms/force/value']
        forces_unit = forces.attrs['unit']
        sim_time = h5file['particles/atoms/id/time']
        print(f"last frame: {sim_time[-2]:.3f} {sim_time.attrs['unit'].decode('utf8')}")

To read with the python module ``pandas`` (documentation: `HDFStore: PyTables
<https://pandas.pydata.org/docs/reference/io.html#hdfstore-pytables-hdf5>`_)::

    import pandas
    with pandas.HDFStore("sample.h5", mode='r') as h5file:
        positions = h5file.root.particles.atoms.position.value
        positions.attrs['unit']
        forces = h5file.root.particles.atoms.force.value
        forces_unit = forces.attrs['unit']
        sim_time = h5file.root.particles.atoms.id.time
        print(f"last frame: {sim_time[-2]:.3f} {sim_time.attrs['unit'].decode('utf8')}")

To read from the command line with
`h5dump <https://support.hdfgroup.org/HDF5/doc/RM/Tools/h5dump.htm>`__
(Ubuntu package ``hdf5-tools``):

.. code-block:: sh

    # show metadata only
    h5dump --header sample.h5 | less
    # show metadata + data
    h5dump sample.h5 | less

H5MD files can also be inspected with the GUI tool
`HDFView <https://www.hdfgroup.org/downloads/hdfview>`__ (Ubuntu package
``hdfview``) or visually with the H5MD VMD plugin (GitHub project
`h5md/VMD-h5mdplugin <https://github.com/h5md/VMD-h5mdplugin>`__).

For an example involving ``h5py``, coordinates resorting and reconstruction
of the unfolded coordinates, see :file:`/samples/h5md_trajectory.py`.

.. _Writing MPI-IO binary files:

Writing MPI-IO binary files
---------------------------

This method outputs binary data in parallel and is, thus, also suitable for
large-scale simulations. Generally, H5MD is the preferred method because the
data is easily accessible. In contrast to H5MD, the MPI-IO functionality
outputs data in a *machine-dependent format*, but has write and read
capabilities. The usage is quite simple:

.. code-block:: python

    import espressomd
    import espressomd.io
    system = espressomd.System(box_l=[1, 1, 1])
    # ... add particles here
    mpiio = espressomd.io.mpiio.Mpiio()
    mpiio.write("/tmp/mydata", positions=True, velocities=True, types=True, bonds=True)

Here, :file:`/tmp/mydata` is the prefix used to generate several files.
The call will output particle positions, velocities, types and their bonds
to the following files in folder :file:`/tmp`:

- :file:`mydata.head`
- :file:`mydata.id`
- :file:`mydata.pos`
- :file:`mydata.pref`
- :file:`mydata.type`
- :file:`mydata.vel`
- :file:`mydata.boff`
- :file:`mydata.bond`

Depending on the chosen output, not all of these files might be created.
To read these in again, simply call :meth:`espressomd.io.mpiio.Mpiio.read`.
It has the same signature as :meth:`espressomd.io.mpiio.Mpiio.write`.
When writing files, make sure the prefix hasn't been used before
(e.g. by a different simulation script), otherwise the write operation
will fail to avoid accidentally overwriting pre-existing data. Likewise,
reading incomplete data (or complete data but with the wrong number of MPI
ranks) will throw an error.

*WARNING*: Do not attempt to read these binary files on a machine
with a different architecture! This will read malformed data without
necessarily throwing an error.

In case of read failure or write failure, the simulation will halt.
On 1 MPI rank, the simulation will halt with a python runtime error.
This exception can be recovered from; in case of a write operation,
any written file must be deleted before attempting to write again
(since the prefix argument must be unique). On more than 1 MPI rank,
the simulation will halt with a call to ``MPI_Abort`` and will send
the ``SIGABRT`` signal.

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

.. code-block:: python

    import espressomd
    import espressomd.io.writer.vtf
    system = espressomd.System(box_l=[100.0, 100.0, 100.0])
    fp = open('trajectory.vtf', mode='w+t')

    # ... add particles here

    # write structure block as header
    espressomd.io.writer.vtf.writevsf(system, fp)
    # write initial positions as coordinate block
    espressomd.io.writer.vtf.writevcf(system, fp)

    # integrate and write the frame
    for n in num_steps:
        system.integrator.run(100)
        espressomd.io.writer.vtf.writevcf(system, fp)
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

.. code-block:: python

    import espressomd
    import espressomd.io.writer.vtf
    system = espressomd.System(box_l=[100.0, 100.0, 100.0])
    # ... add particles here
    fp = open('trajectory.vsf', mode='w+t')
    espressomd.io.writer.vtf.writevsf(system, fp, types='all')

The output of this command can be
used for a standalone VSF file, or at the beginning of a VTF file that
contains a trajectory of a whole simulation.

.. _writevcf\: Writing the coordinates:

``writevcf``: Writing the coordinates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:meth:`espressomd.io.writer.vtf.writevcf`

Writes a coordinate (or timestep) block that contains all coordinates of
the system's particles.

.. code-block:: python

    import espressomd
    import espressomd.io.writer.vtf
    system = espressomd.System(box_l=[100.0, 100.0, 100.0])
    # ... add particles here
    fp = open('trajectory.vcf', mode='w+t')
    espressomd.io.writer.vtf.writevcf(system, fp, types='all')

.. _vtf_pid_map\: Going back and forth between ESPResSo and VTF indexing:

``vtf_pid_map``: Going back and forth between |es| and VTF indexing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:meth:`espressomd.io.writer.vtf.vtf_pid_map`

Generates a dictionary which maps |es| particle ``id`` to VTF indices.
This is motivated by the fact that the list of |es| particle ``id`` is allowed to contain *holes* but VMD
requires increasing and continuous indexing. The |es| ``id`` can be used as *key* to obtain the VTF index as the *value*, for example:

.. code-block:: python

    import espressomd
    import espressomd.io.writer.vtf
    system = espressomd.System(box_l=[100.0, 100.0, 100.0])
    system.part.add(id=5, pos=[0, 0, 0])
    system.part.add(id=3, pos=[0, 0, 0])
    vtf_index = espressomd.io.writer.vtf.vtf_pid_map(system)
    vtf_index[3]

Note that the |es| particles are ordered in increasing order, thus ``id=3`` corresponds to the zeroth VTF index.

.. _Reading VTK files:

Reading VTK files
-----------------

The waLBerla library writes VTK multi-piece uniform grids in XML format.
Each piece contains information about its spatial extent, from which it is
possible to deduce the grid dimensions. Each piece may contain one or more
array, which are uniquely identified by name. While the Python package ``vtk``
provides tools to read VTK files as numpy arrays, it doesn't automatically
reconstruct the 3D grids using the topology information of each piece; this
functionality is provided by the wrapper :class:`~espressomd.io.vtk.VTKReader`:

.. code-block:: python

    import espressomd.io.vtk
    vtk_reader = espressomd.io.vtk.VTKReader()
    vtk_grids = vtk_reader.parse("simulation_step_0.vtu")
    vtk_density = vtk_grids["density"]
    print(vtk_density.shape)

For a self-contained example, please refer to :ref:`LB VTK output`.
