Running the simulation
======================

``integrate``: Running the simulation
-------------------------------------

integrate integrate set integrate set npt\_isotropic

uses the Velocity Verlet algorithm for the integration of the equations
of motion. The command ``integrate`` with an integer ``steps`` as
parameter integrates the system for ``steps`` time steps.

Note that this implementation of the Velocity Verlet algorithm reuses
forces, that is, they are computed once in the middle of the time step,
but used twice, at the beginning and end. However, in the first time
step after setting up, there are no forces present yet. Therefore, has
to compute them before the first time step. That has two consequences:
first, random forces are redrawn, resulting in a narrower distribution
of the random forces, which we compensate by stretching. Second,
coupling forces of e.g. the Lattice Boltzmann fluid cannot be computed
and are therefore lacking in the first half time step. In order to
minimize these effects, has a quite conservative heuristics to decide
whether a change makes it necessary to recompute forces before the first
time step. Therefore, calling hundred times ``integrate 1`` does the
same as ``integrate 100``, apart from some small calling overhead.

However, for checkpointing, there is no way for to tell that the forces
that you read back in actually match the parameters that are set.
Therefore, would recompute the forces before the first time step, which
makes it essentially impossible to checkpoint LB simulations, where it
is vital to keep the coupling forces. To work around this, ``integrate``
has an additional parameter , which tells integrate to not recalculate
the forces for the first time step, but use that the values still stored
with the particles. Use this only if you are absolutely sure that the
forces stored match your current setup!

The opposite problem occurs when timing interactions: In this case, one
would like to recompute the forces, despite the fact that they are
already correctly calculated. To this aim, the option can be used to
enforce force recalculation.

Two methods for the integration can be set: For an NVT ensemble
(thermostat) and for an NPT isotropic ensemble (barostat). The current
method can be detected with the command ``integrate set`` without any
parameters.

The NVT integrator is set without parameters (the temperature can be set
with the thermostat). For the NPT ensemble, the parameters that can be
added are:

-  The external pressure as float variable. This parameter is required.

-  The mass of the applied piston as float variable. This parameter is
   required.

-  :: Three integers to set the box geometry for non-cubic boxes. This
   parameter is optional.

-  ``-cubic_box`` If this optional parameter is added, a cubic box is
   assumed.

``time_integration``: Runtime of the integration loop
-----------------------------------------------------

time\_integration time\_integration

This command runs the integration as would the ``integrate`` command and
returns the wall runtime in seconds.

``minimize_energy``: Run steepest descent minimization
------------------------------------------------------

minimize\_energy

In Python the minimize\_energy functionality can be imported from
``espressomd`` as class ``MinimizeEnergy``. Alternatively it is already
part of the System class object and can be called from there (second
variant).

This command runs a steepest descent energy minimization on the system.
Please note that the behaviour is undefined if either a thermostat,
Maggs electrostatics or Lattice-Boltzmann is activated. It runs a simple
steepest descent algorithm:

Iterate

.. math:: p_i = p_i + \min(\text{\var{gamma}} \times F_i, \text{\var{maxdisplacement}}),

 while the maximal force is bigger than or for at most times. The energy
is relaxed by , while the change per coordinate per step is limited to .
The combination of and can be used to get an poor man’s adaptive update.
Rotational degrees of freedom are treated similarly: each particle is
rotated around an axis parallel to the torque acting on the particle.
Please be aware of the fact that this needs not to converge to a local
minimum in periodic boundary conditions. Translational and rotational
coordinates that are fixed using the “fix“ command or the
ROTATION\_PER\_PARTICLE feature are not altered.

``tune_skin``: Tune the skin
----------------------------

tune\_skin

Determines the fastest skin between and with tolerance by bisection. The
integration time is determined by timing intergration. You should chose
big engough so that multiple verlet updates occure even for the skin,
otherwise the timings are not meaningful. Please be aware that this
command runs actual integrations and propagates the system. In a typical
MD simulation it should be used after warmup and equilibration, in the
same conditions where sampling is done.

``change_volume``: Changing the box volume
------------------------------------------

change\_volume change\_volume

Changes the volume of either a cubic simulation box to the new volume or
its given x-/y-/z-/xyz-extension to the new box-length , and
isotropically adjusts the particles coordinates as well. The function
returns the new volume of the deformed simulation box.

``rotate_system``: Rotating the system around its center of mass
----------------------------------------------------------------

rotate\_system

Rotates the particle coordinates around the system’s center of mass.
This only makes sense for non-periodic boundaries, but no check is
performed. In addition to the particle positions, the command also
rotates the particles themselves, if the ROTATION feature is activated.
Hence, dipole moments as well as virtual sites based on the VS\_RELATIVE
method will also be affected. The command only works on a single cpu.

The rotation axis is given by the parameters and as

.. math:: \vec{a} = (\sin \theta \cos \phi ; \sin \theta \sin \phi ; \cos \theta),

 and denotes the rotation angle.

``lees_edwards_offset``: Applying shear between periodic images
---------------------------------------------------------------

[sec:lees-edwards]

lees\_edwards\_offset

Lees-Edwards Periodic Boundary Conditions are used to impose a shear
flow of speed :math:`\dot{\gamma}` on the system relative to its
periodic images by moving the PBC wrap such that:
:math:`v\_x_{unfolded} =  v\_x_{folded} + \dot{\gamma} y_{img}` (where
:math:`v\_x_{unfolded}` is the :math:`x`-component of the velocity of an
image particle outside the main simulation box, and :math:`y_{img}` is
the count of PBC boundaries crossed in the :math:`y`-direction). The
absolute value of the shear offset is set using this command; with the
shear flow rate :math:`\dot{\gamma}` then determined internally as the
difference between successive offsets. A typical usage would be to
integrate by 1 MD timestep and then to increase the offset to a new
value using this command; this usage pattern is intended to allow for
arbitrary shear flow time profiles, such as an oscillatory shear. A
common calculation to make using Lees-Edwards boundary conditions is to
find the shear viscosity (or kinematic viscosity) by plotting shear
stress (or shear stress/density) against the applied strain for
different values of constant :math:`\dot{\gamma}`.

Lees-Edwards differs from the NEMD approach (see ) in that the shear
imposed is homogenous across the system (but only on average: symmetry
breaking effects are not ruled out) rather than reversing direction with
a periodicity of the box length. Accordingly the transport properties
measured using Lees-Edwards are likely to be different to (and arguably
more physical than) those measured using NEMD or those from equilibrium
simulations by a Green-Kubo type approach.

When the shear flow rate :math:`\dot{\gamma}` is non-zero, the Langevin
thermostat will treat :math:`v\_x_{folded}` as being relative to a flow
field which changes smoothly from :math:`-\dot{\gamma}/2` at the bottom
of the periodic box to :math:`\dot{\gamma}/2` at the top. This ‘laminar’
thermostatting is provided mostly because it gives quite convenient
equilibration of a flowing system. In order to correctly observe
transport properties, symmetry-breaking or entropy production in
relation to shear flow is probably better to use the DPD thermostat (see
) once the initial heat-up has been carried out. The DPD thermostat
removes kinetic energy from the system based on a frictional term
defined relative to a local reference frame of a given particle-pair,
without enforcing any specific flow pattern *a priori*. At high rates of
dissipation, this can however lead to an artefactual shear-banding type
effect at the periodic boundaries, such that the bulk fluid is nearly
stationary. This effect is removed using the modification proposed to
the DPD thermostat by Chatterjee :raw-latex:`\cite{chatterjee2007}` to
allow treatment of systems with high dissipation rates, which is applied
automatically if is compiled in. Chatterjee’s modification is just to
skip calculation of DPD forces (both dissipative and random) for
particle pairs which cross a boundary in *y*.

The function returns the old value of the offset.

If is compiled in, then coordinates are folded into the primary
simulation box as the integration progresses, to prevent a numerical
overflow.

Stopping particles
------------------

Use the following functions, also see Section [sec:Galilei]:

-  ``kill_particle_motion``: halts all particles in the current
   simulation, setting their velocities to zero, as well as their
   angular momentum if the feature ROTATION has been compiled in.

-  ``kill_particle_forces``: sets all forces on the particles to zero,
   as well as all torques if the feature ROTATION has been compiled in.

``velocities``: Setting the velocities
--------------------------------------

velocities

Sets the velocities of the particles with particle IDs between and
:math:`\var{pid}+\var{N}` to a random vector with a length less than ,
and returns the absolute value of the total velocity assigned. By
default, all particles are affected.

Fixing the particle sorting
---------------------------

sort\_particles

Resorts the particles, making sure that

-  the domain decomposition is strictly fullfilled, each particle is on
   the processor and in the cell that its position belongs to

-  the particles within each cell are ordered with ascending identity.

Both conditions together form a unique particle ordering. This is
important when doing checkpointing, because this makes sure that random
numbers are applied in a specific order. Therefore, after writing or
reading a checkpoint, you should call ``sort_particles``.

Parallel tempering
------------------

parallel\_tempering::main -rounds -swap -perform

This command can be used to run a parallel tempering simulation. Since
the simulation routines and the calculation of the swap probabilities
are provided by the user, the method is not limited to sampling in the
temperature space. However, we assume in the following that the sampled
values are temperatures, and call them accordingly. It is possible to
use multiple processors via TCP/IP networking, but the number of
processors can be smaller than the number of temperatures.

specifies the name of the routine calculating the swap probability for a
system. The routine has to accept three parameters: the of the system to
evaluate, and two temperatures and . The routine should return a list
containing the energy of the system at temperatures and , respectively.

specifies the name of the routine performing the simulation between two
swap tries. The routine has to accept two parameters: the of the system
to propagate and the temperature at which to run it. Return values are
ignored.

specifies the name of a routine initializing a system. This routine can
for example create the particles, perform some intial equilibration or
open output files. The routine has to accept two parameters: the of the
system to initialize and its initial temperature . Return values are
ignored.

specifies the number of swap trial rounds; in each round, neighboring
temperatures are tried for swapping alternatingly, i.e. with four
temperatures, The first swap trial round tries to swap
:math:`1\leftrightarrow 2` and :math:`3\leftrightarrow 4`, the second
round :math:`2\leftrightarrow 3`, and so on.

the name of the host on which the parallel\_tempering master node is
running.

the TCP/IP port on which the parallel\_tempering master should listen.
This defaults to 12000.

specifies how many systems to run per -instance. If this is more than 1,
it is the user’s responsibility to manage the storage of configurations,
see below for examples. This defaults to 1.

specifies after how many swap trial rounds to reset the counters for the
acceptance rate statistics. This defaults to 10.

specifies which output the parallel tempering code should produce:

-  parallel tempering will be totally quiet, except for fatal errors

-  information on client activities, such as connecting, is printed to
   stderr

-  print lots of information on swap energies and probabilities to
   stdout. This is useful for debugging and quickly checking the
   acceptance rates.

This defaults to ``all``.

Introduction
^^^^^^^^^^^^

The basic idea of parallel tempering is to run :math:`N` simulations
with configurations :math:`C_i` in parallel at different temperatures
:math:`T_1<T_2<\hdots<T_N`, and exchange configurations between
neighboring temperatures. This is done according to the Boltzmann rule,
the swap probability for two configurations A and B at two different
parameters :math:`T_1` and :math:`T_2` is given by

.. math::

   \label{eq:ptacceptance}
     \min\left(1,\exp -\left[\beta(T_2)U_\text{A}(T_2) + \beta(T_1)U_\text{B}(T_1) -
       \beta(T_1)U_\text{A}(T_1) - \beta(T_2)U_\text{B}(T_2)\right]\right),

 where :math:`U_C(T)` denotes the potential energy of configuration
:math:`C` at parameter :math:`T` and :math:`\beta(T)` the corresponding
inverse temperature. If :math:`T` is the temperature, :math:`U_C` is
indepedent of :math:`T`, and :math:`\beta(T)=1/(k_BT)`. In this case,
the swap probability reduces to the textbook result

.. math:: \min(1,\exp -\left[\left(1/T_2 - 1/T_1\right)\left(U_\text{A} - U_\text{B}\right)/k_B\right].

 However, :math:`T` can also be chosen to be any other parameter, for
example the Bjerrum length, the the strength of the electrostatic
interaction. In this case, :math:`\beta(T)=\beta` is a constant, but the
energy :math:`U_C(T)` of a configuration :math:`C` depends on :math:`T`,
and one needs the full expression . always uses this expression.

In practice, one does not swap configurations, but temperatures, simply
because exchanging temperatures requires much less communication than
exchanging the properties of all particles.

Th implementation of parallel tempering repeatedly propagates all
configurations :math:`C_i` and tries to swap neighboring temperatures.
After the first propagation, the routine attempts to swap temperatures
:math:`T_1` and :math:`T_2`, :math:`T_3` and :math:`T_4`, and so on.
After the second propagation, swaps are attempted between temperatures
:math:`T_2` and :math:`T_3`, :math:`T_4` and :math:`T_5`, and so on. For
the propagation, parallel tempering relies on a user routine; typically,
one will simply propagate the configuration by a few 100 MD time steps.

Details on usage and an example
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The parallel tempering code has to be loaded explicitely by
``source scripts/parallel_tempering.tcl`` from the Espresso directory.
To make use of the parallel tempering tool, one needs to implement three
methods: the propagation, the energy calculation and an initialization
routine for a configuration. A typical initialization routine will look
roughly like this:

::

    proc init {id temp} {
      # create output files for temperature temp
      set f [open "out-$temp.dat" w]; close $f
      init_particle_positions
      thermostat langevin $temp 1.0
      equilibration_integration
      global config
      set config($id) "{[part]} [setmd time]"
    }

The last two lines are only necessary if each instance of handles more
than one configuration, if you have 300 temperatures, but only 10
processes (``-load 30``). In this case, all user provided routines need
to save and restore the configurations. Saving the time is not necessary
because the simulation tine across swaps is not meaningful anyways; it
is however convenient for investigating the (temperature-)history of
individual configurations.

A typical propagation routine accordingly looks like this

::

    proc perform {id temp} {
      global config
      particle delete
      foreach p [lindex $config($id) 0] { eval part $p }
      setmd time [lindex $config($id) 1]
      thermostat langevin $temp 1.0
      set f [open "out-$temp.dat" a];
      integrate 1000
      puts $f "[setmd time] [analyze energy]"
      close $f
      set config($id) "{[part]} [setmd time]"
    }

Again, the saving and storing of the current particle properties in the
config array are only necessary if there is more than one configuration
per process. In practice, one will rescale the velocities at the
beginning of perform to match the current temperature, otherwise the
thermostat needs a short time to equilibrate. The energies necessary to
determine the swap probablility are calculated like this:

::

    proc swap {id temp1 temp2} {
      global config
      particle delete
      foreach p $config($id) { eval part $p }
      set epot [expr [analyze energy total] - [analyze energy kinetic]]
      return "[expr $epot/$temp1] [expr $epot/$temp2]"
    }

Note that only the potential energy is taken into account. The
temperature enters only indirectly through the inverse temperature
prefactor, see Eqn. .

The simulation is then started as follows. One of the processes runs the
command

::

    for {set T 0} {$T < 3} {set T [expr $T + 0.01]} {
        lappend temperatures $T }
    parallel_tempering::main -load 30 -values $temperatures -rounds 1000 \
        -init init -swap swap -perform perform

This command turns the instance executing it into the master part of the
parallel tempering simulation. It waits until a sufficient number of
clients has connected. This are additional instances, which are
identical to the master script, except that they execute

::

    parallel_tempering::main -connect $host -load 30 \
        -init init -swap swap -perform perform

Here, ``host`` is a variable containing the TCP/IP hostname of the
computer running the master process. Note that the master process waits
until enough processes have connected to start the simulation. In the
example, there are 300 temperatures, and each process, including the
master process, will deal with 30 of them. Therefore, 1 master and 9
slave processes are required. For a typical queueing system, a starting
routine could look like this:

::

    master=
    for h in $HOSTS; do
      if [ "$master" == "" ]; then
        ssh $h "cd run; ./pt_test.tcl"
        master=$h;
      else
        ssh $h "cd run; ./pt_test.tcl -connect $host"
      fi
    done

where ``pt_test.tcl`` passes the ``-connect`` option on to
``parallel_tempering::main``.

Sharing data
^^^^^^^^^^^^

parallel\_tempering::set\_shareddata

can be used at any time *by the master process* to specify additional
data that is available on all processes of the parallel\_tempering
simulation. The data is accessible from all processes as
``parallel_tempering::shareddata``.

Metadynamics
------------

metadynamics metadynamics set off metadynamics set distance metadynamics
set relative\_z metadynamics print\_stat current\_coord metadynamics
print\_stat coord\_values metadynamics print\_stat profile metadynamics
print\_stat force metadynamics load\_stat

Required features: METADYNAMICS

Performs metadynamics sampling. Metadynamics is an efficient scheme to
calculate the potential of mean force of a system as a function of a
given reaction coordinate from a canonical simulation. The user first
chooses a reaction coordinate (``distance``) between two particles ( and
). As the system samples values along this reaction coordinate (here the
distance between and ), an iterative biased force pulls the system away
from the values of the reaction coordinate most sampled. Ultimately, the
system is driven in such a way that it self-diffuses along the reaction
coordinate between the two boundaries (here and ). The potential of mean
force (or free energy profile) can be extracted by reading the
``profile``.

ID of the first particle involved in the metadynamics scheme.

ID of the second particle involved in the metadynamics scheme.

: minimum value of the reaction coordinate. While must be positive (it’s
a distance), can be negative since it’s the relative height of with
respect to .

: maximum value of the reaction coordinate.

height of the bias function.

width of the bias function.

strength of the ramping force at the boundaries of the reaction
coordinate interval.

: number of bins of the reaction coordinate. This is only used for the
numerical evaluation of the bias function.

number of relaxation steps before setting a new hill.

Tcl list of a previous metadynamics profile.

Tcl list of a previous metadynamics force.

Details on usage
^^^^^^^^^^^^^^^^

Variant returns the status of the metadynamics routine. Variant turns
metadynamics off (default value). Variant sets a metadynamics scheme
with the reaction coordinate ``distance``, which corresponds to the
distance between any two particles of the system (calculate the
potential of mean force of the end-to-end distance of a polymer).
Variant sets a metadynamics scheme with the reaction coordinate
``relative_z``: relative height (z coordinate) of particle with respect
to (calculate the potential of mean force of inserting one particle
through an interface with center of mass ). Variant prints the current
value of the reaction coordinate. Variant prints a list of the binned
values of the reaction coordinate ( values between and ). Variant prints
the current potential of mean force for all values of the reaction
coordinate considered. Variant prints the current force (norm rather
than vector) for all values of the reaction coordinate considered.
Variant loads a previous metadynamics sampling by reading a Tcl list of
the potential of mean force and applied force. This is especially useful
to restart a simulation.

Note that the metadynamics scheme works seamlessly with the
VIRTUAL\_SITES feature, allowing to define centers of mass of groups of
particles as end points of the reaction coordinate. One can therefore
measure the potential of mean force of the distance between a particle
and a *molecule* or *interface*.

The metadynamics scheme has (as of now) only been implemented for one
processor: MPI usage is *not* supported. However, one can speed up
sampling by communicating the ``profile`` and ``force`` between
independent simulations (denoted *walkers*). The ``print_stat`` and
``load_stat`` can be used to input/output metadynamics information
between walkers at regular intervals. Warning: the information extracted
from ``print_stat`` contains the entire history of the simulation, while
only the *last* increment of sampling should be communicated between
walkers in order to avoid counting the same samples multiple times.

Details on implementation
^^^^^^^^^^^^^^^^^^^^^^^^^

As of now, only two reaction coordinates have been implemented:
``distance`` and ``relative_z``. Many different reaction coordinates can
be set up, and it is rather easy to implement new ones. See the code in
``metadynamics.{h,c}`` for further details.

The bias functions that are applied to the potential of mean force and
the biased force are not gaussian function (as in many metadynamics
codes) but so-called Lucy functions. See :raw-latex:`\cite{marsili09}`
for more details. These avoid the calculation of exponentials.

``integrate_sd``: Running a stokesian dynamics simulation
---------------------------------------------------------

integrate\_sd

uses in the stokesian dynamics algorithm the euler integrater for the
equations of motion. The motion are overdamped. The velocities of the
particles are related to the displacement in the last timestep. This
are, at least if the system is thermalized, however not usefull, as the
displacements don’t scale linear with the timestep. The command
``integrate_sd`` with an integer as parameter integrates the system for
time steps. This is implemented using CUDA, so has to be available in
the system.

Currently there is no parallel implementation of this integrator, so all
particles have to be in a single process.

Setting up the system
~~~~~~~~~~~~~~~~~~~~~

Before running a stokesian dynamics simulation, you have to set a view
constants:

The hydrodynamic particle radius of the (sperical) particles. Only one
particle size is supported.

The viscisity of the fluid. Remember this is only a scaling of the
timestep, so you can set it without problems to one.

(int[2]) seed of the Stokes Dynamics random number generator. As the
generator used for the SD runes on the GPU (cuRAND) it has its own seed
and own state.

(int[2]) offset of the random number generator. Together with the seed,
the state of the random number generator is well defined.

(double) precision used for the approximation of the square root of the
mobility. Sometimes higher accuracy can speedup the simulation.

Make sure to only use the stokesian dynamics thermostat. If you made a
warmup integration without SD, set it with:

thermostat off thermostat sd $temp

Periodicity
~~~~~~~~~~~

The Code uses the ewald-summation (as derived in
:raw-latex:`\cite{beenakker86a}`) of the Rotne-Prager tensor to be
usefull in periodic boundary conditions. To disable the peridoicity of
the hydrodynamic, use the feature .

The ewald summation is required if the system is periodic, as otherwise
an cutoff is introduced, which can lead to negative eigenvalues of the
mobility. This leads to a break down of the code, as it is not possible
to thermalize such a system.

Floatingpoint precision
~~~~~~~~~~~~~~~~~~~~~~~

It is possible to switch between double and float as datatype in
Stokesian Dynamics. As GPUs have more single than double precision
floating point units, the single precision is faster. However this can
lead to problems, if the eigenvalues are not precise enough and get
negative. Therfore the default is double precision. If you want to use
single precision, use the feature .

Farfield only
~~~~~~~~~~~~~

The mobility matrix consists of two contribution, the farfield and the
nearfield lubrication correction. If the nearfield is neglegible, than
the feature can be used. This should speedup the simulation, as the
nearfield doesn’t need to be inverted. Additional the thermal
displacements can be directly calculated.

Multi-timestepping
------------------

| setmd smaller\_time\_step 0.001
| part :math:`i` smaller\_timestep 1

Required features: MULTI\_TIMESTEP

The multi-timestepping integrator allows to run two concurrent
integration time steps within a simulation, associating beads with
either the large or the other . Setting to a positive value turns on the
multi-timestepping algorithm. The ratio / *must* be an integer. Beads
are by default associated with , corresponding to the particle property
0. Setting to 1 associate the particle to the integration. The
integrator can be used in the NVE ensemble, as well as with the Langevin
thermostat and the modified Andersen barostat for NVT and NPT
simulations, respectively. See :raw-latex:`\cite{bereau15}` for more
details.
