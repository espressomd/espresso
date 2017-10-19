.. _Analysis in the core:

Analysis in the core
====================

Analysis in the core is a new concept introduced in since version 3.1.
It was motivated by the fact, that sometimes it is desirable that the
analysis functions do more than just return a value to the scripting
interface. For some observables it is desirable to be sampled every few
integrations steps. In addition, it should be possible to pass the
observable values to other functions which compute history-dependent
quantities, such as correlation functions. All this should be done
without the need to interrupt the integration by passing the control to
the script level and back, which produces a significant overhead when
performed too often.

Some observables in the core have their corresponding counterparts in
the Tcl observables of the ``analyze`` command described in
Chapter [chap:analysis]. However, only the core-observables can be used
on the fly with the toolbox of the correlator and on the fly analysis of
time series. Similarly, some special cases of using the correlator have
their redundant counterparts in the analysis in Tcl
(Chapter [chap:analysis]), but the correlator provides a general and
versatile toolbox which can be used with any implemented
core-observables. The only trick to bridge the gap between Tcl based
analysis and core analysis is the tclcommand observable that allows use
the return value of arbitrary Tcl functions (also self-written) as input
for the core analysis. See more below.

Observables
-----------

Introduction
~~~~~~~~~~~~

The first step of the core analysis is to tell to create an observable.
An observable in the sense of the core analysis can be considered as a
rule how to compute a certain set of numbers from a given state of the
system or a role how to collect data from other observables. Any
observable is represented as a single array of double values. Any more
complex shape (tensor, complex number, …) must be compatible to this
prerequisite. Every observable however documents the storage order.

Creating an observable means just allocating the corresponding memory,
assigning a function to compute the observable value and reserving an
which will be used to refer to the observable. In addition to the
possibility to print the observable value (return the observable value
to the script interface), the of a core-observable can be passed to
another analysis function. The observable value is computed from the
current state of the system at the moment when it is needed, when
requested explicitly by the user calling the ``observable print``
function or when requested automatically by some other analysis
function. Updating is an orthogonal concept: Observables that collect
data over time (e.g. the average observable) need to be updated
regularly, even though their current value is not of interest.

Not all observables are implemented in parallel. When performing a
parallel computation, too frequent updates to observables which are not
implemented in parallel may produce a significant slowdown.

Creating an observable
~~~~~~~~~~~~~~~~~~~~~~

Tcl
^^^

To create a new observable, use

observable new

Upon this call, allocates the necessary amount of memory and returns an
integer which will be used later to refer to the observable. The
parameter and further arguments have to correspond to one of the
observables described below.

Python
^^^^^^

The observables are represented as Python classes. They are contained in
the ``espressomd.observables`` module. An observable is instanced as
follows

::

    from espressomd.observables import *
    part_pos=ParticlePositions(ids=(1,2,3,4,5))

Here, the keyword argument ``ids`` specifies the ids of the particles,
which the observable should take into account.

Available observables
^^^^^^^^^^^^^^^^^^^^^

Tcl
^^^

Currently the following observables are implemented. Particle
specifications (see section [sec:PartSpec] below) define a group of
particles, from which the observable should be calculated. They are
generic to all observables and are described after the list of
observables.

Here are the observables, that only depend on the current state of the
simulation system:

-  | 
   | Positions of the particles, in the format
     :math:`x_1,\ y_1,\ z_1,\ x_2,\ y_2,\ z_2,\ \dots\ x_n,\ y_n,\ z_n`.
     The particles are ordered ascending according to their ids.

-  | 
   | Velocities of the particles, in the format
   | :math:`v^x_1,\ v^y_1,\ v^z_1,\ v^x_2,\ v^y_2,\ v^z_2,\ 
               \dots\ v^x_n,\ v^y_n,\ v^z_n`. The particles are ordered
     ascending according to their ids.

-  | 
   | Velocities of the particles in the body frame, in the format
   | :math:`v^x_1,\ v^y_1,\ v^z_1,\ v^x_2,\ v^y_2,\ v^z_2,\ 
               \dots\ v^x_n,\ v^y_n,\ v^z_n`. The particles are ordered
     ascending according to their ids. This command only produces a
     meaningful result when ``ROTATIONS`` is compiled in.

-  | 
   | Forces on the particles, in the format
   | :math:`f^x_1,\ f^y_1,\ f^z_1,\ f^x_2,\ f^y_2,\ f^z_2,\ 
               \dots\ f^x_n,\ f^y_n,\ f^z_n`. The particles are ordered
     ascending according to their ids.

-  | 
   | Angular momenta (omega) of the particles, in the format
   | :math:`\omega^x_1,\ \omega^y_1,\ \omega^z_1,\ \omega^x_2,\ \omega^y_2,\ \omega^z_2,\ 
               \dots\ \omega^x_n,\ \omega^y_n,\ \omega^z_n`. The
     particles are ordered ascending according to their ids and the
     angular velocity/momentum is specified in the laboratory frame.

-  | 
   | Angular momenta (omega) of the particles, in the format
   | :math:`\omega^x_1,\ \omega^y_1,\ \omega^z_1,\ \omega^x_2,\ \omega^y_2,\ \omega^z_2,\ 
               \dots\ \omega^x_n,\ \omega^y_n,\ \omega^z_n`. The
     particles are ordered ascending according to their ids and the
     angular velocity/momentum is specified in the body (co-rotating)
     frame.

-  | 
   | Position of the centre of mass. If is specified, the particles are
     subdivided into blocks of size and the centre of mass position is
     calculated for each block separately.

-  | 
   | Velocity of the centre of mass. If is specified, the particles are
     subdivided into blocks of size and the centre of mass velocity is
     calculated for each block separately.

-  | 
   | Total force on the specified particles. If is specified, the
     particles are subdivided into blocks of size and the total force is
     calculated for each block separately.

-  | 
   | The stress tensor. It only works with all particles. It is returned
     as a 9-dimensional array:
   | :math:` \{\ \sigma_{xx},\ \sigma_{xy},\ \sigma_{xz},\ \sigma_{yx},\ \sigma_{yy},\
               \sigma_{yz},\ \sigma_{zx},\ \sigma_{zy},\ \sigma_{zz}\ \} `

-  | 
   | The observable for computation of the Stress tensor autocorrelation
     function. Similarly to the stress tensor, it only works with all
     particles. It is returned as a 6-dimensional array:
   | :math:` \{\ \sigma_{xy},\ \sigma_{yz},\ \sigma_{zx},\ 
               ( \sigma_{xx} - \sigma_{yy}),\ 
               ( \sigma_{xx} - \sigma_{zz}),\ 
               ( \sigma_{yy} - \sigma_{zz})\  
               \} `
   | where :math:`\sigma_{ij}` are the components of the stress tensor.

-  | 
   | Electric currents due to individual particles. For a particle
     :math:`i`: :math:`j_i^x = q_i v_i^x / \Delta t` where
     :math:`\Delta t` is the simulation time step. Required feature:

-  | 
   | Electric currents summed over all particles:
     :math:`j^x = \sum_i q_i v_i^x / \Delta t` where :math:`\Delta t` is
     the simulation time step. Required feature:

-  | 
   | The dipole moment of the specified group of particles:
     :math:`\mu^x = \sum_i q_i r_i^x` Required feature:

-  | 
   | Total dipole moment of the specified group of particles:
     :math:`\mu^x = \sum_i d_i^x` where :math:`d_i^x` is an intrinsic
     dipole moment of the individual :math:`i`-th particle. Required
     feature:

-  | 
   | For each particle belonging to the observable is unity if a
     neighbour of a type from is found within the distance defined by
     the . If no such neighbour is found, the observable is zero. The
     observable has one dimension per each particle of

-  | 
   | Compute the density profile within the specified cube. For profile
     specifications, see section [sec:DensProfSpec].

-  | 
   | Compute the force density profile within the specified cube. For
     profile specifications, see section [sec:DensProfSpec].

-  | 
   | Compute the Lattice-Boltzmann velocity profile within the specified
     cube. For profile specifications, see section [sec:DensProfSpec].

-  | 
   | Compute the flux density within the specified cube. For profile
     specifications, see section [sec:DensProfSpec].

-  Compute the density profile in cylindrical coordinates. For profile
   specifications, see section [sec:DensProfSpec].

-  | 
   | Compute the flux density profile in cylindrical coordinates. For
     profile specifications, see section [sec:DensProfSpec].

-  | 
   | Compute the Lattice-Boltzmann velocity profile in cylindrical
     coordinates. For profile specifications, see
     section [sec:DensProfSpec].

-  | 
   | Compute the radial distribution function, see [analyze:rdf].

-  | 
   | Compute the structure factor. Remember it scales as
     order\ :math:`^3`, see [analyze:structurefactor].

-  | ``\``
   | ``\``
   | Computes the radial density distribution for particles of the given
     type around the axis given either as fixed positions or as particle
     ids. The binning is done between ``minr`` and ``maxr`` with
     ``rbins``.

-  | 
   | Calculates the mean charge along ``Npoly`` weak polyelectrolytes of
     the *same* length. The particles can be specified as a list of
     particle ids or through the particle type. The result will be the
     average charge during the simulation in dependence of the monomer
     ranking number, its position on the chain.

-  | 
   | Calculates the persistence length based on the bond vector angle
     correlation of the given polymer specified through the ``id_list``.
     ``max_d`` is the maximum distance (in terms of particles) for which
     the correlation is computed. With ``cut_off`` the number of
     particles at the polymer ends that are ignored for the calculation,
     can be specified. The observable will not calculate the actual
     persistence length but the correlation function of the bond angles.

     .. math::

        \mathrm{Obs}(i) =\mathrm{Norm\_const} \langle \sum_j \frac{(\vec{r}_j - \vec{r}_{j+1})
                    (\vec{r}_{j+i}-\vec{r}_{j+i+1})}{|\vec{r}_j-\vec{r}_{j+1}||\vec{r}_{i+j}-\vec{r}_{i+j+1}|}\rangle

-  | ``\``
   | Computes the pair correlation of particles on a polymer chain given
     through the ``id_list``, here ``k`` is the distance (in terms of
     particles) between the mononers on the chain for which the pair
     correlation is computed. This is averaged over the whole chain and
     all the given chains. The number of polymers for which this
     distribution is computed has to be given as ``Npoly`` and their
     length through ``poly_len``. The distribution is calculated for
     distances between ``minr`` and ``maxr`` with ``rbins``.

     .. math::

        \mathrm{Obs}(r) = \mathrm{Norm\_const} \langle \sum_{j=0,
                    j+=k}^{j<\mathrm{poly\_len} - k} \delta(r - |\vec{r}_j -
                    \vec{r}_{j+k}|)\rangle

The tclcommand observable is a helpful tool, that allows to make the
analysis framework much more versatile, by allowing the evaluation of
arbitrary tcl commands.

-  | 
   | An arbitrary Tcl function that returns a list of floating point
     numbers of fixed size can be specified. Although its execution
     might be slow, it allows to prototype new observables without a lot
     of trouble. Many existing analysis commands can be made to
     cooperate with the core analysis that way.

The following commands allow to collect data automatically over time
once their autoupdate feature is enabled.

-  | 
   | The running average of the reference observable with id . It can be
     resetted by

Printing an observable
~~~~~~~~~~~~~~~~~~~~~~

observable print

Prints the value of the observable with a given :math:`id`. If the
observable refers to the current state of the system, its value is
updated before printing, except if is given.

Python
^^^^^^

The following observables are available

-  Observables working on a given set of particles specified as follows

   ::

       part_vel=ParticleVelocities(ids=(1,2,3,4,5))

   -  ParticlePositions

   -  ParticleVelocities

   -  ParticleForces

   -  ParticleBodyVelocities

   -  ParticleAngularMomentum

   -  ParticleBodyAngularMomentum

   -  ParticleCurrent

   -  Current

   -  DipoleMoment

   -  MagneticDipoleMoment

   -  ComPosition

   -  ComVelocity

   -  ComForce

-  Profile observables sampling the spacial profile of various
   quantities

   ::

       dp =DensityProfile(
          xbins=50, ybins=50, zbins=50, 
          minx=0, miny=0, minz=0,
          maxx=10, maxy=10, maxz=10,
          ids=(1,2,3,4,5))

   -  DensityProfile

   -  FluxDensityProfile

   -  ForceDensityProfile

   -  LBVelocityProfile

Checkpointing observables
~~~~~~~~~~~~~~~~~~~~~~~~~

observable write\_checkpoint

Writes the current state of the observable to . If is given then it does
so in binary form. This is useful to save the state of statful
observables as the average. This checkpoing mechanism is not portable,
rereading is only guaranteed to work on the same system with the same
Espresso binary.

observable read\_checkpoint

Reads the state of the observable from a checkpoint file. The observable
has to be configured beforehand in the script in exactly the same way as
it was when the checkpoint was written. That means that the
configuration of the observable – including possible particle and
profile specifications – is not saved in the file but has to be provided
by the user. Please be aware that the value of observables that refer to
the current state of the system are overwritten by the command except if
the option is given.

Passing an observable to an analysis function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Currently the only analysis function which uses the core observables is
the correlator (section [sec:Correlations]).

Deleting an observable to an analysis function
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

observable delete

Deletes the observable, frees the allocated memory and makes the
:math:`id` free for a new observable.

Particle specifications
~~~~~~~~~~~~~~~~~~~~~~~

You can specify from which particles the observable should be computed
in one of the following ways. In all cases, particle specifications
refer to the current state of espresso. Any later changes to particles
(additions, deletions, changes of types) will not be automatically
reflected in the observable.

-  | 
   | Requests observable calculation based on all particles in the
     system.

-  | 
   | Restricts observable calculation to a given particle type(s). The
     type list is a tcl list of existing particle types.

-  | 
   | Restricts observable calculation to a given list of particle id(s).
     The id list is a tcl list of existing particle ids.

Profile specifications
~~~~~~~~~~~~~~~~~~~~~~

Profiles are specified by giving the spacial area that is to be profiled
and the number of bins in each spacial direction. The area to be
analyzed is characterized by / / and /. The defaults correspond to the
box size when the observable is created. The bin size in each direction
defaults to 1, and can be change with the parameter //. Changing one,
two or three of them to a value :math:`>1` with thus create a one-, two-
or three-dimensional map of the desired quantity. The full syntax thus
reads as:

observable new needs\_profile\_specs

Radial profiles allow to do the same as usual profiles, except the
coordinate system is a cylindrical one and the binning is done in the
cylindrical coordinates (defined with the axis in z-direction). This is
very helpful if the symmetry of the system is cylindrical. The spacial
are is characterized by a center (default to the center of the box) a
maximum radial position (defaults to the smaller value of the box
lengths in x and y directions) and a minimum and maximum value of
:math:`z`. It is possible to also resolve different polar angles, thus
using it as a full 3D mapping tool, but this will only rarely be used.
The full syntax is:

observable new needs\_radial\_profile\_specs

Correlations
------------

[sec:Correlations]

Introduction
~~~~~~~~~~~~

Time correlation functions are ubiquitous in statistical mechanics and
molecular simulations when dynamical properties of many-body systems are
concerned. A prominent example is the velocity autocorrelation function,
:math:` \left< \mathbf{v}(t) \cdot \mathbf{v}(t+\tau) \right> ` which is
used in the Green-Kubo relations. In general, time correlation functions
are of the form

.. math::

   C(\tau) = \left<A\left(t\right) \otimes B\left(t+\tau\right)\right>\,,
   \label{eq:corr.def}

where :math:`t` is time, :math:`\tau` is the lag time (time difference)
between the measurements of (vector) observables :math:`A` and
:math:`B`, and :math:`\otimes` is an operator which produces the vector
quantity :math:`C` from :math:`A` and :math:`B`. The ensemble average
:math:`\left< \cdot \right>` is taken over all time origins \ :math:`t`.
Correlation functions describing dynamics of large and complex molecules
such as polymers span many orders of magnitude, ranging from MD time
step up to the total simulation time.

uses a fast correlation algorithm (see section [sec:multipleTau]) which
enables efficient computation of correlation functions spanning many
orders of magnitude in the lag time.

The generic correlation interface of may process either observables
defined in the kernel, or data which it reads from an external file or
values entered through the scripting interface. Thus, apart from data
processing on the fly, it can also be used as an efficient correlator
for stored data. In all cases it produces a matrix of :math:`n+2`
columns. The first two columns are the values of lag times :math:`\tau`
and the number of samples taken for a particular value of :math:`\tau`.
The remaining ones are the elements of the :math:`n`-dimensional vector
:math:`C(\tau)`.

The command for computing averages and error estimates of a time series
of observables relies on estimates of autocorrelation functions and the
respective autocorrelation times. The correlator provides the same
functionality as a by-product of computing the correlation function (see
section [ssec:CorrError].

An example of the usage of observables and correlations is provided in
the script in the samples directory.

Creating a correlation
~~~~~~~~~~~~~~~~~~~~~~

Correlation first has to be defined by saying which observables are to
be correlated, what should be the correlation operation, sampling
frequency, etc. When a correlation is defined, its id is returned which
is used further to do other operations with the correlation. The
correlation can be either updated automatically on the fly without
direct user intervention, or by an explicit user call for an update.

correlation new obs1 corr\_operation dt tau\_max

Defines a new correlation and returns an integer which has been assigned
to it. Its further arguments are described below.

| and
| are ids of the observables A and B that are to correlated. The ids
  have to refer to existing observables which have been previously
  defined by the command. Some observables are already implemented, and
  others can be easily added. This can be done with very limited
  knowledge just by following the implementations that are already in.
  If is omitted, autocorrelation of is calculated by default.

| 
| The operation that is performed on :math:`A(t)` and :math:`B(t+\tau)`
  to obtain :math:`C(\tau)`. The following operations are currently is
  available:

-  | 
   | Scalar product of :math:`A` and :math:`B`,
     :math:`C=\sum\limits_{i} A_i B_i`

-  | 
   | Comnponentwise product of :math:`A` and :math:`B`,
     :math:`C_i = A_i B_i`

-  | 
   | Each component of the correlation vector is the square of the
     difference between the corresponding components of the observables,
     :math:`C_i = (A_i-B_i)^2`. Example: when :math:`A` is , it produces
     the mean square displacement (for each component separately).

-  | 
   | Tensor product of :math:`A` and :math:`B`,
     :math:`C_{i \cdot l_B + j} = A_i B_j`, with :math:`l_B` the length
     of :math:`B`.

-  
-  | 
   | Fluorescence Correlation Spectroscopy (FCS) autocorrelation
     function,

     .. math::

        G_i(\tau) = \frac{1}{N} \Bigl< \exp \Bigl( - \frac{\Delta x_i^2(\tau) }{w_x^2} - \frac{\Delta y_i^2(\tau)}{w_y^2} - \frac{\Delta z_i^2(\tau)}{w_z^2} \Bigr) \Bigr>\,,
            \label{eq:Gtau}

     where
     :math:`\Delta x_i^2(\tau) = \bigl(x_i(0) - x_i(\tau) \bigr)^2` is
     the square discplacement of particle :math:`i` in the :math:`x`
     direction, and :math:`w_x` is the beam waist of the intensity
     profile of the exciting laser beam,

     .. math:: W(x,y,z) = I_0 \exp \Bigl( - \frac{2x^2}{w_x^2} - \frac{2y^2}{w_y^2} - \frac{2z^2}{w_z^2} \Bigr)\,.

     Equation  is a generalization of the formula presented by Höfling
     :cite:`hofling11a`. For more information, see references therein. Per each
     3 dimensions of the observable, one dimension of the correlation output is
     produced. If is used with other observables than , the physical meaning of
     the result is unclear.

| 
| The time interval of sampling data points. When autoupdate is used,
  has to be a multiple of timestep. It is also used to produce time axis
  in real units. *Warning: if is close to the timestep, autoupdate is
  strongly recommended. Otherwise cpu time is wasted on passing the
  control between the script and kernel.*

| 
| This is the maximum value of :math:`\tau` for which the correlation
  should be computed. *Warning: Unless you are using the multiple tau
  correlator, choosing of more than 100 will result in a huge
  computational overhead. In a multiple tau correlator with reasonable
  parameters, can span the entire simulation without too much additional
  cpu time.*

| 
| The number of data-points for which the results are linearly spaced in
  tau. This is a parameter of the multiple tau correlator. If you want
  to use it, make sure that you know how it works. By default, it is set
  equal to which results in the trivial linear correlator. By setting
  :math:`<` the multiple tau correlator is switched on. In many cases,
  =16 is a good choice but this may strongly depend on the observables
  you are correlating. For more information, we recommend to read
  Ref. :cite:`ramirez10a` or to perform your own tests.

| and
| Are functions used to compress the data when going to the next level
  of the multiple tau correlator. Different compression functions for
  different observables can be specified if desired, otherwise the same
  function is used for both. Default is which takes one of the
  observable values and discards the other one. This is safe for all
  observables but produces poor statistics in the tail. For some
  observables, compression can be used which makes an average of two
  neighbouring values but produces systematic errors. Depending on the
  observable, the systematic error can be anything between harmless and
  disastrous. For more information, we recommend to read
  Ref. :cite:`ramirez10a` or to perform your own tests.

Python
^^^^^^

Each correlator is represented by an instance of the Correlator class,
which is defined in the ``espressomd.correlators`` module.

The meaning of the arguments is as described for TCL. The only
exceptions are the ``obs1`` and ``obs2`` arguments, which take instances
of the Observable class.

Correlators can be registered for automatic updating during the
integration by adding them to ``system.auto_update_correlators``.

::

    system.auto_update_correlators.add(corr)

Inquiring about already existing correlations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

correlation correlation n\_corr

Variant returns a tcl list of the defined correlations including their
parameters.

Variant returns the number of currently defined correlations.

Collecting time series data for the correlation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

correlation autoupdate { start \| stop} correlation update correlation
finalize

Variant is the recommended way of updating the correlations. By
specifying or it starts or stops automatically updating the correlation
estimates. The automatic updates are done within the integration loop
without further user intervention. The update frequency is adjusted
based on the value of provided when defining the correlation. Note that
autoupdate has to be started setting the sim-time (e.g. after ).

Variant is an explicit call for an instantaneous update of the
correlation estimates, using the current system state. It is only
possible to use if the correlation is not being autoupdated. However, it
is possible to use it after autoupdate has been stopped. When updating
by an explicit call, does not check if the lag time between two updates
corresponds the value of specified when creating the correlation.

Variant correlates all data from history which are left in the buffers.
Once this has been done, the history is lost and no further updates are
possible. When a new observable value is passed to a correlation, level
0 of the compression buffers of the multiple tau correlator (see
section [sec:multipleTau] for details) is updated immediately. Higher
levels are updated only when the lower level buffers are filled and
there is a need to push some values one level up. When the updating is
stopped, a number of observable values have not reached the higher
level, especially when is comparable to the total simulation time and if
there are many compression levels. In such case, variant is very useful.
If is much shorter, it does not have a big effect.

Printing out the correlation and related quantities
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

correlation write\_to\_file correlation print correlation print
correlation print

Variant writes the current status of the correlation estimate to the
specified filename. If the file exists, its contents will be
overwritten.

The output looks as follows:

tau1 n\_samples C1 C2 ... Cn tau2 n\_samples C1 C2 ... Cn

Where each line corresponds to a given value of , is the number of
samples which contributed to the correlation at this level and
:math:`C_i` are the individual components of the correlation.

Variant returns the current status of the correlation estimate as a Tcl
variable. The output looks as follows:

 tau1 n\_samples C1 C2 ... Cn tau2 n\_samples C1 C2 ... Cn

| Variants and return the corresponding estimate of the statistical
  property as a Tcl variable.
|  prints the average of observable1.
|  prints the variance of observable1.
|  prints the estimate of the correlation time.
|  prints the estimate of the error of the average based on the method
  according to :cite:`wolff04a` (same as used by the
  command).

The correlation algorithm: multiple tau correlator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here we briefly describe the multiple tau correlator which is
implemented in . For a more detailed description and discussion of its
behaviour with respect to statistical and systematic errors, please read
the cited literature. This type of correlator has been in use for years
in the analysis of dynamic light
scattering :cite:`schatzel88a`. About a decade later it
found its way to the Fluorescence Correlation Spectroscopy
(FCS) :cite:`magatti01a`. The book of Frenkel and
Smit :cite:`frenkel02b` describes its application for the
special case of the velocity autocorrelation function.

.. figure:: figures/correlator_scheme.pdf
   :alt: Schematic representation of buffers in the correlator.

   Schematic representation of buffers in the correlator.

Let us consider a set of :math:`N` observable values as schematically
shown in figure [fig:dataSet], where a value of index :math:`i` was
measured in time :math:`i\delta t`. We are interested in computing the
correlation function according to equation  for a range lag times
:math:`\tau = (i-j)\delta t` between the measurements :math:`i` and
:math:`j`. To simplify the notation, we further drop :math:`\delta t`
when referring to observables and lag times.

The trivial implementation takes all possible pairs of values
corresponding to lag times
:math:`\tau \in [{\tau_{\mathrm{min}}}:{\tau_{\mathrm{max}}}]`. Without
loss of generality, let us further consider
:math:`{\tau_{\mathrm{min}}}=0`. The computational effort for such an
algorithm scales as
:math:`{\cal O} \bigl({\tau_{\mathrm{max}}}^2\bigr)`. As a rule of
thumb, this is feasible if :math:`{\tau_{\mathrm{max}}}< 10^3`. The
multiple tau correlator provides a solution to compute the correlation
functions for arbitrary range of the lag times by coarse-graining the
high :math:`\tau` values. It applies the naive algorithm to a relatively
small range of lag times :math:`\tau \in [0:p-1]`. This we refer to as
compression level 0. To compute the correlations for lag times
:math:`\tau \in [p:2(p-1)]`, the original data are first coarse-grained,
so that :math:`m` values of the original data are compressed to produce
a single data point in the higher compression level. Thus the lag time
between the neighbouring values in the higher compression level
increases by a factor of :math:`m`, while the number of stored values
decreases by the same factor and the number of correlation operations at
this level reduces by a factor of :math:`m^2`. Correlations for lag
times :math:`\tau \in [2p:4(p-1)]` are computed at compression level 2,
which is created in an analogous manner from level 1. This can continue
hierarchically up to an arbitrary level for which enough data is
available. Due to the hierarchical reduction of the data, the algorithm
scales as
:math:`{\cal O} \bigl( p^2 \log({\tau_{\mathrm{max}}}) \bigr)`. Thus an
additional order of magnitude in :math:`{\tau_{\mathrm{max}}}` costs
just a constant extra effort.

The speedup is gained at the expense of statistical accuracy. The loss
of accuracy occurs at the compression step. In principle one can use any
value of :math:`m` and :math:`p` to tune the algorithm performance.
However, it turns out that using a high :math:`m` dilutes the data at
high :math:`\tau`. Therefore :math:`m=2` is hard-coded in the correlator
and cannot be modified by user. The value of :math:`p` remains an
adjustable parameter which can be modified by user by setting when
defining a correlation. In general, one should choose :math:`p \gg m` to
avoid loss of statistical accuracy. Choosing :math:`p=16` seems to be
safe but it may depend on the properties of the analyzed correlation
functions. A detailed analysis has been performed in
Ref. :cite:`ramirez10a`.

The choice of the compression function also influences the statistical
accuracy and can even lead to systematic errors. The default compression
function is which discards the second for the compressed values and
pushes the first one to the higher level. This is robust and can be
applied universally to any combination of observables and correlation
operation. On the other hand, it reduces the statistical accuracy as the
compression level increases. In many cases, the compression operation
can be applied, which averages the two neighbouring values and the
average then enters the higher level, preserving almost the full
statistical accuracy of the original data. In general, if averaging can
be safely used or not, depends on the properties of the difference

.. math::

   \frac{1}{2} (A_i \otimes B_{i+p} + A_{i+1} \otimes B_{i+p+1} ) - 
   \frac{1}{2} (A_i + A_{i+1} ) \otimes \frac{1}{2} (B_{i+p} +  B_{i+p+1})
   \label{eq:difference}

For example in the case of velocity autocorrelation function, the
above-mentioned difference has a small value and a random sign,  
different contributions cancel each other. On the other hand, in the of
the case of mean square displacement the difference is always positive,
resulting in a non-negligible systematic error. A more general
discussion is presented in Ref. :cite:`ramirez10a`.

Checkpointing the correlator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

It is possible to checkpoint the correlator. Thereby the data is written
directly to a file. It may be usefull to write to a binary file, as this
preserves the full bit-value of the variables, whereas the text file
representation has a lower accuracy.

correlation write\_checkpoint\_binary correlation
write\_checkpoint\_ascii

In order to load a checkpoint, the correlator has to be initialized.
Therefore the observable(s) have to be created. Make sure that the
correlator is exactly initilized as it was when the checkpoint was
created. If this is not fullfilled, and e.g. the size of an observable
has changed, loading the checkpoint failes.

correlation read\_checkpoint\_binary correlation read\_checkpoint\_ascii

Depending on whether the checkpoint was written as binary or as text,
the corresponding variant for reading the checkpoint has to be used.

An simple example for checkpointing::

    set pp [observable new particle_positions all]
    set cor1 [correlation new obs1 pp corr_operation square_distance_componentwise \
    dt 0.01 tau_max 1000 tau_lin 16]
    integrate 1000
    correlation cor1 write_checkpoint_binary “cor1.bin”

And then to continue the simulation::

    set pp [observable new particle\_positions all] set cor1 [correlation
    new obs1 pp corr_operation square_distance_componentwise \
    dt 0.01 tau_max 1000 tau_lin 16]
    correlation `\ cor1 read\_checkpoint\_binary “cor1.bin”
