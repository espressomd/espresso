.. _Analysis in the core:

Analysis in the core
====================

Analysis in the core is a new concept introduced in since version 3.1.
It was motivated by the fact, that sometimes it is desirable that the
analysis functions do more than just return a value to the scripting
interface. For some observables it is desirable to be sampled every few
integration steps. In addition, it should be possible to pass the
observable values to other functions which compute history-dependent
quantities, such as correlation functions. All this should be done
without the need to interrupt the integration by passing the control to
the script level and back, which produces a significant overhead when
performed too often.

Some observables in the core have their corresponding counterparts in
the :mod:`espressomd.analysis` module. However, only the core-observables can be used
on the fly with the toolbox of the correlator and on the fly analysis of
time series. 
Similarly, some special cases of using the correlator have
their redundant counterparts in :mod:`espressomd.analysis`,
but the correlator provides a general and
versatile toolbox which can be used with any implemented
core-observables. 

Observables
-----------

Introduction
~~~~~~~~~~~~

The first step of the core analysis is to create an observable.
An observable in the sense of the core analysis can be considered as a
rule how to compute a certain set of numbers from a given state of the
system or a role how to collect data from other observables. Any
observable is represented as a single array of double values. Any more
complex shape (tensor, complex number, …) must be compatible to this
prerequisite. Every observable however documents the storage order.

The observables can be used in parallel simulations. However,
not all observables carry out their calculations in parallel. 
Instead, the entire particle configuration is collected on the head node, and the calculations are carried out there.
This is only performance-relevant if the number of processor cores is large and/or interactions are calculated very frequently.


Creating an observable
~~~~~~~~~~~~~~~~~~~~~~

The observables are represented as Python classes derived from :class:`espressomd.observables.Observable`. They are contained in
the ``espressomd.observables`` module. An observable is instanced as
follows

::

    from espressomd.observables import ParticlePositions
    part_pos=ParticlePositions(ids=(1,2,3,4,5))

Here, the keyword argument ``ids`` specifies the ids of the particles,
which the observable should take into account.

The current value of an observable can be obtained using its calculate()-method::
    
    print(par_pos.calculate())

Available observables
^^^^^^^^^^^^^^^^^^^^^
The following observables are available:

- Observables working on a given set of particles specified as follows

   - ParticlePositions: Positions of the particles, in the format
     :math:`x_1,\ y_1,\ z_1,\ x_2,\ y_2,\ z_2,\ \dots\ x_n,\ y_n,\ z_n`.
     The particles are ordered according to the list of ids passed to the observable.
   - ParticleVelocities: Velocities of the particles in the form
     :math:`v_{x1},\ v_{y1},\ v_{z1},\ v_{x2},\ v_{y2},\ v_{z2},\ \dots\ v_{xn},\ v_{yn},\ v_{zn}`.
     The particles are ordered according to the list of ids passed to the observable.
   - ParticleForces: Forces on the particles in the form
     :math:`f_{x1},\ f_{y1},\ f_{z1},\ f_{x2},\ f_{y2},\ f_{z2},\ \dots\ f_{xn},\ f_{yn},\ f_{zn}`.
   - ParticleBodyVelocities: the particles' velocity in their respective body-fixed frames.
     :math:`v_{x1},\ v_{y1},\ v_{z1},\ v_{x2},\ v_{y2},\ v_{z2},\ \dots\ v_{xn},\ v_{yn},\ v_{zn}`.
     The particles are ordered according to the list of ids passed to the observable.
   - ParticleAngularVelocities: The particles' angular velocities in the space-fixed frame
     :math:`\omega^x_1,\ \omega^y_1,\ \omega^z_1,\ \omega^x_2,\ \omega^y_2,\ \omega^z_2, \dots\ \omega^x_n,\ \omega^y_n,\ \omega^z_n`. 
     The particles are ordered according to the list of ids passed to the observable.
   - ParticleBodyAngularVelocities: As above, but in the particles' body-fixed frame.
   - ParticleCurrent: Product of the particles' velocity and charge
     :math:`m_1 v^x_1, m_1 v^y_1, m_1 v^z_1, \ldots` 
     The particles are ordered according to the list of ids passed to the observable.
   - Current: Total current of the system
     :math:`\sum_i m_i v^x_i, \sum_i m_i v^y_i, \sum_i m_i v^z_i, \ldots` 
   - DipoleMoment: Total electric dipole moment of the system obtained based on unfolded positions
     :math:`\sum_i q_i r^x_i, \sum_i q_i r^y_i, \sum_i q_i r^z_i` 
   - MagneticDipoleMoment: Total magnetic dipole moment of the system based on the :attr:`espressomd.particle_data.ParticleHandle.dip` property.
     :math:`\sum_i \mu^x_i, \sum_i \mu^y_i, \sum_i \mu^z_i` 
   - ComPosition: The system's center of mass based on unfolded coordinates
     :math:`\frac{1}{\sum_i m_i} \left( \sum_i m_i r^x_i, \sum_i m_i r^y_i, \sum_i m_i r^z_i\right)` 
   - ComVelocity: Velocity of the center of mass
     :math:`\frac{1}{\sum_i m_i} \left( \sum_i m_i v^x_i, \sum_i m_i v^y_i, \sum_i m_i v^z_i\right)` 
   - ComForce: Sum of the forces on the particles
     :math:`\sum_i f^x_i, \sum_i f^y_i, \sum_i f^z_i` 

- Profile observables sampling the spacial profile of various
   quantities

   -  DensityProfile

   -  FluxDensityProfile

   -  ForceDensityProfile

   -  LBVelocityProfile


Correlations
------------

Introduction
~~~~~~~~~~~~

Time correlation functions are ubiquitous in statistical mechanics and
molecular simulations when dynamical properties of many-body systems are
concerned. A prominent example is the velocity autocorrelation function,
:math:`\left< \mathbf{v}(t) \cdot \mathbf{v}(t+\tau) \right>` which is
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

A correlator takes one or two observables, obtains values from them during the simulation and 
finally uses a fast correlation algorithm which enables efficient computation 
of correlation functions spanning many orders of magnitude in the lag time.

The implementation for computing averages and error estimates of a time series
of observables relies on estimates of autocorrelation functions and the
respective autocorrelation times. The correlator provides the same
functionality as a by-product of computing the correlation function.

An example of the usage of observables and correlations is provided in
the script in the samples directory.

Creating a correlation
~~~~~~~~~~~~~~~~~~~~~~

Each correlator is represented by an instance of the :class:`espressomd.correlators.Correlator`. Please see its documentation for an explanation of the arguments that have to be passed to the constructor.

Correlators can be registered for automatic updating during the
integration by adding them to :attr:`espressomd.system.System.auto_update_correlators`.

::

    system.auto_update_correlators.add(corr)

Alternatively, an update can triggered by calling the update()-method of the correlator instance. In that case, one has to make sure to call the update in the correct time intervals.


The current on-the-fly correlation result can of a correlator can be obtained using its result()-method.
The final result (including the latest data in the buffers) is obtained using the finalize()-method. After this, no further update of the correlator is possible.

Example: Calculating a particle's diffusion coefficient
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For setting up an observable and correlator to obtain the mean square displacement of particle 0, use::

    pos_obs=ParticlePositions(ids=(0,))
    c_pos = Correlator(obs1=pos_obs, tau_lin=16, tau_max=100., dt=10*dt,
        corr_operation="square_distance_componentwise", compress1="discard1")

To obtain the velocity auto-correlation function of particle 0, use::

    obs=ParticleVelocities(ids=(0,))
    c_vel = Correlator(obs1=vel_obs, tau_lin=16, tau_max=20., dt=dt,
        corr_operation="scalar_product", compress1="discard1")

The full example can be found in `samples/diffusion_coefficient.py`.




The correlation algorithm: multiple tau correlator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Here we briefly describe the multiple tau correlator which is
implemented in . For a more detailed description and discussion of its
behavior with respect to statistical and systematic errors, please read
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
between the neighboring values in the higher compression level
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
can be applied, which averages the two neighboring values and the
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

