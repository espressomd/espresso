.. _Under the hood:

Under the hood
==============

-  Implementation issues that are interesting for the user

-  Main loop in pseudo code (for comparison)

.. _Internal particle organization:

Internal particle organization
------------------------------

Since basically all major parts of the main MD integration have to
access the particle data, efficient access to the particle data is
crucial for a fast MD code. Therefore the particle data needs some more
elaborate organization, which will be presented here. A particle itself
is represented by a structure (``Particle``) consisting of several
substructures (e.g. ``ParticlePosition``, ``ParticleForce`` or
``ParticleProperties``), which in turn represent basic physical properties
such as position, force or charge. The particles are organized in one or
more particle lists on each node, called ``CellPList``. The cells are
arranged by several possible systems, as described in :ref:`Cellsystems`.
A cell system defines a way the particles are stored in |es|, i.e.
how they are distributed onto the processor nodes and how they are
organized on each of them. Moreover a cell system also defines
procedures to efficiently calculate the force, energy and pressure for
the short ranged interactions, since these can be heavily optimized
depending on the cell system. For example, the domain decomposition
cellsystem allows an order N interactions evaluation.

Technically, a cell is organized as a dynamically growing array, not as
a list. This ensures that the data of all particles in a cell is stored
contiguously in the memory. The particle data is accessed transparently
through a set of methods common to all cell systems, which allocate the
cells, add new particles, retrieve particle information and are
responsible for communicating the particle data between the nodes.
Therefore most portions of the code can access the particle data safely
without direct knowledge of the currently used cell system. Only the
force, energy and pressure loops are implemented separately for each
cell model as explained above.

The domain decomposition or link cell algorithm is implemented such
that the cells equal the cells, i.e. each cell is a separate particle
list. For an example let us assume that the simulation box has size
:math:`20\times 20\times 20` and that we assign 2 processors to the
simulation. Then each processor is responsible for the particles inside
a :math:`10\times 20\times 20` box. If the maximal interaction range is
1.2, the minimal possible cell size is 1.25 for 8 cells along the first
coordinate, allowing for a small skin of 0.05. If one chooses only 6
boxes in the first coordinate, the skin depth increases to 0.467. In
this example we assume that the number of cells in the first coordinate
was chosen to be 6 and that the cells are cubic. One would then organize
the cells on each node in a :math:`6\times
12\times 12` cell grid embedded at the center of a
:math:`8\times 14 \times
14` grid. The additional cells around the cells containing the particles
represent the ghost shell in which the information of the ghost
particles from the neighboring nodes is stored. Therefore the particle
information stored on each node resides in 1568 particle lists of which
864 cells contain particles assigned to the node, the rest contain
information of particles from other nodes.

Classically, the link cell algorithm is implemented differently. Instead
of having separate particle lists for each cell, there is only one
particle list per node, and the cells actually only contain pointers
to this particle list. This has the advantage that when particles are
moved from one cell to another on the same processor, only the pointers
have to be updated, which is much fewer data (4 rsp. 8 bytes) than the
full particle structure (around 192 bytes, depending on the features
compiled in). The data storage scheme of however requires to always move
the full particle data. Nevertheless, from our experience, the second
approach is 2-3 times faster than the classical one.

To understand this, one has to know a little bit about the architecture
of modern computers. Most modern processors have a clock frequency above
1GHz and are able to execute nearly one instruction per clock tick. In
contrast, the memory runs at a clock speed around 200MHz. Modern
double data rate (DDR) RAM transfers up to 3.2GB/s at this clock speed
(at each edge of the clock signal 8 bytes are transferred). But in
addition to the data transfer speed, DDR RAM has some latency for
fetching the data, which can be up to 50ns in the worst case. Memory is
organized internally in pages or rows of typically 8KB size. The full
:math:`2\times 200` MHz data rate can only be achieved if the access is
within the same memory page (page hit), otherwise some latency has to be
added (page miss). The actual latency depends on some other aspects of
the memory organization which will not be discussed here, but the
penalty is at least 10ns, resulting in an effective memory transfer rate
of only 800MB/s. To remedy this, modern processors have a small amount
of low latency memory directly attached to the processor, the cache.

The processor cache is organized in different levels. The level 1 (L1)
cache is built directly into the processor core, has no latency and
delivers the data immediately on demand, but has only a small size of
around 128KB. This is important since modern processors can issue
several simple operations such as additions simultaneously. The L2 cache
is larger, typically around 1MB, but is located outside the processor
core and delivers data at the processor clock rate or some fraction of
it.

In a typical implementation of the link cell scheme, the order of the
particles is fairly random, determined e.g. by the order in which the
particles are set up or have been communicated across the processor
boundaries. The force loop therefore accesses the particle array in
arbitrary order, resulting in a lot of unfavorable page misses. In the
memory organization of |es|, the particles are accessed in a virtually
linear order. Because the force calculation goes through the cells in a
linear fashion, all accesses to a single cell occur close in time, for
the force calculation of the cell itself as well as for its neighbors.
Using the domain decomposition cell scheme, two cell layers have to be
kept in the processor cache. For 10000 particles and a typical cell grid
size of 20, these two cell layers consume roughly 200 KBytes, which
nearly fits into the L2 cache. Therefore every cell has to be read from
the main memory only once per force calculation.
