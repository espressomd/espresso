Immersed Boundary Method for soft elastic objects
=================================================

Please contact the Biofluid Simulation and Modeling Group at the
University of Bayreuth if you plan to use this feature.

This section describes an alternative way to include soft elastic
objects somewhat different from the previous chapter. In the Immersed
Boundary Method (IBM), soft particles are considered as an infinitely
thin shell filled with liquid. When the shell is deformed by an external flow it responds by elastic restoring
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

To compute the elastic forces, three new bonded interactions are defined ibm\_triel, ibm\_tribend and ibm\_volCons. ibm\_triel is used to compute elastic shear forces, ibm\_tribend computes out-of-plane bending forces and ibm\_volCons is an artificial force required to ensure proper volume conservation. An example python script can be found in samples/immersed_boundary.
For a more detailed description, see e.g. Guckenberger and Gekle, J. Phys. Cond. Mat. (2017) or contact us.
This feature probably does not work with advanced LB features such electro kinetics or Shan-Chen.
