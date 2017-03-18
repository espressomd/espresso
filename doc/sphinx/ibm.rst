Immersed Boundary Method for soft elastic objects
=================================================

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
