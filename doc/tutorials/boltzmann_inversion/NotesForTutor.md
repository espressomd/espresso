# Tutorial: Boltzmann inversion

## Physics learning objectives

* Develop a simple coarse-grained model for charged colloidal particles with explicit counterions and salt using the method of Boltzmann inversion.
* In the coarse-grained model, only the colloidal particles are explicitly modeled, while the effects of the small ions are considered implicitly through an effective coarse-grained potential.
* Students should plot the RDF between colloidal particles in the system with explicit counterions and salt, interpret it physically and calculate the effective potential between colloidal particles from it.
* For consistency, the obtained effective potentials should be compared with analytical prediction by the Debye-Hückel theory, which is a linearized mean field approach.
* Finally, a coarse-grained MD run should be performed where only the colloidal particles interact with an effective coarse-grained potential.

After the tutorial, students should be able to:

* Explain the shorter runtime of a coarse-grained simulation compared to a simulation with explicit salt.
* Explain the discrepancy in the RDFs (especially in the short range) with explicit and implicit salt.
* Limitations of the Boltzmann inversion method, particularly with regards to many-body effects.

## ESPResSo learning objectives

* Setting up a system with colloidal particles of desired size, counterions and salt.
* Calculating the RDF between colloidal particles using observables and accumulators.
* To use the effective potential via a tabulated interaction, while retaining the WCA potential for separations shorter than the tabulated cutoff.
