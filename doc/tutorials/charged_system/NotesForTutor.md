# Part 1: Charged rod

## Physical learning goals

* The physical system to be investigated is a polyelectrolyte gel: A polymer that is linked into a gel and releases ions into solution
* In the cell model simplification, we only look at one polymer segment and not the entire gel structure
* A further simplification is treating the polymer as a fixed rod
* The question is whether the dissociated ions stay condensed onto polymer or diffuse away. How to define these notions is a research question in itself
* Poisson Boltzmann theory (PB) predicts a solution for the density that has a characteristic length (Manning radius). Ions within that radius are considered condensed
* In the simulations, we have excluded volume, correlations and therefore go beyond the meanfield PB description. In the general case, ions are considered condensed if they are inside the radius at which the logarithmically plotted integrated radial counterion distribution function $P(r)$ has an inflection point
* The PB solution only depends on the product of the Manning parameter $\xi = \lambda l_B/e$ and the counterion valency, where $\lambda$ is the rod line charge density, $l_B$ is the Bjerrum length and $e$ the elementary charge. It is a dimensionless measure for the rod charge compared to the strength of electrostatic interactions in a system at nonzero temperature. In the simulations we see that the result does indeed depend on the individual choice of $\xi$ and the valency. Deviations from PB mainly occur when the counterions are not monovalent.
* By adding concentrated salt we go further beyond the applicability of mean-field theories. The effect we observe is an overcharging, where we see ions accumulating near the positively charged rod
* A semi-analytical solution for the salt case exists (see the referenced paper) but it is not as trivial as in the salt-free case


## Espresso learning goals

* Setting up WCA interactions (should be clear from previous tutorial)
* Using particle properties like ``q``, ``fix`` and ``type``
* Setting up electrostatic methods, the concept of ``actors``
* There can only be one ESPResSo system. For multiple runs, the relevant parts of the system have to be reset manually
