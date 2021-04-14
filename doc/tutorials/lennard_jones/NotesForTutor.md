# Hints for Tutors: Tutorial 01 - Lennard-Jones

## Learning objectives (physics)

After the tutorial, students should be able to 

* explain
    * role of Lennard-Jones potential:
      * noble gases
      * first MD simulations by Verlet et al. (liquid argon)
    * role of radial distribution function
      * relation to structure factor
      * LJ agrees well with experiments for noble gases!
    * how to use an auto correlation function to estimate correlation times and how that affects error estimation


## Learning objectives (ESPResSo)

After the tutorial, students should be able to 

* instance a system object and set basic parameters like the time step
* add particles to a simulation
* access and change particle properties
* set up non-bonded interactions
* explain the role of particle type
* remove overlap using steepest descent
* setting up a thermostat
* integrating and collecting data
* explain the basic concept of observables and accumulators

## Points to mention throughout the tutorial

* Re-explain Langevin equation
* Mention efficiency (loops and lists) vs. numpy arrays at the cost of some readability.
* Non-bonded interactions are defined between pairs of particle types.
  Mention that type=0 is implied when creating particles without specifying one.
* Explain Lennard-Jones cutoff and shift - we don't shift here to make energies comparable to historic simulations (Verlet at al.).
* Explain reduced units - here: LJ units!
* Explain steepest descent: In general overlap removal can
  be done in multiple ways (warmup with force capping, 
  small time step, etc.).
  Additionally the steepest descent algorithm can get trapped in local minima and the convergence criterion is system-dependent.
* Mention that accumulators are updated auomatically.
* Statistically independent measurements for t > correlation time.
