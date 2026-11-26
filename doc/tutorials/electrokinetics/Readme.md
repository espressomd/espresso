# Tutorial: electrokinetics

Modelling electrokinetics together with hydrodynamic interactions.

## General Remarks

The tutorial is split into three independent notebooks, one per simulated system.
Each notebook sets up its own `System` object and can be run on its own, but the
basic setup of an electrokinetics simulation (lattice, `EKContainer`, Poisson
solver, `EKSpecies`) is explained in detail only in part 1.

In 90 min one should get through parts 1 and 2. Part 3 does not introduce any new
concepts regarding the coupling of the solvers, only the reaction interface, and
can be treated as a showcase if time runs short.

Part 2 requires ESPResSo to be built with the `WALBERLA_FFT` feature, parts 1 and
3 only need `WALBERLA`.

## Part 1: advection-diffusion in 2D

### Physics learning objectives

After the tutorial, students should be able to:

* Write down the advection-diffusion equation and its fundamental solution
* Explain the meaning of the Péclet number and identify the
  advection-dominated regime
* Explain what numerical diffusion is, why a finite-difference discretization of
  the advection term produces it, and why this rules out the method for pure
  advection problems

### ESPResSo learning objectives

In the course of this tutorial, students should learn to:

* Set up a lattice-Boltzmann fluid with a prescribed, non-fluctuating velocity
  field
* Set up an `EKContainer` with the placeholder Poisson solver `EKNone` for
  uncharged species
* Set up an `EKSpecies` and initialize its density on individual lattice nodes

### Exercises

1. Set up the LB fluid with a constant velocity field (`kT=0`)
2. Set up the uncharged `EKSpecies`, add it to the container, and initialize a
   delta-shaped droplet on a single lattice node

The analysis that follows (mass conservation, center-of-mass drift, and the
Gaussian fit for the effective diffusion coefficient) is generic numpy/scipy
post-processing of the density field, so it is given and explained rather than
left as an exercise.

## Part 2: electroosmotic flow

### Physics learning objectives

After the tutorial, students should be able to:

* Describe the geometry of a slit pore and the origin of the electric double
  layer
* Sketch the counterion density, velocity and shear stress profiles of an
  electroosmotic flow
* Explain why the ion distribution orthogonal to the plates decouples from the
  driving parallel to them
* Contrast electroosmotic flow with pressure-driven (Poiseuille) flow and explain
  the different shapes of the velocity profiles

### ESPResSo learning objectives

In the course of this tutorial, students should learn to:

* Set up the FFT-based Poisson solver `EKFFT` for charged species
* Convert between dynamic and kinematic viscosity when setting up an LB fluid
* Model an immobile surface charge with a second, stationary `EKSpecies`
* Set flux, density and no-slip boundary conditions from `espressomd.shapes`
* Drive a system either with an external electric field on a species or with a
  body force on the fluid

The LB-fluid setup is given rather than exercised here, since it was already an
exercise in part 1; only the dynamic-to-kinematic viscosity conversion is
highlighted.

### Exercises

1. Set up the `EKFFT` Poisson solver and the `EKContainer`
2. Set up the immobile wall charge species (the mobile counterion is given)
3. Set the boundary conditions at the two walls
4. Switch from electroosmotic to pressure-driven driving

## Part 3: reaction in a turbulent flow

### Physics learning objectives

After the tutorial, students should be able to:

* Write down the rate equation of a bulk reaction and relate stoichiometric
  coefficients and partial orders to it
* Describe the Kármán vortex street and the role of the Reynolds number
* Explain how mixing by the flow limits the rate at which a product is formed

### ESPResSo learning objectives

In the course of this tutorial, students should learn to:

* Set up `EKReactant` objects and an `EKBulkReaction`
* Use Dirichlet density boundary conditions as sources and sinks for a species

The thermalized LB fluid (with thermostat) and the rasterization of the cylinder
obstacles are given rather than exercised, since setting up an LB fluid was
covered in part 1 and `add_boundary_from_shape` in part 2. The notebook explains
the two new aspects — fluid thermalization to trigger the instability, and using
fluctuations — inline. The final analysis (transverse-mean density profiles and
the Reynolds number) is generic array post-processing, so it is also given and
explained.

### Exercises

1. Set up the reactants and the bulk reaction
2. Add the educt sources

## Points to mention throughout the tutorial

* The EK solver and the LB solver share the same lattice and are coupled in both
  directions: advection moves the species with the fluid, and the friction force
  density of the species drags the fluid along
* Boundary conditions have to be set separately for the fluid and for every
  species; an impenetrable wall for a species needs both a `FluxBoundary` and a
  `DensityBoundary`
* All quantities are in simulation units; part 2 in particular mixes dynamic
  viscosity (analytical solution) and kinematic viscosity (LB interface)
