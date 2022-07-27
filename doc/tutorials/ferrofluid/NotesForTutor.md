# Hints for Tutors: Tutorial 11 - Ferrofluids

## General Remarks:
Try to put the focus on the first tutorial.
In 1.5h one should theoretically get through the first two tutorials, but key concepts are mainly introduced in the first one, the second one only introduces an external field.
The third tutorial deals mainly with physical concepts, the only new espresso-specific thing is the setup of a 3D system, i.e. no more DLC.

## Learning Objectives (Physics)
 1. Tutorial:
   - Understand the role and importance of the order parameter lambda ($\lambda \approx 1 \implies$ no large scale structure, $\lambda > 2 \implies$ formation of chain- and ring-like structures)
   - Which cluster distribution is expected (probability decreases with size of clusters)
   
 2. Tutorial:
   - alpha-parameter compares external field coupling to thermal energy
   - External force/field introduces ordering of moments
   - The Langevin magnetization curve gives a theoretical prediction for magnetization of non-interacting dipoles
   - Field perpendicular to particle layer shows noticeable deviations from Langevin-curve due to torque from dipole-dipole interactions

## Learning Objectives (Espresso):
 1. Tutorial: 
   - How to setup P3M for dipolar interactions
   - How to set up randomly distributed dipole moments without bias:
     - Randomly picking cartesian vectors with x, y, z \in [-1,1] would introduce bias
     - Randomly picking angles (phi \in [0,2pi], theta \in [0,pi]) introduces bias
     - Correct way shown in solution (accounts for jacobian)
   - Basics of how to perform cluster analysis
   - If more time available or if asked: What changes with different lambdas
     - This is not documented in the tutorial
     - For a visual understanding of the meaning of lambda, re-run the tutorial with a different value, e.g. lambda = 0.5, 7
     - Discuss the resulting video

 2. Tutorial:
   - How to setup an external field
   - How to sample magnetization curves on a quasi-2D system
   - As sampling of the magnetization curves takes a while, depending on remaining time it might be sensible to only sample the perpendicular direction and skip the cells pertaining to the parallel measurements. Go directly to plotting (of course only perpendicular).

## Points to Mention Throughout the Tutorial:
 1. Tutorial:
   - Importance of correct setup for non-biased results (especially important for magnetic systems!)
   - Using DLC (analogous to ELC) for quasi 2D system with 3D moment directions
   - Idea of a quasi-2D-System with 3D moments
 2. Tutorial:
   - Appeal of simulations as "numerical" experiment and benchmark against theoretical models
